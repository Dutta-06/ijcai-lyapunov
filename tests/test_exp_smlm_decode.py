"""Tests of experiments/exp_smlm_decode.py: the axis/pixel convention of DECODE's coordinates, planted detections through the development-chosen operating points, the packing of the candidates, the raw
consistency check (a tampered stored detection raises), the report's determinism and byte identity of `--from-raw` with the report made at ingest, and, when the stored raw is present, that the stored
md is what the raw gives. Everything but the last test is planted (no SMLM data, no GPU run)."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
os.environ.setdefault('DECODE_OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'smlm'))
import smlm_common as C            # noqa: E402
import smlm_decode_ingest as IN    # noqa: E402
import exp_smlm_decode as D        # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, '..'))
NAME = 'MT0.N1.LD'


def fake_gt(n_frames, per_frame=2, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for f in range(n_frames):
        for _ in range(per_frame):
            rows.append((f, rng.uniform(300, 6100), rng.uniform(300, 6100), rng.uniform(-400, 400), rng.uniform(300, 6000)))
    a = np.array(rows)
    return dict(frame=a[:, 0].astype(np.int64), x_nm=a[:, 1], y_nm=a[:, 2], z_nm=a[:, 3], photons=a[:, 4])


def candidates(gt, a, b, rng, offset=(0.0, 0.0), miss=0.2, n_false=1, jitter=20.0):
    """Candidates of a 'detector' (rows frame, x, y, photons, prob; frames relative to a): each true emitter found with probability 1 - miss (prob high), jittered and shifted, plus false candidates of low prob."""
    rows = []
    sel = (gt['frame'] >= a) & (gt['frame'] < b)
    for f, x, y in zip(gt['frame'][sel], gt['x_nm'][sel], gt['y_nm'][sel]):
        if rng.uniform() > miss:
            rows.append((f - a, x + offset[0] + jitter * rng.standard_normal(), y + offset[1] + jitter * rng.standard_normal(), 3000.0, rng.uniform(0.7, 1.0)))
    for f in range(b - a):
        for _ in range(n_false):
            rows.append((f, rng.uniform(0, 6400), rng.uniform(0, 6400), 500.0, rng.uniform(0.05, 0.4)))
    return np.array(rows)


def planted_inputs(seed=3, offset=(12.0, -9.0)):
    gt = fake_gt(1600)
    rng = np.random.default_rng(seed)
    dev = candidates(gt, 0, 1000, rng, offset=offset)
    test = candidates(gt, 1000, 1500, rng, offset=offset)
    return gt, dev, test


def planted_comparator(gt, seed=5, jitter=45.0, miss=0.35):
    """Per-frame (n, 4) detections of a weaker 'comparison method' on the test frames."""
    rng = np.random.default_rng(seed)
    rows = candidates(gt, 1000, 1500, rng, miss=miss, n_false=0, jitter=jitter)
    off, flat = D.pack_cand(rows, 500)
    return [np.c_[flat[off[i]:off[i + 1], :3], np.zeros(off[i + 1] - off[i])] for i in range(500)]


def planted_raw():
    gt, dev, test = planted_inputs()
    meta = dict(tag='DECODE', sets=[NAME], decode_commit='planted', software='planted', gpu='planted GPU', budget_text='planted', protocol_text=['planted protocol'], versions_text=['planted versions'])
    sel = json.dumps(dict(chosen_snapshot='snap_0003.pt', swap_chosen=True, best_dev_E_by_axis_convention={'x=DECODE_y (swap)': 61.0, 'x=DECODE_x (no swap)': -1e9},
                          all=[dict(snap='snap_0003.pt', swap=True, n_candidates=len(dev), best=dict(thr=0.5, J=0.8, rmse=20.0, E=61.0))]))
    tr = json.dumps(dict(seconds=12.0, epochs_done=3, cfg=dict(epochs=4, emitter_av=3.0, photons_mean=4690.0, photons_sd=3100.0, bg_photons=99.9)))
    comps = {NAME: {'ISS-E': planted_comparator(gt), 'R2E+c': planted_comparator(gt, seed=6, jitter=60.0)}}
    sets_in = {NAME: dict(gt=gt, dev=dev, test=test, pmin=1851.0, select=sel, train=tr)}
    return D.build_raw(sets_in, meta, comps), gt, dev, test


# ---------------------------------------------------------------------------------------------------------------------------------
def test_decode_pixel_convention_rows_first_and_corner_origin():
    """An emitter in image row 10, column 40 is DECODE (x, y) = (10, 40) (x runs along the rows): the ground-truth convention gives x_nm = 4050 (columns), y_nm = 1050 (rows)."""
    x_nm, y_nm = D.decode_px_to_nm(10.0, 40.0)
    assert x_nm == pytest.approx(4050.0) and y_nm == pytest.approx(1050.0)
    xs, ys = D.decode_px_to_nm(np.array([0.0, 63.0]), np.array([5.0, 0.0]))
    assert list(xs) == [550.0, 50.0] and list(ys) == [50.0, 6350.0]


def test_candidate_packing_round_trip_and_range_check():
    rng = np.random.default_rng(0)
    rows = np.c_[rng.integers(0, 20, 50), rng.uniform(0, 6400, (50, 2)), rng.uniform(100, 5000, 50), rng.uniform(0.05, 1, 50)]
    off, flat = D.pack_cand(rows, 20)
    assert off[0] == 0 and off[-1] == 50 and flat.shape == (50, 4) and np.all(np.diff(off) >= 0)
    back = D.unpack_cand(off, flat)
    order = np.argsort(rows[:, 0], kind='stable')
    assert np.array_equal(back, rows[order])
    assert np.array_equal(np.diff(off), np.bincount(rows[:, 0].astype(int), minlength=20))
    with pytest.raises(ValueError):
        D.pack_cand(np.array([[25.0, 1.0, 2.0, 3.0, 0.5]]), 20)
    s = D.to_struct(rows)
    assert s['prob'].shape == (50,) and np.array_equal(s['frame'], rows[:, 0])
    off0, flat0 = D.pack_cand(np.zeros((0, 5)), 4)
    assert list(off0) == [0, 0, 0, 0, 0] and flat0.shape == (0, 4)


def test_planted_detections_through_the_development_chosen_operating_points():
    gt, dev, test = planted_inputs()
    res = D.set_result(NAME, dev, test, D.gt_array(gt))
    for op in D.OPS:
        o = res[op]
        assert 0.4 <= o['thr'] <= 0.7                                    # false candidates have prob < 0.4, the true ones > 0.7
        assert abs(o['off'][0] - 12.0) < 3 and abs(o['off'][1] + 9.0) < 3
        t = o['test']
        assert t['precision'] > 0.98 and 0.7 < t['recall'] < 0.9 and t['rmse_nm'] < 40
    # the offset is removed, so the RMSE is the jitter (20 nm per axis, ~28 nm), not the offset (15 nm more in quadrature)
    assert res['E']['test']['rmse_nm'] < 35
    # an uncorrected pixel-convention slip would show as a large bias before the correction
    gt2, dev2, test2 = planted_inputs(seed=4, offset=(50.0, 50.0))
    r = D.set_result(NAME, dev2, test2, D.gt_array(gt2))
    assert max(abs(r['E']['uncorrected_bias'])) > 30


def test_scorer_on_planted_detections_matches_hand_counts():
    """Two frames of a ground truth with three emitters; the detections: one exact, one 100 nm off, one beyond the 250 nm tolerance, one false -> TP 2, FP 2, FN 1, RMSE sqrt((0 + 100^2)/2)."""
    gt = D.gt_array(dict(frame=[0, 0, 1], x_nm=[1000.0, 3000.0, 2000.0], y_nm=[1000.0, 3000.0, 2000.0], z_nm=[0.0, 0.0, 0.0], photons=[1000.0, 1000.0, 1000.0]))
    dets = [np.array([[1000.0, 1000.0, 1.0, 0.0], [3100.0, 3000.0, 1.0, 0.0]]), np.array([[2400.0, 2000.0, 1.0, 0.0], [5000.0, 5000.0, 1.0, 0.0]])]
    r = C.score(dets, gt, np.array([0, 1]))
    assert (r['tp'], r['fp'], r['fn']) == (2, 2, 1)
    assert r['jaccard'] == pytest.approx(2 / 5) and r['rmse_nm'] == pytest.approx(np.sqrt(100.0 ** 2 / 2))
    assert r['efficiency'] == pytest.approx(100 - np.hypot(100 - 40.0, np.sqrt(5000.0)))


def test_raw_check_catches_a_tampered_stored_detection():
    raw, *_ = planted_raw()
    D.check_raw(raw, [NAME])
    bad = dict(raw)
    k = '%s|DECODE-E|test|flat' % NAME
    f = np.array(raw[k])
    f[0, 0] += 1.0
    bad[k] = f
    with pytest.raises(ValueError):
        D.check_raw(bad, [NAME])
    bad2 = dict(raw)
    k2 = '%s|DECODE|dev|cand|flat' % NAME
    f2 = np.array(raw[k2])
    f2[:, 3] = np.minimum(f2[:, 3], 0.1)           # the development candidates all low-probability: the chosen thresholds move and the stored test detections no longer match
    bad2[k2] = f2
    with pytest.raises(ValueError):
        D.check_raw(bad2, [NAME])


def test_report_is_deterministic_has_the_expected_content_and_from_raw_is_byte_identical(tmp_path):
    raw, gt, dev, test = planted_raw()
    text = D.report(raw, boot=60)
    assert text == D.report(raw, boot=60)
    assert 'DECODE-E - ISS-E' in text and 'DECODE-E - R2E+c' in text and 'planted protocol' in text and 'Test frames 1000-1499' in text
    # the DECODE E point beats the weak planted comparators, and the report says so with the verdict
    line = [l for l in text.splitlines() if l.startswith('| DECODE - ISS-E')][0]
    assert line.endswith('win |')
    # a save / load round trip and the command line give the same bytes
    out = tmp_path / 'out'
    out.mkdir()
    C.save_npz_deterministic(str(out / 'decode_raw.npz'), raw)
    first = open(str(out / 'decode_raw.npz'), 'rb').read()
    C.save_npz_deterministic(str(out / 'decode_raw2.npz'), D.load_npz(str(out / 'decode_raw.npz')))
    assert first == open(str(out / 'decode_raw2.npz'), 'rb').read()
    txt1 = D.report(D.load_npz(str(out / 'decode_raw.npz')), boot=60)
    assert txt1 == text


def test_from_raw_command_line_regenerates_the_md_exactly(tmp_path):
    raw, *_ = planted_raw()
    out = tmp_path / 'o'
    out.mkdir()
    C.save_npz_deterministic(str(out / 'decode_raw.npz'), raw)
    env = dict(os.environ, DECODE_OUT=str(out), OMP_NUM_THREADS='1')
    for i in range(2):
        p = subprocess.run([sys.executable, os.path.join(HERE, 'exp_smlm_decode.py'), '--from-raw'], env=env, capture_output=True, text=True, timeout=1500)
        assert p.returncode == 0, p.stderr
        if i == 0:
            first = open(str(out / 'decode.md'), 'rb').read()
    assert first == open(str(out / 'decode.md'), 'rb').read() and b'DECODE-E - R2E+c' in first


def candidates_ej(gt, a, b, rng):
    """Candidates whose E-best and J-best thresholds differ: 70 % of the emitters are found accurately (prob ~0.9, jitter 10 nm), 30 % are found inaccurately and displaced (prob ~0.5, +40 nm in x, jitter 85 nm per axis),
    all with the constant offset (12, -9) nm. Including the second class raises J and ruins the RMSE: J picks a threshold <= 0.5, E a threshold >= 0.6."""
    rows = []
    sel = (gt['frame'] >= a) & (gt['frame'] < b)
    for f, x, y in zip(gt['frame'][sel], gt['x_nm'][sel], gt['y_nm'][sel]):
        if rng.uniform() < 0.7:
            rows.append((f - a, x + 12.0 + 10.0 * rng.standard_normal(), y - 9.0 + 10.0 * rng.standard_normal(), 3000.0, rng.uniform(0.85, 0.95)))
        else:
            rows.append((f - a, x + 52.0 + 85.0 * rng.standard_normal(), y - 9.0 + 85.0 * rng.standard_normal(), 3000.0, rng.uniform(0.45, 0.55)))
    return np.array(rows)


def planted_inputs_ej(seed=7, n_frames=1600):
    gt = fake_gt(n_frames)
    rng = np.random.default_rng(seed)
    return gt, candidates_ej(gt, 0, 1000, rng), candidates_ej(gt, 1000, 1500, rng)


def planted_raw_ej(names=(NAME,)):
    gt, dev, test = planted_inputs_ej()
    meta = dict(tag='DECODE', sets=list(names), decode_commit='planted', software='planted', gpu='planted GPU', budget_text='planted', protocol_text=['planted protocol'], versions_text=['planted versions'])
    sets_in = {n: dict(gt=gt, dev=dev, test=test, pmin=1851.0, select='', train='') for n in names}
    comps = {n: {'ISS-E': planted_comparator(gt), 'R2E+c': planted_comparator(gt, seed=6, jitter=60.0)} for n in names}
    return D.build_raw(sets_in, meta, comps), gt, dev, test


def write_run_dir(path, gt, dev, test, epochs=7):
    """A planted return of the GPU run: candidate CSVs, selection and training records, versions, GPU, timing."""
    os.makedirs(path, exist_ok=True)
    for split, rows in (('dev', dev), ('test', test)):
        D_rows = np.asarray(rows)
        with open(os.path.join(path, '%s_%s.csv' % (NAME, split)), 'w') as f:
            f.write('frame,x_nm,y_nm,photons,prob\n')
            for r in D_rows:
                f.write('%d,%.3f,%.3f,%.3f,%.5f\n' % (int(r[0]), r[1], r[2], r[3], r[4]))
    json.dump(dict(chosen_snapshot='snap_0005.pt', swap_chosen=True, best_dev_E_by_axis_convention={'x=DECODE_y (swap)': 60.0, 'x=DECODE_x (no swap)': -3.0},
                   all=[dict(snap='snap_0005.pt', swap=True, n_candidates=len(dev), best=dict(thr=0.5, J=0.8, rmse=20.0, E=60.0))]), open(os.path.join(path, '%s_select.json' % NAME), 'w'))
    json.dump(dict(seconds=99.0, epochs_done=epochs, cfg=dict(budget_s=100.0, seed=0, emitter_av=3.0, z_hi=600.0, max_targets=40, epochs=100000, photons_mean=4690.0, photons_sd=3100.0, bg_photons=99.9)),
              open(os.path.join(path, '%s_train_summary.json' % NAME), 'w'))
    open(os.path.join(path, 'versions.txt'), 'w').write('decode 0.10.1dev1 torch 1.7.1 numpy 1.19.5 scipy 1.8.0\nspline 0.10.0\n')
    open(os.path.join(path, 'gpu.txt'), 'w').write('name, driver_version, memory.total [MiB]\nTesla T4, 580.1, 15360 MiB\nTesla T4, 580.1, 15360 MiB\n')
    json.dump(dict(total_s=1234.0), open(os.path.join(path, 'timing.json'), 'w'))


def write_learned_raw(path, gt):
    """The few keys of learned_raw.npz that the ingest reads, planted."""
    arrays = {}
    for k, (tag, jit, miss) in enumerate((('ISS-0', 40.0, 0.3), ('ISS-E', 30.0, 0.3), ('R0', 50.0, 0.4), ('R2E+c', 45.0, 0.4))):
        dets = planted_comparator(gt, seed=100 + k, jitter=jit, miss=miss)
        off, flat = C.pack_detections(dets)
        arrays['%s|cmp|%s|off' % (NAME, tag)], arrays['%s|cmp|%s|flat' % (NAME, tag)] = off, flat
    arrays['%s|cmp|ISS-E|config' % NAME] = np.array(json.dumps(dict(cfg='hp5h', kappa=5)))
    arrays['%s|pmin' % NAME] = np.array(1851.0)
    C.save_npz_deterministic(path, arrays)


def test_ingest_end_to_end_on_a_planted_run_directory(tmp_path):
    gt, dev, test = planted_inputs()
    run = str(tmp_path / 'run')
    write_run_dir(run, gt, dev, test)
    lr = str(tmp_path / 'learned_raw.npz')
    write_learned_raw(lr, gt)
    out = str(tmp_path / 'decode_raw.npz')
    msgs = []
    D.ingest(run, [NAME], out_path=out, log=msgs.append, gt_loader=lambda name: D.gt_array(gt), learned_path=lr)
    raw = D.load_npz(out)
    meta = json.loads(str(raw['meta']))
    assert meta['sets'] == [NAME] and 'torch 1.7.1' in meta['software'] and meta['gpu'].startswith('Tesla T4') and '7' in meta['budget_text']
    assert any('4c0f38c33681d1c3005bc18bd99a11f46d30f1c1' in l for l in meta['versions_text']) and any('1234' in l for l in meta['versions_text'])
    # the candidates in the raw are those of the CSVs (to the printed precision), the comparison rows are there, the report is consistent
    back = D.cand_of(raw, NAME, 'test')
    assert back.shape == test.shape and np.allclose(back[:, 1:3], test[np.argsort(test[:, 0], kind='stable')][:, 1:3], atol=1e-3)
    backd = D.cand_of(raw, NAME, 'dev')
    assert backd.shape == dev.shape and np.allclose(backd[:, 1:3], dev[np.argsort(dev[:, 0], kind='stable')][:, 1:3], atol=1e-3)
    lz = D.load_npz(lr)
    for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
        got, want = D.cmp_dets(raw, NAME, tag), C.unpack_detections(lz['%s|cmp|%s|off' % (NAME, tag)], lz['%s|cmp|%s|flat' % (NAME, tag)])
        assert len(got) == 500 and all(np.array_equal(x, y) for x, y in zip(got, want))
    assert float(raw['%s|pmin' % NAME]) == 1851.0
    nrows = int((gt['frame'] < 1500).sum())                     # the ground truth of the development and the test frames, nothing later
    assert len(raw['%s|gt|frame' % NAME]) == nrows and raw['%s|gt|frame' % NAME].max() < 1500
    D.check_raw(raw, [NAME])
    text = D.report(raw, boot=40)
    assert 'DECODE-E - ISS-E' in text and 'Tesla T4' in text and 'snap_0005.pt' in text
    # a missing candidate file is an error, not a silent skip
    os.remove(os.path.join(run, '%s_test.csv' % NAME))
    with pytest.raises(Exception):
        D.ingest(run, [NAME], out_path=out, log=msgs.append, gt_loader=lambda name: D.gt_array(gt), learned_path=lr)


def test_compose_meta_reports_the_range_of_training_times_and_the_epochs(tmp_path):
    for name, secs, ep in (('MT0.N1.LD', 4100.0, 80), ('MT0.N1.HD', 4900.0, 61)):
        os.makedirs(str(tmp_path), exist_ok=True)
        json.dump(dict(seconds=secs, epochs_done=ep, cfg=dict(budget_s=5000.0, seed=0, emitter_av=3.0, z_hi=600.0, max_targets=40)), open(str(tmp_path / ('%s_train_summary.json' % name)), 'w'))
    open(str(tmp_path / 'gpu.txt'), 'w').write('name, driver_version, memory.total [MiB]\nTesla T4, 580.1, 15360 MiB\nTesla T4, 580.1, 15360 MiB\n')
    open(str(tmp_path / 'versions.txt'), 'w').write('decode 0.10.1dev1 torch 1.7.1 numpy 1.19.5 scipy 1.8.0\nh5py                    3.7.0\npytorch-1.7.1-py3.8_cuda11.0.221_cudnn8.0.5_0.json\n')
    json.dump(dict(versions_text=['total GPU seconds of the account: 12345']), open(str(tmp_path / 'run_extra.json'), 'w'))
    meta = D.compose_meta(str(tmp_path), ['MT0.N1.LD', 'MT0.N1.HD'])
    assert meta['budget_text'] == '4100-4900 s (80/61 epochs of 10,000 simulated frames)' and meta['gpu'].startswith('Tesla T4;')
    assert meta['versions_text'][0] == 'DECODE commit 4c0f38c33681d1c3005bc18bd99a11f46d30f1c1 (decode 0.10.1dev1 torch 1.7.1 numpy 1.19.5 scipy 1.8.0); packages: h5py 3.7.0; pytorch-1.7.1-py3.8_cuda11.0.221_cudnn8.0.5_0'
    assert 'GPU: Tesla T4, 580.1, 15360 MiB / Tesla T4, 580.1, 15360 MiB' in meta['versions_text']
    assert meta['versions_text'][-1] == 'total GPU seconds of the account: 12345' and any(l.startswith('MT0.N1.HD: 61 epochs in 4900 s (budget 5000 s)') for l in meta['versions_text'])


def test_compose_meta_handles_a_run_without_training_records(tmp_path):
    meta = D.compose_meta(str(tmp_path), [NAME])
    assert meta['budget_text'] == '?' and meta['sets'] == [NAME] and meta['decode_commit'].startswith('4c0f38c3')


def test_e_and_j_points_are_distinct_stored_separately_and_checked():
    raw, gt, dev, test = planted_raw_ej()
    res = D.check_raw(raw, [NAME])[NAME]
    assert res['E']['thr'] >= 0.6 and res['J']['thr'] <= 0.5
    assert res['E']['test']['rmse_nm'] < 25 and res['J']['test']['rmse_nm'] > 45 and res['J']['test']['jaccard'] > res['E']['test']['jaccard']
    assert not np.allclose(res['E']['off'], res['J']['off'])
    # the method table carries the E point under DECODE-E and the J point under DECODE-J, and the stored arrays are those of the right point
    ms = D.method_table(raw, NAME, res)
    assert ms[0][0] == 'DECODE-E' and ms[0][2] is res['E']['test_dets'] and ms[1][0] == 'DECODE-J' and ms[1][2] is res['J']['test_dets']
    for op in D.OPS:
        off, flat = IN.pack(res[op]['test_dets'])
        assert np.array_equal(raw['%s|DECODE-%s|test|off' % (NAME, op)], off) and np.array_equal(raw['%s|DECODE-%s|test|flat' % (NAME, op)], flat)
    # tampering with either point or with the offsets is caught
    for k in ('%s|DECODE-J|test|flat', '%s|DECODE-E|test|flat'):
        bad = dict(raw)
        f = np.array(raw[k % NAME])
        f[-1, 1] += 1.0
        bad[k % NAME] = f
        with pytest.raises(ValueError):
            D.check_raw(bad, [NAME])
    for k in ('%s|DECODE-J|test|off', '%s|DECODE-E|test|off'):
        bad = dict(raw)
        o = np.array(raw[k % NAME])
        o[10:] += 1
        bad[k % NAME] = o
        with pytest.raises(ValueError):
            D.check_raw(bad, [NAME])
    # the raw keeps the photon threshold of the challenge setting and the ground truth of the development and test frames in the given order
    assert float(raw['%s|pmin' % NAME]) == 1851.0
    assert np.array_equal(raw['%s|gt|frame' % NAME], gt['frame']) and np.array_equal(raw['%s|gt|x_nm' % NAME], gt['x_nm'])
    assert raw['%s|cmp|ISS-E|flat' % NAME].shape[1] == 4


def test_summary_and_threshold_table_use_the_e_point_and_the_mean_row_needs_two_sets():
    raw, gt, dev, test = planted_raw_ej(names=(NAME, 'MT0.N2.LD'))
    res = D.check_raw(raw, [NAME, 'MT0.N2.LD'])
    text = D.report(raw, boot=40)
    idx = np.arange(1000, 1500)
    g = D.gt_array(gt)
    eE = C.score(res[NAME]['E']['test_dets'], g, idx)['efficiency']
    eI = C.score(D.cmp_dets(raw, NAME, 'ISS-E'), g, idx)['efficiency']
    line = [l for l in text.splitlines() if l.startswith('| DECODE - ISS-E')][0]
    assert line.split('|')[3].strip().startswith('%.1f' % (eE - eI))
    assert 'mean E over the 2 sets' in text
    mean_row = [l for l in text.splitlines() if l.startswith('| mean E over the 2 sets')][0]
    eR2 = C.score(D.cmp_dets(raw, NAME, 'R2E+c'), g, idx)['efficiency']
    assert mean_row.split('|')[2].strip() == '%.1f / %.1f / %.1f' % (eE, eI, eR2)
    tally_line = [l for l in text.splitlines() if l.startswith('Tally over the sets')][0]
    parts = []
    for kb in ('ISS-E', 'R2E+c'):
        ci = C.bootstrap_diff_ci(C.score(res[NAME]['E']['test_dets'], g, idx)['per_frame'], C.score(D.cmp_dets(raw, NAME, kb), g, idx)['per_frame'], B=40, seed=0, block=D.BLOCK)['efficiency']
        v = 'win' if ci[0] > 0 else ('loss' if ci[1] < 0 else 'tie')
        parts.append('against %s: %d win(s), %d tie(s), %d loss(es)' % (kb, 2 * (v == 'win'), 2 * (v == 'tie'), 2 * (v == 'loss')))
    assert tally_line.endswith(': ' + '; '.join(parts) + '.')
    one = D.report(planted_raw_ej()[0], boot=40)
    assert 'mean E over' not in one
    # emitters per frame of the test frames: the ground truth outside the test frames does not count
    assert 'Test frames 1000-1499 (500 frames, 2.000 ground-truth emitters per frame)' in one
    # the threshold table: the efficiency at each threshold with the E point's offset removed, the E threshold in bold, the best in hindsight
    lines = D.threshold_lines(raw, NAME, g[(g['frame'] >= 1000) & (g['frame'] < 1500)], idx, res[NAME])
    row = [l for l in '\n'.join(lines).splitlines() if l.startswith('| DECODE, offset of the E point')][0]
    exp = [C.score(IN.shift(IN.per_frame(D.to_struct(D.cand_of(raw, NAME, 'test')), 500, t), res[NAME]['E']['off']), g, idx)['efficiency'] for t in IN.THRESHOLDS]
    cells = [c.strip() for c in row.strip('|').split('|')]
    assert cells[1:1 + len(exp)] == [('**%.1f**' if t == res[NAME]['E']['thr'] else '%.1f') % e for t, e in zip(IN.THRESHOLDS, exp)]
    best = int(np.nanargmax(exp))
    assert cells[-1] == '%.2f (E %.1f)' % (IN.THRESHOLDS[best], exp[best])


def test_verdict_boundaries_and_pack_range_checks():
    assert D.verdict(1.0, (0.0, 2.0)) == 'tie' and D.verdict(-1.0, (-2.0, 0.0)) == 'tie'
    assert D.verdict(1.0, (0.1, 2.0)) == 'win' and D.verdict(-1.0, (-2.0, -0.1)) == 'loss'
    rows = np.array([[19.0, 1.0, 2.0, 3.0, 0.5]])
    D.pack_cand(rows, 20)
    for bad in (-1.0, 20.0):
        with pytest.raises(ValueError):
            D.pack_cand(np.array([[bad, 1.0, 2.0, 3.0, 0.5]]), 20)
    assert D.pack_cand(np.array([[0.0, 1.0, 2.0, 3.0, 0.5]]), 20)[0][1] == 1


def test_comparators_from_a_planted_learned_raw_apply_the_threshold_and_the_offset(tmp_path):
    gt, dev, test = planted_inputs()
    p = str(tmp_path / 'learned_raw.npz')
    write_learned_raw(p, gt)
    arrays = D.load_npz(p)
    rng = np.random.default_rng(11)
    cands = [np.c_[rng.uniform(0, 6400, (6, 2)), np.full(6, 100.0), np.array([0.1, 0.9, 0.5, 0.2, 0.7, 0.6])] for _ in range(500)]
    off, flat = C.pack_detections(cands)
    for key, tag in (('sup', 'sup'), ('sim', 'sim')):
        arrays['%s|%s|off' % (NAME, key)], arrays['%s|%s|flat' % (NAME, key)] = off, flat
        arrays['%s|%s|sel' % (NAME, key)] = np.array(json.dumps(dict(thr=0.55, offset=[10.0, -5.0], tag=tag)))
    C.save_npz_deterministic(p, arrays)
    comps, z = D.comparators_from_learned(NAME, p)
    a, b = comps['CNN-a'], comps['CNN-b']
    for i in (0, 7, 499):
        keep = cands[i][cands[i][:, 3] >= 0.55]
        assert len(keep) == 3 and np.allclose(b[i][:, :2], keep[:, :2])                      # (b): threshold only
        assert np.allclose(a[i][:, 0], keep[:, 0] - 10.0) and np.allclose(a[i][:, 1], keep[:, 1] + 5.0)    # (a): threshold and offset removed
    for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
        assert len(comps[tag]) == 500
    assert not np.array_equal(comps['ISS-E'][0], comps['ISS-0'][0])


def test_comparators_from_learned_have_the_decode_format():
    p = os.path.join(REPO, 'results', 'smlm', 'learned_raw.npz')
    if not os.path.exists(p):
        pytest.skip('no learned_raw.npz')
    comps, z = D.comparators_from_learned(NAME, p)
    for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
        assert len(comps[tag]) == 500 and comps[tag][0].shape[1] == 4


def test_stored_md_is_what_the_stored_raw_gives():
    raw_p, md_p = os.path.join(REPO, 'results', 'smlm', 'decode_raw.npz'), os.path.join(REPO, 'results', 'smlm', 'decode.md')
    if not (os.path.exists(raw_p) and os.path.exists(md_p)):
        pytest.skip('no stored DECODE run')
    assert D.report(D.load_npz(raw_p)) == open(md_p).read()
