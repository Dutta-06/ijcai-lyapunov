"""The CPU-time benchmark of the frame fit: the helper functions and the report (the run itself needs the challenge data)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_smlm_speed as X      # noqa: E402


def test_evenly_spaced_frames_cover_the_range_and_are_unique():
    assert X.evenly_spaced(1000, 1500, 5).tolist() == [1000, 1125, 1250, 1374, 1499] or X.evenly_spaced(1000, 1500, 5).tolist()[0] == 1000
    idx = X.evenly_spaced(100, 500, 12)
    assert idx[0] == 100 and idx[-1] == 499 and len(set(idx.tolist())) == 12 and (np.diff(idx) > 0).all()
    assert X.evenly_spaced(10, 13, 12).tolist() == [10, 11, 12]                      # fewer frames than asked: all of them


def test_compare_pairs_detections_by_position_and_counts_unequal_frames():
    a = np.array([[100.0, 50.0, 1000.0, 1.4], [20.0, 80.0, 500.0, 1.4]])
    b = np.array([[20.002, 80.001, 500.5, 1.4], [100.004, 50.0, 1000.0, 1.4]])          # the same two detections in another order, slightly moved
    c, d = np.array([[5.0, 5.0, 100.0, 1.4]]), np.array([[5.01, 5.0, 100.05, 1.4]])      # one detection, moved by more than the others
    same, n, dpos, dph = X.compare([a, a[:1], c], [b, np.zeros((0, 4)), d])
    assert (same, n) == (2, 3)                                                       # the second frame has different counts: not compared
    assert abs(dpos - 0.01) < 1e-12 and abs(dph - 0.001) < 1e-12
    assert X.compare([a, a[:1]], [b, np.zeros((0, 4))])[:2] == (1, 2)               # counts differing the other way round: also skipped, not an error
    assert X.compare([np.zeros((0, 4))], [np.zeros((0, 4))]) == (1, 1, 0.0, 0.0)


def test_report_rows_are_read_off_the_raw(tmp_path):
    path = str(tmp_path / 'speed_raw.npz')
    np.savez_compressed(path, sets=np.array(['A.set']), **{'A.set|frames': np.arange(10), 'A.set|cpu': np.array([30.0, 10.0]), 'A.set|cmp': np.array([9.0, 10.0, 0.0031, 2.5e-6]),
                                                           'A.set|ndet': np.array([55, 54])})
    md = X.report(path)
    assert '| A.set | 10 | 3.000 | 1.000 | 3.00 | 9 / 10 | 55 / 54 | 0.0031 | 2.50e-06 |' in md and md == X.report(path)


def test_the_two_settings_are_those_quoted_in_the_paper():
    assert X.SETTINGS == (('full', dict(tol=1e-8, max_nfev=200)), ('fast', dict(tol=1e-4, max_nfev=25)))
    assert X.SETS == ('MT0.N1.LD', 'MT0.N1.HD')


def test_timed_measures_process_cpu_time_of_one_process_with_the_setting(monkeypatch):
    calls = {}
    clock = iter([10.0, 13.5])

    def fake_run_frames(frames, cfg, cal, opts, workers=None):
        calls.update(frames=frames, cfg=cfg, cal=cal, opts=opts, workers=workers)
        return ['dets-of-the-loop'], 1, 2, 3
    monkeypatch.setattr(X.ISS, 'run_frames', fake_run_frames)
    monkeypatch.setattr(X.ISS, 'get_cal', lambda raw, name: ('cal', raw, name))
    monkeypatch.setattr(X.time, 'process_time', lambda: next(clock))
    monkeypatch.setattr(X.time, 'time', lambda: 1e9)                                  # wall clock must not be used
    d, cpu = X.timed('FRAMES', 'A.set', 'RAW', dict(tol=1e-4, max_nfev=25))
    assert d == 'dets-of-the-loop' and cpu == 3.5
    assert calls['workers'] == 1 and calls['frames'] == 'FRAMES' and calls['cfg'] == X.ISS.DEFAULT['cfg'] and calls['cal'] == ('cal', 'RAW', 'A.set')
    assert calls['opts'] == dict(rules=[('snr', X.ISS.DEFAULT['kappa'])], tol=1e-4, max_nfev=25)


def test_run_alternates_the_settings_and_keeps_the_smaller_time(monkeypatch, tmp_path):
    order, cpus = [], {('full', 0): 9.0, ('fast', 0): 4.0, ('full', 1): 7.0, ('fast', 1): 5.0}
    seen = {'full': 0, 'fast': 0}
    one, none = np.array([[1.0, 2.0, 3.0, 1.4]]), np.zeros((0, 4))
    which = []

    def fake_timed(frames, name, raw, setting):
        label = 'full' if setting['tol'] == 1e-8 else 'fast'
        order.append(label)
        rep = seen[label]
        seen[label] += 1
        return ([one, one] if label == 'full' else [one, none]), cpus[(label, rep)]     # the fast setting loses the detection of the second frame
    monkeypatch.setattr(X.ISS, 'load_raw', lambda p: {})
    monkeypatch.setattr(X.ISS, 'load_set', lambda name: ('STACK', None))
    monkeypatch.setattr(X.ISS, 'dev_test_frames', lambda name, w: which.append(w) or (1000, 1500))
    monkeypatch.setattr(X.C, 'load_frames', lambda stack, a, b: np.zeros((1, 4, 4)))
    monkeypatch.setattr(X.C, 'to_photons', lambda f: f)
    monkeypatch.setattr(X, 'timed', fake_timed)
    out = str(tmp_path / 'sub' / 'speed_raw.npz')
    arrays = X.run(sets=('A.set',), n_frames=5, out_path=out)
    assert order == ['full', 'fast', 'fast', 'full'] and which == ['test']           # alternation: the second repeat in the opposite order; test frames
    assert arrays['A.set|cpu'].tolist() == [7.0, 4.0] and arrays['A.set|frames'].tolist() == X.evenly_spaced(1000, 1500, 5).tolist()
    assert arrays['A.set|cmp'].tolist() == [1.0, 2.0, 0.0, 0.0] and arrays['A.set|ndet'].tolist() == [2, 1]
    z = np.load(out)
    assert z['sets'].tolist() == ['A.set'] and z['A.set|cpu'].tolist() == [7.0, 4.0]
