"""The SMLM spot check of the cycle defect (experiments/exp_cycle_scan.py): the comparison of detections, the grouping of the guard's records by frame and tile, the selection of the frames that are
fitted again without the guard, write-once running, a real fit of two synthetic frames through the experiments' own function, and the report read off a planted raw."""
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_astro_cyclefix as CF            # noqa: E402
import exp_cycle_scan as E                 # noqa: E402
import exp_smlm_iss as ISS                 # noqa: E402
from src.core import imaging as IM         # noqa: E402


def _rec(entries=5, exits=0, act=0, bars=None):
    return dict(counts=np.array([entries, exits, entries, 1, 0, entries + 1, act, 0]), barred=np.array(bars or [], float).reshape(-1, 3), taus=np.arange(float(entries)))


def test_compare_dets_is_bitwise_and_counts_the_emitters_only_one_run_has():
    a = np.array([[100.0, 200.0, 5000.0, 1.4], [300.0, 400.0, 4000.0, 1.5]])
    assert E.compare_dets(a, a.copy()) == (True, 2, 2, 0, 0)
    assert E.compare_dets(a, a[:1]) == (False, 2, 1, 1, 0) and E.compare_dets(a[:1], a) == (False, 1, 2, 0, 1)
    b = a.copy()
    b[1, 2] += 1e-9                                                                      # the photons differ: not identical, but the same two emitters (positions)
    assert E.compare_dets(a, b) == (False, 2, 2, 0, 0)
    c = a.copy()
    c[0, 0] += 5.0                                                                       # one emitter moved by 5 nm: one only in each
    assert E.compare_dets(a, c) == (False, 2, 2, 1, 1)
    assert E.compare_dets(a, np.zeros((0, 4))) == (False, 2, 0, 2, 0) and E.compare_dets(np.zeros((0, 4)), np.zeros((0, 4))) == (True, 0, 0, 0, 0)
    d = a.copy()
    d[0, 1] += E.TOL_NM * 0.5
    assert E.compare_dets(a, d)[3:] == (0, 0) and E.compare_dets(a, d)[0] is False


def test_per_frame_groups_the_records_of_the_tiles_of_a_frame_and_prefixes_the_tile_of_a_barred_point():
    recs = [_rec(5), _rec(6, 1, 1, [[0.3, 4, 5]]), _rec(7), _rec(8), _rec(9, 2, 2, [[0.5, 1, 2], [0.7, 3, 4]]), _rec(10)]
    out = E.per_frame(recs, 2, 3)
    assert [o[:3] for o in out] == [(5 + 6 + 7, 1, 1), (8 + 9 + 10, 2, 2)]
    assert out[0][3].tolist() == [[1, 0.3, 4, 5]] and out[1][3].tolist() == [[1, 0.5, 1, 2], [1, 0.7, 3, 4]] and out[0][3].shape == (1, 4)
    assert E.per_frame([_rec()], 1, 1)[0][3].shape == (0, 4)
    with pytest.raises(RuntimeError, match='5 path calls for 2 frames x 3 tiles'):
        E.per_frame(recs[:5], 2, 3)


def test_the_tiles_of_a_frame_follow_the_recipe_of_the_set():
    assert E.tiles_per_frame('Tubulins_I', (256, 256)) == 16 and E.tiles_per_frame('Bundled_Tubes_Long_Sequence', (64, 64)) == 1 and E.tiles_per_frame('Bundled_Tubes_High_Density', (64, 64)) == 1


def test_a_real_fit_of_two_synthetic_frames_runs_the_guard_and_returns_one_record_per_frame(monkeypatch):
    monkeypatch.setattr(IM, 'sliding_path', IM.sliding_path)
    rng = np.random.default_rng(1)
    P = IM.PixelPSF(64, 64, 1.4)
    frames = np.stack([100.0 + P.synth(np.array([[20.3, 30.1], [44.0, 12.5]]), np.array([900.0, 500.0])) + rng.normal(0, 1.0, (64, 64)) for _ in range(2)])
    cal = dict(bg=np.full((64, 64), 100.0), a=1.0)
    d, recs, secs = E.fit('2016', 'MT0.N1.LD', frames, 'hp3', cal, 5.0, True)
    assert len(d) == 2 and len(recs) == 2 and secs.shape == (2,) and all(r[2] == 0 and r[0] >= 3 for r in recs) and all(x.shape[1] == 4 for x in d) and sum(len(x) for x in d) >= 2
    assert IM.sliding_path is CF.ORIG_PATH and CF.TRACE == []                           # the frozen function is back in place afterwards
    d2, none, _ = E.fit('2016', 'MT0.N1.LD', frames, 'hp3', cal, 5.0, False)
    assert none is None and all(np.array_equal(x, y) for x, y in zip(d, d2))            # no double exit: the guarded fit is the unguarded one
    with pytest.raises(RuntimeError):
        ISS_calls = []
        monkeypatch.setattr(ISS, 'run_frames', lambda *a, **k: (ISS_calls.append(1) or [[np.zeros((0, 4))] * 2], [], [], np.zeros(2)))
        E.fit('2016', 'MT0.N1.LD', frames, 'hp3', cal, 5.0, True)                         # a fit that does not call the path function once per frame is refused


def test_scan_fits_everything_with_the_guard_and_the_frames_with_an_activation_and_the_first_controls_without(monkeypatch):
    dets = [np.array([[10.0 * f, 5.0, 100.0, 1.4]]) for f in range(8)]
    calls = []

    def fake_fit(year, name, frames, cfg, cal, kappa, guard):
        calls.append((guard, len(frames)))
        if guard:
            recs = [(5 + f, 1 if f in (2, 5) else 0, 1 if f in (2, 5) else 0, np.array([[0, 0.3, 4, 5]]) if f in (2, 5) else np.zeros((0, 4))) for f in range(len(frames))]
            return [d.copy() for d in dets], recs, np.full(len(frames), 2.0)
        out = [d.copy() for d in dets]
        out[5] = np.vstack([out[5], [[99.0, 9.0, 50.0, 1.4]]])                           # frame 5 differs by one emitter when unguarded
        idx = [int(frames[j, 0]) for j in range(len(frames))]
        return [out[i] for i in idx], None, np.full(len(frames), 3.0)
    monkeypatch.setattr(E, 'fit', fake_fit)
    frames = np.arange(8.0)[:, None] * np.ones((8, 3))
    arr = E.scan('2016', 'MT0.N1.LD', 'E', frames, {}, 'hp3', 5.0, (1000, 1008), log=lambda *a: None)
    k = lambda f: E.key('2016', 'MT0.N1.LD', 'E', f)          # noqa: E731
    assert calls == [(True, 8), (False, 5)] and arr[k('act')].tolist() == [0, 0, 1, 0, 0, 1, 0, 0] and arr[k('entries')].tolist() == [5, 6, 7, 8, 9, 10, 11, 12]
    assert arr[k('cmp_frames')].tolist() == [2, 5, 0, 1, 3] and arr[k('cmp_why')].tolist() == ['activation', 'activation', 'control', 'control', 'control']
    assert arr[k('cmp')].tolist() == [[1, 1, 1, 0, 0], [0, 1, 2, 0, 1], [1, 1, 1, 0, 0], [1, 1, 1, 0, 0], [1, 1, 1, 0, 0]]
    assert arr[k('cmp10')].tolist() == [[0, 0], [0, 1], [0, 0], [0, 0], [0, 0]]                                         # the same emitters at the physical tolerance of 10 nm
    assert arr[k('barred')].tolist() == [[2, 0, 0.3, 4, 5], [5, 0, 0.3, 4, 5]] and arr[k('cfg')].tolist() == ['hp3', '5.0'] and arr[k('frames')].tolist() == [1000, 1008] and arr[k('secs')].tolist() == [2.0, 3.0]


def test_scan_without_any_activation_compares_only_the_controls(monkeypatch):
    monkeypatch.setattr(E, 'fit', lambda year, name, frames, cfg, cal, kappa, guard: ([np.zeros((0, 4))] * len(frames), [(4, 0, 0, np.zeros((0, 4)))] * len(frames) if guard else None, np.ones(len(frames))))
    arr = E.scan('2013', 'Tubulins_I', 'default', np.zeros((6, 2, 2)), {}, 'hp5h', 5.0, (200, 206), log=lambda *a: None)
    k = lambda f: E.key('2013', 'Tubulins_I', 'default', f)          # noqa: E731
    assert arr[k('cmp_frames')].tolist() == [0, 1, 2] and arr[k('barred')].shape == (0, 5) and arr[k('cmp')].tolist() == [[1, 0, 0, 0, 0]] * 3 and arr[k('cmp10')].tolist() == [[0, 0]] * 3


def test_run_runs_each_set_and_operating_point_once_and_never_again(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(E, 'setup_2016', lambda name, n, raw=None: ('frames', 'cal', {'E': ('hpE', 6.0), 'default': ('hpD', 5.0)}, (1000, 1000 + n)))
    monkeypatch.setattr(E, 'setup_2013', lambda name, n, raw=None: ('frames13', 'cal13', {'E': ('hpE13', 6.0), 'default': ('hpD13', 5.0)}, (200, 200 + n)))
    monkeypatch.setattr(E, 'scan', lambda year, name, op, frames, cal, cfg, kappa, span, log=print: seen.append((year, name, op, cfg, kappa, span)) or {E.key(year, name, op, 'act'): np.zeros(2)})
    out = str(tmp_path / 'raw.npz')
    E.run(['MT0.N1.LD', 'Tubulins_I'], 7, ('E', 'default'), out, lambda *a: None)
    assert seen == [('2016', 'MT0.N1.LD', 'E', 'hpE', 6.0, (1000, 1007)), ('2016', 'MT0.N1.LD', 'default', 'hpD', 5.0, (1000, 1007)), ('2013', 'Tubulins_I', 'E', 'hpE13', 6.0, (200, 207)),
                   ('2013', 'Tubulins_I', 'default', 'hpD13', 5.0, (200, 207))]
    seen.clear()
    E.run(['MT0.N1.LD', 'Tubulins_I'], 7, ('E', 'default'), out, lambda *a: None)
    assert seen == []                                                                   # results exist: nothing is run
    E.run(['MT0.N1.LD'], 7, ('E',), out, lambda *a: None)
    assert seen == []
    E.run(['Tubulins_II'], 7, ('default',), out, lambda *a: None)
    assert [l[:3] for l in seen] == [('2013', 'Tubulins_II', 'default')] and sorted(E.load_raw(out)) == sorted(
        [E.key('2016', 'MT0.N1.LD', o, 'act') for o in ('E', 'default')] + [E.key('2013', 'Tubulins_I', o, 'act') for o in ('E', 'default')] + [E.key('2013', 'Tubulins_II', 'default', 'act')])


def _planted():
    raw = {}
    for year, name, op, act, chg in (('2016', 'MT0.N1.LD', 'E', [0, 2, 0, 0], False), ('2016', 'MT0.N1.LD', 'default', [0, 0, 0, 0], False), ('2013', 'Bundled_Tubes_High_Density', 'E', [1, 0, 3, 0], True)):
        a = np.array(act)
        comp = list(np.flatnonzero(a > 0)) + list(np.flatnonzero(a == 0)[:E.CONTROL])
        why = ['activation'] * int(np.sum(a > 0)) + ['control'] * (len(comp) - int(np.sum(a > 0)))
        cmp = np.ones((len(comp), 5), int)
        if chg:
            cmp[0] = [0, 40, 41, 0, 1]
        raw.update({E.key(year, name, op, 'cfg'): np.array(['hp5h', '6.0']), E.key(year, name, op, 'frames'): np.array([1000, 1004]), E.key(year, name, op, 'act'): a, E.key(year, name, op, 'entries'): np.full(4, 6),
                    E.key(year, name, op, 'exits'): a.copy(), E.key(year, name, op, 'secs'): np.array([1.5, 2.5]), E.key(year, name, op, 'cmp_frames'): np.array(comp, int), E.key(year, name, op, 'cmp_why'): np.array(why),
                    E.key(year, name, op, 'cmp'): cmp, **({E.key(year, name, op, 'cmp10'): np.array([[30, 20]] + [[0, 0]] * (len(comp) - 1))} if chg else {}), E.key(year, name, op, 'barred'): np.array([[f, 0, 0.25 * (f + 1), 3, 4] for f in np.flatnonzero(a > 0) for _ in range(int(a[f]))], float).reshape(-1, 5)})
    return raw


def test_the_report_counts_activations_and_lists_the_frames_whose_detections_changed():
    raw = _planted()
    md = E.report(raw)
    assert md == E.report(raw)
    assert '| 2016 | MT0.N1.LD | E-optimal | hp5h, kappa 6.0 | 4 | 1 | 2 | 6.0 | 2 | 4 | 0 | 0 | 1.50 |' in md
    assert '| 2016 | MT0.N1.LD | a priori default | hp5h, kappa 6.0 | 4 | 0 | 0 | 6.0 | 0 | 3 | 0 | 0 | 1.50 |' in md
    assert '| 2013 | Bundled_Tubes_High_Density | E-optimal | hp5h, kappa 6.0 | 4 | 2 | 4 | 6.0 | 4 | 4 | 1 | 0 | 1.50 |' in md
    assert '**12 frames scanned (all sets and both operating points), 3 of them with at least one activation, 1 frames whose detections changed (8 control frames compared).**' in md
    assert '| 2013 | Bundled_Tubes_High_Density | E | 1000 | activation | 40 | 41 | 0 | 1 | 30 | 20 |' in md                 # frame (stack index) = first frame + position
    no_fine = {k: v for k, v in raw.items() if not k.endswith('|cmp10')}
    assert '| 2013 | Bundled_Tubes_High_Density | E | 1000 | activation | 40 | 41 | 0 | 1 | - | - |' in E.report(no_fine)        # a raw without the physical tolerance
    assert '| 2016 | MT0.N1.LD | E | 1001 | 2 | 0.5 (3, 4); 0.5 (3, 4) |' in md and 'None among the frames compared' not in md
    assert 'Tubulins_I' not in md                                                                                       # a set without results is not listed
    none = {k: v for k, v in raw.items() if not k.startswith('2013')}
    md2 = E.report(none)
    assert 'None among the frames compared.' in md2 and '2 frames' not in md2 or 'None among' in md2
    assert '**8 frames scanned' in md2


def test_from_raw_regenerates_the_report_byte_for_byte(tmp_path):
    raw = _planted()
    np.savez_compressed(tmp_path / 'cycle_scan_raw.npz', **raw)
    env = dict(os.environ, CYCLE_SCAN_OUT=str(tmp_path), PYTHONDONTWRITEBYTECODE='1')
    first = None
    for _ in range(2):
        r = subprocess.run([sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_cycle_scan.py'), '--from-raw'], env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        data = (tmp_path / 'cycle_scan.md').read_bytes()
        assert first is None or data == first
        first = data
        (tmp_path / 'cycle_scan.md').unlink()
    assert first == E.report(E.load_raw(str(tmp_path / 'cycle_scan_raw.npz'))).encode()


def test_the_command_line(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(E, 'run', lambda sets, nframes, ops: calls.append((sets, nframes, ops)))
    monkeypatch.setattr(sys, 'argv', ['x', '--run', '--sets', 'Tubulins_I', '--frames', '5', '--ops', 'E'])
    E.main()
    monkeypatch.setattr(sys, 'argv', ['x', '--run'])
    E.main()
    assert calls == [(['Tubulins_I'], 5, ('E',)), (None, E.N_FRAMES, E.OPS)]
    monkeypatch.setattr(sys, 'argv', ['x'])
    with pytest.raises(SystemExit):
        E.main()
    monkeypatch.setattr(E, 'MD', str(tmp_path / 'r.md'))
    monkeypatch.setattr(E, 'RAW', str(tmp_path / 'raw.npz'))
    np.savez_compressed(tmp_path / 'raw.npz', a=np.zeros(1))
    monkeypatch.setattr(E, 'report', lambda raw: 'the report of %s\n' % sorted(raw))
    monkeypatch.setattr(E, 'load_raw', lambda path=None: {'a': 1})
    monkeypatch.setattr(sys, 'argv', ['x', '--from-raw'])
    E.main()
    assert (tmp_path / 'r.md').read_text() == "the report of ['a']\n"
