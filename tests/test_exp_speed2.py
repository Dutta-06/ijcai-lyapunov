"""experiments/exp_speed2.py: the comparison functions, the ratio and slope statistics, the table builder and the report (from a small constructed raw), the plan's tolerances, and the stored report."""
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, '..')
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'experiments'))
import exp_speed2 as S            # noqa: E402
import smlm_common as C           # noqa: E402


def det(*rows):
    return np.array(rows, float).reshape(-1, 4)


def test_evenly_spaced():
    assert list(S.evenly_spaced(0, 1000, 5)) == [0, 250, 500, 749, 999]
    assert list(S.evenly_spaced(100, 105, 24)) == [100, 101, 102, 103, 104]            # fewer frames than asked: all of them
    assert S.evenly_spaced(10, 20, 1).tolist() == [10]


def test_compare_dets_pairs_by_assignment_not_by_order():
    a = det((100, 200, 1000, 1.4), (900, 50, 500, 1.5))
    b = a[::-1] + np.array([0.01, -0.02, 0.0, 0.0])
    same, dpos, dph, dsig = S.compare_dets(a, b)
    assert same and np.isclose(dpos, 0.02, atol=0, rtol=1e-12) and dph == 0.0 and dsig == 0.0


def test_compare_dets_pairs_detections_that_share_a_column():
    a = det((100, 200, 1000, 1.4), (100, 900, 500, 1.5))
    b = a[::-1] + np.array([0.0, 0.01, 0.0, 0.0])
    assert S.compare_dets(a, b)[1] == pytest.approx(0.01, rel=1e-9)


def test_compare_dets_measures_each_column_with_its_own_scale():
    a = det((100, 200, 1000, 1.4))
    assert S.compare_dets(a, det((100, 200, 1010, 1.4)))[2] == pytest.approx(0.01, rel=1e-12)          # photons relative to the FIRST argument
    assert S.compare_dets(a, det((100, 200, 990, 1.4)))[2] == pytest.approx(0.01, rel=1e-12)
    assert S.compare_dets(a, det((100, 200, 1000, 1.45)))[3] == pytest.approx(0.05, rel=1e-12)
    assert S.compare_dets(a, det((103, 204, 1000, 1.4)))[1] == pytest.approx(4.0, rel=1e-12)           # the larger of |dx|, |dy|
    assert S.compare_dets(det((100, 200, 0.0, 1.4)), det((100, 200, 1.0, 1.4)))[2] == pytest.approx(1.0 / 1e-9, rel=1e-9)     # a zero photon count is guarded, not divided by


def test_compare_dets_counts_and_empties():
    a = det((1, 2, 3, 4))
    ok, p, f, s = S.compare_dets(a, det((1, 2, 3, 4), (5, 6, 7, 8)))
    assert ok is False and np.isnan(p) and np.isnan(f) and np.isnan(s)
    assert S.compare_dets(det(), det()) == (True, 0.0, 0.0, 0.0)
    assert S.compare_dets(a, det())[0] is False


def test_summarize_compare_takes_maxima_over_equal_count_frames_only():
    rows = [(True, 0.01, 1e-5, 1e-6), (False, np.nan, np.nan, np.nan), (True, 0.05, 1e-6, 2e-6)]
    r = S.summarize_compare(rows)
    assert r.tolist() == [3, 2, 0.05, 1e-5, 2e-6]
    r = S.summarize_compare([(False, np.nan, np.nan, np.nan)])
    assert r[0] == 1 and r[1] == 0 and np.isnan(r[2:]).all()


def test_compare_astro():
    a = (np.array([1.0, 5.0]), np.array([2.0, 6.0]), np.array([100.0, 200.0]))
    b = (np.array([5.0005, 1.0]), np.array([6.0, 2.0]), np.array([201.0, 100.0]))
    ok, dp, df = S.compare_astro(a, b)
    assert ok and dp == pytest.approx(5e-4, rel=1e-9) and df == pytest.approx(0.005, rel=1e-12)
    assert S.compare_astro(a, (b[0][:1], b[1][:1], b[2][:1]))[0] is False
    assert S.compare_astro((np.zeros(0),) * 3, (np.zeros(0),) * 3) == (True, 0.0, 0.0)


def test_compare_traces():
    ev = [(0.5, (3, 4), 'enter'), (0.7, (5, 6), 'enter'), (0.7, (5, 6), 'exit')]
    ta = dict(events=ev, taus=np.array([0.5, 0.7, 0.7]), barred=[(0.7, (5, 6))])
    tb = dict(events=[(e[0] * (1 + 1e-9), e[1], e[2]) for e in ev], taus=ta['taus'] * (1 + 1e-9), barred=[(0.7000001, (5, 6))])
    same, dt, bar = S.compare_traces(ta, tb)
    assert same and bar and dt == pytest.approx(1e-9, rel=1e-3)
    tc = dict(tb, events=ev[:2] + [(0.7, (5, 7), 'exit')])
    assert S.compare_traces(ta, tc)[0] is False
    td = dict(tb, barred=[])
    assert S.compare_traces(ta, td)[2] is False
    te = dict(tb, barred=[(0.7, (9, 9))])
    assert S.compare_traces(ta, te)[2] is False                             # the same number of barred points, another point


def test_boot_ratio_is_the_ratio_of_sums_with_a_reproducible_interval():
    num = np.array([2.0, 4.0, 6.0, 8.0])
    den = np.array([1.0, 1.0, 1.0, 1.0])
    r, lo, hi = S.boot_ratio(num, den)
    assert r == 5.0 and lo < 5.0 < hi and (lo, hi) == S.boot_ratio(num, den)[1:]
    r2, lo2, hi2 = S.boot_ratio(3.0 * den, den)
    assert (r2, lo2, hi2) == (3.0, 3.0, 3.0)
    r3 = S.boot_ratio(np.array([1.0, 100.0]), np.array([1.0, 1.0]))[0]
    assert r3 == pytest.approx(50.5)                                          # sum / sum, not the mean of the per-frame ratios of 1 and 100 either way: here equal; a different case:
    assert S.boot_ratio(np.array([10.0, 10.0]), np.array([1.0, 9.0]))[0] == pytest.approx(2.0)      # 20 / 10, the mean of the ratios would be 5.5


def test_fit_slope():
    x = np.array([1, 2, 4, 8, 16, 32], float)
    s, lo, hi = S.fit_slope(x, 3.0 * x ** 2)
    assert s == pytest.approx(2.0, abs=1e-9) and lo == pytest.approx(2.0, abs=1e-9) and hi == pytest.approx(2.0, abs=1e-9)
    s, lo, hi = S.fit_slope(x, x ** 1.0 * np.array([1, 1.2, 0.8, 1.1, 0.9, 1.05]))
    assert lo < s < hi


def test_per_frame_matrix_round_trip():
    sc = C.score([np.array([[100.0, 200.0]]), np.zeros((0, 2))], np.array([(0, 1, 105.0, 195.0, 0.0, 1000.0)], dtype=C.GT_DTYPE), [0, 1])
    m = S.pf_matrix(sc)
    p = S.pf_struct(m)
    for k in ('frame', 'tp', 'fp', 'fn', 'n_det', 'n_gt'):
        assert p[k].tolist() == sc['per_frame'][k].tolist()
    assert np.allclose(p['sse'], sc['per_frame']['sse']) and np.allclose(p['sx'], sc['per_frame']['sx'])


def test_md_table():
    t = S.md_table(['a', 'b'], [['1', '2'], ['3', '4']])
    assert t == '| a | b |\n|---|---|\n| 1 | 2 |\n| 3 | 4 |'


def test_the_tolerances_are_the_plans():
    plan = open(os.path.join(REPO, 'results', 'smlm', 'speed2_plan.md')).read()
    assert (S.TOL_POS_NM, S.TOL_PH_REL, S.TOL_SIG_PX, S.TOL_E) == (0.1, 1e-3, 1e-3, 1e-3)
    assert (S.TOL_ASTRO_POS, S.TOL_ASTRO_FLUX, S.TOL_ASTRO_TAU) == (1e-3, 1e-3, 1e-6)
    for token in ('0.1 nm', 'relative 1e-3 in photons', '1e-3 px in width', '1e-3 in E', '1e-3 px, fluxes within a relative 1e-3', 'relative 1e-6'):
        assert token in plan
    assert (S.N_EQ, S.N_VAR, S.N_TIME, S.REPEATS, S.BOOT) == ({'LD': 150, 'HD': 40}, 12, 24, 3, 1000)
    assert S.ASTRO_CAPS == (25, 50, 100, 200, 400) and S.S1_K == (1, 2, 4, 8, 12, 16, 24) and S.S2_SIDES == (32, 48, 64, 96, 128)
    assert set(S.VARIANTS) == {'hp3_k7', 'greedy', 'nowidth', 'two_rules'}


# ---- the report --------------------------------------------------------------------------------------------------------------------------------

def _raw(cmp_pos=1e-12, equal=True, e_new=70.0, astro_ok=True):
    name = 'MT0.N1.LD'
    raw = {}
    gt = np.array([(0, 1, 105.0, 195.0, 0.0, 1000.0), (1, 1, 305.0, 295.0, 0.0, 900.0)], dtype=C.GT_DTYPE)
    sc = C.score([np.array([[100.0, 200.0]]), np.array([[300.0, 300.0]])], gt, [0, 1])
    sc2 = C.score([np.array([[100.0, 200.0]]), np.array([[300.0, 300.0]] if e_new == 70.0 else [[300.0, 300.0], [10.0, 10.0]])], gt, [0, 1])
    raw['E|%s|cmp' % name] = np.array([2, 2 if equal else 1, cmp_pos, 1e-15, 1e-15])
    raw['E|%s|pf_old' % name] = S.pf_matrix(sc)
    raw['E|%s|pf_new' % name] = S.pf_matrix(sc2)
    for v in S.VARIANTS:
        raw['V|%s|%s|cmp' % (v, name)] = np.array([12, 12, 1e-12, 1e-15, 1e-15])
    cpu = np.array([[0.30, 0.34, 0.32, 0.31], [0.020, 0.022, 0.021, 0.020], [0.09, 0.10, 0.095, 0.09]])
    raw['T|%s|cpu' % name] = np.stack([cpu[0], cpu[2], cpu[1]])         # engines: old, new, no-early (stored order: old, new, no early)
    raw['T|%s|cpu' % name] = np.array([cpu[0], cpu[1], cpu[2]])
    raw['P|%s|cpu' % name] = np.array([0.0025, 0.0026, 0.0025, 0.0024])
    raw['P|%s|version' % name] = np.array('0.11.3')
    raw['ENV|cpu_model'], raw['ENV|cores'] = np.array('Test CPU @ 2.30GHz'), np.array(4)
    raw['ENV|load'] = np.array([[1.0, 8.0, 7.0, 6.0], [31.0, 10.0, 8.0, 7.0], [61.0, 12.0, 9.0, 8.0]])
    raw['TT|%s|cpu' % name] = np.array([[0.6, 0.5, 0.7], [0.03, 0.04, 0.035], [0.2, 0.2, 0.2]])
    raw['TT|%s|cfg' % name] = np.array('hp5h_r25')
    raw['TT|%s|kappa' % name] = np.array(6.0)
    raw['TT|%s|ndet' % name] = np.array([[1, 1], [2, 2], [0, 1]])
    rows = []
    for sweep, xs in ((1, (1, 2, 4)), (2, (32, 48, 64))):
        for x in xs:
            for fr in range(2):
                for e in range(2):
                    k = x if sweep == 1 else int(round(6 * x * x / 4096)) or 1
                    rows.append((sweep, 64 if sweep == 1 else x, k if sweep == 1 else max(k, 1), fr, e, (0.1 if e == 0 else 0.01) * (x if sweep == 1 else (x / 32.0) ** 2) * (1 + 0.01 * fr)))
    raw['S|rows'] = np.array(rows, float)
    for cap_n in S.ASTRO_CAPS:
        k = 'A|x110_y630|%d|' % cap_n
        raw[k + 'cpu'], raw[k + 'wall'], raw[k + 'n'] = np.array([cap_n * 0.05, cap_n * 0.01]), np.array([1.0, 1.0]), np.array([cap_n, cap_n])
        raw[k + 'trace'] = np.array([1.0, 1e-12, 1.0])
        raw[k + 'kappas'] = np.array([2.0, 3.0])
        raw[k + 'cmp'] = np.array([[1.0, 1e-9, 1e-9], [1.0, 1e-9 if astro_ok else 0.5, 1e-9]])
        raw[k + 'ndet'] = np.array([[5, 5], [4, 4]])
    for xy in S.ASTRO_PATCHES:
        for cap_n in (200, 400):
            k = 'AC|x%d_y%d|%d|' % (xy + (cap_n,))
            raw[k + 'cpu'], raw[k + 'trace'], raw[k + 'cmp'] = np.array([1.0, 1.0]), np.array([0.0, np.nan, 1.0]), np.array([[1.0, 2.5, 0.25]])
            if ('A|x%d_y%d|%d|' % (xy + (cap_n,))) + 'cpu' not in raw:
                k2 = 'A|x%d_y%d|%d|' % (xy + (cap_n,))
                raw[k2 + 'cpu'], raw[k2 + 'wall'], raw[k2 + 'n'] = np.array([1.0, 0.5]), np.array([1.0, 1.0]), np.array([cap_n, cap_n])
                raw[k2 + 'trace'], raw[k2 + 'kappas'], raw[k2 + 'cmp'], raw[k2 + 'ndet'] = np.array([1.0, 1e-12, 1.0]), np.array([2.0]), np.array([[1.0, 1e-9, 1e-9]]), np.array([[5, 5]])
        blk = np.zeros((2, 25, 4))
        blk[0, :, :] = [10, 8, 9, 8]
        blk[1, :, :] = [10, 8, 9, 7]
        for e, scale in (('old', 1.0), ('new', 1.0), ('old_tight', 0.5)):
            raw['AX|x%d_y%d|%s|kappas' % (xy + (e,))] = np.array([2.0, 0.0])
            raw['AX|x%d_y%d|%s|blk' % (xy + (e,))] = blk * (np.array([1.0, 0.5, 1.0, 0.5]) if e == 'old_tight' else 1.0)
            raw['AX|x%d_y%d|%s|n' % (xy + (e,))] = np.array({'old': 1424, 'new': 1420, 'old_tight': 1439}[e])
    k = 'AF|x110_y630|'
    raw[k + 'cpu'], raw[k + 'n'], raw[k + 'trace'], raw[k + 'cmp'] = np.array([1200.0, 60.0]), np.array([1300, 1300]), np.array([1.0, 1e-12, 1.0]), np.array([[1.0, 1e-9, 1e-9]])
    return raw


def test_the_report_states_pass_and_the_ratios():
    md = S.report(_raw())
    assert md == S.report(_raw())                                            # deterministic
    assert '| MT0.N1.LD | 2 | 2 / 2 |' in md and md.count('| pass |') >= 2
    assert 'FAIL' not in md
    assert '15.3 [' in md                                                      # old / new = (0.30 + 0.34 + 0.32 + 0.31) / (0.020 + 0.022 + 0.021 + 0.020)
    assert 'Complete path of patch (110, 630): 1300 entries' in md and 'old / new = 20.0' in md
    assert 'Exponent of time in K:' in md and 'Exponent of time in the number of pixels:' in md
    assert 'Picasso 0.11.3' in md
    assert 'Platform: Test CPU @ 2.30GHz, 4 cores; 1-minute load average sampled every 30 s during the timing parts (3 samples): minimum 8.0, median 10.0, maximum 12.0.' in md
    assert '| MT0.N1.LD | hp5h_r25 | 6 | 3 | 0.6000 | 0.2000 | 0.0350 | 17.1 [' in md and '| 2 / 3 |' in md                  # E-tuned table: old / new = 1.8 / 0.105; one frame with unequal counts


def test_the_post_hoc_tables_are_reported():
    md = S.report(_raw())
    assert 'POST HOC control' in md and '| (110, 630) | 200 | NO | nan | 2.5e+00 | 2.5e-01 |' in md
    assert 'POST HOC (not in the plan): equivalence of the COMPLETE paths' in md
    assert '| 2 | 0.8421 | 0.8421 | 0.4211 | +0.0000 [+0.0000, +0.0000] | -0.4211 [-0.4211, -0.4211] |' in md          # completeness 8/10 and purity 8/9; the tight run halves the matches
    assert '| BIC | 0.7887 | 0.7887 | 0.3944 |' in md                                                                  # threshold 0: 7 matched of 9 detections


def test_the_report_flags_every_kind_of_failure():
    assert 'FAIL' in S.report(_raw(cmp_pos=0.2))                              # position beyond 0.1 nm
    assert 'FAIL' in S.report(_raw(equal=False))                              # a frame with a different count
    assert 'FAIL' in S.report(_raw(e_new=71.0))                               # an extra detection changes E
    assert 'FAIL' in S.report(_raw(astro_ok=False))                           # an astronomy flux/position beyond tolerance
    raw = _raw()
    raw['V|greedy|MT0.N1.LD|cmp'] = np.array([12, 12, 1e-12, 0.5, 1e-15])
    assert 'FAIL' in S.report(raw)


@pytest.mark.parametrize('what', ['width', 'events', 'barred', 'time', 'counts', 'flux', 'position'])
def test_each_single_criterion_alone_fails_the_verdict(what):
    raw = _raw()
    k = 'A|x110_y630|25|'
    if what == 'width':
        raw['E|MT0.N1.LD|cmp'] = np.array([2, 2, 1e-12, 1e-15, 0.5])
    elif what == 'events':
        raw[k + 'trace'] = np.array([0.0, 1e-12, 1.0])
    elif what == 'barred':
        raw[k + 'trace'] = np.array([1.0, 1e-12, 0.0])
    elif what == 'time':
        raw[k + 'trace'] = np.array([1.0, 1e-3, 1.0])
    elif what == 'counts':
        raw[k + 'cmp'] = np.array([[0.0, 1e-9, 1e-9], [1.0, 1e-9, 1e-9]])
    elif what == 'flux':
        raw[k + 'cmp'] = np.array([[1.0, 1e-9, 0.5], [1.0, 1e-9, 1e-9]])
    else:
        raw[k + 'cmp'] = np.array([[1.0, 0.5, 1e-9], [1.0, 1e-9, 1e-9]])
    md = S.report(raw)
    assert md.count('FAIL') == 1
    ok = S.report(_raw())
    assert 'FAIL' not in ok


def test_the_report_without_picasso_or_scaling_still_runs():
    raw = {k: v for k, v in _raw().items() if not k.startswith('P|') and not k.startswith('S|') and not k.startswith('AF|')}
    md = S.report(raw)
    assert 'Pooled over the four sets' not in md and 'Exponent' not in md and 'Complete path' not in md


def test_the_main_requires_a_command():
    r = subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'exp_speed2.py')], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='1'))
    assert r.returncode == 2 and 'give --run or --from-raw' in r.stderr


def test_the_stored_report_regenerates_byte_for_byte():
    raw_path = os.path.join(REPO, 'results', 'smlm', 'speed2_raw.npz')
    md_path = os.path.join(REPO, 'results', 'smlm', 'speed2.md')
    if not (os.path.exists(raw_path) and os.path.exists(md_path)):
        pytest.skip('the report has not been generated yet')
    assert S.report(S.load_raw(raw_path)) == open(md_path).read()


# ---- the frame wrapper and the astronomy engine switch --------------------------------------------------------------------------------------

def _frame_and_worker(rules, **extra):
    import exp_smlm_iss as ISS
    from src.core import imaging as IM
    bg = np.full((64, 64), 100.0)
    P = IM.PixelPSF(64, 64, ISS.SIGMA0, parts=ISS.get_cfg('hp5h')['parts'])
    rng = np.random.default_rng(11)
    frame = IM.simulate_emccd(P, np.array([[20.3, 31.2], [44.6, 17.9]]), np.array([2500.0, 1800.0]), bg, rng)
    opts = dict(rules=rules, tol=1e-8, max_nfev=200, **extra)
    ISS._init('hp5h', bg, 0.9, opts)
    return ISS, frame


@pytest.mark.parametrize('rules,extra', [([('snr', 5.0)], {}), ([('snr', 5.0), ('bic', 4)], {}), ([('snr', 7.0)], dict(order='greedy', max_entries=6))], ids=['snr', 'two_rules', 'greedy'])
def test_the_fast_frame_fit_returns_the_same_detections(rules, extra):
    import smlm_fast as SF
    ISS, frame = _frame_and_worker(rules, **extra)
    a = ISS.run_frame((0, frame))
    b = SF.run_frame_fast((0, frame))
    assert len(a[0]) == len(b[0]) == len(rules)
    for da, db in zip(a[0], b[0]):
        ok, dpos, dph, dsig = S.compare_dets(da, db)
        assert ok and dpos < 1e-3 and dph < 1e-6 and dsig < 1e-6
    assert len(a[0][0]) >= 2
    if len(rules) == 1:                                                       # the early exit: a prefix of the full tail
        n = len(b[1])
        assert n <= len(a[1]) and np.allclose(a[1][:n], b[1], rtol=1e-6) and np.allclose(a[2][:n], b[2], rtol=1e-5)
    else:
        assert len(b[1]) == len(a[1])                                         # more than one rule: the whole tail is kept


def test_the_early_exit_can_be_switched_off_in_the_frame_wrapper():
    import smlm_fast as SF
    ISS, frame = _frame_and_worker([('snr', 5.0)])
    on, off = SF.run_frame_fast((0, frame)), SF.run_frame_fast((0, frame), early=False)
    assert len(on[1]) < len(off[1])
    assert S.compare_dets(on[0][0], off[0][0])[0]


def test_the_astronomy_engine_switch_restores_the_path_function_and_the_cap(monkeypatch):
    import astro_fast as AF
    import exp_astro_m2 as X
    from src.core import imaging as IM
    from src.core import imaging_guard as G
    original, cap0 = IM.sliding_path, X.MAX_ENTRIES
    with AF.engine('old'):
        assert IM.sliding_path is AF._path_old
    with AF.engine('new'):
        assert IM.sliding_path is AF._path_new
    assert IM.sliding_path is original
    with pytest.raises(RuntimeError):
        with AF.engine('new'):
            raise RuntimeError
    assert IM.sliding_path is original
    with AF.cap(7):
        assert X.MAX_ENTRIES == 7
    assert X.MAX_ENTRIES == cap0
    psf = IM.PixelPSF(24, 24, 1.3, parts=((1.0, 1.0),), hp=3.0)
    plain = IM.PixelPSF(24, 24, 1.3)
    y = psf.highpass(plain.synth(np.array([[10.2, 11.1]]), np.array([4000.0])) + 0.5 * np.random.default_rng(0).standard_normal((24, 24)))
    L = IM.Lattice(psf, h=0.5, margin=1.0)
    kw = dict(max_entries=4, tol=1e-8, passes=1, rtrunc=4.0)
    sw = np.full((24, 24), 2.0)
    with AF.engine('old'):
        r0 = IM.sliding_path(psf, y, sw, L, **kw)
    t0 = dict(AF.TRACE)
    with AF.engine('new'):
        r1 = IM.sliding_path(psf, y, sw, L, **kw)
    assert [(e[1], e[2]) for e in r0['events']] == [(e[1], e[2]) for e in r1['events']] == [(e[1], e[2]) for e in AF.TRACE['events']]
    assert t0['barred'] == [] and AF.TRACE['barred'] == [] and np.allclose(AF.TRACE['taus'], r1['taus']) and np.allclose(t0['taus'], r0['taus'], rtol=1e-8)
    assert S.compare_traces(t0, dict(AF.TRACE))[0]


# ---- the cost bands, the decomposition of the speed-up, the guard statement and the post hoc F pattern (referee A+B) ---------------------------------------------------------

def test_span_prints_the_smallest_and_the_largest_value_and_one_number_when_they_print_alike():
    assert S.span([5.4, 97.6, 36.0]) == '5-98'
    assert S.span([36.11, 36.14], '%.1f') == '36.1'
    assert S.span([4.8, 3.6, 4.1], '%.1f') == '3.6-4.8'


def test_pooled_boot_pools_the_means_of_the_sets_and_resamples_frames_within_each_set():
    cpus = [np.array([[1.0, 1.0], [0.1, 0.1], [0.5, 0.5]]), np.array([[3.0, 3.0], [0.3, 0.3], [1.5, 1.5]])]
    pics = [np.array([0.01, 0.01]), np.array([0.03, 0.03])]
    r = S.pooled_boot(cpus, pics, B=50)
    assert r['old/new'] == pytest.approx((10.0, 10.0, 10.0)) and r['new/pic'] == pytest.approx((10.0, 10.0, 10.0))
    assert r['nonearly/pic'] == pytest.approx((50.0, 50.0, 50.0)) and r['old/pic'] == pytest.approx((100.0, 100.0, 100.0))        # the ratio of the means of the sets' means: (1 + 3) / (0.01 + 0.03)
    cpus2 = [np.array([[1.0, 3.0], [0.1, 0.1], [0.5, 0.5]]), np.array([[3.0, 3.0], [0.3, 0.3], [1.5, 1.5]])]
    r2 = S.pooled_boot(cpus2, pics, B=400)
    assert r2['old/new'] == pytest.approx((12.5, 10.0, 15.0))                  # only the first set varies: its mean of two resampled frames is 1, 2 or 3, so the pooled old is 2, 2.5 or 3 over a new of 0.2
    assert r2['new/pic'] == pytest.approx((10.0, 10.0, 10.0))
    assert S.pooled_boot(cpus2, pics, B=400) == r2                              # seeded


def test_the_report_states_the_three_cost_bands_and_the_decomposition_of_the_speed_up():
    raw = _raw()
    md = S.report(raw)
    cpu, pic = raw['T|MT0.N1.LD|cpu'], raw['P|MT0.N1.LD|cpu']
    pb = S.pooled_boot([cpu], [pic])
    assert '**Three cost bands against Picasso**' in md
    assert 'the new evaluation with its early exit at one fixed threshold is %.0fx slower than Picasso per set (pooled %.0f [%.0f, %.0f])' % ((cpu[1].sum() / pic.sum(),) + pb['new/pic']) in md
    assert 'without the early exit (the whole path kept) is %.0fx per set (pooled %.0f [%.0f, %.0f])' % ((cpu[2].sum() / pic.sum(),) + pb['nonearly/pic']) in md
    assert 'the original is %.0fx per set (pooled %.0f [%.0f, %.0f])' % ((cpu[0].sum() / pic.sum(),) + pb['old/pic']) in md
    assert 'The pooled old / new = %.1f [%.1f, %.1f] mixes two effects' % pb['old/new'] in md
    assert 'exists only for a single fixed `snr` threshold' in md and 'this holds for every set' in md


def test_a_set_whose_interval_does_not_reach_ten_is_named_as_missing_the_plans_rule():
    raw = _raw()
    raw['T|MT0.N1.LD|cpu'] = raw['T|MT0.N1.LD|cpu'].copy()
    raw['T|MT0.N1.LD|cpu'][1] *= 3.0                                           # the new engine three times slower: old / new about 5
    md = S.report(raw)
    assert 'this does not hold for MT0.N1.LD (lower end' in md and "does not reach the plan's own 10x rule" in md and 'The pooled lower end is' in md and '(below 10)' in md


def test_the_report_says_the_cycle_guard_never_ran_on_real_data():
    md = S.report(_raw())
    assert '**The limit-cycle guard never triggered on real data.**' in md and '`cycle_guard` off' in md and 'planted cycles only' in md
    assert 'is vacuous on real data' in md


def test_the_post_hoc_f_differences_are_counted_by_sign_and_by_interval_and_the_entry_counts_are_printed():
    md = S.report(_raw())
    assert 'The 2 new-minus-old F differences are not all negative (' in md and ' of 2 intervals exclude 0.' in md
    assert 'The tight-tolerance original moves F in ' in md and 'holds in magnitude, not in the sign pattern' in md
    assert '(110, 630): 1424 (old) / 1420 (new) / 1439 (tight)' in md and 'Entry counts of the complete paths can differ between the engines' in md
    raw = _raw()
    for xy in S.ASTRO_PATCHES:
        raw['AX|x%d_y%d|new|blk' % xy] = raw['AX|x%d_y%d|new|blk' % xy] * np.array([1.0, 0.5, 1.0, 0.5])      # the new engine's F falls at every threshold
    md = S.report(raw)
    assert 'The 2 new-minus-old F differences are all negative (' in md


def test_a_point_estimate_above_ten_with_a_lower_end_below_ten_is_not_reported_as_ten_fold():
    raw = _raw()
    raw['T|MT0.N1.LD|cpu'] = np.array([[1.0, 1.0, 1.0, 1.0], [0.05, 0.15, 0.10, 0.099], [0.5, 0.5, 0.5, 0.5]])           # old / new = 4 / 0.399 = 10.03, the lower end of the frame bootstrap below 10
    s = S.boot_ratio(raw['T|MT0.N1.LD|cpu'][0], raw['T|MT0.N1.LD|cpu'][1])
    assert s[0] > 10.0 > s[1]
    md = S.report(raw)
    assert 'this does not hold for MT0.N1.LD (lower end %.1f)' % s[1] in md and 'The pooled lower end is %.1f (below 10)' % s[1] in md


def test_the_f_difference_summary_counts_signs_and_intervals_one_threshold_at_a_time():
    raw = _raw()
    for xy in S.ASTRO_PATCHES:
        blk = raw['AX|x%d_y%d|new|blk' % xy].copy()
        blk[0] = blk[0] * np.array([1.0, 0.5, 1.0, 0.5])                           # the new engine is worse at the first threshold only
        raw['AX|x%d_y%d|new|blk' % xy] = blk
    md = S.report(raw)
    assert 'The 2 new-minus-old F differences are not all negative (-0.4211 to +0.0000); 1 of 2 intervals exclude 0. The tight-tolerance original moves F in one direction (-0.4211 to -0.3944; 2 of 2 of its intervals exclude 0).' in md
    raw2 = _raw()
    for xy in S.ASTRO_PATCHES:
        raw2['AX|x%d_y%d|new|blk' % xy] = raw2['AX|x%d_y%d|new|blk' % xy] * np.array([1.0, 2.0, 1.0, 2.0])               # better than the original: intervals above 0
    assert 'The 2 new-minus-old F differences are not all negative (+0.7887 to +0.8421); 2 of 2 intervals exclude 0.' in S.report(raw2)
    assert 'The 2 new-minus-old F differences are not all negative (+0.0000 to +0.0000); 0 of 2 intervals exclude 0.' in S.report(_raw())
