"""Tests of experiments/exp_peakcount.py with planted data (no CPSC archive needed): the split rule, the front end, the classical counters, the path
selection rules, scoring, the bootstrap, the verdict, and the byte-identical regeneration of the report from a stored raw."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
import exp_peakcount as E  # noqa: E402

SCRIPT = E.__file__


def bumps(centres, amps=None, width=0.04, n=500):
    t = np.arange(n) / E.EFS
    amps = np.ones(len(centres)) if amps is None else np.asarray(amps, float)
    return sum(a * np.exp(-(t - c) ** 2 / (2 * width ** 2)) for c, a in zip(centres, amps))


# ---------------------------------------------------------------------------------------------------------------------------------
# registries and split
# ---------------------------------------------------------------------------------------------------------------------------------
def test_registries():
    assert (E.N_RECORDS, E.FS, E.DEC, E.EFS) == (2000, 500, 10, 50)
    assert (E.SPLIT_SEED, E.N_DEV, E.N_TEST) == (20261110, 150, 500)
    assert E.N_CFG == 42 == len(E.FP_PROM) * len(E.FP_DIST) == len(E.REL_Q) * len(E.FP_DIST) == len(E.MF_SIG) * len(E.FP_PROM) * len(E.MF_DIST)
    assert len(E.path_cfgs()) == E.N_CFG == len(E.WIDTHS) * (len(E.BIC_C) + len(E.DROP_D))
    for f in E.FAMS:
        assert len(E.all_cfgs(f)) == E.N_CFG
    assert E.WIDTHS == (0.03, 0.04, 0.05) and E.NI_MARGIN == 0.10 and E.N_BOOT == 10000 and E.BOOT_SEED == 20261111
    assert E.TOL_S == 0.10 and E.MAXM == 44 and E.GRID_STEP == 0.1
    assert E.CLASSICAL == ('fp', 'fp_rel', 'mf', 'step') and E.FRONT_ENDS == ('sq', 'abs')
    assert (E.BAND, E.WIN_S, E.P_NORM) == ((5.0, 18.0), 0.12, 95.0)
    assert E.FP_PROM == (0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35) and E.FP_DIST == (0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
    assert E.REL_Q == (0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60) and (E.REL_LOW, E.REL_PCT) == (0.02, 75.0)
    assert E.MF_SIG == (0.03, 0.05, 0.07) and E.MF_DIST == (0.20, 0.30)
    assert E.BIC_C == (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0) and E.DROP_D == (0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08)
    assert (E.CAP_PROM, E.CAP_DIST, E.CAP_MULT, E.CAP_ADD) == (0.05, 0.14, 1.4, 4) and E.PATH_FAMS == ('step', 'omp', 'iss')
    assert E.FAMS == ('fp', 'fp_rel', 'mf', 'step', 'omp', 'iss')


def test_split_rule():
    dev, test = E.split_ids()
    assert len(dev) == 150 and len(test) == 500
    assert len(set(dev) & set(test)) == 0
    assert dev.min() >= 1 and test.max() <= 2000
    assert np.all(np.diff(dev) > 0) and np.all(np.diff(test) > 0)
    p = np.random.default_rng(20261110).permutation(np.arange(1, 2001))
    assert set(dev) == set(p[:150]) and set(test) == set(p[150:650])
    d2, t2 = E.split_ids()
    assert np.array_equal(dev, d2) and np.array_equal(test, t2)


# ---------------------------------------------------------------------------------------------------------------------------------
# front end
# ---------------------------------------------------------------------------------------------------------------------------------
def synthetic_ecg(beats, fs=500, dur=10.0, noise=0.02, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(int(fs * dur)) / fs
    x = np.zeros_like(t)
    for b in beats:
        x += np.exp(-(t - b) ** 2 / (2 * 0.012 ** 2)) - 0.2 * np.exp(-(t - b - 0.02) ** 2 / (2 * 0.015 ** 2))
        x += 0.2 * np.exp(-(t - b - 0.25) ** 2 / (2 * 0.04 ** 2))        # a T wave
    return x + noise * rng.standard_normal(len(t))


@pytest.mark.parametrize('fe', E.FRONT_ENDS)
def test_envelope_shape_scale_and_count(fe):
    beats = np.arange(0.7, 9.6, 0.8)
    x = synthetic_ecg(beats)
    e = E.envelope(x, fe)
    assert e.shape == (500,) and np.all(e >= 0)
    assert np.isclose(np.percentile(e, 95), 1.0)
    pos = E.fp_positions(e, 0.2, 0.3)
    assert len(pos) == len(beats)
    assert np.all(np.abs(pos - beats) < 0.12)


def test_envelope_flat_signal_is_finite():
    e = E.envelope(np.zeros(5000), 'sq')
    assert e.shape == (500,) and np.all(np.isfinite(e)) and e.max() == 0


def test_envelope_front_ends_differ():
    x = synthetic_ecg(np.arange(0.7, 9.6, 0.8))
    assert not np.allclose(E.envelope(x, 'sq'), E.envelope(x, 'abs'))


# ---------------------------------------------------------------------------------------------------------------------------------
# classical counters
# ---------------------------------------------------------------------------------------------------------------------------------
def test_fp_distance_and_prominence():
    e = bumps([1.0, 1.2, 5.0], [1.0, 0.8, 0.15], width=0.03)
    assert len(E.fp_positions(e, 0.05, 0.15)) == 3
    assert len(E.fp_positions(e, 0.05, 0.40)) == 2          # the bumps 0.2 s apart merge when the minimum distance is 0.4 s
    assert len(E.fp_positions(e, 0.30, 0.15)) == 2          # the small bump is below the prominence
    assert np.allclose(E.fp_positions(e, 0.05, 0.15), [1.0, 1.2, 5.0], atol=0.02)


def test_fp_rel_scales_with_the_record():
    e1 = bumps([1.0, 3.0, 5.0, 7.0], [1.0, 1.0, 1.0, 0.2], width=0.03)
    e2 = 10.0 * e1                                           # the same record on another scale: the relative rule gives the same answer
    for q in (0.1, 0.4, 0.9):
        assert len(E.fp_rel_positions(e1, q, 0.2)) == len(E.fp_rel_positions(e2, q, 0.2))
    assert len(E.fp_rel_positions(e1, 0.1, 0.2)) == 4
    assert len(E.fp_rel_positions(e1, 0.5, 0.2)) == 3
    assert len(E.fp_rel_positions(np.zeros(500), 0.1, 0.2)) == 0


def test_mf_counts_and_ignores_a_spike():
    e = bumps([1.0, 3.0, 5.0, 7.0], width=0.04)
    e[200] += 0.8                                            # a one-sample spike: the matched filter smooths it below the bumps
    assert len(E.mf_positions(e, 0.05, 0.4, 0.2)) == 4
    assert len(E.fp_positions(e, 0.3, 0.2)) == 5             # find_peaks counts it
    assert len(E.mf_positions(np.zeros(500), 0.05, 0.3, 0.2)) == 0


# ---------------------------------------------------------------------------------------------------------------------------------
# paths and selection
# ---------------------------------------------------------------------------------------------------------------------------------
def test_select_block_bic_and_drop():
    rss = np.array([10.0, 5.0, 2.0, 1.9, 1.89, np.nan, np.nan])
    cnt = np.array([0.0, 1.0, 2.0, 3.0, 4.0, np.nan, np.nan])
    N = 500
    assert E.select_block(rss, cnt, N, 'bic', 1.0) == 3      # the usual BIC accepts a gain of 25.6 > 12.4 but not 2.6
    assert E.select_block(rss, cnt, N, 'bic', 8.0) == 2      # a multiplier of 8 stops where the drops become small
    assert E.select_block(rss, cnt, N, 'bic', 1e6) == 0
    assert E.select_block(rss, cnt, N, 'drop', 0.05) == 2    # the drop 2.0 -> 1.9 is 0.01 r0 < 0.05 r0
    assert E.select_block(rss, cnt, N, 'drop', 0.5) == 1     # 5 -> 2 is 0.3 r0 < 0.5 r0
    assert E.select_block(rss, cnt, N, 'drop', 0.0001) == 4  # nothing below the threshold: the last block
    assert E.select_block(rss, cnt, N, 'drop', 0.9) == 0


def test_select_block_uses_unit_counts_not_indices():
    rss = np.array([10.0, 2.0, 1.0, np.nan])
    cnt = np.array([0.0, 1.0, 1.0, np.nan])                  # block 2 has the same count as block 1 (a unit left): the penalty is by count
    assert E.select_block(rss, cnt, 500, 'bic', 1.0) == 2


@pytest.mark.parametrize('fam', E.PATH_FAMS)
def test_run_path_recovers_planted_bumps(fam):
    centres = [1.5, 3.2, 5.0, 7.7]
    e = bumps(centres, [1.0, 0.8, 0.9, 0.7], width=0.04) + 0.01 * np.random.default_rng(5).random(500)       # noiseless data make ln(rss) degenerate
    rss, cnt, cen = E.run_path(e, fam, 0.04, 10)
    assert rss.shape == (E.MAXM + 1,) and cnt.shape == (E.MAXM + 1,) and cen.shape == (E.MAXM + 1, E.MAXM)
    assert rss[0] == pytest.approx(float(e @ e))
    k = E.select_block(rss, cnt, len(e), 'bic', 1.0)
    assert cnt[k] == 4
    assert np.allclose(np.sort(cen[k, :4]), centres, atol=0.05)
    assert np.all(np.isnan(rss[12:]))                         # no block beyond the cap
    assert E.path_positions((rss[None], cnt[None], cen[None]), 0, 'bic', 1.0, len(e)) == pytest.approx(centres, abs=0.05)


def test_path_cap_is_label_free_and_capped():
    e = bumps(np.arange(0.5, 9.5, 0.5), width=0.03)
    assert E.path_cap(e) == min(E.MAXM, int(E.CAP_MULT * 18) + E.CAP_ADD)
    assert E.path_cap(np.zeros(500)) == E.CAP_ADD
    e2 = bumps(np.arange(0.2, 9.8, 0.1), width=0.015)
    assert E.path_cap(e2) <= E.MAXM


# ---------------------------------------------------------------------------------------------------------------------------------
# scoring and statistics
# ---------------------------------------------------------------------------------------------------------------------------------
def test_match_f1_cases():
    assert E.match_f1([], [], 0.1) == (0, 0, 0)
    assert E.match_f1([1.0], [], 0.1) == (0, 1, 0)
    assert E.match_f1([], [1.0], 0.1) == (0, 0, 1)
    assert E.match_f1([1.0, 2.0], [1.05, 2.2], 0.1) == (1, 1, 1)
    assert E.match_f1([1.0, 1.05], [1.02], 0.1) == (1, 1, 0)      # one-to-one: the second estimate is a false positive
    assert E.match_f1([1.0], [1.09], 0.1) == (1, 0, 0)
    assert E.match_f1([0.5], [0.75], 0.25) == (1, 0, 0)            # the tolerance is inclusive (exact in binary floating point)
    assert E.match_f1([1.0], [1.11], 0.1) == (0, 1, 1)
    assert E.micro_f1(0, 0, 0) == 0.0 and E.micro_f1(3, 1, 1) == pytest.approx(0.75)


def test_boot_ci_properties():
    m, lo, hi = E.boot_ci(np.full(40, 0.5))
    assert (m, lo, hi) == (0.5, 0.5, 0.5)
    rng = np.random.default_rng(0)
    d = rng.normal(1.0, 1.0, 400)
    m, lo, hi = E.boot_ci(d)
    assert m == pytest.approx(d.mean()) and lo < m < hi
    assert hi - lo == pytest.approx(2 * 1.96 / 20, rel=0.25)
    assert E.boot_ci(d) == (m, lo, hi)                           # seeded
    m0, lo0, hi0 = E.boot_ci(-d)
    assert lo0 == pytest.approx(-hi) and hi0 == pytest.approx(-lo)
    assert E.boot_ci(d, n_boot=200, seed=1) != E.boot_ci(d, n_boot=200, seed=2)


# ---------------------------------------------------------------------------------------------------------------------------------
# verdict, report, regeneration
# ---------------------------------------------------------------------------------------------------------------------------------
def planted_res(err_iss, err_base, err_omp, base='fp'):
    res = {f: dict(cfg=(0.1, 0.2), err=np.zeros(len(err_iss)), f1=0.9, tp=9, fp=1, fn=1) for f in E.FAMS}
    res['iss']['err'], res[base]['err'], res['omp']['err'] = np.asarray(err_iss, float), np.asarray(err_base, float), np.asarray(err_omp, float)
    return res, dict(best_classical=base)


def test_verdict_go_and_kill_cases():
    rng = np.random.default_rng(1)
    n = 300
    base = rng.choice([-2, -1, 0, 1, 2], n)
    ours = np.where(rng.random(n) < 0.6, 0.0, base)                      # clearly better than the baseline
    res, frz = planted_res(ours, base, ours)
    v = E.verdict(res, frz)
    assert v['go'] and v['win'] and v['ni'] and v['base'] == 'fp' and v['win_ci'][1] > 0
    # equal to the baseline: no win
    res, frz = planted_res(base, base, base)
    v = E.verdict(res, frz)
    assert not v['go'] and not v['win'] and v['ni']
    # better than the baseline but clearly worse than the greedy control: not non-inferior
    res, frz = planted_res(ours, base, np.zeros(n))
    v = E.verdict(res, frz)
    assert v['win'] and not v['ni'] and not v['go']
    # worse than the baseline
    res, frz = planted_res(base, ours, base)
    v = E.verdict(res, frz)
    assert not v['win'] and not v['go']
    # uses the mean ABSOLUTE error: signs do not cancel
    res, frz = planted_res(np.where(np.arange(n) % 2 == 0, 1.0, -1.0), np.zeros(n), np.zeros(n))
    assert E.verdict(res, frz)['win_ci'][0] == pytest.approx(-1.0)


def test_verdict_margin_is_strict_upper_bound():
    n = 200
    zero = np.zeros(n)
    res, frz = planted_res(np.full(n, 0.0), np.full(n, 1.0), np.full(n, 0.0))
    assert E.verdict(res, frz)['ni']
    ours = np.zeros(n)
    ours[:E.NI_MARGIN and int(n * 0.2)] = 1.0                            # 0.2 beats/record worse than omp on average: over the 0.10 margin
    res, frz = planted_res(ours, np.full(n, 1.0), zero)
    assert not E.verdict(res, frz)['ni']


def make_planted_store(n=6, fam_widths=(1,)):
    rng = np.random.default_rng(3)
    st = E.Store.__new__(E.Store)
    st.ids, st.fe, st.paths = np.arange(1, n + 1), 'abs', {}
    ref, Es = [], []
    for j in range(n):
        k = 5 + j % 3
        c = np.sort(rng.uniform(0.5, 9.5, k))
        while np.diff(c).min() < 0.6:
            c = np.sort(rng.uniform(0.5, 9.5, k))
        ref.append(c)
        Es.append(bumps(c, rng.uniform(0.6, 1.0, k), width=0.04) + 0.005 * rng.random(500))
    st.ref, st.E = ref, np.array(Es)
    st.K = np.array([len(r) for r in ref])
    return st


def planted_freeze():
    fam = {f: dict(cfg=[0.1, 0.2], dev_mae=1.0, dev_bias=0.0, n_cfg=42) for f in ('fp', 'fp_rel')}
    fam['mf'] = dict(cfg=[0.05, 0.1, 0.3], dev_mae=1.0, dev_bias=0.0, n_cfg=42)
    for f in E.PATH_FAMS:
        fam[f] = dict(cfg=[1, 'bic', 8.0], dev_mae=1.0, dev_bias=0.0, n_cfg=42)
    return dict(front_end='abs', families=fam, best_classical='fp', shift=0.01)


def test_evaluate_and_report_on_planted_store(tmp_path):
    st = make_planted_store()
    E.compute_paths(st, [(f, 1) for f in E.PATH_FAMS])
    frz = planted_freeze()
    res = E.evaluate(st, frz)
    assert set(res) == set(E.FAMS)
    for f in E.FAMS:
        assert res[f]['err'].shape == (6,)
    assert np.all(res['iss']['err'] == 0) and np.all(res['omp']['err'] == 0)          # clean planted bumps: the path families count exactly
    text, v = E.render_report(st, frz)
    assert '## Verdict' in text and ('**GO**' in text or '**KILL**' in text)
    assert not v['win']                                                                  # nothing to win against a perfect baseline-free tie
    raw = tmp_path / 'peakcount_raw.npz'
    E.save_raw(st, np.zeros(6, int), str(raw))
    json.dump(frz, open(tmp_path / 'peakcount_freeze.json', 'w'))
    st2 = E.load_raw(str(raw))
    assert np.array_equal(st2.K, st.K) and st2.fe == 'abs'
    assert all(np.allclose(a, b, equal_nan=True) for a, b in zip(st2.paths[('iss', 1)], st.paths[('iss', 1)]))
    text2, _ = E.render_report(st2, frz)
    assert text2 == text


def test_from_raw_cli_is_byte_identical(tmp_path):
    st = make_planted_store()
    E.compute_paths(st, [(f, 1) for f in E.PATH_FAMS])
    E.save_raw(st, np.zeros(6, int), str(tmp_path / 'peakcount_raw.npz'))
    json.dump(planted_freeze(), open(tmp_path / 'peakcount_freeze.json', 'w'))
    outs = []
    for k in range(2):
        md = tmp_path / ('m%d.md' % k)
        r = subprocess.run([sys.executable, '-I', SCRIPT, '--from-raw', str(tmp_path / 'peakcount_raw.npz'), '--md', str(md)], capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS='1'))
        assert r.returncode == 0, r.stderr
        outs.append(md.read_bytes())
    assert outs[0] == outs[1] and b'Verdict' in outs[0]


class FakeStore:
    """A stand-in with prescribed count errors per configuration (tune needs only K and counts)."""

    def __init__(self, K, errs):
        self.K, self.errs = np.asarray(K), errs

    def counts(self, fam, cfg):
        return self.K + np.asarray(self.errs[cfg])


def test_tune_picks_lowest_absolute_error_first_on_ties(monkeypatch):
    cfgs = [(0.1, 0.1), (0.2, 0.1), (0.3, 0.1), (0.4, 0.1)]
    monkeypatch.setattr(E, 'all_cfgs', lambda fam: cfgs)
    errs = {cfgs[0]: [2, -2, 2, -2], cfgs[1]: [1, 1, 1, 1], cfgs[2]: [-1, -1, -1, -1], cfgs[3]: [3, 3, 3, 3]}
    cfg, mae, bias = E.tune(FakeStore([10, 11, 12, 13], errs), 'fp')
    assert mae == [2.0, 1.0, 1.0, 3.0] and bias == [0.0, 1.0, -1.0, 3.0]
    assert cfg == cfgs[1]                                                # the abs error decides (not the signed mean), the first of the tied


def test_tune_on_planted_store_is_the_argmin():
    st = make_planted_store()
    cfg, mae, bias = E.tune(st, 'fp')
    assert len(mae) == E.N_CFG and E.classical_cfgs('fp').index(cfg) == int(np.argmin(mae)) and min(mae) <= 0.5


def test_position_shift_recovers_planted_offset():
    st = make_planted_store()
    st.ref = [r - 0.03 for r in st.ref]                                   # annotations 30 ms before the envelope peaks
    cfg = (0.2, 0.3)
    assert E.position_shift(st, cfg) == pytest.approx(0.03, abs=0.01)


def test_position_shift_is_a_median_and_ignores_far_peaks():
    st = make_planted_store(n=1)
    c = np.array([1.0, 3.0, 5.0, 7.0, 9.0])
    st.E = np.array([bumps(c, width=0.04)])
    st.ref = [np.array([1.0 - 0.0, 3.0 - 0.0, 5.0 - 0.19])]              # three matched peaks with offsets 0, 0, 0.19; two estimates have no reference within 0.2 s
    got = E.position_shift(st, (0.2, 0.3))
    assert got == pytest.approx(0.0, abs=0.011)                           # median 0 (the mean would be 0.063); the far estimates (distance 2 s) are excluded
    st.E = np.array([bumps([1.0, 3.0, 4.0, 4.6, 5.2], width=0.04)])
    st.ref = [np.array([1.0 - 0.03, 3.0 - 0.03])]
    assert E.position_shift(st, (0.2, 0.3)) == pytest.approx(0.03, abs=0.011)   # estimates 1.0 to 2.2 s away from the references do not count


def test_frozen_cfg_types():
    frz = planted_freeze()
    assert E.frozen_cfg(frz, 'iss') == (1, 'bic', 8.0)
    assert E.frozen_cfg(frz, 'mf') == (0.05, 0.1, 0.3)
    assert E.cfg_name('iss', (1, 'bic', 8.0)) == 'width 0.04 s, bic 8'
    assert E.cfg_name('fp', (0.12, 0.3)) == '0.12, 0.3'


# ---------------------------------------------------------------------------------------------------------------------------------
# the stages, end to end on synthetic records (the archive is replaced by a generator; the split sizes are shrunk)
# ---------------------------------------------------------------------------------------------------------------------------------
def fake_load(i, root=None):
    gap = 0.6 + 0.07 * (int(i) % 4)
    beats = np.arange(0.4, 9.7, gap)
    return synthetic_ecg(beats, seed=int(i)), np.round(beats * E.FS).astype(int)


@pytest.fixture()
def tiny(monkeypatch, tmp_path):
    monkeypatch.setattr(E, 'OUT', str(tmp_path))
    monkeypatch.setattr(E, 'load_record', fake_load)
    monkeypatch.setattr(E, 'N_DEV', 5)
    monkeypatch.setattr(E, 'N_TEST', 4)
    return tmp_path


def test_stage_fe_writes_the_choice(tiny):
    out = E.stage_fe()
    z = json.load(open(tiny / 'peakcount_fe.json'))
    dev, _ = E.split_ids()
    assert z['front_end'] in E.FRONT_ENDS and z == json.loads(json.dumps(out))
    assert z['dev_ids'] == [int(i) for i in dev] and len(dev) == 5
    best = min(E.FRONT_ENDS, key=lambda f: z['dev_fp'][f]['mae'])
    assert z['front_end'] == best
    for f in E.FRONT_ENDS:
        st = E.Store(dev, f)
        assert z['dev_fp'][f]['mae'] == pytest.approx(min(E.tune(st, 'fp')[1]))
        assert z['dev_fp'][f]['cfg'] == list(E.tune(st, 'fp')[0])


def test_dev_freeze_test_pipeline(tiny):
    E.stage_fe()
    E.stage_dev()
    raw = np.load(tiny / 'peakcount_dev_raw.npz')
    assert str(raw['fe']) == json.load(open(tiny / 'peakcount_fe.json'))['front_end'] and len(raw['ids']) == 5
    for f in E.PATH_FAMS:
        for wi in range(len(E.WIDTHS)):
            assert raw['%s_%d_rss' % (f, wi)].shape == (5, E.MAXM + 1)
    E.stage_freeze()
    frz = json.load(open(tiny / 'peakcount_freeze.json'))
    assert set(frz['families']) == set(E.FAMS) and frz['best_classical'] in E.CLASSICAL
    st = E.load_raw(str(tiny / 'peakcount_dev_raw.npz'))
    for f in E.FAMS:
        cfg, mae, _ = E.tune(st, f)
        assert frz['families'][f]['dev_mae'] == pytest.approx(min(mae)) and frz['families'][f]['n_cfg'] == E.N_CFG
    assert frz['best_classical'] == min(E.CLASSICAL, key=lambda f: frz['families'][f]['dev_mae'])
    assert 'Best tuned classical counter' in (tiny / 'peakcount_dev.md').read_text()
    E.stage_test()
    z = np.load(tiny / 'peakcount_raw.npz')
    _, test = E.split_ids()
    assert np.array_equal(z['ids'], test) and len(test) == 4
    need = {(f, int(frz['families'][f]['cfg'][0])) for f in E.PATH_FAMS}
    have = {(f, wi) for f in E.PATH_FAMS for wi in range(len(E.WIDTHS)) if '%s_%d_rss' % (f, wi) in z.files}
    assert have == need                                                   # only the frozen widths are run on the test records
    v = E.from_raw()
    assert (tiny / 'peakcount.md').exists() and isinstance(v['go'], bool)
    first = (tiny / 'peakcount.md').read_bytes()
    E.from_raw()
    assert (tiny / 'peakcount.md').read_bytes() == first


# ---------------------------------------------------------------------------------------------------------------------------------
# pinned behaviour (front end values, edge cases of the counters, selection by count, dispatch, scoring, bootstrap, verdict, report)
# ---------------------------------------------------------------------------------------------------------------------------------
def test_envelope_pinned_values():
    x = synthetic_ecg(np.arange(0.7, 9.6, 0.8))
    e = E.envelope(x, 'sq')
    assert e.sum() == pytest.approx(77.415196, rel=1e-4) and int(e.argmax()) == 435 and e.max() == pytest.approx(1.066822, rel=1e-4)
    assert [e[i] for i in (35, 36, 75, 200, 499)] == pytest.approx([1.038247, 0.995255, 1.039079, 0.050552, 0.000302], rel=1e-3, abs=1e-5)
    e = E.envelope(x, 'abs')
    assert e.sum() == pytest.approx(101.329534, rel=1e-4) and int(e.argmax()) == 435 and e.max() == pytest.approx(1.06384, rel=1e-4)
    assert [e[i] for i in (35, 36, 75, 200, 499)] == pytest.approx([1.047473, 0.994679, 1.048975, 0.166398, 0.014763], rel=1e-3, abs=1e-5)


def spikes(positions, heights):
    e = np.zeros(500)
    for p, h in zip(positions, heights):
        e[p] = h
    return e


def test_minimum_distance_is_rounded_to_samples():
    e = spikes([100, 107], [1.0, 0.8])                                    # 7 samples apart; 0.15 s = 7.5 samples rounds to 8 -> only the taller survives
    assert len(E.fp_positions(e, 0.05, 0.15)) == 1
    assert len(E.fp_rel_positions(e, 0.1, 0.15)) == 1
    assert len(E.fp_positions(e, 0.05, 0.14)) == 2                       # 7 samples: both


def test_fp_rel_threshold_is_inclusive_and_relative_to_the_percentile():
    e = spikes([50, 150, 250, 350], [1.0, 1.0, 1.0, 0.5])                 # prominences 1, 1, 1, 0.5: the 75th percentile is 1
    assert len(E.fp_rel_positions(e, 0.5, 0.2)) == 4                      # 0.5 >= 0.5 x 1
    assert len(E.fp_rel_positions(e, 0.51, 0.2)) == 3
    e = spikes([50, 150, 250, 350], [1.0, 0.8, 0.8, 0.8])                 # the percentile (0.85) is not the maximum
    assert len(E.fp_rel_positions(e, 0.9, 0.2)) == 4                      # threshold 0.765 keeps all four; the maximum would give 0.9 and keep one


def test_mf_edge_cases():
    e = bumps([3.0, 3.2, 7.0], width=0.05)
    assert len(E.mf_positions(e, 0.07, 0.12, 0.1)) == 2                   # the 3-sigma template smooths the two close bumps into one
    for k in (1.0, 10.0):                                                 # the scale of the record does not matter (normalised by the 95th percentile)
        assert len(E.mf_positions(k * bumps([1.0, 3.0, 5.0, 7.0], [1.0, 1.0, 1.0, 0.3], width=0.04), 0.05, 0.9, 0.2)) == 3
    e = bumps([1.0, 1.2], width=0.03)
    assert len(E.mf_positions(e, 0.03, 0.05, 0.3)) == 1 and len(E.mf_positions(e, 0.03, 0.05, 0.1)) == 2
    with np.errstate(all='raise'):                                        # an all-zero record: no division by zero
        assert len(E.mf_positions(np.zeros(500), 0.05, 0.3, 0.2)) == 0
        assert len(E.fp_rel_positions(np.zeros(500), 0.1, 0.2)) == 0


def test_classical_dispatch():
    e = bumps([1.0, 3.0, 5.0], [1.0, 1.0, 0.2], width=0.04)
    assert np.array_equal(E.classical_positions('fp', e, (0.1, 0.2)), E.fp_positions(e, 0.1, 0.2))
    e2 = spikes([50, 150, 250, 350], [0.2, 0.2, 0.2, 1.0])                    # relative and absolute prominence rules differ here
    assert np.array_equal(E.classical_positions('fp_rel', e2, (0.4, 0.2)), E.fp_rel_positions(e2, 0.4, 0.2))
    assert len(E.classical_positions('fp_rel', e2, (0.4, 0.2))) == 4 and len(E.classical_positions('fp', e2, (0.4, 0.2))) == 1
    assert np.array_equal(E.classical_positions('mf', e, (0.05, 0.1, 0.2)), E.mf_positions(e, 0.05, 0.1, 0.2))
    assert len(E.fp_positions(e, 0.1, 0.2)) == 3 and len(E.fp_positions(e, 0.5, 0.2)) == 2


def test_bic_penalty_counts_units_not_blocks():
    rss = np.array([10.0, 2.0, 1.5, np.nan])
    cnt = np.array([0.0, 1.0, 1.0, np.nan])
    assert E.select_block(rss, cnt, 500, 'bic', 20.0) == 2                # same count: no extra penalty for the second block (a gain of 144 against 249 per unit)


def test_run_path_dispatch_and_cap_block():
    from src.core import baselines as B
    from src.core import rich_transit_nd as RN
    e = bumps([2.0, 2.12, 5.0, 7.5], [1.0, 0.9, 0.8, 0.7], width=0.04) + 0.01 * np.random.default_rng(7).random(500)     # a case where the three paths differ
    t = np.arange(500) / E.EFS
    K = RN.BumpsND(t[:, None], 0.04)
    grid = np.arange(0.0, 10.0, E.GRID_STEP)
    ref = {'iss': B.iss_path(K, e, grid, max_entries=5), 'omp': B.omp_slide(K, e, grid, max_units=5), 'step': B.stepwise_bic(K, e, grid, max_units=5)}
    got = {}
    for fam in E.PATH_FAMS:
        got[fam] = E.run_path(e, fam, 0.04, 5)[0]
        assert got[fam][:len(ref[fam]['blocks'])] == pytest.approx([b['rss'] for b in ref[fam]['blocks']])
    assert abs(got['iss'][3] - got['omp'][3]) > 0.1 and abs(got['iss'][2] - got['step'][2]) > 0.01


def test_run_path_full_length_and_grid_start():
    c = np.arange(0.5, 9.8, 0.2)[:E.MAXM]
    e = bumps(c, width=0.03) + 0.01 * np.random.default_rng(2).random(500)
    rss, cnt, cen = E.run_path(e, 'step', 0.04, E.MAXM)
    assert not np.isnan(rss[E.MAXM]) and cnt[E.MAXM] == E.MAXM and np.isnan(cen[E.MAXM]).sum() == 0
    e = bumps([0.2, 9.9], width=0.04) + 0.01 * np.random.default_rng(2).random(500)       # the start grid reaches the ends of the record
    rss, cnt, cen = E.run_path(e, 'step', 0.04, 4)
    assert np.allclose(np.sort(cen[2, :2]), [0.2, 9.9], atol=1e-9)


def test_match_f1_assignment_prefers_valid_pairs():
    assert E.match_f1([0.05, 0.31], [0.46, 0.22], 0.1) == (1, 1, 1)


def test_boot_ci_pinned():
    assert E.boot_ci([0, 0, 0, 1.0]) == (0.25, 0.0, 0.75)
    assert E.boot_ci([0.0, 1.0, 2.0, 0.0, 5.0]) == pytest.approx((1.6, 0.2, 3.4))


def test_verdict_uses_the_interval_ends_and_the_planted_baseline(monkeypatch):
    n = 40
    zero = np.zeros(n)
    # a case with mean > 0 and lower end < 0
    err_b = np.zeros(n)
    err_b[:3] = 4.0                                                           # three records where the baseline is worse by 4: mean 0.3, lower end 0
    res, frz = planted_res(zero, err_b, zero)
    v = E.verdict(res, frz)
    assert v['win_ci'][0] == pytest.approx(0.3) and v['win_ci'][1] <= 0 and not v['win'] and not v['go']
    # the upper end decides non-inferiority: mean below the margin, upper end above it
    err_o = np.zeros(n)
    err_o[:2] = 2.0                                                           # ours worse than omp on two records by 2: mean 0.1; with the margin 0.15 the upper end is above it
    monkeypatch.setattr(E, 'NI_MARGIN', 0.15)
    res, frz = planted_res(err_o, zero + 5.0, zero)
    v = E.verdict(res, frz)
    assert v['ni_ci'][0] == pytest.approx(0.1) and v['ni_ci'][2] > 0.15 and not v['ni']
    # the comparison is strict: an upper end exactly at the margin is not enough (binary-exact values, constant differences)
    monkeypatch.setattr(E, 'NI_MARGIN', 0.125)
    n = 64
    res, frz = planted_res(np.full(n, 1.0), np.full(n, 5.0), np.full(n, 0.875))
    v = E.verdict(res, frz)
    assert v['ni_ci'] == (0.125, 0.125, 0.125) and not v['ni'] and v['win']
    # the baseline named in the freeze is the one compared (not always 'fp')
    res, frz = planted_res(np.full(n, 2.0), np.full(n, 1.0), np.full(n, 2.0), base='mf')
    res['fp']['err'] = np.full(n, 9.0)
    v = E.verdict(res, frz)
    assert v['base'] == 'mf' and v['win_ci'][0] == -1.0 and not v['win']


def test_evaluate_signs_and_shift_on_a_planted_store():
    st = make_planted_store()
    E.compute_paths(st, [(f, 1) for f in E.PATH_FAMS])
    frz = planted_freeze()
    frz['families']['fp']['cfg'] = [5.0, 0.3]                                 # prominence 5: nothing is found, so the error is -K
    res = E.evaluate(st, frz)
    assert np.array_equal(res['fp']['err'], -st.K.astype(float)) and res['fp']['f1'] == 0.0
    assert res['iss']['f1'] > 0.9
    frz['shift'] = 0.5
    assert E.evaluate(st, frz)['iss']['f1'] == 0.0                            # the shared offset is subtracted from the estimates


def test_store_reads_references_in_seconds(monkeypatch):
    monkeypatch.setattr(E, 'load_record', fake_load)
    st = E.Store([3, 8], 'sq')
    for j, i in enumerate([3, 8]):
        assert np.allclose(st.ref[j], fake_load(i)[1] / E.FS) and st.K[j] == len(st.ref[j])
    assert st.E.shape == (2, 500)


def test_render_report_rows_and_bins(monkeypatch):
    st = make_planted_store()
    big = np.linspace(0.4, 9.4, 20)                                           # a record with 20 beats
    st.E = np.vstack([st.E, bumps(big, width=0.04)[None]])
    st.ref = st.ref + [big]
    st.K = np.r_[st.K, len(big)]
    st.ids = np.arange(1, 8)
    E.compute_paths(st, [(f, 1) for f in E.PATH_FAMS])
    frz = planted_freeze()
    text, v = E.render_report(st, frz)
    res = E.evaluate(st, frz)
    assert len(big) == 20 and '| 20-99 | 1 |' in text and '| 0-10 | 6 |' in text
    row = [ln for ln in text.splitlines() if ln.startswith('| ' + E.LABEL['fp'])][0]
    assert row.split('|')[-2].strip() == '%.3f' % res['fp']['f1']
    win_line = [ln for ln in text.splitlines() if ln.startswith('* Win condition')][0]
    assert win_line.endswith('not met.') == (not v['win']) and ('-> met.' in win_line) == v['win']


def test_stage_freeze_records_what_tune_returns(monkeypatch, tmp_path):
    from types import SimpleNamespace
    monkeypatch.setattr(E, 'OUT', str(tmp_path))
    monkeypatch.setattr(E, 'load_raw', lambda path: SimpleNamespace(fe='abs'))
    base = {'fp': 1.0, 'fp_rel': 0.5, 'mf': 2.0, 'step': 0.7, 'omp': 3.0, 'iss': 4.0}

    def fake_tune(st, fam):
        cfg = E.all_cfgs(fam)[1]
        return cfg, [base[fam], base[fam] + 1], [-base[fam] / 10, 0.0]

    monkeypatch.setattr(E, 'tune', fake_tune)
    monkeypatch.setattr(E, 'position_shift', lambda st, cfg: 0.0123)
    E.stage_freeze()
    z = json.load(open(tmp_path / 'peakcount_freeze.json'))
    assert z['best_classical'] == 'fp_rel' and z['shift'] == 0.0123 and z['front_end'] == 'abs'
    for f in E.FAMS:
        assert z['families'][f]['dev_mae'] == base[f] and z['families'][f]['dev_bias'] == pytest.approx(-base[f] / 10)
        assert z['families'][f]['cfg'] == list(E.all_cfgs(f)[1]) and z['families'][f]['n_cfg'] == 2


def test_stage_fe_picks_the_lower_error(monkeypatch, tmp_path):
    from types import SimpleNamespace
    monkeypatch.setattr(E, 'OUT', str(tmp_path))
    monkeypatch.setattr(E, 'Store', lambda ids, fe: SimpleNamespace(fe=fe))
    for lo in ('sq', 'abs'):
        monkeypatch.setattr(E, 'tune', lambda st, fam, lo=lo: ((0.1, 0.2), [1.0 if st.fe == lo else 2.0], []))
        assert E.stage_fe()['front_end'] == lo
        assert json.load(open(tmp_path / 'peakcount_fe.json'))['front_end'] == lo
