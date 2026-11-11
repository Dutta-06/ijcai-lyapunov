"""The protocol-critical code of the astronomy experiment, with planted numbers: the selection of the frozen configuration and of the tuned baselines, the F-measure and the bootstrap
statistics and verdicts, the crowdsource post-filter, golden cases of the scorer, the freeze record and the write-once test part (referee round IS1a: 20 of 24 mutants of this code
survived the earlier tests)."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import astro_common as A      # noqa: E402
import exp_astro_m2 as X      # noqa: E402


def tot(true22, matched22, det22, det22_matched, true_all=0, matched_all=0, det_all=0):
    return np.array([true22, matched22, det22, det22_matched, true_all, matched_all, det_all], float)


# ---- F-measure, completeness and purity ---------------------------------------------------------------------------------------------------------------------------------

def test_f_measure_and_cpf_follow_the_definitions_with_purity_over_detections():
    t = tot(100, 40, 80, 30)
    f, c, p = X.f_measure(t)
    assert abs(c - 0.4) < 1e-15 and abs(p - 30 / 80) < 1e-15                              # purity = matched detections / detections below 22 mag, not over the truth
    assert abs(f - 2 * 0.4 * 0.375 / (0.4 + 0.375)) < 1e-15
    c2, p2, f2 = X.cpf(t)
    assert (c2, p2, f2) == (c, p, f) or (abs(c2 - c) < 1e-15 and abs(p2 - p) < 1e-15 and abs(f2 - f) < 1e-15)
    cs, ps, fs = X.cpf(np.stack([t, tot(10, 5, 5, 5)]))                                  # any leading shape
    assert abs(cs[1] - 0.5) < 1e-15 and abs(ps[1] - 1.0) < 1e-15 and abs(fs[1] - 2 / 3) < 1e-15
    assert X.f_measure(tot(0, 0, 0, 0))[0] == 0.0                                          # empty: no division error


# ---- choose: the frozen configuration -------------------------------------------------------------------------------------------------------------------------------------

def _dev_raw(best, second=None, tie_k=None):
    """Planted development totals: the pooled F of every (hp, rex, kappa) is low except at `best`; the two development patches contribute different halves."""
    raw = {}
    for p in X.DEV:
        for h in X.HPS:
            for r in X.REXS:
                for k in X.KAPPAS:
                    base = tot(100, 20, 100, 20)                                           # F = 0.2
                    raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] = base.copy()
    for cfg, (c, pu) in ((best, (60, 70)), ) + (((second, (60, 69)), ) if second else ()) + (((tie_k, (60, 70)), ) if tie_k else ()):
        for i, p in enumerate(X.DEV):
            raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], cfg[0], cfg[1], cfg[2])] = tot(100, c, 100, pu)
    return raw


def test_choose_takes_the_largest_pooled_f_and_the_larger_kappa_on_a_tie():
    best = (X.HPS[1], X.REXS[0], X.KAPPAS[2])
    hp, rex, kap, f = X.choose(_dev_raw(best, second=(X.HPS[0], X.REXS[1], X.KAPPAS[4])))
    assert (hp, rex, kap) == best
    exp_f = X.f_measure(tot(200, 120, 200, 140))[0]
    assert abs(f - exp_f) < 1e-12
    # a tie (identical totals) goes to the larger kappa
    other = (X.HPS[1], X.REXS[0], X.KAPPAS[5])
    assert X.choose(_dev_raw(best, tie_k=other))[:3] == other
    # the worst configuration is never chosen (a sign error would pick it)
    raw = _dev_raw(best)
    raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (X.DEV[0][0], X.DEV[0][1], X.HPS[0], X.REXS[0], X.KAPPAS[0])] = tot(100, 1, 100, 1)
    assert X.choose(raw)[:3] == best


def test_choose_pools_both_development_patches():
    """Patch 1 alone prefers configuration A, patch 2 alone B; the pooled counts prefer B (so the choice uses both)."""
    A_, B_ = (X.HPS[0], X.REXS[0], X.KAPPAS[0]), (X.HPS[2], X.REXS[1], X.KAPPAS[3])
    raw = _dev_raw(A_)
    for cfg, vals in ((A_, [(100, 60, 100, 60), (100, 10, 100, 10)]), (B_, [(100, 40, 100, 40), (100, 40, 100, 40)])):
        for p, v in zip(X.DEV, vals):
            raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], cfg[0], cfg[1], cfg[2])] = tot(*v)
    # pooled: A -> C 70/200, P 70/200 (F 0.35); B -> 80/200 (F 0.4)
    assert X.choose(raw)[:3] == B_


# ---- choose_baselines ------------------------------------------------------------------------------------------------------------------------------------------------------

def test_choose_baselines_picks_the_largest_pooled_f_per_family_and_the_larger_name_on_a_tie():
    raw = {}
    for p in X.DEV:
        for t in X.PHOT_TAGS:
            raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(100, 10, 100, 10)
        for f in X.CS_FILTERS:
            raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], X.cs_name(f))] = tot(100, 30, 100, 30)
    best_p, best_c = X.PHOT_TAGS[3], X.cs_name(X.CS_FILTERS[7])
    for p in X.DEV:
        raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], best_p)] = tot(100, 20, 100, 20)
        raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], best_c)] = tot(100, 50, 100, 50)
    out = X.choose_baselines(raw)
    assert out == {'photutils': best_p, 'crowdsource': best_c}
    # only one development patch has a variant: it is not eligible (pooled counts need both)
    raw.pop('DB|x%d_y%d|%s|tot' % (X.DEV[1][0], X.DEV[1][1], best_p))
    assert X.choose_baselines(raw)['photutils'] != best_p
    # a tie goes to the larger name
    raw2 = {'DB|x%d_y%d|%s|tot' % (p[0], p[1], t): tot(100, 20, 100, 20) for p in X.DEV for t in X.PHOT_TAGS}
    assert X.choose_baselines(raw2)['photutils'] == max(X.PHOT_TAGS)


# ---- crowdsource post-filter ------------------------------------------------------------------------------------------------------------------------------------------------

def test_cs_dets_applies_every_cut_with_the_stated_inequalities():
    key = 'x110_y630|crowdsource|'
    b = {key + 'flux': np.array([100.0, 100.0, 100.0, 100.0, 100.0, -5.0, 0.0]), key + 'dflux': np.array([10.0, 10.0, 25.0, 10.0, 10.0, 1.0, 1.0]),
         key + 'fracflux': np.array([0.8, 0.5, 0.9, 0.9, 0.9, 0.9, 0.9]), key + 'qf': np.array([0.5, 0.5, 0.5, 0.4, 0.5, 0.9, 0.9]),
         key + 'spread': np.array([0.02, 0.02, 0.02, 0.02, 0.03, 0.0, 0.0]), key + 'x': np.arange(7.0), key + 'y': np.arange(7.0) * 2}
    x, y, f = X.cs_dets(b, (110, 630), (0.0, 0.0, 0.0, np.inf))                           # the default: only non-positive fluxes go
    assert list(x) == [0, 1, 2, 3, 4]
    x, _, _ = X.cs_dets(b, (110, 630), (5.0, 0.0, 0.0, np.inf))                           # flux >= 5 dflux: rows 0, 1, 3, 4 (row 2 has 100 < 5 x 25)
    assert list(x) == [0, 1, 3, 4]
    x, _, _ = X.cs_dets(b, (110, 630), (0.0, 0.8, 0.0, np.inf))                           # fracflux >= 0.8 (inclusive)
    assert list(x) == [0, 2, 3, 4]
    x, _, _ = X.cs_dets(b, (110, 630), (0.0, 0.0, 0.5, np.inf))                           # qf >= 0.5 (inclusive)
    assert list(x) == [0, 1, 2, 4]
    x, _, _ = X.cs_dets(b, (110, 630), (0.0, 0.0, 0.0, 0.02))                             # |spread| <= 0.02 (inclusive)
    assert list(x) == [0, 1, 2, 3]


# ---- bootstrap, intervals and verdicts -------------------------------------------------------------------------------------------------------------------------------------

def test_boot_idx_covers_every_block_and_the_stratified_one_stays_inside_its_patch():
    idx = X.boot_idx(25, B=400, seed=0)
    assert idx.shape == (400, 25) and idx.min() == 0 and idx.max() == 24 and len(np.unique(idx)) == 25
    st = X.boot_idx_strat([25, 25, 10], B=300, seed=1)
    assert st.shape == (300, 60)
    assert st[:, :25].min() >= 0 and st[:, :25].max() <= 24 and st[:, 25:50].min() >= 25 and st[:, 25:50].max() <= 49 and st[:, 50:].min() >= 50 and st[:, 50:].max() <= 59
    assert np.unique(st[:, 50:]).size == 10                                                # every block of the last patch is drawn somewhere


def test_verdict_truth_table():
    assert X.verdict(0.1, (0.05, 0.2)) == 'win'
    assert X.verdict(-0.1, (-0.2, -0.05)) == 'loss'
    assert X.verdict(0.0, (-0.1, 0.1)) == 'tie'
    assert X.verdict(0.1, (-0.05, 0.2)) == 'tie'                                           # an interval that contains zero is a tie whatever the point estimate
    assert X.verdict(-0.1, (-0.2, 0.05)) == 'tie'


def _blocks(rng, comp, pur, n=25):
    t22 = rng.integers(30, 50, n)
    m22 = rng.binomial(t22, comp)
    d22 = rng.integers(30, 50, n)
    dm = rng.binomial(d22, pur)
    return np.c_[t22, m22, d22, dm].astype(float)


def test_boot_diff_has_the_sign_of_a_minus_b_and_paired_intervals_give_the_verdicts():
    rng = np.random.default_rng(3)
    a = _blocks(rng, 0.8, 0.8)
    b = _blocks(rng, 0.4, 0.4)
    idx = X.boot_idx(25, B=500, seed=0)
    dc, dp, df = X.boot_diff(a, b, idx)
    assert dc[0] > 0.3 and dp[0] > 0.3 and df[0] > 0.3
    assert dc[1][0] > 0 and df[1][0] > 0 and X.verdict(df[0], df[1]) == 'win'
    dc2, dp2, df2 = X.boot_diff(b, a, idx)
    assert abs(df2[0] + df[0]) < 1e-12 and df2[1][1] < 0 and X.verdict(df2[0], df2[1]) == 'loss'
    d0 = X.boot_diff(a, a, idx)
    assert all(abs(v[0]) < 1e-15 and v[1][0] <= 0 <= v[1][1] for v in d0) and X.verdict(*d0[2]) == 'tie'


def test_boot_cpf_interval_is_the_95_percent_percentile_interval_of_the_resamples():
    rng = np.random.default_rng(4)
    blk = _blocks(rng, 0.6, 0.7)
    idx = X.boot_idx(25, B=800, seed=0)
    (c, cci), (p, pci), (f, fci) = X.boot_cpf(blk, idx)
    cs, ps, fs = X.cpf(blk[idx].sum(1))
    assert abs(c - blk[:, 1].sum() / blk[:, 0].sum()) < 1e-15
    assert np.allclose(cci, np.percentile(cs, [2.5, 97.5])) and np.allclose(fci, np.percentile(fs, [2.5, 97.5])) and np.allclose(pci, np.percentile(ps, [2.5, 97.5]))
    assert cci[0] < c < cci[1]


# ---- golden cases of the scorer ---------------------------------------------------------------------------------------------------------------------------------------------

def _score(det_xy, det_mag, tr_x, tr_y, tr_mag):
    """Score detections given by (x, y, magnitude) against a planted truth, on a 100 x 100 patch."""
    det_xy = np.atleast_2d(np.asarray(det_xy, float))
    flux = A.mag_to_counts(np.asarray(det_mag, float))
    t = dict(x=np.asarray(tr_x, float), y=np.asarray(tr_y, float), mag=np.asarray(tr_mag, float))
    t['inner'] = (t['x'] >= A.EDGE - 0.5) & (t['x'] <= 100 - 0.5 - A.EDGE) & (t['y'] >= A.EDGE - 0.5) & (t['y'] <= 100 - 0.5 - A.EDGE)
    return A.score(det_xy, flux, t)


def test_score_blocks_assign_a_match_to_the_truth_block_for_the_truth_and_to_the_detection_block_for_the_detection():
    s = _score([[19.7, 50.0]], [20.0], [19.4], [50.0], [20.0])                              # the block border is at x = 19.5 (pixel centres at the integers): the pair straddles it
    b_truth = A.block_of(np.array([[19.4, 50.0]]))[0]
    b_det = A.block_of(np.array([[19.7, 50.0]]))[0]
    assert b_truth != b_det
    blk = s['blk22']
    assert blk[b_truth].tolist() == [1, 1, 0, 0] and blk[b_det].tolist() == [0, 0, 1, 1]


def test_score_detections_are_counted_below_22_by_their_own_magnitude_and_truth_by_the_truth_magnitude():
    s = _score([[50.0, 50.0]], [21.9], [50.0], [50.0], [22.2])                             # |dm| = 0.3: matched; the detection is below 22 mag, the truth is not
    assert (s['true22'], s['matched22'], s['det22'], s['det22_matched']) == (0, 0, 1, 1)
    s = _score([[50.0, 50.0]], [21.8], [50.0, 70.0], [50.0, 70.0], [21.9, 22.3])           # truth at 22.3 is not in the below-22 total (depth 22, not 22.5)
    assert s['true22'] == 1 and s['matched22'] == 1 and s['true_all'] == 2


def test_score_matching_is_one_to_one_of_maximum_size_and_minimum_total_distance():
    det = [[10.45, 50.0], [10.0, 50.0]]                                                    # B first, then A: the identity assignment (B-P, A-Q) is the wrong one
    s = _score(det, [20.0, 20.0], [10.1, 10.5], [50.0, 50.0], [20.0, 20.0])                # A-P 0.1, A-Q 0.5, B-P 0.35, B-Q 0.05: {A-P, B-Q} (0.15) beats {A-Q, B-P} (0.85)
    r, c = A.match(np.array(det), np.array([20.0, 20.0]), np.array([[10.1, 50.0], [10.5, 50.0]]), np.array([20.0, 20.0]))
    assert sorted(zip(r.tolist(), c.tolist())) == [(0, 1), (1, 0)]
    assert s['matched22'] == 2 and s['det22_matched'] == 2


def test_truth_box_and_scored_margin_edges():
    z = dict(px=np.array([96.9, 97.2, 99.4, 99.6, -0.4, 2.0, 1.9]) + 10, py=np.full(7, 50.0) + 20, m606=np.full(7, 20.0))
    t = A.truth_for_patch(z, 10, 20)
    # x in [-0.5, 99.5): 99.6 is outside; scored (inner) iff 2.0 <= x <= 97.0
    assert len(t['x']) == 6 and not (t['x'] > 99.5).any()
    inner = dict(zip(np.round(t['x'], 1).tolist(), t['inner'].tolist()))
    assert inner[96.9] and not inner[97.2] and not inner[99.4] and not inner[-0.4] and inner[2.0] and not inner[1.9]
    zy = dict(px=np.full(7, 50.0) + 10, py=np.array([96.9, 97.2, 99.4, 99.6, -0.4, 2.0, 1.9]) + 20, m606=np.full(7, 20.0))              # the same along y
    ty = A.truth_for_patch(zy, 10, 20)
    assert len(ty['y']) == 6 and not (ty['y'] > 99.5).any()
    iny = dict(zip(np.round(ty['y'], 1).tolist(), ty['inner'].tolist()))
    assert iny[96.9] and not iny[97.2] and not iny[99.4] and not iny[-0.4] and iny[2.0] and not iny[1.9]


def test_score_ignores_nonpositive_flux_and_unscored_edge_detections():
    t = dict(x=np.array([50.0]), y=np.array([50.0]), mag=np.array([20.0]), inner=np.array([True]))
    flux = np.array([A.mag_to_counts(20.0), -3.0, 0.0, A.mag_to_counts(20.0)])
    s = A.score(np.array([[50.0, 50.0], [60.0, 60.0], [70.0, 70.0], [1.0, 50.0]]), flux, t)       # detections: matched, negative flux, zero flux, within 2.5 px of the edge
    assert s['det22'] == 1 and s['det22_matched'] == 1 and s['det_all'] == 1


# ---- the freeze record and the write-once test part -----------------------------------------------------------------------------------------------------------------------

def _fake_env(tmp_path, monkeypatch):
    files = []
    for name, text in (('imaging.txt', 'a'), ('common.txt', 'b'), ('base.txt', 'c'), ('cat.txt', 'd')):
        (tmp_path / name).write_text(text)
        files.append(name)
    (tmp_path / 'truth_hugs_meth1.npz').write_bytes(b'truth')
    monkeypatch.setattr(X, 'REPO', str(tmp_path))
    monkeypatch.setattr(X, 'ROOT', str(tmp_path))
    monkeypatch.setattr(X, 'FREEZE_FILES', tuple(files))
    monkeypatch.setattr(X, 'FREEZE', str(tmp_path / 'freeze.json'))
    best = (X.HPS[1], X.REXS[0], X.KAPPAS[2])
    raw = _dev_raw(best)
    for p in X.DEV:
        for t in X.PHOT_TAGS:
            raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(100, 10, 100, 10)
        for f in X.CS_FILTERS:
            raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], X.cs_name(f))] = tot(100, 30, 100, 30)
    return best, raw


def test_freeze_record_matches_until_something_it_covers_changes(tmp_path, monkeypatch):
    best, raw = _fake_env(tmp_path, monkeypatch)
    with pytest.raises(RuntimeError, match='no freeze record'):
        X.check_freeze(raw)
    rec = X.write_freeze(raw)
    assert rec['config']['kappa'] == best[2] and rec['config']['sigma_b'] == best[0] and rec['tests'] == [list(p) for p in X.TESTS]
    assert rec['tuned_baselines']['photutils'] == max(X.PHOT_TAGS) and set(rec['files']) == set(X.FREEZE_FILES)
    assert X.check_freeze(raw)['frozen_code'] == X.frozen_code_hash()
    # a changed input file
    (tmp_path / 'common.txt').write_text('B')
    with pytest.raises(RuntimeError, match='files'):
        X.check_freeze(raw)
    (tmp_path / 'common.txt').write_text('b')
    X.check_freeze(raw)
    # a changed development result
    raw2 = dict(raw)
    k = 'D|x%d_y%d|hp%g|rex%g|k%g|tot' % (X.DEV[0][0], X.DEV[0][1], X.HPS[0], X.REXS[0], X.KAPPAS[0])
    raw2[k] = raw2[k] + 1.0
    with pytest.raises(RuntimeError, match='dev_keys'):
        X.check_freeze(raw2)
    # a changed constant that decides a test number
    monkeypatch.setattr(X, 'MARGIN', 3.0)
    with pytest.raises(RuntimeError, match='frozen_code'):
        X.check_freeze(raw)


def test_frozen_code_hash_is_deterministic_and_sensitive_to_the_constants(monkeypatch):
    h = X.frozen_code_hash()
    assert h == X.frozen_code_hash()
    monkeypatch.setattr(X, 'DEFAULT', (4.0, 1.2, 5.0))
    assert X.frozen_code_hash() != h
    monkeypatch.undo()
    monkeypatch.setattr(X, 'TESTS', X.TESTS[:2])
    assert X.frozen_code_hash() != h


def test_part_t_refuses_without_a_matching_freeze_runs_each_patch_once_and_stores_the_labelled_rows(tmp_path, monkeypatch):
    best, raw = _fake_env(tmp_path, monkeypatch)
    pxy = (310, 630)
    monkeypatch.setattr(X, 'TESTS', [pxy])
    X.save_part(str(tmp_path / 'raw.npz'), 'DUMMY', {})
    out_path = str(tmp_path / 'raw.npz')
    np.savez_compressed(out_path, **raw)
    with pytest.raises(RuntimeError, match='no freeze record'):
        X.part_T(out_path)
    truth = dict(px=np.array([30.0 + 310, 60.0 + 310]), py=np.array([30.0 + 630, 60.0 + 630]), m606=np.array([20.0, 20.5]))
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': truth)
    monkeypatch.setattr(X, 'load_patch', lambda x0, y0, ds='m2': None)
    calls = []

    def fake_run_path(patch, hp, rex, **kw):
        calls.append((hp, rex))
        def det(k):                                                                         # the catalogue depends on kappa: int(k) + 1 sources
            n = int(k) + 1
            return (30.0 + np.arange(n), 30.0 + np.arange(n), A.mag_to_counts(np.full(n, 20.0)))
        return {k: det(k) for k in list(X.KAPPAS) + [0.0]}, 1.5, 2, dict(bic_index=1, n_blocks=3)

    monkeypatch.setattr(X, 'run_path', fake_run_path)
    sc = A.score(np.array([[30.0, 30.0]]), A.mag_to_counts(np.array([20.0])), A.truth_for_patch(truth, *pxy))
    every = dict({tag: sc for tag in X.PHOT_TAGS}, **{X.cs_name(f): sc for f in X.CS_FILTERS})
    sc_asis = A.score(np.zeros((0, 2)), np.zeros(0), A.truth_for_patch(truth, *pxy))                 # as shipped: nothing matches
    every.update({'daophot': sc, 'daophot_asis': sc_asis})
    monkeypatch.setattr(X, 'baseline_scores', lambda p, t, base=None: every)
    X.write_freeze(raw)                                                                      # after the stand-ins are in place: the freeze hashes the functions it will run
    X.part_T(out_path)
    z = X.load_raw(out_path)
    pre = 'T|x310_y630|'
    assert z[pre + 'cfg'].tolist() == [best[0], best[1], best[2]]
    assert pre + 'ours|k%g|tot' % best[2] in z and pre + 'default|k0|tot' in z and pre + 'ours|bic' in z and z[pre + 'ours|bic'].tolist() == [1, 3, 0]
    for label in ('photutils (default)', 'crowdsource (default)', 'daophot (published)', 'daophot (as shipped)', 'photutils (tuned)', 'crowdsource (tuned)'):
        assert pre + label + '|tot' in z, label
    assert str(z[pre + 'photutils (tuned)|name'][0]) == max(X.PHOT_TAGS)
    assert z[pre + 'ours|det'].shape == (int(best[2]) + 1, 3) and z[pre + 'default|det'].shape == (int(X.DEFAULT[2]) + 1, 3)      # each stored at its own kappa
    assert z[pre + 'daophot (published)|tot'].tolist() != z[pre + 'daophot (as shipped)|tot'].tolist()
    n_calls = len(calls)
    X.part_T(out_path)                                                                       # a patch with a result is not run again
    assert len(calls) == n_calls
    (tmp_path / 'base.txt').write_text('changed')                                            # a changed input: refuse
    with pytest.raises(RuntimeError, match='files'):
        X.part_T(out_path)


# ---- more golden cases (the matching tolerances, the block grid), part DB and part D bookkeeping -------------------------------------------------------------------------------------

def test_match_needs_both_tolerances_and_block_of_has_five_columns():
    det_xy, tr_xy = np.array([[10.0, 10.0]]), np.array([[10.3, 10.0]])
    assert len(A.match(det_xy, np.array([20.0]), tr_xy, np.array([20.2]))[0]) == 1                  # position and magnitude within
    assert len(A.match(det_xy, np.array([20.0]), tr_xy, np.array([20.7]))[0]) == 0                  # position within, magnitude beyond 0.5
    assert len(A.match(det_xy, np.array([20.0]), np.array([[10.7, 10.0]]), np.array([20.1]))[0]) == 0   # magnitude within, position beyond 0.5
    cols = A.block_of(np.array([[0.0, 0.0], [19.4, 0.0], [19.6, 0.0], [79.6, 0.0], [99.0, 0.0], [99.0, 99.0]]))
    assert cols.tolist() == [0, 0, 1, 4, 4, 24]                                                       # 5 x 5 blocks of 20 px, borders at x = 19.5, 39.5, ...


def test_part_db_stores_the_baseline_scores_of_the_development_patches_under_db_keys(tmp_path, monkeypatch):
    truth = dict(px=np.array([30.0 + 110, 30.0 + 210]), py=np.array([30.0 + 630, 30.0 + 630]), m606=np.array([20.0, 20.0]))
    monkeypatch.setattr(X, 'truth_arrays', lambda: truth)
    out_path = str(tmp_path / 'raw.npz')
    sc = {p: A.score(np.array([[30.0, 30.0]]), A.mag_to_counts(np.array([20.0])), A.truth_for_patch(truth, *p)) for p in X.DEV}
    sc[X.DEV[1]] = A.score(np.zeros((0, 2)), np.zeros(0), A.truth_for_patch(truth, *X.DEV[1]))           # the second patch: nothing detected
    monkeypatch.setattr(X, 'baseline_scores', lambda p, t: {'photutils': sc[p], 'crowdsource|snr0|ff0|qf0|spinf': sc[p]})
    X.part_DB(out_path)
    z = X.load_raw(out_path)
    assert z['DB|x110_y630|photutils|tot'][:4].tolist() == [1, 1, 1, 1] and z['DB|x210_y630|photutils|tot'][:4].tolist() == [1, 0, 0, 0]
    assert 'DB|x110_y630|crowdsource|snr0|ff0|qf0|spinf|tot' in z and not any(k.startswith('D|') for k in z)


_SC = A.score(np.array([[30.0, 30.0]]), A.mag_to_counts(np.array([20.0])), dict(x=np.array([30.0]), y=np.array([30.0]), mag=np.array([20.0]), inner=np.array([True])))


def _fake_dev_job(args):                                                                          # module level: it is pickled for the worker pool
    return {k: _SC for k in list(X.KAPPAS) + [0.0]}, 2.5, X.MAX_ENTRIES, dict(bic_index=7, n_blocks=9)


def test_part_d_stores_totals_seconds_and_the_bic_bookkeeping(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'HPS', (5.0,))
    monkeypatch.setattr(X, 'REXS', (1.2,))
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': None)
    monkeypatch.setattr(X, '_dev_job', _fake_dev_job)
    out_path = str(tmp_path / 'raw.npz')
    X.part_D(1, out_path)
    z = X.load_raw(out_path)
    pre = 'D|x110_y630|hp5|rex1.2|'
    assert z[pre + 'secs'].tolist() == [2.5, X.MAX_ENTRIES] and z[pre + 'bic'].tolist() == [7, 9, 1]       # the last entry flags a path stopped by the cap
    assert pre + 'k2|tot' in z and pre + 'k0|tot' in z and pre + 'k8|tot' in z and 'D|x210_y630|hp5|rex1.2|bic' in z


def test_protocol_constants_are_pinned():
    """The grids and patches of the protocol (changing one is a change of the protocol: the development part has to be rerun and the freeze rewritten)."""
    assert X.KAPPAS == (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0) and X.HPS == (3.0, 5.0, 8.0) and X.REXS == (0.8, 1.2)
    assert X.DEFAULT == (5.0, 1.2, 5.0) and X.MAX_ENTRIES == 2500 and X.LAT_H == 0.5 and X.MARGIN == 2.0 and X.FIT_RADIUS == 4.0
    assert X.DEV == [(110, 630), (210, 630)] and X.TESTS == [(310, 630), (310, 730), (360, 1040), (100, 1100)]
    assert len(X.PHOT_TAGS) == 15 and len(X.CS_FILTERS) == 48 and X.BOOT == 2000


def test_baseline_scores_reads_the_catalogue_file_and_scores_the_families_and_both_daophot_rows(tmp_path, monkeypatch):
    patch = (310, 630)
    pre = 'x310_y630|'
    cat = {}
    xs = np.array([30.0, 60.0, 80.0])
    fl = A.mag_to_counts(np.array([20.0, 20.5, 21.0]))
    for tag in ('photutils', X.PHOT_TAGS[3]):
        cat[pre + tag + '|x'], cat[pre + tag + '|y'], cat[pre + tag + '|flux'] = xs, xs, fl
    cs = pre + 'crowdsource|'
    cat[cs + 'x'], cat[cs + 'y'], cat[cs + 'flux'] = xs, xs, fl
    cat[cs + 'dflux'], cat[cs + 'qf'], cat[cs + 'fracflux'], cat[cs + 'spread'] = np.ones(3), np.ones(3), np.ones(3), np.zeros(3)
    cat[pre + 'daophot|x'], cat[pre + 'daophot|y'], cat[pre + 'daophot|flux'] = xs, xs, fl                      # the shifted catalogue matches the truth
    cat[pre + 'daophot_asis|x'], cat[pre + 'daophot_asis|y'], cat[pre + 'daophot_asis|flux'] = xs, xs + 1.0, fl   # as shipped: one pixel off in y, nothing matches
    path = str(tmp_path / 'cat.npz')
    np.savez(path, **cat)
    monkeypatch.setattr(X, 'BASE', path)
    t = dict(x=xs.copy(), y=xs.copy(), mag=np.array([20.0, 20.5, 21.0]), inner=np.ones(3, bool))
    out = X.baseline_scores(patch, t)
    assert set(out) == {'photutils', X.PHOT_TAGS[3], 'daophot', 'daophot_asis'} | {X.cs_name(f) for f in X.CS_FILTERS}
    assert out['daophot']['matched22'] == 3 and out['daophot_asis']['matched22'] == 0
    assert out['photutils']['matched22'] == 3 and out[X.cs_name((0.0, 0.0, 0.0, np.inf))]['matched22'] == 3
    assert X.baseline_scores((110, 630), t) == {}                                                                   # a patch without catalogues: nothing
    monkeypatch.setattr(X, 'BASE', str(tmp_path / 'absent.npz'))
    assert X.baseline_scores(patch, t) == {}


def test_part_t2_runs_a_second_cluster_with_the_same_frozen_configuration_and_its_own_catalogues(tmp_path, monkeypatch):
    best, raw = _fake_env(tmp_path, monkeypatch)
    monkeypatch.setattr(X, 'SECOND', {'m53': dict(title='M53', root=str(tmp_path), prefix='p', truth='t.npz', patches=[(5, 6)])})
    out_path = str(tmp_path / 'raw.npz')
    np.savez_compressed(out_path, **raw)
    truth = dict(px=np.array([30.0 + 5]), py=np.array([30.0 + 6]), m606=np.array([20.0]))
    seen = {}
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': (seen.setdefault('truth', ds), truth)[1])
    monkeypatch.setattr(X, 'load_patch', lambda x0, y0, ds='m2': (seen.setdefault('patch', ds), None)[1])
    monkeypatch.setattr(X, 'run_path', lambda patch, hp, rex, **kw: ({k: (np.array([30.0]), np.array([30.0]), A.mag_to_counts(np.array([20.0]))) for k in list(X.KAPPAS) + [0.0]}, 1.0, 1,
                                                                    dict(bic_index=0, n_blocks=2)))
    sc = A.score(np.array([[30.0, 30.0]]), A.mag_to_counts(np.array([20.0])), A.truth_for_patch(truth, 5, 6))
    monkeypatch.setattr(X, 'baseline_scores', lambda p, t, base=None: (seen.setdefault('base', base), {'photutils': sc})[1])
    with pytest.raises(RuntimeError, match='no freeze record'):
        X.part_T2('m53', out_path)
    X.write_freeze(raw)
    X.part_T2('m53', out_path)
    z = X.load_raw(out_path)
    assert 'T2|m53|x5_y6|cfg' in z and z['T2|m53|x5_y6|cfg'].tolist() == [best[0], best[1], best[2]] and 'T2|m53|x5_y6|ours|k%g|tot' % best[2] in z
    assert 'T2|m53|x5_y6|photutils (default)|tot' in z and seen['truth'] == 'm53' and seen['patch'] == 'm53' and seen['base'] == X.BASE2        # its own truth, patches and catalogues
    assert not any(k.startswith('T|') for k in z)
    seen.clear()
    X.part_T2('m53', out_path)                                                                       # once
    assert 'patch' not in seen
