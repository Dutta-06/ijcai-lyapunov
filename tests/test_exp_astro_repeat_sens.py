"""The post hoc sensitivity analyses of the repeat-frame test (experiments/exp_astro_repeat_sens.py): our scorer against `astro_common`, the shift of the truth, the depth cut, the magnitude offset, the block
sizes, the bootstrap, the star-offset fits, the builder of the raw and the report on a planted raw, and the byte identity of `--from-raw`. Tests that need the stored results and the data skip without them."""
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_common as A                    # noqa: E402
import astro_repeat as R                    # noqa: E402
import exp_astro_m2 as X                    # noqa: E402
import exp_astro_repeat as E                # noqa: E402
import exp_astro_repeat_sens as S           # noqa: E402
from test_astro_report import _synthetic_catalogs   # noqa: E402

ICUT = S.ICUT


def _truth_patch(x, y, mag, inner=None):
    return dict(x=np.asarray(x, float), y=np.asarray(y, float), mag=np.asarray(mag, float), inner=np.ones(len(x), bool) if inner is None else np.asarray(inner, bool))


def _det(x, y, mag):
    return np.asarray(x, float), np.asarray(y, float), A.mag_to_counts(np.asarray(mag, float))


def _pairs(r, c):
    return sorted(zip(np.asarray(r).tolist(), np.asarray(c).tolist()))


# ---- the declared design is pinned ------------------------------------------------------------------------------------------------------

def test_the_declared_constants_are_pinned():
    assert S.ALPHAS == (0.0, 0.25, 0.5, 0.75, 1.0) and S.CUTS == (20.0, 21.0, 22.0) and S.ICUT == 2 and S.CELLS == (10, 20, 25, 50)
    assert S.TOLS == ((0.5, 0.5), (0.35, 0.5), (0.5, 0.3), (0.75, 0.75)) and S.BIG == 1.0e4
    assert S.BOOT == 2000 == X.BOOT and S.BIG_B == 40000 and S.SEEDS == tuple(range(1, 201))
    assert (S.MAG_LO, S.MAG_HI, S.MAG_MIN_PAIRS) == (17.0, 20.0, 5)
    assert S.SELECTIONS == {'A': dict(mlo=15.5, mhi=19.5, margin=15.0, riso=8.0, dmag_iso=2.0), 'B': dict(mlo=15.5, mhi=19.5, margin=8.0, riso=5.0, dmag_iso=2.5)} and S.OWN_SEL == 'B'
    assert S.STAR == dict(half=5, max_offset=1.0, min_stars=3)
    assert S.ISS_KEYS == ('iso|ours|k2', 'adq|ours|k2', 'iso|default|k5', 'adq|default|k5') and S.TUNED_BASE == ('crowdsource (tuned)', 'photutils (tuned)')
    assert S.BASE_EXTRA == ('scl|ours|k2', 'scl|default|k5', 'crowdsource (default)', 'photutils (default)') and len(S.STORED_ROWS) == 13
    assert sorted(S.SCENARIOS) == sorted(['base', 'tol35', 'mag30', 'tol75', 'scan25', 'scan50', 'scan75', 'scan100', 'own25', 'own50', 'own75', 'own100', 'ownp100', 'mo', 'mo+scan100', 'mo+own100'])
    sp = S.SCENARIOS
    assert (sp['tol35']['dpos'], sp['tol35']['dmag']) == (0.35, 0.5) and (sp['mag30']['dpos'], sp['mag30']['dmag']) == (0.5, 0.3) and (sp['tol75']['dpos'], sp['tol75']['dmag']) == (0.75, 0.75)
    assert sp['scan50']['alpha'] == 0.5 and sp['scan50']['src'] == 'scan' and sp['own25']['src'] == 'own' and sp['ownp100']['src'] == 'ownp' and sp['ownp100']['alpha'] == 1.0
    assert sp['mo']['moff'] and sp['mo+scan100']['moff'] and sp['mo+scan100']['src'] == 'scan' and sp['mo+own100']['src'] == 'own' and not sp['scan100']['moff'] and not sp['base']['moff']
    assert S.scn_name('scan', 0.0) == 'base' and S.scn_name('own', 0.0) == 'base' and S.scn_name('scan', 0.75) == 'scan75' and S.scn_name('own', 1.0) == 'own100'
    assert S.PATCHES == E.PATCHES and S.FRAMES == ['m2r1', 'm2r2', 'm2r3']


# ---- the matcher ------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('seed', [0, 1, 2])
@pytest.mark.parametrize('tol', S.TOLS)
def test_the_matcher_finds_the_matching_of_astro_common_for_every_tolerance(seed, tol):
    rng = np.random.default_rng(seed)
    n = 400
    tr = rng.uniform(0, 40, (n, 2))                                     # a crowded field: many detections have several admissible truth sources
    tm = rng.uniform(16, 23, n)
    keep = rng.uniform(size=n) > 0.2
    det = tr[keep] + rng.normal(0, 0.25, (keep.sum(), 2))
    dm = tm[keep] + rng.normal(0, 0.2, keep.sum())
    det = np.r_[det, rng.uniform(0, 40, (60, 2))]
    dm = np.r_[dm, rng.uniform(17, 23, 60)]
    r1, c1 = S.match_pairs(det, dm, tr, tm, *tol)
    r2, c2 = A.match(det, dm, tr, tm, dpos=tol[0], dmag=tol[1])
    assert len(r1) > 100 and _pairs(r1, c1) == _pairs(r2, c2)
    assert len(set(r1.tolist())) == len(r1) and len(set(c1.tolist())) == len(c1)           # one to one
    assert np.array_equal(r1, np.sort(r1))                                                    # sorted by detection


def test_the_matcher_takes_the_maximum_number_of_pairs_and_then_the_least_total_distance():
    tr_xy = np.array([[10.2, 10.0], [10.5, 10.0]])
    tm = np.array([20.0, 20.1])
    # a greedy nearest-first choice (detection 0 to truth 0) leaves detection 1 without a partner; the maximum matching has two pairs
    r, c = S.match_pairs(np.array([[10.0, 10.0], [10.45, 10.0]]), np.array([20.0, 20.1]), tr_xy, tm)
    assert _pairs(r, c) == [(0, 0), (1, 1)]
    r, c = S.match_pairs(np.array([[10.0, 10.0], [10.7, 10.0]]), np.array([20.0, 20.1]), tr_xy, tm)
    assert _pairs(r, c) == [(0, 0), (1, 1)]
    # two matchings of the maximum size {d0-t0, d1-t1} and {d0-t1, d1-t0}: the one of the smaller total distance (0.1 + 0.1 against 0.4 + 0.4)
    r, c = S.match_pairs(np.array([[10.1, 10.0], [10.6, 10.0]]), np.array([20.0, 20.1]), tr_xy, tm)
    assert _pairs(r, c) == [(0, 0), (1, 1)]
    r, c = S.match_pairs(np.array([[10.6, 10.0], [10.1, 10.0]]), np.array([20.0, 20.1]), tr_xy, tm)
    assert _pairs(r, c) == [(0, 1), (1, 0)]
    # cardinality dominates the distance: one pair at 0.0 px against two pairs at 0.7 px each (the weight BIG - d must be large enough)
    r, c = S.match_pairs(np.array([[10.0, 10.0], [9.3, 10.0]]), np.array([20.0, 20.0]), np.array([[10.0, 10.0], [10.7, 10.0]]), np.array([20.0, 20.0]), 0.75, 0.75)
    assert _pairs(r, c) == [(0, 1), (1, 0)]
    # a component with more detections than the truth sources they can match: d1 and d2 compete for t0 and the free truth sources t1 and t2 are admissible to d0 only; the third pair of an assignment
    # of the whole block would be inadmissible and is not returned
    tr3 = np.array([[10.0, 10.0], [10.8, 10.0], [10.4, 10.45]])
    de3 = np.array([[10.4, 10.0], [9.8, 10.0], [10.0, 9.7]])
    r, c = S.match_pairs(de3, np.full(3, 20.0), tr3, np.full(3, 20.0))
    assert _pairs(r, c) == [(0, 1), (1, 0)] and _pairs(r, c) == _pairs(*A.match(de3, np.full(3, 20.0), tr3, np.full(3, 20.0)))
    # two detections of one star: one false positive, the nearer one is the match
    r, c = S.match_pairs(np.array([[10.2, 10.1], [10.2, 10.0]]), np.array([20.0, 20.0]), tr_xy[:1], tm[:1])
    assert _pairs(r, c) == [(1, 0)]


def test_the_matcher_tolerances_are_inclusive_and_both_bind():
    tr, tm = np.array([[10.0, 10.0]]), np.array([20.0])
    assert len(S.match_pairs(np.array([[10.5, 10.0]]), np.array([20.5]), tr, tm)[0]) == 1                 # exactly at both tolerances
    assert len(S.match_pairs(np.array([[10.5, 10.0]]), np.array([20.5]), tr, tm, 0.49, 0.5)[0]) == 0
    assert len(S.match_pairs(np.array([[10.5, 10.0]]), np.array([20.5]), tr, tm, 0.5, 0.49)[0]) == 0
    assert len(S.match_pairs(np.array([[10.5, 10.0]]), np.array([20.5]), tr, tm, 0.5, 0.5)[0]) == 1
    assert len(S.match_pairs(np.array([[10.4, 10.0]]), np.array([20.4]), tr, tm, 0.35, 0.5)[0]) == 0         # 0.35 px
    assert len(S.match_pairs(np.array([[10.3, 10.0]]), np.array([20.4]), tr, tm, 0.5, 0.3)[0]) == 0          # 0.3 mag
    assert len(S.match_pairs(np.array([[10.7, 10.0]]), np.array([20.7]), tr, tm, 0.75, 0.75)[0]) == 1
    # the candidate search is wider than the tolerance by 1e-6 px; the exact test is the tolerance itself
    assert len(S.match_pairs(np.array([[10.5 + 5e-7, 10.0]]), np.array([20.0]), tr, tm)[0]) == 0
    assert len(S.match_pairs(np.array([[10.5 - 5e-7, 10.0]]), np.array([20.0]), tr, tm)[0]) == 1
    # the distance is the Euclidean one in both coordinates
    assert len(S.match_pairs(np.array([[10.4, 10.4]]), np.array([20.0]), tr, tm)[0]) == 0                   # 0.566 px
    assert len(S.match_pairs(np.array([[10.35, 10.35]]), np.array([20.0]), tr, tm)[0]) == 1                  # 0.495 px


def test_the_matcher_handles_empty_inputs_and_nan_magnitudes():
    z2, z1 = np.zeros((0, 2)), np.zeros(0)
    for a in (S.match_pairs(z2, z1, np.array([[1.0, 1.0]]), np.array([20.0])), S.match_pairs(np.array([[1.0, 1.0]]), np.array([20.0]), z2, z1), S.match_pairs(z2, z1, z2, z1)):
        assert len(a[0]) == 0 and len(a[1]) == 0
    r, c = S.match_pairs(np.array([[1.0, 1.0], [1.0, 1.0]]), np.array([np.nan, 20.0]), np.array([[1.0, 1.0]]), np.array([20.0]))
    assert _pairs(r, c) == [(1, 0)]
    assert len(S.match_pairs(np.array([[5.0, 5.0]]), np.array([20.0]), np.array([[1.0, 1.0]]), np.array([20.0]))[0]) == 0


# ---- the scorer against astro_common.score ----------------------------------------------------------------------------------------------

def _edge_catalogue():
    """A truth list and a catalogue with the difficult members: detections in the edge zone, with a non-positive flux, a duplicate and faint ones."""
    rng = np.random.default_rng(3)
    truth, cats = _synthetic_catalogs(rng)
    x, y, f = (np.asarray(v, float) for v in cats['mid'])
    x = np.r_[x, 1.9, 2.0, 97.0, 97.1, 50.0, truth['x'][0] + 0.05]
    y = np.r_[y, 50.0, 50.0, 50.0, 50.0, 50.0, truth['y'][0]]
    f = np.r_[f, A.mag_to_counts(19.0), A.mag_to_counts(19.0), A.mag_to_counts(19.0), A.mag_to_counts(19.0), -5.0, A.mag_to_counts(truth['mag'][0])]
    return truth, (x, y, f)


@pytest.mark.parametrize('which', ['good', 'mid', 'poor', 'bic', 'edge'])
def test_rescore_equals_astro_common_score_at_the_default_scoring(which):
    truth, cats = _synthetic_catalogs(np.random.default_rng(5))
    if which == 'edge':
        truth, det = _edge_catalogue()
    else:
        det = cats[which]
    ref = A.score(np.c_[det[0], det[1]], det[2], truth)
    r = S.rescore(det, truth, cuts=S.CUTS, cells=(10, 20, 25, 50))
    tot = [ref['true22'], ref['matched22'], ref['det22'], ref['det22_matched'], ref['true_all'], ref['matched_all'], ref['det_all']]
    assert r['tot'][ICUT].tolist() == tot and r['tot'].shape == (3, 7)
    assert np.array_equal(r['blk'][20][ICUT], ref['blk22']) and np.array_equal(r['off'], ref['offset'], equal_nan=True)
    assert r['moff'] == 0.0 and ref['true_all'] > 40 and ref['matched_all'] > 20


def test_rescore_edge_rule_and_counts_of_the_declared_vector():
    truth = _truth_patch([50.0, 1.0, 97.0, 98.0], [50.0, 50.0, 50.0, 50.0], [20.0, 20.0, 21.0, 21.0], inner=[True, False, True, False])
    det = _det([50.0, 1.9, 97.0, 97.2, 50.0, 2.0], [50.0, 50.0, 50.0, 50.0, 60.0, 80.0], [20.0, 19.0, 21.0, 21.0, 22.5, 25.0])
    r = S.rescore(det, truth, cuts=(22.0,))
    # scored truth: 2 (the edge-zone ones are not scored); scored detections: x = 1.9 and 97.2 are outside [2, 97], the third scored one (22.5 mag) is fainter than the cut
    assert r['tot'][0].tolist() == [2, 2, 2, 2, 2, 2, 4] and r['off'].shape == (2,) and np.allclose(r['off'], [0.0, 0.0])           # the detection at x = 2.0 exactly is scored (it is 25 mag)
    assert np.isnan(S.rescore(_det([5.0], [5.0], [30.0]), truth, cuts=(22.0,))['off']).all()                      # no bright match: NaN offsets


def test_rescore_counts_at_the_depth_cuts_and_the_matching_does_not_depend_on_them():
    # truth magnitudes 19.9, 20.5, 21.5, 22.5; detections 20.3 (matches 19.9), 20.6 (matches 20.5), 21.2 (matches 21.5), 22.8 (matches 22.5)
    truth = _truth_patch([10.0, 30.0, 50.0, 70.0], [10.0, 30.0, 50.0, 70.0], [19.9, 20.5, 21.5, 22.5])
    det = _det([10.1, 30.1, 50.1, 70.1], [10.0, 30.0, 50.0, 70.0], [20.3, 20.6, 21.2, 22.8])
    r = S.rescore(det, truth, cuts=(20.0, 21.0, 22.0))
    assert r['tot'][:, :4].tolist() == [[1, 1, 0, 0], [2, 2, 2, 2], [3, 3, 3, 3]]                      # [true_c, matched_c (by truth mag), det_c, detm_c (by detection mag)] at 20, 21, 22
    assert (r['tot'][:, 4:] == [4, 4, 4]).all()                                                          # without a cut: 4 truth, 4 matched, 4 detections at every cut
    # the cut is exclusive: truth sources of exactly 20, 21 and 22 mag are not below the cuts 20, 21 and 22
    tex = _truth_patch([10.0, 30.0, 50.0], [10.0, 30.0, 50.0], [20.0, 21.0, 22.0])
    rex = S.rescore(_det([10.1, 30.1, 50.1], [10.0, 30.0, 50.0], [20.0, 21.0, 22.0]), tex, cuts=(20.0, 21.0, 22.0))
    assert rex['tot'][:, 0].tolist() == [0, 1, 2] and rex['tot'][:, 1].tolist() == [0, 1, 2] and rex['tot'][:, 4].tolist() == [3, 3, 3]
    one = S.rescore(det, truth, cuts=(21.0,))
    assert one['tot'][0].tolist() == r['tot'][1].tolist()
    assert np.array_equal(r['blk'][20][1], one['blk'][20][0])


def test_rescore_blocks_count_the_truth_and_its_match_in_the_truth_block_and_detections_in_theirs():
    # truth at (19.7, 5.0) is in block 0 (cells of 20 px start at -0.5: [-0.5, 19.5) is block 0), its match at (19.9, 5.0) is in block 1
    truth = _truth_patch([19.0, 19.7], [5.0, 60.0], [20.0, 20.0])
    det = _det([19.0, 20.2], [5.0, 60.0], [20.0, 20.0])
    r = S.rescore(det, truth, cuts=(22.0,), cells=(20,))
    blk = r['blk'][20][0]
    assert blk.shape == (25, 4) and blk.sum() == 8
    # truth 0 in block 0; truth 1 (x 19.7 -> block column 1, row 3) is block 3*5 + 1 = 16; the second detection at x = 20.2 is in column 1 as well
    assert blk[0, 0] == 1 and blk[16, 0] == 1 and blk[0, 1] == 1 and blk[16, 1] == 1 and blk[0, 2] == 1 and blk[16, 2] == 1 and blk[0, 3] == 1 and blk[16, 3] == 1
    # the truth at x = 19.4 is in column 0, at 19.6 in column 1; a match 0.3 px away keeps the truth's block
    tr2 = _truth_patch([19.4], [5.0], [20.0])
    r2 = S.rescore(_det([19.7], [5.0], [20.0]), tr2, cuts=(22.0,), cells=(20,))
    assert r2['blk'][20][0][0, 0] == 1 and r2['blk'][20][0][0, 1] == 1 and r2['blk'][20][0][1, 2] == 1 and r2['blk'][20][0][1, 3] == 1


def test_block_index_follows_astro_common_for_20_px_and_the_other_sizes_partition_the_patch():
    rng = np.random.default_rng(1)
    pts = np.r_[rng.uniform(-0.5, 99.4, (500, 2)), [[-0.5, -0.5], [99.4, 99.4], [19.49, 19.49], [19.51, 19.51], [0.0, 99.0]]]
    assert np.array_equal(S.block_index(pts, 20), A.block_of(pts, 100))
    for cell in S.CELLS:
        idx = S.block_index(pts, cell)
        assert idx.min() >= 0 and idx.max() == (100 // cell) ** 2 - 1
    assert S.block_index([[-0.5, -0.5]], 10)[0] == 0 and S.block_index([[99.4, 99.4]], 50)[0] == 3 and S.block_index([[24.4, 24.6]], 25)[0] == 4
    assert S.block_index([[99.7, 99.7]], 20)[0] == 24 and S.block_index([[-3.0, -3.0]], 20)[0] == 0 and S.block_index([[-3.0, 99.7]], 20)[0] == 20 and S.block_index([[99.7, -3.0]], 20)[0] == 4
    assert S.block_index([[9.4, 9.6]], 10)[0] == 10 and S.block_index([[9.6, 9.4]], 10)[0] == 1


def test_rescore_blocks_of_every_size_sum_to_the_same_counts():
    truth, cats = _synthetic_catalogs(np.random.default_rng(8))
    r = S.rescore(cats['good'], truth, cuts=S.CUTS, cells=S.CELLS)
    for cell in S.CELLS:
        assert r['blk'][cell].shape == (3, (100 // cell) ** 2, 4)
        for k in range(3):
            assert r['blk'][cell][k].sum(0).tolist() == r['tot'][k][:4].tolist()


def test_the_offset_is_the_median_of_the_matches_brighter_than_21_mag():
    truth = _truth_patch([10.0, 30.0, 50.0, 70.0, 90.0], np.full(5, 40.0), [19.0, 20.0, 20.9, 21.5, 21.8])
    det = _det([10.1, 30.2, 50.3, 70.4, 90.45], np.full(5, 40.0), [19.0, 20.0, 20.9, 21.5, 21.8])
    r = S.rescore(det, truth, cuts=(22.0,))
    assert np.allclose(r['off'], [0.2, 0.0]) and r['tot'][0][5] == 5                                    # the median of 0.1, 0.2, 0.3 (the pairs of 21.5 and 21.8 mag are not bright)
    r = S.rescore(_det([10.1], [40.0], [21.0]), _truth_patch([10.0], [40.0], [21.0]), cuts=(22.0,))
    assert np.isnan(r['off']).all() and r['tot'][0][5] == 1                                               # 21.0 is not brighter than 21


def test_the_shift_moves_the_truth_for_the_matching_only_and_the_scored_set_and_blocks_stay():
    truth = _truth_patch([19.4, 50.0], [50.0, 50.0], [20.0, 20.0])
    det = _det([19.9, 50.7], [50.0, 50.0], [20.0, 20.0])                    # the second detection is 0.7 px from its truth source
    r0 = S.rescore(det, truth, cuts=(22.0,))
    assert r0['tot'][0].tolist() == [2, 1, 2, 1, 2, 1, 2]
    r1 = S.rescore(det, truth, shift=(0.5, 0.0), cuts=(22.0,))                # +0.5 px in x: the truth moves toward the detections
    assert r1['tot'][0].tolist() == [2, 2, 2, 2, 2, 2, 2]
    assert np.allclose(r1['off'], [(19.9 - 19.9 + 50.7 - 50.5) / 2.0, 0.0])                                  # (detection - moved truth) of the two matches: 0.0 and 0.2
    r2 = S.rescore(det, truth, shift=(-0.5, 0.0), cuts=(22.0,))
    assert r2['tot'][0].tolist() == [2, 0, 2, 0, 2, 0, 2]                     # the other way: neither is within 0.5 px
    # the truth source at 19.4 moved to 19.9 stays in block column 0 (its unshifted block); the detection at 19.9 is in column 1
    b = r1['blk'][20][0]
    assert b[10, 0] == 1 and b[11, 0] == 0 and b[10, 1] == 1 and b[11, 1] == 0 and b[11, 2] == 1 and b[11, 3] == 1
    # the offset of the matches is (detection - moved truth)
    r3 = S.rescore(_det([19.9], [50.2], [20.0]), _truth_patch([19.4], [50.0], [20.0]), shift=(0.3, 0.1), cuts=(22.0,))
    assert np.allclose(r3['off'], [19.9 - 19.7, 50.2 - 50.1])
    # the denominators do not depend on the shift
    assert r0['tot'][0][0] == r1['tot'][0][0] == r2['tot'][0][0] and r0['tot'][0][4] == r2['tot'][0][4]


def test_the_magnitude_offset_is_estimated_from_17_to_20_mag_matches_and_subtracted_from_the_detections():
    n = 12
    x = np.linspace(8, 90, n)
    tm = np.linspace(17.2, 19.8, n)
    truth = _truth_patch(x, np.full(n, 30.0), tm)
    # detections 0.45 mag too faint, scattered by -0.1 / +0.1: the pairs with +0.1 (0.55 mag) are outside the 0.5 mag tolerance until the offset is removed
    scatter = np.where(np.arange(n) % 2 == 0, -0.1, 0.1)
    det = _det(x + 0.05, np.full(n, 30.0), tm + 0.45 + scatter)
    plain = S.rescore(det, truth, cuts=(22.0,))
    assert plain['tot'][0][1] == n // 2 and plain['moff'] == 0.0
    auto = S.rescore(det, truth, moff='auto', cuts=(22.0,))
    assert abs(auto['moff'] - 0.35) < 1e-9 and auto['tot'][0][1] == n and auto['tot'][0][3] == n               # the first pass saw the six pairs at 0.35 mag
    # the offset is the MEDIAN of the pass-1 pairs: Delta m of 0.30, 0.32, 0.34, 0.36, 0.38 and 0.49 (median 0.35, mean 0.365)
    skew = np.array([0.30, 0.32, 0.34, 0.36, 0.38, 0.49])
    sk = S.rescore(_det(x[:6] + 0.05, np.full(6, 30.0), tm[:6] + skew), _truth_patch(x[:6], np.full(6, 30.0), tm[:6]), moff='auto', cuts=(22.0,))
    assert abs(sk['moff'] - 0.35) < 1e-9
    given = S.rescore(det, truth, moff=0.35, cuts=(22.0,))
    assert given['tot'].tolist() == auto['tot'].tolist() and given['moff'] == 0.35
    # the first pass has the tolerances of the scenario: an offset beyond the tolerance finds no pair and no offset
    far = S.rescore(_det(x + 0.05, np.full(n, 30.0), tm + 0.45), truth, dmag=0.3, moff='auto', cuts=(22.0,))
    assert far['moff'] == 0.0 and far['tot'][0][1] == 0
    exact = S.rescore(_det(x + 0.05, np.full(n, 30.0), tm + 0.45), truth, moff='auto', cuts=(22.0,))
    assert abs(exact['moff'] - 0.45) < 1e-9
    # the corrected magnitudes enter the depth cut: detections of 22.15 and 22.35 mag (true 21.7, 21.9) are fainter than 22 uncorrected and brighter once 0.45 is removed
    truth2 = _truth_patch(np.r_[x, 20.0, 40.0], np.r_[np.full(n, 30.0), 70.0, 70.0], np.r_[tm, 21.7, 21.9])
    det2 = _det(np.r_[x + 0.05, 20.0, 40.0], np.r_[np.full(n, 30.0), 70.0, 70.0], np.r_[tm + 0.45, 22.15, 22.35])
    p2, q2 = S.rescore(det2, truth2, cuts=(22.0,)), S.rescore(det2, truth2, moff='auto', cuts=(22.0,))
    assert p2['tot'][0][2] == n and q2['tot'][0][2] == n + 2 and q2['tot'][0][3] == n + 2
    # a negative offset (detections too bright) is added to the magnitudes
    det3 = _det(x + 0.05, np.full(n, 30.0), tm - 0.4)
    assert abs(S.rescore(det3, truth, moff='auto', cuts=(22.0,))['moff'] + 0.4) < 1e-9
    # fewer than 5 matched pairs between 17 and 20 mag: no offset
    few = S.rescore(_det(x[:4] + 0.05, np.full(4, 30.0), tm[:4] + 0.45), truth, moff='auto', cuts=(22.0,))
    assert few['moff'] == 0.0
    five = S.rescore(_det(x[:5] + 0.05, np.full(5, 30.0), tm[:5] + 0.45), truth, moff='auto', cuts=(22.0,))
    assert abs(five['moff'] - 0.45) < 1e-9
    # only the pairs between 17 and 20 mag count: matches at 20.2 and 16.5 mag are ignored
    tr3 = _truth_patch(np.r_[x[:5], 70.0, 80.0], np.full(7, 30.0), np.r_[tm[:5], 20.2, 16.5])
    d3 = _det(np.r_[x[:5] + 0.05, 70.0, 80.0], np.full(7, 30.0), np.r_[tm[:5] + 0.1, 20.6, 16.9])
    assert abs(S.rescore(d3, tr3, moff='auto', cuts=(22.0,))['moff'] - 0.1) < 1e-9
    # the boundaries 17 and 20 are exclusive: with the boundary pairs the sample would have five pairs, without them three
    for mags in ([17.0, 17.0, 18.0, 18.5, 19.0], [20.0, 20.0, 18.0, 18.5, 19.0], [17.0, 20.0, 18.0, 18.5, 19.0]):
        tr4 = _truth_patch(x[:5], np.full(5, 30.0), mags)
        d4 = _det(x[:5] + 0.05, np.full(5, 30.0), np.array(mags) + 0.3)
        assert S.rescore(d4, tr4, moff='auto', cuts=(22.0,))['moff'] == 0.0, mags
        assert S.rescore(d4, tr4, moff=0.1, cuts=(22.0,))['moff'] == 0.1
    ok4 = _truth_patch(x[:5], np.full(5, 30.0), [17.01, 19.99, 18.0, 18.5, 19.0])
    assert abs(S.rescore(_det(x[:5] + 0.05, np.full(5, 30.0), np.array([17.01, 19.99, 18.0, 18.5, 19.0]) + 0.3), ok4, moff='auto', cuts=(22.0,))['moff'] - 0.3) < 1e-9


def test_the_magnitude_offset_pass_uses_the_shifted_truth():
    n = 8
    x = np.linspace(10, 80, n)
    tm = np.linspace(17.5, 19.5, n)
    truth = _truth_patch(x, np.full(n, 30.0), tm)
    det = _det(x, np.full(n, 30.6), tm + 0.2)                                   # 0.6 px off in y: not matched unless the truth is moved
    assert S.rescore(det, truth, moff='auto', cuts=(22.0,))['moff'] == 0.0
    a = S.rescore(det, truth, shift=(0.0, 0.5), moff='auto', cuts=(22.0,))
    assert abs(a['moff'] - 0.2) < 1e-9 and a['tot'][0][1] == n


# ---- the bootstrap -----------------------------------------------------------------------------------------------------------------------

def _blocks(seed, nb=50):
    rng = np.random.default_rng(seed)
    t = rng.integers(8, 30, nb)
    m = (t * rng.uniform(0.3, 0.7, nb)).astype(int)
    d = rng.integers(8, 30, nb)
    dm = np.minimum((d * rng.uniform(0.3, 0.7, nb)).astype(int), d)
    return np.c_[t, m, d, dm].astype(float)


def test_the_resampled_sums_equal_the_gathered_sums_exactly():
    blk = _blocks(1)
    idx = X.boot_idx_strat([25, 25], 300, 4)
    assert np.array_equal(S.resample_sums(blk, idx), blk[idx].sum(1))
    idx1 = X.boot_idx(50, 77, 1)
    assert np.array_equal(S.resample_sums(blk, idx1), blk[idx1].sum(1))


@pytest.mark.parametrize('nb_patches', [1, 2, 6])
def test_boot_dF_is_the_repository_bootstrap_at_2000_resamples_seed_0(nb_patches):
    a, b = _blocks(2, 25 * nb_patches), _blocks(3, 25 * nb_patches)
    sizes = [25] * nb_patches
    d, ci = S.boot_dF(a, b, sizes)
    ref = X.boot_diff(a, b, X.boot_idx_strat(sizes))[2]
    assert d == ref[0] and np.array_equal(ci, ref[1])


def test_boot_stat_chunks_do_not_change_the_result_and_seed_and_B_do():
    a, b = _blocks(4), _blocks(5)
    r1 = S.boot_stat(S._dF, [a, b], [25, 25], 1000, 0, chunk=1000)
    r2 = S.boot_stat(S._dF, [a, b], [25, 25], 1000, 0, chunk=77)
    assert r1[0] == r2[0] and np.array_equal(r1[1], r2[1])
    r3 = S.boot_dF(a, b, [25, 25], 1000, 1)
    assert not np.array_equal(r3[1], r1[1]) and r3[0] == r1[0]
    r4 = S.boot_dF(a, b, [25, 25], 3000, 0)
    assert not np.array_equal(r4[1], S.boot_dF(a, b, [25, 25], 1000, 0)[1])
    ref = X.boot_diff(a, b, X.boot_idx_strat([25, 25], 3000, 0))[2]
    assert np.array_equal(r4[1], ref[1])


def test_the_seed_lowers_are_the_lower_ends_of_the_intervals_of_each_seed():
    a, b = _blocks(6), _blocks(7)
    lows = S.seed_lowers(a, b, [25, 25], seeds=(1, 2, 3), B=500)
    assert len(lows) == 3 and len(set(lows.tolist())) == 3
    for s, v in zip((1, 2, 3), lows):
        assert v == X.boot_diff(a, b, X.boot_idx_strat([25, 25], 500, s))[2][1][0]


def test_boot_ratio_is_the_paired_ratio_of_the_F_values():
    a, b = _blocks(8), _blocks(9)
    r, ci = S.boot_ratio(a, b, [25, 25])
    assert np.isclose(r, X.cpf(a.sum(0))[2] / X.cpf(b.sum(0))[2], rtol=0, atol=1e-15) and ci[0] < r < ci[1]
    idx = X.boot_idx_strat([25, 25])
    ref = np.percentile(X.cpf(a[idx].sum(1))[2] / X.cpf(b[idx].sum(1))[2], [2.5, 97.5])
    assert np.array_equal(ci, ref)
    # identical rows: dF = 0 with a degenerate interval and the ratio 1
    d0, ci0 = S.boot_dF(a, a, [25, 25])
    assert d0 == 0.0 and np.array_equal(ci0, [0.0, 0.0])
    r1, ci1 = S.boot_ratio(a, a, [25, 25])
    assert r1 == 1.0 and np.array_equal(ci1, [1.0, 1.0])


# ---- the star offsets --------------------------------------------------------------------------------------------------------------------

def _gauss_stamp(sx=1.3, sy=1.3, n=51):
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    st = np.exp(-0.5 * (((xx - c) / sx) ** 2 + ((yy - c) / sy) ** 2))
    return st / st.sum()


def _star_image(stamp, positions, amps, shape=(100, 100), sky=150.0):
    """The image of the stamp's stars at the (x, y) positions with the amplitudes (sums), by the same cubic-spline rendering the fit uses."""
    from scipy.ndimage import map_coordinates, spline_filter
    coef = spline_filter(stamp, order=3)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
    img = np.full(shape, sky)
    c = (stamp.shape[0] - 1) / 2.0
    for (cx, cy), a in zip(positions, amps):
        img += a * map_coordinates(coef, [yy - cy + c, xx - cx + c], order=3, mode='constant', cval=0.0, prefilter=False)
    return img


def test_fit_star_recovers_a_planted_position_and_flux_and_reports_convergence():
    st = _gauss_stamp()
    from scipy.ndimage import spline_filter
    coef = spline_filter(st, order=3)
    img = _star_image(st, [(40.3, 55.8)], [8000.0])
    cx, cy, a, b, ok = S.fit_star(img, coef, 4.6, 1.0, 40.0, 56.0)
    assert ok and abs(cx - 40.3) < 1e-3 and abs(cy - 55.8) < 1e-3 and abs(a - 8000.0) < 1.0 and abs(b - 150.0) < 0.1
    # image rows are y and columns x: a star displaced in x only is recovered in x only
    cx, cy, *_ = S.fit_star(_star_image(st, [(42.0, 56.0)], [8000.0]), coef, 4.6, 1.0, 40.0, 56.0)
    assert abs(cx - 42.0) < 1e-3 and abs(cy - 56.0) < 1e-3
    # noise: Poisson plus dark, still within 0.05 px for a star of 8000 DN
    rng = np.random.default_rng(0)
    noisy = img + rng.normal(0, np.sqrt(np.maximum(img, 1) / 4.6 + 1.0))
    cx, cy, *_ = S.fit_star(noisy, coef, 4.6, 1.0, 40.0, 56.0)
    assert abs(cx - 40.3) < 0.05 and abs(cy - 55.8) < 0.05


def test_the_weights_of_the_star_fit():
    w = S.pixel_weights(np.array([0.0, 4.6, 46.0]), 4.6, 1.0)
    assert np.allclose(w, [1 / np.sqrt(1 / 4.6 + 1.0), 1 / np.sqrt(2.0), 1 / np.sqrt(11.0)], rtol=0, atol=1e-15)


def test_a_fit_that_did_not_converge_is_dropped(monkeypatch):
    st = _gauss_stamp()
    img = _star_image(st, [(30.0, 30.0)], [9000.0])
    rows = dict(x=np.array([30.0]), y=np.array([30.0]), mag=np.array([17.0]), nf=np.ones(1))
    assert len(S.measure_patch(img, st, 4.6, 1.0, rows, 'B')[0]) == 1
    real = S.fit_star
    monkeypatch.setattr(S, 'fit_star', lambda *a, **k: real(*a, **k)[:4] + (False,))
    stars, nsel = S.measure_patch(img, st, 4.6, 1.0, rows, 'B')
    assert nsel == 1 and len(stars) == 0


def test_measure_patch_returns_fitted_minus_truth_for_the_selected_stars_only():
    st = _gauss_stamp()
    truth_pos = [(30.0, 30.0), (70.0, 40.0), (50.0, 75.0), (5.0, 50.0), (95.0, 95.0)]
    img = _star_image(st, [(x + dx, y + dy) for (x, y), (dx, dy) in zip(truth_pos, [(0.2, -0.1), (0.2, -0.1), (0.2, -0.1), (0.2, -0.1), (0.2, -0.1)])], [9000.0] * 5)
    rows = dict(x=np.array([p[0] for p in truth_pos]), y=np.array([p[1] for p in truth_pos]), mag=np.array([17.0, 17.5, 18.0, 17.0, 17.0]), nf=np.ones(5))
    stars, nsel = S.measure_patch(img, st, 4.6, 1.0, rows, 'B')
    # the stars at x = 5 and (95, 95) are within 8 px of an edge (selection B): 3 are selected and kept
    assert nsel == 3 and stars.shape == (3, 6)
    assert np.allclose(stars[:, 3], 0.2, atol=2e-3) and np.allclose(stars[:, 4], -0.1, atol=2e-3)           # fitted - truth
    assert np.allclose(stars[:, 0], [30.0, 70.0, 50.0]) and np.allclose(stars[:, 2], [17.0, 17.5, 18.0])
    stars_a, nsel_a = S.measure_patch(img, st, 4.6, 1.0, rows, 'A')                                       # selection A needs 15 px from the edge: only (30, 30), (70, 40), (50, 75)
    assert nsel_a == 3
    # an offset above the largest one kept (1 px) is dropped
    far = _star_image(st, [(30.0 + 1.2, 30.0)], [9000.0])
    s2, n2 = S.measure_patch(far, st, 4.6, 1.0, dict(x=np.array([30.0]), y=np.array([30.0]), mag=np.array([17.0]), nf=np.ones(1)), 'B')
    assert n2 == 1 and len(s2) == 0


def test_the_star_selection_rules_one_by_one():
    base = dict(x=[50.0], y=[50.0], mag=[17.0], nf=[1.0])

    def rows(**extra):
        d = {k: list(v) for k, v in base.items()}
        for k, v in extra.items():
            for kk, vv in zip(('x', 'y', 'mag', 'nf'), v):
                d[kk].append(vv)
        return {k: np.array(v, float) for k, v in d.items()}
    assert S.select_stars(rows()).tolist() == [0]
    for mg in (15.49, 19.5, 19.51, 12.0):                                                                                                                # outside [15.5, 19.5)
        assert len(S.select_stars(dict(x=np.array([50.0]), y=np.array([50.0]), mag=np.array([mg]), nf=np.ones(1)))) == 0, mg
    assert len(S.select_stars(dict(x=np.array([50.0]), y=np.array([50.0]), mag=np.array([15.5]), nf=np.ones(1)))) == 1                                  # 15.5 is in
    assert len(S.select_stars(dict(x=np.array([50.0]), y=np.array([50.0]), mag=np.array([19.49]), nf=np.ones(1)))) == 1
    assert len(S.select_stars(dict(x=np.array([50.0]), y=np.array([50.0]), mag=np.array([17.0]), nf=np.array([99.0])))) == 0                           # HST-saturated
    for x, y, expected in ((14.9, 50.0, 0), (15.0, 50.0, 1), (84.0, 50.0, 1), (84.1, 50.0, 0), (50.0, 14.9, 0), (50.0, 84.1, 0), (50.0, 15.0, 1), (50.0, 84.0, 1)):
        assert len(S.select_stars(dict(x=np.array([x]), y=np.array([y]), mag=np.array([17.0]), nf=np.ones(1)), mlo=15.5, mhi=19.5, margin=15.0)) == expected, (x, y)
    # neighbours (selection A: 8 px and 2 mag; B: 5 px and 2.5 mag)
    ra = S.select_stars(rows(n=(57.0, 50.0, 18.5, 1.0)), **S.SELECTIONS['A'])                         # 7 px away and 1.5 mag fainter: a neighbour in selection A (8 px, 2 mag)
    assert 0 not in ra.tolist()
    rb = S.select_stars(rows(n=(57.0, 50.0, 18.5, 1.0)), **S.SELECTIONS['B'])
    assert 0 in rb.tolist()                                                                           # beyond 5 px: allowed in B
    assert 0 not in S.select_stars(rows(n=(54.9, 50.0, 18.5, 1.0)), **S.SELECTIONS['B']).tolist()     # within 5 px and brighter than m + 2.5
    assert 0 in S.select_stars(rows(n=(55.0, 50.0, 18.5, 1.0)), **S.SELECTIONS['B']).tolist()          # at 5 px exactly: allowed (strictly closer than riso is a neighbour)
    assert 0 in S.select_stars(rows(n=(54.0, 50.0, 19.6, 1.0)), **S.SELECTIONS['B']).tolist()          # within 5 px but fainter than m + 2.5 (19.5): allowed
    assert 0 not in S.select_stars(rows(n=(54.0, 50.0, 19.4, 1.0)), **S.SELECTIONS['B']).tolist()
    assert 0 not in S.select_stars(rows(n=(54.0, 50.0, 13.0, 99.0)), **S.SELECTIONS['B']).tolist()     # a saturated neighbour counts at its listed magnitude
    # the neighbour may be at any magnitude and is itself a star only if it passes the rules
    both = S.select_stars(dict(x=np.array([50.0, 53.0]), y=np.array([50.0, 50.0]), mag=np.array([17.0, 17.2]), nf=np.ones(2)), **S.SELECTIONS['B'])
    assert len(both) == 0


def test_frame_offsets_median_and_the_per_patch_fallback():
    st1 = np.zeros((4, 6))
    st1[:, 3], st1[:, 4] = [0.1, 0.2, 0.3, 0.9], [0.5, 0.6, 0.7, 0.0]
    st2 = np.zeros((2, 6))
    st2[:, 3], st2[:, 4] = [0.4, 0.5], [0.1, 0.2]
    fm, per = S.frame_offsets({(1, 1): st1, (2, 2): st2})
    allx, ally = [0.1, 0.2, 0.3, 0.9, 0.4, 0.5], [0.5, 0.6, 0.7, 0.0, 0.1, 0.2]
    assert np.allclose(fm, [np.median(allx), np.median(ally)])
    assert np.allclose(per[(1, 1)], [np.median([0.1, 0.2, 0.3, 0.9]), np.median([0.5, 0.6, 0.7, 0.0])]) and np.allclose(per[(2, 2)], fm)          # 2 stars < 3: the frame's median
    assert S.frame_offsets({(1, 1): st2}) == (None, {})                                                  # fewer than 3 stars in the frame: none
    st3 = np.zeros((3, 6))
    st3[:, 3], st3[:, 4] = [0.1, 0.2, 0.6], [0.3, 0.4, 0.5]
    fm3, per3 = S.frame_offsets({(1, 1): st3})
    assert np.allclose(fm3, [0.2, 0.4]) and np.allclose(per3[(1, 1)], [0.2, 0.4])
    fm4, per4 = S.frame_offsets({(1, 1): st3, (2, 2): st2})
    assert np.allclose(per4[(1, 1)], [0.2, 0.4]) and np.allclose(per4[(2, 2)], fm4) and not np.allclose(fm4, [0.2, 0.4])          # a patch with exactly 3 stars has its own median
    assert S.frame_offsets({(1, 1): np.zeros((0, 6)), (2, 2): st2})[0] is None


def test_robust_sd_and_the_median_interval():
    v = np.array([0.0, 0.1, -0.1, 0.05, -0.05, 3.0])
    assert np.isclose(S.rob_sd(v), 1.4826 * np.median(np.abs(v - np.median(v))), rtol=0, atol=1e-15) and np.isnan(S.rob_sd([]))
    m, ci = S.median_interval(v)
    idx = np.random.default_rng(0).integers(0, 6, size=(S.BOOT, 6))
    assert m == np.median(v) and np.array_equal(ci, np.percentile(np.median(v[idx], axis=1), [2.5, 97.5]))
    m1, ci1 = S.median_interval(np.array([0.3]))
    assert m1 == 0.3 and np.array_equal(ci1, [0.3, 0.3])


def test_the_neighbour_bias_grows_with_flux_and_falls_with_distance():
    t = S.neighbour_bias(_gauss_stamp(1.27, 1.27), 4.6, 1.0, ds=(3.0, 5.0, 8.0), qs=(0.1, 0.4), angles=4)
    assert t.shape == (3, 2) and (t[:, 1] > t[:, 0]).all() and (t[0] > t[1]).all() and (t[1] > t[2]).all() and t[2, 0] < 0.005 and t[0, 1] > 0.1


def test_the_unit_pairs_of_two_catalogues():
    rng = np.random.default_rng(0)
    n = 40
    base = np.c_[rng.uniform(5, 95, (n, 2)), A.mag_to_counts(rng.uniform(17, 21, n))]
    other = base.copy()
    other[:, 0] -= 0.10
    other[:, 1] += 0.05
    other[:5, 2] *= 1.5                                                    # five pairs whose flux ratio (1 / 1.5 = 0.667) is outside [0.8, 1.25]
    n_, dx, dy = S.unit_pair_difference(base, other, 22.0)
    assert n_ == n - 5 and np.isclose(dx, 0.10) and np.isclose(dy, -0.05)
    assert S.unit_pair_difference(base, other, 17.0)[0] == 0 or S.unit_pair_difference(base, other, 17.0)[0] < n - 5
    far = other.copy()
    far[:, 0] += 0.5                                                       # farther than 0.4 px: no pair
    assert S.unit_pair_difference(base, far, 22.0)[0] == 0 and np.isnan(S.unit_pair_difference(base, far, 22.0)[1])
    assert S.unit_pair_difference(np.zeros((0, 3)), base, 22.0)[0] == 0
    # the radius is exclusive and the flux-ratio bounds are inclusive
    a3, b3 = np.array([[10.0, 10.0, 800.0], [30.0, 30.0, 1250.0]]), np.array([[10.5, 10.0, 1000.0], [30.0, 30.4, 1000.0]])
    assert S.unit_pair_difference(a3, b3, 30.0, rad=0.5)[0] == 1                                          # 0.5 px is not closer than 0.5 px; the pair at 0.4 px has ratio 1.25 (in)
    assert S.unit_pair_difference(a3, b3, 30.0, rad=0.51)[0] == 2                                         # 800 / 1000 = 0.8 (in)
    assert S.unit_pair_difference(a3, b3, 30.0, rad=0.51, lo=0.81)[0] == 1 and S.unit_pair_difference(a3, b3, 30.0, rad=0.51, hi=1.24)[0] == 1
    # not mutual: a second unit of a next to the same unit of b is not paired with it
    a2 = np.array([[10.0, 10.0, 1000.0], [10.2, 10.0, 1000.0]])
    b2 = np.array([[10.15, 10.0, 1000.0]])
    assert S.unit_pair_difference(a2, b2, 30.0)[0] == 1


# ---- a planted world: the inputs of the builder in the layout of the real ones ---------------------------------------------------------------

OFFS = {'m2r1': (0.10, 0.10), 'm2r2': (-0.05, 0.40), 'm2r3': (0.0, 0.05)}               # the planted truth-to-image offsets of the three frames (image = truth + offset)
CAPPED = {('m2r2', (133, 281)), ('m2r3', (370, 441)), ('m2r3', (110, 501))}              # the patches whose isotropic default run is planted as capped
GENERAL = {'m2r2'}                                                                        # the frame where the planted rule picks a general kernel
TUNED_CS = 'crowdsource|snr5|ff0|qf0|spinf'


def _catalogue(rng, tx, ty, tm, comp, jit, nfp, off=(0.0, 0.0), mbias=0.0):
    """A catalogue (x, y, flux) in patch pixels: a fraction comp of the truth sources at truth + off + jitter, magnitudes + mbias, and nfp false positives."""
    keep = rng.uniform(size=len(tx)) < comp
    x = tx[keep] + off[0] + rng.normal(0, jit, keep.sum())
    y = ty[keep] + off[1] + rng.normal(0, jit, keep.sum())
    m = tm[keep] + mbias + rng.normal(0, 0.1, keep.sum())
    fx, fy, fm = rng.uniform(3, 96, nfp), rng.uniform(3, 96, nfp), rng.uniform(20.0, 23.0, nfp)
    return np.r_[x, fx], np.r_[y, fy], A.mag_to_counts(np.r_[m, fm])


def _pack_row(arrays, pre, label, det, tp):
    for k, v in X._pack(label, A.score(np.c_[det[0], det[1]], det[2], tp)).items():
        arrays[pre + k] = v


def _planted_world(seed=0, capped=CAPPED):
    """(rep, base, truth): the stored repeat raw, the baselines' catalogues and the truth lists (R.truth_arrays layout) of six planted patches of three frames, self-consistent: every stored score is the
    score of the stored detections."""
    rng = np.random.default_rng(seed)
    rep, base, tr = {}, {}, {}
    for name, pxy in S.PATCHES:
        x0, y0 = pxy
        pre = E.pre_of(name, pxy)
        n = 90
        tx, ty, tm = rng.uniform(-0.4, 99.4, n), rng.uniform(-0.4, 99.4, n), rng.uniform(16.5, 23.0, n)
        tf = tr.setdefault(name, dict(px=[], py=[], m606=[]))
        tf['px'].append(tx + x0), tf['py'].append(ty + y0), tf['m606'].append(tm)
        off, general = OFFS[name], name in GENERAL
        is_cap = (name, pxy) in capped
        tp = A.truth_for_patch(dict(px=tx + x0, py=ty + y0, m606=tm), x0, y0)
        cats = {}
        for kern, extra in (('iso', 0.12 if general else 0.0), ('adq', 0.0), ('scl', 0.0)):
            for cfg, comp, nfp in (('ours', 0.75, 30), ('default', 0.5, 12)):
                if kern == 'adq' and not general:
                    cats[(kern, cfg)] = cats[('iso', cfg)]
                else:
                    cats[(kern, cfg)] = _catalogue(rng, tx, ty, tm, comp, 0.1, nfp, (off[0], off[1] + extra))
        for (kern, cfg), det in cats.items():
            flag = int(is_cap and cfg == 'default' and kern in ('iso', 'adq') and (kern == 'iso' or not general))
            rep[pre + '%s|%s|det' % (kern, cfg)] = np.c_[det[0], det[1], det[2]]
            rep[pre + '%s|%s|bic' % (kern, cfg)] = np.array([100, 2501, 1] if flag else [800, 1200, 0])
            rep[pre + '%s|%s|rsd' % (kern, cfg)] = np.array([1.1 + 0.2 * (kern == 'iso') + 0.5 * flag])
            rep[pre + '%s|%s|secs' % (kern, cfg)] = np.array([12.0, 2500.0 if flag else 1200.0])
            _pack_row(rep, pre, '%s|%s|%s' % (kern, cfg, 'k2' if cfg == 'ours' else 'k5'), det, tp)
            if cfg == 'default':
                _pack_row(rep, pre, '%s|default|k0' % kern, _catalogue(rng, tx, ty, tm, 0.52, 0.1, 12, (off[0], off[1] + extra)), tp)
        pcat = {}
        for tag, comp in (('photutils', 0.15), ('photutils_ext_all_local_i6', 0.35)):
            pcat[tag] = _catalogue(rng, tx, ty, tm, comp, 0.15, 6, off)
            for lab, v in zip(('x', 'y', 'flux'), pcat[tag]):
                base['x%d_y%d|%s|%s' % (x0, y0, tag, lab)] = v
        cs = _catalogue(rng, tx, ty, tm, 0.4, 0.15, 8, off)
        k = 'x%d_y%d|crowdsource|' % pxy
        base[k + 'x'], base[k + 'y'], base[k + 'flux'] = cs
        base[k + 'dflux'] = cs[2] / rng.uniform(2.0, 20.0, len(cs[2]))
        base[k + 'fracflux'], base[k + 'qf'], base[k + 'spread'] = np.ones(len(cs[2])), np.ones(len(cs[2])), np.zeros(len(cs[2]))
        rep[pre + 'crowdsource (tuned)|name'], rep[pre + 'photutils (tuned)|name'] = np.array([TUNED_CS]), np.array(['photutils_ext_all_local_i6'])
        for lab, det in (('crowdsource (default)', X.cs_dets(base, pxy, (0.0, 0.0, 0.0, np.inf))), ('crowdsource (tuned)', X.cs_dets(base, pxy, (5.0, 0.0, 0.0, np.inf))),
                         ('photutils (default)', pcat['photutils']), ('photutils (tuned)', pcat['photutils_ext_all_local_i6'])):
            _pack_row(rep, pre, lab, det, tp)
        rep[pre + 'kind'] = np.array(['general' if general else 'isotropic'])
        rep[pre + 'rms'] = np.array([0.005, 0.0015]) if general else np.array([0.0015, 0.0015])
        rep[pre + 'shape'] = np.array([0.09 if general else 0.01, -40.0, 1.4, 1.4])
        rep[pre + 'fwhm'] = np.array([3.0 if general else 2.3, 2.15])
        rep[pre + 'scale'] = np.array([1.4 if general else 1.06])
        rep[pre + 'kernel|same'] = np.array([0 if general else 1])
        rep[pre + 'cfg'] = np.array([3.0, 1.2, 2.0])
        rep[pre + 'parts5'] = np.ones((4, 5)) if general else np.zeros((0, 5))
    truth = {n: dict(px=np.concatenate(v['px']), py=np.concatenate(v['py']), m606=np.concatenate(v['m606'])) for n, v in tr.items()}
    return rep, base, truth


def _world_inputs(seed=0, capped=CAPPED):
    """The full `inp` of `build_raw` on the planted world (the star fits are planted arrays: their own tests are above)."""
    rep, base, truth = _planted_world(seed, capped)
    rng = np.random.default_rng(seed + 1)
    truth_g, cats = _synthetic_catalogs(rng)
    m2raw = {}
    for w in X.TESTS:
        m2raw['T|x%d_y%d|cfg' % w] = np.array([3.0, 1.2, 2.0])
        for lab, key in (('ours|k2', 'good'), ('crowdsource (tuned)', 'mid'), ('photutils (tuned)', 'poor')):
            for kk, v in X._pack(lab, A.score(np.c_[cats[key][0], cats[key][1]], cats[key][2], truth_g)).items():
                m2raw['T|x%d_y%d|%s' % (w[0], w[1], kk)] = v
    psfraw = {}
    for ds, plist in (('m2', X.TESTS), ('m53', X.SECOND['m53']['patches']), ('m92', X.SECOND['m92']['patches'])):
        for pxy in plist:
            k = 'S|%s|x%d_y%d|' % (ds, pxy[0], pxy[1])
            psfraw[k + 'rms'], psfraw[k + 'shape'] = (np.array([0.001, 0.001]), np.array([0.004, 20.0, 1.3, 1.3])) if ds == 'm2' else (np.array([0.007, 0.0015]), np.array([0.1, 30.0, 1.4, 1.5]))
            psfraw[k + 'kind'], psfraw[k + 'parts5'] = (np.array(['isotropic']), np.zeros((0, 5))) if ds == 'm2' else (np.array(['general']), np.ones((5, 5)))
    meas = {}
    for label, plist in S.all_labels():
        meas[label] = {}
        off = OFFS.get(label, (0.08, 0.02))
        for i, pxy in enumerate(plist):
            nb = 4 + i
            stars = {}
            for sel, k in (('A', 1 + i % 2), ('B', nb)):
                a = np.zeros((k, 6))
                a[:, 3] = off[0] + rng.normal(0, 0.03, k)
                a[:, 4] = off[1] + rng.normal(0, 0.03, k)
                stars[sel] = a
            meas[label][tuple(pxy)] = dict(stars=stars, nsel=dict(A=len(stars['A']), B=len(stars['B']) + 1), gain=4.6, dark=1.0, fwhm=2.2 + 0.1 * len(label))
    scan = {(n, p): OFFS[n] for n, p in S.PATCHES}
    scan_m2 = {('M2', tuple(w)): (0.08, 0.03) for w in X.TESTS}
    plain = {k: (v[0] + 0.01, v[1] - 0.01) for k, v in scan.items()}
    st = _gauss_stamp(1.3, 1.3)
    return dict(rep=rep, base=base, m2raw=m2raw, psfraw=psfraw, truth=truth, scan=scan, scan_m2=scan_m2, plain=plain, meas=meas, sha256={'repeat_raw.npz': 'ab' * 32, 'm2_raw.npz': 'cd' * 32},
                bias_psf=st, bias_gain=4.6, bias_dark=1.0)


@pytest.fixture(scope='module')
def world():
    inp = _world_inputs()
    return inp, S.build_raw(inp)


def test_the_planted_world_is_self_consistent_and_the_builder_runs(world):
    inp, raw = world
    assert len(raw) > 2000
    for name, pxy in S.PATCHES:
        pre = E.pre_of(name, pxy)
        for key in S.STORED_ROWS:
            assert np.array_equal(raw['st|%s|%s|tot' % (S.pk(name, pxy), key)], inp['rep'][pre + key + '|tot'])
            assert np.array_equal(raw['st|%s|%s|blk22' % (S.pk(name, pxy), key)], inp['rep'][pre + key + '|blk22'])
            assert np.array_equal(raw['st|%s|%s|offset' % (S.pk(name, pxy), key)], inp['rep'][pre + key + '|offset'], equal_nan=True)
            bright = [inp['rep'][pre + key + '|n_matched'][:9].sum(), inp['rep'][pre + key + '|n_true'][:9].sum()]
            assert raw['st|%s|%s|bright' % (S.pk(name, pxy), key)].tolist() == bright and bright[1] < inp['rep'][pre + key + '|n_true'][:10].sum()      # the first nine bins: 15.0 to 19.5 mag
        for key in S.SC_KEYS + S.BASE_EXTRA:                                   # the base scenario is the stored scoring
            assert np.array_equal(raw['sc|base|%s|%s|tot' % (S.pk(name, pxy), key)][ICUT], inp['rep'][pre + key + '|tot'])
            assert np.array_equal(raw['sc|base|%s|%s|blk' % (S.pk(name, pxy), key)][ICUT], inp['rep'][pre + key + '|blk22'])
    assert S.general_patches(raw) == [('m2r2', (133, 281)), ('m2r2', (195, 689))] == S.GEN_PATCHES
    # the unit-pair differences of the two kernels exist for the patches of the general kernel only, and are chosen minus isotropic (the planted isotropic run sits 0.12 px further in y)
    assert sorted(k for k in raw if k.startswith('pos|')) == sorted('pos|%s|m%d' % (S.pk(n, p), m) for n, p in S.GEN_PATCHES for m in (19, 20, 21, 22))
    for n, p in S.GEN_PATCHES:
        pre = E.pre_of(n, p)
        for m in (19, 20, 21, 22):
            v = raw['pos|%s|m%d' % (S.pk(n, p), m)]
            assert np.array_equal(v, S.unit_pair_difference(inp['rep'][pre + 'adq|ours|det'], inp['rep'][pre + 'iso|ours|det'], float(m))) and v[0] > 5
        assert -0.2 < raw['pos|%s|m22' % S.pk(n, p)][2] < -0.01 and abs(raw['pos|%s|m22' % S.pk(n, p)][1]) < 0.1                 # chosen minus isotropic: the planted isotropic run is further in +y


# ---- the builder: scenarios, shifts, offsets, defective flags ----------------------------------------------------------------------------

def test_the_builder_stores_every_scenario_row_with_the_declared_shapes(world):
    inp, raw = world
    for scn, spec in S.SCENARIOS.items():
        for name, pxy in S.PATCHES:
            p = S.pk(name, pxy)
            assert raw['shift|%s|%s' % (scn, p)].shape == (2,)
            for key in S.SC_KEYS:
                pre = 'sc|%s|%s|%s|' % (scn, p, key)
                assert raw[pre + 'tot'].shape == (3, 7) and raw[pre + 'blk'].shape == (3, 25, 4) and raw[pre + 'off'].shape == (2,) and raw[pre + 'moff'].shape == (1,)
                assert (pre + 'blk10' in raw) == (scn == 'base' and key in S.ISS_KEYS)
            for key in S.BASE_EXTRA:
                assert ('sc|%s|%s|%s|tot' % (scn, p, key) in raw) == (scn == 'base')
    p = S.pk('m2r2', (133, 281))
    assert raw['sc|base|%s|iso|ours|k2|blk10' % p].shape == (3, 100, 4) and raw['sc|base|%s|adq|ours|k2|blk25' % p].shape == (3, 16, 4) and raw['sc|base|%s|adq|default|k5|blk50' % p].shape == (3, 4, 4)
    assert 'sc|base|%s|crowdsource (tuned)|blk10' % p not in raw
    # the tolerances and the cuts of each scenario
    for scn in ('tol35', 'mag30', 'tol75', 'base'):
        spec = S.SCENARIOS[scn]
        r = S.rescore(S.catalogues(inp['rep'], inp['base'], 'm2r1', (316, 124))['iso|ours|k2'], A.truth_for_patch(inp['truth']['m2r1'], 316, 124), spec['dpos'], spec['dmag'])
        assert np.array_equal(raw['sc|%s|%s|iso|ours|k2|tot' % (scn, S.pk('m2r1', (316, 124)))], r['tot'])


def test_the_tolerance_scenarios_order_the_matched_counts(world):
    inp, raw = world
    tot = {scn: sum(raw['sc|%s|%s|%s|tot' % (scn, S.pk(n, p), 'crowdsource (tuned)')][ICUT][5] for n, p in S.PATCHES) for scn in ('tol35', 'base', 'tol75', 'mag30')}
    assert tot['tol35'] < tot['base'] < tot['tol75'] and tot['mag30'] < tot['base']


def test_the_shift_of_each_scenario_is_alpha_times_the_offset_of_its_source(world):
    inp, raw = world
    own, ownp = S.own_offsets(inp['meas'])
    for name, pxy in S.PATCHES:
        p = S.pk(name, pxy)
        stars = {q: inp['meas'][name][q]['stars']['B'] for q in inp['meas'][name]}
        fm = np.median(np.concatenate(list(stars.values()))[:, 3:5], axis=0)
        assert np.allclose(own[name], fm) and np.allclose(ownp[(name, pxy)], np.median(stars[pxy][:, 3:5], axis=0) if len(stars[pxy]) >= 3 else fm)
        for scn, expected in (('base', (0.0, 0.0)), ('tol35', (0.0, 0.0)), ('mo', (0.0, 0.0)), ('scan25', 0.25 * np.array(OFFS[name])), ('scan50', 0.5 * np.array(OFFS[name])), ('scan75', 0.75 * np.array(OFFS[name])),
                              ('scan100', np.array(OFFS[name])), ('mo+scan100', np.array(OFFS[name])), ('own25', 0.25 * fm), ('own50', 0.5 * fm), ('own75', 0.75 * fm), ('own100', fm), ('mo+own100', fm),
                              ('ownp100', ownp[(name, pxy)])):
            assert np.allclose(raw['shift|%s|%s' % (scn, p)], expected, atol=0, rtol=1e-14), (scn, name)
    # the M2 windows' offsets and the plain scan are stored for the report
    assert np.allclose(raw['scan|M2|x310_y630'], (0.08, 0.03)) and np.allclose(raw['plain|m2r2|x133_y281'], (OFFS['m2r2'][0] + 0.01, OFFS['m2r2'][1] - 0.01))
    # the sign: the shifts of the planted frames are about the planted offsets, which are the offset of the image against the truth (added to the truth)
    assert raw['shift|scan100|%s' % S.pk('m2r2', (133, 281))][1] == 0.40


def test_moving_the_truth_by_the_planted_offset_raises_the_scores_and_removes_the_median_offset(world):
    inp, raw = world
    gen = S.GEN_PATCHES
    f = [S.F_of(raw, gen, 'iso|ours|k2', S.scn_name('scan', a)) for a in S.ALPHAS]
    assert all(b > a for a, b in zip(f, f[1:]))                                                       # the isotropic run sits furthest from the unshifted truth
    d = [S.dF_stats(raw, gen, 'adq|ours|k2', 'iso|ours|k2', S.scn_name('scan', a))[0] for a in S.ALPHAS]
    assert d[0] > d[-1] and d[0] > 0
    for n, p in gen:
        o0 = raw['sc|base|%s|iso|ours|k2|off' % S.pk(n, p)]
        o1 = raw['sc|scan100|%s|iso|ours|k2|off' % S.pk(n, p)]
        assert 0.3 < o0[1] < 0.5 and abs(o1[1] - 0.12) < 0.06       # the planted 0.40 + 0.12 px in y: only the matches within 0.5 px are in the median (truncated); then 0.12 px
        assert abs(raw['sc|scan100|%s|adq|ours|k2|off' % S.pk(n, p)][1]) < 0.06
    # the own source (the planted star offsets) does the same
    assert S.F_of(raw, gen, 'iso|ours|k2', 'own100') > S.F_of(raw, gen, 'iso|ours|k2', 'base')


def test_the_magnitude_offset_scenarios_store_the_offset_used(world):
    inp, raw = world
    name, pxy = 'm2r1', (316, 124)
    cats = S.catalogues(inp['rep'], inp['base'], name, pxy)
    tp = A.truth_for_patch(inp['truth'][name], *pxy)
    for key in S.SC_KEYS:
        r = S.rescore(cats[key], tp, moff='auto')
        assert np.array_equal(raw['sc|mo|%s|%s|moff' % (S.pk(name, pxy), key)], [r['moff']]) and np.array_equal(raw['sc|mo|%s|%s|tot' % (S.pk(name, pxy), key)], r['tot'])
    assert np.array_equal(raw['sc|base|%s|%s|moff' % (S.pk(name, pxy), key)], [0.0])


def test_the_builder_refuses_a_base_row_that_differs_from_the_stored_one_and_a_frame_without_stars():
    inp = _world_inputs()
    bad = dict(inp['rep'])
    k = E.pre_of('m2r1', (316, 124)) + 'adq|ours|k2|tot'
    bad[k] = bad[k] + 1
    with pytest.raises(RuntimeError, match='differs from the stored row'):
        S.build_scenarios(bad, inp['base'], inp['truth'], inp['scan'], {n: (0.0, 0.0) for n in S.FRAMES}, {k_: (0.0, 0.0) for k_ in S.PATCHES}, {})
    bad = dict(inp['rep'])
    k = E.pre_of('m2r2', (195, 689)) + 'photutils (tuned)|blk22'
    bad[k] = bad[k].copy()
    bad[k][0, 0] += 1
    with pytest.raises(RuntimeError, match='differs from the stored row'):
        S.build_scenarios(bad, inp['base'], inp['truth'], inp['scan'], {n: (0.0, 0.0) for n in S.FRAMES}, {k_: (0.0, 0.0) for k_ in S.PATCHES}, {})
    bad = dict(inp['rep'])
    k = E.pre_of('m2r1', (316, 224)) + 'iso|default|k5|offset'
    bad[k] = bad[k] + 0.001
    with pytest.raises(RuntimeError, match='differs from the stored row'):
        S.build_scenarios(bad, inp['base'], inp['truth'], inp['scan'], {n: (0.0, 0.0) for n in S.FRAMES}, {k_: (0.0, 0.0) for k_ in S.PATCHES}, {})
    few = dict(inp)
    few['meas'] = {k_: dict(v) for k_, v in inp['meas'].items()}
    few['meas']['m2r3'] = {p: dict(v, stars={'A': v['stars']['A'][:1], 'B': v['stars']['B'][:1]}) for p, v in inp['meas']['m2r3'].items()}
    with pytest.raises(RuntimeError, match='no own offset'):
        S.build_raw(few)


def test_the_defective_flags_are_read_from_the_stored_bic_of_each_run(world):
    inp, raw = world
    assert S.capped(raw, 'm2r2', (133, 281), 'iso|default|k5') and S.capped(raw, 'm2r2', (133, 281), 'iso|default|k0')
    assert not S.capped(raw, 'm2r2', (133, 281), 'adq|default|k5') and not S.capped(raw, 'm2r2', (195, 689), 'iso|default|k5')
    assert S.capped(raw, 'm2r3', (370, 441), 'adq|default|k5') and S.capped(raw, 'm2r3', (110, 501), 'adq|default|k0')                      # the chosen kernel is the isotropic run there
    for key in ('iso|ours|k2', 'adq|ours|k2', 'scl|default|k5', 'crowdsource (tuned)', 'photutils (default)'):
        assert not any(S.capped(raw, n, p, key) for n, p in S.PATCHES)
    assert S.defective(raw, S.PATCHES, ['iso|default|k5', 'adq|default|k5']) == [('m2r2', (133, 281)), ('m2r3', (370, 441)), ('m2r3', (110, 501))]
    assert S.defective(raw, S.PATCHES, ['adq|default|k5']) == [('m2r3', (370, 441)), ('m2r3', (110, 501))]
    assert S.clean(raw, S.PATCHES, ['iso|default|k5']) == [('m2r1', (316, 124)), ('m2r1', (316, 224)), ('m2r2', (195, 689))]
    assert S.defective(raw, S.PATCHES, ['adq|ours|k2', 'crowdsource (tuned)']) == [] and S.clean(raw, S.GEN_PATCHES, ['adq|default|k0', 'iso|default|k0']) == [('m2r2', (195, 689))]
    assert S.entrants(raw, 'm2r2', (133, 281), 'iso|default|k5') == 2500 and S.entrants(raw, 'm2r2', (133, 281), 'adq|default|k5') == 1200


def test_blocks_reads_the_scenario_rows_and_the_stored_bic_row(world):
    inp, raw = world
    b, sizes = S.blocks(raw, S.GEN_PATCHES, 'adq|ours|k2', 'scan50', 21.0)
    assert b.shape == (50, 4) and sizes == [25, 25]
    assert np.array_equal(b[:25], raw['sc|scan50|%s|adq|ours|k2|blk' % S.pk(*S.GEN_PATCHES[0])][1])
    b10, s10 = S.blocks(raw, S.GEN_PATCHES, 'iso|ours|k2', cell=10)
    assert b10.shape == (200, 4) and s10 == [100, 100] and b10.sum(0).tolist() == S.total(raw, S.GEN_PATCHES, 'iso|ours|k2').tolist()
    k0, sk = S.blocks(raw, S.PATCHES, 'adq|default|k0')
    assert k0.shape == (150, 4) and np.array_equal(k0[:25], raw['st|%s|adq|default|k0|blk22' % S.pk(*S.PATCHES[0])])
    with pytest.raises(AssertionError):
        S.blocks(raw, S.PATCHES, 'adq|default|k0', 'scan50')
    with pytest.raises(AssertionError):
        S.blocks(raw, S.PATCHES, 'adq|default|k0', cut=20.0)
    ba, bb = S.blocks(raw, S.GEN_PATCHES, 'adq|ours|k2')[0], S.blocks(raw, S.GEN_PATCHES, 'iso|ours|k2')[0]
    d, ci = S.dF_stats(raw, S.GEN_PATCHES, 'adq|ours|k2', 'iso|ours|k2')
    ref = X.boot_diff(ba, bb, X.boot_idx_strat([25, 25]))[2]
    assert d == ref[0] and np.array_equal(ci, ref[1])
    assert S.F_of(raw, S.GEN_PATCHES, 'adq|ours|k2') == X.cpf(S.total(raw, S.GEN_PATCHES, 'adq|ours|k2'))[2]


# ---- the report on the planted raw -------------------------------------------------------------------------------------------------------

@pytest.fixture()
def small(monkeypatch):
    """The 200 seeds and the 40000 resamples shortened (the constants themselves are pinned above)."""
    monkeypatch.setattr(S, 'SEEDS', (1, 2, 3))
    monkeypatch.setattr(S, 'BIG_B', 3000)


@pytest.fixture(scope='module')
def md_small(world):
    """The report of the planted raw with the 200 seeds and the 40000 resamples shortened, computed once."""
    inp, raw = world
    with pytest.MonkeyPatch.context() as mp:                                  # the shortened constants act during the call only
        mp.setattr(S, 'SEEDS', (1, 2, 3))
        mp.setattr(S, 'BIG_B', 3000)
        return S.report(raw)


def _rows(md, start):
    """The lines of md that start with `start`, split into cells."""
    return [[c.strip() for c in l.strip().strip('|').split(' | ')] for l in md.split('\n') if l.startswith(start)]


def test_the_report_has_every_section_and_no_nan(world, md_small):
    inp, raw = world
    md = md_small
    for needle in ('# Post hoc sensitivity analyses of the held-out repeat-frame test', '**POST HOC.**', 'all 60 stored rows', '## 1. Truth-to-image offsets', '## 2. (a) P1 against the truth offset',
                   '### 2.1 Tuned configuration (kappa 2), source: scan offsets of each patch', '### 2.2 Tuned configuration (kappa 2), source: own: frame median of the star fits',
                   '### 2.3 Intervals of the pooled dF', '### 2.4 The a priori default and the BIC stop, with and without the defective path', '### 2.5 With the per-method magnitude offsets removed',
                   '## 3. (b) P1 robustness: block size', '## 4. (c) P2 robustness', '### 4.1 (i) Match tolerances', '### 4.2 (ii) Depth cuts', '### 4.3 (iii) The truth moved by alpha', '### 4.4 (iv)',
                   '### 4.5 Every cell', 'Of 408 (setting, comparison, patch) cells', '## 5. (d) Scaled lengths by frame', '## 6. (e) F ratios and the seeing ladder', '## 7. (f) Plan items missing from `repeat.md`',
                   '### 7.1 The default and BIC readouts', '### 7.2 P3', '### 7.3 Stars brighter than 19.5 mag matched', '## 8. Defective (capped) paths', '## 9. Frames, camera columns and stamp fits',
                   '## 10. Positions of the two kernels on m2r2', '## 11. Limits, deviations, what was not done', '**What section 2 shows.**', '**What it does not show.**', 'DILUTED',
                   'repeat_validation.json[<frame>]["<x0>,<y0>"]["masked"]["asTrans"]', '| `repeat_raw.npz` |'):
        assert needle in md, needle
    assert 'nan' not in md.lower().replace('nan,', '')
    assert md.endswith('\n') and not md.endswith('\n\n')


def test_the_alpha_table_reproduces_the_stored_p1_and_shows_the_shifts(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('### 2.1 ')[1].split('### 2.2 ')[0]
    rows = _rows(sec, '| ')
    rows = [r for r in rows if r[0] in ('0', '0.25', '0.5', '0.75', '1')]
    assert [r[0] for r in rows] == ['0', '0.25', '0.5', '0.75', '1']
    gen = S.GEN_PATCHES
    ba = np.concatenate([inp['rep'][E.pre_of(n, p) + 'adq|ours|k2|blk22'] for n, p in gen])
    bb = np.concatenate([inp['rep'][E.pre_of(n, p) + 'iso|ours|k2|blk22'] for n, p in gen])
    d, ci = X.boot_diff(ba, bb, X.boot_idx_strat([25, 25]))[2]
    assert rows[0][1] == '(+0.000, +0.000) \\| (+0.000, +0.000)' and rows[0][4] == X.fmt_iv(d, ci) and rows[0][5] == X.verdict(d, ci)
    assert rows[0][2] == '%.3f' % X.cpf(bb.sum(0))[2] and rows[0][3] == '%.3f' % X.cpf(ba.sum(0))[2]
    assert rows[4][1] == '(-0.050, +0.400) \\| (-0.050, +0.400)' and rows[2][1] == '(-0.025, +0.200) \\| (-0.025, +0.200)'
    # per patch columns are the intervals of one patch's blocks
    d1, ci1 = X.boot_diff(inp['rep'][E.pre_of(*gen[0]) + 'adq|ours|k2|blk22'], inp['rep'][E.pre_of(*gen[0]) + 'iso|ours|k2|blk22'], X.boot_idx(25))[2]
    assert rows[0][6] == X.fmt_iv(d1, ci1)
    own = md.split('### 2.2 ')[1].split('### 2.3 ')[0]
    assert '| 1 (per-patch medians) |' in own and own.count('| 1 |') == 1
    # the shifts of the own scenario are the frame median of the planted star offsets (about the planted 0.40)
    r_own = [r for r in _rows(own, '| ') if r[0] == '1'][0]
    fm = S.own_offsets(inp['meas'])[0]['m2r2']
    assert r_own[1] == ' \\| '.join(['(%+.3f, %+.3f)' % fm] * 2) and 0.35 < fm[1] < 0.45


def test_defective_runs_are_marked_and_every_pooled_default_row_is_given_without_them(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('### 2.3 ')[1].split('### 2.4 ')[0]
    rows = [r for r in _rows(sec, '| ') if len(r) == 8 and r[0] not in ('readout', '---')]
    tuned = [r for r in rows if r[0].startswith('tuned')]
    default = [r for r in rows if r[0].startswith('a priori default (')]
    bic = [r for r in rows if r[0].startswith('a priori default with the BIC')]
    assert len(tuned) == 10 and all(r[2] == 'no' for r in tuned)
    assert len(default) == 11 and len(bic) == 2                                                           # the readout at alpha 0 is also given without the defective patch
    assert sum(r[2] == 'yes' for r in default) == 10 and sum(r[2] == 'no' for r in default) == 1 and 'without the defective patch (m2r2 (195, 689))' in [r[1] for r in default][1] + [r[1] for r in default][0]
    assert sorted(r[2] for r in bic) == ['no', 'yes']
    sec = md.split('### 2.4 ')[1].split('### 2.5 ')[0]
    rows = [r for r in _rows(sec, '| a priori') if len(r) == 9]
    assert len(rows) == 3 * 2 + 3 and [r[3] for r in rows[:3]] == ['yes', 'no', 'yes'] and [r[3] for r in rows[-3:]] == ['yes', 'no', 'yes']
    assert rows[1][2] == 'without the defective patch (m2r2 (195, 689))' + '' or rows[1][2].startswith('without the defective patch (m2r2 (195, 689))')
    # section 7.1: a D marks the patches whose compared run is capped, and the pooled row is given without them
    sec = md.split('### 7.1 ')[1].split('By frame')[0]
    rows = {r[0]: r for r in _rows(sec, '| ') if len(r) == 13}
    iso = rows['isotropic kernel, a priori default vs photutils (default)']
    chosen = rows['chosen kernel, a priori default vs photutils (default)']
    assert [c.endswith(' D') for c in iso[1:7]] == [False, False, True, False, True, True] and [c.endswith(' D') for c in chosen[1:7]] == [False, False, False, False, True, True]
    assert iso[9] == 'yes' and iso[10] == '3' and chosen[10] == '4'
    cl = S.clean(raw, S.PATCHES, ['iso|default|k5', 'photutils (default)'])
    d, ci = S.dF_stats(raw, cl, 'iso|default|k5', 'photutils (default)')
    assert iso[11] == X.fmt_iv(d, ci)
    # sections 5 and 7.3: the pooled rows without the defective runs
    sc = md.split('## 5. (d)')[1].split('## 6. (e)')[0]
    assert sc.count('all patches without a defective run') == 2 and '| tuned configuration (kappa 2) | all six patches | no |' in sc and '| a priori default (kappa 5) | m2r3 | yes |' in sc
    br = md.split('### 7.3 ')[1].split('**What section 7')[0]
    assert 'the patches without a defective run (their number in brackets)' in br and ' D |' in br
    assert 'Flagged runs: m2r2 (133, 281) iso default; m2r3 (370, 441) iso default; m2r3 (370, 441) adq default; m2r3 (110, 501) iso default; m2r3 (110, 501) adq default.' in md


def test_the_p3_verdict_per_readout_follows_the_plans_rule(world, small):
    inp, raw = world
    md = S.report(raw)
    sec = md.split('### 7.2 ')[1].split('### 7.3 ')[0]
    verd = [r for r in _rows(sec, '| ') if len(r) == 5 and r[0] in ('tuned configuration (kappa 2)', 'a priori default (kappa 5)')]
    assert verd == [['tuned configuration (kappa 2)', '2 of 2', 'supported', '2 of 2', 'supported'], ['a priori default (kappa 5)', '2 of 2', 'supported', '1 of 1', 'supported']]
    pr = [r for r in _rows(sec, '| ') if len(r) == 7 and r[1].startswith('m2r2')]
    assert [r[4] for r in pr] == ['yes'] * 4 and [r[5] for r in pr] == ['no', 'no', 'yes', 'no'] and pr[0][2] == '1.30' and pr[0][3] == '1.10'
    # flip the planted residual scales of the chosen kernel above the isotropic one on one patch: not supported
    raw2 = dict(raw)
    raw2['st|%s|adq|ours|rsd' % S.pk('m2r2', (195, 689))] = np.array([1.9])
    sec2 = S.report(raw2).split('### 7.2 ')[1].split('### 7.3 ')[0]
    verd2 = [r for r in _rows(sec2, '| ') if len(r) == 5 and r[0] == 'tuned configuration (kappa 2)']
    assert verd2 == [['tuned configuration (kappa 2)', '1 of 2', 'not supported', '1 of 2', 'not supported']]


def test_the_p2_tables_record_the_wins_the_depth_cut_and_the_exceptions(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('### 4.1 ')[1].split('F of the catalogues under the tolerances')[0]
    rows = [r for r in _rows(sec, '| ') if len(r) == 11 and r[1].endswith('(tuned)')]
    assert len(rows) == 16 and [r[0] for r in rows[:4]] == ['0.5 px / 0.5 mag (stored)'] * 4 and rows[4][0] == '0.35 px / 0.5 mag' and rows[8][0] == '0.5 px / 0.3 mag' and rows[12][0] == '0.75 px / 0.75 mag'
    first = rows[0]
    per, pooled, v = S.pair_cells(raw, 'base', 'adq|ours|k2', 'crowdsource (tuned)')
    assert first[1] == 'chosen kernel vs crowdsource (tuned)' and first[8] == X.fmt_iv(*pooled) and first[9] == X.verdict(*pooled) and first[10] == '%d / %d / %d' % (v.count('win'), v.count('tie'), v.count('loss'))
    assert first[2:8] == [S.vcell(d, ci) for d, ci in per]
    # the pooled difference is read off the pooled blocks of the six patches (stratified), with the pair taken first minus second
    ba = np.concatenate([inp['rep'][E.pre_of(n, p) + 'adq|ours|k2|blk22'] for n, p in S.PATCHES])
    bb = np.concatenate([inp['rep'][E.pre_of(n, p) + 'crowdsource (tuned)|blk22'] for n, p in S.PATCHES])
    assert first[8] == X.fmt_iv(*X.boot_diff(ba, bb, X.boot_idx_strat([25] * 6))[2])
    dep = md.split('### 4.2 ')[1].split('F of the catalogues at the depth cuts')[0]
    drows = [r for r in _rows(dep, '| ') if len(r) == 11 and r[1].endswith('(tuned)')]
    assert [r[0] for r in drows] == ['m < 20'] * 4 + ['m < 21'] * 4 + ['m < 22'] * 4
    per20, pooled20, _ = S.pair_cells(raw, 'base', 'adq|ours|k2', 'crowdsource (tuned)', 20.0)
    assert drows[0][8] == X.fmt_iv(*pooled20) and drows[8][8] == first[8]
    # F per patch and method at the cuts
    f20 = S.F_of(raw, [S.PATCHES[2]], 'crowdsource (tuned)', 'base', 20.0)
    assert '| m < 20 | crowdsource (tuned) | ' in md and ('%.3f' % f20) in md.split('F of the catalogues at the depth cuts')[1]
    # section 4.5 lists exactly the cells that are not a win
    exc = md.split('### 4.5 ')[1].split('**What section 4 shows.**')[0]
    n_not = int(exc.split('cells ')[1].split(' are not a win')[0])
    assert n_not == len([l for l in exc.split('\n') if l.startswith('| ') and not l.startswith('| setting') and '---' not in l])


def test_the_ratio_table_and_the_ladder_are_read_from_the_stored_rows(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('## 6. (e)')[1].split('The seeing ladder')[0]
    rows = {r[0]: r for r in _rows(sec, '| ') if len(r) == 9}
    r1 = rows['m2r1']
    pl = [(n, p) for n, p in S.PATCHES if n == 'm2r1']
    ba, sizes = S.blocks(raw, pl, 'adq|ours|k2')
    bc, _ = S.blocks(raw, pl, 'crowdsource (tuned)')
    r, ci = S.boot_ratio(ba, bc, sizes)
    assert r1[4] == X.fmt_iv(r, ci, 2) and r1[2] == '%.3f' % X.cpf(ba.sum(0))[2] and r1[3] == '%.3f' % X.cpf(bc.sum(0))[2]
    # the M2 row reads m2_raw (the frozen rows), the M2-window column pools the windows the frame repeats
    m2 = rows['M2 (4 test windows)']
    parts = [inp['m2raw']['T|x%d_y%d|ours|k2|blk22' % w] for w in X.TESTS]
    cparts = [inp['m2raw']['T|x%d_y%d|crowdsource (tuned)|blk22' % w] for w in X.TESTS]
    r, ci = S.boot_ratio(np.concatenate(parts), np.concatenate(cparts), [25] * 4)
    assert m2[4] == X.fmt_iv(r, ci, 2) and m2[7] == '-'
    wins = [tuple(R.REPEAT['m2r1']['windows'][i]) for i in range(2)]
    mp = np.concatenate([inp['m2raw']['T|x%d_y%d|ours|k2|blk22' % w] for w in wins])
    mc = np.concatenate([inp['m2raw']['T|x%d_y%d|crowdsource (tuned)|blk22' % w] for w in wins])
    r, ci = S.boot_ratio(mp, mc, [25, 25])
    assert r1[7] == X.fmt_iv(r, ci, 2)
    lad = md.split('The seeing ladder')[1].split('**What section 6')[0]
    l630 = [l for l in lad.split('\n') if l.startswith('| (310, 630) |')][0]
    assert l630.count('ISS') == 3 and l630.count('M2 (') == 1 and 'm2r1 (' in l630 and 'm2r2 (' in l630                         # (310, 630) is in M2, m2r1 and m2r2
    assert [l for l in lad.split('\n') if l.startswith('| (310, 730) |')][0].count('ISS') == 2


def test_the_confound_table_lists_the_camera_columns_and_the_positions_table_the_unit_pairs(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('## 9. Frames')[1].split('## 10.')[0]
    assert '| m2r2 | 7778-6-72 | 6 | (133, 281) |' in sec and '| M2 | 2583-2-136 | 2 | (310, 630) |' in sec and '| M53 | 5390-6-162 | 6 | (760, 590) |' in sec and '| M92 | 4682-6-208 | 6 | (1100, 220) |' in sec
    assert '| m2r1 | 4187-2-105 | 2 | (316, 124) |' in sec and '| m2r3 | 3388-2-75 | 2 | (110, 501) |' in sec
    r1 = [l for l in sec.split('\n') if l.startswith('| m2r1 | 4187-2-105 | 2 | (316, 124) |')][0].split(' | ')
    assert r1[9:12] == ['isotropic', '3', '0.00150 |'] and r1[7] == '0.00150' and r1[8] == '0.010'                      # the isotropic kernel has three parts
    r2 = [l for l in sec.split('\n') if l.startswith('| m2r2 | 7778-6-72 | 6 | (133, 281) |')][0].split(' | ')
    assert r2[7] == '0.00500' and r2[9:12] == ['general', '4', '0.00150 |'] and r2[8] == '0.090'
    m2 = [l for l in sec.split('\n') if l.startswith('| M2 | 2583-2-136 | 2 | (310, 630) |')][0].split(' | ')
    assert m2[9:11] == ['isotropic', '3'] and m2[7] == '0.00100'
    m53 = [l for l in sec.split('\n') if l.startswith('| M53 | 5390-6-162 | 6 | (760, 590) |')][0].split(' | ')
    assert m53[9:12] == ['general', '5', '0.00150 |'] and m53[7] == '0.00700'
    assert 'where the rule picks a general kernel on at least one patch: m2r2, M53, M92 (camera columns 6, 6, 6); where it keeps the isotropic kernel on every patch: m2r1, m2r3, M2 (camera columns 2, 2, 2)' in sec
    row = [l for l in sec.split('\n') if l.startswith('| m2r2 | 7778-6-72 | 6 | (195, 689) |')][0].split(' | ')
    assert row[4] == '4.600' or row[4] == '4.600'
    assert S.frame_id('m2r2') == (7778, 6, 72) and S.frame_id('M2') == (2583, 2, 136) and S.frame_id('M53') == (5390, 6, 162) and S.frame_id('M92') == (4682, 6, 208)
    pos = md.split('## 10.')[1].split('**What section 10')[0]
    for n, p in S.GEN_PATCHES:
        assert '| m2r2 (%d, %d) |' % p in pos
    assert pos.count('| m2r2 (') == 2


def test_the_offsets_table_prints_both_selections_and_the_scan_key(world, md_small):
    inp, raw = world
    md = md_small
    sec = md.split('## 1. Truth-to-image offsets')[1].split('## 2. (a)')[0]
    rows = {(r[0], r[1]): r for r in _rows(sec, '| ') if len(r) == 9}
    r = rows[('m2r2', '(133, 281)')]
    assert r[2] == '(-0.050, +0.400)' and r[3] == '(-0.040, +0.390)'
    stars = inp['meas']['m2r2'][(133, 281)]['stars']
    assert r[4] == '%d (%d)' % (len(stars['A']), inp['meas']['m2r2'][(133, 281)]['nsel']['A']) and r[6] == '%d (%d)' % (len(stars['B']), inp['meas']['m2r2'][(133, 281)]['nsel']['B'])
    assert r[7] == S.xy(np.median(stars['B'][:, 3:5], axis=0)) and rows[('M2', '(310, 630)')][2] == '(+0.080, +0.030)' and rows[('M53', '(760, 590)')][2] == '-'
    frames = {r[0]: r for r in _rows(sec, '| ') if len(r) == 9 and r[0] in ('m2r2', 'M92') and r[1].isdigit()}
    allB = np.concatenate([inp['meas']['m2r2'][q]['stars']['B'] for q in inp['meas']['m2r2']])
    m, c = S.median_interval(allB[:, 3])
    assert frames['m2r2'][5] == '%+.3f [%+.3f, %+.3f]' % (m, c[0], c[1]) and frames['m2r2'][7] == '(-0.050, +0.400)' and frames['M92'][7] == '-'
    al = [l for l in sec.split('\n') if l.startswith('| m2r2 (133, 281) | ISS isotropic, tuned |')][0]
    assert al.count('(') >= 4 and '| bias' not in al
    assert '| 5 px | 0.' in sec and '| distance | flux 0.063 |' in sec


# ---- the command line --------------------------------------------------------------------------------------------------------------------

def test_from_raw_regenerates_the_md_byte_for_byte_in_a_fresh_process(world, small, tmp_path):
    inp, raw = world
    np.savez_compressed(str(tmp_path / 'repeat_sens_raw.npz'), **raw)
    code = ("import sys; sys.path.insert(0, %r); sys.path.insert(0, %r); import exp_astro_repeat_sens as S; S.SEEDS = (1, 2, 3); S.BIG_B = 3000; "
            "sys.argv = ['exp_astro_repeat_sens.py', '--from-raw']; S.main()") % (os.path.join(HERE, '..'), os.path.join(HERE, '..', 'experiments'))
    env = dict(os.environ, REPEAT_SENS_OUT=str(tmp_path), OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='7')
    p = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, cwd=str(tmp_path))
    assert p.returncode == 0 and 'wrote' in p.stdout, p.stderr
    out = (tmp_path / 'repeat_sens.md').read_bytes()
    assert out.decode() == S.report(X.load_raw(str(tmp_path / 'repeat_sens_raw.npz'))) and out.endswith(b'\n')                   # the in-process report has another hash seed
    p = subprocess.run([sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_astro_repeat_sens.py')], env=env, capture_output=True, text=True, cwd=str(tmp_path))
    assert p.returncode != 0 and '--run' in p.stderr                                                         # no option: the parser asks for one


def test_the_interval_table_uses_the_declared_resamples_seeds_and_the_pooled_blocks(world, md_small):
    inp, raw = world
    sec = md_small.split('### 2.3 ')[1].split('### 2.4 ')[0]
    rows = [r for r in _rows(sec, '| ') if len(r) == 8 and r[0].startswith('tuned')]
    r0, r1 = rows[0], rows[4]
    gen = S.GEN_PATCHES
    for row, scn in ((r0, 'base'), (r1, 'scan100')):
        ba, sizes = S.blocks(raw, gen, 'adq|ours|k2', scn)
        bb, _ = S.blocks(raw, gen, 'iso|ours|k2', scn)
        d, ci = S.boot_dF(ba, bb, sizes)
        ci3 = S.boot_dF(ba, bb, sizes, 3000, 0)[1]
        lows = S.seed_lowers(ba, bb, sizes, (1, 2, 3))
        assert row[3] == '%.4f' % d and row[4] == '[%.4f, %.4f]' % tuple(ci) and row[5] == '[%.4f, %.4f]' % tuple(ci3) and row[4] != row[5]
        assert row[6] == '%d of 3' % (lows > 0).sum() and row[7] == '%.4f (%.4f, %.4f)' % (lows.mean(), lows.min(), lows.max())
    assert 'seeds 1 to 3 (2000 resamples each)' in sec and 'against 3000 resamples (seed 0)' in sec


def test_seed_lowers_defaults_to_the_declared_seeds(monkeypatch):
    a, b = _blocks(11), _blocks(12)
    monkeypatch.setattr(S, 'SEEDS', (4, 9, 11, 12, 13, 14, 15))
    low = S.seed_lowers(a, b, [25, 25], B=300)
    assert len(low) == 7 and low[0] == S.boot_dF(a, b, [25, 25], 300, 4)[1][0] and low[1] == S.boot_dF(a, b, [25, 25], 300, 9)[1][0] and low[6] == S.boot_dF(a, b, [25, 25], 300, 15)[1][0]


def test_run_writes_the_raw_that_the_builder_returns(world, tmp_path, monkeypatch):
    inp, raw = world
    seen = []
    monkeypatch.setattr(S, 'load_inputs', lambda cyclefix=None: (seen.append(cyclefix), inp)[1])
    out = str(tmp_path / 'sub' / 'repeat_sens_raw.npz')
    got = S.run(out)
    assert seen == [None]
    back = S.load_raw(out)
    assert set(back) == set(raw) == set(got)
    for k in raw:
        assert back[k].dtype == raw[k].dtype and np.array_equal(back[k], raw[k], equal_nan=raw[k].dtype.kind == 'f'), k


# ---- the optional hook for the guarded re-runs of exp_astro_cyclefix.py ------------------------------------------------------------------

def _cyclefix_file(inp, path, runs=tuple(sorted(CAPPED)), capped_flag=0):
    """A planted cyclefix raw in the key format of exp_astro_cyclefix.py: 'G|' + the stored key of the run + the relative key (k<kappa>|tot, k<kappa>|blk22, bic, secs, det, wall, trace|...)."""
    rng = np.random.default_rng(5)
    arrays = {'meta|files': np.array(['a']), 'C|R|m2r1|x316_y124|iso|ours|k2|tot': np.arange(7), 'G|D|x110_y630|hp5|rex0.8|k2|tot': np.arange(7), 'G|T2|m53|x760_y590|default|k5|tot': np.arange(7),
              'G|T2|m2r1|x316_y124|iso|default|k5|tot': np.arange(7), 'G|T2|m2r1|x316_y124|iso|default|k5|blk22': np.zeros((25, 4))}                 # other families of the cyclefix raw are not read
    for name, pxy in runs:
        tp = A.truth_for_patch(inp['truth'][name], *pxy)
        pre = 'G|R|%s|x%d_y%d|iso|default|' % (name, pxy[0], pxy[1])
        for k, comp in (('k5', 0.8), ('k0', 0.82), ('k2', 0.85)):
            x, y, f = _catalogue(rng, np.array([c for c in tp['x']]), np.array([c for c in tp['y']]), tp['mag'], comp, 0.1, 10, OFFS[name])
            for kk, v in X._pack(k, A.score(np.c_[x, y], f, tp)).items():
                arrays[pre + kk] = v
        arrays[pre + 'bic'] = np.array([900, 1500, capped_flag])
        arrays[pre + 'secs'] = np.array([20.0, 1500.0])
        arrays[pre + 'wall'] = np.array([20.0])
        arrays[pre + 'trace|counts'] = np.arange(8)
        arrays[pre + 'det'] = np.zeros((3, 3))
    np.savez_compressed(path, **arrays)
    return arrays


def test_the_hook_reads_only_the_guarded_repeat_runs_with_their_counts_and_blocks(world, tmp_path):
    inp, raw = world
    assert S.guarded_rows(None) == {} and S.guarded_rows(str(tmp_path / 'absent.npz')) == {}
    arrays = _cyclefix_file(inp, str(tmp_path / 'cyclefix_raw.npz'))
    g = S.guarded_rows(str(tmp_path / 'cyclefix_raw.npz'))
    assert sorted(g) == [(n, p, 'iso', 'default') for n, p in sorted(CAPPED)]
    run = g[('m2r2', (133, 281), 'iso', 'default')]
    assert sorted(k for k in run) == ['bic', 'k0', 'k2', 'k5', 'secs']
    pre = 'G|R|m2r2|x133_y281|iso|default|'
    assert np.array_equal(run['k5']['tot'], arrays[pre + 'k5|tot']) and np.array_equal(run['k0']['blk'], arrays[pre + 'k0|blk22']) and np.array_equal(run['bic'], [900, 1500, 0])
    np.savez_compressed(str(tmp_path / 'none.npz'), **{'G|D|x110_y630|hp5|rex0.8|k2|tot': np.arange(7), 'meta|files': np.array(['a'])})
    assert S.guarded_rows(str(tmp_path / 'none.npz')) == {}


def test_the_hook_copies_the_rows_and_the_report_gains_section_8_1(world, tmp_path, small):
    inp, raw0 = world
    arrays = _cyclefix_file(inp, str(tmp_path / 'cyclefix_raw.npz'))
    inp2 = dict(inp, guarded=S.guarded_rows(str(tmp_path / 'cyclefix_raw.npz')))
    raw = S.build_raw(inp2)
    assert str(raw['meta|cyclefix'][0]) == 'applied' and str(raw0['meta|cyclefix'][0]) == 'not applied'
    assert not [k for k in raw0 if k.startswith('gr|')] and len([k for k in raw if k.startswith('gr|')]) == 3 * (3 * 2 + 2)
    assert np.array_equal(raw['gr|m2r3|x110_y501|iso|default|k5|blk22'], arrays['G|R|m2r3|x110_y501|iso|default|k5|blk22'])
    assert S.gkey(raw, 'm2r3', (110, 501), 'adq|default|k5') == 'gr|m2r3|x110_y501|iso|default|k5|blk22'                  # the chosen kernel is the isotropic run there
    assert S.gkey(raw, 'm2r2', (133, 281), 'adq|default|k5') is None and S.gkey(raw, 'm2r2', (133, 281), 'iso|default|k5') is not None
    assert S.gkey(raw, 'm2r1', (316, 124), 'iso|default|k5') is None and S.gkey(raw, 'm2r2', (133, 281), 'crowdsource (default)') is None and S.gkey(raw, 'm2r2', (133, 281), 'iso|ours|k2') is None
    assert S.still_defective(raw, S.PATCHES, ['iso|default|k5']) == [] and S.still_defective(raw0, S.PATCHES, ['iso|default|k5']) == S.defective(raw0, S.PATCHES, ['iso|default|k5'])
    b, sizes = S.guarded_blocks(raw, S.GEN_PATCHES, 'iso|default|k5')
    assert np.array_equal(b[:25], arrays['G|R|m2r2|x133_y281|iso|default|k5|blk22']) and np.array_equal(b[25:], raw['st|m2r2|x195_y689|iso|default|k5|blk22'])
    md = S.report(raw)
    assert '### 8.1 With the guarded re-runs' in md and '### 8.1' not in S.report(raw0) and 'it was used for this raw: section 8.1.' in md and 'it was not used for this raw (no section 8.1).' in S.report(raw0)
    sec = md.split('### 8.1')[1].split('## 9.')[0]
    # P1 default at alpha 0: stored rows are defective, the guarded rows are not
    ba, sz = S.guarded_blocks(raw, S.GEN_PATCHES, 'adq|default|k5')
    bb, _ = S.guarded_blocks(raw, S.GEN_PATCHES, 'iso|default|k5')
    d, ci = S.boot_dF(ba, bb, sz)
    rows = [r for r in _rows(sec, '| a priori default (kappa 5)') if len(r) == 7]
    assert [r[1] for r in rows] == ['stored rows', 'guarded rows in place of the capped ones'] and [r[2] for r in rows] == ['yes', 'no'] and rows[1][3] == X.fmt_iv(d, ci)
    # the default pairs: pooled over the six patches, stored against guarded; the m2r3 chosen-kernel runs are replaced through the isotropic guarded run
    pair = [r for r in _rows(sec, '| chosen kernel, a priori default vs photutils (default)') if len(r) == 7][0]
    ba, sz = S.guarded_blocks(raw, S.PATCHES, 'adq|default|k5')
    bb, _ = S.blocks(raw, S.PATCHES, 'photutils (default)')
    d, ci = S.boot_dF(ba, bb, sz)
    assert pair[4] == X.fmt_iv(d, ci) and pair[3] == 'yes' and pair[6] == 'no' and pair[1] == X.fmt_iv(*S.dF_stats(raw, S.PATCHES, 'adq|default|k5', 'photutils (default)'))
    # the runs table: entrants and cap flags stored / guarded
    run_row = [r for r in _rows(sec, '| m2r3 (110, 501) | iso | default |')][0]
    assert run_row[3] == '2500 / 1500' and run_row[4] == '1 / 0'
    # a guarded run that is itself capped does not make the row clean
    _cyclefix_file(inp, str(tmp_path / 'cap.npz'), capped_flag=1)
    rawc = S.build_raw(dict(inp, guarded=S.guarded_rows(str(tmp_path / 'cap.npz'))))
    assert S.still_defective(rawc, S.PATCHES, ['iso|default|k5']) == S.defective(rawc, S.PATCHES, ['iso|default|k5'])
    # a partial raw (one run) replaces only that run
    _cyclefix_file(inp, str(tmp_path / 'one.npz'), runs=(('m2r2', (133, 281)),))
    raw1 = S.build_raw(dict(inp, guarded=S.guarded_rows(str(tmp_path / 'one.npz'))))
    assert S.still_defective(raw1, S.PATCHES, ['iso|default|k5']) == [('m2r3', (370, 441)), ('m2r3', (110, 501))] and S.still_defective(raw1, S.GEN_PATCHES, ['iso|default|k5']) == []


# ---- the stored results and the data (skipped without them) ------------------------------------------------------------------------------

HAVE_STORED = os.path.exists(S.REPEAT_RAW) and os.path.exists(S.BASE) and all(os.path.exists(os.path.join(R.REPEAT[n]['root'], R.REPEAT[n]['truth'])) for n in S.FRAMES)
HAVE_DATA = HAVE_STORED and all(os.path.exists(R.patch_path(p[0], p[1], n)) for n, p in S.PATCHES) and os.path.exists(S.VALID)
HAVE_RAW = os.path.exists(S.RAW) and os.path.exists(S.MD)
needs_stored = pytest.mark.skipif(not HAVE_STORED, reason='the stored repeat results, the baselines\' catalogues or the truth lists are not on disk')
needs_data = pytest.mark.skipif(not HAVE_DATA, reason='the repeat patches or repeat_validation.json are not on disk')
needs_raw = pytest.mark.skipif(not (HAVE_RAW and HAVE_STORED), reason='the committed repeat_sens raw and report, or the stored results, are not on disk')


@needs_stored
def test_the_scorer_reproduces_every_stored_row_of_the_repeat_test_and_astro_common_on_a_sample():
    rep = X.load_raw(S.REPEAT_RAW)
    base = np.load(S.BASE, allow_pickle=False)
    n = 0
    for name, pxy in S.PATCHES:
        tp = A.truth_for_patch(R.truth_arrays(name), *pxy)
        pre = E.pre_of(name, pxy)
        cats = S.catalogues(rep, base, name, pxy)
        assert sorted(cats) == sorted(S.SC_KEYS + S.BASE_EXTRA)
        for key, det in cats.items():
            r = S.rescore(det, tp)
            assert np.array_equal(r['tot'][ICUT], rep[pre + key + '|tot']) and np.array_equal(r['blk'][20][ICUT], rep[pre + key + '|blk22']), (name, pxy, key)
            assert np.array_equal(r['off'], rep[pre + key + '|offset'], equal_nan=True)
            n += 1
    assert n == 60
    # and against astro_common.score itself on the two patches of the rule's frame (the stored rows came from it; this checks the function, not the stored numbers)
    for name, pxy in S.GEN_PATCHES:
        tp = A.truth_for_patch(R.truth_arrays(name), *pxy)
        det = S.catalogues(rep, base, name, pxy)['adq|ours|k2']
        ref = A.score(np.c_[det[0], det[1]], det[2], tp)
        r = S.rescore(det, tp)
        assert r['tot'][ICUT].tolist() == [ref['true22'], ref['matched22'], ref['det22'], ref['det22_matched'], ref['true_all'], ref['matched_all'], ref['det_all']]


@needs_stored
def test_the_baselines_are_read_through_the_frozen_post_filter_and_tags():
    rep = X.load_raw(S.REPEAT_RAW)
    base = np.load(S.BASE, allow_pickle=False)
    name, pxy = 'm2r1', (316, 124)
    cats = S.catalogues(rep, base, name, pxy)
    tuned = str(rep[E.pre_of(name, pxy) + 'crowdsource (tuned)|name'][0])
    assert tuned == 'crowdsource|snr5|ff0|qf0|spinf' and str(rep[E.pre_of(name, pxy) + 'photutils (tuned)|name'][0]) == 'photutils_ext_all_local_i6'
    k = 'x316_y124|crowdsource|'
    keep = (base[k + 'flux'] > 0) & (base[k + 'flux'] >= 5.0 * base[k + 'dflux'])
    assert np.array_equal(cats['crowdsource (tuned)'][0], base[k + 'x'][keep]) and len(cats['crowdsource (tuned)'][0]) <= len(cats['crowdsource (default)'][0])
    assert np.array_equal(cats['photutils (tuned)'][0], base['x316_y124|photutils_ext_all_local_i6|x']) and np.array_equal(cats['photutils (default)'][2], base['x316_y124|photutils|flux'])


@needs_data
def test_the_scan_offsets_are_the_masked_astrans_numbers_of_the_validation_json():
    scan, scan_m2, plain = S.read_scan()
    import json
    d = json.load(open(S.VALID))
    for name, pxy in S.PATCHES:
        e = d[name]['%d,%d' % pxy]
        assert scan[(name, pxy)] == (e['masked']['asTrans']['dx'], e['masked']['asTrans']['dy']) and plain[(name, pxy)] == (e['image_shift_px']['asTrans']['dx'], e['image_shift_px']['asTrans']['dy'])
    for w in X.TESTS:
        assert scan_m2[('M2', tuple(w))] == (d['m2']['%d,%d' % w]['masked']['header']['dx'], d['m2']['%d,%d' % w]['masked']['header']['dy'])
    # the numbers of the table of repeat_spec.md section 6.2
    for key, (dx, dy) in {('m2r2', (133, 281)): (-0.090, 0.257), ('m2r2', (195, 689)): (-0.008, 0.240), ('m2r1', (316, 124)): (0.074, 0.125), ('m2r1', (316, 224)): (0.109, 0.139),
                          ('m2r3', (370, 441)): (-0.028, 0.043), ('m2r3', (110, 501)): (-0.014, 0.043)}.items():
        assert np.allclose(scan[key], (dx, dy), atol=6e-4)
    assert np.allclose(scan_m2[('M2', (310, 630))], (0.093, 0.039), atol=6e-4) and np.allclose(plain[('m2r1', (316, 224))], (-0.517, 0.130), atol=6e-4)       # the plain scan fails there


@needs_data
@needs_raw
def test_the_star_fits_on_a_real_patch_are_those_of_the_committed_raw():
    raw = S.load_raw()
    label, plist, loader, tpath = [d for d in S.datasets() if d[0] == 'm2r2'][0]
    z = np.load(tpath)
    truth = {k: z[k] for k in ('px', 'py', 'm606', 'nf606')}
    pxy = (195, 689)
    patch = loader(*pxy)
    rows = S.box_rows(truth, *pxy)
    for sel in ('A', 'B'):
        stars, nsel = S.measure_patch(patch['dn'], patch['psf'], patch['gain'], patch['dark'], rows, sel)
        assert np.array_equal(stars, raw['own|%s|m2r2|x195_y689' % sel]) and nsel == raw['own|%s|m2r2|x195_y689|nsel' % sel][0]
    assert len(raw['own|B|m2r2|x195_y689']) > len(raw['own|A|m2r2|x195_y689']) >= 1
    # the offsets are of the order of the spec's 0.25 px in y for this frame
    st = np.concatenate([raw['own|B|m2r2|x133_y281'], raw['own|B|m2r2|x195_y689']])
    assert 0.1 < np.median(st[:, 4]) < 0.4 and abs(np.median(st[:, 3])) < 0.15


@needs_raw
def test_the_committed_raw_holds_the_stored_rows_unchanged_and_reproduces_the_p1_lines_of_repeat_md():
    raw = S.load_raw()
    rep = X.load_raw(S.REPEAT_RAW)
    for name, pxy in S.PATCHES:
        pre = E.pre_of(name, pxy)
        for key in S.STORED_ROWS:
            assert np.array_equal(raw['st|%s|%s|tot' % (S.pk(name, pxy), key)], rep[pre + key + '|tot']) and np.array_equal(raw['st|%s|%s|blk22' % (S.pk(name, pxy), key)], rep[pre + key + '|blk22'])
        for key in S.SC_KEYS + S.BASE_EXTRA:
            assert np.array_equal(raw['sc|base|%s|%s|tot' % (S.pk(name, pxy), key)][ICUT], rep[pre + key + '|tot'])
        assert np.array_equal(raw['st|%s|iso|default|bic' % S.pk(name, pxy)], rep[pre + 'iso|default|bic'])
    txt = E.p1_text(rep, E.done_patches(rep))
    for (label, ka, kb), line in zip(S.READOUTS, txt.split('\n')[1:]):
        d, ci = S.dF_stats(raw, S.GEN_PATCHES, ka, kb)
        assert X.fmt_iv(d, ci) in line, (label, line)
    # the sha256 of the inputs recorded in the raw are those of the files
    for k in raw:
        if k.startswith('meta|sha256|') and k.endswith('repeat_raw.npz'):
            assert str(raw[k][0]) == S._sha256_file(S.REPEAT_RAW)
    assert 'masked' in str(raw['meta|scan_key'][0]) and str(raw['meta|cyclefix'][0]) in ('applied', 'not applied')


@needs_raw
def test_the_committed_report_is_regenerated_byte_for_byte_from_the_committed_raw():
    raw = S.load_raw()
    assert S.SEEDS == tuple(range(1, 201)) and S.BIG_B == 40000
    md = S.report(raw)
    with open(S.MD) as f:
        assert f.read() == md


# ---- the sentences of the report follow the data -----------------------------------------------------------------------------------------

def test_the_closing_sentences_follow_the_data(world, md_small):
    inp, raw = world
    md = md_small
    # section 9: the planted rule acts in the frames of camera column 6 only
    assert 'it acts in every frame of camera column 6 and in no frame of the other camera column, so the stamp shape cannot be separated from the camera column' in md
    # the same planted data with the rule also acting in a camera-column-2 frame: no coincidence is claimed
    raw2 = dict(raw)
    for p in R.REPEAT['m2r1']['patches']:
        raw2['st|%s|kind' % S.pk('m2r1', tuple(p))] = np.array(['general'])
    txt = '\n'.join(S.sec_confound(raw2))
    assert ('does not coincide with a single camera column in these frames (the camera columns of the frames where it acts: 2, 6; of the others: 2)' in txt and 'coincides with the camera column' not in txt
            and 'where the rule picks a general kernel on at least one patch: m2r1, m2r2, M53, M92 (camera columns 2, 6, 6, 6)' in txt)
    # section 4: the cells with a loss are listed from the data (the planted isotropic kernel loses to the baselines on one patch at the tight tolerance)
    sec4 = '\n'.join(S.sec_p2(raw))
    listed = sec4.split('a baseline ahead) are: ')[1].split('. **What it does not show.**')[0].split('; ')
    per_c, _, v_c = S.pair_cells(raw, 'tol35', 'iso|ours|k2', 'crowdsource (tuned)')
    per_p, _, v_p = S.pair_cells(raw, 'tol35', 'iso|ours|k2', 'photutils (tuned)')
    assert v_c[2] == 'loss' and v_p[2] == 'loss' and S.PATCHES[2] == ('m2r2', (133, 281))
    assert '0.35 px / 0.5 mag, isotropic kernel vs crowdsource (tuned), m2r2 (133, 281)' in listed and '0.35 px / 0.5 mag, isotropic kernel vs photutils (tuned), m2r2 (133, 281)' in listed
    exc_rows = [l for l in sec4.split('### 4.5')[1].split('**What section 4 shows.**')[0].split('\n') if l.endswith('| loss |')]
    assert len(listed) == len(exc_rows) == 2
    raw3 = dict(raw)
    for n, p in S.PATCHES[:1]:
        for scn in S.SCENARIOS:
            for key in ('adq|ours|k2', 'iso|ours|k2'):
                raw3['sc|%s|%s|%s|blk' % (scn, S.pk(n, p), key)] = np.zeros_like(raw['sc|%s|%s|%s|blk' % (scn, S.pk(n, p), key)])
    sec4b = '\n'.join(S.sec_p2(raw3))
    assert 'm2r1 (316, 124)' in sec4b.split('a baseline ahead) are: ')[1].split('**What it does not show.**')[0] and len(sec4b.split('a baseline ahead) are: ')[1].split('. **What it does not show.**')[0].split('; ')) > 10
    # section 5 and 7: the capped patches are named from the flags
    assert 'hit the entry cap on m2r3 (370, 441), m2r3 (110, 501) (where the rule keeps' in md and 'the scaled default runs hit it on no patch' in md
    assert 'A comparison of valid runs for the isotropic default on 3 patches' in md and '* **Defective runs.** 3 a priori default runs' in md
    # section 3: the number of patches that keep the isotropic kernel
    assert '4 of the 6 patches keep the isotropic kernel' in md
    # section 10: the isotropic run is further from the truth in the planted frame
    assert "the isotropic run's matches sit further from the unshifted truth" in md
    raw4 = dict(raw)
    raw4['st|%s|iso|ours|k2|offset' % S.pk('m2r2', (133, 281))] = np.array([0.0, 0.0])
    assert 'not further on every patch' in '\n'.join(S.sec_positions(raw4))


def test_the_box_of_a_patch_is_half_a_pixel_wider_than_the_pixel_grid():
    x0, y0 = 133, 281
    px = np.array([x0 - 0.4, x0 - 0.6, x0 + 99.4, x0 + 99.6, x0 + 50.0, x0 + 50.0, x0 + 50.0, x0 + 50.0])
    py = np.array([y0 + 50.0, y0 + 50.0, y0 + 50.0, y0 + 50.0, y0 - 0.4, y0 - 0.6, y0 + 99.4, y0 + 99.6])
    truth = dict(px=px, py=py, m606=np.arange(8.0) + 17.0, nf606=np.arange(8.0))
    rows = S.box_rows(truth, x0, y0)
    assert rows['mag'].tolist() == [17.0, 19.0, 21.0, 23.0] and rows['nf'].tolist() == [0.0, 2.0, 4.0, 6.0]
    assert np.allclose(rows['x'], [-0.4, 99.4, 50.0, 50.0]) and np.allclose(rows['y'], [50.0, 50.0, -0.4, 99.4])


def test_main_passes_the_cyclefix_option_to_run_and_asks_for_one_action(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(S, 'run', lambda out_path=S.RAW, cyclefix=None: seen.append(cyclefix))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_repeat_sens.py', '--run', '--cyclefix', str(tmp_path / 'c.npz')])
    S.main()
    monkeypatch.setattr(sys, 'argv', ['exp_astro_repeat_sens.py', '--run'])
    S.main()
    assert seen == [str(tmp_path / 'c.npz'), None]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_repeat_sens.py', '--cyclefix', 'x'])
    with pytest.raises(SystemExit):
        S.main()
    assert len(seen) == 2
