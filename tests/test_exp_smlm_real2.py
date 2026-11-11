"""Tests of the pure functions of experiments/exp_smlm_real2.py (the second real-data pore study: registration of raw frames against a published table, the four-dial ISS
post-filters and their coordinate ascent, the pore model of a test set, the merged-frame overlap test, the paired bootstrap on synthetic data, and the reading of the outcomes
under the plan's rules). No downloads, no Picasso and no tifffile needed. About two minutes."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_smlm_real_gt as G      # noqa: E402
import exp_smlm_real2 as R        # noqa: E402

R0 = 87.0


def ring_points(rng, centre, n, R_=R0, sd=10.0):
    th = rng.uniform(0, 2 * np.pi, n)
    r = R_ + sd * rng.standard_normal(n)
    return np.c_[centre[0] + r * np.cos(th), centre[1] + r * np.sin(th)]


# ----------------------------------------------------------------------------------------------------------------------------------
# merged frames
# ----------------------------------------------------------------------------------------------------------------------------------
def test_make_groups_are_frames_seven_apart_inside_blocks_of_twenty_one_and_the_leftover_is_dropped():
    g = R.make_groups(50)
    assert g.shape == (14, 3)                                     # two blocks of 21, 8 frames left over
    assert g[0].tolist() == [0, 7, 14] and g[6].tolist() == [6, 13, 20] and g[7].tolist() == [21, 28, 35]
    assert g.max() == 41 and len(np.unique(g)) == 42
    assert R.make_groups(20).shape == (0, 3)


def test_merge_frames_subtracts_the_background_k_minus_one_times():
    rng = np.random.default_rng(0)
    fr = 5.0 + rng.random((21, 4, 4))
    bg = np.median(fr, axis=0)
    g = R.make_groups(21)
    S = R.merge_frames(fr, g)
    assert S.shape == (7, 4, 4)
    assert np.allclose(S[2], fr[2] + fr[9] + fr[16] - 2 * bg, atol=1e-5)
    assert np.allclose(R.merge_frames(fr, g, bg=np.zeros((4, 4)))[0], fr[0] + fr[7] + fr[14], atol=1e-5)


# ----------------------------------------------------------------------------------------------------------------------------------
# registration
# ----------------------------------------------------------------------------------------------------------------------------------
def planted_pair(shift_px=(3, -2), ny=60, nx=64, sigma=1.5, seed=0):
    rng = np.random.default_rng(seed)
    pix = 100.0
    pts = rng.uniform([6, 6], [nx - 6, ny - 6], (80, 2))                      # table positions in px (x, y)
    tab_xy = pts * pix
    B = R.render_table(tab_xy[:, 0], tab_xy[:, 1], (ny, nx), pix, sigma)
    A = np.roll(np.roll(B, shift_px[1], axis=0), shift_px[0], axis=1)         # the raw image is the table image displaced by +shift (raw = table + shift)
    return A, B, pix


def test_xcorr_shift_recovers_a_planted_displacement_with_its_sign_and_a_subpixel_part():
    A, B, pix = planted_pair()
    dx, dy, q = R.xcorr_shift(A, B, pix)
    assert abs(dx - 300.0) < 15.0 and abs(dy + 200.0) < 15.0 and q > 2.5
    # half a pixel
    ny, nx = B.shape
    from scipy import ndimage as ndi
    A2 = ndi.shift(B, (0.0, 0.5), order=3, mode='wrap')
    dx2, dy2, _ = R.xcorr_shift(A2, B, pix)
    assert abs(dx2 - 50.0) < 15.0 and abs(dy2) < 15.0


def test_render_table_places_a_localization_in_the_pixel_that_contains_it():
    H = R.render_table(np.array([250.0]), np.array([150.0]), (4, 6), 100.0, 0.5)
    iy, ix = np.unravel_index(np.argmax(H), H.shape)
    assert (iy, ix) == (1, 2)


def test_sum_image_removes_the_temporal_median_and_keeps_positive_parts_only():
    fr = np.zeros((5, 6, 6), np.float32)
    fr += 10.0
    fr[1, 2, 3] += 4.0
    fr[3, 2, 3] += 6.0
    fr[0, 0, 0] -= 3.0
    S = R.sum_image(fr)
    assert abs(S[2, 3] - 10.0) < 1e-5 and S[0, 0] == 0.0 and S[5, 5] == 0.0


def test_reg_blocks_partition_the_frames_into_near_equal_blocks_of_about_125():
    for n in (754, 1055, 2000, 100):
        b = R.reg_blocks(n)
        assert b[0][0] == 0 and b[-1][1] == n and all(b[i][1] == b[i + 1][0] for i in range(len(b) - 1))
        sizes = [y - x for x, y in b]
        assert max(sizes) - min(sizes) <= 1
    assert len(R.reg_blocks(2000)) == 16 and len(R.reg_blocks(754)) == 6 and len(R.reg_blocks(1055)) == 8 and len(R.reg_blocks(40)) == 1


def test_fit_shift_curve_recovers_a_cubic_and_ignores_one_outlying_block_but_reports_its_residual():
    n = 2000
    blocks = R.reg_blocks(n)
    t = np.array([(a + b) / (2.0 * n) for a, b in blocks])
    truth = np.c_[300 + 40 * t - 30 * t ** 2 + 20 * t ** 3, 100 - 15 * t]
    sh = truth.copy()
    sh[7, 1] += 40.0
    coef, res = R.fit_shift_curve(blocks, sh, n)
    assert np.abs(res[np.arange(len(blocks)) != 7]).max() < 1.0
    assert abs(res[7, 1] - 40.0) < 2.0
    # evaluated per frame
    raw = {'reg|X|coef': coef}
    R.TESTS['X'] = dict(ds='D1', field=None, frames=(1, n), ref_ranges=[], group='D1')
    try:
        s = R.shift_of(raw, 'X')
        assert s.shape == (n, 2) and abs(s[0, 0] - (300 + 40 * 0.00025)) < 2.0 and abs(s[-1, 1] - (100 - 15)) < 2.0
    finally:
        del R.TESTS['X']


def test_gate_r1_uses_the_largest_deviation_and_the_smallest_quality():
    raw = {'reg|T1|resid': np.array([[1.0, 2.0], [19.9, 0.0]]), 'reg|T1|blocks': np.array([[0, 1, 0, 0, 3.0], [1, 2, 0, 0, 2.5]])}
    assert R.gate_r1(raw, 'T1')['passed']
    raw['reg|T1|resid'][1, 0] = -20.5
    assert not R.gate_r1(raw, 'T1')['passed']
    raw['reg|T1|resid'][1, 0] = 5.0
    raw['reg|T1|blocks'][1, 4] = 2.4
    assert not R.gate_r1(raw, 'T1')['passed']


# ----------------------------------------------------------------------------------------------------------------------------------
# the four ISS post-filters and the ascent
# ----------------------------------------------------------------------------------------------------------------------------------
def test_filter_dets_amplitude_floor_and_width_limits_are_inclusive():
    d = np.array([[0, 0, 100.0, 1.0], [0, 0, 200.0, 1.2], [0, 0, 300.0, 1.4], [0, 0, 50.0, 1.2]])
    f = dict(fw=0.5, smin=1.0, smax=1.4)
    out = R.filter_dets([d], f, w_ref=200.0, s_scale=1.0)[0]
    assert out[:, 2].tolist() == [100.0, 200.0, 300.0]                       # fw x w_ref = 100 kept, 50 dropped; sigma 1.0 and 1.4 kept
    out = R.filter_dets([d], dict(fw=0.0, smin=None, smax=1.2), 200.0, 1.0)[0]
    assert out[:, 3].max() <= 1.2 and len(out) == 3
    out = R.filter_dets([d], dict(fw=0.0, smin=1.1, smax=3.2), 200.0, 1.0)[0]
    assert sorted(out[:, 2].tolist()) == [50.0, 200.0, 300.0]
    assert len(R.filter_dets([np.zeros((0, 4))], f, 200.0, 1.0)[0]) == 0
    out = R.filter_dets([d], dict(fw=0.0, smin=1.0, smax=1.0), 200.0, 2.0)[0]            # widths scale with s: [2, 2]
    assert len(out) == 0


def test_w_ref_is_the_median_amplitude_of_the_kappa_10_detections():
    dets = {10.0: [np.array([[0, 0, 1.0, 1], [0, 0, 9.0, 1]]), np.zeros((0, 4)), np.array([[0, 0, 2.0, 1]])], 5.0: [np.array([[0, 0, 100.0, 1]])]}
    assert R.w_ref_of(dets) == 2.0


def test_ascent_moves_only_on_a_strict_improvement_in_the_order_of_the_dials_and_stops_when_nothing_moves():
    seen = []
    target = dict(kappa=10.0, fw=1.0, smin=1.1, smax=1.6)

    def ev(p):
        seen.append(tuple(p[k] for k in R.ISS_KEYS))
        return -sum(abs((0 if p[k] is None else p[k]) - target[k]) for k in R.ISS_KEYS)
    cur, best, rounds, n_ev = R.ascent_iss(ev, R.ISS_GRIDS)
    assert cur == target and best == 0.0
    assert rounds == 2 and n_ev == len(seen)                                # round 1 reaches the optimum, round 2 confirms
    assert seen[0] == (5.0, 0.0, None, 3.2)

    # ties never move
    cur2, best2, rounds2, n2 = R.ascent_iss(lambda p: 1.0, R.ISS_GRIDS)
    assert cur2 == R.ISS_START and rounds2 == 1
    # at most max_rounds rounds
    cur3, _, rounds3, _ = R.ascent_iss(lambda p: -abs(p['kappa'] - 50.0) if True else 0, {'kappa': (6.0, 8.0), 'fw': (0.0,), 'smin': (None,), 'smax': (3.2,)}, max_rounds=1)
    assert rounds3 == 1 and cur3['kappa'] == 8.0


def test_choose_cfg_takes_the_largest_value_and_the_earlier_configuration_on_a_tie():
    r = [('a', 40.0, {}), ('b', 41.0, {}), ('c', 41.0000001, {}), ('d', 39.0, {})]
    assert R.choose_cfg(r)[0] == 'b'
    assert R.choose_cfg([('x', 1.0), ('y', 2.0)])[0] == 'y'


def test_paths_round_trip_through_the_flat_arrays():
    rng = np.random.default_rng(0)
    dets = {k: [rng.random((rng.integers(0, 4), 4)) for _ in range(5)] for k in R.KAPPAS2}
    raw = {}
    R.put_paths(raw, 'p', dets, 5, np.arange(5.0), np.arange(5))
    back = R.get_paths(raw, 'p')
    for k in R.KAPPAS2:
        for a, b in zip(dets[k], back[k]):
            assert a.shape == b.shape and np.allclose(a, b)
    assert int(raw['p|n_frames']) == 5


def test_iss_frame2_reads_every_requested_kappa_and_finds_a_planted_emitter():
    from src.core import imaging as IM
    import exp_smlm_real as REAL
    rng = np.random.default_rng(6)
    c = np.array([[20.3, 18.6]])
    P = IM.PixelPSF(40, 40, 1.5, parts=((1.0, 1.0),))
    mean = 300.0 + P.synth(c, np.array([3000.0]))
    frame = mean + np.sqrt(0.8 * (mean - 200.0) + 100.0) * rng.standard_normal((40, 40))
    cfg = G.grid_cfg('hp5', 1.5)
    cal = dict(bg=np.full((40, 40), 300.0), a=0.8, b0=100.0, x_ref=200.0)
    dets, secs, ent = R.iss_run2(np.stack([frame]), cal, cfg, 100.0, kappas=(5.0, 30.0), log=lambda s: None)
    assert set(dets) == {5.0, 30.0}
    d = dets[5.0][0]
    k = np.argmin(np.hypot(d[:, 0] - 100.0 * 20.8, d[:, 1] - 100.0 * 19.1))
    assert np.hypot(d[k, 0] - 2080.0, d[k, 1] - 1910.0) < 30.0 and abs(d[k, 2] / 3000.0 - 1) < 0.25 and d.shape[1] == 4


# ----------------------------------------------------------------------------------------------------------------------------------
# the pore model of a test set
# ----------------------------------------------------------------------------------------------------------------------------------
def test_scorer2_scores_pores_inside_the_region_and_beyond_x_min_and_its_areas_add_up():
    sp = G.DS['D1']
    xs = 6500.0 + 1000.0 * np.arange(5)
    cen = np.array([(x, 1500.0 + 1000.0 * j) for x in xs for j in range(5)] + [(1500.0, 3000.0), (6500.0, 150.0)])      # a development pore, a pore at the border
    sc = R.Scorer2(cen, R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0), x_min=6200.0)
    assert sc.n_pores == 25 and sc.scored.sum() == 25
    assert abs(sc.area.sum() / (25 * np.pi * sc.rmax ** 2) - 1.0) < 0.02
    assert np.allclose(sc.area_pore.sum(0), sc.area)
    r, j, keep = sc.assign(np.array([[6500.0 + R0, 1500.0], [1500.0 + R0, 3000.0]]))
    assert list(keep) == [True, False]
    # a convention shift moves the pores and drops those whose disc leaves the region
    cen2 = np.r_[cen[:25], [[6500.0, 250.0]]]
    assert R.Scorer2(cen2, R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0), x_min=6200.0).n_pores == 25 + 1 - 0
    assert R.Scorer2(cen2, R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, -53.0), x_min=6200.0).n_pores == 25
    # a pore exactly at x_min is not scored
    assert R.Scorer2(np.r_[cen, [[6200.0, 6000.0]]], R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0), x_min=6200.0).n_pores == 25
    # without x_min every pore of the region is scored
    assert R.Scorer2(cen, R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0)).n_pores == 26


def test_pick_window_takes_the_densest_window_and_never_one_overlapping_the_excluded_origin():
    sp = G.DS['D2']
    R_ = 53.0
    cen = np.array([(x, y) for x in np.linspace(30000, 36000, 12) for y in np.linspace(20000, 26000, 12)])        # a dense patch near (33000, 23000) nm
    bw, n = R.pick_window(cen, R_, sp)
    assert n >= 100 and abs(bw[0] * 106 + 64 * 106 - 33000) < 4000
    bw2, n2 = R.pick_window(cen, R_, sp, exclude=bw)
    assert abs(bw2[0] - bw[0]) >= 128 or abs(bw2[1] - bw[1]) >= 128
    assert n2 <= n


# ----------------------------------------------------------------------------------------------------------------------------------
# part B
# ----------------------------------------------------------------------------------------------------------------------------------
def test_reference_events_use_the_frames_of_the_group_the_shift_the_border_and_the_bright_threshold():
    tab = dict(frame=np.array([11, 11, 12, 15, 19, 19, 19]), x=np.array([1000.0, 3000.0, 1500.0, 1200.0, 1000.0, 1000.0, 2000.0]),
               y=np.array([1000.0, 1000.0, 1000.0, 1100.0, 150.0, 1000.0, 2000.0]), ph=np.array([500.0, 900.0, 100.0, 700.0, 800.0, 600.0, 650.0]))
    groups = np.array([[0, 4, 8]])                                                # set positions; table frame = a + position, a = 11
    shifts = np.zeros((9, 2))
    shifts[4] = [10.0, -20.0]
    ev = R.reference_events(tab, 11, groups, (0.0, 0.0), np.array([-53.0, -53.0]), shifts, (0.0, 0.0, 4000.0, 4000.0), 600.0)[0]
    # frame 11 (position 0): two events; position 4 (frame 15): one with the shift; position 8 (frame 19): y = 150 (97 after the convention) is within the border, the others are not
    assert len(ev) == 5 and 800.0 not in ev[:, 2].tolist()
    e = ev[ev[:, 3] == 1][0]
    assert abs(e[0] - (1200 - 53 + 10)) < 1e-9 and abs(e[1] - (1100 - 53 - 20)) < 1e-9
    assert ev[ev[:, 2] == 500.0][0, 4] == 0.0 and ev[ev[:, 2] == 900.0][0, 4] == 1.0 and ev[ev[:, 2] == 600.0][0, 4] == 1.0         # a count equal to the threshold is bright
    # the window origin is subtracted
    ev2 = R.reference_events(tab, 11, groups, (400.0, 0.0), np.array([0.0, 0.0]), np.zeros((9, 2)), (400.0, 0.0, 4400.0, 4000.0), 600.0)[0]
    assert abs(ev2[ev2[:, 2] == 500.0][0, 0] - 600.0) < 1e-9


def test_match_group_is_one_to_one_and_inclusive_at_the_radius():
    Rr = np.array([[0.0, 0.0, 1, 0, 1], [100.0, 0.0, 1, 1, 1]])
    D = np.array([[10.0, 0.0], [50.0, 0.0], [300.0, 0.0]])
    mr, md, which = R.match_group(D, Rr, 80.0)
    assert mr.tolist() == [True, True] and md.tolist() == [True, True, False] and which.tolist() == [0, 1, -1]
    # one detection between two references matches only one of them
    mr, md, which = R.match_group(np.array([[50.0, 0.0]]), Rr, 80.0)
    assert mr.sum() == 1 and md.sum() == 1
    # a pair exactly at the radius matches
    mr, md, which = R.match_group(np.array([[80.0, 0.0]]), Rr[:1], 80.0)
    assert mr.tolist() == [True]
    mr, md, which = R.match_group(np.zeros((0, 2)), Rr, 80.0)
    assert mr.sum() == 0 and len(which) == 0


def test_close_flags_pairs_bright_events_of_different_constituents_only():
    Rr = np.array([[0.0, 0.0, 1, 0, 1], [100.0, 0.0, 1, 1, 1], [150.0, 0.0, 1, 1, 1], [0.0, 120.0, 1, 0, 1], [60.0, 0.0, 1, 2, 0]])
    flag, pairs = R.close_flags(Rr, 130.0)
    # (0,1) 100 apart; (0,2) 150 no; (1,2) same constituent; (1,3) 156 no; (3,0) same constituent; the dim event (4) never counts
    assert sorted(pairs) == [(0, 1)] and flag.tolist() == [True, True, False, False, False]
    flag, pairs = R.close_flags(Rr, 100.0)
    assert sorted(pairs) == [(0, 1)]                                              # within the distance, inclusive
    flag, pairs = R.close_flags(Rr, 99.9)
    assert pairs == []


def test_group_counts_tp_fp_fn_dim_events_and_close_pair_recovery():
    # two bright references 100 nm apart (different constituents), one dim reference elsewhere
    Rr = np.array([[0.0, 0.0, 900, 0, 1], [100.0, 0.0, 900, 1, 1], [1000.0, 0.0, 100, 2, 0], [2000.0, 0.0, 900, 2, 1]])
    D = np.array([[5.0, 0.0], [1005.0, 0.0], [3000.0, 3000.0]])
    c = R.group_counts(D, Rr, 70.0, 130.0)
    # TP: the detection at 5 (bright ref 0); the detection at 1005 matches the dim event (neither); the one at 3000 is a FP
    # FN: bright refs 1 and 3 unmatched; close events: refs 0 and 1 (2), matched 1; one pair, resolved 0
    assert c.tolist() == [1, 1, 2, 2, 1, 1, 0, 3]
    D2 = np.array([[5.0, 0.0], [95.0, 0.0]])
    c2 = R.group_counts(D2, Rr, 70.0, 130.0)
    assert c2[0] == 2 and c2[1] == 0 and c2[3] == 2 and c2[4] == 2 and c2[5] == 1 and c2[6] == 1
    c3 = R.group_counts(np.zeros((0, 2)), Rr, 70.0, 130.0)
    assert c3.tolist() == [0, 0, 3, 2, 0, 1, 0, 0]
    c4 = R.group_counts(D, np.zeros((0, 5)), 70.0, 130.0)
    assert c4.tolist() == [0, 3, 0, 0, 0, 0, 0, 3]


def test_jb_and_b_stats():
    c = np.array([[3, 1, 2, 4, 2, 2, 1, 4], [1, 0, 0, 0, 0, 0, 0, 1]], float)
    assert abs(R.jb_of(c) - 4 / 7) < 1e-12
    s = R.b_stats(c)
    assert abs(s['recall'] - 4 / 6) < 1e-12 and abs(s['fp_per_frame'] - 0.5) < 1e-12 and abs(s['close_recall'] - 0.5) < 1e-12 and abs(s['pair_resolution'] - 0.5) < 1e-12
    assert abs(s['det_per_frame'] - 2.5) < 1e-12


def test_a_second_emitter_inside_the_resolution_is_a_close_pair_that_a_single_detection_cannot_resolve():
    rng = np.random.default_rng(3)
    ref = np.array([[500.0, 500.0, 900, 0, 1], [570.0, 500.0, 900, 1, 1]])
    c = R.group_counts(np.array([[535.0, 500.0]]), ref, 40.0, 130.0)
    assert c[0] == 1 and c[2] == 1 and c[4] == 1 and c[5] == 1 and c[6] == 0          # one detection between the two matches one of them only


# ----------------------------------------------------------------------------------------------------------------------------------
# part A statistics and the reading of the rules
# ----------------------------------------------------------------------------------------------------------------------------------
def test_yield_at_floor_interpolates_log_log_and_clamps():
    fl, N = [100.0, 400.0, 1600.0], [50.0, 100.0, 200.0]
    assert abs(R.yield_at_floor(fl, N, 200.0) - 70.71) < 0.1
    assert R.yield_at_floor(fl, N, 10.0) == pytest.approx(50.0) and R.yield_at_floor(fl, N, 1e6) == pytest.approx(200.0)
    assert R.yield_at_floor(fl[::-1], N[::-1], 400.0) == pytest.approx(100.0)
    assert np.isnan(R.yield_at_floor([1.0], [1.0], 1.0)) and np.isnan(R.yield_at_floor(fl, N, 0.0))


def mk_pair(dE_ci, dEp_ci, ratio):
    return {'pairs': {'ISS4-R2E': dict(dE_ci=dE_ci, dE_ci_pores=dEp_ci, off_ring_ratio=ratio)}}


def test_read_a_follows_the_rule():
    assert R.read_A(mk_pair([1.0, 5.0], [0.5, 6.0], 1.2)) == 'ISS better'
    assert R.read_A(mk_pair([1.0, 5.0], [0.5, 6.0], 1.0)) == 'tradeoff'                  # the ratio must exceed 1
    assert R.read_A(mk_pair([1.0, 5.0], [-0.5, 6.0], 1.5)).startswith('ISS better in the frame bootstrap only')
    assert R.read_A(mk_pair([-5.0, -0.1], [-6.0, 1.0], 0.5)) == 'Picasso better'
    assert R.read_A(mk_pair([-2.9, 2.9], [-3.0, 3.0], 1.0)) == 'equivalent'
    assert R.read_A(mk_pair([-3.0, 2.9], [-3.0, 3.0], 1.0)) == 'inconclusive'            # the lower limit must be above -3
    assert R.read_A(mk_pair([-2.9, 3.0], [-3.0, 3.0], 1.0)) == 'inconclusive'
    assert R.read_A(mk_pair([0.0, 5.0], [-1.0, 6.0], 1.0)) == 'inconclusive'


def test_overall_a_needs_both_valid_datasets_for_an_advantage():
    assert R.overall_A({'D1': 'ISS better', 'D2': 'ISS better'}, {'D1': True, 'D2': True}).startswith('ISS better on BOTH')
    t = R.overall_A({'D1': 'ISS better', 'D2': 'ISS better'}, {'D1': False, 'D2': True})
    assert t.startswith('NOT CONCLUDED') and 'D1 void' in t.replace('D1 and D1', 'D1') or 'D1' in t
    assert 'D2 void' not in t
    assert R.overall_A({'D1': 'Picasso better', 'D2': 'equivalent'}, {'D1': True, 'D2': True}).startswith('ISS WORSE')
    assert R.overall_A({'D1': 'equivalent', 'D2': 'equivalent'}, {'D1': True, 'D2': True}).startswith('EQUIVALENT')
    assert R.overall_A({'D1': 'ISS better', 'D2': 'Picasso better'}, {'D1': True, 'D2': True}).startswith('MIXED')
    assert R.overall_A({'D1': 'tradeoff', 'D2': 'ISS better'}, {'D1': True, 'D2': True}).startswith('MIXED')


def mk_b(dJ, dc, est=0.0):
    return {'pairs': {'ISSB-R2B': {'J': dict(est=0, ci=dJ), 'close_recall': dict(est=est, ci=dc)}}}


def test_read_b_and_overall_b_follow_the_rule():
    assert R.read_B(mk_b([0.01, 0.05], [0.0, 0.1])) == 'ISS better'
    assert R.read_B(mk_b([0.01, 0.05], [-0.01, 0.1], est=0.02)) == 'ISS better'             # the point estimate of the close-pair recall difference decides, not its interval
    assert R.read_B(mk_b([0.01, 0.05], [-0.1, 0.1], est=-0.001)).startswith('ISS better in J_B but')
    assert R.read_B(mk_b([-0.05, -0.01], [0.0, 0.1])) == 'Picasso better'
    assert R.read_B(mk_b([-0.029, 0.029], [0.0, 0.1])) == 'equivalent'
    assert R.read_B(mk_b([-0.03, 0.029], [0.0, 0.1])) == 'inconclusive'
    assert R.overall_B({'D1': 'ISS better', 'D2': 'ISS better'}, {'D1': True, 'D2': True}).startswith('ISS better on BOTH')
    assert R.overall_B({'D1': 'ISS better', 'D2': 'equivalent'}, {'D1': True, 'D2': True}).startswith('MIXED')
    assert R.overall_B({'D1': 'ISS better', 'D2': 'ISS better'}, {'D1': True, 'D2': False}).startswith('NOT CONCLUDED')


# ----------------------------------------------------------------------------------------------------------------------------------
# the analysis on a synthetic test set (bootstrap, pooling, the paired differences)
# ----------------------------------------------------------------------------------------------------------------------------------
def synth_set(monkeypatch, n_frames=60, seed=0, better=True):
    """One synthetic test set 'S1' (D1 geometry): pores on a grid; rows R0 / R1E / R2E / ISS4 / curve rows / PUB as Locs with a ring and a floor."""
    rng = np.random.default_rng(seed)
    sp = G.DS['D1']
    cen = np.array([(6500.0 + 1000.0 * i, 1500.0 + 1000.0 * j) for i in range(4) for j in range(4)])
    region = (0.0, 0.0, 12100.0, 12100.0)

    def loc(n_ring, n_floor, w=12.0):
        ring_pts, fr = [], []
        for _ in range(n_ring):
            c = cen[rng.integers(len(cen))]
            th = rng.uniform(0, 2 * np.pi)
            r = R0 + w * rng.standard_normal()
            ring_pts.append(c + r * np.array([np.cos(th), np.sin(th)]))
            fr.append(rng.integers(n_frames))
        for _ in range(n_floor):
            c = cen[rng.integers(len(cen))]
            ring_pts.append(c + rng.uniform(-230, 230, 2))
            fr.append(rng.integers(n_frames))
        fr = np.array(fr)
        o = np.argsort(fr, kind='stable')
        return G.Loc(fr[o], np.array(ring_pts)[o])
    rows = {'R0': loc(300, 300), 'R1E': loc(250, 200), 'R2E': loc(220, 120), 'ISS4': loc(320 if better else 200, 120 if better else 120), 'ISS_OLD': loc(300, 400), 'PUB': loc(250, 250)}
    for k in R.KAPPAS2:
        rows['K%g' % k] = loc(int(200 + 3 * k) if better else 150, int(60 + 10 * k))
    monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, n_frames), ref_ranges=[(1, 1)], group='D1'))
    monkeypatch.setitem(R.GROUPS, 'S1', ['S1'])
    sc = R.Scorer2(cen, R0, region, sp, (0.0, 0.0), x_min=6200.0)
    monkeypatch.setattr(R, 'scorer_of', lambda raw, name, old=None: sc)
    monkeypatch.setattr(R, 'set_rows', lambda raw, old, name, register=True, curve=True: {k: v for k, v in rows.items() if curve or not k.startswith('K')})
    raw = {'ref|gates': np.array(json.dumps({'S1': dict(n_exp=900.0)})), 'cal|S1|sigma0': np.array(json.dumps({'sigma0': 1.5}))}
    return raw, sc, rows


def test_analyse_a_on_a_synthetic_set_is_deterministic_orders_the_rows_and_gives_a_positive_interval_when_the_iss_is_better(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    a = R.analyse_A(raw, None, 'S1', n_boot=40)
    b = R.analyse_A(raw, None, 'S1', n_boot=40)
    assert json.dumps(a, sort_keys=True, default=float) == json.dumps(b, sort_keys=True, default=float)
    pr = a['pairs']['ISS4-R2E']
    assert pr['dE'] > 0 and pr['dE_ci'][0] <= pr['dE'] <= pr['dE_ci'][1] and pr['dE_ci_pores'][1] > 0
    assert a['rows']['R2E']['E'] > a['rows']['R0']['E'] - 40
    assert a['n_pores'] == sc.n_pores and abs(a['n_exp'] - 900.0) < 1e-9
    assert set(a['curve']) == {'K%g' % k for k in R.KAPPAS2}
    assert np.isfinite(pr['off_ring_ratio']) and pr['off_ring_ratio_ci'][0] <= pr['off_ring_ratio_ci'][1]
    # the unregistered / curve-free variant
    c = R.analyse_A(raw, None, 'S1', n_boot=10, with_curve=False)
    assert 'curve' not in c and np.isnan(c['pairs']['ISS4-R2E']['off_ring_ratio'])


def test_analyse_a_pooled_over_two_members_sums_histograms_and_expected_yields(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    monkeypatch.setitem(R.TESTS, 'S2', dict(ds='D1', field=None, frames=(1, 60), ref_ranges=[(1, 1)], group='D1'))
    monkeypatch.setitem(R.GROUPS, 'SS', ['S1', 'S2'])
    raw['ref|gates'] = np.array(json.dumps({'S1': dict(n_exp=900.0), 'S2': dict(n_exp=400.0)}))
    raw['cal|S2|sigma0'] = np.array(json.dumps({'sigma0': 1.5}))
    one = R.analyse_A(raw, None, 'S1', n_boot=5)
    two = R.analyse_A(raw, None, 'SS', n_boot=5)
    assert two['n_frames'] == 120 and two['n_pores'] == 2 * one['n_pores'] and abs(two['n_exp'] - 1300.0) < 1e-9
    assert two['rows']['R2E']['n_loc'] == 2 * one['rows']['R2E']['n_loc']
    # the same data twice: the pooled radial histogram is twice the single one, so the ring yield and the floor double and the width is unchanged
    assert abs(two['rows']['R2E']['fit']['N'] / one['rows']['R2E']['fit']['N'] - 2.0) < 0.02
    assert abs(two['rows']['R2E']['fit']['floor'] / one['rows']['R2E']['fit']['floor'] - 2.0) < 0.02
    assert abs(two['rows']['R2E']['fit']['w'] - one['rows']['R2E']['fit']['w']) < 0.2
    assert two['members_info']['S1']['n_exp'] == 900.0 and two['members_info']['S2']['n_exp'] == 400.0


def test_a_set_that_is_worse_gives_a_negative_difference(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=False)
    a = R.analyse_A(raw, None, 'S1', n_boot=20)
    assert a['pairs']['ISS4-R2E']['dE'] < 0 or a['pairs']['ISS4-R2E']['dJ'] < 0


def test_the_equal_off_ring_ratio_is_the_iss_curve_read_at_the_floor_of_r2e(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    a = R.analyse_A(raw, None, 'S1', n_boot=5)
    fl = [a['curve']['K%g' % k]['floor'] for k in R.KAPPAS2]
    Ns = [a['curve']['K%g' % k]['N'] for k in R.KAPPAS2]
    expect = R.yield_at_floor(fl, Ns, a['rows']['R2E']['fit']['floor']) / a['rows']['R2E']['fit']['N']
    assert abs(a['off_ring_ratio'] - expect) < 1e-9


def test_analyse_b_pools_counts_and_reads_the_bootstrap_with_the_same_blocks(monkeypatch):
    rng = np.random.default_rng(0)

    def counts(p):
        n = 30
        tp = rng.binomial(10, p, n)
        return np.c_[tp, rng.binomial(5, 0.3, n), 10 - tp, np.full(n, 4), rng.binomial(4, p, n), np.full(n, 2), rng.binomial(2, p, n), tp + 2]
    mem = lambda: dict(name='S1', counts={'R0': counts(0.5), 'R2B': counts(0.6), 'ISSB': counts(0.8)}, n_groups=30, rm=80.0, close_nm=200.0, ref_per_frame=12.0, bright_per_frame=6.0)
    monkeypatch.setattr(R, 'b_member', lambda raw, old, name: mem())
    monkeypatch.setitem(R.GROUPS, 'SB', ['S1'])
    a = R.analyse_B({}, None, 'SB', n_boot=30)
    pr = a['pairs']['ISSB-R2B']
    assert pr['J']['est'] > 0 and pr['J']['ci'][0] < pr['J']['est'] < pr['J']['ci'][1] and a['n_groups'] == 30
    assert abs(a['rows']['ISSB']['J'] - R.jb_of(np.concatenate([mem()['counts']['ISSB']])) ) < 1.0
    assert a['rows']['ISSB']['ci']['J'][0] <= a['rows']['ISSB']['J'] <= a['rows']['ISSB']['ci']['J'][1]


# ----------------------------------------------------------------------------------------------------------------------------------
# freezing and the raw
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_freeze_hash_covers_the_reference_registration_calibration_and_development_keys_only():
    raw = {'ref|a': np.arange(3), 'reg|a': np.arange(2), 'cal|a': np.arange(4), 'dev|a': np.arange(5), 'test|a': np.arange(6), 'devpath|a': np.arange(7), 'freeze|x': np.arange(1)}
    h = R.freeze_hash2(raw)
    raw2 = dict(raw)
    raw2['test|a'] = np.arange(60)
    raw2['devpath|a'] = np.arange(70)
    raw2['freeze|x'] = np.arange(10)
    raw2['cal|T1B|law'] = np.arange(3)                              # made during the test run: not part of the freeze
    assert R.freeze_hash2(raw2) == h
    raw2['cal|D1devB|law'] = np.arange(3)                           # made during the development: part of it
    assert R.freeze_hash2(raw2) != h
    del raw2['cal|D1devB|law']
    for k in ('ref|a', 'reg|a', 'cal|a', 'dev|a'):
        raw3 = dict(raw)
        raw3[k] = raw[k] + 1
        assert R.freeze_hash2(raw3) != h


def test_the_test_run_is_refused_when_the_frozen_choices_changed():
    raw = {'ref|a': np.arange(3), 'dev|a': np.arange(5)}
    raw['freeze|sha256'] = np.array(R.freeze_hash2(raw))
    raw['dev|a'] = raw['dev|a'] + 1
    with pytest.raises(RuntimeError):
        R.step_test(raw, {}, log=lambda s: None)
    with pytest.raises(RuntimeError):
        R.step_test({'ref|a': np.arange(3)}, {}, log=lambda s: None)


def test_raw_round_trip_is_deterministic(tmp_path, monkeypatch):
    raw = {'a': np.arange(5.0), 'ref|gates': R.jdump({'x': 1})}
    p1, p2 = str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz')
    R.save_raw(raw, p1)
    R.save_raw(dict(reversed(list(raw.items()))), p2)
    assert open(p1, 'rb').read() == open(p2, 'rb').read()
    back = R.load_raw(p1)
    assert np.array_equal(back['a'], raw['a']) and R.jload(back['ref|gates']) == {'x': 1}


def test_the_grids_are_the_planned_ones():
    assert R.KAPPAS2 == (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0)
    assert R.FW_GRID == (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
    assert R.SMIN_GRID == (None, 1.0, 1.02, 1.05, 1.1, 1.15, 1.2, 1.3) and R.SMAX_GRID == (1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.5, 3.2)
    assert R.ISS_START == dict(kappa=5.0, fw=0.0, smin=None, smax=3.2) and R.ISS_ROUNDS == 3
    assert (R.K_MERGE, R.GAP_MERGE, R.B_BLOCK) == (3, 7, 21) and R.RM_REL == 0.5 and R.CLOSE_REL == 2.0 and R.BORDER_NM == 200.0
    assert R.CFG_SETS['D1'][4:] == ('hp3h_s12', 'hp3h_s17', 'hp3h_r15', 'hp3h_r25') and R.CFG_SETS['D2'][4:] == ('hp5_s12', 'hp5_s17', 'hp5_r15', 'hp5_r25')
    assert R.TESTS['T1']['frames'] == (15001, 17000) and R.TESTS['T2a']['frames'] == (1, 1055) and R.TESTS['T2b']['frames'] == (301, 1054)
    assert R.TESTS['T1']['ref_ranges'] == [(2001, 10000), (12001, 15000), (17001, 20000)]
    assert (R.R1_MAX_DEV_NM, R.R1_MIN_Q, R.R2_MED_NM, R.R2_P90_NM, R.REG_BLOCK) == (20.0, 2.5, 15.0, 30.0, 125)
    assert R.MARGIN_E == 3.0 and R.MARGIN_J == 0.03 and R.BOOT_BLOCK_B == 10


# ----------------------------------------------------------------------------------------------------------------------------------
# more: the ISS call, the development scorer, the gate R2, the validity record
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_iss_calls_the_guarded_path_with_the_requested_order_and_reads_every_requested_kappa(monkeypatch):
    from src.core import imaging as IM
    from src.core import imaging_guard as IG
    orders, kappas = [], []
    orig_path, orig_stop = IG.sliding_path_guarded, IM.stop_snr

    def path(*a, **kw):
        orders.append((kw.get('order'), kw.get('cycle_guard'), kw.get('snr_min')))
        return orig_path(*a, **kw)

    def stop(hist, taus, snr, kappa=4.5):
        kappas.append(kappa)
        return orig_stop(hist, taus, snr, kappa)
    monkeypatch.setattr(IG, 'sliding_path_guarded', path)
    monkeypatch.setattr(IM, 'stop_snr', stop)
    rng = np.random.default_rng(13)
    P = IM.PixelPSF(40, 40, 1.5, parts=((1.0, 1.0),))
    mean = 300.0 + P.synth(np.array([[20.0, 20.0]]), np.array([3000.0]))
    frames = np.stack([mean + np.sqrt(0.8 * (mean - 200.0) + 100.0) * rng.standard_normal((40, 40))])
    cfg = G.grid_cfg('hp5', 1.5)
    cal = dict(bg=np.full((40, 40), 300.0), a=0.8, b0=100.0, x_ref=200.0)
    for order in ('dual', 'greedy'):
        R.iss_run2(frames, cal, cfg, 100.0, kappas=(5.0, 20.0), order=order, log=lambda s: None)
    assert orders == [('dual', True, 4.0), ('greedy', True, 4.0)]
    assert kappas == [5.0, 20.0, 5.0, 20.0]


def test_make_ev_e_scores_a_point_through_the_filters_caches_and_refuses_few_detections(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    r = rows['ISS4']
    nfr = 60
    dets = {k: [np.zeros((0, 4)) for _ in range(nfr)] for k in R.KAPPAS2}
    amp = np.where(np.arange(len(r.frame)) % 2 == 0, 1000.0, 3000.0)
    per = [np.c_[r.xy[r.frame == i], amp[r.frame == i], np.full((r.frame == i).sum(), 1.4)] for i in range(nfr)]
    dets[10.0] = per
    dets[20.0] = [d[:1] if i < 10 else d[:0] for i, d in enumerate(per)]        # ten detections only
    ev, w_ref = R.make_ev_E(dets, sc, 900.0, (0.0, 0.0), 1.2)
    assert w_ref == np.median(amp)
    full = ev(dict(kappa=10.0, fw=0.0, smin=None, smax=3.2))
    half = ev(dict(kappa=10.0, fw=1.25, smin=None, smax=3.2))                 # drops the amplitude-1000 half of the detections (and keeps those at >= 1.25 w_ref only if larger)
    assert np.isfinite(full) and half != full
    assert ev(dict(kappa=20.0, fw=0.0, smin=None, smax=3.2)) == -np.inf
    assert ev(dict(kappa=10.0, fw=0.0, smin=None, smax=1.1)) == -np.inf        # all widths 1.4 > 1.1 x 1.2


def test_gate_r2_measures_the_residual_translation_of_each_row_onto_the_table(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, n_frames=250, seed=4)
    pub = rows['PUB']
    for off, ok in ((10.0, True), (120.0, False)):
        shifted = {k: G.Loc(pub.frame, pub.xy + np.array([off, 0.0])) for k in ('R0', 'R2E', 'ISS4')}
        rr = dict(rows)
        rr.update(shifted)
        monkeypatch.setattr(R, 'set_rows', lambda raw_, old, name, register=True, curve=True, rr=rr: rr)
        monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, 250), ref_ranges=[(1, 1)], group='D1'))
        out = R.gate_r2(raw, None, 'S1')
        assert set(out) == {'R0', 'R2E', 'ISS4'} and all(abs(v['median'] - off) < 21.0 for v in out.values())
        assert all(v['passed'] == ok for v in out.values()) and all(v['n_blocks'] == 1 for v in out.values())


def test_set_valid_requires_v1_v4_r1_and_r2():
    gates = {'S1': dict(V1_abs_diff=1.9, scored_pores=100)}
    raw = {'ref|gates': R.jdump(gates), 'reg|S1|resid': np.zeros((3, 2)), 'reg|S1|blocks': np.array([[0, 1, 0, 0, 5.0]] * 3)}
    r2 = {'R0': dict(passed=True), 'R2E': dict(passed=True), 'ISS4': dict(passed=True)}
    assert R.set_valid(raw, None, 'S1', r2)['valid']
    assert not R.set_valid(raw, None, 'S1', dict(r2, ISS4=dict(passed=False)))['valid']
    gates['S1']['V1_abs_diff'] = 2.1
    raw['ref|gates'] = R.jdump(gates)
    assert not R.set_valid(raw, None, 'S1', r2)['valid']
    gates['S1'].update(V1_abs_diff=1.0, scored_pores=99)
    raw['ref|gates'] = R.jdump(gates)
    v = R.set_valid(raw, None, 'S1', r2)
    assert not v['valid'] and not v['V4'] and v['V1'] and v['R1']
    gates['S1']['scored_pores'] = 100
    raw['ref|gates'] = R.jdump(gates)
    raw['reg|S1|resid'] = np.full((3, 2), 25.0)
    assert not R.set_valid(raw, None, 'S1', r2)['R1']
    assert 'R2' not in R.set_valid(raw, None, 'S1')


def test_common_radius_is_the_mean_of_the_two_d2_fields_and_the_own_radius_for_d1():
    raw = {'ref|gates': R.jdump({'T1': dict(R_ref=87.0), 'T2a': dict(R_ref=53.0), 'T2b': dict(R_ref=54.0)})}
    assert R.common_R(raw, 'T1') == 87.0 and R.common_R(raw, 'T2a') == 53.5 and R.common_R(raw, 'T2b') == 53.5


def test_expected_yield_is_the_reference_ring_yield_scaled_to_the_scored_frames(monkeypatch):
    rng = np.random.default_rng(0)
    sp = G.DS['D1']
    cen = np.array([(6500.0 + 1000.0 * i, 1500.0 + 1000.0 * j) for i in range(3) for j in range(3)])
    sc = R.Scorer2(cen, R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0), x_min=6200.0)
    pts = np.concatenate([ring_points(rng, c, 600, sd=12.0) for c in cen])
    tab = dict(frame=np.arange(len(pts)) % 100 + 1, x=pts[:, 0], y=pts[:, 1], ph=np.ones(len(pts)))
    monkeypatch.setattr(R, 'load_tab', lambda ds, field=None: tab)
    monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, 20), ref_ranges=[(1, 100)], group='D1'))
    n_exp, f = R.expected_yield2({}, 'S1', sc)
    assert abs(f['N'] / len(pts) - 1.0) < 0.05 and abs(n_exp / (f['N'] * 0.2) - 1.0) < 1e-9
    monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, 20), ref_ranges=[(1, 50)], group='D1'))
    n_exp2, f2 = R.expected_yield2({}, 'S1', sc)
    assert abs(n_exp2 / (f2['N'] * 20 / 50) - 1.0) < 1e-9 and f2['N'] < f['N']


def test_the_polynomial_degree_stays_below_the_number_of_blocks_minus_one_so_that_an_outlier_leaves_a_residual():
    n = 500
    blocks = R.reg_blocks(n)
    assert len(blocks) == 4
    sh = np.c_[np.array([0.0, 10.0, 20.0, 30.0]), np.zeros(4)]
    sh[2, 0] += 40.0
    coef, res = R.fit_shift_curve(blocks, sh, n)
    assert coef[0, 0] == 0.0                                  # degree 2 at most with four blocks
    assert np.abs(res[:, 0]).max() > 10.0


def test_shift_of_evaluates_the_curve_at_the_frame_centre():
    raw = {'reg|X|coef': np.array([[0.0, 0.0, 100.0, 7.0], [0.0, 0.0, 0.0, 1.0]])}
    s = R.shift_of(raw, 'X', n=4)
    assert np.allclose(s[:, 0], 7.0 + 100.0 * (np.arange(4) + 0.5) / 4) and np.allclose(s[:, 1], 1.0)


def test_the_bootstrap_draws_every_block_including_the_last(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    last = rows['ISS4']
    keep = last.frame >= 50
    rows['ISS4'] = G.Loc(last.frame[keep], last.xy[keep])                      # all of its localizations are in the last of the six blocks
    a = R.analyse_A(raw, None, 'S1', n_boot=60, with_curve=False)
    d = a['rows']['ISS4']
    assert d['ci']['N'][1] > 0.5 * d['fit']['N']


def test_gate_r2_fails_on_the_ninetieth_percentile_even_when_the_median_passes(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, n_frames=1000, seed=5)
    pub = rows['PUB']
    off = np.where(pub.frame >= 750, 60.0, 3.0)
    row = G.Loc(pub.frame, pub.xy + np.c_[off, np.zeros(len(off))])
    rr = dict(rows)
    rr.update({k: row for k in ('R0', 'R2E', 'ISS4')})
    monkeypatch.setattr(R, 'set_rows', lambda raw_, old, name, register=True, curve=True: rr)
    monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, 1000), ref_ranges=[(1, 1)], group='D1'))
    out = R.gate_r2(raw, None, 'S1')
    for v in out.values():
        assert v['n_blocks'] == 4 and v['median'] < 15.0 and v['p90'] > 30.0 and not v['passed']


def test_reg_loc_subtracts_the_frame_shift_from_the_localizations():
    loc = G.Loc(np.array([0, 1, 1]), np.array([[100.0, 100.0], [200.0, 200.0], [300.0, 300.0]]), np.array([[1.0], [2.0], [3.0]]))
    shifts = np.array([[10.0, -5.0], [20.0, 0.0]])
    out = R.reg_loc(loc, shifts)
    assert np.allclose(out.xy, [[90.0, 105.0], [180.0, 200.0], [280.0, 300.0]]) and out.frame.tolist() == [0, 1, 1] and out.extra.shape == (3, 1)


def test_the_continuity_row_applies_the_old_width_limit_in_units_of_sigma0_over_1_4(monkeypatch):
    monkeypatch.setitem(R.TESTS, 'S1', dict(ds='D1', field=None, frames=(1, 2), ref_ranges=[(1, 1)], group='D1'))
    dets = {k: [np.array([[100.0, 100.0, 5.0, 1.0], [200.0, 200.0, 5.0, 3.0]]), np.zeros((0, 4))] for k in R.KAPPAS2}
    raw = {'dev|D1|iss4': R.jdump(dict(cfg='A')), 'cal|S1|sigma0': R.jdump({'sigma0': 1.5}), 'cal|S1|law': np.arange(9.0) + 1.0, 'cal|S1|bg': np.zeros((2, 2))}
    R.put_paths(raw, 'testpath|S1|iss4', dets, 2, np.zeros(2), np.zeros(2))
    old = {'dev|D1|dual': R.jdump(dict(cfg='A', kappa=10.0, smax_rel=2.0))}
    loc = R.iss_old_loc(raw, old, 'S1')
    assert len(loc.frame) == 1 and abs(loc.extra[0, 1] - 1.0) < 1e-12                  # limit 2 x 1.5 / 1.4 = 2.14: the width 3.0 is dropped, 1.0 kept
    old = {'dev|D1|dual': R.jdump(dict(cfg='B', kappa=10.0, smax_rel=2.0))}                # another configuration and no separate path: no row
    assert R.iss_old_loc(raw, old, 'S1') is None
    old = {'dev|D1|dual': R.jdump(dict(cfg='A', kappa=4.5, smax_rel=2.0))}                 # a kappa outside the grid needs its own path
    assert R.iss_old_loc(raw, old, 'S1') is None


def test_pooled_areas_are_summed_over_the_members(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    monkeypatch.setitem(R.TESTS, 'S2', dict(ds='D1', field=None, frames=(1, 60), ref_ranges=[(1, 1)], group='D1'))
    monkeypatch.setitem(R.GROUPS, 'SS', ['S1', 'S2'])
    sp = G.DS['D1']
    cen = np.array([(6500.0 + 1000.0 * i, 1500.0 + 1000.0 * j) for i in range(4) for j in range(4)])
    sc2 = R.Scorer2(cen[:4], R0, (0.0, 0.0, 12100.0, 12100.0), sp, (0.0, 0.0), x_min=6200.0)
    monkeypatch.setattr(R, 'scorer_of', lambda raw_, name, old=None: sc if name == 'S1' else sc2)
    raw['ref|gates'] = np.array(json.dumps({'S1': dict(n_exp=900.0), 'S2': dict(n_exp=400.0)}))
    raw['cal|S2|sigma0'] = np.array(json.dumps({'sigma0': 1.5}))
    one = R.analyse_A(raw, None, 'S1', n_boot=2, with_curve=False)
    two = R.analyse_A(raw, None, 'SS', n_boot=2, with_curve=False)
    f1, f2 = one['rows']['R0']['fit']['floor'], two['rows']['R0']['fit']['floor']
    # the second member scores the same localizations on a quarter of the pores: its floor is about a quarter of the first one's, so the pooled floor is about 1.25 times
    assert 1.1 < f2 / f1 < 1.4 and two['n_pores'] == sc.n_pores + sc2.n_pores


# ----------------------------------------------------------------------------------------------------------------------------------
# the wording round after the referee report A+B: (void) readings, the transfer of the tuning, the continuity row, the grid edges
# ----------------------------------------------------------------------------------------------------------------------------------
def test_with_void_marks_the_reading_of_a_dataset_whose_gates_failed():
    assert R.with_void('equivalent', False) == 'equivalent (void)' and R.with_void('Picasso better', False) == 'Picasso better (void)'
    assert R.with_void('ISS better', True) == 'ISS better'


def test_the_pair_list_of_part_a_has_the_continuity_row_against_both_picasso_points(monkeypatch):
    raw, sc, rows = synth_set(monkeypatch, better=True)
    a = R.analyse_A(raw, None, 'S1', n_boot=10, with_curve=False)
    assert list(a['pairs']) == ['ISS4-R2E', 'ISS4-R0', 'R2E-R0', 'ISS_OLD-R2E', 'ISS_OLD-R0']
    d = a['pairs']['ISS_OLD-R0']
    assert d['dE'] == pytest.approx(a['rows']['ISS_OLD']['E'] - a['rows']['R0']['E']) and d['dE_ci'][0] <= d['dE'] <= d['dE_ci'][1]
    assert d['dJ'] == pytest.approx(a['rows']['ISS_OLD']['J'] - a['rows']['R0']['J'])


def test_edge_of_names_the_lower_edge_the_top_edge_and_the_interior_of_the_kappa_grid():
    assert R.edge_of(R.KAPPAS2[0]) == 'at the lower edge' and R.edge_of(R.KAPPAS2[-1]) == 'at the top edge' and R.edge_of(R.KAPPAS2[3]) == 'in the interior'


def test_rank_str_orders_largest_first():
    assert R.rank_str({'a': 1.0, 'b': 3.0, 'c': 2.0}) == 'b 3.0 > c 2.0 > a 1.0'


def fake_res(dE, ci, e_rows, extra=None):
    pairs = {k: dict(dE=v, dE_ci=[v - 1.0, v + 1.0]) for k, v in dE.items()}
    return dict(pairs=pairs, rows={k: dict(E=v) for k, v in e_rows.items()})


def fake_dev(raw_old, raw_new, g, r0, r2e, old_iss, iss4):
    raw_old['dev|%s|picasso' % g] = np.array(json.dumps(dict(R0=dict(E=r0), R2E=dict(E=r2e))))
    raw_old['dev|%s|dual' % g] = np.array(json.dumps(dict(E=old_iss)))
    raw_new['dev|%s|iss4' % g] = np.array(json.dumps(dict(E=iss4)))


def test_the_transfer_paragraph_says_that_untuned_picasso_beats_tuned_picasso_and_reverses_the_ranking():
    old, raw = {}, {}
    fake_dev(old, raw, 'D1', 45.1, 50.0, 44.1, 50.3)
    fake_dev(old, raw, 'D2', 12.1, 19.4, 11.7, 18.7)
    res = {'D1': fake_res({'R2E-R0': -4.7, 'ISS4-R0': -4.2, 'ISS4-R2E': 0.5}, None, {'R0': 31.1, 'R2E': 26.4, 'ISS4': 26.9}),
           'D2': fake_res({'R2E-R0': -23.9, 'ISS4-R0': -34.4, 'ISS4-R2E': -10.4, 'ISS_OLD-R2E': 43.5, 'ISS_OLD-R0': 19.6}, None, {'R0': 34.3, 'R2E': 10.3, 'ISS4': -0.1, 'ISS_OLD': 53.8})}
    out = R.transfer_paragraph(raw, old, res)
    assert len(out) == 4
    assert out[0].startswith('**Tuning on the development frames did not transfer.** On the test frames untuned Picasso (R0) beats Picasso tuned on the development frames (R2E) on both datasets (paired Delta E_s, R2E - R0: D1 -4.7 [-5.7, -3.7]; D2 -23.9 [-24.9, -22.9]).')
    assert 'both worse than untuned Picasso (ISS4 - R0: D1 -4.2 [-5.2, -3.2]; D2 -34.4 [-35.4, -33.4]).' in out[0]
    assert out[1] == 'D1: development ranking (dev E_s) ISS4 50.3 > R2E 50.0 > R0 45.1 > first-study ISS 44.1; test ranking (E_s) R0 31.1 > ISS4 26.9 > R2E 26.4; the order of R0, R2E is exactly reversed between the development and the test frames.'
    assert out[2] == ('D2: development ranking (dev E_s) R2E 19.4 > ISS4 18.7 > R0 12.1 > first-study ISS 11.7; test ranking (E_s) first-study ISS 53.8 > R0 34.3 > R2E 10.3 > ISS4 -0.1; '
                      'the order of R0, R2E, first-study ISS is exactly reversed between the development and the test frames.')
    assert out[3].startswith('Continuity row (D2; POST HOC, the first study\'s frozen two-dial point, outside the registered rule and on a void pooled set, so not evidence of an advantage): ISS_OLD - R2E 43.5 [42.5, 44.5] and ISS_OLD - R0 19.6 [18.6, 20.6]')
    assert 'ISS4 - R2E -10.4 [-11.4, -9.4] the sign of the D2 contrast depends on which frozen ISS point is used.' in out[3]


def test_the_transfer_paragraph_does_not_claim_a_reversal_or_a_win_the_numbers_do_not_show():
    old, raw = {}, {}
    fake_dev(old, raw, 'D1', 45.1, 50.0, 44.1, 50.3)
    fake_dev(old, raw, 'D2', 12.1, 19.4, 11.7, 18.7)
    res = {'D1': fake_res({'R2E-R0': -0.5, 'ISS4-R0': -0.5, 'ISS4-R2E': 1.0}, None, {'R0': 20.0, 'R2E': 22.0, 'ISS4': 23.0}),             # intervals [-1.5, 0.5]: they straddle 0
           'D2': fake_res({'R2E-R0': -5.0, 'ISS4-R0': -1.0, 'ISS4-R2E': 4.0}, None, {'R0': 30.0, 'R2E': 25.0, 'ISS4': 29.0})}
    out = R.transfer_paragraph(raw, old, res)
    assert len(out) == 3                                                     # no ISS_OLD row: no continuity sentence
    assert 'untuned Picasso (R0) does not beat Picasso tuned on the development frames (R2E) on both datasets' in out[0] and 'both worse than' not in out[0] and 'not always worse than untuned Picasso' in out[0]
    assert 'the order of R0, R2E is not reversed between' in out[1]


def test_the_grid_edge_sentence_names_the_part_b_optimum_at_the_lower_edge_and_the_part_a_count_at_the_top():
    t = R.grid_edges_line({'D1': 0, 'D2': 6}, {'D1': 4.0, 'D2': 6.0})
    assert t == ('Grid edges of the ISS tuning (development): Part A, the ISS4 optimum is at the top kappa (50) for 0 of 8 configurations of D1 and 6 of 8 of D2; Part B, the ISS-B optimum is at kappa 4 for D1 '
                 '(at the lower edge of the kappa grid 4-50) and at kappa 6 for D2 (in the interior).')
    assert 'for D2 (at the top edge).' in R.grid_edges_line({'D1': 1, 'D2': 1}, {'D1': 4.0, 'D2': 50.0})
