"""Tests of the pure functions of experiments/exp_smlm_real_gt.py (the pre-registered real-data test on nuclear pores with a known geometry): the pore finder, the circle
fit and the 8-fold phase on synthetic rings, the pore model (assignment, areas, radial bins), the radial-profile fit, J_s and E_s, the coordinate ascent, Picasso's
post-filter semantics, the ISS / greedy path on a small synthetic frame, the equal-yield thinning, the paired bootstrap on a synthetic raw, the reading of the outcomes
under the plan's decision rule, and the report's determinism. No downloads, no Picasso and no tifffile needed. About a minute."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_smlm_real_gt as G      # noqa: E402

R0 = 87.0


def ring_points(rng, centre, n, R=R0, sd=10.0, floor=0):
    th = rng.uniform(0, 2 * np.pi, n)
    r = R + sd * rng.standard_normal(n)
    pts = np.c_[centre[0] + r * np.cos(th), centre[1] + r * np.sin(th)]
    if floor:
        pts = np.r_[pts, centre + rng.uniform(-300, 300, (floor, 2))]
    return pts


def synth_raw(n_pores_side=6, spacing=1000.0):
    """A synthetic raw with the keys the Scorer reads for 'D1test': isolated pores on a grid (no two discs overlap), a zero convention shift."""
    xs = 6500.0 + spacing * np.arange(n_pores_side)
    ys = 1500.0 + spacing * np.arange(n_pores_side)
    cen = np.array([(x, y) for x in xs for y in ys])
    left = np.array([(x, y) for x in 1500.0 + spacing * np.arange(3) for y in ys])           # development pores (x < 5800)
    cen = np.r_[left, cen]
    raw = {'pores|D1|centres': cen, 'gates': np.array(json.dumps({'D1': {'R_ref': R0}})), 'cal|D1|convention_offset': np.array([0.0, 0.0, 10.0])}
    return raw, cen


def test_ring_candidates_and_find_pores_recover_planted_rings():
    rng = np.random.default_rng(1)
    centres = np.array([(1000.0, 1000.0), (1800.0, 1100.0), (1200.0, 2000.0), (2300.0, 2200.0)])
    pts = np.concatenate([ring_points(rng, c, 400, sd=8.0) for c in centres])
    sp = G.DS['D1']
    cen, rad, cnt = G.find_pores(pts[:, 0], pts[:, 1], (0.0, 0.0, 3200.0, 3200.0), sp)
    assert len(cen) == 4
    d = np.hypot(*(cen[:, None, :] - centres[None]).transpose(2, 0, 1)).min(axis=1)
    assert d.max() < 8.0
    assert np.all(np.abs(rad - R0) < 4.0)
    assert np.all(cnt >= G.MIN_ANNULUS)


def test_find_pores_merges_close_centres_and_rejects_off_radius_candidates():
    rng = np.random.default_rng(2)
    c = np.array([1500.0, 1500.0])
    pts = np.r_[ring_points(rng, c, 500), ring_points(rng, c + [60.0, 0.0], 500)]       # two overlapping copies: one pore
    cen, rad, cnt = G.find_pores(pts[:, 0], pts[:, 1], (0.0, 0.0, 3000.0, 3000.0), G.DS['D1'])
    assert len(cen) == 1


def test_circle_fit_recovers_centre_and_radius_and_refuses_a_sparse_annulus():
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(3)
    c = np.array([700.0, 900.0])
    pts = ring_points(rng, c, 300, R=60.0, sd=3.0)
    t = cKDTree(pts)
    f = G.circle_fit(t, pts[:, 0], pts[:, 1], c[0] + 15, c[1] - 12, 25.0, 80.0, 90.0)
    assert f is not None and abs(f[0] - c[0]) < 3 and abs(f[1] - c[1]) < 3 and abs(f[2] - 60.0) < 3
    assert G.circle_fit(t, pts[:, 0], pts[:, 1], 5000.0, 5000.0, 25.0, 80.0, 90.0) is None


def test_pore_phases_recover_a_planted_eight_fold_phase_and_have_no_phase_for_a_uniform_ring():
    rng = np.random.default_rng(4)
    sp = G.DS['D2']
    c = np.array([[0.0, 0.0]])
    phi = 0.21                                                                  # radians, inside [0, 45 degrees)
    corners = phi + np.arange(8) * np.pi / 4
    th = rng.choice(corners, 600) + 0.08 * rng.standard_normal(600)
    pts = np.c_[53.0 * np.cos(th), 53.0 * np.sin(th)]
    ph = G.pore_phases(pts[:, 0], pts[:, 1], c, 53.0, sp)
    assert abs(ph[0, 0] - phi) < 0.03 and ph[0, 2] < 1e-6
    u = rng.uniform(0, 2 * np.pi, 600)
    ph2 = G.pore_phases(53.0 * np.cos(u), 53.0 * np.sin(u), c, 53.0, sp)
    assert ph2[0, 1] < 0.12 and ph2[0, 2] > 1e-3


def test_block_split_deals_blocks_of_100_reference_frames_in_order_and_table_select_is_inclusive():
    frames = np.arange(1, 801)
    tab = dict(frame=frames, x=frames * 1.0, y=frames * 2.0)
    ref = [(101, 400), (601, 800)]
    odd = G.block_split(tab, ref, 1)
    even = G.block_split(tab, ref, 0)
    assert len(odd['frame']) == 200 and len(even['frame']) == 300 and set(odd['frame']).isdisjoint(even['frame'])
    assert set(odd['frame']) == set(range(201, 301)) | set(range(601, 701))                     # blocks 1 and 3 of the 500 reference frames
    assert set(even['frame']) == set(range(101, 201)) | set(range(301, 401)) | set(range(701, 801))
    sel = G.table_select(tab, [(10, 12), (20, 20)])
    assert list(sel['frame']) == [10, 11, 12, 20]


def test_disc_inside_requires_the_whole_disc():
    cen = np.array([[300.0, 300.0], [200.0, 300.0], [300.0, 900.0]])
    m = G.disc_inside(cen, 87.0, 150.0, (0.0, 0.0, 1200.0, 1200.0))
    assert list(m) == [True, False, True]


def test_convention_is_the_multiple_of_half_a_pixel_nearest_to_the_measured_offset():
    raw = {'cal|D1|convention_offset': np.array([-3.8, -5.3, 8862.0]), 'cal|D2|convention_offset': np.array([-52.5, -52.7, 2552.0])}
    assert np.allclose(G.convention(raw, 'D1'), [0.0, 0.0])
    assert np.allclose(G.convention(raw, 'D2'), [-53.0, -53.0])
    raw['cal|D2|convention_offset'] = np.array([-80.0, 27.0, 10.0])
    assert np.allclose(G.convention(raw, 'D2'), [-106.0, 53.0])


def test_scorer_keeps_only_scored_pores_and_its_areas_add_up():
    raw, cen = synth_raw()
    sc = G.Scorer(raw, 'D1test')
    assert sc.n_pores == 36 and sc.scored.sum() == 36
    # isolated pores: every cell is the whole disc
    assert abs(sc.area.sum() / (36 * np.pi * sc.rmax ** 2) - 1.0) < 0.02
    assert np.allclose(sc.area_pore.sum(0), sc.area)
    # a localization on a development pore (x < 5800) is not in a scored disc, one on a test pore is
    r, j, keep = sc.assign(np.array([[1500.0 + R0, 1500.0], [6500.0 + R0, 1500.0], [6500.0, 1500.0 + 2 * R0 + 400]]))
    assert list(keep) == [False, True, False] and abs(r[1] - R0) < 1e-9
    d = G.Scorer(raw, 'D1dev')
    assert d.n_pores == 18


def test_scorer_drops_pores_whose_disc_leaves_the_region_after_the_shift():
    raw, cen = synth_raw()
    raw['cal|D1|convention_offset'] = np.array([0.0, 0.0, 10.0])
    base = G.Scorer(raw, 'D1test').n_pores
    c2 = np.r_[cen, [[6500.0, 100.0]]]                                       # a pore 100 nm from the border: its disc (237 nm) is outside
    raw['pores|D1|centres'] = c2
    assert G.Scorer(raw, 'D1test').n_pores == base


def test_fit_profile_recovers_a_planted_ring_and_floor_and_a_warm_start_agrees():
    raw, _ = synth_raw()
    sc = G.Scorer(raw, 'D1test')
    N, mu, w, b = 8000.0, R0 - 3.0, 14.0, 4e-5
    e = N * (G.ndtr((sc.edges[1:] - mu) / w) - G.ndtr((sc.edges[:-1] - mu) / w)) + b * sc.area
    counts = np.random.default_rng(5).poisson(e).astype(float)
    f = G.fit_profile(counts, sc.area, sc.edges, sc.R)
    assert abs(f['N'] / N - 1) < 0.05 and abs(f['mu'] - mu) < 1.0 and abs(f['w'] / w - 1) < 0.07 and abs(f['floor'] / (b * sc.area.sum()) - 1) < 0.2
    g = G.fit_profile(counts, sc.area, sc.edges, sc.R, start=f['x'])
    assert abs(g['w'] - f['w']) < 1e-3


def test_e_score_is_the_challenge_efficiency_with_the_radial_width_and_the_missed_yield():
    f = dict(N=800.0, floor=100.0, w=20.0)
    J, E = G.e_score(f, 1000.0)                                              # missed 200
    assert abs(J - 800.0 / 1100.0) < 1e-12
    assert abs(E - (100.0 - np.hypot(100.0 - 100.0 * J, 20.0))) < 1e-9
    J2, _ = G.e_score(f, 500.0)                                              # more than expected: nothing missed
    assert abs(J2 - 800.0 / 900.0) < 1e-12


def test_coordinate_ascent_moves_only_on_a_strict_improvement_and_keeps_the_earlier_value_on_ties():
    grids = {'a': [None, 1, 2, 3], 'b': [None, 10, 20]}
    calls = []

    def ev(f):
        calls.append(dict(f))
        return -abs((f['a'] or 0) - 2) - abs((f['b'] or 0) - 20) / 10.0
    cur, best, rounds = G.ascent(ev, grids, 3)
    assert cur == {'a': 2, 'b': 20} and best == 0.0 and rounds == 2
    cur2, best2, _ = G.ascent(lambda f: 1.0, grids, 3)                      # a constant objective: nothing moves
    assert cur2 == {'a': None, 'b': None}


def test_picasso_post_filters_are_strict_inequalities_with_both_widths():
    class S:
        nb = 10
        edges = np.arange(11) * 5.0

        def assign(self, xy):
            n = len(xy)
            return np.full(n, 10.0), np.zeros(n, int), np.ones(n, bool)
    a = np.array([[0.0, 0.0, 100.0, 150.0, 150.0, 1000.0],                    # x, y, photons, sx, sy, net gradient
                  [0.0, 0.0, 100.0, 150.0, 250.0, 2000.0],
                  [0.0, 0.0, 300.0, 150.0, 150.0, 2000.0]])
    cand = G.PicCand([a[:1], a[1:2], a[2:]], (0.0, 0.0), S())
    assert list(cand.mask({})) == [True, True, True]
    assert list(cand.mask({'min_net_gradient': 1000.0})) == [False, True, True]
    assert list(cand.mask({'min_photons': 100.0})) == [False, False, True]
    assert list(cand.mask({'max_sigma_nm': 250.0})) == [True, False, True]       # |sy| < 250 is strict
    assert list(cand.mask({'min_sigma_nm': 150.0})) == [False, False, False]


def test_grid_cfg_scales_widths_and_the_variants_change_one_thing_each():
    s0 = 1.6
    base = G.grid_cfg('hp5h', s0)
    s = s0 / 1.4
    assert abs(base['sigma0'] - s0) < 1e-12 and abs(base['r_ex'] - 2.0 * s) < 1e-12 and abs(base['widths'][1] - 3.2 * s) < 1e-12
    assert base['parts'] == [[0.75, 1.0], [0.25, 3.0]] and base['hp'] == 5.0
    assert abs(G.grid_cfg('hp3_s12', s0)['sigma0'] - s0 * 1.2 / 1.4) < 1e-12 and G.grid_cfg('hp3_s12', s0)['hp'] == 3.0
    assert abs(G.grid_cfg('hp5_r25', s0)['r_ex'] - 2.5 * s) < 1e-12 and G.grid_cfg('hp5_r25', s0)['sigma0'] == s0
    assert abs(G.grid_cfg('hp5h_r15', s0)['r_ex'] - 1.5 * s) < 1e-12
    with pytest.raises(ValueError):
        G.grid_cfg('hp5h_x', s0)


def test_choose_op_takes_the_largest_E_then_the_larger_kappa_then_the_smaller_width():
    rows = [(4.0, 2.0, 50.0, .5, 1), (5.0, 2.0, 50.0, .5, 1), (5.0, 2.5, 50.0, .5, 1), (6.0, 3.2, 49.0, .5, 1)]
    assert G.choose_op(rows) == (5.0, 2.0, 50.0, .5, 1)
    assert G.post_width([np.array([[0, 0, 1, 1.0], [0, 0, 1, 2.1]])], 2.0)[0].shape[0] == 1


def synth_frame(rng, c, ny=40, nx=40, sigma=1.5, photons=3000.0, bg=300.0, a=0.8, b0=100.0):
    from src.core import imaging as IM
    P = IM.PixelPSF(ny, nx, sigma, parts=((1.0, 1.0),))
    mean = bg + P.synth(np.atleast_2d(c), np.array([photons]))
    return mean + np.sqrt(a * (mean - 200.0) + b0) * rng.standard_normal((ny, nx))


@pytest.mark.parametrize('order', ['dual', 'greedy'])
def test_iss_and_greedy_paths_find_a_planted_emitter_in_noise_and_nothing_in_an_empty_frame(order):
    rng = np.random.default_rng(6)
    c = np.array([[20.3, 18.6]])
    frames = np.stack([synth_frame(rng, c) for _ in range(2)] + [300.0 + np.sqrt(0.8 * 100.0 + 100.0) * rng.standard_normal((40, 40))])
    cfg = G.grid_cfg('hp5', 1.5)
    cal = dict(bg=np.full((40, 40), 300.0), a=0.8, b0=100.0, x_ref=200.0)
    dets, secs, ent = G.iss_run(frames, cal, cfg, order, 100.0, log=lambda s: None)
    for i in (0, 1):
        d = dets[5.0][i]
        assert len(d) >= 1
        k = np.argmin(np.hypot(d[:, 0] - 100.0 * (c[0, 0] + 0.5), d[:, 1] - 100.0 * (c[0, 1] + 0.5)))
        assert np.hypot(d[k, 0] - 100.0 * (c[0, 0] + 0.5), d[k, 1] - 100.0 * (c[0, 1] + 0.5)) < 30.0          # nm
        assert abs(d[k, 2] / 3000.0 - 1) < 0.25
    assert len(dets[10.0][2]) == 0
    assert set(dets) == set(G.KAPPAS) and len(secs) == 3


def test_eq_yield_thinning_leaves_the_smaller_set_alone_and_matches_the_yield():
    raw, _ = synth_raw()
    sc = G.Scorer(raw, 'D1test')
    rng = np.random.default_rng(7)

    def counts(N, w, seed):
        e = N * (G.ndtr((sc.edges[1:] - R0) / w) - G.ndtr((sc.edges[:-1] - R0) / w)) + 2e-5 * sc.area
        return np.random.default_rng(seed).poisson(e).astype(float)
    Ha, Hb = counts(12000.0, 12.0, 1), counts(3000.0, 25.0, 2)
    fa, fb = (G.fit_profile(H, sc.area, sc.edges, sc.R) for H in (Ha, Hb))
    wa, wb = G.eq_yield_w(sc, fa, Ha, fb, Hb, sc.area, rng)
    assert wb == fb['w'] and abs(wa - 12.0) < 1.5                              # A is the larger set: thinned; its width does not depend on the yield
    wb2, wa2 = G.eq_yield_w(sc, fb, Hb, fa, Ha, sc.area, np.random.default_rng(7))
    assert abs(wa2 - wa) < 1e-9 and abs(wb2 - wb) < 1e-9                       # the order of the arguments does not matter


def test_verdict_readings_follow_the_plan_rule():
    def pr(dE, dw, dEp=None, dwp=None):
        return dict(dE_ci=dE, dw_ci=dw, dE_ci_pores=dEp or dE, dw_ci_pores=dwp or dw)
    assert G.verdict_ds(pr([1.0, 5.0], [-4.0, -1.0])) == 'ISS better'
    assert 'not claimed' in G.verdict_ds(pr([1.0, 5.0], [-4.0, -1.0], [-1.0, 5.0], [-4.0, 1.0]))
    assert G.verdict_ds(pr([-6.0, -1.0], [1.0, 4.0])) == 'Picasso better'
    assert 'more localizations' in G.verdict_ds(pr([1.0, 5.0], [-2.0, 1.0]))
    assert 'narrower' in G.verdict_ds(pr([-1.0, 4.0], [-4.0, -1.0]))
    assert G.verdict_ds(pr([-2.0, 2.0], [-1.0, 1.0])).startswith('equivalent')
    assert G.verdict_ds(pr([-9.0, 9.0], [-1.0, 1.0])).startswith('no detectable')
    b, w, n = 'ISS better', 'Picasso better', 'no detectable difference'
    assert G.overall_verdict(b, b, False).startswith('ISS better on BOTH')
    assert G.overall_verdict(w, n, False).startswith('DROP')
    assert G.overall_verdict(w, w, False).startswith('DROP')
    assert G.overall_verdict(b, w, False).startswith('MIXED: the ISS is better on one dataset and Picasso')
    assert G.overall_verdict(b, n, False).startswith('MIXED: the ISS is better on one dataset only')
    assert G.overall_verdict(n, n, True).startswith('NO DETECTABLE') and 'development gate' in G.overall_verdict(n, n, True)


def synth_test_raw(seed=8):
    """A synthetic raw with test rows for D1: ring localizations with a floor; ISS tighter and more numerous than R2E, the other rows in between."""
    raw, cen = synth_raw()
    rng = np.random.default_rng(seed)
    n_frames = 200
    raw['test|D1|n_frames'] = np.array(n_frames)
    raw['cal|D1test|sigma0'] = np.array(json.dumps(dict(sigma0=1.5)))
    sc = G.Scorer(raw, 'D1test')
    spec = {'R0': (2500, 24.0), 'R1E': (2800, 22.0), 'R2E': (3000, 20.0), 'ISS': (3600, 12.0), 'ISS0': (3300, 14.0), 'GREEDY': (3500, 13.0), 'PUB': (3200, 15.0)}
    for row, (npp, sd) in spec.items():
        pts = np.concatenate([ring_points(rng, c, int(npp / 36 * 1.0) + 1, sd=sd, floor=20) for c in sc.allc[sc.scored]])
        fr = rng.integers(0, n_frames, len(pts))
        G.put_loc(raw, 'test|D1|%s' % row, G.Loc(fr, pts), secs=np.full(n_frames, 0.5) if row != 'PUB' else None)
        if row in ('ISS', 'ISS0', 'GREEDY'):
            raw['test|D1|%s|entries' % row] = np.full(n_frames, 20.0)
    raw['test|D1|ref_fit'] = G.jdump(dict(n_exp=4000.0, fit=dict(N=1.0, mu=R0, w=15.0, b=1e-5, floor=1.0, n=1.0, ok=True, x=[0.0, 0.0, R0, 2.7])))
    return raw


def test_analyse_orders_the_rows_and_is_deterministic():
    raw = synth_test_raw()
    a1 = G.analyse(raw, 'D1', n_boot=40)
    a2 = G.analyse(raw, 'D1', n_boot=40)
    assert json.dumps(a1, sort_keys=True, default=float) == json.dumps(a2, sort_keys=True, default=float)
    w = {r: a1['rows'][r]['fit']['w'] for r in G.ROWS}
    assert w['ISS'] < w['GREEDY'] < w['R2E'] < w['R0']
    assert a1['rows']['ISS']['E'] > a1['rows']['R2E']['E'] > a1['rows']['R0']['E']
    pr = a1['pairs']['ISS-R2E']
    assert pr['dE'] > 0 and pr['dw'] < 0 and pr['dE_ci'][0] <= pr['dE'] <= pr['dE_ci'][1]
    assert pr['dw_ci'] is not None and a1['pairs']['ISS0-R2E']['dw_ci'] is None
    assert a1['n_pores'] == 36 and a1['nblocks'] == 20
    assert abs(a1['rows']['ISS']['fit']['mu'] - R0) < 2.0


def test_the_pore_list_hash_covers_the_pore_and_window_keys_only():
    raw, _ = synth_raw()
    raw['pores|D1|radii'] = np.ones(3)
    h = G.pore_list_hash(raw)
    raw['cal|x'] = np.zeros(2)
    assert G.pore_list_hash(raw) == h
    raw['pores|D1|radii'] = np.ones(3) * 2
    assert G.pore_list_hash(raw) != h


def test_the_test_run_is_refused_when_the_frozen_choices_changed():
    raw = {'freeze|sha256': np.array('0' * 64), 'pores|D1|centres': np.zeros((1, 2)), 'dev|D1|picasso': np.array('{}')}
    with pytest.raises(RuntimeError, match='refused'):
        G.step_test_all(raw)
    raw.pop('freeze|sha256')
    with pytest.raises(RuntimeError, match='refused'):
        G.step_test_all(raw)


def test_raw_round_trip_is_deterministic(tmp_path):
    raw = synth_test_raw()
    p1, p2 = str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz')
    G.save_raw(raw, p1)
    G.save_raw(G.load_raw(p1), p2)
    assert open(p1, 'rb').read() == open(p2, 'rb').read()


def synth_full_raw():
    """A synthetic raw complete enough for the report: D1 (synth_test_raw) and a D2 with 8-fold pores in its window, the development JSONs, the gates, the calibration."""
    raw = synth_test_raw()
    rng = np.random.default_rng(9)
    sp = G.DS['D2']
    R2 = 53.0
    xs = 2000.0 + 900.0 * np.arange(6)
    cen2 = np.array([(x, y) for x in xs for y in xs])
    phase = rng.uniform(0, np.pi / 4, len(cen2))
    for f in (17, 18):
        raw['pores|D2_f%d|centres' % f] = cen2
        raw['pores|D2_f%d|phase' % f] = np.c_[phase, np.full(len(cen2), 0.4), np.full(len(cen2), 1e-6), np.full(len(cen2), 100.0)]
    raw['window|D2dev'] = np.array([10, 10])
    raw['window|D2test'] = np.array([10, 10])
    raw['cal|D2|convention_offset'] = np.array([-52.5, -52.7, 2552.0])
    gates = {'D1': dict(R_ref=R0, n_pores=36, V1_abs_diff=0.07), 'V1': {'D1': True, 'D2_f17': True, 'D2_f18': True},
             'V2': dict(diff=0.7, lit_dev=0.006, lit_test=-0.007, **{'pass': True}), 'V3': {'D2_f17': 0.62, 'D2_f18': 0.64}}
    for f in (17, 18):
        gates['D2_f%d' % f] = dict(R_ref=R2, n_pores=len(cen2), V1_abs_diff=0.02, frac_8fold_significant=0.6)
    for sname in G.SETS:
        gates['V4_' + sname] = dict(n=100, ok=True)
    raw['gates'] = np.array(json.dumps(gates))
    raw['pores|sha256'] = np.array('0' * 64)
    raw['freeze|time_utc'], raw['freeze|sha256'] = np.array('2026-01-01T00:00:00Z'), np.array('1' * 64)
    raw['test|start_utc'], raw['test|end_utc'] = np.array('2026-01-02T00:00:00Z'), np.array('2026-01-03T00:00:00Z')
    raw['cal|N1|ng_median'], raw['cal|N1|ph_median'] = np.array(1000.0), np.array(80.0)
    for sname in G.SETS:
        raw['cal|%s|sigma0' % sname] = np.array(json.dumps(dict(sigma0=1.5, n_sel=100, n_all=300, fallback=False)))
        raw['cal|%s|law' % sname] = np.array([1.0, 100.0, 200.0, 0.0, 1.2, 1.0, 90.0, 120.0, 250.0])
    for ds in ('D1', 'D2'):
        fac = dict(f_ng=1.5, f_ph=4.0, f_sig=1.1)
        pic = dict(fac=fac, R0=dict(E=40.0), R1E=dict(E=41.0), R2E=dict(E=45.0, n_configs=11000, struct='box=7|fit_box=5', filters=dict(min_net_gradient=2000.0, min_photons=None)))
        d = dict(n_frames=300, E=44.0, cfg='hp5h', kappa=7.0, smax_rel=2.0, default_E=42.0, rows={'hp5h': [[7.0, 2.0, 44.0, .5, 1.0], [10.0, 2.0, 43.0, .5, 1.0]]})
        raw['dev|%s|picasso' % ds] = G.jdump(pic)
        raw['dev|%s|dual' % ds] = G.jdump(d)
        raw['dev|%s|greedy' % ds] = G.jdump(dict(d, E=43.0))
    # a D2 test set: ring localizations on the 8-fold corners (with the reference phases), a floor, the rows of D1's recipe
    sc = G.Scorer(raw, 'D2test')
    raw['test|D2|n_frames'] = np.array(200)
    raw['test|D2|ref_fit'] = G.jdump(dict(n_exp=3000.0, fit=dict(N=1.0, mu=R2, w=10.0, b=1e-5, floor=1.0, n=1.0, ok=True, x=[0.0, 0.0, R2, 2.3])))
    spec = {'R0': (1800, 18.0), 'R1E': (2000, 16.0), 'R2E': (2200, 14.0), 'ISS': (2600, 9.0), 'ISS0': (2400, 10.0), 'GREEDY': (2500, 9.5), 'PUB': (2300, 11.0)}
    for row, (n_tot, sd) in spec.items():
        pts = []
        for c, ph in zip(sc.allc[sc.scored], sc.phase_all[sc.scored]):
            k = max(int(n_tot / sc.n_pores), 1)
            th = ph + rng.integers(0, 8, k) * np.pi / 4 + 0.1 * rng.standard_normal(k)
            r = R2 + sd * rng.standard_normal(k)
            pts.append(np.c_[c[0] + r * np.cos(th), c[1] + r * np.sin(th)])
            pts.append(c + rng.uniform(-150, 150, (3, 2)))
        pts = np.concatenate(pts)
        G.put_loc(raw, 'test|D2|%s' % row, G.Loc(rng.integers(0, 200, len(pts)), pts), secs=np.full(200, 0.4) if row != 'PUB' else None)
        if row in ('ISS', 'ISS0', 'GREEDY'):
            raw['test|D2|%s|entries' % row] = np.full(200, 12.0)
    return raw


def test_report_runs_on_a_synthetic_raw_is_deterministic_and_states_the_rule(tmp_path):
    raw = synth_full_raw()
    G.save_raw(raw, str(tmp_path / 'raw.npz'))
    md1 = G.build_md(G.load_raw(str(tmp_path / 'raw.npz')), n_boot=12)
    md2 = G.build_md(G.load_raw(str(tmp_path / 'raw.npz')), n_boot=12)
    assert md1 == md2
    for needle in ('Outcome (registered rule', 'Validity gates', 'Corner contrast C8', 'Development gate', 'It is **not** a measurement of detection accuracy', 'ThunderSTORM was not run',
                   'ISS tuned - greedy tuned (Amendment 1, descriptive)', '12 resamples'):
        assert needle in md1, needle
    assert 'ISS better on BOTH' in md1 or 'MIXED' in md1 or 'DROP' in md1 or 'NO DETECTABLE' in md1


def test_block_shift_recovers_a_planted_translation_and_register_loc_removes_it_per_block():
    rng = np.random.default_rng(10)
    base = rng.uniform(1000, 9000, (4000, 2))
    dx, dy = -234.0, 97.0
    sh = G.block_shift(base + [dx, dy], base, 12000.0)
    assert abs(sh[0] - (-dx)) < 6 and abs(sh[1] - (-dy)) < 6 and sh[2] > 10
    fr = rng.integers(0, 500, 4000)
    ref = G.Loc(fr, base)
    moved = G.Loc(fr, base + np.where(fr[:, None] < 250, [dx, dy], [dx + 40.0, dy]))        # the shift changes between the two blocks
    out, shifts = G.register_loc(moved, ref, 500, 12000.0)
    assert np.abs(out.xy - base).mean() < 10.0 and shifts.shape == (2, 3)
    assert abs(shifts[1, 0] - (-(dx + 40.0))) < 6
    sparse = G.Loc(fr[:10], base[:10] + 100.0)                                              # fewer than 50 localizations: the all-frame shift is used
    assert G.register_loc(sparse, ref, 500, 12000.0)[1].shape == (2, 3)


def test_analyse_with_registration_removes_a_planted_offset_of_the_rows():
    raw = synth_test_raw()
    sc = G.Scorer(raw, 'D1test')
    for row in ('R0', 'ISS'):
        loc = G.get_loc(raw, 'test|D1|%s' % row)
        G.put_loc(raw, 'test|D1|%s' % row, G.Loc(loc.frame, loc.xy + [230.0, 100.0]), secs=np.full(200, 0.5))
    raw['test|D1|PUB|frame'] = np.concatenate([G.get_loc(raw, 'test|D1|%s' % r).frame for r in ('R2E',)])
    raw['test|D1|PUB|xy'] = G.get_loc(raw, 'test|D1|R2E').xy
    plain = G.analyse(raw, 'D1', n_boot=5)
    reg = G.analyse(raw, 'D1', n_boot=5, register=True)
    assert plain['rows']['ISS']['fit']['w'] > reg['rows']['ISS']['fit']['w'] + 3.0
    assert reg['rows']['ISS']['fit']['N'] > 2 * plain['rows']['ISS']['fit']['N'] or plain['rows']['ISS']['fit']['N'] < 50
    assert abs(reg['shifts']['ISS'][0][0] + 230.0) < 15


def test_find_pores_rejects_a_candidate_whose_radius_is_far_from_the_median_and_keeps_pores_closer_than_three_start_radii():
    rng = np.random.default_rng(11)
    centres = [(1000.0, 1000.0), (2400.0, 1000.0), (1000.0, 2400.0), (2400.0, 2400.0), (1700.0, 1700.0)]
    pts = [ring_points(rng, np.array(c), 500, sd=6.0) for c in centres]
    odd = np.array([3300.0, 1700.0])
    pts.append(ring_points(rng, odd, 800, R=112.0, sd=4.0))                              # radius 112 nm: +29 % of the median, outside +-25 %
    pts = np.concatenate(pts)
    cen, rad, cnt = G.find_pores(pts[:, 0], pts[:, 1], (0.0, 0.0, 4000.0, 3000.0), G.DS['D1'])
    assert len(cen) == 5 and np.all(np.hypot(*(cen - odd).T) > 300.0)
    # two pores 190 nm apart (more than 1.5 x R_start = 120 nm, less than 3 x R_start): both are kept
    pts2 = np.r_[ring_points(rng, np.array([1500.0, 1500.0]), 600, sd=5.0), ring_points(rng, np.array([1690.0, 1500.0]), 600, sd=5.0)]
    cen2, _, _ = G.find_pores(pts2[:, 0], pts2[:, 1], (0.0, 0.0, 3000.0, 3000.0), G.DS['D1'])
    assert len(cen2) == 2


def test_pore_phase_is_reported_in_zero_to_a_quarter_pi_and_uses_only_the_ring_window():
    rng = np.random.default_rng(12)
    sp = G.DS['D2']
    phi = -0.3
    th = phi + rng.choice(np.arange(8) * np.pi / 4, 500) + 0.06 * rng.standard_normal(500)
    ring = np.c_[53.0 * np.cos(th), 53.0 * np.sin(th)]
    other = rng.uniform(0, 2 * np.pi, 500)
    contam = np.c_[13.0 * np.cos(other), 13.0 * np.sin(other)]                           # 40 nm inside the ring: outside the ring window (25 nm), inside twice it
    ph = G.pore_phases(np.r_[ring[:, 0], contam[:, 0]], np.r_[ring[:, 1], contam[:, 1]], np.array([[0.0, 0.0]]), 53.0, sp)
    assert abs(ph[0, 0] - (phi + np.pi / 4)) < 0.03 and ph[0, 3] == 500 and 0 <= ph[0, 0] < np.pi / 4


def test_post_width_keeps_a_width_equal_to_the_bound():
    d = G.post_width([np.array([[0, 0, 1, 2.0], [0, 0, 1, 2.0000000001], [0, 0, 1, 2.01]])], 2.0)[0]
    assert d.shape[0] == 2


def test_equal_yield_thinning_thins_the_larger_set_to_the_smaller_yield_n_thin_times(monkeypatch):
    raw, _ = synth_raw()
    sc = G.Scorer(raw, 'D1test')
    Ha = np.full(sc.nb, 400.0)
    Hb = np.full(sc.nb, 100.0)
    fa, fb = dict(N=4000.0, w=10.0, x=[8.0, -9.0, R0, 2.3]), dict(N=1000.0, w=20.0, x=[7.0, -9.0, R0, 3.0])
    seen = []
    orig = G.fit_profile

    def spy(counts, area, edges, R, w0=20.0, start=None):
        seen.append(float(np.sum(counts)))
        return orig(counts, area, edges, R, w0, start)
    monkeypatch.setattr(G, 'fit_profile', spy)
    G.eq_yield_w(sc, fa, Ha, fb, Hb, sc.area, np.random.default_rng(1))
    assert len(seen) == G.N_THIN
    assert abs(np.mean(seen) / Ha.sum() - 0.25) < 0.02                                    # keep probability N_small / N_large


def test_percentile_ci_is_the_central_95_percent_interval():
    ci = G.percentile_ci(np.arange(1001.0))
    assert np.allclose(ci, [25.0, 975.0])
    assert np.allclose(G.percentile_ci([1.0, np.nan, 3.0]), [1.05, 2.95])


def test_prep_row_assigns_frames_to_blocks_of_ten():
    raw, cen = synth_raw()
    sc = G.Scorer(raw, 'D1test')
    c = sc.allc[sc.scored][0]
    loc = G.Loc(np.array([0, 9, 10, 25]), np.tile(c + [R0, 0.0], (4, 1)))
    p = G.prep_row(sc, loc, 3, False)
    assert list(p['H'].sum(axis=1)) == [2.0, 1.0, 1.0] and p['n_keep'] == 4


def test_table_loc_gives_zero_based_frames_inside_the_set_and_adds_the_convention(monkeypatch):
    tab = dict(frame=np.array([10000, 10001, 10002, 12000, 12001]), x=np.arange(5) * 1.0, y=np.arange(5) * 2.0)
    monkeypatch.setattr(G, 'load_table', lambda ds, field=None: dict(tab))
    loc = G.table_loc('D1', 'D1test', np.array([5.0, -5.0]))
    assert list(loc.frame) == [0, 1, 1999] and np.allclose(loc.xy[0], [6.0, -3.0])


def test_verdict_uses_the_three_point_margin_and_a_strictly_negative_upper_limit():
    def pr(dE, dw):
        return dict(dE_ci=dE, dw_ci=dw, dE_ci_pores=dE, dw_ci_pores=dw)
    assert G.verdict_ds(pr([-4.0, 4.0], [-1.0, 1.0])).startswith('no detectable')            # wider than +-3
    assert G.verdict_ds(pr([-2.9, 2.9], [-1.0, 1.0])).startswith('equivalent')
    assert G.verdict_ds(pr([-6.0, 0.5], [-1.0, 1.0])).startswith('no detectable')            # the upper limit is above 0: Picasso is not better


def test_the_iss_calls_the_guarded_path_with_the_requested_order_and_every_kappa_block(monkeypatch):
    from src.core import imaging as IM
    from src.core import imaging_guard as IG
    orders, kappas = [], []
    orig_path, orig_stop = IG.sliding_path_guarded, IM.stop_snr

    def path(*a, **kw):
        orders.append((kw.get('order'), kw.get('cycle_guard')))
        return orig_path(*a, **kw)

    def stop(hist, taus, snr, kappa=4.5):
        kappas.append(kappa)
        return orig_stop(hist, taus, snr, kappa)
    monkeypatch.setattr(IG, 'sliding_path_guarded', path)
    monkeypatch.setattr(IM, 'stop_snr', stop)
    rng = np.random.default_rng(13)
    frames = np.stack([synth_frame(rng, np.array([[20.0, 20.0]]))])
    cfg = G.grid_cfg('hp5', 1.5)
    cal = dict(bg=np.full((40, 40), 300.0), a=0.8, b0=100.0, x_ref=200.0)
    for order in ('dual', 'greedy'):
        G.iss_run(frames, cal, cfg, order, 100.0, log=lambda s: None)
    assert orders == [('dual', True), ('greedy', True)]
    assert sorted(set(kappas)) == sorted(G.KAPPAS) and len(kappas) == 2 * len(G.KAPPAS)
