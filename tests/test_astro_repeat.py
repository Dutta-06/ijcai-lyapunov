"""Tests of experiments/astro_repeat.py (the four M2 test windows in other SDSS runs), data/astro/fetch_repeat.py and experiments/exp_astro_baselines_repeat.py: the window mapping and the
integer origin, the margins and the five rules on a synthetic frame, the selection rule and its tie-breaking, the loaders on a tiny fake patch directory (also against the loaders of
exp_astro_m2), the registry, the fetch script's items and hashes, the baselines' frozen settings, and (skipped when the data are absent) the survey, the patches and the truth lists that
results/astro/repeat_spec.md reports. Data tests read $ASTRO_DATA_3 (fetch_repeat.py, then `astro_repeat.py build`) and need astropy for the frame."""
import importlib.util
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
import astro_common as A      # noqa: E402
import astro_repeat as R      # noqa: E402
import astro_second as S      # noqa: E402


def _astrans(**kw):
    t = dict(NODE=48.377, INCL=0.00036, A=323.2609, B=1.1e-4, C=2.9e-9, D=-0.8463, E=1.1e-9, F=1.1e-4, DROW0=0.03, DROW1=-9e-5, DROW2=7e-8, DROW3=-1.3e-11, DCOL0=-0.16, DCOL1=7.8e-4,
             DCOL2=-8.5e-7, DCOL3=2.6e-10, CSROW=0.017, CSCOL=0.025, CCROW=1.7, CCCOL=2.5, RICUT=100.0)
    t.update(kw)
    return t


# ----------------------------------------------------------------------------- the windows
def test_window_centre_sky_position_and_integer_origin_round_trip():
    t = _astrans()
    ra, dec = R.window_sky(t)
    assert ra.shape == dec.shape == (4,)
    # the centre of the box with origin (x0, y0) is the pixel (x0 + 49.5, y0 + 49.5); mapped to the sky and back it gives the same box
    for (xc, yc, x0, y0), w in zip(R.window_pixels(ra, dec, t), R.M2_WINDOWS):
        assert abs(xc - (w[0] + 49.5)) < 1e-7 and abs(yc - (w[1] + 49.5)) < 1e-7 and (x0, y0) == w
    assert R.window_centre(310, 630) == (359.5, 679.5) and R.window_centre(0, 0, 51) == (25.0, 25.0)
    # the same sky positions in a frame whose asTrans solution is shifted by whole pixels land on boxes shifted by the same whole pixels; the SDSS rows run along mu, the columns along nu
    z = dict(DROW0=0, DROW1=0, DROW2=0, DROW3=0, DCOL0=0, DCOL1=0, DCOL2=0, DCOL3=0)           # without the column distortion a whole-pixel shift of the constants is a whole-pixel shift of the frame
    t = _astrans(**z)
    ra, dec = R.window_sky(t)
    s = _astrans(A=t['A'] - 20.0 * t['B'], D=t['D'] + 13.0 * t['F'], **z)                       # the frame's row 0 starts 20 rows later and its column 0 13 columns earlier
    got = R.window_pixels(ra, dec, s)
    for (xc, yc, x0, y0), w in zip(got, R.M2_WINDOWS):
        assert abs(xc - (w[0] + 49.5 - 13.0)) < 1e-3 and abs(yc - (w[1] + 49.5 + 20.0)) < 1e-3 and (x0, y0) == (w[0] - 13, w[1] + 20)


def test_origin_of_centre_rounds_half_up_and_centres_the_box_to_half_a_pixel():
    for xc, want in ((59.5, 10), (59.49, 10), (59.99, 10), (60.0, 11), (58.99, 9), (49.5, 0), (49.0, 0), (48.99, -1), (-10.2, -60), (2047.9, 1998)):
        assert R.origin_of_centre(xc, xc)[0] == want, xc
    rng = np.random.default_rng(1)
    for xc, yc in rng.uniform(-100, 2100, (200, 2)):
        x0, y0 = R.origin_of_centre(xc, yc)
        assert abs(x0 + 49.5 - xc) <= 0.5 + 1e-12 and abs(y0 + 49.5 - yc) <= 0.5 + 1e-12
    assert R.origin_of_centre(59.5, 59.5) == (10, 10) and R.origin_of_centre(59.5 + 0.5, 59.5 - 0.5 - 1e-9) == (11, 9)        # half-way: up; just below half-way: down


def test_edge_margin_and_the_forty_pixel_rule():
    assert R.edge_margin(310, 630) == 310 and R.edge_margin(100, 1100) == 100 and R.edge_margin(1948, 600) == 0 and R.edge_margin(1950, 600) == -2
    assert R.edge_margin(500, 1389) == 0 and R.edge_margin(500, 1390) == -1 and R.edge_margin(-5, 500) == -5
    ctx = R.frame_context(dict(px=np.array([]), py=np.array([]), m606=np.array([])), dict(px=np.array([]), py=np.array([]), f814=np.array([])), np.zeros((1489, 2048)), 50000.0)
    for x0, y0, ok in ((40, 40, True), (39, 40, False), (40, 39, False), (1908, 40, True), (1909, 40, False), (40, 1349, True), (40, 1350, False)):
        r = R.window_rules(ctx, x0, y0)
        assert r['margin'] == ok and r['margin_px'] == R.edge_margin(x0, y0) and (r['margin_px'] >= S.FRAME_MARGIN) == ok


def _grid_truth(lo, hi, step=7.0, mag=20.0):
    gx, gy = np.meshgrid(np.arange(lo, hi, step), np.arange(lo, hi, step))
    return dict(px=gx.ravel(), py=gy.ravel(), m606=np.full(gx.size, mag), e606=np.zeros(gx.size), nf606=np.ones(gx.size))


def test_window_rules_saturation_bright_star_lost_star_and_footprint_on_a_synthetic_frame():
    shape = (400, 500)
    truth = _grid_truth(0.0, 400.0)                                                           # truth sources everywhere in x, y < 400: the footprint ends at x = 400
    none = dict(px=np.array([]), py=np.array([]), f814=np.array([]))
    x0, y0, sat = 150, 150, 50000.0
    dn = np.full(shape, 300.0)
    base = R.window_rules(R.frame_context(truth, none, dn, sat, shape), x0, y0)
    assert base['eligible'] and base['eligible_literal'] and base['footprint'] and base['saturation'] and base['bright'] and base['lost'] and base['dn_max'] == 300.0
    # a pixel at 0.8 x the saturation level 10 px from the box fails rule (3), 11 px away it is not looked at; exactly 0.8 x level fails (the rule is "below")
    for px_, ok in ((x0 - 10, False), (x0 - 11, True), (x0 + 99 + 10, False), (x0 + 99 + 11, True)):
        d2 = dn.copy()
        d2[y0 + 50, px_] = 0.8 * sat
        r = R.window_rules(R.frame_context(truth, none, d2, sat, shape), x0, y0)
        assert r['saturation'] == ok and r['eligible'] == ok and r['eligible_literal'] == ok
    d2 = dn.copy()
    d2[y0 + 50, x0 + 50] = 0.8 * sat - 1.0
    assert R.window_rules(R.frame_context(truth, none, d2, sat, shape), x0, y0)['saturation']
    # a truth source brighter than 14.5 within 10 px of the box fails rule (4) only: the window is still ELIGIBLE (the refined rule R3') but not eligible_literal
    br = dict(truth, px=np.r_[truth['px'], x0 - 10.4], py=np.r_[truth['py'], y0 + 30.0], m606=np.r_[truth['m606'], 14.0])
    r = R.window_rules(R.frame_context(br, none, dn, sat, shape), x0, y0)
    assert not r['bright'] and r['eligible'] and not r['eligible_literal'] and r['saturation']
    far = dict(truth, px=np.r_[truth['px'], x0 - 10.6], py=np.r_[truth['py'], y0 + 30.0], m606=np.r_[truth['m606'], 14.0])
    assert R.window_rules(R.frame_context(far, none, dn, sat, shape), x0, y0)['bright']
    # a bright source without F606W (F814W < 17) near the box fails rule (5); a fainter one does not
    r = R.window_rules(R.frame_context(truth, dict(px=np.array([x0 + 20.0]), py=np.array([y0 - 8.0]), f814=np.array([16.5])), dn, sat, shape), x0, y0)
    assert not r['lost'] and r['bright'] and r['eligible'] and not r['eligible_literal']
    assert R.window_rules(R.frame_context(truth, dict(px=np.array([x0 + 20.0]), py=np.array([y0 - 8.0]), f814=np.array([17.5])), dn, sat, shape), x0, y0)['lost']
    # the footprint: a box whose grown cells reach beyond the last truth source fails rule (2) only
    r = R.window_rules(R.frame_context(truth, none, dn, sat, shape), 300, 150)
    assert not r['footprint'] and r['eligible'] and not r['eligible_literal']
    # a window whose DN neighbourhood is empty has no maximum and is not eligible
    r = R.window_rules(R.frame_context(truth, none, dn, sat, shape), 600, 150)
    assert np.isnan(r['dn_max']) and not r['saturation'] and not r['eligible'] and not r['margin']


def test_window_rules_agree_with_astro_second_screen_patches_on_the_grid():
    """On the 10 px grid the five rules for an origin are those of `astro_second.screen_patches` for the same origin."""
    rng = np.random.default_rng(5)
    shape = (360, 380)
    n = 9000
    truth = dict(px=rng.uniform(0, 380, n), py=rng.uniform(0, 360, n), m606=rng.uniform(15.0, 25.0, n), e606=np.zeros(n), nf606=np.ones(n))
    truth['m606'][:3] = (13.0, 14.4, 14.6)
    truth['px'][:3], truth['py'][:3] = (150.0, 250.0, 330.0), (100.0, 200.0, 150.0)
    dn = rng.uniform(100, 300, shape)
    dn[200, 150], dn[120, 300] = 45000.0, 43000.0
    missing = dict(px=np.array([200.0, 120.0]), py=np.array([100.0, 250.0]), f814=np.array([16.0, 18.0]))
    ctx = R.frame_context(truth, missing, dn, 52900.0, shape)
    screened = S.screen_patches(truth, dn, 52900.0, shape, missing=missing)
    assert len(screened) > 10
    for q in screened:
        r = R.window_rules(ctx, q['x0'], q['y0'])
        assert (r['margin'], r['footprint'], r['saturation'], r['bright'], r['lost']) == (q['margin'], q['footprint'], q['saturation'], q['bright'], q['lost'])
        assert abs(r['dn_max'] - q['dn_max']) < 1e-12


# ----------------------------------------------------------------------------- the selection
def _f(run, camcol, field, fwhm, n_w):
    return dict(run=run, camcol=camcol, field=field, fwhm=fwhm, n_w=n_w)


def test_quantised_distance_and_the_sort_key():
    assert R.quantised_distance(2.70, 2.7) == 0 and R.quantised_distance(2.73, 2.7) == 1 and R.quantised_distance(2.67, 2.7) == 1 and R.quantised_distance(2.78, 2.7) == 2
    assert R.quantised_distance(3.52, 3.5) == 0 and R.quantised_distance(2.279, 2.7) == 8 and R.quantised_distance(3.047, 3.1) == 1
    a, b = _f(5, 2, 10, 3.51, 2), _f(5, 2, 11, 3.49, 3)
    assert R.select_key(a, 3.5) == (0, -2, 5, 2, 10) and R.select_key(b, 3.5) == (0, -3, 5, 2, 11) and R.select_key(b, 3.5) < R.select_key(a, 3.5)         # equal q: more windows first


def test_selection_ties_are_broken_by_windows_then_run_then_camcol_then_field():
    T = (3.5,)
    # same quantised distance (both within 0.025 px of 3.5): the frame with more eligible windows wins, whatever the exact distance
    assert R.select_frames([_f(100, 2, 10, 3.50, 2), _f(200, 2, 10, 3.52, 3)], T)[0]['frame']['run'] == 200
    # equal windows: the smaller run; then the smaller camcol; then the smaller field
    assert R.select_frames([_f(300, 2, 10, 3.50, 2), _f(200, 2, 10, 3.51, 2)], T)[0]['frame']['run'] == 200
    assert R.select_frames([_f(200, 3, 10, 3.50, 2), _f(200, 2, 10, 3.51, 2)], T)[0]['frame']['camcol'] == 2
    assert R.select_frames([_f(200, 2, 12, 3.50, 2), _f(200, 2, 11, 3.51, 2)], T)[0]['frame']['field'] == 11
    # a clearly closer frame beats one with more windows
    assert R.select_frames([_f(100, 2, 10, 3.50, 2), _f(50, 2, 10, 3.65, 4)], T)[0]['frame']['run'] == 100
    # the order of the list does not matter
    fr = [_f(300, 2, 10, 3.50, 2), _f(200, 2, 10, 3.51, 2), _f(100, 2, 10, 3.52, 2), _f(100, 1, 10, 3.52, 2)]
    assert R.select_frames(fr, T) == R.select_frames(fr[::-1], T) and R.select_frames(fr, T)[0]['frame'] == _f(100, 1, 10, 3.52, 2)


def test_eligible_frames_need_two_windows_a_finite_seeing_and_not_run_2583():
    fr = [_f(100, 2, 10, 3.0, 1), _f(101, 2, 10, 3.0, 2), _f(2583, 2, 136, 2.1, 4), _f(102, 2, 10, float('nan'), 2), _f(103, 2, 10, 3.0, 4)]
    assert [f['run'] for f in R.eligible_frames(fr)] == [101, 103]
    assert [f['run'] for f in R.eligible_frames(fr, min_windows=1)] == [100, 101, 103]
    assert R.select_frames([_f(100, 2, 10, 3.0, 1)]) == [] and R.select_frames([]) == []


def test_closest_pair_first_gives_the_target_to_the_frame_that_is_closer_and_never_repeats_a_frame():
    A_, B_, C_ = _f(1, 2, 1, 3.05, 2), _f(2, 2, 1, 2.28, 2), _f(3, 2, 1, 3.52, 2)
    got = R.select_frames([A_, B_, C_])
    assert [(c['target'], c['frame']['run']) for c in got] == [(3.5, 3), (3.1, 1), (2.7, 2)]               # matched in the order of closeness; the 3.05 frame goes to 3.1, not to 2.7
    assert [c['q'] for c in got] == [0, 1, 8]
    old = R.select_frames_target_order([A_, B_, C_])                                                     # the first declared version: 2.7 takes the 3.05 frame (0.35 px away)
    assert [(c['target'], c['frame']['run']) for c in old] == [(2.7, 1), (3.1, 3), (3.5, 2)]               # and leaves the 2.28 frame for the 3.5 target
    for fr in ([A_, B_, C_], [A_, B_], [A_]):                                                            # fewer eligible frames than targets: what exists, no frame twice
        got = R.select_frames(fr)
        ids = [(c['frame']['run'], c['frame']['camcol'], c['frame']['field']) for c in got]
        assert len(ids) == len(fr) and len(set(ids)) == len(ids)
    # equal quantised distance of two targets for the same frame: the smaller target takes it
    assert R.select_frames([_f(1, 2, 1, 3.30, 2)], targets=(3.5, 3.1))[0]['target'] == 3.1


def test_frame_names_follow_the_seeing():
    chosen = [dict(target=3.5, q=0, frame=_f(3, 2, 75, 3.52, 2)), dict(target=3.1, q=1, frame=_f(7, 6, 72, 3.05, 2)), dict(target=2.7, q=8, frame=_f(4, 2, 105, 2.28, 2))]
    names = R.frame_names(chosen)
    assert list(names) == ['m2r1', 'm2r2', 'm2r3'] and [names[k]['run'] for k in names] == [4, 7, 3]


# ----------------------------------------------------------------------------- the registry and the loaders
def test_registry_is_consistent_and_has_distinct_origins():
    assert list(R.REPEAT) == ['m2r1', 'm2r2', 'm2r3']
    origins = []
    for name, e in R.REPEAT.items():
        assert set(e) >= {'title', 'root', 'prefix', 'truth', 'patches', 'windows', 'frame'}
        assert len(e['patches']) == len(e['windows']) >= R.MIN_WINDOWS and len(set(e['windows'])) == len(e['windows']) and set(e['windows']) <= set(R.M2_WINDOWS)
        run, cam, fld = e['frame']
        assert e['prefix'] == 'patch_r_%06d_%d_%04d' % (run, cam, fld) and e['truth'] == 'truth_hugs_ngc7089_%06d_%d_%04d_meth1.npz' % (run, cam, fld)
        assert run != R.M2_FRAME[0] and all(R.edge_margin(*p) >= S.FRAME_MARGIN for p in e['patches'])
        assert S.patch_name(run, cam, fld, *e['patches'][0]) == '%s_x%d_y%d_100x100' % (e['prefix'], *e['patches'][0])
        origins += e['patches']
    assert len(set(origins)) == len(origins)                                                             # the baselines' keys x<x0>_y<y0> are unique


def _fake_patch(root, run, camcol, field, x0, y0, n=6):
    d = os.path.join(root, S.patch_name(run, camcol, field, x0, y0, 100))
    os.makedirs(d)
    rng = np.random.default_rng(0)
    dn = rng.uniform(100, 300, (n, n))
    psf = np.abs(rng.normal(0, 1, (5, 5)))
    psf /= psf.sum()
    psf25 = psf[1:4, 1:4] / psf[1:4, 1:4].sum()                                                         # a different array from psf51, as in the real files
    np.savez(os.path.join(d, 'patch_arrays.npz'), dn=dn, electrons=dn * 4.6, nmgy=dn * 0.005, sky_dn=np.full((n, n), 150.0), psf51=psf, psf25=psf25, gain=4.6, dark_variance=1.0, x0=x0, y0=y0)
    json.dump(dict(x0=x0, y0=y0, size=n), open(os.path.join(d, 'meta.json'), 'w'))
    t = _grid_truth(0.0, 200.0, 5.0, 19.0)
    np.savez(os.path.join(root, 'truth_hugs_ngc7089_%06d_%d_%04d_meth1.npz' % (run, camcol, field)), px=t['px'], py=t['py'], m606=t['m606'], e606=t['e606'], nf606=t['nf606'])
    return dn, psf, t


def test_loaders_read_a_tiny_fake_patch_directory_like_exp_astro_m2(tmp_path, monkeypatch):
    dn, psf, t = _fake_patch(str(tmp_path), 4187, 2, 105, 16, 24)
    monkeypatch.setitem(R.REPEAT, 'fake', R.registry_entry('fake', 4187, 2, 105, [(16, 24)], [(310, 630)], root=str(tmp_path)))
    assert R.patch_path(16, 24, 'fake') == os.path.join(str(tmp_path), 'patch_r_004187_2_0105_x16_y24_100x100', 'patch_arrays.npz') and os.path.exists(R.patch_path(16, 24, 'fake'))
    p = R.load_patch(16, 24, 'fake')
    assert set(p) == {'dn', 'psf', 'gain', 'dark'} and np.array_equal(p['dn'], dn) and np.array_equal(p['psf'], psf) and p['gain'] == 4.6 and p['dark'] == 1.0
    assert p['dn'].dtype == float and isinstance(p['gain'], float) and isinstance(p['dark'], float) and p['psf'].shape == (5, 5)
    tr = R.truth_arrays('fake')
    assert set(tr) == {'px', 'py', 'm606'} and np.array_equal(tr['px'], t['px'])
    tp = A.truth_for_patch(tr, 16, 24)                                                                    # the truth list works with astro_common.truth_for_patch
    assert len(tp['mag']) == int(((t['px'] >= 15.5) & (t['px'] < 115.5) & (t['py'] >= 23.5) & (t['py'] < 123.5)).sum()) > 0 and tp['inner'].sum() < len(tp['mag'])
    try:
        import exp_astro_m2 as M2
    except Exception as e:                                                                                # pragma: no cover
        pytest.skip('exp_astro_m2 cannot be imported here (%s)' % e)
    monkeypatch.setitem(M2.SECOND, 'fake', dict(title='fake', root=str(tmp_path), prefix='patch_r_004187_2_0105', truth='truth_hugs_ngc7089_004187_2_0105_meth1.npz', patches=[(16, 24)]))
    q = M2.load_patch(16, 24, 'fake')
    assert set(q) == set(p) and all(np.array_equal(np.asarray(q[k]), np.asarray(p[k])) for k in p)         # the same dict as exp_astro_m2.load_patch
    qt = M2.truth_arrays('fake')
    assert set(qt) == set(tr) and all(np.array_equal(qt[k], tr[k]) for k in tr)
    assert M2.patch_path(16, 24, 'fake') == R.patch_path(16, 24, 'fake')


# ----------------------------------------------------------------------------- data/astro/fetch_repeat.py (no network)
def _fetch_repeat():
    spec = importlib.util.spec_from_file_location('fetch_repeat', os.path.join(HERE, '..', 'data', 'astro', 'fetch_repeat.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_fetch_repeat_items_urls_directories_and_groups():
    F = _fetch_repeat()
    items = F.build_items()
    assert items['sdss_frame_r_004187_2_0105'][2] == 'https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/4187/2/frame-r-004187-2-0105.fits.bz2'
    assert items['sdss_psfield_007778_6_0072'][2] == 'https://data.sdss.org/sas/dr17/eboss/photo/redux/301/7778/objcs/6/psField-007778-6-0072.fit'
    assert items['sdss_photofield_003388_2'][2] == 'https://data.sdss.org/sas/dr17/eboss/photoObj/301/3388/photoField-003388-2.fits'
    assert items['photorunall'][2] == 'https://data.sdss.org/sas/dr17/eboss/photoObj/photoRunAll-dr17.fits' and 'mast:HLSP/hugs/ngc7089/' in items['hugs_ngc7089'][2].replace('%3A', ':').replace('%2F', '/')
    assert len({os.path.dirname(F.dest_path(items, k)) for k in items}) == len(items)                      # every download has its own directory
    assert {v[0] for v in items.values()} == set(F.GROUPS)
    assert all(('sdss_photofield_%06d_%d' % (r, c)) in items for r in F.RUNS for c in range(1, 7))
    for e in R.REPEAT.values():                                                                           # the chosen frames are among the fetched ones, with their photoField
        run, cam, fld = e['frame']
        assert (run, cam, fld) in F.FRAMES and run in F.RUNS and 'sdss_photofield_%06d_%d' % (run, cam) in items
    assert F.build_items((), [(1, 2, 3)])['sdss_frame_r_000001_2_0003'][0] == 'frames'
    assert F.parse_frames(['4187:2:105']) == [(4187, 2, 105)]
    if os.path.exists(F.CHECKSUMS):
        assert set(json.load(open(F.CHECKSUMS))['files']) <= set(items)


def test_fetch_repeat_checks_a_first_download_against_the_earlier_record_and_refuses_a_corrupted_file(tmp_path, monkeypatch):
    F = _fetch_repeat()
    data = np.random.default_rng(0).bytes(5000)
    import hashlib
    digest = hashlib.sha256(data).hexdigest()
    second = tmp_path / 'second.json'
    second.write_text(json.dumps({'files': {'k': {'bytes': len(data), 'sha256': digest}}}))
    monkeypatch.setattr(F, 'CHECKSUMS', str(tmp_path / 'sums.json'))
    items = {'k': ('survey', 'd', 'http://127.0.0.1:1/file.bin', 'file.bin')}
    path = tmp_path / 'file.bin'
    path.write_bytes(data)
    sums = F.load_sums(F.CHECKSUMS)
    monkeypatch.setattr(F.F2, 'head_lastmod', lambda url: 'Wed, 01 Jan 2020 00:00:00 GMT')
    F.check_or_record(items, 'k', str(path), sums, 'Wed, 01 Jan 2020 00:00:00 GMT', second_path=str(second))
    rec = json.load(open(F.CHECKSUMS))['files']['k']
    assert rec['sha256'] == digest and rec['bytes'] == len(data) and rec['matches_earlier_record'] is True
    F.check_or_record(items, 'k', str(path), F.load_sums(F.CHECKSUMS), second_path=str(second))             # verified against its own record
    path.write_bytes(data[:100] + b'\x00\x01' + data[102:])                                                 # a corrupted copy: refused and moved aside
    with pytest.raises(SystemExit):
        F.check_or_record(items, 'k', str(path), F.load_sums(F.CHECKSUMS), second_path=str(second))
    assert os.path.exists(str(path) + '.bad') and not path.exists()
    # a first download that differs from the earlier record of the same key is refused too
    path2 = tmp_path / 'file2.bin'
    path2.write_bytes(data[::-1])
    with pytest.raises(SystemExit):
        F.check_or_record(items, 'k', str(path2), F.load_sums(str(tmp_path / 'other.json')), second_path=str(second))
    assert os.path.exists(str(path2) + '.bad')
    assert set(json.load(open(F.CHECKSUMS))['files']) == {'k'}                                                # the real checksums_repeat.json was never touched: only the temporary file was written


# ----------------------------------------------------------------------------- experiments/exp_astro_baselines_repeat.py (frozen settings; no photutils needed)
def test_baselines_repeat_frozen_settings_and_crowdsource_filter():
    import exp_astro_baselines_repeat as BR
    tag, cs = BR.frozen_choices()
    assert tag in __import__('exp_astro_baselines').EXT_TAGS and cs.startswith('crowdsource|') and BR.parse_cs_name(cs) == (5.0, 0.0, 0.0, np.inf)
    assert BR.parse_cs_name('crowdsource|snr3|ff0.8|qf0.5|sp0.02') == (3.0, 0.8, 0.5, 0.02)
    with pytest.raises(ValueError):
        BR.parse_cs_name('photutils|snr5|ff0|qf0|spinf')
    rng = np.random.default_rng(3)
    n = 50
    k = 'x1_y2|crowdsource|'
    b = {k + 'x': rng.uniform(0, 100, n), k + 'y': rng.uniform(0, 100, n), k + 'flux': rng.uniform(-10, 1000, n), k + 'dflux': rng.uniform(1, 100, n), k + 'qf': rng.uniform(0, 1, n),
         k + 'fracflux': rng.uniform(0, 1, n), k + 'spread': rng.normal(0, 0.05, n)}
    x, y, f = BR.crowdsource_filtered(b, 1, 2, BR.parse_cs_name(cs))
    m = (b[k + 'flux'] > 0) & (b[k + 'flux'] >= 5 * b[k + 'dflux'])
    assert np.array_equal(x, b[k + 'x'][m]) and np.array_equal(f, b[k + 'flux'][m]) and 0 < len(x) < n
    x2, _, _ = BR.crowdsource_filtered(b, 1, 2, (5.0, 0.5, 0.5, 0.02))
    assert len(x2) < len(x)
    for f_ in ((0.0, 0.8, 0.0, np.inf), (0.0, 0.0, 0.8, np.inf), (0.0, 0.0, 0.0, 0.02)):                   # each column of the post-filter acts on its own array
        mm = (b[k + 'flux'] > 0) & (b[k + 'fracflux'] >= f_[1]) & (b[k + 'qf'] >= f_[2]) & (np.abs(b[k + 'spread']) <= f_[3])
        assert np.array_equal(BR.crowdsource_filtered(b, 1, 2, f_)[0], b[k + 'x'][mm])
    try:
        import exp_astro_m2 as M2
    except Exception as e:                                                                                # pragma: no cover
        pytest.skip('exp_astro_m2 cannot be imported here (%s)' % e)
    for f_ in ((5.0, 0.0, 0.0, np.inf), (3.0, 0.5, 0.5, 0.02)):                                           # the same selection as exp_astro_m2.cs_dets
        a, c = BR.crowdsource_filtered(b, 1, 2, f_), M2.cs_dets(b, (1, 2), f_)
        assert all(np.array_equal(u, v) for u, v in zip(a, c))
    full = dict(b)
    for t_ in ('photutils', tag):
        for c_ in ('x', 'y', 'flux'):
            full['x1_y2|%s|%s' % (t_, c_)] = rng.uniform(0, 100, 7)
    cats = BR.frozen_catalogues(full, 1, 2)
    assert len(cats) == 3 and all(len(v) == 3 for v in cats.values()) and [len(v[0]) for v in cats.values()][:2] == [7, 7]


# ----------------------------------------------------------------------------- data tests (skipped when the data are absent)
def _have(*paths):
    return all(os.path.exists(p) for p in paths)


def test_survey_selection_reproduces_the_registry():
    path = os.path.join(R.ROOT3, 'repeat_survey.json')
    if not _have(path):
        pytest.skip('the survey JSON is not in $ASTRO_DATA_3 (python -I experiments/astro_repeat.py survey)')
    d = json.load(open(path))
    frames = [dict(f, fwhm=float('nan') if f['fwhm'] is None else f['fwhm']) for f in d['frames'] if not f.get('reference')]
    chosen = R.select_frames(frames)
    names = R.frame_names(chosen)
    assert list(names) == list(R.REPEAT)
    for name, f in names.items():
        e = R.REPEAT[name]
        assert (f['run'], f['camcol'], f['field']) == tuple(e['frame'])
        el = [q for q in f['windows'] if q['eligible']]
        assert [(q['x0'], q['y0']) for q in el] == list(e['patches']) and [tuple(q['window']) for q in el] == list(e['windows'])
    assert [f['fwhm'] for f in names.values()] == sorted(f['fwhm'] for f in names.values())
    assert sum(1 for f in frames if f['n_w_literal'] > 0) == 0                                         # the literal rule (all five) leaves nothing, section 2 of the spec
    ref = [f for f in d['frames'] if f.get('reference')][0]
    assert [tuple(q['window']) for q in ref['windows']] == list(R.M2_WINDOWS) and all((q['x0'], q['y0']) == tuple(q['window']) for q in ref['windows'])
    assert [q['saturation'] for q in ref['windows']] == [True, True, False, True] and not any(q['bright'] for q in ref['windows'])


# reported in results/astro/repeat_spec.md: patch -> (F606W < 22 in the box, scored, all valid F606W in the box)
REPORTED = {'m2r1': {(316, 124): (1065, 971, 5381), (316, 224): (2294, 2047, 9354)}, 'm2r2': {(133, 281): (1062, 969, 5361), (195, 689): (1011, 923, 5195)}, 'm2r3': {(370, 441): (1016, 915, 5193), (110, 501): (1630, 1450, 7121)}}


@pytest.mark.parametrize('name', sorted(REPORTED))
def test_repeat_patches_on_disk_match_the_spec(name):
    e = R.REPEAT[name]
    paths = [os.path.join(e['root'], e['truth'])] + [R.patch_path(x0, y0, name) for (x0, y0) in e['patches']]
    if not _have(*paths):
        pytest.skip('the patches of %s are not in $ASTRO_DATA_3 (python -I experiments/astro_repeat.py build)' % name)
    assert sorted(REPORTED[name]) == sorted(e['patches'])
    tr_file = np.load(paths[0])
    assert set(tr_file.files) == {'px', 'py', 'm606', 'e606', 'nf606'} and len(tr_file['px']) == 287653 and _have(paths[0].replace('_meth1.npz', '_meth1_header.npz'))
    tr = R.truth_arrays(name)
    old_keys = {'dn', 'electrons', 'nmgy', 'sky_dn', 'psf51', 'psf25', 'gain', 'dark_variance', 'x0', 'y0'}
    for (x0, y0), (n22, scored, n_all) in REPORTED[name].items():
        z = np.load(R.patch_path(x0, y0, name))
        assert set(z.files) == old_keys and z['dn'].shape == (100, 100) and z['psf51'].shape == (51, 51) and z['psf25'].shape == (25, 25) and int(z['x0']) == x0 and int(z['y0']) == y0
        assert abs(z['psf51'].sum() - 1) < 1e-12 and abs(z['psf25'].sum() - 1) < 1e-12 and np.allclose(z['electrons'], z['dn'] * float(z['gain']))
        p = R.load_patch(x0, y0, name)
        assert set(p) == {'dn', 'psf', 'gain', 'dark'} and np.array_equal(p['dn'], z['dn'])
        t = A.truth_for_patch(tr, x0, y0)
        assert int((t['mag'] < 22).sum()) == n22 and int(((t['mag'] < 22) & t['inner']).sum()) == scored and len(t['mag']) == n_all
        assert sorted(os.listdir(os.path.dirname(R.patch_path(x0, y0, name)))) == ['meta.json', 'patch_arrays.npz', 'patch_dn.fits']
        # the saturation condition of the registry: no pixel above 0.8 x the saturation level in the patch (the box + 10 px was checked by the survey)
        assert float(z['dn'].max()) < S.SAT_FRACTION * (52900.0 if name == 'm2r2' else 50631.0)


def test_m2_windows_map_to_themselves_in_the_m2_frame_and_the_chosen_frame_rebuilds_a_patch():
    pytest.importorskip('astropy')
    P = R.frame_paths(*R.M2_FRAME)
    if not _have(P['frame'], P['psfield']):
        pytest.skip('the M2 frame is not in $ASTRO_DATA_3 (python -I data/astro/fetch_repeat.py fetch --group m2)')
    fr = S.read_frame(P['frame'])
    ra, dec = R.window_sky(fr['astrans'])
    assert abs(ra[0] - 323.335696) < 2e-6 and abs(dec[0] + 0.807082) < 2e-6                              # the benchmark patch centre of data_access.md (323.3357, -0.8071)
    assert [(x0, y0) for (_, _, x0, y0) in R.window_pixels(ra, dec, fr['astrans'])] == list(R.M2_WINDOWS)
    e = R.REPEAT['m2r1']
    run, cam, fld = e['frame']
    P = R.frame_paths(run, cam, fld)
    if not _have(P['frame'], P['psfield'], R.patch_path(*e['patches'][0], 'm2r1')):
        pytest.skip('the frame of m2r1 or its patches are not on disk')
    f2, ps = S.read_frame(P['frame']), S.read_psfield(P['psfield'])
    x0, y0 = e['patches'][0]
    arrays, _ = S.extract_patch(f2, ps, x0, y0)
    z = np.load(R.patch_path(x0, y0, 'm2r1'))
    assert all(np.array_equal(np.asarray(arrays[k]), z[k]) for k in z.files)
    # the window of the repeat frame is centred on the sky position of the M2 window: the centre of the box maps to within 0.5 px of the pixel position of that sky position
    xc, yc = S.astrans_world2pix(ra[0], dec[0], f2['astrans'])
    assert abs(x0 + 49.5 - float(xc)) <= 0.5 and abs(y0 + 49.5 - float(yc)) <= 0.5
