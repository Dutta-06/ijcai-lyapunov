"""The astronomy controls (experiments/exp_astro_diag.py): the renderer's flux and position, the noise model, the injected truth, the selection rules, the stack of planted stars against the stamp, the splitting and
flux-bias statistics on planted catalogues, the run bookkeeping and the report on a planted raw. Everything is planted: no astronomy data is read except in the tests marked needs_data."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import astro_common as A                # noqa: E402
import astro_repeat as R                # noqa: E402
import exp_astro_diag as E              # noqa: E402
import exp_astro_m2 as X                # noqa: E402


def gauss_stamp(sx, sy=None, n=51):
    """A point-sampled Gaussian stamp of unit sum (widths sx along the columns, sy along the rows; the middle pixel is the centre) and the unnormalised sum."""
    sy = sx if sy is None else sy
    yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2.0
    g = np.exp(-0.5 * (xx / sx) ** 2 - 0.5 * (yy / sy) ** 2)
    return g / g.sum()


def planted_truth(rows, x0=0, y0=0):
    """A truth dict (frame pixel coordinates, as the truth files) from rows (x, y, mag) given in patch coordinates."""
    a = np.array(rows, float).reshape(-1, 3)
    return dict(px=a[:, 0] + x0, py=a[:, 1] + y0, m606=a[:, 2])


# ----------------------------------------------------------------------------------------------------------------------------------
# the renderer
# ----------------------------------------------------------------------------------------------------------------------------------
def test_a_planted_source_is_rendered_with_its_flux_at_its_position_and_with_the_stamps_shape():
    sx, sy = 1.6, 1.1                                                  # elliptical: a swap of the axes or of the sign of the shift moves the centroid by tenths of a pixel
    st = gauss_stamp(sx, sy)
    x, y, F = 40.3, 55.7, 12345.0
    img = E.render_sources(st, [x], [y], [F], (100, 100), margin=20)
    yy, xx = np.mgrid[0:100, 0:100]
    assert abs(img.sum() / F - 1) < 2e-3
    assert abs((img * xx).sum() / img.sum() - x) < 0.01 and abs((img * yy).sum() / img.sum() - y) < 0.01
    mxx, myy = (img * (xx - x) ** 2).sum() / img.sum(), (img * (yy - y) ** 2).sum() / img.sum()
    assert abs(np.sqrt(mxx) / sx - 1) < 0.01 and abs(np.sqrt(myy) / sy - 1) < 0.01
    g = np.exp(-0.5 * ((xx - x) / sx) ** 2 - 0.5 * ((yy - y) / sy) ** 2)
    ref = F * g / (2 * np.pi * sx * sy)                                # the point-sampled Gaussian of the stamp, at the sub-pixel position
    assert np.abs(img - ref).max() / ref.max() < 0.03
    # sources add, and an integer position gives the stamp itself
    half = E.render_sources(st, [30.0], [35.4], [F], (100, 100))                  # only the row is fractional
    assert abs((half * xx).sum() / half.sum() - 30.0) < 0.01 and abs((half * yy).sum() / half.sum() - 35.4) < 0.01
    two = E.render_sources(st, [x, 30.0], [y, 35.0], [F, 1000.0], (100, 100))
    one = E.render_sources(st, [30.0], [35.0], [1000.0], (100, 100))
    assert np.allclose(two, img + one) and np.allclose(one[35 - 25:35 + 26, 30 - 25:30 + 26], 1000.0 * st)


def test_sources_outside_the_image_contribute_their_wings_and_the_margin_selects_the_truth():
    st = gauss_stamp(4.0)                                              # wide, so that the wings reach 10 px
    img = E.render_sources(st, [-10.0], [50.0], [1.0e5], (100, 100), margin=20)
    assert img[:, :10].sum() > 0 and abs(img.sum() / 1.0e5 - st[:, 25 + 10:].sum()) < 1e-6        # exactly the part of the stamp that lies in the image
    truth = planted_truth([(-0.4, 10, 20.0), (-0.6, 10, 20.0), (99.4, 10, 20.0), (99.6, 10, 20.0), (-20.4, 50, 21.0), (-20.6, 50, 21.0), (50, 119.4, 22.0), (50, 119.6, 22.0), (50, 50, 30.0), (50, 51, 29.99)])
    x, y, m = E.sources_near(truth, 0, 0, 0)
    assert sorted(zip(x, y)) == [(-0.4, 10), (50, 51), (99.4, 10)] and sorted(m) == [20.0, 20.0, 29.99]
    x, y, m = E.sources_near(truth, 0, 0, 20)
    assert len(x) == 7 and -20.4 in x and -20.6 not in x and 119.4 in y and 119.6 not in y and 30.0 not in m and 29.99 in m           # the box grown by 20 px, F606W < 30
    x, y, m = E.sources_near(truth, 100, 200, 0)                                                                         # the patch origin is subtracted
    assert len(x) == 0
    x, y, m = E.sources_near(planted_truth([(5, 6, 20.0)], 100, 200), 100, 200, 0)
    assert (x[0], y[0], m[0]) == (5.0, 6.0, 20.0)


def test_a_stamp_that_does_not_vanish_at_its_edge_is_placed_whole_and_loses_its_edge_when_shifted():
    flat = np.full((5, 5), 1.0 / 25.0)                                 # the pixels at the stamp's edge matter here (a psField stamp is ~0 there, which hides the edge handling)
    F = 1000.0
    whole = E.render_sources(flat, [30.0], [40.0], [F], (80, 80), margin=20)
    assert np.allclose(whole[38:43, 28:33], F / 25.0) and abs(whole.sum() / F - 1.0) < 1e-12 and np.count_nonzero(whole) == 25      # all 25 pixels, at the pixels 28..32 and 38..42
    for x, y, kept in ((30.3, 40.0, 0.8), (30.7, 40.0, 0.8), (30.0, 40.4, 0.8), (30.0, 39.6, 0.8), (30.3, 40.4, 0.64)):
        img = E.render_sources(flat, [x], [y], [F], (80, 80), margin=20)
        assert abs(img.sum() / F - kept) < 1e-9, (x, y)                 # zero outside the stamp: a shift by a fraction of a pixel leaves one column (row) of five without a source pixel


def test_flux_follows_the_magnitude_law_the_scorer_inverts():
    st = gauss_stamp(1.5)
    truth = planted_truth([(50.0, 50.0, 19.0), (20.0, 70.0, 21.5)])
    w = E.render_world(truth, 0, 0, st, 4.6, 1.0, np.zeros((100, 100)), 3, noise=False)
    assert np.allclose(w['clean'], E.render_sources(st, w['x'], w['y'], A.mag_to_counts(w['mag']), (100, 100)))
    # the aperture flux of the brighter source is its counts (10^(-0.4 (m - 22.5)) / 0.00546689) and the scorer reads the magnitude back
    flux = w['clean'][50 - 12:50 + 13, 50 - 12:50 + 13].sum()
    assert abs(flux / A.mag_to_counts(19.0) - 1) < 1e-3 and abs(A.counts_to_mag(flux) - 19.0) < 2e-3
    assert abs(A.mag_to_counts(19.0) - 10 ** (-0.4 * (19.0 - 22.5)) / 0.00546689) < 1e-6 * A.mag_to_counts(19.0)


def test_the_sky_is_added_per_pixel_and_the_noise_follows_the_gain_and_the_dark_variance():
    sky = np.full((100, 100), 120.0)
    sky[:, 50:] = 130.0
    w = E.render_world(planted_truth([]), 0, 0, gauss_stamp(1.5), 4.6, 1.0, sky, 1, noise=False)
    assert np.array_equal(w['clean'], sky) and np.array_equal(w['dn'], sky)
    # statistics of the noise: a flat frame, variance DN / gain + dark variance (the dark term dominates the second frame)
    for level, gain, dark in ((400.0, 4.6, 1.0), (10.0, 4.6, 9.0), (50.0, 1.7, 0.25)):
        clean = np.full((400, 400), level)
        dn = E.add_noise(clean, gain, dark, np.random.default_rng(5))
        assert abs(dn.mean() - level) < 4 * np.sqrt((level / gain + dark) / dn.size) + 1e-9
        assert abs(dn.var() / (level / gain + dark) - 1) < 0.03, (level, gain, dark)
    neg = E.add_noise(np.full((10, 10), -3.0), 4.6, 1.0, np.random.default_rng(1))          # a negative model pixel is a pixel without signal, not an error
    assert np.isfinite(neg).all()
    assert abs(E.add_noise(np.full((300, 300), 0.0), 4.6, 4.0, np.random.default_rng(2)).var() - 4.0) < 0.2


def test_frame_statistics_are_the_median_the_two_tails_and_the_maximum():
    assert E.frame_stats(np.arange(101.0)).tolist() == [50.0, 5.0, 95.0, 100.0]
    assert np.allclose(E.frame_stats(np.array([[1.0, 2.0], [3.0, 10.0]])), [2.5, 1.15, 8.95, 10.0], atol=1e-12)


def test_the_world_renders_the_sources_up_to_twenty_pixels_outside_the_box():
    st = gauss_stamp(4.0)
    sky = np.full((100, 100), 100.0)
    near = E.render_world(planted_truth([(-15.0, 50.0, 17.0)]), 0, 0, st, 4.6, 1.0, sky, 5, noise=False)['clean']
    far = E.render_world(planted_truth([(-25.0, 50.0, 17.0)]), 0, 0, st, 4.6, 1.0, sky, 5, noise=False)['clean']
    assert E.MARGIN == 20 and E.MAG_MAX == 30.0
    assert near[:, :10].sum() > 100.0 * 100 * 10 + 1.0 and np.array_equal(far, sky)                   # a source 15 px outside leaves its wings in the image, one 25 px outside is not rendered


def test_the_header_of_a_world_describes_the_stamp_the_frames_and_the_sources(monkeypatch):
    truth = planted_truth([(30.3, 30.7, 19.0), (60.2, 61.4, 20.5), (75.5, 25.5, 21.9), (-10.0, 20.0, 18.0), (45.0, 80.0, 22.1), (1.0, 50.0, 19.0)])
    patch = dict(dn=np.full((100, 100), 150.0), psf=gauss_stamp(1.5), gain=4.6, dark=1.0)
    monkeypatch.setattr(E, 'load_inputs', lambda name, pxy: (patch, np.full((100, 100), 150.0), truth, 5.5))
    h = E.world_header('m92', (0, 0))
    assert str(h['kind'][0]) == 'isotropic' and h['stats'].shape == (2, 4) and h['stats'][0].tolist() == [150.0, 150.0, 150.0, 150.0] and 140 < h['stats'][1][0] < 160 and h['stats'][1][3] > 300
    assert h['counts'].tolist() == [6, 5, 3] and abs(h['fwhm'][0] / (2.3548 * 1.5) - 1) < 0.01 and 'rms' in h and 'parts5' in h


def test_a_world_is_reproducible_and_its_noise_differs_by_patch_and_set():
    truth = planted_truth([(30.5, 40.5, 20.0), (70.2, 20.7, 21.0)])
    st = gauss_stamp(1.5)
    a = E.render_world(truth, 0, 0, st, 4.6, 1.0, np.full((100, 100), 100.0), E.world_seed('m92', (1100, 220)))
    b = E.render_world(truth, 0, 0, st, 4.6, 1.0, np.full((100, 100), 100.0), E.world_seed('m92', (1100, 220)))
    c = E.render_world(truth, 0, 0, st, 4.6, 1.0, np.full((100, 100), 100.0), E.world_seed('m92', (1430, 550)))
    assert np.array_equal(a['dn'], b['dn']) and not np.array_equal(a['dn'], c['dn']) and np.array_equal(a['clean'], c['clean'])
    seeds = {tuple(E.world_seed(n, p)) for n, p in E.SYNTH + [(n, p) for n in E.STAMP_SETS for p in E.set_patches(n)]}
    assert len(seeds) == len(set(E.SYNTH + [(n, p) for n in E.STAMP_SETS for p in E.set_patches(n)])) and E.world_seed('m2', (310, 630))[0] == 1108
    assert E.world_seed('m2', (310, 630)) != E.world_seed('m53', (310, 630)) and E.world_seed('m2', (310, 630)) != E.world_seed('m2', (310, 730))
    import zlib
    assert E.world_seed('m2', (310, 630)) == [1108, zlib.crc32(b'm2'), 310, 630] and E.world_seed('m92', (1430, 550)) == [1108, zlib.crc32(b'm92'), 1430, 550]         # the declared seed: [1108, crc32(set), x0, y0]


def test_the_injected_truth_of_a_world_is_the_truth_the_scorer_uses():
    rng = np.random.default_rng(4)
    rows = np.c_[rng.uniform(-30, 130, 400), rng.uniform(-30, 130, 400), rng.uniform(15, 31, 400)]
    truth = planted_truth(rows, 300, 700)
    for margin_free in ((300, 700), (310, 710)):
        x, y, m = E.sources_near(truth, margin_free[0], margin_free[1], 0)
        t = A.truth_for_patch(truth, margin_free[0], margin_free[1])
        assert np.array_equal(x, t['x']) and np.array_equal(y, t['y']) and np.array_equal(m, t['mag'])
    # the ten worlds are the patches whose real results the frozen raws hold, and the sets are those of the stamp part
    assert len(E.SYNTH) == 10 and [p for p in E.SYNTH if p[0] == 'm92'] == [('m92', (1100, 220)), ('m92', (1430, 550))]
    assert sum(len(E.set_patches(n)) for n in E.STAMP_SETS) == 18 and E.set_patches('m2') == [(310, 630), (310, 730), (360, 1040), (100, 1100), (110, 630), (210, 630)]
    assert all(tuple(p) in [tuple(q) for q in R.REPEAT[n]['patches']] for n, p in E.SYNTH if n in R.REPEAT)


# ----------------------------------------------------------------------------------------------------------------------------------
# item 1: the stack of isolated bright stars against the stamp
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_isolation_rule_of_the_stack_on_a_planted_truth_list():
    rows = [(50, 50, 19.0), (50, 53, 22.0),                          # A: a faint neighbour at 3 px (22 is fainter than 19 + 1.5): selected
            (30, 30, 19.0), (34, 30, 20.4),                          # B: a neighbour at 4 px brighter than 19 + 1.5: excluded
            (70, 30, 19.0), (74, 30, 20.6),                          # C: a neighbour at 4 px fainter than 19 + 1.5: selected
            (30, 70, 19.0), (34.6, 70, 15.0),                        # D: a bright neighbour at 4.6 px, beyond 4.5: selected
            (50, 80, 20.0), (60, 80, 19.99),                         # E: m = 20.0 is not below 20: excluded; F: 19.99: selected
            (12.0, 50, 19.0), (12.1, 40, 19.0), (86.9, 40, 19.0), (87.0, 50, 19.0),            # the edge: more than 12 px from the edge, i.e. 12 < x < 87
            (50, 12.0, 19.0), (60, 12.1, 19.0), (70, 86.9, 19.0), (80, 87.0, 19.0),
            (60, 60, 19.0), (60, 60, 19.5),                          # J: a second source at the same position, brighter than m + 1.5: excluded (a neighbour is not 'itself')
            (30, 50, 19.0), (34.3, 50, 15.0),                        # K: a bright neighbour at 4.3 px (inside 4.5): excluded; the bright one has K fainter than 15 + 1.5: selected
            (50, 20, 19.0), (54, 20, 20.5)]                          # L: a neighbour at exactly m + 1.5 is not brighter than m + 1.5: selected
    a = np.array(rows)
    idx = E.stack_stars(a[:, 0], a[:, 1], a[:, 2])
    chosen = sorted((a[i, 0], a[i, 1]) for i in idx)
    # D's bright neighbour (34.6, 70) is a star too (its own neighbour D is 4.6 px away); both members of the duplicate pair are excluded (each has the other at distance 0)
    assert chosen == sorted([(50, 50), (70, 30), (30, 70), (34.6, 70), (60, 80), (12.1, 40), (86.9, 40), (60, 12.1), (70, 86.9), (34.3, 50), (50, 20)]), chosen
    assert E.STACK_MAG == 20.0 and E.STACK_RAD == 4.5 and E.STACK_DM == 1.5 and E.HALF + E.PAD == 12
    # a different patch size moves the upper edge
    idx = E.stack_stars(np.array([40.0, 60.0]), np.array([30.0, 30.0]), np.array([19.0, 19.0]), n=60)
    assert list(idx) == [0]


def test_a_star_is_centred_on_its_flux_centroid_and_normalised_within_six_pixels():
    s, n = 1.3, 25
    yy, xx = np.mgrid[0:n, 0:n]
    for dx, dy in ((0.0, 0.0), (0.37, -0.29), (-0.45, 0.4)):
        sub = np.exp(-0.5 * ((xx - 12 - dx) ** 2 + (yy - 12 - dy) ** 2) / s ** 2)
        a = E.aligned(sub)
        assert a.shape == (19, 19)
        y9, x9 = np.mgrid[-9:10, -9:10]
        w = a * (np.hypot(y9, x9) <= 3)
        assert abs((w * x9).sum() / w.sum()) < 0.05 and abs((w * y9).sum() / w.sum()) < 0.05            # the centroid is the middle pixel (the iteration is cut at 6 steps)
        assert abs(a[np.hypot(y9, x9) <= 6].sum() - 1.0) < 1e-12                                          # unit sum within 6 px
        # the middle pixel against an independent sum (a Gaussian centred on the middle pixel)
        r2 = (y9 ** 2 + x9 ** 2)
        expect = 1.0 / np.exp(-0.5 * r2[r2 <= 36] / s ** 2).sum()
        assert abs(a[9, 9] / expect - 1) < 0.02, (dx, dy)
    c = E.centroid(np.exp(-0.5 * ((xx - 13.2) ** 2 + (yy - 11.1) ** 2) / s ** 2))
    assert abs(c[0] - 13.2) < 0.05 and abs(c[1] - 11.1) < 0.05
    flat = np.zeros((25, 25))
    assert E.centroid(flat) == (12.0, 12.0)                                                              # nothing to centre on: the middle pixel


def test_the_fwhm_of_a_centred_cutout_is_the_gaussian_one():
    for s in (1.0, 1.5, 2.2):
        n = 19
        yy, xx = np.mgrid[0:n, 0:n] - 9
        g = np.exp(-0.5 * (xx ** 2 + yy ** 2) / s ** 2)
        assert abs(E.fwhm(g / g.sum()) / (2.3548 * s) - 1) < 0.01, s
    assert np.isnan(E.fwhm(np.ones((19, 19))))                                                           # no half maximum: not a number


def test_the_centroid_ignores_negative_pixels_and_the_cutout_takes_its_sky_beyond_8_5_px():
    n = 25
    yy, xx = np.mgrid[0:n, 0:n]
    img = np.exp(-0.5 * ((xx - 12) ** 2 + (yy - 12) ** 2) / 1.5 ** 2)
    zeroed = img.copy()
    zeroed[12, 14] = 0.0
    img[12, 14] = -3.0                                                   # a negative pixel inside the window, 2 px to the right of the star
    assert E.centroid(img) == E.centroid(zeroed) and abs(E.centroid(img)[0] - 12) > 0.05                 # counted as zero (the star loses a pixel on one side: its centroid moves away from it)
    # the sky of a cut-out is the median of the pixels beyond 8.5 px: a broad star (sd 6 px) puts light in that ring, and the recipe is the same one written out here
    frame = 100.0 + 5.0e4 * np.exp(-0.5 * ((np.arange(100)[None, :] - 50.4) ** 2 + (np.arange(100)[:, None] - 50.7) ** 2) / 6.0 ** 2) / (2 * np.pi * 36.0)
    cut = E.cutouts(frame, np.array([50.4]), np.array([50.7]), [0])[0]
    sub = frame[51 - 12: 51 + 13, 50 - 12: 50 + 13].copy()                  # around the nearest pixel (50, 51) of the position (50.4, 50.7)
    yy, xx = np.mgrid[-12:13, -12:13]
    sub -= np.median(sub[np.hypot(yy, xx) > 8.5])
    assert np.allclose(cut, E.aligned(sub), rtol=0, atol=1e-12) and abs(np.median(frame[51 - 12: 51 + 13, 50 - 12: 50 + 13][np.hypot(yy, xx) > 8.5]) - 100.0) > 1.0
    # the cut-outs come in the order of the list of stars
    fr = np.full((100, 100), 10.0)
    yy, xx = np.mgrid[0:100, 0:100]
    for cx, cy, sd in ((30, 40, 1.0), (70, 60, 2.0)):
        fr += 1000.0 * np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / sd ** 2)
    two = E.cutouts(fr, np.array([30.0, 70.0]), np.array([40.0, 60.0]), [1, 0])
    assert two[0][9, 9] < 0.5 * two[1][9, 9]                              # the broad star (index 1) first: its middle pixel is the lower one (0.040 against 0.159)
    # the star is found at the nearest pixel (a position at x.6 is cut out around x + 1)
    frame2 = np.zeros((100, 100)) + 50.0
    frame2[40:61, 40:61] += 1000.0 * np.exp(-0.5 * ((np.arange(40, 61)[None, :] - 52.0) ** 2 + (np.arange(40, 61)[:, None] - 48.0) ** 2) / 1.2 ** 2)
    c2 = E.cutouts(frame2, np.array([51.6]), np.array([48.4]), [0])[0]
    assert abs(c2[9, 9] - E.cutouts(frame2, np.array([52.0]), np.array([48.0]), [0])[0][9, 9]) < 0.01


def test_the_paired_statistic_is_the_median_of_the_ratios_of_the_stars_central_pixels():
    ref = E.stamp_core(gauss_stamp(1.3))
    base = E.aligned(np.pad(gauss_stamp(1.3)[13:38, 13:38], 0))
    n = 5
    f = np.array([1.0, 1.1, 1.2, 1.3, 50.0])
    g = np.array([1.0, 1.0, 1.0, 1.0, 0.01])
    real, syn = np.repeat(base[None], n, 0), np.repeat(base[None], n, 0)
    real[:, 9, 9] = base[9, 9] * f
    syn[:, 9, 9] = base[9, 9] * g
    assert np.array_equal(ref, base)                                     # the stamp's core is the 25 x 25 pixels around its middle pixel (13..37), centred and normalised as a star is
    est = E.stack_estimates(real, syn, ref)
    assert abs(est[8] - 0.2) < 1e-12                                     # the median of f / g = 1.0, 1.1, 1.2, 1.3, 5000 is 1.2: the mean would be 1001
    assert abs(est[3] - 0.2) < 1e-12 and abs(est[0] - (1.2 * base[9, 9] / ref[9, 9] - 1.0)) < 1e-12 and abs(est[1] - (base[9, 9] / ref[9, 9] - 1.0)) < 1e-12
    # a patch without stars in the list of sizes is no stratum
    truth, w1 = planted_field(1.3, seed=31)
    _, w2 = planted_field(1.3, seed=32)
    idx = E.stack_stars(w1['x'], w1['y'], w1['mag'])
    r1, r2 = E.cutouts(w1['dn'], w1['x'], w1['y'], idx), E.cutouts(w2['dn'], w2['x'], w2['y'], idx)
    a = E.stack_bootstrap(r1, r2, ref, [8, 0, 8], B=10, seed=2)
    b = E.stack_bootstrap(r1, r2, ref, [8, 8], B=10, seed=2)
    assert all(np.array_equal(u, v) for u, v in zip(a[:3], b[:3]))


def planted_field(stamp_sigma, seed, n_side=4, flux=3.0e4, shift=0.0):
    """A 100 x 100 frame with n_side^2 stars on a grid (17 px apart, the sub-pixel positions fixed by the grid), rendered with a Gaussian stamp of width stamp_sigma; the truth is the same for every sigma."""
    xs = 20.0 + 17.0 * np.arange(n_side)
    rows = [(x + 0.31 * ((i + j) % 3) - 0.3 + shift, y + 0.27 * ((i * j) % 4) - 0.4, A.counts_to_mag(flux)) for i, x in enumerate(xs) for j, y in enumerate(xs)]
    truth = planted_truth(rows)
    st = gauss_stamp(stamp_sigma)
    w = E.render_world(truth, 0, 0, st, 4.6, 1.0, np.full((100, 100), 150.0), seed, margin=20)
    return truth, w


def test_the_stack_of_planted_stars_recovers_the_stamps_central_pixel_and_fwhm():
    truth, w1 = planted_field(1.3, seed=11)
    _, w2 = planted_field(1.3, seed=12)                                                                  # the same stars, another noise realisation: PSF = stamp in both
    idx = E.stack_stars(w1['x'], w1['y'], w1['mag'])
    assert len(idx) == 16
    ref = E.stamp_core(gauss_stamp(1.3))
    r1, r2 = E.cutouts(w1['dn'], w1['x'], w1['y'], idx), E.cutouts(w2['dn'], w2['x'], w2['y'], idx)
    est = E.stack_estimates(r1, r2, ref)
    assert est.shape == (9,) and abs(est[8]) < 0.02 and abs(est[0]) < 0.02 and abs(est[1]) < 0.02 and abs(est[2]) < 0.02 and abs(est[3]) < 0.02            # the stack of the stars is the stamp, centring cancels
    assert abs(est[4] - E.fwhm(ref)) < 0.05 and abs(E.fwhm(ref) / (2.3548 * 1.3) - 1) < 0.01 and abs(est[6]) < 0.05
    # stars 10 % broader than the stamp: the central pixel falls by the ratio of the areas, the FWHM grows by 10 %
    _, wb = planted_field(1.43, seed=13)
    rb = E.cutouts(wb['dn'], wb['x'], wb['y'], idx)
    est = E.stack_estimates(rb, r2, ref)
    assert -0.20 < est[3] < -0.14 and -0.20 < est[0] < -0.14 and abs(est[2] - est[0]) < 0.03 and -0.20 < est[8] < -0.14
    assert abs(est[6] - 0.306) < 0.08 and abs(est[7] - 0.10) < 0.03
    est = E.stack_estimates(r2, rb, ref)                                                                     # the roles exchanged: the synthetic stars are the broad ones
    assert abs(est[0]) < 0.02 and -0.20 < est[1] < -0.14 and 0.17 < est[3] < 0.25 and est[6] < -0.2 and 0.17 < est[8] < 0.25
    # more peaked stars: the other sign
    _, wn = planted_field(1.17, seed=14)
    est = E.stack_estimates(E.cutouts(wn['dn'], wn['x'], wn['y'], idx), r2, ref)
    assert 0.15 < est[3] < 0.30 and est[6] < -0.2 and 0.15 < est[8] < 0.30


def test_the_bootstrap_over_the_stars_brackets_the_estimate_and_is_paired_and_stratified():
    truth, w1 = planted_field(1.3, seed=21)
    _, w2 = planted_field(1.3, seed=22)
    idx = E.stack_stars(w1['x'], w1['y'], w1['mag'])
    ref = E.stamp_core(gauss_stamp(1.3))
    r1, r2 = E.cutouts(w1['dn'], w1['x'], w1['y'], idx), E.cutouts(w2['dn'], w2['x'], w2['y'], idx)
    est, ci, sd, nbad = E.stack_bootstrap(r1, r2, ref, [8, 8], B=300, seed=3)
    assert est.shape == (9,) and ci.shape == (9, 2) and sd.shape == (9,) and nbad == 0
    assert (ci[:, 0] <= est + 1e-9).all() and (est <= ci[:, 1] + 1e-9).all() and (sd > 0).all() and (ci[:, 1] > ci[:, 0]).all()
    assert abs(ci[2].mean()) < 0.05                                                                      # real minus synthetic: both stacks are the stamp
    rs = X.boot_idx_strat([8, 8], 300, 3)                                                                    # the intervals are the 2.5 and 97.5 percentiles and the sd the standard deviation of the resamples
    draws = np.array([E.stack_estimates(r1[i], r2[i], ref) for i in rs])
    assert np.allclose(ci, np.percentile(draws, [2.5, 97.5], axis=0).T) and np.allclose(sd, draws.std(0))
    again = E.stack_bootstrap(r1, r2, ref, [8, 8], B=300, seed=3)
    assert all(np.array_equal(u, v) for u, v in zip((est, ci, sd), again[:3]))
    other = E.stack_bootstrap(r1, r2, ref, [8, 8], B=300, seed=4)
    assert not np.array_equal(ci, other[1])
    # paired: if the 'synthetic' stack is the real stack, real minus synthetic is zero in every resample; stratified: with the stars of a patch all equal, resampling within the patch changes nothing
    est, ci, sd, nbad = E.stack_bootstrap(r1, r1, ref, [16], B=100, seed=1)
    assert np.allclose(ci[2], 0.0, atol=1e-12) and np.allclose(ci[3], 0.0, atol=1e-12) and np.allclose(ci[6], 0.0, atol=1e-12)
    same = np.repeat(r1[:1], 16, axis=0)
    est, ci, sd, nbad = E.stack_bootstrap(same, same, ref, [8, 8], B=50, seed=1)
    assert np.allclose(sd, 0.0, atol=1e-12)
    # a patch with the stars of one kind and another with the stars of a second kind: stratification keeps the 8 + 8 split, so the pooled median stack never leaves the mix of the two
    a, b = np.repeat(r1[:1], 8, axis=0), np.repeat(r1[8:9] * 1.2, 8, axis=0)
    est, ci, sd, nbad = E.stack_bootstrap(np.concatenate([a, b]), np.concatenate([a, b]), ref, [8, 8], B=60, seed=2)
    assert np.allclose(sd, 0.0, atol=1e-12) and nbad == 0                                                 # always 8 + 8 stars, so the median stack is always the mean of the two kinds
    arrays = E.stack_arrays(np.concatenate([a, b]), np.concatenate([a, b]), ref, [8, 8], seed=2, B=20)
    assert sorted(arrays) == ['ci', 'est', 'n', 'nan', 'real', 'ref', 'sd', 'sizes', 'syn'] and arrays['n'][0] == 16 and arrays['real'].shape == (19, 19)


# ----------------------------------------------------------------------------------------------------------------------------------
# item 3: the bright isolated stars, their units, the bright pairs
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_bright_isolated_rule_on_a_planted_truth_list():
    rows = [(2.0, 50, 19.0), (1.9, 60, 19.0), (97.0, 50, 19.0), (97.1, 60, 19.0),                # the scored area: at least 2.5 px from the edge, x >= 2.0 and x <= 97.0
            (50, 2.0, 19.0), (60, 1.9, 19.0), (70, 97.0, 19.0), (80, 97.1, 19.0),
            (20, 20, 20.49), (30, 20, 20.5),                                                       # F606W < 20.5
            (20, 40, 19.0), (22.9, 40, 21.4),                                                      # a neighbour at 2.9 px brighter than 19 + 2.5: excluded
            (40, 40, 19.0), (42.9, 40, 21.6),                                                      # fainter than 19 + 2.5: kept
            (60, 40, 19.0), (63.1, 40, 15.0),                                                      # a bright neighbour beyond 3 px: kept
            (20, 80, 19.0), (20, 80, 25.0),                                                        # a faint source at the same position: kept
            (80, 20, 19.0), (82.9, 20, 21.5)]                                                      # a neighbour at exactly m + 2.5 is not brighter than m + 2.5: kept
    a = np.array(rows)
    idx = E.bright_isolated(a[:, 0], a[:, 1], a[:, 2])
    chosen = sorted((a[i, 0], a[i, 1]) for i in idx)
    # (22.9, 40) and (42.9, 40) are too faint to be stars; (63.1, 40) is a star: its neighbour (60, 40) is beyond 3 px; (20, 40) has a neighbour at 2.9 px brighter than 19 + 2.5
    assert chosen == sorted([(2.0, 50), (97.0, 50), (50, 2.0), (70, 97.0), (20, 20), (40, 40), (60, 40), (63.1, 40), (20, 80), (80, 20)]), chosen
    assert E.BRIGHT_MAG == 20.5 and E.BRIGHT_RAD == 3.0 and E.BRIGHT_DM == 2.5


def counts(m):
    return A.mag_to_counts(m)


def planted_catalogue():
    """Seven isolated bright stars in a 100 x 100 patch and a catalogue with a known fate for each; returns (universe, box, det)."""
    truth = [(30.0, 30.0, 19.0), (60.0, 30.0, 19.0), (30.0, 70.0, 19.0), (60.0, 70.0, 19.0), (80.0, 50.0, 20.4), (10.0, 90.0, 19.3), (20.0, 20.0, 19.2), (2.1, 50.0, 19.0), (1.9, 40.0, 19.0)]
    det = [(30.4, 30.45, counts(19.0)),                                         # s0: one unit, the right flux, 0.6 px off
           (60.0, 30.0, 0.6 * counts(19.0)), (61.5, 30.0, 0.4 * counts(19.0)),  # s1: split in two significant units (1.5 px apart)
           (30.0, 70.0, 0.9 * counts(19.0)), (30.0, 71.5, 0.1 * counts(19.0)),     # s2: a second unit with exactly 10 % of the flux is significant
           (62.5, 70.0, counts(19.0)),                                          # s3: the only unit is 2.5 px away: missed
           (80.5, 50.0, 1.3 * counts(20.4)),                                    # s4: faint star (not a bright pair), too much flux
           (10.3, 89.8, counts(19.3) * 10 ** (-0.4 * 0.25)),                    # s5: 0.3 px, -0.2 px off and 0.25 mag too faint
           (20.0, 20.0, counts(19.2)), (20.0, 21.0, 0.05 * counts(19.2)),       # s6: a second unit with 5 % of the flux is not significant
           (1.9, 50.0, counts(19.0)),                                           # s7: a star at the edge of the scored area, its unit just outside it: a unit, no pair (the detection is not scored)
           (2.1, 40.0, counts(19.0))]                                           # s8: a truth source outside the scored area (not a star) with a detection inside it: no pair (the truth is not scored)
    t = np.array(truth)
    return (t[:, 0], t[:, 1], t[:, 2]), (t[:, 0], t[:, 1], t[:, 2]), tuple(np.array(det).T)


def test_units_flux_ratio_and_pairs_of_planted_catalogues():
    universe, box, det = planted_catalogue()
    a = E.bright_arrays(det, universe, box)
    assert a['nsig'].tolist() == [1, 2, 2, 0, 1, 1, 1, 1] and a['nunits'].tolist() == [1, 2, 2, 0, 1, 1, 2, 1]
    assert np.allclose(a['fratio'], [1.0, 1.0, 1.0, 0.0, 1.3, 10 ** (-0.1), 1.05, 1.0])
    p = a['pairs']                                      # truth F606W < 19.5: s0, s1 (the nearer unit), s2, s5, s6 (s4 is fainter than 19.5, s3 is missed)
    assert p.shape == (5, 3)
    expect = [(0.4, 0.45, 0.0), (0.0, 0.0, -2.5 * np.log10(0.6)), (0.0, 0.0, -2.5 * np.log10(0.9)), (0.3, -0.2, 0.25), (0.0, 0.0, 0.0)]
    assert np.allclose(sorted(map(tuple, p)), sorted(expect), atol=1e-9)
    # the neighbours of the universe decide who is a bright isolated star; the box decides the pairs
    far = (np.r_[universe[0], 30.9], np.r_[universe[1], 30.45], np.r_[universe[2], 19.2])      # a source 0.9 px from s0 that is not in the box list: s0 is no longer isolated
    b = E.bright_arrays(det, far, box)
    assert b['nsig'].tolist() == [2, 2, 0, 1, 1, 1, 1] and b['pairs'].shape == (5, 3) and np.allclose(sorted(map(tuple, b['pairs'])), sorted(expect), atol=1e-9)
    none = E.bright_arrays((np.zeros(0), np.zeros(0), np.zeros(0)), universe, box)
    assert none['nsig'].tolist() == [0] * 8 and none['pairs'].shape == (0, 3)
    assert E.UNIT_RAD == 2.0 and E.UNIT_FRAC == 0.1 and E.PAIR_RAD == 0.7 and E.PAIR_MAG == 19.5
    # a detection of non-positive flux is no unit and no pair
    d = tuple(np.vstack([np.array(det).T, [30.0, 30.0, -5.0], [30.0, 30.0, 0.0]]).T)
    assert E.bright_arrays(d, universe, box)['nunits'].tolist() == a['nunits'].tolist() and E.bright_arrays(d, universe, box)['pairs'].shape == (5, 3)


def test_bootstrap_statistics_are_stratified_by_patch_and_reproducible():
    f, ci = E.boot_fraction([np.ones(10), np.zeros(30)], lambda v: v == 1)
    assert f == 0.25 and np.allclose(ci, 0.25)                                  # every stratified resample keeps 10 + 30: unstratified resamples would vary
    f, ci = E.boot_fraction([np.array([1, 0, 1, 1]), np.array([], int), np.array([0, 0])], lambda v: v >= 1)
    assert f == 0.5 and ci[0] <= f <= ci[1] and (f, ci[0], ci[1]) == (E.boot_fraction([np.array([1, 0, 1, 1]), np.array([0, 0])], lambda v: v >= 1)[0],) + tuple(
        E.boot_fraction([np.array([1, 0, 1, 1]), np.array([], int), np.array([0, 0])], lambda v: v >= 1)[1])
    m, ci = E.boot_median([np.full(10, 2.0), np.full(30, 4.0)])
    assert m == 4.0 and np.allclose(ci, 4.0)
    m, ci = E.boot_median([np.arange(21.0), np.arange(21.0)])
    assert m == 10.0 and ci[0] < 10 < ci[1] and np.array_equal(E.boot_median([np.arange(21.0), np.arange(21.0)])[1], ci)
    assert np.isnan(E.boot_median([np.zeros(0)])[0]) and np.isnan(E.boot_fraction([np.zeros(0)], lambda v: v > 0)[0])
    # nothing to resample is skipped quietly: not a number, and no warning (a median of nothing would warn)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        m, ci = E.boot_median([np.zeros(0), np.zeros(0)])
        f, cf = E.boot_fraction([np.zeros(0, int), np.zeros(0, int)], lambda v: v > 0)
    assert np.isnan(m) and np.isnan(ci).all() and np.isnan(f) and np.isnan(cf).all()


# ----------------------------------------------------------------------------------------------------------------------------------
# the runs on a world and the bookkeeping
# ----------------------------------------------------------------------------------------------------------------------------------
def fake_dets(truth_patch, shift=0.1):
    """A catalogue per readout: every scored truth source detected at its position + shift with its flux, and a few fakes."""
    x, y, m = truth_patch['x'], truth_patch['y'], truth_patch['mag']
    f = A.mag_to_counts(m)
    dets = {}
    for k in X.KAPPAS + (0.0,):
        keep = m < 24 - 0.2 * (k if k else 3)
        dets[k] = (np.r_[x[keep] + shift, 50.0], np.r_[y[keep], 51.0], np.r_[f[keep], 100.0])
    return dets


def test_iss_arrays_hold_the_scores_the_catalogue_at_the_reported_kappa_and_the_run_records():
    truth = planted_truth([(30.0, 30.0, 19.0), (60.0, 60.0, 21.0), (70.0, 20.0, 22.5), (40.0, 80.0, 23.3)])          # the last one is in the catalogue of kappa 2 and not in that of kappa 5
    t = A.truth_for_patch(truth, 0, 0)
    dets = fake_dets(t)
    assert len(dets[2.0][0]) == len(dets[5.0][0]) + 1
    arrays = E.iss_arrays(dets, 12.5, 700, dict(bic_index=650, n_blocks=701), t, 2.0, 1.21)
    assert sorted(k for k in arrays if k.endswith('|tot')) == sorted('k%g|tot' % k for k in X.KAPPAS + (0.0,))
    s = X.score_dets(dets[2.0], t)
    assert arrays['k2|tot'].tolist() == [s['true22'], s['matched22'], s['det22'], s['det22_matched'], s['true_all'], s['matched_all'], s['det_all']]
    seen = set()
    for k in X.KAPPAS + (0.0,):                                                                          # every readout is scored on its own catalogue (the catalogues differ in the faint end)
        sk = X.score_dets(dets[k], t)
        assert arrays['k%g|tot' % k].tolist() == [sk['true22'], sk['matched22'], sk['det22'], sk['det22_matched'], sk['true_all'], sk['matched_all'], sk['det_all']], k
        seen.add(tuple(arrays['k%g|tot' % k]))
    assert len(seen) > 1
    assert np.array_equal(arrays['k2|blk22'], s['blk22']) and arrays['k2|blk22'].shape == (25, 4)
    assert np.array_equal(arrays['k2|det'], np.c_[dets[2.0][0], dets[2.0][1], dets[2.0][2]]) and 'k5|det' not in arrays
    assert arrays['secs'].tolist() == [12.5, 700] and arrays['bic'].tolist() == [650, 701, 0] and arrays['rsd'].tolist() == [1.21]
    assert E.iss_arrays(dets, 1.0, X.MAX_ENTRIES, dict(bic_index=1, n_blocks=3), t, 5.0, 1.0)['bic'].tolist() == [1, 3, 1]
    a5 = E.iss_arrays(dets, 1.0, 5, dict(bic_index=1, n_blocks=3), t, 5.0, 1.0)
    assert np.array_equal(a5['k5|det'], np.c_[dets[5.0][0], dets[5.0][1], dets[5.0][2]]) and 'k2|det' not in a5 and len(a5['k5|det']) == len(arrays['k2|det']) - 1       # the catalogue is the one of the reported kappa


def test_crowdsource_default_and_frozen_post_filter_on_a_world_catalogue():
    truth = planted_truth([(30.0, 30.0, 19.0), (60.0, 60.0, 21.0)])
    t = A.truth_for_patch(truth, 0, 0)
    f = A.mag_to_counts(np.array([19.0, 21.0, 22.0, 22.5, 21.5, 21.7, 21.8]))
    nan = np.nan
    cat = dict(x=np.array([30.0, 60.0, 10.0, 20.0, 40.0, 50.0, 55.0]), y=np.array([30.0, 60.0, 10.0, 20.0, 40.0, 50.0, 55.0]), flux=f, dflux=f / np.array([50, 20, 4.0, nan, 9.0, 20.0, 20.0]),
               qf=np.array([0.9, 0.8, 0.7, 0.7, nan, -0.5, 0.7]), fracflux=np.array([1.0, 0.9, 0.5, 0.5, 0.5, 0.5, -0.2]), spread=np.zeros(7), rchi2=np.ones(7))
    arrays = E.cs_arrays(cat, (7, 8), t, 'crowdsource|snr5|ff0|qf0|spinf', 3.5)
    # the default filter (0, 0, 0, inf) drops a source with a NaN error or quality flag and one with a negative quality flag or fraction of flux; the frozen post-filter also needs flux >= 5 dflux
    assert arrays['default|det'][:, 0].tolist() == [30.0, 60.0, 10.0] and arrays['tuned|det'][:, 0].tolist() == [30.0, 60.0]
    assert arrays['secs'].tolist() == [3.5, 7]
    assert arrays['default|tot'][2] >= arrays['tuned|tot'][2] and arrays['default|det'].shape[0] > arrays['tuned|det'].shape[0]
    assert E.cs_filter('crowdsource|snr5|ff0|qf0|spinf') == (5.0, 0.0, 0.0, np.inf)
    b = E.cs_dict(cat, (7, 8))
    assert sorted(b)[0].startswith('x7_y8|crowdsource|') and X.cs_dets(b, (7, 8), (0.0, 0.0, 0.0, np.inf))[0].shape[0] == arrays['default|det'].shape[0]
    p = E.phot_arrays((cat['x'], cat['y'], cat['flux']), t, 9.0)
    assert p['secs'].tolist() == [9.0, 7] and p['tuned|det'].shape == (7, 3) and 'tuned|blk22' in p


def test_world_job_scores_against_the_injected_truth_and_runs_the_frozen_functions(monkeypatch):
    rows = [(30.3, 30.7, 19.0), (60.2, 61.4, 20.5), (75.5, 25.5, 22.0), (45.0, 80.0, 21.0), (2.4, 50.0, 19.0)]          # the last one is 0.4 px inside the scored area: a box moved by a pixel loses it
    truth = planted_truth(rows, 10, 20)
    patch = dict(dn=np.zeros((100, 100)), psf=gauss_stamp(1.5), gain=4.6, dark=1.0)
    ramp = 150.0 + np.tile(np.linspace(0.0, 10.0, 100), (100, 1)) ** 2 / 10.0          # a sky that is not flat: mean 153.3, median 152.5, maximum 160
    monkeypatch.setattr(E, 'load_inputs', lambda name, pxy: (patch, ramp, truth, 5.5))
    seen = {}

    def fake_run_path(world, h, r, **kw):
        seen['iso'] = (h, r, world['dn'].shape, float(world['dn'].mean()))
        seen['iso_dn'] = world['dn'].copy()
        t = A.truth_for_patch(truth, 10, 20)
        return fake_dets(t), 1.5, 100, dict(bic_index=90, n_blocks=101)

    def fake_run_psf(world, h, r, fit):
        seen['adq'] = (h, r, fit['kind'])
        t = A.truth_for_patch(truth, 10, 20)
        return fake_dets(t, 0.2), 2.5, 120, dict(bic_index=100, n_blocks=121, rsd={k: 1.1 + 0.01 * i for i, k in enumerate(X.KAPPAS + (0.0,))})

    def fake_cs(dn, sky, psf, gain, dark):
        seen['cs'] = (dn.shape, sky, gain, dark)
        return dict(x=np.array([30.3, 60.2]), y=np.array([30.7, 61.4]), flux=A.mag_to_counts(np.array([19.0, 20.5])), dflux=np.ones(2), qf=np.ones(2), fracflux=np.ones(2), spread=np.zeros(2), rchi2=np.ones(2)), 4.0

    def fake_phot(dn, psf, gain, dark, skysig, spec):
        seen['phot'] = (skysig, spec['tag'])
        return dict(x=np.array([30.3]), y=np.array([30.7]), flux=A.mag_to_counts(np.array([19.0]))), 6.0, 2.5
    monkeypatch.setattr(X, 'run_path', fake_run_path)
    monkeypatch.setattr(E.PS, 'run_path_psf', fake_run_psf)
    monkeypatch.setattr(E.PS, 'residual_sd', lambda world, P, det: 1.07)
    monkeypatch.setattr(X.AB, 'run_crowdsource', fake_cs)
    monkeypatch.setattr(X.AB, 'run_photutils_ext', fake_phot)
    cfg = dict(hp=3.0, rex=1.2, kap=2.0, tuned=dict(crowdsource='crowdsource|snr5|ff0|qf0|spinf', photutils='photutils_ext_all_local_i6'))
    kind, name, pxy, arrays = E.world_job(('iso', 'm92', (10, 20), cfg))
    assert (kind, name, pxy) == ('iso', 'm92', (10, 20)) and seen['iso'][:3] == (3.0, 1.2, (100, 100)) and seen['iso'][3] > 140.0           # the frame is the world (sky 150 DN), not the patch
    assert 'k2|det' in arrays and arrays['rsd'].tolist() == [1.07] and arrays['secs'].tolist() == [1.5, 100]
    w = E.render_world(truth, 10, 20, patch['psf'], 4.6, 1.0, ramp, E.world_seed('m92', (10, 20)))
    assert np.array_equal(seen['iso_dn'], w['dn']) and not np.array_equal(w['dn'], w['clean'])                                                     # the path runs on the noisy world of the patch's seed
    t = A.truth_for_patch(truth, 10, 20)
    assert arrays['k2|tot'].tolist()[0] == int(((t['mag'] < 22) & t['inner']).sum()) == 4                                                          # scored against the injected truth (the box of the patch)
    assert E.world_job(('iso_def', 'm92', (10, 20), cfg))[3]['k5|det'].shape[1] == 3 and seen['iso'][:2] == (X.DEFAULT[0], X.DEFAULT[1])
    arr = E.world_job(('adq', 'm92', (10, 20), cfg))[3]
    assert seen['adq'][:2] == (3.0, 1.2) and 'k2|det' in arr and arr['rsd'].tolist() == [1.1 + 0.01 * X.KAPPAS.index(2.0)] and arr['secs'].tolist() == [2.5, 120]
    assert E.world_job(('adq_def', 'm92', (10, 20), cfg))[3]['rsd'].tolist() == [1.1 + 0.01 * X.KAPPAS.index(5.0)] and seen['adq'][:2] == (X.DEFAULT[0], X.DEFAULT[1])
    arr = E.world_job(('cs', 'm92', (10, 20), cfg))[3]
    assert seen['cs'][0] == (100, 100) and abs(seen['cs'][1] - ramp.mean()) < 1e-9 and seen['cs'][2:] == (4.6, 1.0) and arr['default|det'].shape == (2, 3)
    arr = E.world_job(('phot', 'm92', (10, 20), cfg))[3]
    assert seen['phot'] == (5.5, 'photutils_ext_all_local_i6') and arr['tuned|det'].shape == (1, 3)
    with pytest.raises(ValueError):
        E.world_job(('nothing', 'm92', (10, 20), cfg))


def test_the_frozen_configuration_is_read_after_the_freeze_is_checked(monkeypatch):
    events = []
    monkeypatch.setattr(X, 'load_raw', lambda p: events.append(('load', p)) or {'frozen': 1})
    monkeypatch.setattr(X, 'check_freeze', lambda raw: events.append(('freeze', raw)))
    monkeypatch.setattr(X, 'choose', lambda raw: (3.0, 1.2, 2.0, 'f'))                                 # sigma_b, r_ex, kappa, the F
    monkeypatch.setattr(X, 'choose_baselines', lambda raw: dict(crowdsource='c', photutils='p'))
    cfg = E.run_cfg('frozen.npz')
    assert cfg == dict(hp=3.0, rex=1.2, kap=2.0, tuned=dict(crowdsource='c', photutils='p')) and events == [('load', 'frozen.npz'), ('freeze', {'frozen': 1})]


def test_jobs_are_run_in_the_declared_order_once_and_the_adequate_run_is_skipped_where_the_rule_keeps_the_isotropic_kernel():
    cur = {}
    todo = E.jobs_todo(cur)
    stage1 = [(k, n, p) for k, n, p in todo[:12]]
    assert stage1[:6] == [('iso', 'm2', (310, 630)), ('adq', 'm2', (310, 630)), ('cs', 'm2', (310, 630)), ('iso', 'm53', (660, 420)), ('adq', 'm53', (660, 420)), ('cs', 'm53', (660, 420))]
    assert len(todo) == 10 * 3 + 10 + 10 * 2 and todo[-1][0] == 'adq_def' and [t[0] for t in todo[30:40]] == ['phot'] * 10
    assert [t for t in todo if t[1] in R.REPEAT][0][0] == 'iso' and todo.index(('iso', 'm2r1', (316, 124))) > todo.index(('cs', 'm92', (1430, 550)))        # the referee's four first, then the repeat patches
    cur = {'WH|m2|x310_y630|kind': np.array(['isotropic']), 'WH|m53|x660_y420|kind': np.array(['general']), 'W|m2|x310_y630|iso|ours|secs': np.array([1.0, 2]), 'W|m53|x660_y420|cs|secs': np.array([1.0, 2])}
    todo = E.jobs_todo(cur)
    assert ('iso', 'm2', (310, 630)) not in todo and ('adq', 'm2', (310, 630)) not in todo and ('adq_def', 'm2', (310, 630)) not in todo and ('cs', 'm2', (310, 630)) in todo
    assert ('adq', 'm53', (660, 420)) in todo and ('cs', 'm53', (660, 420)) not in todo and ('iso', 'm53', (660, 420)) in todo
    assert E.jobs_todo(cur, kinds=['cs']) == [t for t in todo if t[0] == 'cs'] and all(t[1] == 'm92' for t in E.jobs_todo(cur, sets=['m92']))
    assert E.jobs_todo(cur, kinds=['phot'], sets=['m2']) == [('phot', 'm2', (310, 630))]


def test_part_synth_checks_the_freeze_stores_every_result_under_its_key_and_does_not_run_a_result_again(monkeypatch, tmp_path):
    calls, events = [], []
    monkeypatch.setattr(E, 'SYNTH', [('m2', (1, 2)), ('m92', (3, 4))])
    monkeypatch.setattr(E, 'STAGES', [(('iso', 'adq', 'cs'), [('m2', (1, 2)), ('m92', (3, 4))])])
    monkeypatch.setattr(E, 'run_cfg', lambda m2raw: events.append('freeze') or dict(hp=3.0, rex=1.2, kap=2.0, tuned=dict(crowdsource='c', photutils='p')))
    monkeypatch.setattr(E, 'meta_arrays', lambda cfg: dict(cfg=np.array([cfg['hp']])))
    monkeypatch.setattr(E, 'part_real', lambda out, *a: events.append('real'))
    headers = []
    monkeypatch.setattr(E, 'world_header', lambda name, pxy: headers.append((name, pxy)) or dict(kind=np.array(['isotropic' if name == 'm2' else 'general']), marker=np.array([pxy[0]])))

    def fake_job(args):
        kind, name, pxy, cfg = args
        calls.append((kind, name, pxy))
        return kind, name, pxy, {'secs': np.array([1.0 + len(calls), 5]), 'k2|det': np.zeros((2, 3))}
    monkeypatch.setattr(E, 'world_job', fake_job)

    class FakePool:
        made = []

        def __init__(self, n):
            FakePool.made.append(n)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def imap_unordered(self, fn, jobs, chunksize=1):
            return iter([fn(j) for j in jobs])
    import multiprocessing
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    out = str(tmp_path / 'diag_raw.npz')
    E.part_synth(out, workers=1)
    assert events == ['freeze', 'real'] and calls == [('iso', 'm2', (1, 2)), ('cs', 'm2', (1, 2)), ('iso', 'm92', (3, 4)), ('adq', 'm92', (3, 4)), ('cs', 'm92', (3, 4))]
    raw = X.load_raw(out)
    assert raw['META|cfg'].tolist() == [3.0] and raw['WH|m92|x3_y4|marker'].tolist() == [3] and raw['WH|m2|x1_y2|kind'].tolist() == ['isotropic']
    assert 'W|m2|x1_y2|iso|ours|secs' in raw and 'W|m92|x3_y4|adq|ours|k2|det' in raw and 'W|m92|x3_y4|cs|secs' in raw and 'W|m2|x1_y2|adq|ours|secs' not in raw
    calls.clear()
    assert headers == [('m2', (1, 2)), ('m92', (3, 4))]
    E.part_synth(out, workers=1)                                            # nothing is run again
    assert calls == [] and events.count('freeze') == 2 and len(headers) == 2                 # nor a header
    E.part_synth(out, workers=3, kinds=['cs'], sets=['m92'])                # nothing to do for this filter either, and no pool is made for it
    assert FakePool.made == []
    raw = X.load_raw(out)
    del raw['W|m92|x3_y4|cs|secs']
    np.savez_compressed(out, **raw)
    E.part_synth(out, workers=3)
    assert calls == [('cs', 'm92', (3, 4))] and FakePool.made == []        # one job: run in this process, not in a pool
    raw = X.load_raw(out)
    for k in [k for k in raw if k.startswith('W|m92|x3_y4|')]:
        del raw[k]
    np.savez_compressed(out, **raw)
    calls.clear()
    E.part_synth(out, workers=2)
    assert sorted(calls) == [('adq', 'm92', (3, 4)), ('cs', 'm92', (3, 4)), ('iso', 'm92', (3, 4))] and FakePool.made == [2]            # three jobs, two workers
    raw = X.load_raw(out)
    for k in [k for k in raw if k.startswith('W|m92|x3_y4|')]:
        del raw[k]
    np.savez_compressed(out, **raw)
    E.part_synth(out, workers=8)
    assert FakePool.made == [2, 3]                                                                                               # eight workers, three jobs: a pool of three


def planted_base(pxy, det, flux_cut=0.0, low=()):
    """The baseline catalogue arrays of a patch in the key format of baseline_catalogs.npz (crowdsource with its quality columns, photutils' tuned variant): `det` = (x, y, flux); the detections
    listed in `low` have a signal to noise of 1 (the frozen post-filter, snr >= 5, removes them; the default keeps them), the others 50."""
    x, y, f = [np.asarray(a, float) for a in det]
    k = 'x%d_y%d|' % pxy
    n = len(x)
    keep = f >= flux_cut
    df = f / 50.0
    df[list(low)] = f[list(low)]
    return {k + 'crowdsource|x': x[keep], k + 'crowdsource|y': y[keep], k + 'crowdsource|flux': f[keep], k + 'crowdsource|dflux': df[keep], k + 'crowdsource|fracflux': np.ones(keep.sum()),
            k + 'crowdsource|qf': np.ones(keep.sum()), k + 'crowdsource|spread': np.zeros(keep.sum()), k + 'photutils_ext_all_local_i6|x': x[:2], k + 'photutils_ext_all_local_i6|y': y[:2],
            k + 'photutils_ext_all_local_i6|flux': f[:2]}


def test_the_bright_star_part_reads_the_real_catalogues_of_each_family_and_the_worlds_and_aliases_the_isotropic_kernel(monkeypatch, tmp_path):
    universe, box, det = planted_catalogue()
    x, y, f = det
    rows = list(zip(*universe)) + [(96.9, 50.0, 19.0), (99.6, 50.0, 16.0)]                              # a star in the scored area with a much brighter neighbour 2.7 px away, outside the box: no bright isolated star
    one = (x[[0, 1, 3, 6, 7, 8, 10]], y[[0, 1, 3, 6, 7, 8, 10]], f[[0, 1, 3, 6, 7, 8, 10]])        # one unit per star, no split, s3 missed (s1's first unit, s2's first, s4's, s5's, s6's first, s7's)
    first = (x[:2], y[:2], f[:2])
    patches = {'m92': (5, 6), 'm2': (9, 10), 'm2r2': (7, 8)}
    monkeypatch.setattr(E, 'BRIGHT_SETS', {n: [p] for n, p in patches.items()})
    monkeypatch.setattr(E, 'truth_of', lambda name: planted_truth(rows, *patches[name]))

    def d(a):
        return np.c_[a[0], a[1], a[2]]
    ok, cap = np.array([100, 1200, 0]), np.array([100, 2501, 1])
    m2raw = {'T2|m92|x5_y6|ours|det': d(det), 'T2|m92|x5_y6|default|det': d(first), 'T|x9_y10|ours|det': d(det), 'T|x9_y10|default|det': d(first), 'T2|m92|x5_y6|ours|bic': ok, 'T2|m92|x5_y6|default|bic': cap,
             'T|x9_y10|ours|bic': ok, 'T|x9_y10|default|bic': ok}
    psf = {'S|m92|x5_y6|kind': np.array(['general']), 'S|m2|x9_y10|kind': np.array(['isotropic']), 'P|m92|x5_y6|ours|det': d(one), 'P|m92|x5_y6|default|det': d(first), 'P|m92|x5_y6|ours|bic': ok, 'P|m92|x5_y6|default|bic': cap}
    rep = {'R|m2r2|x7_y8|kind': np.array(['general']), 'R|m2r2|x7_y8|iso|ours|det': d(first), 'R|m2r2|x7_y8|iso|default|det': d(first), 'R|m2r2|x7_y8|adq|ours|det': d(one), 'R|m2r2|x7_y8|adq|default|det': d(one),
           'R|m2r2|x7_y8|iso|ours|bic': ok, 'R|m2r2|x7_y8|iso|default|bic': ok, 'R|m2r2|x7_y8|adq|ours|bic': cap, 'R|m2r2|x7_y8|adq|default|bic': ok}
    src = dict(m2raw=m2raw, psf=psf, rep=rep)
    tuned = dict(crowdsource='crowdsource|snr5|ff0|qf0|spinf', photutils='photutils_ext_all_local_i6')
    bases = dict(m2=planted_base((9, 10), first), second=planted_base((5, 6), one, low=(0,)), repeat=planted_base((7, 8), det), tuned=tuned)
    monkeypatch.setattr(E, 'real_sources', lambda *a: src)
    out = str(tmp_path / 'diag_raw.npz')
    cur = {'W|m92|x5_y6|iso|ours|k2|det': d(first), 'W|m92|x5_y6|cs|default|det': d(one), 'W|m92|x5_y6|cs|tuned|det': d(one), 'W|m92|x5_y6|phot|tuned|det': d(first),
           'W|m2|x9_y10|iso|ours|k2|det': d(one), 'W|m92|x5_y6|iso|ours|bic': ok, 'W|m2|x9_y10|iso|ours|bic': cap,
           'W|m92|x5_y6|adq|ours|k2|det': d(first), 'W|m92|x5_y6|adq|ours|bic': cap}
    np.savez_compressed(out, **cur)
    E.part_bright(out, bases=bases)
    raw = X.load_raw(out)
    assert all(k in raw for k in cur)                                                                    # the worlds' results stay
    nsig = lambda w, n, p, s: raw['B|%s|%s|x%d_y%d|%s|nsig' % (w, n, p[0], p[1], s)].tolist()
    # real, frame family M92 (general kernel): the tuned and default isotropic catalogues, the adequate ones from the kernel study, crowdsource and photutils from the baseline catalogue
    assert nsig('real', 'm92', (5, 6), 'iso|ours|k2') == [1, 2, 2, 0, 1, 1, 1, 1] and nsig('real', 'm92', (5, 6), 'iso|default|k5') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('real', 'm92', (5, 6), 'adq|ours|k2') == [1, 1, 1, 0, 1, 1, 1, 1] and nsig('real', 'm92', (5, 6), 'adq|default|k5') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('real', 'm92', (5, 6), 'cs|default') == [1, 1, 1, 0, 1, 1, 1, 1] and nsig('real', 'm92', (5, 6), 'phot|tuned') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('real', 'm92', (5, 6), 'cs|tuned') == [0, 1, 1, 0, 1, 1, 1, 1]                           # the frozen post-filter drops the detection of s0 (signal to noise 1), the default keeps it
    assert str(raw['B|real|m92|x5_y6|kind'][0]) == 'general' and raw['B|real|m92|x5_y6|iso|ours|k2|pairs'].shape == (5, 3)
    # real, M2 (isotropic kernel, no adequate catalogue), repeat family (all four ISS catalogues from the repeat raw)
    assert 'B|real|m2|x9_y10|adq|ours|k2|nsig' not in raw and nsig('real', 'm2', (9, 10), 'iso|ours|k2') == [1, 2, 2, 0, 1, 1, 1, 1] and nsig('real', 'm2', (9, 10), 'cs|default') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('real', 'm2r2', (7, 8), 'adq|ours|k2') == [1, 1, 1, 0, 1, 1, 1, 1] and nsig('real', 'm2r2', (7, 8), 'iso|default|k5') == [1, 1, 0, 0, 0, 0, 0, 0] and nsig('real', 'm2r2', (7, 8), 'cs|default') == [1, 2, 2, 0, 1, 1, 1, 1]
    # synthetic: the worlds' catalogues where they exist; the adequate catalogue of an isotropic-kernel patch is its isotropic one
    assert nsig('synth', 'm92', (5, 6), 'iso|ours|k2') == [1, 1, 0, 0, 0, 0, 0, 0] and nsig('synth', 'm92', (5, 6), 'cs|default') == [1, 1, 1, 0, 1, 1, 1, 1] and nsig('synth', 'm92', (5, 6), 'adq|ours|k2') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('synth', 'm92', (5, 6), 'phot|tuned') == [1, 1, 0, 0, 0, 0, 0, 0]
    assert nsig('synth', 'm2', (9, 10), 'iso|ours|k2') == nsig('synth', 'm2', (9, 10), 'adq|ours|k2') == [1, 1, 1, 0, 1, 1, 1, 1]
    assert not any(k.startswith('B|synth|m2r2|') for k in raw)
    # the stalled paths (entry cap reached) are flagged beside the statistics, for each family of frozen results and for the worlds
    flag = lambda w, n, p, s: raw['B|%s|%s|x%d_y%d|%s|stalled' % (w, n, p[0], p[1], s)].tolist()
    assert flag('real', 'm92', (5, 6), 'iso|ours|k2') == [0] and flag('real', 'm92', (5, 6), 'iso|default|k5') == [1] and flag('real', 'm92', (5, 6), 'adq|ours|k2') == [0] and flag('real', 'm92', (5, 6), 'adq|default|k5') == [1]
    assert flag('real', 'm2', (9, 10), 'iso|ours|k2') == [0] and flag('real', 'm2', (9, 10), 'iso|default|k5') == [0]
    assert flag('real', 'm2r2', (7, 8), 'iso|ours|k2') == [0] and flag('real', 'm2r2', (7, 8), 'adq|ours|k2') == [1] and flag('real', 'm2r2', (7, 8), 'adq|default|k5') == [0]
    assert flag('synth', 'm92', (5, 6), 'iso|ours|k2') == [0] and flag('synth', 'm92', (5, 6), 'adq|ours|k2') == [1]                      # a general-kernel patch: the two runs have their own records
    assert flag('synth', 'm2', (9, 10), 'iso|ours|k2') == [1] and flag('synth', 'm2', (9, 10), 'adq|ours|k2') == [1]                      # an isotropic-kernel patch: the adequate run is the isotropic one
    assert flag('real', 'm92', (5, 6), 'cs|default') == [0] and flag('synth', 'm92', (5, 6), 'cs|default') == [0] and flag('real', 'm92', (5, 6), 'phot|tuned') == [0] and flag('synth', 'm92', (5, 6), 'phot|tuned') == [0]      # no path, no flag
    # a repeat set is skipped where the frozen results of the repeat patches are absent; a patch the kernel study does not hold is skipped
    monkeypatch.setattr(E, 'real_sources', lambda *a: dict(src, rep=None))
    out2 = str(tmp_path / 'raw2.npz')
    np.savez_compressed(out2, **cur)
    E.part_bright(out2, bases=dict(bases, repeat=None))
    assert not any(k.startswith('B|real|m2r2|') for k in X.load_raw(out2)) and any(k.startswith('B|real|m92|') for k in X.load_raw(out2))
    # the 'real_dets' of a repeat patch whose baseline catalogue file is absent is empty
    assert E.real_dets(src, dict(bases, repeat=None), 'm2r2', (7, 8)) == {}


def test_a_path_that_reached_the_entry_cap_is_stalled_and_the_run_is_found_from_the_method_suffix():
    cap = X.MAX_ENTRIES
    assert E.is_stalled(np.array([293, cap + 1, 1])) and E.is_stalled(np.array([5, cap + 1, 0])) and E.is_stalled(np.array([5, 100, 1]))          # the flag, or the blocks of the path ending in the cap + 1
    assert not E.is_stalled(np.array([1338, 1339, 0])) and not E.is_stalled(np.array([5, cap, 0])) and not E.is_stalled(None)
    assert E.path_tag('iso|ours|k2') == 'iso|ours' and E.path_tag('adq|default|k5') == 'adq|default' and E.path_tag('iso|default') == 'iso|default'
    assert E.path_tag('cs|default') is None and E.path_tag('cs|tuned') is None and E.path_tag('phot|tuned') is None
    src = planted_src()
    assert E.real_bic(src, 'm92', (1100, 220), 'iso|default|k5').tolist() == [293, 2501, 1] and E.real_bic(src, 'm92', (1100, 220), 'iso|ours|k2').tolist() == [900, 1291, 0]
    assert E.real_bic(src, 'm92', (1100, 220), 'adq|ours|k2').tolist() == [1100, 1401, 0] and E.real_bic(src, 'm2', (310, 630), 'iso|default|k5').tolist() == [700, 1500, 0]
    assert E.real_bic(src, 'm2r2', (133, 281), 'iso|default|k5').tolist() == [508, 2501, 1] and E.real_bic(src, 'm2r2', (133, 281), 'adq|ours|k2').tolist() == [600, 1201, 0]
    assert E.real_bic(src, 'm92', (1100, 220), 'cs|default') is None and E.real_bic(dict(src, rep=None), 'm2r2', (133, 281), 'iso|ours|k2') is None and E.real_bic(src, 'm92', (0, 0), 'iso|ours|k2') is None


def test_bright_statistics_of_a_set_of_mixed_kernels_read_the_isotropic_rows_for_the_isotropic_patch_and_need_every_patch():
    def put(raw, world, name, pxy, s, nsig, kind=None):
        k = 'B|%s|%s|x%d_y%d|%s|' % (world, name, pxy[0], pxy[1], s)
        raw[k + 'nsig'], raw[k + 'fratio'], raw[k + 'pairs'] = np.array(nsig), np.ones(len(nsig)), np.zeros((2, 3))
        if kind:
            raw['WH|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1])] = np.array([kind])
    raw = {}
    put(raw, 'real', 'm92', (1, 1), 'iso|ours|k2', [1, 1, 1], 'isotropic')                                  # an isotropic-kernel patch: only the isotropic catalogue
    put(raw, 'real', 'm92', (2, 2), 'iso|ours|k2', [1, 1], 'general')
    put(raw, 'real', 'm92', (2, 2), 'adq|ours|k2', [2, 2, 2, 2, 0], None)                                    # a general-kernel patch: both
    st = E.bright_stats(raw, 'real', 'm92', [(1, 1), (2, 2)], 'adq|ours|k2')
    assert st['n'] == 3 + 5 and abs(st['split'][0] - 4 / 8) < 1e-12 and abs(st['none'][0] - 1 / 8) < 1e-12                  # the isotropic patch's isotropic stars and the other patch's adequate ones
    st = E.bright_stats(raw, 'real', 'm92', [(1, 1), (2, 2)], 'iso|ours|k2')
    assert st['n'] == 3 + 2 and st['split'][0] == 0.0
    assert E.bright_stats(raw, 'real', 'm92', [(1, 1), (2, 2), (3, 3)], 'iso|ours|k2') is None              # a patch without the catalogue: no statistic, not a pool over the rest
    assert E.bright_stats(raw, 'real', 'm92', [], 'iso|ours|k2') is None


def test_two_results_of_bright_stats_are_the_same_only_if_every_number_is():
    a = dict(n=10, pairs=5, excluded=[], split=(0.5, np.array([0.2, 0.8])), none=(0.0, np.array([0.0, 0.0])), fratio=(1.0, np.array([0.9, 1.1])), dmag=(np.nan, np.array([np.nan, np.nan])),
             dx=(0.1, np.array([0.0, 0.2])), dy=(0.1, np.array([0.0, 0.2])))
    import copy
    b = copy.deepcopy(a)
    assert E.same_stats(a, b) and not E.same_stats(a, None) and not E.same_stats(None, b)
    for k, v in (('n', 11), ('pairs', 6), ('split', (0.5, np.array([0.2, 0.7]))), ('dy', (0.2, np.array([0.0, 0.2]))), ('fratio', (1.0, np.array([0.9, 1.2]))), ('excluded', [(1, 2)])):
        c = copy.deepcopy(a)
        c[k] = v
        assert not E.same_stats(a, c), k
    assert not E.same_stats(dict(n=None, excluded=[(1, 2)]), dict(n=None, excluded=[(1, 2)]))                  # two results without a statistic are not 'the same'
    assert not E.same_stats(a, dict(n=None, excluded=[(1, 2)]))


def planted_src():
    """Frozen-result dicts with distinguishable blocks: F = (matched / 10) per method, the key layout of m2_raw, psf_raw and repeat_raw."""
    def pack(name, f10, mult=1):
        blk = np.zeros((25, 4))
        blk[:, 0], blk[:, 1], blk[:, 2], blk[:, 3] = 10, f10, 10, f10
        return {name + '|blk22': blk * mult, name + '|tot': np.array([250, 25 * f10, 250, 25 * f10, 300, 25 * f10, 260]) * mult}
    m2raw, psf, rep = {}, {}, {}
    for pre, vals in (('T|x310_y630|', dict(iso=6, dflt=5, cs=4, csT=3, ph=2)), ('T2|m92|x1100_y220|', dict(iso=2, dflt=1, cs=3, csT=3, ph=2)), ('T2|m53|x660_y420|', dict(iso=3, dflt=2, cs=2, csT=2, ph=1))):
        for k, v in (('ours|k2', vals['iso']), ('default|k5', vals['dflt']), ('crowdsource (default)', vals['cs']), ('crowdsource (tuned)', vals['csT']), ('photutils (tuned)', vals['ph'])):
            m2raw.update({pre + kk: vv for kk, vv in pack(k, v).items()})
        m2raw[pre + 'ours|det'] = np.c_[np.arange(5.0), np.arange(5.0) + 1, np.full(5, 1000.0)]
        m2raw[pre + 'default|det'] = np.c_[np.arange(4.0), np.arange(4.0) + 1, np.full(4, 900.0)]
        m2raw[pre + 'ours|bic'] = np.array([900, 1291, 0])
        m2raw[pre + 'default|bic'] = np.array([293, 2501, 1]) if pre == 'T2|m92|x1100_y220|' else np.array([700, 1500, 0])
    psf['S|m2|x310_y630|kind'], psf['S|m92|x1100_y220|kind'], psf['S|m53|x660_y420|kind'] = np.array(['isotropic']), np.array(['general']), np.array(['general'])
    psf['R|m2|x310_y630|ours|rsd'], psf['R|m2|x310_y630|default|rsd'] = np.array([0.96]), np.array([1.01])
    for ds, p in (('m92', 'x1100_y220'), ('m53', 'x660_y420')):
        for k, v in (('ours|k2', 7), ('default|k5', 6)):
            psf.update({'P|%s|%s|%s' % (ds, p, kk): vv for kk, vv in pack(k, v).items()})
        psf['P|%s|%s|ours|rsd' % (ds, p)], psf['P|%s|%s|default|rsd' % (ds, p)] = np.array([1.84, 1.22]), np.array([2.56, 1.60])
        psf['P|%s|%s|ours|bic' % (ds, p)], psf['P|%s|%s|default|bic' % (ds, p)] = np.array([1100, 1401, 0]), np.array([1200, 1801, 0])
        psf['P|%s|%s|ours|det' % (ds, p)], psf['P|%s|%s|default|det' % (ds, p)] = np.c_[np.arange(6.0), np.arange(6.0), np.full(6, 800.0)], np.c_[np.arange(3.0), np.arange(3.0), np.full(3, 700.0)]
    for pre, kind in (('R|m2r1|x316_y124|', 'isotropic'), ('R|m2r2|x133_y281|', 'general')):
        rep[pre + 'kind'] = np.array([kind])
        for k, v in (('iso|ours|k2', 4), ('iso|default|k5', 3), ('adq|ours|k2', 5), ('adq|default|k5', 4), ('crowdsource (default)', 2), ('crowdsource (tuned)', 2), ('photutils (tuned)', 1)):
            rep.update({pre + kk: vv for kk, vv in pack(k, v).items()})
        for kern, lab, val in (('iso', 'ours', 1.07), ('iso', 'default', 1.5), ('adq', 'ours', 1.05), ('adq', 'default', 1.2)):
            rep[pre + '%s|%s|rsd' % (kern, lab)] = np.array([val])
            rep[pre + '%s|%s|bic' % (kern, lab)] = np.array([508, 2501, 1]) if (pre, kern, lab) == ('R|m2r2|x133_y281|', 'iso', 'default') else np.array([600, 1201, 0])
    return dict(m2raw=m2raw, psf=psf, rep=rep)


def test_the_real_results_are_copied_under_the_keys_of_the_worlds_with_the_adequate_rows_only_where_the_kernel_differs():
    src = planted_src()
    a = E.real_arrays(src, 'm92', (1100, 220))
    assert str(a['kind'][0]) == 'general' and a['iso|ours|k2|tot'][1] == 2 * 25 and a['adq|ours|k2|tot'][1] == 7 * 25 and a['adq|default|k5|tot'][1] == 6 * 25
    assert a['cs|default|tot'][1] == 3 * 25 and a['phot|tuned|tot'][1] == 2 * 25 and a['iso|default|k5|tot'][1] == 1 * 25 and a['cs|tuned|blk22'].shape == (25, 4)
    assert a['iso|ours|rsd'].tolist() == [1.84] and a['adq|ours|rsd'].tolist() == [1.22] and a['iso|default|rsd'].tolist() == [2.56] and a['adq|default|rsd'].tolist() == [1.6]
    assert a['iso|ours|bic'].tolist() == [900, 1291, 0] and a['iso|default|bic'].tolist() == [293, 2501, 1] and a['adq|ours|bic'].tolist() == [1100, 1401, 0] and a['adq|default|bic'].tolist() == [1200, 1801, 0]
    m2 = E.real_arrays(src, 'm2', (310, 630))
    assert str(m2['kind'][0]) == 'isotropic' and not any(k.startswith('adq') for k in m2) and m2['iso|ours|rsd'].tolist() == [0.96] and m2['iso|default|rsd'].tolist() == [1.01]
    assert m2['iso|ours|bic'].tolist() == [900, 1291, 0] and m2['iso|default|bic'].tolist() == [700, 1500, 0] and 'adq|ours|bic' not in m2
    assert m2['iso|ours|k2|tot'][1] == 6 * 25 and m2['cs|default|tot'][1] == 4 * 25 and m2['cs|tuned|tot'][1] == 3 * 25 and m2['phot|tuned|tot'][1] == 2 * 25
    r1 = E.real_arrays(src, 'm2r1', (316, 124))
    assert str(r1['kind'][0]) == 'isotropic' and not any(k.startswith('adq') for k in r1) and r1['iso|ours|k2|tot'][1] == 4 * 25 and r1['iso|default|rsd'].tolist() == [1.5]
    r2 = E.real_arrays(src, 'm2r2', (133, 281))
    assert str(r2['kind'][0]) == 'general' and r2['adq|ours|k2|tot'][1] == 5 * 25 and r2['adq|default|rsd'].tolist() == [1.2] and r2['iso|ours|rsd'].tolist() == [1.07]
    assert r2['iso|default|bic'].tolist() == [508, 2501, 1] and r2['adq|default|bic'].tolist() == [600, 1201, 0] and r2['iso|ours|bic'].tolist() == [600, 1201, 0]
    assert 'iso|default|bic' in r1 and 'adq|default|bic' not in r1
    assert E.real_arrays(dict(src, rep=None), 'm2r2', (133, 281)) == {} and E.real_arrays(src, 'm2r2', (0, 0)) == {}
    assert E.alias_suffix('isotropic', 'adq|ours|k2') == 'iso|ours|k2' and E.alias_suffix('general', 'adq|ours|k2') == 'adq|ours|k2' and E.alias_suffix('isotropic', 'cs|default') == 'cs|default'
    assert E.alias_suffix('isotropic', 'iso|ours|k2') == 'iso|ours|k2' and E.alias_suffix('', 'adq|ours|k2') == 'adq|ours|k2'


def test_the_frozen_sources_are_read_and_the_repeat_raw_may_be_absent(tmp_path):
    a, b, c = str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz'), str(tmp_path / 'c.npz')
    np.savez_compressed(a, x=np.arange(3))
    np.savez_compressed(b, y=np.arange(2))
    src = E.real_sources(a, b, c)
    assert src['m2raw']['x'].tolist() == [0, 1, 2] and src['psf']['y'].tolist() == [0, 1] and src['rep'] is None
    np.savez_compressed(c, z=np.arange(4))
    assert E.real_sources(a, b, c)['rep']['z'].tolist() == [0, 1, 2, 3]


def test_part_real_stores_the_copy_and_the_hashes_of_its_sources(monkeypatch, tmp_path):
    src = planted_src()
    monkeypatch.setattr(E, 'real_sources', lambda *a: src)
    monkeypatch.setattr(E, 'SYNTH', [('m2', (310, 630)), ('m92', (1100, 220)), ('m92', (1430, 550)), ('m2r2', (133, 281))])
    srcs = []
    for nm in ('m2_raw.npz', 'psf_raw.npz', 'repeat_raw.npz'):
        (tmp_path / nm).write_bytes(nm.encode())
        srcs.append(str(tmp_path / nm))
    out = str(tmp_path / 'diag_raw.npz')
    E.part_real(out, *srcs)
    raw = X.load_raw(out)
    assert 'RF|m92|x1100_y220|adq|ours|k2|blk22' in raw and 'RF|m2|x310_y630|iso|ours|k2|blk22' in raw and 'RF|m2r2|x133_y281|adq|ours|k2|tot' in raw
    assert not any(k.startswith('RF|m92|x1430_y550|') for k in raw)                                    # the patch the frozen results do not hold is not copied
    import hashlib
    assert raw['RFsrc|files'].tolist() == ['%s:%s' % (n, hashlib.sha256(n.encode()).hexdigest()) for n in ('m2_raw.npz', 'psf_raw.npz', 'repeat_raw.npz')]


def test_the_provenance_hashes_the_code_the_runs_use():
    m = E.meta_arrays(dict(hp=3.0, rex=1.2, kap=2.0, tuned=dict(crowdsource='c', photutils='p')))
    assert m['cfg'].tolist() == [3.0, 1.2, 2.0] and m['tuned'].tolist() == ['c', 'p'] and len(m['files']) == len(E.CODE_FILES) == 6
    import hashlib
    f = [x for x in m['files'] if x.startswith('src/core/imaging.py:')][0]
    assert f.split(':')[1] == hashlib.sha256(open(os.path.join(E.REPO, 'src/core/imaging.py'), 'rb').read()).hexdigest()


# ----------------------------------------------------------------------------------------------------------------------------------
# the report on a planted raw
# ----------------------------------------------------------------------------------------------------------------------------------
REAL_F = {'iso|ours|k2': 2, 'adq|ours|k2': 3, 'cs|default': 4, 'iso|default|k5': 2.5, 'adq|default|k5': 3.5, 'phot|tuned': 1}          # F of the planted blocks in units of 0.1
SYNTH_F = {'iso|ours|k2': 5, 'adq|ours|k2': 6.5, 'cs|default': 4.5, 'iso|default|k5': 5.5, 'adq|default|k5': 7.5, 'phot|tuned': 1.5}
BUMP = {(1430, 550): 1.0}                                                                                                                      # the second M92 patch is 0.1 higher in every column
RUNS = ('iso|ours', 'iso|default', 'adq|ours', 'adq|default')


def tot_for(f10):
    """Planted totals [true22, matched22, det22, det22_matched, true_all, matched_all, det_all] of 25 blocks: completeness f10 / 10, purity f10 / 8 (they differ, so that a table cannot swap them unseen)."""
    return np.array([250, 25 * f10, 200, 25 * f10, 400, 30 * f10, 500])


def cell(v10, same=False, stalled=False):
    v = v10 / 10.0
    return '%.3f [%.3f, %.3f]%s%s' % (v, v, v, ' (same)' if same else '', ' STALLED' if stalled else '')


def planted_raw(stalled=()):
    """A raw as the parts write it, with distinguishable values (REAL_F, SYNTH_F; every block of a patch equal, so that a stratified bootstrap has no spread and an unstratified one has). `stalled`: the runs
    whose 'bic' record has reached the entry cap, as 'W|<set>|x<x0>_y<y0>|<kernel>|<configuration>' (a world) or 'RF|...' (the frozen results)."""
    def blk(f10):
        b = np.zeros((25, 4))
        b[:, 0], b[:, 1], b[:, 2], b[:, 3] = 10, f10, 10, f10
        return b

    def bic(world, tag, run, ok):
        return np.array([ok[0], X.MAX_ENTRIES + 1, 1]) if '%s|%s|%s' % (world, tag, run) in stalled else np.array(ok)
    raw = {}
    patches = [('m92', (1100, 220)), ('m92', (1430, 550)), ('m2', (310, 630))]
    real_f, synth_f = REAL_F, SYNTH_F
    for name, pxy in patches:
        kind = 'isotropic' if name == 'm2' else 'general'
        tag = '%s|x%d_y%d' % (name, pxy[0], pxy[1])
        raw['WH|%s|kind' % tag] = np.array([kind])
        raw['WH|%s|stats' % tag] = np.array([[300.0, 250.0, 450.0, 3000.0], [290.0, 245.0, 440.0, 2800.0]])
        raw['WH|%s|counts' % tag] = np.array([6000, 3141, 600])
        raw['WH|%s|fwhm' % tag] = np.array([3.0])
        raw['RF|%s|kind' % tag] = np.array([kind])
        for suffix in real_f:
            if kind == 'isotropic' and suffix.startswith('adq'):
                continue
            raw['RF|%s|%s|blk22' % (tag, suffix)] = blk(real_f[suffix] + BUMP.get(pxy, 0.0))
            raw['W|%s|%s|blk22' % (tag, suffix)] = blk(synth_f[suffix] + BUMP.get(pxy, 0.0))
            raw['RF|%s|%s|tot' % (tag, suffix)] = tot_for(real_f[suffix] + BUMP.get(pxy, 0.0))
            raw['W|%s|%s|tot' % (tag, suffix)] = tot_for(synth_f[suffix] + BUMP.get(pxy, 0.0))
        raw['W|%s|cs|secs' % tag] = np.array([5.0, 900])
        for run in RUNS:
            if kind == 'isotropic' and run.startswith('adq'):
                continue
            raw['W|%s|%s|secs' % (tag, run)] = np.array([200.0 if run.startswith('adq') else 100.0, 1200 if run.startswith('adq') else 900])
            raw['W|%s|%s|bic' % (tag, run)] = bic('W', tag, run, [1100, 1201, 0] if run.startswith('adq') else [800, 901, 0])
            raw['RF|%s|%s|bic' % (tag, run)] = bic('RF', tag, run, [700, 1001, 0])
        for lab, (rv, sv) in (('ours', (1.84, 2.19)), ('default', (2.5, 3.0))):
            raw['RF|%s|iso|%s|rsd' % (tag, lab)], raw['W|%s|iso|%s|rsd' % (tag, lab)] = np.array([rv]), np.array([sv])
            if kind == 'general':
                raw['RF|%s|adq|%s|rsd' % (tag, lab)], raw['W|%s|adq|%s|rsd' % (tag, lab)] = np.array([rv - 0.6]), np.array([sv - 1.0])
    raw['META|cfg'], raw['META|tuned'], raw['META|files'] = np.array([3.0, 1.2, 2.0]), np.array(['crowdsource|snr5|ff0|qf0|spinf', 'photutils_ext_all_local_i6']), np.array(['src/core/imaging.py:' + 'ab' * 32])
    raw['RFsrc|files'] = np.array(['m2_raw.npz:' + 'cd' * 32])
    # item 1: two patches of M92 and the pooled stack, stars broader than the stamp by 11.5 %
    rng = np.random.default_rng(3)
    ref = E.stamp_core(gauss_stamp(1.3))
    for pxy in ((1100, 220), (1430, 550)):
        real = np.array([E.aligned(gauss_stamp(1.45)[13:38, 13:38] + 0.0005 * rng.normal(size=(25, 25))) for _ in range(9)])
        syn = np.array([E.aligned(gauss_stamp(1.3)[13:38, 13:38] + 0.0005 * rng.normal(size=(25, 25))) for _ in range(9)])
        for k, v in E.stack_arrays(real, syn, ref, [9], seed=1, B=40).items():
            raw['S|m92|x%d_y%d|%s' % (pxy + (k,))] = v
    for k, v in E.stack_arrays(np.concatenate([real, real]), np.concatenate([syn, syn]), ref, [9, 9], seed=2, B=40).items():
        raw['S|m92|all|%s' % k] = v
    # item 3: six bright stars per patch; the M2 test patches other than the world's (real only)
    patches_b = patches + [('m2', (310, 730)), ('m2', (360, 1040)), ('m2', (100, 1100))]
    for world in ('real', 'synth'):
        for name, pxy in patches_b:
            if world == 'synth' and (name, pxy) not in patches:
                continue
            tag = '%s|x%d_y%d' % (name, pxy[0], pxy[1])
            raw['B|real|%s|kind' % tag] = np.array(['isotropic' if name == 'm2' else 'general'])
            nsig = {'real': [2, 2, 1, 1, 0, 2], 'synth': [1, 1, 1, 2, 1, 1]}[world] if (name, pxy) in patches else [0, 0, 1, 1, 2, 2, 2, 2]
            for s_ in ('iso|ours|k2', 'adq|ours|k2', 'cs|default', 'cs|tuned'):
                ns_ = [v + 1 for v in nsig] if (s_ == 'cs|tuned' and name == 'm2') else nsig           # the frozen post-filter differs from the default on M2 only
                k = 'B|%s|%s|%s|' % (world, tag, s_)
                raw[k + 'nsig'] = np.array(ns_)
                raw[k + 'nunits'] = np.array(ns_)
                raw[k + 'fratio'] = np.array([1.0, 1.1, 0.9, 1.0, 0.0, 1.2][:len(nsig)] + [1.0] * (len(nsig) - 6))
                raw[k + 'pairs'] = np.array([[0.2, -0.1, 0.3], [0.3, -0.2, 0.4], [0.1, 0.0, 0.2]]) if world == 'real' else np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.1]])
                run = E.path_tag(E.alias_suffix(raw['B|real|%s|kind' % tag][0], s_))
                raw[k + 'stalled'] = np.array([int(run is not None and '%s|%s|%s' % ('RF' if world == 'real' else 'W', tag, run) in stalled)])
    return raw


MAIN = ('iso|ours|k2', 'adq|ours|k2', 'cs|default')
DFLT = ('iso|default|k5', 'adq|default|k5', 'cs|default')
PHOT = ('iso|ours|k2', 'phot|tuned', 'cs|default')


def row_of(key, bump, same_m2=False):
    """The six F cells of a table row: the real columns then the synthetic ones."""
    return ' | '.join([cell(REAL_F[k] + bump, same_m2 and k.startswith('adq')) for k in key] + [cell(SYNTH_F[k] + bump, same_m2 and k.startswith('adq')) for k in key])


def test_the_report_reads_the_right_worlds_columns_and_keys_and_is_deterministic():
    raw = planted_raw()
    md = E.report(raw)
    assert md == E.report(raw) and md.endswith('\n') and 'nan' not in md.lower().replace('nan,', '')
    sec2 = md.split('## 2.')[1].split('## 3.')[0]
    for needle in ('POST HOC DIAGNOSTICS, NOT TESTS', '## 1. The stamp against the stars', '### By patch', '### Each set pooled', '## 2. Exact-stamp synthetic worlds', '## 3. Splitting and flux bias of bright isolated stars',
                   '## Limitations', 'Configuration read from `freeze.json` (checked by `exp_astro_m2.check_freeze` before the runs): sigma_b 3, r_ex 1.2, kappa 2; the a priori default (5.0, 1.2, 5.0)', 'src/core/imaging.py abababababab', 'm2_raw.npz cdcdcdcdcdcd'):
        assert needle in md, needle
    for needle in ('| set | patch | kernel chosen | real iso | real adequate | real crowdsource | synthetic iso | synthetic adequate | synthetic crowdsource |',
                   '| m92 | (1100, 220) | general | ' + row_of(MAIN, 0.0) + ' |', '| m92 | (1430, 550) | general | ' + row_of(MAIN, 1.0) + ' |',
                   '| m2 | (310, 630) | isotropic | ' + ' | '.join([cell(2), cell(2, True), cell(4), cell(5), cell(5, True), cell(4.5)]) + ' |',
                   '| M92 (two windows of one frame) | 2 of 2 | ' + row_of(MAIN, 0.5) + ' | - |',                                           # the pooled F of equal patches is the mean of the two
                   '| M2 (310, 630) | 1 of 1 | ' + ' | '.join([cell(2), cell(2, True), cell(4), cell(5), cell(5, True), cell(4.5)]) + ' | - |',
                   '| pool | patches | real iso | real adequate | real crowdsource | synthetic iso | synthetic adequate | synthetic crowdsource | left out (stalled path) |',
                   '| M92 (two windows of one frame) | synthetic iso - real iso | 0.300 [0.300, 0.300] | synthetic higher | 2 of 2 |', '| M92 (two windows of one frame) | synthetic adequate - real adequate | 0.350 [0.350, 0.350] | synthetic higher | 2 of 2 |',
                   '| M92 (two windows of one frame) | synthetic crowdsource - real crowdsource | 0.050 [0.050, 0.050] | synthetic higher | 2 of 2 |', '| M92 (two windows of one frame) | real iso - real crowdsource | -0.200 [-0.200, -0.200] | loss | 2 of 2 |',
                   '| M92 (two windows of one frame) | synthetic iso - synthetic crowdsource | 0.050 [0.050, 0.050] | win | 2 of 2 |', '| M92 (two windows of one frame) | real adequate - real crowdsource | -0.100 [-0.100, -0.100] | loss | 2 of 2 |',
                   '| M92 (two windows of one frame) | synthetic adequate - synthetic crowdsource | 0.200 [0.200, 0.200] | win | 2 of 2 |',
                   '| M92 (two windows of one frame) | lead of iso over crowdsource, real minus synthetic | -0.250 [-0.250, -0.250] | real lead smaller | 2 of 2 |',
                   '| M92 (two windows of one frame) | lead of adequate over crowdsource, real minus synthetic | -0.300 [-0.300, -0.300] | real lead smaller | 2 of 2 |',
                   '| M2 (310, 630) | lead of iso over crowdsource, real minus synthetic | -0.250 [-0.250, -0.250] | real lead smaller | 1 of 1 |'):
        assert needle in sec2, needle
    # the completeness / purity / unmatched table: the totals of the pooled patches
    def cpu(fs):
        t = sum(tot_for(f) for f in fs)
        return '%.3f / %.3f / %d' % (t[1] / t[0], t[3] / t[2], t[6] - t[5])
    cpt = sec2.split('### Completeness / purity')[1].split('### ')[0]
    assert '| M92 (two windows of one frame) | 2 of 2 | ' + ' | '.join(cpu([REAL_F[k], REAL_F[k] + 1.0]) for k in MAIN) + ' | ' + ' | '.join(cpu([SYNTH_F[k], SYNTH_F[k] + 1.0]) for k in MAIN) + ' | - |' in cpt
    assert '| M2 (310, 630) | 1 of 1 | ' + ' | '.join(cpu([REAL_F[k if k != 'adq|ours|k2' else 'iso|ours|k2']]) for k in MAIN) + ' | ' + ' | '.join(cpu([SYNTH_F[k if k != 'adq|ours|k2' else 'iso|ours|k2']]) for k in MAIN) + ' | - |' in cpt
    assert '| M2 (310, 630) | synthetic adequate - real adequate |' not in md and 'M2 (310, 630) | lead of adequate' not in md                  # an isotropic-kernel patch has no adequate comparison of its own
    # the a priori default and photutils: tables of their own, the right keys in the right worlds
    assert '| m92 | (1100, 220) | general | ' + row_of(DFLT, 0.0) + ' |' in sec2 and '| M92 (two windows of one frame) | 2 of 2 | ' + row_of(DFLT, 0.5) + ' | - |' in sec2
    assert '| m92 | (1100, 220) | general | ' + row_of(PHOT, 0.0) + ' |' in sec2 and '| M92 (two windows of one frame) | 2 of 2 | ' + row_of(PHOT, 0.5) + ' | - |' in sec2
    assert 'F at mag < 22, a priori default (sigma_b 5, r_ex 1.2, kappa 5): pooled' in sec2 and 'photutils (the variant tuned on the development patches) beside the ISS and crowdsource: pooled' in sec2
    # the residual-scale table pairs the real and the synthetic value of each kernel and shows the isotropic one for the adequate kernel of an isotropic patch
    assert '| m92 | (1100, 220) | general | 1.84 | 2.19 | 1.24 | 1.19 | 2.50 | 3.00 | 1.90 | 2.00 |' in sec2 and '| m2 | (310, 630) | isotropic | 1.84 | 2.19 | 1.84 | 2.19 | 2.50 | 3.00 | 2.50 | 3.00 |' in sec2
    assert '| m92 | (1100, 220) | 6000 | 3141 | 600 | 300.0 / 290.0 | 250 - 450 / 245 - 440 | 3000 / 2800 | 3.00 |' in sec2
    assert '| m92 | (1100, 220) | ISS, chosen kernel, tuned (kappa 2) | 1200 | 200 | 1100 of 1201 | - |' in sec2 and '| m92 | (1100, 220) | crowdsource (default) | 900 | 5 | - | - |' in sec2
    assert '| m92 | (1100, 220) | ISS, isotropic kernel, tuned (kappa 2) | 900 | 100 | 800 of 901 | - |' in sec2 and '| m92 | (1100, 220) | ISS, isotropic kernel, a priori default (kappa 5) | 900 | 100 | 800 of 901 | - |' in sec2
    # nothing stalled: the section says so, with the number of path runs of the worlds (10) and of the frozen results copied (10)
    assert 'Synthetic worlds: none of the 10 ISS path runs reached the cap.' in sec2 and 'Real data (the frozen results of the same patches): none of the 10 ISS path runs copied reached the cap.' in sec2
    assert 'STALLED' not in sec2.replace('marked STALLED', '').replace('STALLED: the path reached', '').replace('stalled path', '')
    # section 1: the stamp rows (stars 9 + 9 = 18; stars 11.5 % broader than the stamp: central pixel about -20 %, FWHM about +0.35 px)
    sec1 = md.split('## 1.')[1].split('## 2.')[0]
    assert '| m92 | (1100, 220) | 9 | 3.07 | 3.4' in sec1 and '| m92 | all | 18 | 3.07 | 3.4' in sec1 and '| m2 | ' not in sec1.split('### Each set pooled')[1]
    row = [l for l in sec1.split('\n') if l.startswith('| m92 | all |')][0].split(' | ')
    assert row[6].startswith('+0.3') and row[7].startswith('-19.') and row[8].startswith('-0.') and row[9].startswith('-19.') and row[10].startswith('-19.'), row
    # section 3: real splits 3 of 6 = 50 %, synthetic 1 of 6 = 17 %, no unit 17 % and 0 %; the median dmag of the planted pairs
    sec3 = md.split('## 3.')[1].split('## Limitations')[0]
    assert '| ISS, isotropic kernel, tuned (kappa 2) | real, the patches of the worlds | 12 | 50 [25, 75] | 17 [0, 42] | 1.00 [0.90, 1.10] | 6 | +0.30 [+0.20, +0.40] | +0.20 [+0.10, +0.30] | -0.10 [-0.20, +0.00] | - |' in sec3
    assert '| ISS, isotropic kernel, tuned (kappa 2) | synthetic | 12 | 17 [0, 42] | 0 [0, 0] | 1.00 [0.90, 1.10] | 4 | +0.05 [+0.00, +0.10] | +0.05 [+0.00, +0.10] | +0.00 [+0.00, +0.00] | - |' in sec3
    # the M2 table also has the real rows of all four test patches (6 + 3 x 8 stars), the world's patch being one of them
    m2tab = sec3.split('### m2:')[1].split('### ')[0]
    assert '| ISS, isotropic kernel, tuned (kappa 2) | real, all 4 patches | 30 | ' in m2tab and '| ISS, isotropic kernel, tuned (kappa 2) | real, the patches of the worlds | 6 | ' in m2tab
    assert '| ISS, isotropic kernel, tuned (kappa 2) | synthetic | 6 | ' in m2tab
    # an M2 set whose rule keeps the isotropic kernel has no rows of the chosen kernel; the post-filter's rows are printed where they differ from the default's (M2) and not where they do not (M92)
    m92tab = sec3.split('### m92:')[1].split('### ')[0]
    assert 'chosen kernel' not in m2tab and 'chosen kernel' in m92tab and 'crowdsource (frozen post-filter)' in m2tab and 'crowdsource (frozen post-filter)' not in m92tab
    assert 'real, all' not in m92tab and 'real, all 4 patches' in m2tab                                                     # a set whose patches all have a world has one real row, M2 has the extra one for its four test patches
    assert '| crowdsource (frozen post-filter) | synthetic | 6 | ' in m2tab


def test_each_statistic_of_a_stamp_row_is_printed_in_its_own_column():
    ref = E.stamp_core(gauss_stamp(1.3))
    est = np.array([-0.11, -0.12, 0.013, 0.021, 3.31, 3.12, 0.19, 0.94, -0.071])         # real / stamp - 1, synthetic / stamp - 1, real minus synthetic, ratio - 1, FWHM real, synthetic, difference, ratio, paired
    ci = np.stack([est - 0.01 * np.arange(1, 10), est + 0.02 * np.arange(1, 10)], axis=1)
    raw = {'S|m92|all|est': est, 'S|m92|all|ci': ci, 'S|m92|all|n': np.array([18]), 'S|m92|all|ref': ref}
    row = [l for l in '\n'.join(E.stamp_section(raw)).split('\n') if l.startswith('| m92 | all |')][0]
    cells = [c.strip() for c in row.strip().strip('|').split(' | ')]
    assert cells == ['m92', 'all', '18', '%.2f' % E.fwhm(ref), E.num(est[4], ci[4], sign=False), E.num(est[5], ci[5], sign=False), E.num(est[6], ci[6]),
                     E.pci(est[0], ci[0]), E.pci(est[1], ci[1]), E.pci(est[2], ci[2]), E.pci(est[3], ci[3]), E.pci(est[8], ci[8])], cells
    assert len(set(cells[7:])) == 5                                                       # the five statistics print differently here: a column cannot be read from another


STALLED_RUNS = ('W|m92|x1430_y550|adq|ours', 'RF|m92|x1100_y220|iso|default', 'W|m2|x310_y630|iso|default')


def test_stalled_paths_are_flagged_everywhere_left_out_of_every_pool_and_named():
    raw = planted_raw(stalled=STALLED_RUNS)
    md = E.report(raw)
    sec2 = md.split('## 2.')[1].split('## 3.')[0]
    # the section lists the defective runs of the worlds and of the frozen results
    assert 'Synthetic worlds: 2 of the 10 ISS path runs reached the cap and are defective: m2 (310, 630) ISS, isotropic kernel, a priori default (kappa 5); m92 (1430, 550) ISS, chosen kernel, tuned (kappa 2).' in sec2
    assert 'Real data (the frozen results of the same patches): 1 of the 10 ISS path runs copied reached the cap and are defective: m92 (1100, 220) ISS, isotropic kernel, a priori default (kappa 5).' in sec2
    assert 'defective run, not a censored valid one' in sec2
    # the tuned table: the stalled synthetic adequate run of M92 (1430, 550) is marked, the pool takes the other patch only and names the one left out; the pools of the other sets are untouched
    assert '| m92 | (1430, 550) | general | ' + ' | '.join([cell(REAL_F[k] + 1.0) for k in MAIN] + [cell(SYNTH_F['iso|ours|k2'] + 1.0), cell(SYNTH_F['adq|ours|k2'] + 1.0, stalled=True), cell(SYNTH_F['cs|default'] + 1.0)]) + ' |' in sec2
    assert '| m92 | (1100, 220) | general | ' + row_of(MAIN, 0.0) + ' |' in sec2
    assert '| M92 (two windows of one frame) | 1 of 2 | ' + row_of(MAIN, 0.0) + ' | m92 (1430, 550): synthetic adequate |' in sec2
    assert '| M2 (310, 630) | 1 of 1 | ' + ' | '.join([cell(2), cell(2, True), cell(4), cell(5), cell(5, True), cell(4.5)]) + ' | - |' in sec2
    # the completeness / purity table and the differences leave the patch out as well and say so
    def cpu(fs):
        t = sum(tot_for(f) for f in fs)
        return '%.3f / %.3f / %d' % (t[1] / t[0], t[3] / t[2], t[6] - t[5])
    cpt = sec2.split('### Completeness / purity')[1].split('### ')[0]
    assert '| M92 (two windows of one frame) | 1 of 2 | ' + ' | '.join(cpu([REAL_F[k]]) for k in MAIN) + ' | ' + ' | '.join(cpu([SYNTH_F[k]]) for k in MAIN) + ' | m92 (1430, 550): synthetic adequate |' in cpt
    dif = sec2.split('### Paired differences')[1].split('### ')[0]
    assert '| M92 (two windows of one frame) | synthetic iso - real iso | 0.300 [0.300, 0.300] | synthetic higher | 1 of 2; left out (stalled path): m92 (1430, 550): synthetic adequate |' in dif
    assert '| M92 (two windows of one frame) | lead of adequate over crowdsource, real minus synthetic | -0.300 [-0.300, -0.300] | real lead smaller | 1 of 2; left out (stalled path): m92 (1430, 550): synthetic adequate |' in dif
    assert '| M2 (310, 630) | synthetic iso - real iso | 0.300 [0.300, 0.300] | synthetic higher | 1 of 1 |' in dif
    # the a priori default: the real isotropic run of M92 (1100, 220) and the synthetic isotropic run of M2 (and, for the isotropic kernel, its adequate alias) are marked; a pool with no valid patch is a row of dashes
    dft = sec2.split('a priori default (sigma_b 5, r_ex 1.2, kappa 5): by patch')[1]
    assert 'STALLED: the path reached the entry cap (a defective run, see above); its F is not a valid row.' in dft.split('| set |')[0] and 'censored' not in dft.split('| set |')[0]               # the note says defective, not censored
    assert '| m92 | (1100, 220) | general | ' + ' | '.join([cell(REAL_F['iso|default|k5'], stalled=True), cell(REAL_F['adq|default|k5']), cell(REAL_F['cs|default'])] + [cell(SYNTH_F[k]) for k in DFLT]) + ' |' in dft
    assert '| m2 | (310, 630) | isotropic | ' + ' | '.join([cell(REAL_F['iso|default|k5']), cell(REAL_F['iso|default|k5'], True), cell(REAL_F['cs|default']), cell(SYNTH_F['iso|default|k5'], stalled=True),
                                                         cell(SYNTH_F['iso|default|k5'], True, True), cell(SYNTH_F['cs|default'])]) + ' |' in dft
    assert '| M92 (two windows of one frame) | 1 of 2 | ' + ' | '.join([cell(REAL_F[k] + 1.0) for k in DFLT] + [cell(SYNTH_F[k] + 1.0) for k in DFLT]) + ' | m92 (1100, 220): real iso |' in dft
    assert '| M2 (310, 630) | 0 of 1 | - | - | - | - | - | - | m2 (310, 630): synthetic iso, synthetic adequate |' in dft
    # the residual scale and the run records carry the mark
    assert '| m92 | (1100, 220) | general | 1.84 | 2.19 | 1.24 | 1.19 | 2.50 STALLED | 3.00 | 1.90 | 2.00 |' in sec2
    assert '| m92 | (1430, 550) | general | 1.84 | 2.19 | 1.24 | 1.19 STALLED | 2.50 | 3.00 | 1.90 | 2.00 |' in sec2
    assert '| m2 | (310, 630) | isotropic | 1.84 | 2.19 | 1.84 | 2.19 | 2.50 | 3.00 STALLED | 2.50 | 3.00 STALLED |' in sec2
    assert '| m92 | (1430, 550) | ISS, chosen kernel, tuned (kappa 2) | 1200 | 200 | 1100 of 2501 | STALLED: the entry cap (2500) was reached, a defective run |' in sec2
    assert '| m2 | (310, 630) | ISS, isotropic kernel, a priori default (kappa 5) | 900 | 100 | 800 of 2501 | STALLED: the entry cap (2500) was reached, a defective run |' in sec2
    assert '| m92 | (1100, 220) | ISS, chosen kernel, tuned (kappa 2) | 1200 | 200 | 1100 of 1201 | - |' in sec2
    # section 3: the stalled synthetic adequate catalogue of M92 (1430, 550) is left out of the pool (6 stars of the other patch) and named; the real rows are untouched
    sec3 = md.split('## 3.')[1].split('## Limitations')[0]
    m92tab = sec3.split('### m92:')[1].split('### ')[0]
    assert '| ISS, chosen kernel, tuned (kappa 2) | synthetic | 6 | 17 [0, 50] | 0 [0, 0] | 1.00 [0.45, 1.15] | 2 | +0.05 [+0.00, +0.10] | +0.05 [+0.00, +0.10] | +0.00 [+0.00, +0.00] | (1430, 550) |' in m92tab
    assert '| ISS, chosen kernel, tuned (kappa 2) | real, the patches of the worlds | 12 | 50 [25, 75] | 17 [0, 42] | 1.00 [0.90, 1.10] | 6 | +0.30 [+0.20, +0.40] | +0.20 [+0.10, +0.30] | -0.10 [-0.20, +0.00] | - |' in m92tab
    assert '| ISS, isotropic kernel, tuned (kappa 2) | synthetic | 12 | 17 [0, 42] |' in m92tab
    # every path stalled: a row with no statistic, the patches named
    for pxy in ((1100, 220), (1430, 550)):
        raw['B|synth|m92|x%d_y%d|adq|ours|k2|stalled' % pxy] = np.array([1])
    m92tab = E.report(raw).split('## 3.')[1].split('### m92:')[1].split('### ')[0]
    assert '| ISS, chosen kernel, tuned (kappa 2) | synthetic | STALLED: every path stalled | - | - | - | - | - | - | - | (1100, 220); (1430, 550) |' in m92tab
    # two patches left out of one pool are separated by semicolons, the columns of one patch by commas
    md2 = E.report(planted_raw(stalled=STALLED_RUNS + ('W|m92|x1100_y220|iso|ours', 'W|m92|x1100_y220|iso|default')))
    sec2b = md2.split('## 2.')[1].split('## 3.')[0]
    assert '| M92 (two windows of one frame) | 0 of 2 | - | - | - | - | - | - | m92 (1100, 220): synthetic iso; m92 (1430, 550): synthetic adequate |' in sec2b
    left = 'm92 (1100, 220): synthetic iso; m92 (1430, 550): synthetic adequate'
    assert '| M92 (two windows of one frame) | 0 of 2 | - | - | - | - | - | - | ' + left + ' |' in sec2b.split('### Completeness / purity')[1].split('### ')[0]                      # a pool with no valid patch is a row of dashes
    assert '| M92 (two windows of one frame) | - | - | - | 0 of 2; left out (stalled path): ' + left + ' |' in sec2b.split('### Paired differences')[1].split('### ')[0]
    assert E.excluded_text([(('m92', (1100, 220)), ['real iso', 'synthetic iso']), (('m92', (1430, 550)), ['synthetic adequate'])]) == 'm92 (1100, 220): real iso, synthetic iso; m92 (1430, 550): synthetic adequate' and E.excluded_text([]) == '-'
    # the limitation is printed
    assert 'A path that reaches the entry cap has stalled' in md.split('## Limitations')[1]


def test_the_report_of_an_empty_raw_says_what_is_missing_and_a_partial_raw_has_dashes():
    md = E.report({})
    assert 'No world has been run.' in md and '## 3.' in md and md == E.report({})
    raw = planted_raw()
    for k in [k for k in raw if k.startswith('W|m92|x1430_y550|cs|')]:
        del raw[k]
    md = E.report(raw)
    row = [l for l in md.split('\n') if l.startswith('| m92 | (1430, 550) | general |')][0]
    assert row.split(' | ')[-1].startswith('-') and '0.600 [0.600, 0.600]' in row
    assert '| M92 (two windows of one frame) | 1 of 2 |' in md                                      # a pool counts the patches that have every column
    assert 'F at mag < 22, a priori default (sigma_b 5, r_ex 1.2, kappa 5): by patch' in md
    for k in [k for k in raw if k.startswith('W|') and '|default|' in k and 'cs|' not in k]:                   # no a priori default run in any world: that table is not printed
        del raw[k]
    md = E.report(raw)
    assert 'F at mag < 22, a priori default' not in md and 'F at mag < 22, tuned configuration (kappa 2): by patch' in md
    for k in [k for k in raw if k.startswith('W|') and k.endswith('|blk22')]:
        del raw[k]
    md = E.report(raw)
    assert 'No world has been run.' not in md                                                      # the runs' records are still there
    assert '### F at mag < 22, tuned configuration (kappa 2): by patch' not in md                  # but there is no synthetic F to tabulate
    # the runs of the worlds are counted where they exist, and so are those of the frozen results
    raw = planted_raw()
    del raw['W|m92|x1100_y220|iso|default|bic']
    sec2 = E.report(raw).split('## 2.')[1].split('## 3.')[0]
    assert 'Synthetic worlds: none of the 9 ISS path runs reached the cap.' in sec2 and 'Real data (the frozen results of the same patches): none of the 10 ISS path runs copied reached the cap.' in sec2
    # a patch whose world has no bright-star catalogues is not one of the worlds' patches: the real rows of the set then also show all of its patches
    raw = planted_raw()
    for k in [k for k in raw if k.startswith('B|synth|m92|x1430_y550|')]:
        del raw[k]
    m92tab = E.report(raw).split('## 3.')[1].split('### m92:')[1].split('### ')[0]
    assert ('| ISS, isotropic kernel, tuned (kappa 2) | real, all 2 patches | 12 |' in m92tab and '| ISS, isotropic kernel, tuned (kappa 2) | real, the patches of the worlds | 6 |' in m92tab
            and '| ISS, isotropic kernel, tuned (kappa 2) | synthetic | 6 |' in m92tab)


def test_bootstrap_intervals_are_the_2_5_and_97_5_percentiles_of_the_stratified_resamples():
    a, b = np.arange(10.0), np.arange(10.0, 25.0)
    idx = X.boot_idx_strat([10, 15], 200, 5)
    v = np.concatenate([a, b])
    f, ci = E.boot_fraction([a, b], lambda u: u > 12, B=200, seed=5)
    assert f == (v > 12).mean() and np.allclose(ci, np.percentile((v > 12)[idx].mean(1), [2.5, 97.5]))
    m, ci = E.boot_median([a, b], B=200, seed=5)
    assert m == np.median(v) and np.allclose(ci, np.percentile(np.median(v[idx], axis=1), [2.5, 97.5]))


def test_a_by_patch_and_a_pooled_row_use_the_blocks_of_the_patch_with_the_repositorys_bootstrap():
    rng = np.random.default_rng(8)
    raw = {}
    for pre in ('RF', 'W'):
        for s_ in ('iso|ours|k2', 'cs|default'):
            b = np.zeros((25, 4))
            b[:, 0] = rng.integers(8, 14, 25)
            b[:, 1] = rng.integers(2, 8, 25)
            b[:, 2] = rng.integers(8, 14, 25)
            b[:, 3] = rng.integers(2, 8, 25)
            raw['%s|m2|x310_y630|%s|blk22' % (pre, s_)] = b
    raw['WH|m2|x310_y630|kind'] = np.array(['isotropic'])
    md = E.report(raw)
    idx = X.boot_idx(25)
    cells = [X.fmt_iv(*X.boot_cpf(raw['%s|m2|x310_y630|%s|blk22' % (pre, s_)], idx)[2]) for pre, s_ in (('RF', 'iso|ours|k2'), ('RF', 'iso|ours|k2'), ('RF', 'cs|default'), ('W', 'iso|ours|k2'), ('W', 'iso|ours|k2'), ('W', 'cs|default'))]
    cells[1] += ' (same)'
    cells[4] += ' (same)'
    assert '| m2 | (310, 630) | isotropic | ' + ' | '.join(cells) + ' |' in md
    idx = X.boot_idx_strat([25])
    pooled = [X.fmt_iv(*X.boot_cpf(raw['%s|m2|x310_y630|%s|blk22' % (pre, s_)], idx)[2]) for pre, s_ in (('RF', 'iso|ours|k2'), ('RF', 'iso|ours|k2'), ('RF', 'cs|default'), ('W', 'iso|ours|k2'), ('W', 'iso|ours|k2'), ('W', 'cs|default'))]
    pooled[1] += ' (same)'
    pooled[4] += ' (same)'
    assert '| M2 (310, 630) | 1 of 1 | ' + ' | '.join(pooled) + ' |' in md
    # the paired difference of the two worlds on the same resamples
    b_s, b_r = raw['W|m2|x310_y630|iso|ours|k2|blk22'], raw['RF|m2|x310_y630|iso|ours|k2|blk22']
    df = X.boot_diff(b_s, b_r, idx)[2]
    assert '| M2 (310, 630) | synthetic iso - real iso | ' + X.fmt_iv(df[0], df[1]) + ' | ' in md
    # the change of the lead of the ISS over crowdsource between the real data and the world: the two leads on the same resamples of the blocks
    blk = lambda pre, s_: raw['%s|m2|x310_y630|%s|blk22' % (pre, s_)]
    def lead(pre, ix=None):                                                       # F of the ISS minus F of crowdsource: of the blocks as they are, or of each resample
        f = lambda b: X.cpf(b.sum(0) if ix is None else b[ix].sum(1))[2]
        return f(blk(pre, 'iso|ours|k2')) - f(blk(pre, 'cs|default'))
    point = lead('RF') - lead('W')
    ci = np.percentile(lead('RF', idx) - lead('W', idx), [2.5, 97.5])
    assert '| M2 (310, 630) | lead of iso over crowdsource, real minus synthetic | ' + X.fmt_iv(point, ci) + ' | ' in md


def test_part_stamps_stores_each_patch_and_the_set_pooled_with_the_stars_weighted_stamp(monkeypatch, tmp_path):
    truths, patches = {}, {}
    for pxy, side, ssig, keep in (((0, 0), 4, 1.2, 16), ((5, 5), 3, 1.5, 9), ((7, 7), 1, 1.3, 2), ((9, 9), 2, 1.4, 3)):
        xs = 20.0 + 17.0 * np.arange(side)
        rows = [(x + 0.31 * ((i + j) % 3) - 0.3, y + 0.27 * ((i * j) % 4) - 0.4, A.counts_to_mag(3.0e4)) for i, x in enumerate(xs) for j, y in enumerate(xs)]
        truths[pxy] = planted_truth(rows[:keep], *pxy)
        real = E.render_world(truths[pxy], pxy[0], pxy[1], gauss_stamp(1.1 * ssig), 4.6, 1.0, np.full((100, 100), 150.0), 99)['dn']
        patches[pxy] = (dict(dn=real, psf=gauss_stamp(ssig), gain=4.6, dark=1.0), np.full((100, 100), 150.0), truths[pxy], 5.5)
    monkeypatch.setattr(E, 'STAMP_SETS', ('m92',))
    monkeypatch.setattr(E, 'set_patches', lambda name: [(0, 0), (5, 5), (7, 7), (9, 9)])
    monkeypatch.setattr(E, 'load_inputs', lambda name, pxy: patches[pxy])
    seen = dict(seeds=[], bsize=[], worlds=[])
    orig_arrays, orig_world = E.stack_arrays, E.render_world
    monkeypatch.setattr(E, 'stack_arrays', lambda real, syn, ref, sizes, seed, B=None: (seen['seeds'].append(seed), seen['bsize'].append(B), orig_arrays(real, syn, ref, sizes, seed, B))[2])
    monkeypatch.setattr(E, 'render_world', lambda truth, x0, y0, stamp, gain, dark, sky, seed, **kw: (seen['worlds'].append(seed), orig_world(truth, x0, y0, stamp, gain, dark, sky, seed, **kw))[1])
    out = str(tmp_path / 'diag_raw.npz')
    E.part_stamps(out, B=20)
    raw = X.load_raw(out)
    assert seen['seeds'] == [0, 1, 3, 99] and seen['bsize'] == [20] * 4                                                                       # one seed per patch (its place in the set), another for the pool
    assert seen['worlds'] == [E.world_seed('m92', p) for p in ((0, 0), (5, 5), (7, 7), (9, 9))]
    assert raw['S|m92|x0_y0|n'].tolist() == [16] and raw['S|m92|x5_y5|n'].tolist() == [9] and 'S|m92|x7_y7|n' not in raw and raw['S|m92|x9_y9|n'].tolist() == [3]                # three stars are enough
    assert raw['S|m92|all|n'].tolist() == [28] and raw['S|m92|all|sizes'].tolist() == [16, 9, 3]
    for k in ('x0_y0', 'x5_y5', 'x9_y9', 'all'):
        assert raw['S|m92|%s|est' % k].shape == (9,) and raw['S|m92|%s|ci' % k].shape == (9, 2) and -0.21 < raw['S|m92|%s|est' % k][8] < -0.14, k             # stars 10 % broader than the stamp
    ref0, ref1, ref3 = E.stamp_core(gauss_stamp(1.2)), E.stamp_core(gauss_stamp(1.5)), E.stamp_core(gauss_stamp(1.4))
    assert np.allclose(raw['S|m92|x0_y0|ref'], ref0) and np.allclose(raw['S|m92|x5_y5|ref'], ref1)
    assert np.allclose(raw['S|m92|all|ref'], (16 * ref0 + 9 * ref1 + 3 * ref3) / 28) and not np.allclose(raw['S|m92|all|ref'], (ref0 + ref1 + ref3) / 3)
    # the stack of the pooled stars is the median of the patches' stars, the bootstrap seeds differ by patch
    assert not np.array_equal(raw['S|m92|x0_y0|ci'], raw['S|m92|x5_y5|ci']) and raw['S|m92|all|real'].shape == (19, 19)
    again = str(tmp_path / 'again.npz')
    E.part_stamps(again, B=20)
    assert all(np.array_equal(raw[k], X.load_raw(again)[k]) for k in raw)                                                                                       # reproducible


def test_the_command_line_runs_one_part_per_flag_and_regenerates_the_report_from_the_raw(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(E, 'part_stamps', lambda *a, **k: calls.append(('stamps', a, k)))
    monkeypatch.setattr(E, 'part_synth', lambda *a, **k: calls.append(('synth', a, k)))
    monkeypatch.setattr(E, 'part_real', lambda *a, **k: calls.append(('real', a, k)))
    monkeypatch.setattr(E, 'part_bright', lambda *a, **k: calls.append(('bright', a, k)))
    raw_path, md_path = str(tmp_path / 'diag_raw.npz'), str(tmp_path / 'diag.md')
    np.savez_compressed(raw_path, marker=np.array([7]))
    monkeypatch.setattr(E, 'RAW', raw_path)
    monkeypatch.setattr(E, 'MD', md_path)
    monkeypatch.setattr(E, 'report', lambda raw: 'REPORT %d\n' % raw['marker'][0])

    def run(*argv):
        monkeypatch.setattr(sys, 'argv', ['exp_astro_diag.py'] + list(argv))
        E.main()
    run('--stamps')
    run('--synth')
    run('--synth', '--workers', '3', '--only', 'iso,cs', '--sets', 'm92,m53')
    run('--real')
    run('--bright')
    assert [c[0] for c in calls] == ['stamps', 'synth', 'synth', 'real', 'bright']
    assert calls[1][2] == dict(workers=1, kinds=None, sets=None) and calls[2][2] == dict(workers=3, kinds=['iso', 'cs'], sets=['m92', 'm53'])
    assert not os.path.exists(md_path)
    run('--from-raw')
    assert open(md_path).read() == 'REPORT 7\n' and len(calls) == 5                                 # the report is made from the raw alone: no part is run
    with pytest.raises(SystemExit):
        run()


needs_data = pytest.mark.skipif(not all(os.path.exists(E.patch_file(n, p)) for n, p in (('m2', (310, 630)), ('m92', (1100, 220)))), reason='the astronomy patches are not on disk')


@needs_data
def test_on_a_real_patch_the_world_has_the_real_frames_scale_and_the_synthetic_stack_is_the_stamps():
    patch, sky, truth, skysig = E.load_inputs('m92', (1100, 220))
    w = E.render_world(truth, 1100, 220, patch['psf'], patch['gain'], patch['dark'], sky, E.world_seed('m92', (1100, 220)))
    st = E.frame_stats(patch['dn']), E.frame_stats(w['dn'])
    assert all(abs(st[1][i] / st[0][i] - 1) < 0.25 for i in range(3))      # the world is as bright as the frame to within a quarter (it lacks what the truth list lacks)
    idx = E.stack_stars(w['x'], w['y'], w['mag'])
    assert 10 < len(idx) < 60
    ref = E.stamp_core(patch['psf'])
    syn = E.cutouts(w['dn'], w['x'], w['y'], idx)
    est = E.stack_estimates(syn, syn, ref)
    assert -0.3 < est[0] < 0.1 and est[2] == 0.0 and abs(est[7]) < 1e-12 and est[8] == 0.0                          # centring costs a few per cent of the peak, a stack against itself has no difference
    assert abs(est[4] / E.fwhm(ref) - 1) < 0.1


@needs_data
def test_the_world_of_a_real_patch_is_deterministic_and_scored_truth_is_injected_truth():
    patch, sky, truth, skysig = E.load_inputs('m2', (310, 630))
    a = E.render_world(truth, 310, 630, patch['psf'], patch['gain'], patch['dark'], sky, E.world_seed('m2', (310, 630)))
    b = E.render_world(truth, 310, 630, patch['psf'], patch['gain'], patch['dark'], sky, E.world_seed('m2', (310, 630)))
    assert np.array_equal(a['dn'], b['dn']) and skysig > 0
    t = A.truth_for_patch(truth, 310, 630)
    x, y, m = E.sources_near(truth, 310, 630, 0)
    assert np.array_equal(x, t['x']) and np.array_equal(m, t['mag']) and len(a['x']) > len(x)
