"""The in-situ PSF study (experiments/exp_astro_insitu.py): detection and the selection rules of the stars, the sky, the aligned cut-outs and their sub-pixel offsets, the flux-weighted stack with outlier rejection, the
hybrid stamp, the alignment-bias correction, the recovery of a PLANTED stamp from planted stars, my own matcher check, the raw bookkeeping, the runs (with the run functions replaced), Task B's bookkeeping and the report.
Everything is planted: no astronomy data is read except in the tests marked needs_data."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import astro_common as A                # noqa: E402
import exp_astro_diag as D              # noqa: E402
import exp_astro_insitu as E            # noqa: E402
import exp_astro_m2 as X                # noqa: E402

GAIN, DARK, SKY = 4.6, 1.0, 300.0
REAL_TARGETS = list(E.TARGETS)                                         # the module-scoped planted fixture replaces E.TARGETS for the whole module


def lumpy_stamp(sx=1.5, sy=1.1, lump=0.12, n=51):
    """A point-sampled elliptical Gaussian with an off-centre lump (so that a symmetric stack would be wrong), unit sum; the middle pixel is the centre."""
    yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2.0
    g = np.exp(-0.5 * (xx / sx) ** 2 - 0.5 * (yy / sy) ** 2) + lump * np.exp(-0.5 * ((xx - 1.2) ** 2 + (yy + 0.9) ** 2) / 0.7 ** 2)
    return g / g.sum()


def planted_frame(stamp, stars, seed=0, sky=SKY, noise=True, n=100):
    """A 100 x 100 frame: the stamp rendered (`exp_astro_diag.render_sources`) at the (x, y, flux) of `stars` on a flat sky, with the noise of the methods' model."""
    x, y, f = (np.array(v, float) for v in zip(*stars))
    clean = E.D.render_sources(stamp, x, y, f, (n, n), 20) + sky
    return D.add_noise(clean, GAIN, DARK, np.random.default_rng(seed)) if noise else clean


STARS = [(25.3, 25.7, 1.0e5), (60.6, 27.2, 8.0e4), (27.1, 62.4, 6.0e4), (62.8, 63.3, 9.0e4), (45.5, 45.2, 1.2e5), (75.4, 40.1, 5.0e4), (40.2, 76.9, 7.0e4)]


@pytest.fixture(scope='module')
def planted():
    st = lumpy_stamp()
    dn = planted_frame(st, STARS, seed=1)
    wrong = lumpy_stamp(2.0, 1.8, 0.0)                                  # a 'psField' stamp that does not match the stars: broader, round, no lump
    return st, dn, wrong


# ----------------------------------------------------------------------------------------------------------------------------------
# detection
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_maxima_are_the_planted_stars_at_their_pixels_in_column_row_order_and_the_noise_unit_is_the_models(planted):
    st, dn, _ = planted
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    for sx, sy, f in STARS:
        j = np.flatnonzero((np.abs(x - np.round(sx)) <= 1) & (np.abs(y - np.round(sy)) <= 1))
        assert len(j) == 1, (sx, sy)
        assert pk[j[0]] > 0.01 * f                                      # the smoothed, high-passed peak is a fixed fraction of the flux
    assert noise == pytest.approx(float(np.median(1.0 / E.HP.iss_weights(dn, GAIN, DARK))), rel=1e-12)
    assert 7.0 < noise < 10.0                                          # sqrt(300 / 4.6 + 1) = 8.1
    big = np.argmax(pk)                                                 # the brightest star is (45.5, 45.2)
    assert abs(x[big] - 45.5) <= 1 and abs(y[big] - 45.2) <= 1


def test_a_star_in_a_corner_of_the_frame_is_found_with_its_column_and_row_not_exchanged():
    st = lumpy_stamp()
    dn = planted_frame(st, [(20.0, 60.0, 1.0e5)], seed=3)             # column 20, row 60
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    b = np.argmax(pk)
    assert abs(x[b] - 20) <= 1 and abs(y[b] - 60) <= 1


# ----------------------------------------------------------------------------------------------------------------------------------
# the selection of stars
# ----------------------------------------------------------------------------------------------------------------------------------
def sel(x, y, pk, dn=None, noise=1.0, **kw):
    dn = np.full((100, 100), 100.0) if dn is None else dn
    xs, ys, pks, idx = E.select_stars(np.array(x, float), np.array(y, float), np.array(pk, float), noise, dn, **kw)
    return sorted((float(xs[i]), float(ys[i])) for i in idx), idx, (xs, ys, pks)


def test_a_star_needs_twenty_sigma_and_a_neighbour_counts_from_five_sigma():
    got, idx, _ = sel([50, 20], [50, 20], [19.99, 20.0])
    assert got == [(20.0, 20.0)]
    got, _, _ = sel([50], [50], [20.0])
    assert got == [(50.0, 50.0)]                                       # exactly 20 sigma is a star (>=)
    got, _, _ = sel([50], [50], [19.999])
    assert got == []
    # a neighbour of 4.99 sigma is not a maximum at all; of 5.0 sigma it is one, and (being 0.25 of the star) inside 5 px it excludes the star
    got, _, _ = sel([50, 53], [50, 50], [100.0, 4.99])
    assert got == [(50.0, 50.0)]
    got, _, _ = sel([50, 53], [50, 50], [100.0, 5.0])
    assert got == [(50.0, 50.0)]                                       # 5 sigma is 0.05 of the star: below the 0.1 inner fraction
    got, _, _ = sel([50, 53], [50, 50], [100.0, 10.0])
    assert (50.0, 50.0) not in got                                      # 0.1 of the star at 3 px (the inner rule is >=); the faint one is not a star (10 < 20 sigma)
    assert got == []


def test_the_inner_neighbour_rule_is_a_fraction_within_five_px_inclusive_and_the_outer_rule_a_larger_fraction_within_twelve():
    # inner: fraction 0.1, radius 5 (distance <= 5 counts)
    got, _, _ = sel([50, 55], [50, 50], [1000.0, 100.0])               # exactly 0.1 at exactly 5 px: excluded
    assert (50.0, 50.0) not in got
    got, _, _ = sel([50, 55.0001], [50, 50], [1000.0, 100.0])          # 5.0001 px: not inside the inner radius, and 0.1 < 0.5: kept
    assert (50.0, 50.0) in got
    got, _, _ = sel([50, 55], [50, 50], [1000.0, 99.9])                # below 0.1 at 5 px: kept
    assert (50.0, 50.0) in got
    # outer: fraction 0.5 within 12 px (inclusive)
    got, _, _ = sel([50, 62], [50, 50], [1000.0, 500.0])
    assert (50.0, 50.0) not in got
    got, _, _ = sel([50, 62], [50, 50], [1000.0, 499.9])
    assert (50.0, 50.0) in got
    got, _, _ = sel([50, 62.01], [50, 50], [1000.0, 900.0])            # beyond 12 px: not a neighbour
    assert (50.0, 50.0) in got
    # a brighter neighbour anywhere inside 12 px excludes the star, a fainter one at 8 px does not
    got, _, _ = sel([50, 58], [50, 50], [1000.0, 300.0])
    assert (50.0, 50.0) in got and (58.0, 50.0) not in got


def test_the_edge_and_saturation_cuts_and_the_ordering_by_brightness_and_the_cap():
    for x, y, ok in ((12, 50, False), (13, 50, True), (86, 50, True), (87, 50, False), (50, 12, False), (50, 13, True), (50, 86, True), (50, 87, False)):
        assert (sel([x], [y], [100.0])[0] == [(float(x), float(y))]) == ok, (x, y)           # more than 12 px from every edge of the 100 px patch
    assert sel([50], [50], [100.0], n=100, edge=49)[0] == []              # the edge distance and the patch size are arguments
    assert sel([50], [50], [100.0], n=100, edge=48)[0] == [(50.0, 50.0)]
    assert sel([50], [50], [100.0], n=200, edge=49)[0] == [(50.0, 50.0)]
    # saturation: any pixel above the limit in the 7 x 7 box (3 px each side), not the 9 x 9
    dn = np.full((100, 100), 100.0)
    dn[50 + 3, 50 - 3] = E.SAT_DN + 1
    assert sel([50], [50], [100.0], dn=dn)[0] == []
    dn[50 + 3, 50 - 3] = E.SAT_DN
    assert sel([50], [50], [100.0], dn=dn)[0] == [(50.0, 50.0)]
    dn[50 + 3, 50 - 3] = 100.0
    dn[50 + 4, 50] = E.SAT_DN + 1
    assert sel([50], [50], [100.0], dn=dn)[0] == [(50.0, 50.0)]
    # a star at the very edge of the frame array does not index out of range on the low side
    dn2 = np.full((100, 100), 100.0)
    assert sel([50], [50], [100.0], dn=dn2, edge=0)[0] == [(50.0, 50.0)]
    # ordered by peak, brightest first, and capped
    xs = [20 + 15 * i for i in range(5)]
    got, idx, (sx, sy, spk) = sel(xs, [50] * 5, [30.0, 90.0, 50.0, 70.0, 60.0])
    assert list(spk[idx]) == [90.0, 70.0, 60.0, 50.0, 30.0]
    got, idx, (sx, sy, spk) = sel(xs, [50] * 5, [30.0, 90.0, 50.0, 70.0, 60.0], n_max=2)
    assert list(spk[idx]) == [90.0, 70.0]


def test_the_maxima_returned_are_the_neighbours_of_five_sigma_and_a_fainter_maximum_never_makes_a_star():
    xs, ys, pks, idx = E.select_stars(np.array([10.0, 50.0, 80.0]), np.array([10.0, 50.0, 80.0]), np.array([4.0, 6.0, 30.0]), 1.0, np.full((100, 100), 100.0))
    assert list(pks) == [6.0, 30.0] and list(xs) == [50.0, 80.0]
    assert list(idx) == [1]                                             # indices into the returned (filtered) arrays


# ----------------------------------------------------------------------------------------------------------------------------------
# the sky and the cut-outs
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_clipped_mean_ignores_outliers_but_is_a_mean_not_a_median():
    rng = np.random.default_rng(0)
    v = rng.normal(100.0, 3.0, 400)
    assert abs(E.sky_clip(v) - v.mean()) < 0.3                         # no outlier: about the mean
    w = np.r_[v, [900.0, 1500.0, 700.0]]
    assert abs(E.sky_clip(w) - v.mean()) < 0.5 and w.mean() > 105
    skew = np.r_[np.zeros(300), np.full(100, 40.0)]                    # a skewed ring (crowd): the mean (10) is above the median (0); clipping at 3 sigma keeps it
    assert E.sky_clip(skew) == pytest.approx(10.0, abs=1e-9) and np.median(skew) == 0.0
    assert E.sky_clip(np.full(50, 7.0)) == 7.0                          # zero spread


def test_a_cut_out_is_centred_on_the_star_and_its_centroid_offset_is_the_subpixel_phase_with_the_axes_the_right_way_round(planted):
    st, dn, _ = planted
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    xs, ys, pks, idx = E.select_stars(x, y, pk, noise, dn)
    assert len(idx) == len(STARS)
    c, off = E.star_cores(dn, xs, ys, idx)
    assert c.shape == (len(idx), E.CORE, E.CORE) and off.shape == (len(idx), 2)
    for k, i in enumerate(idx):
        sx, sy, f = min(STARS, key=lambda s: (s[0] - xs[i]) ** 2 + (s[1] - ys[i]) ** 2)
        # the centroid of the cut-out measures the star's centroid (the stamp has its lump: not the position; the offset is the star's position relative to the maximum pixel, plus the stamp's own centroid offset)
        cx0, cy0 = E.D.centroid(st[25 - 12:25 + 13, 25 - 12:25 + 13])
        assert off[k, 0] == pytest.approx(sx - np.round(xs[i]) + (cx0 - 12), abs=0.12)
        assert off[k, 1] == pytest.approx(sy - np.round(ys[i]) + (cy0 - 12), abs=0.12)
        cy, cx = np.unravel_index(np.argmax(c[k]), c[k].shape)
        assert (cy, cx) in ((9, 9), (9, 8), (9, 10), (8, 9), (10, 9), (8, 8), (10, 10), (8, 10), (10, 8))
        assert abs(E.D.centroid(np.pad(c[k], 3))[0] - 12) < 0.15 and abs(E.D.centroid(np.pad(c[k], 3))[1] - 12) < 0.15         # centred on the middle pixel to a tenth of a pixel
    # the sky is subtracted: the ring of the aligned core is about zero for a star on a flat sky
    assert abs(np.median(c[0][:, 0])) < 40.0


def test_the_cut_out_amplitude_follows_the_flux_and_the_sky_is_removed_even_when_the_background_is_a_gradient_free_pedestal():
    st = lumpy_stamp()
    for sky in (0.0, 300.0, 3000.0):
        dn = planted_frame(st, [(50.0, 50.0, 1.0e5), (25.0, 25.0, 5.0e4)], seed=2, sky=sky, noise=False)
        c, off = E.star_cores(dn, np.array([50.0, 25.0]), np.array([50.0, 25.0]), np.array([0, 1]))
        assert c[0][E.DISK].sum() / c[1][E.DISK].sum() == pytest.approx(2.0, rel=0.02)
        assert c[0][9, 9] > 0.05 * 1.0e5                                                                                  # the sky does not eat the star


def test_a_noise_free_star_at_an_integer_position_has_zero_offset_and_the_stamp_back():
    st = lumpy_stamp()
    dn = planted_frame(st, [(50.0, 50.0, 1.0e5)], noise=False, sky=0.0)
    c, off = E.star_cores(dn, np.array([50.0]), np.array([50.0]), np.array([0]))
    cx0, cy0 = E.D.centroid(st[25 - 12:25 + 13, 25 - 12:25 + 13])
    assert off[0, 0] == pytest.approx(cx0 - 12, abs=1e-9) and off[0, 1] == pytest.approx(cy0 - 12, abs=1e-9)
    assert np.abs(c[0] / 1.0e5 - st[25 - 9:25 + 10, 25 - 9:25 + 10]).max() < 0.05 * st.max()      # centred: the stamp (its own centroid is the lump's offset, a sub-pixel shift)


# ----------------------------------------------------------------------------------------------------------------------------------
# the stack
# ----------------------------------------------------------------------------------------------------------------------------------
def shape_core(sx=1.4, sy=1.4):
    yy, xx = np.mgrid[-9:10, -9:10]
    g = np.exp(-0.5 * (xx / sx) ** 2 - 0.5 * (yy / sy) ** 2)
    return g / g[E.DISK].sum()


def test_the_stack_recovers_the_common_shape_and_the_amplitudes_and_normalises_to_the_six_pixel_disc():
    s0 = shape_core(1.4, 1.1)
    amps = np.array([5.0, 1.0, 12.0, 3.0, 0.5])
    rng = np.random.default_rng(1)
    c = amps[:, None, None] * s0 + 0.0005 * rng.standard_normal((5, 19, 19))
    s, f, e, use = E.robust_stack(c)
    assert use.all() and np.allclose(f, amps, rtol=0.03)
    assert s[E.DISK].sum() == pytest.approx(1.0, abs=1e-12)
    assert np.abs(s - s0).max() < 0.01 * s0.max()
    assert np.all(e < 0.1) and e.shape == (5,)


def test_a_star_with_a_bright_companion_in_its_cutout_is_rejected_and_does_not_change_the_stack():
    s0 = shape_core(1.4, 1.1)
    amps = np.array([5.0, 3.0, 8.0, 4.0, 6.0])
    c = amps[:, None, None] * s0
    base, _, _, _ = E.robust_stack(c)
    bad = c.copy()
    bad[2, 11:14, 5:8] += 1.0 * amps[2] * s0.max()                         # a neighbour of the third star's own peak height, 3 to 4.5 px from its centre (inside the disc)
    s, f, e, use = E.robust_stack(bad)
    assert list(use) == [True, True, False, True, True]
    assert np.abs(s - base).max() < 1e-6 * s0.max()                         # the rest are exactly proportional to the shape: the stack is the shape whether or not the star is there
    assert e[2] > E.KEEP_FLOOR and np.all(e[[0, 1, 3, 4]] < 1e-6)
    far = c.copy()
    far[2, 3:6, 14:17] += 1.0 * amps[2] * s0.max()                         # the same neighbour 7 to 8 px away: outside the 6 px disc that the misfit is measured in, so kept
    assert E.robust_stack(far)[3].all()


def test_a_star_with_a_negative_amplitude_is_dropped_and_fewer_than_three_stars_give_nothing():
    s0 = shape_core()
    c = np.array([3.0, 2.0, 4.0, 5.0])[:, None, None] * s0
    c[1] = -c[1]
    s, f, e, use = E.robust_stack(c)
    assert list(use) == [True, False, True, True] and f[1] < 0
    s2, f2, e2, use2 = E.robust_stack(c[:3])
    assert np.isnan(s2).all() and not use2.any()                       # only two stars with a positive sum remain
    s3, f3, e3, use3 = E.robust_stack(np.array([3.0, 2.0])[:, None, None] * s0)
    assert np.isnan(s3).all() and not use3.any() and len(f3) == 2


def flip_point(median_misfit, n=8):
    """Bisects the height of a blob (3 x 3 px, inside the disc) added to the last of n stars (the others carry noise at `median_misfit` of the peak) between kept and rejected: returns (e of the star at the
    largest height still kept, its e at the smallest height rejected, the rule's threshold (3 x median of the others' e, floor 0.10) at each)."""
    s0 = shape_core()
    rng = np.random.default_rng(5)
    noise = rng.standard_normal((n, 19, 19))
    noise = noise / noise[:, E.DISK].std(1)[:, None, None]
    pk = 4.0 * s0.max()

    def run(h):
        c = 4.0 * s0 + median_misfit * pk * noise
        c[-1, 11:14, 5:8] += h * pk
        s, f, e, use = E.robust_stack(c)
        return use[-1], e, use
    lo, hi = 0.0, 4.0
    assert run(lo)[0] and not run(hi)[0]
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if run(mid)[0] else (lo, mid)
    (k1, e1, u1), (k2, e2, u2) = run(lo), run(hi)
    assert k1 and not k2 and u1[:-1].all() and u2[:-1].all()
    thr = lambda e: max(3.0 * float(np.median(e[:-1])), 0.10)                  # noqa: E731
    return e1[-1], e2[-1], thr(e1), thr(e2)


def test_the_threshold_of_rejection_is_three_times_the_median_misfit_with_a_floor_of_a_tenth_of_the_peak():
    e_keep, e_rej, t_keep, t_rej = flip_point(0.02)                     # small misfits everywhere: the floor binds (3 x 0.02 < 0.10)
    assert t_keep == t_rej == 0.10 and 0.08 < e_keep <= 0.10 < e_rej < 0.13
    # a median misfit of a tenth: 3 x the median (0.33, measured against the first, median stack) binds, the floor does not; the planted blob flips the star between 1.10 and 1.12 of the peak
    s0 = shape_core()
    rng = np.random.default_rng(5)
    noise = rng.standard_normal((8, 19, 19))
    noise = noise / noise[:, E.DISK].std(1)[:, None, None]
    pk = 4.0 * s0.max()
    for h, kept in ((0.5, True), (1.10, True), (1.12, False), (2.0, False)):
        c = 4.0 * s0 + 0.10 * pk * noise
        c[-1, 11:14, 5:8] += h * pk
        use = E.robust_stack(c)[3]
        assert use[-1] == kept and use[:-1].all(), h


# ----------------------------------------------------------------------------------------------------------------------------------
# the hybrid stamp and the forward model
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_hybrid_stamp_is_the_core_inside_the_psfield_wings_outside_unit_sum_and_non_negative():
    psf = lumpy_stamp(1.6, 1.6, 0.0)
    core = shape_core(1.1, 1.1)
    h = E.hybrid_stamp(core, psf)
    assert h.shape == psf.shape and h.sum() == pytest.approx(1.0, abs=1e-12) and h.min() >= 0.0
    r = np.hypot(*np.mgrid[-25:26, -25:26])
    far = r > E.R_OUT + 1
    assert np.allclose(h[far] / psf[far], h[far][0] / psf[far][0], rtol=1e-9)      # the wings are the psField stamp (up to the one renormalisation)
    inner = r <= E.R_FULL
    a = psf[r <= E.NORM_RAD].sum()
    assert np.allclose(h[inner] * (1.0 / h.sum()) * h.sum(), h[inner])
    # the core carries the psField flux within 6 px: (the unit-sum renormalisation aside) the stamp's 6-px sum over the stamp's 6-px sum of psField is within a few per cent of 1
    assert h[r <= E.NORM_RAD].sum() == pytest.approx(a, rel=0.05)
    # a core 20 % higher in the middle is 20 % higher in the middle of the stamp (relative to a), the shape inside 6.5 px follows the core
    narrow = E.hybrid_stamp(shape_core(0.9, 0.9), psf)
    assert narrow[25, 25] > h[25, 25] > psf[25, 25] * 0.5
    ratio = narrow[25, 25] / h[25, 25]
    assert ratio == pytest.approx(shape_core(0.9, 0.9)[9, 9] / shape_core(1.1, 1.1)[9, 9], rel=0.05)
    # negative pixels of the core are clipped
    core2 = core.copy()
    core2[9, 9] = -1.0
    assert E.hybrid_stamp(core2, psf).min() == 0.0


def test_the_blend_is_linear_between_the_full_radius_and_the_outer_radius():
    psf = np.full((51, 51), 1.0 / 51 ** 2)
    core = np.full((19, 19), 0.0)
    core[E.DISK] = 1.0 / E.DISK.sum()
    h = E.hybrid_stamp(core, psf)
    r = np.hypot(*np.mgrid[-25:26, -25:26])
    a = psf[r <= 6.0].sum()
    # the weight of the core at radius r is clip((9 - r) / 2.5, 0, 1): the stamp before renormalisation is w a core + (1 - w) psf
    raw = np.zeros((51, 51))
    raw[25 - 9:25 + 10, 25 - 9:25 + 10] = a * core
    w = np.clip((9.0 - r) / 2.5, 0, 1)
    exp = w * raw + (1 - w) * psf
    exp /= exp.sum()
    assert np.allclose(h, exp, atol=1e-15)
    assert w[25, 25 + 7] == pytest.approx(0.8) and w[25, 25 + 6] == 1.0 and w[25, 25 + 9] == 0.0


def test_the_forward_model_of_the_cutting_procedure_reproduces_the_aligned_normalised_core_of_a_noise_free_star():
    st = lumpy_stamp()
    dn = planted_frame(st, [(50.3, 50.6, 1.0e5)], noise=False, sky=0.0)
    c, off = E.star_cores(dn, np.array([50.0]), np.array([51.0]), np.array([0]))
    obs = c[0] / c[0][E.DISK].sum()
    # the star's centroid offset from its maximum pixel, in the units of the shifts of `forward_core`
    t = E.forward_core(st, off, np.array([1.0]))
    assert np.abs(t - obs).max() < 0.01 * obs.max()
    # weights: two identical stars, any weights, give the same core; two different offsets are averaged with the weights
    t2 = E.forward_core(st, np.vstack([off, off]), np.array([0.3, 0.7]))
    assert np.allclose(t2, t, atol=1e-12)
    a = E.forward_core(st, np.array([[0.4, -0.3]]), np.array([1.0]))
    b = E.forward_core(st, np.array([[-0.2, 0.45]]), np.array([1.0]))
    ab = E.forward_core(st, np.array([[0.4, -0.3], [-0.2, 0.45]]), np.array([1.0, 3.0]))
    mix = 0.25 * a + 0.75 * b
    assert np.allclose(ab, mix / mix[E.DISK].sum(), atol=1e-12) and ab[E.DISK].sum() == pytest.approx(1.0)


def test_the_forward_model_uses_the_offset_with_x_along_the_columns_and_y_along_the_rows():
    st = lumpy_stamp(2.0, 1.0, 0.0)
    ox_only = E.forward_core(st, np.array([[0.45, 0.0]]), np.array([1.0]))
    oy_only = E.forward_core(st, np.array([[0.0, 0.45]]), np.array([1.0]))
    base = E.forward_core(st, np.array([[0.0, 0.0]]), np.array([1.0]))
    # the stamp is wider along x (columns): the aligned core of a star shifted along x differs from the unshifted one in a different way from one shifted along y
    assert np.abs(ox_only - base).max() != pytest.approx(np.abs(oy_only - base).max(), rel=1e-3)
    # a symmetric stamp shifted by +dx and by -dx and then centred gives the cores that are mirror images of each other along x (the centring undoes the shift, the cubic spline's rounding is mirrored)
    plus = E.forward_core(st, np.array([[0.45, 0.0]]), np.array([1.0]))
    minus = E.forward_core(st, np.array([[-0.45, 0.0]]), np.array([1.0]))
    assert np.allclose(plus, minus[:, ::-1], atol=2e-3 * plus.max())
    assert np.allclose(plus.sum(), minus.sum())


# ----------------------------------------------------------------------------------------------------------------------------------
# the estimator on planted stars
# ----------------------------------------------------------------------------------------------------------------------------------
def test_a_planted_stamp_is_recovered_from_planted_stars_where_the_psfield_stamp_is_wrong(planted):
    st, dn, wrong = planted
    ref = E.stamp_core(st)
    est = E.insitu_stamp([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    assert est['ok'] and len(est['stars']['c']) == len(STARS) and est['use'].all()
    core = E.stamp_core(est['stamp'])
    v = E.compare_cores(core, ref)
    assert abs(v[0]) < 0.03 and abs(v[1] - v[2]) < 0.06 and v[3] < 0.015, v
    vw = E.compare_cores(E.stamp_core(wrong), ref)
    assert abs(vw[0]) > 0.15 and vw[3] > 0.04                          # the psField stamp is far from the stars; the in-situ stamp is not
    # the elliptical, lumpy shape is recovered, not a round one
    assert np.abs(est['stamp'][25 - 6:25 + 7, 25 - 6:25 + 7] - st[25 - 6:25 + 7, 25 - 6:25 + 7]).max() < 0.05 * st.max()
    assert est['stamp'].sum() == pytest.approx(1.0)
    assert est['resid'] < 0.02


def test_the_cubic_spline_alignment_bias_is_removed_by_the_fixed_point_and_shows_without_it():
    st = lumpy_stamp(1.0, 0.9, 0.0)                                     # narrow: undersampled, so the spline round trip changes the peak
    rng = np.random.default_rng(4)
    stars = [(20.0 + 12.0 * (i % 6) + rng.uniform(-0.5, 0.5), 20.0 + 14.0 * (i // 6) + rng.uniform(-0.5, 0.5), 1.0e5) for i in range(12)]
    dn = planted_frame(st, stars, noise=False, sky=0.0)
    wrong = lumpy_stamp(1.4, 1.4, 0.0)
    ref = E.stamp_core(st)
    est0 = E.estimate(*E.star_cores(dn, *np.array([[s[0] for s in stars], [s[1] for s in stars]]), np.arange(12)), wrong, n_fix=0)
    est3 = E.estimate(*E.star_cores(dn, *np.array([[s[0] for s in stars], [s[1] for s in stars]]), np.arange(12)), wrong, n_fix=3)
    b0 = abs(E.compare_cores(E.stamp_core(est0['stamp']), ref)[0])
    b3 = abs(E.compare_cores(E.stamp_core(est3['stamp']), ref)[0])
    assert b0 > 0.02 and b3 < 0.01 and b3 < 0.5 * b0, (b0, b3)


def test_too_few_stars_give_no_stamp_and_a_nan_stamp_is_reported_as_not_ok(planted):
    st, dn, wrong = planted
    few = dict(dn=planted_frame(st, STARS[:2], seed=2), gain=GAIN, dark=DARK)
    est = E.insitu_stamp([few], wrong)
    assert not est['ok'] and np.isnan(est['stamp']).all() and len(est['stars']['c']) == 2
    est0 = E.insitu_stamp([dict(dn=np.full((100, 100), 300.0) + np.random.default_rng(0).normal(0, 3, (100, 100)), gain=GAIN, dark=DARK)], wrong)
    assert not est0['ok'] and len(est0['stars']['c']) == 0 and est0['use'].shape == (0,)


def test_several_frames_pool_their_stars_and_the_star_table_follows_the_frames(planted):
    st, dn, wrong = planted
    dn2 = planted_frame(st, [(30.0, 30.0, 8.0e4), (70.0, 70.0, 9.0e4), (30.0, 70.0, 7.0e4)], seed=9)
    est = E.insitu_stamp([dict(dn=dn, gain=GAIN, dark=DARK), dict(dn=dn2, gain=GAIN, dark=DARK)], wrong)
    s = est['stars']
    assert list(s['n_by_frame']) == [len(STARS), 3] and len(s['c']) == len(STARS) + 3
    assert list(s['frame']) == [0] * len(STARS) + [1] * 3 and s['snr'].shape == (len(STARS) + 3,)

    second = {(round(a), round(b)) for a, b in zip(s['x'][len(STARS):], s['y'][len(STARS):])}
    assert second == {(30, 30), (70, 70), (30, 70)}
    v = E.compare_cores(E.stamp_core(est['stamp']), E.stamp_core(st))
    assert abs(v[0]) < 0.03


def test_the_stamp_core_and_the_comparison_statistics_have_the_documented_meaning():
    a = shape_core(1.4)
    b = shape_core(1.4) * 0.9
    b[9, 9] = a[9, 9] * 0.9
    v = E.compare_cores(a, b)
    assert v[0] == pytest.approx(1 / 0.9 - 1)
    assert v[3] == pytest.approx(np.sqrt(np.mean((a - b) ** 2)) / b[9, 9])
    assert v[1] == pytest.approx(E.D.fwhm(a)) and v[2] == pytest.approx(E.D.fwhm(b))
    st = lumpy_stamp(1.5, 1.5, 0.0)
    core = E.stamp_core(st)
    assert core.shape == (19, 19) and core[E.DISK].sum() == pytest.approx(1.0)
    sv = E.stat_vector(st, core)
    assert sv[0] == pytest.approx(0.0, abs=1e-12) and sv[2] == pytest.approx(0.0, abs=1e-12) and sv[1] == pytest.approx(E.D.fwhm(core))
    assert np.isnan(E.stat_vector(np.full((51, 51), np.nan), core)).all()
    sh = E.shape_of(lumpy_stamp(2.0, 1.0, 0.0))
    assert sh[0] > 0.3 and abs(sh[1]) < 5.0                             # elongated along x: orientation 0 degrees
    sh2 = E.shape_of(lumpy_stamp(1.0, 2.0, 0.0))
    assert sh2[0] > 0.3 and abs(abs(sh2[1]) - 90.0) < 5.0


def test_the_star_bootstrap_is_reproducible_stratified_and_brackets_the_estimate(planted):
    st, dn, wrong = planted
    est = E.insitu_stamp([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    ref = E.stamp_core(st)
    ci, nf = E.star_bootstrap(est['stars'], [ref], [wrong], est['stars']['n_by_frame'], B=60, seed=0)
    ci2, nf2 = E.star_bootstrap(est['stars'], [ref], [wrong], est['stars']['n_by_frame'], B=60, seed=0)
    assert ci.shape == (1, 3, 2) and np.array_equal(ci, ci2) and nf == nf2 == 0
    pt = E.stat_vector(est['stamp'], ref)
    assert np.all(ci[0, :, 0] <= ci[0, :, 1])
    assert ci[0, 0, 0] - 0.05 < pt[0] < ci[0, 0, 1] + 0.05
    ci3, _ = E.star_bootstrap(est['stars'], [ref], [wrong], est['stars']['n_by_frame'], B=60, seed=1)
    assert not np.array_equal(ci, ci3)
    # two strata: the resamples stay within their stratum (a star never moves to the other)
    s2 = dict(est['stars'])
    idx = X.boot_idx_strat([3, 4], 5, 0)
    assert (idx[:, :3] < 3).all() and (idx[:, 3:] >= 3).all()


# ----------------------------------------------------------------------------------------------------------------------------------
# my own matcher check
# ----------------------------------------------------------------------------------------------------------------------------------
def brute_nearest(x, y, tx, ty, tm, rad):
    mag, ind = [], []
    for xi, yi in zip(x, y):
        best, bj = None, -1
        for j in range(len(tx)):
            d = ((tx[j] - xi) ** 2 + (ty[j] - yi) ** 2) ** 0.5
            if d <= rad and (best is None or d < best):
                best, bj = d, j
        mag.append(tm[bj] if bj >= 0 else np.nan)
        ind.append(bj)
    return np.array(mag), np.array(ind)


def test_the_nearest_truth_matcher_agrees_with_a_brute_force_loop_including_ties_and_the_boundary():
    rng = np.random.default_rng(3)
    tx, ty, tm = rng.uniform(0, 30, 60), rng.uniform(0, 30, 60), rng.uniform(15, 25, 60)
    x, y = rng.uniform(0, 30, 40), rng.uniform(0, 30, 40)
    for rad in (0.5, 1.0, 2.5):
        m1, i1 = E.nearest_truth(x, y, tx, ty, tm, rad)
        m2, i2 = brute_nearest(x, y, tx, ty, tm, rad)
        assert np.array_equal(i1, i2) and np.array_equal(m1, m2, equal_nan=True)
    # boundary: exactly at rad counts; the nearest wins; no truth at all
    m, i = E.nearest_truth(np.array([0.0]), np.array([0.0]), np.array([1.0, 0.6]), np.array([0.0, 0.0]), np.array([20.0, 21.0]), 1.0)
    assert i[0] == 1 and m[0] == 21.0
    m, i = E.nearest_truth(np.array([0.0]), np.array([0.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), 1.0)
    assert i[0] == 0
    m, i = E.nearest_truth(np.array([0.0]), np.array([0.0]), np.array([1.0001]), np.array([0.0]), np.array([20.0]), 1.0)
    assert i[0] == -1 and np.isnan(m[0])
    m, i = E.nearest_truth(np.array([0.0]), np.array([0.0]), np.zeros(0), np.zeros(0), np.zeros(0), 1.0)
    assert i[0] == -1


def test_truth_isolation_in_diags_sense_is_four_and_a_half_pixels_and_one_and_a_half_magnitudes():
    tx = np.array([10.0, 14.5, 10.0, 30.0, 10.0])
    ty = np.array([10.0, 10.0, 14.6, 30.0, 6.0])
    tm = np.array([18.0, 19.4, 17.0, 17.0, 19.5])
    # source 0: neighbours: 1 at exactly 4.5 px, mag 19.4 < 19.5 (brighter than m + 1.5): not isolated; 2 at 4.6 px is too far; 4 at 4.0 px mag 19.5 is not brighter than 19.5 (strict)
    got = E.truth_isolated_mask(np.array([0, 3, -1]), tx, ty, tm)
    assert list(got) == [False, True, False]
    tm2 = tm.copy()
    tm2[1] = 19.5                                                       # now not brighter than m + 1.5 (strict): isolated
    assert list(E.truth_isolated_mask(np.array([0]), tx, ty, tm2)) == [True]
    tx3 = tx.copy()
    tx3[1] = 14.6
    assert list(E.truth_isolated_mask(np.array([0]), tx3, ty, tm)) == [True]


# ----------------------------------------------------------------------------------------------------------------------------------
# the held-out fit in the high-passed image
# ----------------------------------------------------------------------------------------------------------------------------------
def test_a_star_is_fitted_perfectly_by_its_own_stamp_and_worse_by_a_wrong_one_and_the_misfit_columns_are_in_order(planted):
    st, dn, wrong = planted
    pos = [(s[0], s[1]) for s in STARS[:3]]
    r_true = E.hp_fit_rms(dn, GAIN, DARK, st, pos)
    r_wrong = E.hp_fit_rms(dn, GAIN, DARK, wrong, pos)
    assert r_true.shape == (3, 2) and r_wrong.shape == (3, 2)
    assert np.all(r_true[:, 0] < 2.0) and np.all(r_true[:, 1] < 0.02)       # standardized rms near 1 (noise only), misfit under 2 % of the peak
    assert np.all(r_wrong[:, 1] > 3 * r_true[:, 1]) and np.all(r_wrong[:, 0] > 2 * r_true[:, 0])
    out = E.hp_fit_rms(dn, GAIN, DARK, st, [(2.0, 50.0), (50.0, 99.0), (50.0, 50.0)])
    assert np.isnan(out[0]).all() and np.isnan(out[1]).all() and np.isfinite(out[2]).all()       # boxes outside the frame are NaN


# ----------------------------------------------------------------------------------------------------------------------------------
# the raw
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_raw_is_write_once_records_the_hashes_and_refuses_a_changed_file(tmp_path, monkeypatch):
    path = str(tmp_path / 'insitu_raw.npz')
    E.save_with_meta({'a|x': np.arange(3)}, path)
    raw = E.load_raw(path)
    assert 'meta|files' in raw and 'meta|sha256' in raw and 'meta|versions' in raw and list(raw['a|x']) == [0, 1, 2]
    E.save_with_meta({'b|x': np.arange(2)}, path)                          # the hashes are recorded once
    raw2 = E.load_raw(path)
    assert np.array_equal(raw2['meta|sha256'], raw['meta|sha256']) and set(raw2) - set(raw) == {'b|x'}
    E.save_with_meta({'a|x': np.arange(3)}, path)                          # the same value again: allowed
    with pytest.raises(RuntimeError, match='write-once'):
        E.save_with_meta({'a|x': np.arange(4)}, path)
    E.check_hashes(path)                                                   # unchanged: no error
    bad = dict(raw2)
    bad['meta|sha256'] = np.array(['0' * 64] + [str(h) for h in raw2['meta|sha256'][1:]])
    np.savez_compressed(path, **bad)
    with pytest.raises(RuntimeError, match='differ'):
        E.check_hashes(path)
    assert E.load_raw(str(tmp_path / 'nothing.npz')) == {}
    E.check_hashes(str(tmp_path / 'nothing.npz'))                          # no raw yet: nothing to check


# ----------------------------------------------------------------------------------------------------------------------------------
# --stamps on planted patches
# ----------------------------------------------------------------------------------------------------------------------------------
def planted_world(monkeypatch, sets=(('m2', (0, 0)), ('m2', (1000, 0)))):
    """Two planted patches of one 'frame' (origins (0, 0) and (1000, 0)); E.inputs and E.TARGETS are replaced. Returns the stamps and the truth."""
    st = lumpy_stamp()
    wrong = lumpy_stamp(2.0, 1.8, 0.0)
    cache = {}
    for k, (name, pxy) in enumerate(sets):
        stars = [(x + (k * 3 % 5) * 0.1, y, f) for x, y, f in STARS]
        dn = planted_frame(st, stars, seed=10 + k)
        truth = dict(px=np.array([s[0] for s in stars]) + pxy[0], py=np.array([s[1] for s in stars]) + pxy[1], m606=A.counts_to_mag(np.array([s[2] for s in stars])))
        cache[(name, pxy)] = (dict(dn=dn, psf=wrong, gain=GAIN, dark=DARK), np.full((100, 100), SKY), truth)
    monkeypatch.setattr(E, 'TARGETS', list(sets))
    monkeypatch.setattr(E, '_CACHE', {})
    monkeypatch.setattr(E, 'inputs', lambda name, pxy: cache[(name, tuple(pxy))])
    return st, wrong, cache


def test_part_stamps_stores_both_scopes_the_validation_and_does_not_recompute(monkeypatch, tmp_path):
    st, wrong, cache = planted_world(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_stamps(path, B=6)
    raw = E.load_raw(path)
    ref = E.stamp_core(st)
    for (name, pxy) in E.TARGETS:
        k = E.pkey(name, pxy)
        for scope in ('patch', 'frame'):
            p = 'ST|%s|%s|' % (scope, k)
            assert raw[p + 'stamp'].shape == (51, 51) and raw[p + 'n_sel'][0] == len(STARS) * (2 if scope == 'frame' else 1) and raw[p + 'n_used'][0] == raw[p + 'n_sel'][0]
            assert raw[p + 'vs_psf'].shape == (4,) and raw[p + 'vs_psf_ci'].shape == (1, 3, 2) and raw[p + 'boot_fail'][0] == 0
            # the in-situ stamp is the planted stamp: its distance from the planted core is small, from the wrong psField core large
            assert abs(E.compare_cores(E.stamp_core(raw[p + 'stamp']), ref)[0]) < 0.04
            assert abs(raw[p + 'vs_psf'][0]) > 0.1
            assert list(raw[p + 'n_by_frame']) == ([len(STARS)] * 2 if scope == 'frame' else [len(STARS)])
        assert np.array_equal(raw['ST|psf|%s|stamp' % k], wrong) and raw['ST|psf|%s|core' % k].shape == (19, 19)
        assert raw['ST|patch|%s|truth|n' % k][0] == len(STARS)                                   # the planted stars are the truth stars (isolated, m < 20)
        assert raw['ST|patch|%s|truth|vs_psf' % k].shape == (4,) and np.isfinite(raw['ST|patch|%s|vs_truth' % k]).all()
        assert abs(raw['ST|patch|%s|vs_truth' % k][0]) < 0.05                                      # the in-situ stamp and the truth-star stamp are the same stamp
        tm = raw['ST|patch|%s|truthmatch' % k]
        assert tm[0] == len(STARS) and tm[1] == len(STARS)                                         # every in-situ star is a truth source and truth-isolated
        assert raw['ST|patch|%s|world|n_sel' % k][0] >= 3 and np.isfinite(raw['ST|patch|%s|world|vs_planted' % k]).all()
        h = raw['ST|patch|%s|held' % k]
        assert h.ndim == 2 and h.shape[1] == 4
        assert raw['ST|frame|%s|vs_truth' % k].shape == (4,)
    for name in ('m2',):
        assert raw['ST|frame|%s|all|truth_n' % name][0] == 2 * len(STARS) and abs(raw['ST|frame|%s|all|insitu_vs_truth' % name][0]) < 0.05
    n_keys = len(raw)
    E.part_stamps(path, B=6)                                                                        # a second call computes nothing and changes nothing
    assert len(E.load_raw(path)) == n_keys
    assert 'meta|files' in raw


def test_the_stars_without_a_truth_source_or_with_a_close_brighter_one_are_counted_in_the_truth_match(monkeypatch):
    st, wrong, cache = planted_world(monkeypatch, sets=(('m2', (0, 0)),))
    patch, sky, truth = cache[('m2', (0, 0))]
    truth = dict(truth)
    truth['px'] = np.r_[truth['px'], truth['px'][0] + 2.0]                                         # a brighter neighbour 2 px from the first star
    truth['py'] = np.r_[truth['py'], truth['py'][0]]
    truth['m606'] = np.r_[truth['m606'], truth['m606'][0] - 3.0]
    truth['px'][2] += 30.0                                                                          # the third star has no truth source any more
    cache[('m2', (0, 0))] = (patch, sky, truth)
    est = E.insitu_stamp([dict(dn=patch['dn'], gain=GAIN, dark=DARK)], wrong)
    sx, sy = est['stars']['x'], est['stars']['y']
    tx, ty, tm = D.sources_near(truth, 0, 0, D.MARGIN)
    mag, ind = E.nearest_truth(sx, sy, tx, ty, tm)
    iso = E.truth_isolated_mask(ind, tx, ty, tm)
    assert np.isnan(mag).sum() == 1 and (ind < 0).sum() == 1
    assert iso.sum() == len(sx) - 2                                                                  # one with no truth source, one with a brighter truth neighbour within 4.5 px


# ----------------------------------------------------------------------------------------------------------------------------------
# --ref
# ----------------------------------------------------------------------------------------------------------------------------------
def fake_scores(seed, nblk=25):
    rng = np.random.default_rng(seed)
    blk = rng.integers(1, 6, size=(nblk, 4)).astype(float)
    blk[:, 1] = np.minimum(blk[:, 1], blk[:, 0])
    blk[:, 3] = np.minimum(blk[:, 3], blk[:, 2])
    tot = np.array([blk[:, 0].sum(), blk[:, 1].sum(), blk[:, 2].sum(), blk[:, 3].sum(), 10, 5, 8])
    return blk, tot


def fake_frozen(tmp_path, targets):
    cyc, hp2, m2 = {}, {}, {}
    seed = 0
    for name, pxy in targets:
        src = E.ref_sources(name, pxy)
        for meth, (rn, pre) in src.items():
            seed += 1
            blk, tot = fake_scores(seed)
            d = dict(cyc=cyc, hp2=hp2, m2=m2)[rn]
            if meth in ('iso', 'adq'):
                for k in ('k2', 'k0'):
                    d['%s|%s|tot' % (pre, k)] = tot + (0 if k == 'k2' else 1)
                    d['%s|%s|blk22' % (pre, k)] = blk
                d[pre + '|secs'] = np.array([100.0 + seed, 1000.0 + seed])
                d[pre + '|det'] = np.zeros((seed, 3))
                d[pre + '|trace|barred'] = np.zeros((seed % 3, 3))
                d[pre + '|bic'] = np.array([5, 50, 0])
            else:
                d[pre + '|tot'] = tot
                d[pre + '|blk22'] = blk
                d[pre + '|secs'] = np.array([3.0])
    paths = []
    for nm, d in (('cyclefix_raw.npz', cyc), ('baselines_hp2_raw.npz', hp2), ('m2_raw.npz', m2)):
        np.savez_compressed(str(tmp_path / nm), **d)
        paths.append(str(tmp_path / nm))
    return paths


def test_the_psfield_rows_are_copied_from_the_right_frozen_keys_and_the_m2_adequate_row_is_the_isotropic_one(monkeypatch, tmp_path):
    targets = [('m2', (310, 630)), ('m53', (760, 590)), ('m92', (1100, 220))]
    monkeypatch.setattr(E, 'TARGETS', targets)
    src = E.ref_sources('m2', (310, 630))
    assert src['iso'] == src['adq'] == ('cyc', 'C|T|x310_y630|ours') and src['cs_def'] == ('m2', 'T|x310_y630|crowdsource (default)')
    assert E.ref_sources('m53', (760, 590))['iso'] == ('cyc', 'C|T2|m53|x760_y590|ours') and E.ref_sources('m53', (760, 590))['adq'] == ('cyc', 'C|P|m53|x760_y590|ours')
    assert E.ref_sources('m92', (1100, 220))['cs_def'] == ('m2', 'T2|m92|x1100_y220|crowdsource (default)')
    assert E.ref_sources('m92', (1100, 220))['cs_hp'] == ('hp2', 'T|m92|x1100_y220|cs2') and E.ref_sources('m92', (1100, 220))['phu_hp'] == ('hp2', 'T|m92|x1100_y220|phu2')
    a, b, c = fake_frozen(tmp_path, targets)
    monkeypatch.setattr(E, 'CYCLEFIX_RAW', a)
    monkeypatch.setattr(E, 'HP2_RAW', b)
    monkeypatch.setattr(E, 'M2_RAW', c)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_ref(path)
    raw = E.load_raw(path)
    cyc, hp2, m2 = (dict(np.load(p)) for p in (a, b, c))
    for name, pxy in targets:
        k = 'REF|%s' % E.pkey(name, pxy)
        s = E.ref_sources(name, pxy)
        assert np.array_equal(raw[k + '|iso|k2|blk22'], cyc[s['iso'][1] + '|k2|blk22']) and np.array_equal(raw[k + '|adq|k0|tot'], cyc[s['adq'][1] + '|k0|tot'])
        assert np.array_equal(raw[k + '|cs_hp|k2|blk22'], hp2[s['cs_hp'][1] + '|blk22']) and np.array_equal(raw[k + '|phu_hp|k2|tot'], hp2[s['phu_hp'][1] + '|tot'])
        assert np.array_equal(raw[k + '|cs_def|k2|tot'], m2[s['cs_def'][1] + '|tot'])
        assert np.array_equal(raw[k + '|iso|bic'], [5, 50, 0]) and raw[k + '|iso|barred'][0] == len(cyc[s['iso'][1] + '|trace|barred'])
        assert raw[k + '|iso|det'].shape == cyc[s['iso'][1] + '|det'].shape and np.array_equal(raw[k + '|iso|secs'], cyc[s['iso'][1] + '|secs'])
    assert list(raw['REF|sources|names']) == ['cyclefix_raw.npz', 'baselines_hp2_raw.npz', 'm2_raw.npz'] and len(raw['REF|sources|sha256'][0]) == 64
    assert raw['REF|sources|sha256'][0] == E.sha256_file(a) and raw['REF|sources|sha256'][2] == E.sha256_file(c)
    E.part_ref(path)                                                                                  # idempotent (write-once, same values)


# ----------------------------------------------------------------------------------------------------------------------------------
# the guard and the runs
# ----------------------------------------------------------------------------------------------------------------------------------
def tiny_path_call():
    P = E.IM.PixelPSF(24, 24, 1.2)
    L = E.IM.Lattice(P, h=0.5, margin=2.0)
    y = P.synth(np.array([[12.0, 12.0]]), np.array([100.0]))
    return E.IM.sliding_path(P, y, np.ones((24, 24)), L, max_entries=3, tol=1e-6, passes=1, rtrunc=6.0, r_ex=1.2, fit_radius=4.0, order='dual')


def test_with_guard_replaces_the_sliding_path_during_the_call_restores_it_and_expects_exactly_one_path():
    orig = E.IM.sliding_path
    seen = []

    def fn():
        seen.append(E.IM.sliding_path is not orig)
        return tiny_path_call()
    res, tr = E.with_guard(fn)
    assert seen == [True] and E.IM.sliding_path is orig
    assert set(tr) == {'counts', 'barred', 'taus'} and tr['counts'][0] == len(res['taus']) and len(tr['barred']) == len(res['barred'])
    with pytest.raises(RuntimeError, match='expected one'):
        E.with_guard(lambda: (tiny_path_call(), tiny_path_call()))
    assert E.IM.sliding_path is orig
    with pytest.raises(RuntimeError, match='expected one'):
        E.with_guard(lambda: 0)
    assert E.IM.sliding_path is orig

    def boom():
        raise ValueError('x')
    with pytest.raises(ValueError):
        E.with_guard(boom)
    assert E.IM.sliding_path is orig                                                                  # restored when the run raises


def small_patch(stamp):
    stars = [(14.0, 14.0, 6.0e4), (31.5, 16.2, 3.0e4), (20.3, 33.1, 2.0e4)]
    dn = planted_frame(stamp, stars, seed=11, n=48)
    truth = dict(px=np.array([s[0] for s in stars]), py=np.array([s[1] for s in stars]), m606=A.counts_to_mag(np.array([s[2] for s in stars])))
    return dict(dn=dn, psf=stamp, gain=GAIN, dark=DARK), truth


def test_run_iss_runs_the_frozen_functions_under_the_guard_and_packs_the_scores_of_the_two_readouts(monkeypatch):
    st = lumpy_stamp(1.4, 1.4, 0.0)
    patch, truth = small_patch(st)
    t = dict(x=truth['px'], y=truth['py'], mag=truth['m606'], inner=np.ones(3, bool))
    calls = []
    orig = E.X.run_path

    def spy(*a, **k):
        calls.append(('iso', E.IM.sliding_path is not E.CF.ORIG_PATH))
        return orig(*a, **k)
    monkeypatch.setattr(E.X, 'run_path', spy)
    arr = E.run_iss('iso', patch, t, (3.0, 1.2, 2.0))
    assert calls == [('iso', True)]
    for k in ('k2', 'k0'):
        assert arr[k + '|tot'].shape == (7,) and arr[k + '|blk22'].shape == (25, 4)
    assert arr['det'].shape[1] == 3 and arr['secs'].shape == (2,) and arr['secs'][1] == arr['trace|counts'][0]
    assert arr['bic'][2] == 0 and arr['trace|barred'].shape[1] == 3
    # the two stars brighter than 22 are found; the catalogue at kappa 2 is the stored one
    assert arr['k2|tot'][0] == 3 and arr['k2|tot'][1] >= 2
    fit = E.PG.fit_adequate(st)
    arr2 = E.run_iss('adq', patch, t, (3.0, 1.2, 2.0), fit)
    assert arr2['k2|tot'][0] == 3 and arr2['det'].shape[1] == 3
    if fit['kind'] == 'isotropic':
        assert np.array_equal(arr2['k2|blk22'], arr['k2|blk22'])                                        # the same program on the same input


def test_run_baseline_scores_the_catalogue_with_the_frozen_scorer_and_the_crowdsource_default_filter(monkeypatch):
    st = lumpy_stamp()
    patch, truth = small_patch(st)
    t = dict(x=truth['px'], y=truth['py'], mag=truth['m606'], inner=np.ones(3, bool))
    cat = dict(x=np.array([14.0, 31.5, 5.0]), y=np.array([14.0, 16.2, 5.0]), flux=np.array([6.0e4, 3.0e4, -5.0]))
    monkeypatch.setattr(E.H2, 'run_crowdsource_hp2', lambda dn, psf, gain, dark, spec: (cat, 1.5))
    monkeypatch.setattr(E.H2, 'run_photutils_hp2', lambda dn, psf, gain, dark, spec: (cat, 2.5))
    for meth, secs in (('cs_hp', 1.5), ('phu_hp', 2.5)):
        arr = E.run_baseline(meth, patch, 300.0, t, (None, None))
        assert arr['secs'][0] == secs and arr['det'].shape == (3, 3)
        want = E.X._pack('k2', E.X.score_dets((cat['x'], cat['y'], cat['flux']), t))
        assert np.array_equal(arr['k2|blk22'], want['k2|blk22']) and np.array_equal(arr['k2|tot'], want['k2|tot'])
    full = dict(x=np.array([14.0, 31.5, 5.0, 40.0]), y=np.array([14.0, 16.2, 5.0, 40.0]), flux=np.array([6.0e4, 3.0e4, -5.0, 100.0]), dflux=np.ones(4), qf=np.ones(4), fracflux=np.ones(4),
                spread=np.array([0.0, 0.0, 0.0, np.nan]))
    seen = {}

    def fake_cs(dn, sky, psf, gain, dark):
        seen['args'] = (sky, gain, dark)
        return full, 1.0
    monkeypatch.setattr(E.HP.AB, 'run_crowdsource', fake_cs)
    arr = E.run_baseline('cs_def', patch, 123.0, t, (None, None))
    assert seen['args'] == (123.0, GAIN, DARK)
    assert arr['det'].shape == (2, 3) and np.array_equal(arr['det'][:, 2], [6.0e4, 3.0e4])           # the default filter: positive flux, a finite spread (NaN dropped), as `exp_astro_m2.cs_dets`
    with pytest.raises(ValueError):
        E.run_baseline('nope', patch, 1.0, t, (None, None))


def test_part_run_runs_each_method_once_per_patch_with_its_stamp_skips_the_nan_stamp_and_the_isotropic_adequate_kernel(monkeypatch, tmp_path):
    st, wrong, cache = planted_world(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    monkeypatch.setattr(E, 'frozen_config', lambda: (3.0, 1.2, 2.0))
    monkeypatch.setattr(E, 'round2_specs', lambda: ('cs-spec', 'phu-spec'))
    E.part_stamps(path, B=4)
    calls = []

    def fake_iss(kind, patch, t, cfg, fit=None):
        calls.append((kind, float(patch['psf'][25, 25]), cfg))
        blk, tot = fake_scores(len(calls))
        return {'k2|tot': tot, 'k2|blk22': blk, 'k0|tot': tot, 'k0|blk22': blk, 'det': np.zeros((2, 3)), 'secs': np.array([5.0, 100.0]), 'bic': np.array([3, 30, 0]),
                'trace|counts': np.array([100, 0, 100, 1, 0, 101, 0, 0]), 'trace|barred': np.zeros((0, 3))}

    def fake_base(meth, patch, sky_mean, t, specs):
        calls.append((meth, specs, sky_mean))
        blk, tot = fake_scores(len(calls))
        return {'k2|tot': tot, 'k2|blk22': blk, 'det': np.zeros((4, 3)), 'secs': np.array([2.0])}
    monkeypatch.setattr(E, 'run_iss', fake_iss)
    monkeypatch.setattr(E, 'run_baseline', fake_base)
    kinds = {'isotropic': dict(kind='isotropic', rms_isotropic=1e-3, rms=1e-3), 'general': dict(kind='general', rms_isotropic=5e-3, rms=1e-3)}
    which = ['isotropic', 'general']
    monkeypatch.setattr(E.PG, 'fit_adequate', lambda stamp: kinds[which[0]])
    E.part_run('patch', ['iso', 'adq', 'cs_hp', 'phu_hp', 'cs_def'], ['m2'], path)
    raw = E.load_raw(path)
    names = [c[0] for c in calls]
    assert names == ['iso', 'cs_hp', 'phu_hp', 'cs_def'] * 2                                          # adq is 'same' (isotropic): not run; per patch, in the order given
    stamp0 = raw['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))]
    assert calls[0][1] == float(stamp0[25, 25]) and calls[0][2] == (3.0, 1.2, 2.0)                    # the in-situ stamp, not the psField one, and the frozen configuration
    assert calls[1][1:] == (('cs-spec', 'phu-spec'), float(cache[('m2', (0, 0))][1].mean()))
    k = 'RUN|patch|%s|' % E.pkey('m2', (0, 0))
    assert raw[k + 'adq|same'][0] == 1 and str(raw[k + 'adq|kernel'][0]) == 'isotropic' and k + 'adq|k2|tot' not in raw
    assert raw[k + 'iso|same'][0] == 0 and str(raw[k + 'iso|kernel'][0]) == 'isotropic' and raw[k + 'iso|wall'].shape == (1,)
    assert raw[k + 'cs_def|k2|blk22'].shape == (25, 4) and raw[k + 'cs_def|wall'].shape == (1,)
    n = len(calls)
    E.part_run('patch', ['iso', 'adq', 'cs_hp', 'phu_hp', 'cs_def'], ['m2'], path)                    # nothing is run again
    assert len(calls) == n
    # the general kernel is run (adq) when the rule rejects the isotropic one
    which[0] = 'general'
    E.part_run('frame', ['adq'], ['m2'], path)
    assert [c[0] for c in calls[n:]] == ['adq', 'adq']
    k2 = 'RUN|frame|%s|' % E.pkey('m2', (1000, 0))
    assert raw.get(k2 + 'adq|same') is None
    raw2 = E.load_raw(path)
    assert raw2[k2 + 'adq|same'][0] == 0 and str(raw2[k2 + 'adq|kernel'][0]) == 'general' and list(raw2[k2 + 'adq|rms']) == [5e-3, 1e-3]
    # a patch without an in-situ stamp is skipped
    bad = dict(raw2)
    bad['ST|patch|%s|stamp' % E.pkey('m2', (1000, 0))] = np.full((51, 51), np.nan)
    p2 = str(tmp_path / 'bad.npz')
    np.savez_compressed(p2, **{k: v for k, v in bad.items() if not k.startswith('RUN|')})
    calls.clear()
    E.part_run('patch', ['iso'], ['m2'], p2)
    assert len(calls) == 1
    E.part_run('patch', ['iso'], ['m53'], p2)                                                          # no patches of that set among the targets
    assert len(calls) == 1


# ----------------------------------------------------------------------------------------------------------------------------------
# Task B: crowdsource options on the development patches
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_declared_set_of_variants_is_the_base_and_one_option_at_a_time():
    v = E.cs3_variants()
    assert len(v) == 24 and v[0][0] == 'cs3_base' and v[0][3] == E.CS_BASE and len({t for t, _, _, _ in v}) == 24
    levels = dict(E.CS_LEVELS)
    assert levels['blendthreshu'] == (1, 3, 4) and levels['titer_thresh'] == (1, 3, 5) and levels['psfderiv'] == (False,) and levels['maxiter'] == (30, 40) and levels['miniter'] == (4, 12)
    assert levels['fewstars'] == (5,) and levels['psfvalsharpcutfac'] == (0.5, 0.9) and levels['psfsharpsat'] == (0.5, 1.0) and levels['thr'] == (2.5, 4.0) and levels['blist'] == (5, 15)
    assert levels['weights'] == ('flat', 'smooth3', 'satmask')
    for tag, fam, lv, spec in v[1:]:
        diff = [k for k in spec if spec[k] != E.CS_BASE[k]]
        assert diff == [fam] and spec[fam] == lv and tag == E.cs3_tag(fam, lv)
    assert E.CS_BASE == dict(sb=2.0, thr=3.0, miniter=8, maxiter=20, fewstars=20, psfderiv=True, blendthreshu=2, titer_thresh=2, psfvalsharpcutfac=0.7, psfsharpsat=0.7, weights='iss', blist=0)
    assert E.cs3_tag('thr', 2.5) == 'cs3_thr2p5' and E.cs3_tag('psfderiv', False) == 'cs3_psfderivFalse'
    h2 = next(d for d in E.H2.CS2 if d['tag'] == E.CS_BASE_TAG)
    assert (h2['sb'], h2['thr'], E.H2.PASSES[h2['passes']]) == (E.CS_BASE['sb'], E.CS_BASE['thr'], dict(miniter=8, maxiter=20, fewstars=20))        # the base is the round-2 selection
    assert E.COMBO_GAIN == 0.003 and E.FINDING_GAIN == 0.010


def test_the_weights_of_the_variants_are_the_iss_weights_a_constant_a_wider_smoothing_and_a_saturation_mask():
    st = lumpy_stamp()
    dn = planted_frame(st, STARS, seed=2)
    dn[40, 40] = E.SAT_DN + 1
    iss = E.HP.iss_weights(dn, GAIN, DARK)
    assert np.array_equal(E.cs3_weights(dn, GAIN, DARK, 'iss'), iss)
    flat = E.cs3_weights(dn, GAIN, DARK, 'flat')
    assert flat.shape == dn.shape and np.all(flat == np.median(iss))
    sm = E.cs3_weights(dn, GAIN, DARK, 'smooth3')
    exp = 1.0 / np.sqrt(np.maximum(E.gaussian_filter(dn, 3.0, mode='reflect'), np.percentile(dn, 1)) / GAIN + DARK)
    assert np.allclose(sm, exp) and not np.allclose(sm, iss)
    mk = E.cs3_weights(dn, GAIN, DARK, 'satmask')
    assert mk[40, 40] == 0.0 and (mk == 0).sum() == (dn > E.SAT_DN).sum() and np.array_equal(mk[mk > 0], iss[mk > 0])
    dn[40, 40] = E.SAT_DN                                                                                  # at the limit: not masked
    assert (E.cs3_weights(dn, GAIN, DARK, 'satmask') == 0).sum() == 0
    with pytest.raises(ValueError):
        E.cs3_weights(dn, GAIN, DARK, 'nope')


def test_the_bright_star_list_is_the_brightest_maxima_with_the_first_axis_first():
    st = lumpy_stamp()
    dn = planted_frame(st, [(20.0, 60.0, 1.0e5), (70.0, 25.0, 6.0e4), (45.0, 45.0, 8.0e4)], seed=4)
    bx, by, bm = E.cs3_blist(dn, GAIN, DARK, 2)
    assert len(bx) == len(by) == len(bm) == 2 and bx.dtype == np.float32
    assert (abs(bx[0] - 60) <= 1 and abs(by[0] - 20) <= 1) and (abs(bx[1] - 45) <= 1 and abs(by[1] - 45) <= 1)          # crowdsource's x runs along the rows; brightest first
    assert len(E.cs3_blist(dn, GAIN, DARK, 15)[0]) > 2


def fake_dev_raw(scores, truth_total=(100, 0, 0, 0, 200, 0, 0)):
    """A raw with the development scores of the variants: scores = {tag: F}, built from totals with true22 = 100 and completeness = purity = F."""
    raw = {}
    for tag, F in scores.items():
        for p in E.cs_dev_patches():
            n = int(round(F * 100000))
            raw['CSD|%s|x%d_y%d|tot' % (tag, p[0], p[1])] = np.array([100000, n, 100000, n, 200000, n, 100000])
    return raw


def test_the_development_f_pools_the_counts_of_the_two_patches_and_is_nan_without_both():
    raw = {}
    pa, pb = E.cs_dev_patches()
    raw['CSD|t|x%d_y%d|tot' % pa] = np.array([100, 40, 80, 40, 150, 40, 90])
    assert np.isnan(E.dev_f(raw, 't'))
    raw['CSD|t|x%d_y%d|tot' % pb] = np.array([300, 90, 220, 100, 400, 90, 250])
    assert E.dev_f(raw, 't') == pytest.approx(E.X.f_measure(np.array([400, 130, 300, 140, 550, 130, 340]))[0])
    assert not np.isnan(E.dev_f(raw, 't')) and np.isnan(E.dev_f(raw, 'other'))
    assert E.dev_f({'CSTEST|t|x%d_y%d|tot' % pa: raw['CSD|t|x%d_y%d|tot' % pa], 'CSTEST|t|x%d_y%d|tot' % pb: raw['CSD|t|x%d_y%d|tot' % pb]}, 't', 'CSTEST') == E.dev_f(raw, 't')


def test_the_combination_takes_the_best_level_of_every_option_that_gains_three_thousandths():
    sc = {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 1): 0.442, E.cs3_tag('blendthreshu', 3): 0.4435, E.cs3_tag('blendthreshu', 4): 0.450, E.cs3_tag('titer_thresh', 1): 0.4429,
          E.cs3_tag('maxiter', 30): 0.443, E.cs3_tag('maxiter', 40): 0.443, E.cs3_tag('thr', 2.5): 0.430, E.cs3_tag('weights', 'flat'): 0.4431}
    spec, picked = E.combo_spec(fake_dev_raw(sc))
    assert spec['blendthreshu'] == 4 and spec['maxiter'] == 30 and spec['weights'] == 'flat'          # 0.443 is the base + 0.003 (counts as a gain); a tie goes to the first level
    assert spec['titer_thresh'] == 2 and spec['thr'] == 3.0 and spec['psfderiv'] is True and E.CS_BASE['blendthreshu'] == 2
    assert picked == ['blendthreshu=4', 'maxiter=30', 'weights=flat']
    sc2 = {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 3): 0.4429}
    spec2, picked2 = E.combo_spec(fake_dev_raw(sc2))
    assert picked2 == [] and spec2 == E.CS_BASE


def run_part_cs_dev_with(monkeypatch, tmp_path, scores, fname='insitu_raw.npz'):
    path = str(tmp_path / fname)
    monkeypatch.setattr(E, 'check_hashes', lambda path=None: None)
    monkeypatch.setattr(E.X, 'truth_arrays', lambda ds='m2': 'truth')
    ran = []

    def fake_variant(tag, spec, patches, truth, prefix, path=None):
        ran.append((tag, prefix))
        F = scores.get(tag, 0.3)
        E.save_with_meta({'%s|%s|x%d_y%d|tot' % (prefix, tag, p[0], p[1]): np.array([100000, int(round(F * 100000)), 100000, int(round(F * 100000)), 200000, 0, 100000]) for p in patches}, path)
    monkeypatch.setattr(E, 'part_cs_variant', fake_variant)
    return path, ran


def test_stage_one_runs_every_declared_variant_and_stage_two_the_combination_and_the_selection(monkeypatch, tmp_path):
    sc = {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 4): 0.4504, E.cs3_tag('maxiter', 30): 0.4431}
    path, ran = run_part_cs_dev_with(monkeypatch, tmp_path, dict(sc, cs3_combo=0.4502))
    E.part_cs_dev(1, path)
    assert [t for t, _ in ran] == [t for t, _, _, _ in E.cs3_variants()] and all(p == 'CSD' for _, p in ran)
    ran.clear()
    E.part_cs_dev(2, path)
    raw = E.load_raw(path)
    assert [t for t, _ in ran] == ['cs3_combo'] and list(raw['CSD|SEL|combo']) == ['blendthreshu=4', 'maxiter=30']
    assert str(raw['CSD|SEL|best'][0]) == E.cs3_tag('blendthreshu', 4) and raw['CSD|SEL|best_F'][0] == pytest.approx(0.4504, abs=1e-6) and raw['CSD|SEL|base_F'][0] == pytest.approx(0.440, abs=1e-6)
    assert raw['CSD|SEL|finding'][0] == 1                                                                 # 0.0104 > 0.010
    with pytest.raises(RuntimeError, match='never made twice'):
        E.part_cs_dev(2, path)


def test_a_gain_of_exactly_the_threshold_is_not_a_finding_a_tie_goes_to_the_base_and_no_combination_is_run_without_gains(monkeypatch, tmp_path):
    path, ran = run_part_cs_dev_with(monkeypatch, tmp_path, {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 4): 0.450, E.cs3_tag('maxiter', 40): 0.450})
    E.part_cs_dev(1, path)
    ran.clear()
    E.part_cs_dev(2, path)
    raw = E.load_raw(path)
    assert str(raw['CSD|SEL|best'][0]) == E.cs3_tag('blendthreshu', 4)       # two variants tie at 0.450: the first in the declared order
    assert raw['CSD|SEL|finding'][0] == 0                                                                  # 0.450 - 0.440 = 0.010 is not MORE than 0.010
    path3, ran3 = run_part_cs_dev_with(monkeypatch, tmp_path, {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 4): 0.440}, 'raw3.npz')
    E.part_cs_dev(1, path3)
    ran3.clear()
    E.part_cs_dev(2, path3)
    raw3 = E.load_raw(path3)
    assert ran3 == [] and list(raw3['CSD|SEL|combo']) == ['none'] and str(raw3['CSD|SEL|best'][0]) == 'cs3_base' and raw3['CSD|SEL|finding'][0] == 0


def test_the_test_patches_are_scored_only_for_a_finding_and_only_the_best_and_the_base(monkeypatch, tmp_path):
    path, ran = run_part_cs_dev_with(monkeypatch, tmp_path, {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 4): 0.4504})
    with pytest.raises(RuntimeError, match='not a finding'):
        E.part_cs_test(path)
    E.part_cs_dev(1, path)
    E.part_cs_dev(2, path)
    ran.clear()
    E.part_cs_test(path)
    assert ran == [('cs3_base', 'CSTEST'), (E.cs3_tag('blendthreshu', 4), 'CSTEST')]
    path2, ran2 = run_part_cs_dev_with(monkeypatch, tmp_path, {'cs3_base': 0.440}, 'raw2.npz')
    p2 = path2
    E.part_cs_dev(1, p2)
    E.part_cs_dev(2, p2)
    ran2.clear()
    with pytest.raises(RuntimeError, match='not a finding'):
        E.part_cs_test(p2)
    assert ran2 == []


def test_part_cs_variant_scores_a_variant_per_patch_records_a_failure_and_does_not_run_again(monkeypatch, tmp_path):
    path = str(tmp_path / 'insitu_raw.npz')
    st = lumpy_stamp()
    dn = planted_frame(st, STARS, seed=1)
    truth = dict(px=np.array([s[0] for s in STARS]) + 110, py=np.array([s[1] for s in STARS]) + 630, m606=A.counts_to_mag(np.array([s[2] for s in STARS])))
    monkeypatch.setattr(E.HP, 'load_arrays', lambda ds, x0, y0: dict(dn=dn, psf51=st, gain=GAIN, dark=DARK, sky=SKY))
    calls = []

    def fake_run(dn_, psf, gain, dark, spec):
        calls.append(spec['tag'])
        if spec['tag'] == 'bad' and calls.count('bad') == 1:
            raise ValueError('boom')
        return dict(x=np.array([s[0] for s in STARS]), y=np.array([s[1] for s in STARS]), flux=np.array([s[2] for s in STARS])), 3.5
    monkeypatch.setattr(E, 'run_cs3', fake_run)
    pats = E.cs_dev_patches()
    E.part_cs_variant('good', dict(tag='good'), pats, dict(truth, **{}), 'CSD', path)
    raw = E.load_raw(path)
    for p in pats:
        k = 'CSD|good|x%d_y%d|' % p
        t = A.truth_for_patch(truth, p[0], p[1])
        want = E.X._pack('', E.X.score_dets((np.array([s[0] for s in STARS]), np.array([s[1] for s in STARS]), np.array([s[2] for s in STARS])), t))
        assert np.array_equal(raw[k + 'tot'], want['|tot']) and np.array_equal(raw[k + 'blk22'], want['|blk22']) and raw[k + 'secs'][0] == 3.5 and raw[k + 'ndet'][0] == len(STARS)
    n = len(calls)
    E.part_cs_variant('good', dict(tag='good'), pats, truth, 'CSD', path)
    assert len(calls) == n                                                                                 # not run again
    E.part_cs_variant('bad', dict(tag='bad'), pats, truth, 'CSD', path)
    raw = E.load_raw(path)
    assert 'CSD|bad|x%d_y%d|failed' % pats[0] in raw and 'ValueError: boom' in str(raw['CSD|bad|x%d_y%d|failed' % pats[0]][0]) and 'CSD|bad|x%d_y%d|tot' % pats[1] in raw
    assert np.isnan(E.dev_f(raw, 'bad'))
    calls.clear()
    E.part_cs_variant('bad', dict(tag='bad'), pats, truth, 'CSD', path)                                   # a failure is recorded, not retried
    assert calls == []


def test_the_crowdsource_variant_runs_on_a_planted_patch_when_the_package_is_installed():
    pytest.importorskip('crowdsource')
    st = lumpy_stamp(1.4, 1.4, 0.0)
    patch, truth = small_patch(st)
    cat, secs = E.run_cs3(patch['dn'], st, GAIN, DARK, dict(E.CS_BASE))
    assert len(cat['x']) >= 3 and secs >= 0
    near = [np.min(np.hypot(cat['x'] - tx, cat['y'] - ty)) for tx, ty in zip(truth['px'], truth['py'])]
    assert max(near) < 1.0                                                                                  # column and row not exchanged
    cat2, _ = E.run_cs3(patch['dn'], st, GAIN, DARK, dict(E.CS_BASE, blist=2, psfderiv=False, weights='satmask'))
    assert len(cat2['x']) >= 3


# ----------------------------------------------------------------------------------------------------------------------------------
# the report on a planted raw
# ----------------------------------------------------------------------------------------------------------------------------------
SETS4 = (('m2', (0, 0)), ('m2', (1000, 0)), ('m53', (0, 0)), ('m92', (0, 0)))


def blocks_for(F, seed, noise=0.03):
    """25 blocks [true22, matched22, det22, matched det22] whose pooled completeness and purity are about F."""
    rng = np.random.default_rng(seed)
    t, d = rng.integers(8, 13, 25), rng.integers(8, 13, 25)
    c, p = np.clip(F + noise * rng.standard_normal(25), 0, 1), np.clip(F + noise * rng.standard_normal(25), 0, 1)
    return np.c_[t, np.round(c * t), d, np.round(p * d)].astype(float)


def tot_of(blk):
    return np.array([blk[:, 0].sum(), blk[:, 1].sum(), blk[:, 2].sum(), blk[:, 3].sum(), 400, 200, 300])


LEVEL = {   # F per (set, stamp, method): the M2 control keeps its lead; M53: the loss reverses; M92: it closes to a tie
    ('m2', 'psf'): dict(iso=0.60, adq=0.60, cs_hp=0.40, phu_hp=0.25, cs_def=0.30), ('m2', 'patch'): dict(iso=0.58, adq=0.58, cs_hp=0.41, phu_hp=0.25, cs_def=0.30),
    ('m2', 'frame'): dict(iso=0.59, adq=0.59, cs_hp=0.40, phu_hp=0.25, cs_def=0.30),
    ('m53', 'psf'): dict(iso=0.20, adq=0.30, cs_hp=0.40, phu_hp=0.15, cs_def=0.25), ('m53', 'patch'): dict(iso=0.50, adq=0.55, cs_hp=0.40, phu_hp=0.15, cs_def=0.25),
    ('m53', 'frame'): dict(iso=0.40, adq=0.42, cs_hp=0.40, phu_hp=0.15, cs_def=0.25),
    ('m92', 'psf'): dict(iso=0.15, adq=0.25, cs_hp=0.40, phu_hp=0.20, cs_def=0.30), ('m92', 'patch'): dict(iso=0.40, adq=0.40, cs_hp=0.40, phu_hp=0.20, cs_def=0.30),
    ('m92', 'frame'): dict(iso=0.40, adq=0.30, cs_hp=0.40, phu_hp=0.20, cs_def=0.30),
    ('m2', 'scr'): dict(iso=0.60, adq=0.60, cs_hp=0.30, phu_hp=0.25, cs_def=0.28), ('m53', 'scr'): dict(iso=0.20, adq=0.30, cs_hp=0.30, phu_hp=0.15, cs_def=0.22),
    ('m92', 'scr'): dict(iso=0.15, adq=0.25, cs_hp=0.30, phu_hp=0.20, cs_def=0.28)}


@pytest.fixture(scope='module')
def planted_raw(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    st, wrong, cache = planted_world(mp, SETS4)
    path = str(tmp_path_factory.mktemp('r') / 'insitu_raw.npz')
    E.part_stamps(path, B=4)
    arrays = {}
    seed = 0
    for name, pxy in SETS4:
        for stamp in ('psf', 'patch', 'frame', 'scr'):
            for meth in E.METHODS:
                seed += 1
                F = LEVEL[(name, stamp)][meth]
                blk = blocks_for(F, seed if not (name == 'm92' and stamp == 'patch' and meth in ('iso', 'cs_hp')) else 999)         # the m92 tie: identical blocks for iso and cs_hp
                pre = 'REF|%s|%s' % (E.pkey(name, pxy), meth) if stamp == 'psf' else 'RUN|%s|%s|%s' % (stamp, E.pkey(name, pxy), meth)
                if stamp != 'psf' and meth == 'adq' and name == 'm2':
                    arrays[pre + '|same'] = np.array([1])
                    arrays[pre + '|kernel'] = np.array(['isotropic'])
                    arrays[pre + '|rms'] = np.array([1e-3, 1e-3])
                    continue
                for k in (('k2', 'k0') if meth in ('iso', 'adq') else ('k2',)):
                    arrays['%s|%s|blk22' % (pre, k)] = blk
                    arrays['%s|%s|tot' % (pre, k)] = tot_of(blk)
                arrays[pre + '|secs'] = np.array([50.0 + seed, 1000.0 + seed]) if meth in ('iso', 'adq') else np.array([2.0])
                arrays[pre + '|det'] = np.zeros((300 + seed, 3))
                if meth in ('iso', 'adq'):
                    arrays[pre + '|bic'] = np.array([5, 50, 0])
                    if stamp == 'psf':
                        arrays[pre + '|barred'] = np.array([seed % 3])
                    else:
                        arrays[pre + '|trace|counts'] = np.array([1000 + seed, 10, 1000, 2, 0, 1001, seed % 3, 0])
                        arrays[pre + '|kernel'] = np.array(['general' if meth == 'adq' else 'isotropic'])
                        arrays[pre + '|same'] = np.array([0])
                        arrays[pre + '|rms'] = np.array([5e-3, 1e-3])
                if stamp != 'psf':
                    arrays[pre + '|wall'] = np.array([60.0 + seed])
    E.save_with_meta(arrays, path)
    # Task B
    cs = {}
    for k, (tag, F) in enumerate([('cs3_base', 0.4480), ('cs3_blendthreshu1', 0.4400), ('cs3_blendthreshu4', 0.4520), ('cs3_maxiter30', 0.4485)]):
        for p in E.cs_dev_patches():
            n = int(round(F * 100000))
            cs['CSD|%s|x%d_y%d|tot' % (tag, p[0], p[1])] = np.array([100000, n, 100000, n, 200000, n, 100000])
            cs['CSD|%s|x%d_y%d|secs' % (tag, p[0], p[1])] = np.array([3.0])
    cs['CSD|cs3_psfderivFalse|x110_y630|failed'] = np.array(['ValueError: boom'])
    cs.update({'CSD|SEL|best': np.array(['cs3_blendthreshu4']), 'CSD|SEL|best_F': np.array([0.4520]), 'CSD|SEL|base_F': np.array([0.4480]), 'CSD|SEL|combo': np.array(['none']), 'CSD|SEL|finding': np.array([0])})
    E.save_with_meta(cs, path)
    raw = E.load_raw(path)
    yield raw, path
    mp.undo()


def table_rows(text, start):
    """The rows of the first markdown table after the line that starts with `start`."""
    lines = text.split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith(start))
    j = i + 1
    while not lines[j].startswith('|'):
        j += 1
    rows = []
    while j < len(lines) and lines[j].startswith('|'):
        rows.append([c.strip() for c in lines[j].strip('|').split('|')])
        j += 1
    return rows


def test_the_report_has_every_section_and_is_deterministic(planted_raw):
    raw, path = planted_raw
    text = E.report(raw)
    assert text == E.report(E.load_raw(path)) and text.endswith('\n') and not text.endswith('\n\n') and '\n\n\n' not in text
    for h in ('# An in-situ PSF for every method', '## 1. The in-situ stamp', '### 1.a Scope: PATCH', '### 1.b Scope: FRAME', '### 1.c Against the stars selected with the truth', '### 1.d Recovery in the exact-stamp world',
              '## 2. F at mag < 22 with each stamp', '### 2.1 Pooled by group', '### 2.2 Per patch', '### 2.3 Paired differences', '### 2.4 The effect of the stamp', '### 2.5 The ISS-minus-crowdsource gap',
              '## 3. The answer, as the tables give it', '## 4. Run records', '## 5. Task B', '## 6. Limitations'):
        assert h in text, h
    assert 'POST HOC' in text and 'Code the runs used' in text and 'src/core/imaging_guard.py' in text


def test_the_pooled_f_cells_are_the_block_bootstrap_of_the_pooled_blocks_with_the_stratified_seed_zero_resamples(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    rows = table_rows(text, '### 2.1 Pooled by group')
    group, by = '', {}
    for r in rows[2:]:
        group = r[0] or group
        by[(group, r[1])] = r[2:]
    row = by[('M53 (NGC 5024)', E.LABEL['cs_hp'])]
    blk = raw['REF|%s|cs_hp|k2|blk22' % E.pkey('m53', (0, 0))]
    f, ci = E.X.boot_cpf(blk, E.X.boot_idx_strat([25]))[2]
    assert row[0] == '%s (1)' % E.X.fmt_iv(f, ci)
    row2 = by[('M2 test patches (control)', E.LABEL['iso'])]
    blk2 = np.concatenate([raw['REF|%s|iso|k2|blk22' % E.pkey('m2', p)] for p in ((0, 0), (1000, 0))])
    f2, ci2 = E.X.boot_cpf(blk2, E.X.boot_idx_strat([25, 25]))[2]
    assert row2[0] == '%s (2)' % E.X.fmt_iv(f2, ci2)
    assert 'Left out of a pool because the ISS path reached the entry cap (a defective run): none.' in text and 'Rows with no run (patches missing from the pool): none.' in text


def test_the_paired_differences_use_the_same_stamp_for_both_methods_and_the_gap_sentences_follow_the_verdicts(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    section = text.split('### 2.5 ')[1]
    m53 = section.split('**M53 (NGC 5024)**')[1].split('**M92 (NGC 6341)**')[0]
    assert 'ISS iso minus crowdsource+HP | in-situ stamp (patch stars) | 1 |' in m53 and 'the loss REVERSES (the ISS wins)' in m53
    m92 = section.split('**M92 (NGC 6341)**')[1]
    line = next(l for l in m92.split('\n') if l.startswith('| ISS iso minus crowdsource+HP | in-situ stamp (patch stars)'))
    assert 'the loss closes (the interval now contains 0)' in line and 'the change of the gap is not distinguishable from 0' not in line
    m2 = section.split('**M2 test patches (control)**')[1].split('**M53 (NGC 5024)**')[0]
    line = next(l for l in m2.split('\n') if l.startswith('| ISS iso minus crowdsource+HP | in-situ stamp (patch stars)'))
    assert 'the ISS still wins' in line
    # the paired difference of 2.3 is the bootstrap of the two rows' blocks
    d23 = table_rows(text, '**M53 (NGC 5024)**')
    row = next(r for r in d23 if r[0] == 'ISS iso minus crowdsource+HP')
    a, b = raw['RUN|patch|%s|iso|k2|blk22' % E.pkey('m53', (0, 0))], raw['RUN|patch|%s|cs_hp|k2|blk22' % E.pkey('m53', (0, 0))]
    df = E.X.boot_diff(a, b, E.X.boot_idx_strat([25]))[2]
    assert row[2] == '%s (1)' % E.diff_cell(df[0], df[1])
    # the effect of the stamp (2.4) is in-situ minus psField, paired
    e24 = table_rows(text.split('### 2.4 ')[1], '**M53 (NGC 5024)**')
    row = next(r for r in e24 if r[0] == E.LABEL['iso'])
    a0 = raw['REF|%s|iso|k2|blk22' % E.pkey('m53', (0, 0))]
    df = E.X.boot_diff(a, a0, E.X.boot_idx_strat([25]))[2]
    assert row[1] == '%s (1)' % E.diff_cell(df[0], df[1]) and 'win' in row[1]


def test_the_gap_sentences_cover_the_nine_verdict_pairs_and_the_three_changes():
    seen = set()
    for v0 in ('loss', 'tie', 'win'):
        for v1 in ('loss', 'tie', 'win'):
            for vd in ('loss', 'tie', 'win'):
                t = E.gap_sentence(v0, v1, vd)
                assert t not in seen
                seen.add(t)
    assert len(seen) == 27
    assert E.gap_sentence('loss', 'win', 'win').startswith('the loss REVERSES') and E.gap_sentence('loss', 'tie', 'tie').startswith('the loss closes')
    assert E.gap_sentence('win', 'loss', 'loss').startswith('a win becomes a LOSS') and 'AGAINST the ISS' in E.gap_sentence('win', 'loss', 'loss')
    assert E.gap_sentence('tie', 'tie', 'win').startswith('still a tie') and 'in the ISS\'s favour' in E.gap_sentence('tie', 'tie', 'win')


def test_the_answer_section_states_each_group_and_scope(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    ans = text.split('## 3. The answer')[1].split('## 4.')[0]
    lines = [l for l in ans.split('\n') if l.startswith('* ')]
    assert len(lines) == 3 * 2 * 2                                                                         # three groups, two scopes, ISS iso and adq against crowdsource+HP
    l = next(l for l in lines if l.startswith('* M53 (NGC 5024), patch stars stamp, ISS iso minus crowdsource+HP'))
    assert 'psField verdict loss, in-situ verdict win' in l and 'REVERSES' in l
    l = next(l for l in lines if l.startswith('* M2 test patches (control), frame stars stamp, ISS adq minus crowdsource+HP'))
    assert 'psField verdict win' in l


def test_a_stalled_iss_path_is_left_out_of_every_pool_and_named(planted_raw, tmp_path):
    raw, _ = planted_raw
    bad = dict(raw)
    k = 'RUN|patch|%s|iso|bic' % E.pkey('m2', (1000, 0))
    bad[k] = np.array([5, 2501, 1])
    text = E.report(bad)
    assert 'm2 (1000, 0), in-situ stamp (patch stars), ISS, isotropic kernel (frozen config)' in text.split('Left out of a pool')[1].split('\n')[0]
    assert 'm2 (1000, 0), in-situ stamp (patch stars), ISS, kernel from the stamp (rule)' in text.split('Left out of a pool')[1].split('\n')[0]       # the adequate row is the isotropic run
    assert 'Rows with no run (patches missing from the pool): none.' in text
    rows = table_rows(text, '### 2.1 Pooled by group')
    row = next(r for r in rows if r[1] == E.LABEL['iso'] and r[0].startswith('M2 test'))
    assert row[3].endswith('(1)') and row[2].endswith('(2)')
    blk = bad['RUN|patch|%s|iso|k2|blk22' % E.pkey('m2', (0, 0))]
    f, ci = E.X.boot_cpf(blk, E.X.boot_idx_strat([25]))[2]
    assert row[3] == '%s (1)' % E.X.fmt_iv(f, ci)
    runs = table_rows(text, '## 4. Run records')
    assert any(r[1] == 'm2' and r[2] == '(1000, 0)' and r[3] == E.LABEL['iso'] and r[0] == 'patch' and r[7] == 'STALLED' for r in runs)
    assert not any(r[7] == 'STALLED' for r in table_rows(E.report(raw), '## 4. Run records'))
    # the adequate row of a patch where the rule kept the isotropic kernel is the isotropic run: a stalled iso path stalls adq too
    k2 = 'RUN|patch|%s|adq|k2|blk22' % E.pkey('m2', (1000, 0))
    assert E.is_stalled(bad, 'patch', 'm2', (1000, 0), 'adq') and k2 not in bad and E.blk_of(bad, 'patch', 'm2', (1000, 0), 'adq') is not None
    assert not E.is_stalled(bad, 'patch', 'm2', (1000, 0), 'cs_hp') and not E.is_stalled(bad, 'psf', 'm2', (0, 0), 'iso')


def test_the_run_records_list_kernel_entries_guard_activations_and_seconds(planted_raw):
    raw, _ = planted_raw
    runs = table_rows(E.report(raw), '## 4. Run records')
    row = next(r for r in runs if r[0] == 'patch' and r[1] == 'm53' and r[3] == E.LABEL['adq'])
    pre = 'RUN|patch|%s|adq' % E.pkey('m53', (0, 0))
    assert row[4] == 'general' and row[5] == str(raw[pre + '|trace|counts'][0]) and row[6] == str(raw[pre + '|trace|counts'][6]) and row[8] == '5 of 50' and row[9] == '%.0f' % raw[pre + '|secs'][0]
    row = next(r for r in runs if r[0] == 'patch' and r[1] == 'm2' and r[3] == E.LABEL['adq'])
    assert row[4] == 'isotropic (same as iso)' and row[5] == '-'
    row = next(r for r in runs if r[0] == 'frame' and r[1] == 'm92' and r[3] == E.LABEL['cs_hp'])
    assert row[5] == str(raw['RUN|frame|%s|cs_hp|det' % E.pkey('m92', (0, 0))].shape[0]) and row[4] == '-'


def test_the_stamp_tables_print_the_stored_numbers(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    rows = table_rows(text, '### 1.a Scope: PATCH')
    k = 'ST|patch|%s|' % E.pkey('m53', (0, 0))
    row = next(r for r in rows if r[0] == 'm53')
    v, ci = raw[k + 'vs_psf'], raw[k + 'vs_psf_ci']
    assert row[1] == '(0, 0)' and row[2] == '%d / %d' % (raw[k + 'n_sel'][0], raw[k + 'n_used'][0]) and row[3] == E.pct_ci(v[0], ci[0, 0])
    assert row[4] == '%s [%s, %s]' % (E.num(v[1]), E.num(ci[0, 1, 0]), E.num(ci[0, 1, 1])) and row[5] == E.num(v[2]) and row[6] == '%s [%s, %s]' % (E.num(100 * v[3], 1), E.num(100 * ci[0, 2, 0], 1), E.num(100 * ci[0, 2, 1], 1))
    assert row[7] == '%s / %s' % (E.num(raw[k + 'shape'][0], 3), E.num(raw['ST|psf|%s|shape' % E.pkey('m53', (0, 0))][0], 3))
    assert 'Mean over the 2 m2 patches with a stamp: central pixel' in text and 'Mean over the 1 m92 patches with a stamp' in text
    c = table_rows(text, '### 1.c Against the stars')
    row = next(r for r in c if r[0] == 'm2' and r[1] == '(1000, 0)')
    assert row[2].startswith('%d / ' % len(STARS)) and row[9].startswith('%d of %d' % (len(STARS), len(STARS)))
    d = table_rows(text, '### 1.d Recovery')
    row = next(r for r in d if r[0] == 'm92')
    h = raw['ST|patch|%s|held' % E.pkey('m92', (0, 0))]
    assert row[8] == str(int(np.isfinite(h).all(1).sum())) and row[2] == str(raw['ST|patch|%s|world|n_sel' % E.pkey('m92', (0, 0))][0])
    frame = table_rows(text, 'Pooled by frame')
    assert [r[0] for r in frame[2:]] == ['m2', 'm53', 'm92']
    assert frame[2][1] == '%d / %d' % (raw['ST|frame|m2|all|truth_n'][0], raw['ST|frame|m2|all|truth_n_used'][0])


def test_task_b_section_lists_the_variants_the_failure_the_selection_and_whether_it_is_a_finding(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    rows = table_rows(text, 'Base: the round-2 selection')
    base = next(r for r in rows if r[0] == 'cs3_base')
    assert base[1] == 'base' and base[4] == '0.448' and base[5] == '+0.000'
    assert next(r for r in rows if r[0] == 'cs3_blendthreshu4')[5] == '+0.004' and next(r for r in rows if r[0] == 'cs3_blendthreshu1')[5] == '-0.008'
    assert next(r for r in rows if r[0] == 'cs3_blendthreshu4')[1] == 'blendthreshu = 4'
    assert 'FAILED: ValueError: boom' in next(r for r in rows if r[0] == 'cs3_psfderivFalse')[6]
    assert 'Selected (largest pooled development F; a tie to the base): **cs3_blendthreshu4**, F 0.452 against the round-2 selection 0.448 (difference +0.004)' in text
    assert 'Not a finding' in text and 'the test patches were not scored' in text
    bad = dict(raw)
    bad['CSD|SEL|finding'] = np.array([1])
    for tag in ('cs3_base', 'cs3_blendthreshu4'):
        for p in E.X.TESTS:
            bad['CSTEST|%s|x%d_y%d|tot' % (tag, p[0], p[1])] = np.array([100, 40 if tag == 'cs3_base' else 50, 100, 40 if tag == 'cs3_base' else 50, 0, 0, 0])
    t2 = E.report(bad)
    assert '**A finding**' in t2 and 'pooled counts | 0.400 | 0.500' in t2


def test_an_empty_raw_says_what_is_missing_and_a_partial_raw_has_dashes(planted_raw):
    t = E.report({})
    assert '(no stamps in the raw yet)' in t and '(no runs in the raw yet)' in t and '(not run yet)' in t and '## 4. Run records' in t
    raw, _ = planted_raw
    part = {k: v for k, v in raw.items() if not k.startswith('RUN|frame|') and not k.startswith('RUN|scr|') and not k.startswith('CSD|')}
    t2 = E.report(part)
    assert '(not run yet)' in t2
    rows = table_rows(t2, '### 2.1 Pooled by group')
    assert len(rows[0]) == 4 and rows[0][3].startswith('in-situ stamp (patch stars)')                    # no frame column without frame runs
    one = {k: v for k, v in part.items() if not k.startswith('RUN|patch|%s|cs_hp' % E.pkey('m53', (0, 0)))}
    t3 = E.report(one)
    rows = table_rows(t3, '### 2.2 Per patch')
    assert rows[0][:2] == ['set', 'patch']
    assert 'Rows with no run (patches missing from the pool): ' in t3 and 'in-situ stamp (patch stars), crowdsource + HP (round 2): 1' in t3 and 'Rows with no run (patches missing from the pool): none.' not in t3


def test_the_command_line_runs_one_part_per_flag_and_regenerates_the_report_byte_for_byte(planted_raw, monkeypatch, tmp_path):
    raw, path = planted_raw
    out = tmp_path / 'out'
    monkeypatch.setattr(E, 'OUT', str(out))
    monkeypatch.setattr(E, 'RAW', path)
    monkeypatch.setattr(E, 'MD', str(out / 'insitu.md'))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_insitu.py', '--from-raw'])
    E.main()
    first = (out / 'insitu.md').read_bytes()
    assert first == E.report(E.load_raw(path)).encode()
    E.main()
    assert (out / 'insitu.md').read_bytes() == first
    calls = []
    for flag, fn, args in (('--stamps', 'part_stamps', ()), ('--scr', 'part_scr', ()), ('--ref', 'part_ref', ()), ('--cs-test', 'part_cs_test', ())):
        monkeypatch.setattr(E, fn, lambda *a, _f=fn: calls.append(_f))
        monkeypatch.setattr(sys, 'argv', ['exp_astro_insitu.py', flag])
        E.main()
    monkeypatch.setattr(E, 'part_cs_dev', lambda stage: calls.append(('cs_dev', stage)))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_insitu.py', '--cs-dev', '--stage', '2'])
    E.main()
    monkeypatch.setattr(E, 'part_run', lambda scope, methods, sets: calls.append(('run', scope, methods, sets)))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_insitu.py', '--run', '--scope', 'frame', '--methods', 'iso,cs_hp', '--sets', 'm53'])
    E.main()
    assert calls == ['part_stamps', 'part_scr', 'part_ref', 'part_cs_test', ('cs_dev', 2), ('run', 'frame', ['iso', 'cs_hp'], ['m53'])]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_insitu.py'])
    with pytest.raises(SystemExit):
        E.main()


def test_run_cs3_passes_the_declared_options_to_fit_im_with_the_front_end_of_round_two(monkeypatch):
    import types
    calls = {}
    cs, psfmod, base = types.ModuleType('crowdsource'), types.ModuleType('crowdsource.psf'), types.ModuleType('crowdsource.crowdsource_base')

    def central_stamp(a, n):
        c, h = a.shape[0] // 2, n // 2
        return a[c - h: c + h + 1, c - h: c + h + 1]

    class SimplePSF:
        def __init__(self, stamp, normalize=-1):
            calls['psf'] = (np.array(stamp), normalize)

    def fit_im(im, psf, **kw):
        calls['im'], calls['kw'] = im, kw
        print('crowdsource prints its centroids')
        return {'stars': {'x': np.array([1.0, 2.0]), 'y': np.array([10.0, 20.0]), 'flux': np.array([5.0, 6.0])}}
    psfmod.central_stamp, psfmod.SimplePSF, base.fit_im = central_stamp, SimplePSF, fit_im
    cs.psf, cs.crowdsource_base = psfmod, base
    for name, mod in (('crowdsource', cs), ('crowdsource.psf', psfmod), ('crowdsource.crowdsource_base', base)):
        monkeypatch.setitem(sys.modules, name, mod)
    st = lumpy_stamp()
    dn = planted_frame(st, STARS, seed=1)
    spec = dict(E.CS_BASE, thr=4.0, miniter=12, maxiter=40, fewstars=5, psfderiv=False, blendthreshu=3, titer_thresh=5, psfvalsharpcutfac=0.9, psfsharpsat=1.0, weights='flat', sb=3.0)
    cat, secs = E.run_cs3(dn, st, GAIN, DARK, spec)
    kw = calls['kw']
    assert {k: v for k, v in kw.items() if k not in ('weights', 'dq')} == dict(maxstars=10000000, fewstars=5, miniter=12, maxiter=40, ntilex=1, ntiley=1, refit_psf=False, verbose=False, threshold=4.0,
                                                                          psfderiv=False, blendthreshu=3, titer_thresh=5, psfvalsharpcutfac=0.9, psfsharpsat=1.0)
    assert np.array_equal(kw['weights'], E.cs3_weights(dn, GAIN, DARK, 'flat')) and kw['weights'].dtype == np.float64 and kw['dq'].dtype == np.int64 and not kw['dq'].any() and kw['dq'].shape == dn.shape
    assert np.array_equal(calls['im'], E.HP.highpass_frame(dn, 3.0))                                       # the frame is high-passed with the variant's width
    stamp, norm = calls['psf']
    b = E.HP.cs_stamp(st)
    assert norm == 0 and np.allclose(stamp, E.HP.highpass_stamp(b / central_stamp(b, 19).sum(), 3.0))     # the high-passed stamp in the default's flux convention
    assert list(cat['x']) == [10.0, 20.0] and list(cat['y']) == [1.0, 2.0] and list(cat['flux']) == [5.0, 6.0] and secs >= 0       # crowdsource's x runs along the rows
    assert 'blist' not in kw
    E.run_cs3(dn, st, GAIN, DARK, dict(E.CS_BASE, blist=3))
    bl = calls['kw']['blist']
    assert len(bl[0]) == 3 and np.array_equal(bl[0], E.cs3_blist(dn, GAIN, DARK, 3)[0])


def test_the_committed_report_regenerates_from_the_committed_raw_byte_for_byte(monkeypatch):
    monkeypatch.setattr(E, 'TARGETS', REAL_TARGETS)
    if not (os.path.exists(E.RAW) and os.path.exists(E.MD)):
        pytest.skip('no committed raw / report yet')
    assert open(E.MD, 'rb').read() == E.report(E.load_raw(E.RAW)).encode()


# ----------------------------------------------------------------------------------------------------------------------------------
# Amendment 1: the control stamps
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_control_stamp_is_the_psfield_stamp_plus_the_in_situ_difference_in_a_random_order():
    rng = np.random.default_rng(0)
    psf = 0.01 + 0.0001 * rng.random((51, 51))
    psf /= psf.sum()
    insitu = psf + 1e-5 * rng.standard_normal((51, 51))
    d = (insitu - psf).ravel()
    s = E.scrambled_stamp(psf, insitu, [1, 2, 3])
    assert s.shape == (51, 51) and s.sum() == pytest.approx(1.0) and s.min() > 0
    back = np.sort((s * (psf.sum() + d.sum()) - psf).ravel())                                            # no pixel was clipped: the multiset of the added pixels is the difference's
    assert np.allclose(back, np.sort(d), atol=1e-15)
    assert not np.allclose(s * (psf.sum() + d.sum()), insitu)                                              # in a different order
    assert np.array_equal(s, E.scrambled_stamp(psf, insitu, [1, 2, 3])) and not np.array_equal(s, E.scrambled_stamp(psf, insitu, [1, 2, 4]))
    clipped = E.scrambled_stamp(np.zeros((51, 51)) + 1e-9, np.full((51, 51), -1.0) + np.arange(2601).reshape(51, 51) / 1300.0, 7)
    assert clipped.min() >= 0 and clipped.sum() == pytest.approx(1.0)                                    # negative pixels clipped


def test_the_control_seed_depends_on_the_set_and_the_origin():
    a, b, c, d = E.scr_seed('m2', (310, 630)), E.scr_seed('m53', (310, 630)), E.scr_seed('m2', (310, 730)), E.scr_seed('m2', (630, 310))
    assert a == E.scr_seed('m2', (310, 630)) and len({tuple(x) for x in (a, b, c, d)}) == 4 and a[0] == E.SCR_SEED and a[2:] == [310, 630]


def test_part_scr_stores_one_control_stamp_per_patch_with_a_stamp_and_none_for_a_missing_or_nan_stamp(monkeypatch, tmp_path):
    st, wrong, cache = planted_world(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_stamps(path, B=4)
    raw = E.load_raw(path)
    bad = dict(raw)
    bad['ST|patch|%s|stamp' % E.pkey('m2', (1000, 0))] = np.full((51, 51), np.nan)
    p2 = str(tmp_path / 'p2.npz')
    np.savez_compressed(p2, **bad)
    E.part_scr(p2)
    r2 = E.load_raw(p2)
    assert 'ST|scr|%s|stamp' % E.pkey('m2', (1000, 0)) not in r2
    E.part_scr(path)
    r = E.load_raw(path)
    k = 'ST|scr|%s|' % E.pkey('m2', (0, 0))
    ins, psf = raw['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))], wrong
    assert np.array_equal(r[k + 'stamp'], E.scrambled_stamp(psf, ins, E.scr_seed('m2', (0, 0)))) and r[k + 'vs_psf'].shape == (4,)
    assert r[k + 'rms_to_psf'][0] == pytest.approx(np.sqrt(np.mean((r[k + 'stamp'] - psf) ** 2)) / psf.max()) and r[k + 'rms_to_psf'][1] == pytest.approx(np.sqrt(np.mean((ins - psf) ** 2)) / psf.max())
    assert np.array_equal(r[k + 'stamp'], r['ST|scr|%s|stamp' % E.pkey('m2', (0, 0))]) and not np.array_equal(r[k + 'stamp'], r['ST|scr|%s|stamp' % E.pkey('m2', (1000, 0))])
    E.part_scr(path)                                                                                      # idempotent


def test_the_control_column_appears_in_the_tables_with_the_psfield_pairing(planted_raw):
    raw, _ = planted_raw
    text = E.report(raw)
    rows = table_rows(text, '### 2.1 Pooled by group')
    assert rows[0][-1].startswith('control: psField + scrambled in-situ difference')
    e24 = table_rows(text.split('### 2.4 ')[1], '**M53 (NGC 5024)**')
    assert e24[0][-1].startswith('control: psField + scrambled in-situ difference minus psField')
    row = next(r for r in e24 if r[0] == E.LABEL['cs_hp'])
    a, b = raw['RUN|scr|%s|cs_hp|k2|blk22' % E.pkey('m53', (0, 0))], raw['REF|%s|cs_hp|k2|blk22' % E.pkey('m53', (0, 0))]
    df = E.X.boot_diff(a, b, E.X.boot_idx_strat([25]))[2]
    assert row[-1] == '%s (1)' % E.diff_cell(df[0], df[1]) and 'loss' in row[-1]
    assert 'psField stamp plus a random permutation of the in-situ perturbation' in text and '* SENSITIVITY.' in text


def test_the_base_check_compares_the_development_totals_with_the_stored_round_two_selection(tmp_path, monkeypatch):
    path, hp2 = str(tmp_path / 'insitu_raw.npz'), str(tmp_path / 'hp2.npz')
    pa, pb = E.cs_dev_patches()
    tot = np.array([100, 40, 80, 40, 150, 40, 90])
    E.save_with_meta({'CSD|cs3_base|x%d_y%d|tot' % pa: tot, 'CSD|cs3_base|x%d_y%d|tot' % pb: tot}, path)
    np.savez_compressed(hp2, **{'DEV|x%d_y%d|%s|tot' % (pa[0], pa[1], E.CS_BASE_TAG): tot, 'DEV|x%d_y%d|%s|tot' % (pb[0], pb[1], E.CS_BASE_TAG): tot + 1})
    E.part_cs_check(path, hp2)
    raw = E.load_raw(path)
    assert raw['CSD|check|x%d_y%d' % pa][0] == 1 and raw['CSD|check|x%d_y%d' % pb][0] == 0
    E.part_cs_check(str(tmp_path / 'other.npz'), str(tmp_path / 'missing.npz'))                    # no stored file: nothing recorded, no error
    assert 'CSD|check|x%d_y%d' % pa not in E.load_raw(str(tmp_path / 'other.npz'))
    text = E.report(dict(raw, **{'CSD|cs3_base|x%d_y%d|secs' % pa: np.array([1.0])}))
    assert 'Check: the base reproduces the stored development totals of `cs2_hp2_t3_pm` (`baselines_hp2_raw.npz`) on 1 of 2 patches.' in text


# ----------------------------------------------------------------------------------------------------------------------------------
# model-based tests: independent implementations of the documented procedure
# ----------------------------------------------------------------------------------------------------------------------------------
def ref_sky_clip(v, k=3.0, passes=5):
    v = np.asarray(v, float)
    keep = np.ones(len(v), bool)
    for _ in range(passes):
        mu, sd = v[keep].mean(), v[keep].std()
        new = np.abs(v - mu) <= k * sd
        if new.sum() == keep.sum() or new.sum() < 3:
            break
        keep = new
    return v[keep].mean()


def test_the_clipped_mean_agrees_with_an_independent_implementation_on_skewed_and_tiny_samples():
    rng = np.random.default_rng(11)
    samples = [rng.lognormal(0, 1.0, 300) for _ in range(5)] + [rng.exponential(5.0, 40) for _ in range(5)] + [rng.normal(0, 1, n) for n in (3, 4, 5, 6)] + [np.r_[rng.normal(0, 1, 5), 40.0]]
    for v in samples:
        assert E.sky_clip(v) == pytest.approx(ref_sky_clip(v), rel=1e-12, abs=1e-12)
    assert E.sky_clip(samples[0], k=2.0, passes=2) == pytest.approx(ref_sky_clip(samples[0], k=2.0, passes=2), rel=1e-12)


def test_the_aligned_cut_out_is_the_cubic_spline_shift_to_the_flux_centroid_with_edge_values_repeated():
    rng = np.random.default_rng(2)
    st = lumpy_stamp()
    sub = st[25 - 12:25 + 13, 25 - 12:25 + 13] * 1e5 + 5.0 + rng.normal(0, 3.0, (25, 25)) + np.linspace(0, 20, 25)[None, :]
    sh, ox, oy = E.align_unnorm(sub)
    cx, cy = E.D.centroid(sub)
    want = E.nd_shift(sub, (12 - cy, 12 - cx), order=3, mode='nearest')[3:-3, 3:-3]
    assert np.array_equal(sh, want) and sh.shape == (19, 19) and ox == pytest.approx(cx - 12) and oy == pytest.approx(cy - 12)


def test_star_cores_subtract_the_clipped_mean_of_a_skewed_ring_and_cut_the_right_rows_and_columns():
    rng = np.random.default_rng(6)
    st = lumpy_stamp()
    dn = rng.exponential(20.0, (100, 100)) + 100.0
    dn += E.D.render_sources(st, np.array([30.0, 70.0]), np.array([62.0, 35.0]), np.array([8e4, 6e4]), (100, 100), 20)
    x, y = np.array([30.0, 70.0]), np.array([62.0, 35.0])
    c, off = E.star_cores(dn, x, y, np.array([1, 0]))                                       # the order of idx is the order of the output
    for k, i in enumerate((1, 0)):
        ix, iy = int(round(x[i])), int(round(y[i]))
        sub = dn[iy - 12:iy + 13, ix - 12:ix + 13].copy()
        ring = np.hypot(*np.mgrid[-12:13, -12:13]) > 8.5
        sub = sub - ref_sky_clip(sub[ring])
        want, ox, oy = E.align_unnorm(sub)
        assert np.allclose(c[k], want, atol=1e-9) and np.allclose(off[k], (ox, oy))


def ref_robust(c, n_it=4, kk=3.0, floor=0.10):
    disk = E.DISK
    n = len(c)
    ssum = c[:, disk].sum(1)
    ok = ssum > 0
    s = np.median(c[ok] / ssum[ok][:, None, None], 0)
    s = s / s[disk].sum()
    use = ok.copy()
    for _ in range(n_it):
        f = np.array([(ci * s).sum() / (s * s).sum() for ci in c])
        e = np.array([np.sqrt(np.mean((ci - fi * s)[disk] ** 2)) / abs(fi * s.max()) for ci, fi in zip(c, f)])
        thr = max(kk * np.median(e[use]), floor)
        use = (f > 0) & (e <= thr)
        s = sum(f[i] * c[i] for i in np.flatnonzero(use)) / sum(f[i] ** 2 for i in np.flatnonzero(use))
        s = s / s[disk].sum()
    return s, f, e, use


def test_the_stack_agrees_with_an_independent_implementation_of_the_documented_iteration():
    rng = np.random.default_rng(8)
    s0 = shape_core(1.5, 1.1)
    for trial in range(4):
        amps = rng.uniform(1.0, 20.0, 7)
        c = amps[:, None, None] * s0 * (1 + 0.05 * rng.standard_normal((7, 1, 1))) + (0.002 + 0.01 * rng.random((7, 1, 1))) * rng.standard_normal((7, 19, 19))
        if trial >= 2:
            c[3, 10:13, 4:8] += 0.8 * amps[3] * s0.max()                                    # a contaminated star
            c[5] = -c[5]                                                                    # a negative one
        s, f, e, use = E.robust_stack(c)
        rs, rf, re_, ru = ref_robust(c)
        assert np.array_equal(use, ru) and np.allclose(s, rs, rtol=1e-10, atol=1e-14) and np.allclose(f, rf, rtol=1e-10) and np.allclose(e, re_, rtol=1e-10)
    z = np.zeros((4, 19, 19))
    z[:3] = np.array([3.0, 2.0, 4.0])[:, None, None] * s0
    s, f, e, use = E.robust_stack(z)                                                         # a star of zero flux is not a candidate for the initial shape and is not kept
    assert list(use) == [True, True, True, False] and np.isfinite(s).all() and f[3] == 0


def test_the_forward_model_agrees_with_a_direct_computation_for_two_offsets_and_weights():
    st = lumpy_stamp(1.2, 1.0, 0.1)
    offs = np.array([[0.31, -0.22], [-0.4, 0.17]])
    w = np.array([0.3, 0.9])
    acc = np.zeros((19, 19))
    ring = np.hypot(*np.mgrid[-12:13, -12:13]) > 8.5
    for (ox, oy), wi in zip(offs, w):
        sub = E.nd_shift(st, (oy, ox), order=3, mode='constant', cval=0.0)[13:38, 13:38]
        sub = sub - ref_sky_clip(sub[ring])
        core, _, _ = E.align_unnorm(sub)
        acc += wi * core / core[E.DISK].sum()
    assert np.allclose(E.forward_core(st, offs, w), acc / acc[E.DISK].sum(), atol=1e-12, rtol=0)
    wide = lumpy_stamp(3.0, 3.0, 0.0)                                                         # a wide stamp has wings in the ring: skipping the ring subtraction would change the core
    sub = E.nd_shift(wide, (0.2, 0.3), order=3, mode='constant', cval=0.0)[13:38, 13:38]
    assert abs(ref_sky_clip(sub[ring])) > 1e-5


def test_estimate_defaults_normalisation_and_residual_follow_the_documented_definitions(planted):
    st, dn, wrong = planted
    stars = E.collect([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    est = E.estimate(stars['c'], stars['off'], wrong)
    e3 = E.estimate(stars['c'], stars['off'], wrong, n_fix=3)
    e0 = E.estimate(stars['c'], stars['off'], wrong, n_fix=0)
    e1 = E.estimate(stars['c'], stars['off'], wrong, n_fix=1)
    assert np.array_equal(est['stamp'], e3['stamp']) and not np.allclose(est['stamp'], e0['stamp'], atol=1e-8) and not np.allclose(est['stamp'], e1['stamp'], atol=1e-8)
    assert est['core'][E.DISK].sum() == pytest.approx(1.0, abs=1e-12) and est['ok']
    w = est['f'][est['use']] ** 2 / (est['f'][est['use']] ** 2).sum()
    t = E.forward_core(est['stamp'], stars['off'][est['use']], w)
    assert est['resid'] == pytest.approx(np.sqrt(np.mean((est['core_obs'] - t) ** 2)) / est['core_obs'].max(), rel=1e-12)
    # exactly N_MIN kept stars still give a stamp
    three = E.estimate(stars['c'][:3], stars['off'][:3], wrong)
    assert three['ok'] and three['use'].sum() == 3 and np.isfinite(three['stamp']).all()
    # the weights of the correction are the squared amplitudes: a star with a tiny amplitude barely moves the correction
    c2, o2 = stars['c'].copy(), stars['off'].copy()
    c2[0] = c2[0] * 1e-3
    a, b = E.estimate(stars['c'], o2, wrong), E.estimate(c2, o2, wrong)
    assert np.isfinite(b['stamp']).all()


def test_more_than_twenty_isolated_stars_are_all_kept_by_default_up_to_forty():
    st = lumpy_stamp(1.4, 1.2, 0.05)
    pos = [(14.0 + 18 * i + 0.3 * (j % 2), 14.0 + 18 * j + 0.2 * (i % 2), 6e4) for i in range(5) for j in range(5)]
    dn = planted_frame(st, pos, seed=3)
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    xs, ys, pks, idx = E.select_stars(x, y, pk, noise, dn)
    assert len(idx) == 25
    assert len(E.select_stars(x, y, pk, noise, dn, n_max=10)[3]) == 10


def test_detection_is_the_high_passed_smoothed_frame_at_its_local_maxima(planted):
    st, dn, _ = planted
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    sm = E.gaussian_filter(E.HP.highpass_frame(dn, 3.0), 1.0)
    assert np.array_equal(pk, sm[y.astype(int), x.astype(int)])
    mx = E.maximum_filter(sm, size=5)
    assert np.all(sm[y.astype(int), x.astype(int)] == mx[y.astype(int), x.astype(int)]) and len(x) == int((sm == mx).sum())
    # two peaks 3 px apart: with a 5 x 5 footprint only the higher is a maximum
    two = planted_frame(lumpy_stamp(1.3, 1.3, 0.0), [(50.0, 50.0, 1e5), (53.0, 50.0, 6e4)], seed=1, noise=False)
    xx, yy, pp, nz = E.detect_maxima(two, GAIN, DARK)
    near = [(a, b) for a, b in zip(xx, yy) if abs(a - 51) < 4 and abs(b - 50) < 4]
    assert len(near) == 1


def test_selection_boundary_cases_of_the_neighbour_rules_the_saturation_box_and_the_array_edge():
    # a neighbour of exactly five sigma counts (>=): at 0.1 of a 50 sigma star, 3 px away, it excludes the star
    assert sel([50, 53], [50, 50], [50.0, 5.0])[0] == []
    assert sel([50, 53], [50, 50], [50.0, 4.99])[0] == [(50.0, 50.0)]
    # a faint neighbour straight above (dx = 0, 8 px) is not inside the 5 px inner radius
    got = sel([50, 50], [50, 58], [1000.0, 300.0])[0]
    assert (50.0, 50.0) in got
    # the saturation box is read at the row and column of the star: a saturated pixel 3 px to the right of a star at (x, y) = (30, 60)
    dn = np.full((100, 100), 100.0)
    dn[60, 33] = E.SAT_DN + 1
    assert sel([30], [60], [100.0], dn=dn)[0] == []
    dn = np.full((100, 100), 100.0)
    dn[33, 60] = E.SAT_DN + 1                                                                    # the transposed pixel is not in the box
    assert sel([30], [60], [100.0], dn=dn)[0] == [(30.0, 60.0)]
    # a star near the array edge (edge cut off): the box is clipped, not wrapped
    dn = np.full((100, 100), 100.0)
    dn[-1, 1] = E.SAT_DN + 1                                                                     # would be wrapped into the box by a negative start
    assert sel([1], [1], [100.0], dn=dn, edge=0)[0] == [(1.0, 1.0)]


# ----------------------------------------------------------------------------------------------------------------------------------
# second round of tests (written after the first mutation run left survivors)
# ----------------------------------------------------------------------------------------------------------------------------------
def test_collect_reports_the_peak_in_noise_units_and_three_stars_are_enough_for_a_stamp(planted):
    st, dn, wrong = planted
    s = E.collect([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    x, y, pk, noise = E.detect_maxima(dn, GAIN, DARK)
    for xi, yi, sn in zip(s['x'], s['y'], s['snr']):
        j = int(np.flatnonzero((x == xi) & (y == yi))[0])
        assert sn == pytest.approx(pk[j] / noise, rel=1e-12)
    three = E.insitu_stamp([dict(dn=planted_frame(st, STARS[:3], seed=5), gain=GAIN, dark=DARK)], wrong)
    assert len(three['stars']['c']) == 3 and three['ok'] and np.isfinite(three['stamp']).all()


def test_the_comparison_statistics_keep_their_entries_apart_for_stamps_of_different_widths():
    a, b = shape_core(1.2), shape_core(1.8)
    v = E.compare_cores(a, b)
    assert v[1] == pytest.approx(E.D.fwhm(a)) and v[2] == pytest.approx(E.D.fwhm(b)) and v[1] < v[2] - 0.5
    st_a, st_b = lumpy_stamp(1.2, 1.2, 0.0), lumpy_stamp(1.8, 1.8, 0.0)
    ref = E.stamp_core(st_b)
    sv = E.stat_vector(st_a, ref)
    full = E.compare_cores(E.stamp_core(st_a), ref)
    assert sv[0] == pytest.approx(full[0]) and sv[1] == pytest.approx(full[1]) and sv[2] == pytest.approx(full[3]) and sv[1] != pytest.approx(full[2])
    nanst = st_a.copy()
    nanst[10, 10] = np.nan
    assert np.isnan(E.stat_vector(nanst, ref)).all()


def test_truth_stars_exclude_saturation_at_the_limit_and_the_unmatched_are_not_isolated():
    dn = np.full((100, 100), 100.0)
    truth = dict(px=np.array([50.0]), py=np.array([50.0]), m606=np.array([17.0]))
    dn[50, 50] = E.SAT_DN
    assert len(E.truth_stars(truth, 0, 0, dn)[3]) == 1
    dn[50, 50] = E.SAT_DN + 1
    assert len(E.truth_stars(truth, 0, 0, dn)[3]) == 0
    tx, ty, tm = np.array([10.0, 30.0, 20.0]), np.array([10.0, 30.0, 20.0]), np.array([18.0, 18.0, 18.0])
    assert list(E.truth_isolated_mask(np.array([-1, 2]), tx, ty, tm)) == [False, True]


def test_the_star_bootstrap_equals_an_explicit_resampling_of_the_stars_with_the_stated_seed_and_strata(planted):
    st, dn, wrong = planted
    est = E.insitu_stamp([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    s, ref = est['stars'], E.stamp_core(st)
    sizes = [4, 3]
    ci, nf = E.star_bootstrap(s, [ref], [wrong], sizes, B=7, seed=3)
    idx = E.X.boot_idx_strat(sizes, 7, 3)
    draws = np.array([E.stat_vector(E.estimate(s['c'][i], s['off'][i], wrong)['stamp'], ref) for i in idx])
    assert np.allclose(ci[0], np.percentile(draws, [2.5, 97.5], axis=0).T) and nf == 0
    ci_b, _ = E.star_bootstrap(s, [ref], [wrong], sizes, B=7, seed=4)
    assert not np.allclose(ci, ci_b)
    ci_c, _ = E.star_bootstrap(s, [ref], [wrong], [7], B=7, seed=3)
    assert not np.allclose(ci, ci_c)
    bad = dict(s, c=s['c'][:3], off=s['off'][:3])
    bad['c'] = -np.abs(bad['c'])                                                                   # all amplitudes negative: no resample has an estimate
    ci_n, nf_n = E.star_bootstrap(bad, [ref], [wrong], [3], B=5, seed=0)
    assert nf_n == 5 and np.isnan(ci_n).all()


def test_the_held_out_fit_is_exact_for_the_true_stamp_has_a_standardized_rms_near_one_and_respects_the_frame_edges(planted):
    st, dn, wrong = planted
    noise_free = planted_frame(st, STARS[:3], noise=False, sky=0.0)
    r = E.hp_fit_rms(noise_free, GAIN, DARK, st, [(s[0], s[1]) for s in STARS[:3]])
    assert np.all(r[:, 1] < 3e-3)                                                                  # misfit of a noise-free star fitted by its own stamp, shifts included
    rw = E.hp_fit_rms(noise_free, GAIN, DARK, wrong, [(s[0], s[1]) for s in STARS[:3]])
    assert np.all(rw[:, 1] > 10 * r[:, 1])
    rn = E.hp_fit_rms(dn, GAIN, DARK, st, [(s[0], s[1]) for s in STARS[:3]])
    assert np.all((rn[:, 0] > 0.7) & (rn[:, 0] < 1.5))                                              # noise only, weights right: standardized rms near 1
    img = planted_frame(st, [(94.0, 50.0, 1e5), (50.0, 94.0, 1e5), (95.0, 50.0, 1e5), (50.0, 95.0, 1e5)], noise=False, sky=0.0)
    r = E.hp_fit_rms(img, GAIN, DARK, st, [(94.0, 50.0), (50.0, 94.0), (95.0, 50.0), (50.0, 95.0)])
    assert np.isfinite(r[:2]).all() and np.isnan(r[2:]).all()                                       # the 11 x 11 box must lie inside the frame (index 94 is the last that does)


def test_the_two_run_functions_get_the_stamp_the_frozen_configuration_and_the_kernel_fit(monkeypatch):
    seen = []

    def fake_path(patch, hp, rex, *a, **k):
        seen.append(('iso', patch['tag'], hp, rex))
        E.IM.sliding_path(*tiny_args())
        return {0.0: (np.zeros(0), np.zeros(0), np.zeros(0)), 2.0: (np.array([1.0]), np.array([2.0]), np.array([3.0]))}, 1.5, 7, dict(bic_index=4, n_blocks=9)

    def fake_psf(patch, hp, rex, fit, *a, **k):
        seen.append(('adq', patch['tag'], hp, rex, fit))
        E.IM.sliding_path(*tiny_args())
        return {0.0: (np.zeros(0), np.zeros(0), np.zeros(0)), 2.0: (np.array([1.0]), np.array([2.0]), np.array([3.0]))}, 2.5, 8, dict(bic_index=5, n_blocks=11)

    def tiny_args():
        P = E.IM.PixelPSF(24, 24, 1.2)
        L = E.IM.Lattice(P, h=0.5, margin=2.0)
        return P, P.synth(np.array([[12.0, 12.0]]), np.array([100.0])), np.ones((24, 24)), L

    monkeypatch.setattr(E.X, 'run_path', fake_path)
    monkeypatch.setattr(E.PS, 'run_path_psf', fake_psf)
    monkeypatch.setattr(E.X, 'score_dets', lambda d, t: E.A.score(np.c_[d[0], d[1]], d[2], t))
    t = dict(x=np.array([2.0]), y=np.array([2.0]), mag=np.array([20.0]), inner=np.array([True]))
    orig_call = E.IM.sliding_path
    arr = E.run_iss('iso', dict(tag='P'), t, (3.0, 1.2, 2.0))
    arr2 = E.run_iss('adq', dict(tag='P'), t, (3.0, 1.2, 2.0), 'FIT')
    assert seen == [('iso', 'P', 3.0, 1.2), ('adq', 'P', 3.0, 1.2, 'FIT')] and E.IM.sliding_path is orig_call
    assert arr['secs'].tolist() == [1.5, 7] and arr2['secs'].tolist() == [2.5, 8] and arr['bic'].tolist() == [4, 9, 0] and arr2['bic'].tolist() == [5, 11, 0]
    assert arr['det'].tolist() == [[1.0, 2.0, 3.0]]                                                  # the catalogue at kappa 2, not the BIC readout
    assert 'k2|tot' in arr and 'k0|tot' in arr


def test_pack_iss_flags_exactly_the_cap_and_stores_the_two_readouts(monkeypatch):
    monkeypatch.setattr(E.X, 'score_dets', lambda d, t: E.A.score(np.c_[d[0], d[1]], d[2], t))
    t = dict(x=np.array([2.0]), y=np.array([2.0]), mag=np.array([20.0]), inner=np.array([True]))
    dets = {2.0: (np.array([2.0]), np.array([2.0]), np.array([1e4])), 0.0: (np.array([3.0, 4.0]), np.array([3.0, 4.0]), np.array([1e3, 1e3]))}
    tr = dict(counts=np.arange(8), barred=np.zeros((0, 3)), taus=np.zeros(0))
    for n, flag in ((E.X.MAX_ENTRIES - 1, 0), (E.X.MAX_ENTRIES, 1), (E.X.MAX_ENTRIES + 1, 1)):
        a = E.pack_iss(dets, 3.0, n, dict(bic_index=2, n_blocks=40), t, 2.0, tr)
        assert a['bic'].tolist() == [2, 40, flag] and a['secs'].tolist() == [3.0, n]
    assert a['det'].shape == (1, 3) and a['k2|tot'][2] == 1 and a['k0|tot'][2] == 2 and a['trace|counts'].tolist() == list(range(8))


def test_run_baseline_passes_the_crowdsource_spec_to_crowdsource_and_the_photutils_spec_to_photutils(monkeypatch):
    got = {}
    cat = dict(x=np.array([14.0]), y=np.array([14.0]), flux=np.array([6.0e4]))
    monkeypatch.setattr(E.H2, 'run_crowdsource_hp2', lambda dn, psf, gain, dark, spec: (got.setdefault('cs', spec), (cat, 1.0))[1])
    monkeypatch.setattr(E.H2, 'run_photutils_hp2', lambda dn, psf, gain, dark, spec: (got.setdefault('phu', spec), (cat, 1.0))[1])
    st = lumpy_stamp()
    patch, truth = small_patch(st)
    t = dict(x=truth['px'], y=truth['py'], mag=truth['m606'], inner=np.ones(3, bool))
    E.run_baseline('cs_hp', patch, 1.0, t, ('CS', 'PHU'))
    E.run_baseline('phu_hp', patch, 1.0, t, ('CS', 'PHU'))
    assert got == {'cs': 'CS', 'phu': 'PHU'}


def test_round2_specs_are_the_registered_variants_of_the_verified_selection(monkeypatch):
    monkeypatch.setattr(E.CF, 'load_raw', lambda p: {})
    monkeypatch.setattr(E.HS, 'verify_selection2', lambda raw: dict(crowdsource='cs2_hp3_t3_pd', photutils='phu2_hp5_new_i3_b2'))
    cs, phu = E.round2_specs()
    assert cs['tag'] == 'cs2_hp3_t3_pd' and phu['tag'] == 'phu2_hp5_new_i3_b2' and cs['sb'] == 3.0 and phu['bounds'] == 2.0


def test_the_hash_record_covers_the_freeze_file_when_it_exists(monkeypatch, tmp_path):
    root = tmp_path / 'repo'
    (root / 'results' / 'astro').mkdir(parents=True)
    (root / 'results' / 'astro' / 'freeze.json').write_text('{"a": 1}')
    monkeypatch.setattr(E, 'REPO', str(root))
    m = E.meta_arrays()
    assert 'results/astro/freeze.json' in list(m['meta|files']) and len(m['meta|files']) == 1
    h = dict(zip(m['meta|files'], m['meta|sha256']))
    assert h['results/astro/freeze.json'] == E.sha256_file(str(root / 'results' / 'astro' / 'freeze.json'))


def test_part_run_second_round_sky_scope_sets_hashes_and_the_best_variant_spec(monkeypatch, tmp_path):
    st, wrong, cache = planted_world(monkeypatch)
    for key in list(cache):
        patch, sky, truth = cache[key]
        cache[key] = (patch, SKY + np.arange(10000.0).reshape(100, 100) * 0.001, truth)            # a sky whose mean differs from its maximum
    path = str(tmp_path / 'insitu_raw.npz')
    monkeypatch.setattr(E, 'frozen_config', lambda: (3.0, 1.2, 2.0))
    monkeypatch.setattr(E, 'round2_specs', lambda: ('cs-spec', 'phu-spec'))
    E.part_stamps(path, B=4)
    calls = []

    def fake_iss(kind, patch, t, cfg, fit=None):
        calls.append((kind, float(patch['psf'][25, 25])))
        blk, tot = fake_scores(len(calls))
        return {'k2|tot': tot, 'k2|blk22': blk, 'det': np.zeros((2, 3)), 'secs': np.array([5.0, 100.0])}

    def fake_base(meth, patch, sky_mean, t, specs):
        calls.append((meth, sky_mean, float(patch['psf'][25, 25])))
        blk, tot = fake_scores(len(calls))
        return {'k2|tot': tot, 'k2|blk22': blk, 'det': np.zeros((4, 3)), 'secs': np.array([2.0])}
    monkeypatch.setattr(E, 'run_iss', fake_iss)
    monkeypatch.setattr(E, 'run_baseline', fake_base)
    monkeypatch.setattr(E.PG, 'fit_adequate', lambda stamp: dict(kind='general', rms_isotropic=5e-3, rms=1e-3))
    E.part_run('frame', ['iso', 'cs_def'], ['m2'], path)
    raw = E.load_raw(path)
    fr0 = raw['ST|frame|%s|stamp' % E.pkey('m2', (0, 0))]
    assert calls[0] == ('iso', float(fr0[25, 25])) and calls[1][0] == 'cs_def' and calls[1][1] == pytest.approx(float(cache[('m2', (0, 0))][1].mean()), rel=1e-12) and calls[1][1] != cache[('m2', (0, 0))][1].max()
    assert calls[1][2] == float(fr0[25, 25])                                                         # the frame-scope stamp, not the patch-scope one
    assert not np.isclose(raw['ST|frame|%s|stamp' % E.pkey('m2', (0, 0))][25, 25], raw['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))][25, 25], atol=1e-12)
    assert str(raw['RUN|frame|%s|iso|kernel' % E.pkey('m2', (0, 0))][0]) == 'isotropic'              # the isotropic run is the isotropic kernel whatever the rule says
    assert list(raw['RUN|frame|%s|iso|rms' % E.pkey('m2', (0, 0))]) == [5e-3, 1e-3]
    # a set not in the list is not run; a changed hash refuses
    n = len(calls)
    E.part_run('frame', ['iso'], ['m53'], path)
    assert len(calls) == n
    bad = dict(raw)
    bad['meta|sha256'] = np.array(['0' * 64] + [str(h) for h in raw['meta|sha256'][1:]])
    p2 = str(tmp_path / 'bad.npz')
    np.savez_compressed(p2, **bad)
    with pytest.raises(RuntimeError, match='differ'):
        E.part_run('frame', ['iso'], ['m2'], p2)
    # a stamp with a NaN anywhere (not only in its first pixel) is skipped
    nan = dict(raw)
    s = nan['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))].copy()
    s[30, 31] = np.nan
    nan['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))] = s
    p3 = str(tmp_path / 'nan.npz')
    np.savez_compressed(p3, **{k: v for k, v in nan.items() if not k.startswith('RUN|')})
    calls.clear()
    E.part_run('patch', ['iso'], ['m2'], p3)
    assert len(calls) == 1 and calls[0][1] == float(raw['ST|patch|%s|stamp' % E.pkey('m2', (1000, 0))][25, 25])


def test_part_cs_test_runs_the_combination_with_its_declared_spec_on_the_test_patches(monkeypatch, tmp_path):
    path, ran = run_part_cs_dev_with(monkeypatch, tmp_path, {'cs3_base': 0.440, E.cs3_tag('blendthreshu', 4): 0.4504, E.cs3_tag('maxiter', 30): 0.4431, 'cs3_combo': 0.470}, 'combo.npz')
    E.part_cs_dev(1, path)
    E.part_cs_dev(2, path)
    raw = E.load_raw(path)
    assert str(raw['CSD|SEL|best'][0]) == 'cs3_combo' and raw['CSD|SEL|finding'][0] == 1
    seen = []
    monkeypatch.setattr(E, 'part_cs_variant', lambda tag, spec, patches, truth, prefix, path=None: seen.append((tag, spec, [tuple(p) for p in patches], prefix)))
    E.part_cs_test(path)
    assert [s[0] for s in seen] == ['cs3_base', 'cs3_combo'] and all(s[3] == 'CSTEST' and s[2] == [tuple(p) for p in E.X.TESTS] for s in seen)
    assert seen[0][1] == E.CS_BASE and seen[1][1] == dict(E.CS_BASE, blendthreshu=4, maxiter=30)


def test_run_cs3_normalises_the_stamp_to_the_default_flux_convention_for_a_wide_stamp(monkeypatch):
    import types
    calls = {}
    cs, psfmod, base = types.ModuleType('crowdsource'), types.ModuleType('crowdsource.psf'), types.ModuleType('crowdsource.crowdsource_base')

    def central_stamp(a, n):
        c, h = a.shape[0] // 2, n // 2
        return a[c - h: c + h + 1, c - h: c + h + 1]

    class SimplePSF:
        def __init__(self, stamp, normalize=-1):
            calls['psf'] = np.array(stamp)
    base.fit_im = lambda im, psf, **kw: {'stars': {'x': np.zeros(1), 'y': np.zeros(1), 'flux': np.ones(1)}}
    psfmod.central_stamp, psfmod.SimplePSF = central_stamp, SimplePSF
    cs.psf, cs.crowdsource_base = psfmod, base
    for name, mod in (('crowdsource', cs), ('crowdsource.psf', psfmod), ('crowdsource.crowdsource_base', base)):
        monkeypatch.setitem(sys.modules, name, mod)
    wide = lumpy_stamp(7.0, 7.0, 0.0)
    b = E.HP.cs_stamp(wide)
    s19 = central_stamp(b, 19).sum()
    assert 0.5 < s19 < 0.95
    E.run_cs3(np.full((60, 60), 100.0), wide, GAIN, DARK, dict(E.CS_BASE))
    assert np.allclose(calls['psf'], E.HP.highpass_stamp(b / s19, 2.0), rtol=1e-9, atol=1e-15) and not np.allclose(calls['psf'], E.HP.highpass_stamp(b, 2.0), rtol=1e-3, atol=1e-12)


def world_with_pairs(monkeypatch):
    """Two planted patches of one frame with DIFFERENT psField stamps and, in the first, an extra pair of truth stars 8 px apart (comparable fluxes: the in-situ rules exclude both, the truth rule keeps both)."""
    st = lumpy_stamp()
    wrongs = [lumpy_stamp(2.0, 1.8, 0.0), lumpy_stamp(2.4, 1.6, 0.0)]
    cache = {}
    sets = (('m2', (0, 0)), ('m2', (1000, 0)))
    for k, (name, pxy) in enumerate(sets):
        stars = [(x + 0.1 * k, y, f) for x, y, f in STARS]
        if k == 0:
            stars += [(76.0, 84.0, 1.0e5), (84.0, 84.0, 7.0e4)]
        dn = planted_frame(st, stars, seed=20 + k)
        truth = dict(px=np.array([s[0] for s in stars]) + pxy[0], py=np.array([s[1] for s in stars]) + pxy[1], m606=A.counts_to_mag(np.array([s[2] for s in stars])))
        sky = SKY + np.arange(10000.0).reshape(100, 100) * 0.0005
        cache[(name, pxy)] = (dict(dn=dn, psf=wrongs[k], gain=GAIN, dark=DARK), sky, truth)
    monkeypatch.setattr(E, 'TARGETS', list(sets))
    monkeypatch.setattr(E, '_CACHE', {})
    monkeypatch.setattr(E, 'inputs', lambda name, pxy: cache[(name, tuple(pxy))])
    return st, wrongs, cache, sets


def test_part_stamps_wiring_second_round(monkeypatch, tmp_path):
    st, wrongs, cache, sets = world_with_pairs(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_stamps(path, B=5)
    raw = E.load_raw(path)
    for k, (name, pxy) in enumerate(sets):
        patch, sky, truth = cache[(name, pxy)]
        key = E.pkey(name, pxy)
        ref = E.stamp_core(patch['psf'])
        fr = dict(dn=patch['dn'], gain=GAIN, dark=DARK)
        ep = E.insitu_stamp([fr], patch['psf'])
        assert np.array_equal(raw['ST|patch|%s|stamp' % key], ep['stamp'])
        # patch scope: bootstrap with seed 0 over the selected stars, the stored comparison with the psField core
        ci, nf = E.star_bootstrap(ep['stars'], [ref], [patch['psf']], ep['stars']['n_by_frame'], 5, 0)
        assert np.array_equal(raw['ST|patch|%s|vs_psf_ci' % key], ci) and np.array_equal(raw['ST|patch|%s|vs_psf' % key], E.compare_cores(E.stamp_core(ep['stamp']), ref))
        # frame scope: the pooled stars through THIS patch's wings, bootstrap stratified by patch
        frames = [dict(dn=cache[s][0]['dn'], gain=GAIN, dark=DARK) for s in sets]
        pooled = E.collect(frames, patch['psf'])
        ef = E.estimate(pooled['c'], pooled['off'], patch['psf'])
        assert np.array_equal(raw['ST|frame|%s|stamp' % key], ef['stamp'])
        cif, _ = E.star_bootstrap(pooled, [ref], [patch['psf']], pooled['n_by_frame'], 5, 0)
        assert np.array_equal(raw['ST|frame|%s|vs_psf_ci' % key], cif)
        assert np.array_equal(raw['ST|frame|%s|vs_psf' % key], E.compare_cores(E.stamp_core(ef['stamp']), ref))
        assert np.array_equal(raw['ST|psf|%s|core' % key], ref) and np.array_equal(raw['ST|psf|%s|shape' % key], E.shape_of(patch['psf']))
        # the world of the patch is rendered with the seed of the patch
        w = E.D.render_world(truth, pxy[0], pxy[1], patch['psf'], GAIN, DARK, sky, E.D.world_seed(name, pxy))
        ew = E.insitu_stamp([dict(dn=w['dn'], gain=GAIN, dark=DARK)], patch['psf'])
        assert raw['ST|patch|%s|world|n_sel' % key][0] == len(ew['stars']['c']) and np.array_equal(raw['ST|patch|%s|world|vs_planted' % key], E.compare_cores(E.stamp_core(ew['stamp']), ref))
        sx, sy = ew['stars']['x'], ew['stars']['y']
        assert np.array_equal(raw['ST|patch|%s|world|plain' % key], E.plain_stack_stat(w['dn'], sx, sy, np.arange(len(sx)), ref))
        # the plain stack of the in-situ stars and the truth match
        px, py = ep['stars']['x'], ep['stars']['y']
        assert np.array_equal(raw['ST|patch|%s|plain' % key], E.plain_stack_stat(patch['dn'], px, py, np.arange(len(px)), ref))
        tx, ty, tm = E.D.sources_near(truth, pxy[0], pxy[1], E.D.MARGIN)
        mag, ind = E.nearest_truth(px, py, tx, ty, tm)
        iso = E.truth_isolated_mask(ind, tx, ty, tm)
        assert raw['ST|patch|%s|truthmatch' % key].tolist()[:2] == [int(np.isfinite(mag).sum()), int(iso.sum())]
        # vs the truth-star stamp (same estimator) and the truth stars' own comparison
        txs, tys, tms, tidx = E.truth_stars(truth, pxy[0], pxy[1], patch['dn'])
        c_t, off_t = E.star_cores(patch['dn'], txs, tys, tidx)
        et = E.estimate(c_t, off_t, patch['psf'])
        assert np.array_equal(raw['ST|patch|%s|truth|vs_psf' % key], E.compare_cores(E.stamp_core(et['stamp']), ref))
        assert np.array_equal(raw['ST|patch|%s|vs_truth' % key], E.compare_cores(E.stamp_core(ep['stamp']), E.stamp_core(et['stamp'])))
        assert np.array_equal(raw['ST|patch|%s|truth|plain' % key], E.plain_stack_stat(patch['dn'], txs, tys, tidx, ref))
        assert raw['ST|patch|%s|truth|n_used' % key][0] == int(et['use'].sum()) and raw['ST|patch|%s|n_used' % key][0] == int(ep['use'].sum())
    # held-out stars: the pair 8 px apart in the first patch, in the columns [std psField, frac psField, std in-situ, frac in-situ]
    h0 = raw['ST|patch|%s|held' % E.pkey('m2', (0, 0))]
    assert h0.shape == (2, 4) and np.all(h0[:, 1] > h0[:, 3]) and np.all(h0[:, 0] > 0)
    patch, sky, truth = cache[('m2', (0, 0))]
    ep = E.insitu_stamp([dict(dn=patch['dn'], gain=GAIN, dark=DARK)], patch['psf'])
    want_psf = E.hp_fit_rms(patch['dn'], GAIN, DARK, patch['psf'], [(76.0, 84.0), (84.0, 84.0)])
    want_ins = E.hp_fit_rms(patch['dn'], GAIN, DARK, ep['stamp'], [(76.0, 84.0), (84.0, 84.0)])
    assert np.allclose(h0, np.hstack([want_psf, want_ins]), rtol=1e-6)
    assert raw['ST|patch|%s|held' % E.pkey('m2', (1000, 0))].shape == (0, 4)
    # frame-level truth comparison: the mean of the two psField cores is the reference
    ref_mean = np.mean([E.stamp_core(cache[s][0]['psf']) for s in sets], 0)
    assert np.array_equal(raw['ST|frame|m2|all|psf_mean_core'], ref_mean)
    tr = raw['ST|frame|m2|all|truth_vs_psfmean']
    assert tr.shape == (4,) and np.isfinite(tr).all() and raw['ST|frame|m2|all|insitu_vs_psfmean'][0] != tr[0]
    # a patch with only one of its two scopes stored is completed on the next call
    cut = {k: v for k, v in raw.items() if not k.startswith('ST|frame|%s|' % E.pkey('m2', (0, 0)))}
    p2 = str(tmp_path / 'cut.npz')
    np.savez_compressed(p2, **cut)
    E.part_stamps(p2, B=5)
    assert 'ST|frame|%s|stamp' % E.pkey('m2', (0, 0)) in E.load_raw(p2)


def test_the_stamp_arrays_count_the_stars_that_were_kept_not_the_stars_selected(planted):
    st, dn, wrong = planted
    stars = E.collect([dict(dn=dn, gain=GAIN, dark=DARK)], wrong)
    c = stars['c'].copy()
    c[2, 11:14, 5:8] += 0.8 * c[2].max()
    est = E.estimate(c, stars['off'], wrong)
    est['stars'] = dict(stars, c=c)
    assert est['use'].sum() < len(c)
    arr = E.stamp_arrays('P', est, E.stamp_core(wrong))
    assert arr['P|n_sel'][0] == len(c) and arr['P|n_used'][0] == int(est['use'].sum()) < arr['P|n_sel'][0]
    assert np.array_equal(arr['P|vs_psf'], E.compare_cores(E.stamp_core(est['stamp']), E.stamp_core(wrong)))
    assert list(arr['P|star_use']) == [int(u) for u in est['use']]


def test_the_control_comparison_is_made_against_the_psfield_core(monkeypatch, tmp_path):
    st, wrong, cache = planted_world(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_stamps(path, B=4)
    E.part_scr(path)
    r = E.load_raw(path)
    k = 'ST|scr|%s|' % E.pkey('m2', (0, 0))
    assert np.array_equal(r[k + 'vs_psf'], E.compare_cores(E.stamp_core(r[k + 'stamp']), E.stamp_core(wrong)))
    assert not np.array_equal(r[k + 'vs_psf'], E.compare_cores(E.stamp_core(r[k + 'stamp']), E.stamp_core(r['ST|patch|%s|stamp' % E.pkey('m2', (0, 0))])))


def test_the_report_skips_patches_without_a_stamp_in_the_mean_and_the_finding_table_pools_the_counts_of_every_patch(planted_raw):
    raw, _ = planted_raw
    k = 'ST|patch|%s|vs_psf' % E.pkey('m2', (0, 0))
    bad = dict(raw)
    bad[k] = np.full(4, np.nan)
    t = E.report(bad)
    other = raw['ST|patch|%s|vs_psf' % E.pkey('m2', (1000, 0))]
    assert 'Mean over the 1 m2 patches with a stamp: central pixel %s %%' % E.pct(other[0]) in t
    pooled = dict(raw)
    pooled['CSD|SEL|finding'] = np.array([1])
    for i, p in enumerate(E.X.TESTS):
        pooled['CSTEST|cs3_base|x%d_y%d|tot' % p] = np.array([100, 10 * (i + 1), 100, 10 * (i + 1), 0, 0, 0])
        pooled['CSTEST|cs3_blendthreshu4|x%d_y%d|tot' % p] = np.array([100, 20 * (i + 1), 100, 20 * (i + 1), 0, 0, 0])
    t2 = E.report(pooled)
    assert 'pooled counts | 0.250 | 0.500' in t2


def test_the_held_out_misfit_is_the_residual_rms_over_the_fitted_peak_of_the_model(monkeypatch):
    st = lumpy_stamp(1.4, 1.2, 0.0)
    img = planted_frame(st, [(50.0, 50.0, 1e5)], noise=False, sky=0.0)

    class Sol:
        x = np.array([2.0e4, 0.1, -0.2])
        fun = np.full(121, 3.0)
    monkeypatch.setattr(E, 'least_squares', lambda f, x0, x_scale=None: Sol)
    r = E.hp_fit_rms(img, GAIN, DARK, st, [(50.0, 50.0)])
    yh = E.HP.highpass_frame(img, 3.0)
    d = yh[45:56, 45:56]
    sh = E.nd_shift(st, (-0.2, 0.1), order=3, mode='constant', cval=0.0)
    m = 2.0e4 * E.HP.highpass_stamp(sh, 3.0)[20:31, 20:31]
    assert r[0, 0] == pytest.approx(3.0) and r[0, 1] == pytest.approx(np.sqrt(np.mean((d - m) ** 2)) / m.max(), rel=1e-12)
    assert r[0, 1] != pytest.approx(np.sqrt(np.mean((d - m) ** 2)) / d.max(), rel=1e-3)


def test_the_plain_stack_statistic_is_the_central_pixel_ratio_the_fwhm_and_the_star_count():
    st = lumpy_stamp(1.4, 1.4, 0.0)
    dn = planted_frame(st, [(30.0, 30.0, 1e5), (70.0, 60.0, 8e4)], noise=False, sky=0.0)
    core = np.median(E.D.cutouts(dn, np.array([30.0, 70.0]), np.array([30.0, 60.0]), np.array([0, 1])), 0)
    ref = core * 2.0
    v = E.plain_stack_stat(dn, np.array([30.0, 70.0]), np.array([30.0, 60.0]), np.array([0, 1]), ref)
    assert v[0] == pytest.approx(-0.5) and v[1] == pytest.approx(E.D.fwhm(core)) and v[2] == 2
    assert np.isnan(E.plain_stack_stat(dn, np.array([30.0]), np.array([30.0]), np.array([], int), ref)).all()


def test_the_frame_level_truth_comparison_compares_the_in_situ_frame_stamp_with_the_truth_star_stamp(monkeypatch, tmp_path):
    st, wrongs, cache, sets = world_with_pairs(monkeypatch)
    path = str(tmp_path / 'insitu_raw.npz')
    E.part_stamps(path, B=4)
    raw = E.load_raw(path)
    frames = [dict(dn=cache[s][0]['dn'], gain=GAIN, dark=DARK) for s in sets]
    wings = cache[sets[0]][0]['psf']
    ef = E.estimate(*[E.collect(frames, wings)[k] for k in ('c', 'off')], wings)
    cs, offs = [], []
    for s in sets:
        patch, sky, truth = cache[s]
        tx, ty, tm, tidx = E.truth_stars(truth, s[1][0], s[1][1], patch['dn'])
        c, o = E.star_cores(patch['dn'], tx, ty, tidx)
        cs.append(c)
        offs.append(o)
    et = E.estimate(np.concatenate(cs), np.concatenate(offs), wings)
    ref = np.mean([E.stamp_core(cache[s][0]['psf']) for s in sets], 0)
    k = 'ST|frame|m2|all|'
    assert np.array_equal(raw[k + 'insitu_vs_truth'], E.compare_cores(E.stamp_core(ef['stamp']), E.stamp_core(et['stamp'])))
    assert np.array_equal(raw[k + 'truth_vs_psfmean'], E.compare_cores(E.stamp_core(et['stamp']), ref)) and np.array_equal(raw[k + 'insitu_vs_psfmean'], E.compare_cores(E.stamp_core(ef['stamp']), ref))
    assert raw[k + 'truth_n'][0] == len(np.concatenate(cs)) and raw[k + 'truth_n_used'][0] == int(et['use'].sum())


def test_part_run_runs_only_the_patches_of_the_sets_asked_for(monkeypatch, tmp_path):
    sets = (('m2', (0, 0)), ('m53', (0, 0)))
    st, wrong, cache = planted_world(monkeypatch, sets)
    path = str(tmp_path / 'insitu_raw.npz')
    monkeypatch.setattr(E, 'frozen_config', lambda: (3.0, 1.2, 2.0))
    monkeypatch.setattr(E, 'round2_specs', lambda: (None, None))
    E.part_stamps(path, B=4)
    seen = []

    def fake_base(meth, patch, sky_mean, t, specs):
        seen.append(float(patch['psf'][25, 25]))
        blk, tot = fake_scores(len(seen))
        return {'k2|tot': tot, 'k2|blk22': blk, 'det': np.zeros((1, 3)), 'secs': np.array([1.0])}
    monkeypatch.setattr(E, 'run_baseline', fake_base)
    E.part_run('patch', ['cs_def'], ['m53'], path)
    raw = E.load_raw(path)
    assert len(seen) == 1 and 'RUN|patch|%s|cs_def|wall' % E.pkey('m53', (0, 0)) in raw and 'RUN|patch|%s|cs_def|wall' % E.pkey('m2', (0, 0)) not in raw
