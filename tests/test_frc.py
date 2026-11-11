"""Tests of experiments/frc.py (Fourier ring correlation resolution of a localization set) on synthetic data with known answers: the ring
statistics against a brute-force loop, the frequency -> resolution mapping against an analytic signal + noise model, the resolution of
localization sets blurred by known Gaussian errors (monotone in the error, of the right order), identical halves (no crossing), uniform
random points (lowest-frequency crossing or none), the edge window, the splits and the block bootstrap (determinism, no localization in
both halves). Fast: fields of 3.2 um at 10 nm, one 512 x 512 image."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
import frc as F        # noqa: E402

FIELD = 3200.0


def structure(rng, n):
    """n points on a known structure in the 3.2 um field: two lines and a ring (uniform in arc length), no localization error yet."""
    kind = rng.integers(0, 3, n)
    t = rng.uniform(0, 1, n)
    x, y = np.empty(n), np.empty(n)
    m = kind == 0
    x[m], y[m] = 300 + 2600 * t[m], 400 + 600 * t[m]
    m = kind == 1
    x[m], y[m] = 400 + 2000 * t[m], 2800 - 300 * t[m]
    m = kind == 2
    ang = 2 * np.pi * t[m]
    x[m], y[m] = 2200 + 500 * np.cos(ang), 1500 + 500 * np.sin(ang)
    return np.c_[x, y]


def halves(rng, sigma, n=20000):
    """Two independent half sets: n draws from the structure each, with a Gaussian localization error of sd sigma nm per axis."""
    s1, s2 = structure(rng, n), structure(rng, n)
    return s1 + rng.normal(0, sigma, s1.shape), s2 + rng.normal(0, sigma, s2.shape)


# ------------------------------------------------------------------ rendering and the ring statistics
def test_render_puts_each_localization_in_its_pixel_and_ignores_the_outside():
    xy = np.array([[5.0, 5.0], [15.0, 25.0], [15.0, 25.0], [3199.9, 0.0], [3200.0, 10.0], [-0.1, 10.0], [np.nan, 10.0], [10.0, 3200.0]])
    img = F.render(xy, FIELD, 10.0)
    assert img.shape == (320, 320) and img.sum() == 4
    assert img[0, 0] == 1 and img[2, 1] == 2 and img[0, 319] == 1                 # (row, column) = (y, x) pixel indices
    w = F.render(xy[:3], FIELD, 10.0, weights=[1.0, 2.0, 3.0])
    assert w[0, 0] == 1 and w[2, 1] == 5
    assert F.render(np.zeros((0, 2)), FIELD, 10.0).sum() == 0 and F.render([], FIELD, 10.0).shape == (320, 320)
    assert F.render(np.c_[xy, np.arange(len(xy))], FIELD, 10.0).sum() == 4        # extra columns are ignored
    assert list(F.inside_mask(xy, FIELD)) == [True, True, True, True, False, False, False, False]


def _brute_frc(a, b, alpha=0.25, subtract_mean=True):
    n = a.shape[0]
    w = F.window2d(n, alpha)
    if subtract_mean:
        a, b = a - (w * a).sum() / w.sum(), b - (w * b).sum() / w.sum()
    fa, fb = np.fft.fft2(a * w), np.fft.fft2(b * w)
    k = np.fft.fftfreq(n) * n
    ring = np.rint(np.sqrt(k[:, None] ** 2 + k[None, :] ** 2)).astype(int)
    out = []
    for r in range(1, n // 2 + 1):
        m = ring == r
        out.append((fa[m] * np.conj(fb[m])).real.sum() / np.sqrt((abs(fa[m]) ** 2).sum() * (abs(fb[m]) ** 2).sum()))
    return np.array(out)


@pytest.mark.parametrize('n', [16, 17, 33])
def test_ring_correlation_equals_a_brute_force_loop_on_the_full_fft(n):
    rng = np.random.default_rng(n)
    a = rng.random((n, n))
    b = a + rng.random((n, n))
    c = F.frc_curve(a, b, 10.0)
    assert np.abs(c['frc'] - _brute_frc(a, b)).max() < 1e-12
    assert np.allclose(c['freq'], np.arange(1, n // 2 + 1) / (n * 10.0))              # cycles per nm; the last ring is the Nyquist frequency
    assert np.abs(F.frc_curve(a, b, 10.0, alpha=0.0, subtract_mean=False)['frc'] - _brute_frc(a, b, 0.0, False)).max() < 1e-12
    assert c['npix'].sum() <= n * n and (c['npix'] > 0).all()


def test_frequency_to_resolution_mapping_against_an_analytic_signal_plus_noise_model():
    """Signal with the Gaussian power spectrum P_s(f) = h0^2 exp(-f^2 / (2 f0^2)) shared by both images, independent white noise of variance 1:
    the expected ring correlation is P_s / (P_s + P_n), which falls through 1/7 at f* = f0 sqrt(2 ln(6 h0^2)); the resolution is 1 / f*."""
    n = 512
    kx, ky = np.fft.fftfreq(n)[None, :], np.fft.fftfreq(n)[:, None]
    f2 = kx ** 2 + ky ** 2
    for seed, f0, snr0 in ((1, 0.03, 20.0), (2, 0.02, 20.0), (3, 0.04, 5.0)):
        rng = np.random.default_rng(seed)
        H = np.sqrt(snr0) * np.exp(-f2 / (4 * f0 ** 2))
        s = np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * H).real
        a, b = s + rng.standard_normal((n, n)), s + rng.standard_normal((n, n))
        pixel = 10.0
        r = F.frc_images(a, b, pixel)
        predicted = pixel / (f0 * np.sqrt(2 * np.log(6 * snr0)))
        assert r['crossed'] and abs(r['resolution_nm'] / predicted - 1) < 0.08, (seed, r['resolution_nm'], predicted)
        assert r['frc'][0] > 0.8 and r['frc'][-1] < 0.1                                 # high at the lowest rings, ~0 at the Nyquist frequency


def test_threshold_interpolation_smoothing_and_the_0143_variant():
    k = np.arange(1, 101)
    curve = dict(ring=k, freq=k / 1000.0, frc=1.0 - (k - 1) / 99.0)                     # straight line from 1 at ring 1 to 0 at ring 100
    r = F.frc_resolution(curve, pixel_nm=5.0, smooth=1)
    k_star = 1 + 99.0 * (1.0 - 1.0 / 7.0)                                               # the exact ring where the line crosses 1/7
    assert r['crossed'] and abs(r['resolution_nm'] - 1000.0 / k_star * 1.0) < 1e-9
    assert abs(r['resolution_nm'] - 1000.0 / (k_star / 1000.0) / 1000.0) < 1e-9
    r143 = F.frc_resolution(curve, pixel_nm=5.0, smooth=1, threshold=0.143)
    assert abs(r143['resolution_nm'] / r['resolution_nm'] - 1) < 5e-4                    # "1/7" and "0.143" are the same number in practice
    assert r['limit_nm'] == 10.0
    # a moving average of a line does not move it (except at the ends)
    assert abs(F.frc_resolution(curve, 5.0, smooth=7)['resolution_nm'] - r['resolution_nm']) < 1e-9
    assert np.allclose(F.smooth_curve([0, 0, 6, 0, 0], 3), [0, 2, 2, 2, 0]) and np.allclose(F.smooth_curve([1, 2, 3], 1), [1, 2, 3])
    # first and last crossing of a wandering curve
    y = np.array([1.0, 0.9, 0.5, 0.05, 0.3, 0.05, 0.02])
    w = F.frc_resolution(dict(freq=np.arange(1, 8) / 100.0, frc=y), 10.0, smooth=1)
    assert w['crossed'] and w['resolution_nm'] > w['last_nm'] > 0          # first crossing at a lower frequency (larger nm) than the last one


def test_degenerate_curves_and_inputs():
    k = np.arange(1, 11)
    below = F.frc_resolution(dict(freq=k / 100.0, frc=np.full(10, 0.01)), 10.0, smooth=1)
    assert below['crossed'] and below['at_first_ring'] and abs(below['resolution_nm'] - 100.0) < 1e-9      # no correlation at all: the field size
    never = F.frc_resolution(dict(freq=k / 100.0, frc=np.full(10, 0.9)), 10.0, smooth=1)
    assert not never['crossed'] and np.isnan(never['resolution_nm']) and never['limit_nm'] == 20.0
    empty = F.frc_resolution(dict(freq=k / 100.0, frc=np.full(10, np.nan)), 10.0)
    assert not empty['crossed'] and np.isnan(empty['resolution_nm'])
    r = F.frc_localizations(np.zeros((0, 2)), np.zeros(0, int), FIELD, 10.0, split='random')          # an empty set has no curve
    assert not r['crossed'] and np.isnan(r['resolution_nm']) and r['n_a'] == r['n_b'] == 0
    with pytest.raises(ValueError):
        F.frc_curve(np.zeros((8, 9)), np.zeros((8, 9)))
    with pytest.raises(ValueError):
        F.frc_curve(np.zeros((8, 8)), np.zeros((9, 9)))
    with pytest.raises(ValueError):
        F.grid_size(30.0, 10.0)


# ------------------------------------------------------------------ the four properties asked for
def test_resolution_is_monotone_in_the_localization_error_and_of_the_right_order():
    """(a) Two halves of a localization set drawn from lines and a ring, blurred by Gaussian errors of sd sigma: the FRC resolution grows with sigma
    and is a few sigma (a few FWHM / 1.5 ... 3 sigma, the usual relation), here 2 .. 4 sigma at 20,000 localizations per half."""
    rng = np.random.default_rng(1)
    sigmas = [10.0, 15.0, 25.0, 40.0, 60.0]
    res = []
    for s in sigmas:
        a, b = halves(rng, s)
        r = F.frc_images(F.render(a, FIELD, 10.0), F.render(b, FIELD, 10.0), 10.0)
        assert r['crossed']
        res.append(r['resolution_nm'])
    assert all(x < y for x, y in zip(res, res[1:])), res
    for s, x in zip(sigmas, res):
        assert 1.8 * s < x < 4.5 * s, (s, x)
    # more localizations never make it worse at the same error (same structure, same sigma): 5,000 -> 80,000 per half
    a5, b5 = halves(np.random.default_rng(5), 15.0, 5000)
    a8, b8 = halves(np.random.default_rng(6), 15.0, 80000)
    r5 = F.frc_images(F.render(a5, FIELD, 10.0), F.render(b5, FIELD, 10.0), 10.0)['resolution_nm']
    r8 = F.frc_images(F.render(a8, FIELD, 10.0), F.render(b8, FIELD, 10.0), 10.0)['resolution_nm']
    assert r8 < r5


def test_identical_halves_never_cross_the_threshold():
    """(b) The same image twice: the correlation is 1 in every ring; no crossing, resolution NaN, the pixel limit reported."""
    rng = np.random.default_rng(2)
    a, _ = halves(rng, 10.0)
    img = F.render(a, FIELD, 10.0)
    r = F.frc_images(img, img.copy(), 10.0)
    assert np.allclose(r['frc'], 1.0, atol=1e-9)
    assert not r['crossed'] and np.isnan(r['resolution_nm']) and r['limit_nm'] == 20.0
    xy = np.c_[a, np.arange(len(a)) // 7]
    both = F.frc_images(F.render(xy, FIELD, 10.0), F.render(xy, FIELD, 10.0), 10.0)
    assert not both['crossed']


def test_uniform_random_points_cross_at_the_lowest_frequencies_or_not_at_all():
    """(c) No structure: the halves correlate at no frequency, so the curve is at noise level (mean ~ 0) and, if the smoothed curve falls below 1/7 at
    all, it does so within the first rings (a 'resolution' of the order of the field size), never at a fine scale."""
    rng = np.random.default_rng(3)
    for n in (5000, 40000):
        u1, u2 = rng.uniform(0, FIELD, (n, 2)), rng.uniform(0, FIELD, (n, 2))
        r = F.frc_images(F.render(u1, FIELD, 10.0), F.render(u2, FIELD, 10.0), 10.0)
        assert abs(np.nanmean(r['frc'])) < 0.02
        assert (not r['crossed']) or r['resolution_nm'] >= FIELD / 4, r['resolution_nm']


def test_the_split_is_a_function_of_the_seed_alone():
    """(d) Determinism: the same seed gives the same halves and the same numbers; another seed another split; odd/even does not use the seed."""
    rng = np.random.default_rng(4)
    xy = structure(rng, 30000) + rng.normal(0, 20.0, (30000, 2))
    frames = rng.integers(0, 200, 30000)
    for mode in F.SPLITS:
        m1, m2, m3 = (F.split_mask(len(xy), mode, s, frames, 10) for s in (7, 7, 8))
        assert (m1 == m2).all()
        assert abs(int(m1.sum()) - len(xy) / 2) < 0.12 * len(xy) or mode != 'random'
        if mode == 'oddeven':
            assert (m1 == m3).all() and (m1 == (frames % 2 == 1)).all()
        elif mode == 'alternating':
            assert (m1 == m3).all() and (m1 == ((frames // 10) % 2 == 1)).all()                # no randomness: blocks of 10 frames alternately to the halves
        else:
            assert not (m1 == m3).all()
    m = F.split_mask(len(xy), 'random', 7)
    assert m.sum() == len(xy) // 2                                                  # exactly balanced
    a = F.frc_localizations(xy, frames, FIELD, 10.0, 'random', 7)
    b = F.frc_localizations(xy, frames, FIELD, 10.0, 'random', 7)
    c = F.frc_localizations(xy, frames, FIELD, 10.0, 'random', 8)
    assert a['resolution_nm'] == b['resolution_nm'] and np.array_equal(a['frc'], b['frc'])
    assert a['resolution_nm'] != c['resolution_nm'] and not np.array_equal(a['frc'], c['frc'])
    assert a['n_a'] + a['n_b'] == len(xy) and a['n_outside'] == 0 or a['n_a'] + a['n_b'] + a['n_outside'] == len(xy)
    for mode in ('oddeven', 'blocks', 'alternating'):
        r = F.frc_localizations(xy, frames, FIELD, 10.0, mode, 3)
        assert r['crossed'] and 20 < r['resolution_nm'] < 400
    with pytest.raises(ValueError):
        F.split_mask(10, 'oddeven')                                                 # needs frame numbers
    with pytest.raises(ValueError):
        F.split_mask(10, 'nonsense')


def test_block_split_keeps_every_frame_of_a_block_in_one_half():
    frames = np.repeat(np.arange(100), 3)
    mask = F.split_mask(len(frames), 'blocks', 5, frames, block=10)
    blocks = frames // 10
    for b in range(10):
        assert len(set(mask[blocks == b])) == 1
    assert sum(mask[blocks == b][0] for b in range(10)) == 5                         # five of the ten blocks in each half


# ------------------------------------------------------------------ the window
def test_the_window_removes_the_edge_ridge_of_a_structure_cut_by_the_border():
    """Two images with the same strong ramp (a structure running off the field) and independent noise: without a window the wrap-around step puts a
    ridge along the frequency axes into both images and the rings correlate at every frequency; with the Tukey window they do not."""
    n = 256
    rng = np.random.default_rng(8)
    ramp = np.linspace(0, 1, n)[None, :] * np.ones((n, 1)) * 2000.0
    a, b = ramp + rng.normal(0, 5, (n, n)), ramp + rng.normal(0, 5, (n, n))
    raw = F.frc_curve(a, b, 10.0, alpha=0.0)['frc'][40:100].mean()
    win = F.frc_curve(a, b, 10.0, alpha=0.25)['frc'][40:100].mean()
    assert raw > 0.9 and win < 0.1
    w = F.window2d(n, 0.25)
    assert w.shape == (n, n) and w[0, 0] == 0 and w[n // 2, n // 2] == 1 and w[n // 2, 0] == 0
    assert np.allclose(F.window2d(n, 0.0), 1.0)


# ------------------------------------------------------------------ the block bootstrap
def test_bootstrap_multiplicities_and_halves_never_share_a_block():
    mult = F.bootstrap_draw(50, seed=3, n_boot=40)
    assert mult.shape == (40, 50) and (mult.sum(axis=1) == 50).all()                # as many draws as blocks
    assert np.array_equal(mult, F.bootstrap_draw(50, seed=3, n_boot=40)) and not np.array_equal(mult, F.bootstrap_draw(50, seed=4, n_boot=40))
    for r in range(40):
        half = F.bootstrap_halves(mult[r], np.random.default_rng([3, r]))
        used = mult[r] > 0
        assert not half[~used].any()                                                # an undrawn block is in neither half
        assert abs(int(half[used].sum()) - int((~half[used]).sum())) <= 1           # halves of (nearly) equal numbers of blocks
    ids, uniq, inv = F.block_ids(np.array([0, 3, 9, 10, 25, 99]), 10)
    assert list(uniq) == [0, 1, 2, 9] and list(inv) == [0, 0, 0, 1, 2, 3]
    # with n_frames the blocks are those of the whole sequence, the same for every localization set of it (paired resamples)
    _, uniq_all, inv_all = F.block_ids(np.array([0, 3, 9, 10, 25, 99]), 10, n_frames=100)
    assert list(uniq_all) == list(range(10)) and list(inv_all) == [0, 0, 0, 1, 2, 9]
    _, uniq_sub, inv_sub = F.block_ids(np.array([0, 3]), 10, n_frames=100)
    assert list(uniq_sub) == list(range(10))
    with pytest.raises(ValueError):
        F.block_ids(np.array([100]), 10, n_frames=100)


def test_bootstrap_resolution_is_deterministic_and_centred_near_the_point_estimate():
    rng = np.random.default_rng(10)
    xy = structure(rng, 40000) + rng.normal(0, 20.0, (40000, 2))
    frames = rng.integers(0, 300, 40000)
    a = F.bootstrap_blocks(xy, frames, n_boot=12, seed=1, block=10, field_nm=FIELD)
    b = F.bootstrap_blocks(xy, frames, n_boot=12, seed=1, block=10, field_nm=FIELD)
    c = F.bootstrap_blocks(xy, frames, n_boot=12, seed=2, block=10, field_nm=FIELD)
    assert a['n_blocks'] == 30 and a['resolution_nm'].shape == (12,)
    assert np.array_equal(a['resolution_nm'], b['resolution_nm']) and not np.array_equal(a['resolution_nm'], c['resolution_nm'])
    assert a['crossed'].all()
    point = F.frc_localizations(xy, frames, FIELD, 10.0, 'blocks', 0)['resolution_nm']
    assert abs(np.median(a['resolution_nm']) / point - 1) < 0.25                     # the replicates use ~63 % of the blocks, so a somewhat larger value is expected


def test_bootstrap_replicates_are_coarser_than_the_estimate_on_the_full_sample_and_paired_sets_share_their_resamples():
    """A replicate holds ~63 % of the distinct blocks, so its halves hold fewer localizations and the resolution is coarser than on the full sample: the percentile
    interval of the replicates does not bracket the full-sample estimate (documented in ``bootstrap_blocks``). Two sets of the same frames resampled with
    ``n_frames`` get identical resamples, so the ratio of their resolutions is a paired quantity."""
    rng = np.random.default_rng(12)
    xy = structure(rng, 60000) + rng.normal(0, 20.0, (60000, 2))
    frames = rng.integers(0, 400, 60000)
    full = F.frc_localizations(xy, frames, FIELD, 10.0, 'alternating')['resolution_nm']
    boot = F.bootstrap_blocks(xy, frames, n_boot=16, seed=2, block=10, field_nm=FIELD, n_frames=400)['resolution_nm']
    assert np.percentile(boot, 5) > full                                                  # every replicate (but at most one) coarser than the estimate on all the data
    # a second set of the same frames (here: the same localizations with a different localization error) is resampled identically; its ratio to the first is stable
    xy2 = xy + rng.normal(0, 20.0, xy.shape)                                              # a method with a larger error
    boot2 = F.bootstrap_blocks(xy2, frames, n_boot=16, seed=2, block=10, field_nm=FIELD, n_frames=400)['resolution_nm']
    ratio = boot2 / boot
    full2 = F.frc_localizations(xy2, frames, FIELD, 10.0, 'alternating')['resolution_nm']
    assert (ratio > 1).all() and abs(np.median(ratio) / (full2 / full) - 1) < 0.12


def test_bootstrap_relative_spread_is_at_least_the_spread_across_independent_datasets():
    """The report applies the relative spread of the bootstrap replicates to the full-sample estimate and calls it conservative: on independent synthetic datasets of the
    same size the spread of the alternating-block estimate is smaller than the relative spread of the bootstrap replicates of one dataset."""
    def dataset(seed, n=30000, nfr=300):
        rng = np.random.default_rng(seed)
        return structure(rng, n) + rng.normal(0, 20.0, (n, 2)), rng.integers(0, nfr, n)
    est = np.array([F.frc_localizations(*dataset(100 + m), FIELD, 10.0, 'alternating')['resolution_nm'] for m in range(12)])
    xy, fr = dataset(100)
    boot = F.bootstrap_blocks(xy, fr, n_boot=24, seed=1, block=10, field_nm=FIELD, n_frames=300)['resolution_nm']
    assert np.std(boot, ddof=1) / np.median(boot) > np.std(est, ddof=1) / np.mean(est)


def test_repeated_localizations_in_both_halves_inflate_the_resolution_and_a_block_split_does_not():
    """A sparse set of fluorophore sites, each localized in 60 consecutive frames with independent errors (a blinking dye that is on for a long time): a
    random split (and the odd/even split) puts the repeats of every site into BOTH halves, so the halves share the sites and the resolution looks
    much finer than it is; a split by blocks of frames puts a site into one half only and reads the structure's own resolution."""
    rng = np.random.default_rng(11)
    n_sites, rep, sigma = 60, 60, 15.0
    sites = structure(rng, n_sites)
    start = rng.integers(0, 100, n_sites) * rep
    xy = np.concatenate([sites + rng.normal(0, sigma, sites.shape) for _ in range(rep)])
    frames = np.concatenate([start + k for k in range(rep)])
    res = {m: np.median([F.frc_localizations(xy, frames, FIELD, 10.0, m, s, block=rep)['resolution_nm'] for s in range(5)]) for m in ('random', 'oddeven', 'blocks')}
    assert res['random'] < 0.6 * res['blocks'] and res['oddeven'] < 0.6 * res['blocks'], res
