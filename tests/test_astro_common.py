"""Tests of experiments/astro_common.py: magnitudes, the truth list of a patch, the one-to-one matching at 0.5 px / 0.5 mag, completeness and purity by magnitude."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
import astro_common as A      # noqa: E402


def test_magnitude_conversions_invert_each_other_and_handle_nonpositive_flux():
    f = np.array([100.0, 1e4, 3e5])
    assert np.allclose(A.mag_to_counts(A.counts_to_mag(f)), f, rtol=1e-12)
    assert abs(A.counts_to_mag(1.0 / 0.00546689) - 22.5) < 1e-12 and abs(A.counts_to_mag(1.0 / 0.00546689 * 100) - 17.5) < 1e-12
    m = A.counts_to_mag(np.array([-1.0, 0.0, 5.0]))
    assert np.isnan(m[0]) and np.isnan(m[1]) and np.isfinite(m[2])


def test_truth_for_patch_box_edges_and_inner_mask():
    truth = dict(px=np.array([309.4, 310.0, 359.5, 409.4, 409.6, 312.4, 312.6]), py=np.array([680.0] * 7), m606=np.array([20.0, 21, 22, 23, 24, 25, 26]))
    t = A.truth_for_patch(truth, 310, 630, 100)
    assert len(t['x']) == 5                                              # 309.4 is outside (< 309.5) and 409.6 is outside (>= 409.5)
    assert np.allclose(t['x'], [0.0, 49.5, 99.4, 2.4, 2.6])
    # patch pixels: x - 310; inner = at least 2.5 px from the edges: x in [2.0, 97.0]
    inner = {round(float(x), 2): bool(i) for x, i in zip(t['x'], t['inner'])}
    assert inner[0.0] is False and inner[49.5] is True and inner[99.4] is False and inner[2.4] is True
    t2 = A.truth_for_patch(truth, 310, 630, 100, mag_max=22.5)
    assert len(t2['x']) == 2                                             # magnitudes 21 and 22 of the five in the box


def test_match_is_one_to_one_within_both_tolerances_and_maximal():
    tr_xy = np.array([[10.0, 10.0], [10.4, 10.0]])
    tr_mag = np.array([20.0, 20.2])
    # two detections near both truths: both can match
    r, c = A.match(np.array([[10.0, 10.0], [10.4, 10.0]]), np.array([20.0, 20.2]), tr_xy, tr_mag)
    assert sorted(zip(r.tolist(), c.tolist())) == [(0, 0), (1, 1)]
    # one detection between them takes only one truth; a second detection of the same star is a false positive
    r, c = A.match(np.array([[10.2, 10.0], [10.2, 10.0]]), np.array([20.1, 20.1]), tr_xy, tr_mag)
    assert len(r) == 2 and len(set(c.tolist())) == 2
    r, c = A.match(np.array([[10.2, 10.0]]), np.array([20.1]), tr_xy, tr_mag)
    assert len(r) == 1
    # the magnitude tolerance and the position tolerance both bind
    assert len(A.match(np.array([[10.0, 10.0]]), np.array([20.6]), tr_xy[:1], tr_mag[:1])[0]) == 0
    assert len(A.match(np.array([[10.6, 10.0]]), np.array([20.0]), tr_xy[:1], tr_mag[:1])[0]) == 0
    assert len(A.match(np.array([[10.5, 10.0]]), np.array([20.5]), tr_xy[:1], tr_mag[:1])[0]) == 1      # inclusive
    assert A.match(np.zeros((0, 2)), np.zeros(0), tr_xy, tr_mag)[0].size == 0
    # maximal: a greedy nearest-first choice would match one pair, the optimum two
    r, c = A.match(np.array([[10.0, 10.0], [10.45, 10.0]]), np.array([20.0, 20.1]), np.array([[10.2, 10.0], [10.5, 10.0]]), np.array([20.0, 20.1]))
    assert len(r) == 2


def test_score_counts_bins_edges_and_purity():
    truth = dict(px=np.array([310 + 20.0, 310 + 40.0, 310 + 60.0, 310 + 1.0]), py=np.array([630 + 20.0, 630 + 40.0, 630 + 60.0, 630 + 50.0]),
                 m606=np.array([18.2, 19.7, 22.8, 18.0]))
    t = A.truth_for_patch(truth, 310, 630, 100)
    assert t['inner'].sum() == 3                                         # the source at x = 1 is within 2.5 px of the edge
    det_xy = np.array([[20.1, 20.0], [40.0, 40.2], [60.0, 60.0], [80.0, 80.0], [1.0, 50.0]])
    det_flux = A.mag_to_counts(np.array([18.3, 19.9, 22.7, 21.0, 18.0]))          # the 4th has no truth; the 5th is in the edge zone (not scored)
    s = A.score(det_xy, det_flux, t)
    assert s['true_all'] == 3 and s['matched_all'] == 3 and s['det_all'] == 4
    assert s['true22'] == 2 and s['matched22'] == 2 and s['det22'] == 3 and s['det22_matched'] == 2
    b = np.searchsorted(A.BIN_EDGES, 18.2, side='right') - 1
    assert s['n_true'][b] == 1 and s['n_matched'][b] == 1 and s['completeness'][b] == 1.0 and s['fp_fraction'][b] == 0.0
    b21 = np.searchsorted(A.BIN_EDGES, 21.0, side='right') - 1
    assert s['n_det'][b21] == 1 and s['fp_fraction'][b21] == 1.0 and np.isnan(s['completeness'][b21])
    assert np.allclose(s['offset'], [0.05, 0.1], atol=1e-9)                # median over the two matched pairs brighter than 21
    assert len(s['completeness']) == 17


def test_score_block_counts_add_up_and_follow_the_positions():
    truth = dict(px=np.array([310 + 5.0, 310 + 45.0, 310 + 45.0, 310 + 85.0, 310 + 85.0]), py=np.array([630 + 5.0, 630 + 25.0, 630 + 25.5, 630 + 85.0, 630 + 60.0]),
                 m606=np.array([19.0, 19.5, 21.0, 21.5, 23.0]))
    t = A.truth_for_patch(truth, 310, 630, 100)
    det_xy = np.array([[5.0, 5.0], [45.0, 25.0], [85.0, 85.0], [60.0, 10.0], [30.0, 90.0]])            # the 4th and 5th have no truth
    det_flux = A.mag_to_counts(np.array([19.1, 19.4, 21.6, 20.0, 21.0]))
    s = A.score(det_xy, det_flux, t)
    b = s['blk22']
    assert b.shape == (25, 4)
    assert np.allclose(b.sum(0), [s['true22'], s['matched22'], s['det22'], s['det22_matched']])
    assert s['true22'] == 4 and s['matched22'] == 3 and s['det22'] == 5 and s['det22_matched'] == 3         # the 23 mag truth is not counted
    cell = lambda x, y: (int(y // 20) * 5 + int(x // 20))
    assert b[cell(5, 5)].tolist() == [1, 1, 1, 1]
    assert b[cell(45, 25)].tolist() == [2, 1, 1, 1]                                                           # two truth sources in the block, one matched (one detection)
    assert b[cell(85, 85)].tolist() == [1, 1, 1, 1]
    assert b[cell(60, 10)].tolist() == [0, 0, 1, 0] and b[cell(30, 90)].tolist() == [0, 0, 1, 0]
    assert A.block_of(np.array([[-0.5, -0.5], [19.4, 19.4], [19.6, 19.6], [99.4, 99.4]])).tolist() == [0, 0, 6, 24]


def test_path_on_a_synthetic_crowded_patch_finds_the_sources_with_the_stopping_time():
    """A 48 x 48 patch whose PSF is a known sum of three Gaussians (the stamp given to the path is that PSF), 25 stars of 17-21.5 mag on a sky of 300 DN with
    Poisson-plus-dark noise at gain 4.6: the packaged path run (`exp_astro_m2.run_path`) recovers most of them with their magnitudes and positions and few
    false ones (the real PSF is not a sum of three Gaussians: its wings leave residuals around bright stars, which is a result of the experiment, not of this test)."""
    import exp_astro_m2 as X
    from src.core import imaging as IM
    parts_true = ((0.55, 1.0), (0.35, 2.0), (0.10, 4.5))
    Q = IM.PixelPSF(48, 48, 0.9, parts=parts_true)
    stamp = IM.PixelPSF(51, 51, 0.9, parts=parts_true).phi(np.array([[25.0, 25.0]]))[:, 0].reshape(51, 51)
    rng = np.random.default_rng(31)
    pos = []
    while len(pos) < 25:
        c = rng.uniform(4, 44, 2)
        if all(np.linalg.norm(c - q) > 3.0 for q in pos):
            pos.append(c)
    pos = np.array(pos)
    mag = rng.uniform(17.0, 21.5, 25)
    flux = A.mag_to_counts(mag)
    sky, gain, dark = 300.0, 4.6, 1.0
    mean = sky + Q.synth(pos, flux)
    dn = (rng.poisson(mean * gain).astype(float) / gain) + rng.normal(0, np.sqrt(dark) / gain, mean.shape)
    patch = dict(dn=dn, psf=stamp, gain=gain, dark=dark)
    out, secs, n, info = X.run_path(patch, 5.0, 1.2, kappas=(3.0, 5.0))
    assert 0 <= info['bic_index'] < info['n_blocks'] and info['n_blocks'] >= n
    truth = dict(x=pos[:, 0], y=pos[:, 1], mag=mag, inner=np.ones(25, bool))
    s3 = A.score(np.c_[out[3.0][0], out[3.0][1]], out[3.0][2], truth, n=48)
    s5 = A.score(np.c_[out[5.0][0], out[5.0][1]], out[5.0][2], truth, n=48)
    assert s3['matched_all'] >= 0.85 * 25 and s3['matched_all'] >= s5['matched_all'] and s5['det_all'] <= s3['det_all']
    assert s3['det_all'] - s3['matched_all'] <= 5 and np.all(np.abs(s3['offset']) < 0.1)
    assert n >= s3['matched_all']
    sb = A.score(np.c_[out[0.0][0], out[0.0][1]], out[0.0][2], truth, n=48)                                     # the BIC block of the same path (nothing to tune)
    assert sb['matched_all'] >= 0.8 * 25 and sb['det_all'] - sb['matched_all'] <= 5 and np.all(np.abs(sb['offset']) < 0.1)


def test_path_variants_on_the_synthetic_patch_change_what_they_claim_to_change():
    """The options of `run_path` behind the ablations: greedy order (the instantaneous correlation) finds the sources too on this sparse patch, 'snr' stops by the entrant's
    instantaneous significance, cleanup=False keeps the units whose own significance is below kappa (so never fewer), and hp = 0 / rex = 0 switch the high-pass and the exclusion off."""
    import exp_astro_m2 as X
    from src.core import imaging as IM
    parts_true = ((0.55, 1.0), (0.35, 2.0), (0.10, 4.5))
    Q = IM.PixelPSF(40, 40, 0.9, parts=parts_true)
    stamp = IM.PixelPSF(51, 51, 0.9, parts=parts_true).phi(np.array([[25.0, 25.0]]))[:, 0].reshape(51, 51)
    rng = np.random.default_rng(5)
    pos = np.array([[8.0, 9.0], [30.0, 12.0], [20.0, 30.0], [10.5, 28.0], [32.0, 33.0]])
    mag = np.array([18.0, 19.0, 17.5, 20.0, 19.5])
    mean = 300.0 + Q.synth(pos, A.mag_to_counts(mag))
    dn = (rng.poisson(mean * 4.6).astype(float) / 4.6) + rng.normal(0, 1.0 / 4.6, mean.shape)
    patch = dict(dn=dn, psf=stamp, gain=4.6, dark=1.0)
    truth = dict(x=pos[:, 0], y=pos[:, 1], mag=mag, inner=np.ones(5, bool))
    base = X.run_path(patch, 5.0, 1.2, kappas=(4.0,))[0][4.0]
    res = {}
    for name, kw in (('greedy', dict(order='greedy')), ('snr', dict(stop='snr')), ('nocleanup', dict(cleanup=False))):
        res[name] = X.run_path(patch, 5.0, 1.2, kappas=(4.0,), **kw)[0][4.0]
    for name, d in list(res.items()) + [('base', base)]:
        s = A.score(np.c_[d[0], d[1]], d[2], truth, n=40)
        assert s['matched_all'] >= 4, name                                                  # every variant finds at least 4 of the 5 bright, separated stars
    assert len(res['nocleanup'][0]) >= len(base[0])
    nohp = X.run_path(patch, 0.0, 1.2, kappas=(4.0,))[0][4.0]
    norex = X.run_path(patch, 5.0, 0.0, kappas=(4.0,))[0][4.0]
    assert len(nohp[0]) >= 4 and len(norex[0]) >= 4


def test_published_daophot_catalogue_loads_with_the_recorded_hash_and_the_documented_shift():
    import pytest
    import exp_astro_baselines as XB
    if not os.path.exists(XB.DAOPHOT):
        pytest.skip('the shipped DAOPHOT catalogue is not in $ASTRO_DATA')
    dp = XB.load_daophot()
    raw = np.loadtxt(XB.DAOPHOT)
    assert len(dp['x']) == 356 and np.allclose(dp['x'], raw[:, 0]) and np.allclose(dp['y'], raw[:, 1] - 1.0) and np.allclose(dp['flux'], raw[:, 2])
    z = np.load(os.path.join(XB.ROOT, 'truth_hugs_meth1.npz'))
    t = A.truth_for_patch(dict(px=z['px'], py=z['py'], m606=z['m606']), 310, 630)
    s = A.score(np.c_[dp['x'], dp['y']], dp['flux'], t)
    assert s['matched_all'] > 200 and abs(s['offset'][0]) < 0.3 and abs(s['offset'][1]) < 0.3                         # aligned: the sources match the truth
    shifted = A.score(np.c_[dp['x'], dp['y'] + 1.0], dp['flux'], t)
    assert shifted['matched_all'] < 20                                                                              # the file's own convention would not
    assert XB.load_daophot(os.path.join(XB.ROOT, 'nonexistent.txt')) is None
