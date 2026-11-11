"""Tests of the pure functions of experiments/exp_smlm_iss.py: the label-free PSF width, the BIC part's keys."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_smlm_iss as X      # noqa: E402


def _frames(n, sigma, rng, density=3, photons=(3000, 8000), bg=100.0):
    """n frames of 64 x 64 with `density` isolated emitters each (pixel-integrated Gaussians of width sigma), Poisson noise on a flat background."""
    out = np.empty((n, 64, 64))
    for i in range(n):
        f = np.full((64, 64), bg)
        pos = []
        while len(pos) < density:
            c = rng.uniform(8, 56, 2)
            if all(np.hypot(*(c - q)) > 14 for q in pos):
                pos.append(c)
        for c in pos:
            yy, xx = np.mgrid[-5:6, -5:6]
            cy, cx = int(round(c[1])), int(round(c[0]))
            f[cy - 5:cy + 6, cx - 5:cx + 6] += rng.uniform(*photons) * X._pix_gauss(c[0] - cx, c[1] - cy, sigma, 5)
        out[i] = rng.poisson(f).astype(float)
    return out


def test_label_free_psf_width_recovers_the_planted_width_and_is_nan_without_emitters():
    rng = np.random.default_rng(3)
    bg = np.full((64, 64), 100.0)
    for sigma in (1.2, 1.6):
        fr = _frames(40, sigma, rng)
        est, n, q25, q75 = X.estimate_sigma0(fr, bg, 1.0, b0=0.0)
        assert n > 40 and abs(est - sigma) < 0.04 * sigma and q25 <= est <= q75
    empty = rng.poisson(np.full((20, 64, 64), 100.0)).astype(float)
    est, n, q25, q75 = X.estimate_sigma0(empty, bg, 1.0, b0=0.0)
    assert np.isnan(est) and n == 0


def test_label_free_psf_width_is_not_corrupted_by_blends_and_skips_border_spots():
    rng = np.random.default_rng(4)
    bg = np.full((64, 64), 100.0)
    sigma = 1.4
    fr = _frames(40, sigma, rng)
    for i in range(40):                                                      # a blend of two emitters 3 px apart (unimodal: one local maximum) in every frame, brighter than the singles
        for cx in (30, 33):
            fr[i][25:36, cx - 5:cx + 6] += rng.poisson(9000 * X._pix_gauss(0.0, 0.0, sigma, 5)).astype(float) if cx < 55 else 0.0
    est, n, q25, q75 = X.estimate_sigma0(fr, bg, 1.0, b0=0.0)
    assert abs(est - sigma) < 0.05 * sigma                                  # the blends fit a Gaussian badly (reduced chi^2) and are rejected, whatever their brightness
    edge = np.empty((20, 64, 64))
    for i in range(20):
        f = np.full((64, 64), 100.0)
        f[0:9, 0:9] += 6000 * X._pix_gauss(0.0, 0.0, 1.5, 5)[2:, 2:]           # a spot centred 3 px from the corner: closer to the border than the window
        edge[i] = rng.poisson(f).astype(float)
    est2, n2, _, _ = X.estimate_sigma0(edge, bg, 1.0, b0=0.0)
    assert n2 == 0 and np.isnan(est2)


def test_the_bic_rule_is_registered_and_names_its_stored_keys():
    assert X._rule_block(('bic', 4), dict(hist=[dict(w=np.zeros(0), chi2=1.0), dict(w=np.zeros(1), chi2=0.5)])) == 1
    arrays = {}
    X._store(arrays, 'dev', [('bic', 4)], [[np.zeros((0, 4))]])
    assert 'dev|bic4|off' in arrays and 'dev|bic4|flat' in arrays


def test_the_ols_variant_has_flat_weights_and_the_default_has_model_dependent_ones():
    assert X.VARIANTS['ols'] == dict(flat_var=True)
    P = X.IM.PixelPSF(64, 64, 1.4, parts=((1.0, 1.0),))
    saved = dict(X._W)
    try:
        X._W.update(Pp=P, a=2.0, bg=100.0)
        c, w = np.array([[30.2, 31.7, 1.4]]), np.array([5000.0])
        X._W['opts'] = dict()
        full0, full1 = X._sw(np.zeros((0, 3)), np.zeros(0)), X._sw(c, w)
        assert np.allclose(full0, 1.0 / np.sqrt(2.0 * 100.0 + X.B0))                 # the empty frame: the background's variance
        assert (full1 <= full0 + 1e-15).all() and full1.min() < 0.9 * full0.min()      # a source raises the variance of its pixels (the weights fall)
        X._W['opts'] = dict(flat_var=True)
        flat0, flat1 = X._sw(np.zeros((0, 3)), np.zeros(0)), X._sw(c, w)
        assert np.allclose(flat0, full0) and np.allclose(flat1, full0)                 # unweighted: the same map whatever the model
    finally:
        X._W.clear()
        X._W.update(saved)
