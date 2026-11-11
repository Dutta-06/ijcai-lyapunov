"""The general kernel family (src/core/psf_general.py): it reproduces the isotropic PixelPSF bit for bit where it should, its derivatives and high-pass partners are right, the fit of
a stamp is adequate or richer by the rule, and the sliding path runs on it."""
import os
import sys

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter
from scipy.special import erf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.core import imaging as IM            # noqa: E402
from src.core import psf_general as G         # noqa: E402

SIG, PARTS = 0.8, ((0.53, 1.0), (0.35, 1.9), (0.12, 4.7))
LUMPY = ((0.55, 1.0, 1.15, 0.0, 0.0), (0.2, 0.9, 0.9, 0.7, -0.5), (0.15, 1.6, 1.3, -0.6, 0.4), (0.1, 3.2, 3.8, 0.1, 0.0))


def _centres():
    return np.array([[7.3, 9.1], [15.0, 4.25], [3.9, 18.8]])


def test_isotropic_parts_reproduce_pixelpsf_bit_for_bit():
    for hp in (None, 3.0):
        P0 = IM.PixelPSF(26, 31, SIG, parts=PARTS, hp=hp)
        P1 = G.GeneralPixelPSF(26, 31, G.isotropic_as_general(SIG, PARTS), hp=hp)
        c = _centres()
        w = np.array([1.5, -0.7, 2.2])
        assert np.array_equal(P0.phi(c), P1.phi(c))
        assert np.array_equal(P0.synth(c, w), P1.synth(c, w))
        sw = np.full((26, 31), 0.37)
        assert np.array_equal(P0.jac(c, w, sw), P1.jac(c, w, sw))
        assert P0.reach() == P1.reach() and np.array_equal(P0.corr(c, sw)[0], P1.corr(c, sw)[0])
        L0, L1 = IM.Lattice(P0, h=0.5), IM.Lattice(P1, h=0.5)
        Rw = np.random.default_rng(1).normal(size=(26, 31))
        assert np.array_equal(L0.corr(Rw), L1.corr(Rw)) and np.array_equal(L0.norm2(sw ** 2), L1.norm2(sw ** 2))


def test_offset_part_is_the_shifted_gaussian():
    P1 = G.GeneralPixelPSF(20, 24, ((1.0, 1.3, 1.3, 0.7, -0.4),))
    P0 = IM.PixelPSF(20, 24, 1.3)
    c = np.array([[8.2, 9.9], [12.0, 5.5]])
    assert np.allclose(P1.phi(c), P0.phi(c + np.array([0.7, -0.4])), atol=1e-15, rtol=0)


def test_widths_by_axis_and_the_pixel_integral():
    kx, ky = 1.1, 1.9
    P = G.GeneralPixelPSF(15, 17, ((1.0, kx, ky, 0.0, 0.0),))
    phi = P.phi(np.array([[8.0, 7.0]]))[:, 0].reshape(15, 17)
    px, py = np.arange(17), np.arange(15)
    ix = 0.5 * (erf((px + 0.5 - 8.0) / (np.sqrt(2) * kx)) - erf((px - 0.5 - 8.0) / (np.sqrt(2) * kx)))
    iy = 0.5 * (erf((py + 0.5 - 7.0) / (np.sqrt(2) * ky)) - erf((py - 0.5 - 7.0) / (np.sqrt(2) * ky)))
    assert np.allclose(phi, np.outer(iy, ix), atol=1e-15, rtol=0)
    assert abs(phi.sum() - 1.0) < 2e-4                                    # the frame cuts the tails of the wider axis


@pytest.mark.parametrize('hp', [None, 3.0])
def test_jacobian_is_the_derivative_of_the_residual(hp):
    P = G.GeneralPixelPSF(22, 22, LUMPY, hp=hp)
    rng = np.random.default_rng(3)
    y = rng.normal(size=(22, 22))
    sw = rng.uniform(0.5, 1.5, size=(22, 22))
    c = np.array([[9.4, 11.2], [14.1, 6.3]])
    w = np.array([2.0, 1.3])

    def res(z):
        k = len(w)
        return (sw * (y - P.synth(np.c_[z[k:2 * k], z[2 * k:3 * k]], z[:k]))).ravel()
    z0 = np.r_[w, c[:, 0], c[:, 1]]
    J = P.jac(c, w, sw)
    eps = 1e-6
    for j in range(len(z0)):
        dz = np.zeros_like(z0)
        dz[j] = eps
        num = (res(z0 + dz) - res(z0 - dz)) / (2 * eps)
        assert np.allclose(J[:, j], num, atol=1e-7, rtol=1e-6), j


def test_highpass_partner_is_the_highpass_of_the_kernel():
    Pn = G.GeneralPixelPSF(24, 28, LUMPY)
    Ph = G.GeneralPixelPSF(24, 28, LUMPY, hp=3.0)
    c, w = _centres(), np.array([1.0, 2.0, -0.5])
    assert np.allclose(Ph.synth(c, w), Ph.highpass(Pn.synth(c, w)), atol=1e-12, rtol=0)


def test_window_is_the_cropped_kernel_of_the_full_frame():
    P = G.GeneralPixelPSF(30, 34, LUMPY, hp=3.0)
    c = np.array([[17.2, 14.6]])
    full = P.phi(c)[:, 0].reshape(30, 34)
    Wd = P.window(8, 21, 10, 26)
    assert isinstance(Wd, G.GeneralPixelPSF)
    assert np.allclose(Wd.phi(c)[:, 0].reshape(13, 16), full[8:21, 10:26], atol=1e-13, rtol=0)
    W2 = Wd.window(2, 9, 3, 12)                                          # a window of a window
    assert np.allclose(W2.phi(c)[:, 0].reshape(7, 9), full[10:17, 13:22], atol=1e-13, rtol=0)


def test_reach_grows_with_offsets_and_widths_and_includes_the_highpass():
    a = G.GeneralPixelPSF(10, 10, ((1.0, 1.0, 1.0, 0.0, 0.0),))
    b = G.GeneralPixelPSF(10, 10, ((1.0, 1.0, 1.0, 0.0, 1.5),))
    c = G.GeneralPixelPSF(10, 10, ((1.0, 1.0, 1.0, 0.0, 0.0),), hp=3.0)
    assert a.reach() == 1.0 and b.reach() == 2.5 and abs(c.reach() - np.hypot(1.0, 3.0)) < 1e-15


def test_base_methods_that_read_the_isotropic_parts_fail_loudly():
    P = G.GeneralPixelPSF(10, 10, LUMPY)
    assert P.parts is None
    with pytest.raises(NotImplementedError):
        P.factors(np.array([[3.0, 3.0]]))
    with pytest.raises(AssertionError):
        P.comps(np.array([[3.0, 3.0, 1.0]]))                             # a width of its own per kernel
    with pytest.raises(AssertionError):
        P.comps(np.array([[3.0, 3.0]]), deriv_s=True)                    # the derivative in that width
    with pytest.raises(AssertionError):
        P.reach(2.0)


def _stamp(parts, n=41):
    s = G.stamp_model(parts, n)
    return s / s.sum()


def test_stamp_model_is_centred_and_normalised():
    s = G.stamp_model(((0.6, 1.2, 1.2, 0.0, 0.0), (0.4, 2.5, 2.5, 0.0, 0.0)), 41)
    assert abs(s.sum() - 1.0) < 1e-6 and np.argmax(s) == 20 * 41 + 20 and np.allclose(s, s[::-1, ::-1], atol=1e-15)
    t = G.stamp_model(((1.0, 1.2, 1.2, 1.0, 0.0),), 41)
    assert np.unravel_index(np.argmax(t), t.shape) == (20, 21)


def test_fit_general_parts_recovers_a_lumpy_stamp():
    st = _stamp(LUMPY)
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    assert rms0 > G.RMS_GOAL                                            # the isotropic family cannot represent it
    rms, P = G.fit_general_parts(st, G.isotropic_as_general(sigma, parts))
    assert rms <= G.RMS_GOAL and abs(sum(p[0] for p in P) - 1.0) < 1e-12
    res = (G.stamp_model(P, st.shape[0]) - st) / st.max()
    assert abs(np.sqrt(np.mean(res ** 2)) - rms) < 1e-12              # the rms returned is the rms of the kernel returned
    assert np.abs(res).max() < 2e-2                                    # no pixel is off by more than 2 % of the peak
    assert rms < rms0


def test_parts_are_added_where_the_residual_is_largest_until_the_goal_or_no_gain():
    hard = ((0.3, 1.0, 1.0, 0.0, 0.0), (0.15, 0.7, 0.7, 1.1, 0.9), (0.12, 0.7, 0.7, -1.2, 0.8), (0.12, 0.7, 0.7, 0.9, -1.3), (0.1, 0.7, 0.7, -0.8, -1.1), (0.11, 1.8, 1.4, 0.0, 0.2),
            (0.1, 3.5, 3.5, 0.0, 0.0))
    st = _stamp(hard)
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    start = G.isotropic_as_general(sigma, parts)
    r3, P3 = G.fit_general_parts(st, start, k_max=3)                      # no part may be added
    r, P = G.fit_general_parts(st, start)
    assert len(P3) == 3 and len(P) > 3 and len(P) <= G.K_MAX and r < r3 < rms0
    r_loose, P_loose = G.fit_general_parts(st, start, rms_goal=10.0)     # already adequate by a loose goal: nothing is added
    assert len(P_loose) == 3 and r_loose == r3


def test_fit_adequate_keeps_the_isotropic_fit_where_it_is_adequate():
    st = _stamp(G.isotropic_as_general(SIG, PARTS))
    f = G.fit_adequate(st)
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    assert f['kind'] == 'isotropic' and f['parts5'] is None and f['sigma'] == sigma and f['parts'] == tuple(parts) and f['rms'] == f['rms_isotropic'] == rms0
    assert isinstance(G.make_psf(f, 12, 13, hp=3.0), IM.PixelPSF) and not isinstance(G.make_psf(f, 12, 13), G.GeneralPixelPSF)
    assert G.make_psf(f, 12, 13, hp=3.0).hp == 3.0 and G.make_psf(f, 12, 13).hp is None and G.make_psf(f, 12, 13, hp=0).hp is None


def test_fit_adequate_takes_the_general_family_where_it_is_not():
    f = G.fit_adequate(_stamp(LUMPY))
    assert f['kind'] == 'general' and f['parts'] is None and f['rms'] < f['rms_isotropic'] and f['rms_isotropic'] > G.RMS_GOAL
    P = G.make_psf(f, 12, 13, hp=3.0)
    assert isinstance(P, G.GeneralPixelPSF) and P.parts5 == f['parts5'] and (P.ny, P.nx, P.hp) == (12, 13, 3.0)


def test_sliding_path_finds_the_stars_with_the_general_kernel():
    rng = np.random.default_rng(7)
    ny = nx = 34
    Pn = G.GeneralPixelPSF(ny, nx, LUMPY)
    truth = np.array([[12.3, 15.8], [22.1, 10.4]])
    flux = np.array([9000.0, 5000.0])
    mean = 200.0 + Pn.synth(truth, flux)
    dn = rng.poisson(mean).astype(float)
    P = G.GeneralPixelPSF(ny, nx, LUMPY, hp=3.0)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    sw = 1.0 / np.sqrt(np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), 1.0))                  # the weights from the smoothed data, as the astronomy runs do
    y = P.highpass(dn - np.median(dn))
    scale = IM.noise_scale(L, sw ** 2)
    R = IM.sliding_path(P, y, sw, L, max_entries=12, tol=1e-6, passes=1, rtrunc=10.0 / P.reach(), r_ex=1.2, fit_radius=5.0, order='dual', tau_max=1.0 / (3.0 * scale))
    i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (6.0 * scale))
    h = R['hist'][i]
    big = np.argsort(-h['w'])[:2]
    found = h['c'][big][np.argsort(h['c'][big][:, 0])]
    want = truth[np.argsort(truth[:, 0])]
    assert K >= 2 and np.abs(found - want).max() < 0.15
    assert np.allclose(np.sort(h['w'][big]), np.sort(flux), rtol=0.06)


def test_the_number_of_parts_is_bounded_by_k_max_and_the_goal_is_a_parameter():
    hard = ((0.3, 1.0, 1.0, 0.0, 0.0), (0.15, 0.7, 0.7, 1.1, 0.9), (0.12, 0.7, 0.7, -1.2, 0.8), (0.12, 0.7, 0.7, 0.9, -1.3), (0.1, 0.7, 0.7, -0.8, -1.1), (0.11, 1.8, 1.4, 0.0, 0.2),
            (0.1, 3.5, 3.5, 0.0, 0.0))
    st = _stamp(hard)
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    start = G.isotropic_as_general(sigma, parts)
    for kmax in (4, 5):
        r, P = G.fit_general_parts(st, start, rms_goal=0.0, k_max=kmax)
        assert len(P) == kmax                                             # the goal 0 is never met: the parts are added up to k_max
    f = G.fit_adequate(_stamp(LUMPY), rms_goal=1e-9, k_max=4)
    assert f['kind'] == 'general' and len(f['parts5']) == 4
    f0 = G.fit_adequate(_stamp(LUMPY), rms_goal=1e-9, k_max=3)
    assert len(f0['parts5']) == 3 and f0['rms'] >= f['rms']
    assert G.fit_adequate(_stamp(LUMPY), rms_goal=1.0)['kind'] == 'isotropic'     # a loose goal accepts the isotropic fit
    assert (G.RMS_GOAL, G.K_MAX, G.W_LO, G.W_HI, G.D_MAX) == (2e-3, 10, 0.35, 12.0, 2.5)


def test_the_largest_residual_is_found_by_its_absolute_value():
    """A stamp with a deficit (a dent in the core) from the isotropic start: the part must go where the dent is, whatever the sign."""
    base = ((0.6, 1.0, 1.0, 0.0, 0.0), (0.4, 2.4, 2.4, 0.0, 0.0))
    st = _stamp(base)
    dent = G.stamp_model(((1.0, 0.8, 0.8, 1.0, 0.0),), 41)
    st = np.clip(st - 0.04 * dent, 1e-12, None)
    st = st / st.sum()
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    start = G.isotropic_as_general(sigma, parts)
    r, P = G.fit_general_parts(st, start, rms_goal=0.0, k_max=4)
    assert len(P) == 4 and r < G.fit_general_parts(st, start, rms_goal=0.0, k_max=3)[0] * 0.8


def test_window_without_a_highpass_uses_its_own_pixel_coordinates():
    P = G.GeneralPixelPSF(30, 34, LUMPY)
    c = np.array([[17.2, 14.6]])
    full = P.phi(c)[:, 0].reshape(30, 34)
    Wd = P.window(8, 21, 10, 26)
    assert (Wd.y0, Wd.x0) == (8, 10)
    assert np.allclose(Wd.phi(c)[:, 0].reshape(13, 16), full[8:21, 10:26], atol=1e-13, rtol=0)


def test_reach_is_the_largest_extent_over_the_axes_and_parts():
    P = G.GeneralPixelPSF(10, 10, ((0.5, 1.0, 1.8, 0.0, 0.0), (0.5, 2.1, 0.7, 0.0, 0.0)))
    assert P.reach() == 2.1
    P = G.GeneralPixelPSF(10, 10, ((0.5, 1.0, 1.8, 0.3, -1.0), (0.5, 0.5, 0.7, 0.0, 0.0)))
    assert abs(P.reach() - 2.8) < 1e-15


def test_fitted_parts_stay_inside_the_bounds():
    """A stamp that asks for a spike narrower than the bound and a satellite farther than the offset bound: the fit keeps every parameter inside."""
    wild = ((0.5, 1.2, 1.2, 0.0, 0.0), (0.2, 0.2, 0.2, 0.0, 0.0), (0.2, 0.6, 0.6, 3.6, -3.2), (0.1, 2.5, 2.5, 0.0, 0.0))
    st = _stamp(wild)
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    r, P = G.fit_general_parts(st, G.isotropic_as_general(sigma, parts), rms_goal=0.0, k_max=6)
    for a, kx, ky, dx, dy in P:
        assert 0.0 <= a <= 1.0 and G.W_LO - 1e-9 <= kx <= G.W_HI + 1e-9 and G.W_LO - 1e-9 <= ky <= G.W_HI + 1e-9 and abs(dx) <= G.D_MAX + 1e-9 and abs(dy) <= G.D_MAX + 1e-9
    assert max(max(abs(p[3]), abs(p[4])) for p in P) > 1.0              # a part did use the room


def test_a_bump_the_isotropic_start_lacks_is_found_where_the_residual_is_largest_in_absolute_value():
    base = ((0.5, 1.0, 1.0, 0.0, 0.0), (0.3, 2.4, 2.4, 0.0, 0.0), (0.2, 4.5, 4.5, 0.0, 0.0))
    st = _stamp(base) + 0.05 * G.stamp_model(((1.0, 0.7, 0.7, -1.6, 1.4),), 41)       # an excess: the model is below the stamp there; three parts cannot hold it, four can
    st = st / st.sum()
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    start = G.isotropic_as_general(sigma, parts)
    r3, _ = G.fit_general_parts(st, start, rms_goal=0.0, k_max=3)
    r4, P4 = G.fit_general_parts(st, start, rms_goal=0.0, k_max=4)
    assert r4 < 0.35 * r3 and len(P4) == 4


def test_no_part_is_added_when_none_lowers_the_residual():
    st = _stamp(G.isotropic_as_general(SIG, PARTS))
    sigma, parts, rms0 = IM.fit_psf_parts(st, n_parts=3, sigma0=1.0)
    r, P = G.fit_general_parts(st, G.isotropic_as_general(sigma, parts), rms_goal=0.0)
    assert len(P) == 3 and r < 1e-5
