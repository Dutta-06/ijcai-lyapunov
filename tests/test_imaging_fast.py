"""The fast evaluation of the sliding path (src/core/imaging_fast.py) held to the original: the kernels, the lattice correlations, the trust-region iteration on the Gram matrix, a block fit, the
whole path (events, entry times, blocks), the limit-cycle guard (planted cycles through both engines), the exact early exit, the incremental model, and the fall-backs. The original functions
(`imaging.PixelPSF`, `Lattice`, `block_fit`, `refine`, `sliding_path`, `imaging_guard.sliding_path_guarded`) are the reference; they are not modified."""
import inspect
import os
import sys

import numpy as np
import pytest
import scipy
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.core import imaging as IM            # noqa: E402
from src.core import imaging_fast as F        # noqa: E402
from src.core import imaging_guard as G       # noqa: E402

PARTS = ((0.75, 1.0), (0.25, 3.0))
KW = dict(max_entries=8, tol=1e-8, passes=1, rtrunc=4.0)
HOLE = (30.0, 30.0)


def _psf(n=40, hp=4.0, parts=PARTS, sigma=1.3):
    return IM.PixelPSF(n, n, sigma, parts=parts, hp=hp)


def _world(seed=3, n=40, stars=((12.3, 11.8, 5000.0), (24.6, 14.2, 3000.0), (15.1, 27.4, 2000.0)), noise=0.5, hp=4.0, parts=PARTS, hole=None):
    psf = _psf(n, hp, parts)
    plain = IM.PixelPSF(n, n, psf.sigma, parts=parts)
    rng = np.random.default_rng(seed)
    y = plain.synth(np.array([[s[0], s[1]] for s in stars]), np.array([s[2] for s in stars]))
    if hole is not None:
        y = y + plain.synth(np.array([[hole[0], hole[1]]]), np.array([hole[2]]))
    y = y + noise * rng.standard_normal((n, n))
    return psf, psf.highpass(y), np.full((n, n), 1.0 / noise), IM.Lattice(psf, h=0.5, margin=1.0)


def _rel(a, b):
    return float(np.abs(np.asarray(a) - np.asarray(b)).max() / max(np.abs(b).max(), 1e-300))


# ---- the kernels ----------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('hp,parts', [(4.0, PARTS), (None, ((1.0, 1.0),)), (3.0, ((1.0, 1.0),)), (None, PARTS)], ids=['hp_halo', 'plain', 'hp_gauss', 'halo'])
@pytest.mark.parametrize('width', [False, True], ids=['fixed_sigma', 'own_sigma'])
def test_kernel_images_jacobian_and_correlations_equal_the_originals(hp, parts, width):
    psf = _psf(36, hp, parts)
    rng = np.random.default_rng(1)
    k = 5
    c = np.c_[rng.uniform(2, 33, k), rng.uniform(2, 33, k)]
    if width:
        c = np.c_[c, rng.uniform(1.0, 2.4, k)]
    w = rng.uniform(50, 500, k)
    sw = rng.uniform(0.2, 1.0, (36, 36))
    Rw = rng.standard_normal((36, 36))
    K = F.Kern(psf)
    assert _rel(K.phiT(c).T, psf.phi(c)) < 1e-13
    assert _rel(K.synth(c, w), psf.synth(c, w)) < 1e-13
    assert _rel(K.jacT(c, w, sw).T, psf.jac(c, w, sw)) < 1e-13
    for a, b in zip(K.corr(c[:, :2] if not width else c, Rw), psf.corr(c, Rw)):
        assert _rel(a, b) < 1e-13
    assert K.synth(c[:0], w[:0]).shape == (36, 36) and not K.synth(c[:0], w[:0]).any()


@pytest.mark.parametrize('hp,parts', [(4.0, PARTS), (None, ((1.0, 1.0),))], ids=['hp_halo', 'plain'])
def test_the_window_kernels_equal_the_originals(hp, parts):
    psf = _psf(40, hp, parts)
    pw = psf.window(5, 31, 8, 37)
    pw2 = pw.window(2, 20, 3, 25)                                        # a window of a window
    rng = np.random.default_rng(2)
    c = np.c_[rng.uniform(8, 36, 4), rng.uniform(5, 30, 4), rng.uniform(1.0, 2.0, 4)]
    w = rng.uniform(50, 500, 4)
    for q in (pw, pw2):
        sw = rng.uniform(0.2, 1.0, (q.ny, q.nx))
        K = F.Kern(q)
        assert _rel(K.phiT(c).T, q.phi(c)) < 1e-13
        assert _rel(K.synth(c, w), q.synth(c, w)) < 1e-13
        assert _rel(K.jacT(c, w, sw).T, q.jac(c, w, sw)) < 1e-13


def test_kernels_and_lattice_on_a_frame_that_is_not_square():
    """The two smoothers differ (different frame sides): x must use S_x and y S_y everywhere."""
    psf = IM.PixelPSF(30, 44, 1.3, parts=PARTS, hp=4.0)
    rng = np.random.default_rng(7)
    c = np.c_[rng.uniform(2, 41, 4), rng.uniform(2, 27, 4), rng.uniform(1.0, 2.2, 4)]
    w = rng.uniform(50, 500, 4)
    sw = rng.uniform(0.2, 1.0, (30, 44))
    Rw = rng.standard_normal((30, 44))
    K = F.Kern(psf)
    assert _rel(K.phiT(c).T, psf.phi(c)) < 1e-13 and _rel(K.synth(c, w), psf.synth(c, w)) < 1e-13 and _rel(K.jacT(c, w, sw).T, psf.jac(c, w, sw)) < 1e-13
    for a, b in zip(K.corr(c, Rw), psf.corr(c, Rw)):
        assert _rel(a, b) < 1e-13
    pw = psf.window(4, 26, 6, 40)
    Kw = F.Kern(pw)
    assert _rel(Kw.phiT(c).T, pw.phi(c)) < 1e-13 and _rel(Kw.jacT(c, w, sw[4:26, 6:40]).T, pw.jac(c, w, sw[4:26, 6:40])) < 1e-13
    L0, L1 = IM.Lattice(psf, h=0.5, margin=1.0), F.FastLattice(psf, h=0.5, margin=1.0)
    assert _rel(L1.corr(Rw), L0.corr(Rw)) < 1e-13


def test_pixel_edges_are_shared_between_neighbouring_pixels_exactly():
    """erf and exp are evaluated once per edge: the factors equal `PixelPSF.axis` to the last bit on every pixel."""
    psf = _psf(30, None, ((1.0, 1.0),))
    p = psf.px
    c = np.array([3.3, 14.7, 28.2])
    sg = np.array([1.2, 1.7, 2.3])
    a = F._axis(p, c, sg, True)
    b = psf.axis(p, c, sg, deriv_s=True)
    for u, v in zip(a, b):
        assert np.array_equal(u, v)


@pytest.mark.parametrize('hp', [4.0, None])
def test_fast_lattice_correlation_equals_the_original(hp):
    psf = _psf(40, hp, PARTS)
    L0, L1 = IM.Lattice(psf, h=0.5, margin=1.0), F.FastLattice(psf, h=0.5, margin=1.0)
    Rw = np.random.default_rng(0).standard_normal((40, 40))
    assert _rel(L1.corr(Rw), L0.corr(Rw)) < 1e-13
    assert L1.shape == L0.shape and np.array_equal(L1.xs, L0.xs)
    assert np.array_equal(L1.near(np.array([[10.0, 12.0]]), 2.0), L0.near(np.array([[10.0, 12.0]]), 2.0))


def test_fast_lattice_on_a_window_falls_back_to_the_original_correlation():
    psf = _psf(40, 4.0, PARTS).window(3, 30, 4, 35)
    L1 = F.FastLattice(psf, h=1.0, margin=1.0)
    Rw = np.random.default_rng(0).standard_normal((27, 31))
    assert np.array_equal(L1.corr(Rw), IM.Lattice(psf, h=1.0, margin=1.0).corr(Rw))


# ---- the trust-region iteration -------------------------------------------------------------------------------------------------------------

def test_trf_gram_matches_scipy_on_a_bounded_nonlinear_problem():
    t = np.linspace(0, 4, 60)
    rng = np.random.default_rng(0)
    y = 2.5 * np.exp(-1.3 * t) + 0.7 * np.exp(-0.2 * t) + 0.01 * rng.standard_normal(60)

    def fun(z):
        return z[0] * np.exp(-z[1] * t) + z[2] * np.exp(-z[3] * t) - y

    def jacT(z):
        return np.stack([np.exp(-z[1] * t), -z[0] * t * np.exp(-z[1] * t), np.exp(-z[3] * t), -z[2] * t * np.exp(-z[3] * t)])
    z0 = np.array([1.0, 0.4, 2.0, 0.9])
    lo, hi = np.array([0.0, 0.0, 0.0, 0.0]), np.array([10.0, 3.0, 10.0, 3.0])
    ref = least_squares(fun, z0, jac=lambda z: jacT(z).T, bounds=(lo, hi), method='trf', x_scale='jac', xtol=1e-12, ftol=1e-12, gtol=1e-12)
    x, nfev, status = F.trf_gram(fun, jacT, z0, lo, hi, 1e-12, 1e-12, 1e-12, 200)
    assert np.allclose(x, ref.x, rtol=1e-7, atol=1e-9) and abs(nfev - ref.nfev) <= 5 and status == ref.status


def test_trf_gram_respects_the_bounds_and_stops_at_a_bound_solution():
    t = np.linspace(0, 1, 30)
    y = 3.0 * t

    def fun(z):
        return z[0] * t - y

    def jacT(z):
        return t[None, :].copy()
    x, nfev, status = F.trf_gram(fun, jacT, np.array([0.5]), np.array([0.0]), np.array([2.0]), 1e-12, 1e-12, 1e-12, 100)
    assert 2.0 - 1e-6 < x[0] <= 2.0 and status != 0


def test_max_nfev_is_respected():
    t = np.linspace(0, 4, 60)
    y = 2.5 * np.exp(-1.3 * t)

    def fun(z):
        return z[0] * np.exp(-z[1] * t) - y

    def jacT(z):
        return np.stack([np.exp(-z[1] * t), -z[0] * t * np.exp(-z[1] * t)])
    x, nfev, status = F.trf_gram(fun, jacT, np.array([0.1, 3.0]), np.zeros(2), np.array([10.0, 10.0]), 1e-14, 1e-14, 1e-14, 3)
    assert nfev == 3 and status == 0


def test_rank_threshold_constants_and_the_scipy_functions_the_iteration_borrows():
    assert F.RANK_RATIO == 1e-6
    from scipy.optimize._lsq import trf as T
    assert list(inspect.signature(T.select_step).parameters) == ['x', 'J_h', 'diag_h', 'g_h', 'p', 'p_h', 'd', 'Delta', 'lb', 'ub', 'theta']
    assert tuple(int(v) for v in scipy.__version__.split('.')[:2]) >= (1, 10)


@pytest.mark.parametrize('width', [False, True], ids=['fixed_sigma', 'own_sigma'])
@pytest.mark.parametrize('hp,parts', [(4.0, PARTS), (None, ((1.0, 1.0),))], ids=['hp_halo', 'plain'])
def test_a_block_fit_agrees_with_scipy_to_the_iteration_tolerance(width, hp, parts):
    psf = _psf(36, hp, parts)
    plain = IM.PixelPSF(36, 36, psf.sigma, parts=parts)
    rng = np.random.default_rng(4)
    ctrue = np.array([[12.2, 13.1, 1.3], [17.4, 14.0, 1.3], [24.0, 24.5, 1.3]])
    wtrue = np.array([800.0, 500.0, 900.0])
    y = psf.highpass(plain.synth(ctrue[:, :2], wtrue) + 0.3 * rng.standard_normal((36, 36)))
    sw = np.full((36, 36), 1 / 0.3)
    c0 = ctrue + rng.uniform(-0.4, 0.4, ctrue.shape)
    c0 = c0 if width else c0[:, :2]
    w0 = np.array([300.0, 300.0, 300.0])
    lo, hi = c0 - 2.0, c0 + 2.0
    if width:
        lo[:, 2], hi[:, 2] = 1.0, 3.2
    a = IM.block_fit(psf, y, sw, c0, w0, lo, hi, max_nfev=200, tol=1e-10)
    b = F.fit_block_fast(psf, y, sw, c0, w0, lo, hi, max_nfev=200, tol=1e-10)
    assert _rel(b[0], a[0]) < 1e-6 and _rel(b[1], a[1]) < 1e-6
    assert abs(b[2]['nfev'] - a[2]['nfev']) <= 3
    j = F.fit_block_fast(psf, y, sw, c0, w0, lo, hi, joint=False)
    assert np.allclose(j[1], IM.block_fit(psf, y, sw, c0, w0, lo, hi, joint=False)[1], rtol=1e-9) and np.array_equal(j[0], c0)
    e = F.fit_block_fast(psf, y, sw, c0[:0], w0[:0], lo[:0], hi[:0])
    assert len(e[1]) == 0 and e[2]['nfev'] == 0


def test_a_failing_eigensolver_falls_back_to_scipys_block_fit_bit_for_bit(monkeypatch):
    psf = _psf(30, 4.0, PARTS)
    plain = IM.PixelPSF(30, 30, psf.sigma, parts=PARTS)
    y = psf.highpass(plain.synth(np.array([[12.2, 13.1]]), np.array([800.0])))
    sw = np.ones((30, 30))
    c0, w0 = np.array([[12.5, 13.0]]), np.array([300.0])

    def boom(*a, **k):
        raise np.linalg.LinAlgError('no convergence')
    monkeypatch.setattr(np.linalg, 'eigh', boom)
    a = IM.block_fit(psf, y, sw, c0, w0, c0 - 1, c0 + 1, tol=1e-8)
    b = F.fit_block_fast(psf, y, sw, c0, w0, c0 - 1, c0 + 1, tol=1e-8)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def test_refine_agrees_with_the_original():
    psf = _psf(36, 4.0, PARTS)
    plain = IM.PixelPSF(36, 36, psf.sigma, parts=PARTS)
    rng = np.random.default_rng(5)
    ct = np.array([[12.2, 13.1], [24.0, 24.5]])
    y = psf.highpass(plain.synth(ct, np.array([800.0, 900.0])) + 0.3 * rng.standard_normal((36, 36)))
    sw = np.full((36, 36), 1 / 0.3)

    def swf(c, w):
        return 1.0 / np.sqrt(0.09 + 0.01 * np.maximum(plain.synth(c, w) if len(w) else 0.0, 0.0))
    c0, w0 = ct + 0.2, np.array([600.0, 600.0])
    a = IM.refine(psf, y, sw, c0, w0, box=1.0, sw_fn=swf, sw_on='params')
    b = F.refine_fast(psf, y, sw, c0, w0, box=1.0, sw_fn=swf, sw_on='params')
    assert _rel(b[0], a[0]) < 1e-6 and _rel(b[1], a[1]) < 1e-6 and _rel(b[2], a[2]) < 1e-6
    e = F.refine_fast(psf, y, sw, c0[:0], w0[:0])
    assert len(e[1]) == 0


# ---- the path -------------------------------------------------------------------------------------------------------------------------------

def _agree(R0, R1, tol=1e-6, ttol=1e-8):
    assert [(e[1], e[2]) for e in R0['events']] == [(e[1], e[2]) for e in R1['events']]
    assert len(R0['hist']) == len(R1['hist'])
    assert np.allclose(R0['taus'], R1['taus'], rtol=ttol, atol=0) and np.allclose(R0['snr'], R1['snr'], rtol=1e-6)
    for a, b in zip(R0['hist'], R1['hist']):
        assert a['ids'] == b['ids'] and a['n_ent'] == b['n_ent']
        if len(a['w']):
            assert _rel(b['c'], a['c']) < tol and _rel(b['w'], a['w']) < tol
        assert abs(a['chi2'] - b['chi2']) <= 1e-8 * max(a['chi2'], 1e-12)


@pytest.mark.parametrize('hp,parts', [(4.0, PARTS), (None, ((1.0, 1.0),))], ids=['hp_halo', 'plain'])
@pytest.mark.parametrize('opts', [dict(), dict(r_ex=1.2), dict(fit_radius=5.0, r_ex=0.8), dict(local=False), dict(tau_max=None, order='greedy')],
                         ids=['default', 'r_ex', 'fit_radius', 'global', 'greedy'])
def test_the_fast_path_is_the_original_path(hp, parts, opts):
    psf, y, sw, lat = _world(hp=hp, parts=parts)
    R0 = IM.sliding_path(psf, y, sw, lat, **dict(KW, **opts))
    R1 = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **dict(KW, **opts))
    _agree(R0, R1)
    assert R1['barred'] == []


def test_the_fast_path_with_an_iterated_variance_and_sliding_widths():
    psf, y, sw, lat = _world()
    plain = IM.PixelPSF(40, 40, psf.sigma, parts=PARTS)

    def swf(c, w):
        return 1.0 / np.sqrt(0.25 + 0.002 * np.maximum(plain.synth(c, w) if len(w) else 0.0, 0.0))
    kw = dict(KW, sw_fn=swf, sw_on='params', sig_bounds=(1.0, 3.2), r_ex=1.0)
    R0 = IM.sliding_path(psf, y, swf(np.zeros((0, 3)), np.zeros(0)), lat, **kw)
    R1 = F.sliding_path_fast(psf, y, swf(np.zeros((0, 3)), np.zeros(0)), F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **kw)
    _agree(R0, R1, tol=1e-5)
    kw = dict(kw, sw_on='model', sw_fn=lambda m: np.full(m.shape, 2.0), passes=2)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    R1 = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **kw)
    _agree(R0, R1, tol=1e-5)


def test_the_incremental_model_equals_the_model_summed_afresh():
    psf, y, sw, lat = _world(stars=((12.3, 11.8, 5000.0), (24.6, 14.2, 3000.0), (15.1, 27.4, 2000.0), (30.0, 30.0, 1500.0), (8.0, 31.0, 1200.0)))
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    kw = dict(KW, max_entries=12, fit_radius=5.0, r_ex=0.8)
    A = F.sliding_path_fast(psf, y, sw, LF, incremental=True, refresh=0, **kw)
    B = F.sliding_path_fast(psf, y, sw, LF, incremental=False, **kw)
    C = F.sliding_path_fast(psf, y, sw, LF, incremental=True, refresh=3, **kw)
    _agree(B, A, tol=1e-8)
    _agree(B, C, tol=1e-8)
    assert any(len(h['w']) >= 3 for h in A['hist'])                      # the incremental branch (a fit of a strict subset of the units) was taken


def test_incremental_model_with_exits_equals_the_model_summed_afresh(monkeypatch):
    """A unit forced to leave (its coefficient 0) is subtracted from the carried model: same path as summing afresh."""
    psf, y, sw, lat, bad, start = _planted_world()
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    held = _plant_fast(monkeypatch, start, when={1})
    A = F.sliding_path_fast(psf, y, sw, LF, incremental=True, **dict(KW, fit_radius=4.0))
    held.clear()
    B = F.sliding_path_fast(psf, y, sw, LF, incremental=False, **dict(KW, fit_radius=4.0))
    assert any(e[2] == 'exit' for e in A['events'])
    _agree(B, A, tol=1e-8)


# ---- the limit-cycle guard through both engines -----------------------------------------------------------------------------------------------

def _planted_world(second=None, noise=0.5):
    stars = ((12.3, 11.8, 5000.0),) + ((second,) if second else ())
    psf, y, sw, lat = _world(hole=HOLE + (-3000.0,), stars=stars, noise=noise)
    bad = IM.sliding_path(psf, y, sw, lat, **dict(KW, max_entries=1))['events'][0][1]
    return psf, y, sw, lat, bad, lat.pos(*bad)


def _stub(real, starts, when):
    starts = np.atleast_2d(starts)
    held = []

    def stub(psf, y, sw, c, w, lo, hi, **kw):
        c0 = 0.5 * (lo + hi)
        j = np.flatnonzero(np.any(np.all(np.isclose(c0[:, None, :], starts[None, :, :]), axis=2), axis=1))
        if len(j):
            held.append(len(held) + 1)
            if when is None or held[-1] in when:
                c2, w2 = c.copy(), w.copy()
                w2[j], c2[j] = 0.0, np.array(HOLE)
                return c2, w2, dict(nfev=0, status=0)
        return real(psf, y, sw, c, w, lo, hi, **kw)
    return stub, held


def _plant_fast(monkeypatch, starts, when=None):
    stub, held = _stub(F.fit_block_fast, starts, when)
    monkeypatch.setattr(F, 'fit_block_fast', stub)
    return held


def _plant_both(monkeypatch, starts, when=None):
    s0, h0 = _stub(IM.block_fit, starts, when)
    s1, h1 = _stub(F.fit_block_fast, starts, when)
    monkeypatch.setattr(IM, 'block_fit', s0)
    monkeypatch.setattr(F, 'fit_block_fast', s1)
    return h0, h1


@pytest.mark.parametrize('r_ex', [0.0, 0.2, 1.2], ids=['no_zone', 'zone_0.2', 'zone_1.2'])
def test_the_planted_cycle_is_guarded_identically(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant_both(monkeypatch, start)
    kw = dict(KW, r_ex=r_ex)
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    Rf = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), **kw)
    assert [b for _, b in Rf['barred']] == [b for _, b in Rg['barred']] and len(Rg['barred']) == 1
    assert np.isclose(Rf['barred'][0][0], Rg['barred'][0][0], rtol=1e-12, atol=0)
    assert [(e[1], e[2]) for e in Rf['events']] == [(e[1], e[2]) for e in Rg['events']]
    assert np.allclose([e[0] for e in Rf['events']], [e[0] for e in Rg['events']], rtol=1e-12, atol=0)
    t0 = Rf['events'][0][0]
    assert [(e[1], e[2]) for e in Rf['events'][:4]] == [(bad, 'enter'), (bad, 'exit'), (bad, 'enter'), (bad, 'exit')] and len({e[0] for e in Rf['events'][:4]}) == 1
    assert all(e[1] != bad for e in Rf['events'][4:])
    assert len(Rf['taus']) == KW['max_entries'] and all(t > t0 for t in Rf['taus'][2:])
    assert np.allclose(Rf['taus'], Rg['taus'], rtol=1e-9)


@pytest.mark.parametrize('r_ex', [0.0, 1.2], ids=['no_zone', 'zone_1.2'])
def test_the_unguarded_fast_path_cycles_exactly_as_the_original(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant_both(monkeypatch, start)
    kw = dict(KW, r_ex=r_ex)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    R1 = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **kw)
    assert len(R1['taus']) == KW['max_entries'] and len(set(R1['taus'])) == 1 and R1['barred'] == []
    assert [(e[1], e[2]) for e in R1['events']] == [(e[1], e[2]) for e in R0['events']]
    assert np.allclose([e[0] for e in R1['events']], [e[0] for e in R0['events']], rtol=1e-12, atol=0)


def test_a_single_exit_followed_by_a_reentry_and_two_exits_at_different_times_are_not_barred(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    _plant_both(monkeypatch, start, when={1})
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **KW)
    Rf = F.sliding_path_fast(psf, y, sw, LF, **KW)
    assert Rf['barred'] == [] == Rg['barred'] and [(e[1], e[2]) for e in Rf['events']] == [(e[1], e[2]) for e in Rg['events']]
    assert sum(1 for e in Rf['events'] if e[2] == 'exit' and e[1] == bad) == 1
    monkeypatch.undo()
    _plant_both(monkeypatch, start, when={1, 3})
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **KW)
    Rf = F.sliding_path_fast(psf, y, sw, LF, **KW)
    ex = [e[0] for e in Rf['events'] if e[2] == 'exit' and e[1] == bad]
    assert len(ex) == 2 and ex[0] != ex[1] and Rf['barred'] == [] == Rg['barred']
    assert [(e[1], e[2]) for e in Rf['events']] == [(e[1], e[2]) for e in Rg['events']]


def test_the_guard_acts_on_a_second_planted_point_and_the_switch_defaults():
    p = inspect.signature(F.sliding_path_fast).parameters
    q = inspect.signature(G.sliding_path_guarded).parameters
    assert p['cycle_guard'].default is True
    for name, par in q.items():
        assert name in p and p[name].default == par.default                 # every option of the guarded path, the same default
    assert p['stop_kappa'].default is None and p['incremental'].default is True and p['refresh'].default == 64


# ---- the exact early exit -------------------------------------------------------------------------------------------------------------------

def _noisy_world():
    return _world(stars=((12.3, 11.8, 1200.0), (24.6, 14.2, 700.0), (15.1, 27.4, 250.0), (30.0, 30.0, 90.0)), noise=3.0)


@pytest.mark.parametrize('kappa', [3.0, 5.0, 8.0])
def test_the_early_exit_is_a_prefix_of_the_full_path_with_the_same_chosen_block(kappa):
    psf, y, sw0, lat = _noisy_world()
    sw = sw0 / 6.0
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    kw = dict(KW, max_entries=14, snr_min=3.0, snr_extra=3, r_ex=1.0)
    full = F.sliding_path_fast(psf, y, sw, LF, **kw)
    early = F.sliding_path_fast(psf, y, sw, LF, stop_kappa=kappa, **kw)
    ia, Ka = IM.stop_snr(full['hist'], full['taus'], full['snr'], kappa)
    ib, Kb = IM.stop_snr(early['hist'], early['taus'], early['snr'], kappa)
    assert Ka == Kb and len(early['taus']) <= len(full['taus']) and len(early['hist']) <= len(full['hist'])
    n = len(early['taus'])
    assert np.array_equal(early['taus'], full['taus'][:n]) and np.array_equal(early['snr'], full['snr'][:n])
    assert np.array_equal(early['hist'][ib]['c'], full['hist'][ia]['c']) and np.array_equal(early['hist'][ib]['w'], full['hist'][ia]['w'])
    assert early['hist'][ib]['n_ent'] == full['hist'][ia]['n_ent'] and early['hist'][ib]['ids'] == full['hist'][ia]['ids']
    if len(early['taus']) < len(full['taus']):                               # it did stop early, at the first entrant below kappa
        assert early['snr'][-1] < kappa and np.all(early['snr'][:-1] >= kappa)


def test_the_early_exit_actually_saves_fits_and_is_off_below_the_paths_own_threshold():
    psf, y, sw0, lat = _noisy_world()
    sw = sw0 / 6.0
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    kw = dict(KW, max_entries=14, snr_min=3.0, snr_extra=3, r_ex=1.0)
    full = F.sliding_path_fast(psf, y, sw, LF, **kw)
    early = F.sliding_path_fast(psf, y, sw, LF, stop_kappa=8.0, **kw)
    assert len(early['hist']) < len(full['hist'])
    off = F.sliding_path_fast(psf, y, sw, LF, stop_kappa=2.0, **kw)               # kappa below snr_min: ignored
    assert len(off['hist']) == len(full['hist']) and np.array_equal(off['taus'], full['taus'])
    none = F.sliding_path_fast(psf, y, sw, LF, stop_kappa=8.0, **dict(kw, snr_min=None))
    assert len(none['hist']) == len(F.sliding_path_fast(psf, y, sw, LF, **dict(kw, snr_min=None))['hist'])


def test_the_early_exit_of_the_fast_path_gives_the_original_paths_chosen_block():
    psf, y, sw0, lat = _noisy_world()
    sw = sw0 / 6.0
    kw = dict(KW, max_entries=14, snr_min=3.0, snr_extra=3, r_ex=1.0)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    R1 = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, stop_kappa=5.0, **kw)
    i0, K0 = IM.stop_snr(R0['hist'], R0['taus'], R0['snr'], 5.0)
    i1, K1 = IM.stop_snr(R1['hist'], R1['taus'], R1['snr'], 5.0)
    assert K0 == K1 and R0['hist'][i0]['ids'] == R1['hist'][i1]['ids']
    if K0:
        assert _rel(R1['hist'][i1]['c'], R0['hist'][i0]['c']) < 1e-6 and _rel(R1['hist'][i1]['w'], R0['hist'][i0]['w']) < 1e-6


def test_an_empty_frame_returns_an_empty_path():
    psf = _psf(24, 4.0, PARTS)
    lat = F.FastLattice(psf, h=1.0, margin=1.0)
    R = F.sliding_path_fast(psf, np.zeros((24, 24)), np.ones((24, 24)), lat, **KW)
    assert len(R['hist']) == 1 and len(R['taus']) == 0 and R['barred'] == []


# ---- the path's own stopping options ----------------------------------------------------------------------------------------------------------

def test_the_stopping_options_stop_where_the_original_stops():
    psf, y, sw0, lat = _noisy_world()
    sw = sw0 / 6.0
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    base = dict(KW, max_entries=14, r_ex=1.0)
    full = IM.sliding_path(psf, y, sw, lat, **base)
    ts = np.sort(full['taus'])
    tmid = float(0.5 * (ts[len(ts) // 2] + ts[len(ts) // 2 + 1]))                       # strictly between two entry times (a threshold equal to an entry time is decided by the last digit)
    cases = [dict(tau_max=tmid), dict(chi_stop=float(full['hist'][3]['chi2']) * 1.0001), dict(snr_min=5.0, snr_extra=1), dict(snr_min=3.0, snr_extra=2), dict(max_entries=3),
             dict(jump=1.5, extra=1)]
    for extra in cases:
        kw = dict(base, **extra)
        R0 = IM.sliding_path(psf, y, sw, lat, **kw)
        R1 = F.sliding_path_fast(psf, y, sw, LF, cycle_guard=False, **kw)
        _agree(R0, R1, tol=1e-5)
        assert len(R0['taus']) < len(full['taus']) or extra.get('snr_min') == 3.0 or 'jump' in extra       # the option did cut the path (the world's tail permitting)
    R0 = IM.sliding_path(psf, y, sw, lat, **dict(base, tau_max=tmid))
    assert len(R0['taus']) < len(full['taus'])


def test_the_early_exit_on_an_equal_threshold_saves_the_extra_entries():
    psf, y, sw0, lat = _noisy_world()
    sw = sw0 / 6.0
    LF = F.FastLattice(psf, h=0.5, margin=1.0)
    kw = dict(KW, max_entries=14, snr_min=5.0, snr_extra=3, r_ex=1.0)
    full = F.sliding_path_fast(psf, y, sw, LF, **kw)
    early = F.sliding_path_fast(psf, y, sw, LF, stop_kappa=5.0, **kw)          # stop_kappa == snr_min: the early exit is on
    assert early['snr'][-1] < 5.0 and len(early['taus']) < len(full['taus'])


# ---- added after the first mutation run: the iteration step for step, and the options the first tests left open --------------------------------------------

def _two_exp():
    t = np.linspace(0, 4, 60)
    y = 2.5 * np.exp(-1.3 * t) + 0.7 * np.exp(-0.2 * t) + 0.01 * np.random.default_rng(0).standard_normal(60)

    def fun(z):
        return z[0] * np.exp(-z[1] * t) + z[2] * np.exp(-z[3] * t) - y

    def jacT(z):
        return np.stack([np.exp(-z[1] * t), -z[0] * t * np.exp(-z[1] * t), np.exp(-z[3] * t), -z[2] * t * np.exp(-z[3] * t)])
    return fun, jacT, np.zeros(4), np.array([10.0, 3.0, 10.0, 3.0])


@pytest.mark.parametrize('z0', [[1, .4, 2, .9], [0, .4, 2, .9], [9.9, 2.9, .1, .1], [1, 2.95, 1, .05], [5, 1, .5, .1]], ids=['interior', 'zero_start', 'upper_corner', 'near_bound', 'far'])
@pytest.mark.parametrize('tol', [1e-8, 1e-10])
def test_the_gram_iteration_follows_scipys_iteration_step_for_step(z0, tol):
    fun, jacT, lo, hi = _two_exp()
    z0 = np.array(z0, float)
    ref = least_squares(fun, z0, jac=lambda z: jacT(z).T, bounds=(lo, hi), method='trf', x_scale='jac', xtol=tol, ftol=tol, gtol=tol)
    x, nfev, status = F.trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, 200)
    assert nfev == ref.nfev and status == ref.status
    assert np.allclose(x, ref.x, rtol=1e-12, atol=1e-13)


def test_a_variable_that_does_not_enter_the_residual_stays_and_the_rest_converges():
    fun0, jac0, lo, hi = _two_exp()

    def fun(z):
        return fun0(z[:4])

    def jacT(z):
        return np.vstack([jac0(z[:4]), np.zeros((1, 60))])                  # a zero column
    z0 = np.array([1.0, .4, 2.0, .9, 0.5])
    lo5, hi5 = np.r_[lo, 0.0], np.r_[hi, 1.0]
    ref = least_squares(fun, z0, jac=lambda z: jacT(z).T, bounds=(lo5, hi5), method='trf', x_scale='jac', xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=40)
    x, nfev, status = F.trf_gram(fun, jacT, z0, lo5, hi5, 1e-10, 1e-10, 1e-10, 40)
    c_mine, c_ref = 0.5 * np.sum(fun(x) ** 2), 0.5 * np.sum(fun(ref.x) ** 2)           # the rank-deficient steps are not scipy's (the null direction is arbitrary): the same cost after the same 40 evaluations (this valley is crept along by both solvers for hundreds of evaluations: the damped steps are not scipy's later)
    assert np.isfinite(x).all() and abs(c_mine - c_ref) <= 1e-6 * c_ref and abs(x[4] - 0.5) < 1e-6 and nfev == ref.nfev == 40


def test_the_default_evaluation_cap_is_a_hundred_per_variable():
    fun, jacT, lo, hi = _two_exp()
    x, nfev, status = F.trf_gram(fun, jacT, np.array([0.1, 3.0, 9.0, 0.1]), lo, hi, 1e-12, 1e-12, 1e-12, None)
    assert 4 < nfev <= 400 and status != 0


def test_a_non_finite_residual_shrinks_the_trust_region_and_the_iteration_goes_on():
    t = np.linspace(0, 1, 30)
    y = 3.0 * t

    def fun(z):
        return np.where(z[0] > 1.4, np.nan, z[0] * t - y)

    def jacT(z):
        return t[None, :].copy()
    x, nfev, status = F.trf_gram(fun, jacT, np.array([0.5]), np.array([0.0]), np.array([10.0]), 1e-10, 1e-10, 1e-10, 100)
    assert np.isfinite(x).all() and 1.3 < x[0] <= 1.4 and np.isfinite(fun(x)).all()


def test_a_start_outside_the_box_is_clipped_into_it_and_the_evaluation_cap_is_passed_on():
    psf = _psf(30, 4.0, PARTS)
    plain = IM.PixelPSF(30, 30, psf.sigma, parts=PARTS)
    y = psf.highpass(plain.synth(np.array([[12.2, 13.1]]), np.array([800.0])))
    sw = np.ones((30, 30))
    anchor = np.array([[12.0, 13.0]])
    c_out, w0 = np.array([[16.0, 9.0]]), np.array([300.0])                         # four pixels from the anchor, the box is 1 px
    lo, hi = anchor - 1.0, anchor + 1.0
    a = IM.block_fit(psf, y, sw, c_out, w0, lo, hi, tol=1e-9)
    b = F.fit_block_fast(psf, y, sw, c_out, w0, lo, hi, tol=1e-9)
    assert np.all(b[0] >= lo - 1e-12) and np.all(b[0] <= hi + 1e-12) and _rel(b[0], a[0]) < 1e-6 and _rel(b[1], a[1]) < 1e-5
    few = F.fit_block_fast(psf, y, sw, anchor + 0.3, w0, lo, hi, tol=1e-14, max_nfev=3)
    assert few[2]['nfev'] <= 3


def test_refine_keeps_the_centres_inside_its_box():
    psf = _psf(36, 4.0, PARTS)
    plain = IM.PixelPSF(36, 36, psf.sigma, parts=PARTS)
    y = psf.highpass(plain.synth(np.array([[12.2, 13.1]]), np.array([900.0])))
    sw = np.ones((36, 36))
    c0, w0 = np.array([[14.0, 15.0]]), np.array([600.0])                            # 1.8 px away: the optimum is outside the box
    a = IM.refine(psf, y, sw, c0, w0, box=1.0, rounds=1)               # one round: the box is about the start (a second round re-centres it)
    b = F.refine_fast(psf, y, sw, c0, w0, box=1.0, rounds=1)
    assert np.all(np.abs(b[0] - c0) <= 1.0 + 1e-9) and _rel(b[0], a[0]) < 1e-5 and _rel(b[1], a[1]) < 1e-4


def test_two_passes_on_a_cluster_use_the_updated_image_of_the_units_just_fitted():
    psf, y, sw, lat = _world(stars=((12.3, 11.8, 5000.0), (24.6, 14.2, 3000.0), (15.1, 27.4, 2000.0), (30.0, 30.0, 1500.0), (8.0, 31.0, 1200.0)))

    def swf(model):
        return np.full(model.shape, 2.0) / (1.0 + 1e-4 * np.maximum(model, 0.0))
    kw = dict(KW, max_entries=10, sw_fn=swf, sw_on='model', passes=2, fit_radius=5.0, r_ex=0.8)
    R0 = IM.sliding_path(psf, y, swf(np.zeros((40, 40))), lat, **kw)
    R1 = F.sliding_path_fast(psf, y, swf(np.zeros((40, 40))), F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **kw)
    _agree(R0, R1, tol=1e-5, ttol=1e-5)                                    # flat fits: the Gram iteration and the SVD iteration agree to the tolerance of the fits (2e-7 here), not to the last digit
    assert any(len(h['w']) >= 3 for h in R1['hist'])


def test_tied_entrants_are_all_seeds_of_the_local_refit():
    """Two mirror-image stars enter together (a tie of the duals); with a fit radius both neighbourhoods are refitted."""
    n = 40
    psf = _psf(n, 4.0, PARTS)
    plain = IM.PixelPSF(n, n, psf.sigma, parts=PARTS)
    y = psf.highpass(plain.synth(np.array([[10.0, 20.0], [29.0, 20.0]]), np.array([4000.0, 4000.0])))          # mirror images about x = 19.5
    sw = np.ones((n, n))
    lat = IM.Lattice(psf, h=0.5, margin=1.0)
    kw = dict(KW, max_entries=2, fit_radius=3.0)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    R1 = F.sliding_path_fast(psf, y, sw, F.FastLattice(psf, h=0.5, margin=1.0), cycle_guard=False, **kw)
    first = [e for e in R1['events'] if e[0] == R1['events'][0][0]]
    assert len(first) >= 2 and all(e[2] == 'enter' for e in first)          # the tie
    _agree(R0, R1, tol=1e-6)
    h = R1['hist'][1]
    assert len(h['w']) >= 2 and h['w'].min() > 100.0                        # both entrants were fitted in the first step


def test_an_unbounded_problem_started_at_zero_takes_scipys_unit_initial_radius():
    t = np.linspace(0, 1, 40)
    y = 3.0 * t + 2.0 + 0.01 * np.random.default_rng(3).standard_normal(40)

    def fun(z):
        return z[0] * t + z[1] - y

    def jacT(z):
        return np.stack([t, np.ones_like(t)])
    lo, hi = np.full(2, -np.inf), np.full(2, np.inf)
    ref = least_squares(fun, np.zeros(2), jac=lambda z: jacT(z).T, bounds=(lo, hi), method='trf', x_scale='jac', xtol=1e-10, ftol=1e-10, gtol=1e-10)
    x, nfev, status = F.trf_gram(fun, jacT, np.zeros(2), lo, hi, 1e-10, 1e-10, 1e-10, 100)
    assert nfev == ref.nfev and status == ref.status and np.allclose(x, ref.x, rtol=1e-9) and nfev > 2          # the first step is cut by the radius 1 (the solution is 3.6 away)


@pytest.mark.parametrize('z0', [[1, .4, 2, .9], [5, 1, .5, .1]])
def test_the_step_tolerance_stops_the_iteration_where_scipy_stops_it(z0):
    fun, jacT, lo, hi = _two_exp()
    z0 = np.array(z0, float)
    ref = least_squares(fun, z0, jac=lambda z: jacT(z).T, bounds=(lo, hi), method='trf', x_scale='jac', xtol=1e-2, ftol=1e-15, gtol=1e-15)
    x, nfev, status = F.trf_gram(fun, jacT, z0, lo, hi, 1e-15, 1e-2, 1e-15, 200)
    assert ref.status == 3 and status == 3 and nfev == ref.nfev and np.allclose(x, ref.x, rtol=1e-9)


def test_a_block_fit_takes_as_many_evaluations_as_scipys():
    psf = _psf(36, 4.0, PARTS)
    plain = IM.PixelPSF(36, 36, psf.sigma, parts=PARTS)
    rng = np.random.default_rng(4)
    ctrue = np.array([[12.2, 13.1], [17.4, 14.0], [24.0, 24.5]])
    y = psf.highpass(plain.synth(ctrue, np.array([800.0, 500.0, 900.0])) + 0.3 * rng.standard_normal((36, 36)))
    sw = np.full((36, 36), 1 / 0.3)
    c0 = ctrue + rng.uniform(-0.4, 0.4, ctrue.shape)
    a = IM.block_fit(psf, y, sw, c0, np.array([300.0, 300.0, 300.0]), c0 - 2.0, c0 + 2.0, tol=1e-10)
    b = F.fit_block_fast(psf, y, sw, c0, np.array([300.0, 300.0, 300.0]), c0 - 2.0, c0 + 2.0, tol=1e-10)
    assert a[2]['nfev'] == b[2]['nfev'] and a[2]['status'] == b[2]['status']
