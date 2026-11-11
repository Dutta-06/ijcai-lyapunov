"""The catalogue generator of the baselines given the ISS's front end (experiments/exp_astro_baselines_hp.py): the front end against independent computations, the wiring of the runners (checked with stand-in
packages so that the tests run without photutils and crowdsource), the guards and the test stage of the command line, and, where photutils and crowdsource are installed (the domain environment), planted
stars through the real packages."""
import hashlib
import os
import sys
import types

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import exp_astro_baselines as AB        # noqa: E402
import exp_astro_baselines_hp as G      # noqa: E402
import exp_astro_hp as E                # noqa: E402
import exp_astro_m2 as X                # noqa: E402
from src.core import imaging as IM      # noqa: E402

PHU, CS = G.FAMILIES['photutils'], G.FAMILIES['crowdsource']


def gauss_stamp(sigma, n=51):
    P = IM.PixelPSF(n, n, sigma)
    s = P.synth(np.array([[(n - 1) / 2, (n - 1) / 2]]), np.array([1.0]))
    return s / s.sum()


def planted_frame(sigma, centres, fluxes, sky=300.0, shape=(100, 100)):
    return sky + IM.PixelPSF(shape[0], shape[1], sigma).synth(np.atleast_2d(centres), np.atleast_1d(fluxes))


# ---- the registries and the plan ---------------------------------------------------------------------------------------------------------

def test_the_registries_have_eight_unique_variants_per_family_and_contain_the_references():
    assert len(PHU) == len(CS) == 8 and len(set(PHU)) == 8 and len(set(CS)) == 8
    assert G.PHOT_REF in PHU and set(G.CS_REFS) <= set(CS) and G.PHOT_REF == 'phu_hp5_new_i3' and G.CS_REFS == ('cs_default', 'cs_hp5_t5')
    assert not set(PHU) & set(CS)


def test_the_patch_origins_of_all_sets_are_distinct_and_there_are_eighteen():
    assert G.check_origins() == 18
    assert [len(G.SETS[k][1]) for k in ('dev', 'm2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3')] == [2, 4, 4, 2, 2, 2, 2]
    assert G.SETS['dev'][1] == [(110, 630), (210, 630)] and G.TEST_SETS == ('m2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3') and G.REPEAT_SETS == ('m2r1', 'm2r2', 'm2r3')


PLAN = os.path.join(HERE, '..', 'results', 'astro', 'baselines_hp_plan.md')


@pytest.mark.skipif(not os.path.exists(PLAN), reason='the plan is not in this copy')
def test_the_registries_are_the_variants_of_the_plan_tables():
    """The two tables of sections 4 and 5 of the plan (declared before the scoring) list exactly the variants and the parameters of the registries."""
    text = open(PLAN).read()
    s4 = text.split('## 4.')[1].split('## 5.')[0]
    s5 = text.split('## 5.')[1].split('## 6.')[0]
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in s4.split('\n') if l.startswith('| `phu_')]
    assert len(rows) == 8
    for r, d in zip(rows, G.PHOT_HP):
        assert r[0].split('`')[1] == d['tag'] and float(r[1]) == d['sb'] and r[2] == d['mode'] and int(r[3]) == d['iters']
        assert (r[4] == 'none') == (d['group'] == 0.0) and (r[4] == 'none' or ('min_separation=%g' % d['group']) in r[4])
        assert r[5].startswith('relaxed') == d['relax']
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in s5.split('\n') if l.startswith('| `cs_')]
    assert len(rows) == 8
    for r, prev in zip(rows, [None] + rows[:-1]):                               # 'same' in the plan's table repeats the cell above
        for j in (1, 2):
            if r[j] == 'same':
                r[j] = prev[j]
    for r, d in zip(rows, G.CS_HP):
        assert r[0].split('`')[1] == d['tag'] and float(r[-1]) == d['thr']
        if d['kind'] == 'hp':
            assert r[1].startswith('`H (dn - median)`') and ('%g px' % d['sb']) in r[1] and r[2].startswith('`H` of the `cs_default` stamp') and 'normalize=0' in r[2]
        else:
            assert r[1].startswith('`dn - sky_dn.mean()`') and not r[2].startswith('`H`')
        assert (r[3].startswith('`sw`')) == (d['kind'] != 'default')
    assert 'n_sig is fixed at 5' in s4 and G.NSIG == 5.0 and G.FIT == 7 and G.APER == 3.0 and G.GROUP == 4.4


# ---- the front end ----------------------------------------------------------------------------------------------------------------------

def test_iss_weights_are_the_inverse_noise_of_the_smoothed_counts():
    rng = np.random.default_rng(1)
    dn = 300.0 + 30.0 * rng.normal(size=(40, 60)) + 500.0 * (rng.uniform(size=(40, 60)) > 0.97)
    sw = G.iss_weights(dn, 4.6, 1.0)
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    assert sw.shape == dn.shape and np.allclose(sw, 1.0 / np.sqrt(expected / 4.6 + 1.0), rtol=1e-14, atol=0)
    assert not np.allclose(sw, 1.0 / np.sqrt(np.clip(dn, 1, None) / 4.6 + 1.0))                  # not the raw-count variance of the default baselines
    assert abs(G.noise_sigma(dn, 4.6, 1.0) - np.median(1.0 / sw)) < 1e-12
    flat = np.full((20, 20), 400.0)
    assert np.allclose(G.iss_weights(flat, 4.6, 1.0), 1 / np.sqrt(400 / 4.6 + 1.0))                # the dark variance is added to the counts / gain
    # the floor: 3 % of the pixels are dark and scattered, so the 1st percentile of the counts is the dark level and the 5th is the sky; the smoothed counts are above the former, below the latter
    rng = np.random.default_rng(3)
    dn2 = np.full((60, 60), 300.0)
    dn2[rng.uniform(size=dn2.shape) < 0.03] = 100.0
    p1 = np.percentile(dn2, 1)
    assert p1 == 100.0 and np.percentile(dn2, 5) == 300.0
    expected2 = np.maximum(gaussian_filter(dn2, 1.5, mode='reflect'), p1)
    assert (expected2 < 299.0).any() and np.allclose(G.iss_weights(dn2, 4.6, 1.0), 1.0 / np.sqrt(expected2 / 4.6 + 1.0), rtol=1e-14, atol=0)


def test_highpass_frame_is_the_isss_border_consistent_operator():
    rng = np.random.default_rng(2)
    dn = 300.0 + 40.0 * rng.normal(size=(60, 90)) + 2000.0 * np.exp(-0.5 * (((np.arange(90)[None, :] - 30) ** 2 + (np.arange(60)[:, None] - 20) ** 2) / 1.2 ** 2))
    for sb in (3.0, 5.0):
        y = dn - np.median(dn)
        Sy, Sx = IM.smoother(60, sb), IM.smoother(90, sb)
        h = G.highpass_frame(dn, sb)
        assert h.shape == dn.shape and np.array_equal(h, y - Sy @ y @ Sx.T)                            # rows are the first axis, columns the second; the very operations of the ISS
        assert np.allclose(h, IM.PixelPSF(60, 90, 0.9, hp=sb).highpass(y), atol=1e-10, rtol=0)
        # inside, 4 sigma_b from the border, it is the frame minus its Gaussian smoothing (the same truncation)
        m = int(np.ceil(4 * sb))
        ref = (y - gaussian_filter(y, sb, truncate=4.0, mode='constant'))[m:-m, m:-m]
        assert np.allclose(h[m:-m, m:-m], ref, atol=1e-9, rtol=0)
        assert np.allclose(G.highpass_frame(dn + 1234.5, sb), h, atol=1e-9, rtol=0)                      # constants and the median do not matter
        assert abs(G.highpass_frame(np.full((60, 90), 700.0), sb)).max() < 1e-9
    assert not np.allclose(G.highpass_frame(dn, 3.0), G.highpass_frame(dn, 5.0), atol=1e-3)


def test_highpass_stamp_is_the_interior_operator_of_a_zero_padded_stamp():
    stamp = gauss_stamp(0.93)
    for sb in (3.0, 5.0):
        pad = int(np.ceil(4 * sb))
        big = np.zeros((51 + 2 * pad,) * 2)
        big[pad:-pad, pad:-pad] = stamp
        ref = (big - gaussian_filter(big, sb, truncate=4.0, mode='constant'))[pad:-pad, pad:-pad]
        h = G.highpass_stamp(stamp, sb)
        assert h.shape == (51, 51) and np.allclose(h, ref, atol=1e-14, rtol=0)
        assert abs(h.sum()) < 1e-3 and h[25, 25] < stamp[25, 25] and h.min() < 0                              # a zero-sum bowl around a lower peak
        assert np.allclose(h, h.T, atol=1e-15, rtol=0)                                                      # a symmetric stamp stays symmetric (the operator is applied on both axes)
    asym = np.zeros((51, 51))
    asym[20:30, 22:27] = 1.0 / 50
    ha = G.highpass_stamp(asym, 3.0)
    assert np.allclose(ha, G.highpass_stamp(asym.T, 3.0).T, atol=1e-14) and not np.allclose(ha, ha.T)


def test_the_high_passed_stamp_is_the_model_of_a_high_passed_star():
    """A star of flux w on an integer pixel of a large frame, high-passed by the frame's operator, is w times the high-passed stamp (the frame and the stamp operators agree in the interior)."""
    stamp = gauss_stamp(0.93)
    shape = (140, 150)
    for sb in (3.0, 5.0):
        frame = np.full(shape, 350.0)
        cy, cx = 70, 75
        frame[cy - 25:cy + 26, cx - 25:cx + 26] += 30000.0 * stamp
        h = G.highpass_frame(frame, sb)
        got = h[cy - 25:cy + 26, cx - 25:cx + 26]
        assert np.allclose(got, 30000.0 * G.highpass_stamp(stamp, sb), atol=1e-6 * 30000, rtol=0)
        assert np.abs(got - 30000.0 * stamp).max() > 100.0                                              # and not the stamp itself


def test_cs_stamp_is_clipped_and_normalised_like_the_default_baseline():
    st = np.array([[0.2, -0.1], [0.0, 0.5]])
    out = G.cs_stamp(st)
    c = np.clip(st, 1e-10, np.inf)
    assert np.allclose(out, c / c.sum(), rtol=1e-15, atol=0) and abs(out.sum() - 1) < 1e-15 and out.min() > 0


def test_check_origins_refuses_two_patches_with_the_same_origin(monkeypatch):
    assert G.check_origins() == 18
    sets = dict(G.SETS)
    sets['m53'] = ('m53', [(310, 630)] + list(G.SETS['m53'][1][1:]))                                  # the origin of an M2 test patch
    monkeypatch.setattr(G, 'SETS', sets)
    with pytest.raises(SystemExit, match='share an origin'):
        G.check_origins()


def test_load_arrays_reads_the_patch_files_of_each_data_set(monkeypatch, tmp_path):
    def make(dirname, dn_value):
        d = tmp_path / dirname
        d.mkdir()
        np.savez_compressed(d / 'patch_arrays.npz', dn=np.full((100, 100), dn_value, dtype=np.float32), psf51=np.full((51, 51), 1.0 / 2601, dtype=np.float32), gain=np.float64(4.6),
                            dark_variance=np.float64(1.0), sky_dn=np.arange(4.0).reshape(2, 2))
        return str(d)
    dirs = {'m2': make('m2', 1.0), 'm53': make('m53', 2.0)}
    monkeypatch.setattr(AB, 'patch_dir', lambda x0, y0, n=100, ds='m2': dirs[ds])
    rep = make('rep', 3.0)
    monkeypatch.setattr(G.R, 'patch_path', lambda x0, y0, name: os.path.join(rep, 'patch_arrays.npz'))
    for ds, value in (('m2', 1.0), ('m53', 2.0), ('m2r1', 3.0)):
        a = G.load_arrays(ds, 5, 6)
        assert set(a) == {'dn', 'psf51', 'gain', 'dark', 'sky'} and a['dn'].dtype == np.float64 and a['psf51'].dtype == np.float64 and a['dn'][0, 0] == value
        assert a['gain'] == 4.6 and a['dark'] == 1.0 and a['sky'] == 1.5                               # the mean of sky_dn


def test_store_writes_the_keys_of_the_catalogue_files(monkeypatch):
    out = {}
    G.store(out, 110, 630, 'cs_w_t3', dict(x=np.arange(2.0), y=np.ones(2), flux=np.full(2, 5.0), dflux=np.zeros(2)), 7.5)
    assert sorted(out) == ['x110_y630|cs_w_t3|dflux', 'x110_y630|cs_w_t3|flux', 'x110_y630|cs_w_t3|secs', 'x110_y630|cs_w_t3|x', 'x110_y630|cs_w_t3|y'] and out['x110_y630|cs_w_t3|secs'].tolist() == [7.5]


# ---- the stand-in packages and the wiring of the runners ----------------------------------------------------------------------------------

class _Recorder:
    def __init__(self):
        self.finder = self.model = self.grouper = self.phot = self.call = None
        self.result = None


def _fake_photutils(monkeypatch, result):
    rec = _Recorder()
    rec.result = result

    class DAOStarFinder:
        def __init__(self, threshold, fwhm, **kw):
            self.threshold, self.fwhm, self.kw = threshold, fwhm, kw
            rec.finder = self

    class ImagePSF:
        def __init__(self, data, **kw):
            self.data, self.kw = data, kw
            rec.model = self

    class SourceGrouper:
        def __init__(self, min_separation):
            self.min_separation = min_separation
            rec.grouper = self

    class IterativePSFPhotometry:
        def __init__(self, model, fit_shape, finder, **kw):
            self.model, self.fit_shape, self.finder, self.kw = model, fit_shape, finder, kw
            rec.phot = self

        def __call__(self, data, error=None):
            rec.call = (data, error)
            return rec.result
    mods = {'photutils': types.ModuleType('photutils'), 'photutils.detection': types.ModuleType('photutils.detection'), 'photutils.psf': types.ModuleType('photutils.psf'),
            'astropy': types.ModuleType('astropy'), 'astropy.utils': types.ModuleType('astropy.utils'), 'astropy.utils.exceptions': types.ModuleType('astropy.utils.exceptions')}
    mods['photutils.detection'].DAOStarFinder = DAOStarFinder
    mods['photutils.psf'].ImagePSF, mods['photutils.psf'].SourceGrouper, mods['photutils.psf'].IterativePSFPhotometry = ImagePSF, SourceGrouper, IterativePSFPhotometry
    mods['astropy.utils.exceptions'].AstropyUserWarning = type('AstropyUserWarning', (UserWarning,), {})
    for name, m in mods.items():
        monkeypatch.setitem(sys.modules, name, m)
    return rec


def _patch(seed=5, shape=(100, 100)):
    rng = np.random.default_rng(seed)
    dn = 320.0 + 25.0 * rng.normal(size=shape) + 3000.0 * np.exp(-0.5 * (((np.arange(shape[1])[None, :] - 40) ** 2 + (np.arange(shape[0])[:, None] - 55) ** 2) / 1.1 ** 2))
    return dn, gauss_stamp(0.95), 4.6, 1.0


def test_run_photutils_hp_wires_the_front_end_the_finder_the_model_and_the_photometry(monkeypatch):
    dn, psf51, gain, dark = _patch()
    n = 5
    result = {'flux_fit': np.array([10.0, -3.0, np.nan, 7.0, 0.0]), 'x_fit': np.arange(n, dtype=float), 'y_fit': np.arange(n, dtype=float) * 10}
    for spec in G.PHOT_HP:
        rec = _fake_photutils(monkeypatch, result)
        cat, secs = G.run_photutils_hp(dn, psf51, gain, dark, spec)
        sb = spec['sb']
        err = 1.0 / G.iss_weights(dn, gain, dark)
        assert np.array_equal(rec.call[0], G.highpass_frame(dn, sb)) and np.array_equal(rec.call[1], err)                      # the high-passed data with the ISS's errors
        assert abs(rec.finder.threshold - 5.0 * np.median(err)) < 1e-12 and abs(rec.finder.fwhm - AB.psf_fwhm(psf51)) < 1e-15           # 5 noise sigmas, the stamp's FWHM
        assert rec.finder.kw == (dict(sharpness_range=(0.0, 2.0), roundness_range=(-2.0, 2.0)) if spec['relax'] else {})
        assert np.array_equal(rec.model.data, G.highpass_stamp(psf51, sb)) and rec.model.kw == dict(flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)
        assert rec.phot.fit_shape == (7, 7) and rec.phot.finder is rec.finder and rec.phot.model is rec.model
        assert rec.phot.kw == dict(grouper=rec.phot.kw['grouper'], maxiters=spec['iters'], mode=spec['mode'], aperture_radius=3.0, n_threads=1)
        assert (rec.phot.kw['grouper'] is None) == (spec['group'] == 0.0) and (spec['group'] == 0.0 or rec.phot.kw['grouper'].min_separation == spec['group'] == 4.4)
        assert np.array_equal(cat['x'], [0.0, 3.0]) and np.array_equal(cat['y'], [0.0, 30.0]) and np.array_equal(cat['flux'], [10.0, 7.0]) and secs >= 0       # finite and positive fluxes only
    rec = _fake_photutils(monkeypatch, None)                                                      # no source at all
    cat, secs = G.run_photutils_hp(dn, psf51, gain, dark, G.PHOT_HP[0])
    assert all(len(v) == 0 for v in cat.values()) and set(cat) == {'x', 'y', 'flux'}


def test_the_photutils_variants_differ_in_the_high_pass_the_mode_and_the_depth_as_declared(monkeypatch):
    by = {d['tag']: d for d in G.PHOT_HP}
    assert [(d['sb'], d['mode'], d['iters']) for d in G.PHOT_HP] == [(5.0, 'new', 3), (3.0, 'new', 3), (3.0, 'new', 10), (5.0, 'new', 10), (3.0, 'new', 30), (5.0, 'new', 30), (3.0, 'all', 6), (5.0, 'all', 6)]
    assert all(by[t]['group'] == 4.4 and by[t]['relax'] for t in ('phu_hp3_all_i6', 'phu_hp5_all_i6')) and not any(d['group'] or d['relax'] for d in G.PHOT_HP if d['mode'] == 'new')
    for d in G.PHOT_HP:
        assert d['tag'] == 'phu_hp%d_%s_i%d' % (d['sb'], d['mode'], d['iters'])


def _fake_crowdsource(monkeypatch, stars):
    rec = _Recorder()
    rec.fit_kw, rec.psfs, rec.images = None, [], []

    def central_stamp(stamp, censize=19):
        trim = (stamp.shape[-1] - censize) // 2
        return stamp[..., trim:stamp.shape[-1] - trim, trim:stamp.shape[-1] - trim]

    class SimplePSF:
        def __init__(self, stamp, normalize=19):
            self.stamp, self.normalize = np.array(stamp), normalize
            rec.psfs.append(self)

    def fit_im(images, psf, weights=None, dq=None, **kw):
        rec.fit_kw = dict(kw, images=images, psf=psf, weights=weights, dq=dq)
        return {'stars': stars}
    psfmod, base = types.ModuleType('crowdsource.psf'), types.ModuleType('crowdsource.crowdsource_base')
    psfmod.SimplePSF, psfmod.central_stamp = SimplePSF, central_stamp
    base.fit_im = fit_im
    top = types.ModuleType('crowdsource')
    top.psf, top.crowdsource_base = psfmod, base
    for name, m in (('crowdsource', top), ('crowdsource.psf', psfmod), ('crowdsource.crowdsource_base', base)):
        monkeypatch.setitem(sys.modules, name, m)
    return rec


def _stars():
    names = ['x', 'y', 'flux', 'dflux', 'qf', 'fracflux', 'spread_model', 'rchi2']
    st = np.zeros(3, dtype=[(n, 'f8') for n in names])
    st['x'], st['y'], st['flux'], st['dflux'], st['qf'], st['fracflux'], st['spread_model'], st['rchi2'] = [1.0, 2.0, 3.0], [10.0, 20.0, 30.0], [100.0, 200.0, 300.0], [1.0, 2.0, 3.0], [.1, .2, .3], [.4, .5, .6], [.7, .8, .9], [1.1, 1.2, 1.3]
    return st


def test_run_crowdsource_hp_wires_each_kind_of_variant(monkeypatch):
    dn, psf51, gain, dark = _patch()
    sky = 311.5
    sw = G.iss_weights(dn, gain, dark)
    base = G.cs_stamp(psf51)
    for spec in G.CS_HP:
        if spec['kind'] == 'default':
            continue
        rec = _fake_crowdsource(monkeypatch, _stars())
        cat, secs = G.run_crowdsource_hp(dn, sky, psf51, gain, dark, spec)
        kw = rec.fit_kw
        assert kw['threshold'] == spec['thr'] and np.array_equal(kw['weights'], sw) and kw['weights'].dtype == np.float64                       # the ISS's weights and the declared threshold
        assert kw['dq'].dtype == np.int64 and not kw['dq'].any() and kw['dq'].shape == dn.shape
        assert (kw['maxstars'], kw['fewstars'], kw['ntilex'], kw['ntiley'], kw['refit_psf'], kw['verbose']) == (10000000, 60, 1, 1, False, False)
        assert set(kw) == {'maxstars', 'fewstars', 'ntilex', 'ntiley', 'refit_psf', 'verbose', 'threshold', 'images', 'psf', 'weights', 'dq'}      # the package's other defaults stand
        assert len(rec.psfs) == 1 and kw['psf'] is rec.psfs[0]
        if spec['kind'] == 'raw':
            assert np.array_equal(kw['images'], dn - sky) and rec.psfs[0].normalize == 19 and np.array_equal(rec.psfs[0].stamp, base)
        else:
            s19 = base[16:35, 16:35].sum()
            assert np.array_equal(kw['images'], G.highpass_frame(dn, spec['sb'])) and rec.psfs[0].normalize == 0
            assert np.array_equal(rec.psfs[0].stamp, G.highpass_stamp(base / s19, spec['sb']))
            assert abs(rec.psfs[0].stamp.sum()) < 1e-3
        assert np.array_equal(cat['x'], [10.0, 20.0, 30.0]) and np.array_equal(cat['y'], [1.0, 2.0, 3.0])                                       # crowdsource's x runs along the first axis: swapped
        assert np.array_equal(cat['flux'], [100.0, 200.0, 300.0]) and np.array_equal(cat['dflux'], [1.0, 2.0, 3.0]) and np.array_equal(cat['qf'], [.1, .2, .3])
        assert np.array_equal(cat['fracflux'], [.4, .5, .6]) and np.array_equal(cat['spread'], [.7, .8, .9]) and np.array_equal(cat['rchi2'], [1.1, 1.2, 1.3]) and secs >= 0
    with pytest.raises(ValueError):
        G.run_crowdsource_hp(dn, sky, psf51, gain, dark, dict(tag='x', kind='other', sb=0.0, thr=5.0))


def test_the_default_crowdsource_variant_is_the_unchanged_baseline_call(monkeypatch):
    dn, psf51, gain, dark = _patch()
    seen = []

    def fake_run(*a):
        seen.append(a)
        print('centroid lines')
        sum(i * i for i in range(300000))                                                           # some CPU time to be measured
        return {'x': np.arange(3.0), 'y': np.zeros(3), 'flux': np.ones(3)}, 12.0
    monkeypatch.setattr(AB, 'run_crowdsource', fake_run)
    cat, secs = G.run_crowdsource_hp(dn, 300.0, psf51, gain, dark, G.CS_HP[0])
    assert len(seen) == 1 and seen[0][0] is dn and seen[0][1] == 300.0 and seen[0][2] is psf51 and seen[0][3:] == (gain, dark)
    assert np.array_equal(cat['x'], np.arange(3.0)) and 0.0 < secs < 12.0                              # the CPU seconds are measured here, not taken from the call
    assert G.CS_HP[0]['tag'] == 'cs_default' and G.CS_HP[0]['kind'] == 'default'


def test_run_variant_dispatches_on_the_tag_and_refuses_an_unknown_one(monkeypatch):
    calls = []
    monkeypatch.setattr(G, 'run_photutils_hp', lambda dn, psf, g, d, spec: (calls.append(('phu', spec['tag'])), ({}, 0.0))[1])
    monkeypatch.setattr(G, 'run_crowdsource_hp', lambda dn, sky, psf, g, d, spec: (calls.append(('cs', spec['tag'], sky)), ({}, 0.0))[1])
    a = dict(dn=0, psf51=1, gain=2, dark=3, sky=4)
    for tag in PHU + CS:
        G.run_variant(tag, a)
    assert calls == [('phu', t) for t in PHU] + [('cs', t, 4) for t in CS]
    with pytest.raises(KeyError):
        G.run_variant('nope', a)


# ---- the command line ---------------------------------------------------------------------------------------------------------------------

def _argv(monkeypatch, *args):
    monkeypatch.setattr(sys, 'argv', ['exp_astro_baselines_hp.py'] + list(args))


def test_the_command_line_refuses_what_the_plan_forbids(monkeypatch, tmp_path):
    out = str(tmp_path / 'c.npz')
    for args in (['--out', out], ['--stage', 'dev', '--out', out, '--variants', 'nope'], ['--stage', 'test', '--out', out]):
        _argv(monkeypatch, *args)
        with pytest.raises(SystemExit):
            G.main()
    (tmp_path / 'there.npz').write_bytes(b'x')
    _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'there.npz'))
    with pytest.raises(SystemExit, match='never overwritten'):
        G.main()


def _selection_raw(tmp_path, monkeypatch, sel_phu, sel_cs):
    plan = tmp_path / 'plan.md'
    plan.write_text('the plan\n')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    raw = {}
    for tags in (PHU, CS):
        for t in tags:
            for p in X.DEV:
                raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = np.array([100, 20, 100, 20, 0, 0, 0.0])
    for t, vals in ((sel_phu, (100, 60, 100, 70)), (sel_cs, (100, 50, 100, 55))):
        for p in X.DEV:
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = np.array(list(vals) + [0, 0, 0.0])
    raw['SEL|photutils'], raw['SEL|crowdsource'] = np.array([sel_phu]), np.array([sel_cs])
    raw['SEL|plan_sha256'] = np.array([hashlib.sha256(b'the plan\n').hexdigest()])
    raw['SEL|dev_hash'] = np.array([E.dev_hash(raw)])
    path = tmp_path / 'raw.npz'
    np.savez_compressed(path, **raw)
    return str(path), raw


def test_the_test_stage_runs_the_selected_variants_and_the_references_on_the_patches_of_the_sets(monkeypatch, tmp_path):
    rawp, raw = _selection_raw(tmp_path, monkeypatch, 'phu_hp3_all_i6', 'cs_hp3_t3')
    runs = []
    monkeypatch.setattr(G, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))

    def run(tag, a):
        runs.append((a['ds'], a['x0'], a['y0'], tag))
        return dict(x=np.array([a['x0'] + 0.5]), y=np.array([a['y0'] + 0.5]), flux=np.array([1000.0])), 2.5
    monkeypatch.setattr(G, 'run_variant', run)
    out = str(tmp_path / 'c.npz')
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm92', 'm53')
    G.main()
    tags = ['phu_hp3_all_i6', 'phu_hp5_new_i3', 'cs_hp3_t3', 'cs_default', 'cs_hp5_t5']
    assert runs == [(ds, p[0], p[1], t) for ds, plist in (('m92', G.SETS['m92'][1]), ('m53', G.SETS['m53'][1])) for p in plist for t in tags]
    z = np.load(out)
    assert len(z.files) == 6 * 5 * 4 and 'x1100_y220|cs_hp3_t3|secs' in z.files and z['x1100_y220|cs_hp3_t3|x'][0] == 1100.5 and z['x760_y590|cs_default|secs'][0] == 2.5
    assert not any('phu_hp3_new' in k or 'cs_w_' in k for k in z.files)                         # no other variant on a test patch
    # when the selection is a reference the tag is run once
    rawp2, _ = _selection_raw(tmp_path, monkeypatch, 'phu_hp5_new_i3', 'cs_hp5_t5')
    runs.clear()
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'c2.npz'), '--raw', rawp2, '--sets', 'm92', '--patches', '1430,550')
    G.main()
    assert runs == [('m92', 1430, 550, t) for t in ('phu_hp5_new_i3', 'cs_hp5_t5', 'cs_default')]


def test_the_test_stage_refuses_the_repeat_frames_without_the_flag_and_a_changed_selection(monkeypatch, tmp_path):
    rawp, raw = _selection_raw(tmp_path, monkeypatch, 'phu_hp3_new_i10', 'cs_hp3_t5')
    runs = []
    monkeypatch.setattr(G, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))
    monkeypatch.setattr(G, 'run_variant', lambda tag, a: (runs.append((a['ds'], tag)), (dict(x=np.zeros(1), y=np.zeros(1), flux=np.ones(1)), 0.1))[1])
    out = str(tmp_path / 'c.npz')
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm2r1')
    with pytest.raises(SystemExit, match='flag'):
        G.main()
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm2r1', '--flag', str(tmp_path / 'absent.flag'))
    with pytest.raises(SystemExit, match='flag'):
        G.main()
    flag = tmp_path / 'done.flag'
    flag.write_text('')
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm2r1', '--flag', str(flag))
    G.main()
    assert sorted({r[0] for r in runs}) == ['m2r1'] and len(runs) == 2 * 5                         # two patches, five distinct tags
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'd.npz'), '--raw', rawp, '--sets', 'm2r7')
    with pytest.raises(SystemExit, match='unknown set'):
        G.main()
    open(E.PLAN, 'w').write('edited\n')                                                         # the plan changed after the selection: nothing is run
    runs.clear()
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'e.npz'), '--raw', rawp, '--sets', 'm92')
    with pytest.raises(RuntimeError, match='plan file changed'):
        G.main()
    assert runs == []


def test_the_test_stage_without_sets_runs_the_m2_and_the_second_clusters_but_not_the_repeat_frames(monkeypatch, tmp_path):
    rawp, raw = _selection_raw(tmp_path, monkeypatch, 'phu_hp3_new_i10', 'cs_hp3_t5')
    runs = []
    monkeypatch.setattr(G, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))
    monkeypatch.setattr(G, 'run_variant', lambda tag, a: (runs.append((a['ds'], a['x0'], a['y0'])), (dict(x=np.zeros(1), y=np.zeros(1), flux=np.ones(1)), 0.1))[1])
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'c.npz'), '--raw', rawp)
    G.main()
    assert sorted({(r[0], r[1], r[2]) for r in runs}) == sorted((ds, p[0], p[1]) for ds in ('m2', 'm53', 'm92') for p in G.SETS[ds][1])


def test_the_command_line_stops_when_two_patches_share_an_origin(monkeypatch, tmp_path):
    def refuse():
        raise SystemExit('two patches of the sets share an origin')
    monkeypatch.setattr(G, 'check_origins', refuse)
    monkeypatch.setattr(G, 'run_variant', lambda tag, a: pytest.fail('a run was started'))
    _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'c.npz'))
    with pytest.raises(SystemExit, match='share an origin'):
        G.main()
    assert not (tmp_path / 'c.npz').exists()


def test_the_development_stage_runs_every_variant_on_the_two_development_patches_or_the_requested_subset(monkeypatch, tmp_path):
    runs = []
    monkeypatch.setattr(G, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))

    def run(tag, a):
        runs.append((a['x0'], a['y0'], tag))
        cat = dict(x=np.zeros(2), y=np.zeros(2), flux=np.ones(2))
        if tag.startswith('cs_'):
            cat.update(dflux=np.ones(2), qf=np.ones(2), fracflux=np.ones(2), spread=np.ones(2), rchi2=np.ones(2))
        return cat, 1.25
    monkeypatch.setattr(G, 'run_variant', run)
    out = str(tmp_path / 'dev.npz')
    _argv(monkeypatch, '--stage', 'dev', '--out', out)
    G.main()
    assert runs == [(x0, y0, t) for (x0, y0) in X.DEV for t in PHU + CS]
    z = np.load(out)
    assert 'x110_y630|cs_hp3_t3|dflux' in z.files and 'x210_y630|phu_hp5_all_i6|secs' in z.files and 'x210_y630|phu_hp5_all_i6|dflux' not in z.files
    assert len(z.files) == 2 * (8 * 4 + 8 * 9)
    runs.clear()
    _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'sub.npz'), '--variants', 'cs_w_t3', 'phu_hp3_new_i3', '--patches', '210,630')
    G.main()
    assert runs == [(210, 630, 'cs_w_t3'), (210, 630, 'phu_hp3_new_i3')]


def test_merge_is_the_union_of_catalogue_files_and_a_key_present_twice_must_be_equal(monkeypatch, tmp_path):
    np.savez_compressed(tmp_path / 'a.npz', a=np.arange(3.0), shared=np.ones(2))
    np.savez_compressed(tmp_path / 'b.npz', b=np.arange(2.0), shared=np.ones(2))
    _argv(monkeypatch, '--merge', str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz'), '--out', str(tmp_path / 'm.npz'))
    G.main()
    assert sorted(np.load(tmp_path / 'm.npz').files) == ['a', 'b', 'shared']
    np.savez_compressed(tmp_path / 'c.npz', shared=np.zeros(2))
    _argv(monkeypatch, '--merge', str(tmp_path / 'a.npz'), str(tmp_path / 'c.npz'), '--out', str(tmp_path / 'm2.npz'))
    with pytest.raises(ValueError):
        G.main()


# ---- planted stars through the real packages (the domain environment) -------------------------------------------------------------------------

try:
    import photutils      # noqa: F401
    import crowdsource    # noqa: F401
    HAVE_DOMAIN = True
except Exception:
    HAVE_DOMAIN = False
domain = pytest.mark.skipif(not HAVE_DOMAIN, reason='photutils and crowdsource are not installed (the domain environment)')


@domain
@pytest.mark.parametrize('sb', [3.0, 5.0])
def test_photutils_on_the_high_passed_patch_recovers_the_flux_and_the_position_of_a_planted_star(sb):
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    rng = np.random.default_rng(10 + int(sb))
    spec = dict(tag='t', sb=sb, mode='new', iters=3, group=0.0, relax=False)
    for _ in range(3):
        c = np.array([50.0, 50.0]) + rng.uniform(-0.5, 0.5, 2)
        dn = planted_frame(sigma, c, 20000.0)
        cat, _ = G.run_photutils_hp(dn, stamp, 4.6, 1.0, spec)
        k = np.argmax(cat['flux'])
        assert abs(cat['flux'][k] / 20000.0 - 1) < 0.015 and abs(cat['x'][k] - c[0]) < 0.02 and abs(cat['y'][k] - c[1]) < 0.02          # the planted-star check of the plan: +0.5 % and 0.004 px
        assert (cat['flux'] > 0.02 * 20000.0).sum() == 1                                                                              # no ghost from the bowl of the high-pass


@domain
def test_photutils_without_the_high_pass_of_the_stamp_would_be_biased():
    """The check has power: the high-passed data fitted with the plain stamp is off by an order of magnitude more."""
    from photutils.psf import PSFPhotometry, ImagePSF
    from astropy.table import Table
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    c = np.array([50.3, 49.8])
    y = G.highpass_frame(planted_frame(sigma, c, 20000.0), 3.0)
    res = {}
    for label, st in (('plain', stamp), ('highpassed', G.highpass_stamp(stamp, 3.0))):
        phot = PSFPhotometry(ImagePSF(st, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1), fit_shape=(7, 7), aperture_radius=3.0)
        r = phot(y, error=np.full(y.shape, 8.0), init_params=Table({'x_init': [50.0], 'y_init': [50.0]}))
        res[label] = abs(float(r['flux_fit'][0]) / 20000.0 - 1)
    assert res['highpassed'] < 0.015 and res['plain'] > 0.05


@domain
def test_the_all_mode_variant_runs_and_recovers_two_close_planted_stars():
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    c = np.array([[40.0, 50.0], [43.2, 50.6]])
    dn = planted_frame(sigma, c, [15000.0, 9000.0])
    spec = [d for d in G.PHOT_HP if d['tag'] == 'phu_hp3_all_i6'][0]
    cat, secs = G.run_photutils_hp(dn, stamp, 4.6, 1.0, dict(spec, iters=3))
    order = np.argsort(-cat['flux'])[:2]
    found = np.c_[cat['x'][order], cat['y'][order]]
    assert np.abs(found - c).max() < 0.1 and np.allclose(cat['flux'][order], [15000.0, 9000.0], rtol=0.03)


@domain
@pytest.mark.parametrize('spec', [d for d in G.CS_HP if d['kind'] != 'default' and d['thr'] == 5.0], ids=lambda d: d['tag'])
def test_crowdsource_recovers_the_flux_and_the_position_of_a_planted_star(spec):
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    rng = np.random.default_rng(21)
    for _ in range(2):
        c = np.array([50.0, 50.0]) + rng.uniform(-0.5, 0.5, 2)
        dn = planted_frame(sigma, c, 20000.0)
        cat, _ = G.run_crowdsource_hp(dn, 300.0, stamp, 4.6, 1.0, spec)
        assert len(cat['x']) == 1 and abs(cat['flux'][0] / 20000.0 - 1) < 0.01 and abs(cat['x'][0] - c[0]) < 0.02 and abs(cat['y'][0] - c[1]) < 0.02


@domain
def test_the_crowdsource_threshold_selects_the_faint_star_between_the_two_levels():
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    dn = planted_frame(sigma, np.array([[30.0, 30.0], [70.0, 70.0]]), [20000.0, 190.0])
    n = {}
    for thr in (3.0, 5.0, 8.0):
        spec = dict(tag='t', kind='raw', sb=0.0, thr=thr)
        n[thr] = len(G.run_crowdsource_hp(dn, 300.0, stamp, 4.6, 1.0, spec)[0]['x'])
    assert n[3.0] == n[5.0] == 2 and n[8.0] == 1                                                  # the 190 DN star has a significance of about 6.5


@domain
def test_the_raw_crowdsource_call_reproduces_the_baseline_call_when_given_its_weights():
    """`cs_fit` is the call of `exp_astro_baselines.run_crowdsource`: with the default weights it gives its catalogue exactly."""
    import crowdsource.psf as psfmod
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    dn = planted_frame(sigma, np.array([[30.2, 31.1], [60.7, 62.4], [61.9, 64.0]]), [20000.0, 6000.0, 3000.0])
    ref, _ = AB.run_crowdsource(dn, 300.0, stamp, 4.6, 1.0)
    sig = np.sqrt(np.clip(dn, 1, None) / 4.6 + 1.0)
    mine = G.cs_fit(dn - 300.0, psfmod.SimplePSF(G.cs_stamp(stamp)), 1.0 / sig, 5.0)
    assert sorted(mine) == sorted(ref) and all(np.array_equal(mine[k], ref[k]) for k in ref)
