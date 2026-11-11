"""The kernel study (experiments/exp_astro_psf.py): the path with the adequate kernel is the frozen path where the kernel is isotropic, the residual scale is a standard deviation, and the
report is read off the raw."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import astro_common as A                # noqa: E402
import exp_astro_m2 as X                # noqa: E402
import exp_astro_psf as E               # noqa: E402
from src.core import imaging as IM      # noqa: E402
from src.core import psf_general as G   # noqa: E402

ISO = ((0.53, 1.0), (0.35, 1.9), (0.12, 4.7))
LUMPY = ((0.55, 1.0, 1.15, 0.0, 0.0), (0.2, 0.9, 0.9, 0.7, -0.5), (0.15, 1.6, 1.3, -0.6, 0.4), (0.1, 3.2, 3.8, 0.1, 0.0))


def _patch(parts_fn, seed=3, noise=1.0, stars=((14.2, 17.6, 40000.0), (27.5, 12.1, 18000.0))):
    """A 40 x 40 patch: sky 150 DN, two stars rendered by `parts_fn` (general parts), Poisson noise (gain 4.6, dark variance 1) times `noise`."""
    rng = np.random.default_rng(seed)
    P = G.GeneralPixelPSF(40, 40, parts_fn)
    c = np.array([[s[0], s[1]] for s in stars])
    w = np.array([s[2] for s in stars])
    mean = 150.0 + P.synth(c, w)
    gain = 4.6
    dn = mean + noise * (rng.poisson(mean * gain) / gain - mean) + noise * rng.normal(0, 1.0 / gain, size=mean.shape)
    stamp = G.stamp_model(parts_fn, 51)
    return dict(dn=dn, psf=stamp / stamp.sum(), gain=gain, dark=1.0), c


def test_the_path_with_an_isotropic_kernel_is_the_frozen_path_bit_for_bit():
    patch, _ = _patch(G.isotropic_as_general(0.8, ISO))
    fit = G.fit_adequate(patch['psf'])
    assert fit['kind'] == 'isotropic'
    ks = (3.0, 5.0)
    a = X.run_path(patch, 3.0, 1.2, kappas=ks)
    b = E.run_path_psf(patch, 3.0, 1.2, fit, kappas=ks)
    assert sorted(a[0]) == sorted(b[0]) == [0.0, 3.0, 5.0] and a[2] == b[2] and a[3]['bic_index'] == b[3]['bic_index'] and a[3]['n_blocks'] == b[3]['n_blocks']
    for k in a[0]:
        for u, v in zip(a[0][k], b[0][k]):
            assert np.array_equal(u, v)


def test_the_path_with_the_general_kernel_finds_the_stars_and_reports_the_residual_scale():
    patch, c = _patch(LUMPY)
    fit = G.fit_adequate(patch['psf'])
    assert fit['kind'] == 'general'
    dets, secs, n, info = E.run_path_psf(patch, 3.0, 1.2, fit, kappas=(4.0,))
    x, y, f = dets[4.0]
    big = np.argsort(-f)[:2]
    found = np.c_[x[big], y[big]]
    found = found[np.argsort(found[:, 0])]
    assert np.abs(found - c[np.argsort(c[:, 0])]).max() < 0.15
    assert sorted(info['rsd']) == [0.0, 4.0] and all(0.7 < v < 1.5 for v in info['rsd'].values())


def test_residual_sd_is_a_standard_deviation_of_the_standardized_residual():
    flat = np.full((40, 40), 150.0)
    rng = np.random.default_rng(5)
    stamp = G.stamp_model(G.isotropic_as_general(0.8, ISO), 51)
    P = IM.PixelPSF(40, 40, 0.8, parts=ISO, hp=3.0)
    none = np.zeros((0, 3))
    sds = []
    for noise in (1.0, 2.0):
        dn = flat + noise * (rng.poisson(flat * 4.6) / 4.6 - flat)
        sds.append(E.residual_sd(dict(dn=dn, psf=stamp / stamp.sum(), gain=4.6, dark=1.0), P, none))
    assert 0.8 < sds[0] < 1.1 and abs(sds[1] / sds[0] - 2.0) < 0.15


def test_stamp_arrays_of_an_adequate_and_of_a_lumpy_stamp():
    a, fa = E.stamp_arrays(dict(psf=G.stamp_model(G.isotropic_as_general(0.8, ISO), 51)))
    b, fb = E.stamp_arrays(dict(psf=G.stamp_model(LUMPY, 51)))
    assert str(a['kind'][0]) == 'isotropic' and a['parts5'].shape == (0, 5) and a['rms'][0] == a['rms'][1]
    assert str(b['kind'][0]) == 'general' and b['parts5'].shape[1] == 5 and b['parts5'].shape[0] >= 3 and b['rms'][1] < b['rms'][0] and abs(b['parts5'][:, 0].sum() - 1.0) < 1e-12


def test_report_reads_the_raw_and_the_frozen_raw_and_is_deterministic(monkeypatch):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from test_astro_report import _synthetic_catalogs
    rng = np.random.default_rng(5)
    truth, cats = _synthetic_catalogs(rng)
    pxy = (1, 2)
    monkeypatch.setattr(X, 'SECOND', {'m53': dict(title='M53 (NGC 5024)', root='/nonexistent', prefix='p', truth='t.npz', patches=[pxy])})
    monkeypatch.setattr(E, 'CLUSTERS', ('m53',))
    monkeypatch.setattr(E, 'PATCHES', [('m53', pxy)])

    def sc(key):
        return A.score(np.c_[cats[key][0], cats[key][1]], cats[key][2], truth)
    froz, new = {}, {}
    pre_f, pre_n = 'T2|m53|x1_y2|', 'P|m53|x1_y2|'
    for name, key in (('ours|k2', 'mid'), ('default|k5', 'mid'), ('default|k0', 'mid')):
        froz.update({pre_f + kk: v for kk, v in X._pack(name, sc(key)).items()})
    for name, key in (('ours|k2', 'good'), ('default|k5', 'good'), ('default|k0', 'bic')):
        new.update({pre_n + kk: v for kk, v in X._pack(name, sc(key)).items()})
    for label, key in (('crowdsource (default)', 'poor'), ('crowdsource (tuned)', 'poor'), ('photutils (default)', 'poor'), ('photutils (tuned)', 'mid')):
        froz.update({pre_f + kk: v for kk, v in X._pack(label, sc(key)).items()})
    new[pre_n + 'cfg'] = np.array([3.0, 1.2, 2.0])
    new[pre_n + 'ours|rsd'] = np.array([1.84, 1.22])
    new[pre_n + 'default|rsd'] = np.array([1.5, 1.1])
    new['S|m53|x1_y2|kind'] = np.array(['general'])
    new['S|m53|x1_y2|rms'] = np.array([0.0073, 0.00098])
    new['S|m53|x1_y2|parts5'] = np.zeros((6, 5))
    new['S|m53|x1_y2|shape'] = np.array([0.112, -64.0, 1.37, 1.46])
    new['R|m2|x310_y630|ours|rsd'] = np.array([0.96])
    new['R|m2|x310_y630|default|rsd'] = np.array([1.01])
    new['S|m2|x310_y630|kind'] = np.array(['isotropic'])
    new['S|m2|x310_y630|rms'] = np.array([0.001, 0.001])
    new['S|m2|x310_y630|parts5'] = np.zeros((0, 5))
    monkeypatch.setattr(E, 'PATCHES', [('m2', (310, 630)), ('m53', pxy)])
    sizes = []
    for name in ('boot_idx', 'boot_idx_strat'):
        orig = getattr(X, name)
        monkeypatch.setattr(X, name, (lambda o, nm: (lambda *a, **k: (sizes.append((nm, a[0])), o(*a, **k))[1]))(orig, name))
    md = E.report(new, froz)
    assert ('boot_idx', 25) in sizes and ('boot_idx_strat', [25]) in sizes
    assert md == E.report(new, froz)
    tot = new[pre_n + 'ours|k2|blk22'].sum(0)
    row = [l for l in md.split('\n') if l.startswith('| ISS with the adequate kernel, the same configuration (kappa 2) | %d / %d |' % (tot[1], tot[0]))]
    assert len(row) >= 1                                                               # the pooled row of the adequate kernel is read from the new raw
    totf = froz[pre_f + 'ours|k2|blk22'].sum(0)
    assert any(l.startswith('| ISS as frozen (isotropic kernel), tuned on the M2 development patches (kappa 2) | %d / %d |' % (totf[1], totf[0])) for l in md.split('\n'))
    import re
    for name in ('adequate kernel vs frozen, tuned configuration', 'adequate kernel (tuned configuration) vs crowdsource (tuned)', 'adequate kernel (tuned configuration) vs photutils (tuned)'):
        assert re.search(r'\| %s \| [^\n]* \| win \|' % re.escape(name), md), name         # the better catalogue minus the worse one: a win, not a loss
    for needle in ('POST HOC, NOT A HELD-OUT TEST', '## The stamp fits (label-free)', '| m53 | (1, 2) | 0.00730 | general | 6 | 0.00098 | 0.112 | -64 |', '| m2 | (310, 630) | 0.00100 | isotropic | 3 | 0.00100 | - | - |',
                   '## M53 (NGC 5024): the 1 test patches', '### Patch (1, 2)', 'ISS with the adequate kernel, the same configuration (kappa 2)', 'ISS as frozen (isotropic kernel), tuned on the M2 development patches (kappa 2)',
                   'adequate kernel vs frozen, tuned configuration', 'adequate kernel (tuned configuration) vs crowdsource (tuned)', '| m53 | (1, 2) | 1.84 | 1.22 | 1.50 | 1.10 |', '| m2 | (310, 630) | 0.96 | - | 1.01 | - |', '## The residual scale of the fitted catalogues (label-free)',
                   '## Stars brighter than 19.5 mag: matched / in the truth', '| m53 | (1, 2) | '):
        assert needle in md, needle
    assert 'nan' not in md.lower().replace('nan,', '')


needs_data = pytest.mark.skipif(not all(os.path.exists(X.patch_path(p[0], p[1], ds)) for ds, p in (('m2', (110, 630)), ('m53', (760, 590)), ('m92', (1100, 220)))), reason='the astronomy patches are not on disk')


@needs_data
def test_the_adequacy_rule_on_the_real_stamps():
    """M2's psField stamps are fitted by the isotropic mixture to 0.1 % of the peak and keep it (the frozen results stand); M53's and M92's are fitted to 0.65 to 0.84 % and take the general kernel."""
    for ds, p, kind in (('m2', (110, 630), 'isotropic'), ('m2', (310, 630), 'isotropic'), ('m53', (760, 590), 'general'), ('m92', (1100, 220), 'general')):
        f = G.fit_adequate(X.load_patch(*p, ds=ds)['psf'])
        assert f['kind'] == kind, (ds, p)
        if kind == 'isotropic':
            assert f['rms'] < 1.2e-3
        else:
            assert f['rms_isotropic'] > 5e-3 and f['rms'] < 2e-3 and 4 <= len(f['parts5']) <= 8


def test_bright_counts_are_the_sums_of_the_nine_brightest_bins():
    raw = {'p|k|n_true': np.arange(17), 'p|k|n_matched': np.arange(17) // 2}
    assert E.bright(raw, 'p|', 'k') == (int((np.arange(9) // 2).sum()), int(np.arange(9).sum())) == (16, 36)
    assert E.bright(raw, 'q|', 'k') is None and E.NB == 9 and A.BIN_EDGES[E.NB] == 19.5


def test_stamp_shape_of_a_round_and_of_an_elliptical_stamp():
    r = E.stamp_shape(G.stamp_model(((1.0, 1.3, 1.3, 0.0, 0.0),), 41))
    assert r[0] < 1e-12 and abs(r[2] - r[3]) < 1e-12
    e = E.stamp_shape(G.stamp_model(((1.0, 1.0, 1.6, 0.0, 0.0),), 41))           # longer along y: orientation 90 degrees, sd_y > sd_x
    assert e[0] > 0.1 and abs(abs(e[1]) - 90.0) < 1e-6 and e[3] > e[2]
    f = E.stamp_shape(G.stamp_model(((1.0, 1.6, 1.0, 0.0, 0.0),), 41))           # longer along x: orientation 0
    assert abs(f[1]) < 1e-6 and f[2] > f[3]


def test_stamp_shape_values_against_an_independent_computation():
    n = 41
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    th = np.radians(30.0)
    u, v = (xx - c) * np.cos(th) + (yy - c) * np.sin(th), -(xx - c) * np.sin(th) + (yy - c) * np.cos(th)
    st = np.exp(-0.5 * (u / 1.6) ** 2 - 0.5 * (v / 1.0) ** 2)
    st = st / st.sum()
    r = E.stamp_shape(st)
    ell = (1.6 ** 2 - 1.0 ** 2) / (1.6 ** 2 + 1.0 ** 2)
    assert abs(r[0] - ell) < 0.02 and abs(r[1] - 30.0) < 1.0                       # tilted: both e1 and e2 are non-zero
    m = np.hypot(yy - c, xx - c) <= 5.0
    w = st * m / (st * m).sum()
    assert abs(r[2] - np.sqrt(np.average((xx - c) ** 2, weights=w))) < 1e-12 and abs(r[3] - np.sqrt(np.average((yy - c) ** 2, weights=w))) < 1e-12
    assert abs(E.stamp_shape(st, radius=3.0)[2] - r[2]) > 1e-3                      # the radius is a parameter


def test_part_run_runs_each_patch_once_checks_the_freeze_and_stores_under_the_part_key(monkeypatch, tmp_path):
    calls, frozen = [], []
    monkeypatch.setattr(X, 'SECOND', {'m53': dict(title='M53', patches=[(1, 2), (3, 4)]), 'm92': dict(title='M92', patches=[(5, 6)])})
    monkeypatch.setattr(E, 'CLUSTERS', ('m53', 'm92'))
    monkeypatch.setattr(X, 'check_freeze', lambda raw: frozen.append(1))
    monkeypatch.setattr(X, 'choose', lambda raw: (3.0, 1.2, 2.0, 0.56))
    m2raw = str(tmp_path / 'm2_raw.npz')
    np.savez_compressed(m2raw, dummy=np.zeros(1))

    def fake(args):
        ds, pxy, hp, rex, kap, path = args
        calls.append((ds, pxy, hp, rex, kap, path))
        return {'cfg': np.array([hp, rex, kap]), 'marker': np.array([pxy[0] * 10 + pxy[1]])}
    monkeypatch.setattr(E, '_psf_patch_arrays', fake)
    out = str(tmp_path / 'psf_raw.npz')
    E.part_run(out, workers=1, m2raw=m2raw)
    assert frozen == [1] and [(c[0], c[1]) for c in calls] == [('m53', (1, 2)), ('m53', (3, 4)), ('m92', (5, 6))] and all(c[2:5] == (3.0, 1.2, 2.0) and c[5] == m2raw for c in calls)
    raw = X.load_raw(out)
    assert raw['P|m53|x1_y2|marker'].tolist() == [12] and raw['P|m92|x5_y6|marker'].tolist() == [56] and raw['P|m53|x3_y4|cfg'].tolist() == [3.0, 1.2, 2.0]
    calls.clear()
    E.part_run(out, workers=1, m2raw=m2raw)                                           # nothing is run again
    assert calls == [] and len(frozen) == 2
