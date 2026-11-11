"""Round 2 of the baselines given the ISS's front end, the catalogue generator (experiments/exp_astro_baselines_hp2.py): the registries against the round-2 plan and against the round-1 variants they extend, the wiring of
the runners (checked against the round-1 runners with stand-in packages, so that the tests run without photutils and crowdsource), the command line and the carrying of a selected round-1 variant, and, where photutils and
crowdsource are installed (the domain environment), planted stars through the real packages and the bit-for-bit reproduction of round 1."""
import hashlib
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import exp_astro_baselines_hp as G1     # noqa: E402
import exp_astro_baselines_hp2 as G2    # noqa: E402
import exp_astro_hp as E1               # noqa: E402
import exp_astro_hp2 as H2              # noqa: E402
import exp_astro_m2 as X                # noqa: E402
from test_exp_astro_baselines_hp import (HAVE_DOMAIN, _argv, _fake_crowdsource, _fake_photutils, _patch, _selection_raw, _stars, gauss_stamp, planted_frame)   # noqa: E402

PLAN2 = os.path.join(HERE, '..', 'results', 'astro', 'baselines_hp2_plan.md')
P1 = {d['tag']: d for d in G1.PHOT_HP}


# ---- the registries ----------------------------------------------------------------------------------------------------------------------

def test_the_photutils_variants_are_the_three_round_1_variants_with_xy_bounds_one_and_two():
    assert [d['tag'] for d in G2.PHOT2] == ['phu2_hp3_all_i6_b1', 'phu2_hp3_all_i6_b2', 'phu2_hp3_new_i30_b1', 'phu2_hp3_new_i30_b2', 'phu2_hp5_new_i3_b1', 'phu2_hp5_new_i3_b2']
    for d in G2.PHOT2:
        parent = P1[d['parent']]
        assert d['parent'] in ('phu_hp3_all_i6', 'phu_hp3_new_i30', 'phu_hp5_new_i3') and d['tag'] == 'phu2_%s_b%d' % (d['parent'][4:], d['bounds'])
        assert all(d[k] == parent[k] for k in ('sb', 'mode', 'iters', 'group', 'relax')) and d['bounds'] in (1.0, 2.0) and isinstance(d['bounds'], float)
    assert sorted({d['bounds'] for d in G2.PHOT2}) == [1.0, 2.0] and [d['mode'] for d in G2.PHOT2] == ['all', 'all', 'new', 'new', 'new', 'new']
    assert [d['iters'] for d in G2.PHOT2] == [6, 6, 30, 30, 3, 3] and [d['sb'] for d in G2.PHOT2] == [3.0, 3.0, 3.0, 3.0, 5.0, 5.0]


def test_the_crowdsource_variants_are_the_width_threshold_and_passes_grid():
    assert [d['tag'] for d in G2.CS2] == ['cs2_hp2_t3_pd', 'cs2_hp2_t3_pm', 'cs2_hp2_t5_pd', 'cs2_hp2_t5_pm', 'cs2_hp3_t3_pd', 'cs2_hp3_t3_pm', 'cs2_hp3_t5_pd', 'cs2_hp3_t5_pm']
    assert {(d['sb'], d['thr'], d['passes']) for d in G2.CS2} == {(sb, thr, ps) for sb in (2.0, 3.0) for thr in (3.0, 5.0) for ps in ('d', 'm')} and all(d['kind'] == 'hp' for d in G2.CS2)
    assert G2.PASSES == {'d': dict(miniter=4, maxiter=10, fewstars=60), 'm': dict(miniter=8, maxiter=20, fewstars=20)}


def test_the_families_have_the_declared_sizes_and_no_tag_of_round_1():
    assert len(G2.FAMILIES2['photutils']) == 6 and len(G2.FAMILIES2['crowdsource']) == 8
    all2 = G2.FAMILIES2['photutils'] + G2.FAMILIES2['crowdsource']
    assert len(set(all2)) == 14 and not set(all2) & set(G1.FAMILIES['photutils'] + G1.FAMILIES['crowdsource'])
    for t in G2.FAMILIES2['photutils']:                                            # the round-1 tags sort above the round-2 tags: a tie goes to round 1
        assert any(t < r for r in G1.FAMILIES['photutils'])
    assert all(t < 'cs_hp3_t5' for t in G2.FAMILIES2['crowdsource']) and all(t < 'phu_hp3_all_i6' for t in G2.FAMILIES2['photutils'])


@pytest.mark.skipif(not os.path.exists(PLAN2), reason='the round-2 plan is not in this copy')
def test_the_registries_are_the_variants_of_the_round_2_plan():
    text = open(PLAN2).read()
    s2 = text.split('## 2.')[1].split('## 3.')[0]
    s3 = text.split('## 3.')[1].split('## 4.')[0]
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in s2.split('\n') if l.startswith('| `phu2_')]
    assert len(rows) == 6
    for r, prev in zip(rows, [None] + rows[:-1]):                                           # 'same' repeats the cell above
        for j in range(1, 8):
            if r[j] == 'same':
                r[j] = prev[j]
    for r, d in zip(rows, G2.PHOT2):
        assert r[0].strip('`') == d['tag'] and r[1].split('`')[1] == d['parent'] and float(r[2]) == d['sb'] and r[3] == d['mode'] and int(r[4]) == d['iters'] and float(r[7]) == d['bounds']
        assert (r[5] == 'none') == (d['group'] == 0.0) and (r[6] == 'relaxed') == d['relax']
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in s3.split('\n') if l.startswith('| `cs2_')]
    assert len(rows) == 8
    for r, d in zip(rows, G2.CS2):
        assert r[0].strip('`') == d['tag'] and float(r[1]) == d['sb'] and float(r[2]) == d['thr'] and r[3] == d['passes']
    assert 'Passes `d`: `miniter 4, maxiter 10, fewstars 60`' in s3 and 'Passes `m`: `miniter 8, maxiter 20, fewstars 20`' in s3


# ---- the runners against the round-1 runners ------------------------------------------------------------------------------------------

def test_the_photutils_runner_is_the_round_1_runner_with_the_bound(monkeypatch):
    dn, psf51, gain, dark = _patch()
    result = {'flux_fit': np.array([10.0, -3.0, np.nan, 7.0, 0.0]), 'x_fit': np.arange(5, dtype=float), 'y_fit': np.arange(5, dtype=float) * 10}
    for spec2 in G2.PHOT2:
        rec1 = _fake_photutils(monkeypatch, result)
        cat1, _ = G1.run_photutils_hp(dn, psf51, gain, dark, P1[spec2['parent']])
        rec2 = _fake_photutils(monkeypatch, result)
        cat2, secs = G2.run_photutils_hp2(dn, psf51, gain, dark, spec2)
        assert np.array_equal(rec1.call[0], rec2.call[0]) and np.array_equal(rec1.call[1], rec2.call[1])
        assert rec1.finder.threshold == rec2.finder.threshold and rec1.finder.fwhm == rec2.finder.fwhm and rec1.finder.kw == rec2.finder.kw
        assert np.array_equal(rec1.model.data, rec2.model.data) and rec1.model.kw == rec2.model.kw and rec1.phot.fit_shape == rec2.phot.fit_shape == (7, 7)
        kw1, kw2 = dict(rec1.phot.kw), dict(rec2.phot.kw)
        g1, g2 = kw1.pop('grouper'), kw2.pop('grouper')
        assert (g1 is None) == (g2 is None) and (g1 is None or g1.min_separation == g2.min_separation)
        assert 'xy_bounds' not in kw1 and kw2.pop('xy_bounds') == spec2['bounds'] and kw1 == kw2
        assert all(np.array_equal(cat1[k], cat2[k]) for k in cat1) and secs >= 0
    rec = _fake_photutils(monkeypatch, None)
    cat, _ = G2.run_photutils_hp2(dn, psf51, gain, dark, G2.PHOT2[0])
    assert set(cat) == {'x', 'y', 'flux'} and all(len(v) == 0 for v in cat.values())


def test_the_crowdsource_runner_is_the_round_1_call_with_the_passes_and_the_width(monkeypatch):
    dn, psf51, gain, dark = _patch()
    import_base = G1.cs_stamp(psf51)
    for spec in G2.CS2:
        rec = _fake_crowdsource(monkeypatch, _stars())
        cat, secs = G2.run_crowdsource_hp2(dn, psf51, gain, dark, spec)
        kw = rec.fit_kw
        p = G2.PASSES[spec['passes']]
        assert kw['threshold'] == spec['thr'] and (kw['miniter'], kw['maxiter'], kw['fewstars']) == (p['miniter'], p['maxiter'], p['fewstars'])
        assert (kw['maxstars'], kw['ntilex'], kw['ntiley'], kw['refit_psf'], kw['verbose']) == (10000000, 1, 1, False, False)
        assert set(kw) == {'maxstars', 'fewstars', 'miniter', 'maxiter', 'ntilex', 'ntiley', 'refit_psf', 'verbose', 'threshold', 'images', 'psf', 'weights', 'dq'}
        assert np.array_equal(kw['images'], G1.highpass_frame(dn, spec['sb'])) and np.array_equal(kw['weights'], G1.iss_weights(dn, gain, dark)) and kw['weights'].dtype == np.float64
        assert kw['dq'].dtype == np.int64 and not kw['dq'].any()
        s19 = import_base[16:35, 16:35].sum()
        assert len(rec.psfs) == 1 and rec.psfs[0].normalize == 0 and np.array_equal(rec.psfs[0].stamp, G1.highpass_stamp(import_base / s19, spec['sb']))
        assert np.array_equal(cat['x'], [10.0, 20.0, 30.0]) and np.array_equal(cat['y'], [1.0, 2.0, 3.0]) and np.array_equal(cat['flux'], [100.0, 200.0, 300.0]) and secs >= 0
        assert set(cat) == {'x', 'y', 'flux', 'dflux', 'qf', 'fracflux', 'spread', 'rchi2'}
    # the cell (3 px, threshold 5, passes d) is the round-1 selected variant: the same call plus the package's own miniter and maxiter
    rec1 = _fake_crowdsource(monkeypatch, _stars())
    G1.run_crowdsource_hp(dn, 300.0, psf51, gain, dark, next(d for d in G1.CS_HP if d['tag'] == 'cs_hp3_t5'))
    rec2 = _fake_crowdsource(monkeypatch, _stars())
    G2.run_crowdsource_hp2(dn, psf51, gain, dark, next(d for d in G2.CS2 if d['tag'] == 'cs2_hp3_t5_pd'))
    k1, k2 = dict(rec1.fit_kw), dict(rec2.fit_kw)
    s1, s2 = k1.pop('psf').stamp, k2.pop('psf').stamp
    assert np.array_equal(s1, s2) and k2.pop('miniter') == 4 and k2.pop('maxiter') == 10
    assert set(k1) == set(k2) and all(np.array_equal(k1[k], k2[k]) for k in k1)


def test_run_variant2_dispatches_on_the_tag_and_refuses_an_unknown_one(monkeypatch):
    calls = []
    monkeypatch.setattr(G2, 'run_photutils_hp2', lambda dn, psf, g, d, spec: (calls.append(('phu', spec['tag'], dn, psf, g, d)), ({}, 0.0))[1])
    monkeypatch.setattr(G2, 'run_crowdsource_hp2', lambda dn, psf, g, d, spec: (calls.append(('cs', spec['tag'], dn, psf, g, d)), ({}, 0.0))[1])
    a = dict(dn=0, psf51=1, gain=2, dark=3, sky=4)
    for tag in G2.FAMILIES2['photutils'] + G2.FAMILIES2['crowdsource']:
        G2.run_variant2(tag, a)
    assert calls == [('phu', t, 0, 1, 2, 3) for t in G2.FAMILIES2['photutils']] + [('cs', t, 0, 1, 2, 3) for t in G2.FAMILIES2['crowdsource']]   # frame, stamp, gain, dark in this order
    with pytest.raises(KeyError):
        G2.run_variant2('phu_hp3_all_i6', a)                                       # a round-1 tag is not a round-2 variant


def test_the_test_stage_tags_are_listed_once(monkeypatch):
    monkeypatch.setattr(H2, 'verify_selection2', lambda raw: dict(photutils='t1', crowdsource='t2'))
    assert G2.test_stage_tags(None) == ['t1', 't2']
    monkeypatch.setattr(H2, 'verify_selection2', lambda raw: dict(photutils='t1', crowdsource='t1'))
    assert G2.test_stage_tags(None) == ['t1']


# ---- the command line --------------------------------------------------------------------------------------------------------------------

def _raw2(tmp_path, monkeypatch, sel_phu, sel_cs, sel1=('phu_hp3_all_i6', 'cs_hp3_t5')):
    """A round-2 raw with a selection: planted development totals (the selected candidates best), the round-1 raw and both plans planted in tmp_path, the module paths monkeypatched."""
    raw1_path, raw1 = _selection_raw(tmp_path, monkeypatch, *sel1)
    monkeypatch.setattr(H2, 'RAW1', raw1_path)
    plan2 = tmp_path / 'plan2.md'
    plan2.write_text('plan two\n')
    monkeypatch.setattr(H2, 'PLAN2', str(plan2))
    raw = {}
    for tags in (G2.FAMILIES2['photutils'] + [sel1[0]], G2.FAMILIES2['crowdsource'] + [sel1[1]]):
        for t in tags:
            for p in X.DEV:
                raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = np.array([100, 20, 100, 20, 0, 0, 0.0])
    for t, vals in ((sel_phu, (100, 60, 100, 70)), (sel_cs, (100, 50, 100, 55))):
        for p in X.DEV:
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = np.array(list(vals) + [0, 0, 0.0])
    raw['SEL|photutils'], raw['SEL|crowdsource'] = np.array([sel_phu]), np.array([sel_cs])
    raw['SEL|r1_photutils'], raw['SEL|r1_crowdsource'] = np.array([sel1[0]]), np.array([sel1[1]])
    raw['SEL|plan_sha256'] = np.array([hashlib.sha256(b'plan two\n').hexdigest()])
    raw['SEL|dev_hash'] = np.array([E1.dev_hash(raw)])
    path = tmp_path / 'raw2.npz'
    np.savez_compressed(path, **raw)
    return str(path), raw


def test_the_test_stage_runs_the_selected_round_2_variants_only_on_the_patches_of_the_sets(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, 'phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm')
    runs = []
    monkeypatch.setattr(G1, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))

    def run(tag, a):
        runs.append((a['ds'], a['x0'], a['y0'], tag))
        return dict(x=np.array([a['x0'] + 0.5]), y=np.array([a['y0'] + 0.5]), flux=np.array([1000.0])), 2.5
    monkeypatch.setattr(G2, 'run_variant2', run)
    out = str(tmp_path / 'c.npz')
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm92', 'm53')
    G2.main()
    assert runs == [(ds, p[0], p[1], t) for ds, plist in (('m92', G1.SETS['m92'][1]), ('m53', G1.SETS['m53'][1])) for p in plist for t in ('phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm')]
    z = np.load(out)
    assert len(z.files) == 6 * 2 * 4 and z['x1100_y220|cs2_hp2_t3_pm|x'][0] == 1100.5 and z['x760_y590|phu2_hp3_new_i30_b1|secs'][0] == 2.5


def test_the_test_stage_defaults_to_the_m2_m53_and_m92_test_patches_and_checks_the_origins(monkeypatch, tmp_path, capsys):
    rawp, raw = _raw2(tmp_path, monkeypatch, 'phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm')
    runs = []
    monkeypatch.setattr(G1, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))
    monkeypatch.setattr(G2, 'run_variant2', lambda tag, a: (runs.append((a['ds'], a['x0'], a['y0'], tag)), (dict(x=np.zeros(1), y=np.zeros(1), flux=np.ones(1)), 0.5))[1])
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'c.npz'), '--raw', rawp)
    G2.main()
    assert sorted(set(r[0] for r in runs)) == ['m2', 'm53', 'm92']
    assert runs == [(name, p[0], p[1], t) for name in ('m2', 'm53', 'm92') for p in G1.SETS[name][1] for t in ('phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm')]
    assert '%d patch origins, all distinct' % G1.check_origins() in capsys.readouterr().out


def test_a_selected_round_1_variant_is_carried_from_the_round_1_catalogues_and_not_run(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, 'phu_hp3_all_i6', 'cs2_hp3_t3_pd')
    runs = []
    monkeypatch.setattr(G1, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))
    monkeypatch.setattr(G2, 'run_variant2', lambda tag, a: (runs.append(tag), (dict(x=np.zeros(1), y=np.zeros(1), flux=np.ones(1)), 0.5))[1])
    c1 = {}
    for p in G1.SETS['m92'][1]:
        for f, v in (('x', np.arange(3.0) + p[0]), ('y', np.zeros(3)), ('flux', np.ones(3)), ('secs', np.array([9.5]))):
            c1['x%d_y%d|phu_hp3_all_i6|%s' % (p[0], p[1], f)] = v
        c1['x%d_y%d|phu_hp5_new_i3|x' % p] = np.zeros(1)                                  # another round-1 variant of the patch: not carried
    np.savez_compressed(tmp_path / 'cats1.npz', **c1)
    out = str(tmp_path / 'c.npz')
    _argv(monkeypatch, '--stage', 'test', '--out', out, '--raw', rawp, '--sets', 'm92', '--cats1', str(tmp_path / 'cats1.npz'))
    G2.main()
    assert runs == ['cs2_hp3_t3_pd'] * 2                                                 # the carried variant is not run
    z = np.load(out)
    assert sorted(z.files) == sorted(['x%d_y%d|phu_hp3_all_i6|%s' % (p[0], p[1], f) for p in G1.SETS['m92'][1] for f in ('x', 'y', 'flux', 'secs')] +
                                     ['x%d_y%d|cs2_hp3_t3_pd|%s' % (p[0], p[1], f) for p in G1.SETS['m92'][1] for f in ('x', 'y', 'flux', 'secs')])
    assert np.array_equal(z['x1100_y220|phu_hp3_all_i6|x'], np.arange(3.0) + 1100) and z['x1430_y550|phu_hp3_all_i6|secs'][0] == 9.5
    np.savez_compressed(tmp_path / 'cats0.npz', a=np.zeros(1))
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'd.npz'), '--raw', rawp, '--sets', 'm92', '--cats1', str(tmp_path / 'cats0.npz'))
    with pytest.raises(SystemExit, match='no arrays'):
        G2.main()
    # both selections are round-1 variants: nothing is run
    rawp2, _ = _raw2(tmp_path, monkeypatch, 'phu_hp3_all_i6', 'cs_hp3_t5')
    runs.clear()
    c1.update({'x%d_y%d|cs_hp3_t5|%s' % (p[0], p[1], f): np.zeros(2) for p in G1.SETS['m92'][1] for f in ('x', 'y', 'flux', 'secs')})
    np.savez_compressed(tmp_path / 'cats1.npz', **c1)
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'e.npz'), '--raw', rawp2, '--sets', 'm92', '--cats1', str(tmp_path / 'cats1.npz'))
    G2.main()
    assert runs == []


def test_the_test_stage_refuses_a_changed_plan_the_repeat_frames_without_the_flag_and_an_unknown_set(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, 'phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm')
    monkeypatch.setattr(G1, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))
    runs = []
    monkeypatch.setattr(G2, 'run_variant2', lambda tag, a: (runs.append(tag), (dict(x=np.zeros(1), y=np.zeros(1), flux=np.ones(1)), 0.5))[1])
    for extra, msg in ((['--sets', 'm2r1'], 'flag'), (['--sets', 'm2r1', '--flag', str(tmp_path / 'absent.flag')], 'flag'), (['--sets', 'm2r9'], 'unknown set')):
        _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'x.npz'), '--raw', rawp, *extra)
        with pytest.raises(SystemExit, match=msg):
            G2.main()
    flag = tmp_path / 'done.flag'
    flag.write_text('')
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'f.npz'), '--raw', rawp, '--sets', 'm2r1', '--flag', str(flag), '--patches', '316,124')
    G2.main()
    assert runs == ['phu2_hp3_new_i30_b1', 'cs2_hp2_t3_pm']
    open(H2.PLAN2, 'w').write('edited\n')
    runs.clear()
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'g.npz'), '--raw', rawp, '--sets', 'm92')
    with pytest.raises(RuntimeError, match='plan file changed'):
        G2.main()
    assert runs == []
    _argv(monkeypatch, '--stage', 'test', '--out', str(tmp_path / 'h.npz'))
    with pytest.raises(SystemExit, match='--raw'):
        G2.main()


def test_the_development_stage_runs_every_round_2_variant_on_the_development_patches_or_the_subset(monkeypatch, tmp_path):
    runs = []
    monkeypatch.setattr(G1, 'load_arrays', lambda ds, x0, y0: dict(ds=ds, x0=x0, y0=y0))

    def run(tag, a):
        runs.append((a['x0'], a['y0'], tag))
        cat = dict(x=np.zeros(2), y=np.zeros(2), flux=np.ones(2))
        if tag.startswith('cs2_'):
            cat.update(dflux=np.ones(2), qf=np.ones(2), fracflux=np.ones(2), spread=np.ones(2), rchi2=np.ones(2))
        return cat, 1.25
    monkeypatch.setattr(G2, 'run_variant2', run)
    out = str(tmp_path / 'dev.npz')
    _argv(monkeypatch, '--stage', 'dev', '--out', out)
    G2.main()
    assert runs == [(x0, y0, t) for (x0, y0) in X.DEV for t in G2.FAMILIES2['photutils'] + G2.FAMILIES2['crowdsource']]
    z = np.load(out)
    assert len(z.files) == 2 * (6 * 4 + 8 * 9) and 'x110_y630|cs2_hp2_t3_pm|dflux' in z.files and 'x210_y630|phu2_hp5_new_i3_b2|secs' in z.files
    runs.clear()
    _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'sub.npz'), '--variants', 'cs2_hp3_t5_pd', 'phu2_hp3_all_i6_b1', '--patches', '210,630')
    G2.main()
    assert runs == [(210, 630, 'cs2_hp3_t5_pd'), (210, 630, 'phu2_hp3_all_i6_b1')]
    for args in (['--variants', 'phu_hp3_all_i6'], ['--variants', 'nope']):                 # a round-1 tag is not a round-2 variant
        _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'u.npz'), *args)
        with pytest.raises(SystemExit, match='unknown'):
            G2.main()
    (tmp_path / 'there.npz').write_bytes(b'x')
    _argv(monkeypatch, '--stage', 'dev', '--out', str(tmp_path / 'there.npz'))
    with pytest.raises(SystemExit, match='never overwritten'):
        G2.main()
    _argv(monkeypatch, '--out', str(tmp_path / 'z.npz'))
    with pytest.raises(SystemExit):
        G2.main()


def test_merge_is_the_union_of_catalogue_files(monkeypatch, tmp_path):
    np.savez_compressed(tmp_path / 'a.npz', a=np.arange(3.0), shared=np.ones(2))
    np.savez_compressed(tmp_path / 'b.npz', b=np.arange(2.0), shared=np.ones(2))
    _argv(monkeypatch, '--merge', str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz'), '--out', str(tmp_path / 'm.npz'))
    G2.main()
    assert sorted(np.load(tmp_path / 'm.npz').files) == ['a', 'b', 'shared']


# ---- the real packages (the domain environment) -----------------------------------------------------------------------------------------

domain = pytest.mark.skipif(not HAVE_DOMAIN, reason='photutils and crowdsource are not installed (the domain environment)')


@domain
def test_the_round_2_crowdsource_runner_reproduces_the_round_1_variant_bit_for_bit_and_recovers_a_planted_star():
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    dn = planted_frame(sigma, np.array([[30.2, 31.1], [60.7, 62.4], [61.9, 64.0]]), [20000.0, 6000.0, 3000.0])
    r1, _ = G1.run_crowdsource_hp(dn, 300.0, stamp, 4.6, 1.0, next(d for d in G1.CS_HP if d['tag'] == 'cs_hp3_t5'))
    r2, _ = G2.run_crowdsource_hp2(dn, stamp, 4.6, 1.0, next(d for d in G2.CS2 if d['tag'] == 'cs2_hp3_t5_pd'))
    assert sorted(r1) == sorted(r2) and all(np.array_equal(r1[k], r2[k]) for k in r1)
    for spec in G2.CS2:
        c = np.array([50.2, 49.7])
        cat, _ = G2.run_crowdsource_hp2(planted_frame(sigma, c, 20000.0), stamp, 4.6, 1.0, spec)
        assert len(cat['x']) == 1 and abs(cat['flux'][0] / 20000.0 - 1) < 0.01 and abs(cat['x'][0] - c[0]) < 0.02 and abs(cat['y'][0] - c[1]) < 0.02, spec['tag']


@domain
def test_the_round_2_photutils_runner_without_a_bound_is_the_round_1_runner_and_with_a_bound_recovers_a_planted_star():
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    dn = planted_frame(sigma, np.array([[30.2, 31.1], [60.7, 62.4], [63.9, 64.0]]), [20000.0, 6000.0, 3000.0])
    for parent in ('phu_hp5_new_i3', 'phu_hp3_new_i30'):
        r1, _ = G1.run_photutils_hp(dn, stamp, 4.6, 1.0, P1[parent])
        r2, _ = G2.run_photutils_hp2(dn, stamp, 4.6, 1.0, dict(P1[parent], bounds=None))
        assert len(r1['x']) > 0 and all(np.array_equal(r1[k], r2[k]) for k in r1)
    c = np.array([50.2, 49.7])
    for spec in G2.PHOT2:
        cat, _ = G2.run_photutils_hp2(planted_frame(sigma, c, 20000.0), stamp, 4.6, 1.0, dict(spec, iters=min(spec['iters'], 3)))
        k = int(np.argmax(cat['flux']))
        assert abs(cat['flux'][k] / 20000.0 - 1) < 0.015 and abs(cat['x'][k] - c[0]) < 0.02 and abs(cat['y'][k] - c[1]) < 0.02, spec['tag']
        assert (cat['flux'] > 0.02 * 20000.0).sum() == 1


@domain
def test_a_tight_bound_keeps_the_fitted_positions_at_the_finder_positions():
    sigma = 0.93
    stamp = gauss_stamp(sigma)
    dn = planted_frame(sigma, np.array([[30.3, 31.2], [60.7, 62.4]]), [20000.0, 9000.0])
    spec = dict(P1['phu_hp5_new_i3'], tag='t', parent='phu_hp5_new_i3', iters=1)
    free, _ = G2.run_photutils_hp2(dn, stamp, 4.6, 1.0, dict(spec, bounds=None))
    tight, _ = G2.run_photutils_hp2(dn, stamp, 4.6, 1.0, dict(spec, bounds=1e-3))
    assert len(free['x']) == len(tight['x']) == 2
    moved = np.hypot(free['x'] - tight['x'], free['y'] - tight['y']).max()
    assert moved > 1e-3                                                                  # the bound changes the fit (the free fit leaves the finder's position by more than the bound)
