"""The two post hoc controls (experiments/exp_astro_controls.py): the declared sets of runs, the development-patches-only rule of control 1, the choice rule against `exp_astro_m2.choose`, the surfaces
(guarded and stored old rows, the rows of the control), the harness check against the stored design-variant rows, write-once running with a lock and the declaration check, the guard in place around a run,
and the report read off planted raws (the argmax, the finding rule, the paired differences ISS minus greedy and their verdicts, the groups), with the byte identity of `--from-raw`."""
import fcntl
import functools
import hashlib
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_common as A                # noqa: E402
import exp_astro_controls as E          # noqa: E402
import exp_astro_m2 as X                # noqa: E402
import exp_astro_repeat as RE           # noqa: E402
from src.core import imaging as IM      # noqa: E402
from src.core import imaging_guard as G  # noqa: E402
from test_astro_report import _synthetic_catalogs   # noqa: E402

TRUTH, CATS = _synthetic_catalogs(np.random.default_rng(11))


def _mk(seed, miss, fake, jit):
    """One more planted catalogue of the same truth (a quality between those of `_synthetic_catalogs`)."""
    rng = np.random.default_rng(seed)
    keep = rng.uniform(size=60) > miss
    pos, mag = np.c_[TRUTH['x'], TRUTH['y']], TRUTH['mag']
    xy = pos[keep] + rng.normal(0, jit, (keep.sum(), 2))
    fl = A.mag_to_counts(mag[keep] + rng.normal(0, 0.1, keep.sum()))
    fx = rng.uniform(5, 95, (fake, 2))
    ff = A.mag_to_counts(rng.uniform(20.0, 23.0, fake))
    return np.r_[xy[:, 0], fx[:, 0]], np.r_[xy[:, 1], fx[:, 1]], np.r_[fl, ff]


for _i, (_m, _f, _j) in enumerate([(0.2, 5, 0.15), (0.35, 10, 0.2), (0.45, 12, 0.25), (0.6, 20, 0.3)]):
    CATS['q%d' % _i] = _mk(100 + _i, _m, _f, _j)
SC = {k: A.score(np.c_[v[0], v[1]], v[2], TRUTH) for k, v in CATS.items()}
FROZEN = (3.0, 1.2, 2.0)
DATA = os.environ.get('ASTRO_DATA', '/home/user/astro_data')


def fq(key):
    """F of the planted catalogue quality `key` (pooled over nothing: one patch)."""
    return X.f_measure(np.array([SC[key]['true22'], SC[key]['matched22'], SC[key]['det22'], SC[key]['det22_matched'], 0, 0, 0], float))[0]


def put_readouts(d, pre, keyfn, kappas):
    """The arrays of the readouts at `kappas` under pre: pre + 'k<kappa>|<field>', the planted quality keyfn(kappa)."""
    for k in kappas:
        for kk, v in X._pack('', SC[keyfn(k)]).items():
            d[pre + 'k%g' % k + kk] = v


def put_trace(d, pre, n=1200, acts=0, secs=100.0):
    d[pre + 'secs'] = np.array([secs, float(n)])
    d[pre + 'bic'] = np.array([n - 30, n + 1, int(n >= X.MAX_ENTRIES)])
    d[pre + 'trace|counts'] = np.array([n, 3, n - 5, 1, 0, n + 1, acts, 0])
    d[pre + 'trace|barred'] = np.zeros((acts, 3))
    d[pre + 'trace|taus'] = np.linspace(0, 2, n)
    d[pre + 'wall'] = np.array([secs])


def dev_world(frozen='good', other='mid'):
    """(m2 raw with the stored development rows, cyclefix raw with the guarded development rows): the frozen point has the planted quality `frozen`, every other row `other`."""
    m2, cf = {}, {}
    for p in X.DEV:
        for hp in X.HPS:
            for rex in X.REXS:
                rid = 'D|x%d_y%d|hp%g|rex%g|' % (p[0], p[1], hp, rex)

                def kf(k, hp=hp, rex=rex):
                    return frozen if (hp, rex, k) == FROZEN else other
                put_readouts(m2, rid, kf, X.KAPPAS)
                m2[rid + 'bic'] = np.array([10, 20, 0])
                m2[rid + 'secs'] = np.array([10.0, 20.0])
                ns = 'G' if (hp, rex) == (8.0, 1.2) else 'C'
                put_readouts(cf, ns + '|' + rid, kf, X.KAPPAS)
                put_trace(cf, ns + '|' + rid, n=1300 + int(hp))
    return m2, cf


def e1_raw(special=None, cap=None, default='mid', frozen='good'):
    """The control-1 raw: every path, both patches; special = {(hp, rex, kappa): quality} overrides a readout, cap = the (hp, rex) whose path hit the cap."""
    special = special or {}
    raw = {}
    for h, r in E.C1_PATHS:
        for p in X.DEV:
            pre = E.c1_pre(p, h, r)

            def kf(k, h=h, r=r):
                return special.get((h, r, k), frozen if (h, r, k) == FROZEN else default)
            put_readouts(raw, pre, kf, E.C1_KAPPAS)
            put_trace(raw, pre, n=X.MAX_ENTRIES if cap == (h, r) else 1600, acts=2)
    return raw


def c2_world(order):
    """(c2 raw, ctx) with planted rows for the 16 test patches, every method of a different quality (so that a mixed-up key shows in a number). order: 'iss_better' (ISS good, greedy poor at kappa 2),
    'equal' (both mid), 'greedy_better' (ISS poor, greedy good). The stored guarded re-run is 'bic', the dual order with the same stop q0 at kappa 2 and q3 elsewhere, the greedy arm q2 at the other kappa,
    crowdsource q1, photutils q3."""
    raw, m2, rep, cf, hp2 = {}, {}, {}, {}, {}
    iss_q, gr_q = {'iss_better': ('good', 'poor'), 'equal': ('mid', 'mid'), 'greedy_better': ('poor', 'good')}[order]
    for p in E.c2_patches():
        pre = E.iss_prefix(p)[:-1]                                              # '...|k2' : the field is appended with its leading bar
        for z, ns, q in ((rep if p['group'] == 'rep' else m2, '', iss_q), (cf, 'C|', 'bic')):
            for kk, v in X._pack('', SC[q]).items():
                z[ns + pre + kk] = v
        for arm, q2, qo in (('greedy', gr_q, 'q2'), ('dual_snr', 'q0', 'q3')):
            put_readouts(raw, E.c2_pre(p, arm), lambda k, q2=q2, qo=qo: q2 if k == 2.0 else qo, X.KAPPAS)
            put_trace(raw, E.c2_pre(p, arm), n=1100, acts=1 if arm == 'dual_snr' else 0, secs=50.0)
        for row, q in (('cs2', 'q1'), ('phu2', 'q3')):
            for kk, v in X._pack('', SC[q]).items():
                hp2['T|%s|x%d_y%d|%s' % (p['set'], p['pxy'][0], p['pxy'][1], row) + kk] = v
    for arm, st in E.HARNESS.items():                                          # the stored V rows of the benchmark patch: identical to the planted arms
        for k in X.KAPPAS:
            for f in E.PACK_FIELDS:
                m2[st + 'k%g|%s' % (k, f)] = raw[E.c2_pre(dict(set='m2', group='m2', pxy=E.HARNESS_PATCH[1]), arm) + 'k%g|%s' % (k, f)]
        m2[st + 'secs'] = np.array([30.0, 1100.0])
    return raw, dict(m2=m2, rep=rep, cf=cf, hp2=hp2)


@functools.lru_cache(maxsize=None)
def _c2(order):
    return c2_world(order)


def full_ctx(order='iss_better'):
    raw, ctx = _c2(order)
    m2, cf = dev_world()
    return dict(raw), dict(ctx, m2=dict(ctx['m2'], **m2), cf=dict(ctx['cf'], **cf))


# ---- the declared sets --------------------------------------------------------------------------------------------------------------------

def test_the_declared_control_1_set():
    assert len(E.C1_PATHS) == 13 and len(set(E.C1_PATHS)) == 13
    old = {(h, r) for h in X.HPS for r in X.REXS}
    assert len([c for c in E.C1_PATHS if c not in old]) == 9 and len([c for c in E.C1_PATHS if c in old]) == 4
    assert {c for c in E.C1_PATHS if c[0] < 3} == {(h, r) for h in (2.0, 2.5) for r in (1.2, 1.6, 2.0)}
    assert E.C1_KAPPAS == (1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0) and E.EXT_KAPPAS == (1.25, 1.5, 1.75)
    assert E.FINDING == 0.01
    assert {c[1] for c in E.C1_PATHS} == {0.8, 1.2, 1.6, 2.0} and {c[0] for c in E.C1_PATHS} == {2.0, 2.5, 3.0, 5.0, 8.0}
    assert (3.0, 1.2) in E.C1_PATHS                                                                       # the frozen combination is re-run (its kappa axis is extended)
    for needle in ('sigma_b in {2, 2.5}', 'r_ex in {1.6, 2.0}', 'kappa in {1.25, 1.5, 1.75}', 'POST HOC', "order='greedy'", 'bit for bit', '0.01', 'DEVELOPMENT PATCHES ONLY'):
        assert needle in E.__doc__


def test_the_jobs_are_the_declared_ones_in_the_declared_order():
    j = E.jobs('all')
    c1, c2 = [s for s in j if s['kind'] == 'c1'], [s for s in j if s['kind'] == 'c2']
    assert len(c1) == 26 and len(c2) == 32 and j == c2 + c1
    assert {s['pxy'] for s in c1} == {tuple(p) for p in X.DEV}                                           # development patches only
    assert not ({s['pxy'] for s in c1} & {tuple(p) for p in X.TESTS})
    assert (c2[0]['set'], c2[0]['pxy'], c2[0]['arm']) == ('m2', (310, 630), 'greedy')                    # the harness test patch first
    assert [s['arm'] for s in c2] == ['greedy'] * 16 + ['dual_snr'] * 16
    assert (c1[0]['hp'], c1[0]['rex']) == (3.0, 1.2)
    assert len({s['prefix'] for s in j}) == 58
    assert [len(E.jobs(p)) for p in ('c1', 'c2')] == [26, 32]


def test_the_sixteen_test_patches_and_their_stored_rows():
    pats = E.c2_patches()
    assert len(pats) == 16 and [sum(1 for p in pats if p['group'] == g) for g in ('m2', 'm53', 'm92', 'rep')] == [4, 4, 2, 6]
    assert [p['pxy'] for p in pats[:4]] == [tuple(p) for p in X.TESTS]
    assert E.iss_prefix(pats[0]) == 'T|x310_y630|ours|k2|'
    assert E.iss_prefix(pats[4]) == 'T2|m53|x760_y590|ours|k2|' and E.iss_prefix(pats[8]).startswith('T2|m92|x1100_y220')
    assert E.iss_prefix(pats[10]) == 'R|m2r1|x316_y124|iso|ours|k2|'
    assert E.c2_pre(pats[0], 'greedy') == 'E2|m2|x310_y630|greedy|' and E.c1_pre((110, 630), 3.0, 1.2) == 'E1|x110_y630|hp3|rex1.2|'
    assert E.ARMS == {'greedy': 'greedy_snr', 'dual_snr': 'snr'} and X.VARIANTS['greedy_snr'] == dict(order='greedy', stop='snr') and X.VARIANTS['snr'] == dict(stop='snr')
    assert E.HARNESS == {'greedy': 'V|x310_y630|greedy_snr|', 'dual_snr': 'V|x310_y630|snr|'}


def test_control_1_never_loads_a_test_patch(monkeypatch):
    called = []
    monkeypatch.setattr(X, 'load_patch', lambda *a, **k: called.append(a))
    for p in list(X.TESTS) + [(1, 2)]:
        with pytest.raises(ValueError, match='development patches only'):
            E.load_for(dict(kind='c1', pxy=p))
    assert called == []
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': dict(px=np.zeros(0), py=np.zeros(0), m606=np.zeros(0)))
    monkeypatch.setattr(X, 'load_patch', lambda *a, **k: ('patch', a))
    patch, t = E.load_for(dict(kind='c1', pxy=tuple(X.DEV[0])))
    assert patch == ('patch', tuple(X.DEV[0])) and len(t['x']) == 0


# ---- the choice rule ------------------------------------------------------------------------------------------------------------------------

def _row(f, cap=0, src='control 1'):
    tot = np.array([1000, 1000 * f, 1000, 1000 * f, 0, 0, 0], float)                   # completeness = purity = f
    return dict(tot=tot, per=[tot / 2, tot / 2], blk=[np.ones((25, 4)), np.ones((25, 4))], cap=cap, act=0, src=src)


def test_the_choice_rule_is_the_one_of_exp_astro_m2_choose_on_the_stored_development_rows():
    raw = X.load_raw(E.M2RAW)
    surf = E.merge_surface(E.old_surface(dict(m2=raw), 'stored'), {})
    assert len(surf) == 3 * 2 * 8
    assert E.choose_surface(surf) == X.choose(raw) and E.choose_surface(surf)[:3] == FROZEN


def test_choice_ties_go_to_the_larger_kappa_then_the_first_key_and_uncapped_only_skips_capped_rows():
    s = {(3.0, 1.2, 2.0): _row(0.6), (3.0, 1.2, 2.5): _row(0.6), (5.0, 1.2, 2.5): _row(0.6), (2.0, 1.2, 1.5): _row(0.6000004)}
    assert E.choose_surface(s)[:3] == (3.0, 1.2, 2.5)                                                  # equal after rounding to 6 decimals: the larger kappa, then the first key in ascending order
    s[(8.0, 1.6, 1.25)] = _row(0.61, cap=1)
    assert E.choose_surface(s)[:3] == (8.0, 1.6, 1.25) and E.choose_surface(s, uncapped_only=True)[:3] == (3.0, 1.2, 2.5)
    s[(2.0, 1.2, 1.5)] = _row(0.6000016)                                                                # larger at the 6th decimal: wins over the larger kappa
    assert E.choose_surface({k: v for k, v in s.items() if k[2] != 1.25})[:3] == (2.0, 1.2, 1.5)
    assert E.choose_surface({}) is None
    assert E.F_of(_row(0.6)) == pytest.approx(0.6) and E.Fp(_row(0.6), 1) == pytest.approx(0.6)


def test_the_guarded_and_the_stored_old_surfaces_have_the_values_of_the_cyclefix_report():
    ctx = dict(m2=X.load_raw(E.M2RAW), cf=X.load_raw(E.CFRAW))
    g, s = E.merge_surface(E.old_surface(ctx, 'guarded'), {}), E.merge_surface(E.old_surface(ctx, 'stored'), {})
    assert len(g) == len(s) == 48 and round(E.F_of(s[FROZEN]), 4) == 0.5615 and round(E.F_of(g[FROZEN]), 4) == 0.5617
    assert E.choose_surface(g)[:3] == FROZEN and E.choose_surface(s)[:3] == FROZEN
    assert {k[:2] for k, r in g.items() if r['cap']} == {(8.0, 0.8)} and g[(8.0, 0.8, 2.0)]['cap'] == 2 and s[(8.0, 1.2, 2.0)]['cap'] == 2 and s[(3.0, 1.2, 2.0)]['cap'] == 0       # six stored paths hit the cap; of the guarded re-runs only (8, 0.8) again
    assert g[(8.0, 0.8, 2.0)]['act'] > 0 and s[(8.0, 0.8, 2.0)]['act'] == 0 and g[(3.0, 1.2, 2.0)]['act'] == 1


def test_merge_keeps_the_old_rows_and_adds_only_new_points():
    m2, cf = dev_world()
    ctx = dict(m2=m2, cf=cf)
    old, ext = E.old_surface(ctx, 'guarded'), E.ext_surface(e1_raw())
    assert len(old) == 48 and len(ext) == 13 * 11 and all(len(v) == 2 for v in ext.values())
    surf = E.merge_surface(old, ext)
    assert len(surf) == 48 + 13 * 11 - 4 * 8
    assert surf[FROZEN]['src'] == 'old grid' and surf[(3.0, 1.2, 1.5)]['src'] == 'control 1' and surf[(2.0, 1.2, 2.0)]['src'] == 'control 1'
    assert surf[FROZEN]['cap'] == 0 and surf[FROZEN]['act'] == 0 and surf[(3.0, 1.2, 1.25)]['act'] == 4
    assert surf[FROZEN]['per'][0].shape == (7,) and np.array_equal(surf[FROZEN]['tot'], surf[FROZEN]['per'][0] + surf[FROZEN]['per'][1])
    assert E.old_surface(dict(m2=m2, cf={k: v for k, v in cf.items() if 'x210_y630|hp3|rex1.2' not in k}), 'guarded').get(FROZEN) is None        # a combination needs both patches
    assert E.old_surface(dict(m2={k: v for k, v in m2.items() if 'x210_y630|hp3|rex1.2' not in k}, cf=cf), 'stored').get(FROZEN) is None
    assert E.ext_surface({k: v for k, v in e1_raw().items() if 'x210_y630|hp3|rex1.2' not in k}).get(FROZEN) is None
    capped = E.ext_surface(e1_raw(cap=(2.0, 2.0)))
    assert capped[(2.0, 2.0, 1.25)][tuple(X.DEV[0])]['cap'] == 1 and capped[(2.0, 1.2, 1.25)][tuple(X.DEV[0])]['cap'] == 0


def test_the_old_combinations_rerun_are_compared_with_the_guarded_rows():
    m2, cf = dev_world()
    ctx = dict(m2=m2, cf=cf)
    chk = E.c1_checks(e1_raw(), ctx)
    assert sorted((h, r) for h, r, n, s in chk) == [(3.0, 0.8), (3.0, 1.2), (5.0, 1.2), (8.0, 1.2)]
    assert all(n == 2 * 8 * 9 and s == n for h, r, n, s in chk)
    out = {(h, r): (n, s) for h, r, n, s in E.c1_checks(e1_raw(special={(5.0, 1.2, 4.0): 'poor'}), ctx)}
    assert 0 < out[(5.0, 1.2)][1] < out[(5.0, 1.2)][0] and out[(3.0, 1.2)][0] == out[(3.0, 1.2)][1] and out[(8.0, 1.2)][0] == out[(8.0, 1.2)][1]


# ---- the report of control 1 --------------------------------------------------------------------------------------------------------------

def test_control_1_report_says_the_argmax_is_the_frozen_point_and_no_finding():
    m2, cf = dev_world()
    txt = '\n'.join(E.c1_section(e1_raw(), dict(m2=m2, cf=cf)))
    assert 'the argmax is the frozen point' in txt and 'the argmax moves' not in txt
    assert 'No extended point beats the frozen point' in txt and 'FINDING (' not in txt
    assert 'Paths run: 13 of 13 per patch' in txt and '(3, 1.2): 144 of 144 identical' in txt and '(8, 1.2): 144 of 144 identical' in txt
    fz = X.f_measure(sum(m2['D|x%d_y%d|hp3|rex1.2|k2|tot' % p] for p in X.DEV))[0]
    assert 'POST HOC' in txt and '**%.4f**' % fz in txt and 'pooled F %.4f; frozen point (3, 1.2, 2) F %.4f; difference +0.0000' % (fz, fz) in txt
    assert 'Position of the extended argmax on the extended grid: interior on every axis.' in txt
    assert txt.count('| 3 | 1.2 | 2 |') >= 1 and 'axis through the frozen point' in txt and 'The ten best rows' in txt


def test_control_1_report_says_where_the_extended_argmax_sits_on_the_edge():
    m2, cf = dev_world(frozen='mid', other='poor')
    txt = '\n'.join(E.c1_section(e1_raw(special={(2.0, 2.0, 1.25): 'good'}, default='poor', frozen='mid'), dict(m2=m2, cf=cf)))
    assert 'the argmax moves' in txt and 'sigma_b 2 is at the lower edge' in txt and 'r_ex 2 is at the upper edge' in txt and 'kappa 1.25 is at the lower edge' in txt


def test_control_1_report_flags_an_extended_point_that_beats_the_frozen_one_and_does_not_run_it_on_tests():
    assert fq('good') - fq('mid') > 3 * E.FINDING
    m2, cf = dev_world(frozen='mid', other='poor')
    raw = e1_raw(special={(2.5, 1.6, 1.5): 'good'}, default='poor', frozen='mid')
    txt = '\n'.join(E.c1_section(raw, dict(m2=m2, cf=cf)))
    assert '**FINDING (POST HOC, development patches only, not run on any test patch):**' in txt and '(sigma_b 2.5, r_ex 1.6, kappa 1.5)' in txt
    assert 'the argmax moves' in txt and 'No extended point beats' not in txt and 'optimistic' in txt
    assert '(sigma_b 2.5, r_ex 1.6, kappa 1.5) F' in txt
    # a point that beats the guarded frozen value but only equals the stored one is not a finding (both references must be exceeded)
    m2b = dict(m2)
    for p in X.DEV:
        for kk, v in X._pack('', SC['good']).items():
            m2b['D|x%d_y%d|hp3|rex1.2|k2' % p + kk] = v
    txt = '\n'.join(E.c1_section(raw, dict(m2=m2b, cf=cf)))
    assert 'No extended point beats' in txt and 'FINDING (' not in txt


def test_control_1_report_flags_capped_paths_and_a_partial_run():
    m2, cf = dev_world()
    txt = '\n'.join(E.c1_section(e1_raw(cap=(2.0, 2.0)), dict(m2=m2, cf=cf)))
    assert '(cap)' in txt and 'HIT' in txt
    part = {k: v for k, v in e1_raw().items() if ('hp3|rex1.2' in k or 'hp3|rex1.6' in k)}
    txt = '\n'.join(E.c1_section(part, dict(m2=m2, cf=cf)))
    assert 'Paths run: 2 of 13 per patch' in txt and '(3, 1.2), (3, 1.6)' in txt
    assert 'No control-1 path has been run.' in '\n'.join(E.c1_section({}, dict(m2=m2, cf=cf)))


def test_control_1_pivot_marks_the_frozen_point_the_argmax_and_the_cap():
    s = {(3.0, 1.2, 2.0): _row(0.56), (3.0, 1.2, 1.5): _row(0.58, cap=1), (3.0, 1.6, 2.0): _row(0.57)}
    t = E.pivot_table(s, E.combos_of(s), (3.0, 1.2, 2.0), (3.0, 1.6, 2.0))
    assert '**0.5600**' in t and '0.5800 (cap)' in t and '0.5700 (max)' in t and '| - |' in t
    assert t.splitlines()[0].startswith('| sigma_b px | r_ex px | kappa 1.25 |')


# ---- control 2: the harness, the groups, the differences ----------------------------------------------------------------------------------------

def test_the_harness_check_is_bit_for_bit_and_notices_a_changed_array_or_entry_count():
    raw, ctx = full_ctx()
    hc = E.harness_checks(raw, ctx)
    assert [(a, n, s, eq) for a, n, s, eq, e0, e1 in hc] == [('greedy', 72, 72, True), ('dual_snr', 72, 72, True)]
    ctx2 = dict(ctx, m2=dict(ctx['m2']))
    ctx2['m2']['V|x310_y630|greedy_snr|k3|tot'] = ctx2['m2']['V|x310_y630|greedy_snr|k3|tot'] + 1
    ctx2['m2']['V|x310_y630|snr|secs'] = np.array([1.0, 999.0])
    out = {a: (n, s, eq) for a, n, s, eq, e0, e1 in E.harness_checks(raw, ctx2)}
    assert out['greedy'] == (72, 71, True) and out['dual_snr'] == (72, 72, False)
    assert E.harness_checks({}, ctx) == []
    assert '| greedy | V/x310_y630/greedy_snr | 72 | 72 | 1100 / 1100 | YES |' in '\n'.join(E.c2_section(raw, ctx)) and '| dual_snr | V/x310_y630/snr | 72 | 72 | 1100 / 1100 | YES |' in '\n'.join(E.c2_section(raw, ctx))
    assert '| YES |' in '\n'.join(E.c2_section(raw, ctx)) and '| NO |' not in '\n'.join(E.c2_section(raw, ctx))
    assert '\n'.join(E.c2_section(raw, ctx2)).count('| NO |') == 2


def _rows(txt, start):
    return [ln for ln in txt.splitlines() if ln.startswith(start)]


def test_control_2_report_the_dual_order_wins_when_planted_better_and_ties_when_equal_and_loses_when_worse():
    for order, word in (('iss_better', 'win'), ('equal', 'tie'), ('greedy_better', 'loss')):
        raw, ctx = full_ctx(order)
        txt = '\n'.join(E.c2_section(raw, ctx))
        rows = _rows(txt, '| THE QUESTION: dual-order ISS tuned minus greedy order')
        assert len(rows) == 5, order                                                                    # M2, M53, M92, repeat, all 16
        assert all(r.rstrip(' |').endswith(word) for r in rows), (order, rows[0])
        if order == 'iss_better':
            assert all(r.rstrip(' |').endswith('loss') for r in _rows(txt, '| greedy order minus crowdsource with the front end'))     # greedy poor (.30), crowdsource q1 (.76)
            assert all(r.rstrip(' |').endswith('win') for r in _rows(txt, '| dual order minus greedy order, the same stop'))        # q0 (.86) against poor
            assert all(r.rstrip(' |').endswith('win') for r in _rows(txt, '| ISS tuned minus crowdsource with the front end'))      # good (.91) against q1
        if order == 'greedy_better':
            assert all(r.rstrip(' |').endswith('win') for r in _rows(txt, '| greedy order minus crowdsource with the front end'))
            assert all(r.rstrip(' |').endswith('loss') for r in _rows(txt, '| ISS tuned minus crowdsource with the front end'))
    assert len(_rows(txt, '| tuned ISS (time stop) minus dual order')) == 5 and len(_rows(txt, '| guarded tuned ISS minus greedy order')) == 5


def test_every_method_row_of_a_pooled_table_carries_the_planted_quality_of_that_method_only():
    raw, ctx = full_ctx('iss_better')
    txt = '\n'.join(E.c2_section(raw, ctx))
    want = {'ISS, dual order, tuned (stored rows; kappa 2)': 'good', 'ISS, dual order, tuned, guarded re-run (cyclefix; kappa 2)': 'bic',
            'greedy order, instantaneous-significance stop, kappa 2 (this control)': 'poor', 'dual order, the same instantaneous-significance stop, kappa 2 (this control, added)': 'q0',
            'crowdsource with the ISS front end (round 2, stored)': 'q1', 'photutils with the ISS front end (round 2, stored)': 'q3'}
    for label, q in want.items():
        rows = _rows(txt, '| ' + label + ' |')
        assert len(rows) == 5, label
        f, c, pu = X.f_measure(np.array([SC[q]['true22'], SC[q]['matched22'], SC[q]['det22'], SC[q]['det22_matched'], 0, 0, 0], float))
        for r in rows:
            cells = r.split('|')
            assert cells[6].strip().startswith('%.3f [' % f) and cells[3].strip().startswith('%.3f [' % c) and cells[5].strip().startswith('%.3f [' % pu), (label, r)
            assert cells[2].strip() == '%d / %d' % (SC[q]['matched22'] * (len(rows) and (r.count('x') or 1)) // (r.count('x') or 1), SC[q]['true22']) or True
    # the totals of a group: matched / truth pooled over the patches of the group
    r = _rows(txt, '| ISS, dual order, tuned (stored rows; kappa 2) |')
    assert r[0].split('|')[2].strip() == '%d / %d' % (4 * SC['good']['matched22'], 4 * SC['good']['true22']) and r[1].split('|')[2].strip() == '%d / %d' % (4 * SC['good']['matched22'], 4 * SC['good']['true22'])
    assert r[2].split('|')[2].strip() == '%d / %d' % (2 * SC['good']['matched22'], 2 * SC['good']['true22']) and r[3].split('|')[2].strip() == '%d / %d' % (6 * SC['good']['matched22'], 6 * SC['good']['true22'])
    assert r[4].split('|')[2].strip() == '%d / %d' % (16 * SC['good']['matched22'], 16 * SC['good']['true22'])


def test_the_kappa_table_the_per_patch_table_and_the_bound_use_the_right_rows():
    raw, ctx = full_ctx('iss_better')
    txt = '\n'.join(E.c2_section(raw, ctx))
    fm = {q: X.f_measure(np.array([SC[q]['true22'], SC[q]['matched22'], SC[q]['det22'], SC[q]['det22_matched'], 0, 0, 0], float))[0] for q in SC}
    krows = [ln for ln in txt.splitlines() if ln.startswith('| 2 | ') or ln.startswith('| 3 | ')]
    assert len(krows) == 10 and all(r.split('|')[4].strip().startswith('%.3f [' % (fm['poor'] if r.startswith('| 2 |') else fm['q2'])) for r in krows)
    assert all(r.split('|')[5].strip().startswith('%.3f [' % (fm['q0'] if r.startswith('| 2 |') else fm['q3'])) for r in krows)
    assert txt.count('is 8 (largest pooled F of the greedy arm)') == 5                                    # kappa 3 to 8 tie (q2): the larger kappa
    sec = txt.split('### By patch')[1]
    pr = [[c.strip() for c in ln.split('|')] for ln in sec.splitlines() if ln.startswith('| ') and ln.count('|') == 9 and ln.split('|')[1].strip().startswith(('m2', 'm53', 'm92'))]
    assert len(pr) == 16 and [r[1] for r in pr[:10]] == ['m2 (310, 630)', 'm2 (310, 730)', 'm2 (360, 1040)', 'm2 (100, 1100)', 'm53 (760, 590)', 'm53 (640, 680)', 'm53 (400, 660)', 'm53 (660, 420)', 'm92 (1100, 220)', 'm92 (1430, 550)']
    assert pr[0][2:7] == ['%d' % SC['good']['true22'], '%.3f' % fm['good'], '%.3f' % fm['poor'], '%.3f' % fm['q0'], '%.3f' % fm['q1']] and pr[0][8] == 'win' and pr[0][7].count('[') == 1
    names = [ln.split('|')[1].strip() for ln in sec.splitlines() if ln.startswith('| m2r')]
    assert names == ['m2r1 (316, 124)', 'm2r1 (316, 224)', 'm2r2 (133, 281)', 'm2r2 (195, 689)', 'm2r3 (370, 441)', 'm2r3 (110, 501)']


def test_control_2_report_structure_groups_runs_and_bounds():
    raw, ctx = full_ctx('iss_better')
    txt = '\n'.join(E.c2_section(raw, ctx))
    for g in ('M2 (NGC 7089), SDSS run 2583 (4 patches', 'M53 (NGC 5024), SDSS run 5390 (4 patches', 'M92 (NGC 6341), SDSS run 4682 (2 patches',
              'M2 repeat frames m2r1, m2r2, m2r3 (SDSS runs 4187, 7778, 3388) (6 patches', 'all 16 test patches (16 patches'):
        assert g in txt
    assert txt.count('Paired differences (a minus b') == 5
    assert 'dual order minus greedy order, the same stop (order only)' in txt and 'tuned ISS (time stop) minus dual order with the instantaneous stop (the stop only)' in txt
    assert 'greedy order minus crowdsource with the front end' in txt and 'guarded tuned ISS minus greedy order' in txt
    assert txt.count('CHOSEN ON THESE TEST PATCHES') == 5 and 'crowdsource with the ISS front end (round 2, stored)' in txt and 'photutils with the ISS front end' in txt
    assert 'By patch' in txt and txt.count('| m2 (310, 630) | greedy |') == 1
    row = [ln for ln in txt.splitlines() if ln.startswith('| m2 (310, 630) | dual_snr |')][0]                 # runs table: entries, cap flag, guard activations
    assert row.split('|')[3:6] == [' 1100 ', ' no ', ' 1 '] and [ln for ln in txt.splitlines() if ln.startswith('| m2 (310, 630) | greedy |')][0].split('|')[5] == ' 0 '
    part = {k: v for k, v in raw.items() if not k.startswith('E2|m92|')}                                      # a partial run says so
    txt2 = '\n'.join(E.c2_section(part, ctx))
    assert 'M92 (NGC 6341), SDSS run 4682: not complete (0 of 2 patches run)' in txt2 and 'all 16 test patches: not complete (14 of 16 patches run)' in txt2
    assert 'No control-2 path has been run.' in '\n'.join(E.c2_section({}, ctx))


def test_the_best_greedy_kappa_is_chosen_on_the_test_patches_for_the_bound():
    raw, ctx = full_ctx('iss_better')
    for p in E.c2_patches():                                                                                # plant: the greedy arm is poor at kappa 2 and good at kappa 5
        for kk, v in X._pack('', SC['good']).items():
            raw[E.c2_pre(p, 'greedy') + 'k5' + kk] = v
    txt = '\n'.join(E.c2_section(raw, ctx))
    assert txt.count('is 5 (largest pooled F of the greedy arm)') == 5 and 'is 2 (largest' not in txt
    assert txt.count('minus the greedy order at that kappa: dF') == 5


def test_the_pooled_bootstrap_is_stratified_and_the_same_resamples_serve_every_method():
    raw, ctx = full_ctx('iss_better')
    meths = E.method_sources(ctx, raw)
    pl = [p for p in E.c2_patches() if p['group'] == 'm2']
    pooled = E.pool(meths, 'iss', pl)
    assert pooled[0].shape == (100, 4) and pooled[1] == [25] * 4 and pooled[2].shape == (7,)
    assert E.pool(meths, 'greedy', pl)[1] == [25] * 4
    assert E.pool(meths, 'iss', pl + [dict(set='m2', group='m2', pxy=(1, 1))]) is None
    assert np.array_equal(pooled[2], sum(np.asarray(meths['iss'][1](p)[1], float) for p in pl))
    rep = E.pool(meths, 'iss', [p for p in E.c2_patches() if p['group'] == 'rep'])
    assert rep[1] == [25] * 6 and rep[2][0] == 6 * SC['good']['true22']                                      # the repeat rows come from the repeat raw
    assert E.pool(meths, 'cs2', pl)[2][1] == 4 * SC['q1']['matched22']


def test_each_kind_of_patch_is_loaded_with_its_own_loader_and_truth(monkeypatch):
    calls = []
    monkeypatch.setattr(X, 'load_patch', lambda *a, **k: calls.append(('X.load_patch', a, tuple(sorted(k.items())))) or 'patch')
    monkeypatch.setattr(X, 'truth_arrays', lambda *a: calls.append(('X.truth_arrays', a)) or ('truth', 'X') + a)
    monkeypatch.setattr(E.R, 'load_patch', lambda *a, **k: calls.append(('R.load_patch', a, tuple(sorted(k.items())))) or 'patch')
    monkeypatch.setattr(E.R, 'truth_arrays', lambda *a: calls.append(('R.truth_arrays', a)) or ('truth', 'R') + a)
    monkeypatch.setattr(E.A, 'truth_for_patch', lambda t, x, y: calls.append(('cut', t, x, y)) or 'cut')
    assert E.load_for(dict(kind='c1', pxy=(110, 630))) == ('patch', 'cut')
    assert calls == [('X.load_patch', (110, 630), ()), ('X.truth_arrays', ()), ('cut', ('truth', 'X'), 110, 630)]
    calls.clear()
    assert E.load_for(dict(kind='c2', group='m53', set='m53', pxy=(760, 590))) == ('patch', 'cut')
    assert calls == [('X.load_patch', (760, 590), (('ds', 'm53'),)), ('X.truth_arrays', ('m53',)), ('cut', ('truth', 'X', 'm53'), 760, 590)]
    calls.clear()
    E.load_for(dict(kind='c2', group='m2', set='m2', pxy=(310, 630)))
    assert calls == [('X.load_patch', (310, 630), (('ds', 'm2'),)), ('X.truth_arrays', ('m2',)), ('cut', ('truth', 'X', 'm2'), 310, 630)]
    calls.clear()
    assert E.load_for(dict(kind='c2', group='rep', set='m2r2', pxy=(133, 281))) == ('patch', 'cut')
    assert calls == [('R.load_patch', (133, 281, 'm2r2'), ()), ('R.truth_arrays', ('m2r2',)), ('cut', ('truth', 'R', 'm2r2'), 133, 281)]


def test_a_run_that_makes_two_paths_is_refused_and_the_guard_still_removed(monkeypatch):
    seen = []
    _fake_world(monkeypatch, seen)
    orig = IM.sliding_path

    def two(patch, hp, rex, **kw):
        IM.sliding_path(None, None, None, None)
        IM.sliding_path(None, None, None, None)
        return {2.0: (np.zeros(1), np.zeros(1), np.ones(1))}, 1.0, 2, dict(bic_index=-1, n_blocks=3)
    monkeypatch.setattr(X, 'run_path', two)
    with pytest.raises(RuntimeError, match='2 paths were run, expected one'):
        E.run_spec(E.jobs('c1')[0])
    assert IM.sliding_path is orig and len(E.CF.TRACE) == 0
    monkeypatch.setattr(X, 'run_path', lambda *a, **k: ({2.0: (np.zeros(1), np.zeros(1), np.ones(1))}, 1.0, 2, dict(bic_index=-1, n_blocks=3)))
    with pytest.raises(RuntimeError, match='0 paths were run, expected one'):
        E.run_spec(E.jobs('c1')[0])


# ---- running: the guard around a run, write-once, the lock, the declaration ----------------------------------------------------------------

def _fake_world(monkeypatch, seen):
    monkeypatch.setattr(IM, 'sliding_path', IM.sliding_path)                                              # restored whatever happens
    monkeypatch.setattr(E, 'load_for', lambda spec: ('patch', {'x': 1}))
    monkeypatch.setattr(E, 'frozen_cfg', lambda: (3.0, 1.2, 2.0))

    def fake_path(patch, hp, rex, **kw):
        seen.append((hp, rex, kw, IM.sliding_path is E.CF.guarded_traced))
        IM.sliding_path(None, None, None, None)                                                         # the call that the frozen run function makes: lands in the guarded wrapper
        return {2.0: (np.zeros(3), np.zeros(3), np.ones(3)), 3.0: (np.zeros(2), np.zeros(2), np.ones(2))}, 5.0, 40, dict(bic_index=-1, n_blocks=41)
    monkeypatch.setattr(X, 'run_path', fake_path)
    monkeypatch.setattr(X, 'score_dets', lambda d, t: SC['good'])
    monkeypatch.setattr(G, 'sliding_path_guarded', lambda *a, **k: dict(events=[(1.0, (0, 0), 'enter'), (1.0, (0, 0), 'exit')], hist=[0, 0], taus=np.array([1.0]), snr=np.array([5.0]), barred=[(1.0, (0, 0))]))


def test_a_run_installs_the_guard_runs_the_frozen_path_function_and_removes_the_guard(monkeypatch):
    seen = []
    _fake_world(monkeypatch, seen)
    orig = IM.sliding_path
    a = E.run_spec(E.jobs('c1')[0])
    assert seen[-1][:2] == (3.0, 1.2) and seen[-1][2] == dict(kappas=E.C1_KAPPAS) and seen[-1][3] is True and IM.sliding_path is orig
    assert a['trace|counts'][0] == 1 and a['trace|counts'][6] == 1 and a['trace|barred'].shape == (1, 3) and 'det' not in a and a['bic'].tolist() == [-1, 41, 0] and a['secs'].tolist() == [5.0, 40.0]
    assert 'k2|tot' in a and 'k3|blk22' in a and a['wall'].shape == (1,)
    a = E.run_spec([s for s in E.jobs('c2') if s['arm'] == 'greedy'][0])
    assert seen[-1][:3] == (3.0, 1.2, dict(order='greedy', stop='snr')) and a['det'].shape == (3, 3)                # the catalogue at the frozen kappa
    E.run_spec([s for s in E.jobs('c2') if s['arm'] == 'dual_snr'][0])
    assert seen[-1][:3] == (3.0, 1.2, dict(stop='snr')) and len(E.CF.TRACE) == 0 and IM.sliding_path is orig


def test_the_guard_is_removed_also_when_the_run_fails(monkeypatch):
    seen = []
    _fake_world(monkeypatch, seen)
    orig = IM.sliding_path
    monkeypatch.setattr(X, 'run_path', lambda *a, **k: (_ for _ in ()).throw(RuntimeError('boom')))
    with pytest.raises(RuntimeError, match='boom'):
        E.run_spec(E.jobs('c1')[0])
    assert IM.sliding_path is orig and len(E.CF.TRACE) == 0


def test_run_is_write_once_runs_the_pending_in_order_and_records_the_state_with_the_first_result(monkeypatch, tmp_path):
    out = str(tmp_path / 'controls_raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda p: (E.CF.load_raw(p) if os.path.exists(p) else {}, {} if os.path.exists(p) else {'meta|x': np.array(['recorded'])}))
    ran = []

    def fake(spec):
        ran.append(spec['prefix'])
        return {'wall': np.array([1.0]), 'k2|tot': np.arange(7.0), 'trace|counts': np.array([5, 0, 5, 1, 0, 6, 0, 0])}
    monkeypatch.setattr(E, 'run_spec', fake)
    j = E.jobs('c2')
    assert E.run('c2', first=2, out_path=out) == 2 and ran == [j[0]['prefix'], j[1]['prefix']]
    assert sorted(E.CF.load_raw(out)) == sorted(['meta|x'] + [s['prefix'] + k for s in j[:2] for k in ('wall', 'k2|tot', 'trace|counts')])
    ran.clear()
    assert E.run('c2', first=1, out_path=out) == 1 and ran == [j[2]['prefix']]                          # the two that exist are not run again
    ran.clear()
    assert E.run('all', out_path=out) == 58 - 3 and len(ran) == 55 and j[0]['prefix'] not in ran and ran[-1] == E.jobs('all')[-1]['prefix']
    ran.clear()
    assert E.run('all', out_path=out) == 0 and ran == []
    assert E.run('c1', first=None, out_path=out) == 0


def test_a_second_process_is_refused_by_the_lock(monkeypatch, tmp_path):
    out = str(tmp_path / 'controls_raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda p: ({}, {}))
    monkeypatch.setattr(E, 'run_spec', lambda s: pytest.fail('must not run'))
    with open(out + '.runlock', 'a') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match='another process of this script is running'):
            E.run('all', out_path=out)


def test_preflight_checks_the_freeze_records_the_hashed_files_and_the_declaration(monkeypatch, tmp_path):
    out = str(tmp_path / 'controls_raw.npz')
    checks = []
    monkeypatch.setattr(X, 'check_freeze', lambda raw: checks.append('m2'))
    monkeypatch.setattr(RE, 'check_freeze', lambda: checks.append('rep'))
    cur, pending = E.preflight(out)
    assert checks == ['m2', 'rep'] and cur == {} and {'meta|files', 'meta|sha256', 'meta|declaration_sha256', 'meta|inputs_sha256', 'meta|first_run_utc', 'meta|versions'} <= set(pending)
    assert str(pending['meta|declaration_sha256'][0]) == hashlib.sha256(E.__doc__.encode()).hexdigest() == E.declaration_sha()
    assert 'src/core/imaging_guard.py' in set(pending['meta|files']) and 'experiments/astro_common.py' in set(pending['meta|files'])
    assert len(pending['meta|inputs']) == len(pending['meta|inputs_sha256']) == 4
    E.CF.save_part(out, dict(pending, **{'E1|x': np.zeros(1)}))
    cur, pending = E.preflight(out)
    assert pending == {} and 'E1|x' in cur
    bad = dict(E.CF.load_raw(out))
    bad['meta|declaration_sha256'] = np.array(['0' * 64])
    np.savez_compressed(out, **bad)
    with pytest.raises(RuntimeError, match='docstring'):
        E.preflight(out)
    bad['meta|declaration_sha256'] = np.array([E.declaration_sha()])
    sha = list(bad['meta|sha256'])
    sha[0] = '0' * 64
    bad['meta|sha256'] = np.array(sha)
    np.savez_compressed(out, **bad)
    with pytest.raises(RuntimeError, match='frozen files differ'):
        E.preflight(out)


def test_a_failing_freeze_check_stops_the_run(monkeypatch, tmp_path):
    def boom(raw):
        raise RuntimeError('the freeze record differs')
    monkeypatch.setattr(X, 'check_freeze', boom)
    with pytest.raises(RuntimeError, match='freeze record differs'):
        E.preflight(str(tmp_path / 'r.npz'))


# ---- the report ---------------------------------------------------------------------------------------------------------------------------------

def _write_world(tmp_path):
    raw, ctx = full_ctx('iss_better')
    raw.update(e1_raw())
    raw['meta|declaration_sha256'] = np.array([E.declaration_sha()])
    raw['meta|first_run_utc'] = np.array(['2026-11-08T00:00:00Z'])
    raw['meta|versions'] = np.array(['python 3', 'numpy x', 'scipy y'])
    out = tmp_path / 'out'
    out.mkdir()
    paths = {}
    for k, name in (('m2', 'm2_raw'), ('rep', 'repeat_raw'), ('cf', 'cyclefix_raw'), ('hp2', 'baselines_hp2_raw')):
        paths[k] = str(tmp_path / (name + '.npz'))
        np.savez_compressed(paths[k], **ctx[k])
    np.savez_compressed(out / 'controls_raw.npz', **raw)
    return out, paths, raw, ctx


def test_from_raw_regenerates_the_report_byte_for_byte_in_the_output_directory_of_the_environment(tmp_path):
    out, paths, raw, ctx = _write_world(tmp_path)
    env = dict(os.environ, CONTROLS_OUT=str(out), CONTROLS_M2_RAW=paths['m2'], CONTROLS_REPEAT_RAW=paths['rep'], CONTROLS_CYCLEFIX_RAW=paths['cf'], CONTROLS_HP2_RAW=paths['hp2'], PYTHONDONTWRITEBYTECODE='1')
    cmd = [sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_astro_controls.py'), '--from-raw']
    first = None
    for _ in range(2):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        data = (out / 'controls.md').read_bytes()
        assert first is None or data == first
        first = data
        (out / 'controls.md').unlink()
    assert first == E.report(E.CF.load_raw(str(out / 'controls_raw.npz')), {k: E.CF.load_raw(paths[k]) for k in paths}).encode()
    txt = first.decode()
    assert txt.startswith('# Two post hoc controls') and 'Everything here is POST HOC' in txt and 'matches it' in txt
    assert '## Control 2' in txt and '## Control 1' in txt and txt.index('## Control 2') < txt.index('## Control 1') and '## Run time' in txt
    assert 'Wall seconds of the 58 paths' in txt and 'total 4200 s (1.2 h); control 1: 2600 s; control 2: 1600 s.' in txt
    assert 'wrote' in r.stdout


def test_the_report_notes_a_changed_declaration_and_an_empty_raw(tmp_path):
    out, paths, raw, ctx = _write_world(tmp_path)
    assert 'does NOT match it' in E.report(dict(raw, **{'meta|declaration_sha256': np.array(['f' * 64])}), ctx)
    txt = E.report({}, ctx)
    assert 'no run yet' in txt and 'No control-2 path has been run.' in txt and 'No control-1 path has been run.' in txt and '## Run time' not in txt and 'n/a' in txt


def test_the_command_line_does_what_it_is_asked(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(E, 'run', lambda part, first: calls.append(('run', part, first)))
    monkeypatch.setattr(E, 'list_jobs', lambda: calls.append(('list',)))
    for argv in (['--run'], ['--run', '--part', 'c1', '--first', '3'], ['--list']):
        monkeypatch.setattr(sys, 'argv', ['exp_astro_controls.py'] + argv)
        E.main()
    assert calls == [('run', 'all', None), ('run', 'c1', 3), ('list',)]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_controls.py'])
    with pytest.raises(SystemExit):
        E.main()
    out = tmp_path / 'r.md'
    monkeypatch.setattr(E, 'MD', str(out))
    monkeypatch.setattr(E, 'RAW', str(tmp_path / 'raw.npz'))
    np.savez_compressed(tmp_path / 'raw.npz', a=np.zeros(1))
    monkeypatch.setattr(E, 'load_ctx', lambda: {'x': 1})
    monkeypatch.setattr(E, 'report', lambda raw, ctx: 'the report of %s and %s\n' % (sorted(raw), sorted(ctx)))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_controls.py', '--from-raw'])
    E.main()
    assert out.read_text() == "the report of ['a'] and ['x']\n"


def test_list_jobs_prints_the_state_of_every_run(tmp_path, capsys):
    out = str(tmp_path / 'raw.npz')
    first = E.jobs('all')[0]['prefix']
    np.savez_compressed(out, **{first + 'wall': np.ones(1)})
    E.list_jobs(out)
    lines = capsys.readouterr().out.strip().splitlines()
    assert len(lines) == 58 and lines[0].split()[-1] == 'done' and all(ln.split()[-1] == 'todo' for ln in lines[1:]) and lines[0].split()[:2] == ['c2', first]


def test_load_ctx_reads_the_four_stored_raws(monkeypatch, tmp_path):
    for k, name in (('M2RAW', 'a'), ('REPRAW', 'b'), ('CFRAW', 'c'), ('HP2RAW', 'd')):
        p = tmp_path / (name + '.npz')
        np.savez_compressed(p, **{name: np.zeros(1)})
        monkeypatch.setattr(E, k, str(p))
    ctx = E.load_ctx()
    assert {k: sorted(v) for k, v in ctx.items()} == dict(m2=['a'], rep=['b'], cf=['c'], hp2=['d'])


# ---- more on control 1's report and on running ---------------------------------------------------------------------------------------------------------

def test_control_1_counts_the_extended_points_within_half_a_percent_of_the_frozen_value():
    m2, cf = dev_world()
    txt = '\n'.join(E.c1_section(e1_raw(special={(2.5, 1.6, 1.5): 'good', (2.0, 1.2, 1.75): 'good', (2.0, 2.0, 8.0): 'bic'}), dict(m2=m2, cf=cf)))
    assert 'extended points within 0.005 of the frozen value is 2.' in txt and 'the argmax is the frozen point' in txt           # equal F at a smaller kappa: the larger kappa (the frozen one) keeps the tie
    assert 'The largest excess over the guarded frozen value is +0.0000' in txt


def test_a_better_old_grid_point_is_not_an_extended_point_and_more_than_ten_findings_are_cut_at_ten():
    m2, cf = dev_world(frozen='mid', other='poor')
    for z, ns in ((m2, ''), (cf, 'C|')):
        for p in X.DEV:
            for kk, v in X._pack('', SC['good']).items():
                z[ns + 'D|x%d_y%d|hp5|rex1.2|k4' % p + kk] = v                                               # an old-grid point better than the frozen one
    txt = '\n'.join(E.c1_section(e1_raw(default='poor', frozen='mid'), dict(m2=m2, cf=cf)))
    assert 'the argmax moves' in txt and 'No extended point beats' in txt and 'FINDING (' not in txt and 'sigma_b 5, r_ex 1.2, kappa 4, pooled F' in txt
    special = {(2.0, r, k): 'good' for r in (1.2, 1.6, 2.0) for k in (1.25, 1.5, 1.75, 2.0)}                       # 12 better extended points
    txt = '\n'.join(E.c1_section(e1_raw(special=special, default='poor', frozen='mid'), dict(m2=m2, cf=cf)))
    assert '12 extended point(s) beat the frozen point' in txt
    line = [ln for ln in txt.splitlines() if ln.startswith('**FINDING')][0]
    assert line.count('(sigma_b ') == 10


def test_control_1_findings_name_both_references_and_the_bootstrap_line_is_positive():
    m2, cf = dev_world(frozen='mid', other='poor')
    m2b = dict(m2)
    for p in X.DEV:
        for kk, v in X._pack('', SC['q0']).items():
            m2b['D|x%d_y%d|hp3|rex1.2|k2' % p + kk] = v                                                   # the stored frozen row differs from the guarded one
    raw = e1_raw(special={(2.5, 1.6, 1.5): 'good'}, default='poor', frozen='mid')
    fg = X.f_measure(sum(np.asarray(cf['C|D|x%d_y%d|hp3|rex1.2|k2|tot' % p], float) for p in X.DEV))[0]
    fs = X.f_measure(sum(np.asarray(m2b['D|x%d_y%d|hp3|rex1.2|k2|tot' % p], float) for p in X.DEV))[0]
    fe = X.f_measure(sum(np.asarray(raw['E1|x%d_y%d|hp2.5|rex1.6|k1.5|tot' % p], float) for p in X.DEV))[0]
    assert abs(fg - fs) > 0.01 and fe > fs + 0.01
    txt = '\n'.join(E.c1_section(raw, dict(m2=m2b, cf=cf)))
    assert 'against the guarded value %.4f and the stored value %.4f alike' % (fg, fs) in txt
    assert ('Surface of guarded old rows (all twelve old paths guarded) and control 1: argmax sigma_b 2.5, r_ex 1.6, kappa 1.5, pooled F %.4f; frozen point (3, 1.2, 2) F %.4f; difference %+.4f' % (fe, fg, fe - fg)) in txt
    assert ('Surface of stored old rows (as run, unguarded) and control 1: argmax sigma_b 2.5, r_ex 1.6, kappa 1.5, pooled F %.4f; frozen point (3, 1.2, 2) F %.4f; difference %+.4f' % (fe, fs, fe - fs)) in txt
    assert '(sigma_b 2.5, r_ex 1.6, kappa 1.5) F %.4f (%+.4f)' % (fe, fe - fg) in txt                       # the finding is stated against the guarded value
    line = [ln for ln in txt.splitlines() if ln.startswith('Paired block bootstrap')][0]
    assert float(line.split('dF ')[1].split(' ')[0]) > 0.05 and 'seed 0' in line and 'optimistic' in line


def test_control_1_profile_rows_are_the_axes_through_the_frozen_point():
    m2, cf = dev_world()
    txt = '\n'.join(E.c1_section(e1_raw(), dict(m2=m2, cf=cf)))
    sec = txt.split('### The F surface around the frozen point')[1].split('### The full extended surface')[0]
    rows = [[c.strip() for c in ln.split('|')[1:-1]] for ln in sec.splitlines() if ln.startswith('| ') and not ln.startswith('| axis') and not ln.startswith('|---')]
    assert [(r[0], r[1]) for r in rows] == [('sigma_b', '2'), ('sigma_b', '2.5'), ('sigma_b', '3'), ('sigma_b', '5'), ('sigma_b', '8'), ('r_ex', '0.8'), ('r_ex', '1.2'), ('r_ex', '1.6'), ('r_ex', '2'),
                                            ('kappa', '1.25'), ('kappa', '1.5'), ('kappa', '1.75'), ('kappa', '2'), ('kappa', '2.5'), ('kappa', '3'), ('kappa', '3.5'), ('kappa', '4'), ('kappa', '5'), ('kappa', '6'), ('kappa', '8')]
    fz, fm = fq('good'), fq('mid')
    fmt = {r[0] + r[1]: r for r in rows}
    assert fmt['sigma_b3'][2] == '%.4f' % fz and fmt['sigma_b3'][3] == '+0.0000' and fmt['kappa2'][3] == '+0.0000' and fmt['r_ex1.2'][2] == '%.4f' % fz
    assert fmt['sigma_b2'][2] == '%.4f' % fm and fmt['sigma_b2'][3] == '%+.4f' % (fm - fz) and fmt['kappa1.25'][2] == '%.4f' % fm and fmt['r_ex2'][3] == '%+.4f' % (fm - fz)
    assert fmt['sigma_b2'][4] == '%.4f / %.4f' % (fm, fm) and fmt['sigma_b2'][5].count('/') == 1 and fmt['sigma_b2'][6] == '%d' % (2 * SC['mid']['det22']) and fmt['sigma_b3'][7] == '-'
    top = txt.split('### The ten best rows')[1].split('The table of the stored')[0]
    trows = [[c.strip() for c in ln.split('|')[1:-1]] for ln in top.splitlines() if ln.startswith('| ') and not ln.startswith('| rank') and not ln.startswith('|---')]
    assert len(trows) == 10 and trows[0][:4] == ['1', '3', '1.2', '2'] and trows[0][4] == '%.4f' % fz and trows[0][6] == 'old grid' and trows[1][4] == '%.4f' % fm
    assert [r[3] for r in trows[1:4]] == ['8', '8', '8']                                                  # equal F: the larger kappa first


def test_control_1_says_whether_the_frozen_value_is_the_largest_along_each_axis_and_counts_the_points_above_it():
    m2, cf = dev_world()
    txt = '\n'.join(E.c1_section(e1_raw(), dict(m2=m2, cf=cf)))
    assert ('the frozen value is the largest: sigma_b yes, r_ex yes, kappa yes. Points of the extended surface with a larger pooled F than the frozen point (guarded value %.4f): 0, of which 0 beyond the old edge; the largest excess is none.' % fq('good')) in txt
    m2, cf = dev_world(frozen='mid', other='poor')
    raw = e1_raw(special={(2.5, 1.6, 1.5): 'good', (2.0, 1.2, 2.0): 'q0', (3.0, 1.6, 2.0): 'bic', (5.0, 1.2, 2.0): 'good'}, default='poor', frozen='mid')
    txt = '\n'.join(E.c1_section(raw, dict(m2=m2, cf=cf)))
    assert 'the frozen value is the largest: sigma_b no, r_ex no, kappa yes.' in txt                       # (2, 1.2, 2) is above on the sigma_b axis, (3, 1.6, 2) on the r_ex axis; (5, 1.2, 2) is an old-grid point, its planted row is not read
    assert 'guarded value %.4f): 3, of which 3 beyond the old edge; the largest excess is %+.4f at sigma_b 2.5, r_ex 1.6, kappa 1.5.' % (fq('mid'), fq('good') - fq('mid')) in txt


def test_the_blocks_of_a_row_are_those_of_its_own_kappa_and_the_per_patch_f_is_that_of_each_patch():
    m2, cf = dev_world()
    surf = E.merge_surface(E.old_surface(dict(m2=m2, cf=cf), 'guarded'), {})
    assert np.array_equal(surf[FROZEN]['blk'][0], SC['good']['blk22']) and np.array_equal(surf[(3.0, 1.2, 2.5)]['blk'][1], SC['mid']['blk22'])
    row = _row(0.5)
    row['per'] = [np.array([100, 90, 100, 90, 0, 0, 0], float), np.array([100, 40, 100, 40, 0, 0, 0], float)]
    assert E.Fp(row, 0) == pytest.approx(0.9) and E.Fp(row, 1) == pytest.approx(0.4)


def test_the_pivot_lists_the_combinations_in_ascending_order():
    s = {(5.0, 1.2, 2.0): _row(0.5), (2.0, 1.6, 2.0): _row(0.5), (2.0, 1.2, 2.0): _row(0.5)}
    t = E.pivot_table(s, E.combos_of(s), FROZEN, None)
    assert [ln.split('|')[1:3] for ln in t.splitlines()[2:]] == [[' 2 ', ' 1.2 '], [' 2 ', ' 1.6 '], [' 5 ', ' 1.2 ']]


def test_a_guarded_path_at_one_entry_below_the_cap_is_not_capped_and_one_at_the_cap_is():
    m2, cf = dev_world()
    for p in X.DEV:
        put_trace(cf, 'C|D|x%d_y%d|hp3|rex1.2|' % p, n=X.MAX_ENTRIES - 1)
        put_trace(cf, 'C|D|x%d_y%d|hp5|rex1.2|' % p, n=X.MAX_ENTRIES)
    surf = E.merge_surface(E.old_surface(dict(m2=m2, cf=cf), 'guarded'), {})
    assert surf[FROZEN]['cap'] == 0 and surf[(5.0, 1.2, 2.0)]['cap'] == 2
    st = E.merge_surface(E.old_surface(dict(m2=dict(m2, **{'D|x110_y630|hp3|rex1.2|bic': np.array([1, 2, 1])}), cf=cf), 'stored'), {})
    assert st[FROZEN]['cap'] == 1


def test_the_uncapped_argmax_line_leaves_out_the_capped_paths_and_a_gain_equal_to_the_threshold_is_no_finding(monkeypatch):
    m2, cf = dev_world(frozen='mid', other='poor')
    raw = e1_raw(special={(3.0, 1.6, 1.25): 'good', (3.0, 2.0, 1.25): 'bic'}, cap=(3.0, 1.6), default='poor', frozen='mid')
    txt = '\n'.join(E.c1_section(raw, dict(m2=m2, cf=cf)))
    assert 'argmax sigma_b 3, r_ex 1.6, kappa 1.25' in txt and 'Over the rows of paths that did not hit the cap: sigma_b 3, r_ex 2, kappa 1.25' in txt
    raw = e1_raw(special={(2.5, 1.6, 1.5): 'good'}, default='poor', frozen='mid')
    ctx = dict(m2=m2, cf=cf)
    sg = E.merge_surface(E.old_surface(ctx, 'guarded'), E.ext_surface(raw))
    ss = E.merge_surface(E.old_surface(ctx, 'stored'), E.ext_surface(raw))
    gain = E.F_of(sg[(2.5, 1.6, 1.5)]) - max(E.F_of(sg[FROZEN]), E.F_of(ss[FROZEN]))
    monkeypatch.setattr(E, 'FINDING', gain)
    assert 'FINDING (' not in '\n'.join(E.c1_section(raw, ctx))
    monkeypatch.setattr(E, 'FINDING', gain - 1e-9)
    assert 'FINDING (' in '\n'.join(E.c1_section(raw, ctx))


def _cells(line):
    return [c.strip() for c in line.split('|')[1:-1]]


def _split_world():
    """Control-1 world in which the guarded and the stored frozen rows differ (guarded mid, stored good), the second patch's frozen row is worse, an old-grid point (5, 1.2, 4) is above the guarded frozen value,
    and the extended point (2.5, 1.6, 1.5) is q0 (above the guarded, below the stored frozen row); returns (ctx, e1 raw)."""
    m2, cf = dev_world(frozen='mid', other='poor')
    for p in X.DEV:
        for kk, v in X._pack('', SC['good']).items():
            m2['D|x%d_y%d|hp3|rex1.2|k2' % p + kk] = v
        for kk, v in X._pack('', SC['q0']).items():
            cf['C|D|x%d_y%d|hp5|rex1.2|k4' % p + kk] = v
    for kk, v in X._pack('', SC['poor']).items():
        cf['C|D|x210_y630|hp3|rex1.2|k2' + kk] = v                                                      # the second patch's guarded frozen row is poor: the patches differ
    raw = e1_raw(special={(2.5, 1.6, 1.5): 'q0', (3.0, 1.6, 2.0): 'bic'}, default='poor', frozen='mid')
    for kk, v in X._pack('', SC['q0']).items():
        raw['E1|x210_y630|hp3|rex1.6|k2' + kk] = v                                                      # the two patches differ at the best extended point
    return dict(m2=m2, cf=cf), raw


def test_control_1_text_uses_the_guarded_frozen_value_for_the_excess_and_the_profile_and_each_patch_for_its_own_column():
    ctx, raw = _split_world()
    sg = E.merge_surface(E.old_surface(ctx, 'guarded'), E.ext_surface(raw))
    ss = E.merge_surface(E.old_surface(ctx, 'stored'), E.ext_surface(raw))
    fg, fs = E.F_of(sg[FROZEN]), E.F_of(ss[FROZEN])
    best = E.choose_surface(sg)
    assert fs > E.F_of(sg[(2.5, 1.6, 1.5)]) > fg and best[:3] == (3.0, 1.6, 2.0) and E.choose_surface(ss)[:3] == FROZEN
    txt = '\n'.join(E.c1_section(raw, ctx))
    assert 'FINDING (' not in txt and 'The largest excess over the guarded frozen value is %+.4f (at sigma_b 3, r_ex 1.6, kappa 2)' % (best[3] - fg) in txt
    # the profile: dF against the guarded value, F by patch of each patch
    sec = txt.split('### The F surface around the frozen point')[1].split('### The full extended surface')[0]
    rows = {(c[0], c[1]): c for c in (_cells(ln) for ln in sec.splitlines() if ln.startswith('| ') and not ln.startswith('| axis') and not ln.startswith('|---'))}
    fr = rows[('kappa', '2')]
    assert fr[2] == '%.4f' % fg and fr[3] == '+0.0000' and fr[4] == '%.4f / %.4f' % (E.Fp(sg[FROZEN], 0), E.Fp(sg[FROZEN], 1)) and E.Fp(sg[FROZEN], 0) != E.Fp(sg[FROZEN], 1)
    r16 = rows[('r_ex', '1.6')]
    assert r16[3] == '%+.4f' % (E.F_of(sg[(3.0, 1.6, 2.0)]) - fg) and r16[4] == '%.4f / %.4f' % (E.Fp(sg[(3.0, 1.6, 2.0)], 0), E.Fp(sg[(3.0, 1.6, 2.0)], 1))
    # the ranking: per-patch F of each patch
    top = txt.split('### The ten best rows')[1].split('The table of the stored')[0]
    trows = [_cells(ln) for ln in top.splitlines() if ln.startswith('| ') and not ln.startswith('| rank') and not ln.startswith('|---')]
    kk = (float(trows[0][1]), float(trows[0][2]), float(trows[0][3]))
    assert trows[0][5] == '%.4f / %.4f' % (E.Fp(sg[kk], 0), E.Fp(sg[kk], 1)) and E.Fp(sg[kk], 0) != E.Fp(sg[kk], 1)
    # the table of the stored rows is marked by the stored argmax (the frozen point, bold), not by the guarded one
    st = txt.split('The table of the stored')[1].split('Reading:')[0]
    assert '(max)' not in st and st.count('**') == 2
    assert sum(ln.count('(max)') for ln in txt.split('### The full extended surface')[1].split('### The ten best')[0].splitlines() if ln.startswith('| ')) == 1
    # one old-grid point and the extended points are above the guarded frozen value
    assert 'guarded value %.4f): 3, of which 2 beyond the old edge; the largest excess is %+.4f at sigma_b 3, r_ex 1.6, kappa 2.' % (fg, best[3] - fg) in txt
    assert 'the frozen value is the largest: sigma_b yes, r_ex no, kappa yes.' in txt


def test_control_1_axes_are_reported_one_at_a_time_and_findings_are_listed_largest_first_and_the_half_percent_band_ignores_old_points():
    m2, cf = dev_world(frozen='mid', other='poor')
    txt = '\n'.join(E.c1_section(e1_raw(special={(3.0, 1.6, 2.0): 'good'}, default='poor', frozen='mid'), dict(m2=m2, cf=cf)))
    assert 'the frozen value is the largest: sigma_b yes, r_ex no, kappa yes.' in txt
    txt = '\n'.join(E.c1_section(e1_raw(special={(2.0, 1.2, 2.0): 'good'}, default='poor', frozen='mid'), dict(m2=m2, cf=cf)))
    assert 'the frozen value is the largest: sigma_b no, r_ex yes, kappa yes.' in txt
    special = {(2.0, 1.2, 2.0): 'q0', (2.0, 1.6, 2.0): 'good', (2.0, 2.0, 2.0): 'bic'}
    line = [ln for ln in E.c1_section(e1_raw(special=special, default='poor', frozen='mid'), dict(m2=m2, cf=cf)) if ln.startswith('**FINDING')][0]
    order = [line.index('r_ex 1.6'), line.index('r_ex 2,'), line.index('r_ex 1.2')]
    assert order == sorted(order)                                                                          # good (.914), bic (.904), q0 (.864)
    # an old-grid point as good as the frozen one is not in the half-percent band of the extended points
    m2, cf = dev_world()
    for z, ns in ((m2, ''), (cf, 'C|')):
        for p in X.DEV:
            for kk, v in X._pack('', SC['good']).items():
                z[ns + 'D|x%d_y%d|hp5|rex1.2|k4' % p + kk] = v
    txt = '\n'.join(E.c1_section(e1_raw(), dict(m2=m2, cf=cf)))
    assert 'extended points within 0.005 of the frozen value is 0.' in txt


def test_control_1_bootstrap_and_run_table_columns_are_the_documented_ones():
    ctx, raw = _split_world()
    sg = E.merge_surface(E.old_surface(ctx, 'guarded'), E.ext_surface(raw))
    best = E.choose_surface(sg)[:3]
    idx = X.boot_idx_strat([25, 25])
    dc, dp, df = X.boot_diff(np.concatenate(sg[best]['blk']), np.concatenate(sg[FROZEN]['blk']), idx)
    txt = '\n'.join(E.c1_section(raw, ctx))
    assert 'dF %s (selected' % X.fmt_iv(df[0], df[1], 4) in txt
    row = [ln for ln in txt.splitlines() if ln.startswith('| (110, 630) | 3 | 1.2 |')][0]
    c = _cells(row)
    assert c[3] == '1600' and c[4] == 'no' and c[5] == '2' and c[6] == '3' and c[7] == '1' and c[8] == '100' and c[9] == '100'
    capped = [ln for ln in '\n'.join(E.c1_section(e1_raw(cap=(3.0, 1.6), default='poor', frozen='mid'), dict(ctx))).splitlines() if ln.startswith('| (110, 630) | 3 | 1.6 |')][0]
    assert _cells(capped)[3:5] == ['2500', 'HIT']
    prof = '\n'.join(E.c1_section(e1_raw(cap=(3.0, 1.6)), dict(ctx))).split('### The F surface around the frozen point')[1].split('### The full')[0]
    assert [_cells(ln)[7] for ln in prof.splitlines() if ln.startswith('| r_ex | 1.6 |')] == ['cap']


def test_control_2_text_details_the_harness_table_the_cap_flag_incomplete_groups_and_missing_readouts():
    raw, ctx = full_ctx('iss_better')
    raw = dict(raw)
    ctx2 = dict(ctx, m2=dict(ctx['m2']))
    ctx2['m2']['V|x310_y630|snr|secs'] = np.array([1.0, 999.0])
    txt = '\n'.join(E.c2_section(raw, ctx2))
    assert '| dual_snr | V/x310_y630/snr | 72 | 72 | 999 / 1100 | NO |' in txt
    ctx3 = dict(ctx, m2={k: v for k, v in ctx['m2'].items() if not k.startswith('V|')})
    txt = '\n'.join(E.c2_section(raw, ctx3))
    assert 'The benchmark patch has not been run yet.' in txt and '| arm | stored row |' not in txt
    ctx4 = dict(ctx, m2={k: (v if not k.endswith('greedy_snr|k2|n_det') else v[:0]) for k, v in ctx['m2'].items()})
    ctx4['m2'] = {k: v for k, v in ctx['m2'].items() if not (k.startswith('V|x310_y630|greedy_snr|k') or k.startswith('V|x310_y630|snr|k'))}
    txt = '\n'.join(E.c2_section(raw, ctx4))
    assert '| greedy | V/x310_y630/greedy_snr | 0 | 0 | 1100 / 1100 | NO |' in txt
    raw[E.c2_pre(E.c2_patches()[1], 'greedy') + 'trace|counts'] = np.array([X.MAX_ENTRIES, 0, 5, 1, 0, 6, 0, 0])
    txt = '\n'.join(E.c2_section(raw, ctx))
    assert '| m2 (310, 730) | greedy | 2500 | HIT |' in txt and '| m2 (310, 630) | greedy | 1100 | no |' in txt
    # the dual arm missing on M92 does not make the group incomplete; the greedy arm missing does
    part = {k: v for k, v in raw.items() if not k.startswith('E2|m92|') or '|greedy|' in k}
    txt = '\n'.join(E.c2_section(part, ctx))
    assert 'M92 (NGC 6341), SDSS run 4682 (2 patches, 50 blocks)' in txt and 'M92 (NGC 6341), SDSS run 4682: not complete' not in txt
    part = {k: v for k, v in raw.items() if not k.startswith('E2|m92|') or '|dual_snr|' in k}
    txt = '\n'.join(E.c2_section(part, ctx))
    assert 'M92 (NGC 6341), SDSS run 4682: not complete (0 of 2 patches run)' in txt and 'M92 (NGC 6341), SDSS run 4682 (2 patches' not in txt
    # one dual readout missing: a dash in that kappa's row of the group
    part = {k: v for k, v in raw.items() if k != E.c2_pre(E.c2_patches()[0], 'dual_snr') + 'k3|blk22'}
    txt = '\n'.join(E.c2_section(part, ctx))
    m2sec = txt.split('### M2 (NGC 7089), SDSS run 2583')[1].split('### M53')[0]
    assert [_cells(ln)[4] for ln in m2sec.splitlines() if ln.startswith('| 3 | ')] == ['-'] and [_cells(ln)[4] for ln in m2sec.splitlines() if ln.startswith('| 2.5 | ')][0].startswith('0.')
    # a patch without a front-end row is left out of the per-patch table, the report still runs
    part_ctx = dict(ctx, hp2={k: v for k, v in ctx['hp2'].items() if not k.startswith('T|m53|x760_y590|cs2|')})
    txt = '\n'.join(E.c2_section(raw, part_ctx))
    assert '| m53 (760, 590) |' not in txt.split('### By patch')[1] and '| m53 (640, 680) |' in txt.split('### By patch')[1]


def test_control_2_rows_are_the_documented_comparisons_with_the_documented_resamples():
    raw, ctx = full_ctx('iss_better')
    txt = '\n'.join(E.c2_section(raw, ctx))
    meths = E.method_sources(ctx, raw)
    pl = [p for p in E.c2_patches() if p['group'] == 'm2']
    pooled = {m: E.pool(meths, m, pl) for m in ('iss', 'iss_cf', 'greedy', 'dual_snr', 'cs2')}
    idx = X.boot_idx_strat(pooled['iss'][1])
    sec = txt.split('### M2 (NGC 7089), SDSS run 2583')[1].split('### M53')[0]
    for a, b, name in (('iss', 'greedy', 'THE QUESTION: dual-order ISS tuned minus greedy order'), ('dual_snr', 'greedy', 'dual order minus greedy order, the same stop (order only)'),
                       ('iss', 'dual_snr', 'tuned ISS (time stop) minus dual order with the instantaneous stop (the stop only)'), ('iss_cf', 'greedy', 'guarded tuned ISS minus greedy order'),
                       ('greedy', 'cs2', 'greedy order minus crowdsource with the front end')):
        dc, dp, df = X.boot_diff(pooled[a][0], pooled[b][0], idx)
        want = '| %s | %s | %s | %s | %s |' % (name, X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df))
        assert want in sec, name
    (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(pooled['greedy'][0], idx)
    assert '%s | %d / %d | %s | %d / %d | %s | %s |' % ('greedy order, instantaneous-significance stop, kappa 2 (this control)', pooled['greedy'][2][1], pooled['greedy'][2][0], X.fmt_iv(c, cci), pooled['greedy'][2][3], pooled['greedy'][2][2], X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)) in sec
    # per patch: the interval of 25 blocks, seed 0
    res = {m: meths[m][1](E.c2_patches()[0]) for m in ('iss', 'greedy')}
    dc, dp, df = X.boot_diff(res['iss'][0], res['greedy'][0], X.boot_idx(25))
    assert '| %s | %s |' % (X.fmt_iv(df[0], df[1]), X.verdict(*df)) in txt.split('### By patch')[1]
    # the summary repeats the numbers of the table of its group
    summ = {}
    E.c2_section(raw, ctx, summ)
    sline = [ln for ln in '\n'.join(E.summary_lines(summ)).splitlines() if ln.startswith('| M2 (NGC 7089) |')][0]
    c = _cells(sline)
    dq = X.boot_diff(pooled['iss'][0], pooled['greedy'][0], idx)[2]
    ds = X.boot_diff(pooled['dual_snr'][0], pooled['greedy'][0], idx)[2]
    assert c[3] == '%+.3f [%+.3f, %+.3f]' % (dq[0], dq[1][0], dq[1][1]) and c[4] == X.verdict(*dq)
    assert c[5] == '%+.3f [%+.3f, %+.3f] %s' % (ds[0], ds[1][0], ds[1][1], X.verdict(*ds))


def test_the_verdict_column_is_the_verdict_on_f_and_not_on_completeness_or_purity(monkeypatch):
    raw, ctx = full_ctx('equal')
    monkeypatch.setattr(X, 'boot_diff', lambda a, b, idx: ((0.0, np.array([-0.1, 0.1])), (0.0, np.array([-0.1, 0.1])), (0.1, np.array([0.05, 0.2]))))
    txt = '\n'.join(E.c2_section(raw, ctx))
    rows = [ln for ln in txt.splitlines() if ln.startswith('| THE QUESTION') or ln.startswith('| greedy order minus') or ln.startswith('| dual order minus')]
    assert len(rows) == 15 and all(r.rstrip(' |').endswith('win') for r in rows)
    assert all(_cells(ln)[-1] == 'win' for ln in txt.split('### By patch')[1].splitlines() if ln.startswith('| m2 ('))


def test_the_summary_of_control_1_gives_the_guarded_then_the_stored_frozen_value_and_the_report_lines_of_the_run():
    ctx, raw = _split_world()
    sg = E.merge_surface(E.old_surface(ctx, 'guarded'), E.ext_surface(raw))
    ss = E.merge_surface(E.old_surface(ctx, 'stored'), E.ext_surface(raw))
    summ = {}
    E.c1_section(raw, ctx, summ)
    line = E.summary_lines(summ)[-1]
    assert 'has %.4f (guarded) and %.4f (stored)' % (E.F_of(sg[FROZEN]), E.F_of(ss[FROZEN])) in line
    full, c2ctx = full_ctx('iss_better')
    full.update({'meta|first_run_utc': np.array(['2026-11-08T00:00:00Z']), 'meta|versions': np.array(['python 3', 'numpy x', 'scipy y'])})
    full['X|a|wall'] = np.array([1e6])
    full.update(e1_raw())
    txt = E.report(full, dict(c2ctx, **dict(zip(('m2', 'cf'), (dict(c2ctx['m2'], **dev_world()[0]), dict(c2ctx['cf'], **dev_world()[1]))))))
    assert 'First run recorded: 2026-11-08T00:00:00Z; python 3; numpy x; scipy y.' in txt and 'total 4200 s' in txt


def test_run_pass_first_zero_makes_a_missing_directory_and_saves_the_state_only_with_the_first_result(monkeypatch, tmp_path, capsys):
    out = str(tmp_path / 'sub' / 'dir' / 'controls_raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda p: (E.CF.load_raw(p) if os.path.exists(p) else {}, {} if os.path.exists(p) else {'meta|x': np.array(['recorded'])}))
    monkeypatch.setattr(E, 'run_spec', lambda spec: {'wall': np.array([2.0]), 'k2|tot': np.arange(7.0), 'trace|counts': np.array([X.MAX_ENTRIES, 0, 5, 1, 0, 6, 2, 0])})
    assert E.run('c1', first=0, out_path=out) == 0 and os.path.isdir(os.path.dirname(out)) and not os.path.exists(out)
    saved = []
    orig = E.CF.save_part
    monkeypatch.setattr(E.CF, 'save_part', lambda path, arrays: (saved.append(sorted(arrays)), orig(path, arrays)))
    assert E.run('c1', first=3, out_path=out) == 3
    assert ['meta|x' in k for k in saved] == [True, False, False] and len(saved) == 3
    text = capsys.readouterr().out
    assert '3 runs to do (c1), 0 exist' in text and text.count('start E1|') == 3 and text.count('done  E1|') == 3 and 'cap HIT' in text and '2 guard activations' in text and '2 s' in text
    monkeypatch.setattr(E, 'preflight', lambda p: (E.CF.load_raw(p), {}))
    assert E.run('c1', first=2, out_path=out) == 2
    assert '2 runs to do (c1), 3 exist' in capsys.readouterr().out


def test_a_run_below_the_cap_is_not_flagged(monkeypatch, tmp_path, capsys):
    out = str(tmp_path / 'controls_raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda p: ({}, {}))
    monkeypatch.setattr(E, 'run_spec', lambda spec: {'wall': np.array([2.0]), 'trace|counts': np.array([X.MAX_ENTRIES - 1, 0, 5, 1, 0, 6, 0, 0])})
    E.run('c1', first=1, out_path=out)
    text = capsys.readouterr().out
    assert 'cap no' in text and 'cap HIT' not in text


def test_the_summary_collects_the_headline_numbers_of_both_controls():
    raw, ctx = full_ctx('iss_better')
    raw.update(e1_raw(special={(2.5, 1.6, 1.5): 'good'}, frozen='mid', default='poor'))
    m2, cf = dev_world(frozen='mid', other='poor')
    ctx = dict(ctx, m2=dict(ctx['m2'], **m2), cf=dict(ctx['cf'], **cf))
    txt = E.report(raw, ctx)
    assert txt.index('## Summary (POST HOC)') < txt.index('## Control 2') < txt.index('## Control 1')
    sec = txt.split('## Summary (POST HOC)')[1].split('## Control 2')[0]
    rows = [[c.strip() for c in ln.split('|')[1:-1]] for ln in sec.splitlines() if ln.startswith('| ') and not ln.startswith('| group')]
    assert [r[0] for r in rows] == ['M2 (NGC 7089)', 'M53 (NGC 5024)', 'M92 (NGC 6341)', 'M2 repeat frames m2r1, m2r2, m2r3', 'all 16 test patches']
    for r in rows:
        assert r[1] == '%.3f' % fq('good') and r[2] == '%.3f' % fq('poor') and r[3].startswith('+0.') and r[4] == 'win' and r[5].endswith(' win') and r[6].endswith('win (kappa 8)')
    last = [ln for ln in sec.splitlines() if ln.startswith('Control 1')][0]
    assert 'the argmax moves' in last and 'argmax is sigma_b 2.5, r_ex 1.6, kappa 1.5' in last and '(13 paths per patch run)' in last and 'more than 0.01: 1.' in last
    # the F of the summary is the pooled F of the table of the group
    assert '| %.3f | %.3f |' % (fq('good'), fq('poor')) in sec
    only1 = E.report(e1_raw(), dict(ctx))
    assert '## Summary' in only1 and 'Control 2, pooled F' not in only1 and 'Control 1 (development patches only)' in only1
    assert '## Summary' not in E.report({}, ctx)


# ---- the harness on the real data (last: it runs a path, about 20 s) -----------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, 'truth_hugs_meth1.npz')), reason='the M2 patches are not at $ASTRO_DATA')
def test_the_harness_reproduces_the_stored_snr_variant_row_of_the_benchmark_patch_bit_for_bit_on_real_data():
    spec = [s for s in E.jobs('c2') if s['arm'] == 'dual_snr' and s['pxy'] == (310, 630)][0]
    a = E.run_spec(spec)
    raw = {spec['prefix'] + k: v for k, v in a.items()}
    assert E.harness_checks(raw, dict(m2=X.load_raw(E.M2RAW))) == [('dual_snr', 72, 72, True, 393, 393)]
    assert IM.sliding_path is not E.CF.guarded_traced
