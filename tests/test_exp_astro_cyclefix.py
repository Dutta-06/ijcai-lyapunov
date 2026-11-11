"""The guarded re-runs of the capped astronomy paths and the check of every other path (experiments/exp_astro_cyclefix.py): the inventory of the stored runs and the scan for the entry cap, the
replacement of `imaging.sliding_path` at every call site of the frozen modules, the layout of a re-run in the frozen scripts' own keys, the comparison with the stored run, write-once running,
the freeze check, and the report read off a planted raw (labels, the as-run and the guarded rows side by side, the development choice, the corrected rows, the check-mode table), with
the byte identity of `--from-raw`."""
import functools
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_common as A                # noqa: E402
import astro_repeat as R                # noqa: E402
import exp_astro_cyclefix as E          # noqa: E402
import exp_astro_m2 as X                # noqa: E402
import exp_astro_psf as PS              # noqa: E402
import exp_astro_repeat as RE           # noqa: E402
from src.core import imaging as IM      # noqa: E402
from src.core import imaging_guard as G  # noqa: E402
from test_astro_report import _synthetic_catalogs   # noqa: E402

CFG = dict(hp=3.0, rex=1.2, kap=2.0, default=(5.0, 1.2, 5.0))
STAGE1 = ['D|x110_y630|hp5|rex0.8', 'D|x110_y630|hp8|rex0.8', 'D|x110_y630|hp8|rex1.2', 'D|x210_y630|hp5|rex0.8', 'D|x210_y630|hp8|rex0.8', 'D|x210_y630|hp8|rex1.2',
          'T2|m92|x1100_y220|default', 'T2|m92|x1430_y550|default', 'V|x310_y630|norex', 'R|m2r2|x133_y281|iso|default', 'R|m2r3|x370_y441|iso|default', 'R|m2r3|x110_y501|iso|default']
SAME = {('m2r1', (316, 124)): 1, ('m2r1', (316, 224)): 1, ('m2r2', (133, 281)): 0, ('m2r2', (195, 689)): 0, ('m2r3', (370, 441)): 1, ('m2r3', (110, 501)): 1}
TRUTH, CATS = _synthetic_catalogs(np.random.default_rng(11))
SC = {k: A.score(np.c_[v[0], v[1]], v[2], TRUTH) for k, v in CATS.items()}


def _run(key, n=1200, cap=False, det=None, rsd=None, bic=True, k0=True, secs=100.0, best=None, bic_key=None):
    """The arrays of a stored run (keys relative to the run's prefix) in the layout of the frozen scripts: the scores of every readout (those of the planted catalogue `key`, those of kappa 2 of the
    catalogue `best` if given), [seconds, entries], [BIC block, blocks, cap], the catalogue."""
    out = {}
    for k in X.KAPPAS + ((0.0,) if k0 else ()):
        out.update(X._pack('k%g' % k, SC[best if (best and k == 2.0) else (bic_key if (bic_key and k == 0.0) else key)]))
    out['secs'] = np.array([secs, float(X.MAX_ENTRIES if cap else n)])
    if bic:
        out['bic'] = np.array([n - 30, X.MAX_ENTRIES + 1 if cap else n + 1, int(cap)])
    if det:
        out['det'] = np.c_[CATS[det][0], CATS[det][1], CATS[det][2]]
    if rsd is not None:
        out['rsd'] = np.array(rsd)
    return out


def _put(raw, rid, arrays):
    for f, v in arrays.items():
        raw[rid + '|' + f] = v


def _baselines(raw, pre):
    for label, key in (('photutils (default)', 'poor'), ('photutils (tuned)', 'mid'), ('crowdsource (default)', 'poor'), ('crowdsource (tuned)', 'mid')):
        for kk, v in X._pack(label, SC[key]).items():
            raw[pre + kk] = v
        if '(tuned)' in label:
            raw[pre + label.split(' ')[0] + ' (tuned)|name'] = np.array(['variant'])


def world():
    """The three frozen raws as the real constants lay them out (every run of the real inventory), with the 12 capped runs capped; a fresh shallow copy of a template; dict(m2, psf, rep)."""
    t = _world_template()
    return {k: dict(v) for k, v in t.items()}


@functools.lru_cache(maxsize=1)
def _world_template():
    m2, psf, rep = {}, {}, {}
    for p in X.DEV:
        for hp in X.HPS:
            for rex in X.REXS:
                rid = 'D|x%d_y%d|hp%g|rex%g' % (p[0], p[1], hp, rex)
                _put(m2, rid, _run('mid', 1500 + int(hp * 10), cap=rid in STAGE1, secs=300.0, best='good' if (hp, rex) == (3.0, 1.2) else None))        # the frozen point is the best row
    for ds, kfmt, plist in E.test_groups():
        for p in plist:
            pre = kfmt % tuple(p)
            m2[pre + 'cfg'] = np.array([3.0, 1.2, 2.0])
            _put(m2, pre + 'ours', _run('good', 1300, det='good'))
            _put(m2, pre + 'default', _run('mid', 1400, cap=(pre + 'default') in STAGE1, det='mid', bic_key='bic'))             # the BIC readout differs from the others
            _baselines(m2, pre)
    for name in X.VARIANTS:
        rid = 'V|x310_y630|%s' % name
        a = _run('mid', 1100, cap=rid in STAGE1, bic=False, k0=False, secs=200.0 + 10.0 * list(X.VARIANTS).index(name))
        _put(m2, rid, a)
    for ds in PS.CLUSTERS:
        for p in X.SECOND[ds]['patches']:
            psf['P|%s|x%d_y%d|cfg' % (ds, p[0], p[1])] = np.array([3.0, 1.2, 2.0])
            for label in ('ours', 'default'):
                _put(psf, 'P|%s|x%d_y%d|%s' % (ds, p[0], p[1], label), _run('good', 1350, det='good', rsd=[1.8, 1.2]))
    for name, p in E.RE.PATCHES:
        pre = RE.pre_of(name, p)
        same = SAME[(name, tuple(p))]
        rep[pre + 'cfg'] = np.array([3.0, 1.2, 2.0])
        rep[pre + 'kernel|same'] = np.array([same])
        rep[pre + 'scale'] = np.array([1.3])
        for lab in ('ours', 'default'):
            iso = _run('mid', 1250, cap=(pre + 'iso|' + lab) in STAGE1, det='mid', rsd=[1.2])
            _put(rep, pre + 'iso|' + lab, iso)
            _put(rep, pre + 'adq|' + lab, iso if same else _run('good', 1260, det='good', rsd=[1.1]))
            _put(rep, pre + 'scl|' + lab, _run('good', 1000, det='good', rsd=[1.1]))
        _baselines(rep, pre)
    return dict(m2=m2, psf=psf, rep=rep)


def cells(row):
    """The cells of a markdown table row (pipes escaped as \\| do not split)."""
    return [c.strip().replace('\\|', '|') for c in re.split(r'(?<!\\)\|', row.strip())[1:-1]]


def _trace(n, act=0, bars=None, most=1, taus=True):
    cnt = np.array([n, 0, n, most, 0, n + 1, act, 0])
    out = {'trace|counts': cnt, 'trace|barred': np.array(bars, float).reshape(-1, 3) if bars else np.zeros((0, 3)), 'wall': np.array([12.0])}
    if taus:
        out['trace|taus'] = np.linspace(0.01, 1.0, n)
    return out


def guarded(fro, rids=STAGE1, key='good', cap=False, act=0):
    """The cyclefix raw with the guarded re-runs of `rids`: scores of the planted catalogue `key`, entries 1700 (the cap if cap), `act` activations."""
    cf = {'meta|files': np.array(['src/core/imaging.py']), 'meta|sha256': np.array(['0' * 64]), 'meta|plan_sha256': np.array(['1' * 64]), 'meta|versions': np.array(['python 3.13.0', 'numpy 2.5.4', 'scipy 1.18.1'])}
    for rid in rids:
        arr = _run(key, 1700, cap=cap, det=key if rid.split('|')[0] in ('T2', 'R') else None, rsd=[1.3] if rid.split('|')[0] in ('T2', 'R') else None, bic=rid != 'V|x310_y630|norex' or True,
                   bic_key='poor' if key != 'poor' else 'good')
        arr.update(_trace(X.MAX_ENTRIES if cap else 1700, act=act, bars=[[0.5 + i, 3, 4 + i] for i in range(act)]))
        for f, v in arr.items():
            cf['G|' + rid + '|' + f] = v
    return cf


def checked(fro, cf, rids=None, differ=(), act=None):
    """Adds check-mode runs: copies of the stored runs (identical), except those in `differ` (another catalogue and one more entry); `act`: {rid: activations}."""
    aff, chk, alias = E.split_inventory(fro, CFG)
    for s in chk:
        if rids is not None and s['rid'] not in rids:
            continue
        fields = E.stored_fields(fro[s['raw']], s['rid'])
        fields = {f: v.copy() for f, v in fields.items()}
        if s['rid'] in differ:
            n = int(fields['secs'][1])
            fields.update(_run('poor', n + 1, det='poor' if 'det' in fields else None, bic='bic' in fields, k0='k0|tot' in fields))
            if 'rsd' in fields:
                fields['rsd'] = fields['rsd'] + 0.1
        a = (act or {}).get(s['rid'], 0)
        fields.update(_trace(int(fields['secs'][1]), act=a, bars=[[0.4, 1, 2]] * a, taus=False))
        for f, v in fields.items():
            cf['C|' + s['rid'] + '|' + f] = v
    return cf


def test_the_hashed_files_are_the_frozen_ones_the_guard_and_every_file_the_freeze_records_hash():
    assert 'src/core/imaging_guard.py' in E.HASHED and set(E.FROZEN) <= set(E.HASHED) and len(E.HASHED) == len(E.FROZEN) + 1 and 'src/core/imaging.py' in E.FROZEN
    for rec in ('freeze.json', 'freeze_repeat.json'):
        path = os.path.join(E.ASTRO, rec)
        if os.path.exists(path):
            files = json.load(open(path))['files']
            assert all(f in E.FROZEN for f in files), [f for f in files if f not in E.FROZEN]


# ---- the inventory and the scan ---------------------------------------------------------------------------------------------------------

def test_the_inventory_lists_every_stored_run_once_and_the_aliased_runs_are_not_listed():
    fro = world()
    inv = E.inventory(fro, CFG)
    rids = [s['rid'] for s in inv]
    assert len(rids) == len(set(rids)) == 77 and {f: sum(1 for s in inv if s['fam'] == f) for f in E.FAMILIES} == dict(D=12, T=20, V=5, P=12, R=28)
    assert not any('|adq|' in r and ('m2r1' in r or 'm2r3' in r) for r in rids)                                            # the rule keeps the isotropic kernel there: one run, listed as 'iso'
    assert 'R|m2r2|x133_y281|adq|default' in rids and 'R|m2r2|x195_y689|adq|ours' in rids
    s = {x['rid']: x for x in inv}
    assert (s['T2|m92|x1100_y220|ours']['h'], s['T2|m92|x1100_y220|ours']['r'], s['T2|m92|x1100_y220|ours']['dk']) == (3.0, 1.2, 2.0)
    assert (s['T2|m92|x1100_y220|default']['h'], s['T2|m92|x1100_y220|default']['r'], s['T2|m92|x1100_y220|default']['dk']) == (5.0, 1.2, 5.0)
    assert (s['D|x110_y630|hp8|rex0.8']['h'], s['D|x110_y630|hp8|rex0.8']['r'], s['D|x110_y630|hp8|rex0.8']['dk']) == (8.0, 0.8, None)
    assert (s['V|x310_y630|nohp']['h'], s['V|x310_y630|nohp']['r']) == (3.0, 1.2)
    assert (s['R|m2r2|x133_y281|scl|ours']['h'], s['R|m2r2|x133_y281|scl|ours']['r']) == (3.0 * 1.3, 1.2 * 1.3)           # the stored scale (the seeing) multiplies both lengths
    assert (s['R|m2r2|x133_y281|scl|default']['h'], s['R|m2r2|x133_y281|scl|default']['dk']) == (5.0 * 1.3, 5.0)
    assert s['P|m53|x760_y590|default']['froz'] == 'T2|m53|x760_y590|default|det' and s['P|m53|x760_y590|default']['raw'] == 'psf' and s['R|m2r1|x316_y124|iso|ours']['raw'] == 'rep'
    assert s['D|x110_y630|hp3|rex0.8']['cost'] == 300.0 and s['V|x310_y630|nohp']['cost'] == 210.0


def test_the_scan_finds_the_runs_that_hit_the_cap_by_the_bic_flag_and_by_the_entry_count():
    fro = world()
    cap = E.capped(fro)
    assert sorted(cap) == sorted(STAGE1 + ['R|m2r3|x370_y441|adq|default', 'R|m2r3|x110_y501|adq|default'])
    assert {r: h for r, h in cap.items() if h == 'entries'} == {'V|x310_y630|norex': 'entries'}                           # no `bic` key: found by the stored entry count
    aff, chk, alias = E.split_inventory(fro, CFG)
    assert sorted(s['rid'] for s in aff) == sorted(STAGE1) and len(chk) == 65 and alias == ['R|m2r3|x110_y501|adq|default', 'R|m2r3|x370_y441|adq|default']
    assert {s['rid']: s['how'] for s in aff}['D|x110_y630|hp5|rex0.8'] == 'bic'
    fro['m2']['D|x110_y630|hp3|rex0.8|bic'] = np.array([5, 2501, 1])                                                     # one more flag: one more affected path
    assert 'D|x110_y630|hp3|rex0.8' in E.capped(fro)
    fro['m2']['V|x310_y630|nohp|secs'] = np.array([10.0, 2499.0])                                                         # one entry below the cap is not capped
    assert 'V|x310_y630|nohp' not in E.capped(fro)
    fro['m2']['V|x310_y630|nohp|secs'] = np.array([10.0, 2500.0])
    assert E.capped(fro)['V|x310_y630|nohp'] == 'entries'
    fro['m2']['V|x310_y630|nohp|bic'] = np.array([5, 2501, 0])                                                            # a record that says 'not capped' is trusted over the entry count
    assert 'V|x310_y630|nohp' not in E.capped(fro)


@pytest.mark.skipif(not os.path.exists(E.M2RAW) or not os.path.exists(E.REPRAW), reason='the frozen raws are not there')
def test_the_scan_of_the_real_raws_finds_the_declared_paths():
    fro = E.frozen_raws()
    aff, chk, alias = E.split_inventory(fro, E.cfg_from_raw(fro['m2']))
    assert sorted(s['rid'] for s in aff) == sorted(STAGE1) and len(chk) == 65 and len(E.inventory(fro, E.cfg_from_raw(fro['m2']))) == 77
    assert E.cfg_from_raw(fro['m2']) == dict(CFG, default=tuple(X.DEFAULT))
    assert {s['rid'] for s in aff if s['how'] == 'entries'} == {'V|x310_y630|norex'} and not any(s['fam'] == 'P' for s in aff)       # psf_raw has none; the scan adds the variant


# ---- the replacement of the path function at every call site -----------------------------------------------------------------------------

class _Hit(Exception):
    pass


def _tiny_patch(n=40):
    yy, xx = np.mgrid[0:51, 0:51]
    stamp = np.exp(-0.5 * ((yy - 25) ** 2 + (xx - 25) ** 2) / 1.5 ** 2)
    return dict(dn=100.0 + np.random.default_rng(0).normal(0, 3, (n, n)), psf=stamp / stamp.sum(), gain=4.0, dark=1.0)


def test_every_call_site_of_the_frozen_modules_looks_the_path_function_up_at_call_time(monkeypatch):
    """`exp_astro_m2.run_path`, `exp_astro_psf.run_path_psf` (which the repeat code calls) and `imaging.fit_frame` (the SMLM frames) reach the replacement."""
    seen = []

    def sentinel(*a, **k):
        seen.append(len(a))
        raise _Hit

    monkeypatch.setattr(IM, 'sliding_path', sentinel)
    patch = _tiny_patch()
    with pytest.raises(_Hit):
        X.run_path(patch, 3.0, 1.2)
    fit = E.PG.fit_adequate(patch['psf'])
    with pytest.raises(_Hit):
        PS.run_path_psf(patch, 3.0, 1.2, fit)
    psf = IM.PixelPSF(24, 24, 1.3)
    with pytest.raises(_Hit):
        IM.fit_frame(np.full((24, 24), 50.0), psf, IM.Lattice(psf, h=0.5, margin=1.0), bg=50.0, var_a=1.0)
    assert len(seen) == 3


def test_install_puts_the_guarded_wrapper_in_place_and_uninstall_restores_the_frozen_function(monkeypatch):
    monkeypatch.setattr(IM, 'sliding_path', IM.sliding_path)
    orig = E.ORIG_PATH
    assert orig is not E.guarded_traced
    E.install()
    assert IM.sliding_path is E.guarded_traced
    E.uninstall()
    assert IM.sliding_path is orig and orig.__module__ == 'src.core.imaging'


def test_the_wrapper_calls_the_guarded_function_with_the_arguments_it_got_and_keeps_its_trace(monkeypatch):
    got = []

    def fake(*a, **k):
        got.append((a, k))
        return dict(events=[(0.1, (1, 2), 'enter'), (0.1, (1, 2), 'exit'), (0.1, (1, 2), 'enter'), (0.1, (1, 2), 'exit')], hist=[0, 1, 2], taus=np.array([0.1, 0.1]), snr=np.array([5.0, 5.0]), barred=[(0.1, (1, 2))])

    monkeypatch.setattr(G, 'sliding_path_guarded', fake)
    E.TRACE.clear()
    res = E.guarded_traced('psf', 'y', max_entries=7)
    assert got == [(('psf', 'y'), {'max_entries': 7})] and res['barred'] == [(0.1, (1, 2))]
    t = E.TRACE[-1]
    assert t['counts'].tolist() == [2, 2, 1, 2, 0, 3, 1, 1] and t['barred'].tolist() == [[0.1, 1, 2]] and t['taus'].tolist() == [0.1, 0.1]
    E.TRACE.clear()


def test_trace_of_counts_the_events_the_points_and_the_reentries_at_the_exit_time():
    ev = [(0.1, (1, 1), 'enter'), (0.1, (1, 1), 'exit'), (0.1, (1, 1), 'enter'), (0.2, (2, 2), 'enter'), (0.2, (1, 1), 'exit'), (0.3, (1, 1), 'enter'), (0.4, (3, 3), 'enter'), (0.4, (3, 3), 'exit'),
          (0.4, (3, 3), 'enter'), (0.4, (3, 3), 'exit'), (0.4, (3, 3), 'enter'), (0.9, (4, 4), 'enter')]
    taus = np.array([e[0] for e in ev if e[2] == 'enter'])
    t = E.trace_of(dict(events=ev, hist=list(range(8)), taus=taus, snr=taus, barred=[(0.4, (3, 3))]))
    # entries 8; exit events 4; distinct points 4; the most entries of one point: 3; points entered 3 or more times: (1,1) and (3,3); blocks 8; activations 1; re-entries at the exit time:
    # (1,1) at 0.1 (after its exit at 0.1), (3,3) at 0.4 twice -> 3 (the entry of (1,1) at 0.3 follows its exit at 0.2: another time)
    assert t['counts'].tolist() == [8, 4, 4, 3, 2, 8, 1, 3] and t['barred'].tolist() == [[0.4, 3, 3]] and len(E.COUNT_FIELDS) == 8
    assert t['taus'].tolist() == taus.tolist()
    t0 = E.trace_of(dict(events=[], hist=[0], taus=np.zeros(0), snr=np.zeros(0), barred=[]))
    assert t0['counts'].tolist() == [0, 0, 0, 0, 0, 1, 0, 0] and t0['barred'].shape == (0, 3)


# ---- the layout of a re-run -------------------------------------------------------------------------------------------------------------

def _stub_world(monkeypatch):
    """Planted data and run functions: the frozen run functions return the catalogues of the planted world and record their arguments."""
    calls = []
    monkeypatch.setattr(X, 'load_patch', lambda x0, y0, ds='m2': dict(tag=('m2patch', ds, x0, y0), dn=np.full((100, 100), 150.0), psf=np.full((51, 51), 1.0 / 2601), gain=4.6, dark=1.0))
    monkeypatch.setattr(R, 'load_patch', lambda x0, y0, name: dict(tag=('rpatch', name, x0, y0), dn=np.full((100, 100), 150.0), psf=np.full((51, 51), 1.0 / 2601), gain=4.6, dark=1.0))
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': dict(px=TRUTH['x'] + 100.0, py=TRUTH['y'] + 630.0, m606=TRUTH['mag']))
    monkeypatch.setattr(R, 'truth_arrays', lambda name: dict(px=TRUTH['x'] + 100.0, py=TRUTH['y'] + 630.0, m606=TRUTH['mag']))
    monkeypatch.setattr(A, 'truth_for_patch', lambda truth, x0, y0, n=100, mag_max=30.0: TRUTH)
    monkeypatch.setattr(IM, 'fit_psf_parts', lambda stamp, n_parts=3, sigma0=1.0: (0.9, ((0.53, 1.0), (0.35, 1.9), (0.12, 4.7)), 0.001))
    monkeypatch.setattr(PS, 'residual_sd', lambda patch, P, det: 1.05 + 0.01 * len(det) + 0.001 * (P.hp or 0.0))
    monkeypatch.setattr(E.PG, 'fit_adequate', lambda stamp: dict(kind='general', sigma=0.9, parts=None, parts5=((1.0, 1.0, 1.0, 0.0, 0.0),), rms_isotropic=0.004, rms=0.001))

    def run_path(patch, h, r, **kw):
        calls.append(('run_path', patch['tag'], h, r, tuple(sorted(kw.items()))))
        return _dets(0), 1.5, 40, dict(bic_index=3, n_blocks=9)

    def run_path_psf(patch, h, r, fit, **kw):
        calls.append(('run_path_psf', patch['tag'], h, r, fit['kind']))
        return _dets(1), 2.5, 41, dict(bic_index=4, n_blocks=9, rsd={k: 1.1 + 0.1 * i for i, k in enumerate(X.KAPPAS + (0.0,))})
    monkeypatch.setattr(X, 'run_path', run_path)
    monkeypatch.setattr(PS, 'run_path_psf', run_path_psf)
    return calls


def _dets(shift):
    """A different catalogue for every readout (so that a wrong reported kappa shows): {kappa: (x, y, flux)}; `shift` rotates the choice (the two run functions return different ones)."""
    names = ['good', 'mid', 'poor', 'bic']
    return {k: (CATS[names[(i + shift) % 4]][0], CATS[names[(i + shift) % 4]][1], CATS[names[(i + shift) % 4]][2]) for i, k in enumerate(X.KAPPAS + (0.0,))}


def _spec(rid, **kw):
    fro = world()
    s = next(x for x in E.inventory(fro, CFG) if x['rid'] == rid)
    return dict(s, **kw)


def _expected(shift, n, secs=1.5, bic=(3, 9, 0)):
    """What the frozen scripts store for a run whose readouts are the catalogues of `_dets(shift)` (no catalogue at a reported kappa: see `_det_of`)."""
    names = ['good', 'mid', 'poor', 'bic']
    out = {}
    for i, k in enumerate(X.KAPPAS + (0.0,)):
        out.update(X._pack('k%g' % k, SC[names[(i + shift) % 4]]))
    out['secs'] = np.array([secs, n])
    out['bic'] = np.array(bic)
    return out


def _det_of(shift, dk):
    d = _dets(shift)[dk]
    return np.c_[d[0], d[1], d[2]]


def test_a_development_run_is_stored_as_part_D_stores_it(monkeypatch):
    calls = _stub_world(monkeypatch)
    arr = E.run_path_spec(_spec('D|x210_y630|hp8|rex1.2'))
    assert calls == [('run_path', ('m2patch', 'm2', 210, 630), 8.0, 1.2, ())]
    exp = _expected(0, 40)
    assert sorted(arr) == sorted(exp) and all(np.array_equal(arr[k], exp[k], equal_nan=True) for k in exp)
    assert 'det' not in arr and np.array_equal(arr['bic'], [3, 9, 0])                                                      # part D stores no catalogue


def test_a_test_patch_run_is_stored_as_the_test_part_stores_it_and_gets_the_isotropic_residual_scale(monkeypatch):
    calls = _stub_world(monkeypatch)
    arr = E.run_path_spec(_spec('T2|m92|x1100_y220|default'))
    assert calls == [('run_path', ('m2patch', 'm92', 1100, 220), 5.0, 1.2, ())]
    exp = _expected(0, 40)
    assert sorted(set(arr) - {'det', 'rsd'}) == sorted(exp) and all(np.array_equal(arr[k], exp[k], equal_nan=True) for k in exp)
    assert np.array_equal(arr['det'], _det_of(0, 5.0))                                                                       # the catalogue at the reported kappa of the a priori default, 5
    assert arr['rsd'].tolist() == [1.05 + 0.01 * len(arr['det']) + 0.001 * 5.0]                                             # extra field: the isotropic kernel's residual scale of that catalogue (at the high-pass of the run)
    arr = E.run_path_spec(_spec('T|x310_y630|ours'))
    assert calls[-1] == ('run_path', ('m2patch', 'm2', 310, 630), 3.0, 1.2, ()) and np.array_equal(arr['det'], _det_of(0, 2.0))     # and kappa 2 for the tuned configuration


def test_a_variant_run_passes_its_overrides_and_stores_no_catalogue(monkeypatch):
    calls = _stub_world(monkeypatch)
    for name, expect in (('nohp', (0.0, 1.2, ())), ('norex', (3.0, 0.0, ())), ('greedy_snr', (3.0, 1.2, (('order', 'greedy'), ('stop', 'snr')))), ('snr', (3.0, 1.2, (('stop', 'snr'),))),
                         ('nocleanup', (3.0, 1.2, (('cleanup', False),)))):
        arr = E.run_path_spec(_spec('V|x310_y630|%s' % name))
        assert calls[-1] == ('run_path', ('m2patch', 'm2', 310, 630), expect[0], expect[1], expect[2]), name
        assert 'det' not in arr and 'secs' in arr and 'k2|tot' in arr


def test_a_kernel_study_run_stores_the_residual_scale_of_the_frozen_catalogue_and_of_the_new_one(monkeypatch):
    calls = _stub_world(monkeypatch)
    m2 = {'T2|m53|x760_y590|ours|det': np.zeros((7, 3))}
    monkeypatch.setattr(E, '_m2raw', lambda: m2)
    arr = E.run_path_spec(_spec('P|m53|x760_y590|ours'))
    assert calls == [('run_path_psf', ('m2patch', 'm53', 760, 590), 3.0, 1.2, 'general')]
    ks = X.KAPPAS + (0.0,)
    assert arr['rsd'].tolist() == [1.05 + 0.01 * 7 + 0.001 * 3.0, 1.1 + 0.1 * ks.index(2.0)]                                               # [frozen catalogue with the isotropic kernel, the new run at the reported kappa]
    assert np.array_equal(arr['det'], _det_of(1, 2.0)) and arr['bic'].tolist() == [4, 9, 0] and np.array_equal(arr['k5|tot'], _expected(1, 41)['k5|tot'])


def test_a_repeat_run_uses_the_kernel_the_lengths_and_the_residual_scale_of_the_frozen_script(monkeypatch):
    calls = _stub_world(monkeypatch)
    arr = E.run_path_spec(_spec('R|m2r2|x133_y281|iso|default'))
    assert calls[-1] == ('run_path', ('rpatch', 'm2r2', 133, 281), 5.0, 1.2, ()) and arr['rsd'].tolist() == [1.05 + 0.01 * len(arr['det']) + 0.001 * 5.0] and np.array_equal(arr['det'], _det_of(0, 5.0))
    arr = E.run_path_spec(_spec('R|m2r2|x133_y281|adq|ours'))
    assert calls[-1] == ('run_path_psf', ('rpatch', 'm2r2', 133, 281), 3.0, 1.2, 'general') and arr['rsd'].tolist() == [1.1 + 0.1 * (X.KAPPAS + (0.0,)).index(2.0)] and np.array_equal(arr['det'], _det_of(1, 2.0))
    arr = E.run_path_spec(_spec('R|m2r2|x133_y281|scl|default'))
    assert calls[-1] == ('run_path_psf', ('rpatch', 'm2r2', 133, 281), 5.0 * 1.3, 1.2 * 1.3, 'general') and arr['rsd'].tolist() == [1.1 + 0.1 * (X.KAPPAS + (0.0,)).index(5.0)] and np.array_equal(arr['det'], _det_of(1, 5.0))


def test_store_run_flags_a_path_that_reached_the_cap_and_stores_no_catalogue_unless_asked():
    dets = _dets(0)
    arr = E.store_run({}, dets, 3.0, X.MAX_ENTRIES, dict(bic_index=7, n_blocks=X.MAX_ENTRIES + 1), TRUTH, 2.0)
    assert arr['bic'].tolist() == [7, X.MAX_ENTRIES + 1, 1] and arr['secs'].tolist() == [3.0, float(X.MAX_ENTRIES)] and np.array_equal(arr['det'], _det_of(0, 2.0))
    arr = E.store_run({}, dets, 3.0, X.MAX_ENTRIES - 1, dict(bic_index=7, n_blocks=X.MAX_ENTRIES), TRUTH, 2.0, det=False, bic=False)
    assert 'det' not in arr and 'bic' not in arr and len([k for k in arr if k.endswith('|tot')]) == len(X.KAPPAS) + 1


def test_an_unknown_family_is_refused():
    with pytest.raises(ValueError):
        E.run_path_spec(dict(rid='x', fam='Z', pxy=(0, 0)))


def test_a_job_installs_the_guard_runs_the_path_once_and_returns_the_arrays_with_the_trace(monkeypatch, capsys):
    monkeypatch.setattr(IM, 'sliding_path', IM.sliding_path)

    def fake_spec(spec):
        psf = IM.PixelPSF(20, 20, 1.3)
        y = psf.synth(np.array([[8.0, 9.0]]), np.array([900.0])) + np.random.default_rng(0).normal(0, 1, (20, 20))
        IM.sliding_path(psf, y, np.ones((20, 20)), IM.Lattice(psf, h=0.5, margin=1.0), max_entries=3, tol=1e-8, passes=1, rtrunc=4.0)
        return {'k2|tot': np.arange(7.0), 'secs': np.array([1.0, 3.0])}
    monkeypatch.setattr(E, 'run_path_spec', fake_spec)
    E.TRACE.append(dict(counts=np.zeros(8), barred=np.zeros((0, 3)), taus=np.zeros(0)))                                     # a leftover record of an earlier path
    prefix, arrays, wall = E.run_job((dict(rid='D|x1_y2|hp3|rex1.2', fam='D'), 'G'))
    assert prefix == 'G|D|x1_y2|hp3|rex1.2' and arrays['trace|counts'][0] == 3 and arrays['trace|counts'][6] == 0 and arrays['trace|barred'].shape == (0, 3)
    assert 'trace|taus' in arrays and arrays['trace|taus'].shape == (3,) and arrays['wall'][0] == wall > 0                    # stage-1 runs keep their entry times
    prefix, arrays, wall = E.run_job((dict(rid='D|x1_y2|hp3|rex1.2', fam='D'), 'C'))
    assert prefix == 'C|D|x1_y2|hp3|rex1.2' and 'trace|taus' not in arrays                                                    # a check run without activation or cap does not
    assert 'start' in capsys.readouterr().out
    monkeypatch.setattr(E, 'run_path_spec', lambda spec: (IM.sliding_path(IM.PixelPSF(8, 8, 1.3), np.zeros((8, 8)), np.ones((8, 8)), IM.Lattice(IM.PixelPSF(8, 8, 1.3)), max_entries=1), IM.sliding_path(
        IM.PixelPSF(8, 8, 1.3), np.zeros((8, 8)), np.ones((8, 8)), IM.Lattice(IM.PixelPSF(8, 8, 1.3)), max_entries=1), {})[2])
    with pytest.raises(RuntimeError, match='2 paths'):
        E.run_job((dict(rid='x', fam='D'), 'G'))
    E.uninstall()


def test_a_check_run_keeps_its_entry_times_when_the_guard_acted_or_the_cap_was_reached(monkeypatch):
    monkeypatch.setattr(IM, 'sliding_path', IM.sliding_path)
    monkeypatch.setattr(E, 'run_path_spec', lambda spec: (E.TRACE.append(dict(counts=np.array([5, 0, 5, 1, 0, 6, 1, 0]), barred=np.array([[0.1, 1, 2]]), taus=np.arange(5.0))), {'secs': np.array([1.0, 5.0])})[1])
    monkeypatch.setattr(E, 'install', lambda: None)
    _, arrays, _ = E.run_job((dict(rid='r', fam='D'), 'C'))
    assert arrays['trace|taus'].tolist() == [0.0, 1.0, 2.0, 3.0, 4.0] and arrays['trace|barred'].tolist() == [[0.1, 1, 2]]
    monkeypatch.setattr(E, 'run_path_spec', lambda spec: (E.TRACE.append(dict(counts=np.array([X.MAX_ENTRIES, 0, 5, 1, 0, 6, 0, 0]), barred=np.zeros((0, 3)), taus=np.arange(3.0))), {'secs': np.array([1.0, 5.0])})[1])
    _, arrays, _ = E.run_job((dict(rid='r', fam='D'), 'C'))
    assert 'trace|taus' in arrays


# ---- comparison with the stored run, overlay ---------------------------------------------------------------------------------------------

def test_compare_run_demands_equality_of_the_entry_count_the_bic_record_and_every_score_array():
    old = _run('good', 1300, det='good', rsd=[1.0, 2.0])
    new = {k: v.copy() for k, v in old.items()}
    new['secs'] = np.array([999.0, 1300.0])                                                                              # another run time is not a difference
    new['trace|counts'] = np.zeros(8)
    new['wall'] = np.array([1.0])
    new['extra_field'] = np.zeros(2)
    ok, diff, missing, extra = E.compare_run(new, old)
    assert ok and diff == [] and missing == [] and extra == ['extra_field']                                                # trace and wall are not "extra"; a field the frozen script did not store is
    for f, change in (('secs', lambda v: v + np.array([0.0, 1.0])), ('bic', lambda v: v + np.array([0, 0, 1])), ('k2|tot', lambda v: v + 1), ('k0|blk22', lambda v: v + 1.0), ('det', lambda v: v + 1e-12),
                      ('rsd', lambda v: v * 1.0000001), ('k5|offset', lambda v: v + 0.001)):
        bad = {k: v.copy() for k, v in old.items()}
        bad[f] = change(bad[f])
        ok, diff, missing, extra = E.compare_run(bad, old)
        assert not ok and diff == [f], f
    gone = {k: v for k, v in old.items() if k != 'det'}
    ok, diff, missing, extra = E.compare_run(gone, old)
    assert not ok and missing == ['det'] and diff == []
    shape = dict(new, det=new['det'][:-1])
    assert E.compare_run(shape, old)[1] == ['det']
    nan_old = {'k2|offset': np.array([np.nan, 1.0])}
    assert E.compare_run({'k2|offset': np.array([np.nan, 1.0])}, nan_old)[0]                                                 # NaN equals NaN
    assert not E.compare_run({'k2|offset': np.array([np.nan, 2.0])}, nan_old)[0]
    assert E.compare_run({'name': np.array(['a'])}, {'name': np.array(['a'])})[0] and not E.compare_run({'name': np.array(['b'])}, {'name': np.array(['a'])})[0]


def test_stored_fields_and_overlay_replace_only_the_namespace_asked_for():
    base = {'D|a|k2|tot': np.arange(3.0), 'D|a|secs': np.array([1.0, 2.0]), 'D|b|k2|tot': np.zeros(3), 'D|ab|k2|tot': np.ones(3)}
    assert sorted(E.stored_fields(base, 'D|a')) == ['k2|tot', 'secs']                                                      # the prefix is the run, not a longer name
    cf = {'G|D|a|k2|tot': np.full(3, 9.0), 'G|D|a|trace|counts': np.zeros(8), 'G|D|a|wall': np.array([3.0]), 'C|D|b|k2|tot': np.full(3, 7.0)}
    out = E.overlay(base, cf, 'G')
    assert out['D|a|k2|tot'].tolist() == [9.0] * 3 and out['D|b|k2|tot'].tolist() == [0.0] * 3 and 'D|a|trace|counts' not in out and 'D|a|wall' not in out and out['D|a|secs'].tolist() == [1.0, 2.0]
    assert base['D|a|k2|tot'].tolist() == [0.0, 1.0, 2.0]                                                                  # the frozen dict is not modified
    assert E.overlay(base, cf, 'C')['D|b|k2|tot'].tolist() == [7.0] * 3


def test_alias_repeat_copies_the_isotropic_runs_over_the_chosen_kernel_only_where_the_rule_keeps_it():
    fro = world()
    rep = fro['rep']
    rep['R|m2r3|x110_y501|iso|default|k5|tot'] = np.full(7, 5.0)                                                            # a replaced isotropic run
    rep['R|m2r2|x133_y281|iso|default|k5|tot'] = np.full(7, 6.0)
    out = E.alias_repeat(rep)
    assert out['R|m2r3|x110_y501|adq|default|k5|tot'].tolist() == [5.0] * 7                                                 # the rule keeps the isotropic kernel: the same run
    assert out['R|m2r2|x133_y281|adq|default|k5|tot'].tolist() != [6.0] * 7                                                 # a general kernel there: its own run
    assert out['R|m2r1|x316_y124|adq|ours|k2|tot'].tolist() == rep['R|m2r1|x316_y124|iso|ours|k2|tot'].tolist()
    assert rep['R|m2r3|x110_y501|adq|default|k5|tot'].tolist() != [5.0] * 7                                                 # the input is not modified


# ---- raw files, write-once running, the freeze ------------------------------------------------------------------------------------------

def test_two_processes_saving_into_one_raw_lose_nothing(tmp_path):
    p = str(tmp_path / 'raw.npz')
    code = ("import sys, numpy as np; sys.path.insert(0, %r); sys.path.insert(0, %r); import exp_astro_cyclefix as E\n"
            "rng = np.random.default_rng(0)\nfor i in range(25):\n    E.save_part(%r, {'%%s%%d' %% (sys.argv[1], i): rng.normal(size=60000)})\n") % (os.path.join(HERE, '..', 'experiments'), os.path.join(HERE, '..'), p)
    procs = [subprocess.Popen([sys.executable, '-c', code, tag], env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1')) for tag in ('a', 'b')]
    assert [pr.wait() for pr in procs] == [0, 0]
    assert sorted(E.load_raw(p)) == sorted(['a%d' % i for i in range(25)] + ['b%d' % i for i in range(25)])                     # every save of both processes is there


def test_an_array_with_nan_saved_again_is_the_same_array(tmp_path):
    p = str(tmp_path / 'raw.npz')
    a = np.array([np.nan, 1.0, 2.0])
    E.save_part(p, {'k2|completeness': a})
    E.save_part(p, {'k2|completeness': a.copy(), 'other': np.ones(2)})                                                      # a path saved twice with the same scores (NaN where a bin is empty) is not a clash
    assert E.same_array(a, a.copy()) and not E.same_array(a, np.array([np.nan, 1.0, 3.0])) and not E.same_array(a, a[:2]) and E.same_array(np.array(['x']), np.array(['x'])) and not E.same_array(np.array(['x']), np.array(['y']))
    with pytest.raises(RuntimeError, match='write-once'):
        E.save_part(p, {'k2|completeness': np.array([np.nan, 1.0, 9.0])})


def test_save_part_merges_atomically_and_refuses_to_change_a_stored_value(tmp_path):
    p = str(tmp_path / 'raw.npz')
    E.save_part(p, {'a': np.arange(3.0)})
    E.save_part(p, {'b': np.ones(2), 'a': np.arange(3.0)})                                                                 # the same value again is fine
    assert sorted(E.load_raw(p)) == ['a', 'b'] and sorted(os.listdir(tmp_path)) == ['raw.npz', 'raw.npz.lock']                  # nothing left beside the raw but the lock file
    with pytest.raises(RuntimeError, match='write-once'):
        E.save_part(p, {'a': np.zeros(3)})
    assert E.load_raw(p)['a'].tolist() == [0.0, 1.0, 2.0]


_CALLS = []


def _stub_job(job):
    """A module-level stand-in for `run_job` (a pool pickles its function by name): records the call and returns arrays like a job does."""
    spec, ns = job
    _CALLS.append((spec['rid'], ns))
    return ns + '|' + spec['rid'], {'secs': np.array([1.0, 2.0]), 'wall': np.array([1.0]), 'x': np.array([float(len(spec['rid']))])}, 1.0


@pytest.mark.parametrize('workers', [1, 2])
def test_run_jobs_runs_each_missing_path_once_longest_first_and_never_overwrites(tmp_path, monkeypatch, workers):
    specs = [dict(rid='A|1', cost=10.0), dict(rid='A|22', cost=300.0), dict(rid='A|333', cost=50.0)]
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: ({}, {'meta|sha256': np.array(['h']), 'meta|files': np.array(['f']), 'meta|plan_sha256': np.array(['p']), 'meta|versions': np.array(['v'])}))
    monkeypatch.setattr(E, 'run_job', _stub_job)
    out = str(tmp_path / 'raw.npz')
    _CALLS.clear()
    assert E.run_jobs(specs, 'G', out, workers=workers) == 3
    if workers == 1:
        assert [c[0] for c in _CALLS] == ['A|22', 'A|333', 'A|1']                                                          # the longest stored run first
    raw = E.load_raw(out)
    assert sorted(k for k in raw if k.endswith('|wall')) == ['G|A|1|wall', 'G|A|22|wall', 'G|A|333|wall'] and raw['G|A|22|x'][0] == 4.0 and raw['meta|sha256'][0] == 'h'
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: (E.load_raw(out), {}))
    _CALLS.clear()
    assert E.run_jobs(specs + [dict(rid='B|5', cost=5.0)], 'R', out, workers=1, reverse=True) == 4 and [c[0] for c in _CALLS] == ['B|5', 'A|1', 'A|333', 'A|22']      # the shortest first
    _CALLS.clear()
    assert E.run_jobs(specs, 'G', out, workers=1) == 0 and _CALLS == []                                                    # nothing to do; the record of the first run is kept
    assert E.load_raw(out)['meta|sha256'][0] == 'h'
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: (E.load_raw(out), {'meta|files_late': np.array(['x.json']), 'meta|sha256_late': np.array(['hx'])}))
    assert E.run_jobs(specs + [dict(rid='A|4444', cost=1.0)], 'C', out, workers=1) == 4                                      # the other namespace is a different run
    assert sorted(c[1] for c in _CALLS) == ['C'] * 4
    raw = E.load_raw(out)
    assert raw['meta|files_late'].tolist() == ['x.json'] and raw['meta|sha256_late'].tolist() == ['hx'] and raw['meta|sha256'][0] == 'h'          # a file added to the record later is recorded once


def test_run_jobs_keeps_the_result_of_another_process_that_finished_the_same_path_first(tmp_path, monkeypatch, capsys):
    out = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: ({}, {}))

    def rival_job(job):
        spec, ns = job
        E.save_part(out, {'%s|%s|wall' % (ns, spec['rid']): np.array([99.0]), '%s|%s|x' % (ns, spec['rid']): np.array([7.0])})     # the other process saves first
        return ns + '|' + spec['rid'], {'wall': np.array([1.0]), 'x': np.array([1.0])}, 1.0
    monkeypatch.setattr(E, 'run_job', rival_job)
    assert E.run_jobs([dict(rid='A|1', cost=1.0)], 'C', out, workers=1) == 1
    raw = E.load_raw(out)
    assert raw['C|A|1|x'].tolist() == [7.0] and raw['C|A|1|wall'].tolist() == [99.0]                                          # not overwritten, and no clash error
    assert 'already saved by another process' in capsys.readouterr().out


def test_run_jobs_writes_the_record_of_the_frozen_state_with_the_first_result_only(tmp_path, monkeypatch):
    out = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: ({}, {'meta|sha256': np.array(['h']), 'meta|files': np.array(['f'])}))
    monkeypatch.setattr(E, 'run_job', _stub_job)
    calls, real = [], E.save_part
    monkeypatch.setattr(E, 'save_part', lambda path, arrays: (calls.append(sorted(k for k in arrays if k.startswith('meta|'))), real(path, arrays))[1])
    _CALLS.clear()
    assert E.run_jobs([dict(rid='A|1', cost=3.0), dict(rid='A|2', cost=2.0), dict(rid='A|3', cost=1.0)], 'G', out, workers=1) == 3
    assert calls == [['meta|files', 'meta|sha256'], [], []]                                                                # the record goes in with the first result, once


def test_run_jobs_does_not_start_a_path_that_another_process_saved_meanwhile(tmp_path, monkeypatch, capsys):
    out = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: ({}, {}))
    ran = []

    def job_with_rival(job):
        spec, ns = job
        ran.append(spec['rid'])
        if spec['rid'] == 'A|1':                                                                                           # while this path runs, the other process saves the next one
            E.save_part(out, {'%s|A|2|wall' % ns: np.array([99.0]), '%s|A|2|x' % ns: np.array([7.0])})
        return ns + '|' + spec['rid'], {'wall': np.array([1.0]), 'x': np.array([1.0])}, 1.0
    monkeypatch.setattr(E, 'run_job', job_with_rival)
    assert E.run_jobs([dict(rid='A|1', cost=9.0), dict(rid='A|2', cost=5.0), dict(rid='A|3', cost=1.0)], 'C', out, workers=1) == 3
    assert ran == ['A|1', 'A|3']                                                                                           # A|2 is not run again
    raw = E.load_raw(out)
    assert raw['C|A|2|x'].tolist() == [7.0] and raw['C|A|3|x'].tolist() == [1.0]
    assert 'saved meanwhile by another process, not run' in capsys.readouterr().out


def test_preflight_refuses_to_start_when_a_freeze_record_a_hashed_file_or_the_plan_differs(tmp_path, monkeypatch):
    plan = tmp_path / 'plan.md'
    plan.write_text('the plan')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    monkeypatch.setattr(E, 'hashes', lambda files=E.HASHED: {'a.py': 'h1', 'b.py': 'h2'})
    monkeypatch.setattr(E, 'load_raw', lambda p: {'m2': 1} if p == E.M2RAW else np.load(p, allow_pickle=False))
    checks = []
    monkeypatch.setattr(X, 'check_freeze', lambda raw, path=None: checks.append('m2'))
    monkeypatch.setattr(RE, 'check_freeze', lambda path=None: checks.append('repeat'))
    out = str(tmp_path / 'raw.npz')
    cur, now = E.preflight(out)
    assert cur == {} and checks == ['m2', 'repeat'] and now['meta|files'].tolist() == ['a.py', 'b.py'] and now['meta|plan_sha256'][0] == hashlib.sha256(b'the plan').hexdigest()
    assert now['meta|versions'][1].startswith('numpy ') and sorted(now) == ['meta|files', 'meta|plan_sha256', 'meta|sha256', 'meta|versions']
    np.savez_compressed(out, **now)
    cur, again = E.preflight(out)                                                                                          # the same state: fine, nothing new to record
    assert again == {} and 'meta|sha256' in cur
    monkeypatch.setattr(E, 'hashes', lambda files=E.HASHED: {'a.py': 'h1', 'b.py': 'CHANGED'})
    with pytest.raises(RuntimeError, match='b.py'):
        E.preflight(out)
    monkeypatch.setattr(E, 'hashes', lambda files=E.HASHED: {'a.py': 'h1', 'b.py': 'h2', 'c.py': 'h3'})
    cur, late = E.preflight(out)                                                                                           # a file added to the list later: recorded once, then checked like the others
    assert late['meta|files_late'].tolist() == ['c.py'] and late['meta|sha256_late'].tolist() == ['h3']
    np.savez_compressed(out, **dict(cur, **late))
    assert E.preflight(out)[1] == {}
    monkeypatch.setattr(E, 'hashes', lambda files=E.HASHED: {'a.py': 'h1', 'b.py': 'h2', 'c.py': 'other'})
    with pytest.raises(RuntimeError, match='c.py'):
        E.preflight(out)
    monkeypatch.setattr(E, 'hashes', lambda files=E.HASHED: {'a.py': 'h1', 'b.py': 'h2', 'c.py': 'h3'})
    plan.write_text('another plan')
    with pytest.raises(RuntimeError, match='the plan'):
        E.preflight(out)
    plan.write_text('the plan')

    def refuse(raw, path=None):
        raise RuntimeError('the freeze record differs')
    monkeypatch.setattr(X, 'check_freeze', refuse)
    with pytest.raises(RuntimeError, match='freeze'):
        E.preflight(out)
    E.preflight(out, freeze_checks=False)


def test_stage_one_runs_the_capped_paths_and_check_mode_the_others_in_the_declared_order(tmp_path, monkeypatch):
    fro = world()
    monkeypatch.setattr(E, 'frozen_raws', lambda: fro)
    monkeypatch.setattr(E, 'frozen_config', lambda: CFG)
    monkeypatch.setattr(E, 'preflight', lambda out_path=None, freeze_checks=True: ({}, {'meta|sha256': np.array(['h'])}))
    monkeypatch.setattr(E, 'run_job', _stub_job)
    out = str(tmp_path / 'raw.npz')
    _CALLS.clear()
    assert E.stage1(1, out) == 12 and {c[1] for c in _CALLS} == {'G'} and sorted(c[0] for c in _CALLS) == sorted(STAGE1)
    _CALLS.clear()
    assert E.check_mode(1, out_path=out) == 65 and {c[1] for c in _CALLS} == {'C'}
    fam = [c[0].split('|')[0] for c in _CALLS]
    first = {f: fam.index(f) for f in ('T', 'T2', 'D', 'P', 'R', 'V')}
    assert max(first['T'], first['T2']) < first['D'] < first['P'] < first['R'] < first['V'] or fam[:18].count('D') == 0                  # tiers: test patches, then development paths, kernel study, repeat patches, variants
    assert sorted(set(fam[:18])) == ['T', 'T2'] and set(fam[18:24]) == {'D'} and set(fam[24:36]) == {'P'} and set(fam[36:61]) == {'R'} and set(fam[61:]) == {'V'}
    _CALLS.clear()
    assert E.check_mode(1, only=['V'], out_path=str(tmp_path / 'other.npz')) == 4 and all(c[0].startswith('V|') for c in _CALLS)
    order = [c[0] for c in _CALLS]
    assert order == ['V|x310_y630|nocleanup', 'V|x310_y630|snr', 'V|x310_y630|nohp', 'V|x310_y630|greedy_snr']              # the longest stored run first
    _CALLS.clear()
    assert E.check_mode(1, only=['V'], out_path=str(tmp_path / 'third.npz'), reverse=True) == 4
    assert [c[0] for c in _CALLS] == order[::-1]                                                                           # and the shortest first for a second process


# ---- the report -------------------------------------------------------------------------------------------------------------------------

def _report(cf, fro=None):
    return E.report(cf, fro or world())


def test_the_report_is_deterministic_and_carries_the_labels_and_the_sections():
    fro = world()
    cf = checked(fro, guarded(fro))
    md = E.report(cf, fro)
    assert md == E.report(cf, fro)
    heads = [md.index(h) for h in ('## 1. ', '## 2. ', '## 3. ', '## 4. ', '### 4b. ', '## 5. ', '### 5.1 ', '### 5.2 ', '### 5.3 ', '### 5.4 ', '### 5.5 ', '## 6. ', '## 7. ')]
    assert heads == sorted(heads)                                                                       # the sections come in the declared order
    fam = [l.split('|')[1].strip() for l in md.split('## 6. ')[1].split('## 7. ')[0].split('\n') if l.startswith('| ') and l.split('|')[1].strip() in 'TDPRV' and len(l.split('|')[1].strip()) == 1]
    assert fam == sorted(fam, key=lambda f: E.CHECK_ORDER.index(f)) and set(fam) == set('TDPRV')           # the check table lists the families in the order of check mode
    for needle in ('# The astronomy paths that an enter/exit limit cycle', E.LABEL, E.ASRUN, '## 1. The defect and the guard', '## 2. Frozen state', '## 3. The affected paths: as run and guarded, side by side',
                   '## 4. The development table recomputed', '## 5. The corrected rows that the paper uses', '### 5.1 M92 (NGC 6341)', '### 5.2 The repeat frames', '### 5.3 P1, the secondary lines',
                   '### 5.4 Added after the plan', '### 5.5 The kernel study', '## 6. Check mode', '## 7. Limitations', 'Frozen configuration: sigma_b 3, r_ex 1.2, kappa 2; a priori default 5, 1.2, 5.', 'sha256 ' + '1' * 64):
        assert needle in md, needle
    assert not re.search(r'(\|\s*nan\s*\||\bnan \[|\(nan,|\[nan,)', md.lower())


def test_the_affected_table_shows_the_stored_and_the_guarded_run_side_by_side_with_the_activations():
    fro = world()
    rid = 'D|x110_y630|hp5|rex0.8'
    cf = guarded(fro, act=0)
    cf.update({'G|%s|%s' % (rid, f): v for f, v in _trace(1700, act=2, bars=[[2.19, 110, 54], [2.41, 179, 88]], most=3).items()})
    cf['G|V|x310_y630|norex|trace|counts'] = np.array([X.MAX_ENTRIES, 4, 1700, 7, 2, X.MAX_ENTRIES + 1, 0, 1])
    md = E.report(cf, fro)
    sec = md.split('## 3. ')[1].split('## 4. ')[0]
    row = [l for l in sec.split('\n') if l.startswith('| development patch (110, 630), sigma_b 5, r_ex 0.8 `%s`' % rid.replace('|', '\\|'))][0]
    c = cells(row)
    assert c[1:] == ['2500', '1520 of 2501', '300', '1700', '1670 of 1701', '100', '2', '0', '1700', '3', 'ended before the cap']
    assert c[0].endswith('`%s`' % rid)
    assert '* `%s`: guard activations (entry time, lattice point (b, a)): 2.19 (110, 54); 2.41 (179, 88).' % rid.replace('|', '\\|') in sec
    norex = [l for l in sec.split('\n') if 'V\\|x310_y630\\|norex' in l and l.startswith('| design variant')][0]
    assert 'CAP REACHED, no activation: censoring of a long path; possible residual cycling (a point entered 7 times)' in norex and 'added by the scan of the entry counts' in norex
    assert cells(norex)[1:4] == ['2500', '-', '%.0f' % (200.0 + 10.0 * list(X.VARIANTS).index('norex'))]                                                                       # (the stored variant has no BIC block)
    assert sec.count('not re-run') == 0
    assert 'Cap reached again' not in sec                                                                                    # (only for a run that reached the cap: none of these did)


def test_the_affected_section_compares_the_guarded_runs_with_the_check_runs_when_there_are_both():
    fro = world()
    cf = guarded(fro)
    sec = E.report(cf, fro).split('## 3. ')[1].split('## 4. ')[0]
    assert 'For comparison' not in sec                                                                                  # no check-mode run yet
    cf = checked(fro, cf, rids={'T|x310_y630|ours', 'T|x310_y630|default'})
    cf['C|T|x310_y630|default|trace|counts'] = np.array([1000, 80, 900, 4, 0, 1001, 0, 0])
    cf['C|T|x310_y630|ours|trace|counts'] = np.array([1300, 26, 1300, 2, 0, 1301, 0, 0])
    sec = E.report(cf, fro).split('## 3. ')[1].split('## 4. ')[0]
    assert 'the 2 check-mode paths run so far (most entries of one point: median 3, maximum 4; exit events per entry: median 0.050, maximum 0.080)' in sec
    assert '(most entries of one point: median 1, maximum 1; exit events per entry: median 0.000, maximum 0.000)' in sec


def test_a_guarded_run_that_reaches_the_cap_with_activations_is_said_to_and_shows_its_entry_times():
    fro = world()
    cf = guarded(fro, cap=True, act=1)
    sec = E.report(cf, fro).split('## 3. ')[1].split('## 4. ')[0]
    q = np.percentile(np.linspace(0.01, 1.0, X.MAX_ENTRIES), [0, 25, 50, 75, 100])
    assert 'CAP REACHED with 1 activations' in sec
    assert 'Cap reached again. Entry times (of 2500 entries): min %.5g, quartiles %.5g, %.5g, %.5g, max %.5g.' % tuple(q[[0, 1, 2, 3, 4]]) in sec


def test_runs_not_yet_done_are_listed_as_such_and_no_outcome_is_declared():
    fro = world()
    cf = guarded(fro, rids=['D|x110_y630|hp5|rex0.8', 'R|m2r2|x133_y281|iso|default'])
    md = E.report(cf, fro)
    assert md.count('| not re-run |') == 10 and 'INCOMPLETE: 5 of the 6 capped development paths are not re-run yet; no outcome is declared.' in md
    assert 'Not all M92 default paths are re-run yet.' in md and 'Not all three repeat default paths are re-run yet.' in md and 'The isotropic default path of m2r2 (133, 281) is not re-run yet.' not in md
    assert md.count('not run') >= 65 and 'Summary: 0 paths identical, 0 different, 65 not run (of 65).' in md
    empty = E.report({}, fro)
    assert 'No run yet.' in empty and 'INCOMPLETE: 6 of the 6' in empty


def test_the_development_table_declares_outcome_A_when_the_frozen_configuration_is_still_the_argmax():
    fro = world()
    cf = guarded(fro, key='poor')                                                                                          # the guarded rows are worse than the frozen point's rows
    sec = E.report(cf, fro).split('## 4. ')[1].split('## 5. ')[0]
    assert '**Outcome (A)**' in sec and 'Outcome (B)' not in sec
    fstored = X.choose(fro['m2'])[3]
    assert 'Corrected pooled F of the frozen point: %.4f (stored: %.4f).' % (fstored, fstored) in sec
    h, r, k, f = X.choose(overlay_m2(fro, cf))
    assert (h, r, k) == (3.0, 1.2, 2.0)
    frozen_row = [cells(l) for l in sec.split('\n') if l.startswith('| 3 | 1.2 | 2 |')][0]
    assert frozen_row[5] == '**%.3f**' % f and frozen_row[9] == '**%.3f**' % f and frozen_row[7:9] == ['same', 'same']      # the argmax is marked in both blocks, though the path was not re-run
    assert '* Argmax of the stored table: sigma_b 3, r_ex 1.2, kappa 2, pooled F %.4f (the frozen configuration).' % X.choose(fro['m2'])[3] in sec
    assert '* Argmax of the corrected table: sigma_b 3, r_ex 1.2, kappa 2, pooled F %.4f.' % f in sec
    assert 'Robustness line' in sec and '(differs from the argmax' not in sec
    top = [cells(l) for l in sec.split('The ten best rows')[1].split('\n') if l.startswith('| ') and l.split('|')[1].strip().isdigit()]
    assert len(top) == 10 and top[0][1:4] == ['3', '1.2', '2'] and top[0][4] == '%.4f' % f and [float(r[4]) for r in top] == sorted((float(r[4]) for r in top), reverse=True)
    assert all(r[5] == '-' for r in top[:1])


def overlay_m2(fro, cf):
    return E.overlay(fro['m2'], cf, 'G')


def test_the_development_table_declares_outcome_B_when_another_configuration_wins_and_flags_a_capped_winner():
    fro = world()
    cf = guarded(fro, key='good')                                                                                          # the guarded rows are better than every stored row
    sec = E.report(cf, fro).split('## 4. ')[1].split('## 5. ')[0]
    h, r, k, f = X.choose(overlay_m2(fro, cf))
    assert (h, r, k) != (3.0, 1.2, 2.0)
    assert '**Outcome (B)**' in sec and 'sigma_b %g, r_ex %g, kappa %g with pooled F %.4f' % (h, r, k, f) in sec and 'on the 16 test patches' in sec and 'Outcome (A)' not in sec
    assert '| %g | %g | %g |' % (h, r, k) in sec and '**%.3f**' % f in sec
    fz = X.f_measure(sum(fro['m2']['D|x%d_y%d|hp3|rex1.2|k2|tot' % p] for p in X.DEV))[0]
    assert 'against %.4f for the frozen point (stored: %.4f)' % (fz, fz) in sec
    row = [cells(l) for l in sec.split('\n') if l.startswith('| 5 | 0.8 | 2 |')][0]                                          # a re-run row: the stored and the guarded block side by side
    tot_old = sum(fro['m2']['D|x%d_y%d|hp5|rex0.8|k2|tot' % p] for p in X.DEV)
    tot_new = sum(cf['G|D|x%d_y%d|hp5|rex0.8|k2|tot' % p] for p in X.DEV)
    (f0, c0, p0), (f1, c1, p1) = X.f_measure(tot_old), X.f_measure(tot_new)
    assert row[3:8] == ['%.3f' % c0, '%.3f' % p0, '%.3f' % f0, '2 of 2', '%.3f' % c1] and row[8:10] == ['%.3f' % p1, '%.3f' % f1]
    assert [cells(l) for l in sec.split('\n') if l.startswith('| 3 | 0.8 | 2 |')][0][7:10] == ['same', 'same', 'same']                 # a path that was not re-run is shown as it was
    cf = guarded(fro, key='good', cap=True)                                                                                # capped again: the winner's row is flagged, and the robustness line excludes it
    sec = E.report(cf, fro).split('## 4. ')[1].split('## 5. ')[0]
    assert '| 2 of 2 |' in sec and '(differs from the argmax of the full corrected table)' in sec


DEV_UNCAPPED = ['D|x110_y630|hp3|rex0.8', 'D|x110_y630|hp3|rex1.2', 'D|x110_y630|hp5|rex1.2', 'D|x210_y630|hp3|rex0.8', 'D|x210_y630|hp3|rex1.2', 'D|x210_y630|hp5|rex1.2']


def test_the_fully_guarded_development_table_needs_all_twelve_paths_and_keeps_the_frozen_argmax_when_check_runs_are_identical():
    fro = world()
    cf = guarded(fro, key='poor')
    sec = E.report(cf, fro).split('### 4b. ')[1].split('## 5. ')[0]
    assert 'added after the plan' in sec and 'INCOMPLETE: 6 of the 12 development paths have a guarded run so far.' in sec
    half = checked(fro, cf, rids=set(DEV_UNCAPPED[:3]))
    assert 'INCOMPLETE: 9 of the 12' in E.report(half, fro).split('### 4b. ')[1].split('## 5. ')[0]
    cf = checked(fro, cf, rids=set(DEV_UNCAPPED), act={'D|x110_y630|hp3|rex1.2': 2, 'D|x210_y630|hp5|rex1.2': 1})
    sec = E.report(cf, fro).split('### 4b. ')[1].split('## 5. ')[0]
    assert '**the frozen configuration is still the argmax** (rank 1; the second row has F' in sec and 'not the frozen configuration' not in sec
    allg = E.overlay(E.overlay(fro['m2'], cf, 'G'), {k: v for k, v in cf.items() if k.startswith('C|D|')}, 'C')
    h, r, k, f = X.choose(allg)
    assert (h, r, k) == (3.0, 1.2, 2.0) and 'pooled F %.4f' % f in sec
    top = [cells(l) for l in sec.split('The ten best rows')[1].split('\n') if l.startswith('| ') and l.split('|')[1].strip().isdigit()]
    assert len(top) == 10 and top[0][1:4] == ['3', '1.2', '2']
    assert 'Guard activations in the twelve development paths:' in sec and '`D\\|x110_y630\\|hp3\\|rex1.2` 2' in sec and '`D\\|x210_y630\\|hp5\\|rex1.2` 1' in sec and '`D\\|x110_y630\\|hp5\\|rex0.8` 0' in sec


def test_the_fully_guarded_development_table_names_the_new_argmax_when_a_check_run_beats_the_frozen_point():
    fro = world()
    cf = guarded(fro, key='poor')
    cf = checked(fro, cf, rids=set(DEV_UNCAPPED))
    for p in ((110, 630), (210, 630)):                                                                                       # the unaffected path hp3 / r_ex 0.8 is better at every readout when guarded
        rid = 'C|D|x%d_y%d|hp3|rex0.8' % p
        for f, v in _run('good', 1500).items():
            cf[rid + '|' + f] = v
    sec = E.report(cf, fro).split('### 4b. ')[1].split('## 5. ')[0]
    allg = E.overlay(E.overlay(fro['m2'], cf, 'G'), {k: v for k, v in cf.items() if k.startswith('C|D|')}, 'C')
    h, r, k, f = X.choose(allg)
    assert (h, r) == (3.0, 0.8) and k == 8.0
    assert '**not the frozen configuration**, which has rank' in sec and 'Argmax of this table: sigma_b 3, r_ex 0.8, kappa 8, pooled F %.4f' % f in sec and 'still the argmax' not in sec
    assert 'Outcome (B)' not in sec and 'section 4' in sec                                                                  # the plan's outcome is decided by the table of section 4
    rank = int(sec.split('which has rank ')[1].split(' ')[0])
    assert rank == 8                                                                                                        # the seven rows of the better path at kappa 2.5 to 8 rank above the frozen point


def test_the_m92_section_prints_the_stored_and_the_corrected_report_and_a_side_by_side_table():
    fro = world()
    cf = guarded(fro, key='bic')
    md = E.report(cf, fro)
    sec = md.split('### 5.1 ')[1].split('### 5.2 ')[0]
    assert '#### M92 report as `m2.md` prints it, stored arrays (as run)' in sec and '#### M92 report as `m2.md` prints it, **%s**' % E.LABEL in sec
    assert sec.count('#### M92 (NGC 6341), SDSS run 4682 camcol 6 field 208: test patch (1100, 220)') == 2 and sec.count('all 2 test patches together') == 2
    stored = X.report({k: v for k, v in fro['m2'].items() if k.startswith('T2|m92|')}).split('\n## ', 1)[1]
    first = sec.split('stored arrays (as run)')[1].split('#### M92 report as')[0]
    assert all(line in first for line in ('## ' + stored).split('\n')[2:12] if line.strip() and not line.startswith('## '))   # the stored section is X.report's own text
    corr = {k: v for k, v in E.overlay(fro['m2'], cf, 'G').items() if k.startswith('T2|m92|')}
    assert X.report(corr).split('\n## ', 1)[1].split('\n')[2] in sec
    # the side-by-side table: the pooled default row, as run and guarded, equals what X.report computes for each
    pooled = sec.split('| pooled | ISS a priori default (kappa 5) |')[1].split('\n')[0]
    idx = X.boot_idx_strat([25, 25])
    for raw, at in ((fro['m2'], 0), (E.overlay(fro['m2'], cf, 'G'), 1)):
        blk = np.concatenate([raw['T2|m92|x%d_y%d|default|k5|blk22' % tuple(p)] for p in X.SECOND['m92']['patches']])
        (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk, idx)
        cc = [x.strip() for x in pooled.strip('|').split('|')]
        assert cc[at * 4] == X.fmt_iv(c, cci) and cc[at * 4 + 2] == X.fmt_iv(ff, fci)
    assert 'pooled | default ISS vs crowdsource (default) |' in sec and 'BIC-stop default ISS vs photutils (tuned)' in sec
    assert 'M53 (NGC 5024)' not in sec                                                                                      # only the M92 sections
    for z, at in ((fro['m2'], 0), (E.overlay(fro['m2'], cf, 'G'), 1)):                                                      # the paired rows and the BIC row, as run and guarded, from the frozen functions
        for row_name, a_key, b_key in (('default ISS vs crowdsource (default)', 'default|k5', 'crowdsource (default)'), ('BIC-stop default ISS vs photutils (tuned)', 'default|k0', 'photutils (tuned)')):
            ba = np.concatenate([z['T2|m92|x%d_y%d|%s|blk22' % (p[0], p[1], a_key)] for p in X.SECOND['m92']['patches']])
            bb = np.concatenate([z['T2|m92|x%d_y%d|%s|blk22' % (p[0], p[1], b_key)] for p in X.SECOND['m92']['patches']])
            df = X.boot_diff(ba, bb, idx)[2]
            prow = [cells(l) for l in sec.split('\n') if l.startswith('| pooled | %s |' % row_name)][0]
            assert prow[2 + 2 * at] == X.fmt_iv(df[0], df[1]) and prow[3 + 2 * at] == X.verdict(*df)
        bic = np.concatenate([z['T2|m92|x%d_y%d|default|k0|blk22' % tuple(p)] for p in X.SECOND['m92']['patches']])
        (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(bic, idx)
        brow = [cells(l) for l in sec.split('\n') if l.startswith('| pooled | ISS a priori default, BIC stop |')][0]
        assert brow[2 + 4 * at + 2] == X.fmt_iv(ff, fci)
    for z, title in ((fro['m2'], 'stored arrays (as run)'), (E.overlay(fro['m2'], cf, 'G'), '**%s**' % E.LABEL)):               # the whole M92 report of each version is in its block
        md = X.report({k: v for k, v in z.items() if k.startswith('T2|m92|')})
        body = '## ' + md.split('\n## ', 1)[1]
        body = '\n'.join(('#### ' + l[3:]) if l.startswith('## ') else l for l in body.split('\n')).rstrip('\n')
        block = sec.split('#### M92 report as `m2.md` prints it, %s' % title)[1]
        assert body in block


def test_the_repeat_section_replaces_the_capped_isotropic_runs_and_keeps_the_chosen_kernel_of_the_general_patches():
    fro = world()
    cf = guarded(fro, key='bic')
    md = E.report(cf, fro)
    sec = md.split('### 5.2 ')[1].split('### 5.3 ')[0]
    for needle in ('#### All 6 patches, as run', '#### All 6 patches, **%s**' % E.LABEL, '#### Frame m2r1, as run', '#### Frame m2r2, **', '#### Frame m2r3, as run', '(added after the plan)', '#### Record by patch',
                   'adequate kernel vs isotropic, a priori default', 'isotropic kernel (a priori default) vs photutils (default)', 'ISS, isotropic kernel, a priori default with the BIC stop', 'photutils (default)'):
        assert needle in sec, needle
    as_run, corr = RE.method_table(fro['rep'], RE.done_patches(fro['rep']), E.REPEAT_DEFAULT_ROWS)[0], RE.method_table(E.repeat_corrected(fro, cf), RE.done_patches(fro['rep']), E.REPEAT_DEFAULT_ROWS)[0]
    assert as_run in sec and corr in sec and as_run != corr
    rows = [l for l in sec.split('#### All 6 patches, **')[1].split('####')[0].split('\n') if l.startswith('| ISS, isotropic kernel, a priori default (5, 1.2, kappa 5) |')]
    assert len(rows) == 1
    fix = E.repeat_corrected(fro, cf)
    assert np.array_equal(fix['R|m2r3|x110_y501|adq|default|k5|tot'], fix['R|m2r3|x110_y501|iso|default|k5|tot']) and np.array_equal(fix['R|m2r3|x110_y501|iso|default|k5|tot'], cf['G|R|m2r3|x110_y501|iso|default|k5|tot'])
    assert not np.array_equal(fix['R|m2r2|x133_y281|adq|default|k5|tot'], fix['R|m2r2|x133_y281|iso|default|k5|tot'])
    for fname in ('m2r1', 'm2r2', 'm2r3'):                                                                                  # a frame's tables pool the patches of that frame only
        fl = [np_ for np_ in RE.done_patches(fro['rep']) if np_[0] == fname]
        for title, raw in (('as run', fro['rep']), ('**%s**' % E.LABEL, fix)):
            block = sec.split('#### Frame %s, %s' % (fname, title))[1].split('####')[0]
            assert RE.method_table(raw, fl, E.REPEAT_DEFAULT_ROWS)[0] in block
    # the record by patch: one row per pair and per version
    rec = sec.split('#### Record by patch')[1]
    assert rec.count('adequate kernel vs isotropic, a priori default (as run)') == 1 and rec.count('adequate kernel vs isotropic, a priori default (guarded)') == 1


def test_the_p1_lines_are_the_function_of_the_repeat_report_with_two_and_with_forty_thousand_resamples(monkeypatch):
    fro = world()
    plist = RE.done_patches(fro['rep'])
    gen, a2000 = E.p1_lines(fro['rep'], plist)
    assert gen == [('m2r2', (133, 281)), ('m2r2', (195, 689))] and [l[0] for l in a2000] == ['tuned configuration (kappa 2)', 'a priori default (kappa 5)', 'BIC stop']
    text = RE.p1_text(fro['rep'], plist)
    for label, d, ci, word in a2000:
        assert '* %s: dF = %s, %s' % (label, X.fmt_iv(d, ci), word) in text                                                  # the same numbers and wording as `repeat.md` prints
    _, a40k = E.p1_lines(fro['rep'], plist, B=40000)
    assert [x[1] for x in a40k] == [x[1] for x in a2000] and [tuple(x[2]) for x in a40k] != [tuple(x[2]) for x in a2000]                    # the point estimate does not depend on B, the interval does
    cf = guarded(fro, key='bic')
    sec = E.report(cf, fro).split('### 5.3 ')[1].split('### 5.4 ')[0]
    assert 'identical to the lines of `repeat.md`' in sec and 'NOT identical' not in sec
    monkeypatch.setattr(RE, 'p1_text', lambda raw, plist: 'another text')
    assert 'NOT identical to `repeat.md`: a defect of this report' in E.report(cf, fro).split('### 5.3 ')[1].split('### 5.4 ')[0]
    monkeypatch.undo()
    fix = E.repeat_corrected(fro, cf)
    _, g2000 = E.p1_lines(fix, plist)
    assert g2000[0][1] == a2000[0][1] and g2000[1][1] != a2000[1][1]                                                         # only the default readout's isotropic arm changed
    for (label, d0, ci0, w0), (_, d1, ci1, w1) in zip(a2000, g2000):
        assert label in sec
    row = [l for l in sec.split('\n') if l.startswith('| a priori default (kappa 5) |')][0]
    assert X.fmt_iv(g2000[1][1], g2000[1][2]) in row and X.fmt_iv(a2000[1][1], a2000[1][2]) in row and X.fmt_iv(a40k[1][1], a40k[1][2]) in row
    assert 'Only the isotropic a priori default arm of m2r2 (133, 281) is capped; the tuned line does not change' in sec


def test_the_extras_give_the_residual_scale_of_the_corrected_catalogues():
    fro = world()
    cf = guarded(fro, key='bic')
    sec = E.report(cf, fro).split('### 5.4 ')[1].split('## 6. ')[0]
    assert '| repeat frame m2r2 patch (133, 281), isotropic kernel, a priori default | 1.20 | 1.30 |' in sec
    fro['psf']['P|m92|x1100_y220|default|rsd'] = np.array([2.56, 1.6])
    sec = E.report(cf, fro).split('### 5.4 ')[1].split('## 6. ')[0]
    assert '| M92 test patch (1100, 220), a priori default | 2.56 | 1.30 |' in sec


def test_the_kernel_study_section_prints_the_m92_part_of_psf_md_with_the_stored_and_the_corrected_frozen_rows():
    fro = world()
    cf = guarded(fro, key='bic')
    sec = E.report(cf, fro).split('### 5.5 ')[1].split('## 6. ')[0]
    corr = E.overlay(fro['m2'], cf, 'G')
    for froz, label in ((fro['m2'], '#### stored frozen rows (as run)'), (corr, '#### **%s**' % E.LABEL)):
        md = PS.report(fro['psf'], froz)
        body = md[md.find('\n## M92') + 1:].split('\n## ')[0]
        first_rows = [l for l in body.split('\n') if l.startswith('| ISS as frozen, a priori default (sigma_b 5')]
        assert label in sec and len(first_rows) == 3 and all(r in sec for r in first_rows)                                  # two patches and the pooled table, printed as psf.md prints them
    stored_rows = [l for l in PS.report(fro['psf'], fro['m2'])[PS.report(fro['psf'], fro['m2']).find('\n## M92'):].split('\n') if l.startswith('| ISS as frozen, a priori default (sigma_b 5')]
    corr_rows = [l for l in PS.report(fro['psf'], corr)[PS.report(fro['psf'], corr).find('\n## M92'):].split('\n') if l.startswith('| ISS as frozen, a priori default (sigma_b 5')]
    assert stored_rows != corr_rows
    assert 'No M92 section' not in sec and 'Not all M92' not in sec and 'M53 (NGC 5024)' not in sec and 'The residual scale of the fitted catalogues' not in sec and 'Stars brighter than 19.5 mag' not in sec
    half = guarded(fro, rids=['T2|m92|x1100_y220|default'])
    assert 'Not all M92 default paths are re-run yet.' in E.report(half, fro).split('### 5.5 ')[1].split('## 6. ')[0]
    fro['psf'] = {k: v for k, v in fro['psf'].items() if not k.endswith('|cfg')}
    assert 'No M92 section in the kernel-study raw.' in E.report(cf, fro).split('### 5.5 ')[1].split('## 6. ')[0]


def test_check_mode_reports_identical_different_and_not_run_paths_with_f_before_and_after():
    fro = world()
    rids = {s['rid'] for s in E.split_inventory(fro, CFG)[1] if s['fam'] in ('T', 'D')}
    differ = {'T|x310_y630|ours', 'D|x110_y630|hp3|rex1.2'}
    cf = checked(fro, guarded(fro), rids=rids, differ=differ, act={'T|x310_y630|ours': 1, 'D|x110_y630|hp3|rex1.2': 2})
    sec = E.report(cf, fro).split('## 6. ')[1].split('## 7. ')[0]
    assert 'Summary: 22 paths identical, 2 different, 41 not run (of 65).' in sec and 'T 17 / 1 / 0; D 5 / 1 / 0; P 0 / 0 / 12; R 0 / 0 / 25; V 0 / 0 / 4' in sec
    assert 'Paths with a guard activation: `T\\|x310_y630\\|ours` (1), `D\\|x110_y630\\|hp3\\|rex1.2` (2).' in sec
    row = [l for l in sec.split('\n') if 'T\\|x310_y630\\|ours' in l and l.startswith('| T |')][0]
    assert cells(row)[2:] == ['1300', '1301', '1', 'DIFFERENT']
    ok = [l for l in sec.split('\n') if 'T\\|x360_y1040\\|default' in l][0]
    assert cells(ok)[2:] == ['1400', '1400', '0', 'identical']
    assert [l for l in sec.split('\n') if 'R\\|m2r1\\|x316_y124\\|iso\\|ours' in l and 'not run' in l] and sec.count('| not run |') == 41
    diff = sec.split('### Paths that differ: F before and after')[1]
    good, poor = X.f_measure(X._pack('x', SC['good'])['x|tot'])[0], X.f_measure(X._pack('x', SC['poor'])['x|tot'])[0]
    assert '| M2 test patch (310, 630), tuned configuration | kappa 2 | %.3f' % good in diff and ' | %.3f (' % poor in diff
    assert diff.count('M2 test patch (310, 630), tuned configuration') == 2                                                  # kappa 2 and the BIC block
    assert 'development patch (110, 630), sigma_b 3, r_ex 1.2 | kappa 5 |' in diff
    assert 'T\\|x360_y1040\\|default' not in diff                                                                              # identical paths are not listed there


def test_a_path_that_is_identical_but_has_an_activation_is_not_called_different():
    fro = world()
    rids = {'T|x310_y630|ours'}
    cf = checked(fro, guarded(fro), rids=rids, act={'T|x310_y630|ours': 1})
    sec = E.report(cf, fro).split('## 6. ')[1].split('## 7. ')[0]
    row = [l for l in sec.split('\n') if 'T\\|x310_y630\\|ours' in l and l.startswith('| T |')][0]
    assert row.rstrip().endswith('| identical |') and 'Paths with a guard activation: `T\\|x310_y630\\|ours` (1).' in sec                  # the result is the comparison; the activation is listed separately


def test_the_limitations_name_the_check_coverage_and_the_unrun_paths():
    fro = world()
    cf = checked(fro, guarded(fro), rids={'T|x310_y630|ours'})
    lim = E.report(cf, fro).split('## 7. ')[1]
    assert '**Check mode: 1 of 65 paths were run, 64 were not.**' in lim and 'Not run: `D|x110_y630|hp3|rex0.8`' in lim and 'greedy' in lim and 'post hoc' in lim
    cf = checked(fro, guarded(fro))
    lim = E.report(cf, fro).split('## 7. ')[1]
    assert '**Check mode: 65 of 65 paths were run.** All of them were run.' in lim


def test_the_frozen_state_table_prints_the_recorded_hashes():
    fro = world()
    cf = guarded(fro)
    sec = E.report(cf, fro).split('## 2. ')[1].split('## 3. ')[0]
    assert '| src/core/imaging.py | ' + '0' * 64 + ' |' in sec and 'python 3.13.0, numpy 2.5.4, scipy 1.18.1' in sec and 'did not record' not in sec
    cf.update({'meta|files_late': np.array(['data/astro/checksums_repeat.json']), 'meta|sha256_late': np.array(['f' * 64])})
    sec = E.report(cf, fro).split('## 2. ')[1].split('## 3. ')[0]
    assert 'did not record' in sec and '| data/astro/checksums_repeat.json | ' + 'f' * 64 + ' |' in sec


def test_f_txt_readouts_and_run_info_read_the_stored_arrays():
    z = {}
    _put(z, 'D|a', _run('good', 1300, cap=True))
    _put(z, 'V|b', _run('mid', 1100, bic=False, k0=False))
    f, c, p = X.f_measure(X._pack('x', SC['good'])['x|tot'])
    assert E.f_txt(z, 'D|a|', 2.0) == '%.3f (%.3f, %.3f)' % (f, c, p) and E.f_txt(z, 'D|a|', 7.0) == '-' and E.f_txt(z, 'V|b|', 0.0) == '-'
    assert E.readouts('D|a', z, 'D|a|') == [2.0, 5.0, 0.0] and E.readouts('T2|m92|x1_y2|default', {'p|k5|tot': 1, 'p|k0|tot': 1, 'p|k2|tot': 1}, 'p|') == [5.0, 0.0]
    assert E.readouts('R|m2r1|x1_y2|iso|ours', {'p|k5|tot': 1, 'p|k0|tot': 1, 'p|k2|tot': 1}, 'p|') == [2.0, 0.0] and E.readouts('V|x_y|snr', {'p|k2|tot': 1}, 'p|') == [2.0]
    assert E.run_info(z, 'D|a|') == (X.MAX_ENTRIES, 100.0, 1270, X.MAX_ENTRIES + 1, True)
    assert E.run_info(z, 'V|b|') == (1100, 100.0, None, None, False)
    z['V|b|secs'] = np.array([5.0, 2500.0])
    assert E.run_info(z, 'V|b|')[4] is True


def test_the_configuration_is_read_from_the_first_test_patch_of_the_frozen_raw():
    assert E.cfg_from_raw({'T|x1_y2|cfg': np.array([4.0, 1.5, 3.0])}) == dict(hp=4.0, rex=1.5, kap=3.0, default=tuple(X.DEFAULT))
    assert E.cfg_from_raw({'T2|m92|x1_y2|cfg': np.array([3.0, 1.2, 2.0]), 'x': 1}) == dict(hp=3.0, rex=1.2, kap=2.0, default=tuple(X.DEFAULT))
    with pytest.raises(ValueError):
        E.cfg_from_raw({'D|x|k2|tot': 1})


def test_describe_names_every_family():
    assert E.describe('D|x110_y630|hp8|rex1.2') == 'development patch (110, 630), sigma_b 8, r_ex 1.2'
    assert E.describe('T|x310_y630|ours') == 'M2 test patch (310, 630), tuned configuration'
    assert E.describe('T2|m92|x1100_y220|default') == 'M92 test patch (1100, 220), a priori default'
    assert E.describe('V|x310_y630|norex') == 'design variant "norex", M2 test patch (310, 630), tuned configuration'
    assert E.describe('P|m53|x760_y590|ours') == 'kernel study, M53 patch (760, 590), tuned configuration'
    assert E.describe('R|m2r2|x133_y281|scl|default') == 'repeat frame m2r2 patch (133, 281), chosen kernel, lengths scaled by the seeing, a priori default'
    assert E.describe('R|m2r3|x110_y501|iso|ours').endswith('isotropic kernel, tuned configuration')


def test_status_and_bars_text():
    cap = X.MAX_ENTRIES
    assert E.status_of(np.array([100, 0, 100, 1, 0, 101, 0, 0]), None) == 'ended before the cap'
    assert E.status_of(np.array([cap, 0, 100, 1, 0, 101, 0, 0]), None) == 'CAP REACHED, no activation: censoring of a long path'
    assert E.status_of(np.array([cap, 0, 100, 1, 0, 101, 3, 0]), None) == 'CAP REACHED with 3 activations'
    assert E.status_of(np.array([100, 0, 100, E.CYCLING_FLAG, 0, 101, 0, 0]), None).endswith('possible residual cycling (a point entered %d times)' % E.CYCLING_FLAG)
    assert E.status_of(np.array([100, 0, 100, E.CYCLING_FLAG - 1, 0, 101, 0, 0]), None) == 'ended before the cap'
    assert E.CYCLING_FLAG == 5                                                                                              # the declared value of the plan: a point entered 5 times or more is flagged
    assert E.status_of(np.array([cap - 1, 0, 100, 1, 0, 101, 0, 0]), None) == 'ended before the cap'                         # one entry short of the cap is not the cap
    assert E.bars_txt(np.zeros((0, 3))) == '-'
    assert E.bars_txt(np.array([[2.5, 1, 2], [3.0, 4, 5]])) == '2.5 (1, 2); 3 (4, 5)'
    long = np.array([[float(i), i, i] for i in range(15)])
    assert E.bars_txt(long, limit=3) == '0 (0, 0); 1 (1, 1); 2 (2, 2); ... (15 in all)'


# ---- the command line -------------------------------------------------------------------------------------------------------------------

def test_from_raw_regenerates_the_report_byte_for_byte(tmp_path):
    fro = world()
    cf = checked(fro, guarded(fro, key='bic'), rids={s['rid'] for s in E.split_inventory(fro, CFG)[1] if s['fam'] in ('T', 'D')}, differ={'T|x310_y630|ours'})
    outdir = tmp_path / 'out'
    outdir.mkdir()
    np.savez_compressed(outdir / 'cyclefix_raw.npz', **cf)
    paths = {}
    for k in ('m2', 'psf', 'rep'):
        paths[k] = str(tmp_path / (k + '_raw.npz'))
        np.savez_compressed(paths[k], **fro[k])
    env = dict(os.environ, CYCLEFIX_OUT=str(outdir), CYCLEFIX_M2_RAW=paths['m2'], CYCLEFIX_PSF_RAW=paths['psf'], CYCLEFIX_REPEAT_RAW=paths['rep'], PYTHONDONTWRITEBYTECODE='1')
    cmd = [sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_astro_cyclefix.py'), '--from-raw']
    first = None
    for _ in range(2):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        data = (outdir / 'cyclefix.md').read_bytes()
        assert first is None or data == first
        first = data
        (outdir / 'cyclefix.md').unlink()
    assert first == E.report(E.load_raw(str(outdir / 'cyclefix_raw.npz')), {k: E.load_raw(paths[k]) for k in paths}).encode()
    assert 'wrote' in r.stdout


def test_the_command_line_runs_the_stage_asked_for(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(E, 'stage1', lambda workers: calls.append(('run', workers)))
    monkeypatch.setattr(E, 'check_mode', lambda workers, only=None, reverse=False: calls.append(('check', workers, only, reverse)))
    monkeypatch.setattr(E, 'list_inventory', lambda: calls.append(('list',)))
    for argv in (['--run', '--workers', '3'], ['--check', '--workers', '1'], ['--check', '--only', 'V,T'], ['--check', '--reverse'], ['--list']):
        monkeypatch.setattr(sys, 'argv', ['exp_astro_cyclefix.py'] + argv)
        E.main()
    assert calls == [('run', 3), ('check', 1, None, False), ('check', 2, ['V', 'T'], False), ('check', 2, None, True), ('list',)]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_cyclefix.py'])
    with pytest.raises(SystemExit):
        E.main()
    out = tmp_path / 'r.md'
    monkeypatch.setattr(E, 'MD', str(out))
    monkeypatch.setattr(E, 'RAW', str(tmp_path / 'raw.npz'))
    np.savez_compressed(tmp_path / 'raw.npz', a=np.zeros(1))
    monkeypatch.setattr(E, 'frozen_raws', lambda: {'x': 1})
    monkeypatch.setattr(E, 'report', lambda cf, fro: 'the report of %s and %s\n' % (sorted(cf), sorted(fro)))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_cyclefix.py', '--from-raw'])
    E.main()
    assert out.read_text() == "the report of ['a'] and ['x']\n"


def test_list_inventory_prints_the_affected_paths_and_the_check_families(monkeypatch, capsys):
    fro = world()
    monkeypatch.setattr(E, 'frozen_raws', lambda: fro)
    monkeypatch.setattr(E, 'frozen_config', lambda: CFG)
    E.list_inventory()
    out = capsys.readouterr().out
    assert 'stage 1 (12 paths' in out and 'D|x110_y630|hp5|rex0.8' in out and 'flagged by entries' in out and 'check mode (65 paths)' in out and 'aliases of an isotropic run (not run): R|m2r3|x110_y501|adq|default' in out
    assert 'T: 18 paths' in out and 'R: 25 paths' in out and 'V: 4 paths' in out
