"""The held-out repeat-frame test (experiments/exp_astro_repeat.py): the freeze record and its check, the per-patch bookkeeping (the isotropic run stored under both kernel names where the rule keeps
the isotropic kernel, the lengths scaled by the seeing, the baselines' names), write-once running, and the report read off a planted raw."""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_common as A                # noqa: E402
import astro_repeat as R                # noqa: E402
import astro_second as S2               # noqa: E402
import exp_astro_m2 as X                # noqa: E402
import exp_astro_psf as PS              # noqa: E402
import exp_astro_repeat as E            # noqa: E402
from src.core import imaging as IM      # noqa: E402
from test_astro_report import _synthetic_catalogs   # noqa: E402

ISO = ((0.53, 1.0), (0.35, 1.9), (0.12, 4.7))
NAMES = {'photutils': 'photutils', 'photutils_tuned': 'photutils_ext_all_local_i6', 'cs_default': X.cs_name((0.0, 0.0, 0.0, np.inf)), 'cs_tuned': 'crowdsource|snr5|ff0|qf0|spinf'}


def _score(cats, truth, key):
    return A.score(np.c_[cats[key][0], cats[key][1]], cats[key][2], truth)


# ---- the frozen configuration and the freeze record -------------------------------------------------------------------------------------

def test_load_frozen_reads_the_m2_freeze_and_refuses_another_default(tmp_path, monkeypatch):
    rec = dict(config=dict(sigma_b=4.0, r_ex=1.5, kappa=3.0), default=list(X.DEFAULT), tuned_baselines=dict(photutils='p', crowdsource='c'))
    path = tmp_path / 'freeze.json'
    path.write_text(json.dumps(rec))
    fr = E.load_frozen(str(path))
    assert (fr['hp'], fr['rex'], fr['kap']) == (4.0, 1.5, 3.0) and fr['default'] == tuple(X.DEFAULT) and fr['tuned'] == dict(photutils='p', crowdsource='c')
    monkeypatch.setattr(X, 'DEFAULT', (5.0, 1.2, 6.0))
    with pytest.raises(RuntimeError):
        E.load_frozen(str(path))


def test_the_real_frozen_configuration_is_the_m2_one():
    fr = E.load_frozen()
    assert (fr['hp'], fr['rex'], fr['kap']) == (3.0, 1.2, 2.0) and fr['default'] == (5.0, 1.2, 5.0)
    assert fr['tuned'] == dict(photutils='photutils_ext_all_local_i6', crowdsource='crowdsource|snr5|ff0|qf0|spinf')


def test_the_freeze_is_written_checked_and_refuses_any_difference(tmp_path, monkeypatch):
    rec = dict(config=dict(sigma_b=3.0), files={'a': 'x'}, frozen_code='abc', data_files={'d': 'y'}, fwhm_ref=2.15)
    monkeypatch.setattr(E, 'freeze_record', lambda: json.loads(json.dumps(rec)))
    path = str(tmp_path / 'freeze_repeat.json')
    with pytest.raises(RuntimeError, match='no freeze record'):
        E.check_freeze(path)
    E.write_freeze(path)
    assert E.check_freeze(path) == rec
    for key, val in (('frozen_code', 'abd'), ('files', {'a': 'z'}), ('data_files', {'d': 'w'}), ('fwhm_ref', 2.16)):
        changed = dict(rec)
        changed[key] = val
        monkeypatch.setattr(E, 'freeze_record', lambda changed=changed: json.loads(json.dumps(changed)))
        with pytest.raises(RuntimeError, match=key):
            E.check_freeze(path)


def test_the_code_hash_covers_the_functions_and_constants_that_decide_a_number(monkeypatch):
    h0 = E.frozen_code_hash()
    assert h0 == E.frozen_code_hash() and len(h0) == 64
    monkeypatch.setattr(E, 'PATCHES', E.PATCHES[:-1])
    assert E.frozen_code_hash() != h0
    monkeypatch.undo()

    def other_cs_name(f):
        return 'x'
    monkeypatch.setattr(X, 'cs_name', other_cs_name)
    assert E.frozen_code_hash() != h0
    monkeypatch.undo()
    monkeypatch.setattr(X, 'KAPPAS', X.KAPPAS[1:])
    assert E.frozen_code_hash() != h0


needs_data = pytest.mark.skipif(not all(os.path.exists(R.patch_path(p[0], p[1], n)) for n, p in E.PATCHES) or not os.path.exists(X.patch_path(110, 630)),
                                reason='the repeat and M2 patches are not on disk (data/astro/fetch_repeat.py, data/astro/fetch.py)')


@needs_data
def test_the_real_freeze_record_names_every_input():
    rec = E.freeze_record()
    assert set(rec['files']) == set(E.FILES) and len(rec['data_files']) == len(E.PATCHES) + len(R.REPEAT) and rec['frozen_code'] == E.frozen_code_hash()
    assert rec['config'] == dict(sigma_b=3.0, r_ex=1.2, kappa=2.0) and rec['default'] == [5.0, 1.2, 5.0] and 2.0 < rec['fwhm_ref'] < 2.3
    assert rec['patches'] == [[n, list(p)] for n, p in E.PATCHES]


# ---- the bookkeeping of one patch -------------------------------------------------------------------------------------------------------

def _stub_world(monkeypatch, kind, rng_seed=7):
    """Replaces the data, the runs and the scorers of `patch_arrays` by planted ones; returns the list in which every call of a run is recorded."""
    rng = np.random.default_rng(rng_seed)
    truth, cats = _synthetic_catalogs(rng)
    px, py = truth['x'] + 316.0, truth['y'] + 124.0                                  # the truth of the patch with origin (316, 124)
    calls = []
    monkeypatch.setattr(R, 'load_patch', lambda x0, y0, name: dict(dn=np.full((100, 100), 150.0), psf=np.full((51, 51), 1.0 / 2601), gain=4.6, dark=1.0))
    monkeypatch.setattr(R, 'truth_arrays', lambda name: dict(px=px, py=py, m606=truth['mag']))
    fit = dict(kind=kind, sigma=0.9, parts=ISO if kind == 'isotropic' else None, parts5=None if kind == 'isotropic' else ((1.0, 1.0, 1.0, 0.0, 0.0),), rms_isotropic=0.001, rms=0.001)
    monkeypatch.setattr(PS, 'stamp_arrays', lambda patch: (dict(kind=np.array([kind]), rms=np.array([0.001, 0.001]), sigma=np.array([0.9]), parts5=np.zeros((0, 5)), shape=np.zeros(4)), fit))
    monkeypatch.setattr(S2, 'psf_fwhm', lambda psf: 3.0)
    monkeypatch.setattr(IM, 'fit_psf_parts', lambda stamp, n_parts=3, sigma0=1.0: (0.9, ISO, 0.001))
    monkeypatch.setattr(PS, 'residual_sd', lambda patch, P, det: 1.05)

    def dets(key):
        return {k: (cats[key][0], cats[key][1], cats[key][2]) for k in (2.0, 5.0, 0.0)}

    def run_path(patch, h, r, **kw):
        calls.append(('iso', h, r))
        return dets('mid'), 1.5, 40, dict(bic_index=3, n_blocks=9)

    def run_path_psf(patch, h, r, fit_, **kw):
        calls.append(('psf', h, r))
        return dets('good'), 2.5, 41, dict(bic_index=4, n_blocks=9, rsd={2.0: 1.1, 5.0: 1.2, 0.0: 1.3})
    monkeypatch.setattr(X, 'run_path', run_path)
    monkeypatch.setattr(PS, 'run_path_psf', run_path_psf)

    def baseline_scores(pxy, t, base=None):
        return {NAMES['photutils']: _score(cats, truth, 'poor'), NAMES['photutils_tuned']: _score(cats, truth, 'mid'), NAMES['cs_default']: _score(cats, truth, 'poor'), NAMES['cs_tuned']: _score(cats, truth, 'bic')}
    monkeypatch.setattr(X, 'baseline_scores', baseline_scores)
    return calls, truth, cats


def _rounded(calls):
    return [(c[0], round(c[1], 9), round(c[2], 9)) for c in calls]


ARGS = ('m2r1', (316, 124), 3.0, 1.2, 2.0, dict(photutils=NAMES['photutils_tuned'], crowdsource=NAMES['cs_tuned']), 2.0, '/nonexistent')


def test_patch_arrays_stores_one_run_under_both_kernel_names_where_the_rule_keeps_the_isotropic_kernel(monkeypatch):
    calls, truth, cats = _stub_world(monkeypatch, 'isotropic')
    name, pxy, arr = E.patch_arrays(ARGS)
    assert (name, pxy) == ('m2r1', (316, 124)) and arr['kernel|same'][0] == 1
    assert calls[:2] == [('iso', 3.0, 1.2), ('iso', 5.0, 1.2)]                                   # the frozen runs: the tuned configuration and the a priori default
    assert sorted(_rounded(calls[2:])) == [('psf', 4.5, 1.8), ('psf', 7.5, 1.8)] and len(calls) == 4      # no adequate-kernel run; the scaled runs at s = 3.0 / 2.0
    for lab in ('ours', 'default'):
        for k in ('k2', 'k5', 'k0'):
            assert np.array_equal(arr['iso|%s|%s|tot' % (lab, k)], arr['adq|%s|%s|tot' % (lab, k)]) and np.array_equal(arr['iso|%s|%s|blk22' % (lab, k)], arr['adq|%s|%s|blk22' % (lab, k)])
        assert np.array_equal(arr['iso|%s|det' % lab], arr['adq|%s|det' % lab]) and arr['iso|%s|rsd' % lab][0] == arr['adq|%s|rsd' % lab][0] == 1.05
    assert arr['scale'][0] == 1.5 and np.array_equal(arr['fwhm'], [3.0, 2.0]) and np.array_equal(arr['cfg'], [3.0, 1.2, 2.0])
    assert np.array_equal(arr['scl|ours|rsd'], [1.1]) and np.array_equal(arr['scl|default|rsd'], [1.2])      # the residual scale at the reported significance of each configuration


def test_patch_arrays_runs_the_chosen_kernel_where_the_isotropic_one_is_inadequate(monkeypatch):
    calls, truth, cats = _stub_world(monkeypatch, 'general')
    name, pxy, arr = E.patch_arrays(ARGS)
    assert arr['kernel|same'][0] == 0
    assert _rounded(calls) == [('iso', 3.0, 1.2), ('psf', 3.0, 1.2), ('iso', 5.0, 1.2), ('psf', 5.0, 1.2), ('psf', 4.5, 1.8), ('psf', 7.5, 1.8)]
    mid, good = _score(cats, truth, 'mid'), _score(cats, truth, 'good')
    assert np.array_equal(arr['iso|ours|k2|tot'], X._pack('x', mid)['x|tot']) and np.array_equal(arr['adq|ours|k2|tot'], X._pack('x', good)['x|tot'])
    assert np.array_equal(arr['adq|ours|rsd'], [1.1]) and np.array_equal(arr['adq|default|rsd'], [1.2])          # info['rsd'] at kappa 2 and at the default's kappa 5
    assert np.array_equal(arr['adq|ours|bic'], [4, 9, 0]) and np.array_equal(arr['iso|ours|bic'], [3, 9, 0]) and np.array_equal(arr['adq|ours|secs'], [2.5, 41])


def test_patch_arrays_marks_a_capped_run_and_stores_the_baselines_by_their_frozen_names(monkeypatch):
    calls, truth, cats = _stub_world(monkeypatch, 'isotropic')
    monkeypatch.setattr(X, 'run_path', lambda patch, h, r, **kw: ({k: (cats['mid'][0], cats['mid'][1], cats['mid'][2]) for k in (2.0, 5.0, 0.0)}, 1.0, X.MAX_ENTRIES, dict(bic_index=3, n_blocks=9)))
    _, _, arr = E.patch_arrays(ARGS)
    assert np.array_equal(arr['iso|ours|bic'], [3, 9, 1])
    poor, mid, bic = (_score(cats, truth, k) for k in ('poor', 'mid', 'bic'))
    for label, sc in (('photutils (default)', poor), ('photutils (tuned)', mid), ('crowdsource (default)', poor), ('crowdsource (tuned)', bic)):
        assert np.array_equal(arr[label + '|tot'], X._pack('x', sc)['x|tot']), label
    assert str(arr['photutils (tuned)|name'][0]) == NAMES['photutils_tuned'] and str(arr['crowdsource (tuned)|name'][0]) == NAMES['cs_tuned']


def test_the_scaled_run_uses_the_seeing_of_the_patch_over_the_reference(monkeypatch):
    calls, _, _ = _stub_world(monkeypatch, 'isotropic')
    monkeypatch.setattr(S2, 'psf_fwhm', lambda psf: 4.2)
    E.patch_arrays(ARGS[:6] + (3.5, '/x'))
    assert sorted(_rounded([c for c in calls if c[0] == 'psf'])) == [('psf', 3.6, 1.44), ('psf', 6.0, 1.44)]


# ---- write-once running -----------------------------------------------------------------------------------------------------------------

_CALLS = []


def _stub_pa(args):
    """A module-level stand-in for `patch_arrays` (a pool pickles its function by name); it records the call in this process."""
    name, pxy = args[0], args[1]
    _CALLS.append((name, pxy))
    return name, pxy, {'cfg': np.array([1.0, 2.0, 3.0]), 'x': np.arange(3.0) + pxy[0]}


@pytest.mark.parametrize('workers', [1, 2])
def test_part_run_runs_each_missing_patch_once_and_never_overwrites(tmp_path, monkeypatch, workers):
    patches = [('m2r1', (316, 124)), ('m2r2', (133, 281)), ('m2r3', (370, 441))]
    monkeypatch.setattr(E, 'PATCHES', patches)
    monkeypatch.setattr(E, 'check_freeze', lambda path=None: dict(fwhm_ref=2.15))
    monkeypatch.setattr(E, 'load_frozen', lambda path=None: dict(hp=3.0, rex=1.2, kap=2.0, default=(5.0, 1.2, 5.0), tuned={}))
    monkeypatch.setattr(E, 'patch_arrays', _stub_pa)
    _CALLS.clear()
    out = str(tmp_path / 'repeat_raw.npz')
    E.part_run(out_path=out, workers=workers)
    raw = X.load_raw(out)
    assert sorted(k for k in raw if k.endswith('|cfg')) == ['R|m2r1|x316_y124|cfg', 'R|m2r2|x133_y281|cfg', 'R|m2r3|x370_y441|cfg']
    assert np.array_equal(raw['R|m2r2|x133_y281|x'], np.arange(3.0) + 133)
    if workers == 1:
        assert sorted(_CALLS) == sorted(patches)
    # a second run does nothing; with one result removed only that patch is run, and the others are untouched
    _CALLS.clear()
    E.part_run(out_path=out, workers=1)
    assert _CALLS == []
    keep = {k: v for k, v in raw.items() if not k.startswith('R|m2r2|')}
    keep['R|m2r1|x316_y124|x'] = np.full(3, -1.0)                                           # a stored value is never recomputed
    np.savez_compressed(out, **keep)
    E.part_run(out_path=out, workers=1)
    assert _CALLS == [('m2r2', (133, 281))]
    raw = X.load_raw(out)
    assert np.array_equal(raw['R|m2r1|x316_y124|x'], np.full(3, -1.0)) and 'R|m2r2|x133_y281|cfg' in raw


def test_part_run_refuses_to_start_when_the_freeze_does_not_match(tmp_path, monkeypatch):
    def refuse(path=None):
        raise RuntimeError('the freeze record differs')
    monkeypatch.setattr(E, 'check_freeze', refuse)
    _CALLS.clear()
    monkeypatch.setattr(E, 'patch_arrays', _stub_pa)
    with pytest.raises(RuntimeError, match='freeze'):
        E.part_run(out_path=str(tmp_path / 'raw.npz'), workers=1)
    assert _CALLS == [] and not (tmp_path / 'raw.npz').exists()


# ---- the report -------------------------------------------------------------------------------------------------------------------------

def test_purity_over_all_magnitudes_and_the_unmatched_detections():
    tot = np.array([100, 60, 80, 55, 300, 140, 400.0])
    p, un = E.purity_all(tot)
    assert p == 0.35 and un == 260 and isinstance(un, int)


def _planted(monkeypatch, same=(False, False), adq_key='good', iso_key='mid', seed=11):
    """A raw with two patches (frames m2r1 and m2r2) in the layout `patch_arrays` stores, and a frozen raw with the four M2 windows."""
    rng = np.random.default_rng(seed)
    truth, cats = _synthetic_catalogs(rng)
    plist = [('m2r1', (316, 124)), ('m2r2', (133, 281))]
    monkeypatch.setattr(E, 'PATCHES', plist)
    raw = {}
    for (name, pxy), is_same in zip(plist, same):
        pre = E.pre_of(name, pxy)
        raw[pre + 'cfg'] = np.array([3.0, 1.2, 2.0])
        raw[pre + 'kind'] = np.array(['isotropic' if is_same else 'general'])
        raw[pre + 'rms'] = np.array([0.0012, 0.0012]) if is_same else np.array([0.0072, 0.0011])
        raw[pre + 'parts5'] = np.zeros((0, 5)) if is_same else np.zeros((6, 5))
        raw[pre + 'shape'] = np.array([0.11, -60.0, 1.4, 1.5])
        raw[pre + 'fwhm'] = np.array([2.28 if name == 'm2r1' else 3.05, 2.15])
        raw[pre + 'scale'] = np.array([1.06 if name == 'm2r1' else 1.42])
        raw[pre + 'kernel|same'] = np.array([int(is_same)])
        for kern, key in (('iso', iso_key), ('adq', iso_key if is_same else adq_key), ('scl', adq_key)):
            for lab in E.LABELS:
                for k in ('k2', 'k5', 'k0'):
                    for kk, v in X._pack('%s|%s|%s' % (kern, lab, k), _score(cats, truth, key)).items():
                        raw[pre + kk] = v
                raw[pre + '%s|%s|secs' % (kern, lab)] = np.array([120.0, 900.0])
                raw[pre + '%s|%s|bic' % (kern, lab)] = np.array([300, 901, 0])
                raw[pre + '%s|%s|rsd' % (kern, lab)] = np.array([1.0 + 0.3 * (kern == 'iso')])
        raw[pre + 'iso|default|bic'] = np.array([300, 901, 1])
        for label, key in (('crowdsource (default)', 'poor'), ('crowdsource (tuned)', 'mid'), ('photutils (default)', 'poor'), ('photutils (tuned)', 'poor')):
            for kk, v in X._pack(label, _score(cats, truth, key)).items():
                raw[pre + kk] = v
    froz = {}
    for w in X.TESTS:
        for label, key in (('ours|k2', 'good'), ('crowdsource (tuned)', 'mid'), ('photutils (tuned)', 'poor')):
            for kk, v in X._pack(label, _score(cats, truth, key)).items():
                froz['T|x%d_y%d|%s' % (w[0], w[1], kk)] = v
    return raw, froz, truth, cats


def test_the_report_is_deterministic_and_reads_the_raw(monkeypatch):
    raw, froz, truth, cats = _planted(monkeypatch)
    md = E.report(raw, froz)
    assert md == E.report(raw, froz)
    totals = sum(raw[E.pre_of(*np_) + 'adq|ours|k2|tot'] for np_ in E.PATCHES)
    row = [l for l in md.split('\n') if l.startswith('| ISS, kernel chosen from the stamp, tuned configuration (sigma_b 3, r_ex 1.2, kappa 2) | %d / %d |' % (totals[1], totals[0]))]
    assert len(row) >= 1                                                                     # the pooled row is read from the chosen-kernel arrays, not the isotropic ones
    pa, un = E.purity_all(totals)
    assert '%.3f (%d)' % (pa, un) in row[0]
    for needle in ('## The patches and their stamp fits (label-free)', '| m2r1 | (316, 124) | (310, 630) | 2.28 |', '| m2r2 | (133, 281) | (310, 630) | 3.05 |', '## All 2 patches together', '### Primary outcome P1',
                   '## Frame m2r1', '## Frame m2r2', '## Record by patch', '## P3: the residual scale', '| m2r1 | (316, 124) | 1.30 | 1.00 | 1.30 | 1.00 |', '## Run records',
                   '## Seeing ladder on the M2 windows', '## Exploratory: lengths scaled by the seeing', 'm2r1 1.06', 'm2r2 1.42'):
        assert needle in md, needle
    assert 'nan' not in md.lower().replace('nan,', '')


def test_the_report_flags_capped_runs_and_skips_the_second_kernel_where_the_rule_keeps_the_isotropic_one(monkeypatch):
    raw, froz, *_ = _planted(monkeypatch, same=(True, False))
    md = E.report(raw, froz)
    cap_rows = [l for l in md.split('\n') if l.startswith('| m2r') and l.rstrip().endswith('| capped |')]
    assert len(cap_rows) == 2 and all('isotropic kernel (as frozen)' in l and ' default ' in l for l in cap_rows)      # only the isotropic default runs were capped
    run_rows = [l for l in md.split('\n') if l.startswith('| m2r1 | (316, 124) |') and 'BIC block' not in l and '| 900 |' in l]
    assert len(run_rows) == 2                                                                  # m2r1 keeps the isotropic kernel: its two configurations are listed once, not under the chosen kernel
    assert [l for l in md.split('\n') if l.startswith('| m2r1 | (316, 124) | 1.30 | same | 1.30 | same |')]


def test_p1_confirms_a_better_chosen_kernel_and_has_nothing_to_test_when_the_rule_keeps_the_isotropic_one(monkeypatch):
    raw, froz, *_ = _planted(monkeypatch, adq_key='good', iso_key='poor')
    txt = E.p1_text(raw, E.done_patches(raw))
    assert 'the 2 patches where the rule picks a general kernel' in txt and txt.count('confirmed on held-out frames') == 3 and 'not confirmed' not in txt
    raw, *_ = _planted(monkeypatch, adq_key='poor', iso_key='good')
    txt = E.p1_text(raw, E.done_patches(raw))
    assert txt.count('contradicted') == 3
    raw, *_ = _planted(monkeypatch, same=(True, True))
    assert 'dF = 0 by construction' in E.p1_text(raw, E.done_patches(raw))
    raw, *_ = _planted(monkeypatch, same=(True, False), adq_key='good', iso_key='mid')
    assert 'the 1 patches where the rule picks a general kernel (m2r2 (133, 281))' in E.p1_text(raw, E.done_patches(raw))


def test_paired_differences_are_the_first_catalogue_minus_the_second(monkeypatch):
    raw, froz, *_ = _planted(monkeypatch, adq_key='good', iso_key='poor')
    md = E.report(raw, froz)
    import re
    for name in ('adequate kernel vs isotropic, tuned configuration', 'adequate kernel (tuned configuration) vs crowdsource (tuned)'):
        assert re.search(r'\| %s \| [^\n]* \| win \|' % re.escape(name), md), name
    assert re.search(r'\| isotropic kernel \(tuned configuration\) vs crowdsource \(tuned\) \| [^\n]* \| loss \|', md)


def test_the_ladder_lists_the_frames_that_contain_each_window(monkeypatch):
    raw, froz, *_ = _planted(monkeypatch)
    md = E.report(raw, froz)
    lines = [l for l in md.split('\n') if l.startswith('| (310, 630) |')]
    assert len(lines) == 1 and lines[0].count('ISS') == 3 and 'M2: ISS' in lines[0] and 'm2r1: ISS' in lines[0] and 'm2r2: ISS' in lines[0]       # (310, 630) is in M2, m2r1 and m2r2
    lines = [l for l in md.split('\n') if l.startswith('| (360, 1040) |')]
    assert len(lines) == 1 and lines[0].count('ISS') == 1                                      # only M2 among the two planted frames contains (360, 1040)
