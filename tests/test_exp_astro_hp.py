"""The baselines given the ISS's front end (experiments/exp_astro_hp.py and the registries of experiments/exp_astro_baselines_hp.py): the selection rule on planted scores, the freeze of the selection, the
scoring of the development and test patches from planted catalogues, the frozen rows and the paired differences, the report read off a planted raw, and the byte-identity of `--from-raw`."""
import hashlib
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
import exp_astro_baselines_hp as G      # noqa: E402
import exp_astro_hp as E                # noqa: E402
import exp_astro_m2 as X                # noqa: E402
from test_astro_report import _synthetic_catalogs   # noqa: E402

PHU, CS = G.FAMILIES['photutils'], G.FAMILIES['crowdsource']


def tot(true22, matched22, det22, det22_matched, true_all=0, matched_all=0, det_all=0):
    return np.array([true22, matched22, det22, det22_matched, true_all, matched_all, det_all], float)


def _dev_raw(best_phu=None, best_cs=None, base=(100, 20, 100, 20)):
    """Planted development totals: every variant of both families has the pooled F of `base` (0.2) except the given ones, which have a better F on both development patches."""
    raw = {}
    for tags in (PHU, CS):
        for t in tags:
            for p in X.DEV:
                raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*base)
    for t, vals in ((best_phu, (100, 60, 100, 70)), (best_cs, (100, 50, 100, 55))):
        if t:
            for p in X.DEV:
                raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*vals)
    return raw


# ---- the selection rule -----------------------------------------------------------------------------------------------------------------------

def test_choose_hp_takes_the_largest_pooled_f_of_each_family():
    for tp, tc in ((PHU[3], CS[5]), (PHU[0], CS[0]), (PHU[7], CS[7])):
        assert E.choose_hp(_dev_raw(tp, tc)) == {'photutils': tp, 'crowdsource': tc}
    # the F is the F-measure of completeness and purity from the summed counts: (100, 60, 100, 70) on two patches
    raw = _dev_raw(PHU[2], CS[2])
    assert abs(E.dev_F(raw, PHU[2]) - X.f_measure(tot(200, 120, 200, 140))[0]) < 1e-12 and abs(E.dev_F(raw, PHU[0]) - 0.2) < 1e-12


def test_choose_hp_chooses_inside_each_family_only():
    raw = _dev_raw(PHU[1], CS[1])
    for p in X.DEV:                                       # a crowdsource variant much better than every photutils variant does not become the photutils choice
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], CS[4])] = tot(100, 95, 100, 95)
    assert E.choose_hp(raw) == {'photutils': PHU[1], 'crowdsource': CS[4]}


def test_choose_hp_goes_to_the_larger_tag_on_a_tie_after_rounding_to_six_decimals():
    raw = _dev_raw()
    a, b = sorted([PHU[1], PHU[6]])
    assert a < b
    for p in X.DEV:
        for t in (a, b):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(100, 60, 100, 70)
    assert E.choose_hp(raw)['photutils'] == b                                 # identical totals: the larger tag
    # the smaller tag is better by 1e-9 in F: still a tie at six decimals, so the larger tag
    for p in X.DEV:
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], a)] = tot(10 ** 9, 600000000, 10 ** 9, 700000001)
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], b)] = tot(10 ** 9, 600000000, 10 ** 9, 700000000)
    assert E.dev_F(raw, a) > E.dev_F(raw, b) and round(E.dev_F(raw, a), 6) == round(E.dev_F(raw, b), 6) and E.choose_hp(raw)['photutils'] == b
    # better at the sixth decimal: the smaller tag wins
    for p in X.DEV:
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], a)] = tot(10 ** 9, 600000000, 10 ** 9, 700010000)
    assert E.choose_hp(raw)['photutils'] == a
    # the worst variant is never chosen
    for p in X.DEV:
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], PHU[7])] = tot(100, 1, 100, 1)
    assert E.choose_hp(raw)['photutils'] == a


def test_choose_hp_pools_the_counts_of_both_development_patches_and_needs_both():
    A_, B_ = PHU[2], PHU[5]
    raw = _dev_raw()
    for t, vals in ((A_, [(100, 60, 100, 60), (100, 10, 100, 10)]), (B_, [(100, 40, 100, 40), (100, 40, 100, 40)])):
        for p, v in zip(X.DEV, vals):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*v)
    assert E.choose_hp(raw)['photutils'] == B_                                  # pooled 70 / 200 against 80 / 200
    for t_, vals in ((A_, [(100, 10, 100, 10), (100, 60, 100, 60)]), (B_, [(100, 40, 100, 40), (100, 40, 100, 40)])):         # the same with the patches exchanged: the first patch alone prefers B, the second A
        for p, v in zip(X.DEV, vals):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t_)] = tot(*v)
    assert E.choose_hp(raw)['photutils'] == B_
    # a variant scored on one development patch only is not eligible, however good
    raw['DEV|x%d_y%d|%s|tot' % (X.DEV[0][0], X.DEV[0][1], PHU[0])] = tot(100, 100, 100, 100)
    raw.pop('DEV|x%d_y%d|%s|tot' % (X.DEV[1][0], X.DEV[1][1], PHU[0]))
    assert E.choose_hp(raw)['photutils'] == B_
    # a family without an eligible variant is absent from the answer
    assert 'crowdsource' not in E.choose_hp({k: v for k, v in raw.items() if not any(('|%s|' % t) in k for t in CS)})


# ---- the freeze of the selection ---------------------------------------------------------------------------------------------------------

def _frozen_raw(tmp_path, sel=None):
    plan = tmp_path / 'plan.md'
    plan.write_text('the plan\n')
    raw = _dev_raw(PHU[3], CS[5])
    raw['SEL|photutils'] = np.array([PHU[3]])
    raw['SEL|crowdsource'] = np.array([CS[5]])
    raw['SEL|plan_sha256'] = np.array([hashlib.sha256(b'the plan\n').hexdigest()])
    raw['SEL|dev_hash'] = np.array([E.dev_hash(raw)])
    return raw, str(plan)


def test_verify_selection_returns_the_stored_selection_and_refuses_every_change(tmp_path):
    raw, plan = _frozen_raw(tmp_path)
    assert E.verify_selection(raw, plan) == {'photutils': PHU[3], 'crowdsource': CS[5]}
    with pytest.raises(RuntimeError, match='no selection'):
        E.verify_selection({k: v for k, v in raw.items() if not k.startswith('SEL|')}, plan)
    for key in ('SEL|photutils', 'SEL|crowdsource'):                              # the selection of each family is needed
        with pytest.raises(RuntimeError, match='no selection'):
            E.verify_selection({k: v for k, v in raw.items() if k != key}, plan)
    wrong = dict(raw)
    wrong['SEL|photutils'] = np.array([PHU[4]])                                  # not what the rule gives
    with pytest.raises(RuntimeError, match='not what the rule gives'):
        E.verify_selection(wrong, plan)
    moved = dict(raw)                                                            # a development score changed after the selection (a different variant now wins)
    for p in X.DEV:
        moved['DEV|x%d_y%d|%s|tot' % (p[0], p[1], PHU[0])] = tot(100, 90, 100, 90)
    with pytest.raises(RuntimeError):
        E.verify_selection(moved, plan)
    quiet = dict(raw)                                                            # a change that leaves the selection as it is: caught by the hash of the development scores
    quiet['DEV|x%d_y%d|%s|tot' % (X.DEV[0][0], X.DEV[0][1], PHU[1])] = tot(100, 21, 100, 20)
    assert E.choose_hp(quiet) == E.choose_hp(raw)
    with pytest.raises(RuntimeError, match='development scores changed'):
        E.verify_selection(quiet, plan)
    open(plan, 'w').write('the plan, edited\n')
    with pytest.raises(RuntimeError, match='plan file changed'):
        E.verify_selection(raw, plan)


def test_dev_hash_depends_on_every_development_key_and_nothing_else():
    raw = _dev_raw(PHU[0], CS[0])
    h = E.dev_hash(raw)
    assert len(h) == 64 and h == E.dev_hash(dict(raw, other=np.zeros(2), **{'SEL|photutils': np.array(['x'])}))
    k = 'DEV|x%d_y%d|%s|tot' % (X.DEV[1][0], X.DEV[1][1], CS[7])
    assert E.dev_hash(dict(raw, **{k: raw[k] + 1})) != h
    assert E.dev_hash({a: b for a, b in raw.items() if a != k}) != h
    shuffled = dict(reversed(list(raw.items())))
    assert list(shuffled) != list(raw) and E.dev_hash(shuffled) == h                                 # the keys are hashed in sorted order
    renamed = {(a + 'x' if a == k else a): b for a, b in raw.items()}                                 # one key renamed to its neighbour in the sorted order: the values and their order are the same
    assert sorted(renamed)[sorted(raw).index(k)] == k + 'x' and E.dev_hash(renamed) != h                 # the names are hashed with the values


def test_test_tags_are_the_selected_variants_and_the_references():
    t = E.test_tags({'photutils': PHU[6], 'crowdsource': CS[2]})
    assert t == {'phu_sel': PHU[6], 'phu_ref': 'phu_hp5_new_i3', 'cs_sel': CS[2], 'cs_hp': 'cs_hp5_t5', 'cs_default': 'cs_default'} and tuple(t) == E.NEW_ROWS


# ---- scoring the planted catalogues ------------------------------------------------------------------------------------------------------

QUAL = ['good', 'mid', 'poor', 'bic']
TRUTH_CALLS = []                                  # which truth list was asked for, by whom (X: exp_astro_m2, R: astro_repeat)
TFP_CALLS = []                                    # the patch origins for which the truth of a patch was cut out


def _planted_world(monkeypatch, tmp_path, seed=3, only_tags=None):
    """The truth of every patch is one planted truth list (patch coordinates), the catalogue file holds a planted catalogue of one of four qualities for every (patch, tag) of the sets used."""
    rng = np.random.default_rng(seed)
    truth, cats = _synthetic_catalogs(rng)
    monkeypatch.setattr(A, 'truth_for_patch', lambda tr, x0, y0, n=100, mag_max=30.0: TFP_CALLS.append((x0, y0)) or truth)
    monkeypatch.setattr(X, 'truth_arrays', lambda ds='m2': TRUTH_CALLS.append(('X', ds)) or dict(px=np.zeros(1), py=np.zeros(1), m606=np.zeros(1)))
    monkeypatch.setattr(R, 'truth_arrays', lambda name: TRUTH_CALLS.append(('R', name)) or dict(px=np.zeros(1), py=np.zeros(1), m606=np.zeros(1)))
    TRUTH_CALLS.clear()
    TFP_CALLS.clear()
    out = {}
    plist = [p for s in ('dev', 'm2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3') for p in G.SETS[s][1]]
    tags = PHU + CS
    for pi, p in enumerate(plist):
        for ti, tag in enumerate(tags):
            if only_tags and tag not in only_tags:
                continue
            q = cats[QUAL[(pi + 3 * ti) % 4]]
            k = 'x%d_y%d|%s|' % (p[0], p[1], tag)
            out[k + 'x'], out[k + 'y'], out[k + 'flux'], out[k + 'secs'] = q[0], q[1], q[2], np.array([1.5 + ti])
    path = tmp_path / 'cats.npz'
    np.savez_compressed(path, **out)
    return truth, cats, str(path)


def _score_planted(cats, truth, qual):
    return A.score(np.c_[cats[qual][0], cats[qual][1]], cats[qual][2], truth)


def test_part_dev_scores_every_variant_selects_by_the_rule_and_refuses_a_second_selection(monkeypatch, tmp_path):
    truth, cats, cpath = _planted_world(monkeypatch, tmp_path)
    plan = tmp_path / 'plan.md'
    plan.write_text('the plan\n')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    out = str(tmp_path / 'raw.npz')
    sel = E.part_dev(cpath, out)
    raw = X.load_raw(out)
    for pi, p in enumerate(X.DEV):                                # the stored scores are the scores of the planted catalogues, in the order of the registries
        for ti, tag in enumerate(PHU + CS):
            expected = X._pack('x', _score_planted(cats, truth, QUAL[(pi + 3 * ti) % 4]))
            for field in ('tot', 'blk22', 'n_true', 'n_det', 'completeness'):
                assert np.array_equal(raw['DEV|x%d_y%d|%s|%s' % (p[0], p[1], tag, field)], expected['x|' + field], equal_nan=True)
            assert raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)][0] == 1.5 + ti
    # the selection: recomputed independently from the planted catalogues
    def tot_of(pi, ti):
        return X._pack('x', _score_planted(cats, truth, QUAL[(pi + 3 * ti) % 4]))['x|tot']
    for fam, tags in G.FAMILIES.items():
        pooled = {tag: X.f_measure(tot_of(0, ti) + tot_of(1, ti))[0] for ti, tag in enumerate(PHU + CS) if tag in tags}
        best = max(pooled, key=lambda t: (round(pooled[t], 6), t))
        assert sel[fam] == best == str(raw['SEL|' + fam][0])
    assert str(raw['SEL|plan_sha256'][0]) == hashlib.sha256(b'the plan\n').hexdigest() and str(raw['SEL|dev_hash'][0]) == E.dev_hash(raw)
    assert E.verify_selection(raw, str(plan)) == sel
    with pytest.raises(RuntimeError, match='never made twice'):
        E.part_dev(cpath, out)


def test_part_dev_refuses_a_family_without_a_variant_scored_on_both_patches(monkeypatch, tmp_path):
    truth, cats, cpath = _planted_world(monkeypatch, tmp_path)
    plan = tmp_path / 'plan.md'
    plan.write_text('p\n')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    z = {k: v for k, v in np.load(cpath).items() if not any(k.startswith('x%d_y%d|%s|' % (X.DEV[1][0], X.DEV[1][1], t)) for t in PHU)}      # no photutils variant on the second patch
    np.savez_compressed(cpath, **z)
    with pytest.raises(RuntimeError, match='no variant'):
        E.part_dev(cpath, str(tmp_path / 'raw.npz'))


def test_part_dev_leaves_a_variant_without_a_catalogue_out_of_the_scores(monkeypatch, tmp_path):
    truth, cats, cpath = _planted_world(monkeypatch, tmp_path)
    plan = tmp_path / 'plan.md'
    plan.write_text('p\n')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    z = dict(np.load(cpath))
    drop = 'x%d_y%d|%s|' % (X.DEV[1][0], X.DEV[1][1], PHU[0])
    z = {k: v for k, v in z.items() if not k.startswith(drop)}
    np.savez_compressed(cpath, **z)
    out = str(tmp_path / 'raw.npz')
    E.part_dev(cpath, out)
    raw = X.load_raw(out)
    assert 'DEV|x%d_y%d|%s|tot' % (X.DEV[0][0], X.DEV[0][1], PHU[0]) in raw and 'DEV|x%d_y%d|%s|tot' % (X.DEV[1][0], X.DEV[1][1], PHU[0]) not in raw
    assert str(raw['SEL|photutils'][0]) != PHU[0]                                         # not eligible


def _selected_raw(monkeypatch, tmp_path):
    truth, cats, cpath = _planted_world(monkeypatch, tmp_path)
    plan = tmp_path / 'plan.md'
    plan.write_text('the plan\n')
    monkeypatch.setattr(E, 'PLAN', str(plan))
    out = str(tmp_path / 'raw.npz')
    sel = E.part_dev(cpath, out)
    return truth, cats, cpath, out, sel


def test_part_test_scores_the_selected_variants_and_the_references_on_the_requested_sets_only(monkeypatch, tmp_path):
    truth, cats, cpath, out, sel = _selected_raw(monkeypatch, tmp_path)
    E.part_test(['m2', 'm92'], cpath, out)
    raw = X.load_raw(out)
    tags = E.test_tags(sel)
    plist = [p for s in ('dev', 'm2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3') for p in G.SETS[s][1]]
    ti_of = {t: i for i, t in enumerate(PHU + CS)}
    for name, patches in (('m2', G.SETS['m2'][1]), ('m92', G.SETS['m92'][1])):
        for p in patches:
            pi = plist.index(p)
            for row, tag in tags.items():
                pre = 'T|%s|x%d_y%d|%s|' % (name, p[0], p[1], row)
                exp = X._pack('x', _score_planted(cats, truth, QUAL[(pi + 3 * ti_of[tag]) % 4]))
                assert np.array_equal(raw[pre + 'tot'], exp['x|tot']) and np.array_equal(raw[pre + 'blk22'], exp['x|blk22']) and str(raw[pre + 'name'][0]) == tag
                assert raw[pre + 'secs'][0] == 1.5 + ti_of[tag]
    assert not any(k.startswith('T|m53|') or k.startswith('T|m2r') for k in raw)           # nothing outside the requested sets
    assert sorted({k.split('|')[2] for k in raw if k.startswith('T|m2|')}) == sorted('x%d_y%d' % p for p in G.SETS['m2'][1])
    assert sorted({k.split('|')[3] for k in raw if k.startswith('T|m2|')}) == sorted(E.NEW_ROWS)


def test_part_test_never_scores_a_patch_twice_and_checks_the_selection_first(monkeypatch, tmp_path):
    truth, cats, cpath, out, sel = _selected_raw(monkeypatch, tmp_path)
    E.part_test(['m92'], cpath, out)
    raw = X.load_raw(out)
    k = 'T|m92|x1100_y220|phu_sel|tot'
    raw[k] = raw[k] + 1000                                                      # a stored value is never recomputed
    np.savez_compressed(out, **raw)
    E.part_test(['m92'], cpath, out)
    assert np.array_equal(X.load_raw(out)[k], raw[k])
    # a changed plan stops the test stage
    open(E.PLAN, 'w').write('edited\n')
    with pytest.raises(RuntimeError, match='plan file changed'):
        E.part_test(['m53'], cpath, out)
    assert not any(k.startswith('T|m53|') for k in X.load_raw(out))


def test_part_test_scores_the_repeat_frames_only_when_the_flag_exists(monkeypatch, tmp_path):
    truth, cats, cpath, out, sel = _selected_raw(monkeypatch, tmp_path)
    with pytest.raises(RuntimeError, match='flag'):
        E.part_test(['m2r1'], cpath, out)
    with pytest.raises(RuntimeError, match='flag'):
        E.part_test(['m2r1'], cpath, out, flag=str(tmp_path / 'absent.flag'))
    flag = tmp_path / 'done.flag'
    flag.write_text('')
    E.part_test(['m2r1'], cpath, out, flag=str(flag))
    raw = X.load_raw(out)
    assert sorted({k.split('|')[2] for k in raw if k.startswith('T|m2r1|')}) == ['x316_y124', 'x316_y224']
    with pytest.raises(RuntimeError, match='unknown set'):
        E.part_test(['m2r9'], cpath, out, flag=str(flag))


def test_part_dev_asks_for_the_m2_truth_and_part_test_for_the_truth_of_the_frame_of_each_set(monkeypatch, tmp_path):
    truth, cats, cpath, out, sel = _selected_raw(monkeypatch, tmp_path)
    assert TRUTH_CALLS == [('X', 'm2')] and TFP_CALLS == list(X.DEV)                 # the truth of each development patch is cut at its own origin
    flag = tmp_path / 'done.flag'
    flag.write_text('')
    TRUTH_CALLS.clear()
    TFP_CALLS.clear()
    E.part_test(['m2', 'm53', 'm92', 'm2r1', 'm2r3'], cpath, out, flag=str(flag))
    assert TRUTH_CALLS == [('X', 'm2'), ('X', 'm53'), ('X', 'm92'), ('R', 'm2r1'), ('R', 'm2r3')]
    assert TFP_CALLS == [tuple(p) for s_ in ('m2', 'm53', 'm92', 'm2r1', 'm2r3') for p in G.SETS[s_][1]]


def test_part_old_scores_the_default_and_the_tuned_old_baselines_by_the_names_of_the_freeze(monkeypatch, tmp_path):
    truth, cats, cpath = _planted_world(monkeypatch, tmp_path)
    names = {'photutils': 'bic', X.cs_name((0.0, 0.0, 0.0, np.inf)): 'good', 'photutils_ext_all_local_i6': 'poor', 'crowdsource|snr5|ff0|qf0|spinf': 'mid', 'photutils_local_t3_i3': 'good'}
    seen = []
    monkeypatch.setattr(X, 'baseline_scores', lambda p, t, base=None: (seen.append(p), {n: _score_planted(cats, truth, q) for n, q in names.items()})[1])
    fz = tmp_path / 'freeze.json'
    fz.write_text('{"tuned_baselines": {"photutils": "photutils_ext_all_local_i6", "crowdsource": "crowdsource|snr5|ff0|qf0|spinf"}, "config": {"sigma_b": 3.0}}')
    out = str(tmp_path / 'raw.npz')
    E.part_old(out, freeze=str(fz))
    raw = X.load_raw(out)
    assert seen == list(X.DEV)
    for p in X.DEV:
        for label, q in (('photutils (default)', 'bic'), ('crowdsource (default)', 'good'), ('photutils (tuned)', 'poor'), ('crowdsource (tuned)', 'mid')):
            assert np.array_equal(raw['OLD|x%d_y%d|%s|tot' % (p[0], p[1], label)], X._pack('x', _score_planted(cats, truth, q))['x|tot'])
    assert str(raw['OLD|photutils (tuned)|name'][0]) == 'photutils_ext_all_local_i6' and str(raw['OLD|crowdsource (default)|name'][0]) == X.cs_name((0.0, 0.0, 0.0, np.inf))
    assert str(raw['OLD|crowdsource (tuned)|name'][0]) == 'crowdsource|snr5|ff0|qf0|spinf' and str(raw['OLD|photutils (default)|name'][0]) == 'photutils'


def test_describe_names_the_settings_of_each_variant():
    assert E.describe('phu_hp5_new_i3') == 'photutils, high-pass 5 px, mode new, 3 iterations'
    assert E.describe('phu_hp3_all_i6') == 'photutils, high-pass 3 px, mode all, 6 iterations, grouper 4.4 px, relaxed finder cuts'
    assert E.describe('cs_default') == 'crowdsource, raw image, the default weights, threshold 5'
    assert E.describe('cs_w_t8') == 'crowdsource, raw image, the ISS weights, threshold 8'
    assert E.describe('cs_hp3_t3') == 'crowdsource, image and stamp high-passed 3 px, the ISS weights, threshold 3'
    with pytest.raises(KeyError):
        E.describe('nope')
    for tag in PHU + CS:
        E.describe(tag)


def test_the_development_table_row_of_a_variant_is_computed_from_its_pooled_counts():
    raw = {}
    for i, tag in enumerate(PHU + CS):
        for j, p in enumerate(X.DEV):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] = tot(100 + i, 30 + i + j, 120 + 2 * i, 40 + i, 400, 150 + 3 * i + j, 500 + 5 * i)
            raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)] = np.array([10.0 * (i + 1) + j])
    raw['SEL|photutils'], raw['SEL|crowdsource'] = np.array([PHU[2]]), np.array([CS[3]])
    for label in ('photutils (default)', 'crowdsource (tuned)'):
        for p in X.DEV:
            raw['OLD|x%d_y%d|%s|tot' % (p[0], p[1], label)] = tot(100, 30, 120, 40, 400, 150, 500)
        raw['OLD|%s|name' % label] = np.array(['photutils' if label.startswith('photutils') else 'crowdsource|snr5|ff0|qf0|spinf'])      # the names of the old baselines: the crowdsource ones contain pipes
    table_lines = E.dev_table(raw).split('\n')
    assert {l.count('|') for l in table_lines} == {table_lines[0].count('|')}                         # every row has the cells of the header (a pipe in a name would add some)
    rows = {l.split('|')[2].strip(): [c.strip() for c in l.strip().strip('|').split('|')] for l in table_lines[2:]}
    for i, tag in enumerate(PHU + CS):
        s = sum((raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] for p in X.DEV))
        f, c, pu = X.f_measure(s)
        mark = '**' if tag in (PHU[2], CS[3]) else ''
        cells = rows['%s%s%s%s' % (mark, tag, mark, ' (reference)' if tag in (G.PHOT_REF,) + G.CS_REFS else '')]
        assert cells[0] == ('photutils' if tag in PHU else 'crowdsource') and cells[2] == E.describe(tag)
        assert cells[3:] == ['%.3f' % c, '%.3f' % pu, '%s%.3f%s' % (mark, f, mark), '%d' % s[2], '%d' % (s[6] - s[5]), '%.0f' % (20.0 * (i + 1) + 1)]
    assert rows['photutils (default)'][0] == 'old family' and rows['photutils (default)'][2] == 'photutils' and rows['photutils (default)'][3] == '%.3f' % (60 / 200) and rows['photutils (default)'][-1] == '-'
    assert rows['crowdsource (tuned)'][:3] == ['old family', 'crowdsource (tuned)', 'crowdsource/snr5/ff0/qf0/spinf'] and str(raw['OLD|crowdsource (tuned)|name'][0]) == 'crowdsource|snr5|ff0|qf0|spinf'       # printed with slashes, stored as it is
    assert 'crowdsource (default)' not in rows and 'photutils (tuned)' not in rows                    # an old row appears only where it was scored on both development patches
    raw['OLD|x%d_y%d|photutils (tuned)|tot' % (X.DEV[0][0], X.DEV[0][1])] = tot(100, 30, 120, 40, 400, 150, 500)         # ... on one patch only: not shown, and no failure
    assert 'photutils (tuned)' not in {l.split('|')[2].strip() for l in E.dev_table(raw).split('\n')[2:]}


def test_the_patch_table_cells_are_the_intervals_of_the_stratified_bootstrap_of_the_blocks_and_only_the_declared_comparisons_are_listed():
    truth, cats = _synthetic_catalogs(np.random.default_rng(12))
    fz, sc = _planted_frozen(cats, truth)
    raw = {}
    _planted_new(raw, 'm53', (760, 590), cats, truth, dict(phu_sel='good', phu_ref='mid', cs_sel='poor', cs_hp='poor', cs_default='good'))
    _planted_new(raw, 'm53', (640, 680), cats, truth, dict(phu_sel='mid', phu_ref='mid', cs_sel='bic', cs_hp='poor', cs_default='good'))
    md = E.report(raw, fz)
    sect = md.split('### Patch (640, 680)')[0].split('### Patch (760, 590)')[1]
    blk = np.concatenate([X._pack('x', sc[q])['x|blk22'] for q in ['good']])
    (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk, X.boot_idx_strat([25]))
    tot_ = X._pack('x', sc['good'])['x|tot']
    row = [l for l in sect.split('\n') if l.startswith('| tuned ISS, isotropic kernel (as frozen) |')]
    assert len(row) == 1
    cells = [x.strip() for x in row[0].strip().strip('|').split('|')]
    assert cells[1:] == ['%d / %d' % (tot_[1], tot_[0]), X.fmt_iv(c, cci), '%d / %d' % (tot_[3], tot_[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci), '%.3f (%d)' % E.purity_all(tot_)]
    names = [l.split('|')[1].strip() for l in sect.split('Paired differences')[1].split('\n') if l.startswith('| ') and 'minus' in l and not l.startswith('| comparison')]
    assert names == ['tuned ISS minus photutils HP selected', 'tuned ISS minus crowdsource selected', 'tuned ISS (kernel from the stamp) minus photutils HP selected',
                     'tuned ISS (kernel from the stamp) minus crowdsource selected', 'photutils HP selected minus photutils old tuned', 'crowdsource selected minus crowdsource old tuned']
    pooled = md.split('### All 2 test patches together')[1].split('\n## ')[0]
    blk2 = np.concatenate([X._pack('x', sc['good'])['x|blk22']] * 2)
    (c2, cci2), (pu2, pci2), (ff2, fci2) = X.boot_cpf(blk2, X.boot_idx_strat([25, 25]))
    rowp = [l for l in pooled.split('\n') if l.startswith('| tuned ISS, isotropic kernel (as frozen) |')]
    assert rowp and X.fmt_iv(ff2, fci2) in rowp[0] and X.fmt_iv(c2, cci2) in rowp[0] and X.fmt_iv(pu2, pci2) in rowp[0]
    pnames = [l.split('|')[1].strip() for l in pooled.split('Paired differences')[1].split('\n') if l.startswith('| ') and 'minus' in l and not l.startswith('| comparison')]
    assert len(pnames) == 10 and 'tuned ISS (kernel from the stamp)' not in ' '.join(pnames) and pnames[0] == 'tuned ISS minus photutils HP selected' and pnames[-1] == 'crowdsource selected minus crowdsource HP reference'


# ---- the frozen rows (read only) and the groups ------------------------------------------------------------------------------------------

def _put(d, prefix, score):
    for k, v in X._pack(prefix, score).items():
        d[k] = v


def _planted_frozen(cats, truth, kappa=2.0):
    """Frozen raws with a distinct catalogue quality in every source so that a row read from the wrong key is caught."""
    sc = {q: _score_planted(cats, truth, q) for q in QUAL}
    m2, psf, rep = {}, {}, {}
    for pre in ['T|x310_y630|', 'T|x310_y730|', 'T2|m53|x760_y590|', 'T2|m53|x640_y680|', 'T2|m92|x1100_y220|']:
        m2[pre + 'cfg'] = np.array([3.0, 1.2, kappa])
        _put(m2, pre + 'ours|k%g' % kappa, sc['good'])
        _put(m2, pre + 'ours|k5', sc['bic'])                                                 # the default's row: not the tuned one
        for label, q in (('photutils (tuned)', 'poor'), ('crowdsource (tuned)', 'mid'), ('photutils (default)', 'bic'), ('crowdsource (default)', 'good')):
            _put(m2, pre + label, sc[q])
        m2[pre + 'photutils (tuned)|name'] = np.array(['photutils_ext_all_local_i6'])
        m2[pre + 'crowdsource (tuned)|name'] = np.array(['crowdsource|snr5|ff0|qf0|spinf'])
    for pre in ['P|m53|x760_y590|']:
        psf[pre + 'cfg'] = np.array([3.0, 1.2, kappa])
        _put(psf, pre + 'ours|k%g' % kappa, sc['mid'])
    for pre in ['R|m2r1|x316_y124|', 'R|m2r1|x316_y224|', 'R|m2r2|x133_y281|', 'R|m2r3|x370_y441|']:
        rep[pre + 'cfg'] = np.array([3.0, 1.2, kappa])
        _put(rep, pre + 'iso|ours|k%g' % kappa, sc['good'])
        _put(rep, pre + 'adq|ours|k%g' % kappa, sc['mid'])
        for label, q in (('photutils (tuned)', 'poor'), ('crowdsource (tuned)', 'poor'), ('photutils (default)', 'bic'), ('crowdsource (default)', 'bic')):
            _put(rep, pre + label, sc[q])
        rep[pre + 'photutils (tuned)|name'] = np.array(['photutils_ext_all_local_i6'])
        rep[pre + 'crowdsource (tuned)|name'] = np.array(['crowdsource|snr5|ff0|qf0|spinf'])
    return dict(m2=m2, psf=psf, repeat=rep, sha=dict(m2='a' * 64, psf='b' * 64, repeat=None)), sc


def test_frozen_rows_read_each_cluster_from_its_own_keys_and_the_reported_significance(monkeypatch):
    truth, cats = _synthetic_catalogs(np.random.default_rng(4))
    for kappa in (2.0, 3.5):
        fz, sc = _planted_frozen(cats, truth, kappa)
        tot_of = {q: X._pack('x', sc[q])['x|tot'] for q in QUAL}
        for name, pxy, want in (('m2', (310, 630), {'iss': 'good', 'phu_old_t': 'poor', 'cs_old_t': 'mid', 'phu_old_d': 'bic', 'cs_old_d': 'good'}),
                                ('m53', (760, 590), {'iss': 'good', 'iss_adq': 'mid', 'phu_old_t': 'poor', 'cs_old_t': 'mid', 'phu_old_d': 'bic', 'cs_old_d': 'good'}),
                                ('m92', (1100, 220), {'iss': 'good', 'phu_old_t': 'poor', 'cs_old_t': 'mid', 'phu_old_d': 'bic', 'cs_old_d': 'good'}),
                                ('m2r1', (316, 124), {'iss': 'good', 'iss_adq': 'mid', 'phu_old_t': 'poor', 'cs_old_t': 'poor', 'phu_old_d': 'bic', 'cs_old_d': 'bic'})):
            rows = E.frozen_rows(name, pxy, fz)
            assert [r[0] for r in rows] == list(want), name
            for row, label, src, pre in rows:
                assert np.array_equal(src[pre + '|tot'], tot_of[want[row]]), (name, row)
            labels = {r[0]: r[1] for r in rows}
            assert 'isotropic kernel' in labels['iss'] and ('[photutils_ext_all_local_i6]' in labels['phu_old_t']) and ('[crowdsource/snr5/ff0/qf0/spinf]' in labels['cs_old_t']) and '|' not in ''.join(labels.values())
            assert '[' not in labels['phu_old_d']
    fz, _ = _planted_frozen(cats, truth)
    assert E.frozen_rows('m2', (310, 630), dict(fz, m2=None)) == [] and E.frozen_rows('m2r1', (316, 124), dict(fz, repeat=None)) == []
    assert E.frozen_rows('m2', (360, 1040), fz) == [] and E.frozen_rows('m2r2', (195, 689), fz) == []      # a patch the frozen raw does not have
    nopsf = E.frozen_rows('m53', (760, 590), dict(fz, psf=None))
    assert 'iss_adq' not in [r[0] for r in nopsf] and 'iss' in [r[0] for r in nopsf]
    assert 'iss_adq' not in [r[0] for r in E.frozen_rows('m53', (640, 680), fz)]                 # the post hoc kernel run has no row for this patch


def _planted_new(raw, name, pxy, cats, truth, qual):
    for row, q in qual.items():
        _put(raw, 'T|%s|x%d_y%d|%s' % (name, pxy[0], pxy[1], row), _score_planted(cats, truth, q))
        raw['T|%s|x%d_y%d|%s|name' % (name, pxy[0], pxy[1], row)] = np.array(['tag_' + row])


def test_group_data_pools_the_patches_in_the_raw_and_only_the_rows_present_on_all_of_them():
    truth, cats = _synthetic_catalogs(np.random.default_rng(6))
    fz, sc = _planted_frozen(cats, truth)
    raw = {}
    _planted_new(raw, 'm2', (310, 630), cats, truth, dict(phu_sel='good', phu_ref='mid', cs_sel='poor', cs_hp='mid', cs_default='good'))
    _planted_new(raw, 'm2', (310, 730), cats, truth, dict(phu_sel='good', phu_ref='mid', cs_sel='poor', cs_hp='mid', cs_default='good'))
    g = E.group_data(raw, [('m2', G.SETS['m2'][1])], fz)
    assert g['done'] == [('m2', (310, 630)), ('m2', (310, 730))] and g['sizes'] == [25, 25]
    assert g['rows'] == ['iss', 'phu_old_t', 'cs_old_t', 'phu_old_d', 'cs_old_d', 'phu_sel', 'phu_ref', 'cs_sel', 'cs_hp', 'cs_default']           # the table order
    one = X._pack('x', sc['good'])
    assert g['blk']['iss'].shape == (50, 4) and np.array_equal(g['blk']['iss'], np.concatenate([one['x|blk22']] * 2)) and np.array_equal(g['tot']['iss'], 2 * one['x|tot'])
    assert np.array_equal(g['tot']['phu_old_t'], 2 * X._pack('x', sc['poor'])['x|tot'])
    fz['m2']['T|x310_y730|photutils (tuned)|name'] = np.array(['other_name'])                   # the label of a pooled row is that of the first patch
    assert '[photutils_ext_all_local_i6]' in E.group_data(raw, [('m2', G.SETS['m2'][1])], fz)['label']['phu_old_t']
    del fz['m2']['T|x310_y730|crowdsource (tuned)|tot']                                         # missing on one patch: dropped from the pooled rows
    g2 = E.group_data(raw, [('m2', G.SETS['m2'][1])], fz)
    assert 'cs_old_t' not in g2['rows'] and 'phu_old_t' in g2['rows']
    assert E.group_data(raw, [('m53', G.SETS['m53'][1])], fz)['done'] == []
    gp = E.group_data(raw, [('m2', [(310, 730)])], fz)
    assert gp['done'] == [('m2', (310, 730))] and gp['sizes'] == [25]


def test_a_difference_row_shows_completeness_purity_and_f_each_with_its_own_interval_and_the_verdict_of_f():
    # catalogue a: complete but impure; b: pure but less complete; per block the counts are fixed, so the intervals are points
    a = np.tile([10.0, 8.0, 40.0, 8.0], (25, 1))
    b = np.tile([10.0, 4.0, 8.0, 4.0], (25, 1))
    g = dict(blk=dict(iss=a, phu_sel=b), rows=['iss', 'phu_sel'], sizes=[25])
    idx = X.boot_idx_strat([25])
    tbl, dfs = E.diff_table((('iss', 'phu_sel'), ('phu_sel', 'iss'), ('iss', 'cs_sel')), g, idx)
    dc, dp, df = X.boot_diff(a, b, idx)
    assert dc[0] > 0 and dp[0] < 0 and df[0] < 0                                                    # a wins on completeness, loses on purity and on F
    rows = [[c_.strip() for c_ in l.strip().strip('|').split('|')] for l in tbl.split('\n')[2:]]
    assert len(rows) == 2                                                                           # a pair with a missing row is skipped
    assert rows[0] == ['tuned ISS minus photutils HP selected', X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), 'loss']
    assert rows[1][4] == 'win' and rows[1][3].startswith('+') is False and X.fmt_iv(*X.boot_diff(b, a, idx)[2]) == rows[1][3]
    assert E.diff_table((('iss', 'cs_sel'),), g, idx) == ('', {})


def test_purity_over_all_magnitudes_and_the_unmatched_detections():
    p, un = E.purity_all(np.array([100, 60, 80, 55, 300, 140, 400.0]))
    assert p == 0.35 and un == 260 and isinstance(un, int)


def test_paired_differences_are_the_first_minus_the_second_with_the_verdict_of_the_interval():
    truth, cats = _synthetic_catalogs(np.random.default_rng(8))
    fz, sc = _planted_frozen(cats, truth)
    raw = {}
    _planted_new(raw, 'm2', (310, 630), cats, truth, dict(phu_sel='good', phu_ref='mid', cs_sel='bic', cs_hp='poor', cs_default='good'))
    g = E.group_data(raw, [('m2', [(310, 630)])], fz)
    idx = X.boot_idx_strat(g['sizes'])
    tbl, dfs = E.diff_table(E.COMPARISONS_POOLED, g, idx)
    # tuned ISS is 'good'; phu_sel is 'good' too: a tie; against cs_hp = poor a win; phu_sel (good) minus the old tuned photutils (poor) a win; cs_sel (bic, better than mid) minus the old tuned crowdsource (mid) a win
    verdict = {k: X.verdict(*v) for k, v in dfs.items()}
    assert dfs[('iss', 'phu_sel')][0] == 0.0 and np.all(dfs[('iss', 'phu_sel')][1] == 0.0)                                       # the same catalogue block for block: exactly no difference
    assert verdict[('iss', 'phu_sel')] == 'tie' and verdict[('iss', 'cs_hp')] == 'win' and verdict[('phu_sel', 'phu_old_t')] == 'win' and verdict[('cs_sel', 'cs_old_t')] == 'win'
    assert verdict[('iss', 'cs_sel')] != 'loss' and verdict[('phu_ref', 'phu_old_d')] != 'win' and ('iss_adq', 'phu_sel') not in dfs               # mid is not better than bic; no kernel-from-the-stamp row on M2
    d = dfs[('iss', 'cs_hp')][0]
    assert d == X.cpf(g['blk']['iss'].sum(0))[2] - X.cpf(g['blk']['cs_hp'].sum(0))[2] > 0                    # first minus second, on F
    rows = [l for l in tbl.split('\n') if l.startswith('| tuned ISS minus crowdsource HP reference |')]
    assert len(rows) == 1 and rows[0].rstrip().endswith('| win |')
    dfs_rev = E.diff_table((('cs_hp', 'iss'),), g, idx)[1][('cs_hp', 'iss')]
    assert abs(dfs_rev[0] + d) < 1e-12 and X.verdict(*dfs_rev) == 'loss'


# ---- the report -------------------------------------------------------------------------------------------------------------------------

def _planted_report_inputs(monkeypatch, seed=11):
    rng = np.random.default_rng(seed)
    truth, cats = _synthetic_catalogs(rng)
    fz, sc = _planted_frozen(cats, truth)
    raw = dict(_dev_raw(PHU[3], CS[5]))
    for p in X.DEV:
        for t in PHU + CS:
            raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], t)] = np.array([10.0])
    raw['SEL|photutils'], raw['SEL|crowdsource'], raw['SEL|plan_sha256'] = np.array([PHU[3]]), np.array([CS[5]]), np.array(['c' * 64])
    for label in ('photutils (default)', 'photutils (tuned)', 'crowdsource (default)', 'crowdsource (tuned)'):
        for p in X.DEV:
            _put(raw, 'OLD|x%d_y%d|%s' % (p[0], p[1], label), sc['mid'])
        raw['OLD|%s|name' % label] = np.array([{'photutils (default)': 'photutils', 'photutils (tuned)': 'photutils_ext_all_local_i6', 'crowdsource (default)': 'crowdsource|snr0|ff0|qf0|spinf',
                                                'crowdsource (tuned)': 'crowdsource|snr5|ff0|qf0|spinf'}[label]])
    quals = dict(phu_sel='good', phu_ref='mid', cs_sel='poor', cs_hp='poor', cs_default='good')
    for name, pxy in (('m2', (310, 630)), ('m2', (310, 730)), ('m53', (760, 590)), ('m53', (640, 680)), ('m92', (1100, 220)), ('m2r1', (316, 124)), ('m2r1', (316, 224)), ('m2r2', (133, 281)), ('m2r3', (370, 441))):
        _planted_new(raw, name, pxy, cats, truth, quals)
    return raw, fz, truth, cats


def test_the_report_is_deterministic_and_reads_the_raw_and_the_frozen_rows(monkeypatch):
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    md = E.report(raw, fz)
    assert md == E.report(raw, fz)
    for needle in ('## Development patches (110, 630), (210, 630): every variant', '| photutils | **%s** (reference)' % PHU[3] if PHU[3] == G.PHOT_REF else '| photutils | **%s**' % PHU[3],
                   'Selected: photutils `%s`, crowdsource `%s`.' % (PHU[3], CS[5]), '| old family | photutils (tuned) | photutils_ext_all_local_i6 |', '| old family | crowdsource (default) | crowdsource/snr0/ff0/qf0/spinf |', 'SHA-256 of the plan stored with the selection: ' + 'c' * 64,
                   'Frozen results read: `m2_raw.npz` ' + 'a' * 16 + '...; `psf_raw.npz` ' + 'b' * 16 + '...; `repeat_raw.npz` absent.', '## M2 (NGC 7089), SDSS run 2583', '### Patch (310, 630)', '### Patch (310, 730)',
                   '### All 2 test patches together', '## M53 (NGC 5024), SDSS run 5390 camcol 6 field 162', '### Patch (760, 590)', '## M92 (NGC 6341), SDSS run 4682 camcol 6 field 208', '### Patch (1100, 220)',
                   '## M2 repeat 1 (NGC 7089)', '## All the repeat patches together', '## Summary: pooled F at mag < 22 by cluster', 'tuned ISS, isotropic kernel (as frozen)',
                   'tuned ISS, kernel chosen from the stamp (post hoc)', 'tuned ISS, kernel chosen from the stamp |', 'photutils, old tuned [photutils_ext_all_local_i6]',
                   'crowdsource, old tuned [crowdsource/snr5/ff0/qf0/spinf]', 'photutils on the high-passed patch, selected [tag_phu_sel]', 'crowdsource, default (re-run)',
                   'purity over all magnitudes (unmatched detections)', 'tuned ISS minus photutils HP selected', 'photutils HP selected minus photutils old tuned', 'verdict on F'):
        assert needle in md, needle
    assert '(1430, 550)' not in md and '(110, 501)' not in md and 'nan' not in md.lower().replace('nan,', '')               # absent patches are not invented
    # the pooled M2 row is read from the new arrays of both patches
    t = raw['T|m2|x310_y630|phu_sel|tot'] + raw['T|m2|x310_y730|phu_sel|tot']
    pa, un = E.purity_all(t)
    rows = [l for l in md.split('\n') if l.startswith('| photutils on the high-passed patch, selected [tag_phu_sel] | %d / %d |' % (t[1], t[0]))]
    assert rows and '%.3f (%d)' % (pa, un) in rows[0]
    assert sum(1 for l in md.split('\n') if l.startswith('### Patch (')) == 9
    assert '### All 1 test patches together' not in md and '### All 2 test patches together' in md                  # a cluster with one patch has no pooled block
    for name, pxy in (('m2', (310, 630)), ('m53', (640, 680)), ('m92', (1100, 220)), ('m2r3', (370, 441))):
        sect = md.split('### Patch (%d, %d)\n' % pxy)[1].split('\n')[1]
        assert sect == 'Truth sources brighter than 22 mag in the scored area: %d.' % raw['T|%s|x%d_y%d|phu_sel|tot' % (name, pxy[0], pxy[1])][0]


def test_shown_prints_a_tag_without_pipes_and_the_table_helper_refuses_a_pipe_in_a_cell():
    assert E.shown('crowdsource|snr5|ff0|qf0|spinf') == 'crowdsource/snr5/ff0/qf0/spinf' and E.shown(np.str_('a|b|c')) == 'a/b/c' and E.shown('plain') == 'plain'
    assert E.table(['a', 'b'], [['1', '2'], ['3', '4']]) == X.table(['a', 'b'], [['1', '2'], ['3', '4']]) == '| a | b |\n|---|---|\n| 1 | 2 |\n| 3 | 4 |'
    for hdr, rows in ((['a|b', 'c'], [['1', '2']]), (['a', 'b'], [['1', '2|3']]), (['a', 'b'], [['1', '2'], ['3', 'x|y']]), (['a', 'b'], [['p|q', '2']]), (['a', 'b'], [['1', '2'], ['3|', '4']])):
        with pytest.raises(ValueError, match='pipe'):
            E.table(hdr, rows)


def test_the_four_tables_of_the_report_refuse_a_pipe_in_a_cell(monkeypatch):
    """The guard is used by every table of the report: a label with a pipe, however it got there, stops the report instead of printing a table that does not render."""
    blk = np.tile([10.0, 8.0, 40.0, 8.0], (25, 1))
    totv = np.array([250.0, 200.0, 1000.0, 200.0, 250.0, 200.0, 1000.0])
    idx = X.boot_idx_strat([25])
    g = dict(rows=['iss'], label={'iss': 'a|b'}, tot={'iss': totv}, blk={'iss': blk}, sizes=[25])
    with pytest.raises(ValueError, match='pipe'):
        E.method_table(g, idx)
    g['label'] = {'iss': 'a/b'}
    assert E.method_table(g, idx).split('\n')[2].startswith('| a/b |')
    g2 = dict(rows=['iss', 'phu_sel'], blk=dict(iss=blk, phu_sel=blk), sizes=[25])
    assert E.diff_table((('iss', 'phu_sel'),), g2, idx)[0].split('\n')[2].startswith('| tuned ISS minus photutils HP selected |')
    monkeypatch.setitem(E.SHORT, 'iss', 'tuned|ISS')
    with pytest.raises(ValueError, match='pipe'):
        E.diff_table((('iss', 'phu_sel'),), g2, idx)
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    monkeypatch.setattr(E, 'describe', lambda tag: 'photutils|x')
    with pytest.raises(ValueError, match='pipe'):
        E.dev_table(raw)
    monkeypatch.undo()
    assert E.summary_table([('M2', dict(done=[1], tot={}), {}, {})]).split('\n')[-1].startswith('| M2 (1 patches) | - |')
    with pytest.raises(ValueError, match='pipe'):
        E.summary_table([('M|2', dict(done=[1], tot={}), {}, {})])


def test_every_table_of_the_report_has_the_cells_of_its_header_whatever_the_tags(monkeypatch):
    """The tags of the old crowdsource baselines contain pipes: printed in a cell they would split it (the 25 rows of the first report). Every table row has as many cells as its header."""
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    raw['T|m2|x310_y630|phu_sel|name'] = np.array(['tag|with|pipes'])                                  # a pipe in the tag of a new row is printed as a slash too
    md = E.report(raw, fz)
    lines = md.split('\n')
    n_tables, i = 0, 0
    while i < len(lines):
        if lines[i].startswith('|'):
            j = i
            while j < len(lines) and lines[j].startswith('|'):
                j += 1
            block = lines[i:j]
            n_tables += 1
            assert block[1].replace('|', '').replace('-', '') == ''                                    # the separator line
            assert {l.count('|') for l in block} == {block[0].count('|')}, block[0][:80]
            i = j
        else:
            i += 1
    assert n_tables == 28                                                                              # the development table, the summary, two tables for each of 9 patches, and the pooled tables of 4 groups
    assert 'crowdsource|snr' not in md and 'tag|with|pipes' not in md and 'photutils on the high-passed patch, selected [tag/with/pipes]' in md
    assert md.count('crowdsource, old tuned [crowdsource/snr5/ff0/qf0/spinf]') == md.count('| crowdsource, old tuned [') > 10          # every table that has the row prints the tag the same way


def test_the_summary_table_has_one_row_per_pooled_group_with_the_signed_differences(monkeypatch):
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    md = E.report(raw, fz)
    summ = md.split('## Summary: pooled F at mag < 22 by cluster')[1].split('\n## ')[0]
    lines = [l for l in summ.split('\n') if l.startswith('| ') and 'patches)' in l]
    heads = [l.split('|')[1].strip() for l in lines]
    assert heads == ['M2 (NGC 7089) (2 patches)', 'M53 (NGC 5024) (2 patches)', 'M92 (NGC 6341) (1 patches)', 'M2 repeat 1 (NGC 7089) (2 patches)', 'M2 repeat 2 (NGC 7089) (1 patches)',
                     'M2 repeat 3 (NGC 7089) (1 patches)', 'All the repeat patches (4 patches)']
    assert all(len(l.split('|')) == 12 for l in lines)

    def expect(rows, sizes):
        """The cells of a summary row from the planted arrays alone: rows = {row: [(tot, blk) per patch]}."""
        blk = {r: np.concatenate([b for t_, b in v]) for r, v in rows.items()}
        tot_ = {r: sum(t_ for t_, b in v) for r, v in rows.items()}
        idx = X.boot_idx_strat(sizes)

        def cell(r):
            f, ci = X.boot_cpf(blk[r], idx)[2]
            return '%.3f [%.3f, %.3f] (%d)' % (f, ci[0], ci[1], tot_[r][6] - tot_[r][5])

        def vd(a, b):
            df = X.boot_diff(blk[a], blk[b], idx)[2]
            return '%+.3f %s' % (df[0], X.verdict(*df))
        return [cell('iss'), cell('phu_old_t'), cell('cs_old_t'), cell('phu_sel'), cell('cs_sel'), vd('iss', 'phu_sel'), vd('iss', 'cs_sel'), vd('phu_sel', 'phu_old_t'), vd('cs_sel', 'cs_old_t')]

    def pack(patches):
        rows = {}
        for src, pre, row in (('m2', 'T|x%d_y%d|ours|k2|', 'iss'), ('m2', 'T|x%d_y%d|photutils (tuned)|', 'phu_old_t'), ('m2', 'T|x%d_y%d|crowdsource (tuned)|', 'cs_old_t')):
            rows[row] = [(fz[src][pre % p + 'tot'], fz[src][pre % p + 'blk22']) for p in patches]
        for row in ('phu_sel', 'cs_sel'):
            rows[row] = [(raw['T|m2|x%d_y%d|%s|tot' % (p[0], p[1], row)], raw['T|m2|x%d_y%d|%s|blk22' % (p[0], p[1], row)]) for p in patches]
        return rows
    m2 = [c.strip() for c in lines[0].split('|')][2:11]
    assert m2 == expect(pack([(310, 630), (310, 730)]), [25, 25])
    assert m2[5].startswith('+0.000') and m2[5].endswith('tie') and re.fullmatch(r'\+0\.\d{3} win', m2[6]) and re.fullmatch(r'\+0\.\d{3} win', m2[7]) and re.fullmatch(r'-0\.\d{3} (loss|tie)', m2[8]), m2
    # the repeat group: the four patches of three frames, each row taken from its own frozen keys
    rep = [('m2r1', (316, 124)), ('m2r1', (316, 224)), ('m2r2', (133, 281)), ('m2r3', (370, 441))]
    rows = {}
    for row, key in (('iss', 'iso|ours|k2|'), ('phu_old_t', 'photutils (tuned)|'), ('cs_old_t', 'crowdsource (tuned)|')):
        rows[row] = [(fz['repeat']['R|%s|x%d_y%d|%s' % (n, p[0], p[1], key) + 'tot'], fz['repeat']['R|%s|x%d_y%d|%s' % (n, p[0], p[1], key) + 'blk22']) for n, p in rep]
    for row in ('phu_sel', 'cs_sel'):
        rows[row] = [(raw['T|%s|x%d_y%d|%s|tot' % (n, p[0], p[1], row)], raw['T|%s|x%d_y%d|%s|blk22' % (n, p[0], p[1], row)]) for n, p in rep]
    assert [c.strip() for c in lines[-1].split('|')][2:11] == expect(rows, [25] * 4)


def test_the_report_contains_only_what_the_raw_and_the_frozen_rows_hold(monkeypatch):
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    md = E.report({k: v for k, v in raw.items() if not k.startswith('T|m2r') and not k.startswith('OLD|')}, dict(fz, repeat=None))
    assert 'M2 repeat' not in md and 'All the repeat patches' not in md and '| old family |' not in md and '## M2 (NGC 7089)' in md
    two = E.report({k: v for k, v in raw.items() if not (k.startswith('T|m2r2|') or k.startswith('T|m2r3|'))}, fz)           # two repeat patches (one frame): pooled over the repeat patches too
    assert '## All the repeat patches together' in two and 'All the repeat patches (2 patches)' in two and 'M2 repeat 2' not in two
    one = E.report({k: v for k, v in raw.items() if not (k.startswith('T|m2r2|') or k.startswith('T|m2r3|') or k.startswith('T|m2r1|x316_y224'))}, fz)     # a single repeat patch: no pooled group of repeat patches
    assert 'M2 repeat 1 (NGC 7089) (1 patches)' in one and 'All the repeat patches' not in one
    md2 = E.report({k: v for k, v in raw.items() if not k.startswith('T|')}, fz)
    assert '## Development patches' in md2 and '### Patch' not in md2 and 'Summary' not in md2


def test_from_raw_regenerates_the_report_byte_for_byte(monkeypatch, tmp_path):
    raw, fz, truth, cats = _planted_report_inputs(monkeypatch)
    outdir = tmp_path / 'out'
    outdir.mkdir()
    np.savez_compressed(outdir / 'baselines_hp_raw.npz', **raw)
    paths = {}
    for k, name in (('m2', 'm2_raw.npz'), ('psf', 'psf_raw.npz'), ('repeat', 'repeat_raw.npz')):
        np.savez_compressed(tmp_path / name, **fz[k])
        paths[k] = str(tmp_path / name)
    env = dict(os.environ, HP_OUT=str(outdir), HP_M2_RAW=paths['m2'], HP_PSF_RAW=paths['psf'], HP_REPEAT_RAW=paths['repeat'])
    cmd = [sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_astro_hp.py'), '--from-raw']
    first = None
    for _ in range(2):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        data = (outdir / 'baselines_hp.md').read_bytes()
        assert first is None or data == first
        first = data
        (outdir / 'baselines_hp.md').unlink()
    expected = E.report(X.load_raw(str(outdir / 'baselines_hp_raw.npz')), E.load_frozen(paths))
    assert first == expected.encode()
    line = [l for l in first.decode().split('\n') if l.startswith('Frozen results read:')][0]
    assert line == 'Frozen results read: %s.' % '; '.join('`%s` %s...' % (n, hashlib.sha256(open(paths[k], 'rb').read()).hexdigest()[:16]) for k, n in (('m2', 'm2_raw.npz'), ('psf', 'psf_raw.npz'), ('repeat', 'repeat_raw.npz')))
    assert len({hashlib.sha256(open(paths[k], 'rb').read()).hexdigest() for k in paths}) == 3


def test_the_command_line_runs_the_stage_asked_for_with_the_declared_defaults(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(E, 'part_dev', lambda: calls.append(('dev',)) or {'photutils': 'a', 'crowdsource': 'b'})
    monkeypatch.setattr(E, 'part_old', lambda: calls.append(('old',)))
    monkeypatch.setattr(E, 'part_test', lambda sets, flag=None: calls.append(('test', list(sets), flag)))
    for argv in (['--dev'], ['--old'], ['--test'], ['--test', '--sets', 'm92', 'm2r1', '--flag', 'f.flag']):
        monkeypatch.setattr(sys, 'argv', ['exp_astro_hp.py'] + argv)
        E.main()
    assert calls == [('dev',), ('old',), ('test', ['m2', 'm53', 'm92'], None), ('test', ['m92', 'm2r1'], 'f.flag')]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_hp.py'])
    with pytest.raises(SystemExit):
        E.main()
    out = tmp_path / 'r.md'
    monkeypatch.setattr(E, 'MD', str(out))
    monkeypatch.setattr(E, 'RAW', str(tmp_path / 'raw.npz'))
    np.savez_compressed(tmp_path / 'raw.npz', a=np.zeros(1))
    monkeypatch.setattr(E, 'report', lambda raw: 'the report of %s\n' % sorted(raw))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_hp.py', '--from-raw'])
    E.main()
    assert out.read_text() == "the report of ['a']\n"
