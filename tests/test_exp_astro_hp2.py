"""Round 2 of the baselines given the ISS's front end, scoring, selection and report (experiments/exp_astro_hp2.py): the selection with the round-1 selection as a candidate, its freeze (plan, development scores, the
round-1 selection), the scoring of the development and test patches from planted catalogues, the report with the round-1 selections and the frozen rows beside the round-2 ones, the verdict counts, and the byte-identity of
`--from-raw`. Round 1's helpers (planted worlds, planted frozen results) are those of tests/test_exp_astro_hp.py."""
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
import exp_astro_baselines_hp as G1     # noqa: E402
import exp_astro_baselines_hp2 as G2    # noqa: E402
import exp_astro_hp as E1               # noqa: E402
import exp_astro_hp2 as H2              # noqa: E402
import exp_astro_m2 as X                # noqa: E402
from test_exp_astro_baselines_hp2 import _raw2   # noqa: E402
from test_exp_astro_hp import (CS, PHU, QUAL, TFP_CALLS, TRUTH_CALLS, _planted_frozen, _planted_new, _planted_report_inputs, _planted_world, _score_planted, tot)   # noqa: E402

SEL1 = ('phu_hp3_all_i6', 'cs_hp3_t5')
P2, C2 = G2.FAMILIES2['photutils'], G2.FAMILIES2['crowdsource']


def _dev_totals(raw, tags, vals, patches=None):
    for t in tags:
        for p in (patches or X.DEV):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*vals)


# ---- the selection rule --------------------------------------------------------------------------------------------------------------

def _base_raw():
    raw = {}
    _dev_totals(raw, P2 + [SEL1[0]] + C2 + [SEL1[1]], (100, 20, 100, 20))
    return raw


def test_choose_hp2_takes_the_best_candidate_of_each_family_round_1_selection_included():
    sel1 = dict(photutils=SEL1[0], crowdsource=SEL1[1])
    assert H2.candidates(sel1) == {'photutils': P2 + [SEL1[0]], 'crowdsource': C2 + [SEL1[1]]}
    for tp, tc in ((P2[3], C2[5]), (P2[0], C2[0]), (SEL1[0], SEL1[1]), (P2[5], SEL1[1]), (SEL1[0], C2[7])):
        raw = _base_raw()
        _dev_totals(raw, [tp], (100, 60, 100, 70))
        _dev_totals(raw, [tc], (100, 50, 100, 55))
        assert H2.choose_hp2(raw, sel1) == {'photutils': tp, 'crowdsource': tc}


def test_choose_hp2_keeps_the_round_1_variant_on_a_tie_and_when_no_round_2_variant_is_better():
    sel1 = dict(photutils=SEL1[0], crowdsource=SEL1[1])
    raw = _base_raw()
    _dev_totals(raw, [P2[2], SEL1[0]], (100, 60, 100, 70))                          # a round-2 variant and the round-1 variant with identical totals: the larger tag, the round-1 one
    _dev_totals(raw, [C2[6], SEL1[1]], (100, 50, 100, 55))
    assert H2.choose_hp2(raw, sel1) == {'photutils': SEL1[0], 'crowdsource': SEL1[1]}
    _dev_totals(raw, [P2[2]], (100, 60, 100, 71))                                   # strictly better by 1/100 in purity: round 2 wins
    assert H2.choose_hp2(raw, sel1)['photutils'] == P2[2]
    raw = _base_raw()
    _dev_totals(raw, [SEL1[0]], (100, 10, 100, 10))                                  # the round-1 variant worse than the rest: the best of round 2 (ties among the rest to the larger tag)
    assert H2.choose_hp2(raw, sel1)['photutils'] == max(P2)
    # F rounded to six decimals: a tie is a tie
    raw = _base_raw()
    for p in X.DEV:
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], P2[4])] = tot(10 ** 9, 600000000, 10 ** 9, 700000001)
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], SEL1[0])] = tot(10 ** 9, 600000000, 10 ** 9, 700000000)
    assert E1.dev_F(raw, P2[4]) > E1.dev_F(raw, SEL1[0]) and H2.choose_hp2(raw, sel1)['photutils'] == SEL1[0]
    for p in X.DEV:                                                                    # better by 4e-6 in F: different at the sixth decimal, equal at the third
        raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], P2[4])] = tot(10 ** 9, 600000000, 10 ** 9, 700010000)
    assert 1e-6 < E1.dev_F(raw, P2[4]) - E1.dev_F(raw, SEL1[0]) < 1e-5 and round(E1.dev_F(raw, P2[4]), 3) == round(E1.dev_F(raw, SEL1[0]), 3)
    assert H2.choose_hp2(raw, sel1)['photutils'] == P2[4]


def test_choose_hp2_pools_both_patches_and_needs_both_scores_and_the_families_are_separate():
    sel1 = dict(photutils=SEL1[0], crowdsource=SEL1[1])
    raw = _base_raw()
    for t, vals in ((P2[0], [(100, 60, 100, 60), (100, 10, 100, 10)]), (P2[1], [(100, 40, 100, 40), (100, 40, 100, 40)])):
        for p, v in zip(X.DEV, vals):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*v)
    assert H2.choose_hp2(raw, sel1)['photutils'] == P2[1]
    for t, vals in ((P2[0], [(100, 10, 100, 10), (100, 60, 100, 60)]), (P2[1], [(100, 40, 100, 40), (100, 40, 100, 40)])):
        for p, v in zip(X.DEV, vals):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], t)] = tot(*v)
    assert H2.choose_hp2(raw, sel1)['photutils'] == P2[1]
    raw['DEV|x%d_y%d|%s|tot' % (X.DEV[0][0], X.DEV[0][1], P2[2])] = tot(100, 100, 100, 100)             # scored on one patch only: not eligible
    raw.pop('DEV|x%d_y%d|%s|tot' % (X.DEV[1][0], X.DEV[1][1], P2[2]))
    assert H2.choose_hp2(raw, sel1)['photutils'] == P2[1]
    raw2 = _base_raw()
    _dev_totals(raw2, [C2[3]], (100, 95, 100, 95))                                                       # a crowdsource variant much better than every photutils one
    assert H2.choose_hp2(raw2, sel1) == {'photutils': SEL1[0], 'crowdsource': C2[3]}
    assert 'crowdsource' not in H2.choose_hp2({k: v for k, v in _base_raw().items() if not any('|%s|' % t in k for t in C2 + [SEL1[1]])}, sel1)


# ---- the freeze ------------------------------------------------------------------------------------------------------------------------

def test_verify_selection2_returns_the_stored_selection_and_refuses_every_change(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, P2[1], C2[4])
    assert H2.verify_selection2(raw) == {'photutils': P2[1], 'crowdsource': C2[4]}
    for key in ('SEL|photutils', 'SEL|crowdsource', 'SEL|r1_photutils', 'SEL|r1_crowdsource', 'SEL|plan_sha256', 'SEL|dev_hash'):
        with pytest.raises(RuntimeError, match='no selection'):
            H2.verify_selection2({k: v for k, v in raw.items() if k != key})
    wrong = dict(raw)
    wrong['SEL|photutils'] = np.array([P2[0]])
    with pytest.raises(RuntimeError, match='not what the rule gives'):
        H2.verify_selection2(wrong)
    quiet = dict(raw)                                                                                    # a change that leaves the selection as it is: the hash of the development scores
    quiet['DEV|x%d_y%d|%s|tot' % (X.DEV[0][0], X.DEV[0][1], P2[5])] = tot(100, 21, 100, 20)
    assert H2.choose_hp2(quiet, dict(photutils=SEL1[0], crowdsource=SEL1[1])) == H2.choose_hp2(raw, dict(photutils=SEL1[0], crowdsource=SEL1[1]))
    with pytest.raises(RuntimeError, match='development scores changed'):
        H2.verify_selection2(quiet)
    other = dict(raw)                                                                                    # made with another round-1 selection
    other['SEL|r1_crowdsource'] = np.array(['cs_hp3_t3'])
    with pytest.raises(RuntimeError, match='round-1 selection'):
        H2.verify_selection2(other)
    open(H2.PLAN2, 'w').write('edited\n')
    with pytest.raises(RuntimeError, match='plan file changed'):
        H2.verify_selection2(raw)


def test_verify_selection2_uses_the_plan_and_the_round_1_raw_it_is_given(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, P2[1], C2[4])
    raw1 = X.load_raw(H2.RAW1)
    other_plan = tmp_path / 'other_plan.md'
    other_plan.write_text('plan two\n')
    open(H2.PLAN2, 'w').write('the default plan was edited\n')
    assert H2.verify_selection2(raw, plan=str(other_plan)) == {'photutils': P2[1], 'crowdsource': C2[4]}      # the plan given, not the default
    with pytest.raises(RuntimeError, match='plan file changed'):
        H2.verify_selection2(raw)
    monkeypatch.setattr(H2, 'RAW1', str(tmp_path / 'nonexistent_raw1.npz'))
    assert H2.verify_selection2(raw, raw1=raw1, plan=str(other_plan)) == {'photutils': P2[1], 'crowdsource': C2[4]}     # the round-1 raw given, not the file
    assert H2.round1_selection(raw1) == {'photutils': SEL1[0], 'crowdsource': SEL1[1]}
    with pytest.raises(Exception):
        H2.verify_selection2(raw, plan=str(other_plan))


def test_verify_selection2_refuses_when_the_round_1_selection_itself_no_longer_verifies(monkeypatch, tmp_path):
    rawp, raw = _raw2(tmp_path, monkeypatch, P2[1], C2[4])
    open(E1.PLAN, 'w').write('the round-1 plan, edited\n')
    with pytest.raises(RuntimeError, match='plan file changed'):
        H2.verify_selection2(raw)


# ---- scoring the planted catalogues ----------------------------------------------------------------------------------------------------

def _world2(monkeypatch, tmp_path):
    """Planted truth and catalogues: the round-1 file (made by round 1's own scoring, so that its selection exists) and a round-2 file with a planted catalogue of one of four qualities for every (patch, tag)."""
    truth, cats, cpath1 = _planted_world(monkeypatch, tmp_path)
    plan1 = tmp_path / 'plan1.md'
    plan1.write_text('plan one\n')
    monkeypatch.setattr(E1, 'PLAN', str(plan1))
    raw1p = str(tmp_path / 'raw1.npz')
    sel1 = E1.part_dev(cpath1, raw1p)
    plan2 = tmp_path / 'plan2.md'
    plan2.write_text('plan two\n')
    monkeypatch.setattr(H2, 'PLAN2', str(plan2))
    monkeypatch.setattr(H2, 'RAW1', raw1p)
    out = {}
    plist = [p for s in ('dev', 'm2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3') for p in G1.SETS[s][1]]
    for pi, p in enumerate(plist):
        for ti, tag in enumerate(P2 + C2):
            q = cats[QUAL[(pi + 2 * ti + 1) % 4]]
            k = 'x%d_y%d|%s|' % (p[0], p[1], tag)
            out[k + 'x'], out[k + 'y'], out[k + 'flux'], out[k + 'secs'] = q[0], q[1], q[2], np.array([7.5 + ti])
    cpath2 = str(tmp_path / 'cats2.npz')
    np.savez_compressed(cpath2, **out)
    return truth, cats, cpath1, cpath2, raw1p, sel1, plist


def test_part_dev2_scores_the_round_2_variants_and_the_round_1_selection_and_selects_by_the_rule(monkeypatch, tmp_path):
    truth, cats, cpath1, cpath2, raw1p, sel1, plist = _world2(monkeypatch, tmp_path)
    out = str(tmp_path / 'raw2.npz')
    TFP_CALLS.clear()
    sel = H2.part_dev2(cpath2, cpath1, raw1p, out)
    assert TFP_CALLS == list(X.DEV)
    raw = X.load_raw(out)
    for pi, p in enumerate(X.DEV):
        for ti, tag in enumerate(P2 + C2):
            exp = X._pack('x', _score_planted(cats, truth, QUAL[(pi + 2 * ti + 1) % 4]))
            for f in ('tot', 'blk22', 'n_true', 'completeness'):
                assert np.array_equal(raw['DEV|x%d_y%d|%s|%s' % (p[0], p[1], tag, f)], exp['x|' + f], equal_nan=True)
            assert raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)][0] == 7.5 + ti
        for tag in sel1.values():                                                                         # the round-1 selected variants, from the round-1 catalogues
            ti1 = (PHU + CS).index(tag)
            assert np.array_equal(raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)], X._pack('x', _score_planted(cats, truth, QUAL[(pi + 3 * ti1) % 4]))['x|tot'])
    def tot_of(pi, tag):
        if tag in P2 + C2:
            return X._pack('x', _score_planted(cats, truth, QUAL[(pi + 2 * (P2 + C2).index(tag) + 1) % 4]))['x|tot']
        return X._pack('x', _score_planted(cats, truth, QUAL[(pi + 3 * (PHU + CS).index(tag)) % 4]))['x|tot']
    for fam, tags in H2.candidates(sel1).items():
        pooled = {t: X.f_measure(tot_of(0, t) + tot_of(1, t))[0] for t in tags}
        assert sel[fam] == max(pooled, key=lambda t: (round(pooled[t], 6), t)) == str(raw['SEL|' + fam][0])
    assert str(raw['SEL|r1_photutils'][0]) == sel1['photutils'] and str(raw['SEL|r1_crowdsource'][0]) == sel1['crowdsource']
    assert str(raw['SEL|plan_sha256'][0]) == hashlib.sha256(b'plan two\n').hexdigest() and str(raw['SEL|dev_hash'][0]) == E1.dev_hash(raw)
    assert H2.verify_selection2(raw) == sel
    with pytest.raises(RuntimeError, match='never made twice'):
        H2.part_dev2(cpath2, cpath1, raw1p, out)


def test_part_dev2_refuses_a_round_1_variant_that_does_not_reproduce_its_round_1_totals(monkeypatch, tmp_path):
    truth, cats, cpath1, cpath2, raw1p, sel1, plist = _world2(monkeypatch, tmp_path)
    z = dict(np.load(cpath1))
    k = 'x%d_y%d|%s|flux' % (X.DEV[1][0], X.DEV[1][1], sel1['crowdsource'])
    z[k] = z[k] * 0.1                                                                                     # the catalogue file is no longer the one that was scored in round 1
    np.savez_compressed(cpath1, **z)
    with pytest.raises(RuntimeError, match='does not reproduce'):
        H2.part_dev2(cpath2, cpath1, raw1p, str(tmp_path / 'raw2.npz'))
    assert not (tmp_path / 'raw2.npz').exists()


def test_part_dev2_refuses_a_family_without_a_candidate_on_both_patches(monkeypatch, tmp_path):
    truth, cats, cpath1, cpath2, raw1p, sel1, plist = _world2(monkeypatch, tmp_path)
    z = {k: v for k, v in np.load(cpath2).items() if not any(k.startswith('x%d_y%d|%s|' % (X.DEV[1][0], X.DEV[1][1], t)) for t in P2)}
    np.savez_compressed(cpath2, **z)
    out = str(tmp_path / 'raw2.npz')
    sel = H2.part_dev2(cpath2, cpath1, raw1p, out)                                                         # the round-1 variant is still a candidate: the family is not empty
    assert sel['photutils'] == sel1['photutils']


def test_part_test2_scores_the_two_selected_variants_only_and_never_twice(monkeypatch, tmp_path):
    truth, cats, cpath1, cpath2, raw1p, sel1, plist = _world2(monkeypatch, tmp_path)
    out = str(tmp_path / 'raw2.npz')
    sel = H2.part_dev2(cpath2, cpath1, raw1p, out)
    if sel['photutils'] in G1.FAMILIES['photutils'] or sel['crowdsource'] in G1.FAMILIES['crowdsource']:              # carried variants live in the round-1 file
        z = dict(np.load(cpath2))
        z.update(dict(np.load(cpath1)))
        np.savez_compressed(cpath2, **z)
    TRUTH_CALLS.clear()
    TFP_CALLS.clear()
    H2.part_test2(['m2', 'm92'], cpath2, out)
    assert TFP_CALLS == [tuple(p) for s_ in ('m2', 'm92') for p in G1.SETS[s_][1]]
    raw = X.load_raw(out)
    tags = {'phu2': sel['photutils'], 'cs2': sel['crowdsource']}
    allt = P2 + C2
    for name in ('m2', 'm92'):
        for p in G1.SETS[name][1]:
            pi = plist.index(p)
            for row, tag in tags.items():
                q = QUAL[(pi + 2 * allt.index(tag) + 1) % 4] if tag in allt else QUAL[(pi + 3 * (PHU + CS).index(tag)) % 4]
                pre = 'T|%s|x%d_y%d|%s|' % (name, p[0], p[1], row)
                assert np.array_equal(raw[pre + 'tot'], X._pack('x', _score_planted(cats, truth, q))['x|tot']) and str(raw[pre + 'name'][0]) == tag
                assert raw[pre + 'secs'][0] == (7.5 + allt.index(tag) if tag in allt else 1.5 + (PHU + CS).index(tag))
    assert sorted({k.split('|')[3] for k in raw if k.startswith('T|')}) == ['cs2', 'phu2'] and not any(k.startswith('T|m53|') for k in raw)
    assert TRUTH_CALLS == [('X', 'm2'), ('X', 'm92')]
    k = 'T|m92|x1100_y220|phu2|tot'
    raw[k] = raw[k] + 1000
    np.savez_compressed(out, **raw)
    H2.part_test2(['m92'], cpath2, out)                                                                    # a stored result is not recomputed
    assert np.array_equal(X.load_raw(out)[k], raw[k])
    open(H2.PLAN2, 'w').write('edited\n')
    with pytest.raises(RuntimeError, match='plan file changed'):
        H2.part_test2(['m53'], cpath2, out)


def test_part_test2_scores_the_repeat_frames_only_with_the_flag_and_refuses_unknown_sets(monkeypatch, tmp_path):
    truth, cats, cpath1, cpath2, raw1p, sel1, plist = _world2(monkeypatch, tmp_path)
    out = str(tmp_path / 'raw2.npz')
    sel = H2.part_dev2(cpath2, cpath1, raw1p, out)
    z = dict(np.load(cpath2))
    z.update(dict(np.load(cpath1)))
    np.savez_compressed(cpath2, **z)
    for flag in (None, str(tmp_path / 'absent.flag')):
        with pytest.raises(RuntimeError, match='flag'):
            H2.part_test2(['m2r1'], cpath2, out, flag=flag)
    f = tmp_path / 'done.flag'
    f.write_text('')
    TRUTH_CALLS.clear()
    H2.part_test2(['m2r1'], cpath2, out, flag=str(f))
    assert TRUTH_CALLS == [('R', 'm2r1')] and sorted({k.split('|')[2] for k in X.load_raw(out) if k.startswith('T|m2r1|')}) == ['x316_y124', 'x316_y224']
    with pytest.raises(RuntimeError, match='unknown set'):
        H2.part_test2(['m2r9'], cpath2, out, flag=str(f))


# ---- the report ----------------------------------------------------------------------------------------------------------------------

def _inputs2(monkeypatch):
    """The round-1 raw and the frozen results of round 1's planted report, and a round-2 raw: development totals, a selection, the two selected rows on nine patches."""
    raw1, fz, truth, cats = _planted_report_inputs(monkeypatch)
    raw = {}
    _dev_totals(raw, P2 + C2 + list(SEL1), (100, 20, 100, 20))
    _dev_totals(raw, [P2[2]], (100, 60, 100, 70))
    _dev_totals(raw, [C2[6]], (100, 50, 100, 55))
    for p in X.DEV:
        for t in P2 + C2 + list(SEL1):
            raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], t)] = np.array([12.0])
    raw['SEL|photutils'], raw['SEL|crowdsource'] = np.array([P2[2]]), np.array([C2[6]])
    raw['SEL|r1_photutils'], raw['SEL|r1_crowdsource'] = np.array([SEL1[0]]), np.array([SEL1[1]])
    raw['SEL|plan_sha256'] = np.array(['d' * 64])
    patches = (('m2', (310, 630)), ('m2', (310, 730)), ('m53', (760, 590)), ('m53', (640, 680)), ('m92', (1100, 220)), ('m2r1', (316, 124)), ('m2r1', (316, 224)), ('m2r2', (133, 281)), ('m2r3', (370, 441)))
    for name, pxy in patches:
        _planted_new(raw, name, pxy, cats, truth, dict(phu2='mid', cs2='good'))
    return raw, raw1, fz, truth, cats, patches


def test_the_report_is_deterministic_and_shows_the_rows_of_both_rounds_and_the_frozen_rows(monkeypatch):
    raw, raw1, fz, truth, cats, patches = _inputs2(monkeypatch)
    md = H2.report2(raw, fz, raw1)
    assert md == H2.report2(raw, fz, raw1)
    for needle in ('# The baselines given the ISS\'s front end, round 2', 'SHA-256 of the plan stored with the selection: ' + 'd' * 64, 'On the repeat frames these baselines are added later',
                   'Frozen results read: `m2_raw.npz` ' + 'a' * 16 + '...', '## Development patches (110, 630), (210, 630): every candidate',
                   'Selected: photutils `%s` (round 1: `%s`), crowdsource `%s` (round 1: `%s`).' % (P2[2], SEL1[0], C2[6], SEL1[1]), '| photutils | **%s** |' % P2[2], '| crowdsource | **%s** |' % C2[6],
                   '| photutils | %s (round 1 selection) |' % SEL1[0], '| crowdsource | %s (round 1 selection) |' % SEL1[1], '| old family | photutils (tuned) | photutils_ext_all_local_i6 |',
                   '## M2 (NGC 7089), SDSS run 2583', '### Patch (310, 630)', '### All 2 test patches together', '## M53 (NGC 5024)', '### Patch (760, 590)', '## All the repeat patches together',
                   '## Summary: pooled F at mag < 22 by cluster', '## Verdict counts over the test patches', 'tuned ISS, isotropic kernel (as frozen)', 'tuned ISS, kernel chosen from the stamp (post hoc)',
                   'photutils, round 1 selection [tag_phu_sel]', 'crowdsource, round 1 selection [tag_cs_sel]', 'photutils, round 2 selection [tag_phu2]', 'crowdsource, round 2 selection [tag_cs2]',
                   'crowdsource, old tuned [crowdsource/snr5/ff0/qf0/spinf]', 'tuned ISS minus photutils round 2', 'tuned ISS (kernel from the stamp) minus crowdsource round 2',
                   'photutils round 2 minus photutils round 1', 'crowdsource round 2 minus crowdsource old tuned', 'tuned ISS minus crowdsource round 1', 'Paired differences (first minus second, same blocks), 95 % intervals:'):
        assert needle in md, needle
    assert '95 %%' not in md and '(1430, 550)' not in md and 'nan' not in md.lower().replace('nan,', '') and 'photutils, round 1 selection' in md
    assert sum(1 for l in md.split('\n') if l.startswith('### Patch (')) == 9
    t = raw['T|m2|x310_y630|cs2|tot'] + raw['T|m2|x310_y730|cs2|tot']
    pa, un = E1.purity_all(t)
    rows = [l for l in md.split('\n') if l.startswith('| crowdsource, round 2 selection [tag_cs2] | %d / %d |' % (t[1], t[0]))]
    assert rows and '%.3f (%d)' % (pa, un) in rows[0]


def test_the_rows_print_the_tags_without_pipes_and_the_patch_tables_use_the_round_1_bootstrap(monkeypatch):
    raw, raw1, fz, truth, cats, patches = _inputs2(monkeypatch)
    raw['T|m2|x310_y630|phu2|name'] = np.array(['r2|tag'])
    raw1['T|m2|x310_y630|cs_sel|name'] = np.array(['r1|tag'])
    md = H2.report2(raw, fz, raw1)
    assert 'photutils, round 2 selection [r2/tag]' in md and 'crowdsource, round 1 selection [r1/tag]' in md and 'r2|tag' not in md and 'r1|tag' not in md
    sect = md.split('### Patch (310, 630)\n')[1].split('### Patch (310, 730)')[0]
    assert sect.split('\n')[1] == 'Truth sources brighter than 22 mag in the scored area: %d.' % raw['T|m2|x310_y630|phu2|tot'][0]
    row = [l for l in sect.split('\n') if l.startswith('| tuned ISS, isotropic kernel (as frozen) |')][0]
    (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(fz['m2']['T|x310_y630|ours|k2|blk22'], X.boot_idx_strat([25]))
    assert [x.strip() for x in row.strip().strip('|').split('|')][2:6] == [X.fmt_iv(c, cci), '%d / %d' % tuple(fz['m2']['T|x310_y630|ours|k2|tot'][[3, 2]]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)]
    seen = []
    orig = X.boot_diff
    monkeypatch.setattr(X, 'boot_diff', lambda a, b, idx: (seen.append(idx), orig(a, b, idx))[1])
    H2.patch_verdicts(raw, raw1, fz)
    assert seen and all(np.array_equal(i, X.boot_idx_strat([25])) for i in seen)


def test_the_check_of_the_round_1_variant_is_generated_from_the_development_totals(monkeypatch):
    raw, raw1, fz, *_ = _inputs2(monkeypatch)
    assert C2[6] == 'cs2_hp3_t5_pd'
    assert '(`cs2_hp3_t5_pd` against `cs_hp3_t5`, development totals of both patches): DIFFERENT.' in H2.report2(raw, fz, raw1)       # here its totals differ from those of the round-1 variant
    raw = dict(raw)
    for p in X.DEV:
        raw['DEV|x%d_y%d|cs2_hp3_t5_pd|tot' % p] = raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], SEL1[1])]
    assert '(`cs2_hp3_t5_pd` against `cs_hp3_t5`, development totals of both patches): identical.' in H2.report2(raw, fz, raw1)
    raw['DEV|x%d_y%d|cs2_hp3_t5_pd|tot' % X.DEV[1]] = raw['DEV|x%d_y%d|cs2_hp3_t5_pd|tot' % X.DEV[1]] + 1                          # the second patch alone differs
    assert 'development totals of both patches): DIFFERENT.' in H2.report2(raw, fz, raw1)
    assert 'Check:' not in H2.report2({k: v for k, v in raw.items() if not k.startswith('SEL|')}, fz, raw1)


def test_every_table_of_the_report_has_the_cells_of_its_header(monkeypatch):
    raw, raw1, fz, *_ = _inputs2(monkeypatch)
    md = H2.report2(raw, fz, raw1)
    lines = md.split('\n')
    n, i = 0, 0
    while i < len(lines):
        if lines[i].startswith('|'):
            j = i
            while j < len(lines) and lines[j].startswith('|'):
                j += 1
            n += 1
            assert block_ok(lines[i:j])
            i = j
        else:
            i += 1
    assert n == 1 + 2 + 1 + 2 * 9 + 2 * 4                                                                    # dev table, two summary tables, the verdict table, two tables per patch, two per pooled group
    g = dict(blk=dict(iss=np.ones((25, 4)), phu2=np.ones((25, 4))))
    monkeypatch.setitem(H2.SHORT2, 'iss', 'a|b')
    with pytest.raises(ValueError, match='pipe'):
        H2.diff_table2((('iss', 'phu2'),), g, X.boot_idx_strat([25]))


def block_ok(block):
    return block[1].replace('|', '').replace('-', '') == '' and {l.count('|') for l in block} == {block[0].count('|')}


def test_the_summary_tables_and_the_verdict_counts_are_computed_from_the_planted_arrays(monkeypatch):
    raw, raw1, fz, truth, cats, patches = _inputs2(monkeypatch)
    md = H2.report2(raw, fz, raw1)
    summ = md.split('## Summary: pooled F at mag < 22 by cluster')[1].split('## Verdict counts')[0]
    t1 = [l for l in summ.split('\n') if l.startswith('| ') and 'patches)' in l]
    assert [l.split('|')[1].strip() for l in t1] == ['M2 (NGC 7089) (2 patches)', 'M53 (NGC 5024) (2 patches)', 'M92 (NGC 6341) (1 patches)', 'M2 repeat 1 (NGC 7089) (2 patches)',
                                                       'M2 repeat 2 (NGC 7089) (1 patches)', 'M2 repeat 3 (NGC 7089) (1 patches)', 'All the repeat patches (4 patches)']
    # M2 row, from the arrays
    pts = [(310, 630), (310, 730)]
    blk = {'iss': [fz['m2']['T|x%d_y%d|ours|k2|blk22' % p] for p in pts], 'phu_old_t': [fz['m2']['T|x%d_y%d|photutils (tuned)|blk22' % p] for p in pts],
           'cs_old_t': [fz['m2']['T|x%d_y%d|crowdsource (tuned)|blk22' % p] for p in pts], 'phu1': [raw1['T|m2|x%d_y%d|phu_sel|blk22' % p] for p in pts], 'cs1': [raw1['T|m2|x%d_y%d|cs_sel|blk22' % p] for p in pts],
           'phu2': [raw['T|m2|x%d_y%d|phu2|blk22' % p] for p in pts], 'cs2': [raw['T|m2|x%d_y%d|cs2|blk22' % p] for p in pts]}
    tots = {'iss': [fz['m2']['T|x%d_y%d|ours|k2|tot' % p] for p in pts], 'phu_old_t': [fz['m2']['T|x%d_y%d|photutils (tuned)|tot' % p] for p in pts],
            'cs_old_t': [fz['m2']['T|x%d_y%d|crowdsource (tuned)|tot' % p] for p in pts], 'phu1': [raw1['T|m2|x%d_y%d|phu_sel|tot' % p] for p in pts],
            'cs1': [raw1['T|m2|x%d_y%d|cs_sel|tot' % p] for p in pts], 'phu2': [raw['T|m2|x%d_y%d|phu2|tot' % p] for p in pts], 'cs2': [raw['T|m2|x%d_y%d|cs2|tot' % p] for p in pts]}
    idx = X.boot_idx_strat([25, 25])
    B = {r: np.concatenate(v) for r, v in blk.items()}
    cells = []
    for r in ('iss', None, 'phu_old_t', 'cs_old_t', 'phu1', 'cs1', 'phu2', 'cs2'):
        if r is None:
            cells.append('-')                                                                              # no stamp-kernel row on M2
            continue
        f, ci = X.boot_cpf(B[r], idx)[2]
        T = sum(tots[r])
        cells.append('%.3f [%.3f, %.3f] (%d)' % (f, ci[0], ci[1], T[6] - T[5]))
    assert [c.strip() for c in t1[0].split('|')][2:10] == cells
    t2 = [l for l in summ.split('Pooled paired differences')[1].split('\n') if l.startswith('| M2 (NGC 7089) |')][0]
    def vd(a, b):
        df = X.boot_diff(B[a], B[b], idx)[2]
        return '%+.3f %s' % (df[0], X.verdict(*df))
    assert [c.strip() for c in t2.split('|')][2:10] == [vd('iss', 'phu2'), vd('iss', 'cs2'), '-', '-', vd('phu2', 'phu1'), vd('cs2', 'cs1'), vd('phu2', 'phu_old_t'), vd('cs2', 'cs_old_t')]
    # the verdict counts: every patch on its own blocks, the ISS row of the patch
    vtab = md.split('## Verdict counts over the test patches')[1].split('\n## ')[0]
    got = {l.split('|')[1].strip(): [c.strip() for c in l.split('|')[2:-1]] for l in vtab.split('\n') if l.startswith('| ') and ' minus ' in l.split('|')[1]}
    def iss_blocks(name, p, kind):
        if name == 'm2':
            return fz['m2']['T|x%d_y%d|ours|k2|blk22' % p]
        if name in ('m53', 'm92'):
            if kind == 'stamp' and name == 'm53' and p == (760, 590):
                return fz['psf']['P|m53|x760_y590|ours|k2|blk22']
            return fz['m2']['T2|%s|x%d_y%d|ours|k2|blk22' % (name, p[0], p[1])]
        return fz['repeat']['R|%s|x%d_y%d|%s|ours|k2|blk22' % (name, p[0], p[1], 'iso' if kind == 'iso' else 'adq')]
    new_blk = {'phu2': lambda n, p: raw['T|%s|x%d_y%d|phu2|blk22' % (n, p[0], p[1])], 'cs2': lambda n, p: raw['T|%s|x%d_y%d|cs2|blk22' % (n, p[0], p[1])],
               'phu1': lambda n, p: raw1['T|%s|x%d_y%d|phu_sel|blk22' % (n, p[0], p[1])], 'cs1': lambda n, p: raw1['T|%s|x%d_y%d|cs_sel|blk22' % (n, p[0], p[1])]}
    groups = ['m2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3']
    assert len(got) == 8
    for kind, kname in (('iso', 'tuned ISS, isotropic kernel as frozen'), ('stamp', 'tuned ISS, kernel from the stamp (M2: the isotropic one)')):
        for base in ('phu2', 'cs2', 'phu1', 'cs1'):
            per = {g: [X.verdict(*X.boot_diff(iss_blocks(n, p, kind), new_blk[base](n, p), X.boot_idx_strat([25]))[2]) for n, p in patches if n == g] for g in groups}
            allv = [v for g in groups for v in per[g]]
            exp = ['%d / %d / %d (%d patches)' % (tuple(allv.count(w) for w in ('win', 'tie', 'loss')) + (len(allv),))] + ['%d / %d / %d' % tuple(per[g].count(w) for w in ('win', 'tie', 'loss')) for g in groups]
            assert got['%s minus %s' % (kname, H2.SHORT2[base])] == exp, (kind, base)


def test_the_report_leaves_out_what_the_raw_and_the_frozen_rows_do_not_hold(monkeypatch):
    raw, raw1, fz, truth, cats, patches = _inputs2(monkeypatch)
    md = H2.report2({k: v for k, v in raw.items() if not k.startswith('T|m2r')}, dict(fz, repeat=None), raw1)
    assert 'M2 repeat' not in md and 'All the repeat patches' not in md and '## M2 (NGC 7089)' in md
    md = H2.report2({k: v for k, v in raw.items() if not (k.startswith('T|m2r2|') or k.startswith('T|m2r3|'))}, fz, raw1)
    assert '## All the repeat patches together' in md and 'All the repeat patches (2 patches)' in md
    md = H2.report2({k: v for k, v in raw.items() if not k.startswith('T|')}, fz, raw1)
    assert '## Development patches' in md and '### Patch' not in md and 'Summary' not in md and 'Verdict counts' not in md
    md = H2.report2(raw, fz, {k: v for k, v in raw1.items() if not k.startswith('OLD|')})
    assert '| old family |' not in md


def test_a_difference_row_shows_completeness_purity_and_f_each_with_its_interval_and_the_verdict_of_f():
    a = np.tile([10.0, 8.0, 40.0, 8.0], (25, 1))                                                          # complete but impure
    b = np.tile([10.0, 4.0, 8.0, 4.0], (25, 1))                                                           # pure but less complete
    g = dict(blk=dict(iss=a, phu2=b), rows=['iss', 'phu2'], sizes=[25])
    idx = X.boot_idx_strat([25])
    tbl, dfs = H2.diff_table2((('iss', 'phu2'), ('phu2', 'iss'), ('iss', 'cs2')), g, idx)
    dc, dp, df = X.boot_diff(a, b, idx)
    assert dc[0] > 0 and dp[0] < 0 and df[0] < 0
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in tbl.split('\n')[2:]]
    assert len(rows) == 2 and rows[0] == ['tuned ISS minus photutils round 2', X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), 'loss']
    assert rows[1][0] == 'photutils round 2 minus tuned ISS' and rows[1][4] == 'win' and rows[1][3] == X.fmt_iv(*X.boot_diff(b, a, idx)[2])
    assert H2.diff_table2((('iss', 'cs2'),), g, idx) == ('', {})


def test_the_development_table_row_of_a_candidate_is_computed_from_its_pooled_counts_and_old_rows_need_both_patches():
    raw, raw1 = {}, {}
    tags = P2 + [SEL1[0]] + C2 + [SEL1[1]]
    for i, tag in enumerate(tags):
        for j, p in enumerate(X.DEV):
            raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] = tot(100 + i, 30 + i + j, 120 + 2 * i, 40 + i, 400, 150 + 3 * i + j, 500 + 5 * i)
            raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)] = np.array([10.0 * (i + 1) + j])
    raw['SEL|photutils'], raw['SEL|crowdsource'] = np.array([P2[1]]), np.array([C2[2]])
    raw['SEL|r1_photutils'], raw['SEL|r1_crowdsource'] = np.array([SEL1[0]]), np.array([SEL1[1]])
    for label in ('photutils (default)', 'crowdsource (tuned)'):
        for p in X.DEV:
            raw1['OLD|x%d_y%d|%s|tot' % (p[0], p[1], label)] = tot(100, 30, 120, 40, 400, 150, 500)
        raw1['OLD|%s|name' % label] = np.array(['photutils' if label.startswith('photutils') else 'crowdsource|snr5|ff0|qf0|spinf'])
    raw1['OLD|x%d_y%d|photutils (tuned)|tot' % X.DEV[0]] = tot(100, 30, 120, 40, 400, 150, 500)            # one patch only
    raw1['OLD|photutils (tuned)|name'] = np.array(['x'])
    lines = H2.dev_table2(raw, raw1).split('\n')
    assert {l.count('|') for l in lines} == {lines[0].count('|')}
    rows = {l.split('|')[2].strip(): [c.strip() for c in l.strip().strip('|').split('|')] for l in lines[2:]}
    for i, tag in enumerate(tags):
        s = sum(raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] for p in X.DEV)
        f, c, pu = X.f_measure(s)
        mark = '**' if tag in (P2[1], C2[2]) else ''
        cells = rows['%s%s%s%s' % (mark, tag, mark, ' (round 1 selection)' if tag in SEL1 else '')]
        assert cells[0] == ('photutils' if tag in P2 + [SEL1[0]] else 'crowdsource') and cells[2] == H2.describe2(tag)
        assert cells[3:] == ['%.3f' % c, '%.3f' % pu, '%s%.3f%s' % (mark, f, mark), '%d' % s[2], '%d' % (s[6] - s[5]), '%.0f' % (20.0 * (i + 1) + 1)]
    assert rows['photutils (default)'][:3] == ['old family', 'photutils (default)', 'photutils'] and rows['photutils (default)'][3] == '%.3f' % (60 / 200)
    assert rows['crowdsource (tuned)'][2] == 'crowdsource/snr5/ff0/qf0/spinf' and 'photutils (tuned)' not in rows and 'crowdsource (default)' not in rows
    assert H2.dev_table2(raw, None).count('old family') == 0


def test_describe2_names_the_settings_of_each_variant():
    assert H2.describe2('phu2_hp3_all_i6_b2') == 'photutils, high-pass 3 px, mode all, 6 iterations, grouper 4.4 px, relaxed finder cuts, xy_bounds 2 px'
    assert H2.describe2('phu2_hp5_new_i3_b1') == 'photutils, high-pass 5 px, mode new, 3 iterations, xy_bounds 1 px'
    assert H2.describe2('cs2_hp2_t3_pm') == 'crowdsource, image and stamp high-passed 2 px, the ISS weights, threshold 3, passes m (miniter 8, maxiter 20, fewstars 20)'
    assert H2.describe2('cs2_hp3_t5_pd') == 'crowdsource, image and stamp high-passed 3 px, the ISS weights, threshold 5, passes d (miniter 4, maxiter 10, fewstars 60)'
    assert H2.describe2('cs_hp3_t5') == E1.describe('cs_hp3_t5') and H2.describe2('phu_hp3_all_i6') == E1.describe('phu_hp3_all_i6')


def test_from_raw_regenerates_the_report_byte_for_byte(monkeypatch, tmp_path):
    raw, raw1, fz, *_ = _inputs2(monkeypatch)
    outdir = tmp_path / 'out'
    outdir.mkdir()
    np.savez_compressed(outdir / 'baselines_hp2_raw.npz', **raw)
    np.savez_compressed(tmp_path / 'raw1.npz', **raw1)
    paths = {}
    for k, name in (('m2', 'm2_raw.npz'), ('psf', 'psf_raw.npz'), ('repeat', 'repeat_raw.npz')):
        np.savez_compressed(tmp_path / name, **fz[k])
        paths[k] = str(tmp_path / name)
    env = dict(os.environ, HP2_OUT=str(outdir), HP2_RAW1=str(tmp_path / 'raw1.npz'), HP_M2_RAW=paths['m2'], HP_PSF_RAW=paths['psf'], HP_REPEAT_RAW=paths['repeat'])
    cmd = [sys.executable, os.path.join(HERE, '..', 'experiments', 'exp_astro_hp2.py'), '--from-raw']
    first = None
    for _ in range(2):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        data = (outdir / 'baselines_hp2.md').read_bytes()
        assert first is None or data == first
        first = data
        (outdir / 'baselines_hp2.md').unlink()
    assert first == H2.report2(X.load_raw(str(outdir / 'baselines_hp2_raw.npz')), E1.load_frozen(paths), X.load_raw(str(tmp_path / 'raw1.npz'))).encode()
    line = [l for l in first.decode().split('\n') if l.startswith('Frozen results read:')][0]
    assert line == 'Frozen results read: %s.' % '; '.join('`%s` %s...' % (n, hashlib.sha256(open(paths[k], 'rb').read()).hexdigest()[:16]) for k, n in (('m2', 'm2_raw.npz'), ('psf', 'psf_raw.npz'), ('repeat', 'repeat_raw.npz')))


def test_the_command_line_runs_the_stage_asked_for(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(H2, 'part_dev2', lambda: calls.append(('dev',)) or {'photutils': 'a', 'crowdsource': 'b'})
    monkeypatch.setattr(H2, 'part_test2', lambda sets, flag=None: calls.append(('test', list(sets), flag)))
    for argv in (['--dev'], ['--test'], ['--test', '--sets', 'm92', 'm2r1', '--flag', 'f.flag']):
        monkeypatch.setattr(sys, 'argv', ['exp_astro_hp2.py'] + argv)
        H2.main()
    assert calls == [('dev',), ('test', ['m2', 'm53', 'm92'], None), ('test', ['m92', 'm2r1'], 'f.flag')]
    monkeypatch.setattr(sys, 'argv', ['exp_astro_hp2.py'])
    with pytest.raises(SystemExit):
        H2.main()
    np.savez_compressed(tmp_path / 'raw.npz', a=np.zeros(1))
    monkeypatch.setattr(H2, 'MD2', str(tmp_path / 'r.md'))
    monkeypatch.setattr(H2, 'RAW2', str(tmp_path / 'raw.npz'))
    monkeypatch.setattr(H2, 'report2', lambda raw: 'the report of %s\n' % sorted(raw))
    monkeypatch.setattr(sys, 'argv', ['exp_astro_hp2.py', '--from-raw'])
    H2.main()
    assert (tmp_path / 'r.md').read_text() == "the report of ['a']\n"
