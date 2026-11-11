"""Tests of experiments/exp_smlm_controls.py (post hoc controls of the SMLM comparison): the forced order of entry and the path statistics on planted frames, the claim that the truncated stop
changes no detection, the offset-corrected coordinate ascent against `smlm_common.tune_postfilters` on planted candidates (a world where the raw and the corrected objective choose
differently), the search and its resumption with a stub Picasso session, the Picasso part end to end with a fake worker, the greedy part end to end on planted frames, the report and its byte identity
on a planted raw, and a scorer re-derived independently (exhaustive pairing) checked on planted detections and, where the raw is present, on stored rows. Nothing needs the SMLM data or Picasso."""
import itertools
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import smlm_common as C                  # noqa: E402
import exp_smlm_controls as X            # noqa: E402
import exp_smlm_iss as ISS               # noqa: E402
from src.core import imaging as IM       # noqa: E402

SET = 'MT0.N1.LD'
SET13 = 'Tubulins_I'


# ----------------------------------------------------------------------------------------------------------------------------------
# planted worlds
# ----------------------------------------------------------------------------------------------------------------------------------
def pair_world(seed, sep=3.0, ph=2000.0):
    """24 x 24 noisy frame of two equal emitters `sep` widths apart (the regime where the order of entry matters: `order_pairs.md`), the PSF, the lattice and the weights."""
    sig = 1.4
    P = IM.PixelPSF(24, 24, sig)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    rng = np.random.default_rng(seed)
    y = P.synth(np.array([[11.0, 12.0], [11.0 + sep * sig, 12.3]]), np.array([ph, ph])) + rng.normal(0.0, 6.0, (24, 24))
    return P, L, y, np.full((24, 24), 1.0 / 6.0)


def isolated_frames(n, rng, sigma=1.4, bg=100.0, photons=(4000, 8000)):
    """n frames of 64 x 64 photon counts with two isolated emitters (more than 14 px apart) each, Poisson noise on a flat background."""
    out = np.empty((n, 64, 64))
    for i in range(n):
        f = np.full((64, 64), bg)
        pos = []
        while len(pos) < 2:
            c = rng.uniform(10, 54, 2)
            if all(np.hypot(*(c - q)) > 14 for q in pos):
                pos.append(c)
        for c in pos:
            cx, cy = int(round(c[0])), int(round(c[1]))
            f[cy - 5:cy + 6, cx - 5:cx + 6] += rng.uniform(*photons) * ISS._pix_gauss(c[0] - cx, c[1] - cy, sigma, 5)
        out[i] = rng.poisson(f).astype(float)
    return out


def ascent_world(seed=3, nfr=60):
    """Candidates of a Picasso-like detector: 60 % of the emitters are found with a high net gradient and a +60 nm position bias (spread 3 nm), 40 % with a low net gradient, no bias and 50 nm spread.
    Keeping everything is the better RAW choice (E 36 against 28), the net gradient filter is the better choice AFTER the offset is removed (E 60 against 47)."""
    rng = np.random.default_rng(seed)
    gts, cands = [], []
    for f in range(nfr):
        pos = rng.uniform(250, 6150, (10, 2))
        rows = []
        for i, (x, y) in enumerate(pos):
            gts.append((len(gts), f, x, y, 0.0, 3000.0))
            if i < 6:
                rows.append([x + 60 + rng.normal(0, 3), y + rng.normal(0, 3), 2000.0, 140.0, 150.0, 900.0])
            else:
                rows.append([x + rng.normal(0, 50), y + rng.normal(0, 50), 2000.0, 140.0, 150.0, 300.0])
        cands.append(np.array(rows))
    return cands, np.array(gts, dtype=C.GT_DTYPE), np.arange(nfr)


GRID = {'min_net_gradient': [None, 500.0]}


def brute_score(dets, gt, idx, tol=250.0):
    """E, J, RMSE by exhaustive search of the pairing of each frame (the maximum number of pairs within the tolerance, then the minimum total distance), pooled; written independently of
    `smlm_common.match` (no assignment solver)."""
    tp = fp = fn = 0
    sse = 0.0
    for d, f in zip(dets, idx):
        g = np.array([[r['x_nm'], r['y_nm']] for r in gt if r['frame'] == f]).reshape(-1, 2)
        d = np.asarray(d)[:, :2] if len(d) else np.zeros((0, 2))
        best = (-1, 0.0)

        def rec(i, used, cnt, tot):
            nonlocal best
            if i == len(d):
                if cnt > best[0] or (cnt == best[0] and tot < best[1] - 1e-12):
                    best = (cnt, tot)
                return
            rec(i + 1, used, cnt, tot)
            for j in range(len(g)):
                if j not in used:
                    dist = float(np.hypot(*(d[i] - g[j])))
                    if dist <= tol:
                        rec(i + 1, used | {j}, cnt + 1, tot + dist)
        rec(0, frozenset(), 0, 0.0)
        # the squared error of the optimal pairing: enumerate again for the pairing that attains it
        bestp = None

        def rec2(i, used, cnt, tot, pairs):
            nonlocal bestp
            if i == len(d):
                if cnt == best[0] and abs(tot - best[1]) < 1e-9 and bestp is None:
                    bestp = list(pairs)
                return
            rec2(i + 1, used, cnt, tot, pairs)
            for j in range(len(g)):
                if j not in used:
                    dist = float(np.hypot(*(d[i] - g[j])))
                    if dist <= tol:
                        rec2(i + 1, used | {j}, cnt + 1, tot + dist, pairs + [(i, j)])
        rec2(0, frozenset(), 0, 0.0, [])
        k = best[0]
        tp += k
        fp += len(d) - k
        fn += len(g) - k
        sse += sum(float(((d[i] - g[j]) ** 2).sum()) for i, j in bestp)
    jac = tp / (tp + fp + fn)
    rmse = np.sqrt(sse / tp)
    return dict(jaccard=jac, rmse_nm=rmse, efficiency=100.0 - np.sqrt((100.0 - 100.0 * jac) ** 2 + rmse ** 2))


# ----------------------------------------------------------------------------------------------------------------------------------
# the forced order of entry
# ----------------------------------------------------------------------------------------------------------------------------------
def stub_path(events, taus, barred=()):
    calls = []

    def f(*a, **kw):
        calls.append(kw)
        return dict(events=events, taus=np.array(taus), hist=[], snr=np.zeros(len(taus)), barred=list(barred))
    return f, calls


def test_the_forced_order_and_the_overrides_reach_the_path_and_the_function_is_restored(monkeypatch):
    ev = [(1.0, (0, 0), 'enter'), (1.0, (0, 0), 'exit'), (1.0, (0, 0), 'enter'), (1.0, (0, 0), 'exit')]
    f, calls = stub_path(ev, [1.0, 1.0])
    monkeypatch.setattr(IM, 'sliding_path', f)
    with X.OrderForcing('greedy', dict(snr_min=4.0, snr_extra=0)) as of:
        assert IM.sliding_path is not f
        IM.sliding_path(1, 2, max_entries=7, order='dual', snr_min=3.0, snr_extra=3, tol=1e-8)
    assert IM.sliding_path is f
    assert calls[0]['order'] == 'greedy' and calls[0]['snr_min'] == 4.0 and calls[0]['snr_extra'] == 0 and calls[0]['max_entries'] == 7 and calls[0]['tol'] == 1e-8
    assert of.stats.tolist() == [1, 2, 2, 1, 1, 1, 0] and of.barred == 0
    with pytest.raises(ValueError):
        with X.OrderForcing('dual'):
            raise ValueError('inside')
    assert IM.sliding_path is f
    with X.OrderForcing('dual') as of2:                                  # no override: only the order is touched
        IM.sliding_path(1, 2, snr_min=3.0)
    assert calls[-1]['order'] == 'dual' and calls[-1]['snr_min'] == 3.0 and 'snr_extra' not in calls[-1]
    assert of2.stats[0] == 1 and of2.stats[6] == 0
    with X.OrderForcing('dual') as of3:                                  # the default cap of the 2016 pipeline is 40 entries; a path of 2 entries is not at the cap with max_entries 3, is with 2
        IM.sliding_path(max_entries=3)
        IM.sliding_path(max_entries=2)
    assert of3.stats[0] == 2 and of3.stats[6] == 1 and of3.stats[1] == 4
    f40, _ = stub_path([(float(i), (i, i), 'enter') for i in range(40)], [float(i) for i in range(40)])
    monkeypatch.setattr(IM, 'sliding_path', f40)
    with X.OrderForcing('dual') as of4:                                  # no max_entries given: the 2016 pipeline's default cap of 40 entries
        IM.sliding_path()
    assert of4.stats.tolist() == [1, 40, 0, 0, 0, 0, 1]


def test_the_guarded_variant_calls_the_guarded_path_and_counts_its_barred_points(monkeypatch):
    f, calls = stub_path([(1.0, (0, 0), 'enter')], [1.0])
    g, gcalls = stub_path([(1.0, (0, 0), 'enter')], [1.0], barred=[(1.0, (0, 0)), (2.0, (1, 1))])
    monkeypatch.setattr(IM, 'sliding_path', f)
    monkeypatch.setattr(X.IG, 'sliding_path_guarded', g)
    with X.OrderForcing('greedy', guarded=True) as of:
        IM.sliding_path(max_entries=5)
    assert not calls and gcalls[0]['order'] == 'greedy' and gcalls[0]['cycle_guard'] is True and of.barred == 2 and of.stats[0] == 1
    with X.OrderForcing('greedy', guarded=False) as of:
        IM.sliding_path(max_entries=5)
    assert len(calls) == 1 and of.barred == 0


def test_path_statistics_count_exits_reentries_and_the_events_the_guard_bars_on():
    ev = [(1.0, (0, 0), 'enter'), (1.0, (0, 0), 'exit'), (1.0, (0, 0), 'enter'), (1.0, (0, 0), 'exit'), (2.0, (3, 3), 'enter'), (2.5, (3, 3), 'exit'), (3.0, (3, 3), 'enter'), (4.0, (3, 3), 'exit'),
          (5.0, (1, 1), 'enter')]
    assert X.path_stats(ev, 5, 40).tolist() == [1, 5, 4, 2, 1, 1, 0]
    assert X.path_stats(ev, 5, 5).tolist() == [1, 5, 4, 2, 1, 1, 1]                                  # at the cap
    assert X.path_stats([], 0, 40).tolist() == [1, 0, 0, 0, 0, 0, 0]
    ev2 = [(1e6, (0, 0), 'exit'), (1e6 * (1 + 5e-13), (0, 0), 'exit')]                                  # equal within the relative tolerance of the guard
    assert X.path_stats(ev2, 1, 40)[4] == 1
    ev3 = [(1e6, (0, 0), 'exit'), (1e6 * (1 + 5e-11), (0, 0), 'exit')]
    assert X.path_stats(ev3, 1, 40)[4] == 0
    ev4 = [(0.0, (0, 0), 'exit'), (1e-13, (0, 0), 'exit')]                                              # absolute scale near 0: the tolerance is relative to max(1, |tau|)
    assert X.path_stats(ev4, 1, 40)[4] == 1
    assert len(X.STAT_NAMES) == len(X.path_stats([], 0, 1)) == 7


def test_on_a_planted_close_pair_the_two_orders_enter_different_second_candidates_and_greedy_the_more_significant():
    P, L, y, sw = pair_world(0)
    out = {}
    for order in ('dual', 'greedy'):
        with X.OrderForcing(order) as of:
            out[order] = IM.sliding_path(P, y, sw, L, max_entries=4, r_ex=0.0, tol=1e-8, passes=1, rtrunc=4.0)
        assert of.stats[0] == 1 and of.stats[1] == 4
    pos = {o: np.array([L.pos(*e[1]) for e in out[o]['events'] if e[2] == 'enter']) for o in out}
    assert np.array_equal(pos['dual'][0], pos['greedy'][0])                                              # the first entrant is the same (the correlations are constant until it is fitted)
    assert not np.array_equal(pos['dual'][1], pos['greedy'][1]) and not np.array_equal(pos['dual'][2], pos['greedy'][2])
    e1, e2 = np.array([11.0, 12.0]), np.array([11.0 + 3.0 * 1.4, 12.3])
    assert np.hypot(*(pos['dual'][1] - pos['greedy'][1])) > 1.3 and np.hypot(*(pos['greedy'][1] - e2)) < 1.0 and np.hypot(*(pos['dual'][1] - e2)) < 1.0      # both second entrants sit on the second emitter, 1.5 px apart
    assert np.hypot(*(pos['dual'][2] - pos['dual'][0])) < 1.0 and out['dual']['snr'][2] < 1.0                                   # the dual spends its third entry on a copy of the first emitter, at no significance
    assert min(np.hypot(*(pos['greedy'][2] - e1)), np.hypot(*(pos['greedy'][2] - e2))) > 4.0
    assert out['greedy']['snr'][1] > 1.3 * out['dual']['snr'][1]
    assert out['greedy']['taus'].tolist() == [1.0, 2.0, 3.0, 4.0] and not np.allclose(out['dual']['taus'], [1, 2, 3, 4])
    Q, LQ, yq, swq = pair_world(0, sep=8.0)                                                              # far apart (11 px): the same two entrants in both orders
    ev = {}
    for order in ('dual', 'greedy'):
        with X.OrderForcing(order):
            ev[order] = IM.sliding_path(Q, yq, swq, LQ, max_entries=2, r_ex=0.0, tol=1e-8, passes=1, rtrunc=4.0)['events']
    assert sorted(e[1] for e in ev['dual']) == sorted(e[1] for e in ev['greedy'])


def test_the_guarded_greedy_path_is_the_unguarded_greedy_path_with_no_barred_point():
    P, L, y, sw = pair_world(1)
    res = {}
    for guarded in (False, True):
        with X.OrderForcing('greedy', guarded=guarded) as of:
            res[guarded] = IM.sliding_path(P, y, sw, L, max_entries=4, r_ex=0.0, tol=1e-8, passes=1, rtrunc=4.0)
        assert of.barred == 0
    assert res[False]['events'] == res[True]['events'] and np.array_equal(res[False]['taus'], res[True]['taus'])
    assert np.array_equal(res[False]['hist'][-1]['w'], res[True]['hist'][-1]['w'])


# ----------------------------------------------------------------------------------------------------------------------------------
# the 2016 pipeline through the harness
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_truncated_stop_changes_no_detection_and_greedy_equals_dual_on_isolated_emitters():
    rng = np.random.default_rng(11)
    frames = isolated_frames(3, rng)
    cal = dict(bg=np.full((64, 64), 100.0), a=1.0)
    with X.OrderForcing('dual') as of_full:
        full, _, _, _ = ISS.run_frames(frames, 'hp3', cal, dict(rules=ISS.rules_dev()), workers=1)
    arr_full = {}
    ISS._store(arr_full, 'hp3', ISS.rules_dev(), full)
    tr, st_tr, secs, barred = X.run_cfg16('x', frames, 'hp3', cal, 'dual', X.TRUNC)
    assert barred == 0 and len(secs) == 3 and X.TRUNC == dict(snr_min=4.0, snr_extra=0)
    assert set(tr) == set(arr_full) and all(np.array_equal(tr[k], arr_full[k]) for k in tr)               # every rule, off and flat, bit for bit
    assert st_tr[0] == 3 and st_tr[1] < of_full.stats[1] and st_tr[6] == 0
    gr, st_gr, _, _ = X.run_cfg16('x', frames, 'hp3', cal, 'greedy', X.TRUNC)
    for k in (4.0, 5.0, 10.0):
        a, b = ISS.unpack(tr['hp3|snr%g|off' % k], tr['hp3|snr%g|flat' % k]), ISS.unpack(gr['hp3|snr%g|off' % k], gr['hp3|snr%g|flat' % k])
        assert sum(len(x) for x in a) >= 4 and [len(x) for x in a] == [len(x) for x in b]
        assert all(np.allclose(x, y, atol=1e-6, rtol=0) for x, y in zip(a, b))
    assert st_gr[0] == 3 and st_gr[1] == st_tr[1]
    # a kappa above the significance of every entrant removes them all: the rule list covers the seven levels
    assert len(ISS.rules_dev()) == 7 and all('hp3|snr%g|off' % k in tr for k in ISS.KAPPAS)


def test_validate16_compares_the_dual_path_with_the_stored_part(monkeypatch):
    rng = np.random.default_rng(12)
    frames = isolated_frames(2, rng)
    cal = dict(bg=np.full((64, 64), 100.0), a=1.0)
    ops = {'E': dict(cfg='hp3')}
    dets, _, _, _ = ISS.run_frames(frames, 'hp3', cal, dict(rules=ISS.rules_dev()), workers=1)
    iss = {}
    for r, rule in enumerate(ISS.rules_dev()):
        off, flat = ISS.pack(dets[r])
        iss['T|S|hp3|snr%g|off' % rule[1]], iss['T|S|hp3|snr%g|flat' % rule[1]] = off, flat
    logs = []
    v = X.validate16('S', frames, ops, cal, iss, logs.append)
    assert v['dual_original_stop']['identical'] and v['dual_truncated_stop']['identical'] and v['dual_truncated_stop']['rules_identical'] == 7 and v['dual_original_stop']['frames'] == 2
    assert v['dual_original_stop']['cfg'] == 'hp3' and len(logs) == 2 and '7 of 7' in logs[0]
    iss['T|S|hp3|snr5|flat'] = iss['T|S|hp3|snr5|flat'] + 1e-9                                           # one rule differs in the last bits
    v2 = X.validate16('S', frames, ops, cal, iss, logs.append)
    assert not v2['dual_original_stop']['identical'] and v2['dual_original_stop']['rules_identical'] == 6 and not v2['dual_truncated_stop']['identical']
    g = X.validate_guard('S', frames, ops, cal, 2, logs.append)
    assert g == dict(set='S', cfg='hp3', frames=2, identical=True, barred=0)
    # which order and which stop each check runs, and what a guarded path that differs would be recorded as
    seen = []
    real = X.run_cfg16

    def spy(name, fr, cfg, cal_, order, override=None, guarded=False):
        seen.append((order, dict(override or {}), guarded, len(fr)))
        arr, st, secs, barred = real(name, fr, cfg, cal_, order, override, guarded)
        if guarded:
            arr = {k: v + (1.0 if k.endswith('flat') and len(v) else 0.0) for k, v in arr.items()}
            barred = 3
        return arr, st, secs, barred
    monkeypatch.setattr(X, 'run_cfg16', spy)
    iss['T|S|hp3|snr5|flat'] = iss['T|S|hp3|snr5|flat'] - 1e-9
    X.validate16('S', frames, ops, cal, iss, logs.append)
    assert seen == [('dual', X.TRUNC, False, 2)]                                                         # the original-stop check goes through ISS.run_frames, the truncated one through run_cfg16
    seen.clear()
    g2 = X.validate_guard('S', frames, ops, cal, 1, logs.append)
    assert seen == [('greedy', X.TRUNC, False, 1), ('greedy', X.TRUNC, True, 1)] and g2 == dict(set='S', cfg='hp3', frames=1, identical=False, barred=3)


# ----------------------------------------------------------------------------------------------------------------------------------
# the offset-corrected objective and the ascent
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_offset_corrected_metrics_remove_the_constant_offset_and_report_it():
    cands, gt, idx = ascent_world()
    rng = np.random.default_rng(5)
    dets = []
    for f in idx:
        g = gt[gt['frame'] == f]
        dets.append(np.c_[g['x_nm'] + 40.0 + rng.normal(0, 4, len(g)), g['y_nm'] - 25.0 + rng.normal(0, 4, len(g)), np.full(len(g), 1000.0), np.full(len(g), 140.0)])
    raw = X.plain_metrics(dets, gt, idx)
    ref0 = C.score(dets, gt, idx)
    assert all(raw[k] == ref0[k] for k in C._TUNE_KEYS)                                                  # no border exclusion in the metrics of the search
    cor = X.offset_corrected_metrics(dets, gt, idx)
    assert cor['off_x'] == ref0['bias_x_nm'] and cor['off_y'] == ref0['bias_y_nm']                      # the offset is that of the pairing of all emitters
    assert raw['rmse_nm'] > 45 and cor['efficiency'] > 93.5 and abs(cor['off_x'] - 40) < 1.0 and abs(cor['off_y'] + 25) < 1.0
    assert 4.5 < cor['rmse_nm'] < 7.0 and cor['jaccard'] == 1.0 and set(X.plain_metrics(dets, gt, idx)) == set(C._TUNE_KEYS)
    empty = [np.zeros((0, 6)) for _ in idx]
    e = X.offset_corrected_metrics(empty, gt, idx)
    assert e['off_x'] == 0.0 and e['off_y'] == 0.0 and e['tp'] == 0 and np.isnan(e['efficiency'])      # nothing matched: the offset is 0, not NaN


def test_the_ascent_with_the_plain_evaluator_is_tune_postfilters_and_the_corrected_objective_chooses_differently():
    cands, gt, idx = ascent_world()
    grids = {'min_net_gradient': [None, 200.0, 500.0, 800.0], 'min_photons': [None, 1500.0]}
    for obj in ('efficiency', 'jaccard'):
        ref = C.tune_postfilters(cands, gt, idx, grids, objective=obj, max_rounds=3)
        mine = X.ascent_offset(cands, gt, idx, grids, obj, {}, {}, 3, evaluator=X.plain_metrics)
        assert mine['filters'] == ref['filters'] and mine['rounds'] == ref['rounds'] and mine['converged'] == ref['converged'] and mine['n_eval'] == ref['n_eval']
        assert abs(mine['value'] - ref['value']) < 1e-12 and mine['first']['filters'] == ref['first']['filters']
    plain = X.ascent_offset(cands, gt, idx, GRID, 'efficiency', {}, {}, 3, evaluator=X.plain_metrics)
    corr = X.ascent_offset(cands, gt, idx, GRID, 'efficiency', {}, {}, 3)
    assert plain['filters'] == {'min_net_gradient': None} and corr['filters'] == {'min_net_gradient': 500.0}
    assert 34 < plain['metrics']['efficiency'] < 39 and 57 < corr['metrics']['efficiency'] < 62 and abs(corr['metrics']['off_x'] - 60) < 3
    corr_j = X.ascent_offset(cands, gt, idx, GRID, 'jaccard', {}, {}, 3)                           # J: keeping everything is the better J whatever the offset
    assert corr_j['filters'] == {'min_net_gradient': None} and corr_j['metrics']['jaccard'] == 1.0
    # the cache is shared between objectives, the detection digest skips filters that keep the same detections
    cache, dcache = {}, {}
    X.ascent_offset(cands, gt, idx, {'min_net_gradient': [None, 100.0, 200.0, 500.0], 'min_photons': [None, 100.0]}, 'efficiency', cache, dcache, 3)
    n1 = len(cache)
    X.ascent_offset(cands, gt, idx, {'min_net_gradient': [None, 100.0, 200.0, 500.0], 'min_photons': [None, 100.0]}, 'jaccard', cache, dcache, 3)
    n_filter, calls = [], []
    real_filter = C.filter_detections

    def counting_filter(dets, f):
        n_filter.append(1)
        return real_filter(dets, f)
    C.filter_detections = counting_filter
    try:
        c2, d2 = {}, {}
        X.ascent_offset(cands, gt, idx, {'min_net_gradient': [None, 100.0, 200.0, 500.0], 'min_photons': [None, 100.0]}, 'efficiency', c2, d2, 3,
                        evaluator=lambda dets, g, i: (calls.append(1), X.plain_metrics(dets, g, i))[1])
    finally:
        C.filter_detections = real_filter
    assert len(n_filter) == len(c2) > 2 and len(calls) == len(d2) == 2                                  # one filtering per filter combination, one scoring per distinct set of detections
    assert n1 > 2 and len(cache) <= 8 and len(dcache) == 2                                                  # 8 combinations, 2 distinct sets of detections (everything / the high-gradient ones)


def test_the_ascent_takes_any_strict_improvement_and_records_its_first_line_search():
    cands, gt, idx = ascent_world()
    small = lambda dets, g, i: dict(efficiency=10.0 if sum(len(d) for d in dets) == 600 else 10.5, jaccard=0.0)
    r = X.ascent_offset(cands, gt, idx, GRID, 'efficiency', {}, {}, 3, evaluator=small)
    assert r['filters'] == {'min_net_gradient': 500.0} and r['value'] == 10.5                             # half a point is taken
    tiny = lambda dets, g, i: dict(efficiency=10.0 if sum(len(d) for d in dets) == 600 else 10.0 + 1e-13, jaccard=0.0)
    assert X.ascent_offset(cands, gt, idx, GRID, 'efficiency', {}, {}, 3, evaluator=tiny)['filters'] == {'min_net_gradient': None}       # below 1e-12 is a tie
    # 'first' is the state after the first filter's line search of the FIRST round, even when a later round moves that filter
    r2 = X.ascent_offset(cands, gt, idx, {'min_photons': [None, 1500.0], 'min_net_gradient': [None, 500.0]}, 'efficiency', {}, {}, 3)
    assert r2['first']['filters'] == {'min_photons': None, 'min_net_gradient': None} and r2['filters']['min_net_gradient'] == 500.0 and r2['rounds'] == 2


def test_the_ascent_stops_when_a_round_changes_nothing_and_never_moves_on_a_tie():
    cands, gt, idx = ascent_world()
    r = X.ascent_offset(cands, gt, idx, {'min_net_gradient': [None, 250.0, 500.0, 800.0]}, 'efficiency', {}, {}, 5)
    assert r['filters'] == {'min_net_gradient': 500.0}                                                      # 250 keeps everything (a tie with None, not taken); 500 and 800 keep the same detections: the first is kept
    assert r['rounds'] == 2 and r['converged'] and r['first']['filters'] == {'min_net_gradient': 500.0} and 'first' in r
    r1 = X.ascent_offset(cands, gt, idx, {'min_net_gradient': [None, 250.0, 500.0, 800.0]}, 'efficiency', {}, {}, 1)
    assert r1['rounds'] == 1 and not r1['converged']
    nan_eval = lambda dets, g, i: dict(efficiency=float('nan'), jaccard=1.0)
    rn = X.ascent_offset(cands, gt, idx, GRID, 'efficiency', {}, {}, 3, evaluator=nan_eval)
    assert rn['filters'] == {'min_net_gradient': None} and rn['value'] == -np.inf


def test_digests_of_detections_distinguish_values_counts_and_frames():
    a = [np.array([[1.0, 2.0, 3.0]]), np.zeros((0, 3)), np.array([[4.0, 5.0, 6.0]])]
    b = [x.copy() for x in a]
    assert X.digest_dets(a) == X.digest_dets(b)
    assert X.digest_dets(a) != X.digest_dets([a[0], a[2], a[1]])
    c = [x.copy() for x in a]
    c[2][0, 2] = 6.0 + 1e-9
    assert X.digest_dets(a) != X.digest_dets(c)
    assert X.digest_dets([a[0], np.zeros((0, 3))]) != X.digest_dets([np.zeros((0, 3)), a[0]])


def test_the_structural_choice_takes_the_best_corrected_value_and_the_earlier_setting_on_ties():
    def entry(e, j):
        return dict(E=dict(metrics=dict(efficiency=e, jaccard=0.0)), J=dict(metrics=dict(efficiency=0.0, jaccard=j)))
    res = {'a': entry(50.0, 0.5), 'b': entry(60.0, 0.5), 'c': entry(60.0, 0.7), 'd': entry(55.0, 0.7)}
    assert X.choose_offset(res) == ('b', 'c')
    assert X.choose_offset({'only': entry(1.0, 0.1)}) == ('only', 'only')


def test_same_detections_distinguishes_rounding_from_bits_and_counts():
    a = [np.array([[1.0, 2.0, 3.0, 4.0, 9.0, 9.0]]), np.zeros((0, 6))]
    b = [np.array([[1.0, 2.0, 3.0, 4.0]]), np.zeros((0, 4))]
    assert X.same_detections(a, b) == (True, True)                                                       # the first four columns only
    c = [np.array([[1.0, 2.0 + 1e-7, 3.0, 4.0]]), np.zeros((0, 4))]
    assert X.same_detections(a, c) == (True, False)
    d = [np.array([[1.0, 2.5, 3.0, 4.0]]), np.zeros((0, 4))]
    assert X.same_detections(a, d) == (False, False)
    assert X.same_detections(a, [b[0]]) == (False, False) and X.same_detections(a, [np.zeros((0, 4)), np.zeros((0, 4))]) == (False, False)
    f32 = [np.array([[1.1, 2.0, 3.0, 4.0]]), np.zeros((0, 4))]
    ref32 = [f32[0].astype(np.float32).astype(np.float64), np.zeros((0, 4))]
    assert X.same_detections(f32, ref32, float32=True) == (True, True) and X.same_detections(f32, ref32)[1] is False


# ----------------------------------------------------------------------------------------------------------------------------------
# the search with a stub session, and the Picasso part with a fake worker
# ----------------------------------------------------------------------------------------------------------------------------------
class StubSession:
    """A Picasso session that returns the planted candidates filtered by the post-filters of the request (what a real run does with the net gradient)."""

    def __init__(self, cands, calls=None):
        self.cands, self.calls = cands, calls if calls is not None else []

    def run(self, params=None, info=None):
        self.calls.append(dict(params))
        return C.filter_detections(self.cands, params)


def test_the_search_runs_every_setting_once_chooses_nothing_itself_and_resumes_from_its_cache(tmp_path):
    cands, gt, idx = ascent_world()
    structs = X.PIC.struct_list()[:3]
    calls = []
    cache_file = str(tmp_path / 'sub' / 'x.json')
    out = X.search_offset(StubSession(cands, calls), gt, idx, {}, 100.0, GRID, structs, cache_file, 'sig1', None, lambda s: None)
    assert list(out) == [X.PIC.skey(s) for s in structs] and len(calls) == 3 and all(c['min_net_gradient'] == 100.0 for c in calls)
    assert all(set(c) >= set(structs[0]) for c in calls)
    r = out[X.PIC.skey(structs[0])]
    assert r['struct'] == dict(structs[0]) and r['n_cand'] == sum(len(d) for d in cands) and r['n_eval'] == 2 and r['E']['filters'] == {'min_net_gradient': 500.0} and 'first' not in r['E']
    assert r['J']['filters'] == {'min_net_gradient': None}
    assert r['E']['rounds'] == 2 and r['E']['converged'] and r['J']['rounds'] == 1 and X.PIC.MAX_ROUNDS == 3
    calls2 = []
    out2 = X.search_offset(StubSession(cands, calls2), gt, idx, {}, 100.0, GRID, structs, cache_file, 'sig1', None, lambda s: None)
    assert calls2 == [] and out2 == json.loads(json.dumps(out, default=C._json_default))                   # everything read from the cache
    out3 = X.search_offset(StubSession(cands, calls2), gt, idx, {}, 100.0, GRID, structs[:2], cache_file, 'other', None, lambda s: None)       # another signature: recomputed
    assert len(calls2) == 2 and len(out3) == 2
    calls3 = []
    out4 = X.search_offset(StubSession(cands, calls3), gt, idx, {'baseline': 100.0}, 100.0, GRID, structs, None, '', 2, lambda s: None)          # max_structs, the camera is passed on
    assert len(out4) == 2 and calls3[0]['baseline'] == 100.0
    best = X.choose_offset(out)
    assert best[0] == X.PIC.skey(structs[0]) and best[1] == X.PIC.skey(structs[0])


def test_the_picasso_part_selects_with_the_corrected_objective_and_stores_the_corrected_test_detections(tmp_path, monkeypatch):
    cands, gt, idx = ascent_world(seed=3, nfr=60)
    dev_c, dev_gt, dev_idx = cands[:30], gt[gt['frame'] < 30], np.arange(30)
    test_c = cands[30:]
    gt_all = gt.copy()

    class Fake(StubSession):
        def __init__(self, frames, py=None):
            super().__init__(frames)                                                                     # the "frames" are the candidate lists

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def run(self, params=None, info=None):
            Fake.seen.append(dict(params))
            return C.filter_detections(self.cands, dict(params))
    Fake.seen = []
    stored = {'R2E': dict(min_net_gradient=None, box=7, fit_box=7), 'R2J': dict(min_net_gradient=500.0, box=7, fit_box=7)}
    stored_dev = {r: Fake(dev_c).run(stored[r]) for r in stored}
    Fake.seen = []
    P = dict(dev=dev_c, test=test_c, d=(0, 30), t=(30, 60), gt=gt_all, grids=GRID, ng_low=100.0, cam={'baseline': 100.0}, stored=stored, stored_dev=stored_dev)
    monkeypatch.setattr(X, 'pic_setup', lambda name: P)
    monkeypatch.setattr(X.C, 'PicassoSession', Fake)
    monkeypatch.setattr(X.PIC, 'struct_list', lambda: [dict(box=5, fit_box=5, gauss_sigma=None, pad=0, fit='gausslq'), dict(box=7, fit_box=7, gauss_sigma=None, pad=0, fit='gausslq')])
    raw_path = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(X, 'RAW', raw_path)
    logs = []
    X.part_ps(SET, None, str(tmp_path / 'cache'), None, logs.append)
    raw = X.load_raw(raw_path)
    s = X.jload(raw['ps|%s|summary' % SET])
    assert s['structs'] == 2 and s['info']['R2Eo']['key'].startswith('box=5') and s['info']['R2Jo']['key'].startswith('box=5')                  # ties: the earlier setting
    assert X.jload(raw['ps|%s|validation' % SET]) == dict(R2E=dict(identical_to_rounding=True, bit_identical=True, detections=int(sum(len(d) for d in stored_dev['R2E']))),
                                                          R2J=dict(identical_to_rounding=True, bit_identical=True, detections=int(sum(len(d) for d in stored_dev['R2J']))))
    cfgE, cfgJ = X.jload(raw['ps|%s|R2Eo|config' % SET]), X.jload(raw['ps|%s|R2Jo|config' % SET])
    assert cfgE['min_net_gradient'] == 500.0 and cfgJ['min_net_gradient'] is None
    off = raw['ps|%s|R2Eo|offset' % SET]
    assert abs(off[0] - 60.0) < 3.0 and abs(off[1]) < 3.0
    # the stored detections are the chosen configuration's, with the development offset subtracted, on the dev and the test frames
    dev_d, test_d = X.unpack(raw, 'ps|%s|R2Eo|dev' % SET), X.unpack(raw, 'ps|%s|R2Eo|test' % SET)
    direct_dev, direct_test = C.filter_detections(dev_c, dict(min_net_gradient=500.0)), C.filter_detections(test_c, dict(min_net_gradient=500.0))
    assert len(dev_d) == 30 and len(test_d) == 30
    assert all(np.allclose(a[:, 0], b[:, 0] - off[0]) and np.allclose(a[:, 1], b[:, 1] - off[1]) and a.shape[1] == 4 for a, b in zip(test_d, direct_test) if len(a))
    assert all(np.allclose(a[:, 0], b[:, 0] - off[0]) for a, b in zip(dev_d, direct_dev) if len(a))
    sc = X.jload(raw['ps|%s|stored_corrected' % SET])
    assert set(sc) == {'R2E', 'R2J'} and abs(sc['R2E']['off_x'] - 36.0) < 4.0 and sc['R2E']['efficiency'] < 50
    assert any('validation R2E' in m for m in logs) and any(m.startswith('PS ') for m in logs)
    assert Fake.seen and all(c.get('baseline') == 100.0 for c in Fake.seen)                              # the camera reaches every run, the stored ones included
    assert [c.get('min_net_gradient') for c in Fake.seen[:2]] == [None, 500.0]
    tj = X.unpack(raw, 'ps|%s|R2Jo|test' % SET)
    offj = raw['ps|%s|R2Jo|offset' % SET]
    assert all(np.allclose(a[:, 0], b[:, 0] - offj[0]) and a.shape == (len(b), 4) for a, b in zip(tj, test_c) if len(a)) and sum(len(a) for a in tj) == sum(len(b) for b in test_c)      # the J row is run with its own (unfiltered) configuration
    assert np.array_equal(raw['ps|%s|R2Eo|test|flat' % SET][:, 3], np.concatenate([d[:, 3] for d in direct_test if len(d)]))               # the fourth column is the width sx, not sy
    X.part_ps(SET, None, str(tmp_path / 'cache_ms'), 1, logs.append)
    assert X.jload(X.load_raw(raw_path)['ps|%s|summary' % SET])['structs'] == 1
    res = X.jload(raw['ps|%s|search' % SET])
    assert s['n_eval'] == sum(r['n_eval'] for r in res.values()) and s['n_cand'] == sum(r['n_cand'] for r in res.values()) and X.jload(raw['ps|%s|stored_config' % SET]) == stored
    # a stored configuration that the worker does not reproduce is recorded as such, and the filtered low-threshold run must equal the direct run
    P['stored_dev'] = {r: [d + 1.0 if len(d) else d for d in stored_dev[r]] for r in stored}
    X.part_ps(SET, None, str(tmp_path / 'cache2'), None, logs.append)
    v = X.jload(X.load_raw(raw_path)['ps|%s|validation' % SET])
    assert not v['R2E']['bit_identical'] and not v['R2E']['identical_to_rounding']
    # a 2013 set: the stored detections were float32, and are compared as such
    P['stored_dev'] = {r: [np.asarray(d, np.float32).astype(np.float64) for d in stored_dev[r]] for r in stored}
    X.part_ps(SET13, None, str(tmp_path / 'cache3'), None, logs.append)
    v13 = X.jload(X.load_raw(raw_path)['ps|%s|validation' % SET13])
    assert v13['R2E']['bit_identical'] and v13['R2E']['identical_to_rounding']
    X.part_ps(SET, None, str(tmp_path / 'cache4'), None, logs.append)                                       # the same float32 detections against a float64 set: equal to rounding only
    v16 = X.jload(X.load_raw(raw_path)['ps|%s|validation' % SET])
    assert v16['R2E']['identical_to_rounding'] and not v16['R2E']['bit_identical']

    class Drift(Fake):                                                                                       # a worker whose low-threshold run is not the filtered direct run
        def run(self, params=None, info=None):
            out = C.filter_detections(self.cands, dict(params))
            return [d + 5.0 if (len(d) and params.get('min_net_gradient') == 100.0) else d for d in out]
    monkeypatch.setattr(X.C, 'PicassoSession', Drift)
    with pytest.raises(RuntimeError, match='differs from the direct run'):
        X.part_ps(SET, None, str(tmp_path / 'cache5'), None, logs.append)


# ----------------------------------------------------------------------------------------------------------------------------------
# the greedy part on planted frames
# ----------------------------------------------------------------------------------------------------------------------------------
def test_the_greedy_part_stores_the_detections_of_the_needed_configurations_and_validates_the_dual(tmp_path, monkeypatch):
    rng = np.random.default_rng(21)
    frames_ph = isolated_frames(3, rng)
    frames_dn = frames_ph * C.DN_PER_PHOTON + C.BASELINE_DN
    cal = dict(bg=np.full((64, 64), 100.0), a=1.0)
    raw_path = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(X, 'RAW', raw_path)
    ops = dict(E=dict(cfg='hp3', kappa=5.0, smax=2.5, off=[1.0, -1.0], m=[0.0] * 7), J=dict(cfg='hp3h', kappa=4.5, smax=3.2, off=[0.0, 0.0], m=[0.0] * 7))
    X.save_part(raw_path, 'ref|' + SET, dict(ops=X.jdump(ops)))
    dets, _, _, _ = ISS.run_frames(frames_ph, 'hp3', cal, dict(rules=ISS.rules_dev()), workers=1)
    fake_iss = {'C|%s|bg' % SET: cal['bg'], 'C|%s|a' % SET: np.array([cal['a']])}
    for r, rule in enumerate(ISS.rules_dev()):
        off, flat = ISS.pack(dets[r])
        fake_iss['T|%s|hp3|snr%g|off' % (SET, rule[1])], fake_iss['T|%s|hp3|snr%g|flat' % (SET, rule[1])] = off, flat
    monkeypatch.setattr(X.ISS, 'load_raw', lambda path=None: fake_iss)
    monkeypatch.setattr(X.C, 'load_frames', lambda path, a, b: frames_dn[:b - a])
    monkeypatch.setattr(X.ISS, 'dev_test_frames', lambda name, which, n=None: (1000, 1000 + (n or 500)) if which == 'test' else (0, 1000))
    calls = []
    real = X.run_cfg16

    def spy(name, fr, cfg, cal_, order, override=None, guarded=False):
        calls.append((cfg, order, dict(override or {}), guarded))
        return real(name, fr, cfg, cal_, order, override, guarded)
    monkeypatch.setattr(X, 'run_cfg16', spy)
    logs = []
    X.part_g16(SET, 3, logs.append)
    assert [(c[0], c[1]) for c in calls[:6]] == [('hp3', 'greedy'), ('hp3h', 'greedy'), ('hp5h', 'greedy'), ('hp3', 'dual'), ('hp3h', 'dual'), ('hp5h', 'dual')]       # the three configurations, greedy then dual
    assert all(c[2] == X.TRUNC and not c[3] for c in calls[:6])
    raw = X.load_raw(raw_path)
    assert raw['g|%s|frames' % SET].tolist() == [1000, 1003]
    for cfg in ('hp3', 'hp3h', 'hp5h'):                                                                   # the E, J and default configurations
        assert 'g|%s|%s|snr5|off' % (SET, cfg) in raw and 'g|%s|%s|snr4.5|flat' % (SET, cfg) in raw and raw['g|%s|%s|stats' % (SET, cfg)][0] == 3 and raw['g|%s|%s|dual_stats' % (SET, cfg)][0] == 3
        assert raw['g|%s|%s|secs' % (SET, cfg)].shape == (3,) and 6 <= raw['g|%s|%s|stats' % (SET, cfg)][1] <= 12      # the truncated stop: two emitters and the first entrant below 4 in each of three paths
    assert SET == 'MT0.N1.LD' and X.jload(raw['val|%s|dual_original_stop' % SET])['identical'] and X.jload(raw['val|%s|dual_truncated_stop' % SET])['identical'] and 'val|%s|guard' % SET not in raw
    g = X.unpack(raw, 'g|%s|hp3|snr5' % SET)
    assert len(g) == 3 and all(len(x) == 2 for x in g)                                                    # two isolated emitters in every frame
    ref = ISS.unpack(*(fake_iss['T|%s|hp3|snr5|%s' % (SET, k)] for k in ('off', 'flat')))
    assert all(np.allclose(x, y, atol=1e-6, rtol=0) for x, y in zip(g, ref))
    X.save_part(raw_path, 'ref|MT0.N1.HD', dict(ops=X.jdump(ops)))                                          # the dense set runs the guard check instead (its camera is that of N1 LD)
    X.part_g16('MT0.N1.HD', 3, logs.append)
    raw_hd = X.load_raw(raw_path)
    assert X.jload(raw_hd['val|MT0.N1.HD|guard']) == dict(set='MT0.N1.HD', cfg='hp3', frames=3, identical=True, barred=0) and 'val|MT0.N1.HD|dual_original_stop' not in raw_hd
    other = str(tmp_path / 'raw2.npz')                                                                       # without validation: no val keys, same detections
    monkeypatch.setattr(X, 'RAW', other)
    X.save_part(other, 'ref|' + SET, dict(ops=X.jdump(ops)))
    X.part_g16(SET, 3, logs.append, validate=False)
    raw2 = X.load_raw(other)
    assert not [k for k in raw2 if k.startswith('val|')] and np.array_equal(raw2['g|%s|hp3|snr5|flat' % SET], raw['g|%s|hp3|snr5|flat' % SET])


def test_the_guard_check_runs_on_at_most_twenty_frames(tmp_path, monkeypatch):
    raw_path = str(tmp_path / 'raw.npz')
    monkeypatch.setattr(X, 'RAW', raw_path)
    ops = dict(E=dict(cfg='hp3', kappa=5.0, smax=2.5, off=[0.0, 0.0], m=[0.0] * 7), J=dict(cfg='hp3', kappa=5.0, smax=2.5, off=[0.0, 0.0], m=[0.0] * 7))
    X.save_part(raw_path, 'ref|MT0.N1.HD', dict(ops=X.jdump(ops)))
    fake_iss = {'C|MT0.N1.LD|bg': np.zeros((64, 64)), 'C|MT0.N1.LD|a': np.array([1.0])}
    monkeypatch.setattr(X.ISS, 'load_raw', lambda path=None: fake_iss)
    monkeypatch.setattr(X.C, 'load_frames', lambda path, a, b: np.zeros((b - a, 64, 64)))
    monkeypatch.setattr(X.ISS, 'dev_test_frames', lambda name, which, n=None: (100, 100 + n))
    monkeypatch.setattr(X, 'run_cfg16', lambda name, fr, cfg, cal, order, override=None, guarded=False: ({'hp3|snr5|off': np.zeros(len(fr) + 1, np.int64), 'hp3|snr5|flat': np.zeros((0, 4))}, np.zeros(7, np.int64), np.zeros(len(fr)), 0))
    seen = []
    monkeypatch.setattr(X, 'validate_guard', lambda name, fr, ops_, cal, n, log: (seen.append((len(fr), n)), dict(set=name, cfg='hp3', frames=int(n), identical=True, barred=0))[1])
    for n_frames, want in ((25, 20), (7, 7)):
        seen.clear()
        X.part_g16('MT0.N1.HD', n_frames, lambda s: None)
        assert seen == [(n_frames, want)]
    assert X.jload(X.load_raw(raw_path)['val|MT0.N1.HD|guard'])['frames'] == 7


def test_the_frames_that_differ_are_counted_with_a_tolerance_of_a_millionth_of_a_nanometre(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    ref = X.unpack(raw, 'ref|%s|ISS-0|test' % SET)
    mod = []
    for i, d in enumerate(ref):
        d = d.copy()
        if len(d) and i < 5:
            d[:, :2] += 5e-8                                                                            # below the tolerance
        elif len(d) and i < 10:
            d[:, :2] += 0.5                                                                             # above it
        mod.append(d)
    n_shifted = sum(1 for d in ref[5:10] if len(d))
    off, flat = X.pack(mod)
    raw['g|%s|hp5h|snr5|off' % SET], raw['g|%s|hp5h|snr5|flat' % SET] = off, flat
    md = X.report(raw, 60)
    c = md_row(md, '| N1 LD | a priori default | 60 |')
    assert c[3:] == [str(n_shifted), '0', str(sum(len(d) for d in mod)), str(sum(len(d) for d in ref))] and n_shifted > 0


def test_means_minima_and_maxima_over_two_sets_and_the_status_rows(tmp_path):
    path, _ = planted_raw(str(tmp_path / 'r.npz'), seed=7)
    planted_raw(path, name='MT0.N2.LD', seed=9)
    raw = X.load_raw(path)
    md = X.report(raw, 60)
    evs = [X.SetEval(raw, s, 60) for s in (SET, 'MT0.N2.LD')]
    dE = [e.diff('GR-E', 'ISS-E')['E'][0] for e in evs]
    dJ = [100 * e.diff('GR-E', 'ISS-E')['J'][0] for e in evs]
    dR = [e.diff('GR-E', 'ISS-E')['R'][0] for e in evs]
    assert min(dE) < max(dE)
    c = md_row(md, '| E-tuned operating point | 2 |')
    assert c == ['E-tuned operating point', '2', '%+.2f' % np.mean(dE), '%+.2f / %+.2f' % (min(dE), max(dE)), '%+.2f' % np.mean(dJ), '%+.2f' % np.mean(dR)]
    g = [e.res('R2Eo')[0]['efficiency'] - X.score_set(s, X.dev_dets_of(raw, s, 'R2Eo'), raw['ref|%s|gt_dev' % s], np.arange(30))['efficiency'] for e, s in zip(evs, (SET, 'MT0.N2.LD'))]
    c = md_row(md, '| %s | 2 |' % X.LABEL['R2Eo'], ncells=6)
    assert c[2] == '%+.2f' % np.mean(g) and c[3] == '%+.2f / %+.2f' % (min(g), max(g))
    pe = [e.diff('R2Eo', 'R2E+c')['E'][0] for e in evs]
    c = md_row(md, '| E-tuned | %+.2f |' % np.mean([e.diff('ISS-E', 'R2E+c')['E'][0] for e in evs]))
    assert c[3] == '%+.2f' % np.mean(pe) and c[4] == '%+.2f / %+.2f' % (min(pe), max(pe))
    assert 'Of the 12 verdicts' in md
    path3, _ = planted_raw(str(tmp_path / 'r3.npz'), with_ps=True, with_gb=False, with_g=False)
    assert '| N1 LD | 30-89 (60 frames) | no | no | yes |' in X.report(X.load_raw(path3), 60)
    path4, _ = planted_raw(str(tmp_path / 'r4.npz'), with_ps=False, with_gb=True, with_g=False)
    assert '| N1 LD | 30-89 (60 frames) | no | yes | no |' in X.report(X.load_raw(path4), 60)


def test_resamples_without_a_match_do_not_enter_the_gap_interval(tmp_path):
    path, _ = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    dev = X.unpack(raw, 'ps|%s|R2Jo|dev' % SET)
    sparse = [d if i == 0 else np.zeros((0, 4)) for i, d in enumerate(dev)]                             # matches in one development frame only: most block resamples miss it
    raw['ps|%s|R2Jo|dev|off' % SET], raw['ps|%s|R2Jo|dev|flat' % SET] = X.pack(sparse)
    md = X.report(raw, 100)
    row = md_row(md, '| N1 LD | %s |' % X.LABEL['R2Jo'], ncells=9)
    assert 'nan' not in ' '.join(row)


def test_the_default_bootstrap_is_a_thousand_and_a_part_does_not_replace_a_longer_prefix(tmp_path, monkeypatch):
    path, _ = planted_raw(str(tmp_path / 'controls_raw.npz'), with_ps=False, with_gb=False, with_g=False)
    md_path = str(tmp_path / 'controls.md')
    monkeypatch.setattr(X, 'RAW', path)
    monkeypatch.setattr(X, 'MD', md_path)
    X.main(['--from-raw'])
    assert open(md_path).read() == X.report(X.load_raw(path), 1000) and X.report(X.load_raw(path), 1000) != X.report(X.load_raw(path), 100)
    p = str(tmp_path / 'p.npz')
    X.save_part(p, 'x|1', dict(a=np.arange(2)))
    X.save_part(p, 'x|10', dict(a=np.arange(3)))
    X.save_part(p, 'x|1', dict(b=np.arange(1)))
    assert sorted(X.load_raw(p)) == ['README', 'x|10|a', 'x|1|b']


# ----------------------------------------------------------------------------------------------------------------------------------
# the report on a planted raw
# ----------------------------------------------------------------------------------------------------------------------------------
def planted_dets(gt, idx, rng, bias=(0.0, 0.0), sd=10.0, keep=1.0, width=None):
    out = []
    for f in idx:
        g = gt[gt['frame'] == f]
        k = rng.random(len(g)) < keep
        n = int(k.sum())
        cols = [g['x_nm'][k] + bias[0] + rng.normal(0, sd, n), g['y_nm'][k] + bias[1] + rng.normal(0, sd, n), np.full(n, 1500.0)]
        cols.append(np.full(n, 1.4) if width is None else rng.uniform(*width, n))
        out.append(np.c_[tuple(cols)] if n else np.zeros((0, 4)))
    return out


def planted_raw(path, with_ps=True, with_gb=True, with_g=True, n_test=60, n_dev=30, name=SET, seed=7, f32=False):
    """A raw file for one set: planted ground truth, planted detections of every reference method, greedy detections (before the width filter and the offset), path statistics and validation."""
    rng = np.random.default_rng(seed)
    nf = n_dev + n_test
    rows = []
    for f in range(nf):
        for _ in range(4):
            rows.append((len(rows), f, rng.uniform(300, 6100), rng.uniform(300, 6100), 0.0, rng.uniform(1000, 5000)))
    gt = np.array(rows, dtype=C.GT_DTYPE)
    dev_idx, test_idx = np.arange(n_dev), np.arange(n_dev, nf)
    gt_dev, gt_test = gt[gt['frame'] < n_dev], gt[gt['frame'] >= n_dev]
    ops = dict(E=dict(cfg='cfgE', kappa=5.0, smax=2.5, off=[3.0, -2.0], m=[0.7, 0.99, 0.7, 20.0, 60.0, 0.8, 70.0]),
               J=dict(cfg='cfgJ', kappa=4.5, smax=3.2, off=[1.0, 1.0], m=[0.8, 0.99, 0.8, 25.0, 65.0, 0.85, 66.0]))
    pmin = C.photon_threshold(gt, 0.75)
    A = dict(ops=X.jdump(ops), kind=np.array('snr'), pmin=np.array(pmin), frames=np.array([0, n_dev, n_dev, nf]), gt_dev=gt_dev, gt_test=gt_test)
    spec = {'ISS-E': dict(sd=12.0, keep=0.9), 'ISS-J': dict(sd=18.0, keep=0.95), 'ISS-0': dict(bias=(8.0, -6.0), sd=22.0, keep=0.8), 'R0': dict(bias=(5.0, 5.0), sd=40.0, keep=0.5),
            'R2E+c': dict(sd=30.0, keep=0.75), 'R2J+c': dict(sd=35.0, keep=0.8)}
    for m, kw in spec.items():
        A['%s|test|off' % m], A['%s|test|flat' % m] = X.pack(planted_dets(gt, test_idx, rng, **kw))
        if m in ('ISS-E', 'ISS-J', 'R2E+c', 'R2J+c'):
            A['%s|dev|off' % m], A['%s|dev|flat' % m] = X.pack(planted_dets(gt, dev_idx, rng, **dict(kw, sd=kw['sd'] * 0.8)))
    X.save_part(path, 'ref|' + name, A)
    X.save_part(path, 'meta', dict(sha=X.jdump({'iss_raw.npz': 'a' * 64, 'picasso_ladder_raw.npz': 'b' * 64}), versions=X.jdump(dict(python='3', numpy='2'))))
    if with_g:
        G = {'frames': np.array([n_dev, nf])}
        pk = X.E13.pack if f32 else X.pack                                                                  # the 2013 pipeline stores float32
        for cfg, kap, sd in (('cfgE', 5.0, 40.0), ('cfgJ', 4.5, 6.0), ('hp5h', 5.0, 60.0)):
            for k in X.ISS.KAPPAS:
                ds = planted_dets(gt, test_idx, np.random.default_rng(int(k * 10) + len(cfg)), bias=(3.0, -2.0), sd=sd, keep=0.85, width=(0.9, 3.0))
                if f32 and cfg == 'cfgJ':
                    ds[0] = np.vstack([ds[0], [[1000.0, 1000.0, 1500.0, 3.2]]])                              # a width at the bound, rounded to float32 (3.2000000477): kept by the 2013 tolerance only
                G['%s|snr%g|off' % (cfg, k)], G['%s|snr%g|flat' % (cfg, k)] = pk(ds)
            G['%s|secs' % cfg], G['%s|dual_secs' % cfg] = np.full(60, 0.5), np.full(60, 0.25)
            G['%s|stats' % cfg] = np.array([60, 200, 4, 2, 0, 0, 1])
            G['%s|dual_stats' % cfg] = np.array([60, 190, 3, 1, 0, 0, 0])
        X.save_part(path, 'g|' + name, G)
        X.save_part(path, 'val|' + name, dict(dual_original_stop=X.jdump(dict(set=name, cfg='cfgE', frames=60, rules=7, identical=True, rules_identical=7)),
                                              dual_truncated_stop=X.jdump(dict(set=name, cfg='cfgE', frames=60, rules=7, identical=False, rules_identical=6)),
                                              guard=X.jdump(dict(set=name, cfg='cfgE', frames=20, identical=True, barred=0))))
    if with_gb:
        ops_g = dict(E=dict(cfg='cfgG', kappa=5.0, smax=3.2, off=[2.0, 0.0], m=[0.7, 0.99, 0.7, 20.0, 60.0, 0.8, 70.0]), J=dict(cfg='cfgG', kappa=4.5, smax=3.2, off=[0.0, 0.0], m=[0.8, 0.99, 0.8, 25.0, 65.0, 0.85, 66.0]))
        B = dict(ops=X.jdump(ops_g), frames=np.array([0, n_dev, n_dev, nf]))
        for k in (5.0, 4.5):
            B['T|cfgG|snr%g|off' % k], B['T|cfgG|snr%g|flat' % k] = X.pack(planted_dets(gt, test_idx, np.random.default_rng(int(k * 7)), bias=(2.0, 0.0), sd=16.0, keep=0.9, width=(0.9, 3.0)))
            B['D|cfgG|snr%g|off' % k], B['D|cfgG|snr%g|flat' % k] = X.pack(planted_dets(gt, dev_idx, np.random.default_rng(int(k * 9)), bias=(2.0, 0.0), sd=12.0, keep=0.9, width=(0.9, 3.0)))
        X.save_part(path, 'gb|' + name, B)
    if with_ps:
        P = {}
        for r, sd in (('R2Eo', 24.0), ('R2Jo', 28.0)):
            for split, idx in (('dev', dev_idx), ('test', test_idx)):
                P['%s|%s|off' % (r, split)], P['%s|%s|flat' % (r, split)] = X.pack(planted_dets(gt, idx, np.random.default_rng(len(r) + len(split)), sd=sd * (0.8 if split == 'dev' else 1.0), keep=0.85))
            P['%s|config' % r] = X.jdump(C._merge_params(dict(min_net_gradient=1200.0 if r == 'R2Eo' else 800.0, box=7)))
            P['%s|offset' % r], P['%s|key' % r] = np.array([5.0, -4.0]), np.array('box=7|fit_box=7|gauss_sigma=None|pad=0|fit=gausslq')
        res = {'box=7|fit_box=7|gauss_sigma=None|pad=0|fit=gausslq': dict(struct={}, n_cand=1000, n_eval=120, E=dict(metrics=dict(efficiency=55.5, jaccard=0.66)), J=dict(metrics=dict(efficiency=40.0, jaccard=0.71)))}
        P.update(search=X.jdump(res), stored_config=X.jdump({'R2E': C._merge_params(dict(min_net_gradient=1000.0, box=7)), 'R2J': C._merge_params(dict(min_net_gradient=800.0, box=7, fit_box=9))}),
                 stored_corrected=X.jdump({'R2E': dict(efficiency=54.0, jaccard=0.6), 'R2J': dict(efficiency=40.0, jaccard=0.69)}),
                 validation=X.jdump({'R2E': dict(identical_to_rounding=True, bit_identical=True, detections=5), 'R2J': dict(identical_to_rounding=True, bit_identical=False, detections=5)}),
                 summary=X.jdump(dict(info={'R2Eo': dict(key='box=7|fit_box=7|gauss_sigma=None|pad=0|fit=gausslq', objective='E', offset=[5.0, -4.0]),
                                            'R2Jo': dict(key='box=7|fit_box=7|gauss_sigma=None|pad=0|fit=gausslq', objective='J', offset=[5.0, -4.0])}, structs=1, n_eval=120, n_cand=1000, secs=1.0, max_structs=None)))
        X.save_part(path, 'ps|' + name, P)
    return path, gt


def test_the_dets_of_a_method_are_filtered_by_the_width_bound_and_shifted_by_the_operating_points_offset(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    base = X.unpack(raw, 'g|%s|cfgE|snr5' % SET)
    got = X.dets_of(raw, SET, 'GR-E')
    assert len(got) == 60 and any(len(b) > len(g) for b, g in zip(base, got))                           # the width bound (2.5 px) removes some
    for b, g in zip(base, got):
        keep = b[:, 3] <= 2.5 + 1e-9
        assert np.allclose(g[:, 0], b[keep, 0] - 3.0) and np.allclose(g[:, 1], b[keep, 1] + 2.0) and np.allclose(g[:, 2:], b[keep, 2:])
    base_j = X.unpack(raw, 'g|%s|cfgJ|snr4.5' % SET)
    gj = X.dets_of(raw, SET, 'GR-J')
    assert all(np.allclose(g[:, 0], b[b[:, 3] <= 3.2 + 1e-9, 0] - 1.0) for b, g in zip(base_j, gj) if len(g))
    d0 = X.dets_of(raw, SET, 'GR-0')
    b0 = X.unpack(raw, 'g|%s|hp5h|snr5' % SET)
    assert all(np.array_equal(a, b) for a, b in zip(d0, b0))                                            # the a priori default: no width bound, no offset
    gt_e = X.dets_of(raw, SET, 'GT-E')
    bg = X.unpack(raw, 'gb|%s|T|cfgG|snr5' % SET)
    assert all(np.allclose(g[:, 0], b[b[:, 3] <= 3.2 + 1e-9, 0] - 2.0) for b, g in zip(bg, gt_e) if len(g))
    gt_j = X.dets_of(raw, SET, 'GT-J')
    bj = X.unpack(raw, 'gb|%s|T|cfgG|snr4.5' % SET)
    assert all(np.allclose(g[:, 0], b[b[:, 3] <= 3.2 + 1e-9, 0]) for b, g in zip(bj, gt_j) if len(g)) and sum(len(g) for g in gt_j) == sum(int((b[:, 3] <= 3.2 + 1e-9).sum()) for b in bj)      # the J point: kappa 4.5, no offset
    assert sum(len(g) for g in gt_j) != sum(len(g) for g in gt_e)
    dev_e = X.dev_dets_of(raw, SET, 'GT-E')
    bd = X.unpack(raw, 'gb|%s|D|cfgG|snr5' % SET)
    assert all(np.allclose(g[:, 0], b[b[:, 3] <= 3.2 + 1e-9, 0] - 2.0) for b, g in zip(bd, dev_e) if len(g)) and len(dev_e) == 30
    assert X.dets_of(raw, SET, 'ISS-E') is not None and X.dets_of(raw, SET, 'R2Eo') is not None and X.dev_dets_of(raw, SET, 'R2Jo') is not None and X.dev_dets_of(raw, SET, 'ISS-0') is None
    assert X.dets_of(raw, 'Tubulins_II', 'GR-E') is None and X.dets_of(raw, SET, 'nothing') is None
    # a raw without the greedy and Picasso parts has no such methods
    path2, _ = planted_raw(str(tmp_path / 'r2.npz'), with_ps=False, with_gb=False, with_g=False)
    raw2 = X.load_raw(path2)
    assert X.dets_of(raw2, SET, 'GR-E') is None and X.dets_of(raw2, SET, 'GT-E') is None and X.dets_of(raw2, SET, 'R2Eo') is None and X.dev_dets_of(raw2, SET, 'GT-E') is None


def test_a_2013_set_widens_float32_detections_and_keeps_a_width_at_the_bound(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r13.npz'), name=SET13, f32=True)
    raw = X.load_raw(path)
    assert raw['g|%s|cfgE|snr5|flat' % SET13].dtype == np.float32 and raw['ref|%s|ISS-E|test|flat' % SET13].dtype == np.float64
    got = X.dets_of(raw, SET13, 'GR-J')
    base = X.base_dets(raw, SET13, 'g|%s|cfgJ|snr4.5' % SET13)
    assert base[0].dtype == np.float64 and got[0].dtype == np.float64
    keep_13 = sum(int((b[:, 3] <= 3.2 + 1e-4).sum()) for b in base)
    keep_16 = sum(int((b[:, 3] <= 3.2 + 1e-9).sum()) for b in base)
    assert keep_13 == keep_16 + 1 and sum(len(g) for g in got) == keep_13                               # the 2013 tolerance keeps the float32 width 3.2000000477; the 2016 one would drop it
    assert any(abs(g[:, 0] - (1000.0 - 1.0)).min() < 1e-9 for g in got[:1] if len(g)) and X.post_of(SET13) is X.E13.post and X.post_of(SET) is X.ISS.post
    md = X.report(raw, 60)
    assert '| Tubulins I |' in md and '## A1a.' in md and 'Tubulins I | cfgJ kappa 4.5 width <= 3.2' in md


def md_row(md, prefix, nth=0, ncells=None):
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in md.split('\n') if l.startswith(prefix)]
    return [r for r in rows if ncells is None or len(r) == ncells][nth]


def test_the_tables_of_the_report_are_the_computed_ones(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    md = X.report(raw, 100)
    ev = X.SetEval(raw, SET, 100)
    cv = lambda a, b, key, hb=True, scale=1.0: X.dcell(ev.diff(a, b)[key], scale) + ' ' + X.verdict(ev.diff(a, b)[key], hb)
    # the J-tuned and default tables: all ten cells
    c = md_row(md, '| N1 LD | cfgJ kappa 4.5 width <= 3.2 |')
    assert c[2:] == [X.metric_cell(ev, 'ISS-J', 'efficiency'), X.metric_cell(ev, 'GR-J', 'efficiency'), cv('GR-J', 'ISS-J', 'E'), cv('GR-J', 'ISS-J', 'J', True, 100.0), cv('GR-J', 'ISS-J', 'R', False),
                     cv('GR-J', 'ISS-J', 'Ec'), cv('GR-J', 'R2J+c', 'E'), cv('ISS-J', 'R2J+c', 'E')]
    c = md_row(md, '| N1 LD | hp5h kappa 5 |')
    assert c[2:] == [X.metric_cell(ev, 'ISS-0', 'efficiency'), X.metric_cell(ev, 'GR-0', 'efficiency'), cv('GR-0', 'ISS-0', 'E'), cv('GR-0', 'ISS-0', 'J', True, 100.0), cv('GR-0', 'ISS-0', 'R', False),
                     cv('GR-0', 'ISS-0', 'Ec'), cv('GR-0', 'R0', 'E'), cv('ISS-0', 'R0', 'E')]
    c = md_row(md, '| N1 LD | cfgE kappa 5 width <= 2.5 |')
    assert c[9] == cv('ISS-E', 'R2E+c', 'E') and c[7] == cv('GR-E', 'ISS-E', 'Ec')
    # the record over the sets, from the harness's tally
    for w, p in (('E', 'R2E+c'), ('J', 'R2J+c'), ('0', 'R0')):
        t = X.op_title(w)
        for lab, a, b in (('greedy vs dual', 'GR-' + w, 'ISS-' + w), ('greedy vs ' + p, 'GR-' + w, p), ('dual vs ' + p, 'ISS-' + w, p)):
            c = md_row(md, '| %s | %s |' % (t, lab))
            assert c[2:] == ['%d / %d / %d' % X.tally([(ev, a, b)], key, hb) for key, hb in (('E', True), ('J', True), ('R', False))], (w, lab)
        dE, dJ, dR = (ev.diff('GR-' + w, 'ISS-' + w)[k][0] for k in ('E', 'J', 'R'))
        c = md_row(md, '| %s | 1 | %+.2f' % (t, dE))
        assert c == [t, '1', '%+.2f' % dE, '%+.2f / %+.2f' % (dE, dE), '%+.2f' % (100 * dJ), '%+.2f' % dR]
        g, d = X.dets_of(raw, SET, 'GR-' + w), X.dets_of(raw, SET, 'ISS-' + w)
        ndiff = sum(1 for x, y in zip(g, d) if not (x.shape == y.shape and np.abs(x[:, :2] - y[:, :2]).max(initial=0.0) <= 1e-6))
        ncnt = sum(1 for x, y in zip(g, d) if len(x) != len(y))
        c = md_row(md, '| N1 LD | %s | 60 |' % t)
        assert c[3:] == [str(ndiff), str(ncnt), str(sum(len(x) for x in g)), str(sum(len(y) for y in d))] and ndiff >= ncnt > 0
    # the timing table is empty for the planted raw (no seconds were stored), A2 rows
    assert md_row(md, '| N1 LD | cfgE | 0.500 ') == ['N1 LD', 'cfgE', '0.500', '', '0.250', '60'] and 'iss_raw.npz aaaaaaaaaaaa' in md and 'picasso_ladder_raw.npz bbbbbbbbbbbb' in md
    for w in ('E', 'J'):
        c = md_row(md, '| N1 LD | %s-tuned | ' % w)
        assert c[2:] == [cv('ISS-' + w, X.PIC_OF[w], 'E'), cv('ISS-' + w, X.PICO_OF[w], 'E'), cv('ISS-' + w, X.PIC_OF[w], 'J', True, 100.0), cv('ISS-' + w, X.PICO_OF[w], 'J', True, 100.0),
                         cv('ISS-' + w, X.PIC_OF[w], 'R', False), cv('ISS-' + w, X.PICO_OF[w], 'R', False), cv('ISS-' + w, X.PIC_OF[w], 'Ec'), cv('ISS-' + w, X.PICO_OF[w], 'Ec')]
        c = md_row(md, '| N1 LD | %s | ' % w, ncells=5)
        assert c[2:] == [cv(X.PICO_OF[w], X.PIC_OF[w], 'E'), cv(X.PICO_OF[w], X.PIC_OF[w], 'J', True, 100.0), cv(X.PICO_OF[w], X.PIC_OF[w], 'R', False)]
        for key, lab, hb in (('E', 'E', True), ('J', 'J', True), ('R', 'RMSE', False), ('Ec', 'E (challenge)', True)):
            c = md_row(md, '| %s-tuned ISS | %s |' % (w, lab))
            assert c[2:] == ['%d / %d / %d' % X.tally([(ev, 'ISS-' + w, X.PIC_OF[w])], key, hb), '%d / %d / %d' % X.tally([(ev, 'ISS-' + w, X.PICO_OF[w])], key, hb)]
    a = ev.diff('ISS-E', 'R2E+c')['E'][0]
    b = ev.diff('ISS-E', 'R2Eo')['E'][0]
    pe = ev.diff('R2Eo', 'R2E+c')['E'][0]
    assert md_row(md, '| E-tuned | %+.2f |' % a) == ['E-tuned', '%+.2f' % a, '%+.2f' % b, '%+.2f' % pe, '%+.2f / %+.2f' % (pe, pe)]
    # the selection table: validation, configurations, dev values
    c = md_row(md, '| N1 LD | R2Eo |')
    assert c[2] == 'yes' and c[4] == '54.0' and c[6] == '55.5' and c[7] == '(5.0, -4.0)' and c[8] == 'different' and 'net gradient > 1000' in c[3] and 'net gradient > 1200' in c[5]
    c = md_row(md, '| N1 LD | R2Jo |')
    assert c[2] == 'to rounding' and c[4] == '69.0' and c[6] == '71.0'
    # the gap table and its means
    rd = {m: X.score_set(SET, X.dev_dets_of(raw, SET, m), raw['ref|%s|gt_dev' % SET], np.arange(30)) for m in ('ISS-E', 'ISS-J', 'R2E+c', 'R2J+c', 'R2Eo', 'R2Jo', 'GT-E', 'GT-J')}
    es, js = {}, {}
    for m, r in rd.items():
        rt = ev.res(m)[0]
        es[m], js[m] = rt['efficiency'] - r['efficiency'], 100 * (rt['jaccard'] - r['jaccard'])
        c = md_row(md, '| N1 LD | %s |' % X.LABEL[m], ncells=9)
        assert c[2:4] == ['%.1f' % r['efficiency'], '%.1f' % rt['efficiency']] and c[5:7] == ['%.1f' % (100 * r['jaccard']), '%.1f' % (100 * rt['jaccard'])]
        c = md_row(md, '| %s | 1 |' % X.LABEL[m], ncells=6)
        assert c == [X.LABEL[m], '1', '%+.2f' % es[m], '%+.2f / %+.2f' % (es[m], es[m]), '%+.2f' % js[m], '%+.2f / %+.2f' % (js[m], js[m])]
    # A1b rows
    c = md_row(md, '| N1 LD | E | cfgE kappa 5 width <= 2.5 |')
    assert c[3] == 'cfgG kappa 5 width <= 3.2' and c[6] == cv('GT-E', 'ISS-E', 'E') and c[7] == cv('GT-E', 'ISS-E', 'J', True, 100.0) and c[8] == cv('GT-E', 'ISS-E', 'R', False) and c[9] == cv('GT-E', 'R2E+c', 'E')
    c = md_row(md, '| N1 LD | J | cfgJ kappa 4.5 width <= 3.2 |')
    assert c[6] == cv('GT-J', 'ISS-J', 'E') and c[9] == cv('GT-J', 'R2J+c', 'E') and c[3] == 'cfgG kappa 4.5 width <= 3.2'
    for w in ('E', 'J'):
        c = md_row(md, '| %s-tuned | greedy (own tuning) vs Picasso R2*+c |' % w)
        assert c[2:] == ['%d / %d / %d' % X.tally([(ev, 'GT-' + w, X.PIC_OF[w])], key, hb) for key, hb in (('E', True), ('J', True), ('R', False))]
    # the verdicts that change between the stored and the symmetric Picasso row
    ch = []
    for w in ('E', 'J'):
        for key, nm, hb in (('E', 'E', True), ('J', 'J', True), ('R', 'RMSE', False)):
            va, vb = X.verdict(ev.diff('ISS-' + w, X.PIC_OF[w])[key], hb), X.verdict(ev.diff('ISS-' + w, X.PICO_OF[w])[key], hb)
            if va != vb:
                ch.append('N1 LD %s-tuned %s: %s -> %s' % (w, nm, va, vb))
    assert ('1 change' if len(ch) == 1 else '%d change' % len(ch)) in md and all(ch) and (('; '.join(ch) + '.') in md if ch else 'Of the 6 verdicts (dual ISS against Picasso: 1 sets x E-tuned / J-tuned x E, J, RMSE), 0 change' in md)
    # a gap interval, re-derived from the two bootstraps
    xd = X.boot_reps(X.score_set(SET, X.dev_dets_of(raw, SET, 'ISS-E'), raw['ref|%s|gt_dev' % SET], np.arange(30))['per_frame'], 100, 1)
    xt = X.boot_reps(ev.res('ISS-E')[0]['per_frame'], 100, 2)
    dlt = xt['efficiency'] - xd['efficiency']
    dlt = dlt[~np.isnan(dlt)]
    gl = md_row(md, '| N1 LD | dual ISS, E-tuned |', ncells=9)
    assert gl[4] == '%+.1f [%+.1f, %+.1f]' % (ev.res('ISS-E')[0]['efficiency'] - X.score_set(SET, X.dev_dets_of(raw, SET, 'ISS-E'), raw['ref|%s|gt_dev' % SET], np.arange(30))['efficiency'], np.percentile(dlt, 2.5), np.percentile(dlt, 97.5))
    c = md_row(md, '| E-tuned | greedy (own tuning) vs dual |')
    assert c[2:] == ['%d / %d / %d' % X.tally([(ev, 'GT-E', 'ISS-E')], key, hb) for key, hb in (('E', True), ('J', True), ('R', False))]


def test_the_set_evaluation_scores_the_test_frames_and_the_pairs_with_the_harness(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    ev = X.SetEval(raw, SET, 200)
    gt_test = raw['ref|%s|gt_test' % SET]
    assert ev.idx.tolist() == list(range(30, 90)) and ev.has('GR-E') and not ev.has('nothing')
    for m in ('ISS-E', 'GR-E', 'R2Eo'):
        r, r2 = ev.res(m)
        ref = C.score(X.dets_of(raw, SET, m), gt_test, ev.idx)
        assert r['efficiency'] == ref['efficiency'] and r['tp'] == ref['tp']
        ref2 = C.score(X.dets_of(raw, SET, m), gt_test, ev.idx, border_nm=450.0, photon_min=float(raw['ref|%s|pmin' % SET]))
        assert r2['efficiency'] == ref2['efficiency'] and r2['tp'] == ref2['tp'] and r2['tp'] < r['tp']
    d = ev.diff('GR-E', 'ISS-E')
    a, b = ev.res('GR-E')[0], ev.res('ISS-E')[0]
    assert d['E'][0] == a['efficiency'] - b['efficiency'] and d['J'][0] == a['jaccard'] - b['jaccard'] and d['R'][0] == a['rmse_nm'] - b['rmse_nm']
    ci = C.bootstrap_diff_ci(a['per_frame'], b['per_frame'], B=200, seed=0, block=10)
    assert d['E'][1] == ci['efficiency'] and d['R'][1] == ci['rmse_nm'] and ev.diff('GR-E', 'ISS-E') is d
    assert ev.ci('ISS-E') == C.bootstrap_ci(ev.res('ISS-E')[0], B=200, seed=0, block=10)
    a2, b2 = ev.res('GR-E')[1], ev.res('ISS-E')[1]
    ci2 = C.bootstrap_diff_ci(a2['per_frame'], b2['per_frame'], B=200, seed=0, block=10)
    assert d['Ec'][0] == a2['efficiency'] - b2['efficiency'] and d['Ec'][1] == ci2['efficiency'] and d['J'][1] == ci['jaccard']
    r_e = ev.res('ISS-E')[0]
    assert X.metric_cell(ev, 'ISS-E', 'efficiency') == X.fmt_ci(r_e['efficiency'], ev.ci('ISS-E')['efficiency'], 1.0) and X.metric_cell(ev, 'ISS-E', 'jaccard', 100.0) == X.fmt_ci(r_e['jaccard'], ev.ci('ISS-E')['jaccard'], 100.0)
    n_calls = []
    real = X.dets_of
    X.dets_of = lambda *a: (n_calls.append(1), real(*a))[1]
    try:
        ev2 = X.SetEval(raw, SET, 50)
        assert ev2.has('GR-E') and ev2.has('GR-E') and not ev2.has('nothing') and not ev2.has('nothing')
    finally:
        X.dets_of = real
    assert len(n_calls) == 2                                                                              # the existence of a method is looked up once
    with pytest.raises(ValueError, match='59 detection lists for 60 test frames'):                      # a method with the wrong number of frames is refused
        raw_bad = dict(raw)
        raw_bad['ref|%s|ISS-E|test|off' % SET] = raw['ref|%s|ISS-E|test|off' % SET][:-1]
        X.SetEval(raw_bad, SET, 200).res('ISS-E')
    assert X.verdict((1.0, (0.5, 1.5))) == 'win' and X.verdict((1.0, (-0.5, 1.5))) == 'tie' and X.verdict((-1.0, (-1.5, -0.5))) == 'loss'
    assert X.verdict((-1e-14, (-1e-14, -1e-14))) == 'tie' and X.verdict((1e-14, (1e-14, 1e-14)), False) == 'tie' and X.verdict((2e-4, (1.5e-4, 3e-4))) == 'win' and X.verdict((-2e-4, (-3e-4, -1.5e-4))) == 'loss' and X.verdict((3e-5, (2e-5, 4e-5))) == 'tie'
    assert X.verdict((1.0, (0.5, 1.5)), False) == 'loss' and X.verdict((-1.0, (-1.5, -0.5)), False) == 'win'
    assert sum(X.tally([(ev, 'GR-E', 'ISS-E'), (ev, 'ISS-E', 'R0'), (ev, 'nothing', 'R0')], 'E')) == 2 and X.tally([(ev, 'nothing', 'R0')], 'E') == (0, 0, 0)


class StubEv:
    def has(self, m):
        return m != 'missing'

    def diff(self, A, B):
        return {'E': {'win': (1.0, (0.5, 1.5)), 'tie': (0.1, (-0.5, 0.7)), 'loss': (-1.0, (-1.5, -0.5))}[B]}


def test_the_tally_counts_wins_ties_and_losses_of_existing_methods():
    ev = StubEv()
    assert X.tally([(ev, 'x', 'win'), (ev, 'x', 'tie'), (ev, 'x', 'tie'), (ev, 'x', 'loss'), (ev, 'x', 'loss'), (ev, 'x', 'loss')], 'E') == (1, 2, 3)
    assert X.tally([(ev, 'x', 'win'), (ev, 'missing', 'win'), (ev, 'x', 'missing')], 'E') == (1, 0, 0)
    assert X.tally([(ev, 'x', 'win'), (ev, 'x', 'loss')], 'E', False) == (1, 0, 1)


def test_the_gap_intervals_come_from_independent_bootstraps_of_the_two_ranges():
    cands, gt, idx = ascent_world()
    r = C.score(C.filter_detections(cands, {}), gt, idx)
    reps = X.boot_reps(r['per_frame'], 300, 1)
    ci = C.bootstrap_ci(r, B=300, seed=1, block=10)
    assert np.allclose(np.percentile(reps['efficiency'][~np.isnan(reps['efficiency'])], [2.5, 97.5]), ci['efficiency']) and np.allclose(np.percentile(reps['jaccard'], [2.5, 97.5]), ci['jaccard'])
    assert set(reps) == {'jaccard', 'efficiency', 'rmse_nm'} and not np.array_equal(reps['efficiency'], X.boot_reps(r['per_frame'], 300, 2)['efficiency']) and X.GAP_SEEDS == (1, 2)


def test_the_report_has_every_section_and_its_numbers_are_the_harnesss(tmp_path):
    path, gt = planted_raw(str(tmp_path / 'r.npz'))
    raw = X.load_raw(path)
    md = X.report(raw, 100)
    for head in ('## What is in this report', '## Harness validation', '## A1a. Greedy control', '### E-tuned operating point', '### J-tuned operating point', '### A priori default', '### Record over the sets',
                 '### How often the two paths give different detections', '### Guard activations and cycles', '## A1b. Greedy tuned on its own development frames', '## A2. Picasso objective symmetry',
                 '### Dual ISS minus Picasso, with the stored selection and with the symmetric one', '### Did the conclusions change?', '## A3. Development-versus-test gap', '## What this does not establish'):
        assert head in md, head
    ev = X.SetEval(raw, SET, 100)
    # a row of the E-tuned table: the E of the dual and of the greedy path, their paired difference
    line = [l for l in md.split('\n') if l.startswith('| N1 LD | cfgE kappa 5 width <= 2.5 |')][0]
    cells = [c.strip() for c in line.strip('|').split('|')]
    assert cells[2] == X.metric_cell(ev, 'ISS-E', 'efficiency') and cells[3] == X.metric_cell(ev, 'GR-E', 'efficiency')
    d = ev.diff('GR-E', 'ISS-E')
    assert cells[4] == X.dcell(d['E']) + ' ' + X.verdict(d['E']) and cells[5] == X.dcell(d['J'], 100.0) + ' ' + X.verdict(d['J']) and cells[6] == X.dcell(d['R']) + ' ' + X.verdict(d['R'], False)
    dp = ev.diff('GR-E', 'R2E+c')
    assert cells[8] == X.dcell(dp['E']) + ' ' + X.verdict(dp['E'])
    # the validation table, the path statistics and the symmetric Picasso rows
    assert '| N1 LD | dual path through this harness, original stop | cfgE | 60 | 7 of 7 | yes |' in md and '| N1 LD | dual path through this harness, truncated stop | cfgE | 60 | 6 of 7 | NO |' in md
    assert '| N1 LD | cfgE | greedy | 60 | 200 | 4 | 2 | 0 | 0 | 1 |' in md and '| N1 LD | cfgE | dual | 60 | 190 | 3 | 1 | 0 | 0 | 0 |' in md
    assert 'Picasso R2Eo (offset-corrected objective)' in md and '| N1 LD | R2Eo | yes |' in md and '| N1 LD | R2Jo | to rounding |' in md
    sec = md[md.index('### Did the conclusions change?'):]
    assert 'Of the 6 verdicts' in sec
    # A3: the rows of the gap table are the harness's dev and test values
    gline = [l for l in md.split('\n') if l.startswith('| N1 LD | dual ISS, E-tuned |')][0]
    gc = [c.strip() for c in gline.strip('|').split('|')]
    rd = C.score(X.dev_dets_of(raw, SET, 'ISS-E'), raw['ref|%s|gt_dev' % SET], np.arange(30))
    rt = ev.res('ISS-E')[0]
    assert gc[2] == '%.1f' % rd['efficiency'] and gc[3] == '%.1f' % rt['efficiency'] and gc[5] == '%.1f' % (100 * rd['jaccard']) and gc[6] == '%.1f' % (100 * rt['jaccard']) and gc[8] == '30 / 60'
    assert gc[4].startswith('%+.1f [' % (rt['efficiency'] - rd['efficiency']))


def test_a_raw_without_the_optional_parts_gives_a_report_that_says_so(tmp_path):
    path, _ = planted_raw(str(tmp_path / 'r.npz'), with_ps=False, with_gb=False, with_g=False)
    md = X.report(X.load_raw(path), 100)
    assert 'The greedy part is not in the raw.' in md and 'Not run for any set (time-box).' in md and '| N1 LD | 1000-1019' not in md
    assert '| N1 LD | 30-89 (60 frames) | no | no | no |' in md and 'No validation run in the raw.' in md
    assert X.report({}, 100).startswith('# Controls') and 'nothing to report' in X.report({}, 100)
    path3, _ = planted_raw(str(tmp_path / 'r3.npz'), with_ps=True, with_gb=False, with_g=False)
    md3 = X.report(X.load_raw(path3), 100)
    assert '## A2.' in md3 and 'dE stored' in md3 and 'Not run for any set (time-box).' in md3                    # A2 present, A1b absent


def test_the_report_regenerates_byte_for_byte_from_the_raw_and_the_command_line_uses_the_output_variable(tmp_path, monkeypatch):
    path, _ = planted_raw(str(tmp_path / 'controls_raw.npz'))
    md_path = str(tmp_path / 'controls.md')
    monkeypatch.setattr(X, 'RAW', path)
    monkeypatch.setattr(X, 'MD', md_path)
    X.main(['--from-raw', '--boot', '100'])
    first = open(md_path, 'rb').read()
    X.main(['--from-raw', '--boot', '100'])
    assert open(md_path, 'rb').read() == first and first == X.report(X.load_raw(path), 100).encode()
    with pytest.raises(SystemExit):
        X.main([])
    assert os.path.basename(X.__dict__['MD']) == 'controls.md'.replace('controls.md', os.path.basename(md_path))
    src = open(X.__file__).read()
    assert "os.environ.get('CONTROLS_OUT'" in src and "RAW = os.path.join(OUT, 'controls_raw.npz')" in src and "MD = os.path.join(OUT, 'controls.md')" in src


def test_the_raw_file_is_the_same_bytes_for_the_same_arrays_and_parts_replace_their_own_keys(tmp_path):
    p1, p2 = str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz')
    X.save_part(p1, 'x|1', dict(a=np.arange(3), b=np.zeros(2)))
    X.save_part(p1, 'y|1', dict(a=np.ones(2)))
    X.save_part(p2, 'y|1', dict(a=np.ones(2)))
    X.save_part(p2, 'x|1', dict(b=np.zeros(2), a=np.arange(3)))
    assert open(p1, 'rb').read() == open(p2, 'rb').read()
    X.save_part(p1, 'x|1', dict(c=np.arange(4)))                                                         # replaces the keys of its prefix, keeps the others
    raw = X.load_raw(p1)
    assert sorted(raw) == ['README', 'x|1|c', 'y|1|a']
    X.save_part(p1, 'x|1', dict(d=np.arange(2)), replace=False)
    assert sorted(X.load_raw(p1)) == ['README', 'x|1|c', 'x|1|d', 'y|1|a'] and X.load_raw(str(tmp_path / 'nothing.npz')) == {}
    assert X.jload(X.jdump({'b': np.float64(1.5), 'a': [np.int64(2)]})) == {'a': [2], 'b': 1.5} and str(X.jdump({'b': 1, 'a': 2})) == '{"a": 2, "b": 1}'


def test_small_helpers_cut_the_frames_name_the_families_and_pick_the_configurations():
    assert X.cut([1, 2, 3], 2) == [1, 2] and X.cut([1, 2, 3], None) == [1, 2, 3]
    assert X.fam('MT0.N2.HD') == '16' and X.fam('Tubulins_II') == '13' and len(X.SETS) == 8 and X.SETS16 + X.SETS13 == X.SETS
    assert X.cfgs_needed({'E': {'cfg': 'b'}, 'J': {'cfg': 'a'}}, 'hp5h') == ['a', 'b', 'hp5h'] and X.cfgs_needed({'E': {'cfg': 'hp5h'}, 'J': {'cfg': 'hp5h'}}, 'hp5h') == ['hp5h']
    o = dict(cfg='c', kappa=np.float64(5.0), smax=2, off=np.array([1.0, 2.0]), m=np.arange(7))
    assert X.ops_json({'E': o}) == {'E': dict(cfg='c', kappa=5.0, smax=2.0, off=[1.0, 2.0], m=[float(i) for i in range(7)])}
    assert X.equal_dets([np.ones((1, 2))], [np.ones((1, 2))]) and not X.equal_dets([np.ones((1, 2))], [np.ones((1, 2)) * 2]) and not X.equal_dets([np.ones((1, 2))], [])
    assert X.PIC_OF == {'E': 'R2E+c', 'J': 'R2J+c', '0': 'R0'} and X.PICO_OF == {'E': 'R2Eo', 'J': 'R2Jo'} and X.OPS == ('E', 'J', '0') and X.BLOCK == 10 and X.BOOT == 1000
    assert set(X.LABEL) >= set(X.PIC_OF.values()) | set(X.PICO_OF.values()) | {'ISS-' + w for w in X.OPS} | {'GR-' + w for w in X.OPS} | {'GT-E', 'GT-J'}
    assert X.score_set(SET, [np.zeros((0, 2))], np.zeros(0, dtype=C.GT_DTYPE), np.array([0]))['n_det'] == 0
    assert X.op_title('0') == 'a priori default' and X.op_title('E').startswith('E-tuned')


# ----------------------------------------------------------------------------------------------------------------------------------
# the scorer, re-derived independently
# ----------------------------------------------------------------------------------------------------------------------------------
def test_an_exhaustive_pairing_reproduces_the_harness_scores_on_planted_detections():
    rng = np.random.default_rng(31)
    rows = []
    for f in range(12):
        for _ in range(rng.integers(1, 5)):
            rows.append((len(rows), f, rng.uniform(1000, 2200), rng.uniform(1000, 2200), 0.0, 2000.0))
    gt = np.array(rows, dtype=C.GT_DTYPE)
    idx = np.arange(12)
    dets = planted_dets(gt, idx, rng, bias=(10.0, -5.0), sd=70.0, keep=0.8)
    for f in (3, 7):
        dets[f] = np.vstack([dets[f], [[1500.0, 1500.0, 1.0, 1.4], [1560.0, 1500.0, 1.0, 1.4]]])          # two extra detections close together (contested emitters, false positives)
    mine = brute_score(dets, gt, idx)
    ref = C.score(dets, gt, idx)
    assert abs(mine['jaccard'] - ref['jaccard']) < 1e-12 and abs(mine['rmse_nm'] - ref['rmse_nm']) < 1e-9 and abs(mine['efficiency'] - ref['efficiency']) < 1e-9
    assert ref['fp'] > 0 and ref['fn'] > 0


def stored_raw_path():
    return os.path.join(os.environ.get('CONTROLS_OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'smlm')), 'controls_raw.npz')


@pytest.mark.skipif(not os.path.exists(stored_raw_path()), reason='no controls_raw.npz')
def test_stored_rows_rescored_by_an_exhaustive_pairing_agree_with_the_reports_scorer():
    raw = X.load_raw(stored_raw_path())
    done = 0
    for name in X.SETS16:
        if 'ref|%s|ops' % name not in raw:
            continue
        a, b = (int(v) for v in raw['ref|%s|frames' % name][2:4])
        idx = np.arange(a, a + 25)
        gt = raw['ref|%s|gt_test' % name]
        for m in ('ISS-E', 'R2E+c', 'GR-E', 'R2Eo'):
            ds = X.dets_of(raw, name, m)
            if ds is None:
                continue
            sub = ds[:25]
            ng = max((int((gt['frame'] == f).sum()) for f in idx), default=0)
            if max((len(d) for d in sub), default=0) > 5 or ng > 8:
                continue
            ref = C.score(sub, gt, idx)
            mine = brute_score(sub, gt, idx)
            assert abs(mine['jaccard'] - ref['jaccard']) < 1e-12 and abs(mine['rmse_nm'] - ref['rmse_nm']) < 1e-9 and abs(mine['efficiency'] - ref['efficiency']) < 1e-9, (name, m)
            done += 1
    assert done >= 4


@pytest.mark.skipif(not os.path.exists(stored_raw_path()), reason='no controls_raw.npz')
def test_stored_rows_are_consistent_with_their_operating_points_and_the_frozen_references():
    raw = X.load_raw(stored_raw_path())
    n = 0
    for name in X.SETS:
        if 'ref|%s|ops' % name not in raw:
            continue
        ops = X.jload(raw['ref|%s|ops' % name])
        a, b = (int(v) for v in raw['ref|%s|frames' % name][2:4])
        d0, d1 = (int(v) for v in raw['ref|%s|frames' % name][0:2])
        gt_dev = raw['ref|%s|gt_dev' % name]
        assert gt_dev['frame'].min() >= d0 and gt_dev['frame'].max() < d1 and raw['ref|%s|gt_test' % name]['frame'].min() >= a and raw['ref|%s|gt_test' % name]['frame'].max() < b
        for w in ('E', 'J'):
            o = ops[w]
            ds = X.dev_dets_of(raw, name, 'ISS-' + w)
            r = X.score_set(name, ds, gt_dev, np.arange(d0, d1))
            assert abs(r['efficiency'] - o['m'][4]) < 1e-6 and abs(r['jaccard'] - o['m'][0]) < 1e-9 and abs(r['rmse_nm'] - o['m'][3]) < 1e-6, (name, w)     # the development score of the operating point
            n += 1
        for m in ('GR-E', 'GR-J', 'GR-0'):
            ds = X.dets_of(raw, name, m)
            if ds is not None:
                assert len(ds) == b - a
        for m in ('R2Eo', 'R2Jo'):
            if X.dets_of(raw, name, m) is not None:
                dd = X.dev_dets_of(raw, name, m)
                r = X.score_set(name, dd, gt_dev, np.arange(d0, d1))
                s = X.jload(raw['ps|%s|summary' % name])['info'][m]
                res = X.jload(raw['ps|%s|search' % name])[s['key']][s['objective']]['metrics']
                key = 'efficiency' if s['objective'] == 'E' else 'jaccard'
                assert abs(r[key] - res[key]) < 1e-3, (name, m)                                       # the stored corrected development detections score what the search recorded
                n += 1
    assert n >= 8
