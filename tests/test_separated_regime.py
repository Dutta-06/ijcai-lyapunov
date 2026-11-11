"""The separated-regime experiment: the decoupled prediction, the target generators, the path on one target and the report."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_separated_regime as S      # noqa: E402


def test_gamma_is_the_overlap_of_two_kernels():
    assert S.gamma_of(0.0, 0.12) == 1.0
    assert abs(S.gamma_of(2 * 0.12, 0.12) - np.exp(-1.0)) < 1e-15           # d = 2 sigma: exp(-d^2 / (4 sigma^2)) = exp(-1)
    assert abs(S.gamma_of(6 * 0.12, 0.12) - np.exp(-9.0)) < 1e-15


def test_decoupled_picks_the_best_kernel_with_its_margin_and_hitting_time():
    xs = np.linspace(-1, 1, 41)
    c0 = np.array([-0.5, -0.3, 0.0, 0.3, 0.5])
    sig = 0.2
    Phi0 = np.exp(-(xs[:, None] - c0[None, :]) ** 2 / (2 * sig ** 2))
    bump = lambda c: np.exp(-(xs - c) ** 2 / (2 * sig ** 2))
    dec = S.decoupled(Phi0, bump, [0.28, -0.52], [2.0, -0.5])
    assert [d[0] for d in dec] == [3, 0]
    N = len(xs)
    a = np.abs(Phi0.T @ (2.0 * bump(0.28))) / N
    j2 = np.argsort(-a)[1]
    assert abs(dec[0][1] - 1.0 / a[3]) < 1e-12                               # T = 1 / |a_best|
    assert abs(dec[0][2] - (1.0 - a[j2] / a[3])) < 1e-12                     # margin against the second-best kernel
    assert dec[0][2] > 0
    order, Ts, mu0, g = S.decoupled_summary(dec)
    assert order == [j for j, T, m in sorted(dec, key=lambda x: x[1])]
    assert list(Ts) == sorted(Ts)
    assert mu0 == min(d[2] for d in dec)
    assert abs(g - (Ts[1] - Ts[0]) / Ts[1]) < 1e-12                          # relative gap between the two hitting times
    assert S.decoupled_summary(dec[:1])[3] == 1.0                            # one feature: no gap to speak of


def test_summary_gap_is_the_smallest_difference_of_three():
    dec = [(4, 1.0, 0.3), (9, 1.2, 0.5), (2, 2.0, 0.4)]
    order, Ts, mu0, g = S.decoupled_summary(dec)
    assert order == [4, 9, 2] and list(Ts) == [1.0, 1.2, 2.0]
    assert mu0 == 0.3
    assert abs(g - 0.2 / 2.0) < 1e-12                                        # min diff 0.2, not the largest (0.8), relative to the largest hitting time


def _fake_path(act_idx, events, rms, cs, ws):
    """A fake sliding_iss_nd result on the 21 starting kernels: `act_idx` active, `events` = (tau, kernel), centres `cs` and coefficients `ws` of the active ones."""
    n = 21
    act = np.zeros(n, bool)
    act[list(act_idx)] = True
    c = np.zeros((n, 1))
    w = np.zeros(n)
    c[list(act_idx), 0] = cs
    w[list(act_idx)] = ws
    return dict(events=[(t, j, 'enter', c[j]) for t, j in events], act=act, hist=[(0.0, c, w, rms)], c=c, w=w)


def test_one_scores_a_fake_path_field_by_field(monkeypatch):
    import exp_transit_certificate as XC
    L, K, c0, Phi0 = S.setup()
    import exp_mirror_transit as X
    cs, As = np.array([0.35, -0.45]), np.array([0.9, -0.5])                  # the feature at +0.35 (larger amplitude) enters first and has the larger kernel index
    dec = S.decoupled(Phi0, lambda c: X.bumps_target(K, [c], [1.0]), cs, As)
    order, Ts, mu0, g = S.decoupled_summary(dec)
    assert order[0] > order[1]                                               # descending kernel indices: an unsorted comparison with the active set would fail
    act_idx = sorted(order)

    def run(events, rms, centres=None, coefs=None):
        c_act = np.array([cs[1], cs[0]]) + (0.0 if centres is None else np.array(centres))     # active kernels ascending in index = ascending in centre
        w_act = np.array([As[1], As[0]]) + (0.0 if coefs is None else np.array(coefs))
        monkeypatch.setattr(XC, 'path_1d', lambda L_, K_, f, c0_, **kw: _fake_path(act_idx, events, rms, c_act, w_act))
        return dict(zip(S.COLS, S.one(cs, As, K.sig)))

    ev = [(Ts[0] * 1.01, order[0]), (Ts[1] * 0.98, order[1])]
    r = run(ev, 1e-12, centres=[2e-3, 0.0], coefs=[0.0, -3e-3])
    assert r['exact'] == 1.0 and r['units'] == 2.0 and r['exits'] == 0.0
    assert r['same_set'] == 1.0 and r['same_order'] == 1.0
    assert abs(r['terr'] - 0.02) < 1e-12                                     # the largest relative error of the two entry times
    assert abs(r['cerr'] - 2e-3) < 1e-12 and abs(r['aerr'] - 3e-3) < 1e-12  # compared with the truth sorted by centre
    assert r['mu0'] == mu0 and r['gap'] == g
    assert abs(r['sepmin'] - 0.8 / K.sig) < 1e-12                          # the actual minimal separation in widths
    assert run(ev, 1e-3)['exact'] == 0.0                                     # a residual above the tolerance is not exact
    swapped = [(Ts[0] * 1.01, order[1]), (Ts[1] * 0.98, order[0])]
    r2 = run(swapped, 1e-12)
    assert r2['same_set'] == 1.0 and r2['same_order'] == 0.0                 # the right kernels in the wrong order
    assert np.isnan(r2['terr'])                                              # the entry-time error is defined only where the entries are the decoupled prediction's
    assert np.isnan(run([(Ts[0] * 1.01, order[0])], 1e-12)['terr'])           # (one entry of two: not the prediction)


def test_one_three_features_sepmin_is_the_smallest_pair_distance(monkeypatch):
    import exp_transit_certificate as XC
    import exp_mirror_transit as X
    L, K, c0, Phi0 = S.setup()
    cs, As = np.array([-0.5, 0.0, 0.6]), np.array([0.9, -0.6, 0.7])
    dec = S.decoupled(Phi0, lambda c: X.bumps_target(K, [c], [1.0]), cs, As)
    order, Ts, mu0, gap = S.decoupled_summary(dec)
    act_idx = sorted(order)
    monkeypatch.setattr(XC, 'path_1d', lambda L_, K_, f, c0_, **kw: _fake_path(act_idx, [(T, j) for T, j in zip(Ts, order)], 1e-12, np.sort(cs), As[np.argsort(cs)]))
    r = dict(zip(S.COLS, S.one(cs, As, K.sig)))
    assert abs(r['sepmin'] - 0.5 / K.sig) < 1e-12                            # the pairs are 0.5 and 0.6 apart: the smaller, in widths
    assert r['same_set'] == 1.0 and r['same_order'] == 1.0 and r['terr'] < 1e-12


def test_target_generators():
    sig = 0.12
    cs, As = S.target_a(7, 3.5, sig)
    assert abs(cs[1] - cs[0] - 3.5 * sig) < 1e-12 and -0.7 <= cs[0] and cs[1] <= 0.7
    assert (np.abs(As) >= 0.4).all() and (np.abs(As) <= 1.0).all()
    assert S.target_a(7, 3.5, sig)[0].tolist() == cs.tolist()                # deterministic
    cb, Ab = S.target_b(11, 4.0, sig)
    assert len(cb) == 3 and min(abs(cb[i] - cb[j]) for i in range(3) for j in range(i)) >= 4.0 * sig
    assert (np.abs(Ab) >= 0.4).all()
    for seed in range(60):                                                   # the pair stays inside the domain and both amplitude signs occur
        c2, a2 = S.target_a(seed, 8.0, sig)
        assert -0.7 <= c2[0] and c2[1] <= 0.7 + 1e-12
    signs = np.concatenate([np.sign(S.target_b(seed, 3.0, sig)[1]) for seed in range(30)])
    assert (signs > 0).any() and (signs < 0).any()


def test_slope_recovers_a_power_law_and_ignores_nonpositive():
    x = np.array([1e-2, 1e-3, 1e-4, 1e-6])
    assert abs(S.slope(x, 3.0 * x ** 1.5) - 1.5) < 1e-9
    assert abs(S.slope(np.r_[x, 0.0], np.r_[3.0 * x, 1.0]) - 1.0) < 1e-9
    assert np.isnan(S.slope(np.array([1.0]), np.array([1.0])))


def test_one_target_at_six_sigma_is_the_decoupled_prediction():
    sig = 0.12
    cs, As = S.target_a(3, 6.0, sig)
    r = dict(zip(S.COLS, S.one(cs, As, sig)))
    assert r['exact'] == 1.0 and r['units'] == 2.0 and r['exits'] == 0.0
    assert r['same_set'] == 1.0 and r['same_order'] == 1.0
    h = 0.1
    assert r['terr'] < 3 * S.gamma_of(6.0 * sig - h / 2, sig)                # entry times within the first-order cross-talk of a kernel half a cell from its feature
    assert r['cerr'] < 1e-8 and r['aerr'] < 1e-8                             # the final block is the target
    assert r['mu0'] > 0 and r['gap'] > 0 and abs(r['sepmin'] - 6.0) < 1e-9


def _gam(sep_sig, sig=0.12, h=0.1):
    return float(np.exp(-max(sep_sig * sig - h / 2, 0.0) ** 2 / (4 * sig ** 2)))


def _synthetic_raw():
    """Six rows per separation with known counts: exact on 5, one exit, entry set right on 3, order right on 3 (rows 0-2; only they have an entry-time error, 0.5, 1 and 3 times
    Gamma(Delta - h/2) from 4 sigma, a constant 3 times those below), and margins such that the gated rows (min(mu0, gap) > Gamma(Delta - h/2)) differ from those given by mu0 alone;
    in part B the actual separation is 0.5 sigma above the row's lower bound."""
    raw = {'sig': np.array([0.12]), 'h': np.array([0.1])}
    for part, seps in (('A', S.SEPS_A), ('B', S.SEPS_B)):
        for s in seps:
            sep = s + (0.5 if part == 'B' else 0.0)
            G = _gam(sep)
            d = np.zeros((6, len(S.COLS)))
            d[:, 0] = [1, 1, 1, 1, 1, 0]
            d[:, 1] = [2, 2, 2, 2, 2, 3] if part == 'A' else [3, 3, 3, 3, 3, 2]
            d[:, 2] = [1, 0, 0, 0, 0, 0]
            d[:, 3] = [1, 1, 1, 0, 0, 0]
            d[:, 4] = [1, 1, 1, 0, 0, 0]
            d[:, 5] = [np.nan] * 6
            d[:3, 5] = (G if s >= 4.0 else 3.0) * np.array([0.5, 1.0, 6.0 if (s >= 6.0 and part == 'A') else 3.0])      # the largest error is 3 Gamma, 6 Gamma from 6 sigma (so that its slope is not 1)
            d[:, 6] = np.array([1, 2, 3, 4, 5, 6]) * 1e-15
            d[:, 7] = np.array([5, 4, 3, 2, 1, 0.5]) * 1e-15
            d[:, 8] = [0.9, 0.9, 0.9, 0.001 if part == 'A' else 0.003, 0.9, 0.9]
            d[:, 9] = [0.9, 0.0005, 0.9, 0.9, 0.9, 0.9]
            d[:, 10] = sep
            raw['%s|%g' % (part, s)] = d
    return raw


def test_report_rows_have_the_counts_of_the_raw(tmp_path, monkeypatch):
    raw = _synthetic_raw()
    md = S.report(raw)
    assert md == S.report(raw)
    G5 = _gam(5.0)
    # s = 5: Gamma(4.58 sigma) = 5.2e-3; margins min(mu0, gap) = 0.9, 0.0005, 0.9, 0.001, 0.9, 0.9: rows 0, 2, 4, 5 pass the gate (same set on rows 0, 2: 2 / 4); rows 1, 3 do not (row 1 has the right set: 1 / 2)
    row = '| 5 | %.2e | 5 / 6 | 1 | 3 / 6 | 3 / 6 | 2 / 4 | 1 / 2 | 3 | 1.00 | 3.00 | 6.0e-15 | 5.0e-15 | 4 / 6 |' % G5
    assert row in md, row
    G2 = _gam(2.0)                                                           # Gamma(1.58 sigma) = 0.53 < 0.9 but > 0.001: rows 0, 2, 4, 5 pass the gate again
    assert '| 2 | %.2e | 5 / 6 | 1 | 3 / 6 | 3 / 6 | 2 / 4 | 1 / 2 | 3 | %.2f | %.2f | 6.0e-15 | 5.0e-15 | 4 / 6 |' % (G2, 3.0 / G2, 9.0 / G2) in md       # below 4 sigma the errors are constant, not proportional to Gamma
    G8 = _gam(8.0)                                                           # Gamma(7.58 sigma) = 5.6e-7 < every margin: all six rows pass the gate, 3 of 6 have the right set
    assert '| 8 | %.2e | 5 / 6 | 1 | 3 / 6 | 3 / 6 | 3 / 6 | 0 / 0 | 3 | 1.00 | 6.00 | 6.0e-15 | 5.0e-15 | 4 / 6 |' % G8 in md
    # part B: the gate uses the actual separation (row + 0.5 sigma): row 3 (mu0 = 0.003) passes at 5.5 sigma (Gamma = 1.6e-3) and would fail at the nominal 5 sigma (5.2e-3)
    partb = md.split('## K = 3')[1]
    assert '| 5 | %.2e | 5 / 6 | 1 | 3 / 6 | 3 / 6 | 2 / 5 | 1 / 1 | 3 | 1.00 | 3.00 | 6.0e-15 | 5.0e-15 | 4 / 6 |' % _gam(5.0) in partb
    seps_big = [s for s in S.SEPS_A if s >= 4.0]                             # the slope of the largest error (3 Gamma, then 6 Gamma) against Gamma, part A only
    exp_slope = np.polyfit(np.log([_gam(s) for s in seps_big]), np.log([(6.0 if s >= 6.0 else 3.0) * _gam(s) for s in seps_big]), 1)[0]
    assert abs(exp_slope - 1.0) > 0.05
    assert md.count('log-log slope of the largest entry-time error') == 1
    assert 'over the separations of at least 4 sigma: %.2f (first-order' % exp_slope in md
    # totals, recomputed here from the definitions
    n_gate = n_ok = n_non = n_non_ok = 0
    for part, seps in (('A', S.SEPS_A), ('B', S.SEPS_B)):
        for s in seps:
            d = raw['%s|%g' % (part, s)]
            gate = np.minimum(d[:, 8], d[:, 9]) > _gam(d[0, 10])
            ss = d[:, 3] > 0
            n_gate += gate.sum(); n_ok += (gate & ss).sum(); n_non += (~gate).sum(); n_non_ok += (~gate & ss).sum()
    assert ('(min(mu0, gap) > Gamma(Delta - h/2): %d targets) the path\'s entry set is the decoupled prediction\'s on %d; elsewhere (%d targets) on %d.' % (n_gate, n_ok, n_non, n_non_ok)) in md
    assert 0 < n_gate < 72 and 0 < n_ok < n_gate
    assert 'is 6.00.' in md                                                  # the largest ratio over the rows at 5 sigma and above
    for needle in ('K = 2, separation exactly', 'K = 3, minimal separation', 'exact with exactly K units and no exit on 60', 'Over all 72 targets'):
        assert needle in md, needle
    assert 'nan' not in md.lower().replace('nan,', '')
    monkeypatch.setattr(S, 'OUT', str(tmp_path))
    monkeypatch.setattr(S, 'RAW', str(tmp_path / 'r.npz'))
    np.savez_compressed(S.RAW, **raw)
    assert all(np.array_equal(S.load_raw()[k], raw[k], equal_nan=True) for k in raw)


# ---- part C: K = 2 at 10, 12 and 14 sigma in a wider window (added after the full-draft review) ------------------------------------------

def test_target_c_pairs_are_exactly_sep_sigma_apart_inside_the_wider_window():
    sig = 0.12
    for seed in range(40):
        for sep in S.SEPS_C:
            cs, As = S.target_c(seed, sep, sig)
            assert abs(cs[1] - cs[0] - sep * sig) < 1e-12 and -S.WIN_C <= cs[0] and cs[1] <= S.WIN_C + 1e-12
            assert (np.abs(As) >= 0.4).all() and (np.abs(As) <= 1.0).all()
    assert S.target_c(5, 12.0, sig)[0].tolist() == S.target_c(5, 12.0, sig)[0].tolist()
    assert S.target_c(5, 12.0, sig)[0].tolist() != S.target_c(6, 12.0, sig)[0].tolist()
    assert max(S.target_c(seed, 14.0, sig)[0][0] for seed in range(200)) <= S.WIN_C - 14.0 * sig + 1e-12
    assert max(S.SEPS_C) * sig <= 2 * S.WIN_C                                # the widest pair fits in the window; part A's window (1.4) would not hold it


def test_job_draws_part_c_from_its_own_seed_and_generator(monkeypatch):
    seen = []
    monkeypatch.setattr(S, 'setup', lambda: (None, type('K', (), dict(sig=0.12))(), None, None))
    monkeypatch.setattr(S, 'one', lambda cs, As, sig: seen.append((list(cs), list(As), sig)) or np.zeros(len(S.COLS)))
    S.job(('C', 12.0, 7))
    S.job(('A', 5.0, 7))
    S.job(('B', 4.0, 7))
    cs_c, As_c = S.target_c(S.SEED_C + 7, 12.0, 0.12)
    assert seen[0] == (list(cs_c), list(As_c), 0.12)
    assert seen[1][0] == list(S.target_a(S.SEED_A + 7, 5.0, 0.12)[0]) and seen[2][0] == list(S.target_b(S.SEED_B + 7, 4.0, 0.12)[0])
    assert len({S.SEED_A, S.SEED_B, S.SEED_C}) == 3


def test_run_extra_adds_part_c_and_leaves_the_keys_of_a_and_b_alone(tmp_path, monkeypatch):
    raw = _synthetic_raw()
    monkeypatch.setattr(S, 'OUT', str(tmp_path))
    monkeypatch.setattr(S, 'RAW', str(tmp_path / 'r.npz'))
    np.savez_compressed(S.RAW, **raw)
    calls = []

    def fake_job(args):
        calls.append(args)
        return np.full(len(S.COLS), float(args[1]))
    monkeypatch.setattr(S, 'job', fake_job)
    out = S.run_extra(1)
    assert sorted(a for a in calls) == sorted(('C', s, i) for s in S.SEPS_C for i in range(S.N_PER))
    again = S.load_raw()
    for k in raw:
        assert np.array_equal(again[k], raw[k], equal_nan=True), k
    for s in S.SEPS_C:
        assert again['C|%g' % s].shape == (S.N_PER, len(S.COLS)) and (again['C|%g' % s] == s).all()
    assert set(out) == set(again)


def _raw_with_c():
    raw = _synthetic_raw()
    for s in S.SEPS_C:
        d = np.zeros((6, len(S.COLS)))
        d[:, 0] = 1.0
        d[:, 1] = 2.0
        d[:, 3] = 1.0
        d[:, 4] = 1.0
        d[:, 5] = np.array([1e-11, 2e-11, 3e-11, 4e-11, 5e-11, 6e-11]) * (s / 10.0)
        d[:, 6] = 1e-13
        d[:, 7] = 2e-13
        d[:, 8] = 0.9
        d[:, 9] = 0.9
        d[:, 10] = s
        raw['C|%g' % s] = d
    return raw


def test_report_with_part_c_has_its_table_its_note_and_unchanged_other_numbers():
    base, withc = S.report(_synthetic_raw()), S.report(_raw_with_c())
    assert '## K = 2, separation exactly Delta, wider window' not in base and 'solver' not in base
    assert '## K = 2, separation exactly Delta, wider window' in withc and 'Over all 90 targets of A, B and C' in withc and 'Over all 72 targets of A and B' in base
    for s in S.SEPS_C:                                                       # exact on 6 of 6, no exit, same set 6 of 6, the ratio columns are not given
        assert '| %g | %.2e | 6 / 6 | 0 | 6 / 6 | 6 / 6 | 6 / 6 | 0 / 0 | 6 | - | - | 1.0e-13 | 2.0e-13 | 6 / 6 |' % (s, _gam(s)) in withc, s
    assert "largest relative entry-time error over the targets with the decoupled entries is 8.4e-11" in withc            # 6e-11 x 14 / 10
    assert 'is 6.00.' in withc and 'is 6.00.' in base                      # the largest ratio at 5 sigma and above does not read part C
    assert 'Over all 90 targets of A, B and C the path is exact with exactly K units and no exit on 78.' in withc      # 9 x 5 + 3 x 5 of A and B, 3 x 6 of C
    assert 'Over all 72 targets of A and B the path is exact with exactly K units and no exit on 60.' in base
    assert 'log-log slope' in withc and withc.count('log-log slope') == 1
    # the lines of A and B are the same text with and without C
    a_b = lambda md: md.split('## K = 2, separation exactly Delta, wider window')[0]
    assert a_b(withc).split('Over all')[0] == base.split('Over all')[0]
