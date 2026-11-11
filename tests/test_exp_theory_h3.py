"""The coupled-entry-layer experiment: the box conditions, the inner problem, the certificate and the report (deterministic, fast)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_theory_h3 as T      # noqa: E402


def test_gamma_psi_are_the_overlap_and_its_slope():
    assert abs(float(T.gam(2.0)) - np.exp(-1.0)) < 1e-15                      # Gamma(d) = exp(-d^2 / 4)
    assert abs(float(T.psi(2.0)) - np.exp(-1.0)) < 1e-15                      # psi(d) = (d/2) Gamma(d)
    x = np.linspace(np.sqrt(2.0), 8.0, 400)
    assert np.all(np.diff(T.psi(x)) < 0)                                      # psi is decreasing beyond sqrt 2
    y = np.linspace(np.sqrt(6.0), 10.0, 400)
    assert np.all(np.diff(T.g2(y)) < 0) and np.all(T.g2(y) > 0)               # |Gamma''| decreasing beyond sqrt 6
    d = np.linspace(0.0, np.sqrt(2.0), 400)
    assert np.all(np.diff(T.psi(d)) > 0)                                      # psi increasing on [0, sqrt 2]: the inverse used for d_*


def test_box_holds_far_and_fails_near_and_is_monotone_in_the_separation():
    assert not T.box(2, 1.0, 3.0, 0.41, use_arb=False)['ok']
    assert T.box(2, 1.0, 6.0, 0.41, use_arb=False)['ok']
    ok = [T.box_exists(2, 1.0, D, 0.41, n_r=24, use_arb=False) is not None for D in np.arange(4.0, 9.0, 0.1)]
    first = ok.index(True)
    assert all(ok[first:]) and not any(ok[:first])                            # once it holds it holds for every larger separation
    thr = T.entry_threshold(2, 1.0, 0.41)
    assert 4.5 < thr < 4.8                                                    # the table's K = 2, R = 1, d_0 = 0.41 entry
    assert T.entry_threshold(2, 2.4, 0.41) > thr                              # a larger amplitude ratio needs more room
    assert T.entry_threshold(3, 1.0, 0.41) > thr                              # and so do more features
    assert T.entry_threshold(2, 1.0, 0.8) > thr                               # and a larger starting offset


def test_arithmetic_agrees_with_floats_and_dormant_mass_shrinks_the_box():
    a, b = T.box(2, 2.4, 6.5, 0.41, use_arb=False), T.box(2, 2.4, 6.5, 0.41, use_arb=T.HAVE_ARB)
    assert a['ok'] and b['ok']
    for k in ('Ev', 'Ec', 'vbar', 'vlo', 'Cmargin', 'dstar', 'lam', 'Hc'):
        assert abs(a[k] - b[k]) < 1e-9 * max(1.0, abs(a[k]))
    assert T.box(2, 1.0, 6.0, 0.41, use_arb=False, Omega=0.0)['ok']
    assert not T.box(2, 1.0, 6.0, 0.41, use_arb=False, Omega=0.5)['ok']
    assert T.box(2, 1.0, 6.0, 0.41, use_arb=False, Omega=0.0)['Nmargin'] > T.box(2, 1.0, 6.0, 0.41, use_arb=False, Omega=0.01)['Nmargin']


def test_box_conditions_are_the_faces_signs_of_the_vector_field():
    """On the faces of a box that satisfies (L), (C), the field of the two-unit fit points inward, for the worst couplings, on a grid of the box."""
    D, R, r = 6.0, 1.7, 0.45
    b = T.box(2, R, D, r, use_arb=False)
    assert b['ok']
    cs, A = np.array([0.0, D]) * T.SIG, np.array([1.0, R])
    rng = np.random.default_rng(3)
    for _ in range(4000):
        d = rng.uniform(-r, r, 2)
        v = np.array([rng.uniform(b['vlo'], b['vbar']), rng.uniform(b['vlo'], b['vbar'])])
        c = cs + d * T.SIG
        g, h = T._gh(c, v * np.sign(A), cs, A)
        gv = g / T.Z_CONT                                                    # d f / d v_i, negated, for positive coefficients
        # dv/ds ~ +gv: inward on the lower face and outward-forbidden on the upper face
        assert gv[0] > 0 or v[0] > b['vlo'] + 1e-9                           # unit 1 at its lower face rises
        assert gv[1] > 0 or v[1] > b['vlo'] + 1e-9
    for sgn in (-1.0, 1.0):                                                   # the centres' faces
        for _ in range(500):
            d = np.array([sgn * r, rng.uniform(-r, r)])
            v = np.array([rng.uniform(b['vlo'], b['vbar']), rng.uniform(b['vlo'], b['vbar'])])
            c = cs + d * T.SIG
            g, h = T._gh(c, v, cs, A)
            assert np.sign(v[0] * h[0]) == -sgn                              # the centre moves inward


def test_inner_problem_reaches_the_target_inside_the_box():
    run = T.rip_run(6.0, 0.41, (0.7, 0.9), 1e-2, t_chunk=3000.0, max_chunks=30)
    assert run['still']
    assert np.allclose(run['w'], [0.7, 0.9], atol=1e-7) and np.allclose(run['c'], [0.0, 6.0], atol=1e-7)
    r = T.box_exists(2, 0.9 / 0.7, 6.0, 0.41)
    assert r is not None and run['dmax'].max() <= r


def test_certificate_on_a_target_reproduces_the_path(monkeypatch):
    S, L, K, c0, Phi0 = T.setup()
    cs, As = S.target_a(1000 + 3, 6, K.sig)
    o = np.argsort(cs)
    cs, As = cs[o], As[o]
    ce = T.certify(cs, As, c0, T.Inner('samples'))
    assert ce['dual_ok'] and ce['box_ok'] and ce['certified']
    ev, nex, rms, act = T.path_of(cs, As)
    assert [j for t, j in ev] == ce['entries'] and nex == 0 and rms < 1e-6
    assert max(abs(t - tt) / tt for (t, j), tt in zip(ev, ce['times'])) < 1e-3
    assert abs(ce['sep'] - 6.0) < 1e-9
    # at 2 sigma the dual part still predicts the path; the analytic box does not hold
    cs2, As2 = S.target_a(1000 + 3, 2, K.sig)
    o = np.argsort(cs2)
    ce2 = T.certify(cs2[o], As2[o], c0, T.Inner('samples'))
    ev2 = T.path_of(cs2[o], As2[o])[0]
    assert ce2['dual_ok'] and not ce2['box_ok'] and [j for t, j in ev2] == ce2['entries']


def test_certificate_rejects_a_sign_mismatch_and_a_tie():
    S, L, K, c0, Phi0 = T.setup()
    inner = T.Inner('continuum')
    # two features exactly between kernels: a tie of hitting times between features
    h = 0.1
    cs = np.array([-0.45, 0.45])
    ce = T.certify(cs, np.array([1.0, 1.0]), c0, inner)
    assert not ce['dual_ok']                                                 # equal offsets, equal amplitudes: the two entries tie (gap 0)
    # opposite signs close together: the first entrant's dual sign is the neighbour's
    ce = T.certify(np.array([0.0, 0.5 * 0.12]), np.array([1.0, -3.0]), c0, inner)
    assert not ce['certified']


def test_report_is_regenerated_byte_for_byte_from_the_raw(tmp_path, monkeypatch):
    monkeypatch.setattr(T, 'OUT', str(tmp_path))
    monkeypatch.setattr(T, 'RAW', str(tmp_path / 'theory_h3_raw.npz'))
    monkeypatch.setattr(T, 'MD', str(tmp_path / 'theory_h3.md'))
    raw = {'have_arb': np.array([float(T.HAVE_ARB)]), 'A_thr': np.full((2, 5, 5), 4.5), 'A_om': np.full((3, 5), 0.01)}
    T.save_raw(raw)
    a = T.report(T.load_raw())
    b = T.report(T.load_raw())
    assert a == b and a.endswith('\n') and 'K = 2' in a and 'K = 3' in a and '4.50' in a
    monkeypatch.setattr(sys, 'argv', ['exp_theory_h3.py', '--from-raw'])
    T.main()
    first = (tmp_path / 'theory_h3.md').read_bytes()
    T.main()
    assert first == (tmp_path / 'theory_h3.md').read_bytes() == a.encode()


def test_invariance_needs_only_L_and_C_and_tolerates_more_dormant_mass():
    # (L), (C) alone tolerate a much larger Omega than the full set with (N)
    n_, lc = T.omega_max(2, 1.0, 6.0, T.D0_BENCH, True), T.omega_max(2, 1.0, 6.0, T.D0_BENCH, False)
    assert 0.2 < n_ < 0.3 and 0.6 < lc < 0.8 and lc > 2 * n_
    assert T.box_exists(2, 1.0, 6.0, T.D0_BENCH, Omega=T.OMEGA_WIN) is not None
    assert abs(T.D0_BENCH - 0.05 / 0.12) < 1e-3


def test_window_departure_decays_with_the_distance_to_the_edge():
    W = T.part_W()['W_dev']
    assert W.shape == (len(T.W_EDGE), 4)
    assert np.all(np.diff(W[:, 0]) < 0) and W[0, 0] > 0.1 and W[T.W_EDGE.index(3.5), 0] < 1e-4      # <phi^2>: large at 1 width, negligible at 3.5
    assert W[T.W_EDGE.index(3.5)].max() < 1e-4 and W[0, 1] > 0.1


def test_interior_flag_and_the_budgeted_box():
    S, L, K, c0, Phi0 = T.setup()
    cs, As = S.target_a(1000 + 3, 6, K.sig)                                   # centres in [-0.7, 0.7]; here -0.34 and 0.37: interior (|c| <= 0.58)
    o = np.argsort(cs)
    ce = T.certify(cs[o], As[o], c0, T.Inner('samples'))
    assert ce['interior'] and ce['certified_int'] and ce['box_ok_rob'] and ce['certified']
    cs2 = np.array([0.7, -0.1])                                               # a feature 2.5 widths from the edge: not interior
    ce2 = T.certify(np.sort(cs2), np.array([0.8, 0.9]), c0, T.Inner('samples'))
    assert not ce2['interior'] and not ce2['certified_int']
    ce3 = T.certify(np.sort(cs2), np.array([0.8, 0.9]), c0, T.Inner('continuum'))
    assert ce3['interior']                                                    # in the continuum there is no edge


def test_unpolished_path_agrees_with_the_certificate_on_a_target():
    import exp_transit_certificate as XC
    import exp_mirror_transit as X
    S, L, K, c0, Phi0 = T.setup()
    cs, As = S.target_a(1000 + 3, 6, K.sig)
    o = np.argsort(cs)
    cs, As = cs[o], As[o]
    R = XC.path_1d(L, K, X.bumps_target(K, cs, As), c0, polish=False)
    ev = [int(j) for t, j, k, c in R['events'] if k == 'enter']
    ce = T.certify(cs, As, c0, T.Inner('samples'))
    assert ev == ce['entries'] and R['hist'][-1][3] < 1e-8 and all(R['relaxed'])


def test_wide_window_targets_are_interior_deterministic_and_distinct():
    for K, sep, R, mode in T.i_cells()[::3]:
        cs, As = T.i_target(K, sep, R, mode, 123)
        cs2, As2 = T.i_target(K, sep, R, mode, 123)
        assert np.array_equal(cs, cs2) and np.array_equal(As, As2)
        assert np.abs(cs).max() <= T.I_HALF - T.WIN_EDGE * T.SIG + 1e-9                  # interior of the wider window
        assert len(cs) == K and np.all(np.diff(cs) > 0) and abs(abs(As).max() / abs(As).min() - R) < 1e-9
        if mode == 'worst':                                                              # every feature 0.9-1.1 x 0.045 from a lattice point (never exactly equal: no tie)
            off = np.abs(cs - np.round(cs / 0.1) * 0.1)
            assert np.all((off > 0.0405 - 1e-9) & (off < 0.0495 + 1e-9)) and len(set(np.round(off, 9))) == K


def test_part_I_runs_both_paths_on_an_interior_certified_target(monkeypatch):
    monkeypatch.setattr(T, 'i_cells', lambda: [(3, 7.0, 2.4, 'random')])
    monkeypatch.setattr(T, 'I_N', 1)
    r = T.part_I()['I_rows']
    c = {n: i for i, n in enumerate(T.I_COLS)}
    assert r.shape == (1, len(T.I_COLS)) and r[0, c['interior']] == 1 and r[0, c['certified_int']] == 1
    assert r[0, c['ent_T']] == r[0, c['ent_F']] == r[0, c['exact_T']] == r[0, c['exact_F']] == 1
    assert r[0, c['diff_TF']] < 1e-8 and r[0, c['cerr_F']] < 1e-8
