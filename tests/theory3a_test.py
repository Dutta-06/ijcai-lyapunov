"""Theorem T1, the small-eps tail (experiments/theory3a_tail.py, theory3a_chain.py, theory3a_run.py): the h = 1/kappa parametrisation against the existing unit, the envelope of the dormant coefficient,
the tail verifier (accepts, encloses the point verifier, rejects when it should), the sound displacement accounting, the integrated network inside the tail, and the deterministic report.
Run explicitly:  pytest tests/theory3a_test.py   (about 3-4 min; the file name keeps it out of the default collection)."""
import os
import sys

import numpy as np
import pytest

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
from flint import arb                  # noqa: E402
import exp_theory_t1 as T              # noqa: E402
import theory_t1b_rigorous as RG       # noqa: E402
import theory_t1b_wide as W            # noqa: E402
import theory3a_tail as TL             # noqa: E402
import theory3a_chain as CH            # noqa: E402
import theory3a_run as RUN             # noqa: E402

RAW3 = os.path.join(os.environ.get('TH3A_OUT', os.path.join(REPO, 'results', 'synthetic')), 'theory3a_raw.npz')


def _target(n):
    K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
    prep = T.prepare(*T.make_target(K, sep, R, seed))
    return prep, [e['kernel'] for e in prep['L']['events']], len(prep['c0'])


def _cands(prep, ent, J, e, fac):
    o, pi = T.best_cert(prep, 10.0 ** e)
    assert pi >= 0
    return T.PARS[pi], [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(J)], o


def test_parametrisation_h_L_matches_the_existing_unit():
    for eps in (1e-4, 1e-9, 1e-14):
        U0, U1, F = RG.RU(arb(eps)), TL.UH(TL.h_of_eps(eps)), T.Unit(eps)
        assert abs(float(1 / U1.L.mid()) / F.lam - 1) < 1e-12 and abs(float(U1.eps.mid()) / eps - 1) < 1e-12      # lambda = 1/L, eps = 2 v_s h^2/(1 + h^2)
        for x in (0.3, 1.3, 2.0):
            xa = arb(x)
            assert abs(float(U1.s(xa).mid()) - float(U0.s(xa).mid())) < 1e-3 * max(1.0, abs(float(U0.s(xa).mid())))   # RU's series ball is wider; same value
            assert abs(float(U1.s(xa).mid()) - F.s_of_x(x)) < 1e-9 * max(1.0, abs(F.s_of_x(x)))
            assert abs(float(U1.w(xa).mid()) / float(U0.w(xa).mid()) - 1) < 1e-12
            assert abs(float(U1.M(xa).mid()) / float(U0.M(xa).mid()) - 1) < 1e-12
            assert abs(float(U1.Psi(xa).mid()) / float(U0.Psi(xa).mid()) - 1) < 1e-12
            assert float(U1.deta_dx(xa).upper()) >= float(U0.deta_dx(xa).lower()) * (1 - 1e-12)          # the upper-bound form of d eta/dx
            ratio = float(U1.deta_dx(xa).mid()) / float(U0.deta_dx(xa).mid())
            h2, t2 = float(U1.h2.mid()), float(np.tanh(x) ** 2)
            assert 1 - 1e-12 <= ratio <= 1 + h2 / t2 * 1.0001                                          # d eta/dx at eps = 0 against the exact value: ratio = 1 + h^2 (1 - t^2)/(t^2 (1 + h^2))
    # h is a ball [0, h_hi]: the functions enclose the point values for every h in it
    hb = arb(0).union(TL.h_of_eps(1e-9))
    U = TL.UH(hb)
    for e in (1e-9, 1e-12, 1e-18):
        F = T.Unit(e)
        sv, wv = float(F.s_of_x(1.5)), float(F.w_of_x(1.5))                      # the floating-point values carry roundoff ~1e-16 relative
        assert U.s(arb(1.5)).overlaps(arb(sv, 1e-13 * abs(sv))) and U.w(arb(1.5)).overlaps(arb(wv, 1e-13 * abs(wv)))
    # the inversions enclose the point inversions
    xs_, xw_ = float(T.Unit(1e-12).x_of_s(-15.0)), float(T.Unit(1e-12).x_of_w(0.1))
    assert U.x_of_s(-15).overlaps(arb(xs_, 1e-7)) and U.x_of_w(0.1).overlaps(arb(xw_, 1e-7))          # Unit's tables interpolate at ~1e-9


def test_g_atan_encloses_atan_over_y():
    for lo_, hi_ in [(0.0, 1e-6), (0.0, 0.3), (0.2, 0.4), (0.0, 0.7), (0.6, 1.2), (1e-9, 1e-8), (0.0, 1.5)]:
        g = TL.g_atan(arb(lo_).union(arb(hi_)))
        for y in np.linspace(max(lo_, 1e-12), hi_, 25):
            assert g.contains(arb(float(np.arctan(y) / y))), (lo_, hi_, y)
    assert TL.g_atan(arb(0)).contains(arb(1))
    for y in (0.1, 0.3, 0.45):                                                       # a point: only the series' error bar covers the truncation (7 terms: y^14/15)
        g = TL.g_atan(arb(y))
        assert g.overlaps(arb(float(np.arctan(y) / y), 1e-13)) and float(g.rad()) < 1e-5 and float(g.rad()) > 0.5 * y ** 14 / 15 - 1e-300


def test_nonnegative_squares_and_nan_guard():
    z = arb(0).union(arb(2e-8))
    assert float(TL.sqn(z).lower()) >= 0 and float(TL.sqn(z).upper()) >= 4e-16
    assert float(TL.sqn(arb(-3)).lower()) > 8.9                                           # (-3)^2
    with pytest.raises(ValueError):
        TL.bmax_up(arb('nan'), arb(0))
    assert float(TL.bmax_up(arb(1), arb(2)).upper()) == 2.0


def test_envelope_of_the_dormant_coefficient_holds_in_ball_arithmetic_and_against_floats():
    raw = RUN.part_U(quick=True)
    assert np.all(raw['U_rows'][:, 2] <= 1.0) and np.all(raw['U_rows'][:, 1] > 50)
    for eps in (1e-6, 1e-11):
        F = T.Unit(eps)
        s = np.geomspace(26.0, 0.999 * F.lam, 600)
        env = float(TL.C_ENV.mid()) / (s * (s * s - 400.0))
        assert (F.B(s) / env).max() <= 1 + 1e-6
    assert abs(float(TL.C_ENV.mid()) - 4 / (T.LAW['G'] * T.LAW['gamma'] ** 2 * T.LAW['m'])) < 1e-12        # C = 4/(G gamma^2 m)


def test_tail_verifier_accepts_and_encloses_the_point_verifier():
    prep, ent, J = _target(0)
    par, Dv, o = _cands(prep, ent, J, -10.0, 1.5)
    ok, info = TL.verify_tail(prep, par, Dv, TL.h_of_eps(1e-10))
    assert ok and info['worst_ratio'] < 1 and info['cond']['order'] and info['cond']['dc'] and info['cond']['theta']
    assert info['eps_hi'] >= 1e-10 * (1 - 1e-12) and info['eps_hi'] < 1.01e-10                  # the ball [0, h_hi] reaches eps = 1e-10 (h_hi is rounded up)
    for e_ in (1e-4, 1e-9):
        assert abs(float(TL.eps_of_h(TL.h_of_eps(e_)).mid()) / e_ - 1) < 1e-12 or abs(float(TL.eps_of_h(TL.h_of_eps(e_)).upper()) / e_ - 1) < 1e-12      # eps(h(eps)) = eps
        assert abs(float(TL.eps_of_h(arb(0.3)).mid()) - float(2 * (1 / 40) * 0.09 / 1.09)) < 1e-12                                                     # the 1 + h^2
    for e in (-10.0, -13.0):                                         # the tail ball contains every smaller eps: its numbers dominate the point verifier's (sound displacement accounting)
        ok2, i2 = CH.verify_fixed(prep, arb(10.0 ** e), par, Dv)
        assert ok2
        assert info['T1'] >= i2['T1'] * (1 - 1e-12)                                                    # an upper Riemann sum: only the upper end means anything
        assert info['prelo'] <= i2['pre'] * (1 + 1e-12) and info['pre'] >= i2['pre'] * (1 - 1e-12)
        assert info['Atr'] >= i2['Atr'] * (1 - 1e-12) and info['T2'] >= i2['T2'] * (1 - 1e-12) and info['T3'] >= i2['T3'] * (1 - 1e-12) and info['Tlay'] >= i2['Tlay'] * (1 - 1e-12)
        assert info['kappa'] <= i2['kappa'] * (1 + 1e-12)                                             # kappa enters as a divisor: the tail value is the smaller one
        assert info['De'][0] >= i2['De'][0] * (1 - 1e-9)
        assert abs(info['omega_T'] - info['omega'] * info['tau_end'] / info['Lhi']) < 1e-6 * info['omega_T']      # omega_pk x Tfin = omega_pk tau_end / L_hi
        # the dormant terms dominate the floating-point (exact B, mass identity) values of `cert` at the same eps; the ball verifier's own mass is conservative and not a lower bound of ours
        ofl, _ = T.best_cert(prep, 10.0 ** e)
        assert info['mass'] >= ofl['mass_ne'] * (1 - 1e-6) and info['omega'] >= ofl['omega_peak'] * (1 - 1e-6)
        if e == -10.0:                                                                                  # tightness at the top of the ball only (smaller eps: the tail value is the one at eps_hi)
            assert info['mass'] <= 1.2 * ofl['mass_ne'] and info['omega'] <= 5 * ofl['omega_peak']      # the envelope bounds are not loose (sum of peaks against the peak of the sum: a factor ~2)
            assert info['Vtail'] >= ofl['Vtail'] * (1 - 1e-3) and info['Vtail'] <= ofl['Vtail'] * 1.2
            assert info['De'][0] <= i2['De'][0] * 1.005                                                 # pre-layer drift is the same Riemann sum
        for j, v in ofl['Etot'].items() if isinstance(ofl['Etot'], dict) else enumerate(ofl['Etot']):
            if j in info['Etot']:
                assert info['Etot'][j] >= v * (1 - 3e-3)
        assert info['dcmax'] >= i2['dcmax'] * (1 - 1e-9)
        if e == -10.0:
            assert max(info['Etot'].values()) <= 1.3 * max(i2['Etot'].values()) + 1.0                   # the drift bound without caps and gating costs a few percent, not a factor
    # the enclosure widths shrink with eps_hi (the remainder is O(eps_hi))
    ok3, i3 = TL.verify_tail(prep, par, Dv, TL.h_of_eps(1e-14))
    assert ok3 and (i3['T1'] - i3['T1lo']) < 1e-3 * (info['T1'] - info['T1lo']) + 1e-9 and (i3['pre'] - i3['prelo']) < 1e-3 * (info['pre'] - info['prelo']) + 1e-9
    assert abs(i3['T1'] - 44.62521049) < 1e-6                                   # the eps -> 0 value of the layer integral (Part S of theory3a.md)


def test_tail_verifier_rejects_what_it_must_reject():
    prep, ent, J = _target(0)
    par, Dv, o = _cands(prep, ent, J, -10.0, 1.5)
    h = TL.h_of_eps(1e-10)
    non = [j for j in range(J) if j not in ent]
    small = [d if j in ent else 0.5 * o['Etot'][j] for j, d in enumerate(Dv)]
    ok, info = TL.verify_tail(prep, par, small, h)
    assert not ok and not info['cond']['fixed point E<=D']                      # D below the proved deviation
    big = [d if j in ent else 1e9 for j, d in enumerate(Dv)]
    ok, info = TL.verify_tail(prep, par, big, h)
    assert not ok and 'dist' in info['why']                                    # D larger than the margin lambda * mu
    ok, info = TL.verify_tail(prep, par, Dv, TL.h_of_eps(1e-3))
    assert not ok                                                              # eps_hi far above the threshold
    tau = 1.05 * prep['L']['events'][-1]['tau'] + 1.0
    okb, ib = TL.verify_tail(prep, par, Dv, h)
    assert okb
    j0 = min(non, key=lambda j: float(TL.kernel_margin(prep, j, tau).mid()))
    c0m = float(TL.kernel_margin(prep, j0, tau).mid())
    D_edge = list(Dv)
    D_edge[j0] = (c0m - 15.0 * ib['Lhi']) / ib['Lhi']                           # lambda mu - D = 15 < 26: the envelope (s >= 25) is not available
    ok, info = TL.verify_tail(prep, par, D_edge, h)
    assert not ok and 'dist' in info['why']
    D_edge[j0] = (c0m - 40.0 * ib['Lhi']) / ib['Lhi']                           # lambda mu - D = 40: accepted by the distance test (the deviation test may still be decided either way)
    ok, info = TL.verify_tail(prep, par, D_edge, h)
    assert 'dist' not in info.get('why', '')
    p2 = dict(par, dc=1e-6)
    ok, info = TL.verify_tail(prep, p2, Dv, h)
    assert not ok and not info['cond']['dc']                                   # the displacement bound exceeds a tiny d_c
    p3 = dict(par, Om=1e-12)
    ok, info = TL.verify_tail(prep, p3, Dv, h)
    assert not ok                                                              # dormant budget Omega too small (or box)
    assert len(non) > 10


def test_strict_bootstrap_factor_and_solve_D_iteration():
    prep, ent, J = _target(0)
    par, Dv, o = _cands(prep, ent, J, -10.0, 1.5)
    h = TL.h_of_eps(1e-10)
    ok, info = TL.verify_tail(prep, par, Dv, h)
    assert ok
    non = [j for j in range(J) if j not in ent]
    js = max(non, key=lambda j: info['Etot'][j])
    def fixed_point(fac):                                                          # D_js = fac * E_js(D) (E depends weakly on D through the dormant terms)
        D_ = list(Dv)
        for _ in range(12):
            _, inf_ = TL.verify_tail(prep, par, D_, h)
            D_[js] = inf_['Etot'][js] * fac
        return D_
    D1 = fixed_point(1.00005)                                                      # inside the strictness factor 1.0001: must be refused
    ok1, i1 = TL.verify_tail(prep, par, D1, h)
    assert not ok1 and not i1['cond']['fixed point E<=D'] and 1.00004 < D1[js] / i1['Etot'][js] < 1.00006
    D2 = fixed_point(1.0002)
    ok2, i2 = TL.verify_tail(prep, par, D2, h)
    assert i2['cond']['fixed point E<=D']
    # solve_D: the first starts are all too small (0.3 E); the iteration D <- 1.0003 E(D) repairs it in a step or two
    ok3, i3, D3 = TL.solve_D(prep, par, h, [0.3 * o['Etot'][j] if j not in ent else 0.0 for j in range(J)], ladder=(1.0, 1.1), iters=4)
    assert ok3 and i3['worst_ratio'] < 1 / 1.0001
    assert all(D3[j] <= 1.001 * i3['Etot'][j] * 1.0003 + 1e-12 for j in non)


def test_window_order_condition_uses_the_whole_window_of_the_previous_entry():
    prep, ent, J = _target(3)
    o, pi = T.best_cert(prep, 1e-9)
    ok, info, Dv = TL.solve_D(prep, T.PARS[pi], TL.h_of_eps(1e-9), o['Etot'], ladder=(1.1, 1.3, 1.6, 2.0), iters=3)
    assert ok
    ev = prep['L']['events']
    d1 = T.PARS[pi]['delta1']
    p0, p1 = ev[0]['p_e'], ev[1]['p_e']
    lo_ = max((d1 + info['De'][1]) / p1, 0.0)
    full = lo_ + max((info['De'][0] - d1) / p0, 0.0) + info['Tlay']
    short = lo_ + info['Tlay']
    assert full > short * 1.05
    gap = 0.5 * (full + short) * info['Lhi']                                          # between the correct and the (wrong) shortened requirement
    L2 = dict(prep['L'])
    L2['tau_br'] = np.array(prep['L']['tau_br'], float).copy()
    L2['tau_br'][2] = L2['tau_br'][1] + gap
    prep2 = dict(prep)
    prep2['L'] = L2
    ok2, i2 = TL.verify_tail(prep2, T.PARS[pi], Dv, TL.h_of_eps(1e-9))
    assert not ok2 and not i2['cond']['order']
    L2['tau_br'][2] = L2['tau_br'][1] + 1.05 * full * info['Lhi']
    ok3, i3 = TL.verify_tail(prep2, T.PARS[pi], Dv, TL.h_of_eps(1e-9))
    assert i3['cond']['order']


def test_golden_values_of_the_tail_verifier():
    # regression values of this implementation (floating-point candidates from `cert` at 1e-9, factor 1.5, parameters pi of best_cert): they pin the window gating, the pre-layer drift and the layer constants
    gold = {0: dict(drift=36.82093, pre_drift=[10.930649], dcmax=1.25674e-3, mass=2.33998e-4, T1=44.62521088, Tlay=117.4243647),
            3: dict(drift=5.143848, pre_drift=[11.019120, 11.011540], dcmax=1.04529e-3, mass=1.88851e-4, T1=44.67928226, Tlay=117.5518068)}
    for n, g in gold.items():
        prep, ent, J = _target(n)
        par, Dv, o = _cands(prep, ent, J, -9.0, 1.5)
        ok, info = TL.verify_tail(prep, par, Dv, TL.h_of_eps(1e-9))
        assert ok
        assert abs(info['drift_max'] / g['drift'] - 1) < 2e-3 and abs(info['dcmax'] / g['dcmax'] - 1) < 2e-3 and abs(info['mass'] / g['mass'] - 1) < 2e-3
        assert abs(info['T1'] / g['T1'] - 1) < 1e-7 and abs(info['Tlay'] / g['Tlay'] - 1) < 1e-7
        for a_, b_ in zip(info['pre_drift'], g['pre_drift']):
            assert abs(a_ / b_ - 1) < 1e-6


def test_kernel_margin_uses_the_rate_of_the_last_stage():
    prep, ent, J = _target(0)
    L = prep['L']
    K = prep['K']
    tau = 1.05 * L['events'][-1]['tau'] + 1.0
    j = max([j for j in range(J) if j not in ent], key=lambda j: abs(L['U_br'][K][j]))              # a kernel with a large dual at the last entry: the last stage's rate can reduce |u|
    sg = float(np.sign(L['U_br'][K][j]))
    for rate in (-0.02 * sg, 0.02 * sg, 0.0):
        L2 = dict(L)
        L2['P'] = np.array(L['P'], float).copy()
        L2['P'][K][j] = rate
        p2 = dict(prep)
        p2['L'] = L2
        u0 = abs(L['U_br'][K][j])
        uend = abs(L['U_br'][K][j] + rate * (tau - L['tau_br'][K]))
        want = 1 - max(u0, uend, max(abs(L['U_br'][i][j]) for i in range(K + 1)))
        assert abs(float(TL.kernel_margin(p2, j, tau).mid()) - want) < 1e-12


def test_kernel_margin_and_pieces():
    for n in (0, 3):
        prep, ent, J = _target(n)
        L = prep['L']
        K = prep['K']
        tau = 1.05 * L['events'][-1]['tau'] + 1.0
        non = [j for j in range(J) if j not in ent]
        cm = {j: float(TL.kernel_margin(prep, j, tau).mid()) for j in non}
        assert abs(min(cm.values()) - min(prep['mu_min'], min(1 - abs(L['U_br'][K][j] + L['P'][K][j] * (tau - L['tau_br'][K])) for j in non))) < 1e-12
        for j in non:
            ps = TL.dormant_pieces(prep, j, tau)
            assert abs(float(ps[0][1].mid())) == 0.0 and abs(float(ps[-1][2].mid()) - tau) < 1e-12 and all(float(p_[3].mid()) > 0 and float(p_[4].mid()) > 0 for p_ in ps)
            assert abs(min(min(float(p_[3].mid()), float(p_[4].mid())) for p_ in ps) - cm[j]) < 1e-9               # c_min is the smallest c at a piece end (independent of the formula)
            assert abs(cm[j] - (1 - max([abs(L['U_br'][i][j]) for i in range(K + 1)] + [abs(L['U_br'][K][j] + L['P'][K][j] * (tau - L['tau_br'][K]))]))) < 1e-12
            for a_, b_ in zip(ps[:-1], ps[1:]):
                assert abs(float(a_[2].mid()) - float(b_[1].mid())) < 1e-12 and abs(float(a_[4].mid()) - float(b_[3].mid())) < 1e-9     # contiguous
            # the exact integral of c^-3 over a linear piece against quadrature
            for (st, ta, tb_, ca, cb) in ps[:3]:
                t_ = np.linspace(float(ta.mid()), float(tb_.mid()), 20001)
                c_ = float(ca.mid()) + (float(cb.mid()) - float(ca.mid())) * (t_ - t_[0]) / (t_[-1] - t_[0])
                q_ = float(np.trapezoid(c_ ** -3.0, t_))
                f_ = float(((tb_ - ta) * (ca + cb) / (2 * ca * ca * cb * cb)).mid())
                assert abs(q_ / f_ - 1) < 1e-4


def test_tail_two_kernels_target_and_solve_D():
    prep, ent, J = _target(3)
    o, pi = T.best_cert(prep, 1e-9)
    ok, info, Dv = TL.solve_D(prep, T.PARS[pi], TL.h_of_eps(1e-9), o['Etot'], ladder=(1.1, 1.3, 1.6, 2.0), iters=3)
    assert ok and info['worst_ratio'] < 1 and len(info['De']) == 2
    ok2, i2 = TL.verify_tail(prep, T.PARS[pi], Dv, TL.h_of_eps(1e-9))
    assert ok2 and abs(i2['worst_ratio'] - info['worst_ratio']) < 1e-12
    # a window order violation: the second entry cannot be separated if its window starts before the first ends
    prep2 = dict(prep)
    L2 = dict(prep['L'])
    L2['tau_br'] = np.array(prep['L']['tau_br'], float).copy()
    L2['tau_br'][2] = L2['tau_br'][1] * (1 + 1e-9)
    prep2['L'] = L2
    ok3, i3 = TL.verify_tail(prep2, T.PARS[pi], Dv, TL.h_of_eps(1e-9))
    assert not ok3 and not i3['cond']['order']
    parz = dict(T.PARS[pi], delta1=0.0)                                         # with delta_1 = -15 the first lead (delta_1 + De)/p is negative and the first condition is vacuous
    oz = T.cert(prep, 1e-9, parz)
    okz, iz, Dz = TL.solve_D(prep, parz, TL.h_of_eps(1e-9), oz['Etot'], ladder=(1.1, 1.3, 1.6, 2.0), iters=3)
    assert okz and iz['cond']['order']
    L3 = dict(prep['L'])
    L3['tau_br'] = np.array(prep['L']['tau_br'], float).copy()
    L3['tau_br'][1] = 1e-12                                                    # the first entry comes before its lead time
    prep3 = dict(prep)
    prep3['L'] = L3
    ok4, i4 = TL.verify_tail(prep3, parz, Dz, TL.h_of_eps(1e-9))
    assert not ok4 and not i4['cond']['order']
    # K = 2: the tail dominates the sound point verifier on the dual deviations and on the displacement bound
    for e in (-9.0, -12.0):
        okp, ip = CH.verify_fixed(prep, arb(10.0 ** e), T.PARS[pi], Dv)
        ofl, _ = T.best_cert(prep, 10.0 ** e)
        assert okp and info['mass'] >= ofl['mass_ne'] * (1 - 1e-6)
        for j in info['Etot']:
            assert info['Etot'][j] >= ofl['Etot'][j] * (1 - 3e-3)          # the float certificate (0.1 % safety factors on its quadratures, legacy displacement accounting) is not above the proved bound
        assert info['dcmax'] >= ip['dcmax'] * 0.9


def test_displacement_accounting_is_sound_and_changes_only_the_drift():
    prep, ent, J = _target(3)
    par, Dv, o = _cands(prep, ent, J, -8.0, 1.5)
    ball = arb(10.0 ** -8.03).union(arb(10.0 ** -8.0))
    a, ia = RG.verify(prep, ball, par, Dv)
    b, ib = CH.verify_fixed(prep, ball, par, Dv)
    assert ib['dcmax'] > ia['dcmax'] * 1.05                                    # the legacy bound forgets the displacement accumulated in earlier stages
    assert ib['T1'] == ia['T1'] and ib['mass'] == ia['mass'] and ib['omega'] == ia['omega']
    for j in ia['Etot']:
        assert ib['Etot'][j] >= ia['Etot'][j]
    prep1, ent1, J1 = _target(0)
    par1, Dv1, o1 = _cands(prep1, ent1, J1, -9.0, 1.5)
    ball = arb(10.0 ** -9.03).union(arb(10.0 ** -9.0))
    a1, ia1 = RG.verify(prep1, ball, par1, Dv1)
    b1, ib1 = CH.verify_fixed(prep1, ball, par1, Dv1)
    assert a1 and b1 and ib1['dcmax'] >= ia1['dcmax']
    for old, new in CH.PATCH:                                                  # each anchor occurs exactly once in the verifier's source
        import inspect
        assert inspect.getsource(RG.verify).count(old) == 1 and old != new


def test_integrated_network_follows_inside_the_tail():
    for n, e_ok, eps_list in ((0, -9.0, (1e-9, 1e-11)), (3, -9.0, (1e-9, 1e-11))):
        prep, ent, J = _target(n)
        par, Dv, o = _cands(prep, ent, J, e_ok, 1.5)
        ok, info = TL.verify_tail(prep, par, Dv, TL.h_of_eps(10.0 ** e_ok))
        assert ok                                                              # the tail certificate holds for all eps <= 1e-9 on this target
        for eps in eps_list:
            sol, U = T.simulate(prep, eps)
            fo = T.follows(prep, eps, sol, U)
            assert fo[0] and fo[1]                                             # the network forms exactly the limit object's units with the exact fit


def test_sgn_dir_and_report_on_planted_rows():
    assert RUN.sgn_dir([3, 2, 1]) == 'non-increasing' and RUN.sgn_dir([1, 2, 3]) == 'non-decreasing' and RUN.sgn_dir([1, 1, 1]) == 'constant' and RUN.sgn_dir([1, 3, 2]) == 'non-monotone'
    assert RUN.sgn_dir([1.0, 1.0 + 1e-12, 1.0], 1e-9) == 'constant'
    raw = dict(U_rows=np.array([[-2.0, 100.0, 0.5], [-4.0, 100.0, 1.2]]))
    md = RUN.report(raw)
    assert 'All ratios are <= 1: NO' in md
    raw['U_rows'][1, 2] = 0.9
    assert 'All ratios are <= 1: yes' in RUN.report(raw)
    row = [0, 1, 0.0, 1.0, 1, -7.0, -7.5, -8.5, 0, 1e-4, 1e-5, 160.0, 170.0, 0.94, 70.0, 1e-3, 0.1, 1e-9, 0.06, 44.0, 44.5, 0.59, 0.6, 117.0, 85.0]
    raw2 = dict(T_rows=np.array([row]), T_extra=np.array([[0, -12.0, 1.0, 1e-9, 1e-10, 90.0, -24.0, 1.0, 1e-12, 1e-12, 90.0]]), C_rows=np.array([[0, -7.5, -8.5, 1.0, 10.0, 30.0, 0.0625, -8.5, 5.0]]))
    md2 = RUN.report(raw2)
    assert 'EVERY real eps in (0, eps_top]' in md2 and '| 1 | - | 1 | 1 | -7.00 | -7.50 | -8.50 |' in md2
    raw2['C_rows'][0, 3] = 0.0
    assert 'For 0 of 1 targets' in RUN.report(raw2)


@pytest.mark.skipif(not os.path.exists(RAW3), reason='the stored raw of theory3a_run.py is generated by --run')
def test_stored_report_regenerates_byte_identically_and_states_the_result():
    z = np.load(RAW3, allow_pickle=False)
    raw = {k: z[k] for k in z.files}
    md = RUN.report(raw)
    with open(os.path.join(os.path.dirname(RAW3), 'theory3a.md')) as fh:
        assert fh.read() == md
    assert raw['U_rows'][:, 2].max() <= 1.0
    T_ = raw['T_rows']
    assert len(T_) == 8 and np.all(T_[:, 13] < 1 / 1.0001)                     # every target has a tail with E_j/D_j below 1/1.0001
    C_ = raw['C_rows']
    assert np.all(C_[:, 3] > 0) and np.all(C_[:, 7] <= T_[:, 7] + 1e-9)        # every chain closes down to the tail
    # the stored bounds re-verify (one target): the tail verdict is reproduced from D_j
    n = int(T_[0, 0])
    prep, ent, J = _target(n)
    par = T.PARS[int(T_[0, 8])]
    ok, info = TL.verify_tail(prep, par, [float(x) for x in raw['T_D%d' % n]], TL.h_of_eps(10.0 ** T_[0, 7]))
    assert ok
