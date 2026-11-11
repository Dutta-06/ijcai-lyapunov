"""Round theory3, direction N: the corrected, merged verifications of Theorem T1 (experiments/theory3n_kernels.py, theory3n_system.py, theory3n_tail.py, theory3n_run.py).
Run explicitly:  pytest tests/theory3n_test.py   (a few minutes; the file name keeps it out of the default collection)."""
import inspect
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from flint import arb, ctx

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, os.path.join(REPO, 'tests'))
sys.path.insert(0, REPO)
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_wide as W        # noqa: E402
import theory_t1b_planar as PL     # noqa: E402
import theory3c_limit as LM        # noqa: E402
import theory3c_rigorous as R3     # noqa: E402
import theory3e_kernels as TK      # noqa: E402
import theory3n_kernels as NK      # noqa: E402
import theory3n_system as S        # noqa: E402
import theory3n_tail as TL         # noqa: E402
import theory3n_run as RUN         # noqa: E402

ctx.prec = 100


def clone(LR):
    out = dict(LR)
    out['events'] = [dict(e, w=list(e['w']), c=list(e['c'])) for e in LR['events']]
    out['tau_br'] = list(LR['tau_br'])
    out['U_br'] = [list(r) for r in LR['U_br']]
    out['P'] = [list(r) for r in LR['P']]
    return out


_CACHE = {}


def setup(n):
    """prep, limit object (balls), entrant list of the Q-row n (the 8 verified targets)"""
    if n not in _CACHE:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        LR = LM.limit_rig(cs, As, c0)
        _CACHE[n] = (prep, LR, [e['kernel'] for e in LR['events']])
    return _CACHE[n]


def Dvec(prep, ent, eps, fac):
    o, pi = S.best_cert_n(prep, eps)
    return o, pi, [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]


# ---------------------------------------------------------------------------------------------------------------------------------- the patches
def test_anchored_patches_reproduce_the_old_functions_when_switched_off():
    cs, As, c0 = T.make_target(2, 8.0, 1.0, 13)
    prep = T.prepare(cs, As, c0)
    leg = S.build_cert(True, True, True)
    for e in (-8.0, -9.0):
        a, b = T.cert(prep, 10.0 ** e, T.PARS[0]), leg(prep, 10.0 ** e, T.PARS[0])
        assert a['ok'] == b['ok'] and a['dcmax'] == b['dcmax'] and a['Dnew'] == b['Dnew'] and np.array_equal(a['Etot'], b['Etot'])
    pl = S.build_cert(True, True, True, planar=True)
    prep2 = PL.prepare2(*PL.make_target2(2, 8.0, 1.0, 1))
    a, b = PL.cert2d(prep2, 1e-8, T.PARS[0]), pl(prep2, 1e-8, T.PARS[0])
    assert a['ok'] == b['ok'] and a['dcmax'] == b['dcmax'] and a['Dnew'] == b['Dnew']
    prep, LR, ent = setup(3)
    o, pi, Dv = Dvec(prep, ent, 1e-9, 2.0)
    ball = R3.eps_ball(-9.01, -9.0)
    legr = S.build_rig(True, True, True, True, True)
    ok0, i0 = R3.verify_rig(prep, LR, ball, T.PARS[pi], Dv)
    ok1, i1 = legr(prep, LR, ball, T.PARS[pi], Dv)
    assert ok0 and ok1 and i0['dcmax'] == i1['dcmax'] and i0['Etot'] == i1['Etot'] and i0['cum'] == i1['cum']


def test_anchors_fail_loudly_if_the_source_changes():
    with pytest.raises(RuntimeError):
        S._apply('x = 1\n', S.F_FLOAT)
    with pytest.raises(RuntimeError):
        S._apply(S.F_FLOAT[0][0] * 2, S.F_FLOAT[:1])


def test_repair_f_the_displacement_integral_is_stage_resolved_float_and_rig():
    cs, As, c0 = T.make_target(2, 8.0, 1.0, 13)
    prep = T.prepare(cs, As, c0)
    leg, rf = S.build_cert(True, True, True), S.build_cert(False, True, True)
    o0, o1 = leg(prep, 1e-8, T.PARS[0]), rf(prep, 1e-8, T.PARS[0])
    assert o1['dcmax'] > 1.5 * o0['dcmax'] and o1['Dnew'] > o0['Dnew']                      # the old dcmax forgot the earlier stage's mass (factor 1.65 here)
    o0, o1 = leg(prep, 10 ** -7.5, T.PARS[0]), rf(prep, 10 ** -7.5, T.PARS[0])
    assert o0['ok'] and o0['dcmax'] < 0.3 < o1['dcmax']                                      # and it is the flip of the referee: S2 (dc = 0.3) fails at 1e-7.5
    # K = 1: the claim 'legacy and repaired coincide' of the earlier note is false (window 1 starts at t_1; stage 0 has no hd term)
    cs1, As1, c01 = T.make_target(1, 0.0, 1.0, 1)
    p1 = T.prepare(cs1, As1, c01)
    a, b = leg(p1, 1e-8, T.PARS[0]), rf(p1, 1e-8, T.PARS[0])
    assert b['dcmax'] > 1.3 * a['dcmax']
    # rig: same comparison on a ball
    prep, LR, ent = setup(4)
    o, pi, Dv = Dvec(prep, ent, 10 ** -8.0, 2.0)
    ball = R3.eps_ball(-8.01, -8.0)
    okL, iL = S.build_rig(True, True, True, True, True)(prep, LR, ball, T.PARS[pi], Dv)
    okF, iF = S.build_rig(False, True, True, True, True)(prep, LR, ball, T.PARS[pi], Dv)
    assert okL and iF["dcmax"] > 1.1 * iL["dcmax"] and max(iF["drift"]) > max(iL["drift"])


def test_repair_g_pre_layer_term_of_windows_not_started():
    cs, As, c0 = T.make_target(2, 8.0, 1.0, 13)
    prep = T.prepare(cs, As, c0)
    f_, fg = S.build_cert(False, True, True), S.build_cert(False, False, True)
    a, b = f_(prep, 1e-8, T.PARS[0]), fg(prep, 1e-8, T.PARS[0])
    assert b['dcmax'] > a['dcmax'] and b['Dnew'] > a['Dnew']
    assert b['dcmax'] < 1.1 * a['dcmax'] and b['Dnew'] < 1.01 * a['Dnew']                    # small but present
    assert b['comp'][5] > a['comp'][5]                                                          # in the drift component of E
    # the same in ball arithmetic
    prep, LR, ent = setup(4)
    o, pi, Dv = Dvec(prep, ent, 10 ** -8.0, 2.0)
    ball = R3.eps_ball(-8.01, -8.0)
    okF, iF = S.build_rig(False, True, True, True, True)(prep, LR, ball, T.PARS[pi], Dv)
    okG, iG = S.build_rig(False, False, True, True, True)(prep, LR, ball, T.PARS[pi], Dv)
    assert max(iG['drift']) > max(iF['drift']) and iG['dcmax'] >= iF['dcmax'] and max(iG['drift']) < 1.1 * max(iF['drift'])
    # the added term is psi_sup(X_li) pre per unstarted window: it is absent in the stage after the last window
    K1 = T.prepare(*T.make_target(1, 0.0, 1.0, 1))
    a, b = f_(K1, 1e-8, T.PARS[0]), fg(K1, 1e-8, T.PARS[0])
    assert b['Dnew'] >= a['Dnew'] and b['dcmax'] > a['dcmax']                                 # K = 1: the window is not started on the first piece


def test_repair_6_stage_maximum_of_the_old_units_coefficient():
    class E_(dict):
        pass
    ev = [dict(w=np.array([1.0])), dict(w=np.array([0.9, 2.0])), dict(w=np.array([0.95, 2.1, 3.0]))]
    assert np.array_equal(S._stage_max_vhat(ev, 3), [1.0, 2.1, 3.0])
    evb = [dict(w=[arb(1)]), dict(w=[arb(-0.9), arb(2)]), dict(w=[arb(0.95), arb(2.1), arb(-3)])]
    m = S._stage_max_vhat_rig(evb, 3)
    assert [round(float(R3.up(x)), 9) for x in m] == [1.0, 2.1, 3.0]
    prep, LR, ent = setup(3)
    ev_f = prep['L']['events']
    mf = S._stage_max_vhat(ev_f, 2)
    assert (mf >= np.abs(ev_f[-1]['w']) - 1e-18).all() and (mf > np.abs(ev_f[-1]['w'])).any()          # the stage value exceeds the final one (by ~5e-5 relative)
    mr = S._stage_max_vhat_rig(LR['events'], 2)
    assert all(R3.up(a) >= R3.up(abs(b)) for a, b in zip(mr, LR['events'][-1]['w']))
    cs, As, c0 = T.make_target(2, 6.0, 1.0, 11)
    p = T.prepare(cs, As, c0)
    a, b = S.build_cert(False, False, True)(p, 1e-8, T.PARS[0]), S.cert_n(p, 1e-8, T.PARS[0])
    assert b['Dnew'] >= a['Dnew'] and b['Dnew'] > a['Dnew'] * (1 + 1e-12)


def test_kappa_is_selected_by_the_lower_end_and_inv_of_has_no_float_gate():
    k1, k2 = arb(1.0, 0.5), arb(1.2, 0.01)                                  # [0.5, 1.5] and [1.19, 1.21]
    assert TL.min_lower(k1, k2) is k1 and TL.min_lower(k2, k1) is k1
    assert TL.min_lower(arb(2), arb(3)) == arb(2) and TL.min_lower(arb(3), arb(2)) == arb(2)
    # inv_of: certain ratio, None for a flat or undecided piece
    x, y, dt = arb(0.5), arb(0.2), arb(3)
    assert abs(float((S._inv_of(x, y, dt) - dt / arb(0.3)).abs_upper())) < 1e-20
    assert abs(float((S._inv_of(y, x, dt) - dt / arb(-0.3)).abs_upper())) < 1e-20
    assert S._inv_of(arb(0.5), arb(0.5), dt) is None
    assert S._inv_of(arb(0.5, 0.1), arb(0.5, 0.1), dt) is None
    assert S._inv_of(arb(0.5), arb(0.5) - arb(1e-8), dt) is None             # |difference| below 1e-6 (|x| + |y|)
    assert S._inv_of(arb(0.5), arb(0.5) - arb(1e-5), dt) is not None
    src = inspect.getsource(S.build_rig)
    assert 'float(' not in inspect.getsource(S._inv_of)


# ---------------------------------------------------------------------------------------------------------------------------------- flips of the float certificate
def test_the_two_flips_found_by_the_referee():
    prep = T.prepare(*T.make_target(2, 8.0, 1.0, 13))
    o_old, p_old = T.best_cert(prep, 10 ** -7.5)
    o_new, p_new = S.best_cert_n(prep, 10 ** -7.5)
    assert o_old['ok'] and not o_new['ok'] and p_new == -1
    assert S.best_cert_n(prep, 10 ** -8.0)[0]['ok'] and S.best_cert_n(prep, 10 ** -8.5)[0]['ok']
    prep = T.prepare(*W.make_target_k(6, 8.0, 1.5, 3))
    assert T.cert(prep, 10 ** -8.5, T.PARS[1])['ok'] and not S.cert_n(prep, 10 ** -8.5, T.PARS[1])['ok']
    o, p = S.best_cert_n(prep, 10 ** -9.0)
    assert o['ok']
    assert abs(S.cert_n(prep, 10 ** -8.5, T.PARS[1])['dcmax'] - 0.3253) < 2e-3        # legacy 0.299


def test_planar_certificate_uses_the_repairs():
    prep = PL.prepare2(*PL.make_target2(1, 0.0, 1.0, 1))
    a, b = PL.cert2d(prep, 1e-6, T.PARS[0]), S.cert2d_n(prep, 1e-6, T.PARS[0])
    assert a['ok'] and b['ok'] and b['dcmax'] > 1.3 * a['dcmax']


# ---------------------------------------------------------------------------------------------------------------------------------- g_atan, the unit in h
def test_g_atan_encloses_atan_over_y_with_certain_case_distinction():
    for y in (0.0, 1e-9, 0.01, 0.3, 0.4999, 0.5, 0.5001, 0.7, 2.0, 40.0):
        for r in (0.0, 1e-3):
            yb = arb(y).union(arb(y + r))
            g = TL.g_atan(yb)
            for yy in (y, y + r, y + r / 2):
                v = 1.0 if yy == 0 else math.atan(yy) / yy
                assert g.lower() <= v + 1e-15 and g.upper() >= v - 1e-15, (y, r, yy)
    for a_, b_ in ((0.4, 1.0), (0.3, 1.5), (0.0, 3.0), (0.45, 0.55)):                      # wide balls: the series is not used above 1/2 (it would not enclose)
        g = TL.g_atan(arb(a_).union(arb(b_)))
        hi_v, lo_v = (1.0 if a_ == 0 else math.atan(a_) / a_), math.atan(b_) / b_
        assert g.lower() <= lo_v and g.upper() >= hi_v
        assert float(g.upper() - g.lower()) < 1.01 * (hi_v - lo_v) + 1e-9                    # and tight (the monotone enclosure, not a blown-up series)
    assert float(TL.g_atan(arb(0).union(arb(1e-3))).rad()) < 1e-6
    # the branch is taken by a certain comparison: a ball straddling 1/2 uses the monotone enclosure
    from_zero = TL.g_atan(arb(0).union(arb(1.0)))                                           # a ball [0, 1] straddling both the branch and y = 0: g(0) = 1 is the upper end
    assert from_zero.upper() >= 1.0 and from_zero.lower() <= math.atan(1.0)
    wide = TL.g_atan(arb(0.3).union(arb(0.7)))
    assert wide.lower() <= math.atan(0.7) / 0.7 and wide.upper() >= math.atan(0.3) / 0.3


def test_h_parametrisation_and_join():
    for e in (1e-4, 1e-9, 1e-12, 3.3e-14):
        h = TL.h_of_eps(e)
        assert R3.lo(TL.eps_of_h(h)) >= arb(e)                                  # the tail ball (0, eps_of_h(h)] contains eps = e, decidably
        assert R3.up(TL.eps_of_h(h)) < arb(e) * (1 + arb(1e-9))
    assert abs(float(TL.eps_of_h(arb(0.3)).mid()) - 2 * (1 / 40) * 0.09 / 1.09) < 1e-12
    U = TL.UH(arb(0).union(TL.h_of_eps(1e-12)))
    assert U.h.lower() >= 0 and U.L.lower() >= 0 and U.L.upper() > 0


# ---------------------------------------------------------------------------------------------------------------------------------- the tail on the enclosed limit object
@pytest.mark.parametrize('n', [0, 3])
def test_tail_accepts_and_dominates_the_chain_verifier(n):
    prep, LR, ent = setup(n)
    eh = 1e-12
    o, pi, _ = Dvec(prep, ent, eh, 1.5)
    par = T.PARS[pi]
    hh = TL.h_of_eps(eh)
    ok, info, Dv = TL.solve_D_rig(prep, LR, par, hh, o['Etot'])
    assert ok and info['worst_ratio'] < 1 / 1.0001 and all(info['cond'].values())
    assert {'r>d0', 'sign', 'crit in D', 'unique in D', 'theta<v#'} <= set(info['cond'])
    assert info['eps_hi'] >= eh and info['eps_hi'] < 1.0001e-12
    # every eps <= eps_hi: the chain verifier on narrow balls with the SAME bounds D_j gives smaller numbers
    for e in (-12.0, -14.0, -20.0):
        ok2, i2 = S.verify_rig_n(prep, LR, R3.eps_ball(e - 0.01, e), par, Dv)
        assert ok2
        for j, v in i2['Etot'].items():                                   # both are upper bounds of the true E_j (envelope against exact dormant terms): no order, but the same size
            assert info['Etot'][j] >= v * (1 - 3e-3) - 1e-5 and info['Etot'][j] <= v * 1.3 + 1.0, (e, j)
        assert info['T1'] >= i2['T1'] * (1 - 1e-7) and info['Tlay'] >= i2['Tlay'] * (1 - 1e-7) and info['kappa'] <= i2['kappa'] * (1 + 1e-9)
        assert info['De'][0] >= i2['De'][0] * (1 - 1e-9)
    # the floating-point certificate (exact dormant mass) at eps <= eps_hi is below the tail's envelope numbers, and the envelope is not loose
    for e in (-12.0, -14.0):
        ofl = S.cert_n(prep, 10.0 ** e, par)
        assert info['mass'] >= ofl['mass_ne'] * (1 - 1e-6) and info['omega'] >= ofl['omega_peak'] * (1 - 1e-6)
        assert info['mass'] <= 1.01 * S.cert_n(prep, 1e-12, par)['mass_ne'] and info['omega'] <= 1.01 * S.cert_n(prep, 1e-12, par)['omega_peak']
        assert info['dcmax'] >= ofl['dcmax']
        for j, v in info['Etot'].items():
            assert v >= ofl['Etot'][j] * (1 - 3e-3) - 1e-5


def test_tail_layer_uses_the_stage_maximum_of_the_old_units_coefficient():
    prep, LR, ent = setup(3)
    o, pi, _ = Dvec(prep, ent, 1e-12, 1.5)
    ly = TL.layer(prep, LR, T.PARS[pi], arb(0).union(TL.h_of_eps(1e-12)))
    mr = S._stage_max_vhat_rig(LR['events'], 2)
    mf = S._stage_max_vhat(prep['L']['events'], 2)
    for a, b, c in zip(ly['vhat_fin'], mr, mf):
        assert R3.up(a) == R3.up(b) and R3.up(a) >= abs(arb(float(c))).upper() * (1 - 1e-12)
    assert R3.up(ly['vhat_fin'][0]) > R3.up(abs(LR['events'][-1]['w'][0]))                # strictly above the final block's value
    assert all(R3.up(ly['vhat'][k]) == R3.up(abs(LR['events'][k]['w'][-1])) for k in range(2))


def test_tail_rejects_what_it_must_reject():
    prep, LR, ent = setup(3)
    o, pi, _ = Dvec(prep, ent, 1e-12, 1.5)
    par = T.PARS[pi]
    hh = TL.h_of_eps(1e-12)
    ok, info, Dv = TL.solve_D_rig(prep, LR, par, hh, o['Etot'])
    assert ok
    non = [j for j in range(len(Dv)) if j not in ent]
    jmax = max(non, key=lambda j: info['Etot'][j])
    # (1) a bound D just below E x 1.0001 is refused (strict)
    Dbad = list(Dv)
    Dbad[jmax] = info['Etot'][jmax] * 1.00005
    ok1, i1 = TL.verify_tail_rig(prep, LR, par, Dbad, hh)
    assert not ok1 and not i1['cond']['fixed point E<D/1.0001']
    # (2) a D too large for the envelope (distance - D < 26)
    Dbig = list(Dv)
    Dbig[non[0]] = 1.0 / info['Lhi']                                           # L D = 1: the margin is gone
    ok2, i2 = TL.verify_tail_rig(prep, LR, par, Dbig, hh)
    assert not ok2 and 'why' in i2
    # (3) a much larger tail ball (eps_hi = 1e-4): the envelope needs lambda >= 25 and the layer/Omega conditions fail or the distance test fails
    ok3, i3 = TL.verify_tail_rig(prep, LR, par, Dv, TL.h_of_eps(1e-4))
    assert not ok3
    # (4) tampered limit data: hypotheses fail
    bad = clone(LR)
    bad['events'][0]['c'] = [x + arb(0.3) for x in bad['events'][0]['c']]
    ok4, i4 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok4 and not i4['cond']['crit in D']
    bad = clone(LR)
    bad['events'][1]['d0'] = arb(float(par['r']) + 0.01)
    ok4, i4 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok4 and not i4['cond']['r>d0']
    bad = clone(LR)
    bad['events'][0]['sign_ok'] = False
    ok4, i4 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok4 and not i4['cond']['sign']
    bad = clone(LR)
    bad['As'] = [x * arb(1.5) for x in bad['As']]
    ok4, i4 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok4 and not i4['cond']['unique in D']
    bad = clone(LR)
    bad['ok'] = False
    bad['why'] = 'test'
    ok4, i4 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok4 and 'why' in i4
    # (5) windows too close: the ordering condition (against the whole previous window, in L units)
    bad = clone(LR)
    bad['tau_br'][2] = bad['tau_br'][1] * arb(1.00000001)                                  # in tau units the windows are ~120 L = 3e-5 long at eps_hi = 1e-12
    ok5, i5 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok5 and 'cond' in i5 and not i5['cond']['order']
    # (6) a vanishing entry rate
    bad = clone(LR)
    bad['events'][0]['p_e'] = arb(1e-8)
    ok6, i6 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok6 and 'cond' in i6 and not i6['cond']['order']
    # (7) the window of the first entrant opens before t = 0 (K = 1: only the first ordering condition exists)
    p1, L1, e1 = setup(0)
    par1 = dict(T.PARS[0], delta1=0.0)                                                        # d_1 = 0: the window opens at tau_1 - (De_1 + 0)/p_1 (with d_1 = -15 the condition is void for K = 1)
    o1 = S.cert_n(p1, 1e-12, par1)
    ok7, i7, D7 = TL.solve_D_rig(p1, L1, par1, TL.h_of_eps(1e-12), o1['Etot'])
    assert ok7 and i7['cond']['order']
    bad = clone(L1)
    bad['events'][0]['p_e'] = arb(1e-8)
    ok7, i7 = TL.verify_tail_rig(p1, bad, par1, D7, TL.h_of_eps(1e-12))
    assert not ok7 and not i7['cond']['order']
    # (8) the gap between two entries is shorter than the previous window INCLUDING the layer duration (the layer is 117 / lambda long): between the two requirements
    d1 = float(par['delta1'])
    pe = [float(e['p_e'].mid()) for e in LR['events']]
    X_no_layer = max(0.0, (d1 + info['De'][1]) / pe[1]) + max(0.0, (info['De'][0] - d1) / pe[0])
    bad = clone(LR)
    bad['tau_br'][2] = bad['tau_br'][1] + arb(info['Lhi'] * (X_no_layer + 0.5 * info['Tlay']))
    ok8, i8 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert not ok8 and not i8['cond']['order']
    bad['tau_br'][2] = bad['tau_br'][1] + arb(info['Lhi'] * (X_no_layer + 1.5 * info['Tlay']))
    ok8, i8 = TL.verify_tail_rig(prep, bad, par, Dv, hh)
    assert i8['cond']['order']


def test_tail_keeps_the_pre_layer_term_of_windows_not_started():
    prep, LR, ent = setup(4)
    o, pi, _ = Dvec(prep, ent, 1e-8, 1.5)
    par = T.PARS[pi]
    hh = TL.h_of_eps(10 ** -8.5)
    ok, info, Dv = TL.solve_D_rig(prep, LR, par, hh, o['Etot'])
    assert ok
    src = inspect.getsource(TL._verify_tail)
    a = "hpre[i] if arb(gate[i]) >= eu else hd[i]"
    assert src.count(a) == 1
    ns = dict(TL.__dict__)
    exec(compile(src.replace(a, "arb(0) if arb(gate[i]) >= eu else hd[i]"), 'tail_no_g', 'exec'), ns)
    ok0, i0 = ns['_verify_tail'](prep, LR, par, Dv, hh, 1500, None)
    assert i0['drift_max'] < info['drift_max'] < 1.05 * i0['drift_max'] and i0['worst_ratio'] < info['worst_ratio']
    assert i0['dcmax'] <= info['dcmax']
    # without any gating (all windows counted from the start with hd) the bound is larger still
    ns2 = dict(TL.__dict__)
    exec(compile(src.replace(a, "hd[i]"), 'tail_all_hd', 'exec'), ns2)
    ok2, i2 = ns2['_verify_tail'](prep, LR, par, Dv, hh, 1500, None)
    assert i2['drift_max'] >= info['drift_max']


def test_tail_dormant_chord_is_a_lower_bound_of_the_margin_across_a_zero_crossing():
    prep, LR, ent = setup(3)
    non = [j for j in range(len(prep['c0'])) if j not in ent]
    TE = arb(R3.fup(arb(105) / 100 * LR['events'][-1]['tau'] + 1))
    for j in non:
        pcs = TL.dormant_pieces(LR, j, TE)
        cmin = TL.kernel_margin(LR, j, TE)
        assert cmin.rad() == 0 and 0 < cmin.upper() < 1
        K = LR['K']
        for (st, ta, tb_, ca, cb) in pcs:
            assert R3.lo(ca) > 0 and R3.lo(cb) > 0 and R3.lo(ca) >= cmin - arb(1e-25) and R3.lo(cb) >= cmin - arb(1e-25)
        # pieces tile [0, tau_end] contiguously
        assert abs(float(pcs[0][1].mid())) < 1e-30 and abs(float((pcs[-1][2] - TE).mid())) < 1e-30
        for p1, p2 in zip(pcs[:-1], pcs[1:]):
            assert (p1[2] - p2[1]).abs_upper() < 1e-30
    # a synthetic dual crossing zero at the quarter of a piece: u from 0.2 to -0.6
    fake = clone(LR)
    j = non[0]
    t1, t2 = fake['tau_br'][0], fake['tau_br'][1]
    fake['U_br'][0][j] = arb(0.2)
    fake['P'][0][j] = -(arb(0.8)) / (t2 - t1)
    pcs = TL.dormant_pieces(fake, j, TE)
    assert len(pcs) == len(TL.dormant_pieces(LR, j, TE)) + 1
    k0 = [p for p in pcs if p[0] == 0]
    assert len(k0) == 2 and abs(float((k0[0][2] - t1 - (t2 - t1) / 4).abs_upper())) < 1e-25       # the crossing at |u_a|/(|u_a| + |u_b|) = 1/4 of the piece
    assert abs(float(k0[0][3].mid()) - 0.8) < 1e-12 and abs(float(k0[0][4].mid()) - 1.0) < 1e-12 and abs(float(k0[1][3].mid()) - 1.0) < 1e-12 and abs(float(k0[1][4].mid()) - 0.4) < 1e-12
    # the margin is 1 - the largest |u| over ALL breakpoints and the horizon
    fake = clone(LR)
    fake['U_br'][LR['K']][j] = arb(0.9)                                         # largest at the last breakpoint
    assert abs(float(TL.kernel_margin(fake, j, TE).mid()) - 0.1) < 1e-12 and TL.kernel_margin(fake, j, TE).upper() <= arb(0.1)
    fake = clone(LR)
    fake['U_br'][LR['K']][j] = arb(0.2)
    fake['P'][LR['K']][j] = arb(0.7) / (TE - fake['tau_br'][LR['K']])            # |u| reaches 0.9 at the horizon
    assert abs(float(TL.kernel_margin(fake, j, TE).mid()) - 0.1) < 1e-12
    fake = clone(LR)
    fake['U_br'][LR['K']][j] = arb(0.9)
    fake['P'][LR['K']][j] = -arb(0.5) / (TE - fake['tau_br'][LR['K']])           # the last breakpoint is the maximum, the horizon value is 0.4
    assert abs(float(TL.kernel_margin(fake, j, TE).mid()) - 0.1) < 1e-12
    fake = clone(LR)
    fake['U_br'][1][j] = arb(0.6)
    assert float(TL.kernel_margin(fake, j, TE).mid()) <= 0.4 + 1e-12                            # an interior breakpoint also counts


# ---------------------------------------------------------------------------------------------------------------------------------- the chain and the covers
def test_cover_n_chain_contiguous_outward_rounded_and_every_ball_verified():
    prep, LR, ent = setup(0)
    segs, complete, cur, calls = S.cover_n(prep, LR, -7.5, -9.0)
    assert complete and abs(cur + 9.0) < 1e-12 and calls >= len(segs) > 0
    assert abs(segs[0][1] + 7.5) < 1e-12 and all(abs(a[0] - b[1]) < 1e-12 for a, b in zip(segs[:-1], segs[1:]))
    for lo_, hi in segs:
        b = R3.eps_ball(lo_, hi)
        assert b.lower() <= arb(10) ** arb(lo_) and b.upper() >= arb(10) ** arb(hi)
        o, pi = S.best_cert_n(prep, 10.0 ** hi)
        assert pi >= 0
    # an impossible top: nothing is accepted
    segs2, complete2, cur2, _ = S.cover_n(prep, LR, -3.0, -4.0)
    assert not complete2


GOLD_CHAIN = [   # (row, log10 eps, D factor): max E, dcmax, mass, omega, max drift, max cum, first window, V_tail, Env, T_1, T_2, T_3  (values of the corrected verifier; guards against silent changes)
    (0, -8.0, 1.5, 92.24347944312429, 0.017170794746565524, 0.0033719830783488602, 5.9066998320216036e-08, 6.762980005331739, 85.4771274509889, 156.68715722915118, 2.2469275251502636, 13.36926153778153, 44.62521286468, 0.03072296461323276, 72.76842982747118),
    (3, -9.0, 2.0, 108.99254861677235, 0.001107847708971404, 0.0013117450613706469, 1.2017145394742733e-09, 1.3334993611087516, 108.87659993185449, 155.57279975564188, 2.332487692865433, 13.559877902759812, 44.67928202788582, 0.030600102472817525, 72.84192296259053),
    (4, -9.0, 2.0, 107.52518977597973, 0.00885075373507416, 0.0021413182057769724, 1.2185037744850759e-08, 6.642845652207581, 107.41486016741928, 156.8200085649414, 2.2399221554905444, 13.370342803538197, 44.62551540152704, 0.03126786777283996, 72.76882982436409),
    (6, -9.5, 2.0, 130.5665976141946, 0.004986996663530392, 0.001288273714307359, 4.235697389873239e-09, 6.101750771999566, 130.44409434659624, 157.78862388285137, 2.245843277128597, 13.390509817662156, 44.631219863895765, 0.03137708269059658, 72.77647532871707),
    (8, -9.0, 2.0, 449.3410681316393, 0.0008679341283372708, 0.0055282611034937566, 1.3798234544075617e-11, 3.821838469482642, 445.51370116263456, 134.52810615976122, 25.937248695787, 84.98132527071982, 44.68688977713616, 0.02303258076426573, 72.85116196467573),
]


@pytest.mark.parametrize('case', GOLD_CHAIN)
def test_golden_values_of_the_corrected_chain_verifier(case):
    n, ex, f = case[:3]
    prep, LR, ent = setup(n)
    o, pi, Dv = Dvec(prep, ent, 10.0 ** ex, f)
    ball = R3.eps_ball(ex - 0.01, ex)
    ok, info = S.verify_rig_n(prep, LR, ball, T.PARS[pi], Dv)
    assert ok and info['worst_ratio'] < 1 and all(info['cond'].values())
    got = (max(info['Etot'].values()), info['dcmax'], info['mass'], info['omega'], max(info['drift']), max(info['cum']), info['wk'][0], info['vtail'], info['env'], info['T1'], info['T2'], info['T3'])
    for a, b in zip(got, case[3:]):
        assert abs(a / b - 1) < 1e-6, (got, case)
    okL, iL = R3.verify_rig(prep, LR, ball, T.PARS[pi], Dv)                                    # the legacy and the repaired quantities differ
    assert info['dcmax'] > iL['dcmax'] and max(info['drift']) > max(iL['drift'])
    assert max(info['Etot'].values()) >= max(iL['Etot'].values()) * (1 - 1e-9)


GOLD_TAIL = [   # (row, log10 eps_hi): worst E/D, dcmax, mass, omega, max drift, first De, Tlay, kappa, V_tail, max cum
    (0, -12.0, 0.9637653785255286, 1.0932939383851347e-06, 2.0874991177260713e-07, 3.9102040630949516e-14, 1.0378866916680205, 11.529572198558403, 117.42436278635981, 0.06312937804482936, 2.224785105037896, 85.47573965581715),
    (4, -8.5, 0.9998842492058658, 0.031739803163285826, 0.007054725537762063, 7.468167151308729e-08, 37.46791697203046, 12.527441864304105, 117.42561869950298, 0.06312259301546562, 2.274780980962834, 107.50012206128979),
    (6, -10.0, 0.9997030973990068, 0.00124203739537686, 0.00031981581887263977, 6.022207801777433e-10, 7.013931385981046, 12.77088507819087, 117.43907255724, 0.06300693861588047, 2.2390855857785654, 130.42884783768508),
    (8, -9.0, 0.9554406344976896, 3.942138048108153e-05, 3.2891419028195435e-06, 1.2265857282092263e-11, 1.6674222513146248, 13.20486502514624, 117.56108639114503, 0.059595080399885356, 25.93724495272906, 445.5170573475778),
]


@pytest.mark.parametrize('case', GOLD_TAIL)
def test_golden_values_of_the_tail_verifier(case):
    n, e = case[:2]
    prep, LR, ent = setup(n)
    o, pi, _ = Dvec(prep, ent, 10.0 ** e, 1.5)
    ok, info, Dv = TL.solve_D_rig(prep, LR, T.PARS[pi], TL.h_of_eps(10.0 ** e), o['Etot'])
    assert ok and info['worst_ratio'] < 1 / 1.0001 and all(info['cond'].values())
    got = (info['worst_ratio'], info['dcmax'], info['mass'], info['omega'], info['drift_max'], info['De'][0], info['Tlay'], info['kappa'], info['Vtail'], info['cum_max'])
    for a, b in zip(got, case[2:]):
        assert abs(a / b - 1) < 1e-9, (got, case)                                        # deterministic ball arithmetic: the x < 1e-9 part of the pre-layer drift is 7e-9 relative
    # the same bounds D_j, re-verified, give the same numbers (deterministic)
    ok2, i2 = TL.verify_tail_rig(prep, LR, T.PARS[pi], Dv, TL.h_of_eps(10.0 ** e))
    assert ok2 and i2['worst_ratio'] == info['worst_ratio']
    # parameters: a too small dc is refused by its own condition (the displacement bound dcmax is positive and larger)
    ok3, i3 = TL.verify_tail_rig(prep, LR, dict(T.PARS[pi], dc=info['dcmax'] * 1e-3), Dv, TL.h_of_eps(10.0 ** e))
    assert not ok3 and not i3.get('cond', {'dc': True})['dc']
    ok4, i4 = TL.verify_tail_rig(prep, LR, dict(T.PARS[pi], Om=info['omega'] * 0.5), Dv, TL.h_of_eps(10.0 ** e))
    assert not ok4 and not i4.get('cond', {'Omega': True})['Omega']


# ---------------------------------------------------------------------------------------------------------------------------------- the Gram tail below the table start (E-2)
def test_gram_tail_below_the_table_start_is_the_sup_over_the_true_range():
    g = TK.family('gauss', 1)
    x = np.linspace(0.0, 0.4, 4001)
    d2 = np.abs(g.fd2G(x))
    d1 = np.abs(g.fdG(x))
    G0 = np.abs(g.fG(x))
    assert NK.tail_old(g, 0, 0.0) < 0.97 < 1.0 <= g.tail(0, 0.0) < 1.0 + 1e-9                  # the old value was 0.961 at l = 0
    assert NK.tail_old(g, 2, 0.0) < 0.45 and g.tail(2, 0.0) >= 0.5 - 1e-12 and g.tail(2, 0.0) < 0.5002      # |Gamma''(0)| = 1/2
    for l in (0.0, 0.05, 0.2, 0.39):
        sel = x >= l
        assert g.tail(0, l) >= G0[sel].max() - 1e-12 and g.tail(1, l) >= d1[sel].max() - 1e-12 and g.tail(2, l) >= d2[sel].max() - 1e-12
    for o in (0, 1, 2):
        vals = [g.tail(o, l) for l in (0.0, 0.1, 0.2, 0.3, 0.399, 0.4, 0.5, 1.0, 2.0)]
        assert all(a >= b - 1e-12 for a, b in zip(vals[:-1], vals[1:]))                       # non-increasing in l
    for l in (0.4, 0.4001, 0.7, 5.0):
        for o in (0, 1, 2):
            assert g.tail(o, l) == NK.tail_old(g, o, l)                                       # unchanged from the table start on
    # the other kernels of the family (n = 1): Gamma(0) = 1, |Gamma''(0)| = 1/2 in sigma-equivalent units
    for kind in ('exp', 'mat32', 'mat52'):
        k = TK.family(kind, 1)
        assert k.tail(0, 0.0) >= 1.0 and NK.tail_old(k, 0, 0.0) < 1.0
        d2k = np.abs(k.fd2G(x))
        assert k.tail(2, 0.0) >= d2k.max() - 1e-9 and k.tail(2, 0.0) > NK.tail_old(k, 2, 0.0)
    # planar Bessel Matern Grams (raw formula singular at 0: the global fallback max(table, 1, M1, M2))
    for kind in ('exp',):
        kb = TK.family(kind, 2)
        assert kb.tail(0, 0.0) >= 1.0 and kb.tail(0, 0.0) <= 1.0 + 1e-9
        assert kb.tail(2, 0.0) >= max(NK.tail_old(kb, 2, 0.0), kb.M2()) - 1e-12 and kb.tail(1, 0.1) >= kb.M1() - 1e-12
        assert kb.tail('H', 0.2, 2) >= kb.M2(2) - 1e-12
    # Hessian norm for a planar Gaussian: the tangential eigenvalue Gamma/2 -> 1/2 at 0
    g2 = TK.family('gauss', 2)
    assert g2.tail('H', 0.0, 2) >= 0.5 - 1e-12 and NK.tail_old(g2, 'H', 0.0, 2) < 0.5
    # the certificate built for a Gram sees the corrected value at X = 0
    import theory3e_t1gen as G
    mod = G.build(g)
    assert mod.G2_sup(np.array([0.0]))[0] >= 0.5 - 1e-12 and mod.TAIL0(np.array([0.0]))[0] >= 1.0
    # the generic certificate of the Gaussian Gram now reproduces the repository's own Gaussian certificate (the old tail made D 3.6 % too small)
    prep = T.prepare(*T.make_target(2, 8.0, 1.0, 13))
    D_gen, D_rep = mod.cert(prep, 1e-8, T.PARS[0])['Dnew'], T.cert(prep, 1e-8, T.PARS[0])['Dnew']
    assert abs(D_gen / D_rep - 1) < 2e-3
    try:
        TK.Gram.tail = NK.tail_old
        mod_old = G.build(g)
        D_old = mod_old.cert(prep, 1e-8, T.PARS[0])['Dnew']
    finally:
        NK.install()
    assert D_old < 0.97 * D_rep


# ---------------------------------------------------------------------------------------------------------------------------------- Krawczyk refusals carried over from round theory3 C (the limit object is unchanged)
import theory3c_test as C3T        # noqa: E402


def test_carried_over_krawczyk_refusals():
    C3T.test_krawczyk_certifies_the_true_root_and_rejects_a_perturbed_non_root()
    C3T.test_krawczyk_does_not_certify_a_box_with_several_roots()
    C3T.test_floating_point_guess_far_from_any_root_fails_with_a_reason()
    C3T.test_tied_entry_times_are_refused()
    C3T.test_nearest_feature_tie_and_sign_test_are_refused()
    C3T.test_sign_test_sgn_p_equals_sgn_A_is_checked()
    C3T.test_critical_point_is_unique_in_the_whole_relaxation_set_dhat()


def test_hypotheses_fail_on_tampered_limit_data_in_the_corrected_chain_verifier():
    prep, LR, ent = setup(4)
    o, pi, Dv = Dvec(prep, ent, 1e-9, 2.0)
    par = T.PARS[pi]
    ball = R3.eps_ball(-9.01, -9.0)
    ok, info = S.verify_rig_n(prep, LR, ball, par, Dv)
    assert ok
    bad = clone(LR)
    bad['events'][0]['c'] = [x + arb(0.3) for x in bad['events'][0]['c']]
    ok, info = S.verify_rig_n(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    bad = clone(LR)
    bad['As'] = [x * arb(1.5) for x in bad['As']]
    ok, info = S.verify_rig_n(prep, bad, ball, par, Dv)
    assert not ok and info['cond']['crit in D'] and not info['cond']['unique in D']
    bad = clone(LR)
    bad['events'][0]['sign_ok'] = False
    ok, info = S.verify_rig_n(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['sign']
    bad = clone(LR)
    bad['ok'] = False
    bad['why'] = 'x'
    ok, info = S.verify_rig_n(prep, bad, ball, par, Dv)
    assert not ok and 'why' in info
    ok, info = S.verify_rig_n(prep, LR, R3.eps_ball(-9.0, -3.0), par, Dv)                   # a wide ball is refused, not unsoundly accepted
    assert not ok


# ---------------------------------------------------------------------------------------------------------------------------------- the report
def _raw_path():
    return os.path.join(REPO, 'results', 'synthetic', 'theory3n_raw.npz')


def test_report_is_regenerated_from_the_raw_byte_for_byte(tmp_path):
    raw = _raw_path()
    if not os.path.exists(raw):
        return
    shutil.copy(raw, tmp_path / 'theory3n_raw.npz')
    env = dict(os.environ, T3N_OUT=str(tmp_path), OMP_NUM_THREADS='1')
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory3n_run.py'), '--from-raw'], check=True, env=env, capture_output=True)
    assert open(tmp_path / 'theory3n.md', 'rb').read() == open(os.path.join(REPO, 'results', 'synthetic', 'theory3n.md'), 'rb').read()


def test_report_on_planted_rows_states_a_counterexample():
    ex = np.array([-6.0, -7.0])
    raw = dict(FC_ex=ex, FC_targets=np.array([[2, 8.0, 1.0, 13]]), FC_okc_old=np.array([[1.0, 1.0]]), FC_okc_new=np.array([[0.0, 1.0]]), FC_okm=np.array([[0.0, 0.0]]),
               FC_pi_old=np.zeros((1, 2)), FC_pi_new=np.zeros((1, 2)), FC_dc_old=np.full((1, 2), np.nan), FC_dc_same=np.full((1, 2), np.nan), FC_E_old=np.full((1, 2), np.nan),
               FC_E_same=np.full((1, 2), np.nan), FC_D_new=np.array([[np.nan, 100.0]]), FC_dc_new=np.full((1, 2), np.nan), F_net=np.array([[0, 0, -7.0, 1.0, 1.0, 0.0, 0.0]]))
    md = RUN.report(raw)
    assert 'violations among the corrected informative pairs: **1**' in md and 'NO - SEE ROWS' in md
    assert 'FK_ex' not in md and '| 2 | 8 | 1 | 13 | 1e-6.0 | 1e-7.0 | none |' in md


def test_stored_report_states_the_result():
    p = os.path.join(REPO, 'results', 'synthetic', 'theory3n_raw.npz')
    if not os.path.exists(p):
        return
    z = np.load(p, allow_pickle=False)
    Tr = z['E_rows']
    assert len(Tr) == 8
    for t in Tr:                                     # chain complete, tail accepted, join certain, conditions all true
        assert t[11] > 0 and t[14] > 0 and t[22] > 0 and t[23] > 0 and t[18] < 1 / 1.0001
    for key in 'CKP':
        oo, nn, om = z['F%s_okc_old' % key] > 0.5, z['F%s_okc_new' % key] > 0.5, z['F%s_okm' % key] > 0.5
        assert int((nn & ~om).sum()) == 0                                       # no violation among the corrected informative pairs
    net = z['F_net']
    assert all((r[4] < 0.5) or (r[6] > 0.5) for r in net)
    # the chains are contiguous from the top down to 10^ehi, and the tail ball starts at the chain's lowest ball end
    segs = z['E_segs']
    ehi = float(z['E_ehi'][0])
    for t in Tr:
        sg = segs[segs[:, 0] == t[0]]
        sg = sg[np.argsort(-sg[:, 2])]
        assert abs(sg[0, 2] - t[6]) < 1e-12 and abs(sg[-1, 1] - ehi) < 1e-12 and all(abs(a[1] - b[2]) < 1e-12 for a, b in zip(sg[:-1], sg[1:]))
        assert t[13] == R3.fdown(R3.pow10(ehi)) and t[17] >= t[13] and t[16] > 0
        assert abs(t[6] - t[5]) < 1e-12                                             # the tops of theory3c.md are confirmed
    if 'G_rows' in z.files:
        assert all(r[6] > 0 and abs(r[5] + 24.0) < 1e-12 for r in z['G_rows'])
    B = z['B_rows']
    assert B[:, 1].sum() == 364 and B[:, 2].sum() == 364 and B[:, 3].sum() == 362
    ex = z['FC_ex']
    assert RUN.thr_first(z['FC_okc_new'][7] > 0.5, ex) == -8.0 and RUN.thr_first(z['FK_okc_new'][14] > 0.5, z['FK_ex']) == -9.0     # the two flips of the referee
    n_old = sum(int((z['F%s_okc_old' % k] > 0.5).sum()) for k in 'CKP')
    n_new = sum(int((z['F%s_okc_new' % k] > 0.5).sum()) for k in 'CKP')
    assert (n_old, n_new) == (473, 471)
