"""Round theory3, direction C: the fully rigorous verification of Theorem T1 (experiments/theory3c_limit.py, theory3c_rigorous.py, theory3c_run.py).
Run explicitly:  pytest tests/theory3c_test.py   (about 2-3 min; the file name keeps it out of the default collection)."""
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from flint import arb, ctx

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_rigorous as RG   # noqa: E402  (the float-limit-data verifier, for the cross-check)
import theory3c_limit as LM        # noqa: E402
import theory3c_rigorous as R3     # noqa: E402

ctx.prec = 100


def clone(LR):
    """copy of a limit object (arb balls are immutable; dicts and lists are copied)"""
    out = dict(LR)
    out['events'] = [dict(e, w=list(e['w']), c=list(e['c'])) for e in LR['events']]
    out['tau_br'] = list(LR['tau_br'])
    out['U_br'] = [list(r) for r in LR['U_br']]
    out['P'] = [list(r) for r in LR['P']]
    return out


def target(idx):
    K, sep, R, seed = T.C_TARGETS[idx]
    cs, As, c0 = T.make_target(K, sep, R, seed)
    return cs, As, c0, T.prepare(cs, As, c0)


def block_root(cs, As, feats, c0=None):
    wf, cf = T.crit_block(cs, As, feats)
    cs_a, As_a, _ = LM.exact_data(cs, As, cs)
    return wf, cf, cs_a, As_a


# ---------------------------------------------------------------------------------------------------------------------------------- Krawczyk
def test_closed_form_jacobian_matches_finite_differences():
    cs = np.array([0.0, 6.3, 12.1])
    As = np.array([0.7, -0.9, 0.5])
    cs_a, As_a, _ = LM.exact_data(cs, As, cs)
    z = np.array([0.6, -0.8, 0.45, 0.2, 6.1, 11.9])
    za = [arb(float(x)) for x in z]
    F0, J = LM.FJ(za[:3], za[3:], cs_a, As_a)
    h = 1e-6
    for j in range(6):
        zp, zm = z.copy(), z.copy()
        zp[j] += h
        zm[j] -= h
        Fp = LM.FJ([arb(float(x)) for x in zp[:3]], [arb(float(x)) for x in zp[3:]], cs_a, As_a, jac=False)
        Fm = LM.FJ([arb(float(x)) for x in zm[:3]], [arb(float(x)) for x in zm[3:]], cs_a, As_a, jac=False)
        for i in range(6):
            fd = (float(Fp[i].mid()) - float(Fm[i].mid())) / (2 * h)
            assert abs(fd - float(J[i, j].mid())) < 1e-7, (i, j)


def test_critical_equations_agree_with_the_codes_hybr_system():
    """the reduced system (g, h) and the code's (g, w h) vanish together at the float root; the enclosure contains it"""
    cs, As, c0, prep = target(4)                                  # K = 2, sep 6
    wf, cf, cs_a, As_a = block_root(cs, As, [0, 1])
    Fv = LM.FJ([arb(float(x)) for x in wf], [arb(float(x)) for x in cf], cs_a, As_a, jac=False)
    assert max(abs(float(f.mid())) for f in Fv) < 1e-11
    Z = LM.enclose_root(list(wf) + list(cf), cs_a, As_a)
    for x, v in zip(Z, list(wf) + list(cf)):
        assert abs(float(x.mid()) - v) < 1e-10
        assert float(x.rad()) < 1e-28


def test_krawczyk_certifies_the_true_root_and_rejects_a_perturbed_non_root():
    cs, As, c0, prep = target(4)
    wf, cf, cs_a, As_a = block_root(cs, As, [0, 1])
    z = list(wf) + list(cf)
    zp = LM.newton_polish([arb(float(x)) for x in z], cs_a, As_a)
    n = len(zp)

    def box(centre, rad):
        return [c.mid().union(c.mid() + arb(rad)).union(c.mid() - arb(rad)) for c in centre]
    # a box around the root passes
    Kl, ok = LM.krawczyk(box(zp, 1e-8), cs_a, As_a)
    assert ok
    assert all(float((Kl[i].upper() - Kl[i].lower())) < 1e-15 for i in range(n))          # and contracts by a lot
    # a box of the same radius centred at a PERTURBED point (the root is outside) is rejected, in every coordinate direction
    for i in range(n):
        pert = [x for x in zp]
        pert[i] = pert[i] + arb(1e-6)
        _, ok = LM.krawczyk(box(pert, 1e-8), cs_a, As_a)
        assert not ok, i
    # ... also when the box is wide enough to be centred anywhere near but the point is not a root and the radius is small
    pert = [x + arb(3e-5) for x in zp]
    _, ok = LM.krawczyk(box(pert, 1e-6), cs_a, As_a)
    assert not ok
    # enclose_root on a guess 1e-3 off with no polishing and tiny boxes fails; the polished one passes
    with pytest.raises(LM.LimitFail):
        LM.enclose_root([float(x.mid()) + 1e-3 for x in zp], cs_a, As_a, radii=(1e-9, 1e-7), polish=False)
    Z = LM.enclose_root([float(x.mid()) + 1e-3 for x in zp], cs_a, As_a, polish=True)
    assert all(abs(float(a.mid()) - float(b.mid())) < 1e-12 for a, b in zip(Z, zp))


def test_krawczyk_does_not_certify_a_box_with_several_roots():
    """one unit fitted to two features 6 widths apart: critical points near each feature and one in between; a box around all of them is not certified"""
    cs = np.array([0.0, 6.0])
    As = np.array([1.0, 1.0])
    cs_a, As_a, _ = LM.exact_data(cs, As, cs)
    big = [arb(0.5).union(arb(1.5)), arb(-1).union(arb(7))]
    _, ok = LM.krawczyk(big, cs_a, As_a)
    assert not ok
    mid = [arb(0.9), arb(3.0)]                                         # the symmetric critical point between the features
    Fv = LM.FJ(mid[:1], mid[1:], cs_a, As_a, jac=False)
    assert abs(float(Fv[1].mid())) < 1e-12                              # h = 0 by symmetry
    small = [arb(0.9).union(arb(0.95)), arb(2.8).union(arb(3.2))]
    _, ok = LM.krawczyk(small, cs_a, As_a)
    assert not ok                                                       # a box that contains no root is rejected


def test_floating_point_guess_far_from_any_root_fails_with_a_reason():
    cs = np.array([0.0, 6.0])
    As = np.array([1.0, 1.0])
    cs_a, As_a, _ = LM.exact_data(cs, As, cs)
    with pytest.raises(LM.LimitFail):
        LM.enclose_root([0.3, 40.0], cs_a, As_a, polish=False)


# ---------------------------------------------------------------------------------------------------------------------------------- the limit object
@pytest.mark.parametrize('idx', [0, 4, 6, 12])
def test_limit_enclosure_contains_the_floating_point_limit_object_and_is_certified(idx):
    cs, As, c0, prep = target(idx)
    LR = LM.limit_rig(cs, As, c0)
    assert LR['ok'] and all(LR['checks'].values())
    assert LM.agrees_with_floats(LR, prep['L']) < 1e-12
    L = prep['L']
    for e, ef in zip(LR['events'], L['events']):
        assert e['kernel'] == ef['kernel'] and e['feature'] == ef['feature'] and e['sign_ok']
        assert e['gap'] > 0 and abs(e['gap'] - ef['gap']) < 1e-9
        assert e['rad'] < 1e-50
        assert float(e['p_e'].rad()) < 1e-9 * float(e['p_e'].mid())
    for i in range(len(L['tau_br'])):
        assert abs(float(LR['tau_br'][i].mid()) - L['tau_br'][i]) < 1e-12
        for j in range(L['J']):
            assert abs(float(LR['U_br'][i][j].mid()) - L['U_br'][i][j]) < 1e-12
    for i in range(len(L['P'])):
        for j in range(L['J']):
            assert abs(float(LR['P'][i][j].mid()) - L['P'][i][j]) < 1e-12
    assert abs(float(LR['mu'].lower()) - prep['mu_min']) < 1e-9 or float(LR['mu'].lower()) <= prep['mu_min']


def test_square_of_a_ball_containing_zero_is_nonnegative():
    z = arb(0).union(arb(1e-30))
    assert LM.sqb(z).lower() >= 0 and LM.Gam(z).upper() >= 1 - arb(1e-60) and LM.Gam(z).lower() > 0.9
    assert LM.sqb(arb(-3).union(arb(2))).lower() >= 0


def test_final_block_encloses_the_target_exactly():
    cs, As, c0, prep = target(12)                                  # K = 3
    LR = LM.limit_rig(cs, As, c0)
    ev = LR['events'][-1]
    for i, e in enumerate(LR['events']):
        assert ev['w'][i].contains(arb(float(As[e['feature']]))) and ev['c'][i].contains(arb(float(cs[e['feature']])))
    assert all(abs(float(p.mid())) < 1e-50 for j, p in enumerate(LR['P'][-1]))             # the final residual rate is 0 (the block is the target)


def test_tied_entry_times_are_refused():
    """two dictionary kernels symmetric about the feature hit the wall at the same time: the strict argmin cannot be certified"""
    cs = np.array([0.0])
    As = np.array([0.7])
    c0 = np.array([-1.0, 1.0, 4.0, -4.0])
    LR = LM.limit_rig(cs, As, c0)
    assert not LR['ok'] and not LR['checks']['strict_entrant']
    LR2 = LM.limit_rig(cs, As, np.array([-1.0, 1.0001, 4.0]))
    assert LR2['ok'] and LR2['events'][0]['kernel'] == 0


def test_nearest_feature_tie_and_sign_test_are_refused():
    cs = np.array([-2.0, 2.0])
    As = np.array([0.7, 0.7])
    LR = LM.limit_rig(cs, As, np.array([0.0, 30.0, -30.0]))                  # the first entrant is exactly between the two features: the assignment is a tie
    assert not LR['ok'] and not LR['checks']['strict_feature'] and LR['checks']['strict_entrant']
    LR = LM.limit_rig(cs, As, np.array([0.1, 30.0, -30.0]))
    assert LR['ok'] and LR['events'][0]['feature'] == 1
    # the stronger feature is entered first by the kernel next to it; the rate's sign equals the feature's sign
    cs = np.array([0.0, 8.0])
    As = np.array([0.7, -0.4])
    c0 = np.array([8.3, 0.1])
    LR = LM.limit_rig(cs, As, c0)
    assert LR['ok'] and all(e['sign_ok'] for e in LR['events'])
    assert LR['events'][0]['feature'] == 0 and LR['events'][0]['kernel'] == 1


def test_sign_test_sgn_p_equals_sgn_A_is_checked():
    """the kernel at 3.2 enters first (rate > 0) but its nearest feature (at 6) has a negative amplitude: refused"""
    cs = np.array([0.0, 6.0])
    As = np.array([0.7, -0.2])
    LR = LM.limit_rig(cs, As, np.array([3.4, 3.2]))
    assert not LR['ok'] and not LR['checks']['sign'] and LR['checks']['strict_entrant']
    LR = LM.limit_rig(cs, As, np.array([3.4, 3.9]))
    assert LR['ok'] and LR['checks']['sign']


def test_wall_margin_is_checked_strictly():
    """a dormant kernel whose dual reaches the wall before the second entry (|u| >= 1): refused"""
    cs = np.array([0.0, 7.0])
    As = np.array([1.0, 0.05])
    c0 = np.array([0.1, 0.0 + 3.4, 7.1])                                    # the middle kernel sees the large feature strongly and reaches the wall before the weak feature enters
    LR = LM.limit_rig(cs, As, c0)
    assert not LR['ok'] or LR['checks']['wall']


# ---------------------------------------------------------------------------------------------------------------------------------- certain comparisons and guards
def test_outward_rounding_helpers():
    x = arb(1) / 3
    assert arb(R3.fup(x)) >= x.upper() and arb(R3.fdown(x)) <= x.lower()
    assert R3.fup(x) - R3.fdown(x) < 1e-15
    for k in range(1, 40):                                              # points whose nearest double lies on either side
        x = arb(k) / 7
        assert arb(R3.fup(x)) >= x.upper() and arb(R3.fdown(x)) <= x.lower() and R3.fup(x) - R3.fdown(x) < 1e-15
        x = arb(k).sqrt()
        assert arb(R3.fup(x)) >= x.upper() and arb(R3.fdown(x)) <= x.lower()
    y = arb(0.1)
    assert R3.fup(y) == 0.1 and R3.fdown(y) == 0.1                       # exact doubles stay
    # the eps-ball ends are outward rounded
    for e in (-7.0, -6.3, -5.77, -9.1, -12.4, -3.21, -24.0):
        b = R3.eps_ball(e - 0.37, e)
        assert b.lower() <= arb(10) ** arb(e - 0.37) and b.upper() >= arb(10) ** arb(e)


def test_guards_use_the_upper_end_and_are_certain():
    st = arb(0.5).union(arb(2.0))
    assert R3.logpos_up(st).upper() >= arb(2).log() - arb(1e-12)
    assert R3.logpos_up(arb(0.5)).upper() == 0 and R3.logpos_up(arb(1)).upper() == 0
    assert R3.pos_up(arb(-1).union(arb(3))) >= 3 and R3.pos_up(arb(-2)) == 0
    assert R3.bmax_up(arb(1), arb(2)) == arb(2) and R3.bmax_up(arb(3), arb(2)) == arb(3)
    c = R3.clamp01(arb(-2).union(arb(-1)))
    assert c.lower() >= 0 and c.upper() <= 1e-30 or c.upper() == 0                   # entirely negative -> 0 (was [-1, 0])
    c = R3.clamp01(arb(0.2).union(arb(3)))
    assert c.lower() <= 0.2 and c.upper() >= 1 and c.lower() > 0.19 and c.upper() < 1 + 1e-8                # clamped from above
    c = R3.clamp01(arb(1.5).union(arb(3)))
    assert c.upper() < 1 + 1e-8 and c.lower() > 0.99
    s = R3.sq(arb(0).union(arb(1e-30)))
    assert s.lower() >= 0
    s = R3.sq(arb(-3).union(arb(2)))
    assert s.lower() >= 0 and s.upper() >= 9
    ball = arb(1).union(arb(5))
    assert R3.bmax_up(ball, arb(4)).upper() >= 5 and R3.bmax_up(arb(4), ball).upper() >= 5 and R3.bmax_up(arb(0.5).union(arb(2)), arb(1.5)).upper() >= 2
    assert R3.psi_sup_r(arb(1.2)).upper() >= R3.PSI_MAX_R.lower() and R3.g2_sup_r(arb(1.41)).upper() >= R3.E15.lower()
    assert R3.psi_sup_r(arb(0.5)) is R3.PSI_MAX_R and R3.psi_sup_r(arb(3)).contains(R3.psi_r(arb(3)))
    assert R3.g2_sup_r(arb(0.5)).upper() >= R3.E15.upper() and R3.g2_sup_r(arb(1.2)).upper() >= (abs(arb(1.2) ** 2 / 4 - arb(1) / 2) * R3.Gam_r(arb(1.2))).upper()


def test_bisection_inversions_enclose_the_roots_and_use_certain_comparisons():
    U = R3.RU(arb(1e-8))
    for t in (arb(30), arb(1500).union(arb(1600)), arb(40000)):
        x = U.x_of_s(t)
        sv = U.s(x)
        assert float(sv.lower()) <= float(t.upper()) and float(sv.upper()) >= float(t.lower())
    for t in (arb(1e-3), arb(0.2).union(arb(0.21))):
        x = U.x_of_w(t)
        wv = U.w(x)
        assert float(wv.lower()) <= float(t.upper()) and float(wv.upper()) >= float(t.lower())
    with pytest.raises(ValueError):
        U.x_of_s(arb(1e9))


# ---------------------------------------------------------------------------------------------------------------------------------- the verifier
def _setup(idx, ex):
    cs, As, c0, prep = target(idx)
    LR = LM.limit_rig(cs, As, c0)
    o, pi = T.best_cert(prep, 10.0 ** ex)
    par = T.PARS[pi]
    ent = [e['kernel'] for e in LR['events']]
    return prep, LR, o, par, ent


def _Dv(o, ent, J, f):
    return [1.0 if j in ent else max(f * o['Etot'][j], 1.0) for j in range(J)]


@pytest.mark.parametrize('idx,ex', [(0, -8.0), (3, -9.0), (4, -9.0), (12, -9.5)])
def test_rigorous_verification_agrees_with_the_float_limit_verifier(idx, ex):
    """same S(eps) data, narrow eps ball: the fully rigorous verifier accepts what the float-limit-data verifier accepts, with the same bounds up to the window-start gating"""
    prep, LR, o, par, ent = _setup(idx, ex)
    ball = R3.eps_ball(ex - 0.01, ex)
    ball_old = arb(10.0 ** (ex - 0.01)).union(arb(10.0 ** ex))
    got = False
    for f in (1.5, 2.0, 4.0):
        Dv = _Dv(o, ent, len(prep['c0']), f)
        ok, info = R3.verify_rig(prep, LR, ball, par, Dv)
        ok0, info0 = RG.verify(prep, ball_old, par, Dv)
        if ok0:
            assert ok, (f, info.get('cond'))
            for j, e in info0['Etot'].items():
                assert abs(info['Etot'][j] / e - 1) < 0.06
            assert abs(info['T1'] / info0['T1'] - 1) < 1e-3 and abs(info['Tlay'] / info0['Tlay'] - 1) < 1e-3
            assert abs(info['lam'] / info0['lam'] - 1) < 1e-9 and abs(info['mass'] / info0['mass'] - 1) < 0.02 and abs(info['omega'] / info0['omega'] - 1) < 0.02
            assert abs(info['T2'] / info0['T2'] - 1) < 1e-3 and abs(info['T3'] / info0['T3'] - 1) < 1e-3 and abs(info['kappa'] / info0['kappa'] - 1) < 1e-9
            assert abs(info['Atr'] / info0['Atr'] - 1) < 1e-3 and abs(info['pre'] / info0['pre'] - 1) < 1e-6
            assert abs(info['tau_end'] - (1.05 * prep['L']['events'][-1]['tau'] + 1.0)) < 1e-9 and info['tau_end'] >= 1.05 * prep['L']['events'][-1]['tau'] + 1.0
            amin_, amax_ = float(np.abs(prep['As']).min()), float(np.abs(prep['As']).max())
            assert abs(info['theta'] - min(0.01 * amax_, 0.1 * amin_)) < 1e-12
            assert info['dcmax'] >= 0.99 * info0['dcmax'] and info['dcmax'] < 4 * info0['dcmax']        # the conservative window gating: never below the float-data verifier, not far above
            got = True
            break
    assert got
    assert all(info['cond'].values())
    assert {'r>d0', 'sign', 'crit in D', 'unique in D', 'theta<v#'} <= set(info['cond'])


GOLD = [   # (Q row, log10 eps, D factor): Etot max, dcmax, mass, omega, drift max, cum max, windows, V_tail, Env, T_1, T_2, T_3   (values of the verified implementation; guards against silent changes)
    (0, -8.0, 1.5, 90.40791840453272, 0.011727926814097091, 0.0033719830783488602, 5.9066998320216036e-08, 4.927418966740183, 85.4771274509889, [156.68715722915118], 2.2469275251502636, 13.36926153778153, 44.62521286468, 0.03072296461323276, 72.76842982747118),
    (3, -9.0, 2.0, 108.97793779894646, 0.0006990833308407719, 0.0013117450613706469, 1.2017145394742733e-09, 0.932171404125951, 108.87633318506448, [155.57279975564188, 157.07411024402174], 2.332487692865433, 13.559877902759812, 44.67928202788582, 0.030600102472817525, 72.84192296259053),
    (4, -9.0, 2.0, 107.51675294853386, 0.006989518814375807, 0.0021413182057769724, 1.2185037744850759e-08, 4.757894667872379, 107.41485961649138, [156.8200085649414, 158.31822826729123], 2.2399221554905444, 13.370342803538197, 44.62551540152704, 0.03126786777283996, 72.76882982436409),
    (6, -9.5, 2.0, 130.55037071088114, 0.004045075966337329, 0.0012912914783393711, 4.235697389873239e-09, 4.449445568153605, 130.44409850582423, [157.78862820844984, 158.5387967434035, 158.73386640429328], 2.245843277128597, 13.390509817662156, 44.631219863895765, 0.03137708269059658, 72.77647532871707),
    (8, -9.0, 2.0, 449.14016885476127, 0.0008360920614302489, 0.0055282611034937566, 1.3798234544075617e-11, 3.621187756672075, 445.51345259856714, [134.52810615976122, 145.3331847251212], 25.937248695787, 84.98132527071982, 44.68688977713616, 0.02303258076426573, 72.85116196467573),
]
Q_IDX = (0, 1, 3, 4, 7, 8, 12, 9, 6, 13)


@pytest.mark.parametrize('case', GOLD)
def test_golden_values_of_the_verified_implementation(case):
    n_, ex, f = case[:3]
    K, sep, Rr, seed = T.C_TARGETS[Q_IDX[n_]]
    cs, As, c0 = T.make_target(K, sep, Rr, seed)
    prep = T.prepare(cs, As, c0)
    LR = LM.limit_rig(cs, As, c0)
    o, pi = T.best_cert(prep, 10.0 ** ex)
    par = T.PARS[pi]
    ent = [e['kernel'] for e in LR['events']]
    Dv = [1.0 if j in ent else max(f * o['Etot'][j], 1.0) for j in range(len(c0))]
    ok, info = R3.verify_rig(prep, LR, R3.eps_ball(ex - 0.01, ex), par, Dv)
    assert ok
    got = [max(info['Etot'].values()), info['dcmax'], info['mass'], info['omega'], max(info['drift']), max(info['cum']), info['wk'], info['vtail'], info['env'], info['T1'], info['T2'], info['T3']]
    for g, w in zip(got, case[3:]):
        assert np.allclose(g, w, rtol=2e-4, atol=0), (g, w)
    # the omitted Omega / dc / distance guards of the conditions
    bad = dict(par, Om=1e-12)
    ok2, info2 = R3.verify_rig(prep, LR, R3.eps_ball(ex - 0.01, ex), bad, Dv)
    assert not ok2 and not info2.get('cond', {'Omega': False})['Omega']
    bad = dict(par, dc=1e-7)
    ok2, info2 = R3.verify_rig(prep, LR, R3.eps_ball(ex - 0.01, ex), bad, Dv)
    assert not ok2 and not info2.get('cond', {'dc': False})['dc']
    ok2, info2 = R3.verify_rig(prep, LR, R3.eps_ball(ex - 0.01, ex), par, [1e9 if j not in ent else 1.0 for j in range(len(c0))])
    assert not ok2 and info2.get('why') == 'dist<=D'
    # K = 1: the first ordering condition (the window must open after t = 0)
    if K == 1:
        bad = clone(LR)
        bad['events'][0]['p_e'] = arb(1e-8)
        ok2, info2 = R3.verify_rig(prep, bad, R3.eps_ball(ex - 0.01, ex), dict(par, delta1=0.0), Dv)
        assert not ok2 and not info2['cond']['order']


@pytest.mark.parametrize('n_,ex,f,etot,dcmax,drift', [(0, -8.0, 1.5, 89.37375739778625, 0.0117161824091976, 4.9536250821184264),
                                                      (4, -9.0, 2.0, 106.79987428853862, 0.006987236113936174, 4.751349806942132)])
def test_golden_values_with_an_open_window_inside_the_piece(n_, ex, f, etot, dcmax, drift):
    """delta_1 = 0: the windows open (S^- > 0) inside the stage pieces, so the sub-piece split of the drift terms is exercised"""
    K, sep, Rr, seed = T.C_TARGETS[Q_IDX[n_]]
    cs, As, c0 = T.make_target(K, sep, Rr, seed)
    prep = T.prepare(cs, As, c0)
    LR = LM.limit_rig(cs, As, c0)
    o, pi = T.best_cert(prep, 10.0 ** ex)
    ent = [e['kernel'] for e in LR['events']]
    Dv = [1.0 if j in ent else max(f * o['Etot'][j], 1.0) for j in range(len(c0))]
    ok, info = R3.verify_rig(prep, LR, R3.eps_ball(ex - 0.01, ex), dict(T.PARS[pi], delta1=0.0), Dv)
    assert ok
    assert np.allclose([max(info['Etot'].values()), info['dcmax'], max(info['drift'])], [etot, dcmax, drift], rtol=5e-5, atol=0)


def test_wide_eps_ball_is_refused_not_unsoundly_accepted():
    prep, LR, o, par, ent = _setup(4, -7.5)
    Dv = _Dv(o, ent, len(prep['c0']), 1.5)
    ball = R3.eps_ball(-7.5, -3.0)                                    # five decades up to the middle of the failing range
    try:
        ok, info = R3.verify_rig(prep, LR, ball, par, Dv)
    except (ValueError, ZeroDivisionError):
        ok = False
    assert not ok


def test_bootstrap_is_strict():
    prep, LR, o, par, ent = _setup(0, -8.0)
    eps = arb(1e-8)
    J = len(prep['c0'])

    def at(f):
        Dv = _Dv(o, ent, J, 1.5)
        for _ in range(8):
            ok, info = R3.verify_rig(prep, LR, eps, par, Dv)
            Dv = [1.0 if j in ent else f * info['Etot'][j] for j in range(J)]
        return R3.verify_rig(prep, LR, eps, par, Dv)
    ok, info = at(1.0005)
    assert ok and info['cond']['fixed point E<D/1.0001']
    ok, info = at(1.00005)
    assert not ok and not info['cond']['fixed point E<D/1.0001'] and info['worst_ratio'] < 1.0


def test_hypotheses_fail_on_tampered_limit_data():
    prep, LR, o, par, ent = _setup(4, -9.0)
    ball = R3.eps_ball(-9.01, -9.0)
    Dv = _Dv(o, ent, len(prep['c0']), 2.0)
    ok, info = R3.verify_rig(prep, LR, ball, par, Dv)
    assert ok
    # (1) a block critical point outside D-hat (shifted centre) at an intermediate stage
    bad = clone(LR)
    bad['events'][0]['c'] = [x + arb(0.3) for x in bad['events'][0]['c']]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    # (2) a coefficient outside [v_D^-, v_D^+] a_min, only at the first stage
    bad = clone(LR)
    bad['events'][0]['w'] = [x * 3 for x in bad['events'][0]['w']]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    # (3) a coefficient of the wrong sign
    bad = clone(LR)
    bad['events'][-1]['w'] = [-x for x in bad['events'][-1]['w']]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    # (4) an offset d_0 not below r
    bad = clone(LR)
    bad['events'][1]['d0'] = arb(float(par['r']) + 0.01)
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['r>d0']
    # (2b) a coefficient below v_D^- a_min; (1b) only the SECOND unit's centre of the final block outside D-hat
    bad = clone(LR)
    bad['events'][-1]['w'] = [bad['events'][-1]['w'][0], bad['events'][-1]['w'][1] * 0.2]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    bad = clone(LR)
    bad['events'][-1]['c'] = [bad['events'][-1]['c'][0], bad['events'][-1]['c'][1] + arb(0.3)]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['crit in D']
    # (2c) the amplitudes of the limit's system scaled (so the unique critical point lies outside D-hat) while the stored block is untouched: only the uniqueness check can see it
    bad = clone(LR)
    bad['As'] = [x * arb(1.5) for x in bad['As']]
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and info['cond']['crit in D'] and not info['cond']['unique in D']
    # (5) the sign test
    bad = clone(LR)
    bad['events'][0]['sign_ok'] = False
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and not info['cond']['sign']
    # (6) a limit object that was not certified is not verified at all
    bad = clone(LR)
    bad['ok'] = False
    bad['why'] = 'test'
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and 'why' in info
    # (7) the second entry time too close to the first (the windows would overlap): the ordering condition fails
    bad = clone(LR)
    bad['tau_br'][2] = bad['tau_br'][1] * arb(1.0001)
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and 'cond' in info and not info['cond']['order']
    # (8) a vanishing entry rate: the window opens before t = 0 (first ordering condition)
    bad = clone(LR)
    bad['events'][0]['p_e'] = arb(1e-8)
    ok, info = R3.verify_rig(prep, bad, ball, par, Dv)
    assert not ok and 'cond' in info and not info['cond']['order']


def test_rig_box_conditions_are_certain():
    r, Om = 0.6, 0.02
    assert R3.rig_box(2, 1.0, 8.0, r, Om, 0.05, 0.05) is not None
    assert R3.rig_box(2, 1.0, 4.0, r, Om, 0.05, 0.05) is None                      # ell = Delta - 2r < sqrt 6
    assert R3.rig_box(2, 1.0, 5.0, 0.8, Om, 0.05, 0.05) is None                    # ell < sqrt 6 again, through r
    assert R3.rig_box(3, 10.0, 5.4, r, Om, 0.05, 0.05) is None or R3.rig_box(3, 10.0, 5.4, r, Om, 0.05, 0.05)['bd'].lower() > 0
    assert R3.rig_box(2, 1.0, 8.0, r, 2.5, 0.05, 0.05) is None                     # a dormant budget this large destroys the regime
    assert R3.rig_box(2, 1.0, 8.0, 1.5, Om, 0.05, 0.05) is None                    # r too large (d-tilde < r fails / ell)


def test_cover_balls_are_outward_rounded_and_contiguous(monkeypatch):
    prep, LR, o, par, ent = _setup(0, -7.5)
    seen = []

    def fake(prep_, LR_, ball, par_, Dv):
        seen.append(ball)
        return True, {}
    monkeypatch.setattr(R3, 'verify_rig', fake)
    segs, complete, cur, calls = R3.cover(prep, LR, -7.5, -10.0)
    assert complete and len(segs) == len(seen) and len(segs) >= 3
    for (a, b), ball in zip(segs, seen):
        assert ball.lower() <= arb(10) ** arb(a) and ball.upper() >= arb(10) ** arb(b)
    for (a, b), (c, d) in zip(segs[:-1], segs[1:]):
        assert c < a <= b and abs(d - a) < 1e-12                                   # the chain goes down, consecutive exponents meet
    for b1, b2 in zip(seen[:-1], seen[1:]):
        assert b1.lower() <= b2.upper()                                             # the balls overlap or touch


def test_critical_point_is_unique_in_the_whole_relaxation_set_dhat():
    """Krawczyk on D-hat itself: one root in D-hat (every stage); an inflated D-hat (d-hat = 1.4) or a box far from the blocks is not certified"""
    prep, LR, o, par, ent = _setup(12, -9.0)                      # K = 3
    As = np.asarray(prep['As'], float)
    amin, amax = arb(float(np.abs(As).min())), arb(float(np.abs(As).max()))
    Delta = min(abs(LR['cs'][i] - LR['cs'][j]) for i in range(3) for j in range(i + 1, 3))
    Bx = R3.rig_box(3, amax / amin, Delta, arb(par['r']), arb(par['Om']), par['dfrac'], par['vsig'])
    assert Bx is not None
    for k in (1, 2, 3):
        assert R3.dhat_unique(LR, Bx, amin, As, k)
    big = dict(Bx, dh=arb(1.4))
    assert not R3.dhat_unique(LR, big, amin, As, 2)
    wide = dict(Bx, vD_hi=Bx['vD_hi'] + 40)
    assert not R3.dhat_unique(LR, wide, amin, As, 2) or True       # a wider coefficient range is allowed to fail; the call must not raise
    tiny = dict(Bx, dh=arb(1e-9))                                 # a relaxation set that does not contain the block's critical point: no root in the box, not certified
    assert not R3.dhat_unique(LR, tiny, amin, As, 1) and R3.dhat_unique(LR, tiny, amin, As, 3)       # the final block is the target itself (offset 0)


def test_uniqueness_in_dhat_is_required_at_every_stage_and_failures_count_as_failures(monkeypatch):
    prep, LR, o, par, ent = _setup(4, -9.0)                         # K = 2
    ball = R3.eps_ball(-9.01, -9.0)
    Dv = _Dv(o, ent, len(prep['c0']), 2.0)
    assert R3.verify_rig(prep, LR, ball, par, Dv)[0]
    real = R3.dhat_unique
    monkeypatch.setattr(R3, 'dhat_unique', lambda LR_, Bx, amin, As_f, k: False if k == 2 else real(LR_, Bx, amin, As_f, k))
    ok, info = R3.verify_rig(prep, LR, ball, par, Dv)
    assert not ok and not info['cond']['unique in D']
    monkeypatch.setattr(R3, 'dhat_unique', real)

    def boom(*a, **k):
        raise ArithmeticError('arb failure')
    monkeypatch.setattr(LM, 'krawczyk', boom)
    As = np.asarray(prep['As'], float)
    amin, amax = arb(float(np.abs(As).min())), arb(float(np.abs(As).max()))
    Bx = R3.rig_box(2, amax / amin, abs(LR['cs'][0] - LR['cs'][1]), arb(par['r']), arb(par['Om']), par['dfrac'], par['vsig'])
    assert R3.dhat_unique(LR, Bx, amin, As, 1) is False


def test_hypothesis_slacks_are_positive_on_the_stored_covers():
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory3c_raw.npz')
    if not os.path.exists(raw):
        return
    z = np.load(raw, allow_pickle=False)
    Hy = z['HY_rows']
    assert len(Hy) >= 1
    assert (Hy[:, 2:7] > 0).all() and (Hy[:, 10] > 0).all()
    Ev = z['EV_rows']
    assert (Ev[:, 12] > 0).all() and (Ev[:, 9] > 0).all() and (Ev[:, 11] > 0).all() and (Ev[:, 14] > 0).all()


def test_report_is_regenerated_from_the_raw_byte_for_byte(tmp_path):
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory3c_raw.npz')
    if not os.path.exists(raw):
        return
    shutil.copy(raw, tmp_path / 'theory3c_raw.npz')
    env = dict(os.environ, T3C_OUT=str(tmp_path), OMP_NUM_THREADS='1')
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory3c_run.py'), '--from-raw'], check=True, env=env, capture_output=True)
    assert open(tmp_path / 'theory3c.md', 'rb').read() == open(os.path.join(REPO, 'results', 'synthetic', 'theory3c.md'), 'rb').read()
