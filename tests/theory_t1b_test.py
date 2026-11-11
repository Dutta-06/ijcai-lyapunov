"""Theorem T1 widened (experiments/theory_t1b_*.py): Delta_T1(K) grows slowly with K, the planar certificate against the integrated planar network, the ball-arithmetic
verification against the floating-point certificate, the eps-uniform envelope of the dormant coefficient (Lemma 12.4), and the deterministic report.
Run explicitly:  pytest tests/theory_t1b_test.py   (about 1-2 min; the file name keeps it out of the default collection)."""
import os
import subprocess
import sys

import numpy as np

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_planar as PL     # noqa: E402
import theory_t1b_rigorous as RG   # noqa: E402


def test_delta_t1_grows_slowly_with_K_and_R():
    d = {K: [T.delta_t1(K, R, 0.4167, 0.02) for R in (1.0, 2.4)] for K in (2, 4, 8)}
    assert d[2][0] <= d[4][0] <= d[8][0] and d[2][1] <= d[4][1] <= d[8][1]
    assert d[8][0] - d[2][0] < 1.5 and d[8][1] - d[2][1] < 1.5            # logarithmic in K, not linear
    assert all(d[K][1] > d[K][0] for K in d)                              # larger amplitude ratio needs a larger separation
    assert abs(d[2][0] - 4.9) < 0.15 and abs(d[2][1] - 5.8) < 0.15       # the K = 2 values of theory_t1.md


def test_dormant_coefficient_envelope_is_eps_uniform_lemma_12_4():
    t = np.linspace(1e-4, np.pi / 2 - 1e-9, 20001)
    assert (np.sin(t) ** 3 / (t ** 3 * np.cos(t))).min() >= 1 - 1e-12     # Adamovic-Mitrinovic: (sin t / t)^3 >= cos t
    for eps in (1e-4, 1e-8, 1e-12):
        U = T.Unit(eps)
        s = np.geomspace(30.0, 0.999 * U.lam, 800)
        env = 4 / (T.LAW['G'] * T.LAW['gamma'] ** 2 * T.LAW['m']) / (1 - (20.0 / s) ** 2) / s ** 3
        assert (U.B(s) / env).max() <= 1 + 1e-6


def test_planar_certificate_holds_and_the_integrated_planar_network_follows():
    cs, As, c0 = PL.make_target2(1, 0.0, 1.0, 1)
    prep = PL.prepare2(cs, As, c0)
    eps = 1e-7
    o, pi = PL.best_cert2(prep, eps)
    assert o['ok'] and pi >= 0 and o['D'] < 200
    sol, U = PL.simulate2(prep, eps, nlocal=200)
    assert PL.follows2(prep, eps, sol, U)[0]
    row = PL.diagnostics2(prep, eps, T.PARS[pi])
    assert row is not None


def test_ball_arithmetic_verification_agrees_with_the_float_certificate():
    cs, As, c0 = T.make_target(1, 0.0, 1.0, 1)
    prep = T.prepare(cs, As, c0)
    from flint import arb
    eps = arb(10.0 ** -8).union(arb(10.0 ** -8.125))
    o, pi = T.best_cert(prep, 1e-8)
    assert o['ok']
    ent = [e['kernel'] for e in prep['L']['events']]
    Dv = [1.0 if j in ent else max(1.5 * o['Etot'][j], 1.0) for j in range(len(c0))]
    ok, info = RG.verify(prep, eps, T.PARS[pi], Dv)
    assert ok                                                              # the certificate holds for every real eps in the ball
    non = [j for j in range(len(c0)) if j not in ent]
    rel = np.array([(info['Etot'][j] - o['Etot'][j]) / o['Etot'][j] for j in non])
    assert rel.min() > -0.01                                               # the ball bound is not below the float value (up to its own 1 % margins)
    jm = max(non, key=lambda j: o['Etot'][j])
    assert abs(float(info['Etot'][jm]) / o['Etot'][jm] - 1) < 0.15         # and, enclosing a ball of eps of width 0.125 decade, stays within 15 % for the kernel that carries the largest deviation


def _agreement(K, sep, R, seed, log_eps):
    """ball verification at a point eps against the float certificate: (ok, rel. difference of the largest kernel, smallest rel. difference over the dormant kernels)"""
    from flint import arb
    cs, As, c0 = T.make_target(K, sep, R, seed)
    prep = T.prepare(cs, As, c0)
    eps = 10.0 ** log_eps
    o, pi = T.best_cert(prep, eps)
    assert o['ok']
    ent = [e['kernel'] for e in prep['L']['events']]
    non = [j for j in range(len(c0)) if j not in ent]
    Dv = [1.0 if j in ent else max(1.5 * o['Etot'][j], 1.0) for j in range(len(c0))]
    ok, info = RG.verify(prep, arb(eps), T.PARS[pi], Dv)
    jm = max(non, key=lambda j: o['Etot'][j])
    layer = {k: info[k] / o[k] - 1 for k in ('T1', 'T2', 'T3')}                # the layer constants of the ball version against the float ones
    return ok, info['Etot'][jm] / o['Etot'][jm] - 1, min((info['Etot'][j] - o['Etot'][j]) / o['Etot'][j] for j in non), layer


def test_ball_verification_agrees_with_the_float_certificate_for_K2_and_K3():
    """referee round T1b, finding 6: the agreement (about 0.1 %) of the ball bound with the float certificate, for the kernel carrying the largest deviation, for K = 2 and K = 3"""
    for tgt in ((2, 6.0, 1.0, 11, -7.5), (3, 8.0, 1.0, 31, -8.5)):
        ok, rel_max, rel_min, layer = _agreement(*tgt)
        assert ok
        assert layer['T1'] > -0.002 and layer['T1'] < 0.002 and abs(layer['T2']) < 1e-5 and abs(layer['T3']) < 1e-5     # the relaxation times T_2, T_3 are not dropped or clamped wrongly
        assert abs(rel_max) < 0.005                                        # within 0.5 %
        assert rel_min > -0.01                                             # no kernel's ball bound below the float value (beyond its own margins)


def _k1():
    cs, As, c0 = T.make_target(1, 0.0, 1.0, 1)
    prep = T.prepare(cs, As, c0)
    o, pi = T.best_cert(prep, 1e-8)
    return prep, o, T.PARS[pi]


def test_verifier_asserts_the_hypotheses_r_gt_d0_critical_point_in_Dhat_theta_lt_vsharp():
    """referee round T1b, findings 2: each hypothesis is checked inside `verify`; a violated one makes the verification fail"""
    from flint import arb
    prep, o, par = _k1()
    eps = arb(1e-8)
    ent = [e['kernel'] for e in prep['L']['events']]
    Dv = [1.0 if j in ent else max(1.5 * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]
    ok, info = RG.verify(prep, eps, par, Dv)
    assert ok and all(info['cond'][k] for k in ('r>d0', 'crit in D', 'theta<v#'))
    # (i) an entrant offset above the radius
    p1 = dict(prep, d0=[par['r'] + 0.01] * len(prep['d0']))
    ok1, i1 = RG.verify(p1, eps, par, Dv)
    assert not ok1 and not i1['cond']['r>d0'] and i1['cond']['crit in D']
    # (ii) a block critical point outside D-hat: centre displaced beyond d-hat, then coefficient outside [v_-, v_+]
    for key, val in (('c', 5.0), ('w', 0.5)):
        ev = [dict(e) for e in prep['L']['events']]
        ev[-1][key] = ev[-1][key] + val if key == 'c' else ev[-1][key] * val
        p2 = dict(prep, L=dict(prep['L'], events=ev))
        ok2, i2 = RG.verify(p2, eps, par, Dv)
        assert not ok2 and not i2['cond']['crit in D']
    # (iii) theta not below v_#
    Bx = RG.rig_box(prep['K'], prep['R'], 0.0, par['r'], par['Om'], par['dfrac'], par['vsig'])
    h = RG.hypotheses(prep, par, Bx, RG.A(par['r']), Bx['dh'], RG.A(prep['amin']) * 0.001)
    assert not h['theta<v#'] and h['r>d0'] and h['crit in D']
    h = RG.hypotheses(prep, par, Bx, RG.A(par['r']), Bx['dh'], RG.A(prep['amin']) * 0.4)
    assert h['theta<v#']
    h = RG.hypotheses(prep, par, Bx, RG.A(max(prep['d0'])), Bx['dh'], RG.A(prep['amin']) * 0.4)         # r = d_0 is not enough (strict)
    assert not h['r>d0']
    a_ = RG.A(prep['amin'])
    assert RG.hypotheses(prep, par, Bx, RG.A(par['r']), Bx['dh'], a_ * 0.05)['theta<v#']                  # theta = 0.01 a_min < 0.05 a_min
    assert not RG.hypotheses(prep, par, Bx, RG.A(par['r']), Bx['dh'], a_ * 0.005)['theta<v#']
    assert not RG.hypotheses(prep, par, Bx, RG.A(par['r']), Bx['dh'], RG.A(0.005 * prep['amin']).union(RG.A(0.05 * prep['amin'])))['theta<v#']   # undecided: not certain
    h = RG.hypotheses(prep, par, Bx, RG.A(par['r']), RG.A(-0.1), RG.A(prep['amin']) * 0.4)               # a (negative) d-hat below the critical point's offset
    assert not h['crit in D']


def test_hypotheses_are_checked_for_every_entrant_and_every_coefficient_K2():
    """K = 2: r above every d_0 (not only the first or the last), the critical point of the final block in D-hat coordinate by coordinate, both bounds of the coefficients"""
    cs, As, c0 = T.make_target(2, 6.0, 1.0, 11)
    prep = T.prepare(cs, As, c0)
    o, pi = T.best_cert(prep, 1e-7)
    par = T.PARS[pi]
    r = par['r']
    Bx = RG.rig_box(prep['K'], prep['R'], prep['Delta'], r, par['Om'], par['dfrac'], par['vsig'])
    assert Bx is not None
    A = RG.A
    ok = RG.hypotheses(prep, par, Bx, A(r), Bx['dh'], A(prep['amin']) * 0.4)
    assert ok['r>d0'] and ok['crit in D'] and ok['theta<v#']
    for d0 in ([r + 0.01, 0.05], [0.05, r + 0.01], [r, 0.05]):
        assert not RG.hypotheses(dict(prep, d0=d0), par, Bx, A(r), Bx['dh'], A(prep['amin']) * 0.4)['r>d0']
    assert RG.hypotheses(dict(prep, d0=[0.05, r - 0.01]), par, Bx, A(r), Bx['dh'], A(prep['amin']) * 0.4)['r>d0']
    dh, wlo, whi = float(Bx['dh'].lower()), float((A(prep['amin']) * Bx['vD_lo']).upper()), float((A(prep['amin']) * Bx['vD_hi']).lower())
    for i in range(2):
        for key, f in (('c', lambda x: x + dh + 0.01), ('c', lambda x: x - dh - 0.01), ('w', lambda x: np.sign(x) * (wlo - 0.01)), ('w', lambda x: np.sign(x) * (whi + 0.01))):
            ev = [dict(e) for e in prep['L']['events']]
            v = ev[-1][key].copy()
            v[i] = f(v[i]) if key == 'w' else f(float(prep['cs'][ev[i]['feature']]))
            ev[-1][key] = v
            h = RG.hypotheses(dict(prep, L=dict(prep['L'], events=ev)), par, Bx, A(r), Bx['dh'], A(prep['amin']) * 0.4)
            assert not h['crit in D'], (i, key)

    # balls that are not points: the certain end of each ball decides
    rb = A(r - 0.01).union(A(r + 0.01))
    assert not RG.hypotheses(dict(prep, d0=[0.05, r]), par, Bx, rb, Bx['dh'], A(prep['amin']) * 0.4)['r>d0']        # lower end of r is below d_0
    dhw = A(dh).union(A(dh + 0.05))
    ev = [dict(e) for e in prep['L']['events']]
    v = ev[-1]['c'].copy()
    v[0] = float(prep['cs'][ev[0]['feature']]) + dh + 0.02                                                         # inside the upper end of d-hat, outside the lower end
    ev[-1]['c'] = v
    assert not RG.hypotheses(dict(prep, L=dict(prep['L'], events=ev)), par, Bx, A(r), dhw, A(prep['amin']) * 0.4)['crit in D']
    Bw = dict(Bx, vD_lo=A(float(Bx['vD_lo'].mid())).union(A(float(Bx['vD_lo'].mid()) + 0.3)))
    ev = [dict(e) for e in prep['L']['events']]
    v = ev[-1]['w'].copy()
    v[0] = np.sign(v[0]) * prep['amin'] * (float(Bx['vD_lo'].mid()) + 0.1)                                         # above the lower end of v_-, below its upper end
    ev[-1]['w'] = v
    assert not RG.hypotheses(dict(prep, L=dict(prep['L'], events=ev)), par, Bw, A(r), Bx['dh'], A(prep['amin']) * 0.4)['crit in D']


def test_bootstrap_inequality_is_strict_E_lt_D_over_1_0001():
    """referee round T1b, finding 2: the fixed point E_j(D) < D_j/1.0001; with D_j = f E_j(D) the test fails for f = 1.00005 and holds for f = 1.001"""
    from flint import arb
    prep, o, par = _k1()
    eps = arb(1e-8)
    ent = [e['kernel'] for e in prep['L']['events']]
    non = [j for j in range(len(prep['c0'])) if j not in ent]

    def at(f):
        Dv = [1.0 if j in ent else max(1.5 * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]
        for _ in range(8):                                                  # D = f E(D): E depends weakly on D, the iteration contracts
            ok, info = RG.verify(prep, eps, par, Dv)
            Dv = [1.0 if j in ent else f * info['Etot'][j] for j in range(len(prep['c0']))]
        return RG.verify(prep, eps, par, Dv)
    for f in (1.0005, 1.00015):                                             # E/D = 0.9995 and 0.99985: below 1/1.0001
        ok, info = at(f)
        assert ok and info['cond']['fixed point E<=D']
    ok, info = at(1.00005)                                                  # E/D = 0.99995: E <= D holds, E < D/1.0001 does not
    assert not ok and not info['cond']['fixed point E<=D'] and info['worst_ratio'] < 1.0


def test_undecided_balls_are_bounded_above_not_replaced_by_zero():
    """referee round T1b, finding 5: guards on a ball that straddles the threshold use its upper end"""
    from flint import arb
    straddle = arb(0.5).union(arb(2.0))
    assert float(RG.logpos_up(straddle).upper()) >= np.log(2.0) - 1e-12     # was 0
    assert float(RG.logpos_up(arb(0.5)).upper()) == 0.0 and float(RG.logpos_up(arb(1.0)).upper()) == 0.0
    assert float(RG.logpos_up(arb(5.0)).upper()) >= np.log(5.0) - 1e-12
    assert float(RG.pos_up(arb(-1.0).union(arb(3.0))).upper()) >= 3.0       # was 0
    assert float(RG.pos_up(arb(-2.0)).upper()) == 0.0 and float(RG.pos_up(arb(4.0)).upper()) >= 4.0


def test_report_is_regenerated_from_the_raw_byte_for_byte(tmp_path):
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory_t1b_raw.npz')
    if not os.path.exists(raw):
        return
    import shutil
    shutil.copy(raw, tmp_path / 'theory_t1b_raw.npz')
    env = dict(os.environ, TT1B_OUT=str(tmp_path), OMP_NUM_THREADS='1')
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory_t1b_wide.py'), '--from-raw'], check=True, env=env, capture_output=True)
    assert open(tmp_path / 'theory_t1b.md', 'rb').read() == open(os.path.join(REPO, 'results', 'synthetic', 'theory_t1b.md'), 'rb').read()


def test_reverification_report_is_regenerated_from_the_raw_byte_for_byte(tmp_path):
    import shutil
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory_t1b_reverify_raw.npz')
    if not os.path.exists(raw):
        return
    shutil.copy(raw, tmp_path / 'theory_t1b_reverify_raw.npz')
    env = dict(os.environ, TT1B_OUT=str(tmp_path), OMP_NUM_THREADS='1')
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory_t1b_reverify.py'), '--from-raw'], check=True, env=env, capture_output=True)
    assert open(tmp_path / 'theory_t1b_reverify.md', 'rb').read() == open(os.path.join(REPO, 'results', 'synthetic', 'theory_t1b_reverify.md'), 'rb').read()
    assert '8 of 8 reverified covers coincide' in open(tmp_path / 'theory_t1b_reverify.md').read()
