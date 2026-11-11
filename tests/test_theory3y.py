"""The entering orbit of the coupled entry layer and the separated regime (experiments/theory3y_entry.py, theory3y_sep.py, theory3y_run.py): the unit's closed forms and the tail lemma, the layer data and box constants
(golden values), the explicit bounds, the integrated layer flows against the theorems' checks, the aligned / orbit distances, the Hessian bounds of the relaxation set, the closed-form separated-regime condition against the
exact recursion, and the report's regeneration."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3y_entry as Y      # noqa: E402
import theory3y_sep as S        # noqa: E402
import theory3y_run as R        # noqa: E402

TARGET = (2, 6.0, 1.0, 11)


@pytest.fixture(scope='module')
def layer2():
    L = Y.layers_of(TARGET)[1]
    return L, Y.layer_box(L)


# ---- the unit ---------------------------------------------------------------------------------------------------------------------
def test_unit_inverse_and_mobility_closed_form():
    for eps in (0.0, 1e-9, 1e-4):
        for w in (1e-9, 3e-4, 0.3, 0.9):
            assert abs(float(Y.W_of_x(Y.x_of_w(w, eps), eps)) - w) < 1e-12 * max(w, 1e-3)
        assert Y.x_of_w(-0.3, eps) == pytest.approx(-Y.x_of_w(0.3, eps))
    # M = dW/d eta with d eta/dx = a_max/(G q) (THEORY_T1 Lemma 2.1), by finite differences
    eps = 1e-6
    x = np.array([0.3, 1.0, 2.5])
    h = 1e-6
    dW = (Y.W_of_x(x + h, eps) - Y.W_of_x(x - h, eps)) / (2 * h)
    q = eps + (2 * Y.VS - eps) * np.tanh(x) ** 2
    assert np.allclose(dW / (Y.LAW['a_max'] / (Y.LAW['G'] * q)), Y.M_of_x(x, eps), rtol=1e-6)
    assert Y.M_of_x(0.0, eps) == pytest.approx(Y.LAW['G'] * eps ** 2 / Y.LAW['m'])
    assert Y.M_of_x(0.0, 0.0) == 0.0


def test_tail_lemma_constant_and_inequality():
    wh = 0.3359212445442466
    xh = Y.x_of_w(wh, 0.0)
    cM = (1.0 * 12.0 / 1.0) * (2 * Y.VS * 1.0 / (12.0 * np.cosh(xh) ** 4)) ** (2.0 / 3.0)
    assert Y.c_tail(wh) == pytest.approx(cM, rel=1e-12)
    assert Y.c_tail(wh) == pytest.approx(0.000870562, rel=1e-5)
    worst = np.inf
    for eps in (0.0, 1e-12, 1e-6, 1e-3):
        for w in np.geomspace(1e-12, wh, 60):
            worst = min(worst, Y.M_of_w(w, eps) / (Y.c_tail(wh) * w ** (4.0 / 3.0)))
    assert worst > 20                                        # the lemma holds with a margin of about two orders


def test_M_lower_is_a_certain_lower_bound():
    for eps in (0.0, 1e-9, 1e-6):
        for w in (0.05, 0.3, 0.9):
            lo, ex = Y.M_lower(w, eps), Y.M_of_w(w, eps)
            assert lo <= ex and lo > 0.999999 * ex


# ---- layers and the box --------------------------------------------------------------------------------------------------------
def test_layers_of_golden(layer2):
    L, _ = layer2
    Ls = Y.layers_of(TARGET)
    assert [l['k'] for l in Ls] == [0, 1] and len(Ls[0]['wb']) == 0 and len(L['wb']) == 1
    assert Ls[0]['d0'] == pytest.approx(0.018984645911, abs=1e-9) and L['d0'] == pytest.approx(0.376728262852, abs=1e-9)
    assert L['wb'][0] == pytest.approx(0.7000334, abs=1e-6) and L['cb'][0] == pytest.approx(3.15416033, abs=1e-6)
    assert np.allclose(L['what'], [0.7, 0.7]) and np.allclose(L['chat'], [3.15446152, -3.15446152], atol=1e-7)
    assert L['Delta'] == pytest.approx(6.308923, abs=1e-6) and L['amin'] == pytest.approx(0.7) and L['sign_ok']
    assert np.allclose(L['sgn'], [1, 1])


def test_layer_box_golden_and_regime(layer2):
    L, b = layer2
    assert b['r'] == pytest.approx(0.4) and b['r'] > L['d0']
    assert b['vlo'] == pytest.approx(0.95977498, abs=1e-6) and b['vbar'] == pytest.approx(1.00101445, abs=1e-6)
    assert b['Ev'] == pytest.approx(0.00101445, abs=1e-7) and b['Ec'] == pytest.approx(0.00279428, abs=1e-7) and b['dtilde'] == pytest.approx(0.00581663, abs=1e-7)
    assert b['dhat'] == pytest.approx(0.0255258, abs=1e-6) and b['nu'] == pytest.approx(0.47885890, abs=1e-6) and b['Hup'] == pytest.approx(1.02774996, abs=1e-6)
    assert b['dtilde'] < b['dhat'] < b['r']
    # K = 1: no cross terms, v_bar = R exactly, d-tilde = 0
    L1 = Y.layers_of((1, 0.0, 1.0, 1))[0]
    b1 = Y.layer_box(L1)
    assert b1['Ev'] == 0.0 and b1['vbar'] == 1.0 and b1['dtilde'] == 0.0 and b1['r'] > L1['d0']
    # outside the regime: a close pair with a large amplitude ratio has no admissible box
    outs = [Y.layer_box(l) is None for l in Y.stress_layers()]
    assert any(outs) and not all(outs)                       # the stress set straddles the threshold


def test_stress_layers_construction():
    st = Y.stress_layers()
    assert Y.STRESS_D0 == 0.40 and Y.STRESS_DELTAS == (4.9, 5.2, 5.6) and Y.STRESS_RS == (1.0, 2.4)
    assert len(st) == 16 and len({l['target'] for l in st}) == 16
    for l in st:
        assert l['K'] == 2 and l['d0'] == Y.STRESS_D0 and abs(l['ce'] - (l['cs'][1] - Y.STRESS_D0)) < 1e-15 and l['amin'] == 0.7
        assert np.allclose(np.abs(l['what']), np.abs(l['As']))
        assert abs(float(l['wb'][0]) - float(l['As'][0])) < 0.05 * 0.7         # the block is the first feature up to the overlap
    assert any(np.sign(l['As'][0]) != np.sign(l['As'][1]) for l in st)


# ---- the explicit bounds ----------------------------------------------------------------------------------------------------
def test_bounds_golden_and_formulas(layer2):
    L, b = layer2
    bd = Y.bounds(L, b, 0.0)
    assert bd['w_hash'] == pytest.approx(0.5 * 0.7 * b['vlo'], rel=1e-12)
    assert bd['T2'] == pytest.approx(0.0311773, abs=1e-6) and bd['T3'] == pytest.approx(117.5697, abs=1e-3) and bd['That'] == pytest.approx(bd['T2'] + bd['T3'], rel=1e-12)
    assert bd['rho'] == pytest.approx(0.0331450, abs=1e-6) and bd['pref'] == pytest.approx(np.sqrt(b['Hup'] / b['nu']), rel=1e-12)
    # T3 by hand: the lower-face branch with M(w_#)
    Vt = np.exp(-b['dhat'] ** 2 / 4) - b['Ev']
    t3a = np.log((0.7 * Vt - bd['w_hash']) / (0.7 * 0.01)) / bd['mV']
    assert bd['T3'] == pytest.approx(max(t3a, np.log(max(0.7 * (b['vbar'] - 1 - b['Ev']) / (0.7 * 0.01), 1.0)) / bd['mVp']), rel=1e-12)
    assert bd['mV'] == pytest.approx(Y.M_of_w(bd['w_hash'], 0.0), rel=1e-6) and bd['mV'] <= Y.M_of_w(bd['w_hash'], 0.0)
    assert bd['kd'] == pytest.approx(Y.MC * bd['w_hash'] * 0.7 * np.exp(-b['r'] ** 2 / 4) / 2, rel=1e-9)
    assert bd['rho'] == pytest.approx(min(bd['m_lo'], 1.0) * b['nu'], rel=1e-12) and bd['m_lo'] <= Y.MC * 0.49
    # eps > 0 changes only the mobilities
    bde = Y.bounds(L, b, 1e-6)
    assert bde['T2'] == bd['T2'] and abs(bde['T3'] - bd['T3']) < 1e-6 and bde['w_hash'] == bd['w_hash']


# ---- the flows against the theorems' checks -------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def flows(layer2):
    L, b = layer2
    out = {}
    for key, (eps, v0) in {'e6': (1e-6, None), 'z6': (0.0, 1e-6), 'z9': (0.0, 1e-9), 'ref': (0.0, 1e-12)}.items():
        bd = Y.bounds(L, b, eps)
        tr = Y.integrate_layer(L, eps, v0=v0, w_hash=bd['w_hash'], s2=1500.0, dt=0.5)
        out[key] = (eps, bd, tr)
    return L, b, out


def test_flow_checks_all_pass_and_bounds_hold(flows):
    L, b, out = flows
    for key, (eps, bd, tr) in out.items():
        ok, vec = Y.check_layer(L, b, bd, tr, eps)
        assert all(ok.values()), (key, ok)
        assert 0 < vec['t_hat'] < bd['That'] and 0.2 < vec['ratio_hat'] < 0.9
        assert vec['ratio_decay'] < 1 and vec['err_end'] < 1e-8
        assert vec['t_hat'] == pytest.approx(65.5, abs=0.5) and vec['t_star'] == pytest.approx(99.5, abs=0.5) and vec['t_star'] < vec['t_star_bd'] == pytest.approx(260.41, abs=0.05)
        assert vec['margin_star'] == pytest.approx(0.00101445, abs=1e-7) and vec['err_hat'] == pytest.approx(0.01094, abs=1e-4)
        assert vec['dmax_over_r'] < 1 and vec['vmax_over_vbar'] <= 1 + 1e-9
    # the entrant crosses the level at s = 0: |w_e| = w_# there, the flow from the dormant state needs ~eps^(-1/2) to get there
    eps, bd, tr = out['e6']
    i = int(np.argmin(np.abs(tr['s'])))
    assert abs(abs(float(Y.W_of_x(tr['x'][-1][i], eps))) - bd['w_hash']) < 1e-9
    assert 5e3 < tr['t_hash'] < 2e4
    # t_#(v_0) ~ v_0^(-1/3): a factor 10 for three decades
    r = out['z9'][2]['t_hash'] / out['z6'][2]['t_hash']
    assert 8.0 < r < 12.0


def test_aligned_flows_converge_to_one_orbit(flows):
    L, b, out = flows
    ref = out['ref'][2]
    d6 = Y.aligned_distance(out['z6'][2], ref, L, 0.0, 0.0, -0.9 * out['z6'][2]['t_hash'], 1500.0)
    d9 = Y.aligned_distance(out['z9'][2], ref, L, 0.0, 0.0, -0.9 * out['z9'][2]['t_hash'], 1500.0)
    o6 = Y.orbit_distance(out['z6'][2], ref, L, 0.0, 0.0, 1e-6)
    o9 = Y.orbit_distance(out['z9'][2], ref, L, 0.0, 0.0, 1e-9)
    assert d9 < 0.05 * d6 and o9 < 0.05 * o6                  # Cauchy: three decades of v_0 give two orders (v_0^(2/3))
    assert d9 < 5e-4 and o9 < 5e-4
    dn = Y.aligned_distance(out['e6'][2], ref, L, 1e-6, 0.0, -0.9 * out['e6'][2]['t_hash'], 1500.0)
    assert dn < 5e-4                                          # the network's layer flow is close to the limit orbit
    assert Y.aligned_distance(ref, ref, L, 0.0, 0.0, -1000.0, 1000.0) == pytest.approx(0.0, abs=1e-12)


def test_orbit_distance_is_invariant_under_a_time_shift(flows):
    L, b, out = flows
    ref = out['ref'][2]
    sh = dict(ref)
    sh['s'] = ref['s'] + 0.7
    assert Y.aligned_distance(sh, ref, L, 0.0, 0.0, -300.0, 300.0) > 1e-4
    assert Y.orbit_distance(sh, ref, L, 0.0, 0.0, 1e-9) < 1e-9


def test_check_layer_detects_violations(flows):
    L, b, out = flows
    eps, bd, tr = out['z6']
    bad = dict(b)
    bad['r'] = 0.3                                            # smaller than the entrant's offset 0.377
    assert not Y.check_layer(L, bad, bd, tr, eps)[0]['box_d']
    bad = dict(b)
    bad['vbar'] = 0.9
    assert not Y.check_layer(L, bad, bd, tr, eps)[0]['box_vhi']
    bad = dict(b)
    bad['vlo'] = 1.5
    assert not Y.check_layer(L, bad, bd, tr, eps)[0]['box_vlo_block']
    bad_bd = dict(bd)
    bad_bd['rho'] = 5.0                                       # a decay rate far too large for the flow
    assert not Y.check_layer(L, b, bad_bd, tr, eps)[0]['decay']
    bad_bd = dict(bd)
    bad_bd['w_hash'] = 0.9 * bd['w_hash']                      # a different level: the tail check uses the c_M of that level
    assert isinstance(Y.check_layer(L, b, bad_bd, tr, eps)[0]['tail'], bool)
    tr2 = dict(tr)
    xe = tr['x'].copy()
    xe[-1][3000] = 0.5 * xe[-1][3000]                          # a dip of the entrant before it reaches v_lo
    tr2['x'] = xe
    assert not Y.check_layer(L, b, bd, tr2, eps)[0]['mono_entrant']
    # the entrant falls below v_lo after having reached it
    xe = tr['x'].copy()
    xe[-1][-200:] = Y.x_of_w(0.5 * 0.7 * b['vlo'], 0.0)
    tr3 = dict(tr)
    tr3['x'] = xe
    assert not Y.check_layer(L, b, bd, tr3, eps)[0]['entrant_ge_vlo_after']
    # a block unit leaves the relaxation set
    xb = tr['x'].copy()
    xb[0][-50:] = Y.x_of_w(0.7 * 0.9, 0.0)
    tr4 = dict(tr)
    tr4['x'] = xb
    assert not Y.check_layer(L, b, bd, tr4, eps)[0]['block_in_hat']


# ---- the Hessian bounds on the relaxation set (the constants nu and H_up of Theorem E3) -----------------------------------------
def _Fprime(v, d, L):
    """F' = F / a_min^2 as a function of the units' normalised coefficients v and offsets d (continuum, Z = 1)."""
    amin = L['amin']
    feats = L['feats']
    cs, As = L['cs'], L['As']
    c = cs[feats] + d
    w = L['sgn'] * amin * v
    GG = lambda x: np.exp(-np.asarray(x) ** 2 / 4)
    ff = float(As @ GG(cs[:, None] - cs[None, :]) @ As)
    cross = float(w @ (GG(c[:, None] - cs[None, :]) @ As))
    return (0.5 * ff - cross + 0.5 * float(w @ GG(c[:, None] - c[None, :]) @ w)) / amin ** 2


def test_hessian_bounds_on_the_relaxation_set(layer2):
    L, b = layer2
    rng = np.random.default_rng(3)
    n = 2
    dh, R, Ev, vs = b['dhat'], L['R'], b['Ev'], b['varsigma']
    lo, hi = np.exp(-dh ** 2 / 4) - Ev - vs, R + Ev + vs
    assert lo == pytest.approx(b['vD_lo'], abs=1e-12)
    eig_lo, eig_hi = np.inf, -np.inf
    for _ in range(150):
        y = np.r_[rng.uniform(lo, hi, n), rng.uniform(-dh, dh, n)]
        h = 1e-4

        def f(z):
            return _Fprime(z[:n], z[n:], L)
        H = np.zeros((2 * n, 2 * n))
        for i in range(2 * n):
            for j in range(2 * n):
                e_i, e_j = np.eye(2 * n)[i] * h, np.eye(2 * n)[j] * h
                H[i, j] = (f(y + e_i + e_j) - f(y + e_i - e_j) - f(y - e_i + e_j) + f(y - e_i - e_j)) / (4 * h * h)
        ev = np.linalg.eigvalsh(0.5 * (H + H.T))
        eig_lo, eig_hi = min(eig_lo, ev[0]), max(eig_hi, ev[-1])
    assert eig_lo >= b['nu'] - 1e-6 and eig_hi <= b['Hup'] + 1e-6
    assert eig_lo > 0.4 and eig_hi < 1.05                      # the bounds are not vacuous here


# ---- the separated regime ----------------------------------------------------------------------------------------------------------
def test_lemma_phi_and_region():
    assert S.lemma_phi(6.0, 0.8333) and not S.lemma_phi(3.0, 0.8333) and not S.lemma_phi(6.0, 1.5)
    D = np.linspace(0.8333, 3.0, 400)
    assert np.all(np.diff(S.phi(D, 6.0)) <= 1e-15)             # nonincreasing on [h, Delta/2]
    c = S.conditions(6.309, 0.8333, 0.019, 0.3767)
    assert c['ok'] and c['s1'] and c['s2'] and c['s3']
    assert c['P_minus'] == pytest.approx(0.99985913, abs=1e-7) and c['L1'] == pytest.approx(0.96529719, abs=1e-7)
    assert c['q2'] == pytest.approx(0.96514274, abs=1e-7) and c['T_star'] == pytest.approx(0.03622949, abs=1e-7) and c['T_comp'] == pytest.approx(0.05316059, abs=1e-7)
    assert not S.conditions(6.309, 0.8333, 0.3767, 0.019)['ok']          # offsets not ordered
    assert not S.conditions(6.309, 0.8333, 0.019, 0.4166667)['ok']       # a_2 = h/2: two kernels tie
    assert not S.conditions(6.309, 0.8333, 0.2, 0.2)['ok']               # a_1 = a_2: the two entry times tie
    assert not S.conditions(3.0, 0.8333, 0.019, 0.3767)['ok']            # Delta below the box threshold
    s = [np.mean([S.conditions(D, 0.8333, a1, a2)['ok'] for a1 in np.linspace(0.01, 0.4, 25) for a2 in np.linspace(0.01, 0.4, 25) if a1 < a2]) for D in (4.7, 5.6, 8.0)]
    assert s[0] < s[1] < s[2] == 1.0


def test_rho_bar_bounds_the_block_residual():
    b = S.box_for(5.6, 0.8333)
    dt, Ev = b['dtilde'], b['Ev']
    Q = 1 - float(S.g(dt)) + Ev
    rng = np.random.default_rng(5)
    for _ in range(2000):
        x = rng.uniform(dt + 0.01, 9.0) * rng.choice([-1, 1])
        d = rng.uniform(-dt, dt)
        v = rng.uniform(1 - Q, 1 + Q)
        assert abs(float(S.g(x) - v * S.g(x - d))) <= float(S.rho_bar(abs(x), dt, Q)) + 1e-14
    assert S.rho_bar(3.0, dt, Q) < S.rho_bar(2.0, dt, Q)


def test_condition_never_holds_when_the_recursion_fails():
    rng = np.random.default_rng(9)
    n_c = n_o = 0
    for Delta in (4.9, 5.6):
        for _ in range(25):
            a1, a2 = np.sort(rng.uniform(0, 0.4167, 2))
            sg = rng.choice([-1.0, 1.0], 2)
            cs, As, c0, sep = S.build_target(Delta, 0.8333, a1, a2, sg[0], sg[1], float(rng.choice([-1.0, 1.0])), bool(rng.integers(0, 2)))
            c = S.conditions(sep, 0.8333, a1, a2)
            o = S.oracle(cs, As, c0, a1, a2)
            n_c += int(c['ok'])
            n_o += int(o)
            assert (not c['ok']) or o
    assert n_c > 10 and n_o >= n_c
    # the exact tie: equal offsets -> the oracle (unique entrants) is False
    cs, As, c0, sep = S.build_target(6.0, 0.8333, 0.2, 0.2, 1.0, -1.0, 1.0, False)       # mirror-symmetric: offsets +0.2 and -0.2
    assert np.allclose(cs - cs.mean(), [-sep / 2, sep / 2]) and not S.oracle(cs, As, c0, 0.2, 0.2)


# ---- the report ------------------------------------------------------------------------------------------------------------------------
def test_report_regenerates_byte_for_byte(tmp_path, monkeypatch):
    monkeypatch.setattr(R, 'OUT', str(tmp_path))
    monkeypatch.setattr(R, 'RAW', os.path.join(str(tmp_path), 'theory3y_raw.npz'))
    monkeypatch.setattr(R, 'MD', os.path.join(str(tmp_path), 'theory3y.md'))
    raw = R.run(quick=True, parts='A,B,C,D,S')
    R.save_raw(raw)
    md1 = R.report(raw)
    md2 = R.report(R.load_raw())
    assert md1 == md2 and md1.startswith('# The entering orbit')
    assert 'Layers in the regime: **2 of 2**' in md1 and 'Violations of any check or bound over 2 layers x 2 values of eps: **0**' in md1
    assert 'Total violations: **0**' in md1
    # golden rows of the quick report
    assert '| 1e-06 | 2 | 2 | 0.568 | 0.7 | 0.0841 | 6.28e-16 | 1.05e+04 | 1.05e+04 | 0.95 |' in md1
    assert '| (2,4.9,1,101)/2 | S | 1.04e+06 | 43.2 | 122 | 0.354 | 0.307 | 46.8 | 556 | 0.00874 | 0.0545 |' in md1
    assert '| 0.001 | 2 | 2 | 0.568 | 0.7 | 0.0259 | 0.0247 | 0.399 | 0.389 | 103 | 105 |' in md1
    assert '| 1e-06 | 9.16e-05 | 8.92e-05 | 9.13e-05 | 8.89e-05 |' in md1
    assert '| 0.8333 | 4.7 | 1 | 0.2253 | 435 |' in md1 and '| 0.8333 | 5.6 | 30 | 27 | 30 | 27 | **0** |' in md1
    ac = raw['A_num']
    assert ac.shape == (2, len(R.A_COLS)) and np.all(raw['B_ok'][:, :, :]) and np.all(raw['C_ok'])


# ---- further data and a layer whose entrant has the negative sign -----------------------------------------------------------------------
def test_random_targets_golden():
    rt = Y.random_targets(8)
    assert rt[:3] == [(3, 9.03, 1.35, 7001), (3, 5.62, 1.05, 7002), (3, 9.17, 1.98, 7003)]
    assert rt[7] == (2, 6.17, 1.43, 7008) and all(5.5 <= t[1] <= 9.5 and 1.0 <= t[2] <= 2.4 for t in rt) and len(Y.random_targets(2)) == 2
    assert len(Y.D_TARGETS) == 8 and Y.D_TARGETS[3] == (2, 6.0, 1.0, 11)


def test_layer_with_a_negative_entrant():
    L = [l for l in Y.stress_layers() if l['target'][3] == 11 and l['target'][1] == 5.2][0]      # strong first? code 11: opposite signs, strong feature first
    assert L['sgn'][1] == -1 and L['sgn'][0] == 1
    b = Y.layer_box(L)
    assert b is not None
    for eps, v0 in ((0.0, 1e-6), (1e-6, None)):
        bd = Y.bounds(L, b, eps)
        tr = Y.integrate_layer(L, eps, v0=v0, w_hash=bd['w_hash'], s2=1500.0, dt=0.5)
        assert tr['x'][-1][0] < 0 or tr['x'][-1][0] == 0                  # the entrant's dual is negative: x_e < 0 along the flow
        assert np.all(tr['x'][-1][1:] <= 0)
        ok, vec = Y.check_layer(L, b, bd, tr, eps)
        assert all(ok.values()) and vec['ratio_hat'] < 1 and vec['ratio_decay'] < 1
        assert tr['x'][-1][-1] == pytest.approx(Y.x_of_w(-0.7, eps), rel=1e-6)


def test_layer_features_and_signs_follow_the_entry_order():
    Ls = Y.layers_of((3, 5.62, 1.05, 7002))
    assert [l['feats'] for l in Ls] == [[1], [1, 2], [1, 2, 0]]
    assert [list(l['sgn']) for l in Ls] == [[1.0], [1.0, -1.0], [1.0, -1.0, -1.0]]
    assert np.allclose(Ls[2]['a'], [0.735, 0.71728655, 0.7], atol=1e-8) and len(Ls[2]['wb']) == 2 and len(Ls[2]['cb']) == 2
    assert np.allclose(np.abs(Ls[2]['what']), [0.735, 0.71728655, 0.7], atol=1e-6)
    L2 = Y.layers_of((2, 6.0, 2.4, 12))[1]
    assert list(L2['feats']) == [0, 1] and list(L2['sgn']) == [-1.0, -1.0]


# ---- the separated regime: box threshold, target construction, the benchmark rows, the redundancy of the competitor classes --------------------------------------
def test_sep_box_threshold_in_delta():
    h = 0.8333
    assert S.box_for(4.2, h) is None and S.box_for(4.5, h) is None and S.box_for(4.6, h) is None       # (C) fails at 4.2, (N) from 4.3 to 4.6
    b = S.box_for(4.7, h)
    assert b is not None and b['nu'] == pytest.approx(0.0058, abs=5e-4) and b['r'] == pytest.approx(h / 2 + 0.02)
    assert S.box_for(5.6, h)['dtilde'] < S.box_for(4.7, h)['dtilde'] and S.box_for(8.0, h)['Ev'] < 1e-5


def test_sep_build_target_geometry():
    cs, As, c0, sep = S.build_target(6.0, 0.8333, 0.1, 0.3, 1.0, -1.0, -1.0, False)
    assert np.allclose(As, [0.7, -0.7]) and cs[0] == pytest.approx(0.1) and cs[1] == pytest.approx(7 * 0.8333 - 0.3) and sep == pytest.approx(cs[1] - cs[0])
    cs2, As2, c02, sep2 = S.build_target(6.0, 0.8333, 0.1, 0.3, 1.0, -1.0, 1.0, True)
    assert np.allclose(As2, [0.7, 0.7]) and cs2[0] == pytest.approx(-0.3) and cs2[1] == pytest.approx(7 * 0.8333 + 0.1)       # flip: the feature of offset a2 is the left one
    assert c02[0] <= cs2.min() - 8 + 0.8333 and c02[-1] >= cs2.max() + 7 and np.allclose(np.diff(c02), 0.8333)


def test_sep_benchmark_rows_and_the_competitor_classes():
    raw = S.run_part(quick=True)
    b = raw['S_bench']
    assert b.shape[1] == 12 and b[0][:2].tolist() == [6.0, 11.0]
    assert b[0][2] == pytest.approx(0.8333, abs=1e-3) and b[0][3] == pytest.approx(6.3089, abs=1e-3) and b[0][4] == pytest.approx(0.0190, abs=2e-4) and b[0][5] == pytest.approx(0.3767, abs=2e-4)
    assert b[0][6] == 1 and b[0][7] == 1 and b[1][6] == 0 and b[1][7] == 1               # the near-tie of seed 21: not covered, the recursion is exact
    assert raw['S_rows'][:, 6].sum() == 0 and raw['S_reg'].shape == (3, 5) and raw['S_rows'].shape == (3, 7)
    # the four competitor classes of (S3): C1 (the second-nearest kernel of the second feature) is the binding one everywhere on the sampled region
    h, binds = 0.8333, []
    for Delta in (4.9, 6.0, 8.0):
        box = S.box_for(Delta, h)
        dt, Ev = box['dtilde'], box['Ev']
        Q = 1 - float(S.g(dt)) + Ev
        for a1 in np.linspace(0.02, 0.38, 7):
            for a2 in np.linspace(0.03, 0.41, 7):
                if not a1 < a2:
                    continue
                Pm = float(S.g(a1) - S.g(Delta - a1))
                rho = lambda x: float(S.rho_bar(x, dt, Q))
                cl = [(S.phi(h - a2, Delta) / Pm, S.g(h - a2) + rho(Delta - h + a2)), (S.phi(h - a1, Delta) / Pm, S.g(Delta - h + a1) + rho(h - a1)),
                      (S.phi(h, Delta) / Pm, S.g(h) + rho(Delta / 2)), (S.phi(h, Delta) / Pm, S.g(Delta / 2) + rho(h))]
                T = [(1 - u) / r for (u, r) in cl]
                binds.append(int(np.argmin(T)))
                c = S.conditions(Delta, h, a1, a2, box)
                if c['s1'] and c['s2']:
                    assert c['T_comp'] == pytest.approx(min(T), rel=1e-12)
    assert set(binds) == {0}


# ---- the runner: layers outside the regime, the report tables on modified raw data ------------------------------------------------------------------
def test_run_skips_layers_outside_the_regime(monkeypatch):
    out = [l for l in Y.stress_layers() if Y.layer_box(l) is None][0]
    monkeypatch.setattr(R, 'all_layers', lambda quick: [(out, 'S')])
    raw = R.run(quick=True, parts='A,B,C,D')
    assert raw['A_num'][0, R.A_COLS.index('regime')] == 0 and np.isnan(raw['A_num'][0, R.A_COLS.index('r')])
    assert not raw['B_ok'].any() and not raw['C_ok'].any() and np.isnan(raw['B_vec']).all() and np.isnan(raw['C_dist']).all()
    md = R.report(raw)
    assert 'Layers in the regime: **0 of 1** (set D: 0 of 0, set R: 0 of 0, set S: 0 of 1)' in md and '## B.' not in md


@pytest.fixture(scope='module')
def quick_raw():
    return R.run(quick=True, parts='A,B,C,D,S')


def test_report_tables_follow_the_raw_data(quick_raw):
    raw = {k: np.array(v) for k, v in quick_raw.items()}
    md = R.report(raw)
    assert 'set D: 1 of 1, set R: 0 of 0, set S: 1 of 1' in md
    assert '| 1e-06 | 2 | 2 | 0.568 | 0.7 | 0.0841 | 6.28e-16 | 1.05e+04 | 1.05e+04 | 0.95 |' in md and '| 1e-06 | 2 | 2 | 0.568 | 0.7 | 0.0203 |' in md
    bi = {c: i for i, c in enumerate(R.BCOLS)}
    assert raw['B_bd'][0, 0, bi['mV']] == pytest.approx(0.03362575637110985, rel=1e-11) and abs(raw['B_bd'][0, 0, bi['mV']] - 0.0336257563219463) > 1e-11
    # tamper: a failing check in layer 0 at eps = 1e-6 and a smaller t_# in layer 1
    raw['B_ok'][0, 0, 2] = False
    vi = {c: i for i, c in enumerate(R.VEC)}
    raw['B_vec'][1, 0, vi['t_hash']] = 5000.0
    md2 = R.report(raw)
    assert 'Violations of any check or bound over 2 layers x 2 values of eps: **1**.' in md2
    assert '| 1e-06 | 2 | 1 | 0.568 | 0.7 | 0.0841 | 6.28e-16 | 5e+03 | 1.05e+04 | 0.95 |' in md2
    # the first column of the aligned/orbit distance table of part C and the D table keep their order
    assert md.index('max aligned distance | median aligned | max orbit distance') < md.index('| 0.001 | 2 | 2 |')
