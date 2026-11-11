"""Round theory3, direction D: the readouts of depth p (experiments/theory3d_unit.py, theory3d_scan.py, theory3d_readouts.py): the exact scalar unit against the coordinate ODE and the repository's network
flows, ball arithmetic, depth 1 and depth 2, the sliding unit's first integral and overtaking margin, the stagewise path records, the scan helpers, and that the committed report is the report of the committed raw.
Run explicitly:  pytest tests/theory3d_test.py   (about 1 min; the file name keeps it out of the default collection)."""
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pytest

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import theory3d_unit as U            # noqa: E402
import theory3d_scan as SC           # noqa: E402
import theory3d_readouts as R        # noqa: E402
from src.core import factorized_rbf as FR      # noqa: E402
import exp_parametrization as P      # noqa: E402

warnings.filterwarnings('ignore')
OUT = os.environ.get('T3D_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD, SCAN_RAW = os.path.join(OUT, 'theory3d_raw.npz'), os.path.join(OUT, 'theory3d.md'), os.path.join(OUT, 'theory3d_scan_raw.npz')


def same(a, b, rtol=1e-9, atol=1e-300):
    """deep comparison of stored and recomputed results (dicts, lists, arrays, floats); atol absorbs the integrators' platform-dependent residuals"""
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(same(a[k], b[k], rtol, atol) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(same(x, y, rtol, atol) for x, y in zip(a, b))
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return False
        return bool(np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)) if a.dtype.kind in 'fc' else bool(np.array_equal(a, b))
    if isinstance(a, (float, np.floating)) or isinstance(b, (float, np.floating)):
        return bool((np.isnan(a) and np.isnan(b)) or abs(a - b) <= rtol * max(abs(a), abs(b)) + atol)
    return bool(a == b)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the closed forms
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def test_walls_are_the_repositorys():
    for eps in (1e-2, 1e-6, 1e-10):
        assert abs(U.lam_p(3, eps) - FR.Cubic3(1.0).scale(eps)) < 1e-12 * U.lam_p(3, eps)
        for N in (3, 4, 5):
            assert abs(U.lam_p(N, eps, ell=0.5) - FR.PowerN(N, 0.5).scale(eps)) < 1e-12 * U.lam_p(N, eps, ell=0.5)
        assert abs(U.lam_2(eps) - FR.Product2(1.0).scale(eps)) < 1e-14
    L = P.law()
    ro = FR.PaperLaw(L)
    for eps in (1e-4, 1e-8):
        pu = U.PaperUnit(eps)
        assert abs(pu.lam / ro.scale(eps) - np.sqrt(1 - eps / (2 * pu.vs))) < 1e-12      # lam = (transit scale) sqrt(1 - eps/(2 v_s)): the wall is c_u/sqrt(eps) (1 + O(eps))
    assert abs(U.paper_c_u() - 7.0248) < 1e-4 and abs(U.paper_B_const() - 33.3333) < 1e-3


def test_unit_closed_form_against_ode_with_sign_changing_drive():
    A = R.part_A()
    for r in A['rows']:
        assert r['err_eta'] < 1e-9 and r['err_w'] < 1e-8, r              # p = 0 is the paper's unit: err_w is the drift of Q
        assert 0.4 < r['eta_max_over_lam'] < 0.95                       # the drive really goes up to the wall's neighbourhood
    ps = {r['p'] for r in A['rows']}
    assert ps == {0, 2, 3, 4, 5, 6}


def test_wp_is_odd_increasing_and_its_inverse_is_exact():
    for p in (3, 4, 5):
        for eps in (1e-2, 1e-6):
            lam = U.lam_p(p, eps)
            eta = np.array([1e-12, 1e-6, 0.01, 0.5, 0.9, 0.999, 1 - 1e-9]) * lam
            w = U.W_p(eta, p, eps)
            assert np.all(np.diff(w) > 0) and np.allclose(U.W_p(-eta, p, eps), -w, rtol=1e-12, atol=0)
            assert np.max(np.abs(U.eta_of_w_p(w, p, eps) / eta - 1)) < 1e-9
            assert np.max(np.abs(U.eta_of_w_p(-w, p, eps) / eta + 1)) < 1e-9
    e = np.array([1e-6, 0.1, 1.0, 2.0])
    assert np.allclose(U.eta_of_w_2(U.W_2(e, 1e-6), 1e-6), e, rtol=1e-12)


def test_distance_form_and_dormant_bound_are_exact_and_sharp():
    for p in (3, 4, 5):
        q = p / (p - 2.0)
        for eps in (1e-4, 1e-8):
            for mu in (0.5, 0.1, 0.01):
                lam = U.lam_p(p, eps)
                w1 = float(U.W_p((1 - mu) * lam, p, eps))
                assert abs(w1 / float(U.dormant_bound_p(mu, p, eps)) - 1) < 1e-12
                assert abs(float(U.w_of_distance_p(mu * lam, p, eps)) / w1 - 1) < 1e-12
                simple = float(U.dormant_bound_simple_p(mu, p, eps))
                assert simple >= w1 and abs(simple - (eps ** (p / 2.0)) * mu ** (-q)) < 1e-15 * simple
                assert float(U.dormant_bound_p(mu, p, eps)) > 0.99 * simple * (1 - (mu / (2 - mu)) ** q)
    # the exponents themselves: eps^(p/2) and mu^(-p/(p-2))
    assert abs(float(U.dormant_bound_simple_p(0.1, 3, 1e-6)) - 1e-9 * 1e3) < 1e-15
    assert abs(float(U.dormant_bound_simple_p(0.1, 4, 1e-6)) - 1e-12 * 1e2) < 1e-18
    assert abs(float(U.dormant_bound_simple_p(0.001, 5, 1e-4)) / (1e-10 * 1000 ** (5 / 3.0)) - 1) < 1e-12


def test_formation_distance_brackets_and_margin_exponents():
    for p in (3, 4, 5):
        for eps in (1e-4, 1e-6, 1e-10):
            for th in (1e-2, 1e-3):
                dth = U.d_theta_p(th, p, eps)
                d0 = U.d_theta_p_leading(th, p)
                assert dth <= d0 * (1 + 1e-12) and abs(float(U.w_of_distance_p(dth, p, eps)) / th - 1) < 1e-9
                lam = U.lam_p(p, eps)
                assert abs(U.mu_theta_p(th, p, eps) / (th ** (-(p - 2.0) / p) * eps ** ((p - 2) / 2.0)) - 1) < 1e-4
                x0 = d0 / lam
                d_lo = d0 * (1 - (x0 / (2 - x0)) ** (p / (p - 2.0))) ** ((p - 2.0) / p)
                assert d_lo <= dth * (1 + 1e-12) and dth <= d0 * (1 + 1e-12)                         # the bracket d_lo <= d_theta < d_0 of the note
    # the margin at which a depth-p unit is formed scales as eps^((p-2)/2)
    m3 = U.mu_theta_p(0.01, 3, 1e-6) / U.mu_theta_p(0.01, 3, 1e-8)
    m4 = U.mu_theta_p(0.01, 4, 1e-6) / U.mu_theta_p(0.01, 4, 1e-8)
    m5 = U.mu_theta_p(0.01, 5, 1e-6) / U.mu_theta_p(0.01, 5, 1e-8)
    assert abs(m3 - 10) < 1e-3 and abs(m4 - 100) < 1e-2 and abs(m5 - 1000) < 1e-1
    # the inverse near the wall is accurate in the distance (the cancellation branch): d = 1e-6 at a wall of 333
    for dd in (1e-6, 1e-4, 1e-2):
        w = float(U.w_of_distance_p(dd, 3, 1e-6))
        eta = float(U.eta_of_w_p(w, 3, 1e-6)[0])
        assert abs((U.lam_p(3, 1e-6) - eta) / dd - 1) < 1e-6
    # ell cancels
    assert abs(U.mu_theta_p(0.01, 3, 1e-6, ell=7.0) / U.mu_theta_p(0.01, 3, 1e-6) - 1) < 1e-9


def test_depth_two_is_exactly_logarithmic():
    eps, th = 1e-6, 0.01
    f = lambda mu: float(U.dormant_bound_2(mu, eps))
    assert abs(f(0.3) - (eps ** 0.3 - eps ** 1.7)) < 1e-18
    mt = float(U.mu_theta_2(th, eps))
    assert abs(f(mt) / th - 1) < 1e-9                                        # the root of f(mu) = theta
    mu0 = float(U.rule_depth2(th, eps))
    assert abs(mu0 - 1.0 / 3) < 1e-12 and mt < mu0 and mu0 - mt < 1e-7        # the rule is an upper bound, essentially exact
    d = float(U.d_theta_2(th, eps))
    assert abs(float(U.w_of_distance_2(d, eps)) / th - 1) < 1e-9
    # inverse and dual: w = 2 eps sinh(4 eta)
    assert abs(float(U.W_2(0.75, eps)) - 2 * eps * np.sinh(3.0)) < 1e-18
    # the amplitude-dependent wall (and the mass: the wall of b = 1, lam_2(eps), is the unit of the relative dual)
    assert abs(float(U.wall_shift_2(1.0, 1e-6, b=3.0)) - (1.0 - np.log(3.0) / np.log(1e6))) < 1e-12
    assert abs(float(U.wall_shift_2(1.0, 1e-6)) - 1.0) < 1e-12 and abs(float(U.wall_shift_2(1e-3, 1e-6)) - 0.5) < 1e-12
    for A_ in (0.1, 0.6, 3.0):                                                # the dual at which |w| reaches A is lam_2 (1 + ln(A/b)/ln(1/eps)) up to (eps b/A)^2
        assert abs(float(U.wall_shift_2(A_, 1e-6, b=2.0)) * float(U.lam_2(1e-6)) - float(U.eta_of_w_2(A_, 1e-6, b=2.0))) < 1e-9 * float(U.lam_2(1e-6))
    # a unit with margin 0.0664 (the design) and theta = 0.009 needs ln(1/eps) > ln(1/theta)/mu = 70.9
    assert float(U.mu_theta_2(0.009, 1e-30)) > 0.0664 > float(U.mu_theta_2(0.009, 1e-32))


def test_ball_arithmetic_checks_all_certified():
    res = U.ball_checks()
    assert set(res) == {'B1', 'B2', 'B3', 'B4', 'B5', 'B6'}
    for k, (n, bad) in res.items():
        assert n > 100 and bad == 0, (k, n, bad)
    # a planted false statement is rejected: the formation bracket with the wrong exponent fails some case
    from flint import arb
    assert not (arb(2) ** arb(-3) > arb(2) ** arb(-2))


def test_dual_is_the_integral_of_the_drive_in_the_repositorys_flows():
    for name, eps in (('D3', 1e-4), ('D2', 1e-4), ('P0', 1e-4)):
        r = R.run_with_dual(name, eps, tau_end=14.0, n_out=29)
        err = np.abs(r['eta_int'] - r['eta_inv']) / r['lam']
        assert err.max() < 5e-7, (name, err.max())
        assert r['status'] == 0
    r = R.run_with_dual('D3', 1e-4, tau_end=14.0, n_out=29)
    # the flow has reached the first entries: some coefficient is O(1)
    assert np.abs(r['W']).max() > 0.05


def test_parameters_mass_identity_and_paper_unit_forms():
    from numpy.polynomial.legendre import leggauss
    xg, wg = leggauss(80)
    for p in (3, 4, 5):
        eps, b, ell = 1e-4, 1.3, 0.8
        lam = U.lam_p(p, eps, b, ell)
        s = np.sqrt(eps * b)
        for sgm in (0.1, 0.6, 0.9):
            eta = sgm * lam
            u, v = U.uv_p(eta, p, eps, b, ell)
            assert abs((u ** p - v ** p) / U.W_p(eta, p, eps, b, ell) - 1) < 1e-12
            assert abs(u ** (-(p - 2.0)) / (s ** (-(p - 2.0)) * (1 - sgm)) - 1) < 1e-12 and abs(v ** (-(p - 2.0)) / (s ** (-(p - 2.0)) * (1 + sgm)) - 1) < 1e-12
            # Psi_p is the integral of W_p
            z = 0.5 * eta * (xg + 1)
            quad = 0.5 * eta * float((wg * U.W_p(z, p, eps, b, ell)).sum())
            assert abs(quad / float(U.Psi_p(eta, p, eps, b, ell)) - 1) < 1e-9
        mu = 0.05
        assert float(U.mass_bound_p(mu, p, eps, b, ell)) <= eps * b / (2 * ell * p) * mu ** (-2.0 / (p - 2)) and float(U.mass_bound_p(mu, p, eps, b, ell)) > 0.5 * eps * b / (2 * ell * p) * mu ** (-2.0 / (p - 2))
    e2 = np.linspace(0.0, 0.6, 8)
    z = 0.5 * 0.6 * (xg + 1)
    assert abs(0.5 * 0.6 * float((wg * U.W_2(z, 1e-3, 1.2, 0.7)).sum()) / float(U.Psi_2(0.6, 1e-3, 1.2, 0.7)) - 1) < 1e-9
    # the paper's unit
    for eps in (1e-4, 1e-8):
        pu = U.PaperUnit(eps)
        for x in (0.3, 1.0, 2.5):
            w = float(pu.w_of_x(x))
            assert abs(pu.x_of_w(w) / x - 1) < 1e-10
            assert abs(pu.eta_of_w(w) / float(pu.eta_of_x(x)) - 1) < 1e-10 and abs(pu.eta_of_w(-w) / float(pu.eta_of_x(x)) + 1) < 1e-10
            assert abs(pu.dist_of_w(w) - (pu.lam - float(pu.eta_of_x(x)))) < 1e-6 * pu.lam
        # the coefficient at the wall is 0.0079 (so theta = 0.009 is reached just beyond it) and the dormant bound is below Lemma 2.2's forms
        assert abs(pu.w_of_x(pu.x_of_w(0.0079)) - 0.0079) < 1e-12 and pu.mu_theta(0.009) < 0 < pu.mu_theta(0.001)
        assert abs(pu.mu_theta(0.001) - pu.dist_of_w(0.001) / pu.lam) < 1e-14 and 0 < pu.mu_theta(0.001) < 0.1
        for mu in (0.5, 0.1, 0.02):
            ex = pu.dormant_bound(mu)
            if 1.0 / (np.tan(np.pi * mu / 2) * pu.kap) >= 1:              # Lemma 2.2 needs s_eps < 1
                continue
            f1, f2 = pu.dormant_bound_T1(mu)
            assert 0 < ex <= f1 * (1 + 1e-12) <= f2 * (1 + 1e-12)
            lead = 8 / np.pi ** 3 * 1.0 / (12 * np.sqrt(2 * pu.vs)) * eps ** 1.5 * mu ** -3
            if eps == 1e-8 and mu <= 0.1:
                assert 0.8 * lead < ex < 1.3 * lead
    # the distance law of the paper's unit: |w| ~ 33.3/d^3
    pu = U.PaperUnit(1e-8)
    for dd in (400.0, 2000.0):
        w = float(pu.w_of_x(pu.x_of_w(U.paper_B_const() / dd ** 3)))
        assert abs(pu.dist_of_w(w) / dd - 1) < 0.1


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# depth 1
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def test_depth_one_fixed_centres_closed_form_and_certificate():
    D = R.part_D()
    assert D['err_fixed'] < 1e-8 and D['sigma_min'] > 1e-3
    assert D['n_cert'] >= 5 and all(D['cert_ok'])
    g0 = np.array(D['g0'])
    cert, tj, Tj = U.d1_formation_certificate(g0, D['normH'], P.LR_C, D['thr'])
    sim = np.array(D['tj_sim'])
    assert np.all(sim[cert] <= tj[cert] * (1 + 1e-9)) and np.all(tj[cert] <= Tj[cert] * (1 + 1e-12))
    # the leading prediction theta/(ell |g|) is within a factor 2 of the simulated formation time for the strongly driven units
    pred = np.array(D['tj_pred'])
    strong = np.abs(g0) > 0.5 * np.abs(g0).max()
    assert np.all((sim[strong] / pred[strong] < 2.5) & (sim[strong] / pred[strong] > 0.4))
    # a planted H: w(t) = H^-1 (I - exp(-ell H t)) b
    H = np.array([[2.0, 0.5], [0.5, 1.0]])
    b0 = np.array([0.3, -0.2])
    w, g = U.d1_fixed_centres(H, b0, 1.5, np.array([0.0, 1.0, 1e4]))
    assert np.allclose(w[0], 0) and np.allclose(g[0], b0) and np.allclose(w[-1], np.linalg.solve(H, b0), atol=1e-9) and np.allclose(g[-1], 0, atol=1e-9)
    from scipy.linalg import expm
    assert np.allclose(g[1], expm(-1.5 * H) @ b0, atol=1e-12)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the sliding unit
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def test_sliding_unit_first_integral_and_first_order_lead():
    E = R.part_E()
    for r in E['landscape']:
        assert r['drift'] < 1e-8 and 0.79 < r['dc_over_sig'] < 0.81 and r['G_ratio'] > 1.5
    lin = {r['h_rel']: r for r in E['linear']}
    for r in E['linear']:
        assert r['first_integral_err'] < 1e-8 and r['reached']
    assert abs(lin[0.3]['lead'] / lin[0.3]['pred'] - 1) < 5e-3 and abs(lin[1.0]['lead'] / lin[1.0]['pred'] - 1) < 5e-3
    assert abs(lin[3.0]['lead'] / lin[3.0]['pred'] - 1) < 0.05
    assert lin[10.0]['lead'] / lin[10.0]['pred'] > 1.3                           # the first order stops where the displacement is several widths
    assert lin[0.3]['lead'] < lin[1.0]['lead'] < lin[3.0]['lead']


def test_overtaking_margin_root_and_limits():
    for Kc in (1e-8, 1e-6, 1e-4):
        for xa in (0.05, 0.3, 1.0):
            x = float(U.overtake_margin(Kc, xa))
            assert abs(x - Kc * (1 / x - 1 / xa)) < 1e-12 * x
            assert x <= np.sqrt(Kc) * (1 + 1e-12)
        assert abs(float(U.overtake_margin(Kc, 1e9)) / np.sqrt(Kc) - 1) < 1e-6
    assert float(U.overtake_margin(1e-6, 0.5)) < float(U.overtake_margin(1e-4, 0.5))
    assert abs(U.lead_const('D3') - 1 / 6) < 1e-15 and abs(U.lead_const('P0') - 0.33774) < 1e-4
    assert abs(float(U.drift_lead(12.0, 0.02, 1e-6, 0.02, 1.0, 1 / 6.0)) - 12 / 0.02 ** 2 * 1e-6 / 6 * (50 - 1)) < 1e-12


def test_path_record_design_and_scores():
    L, K, c0 = P.setup()
    rec = R.path_record(L, K, c0, np.array(P.CSTAR), np.array(P.ASTAR))
    assert rec['G'].shape == (3, 21) and rec['exits'] == 0 and rec['failed'] == 0
    assert abs(rec['margin'] - 0.06638) < 5e-5 and abs(rec['margin_never'] - 0.06638) < 5e-5
    assert [e[1] for e in rec['ev']] == [5, 16, 11]
    pairs = R.stage_pairs(rec)
    x0, x1, q, ids = pairs
    assert x1.min() > 0 and np.all(x0 >= x1) and np.all(q >= 0) and ids.shape[1] == 2
    # the entrant of a stage is never its own competitor; the active units are not competitors
    for (k, j) in ids:
        assert j not in [e[1] for e in rec['ev'][:k + 1]]
    # the score increases with eps; the critical eps is where it crosses 1
    s = [R.overtake_score(pairs, e, R.A_RO['D3']) for e in (1e-10, 1e-8, 1e-6, 1e-4, 1e-2)]
    assert s == sorted(s) and s[0] < 1 and s[3] < 1 < s[-1]                  # stored runs: D3 forms the path's set at 1e-4, not at 1e-2
    ec = R.critical_eps(pairs, R.A_RO['D3'])
    assert 1e-4 < ec < 1e-2 and abs(R.overtake_score(pairs, ec, R.A_RO['D3']) - 1) < 1e-6
    assert R.rho_stat(pairs) > 0
    assert R.overtake_score(pairs, 1e-6, R.A_RO['D3']) < 0.5


def test_stage_pairs_on_a_planted_record():
    rec = dict(ev=[(1.0, 0, 1), (2.0, 1, 1)], taus=np.array([0.0, 1.0, 2.0]), G=np.array([[0.1, 0.05, 0.02], [0.0, 0.03, 0.02]]), H=np.array([[0.2, 0.1, 0.4], [0.0, 0.1, 0.4]]),
               U0=np.array([[0.0, 0.0, 0.0], [1.0, 0.4, 0.5]]), U1=np.array([[1.0, 0.5, 0.2], [1.0, 1.0, 0.8]]))
    x0, x1, q, ids = R.stage_pairs(rec, mu_c=10.0)
    # stage 0: entrant 0 excluded; competitors 1, 2. stage 1: active {0}, entrant 1 excluded; competitor 2
    assert [tuple(i) for i in ids] == [(0, 1), (0, 2), (1, 2)]
    assert np.allclose(x1, [0.5, 0.8, 0.2]) and np.allclose(x0, [1.0, 1.0, 0.5])
    assert np.allclose(q, [10 * 0.01 / 0.0025, 10 * 0.16 / 0.0004, 10 * 0.16 / 0.0004])
    # rho: the smallest x1/sqrt(q); the score and the critical eps on this record
    qq = np.array([10 * 0.01 / 0.0025, 10 * 0.16 / 0.0004, 10 * 0.16 / 0.0004])
    assert abs(R.rho_stat((x0, x1, q, ids)) - float(np.min(np.array([0.5, 0.8, 0.2]) / np.sqrt(qq)))) < 1e-15
    assert np.isnan(R.rho_stat((x0[:0], x1[:0], q[:0], ids[:0]))) and R.overtake_score((x0[:0], x1[:0], q[:0], ids[:0]), 1e-6, 1 / 6.0) == 0.0
    # inconsistent records exercise the explicit masks: an active unit whose dual is not pinned, an entrant whose dual is not pinned, a competitor at its wall, a competitor with no drive
    rec2 = dict(ev=[(1.0, 0, 1)], taus=np.array([0.0, 1.0]), G=np.array([[0.1, 0.05, 0.02, 0.0]]), H=np.array([[0.2, 0.1, 0.4, 0.1]]),
                U0=np.array([[0.0, 0.0, 0.0, 0.0]]), U1=np.array([[0.7, 0.5, 1.0, 0.3]]))
    x0, x1, q, ids = R.stage_pairs(rec2, mu_c=10.0)
    assert [tuple(i) for i in ids] == [(0, 1)]                                 # entrant 0 (U1 = 0.7), the unit at its wall (2) and the unit without drive (3) are not competitors
    rec3 = dict(ev=[(1.0, 0, 1), (2.0, 1, 1)], taus=np.array([0.0, 1.0, 2.0]), G=np.array([[0.1, 0.05], [0.1, 0.05]]), H=np.array([[0.2, 0.1], [0.2, 0.1]]),
                U0=np.array([[0.0, 0.0], [1.0, 0.4]]), U1=np.array([[1.0, 0.5], [0.6, 1.0]]))
    x0, x1, q, ids = R.stage_pairs(rec3, mu_c=10.0)
    assert [tuple(i) for i in ids] == [(0, 1)]                                 # stage 1: unit 0 is active (its dual 0.6 is not pinned in this record) and unit 1 is the entrant


def test_observed_critical_brackets():
    mk = lambda se, oe: dict(set_equal=se, order_equal=oe)
    cells = {1e-4: mk(False, False), 1e-5: mk(False, False), 1e-6: mk(True, False), 1e-7: mk(True, True), 1e-8: mk(True, True)}
    assert SC_critical(cells) == (1e-6, 1e-7)
    assert SC_critical({1e-4: mk(True, True), 1e-5: mk(True, True)}) == (None, 1e-5)
    assert SC_critical({1e-4: mk(False, True), 1e-5: mk(False, True)}) == (1e-5, None)
    assert R.confusion([1, 1, 0, 0], [1, 0, 1, 0]) == dict(tp=1, fp=1, fn=1, tn=1)
    assert R.formed_theta(np.array([1.0, -0.5])) == 0.01 and R.formed_theta(np.array([1.0, 0.01])) == 0.001


def SC_critical(cells):
    return R.observed_critical(cells)


def test_local_exponent_and_matching():
    assert abs(R.local_exponent(1e-6, 1e-8, 10.0, 1.0) - 0.5) < 1e-12
    assert abs(R.local_exponent(1e-6, 1e-8, 0.0185, 0.0057) - 0.256) < 5e-3                # the paper's 0.26


def test_fixed_centre_limit_rates():
    cells = R.L_CELLS
    R.L_CELLS = (('D3', (1e-4, 1e-8)), ('D2', (1e-6, 1e-24)), ('D1', (1.0,)))
    try:
        out = R.part_L()
    finally:
        R.L_CELLS = cells
    rows = {(r['name'], r['eps']): r for r in out['rows']}
    ratio = rows[('D3', 1e-4)]['sup_dual'] / rows[('D3', 1e-8)]['sup_dual']
    assert 80 < ratio < 150                                                   # sqrt(eps): 100 for four decades
    assert 1.0 < rows[('D3', 1e-8)]['sup_dual'] / np.sqrt(1e-8) < 1.5
    r2 = rows[('D2', 1e-6)]['sup_dual'] / rows[('D2', 1e-24)]['sup_dual']
    assert 2.5 < r2 < 12                                                      # logarithmic: 18 decades of eps gain a factor of a few, not 1e9
    assert rows[('D1', 1.0)]['l2'] > 0.5 and rows[('D1', 1.0)]['n_formed'] >= 15 and rows[('D3', 1e-8)]['n_formed'] == 3
    assert [j for t, j in out['entries']] == [5, 16, 11] and abs(out['entries'][0][0] - 10.7217) < 1e-3
    if os.path.exists(RAW):
        stored = {(r['name'], r['eps']): r for r in R.load(RAW)['L']['rows']}
        for k, r in rows.items():
            assert same(r, stored[k], 1e-6, 1e-8)
        assert out['tau'] == R.L_TAU


def test_tie_family_has_a_sqrt_law_at_fixed_geometry():
    old = R.K_FAMILY
    R.K_FAMILY = {'D3': dict(deltas=(0.0006, 0.005), lo=-10.0, hi=-2.0, nit=8)}
    try:
        K = R.part_K()
    finally:
        R.K_FAMILY = old
    rows = K['D3']
    assert all(r['flip_hi'] and not r['flip_lo'] for r in rows)
    le = np.array([0.5 * (r['log_eps_lo'] + r['log_eps_hi']) for r in rows])
    lm = np.log10([r['margin'] for r in rows])
    slope = (lm[1] - lm[0]) / (le[1] - le[0])
    assert 0.45 < slope < 0.6
    ratio = [r['pred'] / 10 ** (0.5 * (r['log_eps_lo'] + r['log_eps_hi'])) for r in rows]
    assert all(1.5 < x < 3.5 for x in ratio)
    if os.path.exists(RAW):
        stored = {r['delta']: r for r in R.load(RAW)['K']['D3']}
        for r in rows:
            assert same(r, stored[r['delta']], 1e-6, 1e-8)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the stored raws and the report
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def test_scan_helpers_and_targets():
    T = SC.load_targets(R.PARAM_RAW)
    assert len(T) == 200 and T[0]['i'] == 0 and abs(T[0]['ref']['margin'] - 0.00938) < 1e-5
    assert SC.SPEC['D3']['eps'][0] == 1e-4 and abs(SC.SPEC['D3']['eps'][-1] - 1e-11) < 1e-24 and len(SC.SPEC['D3']['eps']) == 15
    assert SC.SPEC['P0']['eps'] == (1e-5, 1e-6, 1e-7, 1e-8, 1e-9)
    n3 = sum(t['ref']['margin'] < SC.SPEC['D3']['margin_max'] for t in T)
    n0 = sum(t['ref']['margin'] < SC.SPEC['P0']['margin_max'] for t in T)
    assert (n3, n0) == (118, 55)


def test_scan_cell_reproduces_a_stored_cell():
    T = SC.load_targets(R.PARAM_RAW)
    t = T[9]
    c = SC.run_cell('D3', t, 1e-6, P.setup())
    stored = [x for x in t['cells'] if x['name'] == 'D3' and x['eps'] == 1e-6][0]
    assert c['set_equal'] == stored['set_equal'] and c['n_formed'] == stored['n_formed'] and c['form_order'] == stored['form_order']


def test_scan_is_resumable_and_records_the_cells():
    old = SC.SPEC['D3']
    SC.SPEC['D3'] = dict(eps=(1e-4, 1e-8), margin_max=0.0015)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'scan.npz')
        try:
            res = SC.scan(('D3',), limit=1, path=path, log=False)
            tab = res['D3']
            assert sorted(tab) == [(9, 1e-8), (9, 1e-4)] and set(tab[(9, 1e-4)]) == set(SC.KEEP)
            first = dict(tab)
            again = SC.scan(('D3',), limit=1, path=path, log=False)['D3']
            assert all(again[k] is not None and again[k]['secs'] == first[k]['secs'] for k in first)           # nothing was run again
            assert SC.load_scan(path)['D3'][(9, 1e-4)]['secs'] == first[(9, 1e-4)]['secs']
        finally:
            SC.SPEC['D3'] = old
    assert not os.path.exists(os.path.join(tmp, 'scan.npz.tmp.npz'))


@pytest.mark.skipif(not os.path.exists(RAW), reason='no stored raw yet')
def test_stored_parts_are_consistent():
    res = R.load(RAW)
    assert set(res) >= {'A', 'B', 'C', 'D', 'E', 'L', 'K', 'P', 'F', 'G', 'H'}
    assert all(v[1] == 0 for v in res['B']['counts'].values())
    for r in res['C']['rows']:
        assert r['err_max'] < 5e-8 and abs(r['ratio_max'] - 1) < 1e-2 and r['ratio_max'] <= 1 + 1e-6
    recs = res['P']['recs']
    assert len(recs) == 1034
    F = res['F']
    assert F['n_matched'] == 200
    by = {(c['name'], c['eps']): c for c in F['cells']}
    # recompute one confusion table from the stored path records and the stored runs
    T = SC.load_targets(R.PARAM_RAW)
    mi = R.match_population(recs, T)
    sc = np.array([R.overtake_score(R.stage_pairs(recs[m]), 1e-6, R.A_RO['D3']) for m in mi])
    anym = np.array([(not R.cells_of(t, 'D3', 1e-6)['set_equal']) or (not R.cells_of(t, 'D3', 1e-6)['order_equal']) for t in T])
    assert R.confusion(sc >= 1, anym) == by[('D3', 1e-6)]['conf_drift_any']
    # the paper's counts of mismatches are reproduced: 34 / 9 for P0 and 24 / 5 for D3
    assert (by[('P0', 1e-6)]['n_set_mismatch'], by[('P0', 1e-8)]['n_set_mismatch'], by[('D3', 1e-6)]['n_set_mismatch'], by[('D3', 1e-8)]['n_set_mismatch']) == (34, 9, 24, 5)
    assert abs(by[('P0', 1e-6)]['max_margin_obs'] - 0.0185) < 1e-4 and abs(by[('P0', 1e-8)]['max_margin_obs'] - 0.0057) < 1e-4


@pytest.mark.skipif(not os.path.exists(RAW), reason='no stored raw yet')
def test_stored_analyses_are_reproduced_from_their_inputs():
    res = R.load(RAW)
    scan = SC.load_scan(SCAN_RAW)
    assert same(R.part_A(), res['A'], 1e-6, 1e-8) and same(R.part_B(), res['B']) and same(R.part_D(), res['D'], 1e-6, 1e-8) and same(R.part_E(), res['E'], 1e-6, 1e-8)
    F, G, H = R.part_F(res['P']), R.part_G(res['P'], scan), R.part_H(res['P'])
    assert same(F, res['F']) and same(G, res['G']) and same(H, res['H'])
    # the scan covers the targets it says it covers, with the grids it says
    assert {i for (i, _) in scan['D3']} >= {i for i, t in enumerate(SC.load_targets(R.PARAM_RAW)) if t['ref']['margin'] < SC.SPEC['D3']['margin_max']}
    assert {i for (i, _) in scan['P0']} == {i for i, t in enumerate(SC.load_targets(R.PARAM_RAW)) if t['ref']['margin'] < SC.SPEC['P0']['margin_max']}
    assert all(len([k for k in scan['P0'] if k[0] == i]) == 5 for i in {i for (i, _) in scan['P0']})
    assert [r['delta'] for r in res['K']['D3']] == list(R.K_FAMILY['D3']['deltas']) and [r['delta'] for r in res['K']['P0']] == list(R.K_FAMILY['P0']['deltas'])
    assert res['L']['tau'] == R.L_TAU and res['L']['n_events'] == 34


@pytest.mark.skipif(not os.path.exists(RAW), reason='no stored raw yet')
def test_part_C_is_reproduced_and_the_runs_have_the_stated_counts():
    res = R.load(RAW)
    C = R.part_C()
    assert same(C, res['C'], 1e-6, 1e-8)
    counts = {r['name']: r['n_formed'] for r in C['rows']}
    assert counts == {'D2': 9, 'D3': 3, 'D4': 3, 'D5': 3, 'P0': 3}
    assert all(r['err_max'] < 5e-8 and abs(r['ratio_max'] - 1) < 1e-3 for r in C['rows'])


@pytest.mark.skipif(not (os.path.exists(RAW) and os.path.exists(MD)), reason='no stored report yet')
def test_report_is_the_report_of_the_raw():
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(RAW, tmp)
        if os.path.exists(SCAN_RAW):
            shutil.copy(SCAN_RAW, tmp)
        env = dict(os.environ, T3D_OUT=tmp, OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory3d_readouts.py'), '--from-raw'], env=env, check=True, capture_output=True)
        assert open(os.path.join(tmp, 'theory3d.md')).read() == open(MD).read()
