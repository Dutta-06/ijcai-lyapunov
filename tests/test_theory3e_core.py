"""Fast core tests of experiments/theory3e_kernels.py (the criterion of the mutation run, tools/mutants/theory3e_kernels_own.py): the box formulas against a hand calculation with a stub Gram
(n = 1 and n >= 2, row-sum and block-norm cross-Hessian, dormant budget, centre-mobility condition number), the Gaussian regression against the paper's table, the directed bounds, the
cusp, the clean corollary, the families.  About 40 s."""
import math
import os
import sys

import numpy as np
import pytest

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import theory3e_kernels as TK      # noqa: E402

D0 = 0.4167


class Stub:
    """A Gram stand-in with fixed numbers (not a real Gram): every method the conditions call."""
    cusp = False

    def __init__(self, n=1, lamps_small=0.45, lamps_r=0.45):
        self.n = n
        self.ls = lamps_small
        self.lr = lamps_r

    def tail(self, o, l, hess_n=None):
        return {0: 0.01, 1: 0.02, 2: 0.03, 'H': 0.05}[o]

    def g_min(self, x, N=None):
        return 0.9 - 0.1 * x

    def psi_at(self, x):
        return 0.3 * x

    def lampsi(self, x, N=None):
        return self.lr if x > 0.2 else self.ls

    def psihat(self, x, N=None):
        return 0.4 * x

    def rho_min(self, x, N=None):
        return 0.45 - 0.05 * x

    def M1(self):
        return 0.43

    def M2(self, n=None):
        return 0.5


def hand(K, R, Delta, r, Om, hc, cond=1.0, n=1, ls=0.45):
    S0, S1 = (K - 1) * 0.01, (K - 1) * 0.02
    S2 = (K - 1) * (0.03 if n == 1 else 0.05)
    vbar = (R * (1 + S0) + Om) / (1 - S0)
    Ev = S0 * (R + vbar) + Om
    Ec = S1 * (R + vbar) + 0.43 * Om
    vlo = (0.9 - 0.1 * r) - Ev
    Cm = 0.3 * r - cond * Ec
    dt = Ec / 0.45
    vlo_d = (0.9 - 0.1 * dt) - Ev
    s = R * 0.4 * dt
    t = vlo_d * (0.45 - 0.05 * dt)
    lam = 0.5 * (1 + t - math.sqrt((1 - t) ** 2 + 4 * s * s))
    if n > 1:
        lam = min(lam, vlo_d * (0.45 if dt > 0.2 else ls))
    Rv = R + vbar
    if hc == 'rowsum':
        rowv = S0 + Rv * S1 + vbar * S1 + 0.43 * Om
        rowd = Rv * S1 + vbar * Rv * S2 + vbar * S1 + vbar * vbar * S2 + 0.43 * Om + vbar * 0.5 * Om
        Hc = max(rowv, rowd)
    else:
        n2 = lambda a, b, c: (a + c + math.sqrt((a - c) ** 2 + 4 * b * b)) / 2
        Hc = n2(0, Rv * S1 + 0.43 * Om, vbar * Rv * S2 + vbar * 0.5 * Om) + n2(S0, vbar * S1, vbar * vbar * S2)
    return dict(vbar=vbar, Ev=Ev, Ec=Ec, vlo=vlo, Cmargin=Cm, dtilde=dt, lam=lam, Hc=Hc, Nmargin=lam - Hc)


@pytest.mark.parametrize('hc', ['norm2', 'rowsum'])
@pytest.mark.parametrize('K,R,Om,cond', [(2, 2.0, 0.1, 1.0), (3, 1.0, 0.0, 1.0), (2, 1.5, 0.05, 2.0), (2, 1.5, 0.05, 1.2)])
def test_box_formulas_equal_the_hand_calculation(hc, K, R, Om, cond):
    b = TK.box(Stub(1), K, R, 10.0, 0.5, Omega=Om, hc=hc, cond_c=cond)
    h = hand(K, R, 10.0, 0.5, Om, hc, cond)
    for k, v in h.items():
        if k in b:
            assert abs(b[k] - v) < 1e-9 * max(1.0, abs(v)), (k, b[k], v)
    assert ('dtilde' in b) == (h['vlo'] > 0 and h['Cmargin'] > 0)
    assert b['ok'] == (h['vlo'] > 0 and h['Cmargin'] > 0 and h['dtilde'] < 0.5 and h['Nmargin'] > 0)


def test_tangential_stiffness_enters_for_n_at_least_two():
    b2 = TK.box(Stub(2, lamps_small=0.02), 2, 1.0, 10.0, 0.5, hess_n=2)
    h2 = hand(2, 1.0, 10.0, 0.5, 0.0, 'norm2', 1.0, n=2, ls=0.02)
    assert abs(b2['lam'] - h2['lam']) < 1e-9 and h2['lam'] < 0.05
    b1 = TK.box(Stub(2, lamps_small=0.02), 2, 1.0, 10.0, 0.5, hess_n=1)
    assert b1['lam'] > 0.2                                         # n = 1 has no tangential direction
    # the hessian-norm tail (n >= 2) enters through S_2
    assert b2['Hc'] > b1['Hc']


def test_failing_conditions_are_reported_not_ok():
    s = Stub(1)
    assert not TK.box(s, 2, 1.0, 10.0, 0.5, cond_c=1000.0)['ok']               # (C) fails with a huge condition number
    assert not TK.box(s, 2, 1.0, 10.0, 0.5, Omega=50.0)['ok']                  # (L) fails with a huge dormant budget
    assert not TK.box(s, 2, 1.0, 1.0, 0.5)['ok']                               # ell <= L_MIN
    assert TK.box(s, 2, 1.0, 10.0, 0.5)['ok']
    assert not TK.box(s, 2, 1.0, 10.0, 0.001)['ok']                            # (C) fails for a tiny radius
    assert not TK.box(s, 2, 60.0, 10.0, 0.5)['ok']                             # a huge amplitude ratio
    assert TK.box(s, 2, 1.0, 10.0, 0.5, sums=(0.01, 0.02, 0.03))['ok'] and not TK.box(s, 2, 1.0, 10.0, 0.5, sums=(1.2, 0.0, 0.0))['ok']
    assert TK.box(s, 1, 1.0, 10.0, 0.5)['ok']                                  # K = 1: no coupling
    # kappa >= 1 stops at once (without the guard the formulas return garbage that happens to pass (L), (C))
    assert 'Ev' not in TK.box(s, 2, 1.0, 10.0, 0.5, sums=(1.2, 0.0, 0.0)) and 'Ev' not in TK.box(s, 2, 1.0, 10.0, 0.5, sums=(1.0, 0.0, 0.0))
    # d-tilde < r is a condition of its own when psi/d is not monotone (lambda_psi(r) < psi(r)/r): here d-tilde = 0.615 > r = 0.5 although (C) and (L) hold
    st = Stub(1, lamps_r=0.1)
    b = TK.box(st, 2, 1.0, 10.0, 0.5, Omega=0.05)
    assert b['Cmargin'] > 0 and b['vlo'] > 0 and b['dtilde'] > 0.5 and not b['ok']
    assert TK.box(st, 2, 1.0, 10.0, 0.7, Omega=0.05)['dtilde'] < 0.7


def test_gaussian_regression_and_block_norm_values():
    g = TK.Gaussian(1)
    paper = {(2, 1.0): 4.66, (2, 2.4): 5.50, (2, 10.0): 6.70, (3, 1.0): 5.08, (3, 2.4): 5.84}
    for (K, R), v in paper.items():
        assert abs(TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=24.0, hc='rowsum') - v) < 1e-9
    blk = {(2, 1.0): 4.56, (2, 2.4): 5.44, (2, 10.0): 6.70, (3, 1.0): 5.00, (3, 2.4): 5.80}
    for (K, R), v in blk.items():
        assert abs(TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=24.0) - v) < 1e-9
    assert TK.entry_threshold(TK.Gaussian(3), 2, 1.0, D0, lo_=3.0, hi=24.0) == 4.56


def test_threshold_order_and_monotonicity():
    th = lambda k, K=2, R=1.0, d0=D0, **kw: TK.entry_threshold(k, K, R, d0, lo_=1.0, hi=60.0, **kw)
    g, e1, m52 = TK.family('gauss', 1), TK.family('exp', 1), TK.family('mat52', 1)
    assert th(g) < th(m52) < th(e1)
    assert th(g, K=3) > th(g) and th(g, R=2.4) > th(g) and th(g, d0=0.8) > th(g) and th(g, d0=0.0) < th(g)
    assert th(g, cond_c=2.0) > th(g) and th(g, cond_c=4.0) > th(g, cond_c=2.0)
    assert th(g, Omega=0.1) > th(g)
    t, r = TK.entry_threshold(g, 2, 1.0, D0, lo_=3.0, hi=24.0, return_r=True)
    assert r > D0 and TK.box(g, 2, 1.0, t, r)['ok'] and not TK.box(g, 2, 1.0, t - 0.1, r)['ok']


def test_families_and_tails():
    assert TK.family('exp', 1).name == 'matern1.5' and TK.family('exp', 2).name == 'matern2' and TK.family('exp', 3).name == 'matern2.5'
    assert TK.family('mat32', 1).name == 'matern3.5' and TK.family('mat52', 3).name == 'matern6.5'
    m = TK.MaternHalf(1, 1)
    assert [float(c) for c in m.P] == [1.0, 1.0] and m.lam0 == 1.0
    m = TK.MaternHalf(2, 1)
    assert abs(float(m.P[2]) - 1 / 3) < 1e-15 and abs(m.lam0 - 1 / 3) < 1e-15
    for k in (TK.family('gauss', 1), TK.family('mat32', 1), TK.family('exp', 2)):
        xs = np.linspace(0.5, 39.9, 20000)
        G, d1, d2 = np.abs(k.fG(xs)), np.abs(k.fdG(xs)), np.abs(k.fd2G(xs))
        H = np.maximum(d2, d1 / xs)
        for l in (0.7, 2.0, 3.3, 6.0, 12.0):
            mk = xs >= l
            assert k.tail(0, l) >= G[mk].max() * (1 - 1e-9) and k.tail(1, l) >= d1[mk].max() * (1 - 1e-9)
            assert k.tail(2, l) >= d2[mk].max() * (1 - 1e-9) and k.tail('H', l, 2) >= H[mk].max() * (1 - 1e-9)
            assert k.tail(0, l) <= G[mk].max() * 1.03 + 1e-12 and k.tail(1, l) <= d1[mk].max() * 1.03 + 1e-12
        assert k.tail(0, 5.0) >= k.tail(0, 6.0) and k.tail(1, 1.0) >= k.tail(1, 3.0) >= k.tail(1, 8.0)       # monotone envelopes
        assert 0.0 < k.tail(0, 40.0) <= k.tail(0, 20.0) and k.tail(2, 40.0) > 0.0                       # the value at L_MAX continues the table beyond it
        assert k.tail(1, 0.7) >= abs(k.fdG(1.4142)) * (1 - 1e-6) if k.name == 'gauss' else True            # the envelope sees the peak of psi to its right


def test_directed_own_constants():
    for k in (TK.family('gauss', 1), TK.family('mat32', 1), TK.MaternBessel(3, 2)):
        xs = np.linspace(1e-6, 1.2, 30001)
        G, d1, d2 = k.fG(xs), k.fdG(xs), k.fd2G(xs)
        for x in (0.1, 0.45, 0.9):
            m = xs <= x
            assert k.g_min(x) <= G[m].min() + 1e-12 and k.g_min(x) >= G[m].min() - 3e-3
            assert k.lampsi(x) <= (-d1[m] / xs[m]).min() + 1e-9 and k.lampsi(x) >= (-d1[m] / xs[m]).min() - 4e-3
            assert k.psihat(x) >= (-d1[m]).max() - 1e-9 and k.psihat(x) <= (-d1[m]).max() + 1.2e-2
            assert k.rho_min(x) <= (-d2[m]).min() + 1e-9 and k.rho_min(x) >= (-d2[m]).min() - 6e-3
            assert abs(k.psi_at(x) + float(k.fdG(x))) < 1e-9
        assert k.M1() >= (-d1).max() - 1e-9 and k.M2(1) >= np.abs(d2).max() - 1e-9
    g = TK.Gaussian(1)
    assert abs(g.M1() - 0.5 * math.sqrt(2) * math.exp(-0.5)) < 1e-2 and abs(g.M2(1) - 0.5) < 1e-3       # sup psi = 0.4289, sup |Gamma''| = 1/2
    assert abs(g.lam0 - 0.5) < 1e-15 and abs(TK.family('exp', 1).s - 1 / math.sqrt(2)) < 1e-15


def test_cusp_conditions_against_a_hand_calculation():
    e = TK.Exp1(1)
    K, R, D, r = 2, 2.0, 9.0, 0.5
    b = TK.box_cusp(e, K, R, D, r)
    ell = D - 2 * r
    S0 = math.exp(-ell)
    vbar = R * (1 + S0) / (1 - S0)
    Ev = S0 * (R + vbar)
    Ec = S0 * (R + vbar)
    assert abs(b['Ev'] - Ev) < 1e-9 and abs(b['Ec'] - Ec) < 1e-9
    assert abs(b['vlo'] - (math.exp(-r) - Ev)) < 1e-9 and abs(b['Cmargin'] - (math.exp(-r) - Ec)) < 1e-9 and b['ok']
    t = [TK.entry_threshold_cusp(e, K_, R_, D0 / 2.4) for (K_, R_) in [(2, 1.0), (2, 2.4), (3, 1.0)]]
    assert t[0] < t[1] and t[0] < t[2]
    tt = TK.entry_threshold_cusp(e, 2, 1.0, D0 / 2.4)
    rr = TK.r_grid(D0 / 2.4, 3.0, 40)
    assert any(TK.box_cusp(e, 2, 1.0, tt, r_)['ok'] for r_ in rr) and not any(TK.box_cusp(e, 2, 1.0, tt - 0.05, r_)['ok'] for r_ in rr)


def test_clean_corollary_constants_and_decay_modulus():
    g = TK.Gaussian(1)
    assert abs(TK.clean_cstar(g, 1.0, 0.45) - 0.0336) < 1.5e-3 and TK.clean_cstar(g, 10.0, 0.45) < TK.clean_cstar(g, 2.4, 0.45) < TK.clean_cstar(g, 1.0, 0.45)
    assert not TK.clean_ok(g, 1.0, 0.45, 0.3) and TK.clean_ok(g, 1.0, 0.45, 0.01) and not TK.clean_ok(g, 1.0, 0.02, 0.01)
    st = Stub(1, lamps_r=0.1)                                                  # lambda_psi(r) = 0.1: d_1 = 2.7 c_1/0.1 must stay below r
    assert TK.clean_ok(st, 1.0, 0.5, 0.018) and not TK.clean_ok(st, 1.0, 0.5, 0.0186)
    assert abs(TK.clean_index(g, 3, 2.0, 8.0, 0.5) - 2 * 4.0 * g.tail(2, 7.0)) < 1e-12
    # extremal case: at the smallest Delta with index <= c_*, the sharp box holds
    for K, R in [(2, 1.0), (3, 2.4)]:
        r = 0.5
        cs = TK.clean_cstar(g, R, r)
        D = next(D_ for D_ in np.arange(2 * r + 2.0, 20.0, 0.02) if TK.clean_index(g, K, R, D_, r) <= cs)
        assert TK.box(g, K, R, float(D), r)['ok']
        assert TK.clean_threshold(g, K, R, D0) >= TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=40.0)
    C = TK.decay_modulus(g, 2, 0.24, l0=3.0)
    assert 9.0 < C < 10.0
    for l in (3.0, 5.0, 9.0):
        m = max(g.tail(0, l), g.tail(1, l), g.tail(2, l))
        assert m <= C * math.exp(-0.24 * l * l) * (1 + 1e-9)
    assert TK.decay_modulus(g, 2, 0.25, l0=3.0) > C                                            # a larger rate needs a larger constant
