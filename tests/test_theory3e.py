"""Theory round 3, direction E: the entry layer and the limit path for a general Gram (experiments/theory3e_*.py).
The Gram classes (closed forms, identification with the autocorrelation of the basis function), the directed bounds, the box conditions (Gaussian regression against the paper's table, monotonicity,
clean corollary implies sharp box), the inner flow inside the box (Proposition G numerically), the cusp, the generic limit path and network against the repository's Gaussian code, and the
deterministic report."""
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import theory3e_kernels as TK      # noqa: E402
import theory3e_wendland as TW     # noqa: E402
import theory3e_path as TP         # noqa: E402
import theory3e_run as RUN         # noqa: E402
import theory3e_t1gen as TG        # noqa: E402
import exp_theory_t1 as T1         # noqa: E402
from flint import arb              # noqa: E402

D0 = 0.4167


def all_kernels():
    return [TK.family('gauss', 1), TK.family('exp', 1), TK.family('exp', 2), TK.family('exp', 3), TK.family('mat32', 1), TK.family('mat52', 2),
            TK.RationalQuadratic(3, 1), TW.Wendland(1), TW.Wendland(3)]


# ----------------------------------------------------------------------------------------------------------------------------------
# the Gram classes
# ----------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('k', all_kernels(), ids=lambda k: '%s-n%d' % (k.name, k.n))
def test_gram_values_derivatives_and_normalisation(k):
    assert abs(float(k.fG(0.0)) - 1.0) < 1e-12
    assert abs(float(k.fd2G(0.0)) + 0.5) < 1e-9                      # sigma-equivalent units: -Gamma''(0) = 1/2
    assert abs(float(k.fdG(0.0))) < 1e-12
    x = np.array([0.15, 0.4, 0.9, 1.6, 2.7])
    h = 1e-5
    assert np.max(np.abs((k.fG(x + h) - k.fG(x - h)) / (2 * h) - k.fdG(x))) < 1e-7
    assert np.max(np.abs((k.fdG(x + h) - k.fdG(x - h)) / (2 * h) - k.fd2G(x))) < 1e-6
    for xv in x:                                                     # the ball (arb) evaluation agrees with the float one
        assert abs(float(k.G(arb(float(xv))).mid()) - float(k.fG(xv))) < 1e-9
        assert abs(float(k.G1(arb(float(xv))).mid()) - float(k.fdG(xv))) < 1e-9
        assert abs(float(k.G2(arb(float(xv))).mid()) - float(k.fd2G(xv))) < 1e-8
    assert np.all(k.fdG(np.linspace(0.01, 20, 400)) <= 1e-12)        # radially decreasing


@pytest.mark.parametrize('k', all_kernels(), ids=lambda k: '%s-n%d' % (k.name, k.n))
def test_gram_matrices_are_positive_semidefinite_in_R_n(k):
    rng = np.random.default_rng(3)
    pts = rng.uniform(-6, 6, (40, k.n))
    r = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    ev = np.linalg.eigvalsh(k.fG(r))
    assert ev.min() > -1e-9


def _acorr(phi, n, d, amax=45.0):
    from scipy import integrate
    if n == 1:
        v = integrate.quad(lambda x: phi(abs(x)) * phi(abs(x - d)), -amax, amax, points=[0, d], limit=400)[0]
        z = integrate.quad(lambda x: phi(abs(x)) ** 2, -amax, amax, points=[0], limit=400)[0]
        return v / z
    f = lambda b, a: a * b * phi(a) * phi(b)
    v = integrate.dblquad(f, 0, amax, lambda a: abs(a - d), lambda a: a + d, epsabs=1e-12, epsrel=1e-9)[0] * 2 * np.pi / d
    z = 4 * np.pi * integrate.quad(lambda a: a * a * phi(a) ** 2, 0, amax)[0]
    return v / z


def test_gram_of_the_basis_function_is_the_stated_matern():
    """The Gram of exp(-|x|) in R^1, R^3 and of the Matern-3/2 basis in R^1, R^3 is the autocorrelation; family() gives it (raw units, spectral parameter 1)."""
    for n in (1, 3):
        g = TK.family('exp', n)
        for d in (0.5, 1.7):
            assert abs(_acorr(lambda a: np.exp(-a), n, d) - g.raw(0, d)) < 1e-6
    for n in (1, 3):
        g = TK.family('mat32', n)
        for d in (0.8, 2.5):
            assert abs(_acorr(lambda a: (1 + a) * np.exp(-a), n, d, 60.0) - g.raw(0, d)) < 1e-6


def test_matern_half_and_bessel_branches_agree_on_the_same_nu():
    a, b = TK.MaternHalf(3, 1), TK.MaternBessel(3.5, 1)                # nu = 7/2 through both branches
    for x in (0.1, 0.5, 1.0):
        assert abs(a.g_min(x) - b.g_min(x)) < 1e-5
        assert abs(a.lampsi(x) - b.lampsi(x)) < 4e-3 * a.lampsi(x)
        assert abs(a.psihat(x) - b.psihat(x)) < 1.5e-2 * a.psihat(x)
    for l in (2.0, 3.0, 5.0):
        for o in (0, 1, 2, 'H'):
            assert abs(a.tail(o, l) - b.tail(o, l)) < 2e-2 * a.tail(o, l)
    assert TK.entry_threshold(a, 2, 1.0, D0, lo_=1, hi=30) == TK.entry_threshold(b, 2, 1.0, D0, lo_=1, hi=30)


def test_wendland_gram_is_exact_piecewise_polynomial_c2_with_compact_support():
    for n in (1, 3):
        w = TW.Wendland(n)
        for b in w.bps[1:-1]:                                          # C^2 at the breakpoints (exact in rationals)
            bf = float(b)
            for o in (0, 1, 2):
                assert abs(float(w.rawf(o, np.array(bf - 1e-7))) - float(w.rawf(o, np.array(bf + 1e-7)))) < 5e-5
        assert float(w.rawf(0, np.array(2.0))) == 0.0 and float(w.rawf(0, np.array(2.5))) == 0.0
        assert abs(w.lam0 - (4.8 if n == 1 else 6.753246753246753)) < 1e-12
        # quadrature check of the autocorrelation
        phi = (lambda a: max(1 - a, 0) ** 3 * (3 * a + 1)) if n == 1 else (lambda a: max(1 - a, 0) ** 4 * (4 * a + 1))
        for d in (0.3, 1.1):
            assert abs(_acorr(phi, n, d, 1.0 + d) - float(w.rawf(0, np.array(d)))) < 1e-7
    assert TW.lam0_n2() > 0


# ----------------------------------------------------------------------------------------------------------------------------------
# directed bounds
# ----------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('k', [TK.family('gauss', 1), TK.family('exp', 1), TK.family('exp', 2), TK.family('mat32', 3), TW.Wendland(1)], ids=lambda k: '%s-n%d' % (k.name, k.n))
def test_tail_envelopes_dominate_the_dense_supremum(k):
    xs = np.linspace(0.5, 39.9, 60000)
    G, d1, d2 = np.abs(k.fG(xs)), np.abs(k.fdG(xs)), np.abs(k.fd2G(xs))
    H = np.maximum(d2, d1 / xs)
    for l in (0.7, 2.0, 3.3, 6.0, 12.0):
        m = xs >= l
        assert k.tail(0, l) >= G[m].max() * (1 - 1e-9)
        assert k.tail(1, l) >= d1[m].max() * (1 - 1e-9)
        assert k.tail(2, l) >= d2[m].max() * (1 - 1e-9)
        assert k.tail('H', l, 2) >= H[m].max() * (1 - 1e-9)
        assert k.tail(0, l) <= max(G[m].max() * (1.25 if k.name.startswith('wendland') else 1.08), 1e-12) + 1e-12      # and not wildly loose (the compact Gram's degree-9 pieces are enclosed less tightly)


@pytest.mark.parametrize('k', [TK.family('gauss', 1), TK.family('exp', 1), TK.family('exp', 2), TK.family('mat32', 2), TW.Wendland(3)], ids=lambda k: '%s-n%d' % (k.name, k.n))
def test_own_constants_bracket_the_dense_values(k):
    xs = np.linspace(1e-6, 1.2, 40001)
    G, d1, d2 = k.fG(xs), k.fdG(xs), k.fd2G(xs)
    for x in (0.1, 0.35, 0.7, 1.0):
        m = xs <= x
        assert k.g_min(x) <= G[m].min() + 1e-12 and k.g_min(x) >= G[m].min() - 3e-3
        assert k.lampsi(x) <= (-d1[m] / xs[m]).min() + 1e-9 and k.lampsi(x) >= (-d1[m] / xs[m]).min() - 3e-3
        assert k.psihat(x) >= (-d1[m]).max() - 1e-9 and k.psihat(x) <= (-d1[m]).max() + 1.2e-2
        assert k.rho_min(x) <= (-d2[m]).min() + 1e-9 and k.rho_min(x) >= (-d2[m]).min() - 6e-3
        assert abs(k.psi_at(x) + float(k.fdG(x))) < 1e-9


# ----------------------------------------------------------------------------------------------------------------------------------
# the box conditions
# ----------------------------------------------------------------------------------------------------------------------------------
def test_gaussian_row_sum_form_reproduces_the_papers_table():
    g = TK.Gaussian(1)
    paper = {(2, 1.0): 4.66, (2, 2.4): 5.50, (2, 10.0): 6.70, (3, 1.0): 5.08, (3, 2.4): 5.84}
    for (K, R), v in paper.items():
        assert abs(TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=24.0, hc='rowsum') - v) < 1e-9
    # the block-norm bound is never looser, and the Gaussian has no n-dependence
    for (K, R), v in list(paper.items())[:3]:
        t1 = TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=24.0)
        assert t1 <= v + 1e-9 and v - t1 < 0.15
        assert TK.entry_threshold(TK.Gaussian(3), K, R, D0, lo_=3.0, hi=24.0) == t1


def test_thresholds_are_monotone_and_ordered():
    g, e1, m32, m52 = TK.family('gauss', 1), TK.family('exp', 1), TK.family('mat32', 1), TK.family('mat52', 1)
    th = lambda k, K=2, R=1.0, d0=D0, **kw: TK.entry_threshold(k, K, R, d0, lo_=1.0, hi=60.0, **kw)
    assert th(g) < th(m52) < th(m32) < th(e1)                        # smoother Gram -> smaller threshold; the Gaussian is the limit
    assert th(g, K=3) > th(g) and th(g, R=2.4) > th(g) and th(g, d0=0.8) > th(g) and th(g, d0=0.0) < th(g)
    assert th(g, cond_c=4.0) > th(g, cond_c=2.0) > th(g)
    # once the box holds it holds for every larger separation (monotone in Delta at fixed r)
    for k in (e1, TK.family('exp', 2)):
        ok = [TK.box(k, 2, 1.0, D, 0.45)['ok'] for D in np.arange(4.0, 12.0, 0.25)]
        f = ok.index(True)
        assert all(ok[f:]) and not any(ok[:f])


def test_dormant_budget_shrinks_the_box_and_arithmetic_is_certain():
    g = TK.Gaussian(1)
    a, b = TK.box(g, 2, 2.4, 6.5, 0.45), TK.box(g, 2, 2.4, 6.5, 0.45, Omega=0.3)
    assert a['ok'] and b['vlo'] < a['vlo'] and b['Ec'] > a['Ec']
    assert TK.box(g, 2, 1.0, 2.0, 0.45)['ok'] is False
    # the threshold with Omega > 0 is larger
    assert TK.entry_threshold(g, 2, 1.0, D0, lo_=3.0, hi=40, Omega=0.1) > TK.entry_threshold(g, 2, 1.0, D0, lo_=3.0, hi=40)


def test_compact_gram_threshold_is_below_support_plus_twice_the_offset():
    """With a compact Gram all tails vanish beyond the support, so the conditions reduce to Gamma(r) > 0, psi(r) > 0, lambda_0 > 0 there (exact decoupling: Delta = supp + 2 d_0); the ball-arithmetic
    threshold, which also uses the tiny tails near the support edge, is far below that (this needs the polynomial pieces expanded about the support edge: in monomials the enclosure is vacuous)."""
    for n in (1, 3):
        w = TW.Wendland(n)
        t = TK.entry_threshold(w, 2, 1.0, D0, lo_=1.0, hi=40.0)
        assert t < w.supp_eq + 2 * D0 - 1.5
        b = TK.box(w, 2, 1.0, w.supp_eq + 2 * D0 + 0.01, D0 * 1.000001, sums=(0.0, 0.0, 0.0))
        assert b['ok'] and b['Ev'] < 1e-300 and b['Hc'] < 1e-300                      # exactly decoupled: E_v = E_c = H_c = 0
        assert TK.entry_threshold(w, 3, 2.4, D0, lo_=1.0, hi=40.0) > t
        assert w.tail(0, w.supp_eq) < 1e-8 and w.tail(1, w.supp_eq + 0.01) < 1e-300
        assert w.tail(0, 3.3) < 1.5 * float(np.abs(w.fG(np.linspace(3.3, w.supp_eq, 4000))).max())


def test_clean_corollary_implies_the_sharp_box():
    rng = np.random.default_rng(5)
    for k in (TK.family('gauss', 1), TK.family('exp', 1), TK.family('mat32', 2)):
        for _ in range(12):
            K, R, r = int(rng.integers(2, 5)), float(rng.choice([1.0, 2.4, 6.0])), float(rng.uniform(0.45, 0.9))
            cs = TK.clean_cstar(k, R, r)
            if cs <= 0:
                continue
            for D in np.arange(2 * r + 3.0, 2 * r + 14.0, 1.0):
                if TK.clean_index(k, K, R, D, r) <= cs:
                    assert TK.box(k, K, R, float(D), r)['ok'], (k.name, K, R, r, D)
    g = TK.Gaussian(1)
    # clean thresholds are conservative
    for (K, R) in [(2, 1.0), (2, 2.4), (3, 1.0)]:
        assert TK.clean_threshold(g, K, R, D0) > TK.entry_threshold(g, K, R, D0, lo_=3.0, hi=40)


def test_cross_hessian_bound_is_below_twelve_times_the_coherence_index():
    g = TK.Gaussian(1)
    rng = np.random.default_rng(1)
    for _ in range(30):
        R = float(rng.choice([1.0, 2.0, 5.0]))
        c = float(rng.uniform(1e-4, 0.02))
        S = [c / R ** 2 * float(rng.uniform(0.2, 1.0)) for _ in range(3)]
        S[int(rng.integers(0, 3))] = c / R ** 2
        b = TK.box(g, 2, R, 8.0, 0.5, sums=S)
        if 'Hc' in b:
            assert b['Hc'] <= 12.0 * R * R * max(S) * (1 + 1e-12)


def test_cusp_conditions_have_no_hessian_part_and_are_monotone():
    e = TK.Exp1(1)
    t = [TK.entry_threshold_cusp(e, K, R, D0 / 2.4) for (K, R) in RUN.COMBOS]
    assert all(np.isfinite(t)) and t[0] < t[1] < t[2] and t[0] < t[3]
    assert TK.box_cusp(e, 2, 1.0, 10.0, 0.3)['ok'] and not TK.box_cusp(e, 2, 1.0, 1.0, 0.3)['ok']
    ok = [TK.box_cusp(e, 2, 1.0, D, 0.3)['ok'] for D in np.arange(1.0, 12.0, 0.25)]
    f = ok.index(True)
    assert all(ok[f:]) and not any(ok[:f])


def _F_fun(k, cs, As, K, n):
    """The fit of K units (one per feature, all features fitted) as a function of z = (v_1..v_K, c_1..c_K): the target is sum_k A_k phi(c*_k) and the units are v_i phi(c_i)."""
    def F(z):
        v, c = z[:K], z[K:].reshape(K, n)
        G0 = TP.gam_mat(k, cs, cs)
        F1 = 0.5 * As @ G0 @ As
        F2 = -(v @ (TP.gam_mat(k, c, cs) @ As))
        F3 = 0.5 * v @ TP.gam_mat(k, c, c) @ v
        return F1 + F2 + F3
    return F


@pytest.mark.parametrize('kind,n,K,R', [('gauss', 1, 2, 1.0), ('exp', 1, 3, 2.4), ('mat32', 2, 2, 1.0), ('exp', 2, 2, 2.4)])
def test_field_and_cross_hessian_bounds_of_the_proof_hold_at_random_points_of_the_box(kind, n, K, R):
    """Steps 1 and 4 of the proof checked on F itself (finite differences of the actual fit, independent of the formulas): |P_i| <= E_v, |C_i| <= E_c, ||H_cross|| <= H_c."""
    k = TK.family(kind, n)
    thr, r = TK.entry_threshold(k, K, R, D0, lo_=1.0, hi=40.0, return_r=True)
    b = TK.box(k, K, R, thr, r)
    assert b['ok']
    rng = np.random.default_rng(11)
    e1 = np.eye(n)[0]
    cs = (np.arange(K) * thr)[:, None] * e1[None, :]
    a = np.ones(K)
    a[1::2] = R
    As = a.copy()
    F = _F_fun(k, cs, As, K, n)
    dim = K * (1 + n)

    def grad(z, h=1e-6):
        g = np.zeros(dim)
        for q in range(dim):
            e = np.zeros(dim)
            e[q] = h
            g[q] = (F(z + e) - F(z - e)) / (2 * h)
        return g

    def hess(z, h=1e-4):
        H = np.zeros((dim, dim))
        for q in range(dim):
            e = np.zeros(dim)
            e[q] = h
            H[:, q] = (grad(z + e) - grad(z - e)) / (2 * h)
        return 0.5 * (H + H.T)
    worst = [0.0, 0.0, 0.0]
    for _ in range(15):
        d = rng.normal(size=(K, n))
        d = d / np.linalg.norm(d, axis=1, keepdims=True) * (r * rng.uniform(0, 1, (K, 1)) ** (1.0 / n))
        v = rng.uniform(b['vlo'], b['vbar'], K)
        z = np.r_[v, (cs + d).ravel()]
        g = grad(z)
        for i in range(K):
            dn = np.linalg.norm(d[i])
            own_v = a[i] * float(k.fG(dn)) - v[i]
            P = -g[i] - own_v
            gc = g[K + i * n: K + (i + 1) * n]
            own_c = v[i] * a[i] * (-float(k.fdG(dn))) * (d[i] / dn if dn > 0 else 0 * d[i])
            C = np.linalg.norm(gc - own_c) / v[i]
            worst[0] = max(worst[0], abs(P) / b['Ev'])
            worst[1] = max(worst[1], C / b['Ec'])
        H = hess(z)
        Hown = np.zeros_like(H)
        for i in range(K):
            dn = np.linalg.norm(d[i])
            dh = d[i] / dn
            psi_, rho_ = -float(k.fdG(dn)), -float(k.fd2G(dn))
            blk = np.zeros((1 + n, 1 + n))
            blk[0, 0] = 1.0
            blk[0, 1:] = blk[1:, 0] = a[i] * psi_ * dh
            blk[1:, 1:] = v[i] * a[i] * (rho_ * np.outer(dh, dh) + (psi_ / dn) * (np.eye(n) - np.outer(dh, dh)))
            idx = [i] + list(range(K + i * n, K + (i + 1) * n))
            Hown[np.ix_(idx, idx)] = blk
        worst[2] = max(worst[2], np.linalg.norm(H - Hown, 2) / b['Hc'])
    assert worst[0] <= 1.0 and worst[1] <= 1.0 and worst[2] <= 1.0, worst
    assert worst[2] > 0.0


# ----------------------------------------------------------------------------------------------------------------------------------
# the inner flow and the limit path
# ----------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('k,K,R', [(TK.family('gauss', 1), 2, 1.0), (TK.family('exp', 1), 2, 2.4), (TK.family('mat32', 2), 2, 1.0), (TW.Wendland(1), 2, 1.0)],
                         ids=['gauss', 'laplace', 'mat32-n2', 'wendland'])
def test_proposition_G_on_random_and_adversarial_starts(k, K, R):
    thr, r = TK.entry_threshold(k, K, R, D0, lo_=1.0, hi=40.0, return_r=True)
    b = TK.box(k, K, R, thr, r)
    assert b['ok']
    for sd in range(6):
        o = TP.box_run(k, K, R, thr, r, b['vlo'], b['vbar'], sd, mob=['1', 'v', 'v43', 'q'][sd % 4], muc=[0.3, 1.0, 3.0][sd % 3], n=k.n, tmax=1200.0, adv=(sd >= 3))
        assert o['in_d'] and o['in_hi'] and o['in_lo'] and o['mono'] and o['converged'], (sd, o)


def test_inner_flow_can_leave_a_too_small_box():
    """Not vacuous: far below the threshold, a start on the boundary of the box does leave it (the conditions are not a theorem about every Delta)."""
    g = TK.Gaussian(1)
    o = TP.box_run(g, 2, 1.0, 2.0, 0.42, 0.8, 1.2, 1, mob='v', muc=1.0, n=1, tmax=800.0, adv=True)
    assert not (o['in_d'] and o['in_lo'] and o['in_hi'])


def test_generic_limit_path_equals_the_repositorys_gaussian_path():
    g = TK.Gaussian(1)
    cs, As, c0 = T1.make_target(2, 8.0, 1.0, 13)
    L1 = T1.limit_path(cs, As, c0)
    L2 = TP.limit_path(g, cs[:, None], As, c0[:, None])
    assert [e['kernel'] for e in L1['events']] == [e['kernel'] for e in L2['events']]
    assert np.allclose([e['tau'] for e in L1['events']], [e['tau'] for e in L2['events']], rtol=1e-9)
    assert np.allclose(L1['events'][-1]['w'], L2['events'][-1]['w'], atol=1e-9)
    assert all(e['block_ok'] for e in L2['events'])
    assert L2['resid'] < 1e-12
    assert abs(L2['final_w'][np.argsort([e['feature'] for e in L2['events']])] - As).max() < 1e-9


def test_planar_path_is_rotation_invariant_and_ends_at_the_target():
    g = TK.Gaussian(2)
    cs, As, c0 = TP.make_target(2, 2, 7.0, 1.0, 5, margin=4.0)
    L = TP.limit_path(g, cs, As, c0)
    th = 0.7
    Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    L2 = TP.limit_path(g, cs @ Rm.T, As, c0 @ Rm.T)
    assert np.allclose([e['tau'] for e in L['events']], [e['tau'] for e in L2['events']], rtol=1e-9)
    feats = [e['feature'] for e in L['events']]
    assert sorted(feats) == [0, 1] and np.abs(L['final_w'] - As[feats]).max() < 1e-9
    assert np.abs(L['final_c'] - cs[feats]).max() < 1e-9


@pytest.mark.parametrize('k', [TK.family('gauss', 1), TK.family('exp', 1), TW.Wendland(1)], ids=['gauss', 'laplace', 'wendland'])
def test_network_follows_the_limit_path_at_small_eps(k):
    thr = TK.entry_threshold(k, 2, 1.0, 0.5, lo_=1.0, hi=60.0)
    cs, As, c0 = TP.make_target(1, 2, float(np.ceil((thr + 0.6) * 10) / 10), 1.0, 5)
    L = TP.limit_path(k, cs, As, c0)
    eps = 1e-7
    sol, U = TP.simulate_network(k, cs, As, c0, eps, L, nglob=1500, nloc=800)
    c = TP.compare(cs, As, c0, eps, L, sol, U)
    assert c['formed_ok'] and c['exact'] and c['order_ok']
    assert max(abs(x) for x in c['shifts']) < 0.05


def test_smoothed_cusp_inner_problem_locks_the_centres():
    r = RUN.part_D(quick=True)
    d = r['D_smooth']
    assert d[0, 2] < 0.05 and d[0, 3] < 1e-6                           # final offset O(delta), coefficient error ~ 0
    assert np.all(np.isfinite(r['D_thr']))


def test_t1_certificate_for_the_gaussian_gram_agrees_with_the_repositorys():
    g = TK.Gaussian(1)
    M = TG.build(g)
    cs, As, c0 = T1.make_target(2, 8.0, 1.0, 13)
    p0, p1 = T1.prepare(cs, As, c0), M.prepare(cs, As, c0)
    for ex in (-8.0, -10.0):
        o0, i0 = T1.best_cert(p0, 10.0 ** ex)
        o1, i1 = M.best_cert(p1, 10.0 ** ex)
        assert o0['ok'] and o1['ok'] and i0 == i1
        assert abs(o1['D'] / o0['D'] - 1) < 0.06                                   # the Gram version uses sup psi = 0.4319 (ball enclosure) for 0.43 and the exact block bounds
    o0, _ = T1.best_cert(p0, 1e-5)
    o1, _ = M.best_cert(p1, 1e-5)
    assert (not o0['ok']) and (not o1['ok'])                                      # both fail above the threshold
    # the box of Lemma 4.1 for the Gaussian Gram: same verdict as the repository's on a grid
    for (K, R, D, r) in [(2, 1.0, 5.0, 0.5), (2, 1.0, 4.6, 0.5), (3, 2.4, 7.5, 0.5), (2, 10.0, 7.5, 0.6)]:
        a, b = T1.box_t1(K, R, D, r, 0.02), TG.box_t1_gen(g, K, R, D, r, 0.02)
        assert bool(a.get('ok')) == bool(b.get('ok')), (K, R, D, r)
        if a.get('ok'):
            assert abs(a['beta_v'] - b['beta_v']) < 0.02 and abs(a['beta_d'] - b['beta_d']) < 0.02


def test_t1_certificate_runs_for_a_non_gaussian_gram():
    k = TK.family('mat32', 1)
    M = TG.build(k)
    d_t1 = M.delta_t1(2, 1.0, 0.5, 0.02, step=0.2, lo=3.0, hi=30.0)
    assert d_t1 == d_t1 and d_t1 > T1.delta_t1(2, 1.0, 0.5, 0.02, step=0.2, lo=3.0, hi=30.0) - 0.5
    cs, As, c0 = M.make_target(2, float(np.ceil(d_t1 + 1.0)), 1.0, 9003)
    prep = M.prepare(cs, As, c0)
    M.PARS = M.PARS[:4]
    o, pi = M.best_cert(prep, 1e-10)
    assert o['ok'] and pi >= 0
    sol, U = M.simulate(prep, 1e-10, nlocal=200)
    assert M.follows(prep, 1e-10, sol, U)[0]


# ----------------------------------------------------------------------------------------------------------------------------------
# the report
# ----------------------------------------------------------------------------------------------------------------------------------
def test_report_regenerates_byte_for_byte_from_the_stored_raw(tmp_path):
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory3e_raw.npz')
    md = os.path.join(REPO, 'results', 'synthetic', 'theory3e.md')
    if not (os.path.exists(raw) and os.path.exists(md)):
        pytest.skip('stored report not present')
    shutil.copy(raw, tmp_path / 'theory3e_raw.npz')
    env = dict(os.environ, T3E_OUT=str(tmp_path))
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'theory3e_run.py'), '--from-raw'], check=True, env=env, cwd=REPO)
    assert open(tmp_path / 'theory3e.md').read() == open(md).read()


def test_report_numbers_are_those_of_the_stored_raw(tmp_path):
    raw = os.path.join(REPO, 'results', 'synthetic', 'theory3e_raw.npz')
    if not os.path.exists(raw):
        pytest.skip('stored raw not present')
    with np.load(raw) as z:
        A = z['A_thr']
        assert np.allclose(z['A_gauss_rowsum'], z['A_paper'], atol=1e-9)
        assert A[0, 0] < A[9, 0] < A[6, 0] < A[3, 0]                       # Gaussian < Matern-5/2 basis < Matern-3/2 basis < Laplace basis (n = 1)
    g = TK.family('gauss', 1)
    with np.load(raw) as z:
        assert abs(z['A_thr'][0, 0] - TK.entry_threshold(g, 2, 1.0, D0, lo_=1.0, hi=60.0)) < 1e-9
