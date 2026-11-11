"""Proposition D‴.58: a feature wider than the kernels as a Gaussian quadrature of its mixing law (`src/core/wide_feature.py`)."""
from math import factorial

import numpy as np

from src.core import wide_feature as W

SIG = 0.12


def ss(r_):
    return SIG * np.sqrt(1 + r_ ** 2)


def test_half_squeeze_error_is_the_mehler_series():
    """(i): the half-squeeze rule's L² error (trapezoid) is the exact Mehler series, at τ/σ = 1/2, 1, 2 and K ≤ 6."""
    for r_ in (0.5, 1.0, 2.0):
        t = W.params(SIG, ss(r_))['t']
        for K in range(1, 7):
            e = W.l2_error(*W.half_squeeze(K, SIG, ss(r_)), SIG, ss(r_))
            s2, f2 = W.mehler_series(K, t)
            assert abs(e / np.sqrt(s2) - 1) < 1e-9
            assert 0.97 < f2 / s2 <= 1.0 and abs(f2 - W.first_term(K, t)) < 1e-12 * f2


def test_gh_defect_vanishes_below_2k_and_its_tail_ratios_are_integer_polynomials():
    """D_{n,K} = 0 for odd n and n < 2K, D_{2K,K} = 2^K K!, D_{2K+2,K} = −(K − 1)D_{2K,K}, D_{2K+4,K} = −(K² − 25K − 6)/2·D_{2K,K}."""
    for K in range(1, 13):
        assert all(W.gh_defect(n, K) == 0 for n in range(0, 2 * K)) and W.gh_defect(2 * K + 1, K) == 0
        D = W.gh_defect(2 * K, K)
        assert D == 2 ** K * factorial(K)
        assert W.gh_defect(2 * K + 2, K) == -(K - 1) * D
        assert 2 * W.gh_defect(2 * K + 4, K) == -(K * K - 25 * K - 6) * D


def test_half_squeeze_rule_is_the_gauss_hermite_rule_of_the_reweighted_law():
    """Nodes τ(σ/σ*)^{1/2}z_k, symmetric; amplitudes positive, A(σ*/σ)^{1/2} for K = 1; moment matching's nodes τz_k."""
    y, a = W.half_squeeze(1, SIG, ss(1.0), A=0.8)
    assert abs(y[0]) < 1e-15 and abs(a[0] - 0.8 * np.sqrt(ss(1.0) / SIG)) < 1e-14
    z, w = W.gh(5)
    y, a = W.half_squeeze(5, SIG, ss(2.0))
    tau = W.params(SIG, ss(2.0))['tau']
    assert np.allclose(y, tau * np.sqrt(SIG / ss(2.0)) * z, rtol=0, atol=1e-15) and np.allclose(y, -y[::-1], atol=1e-15)
    assert np.all(a > 0) and np.allclose(a, a[::-1], rtol=1e-13)
    assert np.allclose(W.moment_matching(5, SIG, ss(2.0))[0], tau * z, rtol=0, atol=1e-15)


def test_rates_per_kernel():
    """At τ/σ = 1 the half-squeeze rule's error falls by 5.64 at K = 8 (→ 1/t = 5.83, the O(1/K) deficit 1/(4Kt)), moment matching's
    by 1/λ = 3 to 0.2%."""
    s_ = ss(1.0)
    P = W.params(SIG, s_)
    e = [W.l2_error(*W.half_squeeze(K, SIG, s_), SIG, s_) for K in (7, 8)]
    m = [W.l2_error(*W.moment_matching(K, SIG, s_), SIG, s_) for K in (7, 8)]
    assert abs(e[0] / e[1] - (1 - 1 / (4 * 8)) / P['t']) < 0.01 / P['t']
    assert abs(m[0] / m[1] * P['lam'] - 1) < 2e-3


def test_optimal_fit_improves_slightly_on_the_rule():
    """The optimal 4-kernel fit at τ/σ = 1: 0.44% below the half-squeeze rule's error, its nodes within 0.5% of the rule's."""
    s_ = ss(1.0)
    y, a = W.half_squeeze(4, SIG, s_)
    e = W.l2_error(y, a, SIG, s_)
    yo, ao, eo = W.optimal_fit(y, SIG, s_)
    assert 0.99 * e < eo < 0.997 * e
    assert np.abs(np.sort(yo) / np.sort(y) - 1).max() < 5e-3


def test_the_fast_defects_and_the_tail_ratio():
    """`gh_defects` (the remainders modulo He_K) equals `gh_defect`; the tail ratio is the series over its first term; the ratios
    D_{2K+2j,K}/D_{2K,K} are polynomials in K of degree j with the leading coefficient binom(1/2, j)(−2)^j = −C_{j−1}/2^{j−1} (j ≤ 5:
    interpolated on K = 1 … j + 1, checked on K = j + 2)."""
    from fractions import Fraction
    for K in range(1, 7):
        assert W.gh_defects(K, 2 * K + 10) == [W.gh_defect(n, K) for n in range(2 * K + 11)]
    s2, f2 = W.mehler_series(6, 0.5, n_terms=60)
    assert abs(W.tail_ratio(6, 0.5, n_terms=60) - s2 / f2) < 1e-12
    catalan = [1, 1, 2, 5, 14]
    for j in range(1, 6):
        pts = [(K, Fraction(W.gh_defects(K, 2 * K + 2 * j)[2 * K + 2 * j], W.gh_defects(K, 2 * K)[2 * K])) for K in range(1, j + 3)]
        # the j-th divided difference over K = 1 … j + 1 is the leading coefficient; the (j + 1)-th vanishes (degree j)
        def dd(ps):
            if len(ps) == 1:
                return ps[0][1]
            return (dd(ps[1:]) - dd(ps[:-1])) / (ps[-1][0] - ps[0][0])
        assert dd(pts[:j + 1]) == Fraction(-catalan[j - 1], 2 ** (j - 1)) and dd(pts) == 0


def test_the_optimal_fit_keeps_the_centre_node_at_odd_k_and_the_free_path_agrees():
    """`optimal_fit` at K = 3 keeps a node at the centre (symmetric: ±y and 0) and improves on the rule; without the symmetry it reaches
    the same error from the same start."""
    s_ = ss(1.0)
    y, a = W.half_squeeze(3, SIG, s_)
    yo, ao, eo = W.optimal_fit(y, SIG, s_)
    assert yo.size == 3 and abs(np.sort(yo)[1]) < 1e-15 and eo < W.l2_error(y, a, SIG, s_)
    y2, a2 = W.half_squeeze(2, SIG, s_)
    e_sym = W.optimal_fit(y2, SIG, s_)[2]
    e_free = W.optimal_fit(y2, SIG, s_, symmetric=False)[2]
    assert abs(e_free / e_sym - 1) < 1e-8


def test_least_squares_amplitudes_are_the_gram_solve():
    """`ls_amplitudes` is the L² projection: G a = b with G_ij = √π σ e^{−(y_i − y_j)²/(4σ²)}, b_i = A√(2π)σσ*/(σ² + σ*²)^{1/2}
    e^{−y_i²/(2(σ² + σ*²))} — on six nodes 0.25σ apart, where the design's condition number is ~10⁴ and a cut-off at 10⁻³ of its
    spectrum would change the amplitudes entirely."""
    s_, A = ss(1.0), 0.8
    y = SIG * 0.25 * (np.arange(6) - 2.5)
    x, h = W.grid(SIG, s_)
    sv = np.linalg.svd(np.exp(-(x[:, None] - y[None, :]) ** 2 / (2 * SIG ** 2)), compute_uv=False)
    assert 1e4 < sv[0] / sv[-1] < 1e7
    G = np.sqrt(np.pi) * SIG * np.exp(-(y[:, None] - y[None, :]) ** 2 / (4 * SIG ** 2))
    b = A * np.sqrt(2 * np.pi) * SIG * s_ / np.sqrt(SIG ** 2 + s_ ** 2) * np.exp(-y ** 2 / (2 * (SIG ** 2 + s_ ** 2)))
    a = W.ls_amplitudes(y, SIG, s_, A)
    assert np.abs(a - np.linalg.solve(G, b)).max() < 1e-6 * np.abs(a).max()


def test_the_fast_defects_beyond_k_21():
    """`gh_defects` against the direct `gh_defect` at K = 25, 30, 40 (Part D's K ≥ 40 rests on them)."""
    for K in (25, 30, 40):
        D = W.gh_defects(K, 2 * K + 12)
        for n in (2 * K, 2 * K + 2, 2 * K + 8, 2 * K + 12):
            assert D[n] == W.gh_defect(n, K)


def test_scaled_rules_and_the_best_scale():
    """`scaled_rule` at v_hs is the half-squeeze rule and at τ² moment matching; the Fock scaling's rate (2 − λ)/λ is (1 − t + t²)/t, slower
    than 1/t; the best variance's error lies between the optimal fit's and the rule's, and at K = 1 it is the optimum."""
    s_ = ss(1.0)
    V = W.scalings(SIG, s_)
    for K in (3, 6):
        assert np.allclose(np.c_[W.scaled_rule(K, SIG, s_, V['hs'])], np.c_[W.half_squeeze(K, SIG, s_)], rtol=1e-13, atol=1e-16)
        assert np.allclose(np.c_[W.scaled_rule(K, SIG, s_, V['mm'])], np.c_[W.moment_matching(K, SIG, s_)], rtol=1e-13, atol=1e-16)
    P = W.params(SIG, s_)
    t, lam = P['t'], P['lam']
    assert abs((2 - lam) / lam - (1 - t + t * t) / t) < 1e-12 and (2 - lam) / lam < 1 / t
    y, a = W.half_squeeze(4, SIG, s_)
    e_hs = W.l2_error(y, a, SIG, s_)
    vb, eb = W.best_scale(4, SIG, s_)
    eo = W.optimal_fit(y, SIG, s_)[2]
    assert eo <= eb * (1 + 1e-9) and eb < e_hs and 0.98 < vb / V['hs'] < 1.0
    e1 = W.optimal_fit(*W.half_squeeze(1, SIG, s_)[:1], SIG, s_)[2]
    assert abs(W.best_scale(1, SIG, s_)[1] / e1 - 1) < 1e-6


def test_the_frame_continues_to_narrower_features():
    """(v): for σ* < σ the Gauss–Hermite cosine sum at the frequencies (2σ²|sinh r|)^{1/2}z_k has the error of the Mehler series at
    |t| = (σ − σ*)/(σ + σ*) — the half-squeeze rule with imaginary nodes in conjugate pairs."""
    for q in (0.5, 0.8):
        s_ = q * SIG
        t = (1 - q) / (1 + q)
        for K in (1, 3, 6):
            e = W.narrow_error(*W.narrow_rule(K, SIG, s_), SIG, s_)
            assert abs(e / np.sqrt(W.mehler_series(K, t)[0]) - 1) < 1e-10
    eta, b = W.narrow_rule(4, SIG, 0.5 * SIG)
    assert np.allclose(eta, -eta[::-1]) and np.allclose(b, b[::-1])


def test_the_scalings_closed_forms_and_the_best_variance():
    """Referee round 54: the Fock scaling's variance v_hs/cosh r is 2τ²σ²/(σ² + σ*²); the Mercer-zero scaling's nodes are the zeros of
    the K-th Mercer eigenfunction of the rule's kernel exp(−(y − y′)²/(4σ²)) (⟨φ_σ(· − y), φ_σ(· − y′)⟩ ∝ it) under the mixing law
    N(0, τ²): Fasshauer–McCourt's H_K(αβy), α² = 1/(2τ²), β = (1 + (2ε/α)²)^{1/4}, ε² = 1/(4σ²). `best_scale`'s variance is the least
    error's: no point of a grid of spacing 10⁻⁴ in log v within 10⁻² of it is lower."""
    from numpy.polynomial.hermite import hermgauss
    for r_ in (0.5, 1.0, 2.0):
        s_ = ss(r_)
        tau2 = s_ ** 2 - SIG ** 2
        V = W.scalings(SIG, s_)
        assert abs(V['fock'] / (2 * tau2 * SIG ** 2 / (SIG ** 2 + s_ ** 2)) - 1) < 1e-14
        al, eps = 1 / np.sqrt(2 * tau2), 1 / (2 * SIG)
        beta = (1 + (2 * eps / al) ** 2) ** 0.25
        for K in (3, 6):
            assert np.allclose(W.scaled_rule(K, SIG, s_, V['mercer'])[0], hermgauss(K)[0] / (al * beta), rtol=1e-13, atol=1e-16)
        for K in (2, 4):
            vb, eb = W.best_scale(K, SIG, s_)
            es = [W.l2_error(*W.scaled_rule(K, SIG, s_, vb * np.exp(h)), SIG, s_) for h in np.linspace(-1e-2, 1e-2, 201)]
            assert eb <= min(es) * (1 + 1e-12)
