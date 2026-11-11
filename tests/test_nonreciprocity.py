"""
Checks for PARADIGM.md Theorem I (reciprocity criterion of the exactly averaged loop at order 1/k_e²), using
src/core/nonreciprocity.py against the generic second-variation code of src/core/averaging.py.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import langevin as L
from src.core import averaging as A
from src.core import nonreciprocity as NR

PER = 2 * np.pi


def _plant_1d():
    f = lambda X: np.sin(2 * X) + 0.3 * X
    J = lambda X: (2 * np.cos(2 * X) + 0.3)[:, :, None]
    xr = lambda t: np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])
    xd = lambda t: np.array([1.5 * np.cos(t) + 1.0 * np.cos(2 * t)])
    return f, J, xr, xd


def _plant_2d(orient=+1):
    def f(X):
        return np.stack([np.sin(X[:, 1]) + 0.3 * X[:, 0], np.cos(X[:, 0]) - 0.2 * X[:, 1]], 1)

    def J(X):
        out = np.zeros((len(X), 2, 2))
        out[:, 0, 0] = 0.3; out[:, 0, 1] = np.cos(X[:, 1]); out[:, 1, 0] = -np.sin(X[:, 0]); out[:, 1, 1] = -0.2
        return out
    xr = lambda t: np.array([1.2 * np.cos(t), orient * 0.8 * np.sin(t)])
    xd = lambda t: np.array([-1.2 * np.sin(t), orient * 0.8 * np.cos(t)])
    return f, J, xr, xd


def _segment_orbit():
    """A segment traversed back and forth: its 1-current vanishes (T_𝒪 = 0), so K ≡ 0 although n = 2."""
    xr = lambda t: np.array([np.cos(t), 0.5 * np.cos(t)])
    xd = lambda t: np.array([-np.sin(t), -0.5 * np.sin(t)])
    return xr, xd


def _setup(n, seed=0, k=5):
    rng = np.random.default_rng(seed)
    net = L.SaturatedLiftedNetwork(n, a_max=3.0, c_max=2.0, s_min=0.3, s_max=1.5)
    return net, rng.normal(0, 1.0, (k, 2 * n + 1)), rng


def test_bilinear_form_B_reproduces_the_antisymmetric_second_variation():
    """A(ϑ, ϑ') = aᵀB(ξ, ξ')a' equals the antisymmetric part of S_even + K computed by averaging.py, for n = 1 and 2."""
    for n in (1, 2):
        net, Th, rng = _setup(n)
        f, J, xr, xd = _plant_1d() if n == 1 else _plant_2d()
        orb = A.Orbit(xr, xd, PER, N=1500, sub=1)
        fm = A.FrozenMeasure(net, Th)
        TE = rng.normal(0, 1.0, (4, 2 * n + 1))
        Aref = NR.A_from_second_variation(f, J, fm, orb, TE)
        for i in range(4):
            for j in range(4):
                a = NR.A_kernel(f, J, fm, orb, TE[i], TE[j])
                assert abs(a - Aref[i, j]) < 1e-10 * (1 + np.abs(Aref).max()), (n, i, j, a, Aref[i, j])
        assert np.abs(Aref).max() > 1e-3                                    # non-trivial instance


def test_zero_feature_collapses_the_cyclic_criterion_to_a_pairwise_one():
    """A(ϑ, ϑ_z) = 0 for the zero feature ϑ_z = (a = 0, ξ), so the cyclic sum over (ϑ₀, ϑ₁, ϑ_z) equals A(ϑ₀, ϑ₁):
    a coboundary A must vanish identically."""
    for n in (1, 2):
        net, Th, rng = _setup(n)
        f, J, xr, xd = _plant_1d() if n == 1 else _plant_2d()
        orb = A.Orbit(xr, xd, PER, N=1500, sub=1)
        fm = A.FrozenMeasure(net, Th)
        th0, th1 = rng.normal(0, 1.0, 2 * n + 1), rng.normal(0, 1.0, 2 * n + 1)
        thz = NR.zero_feature(net, rng.normal(0, 1.0, n + 1))
        Aref = NR.A_from_second_variation(f, J, fm, orb, np.vstack([th0, th1, thz]))
        assert abs(Aref[0, 2]) < 1e-14 and abs(Aref[1, 2]) < 1e-14
        assert abs(NR.cyclic_sum(Aref) - Aref[0, 1]) < 1e-14 * (1 + abs(Aref[0, 1]))
        assert abs(Aref[0, 1]) > 1e-4


def test_n1_reciprocal_iff_residual_vanishes_on_the_orbit():
    """n = 1: with f := f̂_μ (perfect fit on the orbit) A ≡ 0; with a non-zero residual A ≠ 0, and A equals the
    closed form −a a' ((c − c')/s²) ∫ r φ_c φ_c' dρ for equal widths."""
    net, Th, rng = _setup(1)
    f, J, xr, xd = _plant_1d()
    orb = A.Orbit(xr, xd, PER, N=1500, sub=1)
    fm = A.FrozenMeasure(net, Th)
    # perfect fit: the plant IS the frozen network
    f_fit, J_fit = fm.f_hat, fm.jac
    worst = 0.0
    for _ in range(20):
        th, thp = rng.normal(0, 1.0, 3), rng.normal(0, 1.0, 3)
        worst = max(worst, abs(NR.A_kernel(f_fit, J_fit, fm, orb, th, thp)))
    assert worst < 1e-13, worst
    # non-zero residual: A ≠ 0 and matches the equal-width closed form
    found = 0.0
    for _ in range(20):
        th, thp = rng.normal(0, 1.0, 3), rng.normal(0, 1.0, 3); thp[2] = th[2]              # equal widths
        a, c, s = net.physical(th[None]); ap, cp, _ = net.physical(thp[None])
        val = NR.A_kernel(f, J, fm, orb, th, thp)
        cf = NR.A_even_1d_closed_form(f, fm, orb, a[0, 0], c[0, 0], s[0], ap[0, 0], cp[0, 0])
        assert abs(val - cf) < 1e-12 * (1 + abs(cf)), (val, cf)
        found = max(found, abs(val))
    assert found > 1e-3


def test_circulation_kernel_closed_form_and_vanishing_for_a_zero_current_orbit():
    """n = 2, equal widths: κ(ξ, ξ') = −(1/(2s²)) e^{−|c−c'|²/(4s²)} (c − c')·F_s((c + c')/2); κ ≢ 0 on the ellipse
    (F₂((0.5, 0)) > 0 by the positivity argument) and κ ≡ 0 on a back-and-forth segment (T_𝒪 = 0)."""
    _, _, xr, xd = _plant_2d()
    orb = A.Orbit(xr, xd, PER, N=4000, sub=1)
    rng = np.random.default_rng(1)
    for _ in range(20):
        c, cp = rng.uniform(-1.5, 1.5, 2), rng.uniform(-1.5, 1.5, 2); s = rng.uniform(0.4, 1.2)
        k1 = NR.kappa(orb, c, s, cp, s); k2 = NR.kappa_closed_form(orb, c, s, cp)
        assert abs(k1 - k2) < 1e-10 * (1 + abs(k1)), (k1, k2)
    F = NR.current_F(orb, np.array([0.5, 0.0]), 0.7)
    assert F[1] > 1e-3 and abs(F[0]) < 1e-12
    assert abs(NR.kappa(orb, np.array([0.5, 0.2]), 0.7, np.array([0.5, -0.2]), 0.7)) > 1e-5
    xr_s, xd_s = _segment_orbit()
    seg = A.Orbit(xr_s, xd_s, PER, N=4000, sub=1)
    for _ in range(20):
        c, cp = rng.uniform(-1.5, 1.5, 2), rng.uniform(-1.5, 1.5, 2); s = rng.uniform(0.4, 1.2)
        assert abs(NR.kappa(seg, c, s, cp, s)) < 1e-12


def test_perfect_fit_n2_leaves_only_the_circulation_kernel():
    """n = 2 with f := f̂_μ on all of ℝ² (so the residual AND its Jacobian vanish on the orbit): B(ξ, ξ') = 2κ(ξ, ξ') I, so
    A(ϑ, ϑ') = 2(aᵀa')κ(ξ, ξ') ≠ 0 on the ellipse, and A ≡ 0 on the segment orbit although the fit is the same."""
    net, Th, rng = _setup(2)
    f, J, xr, xd = _plant_2d()
    orb = A.Orbit(xr, xd, PER, N=1500, sub=1)
    fm = A.FrozenMeasure(net, Th)
    found = 0.0
    for _ in range(10):
        th, thp = rng.normal(0, 1.0, 5), rng.normal(0, 1.0, 5)
        a, c, s = net.physical(th[None]); ap, cp, sp = net.physical(thp[None])
        B = NR.B_matrix(fm.f_hat, fm.jac, fm, orb, c[0], s[0], cp[0], sp[0])
        k = NR.kappa(orb, c[0], s[0], cp[0], sp[0])
        assert np.abs(B - 2 * k * np.eye(2)).max() < 1e-13 * (1 + abs(k))
        found = max(found, abs(NR.A_kernel(fm.f_hat, fm.jac, fm, orb, th, thp)))
    assert found > 1e-4
    xr_s, xd_s = _segment_orbit()
    seg = A.Orbit(xr_s, xd_s, PER, N=1500, sub=1)
    fm_s = A.FrozenMeasure(net, Th)
    worst = max(abs(NR.A_kernel(fm_s.f_hat, fm_s.jac, fm_s, seg, rng.normal(0, 1.0, 5), rng.normal(0, 1.0, 5))) for _ in range(10))
    assert worst < 1e-13, worst


def test_perfect_fit_on_the_orbit_alone_keeps_the_curl_term():
    """
    Theorem I(iv), corrected after the round-5 review: r_μ = 0 on the orbit constrains only the tangential derivative of
    the residual.  With f := f̂_μ + h e₁ and h = (x₁/1.2)² + (x₂/0.8)² − 1 (zero on the ellipse, normal derivative ≠ 0),
    B_μ − 2κI = −∫φ_ξφ_ξ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) ≠ 0, so A_μ ≠ 2K although sup_𝒪|r_μ| = 0; the generic second-variation
    code agrees with the bilinear form.
    """
    net, Th, rng = _setup(2)
    _, _, xr, xd = _plant_2d()
    orb = A.Orbit(xr, xd, PER, N=2000, sub=1)
    fm = A.FrozenMeasure(net, Th)
    h = lambda X: (X[:, 0] / 1.2) ** 2 + (X[:, 1] / 0.8) ** 2 - 1.0
    f = lambda X: fm.f_hat(X) + np.stack([h(X), np.zeros(len(X))], 1)

    def J(X):
        out = fm.jac(X).copy()
        out[:, 0, 0] += 2 * X[:, 0] / 1.2 ** 2
        out[:, 0, 1] += 2 * X[:, 1] / 0.8 ** 2
        return out
    assert np.max(np.abs(f(orb.xm) - fm.f_hat(orb.xm))) < 1e-12                        # perfect fit on the orbit
    phi_fn = lambda c, s: np.exp(-np.sum((orb.xm - c[None, :]) ** 2, axis=1) / (2 * s ** 2))
    worst_gap, found = 0.0, 0.0
    for _ in range(10):
        c, cp = rng.uniform(-1.5, 1.5, 2), rng.uniform(-1.5, 1.5, 2); s, sp = rng.uniform(0.4, 1.2, 2)
        B = NR.B_matrix(f, J, fm, orb, c, s, cp, sp)
        k = NR.kappa(orb, c, s, cp, sp)
        d2h = 2 * orb.xm[:, 1] / 0.8 ** 2
        coef = -float(np.mean(phi_fn(c, s) * phi_fn(cp, sp) * d2h))
        pred = coef * np.array([[0.0, 1.0], [-1.0, 0.0]])
        worst_gap = max(worst_gap, np.abs(B - 2 * k * np.eye(2) - pred).max())
        found = max(found, abs(coef))
    assert worst_gap < 1e-12 and found > 1e-3, (worst_gap, found)
    TE = rng.normal(0, 1.0, (3, 5))
    Aref = NR.A_from_second_variation(f, J, fm, orb, TE)
    for i in range(3):
        for j in range(3):
            assert abs(NR.A_kernel(f, J, fm, orb, TE[i], TE[j]) - Aref[i, j]) < 1e-10 * (1 + np.abs(Aref).max())


def test_positivity_witness_n1():
    """f ≡ 0, μ = δ_ϑ̄: A(ϑ₀, ϑ₁) = a₀a₁ā((c₀ − c₁)/s²)·∫φ_c̄φ_c₀φ_c₁ dρ with a strictly positive integral; the generic
    code agrees with the closed form."""
    net, _, rng = _setup(1)
    _, _, xr, xd = _plant_1d()
    orb = A.Orbit(xr, xd, PER, N=1500, sub=1)
    s = 0.7; abar, cbar, a0, c0, a1, c1 = 1.2, 0.3, 0.8, -0.5, -0.6, 0.9
    thbar = net.unconstrained(np.array([[abar]]), np.array([[cbar]]), np.array([s]))
    fm = A.FrozenMeasure(net, thbar)
    f0 = lambda X: np.zeros_like(X); J0 = lambda X: np.zeros((len(X), 1, 1))
    th0 = net.unconstrained(np.array([[a0]]), np.array([[c0]]), np.array([s]))[0]
    th1 = net.unconstrained(np.array([[a1]]), np.array([[c1]]), np.array([s]))[0]
    val = NR.A_kernel(f0, J0, fm, orb, th0, th1)
    cf, I = NR.positivity_witness_1d(orb, a0, c0, a1, c1, abar, cbar, s)
    assert I > 0 and abs(val - cf) < 1e-12 * (1 + abs(cf)) and abs(cf) > 1e-3
