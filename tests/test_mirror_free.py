"""Propositions D‴.45–D‴.47: the free mass–amplitude law as a mirror flow, its leak, its branches (src/core/mirror_free.py)."""
import os
import sys
from decimal import Decimal, getcontext

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import mirror_free as MF


def _instance(seed=3, N=5, d=12):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(d); idx = rng.choice(d, 2, replace=False); ws[idx] = rng.choice([-1, 1], 2) * rng.uniform(0.15, 0.3, 2)
    return Phi, Phi @ ws


def test_eta_closed_form_three_cases():
    """η against the quadrature of dα/(Gβa′) for q > c (also q ≫ c), |q| < c, on a branch (q < −c, from a reference point) on both
    sides, and on the critical curve q = −c."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    f = lambda x, q: 1.0 / (L.G * L.beta(x, q) * L.amp(x)[1])
    for q in (3.0, 0.5, 1e4 * L.c, 0.03, 0.0, -0.02):
        for al in (-2.5, -0.4, 0.3, 1.7, 4.0):
            ref = quad(lambda x: f(x, q), 0.0, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
            assert abs(L.eta(al, q)[0] - ref) <= 1e-9 * max(1.0, abs(ref))
    q = -0.05                                                     # a branch: β = 0 at ±α₀
    a0 = float(L.branch_point(q))
    assert a0 > 0 and abs(L.beta(a0, q)) < 1e-12
    for s in (1.0, -1.0):
        for al in (a0 + 0.2, a0 + 1.0, a0 + 3.0):
            ref = quad(lambda x: f(x, q), s * (a0 + 0.1), s * al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
            assert abs((L.eta(s * al, q) - L.eta(s * (a0 + 0.1), q))[0] - ref) <= 1e-9 * max(1.0, abs(ref))
    assert np.isnan(L.eta(0.5 * a0, q)[0])                        # inside the gap: no curve
    qc = -L.c                                                     # the critical curve: η = (α − a_max coth(α/a_max))/(2Gc)
    for al in (1e-3, 0.4, 3.0, -0.7):
        assert abs(L.eta(al, qc)[0] - (al - L.A / np.tanh(al / L.A)) / (2 * L.G * L.c)) < 1e-12 * max(1.0, abs(L.eta(al, qc)[0]))


def test_eta_against_a_60_digit_reference():
    """q > c: η = (α − a_maxκ artanh(κ tanh(α/a_max)))/(2Gc) against 60-digit decimal arithmetic, for κ near 1 and |t| near 1."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    getcontext().prec = 60
    A_, G_, c_ = (Decimal(repr(v)) for v in (L.A, L.G, L.c))

    def ref(al, q):
        al, q = Decimal(repr(al)), Decimal(repr(q))
        e2 = (2 * al / A_).exp(); t = (e2 - 1) / (e2 + 1)
        k = ((q - c_) / (q + c_)).sqrt()
        return float((al - A_ * k * ((1 + k * t) / (1 - k * t)).ln() / 2) / (2 * G_ * c_))
    for q in (0.026, 0.3, 1e4, 1e8, 1e14):
        for al in (1e-8, 0.3, 0.6, 2.0, 25.0, 200.0, -3.0):
            r = ref(al, q)
            assert abs(L.eta(al, q)[0] - r) <= 4e-15 * abs(r)


def test_inverse_round_trips_and_range():
    """alpha_of_eta inverts η at the extremes and on both sides of a branch; it raises instead of clipping."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    for q, als in ((1e4 * L.c, (1e-6, 0.7, 30.0, -12.0)), (0.01, (0.2, -5.0, 40.0))):
        for al in als:
            assert abs(L.alpha_of_eta(L.eta(al, q), q)[0] - al) <= 1e-12 * max(1.0, abs(al))
    q = -2 * L.c
    a0 = float(L.branch_point(q))
    for s in (1.0, -1.0):
        for dl in (1e-12, 1e-3, 2.0):
            al = s * (a0 + dl)
            back = L.alpha_of_eta(L.eta(al, q), q, sgn=s)[0]
            assert abs(back - al) <= 1e-14 * abs(al) + 1e-3 * dl
    with pytest.raises(ValueError):
        L.alpha_of_eta(1e9, 0.5)                                  # beyond |α| = 340
    with pytest.raises(ValueError):
        L.alpha_of_eta(-1e5, q)                                   # closer to the branch point than double precision
    with pytest.raises(ValueError):
        L.alpha_of_w(0.0, q)                                      # a branch's end point


def test_potential_conjugate_and_metric():
    """ψ′ = η, ψ*(η) = (β − β⁰)/(γm) with (ψ*)′ = w, and dw/dη = mD (finite differences)."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    for q in (0.97, 0.01, -0.02):
        for al in (-1.3, 0.2, 2.2):
            hh = 1e-5
            dw = (L.w(al + hh, q) - L.w(al - hh, q)) / (2 * hh)
            de = (L.eta(al + hh, q) - L.eta(al - hh, q))[0] / (2 * hh)
            dpsi = (L.psi(al + hh, q) - L.psi(al - hh, q))[0] / (2 * hh)
            dps = (L.psi_star(al + hh, q) - L.psi_star(al - hh, q)) / (2 * hh)
            assert abs(dpsi / dw - L.eta(al, q)[0]) < 1e-6 * max(1, abs(L.eta(al, q)[0]))
            assert abs(dps / de - L.w(al, q)) < 1e-6
            assert abs(dw / de - L.m * L.relevance(al, q)) < 1e-6 * L.m * L.relevance(al, q)


def test_conservation_for_any_drive_and_conserved_law_shift():
    """Q constant along the loop with a smooth node-correlated noise; for the conserved law dQ_i/dt = −γβ_iλ."""
    Phi, y = _instance()
    L = MF.FreeLaw()
    rng = np.random.default_rng(1)
    om = np.exp(rng.uniform(-6, 1, 12)); ph = rng.uniform(0, 6.28, 12); Mx = 0.3 * rng.standard_normal((5, 12))
    xi = lambda t: Mx @ np.sin(om * t + ph)
    d = Phi.shape[1]
    z0 = np.r_[0.05 * rng.standard_normal(d), np.log(np.full(d, 0.3))]
    sol = solve_ivp(L.rhs(Phi, y, xi), (0, 300), z0, method='DOP853', rtol=1e-11, atol=1e-13)
    q0 = L.Q(z0[:d], np.exp(z0[d:]))
    assert np.abs(L.Q(sol.y[:d, -1], np.exp(sol.y[d:, -1])) - q0).max() < 1e-8
    F = L.rhs(Phi, y, xi, conserved=True)
    z = sol.y[:, -1]
    dz = F(1.0, z)
    al, b = z[:d], np.exp(z[d:])
    a, da = L.amp(al)
    dQ = b * dz[d:] - L.c * (2 / L.A) * np.sinh(2 * al / L.A) * dz[:d]
    r = y - Phi @ (b * a / L.m) + xi(1.0); g = Phi.T @ r / Phi.shape[0]
    lam = b @ (a * g) / b.sum()
    assert np.allclose(dQ, -L.gamma * b * lam, rtol=1e-10, atol=1e-14)


def test_noiseless_limit_is_the_bregman_projection():
    Phi, y = _instance()
    L = MF.FreeLaw()
    d = Phi.shape[1]
    al0 = 0.05 * np.random.default_rng(5).standard_normal(d); b0 = np.full(d, 1.0)
    q = L.Q(al0, b0)
    w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q))
    assert res < 1e-12
    sol = solve_ivp(L.rhs(Phi, y), (0, 4e4), np.r_[al0, np.log(b0)], method='LSODA', rtol=1e-11, atol=1e-13)
    wT = np.exp(sol.y[d:, -1]) * L.amp(sol.y[:d, -1])[0] / L.m
    assert np.abs(wT - w).max() < 1e-9


def test_rich_projection_is_weighted_basis_pursuit():
    """From α = 0 and masses β⁰ → 0 the projection converges (residual recorded) and approaches basis pursuit weighted by 1/√β⁰ —
    also at masses where the earlier solver diverged."""
    Phi, y = _instance(seed=7, N=6, d=14)
    L = MF.FreeLaw()
    d = Phi.shape[1]
    b0 = np.full(d, 1.0); b0[::3] = 4.0                           # unequal masses: weights 1/√β⁰
    wbp = MF.weighted_bp(Phi, y, 1 / np.sqrt(b0))
    dists = []
    for s in (1e-5, 1e-8, 1e-10):
        q = L.Q(np.zeros(d), s * b0)
        w, al, it, res = L.bregman_projection(Phi, y, q, np.zeros(d))
        assert res < 1e-9 * (1 + np.linalg.norm(y))
        dists.append(np.abs(w - wbp).max())
    assert dists[-1] < 1e-4 and dists[-1] < dists[0]


def test_euler_change_of_Q():
    """The Euler step's change of Q: second order in h with the closed-form coefficient."""
    L = MF.FreeLaw()
    rng = np.random.default_rng(2)
    al = rng.uniform(-2, 2, 8); b = rng.uniform(0.1, 3, 8); g = rng.standard_normal(8)
    for h in (1e-2, 1e-3):
        ex = L.dQ_exact(al, b, g, h); pr = L.dQ_euler(al, b, g ** 2, h)
        assert np.all(np.abs(ex / pr - 1) < 50 * h)


def _real_setup(rng):
    L = MF.FreeLaw(m=6.0, G=1.0, gamma=0.5, a_max=1.0)
    d = 6
    cen, wid = np.linspace(-1.2, 1.2, d), np.full(d, 0.4)
    al0, b0 = 0.1 * rng.standard_normal(d), rng.uniform(0.5, 1.5, d)
    q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, d)
    return L, d, cen, wid, al0, b0, q, wstar, L.alpha_of_w(wstar, q)


def test_real_loop_lyapunov():
    """V̇ = −k_ee² from the vector field at random states (the free law, a realisable target); V reads the state's own mass, so a
    wrong mass law (here without the amplitude a) is seen."""
    rng = np.random.default_rng(4)
    L, d, cen, wid, al0, b0, q, wstar, al_star = _real_setup(rng)
    xr, xrd = (lambda t: np.sin(t)), (lambda t: np.cos(t))
    ke = 5.0
    F = L.real_rhs(cen, wid, wstar, xr, xrd, ke)
    worst_mut = 0.0
    for _ in range(5):
        al = al0 + 0.3 * rng.standard_normal(d)
        z = np.r_[rng.uniform(-1, 1), al, np.log(L.beta(al, q))]
        t = rng.uniform(0, 6)
        dz = F(t, z)
        hh = 1e-6
        V = lambda zz, s: L.real_lyapunov(zz + s * dz_, t + s, xr, al_star, q)
        dz_ = dz
        dV = (V(z, hh) - V(z, -hh)) / (2 * hh)
        e = z[0] - xr(t)
        assert abs(dV + ke * e ** 2) < 1e-6 * max(1.0, ke * e ** 2)
        phi = np.exp(-(z[0] - cen) ** 2 / (2 * wid ** 2))
        dz_ = dz.copy(); dz_[1 + d:] = L.gamma * e * phi       # the mass law without a: (log β)˙ = γeφ
        dVm = (V(z, hh) - V(z, -hh)) / (2 * hh)
        worst_mut = max(worst_mut, abs(dVm + ke * e ** 2) / max(1.0, ke * e ** 2))
    assert worst_mut > 1e-3


def test_mirror_step_keeps_Q_and_moves_eta_by_hg():
    """The mirror update on free curves and on both sides of branches, with steps large enough to need many Newton iterations."""
    L = MF.FreeLaw()
    rng = np.random.default_rng(6)
    al = rng.uniform(-2, 2, (4, 7)); q = rng.uniform(0.05, 2.0, 7)
    g = rng.standard_normal((4, 7))
    al1 = L.mirror_step(al, q, g, 0.3)
    assert np.allclose(L.eta(al1, np.broadcast_to(q, al.shape)) - L.eta(al, np.broadcast_to(q, al.shape)), 0.3 * g, atol=1e-11)
    assert np.allclose(L.Q(al1, L.beta(al1, q)), np.broadcast_to(q, al.shape), atol=1e-13)
    qb = np.array([-0.05, -0.05, -0.1, -0.1])                    # branches
    a0 = L.branch_point(qb)
    alb = np.array([1.0, -1.0, -1.0, 1.0]) * (a0 + np.array([0.3, 1e-3, 2.0, 0.05]))
    gb = np.array([-40.0, 25.0, 3.0, 80.0])
    al2 = L.mirror_step(alb, qb, gb, 1.0)
    assert np.all(np.sign(al2) == np.sign(alb)) and np.all(np.abs(al2) > a0)
    assert np.allclose(L.eta(al2, qb) - L.eta(alb, qb), gb, rtol=1e-10, atol=1e-9)


def test_leak_lyapunov_identity_noise_and_equilibrium_limits():
    """Proposition D‴.46: V̇ = −k_ee² + ed − n(f − f̂) − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] with the leak, a disturbance d and a
    measurement noise n in the adaptation; the averaged equilibrium (k_e ≠ 1, a bias d̄) is the ridge at large starting masses, w* + H⁻¹d̄
    as σ → 0, and near the LASSO (same support) at small ones; the solver's residual is returned."""
    rng = np.random.default_rng(4)
    L, d, cen, wid, al0, b0, q, wstar, al_star = _real_setup(rng)
    xr, xrd = (lambda t: np.sin(t)), (lambda t: np.cos(t))
    dist = lambda t: 0.2 * np.sin(3.1 * t)
    noise = lambda t: 0.05 * np.cos(7.3 * t + 0.4)
    ftrue = lambda x: np.exp(-(x - cen) ** 2 / (2 * wid ** 2)) @ wstar
    ke, sig = 5.0, 0.3
    F = L.real_rhs_leak(cen, wid, ftrue, xr, xrd, ke, q, L.eta(al0, q), sig, dist=dist, noise=noise)
    B = lambda a1, a2: L.bregman_div(a1, a2, q)
    for _ in range(3):
        al = al0 + 0.3 * rng.standard_normal(d)
        z = np.r_[rng.uniform(-1, 1), al, np.log(L.beta(al, q))]
        t = rng.uniform(0, 6)
        dz = F(t, z); hh = 1e-6
        V = lambda s: L.real_lyapunov(z + s * dz, t + s, xr, al_star, q)
        dV = (V(hh) - V(-hh)) / (2 * hh)
        e = z[0] - xr(t)
        phi = np.exp(-(z[0] - cen) ** 2 / (2 * wid ** 2))
        ff = phi @ wstar - phi @ L.w(al, q)
        pred = -ke * e ** 2 + e * dist(t) - noise(t) * ff - sig * (B(al_star, al) + B(al, al0) - B(al_star, al0))
        assert abs(dV - pred) < 1e-6 * max(1.0, abs(pred))
        assert abs(dz[1:1 + d] @ (-L.c * 2 * np.sinh(2 * al / L.A) / L.A) + dz[1 + d:] @ L.beta(al, q)) < 1e-10   # Q kept
    Lz = MF.FreeLaw()
    A = rng.standard_normal((30, 10)); H = A.T @ A / 30
    ws = np.zeros(10); ws[[1, 6]] = [0.2, -0.25]
    dbar = 0.02 * rng.standard_normal(10)
    ke = 3.0
    b0 = 100.0
    q = Lz.Q(np.zeros(10), np.full(10, b0))
    w, _, _, res = Lz.leak_equilibrium(H, ws, q, np.zeros(10), 1e-2, ke, dphi=dbar)
    assert res < 1e-10
    ridge = np.linalg.solve(H + ke * 1e-2 * Lz.m / (Lz.G * b0 ** 2) * np.eye(10), H @ ws + dbar)
    assert np.abs(w - ridge).max() < 1e-6 * np.abs(ridge).max()
    w, _, _, res = Lz.leak_equilibrium(H, ws, q, np.zeros(10), 1e-9, ke, dphi=dbar)
    assert res < 1e-10 and np.abs(w - (ws + np.linalg.solve(H, dbar))).max() < 1e-6
    b0 = 1e-8
    q = Lz.Q(np.zeros(10), np.full(10, b0))
    lam1 = np.pi * Lz.A * np.sqrt((2 * Lz.c - b0) / b0) / (4 * Lz.G * Lz.c)
    w, _, _, res = Lz.leak_equilibrium(H, ws, q, np.zeros(10), 0.02 / lam1, 1.0)
    assert res < 1e-9
    ls = MF.lasso(H, ws, np.full(10, 0.02), 1.0)
    assert np.abs(w - ls).max() < 1e-3 * np.abs(ls).max()
    assert np.array_equal(np.abs(w) > 1e-3 * np.abs(w).max(), np.abs(ls) > 1e-12)


def test_leak_comparison_bound_under_a_worst_case_disturbance():
    """V(t) ≤ V(0)e^{−λt} + (C/λ)(1 − e^{−λt}), λ = min(k_e, σ), C = d̄²/(2k_e) + σB(w*, w_p), along trajectories with the disturbance
    d̄·tanh(e/0.01); the trajectory comes within a factor 2 of the ultimate bound C/λ (the bound binds)."""
    rng = np.random.default_rng(4)
    L, d, cen, wid, al0, b0, q, wstar, al_star = _real_setup(rng)
    xr, xrd = (lambda t: np.sin(t)), (lambda t: np.cos(t))
    ftrue = lambda x: np.exp(-(x - cen) ** 2 / (2 * wid ** 2)) @ wstar
    for ke, sig, dbar in ((1.0, 0.5, 0.3), (0.5, 2.0, 0.5)):
        Fb = L.real_rhs_leak(cen, wid, ftrue, xr, xrd, ke, q, L.eta(al0, q), sig)

        def F(t, z):
            dz = Fb(t, z)
            dz[0] += dbar * np.tanh((z[0] - xr(t)) / 0.01)
            return dz
        sol = solve_ivp(F, (0, 40), np.r_[xr(0) + 2.0, al0, np.log(b0)], method='LSODA', rtol=1e-9, atol=1e-11,
                        t_eval=np.linspace(0, 40, 801))
        V = np.array([L.real_lyapunov(sol.y[:, k], t, xr, al_star, q) for k, t in enumerate(sol.t)])
        lam = min(ke, sig)
        C = dbar ** 2 / (2 * ke) + sig * L.bregman_div(al_star, al0, q)
        bound = V[0] * np.exp(-lam * sol.t) + (C / lam) * (1 - np.exp(-lam * sol.t))
        assert np.all(V <= bound * (1 + 1e-6))
        assert V[-1] > 0.5 * C / lam


def test_filtered_gram():
    """H_L = k_e⟨φL[φ]ᵀ⟩ on a scalar periodic orbit: H − ⟨φ̇φ̇ᵀ⟩/k_e² + O(k_e⁻³) (the first-order lag ⟨φφ̇ᵀ⟩ vanishes), its asymmetry
    O(k_e⁻³) and nonzero when ⟨ẋ³(φ_i′φ_j″ − φ_j′φ_i″)⟩ ≠ 0 (x = sin t + 0.4 sin 2t); d̄_L → ⟨φd⟩."""
    n = 4096
    ts = np.arange(n) * 2 * np.pi / n
    x = np.sin(ts) + 0.4 * np.sin(2 * ts); xd = np.cos(ts) + 0.8 * np.cos(2 * ts)
    cen, wid = np.linspace(-1.2, 1.2, 5), 0.4
    P = np.exp(-(x[:, None] - cen) ** 2 / (2 * wid ** 2))
    Pd = P * (-(x[:, None] - cen) / wid ** 2) * xd[:, None]
    H = P.T @ P / n; K = Pd.T @ Pd / n
    dist = 0.3 * np.sin(3 * ts) + 0.1
    err, asym = [], []
    for ke in (20.0, 40.0):
        HL, dL = MF.filtered_gram(P, 2 * np.pi, ke, dist)
        err.append(np.abs(HL - (H - K / ke ** 2)).max()); asym.append(np.abs(HL - HL.T).max())
    assert 5 < err[0] / err[1] < 11 and 5 < asym[0] / asym[1] < 11 and asym[1] > 1e-9
    HL, dL = MF.filtered_gram(P, 2 * np.pi, 1e6, dist)
    assert np.abs(HL - H).max() < 1e-9 and np.abs(dL - P.T @ dist / n).max() < 1e-6


# --- Proposition D‴.47: the branches ------------------------------------------------------------------------------------------------
def test_particle_set_branches_identities():
    """The branch criterion; δ(0) solved from β(0) (masses far below c's precision kept); η's inverse to any depth; the temperature
    1/(γa₀) at a branch point; B(w̄, w) → β/(γm) as w̄ → 0; a ± pair's W(Δ) = 2u sinh(γa₀Δ) at leading order."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    P = MF.ParticleSet(L, np.array([0.3, -0.3, 0.0, 0.5, 0.001]), np.array([1e-3, 1e-3, 1e-3, 1e-200, 1e-3]))
    assert list(P.br) == [True, True, False, True, False]
    assert abs(P.beta(P.x0)[3] / 1e-200 - 1) < 1e-12 and abs(P.beta(P.x0)[0] / 1e-3 - 1) < 1e-12
    for u in (-5000.0, -60.0, 1.0):
        x = P.x0.copy(); x[0] = x[1] = u
        assert np.abs(P.x_of_eta(P.eta(x))[:2] - u).max() < 1e-9 * max(1.0, abs(u))
    xb = P.x0.copy(); xb[0] = -60.0
    a0 = L.A * np.tanh(P.al0[0] / L.A)
    assert abs(P.temperature(xb)[0] * L.gamma * a0 - 1) < 1e-12
    # B(0, w) = β/(γm) on a branch, against the general formula's limit
    Pc = MF.ParticleSet(L, np.array([0.4]), np.array([1e-3]))
    xw = np.array([np.log(0.05)])
    for dbar in (1e-30, 1e-100):
        xbar = np.array([np.log(dbar)])
        wb = Pc.w(xbar)
        gen = float(wb[0] * (Pc.eta(xbar)[0] - Pc.eta(xw)[0]) - (Pc.beta(xbar)[0] - Pc.beta(xw)[0]) / (L.gamma * L.m))
        assert abs(gen - Pc.beta(xw)[0] / (L.gamma * L.m)) < 1e-14
    assert abs(Pc.bdiv(np.zeros(1), xw, xw) - Pc.beta(xw)[0] / (L.gamma * L.m)) < 1e-16
    Pp = MF.ParticleSet(L, np.array([0.3, -0.3]), np.array([1e-7, 1e-7]))
    e0 = Pp.eta(Pp.x0); u = Pp.w(Pp.x0)[0]; a0 = L.A * np.tanh(Pp.al0[0] / L.A)
    for D in (1.0, 10.0):
        W = Pp.w(Pp.x_of_eta(e0 + D)).sum()
        assert abs(W / (2 * u * np.sinh(L.gamma * a0 * D)) - 1) < 1e-4


def test_branch_limit_is_the_projection_onto_the_sign_constrained_fits():
    """Proposition D‴.47(ii): from a start whose frozen signs cannot fit the labels, the loop converges to the Bregman projection of
    its start onto the sign-constrained least-squares set S (computed independently), the coordinates vanishing on S going extinct."""
    rng = np.random.default_rng(1)
    N, d = 6, 12
    Phi = rng.standard_normal((N, d))
    ws = np.zeros(d); ws[[0, 3, 5, 8]] = [0.9, -0.7, 0.6, -0.5]; y = Phi @ ws
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    P = MF.ParticleSet(L, rng.uniform(0.2, 0.4, d), np.full(d, 1e-3))      # all on positive branches
    assert P.br.all()
    yh, Z, _ = MF.sign_ls_set(Phi, y, P.s)
    assert np.linalg.norm(y - yh) > 1e-3 and Z.size > 0
    keep = np.setdiff1d(np.arange(d), Z)
    wbar, it, res = MF.projection(P, Phi, yh, P.eta(P.x0), keep)
    assert res < 1e-11
    sol = solve_ivp(P.rhs(Phi, y), (0, 3e7), P.x0, method='LSODA', rtol=1e-10, atol=1e-12)
    assert np.abs(P.w(sol.y[:, -1]) - wbar).max() < 1e-9
    assert P.beta(sol.y[:, -1])[Z].max() < 1e-30


def test_exact_newton_fallback_at_extreme_masses():
    """Proposition D‴.47(iii) at ε = 10⁻³²: with particles near and far from extinction the float Hessian loses the small curvatures and
    the float Newton stalls; the decimal fallback converges and agrees with the leading order from the vertex's certificate, and leaves a
    converged float iteration unchanged."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((8, 20))
    ws = np.zeros(20); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    sg = rng.choice([-1.0, 1.0], 20); sg[[1, 5, 9]] = np.sign(ws[[1, 5, 9]])
    amp0 = rng.uniform(0.2, 0.6, 20)
    a0 = L.A * np.tanh(amp0 / L.A)
    wl = MF.weighted_bp(Phi, y, 1 / a0, signs=sg)
    S = np.flatnonzero(np.abs(wl) > 1e-9)
    assert S.size == 8
    mu, _ = MF.certificate_slacks(Phi, wl, S, 1 / (L.gamma * a0), frozen=sg)
    off = np.setdiff1d(np.arange(20), S)
    assert np.abs(mu[S]).max() < 1e-12 and mu[off].min() > 0.05
    P = MF.ParticleSet(L, sg * amp0, np.full(20, 1e-4))
    w_f = MF.projection(P, Phi, y, P.eta(P.x0))[0]
    assert np.array_equal(w_f, MF.projection(P, Phi, y, P.eta(P.x0), exact=True)[0])
    P = MF.ParticleSet(L, sg * amp0, np.full(20, 1e-32))
    w, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
    assert res < 1e-11 and it < 200
    w1 = MF.rich_leading_order(P, Phi, y, wl, S)
    assert np.abs(w - w1).max() < 1e-2 * np.abs(w - wl).max()


def test_rich_limit_leading_order_from_alpha_zero():
    """Proposition D‴.45(iii): from α = 0 at the masses εr_i the off-support coefficients are a_max(εr_j)^{3/2}cos(πμ_j/2)/(m(2c)^{1/2}
    sin³(πμ_j/2)) to leading order (μ_j the weighted basis pursuit certificate's slacks): the approach is O(ε^{3/2}), the leading order
    from the certificate exact to higher order."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    rng = np.random.default_rng(0)
    N, d = 5, 12
    Phi = rng.standard_normal((N, d))
    ws = np.zeros(d); ws[rng.choice(d, 3, replace=False)] = rng.choice([-1, 1], 3) * rng.uniform(0.15, 0.3, 3); y = Phi @ ws
    r = np.exp(rng.uniform(np.log(0.25), np.log(4.0), d))
    wl = MF.weighted_bp(Phi, y, r ** -0.5)
    S = np.flatnonzero(np.abs(wl) > 1e-10)
    assert S.size == N
    mu, xi = MF.certificate_slacks(Phi, wl, S, r ** -0.5)
    off = np.setdiff1d(np.arange(d), S)
    dist = {}
    for eps, tol_cf in ((1e-5, 0.15), (1e-7, 0.015)):
        P = MF.ParticleSet(L, np.zeros(d), eps * r)
        w, it, res = MF.projection(P, Phi, y, P.eta(P.x0))
        assert res < 1e-12
        dist[eps] = np.abs(w - wl).max()
        assert np.abs(w - MF.rich_leading_order(P, Phi, y, wl, S)).max() < 1e-4 * dist[eps] + 1e-13          # + rounding
        cf = np.sign(Phi.T @ xi) * L.core_tail(eps * r, mu)
        assert np.max(np.abs(w[off] - cf[off]) / np.abs(cf[off])) < tol_cf
    assert 800 < dist[1e-5] / dist[1e-7] < 1300                                  # ε^{3/2}: 1000


def test_iss_path_and_the_rich_transit():
    """Proposition D‴.48(i): the inverse scale space path ends at weighted basis pursuit (one-sided: the sign-constrained one), and the
    loop from α = 0 at the masses εr_i, in the time t/Λ(ε), follows it: u = η/Λ within 0.02 of the path at ε = 10⁻⁶."""
    rng = np.random.default_rng(3)
    N, d = 6, 14
    Phi = rng.standard_normal((N, d))
    ws = np.zeros(d); ws[[1, 4, 9]] = [0.8, -0.5, 0.4]; y = Phi @ ws
    r = np.exp(rng.uniform(np.log(0.25), np.log(4.0), d))
    rho = r ** -0.5
    path = MF.iss_path(Phi, y, rho)
    assert np.abs(path['w_end'] - MF.weighted_bp(Phi, y, rho)).max() < 1e-12
    s = rng.choice([-1.0, 1.0], d); s[[1, 4, 9]] = np.sign(ws[[1, 4, 9]])
    assert np.abs(MF.iss_path(Phi, y, rho, signs=s)['w_end'] - MF.weighted_bp(Phi, y, rho, signs=s)).max() < 1e-12
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    eps = 1e-6
    Lam = (np.pi * L.A / (4 * L.G)) * np.sqrt(2 / (L.c * eps))
    taus = np.linspace(0, 1.5 * path['taus'][-1], 301)
    b0 = eps * r; q = L.Q(np.zeros(d), b0)
    sol = solve_ivp(L.rhs(Phi, y), (0, Lam * taus[-1]), np.r_[np.zeros(d), np.log(b0)], method='LSODA', rtol=1e-10, atol=1e-13,
                    t_eval=Lam * taus)
    err = max(np.abs(L.eta(sol.y[:d, k], q) / Lam - MF.iss_eval(path, t, Phi, y)[0]).max() for k, t in enumerate(taus))
    assert err < 0.02


def _iss_stepping(Phi, y, rho, signs=None, dtau=2e-4, tau_end=10.0):
    """An independent inverse scale space integrator (round 37): forward steps of the dual u in τ, clipped to the box (one-sided:
    s_iu_i ≤ ρ_i), w the (sign-constrained) least-squares fit on the atoms at their walls; events read off w's support."""
    from scipy.optimize import nnls
    N, d = Phi.shape
    one = signs is not None
    u = np.zeros(d)
    supp = np.zeros(d, bool)
    events = []
    for k in range(int(round(tau_end / dtau))):
        sg = np.asarray(signs, float) if one else np.sign(u)
        wall = (sg * u >= rho * (1 - 1e-12)) & (sg != 0)
        w = np.zeros(d)
        idx = np.flatnonzero(wall)
        if idx.size:
            w[idx] = sg[idx] * nnls(Phi[:, idx] * sg[idx], y, maxiter=50 * d)[0]
        nz = np.abs(w) > 1e-12
        for i in np.flatnonzero(nz & ~supp):
            events.append((k * dtau, int(i), 'enter'))
        for i in np.flatnonzero(supp & ~nz):
            events.append((k * dtau, int(i), 'exit'))
        supp = nz
        u = u + dtau * Phi.T @ (y - Phi @ w) / N
        u = np.minimum(u, rho) if one else np.clip(u, -rho, rho)
        if one:
            u = np.where(sg > 0, np.minimum(u, rho), np.maximum(u, -rho))
    return events, w


@pytest.mark.parametrize('seed', [4, 5])
def test_iss_path_events_against_an_independent_integrator(seed):
    """Proposition D‴.48(i), round 37: on designs with exits, `iss_path`'s entries and exits (atom, kind, time) are those of a forward-stepping
    integrator of the dual clipped to the box (the times to its O(dτ) per event), and both end at weighted basis pursuit."""
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_mirror_transit as X
    Phi, ws, b = X.design(seed); y = Phi @ ws
    rho = b ** -0.5
    path = MF.iss_path(Phi, y, rho)
    assert any(k == 'exit' for t, i, k in path['events'])
    dtau = 2e-4
    ev, w = _iss_stepping(Phi, y, rho, dtau=dtau, tau_end=path['taus'][-1] * 1.05)
    assert [(i, k) for t, i, k in ev] == [(i, k) for t, i, k in path['events']]
    # the clipped Euler's O(dτ) error in u becomes δτ = δu/|du_j/dτ| at an entry — large where the entrant approaches its wall slowly
    # (design 4: |du/dτ| = 0.0015 at τ = 55.5) —, so the times are compared in u: |δτ|·|du_j/dτ| (an exit at an entry's time: the entry's)
    N = Phi.shape[0]
    vel = {}
    for k, t in enumerate(path['taus'][1:]):
        v = Phi.T @ (y - Phi @ path['W'][k]) / N
        for tt, i, kind in path['events']:
            if kind == 'enter' and tt == t:
                vel[t] = abs(v[i])
    for (t1, i1, k1), (t2, i2, k2) in zip(ev, path['events']):
        assert abs(t1 - t2) * vel[t2] < 2 * dtau
    assert np.abs(path['w_end'] - MF.weighted_bp(Phi, y, rho)).max() < 1e-12
    assert np.abs(w - path['w_end']).max() < 1e-6


def test_one_sided_iss_where_the_sign_constraint_binds():
    """Proposition D‴.48(ii), round 37: with frozen signs (Part B's instance) the one-sided path differs from the two-sided one — an atom the
    two-sided path enters at its wrong wall — and matches the independent integrator with one-sided boxes."""
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_mirror_transit as X
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    Phi, ws, b = X.design(X.SEEDS_A[0]); y = Phi @ ws
    rng = np.random.default_rng(11)
    s = rng.choice([-1.0, 1.0], X.DA); s[ws != 0] = np.sign(ws[ws != 0])
    amp = rng.uniform(0.2, 0.6, X.DA)
    rho = 1 / (L.gamma * L.A * np.tanh(amp / L.A))
    one = MF.iss_path(Phi, y, rho, signs=s)
    two = MF.iss_path(Phi, y, rho)
    wrong = [i for t, i, k in two['events'] if k == 'enter' and np.sign(two['W'][-1][i] if two['W'][-1][i] != 0 else 1) != s[i]]
    assert [(i, k) for t, i, k in one['events']] != [(i, k) for t, i, k in two['events']]
    ev, w = _iss_stepping(Phi, y, rho, signs=s, dtau=5e-3, tau_end=one['taus'][-1] * 1.05)
    assert [(i, k) for t, i, k in ev] == [(i, k) for t, i, k in one['events']]
    assert np.abs(one['w_end'] - MF.weighted_bp(Phi, y, rho, signs=s)).max() < 1e-10
    assert np.all(s * one['w_end'] >= 0)


def test_particle_set_free_drive():
    """The free particles' drive (round 37: an a′ mutant survived): α̇ = Gβa′g, β on the curve, against FreeLaw's closed forms, and the
    branch particles' u̇ = sGβa′g/δ against a finite difference of δ along the flow."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    al = np.array([0.0, 0.7, -1.3, 2.4])
    P = MF.ParticleSet(L, al, np.array([0.2, 0.5, 1.0, 3.0]))
    assert not P.br.any()
    g = np.array([0.3, -0.2, 0.5, 0.1])
    a, da = L.amp(al)
    assert np.allclose(P._dstate(P.x0, g), L.G * P.beta(P.x0) * da * g, rtol=1e-13, atol=0)
    assert np.allclose(P._dstate(P.x0, g), L.G * L.beta(al, P.q) * (1 - np.tanh(al / L.A) ** 2) * g, rtol=1e-12, atol=0)
    Pb = MF.ParticleSet(L, np.array([0.4]), np.array([1e-4]))
    assert Pb.br.all()
    x = Pb.x0 + 0.3
    h = 1e-7
    deta = Pb.eta(x + h) - Pb.eta(x - h)
    assert abs(float(Pb._dstate(x, np.array([1.0]))[0]) * float(deta[0]) / (2 * h) - 1) < 1e-6


def test_bregman_selection_with_extinction_on_a_non_singleton_set():
    """Proposition D‴.47(ii)(b), round 36: an inexact sign-constrained fit whose least-squares set S is not a point (each atom duplicated),
    two atoms extinct (Z ≠ ∅), negative branches alive: the flow ends at the Bregman projection onto S — not at another point of S such as
    the equal split of the duplicates — and the extinct masses vanish."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    rng = np.random.default_rng(953)
    N, d = 5, 4
    Phi = rng.standard_normal((N, d)); y = rng.standard_normal(N)
    s = np.sign(np.linalg.lstsq(Phi, y, rcond=None)[0]); s[3] = -s[3]
    Phi2 = np.hstack([Phi, Phi])
    P = MF.ParticleSet(L, np.r_[s, s] * rng.uniform(0.4, 0.8, 2 * d), np.full(2 * d, 1e-2))
    assert P.br.all()
    yh, Z, _ = MF.sign_ls_set(Phi2, y, P.s)
    keep = np.setdiff1d(np.arange(2 * d), Z)
    wbar, it, res = MF.projection(P, Phi2, yh, P.eta(P.x0), keep)
    assert np.linalg.norm(y - yh) > 0.1 and Z.size >= 2 and np.any(wbar[keep] < 0)
    W = wbar[:d] + wbar[d:]
    assert np.abs(wbar - np.r_[W / 2, W / 2]).max() > 0.01
    sol = solve_ivp(P.rhs(Phi2, y), (0, 1e5), P.x0, method='LSODA', rtol=1e-11, atol=1e-13)
    wT = P.w(sol.y[:, -1])
    assert np.abs(wT - wbar).max() < 1e-9 and np.abs(wT[Z]).max() < 1e-15


def test_filtered_gram_is_causal():
    """Round 36: H_L = k_e⟨φL[φ]ᵀ⟩ with L = 1/(s + k_e) the causal filter — for φ = (sin t, cos t), L[cos] = (k cos t + sin t)/(1 + k²), so
    H_12 = k/(2(1 + k²)) = −H_21 and H_11 = k²/(2(1 + k²)) (an anti-causal filter flips the antisymmetric part)."""
    n, k = 512, 3.0
    t = 2 * np.pi * np.arange(n) / n
    H = MF.filtered_gram(np.c_[np.sin(t), np.cos(t)], 2 * np.pi, k)
    assert abs(H[0, 1] - k / (2 * (1 + k * k))) < 1e-12 and abs(H[1, 0] + k / (2 * (1 + k * k))) < 1e-12
    assert abs(H[0, 0] - k * k / (2 * (1 + k * k))) < 1e-12


def test_irrepresentable_and_nonneg_certificate_values():
    """Round 36: closed forms — for the equicorrelated H (ρ off the diagonal), S = {0, 1} with equal signs, 2ρ/(1 + ρ), and with weights
    λ = (1, 2, 4) the weighted value (ρ(1 + 2) /(1 + ρ)·…)/λ₂ computed by hand; the nonnegative certificate's margin 1 on the identity's
    columns and 0 for a column that is the sum of two others."""
    rho = 0.3
    H = np.full((3, 3), rho) + (1 - rho) * np.eye(3)
    assert abs(MF.irrepresentable(H, [0, 1], [1, 1]) - 2 * rho / (1 + rho)) < 1e-14
    lam = np.array([1.0, 2.0, 4.0])
    x = np.linalg.solve(H[:2, :2], lam[:2])
    assert abs(MF.irrepresentable(H, [0, 1], [1, 1], lam) - abs(H[2, :2] @ x) / lam[2]) < 1e-14
    mu, v = MF.nonneg_certificate(np.eye(3), [0])
    assert abs(mu - 1) < 1e-9
    mu2, v2 = MF.nonneg_certificate(np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]), [2])
    assert abs(mu2) < 1e-9


def test_lasso_kkt_at_a_gain():
    """Round 36: the LASSO min (1/(2k_e))(w − w*)ᵀH(w − w*) + Σμ_i|w_i| at k_e = 5 satisfies its KKT conditions."""
    rng = np.random.default_rng(8)
    A = rng.standard_normal((12, 6)); H = A.T @ A / 12
    ws = np.array([0.8, -0.5, 0.0, 0.3, 0.0, 0.0]); mu = np.full(6, 0.02); ke = 5.0
    w = MF.lasso(H, ws, mu, ke)
    g = H @ (ws - w) / ke
    on = np.abs(w) > 1e-12
    assert np.abs(g[on] - mu[on] * np.sign(w[on])).max() < 1e-10 and np.all(np.abs(g[~on]) <= mu[~on] + 1e-10)
    assert on.sum() < 6 and on.sum() >= 2
