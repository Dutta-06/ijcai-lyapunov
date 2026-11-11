"""
Proposition D‴.34 (PARADIGM.md §5d): mass transfer — the replicator weight law with the cohesion leak.

Checked here: (i) the weighted loop is the gradient flow of J(Θ, β) = mL̄(μ_{Θ,β}) + (σ₀m/2)∫∫W dμ_Θdμ_Θ — the positions' drift −G∂_ΘJ, the
weights' γβ_i(u_i − λ) with u_i = −∂_{β_i}J — and with β ≡ 1, γ = 0 it is Proposition D‴.32's drift; (ii) in the real loop at fixed positions
the relative-entropy function ½e² + (1/(γm))Σβ*log(β*/β) has V̇ = −k_ee² exactly and Σβ is conserved; (iv) at fixed positions with two
classes short of mass the replicator reaches the exact fit, transport alone does not, and with Σ|A_k|/a_max > 1 neither does; and at the
lift of an exact representation the weights add only neutral directions (within-class redistributions and one gauge per class), the rest
of the spectrum Proposition D‴.32's.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import cohesion_leak as CO
from src.core import mass_transfer as MT
from src.core.chaos_learning import AveragedSlowLoop

NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
FIVE = [(0.9, -1.6, 0.35), (-0.6, -0.8, 0.3), (0.7, 0.0, 0.35), (-0.5, 0.8, 0.3), (0.8, 1.6, 0.35)]
SIGMA0 = 0.15 * 0.047162567584665846
PER = 2 * np.pi
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])


def _loop(bumps, s_min=0.3, n_quad=128):
    fv = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in bumps)
    net = L.SaturatedLiftedNetwork(1, **dict(NET_KW, s_min=s_min))
    return net, AveragedSlowLoop(net, fv, XR, PER, G, SIGMA0, n_quad=n_quad), fv


def test_the_weighted_loop_is_the_gradient_flow_of_its_energy():
    net, base, fv = _loop(BUMPS)
    rng = np.random.default_rng(3)
    m = 8
    Th = net.unconstrained(rng.uniform(-2, 2, (m, 1)), rng.uniform(-2, 2, (m, 1)), rng.uniform(0.4, 0.9, m))
    beta = rng.uniform(0.5, 1.5, m); beta *= m / beta.sum()
    loop = MT.ReplicatorCohesionLoop(base, 1.5, 'gauss', gamma=0.3)
    dT, db = loop.drift(Th, beta)
    h = 1e-6
    gT = np.zeros_like(Th)
    for i in range(m):
        for j in range(3):
            E = np.zeros_like(Th); E[i, j] = h
            gT[i, j] = (loop.J(Th + E, beta) - loop.J(Th - E, beta)) / (2 * h)
    gb = np.array([(loop.J(Th, beta + h * np.eye(m)[i]) - loop.J(Th, beta - h * np.eye(m)[i])) / (2 * h) for i in range(m)])
    assert np.abs(dT + G[None, :] * gT).max() < 1e-7 * np.abs(dT).max()
    assert np.abs(db + 0.3 * beta * (gb - beta @ gb / beta.sum())).max() < 1e-7 * np.abs(db).max()
    assert abs(db.sum()) < 1e-14 * np.abs(db).max() * m
    d0, b0 = MT.ReplicatorCohesionLoop(base, 1.5, 'gauss', gamma=0.0).drift(Th, np.ones(m))
    assert np.abs(d0 - CO.CohesionLoop(base, 1.5, 'gauss').drift(Th)).max() < 1e-14 and np.abs(b0).max() == 0.0


def test_the_replicator_law_has_the_relative_entropy_lyapunov_function():
    net, base, fv = _loop(BUMPS)
    atoms = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    m, ke, gamma = 3, 25.0, 2.0
    f_vec = lambda X: fv(np.asarray(X)[:, 0])
    rhs = lambda t, z: MT.replicator_real_rhs(t, z, net, f_vec, XR, XRD, ke, atoms, gamma, m)
    sol = solve_ivp(rhs, (0.0, 60.0), np.concatenate([XR(0.0) + 0.05, [1.6, 0.7, 0.7]]), method='DOP853', rtol=1e-11, atol=1e-13,
                    t_eval=np.linspace(0, 60, 601))
    for i, t in enumerate(sol.t):
        z = sol.y[:, i]; zd = rhs(t, z); e = z[0] - XR(t)[0]
        dV = e * (zd[0] - XRD(t)[0]) - float(np.sum(zd[1:] / z[1:])) / (gamma * m)
        assert abs(dV + ke * e ** 2) < 1e-12 + 1e-10 * ke * e ** 2
    V = np.array([MT.kl_lyapunov(sol.y[:, i], XR, t, np.ones(3), gamma, m) for i, t in enumerate(sol.t)])
    assert np.all(np.diff(V) <= 1e-13) and V[-1] < 0.05 * V[0]
    assert np.abs(sol.y[1:].sum(axis=0) - m).max() < 1e-12


def test_at_fixed_positions_the_replicator_reaches_the_exact_fit_iff_the_masses_can_be_representable():
    W0 = np.array([0.15, 0.30, 0.20, 0.25, 0.10])                         # two classes short of mass: |A_k|/a_max = 0.225, 0.2

    class Frozen(MT.ReplicatorCohesionLoop):
        def drift(self, Th, beta):
            dT, db = super().drift(Th, beta)
            dT = dT.copy(); dT[:, 1:] = 0.0
            return dT, db
    fits = {}
    for scale in (1.0, 1.25):
        bumps = [(scale * A, C, S) for A, C, S in FIVE]
        net, base, fv = _loop(bumps, s_min=0.2)
        Th0 = net.unconstrained(np.full((5, 1), 0.5), np.array([[C] for A, C, S in bumps]), np.array([S for A, C, S in bumps]))
        for gamma in (0.0, 1e-2):
            fl = Frozen(base, 1.0, 'compact', sigma0=0.0, gamma=gamma)
            Th, beta = fl.run(Th0, 5 * W0, 20000.0, 0.25)
            fits[(scale, gamma)] = fl.fit_error(Th, beta)
            assert abs(beta.sum() - 5.0) < 1e-9                              # Σβ = m along the run (the balanced λ)
            assert abs(fl.drift(Th, beta)[1].sum()) < 1e-15 * 5
            if scale == 1.0 and gamma > 0:
                assert np.all(beta / 5 > np.array([abs(A) for A, C, S in bumps]) / 4)
    assert fits[(1.0, 1e-2)] < 1e-6 and fits[(1.0, 0.0)] > 0.3
    assert fits[(1.25, 1e-2)] > 0.0683 and fits[(1.25, 0.0)] > 0.3        # above the ℓ¹-ball infimum 0.0684 (round 24, ID-3)


def test_at_an_exact_representation_the_weights_add_only_neutral_directions():
    net, base, fv = _loop(BUMPS, n_quad=256)
    counts = np.array([3, 2, 2]); w = counts / 7
    atoms = net.unconstrained((np.array([A for A, C, S in BUMPS]) / w)[:, None], np.array([[C] for A, C, S in BUMPS]),
                              np.array([S for A, C, S in BUMPS]))
    Th0 = np.repeat(atoms, counts, axis=0)
    ev = {}
    for gamma in (0.0, 1e-2):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)

        def F(x):
            T = x[:21].reshape(7, 3); b = np.exp(x[21:])
            dT, db = loop.drift(T, b)
            return np.concatenate([dT.ravel(), db / b])
        x0 = np.concatenate([Th0.ravel(), np.zeros(7)])
        assert np.abs(F(x0)).max() < 1e-14
        Jm = np.array([(F(x0 + 1e-6 * e) - F(x0 - 1e-6 * e)) / 2e-6 for e in np.eye(28)]).T
        ev[gamma] = np.sort(np.linalg.eigvals(Jm).real)
    zero = {g: int(np.sum(np.abs(v) < 1e-9)) for g, v in ev.items()}
    assert zero[0.0] == 7 and zero[1e-2] == (7 - 3) + 3                   # within-class redistributions and one gauge per class
    sep, coll = CO.predicted_spectrum(base, atoms, counts, SIGMA0)
    nz = ev[0.0][np.abs(ev[0.0]) >= 1e-9]
    assert np.abs(nz - np.sort(np.concatenate([sep, coll]))).max() < 1e-9
    assert np.all(ev[1e-2][np.abs(ev[1e-2]) >= 1e-9] < 0)
    # round 24 (ID-6): the zero modes are the predicted vectors, and the separations keep −σ₀w_k with γ > 0
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)

    def F(x):
        T = x[:21].reshape(7, 3); b = np.exp(x[21:])
        dT, db = loop.drift(T, b)
        return np.concatenate([dT.ravel(), db / b])
    x0 = np.concatenate([Th0.ravel(), np.zeros(7)])
    Jv = lambda v: (F(x0 + 1e-6 * v) - F(x0 - 1e-6 * v)) / 2e-6
    lab = np.repeat(np.arange(3), counts)
    a_k = net.physical(atoms)[0][:, 0]
    da_k = 1.0 - (a_k / net.a_max) ** 2
    for k in range(3):
        idx = np.where(lab == k)[0]
        for j in idx[1:]:                                                   # redistributions inside class k
            v = np.zeros(28); v[21 + idx[0]] = 1.0; v[21 + j] = -1.0
            assert np.abs(Jv(v)).max() < 1e-9
        v = np.zeros(28); v[21 + idx] = 1.0                                 # the gauge of class k: δlog β against δα = −δβ a_k/a′_k
        for i in idx:
            v[3 * i] = -a_k[k] / da_k[k]
        assert np.abs(Jv(v)).max() < 1e-9
        if len(idx) >= 2:                                                   # a separation: −σ₀w_k exactly, no weight response
            v = np.zeros(28); v[3 * idx[0] + 1] = 1.0; v[3 * idx[1] + 1] = -1.0
            out = Jv(v)
            assert np.abs(out[:21] + SIGMA0 * w[k] * v[:21]).max() < 1e-9 and np.abs(out[21:]).max() < 1e-9


def test_the_fitness_against_an_independent_quadrature():
    """Round 24 (ID-6): the weights' drift γβ_i(u_i − λ) at an off-equilibrium configuration against u_i computed independently from the
    features a(α)e^{−(x−c)²/(2s²)} on the loop's quadrature (not through the module's features)."""
    net, base, fv = _loop(BUMPS)
    rng = np.random.default_rng(11)
    m = 6
    Th = net.unconstrained(rng.uniform(-2, 2, (m, 1)), rng.uniform(-2, 2, (m, 1)), rng.uniform(0.4, 0.9, m))
    beta = rng.uniform(0.5, 1.5, m); beta *= m / beta.sum()
    loop = MT.ReplicatorCohesionLoop(base, 1.5, 'gauss', gamma=0.7)
    a, c, s_ = net.physical(Th)
    X = base.X
    Phi = a[:, 0][None, :] * np.exp(-(X[:, None] - c[None, :, 0]) ** 2 / (2 * s_[None, :] ** 2))
    r = fv(X) - Phi @ beta / m
    u = (Phi.T @ r) / X.size
    lam = beta @ u / beta.sum()
    assert np.abs(loop.drift(Th, beta)[1] - 0.7 * beta * (u - lam)).max() < 1e-12
