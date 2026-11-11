"""Theorem D‴.18 (PARADIGM.md §5d): the constants of Theorem J are uniform in m in the per-particle norm |||·||| at a replicated configuration.
Checked on the averaged slow loop at the prior whose centres are shifted by +0.15, at the matching masses q: there the stationary residual is
not zero (R^G_k ≠ 0, λ_sep = 0.80σ₀ < σ₀), so the per-particle terms that vanish at the matching prior are exercised.
  1. The second-order remainder of the drift is bounded by C_N|||u|||² with one C_N for m = 6, 30, 150, 600 (not growing with m); for a
     one-particle displacement the displaced particle's remainder is exactly N_a + N_b/m with N_a m-independent — non-zero at the shifted
     prior, zero at the matching prior (where the residual vanishes) — and every other particle's remainder is exactly O(1/m).
  2. A separation perturbation (two particles of one class moved by ±v) moves only those two particles at first order, with the response
     (GR_k − σ₀I)v — not the leak's −σ₀v — at every m.
  3. The averaged loop contracts |||u||| at the rate min{λ_coll, λ_sep} itself, with one constant C ≤ 1.05 for m = 6, 30, 150."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import prior_design as PD

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup(shift=0.15, n_quad=96):
    """The averaged loop toward the prior with centres shifted by `shift` and its weighted stationary point Θ*(q) (Newton)."""
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    prior = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    loop = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    atoms, res = PD.weighted_stationary_atoms_prior(loop, Q, prior)
    assert res < 1e-12
    return loop, atoms


def _norm3(u, lab, w):
    """|||u||| = max{(Σ_k w_k‖ū_k‖²_{G⁻¹})^{1/2}, max_i‖u_i − ū_{lab(i)}‖_{G⁻¹}}."""
    K = len(w)
    ubar = np.array([u[lab == k].mean(axis=0) for k in range(K)])
    coll = np.sqrt(np.sum(w * np.sum(ubar ** 2 / G[None, :], axis=1)))
    sep = np.sqrt(np.sum((u - ubar[lab]) ** 2 / G[None, :], axis=1)).max()
    return max(coll, sep)


def _gnorm(x):
    return float(np.sqrt(np.sum(np.atleast_2d(x) ** 2 / G[None, :], axis=1)).max())


def _margins(loop, atoms):
    """λ_sep = min_k λ_min(σ₀I − R^G_k) and λ_coll = λ_min of the whitened, weighted collective block (Proposition D‴.5(i)); and R."""
    a, c, s, diff, phi = loop.features(atoms)
    r = loop.F - phi @ (Q * a)
    R = loop.residual_curvature(r, atoms)
    Rw = SQG[None, :, None] * R * SQG[None, None, :]
    lam_sep = min(float(np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k]).min()) for k in range(3))
    d, h = 9, 1e-6
    v = atoms.ravel(); Jm = np.zeros((d, d))
    gf = lambda V: PD.weighted_energy_grad_prior(loop, Q, V.reshape(3, 3)).ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (gf(v + e) - gf(v - e)) / (2 * h)
    sw = np.repeat(np.sqrt(Q), 3); sg = np.tile(SQG, 3)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]
    lam_coll = float(np.linalg.eigvalsh(0.5 * (Hs + Hs.T)).min())
    return lam_sep, lam_coll, R


def test_second_order_remainder_of_the_drift_is_uniform_in_m_in_the_per_particle_norm():
    h = 1e-5
    v1 = 0.2 * SQG * np.array([0.6, -0.5, 0.6]) / _gnorm(SQG * np.array([0.6, -0.5, 0.6]))   # a one-particle displacement of size 0.2
    Na = {}
    for shift in (0.15, 0.0):
        loop, atoms = _setup(shift)
        rem = lambda Th0, d0, u: loop.drift(Th0 + u) - d0 - (loop.drift(Th0 + h * u) - loop.drift(Th0 - h * u)) / (2 * h)
        rng = np.random.default_rng(0)
        CN, N0, mNother = {}, {}, {}
        for m in (6, 30, 150, 600):
            lab = np.repeat(np.arange(3), m // 3)
            Th0 = atoms[lab]
            d0 = loop.drift(Th0)
            assert np.abs(d0).max() < 1e-12                                   # the replicated configuration is stationary
            if shift:
                ratios = []
                for trial in range(4):
                    u = rng.standard_normal((m, 3)) * SQG[None, :]
                    u *= 0.2 / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1)).max()  # every particle displaced by ≤ 0.2 (G⁻¹)
                    N = rem(Th0, d0, u)
                    ratios.append(_gnorm(N) / _norm3(u, lab, Q) ** 2)
                CN[m] = max(ratios)
            u = np.zeros((m, 3)); u[0] = v1
            N = rem(Th0, d0, u)
            N0[m] = N[0]; mNother[m] = m * _gnorm(N[1:])
        if shift:
            vals = np.array(list(CN.values()))
            assert vals.max() < 0.1 and CN[600] <= CN[6]                      # ‖N(u)‖ ≤ C_N|||u|||² with one C_N for every m (not growing)
        # the displaced particle's remainder is N_a + N_b/m exactly (the drift is affine in the empirical measure, whose part from the
        # displaced particle is 1/m): N_a from every pair of sizes is the same vector
        ms = sorted(N0)
        pairs = [(ms[i] * N0[ms[i]] - ms[j] * N0[ms[j]]) / (ms[i] - ms[j]) for i in range(4) for j in range(i + 1, 4)]
        assert max(_gnorm(p - pairs[0]) for p in pairs) < 1e-10
        Na[shift] = _gnorm(pairs[0])
        mo = np.array(list(mNother.values()))
        assert mo.max() - mo.min() < 1e-6 * mo.max()                          # every other particle's remainder is exactly O(1/m)
    # the per-particle coefficient ⟨r∇³Φ̃⟩ is non-zero at a non-zero residual (1.65·10⁻⁵ here) and vanishes at the matching prior
    assert Na[0.15] > 1e-5 and Na[0.0] < 1e-11


def test_a_separation_perturbation_moves_only_its_own_particles_at_first_order_and_is_the_same_at_every_m():
    loop, atoms = _setup(0.15)
    lam_sep, lam_coll, R = _margins(loop, atoms)
    v = np.array([0.1, 0.05, -0.05])
    expected = G * (R[0] @ v) - SIGMA0 * v                                   # (GR₀ − σ₀I)v: the separation block of class 0
    assert np.abs(expected + SIGMA0 * v).max() > 5e-5                        # R₀ ≠ 0: the response is not the leak's alone
    resp = {}
    for m in (6, 60, 300):
        lab = np.repeat(np.arange(3), m // 3)
        Th0 = atoms[lab]
        u = np.zeros((m, 3)); u[0] = v; u[1] = -v                            # zero mean within class 0: a separation
        h = 1e-5
        lin = (loop.drift(Th0 + h * u) - loop.drift(Th0 - h * u)) / (2 * h)
        assert np.abs(lin[2:]).max() < 1e-13                                 # no other particle responds at first order
        assert np.abs(lin[0] + lin[1]).max() < 1e-13                         # the class mean does not move
        resp[m] = lin[0]
    for m in resp:
        assert np.abs(resp[m] - expected).max() < 1e-9 * np.abs(expected).max()


def test_the_averaged_loop_contracts_the_per_particle_norm_at_the_predicted_rate_with_one_constant():
    loop, atoms = _setup(0.15, n_quad=64)
    lam_sep, lam_coll, R = _margins(loop, atoms)
    mu = min(lam_sep, lam_coll)
    assert 0.75 * SIGMA0 < mu < 0.85 * SIGMA0                                # λ_sep = 0.80σ₀ < λ_coll at the shifted prior
    rng = np.random.default_rng(3)
    consts = {}
    for m in (6, 30, 150):
        lab = np.repeat(np.arange(3), m // 3)
        u = rng.standard_normal((m, 3)) * SQG[None, :]
        u *= 0.5 / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1)).max()
        Th = atoms[lab] + u
        n0 = _norm3(u, lab, Q)
        ds, worst = 1.0, 0.0
        for k in range(400):                                                 # midpoint rule in slow time, 400 slow units
            k1 = loop.drift(Th); Th = Th + ds * loop.drift(Th + 0.5 * ds * k1)
            worst = max(worst, _norm3(Th - atoms[lab], lab, Q) * np.exp(mu * (k + 1) * ds) / n0)
        consts[m] = worst
    vals = np.array(list(consts.values()))
    assert vals.max() <= 1.05 and vals.max() < 1.1 * vals.min()              # |||u(s)||| ≤ C e^{−μs}|||u(0)||| at μ = min{λ_coll, λ_sep}, one C
