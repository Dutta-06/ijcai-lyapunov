"""Theorem D‴.21 (PARADIGM.md §5d): a measure-level basin for the barycentric OT leak at commensurate sizes.
  1. The energy identity (E) and the gradient identity (P) along the plan segments hold exactly for every assignment with the counts mq
     (not only the optimal one), at the matching prior (r* = 0) and at a shifted prior (r* ≠ 0); with the optimal plan E is J^{OT}/m − 𝒥*;
     with other counts N both acquire exactly the defects Σ_k(N_k/m − q_k)c_k and −⟨Δf̂, f̂[μ_N] − f̂[μ*]⟩ (the role of the counts).
  2. The Polyak–Łojasiewicz chain |∇J_a|²_G/(mE) ≥ Q²/(DE) ≥ λ holds wherever the plan-segment margin α is positive, and at the matching prior
     the constant 2σ₀ is attained by the separations (it is sharp).
  3. From the lift of q with one particle moved into another cell the margin α is positive, the plan keeps the particle in its class, and the
     excess energy decays below the a-posteriori bound E(0)exp(−∫λ₊), at the rate 2σ₀ at the end; under the nearest-atom leak from the same
     start the energy stays at the floor.
  4. The nearest-atom energy has a critical point at the replicated configuration of an admissible count vector other than mq (a spurious
     minimum, E > 0: the separation margins and the whitened weighted K-atom Hessian positive); under the OT leak that configuration is not
     critical and the PL inequality holds there."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import ot_basin as OB
from src.core import prior_design as PD

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup(shift=0.0, n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    prior = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    atoms, res, (eps, half_gap) = loop.stationary_atoms(prior)
    assert res < 1e-12 and eps < half_gap                                     # μ* stationary, its diagonal plan optimal
    return loop, atoms, OB.BasinDiagnostics(loop, atoms)


def test_the_energy_and_gradient_identities_hold_exactly_for_every_assignment_with_the_counts():
    for shift in (0.0, 0.15):
        loop, atoms, B = _setup(shift)
        assert B.stationarity < 1e-12
        rng = np.random.default_rng(4)
        for m in (12, 30):
            for trial in range(3):
                lab = rng.permutation(np.repeat(np.arange(3), m // 3))         # any assignment with the counts mq, optimal or not
                Th = atoms[lab] + 0.4 * rng.standard_normal((m, 3)) * SQG[None, :]
                t = B.terms(Th, lab)
                assert abs(t['E'] - t['E_ident']) < 1e-12 * abs(t['E'])
                assert abs(t['Q'] - t['Q_ident']) < 1e-12 * abs(t['Q'])
            t = B.terms(Th)                                                   # the optimal plan: E is the OT functional's excess
            assert abs(t['E'] - (loop.J(Th) / m - B.Jstar)) < 1e-12 * max(1.0, abs(t['E']))
        if shift > 0:                                                         # r* ≠ 0: the counts matter
            for m, counts in ((12, (5, 4, 3)), (30, (8, 12, 10))):
                lab = rng.permutation(np.repeat(np.arange(3), counts))
                Th = atoms[lab] + 0.3 * rng.standard_normal((m, 3)) * SQG[None, :]
                t = B.terms(Th, lab, check_counts=False)
                dE, dQ = B.counts_mismatch(t['counts'], t['df'])
                assert abs(dE) > 1e-6 and abs(dQ) > 1e-6                     # the defects are not negligible …
                assert abs(t['E'] - t['E_ident'] - dE) < 1e-12 * max(abs(t['E']), abs(dE))   # … and they are exactly these
                assert abs(t['Q'] - t['Q_ident'] - dQ) < 1e-12 * max(abs(t['Q']), abs(dQ))


def test_the_PL_chain_holds_and_its_constant_two_sigma0_is_attained_by_the_separations():
    loop, atoms, B = _setup(0.0)
    m = 30
    lab0 = np.repeat(np.arange(3), m // 3)
    rng = np.random.default_rng(2)
    u = rng.standard_normal((m, 3)) * SQG[None, :]
    for k in range(3):
        u[lab0 == k] -= u[lab0 == k].mean(axis=0)                             # a pure separation: zero class means
    u *= 1e-3 / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1)).max()
    Th = atoms[lab0] + u
    t = B.terms(Th)
    alpha, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab']))
    assert abs(t['ratio'] / (2 * SIGMA0) - 1) < 1e-4 and abs(lam / (2 * SIGMA0) - 1) < 1e-4   # observed = theorem's = 2σ₀: sharp
    for shift in (0.0, 0.15):
        loop, atoms, B = _setup(shift)
        for trial in range(4):
            Th = atoms[lab0] + 0.05 * rng.standard_normal((m, 3)) * SQG[None, :]
            t = B.terms(Th)
            alpha, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab']))
            assert alpha > 0 and t['E'] > 0
            assert t['ratio'] >= t['lam_eff'] * (1 - 1e-12) and t['lam_eff'] >= lam * (1 - 1e-9) > 0


def test_the_excess_energy_decays_below_the_a_posteriori_bound_from_a_particle_in_another_cell():
    loop, atoms, B = _setup(0.0)
    m = 30
    lab0 = np.repeat(np.arange(3), m // 3)
    Th0 = atoms[lab0].copy()
    Th0[0] = 0.4 * atoms[0] + 0.6 * atoms[1]                                  # one particle of class 0 moved into cell 1
    assert np.argmin(loop.cost_matrix(Th0[:1])[0]) == 1
    out = {}
    for leak in ('ot', 'nearest'):
        Th = Th0.copy()
        t = B.terms(Th)
        if leak == 'ot':
            assert t['lab'][0] == 0                                           # the plan keeps the particle in its class
        alpha, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=5))
        assert alpha > 0                                                      # the margin is positive at the start
        E0, Es, bound = t['E'], [t['E']], 0.0
        ds, lam_prev = 0.5, lam
        for k in range(600):                                                  # 300 slow units, midpoint rule
            T = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
            k1 = loop.drift(Th, T)
            Th = Th + ds * loop.drift(Th + 0.5 * ds * k1, T)
            if k % 10 == 9:
                t = B.terms(Th)
                alpha, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=5))
                bound += 10 * ds * min(lam_prev, lam)                         # λ held at the smaller endpoint value over the interval
                lam_prev = lam
                Es.append(t['E'])
        out[leak] = (E0, np.array(Es), bound, t, Th)
    E0, Es, bound, t, Th = out['ot']
    assert np.all(np.diff(Es) < 0)                                            # a gradient flow: E decreases
    assert Es[-1] <= 1.01 * E0 * np.exp(-bound)                               # E(s) ≤ E(0)exp(−∫λ₊)
    assert abs(t['ratio'] / (2 * SIGMA0) - 1) < 0.01                          # at the end the rate is 2σ₀
    assert t['D'] < 1e-2 and np.argmin(loop.cost_matrix(Th[:1])[0]) == 0      # back in its own cell, near its atom
    E0n, Esn, _, tn, Thn = out['nearest']
    assert Esn[-1] > 1e-3 and Esn[-1] > 100 * Es[-1]                          # the nearest-atom leak keeps the particle in cell 1: the floor


def test_the_nearest_atom_energy_has_a_spurious_minimum_that_the_mass_constraint_removes():
    loop, atoms, B = _setup(0.0)
    m, counts = 30, np.array([11, 9, 10])
    w = counts / m
    Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, atoms, w0=Q, n_cont=4)
    a, c, s, diff, phi = loop.features(Tw)
    margin = SIGMA0 - loop.kappa_R(loop.F - phi @ (w * a), Tw).max()
    lift = np.repeat(Tw, counts, axis=0)
    d = 9; h = 1e-6; Jac = np.zeros((d, d))                                   # the whitened weighted K-atom Hessian of the nearest-atom energy
    gfun = lambda v: PD.weighted_energy_grad_prior(loop.nearest, w, v.reshape(3, 3)).ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jac[:, j] = (gfun(Tw.ravel() + e) - gfun(Tw.ravel() - e)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3) * np.tile(SQG, 3); sv = np.tile(SQG, 3) / np.repeat(np.sqrt(w), 3)
    Hw = sw[:, None] * Jac * sv[None, :]
    assert res < 1e-12 and margin > 0                                         # admissible weights: the separations …
    assert np.linalg.eigvalsh(0.5 * (Hw + Hw.T)).min() > 0                   # … and the collective block positive: a strict minimum
    assert np.abs(loop.nearest.drift(lift)).max() < 1e-12                    # critical for the nearest-atom leak
    t = B.terms(lift)
    alpha, lam = B.theorem_rate(**B.segment_curvatures(lift, t['lab']))
    assert t['E'] > 1e-3                                                      # at a positive excess energy
    assert np.abs(loop.drift(lift)).max() > 1e-2                              # not critical under the OT leak
    assert alpha > 0 and t['ratio'] >= lam > 0                                # and the PL inequality holds there
