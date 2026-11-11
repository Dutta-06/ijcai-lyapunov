"""Proposition D‴.9 (PARADIGM.md §5d): the averaged slow loop with a measure-valued prior (nearest-atom leak): consistency with
the single-prior loop, the gradient of J_m^ν, the coupling identity inside a cell, the exact representation as a stationary
measure with zero residual and maximal margin under the matching prior, and the weight-sensitivity laws at parameter and
function level (both floors are higher with the matching prior: the price of removing the bias)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import prior_design as PD

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846


def _net():
    return L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)


def _exact(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def test_single_atom_prior_is_the_sigma_modification_and_the_gradient_of_J_is_exact():
    net = _net()
    rng = np.random.default_rng(5)
    Th = CL.cold_start(net, 10, rng)
    single = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, SIGMA0, n_quad=200)
    mix1 = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, np.zeros((1, 3)), n_quad=200)
    assert np.abs(single.drift(Th) - mix1.drift(Th)).max() < 1e-14 and np.abs(single.grad_J(Th) - mix1.grad_J(Th)).max() < 1e-14
    prior = _exact(net) + np.array([[0.0, 0.3, 0.0]])
    mix = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=200)
    g = mix.grad_J(Th)
    h = 1e-6
    for (i, j) in ((0, 0), (3, 1), (7, 2), (9, 0)):
        E = np.zeros_like(Th); E[i, j] = h
        fd = (mix.J(Th + E) - mix.J(Th - E)) / (2 * h)
        assert abs(fd - g[i, j]) < 1e-6 * max(1.0, abs(g[i, j]))          # away from the Voronoi boundaries
    assert np.abs(mix.drift(Th) + G[None, :] * g).max() < 1e-13


def test_coupling_identity_inside_the_cells():
    net = _net()
    rng = np.random.default_rng(7)
    prior = _exact(net)
    mix = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=200)
    Thb = CL.mixture_start(prior, 12, rng, 0.05, G)
    Th = Thb + 0.03 * rng.standard_normal(Thb.shape)
    d = mix.decomposition(Th, Thb, fhat_mf=mix.f_hat(CL.mixture_start(prior, 300, rng, 0.05, G)))
    assert d['same_cell'] and abs(d['total'] - d['half_dD']) < 1e-15 * max(1.0, abs(d['half_dD']))
    # the content of (3.1): T_B = −‖h‖² + a cross term bounded by (L₂^G/2)‖h‖D (L₂^G = 11.11 on this network), so T_B < 0 here
    assert abs(d['T_B'] + d['h_norm'] ** 2) <= 11.12 / 2 * d['h_norm'] * d['D'] and d['T_B'] < 0


def test_matching_prior_gives_the_exact_representation_with_maximal_margin_and_the_two_sensitivity_laws():
    net = _net()
    exact = _exact(net)
    w = np.ones(3) / 3
    # the σ-modification's stationary three-atom measure, for comparison
    single = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, SIGMA0, n_quad=256)
    atoms_s, g_s, eig_s = CL.stationary_atoms(single, exact)
    k_s = single.kappa_R_neighbourhood(single.residual(atoms_s), atoms_s, 0.1, n_dirs=100, radii=(1.0, 0.5))
    # the matching prior: the exact representation is stationary with zero residual and zero curvature
    mix = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=256)
    atoms, res = PD.weighted_stationary_atoms_prior(mix, w, exact)
    assert res < 1e-12 and np.abs(atoms - exact).max() < 1e-9
    r = mix.residual(atoms)
    assert np.abs(r).max() < 1e-12 and np.abs(mix.kappa_R(r, atoms)).max() < 1e-10
    k_m = mix.kappa_R_neighbourhood(r, atoms, 0.1, n_dirs=100, radii=(1.0, 0.5))
    assert abs(k_m) < 1e-10 and SIGMA0 - k_m > SIGMA0 - k_s > 0             # margin σ₀ (maximal) against the σ-modification's
    assert mix.cell_margin(atoms).min() > 2 * 0.1                           # the assignment is constant on N_{0.1}(supp μ*)
    # the weight-sensitivity law: parameter level (a_k ∝ 1/w_k at an exact fit) and function level (the weights are invisible)
    S = PD.weight_sensitivity_prior(mix, w, atoms)
    S_s = CL.weight_sensitivity(single, w, atoms_s)
    Fm = PD.function_sensitivity(mix, w, atoms, S)
    Fs = PD.function_sensitivity(single, w, atoms_s, S_s)
    Dm, Ds = CL.predicted_stationary_D(S, w, 64, 2048, G), CL.predicted_stationary_D(S_s, w, 64, 2048, G)
    fm, fs = PD.predicted_stationary_function_D(Fm, w, 64, 2048), PD.predicted_stationary_function_D(Fs, w, 64, 2048)
    assert Dm > Ds and fm > fs                                               # both floors are higher with the matching prior: the exact fit has
    assert 1.4 < Dm / Ds < 1.9 and 1.9 < fm / fs < 2.5                       # larger amplitudes (table 1: 1.61× and 2.18×; bias against variance)
    # the function-level sensitivity against finite differences of the fitted function of the weighted stationary point
    h = 1e-4
    for j in range(3):
        wp = w.copy(); wp[j] += h
        wm = w.copy(); wm[j] -= h
        Tp, _ = PD.weighted_stationary_atoms_prior(mix, wp, atoms)
        Tm, _ = PD.weighted_stationary_atoms_prior(mix, wm, atoms)
        ap, cp, sp, dp, php = mix.features(Tp); am, cm_, sm, dm, phm = mix.features(Tm)
        fd = ((php @ (wp * ap)) - (phm @ (wm * am))) / (2 * h)
        assert np.abs(fd - Fm[:, j]).max() < 1e-5 * max(1.0, np.abs(Fm[:, j]).max())
