"""Proposition D‴.24 (PARADIGM.md §5d): the finite-size floor is the leak's bias at the sampled weights, and an adaptive prior removes it
without a mass constraint. The nearest-atom leak toward an adaptive prior is the gradient flow of L̄ + (σ₀/2)(1/m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹}
in the metric with the cell masses, and that energy is the minimum of the OT energy over the prior's masses (semi-relaxed transport); at a
matching prior the function-level weight-sensitivity of Proposition D‴.9 is the weight-derivative of the leak's bias, Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰,
its floor law the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j, below the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ); and the exact representation of
sampled weights, with the prior at its atoms, is stationary for the adaptive prior while the fixed prior's stationary point is biased; at
the exact representation the splitting margins are σ₀ (R = 0) even where the fixed prior's point is a saddle; the fixed prior's weight
sensitivity is H_L(H_L + σ₀I)⁻¹ times the gauge's, the rest being the bias; and beyond the interior of the cells the energy stays
non-increasing across a cell crossing, while an empty cell keeps its atom frozen and stays empty."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB

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
    rep = lambda sh: net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + sh] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    exact, prior = rep(0.0), rep(shift)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    return loop, exact, prior


def test_the_nearest_atom_adaptive_flow_is_the_gradient_flow_of_the_free_mass_energy():
    loop, exact, prior = _setup(0.15)
    eta = 2 * SIGMA0
    apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest')
    rng = np.random.default_rng(3)
    counts = np.array([5, 4, 3]); m = int(counts.sum())
    Th = np.repeat(exact, counts, axis=0) + 0.05 * rng.standard_normal((m, 3)) * SQG[None, :]
    Z = prior + 0.02 * rng.standard_normal((3, 3))
    gamma, _ = apl.plan(Th, Z)
    w = gamma.sum(axis=0)
    assert np.allclose(w, counts / m)                                         # the cell masses are the sampled weights
    dTh, dZ = apl.drift(Th, Z, gamma)
    J = lambda T_, Z_: apl.energy(T_, Z_)
    h = 1e-6; gT = np.zeros_like(Th); gZ = np.zeros_like(Z)
    for i in range(m):
        for c in range(3):
            E = np.zeros_like(Th); E[i, c] = h
            gT[i, c] = (J(Th + E, Z) - J(Th - E, Z)) / (2 * h)
    for j in range(3):
        for c in range(3):
            E = np.zeros_like(Z); E[j, c] = h
            gZ[j, c] = (J(Th, Z + E) - J(Th, Z - E)) / (2 * h)
    assert np.abs(dTh + m * G[None, :] * gT).max() < 1e-7                      # dϑ_i/ds = −mG∇_{ϑ_i}𝒥_∘
    assert np.abs(dZ + (eta / (SIGMA0 * w[:, None])) * G[None, :] * gZ).max() < 1e-7   # dz_j/ds = −(η/(σ₀w_j))G∇_{z_j}𝒥_∘, w the cell masses


def test_the_free_mass_energy_is_the_minimum_of_the_OT_energy_over_the_masses():
    loop, exact, prior = _setup(0.0)
    rng = np.random.default_rng(5)
    m = 12
    for trial in range(4):
        Th = exact[rng.integers(0, 3, m)] + 1.2 * rng.standard_normal((m, 3)) * SQG[None, :]
        Z = exact + 0.3 * rng.standard_normal((3, 3))
        gN, dist = AP.nearest_plan(Th, Z, G)
        wN = gN.sum(axis=0)
        if np.any(wN == 0):
            continue
        assert abs(AP.plan_to(Th, Z, wN, G)[1] - dist) < 1e-12                 # attained at the cell masses …
        for _ in range(5):
            qq = rng.dirichlet(np.ones(3))
            assert AP.plan_to(Th, Z, qq, G)[1] >= dist - 1e-12               # … and a minimum over the masses
        for a in range(3):                                                    # … where it matters: the masses moved by ±1/m
            for b in range(3):
                if a != b and wN[a] > 0:
                    qq = wN.copy(); qq[a] -= 1.0 / m; qq[b] += 1.0 / m
                    assert AP.plan_to(Th, Z, qq, G)[1] > dist + 1e-6


def test_the_floor_law_is_the_weight_derivative_of_the_leaks_bias():
    loop, exact, prior = _setup(0.0)
    out = FB.bias_floor_law(loop, Q, exact, SIGMA0)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
    F2 = PD.function_sensitivity(loop.nearest, Q, exact, S) @ np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    assert np.linalg.norm(out['F_bias'] - F2) < 2e-6 * np.linalg.norm(F2)       # F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ (exact at q; the FD error)
    numeric = float(np.trace(F2 @ out['C2'] @ F2.T)) / F2.shape[0]
    assert abs(out['law'] - numeric) < 1e-6 * numeric                         # the mode terms sum to the law
    assert out['law'] < out['bound']                                           # each mode's weight ≤ σ₀/4
    # the weight is a band-pass centred at h = σ₀: against the unweighted shares z_jᵀC̃z_j (independent of σ₀) it moves the floor toward the
    # modes with h_j ≈ σ₀, at a small and at a large σ₀
    for s0 in (0.05 * 0.047162567584665846, 0.5 * 0.047162567584665846):
        o = FB.bias_floor_law(loop, Q, exact, s0)
        h = o['h']; wgt = h * s0 ** 2 / (h + s0) ** 2
        un = o['terms'] / wgt; un = un / un.sum(); tw = o['terms'] / o['terms'].sum()
        band = (h / s0 >= 3 - 2 * np.sqrt(2)) & (h / s0 <= 3 + 2 * np.sqrt(2))   # the weight at least half its peak σ₀/4
        assert tw[band].sum() > un[band].sum() + 0.05
        assert abs(np.sum(tw * np.log(h / s0))) < 0.8 * abs(np.sum(un * np.log(h / s0)))


def _margins(loop, w, atoms):
    """Proposition D‴.5's splitting margins σ₀ − λ_max R^G_k at the replicated configuration of the atoms with the weights w."""
    a, c, s_, d_, phi = loop.features(atoms)
    return SIGMA0 - loop.kappa_R(loop.F - phi @ (w * a), atoms)


def test_the_exact_representation_of_sampled_weights_is_stationary_for_the_adaptive_prior():
    loop, exact, prior = _setup(0.0)
    counts = np.array([5, 4, 6]); m = int(counts.sum()); w = counts / m            # a D‴.5-admissible sample of the fixed prior
    te, res, fit = FB.exact_weighted_atoms(loop, w, exact)
    assert res < 1e-12 and fit < 1e-12
    assert np.allclose(loop.features(te)[0], np.array([A for A, C, S in BUMPS]) / w)   # the amplitudes A_k/w_k (the gauge)
    Th = np.repeat(te, counts, axis=0)
    apl = AP.AdaptivePriorLoop(loop, SIGMA0, plan='nearest')
    gamma, _ = apl.plan(Th, te)
    dTh, dZ = apl.drift(Th, te, gamma)
    assert np.abs(dTh).max() < 1e-12 and np.abs(dZ).max() < 1e-12 and loop.fit_error(Th) < 1e-12
    Tw, resw = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)   # the fixed (matching) prior at the same weights:
    assert resw < 1e-12 and np.all(_margins(loop, w, Tw) > 0.003)                      # a strict local minimum (D‴.5-admissible) …
    assert loop.fit_error(np.repeat(Tw, counts, axis=0)) > 1e-3                         # … and biased: the floor
    counts = np.array([6, 5, 4]); w = counts / counts.sum()                            # a sample where the fixed prior's point is a saddle
    Tw, resw = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)
    assert resw < 1e-12 and _margins(loop, w, Tw)[0] < -0.003                         # class 1 can split under the fixed prior …
    te = FB.exact_weighted_atoms(loop, w, exact)[0]
    assert np.abs(_margins(loop, w, te) - SIGMA0).max() < 1e-12                      # … while at ϑ^e(w) R = 0: every margin is σ₀


def test_the_weight_fluctuation_splits_between_the_parameters_and_the_bias():
    # (iii): the fixed prior's weight sensitivity, whitened and weighted, is H_L(H_L + σ₀I)⁻¹ times the gauge's S^e = ∂_wϑ^e; the adaptive
    # prior carries all of S^e, so its parameter-level (coupled-distance) floor is the larger one
    loop, exact, prior = _setup(0.0)
    out = FB.bias_floor_law(loop, Q, exact, SIGMA0)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
    V2 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    Sy = np.column_stack([(np.sqrt(Q)[:, None] * np.einsum('kcj,j->kc', S, V2[:, c]) / SQG[None, :]).ravel() for c in range(2)])
    Sey = -out['dx0']
    split = out['HL'] @ np.linalg.solve(out['HL'] + SIGMA0 * np.eye(9), Sey)
    assert np.linalg.norm(Sy - split) < 1e-5 * np.linalg.norm(Sy)
    assert np.linalg.norm(Sy - Sey) > 0.1 * np.linalg.norm(Sy)                 # the bias part is not negligible
    D_fixed, D_gauge = np.trace(Sy @ out['C2'] @ Sy.T), np.trace(Sey @ out['C2'] @ Sey.T)
    assert D_gauge > 1.5 * D_fixed


def test_beyond_the_cells_the_energy_decreases_across_a_crossing_and_an_empty_cell_stays_empty():
    """(ii) beyond the interior of the cells: from a spread cloud a particle crosses into another cell, and 𝒥_∘ — the minimum over the cells,
    whose active branch's derivative bounds its Dini derivative — stays non-increasing across the change of counts (the vertex-selection
    flow); a cell without particles has a degenerate metric, its atom is frozen, and the cell stays empty."""
    loop, exact, prior = _setup(0.0)
    m = 12
    rng = np.random.default_rng(4)
    Th = exact[rng.integers(0, 3, m)] + 1.0 * rng.standard_normal((m, 3)) * SQG[None, :]
    for eta in (SIGMA0, np.inf):
        apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest'); rec = []
        apl.run(Th.copy(), exact.copy(), 200.0, 0.25, plan_every=4, rec_every=4,
                callback=lambda s_, T_, Z_, ad: rec.append((np.round(apl.plan(T_, Z_)[0].sum(axis=0) * m), apl.energy(T_, Z_))))
        counts = np.array([r[0] for r in rec]); J = np.array([r[1] for r in rec])
        assert np.any(np.any(np.diff(counts, axis=0) != 0, axis=1))                # a crossing happens …
        assert np.all(np.diff(J) < 0)                                             # … and 𝒥_∘ decreases at every refresh
    Th = np.repeat(exact[:2], 6, axis=0) + 0.02 * np.random.default_rng(1).standard_normal((m, 3)) * SQG[None, :]   # no particle near atom 3
    for eta in (SIGMA0, np.inf):
        apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest'); rec = []
        apl.run(Th.copy(), exact.copy(), 300.0, 0.25, plan_every=4, rec_every=4,
                callback=lambda s_, T_, Z_, ad: rec.append((apl.plan(T_, Z_)[0].sum(axis=0) * m, Z_[2].copy(), apl.energy(T_, Z_))))
        assert all(np.allclose(c, [6, 6, 0]) for c, _, _ in rec)                  # the third cell stays empty …
        assert all(np.array_equal(z, exact[2]) for _, z, _ in rec)                 # … its atom frozen
        assert np.all(np.diff([J for _, _, J in rec]) < 0)
