"""The averaged slow loop on ℝⁿ (src/core/planar_loop.py) and Proposition D‴.24 in two dimensions. For n = 1 the loop is the scalar one of
chaos_learning; for n = 2 its learning term is minus m times the gradient of L̄ (finite differences) and the average over the orbit of the
real loop's gradient law at the error r(x); the exact representation of representable sampled weights, lifted with the adaptive prior at its
atoms, is stationary for the nearest-atom adaptive dynamics while the fixed prior's weighted stationary point there is biased; and the
function-level weight-sensitivity of the fixed prior is the weight-derivative of its bias, F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰, with the mode terms
summing to the law and the law below the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import planar_loop as PL

PER = 2 * np.pi
G2 = np.array([1.0, 1.0, 0.25, 0.25, 0.25])
B2 = [(np.array([0.8, 0.3]), np.array([-1.2, 0.8]), 0.7), (np.array([-0.5, 0.6]), np.array([1.0, 1.0]), 0.6),
      (np.array([0.4, -0.7]), np.array([0.8, -1.1]), 0.6), (np.array([-0.6, -0.4]), np.array([-1.0, -0.9]), 0.7)]
A2 = np.array([b[0] for b in B2]); C2 = np.array([b[1] for b in B2]); S2 = np.array([b[2] for b in B2])
SIGMA0 = 0.05
Q4 = np.ones(4) / 4


def f2(X):
    X = np.atleast_2d(X)
    out = np.zeros_like(X)
    for A, C, S in B2:
        out += A[None, :] * np.exp(-np.sum((X - C[None, :]) ** 2, axis=1) / (2 * S ** 2))[:, None]
    return out


def _planar(n_quad=200):
    net = L.SaturatedLiftedNetwork(2, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    loop = PL.AveragedLoopN(net, f2, lambda t: np.array([2.0 * np.sin(t), 1.5 * np.sin(2 * t)]), PER, G2, SIGMA0, n_quad=n_quad)
    exact = lambda w: net.unconstrained(A2 / np.asarray(w, float)[:, None], C2, S2)
    return net, loop, exact


def test_the_planar_loop_reduces_to_the_scalar_loop_for_n_1():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    B = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
    f1 = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in B)
    G1 = np.array([1.0, 0.25, 0.25])
    old = CL.AveragedSlowLoop(net, f1, lambda t: np.array([2 * np.sin(t)]), PER, G1, 0.007, theta0=np.zeros(3), n_quad=128)
    new = PL.AveragedLoopN(net, lambda X: f1(X[:, 0])[:, None], lambda t: np.array([2 * np.sin(t)]), PER, G1, 0.007, n_quad=128)
    Th = np.random.default_rng(0).standard_normal((7, 3))
    assert np.abs(old.drift(Th) - new.drift(Th, np.zeros_like(Th))).max() < 1e-14


def test_the_learning_term_is_the_loss_gradient_and_the_orbit_average_of_the_real_law():
    net, loop, exact = _planar()
    rng = np.random.default_rng(1)
    Th = exact(Q4)[rng.integers(0, 4, 9)] + 0.3 * rng.standard_normal((9, 5))
    m, h = Th.shape[0], 1e-6
    num = np.zeros_like(Th)
    for i in range(m):
        for c in range(5):
            e = np.zeros_like(Th); e[i, c] = h
            num[i, c] = -m * (loop.L_bar(Th + e) - loop.L_bar(Th - e)) / (2 * h)
    lrn = loop.learning(Th)
    assert np.abs(lrn - num).max() < 1e-7 * np.abs(num).max()
    r = loop.residual(Th)
    avg = np.mean([net.grad_g(loop.X[k], r[k], Th) for k in range(loop.N)], axis=0)
    assert np.abs(lrn - avg).max() < 1e-14


def test_the_exact_representation_of_representable_sampled_weights_is_stationary_for_the_adaptive_prior():
    net, loop, exact = _planar()
    counts = np.array([6, 4, 5, 5]); w = counts / counts.sum()                   # representable: w > (0.2, 0.15, 0.175, 0.15)
    assert np.all(w > np.abs(A2).max(axis=1) / net.a_max)
    e = exact(w)
    Th = e[np.repeat(np.arange(4), counts)]
    assert loop.fit_error(Th) < 1e-12                                             # the lift fits f exactly
    lab = PL.cells(Th, e, G2)
    assert np.abs(PL.cell_means(Th, lab, e) - e).max() < 1e-14                     # the adaptive prior's atoms are the cell means …
    assert np.abs(loop.drift(Th, e[lab])).max() < 1e-14                           # … and the particles are stationary
    Th2, Z2 = PL.run(loop, Th, e, 50.0, 0.25, eta=np.inf)
    assert loop.fit_error(Th2) < 1e-12 and np.abs(Z2 - e).max() < 1e-12           # an invariant set of the dynamics
    at, res, own = PL.weighted_stationary(loop, w, exact(Q4))                    # the fixed prior's stationary point at w: biased
    assert res < 1e-12 and own and loop.fit_error(at, w) > 1e-2


def test_the_floor_is_the_weight_derivative_of_the_leaks_bias_in_two_dimensions():
    net, loop, exact = _planar()
    e0 = exact(Q4)
    H, hs, V = PL.gauss_newton(loop, Q4, e0)
    J = PL.whitened_jacobian(loop, Q4, e0)
    d = 1e-4
    F = np.zeros((loop.N * 2, 3)); dx0 = np.zeros((20, 3))
    for j in range(3):
        fh, x0 = [], []
        for sgn in (1, -1):
            w = Q4.copy(); w[j] += sgn * d; w[3] -= sgn * d
            at, res, own = PL.weighted_stationary(loop, w, e0)
            assert res < 1e-12 and own
            fh.append(loop.f_hat(at, w).ravel()); x0.append(PL.whitened(e0 - exact(w), w, G2))
        F[:, j] = (fh[0] - fh[1]) / (2 * d); dx0[:, j] = (x0[0] - x0[1]) / (2 * d)
    Fb = J @ (SIGMA0 * np.linalg.solve(H + SIGMA0 * np.eye(20), dx0))
    assert np.linalg.norm(F - Fb) < 1e-4 * np.linalg.norm(F)                     # the sensitivity is the bias's derivative
    Ct = np.diag(Q4[:3]) - np.outer(Q4[:3], Q4[:3])
    law = np.trace(F @ Ct @ F.T) / loop.N
    z = V.T @ dx0
    terms = hs * SIGMA0 ** 2 / (hs + SIGMA0) ** 2 * np.einsum('jk,kl,jl->j', z, Ct, z)
    assert abs(np.sum(terms) / law - 1) < 1e-3                                    # the mode terms sum to the law …
    assert law < 0.25 * SIGMA0 * np.trace(dx0 @ Ct @ dx0.T)                       # … which is below the bound


def test_the_finite_eta_branch_of_run_is_a_second_order_scheme_of_the_frozen_cell_flow():
    """run() with a finite η (the prior in the state) against DOP853 on the same flow with the cells frozen: errors falling by 4 per halving
    of dt."""
    from scipy.integrate import solve_ivp
    net, loop, exact = _planar(n_quad=96)
    counts = np.array([6, 4, 5, 5]); w = counts / counts.sum()
    rng = np.random.default_rng(3)
    Th0 = exact(w)[np.repeat(np.arange(4), counts)] + 0.05 * rng.standard_normal((20, 5))
    Z0 = exact(Q4) + 0.05 * rng.standard_normal((4, 5))
    eta, S = 10 * SIGMA0, 20.0
    lab = PL.cells(Th0, Z0, G2)

    def f(s, v):
        Th, Z = v[:100].reshape(20, 5), v[100:].reshape(4, 5)
        return np.concatenate([loop.drift(Th, Z[lab]).ravel(), (-eta * (Z - PL.cell_means(Th, lab, Z))).ravel()])
    ref = solve_ivp(f, (0.0, S), np.concatenate([Th0.ravel(), Z0.ravel()]), method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1]
    errs = []
    for dt in (0.5, 0.25):
        Th, Z = PL.run(loop, Th0, Z0, S, dt, eta=eta, cell_every=10 ** 9)
        errs.append(np.abs(np.concatenate([Th.ravel(), Z.ravel()]) - ref).max())
    assert errs[1] < 1e-4 and 3.5 < errs[0] / errs[1] < 4.5                            # second order
