"""Proposition D‴.12 (PARADIGM.md §5d): the barycentric OT leak toward a prior measure — the plan is a vertex with the prescribed
masses and the single-atom prior is the σ-modification; the leak is the G-gradient flow of J_m^{OT} (Danskin, commensurate and
incommensurate m); the replicated weighted configuration Θ*_m(q) is stationary for every commensurate m with the diagonal plan,
so the stationary measure, the learned function and the W₂/function-level floors are m-independent (also for mass-mismatched
masses q, as long as ε < Δ_ν/2); conservation of the cluster masses under the nearest-atom leak against their constraint under
the OT leak on a short run; the forest bound on the split particles; the exact W² between empirical measures."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import prior_design as PD
from src.core import ot_leak as OT

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


def _fd_grad(loop, Th, h=1e-6):
    out = np.zeros_like(Th)
    for i in range(Th.shape[0]):
        for j in range(3):
            E = np.zeros_like(Th); E[i, j] = h
            out[i, j] = (loop.J(Th + E) - loop.J(Th - E)) / (2 * h)
    return out


def test_plan_is_a_vertex_with_the_prescribed_masses_and_the_single_atom_prior_is_the_sigma_modification():
    net = _net()
    exact = _exact(net)
    rng = np.random.default_rng(0)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=200)
    for m in (15, 16, 25):
        Th = CL.mixture_start(exact, m, rng, 0.3, G)
        gamma, W2 = loop.plan(Th)
        assert np.allclose(gamma.sum(axis=1), 1.0 / m) and np.allclose(gamma.sum(axis=0), loop.q)
        n_split = int(loop.split_mask(Th, gamma).sum())
        if loop.commensurate(m):
            assert n_split == 0 and np.allclose(np.sort(gamma.ravel())[::-1][:m], 1.0 / m)         # integral plan
        else:
            assert 1 <= n_split <= loop.K0 - 1                                                     # a vertex: at most K₀ − 1 split rows
        assert np.allclose(loop.leak_target(Th, gamma), m * gamma @ exact)
        # the plan's cost equals the transportation LP's optimum (total unimodularity: the assignment route loses nothing)
        lp = OT.transport_plan(loop.cost_matrix(Th), np.full(m, 1.0 / m), loop.q)
        assert abs(W2 - float(np.sum(lp * loop.cost_matrix(Th)))) < 1e-9
    # a single-atom prior: the OT leak is the σ-modification toward that atom
    single = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, SIGMA0, theta0=exact[0], n_quad=200)
    one = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact[:1], n_quad=200)
    Th = CL.cold_start(net, 9, rng)
    assert np.abs(one.drift(Th) - single.drift(Th)).max() < 1e-14 and abs(one.J(Th) - (9 * single.L_bar(Th) + SIGMA0 * 0.5 * np.sum((Th - exact[0]) ** 2 / G))) < 1e-12


def test_the_OT_leak_is_the_G_gradient_flow_of_J_OT_for_commensurate_and_incommensurate_m():
    net = _net()
    exact = _exact(net)
    rng = np.random.default_rng(1)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, q=(0.5, 0.3, 0.2), n_quad=200)
    for m in (10, 7):                                                     # 10q = (5, 3, 2) commensurate; 7q is not
        Th = CL.mixture_start(exact, m, rng, 0.4, G)
        g = loop.grad_J(Th)
        assert np.abs(g - _fd_grad(loop, Th)).max() < 1e-8 * max(1.0, np.abs(g).max())   # Danskin: the plan is locally constant (error ~10⁻¹⁰)
        assert np.abs(loop.drift(Th) + G[None, :] * g).max() < 1e-13
        assert loop.commensurate(m) == (m == 10)


def test_replicated_weighted_configuration_is_stationary_for_every_commensurate_m_and_the_floor_vanishes():
    net = _net()
    exact = _exact(net)
    for q in (None, (0.4, 0.3, 0.3)):
        loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, q=q, n_quad=256)
        atoms, res, (eps, half) = loop.stationary_atoms(exact, w0=np.ones(3) / 3, n_cont=4)
        assert res < 1e-10 and eps < half                                   # the diagonal plan is the unique optimal plan
        if q is None:
            assert eps < 1e-9 and loop.fit_error(atoms) < 1e-12            # the matching masses: the exact representation
        else:
            assert loop.fit_error(atoms) > 1e-3                             # mass-mismatched masses: the stationary residual is the bias
        reps = {}
        ms = (15, 30, 60) if q is None else (10, 20, 50)                   # the commensurate sizes of each mass vector
        for m in ms:
            R, lab = loop.replicated(atoms, m)
            gamma, _ = loop.plan(R)
            assert np.all(loop.assignment(R, gamma) == lab) and loop.split_mask(R, gamma).sum() == 0
            assert np.abs(loop.drift(R)).max() < 1e-11                     # a critical point of J_m^{OT} for every commensurate m
            reps[m] = R
        assert OT.W2_empirical(reps[ms[0]], reps[ms[2]], G) < 1e-12 and OT.W2_empirical(reps[ms[1]], reps[ms[2]], G) < 1e-12   # round-off
        assert np.abs(loop.f_hat(reps[ms[0]]) - loop.f_hat(reps[ms[2]])).max() < 1e-14      # the learned function does not depend on m
        # the splitting margin of the replicated configuration is Proposition D‴.5's with w = q (nearest-atom loop, diagonal plan)
        mg, Tw, k = PD.splitting_margin_prior(loop.nearest, loop.q, exact)
        assert np.abs(Tw - atoms).max() < 1e-8
        if q is None:
            assert abs(mg - SIGMA0) < 1e-10                                 # zero residual: margin σ₀


def test_conservation_of_the_cluster_masses_against_their_constraint_on_a_short_run():
    net = _net()
    exact = _exact(net)
    rng = np.random.default_rng(0)
    Th0 = CL.mixture_start(exact, 15, rng, 0.02, G)                     # cluster weights (6, 5, 4)/15 for this seed
    lab0, w0 = CL.cluster_weights(Th0, exact, G)
    assert np.abs(w0 - 1.0 / 3).max() > 0.05
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, 10 * SIGMA0, exact, n_quad=128)     # a stronger leak shortens the transit
    rec_n = loop.coupled_run_ot(Th0, [], 40.0, 0.05, rec_every=40, plan_every=5, leak='nearest')
    rec_o = loop.coupled_run_ot(Th0, [], 40.0, 0.05, rec_every=40, plan_every=5, leak='ot')
    W_n, W_o = rec_n['w'][15], rec_o['w'][15]
    assert np.all(W_n == w0[None, :])                                     # conserved under the nearest-atom leak
    assert np.allclose(W_o[-1], 1.0 / 3) and rec_o['transit'][15][0] > 0 and rec_o['transit'][15][-1] == 0.0
    assert np.all(rec_o['n_split'][15] == 0)                              # commensurate: every plan integral
    # the mass in transit at the start is the excess of the over-weighted clusters
    assert abs(rec_o['transit'][15][0] - 0.5 * np.abs(w0 - 1.0 / 3).sum()) < 1e-12
    assert rec_o['fit'][15][-1] < rec_n['fit'][15][-1]                    # the mass constraint removes the weight bias of the fit


def test_forest_bound_and_the_exact_W2_between_empirical_measures():
    net = _net()
    exact = _exact(net)
    rng = np.random.default_rng(3)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, q=(0.5, 0.3, 0.2), n_quad=64)
    for m in (7, 11, 16, 23):
        Th = CL.cold_start(net, m, rng)
        gamma, _ = loop.plan(Th)
        assert int(loop.split_mask(Th, gamma).sum()) <= loop.K0 - 1
    A = CL.cold_start(net, 6, rng); B = CL.cold_start(net, 12, rng); C = CL.cold_start(net, 9, rng)
    assert OT.W2_empirical(A, A, G) < 1e-12 and OT.W2_empirical(A, np.repeat(A, 3, axis=0), G) < 1e-12
    # the replication route (6 | 12) agrees with the transportation LP, and W² is symmetric
    lp = OT.transport_plan(OT._cost_matrix(A, B, G), np.full(6, 1 / 6), np.full(12, 1 / 12))
    assert abs(OT.W2_empirical(A, B, G) - float(np.sum(lp * OT._cost_matrix(A, B, G)))) < 1e-10
    assert abs(OT.W2_empirical(A, C, G) - OT.W2_empirical(C, A, G)) < 1e-10          # 6 ∤ 9: the LP route
    assert np.sqrt(OT.W2_empirical(A, C, G)) <= np.sqrt(OT.W2_empirical(A, B, G)) + np.sqrt(OT.W2_empirical(B, C, G)) + 1e-10


# ---------------------------------------------------------------------- Corollary D‴.14: the collapsed subspace of the real loop
def _real_step(net, Th, x, t, leak_target, ke=25.0, sig=SIGMA0 / 25.0, dt=1e-3):
    """One Euler step of the parameters and one RK4 step of the plant of the real closed loop (CL) with the leak target leak_target(Th)."""
    xr, xrd = X_REF(t), np.array([2.0 * np.cos(t)])
    e = x - xr
    u = xrd - ke * e - net.f_hat(x, Th)
    f = lambda y: F_VEC(y)
    k1 = f(x) + u; k2 = f(x + 0.5 * dt * k1) + u; k3 = f(x + 0.5 * dt * k2) + u; k4 = f(x + dt * k3) + u
    Th_new = Th + dt * (G[None, :] * net.grad_g(x, e, Th) - sig * (Th - leak_target(Th)))
    return Th_new, x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def test_collapsed_subspace_of_the_real_loop_is_invariant_and_carries_the_weighted_K_atom_loop():
    """Corollary D‴.14(i): on a configuration collapsed by class with masses q (OT leak) or any weights inside the cells (nearest-atom
    leak) the real loop (CL) is the weighted K₀-atom loop: the collapsed m-system and the K₀-atom system stay identical as measures."""
    net = _net()
    exact = _exact(net)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=64)
    atoms = exact + np.array([[0.3, -0.2, 0.1], [-0.1, 0.2, 0.0], [0.2, 0.1, -0.2]])          # off the stationary point: a moving orbit
    # OT leak, masses q = 1/3: the 15-particle lift (5 per class) against the 3-atom loop (m = 3 is the weighted loop for uniform q)
    Th15, lab = loop.replicated(atoms, 15)
    Th3 = atoms.copy()
    x15, x3, t = X_REF(0.0) + 0.05, X_REF(0.0) + 0.05, 0.0
    for _ in range(300):
        Th15, x15 = _real_step(net, Th15, x15, t, loop.leak_target)
        Th3, x3 = _real_step(net, Th3, x3, t, loop.leak_target)
        t += 1e-3
    assert np.abs(Th15 - Th3[lab]).max() < 1e-12 and np.abs(x15 - x3).max() < 1e-12       # invariant, and the same dynamics
    assert np.abs(Th3 - atoms).max() > 1e-4                                                 # the atoms moved from their start (a genuine motion)
    # nearest-atom leak, unequal weights w = (6, 5, 4)/15: the collapsed 15-system stays collapsed (any w inside the cells)
    Thw = np.repeat(atoms, [6, 5, 4], axis=0); labw = np.repeat(np.arange(3), [6, 5, 4])
    xw, t = X_REF(0.0) + 0.05, 0.0
    for _ in range(300):
        Thw, xw = _real_step(net, Thw, xw, t, loop.nearest.leak_target)
        t += 1e-3
    for k in range(3):
        assert np.abs(Thw[labw == k] - Thw[labw == k][0]).max() < 1e-12
    # the OT leak with class masses w ≠ q is NOT invariant: the plan assigns one particle of the over-weighted class elsewhere
    Thq = np.repeat(atoms, [6, 5, 4], axis=0)
    T = loop.leak_target(Thq)
    lab_assigned = loop.assignment(Thq)
    assert np.sum(lab_assigned != labw) == 1 and np.abs(T[labw == 0] - T[labw == 0][0]).max() > 1.0


def test_the_lifted_orbit_measure_does_not_depend_on_m():
    """Corollary D‴.14(ii): two commensurate sizes started on their collapsed lifts with masses q follow the same orbit measure under the
    OT leak; with the nearest-atom leak the lifts of two different weight vectors follow different orbits."""
    net = _net()
    exact = _exact(net)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=64)
    atoms = exact + np.array([[0.2, -0.1, 0.1], [-0.1, 0.1, 0.0], [0.1, 0.1, -0.1]])
    Th15, lab15 = loop.replicated(atoms, 15); Th30, lab30 = loop.replicated(atoms, 30)
    x15, x30, t = X_REF(0.0).copy(), X_REF(0.0).copy(), 0.0
    for _ in range(200):
        Th15, x15 = _real_step(net, Th15, x15, t, loop.leak_target)
        Th30, x30 = _real_step(net, Th30, x30, t, loop.leak_target)
        t += 1e-3
    assert OT.W2_empirical(Th15, Th30, G) < 1e-12 and np.abs(x15 - x30).max() < 1e-12            # round-off
    assert np.abs(net.f_hat_batch(np.linspace(-2, 2, 50)[:, None], Th15) - net.f_hat_batch(np.linspace(-2, 2, 50)[:, None], Th30)).max() < 1e-13
    # nearest-atom leak: the lifts of two weight vectors, (6, 5, 4)/15 and (5, 5, 5)/15, are different measures from the start and stay so
    Tha = np.repeat(atoms, [6, 5, 4], axis=0); Thb = np.repeat(atoms, [5, 5, 5], axis=0)
    xa, xb, t = X_REF(0.0).copy(), X_REF(0.0).copy(), 0.0
    for _ in range(200):
        Tha, xa = _real_step(net, Tha, xa, t, loop.nearest.leak_target)
        Thb, xb = _real_step(net, Thb, xb, t, loop.nearest.leak_target)
        t += 1e-3
    da, db = Tha[[0, 6, 11]] - atoms, Thb[[0, 5, 10]] - atoms                                # the class atoms' displacements along the two weighted loops
    assert np.abs(da - db).max() > 1e-8 and np.abs(xa - xb).max() > 1e-9                  # different weights, different orbits (from the same start)


# ---------------------------------------------------------------------- Proposition D‴.15: the transverse spectrum of a lifted orbit
def test_transverse_spectrum_of_the_lifted_orbit_is_the_ghost_spectrum_and_the_collective_block_is_the_K_atom_loop():
    """Proposition D‴.15 on the matching masses: the lift of the 3-atom orbit to m = 6 is a periodic orbit; each class's separation
    subspace is invariant under its monodromy and carries the ghost particle's monodromy (exponents σ₀/k_e + O(k_e⁻²) at zero residual);
    the complement carries the 3-atom loop's monodromy."""
    from src.core import lifted_orbit as LO
    net = _net()
    exact = _exact(net)
    KE = 25.0
    x_ref_dot = lambda t: np.array([2.0 * np.cos(t)])
    clK = LO.ClassLeakLoop(net, F_VEC, X_REF, x_ref_dot, PER, KE, G, SIGMA0, exact)
    zK, DPK, res = clK.shoot(np.concatenate([X_REF(0.0), exact.ravel()]), 3, max_iter=6, tol=1e-10)
    assert res < 1e-10                                                        # the exact representation with e = 0 is the 3-atom orbit
    labels = np.array([0, 0, 1, 1, 2, 2])
    cl6 = LO.ClassLeakLoop(net, F_VEC, X_REF, x_ref_dot, PER, KE, G, SIGMA0, exact[labels])
    z6 = LO.lift(zK, labels, 1, 3)
    assert np.linalg.norm(cl6.poincare(z6, 6) - z6) < 1e-10                    # the lift is a periodic orbit of the 6-loop
    DP6 = LO.monodromy_fd(cl6, z6, 6)
    for k, B in LO.class_separation_bases(labels, 1, 3).items():
        lam, ex, defect = LO.compression_spectrum(DP6, B, PER)
        M, exg, resg = LO.ghost_monodromy(clK, zK, 3, k, exact[k])
        assert defect < 1e-6 and resg < 1e-9
        assert np.abs(ex - exg).max() < 1e-6 and np.abs(ex - SIGMA0 / KE).max() < 1e-5   # σ₀/k_e = 2.83·10⁻⁴ at zero residual
    lamc, exc, dc = LO.compression_spectrum(DP6, LO.collective_basis(labels, 1, 3), PER)
    exK = np.sort(-np.log(np.abs(np.linalg.eigvals(DPK))) / PER)
    assert dc < 1e-6 and np.abs(exc[:-1] - exK[:-1]).max() < 1e-6              # the parameter exponents (the error channel's is at the FD floor)


# ---------------------------------------------------------------------- Proposition D‴.17: J_m^{OT} is a Lyapunov function of the OT dynamics
def test_J_OT_is_a_Lyapunov_function_of_the_OT_dynamics_and_decreases_from_a_tie_along_either_active_plan():
    """Proposition D‴.17: (𝔄_m^{OT}) is the G-gradient flow of J_m^{OT} wherever the plan is unique (the energy decreases at the rate −‖drift‖²_{G⁻¹})
    and J_m^{OT} is non-increasing along a run; at a tie between two optimal plans — a configuration built by symmetry — J_m^{OT} is the minimum
    over the plans (equal for both) and decreases to first order along the drift of either active plan."""
    net = _net()
    exact = _exact(net)
    rng = np.random.default_rng(11)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, 10 * SIGMA0, exact, n_quad=96)
    Th = CL.mixture_start(exact, 15, rng, 0.3, G)                                          # unbalanced classes: the plan reassigns particles
    dt = 0.05
    J_prev = loop.J(Th)
    for k in range(300):
        gamma, _ = loop.plan(Th); T = loop.leak_target(Th, gamma)
        k1 = loop.drift(Th, T); k2 = loop.drift(Th + 0.5 * dt * k1, T)
        Th = Th + dt * k2
        J_new = loop.J(Th)
        assert J_new <= J_prev + 1e-9 * max(1.0, abs(J_prev))                              # non-increasing (up to the step's discretisation error)
        J_prev = J_new
    gamma, _ = loop.plan(Th); T = loop.leak_target(Th, gamma)
    d = loop.drift(Th, T); h = 1e-6
    dJ = (loop.J(Th + h * d) - loop.J(Th - h * d)) / (2 * h)
    n2 = float(np.sum(d ** 2 / G[None, :]))
    assert abs(dJ + n2) < 1e-6 * max(1.0, n2)                                               # dJ/ds = −‖drift‖²_{G⁻¹} where the plan is unique
    # a tie: m = 6 with classes (2, 2, 2) to fill; three particles at atom 1 (two of them mirror images across the plane through atom 1
    # orthogonal, in the G⁻¹ metric, to atom 3 − atom 1), two at atom 2, one at atom 3: the plan must send one of the two mirror images to atom 3
    v = (exact[2] - exact[0]) / G                                                            # ⟨u, ϑ⁰₃ − ϑ⁰₁⟩_{G⁻¹} = u·v: the tie needs u·v = 0
    u = np.array([0.0, 1.0, 0.0]); u = u - (u @ v) / (v @ v) * v                              # a displacement G⁻¹-orthogonal to ϑ⁰₃ − ϑ⁰₁ (equal distances to atoms 1 and 3)
    u = 0.4 * u / np.sqrt(u @ (u / G))
    Th6 = np.vstack([exact[0] + u, exact[0] - u, exact[0] + np.array([0.3, 0.0, 0.0]), exact[1], exact[1] + np.array([0.0, 0.05, 0.0]), exact[2]])
    C = loop.cost_matrix(Th6)
    assert abs(C[0, 2] - C[1, 2]) < 1e-12 and abs(C[0, 0] - C[1, 0]) < 1e-12               # the tie
    pi_a = np.array([0, 2, 0, 1, 1, 2]); pi_b = np.array([2, 0, 0, 1, 1, 2])                # the two optimal integral plans
    J_pi = lambda pi: 6 * loop.L_bar(Th6) + 0.5 * loop.sigma0 * float(np.sum((Th6 - exact[pi]) ** 2 / G[None, :]))
    assert abs(J_pi(pi_a) - J_pi(pi_b)) < 1e-12 and abs(loop.J(Th6) - J_pi(pi_a)) < 1e-10  # J^{OT} is the minimum, attained by both
    for pi in (pi_a, pi_b):
        d_pi = loop.drift(Th6, exact[pi])                                                   # the drift of the frozen plan π
        eps = 1e-4
        assert loop.J(Th6 + eps * d_pi) < loop.J(Th6) - 0.5 * eps * float(np.sum(d_pi ** 2 / G[None, :]))   # decreases along either active plan's drift
