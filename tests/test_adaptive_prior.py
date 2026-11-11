"""Theorem D‴.23 (PARADIGM.md §5d): adaptive priors. The joint system is the gradient flow of 𝒥 = L̄ + (σ₀/2)W²(μ_Θ, ν_Z) in the stated metric
(the particles' drift −mG∇_ϑ𝒥 and the prior's −(η/(σ₀q_j))G∇_{z_j}𝒥, against finite differences); the exact representation with the prior at its
atoms is stationary, while the fixed prior's stationary point is not stationary for the adaptive prior; the linearisation at the exact
representation has the separation rate σ₀ and the collective pairs' rates λ_j^± of pair_rates; the fixed prior's bias is σ₀(H_L + σ₀I)⁻¹x⁰ to
first order; a gated run keeps the prior frozen while mass is in transit, with 𝒥 non-increasing; the joint linearisation is self-adjoint in
the gradient-flow metric (M times the Jacobian is minus the Hessian of 𝒥); the settled gate opens only when the fixed-prior drift has
fallen below its tolerance, after the transit; the limits of (i): a non-replicated stationary point (a class split over two bumps and held
by the leak toward its barycentre, against a non-zero loss gradient — a saddle) and a permuted prior that is stationary at the exact
representation (z_σ(k) = ϑ_k for the plan's bijection σ); and the loss of the a-priori bound: with a target bump beyond the saturation the
fixed prior's K₀-atom solution converges while the adaptive prior's amplitude coordinate runs off; the projected prior (Proposition D‴.27)
keeps it bounded, and at η = ∞ its flow carries a Lloyd-consistent plan, along which the selected energy never increases while the
reduced energy min_Z𝒥 can rise."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup(shift=0.15, n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    rep = lambda sh: net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + sh] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    exact, prior = rep(0.0), rep(shift)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    return loop, exact, prior


def test_the_joint_system_is_the_gradient_flow_of_the_joint_energy():
    loop, exact, prior = _setup()
    eta = 3 * SIGMA0
    apl = AP.AdaptivePriorLoop(loop, eta)
    rng = np.random.default_rng(2)
    m = 12
    Th = np.repeat(exact, m // 3, axis=0) + 0.05 * rng.standard_normal((m, 3)) * np.sqrt(G)[None, :]
    Z = prior + 0.02 * rng.standard_normal((3, 3))
    gamma, _ = AP.plan_to(Th, Z, Q, G)
    dTh, dZ = apl.drift(Th, Z, gamma)
    J = lambda T_, Z_: apl.energy(T_, Z_)
    h = 1e-6
    gT = np.zeros_like(Th); gZ = np.zeros_like(Z)
    for i in range(m):
        for c in range(3):
            E = np.zeros_like(Th); E[i, c] = h
            gT[i, c] = (J(Th + E, Z) - J(Th - E, Z)) / (2 * h)
    for j in range(3):
        for c in range(3):
            E = np.zeros_like(Z); E[j, c] = h
            gZ[j, c] = (J(Th, Z + E) - J(Th, Z - E)) / (2 * h)
    assert np.abs(dTh + m * G[None, :] * gT).max() < 1e-7                          # dϑ_i/ds = −mG∇_{ϑ_i}𝒥
    assert np.abs(dZ + (eta / (SIGMA0 * Q[:, None])) * G[None, :] * gZ).max() < 1e-7  # dz_j/ds = −(η/(σ₀q_j))G∇_{z_j}𝒥


def test_the_exact_representation_is_stationary_and_the_fixed_priors_stationary_point_is_not():
    loop, exact, prior = _setup()
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    m = 15
    Th = np.repeat(exact, m // 3, axis=0)
    gamma, _ = AP.plan_to(Th, exact, Q, G)
    dTh, dZ = apl.drift(Th, exact, gamma)
    assert np.abs(dTh).max() < 1e-14 and np.abs(dZ).max() < 1e-14 and apl.energy(Th, exact) < 1e-28
    star, res, _ = loop.stationary_atoms(prior)                                  # the fixed prior's stationary point (biased)
    Th = np.repeat(star, m // 3, axis=0)
    gamma, _ = AP.plan_to(Th, prior, Q, G)
    dTh, dZ = apl.drift(Th, prior, gamma)
    assert res < 1e-12 and np.abs(dTh).max() < 1e-12                             # the particles are at rest …
    assert np.abs(dZ).max() > 1e-4                                               # … but the prior moves toward the class atoms


def test_the_linearisation_has_the_separation_rate_and_the_collective_pair_rates():
    loop, exact, prior = _setup()
    for eta in (SIGMA0, 10 * SIGMA0):
        apl = AP.AdaptivePriorLoop(loop, eta)
        m = 6
        Th0 = np.repeat(exact, m // 3, axis=0)
        gamma, _ = AP.plan_to(Th0, exact, Q, G)
        def F(v):
            T_ = v[:3 * m].reshape(m, 3); Z_ = v[3 * m:].reshape(3, 3)
            a, b = apl.drift(T_, Z_, gamma)
            return np.concatenate([a.ravel(), b.ravel()])
        v0 = np.concatenate([Th0.ravel(), exact.ravel()]); d = v0.size; h = 1e-6
        Jac = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            Jac[:, j] = (F(v0 + e) - F(v0 - e)) / (2 * h)
        rates = np.sort(-np.linalg.eigvals(Jac).real)
        HL, hs, V = AP.gauss_newton_block(loop, exact, Q)
        lm, lp = AP.pair_rates(hs, SIGMA0, eta)
        pred = np.sort(np.concatenate([lm, lp, np.full(3 * (m - 3), SIGMA0)]))  # 9 pairs and the separations of three classes of two
        assert np.abs(rates - pred).max() < 1e-6 * max(1.0, pred.max())


def test_the_fixed_priors_bias_is_the_first_order_prediction():
    loop, exact, prior = _setup(shift=0.05)
    HL, hs, V = AP.gauss_newton_block(loop, exact, Q)
    star, res, _ = loop.stationary_atoms(prior)
    x0 = AP.whitened(prior - exact, Q, G); y = AP.whitened(star - exact, Q, G)
    pred = SIGMA0 * np.linalg.solve(HL + SIGMA0 * np.eye(9), x0)
    assert res < 1e-12 and np.linalg.norm(y - pred) < 0.05 * np.linalg.norm(y)   # first order in the prior's error (a smaller shift than the experiment's)
    assert np.linalg.norm(y) < 0.2 * np.linalg.norm(x0)                          # the stiff modes are suppressed …
    assert abs((V.T @ y)[0]) > 0.5 * abs((V.T @ x0)[0])                          # … the weakest mode keeps most of the error


def test_a_gated_run_freezes_the_prior_during_the_transit_and_decreases_the_energy():
    loop, exact, prior = _setup()
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    m = 15
    Th = np.repeat(exact, m // 3, axis=0)
    Th[0] = 0.4 * exact[0] + 0.6 * exact[1]                                      # a particle of class 0 inside cell 1: mass in transit
    gamma, _ = AP.plan_to(Th, prior, Q, G)
    assert apl.in_transit(Th, prior, gamma)
    rec = []
    apl.run(Th, prior, 300.0, 0.25, plan_every=4, rec_every=40, callback=lambda s, T_, Z_, ad: rec.append((s, Z_.copy(), ad, apl.energy(T_, Z_))), gate=True)
    Zs = [r[1] for r in rec]; ad = [r[2] for r in rec]; Js = np.array([r[3] for r in rec])
    first = ad.index(True)
    assert first > 0 and all(np.array_equal(Zs[k], prior) for k in range(first))   # frozen while in transit
    assert np.abs(Zs[-1] - prior).max() > 1e-4                                       # adapting afterwards
    assert np.all(np.diff(Js) <= 1e-12)                                              # 𝒥 non-increasing in both phases


def test_the_linearisation_is_self_adjoint_in_the_gradient_flow_metric():
    # (i)–(ii): the joint drift is −M⁻¹∇𝒥 with M = (1/m)G⁻¹ per particle and (σ₀q_j/η)G⁻¹ per prior atom, so M times its Jacobian at the exact
    # representation is minus the Hessian of 𝒥: symmetric and negative semidefinite (the symmetrised pair matrices of (ii))
    loop, exact, prior = _setup()
    for eta in (0.3 * SIGMA0, 3 * SIGMA0):
        apl = AP.AdaptivePriorLoop(loop, eta)
        m = 6
        Th0 = np.repeat(exact, m // 3, axis=0)
        gamma, _ = AP.plan_to(Th0, exact, Q, G)

        def F(v):
            T_ = v[:3 * m].reshape(m, 3); Z_ = v[3 * m:].reshape(3, 3)
            a, b = apl.drift(T_, Z_, gamma)
            return np.concatenate([a.ravel(), b.ravel()])
        v0 = np.concatenate([Th0.ravel(), exact.ravel()]); d = v0.size; h = 1e-6
        Jac = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            Jac[:, j] = (F(v0 + e) - F(v0 - e)) / (2 * h)
        M = np.concatenate([np.tile(1.0 / (m * G), m), np.concatenate([(SIGMA0 * Q[j] / eta) / G for j in range(3)])])
        S = M[:, None] * Jac
        assert np.abs(S - S.T).max() < 1e-7 * np.abs(S).max()
        assert np.linalg.eigvalsh(-0.5 * (S + S.T)).min() > -1e-9 * np.abs(S).max()


def test_the_settled_gate_opens_when_the_fixed_prior_drift_first_falls_below_its_tolerance_after_the_transit():
    loop, exact, prior = _setup()
    m, tol = 15, 1e-3
    Th = np.repeat(exact, m // 3, axis=0)
    Th[0] = 0.4 * exact[0] + 0.6 * exact[1]                                      # mass in transit at the start
    t_transit = AP.AdaptivePriorLoop(loop, SIGMA0).run(Th, prior, 400.0, 0.25, plan_every=4, rec_every=10 ** 9, gate=True)[2]
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    rec = []

    def cb(s, T_, Z_, ad):
        g_, _ = AP.plan_to(T_, Z_, Q, G)
        rec.append((s, Z_.copy(), ad, apl.fixed_prior_drift(T_, Z_, g_), apl.in_transit(T_, Z_, g_), apl.energy(T_, Z_)))
    t_open = apl.run(Th, prior, 1200.0, 0.25, plan_every=4, rec_every=4, callback=cb, gate='settled', settle_tol=tol)[2]
    first = [r[2] for r in rec].index(True)
    assert np.isfinite(t_transit) and t_open > t_transit + 50.0                  # it waits beyond the end of the transit …
    assert all(np.array_equal(r[1], prior) for r in rec[:first + 1])             # … with the prior frozen until it opens
    assert rec[first][3] < tol and not rec[first][4] and rec[first - 1][3] >= tol  # it opens when the fixed-prior drift first falls below tol
    assert np.all(np.diff([r[5] for r in rec]) <= 1e-12)                         # 𝒥 non-increasing in both phases


def _bumps(B):
    return lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in B)


def test_a_non_replicated_stationary_point_is_held_by_the_leak_and_a_permuted_prior_is_stationary():
    from scipy.optimize import least_squares
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    m = 6                                                                        # four bumps, six kernels: one class split over the last two
    Th0 = net.unconstrained(np.array([[3.0], [3.0], [-2.1], [-2.1], [2.4], [2.4]]), np.array([[-1.5], [-1.5], [0.0], [0.0], [1.3], [2.3]]),
                            np.array([0.5, 0.5, 0.4, 0.4, 0.4, 0.4]))
    loop = OT.OTPriorLoop(net, _bumps([(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.4, 1.3, 0.4), (0.4, 2.3, 0.4)]), X_REF, PER, G, SIGMA0, Th0[[0, 2, 4]], n_quad=96)
    lab = np.array([0, 0, 1, 1, 2, 2]); gam = np.zeros((m, 3)); gam[np.arange(m), lab] = 1.0 / m
    assert loop.fit_error(Th0) < 1e-12                                           # the exact six-kernel fit
    F_inf = lambda v: (loop.drift(v.reshape(m, 3), m * (gam @ AP.barycentres(v.reshape(m, 3), gam, Q))) / np.sqrt(G)[None, :]).ravel()
    Th = least_squares(F_inf, Th0.ravel(), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000).x.reshape(m, 3)
    Z = AP.barycentres(Th, gam, Q)
    leak = SIGMA0 * np.sqrt(np.sum((Th - m * (gam @ Z)) ** 2 / G[None, :], axis=1))
    unreg = np.abs(loop.drift(Th, Th) / np.sqrt(G)[None, :]).max()
    assert np.sqrt(np.sum((Th[4] - Th[5]) ** 2 / G)) > 0.5                       # the class stays split …
    assert leak[4:].min() > 1e-3 and leak[:4].max() < 1e-12 and unreg > 1e-3     # … held by the leak against a non-zero loss gradient
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    g_, _ = AP.plan_to(Th, Z, Q, G)
    dT, dZ = apl.drift(Th, Z, g_)
    assert np.array_equal(np.argmax(g_, axis=1), lab) and np.abs(dT).max() < 1e-12 and np.abs(dZ).max() < 1e-12   # stationary for the joint flow

    def Fj(v):
        x, y = apl.drift(v[:3 * m].reshape(m, 3), v[3 * m:].reshape(3, 3), g_)
        return np.concatenate([x.ravel(), y.ravel()])
    v0 = np.concatenate([Th.ravel(), Z.ravel()]); d = v0.size; h = 1e-6; J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (Fj(v0 + e) - Fj(v0 - e)) / (2 * h)
    assert np.linalg.eigvals(J).real.max() > 1e-6                                # a saddle
    loop3, exact, prior = _setup()                                               # the relabelling: a permuted prior at the exact representation
    Th3 = np.repeat(exact, 2, axis=0); Zp = exact[[1, 2, 0]]
    g3, W2 = AP.plan_to(Th3, Zp, Q, G)
    dT3, dZ3 = AP.AdaptivePriorLoop(loop3, SIGMA0).drift(Th3, Zp, g3)
    assert np.abs(Zp - exact).max() > 1.0 and W2 < 1e-24 and np.abs(dT3).max() < 1e-14 and np.abs(dZ3).max() < 1e-14


def test_the_adaptive_prior_gives_up_the_leaks_a_priori_bound():
    """No exact representation (3·1.6 > a_max). The fixed prior converges at finite parameters. Under the adaptive prior (η = σ₀) the force on
    the first amplitude, F = ⟨rφ₁⟩, tends to the positive Karush–Kuhn–Tucker multiplier F_∞ of the loss's minimum on the face a₁ = a_max — a
    bounded limit would be a critical point, where F = 0 — and the amplitude follows the run-off law (a_max − a)s → a_max²(σ₀ + η)/(4ηG_αF_∞)
    of the conservation law d/ds[α + (σ₀/η)z_α] = G_αF∂a/∂α with α − z_α slaved, so α grows like (a_max/2) ln s: the solution is unbounded."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import minimize
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    start = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    loop = OT.OTPriorLoop(net, _bumps([(1.6, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]), X_REF, PER, G, SIGMA0, start, n_quad=96)   # 3·1.6 > a_max
    q = np.ones(3) / 3

    def res_phi(a, c, s):
        phi = np.exp(-(loop.X[:, None] - c[None, :]) ** 2 / (2 * s[None, :] ** 2))
        return loop.F - phi @ (q * a), phi
    face = lambda v: 0.5 * float(np.mean(res_phi(np.concatenate([[net.a_max], v[:2]]), v[2:5], v[5:])[0] ** 2))
    v = minimize(face, [-2.1, 2.4, -1.5, 0.0, 1.5, 0.5, 0.4, 0.5], method='L-BFGS-B', bounds=[(-4, 4)] * 2 + [(-3, 3)] * 3 + [(0.3, 2.0)] * 3,
                 options=dict(ftol=1e-15, gtol=1e-12, maxiter=5000)).x
    r, phi = res_phi(np.concatenate([[net.a_max], v[:2]]), v[2:5], v[5:])
    F_inf = float(np.mean(r * phi[:, 0]))
    assert F_inf > 0.02                                                            # the loss still decreases outward on the face: a KKT point
    gamma = np.eye(3) / 3.0
    out = {}
    for eta in (0.0, SIGMA0):
        apl = AP.AdaptivePriorLoop(loop, eta)

        def rhs(t, v, apl=apl):
            dT, dZ = apl.drift(v[:9].reshape(3, 3), v[9:].reshape(3, 3), gamma)
            return np.concatenate([dT.ravel(), dZ.ravel()])
        S_end = 1e4 if eta == 0.0 else 1e5
        sol = solve_ivp(rhs, (0.0, S_end), np.concatenate([start.ravel(), start.ravel()]), method='LSODA', t_eval=[S_end / 10, S_end], rtol=1e-10, atol=1e-12)
        a, c, s_ = net.physical(sol.y[:9].T.reshape(-1, 3, 3)[-1])
        out[eta] = (sol.y[0], np.abs(rhs(S_end, sol.y[:, -1])).max(), a[0, 0], float(np.mean(res_phi(a[:, 0], c[:, 0], s_)[0] * res_phi(a[:, 0], c[:, 0], s_)[1][:, 0])))
    th, drift, a, F = out[0.0]
    assert abs(th[1] - th[0]) < 1e-6 and drift < 1e-12 and a < 3.7               # the fixed prior: converged at finite parameters
    th, drift, a, F = out[SIGMA0]
    assert abs(F / F_inf - 1) < 0.02                                               # the force tends to the face's multiplier, not to 0
    law = net.a_max ** 2 * (SIGMA0 + SIGMA0) / (4 * SIGMA0 * G[0] * F_inf)
    assert abs((net.a_max - a) * 1e5 / law - 1) < 0.03                             # the run-off law: a_max − a ~ 1/s
    assert (th[1] - th[0]) / np.log(10.0) > 0.45 * net.a_max                       # and α grows like (a_max/2) ln s


def test_a_projected_prior_keeps_the_solutions_bounded_and_is_inactive_inside_its_box():
    from scipy.integrate import solve_ivp
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    lo = np.array([-net.a_max * np.arctanh(0.9), -net.c_max * np.arctanh(0.9), -6.0]); box = (lo, -lo)
    start = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    loop = OT.OTPriorLoop(net, _bumps([(1.6, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]), X_REF, PER, G, SIGMA0, start, n_quad=96)  # 3·1.6 > a_max
    gamma = np.eye(3) / 3.0
    apl = AP.AdaptivePriorLoop(loop, np.inf, box=box)                                # η = ∞: the leak toward the projected class mean
    rhs = lambda t, v: apl.drift(v.reshape(3, 3), v.reshape(3, 3), gamma)[0].ravel()
    sol = solve_ivp(rhs, (0.0, 1e5), start.ravel(), method='LSODA', t_eval=np.geomspace(1.0, 1e5, 60), rtol=1e-10, atol=1e-12)
    Th_t = sol.y.T.reshape(-1, 3, 3)
    dist = np.sqrt(np.sum((Th_t - np.clip(Th_t, lo, -lo)) ** 2 / G[None, None, :], axis=2)).max()
    assert dist < 2.0 and sol.y[0].max() < 8.0                                       # bounded …
    assert np.abs(rhs(1e5, sol.y[:, -1])).max() < 1e-12                              # … and convergent
    free = AP.AdaptivePriorLoop(loop, np.inf)                                        # without the box the prior's coordinate passes 11 by s = 10⁴
    sol0 = solve_ivp(lambda t, v: free.drift(v.reshape(3, 3), v.reshape(3, 3), gamma)[0].ravel(), (0.0, 1e4), start.ravel(), method='LSODA',
                     rtol=1e-10, atol=1e-12)
    assert sol0.y[0, -1] > 11.0
    apl2 = AP.AdaptivePriorLoop(loop, SIGMA0, box=box)                               # η = σ₀: the projected flow, by the projected midpoint rule
    rec = []
    Th2, Z2, _ = apl2.run(start.copy(), start.copy(), 2000.0, 0.125, plan_every=8, rec_every=80,
                          callback=lambda s_, T_, Z_, ad: rec.append((T_[0, 0], Z_.copy(), apl2.energy(T_, Z_))))
    assert all(np.all((Z_ >= lo - 1e-12) & (Z_ <= -lo + 1e-12)) for _, Z_, _ in rec)  # the prior stays in its box (the clip) …
    assert Z2[0, 0] == -lo[0]                       # … and its amplitude atom reaches the face (a velocity projected at the midpoint stalls it inside at this dt)
    assert max(t for t, _, _ in rec) < 8.0 and np.all(np.diff([J for _, _, J in rec]) <= 1e-12)   # the particles bounded, 𝒥 non-increasing
    loop2, exact, prior = _setup()                                                   # the misplaced prior: the box is never reached
    star = loop2.stationary_atoms(prior)[0]
    runs = [AP.AdaptivePriorLoop(loop2, SIGMA0, box=bx).run(star.copy(), prior.copy(), 300.0, 0.25, plan_every=4, rec_every=10 ** 9)[:2] for bx in (None, box)]
    assert np.array_equal(runs[0][0], runs[1][0]) and np.array_equal(runs[0][1], runs[1][1])


def test_at_eta_infinity_the_selected_energy_decreases_while_the_reduced_energy_may_rise():
    """Proposition D‴.27(ii) at η = ∞: the flow carries a Lloyd-consistent plan (optimal for its own projected barycentres), not a jointly
    optimal one. Six particles at zero amplitude (the configuration of exp_adaptive_prior.py Part H(e), found there by a seeded search), a
    small three-bump target, the box: along the run the selected energy 𝒥(Θ, Z) never increases, while the reduced energy min_Z𝒥 (the
    minimum over the 15 balanced pairings) rises before the plan becomes the optimal one."""
    Y = np.array([[1.4296171063502774, -1.8656576987781426], [0.9186217857197763, -1.297377517589764], [1.4527156893995463, 0.16584488099636685],
                  [-0.8011524378504609, -0.3092511152093662], [-1.8867213154181481, -1.5028668940017442], [0.6824976587745213, 0.5887580462970003]])
    pair = [(0, 2), (1, 4), (3, 5)]
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    lo = np.array([-net.a_max * np.arctanh(0.9), -net.c_max * np.arctanh(0.9), -6.0]); hi = -lo
    Th0 = np.column_stack([np.zeros(6), np.sqrt(G[1]) * Y[:, 0], np.sqrt(G[2]) * Y[:, 1]])
    Z0 = np.array([Th0[list(pp)].mean(axis=0) for pp in pair])
    small = lambda X: 0.02 * F_VEC(X)
    loop = OT.OTPriorLoop(net, small, X_REF, PER, G, SIGMA0, Z0, n_quad=96)
    apl = AP.AdaptivePriorLoop(loop, np.inf, plan='ot', box=(lo, hi))

    def pairings(rem):
        if not rem:
            yield []
            return
        for b in rem[1:]:
            for pp in pairings([x for x in rem if x not in (rem[0], b)]):
                yield [(rem[0], b)] + pp
    PAIRS = list(pairings(list(range(6))))

    def reduced(Th):
        best = np.inf
        for pp in PAIRS:
            e = 0.0
            for a, b in pp:
                T = (Th[a] + Th[b]) / 2
                e += (np.sum((Th[a] - T) ** 2 / G) + np.sum((Th[b] - T) ** 2 / G)) / 6 + np.sum((T - np.clip(T, lo, hi)) ** 2 / G) / 3
            best = min(best, loop.L_bar(Th) + 0.5 * SIGMA0 * e)
        return best
    sel, red = [], []
    apl.run(Th0.copy(), Z0.copy(), 200.0, 0.25, plan_every=4, rec_every=40,
            callback=lambda s, Th, Z, ad: (sel.append(apl.energy(Th, Z)), red.append(reduced(Th))))
    sel, red = np.array(sel), np.array(red)
    assert np.all(np.diff(sel) <= 1e-12)                                       # the selected energy: non-increasing
    assert red[0] < sel[0] and red.max() > 1.2 * red[0]                        # the plan not jointly optimal; the reduced energy rises
    assert abs(sel[-1] - red[-1]) < 1e-12                                     # … until the plan is the optimal one
