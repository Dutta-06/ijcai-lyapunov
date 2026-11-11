"""Proposition D‴.20 (PARADIGM.md §5d): the entropic OT leak. The semi-dual Newton plan has the prescribed marginals and agrees with
log-domain Sinkhorn; the leak is the G-gradient flow of J^ε_m = mL̄ + (σ₀m/2)OT_ε (the envelope theorem, against finite differences); its two
limits, tested at their rates (the OT barycentric projection as ε → 0, the barycentre of ν as ε → ∞, the error falling like 1/ε); at a
configuration collapsed by class the entropic plan is class-symmetric, the replicated weighted stationary configuration is a critical point,
and the finite-difference Hessian of J^ε_m leaves each class's separation subspace invariant and acts there as I ⊗ (σ₀(I − (2/ε)Cov^G_k) − R^G_k),
entrywise; where σ₀I − R^G_k ≻ 0 the over-weighted class's margin changes sign at a critical ε. Every weighted stationary point is checked
to a Newton residual of 10⁻⁹."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import entropic_leak as EL

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    exact = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    return net, exact


def test_plan_marginals_and_the_two_solvers_agree():
    rng = np.random.default_rng(0)
    for m, eps in ((15, 4.0), (60, 16.0)):
        C = rng.uniform(0, 40, (m, 3)); a = np.full(m, 1.0 / m); b = np.array([0.5, 0.3, 0.2])
        pi, f, g, it = EL.sinkhorn_plan(C, a, b, eps)                            # semi-dual Newton (three columns)
        assert np.abs(pi.sum(axis=1) - a).max() < 1e-14 and np.abs(pi.sum(axis=0) - b).max() < 1e-12
        # log-domain Sinkhorn on the transposed problem (m columns > 8: the iterative route)
        piT, fT, gT, itT = EL.sinkhorn_plan(C.T.copy(), b, a, eps, tol=1e-14)
        assert np.abs(piT.T - pi).max() < 1e-10


def test_entropic_leak_is_the_gradient_flow_of_J_eps_and_has_the_two_limits():
    net, exact = _setup()
    rng = np.random.default_rng(1)
    loop = EL.EntropicPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, eps=2.0, n_quad=96)
    Th = CL.mixture_start(exact, 9, rng, 0.5, G)
    g = loop.grad_J(Th); h = 1e-6; fd = np.zeros_like(Th)
    for i in range(9):
        for j in range(3):
            E = np.zeros_like(Th); E[i, j] = h
            fd[i, j] = (loop.J(Th + E) - loop.J(Th - E)) / (2 * h)
    assert np.abs(g - fd).max() < 1e-7 * max(1.0, np.abs(g).max())                # the envelope theorem: ∂OT_ε/∂ϑ_i = (2/m)G⁻¹(ϑ_i − T^ε_i)
    ot = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=96)
    gamma, _ = ot.plan(Th)
    lab = np.argmax(gamma, axis=1); C = ot.cost_matrix(Th)
    # the OT plan is unique here: every exchange of two particles' atoms costs strictly more (c-cyclical monotonicity with a margin)
    gaps = [C[i, lab[j]] + C[j, lab[i]] - C[i, lab[i]] - C[j, lab[j]] for i in range(9) for j in range(9) if lab[i] != lab[j]]
    assert min(gaps) > 1.0
    err = lambda eps: np.abs(EL.EntropicPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, eps=eps, n_quad=96).leak_target(Th) - ot.leak_target(Th)).max()
    assert err(0.5) > 1e-3 > 1e-4 > err(0.2)                                       # ε → 0: the OT barycentric projection, exponentially fast
    bar = lambda eps: np.abs(EL.EntropicPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, eps=eps, n_quad=96).leak_target(Th) - (Q @ exact)[None, :]).max()
    b4, b5 = bar(1e4), bar(1e5)
    assert b5 < 1e-3 and 8.0 < b4 / b5 < 12.0                                      # ε → ∞: the barycentre of ν, the error ∝ 1/ε


def test_collapsed_configuration_is_critical_and_its_separation_block_is_the_entropic_formula():
    net, exact = _setup()
    loop = EL.EntropicPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, eps=10.0, n_quad=96)
    counts = np.array([6, 5, 4]); w = counts / counts.sum(); m = int(counts.sum())
    atoms, res = EL.weighted_stationary_atoms_entropic(loop, w, exact, 10.0, Q)
    assert res < 1e-10
    lab = np.repeat(np.arange(3), counts)
    Th = atoms[lab]
    assert np.abs(loop.drift(Th)).max() < 1e-9                                      # the replicated configuration is a critical point of J^ε_m
    pi, C = loop.plan(Th)
    for k in range(3):                                                              # the m-particle plan is class-symmetric
        rows = pi[lab == k]
        assert np.abs(rows - rows[0][None, :]).max() < 1e-12
    # finite-difference Hessian of J^ε_m (whitened), restricted to the separation subspace of each class
    d = 3 * m; hh = 1e-5; H = np.zeros((d, d)); v = Th.ravel()
    gJ = lambda V: loop.grad_J(V.reshape(m, 3)).ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = hh
        H[:, j] = (gJ(v + e) - gJ(v - e)) / (2 * hh)
    sg = np.tile(SQG, m); Hw = sg[:, None] * H * sg[None, :]; Hw = 0.5 * (Hw + Hw.T)
    blocks, mins, cov = EL.separation_blocks_entropic(loop, w, atoms, 10.0, Q)
    for k in range(3):
        idx = np.flatnonzero(lab == k); nk = len(idx)
        Qb, _ = np.linalg.qr(np.column_stack([np.ones(nk), np.eye(nk)[:, :nk - 1]])); V = Qb[:, 1:]
        B = np.zeros((d, (nk - 1) * 3))
        for a_ in range(nk - 1):
            for c_ in range(3):
                B[3 * idx + c_, 3 * a_ + c_] = V[:, a_]
        CB = B.T @ Hw @ B
        assert np.abs(CB - np.kron(np.eye(nk - 1), blocks[k])).max() < 1e-8          # entrywise I ⊗ (σ₀(I − (2/ε)Cov^G_k) − R^G_k)
        assert np.abs(Hw @ B - B @ CB).max() < 1e-8                                   # the separation subspace is invariant (no coupling)
    assert mins[0] < 0 < mins[1] and mins[2] > 0                                    # at ε = 10 the over-weighted class splits


def test_the_over_weighted_class_has_a_critical_regularisation():
    net, exact = _setup()
    loop = EL.EntropicPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, eps=1.0, n_quad=96)
    w = np.array([6, 5, 4]) / 15.0
    margins, atoms = {}, exact
    for eps in (64.0, 20.0, 18.0, 17.0, 16.0):                                     # continuation in ε from above
        atoms, res = EL.weighted_stationary_atoms_entropic(loop, w, atoms, eps, Q)
        assert res < 1e-9
        a, c, s, diff, phi = loop.features(atoms)
        Rw = SQG[None, :, None] * loop.residual_curvature(loop.F - phi @ (w * a), atoms) * SQG[None, None, :]
        assert np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[0]).min() > 0             # σ₀I − R^G₁ ≻ 0: the entropic term decides
        margins[eps] = EL.separation_blocks_entropic(loop, w, atoms, eps, Q)[1][0]
    assert margins[16.0] < margins[17.0] < 0 < margins[18.0] < margins[20.0] < margins[64.0]   # split below ε_c ∈ (17, 18), stable above
