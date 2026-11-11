"""
Proposition D‴.32 (PARADIGM.md §5d): the cohesion leak — the leak on the unexcited subspace, without a prior.

Checked here: (i) the cohesion drift is minus the G-gradient of J^coh_m = mL̄ + (σ₀m/2)∫∫W dμdμ for both kernels, J^coh_m decreases along a
run, the leak is internal (its forces sum to zero) and equals σ₀ times the kernel mass times the mean-shift vector; (ii)–(iii) at the lift of
the exact representation at sampled weights the compact kernel's drift vanishes and the Jacobian's spectrum is −σ₀w_k on the separations
and the unregularised K-atom flow's on the class means; with residual curvature (the best one-atom fit) the separation block is
−(σ₀w_k − R^G_k), so below λ⁺ a class splits; (v) twins merge under the Gaussian kernel and persist under the compact one beyond its range;
(vi) in the real loop the lift is a periodic solution and a separation decays over a period by exactly e^{−σ₀w_kT₀/k_e}.
"""
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import cohesion_leak as CO
from src.core import adaptive_floquet as AF
from src.core.chaos_learning import AveragedSlowLoop

NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
A_B = np.array([b[0] for b in BUMPS]); C_B = np.array([b[1] for b in BUMPS]); S_B = np.array([b[2] for b in BUMPS])
SIGMA0 = 0.15 * 0.047162567584665846
PER = 2 * np.pi


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def setup(n_quad=128):
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    return net, AveragedSlowLoop(net, f_vec, x_ref, PER, G, SIGMA0, n_quad=n_quad)


def exact_rep(net, w):
    return net.unconstrained((A_B / np.asarray(w, float))[:, None], C_B[:, None], S_B)


def test_gradient_flow_internal_force_and_mean_shift():
    net, base = setup()
    rng = np.random.default_rng(3)
    Th = net.unconstrained(rng.uniform(-1, 1, (7, 1)), rng.uniform(-2, 2, (7, 1)), rng.uniform(0.4, 1.0, 7))
    Th[3] = Th[2] + 0.3 * rng.standard_normal(3) * np.sqrt(G)                 # a few pairs within range
    Th[5] = Th[4] + 0.5 * rng.standard_normal(3) * np.sqrt(G)
    for kind, R in (('compact', 1.2), ('gauss', 0.8)):
        loop = CO.CohesionLoop(base, R, kind)
        d = loop.drift(Th)
        h = 1e-6
        grad = np.zeros_like(Th)
        for i in range(Th.shape[0]):
            for j in range(3):
                E = np.zeros_like(Th); E[i, j] = h
                grad[i, j] = (loop.J(Th + E) - loop.J(Th - E)) / (2 * h)
        assert np.abs(d + G[None, :] * grad).max() < 1e-8 * max(1.0, np.abs(d).max())
        F = CO.cohesion_force(Th, G, R, SIGMA0, kind)
        assert np.abs(F.sum(axis=0)).max() < 1e-15                           # internal: the barycentre is not moved by the leak
        D, d2 = CO.pair_d2(Th, G)
        ps = CO.psi(d2 / R ** 2, kind)
        mass = ps.mean(axis=1)
        mshift = (ps @ Th) / ps.sum(axis=1)[:, None] - Th
        assert np.abs(F - SIGMA0 * mass[:, None] * mshift).max() < 1e-15      # σ₀ × kernel mass × mean-shift vector
        Js = [loop.J(Th)]
        T = Th.copy()
        for _ in range(4):
            T = loop.run(T, 5.0, 0.05)
            Js.append(loop.J(T))
        assert np.all(np.diff(Js) < 0)                                          # J^coh_m decreases along the flow


def test_exact_representation_stationary_with_the_predicted_spectrum():
    net, base = setup()
    loop = CO.CohesionLoop(base, 1.0, 'compact')
    for counts in ([3, 2, 2], [4, 3, 3]):                                    # inside W: w_k > |A_k|/a_max = (0.25, 0.175, 0.2)
        counts = np.array(counts); w = counts / counts.sum()
        at = exact_rep(net, w)
        D, d2 = CO.pair_d2(at, G)
        assert np.sqrt(d2[np.triu_indices(3, 1)]).min() > 1.0                  # the classes beyond the range
        Th = np.repeat(at, counts, axis=0)
        assert np.abs(loop.drift(Th)).max() < 1e-14                           # stationary: no leak bias at the sampled weights
        ev = np.sort(np.linalg.eigvals(loop.jacobian(Th)).real)
        sep, coll = CO.predicted_spectrum(base, at, counts, SIGMA0)
        assert np.abs(ev - np.sort(np.concatenate([sep, coll]))).max() < 1e-9
        # the separations at −σ₀w_k exactly (no residual curvature at an exact representation), with multiplicity 3(n_k − 1)
        assert np.abs(sep - np.sort(np.repeat(-SIGMA0 * w, 3 * (counts - 1)))).max() < 1e-12
        # the collective block is the unregularised one: its slowest rate is far below σ₀ (no leak on the class means)
        assert -coll[-1] < 0.2 * SIGMA0


def test_splitting_margin_with_residual_curvature():
    """At the best one-atom fit (a critical point of the unregularised loss with r ≠ 0), the separation block is −(σ₀ − R^G): stable above
    λ⁺ = λ_max(R^G), split below."""
    net, base = setup()
    x0 = net.unconstrained(np.array([[1.0]]), np.array([[-1.5]]), np.array([0.5])).ravel()

    def fun(v):
        return base.L_bar(v.reshape(1, 3))

    def grad(v):
        T = v.reshape(1, 3)
        return -base.grad_weighted(base.residual(T), T).ravel()
    b1 = minimize(fun, x0, jac=grad, method='BFGS', options=dict(gtol=1e-14)).x.reshape(1, 3)
    assert np.abs(grad(b1.ravel())).max() < 1e-10
    lam = float(base.kappa_R(base.residual(b1), b1)[0])
    assert 0 < lam < SIGMA0
    m = 4
    Th = np.repeat(b1, m, axis=0)
    for sig in (SIGMA0, 0.5 * lam):
        loop = CO.CohesionLoop(base, 1.0, 'compact', sigma0=sig)
        ev = np.sort(np.linalg.eigvals(loop.jacobian(Th)).real)
        sep, coll = CO.predicted_spectrum(base, b1, np.array([m]), sig)
        assert np.abs(ev - np.sort(np.concatenate([sep, coll]))).max() < 1e-9
        assert abs(sep[-1] - (lam - sig)) < 1e-9                                # the top separation eigenvalue λ⁺ − σ₀w (w = 1)
    assert sep[-1] > 0                                                          # below λ⁺ the class splits


def test_twins_contract_at_the_predicted_rate_and_persist_beyond_the_compact_range():
    """(v) as restated in round 23 (D3u-1): twins of bump 3 (the same centre and width, its amplitude shared 0.6 : 0.4) with no third class
    within reach contract on their gauge curve at the rate σ₀ψ(Δα²/R²)(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²) (a′ = 1 − (a/a_max)², the fit
    relaxed quasi-statically): the measured ln Δα(40)/Δα(160) is ∫ of the rate along the run to 1% at 10σ₀; and beyond the compact kernel's
    range they are a neutral pair."""
    net, base = setup()
    counts = np.array([2, 2, 1, 1]); w = counts / 6
    amps = np.array([A_B[0] / w[0], A_B[1] / w[1], 0.6 * A_B[2] / w[2], 0.4 * A_B[2] / w[3]])
    at = net.unconstrained(amps[:, None], np.r_[C_B, C_B[2]][:, None], np.r_[S_B, S_B[2]])
    Th0 = np.repeat(at, counts, axis=0)
    assert base.fit_error(Th0) < 1e-12                                          # the twins represent bump 3 exactly
    d0 = np.sqrt(np.sum((at[2] - at[3]) ** 2 / G))
    assert d0 > 1.0
    sig = 10 * SIGMA0
    rec = []

    def cb(s_, T):
        a_, c_, s__ = net.physical(T[4:6])
        rec.append((s_, T[4, 0] - T[5, 0], a_[0, 0], a_[1, 0]))
    CO.CohesionLoop(base, 1.0, 'gauss', sigma0=sig).run(Th0, 160.0, 0.25, rec_every=4, callback=cb)
    r = np.array(rec)
    da = lambda a: 1 - (a / NET_KW['a_max']) ** 2
    kap = sig * np.exp(-0.5 * r[:, 1] ** 2 / G[0]) * (w[2] * da(r[:, 2]) + w[3] * da(r[:, 3])) ** 2 / (w[2] * da(r[:, 2]) ** 2 + w[3] * da(r[:, 3]) ** 2)
    sel = r[:, 0] >= 40.0
    pred = np.sum(0.5 * (kap[sel][1:] + kap[sel][:-1]) * np.diff(r[sel, 0]))
    meas = np.log(abs(r[sel, 1][0]) / abs(r[sel, 1][-1]))
    assert pred > 2.0 and abs(meas / pred - 1) < 0.01
    Tc = CO.CohesionLoop(base, 1.0, 'compact', sigma0=sig).run(Th0, 400.0, 0.25)
    assert abs(np.sqrt(np.sum((Tc[4] - Tc[5]) ** 2 / G)) - d0) < 1e-10        # compact, beyond the range: a neutral, persistent pair


def test_the_kernels_range_and_scale_are_pinned():
    """D3u-5: the compact kernel attracts a pair at 0.99R and not at 1.01R; the Gaussian kernel's pair force at the distance R is e^{−1/2}
    times its linear part; both at a G⁻¹ distance (the whitened coordinates)."""
    R, sig = 1.3, 0.7
    for kind in ('compact', 'gauss'):
        for frac in (0.5, 0.99, 1.01, 2.0):
            d = frac * R
            Th = np.array([[0.0, 0.0, 0.0], [d * np.sqrt(G[0]), 0.0, 0.0]])
            F = CO.cohesion_force(Th, G, R, sig, kind)
            lin = sig / 2 * d * np.sqrt(G[0])                                    # (σ₀/m)·ψ(0)·|Δϑ| for m = 2
            ratio = F[0, 0] / lin                                                  # particle 0 pulled toward particle 1
            if kind == 'compact':
                assert (ratio > 1e-6) == (frac < 1.0) and abs(ratio - max(1 - frac ** 2, 0.0) ** 3) < 1e-12
            else:
                assert abs(ratio - np.exp(-0.5 * frac ** 2)) < 1e-12
            assert abs(F[0, 0] + F[1, 0]) < 1e-15


def test_single_linkage_clusters_on_a_known_configuration():
    """D3u-5: `clusters` counts three groups whose members are within the threshold and whose gaps exceed it, links a chain of steps below
    the threshold into one cluster, and measures distances in the G⁻¹ metric."""
    rng = np.random.default_rng(1)
    centres = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.5]]) * np.sqrt(G)[None, :]
    Th = np.repeat(centres, 4, axis=0) + 0.01 * rng.standard_normal((12, 3)) * np.sqrt(G)[None, :]
    lab = CO.clusters(Th, G, 0.25)
    assert lab.max() + 1 == 3 and all(len(set(lab[4 * k:4 * k + 4])) == 1 for k in range(3))
    chain = np.array([[0.2 * i * np.sqrt(G[1]) * 0 + 0.0, 0.2 * i * np.sqrt(G[1]), 0.0] for i in range(6)])
    assert CO.clusters(chain, G, 0.25).max() == 0
    assert CO.clusters(chain, G, 0.15).max() == 5


def test_real_loop_lift_periodic_and_separation_decay_exact():
    net, base = setup()
    counts = np.array([2, 1, 1]); w = counts / 4; m = 4
    at = exact_rep(net, w)
    ke = 25.0
    cl = CO.CohesionRealLoop(net, f_vec, x_ref, x_ref_dot, PER, ke, G, SIGMA0, 1.0, 'compact')
    z0 = cl.lift(at, counts)
    z1 = AF.period_map(cl, z0, m)
    assert np.abs(z1 - z0).max() < 1e-10                                         # perfect tracking is a periodic solution
    u = np.array([0.3, -0.2, 0.1]) * 1e-4                                         # a separation of class 0 (two particles, ±u)
    dz = np.zeros_like(z0); dz[1:4] = u; dz[4:7] = -u
    z2 = AF.period_map(cl, z0 + dz, m)
    ratio = (z2 - z1)[1:4] / u
    assert np.abs(ratio - np.exp(-SIGMA0 * w[0] * PER / ke)).max() < 1e-6       # e^{−σ₀w_kT₀/k_e}: the separation's exact exponent
