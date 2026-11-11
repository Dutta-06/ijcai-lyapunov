"""Proposition D‴.28 (PARADIGM.md §5d): the adaptive prior's floor removal is a gauge property, so it survives an approximation error. On the
three-bump instance with a fourth bump that three kernels cannot represent (ε = 0.3): the gauge transform G_wϑ^b of the best weighted
approximant is a critical point of the unregularised loss at the sampled weights, with the same function; the adaptive prior's
function-level sensitivity vanishes while the fixed prior's is σ₀Ĵ(Ã + σ₀I)⁻¹g̃ with the full Hessian Ã = H_L − R̃ (not the Gauss–Newton
block); the residual curvature has a zero amplitude row at a critical point and scales by q_k/w_k under the gauge, and the separation
eigenvalues of the adaptive prior's dynamics at the replicated configuration are −σ₀ + (q_k/w_k)eig(G^{1/2}R_kG^{1/2}); and the adaptive
prior's dynamics (η = ∞ and σ₀) converges to the gauge transform."""
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
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3
X_REF = lambda t: np.array([2.0 * np.sin(t)])
EPS = 0.3


def _target(eps):
    return lambda X: (sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
                      + eps * np.exp(-(np.asarray(X, float) - 0.8) ** 2 / (2 * 0.3 ** 2)))


def _setup(n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    th = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    for k in range(1, 13):                                                  # continue the critical point in ε
        loop = OT.OTPriorLoop(net, _target(EPS * k / 12), X_REF, PER, G, SIGMA0, th, n_quad=n_quad)
        th, res, fit = FB.exact_weighted_atoms(loop, Q, th, n_cont=1, w0=Q, fit_tol=1.0)
    loop = OT.OTPriorLoop(net, _target(EPS), X_REF, PER, G, SIGMA0, th, n_quad=n_quad)   # the matching prior ϑ⁰ = ϑ^b(q)
    return net, loop, th, res, fit


def _gauge(net, tb, w):
    a, c, s = net.physical(tb)
    return net.unconstrained(Q[:, None] * a / np.asarray(w, float)[:, None], c, s)


def _fhat(loop, T, w):
    a, c, s, d, phi = loop.features(T)
    return phi @ (np.asarray(w, float) * a)


def test_the_gauge_transform_of_the_best_approximant_is_a_critical_point_with_the_same_function():
    net, loop, tb, res, fit = _setup()
    assert res < 1e-12 and fit > 0.01                                          # a critical point, and not an exact representation
    for w in (np.array([0.42, 0.30, 0.28]), np.array([0.27, 0.38, 0.35])):
        tg = _gauge(net, tb, w)
        assert np.abs(FB.unregularised_grad(loop, w, tg)).max() < 1e-12        # ∇L_w(G_wϑ^b) = 0
        assert np.abs(_fhat(loop, tg, w) - _fhat(loop, tb, Q)).max() < 1e-14   # f̂_w∘G_w = f̂_q (the gauge identity, for every Θ)
        tc = FB.exact_weighted_atoms(loop, w, tb, n_cont=20, w0=Q, fit_tol=1.0)[0]   # the critical point continued in w from ϑ^b(q)
        assert np.abs(tc - tg).max() < 1e-9                                    # is the gauge transform (not started there)


def test_the_fixed_prior_floor_has_the_full_hessian_and_the_adaptive_prior_none():
    net, loop, tb, res, fit = _setup()
    h = 1e-5
    Vd = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    Fa, Ff, gt = [], [], []
    for c in range(2):
        o = {}
        for sg in (+1, -1):
            w = Q + sg * h * Vd[:, c]
            ta = FB.exact_weighted_atoms(loop, w, tb, n_cont=1, w0=w, fit_tol=1.0)[0]       # Newton from ϑ^b(q), not from the gauge transform
            tf = PD.weighted_stationary_atoms_prior(loop.nearest, w, tb)[0]
            o[sg] = (_fhat(loop, ta, w), _fhat(loop, tf, w), _gauge(net, tb, w))
        Fa.append((o[1][0] - o[-1][0]) / (2 * h)); Ff.append((o[1][1] - o[-1][1]) / (2 * h)); gt.append(((o[1][2] - o[-1][2]) / (2 * h)).ravel())
    Fa, Ff = np.array(Fa).T, np.array(Ff).T
    Dh = np.kron(np.sqrt(Q), np.sqrt(1.0 / G))
    gt = -Dh[:, None] * np.array(gt).T                                         # g̃ = −D^{1/2}∂_w(G_wϑ^b)
    gradL = lambda v: (Q[:, None] * FB.unregularised_grad(loop, Q, v.reshape(3, 3))).ravel()
    A = np.zeros((9, 9)); e = 1e-6
    for j in range(9):
        E = np.zeros(9); E[j] = e
        A[:, j] = (gradL(tb.ravel() + E) - gradL(tb.ravel() - E)) / (2 * e)
    At = 0.5 * (A + A.T) / Dh[:, None] / Dh[None, :]                           # Ã, the whitened unregularised Hessian
    Jh = FB.whitened_jacobian(loop, Q, tb); HL = Jh.T @ Jh / loop.N
    F_full = SIGMA0 * Jh @ np.linalg.solve(At + SIGMA0 * np.eye(9), gt)
    F_gn = SIGMA0 * Jh @ np.linalg.solve(HL + SIGMA0 * np.eye(9), gt)
    assert np.linalg.norm(Fa) < 1e-8 * np.linalg.norm(Ff)                      # no function-level sensitivity under the adaptive prior
    assert np.linalg.norm(F_full - Ff) < 1e-6 * np.linalg.norm(Ff)             # the fixed prior's: σ₀Ĵ(Ã + σ₀I)⁻¹g̃ …
    assert np.linalg.norm(F_gn - Ff) > 1e-3 * np.linalg.norm(Ff)               # … not D‴.24's form with the Gauss–Newton block


def test_the_residual_curvature_has_no_amplitude_row_and_scales_with_the_gauge():
    net, loop, tb, res, fit = _setup()
    w = np.array([0.42, 0.30, 0.28])
    tg = _gauge(net, tb, w)
    r = loop.F - _fhat(loop, tb, Q)
    R_q = loop.residual_curvature(r, tb)
    R_w = loop.residual_curvature(loop.F - _fhat(loop, tg, w), tg)
    assert np.abs(R_q[:, 0, :]).max() < 1e-9 * np.abs(R_q).max()             # the amplitude row vanishes at a critical point (to the FD error)
    assert np.abs(R_q).max() > 1e-5                                           # while the residual curvature does not
    assert np.abs(R_w - (Q / w)[:, None, None] * R_q).max() < 1e-9 * np.abs(R_q).max()   # R_k(G_wϑ^b) = (q_k/w_k)R_k(ϑ^b)
    marg = SIGMA0 - loop.kappa_R(loop.F - _fhat(loop, tg, w), tg)
    sqG = np.sqrt(G)
    lam = np.array([np.linalg.eigvalsh(sqG[:, None] * (0.5 * (R_q[k] + R_q[k].T)) * sqG[None, :])[-1] for k in range(3)])
    assert np.all(lam >= -1e-12) and np.abs(marg - (SIGMA0 - (Q / w) * lam)).max() < 1e-10   # the margins σ₀ − (q_k/w_k)λ⁺_k
    # … and the dynamics': the Jacobian of the averaged adaptive drift (η = ∞) at the replicated configuration, counts (3, 2, 2), on the
    # separations (Helmert contrasts per class) has the eigenvalues −σ₀ + (q_k/w_k)eig(G^{1/2}R_k(q)G^{1/2}), each n_k − 1 times
    counts = np.array([3, 2, 2]); m = 7; wc = counts / m; lab = np.repeat(np.arange(3), counts)
    tgc = _gauge(net, tb, wc)
    apl = AP.AdaptivePriorLoop(loop, np.inf, plan='nearest')
    Th = tgc[lab].copy()
    gamma, _ = apl.plan(Th, tgc)
    fun = lambda v: apl.drift(v.reshape(m, 3), tgc, gamma)[0].ravel()
    v0 = Th.ravel(); J = np.zeros((3 * m, 3 * m))
    for j in range(3 * m):
        e = np.zeros(3 * m); hh = 1e-7 * max(1.0, abs(v0[j])); e[j] = hh
        J[:, j] = (fun(v0 + e) - fun(v0 - e)) / (2 * hh)
    cols = []
    for k in range(3):
        idx = np.where(lab == k)[0]
        for a_ in range(1, len(idx)):
            c = np.zeros(m); c[idx[:a_]] = 1.0; c[idx[a_]] = -a_; c /= np.linalg.norm(c)
            for jj in range(3):
                v = np.zeros(3 * m); v[np.arange(m) * 3 + jj] = c; cols.append(v)
    S = np.array(cols).T
    es = np.sort(np.linalg.eigvals(S.T @ J @ S).real)
    ev_q = [np.linalg.eigvalsh(sqG[:, None] * (0.5 * (R_q[k] + R_q[k].T)) * sqG[None, :]) for k in range(3)]
    pred = np.sort(np.concatenate([np.repeat(-SIGMA0 + (Q[k] / wc[k]) * ev_q[k], counts[k] - 1) for k in range(3)]))
    assert np.abs(J @ S - S @ (S.T @ J @ S)).max() < 1e-7 * np.abs(J @ S).max()   # the separations are invariant
    assert np.abs(es - pred).max() < 1e-9                                       # the dynamics' separation eigenvalues


def test_the_adaptive_prior_converges_to_the_gauge_transform():
    from scipy.integrate import solve_ivp
    net, loop, tb, res, fit = _setup()
    counts = np.array([3, 2, 2]); m = int(counts.sum()); w = counts / m
    tg = _gauge(net, tb, w)
    rng = np.random.default_rng(2)
    Th0 = np.repeat(tb, counts, axis=0) + 0.01 * rng.standard_normal((m, 3)) * np.sqrt(G)[None, :]
    # η = σ₀: the prior in the state, started at ϑ^b(q); slower (the weakest collective rate λ^−), at the gauge transform by s = 1.2·10⁵
    aps = AP.AdaptivePriorLoop(loop, SIGMA0, plan='nearest')
    gs, _ = aps.plan(Th0, tb)
    labs = np.argmax(gs, axis=1)

    def rhs_s(t, v):
        d1, d2 = aps.drift(v[:3 * m].reshape(m, 3), v[3 * m:].reshape(3, 3), gs)
        return np.concatenate([d1.ravel(), d2.ravel()])
    sol = solve_ivp(rhs_s, (0.0, 1.2e5), np.concatenate([Th0.ravel(), tb.ravel()]), method='LSODA', rtol=1e-10, atol=1e-12)
    Ths = sol.y[:3 * m, -1].reshape(m, 3)
    assert np.abs(np.array([Ths[labs == k].mean(axis=0) for k in range(3)]) - tg).max() < 1e-6
    assert np.abs(loop.f_hat(Ths) - _fhat(loop, tb, Q)).max() < 1e-6
    apl = AP.AdaptivePriorLoop(loop, np.inf, plan='nearest')
    gamma, _ = apl.plan(Th0, tb)
    lab = np.argmax(gamma, axis=1)
    assert np.array_equal(np.bincount(lab, minlength=3), counts)

    def rhs(t, v):
        Th = v.reshape(m, 3)
        Z = np.array([Th[lab == k].mean(axis=0) for k in range(3)])              # η = ∞: the prior is the class means
        return apl.drift(Th, Z, gamma)[0].ravel()
    sol = solve_ivp(rhs, (0.0, 4e4), Th0.ravel(), method='LSODA', rtol=1e-10, atol=1e-12)
    Th = sol.y[:, -1].reshape(m, 3)
    means = np.array([Th[lab == k].mean(axis=0) for k in range(3)])
    assert np.array_equal(np.argmax(apl.plan(Th, means)[0], axis=1), lab)      # the cells unchanged along the way
    assert np.abs(means - tg).max() < 1e-4 and np.abs(means - tb).max() > 1.0  # at the gauge transform, not at ϑ^b(q)
    assert np.abs(loop.f_hat(Th) - _fhat(loop, tb, Q)).max() < 1e-6          # with the approximant's function
