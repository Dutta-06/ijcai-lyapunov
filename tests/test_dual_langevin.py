"""Proposition D‴.50: under isotropic label noise the free law's dual is a Langevin diffusion in the Bregman projection's dual objective."""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov

from src.core import mirror_free as MF
from src.core import dual_langevin as DL

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)


def small_problem(seed=0, N=3, d=5, eps=1.0):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d))
    y = Phi @ np.r_[0.6, -0.4, np.zeros(d - 2)]
    b = np.exp(rng.uniform(np.log(0.25), np.log(4.0), d))
    return Phi, y, L.Q(np.zeros(d), eps * b), b


def test_fast_free_matches_the_library():
    """η on both sides of q = c, dα/dη = Gβa′, the inverse and dw/dη, against FreeLaw."""
    for eps in (1e-6, 1e-2, 1.0, 30.0):
        q = L.Q(np.zeros(6), eps * np.array([0.25, 0.5, 1.0, 2.0, 3.0, 4.0]))
        F = DL.FastFree(L, q)
        x = np.array([1e-6, 0.01, 0.3, 1.0, 2.5, 7.0])
        assert np.allclose(F.eta(x), L.eta(x, q), rtol=1e-12, atol=0)
        assert np.allclose(F.eta(np.vstack([x, 0.5 * x])), np.vstack([L.eta(x, q), L.eta(0.5 * x, q)]), rtol=1e-12, atol=0)
        a, da = L.amp(x)
        assert np.allclose(F.dal(x), L.G * L.beta(x, q) * da, rtol=1e-9, atol=0)          # FreeLaw loses ~1e-11 (1 − tanh², q + c)
        e = L.eta(x, q) * np.array([1, -1, 1, -1, 1, -1])
        assert np.allclose(F.inv(e), x * np.array([1, -1, 1, -1, 1, -1]), rtol=1e-11, atol=1e-15)
        assert np.allclose(F.inv(e, x0=x[::-1]), F.inv(e), rtol=1e-12, atol=1e-15)
        assert np.allclose(F.dw(x), L.m * L.relevance(x, q), rtol=1e-13, atol=0)


def test_dual_drift_is_the_gradient_of_the_projection_dual():
    """∇F(ζ) = Φw − y for F(ζ) = Ψ*(Φᵀζ) − yᵀζ (finite differences), and dψ*/dη = w."""
    Phi, y, q, b = small_problem()
    F = DL.FastFree(L, q)
    rng = np.random.default_rng(1)
    for _ in range(3):
        z = rng.standard_normal(Phi.shape[0])
        f, g, al = DL.potential(F, Phi, y, z)
        h = 1e-6
        fd = np.array([(DL.potential(F, Phi, y, z + h * e)[0] - DL.potential(F, Phi, y, z - h * e)[0]) / (2 * h) for e in np.eye(z.size)])
        assert np.abs(fd - g).max() < 1e-7 * (1 + np.abs(g).max())
    x = np.array([0.2, 1.3]); qq = q[:2]; Fq = DL.FastFree(L, qq)
    h = 1e-6
    dpsi = (Fq.psi_star(Fq.inv(Fq.eta(x) + h)) - Fq.psi_star(Fq.inv(Fq.eta(x) - h))) / (2 * h)
    assert np.allclose(dpsi, Fq.w(x), rtol=1e-7)


def test_fd_covariance_is_the_linear_loops_stationary_covariance():
    """The stationary covariance of the linearised dual, dδζ = −ΦDΦᵀδζ dt/N + (σ/N)dW, mapped to δw = DΦᵀδζ, is T·D^{1/2}PD^{1/2}; on a
    support with D negligible elsewhere it is T·D_S (uncorrelated, whatever the design)."""
    rng = np.random.default_rng(2)
    N, d, sig = 4, 7, 0.3
    Phi = rng.standard_normal((N, d))
    D = rng.uniform(0.1, 2.0, d)
    T = sig ** 2 / (2 * N)
    A = Phi @ (D[:, None] * Phi.T) / N
    Cz = solve_continuous_lyapunov(-A, -(sig / N) ** 2 * np.eye(N))
    Cw = (D[:, None] * Phi.T) @ Cz @ (Phi * D[None, :])
    assert np.allclose(DL.fd_covariance(Phi, D, T), Cw, rtol=1e-9, atol=1e-14)
    D2 = np.r_[0.3, 1.1, 0.05, np.full(d - 3, 1e-13)]
    C2 = DL.fd_covariance(Phi, D2, T)
    assert np.allclose(C2[:3, :3], T * np.diag(D2[:3]), rtol=1e-6, atol=1e-12)


def test_ito_and_stratonovich_readings_of_q():
    """One step from a fixed state, averaged over the noise: read in Itô's sense (Euler–Maruyama on (α, log β)) Q moves by the drift
    (γβ/2)(γa² − Gβa′(1 + a²/a_max²))σ²|φ|²h/N²; read in Stratonovich's sense it is kept exactly."""
    Phi, y, q, b = small_problem(eps=1.0)
    d = Phi.shape[1]
    al0 = np.array([0.0, 0.4, -0.8, 1.5, -2.2])
    lb0 = np.log(L.beta(al0, q))
    sig, h, M = 2.0, 1e-4, 20000
    dq = np.zeros(d)
    for s in range(M):
        ts, A, B = DL.primal_sde(L, Phi, y, al0, lb0, sig, h, h, seed=s, scheme='ito')
        dq += L.Q(A[-1], np.exp(B[-1])) - q
    pred = DL.ito_q_drift(L, al0, np.exp(lb0), Phi, sig) * h
    assert np.abs(dq / M / pred - 1).max() < 0.06
    ts, A, B = DL.primal_sde(L, Phi, y, al0, lb0, sig, 1e-3, 1.0, seed=0, scheme='stratonovich')
    assert np.abs(L.Q(A, np.exp(B)) - q).max() < 1e-12


def test_mirror_update_is_an_integrator_of_the_residual():
    """Σ_k h(y − Φw_k) = N(ζ_n − ζ_0) − σ√h Σ_k ξ_k exactly: the time-averaged fit is exact up to O(1/t) and the noise's mean."""
    Phi, y, q, b = small_problem(eps=1e-2)
    F = DL.FastFree(L, q)
    N = Phi.shape[0]
    xi = np.random.default_rng(3).standard_normal((2001, N))
    h, sig = 0.05, 0.7
    ts, Z, W, al = DL.mirror_update(F, Phi, y, sig, h, 100.0, noise=lambda k: xi[k])
    lhs = h * (y[None, :] - W[:-1] @ Phi.T).sum(0)
    assert np.allclose(lhs, N * (Z[-1] - Z[0]) - sig * np.sqrt(h) * xi[:2000].sum(0), rtol=0, atol=1e-10)


def well_conditioned(N=2, d=4, eps=2.0):
    """Orthonormal rows scaled by √d (ΦΦᵀ = dI), masses εb_i with b_i ∈ U(1, 2): the dual mixes in a few time units."""
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    Phi = np.sqrt(d) * Q[:N]
    y = Phi @ np.r_[0.6, -0.4, np.zeros(d - 2)]
    b = rng.uniform(1.0, 2.0, d)
    return Phi, y, L.Q(np.zeros(d), eps * b), b


def test_gibbs_law_mirror_update_against_metropolis():
    """Proposition D‴.50(i) on a 2 × 4 dictionary at σ = 0.7 (T = 0.1225, the coefficients' law far from Gaussian): the mirror update's
    long-run means and standard deviations are those of a Metropolis sampler of π ∝ exp(−F/T), and E[Φw] = y."""
    Phi, y, q, b = well_conditioned()
    N, d = Phi.shape
    sig = 0.7
    T = sig ** 2 / (2 * N)
    F = DL.FastFree(L, q)
    ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.05, 1000.0, rec_every=10 ** 9)
    Dm = F.dw(al)
    cov = (2.38 ** 2 / N) * T * np.linalg.inv(Phi @ (Dm[:, None] * Phi.T))
    Zs, Wr, acc = DL.gibbs_rwm(F, Phi, y, T, Z[-1], cov, 100_000, seed=1, thin=5)
    Wr = Wr[len(Wr) // 10:]
    ts, Zm, Wm, al = DL.mirror_update(F, Phi, y, sig, 0.01, 1500.0, seed=2, rec_every=5)
    Wm = Wm[len(Wm) // 10:]
    sd = Wr.std(0)
    assert np.abs(Wm.mean(0) - Wr.mean(0)).max() < 0.25 * sd.min()
    assert np.abs(Wm.std(0) / sd - 1).max() < 0.1
    assert np.abs(Wm.mean(0) - W[-1]).max() > 0.2 * sd.min()          # the noise moves the mean away from the projection …
    assert np.abs(Wm.mean(0) @ Phi.T - y).max() < 0.03                 # … but not the mean fit


def test_heun_follows_the_mirror_update_pathwise():
    """The paradigm's law on its curves read in Stratonovich's sense (stochastic Heun in α) and the mirror update, driven by the same noise,
    stay together on a finite horizon, at a distance of order h (strong order one)."""
    Phi, y, q, b = well_conditioned()
    F = DL.FastFree(L, q)
    err = {}
    for h in (0.01, 0.0025):
        ts, Z, Wm, al = DL.mirror_update(F, Phi, y, 0.7, h, 20.0, seed=5, rec_every=int(round(1 / h)))
        ts2, A, B = DL.primal_sde(L, Phi, y, np.zeros(4), np.log(2.0 * b), 0.7, h, 20.0, seed=5, rec_every=int(round(1 / h)))
        err[h] = np.abs(F.eta(np.abs(A)) * np.sign(A) - Z @ Phi).max()
    assert err[0.0025] < 5e-3 and err[0.01] / err[0.0025] > 2.5


def test_gibbs_identity_pins_the_temperature():
    """The exact identity E_π[∇F∇Fᵀ] = T·E_π[∇²F] — Cov_π(Φw) = T·ΦE_π[D]Φᵀ at every noise level, linear in T — holds on the Metropolis
    sampler's draws to a few per cent and rejects the same draws read at 1.15T: the sampler's temperature is pinned, not only its law's
    shape (Proposition D‴.50(i), round 40)."""
    Phi, y, q, b = well_conditioned()
    N, d = Phi.shape
    T = 0.7 ** 2 / (2 * N)
    F = DL.FastFree(L, q)
    ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.05, 1000.0, rec_every=10 ** 9)
    Dm = F.dw(al)
    cov = (2.38 ** 2 / N) * T * np.linalg.inv(Phi @ (Dm[:, None] * Phi.T))
    Zs, Wr, acc = DL.gibbs_rwm(F, Phi, y, T, Z[-1], cov, 100_000, seed=11, thin=5)
    Wr = Wr[len(Wr) // 10:]
    C = np.cov((Wr @ Phi.T).T)
    ED = np.array([F.dw(F.inv(e)) for e in (Zs[len(Zs) // 10:] @ Phi)]).mean(0)
    pred = T * Phi @ (ED[:, None] * Phi.T)
    rel = lambda P: np.linalg.norm(C - P) / np.linalg.norm(P)
    assert rel(pred) < 0.04
    assert rel(1.15 * pred) > 0.08


def test_linear_response_is_the_projected_noise_on_an_entries_only_path():
    """On an inverse scale space path of entries only the transit's first-order response is δν(τ) = κP_{T(τ)}W(τ) at every grid point —
    the accumulated deviation projected at each entry onto the nested tangent spaces (Proposition D‴.50(iv))."""
    rng = np.random.default_rng(4)
    N, d = 3, 6
    Phi = rng.standard_normal((N, d))
    path = MF.iss_path(Phi, Phi @ np.r_[0.8, -0.5, 0.0, 0.0, 0.0, 0.0], np.ones(d))
    assert all(kind == 'enter' for t, i, kind in path['events']) and len(path['events']) >= 2
    taus = np.linspace(0.0, 1.5 * path['taus'][-1], 3001)
    xi = rng.standard_normal((taus.size - 1, N)) * np.sqrt(np.diff(taus))[:, None]
    kap = 0.37
    dnu = DL.linear_response(path, Phi, taus, xi, kap, return_nu=True)
    P = DL.tangent_projectors(path, Phi)
    Wb = np.vstack([np.zeros(N), np.cumsum(xi, 0)])
    seg = np.clip(np.searchsorted(path['taus'], taus, side='right') - 1, 0, len(P) - 1)
    closed = np.array([kap * P[seg[i]] @ Wb[i] for i in range(taus.size)])
    assert np.abs(dnu[1:] - closed[1:]).max() < 1e-12 * np.abs(closed).max()
    assert np.allclose(DL.linear_response(path, Phi, taus, xi, kap), dnu @ Phi, rtol=0, atol=1e-14)
    exits = dict(path, events=list(path['events']) + [(path['taus'][-1], path['events'][0][1], 'exit')])
    try:
        DL.linear_response(exits, Phi, taus, xi, kap)
        raise AssertionError('a path with exits must be refused')
    except ValueError:
        pass
