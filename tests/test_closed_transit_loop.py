"""Proposition D‴.49: the closed loop in C (`closed_loop`) against a tight Python integration of the same equations, Proposition D‴.46's
leak in it, and the exact averaged drive along the frozen periodic orbit (`rich_transit.FrozenOrbit`)."""
import numpy as np
from scipy.integrate import solve_ivp

from src.core import closed_loop as CL
from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
SIG, KE, OM = 0.12, 10.0, 0.2
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)


def fstar(x):
    return np.exp(-(np.asarray(x, float)[..., None] - CSTAR) ** 2 / (2 * SIG ** 2)) @ ASTAR


def python_loop(eps, T, ts, metric='paradigm', lam=0.0, noise=None, regx=False, a2=0.0, events=None, dense=False):
    """The same loop on (x, α, log β, c) by LSODA at tight tolerances (the leak through mirror_free's η): the measured state x̃ = x + n(t)
    in the error, the control and (regx False) the regressor, x_m = sin ωt + a₂ sin √2ωt."""
    d = C0.size
    q = L.Q(np.zeros(d), np.full(d, eps))
    om2 = np.sqrt(2) * OM
    nz = noise or dict(w=np.zeros(0), ph=np.zeros(0), amp=np.zeros(0))

    def f(t, z):
        x = z[0]
        al, lb, c = z[1:1 + d], z[1 + d:1 + 2 * d], z[1 + 2 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        xm = np.sin(OM * t) + a2 * np.sin(om2 * t)
        dxm = OM * np.cos(OM * t) + a2 * om2 * np.cos(om2 * t)
        xt = x + nz['amp'] @ np.cos(nz['w'] * t + nz['ph'])
        xr = xm if regx else xt
        phc = np.exp(-(xt - c) ** 2 / (2 * SIG ** 2))
        phr = np.exp(-(xr - c) ** 2 / (2 * SIG ** 2))
        e = xt - xm
        u = e * phr - (lam * L.eta(al, q) if lam else 0.0)
        v = L.G * b * a if metric == 'paradigm' else (L.G * a if metric == 'wfr' else 0.0 * a)
        return np.r_[fstar(x) - phc @ (b * a / L.m) + dxm - KE * e, L.G * b * da * u, L.gamma * a * u,
                     v * e * phr * (xr - c) / SIG ** 2]
    sol = solve_ivp(f, (0, T), np.r_[0.0, np.zeros(d), np.log(np.full(d, eps)), C0], method='LSODA', rtol=1e-10, atol=1e-13,
                    t_eval=ts, max_step=1.0 if noise is None else 0.05, dense_output=dense, events=events)
    if ts is None:
        return sol
    W = np.exp(sol.y[1 + d:1 + 2 * d].T) * L.amp(sol.y[1:1 + d].T)[0] / L.m
    return sol.y[0], W, sol.y[1 + 2 * d:].T


def test_c_loop_is_the_python_loop():
    """At ε = 10⁻² through the kernels' entries (t ≤ 1.2·10⁴): the state to 10⁻⁷, under the paradigm's metric, the Wasserstein–Fisher–Rao
    metric and with fixed centres — the centres to 10⁻⁷, under the Wasserstein–Fisher–Rao metric (the inactive kernels' centres move at
    O(1) there) to 10⁻⁶ (measured 1.7·10⁻⁷)."""
    T = 1.2e4
    ts = np.linspace(0, T, 41)
    for metric in ('paradigm', 'wfr', 'fixed'):
        x, W, cc = python_loop(1e-2, T, ts, metric)
        r = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric=metric)
        assert np.abs(x - r['x']).max() < 1e-7 and np.abs(W - r['w']).max() < 1e-7
        assert np.abs(cc - r['c']).max() < (1e-6 if metric == 'wfr' else 1e-7)
        assert np.abs(W[-1]).max() > 0.05                                 # kernels have formed


def test_c_loop_leak_is_the_mirror_coordinate_leak():
    """Proposition D‴.46's leak −λη(α) on the drive (η on the curve from α = 0 with the mass ε, closed form) against mirror_free's η."""
    T = 6000.0
    ts = np.linspace(0, T, 31)
    x, W, cc = python_loop(1e-2, T, ts, 'fixed', lam=0.05)
    r = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric='fixed', lam=0.05)
    assert np.abs(W - r['w']).max() < 1e-7 and np.abs(x - r['x']).max() < 1e-7
    r0 = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric='fixed')
    assert np.abs(r['w'][-1]).sum() < np.abs(r0['w'][-1]).sum()        # the leak shrinks the coefficients


def test_frozen_orbit_is_the_periodic_solution():
    """FrozenOrbit's e_p solves ė = −k_ee + r(x_m + e) periodically (integrated over a period from e_p(0) it returns to e_p and follows
    the samples), and its truncation e_p ≈ L[r(x_m)] (OrbitBumps) holds to O(|r|²): halving r halves e_p up to a quadratic remainder."""
    T = 2 * np.pi / OM
    M = 512
    c = np.array([-0.5, 0.1, 0.6]); w = np.array([0.5, -0.3, 0.4])
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), T, M, SIG, KE, fstar)
    e = O.orbit(c, w)
    ts = T * np.arange(M) / M

    def rhs(t, y):
        x = np.sin(OM * t) + y[0]
        return [-KE * y[0] + fstar(np.array([x]))[0] - np.exp(-(x - c) ** 2 / (2 * SIG ** 2)) @ w]
    sol = solve_ivp(rhs, (0, T), [e[0]], rtol=1e-11, atol=1e-13, t_eval=np.r_[ts, T], method='DOP853')
    assert np.abs(sol.y[0][:-1] - e).max() < 1e-5 and abs(sol.y[0][-1] - e[0]) < 1e-5
    K = RT.OrbitBumps(lambda t: np.sin(OM * t), T, M, SIG, KE)
    lin = K.filt(fstar(K.xs) - K.phi(c) @ w)
    gap = np.abs(e - lin).max()
    assert 0 < gap < 0.2 * np.abs(e).max()                               # the truncation's error, not small at these parameters
    g, h, _ = O.drive(c, w)
    xp = np.sin(OM * ts) + e
    ph = np.exp(-(xp[:, None] - c[None, :]) ** 2 / (2 * SIG ** 2))
    assert np.allclose(g, ph.T @ e / M, atol=1e-12) and np.allclose(h, (ph * (xp[:, None] - c[None, :]) / SIG ** 2).T @ e / M, atol=1e-12)


def test_exact_averaged_loop_tracks_the_real_loop():
    """At ε = 10⁻⁴: the averaged loop along the frozen orbit forms the real loop's five kernels, entering within 0.5% of its times (the
    leading-order averaged loop forms six, D‴.49 A)."""
    eps = 1e-4
    d = C0.size
    Lam = RT.transit_scale(L, eps)
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), 2 * np.pi / OM, 512, SIG, KE, fstar)
    taus = np.linspace(0, 350, 1401)
    sol = solve_ivp(RT.frozen_orbit_rhs(L, O, np.full(d, eps)), (0, Lam * 350), np.r_[np.zeros(d), C0], method='LSODA', rtol=1e-9,
                    atol=1e-12, t_eval=Lam * taus)
    al = sol.y[:d].T
    W = (eps + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    th = 1e-2 * 0.9
    ent_avg = {i: taus[np.flatnonzero(np.abs(W[:, i]) > th)[0]] for i in range(d) if (np.abs(W[:, i]) > th).any()}
    r = CL.simulate(L, eps, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), np.linspace(0, Lam * 350, 11), th=th)
    ent_real = {i: t / Lam for i, t in enumerate(r['entries']) if np.isfinite(t)}
    assert set(ent_avg) == set(ent_real)
    assert len(ent_avg) == 5 and all(abs(ent_avg[i] / ent_real[i] - 1) < 0.005 for i in ent_avg)


def test_c_loop_with_noise_desired_regressor_and_quasi_periodic_reference_is_the_python_loop():
    """At ε = 0.1 to t = 2500 (the kernels forming): measurement noise (16 tones of rms 0.01) with the regressor at x̃ or at x_m, and the
    quasi-periodic reference sin ωt + 0.4 sin √2ωt — the state to 10⁻⁷ (the coefficients) and 10⁻⁶ (the centres)."""
    T, eps = 2500.0, 0.1
    ts = np.linspace(0, T, 21)
    nz = CL.tones(0.01, 10.0, 30.0, seed=3)
    for kw in (dict(noise=nz), dict(noise=nz, regx=True), dict(a2=0.4)):
        ckw = dict(kw)
        a2 = ckw.pop('a2', 0.0)
        x, W, cc = python_loop(eps, T, ts, **kw)
        r = CL.simulate(L, eps, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM, a2, np.sqrt(2) * OM), ts, **ckw)
        assert np.abs(x - r['x']).max() < 1e-7 and np.abs(W - r['w']).max() < 1e-7 and np.abs(cc - r['c']).max() < 1e-6
        assert np.abs(W[-1]).max() > 0.01                                 # kernels are forming


def test_c_loop_error_bookkeeping_and_entry_times():
    """Against a dense Python solution (ε = 0.1, t ≤ 2500, the noise of rms 0.01): the measured error's max over the last period (the C
    loop's accepted steps: to 1%, the max over three periods being 5% larger), the true error's max and rms there (10⁻³), the state's range
    (5·10⁻⁴), and the entry times interpolated within the step (against the events |w_j| = th, to 10⁻³)."""
    T, eps, th = 2500.0, 0.1, 0.009
    per = 2 * np.pi / OM
    nz = CL.tones(0.01, 10.0, 30.0, seed=3)
    d = C0.size

    def mk(j):
        def ev(t, z):
            return abs(np.exp(z[1 + d + j]) * L.A * np.tanh(z[1 + j] / L.A) / L.m) - th
        ev.direction = 1
        return ev
    sol = python_loop(eps, T, None, noise=nz, events=[mk(j) for j in range(d)], dense=True)
    r = CL.simulate(L, eps, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), np.linspace(0, T, 11), noise=nz, th=th)
    noise = lambda t: nz['amp'] @ np.cos(nz['w'][:, None] * t[None, :] + nz['ph'][:, None])
    tt = np.linspace(T - per, T, 200001)
    X = sol.sol(tt)[0]
    em = np.abs(X + noise(tt) - np.sin(OM * tt))
    assert abs(r['e_last'] / em.max() - 1) < 1e-2
    tt3 = np.linspace(T - 3 * per, T, 300001)
    assert np.abs(sol.sol(tt3)[0] + noise(tt3) - np.sin(OM * tt3)).max() > 1.03 * em.max()
    assert abs(r['e_true_last'] / np.abs(X - np.sin(OM * tt)).max() - 1) < 1e-3
    assert abs(r['e_true_rms_last'] / np.sqrt(np.mean((X - np.sin(OM * tt)) ** 2)) - 1) < 1e-3
    XA = sol.sol(np.linspace(0, T, 1000001))[0]
    assert abs(r['x_range'][0] - XA.min()) < 5e-4 and abs(r['x_range'][1] - XA.max()) < 5e-4
    n_ent = 0
    for j in range(d):
        if sol.t_events[j].size:
            assert abs(r['entries'][j] - sol.t_events[j][0]) < 1e-3
            n_ent += 1
    assert n_ent >= 5


def test_frozen_orbit_rhs_metrics():
    """`frozen_orbit_rhs`: the coefficients' rate is the curve mobility times g under every metric; the centres' is Gβa·h (paradigm), Ga·h
    (Wasserstein–Fisher–Rao) and 0 (fixed)."""
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), 2 * np.pi / OM, 256, SIG, KE, fstar)
    d = C0.size
    rng = np.random.default_rng(4)
    al, c = 0.6 * rng.standard_normal(d), C0 + 0.02 * rng.standard_normal(d)
    b0 = np.full(d, 1e-3)
    b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
    a = L.A * np.tanh(al / L.A)
    g, h, _ = O.drive(c, b * a / L.m)
    for metric, speed in (('paradigm', L.G * b * a), ('wfr', L.G * a), ('fixed', 0.0 * a)):
        out = RT.frozen_orbit_rhs(L, O, b0, metric)(0.0, np.r_[al, c])
        assert np.allclose(out[:d], RT.curve_mobility(L, al, b0) * g, rtol=1e-12, atol=1e-15)
        assert np.allclose(out[d:], speed * h, rtol=1e-12, atol=1e-15)


def test_tones_have_the_given_rms():
    """`tones(s, …)`: the time-mean square of the sum of tones is s² (equal amplitudes s(2/k)^{1/2})."""
    nz = CL.tones(0.02, 10.0, 30.0, seed=7)
    t = np.linspace(0, 2e4, 2000001)
    n = np.zeros_like(t)
    for w, ph, amp in zip(nz['w'], nz['ph'], nz['amp']):
        n += amp * np.cos(w * t + ph)
    assert abs(np.sqrt(np.mean(n ** 2)) / 0.02 - 1) < 1e-2
