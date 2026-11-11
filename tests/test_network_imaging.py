"""Tests of experiments/exp_network_imaging.py: the declared frames, the noise laws, the network's right-hand side against finite differences of the loss, its conservation law and
its sparse Jacobian, the tau <-> t mapping, a one-emitter frame on which the network and the path agree, the signed path, the comparison and its criteria on planted cases, the
stop and snapshot rules, the bridge machinery, the raw file and the report. Fast: the networks are on 8 x 8 .. 10 x 10 frames."""
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_network_imaging as X      # noqa: E402
from src.core import imaging as IM   # noqa: E402


# ---- the declared frames ------------------------------------------------------------------------------------------------------------------------------------------
def test_frames_are_the_declared_ones():
    assert X.FRAMES == tuple(['S%02d' % k for k in range(1, 13)] + ['S03n', 'S06n', 'S10n', 'S12n'] + ['A%02d' % k for k in range(1, 7)] + ['A03n', 'A06n'])
    ks = [X.SPECS[n]['K'] for n in X.FRAMES]
    assert ks == [3, 4, 5, 6, 3, 4, 5, 6, 3, 5, 4, 6, 5, 4, 5, 6, 6, 7, 8, 9, 10, 8, 8, 8]
    assert [X.SPECS['S%02d' % k]['level'] for k in range(1, 13)] == list('LLLLMMMMHHHH')
    assert [n for n in X.FRAMES if X.SPECS[n]['pair']] == ['S06', 'S10', 'S06n', 'S10n'] and X.SPECS['S06']['pair'] == 2.5 and X.SPECS['S10']['pair'] == 2.0
    assert [n for n in X.FRAMES if X.SPECS[n]['sat']] == ['A03', 'A06', 'A03n', 'A06n']
    assert [n for n in X.FRAMES if not X.SPECS[n]['noisy']] == ['S03n', 'S06n', 'S10n', 'S12n', 'A03n', 'A06n']
    assert len(X.FRAMES) == 24 and len([n for n in X.FRAMES if X.SPECS[n]['noisy']]) == 18
    for n in X.FRAMES:
        s = X.SPECS[n]
        assert s['source'] == (n if s['noisy'] else n[:-1]) and s['seed'] == (1000 if n[0] == 'S' else 2000) + int(s['source'][1:])


def test_geometry_obeys_the_design_and_is_deterministic():
    for n in X.FRAMES:
        s = X.SPECS[n]
        c, f = X.draw_geometry(s)
        c2, f2 = X.draw_geometry(s)
        assert np.array_equal(c, c2) and np.array_equal(f, f2) and c.shape == (s['K'], 2) and f.shape == (s['K'],)
        hi = 17.0 if s['kind'] == 'smlm' else 21.0
        assert (c >= 3.0).all() and (c <= hi).all() and (f > 0).all()
        D = np.linalg.norm(c[:, None] - c[None], axis=2) + 1e9 * np.eye(s['K'])
        special = s['pair'] or s['sat']
        if special:
            sep = (s['pair'] or X.SAT_SIGMAS) * X.SIGMA
            assert abs(np.linalg.norm(c[0] - c[1]) - sep) < 1e-12
            D[0, 1] = D[1, 0] = 1e9
        assert D.min() >= X.MIN_SEP - 1e-12
        if s['kind'] == 'smlm':
            lo_, hi_ = X.FLUX[s['level']]
            assert (f >= lo_).all() and (f <= hi_).all()
        else:
            assert (f >= 400.0 - 1e-9).all() and (f <= 25000.0 + 1e-9).all() or s['sat']
            if s['sat']:
                assert 15000.0 <= f[0] <= 25000.0 and 0.03 * f[0] <= f[1] <= 0.06 * f[0] + 1e-9
    for n in X.FRAMES:                                        # a noise-free copy has the geometry of its source
        s = X.SPECS[n]
        if not s['noisy']:
            assert all(np.array_equal(a, b) for a, b in zip(X.draw_geometry(s), X.draw_geometry(X.SPECS[s['source']])))
    assert not np.array_equal(X.draw_geometry(X.SPECS['S01'])[0], X.draw_geometry(X.SPECS['S02'])[0])


def test_simulated_frames_follow_the_noise_laws():
    for n in ('S01', 'S09'):
        sim = X.simulate(X.SPECS[n])
        P0 = IM.PixelPSF(20, 20, X.SIGMA)
        sig = P0.synth(sim['c'], sim['f'])
        assert sim['frame'].shape == (20, 20) and sim['bg'] == 100.0
        assert np.allclose(sim['var'], IM.emccd_variance(100.0 + sig), rtol=1e-14, atol=0)
        assert np.array_equal(sim['y0'], sim['frame'] - 100.0)
        z = sim['y0'] - sig
        assert abs(z.mean()) < 3 * z.std() / 20 and 0.8 < (z / np.sqrt(sim['var'])).std() < 1.25                # zero mean, standardised sd near 1
        assert np.isclose(sim['var_a'], 2.0 / 0.9, rtol=1e-14, atol=0) and np.isclose(sim['var_b0'], IM.emccd_variance(0.0), rtol=1e-14, atol=0)
    for n in ('A02',):
        sim = X.simulate(X.SPECS[n])
        P0 = IM.PixelPSF(24, 24, X.SIGMA)
        mean = 300.0 + P0.synth(sim['c'], sim['f'])
        assert np.allclose(sim['var'], mean / 4.6 + 1.0, rtol=1e-14, atol=0)
        assert np.isclose(sim['var_a'], 1 / 4.6, rtol=1e-14, atol=0) and sim['var_b0'] == 1.0 and sim['bg'] == 300.0
        zs = (sim['frame'] - mean) / np.sqrt(sim['var'])
        assert abs(zs.mean()) < 0.2 and 0.9 < zs.std() < 1.1
    free = X.simulate(X.SPECS['S03n'])
    noisy = X.simulate(X.SPECS['S03'])
    P0 = IM.PixelPSF(20, 20, X.SIGMA)
    assert np.array_equal(free['y0'], P0.synth(free['c'], free['f'])) and np.array_equal(free['var'], noisy['var']) and np.array_equal(free['c'], noisy['c'])
    assert not np.array_equal(X.simulate(X.SPECS['S01'])['frame'], X.simulate(X.SPECS['S02'])['frame'])
    assert np.array_equal(X.simulate(X.SPECS['S05'])['frame'], X.simulate(X.SPECS['S05'])['frame'])


def test_astronomy_dark_noise_has_the_declared_variance(monkeypatch):
    monkeypatch.setitem(X.ASTRO, 'dark', 9.0)                                       # a dark variance of 9 DN^2 (noise sd 3 DN), where a mistaken sd = variance would give 9
    sim = X.simulate(X.SPECS['A02'])
    P0 = IM.PixelPSF(24, 24, X.SIGMA)
    mean = 300.0 + P0.synth(sim['c'], sim['f'])
    assert np.allclose(sim['var'], mean / 4.6 + 9.0, rtol=1e-14, atol=0) and sim['var_b0'] == 9.0
    zs = (sim['frame'] - mean) / np.sqrt(sim['var'])
    assert 0.9 < zs.std() < 1.1 and abs(zs.mean()) < 0.2


def test_build_gives_the_lattice_the_highpassed_data_and_the_weights():
    for hp, nm in ((X.HP, 'S02'), (None, 'A01')):
        pr = X.build(X.SPECS[nm], hp)
        n = pr['sim']['nx']
        assert pr['lat'].shape == (n + 2, n + 2) and pr['c0'].shape == ((n + 2) ** 2, 2) and pr['lat'].h == 1.0 and pr['lat'].margin == 1.0
        assert np.array_equal(pr['c0'][0], [-1.0, -1.0]) and np.array_equal(pr['c0'][1], [0.0, -1.0]) and np.array_equal(pr['c0'][n + 2], [-1.0, 0.0])
        assert np.array_equal(pr['y'], pr['P'].highpass(pr['sim']['y0'])) and np.allclose(pr['W'], 1.0 / pr['sim']['var'], rtol=1e-14, atol=0)
        assert (pr['P'].hp == hp) and pr['P'].sigma == 1.2
        if hp:
            assert abs(pr['y'].mean()) < abs(pr['sim']['y0']).mean()


# ---- the network --------------------------------------------------------------------------------------------------------------------------------------------------
def _small_net(eps=1e-3, n=8, d=5, hp=2.0, seed=3):
    rng = np.random.default_rng(seed)
    P = IM.PixelPSF(n, n, 1.2, hp=hp)
    c0 = rng.uniform(1.5, n - 2.5, (d, 2))
    y = rng.normal(0, 1.0, (n, n))
    W = rng.uniform(0.5, 2.0, (n, n))
    net = X.Network(P, y, W, c0, eps)
    z = net.init()
    z[:d] = rng.normal(0, 0.7, d)
    z[d:2 * d] = rng.uniform(0.4, 1.8, d)
    z[2 * d:] += rng.normal(0, 0.05, 2 * d)
    return net, z


def test_pieces_are_the_psf_synthesis_and_correlation():
    net, z = _small_net()
    comps, a, ap, beta, w, model, v, vx, vy = net.pieces(z)
    c = net.centres(z)
    assert np.allclose(model, net.P.synth(c, w), rtol=1e-14, atol=1e-15)
    v2, vx2, vy2 = net.P.corr(c, net.W * (net.y - net.P.synth(c, w)))
    assert np.allclose(v, v2, rtol=1e-13, atol=1e-14) and np.allclose(vx, vx2, rtol=1e-13, atol=1e-14) and np.allclose(vy, vy2, rtol=1e-13, atol=1e-14)
    al, nu, _, _ = net.split(z)
    assert np.allclose(a, np.tanh(al), rtol=1e-15, atol=0) and np.allclose(beta, nu ** 2, rtol=1e-15, atol=0) and np.allclose(w, nu ** 2 * np.tanh(al) / 12.0, rtol=1e-15, atol=0)
    assert np.allclose(ap, 1.0 - np.tanh(al) ** 2, rtol=1e-12, atol=0) and np.allclose(net.w_of(z), w, rtol=1e-15, atol=0)
    big = np.array([9.0, -9.5])                                    # sech^2 without the cancellation of 1 - tanh^2
    assert np.allclose(X.sech2(big), 4 * np.exp(-2 * np.abs(big)) / (1 + np.exp(-2 * np.abs(big))) ** 2, rtol=1e-12, atol=0)
    assert X.sech2(1e4) == X.sech2(300.0) and X.sech2(0.0) == 1.0


def test_rhs_is_minus_the_group_rate_times_the_gradient_of_the_loss():
    net, z = _small_net()
    d = net.d

    def F(z_):
        r = net.y - net.P.synth(net.centres(z_), net.w_of(z_))
        return 0.5 * float(np.sum(net.W * r ** 2))
    grad = np.zeros_like(z)
    for j in range(len(z)):
        h = 1e-5 * max(1.0, abs(z[j]))
        zp, zm = z.copy(), z.copy()
        zp[j] += h
        zm[j] -= h
        grad[j] = (F(zp) - F(zm)) / (2 * h)
    rates = np.r_[np.full(d, 12.0 * 1.0), np.full(d, 0.1 * 12.0 / 4), np.full(2 * d, 12.0 * 1.0)]       # G m, gamma m / 4, G m
    f = net.rhs(0.0, z)
    assert np.allclose(f, -rates * grad, rtol=1e-5, atol=1e-7 * np.abs(f).max())
    assert net.nrhs == 1
    # the equivalent form: alpha' = G beta a' g, beta' = gamma beta a g
    _, a, ap, beta, w, model, v, vx, vy = net.pieces(z)
    nu = z[d:2 * d]
    assert np.allclose(f[:d], beta * ap * v, rtol=1e-13, atol=0) and np.allclose(2 * nu * f[d:2 * d], 0.1 * beta * a * v, rtol=1e-13, atol=0)
    assert np.allclose(f[2 * d:3 * d], 12.0 * w * vx, rtol=1e-13, atol=0) and np.allclose(f[3 * d:], 12.0 * w * vy, rtol=1e-13, atol=0)


def test_sparse_jacobian_has_exact_active_columns_and_the_dormant_blocks():
    net, z = _small_net(eps=1e-4)
    d = net.d
    z[d:2 * d] = [0.5, 1e-2, 1.2, 1e-2, 1e-2]                  # units 0 and 2 active (beta > 1000 eps), the others dormant
    z[:d] = [0.4, 1e-3, -0.8, 2e-3, 1e-3]
    assert list(net.active(z)) == [0, 2]
    J = net.jac(0.0, z).toarray()
    assert J.shape == (4 * d, 4 * d) and net.nfd == 8
    f0 = net.rhs(0.0, z)
    for u in (0, 2):
        for k in range(4):
            j = k * d + u
            h = 1e-6 * max(abs(z[j]), 1.0)
            zp = z.copy()
            zp[j] += h
            zm = z.copy()
            zm[j] -= h
            col = (net.rhs(0.0, zp) - net.rhs(0.0, zm)) / (2 * h)
            assert np.allclose(J[:, j], col, rtol=2e-4, atol=1e-5 * np.abs(col).max()), (u, k)
    for u in (1, 3, 4):                                         # a dormant unit's own block in alpha and nu: the analytic one (drive frozen) against the full derivative
        idx = [k * d + u for k in range(4)]
        for kk in (0, 1):
            j = idx[kk]
            h = 1e-6 * max(abs(z[j]), 1e-3)
            zp = z.copy()
            zp[j] += h
            zm = z.copy()
            zm[j] -= h
            col = (net.rhs(0.0, zp) - net.rhs(0.0, zm)) / (2 * h)
            assert np.allclose(J[idx, j], col[idx], rtol=5e-2, atol=1e-3 * np.abs(col[idx]).max() + 1e-12), (u, kk)
        assert not J[:, idx[2]].any() and not J[:, idx[3]].any()           # the centre columns of a dormant unit are dropped (its rate is w times a tiny number)
    others = np.setdiff1d(np.arange(4 * d), np.r_[[k * d + u for k in range(4) for u in (0, 2)]])
    off = J[np.ix_(others, others)]
    mask = np.ones_like(off, bool)
    for u in (1, 3, 4):
        for a_ in range(4):
            for b_ in range(4):
                mask[np.flatnonzero(others == a_ * d + u)[0], np.flatnonzero(others == b_ * d + u)[0]] = False
    assert np.abs(off[mask]).max() == 0.0                        # nothing between different dormant units, nothing from the dormant columns into the active rows beyond the blocks


def test_active_set_threshold_and_dormant_blocks_on_physical_scales():
    net, z = _small_net(eps=1e-3, d=4)
    net.y = net.y * 1e4                                                               # a large drive against the unit's own coefficient, as for a dormant unit at small eps
    d = net.d
    z[d:2 * d] = np.sqrt(np.array([2e-3, 0.5, 2.0, 0.9]))                              # beta / eps = 2, 500, 2000, 900
    assert list(net.active(z)) == [2]
    z[:d] = [0.35, 0.5, -0.4, 0.3]
    J = net.jac(0.0, z).toarray()
    for u in (0, 1, 3):                                                                # dormant: the own block in alpha and nu against the full derivative, entry by entry
        idx = [k * d + u for k in range(4)]
        for kk in (0, 1):
            j = idx[kk]
            h = 1e-6 * max(abs(z[j]), 1e-2)
            zp, zm = z.copy(), z.copy()
            zp[j] += h
            zm[j] -= h
            col = (net.rhs(0.0, zp) - net.rhs(0.0, zm)) / (2 * h)
            assert np.allclose(J[idx, j], col[idx], rtol=0.1, atol=1e-7 * np.abs(col).max()), (u, kk, J[idx, j], col[idx])


def test_drift_measure_is_relative_to_the_mass_floored_at_eps():
    net, z = _small_net(eps=1e-3, d=3)
    d = net.d
    z0 = net.init()
    zz = z0.copy()
    zz[d:2 * d] *= 0.99                                                                # nu^2 below eps: the drift is relative to eps
    assert np.allclose(net.q_drift(zz), (1.0 - 0.99 ** 2), rtol=1e-9, atol=0)
    zz = z0.copy()
    zz[d:2 * d] *= 1.01
    assert np.allclose(net.q_drift(zz), (1.01 ** 2 - 1.0) / 1.01 ** 2, rtol=1e-9, atol=0)


def test_swapped_unit_factors_equal_a_full_recomputation():
    net, z = _small_net()
    d = net.d
    comps = net.P.comps(net.centres(z))
    zp = z.copy()
    zp[2 * d + 3] += 0.01
    zp[3 * d + 3] -= 0.02
    sw = net._swap_unit(comps, 3, (zp[2 * d + 3], zp[3 * d + 3]))
    ref = net.P.comps(net.centres(zp))
    for (c1, *a1), (c2, *a2) in zip(sw, ref):
        assert c1 == c2 and all(np.allclose(p, q, rtol=1e-13, atol=1e-15) for p, q in zip(a1, a2))
    assert all(np.allclose(p, q, rtol=1e-13, atol=1e-15) for p, q in zip(comps[0][1:], net._swap_unit(comps, 0, net.centres(z)[0])[0][1:]))
    other = net._swap_unit(comps, 3, (zp[2 * d + 3], zp[3 * d + 3]))
    assert all(np.array_equal(p[:, 0], q[:, 0]) for p, q in zip(other[0][1:], comps[0][1:]))                 # the other units' columns are the very same numbers
    assert not np.array_equal(other[0][1][:, 3], comps[0][1][:, 3])


def test_conservation_law_along_the_flow_and_the_drift_measure():
    n, eps = 8, 1e-3
    P = IM.PixelPSF(n, n, 1.2)
    c_true = np.array([[3.0, 4.0], [5.5, 2.5]])
    y = P.synth(c_true, np.array([30.0, 20.0]))
    W = np.full((n, n), 4.0)
    c0 = np.array([[3.0, 4.0], [5.0, 2.0], [1.0, 6.0], [6.0, 6.0]])
    net = X.Network(P, y, W, c0, eps)
    tr = X.Tracker(net.d, 1.0)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, X.lam_of(eps) * 12.0, tr)
    assert status == 'ok' and cnt['steps'] > 10 and abs(t - X.lam_of(eps) * 12.0) < 1e-6
    assert (np.abs(net.w_of(z)) > 1.0).sum() >= 1                  # a unit has formed: the drift below is not the trivial one
    assert np.isfinite(tr.tform).sum() == (np.abs(net.w_of(z)) > 1.0).sum() and (tr.sform[np.isfinite(tr.tform)] == 1.0).all()           # and the tracker saw it grow
    assert net.q_drift(z).max() < 1e-7
    z0 = net.init()
    assert net.q_drift(z0).max() == 0.0 or net.q_drift(z0).max() < 1e-15
    zz = z0.copy()
    zz[net.d:2 * net.d] *= 1.001                                     # a planted violation of Q: nu^2 off by 0.2 %
    assert np.allclose(net.q_drift(zz), (1.001 ** 2 - 1) / 1.001 ** 2, rtol=1e-9, atol=0)
    al, nu, _, _ = net.split(z)
    assert np.allclose(np.abs(nu ** 2 - 0.025 * np.cosh(2 * al) - (eps - 0.025)) / np.maximum(nu ** 2, eps), net.q_drift(z), rtol=1e-12, atol=0)


def test_lambda_and_the_tau_t_mapping():
    law = X.LAW
    assert (law.m, law.G, law.gamma, law.A) == (12.0, 1.0, 0.1, 1.0)
    for eps in (1e-6, 1e-8):
        vs = 0.1 * 1.0 ** 2 / (4 * 1.0)
        assert np.isclose(X.lam_of(eps), (np.pi * 1.0 / (4 * 1.0)) * np.sqrt(2.0 / (vs * eps)), rtol=1e-14, atol=0)
    assert np.isclose(X.lam_of(1e-8) / X.lam_of(1e-6), 10.0, rtol=1e-13, atol=0)
    # a dormant unit under a constant drive g reaches its wall at tau = 1/g: t = Lambda / g.  One unit, data = a multiple of its own kernel.
    n = 8
    P = IM.PixelPSF(n, n, 1.2)
    c0 = np.array([[3.5, 3.5]])
    phi = P.phi(c0)[:, 0].reshape(n, n)
    y = 40.0 * phi
    W = np.full((n, n), 0.5)
    g0 = float(np.sum(W * y * phi))
    for eps, tol in ((1e-6, 0.04), (1e-8, 0.01)):
        net = X.Network(P, y, W, c0, eps)
        lam = X.lam_of(eps)
        tr = X.Tracker(1, 1.0)                                                                            # formed at |w| = 1 (the amplitude is 40)
        z, t, status, _ = X.advance(net, net.init(), 0.0, 2.0 * lam / g0, tr)
        assert status == 'ok' and abs(net.w_of(z)[0] - 40.0) < 1.0
        assert abs(tr.tform[0] * g0 / lam - 1.0) < tol and tr.tform[0] * g0 / lam < 1.0 + tol                      # entered at tau = rho / g (1 + O(eps^(1/3)))
        assert tr.sform[0] == 1.0


def test_solver_fallback_budget_and_failure_are_recorded(monkeypatch):
    n, eps = 8, 1e-3
    P = IM.PixelPSF(n, n, 1.2)
    c0 = np.array([[3.0, 4.0]])
    y = 20.0 * P.phi(c0)[:, 0].reshape(n, n)
    W = np.full((n, n), 4.0)

    class Failing:                                           # a solver that fails at its first step
        def __init__(self, fun, t0, y0, t1, **kw):
            self.status, self.t, self.y, self.nfev, self.njev, self.nlu, self.message = 'running', t0, np.array(y0), 1, 0, 0, 'planted failure'

        def step(self):
            self.status = 'failed'
    net = X.Network(P, y, W, c0, eps)
    t_end = X.lam_of(eps) * 3.0
    monkeypatch.setattr(X, 'Radau', Failing)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, t_end, X.Tracker(1, 1.0))
    assert status == 'ok' and cnt['fallback'] == 1 and t == t_end and cnt['steps'] > 5 and cnt['nfev'] >= 1 and abs(net.w_of(z)[0]) > 1.0     # BDF finished the interval
    monkeypatch.setattr(X, 'BDF', Failing)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, t_end, X.Tracker(1, 1.0))
    assert status == 'failed: planted failure' and t == 0.0 and cnt['steps'] == 0 and cnt['fallback'] == 1 and np.array_equal(z, net.init())
    monkeypatch.undo()
    net = X.Network(P, y, W, c0, eps)
    net.deadline = 0.0                                       # a budget already spent: the run ends with 'timeout' at the first check (the 64th right-hand side)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, t_end, X.Tracker(1, 1.0))
    assert status == 'timeout' and t < t_end
    monkeypatch.setattr(X, 'advance', lambda *a, **k: (_ for _ in ()).throw(RuntimeError('planted')))
    spec = dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None)
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    pr = X.build(spec, None)
    path = X.run_path(pr)
    res = X.run_network(pr, path, 1e-6)
    run = dict(zip(X.RUN_FIELDS, res['run']))
    assert run['status'] == 3.0 and run['steps'] == 0 and res['s_mid'].shape == (len(X.SN_FIELDS),) and res['w_mid'].shape == (pr['c0'].shape[0],)
    assert dict(zip(X.SN_FIELDS, res['s_mid']))['n_formed'] == 0 and dict(zip(X.SN_FIELDS, res['s_mid']))['ok'] == 0.0


def test_a_path_that_stops_at_the_empty_block_is_skipped(tmp_path, monkeypatch, capsys):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    spec = dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None)
    monkeypatch.setitem(X.SPECS, 'T1', spec)
    real = X.run_path(X.build(spec, None))
    empty = dict(real, pv=real['pv'].copy())
    empty['pv'][X.PV_FIELDS.index('K_stop')] = 0
    monkeypatch.setattr(X, 'run_path', lambda prob: empty)
    raw = str(tmp_path / 'r.npz')
    X.part_network('plain', [1e-6], ['T1'], raw)
    assert 'the path stops at the empty block, skipped' in capsys.readouterr().out
    r = X.load_raw(raw)
    assert 'P|plain|T1|pv' in r and 'N|plain|1e-06|T1|run' not in r


def test_one_emitter_frame_the_network_equals_the_path(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    spec = dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None)
    for hp in (None, X.HP):
        prob = X.build(spec, hp)
        path = X.run_path(prob)
        g = dict(zip(X.PV_FIELDS, path['pv']))
        assert g['K_stop'] == 1 and g['exact_ok'] == 1.0 and g['cap_active'] == 1.0 and np.isinf(g['tau_next']) and np.isclose(g['tau_mid'], 2.5 * g['tau_K'], rtol=1e-14, atol=0)
        res = X.run_network(prob, path, 1e-6)
        run = dict(zip(X.RUN_FIELDS, res['run']))
        assert run['status'] == 0.0 and run['fallback'] == 0.0 and run['steps'] > 100
        for key in ('s_q1', 's_mid', 't_q1', 't_mid'):
            s = dict(zip(X.SN_FIELDS, res[key]))
            assert s['ok'] == 1.0 and s['n_formed'] == 1 and s['id_equal'] == 1.0 and s['order_ok'] == 1.0 and s['n_neg'] == 0
            assert s['max_d'] < 1e-5 and s['rel_model'] < 1e-6 and s['amp_max'] < 1e-6 and s['dormant'] < 1e-6 and s['tdev'] < 0.05
        i = int(path['blk_ids'][0])
        assert np.abs(res['w_mid'][i] - path['blk_w'][0]) < 1e-6 * path['blk_w'][0] and np.abs(res['c_mid'][i] - path['blk_c'][0]).max() < 1e-5
        assert np.isclose(res['tform'][i] * 1.0, path['blk_tau'][0], rtol=0.05, atol=0) and res['sform'][i] == 1.0
        assert np.isnan(res['tform']).sum() == res['tform'].size - 1 and np.isnan(run['tau_neg_net']) and run['id_neg_net'] == -1.0


def test_runner_up_margin_of_a_single_emitter_is_one_minus_the_second_best_correlation_ratio(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    spec = dict(X.SPECS['S05'], K=1, seed=78, noisy=False, name='T1', pair=None)
    for hp in (None, X.HP):
        pr = X.build(spec, hp)
        path = X.run_path(pr)
        g = dict(zip(X.PV_FIELDS, path['pv']))
        V = pr['lat'].corr(pr['W'] * pr['y'])
        order = np.argsort(V.ravel())[::-1]
        best, second = V.ravel()[order[0]], V.ravel()[order[1]]
        assert int(path['blk_ids'][0]) == order[0] and g['K_stop'] == 1
        assert np.isclose(g['margin_min'], 1.0 - second / best, rtol=0, atol=1e-8) and int(g['margin_id']) == order[1]
        assert 0 < g['margin_min'] <= g['margin_mid'] + 1e-8 and g['margin_mid'] < 1.0 + 1e-8
        assert np.isclose(g['tau_K'], 1.0 / best, rtol=1e-12, atol=0)
    # dual_margin on the path of a two-emitter frame: positive, below 1, the block's kernels excluded
    pr = X.build(dict(spec, K=2, seed=4244), None)
    path = X.run_path(pr)
    g = dict(zip(X.PV_FIELDS, path['pv']))
    assert 0 < g['margin_min'] < 1 and int(g['margin_id']) not in path['blk_ids'].astype(int).tolist() and g['margin_min'] <= g['margin_mid'] + 1e-12


@pytest.mark.skipif(not os.path.exists(X.RAW), reason='the stored raw is not there')
def test_declared_frames_reproduce_the_stored_paths():
    """The frames' seeds, geometry, noise laws and weights, the path, its stop block and snapshot times, its polish, margins and the signed references (A02 and A04 have negative
    entries of the signed path before tau_mid), recomputed for nine declared frame-variant pairs (plain A06n has the one exit of a path), against the stored raw (everything but the timing)."""
    raw = X.load_raw(X.RAW)
    checked = 0
    for var, hp, fr in (('hp', X.HP, 'S05'), ('hp', X.HP, 'S06n'), ('hp', X.HP, 'A01'), ('hp', X.HP, 'A02'), ('hp', X.HP, 'A04'), ('hp', X.HP, 'A03n'), ('plain', None, 'S09'), ('plain', None, 'A03n'), ('plain', None, 'A06n')):
        stored = X.load_path(raw, var, fr)
        if stored is None:
            continue
        path = X.run_path(X.build(X.SPECS[fr], hp))
        for k in ('taus', 'blk_c', 'blk_w', 'blk_ids', 'blk_tau', 'true_c', 'true_f', 'hist_k', 'hist_tau'):
            assert np.allclose(path[k], stored[k], rtol=1e-9, atol=1e-12, equal_nan=True), (var, fr, k)
        keep = [i for i, f in enumerate(X.PV_FIELDS) if f != 'secs']
        assert np.allclose(path['pv'][keep], stored['pv'][keep], rtol=1e-9, atol=1e-12, equal_nan=True), (var, fr)
        for nm in ('sq1', 'smid'):
            for k in X.SREF_KEYS:
                assert np.allclose(path[nm][k], stored[nm][k], rtol=1e-9, atol=1e-12, equal_nan=True), (var, fr, nm, k)
        checked += 1
    assert checked == 9


# ---- the comparison and its criteria, on planted cases -------------------------------------------------------------------------------------------------------------
def _planted(w_net_extra=None, shift=0.0, amp=1.0, dormant=0.0):
    P = IM.PixelPSF(12, 12, 1.2)
    cp = np.array([[3.0, 4.0], [8.0, 7.5]])
    wp = np.array([1000.0, 500.0])
    c0 = cp.copy()
    c0[0, 0] += shift
    c = np.vstack([c0, [[1.0, 1.0], [10.0, 10.0]]])
    w = np.r_[wp * [amp, 1.0], dormant * 1000.0, 0.0]
    if w_net_extra is not None:
        c = np.vstack([c, w_net_extra[0]])
        w = np.r_[w, w_net_extra[1]]
    return P, c, w, cp, wp


def test_compare_planted_cases():
    thr = X.formed_threshold(np.array([1000.0, 500.0]))
    assert thr == 10.0
    W = np.full((12, 12), 2.0)
    P, c, w, cp, wp = _planted()
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_formed'] == 2 and s['n_unmatched'] == 0 and s['max_d'] == 0.0 and s['rel_model'] == 0.0 and s['amp_max'] == 0.0 and s['dormant'] == 0.0 and s['n_neg'] == 0
    assert (s['c1'], s['c2'], s['c3'], s['c4'], s['ok']) == (1.0, 1.0, 1.0, 1.0, 1.0)
    P, c, w, cp, wp = _planted(shift=0.3)                               # one centre 0.3 px off: C2 fails, C1 holds (within 0.5 px)
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['max_d'], 0.3, rtol=1e-12, atol=0) and (s['c1'], s['c2']) == (1.0, 0.0) and s['ok'] == 0.0 and np.isclose(s['rms_d'], 0.3 / np.sqrt(2.0), rtol=1e-12, atol=0)
    P, c, w, cp, wp = _planted(shift=0.2)
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['max_d'], 0.2, rtol=1e-12, atol=0) and s['c2'] == 1.0 and 0 < s['rel_model'] < 1
    M = P.synth(c, w) - P.synth(cp, wp)
    assert np.isclose(s['rel_model'], np.linalg.norm(M) / np.linalg.norm(P.synth(cp, wp)), rtol=1e-12, atol=0) and np.isclose(s['nsd_model'], np.sqrt(2.0 * np.mean(M ** 2)), rtol=1e-12, atol=0)
    P, c, w, cp, wp = _planted(amp=1.002)                                # amplitude 0.2 % off: C3 fails
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['amp_max'], 0.002, rtol=1e-9, atol=0) and np.isclose(s['amp_rms'], 0.002 * 1000 / np.hypot(1000, 500), rtol=1e-9, atol=0) and (s['c1'], s['c2'], s['c3'], s['ok']) == (1.0, 1.0, 0.0, 0.0)
    P, c, w, cp, wp = _planted(dormant=2e-3)                              # a dormant unit at 0.2 % of the largest (below the formed threshold, 1 %): C4 fails
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['dormant'], 2e-3, rtol=1e-12, atol=0) and (s['c1'], s['c4'], s['ok']) == (1.0, 0.0, 0.0)
    P, c, w, cp, wp = _planted(w_net_extra=(np.array([[5.0, 5.0]]), 200.0))   # an extra formed unit far from the block: count and unmatched
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_formed'] == 3 and s['n_unmatched'] == 1 and s['c1'] == 0.0 and s['c2'] == 0.0 and s['n_neg'] == 0
    P, c, w, cp, wp = _planted(w_net_extra=(np.array([[5.0, 5.0]]), -300.0))  # the same, negative: n_neg and neg_max
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_neg'] == 1 and np.isclose(s['neg_max'], 0.3, rtol=1e-12, atol=0) and s['c1'] == 0.0
    P, c, w, cp, wp = _planted(w_net_extra=(np.array([[5.0, 5.0]]), -5.0))    # a negative unit below the formed threshold is not formed: n_neg 0, neg_max still reports it
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_formed'] == 2 and s['n_neg'] == 0 and np.isclose(s['neg_max'], 0.005, rtol=1e-12, atol=0) and s['c1'] == 1.0 and s['c3'] == 0.0 and s['dormant'] == 0.005
    P, c, w, cp, wp = _planted()
    w2 = w.copy()
    w2[0] = 0.0                                                           # a missing unit
    s = X.compare(P, W, c, w2, cp, wp, thr)
    assert s['n_formed'] == 1 and s['c1'] == 0.0 and not np.isnan(s['max_d']) and s['max_d'] == 0.0 and s['c2'] == 0.0 and s['ok'] == 0.0
    s = X.compare(P, W, c, np.zeros(len(w)), cp, wp, thr)
    assert s['n_formed'] == 0 and np.isnan(s['max_d']) and s['ok'] == 0.0 and np.isnan(s['dormant']) and s['c1'] == 0.0


def test_compare_further_planted_cases():
    thr = 10.0
    W = np.full((12, 12), 2.0)
    P, c, w, cp, wp = _planted(amp=1.5, w_net_extra=(np.array([[5.0, 5.0]]), -300.0))             # the network's largest is 1500: the ratios are to it, not to the block's 1000
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['neg_max'], 300.0 / 1500.0, rtol=1e-12, atol=0)
    P, c, w, cp, wp = _planted(amp=1.5, dormant=1e-3)
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert np.isclose(s['dormant'], 1.0 / 1500.0, rtol=1e-12, atol=0)
    P, c, w, cp, wp = _planted(w_net_extra=(cp[:1] + [0.7, 0.0], 200.0))                          # an extra formed unit 0.7 px from a block unit has no block unit within 0.5 px ...
    assert X.compare(P, W, c, w, cp, wp, thr)['n_unmatched'] == 1
    P, c, w, cp, wp = _planted(w_net_extra=(cp[:1] + [0.4, 0.0], 200.0))                          # ... and one 0.4 px away has
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_unmatched'] == 0 and s['c1'] == 0.0
    # two formed units, one far from the block: the count is right, C1 and C2 fail
    P = IM.PixelPSF(12, 12, 1.2)
    c = np.array([[3.0, 4.0], [5.0, 5.0], [1.0, 1.0]])
    w = np.array([1000.0, 500.0, 0.0])
    s = X.compare(P, W, c, w, cp, wp, thr)
    assert s['n_formed'] == 2 and s['n_unmatched'] == 1 and (s['c1'], s['c2']) == (0.0, 0.0) and s['ok'] == 0.0
    # C3 needs the amplitudes as well as the image: a faint unit 20 % off barely changes the image
    wp2 = np.array([10000.0, 5.0])
    w3 = np.r_[wp2 * [1.0, 1.2], 0.0, 0.0]
    c3 = np.vstack([cp, [[1.0, 1.0], [10.0, 10.0]]])
    s = X.compare(P, W, c3, w3, cp, wp2, X.formed_threshold(wp2))
    assert s['rel_model'] < 1e-3 and np.isclose(s['amp_max'], 0.2, rtol=1e-9, atol=0) and (s['c1'], s['c2'], s['c3'], s['ok']) == (1.0, 1.0, 0.0, 0.0)
    # C4 holds at 5e-4 of the largest and fails at 1.5e-3 (a dormant unit in the corner changes the image by less than 1e-3 of it)
    w4 = np.r_[wp, 0.0, 0.0]
    cc = np.vstack([cp, [[-1.0, -1.0], [10.0, 10.0]]])
    s = X.compare(P, W, cc, np.r_[wp, 0.5, 0.0], cp, wp, thr)
    assert s['dormant'] == 5e-4 and (s['c3'], s['c4'], s['ok']) == (1.0, 1.0, 1.0)
    s = X.compare(P, W, cc, np.r_[wp, 1.5, 0.0], cp, wp, thr)
    assert np.isclose(s['dormant'], 1.5e-3, rtol=1e-12, atol=0) and (s['c1'], s['c2'], s['c3'], s['c4'], s['ok']) == (1.0, 1.0, 1.0, 0.0, 0.0)


def _snapshot_setup():
    n = 10
    P = IM.PixelPSF(n, n, 1.2)
    c0 = np.array([[3.0, 4.0], [6.0, 6.0], [8.0, 3.0]])
    rng = np.random.default_rng(5)
    y = 100.0 * P.phi(c0[:1])[:, 0].reshape(n, n) + rng.normal(0.0, 1.0, (n, n))
    W = rng.uniform(0.5, 2.0, (n, n))
    prob = dict(P=P, W=W, y=y)
    net = X.Network(P, y, W, c0, 1e-6)
    z = net.init()
    from scipy.optimize import brentq
    al = brentq(lambda a: (1e-6 + 2 * 0.025 * np.sinh(a) ** 2) * np.tanh(a) / 12.0 - 100.0, 1.0, 12.0, xtol=1e-14)       # unit 0 on its curve: nu^2 = eps + 2 v_s sinh^2(alpha), w = 100
    z[0] = al
    z[3] = np.sqrt(1e-6 + 2 * 0.025 * np.sinh(al) ** 2)
    ref = dict(c=c0[:1].copy(), w=np.array([100.0]), ids=np.array([0.0]), tau=np.array([1.0]), thr=1.0)
    return prob, net, z, ref


def test_snapshot_stats_on_a_planted_state():
    prob, net, z, ref = _snapshot_setup()
    d = net.d
    tform = np.full(d, np.nan)
    tform[0] = 1.0
    g = dict(zip(X.SN_FIELDS, X.snapshot_stats(prob, ref, net, z, tform, 1.0)))
    M = 100.0 * prob['P'].phi(net.c0[:1])[:, 0].reshape(10, 10)
    assert np.isclose(net.w_of(z)[0], 100.0, rtol=1e-12, atol=0) and g['n_formed'] == 1 and g['ok'] == 1.0 and g['id_equal'] == 1.0 and g['order_ok'] == 1.0
    F = 0.5 * np.sum(prob['W'] * (prob['y'] - (M + prob['P'].synth(net.centres(z)[1:], net.w_of(z)[1:]))) ** 2)
    assert np.isclose(g['F_net'], F, rtol=1e-9, atol=0) and np.isclose(g['F_path'], 0.5 * np.sum(prob['W'] * (prob['y'] - M) ** 2), rtol=1e-12, atol=0) and g['F_net'] != 0.5 * np.sum((prob['y'] - M) ** 2)
    assert g['tdev'] == 0.0 and g['qd_formed'] < 1e-10 and g['qd_all'] < 1e-10
    # the same state against a reference whose lattice unit is another one: the identity fails, the comparison does not
    assert dict(zip(X.SN_FIELDS, X.snapshot_stats(prob, dict(ref, ids=np.array([2.0])), net, z, tform, 1.0)))['id_equal'] == 0.0
    # a drift of Q in a dormant unit shows in qd_all and not in qd_formed
    z2 = z.copy()
    z2[d + 2] *= 1.05
    g2 = dict(zip(X.SN_FIELDS, X.snapshot_stats(prob, ref, net, z2, tform, 1.0)))
    assert g2['qd_formed'] < 1e-10 and np.isclose(g2['qd_all'], (1.05 ** 2 - 1.0) / 1.05 ** 2, rtol=1e-6, atol=0)
    # an empty reference, a reference without a finite threshold, a state that is not finite: all nan
    for bad_ref, zz in ((dict(ref, w=np.zeros(0)), z), (dict(ref, thr=float('nan')), z), (ref, np.where(np.arange(len(z)) == 3, np.nan, z))):
        v = X.snapshot_stats(prob, bad_ref, net, zz, tform, 1.0)
        assert v.shape == (len(X.SN_FIELDS),) and np.isnan(v).all()


def test_run_network_with_a_stubbed_integrator_records_times_thresholds_and_statuses(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    spec = dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None)
    prob = X.build(spec, None)
    path = X.run_path(prob)
    sp = X.SPV_FIELDS.index('thr')
    path = dict(path, sq1=dict(path['sq1'], pv=path['sq1']['pv'].copy()), smid=dict(path['smid'], pv=path['smid']['pv'].copy(), w=path['smid']['w'] * 2.0))
    thr = X.pv_get(path, 'thr')
    path['sq1']['pv'][sp] = 0.25 * thr                                               # the smallest threshold is the signed path's at tau_q1
    path['smid']['pv'][sp] = 0.5 * thr
    calls = []
    out = {'status': 'ok'}

    def stub(net, z, t0, t1, tracker, rtol=X.RTOL, **kw):
        calls.append((t0, t1, tracker.thr, net.deadline, rtol))
        if len(calls) == 1:                                                           # two negative units formed (the earlier is unit 5), one positive
            tracker.tform[[3, 5, 7]] = [30.0, 10.0, 20.0]
            tracker.sform[[3, 5, 7]] = [-1.0, -1.0, 1.0]
        return z, t1, out['status'], dict(nfev=7, njev=3, nlu=5, steps=11, fallback=len(calls) - 1)
    monkeypatch.setattr(X, 'advance', stub)
    lam = X.lam_of(1e-6)
    t_before = time.time()
    res = X.run_network(prob, path, 1e-6, rtol=3e-9, budget=50.0)
    pv = dict(zip(X.PV_FIELDS, path['pv']))
    assert calls[0][:2] == (0.0, lam * pv['tau_q1']) and calls[1][:2] == (lam * pv['tau_q1'], lam * pv['tau_mid']) and calls[0][2] == 0.25 * thr and calls[1][2] == 0.25 * thr
    assert t_before + 40.0 < calls[0][3] < time.time() + 50.0 and calls[0][4] == 3e-9
    run = dict(zip(X.RUN_FIELDS, res['run']))
    assert run['status'] == 0.0 and (run['nfev'], run['njev'], run['nlu'], run['steps'], run['fallback']) == (14, 6, 10, 22, 1) and run['rtol'] == 3e-9 and run['lam'] == lam
    assert np.isclose(run['tau_neg_net'], 10.0 / lam, rtol=1e-14, atol=0) and run['id_neg_net'] == 5.0
    assert np.isclose(res['tform'][5], 10.0 / lam, rtol=1e-14, atol=0) and res['sform'][3] == -1.0 and res['tform'].shape == (prob['c0'].shape[0],)
    fp = X.SN_FIELDS.index('F_path')                                                   # the signed block at tau_mid was doubled in amplitude: its loss differs, the one at tau_q1 is the non-negative block's
    assert res['t_mid'][fp] != res['s_mid'][fp] and res['t_q1'][fp] == res['s_q1'][fp]
    out['status'] = 'timeout'
    calls.clear()
    res = X.run_network(prob, path, 1e-6)
    assert len(calls) == 1 and dict(zip(X.RUN_FIELDS, res['run']))['status'] == 1.0                    # a timeout of the first segment ends the run
    out['status'] = 'failed: planted'
    res = X.run_network(prob, path, 1e-6)
    assert dict(zip(X.RUN_FIELDS, res['run']))['status'] == 2.0


def test_advance_passes_its_tolerances_and_continues_from_the_last_accepted_step(monkeypatch):
    n, eps = 8, 1e-3
    P = IM.PixelPSF(n, n, 1.2)
    c0 = np.array([[3.0, 4.0], [5.0, 2.0]])
    y = 20.0 * P.phi(c0[:1])[:, 0].reshape(n, n)
    W = np.full((n, n), 4.0)
    net = X.Network(P, y, W, c0, eps)
    t_end = 0.3 * X.lam_of(eps) / float(np.sum(W * y * P.phi(c0[:1])[:, 0].reshape(n, n)))        # before the unit enters

    class Stub:
        inst = []

        def __init__(self, fun, t0, y0, t1, **kw):
            self.kw, self.status, self.t, self.y, self.t1 = kw, 'running', t0, np.array(y0), t1
            self.nfev, self.njev, self.nlu, self.message, self.n = 2, 1, 1, 'planted', 0
            Stub.inst.append(self)

        def step(self):
            self.n += 1
            if self.n == 1:
                self.t = self.t + 0.1 * (self.t1 - self.t)
                self.y = self.y.copy()
                self.y[0] += 0.3                                                          # an accepted step off the true trajectory
            else:
                self.status = 'failed'
    monkeypatch.setattr(X, 'Radau', Stub)
    tr = X.Tracker(net.d, 1.0)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, t_end, tr, rtol=2e-9, atol_p=3e-10, atol_c=4e-12)
    kw = Stub.inst[0].kw
    assert kw['rtol'] == 2e-9 and kw['jac'] == net.jac and np.array_equal(kw['atol'], np.r_[np.full(2 * net.d, 3e-10 * np.sqrt(eps)), np.full(2 * net.d, 4e-12)])
    assert status == 'ok' and cnt['fallback'] == 1 and cnt['steps'] >= 2
    monkeypatch.undo()
    t1 = 0.1 * t_end
    y1 = net.init()
    y1[0] += 0.3
    ref = X.BDF(net.rhs, t1, y1, t_end, rtol=2e-9, atol=np.r_[np.full(2 * net.d, 3e-10 * np.sqrt(eps)), np.full(2 * net.d, 4e-12)], jac=net.jac)
    while ref.status == 'running':
        ref.step()
    assert np.allclose(z, ref.y, rtol=1e-6, atol=1e-12) and abs(z[0] - net.init()[0]) > 0.25                    # it continued from the planted accepted state, not from the start


def test_signed_path_stop_time_and_exit_rule(monkeypatch):
    n = 14
    P = IM.PixelPSF(n, n, 1.2)
    lat = IM.Lattice(P, h=1.0, margin=1.0)
    c = np.array([[4.3, 5.1], [9.2, 8.6]])
    y = P.synth(c, np.array([800.0, -600.0]))
    sw = np.full((n, n), 0.05)
    full = X.signed_path(P, y, sw, lat, max_entries=3)
    t_enter = [t for t, ij, kind, s in full['events'] if kind == 'enter']
    cut = X.signed_path(P, y, sw, lat, max_entries=3, tau_stop=0.5 * (t_enter[0] + t_enter[1]))
    assert [e for e in cut['events'] if e[2] == 'enter'] == [e for e in full['events'] if e[2] == 'enter'][:1] and len(cut['hist']) == 2
    # a negative entrant whose fit is forced to zero stays (its correlation pushes it the right way); with the sign dropped from the test it would leave
    monkeypatch.setattr(X, 'signed_block_fit', lambda psf, y_, swv, c_, w_, s_, lo, hi, **kw: (c_, np.zeros_like(w_)))
    ev = X.signed_path(P, -np.abs(y) * 0 + P.synth(c[:1], np.array([-600.0])), sw, lat, max_entries=2)['events']
    assert [e[2] for e in ev][:1] == ['enter'] and ev[0][3] == -1.0 and not any(e[2] == 'exit' for e in ev)


def test_signed_block_fit_starts_from_the_signed_nonnegative_fit():
    n = 14
    P = IM.PixelPSF(n, n, 1.2)
    c = np.array([[4.0, 5.0], [9.0, 8.0]])
    y = P.synth(c, np.array([800.0, -600.0]))
    swv = np.full((n, n), 0.1)
    cs, ws = X.signed_block_fit(P, y, swv, c.copy(), np.zeros(2), np.array([1.0, -1.0]), c - 2.0, c + 2.0, max_nfev=1)       # one evaluation: the start point of the trust-region fit
    assert np.allclose(ws, [800.0, -600.0], rtol=1e-6, atol=0) and np.array_equal(cs, c)


def test_run_path_extends_a_path_whose_bic_minimum_is_not_followed_and_uses_three_parameters(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 12)
    monkeypatch.setitem(X.SMLM, 'nx', 12)
    spec = dict(X.SPECS['S01'], K=2, seed=5038, noisy=True, name='T', level='L')
    pr = X.build(spec, X.HP)
    with X.record_path_calls() as calls:
        path = X.run_path(pr)
    g = dict(zip(X.PV_FIELDS, path['pv']))
    R = IM.sliding_path(pr['P'], pr['y'], pr['sw'], pr['lat'], rho=1.0, local=False, r_ex=0.0, order='dual', max_entries=8, tol=X.PATH_TOL, max_nfev=X.PATH_NFEV, passes=1)
    assert IM.stop_bic(R['hist'], pr['P'].N, n_par=3)[1] == 2 and IM.stop_bic(R['hist'], pr['P'].N, n_par=4)[1] == 1
    assert g['K_stop'] == 2 and len([c for c in calls if c[0] == 'sliding_path']) == 1 and calls[0][2] == 2 + X.EXTRA_NOISY
    real = IM.stop_bic
    seen = []

    def at_the_end(hist, N, n_par=3):                       # the first call reports the minimum at the last block: the path must be extended by EXTRA_NOISY entries
        seen.append(len(hist))
        if len(seen) == 1:
            return len(hist) - 1, len(hist[-1]['w'])
        return real(hist, N, n_par=n_par)
    monkeypatch.setattr(IM, 'stop_bic', at_the_end)
    with X.record_path_calls() as calls:
        path2 = X.run_path(pr)
    sp = [c for c in calls if c[0] == 'sliding_path']
    assert [c[2] for c in sp] == [2 + X.EXTRA_NOISY, 2 + 2 * X.EXTRA_NOISY] and len(seen) == 2 and seen[1] > seen[0]
    assert dict(zip(X.PV_FIELDS, path2['pv']))['interior'] == 1.0


def test_formed_threshold_snapshot_times_and_classification():
    assert X.formed_threshold([1000.0, -500.0, 800.0]) == 10.0
    assert X.formed_threshold([1000.0, 5000.0]) == 50.0
    assert X.formed_threshold([1000.0, 50.0]) == 5.0
    assert X.formed_threshold([1000.0]) == 10.0
    tau_end, q1, mid, capped = X.snapshot_times(10.0, 18.0)
    assert (tau_end, q1, mid, capped) == (18.0, 12.0, 14.0, False)
    tau_end, q1, mid, capped = X.snapshot_times(10.0, 100.0)               # a gap beyond 3 tau_K: the cap 4 tau_K
    assert (tau_end, q1, mid, capped) == (40.0, 17.5, 25.0, True)
    tau_end, q1, mid, capped = X.snapshot_times(10.0, float('inf'))
    assert (tau_end, q1, mid, capped) == (40.0, 17.5, 25.0, True)
    tau_end, q1, mid, capped = X.snapshot_times(10.0, 40.0)
    assert capped is False and tau_end == 40.0
    base = dict(zip(X.SN_FIELDS, np.zeros(len(X.SN_FIELDS))))

    def cl(status=0.0, K=2, **kw):
        g = dict(base, **kw)
        return X.classify([g[f] for f in X.SN_FIELDS], status, K)
    assert cl(ok=1.0, c1=1.0, c2=1.0, c3=1.0, c4=1.0, n_formed=2) == 'ok'
    assert cl(status=2.0, ok=1.0) == 'solver'
    assert cl(n_formed=3, n_neg=1) == 'extra negative unit' and cl(n_formed=3, n_neg=0) == 'extra unit' and cl(n_formed=1) == 'missing unit'
    assert cl(n_formed=2, n_unmatched=1) == 'unmatched unit' and cl(n_formed=2, c1=1.0) == 'position'
    assert cl(n_formed=2, c1=1.0, c2=1.0) == 'model/amplitude' and cl(n_formed=2, c1=1.0, c2=1.0, c3=1.0) == 'dormant'


def test_tracker_interpolates_the_threshold_crossing_in_log_amplitude():
    tr = X.Tracker(3, 10.0)
    tr.update(0.0, 1.0, np.array([1e-6, 1.0, -1e-3]), np.array([1e-3, 5.0, -1.0]))
    assert np.isnan(tr.tform).all()
    tr.update(1.0, 3.0, np.array([1e-3, 5.0, -1.0]), np.array([1e-3, 20.0, -100.0]))
    assert np.isnan(tr.tform[0]) and tr.sform[0] == 0.0
    assert np.isclose(tr.tform[1], 1.0 + 2.0 * np.log(10 / 5.0) / np.log(20 / 5.0), rtol=1e-14, atol=0) and tr.sform[1] == 1.0
    assert np.isclose(tr.tform[2], 1.0 + 2.0 * np.log(10 / 1.0) / np.log(100 / 1.0), rtol=1e-14, atol=0) and tr.sform[2] == -1.0
    tr2 = X.Tracker(1, 10.0)
    tr2.update(0.0, 1.0, np.array([0.5]), np.array([-20.0]))                                   # the sign is the one at the crossing
    assert tr2.sform[0] == -1.0 and 0.0 < tr2.tform[0] <= 1.0
    tr.update(3.0, 4.0, np.array([1e-3, 20.0, -100.0]), np.array([50.0, 1.0, -1.0]))           # later steps do not change a recorded time; unit 0 now crosses
    assert np.isclose(tr.tform[1], 1.0 + 2.0 * np.log(2.0) / np.log(4.0)) and 3.0 < tr.tform[0] <= 4.0


def test_merge_formed_merges_colocated_same_sign_units_only():
    c = np.array([[3.0, 4.0], [3.01, 4.0], [3.0, 4.02], [8.0, 7.0], [8.0, 7.01], [9.0, 9.0], [9.0, 9.0], [1.0, 1.0]])
    w = np.array([600.0, 300.0, 100.0, 500.0, -200.0, 40.0, 40.0, 0.5])
    c2, w2, gone = X.merge_formed(c, w, 10.0)
    assert gone == 3 + 0 + 1 - 1 + 0 or gone == 3                             # {0, 1, 2} -> one, {5, 6} -> one; the opposite-sign pair (3, 4) stays; unit 7 is dormant
    assert gone == 3 and w2[0] == 1000.0 and np.all(w2[[1, 2]] == 0) and w2[5] == 80.0 and w2[6] == 0.0 and w2[3] == 500.0 and w2[4] == -200.0 and w2[7] == 0.5
    assert np.allclose(c2[0], (600 * c[0] + 300 * c[1] + 100 * c[2]) / 1000.0, rtol=1e-14, atol=0) and np.array_equal(c2[3], c[3]) and np.array_equal(c2[7], c[7])
    assert np.array_equal(c[1], [3.01, 4.0]) and w[0] == 600.0                  # the inputs are not modified
    c3, w3, gone3 = X.merge_formed(c, w, 10.0, d_dup=0.001)
    assert gone3 == 1 and w3[5] == 80.0                                         # only the units at the same position
    c4, w4, gone4 = X.merge_formed(c[:1], w[:1], 10.0)
    assert gone4 == 0 and np.array_equal(w4, w[:1])
    assert X.D_DUP == 0.05


def test_formation_against_the_path():
    tform = np.array([np.nan, 0.0, 10.2, 7.1, 5.05])                       # tau units (already divided by Lambda here: lam = 1)
    order, tdev = X.formation_vs_path(tform, 1.0, np.array([4, 3, 2]), np.array([5.0, 7.0, 10.0]))
    assert order == 1.0 and np.isclose(tdev, 0.02, rtol=1e-12, atol=0)
    order, tdev = X.formation_vs_path(tform, 1.0, np.array([4, 2, 3]), np.array([5.0, 7.0, 10.0]))
    assert order == 0.0
    order, tdev = X.formation_vs_path(tform, 1.0, np.array([4, 0]), np.array([5.0, 7.0]))
    assert order == 0.0 and np.isnan(tdev)
    order, tdev = X.formation_vs_path(tform * 100.0, 100.0, np.array([4, 3, 2]), np.array([5.0, 7.0, 10.0]))
    assert order == 1.0 and np.isclose(tdev, 0.02, rtol=1e-12, atol=0)


# ---- the paths ---------------------------------------------------------------------------------------------------------------------------------------------------
@pytest.fixture
def tiny(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 12)
    monkeypatch.setitem(X.SMLM, 'nx', 12)
    return dict(X.SPECS['S05'], K=2, seed=4243, noisy=True, name='T2', pair=None)


def test_path_stop_rules_and_snapshot_times(tiny):
    pr = X.build(tiny, X.HP)
    path = X.run_path(pr)
    g = dict(zip(X.PV_FIELDS, path['pv']))
    R = IM.sliding_path(pr['P'], pr['y'], pr['sw'], pr['lat'], rho=1.0, local=False, r_ex=0.0, order='dual', max_entries=int(g['n_entries']), tol=X.PATH_TOL, max_nfev=X.PATH_NFEV, passes=1)
    i, K = IM.stop_bic(R['hist'], pr['P'].N, n_par=3)
    assert g['i_stop'] == i and g['K_stop'] == K and len(path['blk_w']) == K and g['n_hist'] == len(R['hist'])
    assert np.isclose(g['tau_K'], R['hist'][i]['tau'], rtol=1e-12, atol=0) and g['interior'] == 1.0
    nxt = R['hist'][i + 1]['tau'] if i + 1 < len(R['hist']) else np.inf
    assert np.isclose(g['tau_next'], nxt, rtol=1e-12, atol=0)
    te, q1, mid, capped = X.snapshot_times(g['tau_K'], nxt)
    assert np.isclose(g['tau_mid'], mid, rtol=1e-12, atol=0) and np.isclose(g['tau_q1'], q1, rtol=1e-12, atol=0) and g['cap_active'] == float(capped) and g['tau_K'] < g['tau_q1'] < g['tau_mid'] <= g['tau_end']
    assert g['polished'] == 1.0 and g['box_active'] == 0.0 and g['kkt_pos'] < 1e-5 and g['kkt_amp'] < 1e-9 and g['thr'] == X.formed_threshold(path['blk_w'])
    tk = R['taus'][:R['hist'][i]['n_ent']]
    assert np.isclose(g['min_gap'], np.min(np.diff(tk) / tk[:-1]), rtol=1e-12, atol=0) if len(tk) > 1 else np.isnan(g['min_gap'])
    assert np.isclose(g['gap_next'], (nxt - g['tau_K']) / g['tau_K'], rtol=1e-12, atol=0)
    ids = path['blk_ids'].astype(int)
    nlx = pr['lat'].shape[1]
    assert all(np.allclose(pr['lat'].pos(u // nlx, u % nlx), pr['c0'][u]) for u in ids)                  # the flat index u = b nlx + a is the network's unit index
    assert np.all(path['blk_tau'] <= g['tau_K'] + 1e-12) and np.array_equal(path['true_c'], pr['sim']['c']) and path['hist_k'][i] == K
    assert np.isnan(g['exact_ok']) and g['n_exits'] == 0
    s = path['smid']
    sp = dict(zip(X.SPV_FIELDS, s['pv']))
    assert sp['K'] == len(s['w']) and sp['tau_snap'] <= g['tau_mid'] * (1 + 1e-9) and path['sq1']['pv'][X.SPV_FIELDS.index('tau_snap')] <= g['tau_q1'] * (1 + 1e-9)


def test_run_path_entry_times_of_a_reentering_unit_and_the_exactness_thresholds(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    pr = X.build(dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None), None)
    real = IM.sliding_path

    def fake(*a, **k):                                                    # the same path with one more exit and re-entry of its unit at 1.5 times its entry time
        R = real(*a, **k)
        t, ij, kind = R['events'][0]
        return dict(R, events=list(R['events']) + [(t, ij, 'exit'), (1.5 * t, ij, 'enter')])
    monkeypatch.setattr(IM, 'sliding_path', fake)
    path = X.run_path(pr)
    assert np.isclose(path['blk_tau'][0], 1.5 * path['hist_tau'][1], rtol=1e-14, atol=0) and dict(zip(X.PV_FIELDS, path['pv']))['n_exits'] == 1.0
    monkeypatch.setattr(IM, 'sliding_path', real)
    assert dict(zip(X.PV_FIELDS, X.run_path(pr)['pv']))['exact_ok'] == 1.0
    monkeypatch.setattr(X, 'polish_block', lambda P, y, sw, c, w, **kw: (c + 0.01, w, True, 0.01, 0.0))        # a block 0.01 px from the emitters is not "recovered to 1e-3 px"
    assert dict(zip(X.PV_FIELDS, X.run_path(pr)['pv']))['exact_ok'] == 0.0


def test_noise_free_path_recovers_the_emitters_and_the_signed_path_agrees(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 12)
    monkeypatch.setitem(X.SMLM, 'nx', 12)
    spec = dict(X.SPECS['S05'], K=2, seed=4244, noisy=False, name='T3', pair=None)
    pr = X.build(spec, None)
    path = X.run_path(pr)
    g = dict(zip(X.PV_FIELDS, path['pv']))
    assert g['K_stop'] == 2 and g['exact_ok'] == 1.0 and g['chi2'] < 1e-20
    D = np.linalg.norm(path['blk_c'][:, None] - pr['sim']['c'][None], axis=2)
    assert D.min(1).max() < 1e-4 and D.min(0).max() < 1e-4
    R = X.signed_path(pr['P'], pr['y'], pr['sw'], pr['lat'], max_entries=5)                            # positive data: the signed path is the non-negative one
    Rn = IM.sliding_path(pr['P'], pr['y'], pr['sw'], pr['lat'], rho=1.0, local=False, r_ex=0.0, order='dual', max_entries=5, tol=X.PATH_TOL, max_nfev=X.PATH_NFEV, passes=1)
    en = [(t, ij) for t, ij, kind, s in R['events'] if kind == 'enter']
    assert [ij for t, ij in en[:2]] == [ij for t, ij, kind in Rn['events'][:2]] and np.allclose([t for t, ij in en[:2]], [t for t, ij, kind in Rn['events'][:2]], rtol=1e-9, atol=0)
    assert path['smid']['pv'][X.SPV_FIELDS.index('same_nn')] == 1.0 and path['smid']['pv'][X.SPV_FIELDS.index('n_neg')] == 0 and np.isinf(path['smid']['pv'][X.SPV_FIELDS.index('tau_neg')])
    assert path['smid']['pv'][X.SPV_FIELDS.index('id_neg')] == -1.0


def test_signed_path_enters_a_negative_emitter_with_its_sign_and_block_fits_agree():
    n = 14
    P = IM.PixelPSF(n, n, 1.2)
    lat = IM.Lattice(P, h=1.0, margin=1.0)
    c = np.array([[4.3, 5.1], [9.2, 8.6]])
    f = np.array([800.0, -600.0])
    y = P.synth(c, f)
    sw = np.full((n, n), 0.05)
    R = X.signed_path(P, y, sw, lat, max_entries=3)
    ent = [(t, s) for t, ij, kind, s in R['events'] if kind == 'enter']
    assert [s for t, s in ent[:2]] == [1.0, -1.0] or [s for t, s in ent[:2]] == [-1.0, 1.0]
    h = [q for q in R['hist'] if len(q['w']) == 2][0]
    assert sorted(np.sign(h['w']).tolist()) == [-1.0, 1.0] and h['chi2'] < 1e-18 and np.allclose(np.sort(np.abs(h['w'])), [600.0, 800.0], rtol=1e-6)
    pos = h['c'][np.argmax(h['w'])]
    assert np.linalg.norm(pos - c[0]) < 1e-5
    # the signed block fit with all signs +1 equals `imaging.block_fit` (non-negative)
    c0 = np.array([[4.0, 5.0], [9.0, 9.0]])
    yp = P.synth(c, np.array([800.0, 500.0]))
    swv = np.full((n, n), 0.1)
    lo, hi = c0 - 2.0, c0 + 2.0
    cs, ws = X.signed_block_fit(P, yp, swv, c0.copy(), np.zeros(2), np.ones(2), lo, hi)
    cb, wb, _ = IM.block_fit(P, yp, swv, c0.copy(), np.zeros(2), lo, hi, joint=True, max_nfev=2000, tol=1e-12)
    assert np.allclose(cs, cb, rtol=0, atol=1e-6) and np.allclose(ws, wb, rtol=1e-6, atol=0)
    cn, wn = X.signed_block_fit(P, P.synth(c, np.array([800.0, -600.0])), swv, c0.copy(), np.zeros(2), np.array([1.0, -1.0]), lo, hi)
    assert wn[0] > 0 > wn[1] and np.allclose(wn, [800.0, -600.0], rtol=1e-5, atol=0)
    cx, wx = X.signed_block_fit(P, yp, swv, c0.copy(), np.zeros(2), np.array([1.0, -1.0]), lo, hi)           # the sign constraint binds: the negative unit sits at 0
    assert abs(wx[1]) < 1e-9 and wx[0] > 0


def test_polish_block_moves_to_the_stationary_point_and_refuses_a_large_jump():
    n = 14
    P = IM.PixelPSF(n, n, 1.2, hp=2.0)
    c = np.array([[4.3, 5.1], [9.2, 8.6]])
    f = np.array([600.0, 800.0])
    y = P.highpass(IM.PixelPSF(n, n, 1.2).synth(c, f))
    sw = np.full((n, n), 0.1)
    c1, w1, ok, dc, dw = X.polish_block(P, y, sw, c + 1e-4, f * 1.0001)
    assert ok and np.abs(c1 - c).max() < 1e-8 and np.abs(w1 / f - 1).max() < 1e-8 and 0 < dc < 2e-4 and np.isclose(dw, 1e-4, rtol=1e-3, atol=0)      # dw: relative to the largest amplitude
    c2, w2, ok2, dc2, dw2 = X.polish_block(P, y, sw, c + 0.3, f, max_shift=0.05)                 # the minimum is 0.3 px away: refused
    assert not ok2 and np.array_equal(c2, c + 0.3) and np.array_equal(w2, f) and dc2 > 0.05
    c3, w3, ok3, dc3, dw3 = X.polish_block(P, y, sw, c + 0.3, f * [1, -1], max_shift=2.0)         # an amplitude that would have to change sign: refused
    assert not ok3 or (w3 * (f * [1, -1]) > 0).all()
    ka, kp = X.block_stationarity(P, 1.0 / (sw ** 2) * 0 + sw ** 2, y, c, f)
    assert ka < 1e-9 and kp < 1e-6
    ka2, kp2 = X.block_stationarity(P, sw ** 2, y, c + [0.1, 0.0], f)
    assert ka2 > 1e-3 and kp2 > 1e-2
    assert X.block_stationarity(P, sw ** 2, y, np.zeros((0, 2)), np.zeros(0)) == (0.0, 0.0)


def test_polish_block_rejects_a_worse_point_and_runs_three_rounds(monkeypatch):
    n = 14
    P = IM.PixelPSF(n, n, 1.2, hp=2.0)
    c = np.array([[4.3, 5.1], [9.2, 8.6]])
    f = np.array([800.0, 600.0])
    y = P.highpass(IM.PixelPSF(n, n, 1.2).synth(c, f))
    sw = np.full((n, n), 0.1)
    calls = []

    def worse(res, z, **kw):                                     # a "solver" that returns a point of higher cost, within the shift bound
        calls.append(kw['method'])
        return type('R', (), {'x': z + np.r_[np.full(2, 1.0), np.full(4, 0.001)]})()
    monkeypatch.setattr(X, 'least_squares', worse)
    c1, w1, ok, dc, dw = X.polish_block(P, y, sw, c, f)
    assert calls == ['lm', 'lm', 'lm'] and not ok and np.array_equal(c1, c) and np.array_equal(w1, f) and np.isclose(dc, 0.003, rtol=1e-9, atol=0) and dw > 0


# ---- the bridge ---------------------------------------------------------------------------------------------------------------------------------------------------
def test_bridge_replica_equals_the_packaged_algorithm_and_off4_is_definition_1(tiny):
    spec = dict(tiny, K=2, seed=4246, noisy=True)
    pr = X.build(spec, X.HP)
    path = X.run_path(pr)
    K = int(path['pv'][X.PV_FIELDS.index('K_stop')])
    res = X.bridge_frame(pr, path)
    assert res['replica'][0] == 0.0 and res['secs'][0] > 0 and set(res) == {n for n, _ in X.BR_CONFIGS} | {'replica', 'secs'}
    for n, _ in X.BR_CONFIGS:
        assert res[n].shape == (2, len(X.BR_FIELDS)) and res[n][0][X.BR_FIELDS.index('K_ref')] == K
    for cfg_name in ('full', 'no_reweight'):
        cfg = dict(X.BR_CONFIGS)[cfg_name]
        o = X.packaged_fit(pr['sim'], pr['P'], pr['lat'], cfg, K)
        c, w = X.packaged_reference(pr['sim'], pr['P'], pr['lat'], cfg)
        assert np.array_equal(c, o['bic'][0]) and np.array_equal(w, o['bic'][1])
        assert o['cnt'][1].shape == (K,) or len(o['cnt'][1]) == 0
    # the cumulative chain's end is Definition 1 in the packaged pipeline: noise-free frame, same units to the fit tolerance
    free = dict(tiny, noisy=False, seed=4247)
    prf = X.build(free, X.HP)
    pathf = X.run_path(prf)
    o = X.packaged_fit(prf['sim'], prf['P'], prf['lat'], dict(X.BR_CONFIGS)['off4'], 2)
    b = X.bridge_compare(o['cnt'][0], o['cnt'][1], pathf['blk_c'], pathf['blk_w'])
    assert b[X.BR_FIELDS.index('equal')] == 1.0 and b[X.BR_FIELDS.index('max_d')] < 1e-5 and b[X.BR_FIELDS.index('amp_max')] < 1e-5


def test_bridge_compare_on_planted_blocks_and_replica_error():
    cr = np.array([[3.0, 4.0], [8.0, 7.5]])
    wr = np.array([1000.0, 500.0])
    f = lambda c, w: dict(zip(X.BR_FIELDS, X.bridge_compare(c, w, cr, wr)))
    g = f(cr, wr)
    assert g == dict(n=2, K_ref=2, same_count=1.0, max_d=0.0, amp_max=0.0, equal=1.0)
    g = f(cr[::-1] + [0.1, 0.0], wr[::-1] * [1.04, 1.0])                                         # order is irrelevant, 0.1 px, 4 % amplitude: equal
    assert np.isclose(g['max_d'], 0.1, rtol=1e-12, atol=0) and np.isclose(g['amp_max'], 0.04, rtol=1e-12, atol=0) and g['equal'] == 1.0
    assert f(cr + [0.3, 0.0], wr)['equal'] == 0.0 and f(cr, wr * [1.06, 1.0])['equal'] == 0.0                     # 0.3 px; 6 %
    g = f(np.hstack([cr, [[1.2], [1.1]]]), wr)                                                      # (x, y, sigma) rows: the width column is ignored
    assert g['equal'] == 1.0
    g = f(cr[:1], wr[:1])
    assert g['n'] == 1 and g['same_count'] == 0.0 and g['equal'] == 0.0
    g = f(np.zeros((0, 2)), np.zeros(0))
    assert g['n'] == 0 and g['same_count'] == 0.0 and g['equal'] == 0.0 and np.isnan(g['max_d'])
    g = dict(zip(X.BR_FIELDS, X.bridge_compare(np.zeros((0, 2)), np.zeros(0), np.zeros((0, 2)), np.zeros(0))))
    assert g['same_count'] == 1.0 and g['equal'] == 1.0
    assert X.replica_error((cr, wr), (cr, wr)) == 0.0 and X.replica_error((cr, wr), (cr[:1], wr[:1])) == 1.0 and X.replica_error((np.zeros((0, 2)), np.zeros(0)), (np.zeros((0, 2)), np.zeros(0))) == 0.0
    assert np.isclose(X.replica_error((cr, wr), (cr + 1e-3, wr + 0.5)), 0.5, rtol=1e-12, atol=0)


# ---- the stalled-path check ---------------------------------------------------------------------------------------------------------------------------------------
def test_stall_stats_on_planted_event_lists():
    normal = [(1.0, (1, 2), 'enter'), (2.0, (3, 4), 'enter'), (2.5, (5, 6), 'enter')]
    assert X.stall_stats(normal, 10, 3) == (0.0, 1.0, 0.0, 0.3)
    swap = normal + [(2.5, (7, 8), 'exit')]                                   # one unit enters and another leaves at one tau: two different points
    assert X.stall_stats(swap, 5, 3) == (1.0, 1.0, 0.0, 0.6)
    once = [(1.0, (1, 2), 'enter'), (1.0, (1, 2), 'exit')]                    # a point that enters and leaves at once: 2 events, not a cycle
    assert X.stall_stats(once, 40, 1) == (1.0, 2.0, 0.0, 1 / 40)
    cycle = [(1.0, (1, 2), 'enter')] + [(3.0, (9, 9), k) for k in ('enter', 'exit', 'enter', 'exit', 'enter')]
    s = X.stall_stats(cycle, 6, 4)
    assert s[0] == 2.0 and s[1] == 5.0 and s[2] == 1.0 and np.isclose(s[3], 4 / 6)
    three = [(3.0, (9, 9), 'enter'), (3.0 * (1 + 1e-14), (9, 9), 'exit'), (3.0, (9, 9), 'enter')]       # identical to 12 digits: three events at one tau
    assert X.stall_stats(three, 40, 2)[1:3] == (3.0, 1.0)
    apart = [(3.0, (9, 9), 'enter'), (3.1, (9, 9), 'exit'), (3.2, (9, 9), 'enter')]                      # a re-entry at another tau is no cycle
    assert X.stall_stats(apart, 40, 2)[1:3] == (1.0, 0.0)
    assert X.stall_stats([], 40, 0) == (0.0, 0.0, 0.0, 0.0)
    rows = X.stall_rows({'d1': [('sliding_path', normal, 10, 3), ('sliding_path', cycle, 6, 4)], 'signed': [('signed_path', swap, 5, 3)]})
    assert rows.shape == (3, 5) and list(rows[0]) == [2.0, 2.0, 5.0, 1.0, 4 / 6] and list(rows[1]) == [1.0, 1.0, 1.0, 0.0, 0.6] and not rows[2].any()
    assert X.STALL_KINDS == ('d1', 'signed', 'bridge') and X.STALL_FIELDS == ('calls', 'exits', 'max_repeat', 'stalled', 'fill')


def test_record_path_calls_sees_both_path_functions_and_puts_them_back():
    n = 12
    P = IM.PixelPSF(n, n, 1.2)
    lat = IM.Lattice(P, h=1.0, margin=1.0)
    y = P.synth(np.array([[5.2, 6.1]]), np.array([500.0]))
    sw = np.full((n, n), 0.1)
    real_sp, real_sg = IM.sliding_path, X.signed_path
    with X.record_path_calls() as calls:
        assert IM.sliding_path is not real_sp and X.signed_path is not real_sg
        R = IM.sliding_path(P, y, sw, lat, rho=1.0, max_entries=3, local=False)
        S = X.signed_path(P, y, sw, lat, max_entries=3)
    assert IM.sliding_path is real_sp and X.signed_path is real_sg
    assert [c[0] for c in calls] == ['sliding_path', 'signed_path'] and calls[0][2] == 3 and calls[1][2] == 3
    assert calls[0][1] == [(t, ij, kind) for t, ij, kind in R['events']] and calls[0][3] == len(R['taus'])
    assert calls[1][1] == [(t, ij, kind) for t, ij, kind, s in S['events']] and calls[1][3] == sum(1 for e in S['events'] if e[2] == 'enter') and calls[1][3] >= 1
    with pytest.raises(RuntimeError):
        with X.record_path_calls():
            raise RuntimeError('planted')
    assert IM.sliding_path is real_sp and X.signed_path is real_sg                      # restored after an exception too


def test_stall_section_lists_defective_paths(monkeypatch):
    monkeypatch.setattr(X, 'FRAMES', ('T1', 'T2'))
    good = np.array([[3.0, 0.0, 1.0, 0.0, 0.2], [1.0, 0.0, 1.0, 0.0, 0.1], [13.0, 1.0, 2.0, 0.0, 0.3]])
    bad = np.array([[2.0, 5.0, 5.0, 1.0, 1.0], [1.0, 0.0, 1.0, 0.0, 0.1], [0.0, 0.0, 0.0, 0.0, 0.0]])
    raw = {X.key_stall('hp', 'T1'): good, X.key_stall('plain', 'T1'): good}
    md = '\n'.join(X.sec_stall(raw))
    assert '(2 frame-variant pairs; 34 path calls: 6 of Definition 1, 2 of the signed path, 26 of the bridge\'s packaged pipeline): 2 exit events in all' in md
    assert 'is 2 (a stalled path has more than 2), stalled paths: 0.' in md and '| variant |' not in md
    raw[X.key_stall('plain', 'T2')] = bad
    md = '\n'.join(X.sec_stall(raw))
    assert 'stalled paths: 1 (listed below; defective paths, not path-versus-network results).' in md and 'is 5 (a stalled path has more than 2)' in md
    assert '| variant | frame | path | stalled calls | largest repeat |' in md and '| plain | T2 | of Definition 1 | 1 | 5 |' in md
    assert X.sec_stall({}) == ['## Stalled paths (the enter/exit limit cycle of `imaging.sliding_path`)', '', '(not run)', '']


def test_scaling_section_statistics_on_a_synthetic_raw(monkeypatch):
    names = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
    noisy = {'T1': True, 'T2': True, 'T3': True, 'T4': True, 'T5': False, 'T6': False}
    for n in names:
        monkeypatch.setitem(X.SPECS, n, dict(X.SPECS['S01'], name=n, noisy=noisy[n]))
    monkeypatch.setattr(X, 'FRAMES', names)
    raw = {}

    def put(fr, ok, rel, d, amp, dorm, ide, order, tdev):
        s = dict(zip(X.SN_FIELDS, np.zeros(len(X.SN_FIELDS))))
        s.update(ok=ok, rel_model=rel, max_d=d, amp_max=amp, dormant=dorm, id_equal=ide, order_ok=order, tdev=tdev)
        raw[X.key_net('hp', 1e-8, fr, 's_mid')] = np.array([s[f] for f in X.SN_FIELDS])
        raw[X.key_net('hp', 1e-8, fr, 'run')] = np.zeros(len(X.RUN_FIELDS))
    put('T1', 1.0, 1e-9, 1e-8, 1e-10, 1e-12, 1.0, 1.0, 1e-3)
    put('T2', 1.0, 2e-9, 3e-8, 2e-10, 2e-12, 1.0, 0.0, float('nan'))
    put('T3', 1.0, 9e-9, 2e-8, 9e-10, 9e-12, 1.0, 1.0, 3e-3)
    put('T4', 0.0, 5e-1, 5e-1, 5e-1, 5e-1, 0.0, 0.0, 5e-1)                      # a failure: not in the statistics
    put('T5', 1.0, 1e-13, 1e-12, 1e-12, 1e-13, 1.0, 1.0, 1e-3)
    put('T6', 0.0, 1e-1, 1e-1, 1e-1, 1e-1, 1.0, 1.0, 1e-1)
    L = X.sec_scaling(raw, 'hp', [1e-8])
    row_noisy = '| 1e-08 | noisy | 3 | 2.0e-09 | 9.0e-09 | 2.0e-08 | 3.0e-08 | 2.0e-10 | 9.0e-10 | 2.0e-12 | 9.0e-12 | 3 of 3 | 2 of 3 | 2.0e-03 | %.3g |' % (2e-3 / 1e-8 ** 0.5)
    row_free = '| 1e-08 | noise-free | 1 | 1.0e-13 | 1.0e-13 | 1.0e-12 | 1.0e-12 | 1.0e-12 | 1.0e-12 | 1.0e-13 | 1.0e-13 | 1 of 1 | 1 of 1 | 1.0e-03 | %.3g |' % (1e-3 / 1e-8 ** 0.5)
    assert row_noisy in L and row_free in L and len(L) == 7
    assert X.sec_scaling({}, 'hp', [1e-8])[-2] == '(none passes)'


@pytest.mark.skipif(not (os.path.exists(X.RAW) and os.path.exists(X.MD)), reason='the stored raw and report are not there')
def test_stored_report_is_the_report_of_the_stored_raw():
    assert X.report(X.RAW) == open(X.MD).read()


@pytest.mark.skipif(not os.path.exists(X.RAW), reason='the stored raw is not there')
def test_bridge_reproduces_the_stored_results():
    raw = X.load_raw(X.RAW)
    for fr in ('S02', 'S05', 'S08', 'S10', 'A06'):
        if 'B|%s|full' % fr not in raw:
            continue
        prob = X.build(X.SPECS[fr], X.HP)
        res = X.bridge_frame(prob, X.load_path(raw, 'hp', fr))
        for name, _ in X.BR_CONFIGS:
            assert np.allclose(res[name], raw['B|%s|%s' % (fr, name)], rtol=1e-9, atol=1e-12, equal_nan=True), (fr, name)
        assert res['replica'][0] == raw['B|%s|replica' % fr][0] == 0.0


def _with_a_later_block_of_the_same_size(real, k):
    """`sliding_path` returning, besides the real path, one more block with k units (the centres moved by 0.3 px) after the real ones."""
    def fake(*a, **kw):
        R = real(*a, **kw)
        h = [q for q in R['hist'] if len(q['w']) == k][0]
        return dict(R, hist=list(R['hist']) + [dict(h, c=h['c'] + 0.3, tau=R['hist'][-1]['tau'] * 2.0)])
    return fake


def test_the_first_block_with_the_wanted_number_of_units_is_taken(monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    pr = X.build(dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None), None)
    real = IM.sliding_path
    path = X.run_path(pr)
    monkeypatch.setattr(IM, 'sliding_path', _with_a_later_block_of_the_same_size(real, 1))
    path2 = X.run_path(pr)
    assert np.array_equal(path2['blk_c'], path['blk_c']) and path2['pv'][X.PV_FIELDS.index('i_stop')] == path['pv'][X.PV_FIELDS.index('i_stop')] == 1.0       # the noise-free stop: the first
    monkeypatch.setattr(IM, 'sliding_path', real)
    free = X.build(dict(X.SPECS['S05'], K=1, seed=79, noisy=True, name='T2', pair=None), X.HP)
    cfg = dict(X.BR_CONFIGS)['off4']
    ref = X.packaged_fit(free['sim'], free['P'], free['lat'], cfg, 1)
    monkeypatch.setattr(IM, 'sliding_path', _with_a_later_block_of_the_same_size(real, 1))
    again = X.packaged_fit(free['sim'], free['P'], free['lat'], cfg, 1)
    assert np.array_equal(again['cnt'][0], ref['cnt'][0]) and np.array_equal(again['cnt'][1], ref['cnt'][1])                                                      # the bridge's same-count block: the first


def test_advance_records_the_formation_time_between_accepted_steps(monkeypatch):
    n, eps = 8, 1e-3
    P = IM.PixelPSF(n, n, 1.2)
    net = X.Network(P, np.zeros((n, n)), np.ones((n, n)), np.array([[3.0, 4.0]]), eps)
    states = []
    for nu2 in (0.1, 20.0, 400.0):                                    # accepted states of a one-unit network: alpha 3, nu^2 = 0.1, 20, 400 at t = 1, 2, 3
        z = net.init()
        z[0], z[1] = 3.0, np.sqrt(nu2)
        states.append(z)

    class Stub:
        def __init__(self, fun, t0, y0, t1, **kw):
            self.status, self.t, self.y, self.nfev, self.njev, self.nlu, self.message, self.n = 'running', t0, np.array(y0), 1, 0, 0, '', 0

        def step(self):
            self.t, self.y, self.n = float(self.n + 1), states[self.n].copy(), self.n + 1
            if self.n == 3:
                self.status = 'finished'
    monkeypatch.setattr(X, 'Radau', Stub)
    tr = X.Tracker(1, 10.0)
    z, t, status, cnt = X.advance(net, net.init(), 0.0, 3.0, tr)
    w = [nu2 * np.tanh(3.0) / 12.0 for nu2 in (0.1, 20.0, 400.0)]
    assert w[1] < 10.0 < w[2] and status == 'ok' and cnt['steps'] == 3 and t == 3.0 and np.array_equal(z, states[2])
    assert np.isclose(tr.tform[0], 2.0 + (np.log(10.0) - np.log(w[1])) / (np.log(w[2]) - np.log(w[1])), rtol=1e-14, atol=0)


def test_dual_margin_follows_a_unit_that_enters_and_leaves():
    n = 12
    P = IM.PixelPSF(n, n, 1.2)
    lat = IM.Lattice(P, h=1.0, margin=1.0)
    y = np.random.default_rng(3).normal(size=(n, n)) + P.synth(np.array([[5.0, 6.0]]), np.array([30.0]))
    W = np.full((n, n), 0.5)
    prob = dict(P=P, y=y, W=W, lat=lat)
    b, a = 7, 6                                                           # the best lattice site of the emitter at (x, y) = (5, 6)
    empty = dict(c=np.zeros((0, 2)), w=np.zeros(0), ids=[])
    c1, w1 = lat.pos(b, a)[None, :], np.array([10.0])
    hist = [dict(empty, tau=0.0), dict(c=c1, w=w1, ids=[(b, a)], tau=1.0), dict(empty, tau=2.0)]          # the unit enters at tau = 1 and leaves at tau = 2
    best, mid, idx = X.dual_margin(prob, hist, 3.0, [])
    V0 = lat.corr(W * y)
    V1 = lat.corr(W * (y - P.synth(c1, w1)))
    taken = np.zeros(lat.shape, bool)
    taken[b, a] = True
    U1 = V0.copy()                                                       # tau 0 -> 1: the empty block
    U1[b, a] = 1.0                                                       # the entrant at the wall
    U2 = U1 + np.where(taken, 0.0, V1)                                   # tau 1 -> 2: the others grow with the block's residual
    U3 = U2 + V0                                                         # tau 2 -> 3: the unit has left, it grows from the wall
    assert np.isclose(mid, (1.0 - U3).min(), rtol=1e-12, atol=0) and idx == float(np.argmin(1.0 - U3)) and idx == b * lat.shape[1] + a
    assert np.isclose(best, min((1.0 - V0).min(), np.where(taken, np.inf, 1.0 - U2).min(), (1.0 - U3).min()), rtol=1e-12, atol=0)
    assert mid < 0.0 and best == mid                                     # the unit that left sits beyond the wall: a negative margin


def test_bridge_configurations_are_the_declared_ones():
    full = dict(r_ex=2.0, local=True, widths=True, fixed=False)
    assert X.R_EX == 2.0 and X.WIDTHS_REL == (0.7, 2.3) and X.BR_AMP == 0.05
    assert [n for n, _ in X.BR_CONFIGS] == ['full', 'no_rex', 'no_local', 'no_widths', 'no_reweight', 'off2', 'off3', 'off4']
    assert dict(X.BR_CONFIGS) == {'full': full, 'no_rex': dict(full, r_ex=0.0), 'no_local': dict(full, local=False), 'no_widths': dict(full, widths=False), 'no_reweight': dict(full, fixed=True),
                                  'off2': dict(r_ex=0.0, local=True, widths=False, fixed=False), 'off3': dict(r_ex=0.0, local=False, widths=False, fixed=False),
                                  'off4': dict(r_ex=0.0, local=False, widths=False, fixed=True)}
    assert set(X.BR_LABEL) == {n for n, _ in X.BR_CONFIGS} and X.BR_FIELDS == ('n', 'K_ref', 'same_count', 'max_d', 'amp_max', 'equal')


def test_bridge_replica_error_is_the_largest_over_the_configurations(monkeypatch):
    spec = dict(X.SPECS['S05'], K=1, seed=79, noisy=False, name='T1', pair=None)
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    pr = X.build(spec, X.HP)
    path = X.run_path(pr)
    real = X.packaged_reference
    cfgs = dict(X.BR_CONFIGS)

    def shifted(sim, psf, lat, cfg):
        c, w = real(sim, psf, lat, cfg)
        return (c + 1e-3, w) if cfg == cfgs['no_widths'] else (c, w)
    monkeypatch.setattr(X, 'packaged_reference', shifted)
    res = X.bridge_frame(pr, path)
    assert np.isclose(res['replica'][0], 1e-3, rtol=1e-6, atol=0)
    assert X.replica_error((np.zeros((2, 2)), np.ones(2)), (np.zeros((2, 3)), np.ones(2))) == 1.0


def test_signed_events_with_an_exit_are_counted_as_entries_only(monkeypatch):
    fake = dict(events=[(1.0, (1, 2), 'enter', 1.0), (2.0, (3, 4), 'enter', -1.0), (2.0, (1, 2), 'exit', 1.0)], hist=[])
    monkeypatch.setattr(X, 'signed_path', lambda *a, **k: fake)
    with X.record_path_calls() as calls:
        X.signed_path(None, None, None, None, max_entries=5)
    assert calls == [('signed_path', [(1.0, (1, 2), 'enter'), (2.0, (3, 4), 'enter'), (2.0, (1, 2), 'exit')], 5, 2)]


def test_eps_present_scans_all_variants_and_the_predictor_threshold_is_inclusive(tiny_raw, monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    monkeypatch.setitem(X.SPECS, 'T1', dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None))
    monkeypatch.setattr(X, 'FRAMES', ('T1',))
    monkeypatch.setattr(X, 'RTOL_CHECK', 1e-10)
    raw = dict(X.load_raw(tiny_raw))
    assert X.eps_present({'N|plain|1e-07|T1|run': 0, 'N|hp|1e-06|T1|run': 0, 'P|hp|T1|pv': 0}) == [1e-6, 1e-7]
    # the predictor: a margin of exactly 4 eps^(1/3) is "<= 4"
    pv = raw['P|hp|T1|pv'].copy()
    pv[X.PV_FIELDS.index('margin_min')] = 4.0 * 1e-6 ** (1.0 / 3.0)
    raw['P|hp|T1|pv'] = pv
    md = '\n'.join(X.sec_predictor(raw, 'hp', [1e-6]))
    assert '| 1e-06 | 2 | 0 | 0 | 1 | 0 |' in md and '| 1e-06 | 4 | 1 | 0 | 0 | 0 |' in md
    # an unfinished run in the solver table, a tolerance check against the other state, the largest replica error
    run = raw['N|hp|1e-06|T1|run'].copy()
    run[X.RUN_FIELDS.index('status')] = 1.0
    raw['N|hp|1e-06|T1|run'] = run
    md = '\n'.join(X.sec_solver(raw, [1e-6]))
    assert '| hp | 1e-06 | 1 | 0 | 0 |' in md
    P = X.build(X.SPECS['T1'], X.HP)['P']
    wa = raw['A|hp|1e-06|T1|w_mid'] * 1.5
    raw['A|hp|1e-06|T1|w_mid'] = wa
    Ma, Mb = P.synth(raw['A|hp|1e-06|T1|c_mid'], wa), P.synth(raw['N|hp|1e-06|T1|c_mid'], raw['N|hp|1e-06|T1|w_mid'])
    assert '| T1 | 1e-06 | %.1e |' % (np.linalg.norm(Ma - Mb) / np.linalg.norm(Mb)) in '\n'.join(X.sec_solver(raw, [1e-6]))
    raw['B|T1|replica'] = np.array([0.0])
    raw['B|T2|replica'] = np.array([5e-4])
    for k in [k for k in raw if k.startswith('B|T1|') and k != 'B|T1|replica']:
        raw[k.replace('T1', 'T2')] = raw[k]
    monkeypatch.setattr(X, 'FRAMES', ('T1', 'T2'))
    assert 'itself: 5.0e-04' in '\n'.join(X.sec_bridge(raw))


# ---- the raw, the parts and the report ------------------------------------------------------------------------------------------------------------------------------
def test_raw_roundtrip_merge_and_keys(tmp_path):
    a, b, c = str(tmp_path / 'a.npz'), str(tmp_path / 'b.npz'), str(tmp_path / 'c.npz')
    assert X.load_raw(a) == {}
    X.save_part(a, {'x': np.arange(3.0), 'y': np.array('text')})
    X.save_part(a, {'x': np.arange(4.0), 'z': np.eye(2)})
    r = X.load_raw(a)
    assert sorted(r) == ['x', 'y', 'z'] and np.array_equal(r['x'], np.arange(4.0)) and str(r['y']) == 'text' and not os.path.exists(str(tmp_path / 'a.tmp.npz'))
    X.save_part(b, {'x': np.zeros(2), 'w': np.ones(1)})
    assert X.merge_raws(c, [a, b]) == 4
    r = X.load_raw(c)
    assert sorted(r) == ['w', 'x', 'y', 'z'] and np.array_equal(r['x'], np.zeros(2))
    assert X.eps_tag(1e-6) == '1e-06' and X.eps_tag(1e-8) == '1e-08'
    assert X.key_path('hp', 'S01', 'pv') == 'P|hp|S01|pv' and X.key_net('plain', 1e-8, 'A02', 'run') == 'N|plain|1e-08|A02|run' and X.key_acc('hp', 1e-6, 'S02', 'w_mid') == 'A|hp|1e-06|S02|w_mid'


@pytest.fixture(scope='module')
def tiny_raw(tmp_path_factory):
    """A raw with one tiny noise-free frame (K = 1) through every part (the check part at a loose tolerance); the patches are undone when it is built (the tests patch again)."""
    path = str(tmp_path_factory.mktemp('raw') / 'raw.npz')
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(X.SMLM, 'ny', 10)
        mp.setitem(X.SMLM, 'nx', 10)
        mp.setitem(X.SPECS, 'T1', dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None))
        mp.setattr(X, 'RTOL_CHECK', 1e-10)
        X.part_P(['T1'], path)
        X.part_network('hp', [1e-6], ['T1'], path)
        X.part_network('plain', [1e-6], ['T1'], path)
        X.part_bridge(['T1'], path)
        X.part_check(['T1'], [1e-6], path)
        X.part_stall(['T1'], path)
    return path


def test_parts_skip_what_is_in_the_raw_and_force_recomputes(tiny_raw, monkeypatch, capsys):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    monkeypatch.setitem(X.SPECS, 'T1', dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None))
    monkeypatch.setattr(X, 'RTOL_CHECK', 1e-10)
    capsys.readouterr()
    X.part_P(['T1'], tiny_raw)
    X.part_network('hp', [1e-6], ['T1'], tiny_raw)
    X.part_bridge(['T1'], tiny_raw)
    X.part_check(['T1'], [1e-6], tiny_raw)
    X.part_stall(['T1'], tiny_raw)
    assert capsys.readouterr().out == ''                                           # nothing recomputed
    X.part_P(['T1'], tiny_raw, force=True)
    assert 'P hp T1: K_stop 1 of 1' in capsys.readouterr().out
    X.part_stall(['T1'], tiny_raw, force=True)
    assert 'S hp T1: ' in capsys.readouterr().out.split('P hp T1')[-1]
    raw = X.load_raw(tiny_raw)
    assert 'S|hp|T1' in raw and 'S|plain|T1' in raw and raw['S|hp|T1'].shape == (3, 5) and raw['S|hp|T1'][0, 0] == 1 and raw['S|hp|T1'][1, 0] == 1 and raw['S|hp|T1'][2, 0] >= 8 and raw['S|plain|T1'][2, 0] == 0 and raw['S|hp|T1'][:, 3].sum() == 0
    assert 'P|hp|T1|pv' in raw and 'P|plain|T1|smid_pv' in raw and 'N|hp|1e-06|T1|s_mid' in raw and 'N|plain|1e-06|T1|tform' in raw and 'B|T1|full' in raw and 'A|hp|1e-06|T1|w_mid' in raw
    p = X.load_path(raw, 'hp', 'T1')
    assert set(p) == set(X.PATH_KEYS) | {'sq1', 'smid'} and X.load_path(raw, 'hp', 'nope') is None
    assert set(raw) >= set(X.net_arrays('hp', 1e-6, 'T1', {k: 0 for k in X.NET_KEYS}))
    assert raw['A|hp|1e-06|T1|run'][X.RUN_FIELDS.index('rtol')] == 1e-10 and raw['N|hp|1e-06|T1|run'][X.RUN_FIELDS.index('rtol')] == X.RTOL
    prob = X.build(X.SPECS['T1'], X.HP)
    again = X.bridge_frame(prob, X.load_path(raw, 'hp', 'T1'))
    assert np.array_equal(again['full'], raw['B|T1|full'])
    seen = []
    monkeypatch.setattr(X, 'part_check', lambda frames, eps_list, raw_path, force=False: seen.append(list(frames)))
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'A', '--frames', 'S01,S02,S10'])
    X.main()
    assert seen == [['S02', 'S10']]


def test_report_is_deterministic_and_prints_the_stored_numbers(tiny_raw, monkeypatch):
    monkeypatch.setitem(X.SMLM, 'ny', 10)
    monkeypatch.setitem(X.SMLM, 'nx', 10)
    monkeypatch.setitem(X.SPECS, 'T1', dict(X.SPECS['S05'], K=1, seed=77, noisy=False, name='T1', pair=None))
    monkeypatch.setattr(X, 'FRAMES', ('T1',))
    monkeypatch.setattr(X, 'GROUPS', (('SMLM, noise-free', ['T1']), ('all frames', ['T1'])))
    monkeypatch.setattr(X, 'RTOL_CHECK', 1e-10)
    md = X.report(tiny_raw)
    assert md == X.report(tiny_raw) and md.startswith('# Does the factorised RBF network')
    raw = X.load_raw(tiny_raw)
    s = dict(zip(X.SN_FIELDS, raw['N|hp|1e-06|T1|s_mid']))
    run = dict(zip(X.RUN_FIELDS, raw['N|hp|1e-06|T1|run']))
    assert s['ok'] == 1.0
    assert '| hp | 1e-06 | SMLM, noise-free | 1 of 1 | 1 of 1 | 1 of 1 | 1 of 1 |' in md and '| plain | 1e-06 | all frames | 1 of 1 | 1 of 1 | 1 of 1 | 1 of 1 |' in md
    row = '| T1 | 1 | 1 (0) | 0 | %.2e | %.2e | %.2e | %.1e | %.1e | YYYY | pass | %.1e | %d | %d |' % (s['max_d'], s['rel_model'], s['amp_max'], s['dormant'], s['kkt_pos'], s['tdev'], run['steps'], round(run['secs']))
    assert row in md
    pv = dict(zip(X.PV_FIELDS, raw['P|hp|T1|pv']))
    assert '| T1 | SMLM M, noise-free | 1 | 1 | %.4g | inf | n/a | cap | yes | 1 | %.2e | %.2e | inf | no |' % (pv['tau_K'], pv['margin_min'], pv['margin_mid']) in md
    assert '### Failures at tau_mid, high-passed (hp = 3 px), eps = 1e-06\n\nnone' in md and 'no frame has a negative entry of the signed path before tau_mid' in md
    assert '### The discrepancy against eps, high-passed (hp = 3 px), non-negative path at tau_mid\n\n(needs two eps)' in md
    assert X.eps_present(raw) == [1e-6] and '| 1e-06 | 1 | 1 | 0 | 0 | 0 | 0 |' not in md
    assert '### Failures against the runner-up margin, high-passed (hp = 3 px)' in md
    mu = pv['margin_mid']
    obs = s['dormant'] * np.abs(raw['N|hp|1e-06|T1|w_mid']).max()
    pred = X.LAW.core_tail(1e-6, mu)
    assert 0 < mu < 1 and abs(pred - 12.0 ** -1 * 1e-9 * np.cos(np.pi * mu / 2) / np.sqrt(2 * 0.025) / np.sin(np.pi * mu / 2) ** 3) < 1e-12 * pred
    st = raw['S|hp|T1'] + raw['S|plain|T1']
    assert '## Stalled paths (the enter/exit limit cycle of `imaging.sliding_path`)\n\nCheck of the paths of the declared runs (2 frame-variant pairs; %d path calls: %d of Definition 1, %d of the signed path, %d of the bridge\'s packaged pipeline): %d exit events in all' % (
        st[:, 0].sum(), st[0, 0], st[1, 0], st[2, 0], st[:, 1].sum()) in md and 'stalled paths: 0.' in md
    td = s['tdev']
    assert '### Size of the discrepancy of the runs that pass, high-passed (hp = 3 px)' in md
    assert '| 1e-06 | noise-free | 1 | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | 1 of 1 | 1 of 1 | %.1e | %.3g |' % (
        s['rel_model'], s['rel_model'], s['max_d'], s['max_d'], s['amp_max'], s['amp_max'], s['dormant'], s['dormant'], td, td / 1e-6 ** 0.5) in md
    assert '### Dormant units against the theory, high-passed (hp = 3 px)' in md and '| T1 | %.2e | %.2e | %.2e | %.2g |' % (mu, obs, pred, obs / pred) in md
    assert '### Duplicated units (post hoc), high-passed (hp = 3 px), eps = 1e-06\n\nnone' in md
    mm = pv['margin_min'] / 1e-6 ** (1.0 / 3.0)
    for t in X.PRED_THRESHOLDS:
        assert '| 1e-06 | %g | %d | 0 | %d | 0 |' % (t, int(mm <= t), int(mm > t)) in md
    assert '| T1 | 1e-06 |' in md and 'Tolerance check (rtol 1e-10 against the main runs at rtol %.0e)' % X.RTOL in md
    assert '| equal in |  | 1 of 1 |' in md and 'Largest difference between this script\'s replica of `fit_frame` and `fit_frame` itself: 0.0e+00' in md
    assert '| hp | 1e-06 | 1 | 1 | 0 |' in md
    out = md.split('## The solver')[1]
    assert out.count('| T1 |') >= 1
    assert X._x(float('nan')) == 'n/a' and X._x(float('inf')) == 'inf' and X._x(0.125, '%.1f') == '0.1' and X._x(3, '%d') == '3'
    assert X._tab(['a', 'b'], [[1, 'x']]) == ['| a | b |', '|---|---|', '| 1 | x |']
    assert X._flags(dict(c1=1.0, c2=0.0, c3=1.0, c4=0.0)) == 'YnYn'
    assert X._frame_desc('S06') == 'SMLM M, pair 2.5 sigma' and X._frame_desc('S03n') == 'SMLM L, noise-free' and X._frame_desc('A03') == 'astro, satellite 2.5 sigma' and X._frame_desc('A01') == 'astro'


def test_command_line(tmp_path, monkeypatch, capsys):
    raw = str(tmp_path / 'r.npz')
    md = str(tmp_path / 'r.md')
    monkeypatch.setattr(X, 'RAW', raw)
    monkeypatch.setattr(X, 'MD', md)
    monkeypatch.setattr(sys, 'argv', ['x'])
    with pytest.raises(SystemExit):
        X.main()
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'P', '--frames', 'S01', '--from-raw'])
    X.main()
    assert os.path.exists(md) and 'P|hp|S01|pv' in X.load_raw(raw) and 'P|plain|S01|pv' in X.load_raw(raw)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'P', '--frames', 'S99'])
    with pytest.raises(SystemExit):
        X.main()
    other = str(tmp_path / 'o.npz')
    X.save_part(other, {'k': np.arange(2.0)})
    monkeypatch.setattr(sys, 'argv', ['x', '--merge', other])
    X.main()
    assert 'k' in X.load_raw(raw) and 'merged' in capsys.readouterr().out
