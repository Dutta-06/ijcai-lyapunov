"""
Proposition D‴.37 (PARADIGM.md §5d): the dither's drift along a manifold of equilibria. Checked here: the second-order identities of the
limit map against flows whose limit map is known in closed form (a curved coordinate system on a linear flow; a curved manifold), and
against the limit map computed by integration for a nonlinear gradient flow in a non-orthogonal metric; on Proposition D‴.36's instance,
the diffusion against D‴.36's projection formula, and the drift at small γ against the replicator flow of the log-determinant.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
from src.core import gauge_drift as GDR                                    # noqa: E402
import exp_gauge_dither as GD                                             # noqa: E402
import exp_gauge_drift as X                                               # noqa: E402


def test_curved_coordinates_on_a_linear_flow():
    # x′ = 0, z′ = −λz written in u = x + κz², v = z: Φ(u, v) = (u − κv², 0), so ½∂²Φ_u[Σ] = −κΣ_vv
    lam, kap = 0.7, 0.4
    F = lambda y: np.array([-2 * kap * lam * y[1] ** 2, -lam * y[1]])
    chart = lambda w: np.array([w[0], 0.0])
    coords = lambda y: np.array([y[0]])
    S = np.diag([3e-4, 2e-4])
    d, D, ev = GDR.limit_drift(F, chart, coords, np.array([0.2]), S)
    assert abs(d[0] / (-kap * S[1, 1]) - 1) < 1e-6 and abs(D[0, 0] / S[0, 0] - 1) < 1e-8 and abs(ev[0] + lam) < 1e-8


def test_a_curved_manifold():
    # x′ = 0, z′ = −λz in u = x, v = z + qx²: M = {v = qu²}, Φ = (u, qu²); the drift of v is qΣ_uu (the chart's curvature), of u zero
    lam, q, u0 = 0.5, 0.8, 0.3
    F = lambda y: np.array([0.0, -lam * (y[1] - q * y[0] ** 2)])
    chart = lambda w: np.array([w[0], q * w[0] ** 2])
    S = np.diag([2e-4, 1e-4])
    dv, _, _ = GDR.limit_drift(F, chart, lambda y: np.array([y[1]]), np.array([u0]), S)
    du, _, _ = GDR.limit_drift(F, chart, lambda y: np.array([y[0]]), np.array([u0]), S)
    assert abs(dv[0] / (q * S[0, 0]) - 1) < 1e-6 and abs(du[0]) < 1e-10


def test_against_the_integrated_limit_map():
    # a gradient flow of J = ½(1 + u/2)(v − qu²)² in the metric g: the limit map by integration, its second differences along the noise
    q = 0.6
    g = np.array([[1.0, 0.3], [0.3, 2.0]]); gi = np.linalg.inv(g)

    def F(y):
        u, v = y
        r = v - q * u * u
        grad = np.array([0.25 * r * r - (1 + 0.5 * u) * r * 2 * q * u, (1 + 0.5 * u) * r])
        return -gi @ grad
    chart = lambda w: np.array([w[0], q * w[0] ** 2])
    coords = lambda y: np.array([y[0]])
    u0 = 0.4
    S = np.diag([2e-4, 1e-4])
    d, _, _ = GDR.limit_drift(F, chart, coords, np.array([u0]), S)
    x0 = chart(np.array([u0])); hh = 2e-2
    bf = 0.0
    for i in range(2):
        e = np.eye(2)[i] * hh
        Pp = GDR.limit_map(F, x0 + e, 400.0)[0][0]; Pm = GDR.limit_map(F, x0 - e, 400.0)[0][0]
        bf += 0.5 * S[i, i] * (Pp - 2 * u0 + Pm) / hh ** 2
    assert abs(d[0] / bf - 1) < 2e-3


def test_the_instances_diffusion_is_the_projection_formula():
    w = np.array([0.30, 0.35])
    d, D, ev = X.drift_at(w)
    D3 = GD.gauge_diffusion(np.r_[w, 1 - w.sum()], X.GAMMA, X.TAU)
    assert np.abs(D - D3[:2, :2]).max() < 1e-4 * np.abs(D3).max() and np.all(ev < 0)


def test_the_slow_replicator_limit_is_the_log_determinants_replicator_flow():
    for w in (np.array([1 / 3, 1 / 3]), np.array([0.40, 0.25])):
        d = X.drift_at(w, gamma=1e-5)[0] / (1e-5 * X.TAU)
        law = X.log_det_law(np.r_[w, 1 - w.sum()])[:2]
        assert np.abs(d - law).max() < 1e-2 * np.abs(law).max()
    d2 = X.drift_at(np.array([1 / 3, 1 / 3]), gamma=1e-2)[0] / (1e-2 * X.TAU)    # at γ = 10⁻² the finite-γ correction is O(1)
    assert np.abs(d2 - X.log_det_law(np.ones(3) / 3)[:2]).max() > 0.1


def test_the_log_determinant_does_not_see_the_split_within_a_class():
    # two coincident particles in the first class: the fit's block is rank p along the class's β-weighted mean, its Schur complement the
    # cohesion's on the separation, so det ∇²_ΘJ depends on the class's total mass only
    import exp_cohesion_leak as E
    from src.core import mass_transfer as MT
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    T3 = E.exact_rep(net, np.ones(3) / 3)
    x0 = np.vstack([T3[0], T3[0], T3[1], T3[2]]).ravel()
    out = []
    for t in (0.5, 0.2):
        beta = np.array([t, 1 - t, 1.0, 1.0]); beta *= 4 / beta.sum()
        J = lambda x: loop.J(x.reshape(4, 3), beta)
        h = 1e-4; I = np.eye(12) * h; H = np.zeros((12, 12))
        for i in range(12):
            for j in range(i, 12):
                H[i, j] = H[j, i] = (J(x0 + I[i] + I[j]) - J(x0 + I[i] - I[j]) - J(x0 - I[i] + I[j]) + J(x0 - I[i] - I[j])) / (4 * h * h)
        out.append(np.linalg.slogdet(H))
    assert out[0][0] > 0 and out[1][0] > 0 and abs(out[0][1] - out[1][1]) < 1e-5


def test_the_ito_term_of_a_nonlinear_observable():
    # round 27 (D3z-8): on the first toy, observe w = u²: dw = 2u·(−κΣ_vv) + ½·2·Σ_uu (the Itô term of the observable)
    lam, kap, u0 = 0.7, 0.4, 0.2
    F = lambda y: np.array([-2 * kap * lam * y[1] ** 2, -lam * y[1]])
    chart = lambda w: np.array([np.sqrt(w[0]), 0.0])
    coords = lambda y: np.array([y[0] ** 2])
    S = np.diag([3e-4, 2e-4])
    d, D, _ = GDR.limit_drift(F, chart, coords, np.array([u0 ** 2]), S)
    assert abs(d[0] / (2 * u0 * (-kap * S[1, 1]) + S[0, 0]) - 1) < 1e-5 and abs(D[0, 0] / (4 * u0 ** 2 * S[0, 0]) - 1) < 1e-6


def test_the_instance_at_gamma_1e2_against_the_integrated_limit_map():
    # round 27 (D3z-8): Part A's brute-force value at W = 1/3 (LSODA to 1.2·10⁵, second differences at ±0.02), drift of (W_1, W_2)
    d = X.drift_at(np.array([1 / 3, 1 / 3]))[0]
    assert np.abs(d - np.array([-2.05497e-7, 1.57975e-7])).max() < 2e-3 * 2.05e-7


def test_the_log_determinant_closed_form():
    # round 27 (D3z-6): with one particle per class ℓ = 2Σ log(W_k − W_min,k²/W_k) + const
    Wm = X.W_MIN
    vals = [X.log_det(3 * W) - 2 * np.sum(np.log(W - Wm ** 2 / W)) for W in (np.ones(3) / 3, np.array([0.27, 0.36, 0.37]))]
    assert abs(vals[0] - vals[1]) < 1e-6


def test_the_normal_rates_and_the_leaks_window():
    # referee round 143 (D3el-1, D3el-5): the slowest normal mode at W = 1/3 is a width mode that barely carries mass; the leak's stiffness
    # along M has its eigenvalues in (0, ε) (κ = ε(TᵀgT)⁻¹H_W with 0 ≺ H_W ⪯ TᵀgT); at fixed γ the normal rates do not degenerate at the face,
    # while at small γ the mass channel's rate γ-small sets the slowest one there
    r = X.normal_rates(np.ones(2) / 3, 1e-2)
    assert 2.9e-4 < r['rates'][0] < 3.1e-4 and r['mass'] < 1e-3 and r['pos'][2] > 0.9 and np.all(np.diff(r['rates']) > 0)
    for eps in X.EPS_G:
        ev = np.linalg.eigvals(X.kappa(np.ones(3) / 3, eps)).real
        assert np.all(ev > 0) and np.all(ev < eps)
    near = X.normal_rates(np.array([X.W_MIN[0] + 5e-4, (1 - X.W_MIN[0] - 5e-4) / 2]), 1e-2)
    assert near['rates'][0] > 0.9 * r['rates'][0]
    slow = X.normal_rates(np.array([X.W_MIN[0] + 5e-4, (1 - X.W_MIN[0] - 5e-4) / 2]), 1e-4)
    assert slow['rates'][0] < 0.2 * r['rates'][0]


def test_part_I_on_the_stored_field_and_ensemble(monkeypatch):
    # Part I's ensemble comparison against a stub limit SDE whose paths sit at known points: the means, the combined distance and the drifts
    monkeypatch.setattr(X, 'normal_rates', lambda w, g, n=4: dict(w=np.r_[w, 1 - np.sum(w)], gamma=g, rates=np.ones(n), mass=0.0, pos=np.ones(3) / 3))
    monkeypatch.setattr(X, 'instance', lambda g=X.GAMMA: (None, lambda w: np.zeros(11), None, None))
    monkeypatch.setattr(X, 'drift_at', lambda w, g=X.GAMMA, by_noise=False: (np.array([-2e-7, 1e-7]), np.diag([4e-8, 2e-8]), None))
    Y = np.tile([0.32, 0.34, 0.34], (100, 1)); Y[::2] += [0.01, -0.01, 0.0]
    monkeypatch.setattr(X, 'sde_paths', lambda field, S, t, seed=0: {u: Y for u in t})
    monkeypatch.setattr(X, 'kappa', lambda W, eps: np.diag([2 * eps, eps]))      # eigvals in the diagonal's order: the md lists them sorted
    Ff = dict(s=[1e4], mean=[np.array([0.33, 0.335, 0.335])], se=[np.full(3, 0.002)], path=[np.array([0.331, 0.334, 0.335])])
    out = X.part_I(np.zeros((3, 8)), Ff)
    r = out['sde'][0]
    ys = Y.std(axis=0, ddof=1) / 10
    assert np.allclose(r['sde'], Y.mean(axis=0)) and np.allclose(r['z'], (Ff['mean'][0] - Y.mean(axis=0)) / np.hypot(0.002, ys))
    assert np.allclose(r['rate_loop'], (Ff['mean'][0] - 1 / 3) / (1e4 * X.GAMMA * X.TAU))
    pe = out['pe'][0]
    assert abs(pe['pe'] - 2e-7 * (1 / 3 - X.W_MIN[0]) / 4e-8) < 1e-12 and abs(pe['spread'] - np.sqrt(4e-8 / (pe['gamma'] * X.TAU))) < 1e-12
    assert all(np.array_equal(st['eig'], [st['eps'], 2 * st['eps']]) for st in out['stiff'])
    assert '## I.' in X.report(Ii=out)


def _same(a, b, path='', rtol=1e-9, atol=1e-14):
    """a and b equal as nested dicts and lists of numbers and arrays (to rtol, atol), with the same keys and lengths."""
    if isinstance(b, dict):
        assert isinstance(a, dict) and set(a) == set(b), path
        for k in b:
            _same(a[k], b[k], f'{path}/{k}', rtol, atol)
    elif isinstance(b, (list, tuple)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f'{path}[{i}]', rtol, atol)
    else:
        x, y = np.asarray(a, float), np.asarray(b, float)
        assert x.shape == y.shape and np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True), path


def test_part_I_is_recomputed_from_the_raw():
    # Part I's numbers are stored in the raw (the md's regeneration only formats them): every one recomputed from the stored field E and
    # ensemble F — the normal rates and the slowest mode's shares, the leak's stiffness and spread, the face rates and α_1, the Péclet
    # number and the free spread, the limit SDE's ensemble (its own seed) against F, Katzenberger's drift at 1/3
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(X.OUT, 'gauge_drift_raw.npz'))
    if 'I' not in raw:
        import pytest
        pytest.skip('no Part I in the raw file')
    I = raw['I']
    assert [r['gamma'] for r in I['rates']] == list(X.GAMMAS_I) and [s['eps'] for s in I['stiff']] == list(X.EPS_G)
    assert [(r['gamma'], r['delta']) for r in I['face']] == [(g, d) for g in X.GAMMAS_IF for d in X.DELTAS_I]
    assert [r['gamma'] for r in I['pe']] == list(X.GAMMAS_PE) and [r['s'] for r in I['sde']] == list(raw['F']['s'])
    _same(X.part_I(raw['E']['field'], raw['F']), I)
    md = open(os.path.join(X.OUT, 'gauge_drift.md')).read()
    assert '\n'.join(X.report_I(I)) in md


def test_normal_rates_shares_on_a_planted_linearisation(monkeypatch):
    # normal_rates on a planted Jacobian: two tangent zeros, then the normal rates in order of |λ|; the slowest normal mode's g-norm split
    # into the positions' (Σv²/G by coordinate α, c, s) and the masses' (Σδβ²/(γβ), δβ from δ log(β_k/β_3) projected on Σδβ = 0)
    import src.core.gauge_drift as GDRm
    W = np.array([0.30, 0.35]); Wf = np.r_[W, 1 - W.sum()]; beta = 3 * Wf; g = 1e-2
    rng = np.random.default_rng(5)
    v = rng.standard_normal(11)
    lam = np.r_[0.0, 0.0, -3e-4, -2e-3, -2.5e-3, -1.6e-2, -np.linspace(0.05, 0.5, 5)]
    V = np.linalg.qr(rng.standard_normal((11, 11)))[0]; V[:, 2] = v
    A = V @ np.diag(lam) @ np.linalg.inv(V)
    monkeypatch.setattr(GDRm, '_jac', lambda F, x, h: A)
    monkeypatch.setattr(X, 'instance', lambda gamma=X.GAMMA: (None, lambda w: np.zeros(11), None, None))
    r = X.normal_rates(W, g)
    assert np.allclose(r['rates'], [3e-4, 2e-3, 2.5e-3, 1.6e-2], rtol=1e-8) and np.allclose(r['w'], Wf)
    import exp_cohesion_leak as E
    vv = v / np.linalg.norm(v)                                                # the shares are ratios: the eigenvector's scale drops out
    dz = np.r_[vv[9:], 0.0]; db = beta * (dz - Wf @ dz)
    pos = np.array([sum(vv[3 * k + j] ** 2 / E.G[j] for k in range(3)) for j in range(3)])
    mas = np.sum(db ** 2 / (g * beta))
    assert np.allclose(r['pos'], pos / (pos.sum() + mas), rtol=1e-6) and abs(r['mass'] - mas / (pos.sum() + mas)) < 1e-6 * r['mass']
