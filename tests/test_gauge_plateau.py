"""
Proposition D‴.39 (PARADIGM.md §5d): dissipativity of the dithered gauge. Checked here: the dead zone V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹}
leaves the drift and the energy unchanged inside the box and adds minus its G-gradient, −σ_B(ϑ_i − Π_B(ϑ_i)), outside; J is bounded on
the unconstrained coordinates (so e^{−J/τ} is not integrable) and J + V grows like (σ_B/2)dist²; the box's wall lies inside the saturation,
W♯ is representable in it, and the polytope's vertices have two classes at their walls.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402

BOX = (5.0, 3.0, 4.0)


def test_the_dead_zone_is_inactive_inside_the_box():
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3) + 0.05 * np.random.default_rng(1).standard_normal((3, 3))
    assert np.all(np.abs(Th) < np.array(BOX)[None, :])
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    assert np.array_equal(free.drift(Th, beta)[0], boxed.drift(Th, beta)[0]) and free.J(Th, beta) == boxed.J(Th, beta)


def test_the_dead_zone_force_is_minus_the_g_gradient_of_its_energy():
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3)
    Th[0, 0] = 7.5; Th[2, 1] = -3.4; Th[1, 2] = 4.6                          # three coordinates beyond the box
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    dV = boxed.drift(Th, beta)[0] - free.drift(Th, beta)[0]
    V = lambda T: boxed.J(T, beta) - free.J(T, beta)
    h = 1e-5
    grad = np.array([(V(Th + h * e.reshape(3, 3)) - V(Th - h * e.reshape(3, 3))) / (2 * h) for e in np.eye(9)]).reshape(3, 3)
    assert np.allclose(dV, -E.G[None, :] * grad, atol=1e-8)
    exc = np.zeros((3, 3)); exc[0, 0] = 2.5; exc[2, 1] = -0.4; exc[1, 2] = 0.6
    assert np.allclose(dV, -0.3 * exc, atol=1e-12)


def test_the_energy_is_bounded_on_the_unconstrained_coordinates_and_the_dead_zone_confines():
    # (i): the physical coordinates are bounded and the cohesion energy W ≤ R²/8, so J is bounded — e^{−J/τ} is not integrable; with the dead
    # zone, J + V ≥ (σ_B/2)dist_{G⁻¹}(Θ, B)² − sup|J|
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    m = 3
    sup_r2 = float(np.mean((np.abs(base.F) + net.a_max) ** 2))                # |f̂| ≤ (1/m)Σβ_i a_max = a_max
    bound = 0.5 * m * sup_r2 + 0.5 * free.sigma0 * m * free.R ** 2 / 8
    rng = np.random.default_rng(5)
    for scale in (1.0, 10.0, 1e3, 1e6):
        for _ in range(50):
            Th = scale * rng.standard_normal((3, 3))
            J = free.J(Th, beta)
            assert 0.0 <= J <= bound
            d2 = float(np.sum(boxed.box_excess(Th) ** 2 / E.G[None, :]))
            assert boxed.J(Th, beta) >= 0.15 * d2 - 1e-12


def test_the_box_keeps_the_best_conditioned_fit_and_moves_the_wall_inside_the_saturation():
    import exp_gauge_fixman as Q
    import exp_gauge_plateau as P
    ws = Q.w_sharp()
    assert np.all(P.W_BOX > P.W_MIN) and np.all(ws > P.W_BOX)                # the wall inside the saturation, W♯ representable in the box
    a_max = E.setup()[0].a_max
    assert np.all(np.abs(E.A_BUMP) / ws < a_max * np.tanh(BOX[0] / a_max))   # its amplitudes inside the box
    V = P.vertices()
    assert np.allclose(V.sum(axis=1), 1) and np.all(V >= P.W_BOX[None, :] - 1e-15)
    assert np.allclose(np.sort(V, axis=None)[:6], np.sort(np.r_[P.W_BOX, P.W_BOX]))


def test_a_dithered_run_inside_the_box_is_the_unboxed_run():
    # inside B the dead zone adds exactly zero, so with the same noise the two runs agree to the last bit until a coordinate leaves B
    net, base = E.setup()
    out = []
    for box in (None, BOX):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_A=E.A_BUMP.copy(), box=box)
        rec = []
        loop.run_dithered(E.exact_rep(net, np.ones(3) / 3), np.ones(3), 500.0, 0.25, 1e-4, np.random.default_rng(7), rec_every=200,
                          callback=lambda s_, T, b: rec.append(np.r_[T.ravel(), b]))
        out.append(np.array(rec))
    assert np.all(np.abs(out[0][:, :9]).reshape(-1, 3, 3) < np.array(BOX)) and np.array_equal(out[0], out[1])



def test_part_G_domain_and_arrhenius_on_plants():
    # referee round 143 (D3el-6): the domain's boundary function, and the prefactor fits recover a planted rate law
    import exp_gauge_plateau as P
    net, _ = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    u = np.random.default_rng(2).standard_normal((5, 9)); u /= (np.abs(u) ** P.DOM_P).sum(axis=1, keepdims=True) ** (1 / P.DOM_P)
    pts = Th0[None] + (np.tile(np.asarray(P.DOM_L), 3)[None] * u).reshape(5, 3, 3)
    assert np.allclose(P.dom_phi(pts, Th0), 1.0) and P.dom_phi(Th0[None], Th0)[0] == 0.0
    B, dJ = 6.0e-4, 6.3e-4
    for c in (0.0, 0.5, 1.0):
        rows = [dict(tau=t, n=n, rate=2.0 * t ** (-c) * np.exp(-B / t)) for t, n in zip(P.TAUS_G, (40, 30, 20, 15, 10, 7))]
        f = P.arrhenius(rows, dJ, c)
        assert abs(f['B'] - B) < 1e-12 and f['chi2'] < 1e-15 and abs(f['t4'] - 1e-4 ** c * np.exp(B / 1e-4) / 2.0) < 1e-6 * f['t4']
        g = P.arrhenius(rows, dJ)
        assert abs(g['c'] - c) < 1e-8 and abs(g['B'] - B) < 1e-12 and g['dof'] == len(rows) - 3
        w = [P.arrhenius(rows, dJ, cc)['chi2'] for cc in (0.0, 0.5, 1.0)]
        assert np.argmin(w) == (0.0, 0.5, 1.0).index(c)
    # the standard error is the weighted one: 1/√(Σn(x − x̄)²)
    x = 1 / np.array(P.TAUS_G); n = np.array([40, 30, 20, 15, 10, 7.0]); xb = n @ x / n.sum()
    assert abs(f['B_se'] - 1 / np.sqrt(n @ (x - xb) ** 2)) < 1e-15


def test_part_G_numbers_are_recomputed_from_the_raw():
    import exp_gauge_plateau as P
    raw = E.load(os.path.join(P.OUT, 'gauge_plateau_raw.npz'))
    if 'G' not in raw:
        import pytest
        pytest.skip('no Part G in the raw file')
    G, runs = raw['G'], raw['G_runs']
    assert len(runs) == len(P.TAUS_G) * len(P.SEEDS_G)
    for r in G['rows']:
        rr = [x for x in runs if x['tau'] == r['tau']]
        ex = [x['t'] for x in rr if x['t'] is not None]
        T = sum(x['t'] if x['t'] is not None else x['cap'] for x in rr)
        assert r['n'] == len(ex) and abs(r['T'] - T) < 1e-6 and abs(r['rate'] - len(ex) / T) < 1e-18
        assert all(0 < t <= x['cap'] for x, t in ((x, x['t']) for x in rr if x['t'] is not None))
    for f, c in zip(G['fits'][:3], (0.0, 0.5, 1.0)):                          # an independent weighted fit (np.polyfit, weights √n)
        ok = [r for r in G['rows'] if r['n'] > 0]
        x = np.array([1 / r['tau'] for r in ok]); y = np.log([r['rate'] for r in ok]) + c * np.log([r['tau'] for r in ok])
        assert abs(-np.polyfit(x, y, 1, w=np.sqrt([r['n'] for r in ok]))[0] - f['B']) < 1e-12
    md = open(os.path.join(P.OUT, 'gauge_plateau.md')).read()
    assert '\n'.join(P.report_G(G)) in md


def _same(a, b, path='', rtol=1e-9, atol=1e-14):
    """a and b equal as nested dicts and lists of numbers, arrays and None (numbers to rtol, atol)."""
    if isinstance(b, dict):
        assert isinstance(a, dict) and set(a) == set(b), path
        for k in b:
            _same(a[k], b[k], f'{path}/{k}', rtol, atol)
    elif isinstance(b, (list, tuple)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f'{path}[{i}]', rtol, atol)
    elif b is None:
        assert a is None, path
    else:
        x, y = np.asarray(a, float), np.asarray(b, float)
        assert x.shape == y.shape and np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True), path


def _raw_G():
    import exp_gauge_plateau as P
    raw = E.load(os.path.join(P.OUT, 'gauge_plateau_raw.npz'))
    if 'G' not in raw:
        import pytest
        pytest.skip('no Part G in the raw file')
    return P, raw


def test_part_G_is_recomputed_in_full_from_the_stored_runs():
    # every number of Part G's rows and fits recomputed from the stored domain and exit times; the runs' set-up against the constants
    P, raw = _raw_G()
    _same(P.part_G(raw['G_dom'], raw['G_runs']), raw['G'])
    assert [(r['tau'], r['seed']) for r in raw['G_runs']] == [(t, s) for t in P.TAUS_G for s in P.SEEDS_G]
    assert all(r['cap'] == c for r in raw['G_runs'] for t, c in zip(P.TAUS_G, P.S_G) if r['tau'] == t)
    assert [f['c'] for f in raw['G']['fits'][:3]] == [0.0, 0.5, 1.0] and raw['G']['fits'][3]['c_se'] is not None


def test_domain_job_reproduces_the_stored_domain(monkeypatch):
    # the barrier on ∂𝒪 and the flow back, at the stored set-up: the second L-BFGS start reaches the stored minimum (the other 198 replaced
    # by a stub that keeps the generator's stream, so the 63 random boundary points are the stored ones: their (J − J(exact))/τ); then the
    # deterministic flow from the minimiser alone over the full S_FLOW, its final distance and energy against the stored first entries
    P, raw = _raw_G()
    import scipy.optimize as so
    dom = raw['G_dom']
    real, calls = so.minimize, []

    def stub(f, x0, **kw):
        calls.append(1)
        return real(f, x0, **kw) if len(calls) <= 2 else so.OptimizeResult(fun=np.inf, x=x0)
    monkeypatch.setattr(so, 'minimize', stub)
    monkeypatch.setattr(P, 'S_FLOW', 0.0)
    d = P.domain_job()
    assert len(calls) == P.N_DOM == 200 and len(d['dJ_start']) == P.N_FLOW == 64
    assert abs(d['dJ'] / dom['dJ'] - 1) < 1e-9 and np.allclose(d['at'], dom['at'], rtol=0, atol=1e-6) and abs(d['fit'] - dom['fit']) < 1e-9
    assert np.allclose(d['a'], dom['a'], atol=1e-6) and np.allclose(d['c'], dom['c'], atol=1e-6) and np.allclose(d['s'], dom['s'], atol=1e-6)
    assert np.allclose(d['dJ_start'], dom['dJ_start'], rtol=1e-6)
    monkeypatch.setattr(so, 'minimize', real)
    monkeypatch.setattr(P, 'S_FLOW', 4e4); monkeypatch.setattr(P, 'N_DOM', 2); monkeypatch.setattr(P, 'N_FLOW', 1)
    d = P.domain_job()
    assert abs(d['dist_end'][0] / dom['dist_end'][0] - 1) < 1e-4 and abs(d['dJ_end'][0] / dom['dJ_end'][0] - 1) < 1e-4
    assert abs(P.dom_phi(np.asarray(dom['at'])[None] + P.E.exact_rep(P.E.setup()[0], np.ones(3) / 3)[None],
                         P.E.exact_rep(P.E.setup()[0], np.ones(3) / 3))[0] - 1) < 1e-9


def test_escape_runs_reproduce_the_stored_exits(monkeypatch):
    # escape_runs at its own set-up for the seed that exits first at the two largest τ (capped at 3000 time units): the stored exit times
    P, raw = _raw_G()
    monkeypatch.setattr(P, 'TAUS_G', (2.5e-4, 2e-4)); monkeypatch.setattr(P, 'S_G', (3000.0, 3000.0))
    monkeypatch.setattr(P, 'SEEDS_G', (716,)); monkeypatch.setattr(P, 'CKPT', None)
    out = P.escape_runs(raw['G_dom'])
    ref = {(r['tau'], r['seed']): r['t'] for r in raw['G_runs']}
    assert [(r['tau'], r['seed'], r['cap']) for r in out] == [(2.5e-4, 716, 3000.0), (2e-4, 716, 3000.0)]
    assert all(r['t'] is not None and r['t'] == ref[(r['tau'], 716)] for r in out) and ref[(2.5e-4, 716)] < ref[(2e-4, 716)] < 3000


def _raw_H():
    import exp_gauge_plateau as P
    raw = E.load(os.path.join(P.OUT, 'gauge_plateau_raw.npz'))
    if 'H' not in raw:
        import pytest
        pytest.skip('no Part H in the raw file')
    return P, raw


def test_part_H_fits_on_plants():
    # referee round 150 (D3es-1, -9, -10): with the barrier fixed the prefactor's exponent is recovered from a planted rate law; χ² is
    # quadratic in c (the intervals ĉ ± s.e., ± 2 s.e.); the lag's counting; the extrapolation's s.e. is the fit's prediction error
    import exp_gauge_plateau as P
    dJ, n = 6.3e-4, np.array([40, 30, 20, 15, 10, 7.0])
    for c in (0.0, 0.36, 1.0):
        rows = [dict(tau=t, n=int(k), rate=2.0 * t ** (-c) * np.exp(-dJ / t)) for t, k in zip(P.TAUS_G, n)]
        fb = P.fixed_barrier(rows, dJ)
        assert abs(fb['c'] - c) < 1e-9 and fb['chi2'] < 1e-15 and fb['dof'] == 4 and abs(fb['t4'] - 1e-4 ** c * np.exp(dJ / 1e-4) / 2) < 1e-6 * fb['t4']
        lt = np.log(P.TAUS_G); lb = n @ lt / n.sum()
        assert abs(fb['se'] - 1 / np.sqrt(n @ (lt - lb) ** 2)) < 1e-12
        for f in fb['fixed']:
            assert abs(f['chi2'] - (f['c'] - c) ** 2 / fb['se'] ** 2) < 1e-8 and f['dof'] == 5 and abs(f['z'] - (f['c'] - c) / fb['se']) < 1e-9
        assert np.allclose(fb['int1'], (c - fb['se'], c + fb['se'])) and np.allclose(fb['int2'], (c - 2 * fb['se'], c + 2 * fb['se']))
        e = P.extrap(rows, dJ, c)
        x = 1 / np.array(P.TAUS_G); X = np.c_[np.ones(6), -x]; cov = np.linalg.inv(X.T @ (n[:, None] * X))
        assert abs(e['lse4'] - np.sqrt(np.r_[1, -1e4] @ cov @ np.r_[1, -1e4])) < 1e-12 and abs(e['lseB'] - 1e4 * np.sqrt(cov[1, 1])) < 1e-12
    runs = [dict(tau=P.TAUS_G[0], t=t, cap=100.0) for t in (10.0, 30.0, None, 50.0)] + [dict(tau=t, t=None, cap=100.0) for t in P.TAUS_G[1:]]
    r = P.rows_G(runs, 20.0)[0]
    assert r['n'] == 2 and r['T'] == 0 + 10 + 80 + 30 and r['first'] == 10.0


def test_part_H_is_recomputed_from_the_raw():
    # every number of Part H from the stored runs, Part C, the minima and the exit points; t₀ = 0 gives Part G's rows and fits; the
    # fixed-barrier exponent against an independent weighted fit (np.polyfit, weights √n); the md holds the report
    P, raw = _raw_H()
    H = raw['H']
    _same(P.part_H(raw['G_dom'], raw['G_runs'], raw['C'], raw['H_min'], raw['G_exits']), H)
    g0 = P.rows_G(raw['G_runs'])
    assert all(a['n'] == b['n'] and abs(a['T'] - b['T']) < 1e-6 and a['rate'] == b['rate'] for a, b in zip(g0, raw['G']['rows']))
    _same(H['lags'][0]['fits'], raw['G']['fits'][:3])
    dJ = raw['G_dom']['dJ']
    for g in H['lags']:
        rr = P.rows_G(raw['G_runs'], g['t0'])
        t = np.array([r['tau'] for r in rr]); n = np.array([r['n'] for r in rr], float)
        y = np.log([r['rate'] for r in rr]) + dJ / t
        assert abs(-np.polyfit(np.log(t), y, 1, w=np.sqrt(n))[0] - g['fb']['c']) < 1e-10
        assert g['n'] == [sum(1 for x in raw['G_runs'] if x['tau'] == tau and x['t'] is not None and x['t'] > g['t0']) for tau in P.TAUS_G]
    assert abs(H['extC']['t4'] - P.arrhenius(raw['C']['rows'], raw['C']['dJ'], 0.5)['t4']) < 1e-9 * H['extC']['t4']
    ex = raw['G_exits']
    assert [(e['tau'], e['seed']) for e in ex] == [(r['tau'], r['seed']) for r in raw['G_runs'] if r['t'] is not None]
    Th0 = E.exact_rep(E.setup()[0], np.ones(3) / 3)
    phi = P.dom_phi(np.array([e['dev'] for e in ex]) + Th0[None], Th0)
    assert np.all(phi >= 1) and np.all(phi < 1.05)                           # the first state outside 𝒪, one step past ∂𝒪
    assert sum(r['n'] for r in H['exits']) == len(ex) and all(r['broad'] + r['narrow'] <= r['n'] for r in H['exits'])
    m = H['mins']['minima']
    assert sum(c['n'] for c in m) == P.N_DOM and abs(m[0]['dJ'] - dJ / P.TAU) < 1e-6 and all(0 < 1 - c['phi_max'] for c in m)
    md = open(os.path.join(P.OUT, 'gauge_plateau.md')).read()
    assert '\n'.join(P.report_H(H)) in md


def test_exit_points_reproduce_a_stored_exit(monkeypatch, tmp_path):
    # the stored exit point of the seed that exits first at τ = 2.5·10⁻⁴: escape_runs at its own set-up (capped at 3000 time units) leaves
    # it in its checkpoint, and exit_points reads it back; a checkpoint whose exit step disagrees with the stored exit time is refused;
    # without a checkpoint there are no exit points
    P, raw = _raw_H()
    monkeypatch.setattr(P, 'TAUS_G', (2.5e-4,)); monkeypatch.setattr(P, 'S_G', (3000.0,)); monkeypatch.setattr(P, 'SEEDS_G', (716,))
    monkeypatch.setattr(P, 'CKPT', str(tmp_path))
    out = P.escape_runs(raw['G_dom'])
    ref = [e for e in raw['G_exits'] if (e['tau'], e['seed']) == (2.5e-4, 716)][0]
    ex = P.exit_points(out, str(tmp_path))
    assert len(ex) == 1 and np.allclose(ex[0]['dev'], ref['dev'], rtol=0, atol=1e-10)
    import pytest
    with pytest.raises(AssertionError):
        P.exit_points([dict(out[0], t=out[0]['t'] + P.DT)], str(tmp_path))
    assert P.exit_points(out, None) is None and P.exit_points(out, str(tmp_path / 'none')) is None


def test_boundary_minima_reproduce_the_stored_minima(monkeypatch):
    # the first four L-BFGS starts (the generator 1700) reach the two stored minima of J on ∂𝒪, and the flow from each enters 𝒪 at once
    # and returns to the exact fit (over the full S_FLOW)
    P, raw = _raw_H()
    monkeypatch.setattr(P, 'N_DOM', 4)
    m = P.boundary_minima()
    ref = raw['H_min']['minima']
    assert m['n'] == 4 and len(m['minima']) == len(ref) == 2
    for a, b in zip(m['minima'], ref):
        assert abs(a['dJ'] - b['dJ']) < 1e-6 and np.allclose(a['u'], b['u'], atol=1e-4)
        assert abs(a['phi_max'] - b['phi_max']) < 1e-9 and abs(a['dist_end'] / b['dist_end'] - 1) < 1e-6
    assert ref[0]['u'][1, 2] > 0.9 and ref[1]['u'][1, 2] < -0.9                # the second width broadened, and narrowed


def test_cluster_minima_on_plants():
    # two minima at the same energy but different points are two clusters; nearby results of one minimum one cluster; sorted by energy
    import exp_gauge_plateau as P
    u1, u2 = np.zeros((3, 3)), np.zeros((3, 3)); u1[1, 2], u2[1, 2] = 0.99, -1.0
    vals = [6.4, 6.3, 6.3 + 1e-4, 6.3]
    U = [u2, u1, u1 + 0.01, u2]
    cl = P.cluster_minima(vals, U)
    assert [(c['n'], c['i']) for c in cl] == [(2, 1), (1, 3), (1, 0)] and [c['dJ'] for c in cl] == [6.3, 6.3, 6.4]
