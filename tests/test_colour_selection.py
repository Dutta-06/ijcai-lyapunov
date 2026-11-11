"""
Proposition D‴.43 (PARADIGM.md §5d): the colour of the measurement noise selects the kernels. Checked here: det(I + τK) against its
subset (Cauchy–Binet) expansion; the potential's gradient in the relevances D is (4/N)ω; the potential is the trace for white draws and
(1/τ)log det(I + τK) for a slow Markov noise, the step-h filter tending to it; R is linear in D for white draws and concave for coloured
ones; the support count τR − k log(4τ/N) → log(det(Φ_SᵀΦ_S)Π|w_i|); the projection keeps Φw and matches the direct formula, also when
ΦDΦᵀ is singular; the reduced flow stays on M, its P decays at the rate ω exactly, the white flow ends at weighted basis pursuit and a
slow flow, from the same start, at a sparser exact fit; the gradient flow's interpolant against a small-step descent; one step of
label-noise gradient descent against a hand-coded step with AR(1) draws.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import colour_selection as CS                              # noqa: E402


def _design(N, d, seed):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((N, d)) / np.sqrt(N), rng


def test_the_subset_expansion_and_the_gradient():
    Phi, rng = _design(4, 7, 1)
    D = rng.uniform(0.05, 1.5, 7)
    for tau in (0.3, 10.0, 1e3):
        assert abs(CS.det_subsets(Phi, D, tau) / np.linalg.det(np.eye(4) + tau * CS.gram(Phi, D)) - 1) < 1e-12
        e = 1e-6
        g = np.array([(CS.potential(Phi, D + e * np.eye(7)[i], 0, ('ou', tau)) - CS.potential(Phi, D - e * np.eye(7)[i], 0, ('ou', tau)))
                      / (2 * e) for i in range(7)])
        assert np.abs(g - (4 / 4) * CS.weights(Phi, D, 0, ('ou', tau))).max() < 1e-7 * np.abs(g).max()
    h, rho = 0.25, 0.9                                                    # at a step: ∂R/∂D = (4/N)φ_iᵀf(K)φ_i with the scheme's filter
    g = np.array([(CS.potential(Phi, D + 1e-6 * np.eye(7)[i], h, ('ar1', rho)) - CS.potential(Phi, D - 1e-6 * np.eye(7)[i], h, ('ar1', rho)))
                  / 2e-6 for i in range(7)])
    assert np.abs(g - CS.weights(Phi, D, h, ('ar1', rho))).max() < 1e-6 * np.abs(g).max()


def test_the_potentials_white_and_slow():
    Phi, rng = _design(5, 9, 2)
    D = rng.uniform(0.05, 1.5, 9)
    K = CS.gram(Phi, D)
    assert abs(CS.potential(Phi, D, 0, ('fresh',)) - np.trace(K)) < 1e-12
    assert np.allclose(CS.weights(Phi, D, 0, ('fresh',)), (Phi ** 2).sum(axis=0), rtol=1e-12)
    tc = 10.0
    assert abs(CS.potential(Phi, D, 0, ('ou', tc)) - np.log(np.linalg.det(np.eye(5) + tc * K)) / tc) < 1e-12
    h = 1e-4                                                              # finely sampled AR(1) draws → the continuum filter
    assert abs(CS.potential(Phi, D, h, ('ar1', np.exp(-h / tc))) / CS.potential(Phi, D, 0, ('ou', tc)) - 1) < 1e-3
    for _ in range(10):                                                   # R along a segment in D: linear (white), concave (coloured)
        D1, D2 = rng.uniform(0.05, 1.5, 9), rng.uniform(0.05, 1.5, 9)
        sec = lambda corr: sum(c * CS.potential(Phi, D1 + t * (D2 - D1), 0, corr) for c, t in ((1, 0.0), (-2, 0.5), (1, 1.0)))
        assert abs(sec(('fresh',))) < 1e-12 and sec(('ou', tc)) < 0


def test_the_support_count():
    Phi, rng = _design(10, 40, 3)
    w = np.zeros(40); S = [3, 17, 25]; w[S] = (1.0, -0.6, 0.8)
    ref = np.log(np.linalg.det(Phi[:, S].T @ Phi[:, S]) * np.prod(np.abs(w[S])))
    vals = [tc * CS.potential(Phi, np.abs(w), 0, ('ou', tc)) - 3 * np.log(4 * tc / 10) for tc in (1e2, 1e4, 1e6)]
    assert abs(vals[2] - ref) < 1e-4 and abs(vals[2] - ref) < abs(vals[1] - ref) < abs(vals[0] - ref)


def test_the_projection():
    Phi, rng = _design(6, 15, 4)
    D = rng.uniform(0.1, 1.0, 15); q = rng.standard_normal(15)
    direct = q - D * (Phi.T @ np.linalg.solve((Phi * D) @ Phi.T, Phi @ q))
    p = CS.project(Phi, D, q)
    assert np.abs(p - direct).max() < 1e-12 and np.abs(Phi @ p).max() < 1e-12
    assert np.abs(CS.project(Phi, D, p) - p).max() < 1e-12                # a projection
    D[3:] = 1e-200                                                        # three atoms alive: ΦDΦᵀ singular to working precision;
    q = rng.uniform(0.5, 2, 15) * D * rng.uniform(-1, 1, 15)              # in the flow q = ω⊙w with |w_i| ≤ D_i
    p = CS.project(Phi, D, q)
    assert np.all(np.isfinite(p)) and np.abs(Phi @ p).max() < 1e-12 * np.abs(q).max()
    assert np.abs(p[3:]).max() < 1e-150                                   # the dead atoms do not move


def test_the_reduced_flow():
    Phi, rng = _design(10, 40, 0)
    ws = np.zeros(40); idx = rng.choice(40, 4, replace=False); ws[idx] = rng.choice([-1, 1], 4) * rng.uniform(0.5, 1.5, 4)
    y = Phi @ ws
    w0, P0 = CS.flow_interpolant(Phi, y, np.full(40, 0.5), np.full(40, 0.5))
    assert np.abs(Phi @ w0 - y).max() < 1e-12 and np.allclose(P0, 0.25)
    s, W, P = CS.reduced_flow(Phi, w0, P0, 5.0, 0, ('ou', 10.0), t_eval=np.linspace(0, 5.0, 6))
    assert np.abs(Phi @ W - y[:, None]).max() < 1e-7                      # the flow stays on M
    z = np.r_[w0, np.log(P0)]
    om = CS.weights(Phi, np.sqrt(w0 ** 2 + 4 * P0 ** 2), 0, ('ou', 10.0))
    assert np.allclose(CS.rhs(0.0, z, Phi, 0, ('ou', 10.0))[40:], -om)    # d log P/ds = −ω
    D1 = np.full(40, 1e-300); D1[idx] = np.abs(ws[idx])                   # at a 4-sparse fit K has 6 zero eigenvalues: the
    assert np.all(np.isfinite(CS.weights(Phi, D1, 0.5, ('ar1', 0.9))))    # finite-step weights stay finite there
    with np.errstate(over='ignore'):                                      # an overshooting trial stage returns NaN (rejected)
        assert np.all(np.isnan(CS.rhs(0.0, np.r_[w0, np.full(40, 1e3)], Phi, 0.5, ('fresh',))))
    c = (Phi ** 2).sum(axis=0)
    wbp = CS.basis_pursuit(Phi, y, c)
    _, Ww, _ = CS.reduced_flow(Phi, w0, P0, 2000.0, 0, ('fresh',))
    assert np.abs(Ww[:, -1] - wbp).max() < 1e-5                          # white: weighted basis pursuit (here not the truth)
    assert np.abs(wbp - ws).max() > 0.5
    _, Wt, _ = CS.reduced_flow(Phi, w0, P0, 400.0 + 60 * np.sqrt(101.0), 0, ('ou', 100.0))
    assert np.abs(Wt[:, -1] - ws).max() < 1e-6                            # slow: the 4-sparse truth
    assert (np.abs(Wt[:, -1]) > 1e-3).sum() == 4 < (np.abs(wbp) > 1e-3).sum()


def test_the_flow_interpolant_and_one_noisy_step():
    Phi, rng = _design(8, 20, 5)
    y = Phi @ rng.standard_normal(20)
    u, v = rng.uniform(0.2, 1, 20), rng.uniform(0.2, 1, 20)
    wf, Pf = CS.flow_interpolant(Phi, y, u, v)
    wg, _ = CS.gd_interpolant(Phi, y, u=u, v=v, h=0.002)
    assert np.abs(Phi @ wf - y).max() < 1e-12 and np.allclose(Pf, u * v) and np.abs(wf - wg).max() < 2e-4
    us, vs = np.full(20, 1e-150), np.full(20, 1e-150); us[:3], vs[:3] = u[:3], v[:3]   # three atoms alive, y outside their span:
    W, ok = CS.end_interpolants(Phi, y, np.array([u, us]), np.array([v, vs]))           # the interpolant is at infinity, flagged
    assert ok.tolist() == [True, False] and np.allclose(W[0], wf) and np.all(np.isnan(W[1]))
    ys = Phi[:, :3] @ np.array([0.3, -0.2, 0.1])                                         # within their span it converges, singular Hm
    W, ok = CS.end_interpolants(Phi, ys, us[None, :], vs[None, :])
    assert ok[0] and np.abs(Phi @ W[0] - ys).max() < 1e-9 and np.abs(W[0][3:]).max() < 1e-100

    class Rec:
        def __init__(self):
            self.r, self.out = np.random.default_rng(9), []

        def standard_normal(self, shape):
            self.out.append(self.r.standard_normal(shape)); return self.out[-1]
    h, sig, rho, rec = 0.1, 0.3, 0.8, Rec()
    u1, v1, _ = CS.label_noise_gd(Phi, y, u, v, h, sig, rho, 2, rec, 3)
    xi = sig * rec.out[0]; uu, vv = np.tile(u, (3, 1)), np.tile(v, (3, 1))
    for k in range(2):
        g = ((uu * uu - vv * vv) @ Phi.T - y[None, :] - xi) @ Phi / 8
        uu, vv = uu - 2 * h * g * uu, vv + 2 * h * g * vv
        xi = rho * xi + np.sqrt(1 - rho * rho) * sig * rec.out[k + 1]
    assert np.allclose(u1, uu, rtol=0, atol=1e-14) and np.allclose(v1, vv, rtol=0, atol=1e-14)


def test_the_kernel_dictionary():
    Phi = CS.rbf_dictionary(np.linspace(0, 1, 15), np.linspace(0, 1, 40), 0.05)
    assert Phi.shape == (15, 40) and np.allclose(np.linalg.norm(Phi, axis=0), 1.0)


def test_the_flow_is_the_gradient_flow_of_R_in_its_metric():
    """The reduced flow's w- and P-drifts against the gradient of R in the metric diag(1/(4D), 1/D) (one common factor): a doubled or
    halved w-drift fails."""
    Phi, rng = _design(5, 9, 21)
    w0 = rng.standard_normal(9); y = Phi @ w0
    P = rng.uniform(0.05, 0.6, 9)
    for corr in (('fresh',), ('ou', 30.0)):
        z = np.r_[w0, np.log(P)]
        f = CS.rhs(0.0, z, Phi, 0, corr)
        dw, dP = f[:9], f[9:] * P                                              # rhs is in (w, log P)
        R = lambda w, P: CS.potential(Phi, np.sqrt(w * w + 4 * P * P), 0, corr)
        e = 1e-6
        gw = np.array([(R(w0 + e * np.eye(9)[i], P) - R(w0 - e * np.eye(9)[i], P)) / (2 * e) for i in range(9)])
        gP = np.array([(R(w0, P + e * np.eye(9)[i]) - R(w0, P - e * np.eye(9)[i])) / (2 * e) for i in range(9)])
        D = np.sqrt(w0 ** 2 + 4 * P ** 2)
        pw = -CS.project(Phi, D, 4 * D * gw)
        pP = -D * gP
        c = (dP @ pP) / (pP @ pP)
        assert c > 0
        assert np.allclose(dP, c * pP, rtol=1e-5, atol=1e-10) and np.allclose(dw, c * pw, rtol=1e-5, atol=1e-10)
        assert np.abs(Phi @ dw).max() < 1e-10


def test_uv_round_trip_and_the_singular_projection():
    rng = np.random.default_rng(3)
    w, P = rng.standard_normal(12), rng.uniform(1e-6, 2.0, 12)
    u, v = CS.uv_of(w, P)
    assert np.all(u > 0) and np.all(v > 0)
    assert np.allclose(u * u - v * v, w, atol=1e-12) and np.allclose(u * v, P, rtol=1e-12)
    Phi, rng = _design(6, 10, 4)
    D = np.r_[rng.uniform(0.5, 1.0, 3), np.full(7, 1e-300)]                   # three atoms alive, fewer than N: ΦDΦᵀ singular
    q = D * rng.standard_normal(10)                                             # the flow's form: D times a gradient (dead atoms fixed)
    pq = CS.project(Phi, D, q)
    assert np.all(np.isfinite(pq))
    assert np.abs(Phi @ pq).max() < 1e-8 * np.abs(q).max()
    assert np.allclose(CS.project(Phi, D, pq), pq, atol=1e-10)                 # a projection


def test_the_experiment_helpers_and_the_new_references():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_colour_selection as X
    ws = np.array([1.0, 0.0, -0.5])
    assert X.ok(ws + 0.9e-4, ws) and not X.ok(ws + 1.1e-4, ws)
    Phi = np.c_[np.eye(3), np.array([1.0, 1.0, 0.0]) / np.sqrt(2)]            # three orthonormal atoms and their bisector
    w = np.array([1.0, -1.0, 0.0, 0.0])
    ratio, k = X.certificate(Phi, w, 0.0)
    om = (Phi ** 2).sum(axis=0)
    c = Phi[:, 3] @ np.linalg.pinv(Phi[:, :2]).T @ (om[:2] * np.sign(w[:2]))  # = (1 − 1)/√2 = 0
    assert k == 2 and abs(ratio - abs(c) / om[3]) < 1e-12 and ratio < 1e-12
    w2 = np.array([1.0, 1.0, 0.0, 0.0])                                         # same signs: c = √2 > ω = 1, not certified
    assert abs(X.certificate(Phi, w2, 0.0)[0] - np.sqrt(2)) < 1e-12
    Phi, rng = _design(8, 20, 9)
    ws = np.zeros(20); ws[[3, 11]] = [1.0, -0.7]
    assert np.abs(CS.sbl_em(Phi, Phi @ ws) - ws).max() < 1e-6
    w, P, s, conv = CS.flow_to_convergence(Phi, *X.start(Phi, Phi @ ws, 0.5), 50.0, 0, ('fresh',), S_max=5e3)
    assert conv and P.max() < 1e-8 * np.abs(w).max()
    u0, v0 = CS.uv_of(*X.start(Phi, Phi @ ws, 0.5))
    u, v = CS.label_noise_heun(Phi, Phi @ ws, u0, v0, 0.2, 0.0, 0.0, 50, np.random.default_rng(1), 2)   # noiseless: stays on M
    assert np.abs((u * u - v * v) @ Phi.T - Phi @ ws).max() < 1e-8


def test_part_K_statistics_and_the_vertex_limit():
    """Part K's helpers: the exact binomial interval and McNemar's exact test against closed forms; at an exact fit on N atoms the slow
    certificate's ratio tends, as τ_c → ∞, to the τ_c-free vertex limit."""
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_colour_selection as X
    lo, hi = X.clopper_pearson(0, 40)
    assert lo == 0.0 and abs(hi - (1 - 0.025 ** (1 / 40))) < 1e-12           # x = 0: the upper end solves (1 − p)ⁿ = α/2
    lo, hi = X.clopper_pearson(40, 40)
    assert hi == 1.0 and abs(lo - 0.025 ** (1 / 40)) < 1e-12
    lo, hi = X.clopper_pearson(15, 40)
    assert lo < 15 / 40 < hi
    assert abs(X.mcnemar_exact(15, 0) - 2 * 0.5 ** 15) < 1e-15                # one-sided tail doubled
    assert abs(X.mcnemar_exact(1, 6) - 2 * (1 + 7) / 2 ** 7) < 1e-12
    assert X.mcnemar_exact(0, 0) == 1.0
    Phi, rng = _design(5, 12, 21)
    w = np.zeros(12); S = rng.choice(12, 5, replace=False); w[S] = rng.choice([-1, 1], 5) * rng.uniform(0.5, 1.5, 5)
    lim = X.vertex_limit(Phi, w)
    r6, r8 = X.certificate(Phi, w, 1e6)[0], X.certificate(Phi, w, 1e8)[0]
    assert abs(r8 - lim) < 1e-5 * max(lim, 1.0) and abs(r8 - lim) < abs(r6 - lim) + 1e-12
    w2 = w.copy(); w2[S[0]] = 0.0
    assert np.isnan(X.vertex_limit(Phi, w2))                                   # fewer than N atoms: no τ_c-free limit (the ratio → 0)


RAW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'colour_selection_raw.npz')


def _raw():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_colour_selection as X
    z = np.load(RAW, allow_pickle=True)
    return X, {k: z[k][0] for k in z.files}


def test_part_K1_is_recomputed_from_parts_C_and_G():
    """K1 (mutation testing after round 142): the stored rates, intervals, discordant counts and McNemar p-values are `k_stats` of the
    stored Parts C and G; the sparse Bayesian learning column only where G has every instance of the cell."""
    X, raw = _raw()
    st = X.k_stats(raw['C'], raw['G'])
    S = raw['K']['stats']
    assert len(st['rates']) == len(S['rates']) and len(st['pairs']) == len(S['pairs'])
    for a, b in zip(st['rates'], S['rates']):
        assert {k: a[k] for k in ('N', 'd', 'k', 'method', 'x', 'n')} == {k: b[k] for k in ('N', 'd', 'k', 'method', 'x', 'n')}
        assert np.allclose(a['ci'], b['ci'], rtol=1e-12, atol=0)
    for a, b in zip(st['pairs'], S['pairs']):
        assert {k: a[k] for k in ('N', 'd', 'k', 'other', 'n', 'slow_only', 'other_only')} == \
            {k: b[k] for k in ('N', 'd', 'k', 'other', 'n', 'slow_only', 'other_only')}
        assert np.isclose(a['p'], b['p'], rtol=1e-12, atol=0)
    assert {r['other'] for r in S['pairs']} == {'bp', 'sparse Bayesian learning (EM)', 'slow 1000'}
    G = [r for r in raw['G'] if (r['k'], r['i']) != (3, 0)]                 # one instance of k = 3 missing: no SBL rate for k = 3
    st = X.k_stats(raw['C'], G)
    assert not [r for r in st['rates'] if r['method'].startswith('sparse') and r['k'] == 3]
    assert [r for r in st['rates'] if r['method'].startswith('sparse') and r['k'] == 4]


def test_part_K2_k_job_reproduces_stored_ends():
    """K2 (mutation testing after round 142): the stored failures are Part C's 10 × 40 failures of the slow flow (τ_c = 100), and `k_job`
    rerun on two of them gives the stored end: support, basic, the certificate for each τ_c, the vertex limit, R(end) − R(truth)."""
    X, raw = _raw()
    cp = raw['K']['competitors']
    fails = sorted((r['k'], r['i']) for r in raw['C'] if (r['N'], r['d']) == (10, 40) and r['k'] in X.KS_K and not r['out']['slow 100'][0])
    assert [(r['k'], r['i']) for r in cp] == fails
    for s in (cp[0], cp[-1]):
        r = X.k_job((s['k'], s['i']))
        assert (r['conv'], r['recovers'], r['nnz'], r['basic']) == (s['conv'], s['recovers'], s['nnz'], s['basic'])
        assert np.allclose(r['w'], s['w'], rtol=1e-8, atol=1e-12) and np.isclose(r['s'], s['s'], rtol=1e-12, atol=0)
        assert sorted(r['cert']) == sorted(s['cert']) == sorted(X.TAUS_KC)
        assert all(np.isclose(r['cert'][tc], s['cert'][tc], rtol=1e-8, atol=0) for tc in X.TAUS_KC)
        assert np.isclose(r['limit'], s['limit'], rtol=1e-8, atol=0)
        assert sorted(r['dR']) == sorted(X.TAUS_KR) and all(np.isclose(r['dR'][tc], s['dR'][tc], rtol=1e-8, atol=1e-14) for tc in X.TAUS_KR)


def test_report_K_counts_planted_ends():
    """K2's table on the raw with planted ends (k = 3): one end not converged, one basic end on nine atoms — counted as a failure, not as
    'converged, failing', and outside the N-atom ends' fractions."""
    import copy
    X, raw = _raw()
    raw = copy.deepcopy(raw)
    cp = [r for r in raw['K']['competitors'] if r['k'] == 3]
    cp[0]['conv'] = False
    cp[1]['nnz'] = 9
    md = X.report(*[raw.get(k) for k in 'ABCDEFGHIJK'])
    line = next(ln for ln in md.splitlines() if ln.startswith('| 3 | 10 | '))
    assert line.startswith('| 3 | 10 | 9 | 10 (9–10) | 0.90 | ')


def test_k_job_on_a_planted_end():
    """`k_job` with the flow's end planted (`flow_to_convergence` replaced; mutation testing after round 142): a stored ten-atom end
    plus an eleventh coefficient at 10⁻⁷ of the largest (inside the support threshold 10⁻⁹) and P = 0.05|w| — eleven atoms, so not basic,
    no certificate and no vertex limit; R(end) uses the relevances D = u² + v² = (w² + 4P²)^{1/2}."""
    X, raw = _raw()
    s = raw['K']['competitors'][0]
    k, i = s['k'], s['i']
    w = np.asarray(s['w'], float).copy()
    assert s['nnz'] == 10 and s['basic']
    off = np.flatnonzero(np.abs(w) <= 1e-9 * np.abs(w).max())
    w[off[0]] = 1e-7 * np.abs(w).max()
    P = 0.05 * np.abs(w)
    saved = X.CS.flow_to_convergence
    try:
        X.CS.flow_to_convergence = lambda *a, **kw: (w.copy(), P.copy(), 123.0, True)
        r = X.k_job((k, i))
    finally:
        X.CS.flow_to_convergence = saved
    assert r['nnz'] == 11 and not r['basic'] and r['cert'] == {} and np.isnan(r['limit']) and r['s'] == 123.0 and r['conv']
    Phi, ws = X.instance(10, 40, k, X.seed_of(10, k, i))
    De = np.sqrt(w ** 2 + 4 * P ** 2)
    for tc in X.TAUS_KR:
        exp = CS.potential(Phi, De, 0, ('ou', tc)) - CS.potential(Phi, np.abs(ws), 0, ('ou', tc))
        assert np.isclose(r['dR'][tc], exp, rtol=1e-12, atol=0), tc
