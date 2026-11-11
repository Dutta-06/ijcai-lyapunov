"""Tests of experiments/exp_mirror_leak.py, Part E (round 141): the 1760 patterns of D(2) against the weighted LASSO itself."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments'))
from src.core import mirror_free as MF
import exp_mirror_leak as X

RAW = os.path.join(ROOT, 'results', 'chaos', 'mirror_leak_raw.npz')
MD = os.path.join(ROOT, 'results', 'chaos', 'mirror_leak.md')


def _orbit():
    L = MF.FreeLaw(12.0, 1.0, 0.1, 1.0)
    H, _, _ = X.orbit_gram(np.linspace(-1.4, 1.4, X.D12), np.full(X.D12, X.WID12))
    b0 = 1e-6
    lam1 = np.pi * L.A * np.sqrt(abs(2 * L.c - b0) / b0) / (4 * L.G * L.c)
    q = L.Q(np.zeros(X.D12), np.full(X.D12, b0))
    return L, H, b0, q, X.MU_IC / lam1


def test_lasso_exact_solves_the_kkt_conditions():
    rng = np.random.default_rng(5)
    for rep in range(5):
        A = rng.standard_normal((30, 8)); H = A.T @ A / 30
        ws = np.zeros(8); ws[rng.choice(8, 3, replace=False)] = rng.choice([-1, 1], 3) * rng.uniform(0.1, 0.4, 3)
        for mu, ke in ((1e-3, 1.0), (0.02, 2.0)):
            w, kkt, slack = X.lasso_exact(H, ws, mu, ke)
            c = H @ (ws - w) / (ke * mu)
            on = w != 0
            assert 0 <= kkt < 1e-10                              # the residual is a violation: off S only |c| − 1 above 0 counts
            assert np.allclose(c[on], np.sign(w[on]), atol=1e-10)
            assert np.all(np.abs(c[~on]) <= 1 + 1e-10)
            assert np.allclose(slack[~on], 1 - np.abs(c[~on])) and np.all(np.isnan(slack[on]))
            # the objective is not lowered by the repository's coordinate descent run much longer
            f = lambda v: 0.5 * (v - ws) @ H @ (v - ws) / ke + mu * np.abs(v).sum()
            assert f(w) <= f(MF.lasso(H, ws, np.full(8, mu), ke, it=200000, tol=1e-16)) + 1e-14


def _planted_lasso(rng, n=8, mu=1e-3, ke=2.0, tiny=1e-7):
    """A LASSO whose exact solution is planted: w with support S and signs s_S, one coefficient tiny, and w* = w + k_eμH⁻¹s with
    |s_j| ≤ 0.9 off S (so w is the unique minimiser and its tiny coordinate's KKT multiplier is ≈ 1 + H_jjw_j/(k_eμ) once dropped)."""
    A = rng.standard_normal((40, n)); H = A.T @ A / 40
    w = np.zeros(n); S = rng.choice(n, 3, replace=False)
    w[S] = rng.choice([-1, 1], 3) * np.array([0.3, 0.2, tiny])
    s = rng.uniform(-0.9, 0.9, n); s[S] = np.sign(w[S])
    return H, w + ke * mu * np.linalg.solve(H, s), w, S, mu, ke


@pytest.mark.parametrize('start', ['zeros', 'tiny dropped', 'sign flipped'])
def test_lasso_exact_corrects_any_warm_start(start):
    """The active-set corrections of `lasso_exact` (mutation testing after round 141), with coordinate descent's warm start planted: from
    zero every atom enters as the largest KKT violator with its multiplier's sign; with the tiny coefficient dropped it re-enters at a
    violation ≈ 10⁻⁴ (well below 10⁻³); with a coefficient's sign flipped the flip leaves the support and the atom re-enters. Each ends at
    the planted minimiser with the KKT conditions."""
    rng = np.random.default_rng(11)
    H, ws, w0, S, mu, ke = _planted_lasso(rng)
    j_tiny, j_big = S[2], S[0]
    warm = dict(zeros=np.zeros(len(ws)), **{'tiny dropped': np.where(np.arange(len(ws)) == j_tiny, 0.0, w0),
                                              'sign flipped': np.where(np.arange(len(ws)) == j_big, -w0, w0)})[start]
    c_drop = (H @ (ws - np.where(np.arange(len(ws)) == j_tiny, 0.0, w0)) / (ke * mu))[j_tiny]
    assert 1 + 1e-6 < abs(c_drop) < 1 + 1e-3
    saved = X.MF.lasso
    try:
        X.MF.lasso = lambda *a, **k: warm.copy()
        w, kkt, slack = X.lasso_exact(H, ws, mu, ke)
    finally:
        X.MF.lasso = saved
    assert 0 <= kkt < 1e-10
    assert np.array_equal(np.flatnonzero(w), np.sort(S)) and np.allclose(w, w0, rtol=0, atol=1e-12)


def test_tangent_profile_is_the_core_at_small_eta():
    """At β⁰ → 0 a coordinate whose η sits at (1 − μ)λ₁ has w on the tangent profile (Derived statement 4.4)."""
    L = MF.FreeLaw(12.0, 1.0, 0.1, 1.0)
    b0 = 1e-8
    q = L.Q(np.zeros(1), np.full(1, b0))
    lam1 = np.pi * L.A * np.sqrt(abs(2 * L.c - b0) / b0) / (4 * L.G * L.c)
    for mu in (0.05, 0.3, 0.7):
        al = L.alpha_of_eta(np.array([(1 - mu) * lam1]), q)
        w = float(L.w(al, q)[0])
        tp = float(X.tangent_profile(L, b0, mu, 1.0))
        assert abs(w / tp - 1) < 0.01, (mu, w, tp)


@pytest.mark.parametrize('S, sg, mag', [((1, 5, 8), (1.0, -1.0, 1.0), (0.2, 0.25, 0.18)),
                                        ((2, 6, 9), (1.0, -1.0, 1.0), (0.25, 0.2, 0.3)),
                                        ((0, 1, 2), (1.0, 1.0, 1.0), (0.2, 0.2, 0.2))])
def test_pattern_record_lasso_support_is_the_prediction(S, sg, mag):
    L, H, b0, q, sig = _orbit()
    S_ = np.array(S); sg_ = np.array(sg)
    wt = np.zeros(X.D12); wt[S_] = sg_ * np.array(mag)
    r = X.pattern_record(L, H, wt, S_, sg_, q, sig, X.MU_IC, b0)
    assert r['kkt'] < 1e-10 and r['res'] < 1e-10
    assert r['l0_exact'] == r['pred']                    # the LASSO's exact support is the prediction (the LASSO's theory)
    assert r['rel'] < 0.01                               # the rich equilibrium is the LASSO to O(10⁻³) at β⁰ = 10⁻⁶


def test_raw_part_E_is_consistent_with_part_D_and_regenerates_its_md_section():
    if not os.path.exists(RAW):
        pytest.skip('no raw file')
    z = np.load(RAW, allow_pickle=True)
    raw = {k: z[k][0] for k in z.files}
    assert {'A', 'B', 'C', 'D', 'E'} <= set(raw)
    P, ic = raw['E']['pats'], raw['D']['ic']
    assert len(P) == ic['n'] == 1760
    assert sum(p['eq_exact'] for p in P) == ic['exact']            # the same equilibria as D(2)
    assert sum(p['pred'] for p in P) == ic['pred']
    assert sum(p['eq_exact'] == p['pred'] for p in P) == ic['agree']
    assert all(p['l0_exact'] == p['pred'] for p in P)                # the exact LASSO agrees with the prediction everywhere
    assert max(p['kkt'] for p in P) < 1e-10
    sec = '\n'.join(X.report_E(raw['E']))
    assert sec in open(MD, encoding='utf-8').read()


def test_part_E_replays_the_stored_patterns():
    """`part_E` (mutation testing after round 141) with `pattern_record` recorded: the 1760 patterns in the stored order with the stored
    supports and signs, the magnitudes in [0.15, 0.3] drawn in that order (seed 23), at μ = σλ₁ = MU_IC; the stored β⁰ and core; and seven
    patterns recomputed in full against the raw — the first, the first where the prediction and the equilibrium disagree on each side,
    the first where both fail, the first where the equilibrium's support differs from the LASSO's exact one, and two where the LASSO's
    threshold support differs from its exact one, one its closest call to the threshold (its smallest coefficient 0.998·10⁻³ of its
    largest)."""
    if not os.path.exists(RAW):
        pytest.skip('no raw file')
    z = np.load(RAW, allow_pickle=True)
    Es = z['E'][0]
    P = Es['pats']
    ltd = [i for i, p in enumerate(P) if not p['lt_l0']]
    pick = {0, next(i for i, p in enumerate(P) if p['pred'] and not p['eq_exact']),
            next(i for i, p in enumerate(P) if not p['pred'] and p['eq_exact']),
            next(i for i, p in enumerate(P) if not p['pred'] and not p['eq_exact']),
            next(i for i, p in enumerate(P) if not p['eq_l0']), max(ltd, key=lambda i: P[i]['l_small'])}
    pick.add(next(i for i in range(len(P)) if i not in pick and not P[i]['lt_l0']))
    calls = []
    real = X.pattern_record

    def rec(L, H, wt, S_, sg_, q, sig, mu, b0, thr=1e-3):
        calls.append((wt.copy(), list(S_), list(sg_), sig, mu, b0))
        return real(L, H, wt, S_, sg_, q, sig, mu, b0, thr) if len(calls) - 1 in pick else {}
    try:
        X.pattern_record = rec
        out = X.part_E()
    finally:
        X.pattern_record = real
    assert len(calls) == len(P) == 1760
    assert all(list(c[1]) == list(p['S']) and list(c[2]) == list(p['sg']) for c, p in zip(calls, P))
    L = MF.FreeLaw(12.0, 1.0, 0.1, 1.0)
    lam1 = np.pi * L.A * np.sqrt(abs(2 * L.c - 1e-6) / 1e-6) / (4 * L.G * L.c)
    rng = np.random.default_rng(23)
    for c in calls:
        assert np.array_equal(c[0][c[1]], np.array(c[2]) * rng.uniform(0.15, 0.3, 3)) and np.count_nonzero(c[0]) == 3
        assert c[3] == X.MU_IC / lam1 and c[4] == X.MU_IC == Es['mu'] and c[5] == Es['b0'] == 1e-6
    assert out['b0'] == Es['b0'] and out['mu'] == Es['mu'] and np.isclose(out['core'], Es['core'], rtol=1e-14, atol=0)
    for i in sorted(pick):
        a, b = out['pats'][i], P[i]
        for k, v in b.items():
            if k in ('S', 'sg'):
                continue
            if isinstance(v, (bool, np.bool_, int, np.integer)):
                assert a[k] == v, (i, k)
            elif k == 'kkt':
                assert 0 <= a[k] < 1e-10, (i, k)
            elif k == 'res':
                assert a[k] < 1e-10, (i, k)
            else:
                assert np.isclose(a[k], v, rtol=1e-6, atol=0, equal_nan=True), (i, k, a[k], v)


def test_pattern_record_prediction_and_threshold_support_on_planted_inputs():
    """Two planted cases (mutation testing after round 141; on the stored patterns neither decides, see below): (1) a truth with one
    magnitude 10⁻⁴ ≪ μ: the irrepresentable quantity is that of a stored pattern (< 1) but the shifted w_S − μH_SS⁻¹s flips that sign, so
    the prediction fails — and the exact LASSO drops the atom; (2) the equilibrium planted at half the exact LASSO (`leak_equilibrium`
    replaced): the LASSO's threshold support is 10⁻³ of the LASSO's own largest coefficient, so a coefficient at 0.998·10⁻³ of it stays
    out (with the equilibrium's largest it would be in). On all 1760 stored patterns the sign condition never fails where ic < 1 and no
    LASSO coefficient lies within the equilibrium–LASSO distance of the threshold, so neither rule changes a stored number."""
    L, H, b0, q, sig = _orbit()
    S_, sg_ = np.array([2, 6, 9]), np.array([1.0, -1.0, 1.0])
    wt = np.zeros(X.D12); wt[S_] = sg_ * np.array([0.25, 0.2, 1e-4])
    assert MF.irrepresentable(H, S_, sg_) < 1
    r = X.pattern_record(L, H, wt, S_, sg_, q, sig, X.MU_IC, b0)
    assert not r['pred'] and r['ic'] < 1 and not r['l0_exact'] and r['l0_exact'] == r['pred']
    z = np.load(RAW, allow_pickle=True)
    P = z['E'][0]['pats']
    i = max((i for i, p in enumerate(P) if not p['lt_l0']), key=lambda i: P[i]['l_small'])
    assert 0.99e-3 < P[i]['l_small'] < 1e-3
    import itertools
    rng = np.random.default_rng(23)
    for j, (S, sg) in enumerate((S, sg) for S in itertools.combinations(range(X.D12), 3) for sg in itertools.product((1.0, -1.0), repeat=3)):
        mag = rng.uniform(0.15, 0.3, 3)
        if j == i:
            break
    assert list(S) == list(P[i]['S']) and list(sg) == list(P[i]['sg'])
    wt = np.zeros(X.D12); wt[list(S)] = np.array(sg) * mag
    wl = X.lasso_exact(H, wt, X.MU_IC, 1.0)[0]
    L.leak_equilibrium = lambda *a, **k: (0.5 * wl, None, 0, 0.0)
    r = X.pattern_record(L, H, wt, np.array(S), np.array(sg), q, sig, X.MU_IC, b0)
    assert np.isclose(r['l_small'], P[i]['l_small'], rtol=1e-9, atol=0)
    assert not r['lt_l0'] and r['eq_lt'] and r['n_l0'] == P[i]['n_l0']
