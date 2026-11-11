"""Proposition D‴.59: the transit with moving widths (`src/core/moving_width.py`)."""
import os
import sys

import numpy as np

from src.core import moving_width as MW
from src.core import rich_transit_nd as RN

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_moving_width as MWE  # noqa: E402


def test_bumps_cw_derivatives_and_grams():
    """∂_cφ, ∂_σφ and ⟨r∇²_θφ⟩ against central differences of φ; the Gram blocks against their sums."""
    rng = np.random.default_rng(0)
    K = MW.BumpsCW(np.linspace(-1.5, 1.5, 201), 0.12)
    th = np.c_[rng.uniform(-0.5, 0.5, 4), rng.uniform(0.05, 0.3, 4)]
    r = rng.standard_normal(K.N)
    h = 1e-6
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.phi(th + e) - K.phi(th - e)) / (2 * h)
        assert np.abs(fd - K.dphi(th)[:, :, k]).max() < 1e-6 * np.abs(fd).max()
    g, hh, Q = K.corr(th, r, second=True)
    assert np.allclose(g, K.phi(th).T @ r / K.N, rtol=1e-12, atol=0)
    assert np.allclose(hh, np.einsum('s,sjk->jk', r, K.dphi(th)) / K.N, rtol=1e-12, atol=1e-15)
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.corr(th + e, r)[1] - K.corr(th - e, r)[1]) / (2 * h)
        assert np.abs(fd - Q[:, :, k]).max() < 1e-5 * np.abs(Q).max()
    P, D, H = K.grams(th)
    dP = K.dphi(th)
    assert np.allclose(P, K.phi(th).T @ K.phi(th) / K.N) and np.allclose(D[1, 2], K.phi(th)[:, 1] @ dP[:, 2, :] / K.N)
    assert np.allclose(H[0, 3], dP[:, 0, :].T @ dP[:, 3, :] / K.N)


def test_the_one_kernel_reduced_energy():
    """`overlap` against quadrature; the reduced energy vanishes at (c*, σ*) only and its closed-form gradient matches central differences."""
    x = np.linspace(-12, 12, 240001)
    dx = x[1] - x[0]
    q = np.exp(-(x - 0.3) ** 2 / (2 * 0.2 ** 2)) @ np.exp(-(x + 0.1) ** 2 / (2 * 0.45 ** 2)) * dx
    assert abs(MW.overlap(0.3, 0.2, -0.1, 0.45) / q - 1) < 1e-10
    A, cs, ss = 0.8, 0.03, 0.2
    assert abs(MW.reduced_energy(cs, ss, A, cs, ss)) < 1e-15
    rng = np.random.default_rng(1)
    for _ in range(10):
        c, s = rng.normal(0, 0.3), ss * np.exp(rng.normal(0, 0.7))
        assert MW.reduced_energy(c, s, A, cs, ss) > 0
        h = 1e-6
        fd = np.array([(MW.reduced_energy(c + h, s, A, cs, ss) - MW.reduced_energy(c - h, s, A, cs, ss)) / (2 * h),
                       (MW.reduced_energy(c, s + h, A, cs, ss) - MW.reduced_energy(c, s - h, A, cs, ss)) / (2 * h)])
        g = MW.one_kernel_grad(c, s, A, cs, ss)
        assert np.abs(fd - g).max() < 1e-6 * max(np.abs(g).max(), 1e-12)


def test_a_narrow_feature_is_formed_by_one_kernel_with_moving_widths():
    """At a quarter of the kernels' width the sliding path with moving widths enters one kernel and relaxes it to the feature exactly
    (centre, width, amplitude), within a small budget; with the widths fixed the same budget runs out."""
    import exp_moving_width as M
    L = M.law()
    ss = 0.25 * M.S0
    f = M.feature_sum([M.C_A], [ss], [M.A_A])
    Rm = RN.sliding_iss_nd(L, MW.BumpsCW(M.XS, M.S0), f, M.theta0(), np.ones(M.C0.size), max_events=40, exit_kkt=True, max_rhs=50_000)
    m = M.match(Rm, [M.C_A], [ss], [M.A_A])
    assert m['n'] == 1 and m['exact'] and not m['oob'] and m['entries'] == 1
    Rf = RN.sliding_iss_nd(L, RN.BumpsND(M.XS[:, None], M.S0), f, M.C0[:, None], np.ones(M.C0.size), max_events=40,
                           max_entries=M.MAX_ENT_FIXED, exit_kkt=True, max_rhs=50_000)
    assert Rf['out_of_budget'] and not M.match(Rf, [M.C_A], [ss], [M.A_A])['exact']


def test_part_b_gradient_check_is_well_conditioned():
    """Part B checks the closed-form gradient against differences of the overlap term alone: in the far tail (c = 6.25σ*, where the
    gradient is ~10⁻¹⁵ and differences of J are swamped by its constant) the relative error stays at the differences' own level."""
    import exp_moving_width as M
    B = M.part_B()
    assert max(r['grad_err'] for r in B) < 1e-5
    assert all(len(r['minima']) == 1 and abs(r['minima'][0][0]) < 1e-12 for r in B)


def test_outcome_classes_and_pairs():
    """`pairs_of` finds kernels at one (c, σ) within PAIR_TOL of their width; `outcome` sorts a path's end into exact (idle entrants
    allowed), a coinciding pair, stopped by the budget, or not fitted."""
    import exp_moving_width as M
    th = [[0.35001, 0.2731707], [0.35097, 0.2731707], [-0.6, 0.08], [0.2, 0.1]]
    assert M.pairs_of(th, [0.33, 0.44, 0.6, -0.8]) == [(0.33, 0.44)] and M.pairs_of(th[2:], [0.6, -0.8]) == []
    As = [0.77, 0.6, -0.8]
    base = dict(rms=1e-13, err=[1e-9, 1e-9, 1e-9], oob=False)
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8]), As) == 'exact'
    assert M.outcome(dict(base, theta=th + [[0.9, 0.1]], w=[0.33, 0.44, 0.6, -0.8, 1e-6], err=[0.33, 1e-9, 1e-9]), As) == 'coinciding pair'
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8], oob=True), As) == 'stopped by the budget'
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8], rms=0.2), As) == 'not fitted'


def test_the_g_report():
    """G's tables (referee round 57): the stalled block (kernels, max|w|, distance, residual against the dipole infimum, its polishes);
    the flow's largest |w| and distance; the last interval's exponents; K = v²(E − E_inf) at the last two outputs; Δv³/Δs against
    3μ₀K; the moment's σ-share; the loop's snapshots with the closest pair's distance, coefficients and offset."""
    import exp_moving_width as M
    ws = (0.9, 4.6, 20.0, 92.0, 300.0, 425.0)
    seps = (5.9, 0.29, 0.062, 0.013, 0.0042, 0.0029)
    flow = [dict(s=s_, wmax=w, sep=d, rms=0.2 + 1 / w) for s_, w, d in zip(M.S_G, ws, seps)]
    law = [dict(v=v, E=E, E_inf=0.02, rms_inf=0.2, p=[0.1, -0.35], mu0=0.1, dc=dc, ds=ds)
           for v, E, dc, ds in ((300.0, 0.02 + 1e-5, 2e-4, 6e-4), (425.0, 0.02 + 5e-6, 1.4e-4, 4.2e-4))]
    stalled = dict(n=4, wmax=3910.0, sep=1.9e-4, rms=0.2000001, polish=[dict(status=2, nfev=2900, accepted=True, wmax_ratio=4300.0)], nfev=10)
    cs, ss, As = [-0.3, 0.2, 0.5], [0.1, 0.15, 0.08], [0.9, 0.54, -0.6]
    pair = lambda off, d: [[0.2 + off * 0.15 + 3.96 * d * 0.15 / 0.54, 0.15], [0.2 + off * 0.15 + 4.5 * d * 0.15 / 0.54, 0.15]]
    loops = [dict(eps=1e-8, tau=t, Lam=100.0, theta=pair(o, d) + [[-0.3, 0.1], [0.5, 0.08]], w=[4.5, -3.96, 0.9, -0.6], wmax=4.5, rms=r)
             for t, d, o, r in ((400.0, 4.6e-4, 0.008, 1.1e-11), (1600.0, 3e-4, 0.005, 5e-12), (6400.0, 2e-4, 0.002, 1e-12))]
    G = [dict(dmin=0.4, q=9, cs=cs, ss=ss, As=As, entries=2, tau_d=17.0, stalled=stalled, flow=flow, law=law, status=0, loops=loops)]
    md = M.report(G=G)
    ls = np.log(M.S_G[-1] / M.S_G[-2])
    ex = M.dipole_exponents(flow)
    exs = (f"{ex[0]:+.4f}, {ex[1]:+.4f}, {np.log(1.4e-4 / 2e-4) / ls:+.4f}, {np.log(4.2e-4 / 6e-4) / ls:+.4f}, "
           f"{np.log((flow[-1]['rms'] - 0.2) / (flow[-2]['rms'] - 0.2)) / ls:+.4f}")
    K2 = [300.0 ** 2 * 1e-5, 425.0 ** 2 * 5e-6]
    vv = f"{(425.0 ** 3 - 300.0 ** 3) / (M.S_G[-1] - M.S_G[-2]):.5e}, {3 * 0.1 * K2[1]:.5e}"
    assert (f"| 0.4 | 9 | 2 | 4, 3910, 1.90e-04, 2.000001e-01 (2.000000e-01), 2/2900/a/4.3e+03 | 0.9, 4.6, 20, 92, 300, 425 | "
            f"5.90e+00, 2.90e-01, 6.20e-02, 1.30e-02, 4.20e-03, 2.90e-03 | {exs} | {K2[0]:.5e}, {K2[1]:.5e} | {vv} | "
            f"{0.35 / np.hypot(0.1, 0.35):.3f} |") in md
    lp = M.landing(loops[1]['theta'], loops[1]['w'], cs, ss, As)                         # the dipole pair at τ = 1600 (round 61)
    assert abs(lp['off'] - 0.005) < 1e-12 and abs(lp['dW']) < 1e-12 and lp['w'] == (4.5, -3.96)
    assert (f"| 0.4 | 9 | 1e-08 | 1600 | 4 | {lp['d']:.2e}, +4.500/-3.960, 5.0e-03, {lp['dW']:+.1e} L | +4.500/-3.960 | 5.00 | "
            "5.0e-12 |") in md and abs(lp['d'] - 3e-4) < 1e-12
    ex = M.dipole_exponents([dict(s=1e9, wmax=1.0, sep=1.0), dict(s=1e11, wmax=100 ** (1 / 3), sep=100 ** (-1 / 3))])
    assert abs(ex[0] - 1 / 3) < 1e-12 and abs(ex[1] + 1 / 3) < 1e-12 and M.dipole_exponents(flow[:1] + [dict(flow[1], sep=None)]) is None


def test_the_formed_kernels_pairs_and_outcomes():
    """Referee round 57: `n_formed` counts |w| above 1% of the smallest |amplitude|; `pairs_of` at PAIR_TOL in |Δc|/σ + |Δlog σ|;
    `outcome`'s classes at their thresholds — exact (idle entrants allowed), a coinciding pair, a cancellation (a coefficient beyond 1.5
    times the largest amplitude, or a width below the samples' spacing), not fitted (an error 2·10⁻⁶ or a residual 2·10⁻⁸), a dipole
    run-off or another budget stop (the stalled block)."""
    import exp_moving_width as M
    assert M.n_formed(dict(w=[0.009, 0.0095, -0.5]), [1.0, -0.9]) == 2
    T = M.PAIR_TOL
    assert M.pairs_of([[0.0, 0.1], [0.1 * T * 0.99, 0.1]], [1.0, 2.0]) == [(1.0, 2.0)]
    assert M.pairs_of([[0.0, 0.1], [0.1 * T * 1.01, 0.1]], [1.0, 2.0]) == []
    assert M.pairs_of([[0.0, 0.1], [0.0, 0.1 * np.exp(T * 0.99)]], [1.0, 2.0]) and not M.pairs_of([[0.0, 0.1], [0.0, 0.1 * np.exp(T * 1.01)]], [1, 2])
    th = [[-0.5, 0.1], [0.0, 0.1], [0.5, 0.1]]
    m = lambda **kw: dict(dict(w=[1.0, 1.0, 1.0], theta=th, rms=1e-9, err=[1e-7, 1e-7, 1e-7], oob=False), **kw)
    As = [1.0, 1.0, 1.0]
    assert M.outcome(m(), As) == 'exact'
    assert M.outcome(m(w=[1.0, 1.0, 1.0, 1e-3, 1e-3], theta=th + [[0.9, 0.1], [0.9, 0.1]]), As) == 'exact'
    assert M.outcome(m(err=[2e-6, 1e-7, 1e-7]), As) == 'not fitted' and M.outcome(m(rms=2e-8), As) == 'not fitted'
    assert M.outcome(m(w=[1.0, 0.5, 0.5, 1.0], theta=[th[0], th[1], th[1], th[2]], err=[1e-7, 0.5, 1e-7]), As) == 'coinciding pair'
    assert M.outcome(m(w=[1.6, 1.0, 1.0], err=[0.6, 1e-7, 1e-7]), As) == 'cancellation'
    assert M.outcome(m(theta=[[-0.5, 0.005], th[1], th[2]], err=[0.1, 1e-7, 1e-7]), As) == 'cancellation'
    st = dict(theta=[[0.0, 0.1], [0.0002, 0.1], [0.5, 0.1]], w=[500.0, -499.0, 1.0])
    assert M.outcome(m(oob=True, stalled=st), As) == 'dipole run-off'
    assert M.outcome(m(oob=True, stalled=dict(st, w=[5.0, -4.0, 1.0])), As) == 'stopped by the budget'
    assert M.outcome(m(oob=True), As) == 'stopped by the budget'


def test_the_width_split_cancels_the_third_moment():
    """Referee round 57: a pair of unequal shares Δc apart fits one Gaussian to O(Δc⁴) when its widths split by (w_j − w_i)Δc²/(6Wσ);
    with equal widths the best single Gaussian leaves the third moment, O(Δc³)."""
    import exp_moving_width as M
    from scipy.optimize import least_squares
    x = np.linspace(-1, 1, 801)
    g = lambda c, s: np.exp(-(x - c) ** 2 / (2 * s ** 2))
    wi, wj, sb = 0.3, 0.5, 0.1
    W = wi + wj
    res = []
    for dc in (0.02, 0.01):
        ci, cj = -wj * dc / W, wi * dc / W
        _, _, pred = M.width_split([[ci, sb], [cj, sb]], [wi, wj], 0, 1)
        assert abs(pred - (wj - wi) * dc ** 2 / (6 * W * sb)) < 1e-15
        for d in (0.0, pred):
            f = wi * g(ci, sb - wj * d / W) + wj * g(cj, sb + wi * d / W)
            sol = least_squares(lambda z: f - z[0] * g(z[1], z[2]), [W, 0.0, sb], xtol=1e-15, ftol=1e-15, gtol=1e-15)
            res.append(np.sqrt(np.mean(sol.fun ** 2)))
    (e0, e1), (h0, h1) = res[:2], res[2:]
    assert e1 < 0.2 * e0 and 7 < e0 / h0 < 9 and 14 < e1 / h1 < 18               # O(Δc³) and O(Δc⁴): 8 and 16 per halving
    dc_, ds_, pr_ = M.width_split([[0.01, 0.11], [-0.01, 0.1]], [0.5, 0.3], 0, 1)        # sorted by centre: i the left kernel
    assert dc_ == 0.02 and abs(ds_ - 0.01) < 1e-15 and abs(pr_ - (0.5 - 0.3) * 0.02 ** 2 / (6 * 0.8 * 0.105)) < 1e-15


def test_the_dipole_law_finds_the_infimum():
    """Referee round 57: `dipole_law` on a pair ±v coalescing around a dipole p·∇_θφ beside two kernels: the pair (the largest |w| and its
    nearest kernel of the other sign), v, the moment, and E_inf (the dipole fit's) ≈ 0 for a target that is that dipole exactly (the pair's
    own error the third-order term p³/(24v²)·∇³φ)."""
    import exp_moving_width as M
    L = M.law()
    K = M.MW.BumpsCW(M.XS, M.S0)
    tb, p = np.array([0.1, 0.12]), np.array([0.002, -0.004])
    f = K.synth(tb[None, :], np.array([0.5])) + K.dphi(tb[None, :])[:, 0, :] @ p
    f = f + K.synth(np.array([[0.6, 0.1]]), np.array([0.4]))
    v = 10.0
    th = np.array([tb + p / (2 * v), tb - p / (2 * v), tb, [0.6, 0.1]])
    d = M.dipole_law(L, K, f, th, np.array([v, -v, 0.5, 0.4]))
    assert d['v'] == v and np.allclose(d['p'], p, rtol=1e-6) and d['E_inf'] < 1e-25 < d['E'] and d['mu0'] == 4 * L.G * L.c / L.A


def test_the_outcomes_at_their_thresholds():
    """Referee round 61: `outcome` with the features — exact with an idle null cluster (each feature's nearest kernel exact, the others
    summing to nothing); a coinciding pair only where its net coefficient and centre of mass fit a feature (else a null cluster), among
    formed kernels (|w| above 1% of the smallest amplitude) and only below the residual 1e-8; a budget stop a dipole run-off only for an
    opposite-sign pair closer than 0.05 with |w| beyond 10 times the largest amplitude."""
    import exp_moving_width as M
    cs, ss, As = [-0.5, 0.0, 0.5], [0.1, 0.1, 0.1], [1.0, -0.8, 0.6]
    th = [[c, s] for c, s in zip(cs, ss)]
    m = lambda **kw: dict(dict(w=list(As), theta=th, rms=1e-10, err=[1e-8] * 3, oob=False), **kw)
    oc = lambda mm: M.outcome(mm, As, cs, ss)
    assert oc(m()) == 'exact'
    nul = [[0.9, 0.05], [0.9, 0.05], [0.9, 0.05]]
    assert oc(m(w=list(As) + [-0.05, 0.1, -0.05], theta=th + nul)) == 'exact, null cluster'
    assert oc(m(w=list(As) + [-0.005, 0.005], theta=th + nul[:2])) == 'exact'                     # below 1%: idle
    pr = dict(w=[1.0, -0.3, -0.5, 0.6], theta=[th[0], th[1], [0.0, 0.1 * (1 + 1e-3)], th[2]], err=[1e-8, 0.5, 1e-8])
    assert oc(m(**pr)) == 'coinciding pair'
    assert oc(m(**dict(pr, w=[1.0, -0.3, -0.6, 0.6]))) == 'coinciding pair off the features'      # the pair's W −0.9 against −0.8
    assert oc(m(**dict(pr, theta=[th[0], [0.002, 0.1], [0.002, 0.1], th[2]]))) == 'coinciding pair off the features'   # 0.02σ* off
    assert oc(m(**dict(pr, rms=2e-8))) == 'not fitted'
    assert oc(m(**dict(pr, w=[1.0, -0.795, -0.005, 0.6]))) == 'not fitted'                       # the second kernel unformed
    st = dict(theta=[[0.0, 0.1], [0.004, 0.1], [0.5, 0.1]], w=[50.0, -49.0, 0.6])                # 0.04 apart, 50 = 50·max|A|
    assert oc(m(oob=True, stalled=st)) == 'dipole run-off'
    assert oc(m(oob=True, stalled=dict(st, theta=[[0.0, 0.1], [0.006, 0.1], [0.5, 0.1]]))) == 'stopped by the budget'
    assert oc(m(oob=True, stalled=dict(st, w=[9.9, -8.9, 0.6]))) == 'stopped by the budget'
    st2 = dict(theta=[[0.0, 0.1], [0.001, 0.1], [0.3, 0.1]], w=[50.0, 49.0, -98.0])                # the close pair of one sign
    assert oc(m(oob=True, stalled=st2)) == 'stopped by the budget'


def test_the_closest_pair_counts_the_widths():
    """`closest_pair`'s distance |Δc|/σ + |Δlog σ| (σ the mean width): two kernels at one centre are |log(σ_j/σ_i)| apart."""
    import exp_moving_width as M
    th, w = [[0.0, 0.1], [0.0, 0.2], [0.05, 0.1]], [1.0, -1.0, 1.0]
    d, i, j = M.closest_pair(th, w)
    assert (i, j) == (0, 2) and abs(d - 0.5) < 1e-12
    d, i, j = M.closest_pair(th, w, opposite=True)
    assert (i, j) == (0, 1) and abs(d - np.log(2)) < 1e-12


def test_the_dipole_law_picks_the_nearest_opposite_kernel_and_halves_the_energy():
    """Referee round 61: `dipole_law`'s pair is the largest |w| with the nearest kernel of the other sign (not a farther one); its E and
    E_inf are ½⟨r²⟩ — with a target the dipole model cannot fit (a narrow bump far away), E_inf = ½⟨g²⟩ of that bump."""
    import exp_moving_width as M
    L = M.law()
    K = M.MW.BumpsCW(M.XS, M.S0)
    tb, p, v = np.array([0.1, 0.12]), np.array([0.002, -0.004]), 10.0
    dip = K.synth(tb[None, :], np.array([0.5])) + K.dphi(tb[None, :])[:, 0, :] @ p
    th = np.array([tb + p / (2 * v), tb - p / (2 * v), [-0.9, 0.1], tb])
    d = M.dipole_law(L, K, dip + K.synth(np.array([[-0.9, 0.1]]), np.array([-0.3])), th, np.array([v, -v, -0.3, 0.5]))
    assert abs(d['dc'] - abs(p[0]) / v) < 1e-15
    g = K.synth(np.array([[1.3, 0.03]]), np.array([0.2]))
    d = M.dipole_law(L, K, dip + g, th[[0, 1, 3]], np.array([v, -v, 0.5]))
    assert abs(d['E_inf'] / (0.5 * np.mean(g ** 2)) - 1) < 1e-6 and abs(d['rms_inf'] ** 2 / 2 - d['E_inf']) < 1e-18
    assert abs(d['E'] - 0.5 * np.mean((dip + g - K.synth(th[[0, 1, 3]], np.array([v, -v, 0.5]))) ** 2)) < 1e-18


def test_h_statistics():
    """Referee round 61: H's statistics — the most kernels below −1% of the smallest amplitude in a block (a coefficient −0.005 is not
    one), the most negative coefficient, the end's formed kernels of the other sign, the smallest formed width (not an idle kernel's),
    and over all the end's kernels those narrower than the samples' spacing and those of negative width."""
    import exp_moving_width as M
    As = (1.0, 0.9)
    vis = lambda n: np.c_[np.linspace(-0.2, 0.2, n), np.full(n, 0.06)]
    blocks = [(vis(2), np.array([0.5, -0.005])), (vis(4), np.array([1.0, -0.02, -0.3, 0.9])), (vis(2), np.array([1.0, -0.4]))]
    m = dict(w=[1.0, -0.2, 0.9, 1e-4, -2e-4], theta=[[-0.08, 0.06], [0.0, 0.05], [0.08, 0.07], [0.3, 0.004], [0.4, -0.003]])
    h = M.h_stats(blocks, m, As)
    assert h == dict(neg_path=2, wmin_path=-0.4, neg_end=1, neg_end_all=1, smin=0.05, sub=2, negw=1, sub_formed=0, negw_formed=0,
                     inv_formed=0, wmax=1.0)
    cl = M.clusters(m, [-0.08, 0.08], [0.06, 0.06], As)
    assert [n for n, x in cl] == [2, 1] and abs(cl[0][1] - (0.06 - 0.2 * 0.05) / 0.06) < 1e-12 and abs(cl[1][1] - 0.9 * 0.07 / (0.9 * 0.06)) < 1e-12


def test_the_stalled_block_summary():
    """Referee round 61: G's stalled block — its closest pair of opposite signs (not a closer one of one sign), its largest |w|."""
    import exp_moving_width as M
    K = M.MW.BumpsCW(M.XS, M.S0)
    st = dict(c=[[0.0, 0.1], [0.001, 0.1], [0.3, 0.1]], w=np.array([50.0, 49.0, -98.0]))
    Rc = dict(stalled=st, polish_log=[[dict(status=2, nfev=10, accepted=True)]], nfev=7)
    f = np.zeros(M.XS.size)
    s = M.stalled_summary(Rc, f, K)
    assert s['wmax'] == 98.0 and abs(s['sep'] - 2.99) < 1e-12 and s['n'] == 3 and s['nfev'] == 7 and s['polish'][0]['nfev'] == 10
    assert M.stalled_summary(dict(Rc, stalled=None), f, K) is None


def test_the_loop_phases():
    """Referee round 61: `phases` on a synthetic loop whose dipole pair grows exactly as v³ = 3μ₀KΛ(τ − τ_d) until a third kernel forms
    at τ₃, then lands on the feature with v frozen: the slopes are 1, the landed v the law's (3μ₀KΛ(τ₃ − τ_d))^{1/3}; the collapse's
    exponents from the pair's distance and the residual over 1600 → 6400."""
    import exp_moving_width as M
    assert all(t in M.TAUS_G for t in M.TAUS_G_MD) and 30.0 in M.TAUS_G
    cs, ss, As = [0.0, 0.5], [0.1, 0.1], [0.6, -0.5]
    law = [None, dict(v=100.0, E=0.01 + 1e-6, E_inf=0.01, mu0=0.1)]
    K_, Lam, td, t3 = 1e4 * 1e-6, 50.0, 17.0, 30.0
    rate = 3 * 0.1 * K_ * Lam
    loops = []
    for tau in M.TAUS_G:
        v = (rate * (min(tau, t3) - td)) ** (1 / 3) if tau > td else 0.5
        d = 1e-3 * (tau / 400.0) ** -0.2 if tau >= t3 else 0.01
        W = 0.6 if tau >= t3 else 0.3
        wa, wb = W / 2 + v, W / 2 - v
        th = [[-wb * d * 0.1 / W, 0.1], [wa * d * 0.1 / W, 0.1]]
        ex = [[0.5, 0.1]] if tau >= t3 else []
        loops.append(dict(eps=1e-8, tau=tau, Lam=Lam, theta=th + ex, w=[wa, wb] + ([-0.5] if tau >= t3 else []), wmax=abs(wa),
                          rms=1e-6 * (tau / 400.0) ** -0.8))
    r = dict(cs=cs, ss=ss, As=As, entries=2, tau_d=td, law=law, loops=loops)
    ph = M.phases(r, 1e-8)
    assert ph['t3'] == 30.0 and len(ph['slopes']) >= 4 and max(abs(s - 1) for s in ph['slopes']) < 1e-9
    assert abs(ph['v400'] / ph['vstar'] - 1) < 1e-12 and abs(ph['rate'] - rate) < 1e-12
    assert all(ph['land'][t]['off'] < 1e-9 and abs(ph['land'][t]['dW']) < 1e-12 for t in M.TAUS_G_MD)
    assert abs(ph['collapse'][0] + 0.2) < 1e-9 and abs(ph['collapse'][1] + 0.8) < 1e-9
    md = '\n'.join(M.g_loop_lines([dict(r, dmin=0.4, q=9, loops=loops)]))
    assert (f"| 0.4 | 9 | 1e-08 | 17 | 30 | 1.000 [1.000, 1.000] | {ph['v400']:.4g} | {ph['vstar']:.4g} | yes, yes, yes | -0.200, -0.800 |"
            in md)


def test_a_pair_fits_any_feature_sharing_its_centre():
    """Referee round 61: a pair fits a feature if its net coefficient and centre of mass match any feature's — two features may share a
    centre (D at d = 0: the narrow feature split in two halves on the pedestal), and the nearest by centre alone is then ambiguous."""
    import exp_moving_width as M
    th = [[0.0, 0.06], [2e-4, 0.06]]
    assert M.pair_fits(th, [0.25, 0.25], 0, 1, [0.0, 0.0], [0.3, 0.06], [0.8, 0.5])
    assert not M.pair_fits(th, [0.25, 0.20], 0, 1, [0.0, 0.0], [0.3, 0.06], [0.8, 0.5])
    assert not M.pair_fits(th, [0.25, 0.25], 0, 1, [0.1, 0.1], [0.3, 0.06], [0.8, 0.5])


def test_the_oracle_fit_starts_at_the_truth():
    """Referee round 61 (E4): the oracle fit is the K-kernel least-squares fit of the data started at the truth — without noise it is
    the truth itself (error 0), with noise its error is below the noise's rms."""
    import exp_moving_width as M
    cs, ss, As = [-0.4, 0.3], [0.08, 0.15], [0.7, -0.5]
    K = MW.BumpsCW(M.XS, M.S0)
    f = K.synth(np.c_[cs, ss], np.array(As))
    assert M.oracle_error(f, np.zeros_like(f), cs, ss, As) < 1e-12
    xi = 1e-3 * np.random.default_rng(1).standard_normal(f.size)
    e = M.oracle_error(f, xi, cs, ss, As)
    assert 1e-5 < e < 1e-3


def test_pair_fits_and_pair_index_at_their_boundaries():
    """Referee round 68 (F2, F4, F5, I1): a pair fits a feature only with W of the feature's sign and size (not −A, not 0) and its centre
    of mass — not its midpoint — at the feature; `pair_index` pairs kernels only when both centre and width agree."""
    th = np.array([[0.0, 0.06], [0.02, 0.06]])
    w = np.array([1.5, -0.5])                                          # W = 1, centre of mass −0.01, midpoint +0.01
    assert MWE.pair_fits(th, w, 0, 1, [-0.01], [0.06], [1.0])
    assert not MWE.pair_fits(th, w, 0, 1, [0.01], [0.06], [1.0])      # the midpoint is not the centre of mass
    assert not MWE.pair_fits(th, -w, 0, 1, [-0.01], [0.06], [1.0])    # W = −1: the wrong sign
    assert not MWE.pair_fits(th, np.array([0.5, -0.5]), 0, 1, [0.0], [0.06], [0.0])   # W = 0 fits nothing
    assert MWE.pair_index(np.array([[0.1, 0.06], [0.1, 0.0601]]), [1.0, 1.0]) == [(0, 1)]
    assert MWE.pair_index(np.array([[0.1, 0.3], [0.1, 0.06]]), [1.0, 1.0]) == []      # one centre, two widths


def test_h_stats_counts_negatives_below_one_percent_of_the_smallest_amplitude():
    """Referee round 68 (H1, H7): a block's only negative at −0.5% of the smallest amplitude is not counted, one at −2% is; the
    threshold is the smallest amplitude's, not the largest's."""
    m = dict(w=[1.0, -0.005], theta=[[0.0, 0.06], [0.1, 0.06]])
    th = np.array(m['theta'])
    r = MWE.h_stats([(th, np.array([1.0, -0.005]))], m, (1.0, 0.6))
    assert r['neg_path'] == 0 and r['neg_end'] == 0 and r['wmin_path'] == -0.005
    r = MWE.h_stats([(th, np.array([1.0, -0.008]))], dict(m, w=[1.0, -0.008]), (1.0, 0.6))
    assert r['neg_path'] == 1 and r['neg_end'] == 1                    # −0.008 is beyond 1% of 0.6, not of 1.0


def test_stop_anatomy():
    """Referee round 68: `stop_anatomy` counts the active kernels (all: the dictionary exhausted), the accepted polishes and those stopped
    by their budget, the kernels of negative width and those invisible on the samples."""
    d = MWE.C0.size
    act = np.ones(d, bool)
    c = np.c_[np.linspace(-1, 1, d), np.full(d, 0.06)]
    c[3] = [0.775, -3e-4]                                               # negative width, invisible on the samples
    log = [[dict(status=0, accepted=True), dict(status=1, accepted=False)], [dict(status=2, accepted=True)]]
    r = MWE.stop_anatomy(dict(act=act, c=c, polish_log=log))
    assert r == dict(n_act=d, exhausted=True, pol_acc=2, pol_capped=1, neg_sigma=1, invisible=1)
    act[0] = False
    assert MWE.stop_anatomy(dict(act=act, c=c, polish_log=[]))['exhausted'] is False


def test_fine_phases_with_an_intermediate_phase():
    """Referee round 68: `fine_phases` finds τ₃ and τ₄ as formations, v before τ₃, the intermediate phase between them, the first
    landing after τ₃ with its p, and counts where p ≤ P_LAND agrees with the landing."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    dip = lambda v, c: (np.array([[c - 0.01, 0.06], [c + 0.01, 0.06]]), np.array([1.0 + v, -v]))
    sts = [(10.0, *dip(0.5, 0.2)), (12.0, *dip(0.6, 0.2))]               # the dipole phase (n0 = 2), off the features
    th, w = dip(0.7, 0.2)
    sts.append((14.0, np.r_[th, [[0.5, 0.06]]], np.r_[w, 0.8]))          # τ₃: a third kernel forms
    sts.append((20.0, np.r_[th, [[0.5, 0.06]]], np.r_[w, 0.8]))          # the intermediate phase
    th2, w2 = np.array([[-0.0001, 0.06], [0.0003, 0.06]]), np.array([1.05, -0.05])
    sts.append((26.0, np.r_[th2, [[0.5, 0.06], [0.9, 0.06]]], np.r_[w2, 0.8, 0.01]))   # τ₄, landed
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['v3']) == (14.0, 26.0, 1.1)                # v = |w_a − w_b|/2 at τ = 12
    assert r['inter']['n'] == 3 and r['inter']['off'] > 1 and r['t_land'] == 26.0 and abs(r['v_land'] - 0.55) < 1e-12
    p = 1.1 * 0.0004 / (1.0 * 0.06)
    assert abs(r['p_land'] - p) < 1e-12 and r['agree'][1] >= 1


def test_fine_phases_ignores_an_early_landing_and_an_exit():
    """Referee round 80 (FPf, FPg, FPh): a landing before τ₃ is not the landing; a state with fewer formed kernels after τ₃ (an exit) is
    not τ₄; the intermediate phase is read at the first state from (τ₃ + τ₄)/2, not at τ₃; the intermediate phase's v only at n₃
    formed."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    off = lambda c: (np.array([[c - 0.01, 0.06], [c + 0.01, 0.06]]), np.array([1.5, -0.5]))
    landed = (np.array([[-0.0001, 0.06], [0.0003, 0.06]]), np.array([1.05, -0.05]))
    third = (np.array([[0.5, 0.06]]), np.array([0.8]))
    cat = lambda a, *b: (np.vstack([a[0]] + [x[0] for x in b]), np.concatenate([a[1]] + [x[1] for x in b]))
    sts = [(10.0, *off(0.2)), (12.0, *landed),                            # the dipole phase, landed early at τ = 12
           (14.0, *cat(off(0.2), third)),                                 # τ₃ (n₃ = 3), the pair at 0.2: 3σ* from feature 1
           (16.0, *off(0.2)),                                             # an exit: two formed
           (18.0, *cat(off(0.2), third)), (20.0, *cat(off(0.3), third)),   # the midpoint (14 + 26)/2: the pair at 0.3
           (26.0, *cat(landed, third, (np.array([[0.9, 0.06]]), np.array([0.01]))))]   # τ₄, landed
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['t_land']) == (14.0, 26.0, 26.0)
    assert abs(r['inter']['off'] - abs(0.28 - 0.5) / 0.06) < 1e-9 and r['inter']['n'] == 3
    assert r['dip2'] == [(14.0, 1.0), (20.0, 1.0)]                    # every tenth of 14, 18, 20 (16 has two formed) and the last


RAW_MW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'moving_width_raw.npz')


def test_t3_from_duals_reproduces_gt_and_the_loop():
    """Referee round 68: τ₃ from the sliding path alone — the duals advanced to the dipole block's entry τ_d, then moved by the block's
    residual at the end of its flow — reproduces GT's raw for target 9 and lies within the fine grid's step (0.1) of the loop's τ₃ at the
    smallest mass, for every G target (raw)."""
    raw = MWE.E.load(RAW_MW)
    gt = {x['q']: x['t3_pred'] for x in raw['GT']}
    for q in (9, 20):                                                   # referee round 80 (GT1): gt_job's every field against GT's raw
        g = next(r for r in raw['G'] if r['q'] == q)
        o = MWE.gt_job((q, np.array(g['cs']), np.array(g['ss']), np.array(g['As']), int(g['entries']), g['flow'][-1]['theta'],
                        g['flow'][-1]['w']))
        x = next(r for r in raw['GT'] if r['q'] == q)
        assert set(o) == set(x)
        for k_, v in o.items():
            assert ((v is None and x[k_] is None) or (isinstance(v, bool) and v == x[k_]) or (isinstance(v, list) and v == list(x[k_]))
                    or np.isclose(v, x[k_], rtol=1e-9, atol=0)), k_
    assert next(r for r in raw['GT'] if r['q'] == 20)['exact'] is True
    n = 0
    for r in raw['G']:
        if r['q'] in gt and r.get('fine'):
            ph = r['fine'][min(r['fine'])]
            assert abs(gt[r['q']] - ph['t3']) < 0.1
            n += 1
    assert n == 6


def test_hl_job_on_a_cheap_case(monkeypatch):
    """Referee round 68 (R8): `hl_job` on H's separated target at one mass and time — the formed set is the loop's kernels above 1% of the
    smallest amplitude, with its count of negatives, smallest width, largest |w| over all kernels and the residual."""
    monkeypatch.setattr(MWE, 'EPS_HL', (1e-6,))
    monkeypatch.setattr(MWE, 'TAUS_HL', (60.0,))
    name, cs, ss, As = next(t for t in MWE.H_TARGETS if t[0] == 'σ* 0.06, A₂ 0.6')
    r = MWE.hl_job((name, cs, ss, As))
    L, K = MWE.law(), MWE.MW.BumpsCW(MWE.XS, MWE.S0)
    f = MWE.feature_sum(cs, ss, As)
    _, W, C, Lam, _ = MWE.RN.loop_nd(L, K, f, MWE.theta0(), 1e-6, np.r_[0.0, 60.0], rtol=1e-9, atol=1e-12)
    big = np.flatnonzero(np.abs(W[1]) > 1e-2 * 0.6)
    (q,) = r['loops']
    assert np.allclose(q['w'], W[1][big], rtol=1e-12, atol=0) and np.allclose(q['theta'], C[1][big], rtol=1e-12, atol=0)
    assert q['neg'] == int(np.sum(W[1][big] < 0)) and q['wmax'] == float(np.abs(W[1]).max())
    assert q['smin'] == float(np.abs(C[1][big][:, 1]).min())
    assert np.isclose(q['rms'], np.sqrt(np.mean((f - K.synth(C[1], W[1])) ** 2)), rtol=1e-12, atol=0)


def test_h_job_with_a_tolerance_and_h3_targets():
    """Referee round 68 (R7, R9, R10): `h_job` with a tolerance records the last entry's τ (the path's last event, not its first), the last
    blocks' residuals (the history's last six) and the clusters; H3's targets are equal amplitudes 1 at ±c for every separation."""
    name, cs, ss, As = next(t for t in MWE.H_TARGETS if t[0] == 'σ* 0.06, A₂ 0.6')
    r = MWE.h_job((name, cs, ss, As, 1e-13))
    Rm = MWE.RN.sliding_iss_nd(MWE.law(), MWE.MW.BumpsCW(MWE.XS, MWE.S0), MWE.feature_sum(cs, ss, As), MWE.theta0(),
                               np.ones(MWE.C0.size), max_events=40, exit_kkt=True, max_rhs=MWE.MAX_RHS_MOVING, tol=1e-13)
    ent = [e for e in Rm['events'] if e[1] is not None]
    assert len(ent) >= 2 and r['last_tau'] == float(Rm['events'][-1][0]) and r['last_tau'] > float(Rm['events'][0][0])
    assert r['rms_blocks'] == [float(h[3]) for h in Rm['hist'][-6:]]
    assert r['outcome'] == 'exact'
    for nm, c, s_, A in MWE.h3_targets():
        assert tuple(A) == (1.0, 1.0) and tuple(s_) == (0.06, 0.06) and c[0] == -c[1] and nm.endswith(f"±{c[1]:g}")


def test_oracle_error_against_a_direct_fit():
    """Referee round 68 (E2, E3): the oracle's error is the rms against f of the K-kernel least-squares fit of the data f + ξ (not of
    f + ½ξ, and not measured against the data), started at the truth — checked against a direct fit."""
    from scipy.optimize import least_squares
    cs, ss, As = (-0.3, 0.25), (0.05, 0.08), (1.0, -0.7)
    f = MWE.feature_sum(cs, ss, As)
    xi = 0.01 * np.random.default_rng(5).standard_normal(MWE.XS.size)
    K = MWE.MW.BumpsCW(MWE.XS, MWE.S0)
    sol = least_squares(lambda z: f + xi - K.synth(z[2:].reshape(2, 2), z[:2]), np.r_[As, np.c_[cs, ss].ravel()], method='lm',
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    e = np.sqrt(np.mean((K.synth(sol.x[2:].reshape(2, 2), sol.x[:2]) - f) ** 2))
    assert np.isclose(MWE.oracle_error(f, xi, cs, ss, As), e, rtol=1e-10, atol=0)
    assert abs(e - np.sqrt(np.mean(xi ** 2))) > 1e-3


def test_c_table_counts_a_budget_stop_at_its_stalled_block():
    """Referee round 68 (D3bo-8): C's table counts a budget-stopped path at its stalled block — the residual G's `stalled_summary` records
    for the same target, the formed kernels the stalled block's (raw)."""
    raw = MWE.E.load(RAW_MW)
    g = {(r['dmin'], r['q']): r for r in raw['G']}
    n = 0
    for r in raw['C']:
        if r['moving'].get('stalled') is not None:
            k, rms = MWE.moving_state(r)
            st = g[(r['dmin'], r['q'])]['stalled']
            assert np.isclose(rms, st['rms'], rtol=1e-9, atol=0) and k == int(np.sum(np.abs(r['moving']['stalled']['w']) > 1e-2 * min(abs(a) for a in r['As'])))
            n += 1
    assert n == 6


def test_visibility_and_h_stats_count_only_visible_kernels():
    """Referee round 74 (D3bu-2): a kernel is visible on the samples if its largest φ there is at least VIS — a spike between two
    samples (σ = 3·10⁻⁴, 0.005 from the nearest) is not, one of largest φ ≈ 3.7·10⁻⁶ (σ = 10⁻³) is —; H's negative kernels along the
    path and at the end count the visible only (the end's formed ones of the other sign including the invisible in neg_end_all), the
    formed invisible ones are counted, and the end's sub-sample and negative-width kernels with how many are formed."""
    v = MWE.visibility([[0.005, 3e-4], [0.005, 1e-3], [0.0, 0.06], [0.2, 0.0]])
    assert v[0] < MWE.VIS < v[1] and abs(v[1] - np.exp(-0.005 ** 2 / (2e-6))) < 1e-18 and v[2] == 1.0 and v[3] == 0.0
    th = np.array([[-0.08, 0.06], [0.005, 3e-4], [0.08, 0.06], [0.3, -0.06]])
    w = np.array([1.2, -0.7, 1.1, -0.3])
    r = MWE.h_stats([(th, w)], dict(theta=th.tolist(), w=w.tolist()), (1.0, 1.0))
    assert (r['neg_path'], r['wmin_path'], r['neg_end'], r['neg_end_all'], r['inv_formed']) == (1, -0.3, 1, 2, 1)
    assert (r['sub'], r['sub_formed'], r['negw'], r['negw_formed']) == (1, 1, 1, 1) and r['smin'] == 0.06 and r['wmax'] == 1.2


def test_h_stats_does_not_count_an_unformed_invisible_kernel():
    """Referee round 80 (HS3): `h_stats`' formed invisible kernels are formed (|w| above 1% of the smallest amplitude) and invisible on
    the samples; an unformed invisible one (|w| 10⁻³, a spike between samples) is not counted."""
    m = dict(w=[1.0, 1e-3, 0.5], theta=[[0.0, 0.06], [0.0049, 1e-5], [0.2049, 1e-5]])
    st = MWE.h_stats([], m, [1.0])
    assert st['inv_formed'] == 1 and st['sub'] == 2 and st['sub_formed'] == 1


def test_path_stop_reasons():
    """Referee round 74: why a path ended — its tolerance, its budget, the dictionary exhausted, the events' cap, no entrant."""
    h = lambda rms: [(0.0, None, None, 1.0), (1.0, None, None, rms)]
    act = np.ones(21, bool)
    assert MWE.path_stop(dict(hist=h(1e-10), out_of_budget=True, act=act), 1e-9) == 'tolerance'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=True, act=act), 1e-9) == 'budget'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=False, act=act), 1e-9) == 'dictionary'
    act[3] = False
    assert MWE.path_stop(dict(hist=h(1e-6) * 20, out_of_budget=False, act=act), 1e-9) == 'events'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=False, act=act), 1e-9) == 'no entrant'


def test_h_job_passes_its_polish_and_tolerance_and_clusters_unequal_amplitudes(monkeypatch):
    """Referee round 74 (HJ1, HJ2, HJ5, HJ6): `h_job` with (…, None, False) runs the path without the polish at the tolerance 1e-9, with
    a tolerance at that one with the polish; its stop and its clusters' masses Σw_i|σ_i|/(A_kσ_k) with unequal amplitudes."""
    seen = []
    th = np.array([[-0.08, 0.06], [0.081, 0.061], [0.5, 0.06]])
    w = np.array([1.0, 0.6, 0.0])

    def fake(L, K, f, c0, rho, **kw):
        seen.append(kw)
        act = np.array([True, True] + [False] * (c0.shape[0] - 2))
        c = np.array(c0, float); c[:2] = th[:2]
        ww = np.zeros(c0.shape[0]); ww[:2] = w[:2]
        return dict(act=act, c=c, w=ww, hist=[(0.0, c, np.zeros_like(ww), 1.0), (1.0, c, ww, 2e-10)], events=[(1.0, 0, 'enter', c[0])],
                    out_of_budget=False, nfev=10, polish_log=[], stalled=None)
    monkeypatch.setattr(MWE.RN, 'sliding_iss_nd', fake)
    name, cs, ss, As = 'x', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)
    r = MWE.h_job((name, cs, ss, As, None, False))
    assert seen[-1]['polish'] is False and seen[-1]['tol'] == 1e-9 and r['polish'] is False and r['stop'] == 'tolerance'
    r = MWE.h_job((name, cs, ss, As, 1e-13))
    assert seen[-1]['polish'] is True and seen[-1]['tol'] == 1e-13 and r['stop'] == 'no entrant'
    assert r['clusters'] == [(1, 1.0 * 0.06 / (1.0 * 0.06)), (1, 0.6 * 0.061 / (0.6 * 0.06))]


def test_stop_anatomy_counts_a_faint_kernel_as_visible():
    """Referee round 74 (SA4): a kernel of largest φ ≈ 10⁻⁶ on the samples is visible (VIS = 10⁻¹⁰), one between samples with
    σ = 3·10⁻⁴ is not."""
    c = np.array([[0.005, 1.2e-3], [0.005, 3e-4], [0.0, 0.06]])
    r = MWE.stop_anatomy(dict(act=np.array([True, True, True]), c=c, polish_log=[]))
    assert r['invisible'] == 1 and MWE.visibility(c[:1])[0] > 1e-7


def test_d3phi_and_the_third_order_k():
    """Referee rounds 74, 80: `d3phi` is the third derivative of φ(θ̄ + tp) in t (against differences), and `k_third` is the dipole's K: a
    pair at v, w = W/2 ± v, Δθ = p/v about the optimal dipole fit's θ̄ (`closure_fit`) has v²(E − E_inf) → K, the relative error
    falling like 1/v² — E − E_inf = K/v² + O(v⁻⁴), the pair's sum even in 1/v about the midpoint (D3ca-7)."""
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    tb, p = np.array([0.05, 0.11]), np.array([0.013, -0.007])
    ph = lambda t: K.phi((tb + t * p)[None, :])[:, 0]
    h = 1e-3
    fd = (ph(2 * h) - 2 * ph(h) + 2 * ph(-h) - ph(-2 * h)) / (2 * h ** 3)
    assert np.abs(MWE.d3phi(tb, p) - fd).max() < 1e-4 * np.abs(fd).max()
    f = MWE.feature_sum([0.0, 0.03], [0.06, 0.05], [1.0, -0.6])
    cf = MWE.closure_fit(K, f, np.zeros((0, 2)), np.zeros(0), 0.4, np.array([0.015, 0.06]), np.array([0.015, 0.0]))
    k = MWE.k_third(cf['r'], cf['tb'], cf['p'])
    e_inf = np.mean(cf['r'] ** 2) / 2

    def est(v):
        d = cf['p'] / v
        g = K.synth(np.array([cf['tb'] + d / 2, cf['tb'] - d / 2]), np.array([cf['W'] / 2 + v, cf['W'] / 2 - v]))
        return v ** 2 * (np.mean((f - g) ** 2) / 2 - e_inf)
    e1, e2 = abs(est(30.0) / k - 1), abs(est(300.0) / k - 1)       # beyond v ~ 10³ the pair's cancellation rounds E − E_inf
    assert e2 < 0.05 * e1 and e2 < 1e-3 and k != 0 and np.sqrt(2 * e_inf) < 1e-3
    c30, c100 = (est(30.0) / k - 1) * 30.0 ** 2, (est(100.0) / k - 1) * 100.0 ** 2
    assert abs(c30 / c100 - 1) < 0.01 and c30 < 0                     # O(v⁻²) relative, not O(v⁻¹)


def test_closure_fit_returns_positive_widths():
    """Referee round 80 (D3ca-10): φ is even in σ, so a fit may end at a negative width; `closure_fit` returns |σ| — the other kernel
    started at σ = −0.06 is returned as the fit started at +0.06 is, at a width within 10⁻⁶ of 0.06, with the same residual."""
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    f = MWE.feature_sum([0.0, 0.03, 0.5], [0.06, 0.05, 0.06], [1.0, -0.6, 0.7])
    fit = lambda s_: MWE.closure_fit(K, f, np.array([[0.5, s_]]), np.array([0.7]), 0.4, np.array([0.015, 0.06]), np.array([0.015, 0.0]))
    cn, cp = fit(-0.06), fit(0.06)
    assert np.array_equal(cn['th'], cp['th']) and np.array_equal(cn['r'], cp['r']) and abs(cn['th'][0, 1] - 0.06) < 1e-6


def test_first_wall_and_kappa_from_the_path():
    """Referee rounds 74, 80: `first_wall` — the first inactive dual to reach ±1 at the correlations v —; `kappa_from_path` on target 9
    (G's raw block): K from the third order is G's v²(E − E_inf), τ₃ and τ₄ within the fine grid's step of the loop's at the smallest
    mass, κ within 0.2% of the loop's κ′ (its v carried by the law to the path's τ₃: the grid's offset removed); K′ > 0, so no
    dissolution time; the three-kernel closure not exact."""
    j, t = MWE.first_wall(np.array([0.5, -0.2, 0.9]), np.array([False, True, False]), np.array([0.1, 1.0, -0.5]))
    assert j == 2 and abs(t - 3.8) < 1e-12
    assert MWE.first_wall(np.zeros(2), np.array([True, False]), np.array([1.0, 0.0])) == (None, float('inf'))
    raw = MWE.E.load(RAW_MW)
    g = next(r for r in raw['G'] if r['q'] == 9)
    K, L = MW.BumpsCW(MWE.XS, MWE.S0), MWE.law()
    f = MWE.feature_sum(g['cs'], g['ss'], g['As'])
    R = MWE.RN.sliding_iss_nd(L, K, f, MWE.theta0(), np.ones(MWE.C0.size), max_entries=int(g['entries']), exit_kkt=True)
    o = MWE.kappa_from_path(K, f, R, g['flow'][-1]['theta'], g['flow'][-1]['w'])
    lw = g['law'][-1]
    ph = g['fine'][min(g['fine'])]
    assert abs(o['K_third'] / (lw['v'] ** 2 * (lw['E'] - lw['E_inf'])) - 1) < 1e-4
    assert abs(o['t3_path'] - ph['t3']) < 0.1 and abs(o['t4_path'] - ph['t4']) < 0.2
    assert abs(o['kappa_path'] / MWE.kappa_loop(ph, g['tau_d'], o['t3_path'])[0] - 1) < 2e-3
    assert o['exact'] is False and o['t_diss'] is None and o['Kp'] > 0 and o['rms2'] > MWE.CLOSURE_EXACT
    assert abs(o['kappa_path'] ** 3 - 1 - o['Kp'] * (o['t4_path'] - o['t3_path']) / (o['K_third'] * (o['t3_path'] - g['tau_d']))) < 1e-9
    assert o['t3_path'] == MWE.t3_from_duals(K, f, R, g['flow'][-1]['theta'], g['flow'][-1]['w'])
    assert (o['block'], o['j3'], o['j4']) == ([4, 12], 3, 11)            # referee round 83: the path's block and entrants by index


def test_kappa_from_path_dissolution_and_an_exact_closure(monkeypatch):
    """Referee round 80 (D3ca-1, D3ca-6): `kappa_from_path` with its fits and walls stubbed — for K′ < 0 the dissolution time
    τ₃ + K(τ₃ − τ_d)/|K′| and κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3} while τ₄ comes first; no κ when τ₄ comes after it (the bracket
    negative); for K′ > 0 no dissolution; a closure fitting to below CLOSURE_EXACT is exact: K′ = 0, no τ₄, κ = 1."""
    class K_:
        def synth(self, th, w):
            return np.zeros(3)

        def corr(self, cen, r):
            return (np.zeros(4),)
    th, w = np.array([[0.0, 0.06], [0.01, 0.06], [0.5, 0.06]]), np.array([2.0, -1.5, 0.8])

    def run(k2, d4, r2=np.full(3, 1e-3)):
        ks = iter((2e-3, k2))
        walls = iter(((1, 10.0), (2, d4)))
        fits = iter((dict(r=np.full(3, 1e-2), th=th[2:], w=w[2:], W=0.5, tb=np.zeros(2), p=np.ones(2)),
                     dict(r=r2, th=th[1:], w=w[1:], W=0.5, tb=np.zeros(2), p=np.ones(2))))
        monkeypatch.setattr(MWE, 'k_third', lambda *a: next(ks))
        monkeypatch.setattr(MWE, 'first_wall', lambda *a: next(walls))
        monkeypatch.setattr(MWE, 'closure_fit', lambda *a, **k: next(fits))
        monkeypatch.setattr(MWE, 'path_duals_at', lambda *a: (np.zeros(4), np.array([True, False, False, False]), np.zeros((4, 2)), 10.0))
        return MWE.kappa_from_path(K_(), np.zeros(3), None, th, w)
    o = run(-4e-4, 30.0)                                                # τ₃ = 20, τ₄ = 50, dissolution at 20 + 2e-3·10/4e-4 = 70
    assert (o['t3_path'], o['t4_path']) == (20.0, 50.0) and abs(o['t_diss'] - 70.0) < 1e-12 and o['exact'] is False
    assert (o['block'], o['j3'], o['j4']) == ([0], 1, 2)                # referee round 83: the indices
    assert abs(o['kappa_path'] - 0.4 ** (1 / 3)) < 1e-12 and abs(o['rms2'] - 1e-3) < 1e-15
    o = run(-4e-4, 60.0)                                                # τ₄ = 80 after the dissolution: no κ
    assert o['t4_path'] == 80.0 and abs(o['t_diss'] - 70.0) < 1e-12 and o['kappa_path'] is None
    o = run(1e-4, 30.0)
    assert o['t_diss'] is None and abs(o['kappa_path'] - 1.15 ** (1 / 3)) < 1e-12
    o = run(1e-30, 30.0, r2=np.full(3, 1e-17))                          # an exact closure: its K′ and τ₄ are roundoff
    assert o['exact'] is True and o['Kp'] == 0.0 and o['t4_path'] is None and o['kappa_path'] == 1.0 and o['t_diss'] is None
    assert (o['j3'], o['j4']) == (1, None)


def test_part_gt_takes_the_status_zero_entries(monkeypatch):
    """Referee round 74 (PG1): GT runs on G's entries whose block flow ended (status 0), not the others."""
    got = []
    monkeypatch.setattr(MWE, 'gt_job', lambda a: got.append(a[0]) or dict(q=a[0]))
    G = [dict(q=q, status=st, cs=[0.0], ss=[0.1], As=[1.0], entries=2, flow=[dict(theta=[[0.0, 0.1]], w=[1.0])]) for q, st in ((1, 0), (2, -1), (3, 0))]
    assert [x['q'] for x in MWE.part_GT(G)] == [1, 3] and got == [1, 3]


def test_fine_phases_counts_a_sigma_dipole_with_its_width_difference():
    """Referee round 74 (FP6, FP8, FP11, FP13): p takes |Δθ| — a landed σ-dipole (Δc = 0, Δσ ≠ 0) with 0.28 < p < 1 disagrees with the
    landing —, a state with 0.28 < p < 1 off the features agrees, states with p ≥ 1 are not counted; the dipole phase's v is kept at
    every tenth time with the last one's time; the intermediate phase's v and the sign of its (W − A)/|A|."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    sts = []
    for k in range(25):                                                 # the dipole phase: p ≥ 1, off the features
        v = (0.5 * (k + 1)) ** (1 / 3)
        sts.append((6.0 + k, np.array([[0.2 - 0.05, 0.06], [0.2 + 0.05, 0.06]]), np.array([1.0 + v, -v])))
    third = [[0.5, 0.06]]
    for k in range(30):                                                 # τ₃ = 31 … 60: three formed, the pair at W = 0.5 off the features
        v = (sts[-1][2][0] - sts[-1][2][1]) / 2 + 0.01 * k
        sts.append((31.0 + k, np.r_[np.array([[0.3, 0.06], [0.31, 0.06]]), third], np.array([0.25 + v, 0.25 - v, 0.8])))
    land = np.array([[0.0, 0.06], [0.0001, 0.06], [0.5, 0.06], [0.9, 0.06]])
    sts.append((61.0, land, np.array([1.05, -0.05, 0.8, 0.02])))          # τ₄, landed (p small)
    sts.append((62.0, np.array([[0.0, 0.065], [0.0, 0.055], [0.5, 0.06], [0.9, 0.06]]), np.array([2.0, -1.0, 0.8, 0.02])))  # σ-dipole
    sts.append((63.0, np.array([[0.1, 0.06], [0.106, 0.06], [0.5, 0.06], [0.9, 0.06]]), np.array([2.0, -1.0, 0.8, 0.02])))  # off, p 0.3
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['t_v3']) == (31.0, 61.0, 30.0) and [t for t, v in r['dip']] == [6.0, 16.0, 26.0]
    assert [t for t, v in r['dip2']][:3] == [31.0, 41.0, 51.0] and r['dip2'][-1][0] == 60.0
    assert r['inter']['n'] == 3 and r['inter']['dW'] < 0 and r['t_land'] == 61.0
    p62 = MWE.pair_p(sts[-2][1], sts[-2][2], 0, 1)
    assert abs(p62 - 3 * 0.01 / (1.0 * 0.06)) < 1e-12 and 0.28 < p62 < 1
    assert abs(MWE.pair_p(sts[-1][1], sts[-1][2], 0, 1) - 0.3) < 1e-9
    # counted: τ = 61 (landed, p small: agrees), 62 (landed, p 0.5: disagrees), 63 (not landed, p 0.3: agrees); p ≥ 1 elsewhere
    assert r['agree'] == (2, 3) and r['agree4'] == (1, 2)


def test_gl_job_and_its_constant(monkeypatch):
    """Referee rounds 74, 80 (GJ1, GJ3, GJ4, GJd): `gl_job` reads the loop's states (`loop_nd` monkeypatched) — the dipole pair's Δc,
    w_a, w_b, mean width (of both, at unequal widths) and nearest amplitude, E of all kernels and of the formed ones, Λ, Gm, μ₀ —;
    `gl_constant` on rows on the law
    Δc⁻⁶ = Δc_L⁻⁶ + 48aGm(W²/(w_a² + w_b²))Λτ is 1 (C = 0 at W = A) and its C is Gm W²(W − A)²(1/|w_a| + 1/|w_b|)/(μ₀Δc²(w_a² + w_b²))."""
    cs, ss, As = [0.0], [0.06], [1.0]
    taus = (400.0, 800.0)
    monkeypatch.setattr(MWE, 'TAUS_GL', taus)
    dcs = (2e-3, 1e-3)
    Wt = [np.zeros(3)] + [np.array([3.0, -2.0, 1e-4]) for _ in taus]
    Ct = [np.zeros((3, 2))] + [np.array([[-0.4 * d, 0.05], [0.6 * d, 0.07], [0.0, 0.06]]) for d in dcs]   # the idle kernel on the feature
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a, **k: (None, Wt, Ct, 77.0, None))
    g = MWE.gl_job((5, cs, ss, As))
    K, L = MW.BumpsCW(MWE.XS, MWE.S0), MWE.law()
    f = MWE.feature_sum(cs, ss, As)
    x = g['rows'][1]
    assert g['Lam'] == 77.0 and g['Gm'] == L.G * L.m and g['mu0'] == 4 * L.G * L.c / L.A and g['q'] == 5
    assert abs(x['dc'] - 1e-3) < 1e-15 and (x['wa'], x['wb'], x['A']) == (3.0, -2.0, 1.0) and abs(x['sbar'] - 0.06) < 1e-15
    ef = np.mean((f - K.synth(Ct[2][:2], Wt[2][:2])) ** 2) / 2
    assert np.isclose(x['E'], np.mean((f - K.synth(Ct[2], Wt[2])) ** 2) / 2, rtol=1e-12, atol=0)
    assert abs(x['E_formed'] / ef - 1) < 1e-12 and abs(x['E'] / ef - 1) > 1e-6   # referee round 83 (GJe): the idle kernel left out
    Gm, mu0, Lam, a = 12.0, 0.1, 50.0, 3e6
    wa, wb = 3.0, -2.0
    R = 48 * a * Gm * 1.0 / (wa ** 2 + wb ** 2) * Lam
    row = lambda t, dW=0.0: dict(tau=t, dc=(1e14 + R * t) ** (-1 / 6), wa=wa, wb=wb, sbar=0.06, A=1.0, dW=dW,
                                 E=a * (1e14 + R * t) ** (-8 / 6))
    p, r0, r1, c = MWE.gl_constant(row(100.0), row(200.0), Gm, mu0, Lam)
    assert abs(r0 - 1) < 1e-9 and abs(r1 - 1) < 1e-9 and c == 0.0
    x1 = row(200.0, dW=-0.01)
    p, r0, r1, c = MWE.gl_constant(row(100.0, dW=-0.01), x1, Gm, mu0, Lam)
    x0 = row(100.0, dW=-0.01)
    cc = lambda x: Gm * 1.0 * 0.01 ** 2 * (1 / 3 + 1 / 2) / (mu0 * x['dc'] ** 2 * 13.0)
    rate = lambda x: 48 * x['E'] / x['dc'] ** 8 * Gm * 1.0 / 13.0 * Lam
    meas = (x1['dc'] ** -6 - x0['dc'] ** -6) / 100.0
    assert abs(c - cc(x1)) < 1e-9 * c and abs(r0 - 1) < 1e-9
    assert abs(r1 - meas / (0.5 * (rate(x0) / (1 + cc(x0)) + rate(x1) / (1 + cc(x1))))) < 1e-9 * r1 and r1 > r0
    assert abs(p - 5.0 * x1['dc'] / 0.06) < 1e-12


def test_a2_one_kernel_from_large_coefficients(monkeypatch):
    """Referee rounds 74, 80 (D3bu-5, D3ca-3): from (c*, σ*) with w = 3A (E(start) four times ½⟨f*²⟩, outside (i)'s condition) or 10A
    the one-kernel block collapses its width onto the samples at once — J_N's σ → 0 limit, E/½⟨f*²⟩ = 0.906 at σ*/σ₀ = 1 whatever w —,
    stays there (3A to s ≈ 3·10⁵, 10A to 10¹⁰) and escapes to the exact fit (10A by 3·10¹⁰); at σ*/σ₀ = 0.25 the end is σ/σ* = 1 (A2b);
    the flow's tolerances (A2d); A2's table and the spike's summary."""
    r3, r10, rq = MWE.a2_job((1.0, 3.0)), MWE.a2_job((1.0, 10.0)), MWE.a2_job((0.25, 2.5))
    for r in (r3, r10, rq):
        e = r['rows'][-1]
        assert abs(e['w'] - 1) < 1e-8 and abs(e['sig'] - 1) < 1e-8 and e['E'] < 1e-20 and r['status'] == 0
    assert r3['E0'] == 4.0 and r10['E0'] == 81.0
    e2 = MWE.a2_plateau(1.0)[0]
    sp3, sp10 = MWE.a2_spike(r3['rows'], e2), MWE.a2_spike(r10['rows'], e2)
    assert sp3[0] == 1.0 and np.isclose(sp3[1], 10 ** 5.5, rtol=1e-12) and np.isclose(sp3[3], 1e6, rtol=1e-12)
    assert sp10[0] == 1.0 and np.isclose(sp10[1], 1e10, rtol=1e-12) and np.isclose(sp10[3], 10 ** 10.5, rtol=1e-12)
    assert MWE.a2_spike(rq['rows'], MWE.a2_plateau(0.25)[0]) is None
    assert abs(sp10[2][0] - 0.906) < 1e-3 and abs(sp10[2][1] - sp3[2][1]) < 1e-4    # the plateau's E is J_N's limit, for every w
    assert max(abs(sp3[2][0] - e2), abs(sp3[2][1] - e2)) < 1e-4 and max(abs(sp10[2][0] - e2), abs(sp10[2][1] - e2)) < 1e-8   # E₂
    d3 = MWE.a2_job((0.25, 3.0))                                        # referee round 87 (AS6 inverted): off the sample 3A at 0.25 collapses
    sp_ = MWE.a2_spike(d3['rows'], MWE.a2_plateau(0.25)[0])
    assert min(x['sig'] for x in d3['rows']) < 0.12 and sp_ is not None and sp_[0] == 1.0 and d3['off'] == MWE.A2_OFF == 1e-9
    md = '\n'.join(MWE.a2_lines([r3, r10]))
    e3 = r3['rows'][-1]
    assert (f"| 1 | 3 | 4 | {e3['w']:.4f}, {e3['sig']:.4g}, {e3['dc']:+.1e}, {e3['E']:.2e} (s = 1e+13) | "
            f"{min(x['sig'] for x in r3['rows']):.3g} | 1e+00 … 3e+05, {sp3[2][0]:.3f}–{sp3[2][1]:.3f} | 1e+06 | 0 |") in md
    assert MWE.a2_lines(None) == [] and MWE.a2_lines([]) == []
    pl = MWE.a2_plateau(1.0)
    assert f"E₂ {MWE.a2_plateau(0.25)[0]:.6f}, {pl[0]:.6f}, {MWE.a2_plateau(4.0)[0]:.6f}, E₁" in md
    import scipy.integrate
    seen = {}
    real = scipy.integrate.solve_ivp
    monkeypatch.setattr(scipy.integrate, 'solve_ivp', lambda *a, **k: seen.update(k) or real(*a, **k))
    MWE.a2_job((1.0, 2.5))
    assert (seen['rtol'], seen['atol'], seen['method']) == (1e-10, 1e-14, 'LSODA')


def test_a2_plateau_is_the_two_sample_limit():
    """Referee round 83 (D3cd-4): E₂ is the two samples nearest c* interpolated — J_N minimised over (c, w) at a width an eighth of the
    samples' spacing equals it —, E₁ the one-sample limit at c = c* (a sample); both against direct sums."""
    from scipy.optimize import minimize_scalar
    for ratio in MWE.RATIOS_A2:
        f = MWE.feature_sum([MWE.C_A], [ratio * MWE.S0], [MWE.A_A])
        q = f ** 2 / np.sum(f ** 2)
        j = int(np.argmin(np.abs(MWE.XS - MWE.C_A)))
        e2, e1 = MWE.a2_plateau(ratio)
        assert abs(MWE.XS[j] - MWE.C_A) < 1e-12 and abs(e1 - (1 - q[j])) < 1e-14 and abs(e2 - (1 - q[j] - q[j + 1])) < 1e-14
        h = MWE.XS[1] - MWE.XS[0]
        J = lambda c: 1 - np.dot(f, np.exp(-(MWE.XS - c) ** 2 / (2 * (h / 8) ** 2))) ** 2 / (
            np.sum(np.exp(-(MWE.XS - c) ** 2 / (h / 8) ** 2)) * np.sum(f ** 2))
        m = minimize_scalar(J, bounds=(MWE.XS[j], MWE.XS[j + 1]), method='bounded', options=dict(xatol=1e-12))
        assert abs(m.fun - e2) < 1e-9 and J(MWE.C_A) > e1 - 1e-12


def test_a2e_the_early_collapse(monkeypatch):
    """Referee round 83 (D3cd-4): A2E — over A2's span its outputs at s = 1 and 10 are `a2_job`'s (the terminal event past s = 10 leaves
    the steps alone); the span, the event and the table's first plateau output."""
    r = MWE.a2e_job(1e13)
    a = MWE.a2_job((0.25, 5.0))
    ro = MWE.a2e_job(1e13, off=MWE.A2_OFF)                               # referee round 87: A2 starts off the sample, A2E on it
    ros = {x['s']: x for x in ro['rows']}
    for k_ in (0, 2):
        assert {k: v for k, v in ros[MWE.S_A2[k_]].items() if k != 'jcc'} == a['rows'][k_]
    assert r['off'] == 0.0 and all(len(x) == 5 and 1 <= x[4] <= 12 for x in r['steps']) and all(x[4] <= 5 for x in r['steps'] if x[2] == 2)
    rs = {x['s']: x for x in r['rows']}
    assert r['status'] == 1 and len(r['rows']) == len(MWE.S_A2E) and r['span'] == 1e13
    assert abs(rs[1.0]['dc']) > 0.1 and abs(ros[1.0]['dc']) > 0.1                 # 5A leaves the symmetric state either way
    r10 = MWE.a2e_job(10.0)
    assert r10['status'] == 0 and len(r10['rows']) == len(MWE.S_A2E)
    r3 = {x['s']: x for x in r10['rows']}                               # round 83's AE3: ∂ċ/∂c, not ∂σ̇/∂σ — the saddle's sign
    assert r3[1e-3]['jcc'] > 1e4 and r3[1e-4]['jcc'] < 0 and abs(r3[1e-3]['dc']) < 1e-12
    for r_ in (r, r10):                                                 # round 83's AE4, AE5, AE7: LSODA's own steps
        st = r_['steps']
        assert st and all(MWE.S_A2E_STEPS[0] < x[0] for x in st) and st[-1][0] >= MWE.S_A2E_STEPS[1] > st[-2][0]
        assert {x[2] for x in st} <= {1, 2} and all(abs(b[0] - a[0] - b[1]) < 1e-15 for a, b in zip(st[:-1], st[1:]))
    on = lambda r_: {x[2] for x in r_['steps'] if abs(x[3]) < 1e-6}    # on the symmetric state BDF on the span (0, 10) only
    assert on(r10) == {1, 2} and on(r) == {1}
    assert [x[4] for x in r10['steps'] if abs(x[3]) < 1e-6 and x[2] == 2] == [2, 2, 2, 1, 1, 1, 1, 1, 1]   # round 100: NQU on the job (AE2, AE3)
    e2 = MWE.a2_plateau(0.25)[0]
    md = '\n'.join(MWE.a2_lines([a], [r, r10]))
    first = next(x['s'] for x in r['rows'] if abs(x['E'] - e2) < MWE.PLATEAU_A2)
    x1 = rs[1.0]
    assert f"| 1e+13 | " in md and f"{x1['dc']:+.1e}, {x1['sig']:.4f}, {x1['E']:.3f} |" in md and f" | {first:.2g} |" in md
    near = dict(r, rows=[dict(x, E=e2 + 0.05) if x['s'] < first else x for x in r['rows']])   # round 83's AL6: near E₂ is not on it
    assert f" | {first:.2g} |" in '\n'.join(MWE.a2_lines([a], [near]))
    import scipy.integrate
    seen = {}
    real = scipy.integrate.solve_ivp
    monkeypatch.setattr(scipy.integrate, 'solve_ivp', lambda *a_, **k: seen.update(k, span=a_[1]) or real(*a_, **k))
    MWE.a2e_job(10.0)
    assert (seen['rtol'], seen['atol'], seen['method'], seen['span']) == (1e-10, 1e-14, 'LSODA', (0.0, 10.0))
    assert seen['events'](10.0, None) < 0 < seen['events'](11.0, None) and seen['events'].terminal and tuple(seen['t_eval']) == MWE.S_A2E


def test_loop_row_counts_visible_kernels_and_hll_job(monkeypatch):
    """Referee round 74 (R8, D3bu-2): HL's row counts the visible formed kernels of the other sign and the invisible formed ones;
    `hll_job` keeps the loop's kernels above 1% of the smallest amplitude (one at 0.5% is not formed)."""
    x = dict(tau=1e5, theta=[[-0.08, 0.06], [0.005, 3e-4], [0.08, 0.06], [2.0, 0.06]], w=[1.2, -0.7, -0.3, 0.5], wmax=1.2, rms=2e-6)
    assert MWE.loop_row('t', 1e-8, x) == "| t | 1e-08 | 100000 | 4 | 1 (-0.300) | 1 | 0.060 | 1 | 1.2 | 2.0e-06 |"
    monkeypatch.setattr(MWE, 'TAUS_HLL', (10.0,))
    Wt = [np.zeros(3), np.array([1.0, -0.005, 0.02])]
    Ct = [np.zeros((3, 2)), np.array([[-0.08, 0.06], [0.0, 0.06], [0.08, 0.06]])]
    got, kw = [], []
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a, **k: got.append(a[4]) or kw.append((a[5], k)) or (None, Wt, Ct, 1.0, None))
    r = MWE.hll_job(('σ* 0.06, equal', 1e-6))
    assert r['loops'][0]['w'] == [1.0, 0.02] and r['loops'][0]['wmax'] == 1.0 and r['eps'] == 1e-6 and got == [1e-6]   # HLLb
    r8 = MWE.hll_job(('σ* 0.06, A₂ 0.99', 1e-8))                        # referee round 83 (HJ1–HJ3): the cell's mass, times, tolerances
    assert got[-1] == 1e-8 and r8['eps'] == 1e-8 and list(kw[-1][0]) == [0.0, 10.0] and kw[-1][1] == dict(rtol=1e-9, atol=1e-12)
    md = MWE.report(HLL=[dict(r8, loops=[dict(x, **{k_: None for k_ in ('main_err', 'n_rest')}) for x in r8['loops']])])
    sp8 = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06]], [1.0, 0.02], (-0.08, 0.08), (0.06, 0.06), (1.0, 0.99))
    assert MWE.split_cells(sp8) + ' |' in md and "The others' sum on the samples (referee round 83)" in md   # the split from θ, w
    sp = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06]], [1.0, 0.02], (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0))
    assert {k_: r['loops'][0][k_] for k_ in sp} == sp                   # referee round 80: hll_job's split is `fit_split`'s
    monkeypatch.setattr(MWE, 'EPS_HL', (1e-7,))                         # HLa: hl_job's formed set at 1%, not 0.1%
    monkeypatch.setattr(MWE, 'TAUS_HL', (10.0,))
    got.clear()
    q = MWE.hl_job(('t', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)))
    assert q['loops'][0]['w'] == [1.0, 0.02] and got == [1e-7] and q['loops'][0]['neg'] == 0


def test_fit_split_tells_a_null_cluster_from_a_cancellation():
    """Referee round 80 (D3ca-2): `fit_split` — per feature the nearest formed kernel of its sign (each used once) and its largest
    relative error; the others' count, Σ|w|, sum on the samples over max f*, widths and largest |c| —; `split_cells` and `end_split`."""
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    null = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [0.0, 1.5], [1.9, 0.3]], [1.0, 1.0, 0.4, -0.4], cs, ss, As)
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    foot = np.abs(K.synth(np.array([[0.0, 1.5], [1.9, 0.3]]), np.array([0.4, -0.4]))).max() / np.abs(MWE.feature_sum(cs, ss, As)).max()
    assert null['main_err'] == [0.0, 0.0] and null['n_rest'] == 2 and null['rest_w'] == 0.8 and abs(null['rest_foot'] - foot) < 1e-15
    assert null['rest_sig'] == [0.3, 1.5] and null['rest_c'] == 1.9 and null['main_dw'] == [0.0, 0.0] and null['rest_sw'] == 0.0
    canc = MWE.fit_split([[-0.08, 0.06], [0.0, 0.05], [0.0, 0.04]], [1.0, 3.0, -2.4], cs, ss, As)
    assert canc['main_err'] == [0.0, 2.0] and canc['n_rest'] == 1 and canc['rest_foot'] > 1   # feature 2's: the +3 kernel at 0
    assert canc['main_dw'] == [0.0, 2.0] and canc['rest_sw'] == -2.4    # referee round 83: inflated, the other's of the other sign
    th = [[-0.042, 0.0438], [-0.123, 0.06], [0.08, 0.06]]               # referee round 83 (D3cd-2a): HL's row — nearer in c, worse
    near = MWE.fit_split(th, [0.0706, 1.58, 1.0], cs, ss, As)
    assert abs(near['main_err'][0] - 0.043 / 0.06) < 1e-12 and abs(near['main_dw'][0] - 0.58) < 1e-12 and near['rest_w'] == 0.0706
    tie = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [0.08, 0.06]], [2.0, 2.0, 1.0], cs, ss, As)   # the largest errors tie: the next
    assert tie['main_err'] == [1.0, 0.0] and tie['rest_sw'] == 2.0
    spl = MWE.fit_split([[-0.08, 0.06], [-0.08, 0.06], [0.08, 0.06]], [0.6, 0.4, 1.0], cs, ss, As)   # a split: deficient, the rest
    assert abs(spl['main_dw'][0] + 0.4) < 1e-12 and spl['rest_sw'] == 0.4
    mm = MWE.fit_split([[0.0, 0.06], [0.0, 0.06]], [1.1, 0.2], (0.0, 0.0), (0.06, 0.06), (1.0, 2.0))   # FT1: the largest error first —
    assert np.allclose(mm['main_err'], [0.8, 0.45], rtol=0, atol=1e-12)  # (0.8, 0.45) over (0.1, 0.9), which the smallest would take
    assert MWE.fit_split([[-0.08, -0.06], [0.08, 0.06]], [1.0, 1.0], cs, ss, As)['main_err'] == [0.0, 0.0]   # FT7: |σ|
    assert MWE.fit_split([[0.0, 0.06]], [-1.0], (0.0,), (0.06,), (1.0,))['main_err'] == [None]
    once = MWE.fit_split([[0.0, 0.06], [0.15, 0.06]], [1.0, 1.0], cs, ss, As)    # each kernel main to one feature only (FS1; within SPLIT_CS)
    assert once['main_err'][0] == abs(0.0 + 0.08) / 0.06 and abs(once['main_err'][1] - abs(0.15 - 0.08) / 0.06) < 1e-12
    assert MWE.split_cells(null) == (f"0.0e+00 (+0), 0.0e+00 (+0) | 2, +0.000 (0.800), mass {np.sqrt(2 * np.pi) * (0.4 * 1.5 - 0.4 * 0.3):+.3f} "
                                     f"(window {null['rest_win']:+.1e}), {foot:.1e} (rms {null['rest_rms']:.1e}), 0.300–1.500, 1.90")
    assert MWE.split_cells(canc).startswith("0.0e+00 (+0), 2.0e+00 (+2) | 1, -2.400 (2.400), ")
    r = dict(moving=dict(theta=[[-0.08, 0.06], [0.08, 0.06], [0.5, 0.06]], w=[1.0, 1.0, 0.005]), cs=cs, ss=ss, As=As)
    assert MWE.end_split(r)['n_rest'] == 0                              # the unformed kernel (0.5%) left out


def test_g_fine_lines_and_gl_lines_on_synthetic_results():
    """Referee rounds 74, 80: the fine-grid table evaluates the dipole law at the time of v(τ₃⁻) and the intermediate phase's slopes over
    3μ₀K′Λ, κ and κ′ (v carried by the law to the path's τ₃), † for another sequence; the path's rows (GT: K, K′/K with the closure's
    rms, τ₃, τ₄, the dissolution time, κ, the landing's v by the law, v after τ₄; an exact closure); κ′ against the mass, the grid's
    offsets and the agreement's counts (all, ε ≤ 10⁻⁷, past τ₄); GL's constant table, its ranges and the local slope."""
    Lam, mu0, K_, Kp = 1e3, 0.1, 2e-3, -4e-4
    rate = 3 * mu0 * K_ * Lam
    t_d, t3, t4 = 10.0, 20.0, 30.0
    dip = [(t, (rate * (t - t_d)) ** (1 / 3)) for t in (11.0, 15.0, 19.9)]
    v3 = dip[-1][1]
    dip2 = [(t, (v3 ** 3 + 3 * mu0 * Kp * Lam * (t - t3)) ** (1 / 3)) for t in (20.0, 25.0, 29.9)]
    ph = dict(t3=t3, t4=t4, v3=v3, t_v3=19.9, dip=dip, dip2=dip2, inter=dict(n=3, off=2.5, dW=-0.4), t_land=31.0, v_land=0.97 * v3,
              p_land=0.1, agree=(5, 6), agree4=(2, 3))
    th = [[0.0, 0.06], [0.001, 0.06], [0.5, 0.06]]
    loops = [dict(eps=1e-10, tau=6400.0, Lam=Lam, theta=th, w=[0.97 * v3 + 0.5, -0.97 * v3 + 0.5, 0.8], wmax=1.0, rms=1e-9)]
    lam6 = MWE.RT.transit_scale(MWE.law(), 1e-6)
    dip2x = [(t, (v3 ** 3 - 3 * 3 * mu0 * Kp * lam6 * (t - t3)) ** (1 / 3)) for t in (20.0, 25.0, 29.9)]   # −3× the law: another sequence
    G = [dict(q=4, dmin=0.4, cs=[0.0, 0.5], ss=[0.06, 0.06], As=[1.0, 0.8], entries=2, tau_d=t_d, status=0, loops=loops,
              law=[dict(v=10.0, E=K_ / 100.0 + 0.5, E_inf=0.5, mu0=mu0)], fine={1e-6: dict(ph, agree=(1, 4), dip2=dip2x), 1e-10: ph})]
    G.append(dict(G[0], q=5, fine={1e-10: ph}))
    K3 = K_ * (1 + 1e-5)
    td_ = 19.95 + K3 * (19.95 - t_d) / abs(Kp)
    GT = [dict(q=4, t3_pred=19.95, K_third=K3, Kp=Kp, t3_path=19.95, t4_path=29.95, kappa_path=0.96, t_diss=td_, rms2=1e-3, exact=False),
          dict(q=5, t3_pred=19.95, K_third=K3, Kp=0.0, t3_path=19.95, t4_path=None, kappa_path=1.0, t_diss=None, rms2=4e-17, exact=True)]
    md = '\n'.join(MWE.g_fine_lines(G, GT))
    kp_ = 0.97 / ((19.95 - t_d) / (19.9 - t_d)) ** (1 / 3)
    lo = 0.97 / ((20.0 - t_d) / (19.9 - t_d)) ** (1 / 3)                 # referee round 83: the loop's own crossing's bracket
    assert f"| 4 | 1e-10 | 10 | 20 | 30 | 1.0000 | 1.0000 | 3, 2.5, -0.400, 1.0000 | 31, 0.970, {kp_:.3f}, 0.100 |" in md
    assert f"| 4 | 1e-06 | 10 | 20 | 30 | {(1.0 * Lam / lam6):.4f} |" in md
    assert f"| 3, 2.5, -0.400, -3.0000 | 31, 0.970, {kp_:.3f}†, 0.100 |" in md
    vl = (3 * mu0 * Lam * (K3 * (19.95 - t_d) + Kp * 10.0)) ** (1 / 3)
    vv = (0.97 * v3 + 0.5 - (-0.97 * v3 + 0.5)) / 2
    assert (f"| 4 | {K3:.5e} ({K_:.5e}) | 19.95 (20) | {Kp / K3:+.4f} (1.0e-03) | 29.95 (30) | {td_:.5g} | 0.9600 ({kp_:.4f}; "
            f"[{lo:.4f}, 0.9700]) | "
            f"{vl:.4g} ({0.97 * v3:.4g}) | {vv / (0.97 * v3) - 1:+.1e} |") in md
    vl5 = (3 * mu0 * Lam * K3 * (19.95 - t_d)) ** (1 / 3)
    assert f"| 5 | {K3:.5e} ({K_:.5e}) | 19.95 (20) | 0, exact (4.0e-17) | — (exact) (30) | — | 1.0000 ({kp_:.4f}; [{lo:.4f}, 0.9700]) | {vl5:.4g} (" in md
    assert f"κ′ against the mass (ε = 1e-06, 1e-10; † another sequence): 4: {kp_:.3f}†, {kp_:.3f}; 5: {kp_:.3f}." in md
    off = ((19.95 - t_d) / (19.9 - t_d)) ** (1 / 3) - 1
    assert f"that κ′ removes, at ε = 1e-10: 4: {off:+.2%}, 5: {off:+.2%}." in md
    assert "agree in 11 of 16 (10 of 12 at ε ≤ 10⁻⁷; past τ₄ 6 of 9)." in md
    assert "κ′ − κ (the path's)" not in md and "by index (SQ" not in md
    assert (f"own crossing's bracket at the smallest mass (referee round 83): inside on none; outside on 4 (the bracket {lo / 0.96 - 1:+.2%} … "
            f"{0.97 / 0.96 - 1:+.2%} off it), 5 (the bracket {lo - 1:+.2%} … {0.97 - 1:+.2%} off it).") in md
    # referee round 83: SQ's direct test — 1e-10 on the path's sequence, 1e-6 not (kernel 5 at τ₃) —, † on '—', the rate
    ev = lambda j3: [(10.0, [0, 1], []), (20.0, [j3], []), (30.0, [3], [])]
    SQ = [dict(q=4, seq={1e-6: ev(5), 1e-10: ev(2)})]
    GT[0].update(block=[0, 1], j3=2, j4=3)
    G2 = [dict(G[0], fine={1e-6: dict(ph, agree=(1, 4), dip2=dip2x, v_land=None), 1e-10: ph}), G[1]]
    md = '\n'.join(MWE.g_fine_lines(G2, GT, SQ))
    assert "† another sequence — the loop's formed kernels at its τ₃ or τ₄ not the path's, by index (SQ, referee round 83) —" in md
    assert f"κ′ against the mass (ε = 1e-06, 1e-10; † another sequence): 4: —†, {kp_:.3f}; 5: {kp_:.3f}." in md
    assert ("the path's 4: {0, 1} + 2 + 3. The loop follows it at 1 of 2 cells (SQ's τ₃ and τ₄ are G's at 2 of 2); elsewhere: 4 at "
            "1e-06: τ₃ 20 {0, 1, 5}, τ₄ 30 {0, 1, 3, 5}. The intermediate slope's proxy (beyond 1 ± 0.5 times the path's law, round 80) "
            "marks the same cells.") in md
    GTn = [dict(GT[0], block=None), GT[1]]                              # GF9: without GT's block SQ is not used — the proxy's †
    md = '\n'.join(MWE.g_fine_lines(G2, GTn, SQ))
    assert f"κ′ against the mass (ε = 1e-06, 1e-10; † another sequence): 4: —†, {kp_:.3f}; 5:" in md and "by index (SQ" not in md
    SQ6 = [dict(q=4, seq={1e-6: ev(5), 1e-10: [(10.0, [0, 1], []), (20.0, [2], []), (35.0, [3], [])]})]   # GF6: τ₄ at 35, not G's 30
    assert "(SQ's τ₃ and τ₄ are G's at 1 of 2)" in '\n'.join(MWE.g_fine_lines(G2, GT, SQ6))
    SQ[0]['seq'][1e-6] = ev(2)                                          # on the path's sequence by index, the proxy's † differs
    md = '\n'.join(MWE.g_fine_lines(G2, GT, SQ))
    assert "follows it at 2 of 2 cells" in md and "differs at 4 at 1e-06." in md and "4: —, " in md
    G3 = [dict(G[0], q=6, fine={e: dict(ph, v_land=0.97 * v3 * k_) for e, k_ in ((1e-8, 1.03), (1e-9, 1.01), (1e-10, 1.002))},
               loops=[dict(loops[0], eps=e) for e in (1e-8, 1e-9, 1e-10)])]
    GT3 = [dict(GT[0], q=6, kappa_path=kp_)]
    md = '\n'.join(MWE.g_fine_lines(G3, GT3))
    d8, d9, d10 = (kp_ * k_ - kp_ for k_ in (1.03, 1.01, 1.002))
    assert (f"6: {d8:+.4f}, {d9:+.4f}, {d10:+.4f} (×{d8 / d9:.2f}–{d9 / d10:.2f}, as ε^{{{np.log10(d8 / d9):.2f}–"
            f"{np.log10(d9 / d10):.2f}}}) — over ε^{{1/3}}: {d8 / 1e-8 ** (1 / 3):+.1f}, {d9 / 1e-9 ** (1 / 3):+.1f}, {d10 / 1e-10 ** (1 / 3):+.1f}." in md
            and d10 >= 1e-3)                                            # referee round 87 (D3ch-7): over ε^{1/3}
    G4 = [dict(G3[0], fine={e: G3[0]['fine'][e] for e in (1e-8, 1e-10)})]   # GF4: two decades apart, no ratio
    assert f"6: {d8:+.4f}, {d10:+.4f} — over ε^{{1/3}}: " in '\n'.join(MWE.g_fine_lines(G4, GT3))
    G3[0]['fine'][1e-10]['v_land'] = 0.97 * v3 * 1.0005                 # |κ′ − κ| below 10⁻³: no ratio with it
    md = '\n'.join(MWE.g_fine_lines(G3, GT3))
    assert f"(×{d8 / d9:.2f}, as ε^{{{np.log10(d8 / d9):.2f}}}) — over ε^{{1/3}}: " in md
    rows = [dict(tau=t, d=1e-3 * t ** (-1 / 6), dW=-0.01 * t ** (-1 / 3), rms=1e-6 * t ** (-2 / 3), E=1e-12, E_formed=0.9e-12,
                 dc=1e-3 * t ** (-1 / 6), wa=3.0, wb=-2.0, sbar=0.06, A=1.0) for t in (400.0, 800.0, 1600.0)]
    md = '\n'.join(MWE.gl_lines([dict(q=9, rows=rows, Lam=10.0, Gm=12.0, mu0=0.1)]))
    cons = [MWE.gl_constant(a, b, 12.0, 0.1, 10.0) for a, b in zip(rows[:-1], rows[1:])]
    assert "### GL. The collapse (referee rounds 68, 74): target 9" in md
    assert f"| 800 | {cons[0][0]:.4f} | {cons[0][1]:.4f} | {cons[0][2]:.4f} | {cons[0][3]:.4f} |" in md
    assert (f"Over the 2 intervals the ratio is {min(c[1] for c in cons):.4f}–{max(c[1] for c in cons):.4f}, with C "
            f"{min(c[2] for c in cons):.4f}–{max(c[2] for c in cons):.4f}; the formed kernels' energy over all kernels' 0.900000–0.900000.") in md
    assert "The remainder with C against p" not in md                   # p's range below 0.05: no fit
    # referee rounds 74, 80: E = a₈Δc⁸(1 + kp) along the run — the local slope (d ln E/d ln Δc − 8)/p of the formed kernels' E, the
    # force factor's remainder, the mean width's share p²/36
    k8, sb, wa, wb = 0.8, 0.06, 3.0, -2.0
    pr = lambda dc: abs(wa - wb) * dc / (abs(wa + wb) * sb)
    rows = [dict(tau=t, d=dc, dW=-0.01, rms=1e-6, E=1e3 * dc ** 8 * (1 + k8 * pr(dc)) + 1e-23, E_formed=1e3 * dc ** 8 * (1 + k8 * pr(dc)),
                 dc=dc, wa=wa, wb=wb, sbar=sb, A=1.0) for t, dc in ((400.0, 3e-3), (800.0, 2e-3), (1600.0, 1e-3))]   # E: an idle floor
    for x in rows:                                                      # referee round 87 (GL1): the floor at 1e-21, ≥ 1% of E at the last row (asserted)
        x['E'] = x['E_formed'] + 1e-21
    assert 1e-21 >= 0.01 * max(x['E_formed'] for x in rows[-1:])
    md = '\n'.join(MWE.gl_lines([dict(q=9, rows=rows, Lam=10.0, Gm=12.0, mu0=0.1)]))
    sl = [np.log(b['E_formed'] / a['E_formed']) / np.log(b['dc'] / a['dc']) for a, b in zip(rows[:-1], rows[1:])]
    ks = [(x - 8) / (0.5 * (pr(a['dc']) + pr(b['dc']))) for x, a, b in zip(sl, rows[:-1], rows[1:])]
    assert all(abs(k_ / k8 - 1) < 0.25 for k_ in ks)                   # the local slope sees E's next order, kp/(1 + kp)
    cons = [MWE.gl_constant(a, b, 12.0, 0.1, 10.0) for a, b in zip(rows[:-1], rows[1:])]
    st = [np.log(b['E'] / a['E']) / np.log(b['dc'] / a['dc']) for a, b in zip(rows[:-1], rows[1:])]   # referee round 83: all kernels'
    kt = [(x - 8) / (0.5 * (pr(a['dc']) + pr(b['dc']))) for x, a, b in zip(st, rows[:-1], rows[1:])]
    rem = [c[2] / (x / 8) - 1 for c, x in zip(cons, st)]
    assert f"(p their mean) {min(ks):.3f}–{max(ks):.3f}, median {np.median(ks):.3f}" in md and f"p²/36, at most {pr(3e-3) ** 2 / 36:.5f}." in md
    assert f"(all kernels' energy: {min(kt):.3f}–{max(kt):.3f}, the idle kernels' floor in it)" in md and min(kt) != min(ks)
    assert f"the ratio with C is 1 {min(rem):+.1e} … {max(rem):+.1e}" in md


def test_sq_job_and_sq_check(monkeypatch):
    """Referee round 83 (D3cd-6): `sq_job` reruns G's loop at each mass on G's outputs and records, where the formed set (above 1% of the
    smallest amplitude) changes, the kernels formed and unformed; `sq_check` finds the loop's τ₃ and τ₄ as `fine_phases` does and
    compares the formed sets there with the path's block and entrants (an exact closure: no τ₄)."""
    monkeypatch.setattr(MWE, 'TAUS_G', (10.0, 20.0))
    monkeypatch.setattr(MWE, 'TAUS_FINE', (10.0, 15.0))
    monkeypatch.setattr(MWE, 'EPS_G', (1e-6,))
    monkeypatch.setattr(MWE, 'EPS_GX', (1e-12,))
    Wt = [np.zeros(4), np.array([1.0, 0.0, 0.0, 0.0]), np.array([1.0, -0.5, 0.005, 0.0]), np.array([0.0, -0.5, 0.02, 0.0])]
    got = []
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a, **k: got.append((a[4], list(a[5]), k)) or (None, Wt, None, None, None))
    r = MWE.sq_job((7, [0.0], [0.06], [1.0]))
    ev = [(10.0, [0], []), (15.0, [1], []), (20.0, [2], [0])]           # 0.005 is not formed, 0.02 is
    assert r == dict(q=7, seq={1e-6: ev, 1e-12: ev})
    assert [g[0] for g in got] == [1e-6, 1e-12] and got[0][1] == [0.0, 10.0, 15.0, 20.0] and got[0][2] == dict(rtol=1e-9, atol=1e-12)
    g = dict(block=[0, 1], j3=2, j4=3, exact=False)
    path = [(10.0, [0, 1], []), (15.0, [2], []), (20.0, [3], [])]
    assert MWE.sq_check(path, 12.0, 2, g) == (15.0, 20.0, [0, 1, 2], [0, 1, 2, 3], True)
    assert MWE.sq_check(path, 12.0, 2, dict(g, j4=5))[4] is False        # another τ₄ entrant
    early = [(10.0, [0, 1, 5], []), (20.0, [2], [])]                    # a third kernel before τ_d: τ₃ the first output past it
    assert MWE.sq_check(early, 12.0, 2, g) == (15.0, 20.0, [0, 1, 5], [0, 1, 2, 5], False)
    ex = dict(block=[0, 1], j3=2, j4=None, exact=True)
    assert MWE.sq_check(path[:2], 12.0, 2, ex)[4] is True and MWE.sq_check(path, 12.0, 2, ex)[4] is False
    assert MWE.sq_check(path[:2], 12.0, 2, dict(ex, exact=False))[4] is False
    assert MWE.sq_check(path, 12.0, 2, dict(g, block=None))[4] is False
    assert MWE.sq_check([(10.0, [0, 1], [])], 12.0, 2, g) == (None, None, None, None, False)
    gone = [(10.0, [0, 1, 5], []), (12.5, [], [5]), (15.0, [2], []), (20.0, [3], [])]   # SK6: kernel 5 unformed before τ₃
    assert MWE.sq_check(gone, 12.0, 2, g) == (15.0, 20.0, [0, 1, 2], [0, 1, 2, 3], True)
    two = [(10.0, [0, 1], []), (15.0, [2, 3], []), (17.0, [], [3]), (20.0, [3], [])]   # SK2: two at τ₃, one leaves and comes back
    assert MWE.sq_check(two, 12.0, 2, g)[:2] == (15.0, None)
    assert MWE.sq_check(path, 15.0, 2, g)[:2] == (20.0, None)             # SK1: at an output equal to τ_d, τ₃ is the next one's
    assert [a[0] for a in MWE.sq_args([dict(q=q, status=st, cs=[0.0], ss=[0.1], As=[1.0]) for q, st in ((1, 0), (2, 1), (3, 0))])] == [1, 3]


def test_round87_split_mass_bdf_and_the_md():
    """Referee round 87 (D3ch-2, D3ch-3, D3ch-5): `fit_split` assigns a feature no main kernel beyond the cap and reports the others'
    signed mass √(2π)Σw|σ| and their sum's rms on the samples over max f*; `bdf_bound` is BDF-k's damping threshold Σ_{j≤k}2^j/j; the
    exact-start 3A run of A2E (A2E_MORE) is held on the symmetric state; the md is `report` of the raw, byte for byte."""
    import exp_cohesion_leak as E
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    far = MWE.fit_split([[-0.08, 0.06], [0.0, 0.5]], [1.0, 1.0], cs, ss, As)          # a broad kernel 7 off feature 2: not its main
    assert far['main_err'] == [0.0, None] and far['n_rest'] == 1 and MWE.SPLIT_CS == 1.5
    assert abs(far['rest_mass'] - np.sqrt(2 * np.pi) * 0.5) < 1e-15
    ok = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06]], [1.0, 3.4], cs, ss, As)           # inflated by 2.4: still the main kernel
    assert ok['main_err'][0] == 0.0 and ok['main_err'][1] is not None and abs(ok['main_err'][1] - 2.4) < 1e-12
    sp = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [0.0, 1.6], [0.5, 0.05]], [1.0, 1.0, 0.45, -0.6], cs, ss, As)
    assert abs(sp['rest_mass'] - np.sqrt(2 * np.pi) * (0.45 * 1.6 - 0.6 * 0.05)) < 1e-14 and sp['rest_sw'] < 0 < sp['rest_mass']
    oth = MWE.MW.BumpsCW(MWE.XS, MWE.S0).synth(np.array([[0.0, 1.6], [0.5, 0.05]]), np.array([0.45, -0.6]))
    assert abs(sp['rest_rms'] - np.sqrt(np.mean(oth ** 2)) / np.max(np.abs(MWE.feature_sum(cs, ss, As)))) < 1e-15
    assert ', mass +' in MWE.split_cells(sp)
    assert [round(MWE.bdf_bound(k), 3) for k in range(1, 6)] == [2.0, 4.0, 6.667, 10.667, 17.067]
    assert MWE.A2E_MORE == ((1e13, 0.25, 3.0, 0.0, 10.0),)            # round 100: its steps to s = 10
    raw = E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    x = next(x for x in raw['A2E'] if x['wmul'] == 3.0)
    assert max(abs(y['dc']) for y in x['rows']) < 1e-9 and x['off'] == 0.0
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    keys = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'H3', 'HL', 'HP', 'GL', 'GT', 'A2', 'HLL', 'SQ', 'A2E', 'A2P', 'A2H', 'A2R',
            'A2T')
    assert MWE.KEYS == keys and MWE.report(*[raw.get(k) for k in keys]) == md   # rounds 130, 133, 137: `--from-raw` reports every part
    assert "From 3A at σ*/σ₀ = 0.25 with the centre at c* exactly" in md and "its BDF order-2 steps at hλ 2.41–2.41" in md


def test_round97_the_cap_the_window_mass_the_null_test_and_the_orders():
    """Referee round 97 (D3cr-2 … D3cr-5, D3cr-9): the cap acts on the (c, σ) error alone — a main kernel inflated in w by 3 (or 10)
    with exact (c, σ) stays the main kernel, one 1.6 off in c does not —; main_cs, near_cs (unused kernels only); the others' mass on the
    window Σ oth Δx beside the continuum's, a broad kernel's tail beyond the window counted only by the latter, a negative-width kernel's
    mass by |σ| (FM2); `split_summary`'s gap and null clusters with NULL_TOL; the raw's span-10 BDF orders on the symmetric state;
    `a2_invariance`; the roundoff filter A2E_ROUND."""
    import exp_cohesion_leak as E
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    for wm in (4.0, 11.0):
        sp = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06]], [1.0, wm], cs, ss, As)
        assert sp['main_err'] == [0.0, wm - 1] and sp['main_cs'] == [0.0, 0.0] and sp['main_dw'] == [0.0, wm - 1]
    off = MWE.fit_split([[-0.08, 0.06], [0.08 + 1.6 * 0.06, 0.06]], [1.0, 1.0], cs, ss, As)
    assert off['main_err'][1] is None and abs(off['near_cs'][1] - 1.6) < 1e-12 and off['near_cs'][0] is None
    edge = MWE.fit_split([[-0.08, 0.06], [0.08 + 1.4 * 0.06, 0.06]], [1.0, 1.0], cs, ss, As)
    assert abs(edge['main_cs'][1] - 1.4) < 1e-12
    used = MWE.fit_split([[-0.08, 0.06]], [1.0], cs, ss, As)          # feature 2 nearest to feature 1's kernel, used: no unused kernel
    assert used['main_err'][1] is None and used['near_cs'] == [None, None]
    broad = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [3.0, 0.6]], [1.0, 1.0, 0.5], cs, ss, As)   # outside the window
    oth = MWE.MW.BumpsCW(MWE.XS, MWE.S0).synth(np.array([[3.0, 0.6]]), np.array([0.5]))
    assert broad['rest_win'] == float(np.sum(oth) * (MWE.XS[1] - MWE.XS[0])) and 0 < broad['rest_win'] < 0.02 * broad['rest_mass']
    neg = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [0.5, -0.2]], [1.0, 1.0, 0.3], cs, ss, As)
    assert abs(neg['rest_mass'] - np.sqrt(2 * np.pi) * 0.3 * 0.2) < 1e-15 and neg['rest_mass'] > 0   # FM2: |σ|
    assert MWE.NULL_TOL == 1e-8
    nul = dict(name='n', moving=dict(theta=[[-0.08, 0.06], [0.08, 0.06], [0.5, 0.06], [0.5, 0.06]], w=[1.0, 1.0, 0.3, -0.3], rms=1e-9),
               cs=cs, ss=ss, As=As)
    big = dict(name='b', moving=dict(theta=[[-0.08, 0.06], [0.08, 0.06], [0.5, 0.06], [0.5, 0.07]], w=[1.0, 1.0, 0.3, -0.3], rms=1e-9),
               cs=cs, ss=ss, As=As)
    mid = dict(name='m', moving=dict(theta=[[-0.08, 0.06], [0.08, 0.06], [0.5, 0.06], [0.5, 0.06 + 1e-8]], w=[1.0, 1.0, 0.3, -0.3], rms=1e-9),
               cs=cs, ss=ss, As=As)                                     # K13: the others' peak 3.5e-8, their rms 7.7e-9 — by the rms
    sm = MWE.end_split(mid)
    assert sm['rest_foot'] > MWE.NULL_TOL > sm['rest_rms']
    lines = MWE.split_summary(H=[nul, big, mid])
    assert "over the 3 splits above" in lines[1]
    assert "of max f*; referee round 100): H n, H m." in lines[1]
    gap = MWE.split_summary(H=[dict(name='g', moving=dict(theta=[[-0.08, 0.06], [0.08 + 2.9 * 0.06, 0.06], [-0.08 + 0.5 * 0.06, 0.06]],
                                                          w=[1.0, 1.0, 0.02], rms=1e-9), cs=cs, ss=ss, As=As)])[1]   # the 0.02 kernel 2.17 off feature 2
    assert "1 features unassigned, 1 of them with an unused kernel of their sign, at least 2.17 off in (c, σ) (H g)" in gap
    assert "the main kernels' (c, σ) errors at most 0 (H g)" in gap and MWE.split_summary() == []
    raw = E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    x = next(x for x in raw['A2E'] if x['span'] == 10.0)
    st = [y for y in x['steps'] if abs(y[3]) < 1e-6 and y[2] == 2]
    assert [y[4] for y in st] == [2, 2, 2, 1, 1, 1, 1, 1, 1] and abs(st[0][0] - 0.00264) < 1e-5     # D3cr-5: NQU, not NQCUR or MUSED
    assert MWE.A2E_ROUND == 1e-14
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    assert "its BDF order-1 steps at hλ 3.99–7.56 multiply the offset by -0.335 … -0.152" in md and "1 further BDF step left out, where the offset falls below 1e-14σ*" in md
    n, fu, cr, cy, lg = MWE.a2_invariance(0.25)
    assert n == 6 and fu == 0.0 and cr < 1e-11 and cy < 1e-14 and abs(lg + (1.55 + 0.03) ** 2 / (2 * 0.03 ** 2) / np.log(10)) < 1e-9
    n, fu, cr, cy, lg = MWE.a2_invariance(4.0)
    assert abs(fu - np.exp(-(1.55 + 0.03) ** 2 / (2 * 0.48 ** 2))) < 1e-9 * fu and cr > 1e-5 and cy < 1e-14   # the unpaired sample nearest c*: −1.55
    assert "The problem as posed is a perturbation of the symmetric one" in md and "held on the symmetric state by LSODA's damping" in md



def test_round87_sq_args_sq_job_a2e_order_and_split(monkeypatch):
    """Referee round 87 (D3ch-5: SA3, SJ4, AE12, FS12, SC3): `sq_args` passes G's centres, widths and amplitudes field by field; `sq_job`
    takes the formed threshold from the smallest amplitude (two amplitudes); `part_A2E` runs SPANS_A2E in order, then A2E_MORE;
    `fit_split`'s main kernels' w/A − 1 is over each feature's own amplitude, its error and w/A − 1 apart."""
    G = [dict(q=4, status=0, cs=[-0.08, 0.08], ss=[0.06, 0.05], As=[1.0, 0.4])]
    a = MWE.sq_args(G)[0]
    assert a[0] == 4 and list(a[1]) == [-0.08, 0.08] and list(a[2]) == [0.06, 0.05] and list(a[3]) == [1.0, 0.4]
    monkeypatch.setattr(MWE, 'TAUS_G', (10.0,))
    monkeypatch.setattr(MWE, 'TAUS_FINE', (10.0,))
    monkeypatch.setattr(MWE, 'EPS_G', (1e-6,))
    monkeypatch.setattr(MWE, 'EPS_GX', ())
    Wt = [np.zeros(3), np.array([1.0, 0.004, 0.006])]                   # 1% of the smallest amplitude 0.4 is 0.004
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a_, **k: (None, Wt, None, None, None))
    r = MWE.sq_job((4, [-0.08, 0.08], [0.06, 0.05], [1.0, 0.4]))
    assert r['seq'][1e-6] == [(10.0, [0, 2], [])]                       # 0.006 > 0.004 formed; by the largest (0.01) it would not be
    got = []
    monkeypatch.setattr(MWE, 'a2e_job', lambda *a_: got.append(a_) or a_)
    assert MWE.part_A2E() == [(s_,) for s_ in MWE.SPANS_A2E] + list(MWE.A2E_MORE) and MWE.SPANS_A2E == (10.0, 1e13)
    sp = MWE.fit_split([[-0.08, 0.06], [0.08, 0.05]], [1.1, 0.3], (-0.08, 0.08), (0.06, 0.05), (1.0, 0.4))
    assert abs(sp['main_dw'][0] - 0.1) < 1e-12 and abs(sp['main_dw'][1] + 0.25) < 1e-12 and abs(sp['main_err'][1] - 0.25) < 1e-12
    sp2 = MWE.fit_split([[-0.07, 0.06], [0.08, 0.05]], [1.0, 0.4], (-0.08, 0.08), (0.06, 0.05), (1.0, 0.4))
    assert abs(sp2['main_err'][0] - 0.01 / 0.06) < 1e-12 and sp2['main_dw'][0] == 0.0          # SC3: the error not w/A − 1


def test_round100_the_lemma_the_fixed_point_and_one_null_predicate():
    """Referee round 100 (D3cu-1 … D3cu-4): Lemma D‴.59a's split sums to the exact Jacobian's ∂ċ/∂c, its centre term the overshoot
    (Gm w/(Nσ²))(w − f*(c*)); at the held runs the implicit step's fixed point −d/λ is the held offset to 5%; `null_cluster` is
    `outcome`'s test and `split_summary`'s — the referee's planted ends (a main kernel deficient by 5·10⁻⁷ made up by a pair at its θ: not
    null, the others' rms 7·10⁻⁸; main kernels off by 5·10⁻⁷ with a far cancelling pair: not null, the residual 10⁻⁷) and R19 (a main
    kernel 5·10⁻⁴ off, the others nothing: not null)."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    for r in raw['A2E']:
        sym = [x for x in r['rows'] if abs(x['dc']) < 1e-12 and x['jcc'] > 0]
        if sym:
            x = max(sym, key=lambda x: x['jcc'])
            parts, tot, J, hs = MWE.a2_lemma(r['ratio'], x['w'], x['sig'])
            assert abs(tot - J) < 1e-9 * abs(J) and abs(sum(parts.values()) - tot) < 1e-9 * abs(J) and parts['0'] > 0 > parts['1']
            L, K, ss, z = MWE.a2_state(r['ratio'], x['w'], x['sig'])
            wv, s_ = float(MWE.RT.critical_w(L, z[:1])[0]), float(z[2])
            assert abs(parts['0'] - L.G * L.m * wv / MWE.XS.size / s_ ** 2 * (wv - MWE.A_A)) < 1e-9 * abs(parts['0'])
        last = r['rows'][-1]
        if abs(last['dc']) < 1e-12:
            d, lam, fp, held = MWE.a2_fixed_point(r['ratio'], last['w'], last['sig'], last['dc'])
            assert lam > 100 and d < 0 < fp and abs(fp / held - 1) < 0.05
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    K = MWE.MW.BumpsCW(MWE.XS, MWE.S0)
    f = MWE.feature_sum(cs, ss, As)
    res = lambda th, w: float(np.sqrt(np.mean((f - K.synth(np.array(th), np.array(w))) ** 2)))
    th1, w1 = [[-0.08, 0.06], [0.08, 0.06], [0.08, 0.06], [0.08, 0.06]], [1.0, 1.0 - 5e-7, 0.3, -0.3 + 5e-7]
    assert res(th1, w1) < 1e-8 and not MWE.null_cluster(th1, w1, cs, ss, As, res(th1, w1))
    th2, w2 = [[-0.08 + 5e-7, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05 + 1e-9]], [1.0, 1.0, 0.3, -0.3]
    assert not MWE.null_cluster(th2, w2, cs, ss, As, 1.1e-7)
    th3, w3 = [[-0.08 + 5e-4 * 0.06, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05]], [1.0, 1.0, 0.3, -0.3]
    assert not MWE.null_cluster(th3, w3, cs, ss, As, 1e-9) and MWE.fit_split(th3, w3, cs, ss, As)['rest_rms'] < 1e-8
    th4 = [[-0.08, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05]]
    assert MWE.null_cluster(th4, w3, cs, ss, As, 1e-9)
    m = dict(rms=1e-9, err=[0.0, 0.0, 0.0], theta=th3, w=w3, oob=False)
    assert MWE.outcome(m, As, cs, ss) != 'exact, null cluster' and MWE.outcome(dict(m, theta=th4), As, cs, ss) == 'exact, null cluster'
    r = dict(name='p', moving=dict(theta=th3, w=w3, rms=1e-9), cs=cs, ss=ss, As=As)
    assert MWE.split_summary(H=[r])[1].split(" On the targets symmetric")[0].endswith("referee round 100): none.")


def test_round100_planted_edges():
    """Round 100's own mutants V10, V12, V13: `a2e_job`'s steps_to records LSODA's steps past S_A2E_STEPS[1] (to 0.02 here); the null
    predicate's residual threshold is 10⁻⁸ (a residual of 5·10⁻⁸ is not null) and its w test absolute (a main kernel 5·10⁻⁶ off in w is
    not null)."""
    r = MWE.a2e_job(10.0, ratio=0.25, wmul=3.0, off=0.0, steps_to=0.02)
    assert r['steps_to'] == 0.02 and r['steps'][-1][0] >= 0.02 > MWE.S_A2E_STEPS[1] and r['steps'][-2][0] < 0.02
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    th = [[-0.08, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05]]
    w = [1.0, 1.0, 0.3, -0.3]
    assert MWE.null_cluster(th, w, cs, ss, As, 9e-9) and not MWE.null_cluster(th, w, cs, ss, As, 5e-8)
    assert not MWE.null_cluster(th, [1.0, 1.0 + 5e-6, 0.3, -0.3], cs, ss, As, 1e-9)


def test_round103_the_held_state_the_exit_the_lemma_s_terms_and_mirror_symmetry():
    """Referee round 103 (D3cx-1 … D3cx-7): in double the centre sample lies 8 ulps of c* above c*, and at the held runs the zero of ċ
    along the line lies at 8.1 and 9.9 ulps, the held offsets 8 and 9 ulps; the explicit flow from c* leaves at s ≈ 0.0035 (5A) and 0.048
    (3A) on A's samples and alike on exactly symmetric ones; the lemma's ±1 term turns positive below σ/h ≈ 0.45 at w = 4A and the total
    is stable (negative) at σ/h = 1; `mirror_asym` is 0 on a mirror-symmetric configuration and the polished ends' value is printed;
    the planted ends W26 (an unformed kernel in `outcome`'s null path), W29 (a residual of 1.1·10⁻⁷ through `split_summary`), W30, W31 (an
    HL and an HLL null row) and W36 (a cap change inside (1.5, 2.5))."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    held = [(r, r['rows'][-1]) for r in raw['A2E'] if abs(r['rows'][-1]['dc']) < 1e-12]
    for r, x in held:
        k0, kz, u = MWE.a2_ulps(r['ratio'], x['w'], x['sig'])
        assert k0 == 8.0 and abs(u - np.spacing(MWE.C_A) / (r['ratio'] * MWE.S0)) < 1e-30 and abs(x['dc'] / u - round(x['dc'] / u)) < 1e-6
        assert abs(kz - x['dc'] / u) < 1.0
    assert 0.003 < MWE.a2_exit(0.25, 5.0) < 0.0036 < 0.0039 < MWE.a2_exit(0.25, 5.0, sym=True) < 0.0045
    assert 0.045 < MWE.a2_exit(0.25, 3.0) < 0.05 < MWE.a2_exit(0.25, 3.0, sym=True) < 0.055 and MWE.a2_exit(0.25, 3.0, smax=0.01) is None
    sc = MWE.a2_lemma_scan()
    w4 = {q: (a, b, t) for q, a, b, t, _ in sc['w']}
    assert w4[1.0][2] < 0 < w4[0.6][2] and w4[0.6][1] < 0 < w4[0.45][1] and all(a > 0 for a, _, _ in w4.values())
    ls = {q: (a, b, t, wa) for q, a, b, t, wa in sc['ls']}
    assert ls[0.3][1] > ls[0.3][0] > 0 and abs(ls[0.2][3] - 1) < 1e-4
    assert MWE.mirror_asym([[-0.1, 0.06], [0.1, 0.06]], [1.0, 1.0]) == 0.0 and abs(MWE.mirror_asym([[-0.1, 0.06], [0.12, 0.06]], [1.0, 1.0]) - 0.02) < 1e-15
    assert abs(MWE.mirror_asym([[-0.1, 0.06], [0.1, 0.08]], [1.0, 1.0]) - 0.02) < 1e-15                    # σ counted too
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    assert "the location the problem's, the stability LSODA's step sequence's" in md and "On the targets symmetric about 0" in md
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    m = dict(rms=1e-10, err=[1e-8] * 3, oob=False, theta=[[-0.08, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05], [0.3, 0.06]], w=[1.0, 1.0, 0.3, -0.3, 0.005])
    assert MWE.outcome(m, As, cs, ss) == 'exact, null cluster'                       # W26: the unformed 0.005 kernel left out
    th = [[-0.08, 0.06], [0.08, 0.06], [0.9, 0.05], [0.9, 0.05]]
    r29 = dict(name='q', moving=dict(theta=th, w=[1.0, 1.0, 0.3, -0.3], rms=1.1e-7), cs=cs, ss=ss, As=As)
    assert MWE.split_summary(H=[r29])[1].split(" On the targets symmetric")[0].endswith("referee round 100): none.")
    loop = dict(tau=1.0, eps=1e-8, theta=th, w=[1.0, 1.0, 0.3, -0.3], rms=1e-9)
    hl = dict(name='σ* 0.06, equal', cs=list(cs), ss=list(ss), As=list(As), loops=[loop])
    assert "HL σ* 0.06, equal 1e-08 τ 1" in MWE.split_summary(HL=[hl])[1].split('max f*; referee round 100): ')[1]
    _, cs6, ss6, As6 = next(x for x in MWE.H_TARGETS if x[0] == 'σ* 0.06, equal')
    th6 = [[cs6[0], ss6[0]], [cs6[1], ss6[1]], [0.9, 0.05], [0.9, 0.05]]
    hll = dict(name='σ* 0.06, equal', eps=1e-8, loops=[dict(loop, theta=th6, w=[As6[0], As6[1], 0.3, -0.3])])
    assert "HLL σ* 0.06, equal 1e-08 τ 1" in MWE.split_summary(HLL=[hll])[1].split('max f*; referee round 100): ')[1]
    assert f"the loops' formed kernels (HL, HLL) at most {MWE.mirror_asym(th6, [As6[0], As6[1], 0.3, -0.3]):.1e}" in MWE.split_summary(HLL=[hll])[1]   # HLL's loops counted
    th36 = [[-0.08, 0.06], [0.08 + 2.0 * 0.06, 0.06]]                    # W36: feature 2's kernel 2.0 off: unassigned at 1.5, main at 2.5
    lo, hi = MWE.cap_interval([(th36, [1.0, 1.0], cs, ss, As)])
    assert hi == 2.0 + 0.0 or abs(hi - 2.0) < 1e-12


def test_round106_the_exit_s_conclusion_needs_the_offset_start(monkeypatch):
    """Round 106's Z04: the 'the exit is the data's' clause needs the symmetric problem started at the offset to leave downward too."""
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    mk = lambda: {w: dict(a=[(0.05, -1), (0.05, -1)], imp=[('Radau', (0.06, -1)), ('BDF', (None, 0))], sym=[(0.05, 1), (0.05, 1)], kz=9.0,
                          symoff=(0.05, -1), summ=dict(fsum=((0.05, -1), (None, 0)), rev=((0.05, -1), (0.05, -1))),
                          ld=dict(a=[(0.05, -1), (0.05, -1)], sym=[(0.06, 1), (0.06, -1)], seed=-7e-17, symS=(0.05, -1), ulp=3.47e-18))
                for w in (3.0, 5.0)}
    fake = mk()
    monkeypatch.setattr(MWE, 'a2_exits', lambda *a, **k: fake)
    got = lambda: '\n'.join(MWE.a2_lines(raw['A2'], raw['A2E']))
    assert "the exit is the data's" in got() and "its seed is the arithmetic's" in got()
    fake[3.0]['symoff'] = (0.05, 1)                                     # Z04: the offset start's direction
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[5.0]['symoff'] = None                       # round 109 (R06): a missing offset start fails
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['a'] = [(0.05, 1), (0.05, 1)]          # R09: A's own exits upward
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[5.0]['summ']['fsum'] = ((0.05, -1), (0.07, 1))   # the symmetric problem leaving with exact sums
    assert "its seed is the arithmetic's" not in got() and "the exit is the data's" in got()
    fake.update(mk()); fake[3.0]['summ']['rev'] = ((0.06, -1), (0.05, -1))   # A's time moving by 20% with the order
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['summ']['rev'] = ((0.05, -1), (0.05, 1))    # S06: the reversed symmetric run not flipping
    assert "its seed is the arithmetic's" not in got()
    assert "the integrator's tolerance flips the direction (3A in long double, 5A in long double)" in (fake.update(mk()) or got())
    import math
    late = 0.06 / 0.05 - 1
    assert f"in long double the exit comes {math.floor(late * 100)}–{math.ceil(late * 100)}% later" in got()
    fake.update(mk()); fake[3.0]['summ']['fsum'] = ((0.05, 1), (None, 0))   # round 112 (T03): A's exact-sum exit upward
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[5.0]['summ']['rev'] = ((0.05, 1), (0.05, -1))    # A's reversed exit upward
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['ld']['symS'] = (0.0505, -1)                # round 112: the seeded symmetric run 1% off A's time
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['ld']['symS'] = (0.05, 1)                   # … or the other way
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[5.0]['ld']['a'] = [(0.05, -1), (0.0501, -1)]     # the long-double exit moving with the tolerance
    fake[5.0]['ld']['symS'] = (0.0501, -1)                                    # (the seeded run at the second tolerance's time: V08)
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['sym'] = [(0.05, 1), (0.05, -1)]             # V07: double flipping, long double not
    fake[3.0]['ld']['sym'] = [(0.06, 1), (0.06, 1)]
    assert "the integrator's tolerance flips the direction (3A, 5A in long double)" in got()
    fake.update(mk()); fake[5.0]['ld']['a'] = [(0.0515, -1), (0.0515, -1)]   # A's double runs 3% off the long-double time
    fake[5.0]['ld']['symS'] = (0.0515, -1)
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['ld']['sym'] = [(0.051, 1), (0.051, 1)]     # long double only 2% later: no arithmetic's clause
    assert "its seed is the arithmetic's" not in got()
    fake.update(mk()); fake[3.0]['symoff'] = (0.06, -1)                       # round 112: the offset start's time is no condition
    assert "the exit is the data's" in got()
    fake.update(mk()); fake[3.0]['ld']['a'] = [(0.05, 1), (0.05, 1)]          # round 116 (D3dk-7, X10): A's long-double exits upward
    fake[3.0]['ld']['symS'] = (0.05, 1)                                       # (the seeded run with them: only the direction fails)
    assert "the exit is the data's" not in got()
    fake.update(mk()); fake[3.0]['ld']['symS'] = (0.0500025, -1)              # round 116 (D3dk-2): past the integrator's resolution
    assert "A's time to 5.0e-05" in got() and "the exit is the data's" in got()
    fake.update(mk()); fake[3.0]['a'] = [(0.0505, -1), (0.05, -1)]            # round 116 (D3dk-4): the largest double deviation printed
    assert "in double its time within 1.0% of the long-double one" in got()
    fake[3.0]['a'] = [(0.05051, -1), (0.05, -1)]                              # round 116's W08: 1.02% printed as 1.1%, rounded up
    assert "in double its time within 1.1% of the long-double one" in got()
    fake.update(mk())
    for w in fake:                                                            # round 116 (D3dk-2): the arithmetic's flip and Λ
        fake[w]['ld']['sym'], fake[w]['ld']['lag'] = [(0.06, -1), (0.06, -1)], [7.0, 8.0]
    t = got()
    assert ("the arithmetic flips it at both tolerances (3A: +, + in double and −, − in long double at rtol 1e-10, 1e-12; 5A: +, + in "
            "double and −, − in long double at rtol 1e-10, 1e-12)" in t and "the integrator's tolerance flips" not in t
            and "by Λ = ∫λ ds 7.00 … 8.00 along the route" in t and "the two seeds' ratio e^Λ) —, each in the code's order one draw" in t)
    fake[5.0]['ld']['lag'] = [7.0, None]                                      # a run that does not leave: no Λ
    fake[3.0]['ld']['sym'] = [(0.06, 1), (0.06, -1)]                          # one tolerance flipping with the arithmetic: no such flip
    t = got()
    assert ("by Λ" not in t and "the arithmetic flips it at both tolerances (5A: +, + in double and −, − in long double at rtol 1e-10, "
            "1e-12)" in t)
    fake.update(mk())                                                         # round 125 (D3dt-6): opposite flips at the two tolerances,
    fake[3.0]['sym'], fake[3.0]['ld']['sym'] = [(0.05, 1), (0.05, -1)], [(0.06, -1), (0.06, 1)]   # each tolerance's directions printed
    assert "(3A: +, − in double and −, + in long double at rtol 1e-10, 1e-12)" in got()
    fake.update(mk())                                                         # round 125 (D3dt-7): S in units of σ* (not σ₀) and the
    for w in fake:                                                            # linear response against the exit's level 10⁻²σ*
        fake[w]['ld'].update(ss=0.25 * MWE.S0, Lam=math.log(0.01 * 0.25 * MWE.S0 / 7e-17))
    t = got()
    assert t.count("S = -7.000e-17 (-2.33e-15σ*)") == 2 and t.count("e^{Λ(s_e)}|S| 1.000 times the exit's level") == 2
    assert "on A's samples 0.050000 (−) at rtol 1e-10 and 1e-12 (apart by 0.0e+00 relatively)" in t and "(scipy, rtol 1e-10 unless stated)" in t
    fake[5.0]['ld']['a'] = [(0.05, -1), (0.05000001, 1)]                      # directions differing: both printed, to five digits
    assert "5A: on A's samples 0.050000 (−), 0.050000 (+) (rtol 1e-10, 1e-12)" in got()


def test_round106_the_exit_s_seed_the_exact_zero_the_bottleneck_and_the_lemma_s_grid():
    """Referee round 106 (D3da-1 … D3da-7): on A's samples the exit from c* is downward for DOP853 at two tolerances and for scipy's
    Radau, and the exactly symmetric samples started at the posed problem's offset leave the same way at nearly A's time, while from
    their centre they leave upward (5A); the zero of ċ at the held states in 50 digits (8.12, 9.57 ulps) within an ulp of the held
    offsets; `mirror_asym` a bottleneck matching ({a, a, Ra} 0.2), `mirror_odd` gauge-invariant; the lemma's sign edges on the grid; the
    md's stale clauses gone and the exact-start clause from the data."""
    ex = MWE.a2_exits(0.25, (5.0,), {5.0: 8.118})[5.0]
    assert all(d < 0 for _, d in ex['a']) and ex['imp'][0][1][1] < 0 and ex['symoff'][1] < 0 and all(d > 0 for _, d in ex['sym'])
    assert abs(ex['symoff'][0] / ex['a'][0][0] - 1) < 0.02 and 0.0034 < ex['a'][0][0] < 0.0036
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    held = [(r, r['rows'][-1]) for r in raw['A2E'] if abs(r['rows'][-1]['dc']) < 1e-12]
    z = {r['wmul']: (MWE.a2_ulps_exact(r['ratio'], x['w'], x['sig']), x['dc'] / MWE.a2_ulps(r['ratio'], x['w'], x['sig'])[2]) for r, x in held}
    assert abs(z[5.0][0] - 8.1183) < 1e-3 and abs(z[3.0][0] - 9.5705) < 1e-3 and all(abs(a - b) < 0.75 for a, b in z.values())
    e3 = MWE.a2_exits(0.25, (3.0,), {3.0: z[3.0][0]})[3.0]                  # round 109: 3A too, and the arithmetic varied
    assert 0.047 < e3['a'][0][0] < 0.050 and e3['a'][0][1] < 0 and e3['summ']['fsum'][1][0] is None
    assert e3['summ']['rev'][1][1] == -e3['sym'][0][1] and abs(e3['summ']['fsum'][0][0] / e3['a'][0][0] - 1) < 0.02
    assert MWE.mirror_asym([[0.1, 0.06], [0.1, 0.06], [-0.1, 0.06]], [1.0, 1.0, 1.0]) == 0.2      # the nearest image gave 0
    assert MWE.mirror_asym([[0.1, 0.06], [-0.1, 0.06], [0.0, 0.05]], [1.0, 1.0, 0.5]) == 0.0
    th = [[-0.08, 0.06], [-0.08, 0.06], [0.08, 0.06]]                       # a coinciding pair's gauge: asymmetric parameters
    assert MWE.mirror_asym(th, [0.4, 0.6, 1.0]) > 0.1 and MWE.mirror_odd(th, [0.4, 0.6, 1.0]) < 1e-15
    sg = MWE.a2_lemma_signs()
    assert sg['w'][2] is None and all(abs(a_ - b_) < 1e-3 for a_, b_ in zip(sg['w'][:2] + sg['ls'], (0.8437, 0.4696, 0.7746, 0.6195, 0.5583)))
    for lab, q in (('w', sg['w'][0]), ('ls', sg['ls'][0])):              # round 109: each edge a root of its term
        lo, hi = (next(iter(MWE.a2_lemma_scan(0.25, (x,))[lab]))[3] for x in (q - 1e-6, q + 1e-6))
        assert lo > 0 > hi
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    assert "a null cluster where it fits none" not in md and "seeded by roundoff alike" not in md and "neither damps nor holds" not in md
    assert "on A's samples every exit downward" in md and "a bottleneck matching" in md and "symmetric about 0 to 4.4e-16" in md
    r = dict(raw['A2E'][-1])
    rows = [dict(x) for x in r['rows']]
    rows[-1]['dc'] = 0.1                                                    # a run that leaves: the clause follows the data
    got = '\n'.join(MWE.a2_lines(raw['A2'], list(raw['A2E'][:-1]) + [dict(r, rows=rows)]))
    assert "— leaving the symmetric state" in got
    rows[-1]['dc'], rows[len(rows) // 2]['dc'] = 1e-15, 1e-6                 # S10: held at the end after a transient: still held
    got = '\n'.join(MWE.a2_lines(raw['A2'], list(raw['A2E'][:-1]) + [dict(r, rows=rows)]))
    assert "— leaving the symmetric state" not in got


def test_round112_the_long_double_exit_the_route_seed_and_the_width_split():
    """Referee rounds 112 (D3dg-1, -2) and 116 (D3dk-1, -2, -5): in long double with the offset c − c* as the state and the stored
    double f* as the data, A's exit (3A) is tolerance-free to 10⁻⁴ and downward at 0.04830 (0.04827 with f* re-evaluated in long
    double); the data's drift carried along the route (`a2_seed`) is 2.32–2.35 times c*'s distance below the zero of ċ (2.37 with the
    re-evaluated f*), its linear response e^{Λ(s_e)}|S| the exit's level to 1%, and, planted in the exactly symmetric problem, gives
    A's exit to 3·10⁻⁵; the symmetric samples' long-double exit comes > 10% later than in double, by a Λ that is one draw (A2P's code
    order; round 125, D3dt-1); the mirror
    evidence by width (σ ≥ h) is odd at 1.2·10⁻⁹ while the formed set at a threshold of 0.95% of A is not."""
    import math
    import pytest
    if np.finfo(np.longdouble).eps > 1e-18:
        pytest.skip('no extended long double on this platform')
    ex = [MWE.a2_exit_ld(0.25, 3.0, rtol=q) for q in (1e-10, 1e-12)]
    assert all(0.04829 < t < 0.04830 and d == -1 for t, d in ex) and abs(ex[0][0] / ex[1][0] - 1) < 1e-4
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    e5 = [MWE.a2_exit_ld(0.25, 5.0, rtol=q) for q in (1e-10, 1e-12)]         # round 130 (D3dy-5, -6): "apart by" relative, from the data
    for a, b in (ex, e5):
        assert f"at rtol 1e-10 and 1e-12 (apart by {abs(a[0] / b[0] - 1):.1e} relatively)" in md and abs(a[0] - b[0]) < 1e-6
    S, e, sS, Le = MWE.a2_seed(0.25, 3.0)
    assert e == ex[1] and S < 0 and 2.32 < S / (-9.5705 * float(np.spacing(MWE.C_A))) < 2.35
    assert sS[1] == -1 and abs(sS[0] / e[0] - 1) < 3e-5 and abs(math.exp(Le) * abs(S) / (0.01 * 0.25 * MWE.S0) - 1) < 0.01
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))         # round 125 (D3dt-1): one draw — against A2P's code order,
    p0 = [x['rows'][0] for x in raw['A2P'] if (x['wmul'], x['rtol']) == (3.0, 1e-12)]   # Λ along the line (not 11 ln 2: a draw)
    lag = MWE.a2_sym_lag(0.25, 3.0, 1e-12)
    assert len(p0) == 2 and all(q['d'] == p0[0]['d'] and q['ld'] == p0[0]['ld'] for q in p0)   # one code order, a record per seed
    assert 5 < lag < 13 and all(abs(lag - q['lag']) < 0.02 for q in p0)
    assert (f"(its Λ within {max(abs(lag - q['lag']) for q in p0):.1e} of `a2_sym_lag`'s, along"
            in MWE.a2p_line(raw['A2P'], {(3.0, 1e-12): lag}))
    sy, sd = MWE.a2_exit_ld(0.25, 3.0, sym=True, rtol=1e-10), MWE.a2_exit(0.25, 3.0, sym=True, signed=True)
    assert sy[0] > 1.1 * sd[0]
    r = next(x for x in raw['H'] if x['name'] == 'σ* 0.06, equal')
    th, w = np.atleast_2d(np.asarray(r['moving']['theta'], float)).reshape(-1, 2), np.asarray(r['moving']['w'], float)
    Am, h = max(abs(a) for a in r['As']), MWE.XS[1] - MWE.XS[0]
    wide = (w != 0) & (np.abs(th[:, 1]) >= h)
    assert wide.sum() == 11 and MWE.mirror_odd(th[wide], w[wide], A=Am) < 2e-9
    f95 = np.abs(w) > 0.0095 * min(abs(a) for a in r['As'])
    assert MWE.mirror_odd(th[f95], w[f95], A=Am) > 1e-3                       # the formed set's oddness rests on the threshold
    assert ("the 11 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.041, the widest narrower one 0.001) odd by "
            "1.2e-09 on and between the samples (9 formed, the narrowest of them 0.0496; 2 idle); the formed threshold's margin: the largest "
            "unformed |w| 0.0098 (over max A, against the threshold 0.01 of min A), the smallest formed 0.06") in md
    assert "the formed kernels' distance from their mirror image (`mirror_asym`" in md and "fixed-point subspace (`mirror_asym`" not in md


def test_round125_the_seed_s_distribution_the_no_exit_branch_and_the_md_s_units(monkeypatch):
    """Referee round 125 (D3dt-1, -2, -5, -7) as revised in round 130 (D3dy-1, -3, -4, -6): A2P — the symmetric samples' exit over the
    code's order and K_A2P random orders per fresh seed, in double and in long double — is the code's configuration and reproduces from
    the code (the orders, the exits, ln|S| along the line, the job's first orders at 3A, rtol 1e-12); each cell's exponent Λ/(11 ln 2)
    (the random orders, both seeds) has its 99% interval inside EXP_A2P (rounds 133, 137: D3eb-1, D3ef-2) with each arithmetic's seeds O(1) in
    units of uσ* (not one draw), the median lateness 10–35%; the md's numbers are the data's and its conclusion follows those two
    conditions; `a2_seed` without an exit returns
    four values and the exit line prints without the seed; the BDF entry's tolerance, the mirror block's idle count and narrowest formed
    kernel, and the line's tangential rates and verdicts at A2E's states."""
    import math
    import pytest
    if np.finfo(np.longdouble).eps > 1e-18:
        pytest.skip('no extended long double on this platform')
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    P = raw['A2P']
    # round 130 (D3dy-6): the raw's configuration is the code's — a record per w, rtol and seed (fresh seeds: not round 125's 125 nor
    # round 130's referee's 130 and 2026), each the code's order and K_A2P random orders
    assert [(r['wmul'], r['rtol'], r['k'], r['seed'], r['max_step']) for r in P] == [
        (w, q, MWE.K_A2P, s_, None) for w in MWE.WMUL_A2P for q in MWE.RTOL_A2P for s_ in MWE.SEED_A2P]
    assert MWE.K_A2P == 100 and not {125, 130, 2026} & set(MWE.SEED_A2P) and all(len(r['rows']) == r['k'] + 1 for r in P)
    lnr = 11 * math.log(2)
    assert all(abs(math.log(r['u'][0] / r['u'][1]) - lnr) < 1e-12 for r in P)
    n = len(MWE.SEED_A2P) * MWE.K_A2P
    for c in [(w, q) for w in MWE.WMUL_A2P for q in MWE.RTOL_A2P]:     # each cell pooled over its seeds
        rs = [r for r in P if (r['wmul'], r['rtol']) == c]
        st = MWE.a2p_stats(rs)
        rows = [q for r in rs for q in r['rows'][1:]]                     # round 130: the code's order (row 0) is one draw, not in it
        x = [q['lag'] for q in rows]
        e = MWE.a2p_exp(st, lnr)                                          # round 133 (D3eb-1): the exponent's 99% interval
        assert st['n'] == n == len(x) and 0.9 <= e[1] < e[0] < e[2] <= 1.1 and all(0.05 < g < 20 for g in st['gm'])
        assert abs(e[2] - e[0] - 2.5758293035489 * st['sem'] / lnr) < 1e-9 and abs(e[0] - np.mean(x) / lnr) < 1e-12
        assert f"of {n}), the exponent {e[0]:.3f} (99% interval {e[1]:.3f} … {e[2]:.3f}), the seeds'" in MWE.a2p_line(P)
        assert abs(st['sem'] - np.std(x, ddof=1) / math.sqrt(n)) < 1e-12 and (st['lo'], st['hi']) == (min(x), max(x))
        assert abs(st['t'] - (np.mean(x) - lnr) / st['sem']) < 1e-9
        assert [st['q'][k_] for k_ in ('med', 'q1', 'q3')] == [np.median(x), np.percentile(x, 25), np.percentile(x, 75)]
        assert all(abs(q['lag'] - (q['lnSd'] - q['lnSld'])) < 1e-12 for q in rows)
        assert st['up'] == [(sum(q[k_][1] > 0 for q in rows), n) for k_ in ('d', 'ld')] and all(0.35 < u / m < 0.65 for u, m in st['up'])
        late = [q['ld'][0] / q['d'][0] - 1 for q in rows]                 # round 130 (D3dy-3): the median lateness, not the extremes
        assert st['late']['med'] == np.median(late) and (st['late']['lo'], st['late']['hi']) == (min(late), max(late))
        assert 0.1 < st['late']['med'] < 0.35 and st['late']['q1'] > 0.05
        for i, (k_, n_) in enumerate((('d', 'nd'), ('ld', 'nld'))):       # D3dy-2: DOP853's steps per arithmetic
            assert abs(st['steps'][i][0] - np.mean([q[n_] for q in rows])) < 1e-9
            assert abs(st['steps'][i][1] / np.mean([q[k_][0] / q[n_] for q in rows]) - 1) < 1e-12
        assert abs(st['steps'][1][1] / st['steps'][0][1] - 1) < 0.1       # the pair's mean steps within 10% of each other
    r = next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (3.0, 1e-12, MWE.SEED_A2P[0]))   # the code reproduces the stored draws
    rng = np.random.default_rng(MWE.SEED_A2P[0])
    p1 = rng.permutation(MWE.XS.size)                                        # the first random order
    assert r['rows'][0]['d'] == MWE.a2_exit(0.25, 3.0, sym=True, rtol=1e-12, signed=True)
    assert r['rows'][0]['ld'] == MWE.a2_exit_ld(0.25, 3.0, sym=True, rtol=1e-12)
    assert r['rows'][1]['d'] == MWE.a2_exit(0.25, 3.0, sym=True, rtol=1e-12, signed=True, perm=p1) != r['rows'][0]['d']
    assert r['rows'][1]['ld'] == MWE.a2_exit_ld(0.25, 3.0, sym=True, rtol=1e-12, perm=p1)
    ts, Lam = MWE.a2_line_lam(0.25, 3.0, 1e-12, 1.05 * max(t for q in r['rows'] for t in (q['d'][0], q['ld'][0])))
    for q in r['rows']:
        for k_, t in (('lnSd', q['d'][0]), ('lnSld', q['ld'][0])):
            assert abs(q[k_] - (math.log(0.01 * 0.25 * MWE.S0) - np.interp(t, ts, Lam))) < 1e-12
    j = MWE.a2p_job(0.25, 3.0, 1e-12, k=2)                                  # round 130 (D3dy-6): the job itself at 3A, rtol 1e-12 —
    assert (j['seed'], j['k'], j['max_step'], j['u'], j['ss'], len(j['rows'])) == (MWE.SEED_A2P[0], 2, None, r['u'], r['ss'], 3)
    assert all(a[k_] == b[k_] for a, b in zip(j['rows'], r['rows']) for k_ in ('d', 'ld', 'nd', 'nld'))   # a seed's first orders
    assert all(abs(a[k_] - b[k_]) < 1e-6 for a, b in zip(j['rows'], r['rows']) for k_ in ('lnSd', 'lnSld', 'lag'))
    line = MWE.a2p_line(P)
    x = [q['lag'] for q in r['rows'][1:]]
    assert f"3A at rtol 1e-12: Λ {np.mean(x):.2f} ± {np.std(x, ddof=1) / 10:.2f} and " in line   # per seed, from the data
    st = MWE.a2p_stats([y for y in P if (y['wmul'], y['rtol']) == (3.0, 1e-12)])
    assert (f"pooled {st['mean']:.2f} ± {st['sem']:.2f} (over {n}, {st['t']:+.1f} sem; median {st['q']['med']:.2f}, interquartile "
            f"{st['q']['q1']:.2f} … {st['q']['q3']:.2f}, extremes {st['q']['lo']:.2f} … {st['q']['hi']:.2f} of {n})") in line
    xa = [q['lag'] for y in P for q in y['rows'][1:]]                        # round 130 (D3dy-1): the pooled exponent, from the data
    ma, sa = np.mean(xa), np.std(xa, ddof=1) / math.sqrt(len(xa))
    assert len(xa) == 4 * n and abs(ma / lnr - 1) < 3 * sa / lnr
    assert f"over all {4 * n} draws Λ = {ma:.3f} ± {sa:.3f}, the exponent Λ/(11 ln 2) = {ma / lnr:.3f} ± {sa / lnr:.3f}" in line
    ts_ = [(MWE.a2p_stats(y)['t'], y) for y in P]
    tm, ym = max(ts_, key=lambda z: abs(z[0]))
    assert f"per-seed means' largest offset {tm:+.1f} sem ({ym['wmul']:g}A at rtol {ym['rtol']:g}, seed {ym['seed']})" in line
    assert "u the unit roundoff (2^{-53} in double, 2^{-64} in long double)" in line and "11 ln 2 = 7.62" in line
    assert "so the seed scales with the arithmetic's unit roundoff at DOP853's steps" in line
    assert "(drawn in round 130, after round 125's rule was set: round 125's stored draws, seed 125, were its referee's own)" in line
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    ml = next(x for x in md.splitlines() if x.startswith("The symmetric samples' exit as a draw (A2P"))
    cut = lambda t_: t_.split("; the rule on all of each cell's draws")[0].split('3A at rtol 1e-10:')[1]   # A2R joins the md's rule
    assert cut(ml) == cut(line) and "of `a2_sym_lag`'s, along the" in ml
    bad = [dict(x, rows=[dict(q, lag=q['lag'] + 3.0) for q in x['rows']]) if i == 3 else x for i, x in enumerate(P)]
    assert "so the seed scales" not in MWE.a2p_line(bad)                      # one seed's mean 3 above: its cell's exponent +0.2, none
    bad = [dict(x, rows=[dict(q, lnSld=q['lnSld'] - 10.0) for q in x['rows']]) if i == 0 else x for i, x in enumerate(P)]
    assert "so the seed scales" not in MWE.a2p_line(bad)                      # long-double seeds at e⁻⁵ of their mean: below 0.05
    one = [dict(x, rows=[dict(q, ld=(None, 0), lnSld=None, lag=None) if j_ == 1 else q for j_, q in enumerate(x['rows'])]) for x in P]
    assert MWE.a2p_stats(one[0])['n'] == 99 and MWE.a2p_stats(one[0])['up'][1][1] == 99 and "over 99 and 99" in MWE.a2p_line(one)
    zero = [dict(x, rows=[dict(q, ld=(None, 0), lnSld=None, lag=None) if j_ == 0 else q for j_, q in enumerate(x['rows'])]) for x in P]
    assert MWE.a2p_stats(zero[0])['n'] == 100 and MWE.a2p_stats(zero[0]) == MWE.a2p_stats(P[0])   # the code's order is not a draw of it
    # D3dt-2: no long-double exit — `a2_seed` returns four values and the exit line prints without the seed clause
    monkeypatch.setattr(MWE, 'a2_exit_ld', lambda *a, **k: ((None, 0), None, None) if k.get('dense') else (None, 0))
    calls = []
    monkeypatch.setattr(MWE, 'a2_exit', lambda *a, **k: calls.append(k) or (0.051 if (k.get('method'), k.get('rtol')) == ('BDF', 1e-12)
                                                                             else 0.05, -1))
    assert MWE.a2_seed(0.25, 3.0) == (None, (None, 0), (None, 0), None)
    ex = MWE.a2_exits(0.25, (3.0,))[3.0]                                       # D3dt-7: the BDF entry at its own rtol
    assert ex['imp'][2] == ('BDF at rtol 1e-12', (0.051, -1)) and ex['ld']['lag'] == [None, None] and ex['ld']['Lam'] is None
    t = '\n'.join(MWE.a2_lines(raw['A2'], raw['A2E'], P))
    assert ("(`a2_exit_ld`; referee round 112): 3A: on A's samples —, — (rtol 1e-10, 1e-12), on the symmetric samples —, —; 5A: on A's "
            "samples —, — (rtol 1e-10, 1e-12), on the symmetric samples —, —.") in t
    assert "the data's drift carried" not in t and "started at c* + S" not in t and line in t.splitlines()
    # D3dt-7: the mirror block's idle count and narrowest formed kernel (H σ* 0.06, equal); D3dt-6: no "over min A" at equal amplitudes
    sp = '\n'.join(MWE.split_summary(H=raw['H']))
    assert "(σ ≥ h; the narrowest 0.041, the widest narrower one 0.001) odd by 1.2e-09 on and between the samples (9 formed, the " \
           "narrowest of them 0.0496; 2 idle); the formed threshold's margin: the largest unformed |w| 0.0098 (over max A, against" in sp
    assert "over min A" not in sp
    # D3dt-5: at A2E's held states the normal rate dominates the tangential ones; while the width falls σ contracts faster than c
    held = {r_['wmul']: r_['rows'][-1] for r_ in raw['A2E'] if abs(r_['rows'][-1]['dc']) < 1e-12}
    lam3, tg3 = MWE.a2_tangent(0.25, held[3.0]['w'], held[3.0]['sig'])
    assert abs(lam3 - 628.85) < 0.01 and abs(tg3[0] + 2875.1) < 0.1 and abs(tg3[1] + 6.169e-4) < 1e-6
    neg = [x for r_ in raw['A2E'] if r_['wmul'] == 5.0 for x in r_['rows'] if abs(x['dc']) < 1e-12 and x['jcc'] <= 0]
    assert len(neg) == 8 and all(MWE.a2_tangent(0.25, x['w'], x['sig'])[1][0] < x['jcc'] < 0 for x in neg)
    assert ("λ = 628.9 (the (α, σ) block's eigenvalues -0.00062 and -2875: λ above both, normally hyperbolic, for every r as both are "
            "negative; referee rounds 125, 130)") in t
    assert "λ = 1.05e+04 (the (α, σ) block's eigenvalues -0.0012 and -1.446e+04: λ above both, normally hyperbolic, for every r" in t
    assert ("8 outputs: stable in c; the (α, σ) block's eigenvalues there -886 … -644 and +0.013 … +0.022, the width contracting faster "
            "than c (∂ċ/∂c above the smallest at every output): not normally hyperbolic, referee rounds 125, 130)") in t


def test_round130_the_capped_steps_the_verdicts_and_the_printed_agreement(monkeypatch):
    """Referee round 130 (D3dy-2, -4, -6): A2H — A2P's first orders of its first seed with DOP853's step capped at f·h₀, h₀ the double's
    mean step over A2P's draws — is the code's configuration, steps at its cap in both arithmetics and reproduces from the job (5A, the
    first cap); its line's shifts are of the means over the same orders at DOP853's own steps (round 133, D3eb-3: the sems in
    quadrature, the orders' draws independent); the normal-hyperbolicity verdicts are computed from the printed rates (plants flip
    them); A2P's "within" is the largest of the code order's differences."""
    import math
    import pytest
    if np.finfo(np.longdouble).eps > 1e-18:
        pytest.skip('no extended long double on this platform')
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    P, H = raw['A2P'], raw['A2H']
    h0 = {w: np.mean([q['d'][0] / q['nd'] for x in P if (x['wmul'], x['rtol']) == (w, MWE.RTOL_A2H) for q in x['rows'][1:]])
          for w in MWE.WMUL_A2P}
    assert all(abs(MWE.a2p_h0(P, w, MWE.RTOL_A2H) / h0[w] - 1) < 1e-14 for w in MWE.WMUL_A2P)
    assert [(r['wmul'], r['rtol'], r['k'], r['seed'], r['max_step']) for r in H] == [
        (w, MWE.RTOL_A2H, MWE.K_A2H, MWE.SEED_A2P[0], f * MWE.a2p_h0(P, w, MWE.RTOL_A2H)) for w in MWE.WMUL_A2P for f in MWE.FRAC_A2H]
    assert MWE.FRAC_A2H == (0.25, 0.0625) and MWE.K_A2H <= MWE.K_A2P and MWE.RTOL_A2H in MWE.RTOL_A2P
    line = MWE.a2h_line(H, P)
    for r in H:
        base = next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (r['wmul'], r['rtol'], r['seed']))
        st = MWE.a2p_stats(r)
        assert len(r['rows']) == r['k'] + 1 and all(q['d'][0] is not None and q['ld'][0] is not None for q in r['rows'])
        assert all(abs(h / r['max_step'] - 1) < 0.02 for _, h in st['steps'])        # under the cap both arithmetics step at it
        assert all(q['nd'] > b['nd'] and q['nld'] > b['nld'] for q, b in zip(r['rows'], base['rows']))   # more steps than DOP853's own
        d = []
        for k_ in ('lnSd', 'lnSld'):                                                  # the shift of the means over the same orders
            a_, b_ = (np.array([q[k_] for q in x['rows'][1:]]) for x in (r, base))
            v = a_ - b_[:len(a_)]
            d.append(MWE.a2p_shift(r, base, k_))
            assert len(v) == r['k'] and abs(d[-1][0] - np.mean(v)) < 1e-12 and abs(d[-1][2] - np.corrcoef(a_, b_[:len(a_)])[0, 1]) < 1e-12
            assert abs(d[-1][1] - math.sqrt(np.var(a_, ddof=1) / len(a_) + np.var(b_[:len(a_)], ddof=1) / len(a_))) < 1e-12
        f = r['max_step'] / h0[r['wmul']]
        bz = [(m - 0.5 * math.log(f)) / e for m, e, _ in d]                          # round 133 (D3eb-3): against Brouwer's ½ ln f
        assert (f"its mean moved by {d[0][0]:+.2f} ± {d[0][1]:.2f} and {d[1][0]:+.2f} ± {d[1][1]:.2f} (Brouwer's {0.5 * math.log(f):+.2f}: "
                f"{bz[0]:+.2f} and {bz[1]:+.2f} sem from it), Λ") in line
        assert f"at most {r['max_step'] / h0[r['wmul']]:g}·h₀: mean step {st['steps'][0][1]:.3e} and {st['steps'][1][1]:.3e}" in line
        assert f"Λ {st['mean']:.2f} ± {st['sem']:.2f} ({st['t']:+.1f} sem from 11 ln 2" in line
    r = H[MWE.WMUL_A2P.index(5.0) * len(MWE.FRAC_A2H)]                             # the job with a cap: its first orders (5A, h₀/4)
    j = MWE.a2p_job(0.25, 5.0, MWE.RTOL_A2H, k=1, max_step=r['max_step'])
    assert j['max_step'] == r['max_step'] and all(a[k_] == b[k_] for a, b in zip(j['rows'], r['rows']) for k_ in ('d', 'ld', 'nd', 'nld'))
    sh = [MWE.a2p_shift(r, next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (r['wmul'], r['rtol'], r['seed'])), k_)
          for r in H for k_ in ('lnSd', 'lnSld')]
    big = sum(1 for m, e, _ in sh if abs(m) > 3 * e)
    assert (f"{big} of {len(sh)} beyond 3 sem" in line) == (big > 0) and all(abs(c) < 0.3 for _, _, c in sh)
    lc = [MWE.a2p_stats(r)['mean'] for r in H]                                     # Λ at equal steps (the caps)
    assert (f"there Λ is {min(lc):.2f} … {max(lc):.2f} (the exponent {min(lc) / (11 * math.log(2)):.3f} … "
            f"{max(lc) / (11 * math.log(2)):.3f}, {sum(abs(MWE.a2p_stats(r)['t']) >= 3 for r in H)} of {len(H)} beyond 3 sem") in line
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read()
    assert line in md.splitlines()
    # D3dy-6: the regeneration path — part_A2P's and part_A2H's jobs are the raw's records' configuration
    cfg = lambda xs: [(x['ratio'], x['wmul'], x['rtol'], x['k'], x['seed'], x['max_step']) for x in xs]
    job = lambda ratio, wmul, rtol, k=MWE.K_A2P, seed=MWE.SEED_A2P[0], max_step=None: dict(ratio=ratio, wmul=wmul, rtol=rtol, k=k,
                                                                                           seed=seed, max_step=max_step)
    monkeypatch.setattr(MWE, 'a2p_job', job)
    assert cfg(MWE.part_A2P()) == cfg(P) and cfg(MWE.part_A2H(P)) == cfg(H)
    monkeypatch.undo()
    # A2P's "within": the largest difference of the code order's Λ from `a2_sym_lag`'s (two cells planted 0.01 and 0.03 off)
    p0 = {c: [x['rows'][0]['lag'] for x in P if (x['wmul'], x['rtol']) == c] for c in ((3.0, 1e-12), (5.0, 1e-10))}
    lg = {(3.0, 1e-12): p0[(3.0, 1e-12)][0] + 0.01, (5.0, 1e-10): p0[(5.0, 1e-10)][0] - 0.03}
    dl = max(abs(v - lg[c]) for c, vs in p0.items() for v in vs)
    assert abs(dl - 0.03) < 1e-4 and f"(its Λ within {dl:.1e} of `a2_sym_lag`'s" in MWE.a2p_line(P, lg)
    # D3dy-4: the verdicts from the rates
    assert MWE.nh_verdict(628.9, [-2875.0, -6.2e-4]) == (True, True) and MWE.nh_verdict(628.9, [-2875.0, 1.0]) == (True, False)
    assert MWE.nh_verdict(628.9, [-1.0, 700.0]) == (False, False) and MWE.nh_verdict(-300.0, [-886.0, 0.013]) == (False, False)
    assert MWE.nh_verdict(-300.0, [-200.0, 0.013]) == (True, False) and MWE.nh_verdict(-300.0, [0.1, 0.2]) == (True, True)
    assert MWE.nh_verdict(0.0, [0.1, 0.2]) == (False, False) == MWE.nh_verdict(0.0, [-0.2, -0.1])   # round 133 (D3eb-4): lam = 0
    assert MWE.nh_verdict(-300.0, [0.0, 0.2]) == (True, False) and MWE.nh_verdict(628.9, [-1.0, 0.0]) == (True, False)   # tg at 0
    fake = {w: dict(a=[(0.05, -1), (0.05, -1)], imp=[('Radau', (0.06, -1)), ('BDF', (None, 0))], sym=[(0.05, 1), (0.05, 1)], kz=9.0,
                    symoff=(0.05, -1), summ=dict(fsum=((0.05, -1), (None, 0)), rev=((0.05, -1), (0.05, -1))),
                    ld=dict(a=[(0.05, -1), (0.05, -1)], sym=[(0.06, 1), (0.06, -1)], seed=-7e-17, symS=(0.05, -1), ulp=3.47e-18))
            for w in (3.0, 5.0)}
    monkeypatch.setattr(MWE, 'a2_exits', lambda *a, **k: fake)
    got = lambda: '\n'.join(MWE.a2_lines(raw['A2'], raw['A2E']))
    t = got()
    assert t.count("λ above both, normally hyperbolic, for every r as both are negative;") == 2
    # D3dy-5: "apart by" is the relative difference (planted exits 0.05 and 0.0501 at the two tolerances: 2.0e-03, not 1.0e-04)
    fake[3.0]['ld']['a'] = [(0.05, -1), (0.0501, -1)]
    assert f"(apart by {abs(0.05 / 0.0501 - 1):.1e} relatively)" in got() and "(apart by 2.0e-03 relatively)" in got()
    fake[3.0]['ld']['a'] = [(0.05, -1), (0.05, -1)]
    # the report carries A2H's line after A2P's, and `--from-raw` reports every part (KEYS)
    rp = MWE.report(A2=raw['A2'], A2E=raw['A2E'], A2P=P, A2H=H).splitlines()
    assert rp.index(MWE.a2h_line(H, P)) == rp.index(MWE.a2p_line(P)) + 1
    assert MWE.KEYS == ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'H3', 'HL', 'HP', 'GL', 'GT', 'A2', 'HLL', 'SQ', 'A2E', 'A2P', 'A2H',
                        'A2R', 'A2T')
    assert "the width contracting faster than c (∂ċ/∂c above the smallest at every output): not normally hyperbolic" in t
    tg = MWE.a2_tangent
    monkeypatch.setattr(MWE, 'a2_tangent', lambda *a: (lambda l, e: (l, sorted([e[0], 2 * abs(l)])))(*tg(*a)))   # a tangential rate 2|λ|
    t = got()
    assert t.count("eigenvalues") >= 3 and t.count("a tangential rate above λ: not normally hyperbolic;") == 2
    assert "the width contracting faster than c (∂ċ/∂c above the smallest at every output): not normally hyperbolic" in t
    monkeypatch.setattr(MWE, 'a2_tangent', lambda *a: (lambda l, e: (l, sorted([0.5 * l, e[1]])))(*tg(*a)))     # c the fastest: halved
    t = got()
    assert "c contracting faster than both at every output: normally hyperbolic, referee rounds 125, 130" in t
    assert "λ above both, normally hyperbolic; referee" in t and "for every r as both are negative" not in t   # one rate positive
    neg = [x for r_ in raw['A2E'] if r_['wmul'] == 5.0 for x in r_['rows'] if abs(x['dc']) < 1e-12 and x['jcc'] <= 0]
    first = min(x['sig'] for x in neg)
    monkeypatch.setattr(MWE, 'a2_tangent', lambda ratio, w, sig: (lambda l, e: (l, sorted([0.5 * l, e[1]]) if sig == first else e))(
        *tg(ratio, w, sig)))
    assert "normally hyperbolic at 1 of 8 outputs, referee rounds 125, 130" in got()


def test_round133_the_interval_rule_the_replications_the_shifts_and_the_edges(monkeypatch):
    """Referee round 133 (D3eb-1 … D3eb-6): A2P's conclusion rests on every cell's exponent's 99% interval inside [0.9, 1.1] and every
    geometric mean within 0.05–20 uσ* (both edges of each planted), a cell beyond 3 sem of 11 ln 2 named as a prefactor difference
    without withdrawing it; the pooled exponent's sem over orders; A2R — further seeds at A2P's farthest cell, the double at the long
    double's atol, A2H's farthest cap on the second seed — is the code's configuration, reproduces from the job, and its line's numbers
    and verdicts are the data's (plants flip them); A2H's shifts against Brouwer's law and the exponent's spread over every
    configuration; A2H's base truncated to its orders when K_A2H < K_A2P."""
    import math
    import pytest
    if np.finfo(np.longdouble).eps > 1e-18:
        pytest.skip('no extended long double on this platform')
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    P, H, R = raw['A2P'], raw['A2H'], raw['A2R']
    lnr = 11 * math.log(2)
    cells = [(w, q) for w in MWE.WMUL_A2P for q in MWE.RTOL_A2P]
    st = {c: MWE.a2p_stats([r for r in P if (r['wmul'], r['rtol']) == c]) for c in cells}
    line = MWE.a2p_line(P)
    assert "so the seed scales with the arithmetic's unit roundoff at DOP853's steps" in line
    assert "Λ = 11 ln 2 itself asks the two arithmetics' prefactors to be equal — no cell beyond 3 sem of it (the largest +2.7, 3A" in line
    assert "every cell's pooled mean within 3 sem" not in line and "round 125's rule, here" not in line
    # D3eb-5: the pooled exponent's sem over orders (each seed's order j the same permutation in every cell)
    bl = {}
    for r in P:
        for j, q in enumerate(r['rows'][1:]):
            bl.setdefault((r['seed'], j), []).append(q['lag'])
    v = [np.mean(x) for x in bl.values()]
    assert len(v) == len(MWE.SEED_A2P) * MWE.K_A2P and all(len(x) == len(cells) for x in bl.values())
    assert f"(± {np.std(v, ddof=1) / math.sqrt(len(v)) / lnr:.4f} over the {len(v)} orders, each order's cells averaged" in line
    # D3eb-1, -6: the interval rule's edges — one cell's draws shifted so that its interval's end lies 1e-4 inside or outside
    c = (5.0, 1e-12)
    s = st[c]
    shift = lambda P_, d_: [dict(x, rows=[dict(q, lag=q['lag'] + d_) for q in x['rows']]) if (x['wmul'], x['rtol']) == c else x
                            for x in P_]
    for end, sign in ((MWE.EXP_A2P[1], +1), (MWE.EXP_A2P[0], -1)):    # round 137 (D3ef-2): the band [0.75, 1.25]
        for eps, keep in ((-1e-4, True), (1e-4, False)):
            d_ = (end + sign * eps) * lnr - s['mean'] - sign * MWE.Z99 * s['sem']
            e = MWE.a2p_exp(MWE.a2p_stats([x for x in shift(P, d_) if (x['wmul'], x['rtol']) == c]), lnr)
            assert abs((e[2] if sign > 0 else e[1]) - (end + sign * eps)) < 1e-9
            assert ("so the seed scales" in MWE.a2p_line(shift(P, d_))) == keep
    for t_, named in ((3.5, True), (2.9, False)):                     # a cell 3.5 sem off: named, the conclusion kept (its interval in)
        d_ = lnr + t_ * s['sem'] - s['mean']
        l_ = MWE.a2p_line(shift(P, d_))
        assert MWE.a2p_exp(MWE.a2p_stats([x for x in shift(P, d_) if (x['wmul'], x['rtol']) == c]), lnr)[2] < MWE.EXP_A2P[1]
        assert "so the seed scales" in l_ and ("1 of 4 cells beyond 3 sem of it (5A at rtol 1e-12, +3.5): there they differ" in l_) == named
        assert ("no cell beyond 3 sem of it" in l_) == (not named)
    g = st[c]['gm'][0]                                                 # the geometric mean's upper edge (20 uσ*)
    for f, keep in ((0.99, True), (1.01, False)):
        pl = [dict(x, rows=[dict(q, lnSd=q['lnSd'] + math.log(f * 20 / g)) for q in x['rows']]) if (x['wmul'], x['rtol']) == c else x
              for x in P]
        assert ("so the seed scales" in MWE.a2p_line(pl)) == keep
    # A2R: the code's configuration (further seeds, disjoint from A2P's and round 125's/130's referees'), from the jobs
    assert MWE.CELL_A2R == max(cells, key=lambda c_: abs(st[c_]['t']))
    assert (MWE.CAP_A2R[0], MWE.CAP_A2R[1]) == max(((r['wmul'], round(r['max_step'] / MWE.a2p_h0(P, r['wmul'], r['rtol']), 12)) for r in H),
                                                    key=lambda k_: abs(MWE.a2p_stats(next(
                                                        r for r in H if (r['wmul'], round(r['max_step'] / MWE.a2p_h0(P, r['wmul'], r['rtol']),
                                                                                          12)) == k_))['t']))
    assert not (set(MWE.SEED_A2R) & (set(MWE.SEED_A2P) | {125, 130, 2026})) and MWE.ATOL_A2R == 1e-14
    cfg = lambda xs: [(x['ratio'], x['wmul'], x['rtol'], x['k'], x['seed'], x['max_step'], x['atol_d']) for x in xs]
    w, q = MWE.CELL_A2R
    ss = 0.25 * MWE.S0
    h016 = MWE.CAP_A2R[1] * MWE.a2p_h0(P, MWE.CAP_A2R[0], MWE.RTOL_A2H)
    assert cfg(R) == ([(0.25, w, q, MWE.K_A2P, s_, None, None) for s_ in MWE.SEED_A2R]
                      + [(0.25, w, q, MWE.K_A2P, s_, None, 1e-14 * ss) for s_ in MWE.SEED_A2P + MWE.SEED_A2R]
                      + [(0.25, MWE.CAP_A2R[0], MWE.RTOL_A2H, MWE.K_A2H, MWE.SEED_A2P[1], h016, None)])
    h = next(r for r in H if r['wmul'] == MWE.CAP_A2R[0] and r['max_step'] == h016)
    assert h['seed'] == MWE.SEED_A2P[0] and all(len(r['rows']) == r['k'] + 1 for r in R)
    job = lambda ratio, wmul, rtol, k=MWE.K_A2P, seed=MWE.SEED_A2P[0], max_step=None, atol_d=None: dict(
        ratio=ratio, wmul=wmul, rtol=rtol, k=k, seed=seed, max_step=max_step, atol_d=atol_d)
    monkeypatch.setattr(MWE, 'a2p_job', job)
    assert cfg(MWE.part_A2R(P)) == cfg(R)
    monkeypatch.undo()
    more = [r for r in R if r['max_step'] is None and r['atol_d'] is None]
    mat = [r for r in R if r['atol_d'] is not None]
    cap = [r for r in R if r['max_step'] is not None]
    for r_, kw in ((more[0], {}), (mat[0], dict(atol_d=1e-14 * ss)), (cap[0], dict(max_step=h016))):   # the job's first orders
        j = MWE.a2p_job(0.25, r_['wmul'], r_['rtol'], k=1, seed=r_['seed'], **kw)
        assert all(a[k_] == b[k_] for a, b in zip(j['rows'], r_['rows']) for k_ in ('d', 'ld', 'nd', 'nld'))
    md = {r['seed']: r for r in more}
    for r in mat:                                                      # the matched arm changes the double only
        ref = md.get(r['seed']) or next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (w, q, r['seed']))
        assert all(a['ld'] == b['ld'] and a['nld'] == b['nld'] for a, b in zip(r['rows'], ref['rows']))
        assert np.mean([a['nd'] for a in r['rows'][1:]]) > np.mean([b['nd'] for b in ref['rows'][1:]])
    # A2R's line: its numbers are the data's
    base = [r for r in P if (r['wmul'], r['rtol']) == (w, q)]
    lr = MWE.a2r_line(R, P, H)
    x = np.array([q_['lag'] for r in base + more for q_ in r['rows'][1:]])
    m, e = x.mean(), x.std(ddof=1) / math.sqrt(x.size)
    assert x.size == 600 and f"over all {x.size} Λ = {m:.3f} ± {e:.3f} ({(m - lnr) / e:+.2f} sem) from 11 ln 2" in lr
    assert (f"the exponent {m / lnr:.3f} (99% interval {(m - MWE.Z99 * e) / lnr:.3f} … {(m + MWE.Z99 * e) / lnr:.3f})" in lr
            and 0.9 < (m - MWE.Z99 * e) / lnr and (m + MWE.Z99 * e) / lnr < 1.1)
    lp = [np.array([q_[k_] - math.log(r['u'][i] * r['ss']) for r in base + more for q_ in r['rows'][1:]]) for i, k_ in
          enumerate(('lnSd', 'lnSld'))]
    assert (f"ln(S/uσ*) {lp[0].mean():+.2f} ± {lp[0].std(ddof=1) / math.sqrt(600):.2f} in double and {lp[1].mean():+.2f} ± "
            f"{lp[1].std(ddof=1) / math.sqrt(600):.2f} in long double — their difference Λ − 11 ln 2 = {m - lnr:+.2f} ± {e:.2f}") in lr
    assert abs(lp[0].mean() - lp[1].mean() - (m - lnr)) < 1e-9
    y = np.array([q_['lag'] for r in mat for q_ in r['rows'][1:]])
    my, ey = y.mean(), y.std(ddof=1) / math.sqrt(y.size)
    assert y.size == 600 and f"on the same 600 orders (seeds 131, 132, 133, 134, 135, 136) Λ = {my:.3f} ± {ey:.3f} ({(my - lnr) / ey:+.2f} sem)" in lr
    z = np.array([q_['lag'] for r in (h, cap[0]) for q_ in r['rows'][1:]])
    mz, ez = z.mean(), z.std(ddof=1) / math.sqrt(z.size)
    assert f"over both seeds' {z.size} Λ = {mz:.3f} ± {ez:.3f} ({(mz - lnr) / ez:+.2f} sem)" in lr
    sm_ = MWE.a2p_stats(more)['t']                                    # round 137 (D3ef-4): the verdicts on the further draws alone
    s1_ = MWE.a2p_stats(cap[0])['t']
    for t_, txt in ((sm_, "differ at this cell, by about e^{"), (s1_, ": replicated, the arithmetics' runs' prefactors differ at equal steps")):
        assert (txt in lr) == (abs(t_) >= 3)
    assert ": not the atol" in lr
    # plants flip A2R's verdicts: the further seeds, the matched arm, the cap moved to 11 ln 2
    to = lambda rs, d_: [dict(r, rows=[dict(q_, lag=q_['lag'] + d_) for q_ in r['rows']]) for r in rs]
    Pc = [dict(r, rows=[dict(q_, lag=q_['lag'] + lnr - m) for q_ in r['rows']]) if (r['wmul'], r['rtol']) == (w, q) else r for r in P]
    l2 = MWE.a2r_line(to(more, lnr - m) + mat + cap, Pc, H)
    assert ": on the further draws alone within 3 sem (" in l2 and "differ at this cell" not in l2
    l3 = MWE.a2r_line(more + to(mat, lnr - my) + cap, P, H)          # the matched arm moved by ≈ 4 sem of the paired shift
    assert ": the atol moves Λ" in l3 and ": not the atol" not in l3
    l4 = MWE.a2r_line(more + mat + to(cap, lnr - mz), P, [dict(r, rows=[dict(q_, lag=q_['lag'] + lnr - mz) for q_ in r['rows']])
                                                         if r is h else r for r in H])
    assert ", not replicated beyond 3 sem by itself" in l4 and "differ at equal steps too" not in l4
    # D3eb-2, -3: A2H's line — no "Λ compares the arithmetics at equal steps"; Brouwer's law; the exponent's systematic spread
    lh = MWE.a2h_line(H, P)
    assert "Λ compares the arithmetics at equal steps" not in lh and "per order" not in lh
    xs = [st[c_]['mean'] / lnr for c_ in cells] + [MWE.a2p_stats(r)['mean'] / lnr for r in H]
    al = MWE.a2p_stats(P)
    assert (f"over every configuration (4 cells at DOP853's own steps, A2P; 4 caps) the exponent Λ/(11 ln 2) spans {min(xs):.3f} … "
            f"{max(xs):.3f} — a range of 8 estimates") in lh and "systematic spread" not in lh
    nb = []
    for r in H:
        b_ = next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (r['wmul'], r['rtol'], r['seed']))
        f = r['max_step'] / MWE.a2p_h0(P, r['wmul'], r['rtol'])
        for i, k_ in enumerate(('lnSd', 'lnSld')):
            mm, ee, _ = MWE.a2p_shift(r, b_, k_)
            if abs(mm - 0.5 * math.log(f)) > 3 * ee:
                nb.append(f"{r['wmul']:g}A at {f:g}·h₀ in {('double', 'long double')[i]}")
    assert f"{len(nb)} of 8 beyond 3 sem of Brouwer's law (" + ', '.join(nb) + ")" in lh and len(nb) == 2
    assert "1 of 4 beyond 3 sem of 11 ln 2: 5A at 0.0625·h₀, +4.7 — the two arithmetics' prefactors differ at equal steps too" in lh
    # D3eb-6: K_A2H < K_A2P — the base is A2H's orders, not A2P's record
    Hs = [dict(r, k=50, rows=r['rows'][:51]) for r in H]
    l5 = MWE.a2h_line(Hs, P)
    b0 = next(x for x in P if (x['wmul'], x['rtol'], x['seed']) == (3.0, MWE.RTOL_A2H, MWE.SEED_A2P[0]))
    bt, bf = MWE.a2p_stats(dict(b0, rows=b0['rows'][:51])), MWE.a2p_stats(b0)
    lam = lambda s_: f"Λ {s_['mean']:.2f} ± {s_['sem']:.2f} ({s_['t']:+.1f} sem from 11 ln 2"
    assert lam(bt) != lam(bf) and f"3A (h₀ 8.71e-04): at DOP853's own steps" in l5 and "random orders of seed 131" in l5
    assert l5.split("3A (h₀ 8.71e-04): at DOP853's own steps")[1].split('; at most')[0].endswith(
        f"{lam(bt)}, the exponent {bt['mean'] / lnr:.3f} ± {bt['sem'] / lnr:.3f})")
    d0 = MWE.a2p_shift(Hs[0], dict(b0, rows=b0['rows'][:51]), 'lnSd')
    assert f"its mean moved by {d0[0]:+.2f} ± {d0[1]:.2f} and" in l5
    # the md carries A2R's line after A2H's
    mdl = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read().splitlines()
    ml = next(x for x in mdl if x.startswith("The symmetric samples' exit as a draw (A2P"))     # the md's line carries `a2_sym_lag`
    assert lr in mdl and mdl.index(lr) == mdl.index(lh) + 1 and mdl.index(lh) == mdl.index(ml) + 1
    assert ml.split("; the rule on all of each cell's draws")[0].split('3A at rtol 1e-10:')[1] == line.split(
        "; the rule on all of each cell's draws")[0].split('3A at rtol 1e-10:')[1]
    rp = MWE.report(A2=raw['A2'], A2E=raw['A2E'], A2P=P, A2H=H, A2R=R).splitlines()
    assert rp.index(lr) == rp.index(lh) + 1


def test_round137_the_traps_three_arithmetics_the_gate_s_input_and_the_replications(monkeypatch):
    """Referee round 137 (D3ef-1 … D3ef-7): `_DTBumps` is `_LDBumps` bitwise at long double and `BumpsCW` bitwise at double (c_ref = 0),
    float32 a third arithmetic; A2T — the scans of ċ(δ) and the fixed-step RK4 arms — is the code's configuration and reproduces from
    the code; its line's numbers and verdicts are the data's (plants flip them); A2P's rule is gated on all of each cell's draws (A2P with
    A2R's further seeds: the gate's input asserted, the band's edges planted on the pooled 3A cell, the geometric mean's open band);
    A2R's verdicts on the further draws alone (t to two decimals, the steps to two, the Stouffer combination, the paired atol shift);
    A2H's DerSimonian–Laird heterogeneity; the md carries A2T's line after A2R's."""
    import math
    import pytest
    if np.finfo(np.longdouble).eps > 1e-18:
        pytest.skip('no extended long double on this platform')
    raw = MWE.E.load(os.path.join(MWE.OUT, 'moving_width_raw.npz'))
    P, H, R, T = raw['A2P'], raw['A2H'], raw['A2R'], raw['A2T']
    lnr = 11 * math.log(2)
    # D3ef-3: the dtype-generic right-hand side, bitwise at long double and double
    L, ss, d = MWE.law(), 0.25 * MWE.S0, MWE.XS[1] - MWE.XS[0]
    xs = d * np.arange(-(MWE.XS.size // 2), MWE.XS.size // 2 + 1)
    xp = xs[np.random.default_rng(131).permutation(xs.size)]
    f = MWE.A_A * np.exp(-xp ** 2 / (2 * ss * ss))
    rl = RN.critical_rhs_nd(L, MWE._LDBumps(xp, MWE.S0, 0.0), f.astype(np.longdouble))
    rp = RN.critical_rhs_nd(L, MWE._DTBumps(xp, MWE.S0, 0.0, np.longdouble), f.astype(np.longdouble))
    rd = RN.critical_rhs_nd(L, MW.BumpsCW(xp, MWE.S0), f)
    rpd = RN.critical_rhs_nd(L, MWE._DTBumps(xp, MWE.S0, 0.0, np.float64), f)
    rf = RN.critical_rhs_nd(L, MWE._DTBumps(xp, MWE.S0, 0.0, np.float32), f.astype(np.float32))
    al = float(RT_alpha(L, 3.0))
    for z in ([al, 0.0, ss], [al * 1.1, 3e-18, ss * 0.5], [al * 0.9, -2e-17, ss * 0.2], [al, 1e-4, ss]):
        z = np.array(z)
        assert np.array_equal(rl(0, z), rp(0, z)) and np.array_equal(rd(0, z), rpd(0, z))
        assert not np.array_equal(rd(0, z), rf(0, z)) and np.allclose(rf(0, z)[[0, 2]], rd(0, z)[[0, 2]], rtol=1e-4, atol=0)
    F, z0, s_ = MWE.a2t_rhs(3.0, np.float64, None)
    assert np.array_equal(F(0, z0), RN.critical_rhs_nd(L, MW.BumpsCW(xs, MWE.S0), A_f(xs, ss))(0, z0)) and s_ == ss
    # A2T is the code's configuration (fresh seeds) and its jobs are a2t_jobs'
    assert not ({MWE.SEED_A2T} | {x[3] for x in MWE.RK_A2T}) & (set(MWE.SEED_A2P) | set(MWE.SEED_A2R) | {125, 130, 137, 138, 2026})
    jobs = MWE.a2t_jobs(P)
    sc = [r for r in T if r['kind'] == 'scan']
    rk = [r for r in T if r['kind'] == 'rk']
    assert [(r['wmul'], tuple(x['s'] for x in r['rows'])) for r in sc] == [(w, ts) for w, ts in MWE.S_A2T]
    assert [(r['wmul'], r['frac'], r['k'], r['seed'], r['h']) for r in rk] == [(w, f_, k, s_, f_ * h0) for _, w, f_, k, s_, h0 in jobs[2:]]
    assert all(len(r['rows']) == r['k'] for r in rk) and all(h0 == MWE.a2p_h0(P, w, 1e-10) for _, w, _, _, _, h0 in jobs[2:])
    monkeypatch.setattr(MWE, 'a2t_job', lambda a: a)
    assert MWE.part_A2T(P) == jobs
    monkeypatch.undo()
    # the scans reproduce (one state, every arithmetic, the code's order) and the line's states are the long-double line flow's
    r0 = sc[0]
    sts = MWE.a2t_states(r0['wmul'], [row['s'] for row in r0['rows']])
    assert [(a, sg / ss) for a, sg in sts] == [(row['al'], row['sig']) for row in r0['rows']]
    for x in r0['rows'][0]['scans']:
        if x['order'] == 0:
            y = MWE.a2t_scan(r0['wmul'], r0['rows'][0]['al'], r0['rows'][0]['sig'] * ss, dict(MWE.DT_A2T)[x['arith']])
            assert y == {k_: x[k_] for k_ in y}
    assert {round(x['ulp'], 6) for r in sc for row in r['rows'] for x in row['scans']} == {round(2.0 ** -59 / (2.0 ** -53 * ss), 6)}
    r5 = next(r for r in sc if r['wmul'] == 5.0)                       # the job reproduces the 5A scans whole (states, orders, arithmetics)
    j5 = MWE.a2t_job(('scan', 5.0, dict(MWE.S_A2T)[5.0]))
    assert j5 == {k_: r5[k_] for k_ in j5} and [x['order'] for x in j5['rows'][0]['scans']] == [0, 1, 2] * 3
    # a planted field: the brackets' midpoints and the counts
    g = MWE.a2t_scan(3.0, r0['rows'][0]['al'], r0['rows'][0]['sig'] * ss, np.float64, nt=81, nu=4)
    assert g['att'] and all(abs(a * 10 - math.floor(a * 10) - 0.5) < 1e-9 for a in g['att'] + g['rep'])   # step 0.1 ulp: midpoints
    # the RK4 arms reproduce (the first order of 5A at 0.25·h₀, every arithmetic) and agree with a hand-written RK4 step
    a5 = next(r for r in rk if (r['wmul'], r['frac']) == (5.0, 0.25))
    p0 = np.random.default_rng(a5['seed']).permutation(MWE.XS.size)
    for nm, dt in MWE.DT_A2T:
        (s_e, dr), n_ = MWE.a2t_rk4(5.0, dt, p0, a5['h'])
        assert (s_e, dr, n_) == (a5['rows'][0][nm]['s'], a5['rows'][0][nm]['dir'], a5['rows'][0][nm]['n'])
    F, z, _ = MWE.a2t_rhs(5.0, np.float64, p0)
    h = a5['h']
    k1 = F(0, z); k2 = F(h / 2, z + h / 2 * k1); k3 = F(h / 2, z + h / 2 * k2); k4 = F(h, z + h * k3)
    z1 = z + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    (s1, _), n1 = MWE.a2t_rk4(5.0, np.float64, p0, h, smax=h / 2, level=abs(z1[1]) / ss)
    assert n1 == 1 and abs(s1 - h) < 1e-15
    assert MWE.a2t_rk4(5.0, np.float64, p0, h, smax=2.5 * h) == ((None, 0), 3)
    ts, Lam = MWE.a2_line_lam(0.25, 5.0, 1e-12, 1.05 * max(x['s'] for q in a5['rows'] for x in q.values() if x['s'] is not None))
    x0 = a5['rows'][0]['double']
    assert x0['lnS'] == float(math.log(1e-2 * ss) - np.interp(x0['s'], ts, Lam))
    j2 = MWE.a2t_job(('rk', 5.0, a5['frac'], 2, a5['seed'], a5['h'] / a5['frac']))   # the job's first two orders, its own line
    assert j2['h'] == a5['h'] and j2['u'] == a5['u'] and all(
        (a[nm]['s'], a[nm]['dir'], a[nm]['n']) == (b[nm]['s'], b[nm]['dir'], b[nm]['n']) for a, b in zip(j2['rows'], a5['rows']) for nm in a)
    ts2, Lam2 = MWE.a2_line_lam(0.25, 5.0, 1e-12, 1.05 * max(x['s'] for q in j2['rows'] for x in q.values()))
    assert all(x['lnS'] == float(math.log(1e-2 * ss) - np.interp(x['s'], ts2, Lam2)) for q in j2['rows'] for x in q.values())
    # A2T's line: the numbers are the data's
    lt = MWE.a2t_line(T)
    nms = [nm for nm, _ in MWE.DT_A2T]
    for r in rk:
        ok = [q for q in r['rows'] if all(q[nm]['s'] is not None for nm in nms)]
        X = np.array([math.log(r['u'][nm]) for nm in nms]); Xc = X - X.mean()
        gg = np.array([np.dot(Xc, [q[nm]['lnS'] for nm in nms]) / np.dot(Xc, Xc) for q in ok])
        cnt = [sum(1 for q in r['rows'] if q[nm]['s'] is not None) for nm in nms]
        assert f"at {r['frac']:g}·h₀ ({r['h']:.3e}, {r['k']} orders of seed {r['seed']}): exits {cnt[0]}, {cnt[1]}, {cnt[2]} of {r['k']}" in lt
        if len(ok) > 1:
            m, e = gg.mean(), gg.std(ddof=1) / math.sqrt(gg.size)
            assert (f"over the {len(ok)} orders where all three leave the three-point slope γ = {m:.4f} ± {e:.4f} (99% interval "
                    f"{m - MWE.Z99 * e:.4f} … {m + MWE.Z99 * e:.4f}, median {np.median(gg):.4f})") in lt
            pl = [(q['float32']['lnS'] - q['long double']['lnS']) / (X[0] - X[2]) for q in ok]
            assert f"float32/long double {np.mean(pl):.4f} ± {np.std(pl, ddof=1) / math.sqrt(len(pl)):.4f}" in lt
            cr = [(q['double']['lnS'] - X[1]) - (q['long double']['lnS'] - X[2]) for q in ok]
            assert f"ln(C_double/C_long double) {np.mean(cr):+.2f} ± {np.std(cr, ddof=1) / math.sqrt(len(cr)):.2f}" in lt
        for nm in nms:
            v = [q[nm]['lnS'] - math.log(r['u'][nm] * ss) for q in r['rows'] if q[nm]['s'] is not None]
            if v:
                assert f"{nm} s_e {np.mean([q[nm]['s'] for q in r['rows'] if q[nm]['s'] is not None]):.4f}" in lt
                assert f"ln(S/uσ*) {np.mean(v):+.2f}" in lt
    for r in sc:
        al_ = [x for row in r['rows'] for x in row['scans']]
        assert f"attracting zeros in {sum(1 for x in al_ if x['att'])} of {len(al_)} scans" in lt
    assert "attracting zeros in 27 of 27 scans" in lt and "attracting zeros in 0 of 27 scans" in lt
    assert ("at 3A the computed field's ulp-scale zeros decide: RK4 at 0.25·h₀ holds the line in double and long double to s = 0.5, at "
            "1·h₀ and 0.5·h₀ it leaves: the exit's existence and its seed are the integrator's") in lt
    gm5 = [MWE.a2t_rk_stats(r)['gam'] for r in rk if r['wmul'] == 5.0]
    dg = max(max(abs(g_[2] - 1), abs(g_[3] - 1)) for g_ in gm5)
    assert f"where the scans find no attracting zero (5A) three arithmetics at equal steps put γ's 99% interval inside [0.75, 1.25] at every arm, within {dg:.3f} of 1" in lt
    # plants: the 3A arm at 0.25·h₀ leaving in double and long double removes the 'held' verdict; float32's seeds raised at 5A put γ
    # out of the band; a 5A scan with an attracting zero moves 5A out of the scaling's arms
    lv = lambda q: {nm: (dict(x, s=0.1, dir=1, lnS=-30.0) if x['s'] is None else x) for nm, x in q.items()}
    T1 = [dict(r, rows=[lv(q) for q in r['rows']]) if r['kind'] == 'rk' and r['wmul'] == 3.0 and r['frac'] == 0.25 else r for r in T]
    assert "holds the line" not in MWE.a2t_line(T1) and "the exit's existence" not in MWE.a2t_line(T1)
    lvd = lambda q: dict(q, double=dict(q['double'], s=0.1, dir=1, lnS=-30.0))   # the double alone leaving: not held (both must hold)
    T1d = [dict(r, rows=[lvd(q) for q in r['rows']]) if r['kind'] == 'rk' and r['wmul'] == 3.0 and r['frac'] == 0.25 else r for r in T]
    assert "holds the line" not in MWE.a2t_line(T1d)
    nm3 = [nm for nm, _ in MWE.DT_A2T]                                 # γ planted at 1.17 at both 5A arms: inside [0.75, 1.25] only
    X3 = np.array([math.log(rk[0]['u'][nm]) for nm in nm3]); Xc3 = X3 - X3.mean()
    sh32 = {id(r): (1.17 - MWE.a2t_rk_stats(r)['gam'][0]) * np.dot(Xc3, Xc3) / Xc3[0] for r in rk if r['wmul'] == 5.0}
    upg = lambda q, d_: dict(q, float32=dict(q['float32'], lnS=None if q['float32']['lnS'] is None else q['float32']['lnS'] + d_))
    T4 = [dict(r, rows=[upg(q, sh32[id(r)]) for q in r['rows']]) if id(r) in sh32 else r for r in T]
    g4 = [MWE.a2t_rk_stats(r)['gam'] for r in T4 if r['kind'] == 'rk' and r['wmul'] == 5.0]
    assert all(abs(g_[0] - 1.17) < 1e-9 and 1.1 < g_[2] and g_[3] < 1.25 for g_ in g4)
    assert "three arithmetics at equal steps put γ's 99% interval inside [0.75, 1.25] at every arm" in MWE.a2t_line(T4)
    up = lambda q: dict(q, float32=dict(q['float32'], lnS=None if q['float32']['lnS'] is None else q['float32']['lnS'] + 10.0))
    T2 = [dict(r, rows=[up(q) for q in r['rows']]) if r['kind'] == 'rk' and r['wmul'] == 5.0 else r for r in T]
    assert "γ not inside the band at every arm without traps" in MWE.a2t_line(T2)
    tz = lambda row: dict(row, scans=[dict(x, att=[0.505]) if i == 0 else x for i, x in enumerate(row['scans'])])
    T3 = [dict(r, rows=[tz(row) for row in r['rows']]) if r['kind'] == 'scan' and r['wmul'] == 5.0 else r for r in T]
    l3 = MWE.a2t_line(T3)
    assert "attracting zeros in 3 of 27 scans" in l3 and "no arm with three arithmetics where the scans find no attracting zero" in l3
    # D3ef-2: the gate's input — A2P with A2R's further seeds at the cells they replicate — and its edges on the pooled 3A cell
    ok, gt = MWE.a2p_gate(P, R)
    assert ok and {c: sd for c, (_, _, sd) in gt.items()} == {
        (3.0, 1e-10): [131, 132, 133, 134, 135, 136], (3.0, 1e-12): [131, 132], (5.0, 1e-10): [131, 132], (5.0, 1e-12): [131, 132]}
    assert [s['n'] for s, _, _ in gt.values()] == [600, 200, 200, 200] and MWE.EXP_A2P == (0.75, 1.25) and MWE.EXP_A2P_133 == (0.9, 1.1)
    more = [r for r in R if r['max_step'] is None and r['atol_d'] is None]
    mat = [r for r in R if r['atol_d'] is not None]
    cap = [r for r in R if r['max_step'] is not None]
    c3 = MWE.CELL_A2R
    x6 = np.array([q['lag'] for r in [r for r in P if (r['wmul'], r['rtol']) == c3] + more for q in r['rows'][1:]])
    m6, e6 = x6.mean(), x6.std(ddof=1) / math.sqrt(x6.size)
    assert gt[c3][1] == MWE.a2p_exp(MWE.a2p_stats([r for r in P if (r['wmul'], r['rtol']) == c3] + more), lnr)
    assert abs(gt[c3][1][1] - (m6 - MWE.Z99 * e6) / lnr) < 1e-12 and abs(gt[c3][1][2] - (m6 + MWE.Z99 * e6) / lnr) < 1e-12
    md = open(os.path.join(MWE.OUT, 'moving_width.md'), encoding='utf-8').read().splitlines()
    ml = next(x for x in md if x.startswith("The symmetric samples' exit as a draw (A2P"))
    full = MWE.a2p_line(P, None, R)
    cut = lambda t_: t_.split("; the rule on all of each cell's draws")[1]
    assert cut(full) == cut(ml) and "3A at rtol 1e-10 over 600, 1.015 … 1.073" in ml
    sh = lambda rs, d_: [dict(r, rows=[dict(q, lag=q['lag'] + d_) for q in r['rows']]) for r in rs]
    isc = lambda r: (r['wmul'], r['rtol']) == c3
    for end, sign in ((MWE.EXP_A2P[1], +1), (MWE.EXP_A2P[0], -1)):
        for eps, keep in ((-1e-4, True), (1e-4, False)):
            d_ = (end + sign * eps) * lnr - m6 - sign * MWE.Z99 * e6
            P_ = [sh([r], d_)[0] if isc(r) else r for r in P]
            R_ = sh(more, d_) + mat + cap
            e_ = MWE.a2p_gate(P_, R_)[1][c3][1]
            assert abs((e_[2] if sign > 0 else e_[1]) - (end + sign * eps)) < 1e-9
            assert MWE.a2p_gate(P_, R_)[0] == keep and ("so the seed scales" in MWE.a2p_line(P_, None, R_)) == keep
    # the gate uses the further draws: moving A2R's further seeds alone flips it, A2P's 200 untouched
    d_ = 1.3 * lnr - m6
    assert not MWE.a2p_gate(P, sh(more, 3 * d_) + mat + cap)[0] and MWE.a2p_gate(P, None)[0]
    # R10-type: the geometric-mean band is open — a band whose upper end is a cell's geometric mean itself fails, just above passes
    g0 = max(g for s, _, _ in gt.values() for g in s['gm'])
    for top, keep in ((g0, False), (g0 * (1 + 1e-12), True)):
        monkeypatch.setattr(MWE, 'GM_A2P', (MWE.GM_A2P[0], top))
        assert MWE.a2p_gate(P, R)[0] == keep
    monkeypatch.undo()
    # the old band's margin on A2P's 200 at the 3A cell, from the data
    e2 = MWE.a2p_exp(MWE.a2p_stats([r for r in P if isc(r)]), lnr)
    assert (f"A2P's 200 alone passed at 3A, rtol 1e-10 by {1.1 - e2[2]:.4f} at its upper end, against the interval's half-width "
            f"{e2[2] - e2[0]:.4f}, marginal") in ml
    # R09-type: the largest offset named by |t| — a cell planted at −2.9 sem, beyond the 3A cell's +2.7, is named
    s5 = MWE.a2p_stats([y for y in P if (y['wmul'], y['rtol']) == (5.0, 1e-10)])
    Pn = [sh([r], lnr - 2.9 * s5['sem'] - s5['mean'])[0] if (r['wmul'], r['rtol']) == (5.0, 1e-10) else r for r in P]
    assert "no cell beyond 3 sem of it (the largest -2.9, 5A at rtol 1e-10: A2R, below)" in MWE.a2p_line(Pn, None, R)
    # D3ef-4: A2R's verdicts on the further draws alone, t and steps to two decimals, the Stouffer combination, the paired atol shift
    lr = MWE.a2r_line(R, P, H)
    sm, s1 = MWE.a2p_stats(more), MWE.a2p_stats(cap[0])
    assert f"over their {sm['n']} alone {sm['mean']:.3f} ± {sm['sem']:.3f} ({sm['t']:+.2f} sem)" in lr and f"({sm['t']:+.2f} sem)" == "(+2.97 sem)"
    assert f"on seed {cap[0]['seed']}'s first {cap[0]['k']} orders alone: Λ {s1['mean']:.3f} ± {s1['sem']:.3f} ({s1['t']:+.2f} sem)" in lr
    z = (sm['t'] + s1['t']) / math.sqrt(2)
    assert f"the 2 replications alone combined (Stouffer, Σt/√2) {z:+.2f} sem: the two arithmetics' runs' prefactors differ" in lr
    nd = lambda rs, k_: np.mean([x[k_] for r in rs for x in r['rows'][1:]])
    base = [r for r in P if isc(r)]
    assert (f"DOP853's double taking {nd(mat, 'nd'):.2f} steps against {nd(base + more, 'nd'):.2f} at its own atol (the long double "
            f"{nd(base + more, 'nld'):.2f})") in lr and f"{nd(base + more, 'nd'):.2f}" != f"{nd(more, 'nd'):.2f}"
    bd = {r['seed']: r for r in base + more}
    v = np.array([a['lag'] - b['lag'] for r in mat for a, b in zip(r['rows'][1:], bd[r['seed']]['rows'][1:])])
    assert v.size == 600 and f"on the same 600 orders (the long double's runs the same) {v.mean():+.3f} ± {v.std(ddof=1) / math.sqrt(600):.3f}: not the atol" in lr
    ev = v.std(ddof=1) / math.sqrt(600)
    for k_, txt in ((2.5, ": not the atol"), (3.5, ": the atol moves Λ")):   # the matched arm moved to k_ sem of the paired shift
        lk = MWE.a2r_line(more + sh(mat, k_ * ev - v.mean()) + cap, P, H)
        assert txt in lk and f"{k_ * ev:+.3f} ± {ev:.3f}" in lk
    sxm = MWE.a2p_stats([r for r in mat if r['seed'] in MWE.SEED_A2R])
    assert f"on the further seeds' 400 alone {sxm['mean']:.3f} ± {sxm['sem']:.3f} ({sxm['t']:+.2f} sem)" in lr
    assert "selected as A2P's cell farthest from 11 ln 2 (of 4, on its 200 draws: the selection's)" in lr
    assert "selected as A2H's farthest from 11 ln 2 (of 4, +4.65 sem on seed 131: the selection's)" in lr
    up3 = lambda rs, k: sh(rs, k)
    l5 = MWE.a2r_line(up3(more, 0.1) + mat + up3(cap, 0.1), P, H)      # the further draws beyond 3 sem: both verdicts flip
    assert ": on the further draws alone the two arithmetics' runs' prefactors differ at this cell" in l5
    assert ": replicated, the arithmetics' runs' prefactors differ at equal steps too" in l5
    l6 = MWE.a2r_line(up3(more, lnr - sm['mean']) + mat + up3(cap, lnr - s1['mean']), P, H)
    assert ": not resolved on the further draws alone" in l6
    # D3ef-5: A2H's DerSimonian–Laird, from the data; a homogeneous set gives τ = 0
    cells = [(w, q) for w in MWE.WMUL_A2P for q in MWE.RTOL_A2P]
    ss_ = [MWE.a2p_stats([r for r in P if (r['wmul'], r['rtol']) == c]) for c in cells] + [MWE.a2p_stats(r) for r in H]
    x8, e8 = np.array([s['mean'] / lnr for s in ss_]), np.array([s['sem'] / lnr for s in ss_])
    w8 = 1 / e8 ** 2
    mu = (w8 * x8).sum() / w8.sum()
    Q = (w8 * (x8 - mu) ** 2).sum()
    t2 = (Q - 7) / (w8.sum() - (w8 ** 2).sum() / w8.sum())
    wr = 1 / (e8 ** 2 + t2)
    lh = MWE.a2h_line(H, P)
    assert (f"a range of 8 estimates of sems {e8.min():.3f}–{e8.max():.3f}; their heterogeneity Q = {Q:.1f} on 7 dof, the "
            f"between-configuration sd τ = {math.sqrt(t2):.3f} (DerSimonian–Laird), the random-effects exponent "
            f"{(wr * x8).sum() / wr.sum():.3f} ± {wr.sum() ** -0.5:.3f} over configurations") in lh
    assert MWE.dl_meta([1.0, 1.0, 1.0], [0.1, 0.2, 0.3])[3] == 0.0 and MWE.dl_meta([1.0, 1.0], [0.1, 0.1])[5] == pytest.approx(0.1 / math.sqrt(2))
    # the md carries A2T's line after A2R's
    assert lt in md and md.index(lt) == md.index(lr) + 1 and md.index(lr) == md.index(lh) + 1
    rp = MWE.report(A2=raw['A2'], A2E=raw['A2E'], A2P=P, A2H=H, A2R=R, A2T=T).splitlines()
    assert rp.index(lt) == rp.index(lr) + 1 and ml in rp                # the report's A2P line gates on A2R too


def RT_alpha(L, wmul):
    from src.core import rich_transit as RT
    return RT._alpha_of(L, np.array([wmul * MWE.A_A]), RT.coef_crit)[0]


def A_f(xs, ss):
    return MWE.A_A * np.exp(-xs ** 2 / (2 * ss * ss))
