"""Proposition D‴.52: close features — the second entry's delay ∝ Δ⁻² (a dipole's ∝ Δ⁻¹, the transit of its derivative shape), the
two-kernel block's stiffnesses ∝ Δ² and Δ⁶ in closed form, and the fast flow's rates along those modes."""
import math
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_close_features as CF
import exp_mirror_transit as X
from src.core import rich_transit as RT


def test_block_stiffnesses_are_the_closed_forms():
    """At the exact fit the Gauss–Newton Hessian's two small eigenvalues are λ_split = Δ⁶⟨(P⊥e₃)²⟩/8 — the same for every pair of
    amplitudes — and λ_sep = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²), to O(Δ)."""
    L, K, c0 = CF.setup()
    D = 0.0125
    split = []
    for A1, A2 in ((0.8, 0.6), (0.8, -0.4), (0.8, 0.8)):
        r = CF.b_job((A1, A2, D))
        l1, l2 = r['pred']
        assert abs(r['hess'][0] / l1 - 1) < 5e-3
        assert abs(r['hess'][1] / l2 - 1) < 5e-2
        assert r['hess'][2] > 10 * r['hess'][1]                           # the other two are of order one
        split.append(r['hess'][0] / D ** 6)
    assert np.ptp(split) / np.mean(split) < 1e-3                         # amplitude-free
    vs = CF.taylor_basis(K, CF.CBAR)                                     # the continuum's √π/(48σ⁵ℓ) and √π/(8σ³ℓ), ℓ = Nh (round 55)
    s, ell = K.sig, K.xs.size * (K.xs[1] - K.xs[0])
    assert abs(CF.perp_norm2(vs, 3) / (np.sqrt(np.pi) / (48 * s ** 5 * ell)) - 1) < 1e-4
    assert abs(CF.perp_norm2(vs, 2) / (np.sqrt(np.pi) / (8 * s ** 3 * ell)) - 1) < 1e-4


def test_flow_relaxes_the_split_at_the_mobility_weighted_stiffness():
    """The fast flow's slowest mode is the split, at λ_split·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes' mobility); the separation relaxes at
    Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) while the amplitudes are frozen on its scale; both modes identified by their eigenvectors."""
    r = CF.b_job((0.8, 0.6, 0.0125))
    (rs, os_), (rp, op) = r['split'], r['sep']
    ps, pp = r['pred_flow']
    assert os_ > 0.999 and op > 0.999
    assert abs(rs / ps - 1) < 1e-2
    assert abs(rp / pp - 1) < 5e-2
    assert rs == min(r['rates'])
    L = X.law()
    mu_branch = 4 * L.G * L.c * np.array([0.8, 0.6]) / L.A                 # μ ≈ 4Gc|w|/a on the exponential branch
    assert np.allclose(r['mu'], mu_branch, rtol=1e-3)


def test_second_entry_delay_is_the_leading_order():
    """A same-sign pair Δ = 0.025 apart: the second entry comes after Δτ₂ = min_j (s_j − u_j)/v_j, v_j = ½V⟨ψφ_j⟩, to 0.5%; the path ends at
    the exact two-kernel fit."""
    r = CF.a_job((0.8, 0.6, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    dt, j, wall = r['pred']
    assert j2 == j and wall == 1
    assert abs((t2 - t1) / dt - 1) < 5e-3
    assert r['n_act'] == 2 and np.allclose(r['c_end'], r['cs'], atol=1e-6) and r['rms'] < 1e-10


def test_dipole_is_its_derivative_shape_slowed_by_the_separation():
    """A dipole (A₂ = −A₁) Δ = 0.025 apart: τ₁ and Δτ₂ are the derivative shape −A₁∂_cφ's, divided by Δ, to 0.5% and 1.5%, the same
    kernels entering."""
    r = CF.a_job((0.8, -0.8, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    sh = r['shape']
    assert [j1, j2] == sh['kernels']
    assert abs(t1 / sh['tau1'] - 1) < 5e-3
    assert abs((t2 - t1) / sh['dtau2'] - 1) < 1.5e-2


def test_predicted_delay_scales_as_the_inverse_square_separation():
    """The leading-order delay is ∝ Δ⁻²: halving Δ multiplies it by 4 (to O(Δ))."""
    L, K, c0 = CF.setup()
    d = []
    for D in (0.02, 0.01):
        f = X.bumps_target(K, CF.pair_centres(0.8, 0.6, D), (0.8, 0.6))
        d.append(CF.predicted_delay(K, f, c0, 0.8, 0.6, D)[0])
    assert abs(d[1] / d[0] / 4 - 1) < 2e-2



def test_c_job_stops_at_the_second_entry_and_places_within_a_budget(monkeypatch):
    """C's draws stop at the second entry (the delay only); the first R_D under H1 are relaxed again after it within RHS_D right-hand-side
    evaluations — placed (two kernels, their centres sorted) or out of budget, never unbounded (round 50's runs: a near-coalescent
    block's relaxation exhausted the memory)."""
    monkeypatch.setattr(CF, 'R_C', 3)
    monkeypatch.setattr(CF, 'R_D', 2)
    r = CF.c_job((0.3, 0.003, 1, 7))
    rows = r['rows']
    assert all(np.isfinite(q['dt']) for q in rows) and [('oob' in q) for q in rows] == [True, True, False]
    assert all('c' in q and q['c'][0] < q['c'][1] and not q['oob'] for q in rows[:2])
    assert all(abs(c_ - t) < 0.01 for q in rows[:2] for c_, t in zip(q['c'], CF.pair_centres(*CF.PAIR_C, 0.3)))
    monkeypatch.setattr(CF, 'RHS_D', 50)
    rows = CF.c_job((0.3, 0.003, 1, 7))['rows']
    assert all(q['oob'] and 'c' not in q for q in rows[:2]) and all(np.isfinite(q['dt']) for q in rows)
    h0 = CF.c_job((0.3, 0.003, 0, 7))['rows']
    assert not any('oob' in q for q in h0)


def test_the_opposite_sign_pair_enters_at_the_other_wall():
    """Referee round 55: for (0.8, −0.4) at Δ = 0.0125 the predicted second entrant is kernel 13 at the wall −1 (the flank, the
    minority's sign), the path's too, its delay 0.907 of the leading order (the O(ϑ) correction, ϑ = 25Δ)."""
    r = CF.a_job((0.8, -0.4, 0.0125))
    (t1, j1), (t2, j2) = r['entries'][:2]
    dt, j, wall = r['pred']
    assert (j, wall) == (13, -1) and j2 == 13 and 0.9 < (t2 - t1) / dt < 0.92
    L, K, c0 = CF.setup()                                                # with the T₃ term: within 1.3%
    f = X.bumps_target(K, CF.pair_centres(0.8, -0.4, 0.0125), (0.8, -0.4))
    d3 = CF.predicted_delay(K, f, c0, 0.8, -0.4, 0.0125, order=3)
    assert d3[1:] == (13, -1) and abs((t2 - t1) / d3[0] - 1) < 0.013


def test_the_split_moves_the_centres_by_half_the_separation_per_weight():
    """Referee round 55: the flow's split eigenvector moves each centre by Δ/(2w_k) per unit of weight moved, to 0.2% at Δ = 0.0125."""
    for A in ((0.8, 0.6), (0.8, -0.4)):
        r = CF.b_job(A + (0.0125,))
        assert r['split'][1] > 0.999 and np.allclose(r['split_c'], 1, atol=2e-3)


def test_the_separation_has_two_regimes():
    """Referee round 55: the separation's exact rate (−D·JᵀJ) is the frozen-amplitude formula's at Δ = 0.025 and tends to the
    relaxed-amplitude rate Gm·λ_sep — 2/3 of the frozen one — at Δ = 0.0005."""
    for A in ((0.8, 0.6), (0.8, -0.4)):
        fr, rl, m = CF.sep_regimes(*A, 0.025)
        assert abs(fr - 1) < 5e-3
        fr, rl, m = CF.sep_regimes(*A, 0.0005)
        assert abs(rl - 1) < 1e-2 and abs(fr - 2 / 3) < 1e-2


def test_detection_on_a_synthetic_table():
    """Referee round 55: `detection`'s H0 quantile, power and Δ₅₀ (log-linear interpolation of the power at ½) on synthetic draws."""
    C, pw = [], {}
    Ds = sorted(CF.DELTAS_C)
    for k, D in enumerate(Ds):
        for s in CF.SIGMAS_C:
            p = min(1.0, 0.1 * k * (1 + CF.SIGMAS_C.index(s) % 2))
            n1 = int(round(p * 20))
            C += [dict(D=D, s=s, hyp=0, rows=[dict(dt=float(t)) for t in range(1, 201)]),
                  dict(D=D, s=s, hyp=1, rows=[dict(dt=5.0)] * n1 + [dict(dt=500.0)] * (20 - n1))]
            pw[(D, s)] = n1 / 20
    det, d50 = CF.detection(C)
    q = np.quantile(np.arange(1, 201), CF.ALPHA_C)
    assert all(abs(det[k_]['q0'] - q) < 1e-12 and det[k_]['power'] == pw[k_] for k_ in pw)
    for s in CF.SIGMAS_C:
        p = [pw[(D, s)] for D in Ds]
        k = next(i for i in range(len(Ds) - 1) if p[i] < 0.5 <= p[i + 1])
        t = (0.5 - p[k]) / (p[k + 1] - p[k])
        assert abs(d50[s] - np.exp(np.log(Ds[k]) + t * np.log(Ds[k + 1] / Ds[k]))) < 1e-12


def test_the_knot_law_is_the_paths_first_order():
    """Referee round 55: under one bump with noise sξ, `knot_of` — the first entry by the duals, then their rates ⟨Rξ, φ_j⟩ to the walls —
    is the path's second knot to first order: the relative difference 10⁻⁸ at s = 10⁻³, falling ∝ s; its law's quantile ∝ 1/s up to the
    margins' O(s) noise."""
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, 0.025), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    xi = np.random.default_rng(3).standard_normal(K.xs.size)
    err = []
    for s in (1e-3, 1e-4):
        R = RT.sliding_iss(L, K, w1 * K.phi(np.array([c1]))[:, 0] + s * xi, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [t for t, j, k, c in R['events'] if k == 'enter']
        err.append(abs((ent[1] - ent[0]) / CF.knot_of(K, c0, c1, w1, s * xi[None, :])[0] - 1))
    assert err[0] < 1e-7 and 5 < err[0] / err[1] < 20
    q = [np.quantile(CF.knot_law(K, c0, c1, w1, s, 4000, np.random.default_rng(1)), 0.05) * s for s in (1e-5, 1e-4, 1e-3)]
    d1, d2 = q[1] / q[0] - 1, q[2] / q[0] - 1                    # ∝ 1/s up to the margins' noise, O(s)
    assert abs(d1) < 1e-3 and 5 < d2 / d1 < 15


def test_one_fit_is_stationary():
    """Referee round 55: `one_fit` returns the single-kernel least-squares fit (the gradient in (w, c) vanishes), not one Newton step."""
    L, K, c0 = CF.setup()
    for D in (0.3, 0.025):
        f = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D), CF.PAIR_C)
        c, w = CF.one_fit(K, f)
        ph, dp = K.phi(np.array([c]))[:, 0], K.dphi(np.array([c]))[:, 0]
        r = f - w * ph
        assert abs(ph @ r) < 1e-12 and abs(w * (dp @ r)) < 1e-12


def test_the_delay_draws_are_noisy_bounded_and_the_hypotheses_differ(monkeypatch):
    """Referee round 55: the delay runs stop at the second entry (max_entries = 2); the draws differ (noise); under H0 the delay is far
    longer than under H1 at Δ = 0.3; Δτ₂ is the gap τ₂ − τ₁."""
    monkeypatch.setattr(CF, 'R_C', 4)
    monkeypatch.setattr(CF, 'R_D', 0)
    kws, real = [], RT.sliding_iss

    def spy(*a, **kw):
        kws.append(kw)
        return real(*a, **kw)
    monkeypatch.setattr(RT, 'sliding_iss', spy)
    d1 = [q['dt'] for q in CF.c_job((0.3, 0.003, 1, 7))['rows']]
    d0 = [q['dt'] for q in CF.c_job((0.3, 0.003, 0, 7))['rows']]
    assert all(k.get('max_entries') == 2 for k in kws) and len(kws) == 8
    assert len(set(d1)) == 4 and len(set(d0)) == 4 and max(d1) < 10 < 100 < min(d0)
    fake = dict(events=[(5.0, 10, 'enter', 0.0), (7.5, 11, 'enter', 0.1)], hist=[], act=np.zeros(21, bool), relaxed=[True])
    monkeypatch.setattr(RT, 'sliding_iss', lambda *a, **kw: fake)
    assert CF.c_job((0.3, 0.003, 0, 7))['rows'][0]['dt'] == 2.5


def test_placed_blocks_are_sorted_and_relaxed(monkeypatch):
    """Referee round 55: a placement run's block is recorded with its coefficients in its centres' order, and only if it relaxed."""
    monkeypatch.setattr(CF, 'R_C', 1)
    monkeypatch.setattr(CF, 'R_D', 1)
    act = np.zeros(21, bool); act[[3, 7]] = True
    c = np.zeros(21); c[3], c[7] = 0.5, -0.5
    w = np.zeros(21); w[3], w[7] = 1.0, 2.0
    delay = dict(events=[(5.0, 10, 'enter', 0.0), (7.5, 11, 'enter', 0.1)], hist=[], act=act, relaxed=[True])
    for rel, placed in (([True, True], True), ([True, False], False)):
        place = dict(events=delay['events'], hist=[], act=act, c=c, w=w, relaxed=rel, out_of_budget=False)
        monkeypatch.setattr(RT, 'sliding_iss', lambda *a, **kw: place if kw.get('max_rhs') else delay)
        row = CF.c_job((0.3, 0.003, 1, 7))['rows'][0]
        assert ('c' in row) == placed and not row['oob']
        if placed:
            assert row['c'] == [-0.5, 0.5] and row['w'] == [2.0, 1.0]


def test_the_report_on_synthetic_results():
    """Referee round 55: the report's cells — A's ratio and scaled delays, B's slopes and the split's centre components, C's powers with
    the Neyman–Pearson bound, D's errors over the Cramér–Rao bounds, the blocks near their features and the adjacent second entrants."""
    A = [dict(A=(0.8, 0.6), D=0.025, entries=[(2.0, 10), (1002.0, 11)], pred=(800.0, 11, 1), n_act=2, cs=[0.0, 0.04],
              c_end=[0.0, 0.04], w_end=[0.8, 0.6], rms=1e-14)]
    md = CF.report(A=A)
    L, K, c0 = CF.setup()
    p3 = CF.predicted_delay(K, X.bumps_target(K, CF.pair_centres(0.8, 0.6, 0.025), (0.8, 0.6)), c0, 0.8, 0.6, 0.025, order=3)
    assert f'| 1000 | 800 (11, +1) | 1.2500 | {1000 / p3[0]:.4f} | 0.625 | 25 |' in md and p3[1:] == (11, 1)
    B = [dict(A=(0.8, 0.6), D=D, hess=[1e-6 * (D / 0.05) ** 6, 1e-2 * (D / 0.05) ** 2, 1.0, 2.0], rates=[1.0, 2.0, 3.0, 4.0],
              pred=(1e-6, 1e-2), split=(3e-6 * (D / 0.05) ** 6, 0.9999), sep=(5.0 * (D / 0.05) ** 2, 0.999), pred_flow=(1.0, 1.0),
              split_c=[1.0012, 1.0011], mu=[0.1, 0.1]) for D in (0.05, 0.025)]
    md = CF.report(B=B)
    assert md.count('| 6.000, 2.000 |') == 2 and '; centres 1.0012, 1.0011)' in md
    rows = lambda dts, **kw: [dict(dt=float(t), j=[10, 11], **kw) for t in dts]
    C = []
    for D in CF.DELTAS_C:
        for s in CF.SIGMAS_C:
            C.append(dict(D=D, s=s, hyp=0, rows=rows(range(100, 300))))
            cs = CF.pair_centres(*CF.PAIR_C, D)
            C.append(dict(D=D, s=s, hyp=1, rows=rows([50] * 200, c=list(cs + [0.001, -0.001]), w=[0.8, 0.6], oob=False)))
    md = CF.report(C=C)
    L, K, c0 = CF.setup()
    for D in (0.3, 0.0125):
        cs = CF.pair_centres(*CF.PAIR_C, D)
        f1 = X.bumps_target(K, cs, CF.PAIR_C)
        c1, w1 = CF.one_fit(K, f1)
        npw = CF.np_power(K, f1, w1 * K.phi(np.array([c1]))[:, 0], 0.001)
        assert f"| {D:g} | 1.00 [{npw:.2f}] (50, 200) |" in md
        cp, cf_ = CF.crb(K, cs, np.array(CF.PAIR_C), 0.001)
        assert (f"| {D:g} | 0.001 | 200, 0 (200) | 1.000e-03 | {1e-3 * D ** 2 / 1e-3:.3g} | {1e-3 / cp:.3g} | 0.000e+00 | 0 | 0 | "
                f"200 | 40 |") in md
    assert "Under H0 the second entrant is the first's neighbour in 100.0% of the 8000 draws, in 100.0% of the earliest 5%" in md
    assert 'The largest coefficient of a placed block: Δ = 0.3, s = 0.001, w = (0.8, 0.6) (their sum 1.4)' in md


def test_the_e_report_on_synthetic_draws(monkeypatch):
    """E's rows and summary (rounds 55, 59): the draws' H0 quantile against the law's, the threshold fixed at c̄ = CBAR applied out of
    sample (below CBAR) and beyond, the calibrated test's size with its Clopper–Pearson interval and the powers, on hand-made draws (H0
    delays uniform on [a, 2a], the bootstrap p-values i/100); the second entrants adjacent over all draws and the earliest 5% (round 64:
    the earliest draws' second entrant a neighbour); S and U passed on to E′."""
    rng = np.random.default_rng(3)
    D, s = CF.CELLS_E[0]
    E_ = []
    for cb, a in ((0.0, 1.0), (CF.CBAR, 1.0), (0.049, 0.1)):
        d0 = a * (1 + rng.random(100))
        d1 = a * (0.5 + rng.random(100))
        q0 = np.quantile(d0, CF.ALPHA_C)
        E_.append(dict(D=D, s=s, cbar=cb, rows0=[dict(dt=float(x), p=i / 100, j=[5, 6] if x < q0 else [5, 9]) for i, x in enumerate(d0)],
                       rows1=[dict(dt=float(x), p=0.0, j=[5, 7]) for x in d1], q_law=1.05 * float(np.quantile(d0, CF.ALPHA_C)),
                       ratio=np.full(21, 0.1), c1=cb, np=0.9))
    md = CF.report(E_=E_)
    lo, hi = CF.cp_int(5, 100)
    ref = E_[1]
    qref = np.quantile([q['dt'] for q in ref['rows0']], CF.ALPHA_C)
    fa0 = np.mean(np.array([q['dt'] for q in E_[0]['rows0']]) < qref)
    assert md.count(f'| {D:g}, {s:g} | ') == 3 and f'| 0.050, 0.050 [{lo:.3f}, {hi:.3f}] |' in md and f'| 1.000, 0.050 [{lo:.3f}, {hi:.3f}] |' in md
    assert (f"the false alarms of the threshold fixed at c̄ = {CF.CBAR:g} {fa0:.3f}–{fa0:.3f} at c̄ < {CF.CBAR:g} (out of sample) and "
            "1.000–1.000 beyond, of the calibrated test 0.050–0.050 (nominal 0.05)") in md
    pw = [np.mean(np.array([q['dt'] for q in e['rows1']]) < np.quantile([q['dt'] for q in e['rows0']], CF.ALPHA_C)) for e in E_]
    assert f"power minus the oracle's {1 - max(pw):+.3f} … {1 - min(pw):+.3f}, within 0.05 in 0 of 3" in md
    assert '| 0.100; 0.100, 0.100 |' in md                                       # the fit's kernel j₁, then j₁ ∓ 1
    assert md.count('| 0.05, 1.00 |') == 3                                       # adjacent: 5 of 100 draws, all of the earliest 5%
    seen = []
    monkeypatch.setattr(CF, 'e_prime_lines', lambda E2, S2, U2=None: (seen.append((E2, S2, U2)), ['E′ spy'])[1])
    assert 'E′ spy' in CF.report(E_=E_, S=[], U=['u']) and seen == [(E_, [], ['u'])]


def test_e_job_on_a_few_draws(monkeypatch):
    """Referee round 59: `e_job` draws H0 from the one-bump fit and H1 from the pair, one rng per cell; each draw's p is its knot against
    the law bootstrapped at its own first fit (`knot_law` at (ĉ, ŵ), the noise s, the draw's own rng); the law at the truth and the
    std ratio as defined; the cells' seeds distinct."""
    monkeypatch.setattr(CF, 'R0_E', 2)
    monkeypatch.setattr(CF, 'R1_E', 2)
    monkeypatch.setattr(CF, 'M_E', 64)
    monkeypatch.setattr(CF, 'M_LAW', 400)
    seen = []
    real = RT.sliding_iss
    monkeypatch.setattr(RT, 'sliding_iss', lambda L, K, y, *a, **k: (seen.append(y.copy()), real(L, K, y, *a, **k))[1])
    D, s, cb, seed = CF.e_jobs()[7]
    r = CF.e_job((D, s, cb, seed))
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    ys = [f0 + s * rng.standard_normal(K.xs.size) for _ in range(2)] + [f1 + s * rng.standard_normal(K.xs.size) for _ in range(2)]
    assert len(seen) == 4 and all(np.array_equal(a, b) for a, b in zip(seen, ys))
    for hyp, rows in ((0, r['rows0']), (1, r['rows1'])):
        for k, q in enumerate(rows):
            law = CF.knot_law(K, c0, q['c1'], q['w1'], s, 64, np.random.default_rng([seed, hyp, k]))
            assert q['p'] == np.mean(law <= q['dt'])
    assert r['q_law'] == np.quantile(CF.knot_law(K, c0, c1, w1, s, 400, np.random.default_rng([seed, 9])), CF.ALPHA_C)
    Rm, Phi = CF.tangent_projector(K, c1, w1), K.phi(c0)
    assert np.allclose(r['ratio'], np.sqrt(np.mean((Rm @ Phi) ** 2, axis=0) / np.mean(Phi ** 2, axis=0)))
    assert len({j[3] for j in CF.e_jobs()}) == len(CF.e_jobs()) and len({j[3] for j in CF.c_jobs()}) == len(CF.c_jobs())
    assert len({CF.f_seed(*j[:5]) for j in CF.f_jobs()}) == len(CF.f_jobs())


def test_crb_and_np_power_closed_forms():
    """Referee round 59: `crb` on an orthogonal design (JᵀJ diagonal) and `np_power` against Φ(‖f₁ − f₀‖/s − z₀.₉₅)."""
    from scipy.stats import norm

    class Ortho:
        xs = np.zeros(4)

        def phi(self, cs):
            return np.array([[2.0, 0.0], [0.0, 3.0], [0.0, 0.0], [0.0, 0.0]])

        def dphi(self, cs):
            return np.array([[0.0, 0.0], [0.0, 0.0], [5.0, 0.0], [0.0, 7.0]])
    w, s = np.array([0.8, 0.6]), 0.3
    pos, split = CF.crb(Ortho(), None, w, s)
    W = w.sum()
    assert np.isclose(pos, s * np.sqrt((1 / (25 * w[0] ** 2) + 1 / (49 * w[1] ** 2)) / 2), rtol=1e-12)
    assert np.isclose(split, s * np.sqrt(w[1] ** 2 / W ** 4 / 4 + w[0] ** 2 / W ** 4 / 9), rtol=1e-12)
    L, K, c0 = CF.setup()
    d = np.zeros(K.xs.size)
    d[3] = 1.0
    for a in (0.0, 1.0, 2.5):
        f1 = (norm.ppf(0.95) + a) * 0.01 * d
        assert np.isclose(CF.np_power(K, f1, 0 * d, 0.01), norm.cdf(a), rtol=1e-12)


def test_d50_bootstrap_on_a_deterministic_table():
    """Referee round 59: a power that steps from 0 to 1 between two grid separations crosses ½ at their geometric mean, in every
    resample; `crossing` interpolates in log Δ."""
    Ds = sorted(CF.DELTAS_C)
    C = []
    for s in CF.SIGMAS_C:
        step = Ds[3] if s < 0.01 else Ds[5]
        for D in Ds:
            C.append(dict(D=D, s=s, hyp=0, rows=[dict(dt=float(v)) for v in np.linspace(10.0, 20.0, 20)]))
            C.append(dict(D=D, s=s, hyp=1, rows=[dict(dt=1.0 if D >= step else 50.0) for _ in range(20)]))
    ci, sl = CF.d50_bootstrap(C, n=20)
    for s in CF.SIGMAS_C:
        k = 3 if s < 0.01 else 5
        g = np.sqrt(Ds[k - 1] * Ds[k])
        assert np.isclose(ci[s][0], g) and np.isclose(ci[s][1], g) and ci[s][2] == 0
    assert np.isclose(CF.crossing([1.0, 2.0, 4.0], [0.2, 0.4, 0.6]), 2 * np.sqrt(2))
    assert np.isnan(CF.crossing([1.0, 2.0], [0.1, 0.3]))


def test_the_knot_law_on_several_noise_vectors():
    """Referee round 59: `knot_of` is the path's second knot to first order over several noise vectors — the error ∝ s."""
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, 0.025), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    y0 = w1 * K.phi(np.array([c1]))[:, 0]
    err = {}
    for s in (1e-3, 1e-4):
        e = []
        for k in range(6):
            xi = np.random.default_rng(100 + k).standard_normal(K.xs.size)
            R = RT.sliding_iss(L, K, y0 + s * xi, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
            ent = [t for t, j, kk, c in R['events'] if kk == 'enter']
            e.append(abs((ent[1] - ent[0]) / CF.knot_of(K, c0, c1, w1, s * xi[None, :])[0] - 1))
        err[s] = np.array(e)
    assert err[1e-3].max() < 1e-2 and 4 < np.median(err[1e-3]) / np.median(err[1e-4]) < 25


def test_the_score_statistic():
    """Referee round 59: z = ⟨r, R∂²_cφ⟩/(s‖R∂²_cφ‖) at a one-bump fit — zero on a one-bump target; on a close pair's one-bump fit the
    signal ‖f₁ − f₀‖ (r ∥ R∂²_cφ) with the pair's sign sign(A₁A₂); N(0, 1) under one bump with noise; the width's z_σ equal to it to first
    order; σ̂² − σ² → V/A on a close pair (the pair's variance added)."""
    L, K, c0 = CF.setup()
    y = 0.9 * K.phi(np.array([0.03]))[:, 0]
    assert abs(CF.score_stats(K, y, 0.03, 0.9, 0.01)[0]) < 1e-9
    for A in ((0.8, 0.6), (0.8, -0.4)):
        D = 0.01                                                       # ϑ = 0.12, 0.25: r ∥ R∂²φ to O(ϑ²)
        f1 = X.bumps_target(K, CF.pair_centres(*A, D), A)
        c1, w1 = CF.one_fit(K, f1)
        z = CF.score_stats(K, f1, c1, w1, 1.0)[0]
        assert np.sign(z) == np.sign(A[0] * A[1]) and abs(z) / np.linalg.norm(f1 - w1 * K.phi(np.array([c1]))[:, 0]) > 0.99
        sg, sd = CF.width_fit(K, f1, c1, w1)
        V_A = A[0] * A[1] * D ** 2 / sum(A) ** 2
        assert abs((sg ** 2 - K.sig ** 2) / V_A - 1) < 0.05
    rng = np.random.default_rng(8)
    f0 = 1.4 * K.phi(np.array([0.032]))[:, 0]
    zs, ws = [], []
    for _ in range(600):
        yy = f0 + 0.01 * rng.standard_normal(K.xs.size)
        c, w = CF.one_fit(K, yy)
        zs.append(CF.score_stats(K, yy, c, w, 0.01)[0])
        ws.append(CF.width_z(K, yy, c, w, 0.01))
    zs, ws = np.array(zs), np.array(ws)
    assert abs(zs.mean()) < 0.15 and 0.9 < zs.std() < 1.1
    assert np.abs(ws - zs).max() < 0.05


def test_the_knot_per_wall():
    """Referee round 59: `knot_walls`' minimum is `knot_of`; `walls_at_fit` at the path's own first fit is the path's second knot."""
    L, K, c0 = CF.setup()
    Xi = 0.01 * np.random.default_rng(1).standard_normal((300, K.xs.size))
    a, b = CF.knot_walls(K, c0, 0.032, 1.4, Xi)
    assert np.array_equal(np.minimum(a, b), CF.knot_of(K, c0, 0.032, 1.4, Xi))
    f1 = X.bumps_target(K, CF.pair_centres(0.8, -0.4, 0.025), (0.8, -0.4))
    for k in range(3):
        y = f1 + 0.01 * np.random.default_rng(20 + k).standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [t for t, j, kk, c in R['events'] if kk == 'enter']
        h = R['hist'][1]
        i = int(np.flatnonzero(h[2])[0])
        tA, tB = CF.walls_at_fit(K, c0, y, float(h[1][i]), float(h[2][i]))
        assert abs(min(tA, tB) / (ent[1] - ent[0]) - 1) < 1e-9


def test_f_and_s_jobs_on_a_few_draws(monkeypatch):
    """Referee round 59: `f_job`'s rows — the knot per wall at the path's first fit, the p-values against `knot_walls` at that fit with the
    draw's rng, the law at the truth with its own; `s_job` regenerates C's labels in the job's order (its scores are those of `c_job`'s
    labels at their one-bump fits)."""
    monkeypatch.setattr(CF, 'M_E', 64)
    monkeypatch.setattr(CF, 'M_LAW', 400)
    A, cb, s, D = (0.8, -0.4), 0.032, 0.01, 0.025
    r = CF.f_job((A, cb, s, D, 1, 3))
    L, K, c0 = CF.setup()
    seed = CF.f_seed(A, cb, s, D, 1)
    N = K.xs.size
    for k, q in enumerate(r['rows']):
        assert abs(min(q['tA'], q['tB']) / q['dt'] - 1) < 1e-9
        la, lb = CF.knot_walls(K, c0, q['c1'], q['w1'], s * np.random.default_rng([seed, 1, k]).standard_normal((64, N)))
        assert q['p'] == np.mean(np.minimum(la, lb) <= q['dt']) and q['pA'] == np.mean(la <= q['tA']) and q['pB'] == np.mean(lb <= q['tB'])
    f1 = X.bumps_target(K, CF.pair_centres(*A, D, cb), A)
    c1, w1 = CF.one_fit(K, f1)
    la, lb = CF.knot_walls(K, c0, c1, w1, s * np.random.default_rng([seed, 0]).standard_normal((400, N)))
    assert r['q_law'] == np.quantile(np.minimum(la, lb), CF.ALPHA_C) and r['qB'] == np.quantile(lb, CF.ALPHA_C / 2)
    assert r['qA'] == np.quantile(la, CF.ALPHA_C / 2)                     # referee round 67 (F03): each wall at α/2
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    assert r['np'] == CF.np_power(K, f1, f0, s, CF.ALPHA_C)               # (F04): the bound at the draws' noise s
    seen = []
    monkeypatch.setattr(RT, 'sliding_iss', lambda L_, K_, y, *a, **kw: (seen.append(y.copy()), dict(events=[]))[1])
    Dc, sc, hc, sdc = CF.c_jobs()[3]
    monkeypatch.setattr(CF, 'R_C', 3)
    CF.c_job((Dc, sc, hc, sdc))
    out = CF.s_job(('C', Dc, sc, CF.CBAR, hc, sdc, (3,), None))
    for y, q in zip(seen, out[f'rows{hc}']):
        c, w = CF.one_fit(K, y)
        assert q['z'] == CF.score_stats(K, y, c, w, sc)[0]


def test_the_f_report_on_synthetic_jobs():
    """Referee round 59: F's report — Δ₅₀ per test from the rejections, '<'/'>' at the grid's ends, the false alarms at H0."""
    s, cb = CF.SIGMAS_F[1], CF.CBAR
    Ds = CF.DELTAS_F[s]
    F = []
    for A in CF.PAIRS_F:
        sg = np.sign(A[0] * A[1])
        for D in Ds:
            k = Ds.index(D)
            rows = [dict(dt=1.0 if i < 50 * k else 9.0, p=0.5, pA=0.5, pB=0.5, tA=9.0, tB=9.0, z=sg * (3.0 if i < 40 * k else 0.0),
                         zs=sg * (3.0 if i < 40 * k else 0.0), w1=1.0) for i in range(200)]
            F.append(dict(A=A, cbar=cb, s=s, D=D, hyp=1, rows=rows, q_law=2.0, qA=1.0, qB=1.0, np=min(1.0, 0.2 * (k + 1)), c1=cb,
                          w1=1.0))
        rows0 = [dict(dt=1.0 if i < 20 else 9.0, p=0.01 if i < 10 else 0.5, pA=0.5, pB=0.5, tA=1.0 if i < 10 else 9.0, tB=9.0, z=0.0,
                      zs=0.0, w1=1.0) for i in range(400)]                # the knot's H0 rejections all at the fit's-sign wall
        F.append(dict(A=A, cbar=cb, s=s, D=Ds[CF.H0_F], hyp=0, rows=rows0, q_law=2.0, qA=1.0, qB=1.0, np=0.0, w1=1.0))
    md = '\n'.join(CF.f_lines(F))
    knot = CF.crossing(list(Ds), [min(1.0, 0.25 * j) for j in range(len(Ds))])
    z = CF.crossing(list(Ds), [min(1.0, 0.2 * j) for j in range(len(Ds))])
    npd = CF.crossing(list(Ds), [min(1.0, 0.2 * (j + 1)) for j in range(len(Ds))])
    row = lambda al: (f"| {s:g} | {cb:g} | {knot:.4g}, > {Ds[-1]:g}, > {Ds[-1]:g}, > {Ds[-1]:g}, > {Ds[-1]:g}, {z:.4g}, {z:.4g}, {z:.4g}, "
                      f"{al}, {npd:.4g} | 0.050, 0.025, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |")
    assert md.count(row(f"{z:.4g}")) == 1 and md.count(row(f"> {Ds[-1]:g}")) == 1    # the allocation 1: one-sided toward the fit's sign
    assert 'spends α·1.000 toward the fit\'s sign' in md and 'The noise over the one-bump fit\'s coefficient, s/|ŵ|' in md
    assert '## F. Resolution across the grid cell, for both signs' in md


def test_the_walls_take_the_pairs_sign():
    """Referee round 64: at the path's own first fit the knot falls on the wall of the fit's sign for a same-sign pair (tA = Δτ₂ < tB) and
    on the other wall for an opposite-sign pair (tB = Δτ₂ < tA)."""
    L, K, c0 = CF.setup()
    for A, own in (((0.8, 0.6), True), ((0.8, -0.4), False)):
        f1 = X.bumps_target(K, CF.pair_centres(*A, 0.025), A)
        y = f1 + 0.001 * np.random.default_rng(5).standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [t for t, j, kk, c in R['events'] if kk == 'enter']
        h = R['hist'][1]
        i = int(np.flatnonzero(h[2])[0])
        tA, tB = CF.walls_at_fit(K, c0, y, float(h[1][i]), float(h[2][i]))
        assert (tA < tB) == own and abs(min(tA, tB) / (ent[1] - ent[0]) - 1) < 1e-9


def test_the_score_against_an_independent_formula_and_a_negative_pair():
    """Referee round 64: z = ⟨r, ψ̂⟩/(s‖ψ̂‖) with ψ̂ the projection of ∂²_cφ off span(φ, ŵ∂_cφ), by least squares, to 10⁻¹⁰; for the pair
    (−0.8, −0.6) z < 0 and ŵz > 0 (ŵz has the sign of A₁A₂), and the one-sided rejection uses ŵ's sign; z_σ ≈ sign(ŵ)z at small noise."""
    L, K, c0 = CF.setup()
    rng = np.random.default_rng(11)
    for A in ((0.8, 0.6), (-0.8, -0.6), (0.8, -0.4)):
        f1 = X.bumps_target(K, CF.pair_centres(*A, 0.03), A)
        y = f1 + 0.002 * rng.standard_normal(K.xs.size)
        c, w = CF.one_fit(K, y)
        T = np.c_[K.phi(np.array([c]))[:, 0], w * K.dphi(np.array([c]))[:, 0]]
        d2 = CF.d2phi(K, c)[:, 0]
        psi = d2 - T @ np.linalg.lstsq(T, d2, rcond=None)[0]
        r = y - w * K.phi(np.array([c]))[:, 0]
        z = CF.score_stats(K, y, c, w, 0.002)[0]
        assert abs(z - (r @ psi) / (0.002 * np.linalg.norm(psi))) < 1e-10 * max(1.0, abs(z))
        if A[0] < 0:
            assert w < 0 and z < 0 and w * z > 0
            rows = [dict(z=z, zs=0.0, w1=w)]
            assert CF.rejections(rows, 1.0)['z'][0] == (abs(z) > CF.Z95)
    f0 = -1.4 * K.phi(np.array([0.032]))[:, 0]
    zs, ws = [], []
    for _ in range(200):
        yy = f0 + 0.01 * rng.standard_normal(K.xs.size)
        c, w = CF.one_fit(K, yy)
        zs.append(CF.score_stats(K, yy, c, w, 0.01)[0])
        ws.append(CF.width_z(K, yy, c, w, 0.01))
    assert np.abs(np.array(ws) + np.array(zs)).max() < 0.05                        # ŵ < 0: z_σ ≈ −z


def test_f_job_under_h0_draws_the_null_with_its_own_seed(monkeypatch):
    """Referee round 64: `f_job` with hyp = 0 draws the labels from the pair's one-bump fit with the H0 seed: its first draw's path fit is
    the one-bump fit of f₀ plus that seed's first noise vector."""
    monkeypatch.setattr(CF, 'M_E', 16)
    monkeypatch.setattr(CF, 'M_LAW', 64)
    A, cb, s, D = (0.8, -0.4), 0.032, 0.01, 0.025
    assert CF.f_seed(A, cb, s, D, 0) != CF.f_seed(A, cb, s, D, 1)
    r = CF.f_job((A, cb, s, D, 0, 1))
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*A, D, cb), A)
    c1, w1 = CF.one_fit(K, f1)
    y = w1 * K.phi(np.array([c1]))[:, 0] + s * np.random.default_rng(CF.f_seed(A, cb, s, D, 0)).standard_normal(K.xs.size)
    c, w = CF.one_fit(K, y)
    assert abs(r['rows'][0]['c1'] - c) < 1e-7 and abs(r['rows'][0]['w1'] - w) < 1e-7 and r['hyp'] == 0


def test_part_s_and_part_u_take_e_and_cs_seeds_and_fits(monkeypatch):
    """Referee round 64: `part_S` and `part_U` hand each cell its own job's seed (`c_jobs`, `e_jobs`) and E's fits in the draws' order (H0's
    rows, then H1's)."""
    import multiprocessing
    seen = []

    class FakePool:
        def __init__(self, n):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            seen.append((fn.__name__, list(jobs)))
            return [None] * len(jobs)
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    Dc, sc, hc, sdc = CF.c_jobs()[5]
    De, se, cbe, sde = CF.e_jobs()[7]
    C = [dict(D=Dc, s=sc, hyp=hc, rows=[{}] * 3)]
    E_ = [dict(D=De, s=se, cbar=cbe, rows0=[dict(c1=0.1, w1=1.0), dict(c1=0.2, w1=1.1)], rows1=[dict(c1=0.3, w1=1.2)])]
    CF.part_S(C, E_)
    CF.part_U(E_)
    (n1, j1), (n2, j2) = seen
    assert n1 == 's_job' and j1[0] == ('C', Dc, sc, CF.CBAR, hc, sdc, (3,), None)
    assert j1[1] == ('E', De, se, cbe, 2, sde, (2, 1), [(0.1, 1.0), (0.2, 1.1), (0.3, 1.2)])
    assert n2 == 'u_job' and j2[0] == (De, se, cbe, sde, (2, 1), [(0.1, 1.0), (0.2, 1.1), (0.3, 1.2)])


def test_e_job_calibrates_at_the_paths_fit(monkeypatch):
    """Referee round 64 (D3bf's M08): `e_job` bootstraps each draw's knot law at the path's own first fit, not at the truth."""
    monkeypatch.setattr(CF, 'R0_E', 1)
    monkeypatch.setattr(CF, 'R1_E', 1)
    monkeypatch.setattr(CF, 'M_LAW', 32)
    calls = []
    real = CF.knot_law

    def spy(K, c0, c, w, s, M, rng):
        calls.append((c, w, M))
        return real(K, c0, c, w, s, M, rng)
    monkeypatch.setattr(CF, 'knot_law', spy)
    D, s, cb, seed = CF.e_jobs()[0]
    r = CF.e_job((D, s, cb, seed))
    per_draw = [(c, w) for c, w, M in calls if M == CF.M_E]
    assert per_draw == [(r['rows0'][0]['c1'], r['rows0'][0]['w1']), (r['rows1'][0]['c1'], r['rows1'][0]['w1'])]
    assert [(c, w) for c, w, M in calls if M == 32] == [(r['c1'], r['w1'])]


def test_the_runner_ups_own_knot(monkeypatch):
    """Referee round 64: `knot_of_kernel` is the tie pair's column of `knot_of`'s duals (the runner-up's, or the noiseless first entrant's
    where the runner-up enters first); `u_job` scores E's regenerated draws with it from the stored fits."""
    L, K, c0 = CF.setup()
    N = K.xs.size
    Xi = 0.01 * np.random.default_rng(3).standard_normal((200, N))
    c, w = 0.049, 1.4
    Phi = K.phi(c0)
    v0 = (w * K.phi(np.array([c]))[:, 0][None, :] + Xi) @ Phi / N
    j1 = int(np.argmax(np.abs(w * K.phi(np.array([c]))[:, 0] @ Phi)))
    j2 = j1 + 1 if abs(c0[j1 + 1] - c) < abs(c0[j1 - 1] - c) else j1 - 1
    t = CF.knot_of_kernel(K, c0, c, w, Xi, j2, j1)
    n = Xi @ (CF.tangent_projector(K, c, w) @ Phi) / N
    for m in range(200):
        k = j1 if int(np.argmax(np.abs(v0[m]))) == j2 else j2
        u = v0[m, k] / np.abs(v0[m]).max()
        assert np.isclose(t[m], (1 - np.sign(n[m, k]) * u) / abs(n[m, k]), rtol=1e-12)
    assert np.all(t >= CF.knot_of(K, c0, c, w, Xi) - 1e-12)                      # the dictionary's minimum is at most the pair's
    monkeypatch.setattr(CF, 'M_LAW', 64)
    D, s, cb, seed = CF.e_jobs()[2]
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    ys = [f0 + s * rng.standard_normal(N), f1 + s * rng.standard_normal(N)]
    fits = [CF.one_fit(K, y) for y in ys]
    u = CF.u_job((D, s, cb, seed, (1, 1), fits))
    jj2 = CF.predicted_delay(K, f1, c0, *CF.PAIR_C, D, cbar=cb)[1]
    jj1 = int(np.argmax(np.abs(Phi.T @ f1)))
    for h, (y, (cf, wf)) in enumerate(zip(ys, fits)):
        g = Phi.T @ y / N
        k = jj1 if int(np.argmax(np.abs(g))) == jj2 else jj2
        v = Phi[:, k] @ (y - wf * K.phi(np.array([cf]))[:, 0]) / N
        assert np.isclose(u[f't2_{h}'][0], (1 - np.sign(v) * g[k] / np.abs(g).max()) / abs(v), rtol=1e-12)
    Xi = s * np.random.default_rng([seed, 9]).standard_normal((64, N))            # E's law samples (referee round 71)
    law = CF.knot_of_kernel(K, c0, c1, w1, Xi, jj2, jj1)
    assert u['q2'] == float(np.quantile(law, CF.ALPHA_C)) and u['j2'] == jj2
    assert u['q_full'] == float(np.quantile(CF.knot_of(K, c0, c1, w1, Xi), CF.ALPHA_C))
    assert u['q_full'] == float(np.quantile(CF.knot_law(K, c0, c1, w1, s, 64, np.random.default_rng([seed, 9])), CF.ALPHA_C))


def test_the_e_prime_report_with_the_runner_ups_knot():
    """Referee rounds 64 and 67: E′'s columns — the runner-up's own knot's size and power against its law's quantile, the knot against
    its own law at the truth, the t-version's size (Z02), the free width's power (D02), the margin from the uncensored μ̂ — and the
    minimum's part (the knot against its law less the own knot) per cell, with U's first entrants against S's swaps."""
    D, s_ = CF.CELLS_E[0]
    cb = CF.OFFSETS_E[0]
    rows0 = [dict(dt=1.0 if i < 2 else 3.0, p=0.5, c1=0.0, w1=1.0) for i in range(10)]
    rows1 = [dict(dt=1.0, p=0.5, c1=0.0, w1=1.0) for _ in range(4)]
    E_ = [dict(D=D, s=s_, cbar=cb, rows0=rows0, rows1=rows1, np=0.5, q_law=2.0)]
    row = lambda t=0.0, zs=0.0, mu_u=0.1: dict(z=0.0, t=t, zs=zs, z2=0.0, swap=False, mu=0.1, mu_u=mu_u, dfit=0.0, w1=1.0)
    S = [dict(kind='E', D=D, s=s_, cbar=cb, rows0=[row(t=2.0 if i < 3 else 0.0, mu_u=0.1 + 0.01 * (i % 2)) for i in range(10)],
              rows1=[row(zs=2.0 if i < 2 else -2.0) for i in range(4)], mu_true=0.1, cos2=0.9)]
    U = [dict(D=D, s=s_, cbar=cb, q2=2.0, t2_0=[1.0] + [3.0] * 9, t2_1=[1.0, 1.0, 1.0, 3.0], first_0=[False] * 10, first_1=[False] * 4,
              mu0=0.1)]
    md = '\n'.join(CF.e_prime_lines(E_, S, U))
    assert '| 0.100, 0.750 |' in md and '; 0.200, 1.000 |' in md            # the own knot; the knot against its law at the truth
    assert '| 0.300 |' in md and ', 0.500 |' in md                          # the t size (3 of 10), z_σ's power (2 of 4)
    assert f"0.1000, {0.005 / s_:.3f}, {0.1 / 0.005:.2f} |" in md               # sd(μ̂) of the uncensored margins 0.10/0.11
    assert f"({cb:g}, {D:g}, {s_:g}) +0.2500" in md and '0.0000 at most apart' in md
    assert '| — |' in '\n'.join(CF.e_prime_lines(E_, S))
    U[0]['mu0'] = 0.2                                                      # referee round 75 (EM02): the null's margin, not mu_true
    assert f"0.2000, {0.005 / s_:.3f}, {0.2 / 0.005:.2f} |" in '\n'.join(CF.e_prime_lines(E_, S, U))
    U[0]['mu0'] = 0.1
    U[0].pop('mu0')                                                        # referee round 73: no null margin without U's mu0
    assert f"| —, {0.005 / s_:.3f}, — |" in '\n'.join(CF.e_prime_lines(E_, S, U))


def test_the_walls_per_draw_and_per_law_against_an_independent_definition():
    """Referee round 67 (D01, K01): `draw_stats`' tA, tB are `walls_at_fit`'s in order — the fit's-sign wall first, the knot of a
    same-sign pair there —, and `knot_walls`' per-wall minima equal an independent computation: each kernel but the first entrant moves
    at n_j = ⟨Rξ, φ_j⟩ toward the wall sign(n_j), reaching it at (1 − sign(n_j)u_j)/|n_j|, on the fit's-sign wall iff sign n_j = sign ŵ."""
    L, K, c0 = CF.setup()
    N = K.xs.size
    f1 = X.bumps_target(K, CF.pair_centres(0.8, 0.6, 0.025), (0.8, 0.6))
    y = f1 + 0.003 * np.random.default_rng(8).standard_normal(N)
    c, w = CF.one_fit(K, y)
    d = CF.draw_stats(K, c0, y, c, w, 0.003)
    assert (d['tA'], d['tB']) == CF.walls_at_fit(K, c0, y, c, w) and d['tA'] < d['tB']
    for c_, w_ in ((0.032, 1.4), (0.032, -1.4)):
        Xi = 0.01 * np.random.default_rng(4).standard_normal((50, N))
        la, lb = CF.knot_walls(K, c0, c_, w_, Xi)
        Phi = K.phi(c0)
        R = CF.tangent_projector(K, c_, w_)
        for m in range(50):
            v = (w_ * K.phi(np.array([c_]))[:, 0] + Xi[m]) @ Phi / N
            j1 = int(np.argmax(np.abs(v)))
            n = (R @ Xi[m]) @ Phi / N
            own, oth = np.inf, np.inf
            for j in range(c0.size):
                if j == j1 or n[j] == 0:
                    continue
                t = (1 - np.sign(n[j]) * v[j] / abs(v[j1])) / abs(n[j])
                if np.sign(n[j]) == np.sign(w_):
                    own = min(own, t)
                else:
                    oth = min(oth, t)
            assert np.isclose(la[m], own, rtol=1e-10) and np.isclose(lb[m], oth, rtol=1e-10)


def test_f_jobs_and_the_width_escape():
    """Referee round 67 (F07, S03): F's H0 draws are at the third separation of each (pair, offset, s), the H1 cells at all; `width_z`
    is finite while the fitted width stays in (σ/5, 5σ) and ±∞ beyond."""
    jobs = CF.f_jobs()
    h0 = [j for j in jobs if j[4] == 0]
    assert CF.H0_F == 2 and len(h0) == len(CF.PAIRS_F) * len(CF.OFFSETS_F) * len(CF.SIGMAS_F)
    assert all(j[3] == CF.DELTAS_F[j[2]][2] and j[5] == CF.R0_F for j in h0)
    L, K, c0 = CF.setup()
    s0 = K.sig
    for sg, fin in ((0.21 * s0, True), (0.19 * s0, False), (4.9 * s0, True), (5.1 * s0, False)):
        old = CF.width_fit
        CF.width_fit = lambda K_, y, c, w, sg=sg: (sg, 1.0)
        try:
            z = CF.width_z(K, np.zeros(K.xs.size), 0.0, 1.0, 0.01)
        finally:
            CF.width_fit = old
        assert np.isfinite(z) == fin and (fin or np.sign(z) == np.sign(sg - s0))


def test_u_job_at_a_swapping_tie_and_at_the_centre(monkeypatch):
    """Referee round 67 (U09, U01): at the tie c̄ = 0.049, where the runner-up j₂ often leads, `u_job`'s dual is over the largest
    correlation (j₂'s own when it leads: the first entrant's knot then); at c̄ = 0 the runner-up is j₁ − 1, as `predicted_delay` says."""
    monkeypatch.setattr(CF, 'M_LAW', 64)
    L, K, c0 = CF.setup()
    N = K.xs.size
    Phi = K.phi(c0)
    for D, s, cb, seed in (CF.e_jobs()[14], CF.e_jobs()[0]):
        f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
        c1, w1 = CF.one_fit(K, f1)
        f0 = w1 * K.phi(np.array([c1]))[:, 0]
        rng = np.random.default_rng(seed)
        ys = [f0 + s * rng.standard_normal(N) for _ in range(12)]
        fits = [CF.one_fit(K, y) for y in ys]
        u = CF.u_job((D, s, cb, seed, (12, 0), fits))
        j2 = CF.predicted_delay(K, f1, c0, *CF.PAIR_C, D, cbar=cb)[1]
        j1 = int(np.argmax(np.abs(Phi.T @ f1)))
        assert u['j2'] == j2 and (j2 == j1 - 1) == (cb == 0.0)
        swaps = 0
        for y, (cf, wf), t2, fst in zip(ys, fits, u['t2_0'], u['first_0']):
            g = Phi.T @ y / N
            lead = int(np.argmax(np.abs(g))) == j2
            k = j1 if lead else j2
            v = Phi[:, k] @ (y - wf * K.phi(np.array([cf]))[:, 0]) / N
            assert np.isclose(t2, (1 - np.sign(v) * g[k] / np.abs(g).max()) / abs(v), rtol=1e-12) and fst == lead
            swaps += lead
        assert (swaps > 0) == (cb == 0.049)


def test_the_allocated_score_and_the_paired_comparison():
    """Referee round 67: the score with the knot's allocation spends α·a toward the fit's sign and α·(1 − a) away (a = 1: z one-sided);
    `knot_allocation` is the share of the calibrated knot's H0 rejections at the fit's-sign wall; `paired` compares two tests' Δ₅₀ over
    their common cells."""
    from scipy.stats import norm
    rows = [dict(z=z, zs=0.0, w1=w) for z, w in ((1.70, 1.0), (-1.70, -1.0), (-2.0, 1.0), (-2.7, 1.0), (1.60, 1.0))]
    r = CF.rejections(rows, 1.0, alloc=0.907)
    assert norm.ppf(1 - 0.05 * 0.907) < 1.70 and list(r['alloc']) == [True, True, False, True, False]
    assert list(CF.rejections(rows, 1.0, alloc=1.0)['alloc']) == list(CF.rejections(rows, 1.0)['z'])
    q = lambda p, tA, tB: dict(p=p, tA=tA, tB=tB)
    F = [dict(A=(0.8, 0.6), hyp=0, rows=[q(0.01, 1, 2), q(0.01, 1, 2), q(0.01, 3, 2), q(0.5, 1, 2)]),
         dict(A=(0.8, 0.6), hyp=1, rows=[q(0.01, 3, 2)])]
    assert CF.knot_allocation(F, (0.8, 0.6)) == 2 / 3 and CF.knot_allocation(F, (0.8, -0.4)) == 1.0
    txt = CF.paired({1: 1.1, 2: 0.9, 3: 1.2}, {1: 1.0, 2: 1.0, 4: 1.0}, 'a against b')
    assert txt == 'a against b 1.000 [0.900, 1.100], smaller in 1 of 2' and CF.paired({}, {}, 'x') == 'x —'


def test_the_separations_bound():
    """Referee round 64: Smith's criterion — the separation's bound ∝ s/Δ² with the amplitudes unknown and unequal, ∝ s/Δ for equal or
    known amplitudes."""
    L, K, c0 = CF.setup()
    sl = lambda A, known: np.log(CF.crb_sep(K, CF.pair_centres(*A, 0.00625), np.array(A), 1.0, known)
                                 / CF.crb_sep(K, CF.pair_centres(*A, 0.0125), np.array(A), 1.0, known)) / np.log(2)
    assert abs(sl((0.8, 0.6), False) - 2) < 0.05 and abs(sl((0.8, -0.4), False) - 2) < 0.05
    assert abs(sl((0.8, 0.8), False) - 1) < 0.05 and abs(sl((0.8, 0.6), True) - 1) < 0.05
    assert "Smith's criterion" in CF.smith_line(K)


def test_s_job_on_an_e_cell_and_a_c_null(monkeypatch):
    """Referee round 64 (D3bk-9: `s_job`'s E and H0 paths): on a 2 + 2-draw E cell the labels are the one-bump fit f₀ plus the seed's
    first noise vectors, then the pair f₁ plus the next; the stored fits are used ('dfit' their distance from `one_fit`'s); z₂ is the
    residual's projection on Rφ_{j₂} toward its wall (R by least squares), the swap and the margin 1 − wall·u_{j₂} from Φᵀy, the cosine of
    Rφ_{j₂} with f₁ − f₀, for a pair of each sign (the wall ±1); on a C cell under H0 the law at the truth's fit at the noise s with the
    seed [seed, 9]."""
    monkeypatch.setattr(CF, 'M_LAW', 64)
    L, K, c0 = CF.setup()
    N = K.xs.size
    Phi = K.phi(c0)
    D, s, cb, seed = CF.e_jobs()[13]
    for A in ((0.8, 0.6), (-0.8, -0.6)):                                      # the runner-up's wall +1, then −1
        monkeypatch.setattr(CF, 'PAIR_C', A)
        f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
        c1, w1 = CF.one_fit(K, f1)
        f0 = w1 * K.phi(np.array([c1]))[:, 0]
        rng = np.random.default_rng(seed)
        ys = [f + s * rng.standard_normal(N) for f in (f0, f0, f1, f1)]
        fits = [CF.one_fit(K, y) for y in ys]
        fits = [(c + 1e-3 * (k % 2), w) for k, (c, w) in enumerate(fits)]          # odd draws: a stored fit off `one_fit`'s
        out = CF.s_job(('E', D, s, cb, 2, seed, (2, 2), fits))
        pd = CF.predicted_delay(K, f1, c0, *CF.PAIR_C, D, cbar=cb)
        j2, wall = pd[1], pd[2]
        j1n = int(np.argmax(np.abs(Phi.T @ f1)))

        def rphi(c, w):
            T = np.c_[K.phi(np.array([c]))[:, 0], w * K.dphi(np.array([c]))[:, 0]]
            return Phi[:, j2] - T @ np.linalg.lstsq(T, Phi[:, j2], rcond=None)[0]
        for k, (y, (c, w)) in enumerate(zip(ys, fits)):
            q = (out['rows0'] + out['rows1'])[k]
            assert abs(q['dfit'] - 1e-3 * (k % 2)) < 1e-9
            assert q['z'] == CF.draw_stats(K, c0, y, c, w, s)['z']
            Rp = rphi(c, w)
            r = y - w * K.phi(np.array([c]))[:, 0]
            assert abs(q['z2'] - wall * (r @ Rp) / (s * np.linalg.norm(Rp))) < 1e-9
            g = Phi.T @ y
            assert q['swap'] == (int(np.argmax(np.abs(g))) != j1n) and abs(q['mu'] - (1 - wall * g[j2] / np.abs(g).max())) < 1e-12
            assert abs(q['mu_u'] - (1 - wall * g[j2] / abs(g[j1n]))) < 1e-12        # referee round 73: the uncensored margin
        Rp = rphi(c1, w1)
        g = Phi.T @ f1
        assert out['j2'] == j2 and abs(out['mu_true'] - (1 - wall * g[j2] / np.abs(g).max())) < 1e-12
        assert abs(out['cos2'] - Rp @ (f1 - f0) / (np.linalg.norm(Rp) * np.linalg.norm(f1 - f0))) < 1e-9
    monkeypatch.setattr(CF, 'PAIR_C', (0.8, 0.6))
    Dc, sc, hc, sdc = next(j for j in CF.c_jobs() if j[2] == 0)
    outc = CF.s_job(('C', Dc, sc, CF.CBAR, hc, sdc, (2,), None))
    f1c = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, Dc, CF.CBAR), CF.PAIR_C)
    c1c, w1c = CF.one_fit(K, f1c)
    f0c = w1c * K.phi(np.array([c1c]))[:, 0]
    rngc = np.random.default_rng(sdc)
    for q in outc['rows0']:
        y = f0c + sc * rngc.standard_normal(N)
        assert q['z'] == CF.draw_stats(K, c0, y, *CF.one_fit(K, y), sc)['z'] and q['dfit'] == 0.0
    la, lb = CF.knot_walls(K, c0, c1c, w1c, sc * np.random.default_rng([sdc, 9]).standard_normal((64, N)))
    assert outc['q_law'] == float(np.quantile(np.minimum(la, lb), CF.ALPHA_C))
    assert outc['qA'] == float(np.quantile(la, CF.ALPHA_C / 2)) and outc['qB'] == float(np.quantile(lb, CF.ALPHA_C / 2))


RAW_CF = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'close_features_raw.npz')


def test_draw_stats_t_and_the_width_sign():
    """Referee round 71 (mutants Z02, D02): at a draw with a clear signal `draw_stats`' t-version is z·s/(‖r‖/√(N − 2)), r the residual
    at the fit, and the free-width z_σ has the sign of sign(ŵ)·z (to first order z_σ = sign(w)·z), for a pair of each sign."""
    L, K, c0 = CF.setup()
    N, s = K.xs.size, 0.003
    for A in ((0.8, 0.6), (0.8, -0.4)):
        f1 = X.bumps_target(K, CF.pair_centres(*A, 0.05, CF.CBAR), A)
        y = f1 + s * np.random.default_rng(3).standard_normal(N)
        c, w = CF.one_fit(K, y)
        st = CF.draw_stats(K, c0, y, c, w, s)
        r = y - w * K.phi(np.array([c]))[:, 0]
        assert abs(st['z']) > 3
        assert np.isclose(st['t'], st['z'] * s / (np.linalg.norm(r) / np.sqrt(N - 2)), rtol=1e-12, atol=0)
        assert np.sign(st['zs']) == np.sign(w) * np.sign(st['z'])


def test_s_job_margin_uncensored_on_a_swapping_tie_cell():
    """Referee round 71 (SM01, SM05): `s_job`'s mu_u = 1 − wall·g_{j₂}/|g_{j₁}| with g = Φᵀy of the noisy labels and j₁ the noiseless
    first entrant — uncensored: below 0 exactly on the draws whose first entrant swaps, at E's swapping tie cell (Δ = 0.05, s = 0.03,
    c̄ = 0.049)."""
    D, s, cb, seed = CF.e_jobs()[14]
    n = 60
    r = CF.s_job(('E', D, s, cb, 2, seed, (n, 0), None))
    L, K, c0 = CF.setup()
    N = K.xs.size
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    pd = CF.predicted_delay(K, f1, c0, *CF.PAIR_C, D, cbar=cb)
    Phi = K.phi(c0)
    j1 = int(np.argmax(np.abs(Phi.T @ f1)))
    rng = np.random.default_rng(seed)
    swaps = 0
    for q in r['rows0']:
        g = Phi.T @ (f0 + s * rng.standard_normal(N))
        assert np.isclose(q['mu_u'], 1 - pd[2] * g[pd[1]] / abs(g[j1]), rtol=1e-12, atol=1e-15)
        assert (q['mu_u'] < 0) == q['swap']
        swaps += q['swap']
    assert 5 <= swaps <= n - 5


def test_u_job_shares_e_s_law_samples():
    """Referee round 71: U's own-knot law is drawn from E's law samples, so the full knot's α quantile from them equals E's q_law exactly
    in every cell, and at the tie (c̄ ≥ 0.045) the own knot's equals it — the minimum's part there exactly 0 —; `u_job` reproduces the raw
    at the swapping tie cell, with the null's margin mu0 = 1 − wall·u_{j₂} at f₀ (raw)."""
    raw = CF.E.load(RAW_CF)
    Eq = {(e['D'], e['s'], e['cbar']): e for e in raw['E']}
    for u in raw['U']:
        e = Eq[(u['D'], u['s'], u['cbar'])]
        assert u['q_full'] == e['q_law']
        if u['cbar'] >= 0.045:
            assert u['q2'] == u['q_full']
    D, s, cb, seed = CF.e_jobs()[14]
    e = Eq[(D, s, cb)]
    fits = [(q['c1'], q['w1']) for q in e['rows0']] + [(q['c1'], q['w1']) for q in e['rows1']]
    u = CF.u_job((D, s, cb, seed, (len(e['rows0']), len(e['rows1'])), fits))
    ref = next(q for q in raw['U'] if (q['D'], q['s'], q['cbar']) == (D, s, cb))
    for k in ('q2', 'q_full', 'mu0', 'j2'):
        assert u[k] == ref[k], k
    assert list(u['t2_0']) == list(ref['t2_0']) and list(u['first_0']) == list(ref['first_0'])
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    g0, g1 = K.phi(c0).T @ (w1 * K.phi(np.array([c1]))[:, 0]), K.phi(c0).T @ f1
    pd = CF.predicted_delay(K, f1, c0, *CF.PAIR_C, D, cbar=cb)
    assert np.isclose(u['mu0'], 1 - pd[2] * g0[pd[1]] / abs(g0[int(np.argmax(np.abs(g1)))]), rtol=1e-12, atol=0)


def test_f_alloc_lines_on_synthetic_rows():
    """Referee round 71 (FA01–FA09): `f_alloc_lines` on synthetic H0 rows — the calibrated knot's rejections at the fit's-sign wall per pair
    and per offset, the sizes at α = 0.01 by sign, the calibrated size per offset and the walls test's components."""
    s, D = CF.SIGMAS_F[0], CF.DELTAS_F[CF.SIGMAS_F[0]][CF.H0_F]
    rows = [dict(p=p, tA=a, tB=b, z=z, w1=w) for p, a, b, z, w in
            ((0.01, 0.5, 1.5, 3.0, 1.0), (0.2, 2.0, 0.5, 0.0, 1.0), (0.03, 2.0, 1.0, -2.5, 1.0), (0.5, 2.0, 1.0, 0.1, -1.0))]
    F = [dict(A=A, cbar=cb, s=s, D=D, hyp=0, qA=1.0, qB=1.0, rows=rows) for A in CF.PAIRS_F for cb in CF.OFFSETS_F[:2]]
    keys = ('knot',)
    out = '\n'.join(CF.f_alloc_lines(F, {(A, k): [] for A in CF.PAIRS_F for k in keys},
                                     {(A, k, i): [0, 0] for A in CF.PAIRS_F for k in keys for i in (0, 1)}, dict(knot='knot (law)'), keys))
    for A in CF.PAIRS_F:
        assert (f"Pair {A} under H0: of the calibrated knot's 4 rejections 2 at the wall of the fit's sign (0.500; per offset "
                f"{CF.OFFSETS_F[0]:g}: 1/2, {CF.OFFSETS_F[1]:g}: 1/2).") in out
        assert f"pair {A}, H0: 0.2500, 0.2500" in out
    assert "The sizes at α = 0.01 over the 16 H0 draws" in out
    assert "knot (calibrated) 0.0000 [" in out and "z, one-sided toward the same sign 0.2500 [" in out
    assert "z, one-sided toward the opposite sign 0.2500 [" in out and "\\|z\\| 0.2500 [" in out
    assert f"c̄ = {CF.OFFSETS_F[0]:g}: 0.5000 [" in out and f"c̄ = {CF.OFFSETS_F[1]:g}: 0.5000 [" in out


def test_smith_line_with_a_closed_form_bound(monkeypatch):
    """Referee round 71 (SL01–SL05): `smith_line` with a closed-form bound sd = 2s/Δ³ (amplitudes unknown) and 3s/Δ² (known) prints the
    constants, the powers and the noise at which sd = Δ."""
    monkeypatch.setattr(CF, 'crb_sep', lambda K, cs, w, s, known=False:
                        (3.0 if known else 2.0) * s / abs(cs[1] - cs[0]) ** (2 if known else 3))
    L, K, c0 = CF.setup()
    line = CF.smith_line(K)
    assert "(0.8, 0.6), unknown: 2·s/Δ^3.00, sd = Δ at s ≈ 0.5Δ^4.00 (with sd ∝ s/Δ^3 exactly: s ≈ 0.5Δ^4, 0.5Δ^4" in line
    assert "(0.8, -0.4), known: 3·s/Δ^2.00, sd = Δ at s ≈ 0.333Δ^3.00 (with sd ∝ s/Δ^2 exactly: s ≈ 0.333Δ^3, 0.333Δ^3" in line


def test_equal_size_line_on_synthetic_rows():
    """Referee round 71: `equal_size_line` picks the level α′ whose H0 size is nearest the one-sided z's (here a unique jump: 10 of 200
    p-values at 0.06, z's size 10 of 200 → α′ = 0.0605), uses the fit's-sign wall for a same-sign pair and the other for an opposite-sign
    one, and pairs the Δ₅₀ ratios against z's."""
    s = CF.SIGMAS_F[0]
    Ds, cb = CF.DELTAS_F[s], CF.OFFSETS_F[0]
    for (A, key), wsg in itertools.product(((CF.PAIRS_F[0], 'pA'), (CF.PAIRS_F[1], 'pB')), (1.0, -1.0)):
        sg = float(np.sign(A[0] * A[1]))                                # round 82's ES03: z's size takes the fit's sign (ŵ < 0)
        h0 = [{key: 0.06 if k < 10 else 0.5, 'z': wsg * sg * (3.0 if k < 10 else 0.0), 'w1': wsg} for k in range(200)]
        F = [dict(A=A, cbar=cb, s=s, D=Ds[CF.H0_F], hyp=0, rows=h0)]
        pw, npw = (0.1, 0.3, 0.6, 0.8, 0.9, 1.0), (0.2, 0.4, 0.7, 0.9, 0.95, 1.0)
        for D, p, q in zip(Ds, pw, npw):
            F.append(dict(A=A, cbar=cb, s=s, D=D, hyp=1, np=q, rows=[{key: 0.0 if k < round(100 * p) else 0.9} for k in range(100)]))
        x = CF.crossing(list(Ds), list(pw)) / CF.crossing(list(Ds), list(npw))
        line = CF.equal_size_line(F, A, {(s, cb): 1.0})
        wall = "the fit's-sign wall" if key == 'pA' else 'the other wall'
        assert f"{wall} at p < α′ = 0.0605, its H0 size 0.0500 against the one-sided z's 0.0500" in line
        assert f"Δ₅₀ over the Neyman–Pearson bound's {x:.3f} [{x:.3f}, {x:.3f}]" in line
        assert f"smaller in {int(x < 1)} of 1" in line


def test_u_job_null_margin_with_a_negative_pair(monkeypatch):
    """Referee round 73 (U71-04 … -06; D3bt-5): `u_job`'s null margin mu0 = 1 − wall·u_{j₂} at f₀ with u over the null's own first
    entrant argmax|Φᵀf₀| — for a pair of negative sign (the wall −1) too —; on every E cell that entrant is f₁'s (so the round-73 change
    leaves the stored U unchanged)."""
    monkeypatch.setattr(CF, 'M_LAW', 64)
    L, K, c0 = CF.setup()
    N, Phi = K.xs.size, K.phi(c0)
    for D, s, cb, seed in CF.e_jobs():
        f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
        c1, w1 = CF.one_fit(K, f1)
        assert int(np.argmax(np.abs(Phi.T @ (w1 * K.phi(np.array([c1]))[:, 0])))) == int(np.argmax(np.abs(Phi.T @ f1)))
    D, s, cb, seed = CF.e_jobs()[13]
    for A in ((0.8, 0.6), (-0.8, -0.6)):
        monkeypatch.setattr(CF, 'PAIR_C', A)
        f1 = X.bumps_target(K, CF.pair_centres(*A, D, cb), A)
        c1, w1 = CF.one_fit(K, f1)
        f0 = w1 * K.phi(np.array([c1]))[:, 0]
        rng = np.random.default_rng(seed)
        ys = [f0 + s * rng.standard_normal(N) for _ in range(2)]
        u = CF.u_job((D, s, cb, seed, (2, 0), [CF.one_fit(K, y) for y in ys]))
        pd = CF.predicted_delay(K, f1, c0, *A, D, cbar=cb)
        g0 = Phi.T @ f0
        assert pd[2] == (1.0 if A[0] > 0 else -1.0)
        assert np.isclose(u['mu0'], 1 - pd[2] * g0[pd[1]] / abs(g0[int(np.argmax(np.abs(g0)))]), rtol=1e-12, atol=0)


def test_e_prime_lines_the_tie_and_the_parts_errors_on_the_raw():
    """Referee rounds 73, 78 (D3bt-1, -2, -5; D3by-2, -6): on the raw, E′ prints the minimum's part per cell with its paired standard
    error (a zero count's one-sided bound), the parts pooled per offset with their standard errors and the two lowest offsets' contrast,
    the full knot against the own knot draw by draw (never later, up to round-off), the H1 draws between the two knots at the tie's common
    threshold (none) and the bound on the cells' mean share there (one-sided and two-sided); without U the null's margin is '—'."""
    raw = CF.E.load(RAW_CF)
    md = '\n'.join(CF.e_prime_lines(raw['E'], raw['S'], raw['U']))
    Ue = {(u['D'], u['s'], u['cbar']): u for u in raw['U']}
    betw = n1 = exc_n = exc_all = 0
    exc_max = 0.0
    pool = {}
    for e in raw['E']:
        u = Ue[(e['D'], e['s'], e['cbar'])]
        d0 = np.array([q['dt'] for q in e['rows0']]); d1 = np.array([q['dt'] for q in e['rows1']])
        own = np.array(u['t2_1']) < u['q2']
        dd = (d1 < e['q_law']).astype(float) - own.astype(float)
        part = float(np.mean(d1 < e['q_law']) - np.mean(own))
        se = (f"0 discordant of {dd.size}: ≤ {1 - 0.05 ** (1 / dd.size):.4f} one-sided, {1 - 0.025 ** (1 / dd.size):.4f} two-sided at 95%"
              if not np.any(dd != 0)
              else f"{np.std(dd, ddof=1) / np.sqrt(dd.size):.3f}")
        assert f"({e['cbar']:g}, {e['D']:g}, {e['s']:g}) {part:+.4f} ({se})" in md
        if e['cbar'] <= 0.032:
            pool.setdefault(e['cbar'], []).append((part, np.std(dd, ddof=1) / np.sqrt(dd.size)))
        r = np.r_[d0, d1] / np.r_[u['t2_0'], u['t2_1']] - 1
        exc_n += int(np.sum(r > 0)); exc_all += r.size; exc_max = max(exc_max, float(r.max()))
        if u['q2'] == u['q_full']:
            betw += int(np.sum((d1 < u['q2']) & ~own)); n1 += own.size
    pp = []
    for cb in sorted(pool):
        m_, s_ = np.mean([x for x, _ in pool[cb]]), np.sqrt(sum(y ** 2 for _, y in pool[cb])) / len(pool[cb])
        pp.append((m_, s_))
        assert f"c̄ = {cb:g} {m_:+.4f} ({s_:.4f})" in md
    assert f"the two lowest offsets' pooled parts differ by z = {(pp[0][0] - pp[1][0]) / np.hypot(pp[0][1], pp[1][1]):+.2f}" in md
    assert exc_max < 1e-10 and n1 == 2400
    assert (f"exceeds it in {exc_n} of {exc_all} draws ({exc_all} with a finite own knot), by at most {exc_max:.1e} relative" in md
            and f"the threshold falls between the two knots of no H1 draw ({betw} of {n1} in those cells" in md and betw == 0)
    assert (f"the share at most {1 - 0.05 ** (1 / n1):.4f} (one-sided 95%; {1 - 0.025 ** (1 / n1):.4f} two-sided) for the cells' mean "
            "share") in md
    assert "(0 discordant of 400: ≤ 0.0075 one-sided, 0.0092 two-sided at 95%)" in md
    no_u = '\n'.join(CF.e_prime_lines(raw['E'], raw['S'], None))
    assert no_u.count('| —, ') == len(raw['E']) and 'The minimum' not in no_u


def test_the_size_by_distance_to_the_nearest_kernel():
    """Referee rounds 73, 75, 78, 82 (D3bt-4, D3bv-1, D3by-1, -4, D3cc-4): `size_by_distance` per pair — each H0 null in the interval of
    distances nearest an offset's, by its own centre's distance to the nearest kernel (an offset and its mirror about the tie pooled; a
    null shifted off its c̄ binned by its own centre, here into another bin than its label's) —, the bin printed as its range with its
    nulls' distances, H1 rows left out, a p-value exactly at α no rejection; E's nulls (the strong pair) beside F's and pooled; per line
    the trend over the nulls' own distances; the two nulls in one logistic model; `f_alloc_lines` prints it."""
    s, D = CF.SIGMAS_F[0], CF.DELTAS_F[CF.SIGMAS_F[0]][CF.H0_F]
    a = CF.ALPHA_C
    row = lambda p: dict(p=p, tA=1.0, tB=2.0, z=0.0, w1=1.0)
    A0, A1 = CF.PAIRS_F
    F = [dict(A=A0, cbar=0.02, c1=0.02, w1=1.4, s=s, D=D, hyp=0, qA=1.0, qB=1.0, rows=[row(0.01), row(a), row(0.5), row(0.5)]),
         dict(A=A0, cbar=0.08, c1=0.08, w1=1.4, s=s, D=D, hyp=0, qA=1.0, qB=1.0, rows=[row(0.005), row(0.03), row(0.5), row(0.5)]),
         dict(A=A0, cbar=0.02, c1=0.02, w1=1.4, s=s, D=D, hyp=1, qA=1.0, qB=1.0, rows=[row(0.0)] * 4),
         dict(A=A1, cbar=0.045, c1=0.0547, w1=0.43, s=s, D=D, hyp=0, qA=1.0, qB=1.0, rows=[row(0.001)] * 2 + [row(0.5)] * 2),
         dict(A=A1, cbar=0.032, c1=0.0417, w1=0.40, s=s, D=D, hyp=0, qA=1.0, qB=1.0, rows=[row(0.001)] + [row(0.5)] * 3)]
    out = '\n'.join(CF.size_by_distance(F))
    lo, hi = CF.cp_int(3, 8)
    assert f"Pair {A0}, F (8 H0 draws of 2 nulls, the nulls' |ŵ| 1.4000–1.4000)" in out
    assert f"[0.01, 0.026) (nulls at 0.0200–0.0200): 0.3750 [{lo:.4f}, {hi:.4f}], 0.1250" in out
    lo, hi = CF.cp_int(3, 8)
    assert f"Pair {A1}, F (8 H0 draws of 2 nulls, the nulls' |ŵ| 0.4000–0.4300)" in out
    assert f"[0.0385, 0.047) (nulls at 0.0417–0.0453): 0.3750 [{lo:.4f}, {hi:.4f}], 0.3750" in out and "[0.026, 0.0385)" not in out.split(f"Pair {A1}")[1]
    E_ = [dict(c1=c, rows0=[dict(p=p) for p in ps]) for c, ps in ((0.0, (0.01, 0.5, 0.5, 0.5)), (0.049, (0.02, 0.5, 0.5, 0.5)), (0.032, (0.02, 0.03, 0.5, 0.5)))]
    out = '\n'.join(CF.size_by_distance(F, E_))
    assert f"Pair {A0}, E (12 H0 draws of 3 nulls)" in out and f"Pair {A0}, F+E (20 H0 draws of 5 nulls, the nulls' |ŵ| 1.4000–1.4000)" in out
    assert "[0, 0.01) (nulls at 0.0000–0.0000): 0.2500 [" in out and "[0.047, 0.05] (nulls at 0.0490–0.0490): 0.2500 [" in out
    tie = '\n'.join(CF.size_by_distance(F, E_ + [dict(c1=0.05, rows0=[dict(p=0.5)] * 4)]))   # round 82's S02: h/2 in the last bin
    assert "[0.047, 0.05] (nulls at 0.0490–0.0500): 0.1250 [" in tie
    z, p = CF.trend_z([(0.0, 1, 4), (0.049, 1, 4), (0.032, 2, 4)])
    assert f"the trend over the nulls' own distances (Cochran–Armitage, binomial variance: the draws as the units) z = {z:+.2f} (p = {p:.3f})" in out
    assert "The two nulls in one logistic model of rejection at α over the nulls' own distances" in out
    keys = ('knot',)
    full = '\n'.join(CF.f_alloc_lines(F, {(A, k): [] for A in CF.PAIRS_F for k in keys},
                                      {(A, k, i): [0, 0] for A in CF.PAIRS_F for k in keys for i in (0, 1)}, dict(knot='knot (law)'), keys, E_))
    assert out in full


def test_equal_size_line_below_alpha_and_cells_that_do_not_count():
    """Referee round 73 (ES05, ES12, ES13): `equal_size_line` takes a level α′ below α where the size matches there, and skips a cell
    whose power never crosses ½ and a cell with a separation missing."""
    s = CF.SIGMAS_F[0]
    Ds = CF.DELTAS_F[s]
    A = CF.PAIRS_F[0]
    h0 = [{'pA': 0.04 if k < 10 else 0.5, 'z': 3.0 if k < 10 else 0.0, 'w1': 1.0} for k in range(200)]
    F = [dict(A=A, cbar=CF.OFFSETS_F[0], s=s, D=Ds[CF.H0_F], hyp=0, rows=h0)]
    npw = (0.2, 0.4, 0.7, 0.9, 0.95, 1.0)
    for cb, pw in ((CF.OFFSETS_F[0], (0.1, 0.3, 0.6, 0.8, 0.9, 1.0)), (CF.OFFSETS_F[1], (0.1, 0.1, 0.2, 0.2, 0.3, 0.4))):
        for D, p, q in zip(Ds, pw, npw):
            F.append(dict(A=A, cbar=cb, s=s, D=D, hyp=1, np=q, rows=[{'pA': 0.0 if k < round(100 * p) else 0.9} for k in range(100)]))
    for D, p, q in list(zip(Ds, (0.1, 0.3, 0.6, 0.8, 0.9, 1.0), npw))[:-1]:
        F.append(dict(A=A, cbar=CF.OFFSETS_F[2], s=s, D=D, hyp=1, np=q, rows=[{'pA': 0.0 if k < round(100 * p) else 0.9} for k in range(100)]))
    x = CF.crossing(list(Ds), [0.1, 0.3, 0.6, 0.8, 0.9, 1.0]) / CF.crossing(list(Ds), list(npw))
    line = CF.equal_size_line(F, A, {(s, CF.OFFSETS_F[0]): 1.0, (s, CF.OFFSETS_F[1]): 1.0, (s, CF.OFFSETS_F[2]): 1.0})
    assert "at p < α′ = 0.0405, its H0 size 0.0500 against the one-sided z's 0.0500" in line
    assert f"{x:.3f} [{x:.3f}, {x:.3f}]" in line and "smaller in" in line and "of 1" in line



def test_e_prime_counts_an_infinite_path_knot_and_the_first_entrants():
    """Referee round 75 (D3bv-6, -7): an infinite path knot beside a finite own knot is the worst exceedance, counted (not dropped), the
    denominator every draw; f₁'s and f₀'s first entrants (U's j1, j10) are counted where U records them."""
    D, s_ = CF.CELLS_E[0]
    cb = CF.OFFSETS_E[0]
    rows0 = [dict(dt=1.0 if i < 2 else 3.0, p=0.5, c1=0.0, w1=1.0) for i in range(10)]
    rows1 = [dict(dt=np.inf if i == 0 else 1.0, p=0.5, c1=0.0, w1=1.0) for i in range(4)]
    E_ = [dict(D=D, s=s_, cbar=cb, rows0=rows0, rows1=rows1, np=0.5, q_law=2.0)]
    row = lambda: dict(z=0.0, t=0.0, zs=0.0, z2=0.0, swap=False, mu=0.1, mu_u=0.1, dfit=0.0, w1=1.0)
    S = [dict(kind='E', D=D, s=s_, cbar=cb, rows0=[row() for _ in range(10)], rows1=[row() for _ in range(4)], mu_true=0.1, cos2=0.9)]
    U = [dict(D=D, s=s_, cbar=cb, q2=2.0, q_full=2.0, t2_0=[1.0] + [3.0] * 9, t2_1=[1.0, 1.0, 1.0, 3.0], first_0=[False] * 10,
              first_1=[False] * 4, mu0=0.1, j1=7, j10=7)]
    md = '\n'.join(CF.e_prime_lines(E_, S, U))
    assert "exceeds it in 1 of 14 draws (14 with a finite own knot), by at most inf relative" in md
    U2 = [dict(U[0], t2_1=[1.0, 1.0, 1.0, np.inf])]                     # referee round 78 (X04, X05): an infinite own knot left out
    assert "exceeds it in 1 of 14 draws (13 with a finite own knot), by at most inf relative" in '\n'.join(CF.e_prime_lines(E_, S, U2))
    assert "f₁'s and f₀'s first entrants (S's normalisation of μ̂, U's of μ₀) coincide in 1 of 1 cells" in md
    U[0]['j10'] = 8
    assert "coincide in 0 of 1 cells" in '\n'.join(CF.e_prime_lines(E_, S, U))
    U[0].pop('j10')
    assert "f₁'s and f₀'s first entrants (S's" not in '\n'.join(CF.e_prime_lines(E_, S, U))


def test_smith_line_with_a_non_integer_power_and_the_bound_itself(monkeypatch):
    """Referee round 75 (D3bv-4): with sd = 2s/Δ^{2.6} the line prints the power 2.60 and the integer power's constants at round(2.6) =
    3 from both separations; unpatched, its power and constant are those of `crb_sep` at SMITH_D."""
    L, K, c0 = CF.setup()
    line = CF.smith_line(K)
    for A in CF.SMITH_PAIRS:
        for known in (False, True):
            sd = [CF.crb_sep(K, CF.pair_centres(*A, D), np.array(A), 1.0, known) for D in CF.SMITH_D]
            p = float(np.log(sd[1] / sd[0]) / np.log(CF.SMITH_D[0] / CF.SMITH_D[1]))
            c = sd[0] * CF.SMITH_D[0] ** p
            assert f"({A[0]:g}, {A[1]:g}), {'known' if known else 'unknown'}: {c:.3g}·s/Δ^{p:.2f}, sd = Δ at s ≈ {1 / c:.3g}Δ^{p + 1:.2f}" in line
    monkeypatch.setattr(CF, 'crb_sep', lambda K, cs, w, s, known=False: 2.0 * s / abs(cs[1] - cs[0]) ** 2.6)
    line = CF.smith_line(K)
    ci = [2.0 / D ** 2.6 * D ** 3 for D in CF.SMITH_D]
    assert (f"unknown: 2·s/Δ^2.60, sd = Δ at s ≈ 0.5Δ^3.60 (with sd ∝ s/Δ^3 exactly: s ≈ {1 / ci[0]:.3g}Δ^4, {1 / ci[1]:.3g}Δ^4 from "
            f"Δ = {CF.SMITH_D[0]:g}, {CF.SMITH_D[1]:g})") in line


def test_f_alloc_lines_on_rows_that_differ_by_pair_and_offset():
    """Referee rounds 75, 78 (FA01, FA03–FA07): per pair and offset its own rejections at the fit's-sign wall; H1 rows left out of every H0
    statistic; the sizes at α = 0.01 by the fit's sign (a row with ŵ < 0 flips its z's side) over the H0 draws only; the size per offset
    (both pairs pooled, by the label c̄); the walls test's components at their own thresholds qA ≠ qB, under H0 and under H1 at the two
    smallest separations only."""
    s, D = CF.SIGMAS_F[0], CF.DELTAS_F[CF.SIGMAS_F[0]][CF.H0_F]
    D0, D2 = CF.DELTAS_F[s][0], CF.DELTAS_F[s][2]
    A0, A1 = CF.PAIRS_F
    r = lambda p, a, b, z, w: dict(p=p, tA=a, tB=b, z=z, w1=w)
    F = [dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=D, hyp=0, qA=1.0, qB=2.0, q_law=0.5,
              rows=[r(0.001, 0.5, 1.5, 3.0, 1.0), r(0.002, 2.0, 0.5, 0.0, 1.0)]),
         dict(A=A0, cbar=0.02, c1=0.02, w1=1.4, s=s, D=D, hyp=0, qA=1.0, qB=2.0, q_law=0.5,
              rows=[r(0.003, 0.5, 1.5, -3.0, -1.0), r(0.5, 2.0, 0.5, 0.0, 1.0)]),
         dict(A=A1, cbar=0.0, c1=0.0, w1=0.4, s=s, D=D, hyp=0, qA=1.0, qB=2.0, q_law=0.5,
              rows=[r(0.004, 2.0, 0.5, -3.0, 1.0), r(0.5, 2.0, 0.5, 2.7, 1.0)]),
         dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=D0, hyp=1, qA=1.0, qB=2.0, q_law=0.5,
              rows=[r(0.0, 0.5, 1.5, 9.0, 1.0)] * 4 + [r(0.0, 1.5, 2.5, 9.0, 1.0)] * 4),
         dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=D2, hyp=1, qA=1.0, qB=2.0, q_law=0.5, rows=[r(0.0, 5.0, 5.0, 9.0, 1.0)] * 5)]
    keys = ('knot',)
    out = '\n'.join(CF.f_alloc_lines(F, {(A, k): [] for A in CF.PAIRS_F for k in keys},
                                     {(A, k, i): [0, 0] for A in CF.PAIRS_F for k in keys for i in (0, 1)}, dict(knot='knot (law)'), keys))
    assert f"Pair {A0} under H0: of the calibrated knot's 3 rejections 2 at the wall of the fit's sign (0.667; per offset 0: 1/2, 0.02: 1/1)." in out
    assert f"Pair {A1} under H0: of the calibrated knot's 1 rejections 0 at the wall of the fit's sign (0.000; per offset 0: 0/1)." in out
    assert "The sizes at α = 0.01 over the 6 H0 draws" in out
    k = lambda n, m=6: f"{n / m:.4f} [{CF.cp_int(n, m)[0]:.4f}, {CF.cp_int(n, m)[1]:.4f}]"
    assert f"knot (calibrated) {k(4)}" in out and f"z, one-sided toward the same sign {k(3)}" in out
    assert f"z, one-sided toward the opposite sign {k(1)}" in out and f"\\|z\\| {k(4)}" in out
    assert f"c̄ = 0: {k(3, 4)}, 0.7500" in out and f"c̄ = 0.02: {k(1, 2)}, 0.5000" in out
    assert (f"pair {A0}, H0: 0.5000, 1.0000; pair {A0}, H1: 0.5000, 0.5000; pair {A1}, H0: 0.0000, 1.0000") in out


def test_f_report_prints_the_noise_range_and_the_equal_size_lines():
    """Referee round 75 (D3bv-4): `f_lines` prints the range of s/|ŵ| over F's draws and each pair's equal-size line."""
    s, cb = CF.SIGMAS_F[1], CF.CBAR
    Ds = CF.DELTAS_F[s]
    F = []
    for A in CF.PAIRS_F:
        for D in Ds:
            wf = 0.5 if A == CF.PAIRS_F[0] else -2.0                      # round 82's FR03: |ŵ|, a negative ŵ
            rows = [dict(dt=9.0, p=0.5, pA=0.5, pB=0.5, tA=9.0, tB=9.0, z=0.0, zs=0.0, w1=wf) for i in range(20)]
            F.append(dict(A=A, cbar=cb, s=s, D=D, hyp=1, rows=rows, q_law=2.0, qA=1.0, qB=1.0, np=0.5, c1=cb, w1=wf))
        rows0 = [dict(dt=9.0, p=0.5, pA=0.5, pB=0.5, tA=9.0, tB=9.0, z=0.0, zs=0.0, w1=0.25) for i in range(40)]   # FR01: the extreme in H0
        F.append(dict(A=A, cbar=cb, s=s, D=Ds[CF.H0_F], hyp=0, rows=rows0, q_law=2.0, qA=1.0, qB=1.0, np=0.0, w1=1.0, c1=cb))
    md = '\n'.join(CF.f_lines(F))
    assert f"s/|ŵ|, over F's draws (each at its own fit; referee rounds 67 and 71): {s / 2.0:.4f}–{s / 0.25:.4f}" in md
    assert md.count('at equal size (referee round 71)') == 2


def test_trend_z_and_logit_fit_against_closed_forms():
    """Referee round 78: the Cochran–Armitage statistic on two groups is the pooled two-proportion z with the score's sign; the logistic
    fit of a saturated two-group model recovers the groups' log-odds, and its deviance is zero."""
    z, p = CF.trend_z([(0.0, 30, 100), (1.0, 10, 100)])
    pb = 40 / 200
    zz = (10 / 100 - 30 / 100) / np.sqrt(pb * (1 - pb) * (1 / 100 + 1 / 100))
    assert np.isclose(z, zz, rtol=1e-12) and np.isclose(p, math.erfc(abs(zz) / np.sqrt(2)), rtol=1e-12)
    assert CF.trend_z([(0.0, 0, 10), (1.0, 0, 10)]) == (None, None) and CF.trend_z([(0.0, 3, 10)]) == (None, None)
    b, cov, dv = CF.logit_fit(np.array([[1.0, 0.0], [1.0, 1.0]]), np.array([30.0, 10.0]), np.array([100.0, 100.0]))
    assert np.allclose(b, [np.log(30 / 70), np.log(10 / 90) - np.log(30 / 70)], atol=1e-10) and abs(dv) < 1e-9
    assert np.isclose(cov[1, 1], 1 / 30 + 1 / 70 + 1 / 10 + 1 / 90, rtol=1e-8)


def test_crb_sep_against_a_finite_difference_bound():
    """Referee round 78 (SL05): `crb_sep` is the separation's Cramér–Rao bound s²gᵀ(GᵀG)⁻¹g with G the samples' Jacobian by central
    differences, the amplitudes unknown (G in w₁, w₂, c₁, c₂) or known (G in c₁, c₂)."""
    L, K, c0 = CF.setup()
    for A, D in (((0.8, 0.6), 0.05), ((0.8, -0.4), 0.025)):
        cs, w = CF.pair_centres(*A, D), np.array(A)
        f = lambda th: K.phi(th[2:]) @ th[:2]
        th, h = np.r_[w, cs], 1e-5
        G = np.column_stack([(f(th + h * e) - f(th - h * e)) / (2 * h) for e in np.eye(4)])
        for known, cols, g in ((False, slice(0, 4), np.array([0.0, 0.0, -1.0, 1.0])), (True, slice(2, 4), np.array([-1.0, 1.0]))):
            Gc = G[:, cols]
            ref = 0.01 * np.sqrt(g @ np.linalg.inv(Gc.T @ Gc) @ g)
            assert np.isclose(CF.crb_sep(K, cs, w, 0.01, known), ref, rtol=1e-6)


def test_null_margin_against_its_definition_and_a_monte_carlo():
    """Referee round 82: `null_margin` at a one-bump null — j₂ the kernel other than the first entrant whose dual at the first entry is
    nearest its wall, μ₀ = 1 − |u_{j₂}| — and its delta-method sd(μ̂)/s against the standard deviation of μ̂ over noise draws at a small
    s; for the pair (0.8, 0.6) j₂ is `predicted_delay`'s runner-up; a passed Φ changes nothing."""
    L, K, c0 = CF.setup()
    Phi = K.phi(c0)
    xi = np.random.default_rng(5).standard_normal((20000, K.xs.size))
    for A, cb, D in (((0.8, 0.6), 0.049, 0.025), ((0.8, -0.4), 0.032, 0.05), ((0.8, 0.6), 0.0, 0.0175)):
        f1 = X.bumps_target(K, CF.pair_centres(*A, D, cb), A)
        c1, w1 = CF.one_fit(K, f1)
        f0 = w1 * K.phi(np.array([c1]))[:, 0]
        j2, mu0, sd = CF.null_margin(K, c0, c1, w1)
        g = Phi.T @ f0
        j1 = int(np.argmax(np.abs(g)))
        u = g / abs(g[j1])
        assert j2 == min((j for j in range(c0.size) if j != j1), key=lambda j: 1 - abs(u[j])) and abs(mu0 - (1 - abs(u[j2]))) < 1e-15
        s = 1e-4
        G = (f0[None, :] + s * xi) @ Phi
        muh = 1 - np.sign(u[j2]) * G[:, j2] / np.abs(G[:, j1])
        assert abs(np.std(muh) / s / sd - 1) < 0.03
        if A == (0.8, 0.6):
            assert j2 == CF.predicted_delay(K, f1, c0, *A, D, cbar=cb)[1]
        assert CF.null_margin(K, c0, c1, w1, Phi) == (j2, mu0, sd)


def test_quasi_logit_perm_trend_and_a_separated_design():
    """Referee round 82: `quasi_logit`'s dispersion is the Pearson χ² over m − p and scales the standard errors by √φ; `perm_trend`'s p-value
    counts the permutations of the scores over the groups as defined; `logit_fit` raises on a separated design."""
    import pytest
    x = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05])
    k, n = np.array([20, 25, 18, 10, 6, 9]), np.array([400, 400, 300, 400, 400, 300])
    Xd = np.c_[np.ones_like(x), x]
    b, se, phi, dv = CF.quasi_logit(Xd, k, n)
    b0, cov, dv0 = CF.logit_fit(Xd, k, n)
    p = 1 / (1 + np.exp(-(Xd @ b0)))
    assert np.allclose(b, b0) and np.isclose(dv, dv0)
    assert np.isclose(phi, np.sum((k - n * p) ** 2 / (n * p * (1 - p))) / 4, rtol=1e-12) and np.allclose(se, np.sqrt(np.diag(cov) * phi))
    gr = list(zip(x, k, n))
    z, pp = CF.perm_trend(gr, n_perm=200, seed=3)
    rng = np.random.default_rng(3)
    hits = sum(abs(CF.trend_z([(xx, g[1], g[2]) for xx, g in zip(rng.permutation(x), gr)])[0]) >= abs(z) - 1e-12 for _ in range(200))
    assert z == CF.trend_z(gr)[0] and pp == (1 + hits) / 201 and 0 < hits < 200
    assert CF.perm_trend([(0.0, 0, 10), (1.0, 0, 10)]) == (None, None)
    with pytest.raises(np.linalg.LinAlgError):
        CF.logit_fit(np.c_[np.ones(4), [0.0, 1.0, 2.0, 3.0]], np.array([0, 0, 10, 10]), np.array([10, 10, 10, 10]))
    with pytest.raises(np.linalg.LinAlgError, match='no convergence'):  # round 82's LF1: the check itself, before any singular step
        CF.logit_fit(Xd, k, n, it=2)


def test_trend_z_with_unequal_groups_and_the_logistic_line_against_scipy():
    """Referee round 82 (T01, L03, L06, S06–S10, S15): the Cochran–Armitage z with unequal group sizes against its score form; on a
    non-saturated design of both pairs' nulls the logistic line's three numbers — the interaction's z, its likelihood-ratio p and the
    common slope's z — against an independent maximum-likelihood fit (scipy) with the likelihood ratio from the log-likelihoods."""
    from scipy.optimize import minimize
    from scipy.stats import chi2
    x = np.array([0.0, 0.02, 0.032, 0.045, 0.049]); k = np.array([30, 12, 25, 9, 4]); n = np.array([500, 200, 600, 300, 150])
    z, p = CF.trend_z(list(zip(x, k, n)))
    pb = k.sum() / n.sum()
    zz = np.sum(x * (k - n * pb)) / np.sqrt(pb * (1 - pb) * (np.sum(n * x ** 2) - np.sum(n * x) ** 2 / n.sum()))
    assert np.isclose(z, zz, rtol=1e-12) and np.isclose(p, math.erfc(abs(zz) / np.sqrt(2)), rtol=1e-12)
    s, D = CF.SIGMAS_F[0], CF.DELTAS_F[CF.SIGMAS_F[0]][CF.H0_F]
    A0, A1 = CF.PAIRS_F
    spec = [(A0, 0.0, 0.0, 22, 400), (A0, 0.02, 0.02, 18, 400), (A0, 0.032, 0.032, 21, 400), (A0, 0.045, 0.045, 9, 400),
            (A0, 0.049, 0.049, 12, 400), (A1, 0.0, 0.003, 19, 400), (A1, 0.02, 0.021, 14, 400), (A1, 0.032, 0.0335, 8, 400),
            (A1, 0.045, 0.0465, 5, 400), (A1, 0.049, 0.0493, 9, 400)]
    F = [dict(A=A, cbar=cb, c1=c, w1=1.4 if A == A0 else 0.4, s=s, D=D, hyp=0, qA=1.0, qB=1.0,
              rows=[dict(p=0.01 if i < kk else 0.5, tA=1.0, tB=2.0, z=0.0, w1=1.0) for i in range(nn)]) for A, cb, c, kk, nn in spec]
    out = '\n'.join(CF.size_by_distance(F))
    h = 2.0 / (X.D_C - 1)
    d = np.array([CF.null_distance(c, h) for _, _, c, _, _ in spec]); kk = np.array([t[3] for t in spec], float)
    nn = np.array([t[4] for t in spec], float); g = np.array([1.0 if t[0] == A1 else 0.0 for t in spec])
    if sorted([A0, A1])[1] == A0:
        g = 1 - g
    def fit(Xm):
        nll = lambda b: -np.sum(kk * (Xm @ b) - nn * np.log1p(np.exp(Xm @ b)))
        r = minimize(nll, np.zeros(Xm.shape[1]), method='BFGS', options=dict(gtol=1e-10))
        pr = 1 / (1 + np.exp(-(Xm @ r.x)))
        cov = np.linalg.inv(Xm.T @ ((nn * pr * (1 - pr))[:, None] * Xm))
        return r.x, cov, -r.fun
    bf, cf, lf = fit(np.c_[np.ones_like(d), d, g, g * d])
    br, cr, lr_ = fit(np.c_[np.ones_like(d), d, g])
    lr = 2 * (lf - lr_)
    prs = sorted({A0, A1})
    assert (f"the slopes' difference (pair {prs[1]} against {prs[0]}) z = {bf[3] / np.sqrt(cf[3, 3]):+.2f}, its likelihood-ratio p = "
            f"{chi2.sf(lr, 1):.3f}; the common slope (without the interaction) z = {br[1] / np.sqrt(cr[1, 1]):+.2f}.") in out


def test_per_offset_sizes_count_p_strictly_below_alpha_and_the_components_keep_two_separations():
    """Referee round 82 (F02, F06): a p-value exactly at α or at 0.01 is no rejection in the per-offset line; the walls test's components
    under H1 take the two smallest separations only (an H1 row at the second counted, at the third not)."""
    s = CF.SIGMAS_F[0]
    Ds = CF.DELTAS_F[s]
    A0, A1 = CF.PAIRS_F
    row = lambda p, tA=9.0, tB=9.0: dict(p=p, tA=tA, tB=tB, z=0.0, w1=1.0)
    F = [dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=Ds[CF.H0_F], hyp=0, qA=1.0, qB=1.0,
              rows=[row(CF.ALPHA_C), row(0.01), row(0.009), row(0.5)]),
         dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=Ds[1], hyp=1, qA=1.0, qB=1.0, rows=[row(0.5, 0.5, 9.0)] * 3 + [row(0.5)]),
         dict(A=A0, cbar=0.0, c1=0.0, w1=1.4, s=s, D=Ds[2], hyp=1, qA=1.0, qB=1.0, rows=[row(0.5)] * 4)]
    keys = ('knot',)
    out = '\n'.join(CF.f_alloc_lines(F, {(A, k_): [] for A in CF.PAIRS_F for k_ in keys},
                                     {(A, k_, i): [0, 0] for A in CF.PAIRS_F for k_ in keys for i in (0, 1)}, dict(knot='knot (law)'), keys))
    lo, hi = CF.cp_int(2, 4)
    assert f"c̄ = 0: 0.5000 [{lo:.4f}, {hi:.4f}], 0.2500" in out
    assert f"pair {A0}, H1: 0.7500, 0.0000" in out


def test_the_report_on_the_raw_prints_es_lines_and_the_margin_section():
    """Referee round 82 (R01, R02): from the raw, the report carries E's size lines and the pooled ones beside F's (E_ passed through
    `report` and `f_lines`), and the margin section (F″) pools F's and E's nulls by the band."""
    raw = CF.E.load(os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'close_features_raw.npz'))
    keys = ('A', 'B', 'C', 'E', 'S', 'F', 'U', 'G')
    md = CF.report(*[raw[k] if k in raw else None for k in keys])
    nE = sum(len(e['rows0']) for e in raw['E'])
    assert f"Pair (0.8, 0.6), E ({nE} H0 draws of {len(raw['E'])} nulls" in md and "Pair (0.8, 0.6), F+E (" in md
    ns = CF.margin_nulls(raw['F'], raw['E'])
    band = [x for x in ns if CF.in_band(x) and x['A'] == (0.8, 0.6)]
    assert "### F″. The calibrated knot's size against the runner-up's margin at the null (referee round 82)" in md
    assert f"pair (0.8, 0.6) in the band: {len(band)} nulls, {CF.fmt_rate(*CF.pooled(band))}" in md
    assert any(x['src'] == 'E' for x in band) and any(x['src'] == 'F' for x in band)


def test_margin_lines_on_synthetic_nulls(monkeypatch):
    """Referee round 82: `margin_lines` — per null its μ₀/sd(μ̂) (`null_margin`, here faked), its calibrated size, the law at the truth's
    rejections and how many of them the calibrated test makes, the fitted margins' medians; pooled by the band (its bounds exclusive)
    per pair and for both; G's nulls against the band's and the distance bin's predictions from F's and E's nulls."""
    fake = {0.02: 0.03, 0.0: 0.05, 0.045: 0.005, 0.0335: 0.02, 0.049: 0.001, 0.021: 0.028, 0.0015: 0.05}
    monkeypatch.setattr(CF, 'null_margin', lambda K, c0, c, w, Phi=None: (7, fake[round(c, 4)], 1.0))
    s, D = 0.01, 0.025
    A0, A1 = CF.PAIRS_F
    def null(A, c, k, kt, n=10, src='F'):
        rows = [dict(p=0.01 if i < k else 0.5, dt=0.5 if i < kt else 2.0, c1=c, w1=1.0) for i in range(n)]
        r = dict(A=A, cbar=c, c1=c, w1=1.4, s=s, D=D, hyp=0, q_law=1.0)
        return dict(r, rows0=rows) if src == 'E' else dict(r, rows=rows)
    F = [null(A0, 0.02, 1, 3), null(A0, 0.0, 3, 3), null(A1, 0.045, 2, 1), null(A1, 0.0335, 0, 2)]
    E_ = [null(A0, 0.049, 1, 1, src='E')]
    G = [null(A0, 0.021, 4, 6, n=20), null(A1, 0.0015, 2, 1, n=20)]
    md = '\n'.join(CF.margin_lines(F, E_, G))
    assert "| (0.8, 0.6) | F | 0.01 | 0.02 | 0.0200 | 7 | 0.0300 | 1.0000 | 3.00 | 10 | " in md
    assert f"0.1000 [{CF.cp_int(1, 10)[0]:.4f}, {CF.cp_int(1, 10)[1]:.4f}] | 0.3000 | 1 of 3; 0.0300, 0.0300 |" in md
    assert f"pair (0.8, 0.6) in the band: 1 nulls, {CF.fmt_rate(1, 10)}, against the law at the truth {CF.fmt_rate(3, 10)}" in md
    assert f"pair (0.8, 0.6) outside it: 2 nulls, {CF.fmt_rate(4, 20)}, against the law at the truth {CF.fmt_rate(4, 20)}" in md
    assert f"pair (0.8, -0.4) in the band: 1 nulls, {CF.fmt_rate(0, 10)}" in md and f"both pairs in the band: 2 nulls, {CF.fmt_rate(1, 20)}" in md
    assert f"| (0.8, 0.6) | 0.01 | 0.021 | 0.025 | 0.0210 | 2.80 | yes | 0.1000 | 0.1000 | {CF.fmt_rate(4, 20)} | {CF.fmt_rate(6, 20)} |" in md
    assert f"| (0.8, -0.4) | 0.01 | 0.0015 | 0.025 | 0.0015 | 5.00 | no | 0.2000 | — | {CF.fmt_rate(2, 20)} | {CF.fmt_rate(1, 20)} |" in md
    assert f"G pooled: in the band 1 nulls, {CF.fmt_rate(4, 20)}, against the law at the truth {CF.fmt_rate(6, 20)}; outside it 1 nulls" in md
    assert CF.margin_lines(None) == [] and CF.margin_lines([]) == []


def test_g_jobs_and_f_jobs_seed(monkeypatch, tmp_path):
    """Referee round 82: G's jobs are its cells' H0 nulls with R0_G draws and their own seeds, which `f_job` uses instead of `f_seed`;
    `main --only G` runs `part_G` and keeps the other parts."""
    J = CF.g_jobs()
    assert len(J) == len(CF.G_CELLS) and all(j[4] == 0 and j[5] == CF.R0_G for j in J) and [j[6] for j in J] == [820000 + 10 * i for i in range(len(J))]
    assert [tuple(j[:4]) for j in J] == [(A, cb, s, D) for A, cb, s, D in CF.G_CELLS]
    monkeypatch.setattr(CF, 'f_seed', lambda *a: (_ for _ in ()).throw(AssertionError('f_seed called')))
    r = CF.f_job(((0.8, 0.6), 0.032, 0.03, 0.05, 0, 2, 820000))
    assert len(r['rows']) == 2 and r['hyp'] == 0
    saved, got = [], []
    monkeypatch.setattr(CF.E, 'load', lambda p: dict(F=['f'], E=['e']))
    monkeypatch.setattr(CF.E, 'save', lambda p, **kw: saved.append(kw))
    monkeypatch.setattr(CF, 'report', lambda *a: '')
    monkeypatch.setattr(CF, 'part_G', lambda: got.append(1) or ['g'])
    monkeypatch.setattr(CF, 'OUT', str(tmp_path))
    monkeypatch.setattr(sys, 'argv', ['exp_close_features.py', '--only', 'G'])
    CF.main()
    assert got == [1] and saved[-1] == dict(F=['f'], E=['e'], G=['g'])

def test_bin_index_null_margin_ties_and_signs():
    """Round 82's mutants (BI1, NM1, NM3): `bin_index` puts a distance on an inner edge into the upper bin and closes the last; `null_margin`
    takes the duals' moduli (a negative null gives the positive one's margin) and keeps a tied kernel as the runner-up (its margin 0)."""
    edges = [0.0, 0.01, 0.03, 0.05]
    assert [CF.bin_index(d, edges) for d in (0.0, 0.01, 0.0299, 0.03, 0.05, 0.0501, -1e-9)] == [0, 1, 1, 2, 2, None, None]

    class Kf:
        def phi(self, cs):
            return np.array([[1.0], [2.0], [0.5]])
    Phi = np.array([[1.0, 1.0, 0.2], [2.0, 2.0, 0.1], [0.5, 0.5, 0.9]])   # kernels 0 and 1 alike: a tie
    j2, mu0, sd = CF.null_margin(Kf(), None, 0.0, 1.0, Phi)
    assert (j2, mu0, sd) == (1, 0.0, 0.0)
    assert CF.null_margin(Kf(), None, 0.0, -1.0, Phi) == (j2, mu0, sd)
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(0.8, 0.6, 0.025, 0.049), (0.8, 0.6))
    c1, w1 = CF.one_fit(K, f1)
    assert CF.null_margin(K, c0, c1, -w1) == CF.null_margin(K, c0, c1, w1)


def test_e_prime_discordances_by_sign_and_the_pooled_bound_needs_every_cell():
    """Round 82's mutants (D01, D02, B02, B03 of both lists, P02 of both lists): the minimum's part counts gains and losses alike as
    discordant — a cell with a gain only and one with a loss only print their standard errors, a cell with neither its count and the
    bound 1 − 0.05^{1/n} (one-sided), 1 − 0.025^{1/n} (two-sided) —; the pooled bound on the share between the knots is printed only
    when no cell where the quantiles coincide has a discordant draw."""
    D, s_ = CF.CELLS_E[0]
    row = lambda: dict(z=0.0, t=0.0, zs=0.0, z2=0.0, swap=False, mu=0.1, mu_u=0.1, dfit=0.0, w1=1.0)
    # (c̄, the H1 draws' path knots, their own knots) at q_law = q2 = 2: no discordance, a gain (path below, own above), a loss
    cells = ((0.0, [1.0] * 4, [1.0] * 4), (0.02, [1.0] * 4, [1.0, 1.0, 1.0, 3.0]), (0.032, [3.0, 1.0, 1.0, 1.0], [1.0] * 4))

    def md_of(cells):
        E_ = [dict(D=D, s=s_, cbar=cb, rows0=[dict(dt=3.0, p=0.5, c1=0.0, w1=1.0) for _ in range(10)],
                   rows1=[dict(dt=t, p=0.5, c1=0.0, w1=1.0) for t in dt], np=0.5, q_law=2.0) for cb, dt, _ in cells]
        S = [dict(kind='E', D=D, s=s_, cbar=cb, rows0=[row() for _ in range(10)], rows1=[row() for _ in range(4)], mu_true=0.1, cos2=0.9)
             for cb, _, _ in cells]
        U = [dict(D=D, s=s_, cbar=cb, q2=2.0, q_full=2.0, t2_0=[3.0] * 10, t2_1=t2, first_0=[False] * 10, first_1=[False] * 4, mu0=0.1,
                  j1=7, j10=7) for cb, _, t2 in cells]
        return '\n'.join(CF.e_prime_lines(E_, S, U))
    md = md_of(cells)
    se = np.std([0.0, 0.0, 0.0, 1.0], ddof=1) / 2
    assert (f"(0, {D:g}, {s_:g}) +0.0000 (0 discordant of 4: ≤ {1 - 0.05 ** 0.25:.4f} one-sided, {1 - 0.025 ** 0.25:.4f} two-sided "
            "at 95%)") in md
    assert f"(0.02, {D:g}, {s_:g}) +0.2500 ({se:.3f})" in md and f"(0.032, {D:g}, {s_:g}) -0.2500 ({se:.3f})" in md
    assert "no discordant draw" not in md                               # two of the three cells have one
    md = md_of([(cb, [1.0] * 4, [1.0] * 4) for cb, _, _ in cells])
    assert (f"no discordant draw: the share at most {1 - 0.05 ** (1 / 12):.4f} (one-sided 95%; {1 - 0.025 ** (1 / 12):.4f} two-sided) "
            "for the cells' mean share") in md
