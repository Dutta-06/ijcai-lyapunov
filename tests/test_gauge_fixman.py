"""
Proposition D‴.38 (PARADIGM.md §5d): a Fixman correction for the replicator under dither. Checked here: the loop's correction equals
(τ/2)∂_βℓ along the manifold of exact fits (ℓ = log det ∇²_ΘJ, by finite differences of the Hessian's determinant); it vanishes without
dither and keeps Σβ = m; in the slow-replicator limit it cancels the dither's drift; with κ_F = 2 the corrected field's zero is the unique
maximiser W♯ of ℓ (the D-optimal exact fit); near a saturation the one-dimensional zero-flux density behaves like δ^{κ_F − 1}, and below
the floor's cap the exponent falls back as stated; the tracked average relaxes at its time constant; kappa() is the leak's projected stiffness;
the long runs' and the paired analysis' statistics are right on synthetic records; the correction's general form,
through the exact fit's implicit function and computed from the current point only, reduces to the closed form.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_gauge_dither as GDD                                           # noqa: E402
import exp_gauge_drift as X                                              # noqa: E402
import exp_gauge_fixman as Q                                             # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_correction_is_the_log_determinants_derivative_along_the_manifold():
    net, base = E.setup()
    tau = 1e-4
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau)
    beta = 3 * np.array([0.30, 0.35, 0.35]); hb = 1e-3
    grad = np.array([(X.log_det(beta + hb * e) - X.log_det(beta - hb * e)) / (2 * hb) for e in np.eye(3)])
    c = loop.fixman_correction(E.exact_rep(net, beta / 3), beta)
    assert np.abs(c / (0.5 * tau * grad) - 1).max() < 1e-4


def test_the_correction_vanishes_without_dither_and_keeps_the_mass():
    net, base = E.setup()
    rng = np.random.default_rng(3)
    Th = E.exact_rep(net, np.ones(3) / 3) + 0.05 * rng.standard_normal((3, 3))
    beta = np.array([0.8, 1.1, 1.1])
    d0 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2).drift(Th, beta)
    dz = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0).drift(Th, beta)
    dc = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4).drift(Th, beta)
    assert np.array_equal(d0[1], dz[1]) and np.array_equal(d0[0], dc[0])
    assert abs(dc[1].sum()) < 1e-15 and np.abs(dc[1] - d0[1]).max() > 0


def test_the_slow_limit_cancellation():
    w = np.array([1 / 3, 1 / 3])
    b = Q.drift(w, 1e-5, 0.0)[0]; bc = Q.drift(w, 1e-5, 1.0)[0]
    assert np.abs(bc).max() < 3e-4 * np.abs(b).max()                        # the true ratio at γ = 10⁻⁵ is ≈ 7·10⁻⁵ (referee round 28)


def test_the_gain_selects_the_d_optimal_fit():
    ws = Q.w_sharp()
    f = lambda W, k: (1 + Q.W_MIN[k] ** 2 / W ** 2) / (W - Q.W_MIN[k] ** 2 / W)
    assert abs(ws.sum() - 1) < 1e-12 and np.ptp([f(ws[k], k) for k in range(3)]) < 1e-9     # Lagrange: equal partials of ℓ
    r = Q.b_job((1e-4, 2.0))
    assert r['ier'] == 1 and np.abs(r['zero'] - ws).max() < 1e-4 and np.all(r['eig'].real < 0)


def test_the_density_exponent_near_the_saturation():
    for kF in (0.0, 1.0, 2.0):
        assert abs(Q.c_job((kF, 0.003))['expo'] - (kF - 1)) < 0.03


def test_the_general_form_reduces_to_the_closed_form():
    # the correction's general form at a point of M: (τ/2)[∂_βℓ̂ − ∂_Θℓ̂ (∇²_ΘJ)⁻¹ ∇²_{Θβ}J], ℓ̂(Θ, β) = log det ∇²_ΘJ — the derivative of ℓ
    # along M through the exact fit's implicit function Θ*(β) — computed from the current point only; the closed form is its evaluation
    net, base = E.setup()
    tau = 1e-4
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau)
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3).ravel()
    J = lambda t, b: loop.J(t.reshape(3, 3), b)
    h = 1e-4

    def hess(t, b):
        I = np.eye(9) * h; H = np.zeros((9, 9))
        for i in range(9):
            for j in range(i, 9):
                H[i, j] = H[j, i] = (J(t + I[i] + I[j], b) - J(t + I[i] - I[j], b) - J(t - I[i] + I[j], b) + J(t - I[i] - I[j], b)) / (4 * h * h)
        return H
    lhat = lambda t, b: np.linalg.slogdet(hess(t, b))[1]
    H = hess(Th, beta)
    hm = 1e-4                                                                 # the mixed Hessian ∇²_{Θβ}J by central differences of ∇_ΘJ
    gradT = lambda t, b: np.array([(J(t + hm * e, b) - J(t - hm * e, b)) / (2 * hm) for e in np.eye(9)])
    Hmix = np.array([(gradT(Th, beta + hm * e) - gradT(Th, beta - hm * e)) / (2 * hm) for e in np.eye(3)]).T     # 9 × 3
    hl = 2e-3
    dl_dT = np.array([(lhat(Th + hl * e, beta) - lhat(Th - hl * e, beta)) / (2 * hl) for e in np.eye(9)])
    dl_dB = np.array([(lhat(Th, beta + hl * e) - lhat(Th, beta - hl * e)) / (2 * hl) for e in np.eye(3)])
    general = 0.5 * tau * (dl_dB - dl_dT @ np.linalg.solve(H, Hmix))
    closed = loop.fixman_correction(Th.reshape(3, 3), beta)
    assert np.abs(general / closed - 1).max() < 2e-3


def test_the_masses_form_equals_the_plug_in_on_the_manifold_and_ignores_the_amplitudes_off_it():
    # the masses form takes W_min,k = |Â_k|/a_max from the classes' contributions Â_k: on M (Â_k = W_k a_k) it is the plug-in form; at fixed
    # masses it does not see an amplitude's excursion, which the plug-in form reads as a mass at its saturation
    net, base = E.setup()
    tau = 1e-4
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3)
    A_hat = beta / 3 * base.features(Th)[0]
    assert np.allclose(A_hat, E.A_BUMP, rtol=1e-12)
    plug = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau)
    mass = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau, fixman_A=A_hat)
    assert np.abs(mass.fixman_correction(Th, beta) / plug.fixman_correction(Th, beta) - 1).max() < 1e-12
    Tx = Th.copy(); Tx[0, 0] += 8.0                                           # the first amplitude far up its plateau, the masses fixed
    assert np.array_equal(mass.fixman_correction(Tx, beta), mass.fixman_correction(Th, beta))
    assert plug.fixman_correction(Tx, beta)[0] > 20 * plug.fixman_correction(Th, beta)[0]
    W_min = np.abs(E.A_BUMP) / net.a_max                                    # the closed form's derivative, (τ/2)∂_βℓ(W)
    W = beta / 3
    assert np.allclose(mass.fixman_correction(Tx, beta), 2.0 * (tau / 3) * (1 + W_min ** 2 / W ** 2) / (W - W_min ** 2 / W), rtol=1e-12)


def test_the_tracked_contributions_stay_put_on_the_manifold():
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, W)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP.copy(), fixman_TA=1.0)
    loop.run_dithered(Th, 3 * W, 20.0, GDD.DT, 0.0, np.random.default_rng(0))
    assert np.abs(loop.fixman_A - E.A_BUMP).max() < 1e-6


def test_the_floor_caps_the_correction_and_its_exponent():
    # referee round 28 (D3aa-2): the denominator is floored at 0.01W_k — within δ_cap = W_min(1/√0.99 − 1) of the saturation the correction is
    # κ_F(τ/m)(1 + q)/(0.01W), and the exponent there is the unfloored one minus κ_F(1 − (1 − W_min²/W²)/0.01); without the floor it is not capped
    net, base = E.setup()
    tau, kF = 1e-4, 2.0
    W = np.array([0.2505, 0.37475, 0.37475])                                  # δ_1 = 0.0005 < δ_cap = 0.00126
    q = (Q.W_MIN / W) ** 2
    A_hat = Q.W_MIN * net.a_max * np.sign(E.A_BUMP)
    for floor in (0.01, 0.0):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=kF, fixman_tau=tau, fixman_A=A_hat, fixman_floor=floor)
        c = loop.fixman_correction(E.exact_rep(net, W), 3 * W)
        den = np.maximum(W * (1 - q), floor * W)
        assert np.allclose(c, kF * (tau / 3) * (1 + q) / den, rtol=1e-12)
    assert abs(Q.delta_cap(0) - 0.25 * (1 / np.sqrt(0.99) - 1)) < 1e-15 and 0.00125 < Q.delta_cap(0) < 0.00127
    d = 0.0006
    e1, e0 = Q.c_job((kF, d, 0.01))['expo'], Q.c_job((kF, d, 0.0))['expo']
    cut = (1 - (Q.W_MIN[0] / (Q.W_MIN[0] + d)) ** 2) / 0.01
    assert abs(e1 - (e0 - kF * (1 - cut))) < 0.01 and e1 < 0.1 < 0.9 < e0


def test_the_tracked_average_relaxes_at_its_time_constant():
    # the running average Â ← Â + (dt/T_A)(W a − Â) on M without dither: a perturbation decays as (1 − dt/T_A)^n
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35]); Th = E.exact_rep(net, W)
    dA = np.array([0.05, -0.03, 0.02]); TA = 5.0; dt = GDD.DT
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP + dA, fixman_TA=TA)
    loop.run_dithered(Th, 3 * W, 10 * dt, dt, 0.0, np.random.default_rng(0))
    assert np.allclose(loop.fixman_A - E.A_BUMP, dA * (1 - dt / TA) ** 10, rtol=1e-6, atol=1e-12)


def test_kappa_is_the_leaks_projected_stiffness():
    # kappa() against the Jacobian of ∂W·P·F_leak (limit_drift's F₁ with the leak as the perturbation) at w = 1/3
    w = np.ones(2) / 3; h = 1e-4; gamma, eps = 1e-2, 3e-4
    PF = lambda v: Q.drift(v, gamma, 0.0, eps=eps)[0] - Q.drift(v, gamma, 0.0)[0]
    Jl = np.array([(PF(w + h * e) - PF(w - h * e)) / (2 * h) for e in np.eye(2)]).T
    K = Q.kappa(np.ones(3) / 3, eps, gamma)
    assert np.abs(-Jl - K).max() < 1e-5 * np.abs(K).max()


def test_the_long_runs_and_the_paired_analysis_on_synthetic_records(monkeypatch):
    # summarise_long and part_G's statistics on records with known answers (the predictions stubbed out)
    rng = np.random.default_rng(4)
    n, S = 9, 41
    s = np.arange(S, dtype=float) * 1000
    base = np.zeros((n, S, 8)); base[:, :, 0] = s
    base[:, :, 1:4] = 2.0 + 0.1 * rng.standard_normal((n, S, 3))
    W = 1 / 3 + 0.01 * rng.standard_normal((n, S, 3)); W -= (W.sum(axis=2, keepdims=True) - 1) / 3
    base[:, :, 4:7] = W; base[:, :, 7] = 0.025
    sm = Q.summarise_long(base)
    h = S // 2
    assert np.allclose(sm['half'], W[:, h:].mean(axis=1)) and abs(sm['fit2'] - 0.025 ** 2) < 1e-15 and np.all(sm['below'] == 0)
    shift = np.array([0.002, -0.0015, -0.0005])
    R1 = base.copy(); R1[:, :, 4:7] += shift
    monkeypatch.setattr(Q, 'SEEDS_G', tuple(range(n)))
    monkeypatch.setattr(Q, 'drift', lambda w, g, k, floor=0.01, eps=0.0: (np.zeros(2), np.eye(2), None))
    monkeypatch.setattr(Q, 'second_order', lambda g, e, k, start=None, hH=None: dict(W0=np.zeros(3), parts={}, mean=np.full(3, float(k)), std=np.zeros(3),
                                                                                  resid=0.0, scale=1.0))
    out = {r['kF']: r for r in Q.part_G({0.0: base, 1.0: R1})}
    m0 = W[:, h:].mean(axis=1)
    assert np.allclose(out[0.0]['shift'], m0.mean(axis=0) - 1 / 3) and np.allclose(out[0.0]['se'], m0.std(axis=0, ddof=1) / np.sqrt(n))
    assert np.allclose(out[1.0]['diff'], shift, atol=1e-15) and np.all(out[1.0]['dse'] < 1e-15) and np.allclose(out[1.0]['dso'], 1.0)
    assert out[0.0]['esc'] == [] and out[1.0]['esc'] == []


def test_the_md_needs_the_drift_raw_file(monkeypatch, tmp_path):
    # referee round 143 (D3el-13): G's last table reads gauge_drift_raw.npz from the output directory; without it the script stops
    monkeypatch.setattr(Q, 'OUT', str(tmp_path))
    import pytest
    with pytest.raises(FileNotFoundError):
        Q.drift_raw_G()
    E.save(os.path.join(str(tmp_path), 'gauge_drift_raw.npz'), G=[dict(eps=1.0)])
    assert Q.drift_raw_G()[0]['eps'] == 1.0


def test_the_remainder_is_of_order_gamma_squared_tau():
    # referee round 143 (D3el-10): the corrected drift at κ_F = 1 over γ²τ settles as γ → 0
    r = [Q.a_job((np.ones(2) / 3, g)) for g in (1e-5, 3e-6)]
    big = [x['bc'][np.argmax(np.abs(x['bc']))] / x['gamma'] for x in r]
    assert abs(big[0] / big[1] - 1) < 0.02 and 1 < abs(big[0]) < 5


def test_part_I_on_synthetic_records(monkeypatch):
    # the closest approach, the fractions inside the floor's cap, with a′ < 0.01 and inside the finite-τ threshold
    monkeypatch.setattr(Q, 'a_job', lambda args: dict(w=args[0], gamma=args[1], bc=np.zeros(3), b=np.zeros(3), law=np.zeros(3)))
    Gam = Q.gamma_kk()
    net, base = E.setup()
    phi = base.features(E.exact_rep(net, np.ones(3) / 3))[4]
    assert np.allclose(Gam, (phi ** 2).mean(axis=0)) and np.all(Gam > 0)
    R = np.zeros((2, 4, 8))
    R[:, :, 4:7] = Q.W_MIN + 0.05; R[0, 0, 4] = Q.W_MIN[0] + 0.001; R[1, 2, 6] = Q.W_MIN[2] + 0.0035
    R[:, :, 1:4] = 2.0; R[1, 3, 1] = 13.0                                   # a′ = 1 − tanh²(13/4) < 0.01
    out = Q.part_I({1.0: R, 2.0: R}, None)
    r = [x for x in out['runs'] if x['kF'] == 2.0][0]
    assert np.allclose(r['dmin'], [0.001, 0.05, 0.0035]) and np.allclose(r['in_cap'], [1 / 8, 0, 0]) and np.allclose(r['in_apr'], [1 / 8, 0, 0])
    assert np.allclose(r['in_thr'], [1 / 8, 0, 1 / 8]) and abs(out['thr'][1.0][0] - np.sqrt(Q.TAU / (48 * Gam[0]))) < 1e-15
    assert '## I.' in Q.report(Ii=out)


def test_part_J_exponent_fit_and_layer_law_on_planted_densities():
    # referee round 143 (D3el-3): the weighted log-bin fit recovers a planted power law, and the layer law's own exponent tends to κ_F − 1
    lo, hi = 0.006, 0.06
    groups = Q.window_bins(lo, hi)
    assert Q.EDGES_J[groups[0][0] - 1] >= lo - 1e-15 and Q.EDGES_J[groups[-1][-1]] <= hi + 1e-15 and all(len(g) == 3 for g in groups)
    for e in (-1.0, 0.0, 1.0, 0.37):
        fine = np.zeros(len(Q.EDGES_J) + 1)
        for i in range(1, len(Q.EDGES_J)):
            a, b = Q.EDGES_J[i - 1], Q.EDGES_J[i]
            fine[i] = 1e9 * ((b ** (e + 1) - a ** (e + 1)) / (e + 1) if e != -1 else np.log(b / a))
        assert abs(Q.expo_fit(fine, groups)[0] - e) < 1e-9
    fine = np.ones(len(Q.EDGES_J) + 1)
    for kF in (0.0, 1.0, 2.0):
        p = Q.expo_fit(fine, Q.window_bins(1e-5, 1e-4), Q.layer_density(kF, 0))[0]
        assert abs(p - (kF - 1)) < 2e-3 * max(kF, 1)
        p2 = Q.expo_fit(fine, groups, Q.layer_density(kF, 0))[0]                 # at δ ~ 0.02 the layer law bends below κ_F − 1
        assert kF - 1 - 0.2 < p2 < kF - 1


def test_part_J_verdicts_and_jackknife_on_plants():
    v = Q.verdict_J(0.95, 0.04, 0.90)
    assert v['law'] == 'consistent with the layer law' and abs(v['z'] - 1.25) < 1e-12 and v['norm'] == 'above −1'
    assert Q.verdict_J(0.5, 0.1, 0.9)['law'] == 'not the layer law' and Q.verdict_J(0.5, 0.1, 0.7)['law'] == 'consistent with the layer law'
    assert Q.verdict_J(-1.3, 0.1, -1.0)['norm'] == 'below −1'
    assert Q.verdict_J(-0.9, 0.1, -1.0)['norm'] == 'undecided'
    assert Q.verdict_J(np.nan, np.nan, 0.0)['law'] == 'too few records'
    # the jackknife over seeds: identical seeds give zero spread; the pooled slope is that of the summed counts
    groups = Q.window_bins(0.006, 0.06)
    fine = np.zeros(len(Q.EDGES_J) + 1)
    for i in range(1, len(Q.EDGES_J)):
        fine[i] = 1e6 * (Q.EDGES_J[i] ** 2 - Q.EDGES_J[i - 1] ** 2) / 2
    e, se, nb, nc = Q.expo_jack(np.array([fine] * 5), groups)
    assert abs(e - 1) < 1e-9 and se < 1e-9 and nb == len(groups)


def test_part_J_record_accumulator():
    acc = dict(n=800, cap=np.full(3, 1e-3), thr=np.full(3, 3e-3), hist=np.zeros((Q.N_SEG, 3, len(Q.EDGES_J) + 1), np.int64),
               cnt=np.zeros(Q.N_SEG, np.int64), amax=np.zeros((Q.N_SEG, 3)), dmin=np.full((Q.N_SEG, 3), np.inf),
               in_cap=np.zeros((Q.N_SEG, 3), np.int64), in_thr=np.zeros((Q.N_SEG, 3), np.int64), neg=np.zeros((Q.N_SEG, 3), np.int64), dec=[])
    W = Q.W_MIN + np.array([0.002, -0.001, 0.0005]); W[2] = 1 - W[0] - W[1]
    Th = np.zeros((3, 3)); Th[:, 0] = (5.0, -7.0, 3.0)
    acc = Q.record_J(acc, 700, Th, np.log(3 * W))
    g = 700 * Q.N_SEG // 800
    d = W - Q.W_MIN
    assert acc['cnt'][g] == 1 and acc['hist'][g, 0, np.searchsorted(Q.EDGES_J, d[0], side='right')] == 1 and acc['hist'][g, 1, 0] == 1
    assert np.allclose(acc['dmin'][g], d) and np.array_equal(acc['neg'][g], [0, 1, 0]) and np.array_equal(acc['in_thr'][g], d < 3e-3)
    assert np.allclose(acc['amax'][g], [5, 7, 3]) and acc['dec'] == [] and np.array_equal(acc['in_cap'][g], d < 1e-3)
    assert acc['cnt'].sum() == 1 and acc['in_cap'].sum() == 1 and acc['neg'].sum() == 1 and d[0] < 3e-3 < d[2]
    assert np.allclose(Q.floor_J(1e-4), 0.006) and abs(Q.delta_cap(0, Q.floor_J(1e-4)) / Q.thr_J(1.0, 1e-4, Q.gamma_kk())[0] - 0.25) < 0.01
    assert abs(Q.delta_cap(0, Q.floor_J(5e-5)) / Q.thr_J(1.0, 5e-5, Q.gamma_kk())[0] - 0.25) < 0.01


def test_part_J_numbers_are_recomputed_from_the_raw():
    raw = E.load(os.path.join(Q.OUT, 'gauge_fixman_raw.npz'))
    if 'J' not in raw:
        import pytest
        pytest.skip('no Part J in the raw file')
    J, runs = raw['J'], raw['J_runs']
    Gam = Q.gamma_kk(); ds = Q.layer_J(Q.GAMMA_J)
    for r in J['rows']:
        rr = [x['acc'] for x in runs if x['tau'] == r['tau'] and x['kF'] == r['kF']]
        assert len(rr) == len(Q.SEEDS_J) and np.allclose(r['thr'], Q.thr_J(r['kF'], r['tau'], Gam))
        for k, c in enumerate(r['cls']):
            H = np.array([a['hist'][Q.BURN_J:, k].sum(axis=0) for a in rr])
            groups = Q.window_bins(2 * r['thr'][k], ds[k] / 3)
            fine = H.sum(axis=0)
            lo = np.array([Q.EDGES_J[g[0] - 1] for g in groups]); hi = np.array([Q.EDGES_J[g[-1]] for g in groups])
            cnt = np.array([fine[g].sum() for g in groups], float); ok = cnt > 0
            if ok.sum() >= 3:                                                     # an independent weighted fit (np.polyfit, weights √n)
                x = np.log(np.sqrt(lo * hi))[ok]; y = np.log(cnt[ok] / (hi - lo)[ok])
                assert abs(np.polyfit(x, y, 1, w=np.sqrt(cnt[ok]))[0] - c['e']) < 1e-9
            assert c['nc'] == int(cnt.sum())
        q4 = r['q'][3]
        ix = np.array_split(np.arange(Q.N_SEG), 4)[3]
        assert np.allclose(q4['amax'], np.mean([a['amax'][ix].max(axis=0) for a in rr], axis=0))
        nb = np.sum([a['neg'][ix].sum(axis=0) for a in rr], axis=0) / np.sum([a['cnt'][ix].sum() for a in rr])
        assert np.allclose(q4['neg'], nb) and r['conf']['ok'] == (r['kF'] > 0)
    md = open(os.path.join(Q.OUT, 'gauge_fixman.md')).read()
    assert '\n'.join(Q.report_J(J)) in md


def test_part_J_pooling_and_confinement_verdicts_on_plants():
    cls = [dict(e=1.0, se=0.5, pred=0.8), dict(e=0.5, se=0.25, pred=0.9), dict(e=np.nan, se=np.nan, pred=0.9)]
    p = Q.pool_J(cls, 2.0)
    w = np.array([4.0, 16.0])
    assert p['n'] == 2 and abs(p['e'] - (4 + 8) / 20) < 1e-12 and abs(p['se'] - 1 / np.sqrt(20)) < 1e-12
    assert abs(p['pred'] - (3.2 + 14.4) / 20) < 1e-12 and abs(p['z0'] - (0.6 - 1) * np.sqrt(20)) < 1e-12 and w.sum() == 20
    q = [dict(amax=np.array([5.0, 4.0, 3.0]), neg=np.zeros(3), cap=np.zeros(3)) for _ in range(4)]
    c = Q.confined_J(q)
    assert c['ok'] and c['grow'] == (5.0, 5.0) and c['label'].startswith('confined')
    q[2] = dict(q[2], cap=np.array([0.0, 1e-4, 0.0]))
    assert not Q.confined_J(q)['ok']                                        # inside δ_cap in some quarter: the floor binds
    q[2] = dict(q[2], cap=np.zeros(3)); q[3] = dict(q[3], neg=np.array([0.1, 0.0, 0.0]), amax=np.array([9.0, 4.0, 3.0]))
    c = Q.confined_J(q)
    assert not c['ok'] and c['beyond'] == 0.1 and c['grow'] == (5.0, 9.0) and c['label'].startswith('not confined')


def _same(a, b, path='', rtol=1e-9, atol=1e-14):
    """a and b equal as nested dicts and lists of numbers, arrays, strings and booleans (numbers to rtol, atol)."""
    if isinstance(b, dict):
        assert isinstance(a, dict) and set(a) == set(b), path
        for k in b:
            _same(a[k], b[k], f'{path}/{k}', rtol, atol)
    elif isinstance(b, (list, tuple)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f'{path}[{i}]', rtol, atol)
    elif isinstance(b, (str, bool, np.bool_)) or b is None:
        assert a == b, path
    else:
        x, y = np.asarray(a, float), np.asarray(b, float)
        assert x.shape == y.shape and np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True), path


def _raw(*keys):
    raw = E.load(os.path.join(Q.OUT, 'gauge_fixman_raw.npz'))
    if any(k not in raw for k in keys):
        import pytest
        pytest.skip(f'no {keys} in the raw file')
    return raw


def test_part_I_is_recomputed_from_the_raw():
    # Part I is stored in the raw (the md only formats it): the remainder at smaller γ, Γ_kk, the caps and thresholds, and the stored records'
    # closest approaches and fractions, recomputed from D's and E's stored records
    raw = _raw('I', 'D_rec', 'E_rec')
    I = raw['I']
    assert [(tuple(r['w']), r['gamma']) for r in I['rem']] == [(w, g) for w in Q.W_A for g in Q.GAMMAS_I]
    assert [(r['form'], r['kF']) for r in I['runs']] == [(f, k) for f in ('plug-in', 'masses') for k in Q.GAINS_I]
    _same(Q.part_I(raw['D_rec'], raw['E_rec']), I)
    md = open(os.path.join(Q.OUT, 'gauge_fixman.md')).read()
    assert '\n'.join(Q.report_I(I)) in md


def test_part_J_is_recomputed_in_full_from_the_stored_accumulators():
    # every number and label of Part J (the windows, the exponents and their jackknife errors, the layer law's, the verdicts, the halves, the
    # quarters, the confinement labels, the pooled rows, Part E's window) recomputed from the stored accumulators and Part E's records; the
    # stored accumulators' own set-up (n, the caps at the shrinking floor, the thresholds) against the constants
    raw = _raw('J', 'J_runs', 'E_rec')
    _same(Q.part_J(raw['J_runs'], raw['E_rec']), raw['J'])
    Gam = Q.gamma_kk()
    assert [(r['tau'], r['kF'], r['seed']) for r in raw['J_runs']] == [(t, k, s) for t, k in Q.GROUPS_J for s in Q.SEEDS_J]
    for r in raw['J_runs']:
        a = r['acc']
        assert a['n'] == int(round(Q.S_J[r['tau']] / Q.DT)) and np.allclose(a['thr'], Q.thr_J(r['kF'], r['tau'], Gam), rtol=1e-12)
        assert np.allclose(a['cap'], [Q.delta_cap(k, Q.FLOOR_J * np.sqrt(r['tau'] / 1e-4)) for k in range(3)], rtol=1e-12)
        assert a['cnt'].sum() == a['n'] // Q.REC_J + 1 and np.all(a['hist'].sum(axis=2) == a['cnt'][:, None])
        assert len(a['dec']) == a['n'] // Q.BLOCK_J + 1 and np.allclose(a['dec'][:, 0], Q.BLOCK_J * Q.DT * np.arange(len(a['dec'])))
    assert np.allclose(Q.layer_J(Q.GAMMA_J), 0.5 * np.sqrt(0.1 * 4.0 * np.abs(E.A_BUMP) / (3 * E.G[0])))


def test_j_runs_set_up_and_a_short_run_against_run_dithered(monkeypatch):
    # j_runs' batch at its own set-up (one seed, 2000 time units): its records equal run_dithered's at the same seed, γ, κ_F, τ, Â and
    # floor (the accumulator rebuilt here from the path), and its first rows equal the stored long runs' (the same seed and noise)
    raw = _raw('J_runs')
    monkeypatch.setattr(Q, 'SEEDS_J', (800,))
    monkeypatch.setattr(Q, 'S_J', {1e-4: 2000.0, 5e-5: 2000.0})
    monkeypatch.setattr(Q, 'CKPT', None)
    out = Q.j_runs()
    assert [(r['tau'], r['kF'], r['seed']) for r in out] == [(t, k, 800) for t, k in Q.GROUPS_J]
    net, base = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    A0 = base.features(Th0)[0] / 3
    for r in out:
        st = [x['acc'] for x in raw['J_runs'] if (x['tau'], x['kF'], x['seed']) == (r['tau'], r['kF'], 800)][0]
        a = r['acc']
        assert np.allclose(a['dec'], st['dec'][:len(a['dec'])], rtol=0, atol=1e-12) and len(a['dec']) == 3
        assert np.allclose(a['cap'], st['cap'], rtol=0, atol=0) and np.allclose(a['thr'], st['thr'], rtol=0, atol=0)
    tau, kF = 5e-5, 2.0
    fl = 0.006 * np.sqrt(tau / 1e-4)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.1, fixman=kF, fixman_tau=tau, fixman_A=A0, fixman_floor=fl)
    path = []
    loop.run_dithered(Th0, np.ones(3), 2000.0, Q.DT, tau, np.random.default_rng(800), rec_every=Q.REC_J,
                      callback=lambda s_, T, b: path.append((s_, T[:, 0].copy(), b / 3)))
    a = [x for x in out if (x['tau'], x['kF']) == (tau, kF)][0]['acc']
    n = int(round(2000.0 / Q.DT)); cap = np.array([Q.delta_cap(k, fl) for k in range(3)]); thr = Q.thr_J(kF, tau, Q.gamma_kk())
    hist = np.zeros_like(a['hist']); cnt = np.zeros(Q.N_SEG, int); amax = np.zeros((Q.N_SEG, 3)); dmin = np.full((Q.N_SEG, 3), np.inf)
    inc, inth, neg = (np.zeros((Q.N_SEG, 3), int) for _ in range(3))
    for j, (s_, al, W) in enumerate(path):
        k_ = j * Q.REC_J; g = min(k_ * Q.N_SEG // n, Q.N_SEG - 1); d = W - Q.W_MIN
        for c in range(3):
            hist[g, c, np.searchsorted(Q.EDGES_J, d[c], side='right')] += 1
        cnt[g] += 1; amax[g] = np.maximum(amax[g], np.abs(al)); dmin[g] = np.minimum(dmin[g], d)
        inc[g] += d < cap; inth[g] += d < thr; neg[g] += d < 0
    assert len(path) == n // Q.REC_J + 1 and np.array_equal(a['hist'], hist) and np.array_equal(a['cnt'], cnt)
    assert np.allclose(a['amax'], amax, rtol=0, atol=1e-10) and np.allclose(a['dmin'], dmin, rtol=0, atol=1e-12)
    assert np.array_equal(a['in_cap'], inc) and np.array_equal(a['in_thr'], inth) and np.array_equal(a['neg'], neg)
    dec = [np.r_[s_, al, W] for j, (s_, al, W) in enumerate(path) if (j * Q.REC_J) % Q.BLOCK_J == 0]
    assert np.allclose(a['dec'], dec, rtol=0, atol=1e-12)


def test_j_runs_hands_the_batch_each_runs_parameters(monkeypatch):
    # the batch j_runs builds, read before any step: per run the gain, the correction's τ, the floor 0.006(τ/10⁻⁴)^{1/2}, γ = 0.1, the masses
    # form's Â, the noise √(2τ dt G), the number of steps, the caps at that floor and the thresholds (the floor never binds in the stored
    # runs, so their records cannot tell a wrong floor at τ = 5·10⁻⁵)
    import src.core.gauge_batch as GB
    seen = {}

    def stub(bd, runs, dt, block, **kw):
        seen.update(bd=bd, runs=runs, dt=dt, block=block, kw=kw)
        return {r['key']: dict(done=False) for r in runs}
    monkeypatch.setattr(GB, 'run_batch', stub)
    monkeypatch.setattr(Q, 'CKPT', None)
    assert Q.j_runs() is None
    bd, runs = seen['bd'], seen['runs']
    assert [r['key'] for r in runs] == [(t, k, s) for t, k in Q.GROUPS_J for s in Q.SEEDS_J]
    tau = np.array([r['key'][0] for r in runs]); kF = np.array([r['key'][1] for r in runs])
    assert np.allclose(bd.floor[:, 0], 0.006 * np.sqrt(tau / 1e-4), rtol=1e-14) and np.array_equal(bd.fixman[:, 0], kF)
    assert np.array_equal(bd.tau_F[:, 0], tau) and np.all(bd.gamma == Q.GAMMA_J) and np.allclose(bd.A_hat, E.A_BUMP, rtol=1e-12)
    assert seen['dt'] == Q.DT and seen['block'] == Q.BLOCK_J and seen['kw']['rec_every'] == Q.REC_J and seen['kw']['record'] is Q.record_J
    Gam = Q.gamma_kk()
    for r, t, k in zip(runs, tau, kF):
        assert np.allclose(r['sq'], np.sqrt(2 * t * Q.DT * E.G), rtol=1e-14) and r['n'] == int(round(Q.S_J[t] / Q.DT)) and r['seed'] == r['key'][2]
        fl = 0.006 * np.sqrt(t / 1e-4)
        assert np.allclose(r['acc0']['cap'], np.abs(E.A_BUMP) / 4.0 * (1 / np.sqrt(1 - fl) - 1), rtol=1e-12)
        assert np.allclose(r['acc0']['thr'], np.sqrt(max(k, 1.0) * t / (3 * 16.0 * Gam)), rtol=1e-12)


def test_part_K_cut_density_and_pooling_on_plants():
    # referee round 150 (D3es-5, -6): the cut through each face (the masses sum to one, the other two equal, the direction keeps the sum);
    # the density from local exponents reproduces a planted power law and the layer law's own exponent; the layer law's local exponent is
    # the derivative of its log; a pooled jackknife over identical seeds has no spread
    for k in range(3):
        w, v, a = Q.cut_K(k, 0.01)
        W = np.r_[w, 1 - w.sum()]
        assert abs(W[k] - Q.W_MIN[k] - 0.01) < 1e-15 and np.allclose(np.delete(W, k), (1 - W[k]) / 2)
        assert abs(a @ v - 1) < 1e-15 and abs(a @ w + (k == 2) - W[k]) < 1e-15
    groups = Q.window_bins(0.0063, 0.05)
    counts = np.zeros(len(Q.EDGES_J) + 1); counts[groups[0][0]:groups[-1][-1] + 1] = 100.0
    for s in (-1.0, 0.0, 0.8):
        f = Q.kat_density(Q.GRID_K, np.full(len(Q.GRID_K), s))
        assert abs(Q.expo_fit(counts, groups, f)[0] - s) < 1e-3
    for kF in (0.0, 1.0, 2.0):
        ex = np.array([Q.law_local_K(kF, 0, d) for d in Q.GRID_K])
        e1 = Q.expo_fit(counts, groups, Q.kat_density(Q.GRID_K, ex))[0]; e2 = Q.expo_fit(counts, groups, Q.layer_density(kF, 0))[0]
        assert abs(e1 - e2) < 2e-3
        d, h = 0.02, 1e-6
        lg = lambda x: np.log(x ** (kF - 1) * (x + 2 * Q.W_MIN[0]) ** (kF - 1) * (Q.W_MIN[0] + x) ** (-kF))
        assert abs(Q.law_local_K(kF, 0, d) - d * (lg(d + h) - lg(d - h)) / (2 * h)) < 1e-6
    import pytest
    with pytest.raises(AssertionError):
        Q.kat_density(Q.GRID_K, np.zeros(len(Q.GRID_K)))(0.003, 0.01)


def test_part_K_is_recomputed_from_the_raw():
    # every number of Part K from Part J's stored accumulators, Part J's analysis and the stored grid; the grid's local exponents
    # recomputed at a few points; the jackknife by brute force (each seed removed from every fit, Part J's fits recomputed); the md
    raw = _raw('K', 'K_grid', 'J', 'J_runs')
    K, kg = raw['K'], raw['K_grid']
    _same(Q.part_K(raw['J_runs'], raw['J'], kg), K)
    assert np.allclose(kg['grid'], Q.GRID_K)
    for kF, k, i in ((0.0, 0, 0), (1.0, 1, 17), (2.0, 2, 35)):
        assert abs(Q.local_K(kF, k, Q.GRID_K[i]) - kg['ex'][(kF, k)][i]) < 1e-9
    assert np.allclose(kg['at'][(2.0, 0)], [Q.local_K(2.0, 0, d) for d in Q.DELTAS_K], rtol=0, atol=1e-9)
    Gam = Q.gamma_kk(); ds = Q.layer_J(Q.GAMMA_J); S = len(Q.SEEDS_J)
    for kF in (1.0, 2.0):
        fits = [(r, k) for r in raw['J']['rows'] if r['kF'] == kF for k in range(3)]
        w = np.array([1 / r['cls'][k]['se'] ** 2 for r, k in fits])
        P = []
        for i in range(S):
            e = []
            for r, k in fits:
                rr = [x['acc'] for x in raw['J_runs'] if x['tau'] == r['tau'] and x['kF'] == kF and x['seed'] != Q.SEEDS_J[i]]
                H = np.array([a['hist'][Q.BURN_J:, k].sum(axis=0) for a in rr]).sum(axis=0)
                e.append(Q.expo_fit(H, Q.window_bins(2 * Q.thr_J(kF, r['tau'], Gam)[k], ds[k] / 3))[0])
            P.append(w @ e / w.sum())
        P = np.array(P); p = [x for x in K['pooled'] if x['kF'] == kF][0]
        assert abs(np.sqrt((S - 1) / S * np.sum((P - P.mean()) ** 2)) - p['se']) < 1e-12
        assert abs(p['e'] - [x for x in raw['J']['pooled'] if x['kF'] == kF][0]['e']) < 1e-12
    assert abs(K['diff']['d'] - (K['pooled'][2]['e'] - K['pooled'][1]['e'])) < 1e-15
    for r in K['ratio']:
        assert np.allclose(r['ratio'], np.asarray(r['cap']) / np.asarray(r['thr']))
    md = open(os.path.join(Q.OUT, 'gauge_fixman.md')).read()
    assert '\n'.join(Q.report_K(K)) in md
