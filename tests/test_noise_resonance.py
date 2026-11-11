"""Tests for experiments/exp_noise_resonance.py (Proposition D‴.60): the oscillator law's transfer function, the curve's mobilities, the
averaged loop's first order against the linearised loop, its tone averages against direct time averages, the loop's chunked period
means, the thresholds' bookkeeping and the report."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_noise_resonance as NR

CT, DZ = NR.CT, NR.DZ


def test_transfer_reduces_to_the_high_pass_without_adaptation():
    """At K = 0 and f̂′ = 0 the law's Y = 1/(iω + k_e) and the measured error's iωY = iω/(iω + k_e)."""
    om = np.array([3.0, 17.0])
    o = dict(f1=np.zeros(2), K=np.zeros(2))
    Y, Ee = NR.transfer(o, om, 10.0)
    assert np.allclose(Y, 1 / (1j * om[None, :] + 10.0))
    assert np.allclose(Ee, 1j * om[None, :] / (1j * om[None, :] + 10.0))
    o = dict(f1=np.array([2.0]), K=np.array([200.0]))
    Y, _ = NR.transfer(o, np.array([5.0]), 4.0)
    assert np.isclose(Y[0, 0], (5j - 2.0) / (200.0 - 25.0 + 20j))


def test_curve_inverts_and_its_derivatives():
    """alpha_of inverts w(α); curve's μ agrees with mobilities' and dμ/dw, dv/dw with finite differences along the curve."""
    L = NR.law()
    w = np.array([0.9, -0.6, 0.05])
    al = NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0)
    for metric in NR.METRICS:
        ww, mu, dmu, v, dv = NR.curve(L, NR.EPS, al, metric)
        assert np.allclose(ww, w, rtol=1e-12, atol=1e-14)
        w0, mu0, v0 = NR.mobilities(L, NR.EPS, al, metric)
        assert np.allclose(mu, mu0) and np.allclose(v, v0)
        h = 1e-6
        wp, mup, _, vp, _ = NR.curve(L, NR.EPS, al + h, metric)
        wm, mum, _, vm, _ = NR.curve(L, NR.EPS, al - h, metric)
        assert np.allclose(dmu, (mup - mum) / (wp - wm), rtol=1e-6)
        assert np.allclose(dv, (vp - vm) / (wp - wm), rtol=1e-6)
    _, _, _, v, dv = NR.curve(L, NR.EPS, al, 'paradigm')
    assert np.allclose(v, L.G * L.m * w) and np.allclose(dv, L.G * L.m)              # v = Gβa = Gmw under the paradigm's metric


def test_jitter_filter():
    """The moving average over one period removes a linear trend and keeps a sinusoid of that period."""
    nwin = 101                                                                          # the sinusoid's period in samples (odd: centred)
    t = np.arange(4000) / nwin
    z = 0.3 + 0.02 * t + 0.1 * np.sin(2 * np.pi * t)
    assert np.allclose(NR.jitter(z, nwin), (0.1 * np.sin(2 * np.pi * t))[nwin:-nwin], atol=1e-12)
    assert np.abs(NR.jitter(0.3 + 0.02 * t, nwin)).max() < 1e-12


def _fit():
    L = NR.law()
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    return L, al0, c0, idx


def test_the_fitted_state_is_a_rest_of_the_averaged_loop_without_noise():
    L, al0, c0, idx = _fit()
    nz = dict(w=np.array([20.0]), ph=np.array([0.0]), amp=np.array([0.0]))
    for quasi in (False, True):
        rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, quasi=quasi)
        y = np.zeros(n)
        y[1:4], y[4:7] = CT.ASTAR, CT.CSTAR
        for t in (0.0, 3.0, 11.0):
            assert np.abs(rhs(t, y)).max() < 1e-12


def test_the_averaged_loops_first_order_is_the_linearised_loop():
    """At a tiny amplitude the envelopes' Re(Θc e^{iωt}) is the linearised loop's δc (`variational`) from the same start."""
    L, al0, c0, idx = _fit()
    from scipy.integrate import solve_ivp
    om, amp, ph0 = 19.0, 1e-7, 0.4
    nz = dict(w=np.array([om]), ph=np.array([ph0]), amp=np.array([amp]))
    ts = np.linspace(9.0, 10.5, 301)                                                    # x_m passes the 0.61 kernel near t = 3.3, 12.4
    for metric in NR.METRICS:
        rhs, n = NR.envelope_rhs(L, NR.EPS, metric, nz, CT.KE, CT.CSTAR, CT.ASTAR)
        y0 = np.zeros(n)
        y0[1:4], y0[4:7] = CT.ASTAR, CT.CSTAR
        sol = solve_ivp(rhs, (9.0, 10.5), y0, t_eval=ts, method='DOP853', rtol=1e-10, atol=1e-16)
        Z = sol.y[7:7 + 7] + 1j * sol.y[14:21]                                          # one tone: (E, Θw×3, Θc×3)
        dc_env = np.real(Z[4:7] * np.exp(1j * om * ts)[None, :]).T
        _, dc_var = NR.variational(L, NR.EPS, al0[idx], c0[idx], metric, ts, om, amp, ph0, CT.KE)
        assert np.abs(dc_env - dc_var).max() < 1e-7 * np.abs(dc_var).max()


def test_the_tone_averages_are_time_averages():
    """For hand-built envelopes, the blur and the correlation terms equal their direct averages over one tone period."""
    L, al0, c0, idx = _fit()
    om = 15.0
    nz = dict(w=np.array([om]), ph=np.array([0.0]), amp=np.array([1e-3]))
    rhs_all, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    rhs_no_blur, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'corr'))
    rhs_no_corr, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'blur'))
    rng = np.random.default_rng(3)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    y = np.zeros(n)
    y[1:4], y[4:7] = CT.ASTAR, CT.CSTAR
    y[7:14], y[14:21] = Z.real, Z.imag
    t = 2.0                                                                              # x_m = sin(0.4) ≈ 0.39
    blur = rhs_all(t, y)[0] - rhs_no_blur(t, y)[0]
    xm = np.sin(CT.OM * t)
    u = (xm - CT.CSTAR) / CT.SIG
    ph = np.exp(-u ** 2 / 2)
    d1, d2 = ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2
    s = np.linspace(0, 2 * np.pi, 4001)[:-1]
    rot = np.exp(1j * s)
    tw = np.real(Z[1:4, None] * rot[None, :])
    tc = np.real(Z[4:7, None] * rot[None, :])
    direct = -0.5 * np.mean((2 * d1[:, None] * tw * tc + CT.ASTAR[:, None] * d2[:, None] * tc ** 2).sum(0))
    assert np.isclose(blur, direct, rtol=1e-10, atol=1e-18)
    corr = rhs_all(t, y)[1:7] - rhs_no_corr(t, y)[1:7]
    al = NR.alpha_of(L, NR.EPS, CT.ASTAR, np.sign(CT.ASTAR) * 3.0)
    _, mu, dmu, v, dv = NR.curve(L, NR.EPS, al, 'paradigm')
    D = np.real((Z[0] + nz['amp'][0]) * rot)
    cw = np.mean((dmu[:, None] * ph[:, None] * tw + mu[:, None] * d1[:, None] * tc) * D[None, :], axis=1)
    cc = np.mean((dv[:, None] * d1[:, None] * tw + v[:, None] * d2[:, None] * tc) * D[None, :], axis=1)
    assert np.allclose(corr, np.r_[cw, cc], rtol=1e-10, atol=1e-18)


def test_quasi_static_envelopes_are_the_oscillator_law():
    """quasi=True: the measured error's envelope is A's law iωY·N and the correlation averages vanish."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([1e-3]))
    rhs_q, nq = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, quasi=True)
    rhs_nc, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'blur'), quasi=True)
    assert nq == 7
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (1.0, 2.6, 5.0):
        assert np.allclose(rhs_q(t, y)[1:], rhs_nc(t, y)[1:], atol=1e-15)


def test_sim_means_do_not_depend_on_the_chunks():
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(0.005, dict(freqs=(21.37,)))
    W1, C1 = NR.sim_means(L, 'paradigm', nz, CT.KE, 2, chunk=1, k=400)
    W2, C2 = NR.sim_means(L, 'paradigm', nz, CT.KE, 2, chunk=2, k=400)
    assert W1.shape == (2, 3) and np.allclose(W1, W2, atol=1e-9) and np.allclose(C1, C2, atol=1e-9)


def test_left_periods_and_thresholds():
    dev = np.array([[0.1, 0.0, 0.2], [0.31, 0.0, 0.29], [0.5, 0.1, 0.4]])
    assert NR.left_periods(dev, 0.3) == [2, -1, 3]
    rs = [dict(s=0.02, left=[5, -1, -1]), dict(s=0.01, left=[-1, -1, -1]), dict(s=0.03, left=[4, -1, 7])]
    assert NR.thresholds(rs) == [0.02, None, 0.03]


def test_report_on_synthetic_results():
    A = [dict(metric=mt, om=om, sqrtK=15.0 if mt == 'paradigm' else 4.8, v=[1.0, 2.0, 3.0], mu=[0.1, 0.2, 0.3],
              sim_c=[1.0, 2.0, 3.0], var_c=[1.0, 2.0, 2.0], qs_c=[1.1, 2.0, 2.0] if om >= 19 else [1.5, 2.0, 2.0],
              sim_w=[1.0, 1.0, 1.0], var_w=[1.0, 1.0, 1.0], qs_w=[1.0, 1.0, 1.0]) for mt in NR.METRICS for om in (13.9, 19.1)]
    mk = len(NR.MARKS_B)
    b = lambda mt, fr, s, terms, quasi, f: dict(metric=mt, freqs=fr, s=s, terms=list(terms), quasi=quasi, status=0,
                                                 dw_sim=[[1e-3, -2e-3, 4e-3]] * mk, dc_sim=[[1e-4, 0.0, 2e-4]] * mk,
                                                 dw_env=[[f * 1e-3, -2e-3, 4e-3]] * mk, dc_env=[[1e-4, 0.0, 2e-4]] * mk)
    mt, fr, s = NR.B2_CASE
    B = [b(mt, list(fr), s, NR.TERMS, False, 1.0), b('wfr', None, s, NR.TERMS, False, 1.0)]
    B2 = [b(mt, list(fr), s, ('ctl', 'cross', 'blur'), False, 0.5), b(mt, list(fr), s, NR.TERMS, True, 0.25)]
    C = [dict(m=12.0, ke=10.0, s=x, left=[3 if x >= 0.02 else -1, -1, 5 if x >= 0.03 else -1], sqrtK=15.76) for x in (0.01, 0.02, 0.03)]
    md = NR.report(A, B, B2, C)
    assert '[1.0000, 1.5000]' in md                                                     # the loop over the linearisation, centres
    assert '[1.000, 1.500] at ω ≥ 10, [1.000, 1.100] at ω ≥ 19' in md
    assert "D‴.49's 16 tones" in md and 'one tone 21.37' in md
    assert f"| without {NR.TERM_NAMES['corr']} | +5.00e-04, -2.00e-03, +4.00e-03 | 0.500, 1.000, 1.000 |" in md
    assert '| quasi-static envelopes | +2.50e-04, -2.00e-03, +4.00e-03 | 0.250, 1.000, 1.000 |' in md
    assert '| 12 | 10 | 15.76 | 0.02 (1/1) / 0.02 | — / — | 0.03 (1/1) / 0.03 | 3 (0) | 3 [3, 5] |' in md
    assert '\n\n\n' not in md and md.endswith('\n') and not md.endswith('\n\n')      # referee round 85 (D3cf-9): one blank line between blocks


def test_the_leaving_levels_over_seeds():
    """Referee round 56: per kernel the least level at which a seed left (how many there), the level from which every seed left at
    every level above, and '*' where a seed left at one level and stayed at a higher one."""
    row = lambda x, sd, left: dict(m=12.0, ke=10.0, s=x, opts=dict(seed=sd), left=left, sqrtK=15.76)
    rs = [row(0.01, 7, [-1, -1, -1]), row(0.01, 1, [40, -1, -1]), row(0.02, 7, [30, 50, -1]), row(0.02, 1, [20, -1, -1]),
          row(0.03, 7, [25, -1, -1]), row(0.03, 1, [10, -1, 60])]
    assert NR.c_cells(rs) == ['0.01 (1/2) / 0.02', '0.02 (1/2) / — *', '0.03 (1/2) / —']


def test_the_control_averages_are_time_averages():
    """'ctl' and 'cross' equal the direct averages over one tone period of ½f″e₁² − ½f̂″(e₁ + n)² and −∇_θf̂′·θ₁(e₁ + n), at a mean
    error ē ≠ 0 (f and f̂ at x_m + ē)."""
    L, al0, c0, idx = _fit()
    om = 15.0
    nz = dict(w=np.array([om]), ph=np.array([0.7]), amp=np.array([1e-3]))
    mk = lambda terms: NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=terms)[0]
    rhs_all, rhs_no_ctl, rhs_no_cross = mk(NR.TERMS), mk(('cross', 'blur', 'corr')), mk(('ctl', 'blur', 'corr'))
    rng = np.random.default_rng(5)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    y = np.zeros(7 + 14)
    y[0], y[1:4], y[4:7] = 0.01, CT.ASTAR * 1.01, CT.CSTAR + 0.003
    y[7:14], y[14:21] = Z.real, Z.imag
    t = 2.0
    xs = np.sin(CT.OM * t) + y[0]
    rot = np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    e1 = np.real(Z[0] * rot)
    n = np.real(nz['amp'][0] * np.exp(1j * nz['ph'][0]) * rot)
    u = (xs - CT.CSTAR) / CT.SIG
    g = np.exp(-u ** 2 / 2)
    f2 = (CT.ASTAR * g * (u ** 2 - 1) / CT.SIG ** 2).sum()
    w, c = y[1:4], y[4:7]
    v_ = (xs - c) / CT.SIG
    ph = np.exp(-v_ ** 2 / 2)
    d1, d2 = ph * v_ / CT.SIG, ph * (v_ ** 2 - 1) / CT.SIG ** 2
    fh2 = (w * d2).sum()
    ctl = np.mean(0.5 * f2 * e1 ** 2 - 0.5 * fh2 * (e1 + n) ** 2)
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_ctl(t, y)[0], ctl, rtol=1e-9, atol=1e-18)
    tw, tc = np.real(Z[1:4, None] * rot[None, :]), np.real(Z[4:7, None] * rot[None, :])
    cross = np.mean((d1[:, None] * tw + w[:, None] * d2[:, None] * tc).sum(0) * (e1 + n))      # −∇_θf̂′·θ₁ = Σ(∂_cφ θw + w∂²_cφ θc)
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_cross(t, y)[0], cross, rtol=1e-9, atol=1e-18)


def test_the_envelopes_are_driven_at_the_regressors_point():
    """With ē ≠ 0 the envelopes' drive uses the regressor at x_m, the fit's gradient at x_m + ē."""
    L, al0, c0, idx = _fit()
    om = 15.0
    N = 1e-3 * np.exp(0.2j)
    nz = dict(w=np.array([om]), ph=np.array([0.2]), amp=np.array([1e-3]))
    rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    y = np.zeros(n)
    y[0], y[1:4], y[4:7] = 0.03, CT.ASTAR, CT.CSTAR
    t = 3.0
    dZ = rhs(t, y)[7:14] + 1j * rhs(t, y)[14:21]
    xm = np.sin(CT.OM * t)
    u = (xm - CT.CSTAR) / CT.SIG
    ph = np.exp(-u ** 2 / 2)
    w, mu, v = NR.mobilities(L, NR.EPS, NR.alpha_of(L, NR.EPS, CT.ASTAR, np.sign(CT.ASTAR) * 3.0), 'paradigm')
    assert np.allclose(dZ[1:4], mu * ph * N, rtol=1e-10)
    assert np.allclose(dZ[4:7], v * ph * u / CT.SIG * N, rtol=1e-10)
    us = (xm + 0.03 - CT.CSTAR) / CT.SIG
    fh1 = -(CT.ASTAR * np.exp(-us ** 2 / 2) * us / CT.SIG).sum()
    assert np.isclose(dZ[0], -fh1 * N - CT.KE * N, rtol=1e-10)


def test_quasi_static_blur_is_the_oscillator_laws():
    """quasi=True at the fit: the blur average is −½|Y|²p·Σ_i[2μ_iv_iφ_i(∂_cφ_i)² + v_i²w_i(∂_cφ_i)²∂²_cφ_i] with A's Y."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([1e-3]))
    rq, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('blur',), quasi=True)
    r0, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=(), quasi=True)
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (2.2, 2.6, 3.4):
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.array([np.sin(CT.OM * t)]))
        Y = NR.transfer(o, om, CT.KE)[0][0, 0]
        blur = (2 * o['mu'] * o['v'] * o['ph'][0] * o['d1'][0] ** 2 + o['v'] ** 2 * o['w'] * o['d1'][0] ** 2 * o['d2'][0]).sum()
        assert np.isclose(rq(t, y)[0] - r0(t, y)[0], -0.5 * abs(Y) ** 2 * (1e-3) ** 2 / 2 * blur, rtol=1e-9, atol=1e-18)


def test_sim_means_at_the_noiseless_fit_are_the_fit():
    L, al0, c0, idx = _fit()
    W, C = NR.sim_means(L, 'paradigm', None, CT.KE, 2, chunk=1, k=200)
    assert np.allclose(W, CT.ASTAR[None, :], atol=1e-10) and np.allclose(C, CT.CSTAR[None, :], atol=1e-10)


def test_quasi_static_control_average_has_the_laws_phase():
    """quasi=True at the fit: 'ctl' is ¼f″|D − N|² − ¼f̂″|D|² with A's D = iωY·N — the phase of (iω − f̂′) matters here, not in the blur."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    N = 1e-3
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([N]))
    rq, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl',), quasi=True)
    r0, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=(), quasi=True)
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (2.2, 2.6, 3.4):
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.array([np.sin(CT.OM * t)]))
        D = NR.transfer(o, om, CT.KE)[1][0, 0] * N
        f2 = o['f2'][0]
        assert abs(o['f1'][0]) > 0.5                                                      # f̂′ ≠ 0 there: the sign is tested
        assert np.isclose(rq(t, y)[0] - r0(t, y)[0], 0.25 * f2 * abs(D - N) ** 2 - 0.25 * f2 * abs(D) ** 2, rtol=1e-9, atol=1e-18)


def test_the_second_orders_first_order_is_the_linearised_loop():
    """`second_order`'s first order is the linearised loop (`variational`, the measured error n at the start): their period means agree."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(19.0,)))
    W1, C1, W2, C2, stt = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 1, k=400)
    per = 2 * np.pi / CT.OM
    ts = per * np.arange(401) / 400
    dw, dc = NR.variational(L, NR.EPS, al0[idx], c0[idx], 'paradigm', ts, 19.0, float(nz['amp'][0]), float(nz['ph'][0]), CT.KE)
    mw, mc = np.trapezoid(dw, axis=0) / 400, np.trapezoid(dc, axis=0) / 400
    assert stt == 0 and W1.shape == (1, 3)
    assert np.abs(W1[0] - mw).max() < 1e-7 * np.abs(dw).max() and np.abs(C1[0] - mc).max() < 1e-7 * np.abs(dc).max()


def test_the_loops_even_part_is_the_second_order():
    """Under s·n and −s·n the loop's period means' half-sum less the fit is s²⟨w₂⟩ + O(s⁴) (`second_order`), its half-difference
    s⟨w₁⟩ + O(s³): at s = 0.002, one tone at 21.37, two periods, the even part within 1%; the odd part's remainder grows by ≈ 2³ from s
    to 2s (the first kernel's cubic remainder is 13% of s⟨w₁⟩ at s = 0.002)."""
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    rem = []
    for s in (0.002, 0.004):
        nz = DZ.noise_of(s, dict(freqs=(21.37,)))
        Wp, Cp = NR.sim_means(L, 'paradigm', nz, CT.KE, 2)
        Wm, Cm = NR.sim_means(L, 'paradigm', dict(nz, ph=nz['ph'] + np.pi), CT.KE, 2)
        if s == 0.002:
            even_w, even_c = (Wp + Wm) / 2 - CT.ASTAR, (Cp + Cm) / 2 - CT.CSTAR
            assert np.abs(even_w - s ** 2 * W2).max() < 0.01 * np.abs(s ** 2 * W2).max()
            assert np.abs(even_c - s ** 2 * C2).max() < 0.01 * np.abs(s ** 2 * C2).max()
        rem.append(np.abs((Wp - Wm) / 2 - s * W1).max())
    assert 6.0 < rem[1] / rem[0] < 10.0


def test_the_b3_report_on_synthetic_results():
    """B3's table: the even part over s²⟨w₂⟩ per s, the averaged loop (B) over it, the odd part and the ranges' line."""
    mk = len(NR.MARKS_B)
    w2, c2, w1 = [2.0, -1.0, 4.0], [1.0, 1.0, -1.0], [1e-3, 0.0, 2e-3]
    rows = [dict(s=s, even_w=[[s * s * 2.0 * 1.01, -s * s, s * s * 4.0]] * mk, odd_w=[[1e-6, 0.0, 2e-6]] * mk,
                 even_c=[[s * s, s * s, -s * s]] * mk, odd_c=[[0.0, 0.0, 0.0]] * mk) for s in NR.S_B3]
    B3 = [dict(metric='paradigm', freqs=[21.37], status=0, w1=[w1] * mk, c1=[[0.0] * 3] * mk, w2=[w2] * mk, c2=[c2] * mk, rows=rows)]
    s = NR.S_B3[-1]
    B = [dict(metric='paradigm', freqs=[21.37], s=s, terms=list(NR.TERMS), quasi=False, status=0, dw_sim=[[0.0] * 3] * mk,
              dc_sim=[[0.0] * 3] * mk, dw_env=[[s * s * 2.0 * 1.05, -s * s, s * s * 4.0]] * mk, dc_env=[[0.0] * 3] * mk)]
    md = '\n'.join(NR.b3_lines(B3, B))
    assert '## B3. The exact second order' in md and '| paradigm | one tone 21.37 | +2, -1, +4 |' in md
    assert '1.0100, 1.0000, 1.0000; 1.0100, 1.0000, 1.0000; 1.0100, 1.0000, 1.0000 | 1.050, 1.000, 1.000 |' in md
    assert f"[1.0000, 1.0100] at s = {s:g}" in md and 'paradigm, one tone 21.37: [1.0000, 1.0500]' in md
    assert NR.b3_lines(None) == [] and '## B3' in NR.report(B3=B3)


def test_the_envelopes_couple_through_the_mean_error():
    """Referee round 56: at ē ≠ 0 the envelopes' own terms (μ′φΘw + μ∂_cφΘc)ē and (v′∂_cφΘw + v∂²_cφΘc)ē are the adaptation law
    linearised at the mean state — central differences of μ(w)φ(x_m; c) and v(w)∂_cφ(x_m; c) along (Θw, Θc) —, the regressor at x_m."""
    L, al0, c0, idx = _fit()
    om, e = 15.0, 0.05
    nz = dict(w=np.array([om]), ph=np.array([0.0]), amp=np.array([0.0]))             # N = 0 and E = 0 below: D = 0
    rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    rng = np.random.default_rng(11)
    Z = np.r_[0.0, rng.standard_normal(6)] + 1j * np.r_[0.0, rng.standard_normal(6)]
    w0, c0_ = CT.ASTAR * 1.02, CT.CSTAR + 0.004
    t = 2.0
    out = rhs(t, np.r_[e, w0, c0_, Z.real, Z.imag])
    dZ = out[7:14] + 1j * out[14:21]
    xm = np.sin(CT.OM * t)

    def law_terms(w, c):
        _, mu, v = NR.mobilities(L, NR.EPS, NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0), 'paradigm')
        u = (xm - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        return mu * ph, v * ph * u / CT.SIG
    h = 1e-6
    for part in (np.real, np.imag):
        tw, tc = part(Z[1:4]), part(Z[4:7])
        (gp, kp), (gm, km) = law_terms(w0 + h * tw, c0_ + h * tc), law_terms(w0 - h * tw, c0_ - h * tc)
        assert np.allclose(part(dZ[1:4] + 1j * om * Z[1:4]), e * (gp - gm) / (2 * h), rtol=1e-6, atol=1e-12)
        assert np.allclose(part(dZ[4:7] + 1j * om * Z[4:7]), e * (kp - km) / (2 * h), rtol=1e-6, atol=1e-12)


def test_the_blur_and_the_correlation_at_a_mean_error():
    """Referee round 56: at ē ≠ 0 the blur is the direct tone average of −½θ₁ᵀ∇²_θf̂θ₁ at x_m + ē (where the control evaluates f̂) and
    the regressors' correlation that of ∂_θ(Γg)·θ₁(e₁ + n) with the regressor at x_m."""
    L, al0, c0, idx = _fit()
    om, e = 15.0, 0.04
    nz = dict(w=np.array([om]), ph=np.array([0.3]), amp=np.array([1e-3]))
    mk = lambda terms: NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=terms)[0]
    rhs_all, rhs_no_blur, rhs_no_corr = mk(NR.TERMS), mk(('ctl', 'cross', 'corr')), mk(('ctl', 'cross', 'blur'))
    rng = np.random.default_rng(4)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    w, c = CT.ASTAR * 0.99, CT.CSTAR - 0.002
    y = np.r_[e, w, c, Z.real, Z.imag]
    t = 2.0
    rot = np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    tw, tc = np.real(Z[1:4, None] * rot[None, :]), np.real(Z[4:7, None] * rot[None, :])

    def kern(x):
        u = (x - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        return ph, ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2
    xm = np.sin(CT.OM * t)
    _, d1s, d2s = kern(xm + e)
    blur = -0.5 * np.mean((2 * d1s[:, None] * tw * tc + w[:, None] * d2s[:, None] * tc ** 2).sum(0))
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_blur(t, y)[0], blur, rtol=1e-10, atol=1e-18)
    phr, d1r, d2r = kern(xm)
    _, mu, dmu, v, dv = NR.curve(L, NR.EPS, NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0), 'paradigm')
    D = np.real((Z[0] + nz['amp'][0] * np.exp(1j * nz['ph'][0])) * rot)
    cw = np.mean((dmu[:, None] * phr[:, None] * tw + mu[:, None] * d1r[:, None] * tc) * D[None, :], axis=1)
    cc = np.mean((dv[:, None] * d1r[:, None] * tw + v[:, None] * d2r[:, None] * tc) * D[None, :], axis=1)
    assert np.allclose(rhs_all(t, y)[1:7] - rhs_no_corr(t, y)[1:7], np.r_[cw, cc], rtol=1e-10, atol=1e-18)
    _, d1x, _ = kern(xm + e)                                                            # the regressor point matters here
    assert np.abs(d1x - d1r).max() > 0.1 * np.abs(d1r).max()


def test_b2_drops_each_average_once_and_adds_the_quasi_static_envelopes():
    jobs = NR.b2_jobs()
    assert [NR.variant_label(dict(terms=list(j[3]), quasi=j[4])) for j in jobs] == \
        [f"without {NR.TERM_NAMES[t]}" for t in NR.TERMS] + ['quasi-static envelopes']
    assert all(j[:3] == NR.B2_CASE for j in jobs) and jobs[-1][3] == NR.TERMS


def test_b3_job_on_a_short_horizon(monkeypatch):
    """Referee round 56: `b3_job`'s bookkeeping on two periods — the half-sum of the loop under ±s·n less the fit over s²⟨w₂⟩ within 1%
    at s = 0.002, the half-difference near s⟨w₁⟩, the second order's means at the marks."""
    monkeypatch.setattr(NR, 'NP_B', 2)
    monkeypatch.setattr(NR, 'MARKS_B', (1, 2))
    monkeypatch.setattr(NR, 'S_B3', (0.002, 0.004))
    r = NR.b3_job(('paradigm', (21.37,)))
    w2, w1 = np.array(r['w2']), np.array(r['w1'])
    assert w2.shape == (2, 3) and r['status'] == 0 and [q['s'] for q in r['rows']] == [0.002, 0.004]
    q0 = r['rows'][0]
    assert np.abs(np.array(q0['even_w']) - 0.002 ** 2 * w2).max() < 0.01 * np.abs(0.002 ** 2 * w2).max()
    assert np.abs(np.array(q0['odd_w']) - 0.002 * w1).max() < 0.2 * np.abs(0.002 * w1).max()
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    assert np.allclose(w2, W2) and np.allclose(r['c2'], C2)


def test_the_c_report_uses_each_laws_own_fitted_state():
    """Referee round 56: the per-kernel √K in C's caption — each law's own fitted state (√K ∝ m^{1/2} when the centres' term
    dominates: m = 3 at half m = 12's), the first kernel's equal to the orbit's √K max."""
    C = [dict(m=m, ke=10.0, s=0.01, left=[-1, -1, -1], sqrtK=0.0) for m in NR.M_C]
    md = NR.report(C=C)
    cap = md.split('at its largest: ')[1].split(':\n')[0]
    vals = {float(p.split(':')[0].split('=')[1]): [float(x) for x in p.split(':')[1].split(',')] for p in cap.split('; ')}
    assert np.allclose(np.array(vals[3.0]) * 2, vals[12.0], rtol=1e-2)
    for m in NR.M_C:
        L = NR.law(m)
        al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.sin(np.linspace(0, 2 * np.pi, 2001)))
        assert abs(vals[m][0] - np.sqrt(o['K'].max())) < 0.02


def test_the_monodromy_and_the_floquet_means_are_the_linearised_loops():
    """Referee round 63: the Floquet recursion's period means over three periods equal the period means of `second_order`'s first order
    from the same start (in the first period to the trapezoid's error on the start's fast transient, later to 10⁻⁷), and the unforced
    multipliers lie inside the unit circle with the three slow coefficient modes first."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.4,)))
    a, ph = float(nz['amp'][0]), float(nz['ph'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.4, ph)
    Z = NR.floquet_means(M, P, a * G, a * H, 21.4, 3)
    W1, C1, _, _, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 3)
    for X, Y in ((Z[:, 1:4], W1), (Z[:, 4:7], C1)):
        assert np.abs(X[1:] - Y[1:]).max() < 1e-7 * np.abs(Y).max() and np.abs(X[0] - Y[0]).max() < 1e-3 * np.abs(Y).max()
    rho = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    assert 0.99 < rho[2] <= rho[0] < 0.995 and rho[3] < 0.3


def test_f_job_and_its_report(monkeypatch):
    """F's job at m = 12 (two B4 tones, one commensurate tone, 20 periods): the multipliers and time constants, the Floquet sums against
    a synthetic B4 whose odd part is the sums plus a known error, the saturation at the commensurate tone equal to the recursion's
    limit."""
    monkeypatch.setattr(NR, 'FREQS_B4', ((21.4,), None, (21.3,)))
    monkeypatch.setattr(NR, 'OM_COMM_F', (21.4,))
    monkeypatch.setattr(NR, 'NP_B4', 20)
    r = NR.f_job(('paradigm', 12.0))
    assert len(r['rho']) == 7 and 120 < min(r['tau'][:3]) and max(r['tau'][:3]) < 180 and set(r['mode'][:3]) == {'w1', 'w2', 'w3'}
    fz = np.array(r['sums']['21.4']['w'])
    assert fz.shape == (20, 3)
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.4,)))
    a = float(nz['amp'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.4, float(nz['ph'][0]))
    far = NR.floquet_means(M, P, a * G, a * H, 21.4, 4000)[-1]
    assert np.allclose(r['comm']['21.4']['sat'], far, rtol=1e-3, atol=1e-6)
    B4 = [dict(freqs=[21.4], s=s, odd_w_all=(fz * (1 + 0.01 * (s / 0.00125) ** 2)).tolist()) for s in (0.00125, 0.0025)]
    md = '\n'.join(NR.f_lines([r], B4))
    assert '21.4: 1.00e-02 (s = 0.00125), 4.00e-02 (s = 0.0025), ×4.00' in md and '| 21.4 | 107 |' in md


def test_b4_job_on_two_periods(monkeypatch):
    """Referee round 63: `b4_job`'s bookkeeping on two periods — the even part per s² near ⟨w₂⟩ (`second_order`), the odd part per s
    near ⟨w₁⟩, the full arrays' marks equal to the marked rows, the centre's per-period odd part."""
    monkeypatch.setattr(NR, 'NP_B4', 2)
    monkeypatch.setattr(NR, 'MARKS_B4', (1, 2))
    r = NR.b4_job(((21.37,), 0.002))
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    assert np.abs(np.array(r['even_w']) - W2).max() < 0.02 * np.abs(W2).max()
    assert np.abs(np.array(r['odd_w']) - W1).max() < 0.25 * np.abs(W1).max()                 # 21.37's odd part: its O(s³) is large
    assert np.abs(np.array(r['even_c']) - C2).max() < 0.1 * np.abs(C2).max()                  # the centres' even part (round 63)
    assert np.allclose(r['even_w_all'], r['even_w']) and np.allclose(r['odd_w_all'], r['odd_w']) and len(r['odd_c_per']) == 2


def test_c_noise_and_c_job_pass_their_options(monkeypatch):
    """Referee round 63: `c_noise` and `c_job` hand the seed, the band, the tones, the phases, the metric and k_e to the loop."""
    o = NR.c_noise(0.02, dict(freqs=(18.9, 22.5), ph=(0.1, 0.2)))
    assert np.allclose(o['w'], [18.9, 22.5]) and np.allclose(o['ph'], [0.1, 0.2]) and np.allclose(o['amp'], 0.02)
    a, b = NR.c_noise(0.02, dict(seed=1)), NR.c_noise(0.02, dict(seed=2))
    assert len(a['w']) == 16 and not np.allclose(a['w'], b['w']) and 10 <= a['w'].min() and a['w'].max() <= 30
    c = NR.c_noise(0.02, dict(band=(25.0, 75.0), seed=3))
    assert 25 <= c['w'].min() and c['w'].max() <= 75
    d1, d2 = NR.c_noise(0.04, dict(freqs=(13.1,), seed=1)), NR.c_noise(0.04, dict(freqs=(13.1,), seed=2))
    assert d1['ph'][0] != d2['ph'][0] and np.isclose(d1['amp'][0], 0.04 * np.sqrt(2))
    seen = {}

    def fake(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None):
        seen.update(m=L.m, noise=noise, metric=metric, ke=ke, n=len(ts))
        n = len(ts)
        return dict(c=np.tile(c0, (n, 1)), w=np.zeros((n, c0.size)), al=np.zeros((n, c0.size)), x=np.zeros(n))
    monkeypatch.setattr(NR, 'simulate', fake)
    r = NR.c_job((6.0, 5.0, 0.03, dict(metric='wfr', freqs=(18.9, 22.4), seed=3)))
    assert seen['metric'] == 'wfr' and seen['ke'] == 5.0 and np.isclose(seen['m'], 6.0) and np.allclose(seen['noise']['w'], [18.9, 22.4])
    assert np.allclose(seen['noise']['ph'], NR.c_noise(0.03, dict(freqs=(18.9, 22.4), seed=3))['ph'])
    assert r['left'] == [-1, -1, -1] and r['w_left'] == [None] * 3 and r['opts']['freqs'] == [18.9, 22.4]


def test_c2_cases_and_their_report():
    """Referee round 63: C2's cases (four seeds per case, the draw's own pair with its phases), their grouping over the seeds and the
    jumps counted apart from the leaves."""
    cases = NR.c2_cases()
    labs = {c[0] for c in cases}
    assert labs == {'one tone', 'tone pair', "the draw's pair", 'band 25–75', 'm', 'WFR'}
    fr, ph = NR.draw_pair()
    assert abs(fr[0] - 18.9015) < 1e-3 and abs(fr[1] - 22.5019) < 1e-3 and abs((fr[1] - fr[0]) / CT.OM - 18.0019) < 1e-3
    assert sum(c[0] == "the draw's pair" for c in cases) == 4 and all('ph' in c[4] for c in cases if c[0] == "the draw's pair")
    pair = [c for c in cases if c[0] == 'tone pair' and c[4]['freqs'] == (18.9, 22.5) and c[3] == 0.02]
    assert sorted(c[4]['seed'] for c in pair) == sorted(NR.SEEDS_C)
    rs = [dict(label='tone pair', m=12.0, ke=10.0, s=0.08, opts=dict(freqs=[18.9, 22.5], seed=sd), left=lf, maxdev=md)
          for sd, lf, md in ((7, [124, -1, -1], [0.95, 0.1, 0.1]), (1, [69, -1, -1], [0.36, 0.1, 0.1]), (2, [-1, -1, -1], [0.2, 0.1, 0.1]))]
    md = '\n'.join(NR.c2_lines(rs))
    assert '| tone pair | paradigm | 12 | 10 | tones 18.9/22.5 (difference 18ω) | 0.08 | 3 | 2 (1), 0 (0), 0 (0) | 69–124 | 0.950, 0.100, 0.100 |' in md


def test_leave_stats_and_the_c3_report():
    """Referee round 63: at a leave the coefficient there and the mean offset over the periods leave − 30 … leave − 6; C3's rows."""
    Wp = np.tile(np.array([0.9, -0.6, 0.75]), (60, 1)) + np.linspace(0, 0.1, 60)[:, None] * np.array([1, 0, 0])[None, :]
    Cp = np.tile(CT.CSTAR, (60, 1))
    Cp[:, 0] += np.r_[np.full(40, -0.01), np.full(20, -0.9)]
    left = NR.left_periods(np.abs(Cp - CT.CSTAR), NR.LEFT_C)
    w_left, pre = NR.leave_stats(Wp, Cp, left)
    assert left == [41, -1, -1] and np.isclose(w_left[0], Wp[40, 0]) and np.isclose(pre[0], -0.01) and pre[1] is None
    md = '\n'.join(NR.c3_lines([dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=7), left=left, maxdev=[0.9, 0.0, 0.0], w_left=w_left,
                                      pre=pre, W=Wp.tolist(), C=Cp.tolist())]))
    assert "| D‴.49's 16 tones | 0.02 | 7 | 1 | 41 | 0.900 → 0.968 | -0.010 | 0.90 | -1.42 |" in md


def test_b5_splits_the_second_order_by_the_phase(monkeypatch):
    """Referee round 63: the second order at φ and φ + π/2 — the constant part (their half-sum) does not depend on φ, the 2ω part (their
    half-difference) turns by 2Δφ: at φ + π/4 the half-difference of the shifted pair is the first pair's turned by π/2."""
    monkeypatch.setattr(NR, 'NP_B', 2)
    monkeypatch.setattr(NR, 'MARKS_B', (1, 2))
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.37,)))
    run = lambda dph: NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', dict(nz, ph=np.asarray(nz['ph']) + dph), CT.KE, 2)[2]
    a, b, c, d = run(0.0), run(np.pi / 2), run(np.pi / 4), run(3 * np.pi / 4)
    assert np.allclose((a + b) / 2, (c + d) / 2, rtol=1e-6, atol=1e-9 * np.abs(a).max())
    r = NR.b5_job(('paradigm', 21.37))
    assert np.allclose(r['w2'], a) and np.allclose(r['w2_shift'], b) and len(r['w_env']) == 2
    md = '\n'.join(NR.b5_lines([r]))
    assert '| paradigm | 21.37 | 213.7 |' in md


def test_a2_rows_split_the_roll_off_from_the_resonance(monkeypatch):
    """Referee round 63: at K = 0 the jitter is the roll-off alone and the resonance factor at the smallest law nearer one; the peaks
    within 20% of √K_i."""
    monkeypatch.setattr(NR, 'M_A2', (3.0, 12.0))
    monkeypatch.setattr(NR, 'OM_A2', np.linspace(4.0, 40.0, 145))
    rows = NR.a2_rows()
    assert len(rows) == 4
    f = {(r['m'], tuple(r['band'])): np.array(r['jit']) / np.array(r['jit0']) for r in rows}
    assert np.all(f[(12.0, (10.0, 30.0))] > f[(3.0, (10.0, 30.0))]) and np.all(f[(3.0, (10.0, 30.0))] > 1)
    for r in rows:
        if r['m'] == 12.0:
            assert np.all(np.abs(np.array(r['peak']) / np.array(r['sqrtK']) - 1) < 0.2)
    md = '\n'.join(NR.a2_lines(rows))
    assert 'of which the roll-off' in md


def test_the_floquet_means_at_a_rotating_tone():
    """Referee round 69 (MO3, MO4, FM1–FM3): at 21.37 rad/s, where e^{iωT} is not ±1, the Floquet recursion's period means over three
    periods equal `second_order`'s first order (the rotation, Im G and Im H all enter) — to 2·10⁻³ of the largest: `second_order`'s
    trapezoid (2000 outputs a period) is exact for tones commensurate with the period (10⁻⁸ at 21.4) and errs by ~10⁻⁴ for the others,
    where a conjugated or dropped rotation errs by O(1)."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.37,)))
    a, ph = float(nz['amp'][0]), float(nz['ph'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.37, ph)
    Z = NR.floquet_means(M, P, a * G, a * H, 21.37, 3)
    W1, C1, _, _, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 3)
    for X, Y in ((Z[:, 1:4], W1), (Z[:, 4:7], C1)):
        assert np.abs(X - Y).max() < 2e-3 * np.abs(Y).max()
    assert abs(np.exp(1j * 21.37 * 2 * np.pi / CT.OM).real) < 0.99


def test_leave_stats_window_with_a_varying_offset():
    """Referee round 69 (LS1): the mean offset before a leave is over the periods leave − 30 … leave − 6 exactly."""
    n = 80
    Wp = np.tile(CT.ASTAR, (n, 1))
    Cp = np.tile(CT.CSTAR, (n, 1))
    Cp[:, 0] += np.r_[-0.001 * np.arange(50), np.full(n - 50, -0.9)]
    left = NR.left_periods(np.abs(Cp - CT.CSTAR), NR.LEFT_C)
    _, pre = NR.leave_stats(Wp, Cp, left)
    p = left[0]
    assert p == 51 and np.isclose(pre[0], np.mean(-0.001 * np.arange(p - 1 - 30, p - 1 - 5)))


def test_c_noise_negated_and_the_raised_start():
    """Referee rounds 69, 76: `c_noise`'s neg moves every phase by π (bands, seeds and tones alike); `raised` moves one formed kernel's
    coefficient by Δ in magnitude on its curve, keeping its sign, and nothing else."""
    for o in (dict(seed=3), dict(band=(25.0, 75.0), seed=2), dict(freqs=(13.1,), seed=7)):
        a, b = NR.c_noise(0.02, o), NR.c_noise(0.02, dict(o, neg=True))
        assert np.allclose(b['ph'], np.asarray(a['ph']) + np.pi) and np.allclose(a['w'], b['w']) and np.allclose(a['amp'], b['amp'])
    L, al0, c0, idx = _fit()
    for k in (0, 1, 2):
        al = NR.raised(L, al0, idx, k, 0.05)
        w0, w1 = NR.mobilities(L, NR.EPS, al0, 'paradigm')[0], NR.mobilities(L, NR.EPS, al, 'paradigm')[0]
        assert np.isclose(abs(w1[idx[k]]) - abs(w0[idx[k]]), 0.05, rtol=0, atol=1e-12) and np.sign(w1[idx[k]]) == np.sign(CT.ASTAR[k])
        assert np.array_equal(np.delete(al, idx[k]), np.delete(al0, idx[k]))
    al = NR.raised(L, al0, idx, 1, -0.05)
    assert np.isclose(abs(NR.mobilities(L, NR.EPS, al, 'paradigm')[0][idx[1]]), 0.55, rtol=0, atol=1e-9)


def test_the_basin_edge_brackets(monkeypatch):
    """Referee rounds 69, 76: without noise kernel 1 returns from a raise of 0.07 and is ejected from 0.075 (`edge_run`); with trace the
    anatomy on TRACE_K outputs a period — at 0.075 the centre stays 0.05 off its feature from t/T ≈ 0.594 on (non-returning), x_m below
    0, and ends below −1; at 0.07 it passes 0.05 off on both sides and comes back —; kernel 3 raised detaches at t/T ≈ 0.40, x_m above
    0; the reference's frequency reaches the loop (at ω_ref = 0.4 a raise of 0.075 returns, at 0.1 one of 0.05 is ejected); `edge_job`
    brackets the edge by the grid and bisections, every run on the bracket's side its outcome says, and records ω_ref and the anatomy at
    hi and at lo; EDGE_CASES hold the rows outside the fitted window."""
    p0, x0 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.07)
    p1, x1 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075)
    assert p0 == -1 and x0 < 0.01 and 1 <= p1 <= 5 and x1 > NR.JUMP_C
    t1 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075, NP=3, trace=True)
    t0 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.07, NP=3, trace=True)
    assert 0.55 < t1[2] < 0.65 and t1[3] < 0 and t1[5] < -1.0 and t1[6] < -0.5 and len(t1) == 8
    assert t0[2] is None and t0[3] is None and t0[6] < -0.05 and t0[7] > 0.04 and abs(t0[5] - CT.CSTAR[0]) < 0.01
    t3 = NR.edge_run(12.0, 10.0, 'paradigm', 2, 0.075, NP=3, trace=True)
    assert 0.35 < t3[2] < 0.45 and t3[3] > 0 and t3[0] >= 0
    assert NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075, om=0.4)[0] == -1 and NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05, om=0.1)[0] >= 0
    assert NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05)[0] == -1
    monkeypatch.setattr(NR, 'EDGE_GRID', (0.05, 0.1))
    monkeypatch.setattr(NR, 'BISECT_EDGE', 2)
    e = NR.edge_job((12.0, 10.0, 'paradigm', 0))
    assert np.isclose(e['lo'], 0.0625) and np.isclose(e['hi'], 0.075) and e['om'] == CT.OM
    for x in e['runs']:
        assert (x['left'] >= 0) == (x['d'] >= e['hi'])
    assert e['hi_trace']['t_detach'] is not None and e['lo_trace']['t_detach'] is None and e['t_detach'] == e['hi_trace']['t_detach']
    assert np.isclose(e['hi_trace']['t_detach'], t1[2]) and np.isclose(e['hi_trace']['c_end'], t1[5])
    monkeypatch.setattr(NR, 'EDGE_GRID', (0.1, 0.15))
    monkeypatch.setattr(NR, 'BISECT_EDGE', 0)
    e = NR.edge_job((12.0, 10.0, 'paradigm', 0, 0.4))
    assert e['om'] == 0.4 and e['lo'] == 0.1 and e['hi'] == 0.15
    for c in [(12.0, 2.5, 'paradigm', 0), (12.0, 80.0, 'paradigm', 0), (1.5, 10.0, 'paradigm', 0), (24.0, 10.0, 'paradigm', 0),
              (12.0, 20.0, 'wfr', 0), (12.0, 40.0, 'wfr', 0), (12.0, 10.0, 'paradigm', 0, 0.1), (3.0, 10.0, 'paradigm', 0, 0.4)]:
        assert c in NR.EDGE_CASES


def test_c3_job_with_a_fake_loop_and_a_raised_start(monkeypatch):
    """Referee rounds 69, 76 (C3's mutants): `c3_job` takes the period means of 8 outputs a period of the formed kernels — the coefficients
    as W, the centres as C —, the leave, the largest distance, the coefficient at the leave and the mean offset before it; opts rise starts
    the loop from `raised`."""
    seen = {}
    NPc, k = NR.NP_C, 8

    def fake(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None, om=None):
        seen.update(al0=np.array(al0), n=len(ts))
        n = len(ts)
        c = np.tile(c0, (n, 1)).astype(float)
        w = np.zeros((n, c0.size))
        _, _, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
        w[:, idx] = CT.ASTAR
        per = np.minimum(np.arange(n) // k, NPc - 1)
        c[:, idx[0]] += np.where(per >= 100, 0.5, 0.0)
        w[:, idx[0]] += 0.001 * per
        return dict(c=c, w=w, al=np.zeros((n, c0.size)), x=np.zeros(n))
    monkeypatch.setattr(NR, 'simulate', fake)
    r = NR.c3_job((12.0, 10.0, 0.02, dict(seed=7, rise=[0, 0.04])))
    L, al0, c0, idx = _fit()
    assert np.allclose(seen['al0'], NR.raised(NR.law(12.0), al0, idx, 0, 0.04)) and seen['n'] == NPc * k + 1
    assert r['left'] == [101, -1, -1] and np.isclose(r['maxdev'][0], 0.5) and np.isclose(r['w_left'][0], CT.ASTAR[0] + 0.1)
    assert np.isclose(r['pre'][0], 0.0) and np.allclose(np.asarray(r['W'])[:, 0], CT.ASTAR[0] + 0.001 * np.arange(NPc))
    C = np.asarray(r['C'])
    assert np.allclose(C[:100, 0], CT.CSTAR[0]) and np.allclose(C[100:, 0], CT.CSTAR[0] + 0.5) and np.allclose(C[:, 1:], CT.CSTAR[1:])


def test_w2_job_is_the_even_part(monkeypatch):
    """Referee round 69: `w2_job` is ((W₊ + W₋)/2 − A)/s² of the loop under ±S_W2·n, the −n run's phases moved by π."""
    calls = []

    def fake(L, metric, noise, ke, NP, chunk=10, k=2000):
        calls.append((metric, ke, NP, np.asarray(noise['ph'], float).copy()))
        sgn = 1.0 if len(calls) == 1 else -1.0
        W = np.tile(CT.ASTAR, (NP, 1)) + sgn * 0.01 + 2e-6 * np.arange(1, NP + 1)[:, None]
        return W, np.tile(CT.CSTAR, (NP, 1))
    monkeypatch.setattr(NR, 'sim_means', fake)
    monkeypatch.setattr(NR, 'NP_W2', 3)
    r = NR.w2_job((3.0, 5.0, 'wfr', dict(seed=2)))
    assert [c[:3] for c in calls] == [('wfr', 5.0, 3)] * 2 and np.allclose(calls[1][3], calls[0][3] + np.pi)
    assert np.allclose(np.asarray(r['even_w_all']), 2e-6 * np.arange(1, 4)[:, None] / NR.S_W2 ** 2 * np.ones(3))


def test_the_like_for_like_the_edge_and_the_raised_reports():
    """Referee rounds 69, 76: C4's row — the rise at the leave, the twin's half-sum (not where the twin left by the leave) and the draw's
    own s²⟨w₂⟩ — and its correlation (signed); EDGE's bracket mapped to the level s* = (edge/⟨w₂⟩)^{1/2} against the observed first jump,
    over the testable cases, with the permutation null; and C6's rows, on synthetic results."""
    n = 300
    W = np.tile(CT.ASTAR, (n, 1)); W[:, 0] += np.linspace(0, 0.09, n)
    Wm = np.tile(CT.ASTAR, (n, 1)); Wm[:, 0] += np.linspace(0, 0.07, n)
    C3 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5), left=[150, -1, -1], maxdev=[0.9, 0, 0], W=W.tolist())]
    C4 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5, neg=True), left=[-1, -1, -1], maxdev=[0.1, 0, 0], W=Wm.tolist())]
    ew = np.zeros((n, 3)); ew[:, 0] = np.linspace(100.0, 200.0, n)
    W2 = [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=5), even_w_all=ew.tolist())]
    md = '\n'.join(NR.c4_lines(C3, C4, W2))
    rise, half, own = W[149, 0] - CT.ASTAR[0], (W[149, 0] + Wm[149, 0]) / 2 - CT.ASTAR[0], 0.02 ** 2 * ew[149, 0]
    assert (f"| 0.02 | 5 | 150 | {rise:+.4f} | -1 | {half:+.4f} ({half / rise:.3f}) | {own:+.4f} ({own / (0.02 ** 2 * ew[-1, 0]):.2f}) | "
            f"{rise / own:.3f} |") in md
    C4p = [dict(C4[0], left=[150, -1, -1])]
    assert f"| 0.02 | 5 | 150 | {rise:+.4f} | 150 | — |" in '\n'.join(NR.c4_lines(C3, C4p, W2))
    C3c, W2c = [], []
    for sd, p, e2 in ((1, 100, 300.0), (2, 140, 200.0), (3, 180, 100.0)):
        C3c.append(dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=sd), left=[p, -1, -1], maxdev=[0.9, 0, 0], W=W.tolist()))
        W2c.append(dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=np.full((n, 3), e2).tolist()))
    rs = np.array([W[p - 1, 0] - CT.ASTAR[0] for p in (100, 140, 180)])
    cc = np.corrcoef(rs, 0.02 ** 2 * np.array([300.0, 200.0, 100.0]))[0, 1]
    assert cc < -0.99 and f"across the 3 draws {cc:.2f}" in '\n'.join(NR.c4_lines(C3c, C4, W2c))
    E_ = [dict(m=12.0, ke=10.0, metric='paradigm', k=0, lo=0.07, hi=0.075, runs=[dict(d=0.075, left=2, maxdev=0.9)])]
    W2b = [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=np.full((n, 3), 200.0).tolist()) for sd in NR.W2_SEEDS]
    C = [dict(m=12.0, ke=10.0, s=s_, opts=dict(seed=sd), left=[5 if s_ >= 0.02 else -1, -1, -1], maxdev=[0.9 if s_ >= 0.02 else 0, 0, 0])
         for sd in NR.W2_SEEDS for s_ in (0.015, 0.02)]
    md = '\n'.join(NR.edge_lines(E_, W2b, C, None))
    a, b = np.sqrt(0.07 / 200.0), np.sqrt(0.075 / 200.0)
    assert f"| 12 | 10 | paradigm | 0.2 | 1 | D‴.49's band | (0.07, 0.075] | 2 | — | ({a:.4f}, {b:.4f}], 0.02 |" in md
    assert ("consistent in 4 of the 4 testable draw-configurations (those that jumped, and those that did not where a level above s* "
            "was run), later in 0, earlier in 0; of those that never jumped, s* lay below the largest level run in 0. With the draws' ⟨w₂⟩ "
            "permuted within each configuration (200 permutations) the consistent count is 4.00 on average (4–4; at least 4 in 200)") in md
    W2c = W2b + [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(band=(25.0, 75.0), seed=7), even_w_all=np.full((n, 3), 400.0).tolist())]
    C2 = [dict(m=12.0, ke=10.0, s=s_, opts=dict(band=[25.0, 75.0], seed=7), left=[5 if s_ >= 0.06 else -1, -1, -1],
               maxdev=[0.9 if s_ >= 0.06 else 0, 0, 0]) for s_ in (0.04, 0.06)]
    md = '\n'.join(NR.edge_lines(E_, W2c, C, C2))
    a2, b2 = np.sqrt(0.07 / 400.0), np.sqrt(0.075 / 400.0)
    assert f"| 12 | 10 | paradigm | 0.2 | 1 | 25–75 | (0.07, 0.075] | 2 | — | ({a2:.4f}, {b2:.4f}], 0.06 | —, — |" in md
    assert ("consistent in 4 of the 5 testable draw-configurations (those that jumped, and those that did not where a level above s* was "
            f"run), later in 1 (m = 12, k_e = 10, 25–75 (1; {0.06 / (0.5 * (a2 + b2)):.2f}–{0.06 / (0.5 * (a2 + b2)):.2f} times s*)), "
            "earlier in 0") in md
    C6 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5, rise=[0, 0.05]), left=[40, -1, -1], maxdev=[0.95, 0, 0])]
    assert "| 16 tones, draw 5 | 0.02 | 0.05 | 40 (0.950) | 150 |" in '\n'.join(NR.c6_lines(C6, C3))


def test_f_job_sums_and_saturations_at_a_rotating_tone(monkeypatch):
    """Referee round 69 (FJ4–FJ6): `f_job`'s Floquet sums for a tone are `floquet_means` with the tone's phase and amplitude (rms 1) over
    NP_B4 periods; at a commensurate tone its first three periods and the saturation P(I − M)⁻¹Re(aG) + Re(aH)."""
    monkeypatch.setattr(NR, 'FREQS_B4', ((21.37,),))
    monkeypatch.setattr(NR, 'NP_B4', 5)
    monkeypatch.setattr(NR, 'OM_COMM_F', (21.4,))
    r = NR.f_job(('paradigm', 12.0))
    L, al0, c0, idx = _fit()
    M, P = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE)
    for om, NP in ((21.37, 5), (21.4, 3)):
        nz = DZ.noise_of(1.0, dict(freqs=(om,)))
        a = float(nz['amp'][0])
        _, _, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, om, float(nz['ph'][0]))
        Z = NR.floquet_means(M, P, a * G, a * H, om, NP)
        if om == 21.37:
            assert np.allclose(r['sums']['21.37']['w'], Z[:, 1:4], rtol=1e-10, atol=0) and np.allclose(r['sums']['21.37']['c'], Z[:, 4:7])
        else:
            assert np.allclose(r['comm']['21.4']['first'], Z, rtol=1e-10, atol=0) and len(r['comm']['21.4']['first']) == 3
            sat = P @ np.linalg.solve(np.eye(M.shape[0]) - M, np.real(a * G)) + np.real(a * H)
            assert np.allclose(r['comm']['21.4']['sat'], sat, rtol=1e-10, atol=0)


def test_b5_job_phases_and_the_averaged_loop_per_s2(monkeypatch):
    """Referee round 69 (B5b–B5d): `b5_job` runs the second order at the tone's phase and at φ + π/2 and stores both runs' w and c, and the
    averaged loop at s = S_B[−1] and at S_B5 per s²; the report prints the small level's ratios."""
    calls, env = [], []
    nP = len(NR.MARKS_B)

    def fake_so(L, eps, al, c, metric, noise, ke, NP, k=2000):
        calls.append(float(np.asarray(noise['ph'])[0]))
        q = len(calls)
        return None, None, np.full((NP, 3), float(q)), np.full((NP, 3), 10.0 * q), 0

    def fake_env(L, eps, metric, noise, ke, NP, k=400, terms=NR.TERMS, quasi=False):
        env.append(float(np.asarray(noise['amp'])[0]))
        return np.tile(CT.ASTAR, (NP, 1)) + 0.5, np.tile(CT.CSTAR, (NP, 1)) - 0.25, 0
    monkeypatch.setattr(NR, 'second_order', fake_so)
    monkeypatch.setattr(NR, 'envelope_means', fake_env)
    r = NR.b5_job(('paradigm', 13.1))
    ph = float(DZ.noise_of(1.0, dict(freqs=(13.1,)))['ph'][0])
    assert np.isclose(calls[0], ph) and np.isclose(calls[1], ph + np.pi / 2)
    s = NR.S_B[-1]
    assert np.isclose(env[0], s * np.sqrt(2.0))
    assert np.allclose(r['w2'], 1.0) and np.allclose(r['w2_shift'], 2.0) and np.allclose(r['c2'], 10.0) and np.allclose(r['c2_shift'], 20.0)
    assert np.allclose(r['w_env'], 0.5 / s ** 2) and np.allclose(r['c_env'], -0.25 / s ** 2) and len(r['w2']) == nP
    assert len(env) == 2 and np.isclose(env[1], NR.S_B5 * np.sqrt(2.0))
    assert np.allclose(r['w_env_small'], 0.5 / NR.S_B5 ** 2) and np.allclose(r['c_env_small'], -0.25 / NR.S_B5 ** 2)
    md = '\n'.join(NR.b5_lines([r]))
    q = 0.5 / NR.S_B5 ** 2 / 1.5
    assert f"13.1: {q:.4f}–{q:.4f} (over the exact second order {0.5 / NR.S_B5 ** 2:.4f}–{0.5 / NR.S_B5 ** 2:.4f})" in md


def test_the_sensitivity_and_bodes_integral():
    """Referee round 69 (D3bp-4): A2's factor per tone and orbit point is the adaptation loop's sensitivity — the quasi-static transfer
    with K over that at K = 0, the plant's slope cancelling —; |S|² − 1 = K(2ω² − K)/|K − ω² + ik_eω|², so |S| > 1 exactly above
    √(K/2); |S(√K)| = (1 + K/k_e²)^{1/2}; Bode's integral ∫₀^∞ ln|S(iω)|dω = 0; the report prints them."""
    K, ke = 100.0, 10.0
    om = np.array([5.0, 7.0, 7.08, 10.0, 15.0, 30.0])
    S = NR.sens(K, ke, om)
    assert np.allclose(S ** 2 - 1, K * (2 * om ** 2 - K) / np.abs(K - om ** 2 + 1j * ke * om) ** 2, rtol=1e-12, atol=1e-14)
    assert np.array_equal(S > 1, om > np.sqrt(K / 2))
    assert np.isclose(NR.sens(K, ke, np.sqrt(K)), np.sqrt(1 + K / ke ** 2), rtol=1e-14)
    assert np.allclose(NR.sens(K, ke, [10.0, 15.0, 20.0, 30.0]), [1.41421, 1.38493, 1.24035, 1.11035], atol=1e-5)
    o = dict(f1=np.array([-3.0, 0.0, 2.5]), K=np.array([K, 40.0, 250.0]))
    o0 = dict(f1=o['f1'], K=np.zeros(3))
    for w_ in (4.0, 13.1, 37.0):
        r = np.abs(NR.transfer(o, np.array([w_]), ke)[0][:, 0] / NR.transfer(o0, np.array([w_]), ke)[0][:, 0])
        assert np.allclose(r, NR.sens(o['K'], ke, w_), rtol=1e-12)
    for K_, ke_ in ((248.3, 10.0), (100.0, 5.0), (10.0, 40.0), (1.0, 1.0)):
        assert abs(NR.bode_integral(K_, ke_)) < 1e-8
    rows = [dict(m=m, band=[10.0, 30.0], sqrtK=[sk, 1.0, 1.0], peak=[pk, 4.0, 5.0], jit=[1.0, 1.0, 1.0], jit0=[1.0, 1.0, 1.0])
            for m, sk, pk in ((4.83, 10.0, 8.35), (12.0, 15.76, 14.15))]
    md = '\n'.join(NR.a2_lines(rows))
    assert '|S(√K)| = (1 + K/k_e²)^{1/2} = 1.414, 1.866 at m = 4.83, 12' in md
    assert 'at m = 4.83 (√K₁ = 10.00) |S| = 1.41, 1.38, 1.24, 1.11 at ω = 10, 15, 20, 30' in md
    assert '4.00* (4.000)' in md and '8.35 (0.835)' in md and '5.00 (5.000)' in md and '5.00*' not in md


def test_c_noise_shifted_is_the_draw_from_a_later_start():
    """Referee round 69 (C7): opts shift starts the noise shift periods of the reference into the draw — each phase moved by ωt₀ —,
    so the shifted noise at t is the draw's at t + t₀; with neg and freqs as before."""
    T = 2 * np.pi / CT.OM
    a, b = NR.c_noise(0.02, dict(seed=7)), NR.c_noise(0.02, dict(seed=7, shift=25))
    n = lambda nz, t: float(np.sum(nz['amp'] * np.cos(nz['w'] * t + nz['ph'])))
    for t in (0.0, 1.3, 17.0):
        assert np.isclose(n(b, t), n(a, t + 25 * T), atol=1e-12)
    assert np.allclose(b['w'], a['w']) and np.allclose(b['amp'], a['amp'])
    t1 = NR.c_noise(0.02, dict(freqs=(13.1,), seed=7, shift=3))
    t0 = NR.c_noise(0.02, dict(freqs=(13.1,), seed=7))
    assert np.isclose(t1['ph'][0] - t0['ph'][0], 13.1 * 3 * T)
    assert NR.c_noise(0.0, dict(seed=7, shift=50)) is None and NR.c_noise(0.02, dict(seed=7, shift=0))['ph'] is not None


def test_bracket_edge_on_a_step_and_its_rechecks():
    """Referee round 69: `bracket_edge` goes up the grid to the first ejecting raise, rechecks the next `above` raises, bisects nb times,
    and keeps every run; a step at 0.0731 is bracketed to the grid step over 2^nb; nothing ejecting gives hi None."""
    calls = []

    def run(d):
        calls.append(d)
        return (2 if d > 0.0731 else -1), (0.9 if d > 0.0731 else 0.001)
    lo, hi, runs = NR.bracket_edge(run, grid=(0.02, 0.04, 0.06, 0.08, 0.1, 0.13), nb=6, above=2)
    assert lo < 0.0731 < hi and np.isclose(hi - lo, 0.02 / 64)
    assert calls[:6] == [0.02, 0.04, 0.06, 0.08, 0.1, 0.13] and len(calls) == 12
    assert [x[0] for x in runs] == sorted(calls)
    lo, hi, runs = NR.bracket_edge(lambda d: (-1, 0.0), grid=(0.1, 0.2), nb=3)
    assert hi is None and lo == 0.2 and len(runs) == 2
    lo, hi, runs = NR.bracket_edge(lambda d: (1, 1.0), grid=(0.1, 0.2), nb=2, above=1)
    assert lo == 0.0 and np.isclose(hi, 0.025) and [x[0] for x in runs] == [0.025, 0.05, 0.1, 0.2]


def test_the_edges_scaling_with_the_centre_mobility():
    """Referee rounds 69, 76 (EDGE): `v_of` is a formed kernel's centre mobility at the fit — Gmw under the paradigm's metric (∝ m), G·a
    under WFR's —; on edges following Δ_e = 0.08(k_e/v₁)^{1/2} exactly the window's power law has α = β = 1/2 (kernel 3, WFR and ω_ref
    rows kept out of the fit; the bracket's midpoint), and each row outside the window prints against the law 0.08(k_e/v)^{1/2} and,
    for the paradigm's kernel 1, the power law; the ω_ref slopes; `star_verdict` sorts a first jump against s*'s bracket."""
    v12, v3 = NR.v1_of(12.0, 'paradigm'), NR.v1_of(3.0, 'paradigm')
    assert np.isclose(v12, 12.0 * CT.ASTAR[0], rtol=1e-6) and np.isclose(v12 / v3, 4.0, rtol=1e-9)
    assert 0.9 < NR.v1_of(12.0, 'wfr') < 1.1 and np.isclose(NR.v_of(12.0, 'paradigm', 2), 12.0 * CT.ASTAR[2], rtol=1e-6)
    row = lambda m, ke, d, w=1e-7, **kw: dict(dict(m=m, ke=ke, metric='paradigm', k=0, lo=d - w, hi=d + w, runs=[]), **kw)
    E_ = [row(m, ke, 0.08 * np.sqrt(ke / NR.v1_of(m, 'paradigm')), w=0.004) for m in (3.0, 12.0) for ke in (5.0, 10.0, 20.0)]
    d3 = 0.5 * np.sqrt(10.0 / NR.v_of(12.0, 'paradigm', 2))
    E_ += [row(12.0, 10.0, d3, k=2), row(12.0, 80.0, 1.5 * 0.08 * np.sqrt(80.0 / v12), w=0.01),
           dict(m=12.0, ke=40.0, metric='wfr', k=0, lo=1.0, hi=None, runs=[]),
           row(12.0, 10.0, 0.08 * np.sqrt(10.0 / v12) * 0.5 ** 0.75, om=0.1), row(12.0, 10.0, 0.08 * np.sqrt(10.0 / v12) * 2 ** 0.75, om=0.4),
           row(24.0, 10.0, 0.9 * 0.08 * np.sqrt(10.0 / NR.v1_of(24.0, 'paradigm')), w=0.002)]
    md = '\n'.join(NR.edge_lines(E_, None, None, None))
    assert f"m = 24, k_e = 10: {0.9 * 0.08 * np.sqrt(10.0 / NR.v1_of(24.0, 'paradigm')):.4f}, 0.90 times the law, the power law's 0.90" in md
    assert 'over the paradigm\'s 6 kernel-1 rows in m = 3–12, k_e = 5–40' in md and 'has α = 0.500 ±' in md and 'β = 0.500 ±' in md
    assert 'Δ_e·(v₁/k_e)^{1/2} = 0.0800 … 0.0800 over the window (median 0.0800)' in md
    assert f"m = 12, k_e = 80: {1.5 * 0.08 * np.sqrt(80.0 / v12):.4f}, 1.50 times the law, the power law's 1.50" in md
    assert f"m = 12, k_e = 10, kernel 3: {d3:.4f}, 6.25 times the law" in md
    assert f"m = 12, k_e = 40, wfr: > 1 (the law {0.08 * np.sqrt(40.0 / NR.v1_of(12.0, 'wfr')):.3f})" in md
    assert "the log-slope by least squares (pairwise): m = 12, k_e = 10: 0.75 (0.75, 0.75)" in md
    fit = NR.edge_fit(E_)
    assert len(fit['rows']) == 6 and np.isclose(fit['alpha'], 0.5) and np.isclose(fit['beta'], 0.5)
    assert NR.edge_fit(E_[6:]) is None
    assert NR.star_verdict(0.018, 0.019, 0.02, [0.015, 0.02]) == 'consistent'
    assert NR.star_verdict(0.035, 0.036, 0.08, [0.04, 0.06, 0.08]) == 'later'
    assert NR.star_verdict(0.008, 0.009, 0.0075, [0.005, 0.0075]) == 'earlier'
    assert NR.star_verdict(0.07, 0.072, None, [0.04, 0.06, 0.08]) == 'none, predicted'
    assert NR.star_verdict(0.09, 0.1, None, [0.04, 0.06, 0.08]) == 'untested'


def test_c2_single_tones_rise_on_the_raw():
    """Referee round 69: C2 prints, per level, the first coefficient's rise over the fit by the last period in the single-tone runs (the
    paradigm's metric) whose first kernel did not jump — the drift that never reaches the tones' moving edge."""
    raw = NR.E.load(os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'noise_resonance_raw.npz'))
    md = '\n'.join(NR.c2_lines(raw['C2']))
    rise = {}
    for r in raw['C2']:
        o = r['opts']
        if 'freqs' in o and len(o['freqs']) == 1 and o.get('metric', 'paradigm') == 'paradigm' and r['m'] == 12.0 and r['ke'] == 10.0:
            if not (r['left'][0] >= 0 and r['maxdev'][0] > NR.JUMP_C):
                rise.setdefault(r['s'], []).append(float(r['w_end'][0] - CT.ASTAR[0]))
    assert rise and all(len(v) == 24 for v in rise.values())
    for s_, v in rise.items():
        assert f"s = {s_:g}: {min(v):+.3f} … {max(v):+.3f} ({len(v)} runs)" in md


def test_star_cases_count_and_null():
    """Referee round 76 (s*, ED01–ED05): kernel 3's s* from its own ⟨w₂⟩ and its own jumps; `star_count`'s verdicts with replaced ⟨w₂⟩;
    `star_null` permutes the draws' ⟨w₂⟩ within each configuration only."""
    n = 300
    E_ = [dict(m=12.0, ke=10.0, metric='paradigm', k=0, lo=0.07, hi=0.075, runs=[]),
          dict(m=12.0, ke=10.0, metric='paradigm', k=2, lo=0.074, hi=0.075, runs=[])]
    W2 = []
    for i, sd in enumerate(NR.W2_SEEDS):
        ew = np.zeros((n, 3)); ew[:, 0] = 200.0; ew[:, 2] = 100.0 * (i + 1)
        W2.append(dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=ew.tolist()))
    C = [dict(m=12.0, ke=10.0, s=s_, opts=dict(seed=sd), left=[5 if s_ >= 0.02 else -1, -1, 5 if s_ >= 0.03 else -1],
              maxdev=[0.9 if s_ >= 0.02 else 0, 0, 0.9 if s_ >= 0.03 else 0]) for sd in NR.W2_SEEDS for s_ in (0.015, 0.02, 0.03)]
    cases = NR.star_cases(E_, W2, C, None)
    k3 = [q for q in cases if q['key'][3] == 2]
    assert [q['w2'] for q in k3] == [100.0, 200.0, 300.0, 400.0] and all(q['obs'] == 0.03 for q in k3)
    assert all(q['obs'] == 0.02 for q in cases if q['key'][3] == 0 and q['key'][4] is None)
    assert all(q['obs'] is None and not q['run'] for q in cases if q['key'][4] is not None)   # W2's band rows: no runs here
    rec = NR.star_count(cases)
    v3 = [v for q, v, _ in rec if q['key'][3] == 2]
    sstar = [np.sqrt(0.0745 / x) for x in (100.0, 200.0, 300.0, 400.0)]       # 0.0273, 0.0193, 0.0158, 0.0136
    assert v3 == ['consistent', 'later', 'later', 'later'] and sstar[0] > 0.02 > sstar[1]
    alt = NR.star_count(cases, [q['w2'] if q['key'][3] == 0 else 100.0 for q in cases])
    assert [v for q, v, _ in alt if q['key'][3] == 2] == ['consistent'] * 4
    null = NR.star_null(cases, n=50, seed=1)
    assert len(null) == 50 and all(x == 5 for x in null)
    md = '\n'.join(NR.edge_lines(E_, W2, C, None))
    w1, w3 = 0.005 / (0.5 * (np.sqrt(0.07 / 200) + np.sqrt(0.075 / 200))), 0.01 / (0.5 * (np.sqrt(0.074 / 100) + np.sqrt(0.075 / 100)))
    assert f"are {w1:.2f}–{w3:.2f} of s* wide" in md and (f"{w1:.2f}", f"{w3:.2f}") == ('0.26', '0.37')
    assert "the consistent count is 5.00 on average (5–5; at least 5 in 200)" in md
    assert ("Within a configuration the draws' s* (its bracket's midpoint) differ by at most 100% (median 50%, over the 2 configurations "
            "with two or more draws), against the level grid's steps of 33–50% — by less than the smallest step in 1 of the 2 —: the null "
            "has little power.") in md


def test_edge_fit_standard_errors():
    """Referee round 76: `edge_fit`'s coefficients and standard errors are ordinary least squares' with n − 3 degrees of freedom, over
    the window's rows only."""
    rng = np.random.default_rng(3)
    rows = []
    for m in (3.0, 6.0, 12.0):
        for ke in (5.0, 10.0, 20.0, 40.0):
            d = 0.05 * ke ** 0.5 * m ** -0.5 * np.exp(0.05 * rng.standard_normal())
            rows.append(dict(m=m, ke=ke, metric='paradigm', k=0, lo=d - 1e-9, hi=d + 1e-9, runs=[]))
    fit = NR.edge_fit(rows + [dict(rows[0], m=24.0), dict(rows[1], ke=80.0), dict(rows[2], metric='wfr'), dict(rows[3], k=2)])
    X = np.array([[1.0, np.log(r['ke']), np.log(r['m'])] for r in rows])
    y = np.log([0.5 * (r['lo'] + r['hi']) for r in rows])
    cf = np.linalg.solve(X.T @ X, X.T @ y)
    res = y - X @ cf
    se = np.sqrt(np.diag(res @ res / (len(y) - 3) * np.linalg.inv(X.T @ X)))
    assert len(fit['rows']) == 12 and np.allclose([fit['c'], fit['alpha'], -fit['beta']], cf, rtol=1e-12)
    assert np.allclose(fit['se'], se, rtol=1e-9) and np.allclose(fit['res'], res, atol=1e-12)


def test_bands_of_and_scan_sets():
    """Referee round 76: `bands_of` gives the maximal runs of flagged raises; `scan_sets` the leaving and the jumping bands of a scan."""
    ds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert NR.bands_of(ds, [False, True, False, True, True, True]) == [(0.1, 0.1), (0.3, 0.5)]
    assert NR.bands_of(ds, [False] * 6) == [] and NR.bands_of(ds, [True] * 6) == [(0.0, 0.5)]
    r = dict(runs=[dict(d=d, left=(3 if d in (0.1, 0.3, 0.4, 0.5) else -1), maxdev=(0.9 if d >= 0.4 else 0.4)) for d in reversed(ds)])
    assert NR.scan_sets(r) == ([(0.1, 0.1), (0.3, 0.5)], [(0.4, 0.5)])
    assert NR.fmt_bands([(0.1, 0.1), (0.3, 0.5)]) == '0.1, 0.3–0.5' and NR.fmt_bands([]) == 'none'


def test_c7_job_scans_and_lines_with_a_fake_run(monkeypatch):
    """Referee rounds 69, 76 (C7): `c7_job` runs `edge_run` under the case's noise (the window's shift, the tone's sign) over NP_C7 periods
    at every raise of SCAN_C7 (kernel 1; m = 12, k_e = 10, the paradigm's metric unless given); `c7_lines` prints the noise-free sets, the
    tones' least raises under ±n with the even part, E and the odd part ('>1' kept), and the windows' least raises with their band counts,
    their least and median ('>1' above every raise) and the windows that are not an interval."""
    seen = []
    monkeypatch.setattr(NR, 'SCAN_C7', (0.0, 0.065, 0.0675, 0.07, 0.0725, 0.1, 0.2, 0.5, 1.0))
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: None if s == 0 else dict(
        tag=(tuple(o.get('freqs', ())), bool(o.get('neg')), o.get('shift', 0), s)))

    def fake(m, ke, metric, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None):
        seen.append((m, ke, metric, k, NP, None if noise is None else noise['tag']))
        if noise is None:
            ej, md = d >= 0.0725, 0.9
        else:
            fr, neg, sh, s = noise['tag']
            if fr:
                th = {(False, 0.01): 0.1, (True, 0.01): 0.2, (False, 0.02): 2.0, (True, 0.02): 0.5}[(neg, s)]
                ej, md = d >= th, 0.9
            else:
                ej = sh == 50 and (d == 0.0675 or d >= 0.1)
                md = 0.4 if d == 0.0675 else 0.9
        return (3 if ej else -1), (md if ej else 0.01)
    monkeypatch.setattr(NR, 'edge_run', fake)
    cases = [(dict(seed=7), 0.0), (dict(freqs=(13.1,), seed=7), 0.01), (dict(freqs=(13.1,), seed=7, neg=True), 0.01),
             (dict(freqs=(13.1,), seed=7), 0.02), (dict(freqs=(13.1,), seed=7, neg=True), 0.02),
             (dict(seed=7, shift=0), 0.02), (dict(seed=7, shift=50), 0.02)]
    rs = [NR.c7_job(c) for c in cases]
    assert all(x[:5] == (12.0, 10.0, 'paradigm', 0, NR.NP_C7) for x in seen) and len(seen) == 9 * len(cases)
    assert {x[5] for x in seen if x[5] is not None and x[5][0] == ()} == {((), False, 0, 0.02), ((), False, 50, 0.02)}
    assert rs[6]['opts'] == dict(seed=7, shift=50) and [x['d'] for x in rs[6]['runs']] == list(NR.SCAN_C7)
    md = '\n'.join(NR.c7_lines(rs))
    assert "| 12 | 10 | paradigm | 0.0725–1 | 0.0725–1 |" in md
    assert "| one tone 13.1 | 0.01 | 0.1 (0.1–1) | 0.2 (0.2–1) | 0.15 | +775 | -5.000 |" in md
    assert "| one tone 13.1 | 0.02 | >1 (none) | 0.5 (0.5–1) | >1 | — | — |" in md
    assert "E per tone and level (the even part's rise per s²; '—' where a least raise is beyond the grid): one tone 13.1: +775 at 0.01." in md
    assert "| 12 | 10 | paradigm | 0.02 | >1, 0.0675 [2] | 0.0675, >1 | 1 of 2 | 0.1, >1 |" in md
    assert ("Without noise the leaving set is an interval from the edge to the grid's end at every configuration; under the band's noise it is not an interval in 1 of the 2 windows: m = 12, k_e = 10, s = 0.02, from period 50: 0.0675, 0.1–1.") in md
    rr = NR.c7_job((dict(seed=7, shift=25), 0.04, (12.0, 20.0, 'paradigm')))
    assert seen[-1][:4] == (12.0, 20.0, 'paradigm', 0) and rr['ke'] == 20.0 and rr['s'] == 0.04


def test_c8_verdict_confusion_and_flags():
    """Referee round 76 (C8): a run's first flag against its first leave — hit in the leave's window (p < P ≤ p + STEP_C8), early, miss;
    for a run that never left correct or a false alarm —; the windows' confusion counts up to the leave's window; the rules' flags."""
    st = list(range(0, 300, 10))
    assert NR.c8_verdict(133, [130]) == ('hit', 130) and NR.c8_verdict(130, [120]) == ('hit', 120) and NR.c8_verdict(131, [130]) == ('hit', 130)
    assert NR.c8_verdict(133, [120, 130]) == ('early', 120) and NR.c8_verdict(133, [140, 150]) == ('miss', None)
    assert NR.c8_verdict(-1, []) == ('correct', None) and NR.c8_verdict(-1, [50, 60]) == ('false alarm', 50)
    assert list(NR.c8_confusion(133, [120, 130, 200], st)) == [1, 0, 1, 12]
    assert list(NR.c8_confusion(133, [], st)) == [0, 1, 0, 13] and list(NR.c8_confusion(-1, [50], st)) == [0, 0, 1, 29]
    r = dict(wins=[dict(p=0, d=0.0, left=-1, maxdev=0.01), dict(p=10, d=0.08, left=4, maxdev=0.4), dict(p=20, d=0.06, left=2, maxdev=0.9)])
    assert NR.c8_flags(r, 'leave') == [10, 20] and NR.c8_flags(r, 'jump') == [20]
    assert NR.c8_flags(r, 'the noise-free edge', 0.07) == [10] and NR.c8_flags(r, 'the noise-free edge', None) == []


def test_c8_job_and_lines_with_a_fake_run(monkeypatch):
    """Referee round 76 (C8): `c8_job` runs C's loop (`c_run`, the metric passed) and every STEP_C8 periods p `edge_run` with kernel 1 raised
    by the run's own δ(p) = |W₁| over period p minus its fitted value (0 at p = 0) under the run's noise shifted by p, over NP_C7 periods;
    `c8_lines` prints each run's row and the rules' tallies in and out of sample."""
    L, al0, c0, idx = _fit()
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    Wp = np.tile(CT.ASTAR, (NR.NP_C, 1)) + np.linspace(0.0, 0.15, NR.NP_C)[:, None] * np.array([1.0, 0, 0])[None, :]
    Cp = np.tile(CT.CSTAR, (NR.NP_C, 1)); Cp[132:, 0] -= 0.9
    got = {}

    def fake_run(m, ke, s, opts):
        got.update(m=m, ke=ke, s=s, opts=dict(opts))
        return None, L, al0, c0, idx, Wp, Cp
    calls = []
    monkeypatch.setattr(NR, 'c_run', fake_run)
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: dict(shift=o.get('shift', 0), s=s))
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None: (
        calls.append((m, ke, mt, k, NP, noise['shift'], d)), ((3, 0.9) if d >= 0.055 else (-1, 0.01)))[1])
    r = NR.c8_job((12.0, 10.0, 'wfr', 0.02, dict(seed=3)))
    assert got == dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=3, metric='wfr')) and r['left'] == [133, -1, -1] and np.isclose(r['w0'], w0)
    assert [c[5] for c in calls] == list(range(0, NR.NP_C, NR.STEP_C8)) and all(c[:5] == (12.0, 10.0, 'wfr', 0, NR.NP_C7) for c in calls)
    assert calls[0][6] == 0.0 and np.isclose(calls[13][6], abs(Wp[129, 0]) - w0) and len(r['wins']) == len(calls)
    first = next(x['p'] for x in r['wins'] if x['left'] >= 0)
    rr = [r, dict(r, opts=dict(seed=4)), dict(r, opts=dict(seed=7), left=[-1, -1, -1], maxdev=[0.1, 0, 0])]
    E_ = [dict(m=12.0, ke=10.0, metric='wfr', k=0, lo=0.09, hi=0.11, runs=[])]
    md = '\n'.join(NR.c8_lines(rr, E_))
    assert first == 120 and "| 12 | 10 | wfr | 0.02 | 3 | 133 j |" in md
    assert "120, 130, 140, 150 … (18): early" in md and "| leave | out | 0 | 2 | 0 | 0 | 0 |" in md and "| 12 | 10 | wfr | 0.02 | 4 | 133 j |" in md
    assert "| leave | in | 0 | 0 | 0 | 0 | 1 |" in md and "| the noise-free edge | out |" in md
    assert ("The totals: in sample 1 runs, 0 of them leaving, and 30 earlier windows; out of sample 2 runs, 2 of them leaving, and 26 "
            "earlier windows.") in md


def test_b5_lines_summary_over_every_mark():
    """Referee round 76 (B5L1): B5's summary at the small level spans every mark and kernel, not the last mark's."""
    nP = len(NR.MARKS_B)
    a = np.ones((nP, 3)); es = np.ones((nP, 3)); es[0, 0] = 0.9; es[-1, 2] = 1.1
    r = dict(metric='paradigm', om=13.1, w2=a.tolist(), w2_shift=a.tolist(), w_env=a.tolist(), w_env_small=es.tolist())
    md = '\n'.join(NR.b5_lines([r]))
    assert "13.1: 0.9000–1.1000 (over the exact second order 0.9000–1.1000)" in md


def test_a2_rows_reference_includes_the_plant_slope(monkeypatch):
    """Referee round 76 (AR1): A2's K = 0 reference is the roll-off (iω − f̂′)/(−ω² + ik_eω), the plant's slope f̂′ included."""
    from src.core import closed_loop as CL
    monkeypatch.setattr(NR, 'M_A2', (12.0,))
    monkeypatch.setattr(NR, 'OM_A2', np.linspace(4.0, 40.0, 13))
    rows = NR.a2_rows()
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    xm = np.sin(np.linspace(0, 2 * np.pi, 2001)[:-1])
    o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', xm)
    nz = CL.tones(1.0, *NR.BANDS_A2[0], seed=7)
    om, amp = np.asarray(nz['w'], float), np.asarray(nz['amp'], float)
    Y0 = (1j * om[None, :] - o['f1'][:, None]) / (-om[None, :] ** 2 + 1j * CT.KE * om[None, :])
    g = o['v'][None, :] * o['d1']
    j0 = np.sqrt(np.mean((np.abs(g[:, :, None] * Y0[:, None, :] * amp[None, None, :]) ** 2 / 2).sum(2), axis=0))
    assert np.allclose(rows[0]['jit0'], j0, rtol=1e-12) and np.abs(o['f1']).max() > 0.1


def test_w2_cases_cover_the_band_draws_and_every_configuration():
    """Referee round 76 (W2C): W2's draws are C3's eleven band draws at m = 12, k_e = 10 and W2_SEEDS for every other configuration."""
    cs = NR.w2_cases()
    band = [c for c in cs if c[:3] == (12.0, 10.0, 'paradigm') and 'band' not in c[3]]
    assert sorted(c[3]['seed'] for c in band) == list(range(1, 12))
    others = {(m, ke, mt, b) for m, ke, mt, b in NR.W2_CONFIGS} - {(12.0, 10.0, 'paradigm', None)}
    assert len(cs) == 11 + len(NR.W2_SEEDS) * len(others)


def _c9_runs():
    """Synthetic C3/C6 runs of draw 3 for C9: kernel 1's |w| = w₀ + 0.001·period (+ add), leaving at P; and runs C9 must leave out."""
    L = NR.law(NR.C7_CFG[0])
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    w0 = float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]])

    def run(opts, P, s=0.02, add=0.0):
        W = np.tile(CT.ASTAR, (NR.NP_C, 1)).astype(float)
        W[:, 0] = np.sign(w0) * (abs(w0) + 0.001 * np.arange(NR.NP_C) + add)
        return dict(opts=opts, s=s, left=[P, -1, -1], W=W.tolist())
    C3 = [run(dict(seed=3), 157), run(dict(seed=3), 155, s=0.03), run(dict(freqs=(13.1,), seed=3), 155),
          run(dict(seed=3, neg=True), 155), run(dict(seed=7), 155)]
    C6 = [run(dict(seed=3, rise=[0, 0.05]), 166, add=0.02), run(dict(seed=3, rise=[0, 0.03]), 120, add=0.01),
          run(dict(seed=3, rise=[0, 0.02]), 150, add=0.005), run(dict(seed=3, rise=[0, 0.04]), 160, add=0.015)]
    return C3, C6, abs(w0)


def test_c9_pairs_own_job_and_part(monkeypatch):
    """Referee round 76 (C9D): `c9_pairs` takes, per window (seed, p), the draw's s = 0.02 runs without tones or −n that have not left by
    p (a leave at p itself excluded), with δ(p) = |W₁| over the period before p minus the fitted value; `c9_own_job` runs `edge_run` at
    C7_CFG with kernel 1 raised by δ under the draw's noise from p for NP_C7 periods; `part_C9D` maps it over every pair of the default
    windows."""
    C3, C6, w0 = _c9_runs()
    pr = NR.c9_pairs(C3, C6, [(3, 150), (3, 160), (3, 165), (3, 100)])
    assert [(x[0], x[1], x[2], x[4]) for x in pr] == [(3, 150, 0.0, 157), (3, 150, 0.04, 160), (3, 150, 0.05, 166), (3, 160, 0.05, 166),
                                                      (3, 165, 0.05, 166), (3, 100, 0.0, 157), (3, 100, 0.02, 150), (3, 100, 0.03, 120),
                                                      (3, 100, 0.04, 160), (3, 100, 0.05, 166)]
    assert np.allclose([x[3] for x in pr], [0.149, 0.164, 0.169, 0.179, 0.184, 0.099, 0.104, 0.109, 0.114, 0.119], atol=1e-12)
    calls = []
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: dict(s=s, **o))
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None: (
        calls.append((m, ke, mt, k, d, NP, noise)), (4, 0.7))[1])
    r = NR.c9_own_job((3, 150, 0.05, 0.169))
    assert calls == [(12.0, 10.0, 'paradigm', 0, 0.169, NR.NP_C7, dict(s=0.02, seed=3, shift=150))]
    assert r == dict(seed=3, p=150, rise=0.05, d=0.169, left=4, maxdev=0.7)
    got = []
    monkeypatch.setattr(NR, '_pool_map', lambda fn, jobs: got.append((fn, list(jobs))) or [])
    NR.part_C9D(C3, C6)
    assert got[0][0] is NR.c9_own_job and got[0][1] == [x[:4] for x in NR.c9_pairs(C3, C6)]
    assert {(x[0], x[1]) for x in got[0][1]} <= {(o['seed'], int(o['shift'])) for o, _ in NR.C9_CASES} and len(got[0][1]) > 0


def test_c9_lines_counts_and_the_runs_flagged():
    """Referee round 76 (C9, C9D): the table's cells (δ, whether the nearest raise leaves — a leave without a jump counts —, the loop at
    δ, the run — its leave at the window's end counts —), the 2×2 counts against the run's leave by the nearest raise and by the loop at
    δ, the misses, and the runs flagged in some window containing their leave (an earlier window's flag counts)."""
    C3, C6, _ = _c9_runs()
    sc = lambda sh, lv, md=None: dict(opts=dict(seed=3, shift=sh), s=0.02, runs=[
        dict(d=d, left=l_, maxdev=(md or {}).get(d, 0.9 if l_ >= 0 else 0.01)) for d, l_ in zip((0.14, 0.15, 0.16, 0.17, 0.18, 0.19), lv)])
    C9 = [sc(165, (-1, -1, -1, -1, -1, 1)), sc(160, (-1, -1, -1, -1, 2, 2)), sc(150, (3, 3, -1, 2, -1, -1), {0.17: 0.4})]
    md = '\n'.join(NR.c9_lines(C9, C3, C6))
    assert "| 3 | 150 | 0.14–0.15, 0.17 | 0: +0.1490, yes, 157; 0.04: +0.1640, no, 160; 0.05: +0.1690, yes, — |" in md
    assert "| 3 | 160 | 0.18–0.19 | 0.05: +0.1790, yes, 166 |" in md and "| 3 | 165 | 0.19 | 0.05: +0.1840, no, 166 |" in md
    assert md.index("| 3 | 150 |") < md.index("| 3 | 160 |") < md.index("| 3 | 165 |")
    assert ("Of the 5 (run, window) pairs, by the grid's nearest raise: 2 in the leaving set with the run leaving there, 1 in the set but "
            "staying, 2 leaving though not in the set, 0 neither.") in md and "By the window's loop" not in md
    D = [dict(seed=3, p=150, rise=0.0, d=0.149, left=6, maxdev=0.9), dict(seed=3, p=150, rise=0.04, d=0.164, left=9, maxdev=0.9),
         dict(seed=3, p=150, rise=0.05, d=0.169, left=-1, maxdev=0.01), dict(seed=3, p=160, rise=0.05, d=0.179, left=6, maxdev=0.9),
         dict(seed=3, p=165, rise=0.05, d=0.184, left=-1, maxdev=0.01)]
    md = '\n'.join(NR.c9_lines(C9, C3, C6, D))
    assert "| 3 | 150 | 0.14–0.15, 0.17 | 0: +0.1490, yes, 156, 157; 0.04: +0.1640, no, 159, 160; 0.05: +0.1690, yes, —, — |" in md
    assert "| 3 | 160 | 0.18–0.19 | 0.05: +0.1790, yes, 166, 166 |" in md and "| 3 | 165 | 0.19 | 0.05: +0.1840, no, —, 166 |" in md
    assert ("By the window's loop at δ itself: 3 leaving with the run, 0 leaving while the run stays, 1 staying while the run leaves, "
            "1 neither; the run leaving though the loop at δ stays: draw 3 raised by 0.05 in the window from period 165 (it leaves at "
            "166). Of the 3 runs whose leave a window contains, 3 are flagged — the loop at δ leaving — in at least one of those windows.") in md
    D[3] = dict(D[3], left=-1)
    md = '\n'.join(NR.c9_lines(C9, C3, C6, D))
    assert ("2 leaving with the run, 0 leaving while the run stays, 2 staying while the run leaves, 1 neither; the run leaving though the "
            "loop at δ stays: draw 3 raised by 0.05 in the window from period 160 (it leaves at 166); draw 3 raised by 0.05 in the window "
            "from period 165 (it leaves at 166). Of the 3 runs whose leave a window contains, 2 are flagged — the loop at δ leaving — in at "
            "least one of those windows: not draw 3 raised by 0.05.") in md
    D[2] = dict(D[2], left=2)
    assert "2 leaving with the run, 1 leaving while the run stays, 2 staying while the run leaves, 0 neither" in '\n'.join(
        NR.c9_lines(C9, C3, C6, D))
    assert NR.c9_lines([], C3, C6, D) == [] and NR.c9_lines(None) == []

def test_main_passes_c3_and_c6_to_c9d(monkeypatch, tmp_path):
    """Referee round 76: `main --only C9D` computes C9D from the raw's C3 and C6 and keeps every other part."""
    old = dict(C3=['c3'], C6=['c6'], C9=['c9'])
    saved, got = [], []
    monkeypatch.setattr(NR.E, 'load', lambda p: dict(old))
    monkeypatch.setattr(NR.E, 'save', lambda p, **kw: saved.append(kw))
    monkeypatch.setattr(NR, 'OUT', str(tmp_path))
    monkeypatch.setattr(NR, 'report', lambda **kw: '')
    monkeypatch.setattr(NR, 'part_C9D', lambda C3, C6: got.append((C3, C6)) or ['d'])
    monkeypatch.setattr(sys, 'argv', ['exp_noise_resonance.py', '--only', 'C9D'])
    NR.main()
    assert got == [(['c3'], ['c6'])] and saved[-1] == dict(C3=['c3'], C6=['c6'], C9=['c9'], C9D=['d'])


def test_edge_run_from_an_explicit_start():
    """Referee round 85: `edge_run`'s start — a callable giving the state the loop starts from — replaces the raised fitted state; a start
    returning the raised state reproduces the default run bit for bit."""
    a = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05, NP=2)
    b = NR.edge_run(12.0, 10.0, 'paradigm', 0, 999.0, NP=2, start=lambda L, al0, c0, idx: (NR.raised(L, al0, idx, 0, 0.05), c0))
    assert a == b and a != NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.2, NP=2)
    c = NR.edge_run(12.0, 10.0, 'paradigm', 0, 999.0, NP=2, start=lambda L, al0, c0, idx: (NR.raised(L, al0, idx, 0, 0.05), c0 + 0.03))
    assert c != b                                                       # the start's centres are the loop's (ES2)


def test_c7_job_window_length(monkeypatch):
    """Referee round 85 (D3cf-1): `c7_job`'s fourth argument is the window's length — passed to `edge_run` as NP and recorded as 'np' only
    where it is not NP_C7 —, a None configuration C7_CFG; C7L_CASES are C7's 18 gap windows, four one-band controls and C9's draw 3
    from period 150, over twice the window."""
    calls = []
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, noise=None, **kw: calls.append((m, ke, mt, NP)) or (-1, 0.01))
    monkeypatch.setattr(NR, 'SCAN_C7', (0.0, 0.1))
    r = NR.c7_job((dict(seed=7, shift=25), 0.01, None, 20))
    assert calls == [(12.0, 10.0, 'paradigm', 20)] * 2 and r['np'] == 20 and (r['m'], r['ke'], r['metric']) == NR.C7_CFG
    calls.clear()
    r = NR.c7_job((dict(seed=7), 0.01, (4.0, 10.0, 'paradigm')))
    assert calls == [(4.0, 10.0, 'paradigm', NR.NP_C7)] * 2 and 'np' not in r
    assert len(NR.C7L_CASES) == 23 and all(c[3] == 2 * NR.NP_C7 for c in NR.C7L_CASES)
    assert sum(1 for c in NR.C7L_CASES if c[2] is None and c[1] == 0.02 and c[0].get('shift') in (0, 50, 75, 100)) == 4
    assert (dict(seed=3, shift=150), 0.02, None, 2 * NR.NP_C7) in NR.C7L_CASES


def test_c8_windows_nulls_dense_reach_offsets_and_c9a(monkeypatch):
    """Referee round 85 (D3cf-2 … D3cf-4): `c8_cases_of`, `c8_delta` and `win_run` reproduce C8's stored windows bit for bit; with a fake
    window, C8N's windows (δ = 0 under the coming noise, the run's δ under the draw's noise C8N_OFF periods later), C8D's at every period
    up to the leave's (every window for a run that never left), C8R's starts h periods before the leave (its leave counted from the
    run's start), C8K's start from every formed kernel's offsets (and centres), and C9A's starts; `c8k_cases` pairs C8's runs with C3's."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    cases = NR.c8_cases_of(raw['C8'])
    assert len(cases) == len(raw['C8']) and cases[0][5] == list(raw['C8'][0]['w']) and cases[0][7] == raw['C8'][0]['left'][0]
    for x in raw['C8'][0]['wins'][:2]:
        d = NR.c8_delta(cases[0][5], cases[0][6], x['p'])
        assert d == x['d'] and NR.win_run(*cases[0][:5], d, x['p']) == (x['left'], x['maxdev'])
    assert NR.c8_delta([0.3, -0.4], 0.25, 0) == 0.0 and NR.c8_delta([0.3, -0.4], 0.25, 2) == 0.4 - 0.25
    calls = []

    def fake(m, ke, mt, s, opts, d, shift, NP=NR.NP_C7, start=None):
        calls.append((d, shift, NP, start is not None))
        return ((5 if shift == 40 else -1), 0.1 + d)
    monkeypatch.setattr(NR, 'win_run', fake)
    w, w0 = [0.3 + 0.01 * i for i in range(NR.NP_C)], 0.28
    a = (12.0, 10.0, 'paradigm', 0.02, dict(seed=1), w, w0, 57)
    out = NR.c8n_job(a)
    ps = list(range(0, NR.NP_C, NR.STEP_C8))
    assert calls[:len(ps)] == [(0.0, p, NR.NP_C7, False) for p in ps]
    assert calls[len(ps):] == [(NR.c8_delta(w, w0, p), p + o, NR.NP_C7, False) for o in NR.C8N_OFF for p in ps]
    assert [x['p'] for x in out['zero']] == ps and out['zero'][4]['left'] == 5 and list(out['off']) == [str(o) for o in NR.C8N_OFF]
    calls.clear()
    out = NR.c8d_job(a)
    assert [c[1] for c in calls] == list(range(0, 57)) and [x['d'] for x in out['wins']] == [NR.c8_delta(w, w0, p) for p in range(57)]
    calls.clear()
    out = NR.c8d_job(a[:7] + (-1,))
    assert [c[1] for c in calls] == list(range(0, NR.NP_C - NR.NP_C7 + 1))
    calls.clear()
    out = NR.c8r_job(a[:7] + (60,))
    assert [(c[1], c[2]) for c in calls] == [(60 - h, h + 2) for h in NR.C8R_BACK if h <= 60] and [r['left'] for r in out['rows']] == [45, -1]
    calls.clear()
    out = NR.c8r_job(a[:7] + (30,))
    assert [c[1] for c in calls] == [10] and len(out['rows']) == 1
    ks = NR.c8k_cases(raw['C8'], raw['C3'])
    assert len(ks) == 22 and all(k[0:3] == (12.0, 10.0, 'paradigm') for k in ks)
    q = next(q for q in raw['C3'] if (q['s'], q['opts'].get('seed')) == (ks[0][3], ks[0][4]['seed']) and len(q['opts']) == 1)
    assert ks[0][5] == [list(x) for x in q['W']]
    calls.clear()
    out = NR.c8k_job(ks[0])
    assert [c[1] for c in calls] == list(range(NR.STEP_C8, NR.NP_C, NR.STEP_C8)) * 2 and all(c[3] for c in calls)
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    W = [[0.9, -0.5, 0.7]] * 3
    C = [[0.1, 0.4, -0.3]] * 3
    al, c = NR.offsets_start(W, C, 2, True)(L, al0, c0, idx)
    w_fit = NR.mobilities(L, NR.EPS, al0, 'paradigm')[0]
    w_new = NR.mobilities(L, NR.EPS, al, 'paradigm')[0]
    assert np.allclose(np.abs(w_new[idx]), [0.9, 0.5, 0.7], rtol=1e-9) and np.allclose(c[idx], [0.1, 0.4, -0.3])
    al2, c2 = NR.offsets_start(W, C, 2, False)(L, al0, c0, idx)
    assert np.array_equal(al2, al) and np.array_equal(c2, c0) and np.allclose(np.sign(w_new[idx]), np.sign(w_fit[idx]))
    j9 = NR.c9a_jobs(raw['C3'], raw['C6'])
    runs = {(int(r['opts']['seed']), float(r['opts'].get('rise', [0, 0.0])[1])): r for r in list(raw['C3']) + list(raw['C6'])
            if r['s'] == 0.02 and 'freqs' not in r['opts'] and r['opts'].get('seed') in (3, 7, 10) and not r['opts'].get('neg')}
    assert len(runs) == 12 and len(j9) == sum((r['left'][0] - 1) if r['left'][0] >= 0 else NR.NP_C - NR.NP_C7 for r in runs.values())
    sd, p, rise, d = j9[17]
    w1 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    assert d == abs(float(np.asarray(runs[(sd, rise)]['W'], float)[p - 1, 0])) - w1
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, **kw: (3, 0.9, 1.5, 0.2, 0.8, 1.4, -0.1, 1.5))
    e = NR.e40_job((12.0, 2.5, 'paradigm', 0, 0.123, 40))
    assert e == dict(m=12.0, ke=2.5, metric='paradigm', k=0, d=0.123, np=40, left=3, maxdev=0.9, t_detach=1.5, x_detach=0.2, c_after=0.8,
                     c_end=1.4, off_min=-0.1, off_max=1.5)


def test_c7_gap_lines_on_the_raw_and_gap_devs():
    """Referee round 85 (D3cf-1): `gap_devs` gives a scan's leaving bands at a cut and the largest distance inside each gap; on the raw's C7
    the band's multi-band windows against the cut are 13, 13, 14, 18, 19, 17, 17 at 0.1 … 0.6, and at 0.3 twelve have a gap of the basin
    (0.014–0.074), six the cut's (0.252–0.300)."""
    import exp_cohesion_leak as E
    r = dict(runs=[dict(d=d, left=(0 if x > NR.LEFT_C else -1), maxdev=x) for d, x in ((0.1, 0.5), (0.2, 0.05), (0.3, 0.6), (0.4, 0.2),
                                                                                     (0.5, 0.29), (0.6, 0.7))])
    assert NR.gap_devs(r) == ([(0.1, 0.1), (0.3, 0.3), (0.6, 0.6)], [0.05, 0.29])
    assert NR.gap_devs(r, 0.25) == ([(0.1, 0.1), (0.3, 0.3), (0.5, 0.6)], [0.05, 0.2]) and NR.gap_devs(r, 0.8) == ([], [])
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    s = '\n'.join(NR.c7_gap_lines(raw['C7']))
    assert "13 at 0.1, 13 at 0.2, 14 at 0.25, 18 at 0.3, 19 at 0.35, 17 at 0.4, 17 at 0.6 of 120" in s
    assert "of the 18 such windows 12 have a gap whose raises all stay within 0.1 of the feature (the largest distance in it 0.014–0.074)" in s
    assert "in the other 6 the gaps' raises move the kernel 0.252–0.300, just under the cut" in s and "and 2 empty" in s
    assert "one tone 13.1, s = 0.06: 0.36–0.48, 0.96–1 (0.289: the cut's)" in s and "one tone 21.3, s = 0.02: 0.2–0.205, 0.27–1 (0.013)" in s
    assert NR.GAP_NEAR == 0.1 and NR.CUTS_C7 == (0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6) and NR.c7_gap_lines(None) == []


def test_c7l_lines_against_the_short_windows():
    """Referee round 85 (D3cf-1): `c7l_lines` sets each long window's leaving set beside the short one's (C7's or C9's), counts the windows
    whose set changes and those with a gap of the basin, and reads C9's draw 3 from period 150 against its runs' offsets."""
    run = lambda pairs: [dict(d=d, left=(0 if x > NR.LEFT_C else -1), maxdev=x) for d, x in pairs]
    C7 = [dict(s=0.02, opts=dict(seed=7, shift=25), runs=run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.6)))),
          dict(s=0.02, opts=dict(seed=7, shift=50), runs=run(((0.1, 0.1), (0.2, 0.5), (0.3, 0.6))))]
    C7L = [dict(s=0.02, opts=dict(seed=7, shift=25), np=20, runs=run(((0.1, 0.5), (0.2, 0.4), (0.3, 0.6)))),
           dict(s=0.02, opts=dict(seed=7, shift=50), np=20, runs=run(((0.1, 0.1), (0.2, 0.5), (0.3, 0.6))))]
    s = '\n'.join(NR.c7l_lines(C7L, C7))
    assert "| 12 | 10 | 0.02 | 7, 25 | 0.1, 0.3 (0.050) | 0.1–0.3 |" in s and "| 12 | 10 | 0.02 | 7, 50 | 0.2–0.3 | 0.2–0.3 |" in s
    assert "differs from the 10-period one in 1 of 2 windows; 0 have more than one band" in s and NR.c7l_lines(None) == []
    two = run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.6), (0.4, 0.2), (0.5, 0.6)))   # two gaps, 0.05 and 0.2: a gap of the basin (L75)
    C7L.append(dict(s=0.02, opts=dict(seed=7, shift=75), np=20, runs=two))
    C7.append(dict(s=0.02, opts=dict(seed=7, shift=75), runs=two))
    s = '\n'.join(NR.c7l_lines(C7L, C7))
    assert "1 have more than one band, 1 of them a gap within 0.1 of the feature" in s


def test_c8_more_lines_on_synthetic_runs():
    """Referee round 85 (D3cf-2, D3cf-3): the nulls' counts (out of sample), the rule on δ alone swept, the replay started further back, the
    dense windows' phases, the five-period tiling and the early alarms, and the kernel-offset windows' recovered misses."""
    P = 35
    wins = [dict(p=p, d=0.01 * p, left=(1 if p == 30 else -1), maxdev=0.1) for p in range(0, NR.NP_C, NR.STEP_C8)]
    C8 = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=[P, -1, -1], maxdev=[0.9, 0, 0], wins=wins),
          dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=2), left=[-1, -1, -1], maxdev=[0.1, 0, 0],
               wins=[dict(x, left=-1) for x in wins])]
    z = lambda flag: [dict(p=x['p'], left=(1 if x['p'] in flag else -1), maxdev=0.1) for x in wins]
    C8N = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), zero=z(()), off={str(o): z((10,)) for o in NR.C8N_OFF}),
           dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=2), zero=z((50,)), off={str(o): z(()) for o in NR.C8N_OFF})]
    C8R = [dict(left=P, rows=[dict(h=20, q=15, left=35, maxdev=0.9), dict(h=40, q=-5, left=-1, maxdev=0.1)])]
    dense = [dict(p=p, d=0.0, left=(1 if p in (29, 30, 31, 32) else -1), maxdev=0.1) for p in range(0, P)]
    C8D = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), wins=dense),
           dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=2), wins=[dict(p=p, d=0.0, left=(1 if p == 7 else -1), maxdev=0.1)
                                                                                    for p in range(0, 291)])]
    EDGE = [dict(m=12.0, ke=10.0, metric='paradigm', k=0, lo=0.2, hi=0.3)]
    s = '\n'.join(NR.c8_more_lines(C8, C8N, C8R, C8D, None, None, EDGE))
    assert "δ = 0, its coming noise: the leave's window flagged in 0 of 1, 1 of" in s and "1 false alarms" in s
    assert "its δ, the draw's noise 75 periods later: the leave's window flagged in 0 of 1, 1 of" in s
    assert "20 periods before: of 1 runs 1 leave within one period of the run (out of sample 0 of 0, in sample 1 of 1), 0 earlier, 0 later, 0 never" in s
    assert "40 periods before: of 1 runs 0 leave within one period of the run (out of sample 0 of 0, in sample 0 of 1), 0 earlier, 0 later, 1 never" in s
    assert "the leave's window is flagged in 0–1 of 1" in s and "missed at every phase: none" in s
    assert "of the 1 leaves fall in a flagged window" in s and ("(or whose run never left), 1 are flagged — m = 12, k_e = 10, s = 0.02: 1: 0 of them the replay one period early, 1 isolated starts — a lead other than one period — "
                                                       "(m = 12, k_e = 10, s = 0.02: 1), 1 in runs that never left, in the window's period 1. Of the other dense "
                                                       "windows that cover an isolated alarm's leave period (7), 0 repeat it (leaving within one period of it; referee round 94); C8R's loop, started before such "
                                                       "a period, covers 0 of the 1 and passes 0 of them without leaving.") in s
    sw = NR.c8_sweep(C8, EDGE)
    assert sw[0] == (0.0, 1, 1, 3 + len(wins), 3 + len(wins)) and len(sw) == 201      # f = 0 flags every window
    assert sw[-1] == (2.0, 0, 1, sum(1 for x in wins if x['d'] >= 2.0 * 0.25), 3 + len(wins))  # δ ≥ 2e₀, e₀ = 0.25: only the never-left run's late windows


def test_c9a_lines_and_e40_line():
    """Referee round 85 (D3cf-4, D3cf-5): `c9a_lines` counts the windows that end before their run's leave and are flagged, and the windows
    that contain the leave flagged per run; `e40_line` the ejection's ends over 10 and 40 periods."""
    C3 = [dict(s=0.02, opts=dict(seed=3), left=[30, -1, -1])]
    C6 = [dict(s=0.02, opts=dict(seed=3, rise=[0, 0.05]), left=[25, -1, -1])]
    C9A = ([dict(seed=3, p=p, rise=0.0, d=0.0, left=(5 if p == 12 else -1), maxdev=0.1) for p in range(1, 30)]
           + [dict(seed=3, p=p, rise=0.05, d=0.0, left=(1 if p >= 18 else -1), maxdev=0.1) for p in range(1, 25)])
    s = '\n'.join(NR.c9a_lines(C9A, C3, C6))
    assert ("of the 33 windows that end before their run's leave 1 are flagged — each window with the periods by which its loop leaves "
            "before the run: draw 3 raised by 0 from 12 (13) — 0 of them within one period;") in s
    C9A[2] = dict(C9A[2], left=1)                                        # a second flagged window, its lead the larger: each lead with its window
    s = '\n'.join(NR.c9a_lines(C9A, C3, C6))
    assert "of the 33 windows that end before their run's leave 2 are flagged" in s
    assert "before the run: draw 3 raised by 0 from 3 (26), draw 3 raised by 0 from 12 (13) — 0 of them within one period;" in s
    assert ("0 of 10 starts (draw 3, un-raised); 7 of 10 starts (draw 3 raised by 0.05) — the least 0 of 10 (draw 3, un-raised), the "
            "other 1 runs at 7–7 of 10–10 starts.") in s
    assert NR.c9a_least({(3, 0.0): [4, 10], (7, 0.0): [9, 10], (7, 0.05): [10, 10]}) == (
        " — the least 4 of 10 (draw 3, un-raised), the other 2 runs at 9–10 of 10–10 starts") and NR.c9a_least({(3, 0.0): [1, 2]}) == ''
    assert "4 of 10 starts (draw 3, un-raised)" not in s and "7 of 10 starts (draw 3 raised by 0.05)" in s and "0 of 10 starts (draw 3, un-raised)" in s
    E40 = [dict(ke=2.5, d=0.123, np=10, left=3, c_end=1.383, t_detach=2.918), dict(ke=2.5, d=0.123, np=40, left=3, c_end=1.403, t_detach=2.918)]
    s = '\n'.join(NR.e40_line(E40))
    assert "over 10 periods its kernel leaves at period 3, ends at +1.383, detached for good" in s and "outside the orbit's range [−1, 1]" in s
    assert NR.e40_line(None) == [] and NR.c9a_lines(None) == []


def test_edge_fit_robust_errors_and_curvature():
    """Referee round 85 (D3cf-5): `edge_fit`'s HC3 errors are the sandwich (XᵀX)⁻¹Xᵀdiag(e_i²/(1 − h_i)²)X(XᵀX)⁻¹, its jackknife errors
    ((n − 1)/n)Σ(b_(i) − b̄)², and its curvature test the F of (log k_e)², (log m)² added, against direct computations."""
    rng = np.random.default_rng(5)
    EDGE = [dict(m=m, ke=ke, metric='paradigm', k=0, lo=0.0, hi=2 * float(np.exp(-1.0 + 0.5 * np.log(ke) - 0.5 * np.log(m)
                                                                                  + 0.1 * np.log(ke) ** 2 + rng.normal(0, 0.01))))
            for m in (3.0, 6.0, 12.0) for ke in (5.0, 10.0, 20.0, 40.0)]
    fit = NR.edge_fit(EDGE)
    X = np.array([[1.0, np.log(e['ke']), np.log(e['m'])] for e in fit['rows']])
    y = np.array([np.log(0.5 * (e['lo'] + e['hi'])) for e in fit['rows']])
    b = np.linalg.solve(X.T @ X, X.T @ y)
    e = y - X @ b
    A = np.linalg.inv(X.T @ X)
    h = np.array([x @ A @ x for x in X])
    V = A @ sum(np.outer(x, x) * (ei / (1 - hi)) ** 2 for x, ei, hi in zip(X, e, h)) @ A
    assert np.allclose(fit['se_hc3'], np.sqrt(np.diag(V)), rtol=1e-10)
    bs = np.array([np.linalg.solve(np.delete(X, i, 0).T @ np.delete(X, i, 0), np.delete(X, i, 0).T @ np.delete(y, i)) for i in range(len(y))])
    assert np.allclose(fit['se_jack'], np.sqrt((len(y) - 1) / len(y) * ((bs - bs.mean(0)) ** 2).sum(0)), rtol=1e-10)
    X2 = np.c_[X, X[:, 1] ** 2, X[:, 2] ** 2]
    e2 = y - X2 @ np.linalg.lstsq(X2, y, rcond=None)[0]
    F = ((e @ e - e2 @ e2) / 2) / (e2 @ e2 / (len(y) - 5))
    assert np.isclose(fit['F'], F, rtol=1e-10) and fit['F_df'] == (2, len(y) - 5) and fit['F_p'] < 1e-4


def test_the_case_lists_against_their_stated_contents():
    """Referee round 85 (D3cf-7: S71, S72, S73, C8C, EO3; CL4): SCAN_C7 is the grid 0–0.2 by 0.0025, 0.2–0.4 by 0.005, 0.4–1 by 0.02 (151
    raises); C7_CASES, C7L_CASES, C8_CASES and EDGE_CASES are, one to one and in order, the raw's C7, C7L, C8 and EDGE records; C7_CASES
    holds every tone of TONES_C7 at every level of S_C7 + S_C7_HI under +n and −n, the band's windows from every SHIFT_C7 at S_C7, and
    C7_MORE's configurations; C7L_CASES holds C7's 18 band windows with more than one leaving band, at twice the length; C8_CASES the
    draws 1–11 at s = 0.02, 0.03; EDGE_OUT the rows outside the window, kernel 3's among them."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    g = NR.SCAN_C7
    assert len(g) == 151 and g[:3] == (0.0, 0.0025, 0.005) and g[79:82] == (0.1975, 0.2, 0.205) and g[119:122] == (0.395, 0.4, 0.42)
    assert g[-1] == 1.0 and all(b > a for a, b in zip(g, g[1:]))
    no = lambda o: {k: (tuple(v) if isinstance(v, (list, tuple, np.ndarray)) else v) for k, v in dict(o).items()}
    cfg = lambda r: (r['m'], r['ke'], r['metric'])
    c7 = [(no(c[0]), c[1], c[2] if len(c) > 2 else NR.C7_CFG) for c in NR.C7_CASES]
    assert c7 == [(no(r['opts']), r['s'], cfg(r)) for r in raw['C7']]
    tones = [(o['freqs'][0], s_, bool(o.get('neg'))) for o, s_, _ in c7 if 'freqs' in o]
    assert sorted(tones) == sorted((f, s_, ng) for f in NR.TONES_C7 for s_ in NR.S_C7 + NR.S_C7_HI for ng in (False, True))
    assert NR.S_C7_HI == (0.04, 0.06, 0.08) and NR.TONES_C7 == (10.0, 13.1, 15.0, 21.3) and NR.SHIFT_C7 == tuple(range(0, 200, 25))
    assert [cf for cf, _ in NR.C7_MORE] == [(12.0, 20.0, 'paradigm'), (4.0, 10.0, 'paradigm'), (12.0, 10.0, 'wfr')]
    assert [lv for _, lv in NR.C7_MORE] == [(0.02, 0.04, 0.06, 0.08), (0.04, 0.06, 0.08), (0.08, 0.16, 0.24)]
    c7l = [(no(c[0]), c[1], c[2] if c[2] is not None else NR.C7_CFG, c[3]) for c in NR.C7L_CASES]
    assert c7l == [(no(r['opts']), r['s'], cfg(r), r['np']) for r in raw['C7L']] and all(c[3] == 2 * NR.NP_C7 for c in c7l)
    band = [r for r in raw['C7'] if 'freqs' not in r['opts'] and r['s'] > 0]
    multi = {(no(r['opts'])['shift'], r['s'], cfg(r)) for r in band if len(NR.gap_devs(r)[0]) > 1}
    assert len(multi) == 18 and multi <= {(o['shift'], s_, cf) for o, s_, cf, _ in c7l if o['seed'] == 7}
    assert (dict(seed=3, shift=150), 0.02) == tuple(NR.C7L_CASES[-1][:2])
    c8 = [(c[0], c[1], c[2], c[3], no(c[4])) for c in NR.C8_CASES]
    assert c8 == [(r['m'], r['ke'], r['metric'], r['s'], no(r['opts'])) for r in raw['C8']]
    assert [o['seed'] for m, ke, mt, s_, o in c8 if (m, ke, mt) == (12.0, 10.0, 'paradigm') and s_ == 0.02] == list(range(1, 12))
    ed = [(c[0], c[1], c[2], c[3], c[4] if len(c) > 4 else 0.2) for c in NR.EDGE_CASES]
    assert ed == [(r['m'], r['ke'], r['metric'], r['k'], r['om']) for r in raw['EDGE']]
    assert (12.0, 2.5, 'paradigm', 0) in NR.EDGE_OUT and sum(c[3] == 2 for c in NR.EDGE_OUT) == 3 and NR.TRACE_K == 1024


def test_the_trace_on_a_planted_crossing(monkeypatch):
    """Referee round 85 (D3cf-7: ER2, ER4–ER7, EJ2): `edge_run`'s anatomy on a fake loop whose kernel-1 centre has a planted path — at its
    feature, then 0.045 off (within 0.05), then from sample n1 on 0.2 … 0.29 off for good: the detachment at ts[n1] (the first sample
    after the last within 0.05), x_m = sin(ω·ts[n1]) at the run's own ω (0.4), the centre at the end of the detachment's period and at the
    end, and the least and largest offsets; a path back within 0.05 at the end gives no detachment; `edge_job` traces at its own ω."""
    from src.core import closed_loop as CL  # noqa: F401  (the real simulate is not called)
    kk, NP, om = NR.TRACE_K, 3, 0.4
    n0, n1 = int(1.2 * kk), int(1.6 * kk)
    N = NP * kk + 1
    L = NR.law(12.0)
    al0, c0, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    off = np.zeros(N)
    off[n0:n1] = 0.045
    off[n1:] = 0.2 + 0.09 * np.arange(N - n1) / (N - n1)
    back = off.copy()
    back[-5:] = 0.0

    def fake(path):
        def sim(L_, al, c, astar, ts, noise, metric, ke, om=None):
            C = np.tile(np.asarray(c, float), (len(ts), 1))
            C[:, idx[0]] = NR.CT.CSTAR[0] + path
            return dict(c=C)
        return sim
    monkeypatch.setattr(NR, 'simulate', fake(off))
    out = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05, NP=NP, trace=True, om=om)
    T = 2 * np.pi / om
    ts = np.linspace(0, NP * T, N)
    p = int(ts[n1] // T)
    assert p == 1 and out[2] == float(ts[n1] / T) and out[3] == float(np.sin(om * ts[n1]))
    assert out[4] == float(NR.CT.CSTAR[0] + off[(p + 1) * kk]) and out[4] != float(NR.CT.CSTAR[0] + off[(p + 2) * kk])
    assert out[5:] == (float(NR.CT.CSTAR[0] + off[-1]), 0.0, float(off.max()))
    monkeypatch.setattr(NR, 'simulate', fake(back))
    out = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05, NP=NP, trace=True, om=om)
    assert out[2:5] == (None, None, None) and out[5] == float(NR.CT.CSTAR[0])
    calls = []
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, om=None, **kw:
                        calls.append((NP, trace, om)) or ((2 if d >= 0.1 else -1, 0.5 if d >= 0.1 else 0.01) + ((0.5, 0.1, 0.2, 0.3, 0.0, 0.4) if trace else ())))
    NR.edge_job((12.0, 10.0, 'paradigm', 0, 0.4))
    assert [c for c in calls if c[1]] == [(3, True, 0.4), (3, True, 0.4)] and all(c[2] == 0.4 for c in calls)


def test_win_run_passes_its_window_and_start(monkeypatch):
    """Referee round 85 (WR2, WR3): `win_run` runs `edge_run` on kernel 1 with δ, the window's length NP, the run's noise from period
    shift (`c_noise` with the shift) and the start it is given."""
    calls = []
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=None, noise=None, start=None, **kw: calls.append((m, ke, mt, k, d, NP, noise, start)) or (3, 0.4))
    monkeypatch.setattr(NR, 'c_noise', lambda s, opts: ('noise', s, tuple(sorted(opts.items()))))
    st = lambda L, al0, c0, idx: (al0, c0)
    assert NR.win_run(12.0, 10.0, 'paradigm', 0.02, dict(seed=4), 0.07, 30, NP=7, start=st) == (3, 0.4)
    assert NR.win_run(12.0, 10.0, 'paradigm', 0.02, dict(seed=4), 0.07, 30) == (3, 0.4)
    assert calls == [(12.0, 10.0, 'paradigm', 0, 0.07, 7, ('noise', 0.02, (('seed', 4), ('shift', 30))), st),
                     (12.0, 10.0, 'paradigm', 0, 0.07, NR.NP_C7, ('noise', 0.02, (('seed', 4), ('shift', 30))), None)]


def test_round85_jobs_on_planted_inputs(monkeypatch):
    """Referee round 85 (RR2, PR1, OS3, OS4, KJ1, KC1, 9J7, PE1, E4J1): `c8r_job` starts at q = P − h down to q = 1; `part_C8R` takes only
    the runs that left; `offsets_start` raises each formed kernel by |w_k| over the period before p minus its fitted value and, with
    centres, puts each centre at its period mean there; `c8k_job` runs the coefficients' windows, then the centres'; `c8k_cases` pairs
    every C8 run with C3's run of its own seed; `c9a_jobs` carries each run's rise; `part_E40` runs k_e = 2.5's hi over 10 and 40 periods
    and `e40_job` passes NP."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    calls = []
    monkeypatch.setattr(NR, 'win_run', lambda m, ke, mt, s, opts, d, shift, NP=NR.NP_C7, start=None: calls.append((d, shift, NP, start)) or (-1, 0.1))
    w = [0.3 + 0.01 * i for i in range(NR.NP_C)]
    out = NR.c8r_job((12.0, 10.0, 'paradigm', 0.02, dict(seed=1), w, 0.28, 21))
    assert [(r['h'], r['q']) for r in out['rows']] == [(20, 1)] and calls == [(NR.c8_delta(w, 0.28, 1), 1, 22, None)]
    got = []
    monkeypatch.setattr(NR, '_pool_map', lambda fn, jobs: got.append((fn, list(jobs))) or [])
    C8 = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), w=w, w0=0.28, left=[50, -1, -1]),
          dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=2), w=w, w0=0.28, left=[-1, -1, -1])]
    NR.part_C8R(C8)
    assert got[-1][0] is NR.c8r_job and [j[4]['seed'] for j in got[-1][1]] == [1]
    EDGE = [dict(m=12.0, ke=2.5, metric='paradigm', k=0, om=0.2, lo=0.05875, hi=0.06)]
    NR.part_E40(EDGE)
    assert got[-1][0] is NR.e40_job and got[-1][1] == [(12.0, 2.5, 'paradigm', 0, 0.06, 10), (12.0, 2.5, 'paradigm', 0, 0.06, 40)]
    ecalls = []
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=None, trace=False, **kw: ecalls.append((d, NP, trace)) or tuple(range(8)))
    r = NR.e40_job((12.0, 2.5, 'paradigm', 0, 0.06, 40))
    assert ecalls == [(0.06, 40, True)] and (r['np'], r['left'], r['off_max']) == (40, 0, 7)
    L = NR.law(12.0)
    al0, c0, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    W = [[0.9 + 0.01 * p, -0.6 - 0.01 * p, 0.75 + 0.02 * p] for p in range(20)]
    Cc = [[-0.52 + 0.001 * p, 0.07 - 0.001 * p, 0.61 + 0.002 * p] for p in range(20)]
    rs = []
    monkeypatch.setattr(NR, 'raised', lambda L_, al, idx_, k, delta: rs.append((k, delta)) or al)
    w0 = [abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[k]])) for k in range(len(idx))]
    al, c = NR.offsets_start(W, Cc, 7, True)(L, al0, c0, idx)
    assert rs == [(k, abs(W[6][k]) - w0[k]) for k in range(len(idx))] and [c[idx[k]] for k in range(len(idx))] == [Cc[6][k] for k in range(len(idx))]
    al, c = NR.offsets_start(W, Cc, 7, False)(L, al0, c0, idx)
    assert np.array_equal(c, c0)
    starts = []
    monkeypatch.setattr(NR, 'offsets_start', lambda W_, C_, p, cen: starts.append((p, cen)) or 'start')
    out = NR.c8k_job((12.0, 10.0, 'paradigm', 0.02, dict(seed=1), W, Cc))
    ps = list(range(NR.STEP_C8, NR.NP_C, NR.STEP_C8))
    assert starts == [(p, False) for p in ps] + [(p, True) for p in ps] and [x['p'] for x in out['w']] == ps == [x['p'] for x in out['wc']]
    ks = NR.c8k_cases(raw['C8'], raw['C3'])
    for k_ in ks:
        q = next(q for q in raw['C3'] if (q['m'], q['ke'], q['s'], q['opts'].get('seed')) == (k_[0], k_[1], k_[3], k_[4]['seed']) and len(q['opts']) == 1)
        assert k_[5] == [list(x) for x in q['W']] and k_[6] == [list(x) for x in q['C']], k_[4]
    assert len({(k_[3], k_[4]['seed']) for k_ in ks}) == len(ks)
    C3 = [dict(s=0.02, opts=dict(seed=3), left=[5, -1, -1], W=[[0.95, -0.6, 0.75]] * 10)]
    C6 = [dict(s=0.02, opts=dict(seed=3, rise=[0, 0.05]), left=[4, -1, -1], W=[[0.99, -0.6, 0.75]] * 10)]
    jobs = NR.c9a_jobs(C3, C6)
    assert [(j[0], j[1], j[2]) for j in jobs] == [(3, p, 0.0) for p in range(1, 5)] + [(3, p, 0.05) for p in range(1, 4)]


def test_c8_more_lines_boundaries_and_main_c9a_c8k(monkeypatch, tmp_path):
    """Referee round 85 (MD6, MK2, MK4, RW7, RW8): the five-period tiling counts a leave in a window's last period; C8K counts a miss as
    recovered only when the offsets' windows hit the leave (an early flag is a new early flag), and reads C3's offsets at the leave's
    window — the period before it — also when the leave falls on a multiple of STEP_C8; `main --only C9A,C8K` passes C3 and C6 to C9A and
    C8 and C3 to C8K."""
    P = 35
    C8 = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=[P, -1, -1], maxdev=[0.9, 0, 0],
               wins=[dict(p=p, d=0.0, left=-1, maxdev=0.1) for p in range(0, NR.NP_C, NR.STEP_C8)])]
    C8D = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1),
                wins=[dict(p=p, d=0.0, left=(1 if p == P - NR.NP_C7 else -1), maxdev=0.1) for p in range(0, P)])]
    s = '\n'.join(NR.c8_more_lines(C8, None, None, C8D))
    assert "With a window every 5 periods, 0–1 of the 1 leaves fall in a flagged window" in s
    P = 40
    C8[0]['left'] = [P, -1, -1]
    flag = lambda ps: [dict(p=p, left=(1 if p in ps else -1), maxdev=0.1) for p in range(NR.STEP_C8, NR.NP_C, NR.STEP_C8)]
    C8K = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), w=flag((30,)), wc=flag((20,)))]
    W = [[0.9, 0.6 + 0.001 * p, 0.75 + 0.002 * p] for p in range(NR.NP_C)]
    C3 = [dict(s=0.02, opts=dict(seed=1), left=[P, -1, -1], W=W)]
    s = '\n'.join(NR.c8_more_lines(C8, None, None, None, C8K, C3))
    assert ("the coefficients: misses recovered s = 0.02, draw 1, new early flags or false alarms none; with the centres: misses recovered "
            "none, new early flags or false alarms s = 0.02, draw 1") in s
    L = NR.law(NR.C7_CFG[0])
    al0, _, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    fit = [abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][i])) for i in idx]
    o2, o3 = abs(W[29][1]) - fit[1], abs(W[29][2]) - fit[2]
    assert f"s = 0.02: kernel 2 {o2:+.3f} … {o2:+.3f}, kernel 3 {o3:+.3f} … {o3:+.3f}" in s
    old = dict(C3=['c3'], C6=['c6'], C8=['c8'])
    saved, got = [], []
    monkeypatch.setattr(NR.E, 'load', lambda p: dict(old))
    monkeypatch.setattr(NR.E, 'save', lambda p, **kw: saved.append(kw))
    monkeypatch.setattr(NR, 'OUT', str(tmp_path))
    monkeypatch.setattr(NR, 'report', lambda **kw: '')
    monkeypatch.setattr(NR, 'part_C9A', lambda C3, C6: got.append(('C9A', C3, C6)) or ['a'])
    monkeypatch.setattr(NR, 'part_C8K', lambda C8, C3: got.append(('C8K', C8, C3)) or ['k'])
    monkeypatch.setattr(sys, 'argv', ['exp_noise_resonance.py', '--only', 'C9A,C8K'])
    NR.main()
    assert sorted(got) == [('C8K', ['c8'], ['c3']), ('C9A', ['c3'], ['c6'])] and saved[-1]['C9A'] == ['a'] and saved[-1]['C8K'] == ['k']


def test_round89_alarms_by_kind_and_the_survivors(monkeypatch, tmp_path):
    """Referee round 89 (D3cj-1, D3cj-5: R03, R05, R11, R16, R21, R22, R31, D3cj-6): `c8d_alarms` splits the dense windows' early flags
    into the replay one period early and isolated starts, with the neighbours covering an isolated alarm's leave period (and whether
    those leaving by then are alarms themselves) and C8R's rows started before it; C8D's phases count out-of-sample runs only; C8R's lag
    of −1 is within one period; C9A's lead of exactly 1 is within one period, and its runs are those it ran; a one-band window to the
    grid's end is not "ending below it"; `c8r_job` starts at all three h when P > 81; `main` passes this run's C8 to C8D and EDGE to E40."""
    P = 40
    mk = lambda seed, P_, fl: (dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=seed), left=[P_, -1, -1], maxdev=[0.9, 0, 0],
                                    wins=[dict(p=p, d=0.0, left=-1, maxdev=0.1) for p in range(0, NR.NP_C, NR.STEP_C8)]),
                               dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=seed),
                                    wins=[dict(p=p, d=0.0, left=fl.get(p, -1), maxdev=0.1) for p in range(0, P_)]))
    r1, d1 = mk(1, P, {5: 2, 6: 1, 29: 10})
    r7, d7 = mk(7, P, {30: 1})                                           # in sample: its leave's window flagged
    C8R = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=P,
                rows=[dict(h=20, q=20, left=P - 1, maxdev=0.9), dict(h=37, q=3, left=P, maxdev=0.9)])]
    sp = NR.c8d_alarms([r1, r7], [d1, d7], C8R)['alarms']
    assert [(z['p'], z['kind'], z['lead']) for z in sp] == [(5, 'isolated', 33), (6, 'isolated', 33), (29, 'early', 1)]
    assert [(z['nb'], z['rep'], z['c8r']) for z in sp[:2]] == [(6, 1, True), (6, 1, True)]
    s = '\n'.join(NR.c8_more_lines([r1, r7], None, C8R, [d1, d7]))
    assert ("3 are flagged — m = 12, k_e = 10, s = 0.02: 3: 1 of them the replay one period early (the window's loop leaving in its last "
            "period, the run in the next; m = 12, k_e = 10, s = 0.02: 1), 2 isolated starts — a lead other than one period — (m = 12, k_e = 10, "
            "s = 0.02: 2) whose loops leave "
            "33–33 periods before their runs, in the window's period 1, 2. Of the other dense windows that cover an isolated alarm's leave "
            "period (12), 2 repeat it (leaving within one period of it; referee round 94); C8R's loop, started before such a period, covers 2 of the 2 "
            "and passes 2 of them without leaving (implied by C8R's aggregate: none of its loops leaves more than one period before its run; "
            "referee rounds 102, 105).") in s
    C8Rb = [dict(C8R[0], rows=[dict(h=20, q=20, left=P - 1, maxdev=0.9), dict(h=37, q=3, left=6, maxdev=0.9)])]   # round 102: a loop leaving
    assert [z['c8r'] for z in NR.c8d_alarms([r1, r7], [d1, d7], C8Rb)['alarms'][:2]] == [False, False]            # before the alarm's period
    s_b = '\n'.join(NR.c8_more_lines([r1, r7], None, C8Rb, [d1, d7]))
    assert "covers 2 of the 2 and passes 0 of them without leaving" in s_b and "implied by C8R's aggregate" not in s_b
    C8Re = [dict(C8R[0], rows=[dict(h=20, q=20, left=P - 1, maxdev=0.9), dict(h=37, q=3, left=P, maxdev=0.9)]),
            dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=9), left=P, rows=[dict(h=20, q=20, left=P - 5, maxdev=0.9)])]
    s_e = '\n'.join(NR.c8_more_lines([r1, r7], None, C8Re, [d1, d7]))       # round 105: another loop leaving early breaks the aggregate
    assert "passes 2 of them without leaving" in s_e and "implied by C8R's aggregate" not in s_e
    C8R2 = [C8Re[0], dict(C8Re[1], rows=[dict(h=20, q=20, left=P - 2, maxdev=0.9)])]   # round 111's R07: a loop exactly 2 periods early
    assert "implied by C8R's aggregate" not in '\n'.join(NR.c8_more_lines([r1, r7], None, C8R2, [d1, d7]))
    C8R1 = [C8Re[0], dict(C8Re[1], rows=[dict(h=20, q=20, left=P - 1, maxdev=0.9)])]   # … and one period early: the aggregate holds
    assert "implied by C8R's aggregate" in '\n'.join(NR.c8_more_lines([r1, r7], None, C8R1, [d1, d7]))
    assert "the leave's window is flagged in 0–0 of 1 and 0–1 earlier windows are flagged" in s          # the in-sample run not counted (R31)
    assert "20 periods before: of 1 runs 1 leave within one period of the run (out of sample 1 of 1, in sample 0 of 0), 0 earlier" in s   # a lag of −1 (R05)
    C3 = [dict(s=0.02, opts=dict(seed=3), left=[30, -1, -1]), dict(s=0.02, opts=dict(seed=4), left=[30, -1, -1])]
    C9A = [dict(seed=3, p=p, rise=0.0, d=0.0, left=(9 if p == 10 else -1), maxdev=0.1) for p in range(1, 30)]
    s = '\n'.join(NR.c9a_lines(C9A, C3, []))
    assert "(C9A, 29 windows of 1 runs)" in s and "draw 3 raised by 0 from 10 (11)" in s
    C9A[5] = dict(C9A[5], left=NR.NP_C7)                                  # p = 6: its loop leaves at 16, the run at 30 — 14
    C9A[18] = dict(C9A[18], p=19, left=NR.NP_C7)                          # p = 19: leaves at 29, one period before the run (R11)
    s = '\n'.join(NR.c9a_lines(C9A, C3, []))
    assert "— 1 of them within one period" in s
    band = [dict(s=0.02, opts=dict(seed=7, shift=25), runs=[dict(d=d, left=(0 if d >= 0.5 else -1), maxdev=(0.5 if d >= 0.5 else 0.01))
                                                            for d in NR.SCAN_C7])]
    s = '\n'.join(NR.c7_gap_lines(band))
    assert "1 are one band — 1 of them up to the grid's end." in s and "ending below it" not in s        # R03
    calls = []
    monkeypatch.setattr(NR, 'win_run', lambda m, ke, mt, s_, opts, d, shift, NP=NR.NP_C7, start=None: calls.append((shift, NP)) or (-1, 0.1))
    w = [0.3 + 0.01 * i for i in range(NR.NP_C)]
    out = NR.c8r_job((12.0, 10.0, 'paradigm', 0.02, dict(seed=1), w, 0.28, 100))
    assert [(r['h'], r['q']) for r in out['rows']] == [(20, 80), (40, 60), (80, 20)]                  # R16
    old = dict(C8=['stale8'], EDGE=['staleedge'])
    saved, got = [], []
    monkeypatch.setattr(NR.E, 'load', lambda p: dict(old))
    monkeypatch.setattr(NR.E, 'save', lambda p, **kw: saved.append(kw))
    monkeypatch.setattr(NR, 'OUT', str(tmp_path))
    monkeypatch.setattr(NR, 'report', lambda **kw: '')
    monkeypatch.setattr(NR, 'part_C8', lambda: ['new8'])
    monkeypatch.setattr(NR, 'part_EDGE', lambda: ['newedge'])
    monkeypatch.setattr(NR, 'part_C8D', lambda C8: got.append(('C8D', C8)) or ['d'])
    monkeypatch.setattr(NR, 'part_E40', lambda EDGE: got.append(('E40', EDGE)) or ['e'])
    monkeypatch.setattr(sys, 'argv', ['exp_noise_resonance.py', '--only', 'EDGE,C8,C8D,E40'])
    NR.main()
    assert sorted(got) == [('C8D', ['new8']), ('E40', ['newedge'])]                                    # R21, R22


def test_round89_c7s_and_c7l(monkeypatch):
    """Referee round 89 (D3cj-1, D3cj-4): C7S_CASES are C7's windows with a gap of the basin and C9's draw 3 from 150; `c7s_job` moves the
    start by d and records it; `part_C7S` runs every case at both offsets; `c7s_lines` counts the windows whose neighbours keep a gap of the
    basin; `c7l_lines` counts the long sets that contain the short ones and the gaps of the basin that survive or close; C8K prints the
    recovered misses' own offsets and kernel 3's range over the draws whose kernel 3 has not left."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    cfg = lambda r: (r.get('m', NR.C7_CFG[0]), r.get('ke', NR.C7_CFG[1]), r.get('metric', NR.C7_CFG[2]))
    real = [r for r in raw['C7'] if 'freqs' not in r['opts'] and r['s'] > 0 and NR.gap_devs(r)[1] and min(NR.gap_devs(r)[1]) < NR.GAP_NEAR]
    want = [(None if cfg(r) == NR.C7_CFG else cfg(r), r['s'], r['opts']['seed'], r['opts'].get('shift', 0)) for r in real]
    assert sorted(want, key=str) == sorted(NR.C7S_CASES[:-1], key=str) and NR.C7S_CASES[-1] == (None, 0.02, 3, 150) and NR.C7S_D == (-1, 1)
    assert len(real) == 12
    got = []
    monkeypatch.setattr(NR, 'c7_job', lambda a: got.append(a) or dict(runs=[], s=a[1], opts=dict(a[0])))
    out = NR.c7s_job(((4.0, 10.0, 'paradigm'), 0.08, 7, 100, -1))
    assert got[-1] == (dict(seed=7, shift=99), 0.08, (4.0, 10.0, 'paradigm')) and (out['base'], out['dshift']) == (100, -1)
    jobs = []
    monkeypatch.setattr(NR, '_pool_map', lambda fn, js: jobs.append((fn, list(js))) or [])
    NR.part_C7S()
    assert jobs[-1][0] is NR.c7s_job and jobs[-1][1] == [c + (d,) for c in NR.C7S_CASES for d in (-1, 1)]
    run = lambda pairs: [dict(d=d, left=(0 if x > NR.LEFT_C else -1), maxdev=x) for d, x in pairs]
    gap, flat = run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.6))), run(((0.1, 0.5), (0.2, 0.5), (0.3, 0.6)))
    C7 = [dict(s=0.02, opts=dict(seed=7, shift=25), runs=gap), dict(s=0.02, opts=dict(seed=7, shift=50), runs=gap)]
    C7S = [dict(s=0.02, opts=dict(seed=7, shift=24), runs=gap, base=25, dshift=-1), dict(s=0.02, opts=dict(seed=7, shift=26), runs=flat, base=25, dshift=1),
           dict(s=0.02, opts=dict(seed=7, shift=49), runs=gap, base=50, dshift=-1), dict(s=0.02, opts=dict(seed=7, shift=51), runs=gap, base=50, dshift=1)]
    s = '\n'.join(NR.c7s_lines(C7S, C7))
    assert "| 12 | 10 | 0.02 | 7, 25 | 0.1, 0.3 (0.050) | 0.1, 0.3 (0.050) | 0.1–0.3 | 1 of 2 |" in s
    assert "survives at both neighbours in 1 of the 2 windows, at one in 1, at neither in 0." in s and NR.c7s_lines(None) == []
    C7L = [dict(s=0.02, opts=dict(seed=7, shift=25), np=20, runs=flat), dict(s=0.02, opts=dict(seed=7, shift=50), np=20, runs=run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.1))))]
    s = '\n'.join(NR.c7l_lines(C7L, C7))
    assert ("By construction the 20-period run's first 10 periods are the 10-period run, so the set can only grow (1 of the 2: a check of "
            "determinism, not a measurement; referee round 94); of the 2 windows with a gap of the basin over 10 periods, a gap of the basin "
            "survives over 20 in 0 and closes in m = 12, k_e = 10, s = 0.02, draw 7 from 25; m = 12, k_e = 10, s = 0.02, draw 7 from 50.") in s
    assert NR.KEYS[-3:] == ('C7S', 'C8A', 'C8F')                        # round 134: C8F after C8A, whose 'cold_c' it reads
    P = 40
    C8 = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.03, opts=dict(seed=s_), left=[P, -1, -1], maxdev=[0.9, 0, 0],
               wins=[dict(p=p, d=0.0, left=-1, maxdev=0.1) for p in range(0, NR.NP_C, NR.STEP_C8)]) for s_ in (1, 2)]
    flag = lambda ps: [dict(p=p, left=(1 if p in ps else -1), maxdev=0.1) for p in range(NR.STEP_C8, NR.NP_C, NR.STEP_C8)]
    C8K = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.03, opts=dict(seed=1), w=flag((30,)), wc=flag(()))]
    W = lambda a: [[0.9, 0.6 + a * p, 0.75 + 2 * a * p] for p in range(NR.NP_C)]
    C3 = [dict(s=0.03, opts=dict(seed=1), left=[P, -1, -1], W=W(0.001)), dict(s=0.03, opts=dict(seed=2), left=[P, -1, 12], W=W(0.003))]
    s = '\n'.join(NR.c8_more_lines(C8, None, None, None, C8K, C3))
    L = NR.law(NR.C7_CFG[0])
    al0, _, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    fit = [abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][i])) for i in idx]
    o = lambda a, k: abs(W(a)[29][k]) - fit[k]
    assert (f"kernel 3 {o(0.001, 2):+.3f} … {o(0.001, 2):+.3f} over the 1 of 2 draws whose kernel 3 is still at its feature (draw 2's, gone, "
            f"{o(0.003, 2):+.3f}); the recovered misses' own: s = 0.03, draw 1: kernel 2 {o(0.001, 1):+.3f}, kernel 3 {o(0.001, 2):+.3f}.") in s


def test_round89_md_is_the_report_and_planted_edges():
    """The committed md is `report` of the committed raw, byte for byte (a report-only change fails here: A05, C01, E01, K01, L01); and the
    round's classifications on planted edges (A04, A06, L02, S03): a neighbour that leaves only after the alarm's period does not repeat it;
    C8R passes an alarm only if every row started before it passes; a long window's gap counts as the basin's only within GAP_NEAR; a
    neighbour's gap beyond GAP_NEAR is not one of the basin."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    md = open(os.path.join(NR.OUT, 'noise_resonance.md'), encoding='utf-8').read()
    assert NR.report(**{k: raw.get(k) for k in NR.KEYS}) == md
    P = 40
    r = dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=[P, -1, -1], maxdev=[0.9, 0, 0], wins=[])
    fl = {5: 2, 1: 9}                                                    # p = 1 leaves at 10, after the alarm's period 7: no repeat
    d = dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), wins=[dict(p=p, d=0.0, left=fl.get(p, -1), maxdev=0.1) for p in range(P)])
    C8R = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=P,
                rows=[dict(h=37, q=3, left=P, maxdev=0.9), dict(h=39, q=1, left=6, maxdev=0.9)])]
    z = next(z for z in NR.c8d_alarms([r], [d], C8R)['alarms'] if z['p'] == 5)
    assert (z['nb'], z['rep'], z['c8r']) == (6, 0, False)
    run = lambda pairs: [dict(d=d_, left=(0 if x > NR.LEFT_C else -1), maxdev=x) for d_, x in pairs]
    C7 = [dict(s=0.02, opts=dict(seed=7, shift=25), runs=run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.6))))]
    C7L = [dict(s=0.02, opts=dict(seed=7, shift=25), np=20, runs=run(((0.1, 0.5), (0.2, 0.25), (0.3, 0.6))))]  # a gap, not the basin's
    assert "a gap of the basin survives over 20 in 0 and closes in m = 12, k_e = 10, s = 0.02, draw 7 from 25." in '\n'.join(NR.c7l_lines(C7L, C7))
    C7S = [dict(s=0.02, opts=dict(seed=7, shift=24), runs=run(((0.1, 0.5), (0.2, 0.15), (0.3, 0.6))), base=25, dshift=-1),
           dict(s=0.02, opts=dict(seed=7, shift=26), runs=run(((0.1, 0.5), (0.2, 0.05), (0.3, 0.6))), base=25, dshift=1)]
    assert "survives at both neighbours in 0 of the 1 windows, at one in 1, at neither in 0." in '\n'.join(NR.c7s_lines(C7S, C7))


def test_round94_main_c7s_superset_repeat_and_overlap(monkeypatch, tmp_path):
    """Referee round 94 (D3co-1, D3co-2, D3co-3, D3co-5: Q30): `main --only C7S` runs `part_C7S`; every C7L window's set contains its
    10-period window's (by construction: a test, not a finding); a repeat of an isolated alarm is another dense window leaving within one
    period of its leave period, not merely by then; C7S's line prints the least raise's change and the neighbours' overlap."""
    import exp_cohesion_leak as E
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    key = lambda r: (r.get('m', NR.C7_CFG[0]), r.get('ke', NR.C7_CFG[1]), r['s'], r['opts'].get('seed'), r['opts'].get('shift', 0))
    short = {key(r): r for r in list(raw['C7']) + list(raw['C9']) if 'freqs' not in r['opts']}
    for r in raw['C7L']:
        a = {x['d']: x['left'] for x in short[key(r)]['runs'] if x['left'] >= 0}
        b = {x['d']: x['left'] for x in r['runs']}
        assert all(b[d] == p for d, p in a.items())                     # every short leaver leaves in the same period over 20
    P = 40
    r = dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), left=[P, -1, -1], maxdev=[0.9, 0, 0], wins=[])
    fl = {5: 2, 1: 5, 2: 9}                                             # alarm at 7; window 1 leaves at 6 (a repeat), window 2 at 11 (not)
    d = dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=1), wins=[dict(p=p, d=0.0, left=fl.get(p, -1), maxdev=0.1) for p in range(P)])
    z = next(z for z in NR.c8d_alarms([r], [d])['alarms'] if z['p'] == 5)
    assert z['rep'] == 1
    s = '\n'.join(NR.c7s_lines(raw['C7S'], raw['C7'], raw['C9']))
    assert "the least ejecting raise changes between a window and its neighbour by up to ×4.08" in s and "median 0.975" in s
    got = []
    monkeypatch.setattr(NR.E, 'load', lambda p: {})
    monkeypatch.setattr(NR.E, 'save', lambda p, **kw: None)
    monkeypatch.setattr(NR, 'OUT', str(tmp_path))
    monkeypatch.setattr(NR, 'report', lambda **kw: '')
    monkeypatch.setattr(NR, 'part_C7S', lambda: got.append('C7S') or [])
    monkeypatch.setattr(NR, 'part_C7L', lambda: got.append('C7L') or [])
    monkeypatch.setattr(sys, 'argv', ['exp_noise_resonance.py', '--only', 'C7S'])
    NR.main()
    assert got == ['C7S']


def test_round98_the_anchored_overlap_and_planted_edges():
    """Referee round 98 (D3cs-1, D3cs-3: R01, R04, R06, R14): `c7s_overlap` — a base window whose least raise is 0 kept apart (lr0), not
    dropped; a far pair of two empty leaving sets has overlap 1; far pairs only within one draw; the anchored far pairs each base window
    against C7's windows 25 or more periods away. On the raw the anchored median 0.706 over 84 pairs (round 94's pooled 0.831 kept beside
    it). C8R: a replay lagging 2 periods counts as later and is named."""
    key = lambda r, sh: (12.0, 10.0, 'paradigm', r['s'], r['opts'].get('seed'), sh)
    mk = lambda s, seed, sh, ds: dict(s=s, opts=dict(seed=seed, shift=sh), runs=[dict(d=d, left=(1 if d in ds else -1)) for d in (0.0, 0.1, 0.2, 0.3)])
    b = mk(0.02, 7, 50, {0.0, 0.1})
    q1, q2 = mk(0.02, 7, 49, {0.1, 0.2}), mk(0.02, 7, 51, {0.1})
    f1, f2, od = mk(0.02, 7, 100, set()), mk(0.02, 7, 150, set()), mk(0.02, 3, 100, {0.2})
    ov = NR.c7s_overlap({key(b, 50): {-1: q1, 1: q2}}, {key(b, 50): b}, [b, f1, f2, od], key)
    assert ov['nbj'] == [1 / 3, 1 / 2] and ov['lr'] == [] and ov['lr0'] == 2
    assert sorted(ov['far']) == [0.0, 0.0, 1.0] and ov['anc'] == [0.0, 0.0] and ov['ancr'] == []
    b2 = mk(0.02, 7, 50, {0.1, 0.2})
    ov = NR.c7s_overlap({key(b2, 50): {1: q2}}, {key(b2, 50): b2}, [b2, mk(0.02, 7, 80, {0.3})], key)
    assert ov['lr'] == [(1.0, key(b2, 50), 51)] and ov['lr0'] == 0 and np.allclose(ov['ancr'], [3.0]) and ov['anc'] == [0.0]
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    md = open(os.path.join(NR.OUT, 'noise_resonance.md'), encoding='utf-8').read()
    assert ("against 0.706 for each of these base windows against C7's windows of its configuration, level and draw started 25 or more periods "
            "away (84 pairs; their least raises' ratio median ×1.36 over 84; referee round 98), and 0.831 over all of C7's levels pooled (420 pairs") in md
    C8 = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=3), left=[40, -1, -1], maxdev=[0.9, 0, 0], wins=[dict(p=10, d=0.0, left=-1, maxdev=0.1)])]
    C8R = [dict(m=12.0, ke=10.0, metric='paradigm', s=0.02, opts=dict(seed=3), left=40, rows=[dict(h=20, q=20, left=42, maxdev=0.5)])]
    s = '\n'.join(NR.c8_more_lines(C8, None, C8R))
    assert ("20 periods before: of 1 runs 0 leave within one period of the run (out of sample 0 of 1, in sample 0 of 0), 0 earlier, 1 later, "
            "0 never (m = 12, k_e = 10, s = 0.02, draw 3)") in s


def test_round98_c8a_the_fuller_restarts_and_draw_5(monkeypatch):
    """Referee round 98 (D3cs-2, D3cs-6): `offsets_start`'s `only` applies the chosen kernels' offsets alone; `c8a_cases` takes C3's period
    means where the run is C3's and reruns the others, with draw 5's probe; `c8a_lines` classes each restart against the reduced window
    (kept, delayed, earlier, removed); on the raw 'w' keeps all 12 isolated alarms, 'wc' removes 6, delays 1 (draw 2 from 58: 2 → 5) and
    keeps 5, the reruns leave at C8's periods, and draw 5's kernel 1 leaves only with kernel 3's offset."""
    L = NR.law(NR.C7_CFG[0])
    al0, c0, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    fit = [abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][i])) for i in idx]
    W = [[fit[0] + 0.05, fit[1] + 0.03, fit[2] + 0.02]] * 3
    C = [[0.1, 0.2, 0.3]] * 3
    al_all, _ = NR.offsets_start(W, C, 2, False)(L, al0, c0, idx)
    al3, c3 = NR.offsets_start(W, C, 2, False, only=(2,))(L, al0, c0, idx)
    ref = NR.raised(L, np.array(al0, float), idx, 2, 0.02)
    assert np.allclose(al3, ref, rtol=0, atol=1e-15) and not np.allclose(al3, al_all) and np.array_equal(c3, c0)
    _, cc = NR.offsets_start(W, C, 2, True, only=(0,))(L, al0, c0, idx)
    assert cc[idx[0]] == 0.1 and cc[idx[1]] == c0[idx[1]]
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    cases = NR.c8a_cases(raw['C8'], raw['C8D'], raw['C8R'], raw['C3'])
    assert len(cases) == 7 and sum(len(c[5]) for c in cases) == 12
    assert [c[6] is not None for c in cases] == [c[1] == 10.0 for c in cases] and sum(c[7] is not None for c in cases) == 1
    assert next(c for c in cases if c[7] is not None)[3:5] == (0.03, {'seed': 5}) and NR.C8A_PROBE == (0.03, 5, 150, ((), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)))
    A = raw['C8A']
    rows = [x for r in A for x in r['rows']]
    assert all(x['w'][0] == x['reduced'] for x in rows) and sum(x['wc'][0] < 0 for x in rows) == 6
    assert [(x['p'], x['reduced'], x['wc'][0]) for x in rows if x['wc'][0] > x['reduced']] == [(58, 2, 5)]
    base = {(r['m'], r['ke'], r['s'], r['opts'].get('seed', 7)): r['left'][0] for r in raw['C8']}
    assert all(base[(r['m'], r['ke'], r['s'], r['opts'].get('seed', 7))] == r['left'] for r in A if r['src'] == 'rerun')
    pb = next(r['probe'] for r in A if r.get('probe'))
    assert [(x['only'], x['res'][0]) for x in pb['rows']] == [([], -1), ([0], -1), ([1], -1), ([2], 10), ([0, 1], -1), ([0, 2], 10), ([1, 2], 10), ([0, 1, 2], 10)]
    syn = [dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1), src='C3', left=None, probe=None,
                rows=[dict(p=5, reduced=2, w=[2, 0.5], wc=[1, 0.5]), dict(p=9, reduced=2, w=[-1, 0.1], wc=[3, 0.4])])]
    s = '\n'.join(NR.c8a_lines(syn))
    assert ("('w'): kept 1, removed 1; and the centres at their period means too ('wc'): delayed 1 (m = 12, k_e = 20, s = 0.06, draw 1 from 9: "
            "period 2 → 3), earlier 1 (m = 12, k_e = 20, s = 0.06, draw 1 from 5: period 2 → 1)") in s
    assert NR.c8a_lines(None) == [] and "(0 of the 0 reruns leaving at C8's period)" in s


def test_round98_c8a_job_wiring(monkeypatch):
    """`c8a_job` (rounds 98, 102): C's loop always rerun with the run's (m, k_e, s, draw, metric); C3's period means used as given where
    passed; per window 'w' without and 'wc' with the centres (`offsets_start`) and the one-ingredient variants (`state_start`, the run's
    instantaneous state at the window's start pT, output 8p), under the run's own (m, k_e, metric, s, draw) and noise from the window's
    start; the probe's subsets under the noise and without it (`edge_run`, no noise); round 124: `c8a_pT` at every window from p ≥ 2 with
    the run's (m, k_e, metric, s, draw), the W used (C3's or the rerun's), p and the reduced leave, its fields merged into the row."""
    calls, quiet, runs = [], [], []
    monkeypatch.setattr(NR, 'win_run', lambda m, ke, mt, s, opts, d, p, NP=NR.NP_C7, start=None: (calls.append((m, ke, mt, s, dict(opts), p, d)) or (3, 0.4)))
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=None, start=None, **kw: (quiet.append((m, ke, mt, k, NP, kw.get('noise'))) or (-1, 1e-5)))
    Wp, Cp = np.full((12, 3), 0.5), np.tile(NR.CT.CSTAR, (12, 1))
    Cp[6:, 0] += 1.0
    rr = dict(w=np.arange(97 * 21, dtype=float).reshape(97, 21), c=-np.arange(97 * 21, dtype=float).reshape(97, 21))
    monkeypatch.setattr(NR, 'c_run', lambda m, ke, s, opts: (runs.append((m, ke, s, dict(opts))) or (rr, None, None, None, [2, 5, 9], Wp, Cp)))
    seen = []
    real = NR.state_start
    monkeypatch.setattr(NR, 'state_start', lambda W, C, Wi, Ci, p, spec: (seen.append((p, Wi, Ci, spec)) or real(W, C, Wi, Ci, p, spec)))
    pt = []                                                              # round 124 (D3ds-1): the states at pT, per window
    monkeypatch.setattr(NR, 'c8a_pT', lambda *a: (pt.append(a) or dict(dp=0.25, pT=[a[6]])))
    W, C = [[0.5, 0.4, 0.3]] * 12, [[-0.6, 0.0, 0.6]] * 12
    out = NR.c8a_job((12.0, 10.0, 'wfr', 0.02, {'seed': 7}, [(5, 2), (9, 1)], (W, C), (8, ((2,), (0, 1)))))
    assert pt == [(12.0, 10.0, 'wfr', 0.02, {'seed': 7}, W, 5, 2), (12.0, 10.0, 'wfr', 0.02, {'seed': 7}, W, 9, 1)]
    assert [(x['dp'], x['pT']) for x in out['rows']] == [(0.25, [5]), (0.25, [9])]
    assert runs == [(12.0, 10.0, 0.02, {'seed': 7, 'metric': 'wfr'})] and out['src'] == 'C3' and out['left'] == 7
    assert [x['p'] for x in out['rows']] == [5, 9] and [x['reduced'] for x in out['rows']] == [2, 1]
    assert all(x[t] == [3, 0.4] for x in out['rows'] for t in ('w', 'wc') + tuple(v for v, _ in NR.C8A_VARIANTS))
    assert all(c[:5] == (12.0, 10.0, 'wfr', 0.02, {'seed': 7}) and c[6] == 0.0 for c in calls)
    n_ = 2 + len(NR.C8A_VARIANTS) + 3 + len(NR.C8A_SCAN)                  # round 105: the centre-alone starts, the scan, the warm start
    assert [c[5] for c in calls] == [5] * n_ + [4, 4, 3, 3, 2, 2] + [9] * n_ + [8, 8, 7, 7, 6, 6] + [8, 8]   # round 111: the earlier starts
    assert seen[0][0] == 5 and seen[0][1] == [float(v) for v in rr['w'][40, [2, 5, 9]]] and seen[0][2] == [float(v) for v in rr['c'][40, [2, 5, 9]]]
    assert [x[3] for x in seen[:len(NR.C8A_VARIANTS)]] == [sp for _, sp in NR.C8A_VARIANTS]
    assert out['rows'][0]['c1'] == -0.6 - NR.CT.CSTAR[0]
    assert out['probe']['rows'] == [dict(only=[2], res=[3, 0.4]), dict(only=[0, 1], res=[3, 0.4])]
    assert out['probe']['quiet'] == [dict(only=[2], res=[-1, 1e-5]), dict(only=[0], res=[-1, 1e-5])]
    assert quiet == [(12.0, 10.0, 'wfr', 0, NR.NP_C7, None)] * 2
    out = NR.c8a_job((12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(4, 2)], None, None))
    assert out['src'] == 'rerun' and out['left'] == 7 and out['probe'] is None
    assert pt[-1] == (12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, Wp.tolist(), 4, 2)          # the rerun's W
    pt.clear()
    NR.c8a_job((12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(1, 2)], None, None))           # p = 1: no earlier start, nothing at pT
    assert pt == []


def test_round102_state_start_the_rerun_and_computed_conclusions(monkeypatch):
    """Referee round 102 (D3cw-1 … D3cw-5): `state_start` with every kernel ('mean', 'fit') is `offsets_start` without the centres and
    ('mean', 'mean') with them, 'inst' the run's own values; `c8a_job` rerun on draw 7's job equals the raw's; `c8a_lines` prints a
    conclusion only where the counts imply it (contrary synthetic inputs); the anchored pairs count a zero least raise apart."""
    L = NR.law(NR.C7_CFG[0])
    al0, c0, idx = NR.DZ.fitted_state(L, NR.EPS, NR.CT.CSTAR, NR.CT.ASTAR)
    W, C = [[0.6, 0.45, 0.35]] * 4, [[-0.62, 0.01, 0.6]] * 4
    for cen, spec in ((False, (('mean', 'fit'),) * 3), (True, (('mean', 'mean'),) * 3)):
        a1, c1 = NR.offsets_start(W, C, 2, cen)(L, al0, c0, idx)
        a2, c2 = NR.state_start(W, C, [0, 0, 0], [0, 0, 0], 2, spec)(L, al0, c0, idx)
        assert np.array_equal(a1, a2) and np.array_equal(c1, c2)
    a3, c3 = NR.state_start(W, C, [0.7, 0.4, 0.3], [-0.5, 0.0, 0.5], 2, (('inst', 'inst'), ('fit', 'fit'), ('fit', 'mean')))(L, al0, c0, idx)
    assert c3[idx[0]] == -0.5 and c3[idx[1]] == c0[idx[1]] and c3[idx[2]] == 0.6
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    cases = NR.c8a_cases(raw['C8'], raw['C8D'], raw['C8R'], raw['C3'])
    job = next(c for c in cases if c[1] == 10.0 and c[3] == 0.02)
    ref = next(r for r in raw['C8A'] if r['ke'] == 10.0 and r['s'] == 0.02)
    out = NR.c8a_job(job)
    assert {k: v for k, v in out.items() if k != 'secs'} == {k: v for k, v in ref.items() if k != 'secs'}
    mk = lambda red, **kw: dict(p=10, reduced=red, c1=-0.05, **{t: [kw.get(t, red), 0.4] for t in ('w', 'wc', 'w+c1', 'w+c23', 'i', 'wc+i1c', 'wc+i1w', 'wc+i23')})
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1), src='rerun', left=None, probe=None)
    good = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2, wc=-1, **{'w+c1': -1, 'i': -1}), mk(2, i=-1, **{'wc+i1c': -1})])]))
    assert "no alarm comes from dropping kernels 2 and 3's coefficient offsets" in good and "the centres' effect is kernel 1's centre alone" in good
    assert "the 1 that 'wc' leaves go with the formed kernels' instantaneous state" in good
    bad = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2, w=-1, wc=-1), mk(2, **{'w+c1': -1})])]))
    assert "no alarm comes from" not in bad and "kernel 1's centre alone" not in bad and "go with the formed kernels' instantaneous" not in bad
    pb = dict(p=150, rows=[dict(only=[2], res=[10, 0.5]), dict(only=[0], res=[4, 0.5])], quiet=[dict(only=[2], res=[-1, 1e-5])])
    s = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2)], probe=pb)]))
    assert "needs kernel 3's offset" not in s and "kernel 1's leave with kernel 3's offset needs the noise" in s
    half = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2, wc=-1, **{'w+c1': -1, 'w+c23': -1})])]))   # 'w+c1' as 'wc', 'w+c23' not as 'w'
    assert "kernel 1's centre alone" not in half
    pq = dict(pb, quiet=[dict(only=[2], res=[12, 0.5])])                                        # kernel 3's offset leaves without the noise
    assert "kernel 1's leave with kernel 3's offset needs the noise" not in '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2)], probe=pq)]))
    pn = dict(pb, rows=[dict(only=[2], res=[-1, 0.5]), dict(only=[0], res=[4, 0.5])])   # round 111's R06: kernel 3's offset not leaving
    assert "needs the noise" not in '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2)], probe=pn)]))  # with the noise either: no clause
    key = lambda q, sh: (12.0, 10.0, 'paradigm', q['s'], q['opts'].get('seed'), sh)
    mkw = lambda sh, ds: dict(s=0.02, opts=dict(seed=7, shift=sh), runs=[dict(d=d, left=(1 if d in ds else -1)) for d in (0.0, 0.1, 0.2)])
    b = mkw(50, {0.1, 0.2})
    ov = NR.c7s_overlap({key(b, 50): {1: mkw(51, {0.1})}}, {key(b, 50): b}, [b, mkw(100, {0.0, 0.2})], key)
    assert ov['ancr'] == [] and ov['ancr0'] == 1


def test_round105_the_rerun_path_the_conditions_and_the_centre_coordinate(monkeypatch):
    """Referee round 105 (D3cz-1 … D3cz-6): `c8a_job`'s rerun path passes the rerun's (W, C) in that order to every start; the centre
    clause needs the leave periods equal and some window where 'wc' ≠ 'w'; "needs kernel 3's offset" only when every subset without kernel
    3 was run; the C8R clause only when its aggregate holds; the round-105 sentence's conclusions from the counts (rounds 114, 119: the
    cold restart's only if no earlier start with its own δ leaves near or before the alarm's period and (round 124) every own-δ start's
    state at pT with δ(p) passes, a late δ(p) start's drift only where its arrival is above δ(p) and δ(p) at pT removes it; one coordinate
    only if every scan is monotone with both outcomes, offset 0 the reduced window, and the three samples
    merged into it with the removing set a prefix)."""
    seen = []
    Wp, Cp = 0.3 + np.arange(900).reshape(300, 3) * 1e-4, -0.1 + np.arange(3) * 0.5 + np.arange(900).reshape(300, 3) * 1e-5   # distinct rows
    fake_r = dict(w=np.full((2400, 3), 0.3), c=np.full((2400, 3), -0.1))
    monkeypatch.setattr(NR, 'c_run', lambda m, ke, s, opts: (fake_r, None, None, None, [0, 1, 2], Wp, Cp))
    monkeypatch.setattr(NR, 'left_periods', lambda *a, **k: [5])
    monkeypatch.setattr(NR, 'offsets_start', lambda W, C, p, cen, sub=None: seen.append(('o', W, C)) or None)
    monkeypatch.setattr(NR, 'state_start', lambda W, C, Wi, Ci, p, spec: seen.append(('s', W, C)) or None)
    monkeypatch.setattr(NR, 'win_run', lambda *a, **k: (-1, 0.0))
    seenT = []
    monkeypatch.setattr(NR, 'c8a_pT', lambda m, ke, mt, s, opts, W, p, R: seenT.append(W) or {})   # round 124: at pT, the rerun's W
    NR.c8a_job((12.0, 20.0, 'paradigm', 0.06, dict(seed=1), [(50, 2)], None, None))
    assert seenT == [Wp.tolist()]
    j4 = NR.c8a_job((12.0, 20.0, 'paradigm', 0.06, dict(seed=1), [(4, 2)], None, None))   # round 114 (J04): p − k − 1 = 0 is valid
    assert [y[0] for y in j4['rows'][0]['warm_pk']] == [2, 3] and 'c1i' in j4['rows'][0] and 'c1p' in j4['rows'][0]
    assert seen and all(W[k] == Wp.tolist()[k] and C[k] == Cp.tolist()[k] for _, W, C in seen for k in range(50, 300))   # round 111: every row
    monkeypatch.undo()
    mk = lambda red, **kw: dict(p=10, reduced=red, c1=-0.05, **{t: [kw.get(t, red), 0.4] for t in ('w', 'wc', 'w+c1', 'w+c23', 'i', 'wc+i1c', 'wc+i1w', 'wc+i23')})
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1), src='rerun', left=None, probe=None)
    same = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2), mk(3)])]))              # 'wc' = 'w' everywhere: no centre effect to attribute
    assert "kernel 1's centre alone" not in same
    dl = '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2, wc=5, **{'w+c1': 9})])]))   # both delay, to different periods
    assert "kernel 1's centre alone" not in dl
    pb = dict(p=150, rows=[dict(only=[2], res=[10, 0.5]), dict(only=[0], res=[-1, 0.5])], quiet=[])
    assert "needs kernel 3's offset" not in '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2)], probe=pb)]))   # {} and {2} not run
    pb8 = dict(p=150, rows=[dict(only=list(o), res=[10 if 2 in o else -1, 0.5]) for o in ((), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2))], quiet=[])
    assert "needs kernel 3's offset, not its own" in '\n'.join(NR.c8a_lines([dict(r, rows=[mk(2)], probe=pb8)]))
    row = lambda warm, scan: dict(mk(2), scan=scan, warm=[warm, 0.1], warm_m=[[1, 'own', -1, 0.1]],   # round 124: 'warm_m' passing
                                  **{'r+c1m': [-1, 0.1], 'r+c1i': [-1, 0.1], 'r+c1p': [-1, 0.1]})
    one = ("at the 2 windows with an own-δ start's state at pT (1 where no own-δ start leaves near or before the alarm's period, 1 where one "   # round 128
           "reproduces the alarm): passing at every own-δ k at 2 — there the alarm is the cold restart's")   # (D3dw-2): the reproduction classed too
    two = "at the 2 windows where no own-δ start leaves near or before the alarm's period: passing at every own-δ k at 2 — there the alarm is the cold restart's"
    good = [[-0.06, -1, 0.1], [-0.02, -1, 0.1], [0.0, 2, 0.4], [0.04, 2, 0.4]]
    bad = [[-0.06, 2, 0.4], [-0.02, -1, 0.1], [0.0, 2, 0.4], [0.04, 2, 0.4]]          # removed inward, not outward: not monotone
    s1 = '\n'.join(NR.c8a_lines([dict(r, rows=[row(-1, good), row(-1, good)])]))
    assert two in s1 and "at the 2 kernel 1's centre at the restart is one coordinate against a per-window threshold" in s1
    s2 = '\n'.join(NR.c8a_lines([dict(r, rows=[row(3, good), row(-1, bad)])]))     # a warm start keeping it (2 + 1), a non-monotone scan
    assert one in s2 and "at 1 of 2 windows (not at" in s2                          # round 124 (D3ds-3): the clause from nco alone
    assert "at 1 of 2 windows (not at m = 12, k_e = 20, s = 0.06, draw 1 from 10: not monotone, its largest distance 0.4000) — at the 1 kernel 1's centre at the restart is one coordinate" in s2
    flat = [[-0.06, 2, 0.31], [-0.02, 2, 0.31], [0.0, 2, 0.31], [0.04, 2, 0.31]]       # round 111 (D3df-1): no offset removes it
    allr = [[-0.06, -1, 0.1], [-0.02, -1, 0.1], [0.0, -1, 0.1], [0.04, -1, 0.1]]       # removed everywhere, offset 0 not the alarm
    s3 = '\n'.join(NR.c8a_lines([dict(r, rows=[row(-1, good), row(-1, flat), row(-1, allr)])]))
    assert "at 1 of 3 windows (not at" in s3 and "no offset removes it" in s3 and "at the 1 kernel 1's centre at the restart is one coordinate against" in s3
    assert "every offset removes it" in s3
    off0 = [[-0.06, -1, 0.1], [-0.02, -1, 0.1], [0.0, 3, 0.4], [0.04, 3, 0.4]]        # monotone, both outcomes, offset 0 not the reduced (2)
    s5 = '\n'.join(NR.c8a_lines([dict(r, rows=[row(-1, good), row(-1, off0)])]))   # round 111's mutant T05
    assert "at 1 of 2 windows (not at" in s5 and "offset 0 not the reduced window" in s5 and "at the 1 kernel 1's centre at the restart is one coordinate against" in s5
    s6 = '\n'.join(NR.c8a_lines([dict(r, rows=[row(-1, good), row(-1, bad)])]))
    assert "not monotone" in s6
    near = dict(row(-1, good), warm=[4, 0.4])                                        # the warm run leaving one period after the alarm's
    s4 = '\n'.join(NR.c8a_lines([dict(r, rows=[near, row(-1, good)])]))
    assert one in s4
    deep = dict(row(-1, good), warm_k=[[2, 4, 0.4]])                                 # two periods earlier with its own δ, at the alarm's period
    assert one in '\n'.join(NR.c8a_lines([dict(r, rows=[deep, row(-1, good)])]))
    dp = dict(row(-1, good), warm_pk=[[2, 4, 0.4]])                                  # round 119 (D3dn-1): with δ(p) only — the control's drift
    s_dp = '\n'.join(NR.c8a_lines([dict(r, rows=[dp, row(-1, good)])]))
    assert "with the window's δ(p) instead at 1 more (m = 12, k_e = 20, s = 0.06, draw 1 from 10: 2 periods earlier)" in s_dp and two in s_dp
    # round 114 (D3di-4, -5): both sides of every boundary (the alarm at the window's period 2 = R; k periods earlier it is R + k)
    run = lambda *xs: '\n'.join(NR.c8a_lines([dict(r, rows=list(xs) + [row(-1, good)])]))
    assert "with the window's δ(p) instead at 1 more" in run(dict(row(-1, good), warm_p=[3, 0.4]))           # L01: 'warm_p' at R + 1
    assert one in run(dict(row(-1, good), warm_k=[[3, 5, 0.4]]))                                              # L02: k = 3 at R + 3
    assert one in run(row(2, good))                                                                           # L06: warm at R (R + 1 − 1)
    assert one in run(dict(row(-1, good), warm_k=[[3, 6, 0.4]]))                                              # L09: k = 3 at R + 4
    s7 = run(row(5, good))                                                                                    # L03: warm at R + 3 (late)
    assert two in s7 and "after the alarm's period an earlier start leaves at" in s7
    assert "1 period earlier with its own δ(p − 1): 4 periods after the window's start, the alarm's at 2" in s7 and "own drift" not in s7
    s7n = run(dict(row(-1, good), warm_pk=[[2, 6, 0.3001]]))                                                  # round 119: a δ(p) late leave,
    assert "2 periods earlier with δ(p): 4 periods after" in s7n and "above the cut by 1.0e-04" in s7n and "own drift" not in s7n   # no pT
    pz = lambda k, kd, w1: dict(k=k, kd=kd, w1=w1, c=[-0.05, 0.0, 0.0], b=[4, 0.3001])                         # round 124 (D3ds-1): at pT
    lt = lambda w1, lm, **kw: run(dict(row(-1, good), dp=0.07, pT=[pz(2, 'p', w1)], warm_m=[[1, 'own', -1, 0.1], [2, 'p', lm, 0.1]], **kw))
    s7p = lt(0.071, -1, warm_pk=[[2, 6, 0.3001]])
    assert "2 periods earlier with δ(p): 4 periods after the window's start, the alarm's at 2, its largest distance 0.3001, above the cut by " \
           "1.0e-04, arriving at pT at 0.07100 above δ(p) = 0.07000" in s7p and "the start's own drift, not the alarm" in s7p
    assert "own drift" not in lt(0.069, -1, warm_pk=[[2, 6, 0.3001]]) and "own drift" not in lt(0.071, 5, warm_pk=[[2, 6, 0.3001]])
    s7w = run(dict(row(-1, good), warm_p=[5, 0.3001], dp=0.07, pT=[pz(1, 'p', 0.071)], warm_m=[[1, 'own', -1, 0.1], [1, 'p', -1, 0.1]]))
    assert "1 period earlier with δ(p): 4 periods after" in s7w and "the start's own drift, not the alarm" in s7w   # round 119's Y09
    s7k = run(dict(row(-1, good), warm_k=[[2, 6, 0.3001]]))                                                   # Y10: an own-δ start late
    assert "2 periods earlier with its own δ(p − 2): 4 periods after" in s7k and "own drift" not in s7k
    s8 = run(row(1, good))                                                                                    # L04, P1: warm at R − 1 (early)
    assert "leaves before the alarm's period at 1 (not tested:" in s8 and "the alarms are the cold restart's" not in s8
    assert one.replace("1 where one reproduces the alarm", "1 where one leaves before it") in s8              # (the other window: K01)
    s9 = run(dict(row(-1, good), warm_k=[[3, 3, 0.4]]))                                                       # L05: k = 3 at R + 1 (early)
    assert "leaves before the alarm's period at 1 (not tested:" in s9
    s9p = run(dict(row(-1, good), warm_pk=[[2, 2, 0.4]]))                                                     # round 119's Y03: a δ(p) start
    assert two in s9p and "not tested" not in s9p                             # leaving early tests nothing
    s10 = run(dict(row(-1, good), c1=-0.08, **{'r+c1m': [2, 0.4]}))                                           # P7: a sample against the scan
    assert "at 1 of 2 windows (not at" in s10 and "a sample against the scan" in s10 and "at the 1 kernel 1's centre at the restart is one coordinate" in s10
    for kw in (dict(c1=-0.02, **{'r+c1m': [2, 0.4]}), dict(c1=0.0), dict(c1i=-0.05, **{'r+c1i': [2, 0.4]}),   # round 119 (N01–N05):
               dict(c1p=-0.05, **{'r+c1p': [2, 0.4]}), dict(c1i=-0.015, c1p=-0.005, **{'r+c1i': [2, 0.4]})):  # ends, the samples, inside
        assert "a sample against the scan" in run(dict(row(-1, good), **kw)), kw
    s12 = run(dict(row(-1, good), c1=0.02))                                                                    # K05: the keep side
    assert "a sample against the scan" in s12 and "at the 1 kernel 1's centre at the restart is one coordinate" in s12
    s11 = run(dict(row(-1, good), c1i=-0.01, **{'r+c1i': [-1, 0.1]}))                                          # inside the bracket: no test
    assert "at 2 of 2 windows" in s11
    s11e = run(dict(row(-1, good), c1i=-0.02, **{'r+c1i': [-1, 0.1]}))                                         # round 119's Y07: at the
    assert "(merged into the scan; 0 of the 3 in a bracket inside it" in s11e and "1 of the 3 in a bracket inside" in s11                # bracket's end: outside
    fl = [[-0.06, 2, 0.35], [-0.02, 2, 0.35], [0.0, 2, 0.3104], [0.04, 2, 0.35]]                             # L08: offset 0's distance
    assert "no offset removes it, its largest distance 0.3104" in run(row(-1, fl))


def test_round124_c8a_pT_wiring(monkeypatch):
    """Referee round 124 (D3ds-1, -4, -5): `c8a_pT` runs each earlier start to pT under the noise from p − k with its own δ(p − k) or with
    δ(p), keeps its state there, continues it ('b') and runs it with kernel 1's α set to the cold restart's, the rest of the state kept
    ('warm_m'); runs the cold restart with the start's α ('cold_m') and the own-δ states with one formed centre at its feature ('warm_c');
    runs to p + R + 2 only where that passes the window's end ('ext'); and bisects kernel 1's threshold only between a passing arrival and
    a leaving δ(p), on the cold restart's state and on each own-δ start's. Round 128 (D3dw-3): then 'cold_c', the cold restart (its α, its
    centres, the plant at x_m(0)) with one own-δ start's formed centre at its value at pT, per k and kernel."""
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    off = lambda al: abs(float(NR.mobilities(L, NR.EPS, al, 'paradigm')[0][idx[0]])) - w0
    calls, rets = [], []

    def fake(m, ke, mt, s, opts, shift, NP, al, c, x0=None):
        calls.append((m, ke, mt, s, dict(opts), shift, NP, np.array(al, float), np.array(c, float), x0))
        a = off(al)
        if NP < NR.NP_C7:                                                # a start's own segment to pT: it arrives 0.0015 lower
            alt = NR.raised(L, al0, idx, 0, a - 0.0015)
            alt[idx[1]] += 1e-3                                         # kernel 2 moved too: kept from pT, not reset
            rets.append((0.001 * shift, alt, np.array(c, float) + 0.01))
            return [0.1] * NP, rets[-1]
        return ([0.5] if a > (0.0763 if x0 is None else 0.0762) else [0.1]) * NP, None   # leaves above a threshold: cold, warm
    monkeypatch.setattr(NR, 'seg_run', fake)
    W = [[w0 + 0.07 + 0.001 * (q + 1), 0.0, 0.0] for q in range(12)]       # δ(q) = |W[q − 1][0]| − w₀ = 0.07 + 0.001q
    out = NR.c8a_pT(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, W, 8, 10)
    dp = abs(W[7][0]) - w0
    alc = NR.raised(L, al0, idx, 0, dp)
    assert out['dp'] == dp and all(c_[:5] == (12.0, 20.0, 'paradigm', 0.06, {'seed': 1}) for c_ in calls) and len(calls) == 1 + 24 + 9 + 48 + 9
    assert calls[0][5:7] == (8, 12) and np.array_equal(calls[0][7], alc) and np.array_equal(calls[0][8], c0) and calls[0][9] is None
    st = [c_ for c_ in calls if c_[6] < NR.NP_C7]
    ds = [abs(W[8 - k - 1][0]) - w0 if kd == 'own' else dp for k in (1, 2, 3) for kd in ('own', 'p')]
    assert [(c_[5], c_[6]) for c_ in st] == [(7, 1), (7, 1), (6, 2), (6, 2), (5, 3), (5, 3)]
    assert all(np.array_equal(c_[7], NR.raised(L, al0, idx, 0, d)) and np.array_equal(c_[8], c0) and c_[9] is None for c_, d in zip(st, ds))
    x1, a1, c1 = rets[0]
    b, m_, cm = calls[2], calls[3], calls[4]
    assert b[5:7] == (8, 12) and np.array_equal(b[7], a1) and np.array_equal(b[8], c1) and b[9] == x1          # continued from pT
    am = a1.copy()
    am[idx[0]] = alc[idx[0]]
    assert m_[5:7] == (8, 12) and np.array_equal(m_[7], am) and np.array_equal(m_[8], c1) and m_[9] == x1       # the cold restart's α
    ac = alc.copy()
    ac[idx[0]] = a1[idx[0]]
    assert cm[5:7] == (8, 10) and np.array_equal(cm[7], ac) and np.array_equal(cm[8], c0) and cm[9] is None     # the cold, the arrival's α
    for j, c_ in enumerate(calls[5:8]):
        cj = c1.copy()
        cj[idx[j]] = c0[idx[j]]
        assert c_[5:7] == (8, 10) and np.array_equal(c_[7], a1) and np.array_equal(c_[8], cj) and c_[9] == x1
    assert (calls[8][6], calls[12][6]) == (1, 2)                                                                # no 'warm_c' after δ(p)
    ar = [off(r_[1]) for r_ in rets]
    assert [(z['k'], z['kd']) for z in out['pT']] == [(k, kd) for k in (1, 2, 3) for kd in ('own', 'p')] and [z['w1'] for z in out['pT']] == ar
    assert all(np.allclose(z['c'], 0.01, rtol=0, atol=1e-12) for z in out['pT']) and [z['b'] for z in out['pT']] == [[-1, 0.1], [1, 0.5]] * 3
    assert out['warm_m'] == [[k, kd, 1, 0.5] for k in (1, 2, 3) for kd in ('own', 'p')]
    assert out['cold_m'] == [[k, kd] + ([-1, 0.1] if kd == 'own' else [1, 0.5]) for k in (1, 2, 3) for kd in ('own', 'p')]
    assert out['warm_c'] == [[k, j, -1, 0.1] for k in (1, 2, 3) for j in (0, 1, 2)]
    assert out['ext'] == dict(n=12, cold=[0.5, 0.5], b=[[k, kd] + [0.1 if kd == 'own' else 0.5] * 2 for k in (1, 2, 3) for kd in ('own', 'p')],
                              m=[[k, kd, 0.5, 0.5] for k in (1, 2, 3) for kd in ('own', 'p')])
    lo, hi = out['thr']['cold']                                         # from the largest arrival at which the cold restart passes
    assert lo <= 0.0763 < hi and np.isclose(hi - lo, (dp - ar[0]) / 2 ** NR.NB_C8A, rtol=1e-6, atol=0)
    assert [k for k, _, _ in out['thr']['warm']] == [1, 2, 3]
    for (k, lo, hi), a in zip(out['thr']['warm'], ar[::2]):
        assert lo <= 0.0762 < hi and np.isclose(hi - lo, (dp - a) / 2 ** NR.NB_C8A, rtol=1e-6, atol=0)
    bc, bw = calls[34], calls[46]                                       # the first step of each: the cold, then kernel 1's own-δ k = 1
    assert bc[5:7] == (8, 10) and np.array_equal(bc[8], c0) and bc[9] is None and np.isclose(off(bc[7]), 0.5 * (ar[0] + dp), rtol=0, atol=1e-15)
    assert np.array_equal(np.delete(bc[7], idx[0]), np.delete(alc, idx[0]))
    assert np.array_equal(bw[8], c1) and bw[9] == x1 and np.isclose(off(bw[7]), 0.5 * (ar[0] + dp), rtol=0, atol=1e-15)
    assert np.array_equal(np.delete(bw[7], idx[0]), np.delete(a1, idx[0]))
    for q, c_ in enumerate(calls[82:]):                                 # round 128: 'cold_c', after the bisections
        k, j = q // 3 + 1, q % 3
        cw = c0.copy()
        cw[idx[j]] = rets[2 * (k - 1)][2][idx[j]]                       # the own-δ start k's centre j at pT, the others the fit's
        assert c_[5:7] == (8, 10) and np.array_equal(c_[7], alc) and np.array_equal(c_[8], cw) and c_[9] is None
    assert out['cold_c'] == [[k, j, 1, 0.5] for k in (1, 2, 3) for j in (0, 1, 2)]
    calls.clear()
    rets.clear()
    o2 = NR.c8a_pT(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, W, 8, 8)      # R + 2 = NP_C7: nothing past the window
    assert 'ext' not in o2 and calls[0][6] == NR.NP_C7 and {c_[6] for c_ in calls} == {1, 2, 3, NR.NP_C7}
    o3 = NR.c8a_pT(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, W, 3, 2)       # p − k − 1 ≥ 0: k = 1, 2 only
    assert [(z['k'], z['kd']) for z in o3['pT']] == [(1, 'own'), (1, 'p'), (2, 'own'), (2, 'p')]
    Wl = [[w0 + 0.06 + 0.001 * (q + 1), 0.0, 0.0] for q in range(12)]       # δ(p) = 0.068, below both thresholds: no bracket
    o4 = NR.c8a_pT(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, Wl, 8, 2)
    assert o4['thr'] == {} and all(y[2] < 0 for y in o4['warm_m'])
    Wh = [[w0 + 0.079 + 0.0001 * (q + 1), 0.0, 0.0] for q in range(12)]     # every arrival above both thresholds: no bracket either
    o5 = NR.c8a_pT(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, Wh, 8, 2)
    assert o5['thr'] == {} and all(y[2] >= 0 for y in o5['cold_m']) and all(z['b'][0] >= 0 for z in o5['pT'])


def test_round124_c8a_lines_at_pT_and_the_labels():
    """Referee round 124 (D3ds-1, -3, -4, -5, -6, -7): `c8a_lines` calls a window the cold restart's only where every own-δ start's state
    at pT with δ(p) ('warm_m') passes, the coefficient's where every one reproduces the alarm, neither otherwise, and not run without it;
    prints the arrivals, the cold restart with them ('cold_m') and the thresholds only where they say so; the centre transplants by kernel
    and the merged scan's side; the periods past the window's end; each late leave labelled by its own start (P1, P2) and the drift label
    only where every late start's arrival is above δ(p) and δ(p) at pT removes it (P3); δ(p)-only reproductions exclude the own-δ ones
    (P4); every sample in a bracket counted (P6) and none without one (P7)."""
    mk = lambda red, **kw: dict(p=10, reduced=red, c1=-0.05, **{t: [kw.get(t, red), 0.4] for t in ('w', 'wc', 'w+c1', 'w+c23', 'i', 'wc+i1c', 'wc+i1w', 'wc+i23')})
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1), src='rerun', left=None, probe=None)
    good = [[-0.06, -1, 0.1], [-0.02, -1, 0.1], [0.0, 2, 0.4], [0.04, 2, 0.4]]
    row = lambda **kw: {**mk(2), 'scan': good, 'warm': [-1, 0.1], 'warm_m': [[1, 'own', -1, 0.1]], 'r+c1m': [-1, 0.1], 'r+c1i': [-1, 0.1],
                        'r+c1p': [-1, 0.1], **kw}
    run = lambda *xs: '\n'.join(NR.c8a_lines([dict(r, rows=list(xs) + [row()])]))
    nm = "m = 12, k_e = 20, s = 0.06, draw 1 from 10"
    pz = lambda k, kd, w1, c=(-0.05, 0.0, 0.0): dict(k=k, kd=kd, w1=w1, c=list(c), b=[-1, 0.1])
    co = row(warm_p=[3, 0.31], dp=0.0767, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0776)], warm_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31]],
             cold_m=[[1, 'own', -1, 0.03], [1, 'p', 2, 0.31]], thr=dict(cold=[0.07635, 0.07636], warm=[[1, 0.07633, 0.07634]]))
    s = run(co)
    assert (f"with the window's δ(p) instead at 1 more ({nm}: 1 period earlier, arriving at pT with kernel 1's coefficient offset 0.07760 "
            "above δ(p) = 0.07670)") in s
    assert ("at the 2 windows where no own-δ start leaves near or before the alarm's period: passing at every own-δ k at 1 — there the alarm "
            f"is the cold restart's —, leaving within one period of the alarm's at every own-δ k and with δ(p) too at 1 ({nm}) — there the alarm "
            "is the coefficient's, not the restart's") in s
    assert (f"; at {nm}: the earlier starts arrive at pT with kernel 1's coefficient offset 0.07540 with their own δ (0.00130 below δ(p) = "
            "0.07670 and 0.00093 … 0.00094 below their own thresholds at pT) and 0.07760 with δ(p), each above δ(p) reproducing the alarm and "
            "each below passing; the cold restart with each start's arrival coefficient leaves within one period of the alarm's where it is "
            "above δ(p) and passes where below ('cold_m'); kernel 1's threshold coefficient at pT (a leave in the window; rounded outward) in "
            "(0.0763500, 0.0763600] on the cold restart's state; in (0.0763299, 0.0763401] on the own-δ starts' states (k = 1)") in s   # 128, 134
    s = run(dict(co, pT=[pz(1, 'own', 0.0770), pz(1, 'p', 0.0776)]))         # an own-δ start above δ(p) that passes: no rule
    assert "with their own δ above δ(p) = 0.07670 and" in s and "each above δ(p) reproducing" not in s and "('cold_m')" not in s
    s = run(dict(co, cold_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31]], thr={}))   # the cold restart leaving at a low arrival
    assert "each above δ(p) reproducing the alarm" in s and "('cold_m')" not in s and "threshold coefficient" not in s
    s = run(dict(co, warm_m=[[1, 'own', 2, 0.31], [1, 'p', -1, 0.1]]))
    assert "at every own-δ k at 1 (" in s and "with δ(p) too" not in s
    s = run(row(warm_k=[[2, -1, 0.1]], warm_m=[[1, 'own', 2, 0.31], [2, 'own', -1, 0.1]]), row(warm_m=[[1, 'p', -1, 0.1]]))
    assert f"passing at every own-δ k at 1 — there the alarm is the cold restart's —, mixed over k at 1 ({nm}), not run at 1 ({nm})" in s
    assert "coefficient's" not in s
    k1 = row(pT=[pz(1, 'own', 0.07, (-0.07, 0.0, 0.0))], dp=0.07, warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05]])
    k3 = row(pT=[pz(1, 'own', 0.07, (0.01, 0.09, 0.039))], dp=0.07, warm_c=[[1, 0, -1, 0.05], [1, 1, -1, 0.05], [1, 2, 2, 0.85]])
    kb = row(pT=[pz(1, 'own', 0.07, (-0.01, 0.0, 0.0))], dp=0.07, warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05]])
    kk = row(warm_k=[[2, -1, 0.1]], pT=[pz(1, 'own', 0.07), pz(2, 'own', 0.07)], dp=0.07,
             warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05], [2, 0, -1, 0.05], [2, 1, -1, 0.05], [2, 2, -1, 0.05]])
    s = run(k1, k3, kb, kk)
    assert ("; at 4 of those 5, with one formed centre set to its feature as in the cold restart ('warm_c'), every own-δ start's state at pT "
            "leaves within one period of the alarm's at 3 by a centre — kernel 1's centre at 2 (their kernel 1 centre at pT against the merged "
            f"scan: on its removing side at 1, inside its bracket at 1 ({nm})), kernel 3's centre at 1 ({nm}: their kernel 3 centre at pT +0.0390, "
            f"their kernel 1's +0.0100, on the merged scan's keeping side), by no single centre at 1 ({nm})") in s
    s = run(row(warm_m=[[1, 'own', -1, 0.1], [1, 'p', 2, 0.31]]))              # G31: a δ(p) start's state with δ(p) leaving decides nothing
    assert "passing at every own-δ k at 2 — there the alarm is the cold restart's" in s and "neither" not in s
    s = run(row(warm_m=[[1, 'own', 4, 0.31]]))                                 # G32: two periods after the alarm's: a pass, not the alarm
    assert "passing at every own-δ k at 2 — there the alarm is the cold restart's" in s and "coefficient's" not in s
    s = run(row(warm_m=[[1, 'own', 3, 0.31]]))                                 # G33: one period after it: the alarm, not a pass
    assert "passing at every own-δ k at 1 — there the alarm is the cold restart's —, leaving within one period of the alarm's at every own-δ k at 1" in s
    s = run(dict(co, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0760)]))          # G34: a start below δ(p) reproducing it: no rule
    assert "arriving at pT with kernel 1's coefficient offset 0.07600 below δ(p) = 0.07670" in s and "each above δ(p) reproducing" not in s
    kx = row(warm_k=[[2, -1, 0.1]], pT=[pz(1, 'own', 0.07, (-0.07, 0.0, 0.0)), pz(2, 'own', 0.07, (-0.01, 0.0, 0.0))], dp=0.07,
             warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05], [2, 0, 2, 0.85], [2, 1, -1, 0.05], [2, 2, -1, 0.05]])
    assert f"against the merged scan: on both sides at 1 ({nm}))" in run(kx)   # G39: the starts' centres on either side of the bracket
    s = run(row(ext=dict(n=12, cold=[0.82, 0.83], b=[[1, 'own', 0.007, 0.8]], m=[])))   # G47: crossing the cut past the window
    assert "the earlier starts with their own δ at 0.0070 … 0.8000 (k = 1)" in s and "below the cut" not in s
    p1 = run(row(warm=[5, 0.3001], warm_p=[5, 0.3001]))                       # P1: identical outcomes, each start its own label
    assert "1 period earlier with its own δ(p − 1): 4 periods after" in p1 and "1 period earlier with δ(p): 4 periods after" in p1
    p2 = run(row(warm_k=[[2, 6, 0.3001]], warm_pk=[[2, 6, 0.3001]]))           # P2
    assert "2 periods earlier with δ(p): 4 periods after" in p2 and "2 periods earlier with its own δ(p − 2): 4 periods after" in p2
    assert "own drift" not in p2
    p3 = lambda lm: run(row(warm_pk=[[2, 6, 0.3001]], warm_k=[[3, 7, 0.3002]], dp=0.07, pT=[pz(2, 'p', 0.071), pz(3, 'own', 0.071)],
                            warm_m=[[1, 'own', -1, 0.1], [2, 'p', -1, 0.1], [3, 'own', lm, 0.31]]))   # P3: one late start not removed
    assert "own drift" not in p3(6) and "the start's own drift, not the alarm" in p3(-1)
    p4 = run(row(warm_pk=[[2, 4, 0.4]], warm_k=[[2, 4, 0.4]]))                 # P4: an own-δ and a δ(p) start both near
    assert "with the window's δ(p) instead" not in p4 and f"at 1 of 2 ({nm})" in p4
    assert "(merged into the scan; 1 of the 3 in a bracket inside it" in run(row(c1p=-0.01))   # P6: 'c1p' inside the bracket
    flat = [[-0.06, 10, 0.31], [-0.02, 10, 0.31], [0.0, 10, 0.31], [0.04, 10, 0.31]]
    p7 = {**mk(10), 'scan': flat, 'warm': [-1, 0.1], 'c1i': -0.03, 'c1p': -0.02, 'r+c1m': [10, 0.3], 'r+c1i': [10, 0.3], 'r+c1p': [10, 0.3]}
    assert "(merged into the scan; 0 of the 1 in a bracket inside it" in run(p7)          # P7: no bracket, its samples not counted
    ex = row(ext=dict(n=12, cold=[0.82, 0.83], b=[[1, 'own', 0.007, 0.03], [1, 'p', 0.008, 0.029], [2, 'p', 0.8, 0.81]],
                      m=[[1, 'own', 0.82, 0.83], [1, 'p', 0.821, 0.831]]))
    assert (f"; past the window's end, run to p + R + 2 where that passes it ({nm}, the window's periods 11, 12: the cold restart at 0.8200, "
            "0.8300, the earlier starts with their own δ at 0.0070 … 0.0300 (k = 1: below the cut), with δ(p) at 0.0080 … 0.0290 (k = 1: below "
            "the cut) and 0.8000 … 0.8100 (k = 2), from pT with kernel 1's coefficient set to δ(p) ('warm_m') at 0.8200 … 0.8310)") in run(ex)
    assert "past the window's end" not in run(row())


def test_round124_on_the_raw():
    """Referee round 124 on the raw (D3ds-1 … D3ds-5): the own-δ starts are C8D's windows at p − k run k periods longer — the same leave at
    all 36, the largest distance the same to 1e-12 at 28, to 1e-8 at 2, larger (in the extra periods) at 6, δ(p) C8's —; the two segments
    through pT reproduce every earlier start; at pT with δ(p) every start passes at the 11 windows at k_e = 20 and leaves at draw 7's alarm
    period at every k and either δ; the md prints the one-coordinate clause from the 11, the bracket count of the 33, the classification,
    the thresholds, the centres and the period after the window."""
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    key = lambda r: (r['m'], r['ke'], r['metric'], r['s'], r['opts'].get('seed', 7))
    d8 = {key(r): {x['p']: x for x in r['wins']} for r in raw['C8D']}
    n = [0, 0, 0]
    for r in raw['C8A']:
        for x in r['rows']:
            wins, p = d8[key(r)], x['p']
            assert wins[p]['left'] == x['reduced'] and x['dp'] == wins[p]['d']
            for k, y in [(1, x['warm'])] + [(z[0], z[1:]) for z in x['warm_k']]:
                q = wins[p - k]
                d = y[1] - q['maxdev']
                assert y[0] == q['left'] == -1 and (abs(d) < 1e-8 or d > 0)
                n[0 if abs(d) < 1e-12 else 1 if abs(d) < 1e-8 else 2] += 1
            sts = [(1, 'own', x['warm']), (1, 'p', x['warm_p'])] + [(z[0], 'p', z[1:]) for z in x['warm_pk']] + [(z[0], 'own', z[1:]) for z in x['warm_k']]
            for k, kd, y in sts:
                b = next(z['b'] for z in x['pT'] if (z['k'], z['kd']) == (k, kd))
                assert b[0] == (y[0] - k if y[0] >= 0 else -1) and y[1] >= b[1] - 1e-9
            if key(r) == (12.0, 10.0, 'paradigm', 0.02, 7):
                assert [y[:3] for y in x['warm_m']] == [[k, kd, 10] for k in (1, 2, 3) for kd in ('own', 'p')]
                assert all(abs(y[3] - 0.31036) < 1e-5 for y in x['warm_m']) and 0.0763 < x['thr']['cold'][0] < x['thr']['cold'][1] < 0.0764
                assert [y[0] for y in x['thr']['warm']] == [1, 2, 3] and all(0.0763 < y[1] < y[2] < x['thr']['cold'][0] for y in x['thr']['warm'])
            else:
                assert all(y[2] < 0 for y in x['warm_m']) and x['thr'] == {} and 'ext' not in x
    assert n == [28, 2, 6]
    md = NR.report(**{k: raw.get(k) for k in NR.KEYS})
    assert ("at 11 of 12 windows (not at m = 12, k_e = 10, s = 0.02, draw 7 from 107: no offset removes it, its largest distance 0.3104) — at "
            "the 11 kernel 1's centre at the restart is one coordinate against a per-window threshold") in md
    assert "(merged into the scan; 4 of the 33 in a bracket inside it, untested there)" in md
    assert ("passing at every own-δ k at 11 — there the alarm is the cold restart's —, leaving within one period of the alarm's at every own-δ "
            "k and with δ(p) too at 1 (m = 12, k_e = 10, s = 0.02, draw 7 from 107) — there the alarm is the coefficient's, not the restart's") in md
    assert "each above δ(p) reproducing the alarm and each below passing" in md and "kernel 3's centre at 1 (m = 12, k_e = 20, s = 0.06, draw 2 from 58:" in md
    assert ("kernel 1's centre at 10 (their kernel 1 centre at pT against the merged scan: on its removing side at 9, inside its bracket at 1 "
            "(m = 12, k_e = 20, s = 0.08, draw 1 from 50))") in md
    assert "the window's periods 11, 12: the cold restart at 0.8218, 0.8341" in md and "above the cut by 6.1e-06" in md and "above the cut by 3.6e-04" in md
    assert "the start's own drift, not the alarm" in md and "own drift, not the alarm (referee round 119)" not in md


def _c8a_plant_rows():
    """The synthetic C8A rows of the round-124 tests: a window from p = 10 with the reduced leave at 2, its scan with both outcomes."""
    mk = lambda red, **kw: dict(p=10, reduced=red, c1=-0.05, **{t: [kw.get(t, red), 0.4] for t in ('w', 'wc', 'w+c1', 'w+c23', 'i', 'wc+i1c', 'wc+i1w', 'wc+i23')})
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1), src='rerun', left=None, probe=None)
    good = [[-0.06, -1, 0.1], [-0.02, -1, 0.1], [0.0, 2, 0.4], [0.04, 2, 0.4]]
    row = lambda **kw: {**mk(2), 'scan': good, 'warm': [-1, 0.1], 'warm_m': [[1, 'own', -1, 0.1]], 'r+c1m': [-1, 0.1], 'r+c1i': [-1, 0.1],
                        'r+c1p': [-1, 0.1], **kw}
    run = lambda *xs: '\n'.join(NR.c8a_lines([dict(r, rows=list(xs) + [row()])]))
    pz = lambda k, kd, w1, c=(-0.05, 0.0, 0.0): dict(k=k, kd=kd, w1=w1, c=list(c), b=[-1, 0.1])
    return row, run, pz, "m = 12, k_e = 20, s = 0.06, draw 1 from 10"


def test_round128_brackets_outward_the_matched_arrival_everywhere_and_the_reverse(monkeypatch):
    """Referee round 128 (D3dw-1 … D3dw-3, D3dw-5, D3dw-6): the brackets rounded outward (`bracket_out`); the matched arrival at every
    window that has it, a reproduction whose matched arrival passes printed as the arrival's (plant Q1); "mixed over k" (plant Q2); the
    reverse centre transplant 'cold_c' — necessary and sufficient only where 'warm_c' and 'cold_c' name the same centre at every own-δ k;
    the margins to the own thresholds, the run's own coefficient at pT ('w1i'), no empty "by a centre" clause, the header's rounds;
    `c8a_job` stores 'w1i' from the rerun at pT."""
    assert NR.bracket_out(0.07635111802815617, 0.07635121486619703) == "(0.0763511, 0.0763513]"   # D3dw-1: to nearest it excludes
    assert NR.bracket_out(0.0763370644406115, 0.0763373467359406) == "(0.0763370, 0.0763374]"
    assert NR.bracket_out(0.07635, 0.07636) == "(0.0763500, 0.0763600]" and NR.bracket_out(0.5, 0.75, 1) == "(0.5, 0.8]"
    from decimal import Decimal                                          # round 134 (D3ec-1): against the doubles' exact values
    rng = np.random.default_rng(128)
    for lo in rng.uniform(0.07, 0.08, 200):
        hi = lo + rng.uniform(0, 1e-6)
        plo, phi = (Decimal(v) for v in NR.bracket_out(lo, hi)[1:-1].split(', '))
        assert plo <= Decimal(float(lo)) and Decimal(float(hi)) <= phi and phi - plo < Decimal(hi - lo) + Decimal('2.0001e-7')
    row, run, pz, nm = _c8a_plant_rows()
    q1 = row(warm=[3, 0.4], dp=0.07, pT=[pz(1, 'own', 0.075)], warm_m=[[1, 'own', -1, 0.1]], cold_m=[[1, 'own', 2, 0.4]])   # Q1
    s = run(q1)
    assert (f"leaves within one period of the alarm's at 1 of 2 ({nm}) — the arrival's, not the alarm, at 1 of them ({nm}): from its state "
            "at pT with kernel 1's coefficient set to δ(p) ('warm_m') each reproducing start passes") in s
    assert ("at the 2 windows with an own-δ start's state at pT (1 where no own-δ start leaves near or before the alarm's period, 1 where one "
            "reproduces the alarm): passing at every own-δ k at 2 — there the alarm is the cold restart's") in s
    s = run(dict(q1, warm_m=[[1, 'own', 3, 0.4]]))                       # the matched arrival reproduces it too: the coefficient's
    assert "the arrival's, not the alarm" not in s and f"leaving within one period of the alarm's at every own-δ k at 1 ({nm}) — there the alarm is the coefficient's" in s
    s = run(dict(q1, warm_m=[]))                                         # no matched arrival: neither printed nor counted
    assert "the arrival's" not in s and "at the 1 window where no own-δ start leaves near" in s
    q2 = row(warm_k=[[2, -1, 0.1]], dp=0.07, pT=[pz(1, 'own', 0.069), pz(2, 'own', 0.069)], warm_m=[[1, 'own', 2, 0.4], [2, 'own', -1, 0.1]])
    s = run(q2)                                                           # Q2 (D3dw-5): the code's mixed class, not round 111's "neither"
    assert f"mixed over k at 1 ({nm})" in s and "neither" not in s
    assert "(C8A, referee rounds 98, 102, 105, 111, 114, 119, 124, 128, 134): each of the 2 windows" in s
    cen = lambda **kw: row(dp=0.07, pT=[pz(1, 'own', 0.07, (-0.07, 0.0, 0.0))], warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05]], **kw)
    s = run(cen(cold_c=[[1, 0, -1, 0.05], [1, 1, 2, 0.85], [1, 2, 2, 0.85]]))   # D3dw-3: both directions name kernel 1's centre
    assert ("by a centre — kernel 1's centre at 1 (their kernel 1 centre at pT against the merged scan: on its removing side at 1); in reverse, "
            "with that centre at each own-δ start's value at pT instead of at its feature, the cold restart passes at every own-δ k ('cold_c') "
            "at all 1 — there the centre is necessary and sufficient, no other centre alone making the cold restart pass") in s
    s = run(cen(cold_c=[[1, 0, -1, 0.05], [1, 1, -1, 0.05], [1, 2, 2, 0.85]]))   # kernel 2's centre makes it pass too
    assert "necessary and sufficient" in s and "no other centre alone" not in s
    s = run(cen(cold_c=[[1, 0, 3, 0.85], [1, 1, 2, 0.85], [1, 2, 2, 0.85]]))      # the reverse keeps the alarm: sufficiency only
    assert f"('cold_c') at 0 of those 1 (not at {nm})" in s and "necessary and sufficient" not in s
    two = row(warm_k=[[2, -1, 0.1]], dp=0.07, pT=[pz(1, 'own', 0.07, (-0.07, 0.0, 0.0)), pz(2, 'own', 0.07, (-0.07, 0.0, 0.0))],
              warm_m=[[1, 'own', -1, 0.1], [2, 'own', -1, 0.1]], warm_c=[[k, j, 2 if j == 0 else -1, 0.5] for k in (1, 2) for j in (0, 1, 2)])
    s = run(dict(two, cold_c=[[k, j, -1 if (k, j) == (1, 0) else 2, 0.5] for k in (1, 2) for j in (0, 1, 2)]))   # at k = 1 only
    assert "('cold_c') at 0 of those 1" in s and "necessary and sufficient" not in s
    s = run(dict(two, cold_c=[[k, j, -1 if j == 0 else 2, 0.5] for k in (1, 2) for j in (0, 1, 2)]))
    assert "('cold_c') at all 1 — there the centre is necessary and sufficient" in s
    assert "'cold_c'" not in run(cen())                                  # no reverse run: nothing said about necessity
    z22 = row(dp=0.07, pT=[pz(1, 'own', 0.07)], warm_c=[[1, 0, 7, 0.4], [1, 1, -1, 0.05], [1, 2, -1, 0.05]])   # D3dw-6: no centre
    s = run(z22)
    assert ("('warm_c'), no single centre makes every own-δ start's state at pT leave within one period of the alarm's, at 1 "
            f"({nm})") in s and "by a centre" not in s and "— ," not in s
    co = row(warm_p=[3, 0.31], warm_k=[[2, -1, 0.1]], dp=0.0767, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0776), pz(2, 'own', 0.0755)],
             warm_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31], [2, 'own', 2, 0.31]], cold_m=[[1, 'own', -1, 0.03], [1, 'p', 2, 0.31], [2, 'own', -1, 0.03]],
             cold_c=[[k, j, 2, 0.31] for k in (1, 2) for j in (0, 1, 2)], w1i=0.07613,
             thr=dict(cold=[0.0763511, 0.0763512], warm=[[1, 0.0763377, 0.0763379]]))
    s = run(co)                                                           # the coefficient's: margins (k = 1 only), 'cold_c', 'w1i'
    assert ("with their own δ (0.00120 … 0.00130 below δ(p) = 0.07670 and 0.00094 below their own thresholds at pT (k = 1)) and 0.07760 "
            "with δ(p), each above δ(p) reproducing the alarm and each below passing") in s
    assert ("passes where below ('cold_m'); the cold restart with one own-δ start's formed centre at its value at pT leaves within one period "
            "of the alarm's at every k and kernel ('cold_c'); kernel 1's threshold coefficient at pT (a leave in the window; rounded outward) in "
            "(0.0763511, 0.0763512] on the cold restart's state; in (0.0763376, 0.0763380] on the own-δ starts' states (k = 1); the run's own "
            "kernel 1 coefficient offset at pT 0.07613 ('w1i'), below every threshold") in s
    s = run(dict(co, w1i=0.0764, cold_c=co['cold_c'][:-1] + [[2, 2, -1, 0.03]]))
    assert "0.07640 ('w1i')" in s and "below every threshold" not in s and "('cold_c')" not in s
    s = run(dict(co, thr=dict(cold=[0.0763511, 0.0763512], warm=[[1, 0.0763377, 0.0763379], [2, 0.0763376, 0.0763378]])))
    assert "and 0.00084 … 0.00094 below their own thresholds at pT) and" in s
    L = NR.law(12.0)                                                      # D3dw-6: 'w1i' from the rerun's w at pT, kernel 1
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    rr = dict(w=np.zeros((97, 21)), c=np.zeros((97, 21)))
    rr['w'][40, idx[0]] = -(w0 + 0.0761)                                 # output 8p, p = 5: its magnitude counts
    monkeypatch.setattr(NR, 'win_run', lambda *a, **kw: (3, 0.4))
    monkeypatch.setattr(NR, 'c_run', lambda m, ke, s, opts: (rr, None, None, None, list(idx), np.full((12, 3), 0.5), np.tile(CT.CSTAR, (12, 1))))
    monkeypatch.setattr(NR, 'c8a_pT', lambda *a: {})
    out = NR.c8a_job((12.0, 10.0, 'paradigm', 0.02, {'seed': 7}, [(5, 2)], None, None))
    assert abs(out['rows'][0]['w1i'] - 0.0761) < 1e-12


def test_round128_the_kill_plants_of_round_128s_survivors():
    """Referee round 128 (D3dw-4): one plant per surviving mutant of the round — Z16 ('rel' above with any: two δ(p) reproductions on either
    side of δ(p) are "against" it), Z19 (late from R + k + 1: a δ(p) start leaving at R + k + 1 is a reproduction, not a late leave), Z22
    (a centre explains a window only by a leave near the alarm's period, not by any leave), Z26 (the 'cold_m' rule over every start, the
    δ(p) starts' too)."""
    row, run, pz, nm = _c8a_plant_rows()
    s = run(row(warm_p=[3, 0.31], warm_pk=[[2, 4, 0.31]], dp=0.07, pT=[pz(1, 'p', 0.071), pz(2, 'p', 0.069)]))          # Z16
    assert (f"with the window's δ(p) instead at 1 more ({nm}: 1, 2 periods earlier, arriving at pT with kernel 1's coefficient offset 0.07100, "
            "0.06900 against δ(p) = 0.07000)") in s
    s = run(row(warm_pk=[[3, 6, 0.31]]))                                                                                     # Z19
    assert f"with the window's δ(p) instead at 1 more ({nm}: 3 periods earlier)" in s and "after the alarm's period" not in s
    s = run(row(dp=0.07, pT=[pz(1, 'own', 0.07)], warm_c=[[1, 0, 7, 0.4], [1, 1, -1, 0.05], [1, 2, -1, 0.05]]))            # Z22
    assert "no single centre makes every own-δ start's state at pT leave" in s and "kernel 1's centre at 1" not in s
    s = run(row(warm_p=[-1, 0.1], dp=0.0767, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0776)], warm_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31]],
                cold_m=[[1, 'own', -1, 0.03], [1, 'p', -1, 0.03]]))                                                         # Z26
    assert "there the alarm is the coefficient's, not the restart's" in s and "('cold_m')" not in s


def test_round128_on_the_raw():
    """Referee round 128 on the raw (D3dw-1, -3, -6): every printed threshold bracket contains the stored one; 'cold_c' — the cold restart
    with one own-δ start's formed centre at its value at pT — passes at every k with exactly the centre 'warm_c' names (kernel 1's at 10
    windows, kernel 3's at draw 2 from 58) and leaves at the alarm's period with each other one, and at draw 7 leaves at the alarm's period
    at every k and kernel; 'w1i' is C's loop's own |w₁| − w₀ at pT (0.076126 at draw 7); the md prints the margins to the own thresholds,
    the 2² cells, the reverse's verdict and 'w1i'."""
    import re
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    md = NR.report(**{k: raw.get(k) for k in NR.KEYS})
    key = lambda r: (r['ke'], r['s'], r['opts'].get('seed', 7))
    rows = {key(r) + (x['p'],): x for r in raw['C8A'] for x in r['rows']}
    thr = [x for x in rows.values() if x.get('thr')]
    assert [k_ for k_, x in rows.items() if x.get('thr')] == [(10.0, 0.02, 7, 107)]
    m = re.findall(r"threshold coefficient at pT \(a leave in the window; rounded outward\) in \((0\.\d+), (0\.\d+)\] on the cold restart's state; "
                   r"in ((?:\(0\.\d+, 0\.\d+\](?:, )?)+) on the own-δ starts' states \(k = ([\d, ]+)\)", md)
    assert len(m) == 1
    lo, hi, warm, ks = m[0]
    th = thr[0]['thr']
    assert float(lo) <= th['cold'][0] < th['cold'][1] <= float(hi) and (lo, hi) == ('0.0763511', '0.0763513')   # D3dw-1
    pw = [tuple(float(v) for v in b_.split(', ')) for b_ in re.findall(r"\((0\.\d+, 0\.\d+)\]", warm)]
    assert [int(k_) for k_ in ks.split(', ')] == [y[0] for y in th['warm']] and len(pw) == len(th['warm'])
    assert all(a_ <= y[1] < y[2] <= b_ for (a_, b_), y in zip(pw, th['warm']))
    for k_, x in rows.items():                                           # D3dw-3: both directions, one centre at a time
        R = x['reduced']
        cc, wc = x['cold_c'], x['warm_c']
        assert [(y[0], y[1]) for y in cc] == [(k, j) for k in (1, 2, 3) for j in (0, 1, 2)] == [(y[0], y[1]) for y in wc]
        if k_ == (10.0, 0.02, 7, 107):
            assert all(y[2] == R for y in cc)
            continue
        j0 = 2 if k_ == (20.0, 0.06, 2, 58) else 0
        assert all((y[2] < 0) == (y[1] == j0) for y in cc) and all(y[2] == R for y in cc if y[1] != j0)
        assert all((y[2] == R) == (y[1] == j0) for y in wc) and all(y[2] < 0 for y in wc if y[1] != j0)
    assert abs(rows[(10.0, 0.02, 7, 107)]['w1i'] - 0.076126) < 5e-7 and all(isinstance(x['w1i'], float) for x in rows.values())
    assert "with their own δ (0.00115 … 0.00122 below δ(p) = 0.07667 and 0.00082 … 0.00089 below their own thresholds at pT)" in md
    assert "the run's own kernel 1 coefficient offset at pT 0.07613 ('w1i'), below every threshold" in md
    assert ("at those 11 the cold restart with each own-δ start's arrival coefficient leaves within one period of the alarm's at every own-δ k "
            "('cold_m') and the starts continued from pT pass") in md
    assert ("the cold restart passes at every own-δ k ('cold_c') at all 11 — there the centre is necessary and sufficient, no other centre "
            "alone making the cold restart pass") in md
    assert "the cold restart with one own-δ start's formed centre at its value at pT leaves within one period of the alarm's at every k and kernel ('cold_c')" in md
    assert "(C8A, referee rounds 98, 102, 105, 111, 114, 119, 124, 128, 134)" in md and "by a centre — ," not in md


def test_round128_plants_for_the_rounds_own_mutants():
    """Round 128's own mutants (tools/mutants/round128_D60_own.py) that the other tests leave alive, one plant each: H14 (a reproduction
    is the arrival's only where every reproducing own-δ start's matched arrival passes), H21 ("no other centre alone" only where that
    holds at every such window), H26 ('w1i' below every threshold, the cold restart's too), H27 (the 2² clause only where the continued
    starts pass), H33 (the reverse counted over the windows a centre explains)."""
    row, run, pz, nm = _c8a_plant_rows()
    s = run(row(warm=[3, 0.4], warm_k=[[2, 4, 0.4]], dp=0.07, pT=[pz(1, 'own', 0.075), pz(2, 'own', 0.075)],
                warm_m=[[1, 'own', -1, 0.1], [2, 'own', 2, 0.4]]))      # H14: k = 1 the arrival's, k = 2 the alarm with δ(p) too
    assert f"leaves within one period of the alarm's at 1 of 2 ({nm})" in s and "the arrival's, not the alarm" not in s and f"mixed over k at 1 ({nm})" in s
    cw = lambda j: [[1, q, 2 if q == j else -1, 0.5] for q in (0, 1, 2)]
    a = row(dp=0.07, pT=[pz(1, 'own', 0.07)], warm_c=cw(0), cold_c=[[1, 0, -1, 0.1], [1, 1, 2, 0.5], [1, 2, 2, 0.5]])
    b = dict(a, cold_c=[[1, 0, -1, 0.1], [1, 1, -1, 0.1], [1, 2, 2, 0.5]])   # kernel 2's centre makes it pass too
    s = run(a, b)                                                        # H21
    assert "('cold_c') at all 2 — there the centre is necessary and sufficient" in s and "no other centre alone" not in s
    co = row(warm_p=[3, 0.31], dp=0.0767, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0776)], warm_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31]],
             cold_m=[[1, 'own', -1, 0.03], [1, 'p', 2, 0.31]], thr=dict(cold=[0.0763, 0.07631], warm=[[1, 0.0764, 0.07641]]), w1i=0.07635)
    s = run(co)                                                          # H26: between the cold restart's threshold and the own-δ one
    assert "0.07635 ('w1i')" in s and "below every threshold" not in s
    s = run(row(warm=[3, 0.4], dp=0.07, pT=[dict(pz(1, 'own', 0.075), b=[2, 0.4])], warm_m=[[1, 'own', -1, 0.1]], cold_m=[[1, 'own', 2, 0.4]]))
    assert "there the alarm is the cold restart's" in s and "('cold_m') and the starts continued from pT pass" not in s   # H27
    s = run(a, row(dp=0.07, pT=[pz(1, 'own', 0.07)], warm_c=cw(None), cold_c=[[1, q, 2, 0.5] for q in (0, 1, 2)]))   # H33
    assert f"by no single centre at 1 ({nm})" in s and "('cold_c') at all 1 — there the centre is necessary and sufficient" in s


def test_round134_brackets_from_the_exact_double():
    """Referee round 134 (D3ec-1): `bracket_out` rounds outward from each double's exact value, not from its shortest repr — on the
    7-decimal grid the repr lies inward of the double about half the time (0.0763377 is 0.07633769999…: its lo prints 0.0763376); the
    printed ends contain the exact doubles at every grid point and at their neighbours, and an exactly representable end prints as itself."""
    from decimal import Decimal
    assert NR.bracket_out(0.0763377, 0.0763379) == "(0.0763376, 0.0763380]" and NR.bracket_out(0.07633, 0.07634) == "(0.0763299, 0.0763401]"
    assert NR.bracket_out(0.5, 0.75, 2) == "(0.50, 0.75]" and NR.bracket_out(-0.0763377, -0.0763376) == "(-0.0763377, -0.0763376]"
    inw = 0
    for n in range(763300, 763500):
        for v in (n / 1e7, np.nextafter(n / 1e7, 0.0), np.nextafter(n / 1e7, 1.0)):
            plo, phi = (Decimal(u) for u in NR.bracket_out(v, v)[1:-1].split(', '))
            e = Decimal(float(v))
            assert plo <= e <= phi and phi - plo <= Decimal('2e-7')
        inw += Decimal(repr(n / 1e7)) != Decimal(float(n / 1e7)).quantize(Decimal('1e-7'), rounding='ROUND_FLOOR')
    assert inw > 50                                                      # the shortest repr would have printed lo inward there


def test_round134_kill_plants_and_the_new_clauses():
    """Referee round 134 (D3ec-2, -4, -5): the referee's kill plants P05–P10 of R05–R10 — the 2² clause only where every continued start
    passes (P05) and 'cold_m' leaves at every own-δ k (P06), "necessary and sufficient" only where 'cold_c' names every centre 'warm_c'
    names (P07), the margins only where none is negative (P08), "below every threshold" only below every bracket's lower end (P09, P10)
    —; R04's convention kept: a pass in 'cold_c' is past the alarm's period, a later leave in the window included and counted; the 2²
    clause's contrast in kernel 1's coefficient and 'w1i' against δ(p); the passes' largest distance and its margin to the cut."""
    row, run, pz0, nm = _c8a_plant_rows()
    pz = lambda k, kd, w1, c=(-0.05, 0.0, 0.0), b=(-1, 0.1): dict(k=k, kd=kd, w1=w1, c=list(c), b=list(b))
    two = dict(dp=0.07, warm_k=[[2, -1, 0.1]], warm_m=[[1, 'own', -1, 0.1], [2, 'own', -1, 0.1]])
    s = run(row(**two, pT=[pz(1, 'own', 0.07), pz(2, 'own', 0.07, b=(2, 0.4))], cold_m=[[1, 'own', 2, 0.4], [2, 'own', 2, 0.4]]))   # P05
    assert "passing at every own-δ k at 2 — there the alarm is the cold restart's" in s and "a 2² factor separation" not in s
    s = run(row(**two, pT=[pz(1, 'own', 0.07), pz(2, 'own', 0.07)], cold_m=[[1, 'own', 2, 0.4], [2, 'own', -1, 0.1]]))              # P06
    assert "passing at every own-δ k at 2 — there the alarm is the cold restart's" in s and "a 2² factor separation" not in s
    s = run(row(dp=0.07, pT=[pz(1, 'own', 0.07, (-0.07, 0.05, 0.0))], warm_c=[[1, 0, 2, 0.8], [1, 1, 2, 0.8], [1, 2, -1, 0.05]],
                cold_c=[[1, 0, -1, 0.05], [1, 1, 2, 0.8], [1, 2, 2, 0.8]]))                                                         # P07
    assert f"('cold_c') at 0 of those 1 (not at {nm})" in s and "necessary and sufficient" not in s
    co = row(warm_p=[3, 0.31], warm_k=[[2, -1, 0.1]], dp=0.0767, pT=[pz(1, 'own', 0.0754), pz(1, 'p', 0.0776), pz(2, 'own', 0.0755)],
             warm_m=[[1, 'own', 2, 0.31], [1, 'p', 2, 0.31], [2, 'own', 2, 0.31]], cold_m=[[1, 'own', -1, 0.03], [1, 'p', 2, 0.31], [2, 'own', -1, 0.03]],
             cold_c=[[k, j, 2, 0.31] for k in (1, 2) for j in (0, 1, 2)], w1i=0.07613, thr=dict(cold=[0.0763511, 0.0763512], warm=[[1, 0.0763377, 0.0763379]]))
    s = run(dict(co, thr=dict(cold=[0.0763511, 0.0763512], warm=[[1, 0.0763377, 0.0763379], [2, 0.0754, 0.07541]])))             # P08
    assert "with their own δ (0.00120 … 0.00130 below δ(p) = 0.07670) and 0.07760 with δ(p)" in s and "below their own thresholds" not in s
    for th in (dict(cold=[0.0763, 0.0764], warm=[[1, 0.0765, 0.0766]]), dict(cold=[0.0770, 0.0771], warm=[[1, 0.0763, 0.0764]])):   # P09, P10
        s = run(dict(co, w1i=0.07635, thr=th))
        assert "offset at pT 0.07635 ('w1i')." in s and "below every threshold" not in s
    cen = lambda **kw: row(dp=0.07, pT=[pz(1, 'own', 0.07, (-0.07, 0.0, 0.0))], warm_c=[[1, 0, 2, 0.85], [1, 1, -1, 0.05], [1, 2, -1, 0.05]], **kw)
    s = run(cen(cold_c=[[1, 0, 5, 0.35], [1, 1, 2, 0.85], [1, 2, 2, 0.85]]))   # R04's convention: a leave after the alarm's period passes
    assert ("('cold_c') at all 1 — there the centre is necessary and sufficient, no other centre alone making the cold restart pass (a pass "
            "is past the alarm's period: 1 of the 1 runs with the named centre leave later in the window)") in s
    s = run(cen(cold_c=[[1, 0, -1, 0.05], [1, 1, 2, 0.85], [1, 2, 2, 0.85]]))
    assert ("making the cold restart pass (a pass is past the alarm's period: none of the 1 runs with the named centre leaves later in the "
            f"window, the largest distance of those that do not leave 0.0500 at {nm}, 0.2500 below the cut)") in s
    s = run(cen(cold_c=[[1, 0, 3, 0.35], [1, 1, 2, 0.85], [1, 2, 2, 0.85]]))   # a leave within one period of the alarm's: no pass
    assert "necessary and sufficient" not in s and "(a pass is" not in s
    q = row(dp=0.07, pT=[pz(1, 'own', 0.0695)], warm_m=[[1, 'own', -1, 0.1]], cold_m=[[1, 'own', 2, 0.4]], w1i=0.0712)   # D3ec-4
    s = run(q)
    assert ("not with kernel 1's coefficient (its two levels, each own-δ arrival and δ(p), 5.0e-04 apart; the run's own kernel 1 coefficient "
            "offset at pT ('w1i') 0.00120 above δ(p))") in s
    s = run(q, row(dp=0.07, pT=[pz(1, 'own', 0.06)], warm_m=[[1, 'own', -1, 0.1]]))   # a cold-restart window without 'cold_m': not counted
    assert "(its two levels, each own-δ arrival and δ(p), 5.0e-04 apart;" in s
    q2 = dict(q, pT=[pz(1, 'own', 0.0712)], w1i=0.0690)
    s = run(q, q2)
    assert ("(its two levels, each own-δ arrival and δ(p), 5.0e-04 … 1.2e-03 apart; the run's own kernel 1 coefficient offset at pT ('w1i') "
            "-0.00100 … +0.00120 from δ(p))") in s
    s = run(dict(q, w1i=0.0690), {k: v for k, v in q2.items() if k != 'w1i'})
    assert "('w1i') 0.00100 below δ(p) at 1 of them)" in s
    assert "(its two levels, each own-δ arrival and δ(p), 1.2e-03 apart)" in run({k: v for k, v in q2.items() if k != 'w1i'})


def _c8f_row(p, R, k, j, out, deciding=1, dev=None):
    """A synthetic C8F start: the outcome decided by group `deciding` (a leave at R where its bit is 0, none where 1), the largest
    distance `dev` per corner (default: 0.3 + 0.5·(1 − bit))."""
    import itertools
    cells = []
    for bits in itertools.product((0, 1), repeat=5):
        b = ''.join(map(str, bits))
        cells.append([b, R if b[deciding] == '0' else -1, dev(b) if dev else 0.3 + 0.5 * (b[deciding] == '0')])
    return dict(p=p, reduced=R, k=k, j=j, dp=0.07, cells=cells)


def test_round134_c8f_lines_cases_and_job(monkeypatch):
    """Referee round 134 (D3ec-3): `c8f_lines` names the deciding group per start, the Stein–Alpert terms of the largest distance (the
    deciding group's share of the change, the largest other term), "no single group's" where none decides, and the corners against
    C8A's runs; `c8f_cases` takes the windows with states at pT and the centre 'cold_c' names (kernel 1's where none or several);
    `c8f_job` builds the 32 corners from the cold restart and the own-δ start's state at pT, group by group."""
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1))
    nm = "m = 12, k_e = 20, s = 0.06, draw 1 from 10"
    rows = [_c8f_row(10, 2, k, 0, None) for k in (1, 2, 3)] + [_c8f_row(20, 2, 1, 2, None, deciding=0)]
    s = '\n'.join(NR.c8f_lines([dict(r, rows=rows)]))
    assert ("at each of the 2 windows with states at pT and each own-δ k = 1, 2, 3 (4 starts), the cold restart (C8D's window at p) and "
            "the own-δ start's state at pT crossed in five groups — kernel 1's coefficient; the named formed centre; the other two formed "
            "centres; the plant's x; every other α with the unformed kernels' centres — (the named centre kernel 1's at 1 "
            f"({nm}), kernel 3's at 1 (m = 12, k_e = 20, s = 0.06, draw 1 from 20): the one 'cold_c' names, kernel 1's where it names none)") in s
    assert (f"is the named formed centre's alone at 3 starts at 1 window ({nm}), kernel 1's coefficient's alone at 1 start at 1 window "
            "(m = 12, k_e = 20, s = 0.06, draw 1 from 20); on the window's largest distance") in s
    assert "the deciding group's main term is +1.000 of the change from the cold restart to the start continued and every other main or interaction term at most 0.0e+00 of it" in s
    dev = lambda b: 0.3 + 0.5 * (b[1] == '0') + 0.01 * (b[3] == '1') * (b[4] == '1')   # an interaction of x and the rest, D·E
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_c8f_row(10, 2, 1, 0, None, dev=dev)])]))
    assert "main term is +1.020 of the change" in s and "at most 2.0e-02 of it in magnitude" in s
    bad = _c8f_row(10, 2, 1, 0, None)
    bad['cells'][3][1] = 5                                                # one corner leaving after the alarm's: no group decides
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[bad, _c8f_row(10, 2, 2, 0, None)])]))
    assert "the named formed centre's alone at 1 start at 1 window" in s and "no single group's at 1 start at 1 window" in s
    late = _c8f_row(10, 2, 1, 0, None)
    late['cells'][8][1] = 5                                               # a pass leaving after the alarm's period still passes
    assert "the named formed centre's alone at 1 start at 1 window" in '\n'.join(NR.c8f_lines([dict(r, rows=[late])]))
    assert NR.c8f_lines(None) == [] and NR.c8f_lines([]) == []
    c8a = [dict(r, rows=[dict(p=10, reduced=2, cold_m=[[1, 'own', 2, 0.8]], cold_c=[[1, 0, -1, 0.3], [1, 1, 2, 0.8], [1, 2, 2, 0.8]],
                              warm_c=[[1, 0, 2, 0.8], [1, 1, -1, 0.3], [1, 2, -1, 0.3]], pT=[dict(k=1, kd='own', b=[-1, 0.3 + 1e-13])])])]
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_c8f_row(10, 2, 1, 0, None)])], c8a))
    assert ("the corners that are C8A's runs ('cold_m', 'cold_c', 'warm_c', the start continued) leave at C8A's periods and the cold restart "
            "at the alarm's at 1 of the 1 starts, their largest distances within 1.0e-13 of C8A's") in s
    c0_ = _c8f_row(10, 2, 1, 0, None)
    c0_['cells'][0][1] = 3                                                # the cold restart one period after the alarm's
    assert "at the alarm's at 0 of the 1 starts" in '\n'.join(NR.c8f_lines([dict(r, rows=[c0_])], c8a))
    c8a[0]['rows'][0]['cold_m'] = [[1, 'own', 3, 0.8]]
    assert "at the alarm's at 0 of the 1 starts" in '\n'.join(NR.c8f_lines([dict(r, rows=[_c8f_row(10, 2, 1, 0, None)])], c8a))
    # c8f_cases: the named centre from 'cold_c' (a later leave passes), kernel 1's where none or several; p ≥ 2 only
    cc = lambda js: [[k, j, (-1 if j in js else 2), 0.1] for k in (1, 2) for j in (0, 1, 2)]
    a = dict(r, metric='paradigm', rows=[dict(p=1, reduced=2), dict(p=10, reduced=2, cold_c=cc((2,))), dict(p=20, reduced=2, cold_c=cc(())),
                                        dict(p=30, reduced=2, cold_c=cc((1, 2))), dict(p=40, reduced=2, cold_c=[[1, 1, 5, 0.3], [1, 0, 2, 0.8]]),
                                        dict(p=50, reduced=2, cold_c=[[1, 0, 3, 0.3], [1, 1, -1, 0.1]])])   # R + 1: not a pass
    monkeypatch.setattr(NR, 'c8a_cases', lambda *a_: [(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(1, 2), (10, 2), (20, 2), (30, 2), (40, 2), (50, 2)], None, None),
                                                      (12.0, 10.0, 'paradigm', 0.03, {'seed': 5}, [], None, (150, []))])
    out = NR.c8f_cases(None, None, None, None, [a])
    assert out == [(12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(10, 2), (20, 2), (30, 2), (40, 2), (50, 2)], None, {10: 2, 20: 0, 30: 0, 40: 1, 50: 1})]
    # c8f_job: the corners, group by group
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    calls = []

    def seg(m, ke, mt, s_, opts, shift, n, al, c, x0=None):
        if n < NR.NP_C7:                                                  # the start's k periods: a state distinct in every component
            return [0.0] * n, (7.0, np.asarray(al, float) + 1.0, np.asarray(c, float) + 2.0)
        calls.append((shift, n, np.array(al, float), np.array(c, float), x0))
        return [0.1] * n, None
    monkeypatch.setattr(NR, 'seg_run', seg)
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    W = [[w0 + 0.07, 0.5, 0.5]] * 12
    o = NR.c8f_job((12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(3, 2)], (W, None), {3: 2}))
    assert [(x['k'], len(x['cells'])) for x in o['rows']] == [(1, 32), (2, 32)] and all(x['j'] == 2 and x['p'] == 3 for x in o['rows'])
    assert o['rows'][0]['cells'][0] == ['00000', -1, 0.1] and len(calls) == 64
    alc = NR.raised(L, al0, idx, 0, abs(w0 + 0.07) - w0)
    als = alc + 1.0
    un = [i for i in range(c0.size) if i not in set(int(i_) for i_ in idx)]
    oth = [i for i in range(al0.size) if i != idx[0]]
    for (shift, n, al, c, x0), (b, *_) in zip(calls[:32], o['rows'][0]['cells']):
        A, B, C, D, E_ = (ch == '1' for ch in b)
        assert shift == 3 and n == NR.NP_C7 and (x0 == 7.0 if D else x0 is None)
        assert al[idx[0]] == (als if A else alc)[idx[0]] and np.array_equal(al[oth], (als if E_ else alc)[oth])
        assert c[idx[2]] == c0[idx[2]] + 2.0 * B and all(c[idx[j]] == c0[idx[j]] + 2.0 * C for j in (0, 1))
        assert np.array_equal(c[un], c0[un] + 2.0 * E_)


def test_round134_c8f_on_the_raw():
    """Referee round 134 on the raw (D3ec-3, -4, -5): C8F has the 12 windows' 36 own-δ starts; its corners that are C8A's runs give C8A's
    leave periods ('cold_m', 'cold_c', 'warm_c' to the bit, the continued start to 10⁻¹²: C8A's runs it 12 periods at draw 7) and the
    cold restart leaves at the alarm's period; the outcome is the named centre's alone at the 33 starts at k_e = 20 (kernel 3's at draw 2
    from 58) and kernel 1's coefficient's alone at draw 7's 3; the md prints it, the 2² contrast and the passes' margin."""
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    md = NR.report(**{k: raw.get(k) for k in NR.KEYS})
    key = lambda r: (r['ke'], r['s'], r['opts'].get('seed', 7))
    a = {key(r) + (x['p'],): x for r in raw['C8A'] for x in r['rows']}
    f = [(key(r) + (x['p'],), x) for r in raw['C8F'] for x in r['rows']]
    assert len(f) == 36 and sorted({k_ for k_, _ in f}) == sorted(k_ for k_, x in a.items() if 'pT' in x)
    for k_, x in f:
        y, cl, R, k, j = a[k_], {b: v for b, *v in x['cells']}, x['reduced'], x['k'], x['j']
        assert x['dp'] == y['dp'] and R == y['reduced'] and j == (2 if k_ == (20.0, 0.06, 2, 58) else 0) and len(cl) == 32
        assert cl['10000'] == next(z[2:] for z in y['cold_m'] if (z[0], z[1]) == (k, 'own'))
        assert cl['01000'] == next(z[2:] for z in y['cold_c'] if (z[0], z[1]) == (k, j))
        assert cl['10111'] == next(z[2:] for z in y['warm_c'] if (z[0], z[1]) == (k, j))
        b = next(z['b'] for z in y['pT'] if (z['k'], z['kd']) == (k, 'own'))
        assert cl['11111'][0] == b[0] and abs(cl['11111'][1] - b[1]) < 1e-12 and cl['00000'][0] == R
        q = 0 if k_ == (10.0, 0.02, 7, 107) else 1
        assert all((v[0] == R) if bb[q] == '0' else (v[0] < 0) for bb, v in cl.items())
    assert ("is the named formed centre's alone at 33 starts at 11 windows, kernel 1's coefficient's alone at 3 starts at 1 window (m = 12, "
            "k_e = 10, s = 0.02, draw 7 from 107); on the window's largest distance") in md
    assert ("the cold restart at the alarm's at 36 of the 36 starts, their largest distances within ") in md
    assert "(its two levels, each own-δ arrival and δ(p), 1.2e-06 … 1.7e-03 apart; the run's own kernel 1 coefficient offset at pT ('w1i') 0.00055 … 0.00182 above δ(p))" in md
    assert "(a pass is past the alarm's period: none of the 33 runs with the named centre leaves later in the window, the largest distance of those that do not leave 0.2857 at m = 12, k_e = 20, s = 0.08, draw 1 from 50, 0.0143 below the cut)" in md


def test_round138_c8f_job_start_segment_and_corners(monkeypatch):
    """Referee round 138 (D3eg-2): `c8f_job`'s start segment — each own-δ start runs from (p − k)T for k periods with kernel 1 at
    raised(δ(p − k)) and the fitted centres — and its corners, with a W that varies by period (δ(p) from W[p − 1], not a neighbour), a
    distance list whose last period carries the maximum (the largest distance over every period, the list stored as 'dev'), and the
    named centre jb = 0 (group C then sets both other formed centres)."""
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    W = [[w0 + 0.05 + 0.001 * q, 0.5, 0.5] for q in range(12)]                  # δ(q + 1) = 0.05 + 0.001·q
    starts, calls = [], []
    dv = [0.01 * (q + 1) for q in range(NR.NP_C7 - 1)] + [0.9]                   # the last period the largest

    def seg(m, ke, mt, s_, opts, shift, n, al, c, x0=None):
        if n < NR.NP_C7:                                                       # the start: its state depends on shift and n
            starts.append((shift, n, np.array(al, float), np.array(c, float), x0))
            return [0.0] * n, (7.0 + shift + 100.0 * n, np.asarray(al, float) + 1.0, np.asarray(c, float) + 2.0)
        calls.append((shift, n, np.array(al, float), np.array(c, float), x0))
        return list(dv), None
    monkeypatch.setattr(NR, 'seg_run', seg)
    p = 6
    o = NR.c8f_job((12.0, 20.0, 'paradigm', 0.06, {'seed': 1}, [(p, 2)], (W, None), {p: 0}))
    assert [x['k'] for x in o['rows']] == [1, 2, 3] and len(starts) == 3 and len(calls) == 96
    for (shift, n, al, c, x0), k in zip(starts, (1, 2, 3)):
        assert shift == p - k and n == k and x0 is None and np.array_equal(c, c0)
        assert np.array_equal(al, NR.raised(L, al0, idx, 0, abs(float(W[p - k - 1][0])) - w0))   # δ(p − k), not δ(p)
    alc = NR.raised(L, al0, idx, 0, abs(float(W[p - 1][0])) - w0)
    un = [i for i in range(c0.size) if i not in set(int(i_) for i_ in idx)]
    oth = [i for i in range(al0.size) if i != idx[0]]
    for x, k in zip(o['rows'], (1, 2, 3)):
        als = NR.raised(L, al0, idx, 0, abs(float(W[p - k - 1][0])) - w0) + 1.0   # the start's state at pT (the stub's)
        assert x['j'] == 0 and np.isclose(x['dp'], 0.055) and x['dp'] == abs(float(W[p - 1][0])) - w0
        assert all(cl[1:] == [NR.NP_C7, 0.9] for cl in x['cells'])               # the leave and the largest distance in the last period
        assert x['dev'] == [dv] * 32
        for (shift, n, al, c, x0), (b, *_) in zip(calls[32 * (k - 1): 32 * k], x['cells']):
            A, B, C, D, E_ = (ch == '1' for ch in b)
            assert shift == p and n == NR.NP_C7 and (x0 == 7.0 + (p - k) + 100.0 * k if D else x0 is None)
            assert al[idx[0]] == (als if A else alc)[idx[0]] and np.array_equal(al[oth], (als if E_ else alc)[oth])
            assert c[idx[0]] == c0[idx[0]] + 2.0 * B                             # the named centre, kernel 1's
            assert c[idx[1]] == c0[idx[1]] + 2.0 * C and c[idx[2]] == c0[idx[2]] + 2.0 * C   # C: both other formed centres
            assert np.array_equal(c[un], c0[un] + 2.0 * E_)


def test_round138_c8f_lines_pass_and_alarm_boundaries():
    """Referee round 138 (D3eg-3): a corner on the start's side leaving at R + 2 still passes, at R + 1 it does not; a corner on the cold
    side leaving at R ± 1 is at the alarm, at R ± 2 it is not."""
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1))
    one = "the named formed centre's alone at 1 start at 1 window"
    none_ = "no single group's at 1 start at 1 window"
    cells = lambda row: {b: i for i, (b, *_) in enumerate(row['cells'])}
    for leave, decides in ((6, True), (5, False)):                             # R = 4, the start's side ('1' in bit 1)
        x = _c8f_row(10, 4, 1, 0, None)
        x['cells'][cells(x)['01000']][1] = leave
        s = '\n'.join(NR.c8f_lines([dict(r, rows=[x])]))
        assert (one in s) is decides and (none_ in s) is not decides
    for leave, decides in ((3, True), (5, True), (2, False), (6, False)):      # R = 4, the cold side ('0' in bit 1)
        x = _c8f_row(10, 4, 1, 0, None)
        x['cells'][cells(x)['10111']][1] = leave
        s = '\n'.join(NR.c8f_lines([dict(r, rows=[x])]))
        assert (one in s) is decides and (none_ in s) is not decides


def _pre_row(p, R, f, f0=None):
    """A C8F start decided by the named centre (`_c8f_row`) with per-corner distances 'dev': f(bits) in the period before the alarm's
    (R − 1), f0(bits) (default 0.01) in the others."""
    x = _c8f_row(p, R, 1, 0, None)
    x['dev'] = [[(f(b) if q == R - 2 else (f0(b) if f0 else 0.01)) for q in range(NR.NP_C7)] for b, *_ in x['cells']]
    return x


def test_round138_c8f_pre_alarm_terms():
    """Referee round 138 (D3eg-1): the Möbius terms of kernel 1's distance in the period before the alarm's — windows whose other terms
    are within C8F_TOL pooled, the others' deciding term, largest other main term and largest interaction (named above C8F_TOL), the
    cold restart's farthest pre-alarm period where it is another period —; the largest distance's bimodality; nothing without 'dev'."""
    r = dict(m=12.0, ke=20.0, metric='paradigm', s=0.06, opts=dict(seed=1))
    f = lambda b: 0.2 - 0.1 * (b[1] == '1') - 0.05 * (b[2] == '1') + 0.03 * (b[1] == '1') * (b[2] == '1')   # B, C and B·C
    t, fm = NR.c8f_terms(_pre_row(10, 3, f), 1)
    assert np.isclose(t['01000'], 0.1 / 0.12) and np.isclose(t['00100'], 0.05 / 0.12) and np.isclose(t['01100'], -0.25) and fm == 0.2
    assert all(abs(v) < 1e-12 for S, v in t.items() if S not in ('01000', '00100', '01100')) and np.isclose(sum(t.values()), 1.0)
    quiet = lambda b: 0.2 - 0.1 * (b[1] == '1') + 1e-4 * (b[3] == '1')
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_pre_row(10, 3, f), _pre_row(20, 3, quiet)])]))
    assert ("; the dynamics before the decision do interact: on a response that does not saturate, kernel 1's period-averaged distance "
            "in the period before the alarm's (at most 0.2001 at every corner, below the cut 0.3), the deciding group's term is +1.001 of "
            "the change and every other term at most 1.0e-03 of it at 1 window (1 start), but at m = 12, k_e = 20, s = 0.06, draw 1 from 10 "
            "(period 2) the named centre's +0.833, the other two formed centres' +0.417 and the largest interaction, 01100 (the named centre "
            "with the other two formed centres), -0.250.") in s
    assert ("— bimodal, the leaving corners' 0.800 and the others' 0.300, so its terms restate the outcome — the deciding group's main term "
            "is +1.000") in s
    f2 = lambda b: 0.2 - 0.1 * (b[1] == '1') - 0.05 * (b[4] == '1') + 0.0012 * (b[1] == '1') * (b[4] == '1')   # an interaction below C8F_TOL
    far = lambda b: 0.25 if b == '00000' else 0.24 - 0.04 * (b[1] == '1')                                      # the cold restart farthest at period 1
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_pre_row(10, 3, f2, far)])]))
    assert ("(at most 0.2000 at every corner, below the cut 0.3), at m = 12, k_e = 20, s = 0.06, draw 1 from 10 (period 2) the named "
            "centre's +0.672, the other α and unformed centres' +0.336 and every interaction at most 8.1e-03, and in the period before the "
            "alarm's where the cold restart is farthest from the feature (period 1, 0.2500) the named centre's +1.000, the other α and "
            "unformed centres' +0.200 and the largest interaction, 00011 (the plant's x with the other α and unformed centres), -0.200.") in s
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_pre_row(10, 3, f, lambda b: 0.0)])]))
    assert "farthest" not in s                                              # the period before the alarm's is the farthest
    s = '\n'.join(NR.c8f_lines([dict(r, rows=[_c8f_row(10, 3, 1, 0, None)])]))
    assert "before the decision" not in s and "the deciding group's main term" in s
    bad = _pre_row(10, 3, f)
    bad['cells'][3][1] = 6                                                  # no group decides: no pre-alarm terms
    assert "before the decision" not in '\n'.join(NR.c8f_lines([dict(r, rows=[bad])]))
    assert NR.c8f_pre_text([(r, _pre_row(10, 3, f), 1)], lambda r_, x: 'w').startswith("; the dynamics before the decision")
    assert NR.c8f_pre_text([(r, _pre_row(10, 3, f), None)], lambda r_, x: 'w') == ''


def test_round138_c8f_dev_on_the_raw():
    """Referee round 138 on the raw (D3eg-1): every C8F corner's stored distances give its leave period and largest distance; the Möbius
    terms of one start agree with an independent subset-sum inversion; the md prints the pre-alarm terms and the bimodality."""
    raw = NR.E.load(os.path.join(NR.OUT, "noise_resonance_raw.npz"))
    md = NR.report(**{k: raw.get(k) for k in NR.KEYS})
    n = 0
    for r in raw['C8F']:
        for x in r['rows']:
            assert len(x['dev']) == 32 and all(len(d) == NR.NP_C7 for d in x['dev'])
            for (b, lv, mx), d in zip(x['cells'], x['dev']):
                assert lv == NR.left_periods(np.asarray(d, float)[:, None], NR.LEFT_C)[0] and mx == max(d)
                assert d[x['reduced'] - 2] < NR.LEFT_C                      # the response before the alarm's does not saturate
                n += 1
    assert n == 36 * 32
    x = next(x for r in raw['C8F'] if r['opts'].get('seed') == 2 and r['s'] == 0.06 for x in r['rows'] if x['k'] == 1)
    v = np.array([d[0] for d in x['dev']])                                  # corner i ↔ bits of i (A the highest)
    g = v.copy()
    for q in range(5):                                                      # the inverse zeta transform over the subset lattice
        for i in range(32):
            if i >> q & 1:
                g[i] -= g[i ^ (1 << q)]
    t, _ = NR.c8f_terms(x, 0)
    assert all(np.isclose(t[format(i, '05b')], g[i] / (v[31] - v[0]), rtol=1e-9, atol=1e-12) for i in range(1, 32))
    assert np.isclose(t['01100'], 0.2826, atol=5e-4) and np.isclose(t['01000'], 1.7485, atol=5e-4)
    assert ("— bimodal, the leaving corners' 0.310 … 1.903 and the others' 0.030 … 0.286, so its terms restate the outcome — the deciding "
            "group's main term is +1.000 … +1.001") in md
    assert ("kernel 1's period-averaged distance in the period before the alarm's (at most 0.2581 at every corner, below the cut 0.3), the "
            "deciding group's term is +0.999 … +1.000 of the change and every other term at most 4.0e-03 of it at 10 windows (30 starts), but "
            "at m = 12, k_e = 10, s = 0.02, draw 7 from 107 (period 9) kernel 1's coefficient's +1.018 … +1.112") in md
    assert ("(period 3, 0.0305) kernel 1's coefficient's +1.198 … +1.238, the named centre's -0.223 … -0.198 and every interaction at most "
            "1.3e-03; at m = 12, k_e = 20, s = 0.06, draw 2 from 58 (period 1) the named centre's +1.749 … +1.796, the other two formed "
            "centres' -1.088 … -1.030 and the largest interaction, 01100 (the named centre with the other two formed centres), +0.283 … "
            "+0.292;") in md
