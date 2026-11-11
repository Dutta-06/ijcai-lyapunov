"""Proposition D‴.57: the transit's dead zone — the energy dead zone's frozen block in the sliding limit and its law given the stop, the
gated loop, the report (`experiments/exp_dead_zone.py`), and the dead zone in the closed loop's C integrator
(`closed_loop.simulate(gate=...)`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

from src.core import closed_loop as CL
from src.core import mirror_free as MF
from src.core import rich_transit as RT
import exp_dead_zone as DZ
import exp_mirror_transit as X
import exp_noisy_certificate as NC
import exp_transit_mass as TM


def test_gated_relax_stops_where_the_residuals_rms_reaches_the_threshold():
    """The noiseless path on target 0: the (K − 1)-th block relaxed and the K-th entrant at its seed, gated at half the (K − 1)-th
    relaxed rms: the block's fast flow crosses the threshold and stops on it."""
    cs, As = NC.targets()[0]
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    nK = cs.size
    _, c_prev, w_prev, rms_prev = R['hist'][nK - 1]
    j = [e for e in R['events'] if e[2] == 'enter'][nK - 1][1]
    act = np.flatnonzero(w_prev != 0)
    sg = np.sign(K.phi(c_prev)[:, j] @ (f - K.phi(c_prev) @ w_prev))
    thr = 0.5 * rms_prev
    crossed, cc, wc = DZ.gated_relax(L, K, f, c_prev[np.r_[act, j]], np.r_[w_prev[act], sg * 1e-7], thr)
    assert crossed
    assert abs(np.sqrt(np.mean((f - K.phi(cc) @ wc) ** 2)) / thr - 1) < 1e-8
    assert np.abs(wc[-1]) > 1e-3                                         # the entrant has grown


def test_energy_dead_zone_freezes_the_last_block_on_the_energy_identity():
    """Target 0, s = 0.01 (Part A's first draw): the dead zone stops at the K-th entry; the frozen block's unrelaxed part u obeys |u|² =
    (1 + δ)²s²N − |r_K|² to 10⁻³ (first order) and is the last entrant's amplitude deficit along its kernel off the block's tangent,
    u = δw·P⊥φ_last (to 10⁻³), a shrinkage (δw > 0); the frozen fit is the oracle's displaced by u (an identity of the bookkeeping)."""
    cs, As = NC.targets()[0]
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = cs.size, K.xs.size
    s = 0.01
    xi = s * np.random.default_rng(5700).standard_normal(nx)
    d = DZ.a_draw(L, K, c0, f, sorted(cert['e']), xi, s, nK)
    assert d['exact'] and d['k'] == nK and d['crossed']
    assert abs(d['u2'] / d['u2_pred'] - 1) < 1e-3
    assert d['dw'] > 0 and abs(d['dw'] * d['pnorm'] / np.sqrt(d['u2']) - 1) < 1e-3 and d['cos_last'] > 1 - 1e-3
    assert d['err_dz'] > d['err_or']
    u = np.sqrt(d['u2'] / nx)
    assert abs(d['err_dz'] ** 2 - (d['err_or'] ** 2 + u ** 2 - 2 * d['cross'] * d['err_or'] * u)) < 1e-12 * d['err_dz'] ** 2 + 1e-18


def test_gated_loops_form_the_features_kernels_and_follow_the_sliding_limit():
    """Part B's draw (target 0, s = 0.03, draw 1) at ε = 10⁻⁸ to τ_K + 3W: ungated the loop has formed noise kernels; either gate forms
    the three features' kernels only, the energy dead zone's fit is the limit's frozen fit and the formation dead zone's the oracle's —
    the fits themselves, not only their errors (measured 1.1·10⁻¹¹ rms)."""
    cs, As = NC.targets()[0]
    out = DZ.b_job((0, cs, As, 0.03, 1))
    assert out['exact'] and out['k_lim'] == cs.size
    assert out['none']['formed'] > cs.size and out['energy']['formed'] == out['formation']['formed'] == cs.size
    assert out['energy']['d_dz'] < 1e-8 and out['formation']['d_or'] < 1e-8
    assert out['energy']['err'] > 1.5 * out['formation']['err']


def test_part_B_draws_its_noise_afresh_at_each_level():
    """`b_noise`: Part B's draws at the three levels are not one vector rescaled (referee round 50), and a draw is reproducible."""
    xs = [DZ.b_noise(0, 0, s, 61) / s for s in DZ.LEVELS_B]
    assert all(np.abs(xs[i] - xs[j]).max() > 0.1 for i in range(3) for j in range(i))
    assert np.array_equal(DZ.b_noise(3, 2, 0.03, 61), DZ.b_noise(3, 2, 0.03, 61))


def test_the_stop_law_against_quadrature():
    """`stop_law`: m = E[χ²_{N−2K} | χ²_{N−2K} ≤ (1 + δ)²N] against quadrature of x·p(x) (53.83 at N = 61, K = 3), E|u|²/(s²N) =
    (c − m)/N and E err_dz²/s² = (2K + c − m)/N; the excess's growth with N: the rms ratio over (2 + N(2δ + δ²)/(2K))^{1/2} → 1."""
    from scipy import integrate
    from scipy.stats import chi2
    for N in (61, 181, 541):
        lw = DZ.stop_law(N, 3)
        k, c = N - 6, (1 + DZ.DELTA) ** 2 * N
        num = integrate.quad(lambda x: x * chi2.pdf(x, k), 0, c, limit=400, epsabs=0, epsrel=1e-12)[0]
        assert abs(lw['m'] / (num / chi2.cdf(c, k)) - 1) < 1e-9 and abs(lw['p'] - chi2.cdf(c, k)) < 1e-15
        assert abs(lw['dz'] - (6 + c - lw['m']) / N) < 1e-14 and abs(lw['u'] - (c - lw['m']) / N) < 1e-14
    assert abs(DZ.stop_law(61, 3)['m'] - 53.83) < 5e-3 and abs(DZ.stop_law(61, 3)['dz'] - 0.4260) < 5e-5
    lw = DZ.stop_law(10 ** 5, 3)
    d = 2 * DZ.DELTA + DZ.DELTA ** 2
    assert abs(np.sqrt(lw['dz'] / (6 / 10 ** 5)) / np.sqrt(2 + 10 ** 5 * d / 6) - 1) < 1e-3


def test_gate_stats_and_kernel_events():
    """`gate_stats` on a hand-made filtered energy: shut where z ≤ lo², fully open where z ≥ hi², the first closing at the first shut
    output (not the first below hi²), the last reopening after the last shut stretch, in τ = t/Λ; `kernel_events`: each formed kernel's
    formation, the feature nearest its centre then, and the first output time beyond the range after it."""
    ts = np.arange(10.0)
    z = np.array([9.0, 3.0, 1.0, 0.5, 2.5, 4.0, 0.9, 0.8, 1.5, 0.5])
    g = DZ.gate_stats(z, ts, 1.0, 2.0, 2.0)
    assert g['shut'] == 0.5 and g['open'] == 0.2 and g['first_close'] == 1.0 and g['last_reopen'] == 4.0 and g['shut_end']
    g = DZ.gate_stats(np.full(4, 5.0), ts[:4], 1.0, 2.0, 2.0)
    assert g['shut'] == 0.0 and g['open'] == 1.0 and np.isnan(g['first_close']) and np.isnan(g['last_reopen']) and not g['shut_end']
    C = np.zeros((10, 3))
    C[:, 1] = [-0.4] * 6 + [-1.3] * 4
    C[:, 2] = [0.6] * 8 + [-0.45] * 2                                   # crosses to the other feature after its formation
    ev = DZ.kernel_events(np.array([np.nan, 2.5, 6.0]), ts, C, 2.0, [-0.5, 0.6], 1.2)
    assert [(e['i'], e['formed'], e['feature']) for e in ev] == [(1, 1.25, 0), (2, 3.0, 1)]
    assert ev[0]['ejected'] == 3.0 and np.isnan(ev[1]['ejected'])


def _draw(kept, k, s):
    """A synthetic Part A draw: exact entry set, the order kept or swapped, the stop at entry k (−1: none of the first K)."""
    d = dict(exact=True, exact_order=kept, exits=0, k=k, last=5)
    if k == 3:
        d.update(crossed=True, err_dz=2 * s, err_or=s, u2=0.3 * s * s * 61, u2_pred=0.3 * s * s * 61, cross=0.0, dw=1.8 * s, wj=0.6,
                 pnorm=np.sqrt(0.3 * 61) / 1.8, cos_last=1.0, cos=0.99)
    return d


def test_the_report_on_synthetic_results():
    """A: the draws split by entry order with the stop's count (later includes the draws that stop at none of the first K), the law's
    numbers (E[χ²_{N−2K} | ≤ c] = 53.83); A3: a row per N; B: the fits' distances; C: the gate's shares and times, the timelines."""
    lv = lambda s: dict(s=s, draws=[_draw(True, 3, s), _draw(True, 3, s), _draw(True, -1, s), _draw(False, 3, s),
                                    dict(exact=False, exact_order=False, exits=1)])
    md = DZ.report(A=[dict(q=0, levels=[lv(s) for s in DZ.LEVELS_A])])
    assert '| 0.01 | noiseless | 3 | 2, 1, 0 | 2.000 [2.000, 2.000] |' in md and '| 0.01 | swapped | 1 | 1, 0, 0 |' in md
    assert '= 53.83s² against 55s²' in md
    assert '| 0.01 | noiseless | 2 | 1.0000 [1.0000, 1.0000]; 0.0e+00 |' in md
    A2 = [dict(q=0, N=61, exact0=True, draws=[_draw(True, 3, DZ.S_A2), _draw(True, -1, DZ.S_A2), _draw(False, 3, DZ.S_A2)]),
          dict(q=1, N=61, exact0=False, draws=[])]
    assert '| 61 | 1/2 | 1, 1, 0; 1 |' in DZ.report(A2=A2)
    A2 = [dict(q=0, N=181, exact0=True, draws=[_draw(True, 3, DZ.S_A2)])]
    assert '| 4.0000 (0.2802) | 30.167 | 2.000 (2.907) |' in DZ.report(A2=A2)

    def bm(formed, err, dd):
        return dict(formed=formed, err=err, tau_gate=1.5, rms=0.01, d_dz=dd, d_or=dd)
    B = [dict(q=0, s=0.01, k=0, exact=True, tau_K=1.0, W=1.0, k_lim=3, err_dz_lim=0.02, err_or_lim=0.01,
              none=bm(5, 0.03, 1.0), energy=bm(3, 0.02, 1.3e-8), formation=bm(3, 0.01, 2e-9))]
    assert '| 0.01 | 1, 0 | 5 [5, 5]; 3 [3, 3]; 3 [3, 3] | 3.000, 2.000, 1.000 | 0.667 | 1.3e-08, 2.0e-09 | 1.5000 | — |' in DZ.report(B=B)

    def crun(regx, s, mode, name='base', gate=None, events=()):
        return dict(s=s, regx=regx, mode=mode, name=name, n_act=3, outside=0, w=[0.9, -0.6, 0.75], c=[-0.52, 0.07, 0.61], e_last=0.05,
                    e_true=0.01, fit_win=(0.02, 0.015, 0.07), formed_win=3.0, gate=gate, events=list(events))
    g = dict(shut=0.11, open=0.72, first_close=149.5, last_reopen=float('nan'), shut_end=True)
    ev = [dict(i=4, formed=136.5, feature=0, ejected=164.5), dict(i=5, formed=169.0, feature=2, ejected=float('nan'))]
    md = DZ.report(C=[crun(False, 0.0, 'none'), crun(True, 0.0, 'none'), crun(True, DZ.S_EXTRA_C, 'formation', gate=g, events=ev),
                      crun(True, DZ.S_EXTRA_C, 'none', name='true error')])
    assert ('| x_m | 0.02 | the formation dead zone | 2, 1 | 3, 0 | 0.0200 [0.0150, 0.0700] | 3 | 0.0100 | 0.0500 | 0.11, 0.72 | '
            '149.5, — |') in md
    assert '| base, the formation dead zone | 4: 136.5 (-0.52) → 164.5; 5: 169.0 (+0.61) |' in md
    assert '| true error | ungated | 0, 0 | 3, 0 |' in md and '| true error, ungated | — |' in md


LAW = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
SIG, KE, OM = 0.12, 10.0, 0.2
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)
TS = np.linspace(0, 1.2e4, 41)
OPEN_TOL = 1e-8                                                          # measured: 2.5·10⁻¹¹ on the centres


def test_closed_dead_zone_freezes_the_adaptation_exactly():
    """A gate that never opens (the filtered energy starts at 0 below the thresholds and the filter is slow): every kernel's drive is
    zero — the coefficients stay at 0 and the centres at their start, exactly, under either gate."""
    for mode in ('energy', 'formation'):
        r = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS,
                        gate=dict(mode=mode, Tf=1e12, lo=10.0, hi=20.0, z0=0.0, theta=0.01))
        assert np.abs(r['w']).max() == 0.0 and np.abs(r['c'] - C0).max() == 0.0
        assert 0.0 < r['z'][-1] < 1e-6


def test_open_dead_zone_is_the_ungated_loop():
    """A gate that never closes (thresholds far below the error's energy): the ungated loop, up to the integrator's error control on the
    extra state."""
    r0 = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS)
    assert np.abs(r0['w'][-1]).max() > 0.05                               # kernels have formed
    for mode in ('energy', 'formation'):
        r = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS,
                        gate=dict(mode=mode, Tf=1e12, lo=0.5e-9, hi=1e-9, theta=0.01))
        assert np.abs(r['w'] - r0['w']).max() < OPEN_TOL and np.abs(r['c'] - r0['c']).max() < OPEN_TOL


def test_the_formation_gate_drives_only_the_formed_kernels():
    """The gate shut throughout (the filtered energy starts at 0, slow filter, high thresholds) with kernel 10 formed at the start
    (α = 2, |w| > θ): the formation dead zone drives it (its α and centre move) and keeps every unformed kernel exactly at rest; the
    energy dead zone keeps every kernel at rest."""
    al0 = np.zeros(C0.size)
    al0[10] = 2.0
    ts = np.linspace(0, 200.0, 11)
    rest = np.arange(C0.size) != 10
    for mode in ('energy', 'formation'):
        r = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, al0=al0,
                        gate=dict(mode=mode, Tf=1e12, lo=10.0, hi=20.0, z0=0.0, theta=0.01))
        assert r['w'][0, 10] > 0.01
        assert np.abs(r['al'][:, rest]).max() == 0.0 and np.abs(r['c'][:, rest] - C0[rest]).max() == 0.0
        moved = max(abs(r['al'][-1, 10] - 2.0), abs(r['c'][-1, 10] - C0[10]))
        assert (moved > 0.05) if mode == 'formation' else (moved == 0.0)


def test_noise_of_the_default_and_the_beat_tests():
    """C's noise: by default D‴.49's 16 tones (`closed_loop.tones`, seed 7); with `freqs` the given tones at equal amplitudes, of total rms s;
    none at s = 0."""
    a, b = DZ.noise_of(0.02, {}), CL.tones(0.02, 10.0, 30.0, seed=7)
    assert all(np.array_equal(a[k], b[k]) for k in ('w', 'ph', 'amp'))
    nz = DZ.noise_of(0.02, dict(freqs=(20.0, 20.07)))
    assert np.array_equal(nz['w'], [20.0, 20.07]) and abs(np.sqrt(np.sum(nz['amp'] ** 2) / 2) - 0.02) < 1e-15
    t = np.linspace(0, 2000, 400001)
    x = sum(A * np.sin(w * t + p) for A, w, p in zip(nz['amp'], nz['w'], nz['ph']))
    assert abs(np.sqrt(np.mean(x ** 2)) / 0.02 - 1) < 0.02 and DZ.noise_of(0.0, dict(freqs=(20.0,))) is None


def test_part_D_from_the_fitted_state():
    """D: the fitted state is the noise-free closed loop's rest (each feature's kernel stays at its centre with its amplitude, the others
    at zero); the report's thresholds — the smallest level at which a kernel left — and its coefficient-growth ratio."""
    r = DZ.d_job(('paradigm', 1.0, 0.0))
    assert r['left'] == [-1, -1, -1] and max(r['dev_max']) < 1e-8 and r['others'] < 1e-12
    assert np.allclose(r['w_end'], [0.9, -0.6, 0.75], atol=1e-8)
    D = [dict(metric='paradigm', lam=1.0, s=s, A=[0.9, -0.6, 0.75], left=[-1 if s < 0.02 else 120, -1, -1], dev_max=[0.01] * 3,
              w_end=[a * (1 + 300 * (a * s) ** 2) for a in (0.9, -0.6, 0.75)], others=1e-10) for s in (0.01, 0.02)]
    D.append(dict(metric='wfr', lam=1.0, s=0.04, A=[0.9, -0.6, 0.75], left=[-1] * 3, dev_max=[1e-3] * 3, w_end=[0.9, -0.6, 0.75], others=0.0))
    md = DZ.report(D=D)
    assert 'paradigm, |A| = 0.9: 0.02 (s·|A| = 0.0180); paradigm, |A| = 0.6: kept to s = 0.02' in md
    assert 'wfr, |A| = 0.9: kept to s = 0.04' in md and 'median 300 [5–95% 300, 300] (5 kernels and levels)' in md
    assert '| paradigm | 1 | 0.02 | 120, +0.097 |' in md
