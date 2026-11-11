"""Proposition D‴.58: the experiment's analysis (`experiments/exp_width_mismatch.py`)."""
import os
import sys
import warnings
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_width_mismatch as WM
from src.core import wide_feature as W


def test_flank_peaks_of_the_one_kernel_optimum():
    """B's `flank`: after the optimal one-kernel fit on the samples (centre at the feature, least-squares amplitude) the residual's
    correlation C is nonnegative to rounding and peaks symmetrically ≈ 2σ beyond the node, on both flanks."""
    ss = WM.sigma_star(1.0)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    phi = np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.SIG ** 2))
    w = float(phi @ f / (phi @ phi))
    mn, left, right = WM.flank(1.0, np.array([WM.C_B]), np.array([w]))
    assert abs(mn) < 1e-12 and abs(left - right) < 0.02 and 2.0 < left < 2.5


def test_lead_coefficients_and_the_tail_job():
    """D: the ratios' leading coefficients −C_{j−1}/2^{j−1} (j ≤ 3) with a vanishing next divided difference; `d_job` is the library's
    tail ratio."""
    lc = WM.lead_coeffs(3)
    assert [Fraction(*r['lead']) for r in lc] == [Fraction(-1), Fraction(-1, 2), Fraction(-1, 2)] and all(r['next_dd'] == 0 for r in lc)
    r = WM.d_job((6, 0.5))
    assert r['n_terms'] == 58 and abs(r['ratio'] - W.tail_ratio(6, 0.5, n_terms=58)) < 1e-15


def test_the_discrepancy_stop_prediction():
    """C's prediction: the first K whose rule error e_K satisfies e_K²rf² ≤ s²((1 + δ)² − 1 + 2K/N); larger noise stops earlier."""
    rf, N = 0.4, WM.X_B.size
    ks = [WM.predict_stop(s, rf, N)[0] for s in (1e-4, 1e-3, 1e-2)]
    assert ks[0] > ks[1] > ks[2] >= 1
    K, err = WM.predict_stop(1e-3, rf, N)
    ss = WM.sigma_star(WM.TAU_C)
    e = W.l2_error(*W.half_squeeze(K, WM.SIG, ss), WM.SIG, ss)
    e_prev = W.l2_error(*W.half_squeeze(K - 1, WM.SIG, ss), WM.SIG, ss)
    bound = lambda k: 1e-6 * ((1 + WM.DELTA) ** 2 - 1 + 2 * k / N)
    assert (e * rf) ** 2 <= bound(K) and (e_prev * rf) ** 2 > bound(K - 1)
    assert abs(err - np.sqrt(e ** 2 + 2 * K * 1e-6 / (N * rf ** 2))) < 1e-15


def test_sigma_star_and_the_blocks():
    """σ* = σ(1 + (τ/σ)²)^{1/2}; `blocks`: hist[i + 1] is iteration i's block (its active kernels sorted by centre), relaxed[i] its flag
    (False past the flags' end)."""
    assert abs(WM.sigma_star(1.0) - WM.SIG * np.sqrt(2)) < 1e-15 and abs(WM.sigma_star(0.0) - WM.SIG) < 1e-15
    c = np.array([0.5, -0.2, 0.1])
    R = dict(hist=[(0.0, c, np.zeros(3), 1.0), (1.0, c, np.array([0.0, 0.3, 0.0]), 0.5), (2.0, c, np.array([0.2, 0.3, 0.0]), 0.1)],
             relaxed=[True])
    b = WM.blocks(R)
    assert [x['k'] for x in b] == [1, 2] and [x['relaxed'] for x in b] == [True, False]
    assert np.allclose(b[1]['c'], [-0.2, 0.5]) and np.allclose(b[1]['w'], [0.3, 0.2]) and b[1]['rms'] == 0.1 and b[0]['tau'] == 1.0


def test_sample_fit_sorts_and_projects():
    """`sample_fit` from an unsorted start returns sorted centres, the least-squares coefficients and their rms residual."""
    x = WM.X_B
    phi = lambda c: np.exp(-(x[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * WM.SIG ** 2))
    f = phi([-0.3, 0.2]) @ np.array([0.5, -0.4])
    c, a, rms = WM.sample_fit(phi, f, [0.25, -0.35])
    assert np.all(np.diff(c) > 0) and np.allclose(c, [-0.3, 0.2], atol=1e-8) and np.allclose(a, [0.5, -0.4], atol=1e-8) and rms < 1e-10


def test_the_stop_prediction_uses_2k_over_n():
    """A case where the bound with 2K/N stops at K = 3 and one with K/N would not (N = 4: the variance term dominant)."""
    rf, N = 0.25, 4
    ss_ = WM.sigma_star(WM.TAU_C)
    e3 = W.l2_error(*W.half_squeeze(3, WM.SIG, ss_), WM.SIG, ss_)
    s = e3 * rf / np.sqrt((1 + WM.DELTA) ** 2 - 1 + 1.5 * 3 / N)
    K, err = WM.predict_stop(s, rf, N)
    assert K == 3 and abs(err - np.sqrt(e3 ** 2 + 2 * 3 * s ** 2 / (N * rf ** 2))) < 1e-15


def test_the_tail_job_has_enough_terms_and_the_degree_check_can_fail(monkeypatch):
    """`d_job` at t = 0.8 sums 40/ln(1/t) terms (180), its ratio converged; `lead_coeffs`' (j + 1)-th divided difference is nonzero for
    ratios of degree j + 1 (a fake defect table with P_1(K) = K²)."""
    r = WM.d_job((10, 0.8))
    assert r['n_terms'] == 180 and abs(r['ratio'] - W.tail_ratio(10, 0.8, n_terms=400)) < 1e-12
    fake = lambda K, n: [0] * n + [K * K if n == 2 * K + 2 else 1]
    monkeypatch.setattr(WM.W, 'gh_defects', fake)
    lc = WM.lead_coeffs(1)
    assert lc[0]['next_dd'] != 0 and Fraction(*lc[0]['lead']) == 3


def test_the_report_on_synthetic_results():
    """A2's leading coefficients −1, −1, −3, −15; B's optimal blocks within OPT_TOL and the first non-optimal one; C's stop at the first
    rms within (1 + δ)s; D's Catalan verdict, false for a wrong coefficient."""
    rows = []
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        rows.append(dict(r=r_, K=1, t=P['t'], lam=P['lam'], e_hs=0.1, e_series=0.1, e_first=0.099, e_mm=0.2, e_ls=0.09, y_hs=np.array([0.0])))
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) // W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, 7)}
    share = [dict(r=r_, K=K, share=0.01) for r_ in WM.TAUS_A for K in (5, 10, 20, 30)]
    md = WM.report(A=dict(rows=rows, ratios=ratios, share=share))
    lines = md.splitlines()
    i = next(i for i, l in enumerate(lines) if l.startswith('| j | K = 1, 2, 3, 4, 5'))
    assert [l.rsplit('|', 2)[1].strip() for l in lines[i + 2:i + 6]] == ['-1', '-1', '-3', '-15']
    blk = lambda k, e, e_opt, rel: dict(k=k, tau=float(k), c=WM.C_B + np.linspace(-0.1, 0.1, k), w=np.full(k, 0.8 / k), e=e, e_opt=e_opt,
                                        c_opt=None, relaxed=rel, dc=0.0, dc_hs=0.0, e_hs=e)
    B = [dict(r=0.25, var='base', blocks=[blk(1, 1e-2 * (1 + 5e-4), 1e-2, True), blk(2, 1e-3 * (1 + 2e-3), 1e-3, False)],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.0)], relaxed=[True, False], stiff=0, redone=0, failed=0, rejected=0)]
    md = WM.report(B=B)
    assert '| 0.25 | base | 1 | 2, 1.00e-03, 1.002, False | — | 1/2 | 0, 0, 0 | 0 | 0 |' in md
    rfB = float(np.sqrt(np.mean((WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))) ** 2)))
    B[0]['blocks'][1]['polish'] = [dict(status=0, nfev=5000, accepted=True),
                                   dict(status=1, nfev=12, accepted=False, e0=1.0, e1=1.002e-3 * rfB)]  # its end over the optimum (round 65)
    B[0]['entries'] = [dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5, t_wall=2.0, t_min=2.0),
                       dict(it=1, j=11, u0=0.9, rate=0.1, rate_max=0.4, t_wall=1.0, t_min=1.0)]
    md = WM.report(B=B)
    assert '| 2, 1.00e-03, 1.002, False | 0/5000/a 1/12/r (1.002) | 1/2 |' in md        # the first non-optimal block's polishes (round 54)
    assert '| none (the first entry) | 1 (1) | 1 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |' in md   # the entries by regime (round 60)
    assert '| an optimal block | 1 (1) | 0 | 0.25 (0.250–0.25) | 0.900 (0.900) | 1, 0.58σ (0.58–0.58σ) | — |' in md
    assert '| all | 2 (2) | 1 | 0.62 (0.250–1.00) | 0.450 (0.000) | 1, 0.58σ (0.58–0.58σ) | — |' in md and '| a non-optimal block |' not in md
    x = dict(e_hs=1.0, e_best=0.5, e_opt=0.4)
    assert WM.gap_share(x) == '83.3%' and WM.gap_share(dict(x, e_opt=0.6)) == '100.0%, optimum unconverged'
    assert WM.gap_share(dict(x, e_opt=0.5 * (1 + 1e-12))) == '100.0%'                 # within UNCONV: not flagged
    C = [dict(s=1e-3, q=0, k=[1, 2, 3], rms=[5e-3, 1.05e-3, 0.9e-3], err=[0.02, 0.005, 0.004], at=[1, 2, 3], gaps=[0.5, 0.2], exits=0,
              stiff=0, runaway=False, relaxed=[True] * 3, rf=0.25)]
    md = WM.report(C=C)
    assert '| 0.001 | 2 [2, 2] (1/1 stopped) |' in md and '| 1 [1, 1] (the stop in 0) |' in md
    D = dict(rows=[dict(K=10, t=0.5, ratio=1.0, n_terms=58)], lead=[dict(j=1, lead=(-1, 1), next_dd=0.0), dict(j=2, lead=(-1, 2), next_dd=0.0)])
    assert 'at every j' in WM.report(D=D) and 'NOT' not in WM.report(D=D)
    D['lead'].append(dict(j=3, lead=(1, 2), next_dd=0.0))
    assert 'NOT at every j' in WM.report(D=D)


def test_the_quoted_cells_on_synthetic_results():
    """Referee round 54: the cells the text quotes, on synthetic results. A3: the ratios per kernel at the last K against their rates, the
    best variance at K = 1 and at the last K with an optimum, the gap shares and their range (one optimum unconverged); A4: the errors and
    the last ratio. B: the stall, failure and failed-chunk counts (distinct), the first non-optimal block with its polishes, the summary
    over every width and variant, the base path's table (an exit among the events, ratios only after one more kernel) and the flank of
    its optimal blocks. C: the stop at the first rms within (1 + δ)s, the error there over the oracle, over the bias–variance prediction
    from B's base path (not the shifted grid's) and over the bound, the largest gap's kernel count, the runaway count, a level where no
    draw stops."""
    rows, narrow = [], []
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        x0 = dict(r=r_, t=P['t'], lam=P['lam'], y_hs=np.array([0.0]))
        rows += [dict(x0, K=1, e_hs=0.1, e_series=0.1, e_first=0.099, e_mm=0.5, e_ls=0.09, e_fock=0.2, e_mercer=0.3, e_best=0.09,
                      v_best=0.98765, e_opt=0.08),
                 dict(x0, K=2, e_hs=0.04, e_series=0.04, e_first=0.0399, e_mm=0.25, e_ls=0.035, e_fock=0.05, e_mercer=0.1, e_best=0.03,
                      v_best=0.97531, e_opt=0.035 if r_ == 2.0 else 0.025),
                 dict(x0, K=3, e_hs=0.016, e_series=0.016, e_first=0.0159, e_mm=0.125, e_ls=0.015, e_fock=0.0125, e_mercer=0.0333,
                      e_best=0.012, v_best=0.96)]
    for q in WM.NARROW_A:
        narrow += [dict(q=q, K=K, t=(1 - q) / (1 + q), e=e, e_series=e) for K, e in zip((1, 2, 3), (0.1 * q, 0.03 * q, 0.009 * q * q))]
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) // W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, 7)}
    share = [dict(r=r_, K=K, share=0.01) for r_ in WM.TAUS_A for K in (5, 10, 20, 30)]
    md = WM.report(A=dict(rows=rows, ratios=ratios, share=share, narrow=narrow))
    assert 'over all K with an optimum the shares 50.0%–100.0% (1 optima unconverged):' in md
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        b2 = np.sqrt(1 + 2 * r_ ** 2)
        last = '0.87500 (100.0%, optimum unconverged)' if r_ == 2.0 else '0.62500 (66.7%)'
        assert (f"| {r_:g} | 2.500 ({1 / P['t']:.3f}) | 4.000 ({(2 - P['lam']) / P['lam']:.3f}) | 3.003 ({(b2 + 1) / (b2 - 1):.3f}) | "
                f"0.98765, 0.90000, 0.80000 (50.0%) | 2, 0.97531, 0.75000, {last} |") in md
    for q in WM.NARROW_A:
        assert (f"| {q:g} | {(1 + q) / (1 - q):.4f} | 1: {0.1 * q:.3e} (0e+00), 2: {0.03 * q:.3e} (0e+00), 3: {0.009 * q * q:.3e} (0e+00) | "
                f"{0.03 / (0.009 * q):.4f} |") in md

    def blk(k, tau, e, e_opt, rel, dc=0.0, dc_hs=0.0, e_hs=None, shift=0.0, **kw):
        return dict(k=k, tau=tau, c=WM.C_B + shift + 0.2 * WM.SIG * (np.arange(k) - (k - 1) / 2), w=np.full(k, WM.A_B / k), e=e,
                    e_opt=e_opt, c_opt=None, relaxed=rel, dc=dc, dc_hs=dc_hs, e_hs=e if e_hs is None else e_hs, **kw)
    rfB = float(np.sqrt(np.mean((WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))) ** 2)))
    pol = [dict(status=0, nfev=5000, accepted=True), dict(status=1, nfev=12, accepted=False, e0=1.0, e1=3.006e-8 * rfB)]
    path = lambda r_, var, es: dict(r=r_, var=var, blocks=[blk(k, 10.0 + k, e * (1 + 1e-4), e, True) for k, e in es],
                                    events=[(10.0 + k, 9 + k, 'enter', 0.0) for k, e in es], relaxed=[True] * len(es), stiff=0, redone=0,
                                    failed=0, rejected=0)
    B = [dict(r=0.25, var='base', blocks=[blk(1, 10.72, 1.0005e-2, 1e-2, True, dc=1.2e-8, e_hs=1.1e-2, shift=0.03),
                                          blk(2, 11.87, 2.0002e-3, 2e-3, True, dc=2.4e-8, dc_hs=0.012, e_hs=2.5e-3, shift=-0.02),
                                          blk(3, 12.85, 3.006e-8, 3e-8, False, dc=3.6e-6, dc_hs=0.024, e_hs=5e-4, polish=pol),
                                          blk(2, 13.5, 2.1e-3, 2e-3, True, e_hs=2.5e-3)],
              events=[(10.72, 10, 'enter', 0.0), (11.87, 11, 'enter', 0.1), (12.85, 12, 'enter', 0.2), (13.5, 11, 'exit', 0.1)],
              relaxed=[True, True, False, True], stiff=1, redone=2, failed=3, rejected=4),
         dict(r=0.25, var='grid +0.05', blocks=[blk(1, 10.7, 1.0003e-2, 1e-2, True), blk(2, 11.9, 6e-7, 5e-7, True)],
              events=[(10.7, 10, 'enter', 0.05), (11.9, 11, 'enter', 0.15)], relaxed=[True, True], stiff=0, redone=0, failed=0, rejected=0),
         path(1.0, 'base', ((1, 0.02), (2, 0.005), (3, 0.004))), path(1.0, 'grid +0.05', ((1, 0.03), (2, 0.007), (3, 0.006)))]
    rf, N, s = 0.25, WM.X_B.size, 1e-3
    C = [dict(s=s, k=[1, 2, 3, 4], rms=[5e-3, 1.15e-3, 1.08e-3, 0.95e-3], err=[0.02, 0.006, 0.0042, 0.004], at=[1, 2, 3, 3],
              gaps=[0.5, 0.9, 0.2], runaway=False, rf=rf),
         dict(s=s, k=[1, 2, 3, 4], rms=[4e-3, 1.09e-3, 1.0e-3, 0.9e-3], err=[0.03, 0.0055, 0.005, 0.0049], at=[1, 2, 2, 2],
              gaps=[0.2, 0.9, 0.5], runaway=True, rf=rf),
         dict(s=s, k=[1, 2, 3, 4], rms=[5e-3, 3e-3, 2e-3, 1.2e-3], err=[0.02, 0.01, 0.008, 0.007], at=[1, 1, 1, 1], gaps=[0.4, 0.3, 0.2],
              runaway=False, rf=rf),
         dict(s=1e-2, k=[1, 2], rms=[0.05, 0.02], err=[0.2, 0.1], at=[1, 1], gaps=[0.5], runaway=False, rf=rf)]
    md = WM.report(B=B, C=C)
    assert '| 0.25 | base | 1, 2 | 3, 3.00e-08, 1.002, False | 0/5000/a 1/12/r (1.002) | 3/4 | 1, 2, 3 | 4 | 1 |' in md
    assert '| 0.25 | grid +0.05 | 1 | 2, 5.00e-07, 1.2, True | — | 2/2 | 0, 0, 0 | 0 | 0 |' in md
    assert '| 1 | base | 1, 2, 3 | none | — | 3/3 | 0, 0, 0 | 0 | 0 |' in md
    mg = sorted(WM.cs_margin(b['r'], b['var'], o) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= WM.OPT_TOL and o['relaxed'])
    assert ("Over every width and variant: the optimal blocks' optimum errors down to 2.00e-03; the non-optimal blocks' from 2.00e-03 down "
            f"(1 of 3 failed the stationarity test; for the other 2 its threshold over the Cauchy–Schwarz bound — the residual's cosine with "
            f"the block's directions it allows — is {mg[0]:.1e}–{mg[1]:.2f} (median {np.median(mg):.1e}; above 0.1 in "
            f"{sum(m > 0.1 for m in mg)}): a weak test); non-optimal blocks with an optimum above 10⁻⁷: 2; the non-optimal blocks' errors "
            "3.01e-08 … 2.10e-03.") in md
    assert '### τ/σ = 0.25, base (events: 10+, 11+, 12+, 11−)' in md
    for row in ('| 1 | 10.72 | 1.0005e-02 | 1.0000e-02 | 1.000500 | True | 1.0e-07 | 0.000 | 1.1000e-02 | — |',
                f'| 2 | 11.87 | 2.0002e-03 | 2.0000e-03 | 1.000100 | True | 2.0e-07 | 0.100 | 2.5000e-03 | {1.0005e-2 / 2.0002e-3:.3f} |',
                f'| 3 | 12.85 | 3.0060e-08 | 3.0000e-08 | 1.002000 | False | 3.0e-05 | 0.200 | 5.0000e-04 | {2.0002e-3 / 3.006e-8:.3f} |',
                '| 2 | 13.5 | 2.1000e-03 | 2.0000e-03 | 1.050000 | True | 0.0e+00 | 0.000 | 2.5000e-03 | — |'):
        assert row in md

    def flank(r_, c, w):                                   # C(c) on the fine grid, by a loop
        x = WM.X_B
        r = WM.A_B * np.exp(-(x - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
        for ck, wk in zip(c, w):
            r = r - wk * np.exp(-(x - ck) ** 2 / (2 * WM.SIG ** 2))
        cg = np.linspace(x[0], x[-1], 3201)
        Cg = np.array([np.mean(r * np.exp(-(x - g) ** 2 / (2 * WM.SIG ** 2))) for g in cg])
        lft, rgt = cg < min(c), cg > max(c)
        return (Cg.min() / Cg.max(), (min(c) - cg[lft][np.argmax(Cg[lft])]) / WM.SIG,
                (cg[rgt][np.argmax(Cg[rgt])] - max(c)) / WM.SIG)
    fl = [flank(0.25, o['c'], o['w']) for o in B[0]['blocks'][:2]]
    assert all(f"{x[1]:.2f}" != f"{x[2]:.2f}" for x in fl)                      # off-centre blocks: the two flanks differ
    assert ("after the optimal blocks (K = 1, 2): min C/max C within ±" + f"{max(abs(x[0]) for x in fl):.0e}, its peaks beyond the "
            "outermost nodes at (left, right) " + ', '.join(f"({x[1]:.2f}, {x[2]:.2f})σ" for x in fl) + '.') in md
    fmt = lambda v, f_: f"{np.median(v):{f_}} [{np.min(v):{f_}}, {np.max(v):{f_}}]"
    ks, es = [3, 2], [0.0042, 0.0055]
    eopt = {2: 0.005, 3: 0.004}
    Kp, ep = WM.predict_stop(s, rf, N)
    assert (f"| 0.001 | 2.5 [2, 3] (2/3 stopped) | {Kp} | 4.850e-03 [4.200e-03, 5.500e-03] | {ep:.3e} | 2.5 [2, 3] | 2 [2, 2] (the stop "
            "in 1) |") in md
    r1 = [e / (s * np.sqrt(2 * k / N) / rf) for k, e in zip(ks, es)]
    r2 = [e / np.sqrt(eopt[k] ** 2 + 2 * k * s ** 2 / (N * rf ** 2)) for k, e in zip(ks, es)]
    r3 = [e / (s * np.sqrt((1 + WM.DELTA) ** 2 - 1 + 4 * k / N) / rf) for k, e in zip(ks, es)]
    assert f"| 0.001 | {fmt(r1, '.2f')} | {fmt(r2, '.3f')} | {max(r3):.3f} | 1 |" in md
    Kp, ep = WM.predict_stop(1e-2, rf, N)
    assert f"| 0.01 | — (0/1 stopped) | {Kp} | — | {ep:.3e} | — | — |" in md and '| 0.01 | — | — | — | — |' in md


def test_the_stop_law():
    """Referee round 54: `stop_law`'s K(s) is the root of the first-term formula, its ceiling the predicted stop at the three levels; the
    stop's probabilities sum to one but the chance of no stop by K = 12 and are largest at the predicted K; the report prints them with the expected and observed counts."""
    ss = WM.sigma_star(WM.TAU_C)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    rf, N, t = float(np.sqrt(np.mean(f ** 2))), WM.X_B.size, W.params(WM.SIG, ss)['t']
    for s in WM.LEVELS_C:
        Kr, P = WM.stop_law(s, rf, N)
        rhs = np.log((1 - t * t) ** 0.25 * (np.pi * Kr) ** 0.25 * rf / (s * np.sqrt(2 * WM.DELTA + WM.DELTA ** 2))) / np.log(1 / t)
        Kp = WM.predict_stop(s, rf, N)[0]
        assert abs(Kr - rhs) < 1e-12 and int(np.ceil(Kr)) == Kp and max(P, key=P.get) == Kp
        assert 1 - 1e-3 < sum(P.values()) <= 1 + 1e-12 and min(P.values()) >= 0           # the rest: no stop by K = 12
    C = [dict(s=1e-3, k=[3, 4], rms=[2e-3, 1e-3], err=[0.01, 0.005], at=[3, 4], gaps=[0.5], runaway=False, rf=rf) for _ in range(2)]
    Kr, P = WM.stop_law(1e-3, rf, N)
    md = WM.report(C=C)
    assert (f"| 0.001 | {Kr:.2f} | " + ', '.join(f"{K}: {p:.3f}" for K, p in P.items() if p >= 5e-4) + f" | {P[4]:.3f} | "
            f"{2 * P[4]:.1f} of 2 | 2 |") in md


def test_b2_variable_projection_from_the_blocks():
    """B2 (rounds 54 and 60): from the half-squeeze nodes themselves both starts reach one optimum, at the block's own nodes; the summary
    counts the blocks in the optimum's basin only where the optimum is at least E_FLOOR (below it the solve stops short), the misses
    within 1% (in its basin, unconverged) apart from the far ones with their nodes' distance, and each path's first non-optimal block —
    a path being a width and a variant (B2 in the paths' order)."""
    from src.core import rich_transit as RT
    yh, _ = WM.W.half_squeeze(2, WM.SIG, WM.sigma_star(1.0))
    c = WM.C_B + np.sort(yh)
    q = WM.b2_job((1.0, 'base', 2, c))
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    assert abs(q['e_block'] / q['e_rule'] - 1) < 1e-9 and 0 < q['e_rule'] < 0.1
    assert q['dc'] == float(np.abs(c - WM.sample_fit(RT.Bumps(WM.X_B, WM.SIG).phi, f, c)[0]).max()) and 0 < q['dc'] < 0.01 * WM.SIG
    assert 1 < q['cond_vp'] < q['cond_joint'] < 1e12
    B2 = [dict(r=1.0, var='base', k=10, e_block=1e-6, e_rule=1.5e-6, dc=0.0),
          dict(r=1.0, var='base', k=11, e_block=1.005e-6, e_rule=1e-6, dc=0.001),
          dict(r=1.0, var='grid +0.05', k=10, e_block=3e-6, e_rule=1e-6, dc=0.3),
          dict(r=0.5, var='base', k=8, e_block=5e-12, e_rule=1e-12, dc=0.2),
          dict(r=0.5, var='base', k=9, e_block=2e-11, e_rule=2e-11, dc=0.0),
          dict(r=0.5, var='polish 2σ', k=8, e_block=4e-6 * (1 + 1e-4), e_rule=4e-6, dc=0.0)]
    assert WM.b2_line(B2)[0] == (
        "Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in 2 of the 4 blocks whose "
        "optimum is at least 1e-10; 1 ends within 1% of it (1.005 times: in its basin, unconverged); 1 ends 3 times above it — the "
        "blocks' largest node distances from the optimum's 2.50σ —. The 2 blocks with an optimum below 1e-10, where the solve stops short of the "
        "optimum (A), are not counted (1 reached). From the first non-optimal block of each path: 2 of 3 (the others 3 times above it); on "
        "1 path the first non-optimal block's optimum is below 1e-10.")
    assert WM.b2_line(None) == []
    for q, (cj, cv, sm) in zip(B2, ((1e7, 1e4, 1.5), (2e7, 2e4, 9.0), (4e8, 3e5, 2.0), (9e9, 9e9, 7.0), (9e9, 9e9, 7.0), (5e7, 2e5, 1.0))):
        q.update(cond_joint=cj, cond_vp=cv, smin_ratio=sm, vp_status=3, vp_nfev=77, vp_gap=0.002, vp_cos=3e-4)
    md = WM.b2_line(B2)[0]
    assert md.endswith("on 1 path the first non-optimal block's optimum is below 1e-10.")
    for q, (cx, vx, rx, re) in zip(B2, ((1e6, 2e3, 2.5, 80.0), (2e6, 3e3, 9.0, 9.0), (5e8, 4e4, 1.2, 2000.0), (9e9, 9e9, 7.0, 7.0),
                                        (9e9, 9e9, 7.0, 7.0), (3e6, 1e4, 1.6, 49.0))):
        q.update(cond_joint_x=cx, cond_vp_x=vx, smin_ratio_x=rx, smin_ratio_eq=re)
    md = WM.b2_line(B2)[0]
    assert md.endswith(" There the Jacobian's condition number is 1e+06–5e+08 for the joint fit in (w, c) with the centres in x, the "
                       "polish's parameters (1e+07–4e+08 with the centres in units of σ) and 2e+03–4e+04 for variable projection's in c "
                       "alone (analytic; the same in either unit of the centres); the ratio of their smallest singular values, variable "
                       "projection's over the joint fit's, is at least 1 in any common scaling of the centres for Kaufman's form of the "
                       "former, the centre columns projected off the amplitudes' (referee round 77: the joint fit's smallest singular value "
                       "is at most theirs — a Schur complement), from which the analytic Golub–Pereyra form used here differs by a term "
                       "proportional to the residual; its size depends on the units: 1–2 (median 1.5) with the centres in units of σ, "
                       "1.2–2.5 (median 1.6) in x, 49–2e+03 (median 80) with every column scaled to unit norm, which scales the two "
                       "problems' centre columns differently; TRF itself works in Coleman–Li's variables, the bounded amplitudes scaled by "
                       "v^{1/2}.")
    assert ("where variable projection stalls: its two closest nodes 0.002σ apart, its end not stationary (the largest cosine of its "
            "residual with a Jacobian column 3.0e-04), its stop 3/77") in md          # the far miss (grid +0.05)


def test_the_entries_by_the_block_they_join():
    """Referee round 60: `entry_lines` apart by regime — the first entry, entries joining an optimal block and a non-optimal one — on
    asymmetric data: the largest rate exactly (a rate 0.95 of it is not), medians (not means), beyond the outermost node at any distance
    (0.3σ is beyond) and between the nodes with the duals."""
    blk = lambda k, opt, c: dict(k=k, tau=float(k), c=np.array(c), w=np.full(k, 0.1), e=1e-3 * (1 + (5e-4 if opt else 5e-3)), e_opt=1e-3,
                                 c_opt=None, relaxed=True, dc=0.0, dc_hs=0.0, e_hs=1e-3)
    B = [dict(r=1.0, var='base', blocks=[blk(1, True, [0.0]), blk(2, False, [-0.1, 0.1]), blk(3, True, [-0.2, 0.0, 0.2]),
                                         blk(4, True, [-0.3, -0.1, 0.1, 0.3])],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.036), (3.0, 12, 'enter', 0.05), (4.0, 13, 'enter', 0.5),
                      (5.0, 14, 'enter', 0.42)],
              entries=[dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5), dict(it=1, j=11, u0=0.7, rate=0.38, rate_max=0.4),
                       dict(it=2, j=12, u0=0.97, rate=0.1, rate_max=0.4), dict(it=3, j=13, u0=-0.1, rate=-0.2, rate_max=0.2),
                       dict(it=4, j=14, u0=0.9, rate=0.08, rate_max=0.4)])]
    md = '\n'.join(WM.entry_lines(B))
    assert '| none (the first entry) | 1 (1) | 1 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |' in md
    assert '| an optimal block | 3 (3) | 1 | 0.95 (0.200–1.00) | 0.700 (0.100) | 3, 1.00σ (0.30–2.50σ) | — |' in md
    assert '| a non-optimal block | 1 (1) | 0 | 0.25 (0.250–0.25) | 0.970 (0.970) | — | 1 (0.970–0.970) |' in md
    assert '| all | 5 (5) | 2 | 0.95 (0.200–1.00) | 0.700 (0.000) | 3, 1.00σ (0.30–2.50σ) | 1 (0.970–0.970) |' in md
    B2_ = B + [dict(B[0], var='polish 2σ', entries=[dict(e, u0=e['u0'] + 5e-6) for e in B[0]['entries']]),
               dict(B[0], var='grid +0.05')]                       # a repeat (same events) and another grid (other events)
    md2 = '\n'.join(WM.entry_lines(B2_))
    assert '| an optimal block | 6 (9) |' in md2 and 'the duals of one event agree to 5.0e-06' in md2
    B3_ = B + [dict(B[0], var='polish 2σ', blocks=[B[0]['blocks'][0], blk(2, True, [-0.1, 0.1])] + B[0]['blocks'][2:])]
    md3 = '\n'.join(WM.entry_lines(B3_))                              # one event joins an optimal block in one variant only (round 65)
    assert "1 events join an optimal block in some variants and a non-optimal one in others" in md3 and '| all | 5 (10) |' in md3
    assert '| an optimal block | 4 (7) |' in md3 and 'events join an optimal block' not in md + md2
    assert WM.entry_lines([dict(B[0], entries=None)]) == []


def test_the_stop_laws_moments_are_the_projections():
    """Referee round 60: `stop_law`'s probabilities are the normal law of rms²/s² with the mean b_K² + 1 − 2K/N and the sd
    (2(N − 2K)/N² + 4b_K²/N)^{1/2}, b_K = e_K‖f*‖_N/s (recomputed here from the half-squeeze errors); and these are the exact first two
    moments of |R(f + ξ)|²/(Ns²), R the projection off a 2K-dimensional subspace and Rf of rms b_Ks (Monte Carlo, 4·10⁴ draws)."""
    from scipy.stats import norm
    ss = WM.sigma_star(WM.TAU_C)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    rf, N = float(np.sqrt(np.mean(f ** 2))), WM.X_B.size
    for s in WM.LEVELS_C:
        cdf = [0.0]
        for K in range(1, 13):
            b2 = (W.l2_error(*W.half_squeeze(K, WM.SIG, ss), WM.SIG, ss) * rf / s) ** 2
            cdf.append(norm.cdf(((1 + WM.DELTA) ** 2 - b2 - 1 + 2 * K / N) / np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b2 / N)))
        P = WM.stop_law(s, rf, N)[1]
        assert max(abs(P[K] - (cdf[K] - cdf[K - 1])) for K in range(1, 13)) < 1e-12
    rng = np.random.default_rng(0)
    K, b, s = 4, 0.7, 1.0
    Q = np.linalg.qr(rng.standard_normal((N, 2 * K)))[0]
    bias = rng.standard_normal(N)
    bias -= Q @ (Q.T @ bias)
    bias *= b * s * np.sqrt(N) / np.linalg.norm(bias)
    xi = rng.standard_normal((40000, N)) * s
    res = bias[None, :] + xi - (xi @ Q) @ Q.T
    m = np.sum(res ** 2, axis=1) / (N * s ** 2)
    assert abs(m.mean() - (b * b + 1 - 2 * K / N)) < 4 * m.std() / np.sqrt(m.size)
    assert abs(m.std() / np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b * b / N) - 1) < 0.02


def test_the_stationarity_test_holds_by_cauchy_schwarz_below_an_error():
    """Referee round 60: `cs_vacuous` — a block passes the path's absolute stationarity test at any state once its rms error times the
    largest of rms(φ_i) and |w_i|·rms(∂_cφ_i) is at most 10⁻⁹·max_j|⟨f*φ_j⟩| over the starting grid (shifted in 'grid +0.05')."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    o = dict(c=np.array([WM.C_B - 0.05, WM.C_B + 0.05]), w=np.array([0.5, 0.4]))
    lim = max(np.sqrt(np.mean(K.phi(o['c']) ** 2, axis=0)).max(), (o['w'] * np.sqrt(np.mean(K.dphi(o['c']) ** 2, axis=0))).max())
    for r_ in (0.5, 2.0):
        f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
        rf = float(np.sqrt(np.mean(f ** 2)))
        for var, sh in (('base', 0.0), ('grid +0.05', 0.05)):
            e_ = 1e-9 * float(np.abs(K.phi(np.linspace(-1, 1, 21) + sh).T @ f).max()) / WM.X_B.size / (rf * lim)
            assert WM.cs_vacuous(r_, var, dict(o, e=0.999 * e_)) and not WM.cs_vacuous(r_, var, dict(o, e=1.001 * e_))


def test_the_summary_counts_the_vacuous_stationarity_tests():
    """Referee round 60: B's summary counts the relaxed non-optimal blocks whose stationarity test holds at any state."""
    blk = lambda k, e, q: dict(k=k, tau=float(k), c=WM.C_B + 0.1 * np.arange(k), w=np.full(k, 0.8 / k), e=e * q, e_opt=e, c_opt=None,
                               relaxed=True, dc=0.0, dc_hs=0.0, e_hs=e)
    B = [dict(r=1.0, var='base', blocks=[blk(1, 1e-2, 1.0001), blk(2, 1e-14, 1.5), blk(3, 1e-6, 1.5), blk(4, 1e-7, 1.5)], events=[],
              relaxed=[True] * 4, stiff=0, redone=0, failed=0, rejected=0)]
    md = WM.report(B=B)
    assert ("(0 of 3 failed the stationarity test; 1 of those that passed it pass it at any state by Cauchy–Schwarz, the test being "
            "absolute — their errors up to 1.5e-14; for the other 2 its threshold over the Cauchy–Schwarz bound") in md
    mg = sorted(WM.cs_margin(1.0, 'base', o) for o in B[0]['blocks'][2:])
    assert f"is {mg[0]:.1e}–{mg[1]:.2f} (median {np.median(mg):.1e}; above 0.1 in {sum(m > 0.1 for m in mg)}): a weak test" in md


def test_merge_b_replaces_and_orders_the_paths():
    """Round 60's `--variants`: the new paths replace the raw's of the same width and variant, the others are kept, in the order of
    VARIANTS_B, then TAUS_B."""
    v5 = WM.VARIANTS_B[5][0]
    old = np.array([dict(r=1.0, var=v5, x=0), dict(r=0.5, var='base', x=0), dict(r=0.25, var=v5, x=0)], dtype=object)
    m = WM.merge_B(old, [dict(r=2.0, var=v5, x=1), dict(r=1.0, var=v5, x=1)])
    assert [(b['r'], b['var'], b['x']) for b in m] == [(0.5, 'base', 0), (0.25, v5, 0), (1.0, v5, 1), (2.0, v5, 1)]
    assert [(b['r'], b['x']) for b in WM.merge_B(None, [dict(r=2.0, var=v5, x=1)])] == [(2.0, 1)]
    assert v5 == 'polish 2σ, 5·10⁴, repeated, no gradient test' and WM.VARIANTS_B[5][5] is None
    assert all(v[5] == 1e-15 for v in WM.VARIANTS_B[:5])


def test_b_paths_reproduce_the_raw():
    """Referee round 60: `b_job` reproduces the committed raw on B's two fastest paths at τ/σ = 2 — the base variant (rejected first
    polishes) and the repeated one (rejected repeats): the events, the blocks' kernels, errors and optima, the polish logs, the counts
    and the entries."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'width_mismatch_raw.npz'))
    for var in ('base', 'polish 2σ, 5·10⁴, repeated'):
        ref = next(b for b in raw['B'] if b['r'] == 2.0 and b['var'] == var)
        b = WM.b_job((2.0, var))
        assert [(j, k) for t, j, k, c in b['events']] == [(j, k) for t, j, k, c in ref['events']]
        assert np.allclose([e[0] for e in b['events']], [e[0] for e in ref['events']], rtol=1e-9, atol=0)
        assert [o['k'] for o in b['blocks']] == [o['k'] for o in ref['blocks']]
        for key in ('e', 'e_opt'):
            assert np.allclose([o[key] for o in b['blocks']], [o[key] for o in ref['blocks']], rtol=1e-9, atol=0)
        log = lambda x: [[(p['status'], p['nfev'], p['accepted']) for p in o['polish']] for o in x['blocks']]
        assert log(b) == log(ref) and any(p for o in b['blocks'] for p in o['polish'])
        assert [b[k] for k in ('rejected', 'stiff', 'redone', 'failed')] == [ref[k] for k in ('rejected', 'stiff', 'redone', 'failed')]
        assert [e['j'] for e in b['entries']] == [e['j'] for e in ref['entries']]
        assert max(abs(e['u0'] - f_['u0']) + abs(e['rate'] / f_['rate'] - 1) for e, f_ in zip(b['entries'], ref['entries'])) < 1e-9
    assert ref['rejected'] == 0 and any('e0' in p for o in ref['blocks'] for p in o['polish'])
    assert next(b for b in raw['B'] if b['r'] == 2.0 and b['var'] == 'base')['rejected'] == 13


def test_b_job_passes_the_variants_polish_settings(monkeypatch):
    """Referee round 65 (X01): `b_job` hands `sliding_iss` each variant's polish step, budget, repetitions and gradient test (None: the
    sixth variant), and the starting grid's shift."""
    from src.core import rich_transit as RT
    seen = []

    def stub(L, K, f, c0, rho, **kw):
        seen.append((np.asarray(c0).copy(), kw))
        return dict(hist=[(0.0, c0, np.zeros(len(c0)), 1.0)], relaxed=[], polish_log=[], events=[], stiff_redone=0, solver_redone=0,
                    solver_failed=0, polish_rejected=0, entry_info=[])
    monkeypatch.setattr(RT, 'sliding_iss', stub)
    for name, step, shift, nfev, rep_, gt in WM.VARIANTS_B:
        WM.b_job((0.25, name))
        c0, kw = seen[-1]
        assert (kw['polish_step'], kw['polish_nfev'], kw['repolish'], kw['polish_gtol']) == (step, nfev, rep_, gt)
        assert abs(c0[0] - (-1.0 + shift)) < 1e-12
    assert seen[-1][1]['polish_gtol'] is None


def test_the_polish_defaults_are_pinned():
    """Referee round 65 (R01): `polish_block`'s and `sliding_iss`'s gradient test default to 10⁻¹⁵ — other experiments call the polish
    without it."""
    import inspect
    from src.core import rich_transit as RT
    assert inspect.signature(RT.polish_block).parameters['gtol'].default == 1e-15
    assert inspect.signature(RT.sliding_iss).parameters['polish_gtol'].default == 1e-15


def test_b2_takes_dc_from_the_better_start(monkeypatch):
    """Referee round 65 (X02): `b2_job`'s dc is the block's distance from the better of the two variable-projection ends."""
    ends = iter([(np.array([0.0, 0.1]), np.array([1.0, 1.0]), 2e-6, 3, 50), (np.array([0.02, 0.13]), np.array([1.0, 1.0]), 1e-6, 2, 40)])
    monkeypatch.setattr(WM, 'vp_fit', lambda Phi_of, f, y0, max_nfev=20000: next(ends))
    q = WM.b2_job((1.0, 'base', 2, np.array([0.0, 0.1])))
    assert abs(q['dc'] - 0.03) < 1e-12 and q['vp_status'] == 3 and q['vp_nfev'] == 50 and abs(q['vp_gap'] - 0.1 / WM.SIG) < 1e-12


def test_the_vp_jacobian_and_the_condition_numbers():
    """Referee round 65 (X04, X05, X07): `vp_jacobian` is variable projection's Jacobian (against central differences); `b2_job`'s
    condition numbers are those of the joint Jacobian with the least-squares amplitudes and of the analytic one, the centres in units of
    σ, and the ratio of their smallest singular values."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    c = np.array([-0.15, 0.03, 0.2])

    def vp(cc):
        P = K.phi(cc)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    h = 1e-6
    Jd = np.stack([(vp(c + h * e) - vp(c - h * e)) / (2 * h) for e in np.eye(3)], axis=1)
    J = WM.vp_jacobian(K, f, c)
    assert np.linalg.norm(J - Jd) < 1e-6 * np.linalg.norm(J)
    q = WM.b2_job((1.0, 'base', 3, c))
    wl = np.linalg.lstsq(K.phi(c), f, rcond=None)[0]
    sj = np.linalg.svd(np.c_[K.phi(c), K.dphi(c) * wl[None, :] * WM.SIG], compute_uv=False)
    sp = np.linalg.svd(J * WM.SIG, compute_uv=False)
    assert abs(q['cond_joint'] / (sj[0] / sj[-1]) - 1) < 1e-9 and abs(q['cond_vp'] / (sp[0] / sp[-1]) - 1) < 1e-9
    assert abs(q['smin_ratio'] / (sp[-1] / sj[-1]) - 1) < 1e-9
    Jj = np.c_[K.phi(c), K.dphi(c) * wl[None, :]]
    sx, px = np.linalg.svd(Jj, compute_uv=False), np.linalg.svd(J, compute_uv=False)
    assert abs(q['cond_joint_x'] / (sx[0] / sx[-1]) - 1) < 1e-9 and abs(q['cond_vp_x'] / (px[0] / px[-1]) - 1) < 1e-9
    assert abs(q['smin_ratio_x'] / (px[-1] / sx[-1]) - 1) < 1e-9 and q['cond_joint_x'] != q['cond_joint']
    e1 = np.linalg.svd(Jj / np.linalg.norm(Jj, axis=0), compute_uv=False)
    e2 = np.linalg.svd(J / np.linalg.norm(J, axis=0), compute_uv=False)
    assert abs(q['smin_ratio_eq'] / (e2[-1] / e1[-1]) - 1) < 1e-9


def test_main_variants_merges_the_new_paths(monkeypatch, tmp_path):
    """Referee round 65 (X28–X30): `--variants i,j` runs B for those variants only (in VARIANTS_B's order) and merges them into the raw's B."""
    import exp_cohesion_leak as E
    v = [x[0] for x in WM.VARIANTS_B]
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), B=[dict(r=1.0, var=v[0], x=0), dict(r=1.0, var=v[5], x=0)])
    calls, saved = [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'b_job', lambda a: (calls.append(a), dict(r=a[0], var=a[1], x=1))[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    import multiprocessing

    class FakePool:
        def __init__(self, n):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            return [fn(j) for j in jobs]
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'B', '--variants', '0,5'])
    WM.main()
    assert calls == [(r_, v[0]) for r_ in WM.TAUS_B] + [(r_, v[5]) for r_ in WM.TAUS_B]
    B = saved[-1]['B']
    assert len(B) == 2 * len(WM.TAUS_B) and all(b['x'] == 1 for b in B) and [b['var'] for b in B][:len(WM.TAUS_B)] == [v[0]] * len(WM.TAUS_B)


def test_b2_line_at_its_boundaries():
    """Referee round 65 (X09, X14): a miss exactly 1.01 times the optimum is far, not near."""
    B2 = [dict(r=1.0, var='base', k=10, e_block=0.505, e_rule=0.5, dc=0.0)]
    md = WM.b2_line(B2)[0]
    assert '1 ends 1.01 times above it' in md and 'within 1% of it' not in md


def test_the_entries_at_their_boundaries():
    """Referee round 65 (X23, X24): an entrant at the block's outermost node is between its nodes (d = 0), and a rate within a relative
    10⁻¹² of the largest counts as the largest."""
    blk = lambda k, c: dict(k=k, tau=float(k), c=np.array(c), w=np.full(k, 0.1), e=1e-3, e_opt=1e-3, c_opt=None, relaxed=True, dc=0.0,
                            dc_hs=0.0, e_hs=1e-3)
    B = [dict(r=1.0, var='base', blocks=[blk(1, [0.0]), blk(2, [-0.1, 0.1])],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.0), (3.0, 12, 'enter', 0.1)],
              entries=[dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5), dict(it=1, j=11, u0=0.9, rate=0.4 * (1 - 5e-13), rate_max=0.4),
                       dict(it=2, j=12, u0=0.95, rate=0.1, rate_max=0.4)])]
    md = '\n'.join(WM.entry_lines(B))
    assert '| an optimal block | 2 (2) | 1 |' in md and '| — | 2 (0.900–0.950) |' in md


def test_frame_series_against_the_exact_defects():
    """Referee round 81 (D3cb-5): `gh_nodes` is the Gauss–Hermite rule of N(0, 1) (He_K at its nodes zero, the weights summing to 1, to
    the arithmetic's 110 digits); `frame_series` at rational ρ and v is Σ_{n ≥ 2K} (ρⁿ/n!)D_n(v)² with `gh_defect_v`'s exact defects."""
    from decimal import Decimal, localcontext
    from fractions import Fraction
    from math import factorial
    for K in (1, 3, 6):
        zs, ws = WM.gh_nodes(K)
        with localcontext() as ctx:
            ctx.prec = WM.A5_PREC
            assert abs(sum(ws) - 1) < Decimal(10) ** -100 and all(abs(WM.he_all(z, K)[K]) < Decimal(10) ** -90 for z in zs)
    for K, rho, v in ((1, Fraction(1, 3), Fraction(3, 2)), (2, Fraction(1, 4), Fraction(2)), (3, Fraction(2, 5), Fraction(5, 4))):
        with localcontext() as ctx:
            ctx.prec = WM.A5_PREC
            rd, vd = Decimal(rho.numerator) / rho.denominator, Decimal(v.numerator) / v.denominator
        S, terms, low = WM.frame_series(K, rd, vd)
        ref = [rho ** n / factorial(n) * WM.gh_defect_v(n, K, v) ** 2 for n in range(2 * K, 2 * K + 2 * len(terms), 2)]
        assert all(abs(Fraction(t) / r_ - 1) < Fraction(1, 10 ** 60) for t, r_ in zip(terms, ref) if r_ != 0)
        assert low < Decimal(10) ** -80


def test_the_frames_and_the_first_tail_ratio():
    """Referee round 65: `gh_defect_v` at v = 2 is `gh_defect`, and the first tail ratio is −vK(K + 2) + (v − 1)(K + 1)(2K + 1) exactly;
    the frame at ρ = t is the half-squeeze rule's (v_ρ = τ²σ/σ*, v′ = 2), at ρ = λ v′ = 1, at ρ = q the variance 2τ²/(β² + 1)."""
    from fractions import Fraction
    for K in range(1, 8):
        for n in range(2 * K, 2 * K + 6, 2):
            assert WM.gh_defect_v(n, K, 2) == W.gh_defect(n, K)
        for v in (Fraction(1, 2), Fraction(5, 3)):
            assert WM.gh_defect_v(2 * K + 2, K, v) / WM.gh_defect_v(2 * K, K, v) == -v * K * (K + 2) + (v - 1) * (K + 1) * (2 * K + 1)
    for r_ in (0.25, 1.0, 2.0):
        ss = WM.sigma_star(r_)
        P = W.params(WM.SIG, ss)
        tau2 = P['tau'] ** 2
        v, vp = WM.frame(P['t'], WM.SIG, tau2)
        assert abs(v / (tau2 * WM.SIG / ss) - 1) < 1e-12 and abs(vp - 2) < 1e-12
        assert abs(WM.frame(P['lam'], WM.SIG, tau2)[1] - 1) < 1e-12
        b2 = np.sqrt(1 + 2 * tau2 / WM.SIG ** 2)
        assert abs(WM.frame((b2 - 1) / (b2 + 1), WM.SIG, tau2)[0] / (2 * tau2 / (b2 + 1)) - 1) < 1e-12


def test_b3_restarts_the_polish_from_its_recorded_start(monkeypatch):
    """Referee rounds 65, 70: `b3_job` finds the path's first non-optimal block and restarts its first polish from that polish's recorded
    start, one setting at a time — (the path's budget, its gradient test), (B3_BIG, the test), (the path's budget, none), (B3_BIG, none) —,
    with a variant that has a gradient test; `seed_rule` judges each end from the start; `path_polishes` records the calls in the path's
    order and restores `polish_block`."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    rf = np.sqrt(np.mean(f ** 2))
    n = 4
    cs = np.linspace(-0.3, 0.3, n)
    hist = [(0.0, cs.copy(), np.zeros(n), 1.0)]
    c1 = cs.copy(); c1[1] = 0.03; w1 = np.zeros(n); w1[1] = 0.8
    r1 = np.sqrt(np.mean((f - K.phi(np.array([0.03])) @ np.array([0.8])) ** 2))
    hist.append((1.0, c1, w1, r1))                                     # block 0: one kernel, not the optimum
    R = dict(hist=hist, relaxed=[True], polish_log=[[dict(status=0, nfev=1, accepted=True)]])
    start = dict(w=np.array([0.7]), c=np.array([0.0]), s=np.array([1.0]), max_nfev=5000, gtol=1e-15)
    monkeypatch.setattr(WM, 'path_polishes', lambda r_, var: (R, K, f, [start]))
    calls = []

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append((max_nfev, gtol, tuple(w_), tuple(c_)))
        return w_ * 1.1, c_ + 0.01, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    v = WM.VARIANTS_B[3]
    q = WM.b3_job((1.0, v[0]))
    assert [x[:2] for x in calls] == [(v[3], v[5]), (WM.B3_BIG, v[5]), (v[3], None), (WM.B3_BIG, None)]
    assert all(x[2:] == ((0.7,), (0.0,)) for x in calls)
    e_opt = min(WM.sample_fit(K.phi, f, c1[[1]])[2], WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(1, WM.SIG, WM.sigma_star(1.0))[0])[2]) / rf
    assert np.isclose(q['e_opt'], e_opt, atol=0) and np.isclose(q['e_block'], r1 / rf / e_opt, atol=0)
    assert np.isclose(q['e_start'], np.sqrt(np.mean((f - K.phi(np.array([0.0])) @ np.array([0.7])) ** 2)) / rf / e_opt, atol=0)
    sr = WM.seed_rule(K, f, np.array([0.7]), np.array([0.0]), np.array([1.0]), np.array([0.01]), v[1])
    assert q['path']['move'] == sr['move'] and q['path']['accepted'] == sr['accepted']


def test_path_polishes_records_and_restores(monkeypatch):
    """Referee round 70: `path_polishes` runs the variant's path with `polish_block` wrapped — every call's start and settings recorded
    in order, the real solve's result passed through — and restores `polish_block` even when the path fails."""
    from src.core import rich_transit as RT
    real, seen = RT.polish_block, {}

    def fake_path(L, K, f, c0, rho, **kw):
        seen.update(kw)
        out = [RT.polish_block(K, f, np.array([0.5]), np.array([0.01 * i]), np.array([1.0]), max_nfev=7 + i, info=True, gtol=None)
               for i in range(2)]
        return dict(out=out)
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    v = WM.VARIANTS_B[4]
    R, K, f, calls = WM.path_polishes(1.0, v[0])
    assert RT.polish_block is real and len(calls) == 2
    assert [(c['max_nfev'], c['gtol'], float(c['c'][0])) for c in calls] == [(7, None, 0.0), (8, None, 0.01)]
    assert seen['polish_step'] == v[1] and seen['polish_nfev'] == v[3] and seen['repolish'] == v[4] and seen['polish_gtol'] == v[5]
    assert len(R['out'][0]) == 3 and R['out'][0][2]['nfev'] <= 7

    def boom(*a, **k):
        raise RuntimeError('path failed')
    monkeypatch.setattr(RT, 'sliding_iss', boom)
    try:
        WM.path_polishes(1.0, v[0])
    except RuntimeError:
        pass
    assert RT.polish_block is real


def test_b3_binds_and_lines_on_synthetic_results():
    """Referee rounds 70, 77: `b3_binds` labels from the ends that reach the optimum and are accepted — the budget, the gradient test,
    either, both —, "(and the radius)" where every reaching end is rejected, the radius alone, "neither" with the stop, the
    inconsistency —, judged against its ref; B3's table and summary, below B4_MAX against B4's long-double end."""
    P = lambda ra, acc=True, st=0, nf=5000, mv=0.5: dict(status=st, nfev=nf, ratio=ra, move=mv, accepted=acc)
    mk = lambda r_, var, e_opt, a, b, c, d: dict(r=r_, var=var, k=6, e_opt=e_opt, e_start=17.0, e_block=1.46, n_polish=1,
                                                 **dict(zip(WM.B3_KEYS, (a, b, c, d))))
    B3 = [mk(0.5, 'polish 2σ', 6e-8, P(1.46), P(1.0, st=1, nf=6635), P(1.46), P(1.0, st=3, nf=6654)),
          mk(0.25, 'base', 6e-8, P(1.002, st=1), P(1.0017, st=1), P(1.0000007, st=3), P(1.0, st=3)),
          mk(1.0, 'grid +0.05', 3e-7, P(1.0, acc=False, mv=1.98), P(1.0, acc=False, mv=1.98), P(1.0, acc=False), P(1.0, acc=False)),
          mk(0.75, 'x', 2e-9, P(1.3), P(1.2), P(1.2), P(1.1)),
          mk(1.5, 'y', 4e-9, P(1.0), P(1.0), P(1.0), P(1.0)),
          mk(0.25, 'z', 2.5e-11, P(3.7), P(3.7), P(3.7), P(1.0147))]
    labs = [WM.b3_binds(q) for q in B3]
    assert labs == ['the budget', 'the gradient test', 'the radius', 'neither (the budget spent)',
                    "nothing (inconsistent: the path's own polish reaches it, accepted)", 'neither (the budget spent)']
    alt = lambda i, **kw: dict(B3[i], **kw)
    assert WM.b3_binds(alt(0, budget=P(1.0, acc=False))) == 'the budget and the gradient test'
    assert WM.b3_binds(alt(0, **{'budget, no test': P(1.0, acc=False)})) == 'the budget'
    assert WM.b3_binds(alt(1, budget=P(1.0, acc=False))) == 'the gradient test'
    assert WM.b3_binds(alt(1, budget=P(1.0))) == 'the budget or the gradient test'
    assert WM.b3_binds(alt(1, budget=P(1.0, acc=False), **{'no test': P(1.0, acc=False), 'budget, no test': P(1.0, acc=False)})) == \
        'the budget or the gradient test (and the radius)'
    assert WM.b3_binds(alt(0, budget=P(1.2), **{'budget, no test': P(1.0, acc=False)})) == 'the budget and the gradient test (and the radius)'
    assert WM.b3_binds(alt(3, **{'budget, no test': P(1.1, st=3)})) == "neither (the solve's stop, status 3)"
    assert WM.b3_binds(alt(3, **{'budget, no test': P(1.1, st=1)})) == "neither (the solve's stop, status 1)"
    near = alt(5, **{'budget, no test': P(1.0015)})
    assert WM.b3_binds(near) == 'neither (the budget spent)' and WM.b3_binds(near, ref=1.001) == 'the budget and the gradient test'
    s = '\n'.join(WM.b3_lines(B3 + [dict(r=2.0, var='base', k=None)]))
    assert ("| 0.5 | polish 2σ | 6 | 6.00e-08 | 17 | 1.46 | 0/5000, 1.46000, 0.50σ a | 1/6635, 1.00000, 0.50σ a | 0/5000, 1.46000, 0.50σ a | "
            "3/6654, 1.00000, 0.50σ a | the budget |") in s
    assert "| 0.25 | z | 6 | 2.50e-11† |" in s and s.endswith('\n')
    assert ("\n\nBelow it (1 blocks) the restart with both changed stops at 1.0147 (0/5000: neither (the budget spent)) times the stored "
            "optimum.") in s
    assert "Above 1e-10 (5 blocks: τ/σ (variant)) the end is bounded by" in s
    assert ("- neither (the budget spent): 1 — 0.75 (x);\n- nothing (inconsistent: the path's own polish reaches it, accepted): 1 — 1.5 (y);"
            "\n- the budget: 1 — 0.5 (polish 2σ);\n- the gradient test: 1 — 0.25 (base);\n- the radius: 1 — 1 (grid +0.05).") in s
    assert ("by the optimum (one solver from one kind of start; above B4_MAX a negative value is the end below the stored optimum, "
            "converged only to that): 3e-07: +0.0e+00; 6e-08: +0.0e+00; 4e-09: +0.0e+00; 2e-09: +1.0e-01; 2.5e-11: +1.5e-02.") in s
    B4 = [dict(r=0.25, k=6, start='block', e_opt=2.5e-11, **{'long double': dict(ratio=0.99, rms='2.475e-11')}),
          dict(r=0.25, k=6, start='optimum', e_opt=2.5e-11, **{'long double': dict(ratio=0.995, rms='2.4875e-11')})]
    s4 = '\n'.join(WM.b3_lines([B3[5], alt(5, r=0.75, var='w', **{'budget, no test': P(0.99)})], B4))
    assert "stops at 1.0147 (0/5000: neither (the budget spent)), 0.9900 (0/5000: the budget and the gradient test) times the stored optimum." in s4
    s4 = '\n'.join(WM.b3_lines([B3[5]], B4))
    assert "stops at 1.0147 (0/5000: neither (the budget spent)) times the stored optimum, 1.0249 times B4's best end." in s4
    assert "2.5e-11: +2.5e-02." in s4
    z2 = alt(5, var='z2', e_opt=2.6e-11, **{'budget, no test': P(1.0)})  # referee round 81 (D3cb-6): the reference over each row's own optimum
    s4 = '\n'.join(WM.b3_lines([B3[5], z2], B4))
    assert f"1.0249, {1.0 / (0.99 * 2.5e-11 / 2.6e-11):.4f} times B4's best end." in s4
    Bs = [dict(r=0.25, var='z', blocks=[dict(k=6, e=1.1, e_opt=1.0, polish=[dict(status=0, nfev=5000)])]),
          dict(r=0.25, var='z2', blocks=[dict(k=6, e=1.1, e_opt=1.0, polish=[dict(status=1, nfev=77)])])]
    s5 = '\n'.join(WM.b3_lines([B3[5], z2], B4, Bs))
    assert "The 'path' restart reproduces the first polish B stored — its stop and evaluations — at 1 of 2 blocks (not at 0.25 (z2))." in s5
    q2 = [dict(B3[1]), dict(B3[0])]
    q2[0]['budget, no test'] = P(0.99998, st=3)
    s2 = '\n'.join(WM.b3_lines(q2))
    assert "6e-08: -2.0e-05 … +0.0e+00." in s2
    assert WM.b3_lines([]) == [] and WM.b3_lines(None) == []


def test_seed_rule_with_two_seed_kernels_and_its_default():
    """Referee round 70 (S02, S05–S08): with two kernels near their seeds the rule takes the larger of their moves, each the nearer of its
    start and its ascent target — one decided by the target —; the seed set is |w| ≤ 10⁻²·max|w|, `sliding_iss`'s default."""
    import inspect
    from src.core import rich_transit as RT
    assert (inspect.signature(WM.seed_rule).parameters['wfrac'].default
            == inspect.signature(RT.sliding_iss).parameters['polish_wfrac'].default == 1e-2)
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    c, w = np.array([-0.4, 0.03, 0.5]), np.array([5e-3, 0.8, 7e-3])
    r0 = f - K.phi(c) @ w
    t0, t2 = (RT.ascent_target(K, r0, c[i], 1.0)[0] for i in (0, 2))
    cp = c.copy(); cp[0] = t0 + 0.2 * WM.SIG; cp[2] = c[2] + 0.7 * WM.SIG
    r = WM.seed_rule(K, f, w, c, np.ones(3), cp, 1.0)
    m0 = min(abs(cp[0] - c[0]), abs(cp[0] - t0)) / WM.SIG
    m2 = min(abs(cp[2] - c[2]), abs(cp[2] - t2)) / WM.SIG
    assert abs(cp[0] - c[0]) / WM.SIG > 1.0 and np.isclose(m0, 0.2)
    assert np.isclose(r['move'], max(m0, m2)) and r['accepted']
    assert not WM.seed_rule(K, f, w, c, np.ones(3), cp, 0.1)['accepted']
    assert WM.seed_rule(K, f, np.array([5e-3, 0.8, 0.5]), c, np.ones(3), cp, 1.0)['move'] == m0


def test_part_a5_on_its_own_numbers():
    """Referee rounds 70, 81 (A01–A05, A09; D3cb-5): A5's frames at ρ = t, q, λ and the Mercer frame's rule N(0, 2τ²/(β² + 1)) — not the
    Mercer-zero rule — with its errors at the given K, by hand; the tail-ratio checks all exact; the exact frame series (`a5_series`) is
    the rule's error (its ratio per kernel the actual errors' to 10⁻⁴), its tail above the first term from K = 6 at τ/σ = 2, and its
    lines."""
    A = dict(rows=[dict(r=1.0, K=3), dict(r=1.0, K=4), dict(r=2.0, K=5)])
    A5 = WM.part_A5(A, far=False)
    x = next(r for r in A5['rows'] if r['r'] == 1.0)
    ss = WM.sigma_star(1.0)
    P = W.params(WM.SIG, ss)
    tau2 = P['tau'] ** 2
    b2 = np.sqrt(1 + 2 * tau2 / WM.SIG ** 2)
    vhs, vM = tau2 * WM.SIG / ss, 2 * tau2 / (b2 + 1)
    assert x['Ks'] == [3, 4] and np.isclose(x['vM'], vM / vhs) and np.isclose(x['q'], (b2 - 1) / (b2 + 1))
    assert np.isclose(x['frames']['t']['v'], 1.0) and np.isclose(x['frames']['λ']['vp'], 1.0)
    assert np.allclose(x['e_mframe'], [W.l2_error(*W.scaled_rule(K, WM.SIG, ss, vM), WM.SIG, ss) for K in (3, 4)])
    assert all(c_['ok'] for c_ in A5['check'])
    sr = {y['r']: y for y in A5['series']}
    by = {y['K']: y for y in sr[1.0]['frames']['Mercer']}
    rat = lambda K: float(np.sqrt(by[K - 1]['S'] / by[K]['S']))
    e = [W.l2_error(*W.scaled_rule(K, WM.SIG, ss, vM), WM.SIG, ss) for K in (3, 4)]
    assert abs(rat(4) / (e[0] / e[1]) - 1) < 1e-4                       # referee round 81: the exact series is the rule's error
    assert sr[2.0]['frames']['Mercer'][4]['tail'] < 1 < sr[2.0]['frames']['Mercer'][5]['tail']   # the tail above the first from K = 6
    A5['rows'] = [dict(x, Ks=[3, 4, 5, 6], e_mframe=[1.0, 0.25, 0.05, 0.01]) for x in A5['rows']]
    md = '\n'.join(WM.a5_lines(A5))
    assert f"| 4.0000, 5.0000, 5.0000 at K = 4–6 ({1 / x['t']:.4g}, {1 / x['q']:.4g}) |" in md and "ratio per kernel at the last three K" in md
    assert "2/(q(2 − v′))" not in md and "| 2 | Mercer | " in md and "| 1 | half-squeeze | " in md
    assert (f"| 1 | Mercer | " + ', '.join(f"{by[K]['tail']:.3g}" for K in (5, 10, 20, 40)) + " | 8 | "
            + ', '.join(f"{rat(K):.4f}" for K in (10, 20, 30, 40))) in md


def test_b2_job_where_variable_projection_leaves_the_block():
    """Referee round 70 (B04–B06, B09, V02, V03): `b2_job` from a two-kernel block away from the optimum — variable projection from the
    block and from the half-squeeze nodes, the better end's distance, the stop, the closest nodes and the stationarity cosine — against a
    direct `vp_fit`; `vp_fit` returns sorted centres, their least-squares coefficients, the rms, the status and the evaluations."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    r_ = 1.0
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
    rf = np.sqrt(np.mean(f ** 2))
    c = np.array([0.2, -0.15])
    cv, av, ev, st, nf = WM.vp_fit(K.phi, f, c)
    assert np.all(np.diff(cv) > 0) and np.allclose(av, np.linalg.lstsq(K.phi(cv), f, rcond=None)[0])
    assert np.isclose(ev, np.sqrt(np.mean((f - K.phi(cv) @ av) ** 2)), atol=0) and isinstance(st, int) and isinstance(nf, int) and nf > 0
    q = WM.b2_job((r_, 'base', 2, c))
    yh, _ = W.half_squeeze(2, WM.SIG, WM.sigma_star(r_))
    cr, _, er, _, _ = WM.vp_fit(K.phi, f, WM.C_B + yh)
    assert np.isclose(q['e_block'], ev / rf, atol=0) and np.isclose(q['e_rule'], er / rf, atol=0) and q['vp_status'] == st and q['vp_nfev'] == nf
    assert np.isclose(q['dc'], np.abs(np.sort(c) - (cv if ev <= er else cr)).max()) and q['dc'] > 0.01
    assert np.isclose(q['vp_gap'], np.diff(cv).min() / WM.SIG)
    Je = WM.vp_jacobian(K, f, cv)
    re = f - K.phi(cv) @ av
    assert np.isclose(q['vp_cos'], np.max(np.abs(Je.T @ re) / (np.linalg.norm(Je, axis=0) * np.linalg.norm(re))))


def test_b3_jobs_and_run_b_with_variants(monkeypatch):
    """Referee round 70 (P01, X02): B3's jobs are the paths with a non-optimal block; `run_B` with variants reruns those and merges them
    into the stored paths (none dropped), without variants reruns all."""
    o = lambda e: dict(e=e, e_opt=1.0)
    B = [dict(r=0.5, var='base', blocks=[o(1.0), o(1.0005)]), dict(r=1.0, var='base', blocks=[o(1.0), o(1.2)]),
         dict(r=0.5, var='polish 2σ', blocks=[o(2.0)])]
    assert WM.b3_jobs(B) == [(1.0, 'base'), (0.5, 'polish 2σ')]
    new = [dict(r=0.5, var='polish 2σ', blocks=[o(1.0)])]
    monkeypatch.setattr(WM, 'part_B', lambda idx=None: new if idx is not None else ['all'])
    got = WM.run_B(B, {1})
    assert len(got) == 3 and next(b for b in got if b['var'] == 'polish 2σ')['blocks'] == [o(1.0)]
    assert [(b['var'], b['r']) for b in got] == [('base', 0.5), ('base', 1.0), ('polish 2σ', 0.5)]
    assert WM.run_B(B, None) == ['all']


def test_b3_stops_the_path_at_the_stored_blocks_first_polish(monkeypatch):
    """Referee rounds 70, 77 (B3's cost; FN01, BR08): `first_nonoptimal` finds in B's stored blocks the first polished non-optimal block,
    the number of polish calls before it (unequal counts, a block without polishes) and its first polish's stored stop; `path_polishes`
    with stop_after = n records call n and stops the path before its solve (R None, the solve not called, `polish_block` restored);
    `b3_job` with the block restarts that recorded start, with the stored block's errors over its optimum's, and refuses a rerun whose call
    n is not the block's size or whose 'path' restart is not the stored first polish."""
    from src.core import rich_transit as RT
    o = lambda k, e, n: dict(k=k, e=e, e_opt=0.5, polish=[dict(status=1, nfev=1, accepted=True)] * n)
    b = dict(r=0.5, var='base', blocks=[o(1, 0.5, 3), o(2, 1.5, 0), o(3, 0.75, 3), o(4, 1.0, 1)])
    assert WM.first_nonoptimal(b) == dict(n0=3, k=3, e=0.75, e_opt=0.5, n_polish=3, first=dict(status=1, nfev=1, accepted=True))
    assert WM.first_nonoptimal(dict(blocks=[o(1, 0.5, 2), o(2, 0.50025, 1)])) is None
    assert WM.b3_jobs([b, dict(r=1.0, var='base', blocks=[o(1, 0.5, 1)])]) == [(0.5, 'base', WM.first_nonoptimal(b))]
    solved = []

    def real_stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        solved.append(c_.size)
        return w_, c_, dict(status=1, nfev=1)
    monkeypatch.setattr(RT, 'polish_block', real_stub)

    def fake_path(L, K, f, c0, rho, **kw):
        for k in (1, 1, 1, 3, 3, 4):
            RT.polish_block(K, f, np.full(k, 0.5), np.linspace(-0.1, 0.1, k), np.ones(k), info=True)
        return 'full'
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    R, K, f, calls = WM.path_polishes(0.5, 'base', stop_after=3)
    assert R is None and len(calls) == 4 and calls[3]['c'].size == 3 and solved == [1, 1, 1] and RT.polish_block is real_stub
    R, K, f, calls = WM.path_polishes(0.5, 'base')
    assert R == 'full' and len(calls) == 6 and solved == [1, 1, 1, 1, 1, 1, 3, 3, 4]
    q = WM.b3_job((0.5, 'base', WM.first_nonoptimal(b)))
    assert q['k'] == 3 and q['e_opt'] == 0.5 and q['e_block'] == 1.5 and q['n_polish'] == 3
    assert q['settings'] == dict(nfev=5000, gtol=1e-15) and all(q[k]['nfev'] == 1 for k in WM.B3_KEYS)
    rf = np.sqrt(np.mean(f ** 2))
    st = calls[3]
    assert np.isclose(q['e_start'], np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2)) / rf / 0.5, atol=0)
    assert np.isclose(q['path']['ratio'], np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2)) / rf / 0.5, atol=0)
    for bad, msg in ((dict(WM.first_nonoptimal(b), k=4), 'not the stored block'),
                     (dict(WM.first_nonoptimal(b), first=dict(status=0, nfev=1)), 'is not the stored first polish'),
                     (dict(WM.first_nonoptimal(b), first=dict(status=1, nfev=2)), 'is not the stored first polish')):
        try:
            WM.b3_job((0.5, 'base', bad))
            raise AssertionError('a mismatched block was accepted')
        except RuntimeError as e:
            assert msg in str(e)
    assert WM.b3_job((0.5, 'base', dict(WM.first_nonoptimal(b), first=dict(status=1))))['k'] == 3
    solved.clear()                                                      # referee round 81: acceptance compared, before the larger restarts
    try:
        WM.b3_job((0.5, 'base', dict(WM.first_nonoptimal(b), first=dict(status=1, nfev=1, accepted=False))))
        raise AssertionError('a mismatched acceptance was accepted')
    except RuntimeError as e:
        assert 'is not the stored first polish' in str(e) and solved == [1, 1, 1, 3]


def test_lm_fit_in_double_and_long_double():
    """Referee rounds 70, 77 (B4): `lm_fit` from variable projection's end at τ/σ = 1/2, K = 8 (its optimum ≈ 2.1·10⁻¹⁰): in double the
    solve stalls with nothing gained (its steps lost to the residual's rounding), in long double the same solve goes on 10⁻⁵ below, and so
    does the solve in double with the residual about an anchor; the iteration cap stops it; the anchored solve's marks, and (referee
    round 81) the anchored solve in long double from the anchored double's end, ending no higher; on an exact two-kernel sum from a
    perturbed start both reach the sum; `lm_rms` is the long-double rms; the defaults are B4's."""
    import inspect
    import pytest
    from src.core import rich_transit as RT
    if np.finfo(np.longdouble).eps >= 1e-18:
        pytest.skip('long double is not extended on this platform')
    K = RT.Bumps(WM.X_B, WM.SIG)
    ss = WM.sigma_star(0.5)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    co, ao, eo = WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(8, WM.SIG, ss)[0])
    pd, itd, accd, std = WM.lm_fit(f, ao, co, np.float64)
    pl, itl, accl, stl = WM.lm_fit(f, ao, co, np.longdouble)
    ed, el = WM.lm_rms(f, pd, 8) / eo, WM.lm_rms(f, pl, 8) / eo
    assert pd.dtype == np.float64 and pl.dtype == np.longdouble and std == stl == 'stall' and accd <= 5 and accl > 100
    assert ed > 1 - 1e-6 and el < ed - 5e-6
    pa, ita, acca, sta = WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=1000)
    ea = WM.lm_rms(f, pa, 8) / eo
    assert pa.dtype == np.float64 and ita == 1000 and sta == 'cap' and acca > 100 and abs(ea / el - 1) < 1e-8
    marks = {10: None, 10 ** 9: None}
    WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=50, marks=marks)
    assert np.array_equal(marks[10], WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=10)[0]) and marks[10 ** 9] is None
    pq, itq, accq, stq = WM.lm_fit(f, pa[:8], pa[8:], np.longdouble, anchored=True, max_iter=300)
    assert pq.dtype == np.longdouble and WM.exact_rms(f, pq, 8) <= WM.exact_rms(f, pa, 8) * (1 + 1e-12)
    p, it, acc, st = WM.lm_fit(f, ao, co, np.longdouble, max_iter=50)
    assert it == 50 and st == 'cap' and WM.lm_rms(f, p, 8) / eo > el
    w, c = np.array([0.5, 0.3]), np.array([-0.05, 0.2])
    f2 = K.phi(c) @ w
    for kw in (dict(), dict(anchored=True)):
        p2, it2, acc2, st2 = WM.lm_fit(f2, w * np.array([1.02, 0.97]), c + np.array([0.01, -0.008]), np.float64, **kw)
        assert st2 == 'stall' and acc2 > 3 and np.allclose(p2, np.r_[w, c], rtol=0, atol=1e-12) and WM.lm_rms(f2, p2, 2) < 1e-15
    xl = WM.X_B.astype(np.longdouble)
    pr = np.r_[w, c].astype(np.longdouble)
    rl = (f2 + 1e-3).astype(np.longdouble) - np.exp(-(xl[:, None] - pr[None, 2:]) ** 2 / (2 * np.longdouble(WM.SIG) ** 2)) @ pr[:2]
    assert WM.lm_rms(f2 + 1e-3, np.r_[w, c], 2) == float(np.sqrt(np.mean(rl * rl)))
    sig = inspect.signature(WM.lm_fit).parameters
    assert sig['max_iter'].default == WM.B4_CAP == 400_000 and sig['mu0'].default == 1e-20 and WM.B4_MAX == 2e-9
    assert sig['anchored'].default is False and sig['anchor_every'].default == WM.B4_ANCHOR == 2000 and WM.B4_FLAT == 1e-10
    assert sig['stall_reanchor'].default is True and WM.EXACT_PREC == 60
    assert WM.B4_MARKS == (100_000, 200_000, 300_000) and sig['marks'].default is None
    assert [k for k, _ in WM.B4_PREC] == ['double', 'long double'] and WM.B4_PREC[1][1] is np.longdouble


def test_b4_jobs_job_and_lines(monkeypatch):
    """Referee rounds 70, 77, 81, 84 (B4; BS06, BS07; D3cb-2; D3ce-2, D3ce-4, D3ce-5): the jobs are B3's distinct blocks below B4_MAX (its
    first row for each — not the first path of that width —, the block of that K nearest its optimum), from the block and from the stored
    optimum with least-squares amplitudes; `b4_job` (its solves faked, each call's values distinct) reports the start's and each end's
    exact and long-double error over the optimum's, its exact rms and its Gauss–Newton gap, for both arithmetics, the anchored double
    with its marks, its floors at its end and at its last anchor (from the solve's anchor log) and its anchors, the solve without stall
    re-anchoring with its anchors, the long double's anchored solve from the anchored double's end with its parameters rounded to double,
    and the polish's TRF without the gradient test; `b4_lines` prints the table against each block's best end — the excesses exact —
    and the summaries."""
    from decimal import Decimal, localcontext
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    o = lambda k, e_opt, c: dict(k=k, e=2 * e_opt, e_opt=e_opt, c=np.array(c), w=np.full(len(c), 0.4), c_opt=np.array(c) + 0.01)
    B = [dict(r=0.5, var='v1', blocks=[o(1, 1e-3, [0.0]), o(2, 3e-10, [-0.1, 0.1]), o(2, 5e-10, [-0.2, 0.2])]),
         dict(r=0.5, var='v2', blocks=[o(2, 3e-10, [-0.15, 0.15])]),
         dict(r=1.0, var='v1', blocks=[o(3, 1.5e-9, [-0.1, 0.0, 0.1])]),
         dict(r=1.0, var='v2', blocks=[o(3, 1e-8, [-0.12, 0.0, 0.12])])]
    q = lambda r_, var, k, e_opt: dict(r=r_, var=var, k=k, e_opt=e_opt)
    B3 = [q(0.5, 'v2', 2, 3e-10), q(0.5, 'v1', 2, 5e-10), q(1.0, 'v2', 3, 1e-8), q(1.0, 'v1', 3, 1.5e-9), dict(r=2.0, var='v1', k=None)]
    jobs = WM.b4_jobs(B, B3)
    assert [(j[0], j[1], j[2], j[3]) for j in jobs] == [(0.5, 'v2', 2, 'block'), (0.5, 'v2', 2, 'optimum'), (1.0, 'v1', 3, 'block'),
                                                         (1.0, 'v1', 3, 'optimum')]
    assert np.array_equal(jobs[0][5], [-0.15, 0.15]) and np.allclose(jobs[1][5], [-0.14, 0.16], rtol=0, atol=1e-15) and jobs[0][6] == 3e-10
    assert np.array_equal(jobs[2][5], [-0.1, 0.0, 0.1]) and jobs[2][6] == 1.5e-9
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.5) ** 2))
    assert np.allclose(jobs[1][4], np.linalg.lstsq(K.phi(np.array([-0.14, 0.16])), f, rcond=None)[0])
    assert WM.b4_jobs(B, None) == [] and WM.b4_jobs(B, [q(0.5, 'v1', 2, 2e-9)]) == []
    seen, pol, logs = [], [], {}
    calls = iter(range(10))

    def fake(f_, w_, c_, dtype, max_iter=WM.B4_CAP, mu0=1e-20, anchored=False, anchor_every=WM.B4_ANCHOR, marks=None, stall_reanchor=True,
             kappa=None, anchor_log=None):
        i = next(calls)                                                 # referee round 81: distinct values per call, moved amplitudes
        seen.append((dtype, anchored, marks is not None, stall_reanchor, anchor_log is not None, kappa, np.r_[w_, c_].astype(float)))
        if marks is not None:
            marks[WM.B4_MARKS[0]] = np.r_[w_, c_ + 0.002]
        if anchor_log is not None:                                      # referee round 84: the last anchor, its frozen r₀ planted
            pa_ = np.r_[w_, c_ + 0.0005 * (i + 1)]
            anchor_log.update(stall=3 + i, period=i, p=pa_, r0=WM.direct_double(f_, pa_, len(w_)) + 1e-16 * (i + 1), it=4 + i, kind='stall')
            logs[i] = dict(anchor_log)
        return np.r_[w_ * (1 + 0.001 * i), c_ + 0.001 * (i + 1)].astype(dtype), 7 + i, 5 + i, ('stall', 'cap', 'flat')[i % 3]
    monkeypatch.setattr(WM, 'lm_fit', fake)

    def pb(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        pol.append((max_nfev, gtol, info))
        return w_ * 1.01, c_ + 0.004, dict(status=3, nfev=1234)
    monkeypatch.setattr(RT, 'polish_block', pb)
    rf = float(np.sqrt(np.mean(f ** 2)))
    w0, c0 = jobs[0][4], jobs[0][5]
    ex = lambda w_, c_: WM.exact_rms(f, np.r_[w_, c_], 2) / rf / 3e-10
    out = WM.b4_job(jobs[0])
    assert [x[:6] for x in seen] == [(np.float64, False, False, True, False, None), (np.longdouble, False, False, True, False, None),
                                     (np.float64, True, True, True, True, None), (np.float64, True, False, False, True, None),
                                     (np.longdouble, True, False, True, False, None)] and pol == [(WM.B3_BIG, None, True)]
    pa = np.r_[w0 * 1.002, c0 + 0.003]                                  # the anchored double's end (call 2) starts the long double's
    assert np.allclose(seen[4][6], pa, rtol=0, atol=1e-15) and all(np.array_equal(x[6], np.r_[w0, c0]) for x in seen[:4])
    assert out['start'] == 'block' and out['k'] == 2 and out['e_start'] == ex(w0, c0)
    ends = {}
    for i, k_ in enumerate(('double', 'long double', 'anchored', 'nostall')):
        q = out[k_]
        ends[k_] = np.r_[w0 * (1 + 0.001 * i), c0 + 0.001 * (i + 1)].astype(np.longdouble if i == 1 else np.float64)
        assert (q['it'], q['acc'], q['stop']) == (7 + i, 5 + i, ('stall', 'cap', 'flat')[i % 3]), k_
        assert q['ratio'] == ex(w0 * (1 + 0.001 * i), c0 + 0.001 * (i + 1)), k_
        assert q['ratio_ld'] == WM.lm_rms(f, ends[k_], 2) / rf / 3e-10, k_
    P2 = np.array(out['anchored']['p'])
    pl = np.r_[P2[:2] * (1 + 0.001 * 4), P2[2:] + 0.001 * (4 + 1)].astype(np.longdouble)
    ends['anchored ld'] = pl
    ends['trf'] = np.r_[w0 * 1.01, c0 + 0.004]
    assert out['anchored ld']['ratio'] == WM.exact_rms(f, pl, 2) / rf / 3e-10
    for k_, p_ in ends.items():                                         # referee round 84: the exact rms and the gap at every end
        assert out[k_]['rms'] == str(WM.exact_resid(f, p_, 2, dec=True)[2]) and out[k_]['gap'] == WM.b4_gap(f, p_, 2), k_
    assert out['anchored ld']['rms_d'] == str(WM.exact_resid(f, pl.astype(float), 2, dec=True)[2])
    assert out['double']['p_lo'] == [0.0] * 4 and np.allclose(out['double']['p'], np.r_[w0, c0 + 0.001], rtol=0, atol=1e-15)
    assert len(out['anchored']['marks']) == 1 and out['anchored']['marks'][0][0] == WM.B4_MARKS[0]
    assert out['anchored']['marks'][0][1] == ex(w0, c0 + 0.002)
    seen.clear()
    monkeypatch.setattr(WM, 'lm_rms', lambda f_, p_, k_: 123.0)       # the marks are exact, not lm_rms's (BJ05)
    out2 = WM.b4_job(jobs[0])
    assert out2['anchored']['marks'] == out['anchored']['marks'] and out2['anchored']['ratio_ld'] == 123.0 / rf / 3e-10
    fl = WM.b4_floors(f, ends['anchored'], 2, logs[2]['p'], logs[2]['r0'])
    a = out['anchored']
    assert (a['floor'], a['floor_last'], a['floor_s2'], a['gap']) == (fl['end'], fl['last'], fl['s2'], fl['gap']) and fl['last'] != fl['end']
    assert a['anchors'] == dict(stall=5, period=2, kind='stall', before=3, dist=float(np.abs(ends['anchored'] - logs[2]['p']).max()))
    assert out['nostall']['anchors'] == dict(stall=6, period=3, kind='stall', before=3,
                                             dist=float(np.abs(ends['nostall'] - logs[3]['p']).max()))
    assert out['trf']['status'] == 3 and out['trf']['nfev'] == 1234 and out['trf']['ratio'] == ex(w0 * 1.01, c0 + 0.004)
    assert out['eps'] == float(np.finfo(np.longdouble).eps)
    R = lambda ra, rl, it=100, acc=90, st='stall', gap=1e-17: dict(ratio=ra, ratio_ld=rl, it=it, acc=acc, stop=st, gap=gap)
    T = lambda ra, rl: dict(ratio=ra, ratio_ld=rl, status=3, nfev=2500, gap=1e-9)
    an = lambda st, pe, kind, bef, d: dict(stall=st, period=pe, kind=kind, before=bef, dist=d)

    def cell(r_, k, start, e_opt, e_start, d, ld, a_, ns, al, tr, fl=7e-13, last=6e-13, s2=2e-13, mk=(), anc=None, nsa=None, d_rel=1e-12):
        q = dict(r=r_, var='x', k=k, start=start, e_opt=e_opt, e_start=e_start, eps=1.0842e-19, double=d, trf=tr, nostall=dict(ns, anchors=nsa),
                 anchored=dict(a_, floor=fl, floor_last=last, floor_s2=s2, anchors=anc, marks=[list(x) for x in mk]),
                 **{'long double': ld, 'anchored ld': al})
        for k_, _ in WM.B4_ENDS:
            q[k_]['rms'] = str(Decimal(repr(q[k_]['ratio'])) * Decimal(repr(e_opt)))
        q['anchored ld']['rms_d'] = str(Decimal(q['anchored ld']['rms']) * (1 + Decimal(repr(d_rel))))
        return q
    rows = [cell(0.25, 6, 'block', 2.5e-11, 4.95, R(1.0019, 1.0019), R(0.99997, 0.99996), R(0.99996, 0.99996, 5000, 4000, 'flat', 3e-13),
                 R(0.999962, 0.999962, 900, 800), R(0.9999599, 0.9999599, 60, 50, gap=4e-18), T(1.01, 1.01), mk=((100000, 0.999961),),
                 anc=an(40, 2, 'stall', 3, 2e-9), nsa=an(0, 10, 'period', 5, 1e-6), d_rel=1.5e-12),
            cell(0.25, 6, 'optimum', 2.5e-11, 1.00003, R(0.99997, 0.99997), R(0.99998, 0.999979), R(0.99996, 0.99996, 3000, 2000, 'flat', 4e-13),
                 R(0.99996, 0.99996, 40, 30), R(0.9999599001, 0.9999599001, 70, 60), T(0.99999, 0.99999), fl=2e-13, last=3e-13, s2=1e-13,
                 anc=an(20, 1, 'stall', 7, 5e-10), nsa=an(0, 0, 'start', 40, 1e-7)),
            cell(1.0, 12, 'block', 1.6e-9, 1.47, R(0.999997, 0.999997), R(0.999996, 0.999996), R(0.999996, 0.999996, 8000, 7000, 'stall', 1e-16),
                 R(0.999996, 0.999996), R(0.999996, 0.999996), T(1.0, 1.0), anc=an(30, 3, 'stall', 1, 0.0), nsa=an(0, 75, 'period', 2, 1e-8)),
            cell(1.0, 12, 'optimum', 1.6e-9, 1.0, R(0.999996, 0.999996), R(0.999996, 0.999996), R(0.999996, 0.999996, 9000, 8000, 'flat', 2e-16),
                 R(0.999996, 0.999996), R(0.999996, 0.999996), T(0.999996, 0.999996), anc=an(10, 4, 'period', 28, 4.5e-9), nsa=an(0, 0, 'start', 9, 0.0))]
    s = '\n'.join(WM.b4_lines(rows))
    assert "in double or in long double (its ε 1.08e-19), or about an anchor" in s and "(`polish_block`, 500000 evaluations, no gradient test)" in s
    assert "a step accepted when its predicted decrease is positive and the objective" in s
    exb = lambda ra: WM.b4_excess(Decimal(repr(ra)) * Decimal(repr(2.5e-11)), Decimal(repr(0.9999599)) * Decimal(repr(2.5e-11)))
    ex4 = lambda ra: f"{exb(ra):.1e}"
    assert (f"| 0.25 | x | 6 | 2.50e-11† | block | 4.95 | {0.9999599:.12f} | {ex4(1.0019)}, 100/90, stall | {ex4(0.99997)}, 100/90, stall | "
            f"{ex4(0.99996)}, 5000/4000, flat | {ex4(0.999962)}, 900/800, stall | 0.0e+00, 60/50, stall | {ex4(1.01)}, 3/2500 |") in s
    assert f"{ex4(0.9999599001)}, 70/60, stall" in s and abs(exb(0.9999599001) - 1.0000400e-10) < 1e-15   # exact, not a double ratio's rounding
    ev = [r_[k_]['ratio_ld'] / r_[k_]['ratio'] - 1 for r_ in rows for k_, _ in WM.B4_ENDS]
    evl = [r_['long double']['ratio_ld'] / r_['long double']['ratio'] - 1 for r_ in rows]
    assert (f"against the exact one, minus 1: {min(ev):+.1e} … {max(ev):+.1e} over every end, {min(evl):+.1e} … {max(evl):+.1e} at the "
            "long-double solve's ends (2 of 4 below it") in s
    sp = abs(WM.b4_excess(rows[0]['anchored ld']['rms'], rows[1]['anchored ld']['rms']))
    assert ("The best end is reached by long double at 1 and by long double, anchored, from the anchored double's end at 1 of the 2 blocks; "
            f"that solve's ends from the two starts agree within {sp:.1e} (exactly), and a Gauss–Newton step from the best end's exact residual "
            "predicts at most 1.0e-17 of further decrease") in s
    eg = [WM.b4_excess(r_[k_]['rms'], WM.b4_best(rows)[(r_['r'], r_['k'])]['rms']) / r_[k_]['gap'] for r_ in rows for k_, _ in WM.B4_ENDS]
    eg = [x for x, r_, k_ in zip(eg, [r_ for r_ in rows for _ in WM.B4_ENDS], [k_ for _ in rows for k_, _ in WM.B4_ENDS])
          if WM.b4_excess(r_[k_]['rms'], WM.b4_best(rows)[(r_['r'], r_['k'])]['rms']) >= 1e-14]
    assert (f"The {len(eg)} of the 24 ends whose excess is at least 10³ times that lie at their own Gauss–Newton gap from it: their excess "
            f"over their gap is {min(eg):.3f}–{max(eg):.3f} (median {np.median(eg):.3f})") in s
    r1 = exb(1.0019) / exb(0.99997)
    r2 = exb(0.99997) / exb(0.99998)
    rs = sorted([r1, r2])
    eps = float(np.finfo(float).eps) / 1.0842e-19
    a_ = [np.log(x) / np.log(eps) for x in rs]
    assert (f"where the long double's exceeds that spread (2 of 4 rows): {rs[0]:.2g}, {rs[1]:.2g} — median {np.median(rs):.3g}, as "
            f"(ε_d/ε_ld)^a with a = {min(a_):.2f}–{max(a_):.2f} (median {np.median(a_):.2f}); a floor linear in the rounding (a = 1) gives "
            f"ε_d/ε_ld = {eps:.4g}; at the other 1.6e-09 (τ/σ = 1, K = 12) from the block: 0.0e+00, 1.6e-09 (τ/σ = 1, K = 12) from the "
            "optimum: 0.0e+00.") in s
    with localcontext() as ctx:
        ctx.prec = 60
        d2 = float((Decimal(float(2 * np.float64(WM.SIG) ** 2)) - 2 * Decimal(WM.SIG) ** 2) / (2 * Decimal(WM.SIG) ** 2))
    assert abs(abs(d2) - 4.66e-17) < 1e-19                              # the referee's 4.66·10⁻¹⁷, 2σ² rounded to double
    gl = [3e-13, 4e-13, 1e-16, 2e-16]
    rl = [x / g for x, g in zip([6e-13, 3e-13, 6e-13, 6e-13], gl)]
    re_ = [x / g for x, g in zip([7e-13, 2e-13, 7e-13, 7e-13], gl)]
    r3 = [x / y for x, y in zip([2e-13, 1e-13, 2e-13, 2e-13], [6e-13, 3e-13, 6e-13, 6e-13])]
    assert ("‖P_Je₀‖²/(2‖r‖²) at its end (`b4_floors`, exactly): with e₀ at the end 2.0e-13 … 7.0e-13; with e₀ at its last anchor, the "
            "rounding it froze, 3.0e-13 … 6.0e-13; its gap — the decrease a Gauss–Newton step predicts there — 1.0e-16 … 4.0e-13: the last "
            f"anchor's prediction {min(rl):.2f}–{max(rl):.2f} times the gap, the end's {min(re_):.2f}–{max(re_):.2f}; 2σ² rounded to double "
            f"(relative {d2:.3g}), the same at every anchor, alone predicts {min(r3):.2f}–{max(r3):.2f} of the last anchor's. Its last anchor: "
            "3 stall, 1 period, 1–28 iterations before the end (5.0e-10 … 4.5e-09 from it at the flat stops); per solve 10–40 re-anchors at "
            "stalls and 1–4 by period; without re-anchoring at stalls, from the blocks, 10–75 by period.") in s
    assert ("with its parameters rounded to double, over the best end, minus 1 (a scale, not the double parameters' floor: B5; the anchored double's "
            f"excesses in brackets): 1.6e-09 (τ/σ = 1, K = 12): {WM.b4_excess(rows[2]['anchored ld']['rms_d'], rows[2]['long double']['rms']):.1e} "
            f"(0.0e+00, 0.0e+00); 2.5e-11 (τ/σ = 0.25, K = 6): 1.5e-12 ({ex4(0.99996)}, {ex4(0.99996)}).") in s
    assert "The anchored double's ends from the block and from the stored optimum agree within 0.0e+00; its stops: flat 3, stall 1, after 3000 … 9000 iterations" in s
    gain = [WM.b4_excess(r_['anchored ld']['rms'], r_['anchored']['rms']) for r_ in rows]
    assert f"anchored solve from the anchored double's end gains on it, minus 1, {min(gain):+.1e} … {max(gain):+.1e}." in s
    assert f"The stored optimum over the best end, minus 1: 1.6e-09 (τ/σ = 1, K = 12): {1.0 / 0.999996 - 1:+.1e}; 2.5e-11 (τ/σ = 0.25, K = 6): {1.00003 / 0.9999599 - 1:+.1e}." in s
    assert f"after 100000, 200000, 300000 iterations, over the best end, minus 1: 2.5e-11 (τ/σ = 0.25, K = 6) from the block: {0.999961 / 0.9999599 - 1:+.1e}." in s
    best = WM.b4_best(rows)
    assert {c_: (b['end'], b['start'], b['err']) for c_, b in best.items()} == {(0.25, 6): ('anchored ld', 'block', 0.9999599 * 2.5e-11),
                                                                                (1.0, 12): ('long double', 'block', 0.999996 * 1.6e-9)}
    assert best[(0.25, 6)]['rms'] == Decimal(rows[0]['anchored ld']['rms'])
    assert WM.b4_lines(None) == [] and WM.b4_lines([]) == []


def test_b3_restart_judges_the_seed_from_the_start_with_the_paths_step(monkeypatch):
    """Referee round 77 (BR02, BR03): `b3_restart` judges every end by `seed_rule` from the start's amplitudes — a seed kernel that grows
    by the end is still judged — with the path's polish step: a 2σ path accepts a seed's move between 1σ and 2σ, a 1σ path rejects it."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    w0, c0 = np.array([0.7, 0.005]), np.array([0.0, 0.3])
    end = c0 + np.array([0.0, 1.5 * WM.SIG])
    t = RT.ascent_target(K, f - K.phi(c0) @ w0, c0[1], 1.0)[0]
    assert abs(end[1] - t) / WM.SIG > 2

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        return np.array([0.7, 0.3]), end.copy(), dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    st = dict(w=w0, c=c0, s=np.ones(2), max_nfev=5000, gtol=1e-15)
    for step, ok in ((2.0, True), (1.0, False)):
        q = WM.b3_restart(1.0, 'x', K, f, st, 2, 1e-3, 2.0, 1, step, 5000, 1e-15, 0.0)
        assert all(np.isclose(q[k]['move'], 1.5) and q[k]['accepted'] is ok for k in WM.B3_KEYS)


def test_b3_job_on_a_variant_without_the_gradient_test(monkeypatch):
    """Referee round 77 (BR07): on the variant without a gradient test B3's restarts have none — the path's budget and B3_BIG, then the
    same two again."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.5) ** 2))
    v = WM.VARIANTS_B[5]
    assert v[5] is None
    start = dict(w=np.array([0.7]), c=np.array([0.0]), s=np.array([1.0]), max_nfev=v[3], gtol=None)
    monkeypatch.setattr(WM, 'path_polishes', lambda r_, var, stop_after=None: (None, K, f, [start]))
    calls = []

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append((max_nfev, gtol))
        return w_, c_, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    q = WM.b3_job((0.5, v[0], dict(n0=0, k=1, e=2e-3, e_opt=1e-3, n_polish=1, first=dict(status=0, nfev=v[3]))))
    assert calls == [(v[3], None), (WM.B3_BIG, None), (v[3], None), (WM.B3_BIG, None)] and q['settings'] == dict(nfev=v[3], gtol=None)


def test_path_polishes_passes_the_settings_to_the_real_solve(monkeypatch):
    """Referee round 77 (PP04): `path_polishes`' recorder passes each call's own budget and gradient test to the solve it wraps."""
    from src.core import rich_transit as RT
    got = []

    def solve(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        got.append((max_nfev, gtol))
        return w_, c_, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', solve)

    def fake_path(L, K, f, c0, rho, **kw):
        for i, g in enumerate((None, 1e-12)):
            RT.polish_block(K, f, np.array([0.5]), np.array([0.0]), np.array([1.0]), max_nfev=7 + i, info=True, gtol=g)
        return 'full'
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    R, K, f, calls = WM.path_polishes(1.0, 'base')
    assert R == 'full' and got == [(7, None), (8, 1e-12)] and [(c['max_nfev'], c['gtol']) for c in calls] == got
    assert RT.polish_block is solve


def test_lm_fit_flat_stop(monkeypatch):
    """Referee rounds 77, 81: the anchored solve stops ('flat') once the exact rms over a period of anchor_every accepted steps falls by
    less than B4_FLAT of itself — here a period of 20 and B4_FLAT = ½ —; the period's accepted steps are counted across the re-anchors
    at stalls (the exact rms evaluated once at the start and once per period); without stall re-anchoring the solve stops at its first
    stall, earlier."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    ss = WM.sigma_star(0.5)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    co, ao, eo = WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(8, WM.SIG, ss)[0])
    monkeypatch.setattr(WM, 'B4_FLAT', 0.5)
    p, it, acc, st = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=20)
    assert st == 'flat' and acc == 20 and it < 100
    monkeypatch.setattr(WM, 'B4_FLAT', -1.0)                            # never flat
    real, n_ex, stalls = WM.exact_rms, [], []
    monkeypatch.setattr(WM, 'exact_rms', lambda *a: n_ex.append(1) or real(*a))
    p, it, acc, st = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=30, max_iter=1500)
    assert st in ('cap', 'stall') and acc >= 60 and len(n_ex) == 1 + acc // 30
    pn, itn, accn, stn = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=30, max_iter=1500, stall_reanchor=False)
    assert stn == 'stall' and itn < it and accn < acc
    pd_, itd, accd, std_ = WM.lm_fit(f, ao, co, np.float64, max_iter=1500)
    assert std_ == 'stall' and accd <= accn                              # direct subtraction stalls first: its rounding is not frozen


def test_the_exact_residual_parts_and_the_floors():
    """Referee rounds 81, 84 (D3cb-1; D3ce-4, D3ce-5): `exact_resid` — every double and long double taken exactly, exp in 60-digit
    decimal — agrees with 90 digits to the double's rounding, its decimal rms rounding to its float, the long double entering as its exact
    hi + lo; at B4's 2.55·10⁻¹¹ optimum (τ/σ = 1/4, K = 6) the long-double rms by direct subtraction (`lm_rms`) is off by 10⁻¹³ … 10⁻⁹
    relative, the residual by direct subtraction in double (`direct_double`) within 10⁻¹⁵ of max|f|; `exact_parts`' residual and
    Jacobian are the exact residual and the double Jacobian; `gn_share` is vᵀP_Jv/(2‖r‖²) (a planted v: its part in the Jacobian's range
    and a part orthogonal to it); `b4_floors` with the end as its anchor gives the end's floor twice, its 2σ² part the first order
    −Σ_kw_kφ_k(u_k²/(2σ²))δ (δ the relative rounding of 2σ²); `b4_gap` at the stored optimum predicts its excess over the best end."""
    from decimal import Decimal, localcontext
    from fractions import Fraction
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if j[0] == 0.25 and j[2] == 6 and j[3] == 'optimum')
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))
    p = np.r_[j[4], j[5]]
    r60, e60, d60 = WM.exact_resid(f, p, 6, dec=True)
    r90, e90 = WM.exact_resid(f, p, 6, prec=90)
    assert np.array_equal(r60, r90) and e60 == e90 and float(d60) == e60 and abs(d60 / Decimal(e60) - 1) < Decimal(10) ** -15
    rd = f - np.exp(-(WM.X_B[:, None] - p[None, 6:]) ** 2 / (2 * WM.SIG ** 2)) @ p[:6]
    assert np.array_equal(WM.direct_double(f, p, 6), rd)
    assert np.abs(rd - r60).max() < 1e-15 * np.abs(f).max() and np.abs(rd - r60).max() > 0
    rel = abs(WM.lm_rms(f, p, 6) / e60 - 1)
    assert 1e-13 < rel < 1e-9
    v = np.longdouble(1) / np.longdouble(3)
    hi = float(v)
    assert Fraction(WM._dec(v)) == Fraction(hi) + Fraction(float(v - np.longdouble(hi))) and WM._dec(v) != Decimal(hi)   # exact, in any context
    rx, Jx = WM.exact_parts(f, p, 6)
    u = WM.X_B[:, None] - p[None, 6:]
    P = np.exp(-u ** 2 / (2 * WM.SIG ** 2))
    J = -np.c_[P, P * u / WM.SIG ** 2 * p[None, :6]]
    assert np.array_equal(np.array([float(x) for x in rx]), r60) and np.allclose(np.array([[float(x) for x in row] for row in Jx]), J, rtol=1e-13, atol=1e-300)
    a = J[:, 0] * 1e-15
    b = np.random.default_rng(3).standard_normal(len(WM.X_B))
    b -= J @ np.linalg.lstsq(J, b, rcond=None)[0]                       # orthogonal to the range
    b *= 1e-16 / np.abs(b).max()
    with localcontext() as ctx:
        ctx.prec = WM.GN_PREC
        vs = [Decimal(float(x)) for x in a + b]
    g = WM.gn_share(Jx, rx, [vs, rx])
    assert np.isclose(g[0], (a @ a) / (2 * (r60 @ r60)), rtol=1e-6, atol=0) and g[1] > 0
    fl = WM.b4_floors(f, p, 6, p, rd)
    assert fl['last'] == fl['end'] and fl['gap'] == g[1] == WM.b4_gap(f, p, 6)
    with localcontext() as ctx:
        ctx.prec = 60
        dl = float((Decimal(float(2 * np.float64(WM.SIG) ** 2)) - 2 * Decimal(WM.SIG) ** 2) / (2 * Decimal(WM.SIG) ** 2))
    es = (P * (u ** 2 / (2 * WM.SIG ** 2)) * p[None, :6]).sum(axis=1) * dl  # φ(s2(1 + δ)) − φ(s2) ≈ φ·u²/s2·δ; r falls by Σwφ·u²/s2·δ
    pe = J @ np.linalg.lstsq(J, es, rcond=None)[0]
    assert np.isclose(fl['s2'], (pe @ pe) / (2 * (r60 @ r60)), rtol=1e-6, atol=0)        # atol=0: the floors are ~1e-13
    pa = p.copy()                                                       # an anchor apart from the end: its 2σ² part at the anchor (BF4)
    pa[6:] += 0.02
    pa[:6] *= 1.1
    fla = WM.b4_floors(f, p, 6, pa, WM.direct_double(f, pa, 6))
    ua = WM.X_B[:, None] - pa[None, 6:]
    esa = (np.exp(-ua ** 2 / (2 * WM.SIG ** 2)) * (ua ** 2 / (2 * WM.SIG ** 2)) * pa[None, :6]).sum(axis=1) * dl
    pea = J @ np.linalg.lstsq(J, esa, rcond=None)[0]
    assert np.isclose(fla['s2'], (pea @ pea) / (2 * (r60 @ r60)), rtol=1e-6, atol=0) and not np.isclose(fla['s2'], fl['s2'], rtol=1e-3, atol=0)
    assert fla['end'] == fl['end'] and fla['gap'] == fl['gap'] and fla['last'] != fl['last']
    best = min(min(x[k_]['ratio'] for k_, _ in WM.B4_ENDS) for x in raw['B4'] if (x['r'], x['k']) == (0.25, 6))
    x = next(x for x in raw['B4'] if (x['r'], x['k'], x['start']) == (0.25, 6, 'optimum'))
    assert abs(fl['gap'] / (x['e_start'] / best - 1) - 1) < 2e-3          # the stored optimum's excess, 3.1·10⁻⁵, to first order


def test_b4_job_on_three_cells_against_the_raw():
    """Referee rounds 81, 84 (D3cb-2; D3ce-2): the real `b4_job` — no fakes — on B4's two fastest cells (τ/σ = 1/4, K = 5, from the stored
    optimum and from the path's block) and on one whose anchored solve stops flat after a periodic re-anchor (τ/σ = 3/4, K = 11, from the
    stored optimum) reproduces the raw's row in every field: each solve's exact rms and errors, Gauss–Newton gap, iterations, accepted
    steps, stop and parameters, the anchored double's floors, anchors and marks, the 'nostall' solve's anchors, the long double's anchored
    end rounded to double, TRF's status and evaluations."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    for cell in ((0.25, 5, 'optimum'), (0.25, 5, 'block'), (0.75, 11, 'optimum')):
        j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if (j[0], j[2], j[3]) == cell)
        out = WM.b4_job(j)
        x = next(q for q in raw['B4'] if (q['r'], q['k'], q['start']) == cell)
        assert set(out) == set(x) and all(out[k_] == x[k_] for k_ in out if k_ not in ('secs',) + tuple(k for k, _ in WM.B4_ENDS)), cell
        for k_, _ in WM.B4_ENDS:
            assert out[k_] == x[k_], (cell, k_, [f_ for f_ in out[k_] if out[k_][f_] != x[k_].get(f_)])
    assert x['anchored']['stop'] == 'flat' and x['anchored']['anchors']['period'] >= 1


def test_lm_fit_long_descent_from_the_paths_block():
    """Referee round 77 (LM15): from B4's farthest start — the path's block at τ/σ = 1/4, K = 5, ten times the optimum's error — the
    double solve reaches the optimum within 10⁻⁶ in about 8000 iterations, its damping relaxing after its rejections."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if j[0] == 0.25 and j[2] == 5 and j[3] == 'block')
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    assert WM.lm_rms(f, np.r_[j[4], j[5]], 5) / rf / j[6] > 10
    p, it, acc, st = WM.lm_fit(f, j[4], j[5], np.float64, max_iter=20_000)
    assert st == 'stall' and 5000 < it < 20_000 and acc > 5000 and abs(WM.lm_rms(f, p, 5) / rf / j[6] - 1) < 1e-6


def test_part_b4_and_main_only_b4(monkeypatch, tmp_path):
    """Referee round 77 (PB01, MN01): `part_B4` maps `b4_job` over every job; `main --only B4` passes the raw's B and B3 to it and stores
    its result."""
    import multiprocessing
    import exp_cohesion_leak as E

    class FakePool:
        def __init__(self, n):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            return [fn(j) for j in jobs]
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(WM, 'b4_jobs', lambda B, B3: [('j', i) for i in range(5)])
    monkeypatch.setattr(WM, 'b4_job', lambda j: dict(job=j))
    assert WM.part_B4('B', 'B3') == [dict(job=('j', i)) for i in range(5)]
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), B=[dict(r=1.0, var='base')], B3=[dict(r=1.0, var='base', k=None)])
    got, saved = [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'part_B4', lambda B, B3: (got.append((B, B3)), ['b4'])[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'B4'])
    WM.main()
    assert len(got) == 1 and list(got[0][0]) == [dict(r=1.0, var='base')] and list(got[0][1]) == [dict(r=1.0, var='base', k=None)]
    assert saved[-1]['B4'] == ['b4'] and list(saved[-1]['B']) == [dict(r=1.0, var='base')]


def test_lm_fit_kappa_and_the_anchor_log():
    """Referee round 84 (D3ce-1, D3ce-4): `kappa_blend` at κ = 1 is the double's residual bit for bit (r_d − r_ld exact by Sterbenz's
    lemma), at κ = 0 the long double's rounded to double, at κ = ½ the midpoint; `lm_fit` with kappa = 1 is the double solve, direct
    and anchored (the same ends, iterations, accepted steps and stops), and kappa needs dtype double; anchor_log counts the anchors by
    kind and keeps the last one's parameters, iteration, kind and the r₀ it froze — the double's direct subtraction there."""
    import pytest
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    ss = WM.sigma_star(0.5)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    co, ao, eo = WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(8, WM.SIG, ss)[0])
    p = np.r_[ao, co]
    rd = WM.direct_double(f, p, 8)
    LD = np.longdouble
    q = p.astype(LD)
    rl = np.asarray(f, float).astype(LD) - np.exp(-(WM.X_B.astype(LD)[:, None] - q[None, 8:]) ** 2 / (2 * LD(WM.SIG) ** 2)) @ q[:8]
    assert np.array_equal(WM.kappa_blend(f, p, 8, rd, 1.0), rd) and np.array_equal(WM.kappa_blend(f, p, 8, rd, 0.0), rl.astype(float))
    h = WM.kappa_blend(f, p, 8, rd, 0.5)
    assert np.array_equal(h, (rl + (rd.astype(LD) - rl) / 2).astype(float)) and not np.array_equal(h, rd) and not np.array_equal(h, rl.astype(float))
    for kw in (dict(max_iter=300), dict(anchored=True, max_iter=300)):
        a, b = WM.lm_fit(f, ao, co, np.float64, **kw), WM.lm_fit(f, ao, co, np.float64, kappa=1.0, **kw)
        assert np.array_equal(a[0], b[0]) and a[1:] == b[1:], kw
        c = WM.lm_fit(f, ao, co, np.float64, kappa=0.0, **kw)
        assert not np.array_equal(a[0], c[0]), kw
    with pytest.raises(ValueError):
        WM.lm_fit(f, ao, co, np.longdouble, kappa=1.0)
    log = {}
    p2, it, acc, st = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=20, max_iter=2000, anchor_log=log)
    assert log['start'] == 1 and log.get('stall', 0) + log.get('period', 0) >= 1 and log['kind'] in ('stall', 'period') and 0 < log['it'] <= it
    assert np.array_equal(log['r0'], WM.direct_double(f, log['p'], 8))
    a = WM.anchor_info(log, p2, it)
    assert a == dict(stall=log.get('stall', 0), period=log.get('period', 0), kind=log['kind'], before=it - log['it'],
                     dist=float(np.abs(p2 - log['p']).max()))
    sig = __import__('inspect').signature(WM.lm_fit).parameters
    assert sig['kappa'].default is None and sig['anchor_log'].default is None and WM.GN_PREC == 80


def test_b5_jobs_job_and_lines(monkeypatch):
    """Referee round 84 (D3ce-1): B5's jobs are B4's two starts at its blocks and the starts between them (centres c_opt + θ(c_block −
    c_opt), least-squares amplitudes) for the direct solve, and the stored optimum's for the anchored one; `b5_job` runs `lm_fit` in
    double at every κ of its mode, recording each end's exact rms, ratio, iterations, accepted steps, stop and parameters; `b5_lines`
    prints each end's excess over B4's best end, the slopes of log(excess) against log κ, the geometric means over the starts, the
    anchored rows and whether the ends at κ = 1 are B4's."""
    from decimal import Decimal
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    o = lambda k, e_opt, c: dict(k=k, e=2 * e_opt, e_opt=e_opt, c=np.array(c), w=np.full(len(c), 0.4), c_opt=np.array(c) + 0.01)
    B = [dict(r=0.5, var='v1', blocks=[o(2, 3e-10, [-0.15, 0.15])])]
    B3 = [dict(r=0.5, var='v1', k=2, e_opt=3e-10)]
    monkeypatch.setattr(WM, 'B5_CELLS', ((0.5, 2), (2.0, 7)))
    jobs = WM.b5_jobs(B, B3)
    assert [(j[0], j[1], j[2], j[3], j[4]) for j in jobs] == [(0.5, 'v1', 2, th, 'direct') for th in WM.B5_THETA] + [(0.5, 'v1', 2, 0.0, 'anchored')]
    b4 = {j[3]: j for j in WM.b4_jobs(B, B3)}
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.5) ** 2))
    assert np.array_equal(jobs[0][6], b4['optimum'][5]) and np.array_equal(jobs[0][5], b4['optimum'][4])
    assert np.array_equal(jobs[4][6], b4['block'][5]) and np.array_equal(jobs[4][5], b4['block'][4]) and np.array_equal(jobs[5][6], b4['optimum'][5])
    cm = b4['optimum'][5] + 0.5 * (b4['block'][5] - b4['optimum'][5])
    cq = b4['optimum'][5] + 0.25 * (b4['block'][5] - b4['optimum'][5])  # θ = 1/4 from the optimum's end (B53: not from the block's)
    assert jobs[1][3] == 0.25 and np.allclose(jobs[1][6], cq, rtol=0, atol=1e-15) and not np.allclose(jobs[3][6], cq, rtol=0, atol=1e-6)
    assert np.allclose(jobs[2][6], cm, rtol=0, atol=1e-15) and np.allclose(jobs[2][5], np.linalg.lstsq(K.phi(cm), f, rcond=None)[0])
    assert WM.B5_THETA == (0.0, 0.25, 0.5, 0.75, 1.0) and WM.B5_KAPPA == (0.0, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0, 2.0, 4.0)
    assert WM.B5_KAPPA_ANC == (0.0, 1 / 4, 1 / 2, 1.0, 2.0, 4.0)
    calls = []

    def fake(f_, w_, c_, dtype, anchored=False, kappa=None, **kw):
        calls.append((dtype, anchored, kappa, kw))
        return np.r_[w_, c_ + 1e-4 * (1 + kappa)].astype(dtype), 11, 7, 'stall'
    monkeypatch.setattr(WM, 'lm_fit', fake)
    for j in (jobs[1], jobs[5]):
        calls.clear()
        out = WM.b5_job(j)
        ks = WM.B5_KAPPA if j[4] == 'direct' else WM.B5_KAPPA_ANC
        assert [c_[:3] for c_ in calls] == [(np.float64, j[4] == 'anchored', k_) for k_ in ks] and all(c_[3] == {} for c_ in calls)
        rf = float(np.sqrt(np.mean(f ** 2)))
        for y, k_ in zip(out['rows'], ks):
            pe = np.r_[j[5], j[6] + 1e-4 * (1 + k_)]
            _, e, d = WM.exact_resid(f, pe, 2, dec=True)
            assert y == dict(kappa=k_, rms=str(d), ratio=e / rf / 3e-10, it=11, acc=7, stop='stall', p=pe.tolist())
        assert out['e_start'] == WM.exact_rms(f, np.r_[j[5], j[6]], 2) / rf / 3e-10 and (out['theta'], out['mode']) == (j[3], j[4])
    bst = Decimal('2.0e-10')
    mk = lambda x: str(bst * (1 + Decimal(repr(x))))
    B4 = [dict(r=0.5, k=2, start=s_, e_opt=3e-10, **{k_: dict(rms=str(bst * 2), ratio=2.0, it=1, acc=1, stop='stall', p=[0.0]) for k_, _ in WM.B4_ENDS})
          for s_ in ('block', 'optimum')]
    B4[1]['anchored ld'] = dict(rms=str(bst), ratio=1.0, rms_d=mk(1e-12), it=1, acc=1, stop='stall', p=[0.0])
    B4[1]['double'] = dict(rms=mk(3e-6), ratio=1.0, it=91, acc=37, stop='stall', p=[0.5, 0.1])
    row = lambda k_, x, it=9, p=(0.0,): dict(kappa=k_, rms=mk(x), ratio=1.0, it=it, acc=5, stop='stall', p=list(p))
    D = lambda th, c, a: dict(r=0.5, k=2, theta=th, mode='direct', e_start=1.5, rows=[row(k_, c * k_ ** a if k_ else c * 1e-3) for k_ in WM.B5_KAPPA])
    B5 = [D(0.0, 3e-6, 2.0), D(0.5, 1e-5, 1.5), D(1.0, 4e-6, 1.0),
          dict(r=0.5, k=2, theta=0.0, mode='anchored', e_start=1.0, rows=[row(k_, 1e-13 * k_ ** 2 + 2e-12, 100) for k_ in WM.B5_KAPPA_ANC])]
    B5[0]['rows'][5] = dict(row(1.0, 3e-6, 91, (0.5, 0.1)), acc=37)
    B5[3]['rows'][2] = dict(B5[3]['rows'][2], stop='flat')            # κ = 1/2: one flat stop to count
    B4[0]['double'] = dict(rms=mk(4e-6), ratio=1.0, it=1, acc=1, stop='stall', p=[0.0])   # the block's: the same rms, not the same end (L54)
    s = '\n'.join(WM.b5_lines(B5, B4))
    assert "### B5. The floors' order: the double's rounding scaled (referee round 84)" in s and '| κ = 0 | κ = 1/16 | κ = 1/8 |' in s
    ex = lambda x: WM.b4_excess(mk(x), bst)
    assert f"| 0.5 | 2 | 0 (1.5) | {ex(3e-9):.1e} | {ex(3e-6 / 256):.1e} | " in s and f"| 2.00 | 2.00 |" in s and "| 1.50 | 1.50 |" in s
    assert "| 1.00 | 1.00 |" in s
    g = [float(np.exp(np.mean(np.log([ex(c * k_ ** a if k_ else c * 1e-3) for c, a in ((3e-6, 2.0), (1e-5, 1.5), (4e-6, 1.0))]))))
         for k_ in WM.B5_KAPPA]
    assert f"τ/σ = 0.5, K = 2 (3 starts): 0: {g[0]:.1e}, 1/16: {g[1]:.1e}" in s and "— slopes 1.50 on [1/16, 1], 1.50 on [1/4, 4]" in s
    an = [ex(1e-13 * k_ ** 2 + 2e-12) for k_ in WM.B5_KAPPA_ANC]
    assert f"| 0.5 | 2 | {an[0]:.1e} (100, stall) | {an[1]:.1e} (100, stall) |" in s
    sl = WM.log_slope([k_ for k_ in WM.B5_KAPPA_ANC if k_ >= 0.5], [e for k_, e in zip(WM.B5_KAPPA_ANC, an) if k_ >= 0.5])
    assert f"| {sl:.2f} | " + ', '.join(f"{e / ex(1e-12):.3g}" for e in an) + " |" in s and f"{an[2]:.1e} (100, flat)" in s
    assert "| slope on [1/2, 4] | over the best end rounded to double, κ = 0, 1/4, 1/2, 1, 2, 4 |" in s
    es3 = [[ex(c * k_ ** a if k_ else c * 1e-3) for k_ in WM.B5_KAPPA] for c, a in ((3e-6, 2.0), (1e-5, 1.5), (4e-6, 1.0))]
    sp = [max(col) / min(col) for col in zip(*es3)]
    i_ = int(np.argmax(sp))
    fr = lambda k_: '0' if k_ == 0 else (f"1/{round(1 / k_)}" if k_ < 1 else f"{k_:g}")
    assert (f"Per start (3): the slopes 1.00–2.00 on [1/16, 1] (median 1.50), 1.00–2.00 on [1/4, 4] (median 1.50); the starts' largest "
            f"over least at one κ at most {sp[i_]:.3g} (τ/σ = 0.5, K = 2, κ = {fr(WM.B5_KAPPA[i_])}).") in s
    assert (f"The anchored solves' slopes on [1/2, 4] {sl:.2f}–{sl:.2f}; their flat stops by κ: 0: 0 of 1, 1/4: 0 of 1, 1/2: 1 of 1, 1: 0 of 1, "
            "2: 0 of 1, 4: 0 of 1.") in s
    assert "At κ = 1 the ends are B4's (the double's from the stored optimum and the block, the anchored double's from the stored optimum) in 1 of 3" in s
    assert WM.log_slope([1, 2, 4], [1, 4, 16]) == 2.0 or abs(WM.log_slope([1, 2, 4], [1, 4, 16]) - 2) < 1e-12
    assert WM.log_slope([0.0, 1.0], [1.0, 2.0]) is None and WM.b5_lines(None, B4) == [] and WM.b5_lines(B5, None) == []


def test_b5_job_against_the_raw():
    """Referee round 84 (D3ce-1): the real `b5_job` on B5's fastest job (τ/σ = 1/4, K = 6, the direct solve from the stored optimum)
    reproduces the raw's row in every field, and its end at κ = 1 is B4's double end from that start."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    j = next(j for j in WM.b5_jobs(raw['B'], raw['B3']) if (j[0], j[2], j[3], j[4]) == (0.25, 6, 0.0, 'direct'))
    out = WM.b5_job(j)
    x = next(x for x in raw['B5'] if (x['r'], x['k'], x['theta'], x['mode']) == (0.25, 6, 0.0, 'direct'))
    assert out['e_start'] == x['e_start'] and [dict(y) for y in out['rows']] == [dict(y) for y in x['rows']]
    y = next(y for y in out['rows'] if y['kappa'] == 1.0)
    z = next(q for q in raw['B4'] if (q['r'], q['k'], q['start']) == (0.25, 6, 'optimum'))['double']
    assert (y['rms'], y['it'], y['acc'], y['stop'], list(y['p'])) == (z['rms'], z['it'], z['acc'], z['stop'], list(z['p']))


def test_frame_series_at_k1_in_the_half_squeeze_frame():
    """Referee round 84 (D3ce-6): at K = 1 in the half-squeeze frame (ρ = t, v′ = 2) D_n = (n − 1)!!(1 − (−1)^{n/2}) vanishes for
    n ≡ 0 mod 4, and `frame_series` sums past those zeros to S(1) = 2((1 − t²)^{−1/2} − (1 + t²)^{−1/2}) within its 10⁻⁴⁰."""
    from decimal import Decimal, localcontext
    for r_ in (0.25, 1.0, 2.0):
        with localcontext() as ctx:
            ctx.prec = WM.A5_PREC
            st = (1 + Decimal(r_) ** 2).sqrt()
            t = (st - 1) / (st + 1)
            S, terms, low = WM.frame_series(1, t, Decimal(2))
            ref = 2 * (1 / (1 - t * t).sqrt() - 1 / (1 + t * t).sqrt())
            assert abs(S / ref - 1) < Decimal(10) ** -39 and terms[1] == 0 and terms[0] > 0 and len(terms) > 12, r_


def test_part_a5_against_the_raw(monkeypatch):
    """Referee round 84 (D3ce-2): A5 against the raw in every field — the frames' table and the tail-ratio checks (`part_A5` without its
    series), the exact series at every τ/σ (`a5_series`), and the Mercer frame's rule further out at four cells (`a5_far_row`, τ/σ = 1/4
    in 4K + 60 digits); `a5_far`'s cells, precisions and its check of each τ/σ's last K with A5_FAR_CHECK more digits."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    A5 = WM.part_A5(raw['A'], series=False, far=False)
    eq = lambda a, b: (a.keys() == b.keys() and all(eq(a[k], b[k]) for k in a)) if isinstance(a, dict) else (
        len(a) == len(b) and all(eq(x, y) for x, y in zip(a, b))) if isinstance(a, (list, tuple, np.ndarray)) else a == b
    assert eq(A5['rows'], list(raw['A5']['rows'])) and eq(A5['check'], list(raw['A5']['check'])) and A5['series'] == [] and A5['far'] == []
    for x in raw['A5']['series']:
        assert eq(WM.a5_series(x['r']), x), x['r']
    for r_, K in ((2.0, 60), (1.0, 120), (0.5, 80), (0.25, 40)):
        y = next(y for y in raw['A5']['far'] if (y['r'], y['K']) == (r_, K))
        z = WM.a5_far_row(r_, K, y['prec'])
        assert (z['ratio'], z['low'], z['prec']) == (y['ratio'], y['low'], WM.A5_FAR_PREC if r_ != 0.25 else 4 * K + 60), (r_, K)
    cells = []
    monkeypatch.setattr(WM, 'a5_far_row', lambda r_, K, prec: cells.append((r_, K, prec)) or dict(r=r_, K=K, prec=prec, ratio=K + prec, low=0.0))
    far = WM.a5_far()
    main = [(r_, K, WM.A5_FAR_PREC) for r_ in WM.A5_FAR[1] for K in WM.A5_FAR[0]] + [(0.25, K, 4 * K + 60) for K in WM.A5_FAR_Q]
    last = [(r_, max(K for r2, K, _ in main if r2 == r_)) for r_ in (0.5, 1.0, 2.0, 0.25)]
    assert cells == main + [(r_, K, p + WM.A5_FAR_CHECK) for r_, K in last for r2, K2, p in main if (r2, K2) == (r_, K)]
    assert [y['check'] for y in far if 'check' in y] == [K + p + WM.A5_FAR_CHECK for r_, K in last for r2, K2, p in main if (r2, K2) == (r_, K)]
    assert sum('check' in y for y in far) == 4 and WM.A5_FAR_Q == (40, 60, 80, 100, 120, 160) and WM.A5_FAR_CHECK == 40


def test_a5_far_lines():
    """Referee round 84 (D3ce-3, D3ce-7): the far line prints each τ/σ's ratios, the successive differences' ratios against K⁻²'s, the limit
    by L + c/K² + d/K³ through the last three, the conjecture's (b² + b + 2)/(b(b − 1)) and 1/t, the precision check, and whether the
    limits lie below 1/t."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    A5 = dict(raw['A5'])
    Ks = (60, 80, 100, 120)
    lim = lambda r_: (lambda b: (b * b + b + 2) / (b * (b - 1)))(np.sqrt(1 + 2 * r_ ** 2))
    A5['far'] = [dict(r=r_, K=K, prec=140, ratio=lim(r_) - 0.5 / K ** 2, low=1e-50) for r_ in (2.0, 1.0) for K in Ks]
    A5['far'][-1]['check'] = A5['far'][-1]['ratio'] * (1 + 3e-12)
    s = '\n'.join(WM.a5_series_lines(A5))
    t2 = next(x['t'] for x in A5['series'] if x['r'] == 2.0)
    q2 = [(Ks[i + 1] ** -2.0 - Ks[i + 2] ** -2.0) / (Ks[i] ** -2.0 - Ks[i + 1] ** -2.0) for i in range(2)]
    assert (f"τ/σ = 2 (140 digits): 60: {lim(2.0) - 0.5 / 3600:.6f}, 80: " in s and f"the successive differences' ratios {q2[0]:.3f}, {q2[1]:.3f} "
            f"(K⁻²'s {q2[0]:.3f}, {q2[1]:.3f}); the limit by L + c/K² + d/K³ through the last three {lim(2.0):.6f}, the formula's "
            f"(b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{{1/2}}, {lim(2.0):.6f} (A6), 1/t {1 / t2:.4f}") in s
    assert "τ/σ = 1 (140 digits)" in s and "; the last K with 40 more digits 3.0e-12 apart" in s and s.index("τ/σ = 2 (") < s.index("τ/σ = 1 (")
    A5['far'][-1]['check'] = A5['far'][-1]['ratio']
    assert "; the last K with 40 more digits the same in double" in '\n'.join(WM.a5_series_lines(A5))
    assert "Its limit lies below 1/t at every τ/σ, within " in s and "the low defects at most 1.0e-50" in s
    A5['far'] = [dict(y, ratio=y['ratio'] + 1) for y in A5['far']]
    assert "NOT below 1/t at every τ/σ" in '\n'.join(WM.a5_series_lines(A5))
    assert "it bounds that rounding, not the tail's" in s


def test_a6_scaled_rules_and_the_conjectured_rate(monkeypatch):
    """Referee round 84 (D3ce-3): `a6_error` is `scaled_rule`'s L² error (against the trapezoid at small K, every rule of `a6_scalings`);
    the conjecture `rate_conj` gives λ at moment matching, t at v_hs (both arguments), λ/(2 − λ) at the Fock scaling, q at the Mercer-zero
    scaling and b(b − 1)/(b² + b + 2) for the Mercer frame's rule, at least t everywhere; the real `a6_job` at τ/σ = 2 (the Mercer frame's
    rule) reproduces the raw's row; `a6_lines` fits L + c/K² + d/K³ (at the cusp L + c/K + d/K²) and compares the Mercer frame's rule's
    exact squared errors with A5's series."""
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    for r_ in (0.5, 1.0, 2.0):
        ss = WM.sigma_star(r_)
        for K in (2, 4):
            with localcontext() as ctx:
                ctx.prec = 60
                for name, v in WM.a6_scalings(r_):
                    e_ex = float(WM.a6_error(r_, v, K, 60).sqrt())
                    e_tr = W.l2_error(*W.scaled_rule(K, WM.SIG, ss, float(v) * WM.SIG ** 2), WM.SIG, ss)
                    assert abs(e_ex / e_tr - 1) < 1e-6, (r_, K, name)
    for r_ in (0.25, 0.5, 1.0, 2.0):
        ss = np.sqrt(1 + r_ ** 2)
        t, lam, b = (ss - 1) / (ss + 1), r_ ** 2 / (2 + r_ ** 2), np.sqrt(1 + 2 * r_ ** 2)
        with localcontext() as ctx:
            ctx.prec = 40
            vs = {k_: float(v) for k_, v in WM.a6_scalings(r_)}
        assert [k_ for k_ in vs] == ['moment matching', 'Fock', 'Mercer-zero', 'Mercer frame'] + [f"{x:g}·v_hs" for x in WM.A6_VR]
        for k_, x in (('moment matching', lam), ('1·v_hs', t), ('Fock', lam / (2 - lam)), ('Mercer-zero', (b - 1) / (b + 1)),
                      ('Mercer frame', b * (b - 1) / (b * b + b + 2))):
            assert abs(WM.rate_conj(r_, vs[k_])[0] / x - 1) < 1e-12, (r_, k_)
        assert all(WM.rate_conj(r_, v)[0] >= t * (1 - 1e-12) for v in vs.values())
        vhs = r_ ** 2 / ss                                              # the tanh form: the mis-scaling adds to the squeeze r = ln(σ*/σ)
        assert all(abs(WM.rate_conj(r_, v)[0] - np.tanh((np.log(ss) + abs(np.log(v / vhs))) / 2)) < 1e-14 for v in vs.values())
        assert WM.rate_conj(r_, vs['0.5·v_hs'])[1] == 1 and WM.rate_conj(r_, vs['2·v_hs'])[1] == 2
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    assert WM.A6_K == (20, 40, 60, 80, 120, 160) and WM.A6_MARGIN == 60 and WM.A6_TAUS == (0.25, 0.5, 1.0, 2.0)
    assert [(x['r'], x['name']) for x in raw['A6']] == [(r_, k_) for r_ in WM.A6_TAUS for k_, _ in WM.a6_scalings(r_)]
    for x in raw['A6']:
        assert [y['K'] for y in x['rows']] == list(WM.A6_K), (x['r'], x['name'])
        assert [y['prec'] for y in x['rows']] == [int(np.ceil(2 * K * np.log10(1 / x['t']))) + 60 for K in WM.A6_K], (x['r'], x['name'])
    i = [k_ for k_, _ in WM.a6_scalings(2.0)].index('Mercer frame')
    full = WM.A6_K
    monkeypatch.setattr(WM, 'A6_K', full[:2])                            # the real job at the first two K: the raw's rows there
    out = WM.a6_job((2.0, i))
    monkeypatch.setattr(WM, 'A6_K', full)
    x = next(x for x in raw['A6'] if (x['r'], x['name']) == (2.0, 'Mercer frame'))
    assert {k_: v for k_, v in out.items() if k_ not in ('secs', 'rows')} == {k_: v for k_, v in x.items() if k_ not in ('secs', 'rows')}
    assert [dict(y) for y in out['rows']] == [dict(y) for y in x['rows'][:2]] and [y['K'] for y in out['rows']] == [20, 40]
    Ks = np.array(WM.A6_K, float)
    rows = lambda L, c, p: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=L + c * K ** p) for K in Ks]
    A6 = [dict(r=1.0, name='Mercer frame', v=1.0, v_hs=1.2, rate=0.25, arg=2, t=0.17, rows=rows(4.0, 3.0, -2.0)),
          dict(r=1.0, name='1·v_hs', v=1.0, v_hs=1.0, rate=1 / 5.8, arg=1, t=0.17, rows=rows(5.8, -2.0, -1.0)),
          dict(r=1.0, name='0.5·v_hs', v=1.0, v_hs=0.5, rate=1 / 4.0, arg=1, t=0.17, rows=[dict(K=20, prec=100, e0=None, e1=None, ratio=None)] * 4)]
    s = '\n'.join(WM.a6_lines(A6))
    row = next(x for x in s.split('\n') if x.startswith('| 1 | Mercer frame | 1.200000 | 4.000000 (vσ*²/τ²) | '))
    cells = row.rstrip(' |').split(' | ')
    y1, y0 = (4.0 + 3.0 / 160 ** 2) * 0.25 - 1, (4.0 + 3.0 / 120 ** 2) * 0.25 - 1
    assert cells[-3] == '4.000000' and abs(float(cells[-2])) < 1e-12 and cells[-1] == f"{y1:+.1e} (2.00)" and len(cells) == 9
    assert "| e_{K−1}/e_K·t(v) − 1 at K = 160 (its exponent from K = 120) |" in s and '|---|---|---|---|---|---|---|---|---|' in s
    assert "| 1 | 1·v_hs | 1.000000 | 5.800000 (both arguments, the cusp) | " in s and "| 1 | 0.5·v_hs | 0.500000 | 4.000000 | — | — | — | — | — |" in s
    d = np.diff([4.0 + 3.0 * K ** -2.0 for K in Ks])
    assert ', '.join(f"{d[i + 1] / d[i]:.3f}" for i in range(len(d) - 1)) in s
    assert "(as K⁻²: 0.185, 0.350, 0.714, 0.350; as K⁻¹: 0.333, 0.500, 1.000, 0.500)" in s
    assert "Off the cusp the limit lies within " in s and "at all 1 rules (the farthest: τ/σ = 1, Mercer frame); at v_hs within " in s
    rr_ = [4.0 + 3.0 * K ** -2.0 for K in Ks]
    L1 = float(np.linalg.solve(np.array([[1.0, K ** -2.0, K ** -3.0] for K in Ks[-3:]]), rr_[-3:])[0]) * 0.25 - 1
    L2 = float(np.linalg.solve(np.array([[1.0, K ** -1.0, K ** -2.0] for K in Ks[-3:]]), rr_[-3:])[0]) * 0.25 - 1
    ex_ = -float(np.log((rr_[-1] * 0.25 - 1) / (rr_[-2] * 0.25 - 1)) / np.log(Ks[-1] / Ks[-2]))
    assert ("By the argument attaining t(v), away from the cusp (|ln(v/v_hs)| ≥ 0.05): τ²/v (0 rules): the limit over 1/t(v), minus 1, —; "
            "the ratio's excess at K = 160 —, its exponent —; with a K⁻¹ term —, the two fits bracketing 1/t(v) at 0 of the 0; vσ*²/τ² (1 "
            f"rules): the limit over 1/t(v), minus 1, {WM.out_e(float(cells[-2]), False)} … {WM.out_e(float(cells[-2]), True)}; the ratio's "
            f"excess at K = 160 {WM.out_e(y1, False)} … {WM.out_e(y1, True)}, its "
            f"exponent {np.floor(ex_ * 100) / 100:.2f} … {np.ceil(ex_ * 100) / 100:.2f}; with a K⁻¹ term {WM.out_e(L2, False)} … {WM.out_e(L2, True)}, the two fits bracketing 1/t(v) at {int(L1 * L2 <= 0)} of the 1. "
            "Near it: none.") in s
    assert (f"With a K⁻¹ term instead (L + b/K + c/K² through the same K) off the cusp within {abs(L2):.1e} (the farthest: τ/σ = 1, "
            "Mercer frame): the second fit's largest deviation;") in s
    res = [(1 - 2 / (5.8 * K)) / (1 - 1 / (4 * K)) - 1 for K in Ks]
    assert ("At v_hs the sketch's integral spreads over every angle, so e_K ∝ K^{1/4}t^K and the ratio per kernel is (1 − 1/(4K))/t + "
            "O(K⁻²): the ratio over (1 − 1/(4K))/t, minus 1, at K = 20, 40, 60, 80, 120, 160 (its successive ratios; K⁻²'s 0.250, 0.444, "
            "0.562, 0.444, 0.562): τ/σ = 1: " + ', '.join(f"{v:+.1e}" for v in res) + ' ('
            + ', '.join(f"{res[i + 1] / res[i]:.3f}" for i in range(len(res) - 1)) + ').') in s
    cr = 3.0 / 140 ** 2                                                  # an excess changing sign between K = 120 and 160
    A6b = [A6[0], dict(r=1.0, name='0.5·v_hs', v=1.0, v_hs=0.5, rate=0.25, arg=1, t=0.17, rows=rows(4.0 - cr, 3.0, -2.0)),
           dict(r=1.0, name='1.01·v_hs', v=1.0, v_hs=1.01, rate=0.25, arg=2, t=0.17, rows=rows(4.0, 3.0, -2.0))]
    s2 = '\n'.join(WM.a6_lines(A6b))
    assert "τ²/v (1 rules): the limit over 1/t(v), minus 1, " in s2 and ", its exponent — (1 changing sign); with a K⁻¹ term " in s2
    assert f"Near it: τ/σ = 1, 1.01·v_hs (v/v_hs 1.010000): {y1:+.1e}." in s2 and "at all 3 rules" in s2 and 'At v_hs' not in s2
    A5 = dict(series=[dict(r=1.0, frames=dict(Mercer=[dict(K=K, S=2.0 * 10.0 ** -K) for K in range(19, 41)]))])
    A6[0]['rows'] = [dict(K=K, prec=100, e0=f"{np.sqrt(10.0 ** -(K - 1)):.17e}", e1=f"{np.sqrt(10.0 ** -K):.17e}", ratio=np.sqrt(10.0) + 1 / K ** 2)
                     for K in WM.A6_K]
    s = '\n'.join(WM.a6_lines(A6[:1], A5))
    assert "The Mercer frame's rule's exact squared error over its frame series (A5) is one constant per τ/σ: within " in s
    assert "at τ/σ = 1 (K = 19, 20, 39, 40)" in s and WM.a6_lines(None) == [] and WM.a6_lines([]) == []


def test_main_only_b5_and_a6(monkeypatch, tmp_path):
    """Referee round 84: `main --only B5,A6` passes the raw's B and B3 to `part_B5`, runs `part_A6`, and stores both; the report takes
    them in the order of its keys."""
    import exp_cohesion_leak as E
    import inspect
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), B=[dict(r=1.0, var='base')], B3=[dict(r=1.0, var='base', k=None)])
    got, saved, rep = [], [], []
    params = list(inspect.signature(WM.report).parameters)
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'part_B5', lambda B, B3: (got.append((list(B), list(B3))), ['b5'])[1])
    monkeypatch.setattr(WM, 'part_A6', lambda: ['a6'])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: rep.append(a) or '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'B5,A6'])
    WM.main()
    assert got == [([dict(r=1.0, var='base')], [dict(r=1.0, var='base', k=None)])] and saved[-1]['B5'] == ['b5'] and saved[-1]['A6'] == ['a6']
    assert rep[-1][8:10] == (['b5'], ['a6']) and params == ['A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X', 'A6C']


def test_a6e_the_generating_function(monkeypatch):
    """Referee round 84 (D3ce-3): `dpi` is π (50 digits), `dcos` the cosine; `gen_ratio` at K = 1 — Q_1[e^{ξz}] = 1 — is
    (1 − e^{−ξ²/2})/((ξ²/2)e^{−ξ²/4}) on the real axis and (e^{y²/2} − 1)/((y²/2)e^{y²/4}) on the imaginary; the real `gen_ratio` at
    K = 20 reproduces the raw; `part_A6E` covers A6E_K × both axes × A6E_S in 3K + 60 digits; `a6e_lines` prints the ranges."""
    from decimal import Decimal, localcontext
    import math
    import exp_cohesion_leak as E
    with localcontext() as ctx:
        ctx.prec = 60
        pi = WM.dpi()
        assert abs(pi - Decimal('3.14159265358979323846264338327950288419716939937510582097')) < Decimal(10) ** -55
        for t in (0.3, 2.0, 7.5, 40.0):
            assert abs(float(WM.dcos(Decimal(t), pi)) - math.cos(t)) < 1e-13
    for s in (0.5, 1.0, 2.0):
        x2 = s * s
        assert abs(WM.gen_ratio(1, s, 'real', 60) - (1 - math.exp(-x2 / 2)) / (x2 / 2 * math.exp(-x2 / 4))) < 1e-13
        assert abs(WM.gen_ratio(1, s, 'imag', 60) - (math.exp(x2 / 2) - 1) / (x2 / 2 * math.exp(x2 / 4))) < 1e-13
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    x = next(x for x in raw['A6E'] if (x['K'], x['axis']) == (20, 'imag'))
    assert WM.gen_ratio(20, 1.5, 'imag', 3 * 20 + 60) == x['ratio'][list(x['s']).index(1.5)]
    got = []
    monkeypatch.setattr(WM, 'gen_ratio', lambda K, s, ax, prec: got.append((K, s, ax, prec)) or K + s)
    A6E = WM.part_A6E()
    assert got == [(K, s, ax, 3 * K + 60) for K in WM.A6E_K for ax in ('real', 'imag') for s in WM.A6E_S]
    assert [(y['K'], y['axis']) for y in A6E] == [(K, ax) for K in WM.A6E_K for ax in ('real', 'imag')] and WM.A6E_S[-1] == 2.0
    s = '\n'.join(WM.a6e_lines(A6E))
    assert "K = 20, the real axis: 20.100–22.000 (at s = 0.5, 1, 1.5, 2: 20.500, 21.000, 21.500, 22.000)" in s
    assert "K = 80, the imaginary axis: 80.100–82.000" in s and "Laplace's method" in s and WM.a6e_lines(None) == []
    assert s.startswith("\nThe step of the rate's Laplace sketch (A6E)") and "and D‴.62 (b) for every rule without one, by hand)." in s   # round 132


def test_round88_precision_and_the_md(monkeypatch):
    """Referee round 88 (D3ci-4, P01, P02; D3ci-5g): `a6_error` at τ/σ = 1/4 (the Fock rule), K = 59 and 60, in the job's own precision
    (2K·log₁₀(1/t) + 60 digits) reproduces the raw's e₀ and e₁ to their 26 digits — where `gh_nodes`' default 110 digits do not suffice —;
    `gen_ratio` at K = 80, s = 0.1 on the real axis in 3K + 60 digits reproduces the raw; and the md is `report` of the raw, byte for byte."""
    from decimal import Decimal, localcontext
    import math
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    x = next(x for x in raw['A6'] if (x['r'], x['name']) == (0.25, 'Fock'))
    y = next(y for y in x['rows'] if y['K'] == 60)
    t = x['t']
    assert y['prec'] == int(math.ceil(2 * 60 * math.log10(1 / t))) + 60 and y['prec'] > 200
    with localcontext() as ctx:
        ctx.prec = y['prec']
        v = next(v for nm, v in WM.a6_scalings(0.25) if nm == 'Fock')
        e0, e1 = WM.a6_error(0.25, v, 59, y['prec']), WM.a6_error(0.25, v, 60, y['prec'])
        assert format(e0.sqrt(), '.25e') == y['e0'] and format(e1.sqrt(), '.25e') == y['e1']
    z = next(z for z in raw['A6E'] if (z['K'], z['axis']) == (80, 'real'))
    assert WM.gen_ratio(80, 0.1, 'real', 3 * 80 + 60) == z['ratio'][list(z['s']).index(0.1)]
    md = open(os.path.join(WM.OUT, 'width_mismatch.md'), encoding='utf-8').read()
    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X', 'A6C')
    assert WM.report(*[raw.get(k) for k in keys]) == md
    # D3ci-3: the rules with a Mehler frame, |ρ| < 1, ρ = 1 + 2/τ² − 2/v at σ = 1 (v < τ² and v > τ²/(1 + τ²))
    rho = lambda x: 1 + 2 / x['r'] ** 2 - 2 / float(x['v'])
    fin = [x for x in raw['A6'] if float(x['v']) < x['r'] ** 2 and float(x['v']) > x['r'] ** 2 / (1 + x['r'] ** 2)]
    assert len(fin) == 22 and all(abs(rho(x)) < 1 for x in fin) and [sum(1 for x in fin if x['r'] == r_) for r_ in (0.25, 0.5, 1, 2)] == [4, 4, 6, 8]
    assert sorted({x['name'] for x in raw['A6'] if abs(rho(x) - 1) < 1e-9}) == ['moment matching']
    # D3db-4: Fock and Mercer-zero coincide at τ/σ = 2 (one distinct rule, counted once in the distinct count)
    dst = {(x['r'], round(float(x['v']), 12)) for x in raw['A6']}
    assert len(dst) == 35 and len({(x['r'], round(float(x['v']), 12)) for x in fin}) == 21
    assert ("A Mehler frame (|ρ| < 1, ρ = 1 + 2σ²/τ² − 2σ²/v: v ∈ (τ²σ²/σ*², τ²); D‴.61) exists for 22 of the 36 rules (21 of 35 distinct; 4 at τ/σ = 0.25, 4 "
            "at τ/σ = 0.5, 6 at τ/σ = 1, 8 at τ/σ = 2); ρ = 1 at 4 (moment matching), |ρ| > 1 at 10.") in md
    three = [next(x for x in raw['A6'] if (x['r'], x['name']) == k) for k in ((0.25, 'Fock'), (0.25, 'moment matching'), (0.25, '0.5·v_hs'))]
    assert ("exists for 1 of the 3 rules (1 of 3 distinct; 1 at τ/σ = 0.25); ρ = 1 at 1 (moment matching), |ρ| > 1 at 1.") in '\n'.join(WM.a6_lines(three))


def test_round88_the_second_fit_and_the_censored_ends():
    """Referee round 88, this round's surviving mutants (F24, B02): the K⁻¹ fit's farthest rule is its own, not the first fit's (a rule
    exact in L + c/K² + d/K³ and one exact in L + b/K + c/K²); a direct B5 end is starred at two accepted steps and not at three."""
    from decimal import Decimal
    Ks = np.array(WM.A6_K, float)
    rows = lambda c, p: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0 + c * K ** p) for K in Ks]
    A6 = [dict(r=1.0, name='cubic', v=1.0, v_hs=2.0, rate=0.25, arg=2, t=0.17, rows=rows(3.0, -3.0)),
          dict(r=1.0, name='linear', v=1.5, v_hs=3.0, rate=0.25, arg=2, t=0.17, rows=rows(3.0, -1.0))]   # distinct (r, v): round 107
    fit = lambda rr, p, q: float(np.linalg.solve(np.array([[1.0, K ** p, K ** q] for K in Ks[-3:]]), rr[-3:])[0]) * 0.25 - 1
    rr = {x['name']: [y['ratio'] for y in x['rows']] for x in A6}
    d1, d2 = {k_: fit(v, -2.0, -3.0) for k_, v in rr.items()}, {k_: fit(v, -1.0, -2.0) for k_, v in rr.items()}
    assert abs(d1['cubic']) < 1e-13 < abs(d1['linear']) and abs(d2['linear']) < 1e-13 < abs(d2['cubic'])
    s = '\n'.join(WM.a6_lines(A6))
    assert f"within {WM.up_e(d1['linear'])} of the formula's 1/t(v) at all 2 rules (the farthest: τ/σ = 1, linear)" in s
    assert f"off the cusp within {WM.up_e(d2['cubic'])} (the farthest: τ/σ = 1, cubic): the second fit's largest deviation;" in s
    bst = Decimal('2.0e-10')
    mk = lambda x: str(bst * (1 + Decimal(repr(x))))
    B4 = [dict(r=0.5, k=2, start=s_, e_opt=3e-10, **{k_: dict(rms=str(bst * 2), ratio=2.0, it=1, acc=1, stop='stall', p=[0.0]) for k_, _ in WM.B4_ENDS})
          for s_ in ('block', 'optimum')]
    B4[1]['anchored ld'] = dict(rms=str(bst), ratio=1.0, rms_d=mk(1e-12), it=1, acc=1, stop='stall', p=[0.0])
    acc = dict(zip(WM.B5_KAPPA, (5, 5, 5, 5, 5, 4, 3, 2)))
    B5 = [dict(r=0.5, k=2, theta=0.0, mode='direct', e_start=1.5, rows=[dict(kappa=k_, rms=mk(3e-6 * k_ ** 2 if k_ else 3e-9), ratio=1.0, it=9,
                                                                             acc=acc[k_], stop='stall', p=[0.0]) for k_ in WM.B5_KAPPA])]
    ex = lambda x: WM.b4_excess(mk(x), bst)
    s = '\n'.join(WM.b5_lines(B5, B4))
    assert f"| {ex(3e-6):.1e} | {ex(3e-6 * 4):.1e} | {ex(3e-6 * 16):.1e}* |" in s
    assert "* at most two accepted steps: the end is the start's, censored (the slopes include it)." in s
    B5[0]['rows'][-1]['acc'] = 3
    s = '\n'.join(WM.b5_lines(B5, B4))
    assert f"| {ex(3e-6 * 16):.1e} |" in s and '*' not in s.replace('σ*', '')



def test_round92_frames_and_d61_rstar():
    """Referee round 92 (D3cm-1, D3cm-2): A6's deviations split by frame; the Fock and Mercer-zero thresholds below D‴.61's r*;
    D61_RSTAR is D‴.61's r* from its verified radius, rounded down."""
    import math
    import exp_cohesion_leak as E
    import exp_gh_peano as GP
    assert WM.D61_RSTAR == math.floor(GP.r_star(1.8927) * 1e4) / 1e4 == 0.8052
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    s = '\n'.join(WM.a6_lines(raw['A6'], raw.get('A5')))
    assert "Off the cusp, the 18 rules (17 distinct) with a frame, whose limits D‴.61 (iv) proves" in s and "at all 32 rules (31 distinct)" in s and "the 14 rules without a frame" in s
    rf = next(r_ for r_ in np.linspace(0.01, 10, 100000) if WM.rate_conj(r_, 2 * r_ * r_ / (2 + r_ * r_))[0] >= WM.D61_RSTAR)
    assert f"σ*/σ < {math.sqrt(1 + rf * rf):.3f}" in s or f"σ*/σ < {math.sqrt(1 + (rf - 1e-4) ** 2):.3f}" in s
    assert f"the frame ends above it once σ*/σ > {math.ceil(math.sqrt((1 + 0.8052) / (1 - 0.8052)) * 1000) / 1000:.3f}" in s


def test_round96_outward_rounding_the_split_cells_and_the_thresholds():
    """Referee round 96 (D3cq-2, D3cq-4, D3cq-7): `up_e` rounds an upper bound up at its last digit (7.547e-6 → 7.6e-06, never below),
    `dn_f`/`up_f` round thresholds inward; the split by frame, branch and nearness counts and bounds planted rows per cell; the Fock,
    Mercer-zero and frame-end thresholds against their closed forms (4.18735, 6.51445, 3.04417 at r* = 0.8052)."""
    import math
    assert WM.up_e(7.547e-6) == '7.6e-06' and WM.up_e(4.449e-5) == '4.5e-05' and WM.up_e(3.134e-7) == '3.2e-07'
    assert WM.up_e(-6.0734e-5) == '6.1e-05' and WM.up_e(9.96e-5) == '1.0e-04' and WM.up_e(2.0e-4) == '2.0e-04' and WM.up_e(0) == '0.0e+00'
    for v in np.random.default_rng(96).uniform(1e-9, 1e-3, 200):
        assert float(WM.up_e(v)) >= v and float(WM.up_e(v)) - v <= 1.0001 * 10 ** (math.floor(math.log10(v)) - 1)
    assert WM.dn_f(4.18735) == '4.187' and WM.dn_f(6.5149) == '6.514' and WM.up_f(3.04417) == '3.045' and WM.up_f(3.0440) == '3.044'
    r = WM.D61_RSTAR
    tf = 4.066188                                                       # the Fock rule at t(v) = r*: τ/σ, σ*/σ = (1 + τ²)^{1/2} = 4.18735
    assert abs(WM.rate_conj(tf, 2 * tf * tf / (2 + tf * tf))[0] - r) < 1e-6 and abs(math.sqrt(1 + tf * tf) - 4.18735) < 1e-5
    mz_tau = math.sqrt((((1 + r) / (1 - r)) ** 2 - 1) / 2)              # the Mercer-zero rule's τ/σ at t(v) = r*, closed form
    assert abs(WM.rate_conj(mz_tau, mz_tau ** 2 / math.sqrt(1 + 2 * mz_tau ** 2))[0] - r) < 1e-9
    end = math.sqrt((1 + r) / (1 - r))
    assert abs(mz_tau - 6.51445) < 1e-5 and abs(end - 3.04417) < 1e-5
    A6 = [dict(r=x['r'], name=x['name'], v=x['v'], v_hs=x['v_hs'], rate=x['rate'], arg=x['arg'], t=x['t'], rows=x['rows'])
          for x in WM.E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))['A6']]
    s = '\n'.join(WM.a6_lines(A6))
    assert ("By frame, the branch attaining t(v) and nearness to the cusp (far: |ln(v/v_hs)| ≥ 0.05; referee round 96): with a frame, τ²/v, "
            "far (8, 7 distinct): within 6.1e-05 and 1.4e-04, bracketing at 8;") in s
    assert "without, τ²/v, far (5): within 7.6e-06 and 4.5e-05, bracketing at 5;" in s and "at v_hs within 3.2e-07 of 1/t" in s
    assert "σ*/σ < 4.187 (τ/σ < 4.066), τ/σ < 6.514, the frame ends above it once σ*/σ > 3.045" in s
    rows = lambda c, p: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0 + c * K ** p) for K in WM.A6_K]
    pl = [dict(r=1.0, name='a', v=1.2, v_hs=2.0, rate=0.25, arg=1, t=0.17, rows=rows(3.0, -1.5)),     # ρ = 1 + 2 − 2/1.2: |ρ| > 1
          dict(r=1.0, name='b', v=0.9, v_hs=1.0, rate=0.25, arg=2, t=0.17, rows=rows(3.0, -1.5)),     # near the cusp, ρ = 0.78: a frame
          dict(r=1.0, name='c', v=0.8, v_hs=2.0, rate=0.25, arg=2, t=0.17, rows=rows(-3.0, -1.5))]
    s = '\n'.join(WM.a6_lines(pl))
    assert "without, τ²/v, far (1)" in s and "with a frame, vσ*²/τ², near (1)" in s and "with a frame, vσ*²/τ², far (1)" in s


def test_round96_outward_rounding_planted_edges(monkeypatch):
    """Round 96's own mutants L03, L11: `up_e` stays an upper bound a hair above a printable value (7.5000001e-6 → 7.6e-06, 7.5e-6 itself
    → 7.5e-06); the Fock threshold is printed rounded down, at an r* where rounding to nearest would round it up."""
    import math
    assert WM.up_e(7.5000001e-6) == '7.6e-06' and WM.up_e(7.5e-6) == '7.5e-06' and WM.up_e(1.0000001e-3) == '1.1e-03'
    A6 = [dict(r=x['r'], name=x['name'], v=x['v'], v_hs=x['v_hs'], rate=x['rate'], arg=x['arg'], t=x['t'], rows=x['rows'])
          for x in WM.E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))['A6']]

    def fock_tau(r):
        lo, hi = 1e-3, 50.0
        for _ in range(100):
            m = (lo + hi) / 2
            lo, hi = (m, hi) if WM.rate_conj(m, 2 * m * m / (2 + m * m))[0] < r else (lo, m)
        return lo
    seen = 0
    for r in (0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86):
        monkeypatch.setattr(WM, 'D61_RSTAR', r)
        x = math.sqrt(1 + fock_tau(r) ** 2)
        s = '\n'.join(WM.a6_lines(A6))
        assert f"σ*/σ < {math.floor(x * 1000) / 1000:.3f} (τ/σ < " in s
        seen += f"{x:.3f}" != f"{math.floor(x * 1000) / 1000:.3f}"
    assert seen > 0


def test_round101_outward_ranges_the_cells_the_outlier_and_part_v_s_tier():
    """Referee round 101 (D3cv-1 … D3cv-5): `out_e` rounds a range's lower end down and its upper end up (−3.591e-7 → −3.6e-07 below,
    −3.5e-07 above); the cells count distinct rules (Fock ≡ Mercer-zero at τ/σ = 2); the fast branch's outlier is printed with both fits'
    ratios to the proved rules' largest; the slow branch's rules by t(v), frame or not; a5's t(v) maximum compared with the tiers, not
    fixed text; with D‴.61's part V (r* = 0.9900, an earlier tier since VA: round 121) the Fock, Mercer-zero and frame-end thresholds
    against their closed forms (19.924, 140.712, 14.107), and with VA no threshold; planted cells' maxima and brackets."""
    import math
    import exp_cohesion_leak as E
    assert WM.out_e(-3.591e-7, False) == '-3.6e-07' and WM.out_e(-3.591e-7, True) == '-3.5e-07'
    assert WM.out_e(5.652e-9, False) == '+5.6e-09' and WM.out_e(5.652e-9, True) == '+5.7e-09' and WM.out_e(9.96e-5, True) == '+1.0e-04'
    assert WM.out_e(-9.96e-5, False) == '-1.0e-04' and WM.out_e(0.0, True) == '+0.0e+00' and WM.out_e(2.0e-4, True) == '+2.0e-04'
    for v in np.random.default_rng(101).uniform(-1e-3, 1e-3, 400):
        lo, hi = float(WM.out_e(v, False)), float(WM.out_e(v, True))
        u = 1.0001 * 10 ** (math.floor(math.log10(abs(v))) - 1)
        assert lo <= v <= hi and v - lo <= u and hi - v <= u
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    A6 = [dict(r=x['r'], name=x['name'], v=x['v'], v_hs=x['v_hs'], rate=x['rate'], arg=x['arg'], t=x['t'], rows=x['rows']) for x in raw['A6']]
    s = '\n'.join(WM.a6_lines(A6, raw['A5']))
    assert "with a frame, τ²/v, far (8, 7 distinct): within 6.1e-05 and 1.4e-04, bracketing at 8;" in s
    assert "without, vσ*²/τ², far (8): within 8.8e-08 and 3.6e-07, bracketing at 4;" in s and "without, vσ*²/τ², near (1)" in s
    assert ("the largest deviations without a frame, both at τ/σ = 0.25, 2·v_hs, are 41.4 and 2.0 times the proved rules' largest") in s
    sl = s.split("The slow branch far from the cusp by t(v) (both fits' deviations): ")[1].split('.\n')[0].split('; ')
    ts = [float(x.split(' ')[1]) for x in sl]
    assert len(sl) == 13 and ts == sorted(ts, reverse=True) and sum('(no frame)' in x for x in sl) == 5
    assert "t(v) 0.500 at τ/σ = 2, Fock (a frame): +1.42e-05, -6.15e-05; t(v) 0.500 at τ/σ = 2, Mercer-zero (a frame): +1.42e-05, -6.15e-05" in s
    assert (f"proved by D‴.61 (iv) where a rule has a Mehler frame ({WM.D61_STATUS}); by D‴.62 (b) where it has none, without computer "
            "assistance (referee round 132))") in s and "t(v) < 0.9900" not in s and "conjectural otherwise" not in s
    r = WM.D61_RSTAR_V
    mz_tau, end = math.sqrt((((1 + r) / (1 - r)) ** 2 - 1) / 2), math.sqrt((1 + r) / (1 - r))
    assert abs(mz_tau - 140.71247) < 1e-5 and abs(end - 14.106736) < 1e-6
    lo, hi = 1.0, 100.0
    for _ in range(100):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if WM.rate_conj(m, 2 * m * m / (2 + m * m))[0] < r else (lo, m)
    assert WM.dn_f(math.sqrt(1 + lo * lo)) == '19.924' and WM.dn_f(lo) == '19.899'
    assert ("The Fock rule (ρ = 0) and the Mercer-zero rule have a frame at every τ/σ, so with D‴.61's part VA (its (iv) at every frame, "
            f"r* = 1; {WM.D61_STATUS}; referee rounds 101, 107) their rates — λ/(2 − λ) and q — are proved at every τ/σ, no threshold left. "
            "The earlier tiers (thresholds rounded inward: referee round 96): with part V (r* = 0.9900) while below that r*, σ*/σ < 19.924 "
            "(τ/σ < 19.899) and τ/σ < 140.712 — near the frame interval's ends t(v) → λ, above that r* once σ*/σ > 14.107 —; with the "
            "branch and bound alone (r* = 0.8052): σ*/σ < 4.187") in s and "conjectural beyond" not in s
    a5 = '\n'.join(WM.a5_series_lines(raw['A5']))
    tv = lambda b: b * (b - 1) / (b * b + b + 2)                         # the Mercer frame rule's t(v), b = (1 + 2τ²/σ²)^{1/2}
    tmax = max(tv(math.sqrt(1 + 2 * y['r'] ** 2)) for y in raw['A5']['far'])
    assert (f"the frame's t(v) at most {WM.up_f(tmax, 4)}, below its branch and bound's r* = 0.8052 already, and with VA (iv) holds at "
            f"every frame; {WM.D61_STATUS}; referee rounds 96, 101, 107)") in a5 and WM.up_f(tmax, 4) == '0.4286'
    Ks = np.array(WM.A6_K, float)
    rows = lambda c, p: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0 + c * K ** p) for K in Ks]
    pl = [dict(r=1.0, name='a', v=1.2, v_hs=2.0, rate=0.25, arg=1, t=0.17, rows=rows(3.0, -1.5)),
          dict(r=1.0, name='b', v=0.9, v_hs=1.0, rate=0.25, arg=2, t=0.17, rows=rows(3.0, -1.5)),
          dict(r=1.0, name='c', v=0.8, v_hs=2.0, rate=0.25, arg=2, t=0.17, rows=rows(-3.0, -1.5))]
    rr = lambda c: [4.0 + c * K ** -1.5 for K in Ks]
    f1 = lambda c: float(np.linalg.solve(np.array([[1.0, K ** -2.0, K ** -3.0] for K in Ks[-3:]]), rr(c)[-3:])[0]) * 0.25 - 1
    f2 = lambda c: float(np.linalg.solve(np.array([[1.0, K ** -1.0, K ** -2.0] for K in Ks[-3:]]), rr(c)[-3:])[0]) * 0.25 - 1
    s = '\n'.join(WM.a6_lines(pl))
    assert (f"without, τ²/v, far (1): within {WM.up_e(abs(f1(3.0)))} and {WM.up_e(abs(f2(3.0)))}, bracketing at "
            f"{int(f1(3.0) * f2(3.0) <= 0)}") in s
    assert (f"with a frame, vσ*²/τ², far (1): within {WM.up_e(abs(f1(-3.0)))} and {WM.up_e(abs(f2(-3.0)))}, bracketing at "
            f"{int(f1(-3.0) * f2(-3.0) <= 0)}") in s


def test_round101_a5_label_compared(monkeypatch):
    """Rounds 101 (D3cv-5), 121: a5's tier is computed — the first of D‴.61's steps whose r* exceeds the frame's largest t(v) (its branch
    and bound, V, VA); with the earlier tiers lowered below it the sentence names the later step. Whatever the tiers, a frame rule is
    proved (even at t(v) = 0.995, beyond V's r*) and a rule without a frame is D‴.62 (b)'s, by hand (A6; round 132), A6X's frame/no-frame
    labels unchanged."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    a5 = lambda A5: '\n'.join(WM.a5_series_lines(A5))
    assert "the frame's t(v) at most 0.4286, below its branch and bound's r* = 0.8052 already, and with VA (iv) holds at every frame;" in a5(raw['A5'])
    assert WM.d61_step(0.5) == ('its branch and bound', 0.8052) and WM.d61_step(0.9) == ('V', 0.99) and WM.d61_step(0.995) == ('VA', 1)
    monkeypatch.setattr(WM, 'D61_RSTAR', 0.4)
    assert "the frame's t(v) at most 0.4286, below V's r* = 0.9900 already, and with VA (iv) holds at every frame;" in a5(raw['A5'])
    monkeypatch.setattr(WM, 'D61_RSTAR_V', 0.4)
    assert f"the frame's t(v) at most 0.4286, below VA's r* = 1; {WM.D61_STATUS}" in a5(raw['A5'])
    A5 = dict(raw['A5'])                                                # N13: at τ/σ = 1 t(v) = (3 − √3)/(5 + √3) = 0.188346…, up 0.1884
    lim = lambda r_: (lambda b: (b * b + b + 2) / (b * (b - 1)))(np.sqrt(1 + 2 * r_ ** 2))
    A5['far'] = [dict(r=1.0, K=K, prec=140, ratio=lim(1.0) - 0.5 / K ** 2, low=1e-50) for K in (60, 80, 100, 120)]
    assert "the frame's t(v) at most 0.1884, below its branch and bound's r* = 0.4000 already" in a5(A5)
    Ks = np.array(WM.A6_K, float)                                       # the tiers lowered: a frame at t(v) = 0.995 proved, ρ = 4/3 not
    rows = lambda c: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0 + c * K ** -1.5) for K in Ks]
    pl = [dict(r=1.0, name='a', v=1.2, v_hs=2.0, rate=0.25, arg=1, t=0.17, rows=rows(3.0)),       # ρ = 3 − 2/1.2: no frame
          dict(r=1.0, name='c', v=0.8, v_hs=2.0, rate=0.995, arg=2, t=0.17, rows=rows(-3.0))]     # ρ = 1/2: a frame
    s = '\n'.join(WM.a6_lines(pl))
    assert "Off the cusp, the 1 rules with a frame, whose limits D‴.61 (iv) proves" in s
    assert ("Off the cusp, the 1 rules without a frame, whose limits D‴.62 (b) proves by hand (the fits' deviations their bias too; 'the "
            "proved rules' below are the frame rules, the calibration predating D‴.62)") in s and "t(v) ≥ r*" not in s
    assert "conjectur" not in s                                         # round 132: no rule's rate is conjectural any more
    rx = [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0) for K in WM.A6X_K]
    pl = [dict(x, arg=2) for x in pl]                                   # A6X: the fast branch
    s = '\n'.join(WM.a6x_lines([dict(x, rows=rx) for x in pl], pl))
    assert "(no frame; 100 digits)" in s and "(a frame; 100 digits)" in s and "at the 1 rule without a frame" in s and "at the proved rule" in s
    assert "(A6X extends 1 of the 1 fast-branch rules without a frame)" in s


def test_round101_a6x_the_job_s_k_and_the_lines(monkeypatch):
    """Round 101 (D3cv-1): `a6_job` with a third argument runs at those K (the same record otherwise: A6's first two K reproduce the raw's
    rows); `part_A6X` maps A6X_RULES at A6X_K; `a6x_lines` fits both windows and says whether every deviation shrinks."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    x = next(z for z in raw['A6'] if (z['r'], z['name']) == (2.0, '2·v_hs'))
    out = WM.a6_job((2.0, 8, (20, 40)))
    assert [dict(y) for y in out['rows']] == [dict(y) for y in x['rows'][:2]] and out['name'] == '2·v_hs'
    assert [WM.a6_scalings(r_)[i][0] for r_, i in WM.A6X_RULES] == ['2·v_hs', 'moment matching', 'moment matching', 'moment matching', '2·v_hs',
                                                                     'Mercer frame', 'Mercer frame']     # round 123: the known cases
    seen = []

    class P:
        def __init__(self, n): pass
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def map(self, f, args, chunksize=1): seen.extend(args); return ['ok']
    import multiprocessing
    monkeypatch.setattr(multiprocessing, 'Pool', P)
    assert WM.part_A6X() == ['ok'] and seen == [(r_, i, WM.A6X_K) for r_, i in WM.A6X_RULES]
    Ks = np.array(list(WM.A6_K) + list(WM.A6X_K), float)
    mk = lambda c: [dict(K=int(K), prec=100, e0='1e-10', e1='1e-11', ratio=4.0 * (1 + c / K ** 2 + 5.0 / K ** 4)) for K in Ks]
    a6 = [dict(r=1.0, name='q', v=0.5, v_hs=2.0, rate=0.25, arg=2, t=0.17, rows=mk(0.3)[:len(WM.A6_K)])]
    a6x = [dict(a6[0], rows=mk(0.3)[len(WM.A6_K):])]
    s_ = '\n'.join(WM.a6x_lines(a6x, a6))
    assert "every deviation shrinks as the window moves out" in s_ and "at the 1 rule without a frame" in s_ and "(no frame; 100 digits)" in s_
    ex = [(4.0 * (1 + 0.3 / K ** 2 + 5.0 / K ** 4) * 0.25 - 1) * K ** 2 for K in (160, 200, 240, 320)]
    assert "the excess times K² at K = 160, 200, 240, 320 " + ', '.join(f"{v:+.4f}" for v in ex) in s_
    a6x = [dict(a6[0], rows=[dict(y, ratio=y['ratio'] + 1e-3) for y in mk(0.3)[len(WM.A6_K):]])]
    assert "NOT every deviation shrinks" in '\n'.join(WM.a6x_lines(a6x, a6))
    md = open(os.path.join(WM.OUT, 'width_mismatch.md'), encoding='utf-8').read()          # the run: the referee's extension reproduced
    assert "τ/σ = 0.25, 2·v_hs (no frame; 1225 digits): +8.71e-08 → +1.55e-12 and +3.60e-07 → -3.34e-09" in md
    assert "-0.0928, -0.0927, -0.0927, -0.0926" in md and "every deviation shrinks as the window moves out" in md
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    assert [(x['r'], x['name']) for x in raw['A6X']] == [(r_, WM.a6_scalings(r_)[i][0]) for r_, i in WM.A6X_RULES]


def test_round107_d61_status_pinned_to_its_raw_and_the_a6x_verdict():
    """Round 107 (D3db-1, -2, -5), round 121: `D61_RSTAR_V` is D‴.61's r*(1.995) rounded down (V's tier), and `D61_STATUS` ("three
    steps rigorous in ball arithmetic") holds on D‴.61's raw — every band of V and every radius of D's branch and bound rechecked in ball
    arithmetic, the rigorous bounds negative, no open box; VA's every piece without a failed box, its outer polyline's 76 boxes from where
    V ends, its end piece and hand-off passing and tiled, every soundness check passing; with VA (iv) at every frame, the tails' and the
    mass's c² = 4r ln(1/r)/(1 − r) and 4r below 4 for every r < 1 —; A6X's verdict is `report()`'s on the raw, rule by rule (round 123),
    and a planted offset of the outlier's limit (the referee's 5·10⁻⁸), which still shrinks as the window moves out, fails the transient's
    check; the report's D section names each earlier tier with its own r* (round 123's X12)."""
    import math
    import exp_cohesion_leak as E
    import exp_gh_peano as GP
    assert WM.D61_RSTAR_V == math.floor(GP.r_star(1.995) * 1e4) / 1e4 == 0.99
    gp = E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    bs, db = gp['V']['bands'], list(gp['DB'])
    assert gp['V'].get('ball') and bs and all(b['ball_worst'] < 0 and b['ball_end'] < 0 and b['fails'] == 0 for b in bs)
    assert db and all(x['ball_worst'] < 0 and x['open'] == 0 for x in db) and [x['R'] for x in db] == list(GP.DOM_BB)
    assert max(b['R2'] for b in bs) == 1.995 and "ball arithmetic" in WM.D61_STATUS and "two certificates" in WM.D61_STATUS
    assert WM.D61_STATUS.startswith("by hand since referee round 148") and "now independent cross-checks" in WM.D61_STATUS   # round 148
    by = {r['piece']: r for r in gp['VA']}                              # round 121: part VA on 1.995 ≤ |η| < 2
    assert all(by[p]['failed'] == 0 for p in ('A-J', 'A-FC', 'B-J', 'B-F', 'O')) and by['O']['boxes'] == 76
    assert by['end']['ok'] and by['handoff']['failed'] == 0 and by['handoff']['tiled'] and by['O']['R1'] == max(b['R2'] for b in bs)
    snd = [r for r in gp['VA'] if str(r['piece']).startswith('sound ')]
    assert len(snd) == 4 and all(r['failed'] == 0 for r in snd)
    assert WM.D61_RSTAR_VA == 1 and "VA" in WM.D61_STATUS and "referee rounds 120, 121" in WM.D61_STATUS   # round 129 (D3dx-5): as D58
    assert "Appendix A.0 assumed correct" in WM.D61_STATUS and "V an independent cross-check" in WM.D61_STATUS and "adoption reviewed in round 122" in WM.D61_STATUS and "round 121)" not in WM.D61_STATUS
    rr = np.linspace(1e-6, 1 - 1e-6, 20001)                             # c_need(r)² < 4 on r < 1: the open disc |η| < 2 suffices
    assert np.all(4 * rr * np.log(1 / rr) / (1 - rr) < 4) and np.all(4 * rr < 4)
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X', 'A6C')
    rep = WM.report(*[raw.get(k) for k in keys])
    ok4 = ("(τ/σ, rule) = (0.5, moment matching), (1, moment matching), (2, moment matching), (2, 2·v_hs): decreasing in magnitude through eight "
           "terms, the last within ten times its last change (that, the resolution: 4.1e-18, 9.0e-18, 1.9e-19, 2.3e-18) — 1/t(v) to that resolution; ")
    assert (ok4 + "(1, Mercer frame), (2, Mercer frame): NOT decreasing in magnitude through eight terms — inconclusive: the limit, the premise, "
            "or the window (the expansion not yet asymptotic at the order's smallest K); at (1, Mercer frame), (2, Mercer frame) the window: on "
            "K ≥ 160 the orders decrease (A6C), to -1.1e-21, -7.3e-24 at sixteen terms, the last within ten times its last change, the excess's "
            "coefficients growing 12.4–17.5, 22.7–24.2 times per order (the last four ratios c_k/c_{k−1}; at most 6.3 at A6C's rules where the "
            "check passes), so that at the eight-term order's smallest K, 40, the terms c_k/K^k shrink by only 0.31–0.44, 0.56–0.61 per order: the "
            "expansion is not yet asymptotic there; the known cases, the 3 rules proved: the check passes at 1 of 3 — a failure is inconclusive, "
            "it does not show the limit off;") in rep                                               # round 129 (D3dx-1, -5): the window, from A6C
    assert ("; at the transient's 1 rule ((0.25, 2·v_hs): resolvable at every order's smallest K) the power fits worsen with more terms, B > 0, B_K flat "
            "(within 0.001 of B outside the fit), and the transient fit is within 6.6e-16 (rounded up) of 1/t(v), inside its next order's change (or its "
            "input bound, if larger) 2.5e-15 (that, the resolution; the transient's own q^{2K} term unmodelled) — the transient's own criterion, not "
            "the orders' decrease: its carried orders do NOT decrease in magnitude — numerical, consistent with an expansion in integer powers of "
            "1/K from K⁻² (no K⁻¹ term) and the transient (a non-integer power can pass this check), not a proof;") in rep   # round 129 (D3dx-4)
    assert "(A6X extends 4 of the 9 fast-branch rules without a frame)." in rep
    assert ("B_K = r_K/q^K at K = 60, 80, outside the fit, within 1.2e-04 of B — at the fit's own K it is B by construction —" in rep
            and "its next order (with K⁻⁵, from K = 80) +3.1e-15" in rep)
    assert WM.a6x_resolvable(0.850630) and not WM.a6x_resolvable(0.605) and not WM.a6x_resolvable(0.0)   # round 118 (D3dm-4)
    assert ("at the 3 proved rules the deviations shrink by 20.7–26.2 and 10–12.2 times" in rep
            and "a pure K⁻⁴ and K⁻³ remainder gives 21.11 and 10, 2.533 and 2" in rep)
    assert ("## D. The tail of the exact series (the sharp constant: proved by D‴.61 (iv) at every t < 1 — every σ*/σ —, " + WM.D61_STATUS
            + "; before VA for t < 0.9900 with V, for t < 0.8052 with the branch and bound alone)") in rep    # round 123 (X12): each tier its r*
    assert "q = (h_other/t(v))² = 0.850630 per kernel" in rep and "falling by 0.85219, 0.85069, 0.85062, 0.85026 per kernel" in rep
    assert rep.count("A6 beyond K = 160 on the fast branch (A6X") == 1 and "calibrated by it" not in rep
    assert bs[0]['R1'] == GP.DOM_BB[-1] and all(a['R2'] == b['R1'] for a, b in zip(bs, bs[1:]))   # round 113: V contiguous with DB
    # round 113 (D3dh-3): plants through the 26-digit errors as well as the doubles
    from decimal import Decimal
    pl = lambda xs, key, d: [dict(x, rows=[dict(y, ratio=y['ratio'] * (1 + d), e0=str(Decimal(y['e0']) * (1 + Decimal(repr(d)))))
                                           for y in x['rows']]) if (x['r'], x['name']) == key else x for x in xs]
    s_ = '\n'.join(WM.a6x_lines(pl(raw['A6X'], (0.25, '2·v_hs'), 5e-8), pl(raw['A6'], (0.25, '2·v_hs'), 5e-8)))
    assert "every deviation shrinks as the window moves out (a necessary condition)" in s_
    assert "the transient fit does NOT account for the residual" in s_
    s_ = '\n'.join(WM.a6x_lines(pl(raw['A6X'], (0.5, 'moment matching'), 1e-15), pl(raw['A6'], (0.5, 'moment matching'), 1e-15)))
    assert ("(τ/σ, rule) = (1, moment matching), (2, moment matching), (2, 2·v_hs): decreasing in magnitude through eight terms" in s_   # round 123:
            and "; (0.5, moment matching), (1, Mercer frame), (2, Mercer frame): NOT decreasing in magnitude through eight terms" in s_)  # its own
    assert "; at (1, Mercer frame), (2, Mercer frame) the other branch's transient, resolvable" in s_   # round 129: not at q = 0 (moment matching)
    s0 = '\n'.join(WM.a6x_lines(raw['A6X'], raw['A6']))
    assert s0.count("NOT") == 2 and "(1, Mercer frame), (2, Mercer frame): NOT decreasing" in s0 and "do NOT decrease in magnitude" in s0
    assert ("at (1, Mercer frame), (2, Mercer frame) the other branch's transient, resolvable in the 26-digit errors at some orders' smallest K, "
            "lies below the order's own remainder at each: not carried, it does not account for the stall (a negative result)") in s0
    assert "the window: on K" not in s0 and "the window: on K ≥ 160 the orders decrease (A6C)" in rep  # round 129: the window only from A6C
    key = (0.25, '2·v_hs')                                               # round 118: B_K not flat (K = 60 moved alone)
    a6p = [dict(x, rows=[dict(y, e0=str(Decimal(y['e0']) * (1 + Decimal('1e-5')))) if y['K'] == 60 else y for y in x['rows']])
           if (x['r'], x['name']) == key else x for x in raw['A6']]
    s_ = '\n'.join(WM.a6x_lines(raw['A6X'], a6p))
    assert "; at (0.25, 2·v_hs) the transient fit does NOT account for the residual" in s_ and ok4 in s_


def test_round107_a6x_verdict_plants():
    """Rounds 107, 113, 118, 123 (D3dh-1, -3; D3dm-1; D3dr-3): on mock rules with exact 40-digit errors (ratio 4(1 + 1/(K + 1)²), every
    order of 1/K from K⁻² present) the orders decrease and agree with 1/t(v) at their last change, and a limit offset of 10⁻¹² planted in
    the errors is caught at that rule alone (the orders still decrease, toward the offset: inconclusive); a frame rule
    with t(v) = 0.995 ≥ D61_RSTAR_V is proved since VA (round 121: proved iff a frame), and the coverage counts the same population
    (fast-branch rules without a frame)."""
    from fractions import Fraction
    Ks = list(WM.A6_K) + list(WM.A6X_K)
    def mk(rate=0.25, off=0):
        out = []
        for K in Ks:
            r = (1 / Fraction(rate).limit_denominator(10 ** 6)) * (1 + Fraction(1, (K + 1) ** 2)) * (1 + Fraction(off))
            out.append(dict(K=int(K), prec=100, e0=str(r.numerator * 10 ** 40 // r.denominator) + 'e-40', e1='1',
                            ratio=float(r)))
        return out
    n6 = len(WM.A6_K)
    pr = dict(r=1.0, name='proved', v=0.8, v_hs=2.0, rate=0.25, arg=2, t=0.17)          # ρ = 0.5: a frame, t(v) < r*
    un = dict(r=1.0, name='unproved', v=1.5, v_hs=3.0, rate=0.25, arg=2, t=0.17)        # ρ = 5/3: no frame
    a6 = [dict(pr, rows=mk()[:n6]), dict(un, rows=mk()[:n6])]
    a6x = [dict(pr, rows=mk()[n6:]), dict(un, rows=mk()[n6:])]
    s0 = '\n'.join(WM.a6x_lines(a6x, a6))
    assert ("rule by rule at the 2 rules without a transient resolvable at K = 160, (τ/σ, rule) = (1, proved), (1, unproved): decreasing in "
            "magnitude through eight terms, the last within ten times its last change" in s0 and "NOT" not in s0)
    assert "(A6X extends 1 of the 1 fast-branch rules without a frame)" in s0
    a6 = [dict(pr, rows=mk()[:n6]), dict(un, rows=mk(off=Fraction(1, 10 ** 12))[:n6])]
    a6x = [dict(pr, rows=mk()[n6:]), dict(un, rows=mk(off=Fraction(1, 10 ** 12))[n6:])]
    s_ = '\n'.join(WM.a6x_lines(a6x, a6))
    assert ("(τ/σ, rule) = (1, proved): decreasing in magnitude through eight terms, the last within ten times its last change" in s_
            and "; (1, unproved): decreasing in magnitude through eight terms, NOT within ten times the last change — inconclusive: the limit, "
            "the premise, or the window (the expansion not yet asymptotic at the order's smallest K);" in s_ and "a limit off" not in s_)
    assert "the known cases, the 1 rule proved: the check passes at 1 of 1 — numerical" in s_
    a6 = [dict(pr, rows=mk(off=Fraction(1, 10 ** 12))[:n6]), dict(un, rows=mk()[:n6])]             # round 123: a known case failing
    a6x = [dict(pr, rows=mk(off=Fraction(1, 10 ** 12))[n6:]), dict(un, rows=mk()[n6:])]          # consistency alone
    assert ("the known cases, the 1 rule proved: the check passes at 0 of 1 — a failure is inconclusive, it does not show the limit off"
            in '\n'.join(WM.a6x_lines(a6x, a6)))                                          # round 129 (D3dx-5): rephrased
    hi = dict(pr, rate=0.995)                                                             # a frame, t(v) beyond V's r*: VA's
    a6 = [dict(hi, rows=mk(0.995)[:n6]), dict(un, rows=mk()[:n6])]
    a6x = [dict(hi, rows=mk(0.995)[n6:]), dict(un, rows=mk()[n6:])]
    s1 = '\n'.join(WM.a6x_lines(a6x, a6))
    assert "at the proved rule" in s1 and "at the 1 rule without a frame" in s1 and "(a frame; 100 digits)" in s1 and "unproved rule" not in s1
    assert "(A6X extends 1 of the 1 fast-branch rules without a frame)" in s1            # round 113: the same population


def _a6x_mock(r_, name, fun):
    """Round 123: a mock A6/A6X record pair on a real rule's geometry (its raw A6 fields; `a6_scalings`' Decimal variance, so t(v) and
    q are the code's), the ratios fun(K) (a Decimal function) as exact 60-digit errors e₀ = fun(K), e₁ = 1, and the rule's 1/t(v) and q
    (q as the code takes it, a double)."""
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    x = next(dict(y) for y in raw['A6'] if (y['r'], y['name']) == (r_, name))
    with localcontext() as c:
        c.prec = 60
        vv = dict(WM.a6_scalings(r_))[name]
        r, = (Decimal(repr(r_)),)
        t0 = WM.rate_conj(r, vv)[0]
        hs = [(r * r / vv - 1) / (r * r / vv + 1), (vv * (1 + r * r) / (r * r) - 1) / (vv * (1 + r * r) / (r * r) + 1)]
        ho = hs[0] if abs(hs[1] - t0) < abs(hs[0] - t0) else hs[1]
        q = Decimal(float((Fraction(ho) / Fraction(t0)) ** 2))
        L = 1 / t0
    with localcontext() as c:
        c.prec = 80
        rows = [dict(K=int(K), prec=100, e0=format(fun(Decimal(K), L, q), '.60e'), e1='1', ratio=float(fun(Decimal(K), L, q)))
                for K in list(WM.A6_K) + list(WM.A6X_K)]
    n6 = len(WM.A6_K)
    return '\n'.join(WM.a6x_lines([dict(x, rows=rows[n6:])], [dict(x, rows=rows[:n6])]))


def test_round123_a6x_the_known_cases_and_the_transient_at_each_order():
    """Round 123 (D3dr-1, -2), round 129 (D3dx-1, -2): the proved Mercer frame's rules at τ/σ = 1, 2 in A6X (the known cases, computed by
    `a6_job` and stored in the raw); their orders stall, and the other branch's transient — resolvable in the 26-digit errors at the six-
    to eight-term orders' smallest K (80, 60, 40) — lies below each order's own remainder there (the power it would replace, by the next
    order with Bq^K), so it is not carried, the reason printed per rule (no B printed); on K ≥ 160, where A6C's orders converge, q^160 is
    printed; the verdict names the window there; the round-113 ratio of the proved rules' biases is the raw's. A mock on the (2, Mercer
    frame) geometry with an exact transient 0.24q^K, q = 0.604938 (resolvable at K = 80, 60, 40, not at 160), passes only through the
    transient in those orders (the positive branch; B = 0.24 recovered, in place of the top power on each order's own K); on the outlier's
    geometry an exact transient passes the transient's check at γ = 0 and fails at γ = 0.01 and with B < 0."""
    import exp_cohesion_leak as E
    from decimal import Decimal
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    mer = {x['r']: x for x in raw['A6X'] if x['name'] == 'Mercer frame'}
    assert sorted(mer) == [1.0, 2.0] and all([y['K'] for y in x['rows']] == list(WM.A6X_K) for x in mer.values())
    assert abs(mer[2.0]['rows'][-1]['ratio'] * 3 / 7 - 1) < 1e-6 and mer[2.0]['rows'][-1]['prec'] == 328
    md = open(os.path.join(WM.OUT, 'width_mismatch.md'), encoding='utf-8').read()
    line = '\n'.join(WM.a6x_lines(raw['A6X'], raw['A6'], raw['A6C'])).strip()   # computed from the raw (the stored md holds it)
    assert line.startswith("A6 beyond K = 160 on the fast branch (A6X") and '\n' + line + '\n' in md
    assert ("τ/σ = 2, Mercer frame (a frame; 328 digits): -2.10e-09 → -9.69e-11 and -1.73e-07 → -1.62e-08, the excess times K² at K = 160, 200, "
            "240, 320 -0.0455, -0.0452, -0.0450, -0.0447; by successive orders (four … eight terms) +1.2e-12, +2.5e-13, -1.3e-13, +7.9e-14, -1.4e-13; "
            "the other branch's transient Bq^K, q = (h_other/t(v))² = 0.604938 per kernel, not carried in the six, seven and eight-term orders, "
            "where q^K > 10⁻²⁰ at their smallest K (80, 60, 40): q^K 3.4e-18, 8.0e-14, 1.9e-09 there is below the order's own remainder, its next "
            "change (the power the order would drop for Bq^K, at that K, by the next order with Bq^K) ") in line
    assert ("τ/σ = 1, Mercer frame (a frame; 550 digits): +2.55e-08 → +9.76e-10 and +7.55e-07 → +6.20e-08, the excess times K² at K = 160, 200, "
            "240, 320 +0.0631, +0.0619, +0.0611, +0.0602; by successive orders (four … eight terms) -7.0e-11, +8.5e-12, -2.2e-12, +1.2e-12, -1.7e-12; "
            "the other branch's transient Bq^K, q = (h_other/t(v))² = 0.674644 per kernel, not carried in the six, seven and eight-term orders, "
            "where q^K > 10⁻²⁰ at their smallest K (80, 60, 40): q^K 2.1e-14, 5.6e-11, 1.5e-07 there is below") in line
    assert "; on K ≥ 160, where the orders converge (A6C), q^160 = 4.5e-28;" in line and "; on K ≥ 160, where the orders converge (A6C), q^160 = 1.2e-35" in line
    assert ("at (1, Mercer frame), (2, Mercer frame) the other branch's transient, resolvable in the 26-digit errors at some orders' smallest K, "
            "lies below the order's own remainder at each: not carried, it does not account for the stall, and where the orders converge, on "
            "K ≥ 160 (A6C), q^160 is at most 4.5e-28 (a negative result)") in line
    assert "(B 1" not in line and "changing sign" not in line and "second branch" not in line and "says nothing" not in line   # D3dx-2, -5
    assert ("the proved rules' biases differ (by 29.0 times at τ/σ = 2: -3.35e-12 and -9.69e-11, the first fit through K = 200, 240, 320)"
            in line and "here at most 6.8e-24" in line and "A6X, referee rounds 101, 107, 113, 118, 123, 129" in line)
    assert "not reproduced" not in line and "stall near" not in line and "25-digit" not in md and "a limit off" not in line
    P = lambda K, L: L * (1 - Decimal('0.09') / K ** 2 + Decimal('0.01') / K ** 3 + Decimal('0.001') / K ** 4)
    assert WM.a6x_resolvable(0.6, 40) and WM.a6x_resolvable(0.6, 80) and not WM.a6x_resolvable(0.6, 120)
    s = _a6x_mock(2.0, 'Mercer frame', lambda K, L, q: P(K, L) + Decimal('0.24') * q ** int(K))
    assert ("(τ/σ, rule) = (2, Mercer frame): decreasing in magnitude through eight terms, the last within ten times its last change" in s
            and "at (2, Mercer frame) the orders decrease only with the other branch's transient carried: the power fits' stall is that "
            "unmodelled transient (a positive result)" in s and "the orders still do NOT decrease" not in s)
    assert ("carried in the six, seven and eight-term orders (at their smallest K, 80, 60, 40, q^K 3.4e-18, 8.0e-14, 1.9e-09 above the order's own "
            "remainder there" in s and "on the order's own K): " in s and " (B 0.24, 0.24, 0.24)" in s)          # round 129 (D3dx-2)
    s = _a6x_mock(2.0, 'Mercer frame', lambda K, L, q: P(K, L))                                # no transient: exact, at the floor
    assert "(2, Mercer frame): decreasing in magnitude through eight terms, the last within ten times" in s and "positive result" not in s
    tr = lambda g, B='0.2383': (lambda K, L, q: P(K, L) + Decimal(B) * q ** int(K) * (K / 140) ** Decimal(g))
    s = _a6x_mock(0.25, '2·v_hs', tr('0'))
    assert "at the transient's 1 rule ((0.25, 2·v_hs): resolvable at every order's smallest K) the power fits worsen with more terms, B > 0" in s
    for s in (_a6x_mock(0.25, '2·v_hs', tr('0.01')), _a6x_mock(0.25, '2·v_hs', tr('0', '-0.2383'))):
        assert "; at (0.25, 2·v_hs) the transient fit does NOT account for the residual" in s and "B > 0" not in s


def test_round123_a6x_the_checks_at_their_edges():
    """Round 123 (D3dr-3, -4, -7): `a6x_decreasing` checks the last pair, strictly, an order at its input bound converged;
    `a6x_resolution` is ten times the last change (not once, not 100 times, not the change two back) or the bound; `a6x_transient_ok`
    needs B > 0, flatness below 10⁻³ (2·10⁻³ fails), the deviation within the next order's change or the bound, and the power fits
    worsening; the trigger q^K > 10⁻²⁰ at its K (0.70 and 0.745 not at K = 160, 0.80 yes; the edge strict). On the raw: e₀ at K = 60
    moved by 10⁻⁸ (B_K 1.9·10⁻³ off B) fails flatness; the outlier's offsets within its printed resolution pass, beyond it fail. A
    moment-matching rule at a non-dyadic τ/σ, ρ = 1 to rounding, stays frameless; `d61_step` refuses t ≥ 1."""
    import pytest
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    z = [0.0] * 5
    assert WM.a6x_decreasing([4.0, 3.0, 2.0, 1.0], z[:4]) and not WM.a6x_decreasing([4.0, 3.0, 1.0, 2.0], z[:4])   # the last pair
    assert not WM.a6x_decreasing([4.0, 3.0, 2.0, -2.0], z[:4]) and not WM.a6x_decreasing([4.0, 3.0, 3.0, 1.0], z[:4])
    assert WM.a6x_decreasing([4.0, 3.0, 1e-24, -2e-24], [0, 0, 1e-23, 1e-23]) and not WM.a6x_decreasing([4.0, 3.0, 1e-24, -2e-23], [0, 0, 1e-23, 1e-23])
    assert WM.A6X_CONS == 10 and WM.a6x_resolution([3.0, 1.0, 0.7], z) == pytest.approx(3.0)         # 0.7 within 10 × 0.3, not 1 ×
    assert WM.a6x_resolution([3.0, 1.0, 0.95], z) == pytest.approx(0.5) and WM.a6x_resolution([5.0, 1.0, 1.0], [0, 0, 1e-23]) == 1e-23
    assert WM.a6x_resolution([2.0, 1.01, 1.0], z) == pytest.approx(0.1)                               # the last change, not two back
    T = dict(B=0.2, flat=5e-4, dev=1e-15, dev5=3e-15, bnd=1e-23)
    po = [1e-13, 1e-12, 1e-11]
    assert WM.A6X_FLAT == 1e-3 and WM.a6x_transient_ok(T, po) and not WM.a6x_transient_ok(dict(T, B=-0.2), po)
    assert not WM.a6x_transient_ok(dict(T, flat=2e-3), po) and not WM.a6x_transient_ok(dict(T, dev=2.5e-15), po)
    assert WM.a6x_transient_ok(dict(T, dev=1e-24, dev5=1e-24 + 1e-60), po) and not WM.a6x_transient_ok(T, [1e-13, 1e-12, 1e-12])
    assert not WM.a6x_transient_ok(T, [1e-11, 1e-12, 1e-13])
    assert not WM.a6x_resolvable(0.70) and not WM.a6x_resolvable(0.745) and WM.a6x_resolvable(0.80) and WM.a6x_resolvable(0.70, 120)
    assert not WM.a6x_resolvable(1e-20, 1) and WM.a6x_resolvable(1.01e-20, 1) and not WM.a6x_resolvable(0.0, 1)
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    key = (0.25, '2·v_hs')
    with localcontext() as c:
        c.prec = 60
        a6p = [dict(x, rows=[dict(y, e0=str(Decimal(y['e0']) * (1 + Decimal('1e-8')))) if y['K'] == 60 else y for y in x['rows']])
               if (x['r'], x['name']) == key else x for x in raw['A6']]
        pl = lambda xs, d: [dict(x, rows=[dict(y, ratio=y['ratio'] * (1 + d), e0=str(Decimal(y['e0']) * (1 + Decimal(repr(d)))))
                                          for y in x['rows']]) if (x['r'], x['name']) == key else x for x in xs]
    s = '\n'.join(WM.a6x_lines(raw['A6X'], a6p))
    assert "B_K = r_K/q^K at K = 60, 80, outside the fit, within 1.9e-03 of B" in s and "the transient fit does NOT account" in s
    for d, ok in ((-2e-15, True), (1e-15, True), (2e-15, False), (-3.5e-15, False)):
        s = '\n'.join(WM.a6x_lines(pl(raw['A6X'], d), pl(raw['A6'], d)))
        assert ("the transient fit does NOT account" in s) != ok and ("at the transient's 1 rule" in s) == ok
    x = next(dict(y) for y in raw['A6'] if (y['r'], y['name']) == (2.0, 'moment matching'))
    with localcontext() as c:                                           # τ/σ = 1.1: ρ = 1 − 2.2·10⁻¹⁶ in doubles
        c.prec = 60
        v = float(Decimal('1.1') ** 2)
    assert abs(1 + 2 / 1.1 ** 2 - 2 / v) < 1
    m = dict(x, r=1.1, v=v)
    rows = lambda Ks: [dict(K=int(K), prec=100, e0=str(4 + 1 / K ** 2), e1='1', ratio=4 + 1 / K ** 2) for K in Ks]
    s = '\n'.join(WM.a6x_lines([dict(m, rows=rows(WM.A6X_K))], [dict(m, rows=rows(WM.A6_K))]))
    assert "τ/σ = 1.1, moment matching (no frame; 100 digits)" in s and "(A6X extends 1 of the 1 fast-branch rules without a frame)" in s
    with pytest.raises(AssertionError):
        WM.d61_step(1.0)
    assert WM.d61_step(0.9999) == ('VA', 1)


def test_round123_part_a6x_rules_merged(monkeypatch, tmp_path):
    """Round 123: `--rules i,j` runs A6X for those rules of A6X_RULES only and merges them into the raw's A6X, every record in
    A6X_RULES' order (the proved Mercer frame's rules added beside the stored five)."""
    import exp_cohesion_leak as E
    import multiprocessing
    old = [dict(r=r_, name=WM.a6_scalings(r_)[i][0], x=0) for r_, i in WM.A6X_RULES[:5]]
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), A6X=old)
    calls, saved = [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'a6_job', lambda a: (calls.append(a), dict(r=a[0], name=WM.a6_scalings(a[0])[a[1]][0], x=1))[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))

    class FakePool:
        def __init__(self, n): pass
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def map(self, fn, jobs, chunksize=1): return [fn(j) for j in jobs]
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'A6X', '--rules', '6,5'])
    WM.main()
    assert calls == [(1.0, 3, WM.A6X_K), (2.0, 3, WM.A6X_K)] and WM.A6X_RULES[5:] == ((1.0, 3), (2.0, 3))
    got = saved[-1]['A6X']
    assert [(x['r'], x['name']) for x in got] == [(r_, WM.a6_scalings(r_)[i][0]) for r_, i in WM.A6X_RULES]
    assert [x['x'] for x in got] == [0] * 5 + [1] * 2
    assert WM.part_A6X(got[:2] + got[3:], {2}) == [dict(x, x=1) if k == 2 else x for k, x in enumerate(got)]
    assert WM.part_A6X(got[::-1], {0}) == [dict(x, x=1) if k == 0 else x for k, x in enumerate(got)]      # the rules' order, not old's


def test_round129_a6c_arb_errors_nodes_and_job(monkeypatch):
    """Referee round 129 (D3dx-1): `a6c_nodes` is the Gauss–Hermite rule in Arb (symmetric, the middle node of an odd rule exactly 0, the
    weights summing to 1, exact on z^{2K−2}); `a6c_error` — independent of A6's decimal path — reproduces `a6_error` at small K to the
    latter's 60-digit margin; `a6c_str` rounds the midpoint once; `a6c_job` stores A6C_DIGITS-digit strings at every K of its range with
    the bits per K (⌈2K·log₂(1/t(v))⌉ + A6C_BITS) and restores flint's precision; two stored A6C errors recomputed are the same strings;
    `part_A6C` maps A6C_RULES in order."""
    import math
    import flint
    from flint import arb
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    p0 = flint.ctx.prec
    zs, ws = WM.a6c_nodes(7, 300)
    assert zs[3] == 0 and all(abs(float((a + b).mid())) < 1e-80 for a, b in zip(zs, zs[::-1])) and abs(float((sum(ws) - 1).mid())) < 1e-80
    assert abs(float((sum(w * z ** 12 for w, z in zip(ws, zs)) - 10395).mid())) < 1e-70       # E Z¹² = 11!!, degree 12 ≤ 2K − 1
    flint.ctx.prec = p0
    for r_, i, K in ((2.0, 3, 20), (0.25, 8, 21), (1.0, 3, 19)):
        name, v = WM.a6_scalings(r_)[i]
        rate = WM.rate_conj(r_, float(v))[0]
        prec = int(math.ceil(2 * K * math.log10(1 / rate))) + 60
        s, rr = WM.a6c_str(WM.a6c_error(r_, i, K, int(math.ceil(2 * K * math.log2(1 / rate))) + WM.A6C_BITS))
        flint.ctx.prec = p0
        with localcontext() as c:
            c.prec = prec
            e = WM.a6_error(r_, WM.a6_scalings(r_)[i][1], K, prec).sqrt()
            assert abs(Decimal(s) / e - 1) < Decimal('1e-55') and rr < 1e-100
    flint.ctx.prec = 300
    assert WM.a6c_str(arb(1) / 3, 5)[0] == '3.3333e-1' and WM.a6c_str(arb(2) / 3, 5)[0] == '6.6667e-1'
    flint.ctx.prec = p0
    monkeypatch.setattr(WM, 'A6C_K', (20, 21))
    out = WM.a6c_job((2.0, 3))
    assert flint.ctx.prec == p0 and (out['K0'], out['K1'], out['name'], out['i']) == (19, 21, 'Mercer frame', 3) and len(out['e']) == 3
    assert out['bits'] == [int(math.ceil(2 * K * math.log2(1 / out['rate']))) + 450 for K in (19, 20, 21)] and WM.A6C_BITS == 450
    assert all(len(Decimal(s).as_tuple().digits) == WM.A6C_DIGITS == 70 for s in out['e']) and out['rad'] < 1e-75
    monkeypatch.undo()
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))       # the stored strings recomputed (the proved control, K = 159, 160)
    x = next(x for x in raw['A6C'] if (x['r'], x['name']) == (2.0, '2·v_hs'))
    for K in (159, 160):
        assert WM.a6c_str(WM.a6c_error(2.0, 8, K, x['bits'][K - x['K0']]))[0] == x['e'][K - x['K0']]
    flint.ctx.prec = p0
    seen = []

    class P:
        def __init__(self, n): pass
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def map(self, f, args, chunksize=1): seen.extend(args); return ['ok']
    import multiprocessing
    monkeypatch.setattr(multiprocessing, 'Pool', P)
    assert WM.part_A6C() == ['ok'] and seen == list(WM.A6C_RULES)


def test_round129_a6c_on_the_raw_and_its_lines(monkeypatch):
    """Round 129 (D3dx-1, -4): the raw's A6C — the rules of A6C_RULES at every K = 159 … 400, A6C_DIGITS digits, the balls' radii far below
    them — agrees with the stored 26-digit errors at the K they share, within their rounding; `a6c_analysis`: at the Mercer frame's rules
    (A6X's failing known cases) the orders on K ≥ 160 decrease to the referee's 7.3e-24 and 1.1e-21 at sixteen terms, the residual floor
    falls with N, the excess's coefficients grow above 10 times per order (the last four ratios) against below 7 at the controls, and c₂
    fitted agrees with `a6x_c2`; at the outlier, with q^K(B₀ + B₁/K + B₂/K²), B₀ = 0.238316 > 0 and B₁/B₀ ≈ −0.005; the md prints A6C's
    line (`a6c_lines`) and A6X's window verdict; a planted oscillation at one rule's errors removes the window there alone."""
    import exp_cohesion_leak as E
    from decimal import Decimal
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    C = list(raw['A6C'])
    assert [(x['r'], x['name']) for x in C] == [(r_, WM.a6_scalings(r_)[i][0]) for r_, i in WM.A6C_RULES]
    assert all((x['K0'], x['K1'], len(x['e'])) == (159, 400, 242) and x['rad'] < 1e-85 for x in C) and WM.A6C_K == (160, 400)
    z = {(x['r'], x['name']): WM.a6c_analysis(x) for x in C}
    m2, m1 = z[(2.0, 'Mercer frame')], z[(1.0, 'Mercer frame')]
    assert 7.2e-24 < abs(m2['orders'][-1]) < 7.3e-24 and 1.0e-21 < abs(m1['orders'][-1]) < 1.1e-21 and len(m2['orders']) == 13
    assert all(w['mono'] and w['cons'] for w in (m2, m1)) and m2['floor'][-1] < 1e-40 and m1['floor'][-1] < 3e-36
    assert all(a > b for w in (m2, m1) for a, b in zip(w['floor'], w['floor'][1:]))
    assert m2['g4'][0] > 20 and m1['g4'][0] > 12 and all(z[k]['g4'][1] < 7 for k in ((2.0, '2·v_hs'), (2.0, 'moment matching')))
    assert [k for k, _, _ in m2['coef']] == list(range(2, 12)) and [k for k, _, _ in m1['coef']] == list(range(2, 9))
    assert m2['c2rel'] < 1e-20 and m1['c2rel'] < 1e-17 and abs(m2['coef'][0][1] + 0.044) < 1e-15 and abs(m2['k1']) < 1e-17
    tr = z[(0.25, '2·v_hs')]
    assert tr['tr'] is not None and abs(tr['tr']['J'][2]['B0'][-1] - 0.238316) < 1e-6 and -0.0051 < tr['tr']['J'][2]['b10'] < -0.0049
    assert tr['tr']['J'][1]['B0'][-1] > 0 and tr['c2rel'] < 1e-12 and tr['qK0'] > 1e-20
    md = open(os.path.join(WM.OUT, 'width_mismatch.md'), encoding='utf-8').read()
    lc = '\n'.join(WM.a6c_lines(raw['A6C'], raw['A6'], raw['A6X'])).strip()
    assert lc.startswith("A6X's rules at consecutive K (A6C, referee round 129)") and '\n' + lc + '\n' in md
    assert ("within 2.4e-26 (relative) of A6's and A6X's 26-digit errors at K − 1 and K for K = 160, 200, 240, 320 (their rounding at most "
            "5.0e-26)") in lc and "stored to 70 digits, the balls' relative radii at most 3.5e-86" in lc and "637–2377" in lc
    assert ("(2, Mercer frame; a frame; q = 0.604938, q^160 = 1.2e-35): the orders +7.9e-13, +6.4e-14, -8.0e-15, +7.6e-16, -7.5e-17, +7.6e-18, "
            "-8.4e-19, +1.0e-19, -1.3e-20, +1.8e-21, -2.7e-22, +4.3e-23, -7.3e-24 — decreasing in magnitude, the last within ten times its last "
            "change (5.1e-22; the input bound at most 1.4e-60); the residual floor 6.9e-18, 9.8e-24, 5.3e-29, 6.8e-34, 9.5e-38, 9.4e-41; the K⁻¹ "
            "coefficient +3.1e-18; c_2 … c_11 -0.044, -0.235536, -1.09167, -3.56048, +28.2686, +1174.07, +28436.5, +644380, +1.48467e+07, "
            "+3.56822e+08, the ratios c_k/c_{k−1} +5.4, +4.6, +3.3, -7.9, +41.5, +24.2, +22.7, +23.0, +24.0") in lc
    assert ("with B₂/K² too +7.8e-17, -1.1e-17, -3.1e-20, +2.1e-22, +2.4e-24, +1.4e-26" in lc and "B₁/B₀ -0.0050: a small 1/K term in B) — the "
            "transient's own next term, q^{2K} = 7.9e-29 at K = 200, unmodelled" in lc and "B₁/B₀ -0.0063 at sixteen terms" in lc)
    md = open(os.path.join(WM.OUT, 'width_mismatch.md'), encoding='utf-8').read()
    assert "within 1.2e-04 of B — at the fit's own K it is B by construction —: no power of K beside q^K to that resolution)" in md   # D3dx-4
    assert md.count(WM.D61_STATUS) == 5 and md.count("Appendix A.0 assumed correct") == 5 and "Arb's correctness assumed" not in md and "VA verified by referee round 121" not in md
    from decimal import localcontext
    with localcontext() as c_:                                          # the oscillation planted at A6C's own digits
        c_.prec = 100
        pl = [dict(x, e=[format(Decimal(s) * (1 + Decimal((-1) ** K) * Decimal('1e-15')), '.69e') for K, s in zip(range(x['K0'], x['K1'] + 1), x['e'])])
              if (x['r'], x['name']) == (2.0, 'Mercer frame') else x for x in C]
    assert not WM.a6c_analysis(pl[0])['mono'] and WM.a6c_analysis(pl[0])['g4'] is None
    s = '\n'.join(WM.a6x_lines(raw['A6X'], raw['A6'], pl))
    assert "at (1, Mercer frame) the window: on K ≥ 160 the orders decrease (A6C), to -1.1e-21 at sixteen terms" in s and "(2, Mercer frame) the window" not in s
    assert "q^160 = 4.5e-28" in s and "q^160 = 1.2e-35" not in s and "q^160 is at most 4.5e-28" in s   # q^160 only where A6C's orders converge
    real = WM.a6c_analysis                                              # orders not decreasing, coefficients kept: still no window there
    monkeypatch.setattr(WM, 'a6c_analysis', lambda x: dict(real(x), mono=False) if (x['r'], x['name']) == (2.0, 'Mercer frame') else real(x))
    s = '\n'.join(WM.a6x_lines(raw['A6X'], raw['A6'], raw['A6C']))
    assert WM.a6c_analysis(C[0])['g4'] is not None and "at (1, Mercer frame) the window: on K ≥ 160" in s and "(2, Mercer frame) the window" not in s


def test_round129_the_next_term_c2():
    """Round 129 (D3dx-3): `a6x_c2`, c₂ = s₁/2 in closed form (formal at next order), exact on rational inputs — (2, Mercer) −11/250,
    (½, MM) −13/250, (1, MM) −11/128, (2, MM) −187/4000 — and at 60 digits the referee's table: (1, Mercer) 0.0575470734670659269,
    (2, 2·v_hs) −0.051243224484370181, the outlier −0.0924369845420, the slow branch (2, 0.8·v_hs) 5.7308947080 and (2, 0.5·v_hs)
    21.988001779 (ϱ = −t(v)); None at the cusp; `a6x_branches` picks the larger argument, a tie as the slow branch. The md prints c₂ per
    rule beside its fit and the check with c₂ fixed, which passes at all three known cases."""
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    F = Fraction
    assert WM.a6x_c2(2.0, F(2)) == F(-11, 250) and WM.a6x_c2(0.5, F(1, 4)) == F(-13, 250)
    assert WM.a6x_c2(1.0, F(1)) == F(-11, 128) and WM.a6x_c2(2.0, F(4)) == F(-187, 4000)
    assert WM.a6x_branches(2.0, F(2)) == (F(3, 7), F(1, 3), True) and WM.a6x_branches(0.5, F(1, 4)) == (F(1, 9), F(0), True)
    with localcontext() as c:
        c.prec = 60
        g = lambda r_, i: WM.a6x_c2(r_, WM.a6_scalings(r_)[i][1])
        for r_, i, ref, tol in ((1.0, 3, '0.0575470734670659269', 1e-12), (2.0, 8, '-0.051243224484370181', 1e-12), (0.25, 8, '-0.0924369845420', 1e-12),
                                (2.0, 5, '5.7308947080', 1e-9), (2.0, 4, '21.988001779', 1e-9), (2.0, 3, '-0.044', 1e-15), (0.5, 0, '-0.052', 1e-15)):
            assert abs(float(g(r_, i) / Decimal(ref) - 1)) < tol, (r_, i)
        assert not WM.a6x_branches(2.0, WM.a6_scalings(2.0)[4][1])[2] and g(1.0, 6) is None          # 0.5·v_hs slow; 1·v_hs the cusp
    assert WM.a6x_c2(0.75, F(9, 20)) is None and WM.a6x_c2(0.75, F(9, 20) * (1 + F(1, 10 ** 30))) is None   # σ* = 5/4: the cusp v_hs = 9/20,
    assert WM.a6x_c2(0.75, F(9, 20) * (1 + F(1, 10 ** 12))) is not None                                   # exactly and to 10⁻³⁰ (√q to 10⁻²⁰)
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    s = '\n'.join(WM.a6x_lines(raw['A6X'], raw['A6'], raw['A6C']))
    assert ("(2, Mercer frame) -0.044 and -4.4000000000000e-02 (A6C), 9.4e-22" in s and "(0.5, moment matching) -0.052 and -5.1999999999997e-02 "
            "(the stored errors), 6.1e-14" in s and "(1, Mercer frame) 0.0575470734670659269 and +5.7547073467066e-02 (A6C), 1.2e-18" in s)
    assert ("(2, Mercer frame) -5.14e-13, -4.39e-14, +1.20e-14, -4.41e-15, +2.75e-15, decreasing, the last within ten times its last change "
            "(7.2e-14) — the known cases, the 3 rules proved: with c₂ fixed the check passes at 3 of 3" in s)
    assert "formal at next order — D‴.61 proves only the leading limit, the next-order uniformity and the Laplace remainder are not proved" in s
    off = [dict(x, rows=[dict(y, e0=str(Decimal(y['e0']) * (1 + Decimal('1e-12')))) for y in x['rows']]) if x['name'] == 'Mercer frame' else x
           for x in raw['A6X']]                                         # a limit offset at the Mercer rules' A6X errors: c₂ fixed fails
    s = '\n'.join(WM.a6x_lines(off, raw['A6']))
    assert "with c₂ fixed the check passes at 1 of 3" in s
    for c2m, ok in (('-0.0441', False), ('-0.044', True)):              # a K⁻² coefficient off c₂ by 2.3e-3 is caught with c₂ fixed
        P = lambda K, L: L * (1 + Decimal(c2m) / K ** 2 + Decimal('0.01') / K ** 3 + Decimal('0.001') / K ** 4)
        s = _a6x_mock(2.0, 'Mercer frame', lambda K, L, q: P(K, L))
        assert (", decreasing, NOT within ten times its last change — the known cases, the 1 rule proved: with c₂ fixed the check passes at 0 of 1"
                in s) != ok and ("with c₂ fixed the check passes at 1 of 1" in s) == ok


def _a6x_mock_at(r_, name, fun):
    """Round 129: as `_a6x_mock`, on a geometry not in the raw (r_ any τ/σ; the record's fields from `a6_scalings` and `rate_conj`)."""
    from decimal import Decimal, localcontext
    with localcontext() as c:
        c.prec = 60
        vv = dict(WM.a6_scalings(r_))[name]
        t0, ho, fast = WM.a6x_branches(r_, vv)
        q = Decimal(float((Fraction(ho) / Fraction(t0)) ** 2))
        L = 1 / t0
        tau2 = Decimal(repr(r_)) ** 2
        x = dict(r=r_, name=name, v=float(vv), v_hs=float(vv / (tau2 / (1 + tau2).sqrt())), rate=float(t0), arg=2 if fast else 1, t=0.1)
    with localcontext() as c:
        c.prec = 80
        rows = [dict(K=int(K), prec=100, e0=format(fun(Decimal(K), L, q), '.60e'), e1='1', ratio=float(fun(Decimal(K), L, q)))
                for K in list(WM.A6_K) + list(WM.A6X_K)]
    n6 = len(WM.A6_K)
    return '\n'.join(WM.a6x_lines([dict(x, rows=rows[n6:])], [dict(x, rows=rows[:n6])])), float(q)


def test_round129_the_transient_carry_and_its_plants():
    """Round 129 (D3dx-2, -6): `a6x_carry` — q^K resolvable and above the order's own remainder (None: not carried); the carried order
    replaces the top power on the order's own K (never below its smallest K); `a6x_decreasing` takes each pair's later bound (a bound
    vector unequal across orders: the referee's Y01); a Mercer frame's rule with q ≈ 0.70 (resolvable at the five-term order's K = 120,
    not at 160: the referee's Y07 path) goes through `a6x_lines` without error, judged by its orders, not by the transient's analysis."""
    from decimal import Decimal
    assert WM.a6x_carry(0.6, 40, 1e-9) and not WM.a6x_carry(0.6, 40, 2e-9) and not WM.a6x_carry(0.6, 40, None)
    assert not WM.a6x_carry(0.6, 120, 0.0) and WM.a6x_carry(0.7, 120, 1e-19) and not WM.a6x_carry(0.7, 120, 3e-19)
    assert WM.a6x_decreasing([4.0, 3.0, 1e-24, -2e-24], [0, 0, 0, 1e-23]) and not WM.a6x_decreasing([4.0, 3.0, 1e-24, -2e-24], [0, 0, 1e-23, 0])
    P = lambda K, L: L * (1 - Decimal('0.09') / K ** 2 + Decimal('0.01') / K ** 3 + Decimal('0.001') / K ** 4)
    s, q = _a6x_mock_at(0.88, 'Mercer frame', lambda K, L, q: P(K, L) + Decimal('0.24') * q ** int(K))
    assert 0.68 < q < 0.75 and WM.a6x_resolvable(q, 120) and not WM.a6x_resolvable(q, 160)
    assert "at the transient's" not in s and "rule by rule at the 1 rule without a transient resolvable at K = 160, (τ/σ, rule) = (0.88, Mercer frame): " in s
    assert ("carried in the five, six, seven and eight-term orders (at their smallest K, 120, 80, 60, 40, q^K" in s and "(B 0.24, 0.24, 0.24, 0.24)" in s
            and "(0.88, Mercer frame): decreasing in magnitude through eight terms, the last within ten times its last change" in s
            and "a positive result" in s)
    s, _ = _a6x_mock_at(0.88, 'Mercer frame', lambda K, L, q: P(K, L) + Decimal('0.24') * q ** int(K) + L * Decimal('1e-4') / K ** 7)
    assert "at the transient's" not in s and "(0.88, Mercer frame)" in s


def test_round129_main_only_a6c_and_from_raw(monkeypatch, tmp_path):
    """Round 129: `main --only A6C` runs `part_A6C` once and stores it beside the raw's other parts, the report taking it last; `--from-raw`
    computes nothing (A6C read from the raw, `part_A6C` never called)."""
    import exp_cohesion_leak as E
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), A6X=[dict(r=1.0, name='x')], A6C=['old'])
    calls, saved, rep = [], [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'part_A6C', lambda: (calls.append(1), ['a6c'])[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: rep.append(a) or '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'A6C'])
    WM.main()
    assert calls == [1] and saved[-1]['A6C'] == ['a6c'] and [dict(x) for x in saved[-1]['A6X']] == [dict(r=1.0, name='x')]
    assert len(rep[-1]) == 13 and rep[-1][12] == ['a6c']
    monkeypatch.setattr(WM, 'part_A6C', lambda: (_ for _ in ()).throw(AssertionError('computed by --from-raw')))
    monkeypatch.setattr(sys, 'argv', ['x', '--from-raw'])
    WM.main()
    assert list(rep[-1][12]) == ['old'] and len(saved) == 1
