"""Proposition D‴.61: the Gauss–Hermite Peano kernel (`experiments/exp_gh_peano.py`)."""
import math
import os
import sys
import warnings
from fractions import Fraction
from math import factorial

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_gh_peano as GP
import exp_width_mismatch as WM


def dfact(n):
    return 1 if n <= 0 else n * dfact(n - 2)


def test_dyck_paths_are_the_gauss_hermite_moments_and_defects():
    """`dyck` gives E z^m = (m − 1)!! (odd m: 0) unbounded and Q_K z^m with heights below K; `defect` equals (E − Q_K)[z^m] computed from
    60-digit Gauss–Hermite nodes and weights (`gh_nodes`) for K ≤ 6, m ≤ 2K + 6; D(m) = 0 below 2K and D(2K) = K!."""
    from decimal import Decimal, localcontext
    assert [GP.dyck(m) for m in range(9)] == [1, 0, 1, 0, 3, 0, 15, 0, 105]
    assert GP.dyck(4, 1) == 1 and GP.dyck(6, 2) == 1 + 2 + 2 + 4 and GP.dyck(3, 2) == 0
    for K in range(1, 7):
        with localcontext() as ctx:
            ctx.prec = 60
            zs, ws = WM.gh_nodes(K, 60)
            for m in range(0, 2 * K + 7):
                q = sum((w * (z ** m if m else Decimal(1)) for z, w in zip(zs, ws)), Decimal(0))
                ex = Decimal(dfact(m - 1) if m % 2 == 0 else 0) - q
                assert abs(ex - GP.defect(m, K)) < Decimal(10) ** -40 * (1 + abs(ex)), (K, m)
                assert (GP.defect(m, K) == 0) == (m < 2 * K or m % 2 == 1), (K, m)
        assert GP.defect(2 * K, K) == factorial(K)


def test_p1_is_k2_plus_k_plus_1_and_the_polynomial_fit():
    """P_1(K) = D(2K + 2)/D(2K) = K² + K + 1 (the path count of D‴.61) for K ≤ 30; `poly_fit` recovers a known polynomial exactly and
    `poly_eval` evaluates it; `part_P` returns the polynomials with no point off them, P_1's coefficients 1, 1, 1 and P_2's leading 1/2."""
    for K in range(1, 31):
        assert Fraction(GP.defect(2 * K + 2, K), GP.defect(2 * K, K)) == K * K + K + 1
    c = GP.poly_fit([(x, Fraction(3, 2) * x ** 3 - 2 * x + Fraction(1, 7)) for x in range(4)])
    assert c == [Fraction(1, 7), -2, 0, Fraction(3, 2)] and GP.poly_eval(c, 10) == Fraction(3, 2) * 1000 - 20 + Fraction(1, 7)
    P = GP.part_P()
    assert P['d0_ok'] and all(r['bad'] == [] for r in P['rows']) and [r['j'] for r in P['rows']] == list(range(1, GP.P_J + 1))
    assert P['rows'][0]['coef'] == ['1', '1', '1'] and P['rows'][1]['coef'][-1] == '1/2' and len(P['rows'][2]['coef']) == 7
    for r in P['rows']:                                                 # the leading coefficient 1/j!: P_j ~ K^{2j}/j!
        assert Fraction(r['coef'][-1]) == Fraction(1, factorial(r['j']))


def test_cumulants_moments_and_the_limit_coefficients():
    """`cumulants` of N(0, ½)'s moments is (½, 0, 0, …) and of a ±1 coin's (0, 1, 0, −2, 0, 16); `c_limit` are the coefficients of
    ½log(1 − η²/4) (its series to η¹²); `peano_moments` from the polynomials equals the direct path ratios at K = 7; `part_M`'s variance
    is (K² + K + 1)/((K + 1)(2K + 1)) and its scaled cumulants are the definitions'."""
    g = {2 * j: Fraction(dfact(2 * j - 1), 2 ** j) for j in range(1, 7)}
    k = GP.cumulants(g, 12)
    assert k[2] == Fraction(1, 2) and all(k[i] == 0 for i in range(3, 13))
    k = GP.cumulants({2 * j: 1 for j in range(1, 4)}, 6)
    assert [k[i] for i in range(1, 7)] == [0, 1, 0, -2, 0, 16]
    for j in range(1, 7):                                               # ½log(1 − x/4) = −½Σ (x/4)^n/n, x = η²
        assert GP.c_limit(j) / factorial(2 * j) == -Fraction(1, 2 * j * 4 ** j)
    P = GP.part_P()
    polys = [[Fraction(c) for c in r['coef']] for r in P['rows']]
    K = 7
    mom = GP.peano_moments(K, polys)
    for j in range(1, GP.P_J + 1):
        direct = Fraction(factorial(2 * j) * factorial(2 * K), factorial(2 * K + 2 * j)) * Fraction(GP.defect(2 * K + 2 * j, K), factorial(K))
        assert mom[2 * j] == direct
    M = GP.part_M(P)
    assert [x['K'] for x in M] == list(GP.M_K)
    for x in M[:3]:
        K = x['K']
        assert Fraction(x['var']) == Fraction(K * K + K + 1, (K + 1) * (2 * K + 1))
        k = GP.cumulants(GP.peano_moments(K, polys), 2 * GP.P_J)
        assert x['scaled'][0] == float(K * (k[2] - Fraction(1, 2))) and x['scaled'][1] == float(k[4] * K ** 2)
        assert x['nxt'][2] == float(K * (k[6] * K ** 3 - GP.c_limit(3)))


def test_m_lines_and_p_lines():
    """`m_lines` prints the table, the relative deviations at the last K and the Richardson values 2s(2K) − s(K); `p_lines` the
    polynomials; both empty without data."""
    P = GP.part_P()
    M = [dict(K=K, var=str(Fraction(K * K + K + 1, (K + 1) * (2 * K + 1))), scaled=[float(GP.c_limit(j)) * (1 + 1 / K) for j in range(1, 7)],
              nxt=[float(GP.c_limit(j)) for j in range(1, 7)]) for K in (100, 200)]
    s = '\n'.join(GP.m_lines(M))
    assert '| K | variance | K(κ_2 − ½) | κ_4·K^2 |' in s and f"| 200 | {float(Fraction(40201, 80601)):.10f} | -0.25125 [-0.25] |" in s
    rich = [2 * b - a for a, b in zip(M[0]['scaled'], M[1]['scaled'])]
    assert "At K = 200 the scaled cumulants lie within 5.0e-03, 5.0e-03" in s
    assert "extrapolated from K = 100 and 200 (2s(2K) − s(K), the 1/K term removed) " + ', '.join(f"{v:.6g}" for v in rich) in s
    assert all(abs(v / float(GP.c_limit(j)) - 1) < 1e-12 for j, v in enumerate(rich, 1))
    t = '\n'.join(GP.p_lines(P))
    assert '| 1 | 1 + 1·K + 1·K^2 | 0 |' in t and '| 2 | 3 + 13/2·K + 3·K^2 + 2·K^3 + 1/2·K^4 | 0 |' in t
    assert GP.m_lines(None) == [] and GP.p_lines(None) == [] and GP.l_lines([]) == [] and GP.s_lines(None) == [] and GP.r_lines(None) == []


def test_the_transform_at_k1_and_l_job(monkeypatch):
    """At K = 1 (Q_1[e^{ξz}] = 1) ĝ_1(ξ)e^{−ξ²/4} is (e^{ξ²/2} − 1)/((ξ²/2)e^{ξ²/4}) on both axes (ξ² = ±s²) — `gen_ratio` agrees; `lt_limit`
    is (1 ∓ s²/4)^{1/2}; `l_job` calls `gen_ratio` at 3K + 60 digits for every s; `l_lines` divides by the limit."""
    import math
    for s in (0.5, 1.0, 1.7):
        x2 = s * s
        assert abs(WM.gen_ratio(1, s, 'real', 60) - (math.exp(x2 / 2) - 1) / (x2 / 2 * math.exp(x2 / 4))) < 1e-13
        assert abs(WM.gen_ratio(1, s, 'imag', 60) - (1 - math.exp(-x2 / 2)) / (x2 / 2 * math.exp(-x2 / 4))) < 1e-13
    assert GP.lt_limit(1.0, 'real') == np.sqrt(0.75) and GP.lt_limit(2.0, 'imag') == np.sqrt(2.0)
    calls = []
    monkeypatch.setattr(WM, 'gen_ratio', lambda K, s, ax, prec: calls.append((K, s, ax, prec)) or GP.lt_limit(s, ax) * (1 + s / K))
    out = GP.l_job((20, 'imag'))
    assert calls == [(20, s, 'imag', 120) for s in GP.L_S] and out['s'] == list(GP.L_S) and (out['K'], out['axis']) == (20, 'imag')
    L = [dict(K=K, axis=ax, s=list(GP.L_S), ratio=[GP.lt_limit(s, ax) * (1 + s / K) for s in GP.L_S]) for K in (20, 40) for ax in ('real', 'imag')]
    s = '\n'.join(GP.l_lines(L))
    assert f"| 20 | real | {1 + 0.1 / 20:.4f} | " in s and f" | {1.9 / 20:.2e} |" in s
    assert f"real: {1.9 / 20:.2e}, {1.9 / 40:.2e} (0.500)" in s and "imag: " in s
    assert f"; on s ≤ 1.5: real: {1.5 / 20:.2e}, {1.5 / 40:.2e} (0.500); imag: " in s and GP.L_SMAX == 1.5


def test_gammas_hypergeometric_and_the_predictions():
    """`gammas` at the half-squeeze frame (ρ = t, v′ = 2) are both 1/(2t) (the cusp); `hyp_half` is ₂F₁(−½, −½; 1; x) (4/π at x = 1 by
    its series' sum, 1 at 0); `log_s_pred` approaches `log_s_simple` as K grows, at the cusp and on each branch, with the factor 1 ∓ r;
    `s_prec` grows with K."""
    t = 0.2
    gx, gy = GP.gammas(t, 2.0)
    assert abs(gx - 1 / (2 * t)) < 1e-14 and abs(gy - 1 / (2 * t)) < 1e-14
    assert GP.hyp_half(0.0) == 1.0 and abs(GP.hyp_half(0.25) - (1 + 0.25 / 4 + (1 / 64) * 0.0625 + (1 / 256) * 0.25 ** 3 * 1)) < 2e-5
    assert abs(GP.hyp_half(0.999999) - 4 / np.pi) < 2e-3
    for rho, vp in ((t, 2.0), (t, 1.6), (t, 2.5)):
        d = [GP.log_s_pred(K, rho, vp, n=20000) - GP.log_s_simple(K, rho, vp) for K in (200, 1600)]
        assert abs(d[1]) < abs(d[0]) and abs(d[1]) < 0.02, (vp, d)
    gx, gy = GP.gammas(t, 2.5)
    r = 1 / (2 * gx)
    assert gx < gy and abs(GP.log_s_simple(10, t, 2.5) - (20 * np.log(r) + np.log(1 - r) - np.log(t * 2.5 * np.sqrt(2 * gx * (gy - gx))))) < 1e-12
    gx, gy = GP.gammas(t, 1.6)
    r = 1 / (2 * gy)
    assert gy < gx and abs(GP.log_s_simple(10, t, 1.6) - (20 * np.log(r) + np.log(1 + r) - np.log(t * 1.6 * np.sqrt(2 * gy * (gx - gy))))) < 1e-12
    assert GP.s_prec(20, t, 2.0) > GP.s_prec(10, t, 2.0) > 60


def test_s_frames_and_s_job(monkeypatch):
    """`s_frames` keeps the convergent frames (ρ(v′ − 1) < 1) at τ/σ = 1/2, 1, 2 — ρ = t at v′ = 2, 1.6, 2.5 and ρ = q at v′ = 1 + q, values
    pinned —; `s_job` is `frame_series` in `s_prec` digits at every K of S_K and the last with S_CHECK more (the two Decimals' relative
    difference), beside the two predictions; at K = 10 the exact series is within 6% of the prediction at the cusp of τ/σ = 1; `s_lines`
    summarises by where the mass lies, with the closed forms' range."""
    from decimal import Decimal
    fr = GP.s_frames()
    assert all(rho * (vp - 1) < 1 for _, _, rho, vp in fr) and {r_ for r_, *_ in fr} == {0.5, 1.0, 2.0} and len(fr) == 12
    assert [nm for r_, nm, *_ in fr if r_ == 1.0] == ['half-squeeze'] * 3 + ['Mercer']
    t1, q1 = 3 - 2 * np.sqrt(2), (np.sqrt(3) - 1) / (np.sqrt(3) + 1)
    assert [(nm, round(rho, 12), vp) for r_, nm, rho, vp in fr if r_ == 1.0] == [('half-squeeze', round(t1, 12), 2.0), ('half-squeeze', round(t1, 12), 1.6),
                                                                                 ('half-squeeze', round(t1, 12), 2.5), ('Mercer', round(q1, 12), 1 + q1)]
    assert GP.S_K == (10, 20, 40, 80) and GP.S_CHECK == 40 and GP.S_THETA == 4096
    monkeypatch.setattr(GP, 'S_K', (4, 10))
    job = next(a for a in fr if a[0] == 1.0 and a[3] == 2.0)
    out = GP.s_job(job)
    for y in out['rows']:
        S_, _, _ = WM.frame_series(y['K'], Decimal(repr(float(job[2]))), Decimal(2), prec=y['prec'])
        assert y['logS'] == float(S_.ln()) and y['pred'] == GP.log_s_pred(y['K'], job[2], 2.0) and y['prec'] == GP.s_prec(y['K'], job[2], 2.0)
    assert out['check'] > 60 and abs(np.exp(out['rows'][-1]['logS'] - out['rows'][-1]['pred']) - 1) < 0.06
    real = GP.s_job(next(a for a in fr if a[0] == 1.0 and a[3] == 2.5))
    s = '\n'.join(GP.s_lines([out, real]))
    assert "| 1 | half-squeeze | " in s and " | cusp (" in s and " | real (" in s and "synthetic" in s
    rr = np.exp(out['rows'][-1]['logS'] - out['rows'][-1]['pred'])
    cc = np.exp(out['rows'][-1]['logS'] - out['rows'][-1]['simple'])
    assert (f"cusp (1 frames) {rr:.5f}–{rr:.5f}, K(ratio − 1) {10 * (rr - 1):+.3f} … {10 * (rr - 1):+.3f}, over the closed form "
            f"{cc:.5f}–{cc:.5f}; real (1 frames) ") in s and f"{out['check']:.0f} |" in s and "| digits agreed |" in s
    lo, hi = sorted((out['check'], real['check']))
    assert f"; the last K's series and the same with 40 more digits agree to {lo:.0f}–{hi:.0f} digits." in s


def test_the_cusp_form_and_the_ab_parametrisation():
    """D3ck-9 (S12): at the cusp `log_s_simple` is (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/(ρv′) exactly (₂F₁ by its series, not 1); the
    (a, b) forms equal the (ρ, v′) ones bit for bit; for a rule's own frame 1/(2γx) = h(τ²/v) and 1/(2γy) = h(vσ*²/τ²) (r = t(v)), on
    both sides of v_hs and at ρ ≤ 0; `gammas_ab` is finite at ρ = 0 (the Fock scaling)."""
    t = 0.2
    K = 50
    exp = 0.5 * np.log(np.pi * K) + 2 * K * np.log(t) + np.log(GP.hyp_half(t * t)) + np.log(2 * t / (2 * t))
    assert abs(GP.log_s_simple(K, t, 2.0) - exp) < 1e-12 and GP.hyp_half(t * t) > 1.01
    for rho, vp in ((t, 2.0), (t, 1.6), (0.3, 2.5)):
        assert GP.log_s_simple(K, rho, vp) == GP.log_s_simple_ab(K, rho * vp, rho * (vp - 1))
        assert GP.log_s_pred(K, rho, vp) == GP.log_s_pred_ab(K, rho * vp, rho * (vp - 1))
    h = lambda y: (y - 1) / (y + 1)
    for r_ in (0.5, 1.0, 2.0):
        tau2, ss2 = r_ ** 2, 1 + r_ ** 2
        for vr in (0.5, 0.8, 1.25, 2.0):
            v = vr * tau2 / np.sqrt(ss2)
            rho, a, b, C = GP.rule_frame(r_, v)
            if a is None:
                continue
            gx, gy = GP.gammas_ab(a, b)
            assert abs(1 / (2 * gx) - h(tau2 / v)) < 1e-12 and abs(1 / (2 * gy) - h(v * ss2 / tau2)) < 1e-12, (r_, vr)
            assert abs(1 / (2 * min(gx, gy)) - WM.rate_conj(r_, v)[0]) < 1e-12
    rho, a, b, C = GP.rule_frame(1.0, 2.0 / 3.0)                       # the Fock scaling at τ/σ = 1: v = 2τ²/(2 + τ²)
    assert abs(rho) < 1e-15 and abs(a - 1 / 3) < 1e-15 and np.isfinite(GP.gammas_ab(a, b)).all()


def test_frame_series_any_for_every_frame():
    """`frame_series_any` is `frame_series` for ρ > 0, and for ρ < 0 (v′ < 0) and ρ = 0 (the Fock scaling's limit) it gives the exact
    squared error through the frame constant — e_K² = C·S_K to 10⁻⁴⁰ at 80 digits for τ/σ = 2, 0.5v_hs (ρ = −0.736), τ/σ = 1's Mercer-zero
    rule (ρ = −0.464) and Fock rule (ρ = 0) at K = 2, 3."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 60
        a = GP.frame_series_any(4, Decimal('0.3'), Decimal('0.6'), 60)[0]
        b = WM.frame_series(4, Decimal('0.3'), Decimal('0.6') / Decimal('0.3'), prec=60)[0]
        assert a == b
    for r_, nm in ((2.0, '0.5·v_hs'), (1.0, 'Mercer-zero'), (1.0, 'Fock')):
        for K in (2, 3):
            with localcontext() as ctx:
                ctx.prec = 80
                v = next(v for n_, v in WM.a6_scalings(r_) if n_ == nm)
                tau2 = Decimal(repr(r_)) ** 2
                rh = 1 + 2 / tau2 - 2 / v
                S_, terms, low = GP.frame_series_any(K, rh, v * (1 - rh * rh) / 2, 80)
                C = (1 + tau2).sqrt() * (1 - rh * rh).sqrt() * v / tau2
                assert abs(WM.a6_error(r_, v, K, 80) / (C * S_) - 1) < Decimal(10) ** -38 and low < Decimal(10) ** -50, (r_, nm, K)
                assert (rh < 0) == (nm != 'Fock') and len(terms) > 10


def test_rule_frame_and_r_rules():
    """`rule_frame` returns (ρ, a = ρv′, b = ρ(v′ − 1), C): the half-squeeze rule's frame ρ = t, a = 2t, C = (1 − t²)^{1/2}; the Mercer frame
    rule's ρ = q, a = q(1 + q); the frame's own variance is the rule's; none for |ρ| ≥ 1 (moment matching ρ = 1); `r_rules` splits A6's
    rules into those with a frame and a convergent series and the others; R_K, R_PREC, R_TOL pinned."""
    for r_ in (0.5, 1.0, 2.0):
        ss, bb = np.sqrt(1 + r_ ** 2), np.sqrt(1 + 2 * r_ ** 2)
        t, q = (ss - 1) / (ss + 1), (bb - 1) / (bb + 1)
        rho, a, b, C = GP.rule_frame(r_, r_ ** 2 / ss)
        assert abs(rho - t) < 1e-12 and abs(a - 2 * t) < 1e-12 and abs(b - t) < 1e-12 and abs(C - np.sqrt(1 - t * t)) < 1e-12
        v = 2 * r_ ** 2 / (bb + 1)
        rho, a, b, C = GP.rule_frame(r_, v)
        assert abs(rho - q) < 1e-12 and abs(a - q * (1 + q)) < 1e-12 and abs(2 * r_ ** 2 / (2 + (1 - rho) * r_ ** 2) - v) < 1e-12
        assert GP.rule_frame(r_, r_ ** 2) == (1.0, None, None, None)
    A6 = [dict(r=1.0, name='1·v_hs', v=1 / np.sqrt(2), v_hs=1.0, arg=1, rate=0.17, rows=[]),
          dict(r=1.0, name='moment matching', v=1.0, v_hs=np.sqrt(2), arg=2, rate=1 / 3, rows=[])]
    yes, no = GP.r_rules(A6)
    assert [y[0]['name'] for y in yes] == ['1·v_hs'] and no == [dict(r=1.0, name='moment matching', v_hs=np.sqrt(2), rho=1.0)]
    assert GP.R_K == (2, 3) and GP.R_PREC == 80 and GP.R_TOL == 0.01


def test_r_job_and_r_lines(monkeypatch):
    """`r_job` checks the frame constant at R_K in R_PREC digits and divides A6's e_K² by C·exp(per-angle prediction) and by
    C·exp(closed form) in logs; `r_lines` prints the census of frames (ρ > 0, ρ < 0, ρ = 0, none), both summaries at the last K and the
    rules off by more than R_TOL."""
    A6 = [dict(r=1.0, name='1·v_hs', v=1 / np.sqrt(2), v_hs=1.0, arg=1, rate=0.17, rows=[dict(K=K, e1=None) for K in (20, 40)]),
          dict(r=1.0, name='Mercer-zero', v=1 / np.sqrt(3), v_hs=0.816497, arg=1, rate=0.27, rows=[dict(K=K, e1=None) for K in (20, 40)]),
          dict(r=1.0, name='moment matching', v=1.0, v_hs=np.sqrt(2), arg=2, rate=1 / 3, rows=[])]
    yes, no = GP.r_rules(A6)
    for (x, rho, a, b, C), f_ in zip(yes, (1.01, 1.002)):
        x['rows'][1]['e1'] = repr(float(np.exp(0.5 * (np.log(C) + GP.log_s_pred_ab(40, a, b)))) * f_)
    monkeypatch.setattr(GP, 'R_K', (2,))
    rules = [GP.r_job(y) for y in yes]
    assert rules[0]['rows'][0]['ratio'] is None and abs(rules[0]['rows'][1]['ratio'] - 1.0201) < 1e-9 and abs(rules[1]['rows'][1]['ratio'] - 1.002 ** 2) < 1e-9
    assert all(abs(x['chk'][0]['dev']) < 1e-38 for x in rules) and rules[1]['rho'] < 0
    cl = np.exp(np.log(1.0201) + GP.log_s_pred_ab(40, yes[0][2], yes[0][3]) - GP.log_s_simple_ab(40, yes[0][2], yes[0][3]))
    assert abs(rules[0]['rows'][1]['closed'] - cl) < 1e-9
    s = '\n'.join(GP.r_lines(dict(rules=rules, none=no)))
    assert "Of A6's 3 rules 2 have a frame (1 with ρ < 0, 1 with ρ > 0) and 1 not (|ρ| ≥ 1: τ/σ = 1, moment matching (+1.000))." in s
    assert "over the per-angle prediction 1.00400–1.02010, 1 of the 2 within 0.01 (1.00400–1.00400); the others τ/σ = 1, 1·v_hs (1.02010)" in s
    assert "| 1 | Mercer-zero | 0.816497 | -0." in s and "| —, 1.02010 | " in s and GP.r_lines(None) == []


def test_the_transform_off_the_axes(monkeypatch):
    """`gen_ratio_c` is `gen_ratio` on the axes (θ = 0, π/2) and at K = 1 off them (1 − e^{−ξ²/2})/((ξ²/2)e^{−ξ²/4}); `lc_limit` is
    (1 − ξ²/(4K))^{1/2}; `lc_job`'s precision covers e^{y²/2}; `lc_lines` prints the deviations and their successive ratios; LC_K, LC_TH
    pinned."""
    import cmath, math
    for s in (0.5, 1.3):
        assert GP.gen_ratio_c(20, s, 0, 120)[0] == WM.gen_ratio(20, s, 'real', 120)
        assert abs(GP.gen_ratio_c(20, s, 4, 120)[0] - WM.gen_ratio(20, s, 'imag', 120)) < 1e-14
        for j in (1, 2, 3):
            xi = s * cmath.exp(1j * j * math.pi / 8)
            ref = (1 - cmath.exp(-xi * xi / 2)) / (xi * xi / 2 * cmath.exp(-xi * xi / 4))
            got = complex(*GP.gen_ratio_c(1, s, j, 60))
            assert abs(got - ref) < 1e-13, (s, j)
            assert abs(GP.lc_limit(s, j) - cmath.sqrt(1 - xi * xi / 4)) < 1e-15
    assert GP.LC_K == (20, 40, 80, 160) and GP.LC_TH == (1, 2, 3)
    calls = []
    monkeypatch.setattr(GP, 'gen_ratio_c', lambda K, s, j8, prec, pi=None, nodes=None: calls.append((K, s, j8, prec)) or (lambda z: (z.real, z.imag))(GP.lc_limit(s, j8) * (1 + s / K)))
    out = GP.lc_job((20, 2))
    y2 = (max(GP.L_S) * math.sqrt(20) * math.sin(math.pi / 4)) ** 2
    assert {c[3] for c in calls} == {3 * 20 + 60 + int(math.ceil(y2 / (2 * math.log(10))))} and [c[1] for c in calls] == list(GP.L_S)
    assert all(abs(d - s / 20) < 1e-12 for d, s in zip(out['dev'], GP.L_S)) and max(abs(a) for a in out['arg']) < 1e-12
    LC = [dict(K=K, j8=2, prec=1, s=list(GP.L_S), dev=[s / K for s in GP.L_S], arg=[0.0] * len(GP.L_S)) for K in (20, 40)]
    s = '\n'.join(GP.lc_lines(LC))
    assert "| 20 | 2π/8 | 5.0e-03 |" in s and f"θ = 2π/8: {1.9 / 20:.2e}, {1.9 / 40:.2e} (0.500)" in s and GP.lc_lines(None) == []


def test_report_and_main(monkeypatch, tmp_path):
    """`main --only P,M` runs `part_P` and `part_M` on P's result and saves both; `--from-raw` writes the md from the raw alone, also when
    the md is missing (O6)."""
    import exp_cohesion_leak as E
    monkeypatch.setattr(GP, 'OUT', str(tmp_path))
    monkeypatch.setattr(GP, 'part_P', lambda: dict(rows=[dict(j=1, coef=['1', '1', '1'], bad=[], kchk=3)], d0_ok=True, kchk=3))
    got = []
    monkeypatch.setattr(GP, 'part_M', lambda P: got.append(P) or [dict(K=10, var='31/66', scaled=[-0.2], nxt=[0.5])])
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'P,M'])
    GP.main()
    raw = E.load(str(tmp_path / 'gh_peano_raw.npz'))
    assert set(raw) >= {'P', 'M'} and got[0]['kchk'] == 3
    md1 = open(tmp_path / 'gh_peano.md').read()
    os.remove(tmp_path / 'gh_peano.md')
    monkeypatch.setattr(sys, 'argv', ['x', '--from-raw'])
    GP.main()
    assert open(tmp_path / 'gh_peano.md').read() == md1 and '| 1 | 1 + 1·K + 1·K^2 | 0 |' in md1
    assert GP.KEYS == ('P', 'M', 'L', 'S', 'R', 'LC', 'Q', 'D', 'V', 'DB', 'VC', 'VA', 'VA011') and GP.P_J == 6 and GP.P_KCHK == 20
    assert GP.M_K == (10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120) and GP.L_K == (20, 40, 80, 160) and GP.L_SMAX == 1.5


def test_the_md_is_the_report_of_the_raw():
    """The committed md is `report` of the committed raw, byte for byte (a report-only change fails here, not only in `--from-raw`)."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    md = open(os.path.join(GP.OUT, 'gh_peano.md'), encoding='utf-8').read()
    assert GP.report(*[raw.get(k) for k in GP.KEYS]) == md


def test_round90_survivors_pinned():
    """Referee round 90 (S15, S18, S19, S24) and this round's (LC1): `hyp_half` to 10⁻¹⁵ against the exact series (fractions); `s_prec`'s
    digits 2K·log₁₀(1/r) + 2K + 60; `s_job`'s check two different precisions (it agrees to fewer digits than the series carries);
    `lc_lines`' largest |arg| over every job."""
    from fractions import Fraction
    import math
    x = Fraction(9, 25)
    tot, term = Fraction(1), Fraction(1)
    for n in range(200):
        term *= (Fraction(2 * n - 1, 2 * n + 2)) ** 2 * x
        tot += term
    assert abs(GP.hyp_half(0.36) - float(tot)) < 1e-15
    t = 0.2
    gx, gy = GP.gammas(t, 2.0)
    assert GP.s_prec(20, t, 2.0) == int(math.ceil(40 * math.log10(1 / (1 / (2 * min(gx, gy)))))) + 40 + 60
    out = GP.s_job((1.0, 'half-squeeze', 3 - 2 * math.sqrt(2), 2.0)) if False else None
    import pytest
    mp = pytest.MonkeyPatch()
    mp.setattr(GP, 'S_K', (4,))
    out = GP.s_job((1.0, 'half-squeeze', 3 - 2 * math.sqrt(2), 2.0))
    mp.undo()
    assert 60 < out['check'] < out['rows'][-1]['prec']
    LC = [dict(K=K, j8=2, prec=1, s=[0.5], dev=[0.01], arg=[0.0 if K == 20 else 0.5]) for K in (20, 40)]
    assert "the largest |arg| 5.0e-01." in '\n'.join(GP.lc_lines(LC))


def test_the_corner_ratio_part_q(monkeypatch):
    """Part Q: `corner_ratio` is c_K(u) = Q_K[e^{uz}He_{K−1}(z)]/(u^{K−1}e^{u²/2}) — at K = 1 e^{−u²/2}, at K = 2 sinh(u)/(u e^{u²/2}) (nodes ±1),
    on and off the real axis —; `q_job` divides it by 1 − μ²; `q_lines` prints the deviations and their successive ratios; Q_K, Q_MU,
    Q_TH pinned."""
    import cmath, math
    for j in (0, 2, 4):
        e = cmath.exp(1j * j * math.pi / 8)
        u1, u2 = 0.7 * e, 0.7 * math.sqrt(2) * e
        assert abs(complex(*GP.corner_ratio(1, 0.7, j, 60)) - cmath.exp(-u1 * u1 / 2)) < 1e-14
        assert abs(complex(*GP.corner_ratio(2, 0.7, j, 60)) - cmath.sinh(u2) / (u2 * cmath.exp(u2 * u2 / 2))) < 1e-14
    assert GP.Q_K == (20, 40, 80, 160) and GP.Q_MU == (0.1, 0.3, 0.5, 0.7, 0.9) and GP.Q_TH == (0, 2, 4)
    monkeypatch.setattr(GP, 'corner_ratio', lambda K, m, j8, prec, pi=None, nodes=None, he=None:
                        (lambda z: (z.real, z.imag))((1 - (m * cmath.exp(1j * j8 * math.pi / 8)) ** 2) * (1 + m / K)))
    out = GP.q_job((20, 2))
    assert out['prec'] == 3 * 20 + 80 and all(abs(d - m / 20) < 1e-12 for d, m in zip(out['dev'], GP.Q_MU))
    Q = [dict(K=K, j8=0, prec=1, mu=list(GP.Q_MU), dev=[m / K for m in GP.Q_MU]) for K in (20, 40)]
    s = '\n'.join(GP.q_lines(Q))
    assert "| 20 | 0π/8 | 5.00e-03 |" in s and f"arg μ = 0π/8: {0.9 / 20:.2e}, {0.9 / 40:.2e} (0.500)" in s and GP.q_lines(None) == []
    assert '\\|μ\\| = 0.1' in s


class _SerialPool:
    """A stand-in for multiprocessing.Pool that maps in the calling process (so a patched job is the one called)."""
    def __init__(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def map(self, f, xs, chunksize=1):
        return [f(x) for x in xs]


def test_round91_q_job_unpatched_part_q_and_main(monkeypatch, tmp_path):
    """Referee round 91 (D3cl-5: QJ4, QJ5, QP1, QM1): `q_job` with the real `corner_ratio` at K = 2 against the closed form
    sinh(u)/(ue^{u²/2}) over 1 − μ², off the axis (a complex deviation); `part_Q`'s jobs (K, θ) over Q_K × Q_TH; `main --only Q` saves
    part_Q's output under Q and leaves the other parts."""
    import cmath, math
    import multiprocessing
    import exp_cohesion_leak as E
    out = GP.q_job((2, 2))
    e = cmath.exp(1j * 2 * math.pi / 8)
    for m, d in zip(GP.Q_MU, out['dev']):
        u = m * math.sqrt(2) * e
        assert abs(d - abs(cmath.sinh(u) / (u * cmath.exp(u * u / 2)) / (1 - (m * e) ** 2) - 1)) < 1e-12
    c = complex(*GP.corner_ratio(2, 0.9, 2, 86)) / (1 - (0.9 * e) ** 2) - 1
    assert abs(c.imag) > 1e-3 and abs(out['dev'][-1] - abs(c)) < 1e-12 and abs(out['dev'][-1] - abs(c.real)) > 1e-4
    monkeypatch.setattr(multiprocessing, 'Pool', _SerialPool)
    monkeypatch.setattr(GP, 'q_job', lambda a: ('q', a))
    assert GP.part_Q() == [('q', (K, j)) for K in GP.Q_K for j in GP.Q_TH]
    monkeypatch.setattr(GP, 'OUT', str(tmp_path))
    E.save(str(tmp_path / 'gh_peano_raw.npz'), LC=[dict(K=1, j8=1, prec=1, s=[0.5], dev=[0.25], arg=[0.0])])
    monkeypatch.setattr(GP, 'part_Q', lambda: [dict(K=20, j8=0, prec=1, mu=[0.5], dev=[0.125])])
    monkeypatch.setattr(GP, 'part_LC', lambda: 1 / 0)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'Q'])
    GP.main()
    raw = E.load(str(tmp_path / 'gh_peano_raw.npz'))
    assert list(raw['Q']) == [dict(K=20, j8=0, prec=1, mu=[0.5], dev=[0.125])] and list(raw['LC'])[0]['dev'] == [0.25]
    assert '| 20 | 0π/8 | 1.25e-01 |' in open(tmp_path / 'gh_peano.md').read()


def test_round91_the_dirichlet_form_of_the_corner_ratio():
    """D3cl-6: c_K(u) = e^{−u²/2}E e^{uT}, T the uniform-simplex average of the K nodes (Hermite–Genocchi on the nodes: E e^{uT} =
    (K − 1)!u^{1−K}Σ_k e^{uz_k}/Π_{j≠k}(z_k − z_j)), and Var T = (K − 1)/(K + 1): (c_K(u) − 1)/u² → −1/(K + 1) as u → 0."""
    from decimal import localcontext
    import cmath
    import math
    for K in (3, 5, 8):
        with localcontext() as ctx:
            ctx.prec = 80
            zs = [float(z) for z in WM.gh_nodes(K, 80)[0]]
        for m, j in ((0.4, 0), (0.6, 3), (0.8, 4)):
            u = m * np.sqrt(K) * cmath.exp(1j * j * np.pi / 8)
            dd = sum(cmath.exp(u * z) / np.prod([z - y for y in zs if y != z]) for z in zs)
            ref = cmath.exp(-u * u / 2) * factorial(K - 1) * u ** (1 - K) * dd
            assert abs(complex(*GP.corner_ratio(K, m, j, 80)) - ref) < 1e-10 * abs(ref)
        u = 1e-4
        assert abs((GP.corner_ratio(K, u / math.sqrt(K), 0, 80)[0] - 1) / u ** 2 + 1 / (K + 1)) < 1e-6
        assert abs(sum(z * z for z in zs) / (K * (K + 1)) - (K - 1) / (K + 1)) < 1e-12


def test_round91_part_d_functions():
    """Part D: `lam_plus` is the circle bound (max_θ (2a cos θ − ½cos 2θ) = M(a) + ½) and vanishes on |μ| = 1 inside |Re μ| ≤ 1;
    `dom_gap` on the axes equals its closed forms (the supremum at s = 1/|η|: log(η − 1) + (η − 1)²/2 + ½ − 2 log η + log 4 − η²/4 on
    the real axis, −¼(2 − η)³ to leading order; log(y − 1) − (y − 1)²/2 − ½ − 2 log y + log 4 + y²/4 on the imaginary, −½(2 − y)²),
    −∞ for |η| ≤ 1; `dom_radius` is the crossing; `c_need` and `r_star` invert each other; `alias_rate` vanishes at 1.8589."""
    import math
    th = np.linspace(-np.pi, np.pi, 200001)
    for a in (0.0, 0.3, 0.9, 1.0, 1.3, 2.5):
        M = a * a if a <= 1 else 2 * a - 1
        assert abs(np.max(2 * a * np.cos(th) - 0.5 * np.cos(2 * th)) - (M + 0.5)) < 1e-9
        b = 0.7
        mu = a + 1j * b
        assert abs(GP.lam_plus(mu) - (M - 0.5 - np.log(abs(mu)) - (mu * mu).real / 2)) < 1e-15
    assert all(abs(GP.lam_plus(np.exp(1j * t))) < 1e-15 for t in np.linspace(0, np.pi, 9))
    assert GP.dom_gap(1.0) == -np.inf and GP.dom_gap(0.9j) == -np.inf
    for e in (1.3, 1.7, 1.95, 1.999):
        real = math.log(e - 1) + (e - 1) ** 2 / 2 + 0.5 - 2 * math.log(e) + math.log(4) - e * e / 4
        imag = math.log(e - 1) - (e - 1) ** 2 / 2 - 0.5 - 2 * math.log(e) + math.log(4) + e * e / 4
        assert abs(GP.dom_gap(e) - real) < 1e-12 and abs(GP.dom_gap(1j * e) - imag) < 1e-12 and real < 0 and imag < 0
    assert abs(GP.dom_gap(1.999) / (-0.25 * 0.001 ** 3) - 1) < 0.01 and abs(GP.dom_gap(1.999j) / (-0.5 * 0.001 ** 2) - 1) < 0.01
    r = GP.dom_radius(0.345)
    assert 1.89 < r < 1.895 and GP.gap_closed(r - 1e-9, 0.345) < 0 <= GP.gap_closed(r + 1e-9, 0.345)
    assert GP.dom_radius(0.0) == 1.9999 and GP.dom_radius(np.pi / 2) == 1.9999
    for r_ in (0.3, 0.6345, 0.8):
        assert abs(GP.c_need(r_) ** 2 - 4 * math.log(1 / r_) / (1 / r_ - 1)) < 1e-12
    assert abs(GP.c_need(GP.r_star(1.8927)) - 1.8927) < 1e-9 and GP.c_need(0.5) < GP.c_need(0.8) < 2
    assert abs(GP.alias_rate(2.0) - (5 / 8 * 4 - math.log(2) - 0.5 - 1.5 * math.log(2))) < 1e-15 and GP.alias_rate(1.858) < 0 < GP.alias_rate(1.86)
    assert GP.D_K == (20, 40, 80, 160) and GP.D_ETA == (1.0, 1.5, 1.8, 1.9) and GP.D_TH == (0, 1, 2, 3, 4)
    assert GP.D_QUAD == ((12, 1.5, 2), (12, 1.9, 1), (40, 1.5, 2), (40, 1.9, 4), (40, 1.95, 1)) and (GP.DOM_S, GP.DOM_TH) == (4001, 361)
    assert GP.DOM_MARGIN == (1.7, 1.78, 1.85) and GP.DOM_NR == 120
    assert GP.QX_Y == (1.2, 1.5, 1.8, 1.9, 2.0, 2.2, 2.5)


def test_round91_part_d_jobs_and_lines(monkeypatch):
    """Part D's jobs at small K with the real functions — `duhamel_job` reproduces `gen_ratio_c` by quadrature, `d_job` stores
    `gen_ratio_c` at η and `corner_ratio` at η/2, `qx_job` the closed form at K = 2 (sin(√2y)e^{y²}/(√2y)) —, `part_D`'s job lists, and
    `d_lines` on planted results."""
    import cmath, math
    import multiprocessing
    q = GP.duhamel_job((6, 1.7, 3))
    assert abs(complex(*q['quad']) / complex(*q['exact']) - 1) < 1e-12 and q['n'] == 38 and q['prec'] == GP.d_prec(6, 1.7, 3)
    assert GP.d_prec(20, 1.9, 4) == 4 * 20 + 80 + int(math.ceil((1.9 * math.sqrt(20)) ** 2 / (2 * math.log(10)))) and GP.d_prec(20, 1.9, 0) == 160
    monkeypatch.setattr(GP, 'D_ETA', (1.0, 1.9))
    d = GP.d_job((6, 4))
    assert d['eta'] == [1.0, 1.9] and d['prec'] == GP.d_prec(6, 1.9, 4) != GP.d_prec(6, 1.0, 4)
    for e, a, b in zip(d['eta'], d['g'], d['c']):
        assert tuple(a) == GP.gen_ratio_c(6, e, 4, d['prec']) and tuple(b) == GP.corner_ratio(6, e / 2, 4, d['prec'])
    monkeypatch.setattr(GP, 'QX_Y', (0.5, 1.5))
    x = GP.qx_job(2)
    assert all(abs(c - math.sin(math.sqrt(2) * y) * math.exp(y * y) / (math.sqrt(2) * y)) < 1e-12 * max(1, abs(c)) for y, c in zip(x['y'], x['c']))
    monkeypatch.setattr(multiprocessing, 'Pool', _SerialPool)
    monkeypatch.setattr(GP, 'd_job', lambda a: ('d', a))
    monkeypatch.setattr(GP, 'duhamel_job', lambda a: ('q', a))
    monkeypatch.setattr(GP, 'qx_job', lambda a: ('x', a))
    monkeypatch.setattr(GP, 'part_dom', lambda: 'dom')
    monkeypatch.setattr(GP, 'qxp_job', lambda a: ('p', a))
    monkeypatch.setattr(GP, 'part_dom', lambda: dict(d='dom'))
    monkeypatch.setattr(GP, 'dom_bb', lambda R, ball=False: (R, -1.0, 0) + ((-2.0,) if ball else ()))
    monkeypatch.setattr(GP, 'part_dpath', lambda: 'path')
    assert GP.part_D() == dict(jobs=[('d', (K, j)) for K in GP.D_K for j in GP.D_TH], quad=[('q', a) for a in GP.D_QUAD],
                               qx=[('x', K) for K in GP.D_K], qxp=[('p', K) for K in GP.QXP_K], dpath='path',
                               dom=dict(d='dom', bb=[dict(R=R, n=R, worst=-1.0, open=0, ball_worst=-2.0) for R in GP.DOM_BB]))
    lim = lambda e, j: cmath.sqrt(1 - (e * cmath.exp(1j * j * math.pi / 8)) ** 2 / 4)
    mk = lambda K, j: dict(K=K, j8=j, prec=1, eta=[1.0, 1.5], g=[[(lim(e, j) * (1 + 1 / K)).real, (lim(e, j) * (1 + 1 / K)).imag] for e in (1.0, 1.5)],
                           c=[[(lim(e, j) ** 2).real, (lim(e, j) ** 2).imag] for e in (1.0, 1.5)])
    D = dict(jobs=[mk(K, j) for K in (20, 40) for j in (0, 4)], quad=[dict(K=12, eta=1.5, j8=2, prec=1, n=56, quad=[1.0 + 1e-12, 0.0], exact=[1.0, 0.0])],
             qx=[dict(K=K, prec=1, y=[1.5, 2.0], c=[3.25 * (1 + 1 / K), math.exp(0.25 * K)]) for K in (20, 40)],
             dom=dict(th=[0.0, 0.5], radius=[1.9999, 1.89], axes=[-2.5e-10, -5e-7], viol={'1.95': [0.1, 0.6], '1.999': None},
                      margin={'1.7': [-0.05, 0.0], '1.78': [-0.01, 0.1]}, bb=[dict(R=1.78, n=7, worst=-1e-3, open=0)]),
             dpath=dict(th=[0.0, 1.0], rows=[dict(R=1.95, h=0.01, ok=[True, False]),
                                             dict(R=1.97, h=0.005, ok=[True, True], ma=[-1e-3, -2e-5], mb=[-3e-4, -1e-2], q=[4e-4, 5e-3])]))
    S = [dict(r=2.0, frame='half-squeeze', gx=1 / (2 * 0.5), gy=1.5)]
    R = dict(rules=[dict(r=1.0, name='Fock', gx=2.0, gy=1 / (2 * 0.6))], none=[])
    s = '\n'.join(GP.d_lines(D, S, R))
    assert "| 12 | 1.5 | 2π/8 | 56 | 1.0e-12 |" in s and "| 20 | 0π/8 | 5.0e-02 (" in s and "| 40 | 4π/8 | 2.5e-02 (" in s
    assert "K times the step's deviation at K = 40 (|η| = 1, 1.5), by arg η: 0π/8: 1.000, 1.000; 4π/8: 1.000, 1.000" in s
    assert "(K⁻¹'s 20/40 = 0.5): 0π/8: 0.500, 0.500; 4π/8: 0.500, 0.500." in s
    assert "for |η| < 1.890000 (the least radius, at arg η = 0.500)" in s and ("violated at |η| = 1.95 for at least arg η ∈ [0.100, 0.600] and at "
                                                    "1.999 for at least none; the gap's largest value on the grid (120 radii from 1 to c times the angles) "
                                                    "-5.00e-02 at |η| ≤ 1.7 (at arg η = 0.000), -1.00e-02 at |η| ≤ 1.78 (at arg η = 0.100).") in s
    assert "verifies the condition on 1 < |η| ≤ 1.78, 0 ≤ arg η ≤ π/2 (7 boxes, the largest bound -1.0e-03, 0 unresolved)." in s
    assert "at |η| = 1.95: 1 of 2 angles (grid step 0.01); |η| = 1.97: 2 of 2 angles (grid step 0.005; the seeds' F − F(η/2) -2.0e-05 … " \
        "-1.0e-02, the chord's ratio q at most 5.0e-03) — 2 angles from 0 to 0.31831π" in s
    assert "on the annulus 1.78 < |η| ≤ 1.97 (the branch and bound covers |η| ≤ 1.78)" in s and "chord between the seeds" in s
    rf = math.floor(GP.r_star(1.97) * 1e4) / 1e4                       # D3cp-5: σ*/σ from the floored r, floored
    assert f"(iv) would cover every frame with r < {rf:.4f} (σ*/σ < {math.floor((1 + rf) / (1 - rf) * 100) / 100:.2f} at the half" in s
    assert (f"R's 1: r ≤ 0.6000, |η| ≤ {GP.c_need(0.6):.4f} (τ/σ = 1, Fock); S's 1: r ≤ 0.5000, |η| ≤ {GP.c_need(0.5):.4f} (τ/σ = 2, half-squeeze); "
            f"every frame with r < {math.floor(GP.r_star(1.78) * 1e4) / 1e4:.4f} lies within the verified radius 1.78 (c_need(0.6380) ≤ "
            f"{GP.c_need_ball(0.638):.8f} in ball arithmetic, inside it: referee round 110) — at the half-squeeze rule "
            f"(r = t) σ*/σ < (1 + r)/(1 − r) = ") in s and (f"below r = {math.floor(GP.r_star(1.0) * 1e4) / 1e4:.4f} (|η| ≤ 1, where the "
                                                       "condition is void; c_need(0.0966) ≤ 0.99982566 < 1) no check") in s
    assert "y = 1.5: +5.0e-02, +2.5e-02" in s and "beyond y = 1.8589, where the first alias's rate" in s
    assert GP.d_lines(None) == [] and GP.d_lines(dict(jobs=[])) == []


def test_round91_stored_fields(monkeypatch):
    """Referee round 91 (D3cl-5: SJ1, LJ2, LJ3, raw-only survivors): `lc_job` at K = 20, θ = π/8 and 3π/8 reproduces the raw's LC
    entries (its precision from Im ξ, its signed args); `s_job`'s check is the digits the two series agree to, their relative difference
    taken at the higher precision."""
    import math
    from decimal import Decimal, localcontext
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    for j in (1, 3):
        out = GP.lc_job((20, j))
        x = next(x for x in raw['LC'] if (x['K'], x['j8']) == (20, j))
        assert {k: (list(v) if isinstance(v, (list, tuple, np.ndarray)) else v) for k, v in out.items()} == \
               {k: (list(v) if isinstance(v, (list, tuple, np.ndarray)) else v) for k, v in x.items()}
    assert any(a < 0 for x in raw['LC'] for a in x['arg'])
    monkeypatch.setattr(GP, 'S_K', (4,))
    rho, vp = 3 - 2 * math.sqrt(2), 2.0
    out = GP.s_job((1.0, 'half-squeeze', rho, vp))
    p = out['rows'][-1]['prec']
    S1 = WM.frame_series(4, Decimal(repr(float(rho))), Decimal(repr(float(vp))), prec=p)[0]
    S2 = WM.frame_series(4, Decimal(repr(float(rho))), Decimal(repr(float(vp))), prec=p + GP.S_CHECK)[0]
    with localcontext() as ctx:
        ctx.prec = p + GP.S_CHECK
        assert out['check'] == float(-abs(S2 / S1 - 1).log10())
    with localcontext() as ctx:
        ctx.prec = p
        assert out['check'] != float(-abs(S2 / S1 - 1).log10())


def test_round91_part_dom_and_main_d(monkeypatch, tmp_path):
    """`part_dom` on a small grid against its definition — the radius per angle, the axes' gaps, the violated sectors' first and last
    angles, the margin as the largest gap over DOM_NR radii from 1 to c times the angles —, and `main --only D` saving part_D's output."""
    import exp_cohesion_leak as E
    monkeypatch.setattr(GP, 'DOM_S', 401)
    monkeypatch.setattr(GP, 'DOM_TH', 7)
    monkeypatch.setattr(GP, 'DOM_NR', 9)
    dm = GP.part_dom()
    ths = np.linspace(0, np.pi / 2, 7)
    assert dm['th'] == [float(t) for t in ths] and dm['radius'] == [GP.dom_radius(t) for t in ths]
    assert dm['axes'] == [GP.gap_closed(1.999, 0.0), GP.gap_closed(1.999, np.pi / 2)]
    for r_ in (1.95, 1.999):
        v = [float(t) for t in ths if GP.gap_closed(r_, t) >= 0]
        assert len(v) >= 2 and dm['viol'][str(r_)] == [v[0], v[-1]] and v[0] < v[-1]
    for c in GP.DOM_MARGIN:
        rs = np.linspace(1 + 1e-9, c, 9)
        assert tuple(dm['margin'][str(c)]) == max((GP.gap_closed(r_, t), float(t)) for r_ in rs for t in ths) and dm['margin'][str(c)][0] < 0
    monkeypatch.setattr(GP, 'OUT', str(tmp_path))
    monkeypatch.setattr(GP, 'part_D', lambda: dict(jobs=[], quad=[], qx=[], dom={}))
    monkeypatch.setattr(GP, 'part_Q', lambda: 1 / 0)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'D'])
    GP.main()
    raw = E.load(str(tmp_path / 'gh_peano_raw.npz'))
    assert set(raw) == {'D'} and dict(raw['D']) == dict(jobs=[], quad=[], qx=[], dom={})


def test_round93_closed_form_bb_prony_and_the_path():
    """Referee round 93 (D3cn-1, D3cn-3): `gap_closed` is the exact supremum of `dom_gap`'s grid (within the grid's resolution, never
    below it); `dom_bb` verifies 1 < |η| ≤ 1.78 and finds an unresolved box beyond the least radius; `prony` recovers a planted rate and
    frequency; `path_ok` holds at |η| = 1.95 in the sector where the straight segment fails; the constants pinned."""
    import math
    for R, f in ((1.5, 0.0), (1.9, 0.345), (1.95, 0.35), (1.999, math.pi / 2), (1.7, 1.0)):
        g, c = GP.dom_gap(R * np.exp(1j * f)), GP.gap_closed(R, f)
        assert c >= g - 1e-15 and c - g < 1e-6
    assert GP.gap_closed(1.0, 0.3) == -math.inf and GP.gap_closed(1.95, 0.35) > 0
    n, worst, op = GP.dom_bb(1.78)
    assert op == 0 and n > 100 and worst < -GP.DOM_SLACK
    n2, _, op2 = GP.dom_bb(1.95, max_boxes=20000)
    assert op2 > 0
    k = np.arange(13)
    c = 3.0 * np.exp(0.25 * k) * np.cos(2.1 * k + 0.4)
    lam, fr = GP.prony(c, 0.2)
    assert abs(lam - 0.25) < 1e-9 and abs(fr - 2.1) < 1e-9
    assert GP.path_ok(1.95 * np.exp(0.35j), 0.002) and GP.dom_gap(1.95 * np.exp(0.35j)) > 0   # round 95: h/√2 < 0.2rb
    assert (GP.QXP_K, GP.QXP_N, GP.DOM_BB, GP.DOM_SLACK) == ((100, 150), 13, (1.78, 1.8927), 1e-12)
    assert (GP.DP_R, GP.DP_TH, GP.DP_H) == ((1.9, 1.95, 1.98), 19, (0.004, 0.002, 0.001))
    import pytest
    mp = pytest.MonkeyPatch()
    mp.setattr(GP, 'QX_Y', (0.5,))
    assert GP.qx_job(3)['prec'] == 4 * 3 + 80                          # D3cn-8: qx's stored precision pinned
    mp.undo()
    assert GP.DOM_TH == 361 and GP.DOM_NR == 120 and GP.DOM_MARGIN == (1.7, 1.78, 1.85)


def test_round93_qxp_job_and_the_closed_radius(monkeypatch):
    """Round 93's mutants G04, G17: `qxp_job` stores c_K(i√Ky) − (1 + y²) at QXP_N consecutive K (the y beyond `alias_rate`'s zero only);
    `dom_radius` is the crossing of `gap_closed`, to 10⁻¹², not of the grid's `dom_gap`."""
    import math
    monkeypatch.setattr(GP, 'QX_Y', (1.5, 2.0))
    monkeypatch.setattr(GP, 'QXP_N', 3)
    x = GP.qxp_job(4)
    assert x['y'] == [2.0] and [r['K'] for r in x['rows']] == [4, 5, 6]
    for r in x['rows']:
        assert r['c'][0] == GP.corner_ratio(r['K'], 2.0, 4, r['prec'])[0] - 5.0 and r['prec'] == 4 * r['K'] + 80
    th = 0.3432
    a, b = 1.8, 1.95
    for _ in range(80):
        m = (a + b) / 2
        a, b = (m, b) if GP.gap_closed(m, th) < 0 else (a, m)
    assert abs(GP.dom_radius(th) - a) < 1e-12 and GP.gap_closed(1.8927, th) < 0 < GP.gap_closed(1.89275, th)


def test_round93_the_box_bound_is_an_upper_bound():
    """`box_bound` is not below `gap_closed` anywhere in its box (G09, G18: a bound taken at the wrong corner would verify boxes it should
    not), on 300 random boxes inside 1 < |η| < 2, 0 ≤ arg η ≤ π/2, each sampled at 25 points."""
    rng = np.random.default_rng(3)
    for _ in range(300):
        R1 = 1 + rng.random() * 0.95
        R2 = min(R1 + rng.random() * 0.05, 1.999)
        f1 = rng.random() * (np.pi / 2 - 0.05)
        f2 = f1 + rng.random() * 0.05
        b = GP.box_bound(R1, R2, f1, f2)
        for R in np.linspace(R1, R2, 5):
            for f in np.linspace(f1, f2, 5):
                assert GP.gap_closed(R, f) <= b + 1e-12, (R1, R2, f1, f2, R, f)


def test_round93_box_bound_corners():
    """`box_bound` takes the least R²cos 2φ over the box's corners and the largest R cos φ (G18)."""
    import math
    for R1, R2, f1, f2 in ((1.5, 1.6, 0.2, 0.9), (1.85, 1.9, 0.0, 1.5), (1.2, 1.9, 0.7, 1.2)):
        p_lo = min(R * R * math.cos(2 * f) for R in (R1, R2) for f in (f1, f2))
        assert GP.box_bound(R1, R2, f1, f2) == GP.bound_sup(R1, R2, p_lo, R2 * math.cos(f1))
    assert GP.box_bound(0.5, 1.0, 0.0, 0.1) == -math.inf


def _h_sup(R1, R2, p, B, n=200001):
    """bound_sup's function on a dense grid of w ∈ (0, 1 − 1/R2]."""
    import math
    w = np.linspace(1e-6, 1 - 1 / R2, n)
    a = np.abs(B * (1 - w))
    return float(np.max(np.log(w) + p / 2 * (w * w - 0.5) + np.where(a <= 1, a * a, 2 * a - 1) - 0.5 - math.log(R1) + math.log(4)))


def test_round95_stable_roots_and_bound_sup_against_a_dense_grid():
    """Referee round 95 (D3cp-2, D3cp-3: H02, H06): `quad_roots` keeps both roots when qa is tiny (R²cos(π/2) ≈ 6·10⁻¹⁷) and drops to the
    linear root only below a relative 10⁻¹⁴; `bound_sup` is the supremum of its function — log R1, not log R2 — at random parameters and
    at p = 0 and p ≈ 6·10⁻¹⁷; the box [1.5, 1.8927] × [0, π/4], whose bound the round-93 code put at 0.8458, bounds by 0.9351."""
    import math
    for qa, qb, qc in ((1.0, -3.0, 2.0), (2.0, 1.0, -1.0), (-1.5, 2.0, 1.0), (3e-13, -2.0, 1.0)):
        r = GP.quad_roots(qa, qb, qc)
        assert len(r) == 2 and all(abs(qa * w * w + qb * w + qc) < 1e-9 * max(1.0, w * w) for w in r)
    assert GP.quad_roots(1.0, 0.0, 1.0) == [] and GP.quad_roots(0.0, 0.0, 1.0) == []
    assert GP.quad_roots(6e-17, -2.0, 1.0) == [0.5]
    tiny = 1.5 ** 2 * math.cos(math.pi / 2)
    assert 0 < tiny < 1e-15
    rng = np.random.default_rng(95)
    cases = [(1.5, 1.8927, tiny, 1.8927), (1.5, 1.8927, 0.0, 1.8927), (1.2, 1.9, -2.0, 0.4), (1.9, 1.9, 1.0, 1.9)]
    cases += [(R1, R1 + rng.random() * 0.2, rng.uniform(-4, 4), rng.uniform(0.0, 2.0)) for R1 in 1.05 + 0.9 * rng.random(12)]
    for R1, R2, p, B in cases:
        b, d = GP.bound_sup(R1, R2, p, B), _h_sup(R1, R2, p, B)
        assert d - 1e-9 <= b <= d + 1e-6, (R1, R2, p, B, b, d)
    assert abs(GP.bound_sup(1.5, 1.9, 0.5, 1.2) - GP.bound_sup(1.9, 1.9, 0.5, 1.2) - math.log(1.9 / 1.5)) < 1e-12
    b = GP.box_bound(1.5, 1.8927, 0.0, math.pi / 4)
    assert abs(b - 0.93508) < 1e-4 and b >= _h_sup(1.5, 1.8927, tiny, 1.8927) - 1e-9
    for R, f in zip(1.0 + 0.999 * rng.random(20), rng.random(20) * math.pi / 2):     # H02: the closed form against s-grids, random η
        g, c = GP.dom_gap(R * np.exp(1j * f)), GP.gap_closed(R, f)
        assert c >= g - 1e-15 and c - g < 1e-5, (R, f, c, g)


def test_round95_dom_bb_covers_the_quadrant_and_reports_the_worst(monkeypatch):
    """H10, H12, H13: the verified boxes tile 1 < |η| ≤ Rmax × [0, π/2] exactly; `worst` is the largest verified bound; a bound of
    −5·10⁻¹⁰ is verified with the slack 10⁻¹² (DOM_SLACK), −5·10⁻¹³ is not; a planted violation at arg η > π/4 is not verified."""
    import math
    seen = []
    real = GP.box_bound

    def rec(R1, R2, f1, f2):
        b = real(R1, R2, f1, f2)
        seen.append((R1, R2, f1, f2, b))
        return b
    monkeypatch.setattr(GP, 'box_bound', rec)
    n, worst, op = GP.dom_bb(1.6)
    ver = [x for x in seen if x[4] < -GP.DOM_SLACK]
    assert op == 0 and n == len(ver) and worst == max(x[4] for x in ver) and worst < min(x[4] for x in seen if x[4] >= -GP.DOM_SLACK)
    assert abs(sum((x[1] - x[0]) * (x[3] - x[2]) for x in ver) - 0.6 * math.pi / 2) < 1e-12
    assert min(x[2] for x in ver) == 0.0 and max(x[3] for x in ver) == math.pi / 2 and max(x[1] for x in ver) == 1.6
    monkeypatch.setattr(GP, 'box_bound', lambda *a: -5e-10)
    assert GP.dom_bb(1.6) == (1, -5e-10, 0)
    monkeypatch.setattr(GP, 'box_bound', lambda *a: -5e-13)
    n, _, op = GP.dom_bb(1.6)
    assert n == 0 and op > 0
    monkeypatch.setattr(GP, 'box_bound', lambda R1, R2, f1, f2: 1.0 if f2 > 1.2 and R2 > 1.5 else real(R1, R2, f1, f2))
    n, _, op = GP.dom_bb(1.6)
    assert op > 0


def test_round95_the_seeds_the_bound_f_and_negative_controls(monkeypatch):
    """D3cp-1, D3cp-3 (H15, H18, H21, H22, H23): the seeds lie at PATH_SEED·rb = 1.2rb from η/2 on the steepest-descent line (u²φ″ < 0
    real), the first toward 0; `path_F` carries Λ⁺ from |σ| = 1 on; at 1.95e^{0.35i} the straight segment rises above F(η/2) outside the
    ball (at |σ| > 1, so only with Λ⁺) while the deformed check passes; with Λ⁺ + 50 it fails, and with the sublevel set enlarged by 10
    it fails (one component); the check refuses a grid too coarse for the ball; at |η| = 1.9 on the axis a tolerance of 10⁻⁴ would fail
    the check (margin −9.8·10⁻⁵) and the strict one passes; a bump on the chord between the seeds fails it; part_dpath's angles span
    [0, π/2]."""
    import math
    import pytest
    assert (GP.PATH_SEED, GP.PATH_TOL) == (1.2, 0.0)
    for eta in (1.9 + 0j, 1.95 * np.exp(0.35j), 1.98 * np.exp(1.2j), 1.5j, 1.9 * np.exp(-0.8j), -1.7 - 0.6j):   # the last two need the flip
        a, b, rb = GP.path_seeds(eta)
        assert abs(rb - 0.5 * (1 - abs(eta) / 2)) < 1e-15
        assert abs(abs(a - eta / 2) - 1.2 * rb) < 1e-12 and abs(abs(b - eta / 2) - 1.2 * rb) < 1e-12 and abs(a + b - eta) < 1e-12
        u = (b - eta / 2) / (1.2 * rb)
        q = (2 - 8 / eta ** 2) * u * u
        assert q.real < 0 and abs(q.imag) < 1e-12 * abs(q) and abs(a) < abs(b)
    e = 1.95 * np.exp(0.35j)
    for z in (0.5 + 0.2j, 1.02 * np.exp(0.2j), 1.3 - 0.4j):
        ref = math.log(abs(e - z)) + ((e - z) ** 2).real / 2 + math.log(abs(z)) + (z * z).real / 2 + (0.0 if abs(z) < 1 else float(GP.lam_plus(z)))
        assert abs(float(GP.path_F(e, z)) - ref) < 1e-12
    F0 = 2 * math.log(abs(e) / 2) + (e * e).real / 4
    assert abs(0.6 * e) > 1 and float(GP.path_F(e, 0.6 * e)) > F0 + 1e-3 and abs(0.6 * e - e / 2) > GP.path_seeds(e)[2]
    ok, ma, mb, q = GP.path_check(e, 0.002)
    assert ok and ma < 0 and mb < 0 and 0 < q < 1e-3 and q == GP.chord_q(e, abs(GP.path_seeds(e)[1] - e / 2))
    with pytest.raises(ValueError):
        GP.path_check(e, 0.004)
    lp = GP.lam_plus
    monkeypatch.setattr(GP, 'lam_plus', lambda mu: lp(mu) + 50)
    assert not GP.path_check(e, 0.002)[0]
    monkeypatch.setattr(GP, 'lam_plus', lp)
    cq = GP.chord_q
    monkeypatch.setattr(GP, 'chord_q', lambda eta, D: 1.0)                 # round 99: a chord ratio at 1 alone fails the check
    r = GP.path_check(e, 0.002)
    assert r[0] is False and r[3] == 1.0 and r[1] < 0 and r[2] < 0
    monkeypatch.setattr(GP, 'chord_q', cq)
    monkeypatch.setattr(GP, 'PATH_TOL', 10.0)
    assert not GP.path_check(e, 0.002)[0]
    monkeypatch.setattr(GP, 'PATH_TOL', 1e-4)
    assert not GP.path_check(1.9 + 0j, 0.004)[0]
    monkeypatch.setattr(GP, 'PATH_TOL', 0.0)
    ok, ma, mb, q = GP.path_check(1.9 + 0j, 0.004)
    assert ok and -1e-4 < ma < 0 and abs(q - 5.814e-3) < 1e-6
    monkeypatch.setattr(GP, 'path_check', lambda eta, h: (True, -1.0, -2.0, -3.0))
    monkeypatch.setattr(GP, 'DP_R', (1.9,))
    monkeypatch.setattr(GP, 'DP_H', (0.004,))
    x = GP.part_dpath()
    assert len(x['th']) == GP.DP_TH and x['th'][0] == 0.0 and x['th'][-1] == np.pi / 2 and x['seeds'] == 'descent'
    assert x['rows'] == [dict(R=1.9, h=0.004, ok=[True] * GP.DP_TH, ma=[-1.0] * GP.DP_TH, mb=[-2.0] * GP.DP_TH, q=[-3.0] * GP.DP_TH)]


def test_round95_qxp_uses_the_y_beyond_the_alias_zero(monkeypatch):
    """H26: on the published QX_Y, `qxp_job` fits the y where `alias_rate` > 0 — 1.9, 2.0, 2.2, 2.5, not 1.8 (rate −0.103)."""
    assert GP.QX_Y == (1.2, 1.5, 1.8, 1.9, 2.0, 2.2, 2.5)
    monkeypatch.setattr(GP, 'QXP_N', 1)
    x = GP.qxp_job(3)
    assert x['y'] == [1.9, 2.0, 2.2, 2.5] and -0.2 < GP.alias_rate(1.8) < 0 < GP.alias_rate(1.9)


def test_round99_the_taylor_chord_the_seed_pixel_and_the_sector(monkeypatch):
    """Referee round 99 (D3ct-1, D3ct-4, D3ct-6): `chord_q` bounds F − F(η/2) on the chord by −(|φ″|t²/2)(1 − q(t/D)²) — checked against
    a dense sampling of the chord at 30 points of the annulus, the sampled F − F(η/2) never above the bound —; q grows to ≥ 1 for a chord
    reaching the ball's edge at 2rb only far from the saddle's scale; the seed pixel's guard raises when the rounded pixel falls in the ball
    (K25); `path_ok` returns a bool (K31); `part_dpath` pairs radii and steps in order (K32); `dom_sector` brackets the violated sector by
    bisection and its maximisers' Re(ηs) ≥ 1."""
    import math
    rng = np.random.default_rng(99)
    for R, f in zip(rng.uniform(1.8927, 1.98, 30), rng.uniform(0, math.pi / 2, 30)):
        eta = R * np.exp(1j * f)
        a, b, rb = GP.path_seeds(eta)
        D = abs(b - eta / 2)
        q = GP.chord_q(eta, D)
        assert 0 < q < 0.01
        u = (b - eta / 2) / D
        F0 = 2 * math.log(R / 2) + (eta * eta).real / 4
        ph2 = abs(2 - 8 / eta ** 2)
        for t in np.linspace(-D, D, 201):
            if t == 0:
                continue
            Ft = float(GP.path_F(eta, eta / 2 + t * u)) - F0
            assert Ft <= -(ph2 * t * t / 2) * (1 - q * (t / D) ** 2) + 1e-15
    eta = 1.95 * np.exp(0.35j)
    assert GP.chord_q(eta, 0.012) < GP.chord_q(eta, 0.024) < GP.chord_q(eta, 0.4)
    assert GP.path_ok(eta, 0.002) is True
    monkeypatch.setattr(GP, 'path_seeds', lambda e: (e / 2 - 0.0099, e / 2 + 0.0099, 0.0095))   # seeds a hair outside a 0.0095 ball
    import pytest
    with pytest.raises(ValueError):
        GP.path_check(1.95 + 0j, 0.002)                                  # h/√2 = 0.0014 < 0.0019 = 0.2rb, but a pixel rounds inside (0.009)
    monkeypatch.undo()
    calls = []
    monkeypatch.setattr(GP, 'path_check', lambda e, h: (calls.append((round(abs(e), 6), h)) or (True, -1.0, -1.0, 1e-3)))
    monkeypatch.setattr(GP, 'DP_R', (1.9, 1.98))
    monkeypatch.setattr(GP, 'DP_H', (0.004, 0.001))
    GP.part_dpath()
    assert sorted(set(calls)) == [(1.9, 0.004), (1.98, 0.001)]
    monkeypatch.undo()
    sc = GP.dom_sector()
    assert sc['n'] == GP.DOM_SECT == 20001 and abs(sc['lo'] - 0.0860576) < 1e-6 and abs(sc['hi'] - 0.6175580) < 1e-6
    assert GP.gap_closed(1.95, sc['lo'] - 1e-9) < 0 < GP.gap_closed(1.95, sc['lo'] + 1e-9) and 1.005 < sc['re_min'] < 1.0052 < sc['re_max'] < 1.1
    assert GP.dom_sector(1.85) is None


def test_round99_planted_edges(monkeypatch):
    """Round 99's own mutants T08, T10, T11, T13: the pixel guard covers both seeds (only b's pixel in the ball here); the sector's ends
    to 10⁻¹² (80 bisections); the maximiser's Re(ηs) taken at the sector's ends too (the least there, 1.005108, below the grid's 1.005156);
    `part_dom` stores the sector (`dom_sector`), on a coarse grid."""
    import pytest
    monkeypatch.setattr(GP, 'path_seeds', lambda e: (e / 2 - 0.0111, e / 2 + 0.0099, 0.0095))
    with pytest.raises(ValueError):
        GP.path_check(1.95 + 0j, 0.002)
    monkeypatch.undo()
    sc = GP.dom_sector()
    for end in (sc['lo'], sc['hi']):
        assert GP.gap_closed(1.95, end - 1e-11) * GP.gap_closed(1.95, end + 1e-11) < 0
    assert sc['re_min'] < 1.00512
    monkeypatch.setattr(GP, 'DOM_TH', 5)
    monkeypatch.setattr(GP, 'DOM_NR', 3)
    monkeypatch.setattr(GP, 'DOM_SECT', 201)
    d = GP.part_dom()
    assert d['sect'] == GP.dom_sector() and d['sect']['n'] == 201


def test_part_v_the_valley_certificate(monkeypatch):
    """Part V (referee round 99's recipe): `cert_E` is F_η(σ) − F_η(η/2) at σ = (η/2)(1 + z), inside and outside the unit disc; on the
    straight segment inside the disc it is negative but at the saddle (the analytic part); `cert_seg_ub` bounds `cert_E` from above over
    its segment and η-disc (sampled at random segments and boxes); `cert_end` bounds the last piece; a box verifies, and fails with the
    bound raised by 10⁻³ (a negative control); `cert_band` and `part_V` wire the bands; `v_lines` reports the certified radius."""
    import cmath
    import math
    rng = np.random.default_rng(61)
    for _ in range(400):
        R, f = rng.uniform(1.5, 1.999), rng.uniform(0, math.pi / 2)
        eta = R * cmath.exp(1j * f)
        z = complex(rng.uniform(-1.5, 1.5), rng.uniform(-1, 1))
        if abs(1 - z * z) < 1e-3:
            continue
        F0 = 2 * math.log(R / 2) + (eta * eta).real / 4
        assert abs(GP.cert_E(z, eta) - (float(GP.path_F(eta, eta / 2 * (1 + z))) - F0)) < 1e-12
        tmax = (2 - R) / R
        for x in np.linspace(-0.999, tmax, 50):
            if abs(x) > 1e-9:
                assert GP.cert_E(complex(x), eta) < 0                    # the straight segment inside the disc
    for _ in range(300):
        R2 = rng.uniform(1.89, 1.99)
        R1, f1 = R2 - rng.uniform(0, 0.01), rng.uniform(0, 1.5)
        f2 = f1 + rng.uniform(0, 0.02)
        ec = (R1 + R2) / 2 * cmath.exp(1j * (f1 + f2) / 2)
        d = max(abs(R * cmath.exp(1j * f) - ec) for R in (R1, R2) for f in (f1, f2))
        p = complex(rng.uniform(-0.2, 0.9), rng.uniform(-0.5, 0.5))
        q = p + complex(rng.uniform(-0.05, 0.05), rng.uniform(-0.05, 0.05))
        ub = GP.cert_seg_ub(p, q, ec, d, R2)
        for t in np.linspace(0, 1, 7):
            for R in (R1, R2):
                for f in (f1, f2):
                    assert GP.cert_E(p + t * (q - p), R * cmath.exp(1j * f)) <= ub + 1e-12
    for R2 in (1.9, 1.99):
        assert GP.cert_end(R2) < -1
        eta = R2 * cmath.exp(0.4j)
        for w in (1 - 0.02, 1 - 0.02j, 1 + 0.015 + 0.01j):
            for t in np.linspace(0, 0.999, 30):
                assert GP.cert_E(w + t * (1 - w), eta) <= GP.cert_end(R2)
    ok, n, w = GP.cert_box(1.95, 1.955, 0.3, 0.31)[:3]
    assert ok and n > 10 and w < -GP.CERT_SLACK
    ce = GP.cert_E
    monkeypatch.setattr(GP, 'cert_E', lambda z, eta: ce(z, eta) + 1e-3)
    assert not GP.cert_box(1.95, 1.955, 0.3, 0.31)[0]
    monkeypatch.setattr(GP, 'cert_E', ce)
    b = GP.cert_band(1.9, 1.902, 1, 2, ball=True)
    assert b['boxes'] >= 2 and b['fails'] == 0 and b['end'] == GP.cert_end(1.902) and b['worst'] < -GP.CERT_SLACK
    monkeypatch.setattr(GP, 'CERT_BANDS', ((1.9, 1.902, 1, 2),))
    v = GP.part_V()
    assert v['bands'] == [b] and v['tb'] == GP.CERT_TB == 0.9 and v['eps0'] == GP.CERT_EPS0 == 0.02 and v['slack'] == 1e-12
    assert v['ball']['prec'] == GP.CERT_BALL_PREC == 120 and b['ball_worst'] < 0 and b['ball_end'] < 0
    D = dict(dom=dict(bb=[dict(R=1.9, n=1, worst=-1.0, open=0)]))
    s = '\n'.join(GP.v_lines(v, D))
    rv = math.floor(GP.r_star(1.902) * 1e4) / 1e4
    assert "the limit of (iii)(b) holds on |η| ≤ 1.902" in s and f"(iv) for every frame with r < {rv:.4f}" in s
    s2 = '\n'.join(GP.v_lines(dict(v, bands=[dict(b, fails=1)]), D))
    assert "holds on |η| ≤ 1.9 (" in s2                                  # a failed band certifies nothing
    assert GP.v_lines(None) == [] and GP.KEYS[-5:] == ('V', 'DB', 'VC', 'VA', 'VA011')


def test_part_v_the_bound_s_parts_and_the_wiring(monkeypatch):
    """Part V's own mutants: each Lipschitz constant of `cert_seg_parts` bounds the sampled gradients — L_z over the segment at η_c
    (`cert_grad`, itself checked against differences of `cert_E`), L_η over the segment and the η-disc |η − η_c| ≤ δ, |η| ≤ R2 (by
    differences in η) —, on random segments near z = ±1, outside the disc and near the saddle; the box's disc reaches its corners; b
    lies inside the disc for the box's largest |η|; a bound of −5·10⁻¹³ is not below the slack; the end piece's formula; a path that does
    not reach z = 1 fails; the band's boxes tile [Ra, Rb] × [0, π/2]; failures are counted."""
    import cmath
    import math
    rng = np.random.default_rng(611)
    dE = lambda z, eta, h=1e-7: complex((GP.cert_E(z + h, eta) - GP.cert_E(z - h, eta)) / (2 * h),
                                       (GP.cert_E(z + 1j * h, eta) - GP.cert_E(z - 1j * h, eta)) / (2 * h))
    for z in (0.3 + 0.2j, 0.8 - 0.3j, -0.5 + 0.1j, 1.2 + 0.4j):
        eta = 1.95 * cmath.exp(0.3j)
        assert abs(GP.cert_grad(z, eta) - dE(z, eta)) < 1e-6
    for _ in range(250):
        R2 = rng.uniform(1.89, 1.995)
        R1, f1 = R2 - rng.uniform(0, 0.03), rng.uniform(0, 1.5)
        f2 = f1 + rng.uniform(0, 0.06)
        ec, d = GP.cert_box_disc(R1, R2, f1, f2)
        kind = rng.integers(3)
        p = (complex(rng.uniform(0.6, 0.95), rng.uniform(-0.3, 0.3)) if kind == 0 else
             complex(rng.uniform(-0.3, 1.2), rng.uniform(0.2, 0.8)) if kind == 1 else complex(rng.uniform(-0.05, 0.05), rng.uniform(-0.05, 0.05)))
        q = p + complex(rng.uniform(-0.15, 0.15), rng.uniform(-0.15, 0.15))
        Em, l, Lz, Leta = GP.cert_seg_parts(p, q, ec, R2)
        if not math.isfinite(Lz):
            continue
        assert Em == GP.cert_E((p + q) / 2, ec) and abs(l - abs(q - p) / 2) < 1e-15
        for t in np.linspace(0, 1, 9):
            z = p + t * (q - p)
            assert abs(GP.cert_grad(z, ec)) <= Lz + 1e-9
            for e in [ec + d * cmath.exp(1j * a) for a in np.linspace(0, 2 * math.pi, 8, endpoint=False)] + [R2 * cmath.exp(1j * f2)]:
                if abs(e) > R2 + 1e-12:
                    e = e / abs(e) * R2
                g = complex((GP.cert_E(z, e + 1e-7) - GP.cert_E(z, e - 1e-7)) / 2e-7, (GP.cert_E(z, e + 1e-7j) - GP.cert_E(z, e - 1e-7j)) / 2e-7)
                assert abs(g) <= Leta + 1e-6
        for R in (R1, R2):
            for f in (f1, f2):
                assert abs(R * cmath.exp(1j * f) - ec) <= d + 1e-15
    for R2 in (1.9, 1.98, 1.995):
        assert R2 / 2 * (1 + GP.cert_b(R2)) < 1 and GP.cert_b(R2) > 0.8 * (2 - R2) / R2
    hh = lambda r: (r * r - 1) / 2 - math.log(r)
    assert GP.cert_end(1.99) == math.log(0.02) + math.log(2.02) + 1.99 ** 2 / 4 * 1.02 ** 2 + hh(1.99 / 2 * 2.02)
    monkeypatch.setattr(GP, 'cert_seg_ub', lambda p, q, e, d, R2: -5e-13)
    assert GP.cert_poly([0.1, 0.2], 1.95, 0.0, 1.95)[0] is False
    monkeypatch.undo()
    cf = GP.cert_flow
    full = cf(1.9525 * cmath.exp(0.305j), GP.cert_b(1.955))
    monkeypatch.setattr(GP, 'cert_flow', lambda eta, z0: full[:len(full) * 2 // 3])
    assert abs(1 - full[len(full) * 2 // 3 - 1]) > GP.CERT_EPS0 and GP.cert_box(1.95, 1.955, 0.3, 0.31)[0] is False   # a valid path cut short
    monkeypatch.undo()
    for eta in (1.95 * cmath.exp(0.3j), 1.99 + 0j):                       # near the singularity z = 1
        Em, l, Lz, Leta = GP.cert_seg_parts(0.9 + 0j, 0.99 + 0j, eta, abs(eta))
        assert all(abs(GP.cert_grad(complex(x), eta)) <= Lz for x in np.linspace(0.9, 0.99, 50))
    monkeypatch.undo()
    seen = []
    monkeypatch.setattr(GP, 'cert_box', lambda R1, R2, f1, f2, ball=False: (seen.append((R1, R2, f1, f2)) or (True, 1, -1.0, -1.0, 0.0, False)))
    b = GP.cert_band(1.9, 1.95, 2, 3)
    assert b['boxes'] == 6 and abs(sum((x[1] - x[0]) * (x[3] - x[2]) for x in seen) - 0.05 * math.pi / 2) < 1e-15
    assert min(x[2] for x in seen) == 0.0 and max(x[3] for x in seen) == math.pi / 2
    monkeypatch.setattr(GP, 'cert_box', lambda R1, R2, f1, f2, ball=False: (False, 0, math.inf, math.inf, math.inf, False))
    monkeypatch.setattr(GP, 'CERT_DEPTH', 1)
    assert GP.cert_band(1.9, 1.95, 2, 3)['fails'] == 12 and GP.CERT_BANDS[-1][:2] == (1.99, 1.995)


def test_round104_part_v_recomputed_pinned_and_in_ball_arithmetic(monkeypatch):
    """Referee round 104 (D3cy-2 … D3cy-5): `cert_box` and a thin slice of `cert_band` recomputed and pinned to their values (so an
    unsound change that only makes the certificate cheaper changes them); |∇_σC| ≤ h′(|σ|) where it binds (|Re σ| < 1 < |σ|, and
    beyond |Re σ| = 1); `cert_seg_ball` an upper bound of the float bound within 10⁻¹³ and of `cert_E` on the segment; the rigorous
    disc radius at the true π/2; `cert_box` starts at cert_b(R2) and rejects a path ending beyond CERT_EPS0; forced splits tile the band
    exactly; the bands contiguous with the branch and bound's radius; a gap rejected by `v_lines`; r* floored; σ*/σ in exact arithmetic."""
    import math
    import cmath
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31)[:3] == (True, 62, -2.7704704904874086e-05)
    assert GP.cert_box(1.9944, 1.9945, 0.0058, 0.0059)[:2] == (False, 72)
    assert GP.cert_box(1.99, 1.995, math.pi / 2 * 31 / 32, math.pi / 2)[:2] == (False, 12)
    b = GP.cert_band(1.99, 1.9902, 1, 4)
    assert (b['boxes'], b['segs'], b['worst'], b['fails']) == (317, 105839, -2.1874538491923193e-12, 0)
    rng = np.random.default_rng(104)
    Cf = lambda sg: 0.0 if abs(sg) <= 1 else (abs(sg) ** 2 - 1) / 2 - math.log(abs(sg)) - max(abs(sg.real) - 1, 0.0) ** 2
    hp = lambda r: r - 1 / r if r > 1 else 0.0
    for _ in range(400):                                                 # |∇C| ≤ h′(|σ|), equality where |Re σ| < 1 < |σ|
        sg = complex(*rng.uniform(-2.2, 2.2, 2))
        if abs(abs(sg) - 1) < 1e-3 or abs(abs(sg.real) - 1) < 1e-3:
            continue
        g = abs(complex(Cf(sg + 1e-7) - Cf(sg - 1e-7), Cf(sg + 1e-7j) - Cf(sg - 1e-7j)) / 2e-7)
        assert g <= hp(abs(sg)) + 1e-6
        if abs(sg.real) < 1 < abs(sg):
            assert abs(g - hp(abs(sg))) < 1e-6
    arb, acb = GP._ball()
    for _ in range(60):                                                 # the ball bound bounds the float one and E on the segment
        R1 = rng.uniform(1.8927, 1.99)
        R2, f1 = R1 + 0.002, rng.uniform(0, 1.5)
        f2 = f1 + 0.01
        ec, d = GP.cert_box_disc(R1, R2, f1, f2)
        p = complex(*rng.uniform(-0.3, 0.9, 2)) * 0.8
        q = p + 0.01 * cmath.exp(1j * rng.uniform(0, 2 * math.pi))
        u, ub = GP.cert_seg_ball(p, q, acb(ec.real, ec.imag), arb(d), arb(R2)), GP.cert_seg_ub(p, q, ec, d, R2)
        if math.isfinite(ub):
            assert ub - 1e-13 <= u <= ub + 1e-13 and all(GP.cert_E(p + t * (q - p), ec) <= u for t in np.linspace(0, 1, 9))
    ec_ = 1.9925 * cmath.exp(1j * (1.5 + math.pi / 2) / 2)
    db = GP.cert_box_disc_ball(1.99, 1.995, 1.5, math.pi / 2, ec_)
    assert float(db.upper()) >= GP.cert_box_disc(1.99, 1.995, 1.5, math.pi / 2)[1]
    top = abs(acb(arb(0), arb(1.995)) - acb(ec_.real, ec_.imag))           # the corner at the true π/2, rigorously
    assert top.lower() <= db and top.lower() > abs(acb(1.995 * math.cos(math.pi / 2), 1.995 * math.sin(math.pi / 2)) - acb(ec_.real, ec_.imag)).upper()
    ref = GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)
    cd = GP.cert_box_disc
    monkeypatch.setattr(GP, 'cert_box_disc', lambda *a: (cd(*a)[0], cd(*a)[1] / 2))   # a float radius short: the rigorous one used
    got_ = GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)
    assert got_[:2] == ref[:2] and got_[5] is True and ref[0]
    monkeypatch.undo()
    monkeypatch.setattr(GP, 'cert_seg_ball', lambda *a: 1.0)            # a failed rigorous bound is never accepted
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31)[0] is True and GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)[0] is False
    monkeypatch.undo()
    got = []
    monkeypatch.setattr(GP, 'cert_flow', lambda eta, z0: got.append(z0))
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31)[0] is False and got == [GP.cert_b(1.955)]
    monkeypatch.undo()
    full = GP.cert_flow(GP.cert_box_disc(1.95, 1.955, 0.3, 0.31)[0], GP.cert_b(1.955))
    j = next(i for i, z in enumerate(full) if abs(1 - z) <= GP.CERT_EPS0) - 1
    lo, hi = 0.0, 1.0                                                   # the real path, cut at |1 − z| = 1.5·CERT_EPS0 on its segment
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if abs(1 - (full[j] + m * (full[j + 1] - full[j]))) > 1.5 * GP.CERT_EPS0 else (lo, m)
    cut = full[:j + 1] + [full[j] + lo * (full[j + 1] - full[j])]
    assert GP.CERT_EPS0 < abs(1 - cut[-1]) < 1.6 * GP.CERT_EPS0 and GP.cert_poly(cut, *GP.cert_box_disc(1.95, 1.955, 0.3, 0.31), 1.955)[0]
    monkeypatch.setattr(GP, 'cert_flow', lambda eta, z0: cut)
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31)[0] is False              # rejected by the end check alone
    monkeypatch.undo()
    seen = []

    def fake(R1, R2, f1, f2, ball=False):
        seen.append((R1, R2, f1, f2))
        return ((R2 - R1) < 0.01 and R2 * (f2 - f1) < 0.1, 1, -1.0, -1.0, 0.0, False) + ()
    monkeypatch.setattr(GP, 'cert_box', lambda *a, **k: fake(*a[:4]))
    bb = GP.cert_band(1.9, 1.95, 2, 3)
    ok = [x for x in seen if (x[1] - x[0]) < 0.01 and x[1] * (x[3] - x[2]) < 0.1]
    assert bb['boxes'] == len(ok) > 6 and abs(sum((x[1] - x[0]) * (x[3] - x[2]) for x in ok) - 0.05 * math.pi / 2) < 1e-15
    monkeypatch.undo()
    assert GP.CERT_BANDS[0][0] == GP.DOM_BB[-1] and all(a[1] == c[0] for a, c in zip(GP.CERT_BANDS, GP.CERT_BANDS[1:]))
    band = dict(R1=1.9, R2=1.902, boxes=1, segs=1, worst=-1.0, fails=0, side=1e-3, end=-1.0)
    gap = '\n'.join(GP.v_lines(dict(bands=[band], tb=0.9, eps0=0.02, slack=1e-12), dict(dom=dict(bb=[dict(R=1.9 - 1e-6, n=1, worst=-1.0, open=0)]))))
    rb_, rg_ = (math.floor(GP.r_star(x) * 1e4) / 1e4 for x in (1.902, 1.9 - 1e-6))
    assert rb_ != rg_ and f"r < {rg_:.4f}" in gap and f"r < {rb_:.4f}" not in gap        # the band not contiguous: not counted
    s = '\n'.join(GP.v_lines(dict(bands=[dict(band, R1=1.9, R2=1.99)], tb=0.9, eps0=0.02, slack=1e-12), dict(dom=dict(bb=[dict(R=1.9, n=1, worst=-1.0, open=0)]))))
    assert f"r < {math.floor(GP.r_star(1.99) * 1e4) / 1e4:.4f}" in s and "r < 0.9801" in s
    assert GP.ratio_floor(0.99) == '199.00' and GP.ratio_floor(0.8052) == '9.26' and (1 + 0.99) / (1 - 0.99) < 199
    bb = GP.cert_band(1.99, 1.9902, 1, 4, ball=True)                    # the ball recheck pinned on the thin slice
    assert (bb['boxes'], bb['segs'], bb['ball_short']) == (317, 105839, 108) and bb['ball_worst'] == -2.18741010735331e-12
    assert bb['ball_excess'] == 2.708944180085382e-14 and bb['ball_end'] < -1.35
    cb = GP.cert_b
    monkeypatch.setattr(GP, 'cert_b', lambda R2: 1.2 * (2 - R2) / R2)    # b outside the unit disc: rejected in ball arithmetic
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31)[0] is True and GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)[0] is False
    monkeypatch.setattr(GP, 'cert_b', cb)
    vb = dict(bands=[dict(band, R1=1.9, R2=1.99, ball_worst=-1e-12, ball_excess=0.0, ball_short=0, ball_end=-1.0)], tb=0.9, eps0=0.02,
              slack=1e-12, ball=dict(prec=120, flint='0.9.0'))
    Dd = dict(dom=dict(bb=[dict(R=1.9, n=1, worst=-1.0, open=0)]))
    assert "the valley condition certified in ball arithmetic" in '\n'.join(GP.v_lines(vb, Dd))
    vb['bands'][0]['ball_end'] = 1.0
    assert "modulo both certificates" in '\n'.join(GP.v_lines(vb, Dd))         # round 148: the report says "certificates"


def test_round108_d_in_ball_arithmetic_the_end_piece_and_the_thresholds(monkeypatch):
    """Referee round 108 (D3dc-1 … D3dc-3): D's branch and bound rechecked in ball arithmetic (`box_bound_ball` an upper bound of
    `box_bound` within 10⁻¹², pinned at 1.78 and in the raw at 1.8927; a failed rigorous bound never accepted); the rigorous end piece
    within 10⁻¹⁴ above its float twin and pinned; a rigorous segment bound of +10⁻¹⁵ rejected; b at (R2/2)(1 + t) = 1.0005 rejected;
    contiguity exact (a 10⁻¹³ gap rejected); the floors' c_need inside their radii."""
    import math
    rng = np.random.default_rng(108)
    for _ in range(200):
        R1 = rng.uniform(1.0, 1.89)
        R2, f1 = R1 + rng.uniform(1e-4, 5e-3), rng.uniform(0, 1.5)
        f2 = min(f1 + rng.uniform(1e-4, 5e-2), math.pi / 2)
        u, b = GP.box_bound_ball(R1, R2, f1, f2), GP.box_bound(R1, R2, f1, f2)
        assert u is None or (b - 1e-12 <= u and (b == -math.inf or u <= b + 1e-9 or R1 < 1 + 1e-12))
    assert GP.dom_bb(1.78, ball=True) == (4935, -3.281257330289833e-07, 0, -3.2812573349693787e-07)
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    assert [(x['R'], x['n'], x['open']) for x in raw['DB']] == [(1.78, 4935, 0), (1.8927, 402664, 0)]
    assert raw['DB'][1]['ball_worst'] == -1.9432017877494493e-11
    monkeypatch.setattr(GP, 'box_bound_ball', lambda *a: 1.0)
    n, w, op, bw = GP.dom_bb(1.2, ball=True)
    assert n == 0 and op > 0
    monkeypatch.undo()
    for R in (1.9, 1.98, 1.995):
        assert GP.cert_end(R) <= GP.cert_end_ball(R) <= GP.cert_end(R) + 1e-14
    assert GP.cert_band(1.99, 1.9902, 1, 4, ball=True)['ball_end'] == -1.3566268758611253
    assert GP.cert_band(*GP.CERT_BANDS[0][:2], 1, 1)['R2'] == GP.CERT_BANDS[0][1]
    monkeypatch.setattr(GP, 'cert_seg_ball', lambda *a: 1e-15)          # a rigorous bound just above 0: rejected
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)[0] is False
    monkeypatch.undo()
    t_out = (2 * 1.0005 / 1.955 - 1)                                     # (R2/2)(1 + t) = 1.0005: b just outside the disc
    monkeypatch.setattr(GP, 'cert_b', lambda R2: t_out)
    assert GP.cert_box(1.95, 1.955, 0.3, 0.31, ball=True)[0] is False
    monkeypatch.undo()
    band = dict(R1=1.9 + 1e-13, R2=1.99, boxes=1, segs=1, worst=-1.0, fails=0, side=1e-3, end=-1.0)
    s = '\n'.join(GP.v_lines(dict(bands=[band], tb=0.9, eps0=0.02, slack=1e-12), dict(dom=dict(bb=[dict(R=1.9, n=1, worst=-1.0, open=0)]))))
    assert "holds on |η| ≤ 1.9 (" in s                                   # a 10⁻¹³ gap: the band not counted
    for r_, R in ((0.99, 1.995), (0.8052, 1.8927), (0.0966, 1.0)):
        assert GP.c_need(r_) < R - 1e-7



def test_round110_box_bound_ball_against_a_dense_sup_the_guards_and_c_need_in_balls(monkeypatch):
    """Referee round 110 (D3de-1, -4, -6): `box_bound_ball` is at least the dense maximum of its own h (the box's float parameters,
    4·10⁵ points on (0, wm]) on wide boxes and at |η| ≈ 1.96–1.99, arg η ≈ 0.64–0.78, where the maximiser is often the quadratic piece's
    root (B(1 − w) ≤ 1) — the referee's box among them (0.349033); `part_DB` stores `dom_bb`'s rigorous bound, not the float one; a
    branch and bound with an unresolved box or without `ball_worst` falls back in `d_lines` and `v_lines` (no radius from the grid, no
    "both … in ball arithmetic"); `report` keeps D's own bb when it carries `ball_worst`; c_need at the printed floors in Arb strictly
    inside their radii."""
    import math
    from flint import arb, ctx

    def dense(R1, R2, f1, f2, n=400001):
        p = min(R * R * math.cos(2 * f) for R in (R1, R2) for f in (f1, f2))
        B, wm = R2 * math.cos(f1), 1 - 1 / R2
        w = np.linspace(wm * 1e-6, wm, n)
        a = B * (1 - w)
        h = np.log(w) + p / 2 * (w * w - 0.5) + np.where(np.abs(a) <= 1, a * a, 2 * np.abs(a) - 1) - 0.5 - math.log(R1) + math.log(4)
        k = int(np.argmax(h))
        return float(h[k]), float(a[k])
    boxes = [(1.8080676209435946, 1.821991331786887, 0.9495496783755496, math.pi / 2)]
    rng = np.random.default_rng(110)
    for _ in range(40):
        R1, f1 = rng.uniform(1.96, 1.985), rng.uniform(0.64, 0.77)
        boxes.append((R1, R1 + rng.uniform(1e-4, 5e-3), f1, f1 + rng.uniform(1e-4, 1e-2)))
    for _ in range(20):
        R1, f1 = rng.uniform(1.2, 1.9), rng.uniform(0.0, 1.4)
        boxes.append((R1, R1 + rng.uniform(0.01, 0.08), f1, min(f1 + rng.uniform(0.05, 0.6), math.pi / 2)))
    quad = 0
    for bx in boxes:
        u = GP.box_bound_ball(*bx)
        d, a = dense(*bx)
        assert u is not None and d <= u + 1e-14, (bx, u, d)          # round 115: the float noise is 1.2·10⁻¹⁶
        quad += abs(a) <= 1
    assert quad >= 10 and abs(GP.box_bound_ball(*boxes[0]) - 0.349033) < 1e-6
    monkeypatch.setattr(GP, 'DOM_BB', (1.78,))
    assert GP.part_DB() == [dict(R=1.78, n=4935, worst=-3.281257330289833e-07, open=0, ball_worst=GP.dom_bb(1.78, ball=True)[3])]
    assert GP.part_DB()[0]['ball_worst'] != GP.part_DB()[0]['worst']
    monkeypatch.undo()
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    D = raw['D']
    bbo = [dict(x, open=3) for x in raw['DB']]                          # an unresolved box at every radius
    s = '\n'.join(GP.d_lines(dict(D, dom=dict(D['dom'], bb=bbo)), raw.get('S'), raw.get('R')))
    assert "no radius beyond |η| ≤ 1 is verified" in s and "lies within the verified radius" not in s
    s = '\n'.join(GP.v_lines(raw['V'], dict(D, dom=dict(D['dom'], bb=bbo))))
    assert "gives (iii)(b) no radius here" in s and "holds on |η| ≤" not in s
    bbf = [{k: v for k, v in x.items() if k != 'ball_worst'} for x in raw['DB']]   # float only
    s = '\n'.join(GP.v_lines(raw['V'], dict(D, dom=dict(D['dom'], bb=bbf))))
    assert "modulo the branch and bound, the valley condition certified in ball arithmetic" in s and "both certificates in ball" not in s
    assert ", not interval arithmetic" in '\n'.join(GP.d_lines(dict(D, dom=dict(D['dom'], bb=bbf))))
    bbp = [dict(x, ball_worst=1e-15) for x in raw['DB']]                # a rigorous bound above 0: not "both"
    assert "both certificates in ball" not in '\n'.join(GP.v_lines(raw['V'], dict(D, dom=dict(D['dom'], bb=bbp))))
    s0 = '\n'.join(GP.v_lines(raw['V'], dict(D, dom=dict(D['dom'], bb=list(raw['DB'])))))
    assert "both certificates in ball arithmetic: referee round 108; since round 148 cross-checks of the hand proof" in s0   # round 148
    own = [dict(x, n=x['n'] + 1) for x in raw['DB']]                    # D's own bb in balls: kept, not replaced by DB
    md = GP.report(D=dict(D, dom=dict(D['dom'], bb=own)), DB=raw['DB'])
    assert "(4936 boxes" in md and "(4935 boxes" not in md
    md = GP.report(D=D, DB=raw['DB'])                                    # D's stored bb in float: replaced
    assert "(4935 boxes" in md and "in ball arithmetic -3.2813e-07" in md
    ctx.prec = 200
    for r_, R in (('0.99', 1.995), ('0.8052', 1.8927), ('0.0966', 1.0)):
        uu = 1 / arb(r_)
        assert (4 * uu.log() / (uu - 1)).sqrt() < arb(R)                 # rigorous: the whole ball below R
        assert GP.c_need(float(r_)) <= GP.c_need_ball(float(r_)) < R
    plant = {0.99: 1.9999, 0.8052: 1.8928, 0.0966: 1.5}                  # round 115: just outside each radius, below the wrong ones
    monkeypatch.setattr(GP, 'c_need_ball', lambda r: plant[round(r, 4)])  # a floor outside its radius: said so
    s = '\n'.join(GP.d_lines(D, raw.get('S'), raw.get('R'))) + '\n'.join(GP.v_lines(raw['V'], dict(D, dom=dict(D['dom'], bb=list(raw['DB'])))))
    assert s.count("NOT inside") == 2 and "NOT below 1" in s


def test_round115_part_vc_recomputed_the_regimes_and_dom_bb_s_acceptance(monkeypatch):
    """Referee round 115 (D3dj-1, -2, -6): part VC recomputed equals the raw; φ_dom the root of 8φ³ − 12φ² + 1 to 10⁻⁹; Ψ¹_φ's closed
    ratio against the ball coefficients; the far-out floor positive from φ = 25 at ε = 10⁻³ while the least stays negative; Ψ̃_ψ against
    `cert_E`/ε²; the θ ~ √ε path negative and −ψ² − ψ⁴/2 far out; the regimes' matching limits; `dom_bb` rejecting a positive ball bound."""
    import math
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    assert GP.part_VC() == raw['VC']
    lo, hi = 1.0, 2.0
    for _ in range(80):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if 8 * m ** 3 - 12 * m * m + 1 < 0 else (lo, m)
    assert abs(raw['VC']['phi_dom'] - lo) < 1e-9 and abs((lo - 1) ** 2 * (2 * lo + 1) - 0.75) < 1e-12
    for c in raw['VC']['coef']:
        assert abs(c['ball'][1] / c['closed'] - 1) < 1e-6
        w = complex(c['x'], 2 * c['phi'] * c['x'])
        assert abs(c['closed'] + (0.25 + 5 * c['phi'] ** 2 + 8 * c['phi'] ** 4 - c['x'] ** 2 * ((1 + 2j * c['phi']) ** 4).real / 2)
                   / (1 + 4 * c['phi'] ** 2)) < 1e-9 * (1 + abs(c['closed']))
    far = raw['VC']['far']
    assert min(x['phi'] for x in far if x['floor'] > 0) == 25.0 and all(x['least'] < 0 for x in far)
    assert all(abs(m['vals'][-1] - m['lim']) < 1e-5 * (1 + abs(m['lim'])) for m in raw['VC']['mid'])
    assert raw['VC']['path']['worst'] < 0 and raw['VC']['path']['far'] < 1e-9
    pth = raw['VC']['path']                                              # round 117 (D3dl-4): the grid's range in the raw
    assert (pth['xmin'], pth['xmax'], pth['xmax0'], pth['n_psi'], pth['psi']) == (1e-4, 1602.0, 200.0, 145, (0.02, 20.0))
    pm = [m for m in raw['VC']['mid'] if m['P'] > 0]                     # the (Re σ − 1)₊² branch exercised
    assert len(pm) == 1 and abs(pm[0]['vals'][-1] - pm[0]['lim']) < 1e-6 and abs(pm[0]['lim'] + 0.005825) < 1e-9
    assert all(m['eps'] == list(GP.VC_MID_EPS) for m in raw['VC']['mid']) and all(c['eps'] == ['1e-6', '1e-8'] for c in raw['VC']['coef'])
    for q in (0.3, 4.0):                                                  # inside the disc the path's closed form
        xs = np.linspace(0.01, 0.3, 30)
        y = 2 * q / (1 + 2 * q * q) * xs
        ins = xs - 0.5 + y * y / 2 <= 0
        assert np.allclose(GP.corner_mid_path(q, xs)[ins], (-4 * q * q * xs ** 2 / (1 + 2 * q * q) - y ** 4 / 2)[ins], rtol=1e-12, atol=1e-15)
    from fractions import Fraction as Fr                                 # Lemma VC off the grid, in exact rationals (ψ, x log-uniform)
    def pt_exact(q, x):
        y = min(2 * q * x / (1 + 2 * q * q), q)
        D, P = x - Fr(1, 2) + y * y / 2, x - Fr(1, 2) - q * q / 2 - q * y
        return (1 + 2 * q * q) * y * y - 4 * q * x * y - y ** 4 / 2 + max(D, 0) ** 2 - max(P, 0) ** 2
    rng = np.random.default_rng(115)
    q, xx = np.exp(rng.uniform(np.log(1e-3), np.log(1e2), 3000)), np.exp(rng.uniform(np.log(1e-3), np.log(1e5), 3000))
    assert all(pt_exact(Fr(a), Fr(b)) < 0 for a, b in zip(q, xx))
    assert all(abs(float(GP.corner_mid_path(a, np.array([b]))[0]) - float(pt_exact(Fr(a), Fr(b)))) < 1e-9 * (1 + b * b)
               for a, b in zip(q[:300], xx[:300]) if b < 1e3)
    pinf = lambda x, t: t * t - 4 * x * t + max(x - 0.5, 0.0) * (t + 1) ** 2
    for x, t in ((0.2, 0.3), (0.9, 1.0), (3.0, 0.5)):
        assert abs(GP.corner_psi(x, 1e4 * t, 1e4) / 1e8 - pinf(x, t)) < 1e-3
        assert abs(float(GP.corner_psi_mid(x, 1e-4 * t, 1e-4)) / 1e-8 - pinf(x, t)) < 1e-3
        assert abs(float(GP.corner_psi_mid(x, t / 1e4, 1e4)) - (2 * t * t - 4 * x * t + max(x - 0.5, 0.0) ** 2)) < 1e-3
    monkeypatch.setattr(GP, 'box_bound_ball', lambda *a: 5e-10)          # a rigorous bound in (0, 10⁻⁹): never accepted
    n, w, op, bw = GP.dom_bb(1.2, ball=True, max_boxes=200)
    assert n == 0


def test_part_vc_the_corner_s_blow_up_limit():
    """Part VC (the corner at η = 2): Ψ_φ is `cert_E`'s ε³ limit (the ratio within 1% at ε = 10⁻⁴ on the floor for φ ≤ 2); the floor is
    Ψ's minimiser in y with the stated values; φ_dom ≈ 1.4397 and the measured sector over ε approaches it; the raw holds the run."""
    import math
    import cmath
    for phi in (0.0, 0.7, 2.0):
        for x in (0.1, 0.4, 0.8, 3.0):
            y = GP.corner_floor(x, phi)
            v = GP.corner_psi(x, y, phi)
            exp_ = -(1 + 4 * phi * phi) * x * x if x <= 0.5 else -x * x - phi * phi - (x - 0.5) ** 3 / 3
            assert abs(v - exp_) < 1e-12 and all(GP.corner_psi(x, y + d, phi) > v for d in (-1e-3, 1e-3))
            eps = 1e-4
            E_ = GP.cert_E(eps * complex(x, y), (2 - eps) * cmath.exp(1j * eps * phi))
            assert abs(E_ / (eps ** 3 * v) - 1) < 1e-2
    pd, g0 = GP.corner_phi_dom()
    assert abs(pd - 1.4397) < 1e-3 and abs(g0 + 0.25) < 1e-12
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    lo = [x['lo'] / (2 - x['R']) for x in raw['VC']['sector']]
    assert lo == sorted(lo, reverse=True) and abs(lo[-1] - pd) < 0.01 and all(p_['E'] < 0 for r in raw['VC']['rows'] for p_ in r['pts'])
    s = '\n'.join(GP.vc_lines(raw['VC']))
    assert "the bound negative at every point" in s and "φ_dom = 1.4397" in s and GP.vc_lines(None) == []
    bad = dict(raw['VC'], rows=[dict(r, pts=[dict(p_, E=abs(p_['E'])) for p_ in r['pts']]) for r in raw['VC']['rows']])
    assert "NOT negative at every point" in '\n'.join(GP.vc_lines(bad))     # V05: the negativity checked


def test_part_va_the_last_annulus_certificate(monkeypatch):
    """Part VA (referee round 120, adopted; round 121's fix applied): the certificate recomputed equals the raw (counts, bounds,
    validation), no piece fails, the end piece is negative, the forms agree with `cert_E` and F_η's definition at the validation's
    points; round 121 (D3dp-1, -2, -4): the outer piece's 76 η-boxes and 5193 segments (the connector of the η-box at arg η ≤ π/128
    next to |η| = 2 no longer dropped), the joint checked at every box and the boxes tiling the sector, the true value below the
    bound in every accepted box; negative controls, each run beside its positive twin in the same harness (D3dp-3) — the closed form
    raised by 0.05, the cap's derivative raised —; plants for the joint (a NaN connector, a connector missing the inner path's end)
    and for the soundness check (a bound below the true value)."""
    import math
    import gh_lastann as LA
    from flint import arb
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    va = GP.part_VA()
    assert va == list(raw['VA'])
    by = {r['piece']: r for r in va}
    assert all(by[k]['failed'] == 0 and by[k]['boxes'] > 0 for k in ('A-J', 'A-FC', 'B-J', 'B-F', 'O')) and by['end']['ok']
    assert by['validate']['A_vs_def'] < 1e-20 and by['validate']['B_vs_centred'] < 1e-20 and by['validate']['A0'] < 1e-40
    assert by['validate']['B0'] < 1e-40 and by['validate']['counts']['points'] > 1000
    assert (by['O']['boxes'], by['O']['segs']) == (76, 5193) and by['O']['R1'] == 1.995
    h = by['handoff']
    assert h['failed'] == 0 and h['tiled'] and h['boxes'] == 76 and h['points'] == 76 * 15 and h['min_margin'] >= 0
    assert h['joints'] and h['inside'] and h['disjoint'] and h['sub_boxes'] >= 76 and 'joints' not in by['O']   # round 122
    assert by['meta'] == dict(piece='meta', prec=LA.PREC, delta0=LA.DELTA0, eps0=LA.EPS0)
    assert all(by['sound ' + k]['failed'] == 0 and by['sound ' + k]['worst'] < 0 for k in ('A-J', 'A-FC', 'B-J', 'B-F'))
    assert [by['sound ' + k]['points'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [240, 18326, 3056, 7434]   # round 127: distinct points
    assert [by['sound ' + k]['evals'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [240, 24770, 3056, 7434]    # point–branch evaluations
    assert h['area_exact']
    nan_box = [x for x in by['O']['accepted'] if x[0] == 1.995 and x[2] == 0.0]   # D3dp-1: the box whose κ ball is NaN
    assert len(nan_box) == 1 and nan_box[0][3] == math.pi / 256 and nan_box[0][4] == 0.0 and nan_box[0][5] > 0.0122
    assert LA.handoff_range(1.995, 2.0, 0.0, math.pi / 256, LA.DELTA0)[:2] == (0.0, nan_box[0][5])   # the NaN branch: [0, θ_hi]
    assert LA.handoff_range(1.995, 2.0, 0.0, math.pi / 128, LA.DELTA0)[:2] == (0.0, LA.up(LA.th_ball(0.0, math.pi / 128)))
    O = dict(by['O'], accepted=[list(x) for x in by['O']['accepted']])            # plants for the joint
    O['accepted'][0][4] = float('nan')
    assert LA.handoff_check(O)['failed'] == 1
    O = dict(by['O'], accepted=[list(x) for x in by['O']['accepted']])
    i = next(j for j, x in enumerate(O['accepted']) if x[0] == 1.995 and x[2] == 0.0)
    O['accepted'][i][5] = 1e-4                                                     # a connector short of the inner path's end
    assert LA.handoff_check(O)['failed'] > 0
    O = dict(by['O'], accepted=[x for x in by['O']['accepted'] if not (x[0] == 1.995 and x[2] == 0.0)])
    assert not LA.handoff_check(O)['tiled']                                        # a box missing: the sector not tiled
    acc = [list(x) for x in by['O']['accepted']]                                   # round 122 (D3dq-2, T01): a lower θ-half for
    j = next(j for j, x in enumerate(acc) if x[0] == 1.995 and x[2] == 0.0)        # its upper twin, the area kept: not disjoint
    k = next(k for k, x in enumerate(acc) if x[2] == acc[j][3] and x[0] == acc[j][0])
    acc[k] = list(acc[j][:2]) + [acc[j][2], acc[j][2] + (acc[k][3] - acc[k][2])] + acc[k][4:]
    t = LA.handoff_check(dict(by['O'], accepted=acc))
    assert not t['tiled'] and not t['disjoint'] and abs(t['area_rel']) < 1e-12
    import copy                                                                    # round 122 (D3dq-1): the joint from the pieces
    Oj = LA.piece_outer()                                                          # certified, in memory
    assert LA.handoff_check(Oj)['failed'] == 0 and LA.handoff_check(Oj)['joints']
    plant = lambda f: (lambda O_: (f(O_['joints'][j]), LA.handoff_check(O_))[1])(copy.deepcopy(Oj))
    assert 'connector cover' in str(plant(lambda J: J.update(cert=J['cert'][1:]))['fails'])            # O02: its first piece missing
    assert 'polyline start' in str(plant(lambda J: J.update(pts=[J['pts'][0] + 1e-3j] + J['pts'][1:]))['fails'])   # H10
    assert 'y_c' in str(plant(lambda J: J.update(yc=J['yc'] + 1.0))['fails'])
    def skip(J):                                                                   # a shortcut past one of the polyline's vertices
        pc, v = J['cert'], set(J['pts'][1:-1])
        i = next(i for i in range(len(pc) - 1) if pc[i][1] in v and pc[i][0].real != LA.DELTA0)
        J.update(cert=pc[:i] + [(pc[i][0], pc[i + 1][1])] + pc[i + 2:])
    assert 'polyline cover' in str(plant(skip)['fails'])
    nt = [dict(r, tiled=False) if r['piece'] == 'handoff' else r for r in va]       # round 122 (D3dq-2, W01): no conclusion untiled
    nj = [dict(r, joints=False) if r['piece'] == 'handoff' else r for r in va]      # round 126 (W01): nor with the joints unchecked
    assert "NOT every piece passes" in '\n'.join(GP.va_lines(nj))
    # round 126 (D3du-1 … -4): the failure paths of the new checks
    acc = [list(x) for x in by['O']['accepted']]
    acc[1][2] += 1e-14                                                             # T1: a 1e-14 gap between two θ-boxes
    t = LA.handoff_check(dict(by['O'], accepted=acc))
    assert t['inside'] and t['disjoint'] and t['area_rel'] < 1e-12 and not t['area_exact'] and not t['tiled']
    acc = [list(x) for x in by['O']['accepted']]
    acc[-1][1] = 2.001                                                             # a box outside the sector
    assert not LA.handoff_check(dict(by['O'], accepted=acc))['inside']
    assert LA._y_enclosed(1.995, 2.001, 0.0, 0.01, 0.0, 0.02, LA.DELTA0) is None    # outside the clamps' domain: not proved
    def cut(J):                                                                    # P8: the polyline cut halfway (a chain still)
        m = len(J['pts']) // 2
        pc = J['cert']
        i = next(i for i, x in enumerate(pc) if x[1] == J['pts'][m] and x[0].real != LA.DELTA0)
        J.update(pts=J['pts'][:m + 1], cert=pc[:i + 1])
    assert 'polyline end' in str(plant(cut)['fails'])
    nc = lambda J: sum(1 for x in J['cert'] if x[0].real == LA.DELTA0 and x[1].real == LA.DELTA0 and x[0].imag >= J['lo'] and x[1].imag <= J['hi'])
    assert 'connector cover' in str(plant(lambda J: J.update(cert=J['cert'][:nc(J) - 1] + J['cert'][nc(J):]))['fails'])   # its last piece
    assert 'connector cover' in str(plant(lambda J: J.update(cert=J['cert'][nc(J):]))['fails'])   # round 120's defect: none at all
    def mid(J):                                                                    # a middle half-piece of the polyline dropped
        k = nc(J) + (len(J['cert']) - nc(J)) // 2
        J.update(cert=J['cert'][:k] + J['cert'][k + 1:])
    assert 'polyline cover' in str(plant(mid)['fails'])
    import gh_lastann as LA2
    orig = LA2._y_enclosed
    monkeypatch.setattr(LA2, '_y_enclosed', lambda *a, **k: None)                  # the sub-proof failing: reported
    assert LA.handoff_check(Oj)['failed'] == 76
    monkeypatch.setattr(LA2, '_y_enclosed', orig)
    b2 = next(x for x in by['O']['accepted'] if x[2] > 0 and x[4] == x[2])         # lo = arg η's lower end (the cap):
    assert LA._y_enclosed(*b2[:4], b2[4] + 1e-9, b2[5], LA.DELTA0) is None          # raised above inf Y: not proved
    # round 127 (D3dv-3): plants at the thresholds
    def p7(J):                                                                     # the end 0.03 from z = 1, just past END_R
        pc = list(J['cert'])
        pc[-1] = (pc[-1][0], 0.97 + 0j)
        J.update(cert=pc, pts=J['pts'][:-1] + [0.97 + 0j])
    assert 'polyline end' in str(plant(p7)['fails'])
    acc = [list(x) for x in by['O']['accepted']]                                   # a bottom-row box shifted by 2⁻¹⁴ in |η|: the
    i0 = next(i for i, x in enumerate(acc) if x[0] == 1.995 and x[1] == 2.0)       # area exact, disjoint, but outside the sector
    acc[i0][0], acc[i0][1] = acc[i0][0] - 2.0 ** -14, acc[i0][1] - 2.0 ** -14
    t = LA.handoff_check(dict(by['O'], accepted=acc))
    assert t['area_exact'] and t['disjoint'] and not t['inside'] and not t['tiled']
    for args in ((1.995, 2.0, -1e-3, 0.01), (1.995, 2.0, 1.5, math.pi / 2 + 1e-3), (0.0, 1.0, 0.0, 0.01)):   # outside the clamps'
        assert LA._y_enclosed(*args, 0.0, 10.0, LA.DELTA0) is None                                      # domain, every side
    b0 = max(by['O']['accepted'], key=lambda x: LA._y_enclosed(*x, LA.DELTA0))     # the box needing the most sub-boxes
    assert LA._y_enclosed(*b0, LA.DELTA0) > 4 and LA._y_enclosed(*b0, LA.DELTA0, budget=4) is None    # the budget binds
    # round 127 (D3dv-2): the whole-disc conclusion only with the seams to V and D's branch and bound
    V, D = raw['V'], raw['D']
    D = dict(D, dom=dict(D['dom'], bb=[dict(x) for x in raw['DB']]))
    full = '\n'.join(GP.va_lines(va, V, D))
    assert "on the whole open disc |η| < 2" in full and "meeting VA at exactly its R1" in full
    seamless = lambda VA_, V_, D_: "on the whole open disc |η| < 2" not in '\n'.join(GP.va_lines(VA_, V_, D_))
    bands = [dict(b) for b in V['bands']]
    assert seamless(va, dict(V, bands=bands[:-1] + [dict(bands[-1], fails=1)]), D)       # V's last band failed
    assert seamless(va, dict(V, bands=bands[:-1] + [dict(bands[-1], R2=1.994)]), D)      # V ending at 1.994
    assert seamless(va, V, dict(D, dom=dict(D['dom'], bb=[dict(x, open=1) for x in D['dom']['bb']])))   # D unresolved
    va_r = [dict(r, R1=1.996) if r['piece'] == 'O' else dict(r, eps0=2 - 1.996) if r['piece'] == 'meta' else r for r in va]
    assert seamless(va_r, V, D) and "holds on 1.996 ≤ |η| < 2" in '\n'.join(GP.va_lines(va_r, V, D))   # VA starting at 1.996
    va_r = [dict(r, R1=1.996) if r['piece'] == 'O' else r for r in va]                   # round 131: R1 against meta's ε₀
    assert seamless(va_r, V, D) and "NOT complete (meta)" in '\n'.join(GP.va_lines(va_r, V, D))
    ne = [dict(r, ok=False) if r['piece'] == 'end' else r for r in va]                  # G01: the end piece failing
    assert "NOT every piece passes" in '\n'.join(GP.va_lines(ne, V, D))
    assert seamless(va, None, None)                                                       # no V, no D: the annulus only
    top = next(x for x in by['O']['accepted'] if x[3] == math.pi / 2)              # D3du-4: Y = 0 at the true π/2, so no lo > 0
    assert LA._y_enclosed(top[0], top[1], math.pi / 2 - 1e-12, math.pi / 2, 1e-20, 1.0, LA.DELTA0) is None
    assert LA._y_enclosed(top[0], top[1], math.pi / 2 - 1e-12, math.pi / 2, 0.0, 1.0, LA.DELTA0) is not None
    assert "NOT every piece passes" in '\n'.join(GP.va_lines(nt)) and "NOT every piece passes" not in '\n'.join(GP.va_lines(va))
    assert LA._y_enclosed(1.995, 2.0, 0.0, math.pi / 256, 0.0, 0.01, LA.DELTA0) is None   # a connector short of Y: not proved
    kp = []
    LA.pieceA('J', keep=kp)
    assert LA.soundness('A-J', LA.fA, kp)['failed'] == 0
    low = [(iv, [(br, 'v', v - 1.0) for br, _, v in recs]) for iv, recs in kp]    # a bound below the true value
    assert LA.soundness('A-J', LA.fA, low)['failed'] > 0
    kpf = []
    LA.pieceA('FC', keep=kpf)                                                      # round 122 (S03): the excess there is −2.4e-13,
    low = [(iv, [(br, ro, v - 1e-5 if ro == 'v' else v) for br, ro, v in recs]) for iv, recs in kpf]   # so 1e-5 lower is caught
    assert LA.soundness('A-FC', LA.fA, low)['failed'] > 0
    dkeep = [(iv, [r for r in recs if r[1] == 'd']) for iv, recs in kpf if any(r[1] == 'd' for r in recs)]
    assert dkeep and LA.soundness('A-FC', LA.fA, dkeep)['failed'] == 0              # route 'd' (S02, S06, K02): a form just above 0
    assert LA.soundness('A-FC', lambda *a, **k: arb(1e-300), dkeep)['failed'] > 0    # there is caught
    import functools
    monkeypatch.setattr(LA, 'bb', functools.partial(LA.bb, max_boxes=3000))   # the controls' harness: a failing control stops early
    assert LA.pieceA('J', log=None)['failed'] == 0                      # the positive twin passes in it
    f0 = LA.fA0
    monkeypatch.setattr(LA, 'fA0', lambda *a, **k: f0(*a, **k) + 0.05)
    assert LA.pieceA('J', log=None)['failed'] > 0                       # a closed form 0.05 too high: the junction no longer passes
    monkeypatch.undo()
    monkeypatch.setattr(LA, 'bb', functools.partial(LA.bb, max_boxes=3000))
    assert LA.pieceA('FC', log=None)['failed'] == 0                     # the positive twin
    fa = LA.fA
    monkeypatch.setattr(LA, 'fA', lambda u, q, rho, br, lam=None: fa(u, q, rho, br, lam) + (rho * 10 if br == 'C' else 0))
    assert LA.pieceA('FC', log=None)['failed'] > 0                      # the cap's derivative raised: the matching corner fails
    md = open(os.path.join(GP.OUT, 'gh_peano.md'), encoding='utf-8').read()
    assert "## VA. The last annulus certified in ball arithmetic (referee rounds 120–122, 126, 127, 131)" in md and "pending round 121" not in md
    assert "proved inside it in balls on 484 sub-boxes" in md and "at 128 bits" in md and "— 0 failures —" in md
    assert "the outer piece 76 η-boxes and 5193 segments, 0 failed" in md and "round 120's certificate was not a proof as submitted" in md
    assert "that box now split in two, the lower at arg η ≤ 0.01227 with its connector [0, 0.01227] certified" in md


def test_part_va_series_against_arb():
    """Round 121 (D3dp-4, the mutants S05, G12): `gh_lastann`'s series with rigorous tails — S = sinc√w, Cv = (1 − cos√w)/w,
    Λ = (v − log(1 + v))/v², H₁ = (h(1 + δ) − δ²)/δ³ — contain Arb's own values and derivatives, at the extremes of their use (w up
    to π², |v|, |δ| up to 0.4999)."""
    import gh_lastann as LA
    from flint import arb
    LA.ctx.prec = LA.PREC
    sinc = lambda w: (w.sqrt().sin() / w.sqrt())
    dsinc = lambda w: (w.sqrt() * w.sqrt().cos() - w.sqrt().sin()) / (2 * w * w.sqrt())
    cv = lambda w: (1 - w.sqrt().cos()) / w
    dcv = lambda w: (w.sqrt().sin() / (2 * w.sqrt()) - (1 - w.sqrt().cos()) / w) / w
    lam = lambda v: (v - (1 + v).log()) / (v * v)
    dlam = lambda v: ((1 - 1 / (1 + v)) * v - 2 * (v - (1 + v).log())) / (v * v * v)
    h1 = lambda d: (d + d * d / 2 - (1 + d).log() - d * d) / (d * d * d)
    dh1 = lambda d: ((1 + d - 1 / (1 + d) - 2 * d) * d - 3 * (d + d * d / 2 - (1 + d).log() - d * d)) / (d ** 4)
    X = lambda *v: [arb(float(t)) for t in v]                           # exact doubles: the inputs carry no radius
    cases = [(LA.S, sinc, dsinc, X(0.01, 0.5, 1, 3, 2.4674011, 9.8696044)),   # up to π²/4 and π²
             (LA.Cv, cv, dcv, X(0.01, 1, 2.4674011, 9.8696044)),
             (LA.Lam, lam, dlam, X(-0.4999, -0.1, 0.001, 0.3, 0.4999)),
             (LA.H1, h1, dh1, X(-0.4999, -0.2, 0.01, 0.4999))]
    from flint import ctx

    def ref(f, x):                                                      # the reference at 400 bits: a ball far inside the series'
        p0 = ctx.prec
        ctx.prec = 400
        try:
            return f(x)
        finally:
            ctx.prec = p0
    for ser, ex, dex, xs in cases:
        for x in xs:
            assert ser(x).contains(ref(ex, x)), (ser, x)                # round 122 (D3dq-5): contains, not overlaps
            y = ser(LA.D(x, 1))
            assert y.v.contains(ref(ex, x)) and y.d.contains(ref(dex, x)), (ser, x)
            assert ser(x).rad() < (1e-8 if abs(float(x.mid())) > 0.4 else 1e-15)   # the tails: 5.8e-11 (Λ), 1.0e-9 (H₁) at 0.4999


def test_round120_vc_labels_from_the_raw():
    """Referee round 120 (D3do-6; round 117's mutant Y02): part VC's ε labels are the raw's, not constants — the coefficients' and Ψ̃'s
    ball-arithmetic points printed at the ε stored with them; the (Re σ > 1) label outside the tuple; the floor y = min(sx, ψ) to
    leading order only; two forms, not four charts."""
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    VC = raw['VC'].item() if hasattr(raw['VC'], 'item') else raw['VC']
    VC = dict(VC, coef=[dict(c, eps=['1e-5', '1e-7']) for c in VC['coef']], mid=[dict(m, eps=['1e-3', '1e-5', '1e-7']) for m in VC['mid']])
    t = '\n'.join(GP.vc_lines(VC))
    assert "at ε = 10⁻⁵, 10⁻⁷ against it" in t and "in ball arithmetic at ε = 10⁻³, 10⁻⁵, 10⁻⁷ " in t
    md = open(os.path.join(GP.OUT, 'gh_peano.md'), encoding='utf-8').read()
    assert "(ψ, x, y) = (0.2, 3, 0.1), Re σ > 1:" in md and "y = min(sx, ψ) to leading order (the cap y = ψ exactly" in md
    assert "Two forms suffice (referee round 120)" in md and "four parameter charts" not in md and "Open." not in md


def test_round127_duhamel_identity_moment_by_moment():
    """Referee round 127 (D3dv-1): Duhamel's identity, proved in D61 by telescoping J^m − J₀^m = Σ_j J^{m−1−j}VJ₀^j on the banded Jacobi
    matrices, holds moment by moment: (E − Q_K)[z^m] = (1/(K − 1)!)Σ_{j<m} E[Z^{m−1−j}He_K(Z)]·Q_K[z^jHe_{K−1}], E[Z^nHe_K(Z)] =
    n!/(n − K)!·(n − K − 1)!! (Gaussian integration by parts), at 60 digits for K ≤ 6 and m ≤ 3K + 4 — the identity's ξ^m/m! coefficients."""
    from decimal import Decimal, localcontext
    from math import factorial

    def dfact(n):
        return 1 if n <= 0 else n * dfact(n - 2)

    def eZnHe(n, K):
        return 0 if n < K or (n - K) % 2 else factorial(n) // factorial(n - K) * dfact(n - K - 1)
    for K in range(1, 7):
        with localcontext() as ctx:
            ctx.prec = 60
            zs, ws = WM.gh_nodes(K, 60)
            hk1 = [WM.he_all(z, K)[K - 1] for z in zs]
            qk = lambda j: sum((w * (z ** j if j else Decimal(1)) * h for z, w, h in zip(zs, ws, hk1)), Decimal(0))
            for m in range(0, 3 * K + 5):
                rhs = sum((eZnHe(m - 1 - j, K) * qk(j) for j in range(m)), Decimal(0)) / factorial(K - 1)
                assert abs(Decimal(GP.defect(m, K)) - rhs) < Decimal(10) ** -35 * (1 + abs(rhs)), (K, m)


def _seamed():
    """The raw's VA, V and D with DB substituted for D's float branch and bound (as `report` does)."""
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    D = dict(raw['D'], dom=dict(raw['D']['dom'], bb=[dict(x) for x in raw['DB']]))
    return raw, list(raw['VA']), raw['V'], D


def test_round131_va011_d_and_va_alone_cover_the_disc():
    """Referee round 131 (D3dz-3): part VA at ε₀ = 0.11 (`part_VA011`) recomputed equals the raw, with the referee's counts — the inner
    pieces 54, 2533, 413 (549 skipped), 1189 (1174 skipped), the outer piece 220 η-boxes and 16 083 segments (the largest bound
    −4.2·10⁻⁹), the end piece −1.3314, the joint on 1306 sub-boxes, soundness at 324, 22 258, 3982, 11 634 points —; its seam to D
    (R1 = 1.89 ≤ 1.8927, D in ball arithmetic) and the conclusion printed only with it; plants at the seam's threshold: R1 exactly
    1.8927 (the sound edge) and one ulp above, D's 1.8927 run open or in floats only, VA011 missing a piece."""
    import math
    import gh_lastann as LA
    raw, _, V, D = _seamed()
    va = GP.part_VA011()
    assert va == list(raw['VA011']) and GP.VA011_EPS0 == 0.11
    by = {r['piece']: r for r in va}
    assert [by[k]['boxes'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [54, 2533, 413, 1189]
    assert [by[k]['skipped'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [0, 0, 549, 1174]
    assert all(by[k]['failed'] == 0 for k in ('A-J', 'A-FC', 'B-J', 'B-F', 'O', 'handoff') + tuple('sound ' + k for k in ('A-J', 'A-FC', 'B-J', 'B-F')))
    assert (by['O']['boxes'], by['O']['segs'], by['O']['R1']) == (220, 16083, 1.89) and -5e-9 < by['O']['worst'] < -4e-9
    assert by['O']['R1'] == 2 - 0.11 and by['end']['ok'] and abs(by['end']['bound'] + 1.3314) < 1e-4
    h = by['handoff']
    assert (h['sub_boxes'], h['points'], h['boxes']) == (1306, 3300, 220) and h['tiled'] and h['joints'] and h['area_exact']
    assert [by['sound ' + k]['points'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [324, 22258, 3982, 11634]
    assert [by['sound ' + k]['evals'] for k in ('A-J', 'A-FC', 'B-J', 'B-F')] == [324, 30132, 3982, 11634]
    assert by['meta'] == dict(piece='meta', prec=LA.PREC, delta0=LA.DELTA0, eps0=0.11) and LA.EPS0 == 0.005
    v = by['validate']                                                   # D3dz-4: the bar is ε₀'s — 1.1e-12 here, above VA's 1e-20
    assert 1e-13 < v['B_vs_centred'] < 1e-11 and v['A_vs_def'] < 1e-20 and max(v['A0'], v['B0']) < 1e-40 and v['counts']['skipped'] == 2
    assert v['B_vs_centred'] < LA.validate_tol(0.11) == 1e-10 and not v['B_vs_centred'] < LA.validate_tol(0.005) == 1e-20
    assert LA.validate_tol(0.05) == 1e-10 and LA.validate_tol(0.004) == 1e-20 and LA.validate_tol(0.2) is None
    assert GP.va_complete(va) == (True, []) and GP.va011_seam(va, D) == (1.8927, 1.89, True)
    t = '\n'.join(GP.va011_lines(va, D, V))
    assert "so D and VA011 alone cover 1 < |η| < 2" in t and "V an independent cross-check on 1.8927 ≤ |η| ≤ 1.995." in t
    assert "the outer piece 220 η-boxes and 16083 segments, 0 failed" in t and "the joint proved on 1306 sub-boxes" in t
    alone = lambda VA_, D_: "alone cover 1 < |η| < 2" in '\n'.join(GP.va011_lines(VA_, D_, V))
    at = lambda R1: [dict(r, R1=R1) if r['piece'] == 'O' else dict(r, eps0=2 - R1) if r['piece'] == 'meta' else r for r in va]
    assert alone(at(1.8927), D) and GP.va011_seam(at(1.8927), D)[2]     # R1 = R_bb exactly (2 − ε₀ = 1.8927 too: Sterbenz)
    assert not alone(at(math.nextafter(1.8927, 2)), D)                   # one ulp above D's radius: a gap
    e_ = math.nextafter(2 - 1.8927, 0)                                   # R1 = fl(2 − ε₀) = 1.8927 still, but the inner pieces start
    va_e = [dict(r, R1=1.8927) if r['piece'] == 'O' else dict(r, eps0=e_) if r['piece'] == 'meta' else r for r in va]
    assert 2 - e_ == 1.8927 and Fraction(2) - Fraction(e_) > Fraction(1.8927) and not alone(va_e, D)   # at 2 − ε₀ > 1.8927: a gap
    assert "is not inside D's ball-certified disc" in '\n'.join(GP.va011_lines(at(math.nextafter(1.8927, 2)), D, V))
    assert not alone([dict(r, eps0=0.1) if r['piece'] == 'meta' else r for r in va], D)   # meta's ε₀ against R1
    bb = [dict(x) for x in D['dom']['bb']]
    assert not alone(va, dict(D, dom=dict(D['dom'], bb=[bb[0], dict(bb[1], open=1)])))    # D's 1.8927 run open: D only to 1.78
    assert not alone(va, raw['D']) and alone(va, D)                      # D in floats only (no ball recheck)
    assert not alone(va, dict(D, dom=dict(D['dom'], bb=[bb[0]])))         # D only to 1.78
    for k in GP.VA_PIECES:                                               # every piece required
        assert not alone([r for r in va if r['piece'] != k], D), k
    assert GP.va011_lines(None) == [] and GP.va011_lines([]) == []
    md = open(os.path.join(GP.OUT, 'gh_peano.md'), encoding='utf-8').read()
    assert "Part VA again at ε₀ = 0.11 (`part_VA011`, referee round 131" in md and "so D and VA011 alone cover 1 < |η| < 2" in md
    assert "V an independent cross-check on 1.8927 ≤ |η| ≤ 1.995." in md
    assert ("Since referee round 148 D, V, VA and VA011 are all independent cross-checks: the valley condition on 1 < |η| < 2 is proved "
            "by hand (review/RESEARCH_note_annulus.md, Theorem A).") in md and "and since referee round 148 an independent cross-check" in md


def test_round131_va_seams_complete_and_planted_at_the_thresholds():
    """Referee round 131 (D3dz-5, the survivors Y01, Y05, Y07, Y08, Y09): the whole-disc conclusion needs V's bands in ball arithmetic
    (`ball_worst` and `ball_end`, both defaulting to 0.0), contiguous to the ulp, ending exactly at VA's R1, D's branch and bound
    in balls, and a complete VA — every piece of VA_PIECES present and non-empty, R1 = 2 − meta's ε₀, the validation within its bar
    (`va_complete`); plants for each, beside sound twins (an overlap one ulp past 1.995; a deviation just below the bar); and
    `va_seams` on the raw's numpy array (D3dz-8)."""
    import math
    raw, va, V, D = _seamed()
    whole = lambda VA_, V_, D_: "on the whole open disc |η| < 2" in '\n'.join(GP.va_lines(VA_, V_, D_))
    assert whole(va, V, D) and GP.va_seams(raw['VA'], V, D) == (1.8927, 1.995, 1.995, True)   # the numpy array: no ValueError
    assert GP.va_seams(None, V, D)[3] is False and GP.va_complete(va) == (True, [])
    assert not GP.va_seams([r for r in va if r['piece'] != 'validate'], V, D)[3]   # va_seams itself requires a complete VA
    full = '\n'.join(GP.va_lines(va, V, D))
    assert "the limit (iii)(b) holds on the whole open disc |η| < 2 (D's domination is the valley condition on the straight path, so Ω ⊆ Ω′)" in full
    assert "the valley condition holds on 1.995 ≤ |η| < 2 (with η ↦ η̄, −η the whole annulus)" in full
    bands = [dict(b) for b in V['bands']]
    Vw = lambda bs: dict(V, bands=bs)
    nb = [dict(b) for b in bands]; del nb[1]['ball_worst']
    assert not whole(va, Vw(nb), D)                                      # Y05: a band without its ball recheck
    nb = [dict(b) for b in bands]; del nb[2]['ball_end']
    assert not whole(va, Vw(nb), D)                                      # its end piece not rechecked in balls (round 127's −1.0 passed it)
    assert not whole(va, Vw([bands[0], bands[2]]), D)                    # Y07: the middle band missing
    assert not whole(va, Vw([bands[0], dict(bands[1], R1=math.nextafter(1.98, 2)), bands[2]]), D)   # Y08: a 1-ulp gap at 1.98
    assert not whole(va, Vw(bands[:2] + [dict(bands[2], R2=math.nextafter(1.995, 0))]), D)        # Y09: V one ulp short of R1
    assert not whole(va, Vw(bands[:2] + [dict(bands[2], R2=math.nextafter(1.995, 2))]), D)        # an overlap is refused too (the
    assert whole(va, Vw(bands + [dict(bands[2], R1=1.5, R2=1.995)]), D)                            # seam is exact); an extra band: sound
    assert not whole(va, V, raw['D'])                                    # Y01: D in floats only (no ball_worst)
    assert not whole(va, V, dict(D, dom=dict(D['dom'], bb=[dict(x, ball_worst=1e-12) for x in D['dom']['bb']])))
    for k in GP.VA_PIECES:                                               # a VA raw missing any piece: neither the annulus nor the disc
        t = '\n'.join(GP.va_lines([r for r in va if r['piece'] != k], V, D))
        assert "on the whole open disc" not in t and "So the valley condition holds" not in t and f"missing {k}" in t, k
    put = lambda k, **f: [dict(r, **f) if r['piece'] == k else r for r in va]
    for k, f in (('A-J', dict(boxes=0)), ('A-FC', dict(boxes=0)), ('O', dict(segs=0)), ('O', dict(worst=1e-3)), ('handoff', dict(sub_boxes=0)),
                 ('handoff', dict(points=0)), ('handoff', dict(boxes=75)), ('end', dict(bound=0.5)), ('sound B-F', dict(points=0)),
                 ('sound A-FC', dict(worst=1e-3)), ('meta', dict(eps0=0.004)), ('validate', dict(counts=dict(points=0))),
                 ('validate', dict(A_vs_def=1e-19)), ('validate', dict(B_vs_centred=1.0)), ('validate', dict(A0=1e-30))):
        assert not whole(put(k, **f), V, D) and not GP.va_complete(put(k, **f))[0], (k, f)
        assert "NOT complete" in '\n'.join(GP.va_lines(put(k, **f), V, D)), (k, f)
    assert whole(put('validate', A_vs_def=9e-21, B_vs_centred=9e-21), V, D)       # just below ε₀ = 0.005's bar: passes
    assert not whole(put('meta', eps0=0.11), V, D)                       # meta's ε₀ against R1 = 1.995
    assert not whole(put('meta', eps0=0.005 + 1e-9), V, D)               # … exactly: 10⁻⁹ off is refused
    md = open(os.path.join(GP.OUT, 'gh_peano.md'), encoding='utf-8').read()
    assert "the limit (iii)(b) holds on the whole open disc |η| < 2" in md and "So the valley condition holds on 1.995 ≤ |η| < 2" in md
    assert "(numerically: floating point at 20001 angles and the ends, the printed digits rounded outward; referee round 131)" in md


def test_round131_the_standalone_cli_and_the_tails(tmp_path):
    """Referee round 131 (D3dz-4): `gh_lastann.py all` runs to the end (it crashed at the outer piece's JSON log of the complex joints)
    and writes every piece; the outer piece's log is JSON without the joints, which stay in memory for `handoff_check`; an unknown
    piece is an error. D3dz-1: a series tail that underflows in doubles is floored at TAIL_FLOOR, so the ball still contains the
    sum (a series whose kept terms vanish exactly, its remainder ≈ 10⁻³⁶⁰)."""
    import json
    import subprocess
    import gh_lastann as LA
    from flint import arb, ctx
    LA.ctx.prec = LA.PREC
    logs = []
    O = LA.piece_outer(log=logs.append)
    assert len(logs) == 1 and json.loads(logs[0])['boxes'] == 76 and 'joints' not in json.loads(logs[0]) and len(O['joints']) == 76
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    env = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    out = tmp_path / 'va.json'
    p = subprocess.run([sys.executable, os.path.join('experiments', 'gh_lastann.py'), 'all', '--json', str(out)], cwd=repo, env=env,
                       capture_output=True, text=True, timeout=600)
    assert p.returncode == 0, p.stderr[-3000:]
    js = json.load(open(out))
    by = {r['piece']: r for r in js['results']}
    assert set(by) == set(GP.VA_PIECES) - {'meta'} and js['eps0'] == 0.005 and js['prec'] == LA.PREC
    assert all(by[k]['failed'] == 0 for k in by if 'failed' in by[k]) and by['handoff']['tiled'] and by['end']['ok']
    raw = GP.E.load(os.path.join(GP.OUT, 'gh_peano_raw.npz'))
    rb = {r['piece']: r for r in raw['VA']}
    assert all(by[k]['boxes'] == rb[k]['boxes'] for k in ('A-J', 'A-FC', 'B-J', 'B-F', 'O')) and by['O']['segs'] == rb['O']['segs']
    assert by['validate']['B_vs_centred'] == rb['validate']['B_vs_centred'] and by['handoff']['sub_boxes'] == rb['handoff']['sub_boxes']
    q = subprocess.run([sys.executable, os.path.join('experiments', 'gh_lastann.py'), 'A'], cwd=repo, env=env, capture_output=True, text=True)
    assert q.returncode != 0 and 'unknown piece' in q.stderr
    x = arb(1e-200)                                                      # the kept terms 0 exactly, the remainder x³/(1 − x)
    tl0, tl1 = (lambda r: r ** 3 / (1 - r)), (lambda r: 3 * r ** 2 / (1 - r) ** 2)
    assert tl0(1e-200) == 0.0 and tl1(1e-200) == 0.0 and LA.TAIL_FLOOR == 1e-300   # both double tails underflow
    y = LA._series(x, lambda k: arb(0), 3, tl0, tl1)
    p0 = ctx.prec
    ctx.prec = 600
    true, dtrue = x ** 3 / (1 - x), (3 * x ** 2 - 2 * x ** 3) / (1 - x) ** 2
    ctx.prec = p0
    assert y.contains(true) and y.rad() >= LA.TAIL_FLOOR
    yd = LA._series(LA.D(x, 1), lambda k: arb(0), 3, tl0, tl1)
    assert yd.v.contains(true) and yd.d.contains(dtrue) and yd.d.rad() >= LA.TAIL_FLOOR
    assert LA._series(arb(0.25), lambda k: arb(1), 40, lambda r: r ** 40 / (1 - r), lambda r: 40 * r ** 39 / (1 - r) ** 2).rad() < 1e-20
    nan = LA._series(arb(0.1), lambda k: arb(1), 3, lambda r: float('nan'), lambda r: float('nan'))
    assert not (nan < 10) and not (nan > -10)                            # a NaN tail is not floored away: the ball decides nothing


def test_round131_the_s_form_variance_exactly():
    """Referee round 131 (D3dz-7): the s-form S =_d sT + (1 − s)Z, s ~ Beta(K, K + 1) independent of T and Z, gives Var S = E s²·Var T +
    E(1 − s)² = (K² + K + 1)/((K + 1)(2K + 1)) in one line — exactly in rationals for K ≤ 12, with Var T = Σ_k z_k²/(K(K + 1)) from He_K's
    coefficients (Σ z_k² = K(K − 1)), against (ii)'s E S² = 2(2K)!/(2K + 2)!·P_1(K) from the Dyck-path defects."""
    from math import comb
    for K in range(1, 13):
        Es2 = Fraction(K * (K + 1), (2 * K + 1) * (2 * K + 2))           # Beta(a, b): E s² = a(a + 1)/((a + b)(a + b + 1))
        E1s2 = Fraction((K + 1) * (K + 2), (2 * K + 1) * (2 * K + 2))
        e2 = -comb(K, 2)                                                 # He_K = x^K − C(K, 2)x^{K−2} + …: the second elementary symmetric
        sumsq = 0 - 2 * e2                                               # (Σz)² − 2e₂ with Σz = 0
        assert sumsq == K * (K - 1)
        varT = Fraction(sumsq, K * (K + 1))                              # E λ_kλ_l = (1 + δ_kl)/(K(K + 1)) for Dirichlet(1, …, 1)
        assert varT == Fraction(K - 1, K + 1)
        var = Es2 * varT + E1s2
        p1 = Fraction(GP.defect(2 * K + 2, K), GP.defect(2 * K, K))
        assert var == Fraction(K * K + K + 1, (K + 1) * (2 * K + 1)) == Fraction(2 * factorial(2 * K), factorial(2 * K + 2)) * p1, K


# ---------------------------------------------------------------------------------------------------- round 135 (I4): the stored certificates
CERT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'certificates'))


def _verifier():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tools')))
    import verify_certificates as V
    return V


def test_round135_the_verifier_passes_on_the_stored_certificates():
    """Referee round 135 (I4, review M5): `tools/verify_certificates.py` checks Theorem 3's certificates (D on 1 < |η| ≤ 1.8927, VA at
    ε₀ = 0.11) from results/certificates alone — every lemma's line [ok], the counts of the md (402664 leaves; 54, 2533, 413, 1189 inner
    boxes; 220 η-boxes, 16083 segments, 1306 sub-boxes) and its largest bounds (D −1.9432·10⁻¹¹, A-FC −8.781·10⁻¹⁴, outer −4.228·10⁻⁹, end
    −1.3314); it imports nothing from experiments/ (a fresh interpreter with only tools/ on the path)."""
    import subprocess
    V = _verifier()
    out = []
    assert V.run(CERT, log=out.append) == 0, '\n'.join(out)
    txt = '\n'.join(out)
    assert all(ln.startswith(('[ok]', '[info]')) for ln in out[:-1]) and out[-1].startswith('VERIFIED:'), txt
    assert [ln for ln in out if ln.startswith('[info]')] and all('not a check' in ln for ln in out if ln.startswith('[info]'))  # round 136
    for s in ('402664 leaves, 0 failed, largest bound -1.9432e-11', 'tiled', 'A-FC: 2533 boxes, 0 failed', 'largest bound -8.781e-14',
              'A-J: 54 boxes', 'B-J: 413 boxes', 'B-F: 1189 boxes', '16083 segments on 220 η-boxes, 0 failed, largest bound -4.228e-09',
              '1306 sub-boxes', 'bound -1.3314', "max |Δ| against the proposer's bound 6e-17"):
        assert s in txt, s
    assert sum('Lemma A.8' in ln and 'covered' in ln for ln in out) == 4   # round 136: the paper's labels
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    code = ("import sys; sys.path.insert(0, 'tools'); import verify_certificates as V; "
            "assert not [m for m in sys.modules if m.startswith(('exp_', 'gh_'))]; sys.exit(V.main(['--only', 'VA', '--quiet']))")
    p = subprocess.run([sys.executable, '-I', '-c', code], cwd=repo, capture_output=True, text=True, timeout=600,
                       env=dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1'))
    assert p.returncode == 0, p.stderr[-2000:]


def test_round135_the_verifier_fails_on_planted_corruptions(tmp_path):
    """Round 135 (I4): each corruption is caught by its lemma — D: a leaf removed (the area), a leaf replaced by a copy of an equal one
    elsewhere (the same area, an overlap: the sweep), the worst leaf widened (its bound ≥ 0 or undefined); VA: the worst A-FC box
    widened in q (its bound ≥ 0), an inner box removed (the cover), the floor of a B-F box not provable (widened in p), a polyline vertex moved, a
    segment's end moved (the chain), the connector's end moved, a sub-box removed, an η-box removed (the tiling); the files: a hash
    mismatch, and a leaf removed with SHA256SUMS rewritten (the run fails at D's tiling, exit status 1)."""
    import copy
    import gzip
    import hashlib
    import shutil
    V = _verifier()
    meta, leaves, inner, outer = V.load(CERT)
    quiet = lambda s: None
    big = 10 ** 9                                                       # D2 on leaf 0 only: the tilings are what is tested here
    assert V.verify_D(meta, leaves, sample=big, log=quiet)
    assert not V.tiles((1.0, 1.8927, 0.0, V.PI2), [x[:4] for x in leaves[:100] + leaves[101:]])[0]
    shape = {}
    for i, x in enumerate(leaves):
        shape.setdefault((x[1] - x[0], x[3] - x[2]), []).append(i)
    i, j = next(v[:2] for v in shape.values() if len(v) > 1)
    dup = list(leaves)
    dup[i] = leaves[j]
    ok, why = V.tiles((1.0, 1.8927, 0.0, V.PI2), [x[:4] for x in dup])
    assert not ok and 'overlaps' in why
    w = max(range(len(leaves)), key=lambda k: leaves[k][4])
    R1, R2, f1, f2, H = leaves[w]
    from flint import ctx
    ctx.prec = V.D_PREC
    assert V.d_bound(R1, R2, f1, f2) < 0
    b = V.d_bound(R1, R2 + 4 * (R2 - R1), max(0.0, f1 - 4 * (f2 - f1)), f2)
    assert b is None or not (b < 0)
    ctx.prec = V.VA_PREC
    assert V.verify_inner(meta, inner, log=quiet) and V.verify_outer(meta, outer, log=quiet)
    afc = [k for k, r in enumerate(inner) if r['piece'] == 'A-FC']
    k = max(afc, key=lambda k: max([v for br, route, v in inner[k]['recs'] if route == 'v'] or [-math.inf]))
    bad = copy.deepcopy(inner)
    lo, hi = bad[k]['iv']['q']
    bad[k]['iv']['q'] = (lo, 2 * hi)                                   # its mean value bound becomes +3.0·10⁻¹⁴
    iv = bad[k]['iv']
    ok, route, v = V.mean_value(V.fA, V.fA0, 'C', [V.hull(*iv['u']), V.hull(*iv['q'])], [V.Iv(*iv['u']), V.Iv(*iv['q'])], *iv['r'])
    assert not ok and route == 'v' and v >= 0
    out = []
    assert not V.verify_inner(meta, bad, log=out.append) and any(ln.startswith('[FAIL] Lemmas A.4–A.7 A-FC: 2533 boxes, 1 failed') for ln in out)
    out = []
    assert not V.verify_inner(meta, inner[:afc[7]] + inner[afc[7] + 1:], log=out.append)
    assert any(ln.startswith('[FAIL] Lemma A.8 A-FC') and 'uncovered piece' in ln for ln in out), out
    bf = [k for k, r in enumerate(inner) if r['piece'] == 'B-F']
    k = max(bf, key=lambda k: inner[k]['iv']['p'][1])
    bad = copy.deepcopy(inner)
    bad[k]['iv']['p'] = (0.0, 4.0)
    out = []
    assert not V.verify_inner(meta, bad, log=out.append) and any('B-F: 1189 boxes, 1 failed' in ln for ln in out)

    def outer_fails(mutate, what):
        o = copy.deepcopy(outer)
        mutate(o)
        out = []
        assert not V.verify_outer(meta, o, log=out.append)
        assert any(ln.startswith('[FAIL] ' + what) for ln in out), out

    k = max(range(len(outer)), key=lambda k: len(outer[k]['pts']))
    m = len(outer[k]['pts']) // 2

    def move_vertex(o):
        o[k]['pts'][m] += 1e-9j

    def move_segment_end(o):
        P, Q = o[k]['poly'][5]
        o[k]['poly'][5] = (P, Q + 1e-9)

    def move_connector(o):
        P, Q = o[k]['conn'][-1]
        o[k]['conn'][-1] = (P, Q - 1e-9j)

    outer_fails(move_vertex, 'Lemma A.10')                             # round 136: the paper's labels (VA7, VA8, VA6 in the drafts)
    outer_fails(move_segment_end, 'Lemma A.10')
    outer_fails(move_connector, 'Lemma A.10')
    outer_fails(lambda o: o[k]['sub'].pop(), 'Lemma A.10')
    outer_fails(lambda o: o.pop(k), 'Lemma A.11')
    o = copy.deepcopy(outer)                                            # a polyline piece moved far: its bound is no longer negative
    P, Q = o[k]['poly'][3]
    o[k]['poly'][3] = (P, 1.5 + 0.5j)
    out = []
    assert not V.verify_outer(meta, o, log=out.append) and any(ln.startswith('[FAIL] Lemma A.9') for ln in out)

    d = tmp_path / 'cert'
    shutil.copytree(CERT, d)
    with open(d / 'meta.json', 'a') as fh:
        fh.write(' ')
    out = []
    assert V.run(str(d), log=out.append) == 1 and 'hash mismatch' in out[0] and len(out) == 2
    shutil.copy(os.path.join(CERT, 'meta.json'), d / 'meta.json')
    import array
    a = array.array('d', [v for x in leaves[1:] for v in x])
    with open(d / 'D_leaves.f64le.gz', 'wb') as fh, gzip.GzipFile(filename='', mode='wb', fileobj=fh, mtime=0) as g:
        g.write(a.tobytes())
    sums = ''.join(f"{hashlib.sha256(open(d / f, 'rb').read()).hexdigest()}  {f}\n" for f in V.FILES)
    (d / 'SHA256SUMS').write_text(sums)
    out = []
    assert V.run(str(d), only=('D',), d_sample=big, log=out.append) == 1
    assert out[0].startswith('[ok] SHA256SUMS') and any(ln.startswith('[FAIL] D meta') for ln in out), out


def test_round135_the_verifier_s_lemmas_against_the_programs_and_the_definition():
    """Round 135 (I4): the verifier's independent re-implementation against the proposing programs on random boxes (the same
    inequalities, so the same bounds to rounding: D2's bound vs `box_bound_ball`, VA6's vs `seg_bound`, the forms, the closed forms and
    the Y enclosure vs `fA`, `fB`, `fA0`, `fB0` and `_y_enclosed`) and the forms against E's definition (log|1 − z²| + Re(cz²) + C(σ))
    at path points; the export hook leaves `dom_bb`'s result unchanged."""
    import random
    import gh_lastann as LA
    from flint import acb, arb, ctx
    V = _verifier()
    rnd = random.Random(135)
    ctx.prec = V.D_PREC
    for _ in range(300):
        R1 = rnd.uniform(1.0, 1.85)
        R2, f1 = R1 + rnd.uniform(1e-6, 0.04), rnd.uniform(0, 1.5)
        f2 = min(f1 + rnd.uniform(1e-6, 0.05), math.pi / 2) if rnd.random() < 0.9 else math.pi / 2
        mine, theirs = V.d_bound(R1, R2, f1, f2), GP.box_bound_ball(R1, R2, f1, f2)
        assert (mine is None) == (theirs is None) and (mine is None or abs(float(mine) - theirs) < 1e-15)
    ctx.prec = LA.PREC
    for _ in range(200):
        u, q, lam = rnd.uniform(0, 4), rnd.uniform(0, 2.5), rnd.uniform(0, 1)
        rho = 10 ** rnd.uniform(-8, math.log10(min(0.05, 0.11 / max(u, 1e-9))))
        for br in ('F', 'C', 'J'):
            uu = 4.0 if br == 'J' else u
            a, b = V.fA(br, arb(uu), arb(q), arb(rho), arb(lam)), LA.fA(arb(uu), arb(q), arb(rho), br, arb(lam))
            assert abs(float(a.mid() - b.mid())) <= 1e-30 and a.rad() <= b.rad() * (1 + 1e-6) + 1e-60   # the same value, no looser
            a0, b0 = V.fA0(br, V.Iv(uu), V.Iv(q), V.Iv(lam)), LA.fA0(LA.I(uu), LA.I(q), br, LA.I(lam))
            assert abs(a0.hi - b0.hi) < 1e-30
        p, t = rnd.uniform(0, 0.5), rnd.uniform(1e-6, 0.1)
        for br in ('F', 'J'):
            uu = 4.0 if br == 'J' else u
            if uu * p * t > 0.11:
                continue
            a, b = V.fB(br, arb(uu), arb(p), arb(t), arb(lam)), LA.fB(arb(uu), arb(p), arb(t), br, arb(lam))
            assert abs(float(a.mid() - b.mid())) <= 1e-30 and a.rad() <= b.rad() * (1 + 1e-6) + 1e-60
            assert abs(V.fB0(br, V.Iv(uu), V.Iv(p), V.Iv(lam)).hi - LA.fB0(LA.I(uu), LA.I(p), br, LA.I(lam)).hi) < 1e-30
    for _ in range(100):                                                # the forms against E's definition at path points
        eps, th = 10 ** rnd.uniform(-6, math.log10(0.11)), rnd.uniform(0, 1.5)
        e, t = arb(eps), arb(th)
        om = 1 - e / 2
        eta = (2 - e) * acb(0, t).exp()
        A = om * om * acb(0, 2 * t).exp() - 1
        kap = A.imag / abs(A.real)
        X = e / 4 * rnd.uniform(1, 3)
        z = acb(X, kap * X if kap * X < t else t)
        br = 'F' if kap * X < t else 'C'
        c = eta * eta / 4
        sg = eta / 2 * (1 + z)
        r = abs(sg)
        C = (r * r - 1) / 2 - r.log() - V.pos(abs(sg.real) - 1) ** 2 if r > 1 else arb(0)
        E_ = abs(1 - z * z).log() + (c * z * z).real + C
        q = t * t / X
        if q < 2.5:
            try:
                f = V.fA(br, e / X, q, X) * X * X
            except V.Fail:
                continue
            assert abs(f - E_) <= 1e-10 * abs(E_) + 1e-40, (eps, th, br)
    for o in V.load(CERT)[3][:40]:                                      # VA6's bound and the Y enclosure on stored data
        for P, Q in o['poly'][:20]:
            assert abs(float(V.seg_bound(P, Q, o['box']).upper()) - LA.seg_bound(P, Q, *o['box'])) < 1e-15
        assert all(V.y_enclosed(s, o['lo'], o['hi'], 0.05) for s in o['sub'])
        assert LA._y_enclosed(*o['box'], o['lo'], o['hi'], 0.05) == len(o['sub'])
    keep = []
    assert GP.dom_bb(1.78, ball=True, keep=keep) == GP.dom_bb(1.78, ball=True) and len(keep) == 4935


def test_round135_the_export_reproduces_the_stored_certificates(tmp_path):
    """Round 135 (I4): `exp_gh_peano.py --export-cert DIR` recomputes D and VA011 (asserting equality with the raw's records) and writes
    the same bytes as results/certificates (gzip with mtime 0), with SHA256SUMS over the four data files; nothing else is written."""
    GP.export_cert(str(tmp_path))
    V = _verifier()
    assert sorted(os.listdir(tmp_path)) == sorted(V.FILES + ('SHA256SUMS',))
    for f in V.FILES + ('SHA256SUMS',):
        assert open(tmp_path / f, 'rb').read() == open(os.path.join(CERT, f), 'rb').read(), f
