"""Proposition D‴.62: the scaled Gauss–Hermite rules without a Mehler frame (`experiments/exp_frameless.py`)."""
import math
import os
import shutil
import sys
import warnings
from fractions import Fraction

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_frameless as FL

RAW = os.path.join(FL.OUT, 'frameless_raw.npz')
MD = os.path.join(FL.OUT, 'frameless.md')


@pytest.fixture(scope='module')
def raw():
    return FL.load(RAW)


@pytest.fixture(scope='module')
def F(raw):
    return FL.part_F(raw['S'], raw['A'])


def arb_(prec=600):
    from flint import arb, ctx
    ctx.prec = prec
    return arb


# ------------------------------------------------------------------------------------------------------------ the md and the raw

def test_the_md_is_the_report_of_the_raw(raw):
    """The md is `report` of the raw, byte for byte (a report-only change fails here, not only in `--from-raw`)."""
    md = open(MD, encoding='utf-8').read()
    assert FL.report(*[raw.get(k) for k in FL.KEYS]) == md
    assert FL.KEYS == ('R', 'S', 'A', 'T', 'G') and set(raw) == set(FL.KEYS)


def test_from_raw_regenerates_the_md_and_only_writes(monkeypatch, tmp_path):
    """`--from-raw` writes the md from a copy of the raw alone (the md missing), byte for byte, and leaves the raw unchanged; `main --only A`
    recomputes A alone into the raw (its other parts kept)."""
    shutil.copy(RAW, tmp_path / 'frameless_raw.npz')
    before = open(tmp_path / 'frameless_raw.npz', 'rb').read()
    monkeypatch.setattr(FL, 'OUT', str(tmp_path))
    monkeypatch.setattr(sys, 'argv', ['x', '--from-raw'])
    FL.main()
    assert open(tmp_path / 'frameless.md', encoding='utf-8').read() == open(MD, encoding='utf-8').read()
    assert open(tmp_path / 'frameless_raw.npz', 'rb').read() == before
    calls = []
    monkeypatch.setattr(FL, 'part_A', lambda: calls.append(1) or [])
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'A'])
    FL.main()
    got = FL.load(str(tmp_path / 'frameless_raw.npz'))
    assert calls == [1] and got['A'] == [] and set(got) == set(FL.KEYS)


def test_the_md_s_labels_and_sections():
    md = open(MD, encoding='utf-8').read()
    for h in ('## R. The representation (a) at small K', '## S. The exact errors at consecutive K **[numerical]**', '### F1. ',
              '### F2. The limit and s₁', '### F3. The second ridge', "### F4. A6's stored errors at K = 160", "## T. Lemma T's margin",
              "### F5. Beyond A6's rules", '## G. Assumption (F) on the axes **[numerical]**', 'ρ = −1 (b ≡ 0', 'Pinned for the tests'):
        assert h in md, h
    assert '**[proved]**' in md and '**[formal]**' in md and 'NO' not in md.replace('NOT', '')


# ----------------------------------------------------------------------------------------------------------------- closed forms

def test_moment_matching_constant_and_s1_are_the_rationals():
    """C = (1 + λ)^{3/2}/√2 = 2σ*³/(1 + σ*²)^{3/2} (C² rational), s₁ = −13/125, −11/64, −187/2000 and c₂ = s₁/2 at τ/σ = ½, 1, 2 exactly
    (`mm_s1`, `s1_formal` in Fractions at s = −λ, ϖ = 0, and the unsimplified expansion), and `C_const`, `s1_formal` in arb contain them."""
    arb = arb_()
    want = {Fraction(1, 4): Fraction(-13, 125), Fraction(1): Fraction(-11, 64), Fraction(4): Fraction(-187, 2000)}
    for t2, s1 in want.items():
        lam = FL.mm_lambda(t2)
        assert lam == t2 / (2 + t2) and FL.mm_s1(t2) == s1 == FL.s1_formal(-lam, Fraction(0)) == FL.s1_unsimplified(-lam, Fraction(0))
        assert FL.c2_formal(-lam, Fraction(0)) == s1 / 2 and FL.mm_C2(t2) == (1 + lam) ** 3 / 2
        ss2 = 1 + t2
        assert FL.mm_C2(t2) == Fraction(4) * ss2 ** 3 / (1 + ss2) ** 3          # C² = 4σ*⁶/(1 + σ*²)³
    assert FL.mm_s1(Fraction(1, 16)) == Fraction(8, 4913)
    for tau, t2 in (('0.5', Fraction(1, 4)), ('1', Fraction(1)), ('2', Fraction(4)), ('0.25', Fraction(1, 16))):
        P = FL.params(tau, FL.rule_v(tau, 'mm'))
        assert P['kind'] == 'MM' and P['dom'] == 'b' and float(P['varpi'].mid()) == 0
        C = FL.C_const(tau, P['v'], P)
        assert abs(C * C - FL.A_(FL.mm_C2(t2))) < arb(10) ** -150
        assert abs(C - (1 + FL.A_(FL.mm_lambda(t2))) ** FL.A_('1.5') / arb(2).sqrt()) < arb(10) ** -150
        assert abs(FL.s1_formal(P['s_d'], P['varpi']) - FL.A_(FL.mm_s1(t2))) < arb(10) ** -150
        assert FL.C_other(tau, P['v'], P) is None and FL.B0_formal(tau, P['v'], P) is None


def test_the_outliers_b0_and_the_second_ridge():
    """The outlier (τ/σ = ¼, 2·v_hs): B₀ = 0.2383158459 (10 digits), B₁/B₀ = −0.0050412152, q = ϖ² = 0.8506295388; C_o/C = 0.9411492563;
    and C_o is C with d and o exchanged."""
    arb_()
    P = FL.params('0.25', FL.rule_v('0.25', '2'))
    B0, B1, q = FL.B0_formal('0.25', P['v'], P)
    assert P['kind'] == 'rho>1' and FL.fs(B0, 10) == '0.2383158459' and FL.fs(B1, 8) == '-0.0050412152'
    assert FL.fs(q, 10) == '0.8506295388' and FL.fs(FL.C_other('0.25', P['v'], P) / FL.C_const('0.25', P['v'], P), 10) == '0.9411492563'
    Q = dict(P, s_d=P['s_o'], s_o=P['s_d'], varpi=1 / P['varpi'])
    assert abs(FL.C_const('0.25', P['v'], Q) / FL.C_other('0.25', P['v'], P) - 1) < FL.A_('1e-150')


def test_c_is_d61s_constant_at_the_frames_and_the_algebraic_identity():
    """(e): C(v) with ϖ < 0 equals D‴.61 (iv)'s constant at frames with the mass on either axis, to 10⁻²⁵⁰ at 1000 bits, with r₆₁ = t(v);
    the identity 4(τ² − v)(vσ*² − τ²) = τ⁴(1 − ρ²)v² and c·c_b = −(1 − ρ²)v²/16 hold exactly (Fractions) at random rational τ², v."""
    arb = arb_(1000)
    for tau, w, dom in (('2', '0.5', 'a'), ('2', '2', 'b'), ('1', 'Fock', 'a'), ('1', '1.25', 'b'), ('0.25', 'Mercer-zero', 'a')):
        P = FL.params(tau, FL.rule_v(tau, w))
        assert P['kind'] == 'frame' and P['dom'] == dom and float(P['varpi'].mid()) < 0
        r61, C61 = FL.d61_constant(tau, P['v'])
        assert abs(FL.C_const(tau, P['v'], P) / C61 - 1) < arb(10) ** -250 and abs(r61 - P['t']) < arb(10) ** -250
    rng = np.random.default_rng(62)
    for _ in range(50):
        t2, v = Fraction(int(rng.integers(1, 400)), 100), Fraction(int(rng.integers(1, 900)), 100)
        ss2 = 1 + t2
        rho = 1 + 2 / t2 - 2 / v
        c = (1 - v / t2) / 2
        cb = c - v / 2
        assert 4 * (t2 - v) * (v * ss2 - t2) == t2 ** 2 * (1 - rho ** 2) * v ** 2 and c * cb == -(1 - rho ** 2) * v ** 2 / 16
        assert 1 - rho == 2 * (t2 - v) / (v * t2) and 1 + rho == 2 * (v * ss2 - t2) / (v * t2)
        assert c / (1 - c) == (t2 - v) / (t2 + v) and cb / (1 - cb) == (t2 - v * ss2) / (t2 + v * ss2)


def test_the_dichotomy_the_rate_and_the_duality_on_a_grid():
    """|ρ| < 1 ⟺ s_a s_b < 0 (`params`' kind), t(v) = max(|s_a|, |s_b|) = tanh((r + |ln(v/v_hs)|)/2), the branches (ρ ≥ 1: s_d = s_b < 0,
    ρ ≤ −1: s_d = s_a > 0) and the duality v ↦ v_hs²/v: (s_a, s_b) ↦ (−s_b, −s_a), on a grid of τ/σ and v/v_hs."""
    arb_(200)
    n = 0
    for tau in ('0.1', '0.3', '0.7', '1.5', '4'):
        t2 = float(tau) ** 2
        ss = math.sqrt(1 + t2)
        vhs = t2 / ss
        for x in np.geomspace(0.01, 100, 41):
            v = FL.A_(repr(float(x))) * FL.A_(tau) ** 2 / (1 + FL.A_(tau) ** 2).sqrt()
            P = FL.params(tau, v)
            sa, sb = float(P['s_a'].mid()), float(P['s_b'].mid())
            frame = abs(float(P['rho'].mid())) < 1
            assert frame == (sa * sb < 0) and (P['kind'] in ('frame', 'cusp')) == frame
            assert abs(float(P['t'].mid()) - math.tanh((math.log(ss) + abs(math.log(x))) / 2)) < 1e-12
            if not frame:
                assert (P['kind'] == 'rho>1') == (sb < 0) and P['dom'] == ('b' if sb < 0 else 'a') and 0 <= float(P['varpi'].mid()) < 1
            Q = FL.params(tau, FL.A_(repr(vhs * vhs)) / v)
            assert abs(float(Q['s_a'].mid()) + sb) < 1e-9 and abs(float(Q['s_b'].mid()) + sa) < 1e-9
            n += 1
    assert n == 205


def test_s1_formal_is_round_129s_formula_and_the_unsimplified_expansion():
    """s₁(s, ϖ) = round 129's s₁ at ϱ = −s, √q = −ϖ and = (d)'s unsimplified expansion, exactly (Fractions) at 40 random (s, ϖ) on both
    sides of ϖ = 0; F₁'s definition; c₂ = s₁/2."""
    rng = np.random.default_rng(129)
    for _ in range(40):
        s = Fraction(int(rng.integers(-95, 95)), 100)
        w = Fraction(int(rng.integers(-95, 95)), 100)
        assert FL.s1_formal(s, w) == FL.s1_r129(-s, -w) == FL.s1_unsimplified(s, w) and FL.c2_formal(s, w) == FL.s1_formal(s, w) / 2
    y = Fraction(1, 3)
    assert FL.F1(y) == y * (1 + y) / (1 - y) ** 3 + y / (8 * (1 - y) ** 2) + y / (8 * (1 - y))


def test_stirlings_factor():
    """(K!/(2K)!)²(4K − 1)!! → 2^{−1/2}(1 + 1/(16K) + …), and (K!/(2K)!)K^K(4/e)^K → 2^{−1/2}."""
    for K in (200, 2000):
        lg = math.lgamma
        l1 = 2 * (lg(K + 1) - lg(2 * K + 1)) + lg(4 * K + 1) - lg(2 * K + 1) - 2 * K * math.log(2)
        assert abs(math.exp(l1) * math.sqrt(2) - 1 - 1 / (16 * K)) < 1 / K ** 2
        l2 = lg(K + 1) - lg(2 * K + 1) + K * math.log(K) + K * math.log(4 / math.e)
        assert abs(math.exp(l2) * math.sqrt(2) - 1) < 1 / K


def test_lemma_a_on_a_grid():
    """Lemma A: Φ = 2 + 2 log(|a² − b²|/4) − a²/(2|s_a|) − b²/(2|s_b|) ≤ 2 log t, attained only at (0, ±2|s_b|^{1/2}) or (±2|s_a|^{1/2}, 0),
    and per sector at its own point with 2 log|s|; on a grid (`phi_rate`)."""
    for sa, sb in ((-0.1, -0.35), (0.38, 0.28), (0.0303, 0.01), (-0.5, -0.9)):
        g = np.linspace(-3, 3, 1201)
        a, b = np.meshgrid(g, g)
        F = FL.phi_rate(a, b, sa, sb)
        t = max(abs(sa), abs(sb))
        assert np.nanmax(F) <= 2 * math.log(t) + 1e-12
        for sec, s, pt in ((np.abs(a) >= np.abs(b), sa, (2 * math.sqrt(abs(sa)), 0.0)), (np.abs(b) >= np.abs(a), sb, (0.0, 2 * math.sqrt(abs(sb))))):
            assert np.nanmax(F[sec]) <= 2 * math.log(abs(s)) + 1e-12
            assert abs(FL.phi_rate(np.array(pt[0]), np.array(pt[1]), sa, sb) - 2 * math.log(abs(s))) < 1e-12
            i = np.unravel_index(np.nanargmax(np.where(sec, F, -np.inf)), F.shape)
            assert abs(abs(a[i]) - pt[0]) < 0.02 and abs(abs(b[i]) - pt[1]) < 0.02


def test_lemma_t_reduction_and_rate_l():
    """L(y) = 1 + log(y²/4) − y²/4 ≤ 0, zero only at |y| = 2, increasing on (0, 2]; `t_reduced`'s R1 bound is the least of −Q on the
    corner and below 2 log|s| (log(1 + 1/k) < 1/k); R2's edge value bounds Φ⁺ on a grid of {|y₁| ≥ 2 > |y₂|} in its sector."""
    y = np.linspace(0.01, 2, 500)
    L = FL.L_rate(y)
    assert np.all(L <= 1e-15) and abs(FL.L_rate(np.array(2.0))) < 1e-15 and np.all(np.diff(L) > 0) and FL.L_rate(np.array(3.0)) < 0
    for sa, sb in ((-0.096, -0.126), (0.347, 0.320), (0.17, 0.48)):
        red = FL.t_reduced(sa, sb, 20001)
        ka, kb = abs(sa) / (1 - abs(sa)), abs(sb) / (1 - abs(sb))
        assert red['a'][0] == -2 / ka and red['b'][0] == -2 / kb and -2 / ka < 2 * math.log(abs(sa)) and -2 / kb < 2 * math.log(abs(sb))
        Q = lambda y1, y2: (y1 + y2) ** 2 / (8 * ka) + (y1 - y2) ** 2 / (8 * kb)
        y1, y2 = np.meshgrid(np.linspace(2, 6, 401), np.linspace(-1.999, 1.999, 801))
        Fp = FL.L_rate(y2) - Q(y1, y2)
        assert Fp[y2 > 0].max() <= red['a'][1] + 1e-9 and Fp[y2 < 0].max() <= red['b'][1] + 1e-9
        assert Q(2.0, 2.0) == 2 / ka and Q(2.0, -2.0) == 2 / kb


# ------------------------------------------------------------------------------------------------------------- live recomputation

def test_the_representation_live_at_small_k(monkeypatch):
    """(a) live: moment matching at τ/σ = 1, K = 3 by `r_mm` (acb.integral, 128 bits) and ρ > 1 / ρ < −1 at K = 3 by `r_2d` with a 40 × 40
    and 50 × 50 product rule, against the closed form."""
    monkeypatch.setattr(FL, 'R_BITS', 128)
    r = FL.r_mm('1', 3)
    assert r['overlap'] and r['rel'] < 1e-30
    monkeypatch.setattr(FL, 'R_M', (40, 50))
    for tau, x in (('0.5', '2'), ('0.5', '0.5')):
        r = FL.r_2d(tau, x, 3)
        assert r['kind'] == ('rho>1' if x == '2' else 'rho<-1') and r['rel'] < 1e-12 and r['dM'] < 1e-12, r


def _nosecs(d):
    return {k: v for k, v in d.items() if k != 'secs'}


def test_r_recomputed_where_cheap(raw, monkeypatch):
    """R's moment-matching row at τ/σ = ¼, K = 5 by `r_mm`, every two-dimensional row by `r_2d` (M = 120, 160), the ρ = −1 row by
    `r_rm1` and, by `part_R` (R_MM and R_2D emptied), the ρ = −1 row and the pinned small product rules (R_PIN) recomputed equal the
    stored ones (but their time); the pinned values also numerically to 10⁻⁴⁰."""
    R = raw['R']
    got = FL.r_mm('0.25', 5)
    want = [r for r in R['mm'] if r['tau'] == '0.25'][0]
    assert _nosecs(got) == _nosecs(want)
    for want in R['two']:
        assert _nosecs(FL.r_2d(want['tau'], want['x'], want['K'])) == _nosecs(want), (want['tau'], want['x'])
    assert [(r['tau'], r['x'], r['K'], tuple(r['M'])) for r in R['pin']] == list(FL.R_PIN)
    monkeypatch.setattr(FL, 'R_MM', ())
    monkeypatch.setattr(FL, 'R_2D', ())
    got = FL.part_R()
    assert got['mm'] == [] and got['two'] == [] and [_nosecs(r) for r in got['rm1']] == [_nosecs(r) for r in R['rm1']]
    arb = arb_(400)
    for g, want in zip(got['pin'], R['pin']):
        for a, b in zip(g['vals'] + [g['ex']], want['vals'] + [want['ex']]):
            assert abs(arb(a) / arb(b) - 1) < arb(10) ** -40
        assert _nosecs(g) == _nosecs(want)
    assert len(got['pin']) == len(R['pin']) == 3


def test_s_rows_recomputed_and_the_precision(raw):
    """S's stored rows: `e_row` recomputes e_K at K = 2, 3, 20 of two rules exactly (the 60-digit strings, bits, radius); `s_rule` with K ≤ 4
    gives the stored first rows; every row's bits are `s_bits`'; every rule's K run S_K0 … K_max consecutively."""
    S = {(r['tau'], r['w']): r for r in raw['S']}
    assert [(t, w, k) for t, w, k in FL.S_RULES] == [(r['tau'], r['w'], r['kmax']) for r in raw['S']] and len(S) == 17
    assert FL.S_RULES == FL.S_A6 + FL.S_NEW and len(FL.S_A6) == 14 and sum(FL.s_is_a6(r) for r in raw['S']) == 14
    for key in (('1', '0.5'), ('0.25', '2')):
        rows = {r[0]: r for r in S[key]['rows']}
        for K in (2, 3, 20):
            assert FL.e_row(*key, K)[:4] == tuple(rows[K][:4])
        got = FL.s_rule(*key, 4)
        assert [tuple(r[:4]) for r in got['rows']] == [tuple(rows[K][:4]) for K in (2, 3, 4)]
        assert got['check'][2] == FL.s_bits(*key, 4, FL.S_RECHECK) and got['check'][1][:50] == rows[4][1][:50]
    for rec in raw['S']:
        assert [r[0] for r in rec['rows']] == list(range(FL.S_K0, rec['kmax'] + 1))
        assert all(r[2] == FL.s_bits(rec['tau'], rec['w'], r[0]) for r in rec['rows'])
        assert max(r[3] for r in rec['rows']) < 10.0 ** -(FL.S_DIGITS + 10) and rec['check'][0] == rec['kmax']
        assert rec['check'][1] == rec['rows'][-1][1]
        assert rec['check'][2] == FL.s_bits(rec['tau'], rec['w'], rec['kmax'], FL.S_RECHECK) == rec['rows'][-1][2] + FL.S_RECHECK
        assert all(len(r[1].replace('-', '').replace('.', '').split('e')[0]) >= FL.S_DIGITS - 1 for r in rec['rows'][::50])


def test_a_is_a6s_record(raw):
    """A is width_mismatch_raw.npz's A6 (read-only), string for string; `part_A` recomputes it; `a6_w`/`a6_tau` map its names."""
    if not os.path.exists(FL.A6_RAW):
        pytest.skip('no width_mismatch_raw.npz')
    assert FL.part_A() == raw['A'] and len(raw['A']) == 36
    assert FL.a6_w('moment matching') == 'mm' and FL.a6_w('1.25·v_hs') == '1.25' and FL.a6_w('Fock') == 'Fock' and FL.a6_tau(1.0) == '1'
    arb_(400)
    for x in raw['A']:
        v = FL.rule_v(FL.a6_tau(x['r']), FL.a6_w(x['name']))
        assert abs(float(v.mid()) - x['v']) < 1e-15 * x['v'] and [r['K'] for r in x['rows']] == list(FL.A6_K)


def test_err2_against_a6_at_k20():
    """`err2` at K = 20 against A6's 26-digit e_K of two rules (independent decimal code), in arb."""
    if not os.path.exists(FL.A6_RAW):
        pytest.skip('no width_mismatch_raw.npz')
    arb = arb_()
    A = FL.part_A()
    for tau, w in (('0.5', '0.8'), ('2', 'mm')):
        x = [z for z in A if FL.a6_tau(z['r']) == tau and FL.a6_w(z['name']) == w][0]
        e = FL.err2(tau, w, 20, FL.s_bits(tau, w, 20)).sqrt()
        arb_()
        assert abs(e / arb(x['rows'][0]['e1']) - 1) < arb(10) ** -24


def test_t_recomputed(raw):
    """T's rule records recomputed by `t_rule` at A6's first two rules (ρ < −1) and at one ρ > 1 rule, the sweep's (τ/σ, v/v_hs) by
    `t_sweep_rules`, its (s_a, s_b, t) by `t_sv` from (τ, w) at a sample, one sweep rule's grid (the smallest margin), and the reduction
    at every seventh sweep rule, equal the stored ones."""
    T = raw['T']
    for r in T['rules'][:2] + [T['rules'][6]]:
        assert FL.t_rule(r['tau'], r['w'], *FL.T_GRID) == r
    assert T['rules'][6]['kind'] == 'rho>1' and T['rules'][0]['kind'] == 'rho<-1'
    assert FL.t_sweep_rules() == [(r['tau'], r['w']) for r in T['sweep']]
    for r in T['sweep'][3::11]:
        c, cb, sa, sb = FL.t_sv(r['tau'], r['w'])
        assert (sa, sb, max(abs(sa), abs(sb)), 'rho>1' if c < 0 else 'rho<-1') == (r['sa'], r['sb'], r['t'], r['kind'])
        assert r['two_log_t'] == 2 * math.log(r['t'])
    i = int(np.argmin([FL.t_margins(r)[0] for r in T['sweep']]))
    r = T['sweep'][i]
    assert FL.t_rule(r['tau'], r['w'], *FL.T_SWEEP) == r
    for r in T['sweep'][::7]:
        assert FL.t_reduced(r['sa'], r['sb']) == r['red']
    assert len(T['rules']) == len(FL.T_RULES) == 10 and len(T['sweep']) == FL.T_TAUS[2] * 2 * FL.T_WN
    assert [(r['tau'], r['w']) for r in T['rules']] == list(FL.T_RULES)


def test_part_s_resumes(monkeypatch, tmp_path, raw):
    """`part_S` keeps the rules already in the raw (no `s_rule` call), recomputes a missing one, and saves the raw after each rule."""
    calls = []
    monkeypatch.setattr(FL, 's_rule', lambda tau, w, kmax: calls.append((tau, w, kmax)) or dict(tau=tau, w=w, kmax=kmax, rows=[], secs=0))
    res = dict(raw)
    out = FL.part_S(dict(raw), str(tmp_path / 'r.npz'), res)
    assert calls == [] and [r['kmax'] for r in out] == [k for _, _, k in FL.S_RULES]
    old = dict(raw, S=[r for r in raw['S'] if (r['tau'], r['w']) != ('0.5', '0.8')])
    out = FL.part_S(old, str(tmp_path / 'r.npz'), dict(old))
    assert calls == [('0.5', '0.8', 300)] and len(out) == 17 and len(FL.load(str(tmp_path / 'r.npz'))['S']) == 17


# ----------------------------------------------------------------------------------------------------------- semantics of the raw

def test_r_agrees(raw):
    R = raw['R']
    assert [(r['tau'], r['K']) for r in R['mm']] == list(FL.R_MM) and [(r['tau'], r['x'], r['K']) for r in R['two']] == list(FL.R_2D)
    assert all(r['overlap'] and r['rel'] < 1e-75 and float(r['tail']) < 1e-100 for r in R['mm'])
    assert all(r['rel'] < 1e-40 and r['dM'] < 1e-30 for r in R['two'])
    assert sorted(r['kind'] for r in R['two']) == ['rho<-1'] * 3 + ['rho>1'] * 3
    assert set(R) == {'mm', 'two', 'rm1', 'pin'} and [(r['tau'], r['K']) for r in R['rm1']] == list(FL.R_RM1)


def test_the_limits_and_s1(F):
    """The limit-free fits' best |L − 1| (moment matching and ρ > 1 below 10⁻²⁰, ρ < −1 below 10⁻⁸) and s₁ against the formal value
    (below 10⁻¹⁸, resp. 10⁻⁶)."""
    for r in F['rules']:
        best = min(float(f['L1'].mid()) for f in r['fits'])
        ds = abs(float((r['fits'][-1]['s1'] / r['s1'] - 1).mid()))
        tol = (1e-20, 1e-18) if r['kind'] in ('MM', 'rho>1') else (1e-8, 1e-6)
        assert best < tol[0] and ds < tol[1], (r['tau'], r['w'], best, ds)
        assert max(float(f['L1'].mid()) for f in r['fits']) < 1e-6
    assert sorted(r['kind'] for r in F['rules']) == ['MM'] * 4 + ['rho<-1'] * 5 + ['rho>1'] * 5


def test_the_second_ridge(F):
    """The outlier's second ridge resolved: g₀ = C_o/C within 10⁻¹¹ and g₁/g₀ = s₁′ within 10⁻⁶ (one window); five rules resolved."""
    rid = {(r['tau'], r['w']): r for r in F['rules'] if 'ridge' in r}
    assert len(rid) == 10
    out = rid[('0.25', '2')]
    assert out['verdict'] == 'resolved' and out['best'] < 1e-11
    assert min(abs(float((x['g1g0'] - out['s1o']).mid())) for x in out['ridge']) < 1e-6
    assert min(abs(float((x['b1'] - out['B1']).mid())) for x in out['ridge']) < 1e-6
    assert sum(r['verdict'] == 'resolved' for r in rid.values()) >= 5
    for r in rid.values():
        if r['verdict'] == 'resolved':
            assert min(abs(float((x['g1g0'] / r['s1o'] - 1).mid())) for x in r['ridge']) < 1e-2


def test_a6_and_the_overlap(F):
    """S against A6 within A6's 26-digit rounding; A6 at K = 160: the frame identity at 1000 bits, the cusp excluded, and moment matching
    and ρ > 1 with 160(R − 1) within 1 % of s₁ (10⁻³ where s₁ is small)."""
    w, n = F['overlap']
    assert w < 6e-26 and n == 14 * 12
    a6 = F['a6']
    assert len(a6) == 36 and sum(r['kind'] == 'cusp' for r in a6) == 4 and sum(r['kind'] == 'frame' for r in a6) == 18
    for r in a6:
        if r['kind'] == 'frame':
            assert float(r['ident']) < 1e-250 and float(r['rdiff']) < 1e-250
        if r['kind'] in ('MM', 'rho>1'):                                  # within O(s₂/K): 1 % of s₁, or 10⁻³ where s₁ is small
            assert abs(float((r['KR1'] - r['s1']).mid())) < max(1e-2 * abs(float(r['s1'].mid())), 1e-3)


def test_t_margins(raw):
    """Lemma T's margins positive at A6's rules (by the grid and the reduction, sectors too), the reduction's sup at least the grid's (it is
    the supremum), the grid's maximum at most 2 log t(v) near |y| = 2t^{1/2}; every sweep margin positive."""
    for r in raw['T']['rules'] + raw['T']['sweep']:
        m = FL.t_margins(r)
        assert min(m) > 0 and m[0] <= m[1] + 1e-9, (r['tau'], r['w'], m)
        assert r['grid']['gmax'] <= r['two_log_t'] + 1e-9
    for r in raw['T']['rules']:
        assert abs(r['grid']['gy'] - 2 * math.sqrt(r['t'])) < 0.02 and r['grid']['gmax'] > r['two_log_t'] - 1e-3


# ----------------------------------------------------------------------------------------------------------- round 132's additions

def test_rho_minus_one(raw, F):
    """ρ = −1 (τ/σ = 1, v = τ²/σ*² = ½; moment matching's dual): `params`' kind, s_d = s_a = λ = ⅓, ϖ = 0, no second ridge;
    C² = (1 − λ)²/(2 + τ²) = 4/27 and s₁ = 13/4 exactly (`rm1_C2`, `rm1_s1`, `s1_formal`) and in arb; R's rigorous one-dimensional
    check; S's fits: |L − 1| below 10⁻⁸ in the best window and the fitted s₁ within 10⁻⁶ of 13/4."""
    arb = arb_()
    P = FL.params('1', FL.rule_v('1', 'rho=-1'))
    assert P['kind'] == 'rho=-1' and P['dom'] == 'a' and abs(P['varpi']) < arb(10) ** -150 and FL.C_other('1', P['v'], P) is None
    assert abs(P['s_a'] - arb(1) / 3) < arb(10) ** -150 and abs(P['rho'] + 1) < arb(10) ** -150
    for d, kind in (('1e-14', 'rho=-1'), ('-1e-10', 'rho<-1'), ('1e-10', 'frame')):  # `params`' tolerance 10⁻¹² in ρ, as at ρ = 1
        assert FL.params('1', FL.A_('0.5') * (1 + FL.A_(d)))['kind'] == kind, d
    assert FL.rm1_C2(Fraction(1)) == Fraction(4, 27) and FL.rm1_s1(Fraction(1)) == Fraction(13, 4)
    assert FL.s1_formal(Fraction(1, 3), Fraction(0)) == Fraction(13, 4)
    C = FL.C_const('1', P['v'], P)
    assert abs(C * C - arb(4) / 27) < arb(10) ** -150 and abs(FL.s1_formal(P['s_d'], P['varpi']) - arb(13) / 4) < arb(10) ** -150
    for t2 in (Fraction(1, 4), Fraction(4), Fraction(9, 7)):             # the dual: C² = (1 − λ)²/(1 + σ*²), s₁ at s = +λ
        tau = repr(float(t2) ** 0.5)
        Q = FL.params(tau, FL.A_(t2) / (1 + FL.A_(t2)))
        assert Q['kind'] == 'rho=-1' and abs(FL.C_const(tau, Q['v'], Q) ** 2 / FL.A_(FL.rm1_C2(t2)) - 1) < 1e-12
    r = [x for x in raw['R']['rm1'] if x['tau'] == '1'][0]
    assert r['overlap'] and r['rel'] < 1e-75 and float(r['tail']) < 1e-100
    n = [x for x in F['new'] if x['kind'] == 'rho=-1']
    assert len(n) == 1 and n[0]['tau2'] == 1
    assert min(float(f['L1'].mid()) for f in n[0]['fits']) < 1e-8 and abs(float((n[0]['fits'][-1]['s1'] - FL.A_('3.25')).mid())) < 1e-6


def test_kuo_wozniakowski_threshold_points(F):
    """Moment matching at τ²/(2σ²) = 2αγ² = 1 and 4 (τ/σ = √2, 2√2): λ = 2αγ²/(1 + 2αγ²) = ½, ⅘ exactly, s₁ = −61/432 and −631/11664;
    y_K − 1 → 0 like s₁/K (K_max(y − 1) within 10⁻³ of s₁), the fits' |L − 1| below 10⁻²⁰ and s₁ within 10⁻¹⁵; the cheap bound
    (rate τ²/2 = 2αγ²) does not decay (≥ ½ at K_max) while e²_{K_max} < 10⁻⁵⁰."""
    arb = arb_()
    mm = {x['tau']: x for x in F['new'] if x['kind'] == 'MM'}
    assert set(mm) == {'√2', '2√2'} and mm['√2']['tau2'] == 2 and mm['2√2']['tau2'] == 8
    assert FL.mm_lambda(Fraction(2)) == Fraction(1, 2) and FL.mm_lambda(Fraction(8)) == Fraction(4, 5)
    assert FL.mm_s1(Fraction(2)) == Fraction(-61, 432) and FL.mm_s1(Fraction(8)) == Fraction(-631, 11664)
    for tau, r in mm.items():
        lam = FL.mm_lambda(r['tau2'])
        assert abs(r['t'] - FL.A_(lam)) < arb(10) ** -100 and abs(r['s1'] - FL.A_(FL.mm_s1(r['tau2']))) < arb(10) ** -100
        K, y = r['y1'][-1]
        assert K == r['kmax'] and abs(float((K * y - r['s1']).mid())) < 1e-3
        assert min(float(f['L1'].mid()) for f in r['fits']) < 1e-20 and abs(float((r['fits'][-1]['s1'] / r['s1'] - 1).mid())) < 1e-15
        assert float(r['cheap'].mid()) > 0.5 and float(r['e2max'].mid()) < 1e-50


def test_g_assumption_f_on_the_axes(raw, monkeypatch):
    """G: `part_G` at K = 60 and 120 (all six η) gives the stored strings; extrapolated through all the points D_K is F₁(η²/4) within
    10⁻⁶ on the imaginary axis and within 10⁻⁴ at η = 1 on the real axis."""
    G = raw['G']
    assert G['K'] == list(FL.G_K) and [(g['axis'], g['eta']) for g in G['rows']] == list(FL.G_ETA)
    monkeypatch.setattr(FL, 'G_K', (FL.G_K[0], FL.G_K[2]))
    got = FL.part_G()
    assert got['K'] == [G['K'][0], G['K'][2]] and [dict(g, D=[g['D'][0], g['D'][2]]) for g in G['rows']] == got['rows']
    from flint import arb, ctx
    ctx.prec = FL.F_PREC
    for g in G['rows']:
        e = FL.A_(g['eta'])
        y = e * e / 4 * (1 if g['axis'] == 'real' else -1)
        x6 = FL.g_extrap(G['K'], [arb(x) for x in g['D']])
        if g['axis'] == 'imag':
            assert abs(float((x6 - FL.F1(y)).mid())) < 1e-6
        elif g['eta'] == '1':
            assert abs(float((x6 - FL.F1(y)).mid())) < 1e-4
