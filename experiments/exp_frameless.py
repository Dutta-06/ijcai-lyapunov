"""
Proposition D‴.62: the scaled Gauss–Hermite rules without a Mehler frame — a Bochner representation on one axis, the rate and the sharp
constant
====================================================================================================================================

D‴.58's scaled rule: kernels of width σ = 1, a feature of width σ*, τ² = σ*² − 1, the K-point Gauss–Hermite rule of N(0, v) for the
mixing law reweighted to N(0, τ²); e_K = ‖f* − f_K‖/‖f*‖ (the closed form of referee round 129's arberr.py, as D‴.58's A6). With
c = (1 − v/τ²)/2, c_b = c − v/2, s_a = c/(1 − c) = h(τ²/v), s_b = c_b/(1 − c_b) = −h(vσ*²/τ²), h(y) = (y − 1)/(y + 1): D‴.61 (v)'s ratio
has |ρ| < 1 iff s_a·s_b < 0 (a frame); the rules without a frame are ρ ≥ 1 (v ≥ τ², moment matching ρ = 1) and ρ ≤ −1 (v ≤ τ²/σ*²).
s_d is the one of s_a, s_b of larger modulus, s_o the other, ϖ = s_o/s_d, t(v) = |s_d| (D‴.58's conjectured rate). The research note
(research note, 2026-11-02: review/RESEARCH_note_frameless.md) proves e_K² = σ*(v/τ²)E[φ_K(a + b)φ_K(a − b)], a ~ N(0, c), b ~ N(0, c_b)
(an imaginary Gaussian for a negative variance), φ_K(ξ) = (E − Q_K)[e^{ξz}], and from it e_K² ~ C(v)t(v)^{2K} for every rule without a
frame, C(v) = √2σ*v(1 − s_d)/[(τ² + v)(τ² + vσ*²)(1 + ϖ)]^{1/2}. This script:

R. The representation at small K against the exact closed-form error: at moment matching one-dimensional, e_K² = σ*(K!/(2K)!)²
   E[B^{4K}ĝ_K(iB)²] = σ*E[φ_K(iB)²], B ~ N(0, τ²/2), by a rigorous integral (python-flint's acb.integral, the truncation bounded);
   for ρ > 1 and ρ < −1 the two-dimensional expectation by an M × M Gauss–Hermite product rule (not rigorous; M in R_M).
S. The exact errors e_K at every K = S_K0 … K_max of A6's 14 rules without a frame (python-flint arb at 2K·log₂(1/t(v)) + S_MARGIN
   bits; the nodes by Newton from LAPACK's, as referee round 129's arberr.py — independent of A6's decimal code), stored to S_DIGITS
   digits with the precision and the arb radius; the last K recomputed with S_RECHECK more bits.
A. A6's stored errors (results/chaos/width_mismatch_raw.npz, read-only): its 36 rules' e_{K−1}, e_K at K in A6_K, copied into this raw.
F. (computed by the report from S and A, not stored) the limit-free fits of e_K²/(C t^{2K}) on three windows, the L = 1 fit's s₁ against
   (d)'s (proved since referee round 139); the second ridge's fits (g₀ against C_o/C, g₁/g₀ against the formal s₁′); S against A6 where they overlap; A6 at
   K = 160 for all 36 rules (frames with D‴.61 (iv)'s constant, which C(v) equals there: the identity at F_PREC bits; the cusp excluded).
T. Lemma T's margin (floating point; a sanity check of a proved lemma): sup of Φ⁺ over {max|y_i| ≥ 2} against 2 log t(v), on a grid
   and by the proof's reduction (R1's corners, R2's edge), at A6's 10 two-dimensional rules and on a sweep.
G. Assumption (F) of (d) on the axes (numerical): K[log ĝ_K(√Kη) − Kη²/4 − ½log(1 − η²/4)] at K in G_K, extrapolated in 1/K by the
   report against round 129's F₁(η²/4) (referee round 132's check_F1.py).
Round 132 added to R the rule ρ = −1 (rigorous, one-dimensional on the real axis) and pinned small product rules (R_PIN), and to S
   the rules S_NEW (ρ = −1; moment matching at τ²/(2σ²) = 1 and 4, Kuo–Woźniakowski's 2αγ²), reported in F5.

Usage: python experiments/exp_frameless.py [--only R,S,A,T,G] [--from-raw] → results/chaos/frameless.md (+ _raw.npz). FL_OUT moves both.
--only S resumes: the rules already in the raw's S are kept.
"""
import os
import sys
import math
import time
import warnings
from fractions import Fraction

import numpy as np

OUT = os.environ.get('FL_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
A6_RAW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'width_mismatch_raw.npz')

R_MM = (('1', 6), ('2', 10), ('0.25', 5))                              # R: moment matching (τ/σ, K), one-dimensional, rigorous
R_2D = (('0.25', '2', 6), ('1', '2', 8), ('0.5', '1.25', 6), ('0.25', '0.5', 6), ('1', '0.5', 8), ('0.5', '0.8', 6))  # R: (τ/σ, v/v_hs, K)
R_M = (120, 160)                                                       # R: the product rule's points per axis
R_RM1 = (('1', 6),)                                                    # R: ρ = −1 (τ/σ, K), one-dimensional on the real axis, rigorous
R_PIN = (('0.5', '1.25', 3, (40, 50)), ('1', '2', 4, (24, 32)), ('1', '0.5', 4, (24, 32)))  # R: small product rules pinned at 60 digits
R_BITS = 320                                                           # R: the working precision
R_LCUT = 200                                                           # R: the 1-D integral on |W| ≤ (R_LCUT τ²)^{1/2} + 5
S_A6 = (('0.25', 'mm', 200), ('0.5', 'mm', 300), ('1', 'mm', 400), ('2', 'mm', 400),
           ('0.25', '1.25', 240), ('0.25', '2', 400), ('0.5', '1.25', 300), ('0.5', '2', 360), ('1', '2', 400),
           ('0.25', '0.5', 300), ('0.25', '0.8', 240), ('0.5', '0.5', 360), ('0.5', '0.8', 300), ('1', '0.5', 400))  # S: (τ/σ, v/v_hs or 'mm', K_max)
S_NEW = (('1', 'rho=-1', 300), ('√2', 'mm', 300), ('2√2', 'mm', 300))  # S beyond A6 (round 132): ρ = −1; moment matching at τ²/(2σ²) = 1, 4
S_RULES = S_A6 + S_NEW                                                 # S: A6's 14 rules without a frame, then S_NEW (F5)
S_K0 = 2                                                               # S: the first K
S_MARGIN = 450                                                         # S: bits beyond the cancellation 2K·log₂(1/t(v)) (arberr.py's)
S_DIGITS = 60                                                          # S: the digits stored
S_RECHECK = 128                                                        # S: bits added in the recomputation of e_{K_max}
A6_K = (20, 40, 60, 80, 120, 160)                                      # A: A6's K (e_{K−1}, e_K at each)
TAUS = ('0.25', '0.5', '1', '2')                                       # A6's τ/σ


# ------------------------------------------------------------------------------------------------------------------- arb helpers

def _flint():
    from flint import arb, ctx
    return arb, ctx


def A_(x):
    """x (str, int, Fraction or arb) as an arb at the current precision; a string 'a√b' (a optional) is a·√b."""
    arb, _ = _flint()
    if isinstance(x, arb):
        return x
    if isinstance(x, Fraction):
        return arb(x.numerator) / x.denominator
    if isinstance(x, str) and '√' in x:
        a, b = x.split('√')
        return (arb(a) if a else arb(1)) * arb(b).sqrt()
    return arb(str(x))


def tau2_frac(tau):
    """τ² of a τ/σ string as a Fraction (a decimal, or 'a√b': a²b)."""
    if '√' in tau:
        a, b = tau.split('√')
        return Fraction(a or 1) ** 2 * Fraction(b)
    return Fraction(tau) ** 2


def rule_v(tau, w):
    """The variance v of a scaled rule at τ/σ = tau (σ = 1): w = 'mm' (moment matching, τ²), 'rho=-1' (τ²/σ*², ρ = −1), 'Fock'
    (2τ²/(2 + τ²)), 'Mercer-zero' (τ²/b), 'Mercer frame' (2τ²/(b + 1)), b = (1 + 2τ²)^{1/2}, or a decimal string x (x·v_hs, v_hs =
    τ²/σ*) — A6's rules, exactly."""
    t2 = A_(tau) ** 2
    if w == 'mm':
        return t2
    if w == 'rho=-1':
        return t2 / (1 + t2)
    b = (1 + 2 * t2).sqrt()
    if w == 'Fock':
        return 2 * t2 / (2 + t2)
    if w == 'Mercer-zero':
        return t2 / b
    if w == 'Mercer frame':
        return 2 * t2 / (b + 1)
    return A_(w) * t2 / (1 + t2).sqrt()


def h_(y):
    return (y - 1) / (y + 1)


def t_float(tau, w):
    """t(v) = max(|s_a|, |s_b|) in double (the precision formula's)."""
    arb, ctx = _flint()
    p = ctx.prec
    ctx.prec = 200
    try:
        t2 = A_(tau) ** 2
        v = rule_v(tau, w)
        return max(abs(float(h_(t2 / v).mid())), abs(float(h_(v * (1 + t2) / t2).mid())))
    finally:
        ctx.prec = p


# --------------------------------------------------------------------------------------------------------------- the closed forms

def params(r, v):
    """The model at τ/σ = r and variance v (σ = 1; arb at the current precision): c = (1 − v/τ²)/2, c_b = c − v/2, s_a = c/(1 − c),
    s_b = c_b/(1 − c_b), ρ = 1 + 2/τ² − 2/v (D‴.61 (v)), s_d (the larger modulus; 'dom' names it), s_o, ϖ = s_o/s_d, t = |s_d|; kind
    'frame' (|ρ| < 1), 'cusp' (a frame with ϖ = −1, v = v_hs), 'MM' (ρ = 1), 'rho=-1', 'rho>1' or 'rho<-1' (within 10⁻¹² in double)."""
    t2 = A_(r) ** 2
    v = A_(v)
    ss2 = 1 + t2
    c = (1 - v / t2) / 2
    cb = c - v / 2
    sa, sb = c / (1 - c), cb / (1 - cb)
    rho = 1 + 2 / t2 - 2 / v
    if float(abs(sa).mid()) >= float(abs(sb).mid()):
        sd, so, dom = sa, sb, 'a'
    else:
        sd, so, dom = sb, sa, 'b'
    w = so / sd
    rf = float(rho.mid())
    if abs(rf) < 1 - 1e-12:
        kind = 'cusp' if abs(float((1 + w).mid())) < 1e-12 else 'frame'
    elif abs(rf - 1) <= 1e-12:
        kind = 'MM'
    elif abs(rf + 1) <= 1e-12:
        kind = 'rho=-1'
    else:
        kind = 'rho>1' if rf > 1 else 'rho<-1'
    return dict(c=c, c_b=cb, s_a=sa, s_b=sb, rho=rho, s_d=sd, s_o=so, varpi=w, t=abs(sd), dom=dom, kind=kind, ss=ss2.sqrt(), t2=t2, v=v)


def C_const(r, v, P=None):
    """(b)'s constant, e_K² ~ C(v)t(v)^{2K}: √2σ*v(1 − s_d)/[(τ² + v)(τ² + vσ*²)(1 + ϖ)]^{1/2} (singular at the cusp)."""
    P = P or params(r, v)
    t2, v, ss = P['t2'], P['v'], P['ss']
    return 2 ** A_('0.5') * ss * v * (1 - P['s_d']) / ((t2 + v) * (t2 + v * ss ** 2) * (1 + P['varpi'])).sqrt()


def C_other(r, v, P=None):
    """(c)'s second ridge, J_o ~ C_o(v)|s_o|^{2K}: C(v) with d and o exchanged, √2σ*v(1 − s_o)/[(τ² + v)(τ² + vσ*²)(1 + 1/ϖ)]^{1/2};
    None unless ρ ∉ [−1, 1] (ϖ > 0)."""
    P = P or params(r, v)
    if P['kind'] not in ('rho>1', 'rho<-1'):
        return None
    t2, v, ss = P['t2'], P['v'], P['ss']
    return 2 ** A_('0.5') * ss * v * (1 - P['s_o']) / ((t2 + v) * (t2 + v * ss ** 2) * (1 + 1 / P['varpi'])).sqrt()


def s1_formal(s, w):
    """(d)'s formal 1/K coefficient of a ridge's integral, C t^{2K}(1 + s₁/K + …), s the ridge's signed tilted variance, w = ϖ:
    1/16 + 2s(1 + s)/(1 − s)³ + s/(4(1 − s)²) − ϖ(1 + s)²/(8(1 + ϖ)(1 − s)²) (round 129's formula at ϱ = −s, √q = −ϖ); exact for
    Fractions."""
    one = s * 0 + 1
    return one / 16 + 2 * s * (1 + s) / (1 - s) ** 3 + s / (4 * (1 - s) ** 2) - w * (1 + s) ** 2 / (8 * (1 + w) * (1 - s) ** 2)


def s1_r129(rho_, sq):
    """Round 129's s₁ (review/REVIEW_round129_D3dx.md, its c₂ formula), ϱ = t(v) on the fast branch and −t(v) on the slow, √q signed."""
    one = rho_ * 0 + 1
    return (one / 16 - rho_ / (4 * (1 + rho_) ** 2) - 2 * rho_ * (1 - rho_) / (1 + rho_) ** 3
            + sq * (1 - rho_) ** 2 / (8 * (1 - sq) * (1 + rho_) ** 2))


def s1_unsimplified(s, w):
    """(d)'s expansion before simplification: 1/16 − ϖω/8 + sf′/2 + s²f″/2 + s²f′² + ϖωsf′ + ϖωs²f″ + 2F₁(s), ω = 1/(1 + ϖ),
    f′ = −1/(2(1 − s)), f″ = −1/(2(1 − s)²)."""
    one = s * 0 + 1
    om = one / (1 + w)
    f1, f2 = -one / (2 * (1 - s)), -one / (2 * (1 - s) ** 2)
    return (one / 16 - w * om / 8 + s * f1 / 2 + s * s * f2 / 2 + s * s * f1 * f1 + w * om * s * f1 + w * om * s * s * f2 + 2 * F1(s))


def F1(y):
    """Round 129's F₁(y) = y(1 + y)/(1 − y)³ + y/(8(1 − y)²) + y/(8(1 − y)) (the formal 1/K term of log ĝ_K(√Kη) − Kη²/4, y = η²/4)."""
    return y * (1 + y) / (1 - y) ** 3 + y / (8 * (1 - y) ** 2) + y / (8 * (1 - y))


def c2_formal(s, w):
    """(d): e_{K−1}/e_K = t(v)^{−1}(1 + c₂/K² + O(K⁻³)), c₂ = s₁/2."""
    return s1_formal(s, w) / 2


def B0_formal(r, v, P=None):
    """(c)'s formal transient, e_{K−1}/e_K = R_K^{(d)} + Bq^K(1 + B₁/K + …): (B₀, B₁/B₀, q) = ((C_o/C)(1 − q)/(2tq), s₁′ − s₁, ϖ²),
    s₁′ the other ridge's (d and o exchanged); None without a second ridge."""
    P = P or params(r, v)
    Co = C_other(r, v, P)
    if Co is None:
        return None
    q = P['varpi'] ** 2
    return (Co / C_const(r, v, P) * (1 - q) / (2 * P['t'] * q),
            s1_formal(P['s_o'], 1 / P['varpi']) - s1_formal(P['s_d'], P['varpi']), q)


def rm1_C2(tau2):
    """ρ = −1 (v = τ²/σ*², the dual of moment matching: s_d = s_a = +λ, ϖ = 0): C² = (1 − λ)²/(1 + σ*²), rational for rational τ²."""
    lam = mm_lambda(tau2)
    return (1 - lam) ** 2 / (2 + tau2)


def rm1_s1(tau2):
    """ρ = −1's s₁ = 1/16 + 2λ(1 + λ)/(1 − λ)³ + λ/(4(1 − λ)²) (s = +λ, ϖ = 0): rational for rational τ²."""
    lam = mm_lambda(tau2)
    return s1_formal(lam, Fraction(0))


def mm_lambda(tau2):
    """Moment matching's rate λ = tanh r = τ²/(2 + τ²) (σ = 1), exact for a rational τ²."""
    return tau2 / (2 + tau2)


def mm_C2(tau2):
    """Moment matching's constant squared, C² = (1 + λ)³/2 (C = (1 + λ)^{3/2}/√2 = 2σ*³/(1 + σ*²)^{3/2}): rational for rational τ²."""
    lam = mm_lambda(tau2)
    return (1 + lam) ** 3 / 2


def mm_s1(tau2):
    """Moment matching's s₁ = 1/16 − λ/(4(1 + λ)²) − 2λ(1 − λ)/(1 + λ)³ (s = −λ, ϖ = 0): rational for rational τ²."""
    lam = mm_lambda(tau2)
    return Fraction(1, 16) - lam / (4 * (1 + lam) ** 2) - 2 * lam * (1 - lam) / (1 + lam) ** 3


def d61_constant(r, v):
    """D‴.61 (iv)'s constant at a frame rule off the cusp, e_K² ~ const·r^{2K}: σ*(1 − ρ²)^{1/2}(v/τ²) times (1 − r)/(a(2γx(γy −
    γx))^{1/2}) (γx < γy, r = 1/(2γx)) or (1 + r)/(a(2γy(γx − γy))^{1/2}) (γy < γx, r = 1/(2γy)), a = v(1 − ρ²)/2, γx = A + ½,
    γy = B − ½, A = v/(τ² − v), B = vσ*²/(vσ*² − τ²). Returns (r, the constant)."""
    t2 = A_(r) ** 2
    v = A_(v)
    ss2 = 1 + t2
    rho = 1 + 2 / t2 - 2 / v
    a = v * (1 - rho * rho) / 2                                          # rho * rho: arb's rho ** 2 is nan on a ball about 0 (ρ = 0)
    gx, gy = v / (t2 - v) + A_('0.5'), v * ss2 / (v * ss2 - t2) - A_('0.5')
    if float(gx.mid()) < float(gy.mid()):
        rr, S = 1 / (2 * gx), (1 - 1 / (2 * gx)) / (a * (2 * gx * (gy - gx)).sqrt())
    else:
        rr, S = 1 / (2 * gy), (1 + 1 / (2 * gy)) / (a * (2 * gy * (gx - gy)).sqrt())
    return rr, ss2.sqrt() * (1 - rho * rho).sqrt() * (v / t2) * S


def phi_rate(a, b, sa, sb):
    """Lemma A's exponent in the tilted coordinates (numpy): Φ = 2 + 2 log(|a² − b²|/4) − a²/(2|s_a|) − b²/(2|s_b|)."""
    with np.errstate(divide='ignore'):
        return 2 + 2 * np.log(np.abs(a * a - b * b) / 4) - a * a / (2 * abs(sa)) - b * b / (2 * abs(sb))


# ------------------------------------------------------------------------------------------------------- S: the exact errors (arb)

def gh_nodes(K, bits):
    """The K-point Gauss–Hermite rule of N(0, 1) in arb at `bits` (referee round 129's arberr.py): LAPACK's nodes refined by Newton on
    He_K (He_K′ = KHe_{K−1}) to 2^{−(bits−20)}, the midpoints taken exact (not enclosures), the weights K!/(K He_{K−1}(z_k))²,
    symmetric."""
    from scipy.linalg import eigh_tridiagonal
    arb, ctx = _flint()
    ctx.prec = bits
    z0 = eigh_tridiagonal(np.zeros(K), np.sqrt(np.arange(1, K, dtype=float)), eigvals_only=True)
    zs, ws = [], []
    fK = arb(math.factorial(K))
    tol = arb(2) ** (-(bits - 20))
    for x0 in z0[: (K + 1) // 2]:                                       # the negative half and the middle
        x = arb(float(x0))
        for it in range(40):
            hm, h = arb(1), x
            for n in range(1, K):
                hm, h = h, x * h - n * hm
            dx = h / (K * hm)
            x = arb((x - dx).mid())
            if abs(dx.mid()) < tol * (1 + abs(x.mid())) and it > 2:
                break
        hm, h = arb(1), x
        for n in range(1, K):
            hm, h = h, x * h - n * hm
        zs.append(x)
        ws.append(fK / (K * hm) ** 2)
    if K % 2:
        zs[-1] = arb(0)
        return zs + [-z for z in zs[:-1][::-1]], ws + ws[:-1][::-1]
    return zs + [-z for z in zs[::-1]], ws + ws[::-1]


def err2(tau, w, K, bits, nodes=None):
    """e_K² of the scaled rule (tau, w) in arb at `bits` (arberr.py's closed form): 1 − 2√2(1 + σ*²)^{−1/2}Σ_k a_k e^{−y_k²/(2(1 + σ*²))}
    + σ*⁻¹Σ_{jk}a_ja_k e^{−(y_j − y_k)²/4}, y_k = √v z_k, a_k = σ*w_k(v/τ²)^{1/2}e^{cz_k²}, c = (1 − v/τ²)/2."""
    arb, ctx = _flint()
    zs, ws = nodes if nodes is not None else gh_nodes(K, bits)
    ctx.prec = bits
    tau2 = A_(tau) ** 2
    v = rule_v(tau, w)
    ss2 = 1 + tau2
    ss = ss2.sqrt()
    sv = v.sqrt()
    ys = [sv * z for z in zs]
    g = (1 - v / tau2) / 2
    a = [ss * wk * (v / tau2).sqrt() * (z * z * g).exp() for wk, z in zip(ws, zs)]
    cross = sum((a_ * (-(y * y) / (2 * (ss2 + 1))).exp() for a_, y in zip(a, ys)), arb(0))
    gram = arb(0)
    for j in range(K):
        gram += a[j] * a[j]
        s = arb(0)
        yj = ys[j]
        for k in range(j):
            d = yj - ys[k]
            s += a[k] * (-(d * d) / 4).exp()
        gram += 2 * a[j] * s
    return 1 - 2 * arb(2).sqrt() / (ss2 + 1).sqrt() * cross + gram / ss


def s_bits(tau, w, K, extra=0):
    """S's precision: 2K·log₂(1/t(v)) + S_MARGIN (+ extra) bits."""
    return int(2 * K * math.log2(1 / t_float(tau, w)) + S_MARGIN + extra)


def e_row(tau, w, K, extra=0):
    """(K, e_K to S_DIGITS digits, bits, e_K's relative arb radius, seconds)."""
    t0 = time.time()
    bits = s_bits(tau, w, K, extra)
    E = err2(tau, w, K, bits)
    e = E.sqrt()
    rel = float(e.rad() / abs(e.mid()))
    return (K, e.str(S_DIGITS, radius=False), bits, rel, round(time.time() - t0, 2))


def s_rule(tau, w, kmax):
    """S for one rule: e_K at K = S_K0 … kmax and e_{kmax} again with S_RECHECK more bits."""
    t0 = time.time()
    rows = [e_row(tau, w, K) for K in range(S_K0, kmax + 1)]
    chk = e_row(tau, w, kmax, S_RECHECK)
    return dict(tau=tau, w=w, kmax=kmax, rows=rows, check=chk, secs=round(time.time() - t0, 1))


# ------------------------------------------------------------------------------------------------------- R: the representation (a)

def phi_fns(nodes):
    """φ_K(iw) = e^{−w²/2} − Q_K[cos(wz)] (real) and φ_K(x) = e^{x²/2} − Q_K[e^{xz}] for the rule `nodes`, in arb."""
    arb, _ = _flint()
    zs, ws = nodes

    def phi_i(w):
        return (-(w * w) / 2).exp() - sum((wk * (w * z).cos() for z, wk in zip(zs, ws)), arb(0))

    def phi_r(x):
        return (x * x / 2).exp() - sum((wk * (x * z).exp() for z, wk in zip(zs, ws)), arb(0))
    return phi_i, phi_r


def r_mm(tau, K):
    """R at moment matching: σ*E[φ_K(iB)²], B ~ N(0, τ²/2), by acb.integral on |W| ≤ L = (R_LCUT τ²)^{1/2} + 5 plus the truncation's bound
    (|φ_K(iW)|² ≤ 4 and erfc(L/τ) ≤ e^{−L²/τ²}, added to the radius), against `err2` at the same nodes (the identity holds for any rule)."""
    from flint import acb
    arb, ctx = _flint()
    t0 = time.time()
    nodes = gh_nodes(K, R_BITS)
    ctx.prec = R_BITS
    t2 = A_(tau) ** 2
    ss = (1 + t2).sqrt()
    zc, wc = [acb(z) for z in nodes[0]], [acb(w) for w in nodes[1]]

    def f(W, analytic):
        p = (-(W * W) / 2).exp() - sum((wk * (W * z).cos() for z, wk in zip(zc, wc)), acb(0))
        return p * p * (-(W * W) / t2).exp()
    L = (t2 * R_LCUT).sqrt() + 5
    tail = (4 * (-(L * L) / t2).exp()).upper()
    rep = ss * (acb.integral(f, -L, L).real / (arb.pi() * t2).sqrt() + arb(0, tail))
    ex = err2(tau, 'mm', K, R_BITS, nodes)
    ctx.prec = R_BITS
    rel = (rep - ex) / ex
    return dict(tau=tau, K=K, rep=rep.str(30, radius=False), ex=ex.str(30, radius=False), rel=float(abs(rel).upper()),
                overlap=bool(rep.overlaps(ex)), tail=fs(tail, 2), secs=round(time.time() - t0, 1))


def r_2d(tau, x, K, Ms=None, digits=30):
    """R for ρ > 1 (a = iA, b = iB, φ on the imaginary axis) or ρ < −1 (a, b real): σ*(v/τ²)E[φ_K(a + b)φ_K(a − b)] by the M × M
    Gauss–Hermite product rule of N(0, 1)² (M in Ms, default R_M; not rigorous) against `err2` at the same nodes; the values to `digits`
    digits (R_PIN: 60, the tests' pinned values)."""
    arb, ctx = _flint()
    t0 = time.time()
    nodes = gh_nodes(K, R_BITS)
    ctx.prec = R_BITS
    P = params(tau, rule_v(tau, x))
    assert P['kind'] in ('rho>1', 'rho<-1'), P['kind']
    phi_i, phi_r = phi_fns(nodes)
    phi = phi_i if P['kind'] == 'rho>1' else phi_r
    ka, kb = abs(P['c']).sqrt(), abs(P['c_b']).sqrt()
    vals = []
    for M in (Ms or R_M):
        X, WX = gh_nodes(M, R_BITS)
        ctx.prec = R_BITS
        tot = arb(0)
        for xi, wi in zip(X, WX):
            a = ka * xi
            for yj, wj in zip(X, WX):
                b = kb * yj
                tot += wi * wj * phi(a + b) * phi(a - b)
        vals.append(P['ss'] * (P['v'] / P['t2']) * tot)
    ex = err2(tau, x, K, R_BITS, nodes)
    ctx.prec = R_BITS
    return dict(tau=tau, x=x, K=K, kind=P['kind'], vals=[y.str(digits, radius=False) for y in vals], ex=ex.str(digits, radius=False),
                rel=float(abs((vals[-1] - ex) / ex).upper()), dM=float(abs((vals[-1] - vals[0]) / ex).upper()),
                secs=round(time.time() - t0, 1))


def r_rm1(tau, K):
    """R at ρ = −1 (b ≡ 0, a ~ N(0, c) real, c = τ²/(2σ*²)): σ*(v/τ²)E[φ_K(a)²] = E[φ_K(a)²]/σ* by acb.integral on |x| ≤ L =
    (R_LCUT/β)^{1/2} + 5, β = 1/(2c) − 1 = 1/τ², plus the truncation's bound (0 ≤ φ_K(x) ≤ e^{x²/2}, D‴.61 (i), so the tail is at most
    (2cβ)^{−1/2}e^{−βL²}, added to the radius), against `err2` at the same nodes."""
    from flint import acb
    arb, ctx = _flint()
    t0 = time.time()
    nodes = gh_nodes(K, R_BITS)
    ctx.prec = R_BITS
    t2 = A_(tau) ** 2
    ss2 = 1 + t2
    c = t2 / (2 * ss2)
    beta = 1 / (2 * c) - 1
    zc, wc = [acb(z) for z in nodes[0]], [acb(w) for w in nodes[1]]

    def f(x, analytic):
        p = (x * x / 2).exp() - sum((wk * (x * z).exp() for z, wk in zip(zc, wc)), acb(0))
        return p * p * (-(x * x) / (2 * c)).exp()
    L = (R_LCUT / beta).sqrt() + 5
    tail = ((2 * c * beta).rsqrt() * (-(beta * L * L)).exp()).upper()
    rep = (acb.integral(f, -L, L).real / (2 * arb.pi() * c).sqrt() + arb(0, tail)) / ss2.sqrt()
    ex = err2(tau, 'rho=-1', K, R_BITS, nodes)
    ctx.prec = R_BITS
    rel = (rep - ex) / ex
    return dict(tau=tau, K=K, rep=rep.str(30, radius=False), ex=ex.str(30, radius=False), rel=float(abs(rel).upper()),
                overlap=bool(rep.overlaps(ex)), tail=fs(tail, 2), secs=round(time.time() - t0, 1))


def part_R():
    return dict(mm=[r_mm(tau, K) for tau, K in R_MM], two=[r_2d(tau, x, K) for tau, x, K in R_2D],
                rm1=[r_rm1(tau, K) for tau, K in R_RM1], pin=[dict(r_2d(tau, x, K, Ms, 60), M=list(Ms)) for tau, x, K, Ms in R_PIN])


# ------------------------------------------------------------------------------------------------------------- A: A6's stored errors

def a6_w(name):
    """A6's rule name as this script's w: 'moment matching' → 'mm', 'x·v_hs' → 'x', the named scalings unchanged."""
    if name == 'moment matching':
        return 'mm'
    if name.endswith('·v_hs'):
        return name[:-len('·v_hs')]
    return name


def a6_tau(r):
    """A6's τ/σ (a float) as this script's string (TAUS)."""
    return {float(t): t for t in TAUS}[float(r)]


def part_A():
    """A: A6's 36 records (results/chaos/width_mismatch_raw.npz, read-only), their e_{K−1}, e_K strings at K in A6_K."""
    z = np.load(A6_RAW, allow_pickle=True)
    out = []
    for r in z['A6'][0]:
        out.append(dict(r=float(r['r']), name=str(r['name']), v=float(r['v']),
                        rows=[dict(K=int(y['K']), prec=int(y['prec']), e0=str(y['e0']), e1=str(y['e1'])) for y in r['rows']]))
    return out


# --------------------------------------------------------------------------------------------- T: Lemma T's margin (floating point)

T_GRID = (1601, 7.0)                                                   # T: A6's rules, grid points per axis and the half-width
T_SWEEP = (801, 9.0)                                                   # T: the sweep's grid
T_TAUS = (0.05, 20.0, 25)                                              # T: the sweep's τ/σ (geometric)
T_WN = 8                                                               # T: v/v_hs per branch in the sweep: (0.02, 0.999/σ*) and (1.001σ*, 50)
T_EDGE = 200001                                                        # T: R2's edge, points on [0, 2]
T_RULES = (('0.25', '0.5'), ('0.25', '0.8'), ('0.5', '0.5'), ('0.5', '0.8'), ('1', '0.5'),
           ('0.25', '1.25'), ('0.25', '2'), ('0.5', '1.25'), ('0.5', '2'), ('1', '2'))  # T: A6's two-dimensional frameless rules


def L_rate(y):
    """L(y) = 1 + log(y²/4) − y²/4 (numpy; −∞ at 0): the rate of (−1)^Kφ_K(i√Ky), ≤ 0, increasing on (0, 2], zero only at |y| = 2."""
    with np.errstate(divide='ignore'):
        return 1 + np.log(y * y / 4) - y * y / 4


def t_sv(tau, w):
    """(c, c_b, s_a, s_b, t) in double for T."""
    t2 = float(tau) ** 2
    ss2 = 1 + t2
    v = (t2 if w == 'mm' else float(w) * t2 / math.sqrt(ss2))
    c = (1 - v / t2) / 2
    cb = c - v / 2
    return c, cb, c / (1 - c), cb / (1 - cb)


def t_reduced(sa, sb, n=T_EDGE):
    """Lemma T's proof in floating point, (k_a, k_b) = (|s_a|, |s_b|)/(1 − ·) — (|c|, |c_b|) for ρ > 1, (c̃, c̃_b) for ρ < −1 —,
    Q = (y₁ + y₂)²/(8k_a) + (y₁ − y₂)²/(8k_b): per closed sector (a: y₁y₂ ≥ 0, b: y₁y₂ ≤ 0) R1's bound −2/k (Q's least value at the
    corner) and R2's, the largest of Φ(2, y₂) = L(y₂) − Q(2, y₂) on its half of the edge; returns {sector: (R1, R2)}."""
    ka, kb = abs(sa) / (1 - abs(sa)), abs(sb) / (1 - abs(sb))
    y = np.linspace(0, 2, n)[1:]
    out = {}
    for sec, sgn, k in (('a', 1, ka), ('b', -1, kb)):
        y2 = sgn * y
        e = L_rate(y2) - ((2 + y2) ** 2 / (8 * ka) + (2 - y2) ** 2 / (8 * kb))
        out[sec] = (-2 / k, float(np.max(e)))
    return out


def t_grid(tau, w, n, R):
    """The note's check in the representation's own variables: Φ⁺ = L⁺(y₁) + L⁺(y₂) − Q on an n × n grid of [−R, R]², Q with (|c|, |c_b|),
    L⁺ = L on |y| ≤ 2 and 0 beyond on the imaginary axis (ρ > 1), L + y²/2 and y²/2 beyond on the real axis (ρ < −1, unsubtracted).
    Returns the global maximum, the |y₁| where it is attained, the tail's {max|y_i| ≥ 2} maximum and the tail's per sector."""
    c, cb, sa, sb = t_sv(tau, w)
    g = np.linspace(-R, R, n)
    Y1, Y2 = np.meshgrid(g, g)
    Q = (Y1 + Y2) ** 2 / (8 * abs(c)) + (Y1 - Y2) ** 2 / (8 * abs(cb))
    real = c > 0

    def Lp(y):
        a = np.abs(y)
        out = np.where(a <= 2, L_rate(np.minimum(a, 2)), 0.0)
        return out + (y * y / 2 if real else 0)
    F = Lp(Y1) + Lp(Y2) - Q
    tail = np.maximum(np.abs(Y1), np.abs(Y2)) >= 2
    i = np.unravel_index(np.argmax(F), F.shape)
    return dict(gmax=float(F[i]), gy=float(abs(Y1[i])), tail=float(F[tail].max()),
                tail_a=float(F[tail & (Y1 * Y2 >= 0)].max()), tail_b=float(F[tail & (Y1 * Y2 <= 0)].max()))


def t_rule(tau, w, n, R):
    c, cb, sa, sb = t_sv(tau, w)
    red = t_reduced(sa, sb)
    G = t_grid(tau, w, n, R)
    t = max(abs(sa), abs(sb))
    return dict(tau=tau, w=w, kind='rho>1' if c < 0 else 'rho<-1', sa=sa, sb=sb, t=t, two_log_t=2 * math.log(t), red=red, grid=G)


def t_sweep_rules():
    """The sweep's (τ/σ, v/v_hs) strings: T_TAUS geometric τ/σ, T_WN geometric v/v_hs in (0.02, 0.999/σ*) and in (1.001σ*, 50)."""
    out = []
    for tau in np.geomspace(*T_TAUS):
        ss = math.sqrt(1 + tau ** 2)
        for w in list(np.geomspace(0.02, 0.999 / ss, T_WN)) + list(np.geomspace(1.001 * ss, 50, T_WN)):
            out.append((repr(float(tau)), repr(float(w))))
    return out


def part_T():
    rules = [t_rule(tau, w, *T_GRID) for tau, w in T_RULES]
    sweep = [t_rule(tau, w, *T_SWEEP) for tau, w in t_sweep_rules()]
    return dict(rules=rules, sweep=sweep)

# ------------------------------------------------------------------------------------------------- G: assumption (F) on the axes

G_K = (60, 80, 120, 160, 240, 320)                                     # G: the K
G_ETA = (('real', '1'), ('real', '1.5'), ('real', '1.8'), ('imag', '1'), ('imag', '1.5'), ('imag', '1.9'))  # G: (axis, |η|)
G_MARGIN = 200                                                         # G: bits 2.5K + G_MARGIN
G_DIGITS = 50                                                          # G: the digits stored


def g_D(nodes, K, axis, eta):
    """D_K(η) of G at the rule `nodes` (arb at the current precision): log ĝ_K = log((2K)!/K!) − 2K log ξ + log φ_K(ξ) on the real axis,
    with (−1)^Kφ_K(iξ) on the imaginary, ξ = √K|η|."""
    arb, _ = _flint()
    phi_i, phi_r = phi_fns(nodes)
    e = A_(eta)
    xi = arb(K).sqrt() * e
    lf = arb(math.factorial(2 * K)).log() - arb(math.factorial(K)).log()
    if axis == 'real':
        return K * (lf - 2 * K * xi.log() + phi_r(xi).log() - K * e * e / 4 - (1 - e * e / 4).log() / 2)
    return K * (lf - 2 * K * xi.log() + ((-1) ** K * phi_i(xi)).log() + K * e * e / 4 - (1 + e * e / 4).log() / 2)


def part_G():
    """G: D_K(η) at G_K × G_ETA, stored to G_DIGITS digits."""
    arb, ctx = _flint()
    t0 = time.time()
    D = {g: [] for g in G_ETA}
    for K in G_K:
        bits = int(2.5 * K + G_MARGIN)
        nodes = gh_nodes(K, bits)
        ctx.prec = bits
        for g in G_ETA:
            D[g].append(g_D(nodes, K, *g).str(G_DIGITS, radius=False))
    return dict(K=list(G_K), rows=[dict(axis=a, eta=e, D=D[(a, e)]) for a, e in G_ETA], secs=round(time.time() - t0, 1))


# ------------------------------------------------------------------------------------- F: the fits (computed by the report, not stored)

F_PREC = 1000                                                          # F: bits of the fits and of the closed forms
F_WIN = ((2, 0, 10), (2, 40, 14), (4, 0, 18))                          # F: limit-free windows [max(60 if j = 4, K_max//j − d), K_max], n powers
F_QRES = 1e-28                                                         # F: the fits carry q^K(g₀ + … + g_F_M K^{−F_M}) where q^{K_min} exceeds this
F_M = 3                                                                # F: the q^K terms' powers in the limit-free fits
F_RIDGE = ((18, 3), (22, 4))                                           # F: the second ridge's (n, m) on [k₀, K_max] and [max(20, k₀//2), K_max]
F_K0DIG = 30                                                           # F: k₀ = max(20, min(K_max//4, F_K0DIG/log₁₀(1/q)))
F_RES = (1e-5, 1e-2)                                                   # F: resolved / roughly resolved by the better window's |g₀/(C_o/C) − 1|
F_A6K = 160                                                            # F: A6's comparison at this K


def f_solve(Ks, ys, cols):
    """Least squares of ys (midpoints taken exact) on the columns, by the normal equations in arb (an enclosure of that solution)."""
    from flint import arb_mat
    arb, _ = _flint()
    M = arb_mat([[f(K) for f in cols] for K in Ks])
    b = arb_mat([[arb(y.mid())] for y in ys])
    x = (M.transpose() * M).solve(M.transpose() * b)
    return [x[i, 0] for i in range(len(cols))]


def f_cols(n, m, q):
    """1/K … 1/K^n, then q^K/K^j for j = 0 … m (m < 0: none)."""
    arb, _ = _flint()
    cols = [lambda K, j=j: arb(1) / arb(K) ** j for j in range(1, n + 1)]
    return cols + [lambda K, j=j: q ** K / arb(K) ** j for j in range(0, m + 1)]


def f_rule(rec):
    """F for one rule of S: the closed forms, the limit-free fits of y_K = e_K²/(C t^{2K}) (L + Σ s_jK^{−j} [+ q^K terms]) on F_WIN's
    windows, the L = 1 fit's s₁ and s₂ on the last, and the second ridge's fits of y_K − 1 (F_RIDGE) where one exists."""
    arb, ctx = _flint()
    ctx.prec = F_PREC
    tau, w, kmax = rec['tau'], rec['w'], rec['kmax']
    P = params(tau, rule_v(tau, w))
    C = C_const(tau, P['v'], P)
    t, q = P['t'], P['varpi'] ** 2
    e = {int(r[0]): arb(r[1]) for r in rec['rows']}
    y = {K: e[K] ** 2 / (C * t ** (2 * K)) for K in e}
    Co = C_other(tau, P['v'], P)
    two = Co is not None
    s1 = s1_formal(P['s_d'], P['varpi'])
    out = dict(tau=tau, w=w, kind=P['kind'], kmax=kmax, rho=P['rho'], t=t, C=C, s1=s1, c2=s1 / 2, fits=[])
    for j, d, n in F_WIN:
        kmin = max(60, kmax // j) if j == 4 else kmax // j - d
        m = F_M if (two and float(q.mid()) ** kmin > F_QRES) else -1
        Ks = [K for K in sorted(e) if K >= kmin]
        cols = f_cols(n, m, q)
        x = f_solve(Ks, [y[K] for K in Ks], [lambda K: arb(1)] + cols)
        x1 = f_solve(Ks, [y[K] - 1 for K in Ks], cols)
        out['fits'].append(dict(kmin=kmin, n=n, m=m, L1=abs(x[0] - 1), s1=x1[0], s2=x1[1]))
    if two:
        ratio = Co / C
        B0, B1, _ = B0_formal(tau, P['v'], P)
        s1o = s1_formal(P['s_o'], 1 / P['varpi'])
        k0 = int(max(20, min(kmax // 4, F_K0DIG / -math.log10(float(q.mid())))))
        rid = []
        for (n, m), kmin in zip(F_RIDGE, (k0, max(20, k0 // 2))):
            Ks = [K for K in sorted(e) if K >= kmin]
            x = f_solve(Ks, [y[K] - 1 for K in Ks], f_cols(n, m, q))
            g0, g1 = x[n], x[n + 1]
            rid.append(dict(kmin=kmin, n=n, m=m, g0=abs(g0 / ratio - 1), g1g0=g1 / g0, b1=g1 / g0 - x[0], s1=x[0], s2=x[1],
                            last=abs(x[n - 1]) / arb(kmin) ** n, qk=ratio * q ** kmin))
        best = min(float(r['g0'].mid()) for r in rid)
        verdict = 'resolved' if best < F_RES[0] else 'roughly resolved' if best < F_RES[1] else 'not resolvable'
        out.update(so=abs(P['s_o']), q=q, ratio=ratio, B0=B0, B1=B1, s1o=s1o, ridge=rid, verdict=verdict, best=best)
    return out


def f_overlap(S, A):
    """S against A6 where they overlap: the largest |e_S/e_A6 − 1| over every (rule, K) of the 14 rules, and the number of points."""
    arb, ctx = _flint()
    ctx.prec = F_PREC
    worst, npts = 0.0, 0
    for rec in S:
        e = {int(r[0]): r[1] for r in rec['rows']}
        a = [x for x in A if a6_tau(x['r']) == rec['tau'] and a6_w(x['name']) == rec['w']]
        assert len(a) == 1, (rec['tau'], rec['w'])
        for row in a[0]['rows']:
            for K, s in ((row['K'] - 1, row['e0']), (row['K'], row['e1'])):
                if K in e:
                    worst = max(worst, float(abs(arb(e[K]) / arb(s) - 1).mid()))
                    npts += 1
    return worst, npts


def f_a6(A):
    """A6 at K = F_A6K for all 36 rules: R = e_K²/(C t^{2K}) (frames: C(v) = D‴.61 (iv)'s constant, the identity's |C/D‴.61 − 1|),
    R − 1, K(R − 1) and the formal s₁; the cusp excluded."""
    arb, ctx = _flint()
    ctx.prec = F_PREC
    out = []
    for x in A:
        tau, w = a6_tau(x['r']), a6_w(x['name'])
        P = params(tau, rule_v(tau, w))
        rec = dict(tau=tau, name=x['name'], kind=P['kind'], rho=P['rho'], t=P['t'])
        if P['kind'] != 'cusp':
            C = C_const(tau, P['v'], P)
            row = [z for z in x['rows'] if z['K'] == F_A6K][0]
            R = arb(row['e1']) ** 2 / (C * P['t'] ** (2 * F_A6K))
            rec.update(C=C, R1=R - 1, KR1=F_A6K * (R - 1), s1=s1_formal(P['s_d'], P['varpi']))
            if P['kind'] == 'frame':
                r61, C61 = d61_constant(tau, P['v'])
                rec.update(ident=abs(C / C61 - 1).upper(), rdiff=abs(r61 - P['t']).upper())
        out.append(rec)
    return out


def s_is_a6(rec):
    """Whether an S record is one of A6's 14 rules (S_A6) rather than S_NEW's."""
    return (rec['tau'], rec['w']) in {(t, w) for t, w, _ in S_A6}


F_NEWK = (100, 200)                                                    # F5: y_K − 1 printed at these K and at K_max


def f_new(rec):
    """F5 for one rule of S_NEW: `f_rule`'s fits, y_K − 1 = e_K²/(C t^{2K}) − 1 at F_NEWK and K_max, K_max(y − 1), and at moment matching
    the cheap bound σ*(K!/(2K)!)²(4K − 1)!!(τ²/2)^{2K} (from |ĝ_K(iB)| ≤ 1) at K_max against e²_{K_max}."""
    arb, ctx = _flint()
    out = f_rule(rec)
    ctx.prec = F_PREC
    P = params(rec['tau'], rule_v(rec['tau'], rec['w']))
    e = {int(r[0]): arb(r[1]) for r in rec['rows']}
    ks = list(F_NEWK) + [rec['kmax']]
    out['y1'] = [(K, e[K] ** 2 / (out['C'] * out['t'] ** (2 * K)) - 1) for K in ks]
    out['tau2'] = tau2_frac(rec['tau'])
    if P['kind'] == 'MM':
        K = rec['kmax']
        fac = arb(math.factorial(K)) / arb(math.factorial(2 * K))
        dfac = arb(math.factorial(4 * K)) / (arb(2) ** (2 * K) * arb(math.factorial(2 * K)))
        out['cheap'] = P['ss'] * fac ** 2 * dfac * (P['t2'] / 2) ** (2 * K)
        out['e2max'] = e[K] ** 2
    return out


def part_F(S, A):
    """F from S and A (the report's; nothing stored): A6's 14 rules (F1–F4) and S_NEW's (F5)."""
    S6 = [r for r in S if s_is_a6(r)]
    return dict(rules=[f_rule(r) for r in S6], overlap=f_overlap(S6, A) if A else None, a6=f_a6(A) if A else None,
                new=[f_new(r) for r in S if not s_is_a6(r)])

# ------------------------------------------------------------------------------------------------------------------------- report

def fs(x, d=10):
    """An arb's midpoint to d significant digits, the exponent in Python's form (e-05); a ball containing 0 (ρ = 0 at the Fock rules) as 0."""
    if x.contains(0):
        return '0'
    s = x.str(d, radius=False)
    if 'e' in s:
        m, e = s.split('e')
        s = f"{m}e{int(e):+03d}"
    return s


def fe(x, d=1):
    """|x| (float or arb) in e-notation with d decimals."""
    x = abs(float(x.mid()) if hasattr(x, 'mid') else float(x))
    return f"{x:.{d}e}"


def rule_name(w):
    return 'moment matching' if w == 'mm' else 'v = τ²/σ*² (ρ = −1)' if w == 'rho=-1' else w + '·v_hs'


def r_lines(R):
    if not R:
        return []
    out = ['', '## R. The representation (a) at small K', '',
           "e_K² = σ*(v/τ²)E[φ_K(a + b)φ_K(a − b)] **[proved]**, checked against the exact closed form (`err2`) at the same nodes — the "
           "identity holds for every rule, so the comparison is exact up to the integration. Moment matching: σ*E[φ_K(iB)²], B ~ N(0, τ²/2), "
           f"by python-flint's `acb.integral` on |W| ≤ ({R_LCUT}τ²)^{{1/2}} + 5 at {R_BITS} bits, the truncation's bound 4e^{{−L²/τ²}} added "
           "to the radius (rigorous):", '',
           '| τ/σ | K | σ*E[φ_K(iB)²] | the closed form | relative difference (upper bound) | the balls overlap | truncation bound |',
           '|---|---|---|---|---|---|---|']
    for r in R['mm']:
        out.append(f"| {r['tau']} | {r['K']} | {r['rep']} | {r['ex']} | {fe(r['rel'])} | {'yes' if r['overlap'] else 'NO'} | {r['tail']} |")
    out += ['', "ρ > 1 (a = iA, b = iB: φ_K on the imaginary axis) and ρ < −1 (a, b real): the M × M Gauss–Hermite product rule of N(0, 1)² "
            f"at {R_BITS} bits (not rigorous; a check of the normalisation and of the Fubini step):", '',
            f"| τ/σ | rule | kind | K | M = {R_M[0]} | M = {R_M[1]} | the closed form | relative difference (M = {R_M[1]}) | "
            f"change M = {R_M[0]} → {R_M[1]} |", '|---|---|---|---|---|---|---|---|---|']
    for r in R['two']:
        out.append(f"| {r['tau']} | {rule_name(r['x'])} | {r['kind']} | {r['K']} | {r['vals'][0]} | {r['vals'][1]} | {r['ex']} | "
                   f"{fe(r['rel'])} | {fe(r['dM'])} |")
    out += ['', f"Largest relative difference: one-dimensional {fe(max(r['rel'] for r in R['mm']))}, two-dimensional "
            f"{fe(max(r['rel'] for r in R['two']))}; the largest change M = {R_M[0]} → {R_M[1]} {fe(max(r['dM'] for r in R['two']))}."]
    if R.get('rm1'):
        out += ['', "ρ = −1 (b ≡ 0, a ~ N(0, c) real, c = τ²/(2σ*²)): σ*(v/τ²)E[φ_K(a)²] = E[φ_K(a)²]/σ* by `acb.integral` on |x| ≤ "
                f"({R_LCUT}/β)^{{1/2}} + 5, β = 1/(2c) − 1 = 1/τ², at {R_BITS} bits, the truncation's bound (2cβ)^{{−1/2}}e^{{−βL²}} (from "
                "0 ≤ φ_K(x) ≤ e^{x²/2}, D‴.61 (i)) added to the radius (rigorous):", '',
                '| τ/σ | K | E[φ_K(a)²]/σ* | the closed form | relative difference (upper bound) | the balls overlap | truncation bound |',
                '|---|---|---|---|---|---|---|']
        for r in R['rm1']:
            out.append(f"| {r['tau']} | {r['K']} | {r['rep']} | {r['ex']} | {fe(r['rel'])} | {'yes' if r['overlap'] else 'NO'} | {r['tail']} |")
    if R.get('pin'):
        out += ['', "Pinned for the tests (`r_2d` with small product rules, values to 60 digits; not accuracy checks): "
                + '; '.join(f"τ/σ = {r['tau']}, {rule_name(r['x'])}, K = {r['K']}, M = {r['M'][-1]}: {r['vals'][-1]}" for r in R['pin']) + '.']
    return out


def s_lines(S, F):
    if not S:
        return []
    arb, ctx = _flint()
    ctx.prec = F_PREC
    out = ['', '## S. The exact errors at consecutive K **[numerical]**', '',
           f"e_K for K = {S_K0} … K_max (python-flint arb at 2K·log₂(1/t(v)) + {S_MARGIN} bits; referee round 129's arberr.py: the nodes "
           f"Newton iterates from LAPACK's, independent of A6's decimal code), stored to {S_DIGITS} digits; the radius is the arithmetic's "
           f"at those nodes; e_{{K_max}} recomputed with {S_RECHECK} more bits:", '',
           '| τ/σ | rule | kind | t(v) | K | bits | largest relative radius of e_K | e_{K_max} with more bits: relative difference | seconds |',
           '|---|---|---|---|---|---|---|---|---|']
    for rec in S:
        P = params(rec['tau'], rule_v(rec['tau'], rec['w']))
        rows = rec['rows']
        d = abs(arb(rec['check'][1]) / arb(rows[-1][1]) - 1)
        same = f"the same {S_DIGITS} digits" if rec['check'][1] == rows[-1][1] else fe(d)
        out.append(f"| {rec['tau']} | {rule_name(rec['w'])} | {P['kind']} | {fs(P['t'], 10)} | {rows[0][0]}…{rows[-1][0]} | "
                   f"{rows[0][2]}…{rows[-1][2]} | {fe(max(r[3] for r in rows))} | {same} | {rec['secs']} |")
    S6, SN = [r for r in S if s_is_a6(r)], [r for r in S if not s_is_a6(r)]
    out += ['', f"{len(S6)} rules, {sum(len(r['rows']) for r in S6)} values, the largest relative radius "
            f"{fe(max(x[3] for r in S6 for x in r['rows']))}; {sum(r['secs'] for r in S6):.0f} s in all."]
    if SN:
        out.append(f"Beyond A6's rules (F5): {len(SN)} rules, {sum(len(r['rows']) for r in SN)} values, the largest relative radius "
                   f"{fe(max(x[3] for r in SN for x in r['rows']))}; {sum(r['secs'] for r in SN):.0f} s.")
    if F and F.get('overlap'):
        w, n = F['overlap']
        out.append(f"S against A6's stored 26-digit errors where they overlap (K − 1, K at K in A6_K): the largest |e_S/e_A6 − 1| {fe(w)} "
                   f"over {n} values.")
    return out


def model_lines(F):
    rules = F['rules']
    out = ['', '## F. The constants and the fits', '', "### F1. The model's numbers", '',
           "t(v) and C(v) **[proved]** ((b): e_K² ~ C(v)t(v)^{2K}); |s_o|, q = ϖ² and C_o/C **[proved]** ((c): J_o ~ C_o|s_o|^{2K}); s₁ and "
           "c₂ = s₁/2 **[proved]** ((d), by hand since referee round 139: review/RESEARCH_note_c2.md); s₁′ (s₁ with d and o exchanged), "
           "B₀ = (C_o/C)(1 − q)/(2tq) and B₁/B₀ = s₁′ − s₁ **[formal]** ((c), the transient):", '',
           '| τ/σ | rule | ρ | kind | t(v) | C(v) | s₁ | c₂ | \\|s_o\\| | q | C_o/C | s₁′ | B₀ | B₁/B₀ |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rules:
        row = [r['tau'], rule_name(r['w']), fs(r['rho'], 6), r['kind'], fs(r['t'], 10), fs(r['C'], 13), fs(r['s1'], 12), fs(r['c2'], 12)]
        if 'ratio' in r:
            row += [fs(r['so'], 10), fs(r['q'], 10), fs(r['ratio'], 12), fs(r['s1o'], 10), fs(r['B0'], 10), fs(r['B1'], 8)]
        else:
            row += ['—'] * 6
        out.append('| ' + ' | '.join(row) + ' |')
    mm = [r for r in rules if r['kind'] == 'MM']
    if mm:
        out += ['', "At moment matching C = (1 + λ)^{3/2}/√2 and s₁ = 1/16 − λ/(4(1 + λ)²) − 2λ(1 − λ)/(1 + λ)³, rational with τ² (`mm_C2`, "
                "`mm_s1`): " + '; '.join(f"τ/σ = {r['tau']}: C² = {mm_C2(Fraction(r['tau']) ** 2)}, s₁ = {mm_s1(Fraction(r['tau']) ** 2)}, "
                                          f"c₂ = {mm_s1(Fraction(r['tau']) ** 2) / 2}" for r in mm) + '.']
    return out


def fit_lines(F):
    rules = F['rules']
    j_desc = ', '.join((f"[max(60, K_max//{j}), K_max] with {n} powers" if j == 4 else
                        f"[K_max//{j}{' − ' + str(d) if d else ''}, K_max] with {n} powers") for j, d, n in F_WIN)
    out = ['', '### F2. The limit and s₁ from the exact errors **[numerical]**', '',
           f"Least squares (F_PREC = {F_PREC} bits, the data's midpoints exact) of e_K²/(C(v)t(v)^{{2K}}) on L + Σ_{{j ≤ n}}s_jK^{{−j}}, "
           f"plus q^K(g₀ + … + g_{F_M}K^{{−{F_M}}}) where a second ridge exists and q^{{K_min}} > {F_QRES:.0e}, on three windows: {j_desc}. "
           "|L − 1| is the limit's distance from the proved C(v) (its smallest and largest over the windows); s₁ and s₂ are the last "
           "window's with L = 1 fixed:", '',
           '| τ/σ | rule | kind | K | \\|L − 1\\| (three limit-free fits) | s₁ of (d) (proved) | fitted s₁ (L = 1) | relative difference | fitted s₂ |',
           '|---|---|---|---|---|---|---|---|---|']
    order = {'MM': 0, 'rho>1': 1, 'rho<-1': 2}
    for r in sorted(rules, key=lambda r: (order[r['kind']], float(r['tau']), r['w'])):
        Ls = [float(f['L1'].mid()) for f in r['fits']]
        f3 = r['fits'][-1]
        out.append(f"| {r['tau']} | {rule_name(r['w'])} | {r['kind']} | {S_K0}…{r['kmax']} | {min(Ls):.1e} … {max(Ls):.1e} | "
                   f"{fs(r['s1'], 12)} | {fs(f3['s1'], 12)} | {fe(f3['s1'] / r['s1'] - 1)} | {fs(f3['s2'], 8)} |")
    out.append('')
    for kinds, lab in ((('MM',), 'moment matching'), (('rho>1',), 'ρ > 1'), (('rho<-1',), 'ρ < −1')):
        sub = [r for r in rules if r['kind'] in kinds]
        if sub:
            Lmin = [min(float(f['L1'].mid()) for f in r['fits']) for r in sub]
            Lmax = [max(float(f['L1'].mid()) for f in r['fits']) for r in sub]
            ds = [abs(float((r['fits'][-1]['s1'] / r['s1'] - 1).mid())) for r in sub]
            out.append(f"- {lab} ({len(sub)} rules): the best window's |L − 1| from {min(Lmin):.1e} to {max(Lmin):.1e} (every window: "
                       f"at most {max(Lmax):.1e}); s₁ within {min(ds):.1e} … {max(ds):.1e} of (d)'s proved value.")
    return out


def ridge_lines(F):
    rid = [r for r in F['rules'] if 'ridge' in r]
    if not rid:
        return []
    out = ['', '### F3. The second ridge **[numerical]**', '',
           "Least squares of e_K²/(C t^{2K}) − 1 on Σ_{j ≤ n}s_jK^{−j} + q^K Σ_{j ≤ m}g_jK^{−j}: g₀ against C_o/C (**[proved]** by (c)), g₁/g₀ "
           "against the formal s₁′ (J_o's own 1/K coefficient) and g₁/g₀ − s₁ (the fit's) against s₁′ − s₁ = B₁/B₀ (the ratio's transient, "
           f"(c)), on [k₀, K_max] and [max(20, k₀//2), K_max], k₀ = max(20, min(K_max//4, {F_K0DIG}/log₁₀(1/q))). Resolved when the better "
           f"window's |g₀/(C_o/C) − 1| < {F_RES[0]:.0e}, roughly below {F_RES[1]:.0e}, otherwise not resolvable; in brackets the q^K "
           "term's size at the window's start against the last fitted power's, (C_o/C)q^{K_min} vs |s_n|K_min^{−n}:", '',
           '| τ/σ | rule | kind | q (q^{K_max}) | C_o/C | s₁′ | s₁′ − s₁ | fits: [K_min, K_max] n, m: \\|g₀/(C_o/C) − 1\\|, g₁/g₀, g₁/g₀ − s₁ '
           '| verdict |', '|---|---|---|---|---|---|---|---|---|']
    for r in sorted(rid, key=lambda r: (float(r['tau']), r['w'])):
        cells = '; '.join(f"[{x['kmin']}, {r['kmax']}] {x['n']}, {x['m']}: {fe(x['g0'])}, {fs(x['g1g0'], 8)}, {fs(x['b1'], 6)} "
                          f"({fe(x['qk'])} vs {fe(x['last'])})" for x in r['ridge'])
        verdict = r['verdict']
        if verdict == 'not resolvable':
            x = min(r['ridge'], key=lambda x: float(x['g0'].mid()))
            verdict += (f": at the better window's start (C_o/C)q^{{{x['kmin']}}} = {fe(x['qk'])}, below the last fitted power "
                        f"|s_{x['n']}|·{x['kmin']}^{{−{x['n']}}} = {fe(x['last'])}; fitted s₂ = {fs(r['fits'][-1]['s2'], 4)}")
        out.append(f"| {r['tau']} | {rule_name(r['w'])} | {r['kind']} | {fs(r['q'], 8)} ({fs(r['q'] ** r['kmax'], 2)}) | {fs(r['ratio'], 12)} | "
                   f"{fs(r['s1o'], 10)} | {fs(r['B1'], 8)} | {cells} | {verdict} |")
    out.append('')
    for v in ('resolved', 'roughly resolved', 'not resolvable'):
        sub = [r for r in rid if r['verdict'] == v]
        if sub:
            out.append(f"- {v}: {len(sub)} — " + ', '.join(f"({r['tau']}, {rule_name(r['w'])}) {r['best']:.1e}" for r in sub))
    return out


def a6_lines(F):
    a6 = F.get('a6')
    if not a6:
        return []
    out = ['', f"### F4. A6's stored errors at K = {F_A6K}, all 36 rules", '',
           f"R = e_K²/(C(v)t(v)^{{2K}}) at K = {F_A6K} from A6's 26-digit e_K (this raw's A, a copy of results/chaos/width_mismatch_raw.npz's "
           f"A6): R − 1 and {F_A6K}(R − 1) against s₁ ((d)'s, proved for the rules without a frame; at a frame s₁ is round 129's formula, ϖ < 0, formal). At a frame C(v) is "
           f"D‴.61 (iv)'s constant (e), the identity's |C/C₆₁ − 1| at {F_PREC} bits **[proved; checked]**; the cusp (ϖ = −1) is excluded:", '',
           f'| τ/σ | rule | ρ | kind | t(v) | C(v) | \\|C/C₆₁ − 1\\| | R − 1 | {F_A6K}(R − 1) | s₁ |', '|---|---|---|---|---|---|---|---|---|---|']
    for r in a6:
        if r['kind'] == 'cusp':
            out.append(f"| {r['tau']} | {r['name']} | {fs(r['rho'], 6)} | cusp | {fs(r['t'], 10)} | excluded | | | | |")
            continue
        out.append(f"| {r['tau']} | {r['name']} | {fs(r['rho'], 6)} | {r['kind']} | {fs(r['t'], 10)} | {fs(r['C'], 12)} | "
                   f"{fe(r['ident']) if 'ident' in r else '—'} | {fs(r['R1'], 5)} | {fs(r['KR1'], 6)} | {fs(r['s1'], 8)} |")
    fr = [r for r in a6 if r['kind'] == 'frame']
    nf = [r for r in a6 if r['kind'] not in ('frame', 'cusp')]
    cu = [r for r in a6 if r['kind'] == 'cusp']
    out += ['', f"{len(nf)} rules without a frame, {len(fr)} frame rules off the cusp (the identity's largest |C/C₆₁ − 1| "
            f"{fe(max(float(r['ident']) for r in fr))}, |r₆₁ − t(v)| {fe(max(float(r['rdiff']) for r in fr))}), {len(cu)} at the cusp "
            f"(excluded). Without a frame |R − 1| from {min(abs(float(r['R1'].mid())) for r in nf):.1e} to "
            f"{max(abs(float(r['R1'].mid())) for r in nf):.1e}, and |{F_A6K}(R − 1)/s₁ − 1| from "
            f"{min(abs(float((r['KR1'] / r['s1'] - 1).mid())) for r in nf):.1e} to {max(abs(float((r['KR1'] / r['s1'] - 1).mid())) for r in nf):.1e} "
            f"(the O(1/K) remainder s₂/K)."]
    return out


def new_lines(F):
    new = F.get('new')
    if not new:
        return []
    out = ['', "### F5. Beyond A6's rules: ρ = −1, and moment matching at and beyond Kuo–Woźniakowski's threshold **[numerical]**", '',
           "S_NEW's exact errors against (b): y_K = e_K²/(C(v)t(v)^{2K}), y_K − 1 at a few K, K_max(y_{K_max} − 1) against (d)'s proved s₁, "
           "and the limit-free fits of F2. Kuo–Woźniakowski's parametrisation (the kernel e^{−γ²(x − x′)²}, the weight N(0, α)) is "
           "α = τ², γ² = 1/(4σ²) for D‴.58's kernel e^{−(y − y′)²/(4σ²)} and N(0, τ²), so 2αγ² = τ²/(2σ²) and moment matching's rate is "
           "λ = 2αγ²/(1 + 2αγ²); ρ = −1 (v = τ²/σ*²) is moment matching's dual (u ↦ −u), t(v) = λ with s_d = +λ, on the real axis:", '',
           '| τ/σ | rule | kind | τ²/(2σ²) = 2αγ² | t(v) | C(v) | C², s₁ (exact) | K | y_K − 1 | K_max(y − 1) | \\|L − 1\\| (three fits) '
           '| fitted s₁ (L = 1): relative difference |', '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in new:
        t2 = r['tau2']
        ex = (f"{mm_C2(t2)}, {mm_s1(t2)}" if r['kind'] == 'MM' else f"{rm1_C2(t2)}, {rm1_s1(t2)}" if r['kind'] == 'rho=-1' else '—')
        Ls = [float(f['L1'].mid()) for f in r['fits']]
        y1 = '; '.join(f"K = {K}: {float(y.mid()):.4e}" for K, y in r['y1'])
        K, y = r['y1'][-1]
        out.append(f"| {r['tau']} | {rule_name(r['w'])} | {r['kind']} | {t2 / 2} | {fs(r['t'], 10)} | {fs(r['C'], 13)} | {ex} | "
                   f"{S_K0}…{r['kmax']} | {y1} | {fs(K * y, 8)} (s₁ = {fs(r['s1'], 8)}) | {min(Ls):.1e} … {max(Ls):.1e} | "
                   f"{fs(r['fits'][-1]['s1'], 12)}: {fe(r['fits'][-1]['s1'] / r['s1'] - 1)} |")
    mm = [r for r in new if 'cheap' in r]
    if mm:
        out += ['', "At moment matching the cheap bound e_K² ≤ σ*(K!/(2K)!)²(4K − 1)!!(τ²/2)^{2K} ~ 2^{−1/2}σ*(τ²/2)^{2K} (from |ĝ_K(iB)| ≤ 1; "
                "its rate τ²/2 = 2αγ² is Kuo–Woźniakowski's threshold) at K_max: "
                + '; '.join(f"τ/σ = {r['tau']} (2αγ² = {r['tau2'] / 2}): the bound {fs(r['cheap'], 4)}, e²_{{K_max}} = {fs(r['e2max'], 4)}, "
                            f"λ = {fs(r['t'], 6)}" for r in mm) + '.']
    return out


def g_lines(G):
    if not G:
        return []
    arb, ctx = _flint()
    ctx.prec = F_PREC
    Ks = G['K']
    out = ['', "## G. Assumption (F) on the axes **[numerical]**", '',
           "D_K(η) = K[log ĝ_K(√Kη) − Kη²/4 − ½log(1 − η²/4)] on the real axis and K[log ĝ_K(i√Ky) + Ky²/4 − ½log(1 + y²/4)] on the "
           f"imaginary (η = iy), at K = {', '.join(map(str, Ks))} ({G_DIGITS} digits stored; arb at 2.5K + {G_MARGIN} bits, `gh_nodes`), "
           "extrapolated to K = ∞ by the polynomial in 1/K through all the points and through the last four, against round 129's F₁(η²/4) "
           "((d)'s assumption (F): D_K → F₁). Near η = 2 on the real axis the approach is slow (pre-asymptotic, not a contradiction):", '',
           f'| axis | η | η²/4 | D_K at K = {Ks[-1]} | extrapolated ({len(Ks)} points) | extrapolated (last 4) | F₁(η²/4) | '
           f'\\|extrapolated ({len(Ks)}) − F₁\\| |', '|---|---|---|---|---|---|---|---|']
    for g in G['rows']:
        e = A_(g['eta'])
        y = e * e / 4 * (1 if g['axis'] == 'real' else -1)
        D = [arb(x) for x in g['D']]
        x6, x4 = g_extrap(Ks, D), g_extrap(Ks[-4:], D[-4:])
        f1 = F1(y)
        out.append(f"| {g['axis']} | {g['eta']} | {fs(y, 6)} | {fs(D[-1], 8)} | {fs(x6, 10)} | {fs(x4, 10)} | {fs(f1, 10)} | {fe(x6 - f1)} |")
    out += ['', f"G computed in {G['secs']} s."]
    return out


def g_extrap(Ks, Ds):
    """The value at 1/K = 0 of the polynomial in 1/K through the points (K, D_K)."""
    from flint import arb_mat
    arb, _ = _flint()
    n = len(Ks)
    M = arb_mat([[arb(1) / arb(K) ** j for j in range(n)] for K in Ks])
    return M.solve(arb_mat([[d] for d in Ds]))[0, 0]


def t_lines(T):
    if not T:
        return []
    out = ['', "## T. Lemma T's margin (floating point; a sanity check of a proved lemma) **[numerical]**", '',
           "sup of Φ⁺ = L⁺(y₁) + L⁺(y₂) − Q over the tail {max|y_i| ≥ 2} against 2 log t(v): on a grid in the representation's own variables "
           f"(Q with (|c|, |c_b|), the real axis unsubtracted; {T_GRID[0]}² points on [−{T_GRID[1]:g}, {T_GRID[1]:g}]²) and by the proof's "
           "reduction (R1: −2/k at the corner; R2: the largest Φ(2, y₂) on the edge, "
           f"{T_EDGE} points); the sector forms against 2 log|s| of the sector's own variable:", '',
           '| τ/σ | rule | kind | t(v) | 2 log t(v) | grid max − 2 log t | at \\|y\\| (2t^{1/2}) | margin (grid) | margin (reduction) | R1, R2 '
           '| sector margins (d, o) |', '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in T['rules']:
        out.append(t_row(r))
    sw = T['sweep']
    red = [t_margins(r) for r in sw]
    i = int(np.argmin([m[0] for m in red]))
    j = int(np.argmin([m[1] for m in red]))
    k = int(np.argmin([min(m[2], m[3]) for m in red]))
    gmax = max(r['grid']['gmax'] - r['two_log_t'] for r in sw)
    out += ['', f"The sweep: τ/σ in [{T_TAUS[0]:g}, {T_TAUS[1]:g}] ({T_TAUS[2]} geometric), v/v_hs in {T_WN} geometric values of each "
            f"frameless range ((0.02, 0.999/σ*) and (1.001σ*, 50)), {len(sw)} rules, grids of {T_SWEEP[0]}² points on "
            f"[−{T_SWEEP[1]:g}, {T_SWEEP[1]:g}]²: the smallest margin {red[i][0]:.1e} by the reduction (τ/σ = {float(sw[i]['tau']):.4g}, "
            f"v/v_hs = {float(sw[i]['w']):.4g}, t(v) = {sw[i]['t']:.6f}) and {red[j][1]:.1e} on the grid; the smallest sector margin "
            f"{min(red[k][2], red[k][3]):.1e}; every margin positive: {'yes' if min(min(m) for m in red) > 0 else 'NO'}; the grid's "
            + (f"global maximum never above 2 log t(v) (Lemma A), within {-gmax:.1e} of it at the closest." if gmax <= 0 else
               f"global maximum above 2 log t(v) by {gmax:.1e} (Lemma A: NO).")]
    return out


def t_margins(r):
    red, G = r['red'], r['grid']
    sup = max(max(red['a']), max(red['b']))
    d, o = ('a', 'b') if abs(r['sa']) >= abs(r['sb']) else ('b', 'a')
    sd, so = (r['sa'], r['sb']) if d == 'a' else (r['sb'], r['sa'])
    return (r['two_log_t'] - sup, r['two_log_t'] - G['tail'], 2 * math.log(abs(sd)) - max(red[d]), 2 * math.log(abs(so)) - max(red[o]))


def t_row(r):
    G = r['grid']
    m = t_margins(r)
    red = r['red']
    d = 'a' if abs(r['sa']) >= abs(r['sb']) else 'b'
    return (f"| {r['tau']} | {rule_name(r['w'])} | {r['kind']} | {r['t']:.10f} | {r['two_log_t']:+.6f} | {G['gmax'] - r['two_log_t']:+.1e} | "
            f"{G['gy']:.4f} ({2 * math.sqrt(r['t']):.4f}) | {m[1]:.4f} | {m[0]:.4f} | {max(red[d][0], red['b' if d == 'a' else 'a'][0]):+.4f}, "
            f"{max(red['a'][1], red['b'][1]):+.4f} | {m[2]:.4f}, {m[3]:.4f} |")


def report(R=None, S=None, A=None, T=None, G=None):
    F = part_F(S, A) if S else None
    head = ['# Proposition D‴.62: the scaled Gauss–Hermite rules without a Mehler frame (`exp_frameless.py`)', '',
            "D‴.58's scaled rule (σ = 1): the K-point Gauss–Hermite rule of N(0, v) for the mixing law reweighted to N(0, τ²), e_K = "
            "‖f* − f_K‖/‖f*‖; c = (1 − v/τ²)/2, c_b = c − v/2, s_a = c/(1 − c), s_b = c_b/(1 − c_b), s_d the larger in modulus, s_o the "
            "other, ϖ = s_o/s_d, t(v) = |s_d|, C(v) = √2σ*v(1 − s_d)/[(τ² + v)(τ² + vσ*²)(1 + ϖ)]^{1/2}. Labels: **[proved]** (D‴.62), "
            "**[formal]** (an expansion not proved), **[numerical]**. Sources: research note, 2026-11-02 (review/RESEARCH_note_frameless.md); "
            "the exact-error method: referee round 129's arberr.py."]
    body = r_lines(R) + s_lines(S, F)
    if F:
        body += model_lines(F) + fit_lines(F) + ridge_lines(F) + a6_lines(F) + new_lines(F)
    body += t_lines(T) + g_lines(G)
    return '\n'.join(head + body + ['', 'Numbers are produced by `experiments/exp_frameless.py`; nothing is edited by hand.', ''])

# ------------------------------------------------------------------------------------------------------------------ raw handling

def save(path, **parts):
    np.savez_compressed(path, **{k: np.array([v], dtype=object) for k, v in parts.items() if v is not None})


def load(path):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    out = {k: z[k][0] for k in z.files}
    return {k: (list(v) if isinstance(v, np.ndarray) else v) for k, v in out.items()}   # a stored list of dicts comes back as an array


# ------------------------------------------------------------------------------------------------------------------------ main

KEYS = ('R', 'S', 'A', 'T', 'G')


def part_S(old, raw, res):
    """S: `s_rule` for every rule of S_RULES not yet in the raw's S (resumable: the raw is saved after each rule)."""
    have = {(r['tau'], r['w'], r['kmax']): r for r in (old.get('S') or [])}
    out = []
    for tau, w, kmax in S_RULES:
        r = have.get((tau, w, kmax))
        if r is None:
            r = s_rule(tau, w, kmax)
            print(f"  S {tau} {w} K ≤ {kmax}: {r['secs']}s", flush=True)
        out.append(r)
        res['S'] = list(out)
        save(raw, **res)
    return out


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'frameless_raw.npz')
    md = os.path.join(OUT, 'frameless.md')
    old = load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w', encoding='utf-8').write(report(*[old.get(k) for k in KEYS]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(KEYS)
    res = dict(old)
    t0 = time.time()
    for k in KEYS:
        if k not in parts:
            continue
        if k == 'S':
            res['S'] = part_S(old, raw, res)
        else:
            res[k] = {'R': part_R, 'A': part_A, 'T': part_T, 'G': part_G}[k]()
        save(raw, **res)
        print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w', encoding='utf-8').write(report(*[res.get(k) for k in KEYS]))


if __name__ == '__main__':
    main()
