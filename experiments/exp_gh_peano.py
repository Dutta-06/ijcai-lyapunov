"""
Proposition D‴.61: the Gauss–Hermite rule's Peano kernel — its Laplace transform, the scaled rules' rates and their sharp constants
====================================================================================================================================

The K-point Gauss–Hermite rule Q_K of N(0, 1) integrates polynomials of degree < 2K exactly and its Peano kernel G_K of order 2K is
nonnegative with ∫G_K = K!/(2K)!, so its defect generating function is ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] = (K!/(2K)!)ξ^{2K}e^{−ξ²/2}ĝ_K(ξ),
ĝ_K the two-sided Laplace transform of the probability density g_K = G_K(2K)!/K!. D‴.58's frame series (the scaled rules' squared errors
up to one constant) are Fock integrals of |ε̃_K|². This script:

P. The defects D(m) = (E − Q_K)[z^m] exactly, as weighted Dyck paths (E z^m = e₁ᵀJ^m e₁, Q_K z^m = e₁ᵀJ_K^m e₁, J the Hermite Jacobi matrix):
   the ratios P_j(K) = D(2K + 2j)/D(2K) as exact polynomials in K (interpolated on K ≤ 2j + 1, checked on K ≤ P_KCHK), P_1 = K² + K + 1.
M. g_K's moments and cumulants exactly (fractions) for K in M_K: κ_{2j}·K^j against c_j = −(2j − 1)!/4^j, the coefficients of
   ½log(1 − η²/4) — the conjecture ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} ("conjecture P").
L. ĝ_K(ξ)e^{−ξ²/4} (D‴.58's `gen_ratio`) at |ξ| = s√K on both axes against (1 ∓ s²/4)^{1/2}, K in L_K.
S. The frame series S_K(ρ, v′) exactly (`frame_series`) against conjecture P's Laplace prediction — per angle the radial Laplace
   integral, then the angle by the trapezoid —, at the cusp with ₂F₁(−½, −½; 1; r²), off it with the factor 1 ∓ r.
R. D‴.58's scaled rules (A6's): each rule's own Mehler frame, the frame constant e_K² = (σ*/σ)(1 − ρ²)^{1/2}(v/τ²)·S_K(ρ, v′) checked
   against the exact L² error (`a6_error`), and the rules' exact errors (A6's raw) against the predicted sharp asymptotics.

D. The reduction (referee round 91): Duhamel's identity by quadrature, the Laplace step at s = ½ against the exact transform, the
   domination condition on a grid (`dom_gap`, `lam_plus`) and the radius (iv) needs per frame (`c_need`), and the corner ratio beyond
   Q's disc on the imaginary axis against the first alias's rate (`alias_rate`).

Usage: python experiments/exp_gh_peano.py [--only P,M,L,S,R,LC,Q,D,V,DB,VC,VA,VA011] [--from-raw] → results/chaos/gh_peano.md (+ _raw.npz). GP_WORKERS (default 2).
       python experiments/exp_gh_peano.py --export-cert [DIR] → D's and VA011's certified data in DIR (default results/certificates),
       for tools/verify_certificates.py (round 135, I4); writes nothing else.
"""
import os
import sys
import math
import cmath
import time
import warnings
from fractions import Fraction
from math import comb, factorial

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import exp_cohesion_leak as E
import exp_width_mismatch as WM

OUT = os.environ.get('GP_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GP_WORKERS', '2'))
P_J, P_KCHK = 6, 20                                                    # P: the ratios P_j for j ≤ P_J, checked on K ≤ P_KCHK
M_K = (10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120)                # M: the cumulants at these K
L_K, L_S = (20, 40, 80, 160), tuple(i / 10 for i in range(1, 20))      # L: K, and |ξ|/√K on each axis (below the branch point 2)
L_SMAX = 1.5                                                           # L: the deviation also summarised on s ≤ L_SMAX (away from the branch point)
S_K = (10, 20, 40, 80)                                                 # S: the frame series at these K
S_CHECK = 40                                                           # S: the last K recomputed with this many more digits
S_THETA = 4096                                                         # S: the trapezoid's points in the angle
R_K = (2, 3)                                                           # R: the frame constant checked at these K
R_TOL = 0.01                                                           # R: the rules within this of the prediction, summarised
R_PREC = 80                                                            # R: the digits of the frame constant's check
LC_K, LC_TH = (20, 40, 80, 160), (1, 2, 3)                            # LC (referee round 90): K, and θ = jπ/8 off the axes
Q_K, Q_MU, Q_TH = (20, 40, 80, 160), (0.1, 0.3, 0.5, 0.7, 0.9), (0, 2, 4)  # Q: K, |μ| and θ = jπ/8 of the corner ratio's μ
D_K, D_ETA, D_TH = (20, 40, 80, 160), (1.0, 1.5, 1.8, 1.9), (0, 1, 2, 3, 4)  # D (referee round 91): the Laplace step at K, |η|, arg η = jπ/8
D_QUAD = ((12, 1.5, 2), (12, 1.9, 1), (40, 1.5, 2), (40, 1.9, 4), (40, 1.95, 1))  # D: Duhamel's identity by quadrature at (K, |η|, j)
DOM_S, DOM_TH = 4001, 361                                               # D: the domination grid's points in s and in arg η ∈ [0, π/2]
DOM_MARGIN, DOM_NR = (1.7, 1.78, 1.85), 120                              # D: the gap's largest value on |η| ≤ these radii (DOM_NR radii from 1)
QX_Y = (1.2, 1.5, 1.8, 1.9, 2.0, 2.2, 2.5)                              # D: the corner ratio at u = i√K·y beyond Q's disc
QXP_K, QXP_N = (100, 150), 13                                           # D (referee round 93): Prony fits over QXP_N consecutive K from these
DOM_BB = (1.78, 1.8927)                                                 # D (round 93): the branch and bound on 1 < |η| ≤ these radii
DOM_SLACK = 1e-12                                                       # D: its floating-point slack
DOM_SECT = 20001                                                        # D (round 95): angles for the violated sector at 1.95
CERT_BANDS = ((1.8927, 1.98, 4, 32), (1.98, 1.99, 2, 32), (1.99, 1.995, 2, 32))   # V (round 99's recipe): the annulus bands, first grids
CERT_TB, CERT_EPS0, CERT_H, CERT_STEP = 0.9, 0.02, 2e-3, 10              # V: b at CERT_TB of the disc's edge; the end's radius; the flow's step, decimation
CERT_SLACK, CERT_LMIN, CERT_DEPTH = 1e-12, 1e-7, 20                     # V: the bound's slack, the shortest segment, the deepest box split
CERT_BALL_PREC = 120                                                    # V (referee round 104): the ball-arithmetic recheck's precision (bits)
DP_R, DP_TH, DP_H = (1.9, 1.95, 1.98), 19, (0.004, 0.002, 0.001)       # D: the deformed path's valley check (radii, angles, grid steps)


# ---------------------------------------------------------------------------------------------------------------- P: exact defects

def dyck(m, cap=None):
    """Σ over Dyck paths of length m with heights ≤ cap (None: unbounded) of Π over the up-steps (level + 1): e₁ᵀJ^m e₁ for the Hermite
    Jacobi matrix J (off-diagonals √1, √2, …) — each edge between the levels i − 1 and i traversed up and down contributes i —, and
    e₁ᵀJ_K^m e₁ (cap = K − 1) for its K × K section, i.e. E z^m and Q_K z^m for Z ~ N(0, 1)."""
    v = {0: 1}
    for _ in range(m):
        w = {}
        for h, c in v.items():
            if cap is None or h + 1 <= cap:
                w[h + 1] = w.get(h + 1, 0) + c * (h + 1)
            if h > 0:
                w[h - 1] = w.get(h - 1, 0) + c
        v = w
    return v.get(0, 0)


def defect(m, K):
    """D(m) = (E − Q_K)[z^m] exactly (an integer): the weighted Dyck paths of length m that reach height K."""
    return dyck(m) - dyck(m, K - 1)


def poly_fit(pts):
    """The exact polynomial through the points (x, y) (fractions), as its monomial coefficients (Gauss–Jordan on the Vandermonde system)."""
    n = len(pts)
    A = [[Fraction(x) ** p for p in range(n)] + [Fraction(y)] for x, y in pts]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        A[c] = [a / A[c][c] for a in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[c][n] for c in range(n)]


def poly_eval(coef, x):
    """Σ_p coef[p]·x^p, exactly for a Fraction or integer x."""
    return sum((c * Fraction(x) ** p for p, c in enumerate(coef)), Fraction(0))


def part_P():
    """P: for j ≤ P_J the ratio P_j(K) = D(2K + 2j)/D(2K) as an exact polynomial of degree 2j (interpolated on K = 1 … 2j + 1), checked
    against the path sums on K ≤ P_KCHK; D(2K) = K! checked there too."""
    vals = {}
    for K in range(1, P_KCHK + 1):
        d0 = defect(2 * K, K)
        vals[K] = dict(d0=d0, r=[Fraction(defect(2 * K + 2 * j, K), d0) for j in range(1, P_J + 1)])
    rows = []
    for j in range(1, P_J + 1):
        coef = poly_fit([(K, vals[K]['r'][j - 1]) for K in range(1, 2 * j + 2)])
        bad = [K for K in vals if poly_eval(coef, K) != vals[K]['r'][j - 1]]
        rows.append(dict(j=j, coef=[str(c) for c in coef], bad=bad, kchk=P_KCHK))
    return dict(rows=rows, d0_ok=all(vals[K]['d0'] == factorial(K) for K in vals), kchk=P_KCHK)


# ---------------------------------------------------------------------------------------------------------------- M: cumulants

def peano_moments(K, polys):
    """g_K's even moments M_2j = (2j)!(2K)!/(2K + 2j)!·P_j(K) exactly (its odd moments vanish), j ≤ len(polys)."""
    return {2 * j: Fraction(factorial(2 * j) * factorial(2 * K), factorial(2 * K + 2 * j)) * poly_eval(c, K) for j, c in enumerate(polys, 1)}


def cumulants(mom, n):
    """The cumulants κ_1 … κ_n from the moments (a dict, missing ones zero) by κ_n = m_n − Σ_{i<n} C(n − 1, i − 1)κ_i m_{n−i}, exactly."""
    m = {0: Fraction(1), **{i: Fraction(mom.get(i, 0)) for i in range(1, n + 1)}}
    k = {}
    for i in range(1, n + 1):
        k[i] = m[i] - sum((comb(i - 1, a - 1) * k[a] * m[i - a] for a in range(1, i)), Fraction(0))
    return k


def c_limit(j):
    """c_j = −(2j − 1)!/4^j: Σ_j c_j η^{2j}/(2j)! = ½log(1 − η²/4)."""
    return -Fraction(factorial(2 * j - 1), 4 ** j)


def part_M(P):
    """M: at K in M_K the cumulants' scaled values — K(κ_2 − ½) and κ_{2j}·K^j (2 ≤ j ≤ P_J), the coefficients of η^{2j}/(2j)! in
    log ĝ_K(√Kη) − Kη²/4 — against c_j, and K times their difference from c_j; the variance exactly."""
    polys = [[Fraction(c) for c in r['coef']] for r in P['rows']]
    out = []
    for K in M_K:
        k = cumulants(peano_moments(K, polys), 2 * len(polys))
        sc = [k[2 * j] * Fraction(K) ** j - (Fraction(K, 2) if j == 1 else 0) for j in range(1, len(polys) + 1)]
        out.append(dict(K=K, var=str(k[2]), scaled=[float(v) for v in sc], nxt=[float(K * (v - c_limit(j))) for j, v in enumerate(sc, 1)]))
    return out


# ---------------------------------------------------------------------------------------------------------------- L: the transform

def lt_limit(s, axis):
    """Conjecture P's limit of ĝ_K(ξ)e^{−ξ²/4} at |ξ| = s√K: (1 − s²/4)^{1/2} on the real axis, (1 + s²/4)^{1/2} on the imaginary."""
    return np.sqrt(1 - s * s / 4) if axis == 'real' else np.sqrt(1 + s * s / 4)


def l_job(args):
    """L: at K, axis, ĝ_K(s√K)e^{−Ks²/4} (D‴.58's `gen_ratio`, 3K + 60 digits) for s in L_S."""
    K, ax = args
    t0 = time.time()
    r = [WM.gen_ratio(K, s, ax, 3 * K + 60) for s in L_S]
    print(f"  L K {K} {ax}: " + ', '.join(f"{v / lt_limit(s, ax):.4f}" for v, s in zip(r, L_S)) + f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, axis=ax, s=list(L_S), ratio=r)


def part_L():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(l_job, [(K, ax) for K in L_K for ax in ('real', 'imag')], chunksize=1)


# ---------------------------------------------------------------------------------------------------------------- S: frame series

def gammas_ab(a, b):
    """The Fock integral's exponents in a frame's (a, b) = (ρv′, ρ(v′ − 1)), finite through ρ = 0 (referee round 90, D3ck-1): S_K =
    (πa)⁻¹∫|ε̃_K|²e^{−Ax²−By²}dA, A = (1 − b)/a, B = (1 + b)/a; with |ε̃_K|² ∝ |ξ|^{4K}e^{−(x² − y²)/2}|ĝ_K e^{−ξ²/4}|², the radial exponents
    γx = A + ½ (the real axis) and γy = B − ½ (the imaginary)."""
    return (1 - b) / a + 0.5, (1 + b) / a - 0.5


def gammas(rho, vp):
    """`gammas_ab` at a = ρv′, b = ρ(v′ − 1)."""
    return gammas_ab(rho * vp, rho * (vp - 1))


def log_s_pred_ab(K, a, b, n=S_THETA):
    """log S_K by conjecture P and Laplace's method in the radius at each angle θ (γ(θ) = γx cos²θ + γy sin²θ, the radius² at
    (2K + 1)/γ(θ), the factor |1 − ξ²/(4K)| there), the angle by the trapezoid on n points:
    log[(πa)⁻¹(K!)²/(2K)!²·½(2K)!·∫γ(θ)^{−2K−1}|1 − (2K + 1)e^{2iθ}/(4Kγ(θ))|dθ]."""
    from math import lgamma, log, pi
    gx, gy = gammas_ab(a, b)
    th = np.arange(n) * (2 * np.pi / n)
    g = gx * np.cos(th) ** 2 + gy * np.sin(th) ** 2
    lg = -(2 * K + 1) * np.log(g)
    f = np.abs(1 - (2 * K + 1) * np.exp(2j * th) / (4 * K * g))
    mx = lg.max()
    integ = mx + np.log(np.sum(np.exp(lg - mx) * f) * (2 * np.pi / n))
    return -log(pi * a) + 2 * lgamma(K + 1) - lgamma(2 * K + 1) - log(2) + integ


def log_s_pred(K, rho, vp, n=S_THETA):
    """`log_s_pred_ab` at a = ρv′, b = ρ(v′ − 1)."""
    return log_s_pred_ab(K, rho * vp, rho * (vp - 1), n)


def log_s_simple_ab(K, a, b):
    """The closed forms (conjecture P, Laplace's method to leading order): off the cusp, with γ = min(γx, γy), γ′ the other, r = 1/(2γ),
    S_K ≈ r^{2K}(1 ∓ r)/(a(2γ(γ′ − γ))^{1/2}) (− with the mass on the real axis); at the cusp S_K ≈ (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/a."""
    gx, gy = gammas_ab(a, b)
    real = gx <= gy
    g, g2 = (gx, gy) if real else (gy, gx)
    r = 1 / (2 * g)
    if abs(g2 - g) <= 1e-12 * g:                                        # the cusp (γx = γy up to rounding)
        return 0.5 * np.log(np.pi * K) + 2 * K * np.log(r) + np.log(hyp_half(r * r)) + np.log(2 * r / a)
    return 2 * K * np.log(r) + np.log(1 - r if real else 1 + r) - np.log(a * np.sqrt(2 * g * (g2 - g)))


def log_s_simple(K, rho, vp):
    """`log_s_simple_ab` at a = ρv′, b = ρ(v′ − 1)."""
    return log_s_simple_ab(K, rho * vp, rho * (vp - 1))


def hyp_half(x):
    """₂F₁(−½, −½; 1; x) = Σ_n binom(½, n)²xⁿ, |x| < 1."""
    tot, term, n = 1.0, 1.0, 0
    while abs(term) > 1e-17 * tot:
        term *= ((n - 0.5) / (n + 1)) ** 2 * x
        tot += term
        n += 1
    return tot


def frame_series_any(K, rho, a, prec, tol_exp=40, zero=None):
    """S_K = Σ_{n ≥ 2K even}(ρⁿ/n!)D_n(v′)² exactly in prec digits for any frame ρ ∈ (−1, 1), given (ρ, a = ρv′) as Decimals (referee
    round 90, D3ck-1): ρ > 0 D‴.58's `frame_series` at v′ = a/ρ; ρ < 0 (v′ < 0, √v′ imaginary) with D_n(v′)² = D̃_n(−v′)², D̃_n(u) the
    defect on H̃_n(√u z), H̃_{n+1}(x) = xH̃_n(x) + nH̃_{n−1}(x) (Σ H̃_n sⁿ/n! = e^{xs + s²/2}), E H̃_n(√uZ) = (u + 1)^{n/2}(n − 1)!!, and
    ρⁿ = |ρ|ⁿ at even n; ρ = 0 (zero, or |ρ| < 10^{−30}: the Fock scaling) the limit Σ_n (aⁿ/n!)D(n)², D(n) = (E − Q_K)zⁿ. Summed until the
    last two terms fall below 10^{−tol_exp} of the sum, one of them nonzero. Returns (S, the terms, the largest relative defect below 2K)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        rho, a = +Decimal(rho), +Decimal(a)
        zero = abs(rho) < Decimal(10) ** -30 if zero is None else zero
        if not zero and rho > 0:
            return WM.frame_series(K, rho, a / rho, prec=prec, tol_exp=tol_exp)
        zs, ws = WM.gh_nodes(K, prec)
        u = Decimal(0) if zero else -a / rho
        xs = [z if zero else u.sqrt() * z for z in zs]
        w_ = a if zero else -rho                                      # |ρ|; at ρ = 0 the terms carry aⁿ
        h_prev, h_cur = [Decimal(0)] * K, [Decimal(1)] * K
        S, terms, low = Decimal(0), [], Decimal(0)
        dfact, nfact, wn, tol, n = Decimal(1), Decimal(1), Decimal(1), Decimal(10) ** -tol_exp, 0
        while True:
            if n % 2 == 0:
                E_ = dfact if zero else (u + 1) ** (n // 2) * dfact
                Q = sum(w * h for w, h in zip(ws, h_cur))
                D = E_ - Q
                if n < 2 * K:
                    low = max(low, abs(D) / (abs(E_) + abs(Q) + 1))
                else:
                    terms.append(wn / nfact * D * D)
                    S += terms[-1]
                    if len(terms) > 10 and 0 < max(terms[-2:]) < tol * S:
                        return S, terms, low
            if zero:
                h_prev, h_cur = h_cur, [x_ * hc for x_, hc in zip(xs, h_cur)]
            else:
                h_prev, h_cur = h_cur, [x_ * hc + n * hp for x_, hc, hp in zip(xs, h_cur, h_prev)]
            n += 1
            nfact *= n
            wn *= w_
            if n % 2 == 0:
                dfact *= n - 1


def s_frames():
    """S's frames (ρ, v′) at τ/σ = 1/2, 1, 2 (σ = 1): the half-squeeze frame ρ = t at v′ = 2 (the cusp; the half-squeeze rule's own) and,
    synthetic — no rule's own frame —, at v′ = 1.6 and 2.5 (both branches), and the Mercer frame ρ = q, v′ = 1 + q (the Mercer frame
    rule's own); kept where the series converges (ρ(v′ − 1) < 1)."""
    out = []
    for r_ in (0.5, 1.0, 2.0):
        ss = np.sqrt(1 + r_ ** 2)
        t = (ss - 1) / (ss + 1)
        b = np.sqrt(1 + 2 * r_ ** 2)
        q = (b - 1) / (b + 1)
        for nm, rho, vp in (('half-squeeze', t, 2.0), ('half-squeeze', t, 1.6), ('half-squeeze', t, 2.5), ('Mercer', q, 1 + q)):
            if rho * (vp - 1) < 1:
                out.append((r_, nm, rho, vp))
    return out


def s_prec(K, rho, vp):
    """The digits for S_K: 2K·log₁₀(1/r) for the series' size, 2K for the defects' cancellation, and 60."""
    gx, gy = gammas(rho, vp)
    r = 1 / (2 * min(gx, gy))
    return int(np.ceil(2 * K * np.log10(1 / r))) + 2 * K + 60


def s_job(args):
    """S: at τ/σ = r_ in the frame (ρ, v′) the exact S_K (`frame_series` in s_prec digits) for K in S_K, the last also with S_CHECK more
    digits (referee round 90, D3ck-10: the digits the two Decimals agree to, −log₁₀ of their relative difference), against `log_s_pred`
    and `log_s_simple`."""
    from decimal import Decimal
    r_, nm, rho, vp = args
    t0 = time.time()
    rows = []
    for K in S_K:
        p = s_prec(K, rho, vp)
        S_, _, low = WM.frame_series(K, Decimal(repr(float(rho))), Decimal(repr(float(vp))), prec=p)
        rows.append(dict(K=K, prec=p, logS=float(S_.ln()), low=float(low), pred=log_s_pred(K, rho, vp), simple=log_s_simple(K, rho, vp)))
    K = S_K[-1]
    S2, _, _ = WM.frame_series(K, Decimal(repr(float(rho))), Decimal(repr(float(vp))), prec=rows[-1]['prec'] + S_CHECK)
    gx, gy = gammas(rho, vp)
    from decimal import localcontext
    with localcontext() as ctx:                                        # the digits the two agree to (D3ck-10), at the higher precision
        ctx.prec = rows[-1]['prec'] + S_CHECK
        d = abs(S2 / S_ - 1)
        chk = float(-d.log10()) if d != 0 else float(ctx.prec)
    out = dict(r=r_, frame=nm, rho=rho, vp=vp, gx=gx, gy=gy, rows=rows, check=chk, secs=time.time() - t0)
    print(f"  S tau/sigma {r_} {nm} v' {vp:g}: " + ', '.join(f"{y['K']}: {np.exp(y['logS'] - y['pred']):.5f}" for y in rows)
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_S():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(s_job, s_frames(), chunksize=1)


# ---------------------------------------------------------------------------------------------------------------- R: the rules

def rule_frame(r_, v):
    """The Mehler frame of D‴.58's rule of N(0, v) for the reweighted law at τ/σ = r_ (σ = 1): ρ = 1 + 2/τ² − 2/v, the frame whose own
    Gaussian N(0, v_ρ), v_ρ = 2τ²/(2 + (1 − ρ)τ²), is the rule's — it exists for −1 < ρ < 1, i.e. v ∈ (τ²/σ*², τ²); ρ < 0 below the Fock
    scaling, ρ = 0 at it (as a limit), ρ > 0 above (referee round 90, D3ck-1) —, ℓ² = 2ρ/(1 − ρ²), v′ = v/ℓ²; returned as (ρ, a, b, C)
    with a = ρv′ = v(1 − ρ²)/2, b = ρ(v′ − 1) = a − ρ and the constant C = (σ*/σ)(1 − ρ²)^{1/2}(v/τ²) of e_K²/‖f*‖² = C·S_K; (ρ, None,
    None, None) without a frame."""
    tau2 = r_ ** 2
    rho = 1 + 2 / tau2 - 2 / v
    if not -1 < rho < 1:
        return rho, None, None, None
    a = v * (1 - rho ** 2) / 2
    return rho, a, a - rho, np.sqrt(1 + tau2) * np.sqrt(1 - rho ** 2) * v / tau2


def r_rules(A6):
    """R's rules: A6's rules with a frame (−1 < ρ < 1) and a convergent series (A = (1 − b)/a > 0), and the others."""
    yes, no = [], []
    for x in (A6 if A6 is not None else []):
        rho, a, b, C = rule_frame(x['r'], x['v'])
        if a is not None and 1 - b > 0:
            yes.append((x, rho, a, b, C))
        else:
            no.append(dict(r=x['r'], name=x['name'], v_hs=x['v_hs'], rho=rho))
    return yes, no


def r_job(args):
    """R: the frame constant at K in R_K — e_K² (`a6_error`, R_PREC digits) over C·S_K (`frame_series_any`), minus 1 —, and A6's exact
    e_K² at A6_K over C·exp(`log_s_pred_ab`) (the per-angle prediction) and over C·exp(`log_s_simple_ab`) (the closed forms)."""
    from decimal import Decimal, localcontext
    x, rho, a, b, C = args
    chk = []
    for K in R_K:
        with localcontext() as ctx:
            ctx.prec = R_PREC
            name, v = next((nm, v) for nm, v in WM.a6_scalings(x['r']) if nm == x['name'])
            e2 = WM.a6_error(x['r'], v, K, R_PREC)
            tau2 = Decimal(repr(x['r'])) ** 2
            rh = 1 + 2 / tau2 - 2 / v
            ad = v * (1 - rh * rh) / 2
            Cd = (1 + tau2).sqrt() * (1 - rh * rh).sqrt() * v / tau2
            S_, _, _ = frame_series_any(K, rh, ad, R_PREC)
            chk.append(dict(K=K, dev=float(e2 / (Cd * S_) - 1)))
    le = lambda y: 2 * float(Decimal(y['e1']).ln()) - np.log(C)        # in logs: S_K underflows a double at small τ/σ
    rows = [dict(K=y['K'], ratio=float(np.exp(le(y) - log_s_pred_ab(y['K'], a, b))) if y['e1'] is not None else None,
                 closed=float(np.exp(le(y) - log_s_simple_ab(y['K'], a, b))) if y['e1'] is not None else None) for y in x['rows']]
    gx, gy = gammas_ab(a, b)
    return dict(r=x['r'], name=x['name'], v_hs=x['v_hs'], rho=rho, a=a, b=b, C=C, gx=gx, gy=gy, arg=x['arg'], rate=x['rate'], chk=chk,
                rows=rows)


def part_R(A6):
    from multiprocessing import Pool
    yes, no = r_rules(A6)
    with Pool(N_WORKERS) as pool:
        return dict(rules=pool.map(r_job, yes, chunksize=1), none=no)


# ---------------------------------------------------------------------------------------------------------------- LC: off the axes

def gen_ratio_c(K, s, j8, prec, pi=None, nodes=None):
    """ĝ_K(ξ)e^{−ξ²/4} = ε̃_K(ξ)/((K!/(2K)!)ξ^{2K}e^{−ξ²/4}) at ξ = s√K·e^{iθ}, θ = j8·π/8, in prec-digit decimal arithmetic by real and
    imaginary parts (referee round 90, D3ck-3): Q_K[e^{ξz}] = Σ_k w_k e^{xz_k}(cos yz_k + i sin yz_k), e^{−ξ²/2} = e^{−(x² − y²)/2}(cos xy −
    i sin xy). Returns the complex value as (re, im) floats."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        pi = WM.dpi() if pi is None else pi
        zs, ws = WM.gh_nodes(K, prec) if nodes is None else nodes
        th = pi * j8 / 8
        r = Decimal(repr(s)) * Decimal(K).sqrt()
        c, sn = WM.dcos(th, pi), WM.dcos(th - pi / 2, pi)
        x, y = r * c, r * sn
        qr = sum((w * (x * z).exp() * WM.dcos(y * z, pi) for w, z in zip(ws, zs)), Decimal(0))
        qi = sum((w * (x * z).exp() * WM.dcos(y * z - pi / 2, pi) for w, z in zip(ws, zs)), Decimal(0))
        m = (-(x * x - y * y) / 2).exp()
        er, ei = m * WM.dcos(x * y, pi), -m * WM.dcos(x * y - pi / 2, pi)
        nr, ni = 1 - (er * qr - ei * qi), -(er * qi + ei * qr)       # ε̃ = 1 − e^{−ξ²/2}Q
        mag = Decimal(factorial(K)) / factorial(2 * K) * r ** (2 * K) * (-(x * x - y * y) / 4).exp()
        ph = 2 * K * th - x * y / 2                                    # arg of ξ^{2K}e^{−ξ²/4}
        dr, di = mag * WM.dcos(ph, pi), mag * WM.dcos(ph - pi / 2, pi)
        den = dr * dr + di * di
        return float((nr * dr + ni * di) / den), float((ni * dr - nr * di) / den)


def lc_limit(s, j8):
    """Conjecture P's limit (1 − ξ²/(4K))^{1/2} at ξ = s√K·e^{ijπ/8} (the principal root; |ξ²/(4K)| = s²/4 < 1)."""
    import cmath
    return cmath.sqrt(1 - s * s / 4 * cmath.exp(2j * j8 * np.pi / 8))


def lc_job(args):
    """LC: at K and θ = j8·π/8, `gen_ratio_c` for s in L_S (prec 3K + 60 plus the digits of e^{y²/2}, y = Im ξ at the largest s), over
    conjecture P's limit: |ratio/limit − 1| and arg(ratio/limit)."""
    import math
    K, j8 = args
    t0 = time.time()
    y2 = (max(L_S) * math.sqrt(K) * math.sin(j8 * math.pi / 8)) ** 2
    prec = 3 * K + 60 + int(math.ceil(y2 / (2 * math.log(10))))
    from decimal import localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
    q = [complex(*gen_ratio_c(K, s, j8, prec, pi, nodes)) / lc_limit(s, j8) for s in L_S]
    print(f"  LC K {K} theta {j8}pi/8: " + ', '.join(f"{abs(v - 1):.1e}" for v in q) + f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, j8=j8, prec=prec, s=list(L_S), dev=[float(abs(v - 1)) for v in q], arg=[float(np.angle(v)) for v in q])


def part_LC():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(lc_job, [(K, j8) for K in LC_K for j8 in LC_TH], chunksize=1)


# ---------------------------------------------------------------------------------------------------------------- Q: the corner ratio

def corner_ratio(K, mu, j8, prec, pi=None, nodes=None, he=None):
    """c_K(u) = Q_K[e^{uz}He_{K−1}(z)]/(u^{K−1}e^{u²/2}) at u = μ√K·e^{ijπ/8} — the corner entry of e^{uJ_K} over the same for the infinite
    Jacobi matrix (E[e^{uZ}He_{K−1}(Z)] = u^{K−1}e^{u²/2}) —, in prec-digit decimal arithmetic by real and imaginary parts. Returns (re, im)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        pi = WM.dpi() if pi is None else pi
        zs, ws = WM.gh_nodes(K, prec) if nodes is None else nodes
        he = [WM.he_all(z, K)[K - 1] for z in zs] if he is None else he
        th = pi * j8 / 8
        r = Decimal(repr(mu)) * Decimal(K).sqrt()
        x, y = r * WM.dcos(th, pi), r * WM.dcos(th - pi / 2, pi)
        qr = sum((w * h * (x * z).exp() * WM.dcos(y * z, pi) for w, z, h in zip(ws, zs, he)), Decimal(0))
        qi = sum((w * h * (x * z).exp() * WM.dcos(y * z - pi / 2, pi) for w, z, h in zip(ws, zs, he)), Decimal(0))
        mag = r ** (K - 1) * ((x * x - y * y) / 2).exp()
        ph = (K - 1) * th + x * y                                      # arg of u^{K−1}e^{u²/2}
        dr, di = mag * WM.dcos(ph, pi), mag * WM.dcos(ph - pi / 2, pi)
        den = dr * dr + di * di
        return float((qr * dr + qi * di) / den), float((qi * dr - qr * di) / den)


def q_job(args):
    """Q: at K and θ = j8·π/8, `corner_ratio` for |μ| in Q_MU (3K + 80 digits) over its limit 1 − μ² ((iii)(a), proved in round 91): |c/(1 − μ²) − 1|."""
    import cmath
    from decimal import localcontext
    K, j8 = args
    t0 = time.time()
    prec = 3 * K + 80
    with localcontext() as ctx:
        ctx.prec = prec
        pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
        he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]
    dev = []
    for m in Q_MU:
        c = complex(*corner_ratio(K, m, j8, prec, pi, nodes, he))
        dev.append(float(abs(c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))
    print(f"  Q K {K} theta {j8}pi/8: " + ', '.join(f"{d:.1e}" for d in dev) + f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, j8=j8, prec=prec, mu=list(Q_MU), dev=dev)


def part_Q():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(q_job, [(K, j8) for K in Q_K for j8 in Q_TH], chunksize=1)


# ---------------------------------------------------------------------------------------------------------------- D: the reduction (round 91)

def d_prec(K, eta, j8):
    """The digits of D's jobs at K and η = |η|e^{ij8π/8}: 4K + 80 plus those of e^{y²/2}, y = Im ξ (as LC's)."""
    import math
    return 4 * K + 80 + int(math.ceil((eta * math.sqrt(K) * math.sin(j8 * math.pi / 8)) ** 2 / (2 * math.log(10))))


def d_job(args):
    """D, the Laplace step (referee round 91, D3cl-3): at K and arg η = j8·π/8, for |η| in D_ETA, ĝ_K(ξ)e^{−ξ²/4} (`gen_ratio_c`) and the
    corner ratio at ξ/2 (`corner_ratio`), ξ = √Kη; the step's claim is ĝ_K(ξ)e^{−ξ²/4} ~ (1 − η²/4)^{−1/2}c_K(ξ/2). Returns both, complex
    as (re, im)."""
    from decimal import localcontext
    K, j8 = args
    t0 = time.time()
    prec = d_prec(K, max(D_ETA), j8)
    with localcontext() as ctx:
        ctx.prec = prec
        pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
        he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]
    g = [gen_ratio_c(K, e, j8, prec, pi, nodes) for e in D_ETA]
    c = [corner_ratio(K, e / 2, j8, prec, pi, nodes, he) for e in D_ETA]
    st = [complex(*a) * lc_limit(e, j8) / complex(*b) for a, b, e in zip(g, c, D_ETA)]
    print(f"  D K {K} theta {j8}pi/8: " + ', '.join(f"{abs(v - 1):.1e}" for v in st) + f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, j8=j8, prec=prec, eta=list(D_ETA), g=[list(a) for a in g], c=[list(b) for b in c])


def duhamel_job(args):
    """D, the Duhamel identity by quadrature: ĝ_K(ξ)e^{−ξ²/4} = K·C(2K, K)∫_0^1 s^{K−1}(1 − s)^K e^{ξ²(1 − 2s)²/4}c_K(ξs)ds at
    ξ = √K|η|e^{ij8π/8}, Gauss–Legendre with 3K + 20 points in double, c_K by `corner_ratio` — against `gen_ratio_c`."""
    import cmath
    import math
    from decimal import localcontext
    K, eta, j8 = args
    t0 = time.time()
    prec = d_prec(K, eta, j8)
    with localcontext() as ctx:
        ctx.prec = prec
        pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
        he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]
    xs, ws = np.polynomial.legendre.leggauss(3 * K + 20)
    s, w = (xs + 1) / 2, ws / 2
    xi2 = K * (eta * cmath.exp(1j * j8 * math.pi / 8)) ** 2
    cs = [complex(*corner_ratio(K, eta * float(si), j8, prec, pi, nodes, he)) for si in s]
    lw = math.log(K) + math.lgamma(2 * K + 1) - 2 * math.lgamma(K + 1) + (K - 1) * np.log(s) + K * np.log(1 - s)
    quad = complex(np.sum(w * np.exp(lw + xi2 * (1 - 2 * s) ** 2 / 4) * np.array(cs)))
    exact = complex(*gen_ratio_c(K, eta, j8, prec, pi, nodes))
    print(f"  D quadrature K {K} eta {eta} theta {j8}pi/8: {abs(quad / exact - 1):.1e} [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, eta=eta, j8=j8, prec=prec, n=len(s), quad=[quad.real, quad.imag], exact=[exact.real, exact.imag])


def lam_plus(mu):
    """The rate of the bound on the corner ratio beyond Q's disc (referee round 91, D3cl-1): (1/K)log|c_K(√Kμ)| ≤ Λ⁺(μ) + o(1),
    Λ⁺(μ) = M(Re μ) − ½ − log|μ| − Re(μ²)/2, M(a) = a² for |a| ≤ 1 and 2|a| − 1 beyond — the contour integral on |ζ| = ρ₀ → 1, where
    max_θ Re[μ(ζ + 1/ζ)] − ½Re ζ² = max_θ (2a cos θ − ½cos 2θ) = M(a) + ½. Vectorised over μ."""
    a = np.abs(np.real(mu))
    return np.where(a <= 1, a * a, 2 * a - 1) - 0.5 - np.log(np.abs(mu)) - np.real(mu * mu) / 2


def dom_gap(eta, n=None):
    """The domination condition's gap at η (|η| < 2): sup over s ∈ [1/|η|, 1) — where u = ξs has left Q's disc — of g(s) + Λ⁺(ηs) − g(½),
    g(s) = log(s(1 − s)) + (Re η²/2)((1 − s)² + s²) the s-integrand's rate without c_K; negative: the integrand dominated by its value
    at s = ½ (−∞ when |η| ≤ 1: the segment never leaves the disc)."""
    n = DOM_S if n is None else n
    if abs(eta) <= 1:
        return -np.inf
    s = np.linspace(1 / abs(eta), 1, n)[:-1]
    e2 = (eta * eta).real
    g = np.log(s * (1 - s)) + e2 / 2 * ((1 - s) ** 2 + s ** 2)
    return float(np.max(g + lam_plus(eta * s)) - (np.log(0.25) + e2 / 4))


def dom_radius(th, lo=1.0, hi=1.9999, it=50):
    """The largest |η| < hi up to which the gap (`gap_closed`, its exact supremum over s; referee round 93) is negative along arg η = th
    (the first crossing on a grid of 400 radii, then bisection); hi if it never crosses."""
    rs = np.linspace(lo + 1e-9, hi, 400)
    bad = [r_ for r_ in rs if gap_closed(r_, th) >= 0]
    if not bad:
        return hi
    a, b = max([r_ for r_ in rs if r_ < bad[0]], default=lo), bad[0]
    for _ in range(it):
        m = (a + b) / 2
        a, b = (m, b) if gap_closed(m, th) < 0 else (a, m)
    return a


def c_need(r):
    """The radius (iv) needs P on for a frame with Laplace rate r (referee round 91): c² > 4 ln u/(u − 1), u = 1/r."""
    u = 1 / r
    return float(np.sqrt(4 * np.log(u) / (u - 1)))


def c_need_ball(r):
    """`c_need` at r's 4-decimal value in ball arithmetic (referee round 110, D3de-4: the thresholds are float bisections, so their floors'
    radii are rechecked rigorously): the enclosure's upper end, rounded up at 8 decimals."""
    import math
    arb, _ = _ball()
    u = 1 / (arb(int(round(r * 10 ** 4))) / 10 ** 4)
    return math.ceil(float((4 * u.log() / (u - 1)).sqrt().upper()) * 1e8) / 1e8


def r_star(R_dom, lo=0.01, hi=0.999, it=60):
    """The largest frame rate r whose needed radius `c_need` is below R_dom (c_need increases with r)."""
    for _ in range(it):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if c_need(m) < R_dom else (lo, m)
    return lo


def dom_sector(R=1.95, n=None):
    """The violated sector at |η| = R (referee round 99, D3ct-6): its ends by bisection on `gap_closed` (80 steps from DOM_SECT grid
    angles), and the range of the supremum's maximiser's Re(ηs) (`gap_arg`) over the grid angles inside and the two ends. None if the
    condition holds at every grid angle."""
    n = DOM_SECT if n is None else n
    ph = np.linspace(0, np.pi / 2, n)
    g = [gap_closed(R, float(f)) for f in ph]
    idx = [i for i, v in enumerate(g) if v > 0]
    if not idx:
        return None

    def cross(a, b):                                                    # gap(a) ≤ 0 < gap(b), or the reverse
        sa = gap_closed(R, a) > 0
        for _ in range(80):
            m = (a + b) / 2
            a, b = (m, b) if (gap_closed(R, m) > 0) == sa else (a, m)
        return (a + b) / 2
    lo = cross(float(ph[idx[0] - 1]), float(ph[idx[0]])) if idx[0] > 0 else 0.0
    hi = cross(float(ph[idx[-1]]), float(ph[idx[-1] + 1])) if idx[-1] < n - 1 else float(np.pi / 2)
    re = [gap_arg(R, f)[2] for f in [lo] + [float(ph[i]) for i in idx] + [hi]]
    return dict(R=R, n=n, lo=lo, hi=hi, re_min=min(re), re_max=max(re))


def part_dom():
    """D, the domination grid: `dom_radius` at DOM_TH angles in [0, π/2] (the gap is even in η and in its conjugate), the gap on both axes
    at |η| = 1.999, the angles where it is violated at |η| = 1.95 and 1.999, and the gap's largest value on the grid of |η| ≤ c (DOM_NR
    radii from 1 to c times the DOM_TH angles) for c in DOM_MARGIN — the margin the radius (iv) needs has."""
    ths = np.linspace(0, np.pi / 2, DOM_TH)
    rad = [dom_radius(t) for t in ths]
    viol = {}
    for r_ in (1.95, 1.999):
        v = [float(t) for t in ths if gap_closed(r_, t) >= 0]
        viol[str(r_)] = [min(v), max(v)] if v else None
    margin = {}
    for c in DOM_MARGIN:
        rs = np.linspace(1 + 1e-9, c, DOM_NR)
        margin[str(c)] = max((gap_closed(r_, t), float(t)) for r_ in rs for t in ths)       # (the gap, its angle)
    return dict(th=[float(t) for t in ths], radius=[float(x) for x in rad], axes=[gap_closed(1.999, 0.0), gap_closed(1.999, np.pi / 2)], viol=viol,
                margin=margin, sect=dom_sector())


def qx_job(K):
    """D, the corner ratio beyond Q's disc on the imaginary axis: c_K(i√Ky) (real there) for y in QX_Y, 4K + 80 digits."""
    from decimal import localcontext
    t0 = time.time()
    prec = 4 * K + 80
    with localcontext() as ctx:
        ctx.prec = prec
        pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
        he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]
    c = [corner_ratio(K, y, 4, prec, pi, nodes, he)[0] for y in QX_Y]
    print(f"  D beyond the disc K {K}: " + ', '.join(f"{v:.3g}" for v in c) + f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(K=K, prec=prec, y=list(QX_Y), c=c)


def part_D():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        jobs = pool.map(d_job, [(K, j8) for K in D_K for j8 in D_TH], chunksize=1)
        quad = pool.map(duhamel_job, list(D_QUAD), chunksize=1)
        qx = pool.map(qx_job, list(D_K), chunksize=1)
        qxp = pool.map(qxp_job, list(QXP_K), chunksize=1)
    dom = part_dom()
    dom['bb'] = [dict(R=R, n=n, worst=w, open=o, ball_worst=bw) for R in DOM_BB for n, w, o, bw in [dom_bb(R, ball=True)]]   # round 108
    return dict(jobs=jobs, quad=quad, qx=qx, qxp=qxp, dom=dom, dpath=part_dpath())


# ---------------------------------------------------------------------------------------------------------------- report

def p_lines(P):
    if P is None or len(P) == 0:
        return []
    fmt = lambda coef: ' + '.join((f"{c}" if p == 0 else f"{c}·K" if p == 1 else f"{c}·K^{p}") for p, c in enumerate(coef) if c != '0')
    out = ['', '## P. The defects as weighted paths: the ratios P_j(K) = D(2K + 2j)/D(2K)', '',
           "D(m) = (E − Q_K)[z^m] as the weighted Dyck paths of length m reaching height K (`defect`; D(2K) = K! on K ≤ "
           f"{P['kchk']}: {'yes' if P['d0_ok'] else 'NO'}), each P_j the exact polynomial through K = 1 … 2j + 1, checked on K ≤ {P['kchk']}:", '',
           '| j | P_j(K) | off the polynomial on K ≤ the check |', '|---|---|---|']
    for r in P['rows']:
        out.append(f"| {r['j']} | {fmt(r['coef'])} | {len(r['bad'])} |")
    return out


def m_lines(M):
    if M is None or len(M) == 0:
        return []
    J = len(M[0]['scaled'])
    out = ['', "## M. The normalised Peano kernel's cumulants (exact)", '',
           "g_K = G_K(2K)!/K! (a probability density): its variance exactly, (K² + K + 1)/((K + 1)(2K + 1)) by P_1, and κ_{2j}·K^j against "
           "c_j = −(2j − 1)!/4^j (the coefficients of ½log(1 − η²/4): " + ', '.join(f"{float(c_limit(j)):.6g}" for j in range(1, J + 1))
           + "), with K(κ_{2j}K^j − c_j) in brackets:", '',
           '| K | variance | K(κ_2 − ½) | ' + ' | '.join(f"κ_{2 * j}·K^{j}" for j in range(2, J + 1)) + ' |', '|---|---|' + '---|' * J]
    for x in M:
        out.append(f"| {x['K']} | {float(Fraction(x['var'])):.10f} | " + ' | '.join(f"{a:.6g} [{b:+.4g}]" for a, b in zip(x['scaled'], x['nxt'])) + ' |')
    a, b = (M[-2] if len(M) > 1 else None), M[-1]
    rich = [2 * y - x for x, y in zip(a['scaled'], b['scaled'])] if a is not None and b['K'] == 2 * a['K'] else None
    out += ['', f"At K = {b['K']} the scaled cumulants lie within " + ', '.join(
        f"{abs(v / float(c_limit(j)) - 1):.1e}" for j, v in enumerate(b['scaled'], 1)) + f" of c_1 … c_{J} (relative)"
        + ('' if rich is None else f"; extrapolated from K = {a['K']} and {b['K']} (2s(2K) − s(K), the 1/K term removed) "
           + ', '.join(f"{v:.6g}" for v in rich) + ', within ' + ', '.join(f"{abs(v / float(c_limit(j)) - 1):.1e}" for j, v in enumerate(rich, 1)))
        + '.']
    return out


def l_lines(L):
    if L is None or len(L) == 0:
        return []
    out = ['', "## L. The Laplace transform against conjecture P's limit", '',
           "ĝ_K(ξ)e^{−ξ²/4} (D‴.58's `gen_ratio`) over (1 − s²/4)^{1/2} on the real axis and (1 + s²/4)^{1/2} on the imaginary, |ξ| = s√K:", '',
           '| K | axis | ' + ' | '.join(f"s = {s:g}" for s in L[0]['s']) + ' | largest deviation |', '|---|---|' + '---|' * (len(L[0]['s']) + 1)]
    dev, dev2 = {}, {}
    for x in L:
        q = [v / lt_limit(s, x['axis']) for v, s in zip(x['ratio'], x['s'])]
        dev.setdefault(x['axis'], []).append((x['K'], max(abs(v - 1) for v in q)))
        dev2.setdefault(x['axis'], []).append((x['K'], max(abs(v - 1) for v, s in zip(q, x['s']) if s <= L_SMAX)))
        out.append(f"| {x['K']} | {x['axis']} | " + ' | '.join(f"{v:.4f}" for v in q) + f" | {max(abs(v - 1) for v in q):.2e} |")
    out += ['', "The largest deviation by K (its successive ratios): " + '; '.join(
        f"{ax}: " + ', '.join(f"{d:.2e}" for _, d in v) + ' (' + ', '.join(f"{v[i + 1][1] / v[i][1]:.3f}" for i in range(len(v) - 1)) + ')'
        for ax, v in dev.items()) + f"; on s ≤ {L_SMAX:g}: " + '; '.join(
        f"{ax}: " + ', '.join(f"{d:.2e}" for _, d in v) + ' (' + ', '.join(f"{v[i + 1][1] / v[i][1]:.3f}" for i in range(len(v) - 1)) + ')'
        for ax, v in dev2.items()) + '.']
    return out


def s_lines(S):
    if S is None or len(S) == 0:
        return []
    out = ['', "## S. The frame series against conjecture P's Laplace prediction", '',
           "S_K(ρ, v′) exactly (`frame_series`, 2K·log₁₀(1/r) + 2K + 60 digits) over the prediction — per angle the radius by Laplace's "
           "method with |1 − ξ²/(4K)|, the angle by the trapezoid (`log_s_pred`) —, S_K over the closed form (`log_s_simple`: off the cusp "
           "r^{2K}(1 ∓ r)/(ρv′(2γ(γ′ − γ))^{1/2}), at it (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/(ρv′)), r = 1/(2 min(γx, γy)), and the last K's "
           f"series against the same with {S_CHECK} more digits (the digits they agree to); the half-squeeze frame at v′ = 2 and the Mercer frame "
           "are the half-squeeze and the Mercer frame rule's own frames, v′ = 1.6 and 2.5 at ρ = t synthetic:", '',
           '| τ/σ | frame | ρ | v′ | axis (r) | over the prediction at K = ' + ', '.join(str(K) for K in S_K)
           + ' | K(ratio − 1) at the last two | over the closed form at the last K | digits agreed |', '|---|---|---|---|---|---|---|---|---|']
    by = {}
    for x in S:
        rr = [np.exp(y['logS'] - y['pred']) for y in x['rows']]
        ax = 'cusp' if abs(x['gx'] - x['gy']) <= 1e-12 * x['gx'] else ('real' if x['gx'] < x['gy'] else 'imaginary')
        by.setdefault(ax, []).append((rr[-1], x['rows'][-1]['K'] * (rr[-1] - 1), np.exp(x['rows'][-1]['logS'] - x['rows'][-1]['simple'])))
        r = 1 / (2 * min(x['gx'], x['gy']))
        Ks = [y['K'] for y in x['rows']]
        out.append(f"| {x['r']:g} | {x['frame']} | {x['rho']:.6f} | {x['vp']:g} | {ax} ({r:.6f}) | " + ', '.join(f"{v:.5f}" for v in rr)
                   + f" | {Ks[-2] * (rr[-2] - 1):+.3f}, {Ks[-1] * (rr[-1] - 1):+.3f} | {np.exp(x['rows'][-1]['logS'] - x['rows'][-1]['simple']):.5f} | "
                   f"{x['check']:.0f} |")
    out += ['', f"At K = {S[0]['rows'][-1]['K']}, by where the mass lies: " + '; '.join(
        f"{ax} ({len(v)} frames) {min(a for a, _, _ in v):.5f}–{max(a for a, _, _ in v):.5f}, K(ratio − 1) {min(b for _, b, _ in v):+.3f} … "
        f"{max(b for _, b, _ in v):+.3f}, over the closed form {min(c for _, _, c in v):.5f}–{max(c for _, _, c in v):.5f}"
        for ax, v in sorted(by.items())) + f"; the last K's series and the same with {S_CHECK} more digits agree to "
        f"{min(x['check'] for x in S):.0f}–{max(x['check'] for x in S):.0f} digits."]
    return out


def r_lines(R):
    if R is None or not isinstance(R, dict) or len(R) == 0:            # (a raw from before round 90 kept R as a list)
        return []
    rules, none = list(R['rules']), list(R['none'])
    if not rules:
        return []
    Ks = [y['K'] for y in rules[0]['rows']]
    sgn = lambda x: 'ρ > 0' if x['rho'] > 1e-12 else ('ρ < 0' if x['rho'] < -1e-12 else 'ρ = 0')
    cnt = {}
    for x in rules:
        cnt[sgn(x)] = cnt.get(sgn(x), 0) + 1
    out = ['', "## R. D‴.58's scaled rules: the frame constant and the sharp asymptotics", '',
           "Each rule of A6 with its own Mehler frame (−1 < ρ < 1: ρ = 1 + 2σ²/τ² − 2σ²/v, ℓ² = 2σ²ρ/(1 − ρ²), v′ = v/ℓ², a = ρv′ = v(1 − ρ²)/(2σ²), "
           "b = ρ(v′ − 1); ρ = 0, the Fock scaling, as the limit): the exact squared L² error (`a6_error`, "
           f"{R_PREC} digits) over (σ*/σ)(1 − ρ²)^{{1/2}}(v/τ²)·S_K (`frame_series_any`), minus 1, at K = " + ', '.join(str(K) for K in R_K)
           + ", and A6's exact e_K² over the constant times conjecture P's prediction of S_K — per angle (`log_s_pred_ab`) at K = "
           + ', '.join(str(K) for K in Ks) + ", and the closed forms (`log_s_simple_ab`) at the last —, by where the mass lies:", '',
           '| τ/σ | rule | v/v_hs | ρ | a | axis | the constant: deviation | e_K² over the per-angle prediction | over the closed form |',
           '|---|---|---|---|---|---|---|---|---|']
    for x in rules:
        ax = 'cusp' if abs(x['gx'] - x['gy']) <= 1e-12 * x['gx'] else ('real' if x['gx'] < x['gy'] else 'imaginary')
        out.append(f"| {x['r']:g} | {x['name']} | {x['v_hs']:.6f} | {x['rho']:+.6f} | {x['a']:.6f} | {ax} | "
                   + ', '.join(f"{y['dev']:+.1e}" for y in x['chk']) + ' | '
                   + ', '.join('—' if y['ratio'] is None else f"{y['ratio']:.5f}" for y in x['rows']) + ' | '
                   + ('—' if x['rows'][-1]['closed'] is None else f"{x['rows'][-1]['closed']:.5f}") + ' |')
    dv = max(abs(y['dev']) for x in rules for y in x['chk'])
    summ = lambda key: [(x, x['rows'][-1][key]) for x in rules if x['rows'][-1][key] is not None]
    parts = []
    for key, lab in (('ratio', 'the per-angle prediction'), ('closed', 'the closed forms')):
        last = summ(key)
        near = [v for _, v in last if abs(v - 1) <= R_TOL]
        far = [(x, v) for x, v in last if abs(v - 1) > R_TOL]
        parts.append(f"over {lab} {min(v for _, v in last):.5f}–{max(v for _, v in last):.5f}, {len(near)} of the {len(last)} within "
                     f"{R_TOL:g}" + (f" ({min(near):.5f}–{max(near):.5f})" if near else '')
                     + ('; the others ' + ', '.join(f"τ/σ = {x['r']:g}, {x['name']} ({v:.5f})" for x, v in far) if far else ''))
    out += ['', f"Of A6's {len(rules) + len(none)} rules {len(rules)} have a frame (" + ', '.join(f"{v} with {k}" for k, v in sorted(cnt.items()))
            + f") and {len(none)} not (|ρ| ≥ 1: " + ', '.join(f"τ/σ = {x['r']:g}, {x['name']} ({x['rho']:+.3f})" for x in none)
            + f"). The frame constant holds to {dv:.1e} at all {len(rules)} — the {R_PREC}-digit rounding of a cancelling sum —; at K = {Ks[-1]} "
            "the exact e_K² lie " + '; '.join(parts) + '.']
    return out


def lc_lines(LC):
    if LC is None or len(LC) == 0:
        return []
    out = ['', "## LC. The Laplace transform off the axes (referee round 90)", '',
           "ĝ_K(ξ)e^{−ξ²/4} at ξ = s√K·e^{iθ} (`gen_ratio_c`, complex decimal arithmetic) over conjecture P's (1 − ξ²/(4K))^{1/2}: |ratio/limit − 1|:", '',
           '| K | θ | ' + ' | '.join(f"s = {s:g}" for s in LC[0]['s']) + ' | largest |', '|---|---|' + '---|' * (len(LC[0]['s']) + 1)]
    dev = {}
    for x in LC:
        dev.setdefault(x['j8'], []).append(max(x['dev']))
        out.append(f"| {x['K']} | {x['j8']}π/8 | " + ' | '.join(f"{v:.1e}" for v in x['dev']) + f" | {max(x['dev']):.1e} |")
    out += ['', "The largest deviation by K (its successive ratios): " + '; '.join(
        f"θ = {j}π/8: " + ', '.join(f"{d:.2e}" for d in v) + ' (' + ', '.join(f"{v[i + 1] / v[i]:.3f}" for i in range(len(v) - 1)) + ')'
        for j, v in sorted(dev.items())) + f"; the largest |arg| {max(abs(a) for x in LC for a in x['arg']):.1e}."]
    return out


def q_lines(Q):
    if Q is None or len(Q) == 0:
        return []
    out = ['', "## Q. The corner ratio: its limit, (iii)(a)", '',
           "c_K(u) = Q_K[e^{uz}He_{K−1}(z)]/(u^{K−1}e^{u²/2}) at u = μ√K (`corner_ratio`: the corner entry of e^{uJ_K} over the infinite Jacobi "
           "matrix's, complex decimal arithmetic) against its limit 1 − μ² (|μ| < 1; (iii)(a)): |c/(1 − μ²) − 1|:", '',
           '| K | arg μ | ' + ' | '.join(f"\\|μ\\| = {m:g}" for m in Q[0]['mu']) + ' |', '|---|---|' + '---|' * len(Q[0]['mu'])]
    by = {}
    for x in Q:
        by.setdefault(x['j8'], []).append(max(x['dev']))
        out.append(f"| {x['K']} | {x['j8']}π/8 | " + ' | '.join(f"{d:.2e}" for d in x['dev']) + ' |')
    out += ['', "The largest deviation by K (its successive ratios): " + '; '.join(
        f"arg μ = {j}π/8: " + ', '.join(f"{d:.2e}" for d in v) + ' (' + ', '.join(f"{v[i + 1] / v[i]:.3f}" for i in range(len(v) - 1)) + ')'
        for j, v in sorted(by.items())) + '.']
    return out


def alias_rate(y):
    """The growth rate of c_K(i√Ky) beyond Q's disc (referee round 93, D3cn-1): the first alias — the second saddle t = i(1/ζ − ζ) of
    I_K's integrand —, whose exponent f₁(ζ) = σ(ζ + 1/ζ) + 3 log ζ − ζ² + ζ⁻²/2 (σ = iy) is stationary at ζ = (iy ± (8 − y²)^{1/2})/4, on
    |ζ| = 1/√2: Re f₁ − Re(1 + σ²/2 + log σ) = (5/8)y² − log y − ½ − (3/2)log 2. (Round 91's (11/18)y² − log y − 3/2 took the modulus
    at a real stationary-phase point under a modulus that varies on the same exponential scale: wrong.)"""
    import math
    return 5 / 8 * y * y - math.log(y) - 0.5 - 1.5 * math.log(2)


def prony(c, lam0=0.0):
    """The dominant root of a two-term linear recurrence fitted to c_k (least squares, c_{k+2} = a c_{k+1} + b c_k) after scaling by
    e^{−k·lam0}: (its log-modulus + lam0, its |argument|) — a growth rate and a frequency per step."""
    import math
    d = np.asarray(c, float) * np.exp(-np.arange(len(c)) * lam0)
    a, b = np.linalg.lstsq(np.c_[d[1:-1], d[:-2]], d[2:], rcond=None)[0]
    z = np.roots([1.0, -a, -b])
    z = z[np.argmax(np.abs(z))]
    return float(math.log(abs(z)) + lam0), float(abs(np.angle(z)))


def qxp_job(K0):
    """D (round 93): c_K(i√Ky) − (1 + y²) at QXP_N consecutive K from K0 for the y of QX_Y beyond `alias_rate`'s zero (4K + 80 digits)."""
    from decimal import localcontext
    t0 = time.time()
    ys = [y for y in QX_Y if alias_rate(y) > 0]
    rows = []
    for K in range(K0, K0 + QXP_N):
        prec = 4 * K + 80
        with localcontext() as ctx:
            ctx.prec = prec
            pi, nodes = WM.dpi(), WM.gh_nodes(K, prec)
            he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]
        rows.append(dict(K=K, prec=prec, c=[corner_ratio(K, y, 4, prec, pi, nodes, he)[0] - (1 + y * y) for y in ys]))
    print(f"  D Prony K {K0}: [{time.time() - t0:.0f}s]", flush=True)
    return dict(K0=K0, y=ys, rows=rows)


def gap_closed(R, phi):
    """The domination gap in closed form (referee round 93, D3cn-3): with w = 1 − s ∈ (0, 1 − 1/R], p = R²cos 2φ, the gap is
    log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4 (M as in `lam_plus`); its supremum over w — the endpoint 1 − 1/R, the
    kink R(1 − w)cos φ = 1, and the roots of the quadratics h′ = 0 on each piece. −∞ for R ≤ 1."""
    import math
    if R <= 1:
        return -math.inf
    return bound_sup(R, R, R * R * math.cos(2 * phi), R * math.cos(phi))


def quad_roots(qa, qb, qc):
    """The real roots of qa·w² + qb·w + qc by the stable formula (referee round 95, D3cp-2: an absolute zero test on qa missed qa =
    R²cos(π/2) ≈ 6·10⁻¹⁷ and the cancelling root lost the critical point): q = −(qb + sign(qb)√disc)/2, the roots q/qa and qc/q; qa
    negligible against qb, qc (relative 10⁻¹⁴): the linear root only."""
    import math
    if abs(qa) <= 1e-14 * max(abs(qb), abs(qc)):
        return [-qc / qb] if qb else []
    disc = qb * qb - 4 * qa * qc
    if disc < 0:
        return []
    q = -(qb + math.copysign(math.sqrt(disc), qb)) / 2
    return [q / qa] + ([qc / q] if q != 0 else [])


def bound_sup(R1, R2, p, B, arg=False):
    """sup over w ∈ (0, 1 − 1/R2] of log w + (p/2)(w² − ½) + M(B(1 − w)) − ½ − log R1 + log 4 (`gap_closed` at R1 = R2, p = R²cos 2φ,
    B = R cos φ; on a box an upper bound with R1 the least radius, p the least of R²cos 2φ, B the largest R cos φ: w² < ½ and M increase).
    arg: (the supremum, its w)."""
    import math
    wm = 1 - 1 / R2
    Mf = lambda a: a * a if abs(a) <= 1 else 2 * abs(a) - 1
    h = lambda w: math.log(w) + p / 2 * (w * w - 0.5) + Mf(B * (1 - w)) - 0.5 - math.log(R1) + math.log(4)
    cand = [wm]
    if B > 1 and 0 < 1 - 1 / B < wm:
        cand.append(1 - 1 / B)
    for qa, qb, qc in ((p + 2 * B * B, -2 * B * B, 1.0), (p, -2 * B, 1.0)):   # h′ = 0 on the pieces a ≤ 1 and a > 1
        cand += [w for w in quad_roots(qa, qb, qc) if 0 < w < wm]
    return max((h(w), w) for w in cand) if arg else max(h(w) for w in cand)


def gap_arg(R, phi):
    """The maximiser of `gap_closed` over s ∈ [1/R, 1): (the gap, s, Re(ηs)) — where the straight segment's integrand peaks beyond the
    disc (referee round 95, D3cp-5: the draft's range of Re(ηs), now printed)."""
    g, w = bound_sup(R, R, R * R * math.cos(2 * phi), R * math.cos(phi), arg=True)
    return g, 1 - w, R * (1 - w) * math.cos(phi)


def box_bound(R1, R2, f1, f2):
    """An upper bound of the domination gap over the box R1 ≤ |η| ≤ R2, f1 ≤ arg η ≤ f2 ⊂ [0, π/2] (referee round 93): `bound_sup` with the
    least radius, the least of R²cos 2φ (at a corner: it is monotone in φ and, at fixed sign, in R) and the largest R cos φ (at R2, f1)."""
    import math
    if R2 <= 1:
        return -math.inf
    p_lo = min(R * R * math.cos(2 * f) for R in (R1, R2) for f in (f1, f2))
    return bound_sup(max(R1, 1.0 + 1e-15), R2, p_lo, R2 * math.cos(f1))


def box_bound_ball(R1, R2, f1, f2):
    """`box_bound` in ball arithmetic (referee round 108, D3dc-1, after the referee's port): with R1 exact (no 1 + 10⁻¹⁵), p the rigorous
    least of R²cos 2φ over the corners (the top edge at the true π/2 where it is the float π/2), B = R2 cos f1 and wm = 1 − 1/R2 at
    their upper ends (w² < ½ checked), h(w) = log w + (p/2)(w² − ½) + M(B(1 − w)) − ½ − log R1 + log 4 bounds the gap of every η of the
    box on (0, wm]; its supremum is at wm, the kink B(1 − w) = 1 or a root of a piece's h′ = 0 — (p + 2B²)w² − 2B²w + 1 and pw² − 2Bw
    + 1 —, each enclosed by a ball (the stable formula; a root that cannot be placed, or a candidate ball that may reach 0, returns
    None). Returns the bound's rigorous upper end as a float rounded up, or None."""
    arb, acb = _ball()
    half, log4 = arb(1) / 2, arb(4).log()
    R1a, R2a, f1a = arb(R1), arb(R2), arb(f1)
    f2a = arb.pi() / 2 if f2 == math.pi / 2 else arb(f2)
    p = min((R * R * (2 * f).cos()).lower() for R in (R1a, R2a) for f in (f1a, f2a))
    B = (R2a * f1a.cos()).upper()
    wm = (1 - 1 / R2a).upper()
    if not (wm * wm < half):
        return None
    cand = [wm]
    if B > 1:
        wk = 1 - 1 / B
        if not (wk.upper() <= 0) and not (wk.lower() > wm):
            cand.append(wk)
    for qa, qb in ((p + 2 * B * B, -2 * B * B), (p, -2 * B)):
        if not (qb < 0):
            return None
        disc = qb * qb - 4 * qa
        if disc < 0:
            continue
        if not (disc >= 0):
            return None
        q = (disc.sqrt() - qb) / 2                                      # > 0: the roots 1/q and q/qa
        rs = [1 / q]
        if qa > 0 or qa < 0:
            rs.append(q / qa)
        elif not (q.lower() > abs(qa).upper() * wm):
            return None
        cand += [r for r in rs if not (r.upper() <= 0) and not (r.lower() > wm)]
    best = None
    for w in cand:
        if not (w > 0):
            return None
        a_up = (B * (1 - w)).upper()
        Mu = a_up * a_up if a_up <= 1 else 2 * a_up - 1
        v = (w.log() + p / 2 * (w * w - half) + Mu - half - R1a.log() + log4).upper()
        best = v if best is None or v > best else best
    return math.nextafter(float(best), math.inf)


def dom_bb(Rmax, slack=None, max_boxes=3_000_000, ball=False, keep=None):
    """The domination condition on 1 < |η| ≤ Rmax, 0 ≤ arg η ≤ π/2, by branch and bound on `bound_sup`'s box bound: boxes whose bound is
    below −slack are verified, the others bisected (floating point with slack). With ball (referee round 108) every accepted box is
    rechecked by `box_bound_ball` and bisected further unless that bound is below 0. Returns (verified boxes, the largest bound among
    them, unresolved boxes) and, with ball, the largest rigorous bound. With keep (a list; round 135, I4: the certificate's export, not
    used otherwise), every accepted leaf (R1, R2, f1, f2, its rigorous bound or None) is appended in the order accepted."""
    import math
    slack = DOM_SLACK if slack is None else slack
    stack, n, worst, bw = [(1.0, Rmax, 0.0, math.pi / 2)], 0, -math.inf, -math.inf
    while stack and n < max_boxes:
        R1, R2, f1, f2 = stack.pop()
        b = box_bound(R1, R2, f1, f2)
        if b < -slack:
            u = box_bound_ball(R1, R2, f1, f2) if ball else None
            if not ball or (u is not None and u < 0):
                n += 1
                worst = max(worst, b)
                bw = max(bw, u) if ball else bw
                if keep is not None:
                    keep.append((R1, R2, f1, f2, u))
                continue
        if R2 - R1 > f2 - f1:
            m = (R1 + R2) / 2
            stack += [(R1, m, f1, f2), (m, R2, f1, f2)]
        else:
            m = (f1 + f2) / 2
            stack += [(R1, R2, f1, m), (R1, R2, m, f2)]
        if R2 - R1 < 1e-9 and f2 - f1 < 1e-9:
            return (n, worst, len(stack)) + ((bw,) if ball else ())
    return (n, worst, len(stack)) + ((bw,) if ball else ())


def path_ok(eta, h):
    """D (round 93): the deformed Duhamel path's valley condition at η (|η| < 2) on a grid of step h — with the bound F(σ) = log|η − σ| +
    Re(η − σ)²/2 + log|σ| + Re σ²/2 + C(σ) (C = 0 inside the unit disc, Λ⁺ outside) and F0 = F(η/2), the seeds η/2 ∓ PATH_SEED·rb·u on
    the saddle's steepest-descent line (`path_seeds`; round 95, formerly on the segment), rb = ½(1 − |η|/2), lie in the components of
    {F < F0} ∖ B(η/2, rb) (8-connected) that contain 0 and η, and those two components differ. Then a path 0 → η/2 → η runs in the sublevel set but for the segment inside the ball, where (a) holds."""
    return bool(path_check(eta, h)[0])


PATH_SEED, PATH_TOL = 1.2, 0.0                                          # round 95: the seeds at 1.2·rb from η/2; the sublevel set strict


def path_seeds(eta):
    """The deformed path's seeds (referee round 95, D3cp-1): η/2 ∓ PATH_SEED·rb·u on the steepest-descent line of the saddle — u² ∝
    −conj(φ″), φ″ = 2 − 8/η² the exponent's second derivative at η/2, the sign making the first seed the one toward 0 —, rb = ½(1 − |η|/2)."""
    import cmath
    rb = 0.5 * (1 - abs(eta) / 2)
    ph = 2 - 8 / eta ** 2
    u = cmath.exp(1j * (math.pi - cmath.phase(ph)) / 2)
    if (u * eta.conjugate()).real < 0:
        u = -u
    return eta / 2 - PATH_SEED * rb * u, eta / 2 + PATH_SEED * rb * u, rb


def path_F(eta, S):
    """The valley check's bound F(σ) = log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + C(σ), C = 0 for |σ| < 1 and Λ⁺(σ) (`lam_plus`)
    beyond, on an array (or a number) σ = S: the grid's and the seeds' values from one formula (round 95)."""
    S = np.asarray(S, dtype=complex)
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = np.where(np.abs(S) < 1, 0.0, lam_plus(S))
        return np.log(np.abs(eta - S)) + ((eta - S) ** 2).real / 2 + np.log(np.abs(S)) + (S * S).real / 2 + lam


def chord_q(eta, D):
    """The chord's Taylor ratio (referee round 99, D3ct-1): inside the unit disc F = Re φ, φ(σ) = log(η − σ) + (η − σ)²/2 + log σ + σ²/2,
    with φ′(η/2) = φ‴(η/2) = 0, φ″(η/2) = 2 − 8/η², φ⁗ = −6/(η − σ)⁴ − 6/σ⁴, so on the chord η/2 + tu, |t| ≤ D, u on the steepest-descent
    line (Re φ″u² = −|φ″|), F − F(η/2) ≤ −|φ″|t²/2 + sup|φ⁗|t⁴/24 ≤ −(|φ″|t²/2)(1 − q(t/D)²), sup|φ⁗| ≤ 12/(|η|/2 − D)⁴ on the ball:
    q = D²·12/(|η|/2 − D)⁴/(12|φ″|) = D²/((|η|/2 − D)⁴|φ″|). q < 1: F < F(η/2) on the whole chord but at η/2 (Taylor's theorem)."""
    eta = complex(eta)
    return D * D / ((abs(eta) / 2 - D) ** 4 * abs(2 - 8 / eta ** 2))


def path_check(eta, h):
    """`path_ok`'s check and the seeds' margins F(seed) − F(η/2): (ok, margin toward 0, margin toward η, the chord's Taylor ratio
    `chord_q`). Requires h/√2 < 0.2·rb and the seeds' rounded pixels outside the ball (referee rounds 95, 99). The chord a → b through
    η/2, inside the unit disc (1.2rb < 2rb), stays below F(η/2) but at η/2 where q < 1 (round 99: Taylor's bound in place of round 95's
    sampling)."""
    from scipy import ndimage
    eta = complex(eta)
    a_, b_, rb = path_seeds(eta)
    if not h / math.sqrt(2) < 0.2 * rb:
        raise ValueError(f"the grid step {h} is too coarse for the ball's radius {rb}")
    x0, x1 = min(0, eta.real) - 1.2, max(0, eta.real) + 1.2
    y0, y1 = min(0, eta.imag) - 1.2, max(0, eta.imag) + 1.2
    X, Y = np.meshgrid(np.arange(x0, x1, h), np.arange(y0, y1, h))
    S = X + 1j * Y
    F = path_F(eta, S)
    F0 = 2 * np.log(abs(eta) / 2) + (eta * eta).real / 4
    lab, _ = ndimage.label((F < F0 - 1e-12 + PATH_TOL) & (np.abs(S - eta / 2) > rb), structure=np.ones((3, 3)))
    ix = lambda z: (int(round((z.imag - y0) / h)), int(round((z.real - x0) / h)))
    at = lambda z: lab[ix(z)]
    for z in (a_, b_):                                                  # round 99 (D3ct-4, K25): the seed's own pixel outside the ball
        i, j = ix(z)
        if not abs(complex(X[i, j], Y[i, j]) - eta / 2) > rb:
            raise ValueError(f"a seed's pixel lies in the ball at {eta}")
    l0, le, la, lb = at(0.02 * eta), at(0.98 * eta), at(a_), at(b_)
    fv = lambda z: float(path_F(eta, z)) - F0
    q = chord_q(eta, abs(b_ - eta / 2))
    return bool(l0 and le and l0 == la and le == lb and l0 != le and q < 1), fv(a_), fv(b_), q


def cert_E(z, eta):
    """The valley bound in the saddle's coordinates (round 99's centred form): with σ = (η/2)(1 + z), F_η(σ) − F_η(η/2) = log|1 − z²| +
    Re(cz²) + C(σ), c = η²/4, C = 0 for |σ| ≤ 1 and Λ⁺(σ) = h(|σ|) − (|Re σ| − 1)₊² beyond, h(ρ) = (ρ² − 1)/2 − log ρ (exact: the log and
    quadratic terms of F and F(η/2) combine)."""
    c = eta * eta / 4
    sg = eta / 2 * (1 + z)
    rho = abs(sg)
    C = 0.0 if rho <= 1 else (rho * rho - 1) / 2 - math.log(rho) - max(abs(sg.real) - 1, 0.0) ** 2
    return math.log(abs(1 - z * z)) + (c * z * z).real + C


def cert_grad(z, eta):
    """∇_z of `cert_E` as a complex number (∂x + i∂y): conj(g′) of g = log(1 − z²) + cz², and C's gradient by the chain rule."""
    c = eta * eta / 4
    v = (2 * z * (c - 1 / (1 - z * z))).conjugate()
    sg = eta / 2 * (1 + z)
    rho = abs(sg)
    if rho > 1:
        a = sg.real
        gs = (rho - 1 / rho) * sg / rho - (2 * (abs(a) - 1) * math.copysign(1, a) if abs(a) > 1 else 0.0)
        v += gs * (eta / 2).conjugate()
    return v


def cert_flow(eta, z0):
    """A path from z0 (b) toward z = 1 (σ = η): normalised steepest descent of `cert_E` with step CERT_H until |1 − z| < CERT_EPS0/2,
    decimated to every CERT_STEP-th point (None if it does not arrive in 2·10⁴ steps)."""
    z, pts = complex(z0), [complex(z0)]
    for _ in range(20000):
        if abs(1 - z) < CERT_EPS0 / 2:
            return pts[::CERT_STEP] + ([pts[-1]] if (len(pts) - 1) % CERT_STEP else [])
        g = cert_grad(z, eta)
        z = z - CERT_H * g / abs(g)
        pts.append(z)
    return None


def cert_seg_ub(p, q, eta_c, delta, R2):
    """`cert_seg_parts`' bound E(m; η_c) + ℓ·L_z + δ·L_η (∞ where the segment may reach z = ±1)."""
    Em, l, Lz, Leta = cert_seg_parts(p, q, eta_c, R2)
    return Em + l * Lz + delta * Leta


def cert_seg_parts(p, q, eta_c, R2):
    """An upper bound of `cert_E` over the segment [p, q] × the η-disc |η − η_c| ≤ δ, |η| ≤ R2: E(m; η_c) + ℓ·L_z + δ·L_η, m the
    midpoint, ℓ the half-length; L_z ≥ sup over the segment of |∇_z E(·; η_c)| — |g′(m)| + ℓ·sup|g″|, g″ = 2c − 2(1 + z²)/(1 − z²)²
    bounded with |z| ≤ |m| + ℓ and |1 − z²| ≥ |1 − m²| − ℓ(2|m| + ℓ), and C's |∇_σC| ≤ h′(|σ|) (h′ increasing, 0 inside the disc;
    sharp where |Re σ| ≤ 1 < |σ|) times |dσ/dz| = |η|/2 —; L_η ≥ sup of |∂_ηE| — |η||z|²/2 from Re(cz²) and |1 + z|/2 times C's gradient bound at |σ| ≤
    (R2/2)|1 + z|. Infinite where the segment may reach z = ±1. Floating point; `cert_seg_ball` rechecks it in ball arithmetic.
    Referee round 104 (D3cy-4): |∇_σC| ≤ h′(|σ|) exactly — beyond |Re σ| = 1, with a = |Re σ|, ρ = |σ|, |∇C|² = h′(ρ)² − 4(a − 1)(1 −
    a/ρ²) ≤ h′(ρ)² since a ≤ ρ —, so the 2(|Re σ| − 1)₊ terms of round 99's bound are dropped."""
    m, l = (p + q) / 2, abs(q - p) / 2
    c = eta_c * eta_c / 4
    Zmax = abs(m) + l
    W = abs(1 - m * m) - l * (2 * abs(m) + l)
    if W <= 0:
        return math.inf, l, math.inf, math.inf
    Lg = abs(2 * m * (c - 1 / (1 - m * m))) + l * (2 * abs(c) + 2 * (1 + Zmax ** 2) / W ** 2)
    Rc = abs(eta_c)
    sm = eta_c / 2 * (1 + m)
    hp = lambda r: (r - 1 / r) if r > 1 else 0.0
    rmax = abs(sm) + Rc / 2 * l
    Lz = Lg + Rc / 2 * hp(rmax)
    opz = abs(1 + m) + l
    rB = R2 / 2 * opz
    Leta = R2 * Zmax ** 2 / 2 + opz / 2 * hp(rB)
    return cert_E(m, eta_c), l, Lz, Leta


def _ball(prec=None):
    """python-flint's arb, acb at CERT_BALL_PREC bits (referee round 104: the ball-arithmetic recheck), or at prec."""
    from flint import arb, acb, ctx
    ctx.prec = CERT_BALL_PREC if prec is None else prec
    return arb, acb


def cert_seg_ball(p, q, ec, delta, R2):
    """`cert_seg_parts`' bound in ball arithmetic (referee round 104, D3cy-2): the same formula on the exact binary inputs p, q (floats),
    ec (acb), delta and R2 (arb), every quantity enclosed and every monotone term (h′, |z|'s and |σ|'s bounds) taken at its rigorous
    upper end, E enclosed on the ball (p + q)/2. Returns the bound's rigorous upper end as a float rounded up (∞ where the segment may
    reach z = ±1)."""
    arb, acb = _ball()
    up = lambda x: x.upper()
    hp = lambda r: arb(0) if not (r > 1) else up(r - 1 / r)
    P, Q = acb(p.real, p.imag), acb(q.real, q.imag)
    m, l = (P + Q) / 2, up(abs(Q - P) / 2)
    c = ec * ec / 4
    Zmax = up(abs(m) + l)
    W = (abs(1 - m * m) - l * (2 * abs(m) + l)).lower()
    if not (W > 0):
        return math.inf
    Lg = up(abs(2 * m * (c - 1 / (1 - m * m))) + l * (2 * abs(c) + 2 * (1 + Zmax * Zmax) / (W * W)))
    Rc = abs(ec)
    sm = ec / 2 * (1 + m)
    Lz = up(Lg + Rc / 2 * hp(up(abs(sm) + Rc / 2 * l)))
    opz = up(abs(1 + m) + l)
    Leta = up(R2 * Zmax * Zmax / 2 + opz / 2 * hp(up(R2 / 2 * opz)))
    rho = abs(sm)
    if up(rho) <= 1:
        C = arb(0)
    else:
        h = (rho * rho - 1) / 2 - rho.log()
        if rho.lower() > 1:
            al = abs(sm.real).lower() - 1
            C = up(h - (al * al if al > 0 else arb(0)))
        else:
            C = up(h) if up(h) > 0 else arb(0)
    Em = up(abs(1 - m * m).log() + (c * m * m).real + C)
    u = up(Em + l * Lz + delta * Leta)
    return math.nextafter(float(u), math.inf)


def cert_box_disc_ball(R1, R2, f1, f2, ec):
    """The box's largest distance from ec at a corner, rigorously (an arb upper end), the angle f2 taken as the true π/2 where it is
    the float π/2 (so the top boxes cover the sliver up to the imaginary axis; referee round 104, D3cy-2)."""
    arb, acb = _ball()
    E_ = acb(ec.real, ec.imag)
    best = arb(0)
    for R in (arb(R1), arb(R2)):
        for f in (arb(f1), arb.pi() / 2 if f2 == math.pi / 2 else arb(f2)):
            v = abs(acb(R * f.cos(), R * f.sin()) - E_).upper()
            best = v if v > best else best
    return best


def cert_poly(pts, eta_c, delta, R2, ball=None):
    """`cert_seg_ub` below −CERT_SLACK on every segment of the polyline, bisecting a segment until it is (or shorter than CERT_LMIN: a
    failure); the children share the float midpoint, so the verified segments tile the polyline exactly. With ball = (ec, δ, R2) in
    ball arithmetic (`cert_seg_ball`) every accepted segment is rechecked rigorously and bisected further unless that bound is below 0.
    Returns (ok, segments verified, the largest float bound, the largest rigorous bound, the largest rigorous − float)."""
    n, worst, bw, bx = 0, -math.inf, -math.inf, -math.inf
    stack = [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)][::-1]
    while stack:
        p, q = stack.pop()
        ub = cert_seg_ub(p, q, eta_c, delta, R2)
        if ub < -CERT_SLACK:
            u = None if ball is None else cert_seg_ball(p, q, *ball)
            if u is None or u < 0:
                n, worst = n + 1, max(worst, ub)
                if u is not None:
                    bw, bx = max(bw, u), max(bx, u - ub)
                continue
        if abs(q - p) < CERT_LMIN:
            return False, n, ub, bw, bx
        m = (p + q) / 2
        stack += [(m, q), (p, m)]
    return True, n, worst, bw, bx


def cert_end(R2):
    """The last piece, from a point z_n with |1 − z_n| ≤ CERT_EPS0 straight to z = 1 (σ = η): there log|1 − z²| ≤ log CERT_EPS0 +
    log(2 + CERT_EPS0), Re(cz²) ≤ (R2²/4)(1 + CERT_EPS0)², C ≤ h((R2/2)(2 + CERT_EPS0)): their sum, which must be negative."""
    hh = lambda r: (r * r - 1) / 2 - math.log(r)
    return (math.log(CERT_EPS0) + math.log(2 + CERT_EPS0) + R2 ** 2 / 4 * (1 + CERT_EPS0) ** 2 + hh(R2 / 2 * (2 + CERT_EPS0)))


def cert_b(R2):
    """b's z = t on the straight segment: CERT_TB of the disc's edge (2 − R2)/R2, inside the disc for every |η| ≤ R2."""
    return CERT_TB * (2 - R2) / R2


def cert_box_disc(R1, R2, f1, f2):
    """The box's centre η_c (at the mean |η| and arg η) and its largest distance from it, at a corner (|η − η_c|² is convex in |η|
    and largest at the extreme angle)."""
    ec = (R1 + R2) / 2 * cmath.exp(1j * (f1 + f2) / 2)
    return ec, max(abs(R * cmath.exp(1j * f) - ec) for R in (R1, R2) for f in (f1, f2))


def cert_box(R1, R2, f1, f2, ball=False):
    """The valley certificate on the box R1 ≤ |η| ≤ R2, f1 ≤ arg η ≤ f2: one path for the whole box — from b = CERT_TB·(2 − R2)/R2 (on
    the straight segment, inside the disc for every η of the box) along the box centre's descent flow to |1 − z| < CERT_EPS0, then
    straight to z = 1 —, `cert_poly` with δ the box's largest distance from its centre (at a corner). With ball (referee round 104) also
    in ball arithmetic: δ from the corners rigorously (the top edge at the true π/2), (R2/2)(1 + t) < 1 and |1 − z_n| ≤ CERT_EPS0
    rigorously, every segment by `cert_seg_ball`. Returns (ok, segments, worst, the largest rigorous bound, the largest rigorous − float,
    whether the float δ fell short of the rigorous one)."""
    ec, d = cert_box_disc(R1, R2, f1, f2)
    pts = cert_flow(ec, cert_b(R2))
    if pts is None or abs(1 - pts[-1]) > CERT_EPS0:
        return False, 0, math.inf, math.inf, math.inf, False
    if not ball:
        return cert_poly(pts, ec, d, R2) + (False,)
    arb, acb = _ball()
    db = cert_box_disc_ball(R1, R2, f1, f2, ec)
    short = bool(db > arb(d))
    if not ((arb(R2) / 2 * (1 + arb(cert_b(R2)))).upper() < 1) or not (abs(1 - acb(pts[-1].real, pts[-1].imag)).upper() <= arb(CERT_EPS0)):
        return False, 0, math.inf, math.inf, math.inf, short
    if not (abs(acb(ec.real, ec.imag)) < arb(R2)):                    # round 108 (D3dc-3): L_η's segment [η_c, η] within |η| ≤ R2
        return False, 0, math.inf, math.inf, math.inf, short
    dd = db if short else arb(d)
    return cert_poly(pts, ec, max(d, math.nextafter(float(dd), math.inf)), R2, ball=(acb(ec.real, ec.imag), dd, arb(R2))) + (short,)


def cert_end_ball(R2):
    """`cert_end` in ball arithmetic: its rigorous upper end (a float rounded up; referee round 104)."""
    arb, acb = _ball()
    e, R = arb(CERT_EPS0), arb(R2)
    r = R / 2 * (2 + e)
    return math.nextafter(float((e.log() + (2 + e).log() + R * R / 4 * (1 + e) ** 2 + (r * r - 1) / 2 - r.log()).upper()), math.inf)


def cert_band(Ra, Rb, nR, nT, ball=False):
    """`cert_box` over Ra ≤ |η| ≤ Rb, 0 ≤ arg η ≤ π/2 from an nR × nT grid, a failing box halved along its longer side (in |η| or in
    |η|·arg η) down to CERT_DEPTH levels. Returns the boxes verified, the segments, the largest bound, the boxes failed at the deepest
    level, the smallest box side, and the end piece's bound; with ball also the largest rigorous bound, the largest rigorous − float,
    the boxes whose float δ fell short, and the end piece's rigorous bound (referee round 104)."""
    edge = lambda i: Rb if i == nR else Ra + (Rb - Ra) * i / nR        # round 108 (D3dc-3): the last edge Rb exactly
    stack = [(edge(i), edge(i + 1), (math.pi / 2) * j / nT, (math.pi / 2) * (j + 1) / nT, 0) for i in range(nR) for j in range(nT)]
    done, segs, worst, fails, side = 0, 0, -math.inf, 0, math.inf
    bw, bx, nshort = -math.inf, -math.inf, 0
    while stack:
        R1, R2, f1, f2, dep = stack.pop()
        ok, n, w, bw_, bx_, short = cert_box(R1, R2, f1, f2, ball)
        if ok:
            done, segs, worst, side = done + 1, segs + n, max(worst, w), min(side, R2 - R1, R2 * (f2 - f1))
            bw, bx, nshort = max(bw, bw_), max(bx, bx_), nshort + short
            continue
        if dep >= CERT_DEPTH:
            fails += 1
            continue
        if R2 - R1 > R2 * (f2 - f1):
            m = (R1 + R2) / 2
            stack += [(R1, m, f1, f2, dep + 1), (m, R2, f1, f2, dep + 1)]
        else:
            m = (f1 + f2) / 2
            stack += [(R1, R2, f1, m, dep + 1), (R1, R2, m, f2, dep + 1)]
    out = dict(R1=Ra, R2=Rb, boxes=done, segs=segs, worst=worst, fails=fails, side=side, end=cert_end(Rb))
    if ball:
        out.update(ball_worst=bw, ball_excess=bx, ball_short=nshort, ball_end=cert_end_ball(Rb))
    return out


def part_V():
    """V (referee round 99's recipe, carried out; round 104: the tight gradient bound and the ball-arithmetic recheck): the valley
    condition certified on the bands of CERT_BANDS (`cert_band` with ball)."""
    import flint
    return dict(bands=[cert_band(*b, ball=True) for b in CERT_BANDS], tb=CERT_TB, eps0=CERT_EPS0, slack=CERT_SLACK,
                ball=dict(prec=CERT_BALL_PREC, flint=flint.__version__))


def part_dpath():
    """D (round 93): `path_ok` at DP_TH angles in [0, π/2] for each radius of DP_R (grid step DP_H)."""
    ths = np.linspace(0, np.pi / 2, DP_TH)
    rows = []
    for R, h in zip(DP_R, DP_H):
        res = [path_check(R * np.exp(1j * t), h) for t in ths]
        rows.append(dict(R=R, h=h, ok=[x[0] for x in res], ma=[x[1] for x in res], mb=[x[2] for x in res], q=[x[3] for x in res]))
    return dict(th=[float(t) for t in ths], rows=rows, seeds='descent')


def d_lines(D, S=None, R=None):
    import cmath
    import math
    if D is None or not isinstance(D, dict) or not D.get('jobs'):
        return []
    out = ['', "## D. The reduction: Duhamel's identity, the Laplace step and the domination condition (referee round 91)", '',
           "Duhamel's identity in s: ĝ_K(ξ)e^{−ξ²/4} = K·C(2K, K)∫_0^1 s^{K−1}(1 − s)^K e^{ξ²(1 − 2s)²/4}c_K(ξs)ds, ξ = √Kη, by Gauss–Legendre "
           "quadrature in double (c_K by `corner_ratio`) over `gen_ratio_c`, minus 1:", '',
           '| K | \\|η\\| | arg η | points | quadrature over the exact, minus 1 |', '|---|---|---|---|---|']
    for x in D['quad']:
        q, e = complex(*x['quad']), complex(*x['exact'])
        out.append(f"| {x['K']} | {x['eta']:g} | {x['j8']}π/8 | {x['n']} | {abs(q / e - 1):.1e} |")
    out += ['', "The Laplace step at s = ½ (the phase's stationary point for every η): |ĝ_K(ξ)e^{−ξ²/4}(1 − η²/4)^{1/2}/c_K(ξ/2) − 1| "
            "(`gen_ratio_c`, `corner_ratio`; in brackets the corner ratio's own |c_K(ξ/2)/(1 − η²/4) − 1|):", '',
            '| K | arg η | ' + ' | '.join(f"\\|η\\| = {e:g}" for e in D['jobs'][0]['eta']) + ' |', '|---|---|' + '---|' * len(D['jobs'][0]['eta'])]
    dev = {}
    for x in D['jobs']:
        row = []
        for e, a, b in zip(x['eta'], x['g'], x['c']):
            mu2 = (e * cmath.exp(1j * x['j8'] * math.pi / 8)) ** 2 / 4
            st = abs(complex(*a) * cmath.sqrt(1 - mu2) / complex(*b) - 1)
            dev.setdefault((x['j8'], e), []).append((x['K'], st))
            row.append(f"{st:.1e} ({abs(complex(*b) / (1 - mu2) - 1):.1e})")
        out.append(f"| {x['K']} | {x['j8']}π/8 | " + ' | '.join(row) + ' |')
    Kl = max(x['K'] for x in D['jobs'])
    js = sorted({j for j, _ in dev})
    out += ['', f"K times the step's deviation at K = {Kl} (|η| = " + ', '.join(f"{e:g}" for e in D['jobs'][0]['eta']) + '), by arg η: '
            + '; '.join(f"{j}π/8: " + ', '.join(f"{Kl * dict(dev[(j, e)])[Kl]:.3f}" for e in D['jobs'][0]['eta']) for j in js)
            + f"; the ratio of the last two deviations (K⁻¹'s {Kl // 2}/{Kl} = 0.5): " + '; '.join(
                f"{j}π/8: " + ', '.join(f"{dev[(j, e)][-1][1] / dev[(j, e)][-2][1]:.3f}" for e in D['jobs'][0]['eta']) for j in js) + '.']
    dm = D['dom']
    i = int(np.argmin(dm['radius']))
    R_dom = dm['radius'][i]
    vl = lambda k_: ('none' if dm['viol'][k_] is None else f"arg η ∈ [{dm['viol'][k_][0]:.3f}, {dm['viol'][k_][1]:.3f}]")
    out += ['', "The domination condition — sup over s ∈ [1/|η|, 1) of g(s) + Λ⁺(ηs) below g(½), g(s) = log(s(1 − s)) + (Re η²/2)((1 − s)² + "
            "s²), Λ⁺(μ) = M(Re μ) − ½ − log|μ| − Re(μ²)/2 (`lam_plus`), the corner ratio's bound beyond Q's disc —, the supremum over s in "
            "closed form (`gap_closed`; referee round 93: with w = 1 − s, log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4, η = "
            f"Re^{{iφ}}, p = Re η²): it holds at every one of {len(dm['th'])} angles in [0, π/2] for |η| < {R_dom:.6f} (the least radius, at arg "
            f"η = {dm['th'][i]:.3f}); its gap at |η| = 1.999 on the real axis {dm['axes'][0]:.2e}, on the imaginary {dm['axes'][1]:.2e}; "
            f"violated at |η| = 1.95 for at least {vl('1.95')} and at 1.999 for at least {vl('1.999')}"
            + ('' if not dm.get('margin') else "; the gap's largest value on the grid (" + str(DOM_NR) + " radii from 1 to c times the angles) "
               + ', '.join(f"{v[0]:.2e} at |η| ≤ {c} (at arg η = {v[1]:.3f})" if isinstance(v, (list, tuple)) else f"{v:.2e} at |η| ≤ {c}"
                           for c, v in dm['margin'].items())) + '.']
    if dm.get('bb'):
        out += ['', "A branch and bound on the closed form (`dom_bb`: on a box of (|η|, arg η) the least radius, the least Re η² and the "
                    f"largest |η|cos(arg η) bound the gap from above; floating point with slack {DOM_SLACK:g}"
                    + (", every accepted box rechecked in ball arithmetic — `box_bound_ball`: the same bound with every parameter at its "
                       "rigorous end, the top edge at the true π/2, the supremum's candidates as balls; referee round 108 —"
                       if all('ball_worst' in x for x in dm['bb']) else ", not interval arithmetic") + ") verifies "
                    "the condition on " + '; '.join(f"1 < |η| ≤ {x['R']:g}, 0 ≤ arg η ≤ π/2 ({x['n']} boxes, the largest bound {x['worst']:.1e}"
                                                    + (f", in ball arithmetic {x['ball_worst']:.4e}" if 'ball_worst' in x else '')
                                                    + f", {x['open']} unresolved)" for x in dm['bb']) + '.']
        sc = dm.get('sect')
        if sc:
            fl, cl = (lambda v, d: math.floor(v * 10 ** d) / 10 ** d), (lambda v, d: math.ceil(v * 10 ** d) / 10 ** d)
            out[-1] += (f" The violated sector at |η| = {sc['R']:g} by bisection on the closed form (`dom_sector`): arg η ∈ "
                        f"({fl(sc['lo'], 5):.5f}, {cl(sc['hi'], 5):.5f}), and there the supremum's maximiser s has Re(ηs) between "
                        f"{fl(sc['re_min'], 4):.4f} and {cl(sc['re_max'], 4):.4f} (numerically: floating point at {sc['n']} angles and the ends, the "
                        "printed digits rounded outward; referee round 131) — beyond 1, "
                        "where Λ⁺'s M(Re ηs) is its linear piece (`gap_arg`; referee rounds 95, 99).")
    fr = []
    for lab, xs in (('R', [] if R is None or not isinstance(R, dict) else list(R['rules'])), ('S', [] if S is None else list(S))):
        if xs:
            rs = [(x, 1 / (2 * min(x['gx'], x['gy']))) for x in xs]
            w = max(rs, key=lambda z: z[1])
            fr.append(f"{lab}'s {len(xs)}: r ≤ {w[1]:.4f}, |η| ≤ {c_need(w[1]):.4f} (τ/σ = {w[0]['r']:g}, "
                      f"{w[0].get('name', w[0].get('frame'))})")
    import math as _m
    bbv = [x['R'] for x in dm.get('bb', []) if x['open'] == 0]          # round 93: the radius the branch and bound verified
    R_ver = min(max(bbv), R_dom) if bbv else None                       # round 110 (D3de-4): no fallback to the grid's radius
    r0_ = _m.floor(r_star(1.0) * 1e4) / 1e4                              # round 110: floored, as the others
    if R_ver is not None:
        rs_ = _m.floor(r_star(R_ver) * 1e4) / 1e4                       # referee round 93 (D3cn-6): rounded down
        out += ['', "(iv) needs the limit of (iii) on |η| ≤ c with c² > 4 log u/(u − 1), u = 1/r, r a frame's Laplace rate (`c_need`): "
                + ('; '.join(fr) + '; ' if fr else '') + f"every frame with r < {rs_:.4f} lies within the verified radius {R_ver:g} (c_need({rs_:.4f}) ≤ "
                f"{c_need_ball(rs_):.8f} in ball arithmetic, {'inside' if c_need_ball(rs_) < R_ver else 'NOT inside'} it: referee round 110) — at the half-squeeze rule (r = t) σ*/σ < (1 + r)/(1 − r) = "
                f"{ratio_floor(rs_)}; below r = {r0_:.4f} (|η| ≤ 1, where the condition is void; c_need({r0_:.4f}) ≤ {c_need_ball(r0_):.8f}"
                f"{' < 1' if c_need_ball(r0_) < 1 else ', NOT below 1'}) no "
                "check is needed."]
    else:
        out += ['', "(iv) needs the limit of (iii) on |η| ≤ c with c² > 4 log u/(u − 1), u = 1/r, r a frame's Laplace rate (`c_need`): "
                + ('; '.join(fr) + '; ' if fr else '') + "no branch and bound here resolved every box, so no radius beyond |η| ≤ 1 is "
                f"verified (the angle grid's least radius {R_dom:.6f} is unverified); below r = {r0_:.4f} (|η| ≤ 1, where the "
                "condition is void) no check is needed."]
    if D.get('qx'):
        qx = sorted(D['qx'], key=lambda x: x['K'])
        ys = qx[0]['y']
        lo, hi = 1.0, 3.0
        for _ in range(60):
            m = (lo + hi) / 2
            lo, hi = (lo, m) if alias_rate(m) > 0 else (m, hi)
        inside = [y for y in ys if alias_rate(y) < 0]
        beyond = [y for y in ys if alias_rate(y) >= 0]
        Ks = [x['K'] for x in qx]
        cell = lambda k, y: qx[k]['c'][ys.index(y)]
        out += ['', "Beyond Q's disc on the imaginary axis, c_K(i√Ky) (real there; 4K + 80 digits): over 1 + y², minus 1, at K = "
                + ', '.join(str(K) for K in Ks) + ': ' + '; '.join(f"y = {y:g}: " + ', '.join(f"{cell(k, y) / (1 + y * y) - 1:+.1e}"
                                                                           for k in range(len(qx))) for y in inside)
                + f"; beyond y = {lo:.4f}, where the first alias's rate (5/8)y² − log y − ½ − (3/2)log 2 (`alias_rate`; referee round 93) "
                "vanishes, c_K − (1 + y²) grows" + ('' if not D.get('qxp') else ": Prony fits (`prony`, a two-term recurrence) over "
                f"{QXP_N} consecutive K from K = " + ', '.join(str(x['K0']) for x in D['qxp']) + ' — ' + '; '.join(
                    f"y = {y:g}: the rate " + ', '.join(f"{prony([r['c'][j] for r in x['rows']], alias_rate(y))[0]:.5f}" for x in D['qxp'])
                    + f" (the alias {alias_rate(y):.5f}), the frequency " + ', '.join(
                        f"{prony([r['c'][j] for r in x['rows']], alias_rate(y))[1]:.4f}" for x in D['qxp'])
                    for j, y in enumerate(D['qxp'][0]['y'])) + ' — (two-point slopes oscillate with the frequency and do not measure it)') + '.']
    if D.get('dpath'):
        dp = D['dpath']
        rf_ = math.floor(r_star(max(x['R'] for x in dp['rows'])) * 1e4) / 1e4      # floored, and the ratio from the floored r
        out += ['', "The deformed Duhamel path (`path_ok`): with the bound F(σ) = log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + C(σ) (C = 0 "
                    "inside the unit disc, Λ⁺ outside) and the ball B of radius rb = ½(1 − |η|/2) around η/2, the two seeds η/2 ∓ "
                    f"{PATH_SEED:g}·rb·u on the saddle's steepest-descent line (`path_seeds`: u² ∝ −conj φ″, φ″ = 2 − 8/η²; within 2rb of η/2, so in "
                    "the unit disc) lie in the components of {F < F(η/2)} ∖ B (8-connected on a grid) that contain 0 and η, those two components "
                    "differ, and the chord between the seeds stays below F(η/2) but at η/2 — by Taylor's theorem, φ‴(η/2) = 0, where the ratio "
                    "q = D²/((|η|/2 − D)⁴|φ″|) < 1, D = |b − η/2| (`chord_q`; referee round 99) —, at "
                    + '; '.join(
                        f"|η| = {x['R']:g}: {sum(x['ok'])} of {len(x['ok'])} angles (grid step {x['h']:g}"
                        + (f"; the seeds' F − F(η/2) {max(max(x['ma']), max(x['mb'])):.1e} … {min(min(x['ma']), min(x['mb'])):.1e}"
                           if x.get('ma') else '')
                        + (f", the chord's ratio q at most {max(x['q']):.1e}" if x.get('q') else '') + ")" for x in dp['rows'])
                + f" — {len(dp['th'])} angles from {min(dp['th']):g} to {max(dp['th']) / math.pi:g}π, the sector of the straight segment's "
                  "failure among them; the components a grid check, numerical (referee round 95; part V certifies a valley condition on "
                  "the annulus: below). Were the condition to hold on the annulus "
                  f"{1.0 if R_ver is None else R_ver:g} < |η| ≤ {max(x['R'] for x in dp['rows']):g} (the branch and bound covers |η| ≤ "
                  f"{1.0 if R_ver is None else R_ver:g}), "
                  f"(iv) would cover every frame with r < {rf_:.4f} (σ*/σ < {ratio_floor(rf_)} at the "
                  "half-squeeze rule)."]
    return out


def ratio_floor(r):
    """(1 + r)/(1 − r) in exact arithmetic from r's 4-decimal string, rounded down at 2 decimals (referee round 104, D3cy-5: in floating
    point (1 + 0.99)/(1 − 0.99) = 198.999…83 floored to 198.99)."""
    q = Fraction(int(round(r * 10 ** 4)), 10 ** 4)
    return f"{math.floor((1 + q) / (1 - q) * 100) / 100:.2f}"


def v_lines(V, D=None):
    """Part V's lines: the valley certificate per band and what it gives (iv) with the branch and bound below it."""
    if V is None or not V.get('bands'):
        return []
    bs = sorted(V['bands'], key=lambda b: b['R1'])
    dm = (D or {}).get('dom', {})
    R_bb = max([x['R'] for x in dm.get('bb', []) if x['open'] == 0] or [0.0])
    R_v = R_bb
    for b in bs:                                                         # the certified radius: bands contiguous from the branch and bound's
        if b['fails'] == 0 and b['end'] < 0 and b['R1'] <= R_v:         # round 108 (D3dc-3): no gap at all
            R_v = max(R_v, b['R2'])
    rv = math.floor(r_star(R_v) * 1e4) / 1e4
    out = ['', '## V. The valley condition certified (referee round 99\'s recipe; in ball arithmetic: round 104)', '',
           "In the saddle's coordinates σ = (η/2)(1 + z) the bound is F_η(σ) − F_η(η/2) = log|1 − z²| + Re(cz²) + C(σ), c = η²/4 (`cert_E`): on "
           "the straight segment inside the unit disc, z = x real with (|η|/2)(1 + x) ≤ 1, it is at most log(1 − x²) + x²·Re c < log(1 − x²) + "
           "x² < 0 for x ≠ 0 (|c| < 1), so the path from 0 through η/2 to b = η/2 + t·η/2 (0 < t < (2 − |η|)/|η|) needs no check, Laplace's "
           "point crossed along the real z-axis (Re(c − 1) < 0); from b to η one path per η-box — the box centre's descent flow of the bound "
           f"from b at t = {V['tb']:g}·(2 − R2)/R2, then straight to z = 1 from |1 − z| ≤ {V['eps0']:g} — with an upper bound of the bound over "
           "each segment and the whole box (`cert_seg_ub`: the value at the midpoint and the box centre, plus first-order terms in z and in η "
           "from explicit gradient bounds, C's |∇_σC| ≤ h′(|σ|) exactly (referee round 104); in floating point with slack "
           f"{V['slack']:g}), a segment bisected and a box halved until every bound is below −slack: " + '; '.join(
               f"{b['R1']:g} ≤ |η| ≤ {b['R2']:g}, 0 ≤ arg η ≤ π/2: {b['boxes']} boxes verified ({b['fails']} failed), {b['segs']} segment bounds, "
               f"the largest {b['worst']:.1e}, the smallest box side {b['side']:.1e}, the end piece's bound {b['end']:.3f}" for b in bs)
           + '.']
    ball = V.get('ball') and all('ball_worst' in b for b in bs)
    if ball:                                                            # referee round 104 (D3cy-2): the recheck in ball arithmetic
        out[-1] += (f" Every accepted segment rechecked in ball arithmetic (`cert_seg_ball`: Arb at {V['ball']['prec']} bits, python-flint "
                    f"{V['ball']['flint']}; the same formula on the exact binary inputs, every monotone term at its rigorous upper end): the "
                    "largest rigorous bound " + ', '.join(f"{b['ball_worst']:.4e}" for b in bs) + " (by band), at most "
                    + ', '.join(f"{b['ball_excess']:.1e}" for b in bs) + " above the float one; each box's covering disc from its corners "
                    "rigorously, the top edge at the true π/2 (the float radius short of it at " + ', '.join(str(b['ball_short']) for b in bs)
                    + " boxes, the rigorous one used), b inside the unit disc and the path's end within the end piece's radius rigorously, "
                      "the end piece's rigorous bound " + ', '.join(f"{b['ball_end']:.6f}" for b in bs) + '.')
    ok_b = all(b.get('ball_worst', 0) < 0 and b.get('ball_end', 0) < 0 for b in bs) if ball else False
    if R_bb == 0.0:                                                     # round 110 (D3de-1): no verified branch and bound below it
        out[-1] += (" The bound is even in η ↦ η̄ and η ↦ −η (z ↦ z̄; σ ↦ −σ), so the quadrant covers the annulus. No branch and bound "
                    "below the annulus resolved every box, so the certificate gives (iii)(b) no radius here.")
        return out
    out[-1] += (f" The bound is even in η ↦ η̄ and η ↦ −η (z ↦ z̄; σ ↦ −σ), so the quadrant covers the annulus. With the branch and bound "
                f"on |η| ≤ {R_bb:g}, the limit of (iii)(b) holds on |η| ≤ {R_v:g} ("
                + (("both certificates in ball arithmetic: referee round 108; since round 148 cross-checks of the hand proof" if all(
                    x.get('ball_worst', 0) < 0 and x['open'] == 0 for x in dm.get('bb', [])) and dm.get('bb') else
                    "modulo the branch and bound, the valley condition certified in ball arithmetic") if ok_b else
                   "modulo both certificates") + f"), and (iv) for every frame with r < {rv:.4f} (c_need({rv:.4f}) ≤ {c_need_ball(rv):.8f} in ball arithmetic, "
                f"{'inside' if c_need_ball(rv) < R_v else 'NOT inside'} |η| ≤ {R_v:g}: referee round 110) — "
                "at the half-squeeze rule "
                f"σ*/σ < {ratio_floor(rv)} (in exact arithmetic from r's floor: referee round 104).")
    return out


VC_EPS = (1e-2, 1e-3)                                                   # VC: ε = 2 − |η| for the blow-up's check
VC_PHI = (0.0, 0.5, 1.0, 1.44, 2.0, 5.0)                                # VC: φ = arg η / ε
VC_X = (0.05, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)                            # VC: the floor's points x (w = x + iy)
VC_R = (1.99, 1.995, 1.999)                                             # VC: radii for the straight segment's sector near the real axis


def corner_psi(x, y, phi):
    """The corner's blow-up limit (HANDOFF §7; part VC): with ε = 2 − |η|, θ = arg η = εφ and z = ε(x + iy), `cert_E` = ε³Ψ_φ + o(ε³),
    Ψ_φ(x, y) = −(x² − y²) − 4φxy + (x − ½)₊[(y + φ)² − (x − ½)²/3] — the quadratic part Re((η²/4 − 1)z²) at its leading order
    (η²/4 − 1 ≈ ε(−1 + 2iφ)), and outside the unit disc (x > ½) C = h(|σ|) − (|Re σ| − 1)₊² with h(1 + δ) = δ² − δ³/3 and |σ| − 1,
    Re σ − 1 = ε(x − ½) + ε²((y² − x)/2, −φ²/2 − φy − x/2) + …, whose squares differ by ε³(x − ½)(y + φ)²."""
    u = max(x - 0.5, 0.0)
    return -(x * x - y * y) - 4 * phi * x * y + u * ((y + phi) ** 2 - u * u / 3)


def corner_floor(x, phi):
    """Ψ_φ's valley floor: Ψ_φ is convex in y and its minimiser is y = 2φx inside the disc (x ≤ ½) and y = φ outside, where
    Ψ_φ = −(1 + 4φ²)x² and −x² − φ² − (x − ½)³/3 — negative but at w = 0, for every φ ≥ 0."""
    return min(2 * phi * x, phi)


def corner_phi_dom():
    """The straight segment's limit (y = 0) beyond the disc: max over u = x − ½ ≥ 0 of −(u + ½)² + φ²u − u³/3 (below ½ it is −x², 0 at the
    saddle), negative iff φ < φ_dom, the root of 8φ³ − 12φ² + 1 = 0 ((φ − 1)²(2φ + 1) = ¾);
    φ_dom by bisection on that maximum (on a grid of u, step 10⁻⁴). Returns (φ_dom, the maximum at φ = 0)."""
    def top(phi):
        us = np.linspace(0.0, 4.0, 40001)
        return float(np.max(-(us + 0.5) ** 2 + phi * phi * us - us ** 3 / 3))
    lo, hi = 1.0, 2.0
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if top(m) < 0 else (lo, m)
    return (lo + hi) / 2, top(0.0)


def corner_sector(R, top=0.05):
    """The straight segment's good sector near the real axis at |η| = R: the least arg η ≤ top where `gap_closed` > 0 (bisection from a
    grid of 2001 angles), None if none. Floating point: valid for ε = 2 − R ≥ 10⁻³ (referee round 115: noise below ε ≈ 10⁻⁴)."""
    ph = np.linspace(0.0, top, 2001)
    i = next((k for k, f in enumerate(ph) if gap_closed(R, float(f)) > 0), None)
    if i is None or i == 0:
        return None
    a, b = float(ph[i - 1]), float(ph[i])
    for _ in range(60):
        m = (a + b) / 2
        a, b = (a, m) if gap_closed(R, m) > 0 else (m, b)
    return (a + b) / 2


VC_PHI1 = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0)                               # VC (round 115): φ for the first correction's coefficient
VC_FAIL = (1e-3, (20.0, 25.0, 30.0), (0.05, 0.25, 0.45))                # VC (round 115): ε, φ and x where Ψ_φ's floor is tried far out
VC_MID = ((1.0, 1.0, 0.5), (2.0, 0.3, 0.4), (3.0, 5.0, 3.0), (0.2, 3.0, 0.1))   # VC (rounds 115, 117: the last with Re σ > 1) for Ψ̃_ψ
VC_MID_EPS = ('1e-4', '1e-6', '1e-8')


VC_PREC = 320                                                           # VC: bits for the remainder's coefficient (ε down to 10⁻⁸)


def _cert_E_ball(z, eta):
    """`cert_E` in ball arithmetic (Arb at the current precision; z, η acb): the same centred form, C exact (referee round 115: float
    cannot resolve the ε³ margins below ε ≈ 10⁻⁴)."""
    from flint import arb
    c = eta * eta / 4
    sg = eta / 2 * (1 + z)
    rho = abs(sg)
    C = arb(0)
    if rho > 1:
        a = abs(sg.real) - 1
        C = (rho * rho - 1) / 2 - rho.log() - (a * a if a > 0 else arb(0))
    elif not (rho <= 1):
        raise ValueError('|σ| straddles 1')
    return abs(1 - z * z).log() + (c * z * z).real + C


def corner_psi1(x, y, phi):
    """Ψ_φ's first correction inside the unit disc (referee round 115, D3dj-1): with η²/4 − 1 = ε(−1 + 2iφ) + ε²(¼ − 2φ² − 2iφ) + O(ε³)
    and log(1 − z²) + z² = −z⁴/2 + …, `cert_E`/ε³ = Ψ_φ + εΨ¹_φ + O(ε²), Ψ¹_φ(w) = Re((¼ − 2φ² − 2iφ)w²) − Re(w⁴)/2. On the floor
    y = 2φx its ratio to Ψ_φ is −(¼ + 5φ² + 8φ⁴ − x²Re((1 + 2iφ)⁴)/2)/(1 + 4φ²) ≈ −2φ²: the expansion holds for θ = εφ ≪ √ε only."""
    w = complex(x, y)
    return ((0.25 - 2 * phi * phi - 2j * phi) * w * w).real - (w ** 4).real / 2


def corner_psi_mid(x, y, psi):
    """The limit at θ = √εψ (referee round 115, D3dj-1): with z = εx + i√εy, `cert_E`/ε² → Ψ̃_ψ(x, y) = (1 + 2ψ²)y² − 4ψxy − y⁴/2 +
    D₊² − P₊², D = x − ½ + y²/2 (|σ| − 1 over ε) and P = x − ½ − ψ²/2 − ψy (Re σ − 1 over ε)."""
    D, P = x - 0.5 + y * y / 2, x - 0.5 - psi * psi / 2 - psi * y
    return (1 + 2 * psi * psi) * y * y - 4 * psi * x * y - y ** 4 / 2 + np.maximum(D, 0) ** 2 - np.maximum(P, 0) ** 2


def corner_mid_path(psi, xs):
    """Ψ̃_ψ along the path y = min(2ψx/(1 + 2ψ²), ψ): inside the disc −4ψ²x²/(1 + 2ψ²) − y⁴/2 (the slope minimises the quadratic), far
    out −ψ² − ψ⁴/2 exactly (D − P = (y + ψ)²/2 and the linear terms cancel at y = ψ); negative for every x > 0 and ψ > 0 (proved by
    cases in D61; round 115). Returns the values on xs."""
    return corner_psi_mid(xs, np.minimum(2 * psi / (1 + 2 * psi * psi) * xs, psi), psi)


def part_VC():
    """VC (the corner at η = 2; HANDOFF §7): along the valley floor w = x + i·`corner_floor`(x, φ) the true bound `cert_E` at η = (2 − ε)e^{iεφ},
    z = εw, against ε³Ψ_φ, and the analytic floor values; the straight segment's φ_dom against the measured sector's lower end over ε."""
    rows = []
    for eps in VC_EPS:
        for phi in VC_PHI:
            eta = (2 - eps) * cmath.exp(1j * eps * phi)
            pts = []
            for x in VC_X:
                y = corner_floor(x, phi)
                E_ = cert_E(eps * complex(x, y), eta)
                pts.append(dict(x=x, y=y, E=E_, psi=corner_psi(x, y, phi)))
            rows.append(dict(eps=eps, phi=phi, pts=pts))
    pd, g0 = corner_phi_dom()
    sec = [dict(R=R, lo=corner_sector(R)) for R in VC_R]
    arb, acb = _ball(VC_PREC)                                           # round 115 (D3dj-1): the remainder, far out, and θ ~ √ε
    eta_b = lambda eps, th: (2 - eps) * acb(th.cos(), th.sin())
    coef = []
    for phi in VC_PHI1:
        y = corner_floor(0.05, phi)
        cs = []
        for e in ('1e-6', '1e-8'):
            eps = arb(e)
            E_ = _cert_E_ball(eps * acb(arb('0.05'), arb(y)), eta_b(eps, eps * arb(phi)))
            xa, ya, pa = arb('0.05'), arb(y), arb(phi)                   # round 117 (D3dl-5): Ψ_φ in balls too (inside the disc)
            ps = -(xa * xa - ya * ya) - 4 * pa * xa * ya
            cs.append(float(((E_ / eps ** 3 - ps) / (eps * ps)).mid()))
        coef.append(dict(phi=phi, x=0.05, eps=['1e-6', '1e-8'], ball=cs, closed=corner_psi1(0.05, y, phi) / corner_psi(0.05, y, phi)))
    eps, phis, xs_ = VC_FAIL
    far = []
    for phi in phis:
        eta = (2 - eps) * cmath.exp(1j * eps * phi)
        for x in xs_:
            yf = corner_floor(x, phi)
            far.append(dict(phi=phi, x=x, floor=cert_E(eps * complex(x, yf), eta),
                            least=min(cert_E(eps * complex(x, float(yy)), eta) for yy in np.linspace(0.0, 3 * yf, 3001))))
    mid = []
    for psi, x, y in VC_MID:
        vals = []
        for e in VC_MID_EPS:
            ee = arb(e)
            vals.append(float((_cert_E_ball(acb(ee * x, ee.sqrt() * y), eta_b(ee, ee.sqrt() * psi)) / ee ** 2).mid()))
        mid.append(dict(psi=psi, x=x, y=y, eps=list(VC_MID_EPS), vals=vals, lim=float(corner_psi_mid(x, y, psi)),
                        P=float(x - 0.5 - psi * psi / 2 - psi * y)))
    pg = np.unique(np.r_[np.linspace(0.02, 1.0, 50), np.linspace(1.0, 20.0, 96)])   # round 117: ψ = 1 once
    xq = lambda q: np.r_[np.linspace(1e-4, 1.0, 2001), np.linspace(1.0, max(200.0, 2 + 4 * q * q), 20001)]   # past the switch
    sc = lambda q, x: q * q * np.minimum(x, 1.0) ** 2
    path = dict(n_psi=int(pg.size), psi=(float(pg[0]), float(pg[-1])), xmin=float(min(xq(q)[0] for q in pg)),
                xmax=float(max(xq(q)[-1] for q in pg)), xmax0=float(xq(pg[0])[-1]),   # round 117: the grid used, not a literal
                worst=float(max(np.max(corner_mid_path(q, xq(q)) / sc(q, xq(q))) for q in pg)),
                far=float(max(abs(float(corner_mid_path(q, np.array([2 + 4 * q * q]))[0]) / (q * q + q ** 4 / 2) + 1) for q in pg)))
    return dict(rows=rows, phi_dom=pd, gap0=g0, sector=sec, coef=coef, far=far, mid=mid, path=path)


def vc_lines(VC):
    """Part VC's line: the blow-up limit, its floor, the true bound's ratio to ε³Ψ along it, and the straight segment's sector; round
    115: the first correction and its φ² growth, the floor far out, the limit at θ ~ √ε with its path, and the four regimes."""
    if not VC:
        return []
    rat = lambda r: [p_['E'] / (r['eps'] ** 3 * p_['psi']) for p_ in r['pts'] if p_['psi'] != 0]   # round 115 (D3dj-6): Ψ = 0 skipped
    neg = all(p_['E'] < 0 for r in VC['rows'] for p_ in r['pts'])
    epss = list(dict.fromkeys(r['eps'] for r in VC['rows']))           # round 115 (D3dj-6): from the raw, not the constants
    by = {e: [r for r in VC['rows'] if r['eps'] == e] for e in epss}
    out = ['', '## VC. The corner at η = 2: a blow-up limit (not yet a proof)', '',
           "With ε = 2 − |η|, θ = arg η = εφ and the saddle's coordinate z = ε(x + iy), the centred bound is ε³Ψ_φ(x, y) + o(ε³), Ψ_φ = "
           "−(x² − y²) − 4φxy + (x − ½)₊[(y + φ)² − (x − ½)²/3] (`corner_psi`: the quadratic part's leading order and, outside the unit "
           "disc, C's squares differing at order ε³). Ψ_φ is convex in y; its floor y = min(2φx, φ) (`corner_floor`) gives −(1 + 4φ²)x² "
           "inside the disc and −x² − φ² − (x − ½)³/3 outside: negative but at the saddle for every φ ≥ 0. Along the floor at x = "
           + ', '.join(f"{x:g}" for x in dict.fromkeys(p_['x'] for p_ in VC['rows'][0]['pts'])) + " the true bound (`cert_E`) over ε³Ψ_φ: "
           + '; '.join(f"ε = {e:g}: " + ', '.join(f"φ = {r['phi']:g} {min(rat(r)):.3f}–{max(rat(r)):.3f}" for r in by[e]) for e in epss)
           + (", the bound negative at every point" if neg else ", NOT negative at every point") + ". On the straight segment (y = 0) "
           f"beyond the disc (x ≥ ½; at the saddle Ψ is 0) Ψ's maximum is {VC['gap0']:.4f} at φ = 0 (the real axis' gap −¼ε³, part D's −¼(2 − "
           f"η)³) and turns positive at φ_dom = {VC['phi_dom']:.4f} (`corner_phi_dom`; the root of 8φ³ − 12φ² + 1); the measured sector of the "
           "straight segment near the real axis (`corner_sector`, the domination gap's first violated arg η; floating point, valid for ε ≥ "
           "10⁻³) over ε: " + ', '.join(f"|η| = {x['R']:g}: " + ('—' if x['lo'] is None else f"{x['lo'] / (2 - x['R']):.4f}") for x in VC['sector'])
           + "."]
    if VC.get('coef'):                                                  # round 115 (D3dj-1)
        cf, fr, md_, pt = VC['coef'], VC['far'], VC['mid'], VC['path']
        pos = sorted({x['phi'] for x in fr if x['floor'] > 0})
        out[-1] += (" Its first correction inside the disc (referee round 115): `cert_E`/ε³ = Ψ_φ + εΨ¹_φ + O(ε²), Ψ¹_φ(w) = Re((¼ − 2φ² − "
                    "2iφ)w²) − Re(w⁴)/2 (`corner_psi1`), so that on the floor Ψ¹_φ/Ψ_φ = −(¼ + 5φ² + 8φ⁴ − x²Re((1 + 2iφ)⁴)/2)/(1 + 4φ²) ≈ "
                    f"−2φ²; at x = {cf[0]['x']:g} (E/ε³ − Ψ_φ)/(εΨ_φ) in ball arithmetic (Arb, {VC_PREC} bits) at ε = " + ', '.join(   # round 120 (D3do-6)
                        "10⁻" + str(-int(math.log10(float(e)))).translate(str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')) for e in cf[0].get('eps', ('1e-6', '1e-8')))
                    + " against it: "
                    + '; '.join(f"φ = {c['phi']:g}: {c['ball'][0]:.6g}, {c['ball'][1]:.6g} ({c['closed']:.6g})" for c in cf)
                    + " (Ψ_φ in balls; the midpoints kept)"
                    + f" — the expansion holds for θ = εφ ≪ √ε only. Far out, at ε = {VC_FAIL[0]:g}, the true bound on Ψ_φ's floor and its "
                    "least over y from 0 to three times the floor's: " + '; '.join(
                        f"φ = {x['phi']:g}, x = {x['x']:g}: {x['floor']:+.2e}, {x['least']:+.2e}" for x in fr)
                    + (f" — the floor positive from φ = {pos[0]:g}" if pos else '') + (", the least negative throughout" if all(
                        x['least'] < 0 for x in fr) else ", the least NOT negative throughout") + ": there another limit governs. At θ = √εψ, "
                    "z = εx + i√εy, `cert_E`/ε² → Ψ̃_ψ = (1 + 2ψ²)y² − 4ψxy − y⁴/2 + D₊² − P₊², D = x − ½ + y²/2, P = x − ½ − ψ²/2 − ψy "
                    "(`corner_psi_mid`; the referee's): in ball arithmetic at ε = " + ', '.join(
                        "10⁻" + str(-int(math.log10(float(e)))).translate(str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')) for e in md_[0].get('eps', VC_MID_EPS)) + ' '
                    + '; '.join(f"(ψ, x, y) = ({m['psi']:g}, {m['x']:g}, {m['y']:g})" + (", Re σ > 1" if m.get('P', -1) > 0 else '') + ": "
                                + ', '.join(f"{v:.6f}" for v in m['vals']) + f" against {m['lim']:.6f}" for m in md_)
                    + f". Along y = min(2ψx/(1 + 2ψ²), ψ) (`corner_mid_path`) Ψ̃_ψ is " + ("negative" if pt['worst'] < 0 else "NOT negative")
                    + f" at every grid point for {pt['n_psi']} ψ in [{pt['psi'][0]:g}, {pt['psi'][1]:g}] and {pt.get('xmin', 1e-4):g} ≤ x ≤ max(200, 2 + "
                    f"4ψ²) (at most {math.trunc(pt['worst'] * 1e4) / 1e4:.4f}·ψ²min(x, 1)²; inside the disc −4ψ²x²/(1 + 2ψ²) − y⁴/2, far out −ψ² − ψ⁴/2 exactly, to {pt['far']:.1e}). "
                    "The limits match: Ψ_φ(x, φỹ)/φ² (φ → ∞) and Ψ̃_ψ(x, ψỹ)/ψ² (ψ → 0) share Ψ_∞ = ỹ² − 4xỹ + (x − ½)₊(ỹ + 1)², its floor ỹ = "
                    "min(2x, 1) at −4x² and −1 (ε ≪ θ ≪ √ε), and Ψ̃_ψ(x, ŷ/ψ) (ψ → ∞) tends to 2ŷ² − 4xŷ + (x − ½)₊², its floor ŷ = x at −2x² + "
                    "(x − ½)₊² (θ ≫ √ε, the fixed-θ limit). So the corner has four regimes, each limit with an explicit negative path, all four "
                    "proved (Ψ̃_ψ's by cases, Lemma VC in D61's (iii)(b) Status: on y = sx, s = 2ψ/(1 + 2ψ²), P ≤ −ψ²/2 and D₊² < 2ψxy + y⁴/2; on "
                    "y = ψ a convex quadratic in x negative at both ends, then −ψ² − ψ⁴/2; the grid above a check). The four paths are one family "
                    "(referee round 117): Im z = min(κRe z, θ), κ = Im A/|Re A| ≈ 2θ/(ε + 2θ²), A = η²/4 − 1 — the exact quadratic floor capped at "
                    "Im z = θ —, in Ψ̃'s coordinates y = min(sx, ψ) to leading order (the cap y = ψ exactly; at ε > 0 the floor's slope is "
                    "Im A/|Re A|: referee round 120), tending to Ψ_φ's floor as ψ → 0 and to the fixed-θ floor as ψ → ∞; "
                    "Ψ_φ's own floor, y = min(2ψx, ψ) there, gives Ψ̃ = 4ψ²x²(2ψ² − 1) − 8ψ⁴x⁴ inside the disc, negative near the saddle iff ψ < "
                    "1/√2: valid for θ < √(ε/2) (φ = 20 and 25 at ε = 0.001 are ψ = 0.63 and 0.79). Two forms suffice (referee round 120): in "
                    "the path point's coordinates (u, q) = (ε/X, θ²/X) the four regimes are faces of one square — the ψ chart its rays q/u = ψ², "
                    "the φ and (ε/θ, θ/√ε) charts the face q = 0, the matching region its corner u = q = 0 (E ≈ −εθ² − θ⁴/2 − X³/3), the "
                    "(√ε/θ, θ) chart q → ∞, through p = X/θ² —, the junction from b carried by the forms' derivative where their limit vanishes "
                    "(on the real segment the limits vanish, the margin there V's real-axis bound of order ε³), and part VA certifies the path along "
                    "them (below); the saddle's degeneration φ″ → 0 is harmless on this route and no overlap with the straight segment is needed "
                    "(referee rounds 115, 117).")
    else:
        out[-1] += (" So at the corner the straight segment covers arg η below about φ_dom·ε and the floor every φ — a candidate for P on "
                    "1.995 < |η| < 2 (open).")
    return out


VA011_EPS0 = 0.11                                                       # VA011 (referee round 131, D3dz-3): VA at this ε₀, R1 = 1.89 < 1.8927


def part_VA(eps0=None):
    """VA (referee round 120, adopted; verified in round 121, its defect fixed): the last annulus 1.995 ≤ |η| < 2, 0 ≤ arg η ≤
    π/2, certified in ball arithmetic along part VC's path family by `gh_lastann` (forms A, B, the outer piece, the end piece;
    the forms validated against `cert_E` and F_η's definition), with round 121's checks: the true value below the certified bound
    inside every accepted box (`soundness`) and the joint between the inner and outer pieces (`handoff_check`). Returns the
    pieces' records without their timings. eps0: the annulus 2 − ε₀ ≤ |η| < 2 (default gh_lastann's EPS0 = 0.005; part VA011
    passes VA011_EPS0, round 131), the same ε₀ for every piece, the validation and the meta record."""
    import gh_lastann as LA
    LA.ctx.prec = LA.PREC
    e0 = LA.EPS0 if eps0 is None else eps0
    kp = {k: [] for k in ('A-J', 'A-FC', 'B-J', 'B-F')}
    res = [LA.validate(log=lambda s: None, eps0=e0), LA.pieceA('J', eps0=e0, keep=kp['A-J']), LA.pieceA('FC', eps0=e0, keep=kp['A-FC']),
           LA.pieceB('J', eps0=e0, keep=kp['B-J']), LA.pieceB('F', eps0=e0, keep=kp['B-F']), LA.piece_outer(eps0=e0), LA.piece_end()]
    res += [LA.handoff_check(res[5])] + [LA.soundness(k, LA.fA if k[0] == 'A' else LA.fB, kp[k]) for k in kp]
    res[5] = {k: v for k, v in res[5].items() if k != 'joints'}       # round 122: the certified pieces checked, not stored
    res.append(dict(piece='meta', prec=LA.PREC, delta0=LA.DELTA0, eps0=e0))
    _ball()                                                             # back to CERT_BALL_PREC for the other ball users
    return [{k: v for k, v in r.items() if k != 'seconds'} for r in res]


def part_VA011():
    """VA011 (referee round 131, D3dz-3): part VA at ε₀ = VA011_EPS0 = 0.11 — the annulus 2 − ε₀ ≤ |η| < 2 from R1 = 1.89, below D's
    ball-certified radius 1.8927, so that D's branch and bound and VA011 alone cover 1 < |η| < 2 and V is an independent cross-check."""
    return part_VA(VA011_EPS0)


VA_PIECES = ('validate', 'A-J', 'A-FC', 'B-J', 'B-F', 'O', 'end', 'handoff', 'sound A-J', 'sound A-FC', 'sound B-J', 'sound B-F', 'meta')


def va_complete(VA):
    """Round 131 (D3dz-5): whether a VA record set (part VA or VA011) is a whole certificate, checked where its conclusion is printed —
    every piece of VA_PIECES present; the four inner pieces and the outer one with boxes > 0, 0 failed and their largest bound < 0
    (the outer one also with segments > 0); the end piece ok with its bound < 0; the joint with sub-boxes > 0, sampled points > 0,
    one connector per outer η-box, 0 failures, joints checked and the sector tiled; the four soundness checks with points > 0, 0
    failed and no excess (≤ 10⁻³⁰, `soundness`'s own bar); the outer piece's R1 the float 2 − ε₀ of meta's ε₀; the validation with
    points > 0, the forms within `validate_tol(ε₀)` of `cert_E` and F_η's definition and the closed forms within 10⁻⁴⁰. Returns
    (ok, the failed requirements)."""
    import gh_lastann as LA
    by = {r['piece']: r for r in (list(VA) if VA is not None else [])}
    why = [k for k in VA_PIECES if k not in by]
    if why:
        return False, ['missing ' + ', '.join(why)]
    g = lambda k, f, d=None: by[k].get(f, d)
    for k in ('A-J', 'A-FC', 'B-J', 'B-F', 'O'):
        if not (g(k, 'boxes', 0) > 0 and g(k, 'failed', 1) == 0 and g(k, 'worst', 0.0) < 0):
            why.append(k)
    if not g('O', 'segs', 0) > 0:
        why.append('O segments')
    if not (g('end', 'ok', False) and g('end', 'bound', 0.0) < 0):
        why.append('end')
    h = by['handoff']
    if not (h.get('sub_boxes', 0) > 0 and h.get('points', 0) > 0 and h.get('boxes') == g('O', 'boxes') and h.get('failed', 1) == 0
            and h.get('joints', False) and h.get('tiled', False)):
        why.append('handoff')
    for k in ('A-J', 'A-FC', 'B-J', 'B-F'):
        s = by['sound ' + k]
        if not (s.get('points', 0) > 0 and s.get('failed', 1) == 0 and s.get('worst', 1.0) <= 1e-30):
            why.append('sound ' + k)
    e0 = g('meta', 'eps0')
    if e0 is None or g('O', 'R1') != 2 - e0:
        why.append('meta')
    v, tol = by['validate'], (None if e0 is None else LA.validate_tol(e0))
    if not (v.get('counts', {}).get('points', 0) > 0 and tol is not None
            and max(v.get(k, math.inf) for k in ('A_vs_centred', 'A_vs_def', 'B_vs_centred')) < tol
            and max(v.get(k, math.inf) for k in ('A0', 'B0')) < 1e-40):
        why.append('validate')
    return not why, why


def va_seams(VA, V=None, D=None):
    """Round 127 (D3dv-2): the seams the whole-disc conclusion rests on, from the parts themselves — D's branch and bound resolved
    in ball arithmetic (its largest radius R_bb, no open box), V's bands contiguous from R_bb and passing (its radius R_v, as
    `v_lines` computes it), and VA starting exactly there (its outer piece's R1, the float 2 − ε₀). Round 131 (D3dz-5): V's band
    passes only with its end piece rechecked in balls (`ball_end` defaults to 0.0, as `ball_worst` does), and VA must be complete
    (`va_complete`: every piece present, non-empty and passing, R1 = 2 − meta's ε₀, the validation within its ε₀'s bar). Returns
    (R_bb, R_v, R1, ok)."""
    VA = list(VA) if VA is not None else []                           # round 131 (D3dz-8): the raw's numpy array has no truth value
    by = {r['piece']: r for r in VA}
    R1 = by.get('O', {}).get('R1')
    dm = (D or {}).get('dom', {}) if isinstance(D, dict) else {}
    bb = [x for x in dm.get('bb', []) if x.get('open', 1) == 0 and x.get('ball_worst', 0.0) < 0]
    R_bb = max([x['R'] for x in bb] or [0.0])
    R_v = R_bb
    for b in sorted((V or {}).get('bands', []) if isinstance(V, dict) else [], key=lambda b: b['R1']):
        if b['fails'] == 0 and b['end'] < 0 and b.get('ball_worst', 0.0) < 0 and b.get('ball_end', 0.0) < 0 and b['R1'] <= R_v:
            R_v = max(R_v, b['R2'])
    return R_bb, R_v, R1, bool(R_bb > 1.0 and R1 is not None and R_v == R1 and va_complete(VA)[0])


def va011_seam(VA011, D=None):
    """Round 131 (D3dz-3): the seam between D's branch and bound and VA011 — D resolved in ball arithmetic up to R_bb (as in
    `va_seams`), VA011 complete (`va_complete`), and VA011's annulus starting inside D's disc: both the real 2 − ε₀ (the inner pieces
    cover ε ≤ ε₀, ε₀ the double) and the outer piece's double R1 = fl(2 − ε₀) at most R_bb, in rationals. Returns (R_bb, R1, ok)."""
    by = {r['piece']: r for r in (list(VA011) if VA011 is not None else [])}
    R1, e0 = by.get('O', {}).get('R1'), by.get('meta', {}).get('eps0')
    dm = (D or {}).get('dom', {}) if isinstance(D, dict) else {}
    bb = [x for x in dm.get('bb', []) if x.get('open', 1) == 0 and x.get('ball_worst', 0.0) < 0]
    R_bb = max([x['R'] for x in bb] or [0.0])
    ok = (R_bb > 1.0 and R1 is not None and e0 is not None and va_complete(VA011)[0]
          and Fraction(2) - Fraction(e0) <= Fraction(R_bb) and Fraction(R1) <= Fraction(R_bb))
    return R_bb, R1, bool(ok)


def va_lines(VA, V=None, D=None):
    """Part VA's line: what the certificate covers, piece by piece, and its validation; the annulus only when every piece is present,
    non-empty and passing (`va_complete`; round 131), the whole-disc conclusion only when the seams with D's branch and bound and V
    hold (`va_seams`; round 127)."""
    if VA is None or len(VA) == 0:
        return []
    VA = list(VA)
    by = {r['piece']: r for r in VA}
    R_bb, R_v, R1_va, seams = va_seams(VA, V, D)
    complete, why = va_complete(VA)
    ok = (all(r.get('failed', 0) == 0 for r in VA if 'failed' in r) and by.get('end', {}).get('ok', False)
          and by.get('handoff', {}).get('tiled', False) and by.get('handoff', {}).get('joints', False))   # round 122 (D3dq-2)
    v = by.get('validate', {})
    out = ['', '## VA. The last annulus certified in ball arithmetic (referee rounds 120–122, 126, 127, 131)', '',
           "Along part VC's path family — from b = (1 + ε/4)η/2 (a point of part V's chord, ε = 2 − |η|) straight up to the floor Im z = κRe z, "
           "κ = Im A/|Re A|, A = η²/4 − 1, along it to the cap Im z = arg η, along the cap to Re z = 0.05, then a connector and a descent "
           "polyline per η-box to |1 − z| ≤ 0.02 and part V's end piece — the bound is certified negative (`gh_lastann`, python-flint "
           f"at {by.get('meta', {}).get('prec', 128)} bits) by two scaled forms of E/X² — exact on the path up to two one-sided "
           "simplifications ((Re σ − 1)₊ for (|Re σ| − 1)₊ in A, −P₊² dropped in B: upper bounds) —, A in (u, q, ρ) = (ε/X, θ²/X, X) and "
           "B in (u, p, ϑ) = (ε/X, X/θ², θ²), each "
           "box by the mean value form f(ρ) ≤ f(0) + ρ·sup f′ with f(0) in closed form (Ψ̃'s limit) and sup f′ by forward differentiation "
           "in balls: " + '; '.join(
               f"{lab} {by[k]['boxes']} boxes" + (f" ({by[k]['skipped']} outside its domain: ε > ε₀ or Re z > 0.05)" if by[k].get('skipped') else '')
               + f", {by[k]['failed']} failed" for k, lab in (('A-J', 'the junction (form A)'), ('A-FC', 'the floor and the cap (A)'),
                                                            ('B-J', 'the junction (B)'), ('B-F', 'the floor up to arg η = π/2 (B)')) if k in by)
           + (f"; the outer piece {by['O']['boxes']} η-boxes and {by['O']['segs']} segments, {by['O']['failed']} failed, the largest "
              f"rigorous bound {by['O']['worst']:.1e}" if 'O' in by else '')
           + (f"; the end piece {by['end']['bound']:.4f}" if 'end' in by else '')
           + ('' if 'handoff' not in by else
              f"; the joint (`handoff_check`, referee rounds 121, 122), from the pieces certified: at every η-box the connector finite, the "
              f"inner path's end min(κδ₀, arg η) proved inside it in balls on {by['handoff'].get('sub_boxes', 0)} sub-boxes (and sampled at "
              f"{by['handoff']['points']} η: corners, edges, arg η = 0 and π/2, random), the certified pieces chaining the connector from end to "
              f"end and the polyline through every vertex from its start δ₀ + i·y_c to its end inside |1 − z| < 0.02 (round 126), y_c inside the "
              f"connector — {by['handoff']['failed']} failures —, the η-boxes inside the sector with disjoint interiors and exactly its area "
              f"(in rationals)" + ('' if by['handoff'].get('tiled') else " — NOT tiled"))
           + ('' if not any(k.startswith('sound ') for k in by) else
              "; the forms' value below the certified bound at " + ', '.join(
                  f"{by['sound ' + k]['points']}" for k in ('A-J', 'A-FC', 'B-J', 'B-F') if 'sound ' + k in by)
              + f" distinct seeded points and corners (at the expansion variable's upper end and at its lower end where positive; "
              + ', '.join(f"{by['sound ' + k].get('evals', by['sound ' + k]['points'])}" for k in ('A-J', 'A-FC', 'B-J', 'B-F') if 'sound ' + k in by)
              + " point–branch evaluations) in every accepted box of the four inner pieces (`soundness`), "
              f"{sum(by['sound ' + k]['failed'] for k in ('A-J', 'A-FC', 'B-J', 'B-F') if 'sound ' + k in by)} above it, the largest "
              f"excess {max(by['sound ' + k]['worst'] for k in ('A-J', 'A-FC', 'B-J', 'B-F') if 'sound ' + k in by):.1e}")
           + (f". The forms agree with `cert_E` and with F_η's definition to {max(v.get('A_vs_def', 0), v.get('B_vs_centred', 0)):.1e} at "
              f"{v.get('counts', {}).get('points', 0)} path points, the closed forms with the forms at ρ = 0 to "
              f"{max(v.get('A0', 0), v.get('B0', 0)):.1e}" if v else '')
           + (f". So the valley condition holds on {R1_va:g} ≤ |η| < 2 (with η ↦ η̄, −η the whole annulus)"
              + (f" and, with V (to |η| ≤ {R_v:g}) and D's branch and bound (|η| ≤ {R_bb:g}), both in ball arithmetic and meeting VA at "
                 f"exactly its R1, the limit (iii)(b) holds on the whole open disc |η| < 2 (D's domination is the valley condition on the "
                 f"straight path, so Ω ⊆ Ω′)" if seams else " — NOT joined to V and D's branch and bound (the seams do not meet), so the "
                 "limit not on the whole disc")
              + " — computer-assisted (Arb's, python-flint's and CPython's float arithmetic and "
              "`math.nextafter` assumed correct), and since referee round 148 an independent cross-check: the valley condition on "
              "1 < |η| < 2 is proved by hand (review/RESEARCH_note_annulus.md, Theorem A). Round 121's referee verified "
              "the certificate independently — every identity re-derived, the code read line by line, 444 η with 437 720 path points "
              "at 400 bits — and found one defect: in the η-box nearest η = 2 at arg η ≤ π/128, Re A's ball contains 0, κ's ball and "
              "the connector's float ends were NaN, and the connector was silently dropped (round 120's 75 η-boxes and 5126 segments, "
              "a hole between δ₀ + iY(η) and the polyline); fixed by the referee's patch (a NaN or infinite κ gives the connector "
              "[0, θ_hi]; the connector always certified)" + ''.join(
                  f" — that box now split in two, the lower at arg η ≤ {x[3]:.5f} with its connector [{x[4]:g}, {x[5]:.5f}] certified"
                  for x in by['O'].get('accepted', []) if x[0] == by['O'].get('R1') and x[2] == 0.0)
              + "; the counts above are the fixed certificate's (a negative result: round 120's certificate was not a proof as "
              "submitted). Round 122's referee verified the adopted certificate independently again (1493 η with 1 619 161 path points "
              "at 1200 bits, every η-box's corners and edges, 0 defects) and found that the first joint and tiling checks trusted the "
              "connector `outer_box` reported and an area sum; they now rebuild the cover from the pieces certified (above)." if ok and complete
              and 'handoff' in by and by['handoff']['failed'] == 0 and all(by[k]['failed'] == 0 for k in by if k.startswith('sound '))
              else f". NOT complete ({'; '.join(why)}): not a certificate of the annulus." if why and why[0].startswith('missing')
              else ". NOT every piece passes." if not ok or not all(by[k].get('failed', 0) == 0 for k in by)
              else f". NOT complete ({'; '.join(why)}): not a certificate of the annulus.")]
    return out


def va011_lines(VA011, D=None, V=None):
    """Part VA011's line (referee round 131, D3dz-3): part VA at ε₀ = 0.11, piece by piece; the conclusion that D and VA011 alone cover
    1 < |η| < 2 — V an independent cross-check — only when VA011 is complete and passing (`va_complete`) and starts inside D's disc
    resolved in ball arithmetic (`va011_seam`)."""
    if VA011 is None or len(VA011) == 0:
        return []
    import gh_lastann as LA
    VA011 = list(VA011)
    by = {r['piece']: r for r in VA011}
    R_bb, R1, seam = va011_seam(VA011, D)
    complete, why = va_complete(VA011)
    e0 = by.get('meta', {}).get('eps0')
    R_v = va_seams([], V, D)[1] if isinstance(V, dict) and V.get('bands') else None
    v = by.get('validate', {})
    sk = lambda k: f" ({by[k]['skipped']} outside its domain)" if by[k].get('skipped') else ''
    out = ['', f"Part VA again at ε₀ = {'(no meta record)' if e0 is None else f'{e0:g}'} (`part_VA011`, referee round 131, D3dz-3: VA's "
           "code, every piece at this ε₀, nothing else changed)" + ('' if R1 is None else f", from R1 = 2 − ε₀ = {R1:g}") + ": " + '; '.join(
               f"{lab} {by[k]['boxes']} boxes{sk(k)}, {by[k]['failed']} failed" for k, lab in (
                   ('A-J', 'the junction (form A)'), ('A-FC', 'the floor and the cap (A)'), ('B-J', 'the junction (B)'),
                   ('B-F', 'the floor (B)')) if k in by)
           + (f"; the outer piece {by['O']['boxes']} η-boxes and {by['O']['segs']} segments, {by['O']['failed']} failed, the largest "
              f"rigorous bound {by['O']['worst']:.1e}" if 'O' in by else '')
           + (f"; the end piece {by['end']['bound']:.4f}" if 'end' in by else '')
           + (f"; the joint proved on {by['handoff'].get('sub_boxes', 0)} sub-boxes (sampled at {by['handoff']['points']} η), the certified "
              f"pieces chaining connector and polyline, {by['handoff']['failed']} failures, the η-boxes "
              + ("tiling the sector (in rationals)" if by['handoff'].get('tiled') else "NOT tiling the sector") if 'handoff' in by else '')
           + ('' if not all('sound ' + k in by for k in ('A-J', 'A-FC', 'B-J', 'B-F')) else
              "; the forms' value below the certified bound at " + ', '.join(f"{by['sound ' + k]['points']}" for k in ('A-J', 'A-FC', 'B-J', 'B-F'))
              + " distinct points (" + ', '.join(f"{by['sound ' + k]['evals']}" for k in ('A-J', 'A-FC', 'B-J', 'B-F'))
              + f" evaluations), {sum(by['sound ' + k]['failed'] for k in ('A-J', 'A-FC', 'B-J', 'B-F'))} above it")
           + (f"; the forms against `cert_E` and F_η's definition to {max(v.get('A_vs_def', 0), v.get('B_vs_centred', 0)):.1e} at "
              f"{v.get('counts', {}).get('points', 0)} path points — this ε₀'s bar {LA.validate_tol(e0):g} (`validate_tol`: the forms' "
              "midpoints are truncated series, whose tails grow with ε₀; VA's bar 1e-20 is ε₀ = 0.005's) —, the closed forms to "
              f"{max(v.get('A0', 0), v.get('B0', 0)):.1e}" if v and e0 is not None and LA.validate_tol(e0) is not None else '')]
    if seam:
        out[-1] += (f". VA011 certifies the valley condition at every η with 2 − |η| ≤ ε₀ (ε₀ the double {e0:g}; its outer piece from "
                    f"the double {R1:g}, at most that), which starts inside the disc |η| ≤ {R_bb:g} of D's branch and bound in ball "
                    f"arithmetic (both in rationals): so D and VA011 alone cover 1 < |η| < 2 — with |η| ≤ 1, where the domination "
                    "condition is void, and the reflections, the limit (iii)(b) on the whole open disc without V —; V an independent "
                    + (f"cross-check on {R_bb:g} ≤ |η| ≤ {R_v:g}." if R_v is not None and R_v > R_bb else "cross-check.")
                    + " Since referee round 148 D, V, VA and VA011 are all independent cross-checks: the valley condition on 1 < |η| < 2 "
                      "is proved by hand (review/RESEARCH_note_annulus.md, Theorem A).")
    else:
        out[-1] += ('. NOT joined to D' + (f" (NOT complete: {'; '.join(why)})" if not complete else
                                         f" (R1 = {R1:g} is not inside D's ball-certified disc |η| ≤ {R_bb:g})" if R_bb > 1.0 else
                                         " (no branch and bound resolved in ball arithmetic)") + ': no conclusion.')
    return out


def part_DB():
    """DB (referee round 108, D3dc-1): D's branch and bound on 1 < |η| ≤ DOM_BB's radii with every accepted box rechecked in ball arithmetic
    (`dom_bb` with ball); the report takes it as D's branch and bound."""
    return [dict(R=R, n=n, worst=w, open=o, ball_worst=bw) for R in DOM_BB for n, w, o, bw in [dom_bb(R, ball=True)]]


def report(P=None, M=None, L=None, S=None, R=None, LC=None, Q=None, D=None, V=None, DB=None, VC=None, VA=None, VA011=None):
    if DB is not None and isinstance(D, dict) and D.get('dom') and not (     # round 110 (D3de-6): D's own bb when it is in balls
            D['dom'].get('bb') and all('ball_worst' in x for x in D['dom']['bb'])):
        D = dict(D, dom=dict(D['dom'], bb=[dict(x) for x in DB]))
    head = ['# Proposition D‴.61: the Gauss–Hermite Peano kernel (`exp_gh_peano.py`)', '',
            "The K-point Gauss–Hermite rule of N(0, 1): ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] = (K!/(2K)!)ξ^{2K}e^{−ξ²/2}ĝ_K(ξ), ĝ_K the Laplace "
            "transform of the normalised Peano kernel g_K; conjecture P: ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} for complex |η| < 2."]
    body = p_lines(P) + m_lines(M) + l_lines(L) + lc_lines(LC) + q_lines(Q) + d_lines(D, S, R) + v_lines(V, D) + vc_lines(VC) + va_lines(VA, V, D) + va011_lines(VA011, D, V) + s_lines(S) + r_lines(R)
    return '\n'.join(head + body + ['', 'Numbers are produced by `experiments/exp_gh_peano.py`; nothing is edited by hand.', ''])


CERT_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'certificates'))   # round 135 (I4): the stored certificates
CERT_FORMAT = 1


def _cert_split(cert, delta0, lo, hi):
    """The certified pieces of one η-box (`outer_box`'s, in the order certified) split into the connector's (the first pieces, on
    Re z = δ₀, up to the one ending at δ₀ + i·hi) and the polyline's (the rest)."""
    end = complex(delta0, hi)
    for k, (P, Q) in enumerate(cert):
        if not (P.real == delta0 and Q.real == delta0):
            raise AssertionError(('connector piece off Re z = δ₀', k, P, Q))
        if Q == end:
            return cert[:k + 1], cert[k + 1:]
    raise AssertionError('the connector does not reach δ₀ + i·hi')


def export_cert(outdir=None):
    """Round 135 (I4, review M5): export the certified data of D (`dom_bb(1.8927, ball=True)`) and VA011 (part VA's pieces at
    ε₀ = VA011_EPS0) into outdir (default results/certificates), so that `tools/verify_certificates.py` can check Theorem 3's
    certificates from the data alone, without re-running the search. Recomputes both (the same calls as parts DB and VA011, the
    hooks `dom_bb(keep=)`, `bb(keep=)` and `_y_enclosed(keep=)` only collecting what was accepted), asserts that the recomputation
    equals the stored records of results/chaos/gh_peano_raw.npz, and writes
      meta.json               the parameters, counts and the files' formats (floats as float.hex strings)
      D_leaves.f64le.gz        D's accepted leaves: gzip of n × 5 little-endian binary64 (R1, R2, f1, f2, the proposer's bound)
      VA011_inner.jsonl.gz     one JSON line per accepted box of A-J, A-FC, B-J, B-F (its intervals, branches, routes, bounds)
      VA011_outer.jsonl.gz     one JSON line per accepted η-box (box, connector, y_c, polyline vertices, the certified pieces of
                               the connector and of the polyline in the order certified, the handoff's sub-boxes)
      SHA256SUMS               sha256sum's format, over the four files above
    (gzip with mtime 0: the bytes are a function of the data). Nothing else is written; the md and the raw are untouched."""
    import array
    import gzip
    import hashlib
    import json
    import gh_lastann as LA
    outdir = CERT_DIR if outdir is None else outdir
    os.makedirs(outdir, exist_ok=True)
    old = E.load(os.path.join(OUT, 'gh_peano_raw.npz'))
    hx = float.hex
    names = []

    def put(name, data, zipped=True):
        with open(os.path.join(outdir, name), 'wb') as fh:
            if zipped:
                with gzip.GzipFile(filename='', mode='wb', fileobj=fh, mtime=0, compresslevel=9) as g:
                    g.write(data)
            else:
                fh.write(data)
        names.append(name)

    t0 = time.time()
    RD = DOM_BB[-1]
    leaves = []
    n, w, o, bw = dom_bb(RD, ball=True, keep=leaves)
    rec = [x for x in old['DB'] if x['R'] == RD][0]
    assert (n, w, o, bw) == (rec['n'], rec['worst'], rec['open'], rec['ball_worst']) and len(leaves) == n, ('D differs from the raw', n, bw)
    a = array.array('d', [float(v) for L in leaves for v in L])
    if sys.byteorder != 'little':
        a.byteswap()
    put('D_leaves.f64le.gz', a.tobytes())
    tD = time.time() - t0

    t0 = time.time()
    LA.ctx.prec = LA.PREC
    e0, d0 = VA011_EPS0, LA.DELTA0
    kp = {k: [] for k in ('A-J', 'A-FC', 'B-J', 'B-F')}
    inner = [LA.pieceA('J', eps0=e0, keep=kp['A-J']), LA.pieceA('FC', eps0=e0, keep=kp['A-FC']), LA.pieceB('J', eps0=e0, keep=kp['B-J']),
             LA.pieceB('F', eps0=e0, keep=kp['B-F'])]
    O = LA.piece_outer(eps0=e0)
    end = LA.piece_end()
    ho = LA.handoff_check(O)
    by = {r['piece']: r for r in old['VA011']}
    for r in inner:
        s = by[r['piece']]
        assert all(r[k] == s[k] for k in ('boxes', 'skipped', 'failed', 'worst', 'worst_der', 'n_der', 'evals')), ('VA011 differs', r['piece'])
    assert all(O[k] == by['O'][k] for k in ('boxes', 'segs', 'worst', 'failed', 'R1')) and end == by['end'], 'VA011 outer/end differ'
    assert all(ho[k] == by['handoff'][k] for k in ('boxes', 'sub_boxes', 'failed', 'tiled')), 'VA011 handoff differs'
    lines = []
    for name in kp:
        for iv, recs in kp[name]:
            lines.append(json.dumps(dict(piece=name, iv={k: [hx(float(x)) for x in iv[k]] for k in sorted(iv)},
                                         recs=[[br, route, hx(float(v))] for br, route, v in recs])))
    put('VA011_inner.jsonl.gz', ('\n'.join(lines) + '\n').encode())
    lines, nsub, nseg = [], 0, 0
    cpx = lambda z: [hx(z.real), hx(z.imag)]
    for (a_, b_, f1, f2, lo, hi), J in zip(O['accepted'], O['joints']):
        conn, poly = _cert_split(J['cert'], d0, lo, hi)
        sub = []
        assert LA._y_enclosed(a_, b_, f1, f2, lo, hi, d0, keep=sub) == len(sub)
        nsub, nseg = nsub + len(sub), nseg + len(conn) + len(poly)
        lines.append(json.dumps(dict(box=[hx(x) for x in (a_, b_, f1, f2)], lo=hx(lo), hi=hx(hi), yc=hx(J['yc']),
                                     pts=[cpx(z) for z in J['pts']], conn=[cpx(P) + cpx(Q) for P, Q in conn],
                                     poly=[cpx(P) + cpx(Q) for P, Q in poly], sub=[[hx(x) for x in s] for s in sub])))
    assert nsub == ho['sub_boxes'] and nseg == O['segs']
    put('VA011_outer.jsonl.gz', ('\n'.join(lines) + '\n').encode())
    tVA = time.time() - t0
    _ball()                                                             # back to CERT_BALL_PREC for the other ball users

    meta = dict(
        format=CERT_FORMAT,
        source='experiments/exp_gh_peano.py --export-cert: dom_bb(1.8927, ball=True) and part VA011 (experiments/gh_lastann.py), '
               'recomputed and checked equal to results/chaos/gh_peano_raw.npz',
        floats='every float is a Python float.hex string (exact binary64), except in D_leaves (raw binary64)',
        D=dict(file='D_leaves.f64le.gz', R1=hx(1.0), R2=hx(RD), f1=hx(0.0), f2=hx(math.pi / 2), leaves=n, prec=CERT_BALL_PREC,
               columns=['R1', 'R2', 'f1', 'f2', 'H (the proposer\'s rigorous bound, informational)'],
               layout='gzip of leaves x 5 little-endian IEEE-754 binary64, row-major, in the order accepted',
               proposer_worst=hx(bw)),
        VA011=dict(files=['VA011_inner.jsonl.gz', 'VA011_outer.jsonl.gz'], eps0=hx(e0), delta0=hx(d0), end_r=hx(LA.END_R),
                   QA=hx(LA.QA), PB=hx(LA.PB), UJ=hx(LA.UJ), R1=hx(O['R1']), R2=hx(2.0), prec=LA.PREC, tail_floor=hx(LA.TAIL_FLOOR),
                   terms=dict(S=16, Cv=16, Lam=30, H1=26),
                   pieces={r['piece']: dict(boxes=r['boxes'], skipped=r['skipped'], proposer_worst=hx(r['worst'])) for r in inner},
                   outer=dict(boxes=O['boxes'], segments=O['segs'], proposer_worst=hx(O['worst'])), handoff_sub_boxes=nsub,
                   end=dict(R2=hx(2.0), eps_e=hx(LA.END_R), proposer_bound=hx(end['bound']))))
    put('meta.json', (json.dumps(meta, indent=1) + '\n').encode(), zipped=False)
    sums = ''.join(f"{hashlib.sha256(open(os.path.join(outdir, k), 'rb').read()).hexdigest()}  {k}\n" for k in sorted(names))
    open(os.path.join(outdir, 'SHA256SUMS'), 'w').write(sums)
    print(f"exported to {outdir}: D {n} leaves [{tD:.0f}s], VA011 {sum(len(v) for v in kp.values())} inner boxes, {O['boxes']} η-boxes, "
          f"{nseg} segments, {nsub} sub-boxes [{tVA:.0f}s]", flush=True)


KEYS = ('P', 'M', 'L', 'S', 'R', 'LC', 'Q', 'D', 'V', 'DB', 'VC', 'VA', 'VA011')


def main():
    warnings.filterwarnings('ignore')
    if '--export-cert' in sys.argv:                                     # round 135 (I4): the certificates' data, nothing else written
        i = sys.argv.index('--export-cert')
        export_cert(sys.argv[i + 1] if i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith('--') else None)
        return
    raw = os.path.join(OUT, 'gh_peano_raw.npz')
    md = os.path.join(OUT, 'gh_peano.md')
    old = E.load(raw) if os.path.exists(raw) else {}
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in KEYS]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(KEYS)
    res = dict(old)
    t0 = time.time()
    a6 = lambda: E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))['A6']
    for k, fn in (('P', part_P), ('M', lambda: part_M(res['P'])), ('L', part_L), ('S', part_S), ('R', lambda: part_R(a6())), ('LC', part_LC), ('Q', part_Q), ('D', part_D), ('V', part_V), ('DB', part_DB), ('VC', part_VC), ('VA', part_VA), ('VA011', part_VA011)):
        if k in parts:
            res[k] = fn()
            E.save(raw, **res)
            print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in KEYS]))


if __name__ == '__main__':
    main()
