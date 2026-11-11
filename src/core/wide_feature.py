"""
Proposition D‴.58: the transit under width mismatch — a feature wider than the kernels as a Gaussian quadrature of its mixing law
====================================================================================================================================

A Gaussian feature f* = Aφ_{σ*}(· − c*) wider than the kernels φ_σ(x) = exp(−x²/(2σ²)) (σ* > σ) is the kernel smoothed by its mixing
law: f* = A(σ*/σ)·φ_σ * N(c*, τ²), τ² = σ*² − σ². A K-kernel sum Σ_k a_kφ_σ(· − c* − y_k) is a K-point quadrature of that law and its
L²(ℝ) error the rule's worst-case error in the reproducing kernel Hilbert space of ⟨φ_σ(· − y), φ_σ(· − y′)⟩. The half-squeeze frame —
the dilation to the geometric-mean width (σσ*)^{1/2}, the Bargmann transform, Mehler's formula at ρ = t = tanh(r/2), r = ln(σ*/σ) —
makes the error of the Gauss–Hermite rule of the reweighted law N(0, τ²σ/σ*) (`half_squeeze`) the exact series
‖f*‖²(1 − t²)^{1/2}Σ_{n ≥ 2K}(tⁿ/n!)D_{n,K}², D_{n,K} = E He_n(√2Z) − Σ_k w_k He_n(√2z_k) (`gh_defect`, `mehler_series`): the rate t per
kernel, against λ = tanh r for the Gauss–Hermite rule of the mixing law itself (`moment_matching`).

This module: the rules, the L²(ℝ) error of any K-kernel sum against the feature by a trapezoid rule on a fine grid (spectrally accurate,
free of the cancellation of ‖f*‖² − 2aᵀb + aᵀGa; `l2_error`), the exact D_{n,K} in integer arithmetic, the series, least-squares
amplitudes for given nodes (`ls_amplitudes`) and the optimal K-kernel fit by variable projection (`optimal_fit`).
"""
from fractions import Fraction
from math import factorial, lgamma, log

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import least_squares


def gh(K):
    """The K-point Gauss–Hermite rule of N(0, 1): nodes z_k and weights w_k (Σw_k = 1)."""
    z, w = hermegauss(K)
    return z, w / w.sum()


def params(sig, sig_star):
    """τ, r = ln(σ*/σ), t = tanh(r/2) = (σ* − σ)/(σ* + σ), λ = tanh r = τ²/(σ² + σ*²)."""
    tau = np.sqrt(sig_star ** 2 - sig ** 2)
    return dict(tau=tau, r=np.log(sig_star / sig), t=(sig_star - sig) / (sig_star + sig), lam=tau ** 2 / (sig ** 2 + sig_star ** 2))


def half_squeeze(K, sig, sig_star, A=1.0):
    """The half-squeeze rule: nodes y_k = τ(σ/σ*)^{1/2}z_k (the Gauss–Hermite nodes of N(0, τ²σ/σ*)), amplitudes
    a_k = A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*))), relative to the feature's centre."""
    z, w = gh(K)
    tau = params(sig, sig_star)['tau']
    y = tau * np.sqrt(sig / sig_star) * z
    return y, A * np.sqrt(sig_star / sig) * w * np.exp(y ** 2 / (2 * sig * (sig + sig_star)))


def moment_matching(K, sig, sig_star, A=1.0):
    """The Gauss–Hermite rule of the mixing law N(0, τ²): nodes τz_k, amplitudes A(σ*/σ)w_k."""
    z, w = gh(K)
    return params(sig, sig_star)['tau'] * z, A * (sig_star / sig) * w


def scaled_rule(K, sig, sig_star, v, A=1.0):
    """The Gauss–Hermite rule of N(0, v) for the mixing law reweighted to it: nodes y_k = √v z_k, amplitudes
    A(σ*/σ)w_k(v/τ²)^{1/2}exp(y_k²(1/v − 1/τ²)/2). v = τ² is moment matching, τ²σ/σ* the half-squeeze rule; the natural Fock scaling
    divides the half-squeeze variance by cosh r, the Mercer-zero scaling is τ²(1 + 2τ²/σ²)^{−1/2} (referee round 48's comparisons)."""
    z, w = gh(K)
    tau = params(sig, sig_star)['tau']
    y = np.sqrt(v) * z
    return y, A * (sig_star / sig) * w * np.sqrt(v / tau ** 2) * np.exp(y ** 2 * (1 / v - 1 / tau ** 2) / 2)


def scalings(sig, sig_star):
    """The variances of `scaled_rule`'s named scalings: moment matching, half-squeeze, Fock, Mercer-zero."""
    P = params(sig, sig_star)
    tau2 = P['tau'] ** 2
    vhs = tau2 * sig / sig_star
    return dict(mm=tau2, hs=vhs, fock=vhs / np.cosh(P['r']), mercer=tau2 / np.sqrt(1 + 2 * tau2 / sig ** 2))


def best_scale(K, sig, sig_star, A=1.0):
    """The variance v of least L² error for `scaled_rule` (bounded scalar search in log v around the half-squeeze variance) → (v, error):
    an upper bound on the optimal K-kernel fit's error by one parameter."""
    from scipy.optimize import minimize_scalar
    vhs = scalings(sig, sig_star)['hs']
    f = lambda lv: l2_error(*scaled_rule(K, sig, sig_star, vhs * np.exp(lv), A), sig, sig_star, A)
    res = minimize_scalar(f, bounds=(-1.0, 1.0), method='bounded', options=dict(xatol=1e-10))
    return float(vhs * np.exp(res.x)), float(res.fun)


def narrow_rule(K, sig, sig_star, A=1.0):
    """(v): the half-squeeze frame continued to a feature narrower than the kernels (σ* < σ, r < 0): the nodes iη_k, η_k =
    (2σ²|sinh r|)^{1/2}z_k, in conjugate pairs, the K-term sum real — e^{−x²/(2σ²)}Σ_k b_k cos(η_kx/σ²), b_k = A(σ*/σ)^{1/2}w_k
    exp(η_k²/(2σ²) − η_k²/(2σ(σ + σ*))) (Derevianko–Plonka's Hermite-zero cosine sums) → (η, b)."""
    z, w = gh(K)
    eta = np.sqrt(2 * sig ** 2 * abs(np.sinh(np.log(sig_star / sig)))) * z
    return eta, A * np.sqrt(sig_star / sig) * w * np.exp(eta ** 2 / (2 * sig ** 2) - eta ** 2 / (2 * sig * (sig + sig_star)))


def narrow_error(eta, b, sig, sig_star, A=1.0, per=400):
    """The relative L² error of `narrow_rule`'s cosine sum against Aφ_{σ*}, by the trapezoid rule on ±40σ* with spacing σ*/per."""
    h = sig_star / per
    x = np.arange(-40 * per, 40 * per + 1) * h
    fK = np.exp(-x ** 2 / (2 * sig ** 2)) * (np.cos(np.outer(x, eta) / sig ** 2) @ b)
    return float(np.sqrt(h * np.sum((A * np.exp(-x ** 2 / (2 * sig_star ** 2)) - fK) ** 2)) / (A * np.sqrt(np.sqrt(np.pi) * sig_star)))


def grid(sig, sig_star, span=14.0, per_sig=10):
    """A trapezoid grid on [−span·σ*, span·σ*] with spacing σ/per_sig (the integrands are Gaussians of width ≥ σ/√2: the rule's error is
    of order exp(−π²(σ·per_sig)²/σ²) — far below double precision)."""
    h = sig / per_sig
    n = int(np.ceil(span * sig_star / h))
    return np.arange(-n, n + 1) * h, h


def residual(y, a, sig, sig_star, A=1.0, x=None):
    if x is None:
        x, _ = grid(sig, sig_star)
    return A * np.exp(-x ** 2 / (2 * sig_star ** 2)) - np.exp(-(x[:, None] - np.asarray(y)[None, :]) ** 2 / (2 * sig ** 2)) @ np.asarray(a)


def l2_error(y, a, sig, sig_star, A=1.0, rel=True):
    """‖f* − Σ_k a_kφ_σ(· − y_k)‖ on ℝ (relative to ‖f*‖ = A(π^{1/2}σ*)^{1/2} if rel), by the trapezoid rule on `grid`."""
    x, h = grid(sig, sig_star)
    e = np.sqrt(h * np.sum(residual(y, a, sig, sig_star, A, x) ** 2))
    return e / (A * np.sqrt(np.sqrt(np.pi) * sig_star)) if rel else e


def ls_amplitudes(y, sig, sig_star, A=1.0):
    """The least-squares amplitudes for the nodes y (the L² projection, on `grid`)."""
    x, _ = grid(sig, sig_star)
    P = np.exp(-(x[:, None] - np.asarray(y)[None, :]) ** 2 / (2 * sig ** 2))
    return np.linalg.lstsq(P, A * np.exp(-x ** 2 / (2 * sig_star ** 2)), rcond=None)[0]


def optimal_fit(y0, sig, sig_star, A=1.0, symmetric=True):
    """The K-kernel fit of least L² error by variable projection from the nodes y0 (amplitudes by least squares for each node set);
    symmetric: the nodes kept symmetric (±y, and 0 for odd K). Returns (nodes, amplitudes, relative error)."""
    x, h = grid(sig, sig_star)
    f = A * np.exp(-x ** 2 / (2 * sig_star ** 2))
    y0 = np.sort(np.asarray(y0, float))
    K = y0.size
    if symmetric:
        p0 = y0[y0 > 1e-12]
        full = lambda p: np.r_[-p[::-1], [0.0] * (K % 2), p]
    else:
        p0, full = y0, (lambda p: p)

    def res(p):
        P = np.exp(-(x[:, None] - full(p)[None, :]) ** 2 / (2 * sig ** 2))
        a = np.linalg.lstsq(P, f, rcond=None)[0]
        return (f - P @ a) * np.sqrt(h)
    if p0.size:
        sol = least_squares(res, p0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        p = sol.x
    else:
        p = p0
    y = full(p)
    a = ls_amplitudes(y, sig, sig_star, A)
    return y, a, l2_error(y, a, sig, sig_star, A)


def _he(n):
    """The integer coefficients of He_n in powers of x."""
    a, b = [1], [0, 1]
    if n == 0:
        return a
    for k in range(1, n):
        c = [0] + b
        for i, v in enumerate(a):
            c[i] -= k * v
        a, b = b, c
    return b


def _gauss_moment(m):
    if m % 2:
        return 0
    r = 1
    for k in range(1, m, 2):
        r *= k
    return r


def gh_defect(n, K):
    """D_{n,K} = E He_n(√2Z) − Σ_k w_k He_n(√2z_k), exactly (an integer): the rule integrates the remainder of He_n(√2z) modulo He_K(z)
    — a polynomial of degree < K — exactly, and He_K vanishes at its nodes. Zero for odd n and for n < 2K."""
    if n % 2:
        return 0
    c = _he(n)
    p = [c[m] * 2 ** (m // 2) if m % 2 == 0 else 0 for m in range(len(c))]
    q = _he(K)
    for d in range(len(p) - 1, K - 1, -1):
        f = p[d]
        if f:
            for i in range(K + 1):
                p[d - K + i] -= f * q[i]
    QK = sum(Fraction(p[m]) * _gauss_moment(m) for m in range(min(K, len(p))))
    D = _gauss_moment(n) - QK
    assert D.denominator == 1
    return int(D)


def gh_defects(K, n_max):
    """D_{n,K} for n = 0 … n_max at once (integers): Q_K(He_m) is the He_0 coefficient of He_m modulo He_K — the remainders by
    He_{m+1} = zHe_m − mHe_{m−1} and z·He_h = He_{h+1} + hHe_{h−1} (He_K ≡ 0), O(K) per degree —, and by the multiplication theorem
    He_n(√2z) = Σ_i n!/(2^i i!(n − 2i)!)·2^{(n−2i)/2}He_{n−2i}(z) with E He_m(Z) = 0 for m ≥ 1."""
    q, prev, c = [], None, [1] + [0] * (K - 1)
    for m in range(n_max + 1):
        q.append(c[0])
        zc = [0] * K
        for h, v in enumerate(c):
            if v:
                if h + 1 < K:
                    zc[h + 1] += v
                if h:
                    zc[h - 1] += h * v
        c, prev = [zc[i] - (m * prev[i] if prev is not None else 0) for i in range(K)], c
    return [0 if n % 2 else -sum(factorial(n) // (2 ** i * factorial(i) * factorial(n - 2 * i)) * 2 ** (n // 2 - i) * q[n - 2 * i]
                                 for i in range(n // 2)) for n in range(n_max + 1)]


def tail_ratio(K, t, n_terms=200):
    """The series over its first term, Σ_{j ≥ 0} t^{2j}(2K)!/(2K + 2j)!·(D_{2K+2j,K}/D_{2K,K})² — each term tending to binom(1/2, j)²t^{2j}
    as K → ∞ (the ratios' leading coefficients binom(1/2, j)(−2)^j), so conjecturally the ratio to ₂F₁(−½, −½; 1; t²)."""
    D = gh_defects(K, 2 * K + 2 * n_terms)
    l0 = log(abs(D[2 * K]))
    return float(sum(np.exp(2 * j * log(t) + 2 * (log(abs(D[2 * K + 2 * j])) - l0) + lgamma(2 * K + 1) - lgamma(2 * K + 2 * j + 1))
                     for j in range(n_terms) if D[2 * K + 2 * j]))


def mehler_series(K, t, n_terms=200):
    """(1 − t²)^{1/2}Σ_{n ≥ 2K}(tⁿ/n!)D_{n,K}² (the half-squeeze rule's squared relative L² error) and its first term; the terms by
    logarithms (D_{n,K} grows like a factorial)."""
    terms = []
    for n in range(2 * K, 2 * K + 2 * n_terms, 2):
        D = gh_defect(n, K)
        terms.append(0.0 if D == 0 else np.exp(n * log(t) + 2 * log(abs(D)) - lgamma(n + 1)))
    pre = np.sqrt(1 - t ** 2)
    return pre * float(np.sum(terms)), pre * terms[0]


def first_term(K, t):
    """The series' first term (1 − t²)^{1/2}t^{2K}4^K K!²/(2K)!: the half-squeeze rule's squared relative error to leading order,
    ≈ (1 − t²)^{1/2}(πK)^{1/2}t^{2K}."""
    return np.sqrt(1 - t ** 2) * np.exp(2 * K * log(t) + K * log(4) + 2 * lgamma(K + 1) - lgamma(2 * K + 1))
