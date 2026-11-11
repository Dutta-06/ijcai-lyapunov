"""
Wendland kernels (compactly supported basis functions) as Gram classes for theory3e_kernels.

phi_0(r) = (1 - r)_+^3 (3 r + 1)  (n = 1; Wendland psi_{1,1}, C^2)  and  (1 - r)_+^4 (4 r + 1)  (n = 2, 3; Wendland psi_{3,1}, C^2, positive definite in R^n, n <= 3), support radius 1.
Its Gram Gamma(d) = <phi(.) phi(. - d e)> / ||phi||^2 has support [0, 2] (raw units) and is, on the pieces [0,1], [1,2] (n = 1) and [0,1/2], [1/2,1], [1,2] (n = 3), an EXACT rational
polynomial (n = 1) or polynomial / d (n = 3), computed here in rational arithmetic (the integrals are polynomial in d on each piece: fitted by exact interpolation of exactly integrated values
and verified at further points).  n = 2: the Gram is not piecewise polynomial (arccos terms); it is computed by floating quadrature on a grid and only the support, the curvature
(-Gamma''(0) = ||d_1 phi||^2/||phi||^2, exact rational) and the qualitative facts (Gamma > 0, Gamma' < 0 on (0, 2): Gamma is the autocorrelation of a radially decreasing function) are used.
"""
import math
import os
import sys
from fractions import Fraction

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from flint import arb, fmpq

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import theory3e_kernels as TK     # noqa: E402


# ---- rational polynomials, ascending coefficients
def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fraction(0)) + (b[i] if i < len(b) else Fraction(0)) for i in range(m)]


def pneg(a):
    return [-c for c in a]


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pscale(a, c):
    return [c * x for x in a]


def pcomp_lin(a, alpha, beta):
    """a(alpha x + beta)"""
    out = [Fraction(0)]
    base = [Fraction(1)]
    lin = [Fraction(beta), Fraction(alpha)]
    for c in a:
        out = padd(out, pscale(base, c))
        base = pmul(base, lin)
    return out


def pint(a):
    return [Fraction(0)] + [c / (i + 1) for i, c in enumerate(a)]


def peval(a, x):
    acc = Fraction(0)
    for c in reversed(a):
        acc = acc * x + c
    return acc


def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [Fraction(0)]


def definite(a, lo_, hi_):
    A = pint(a)
    return peval(A, hi_) - peval(A, lo_)


def phi_poly(n):
    """phi_0 as a polynomial in r on [0, 1]"""
    if n == 1:
        return pmul(pmul(pmul([Fraction(1), Fraction(-1)], [Fraction(1), Fraction(-1)]), [Fraction(1), Fraction(-1)]), [Fraction(1), Fraction(3)])
    q = [Fraction(1), Fraction(-1)]
    q4 = pmul(pmul(q, q), pmul(q, q))
    return pmul(q4, [Fraction(1), Fraction(4)])


def gram_exact_value(n, d):
    """Gamma_raw(d) for a rational d in [0, 2] (n = 1 or 3), exactly."""
    d = Fraction(d)
    p = phi_poly(n)
    if n == 1:
        norm = 2 * definite(pmul(p, p), Fraction(0), Fraction(1))
        if d >= 2:
            return Fraction(0)
        if d <= 1:
            # x in [d-1, 0]: |x| = -x, |x-d| = d - x ;  [0, d]: x, d - x ;  [d, 1]: x, x - d
            f1 = pmul(pcomp_lin(p, -1, 0), pcomp_lin(p, -1, d))
            f2 = pmul(p, pcomp_lin(p, -1, d))
            f3 = pmul(p, pcomp_lin(p, 1, -d))
            val = definite(f1, d - 1, Fraction(0)) + definite(f2, Fraction(0), d) + definite(f3, d, Fraction(1))
        else:
            f = pmul(p, pcomp_lin(p, -1, d))
            val = definite(f, d - 1, Fraction(1))
        return val / norm
    # n = 3:  Gamma(d) = I(d) / (2 d int a^2 phi^2),  I(d) = int_a a phi(a) [Q(hi) - Q(lo)] da,  Q(t) = int_0^t b phi(b) db  (t <= 1)
    a_phi = pmul([Fraction(0), Fraction(1)], p)
    Q = pint(a_phi)
    Q1 = peval(Q, Fraction(1))
    nrm = 2 * definite(pmul([Fraction(0), Fraction(0), Fraction(1)], pmul(p, p)), Fraction(0), Fraction(1))
    if d == 0:
        return Fraction(1)
    if d >= 2:
        return Fraction(0)

    def piece(lo_lin, hi_lin, a0, a1):
        # lo_lin, hi_lin = (alpha, beta) giving lo = alpha a + beta ; hi likewise (hi = 'one' -> constant Q(1))
        Qlo = pcomp_lin(Q, lo_lin[0], lo_lin[1])
        Qhi = [Q1] if hi_lin is None else pcomp_lin(Q, hi_lin[0], hi_lin[1])
        integrand = pmul(a_phi, padd(Qhi, pneg(Qlo)))
        return definite(integrand, a0, a1)
    if d < Fraction(1, 2):
        I = piece((-1, d), (1, d), Fraction(0), d) + piece((1, -d), (1, d), d, 1 - d) + piece((1, -d), None, 1 - d, Fraction(1))
    elif d <= 1:
        I = piece((-1, d), (1, d), Fraction(0), 1 - d) + piece((-1, d), None, 1 - d, d) + piece((1, -d), None, d, Fraction(1))
    else:
        I = piece((-1, d), None, d - 1, Fraction(1))
    return I / (nrm * d)


def fit_piece(n, d0, d1, deg):
    """Exact polynomial N(d) on [d0, d1] with Gamma = N(d) (n = 1) or N(d)/d (n = 3), by interpolation at deg+1 interior rational points; verified at 3 more points."""
    pts = [d0 + (d1 - d0) * Fraction(2 * i + 1, 2 * (deg + 1)) for i in range(deg + 1)]
    ys = [gram_exact_value(n, x) * (x if n == 3 else 1) for x in pts]
    # Lagrange -> monomial coefficients by solving the Vandermonde system exactly
    m = deg + 1
    A = [[x ** j for j in range(m)] + [y] for x, y in zip(pts, ys)]
    for c in range(m):
        piv = next(r for r in range(c, m) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [vr - f * vc for vr, vc in zip(A[r], A[c])]
    coef = [A[i][m] for i in range(m)]
    for x in (d0 + (d1 - d0) * Fraction(1, 7), d0 + (d1 - d0) * Fraction(3, 5), d0 + (d1 - d0) * Fraction(6, 7)):
        assert peval(coef, x) == gram_exact_value(n, x) * (x if n == 3 else 1), 'piece is not polynomial (degree too low)'
    return coef


class Wendland(TK.Gram):
    def __init__(self, n):
        assert n in (1, 3)
        self.n = n
        self.name = 'wendland%d' % n
        if n == 1:
            self.bps = [Fraction(0), Fraction(1), Fraction(2)]
            deg = 12
        else:
            self.bps = [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2)]
            deg = 16
        self.pieces = []
        for a, b in zip(self.bps[:-1], self.bps[1:]):
            c = fit_piece(n, a, b, deg)
            while len(c) > 1 and c[-1] == 0:
                c = c[:-1]
            m = 1 if n == 3 else 0
            if m == 1 and a == 0:                       # N(d) = d N1(d) on the first piece: divide exactly
                assert c[0] == 0
                c, m = c[1:], 0
            self.pieces.append((c, m))
        # the same polynomials expanded about the RIGHT endpoint b of their piece (exact): N(b + t); the ball evaluation uses t = x - b (no cancellation between large monomial coefficients near the support edge)
        self.shifted = [(pcomp_lin(c, 1, b), m) for (c, m), b in zip(self.pieces, self.bps[1:])]
        c0, m0 = self.pieces[0]
        assert m0 == 0 and c0[0] == 1 and c0[1] == 0
        self.lam0 = float(-2 * c0[2])
        self.support = 2.0
        super().__init__()
        self.supp_eq = self.support / self.s          # support of the sigma-equivalent Gram

    # raw value of order 0/1/2 on one piece (arb or float)
    def _piece(self, i, order, x):
        if isinstance(x, arb):
            c, m = self.shifted[i]
            xt = x - arb(fmpq(self.bps[i + 1].numerator, self.bps[i + 1].denominator))     # t = x - b
        else:
            c, m = self.pieces[i]
            xt = x
        N = [c, pder(c), pder(pder(c))]
        conv = (lambda q: arb(fmpq(q.numerator, q.denominator))) if isinstance(x, arb) else float

        def hor(C):
            acc = conv(Fraction(0))
            for q in reversed(C):
                acc = acc * xt + conv(q)
            return acc
        if m == 0:
            return hor(N[order])
        n0, n1, n2 = hor(N[0]), hor(N[1]), hor(N[2])
        if order == 0:
            return n0 / x
        if order == 1:
            return n1 / x - n0 / (x * x)
        return n2 / x - 2 * n1 / (x * x) + 2 * n0 / (x * x * x)

    def raw(self, order, x):
        if isinstance(x, arb):
            lo_, hi_ = float(x.lower().mid()), float(x.upper().mid())
            lo_, hi_ = float(x.mid()) - float(x.rad()), float(x.mid()) + float(x.rad())
            if lo_ >= 2.0:
                return arb(0)
            res = None
            for i, (a, b) in enumerate(zip(self.bps[:-1], self.bps[1:])):
                a_, b_ = float(a), float(b)
                if hi_ < a_ or lo_ > b_:
                    continue
                xs = x if (lo_ >= a_ and hi_ <= b_) else arb((max(lo_, a_) + min(hi_, b_)) / 2, (min(hi_, b_) - max(lo_, a_)) / 2 * (1 + 1e-12) + 1e-300)
                v = self._piece(i, order, xs)
                res = v if res is None else res.union(v)
            if hi_ > 2.0:
                res = res.union(arb(0))
            return res
        return float(self.rawf(order, np.asarray(x, float)))

    def rawf(self, order, x):
        """vectorised float evaluation (piecewise polynomial / d^m): float coefficients are cached"""
        if not hasattr(self, '_fc'):
            self._fc = []
            for (c, m) in self.pieces:
                cs_ = [c, pder(c), pder(pder(c))]
                self._fc.append((m, [np.array([float(q) for q in reversed(v)]) for v in cs_]))
        x = np.asarray(x, float)
        out = np.zeros_like(x)
        for i, (a, b) in enumerate(zip(self.bps[:-1], self.bps[1:])):
            msk = (x >= float(a)) & (x < float(b))
            if not msk.any():
                continue
            xv = x[msk]
            m, (n0, n1, n2) = self._fc[i]
            p0, p1, p2 = np.polyval(n0, xv), np.polyval(n1, xv), np.polyval(n2, xv)
            if m == 0:
                out[msk] = (p0, p1, p2)[order]
            elif order == 0:
                out[msk] = p0 / xv
            elif order == 1:
                out[msk] = p1 / xv - p0 / xv ** 2
            else:
                out[msk] = p2 / xv - 2 * p1 / xv ** 2 + 2 * p0 / xv ** 3
        return out

    def ptilde_raw(self, x):
        # psi/d = -Gamma'(d)/d on the first piece: Gamma = c0 + c2 d^2 + ... (c1 = 0): -Gamma'/d = -(2 c2 + 3 c3 d + ...) : polynomial, regular at 0
        c, m = self.pieces[0]
        d = pder(c)                                       # Gamma'
        assert d[0] == 0
        q = [-v for v in d[1:]]                           # -Gamma'/d
        if isinstance(x, arb):
            if float(x.mid()) - float(x.rad()) <= 1.0:
                acc = arb(0)
                for v in reversed(q):
                    acc = acc * x + arb(fmpq(v.numerator, v.denominator))
                return acc
            return -self.raw(1, x) / x
        return -self.raw(1, x) / x if x > 1.0 else float(peval(q, Fraction(x).limit_denominator(10 ** 12)))

    def cell_up(self, a, b):
        ball = TK.cell(a, b)
        v0 = TK.up(abs(self.G(ball)))
        v1 = TK.up(abs(self.G1(ball)))
        v2 = TK.up(abs(self.G2(ball)))
        return v0, v1, v2, max(v2, v1 / a)


def lam0_n2():
    """-Gamma''(0) = ||d_1 phi||^2 / ||phi||^2 = (1/n) int phi'(r)^2 r^{n-1} dr / int phi^2 r^{n-1} dr for n = 2 (exact rational), phi = (1-r)^4 (4r+1)."""
    p = phi_poly(2)
    dp = pder(p)
    num = definite(pmul(dp, dp), Fraction(0), Fraction(1)) if False else definite(pmul(pmul(dp, dp), [Fraction(0), Fraction(1)]), Fraction(0), Fraction(1))
    den = definite(pmul(pmul(p, p), [Fraction(0), Fraction(1)]), Fraction(0), Fraction(1))
    return num / den / 2
