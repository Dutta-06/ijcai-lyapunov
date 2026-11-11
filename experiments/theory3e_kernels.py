"""
Round theory3 / direction E: the coupled entry layer for a GENERAL Gram -- abstract hypotheses, explicit conditions, and their evaluation in ball arithmetic.

The continuum network has kernels phi(x; c) = phi_0(|x - c|) and its fit F = (1/2)<(y - sum w_i phi_i)^2> depends on the kernel only through the Gram
<phi_c phi_c'> = Z Gamma(|c - c'|), Gamma(0) = 1 (the autocorrelation of phi_0).  Proposition G (NOTE.md of the round; the Gaussian Proposition S-A is the case
Gamma(d) = exp(-d^2/4)) states the box conditions (L), (C), (N) of the entry layer through the following numbers of Gamma only (lengths in "sigma-equivalents", the
unit in which -Gamma''(0) = 1/2, as for the Gaussian of width sigma):

  tail envelopes (the coherence of the Gram)   m0(l) = sup_{d>=l}|Gamma|,  m1(l) = sup_{d>=l}|Gamma'|,  m2(l) = sup_{d>=l}||Hess Gamma(|x|)|| = sup max(|Gamma''|, |Gamma'|/d)
                                               (n = 1: m2 = sup |Gamma''|), M1 = m1(0), M2 = m2(0);
  own constants at a radius x                  g(x) = min_{[0,x]} Gamma,  psi(x) = -Gamma'(x),  lampsi(x) = inf_{(0,x]} psi(d)/d,  psihat(x) = sup_{[0,x]} psi,
                                               rho(x) = inf_{[0,x]} (-Gamma'').

`box(...)` evaluates the conditions in ball arithmetic (python-flint): every kernel number is a DIRECTED bound (upper for the tails, lower for the own minima), turned into an exact
ball, and the arithmetic of the conditions is done on balls, so `ok` means the strict inequalities are certainly true.  The kernels: Gaussian; Matern-nu Grams (closed form
polynomial x exp for half-integer nu, Bessel K for integer nu >= 2), which are the Grams of the exponential / Matern basis functions in R^n (nu_Gram = 2 nu_basis + n/2); the Wendland
kernels (compact Gram) are in theory3e_wendland.py.  Directed bounds near d = 0 for the Bessel kernels use monotone brackets (see the methods).
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

L_MAX = 40.0          # the tails are enclosed on [l, L_MAX] by cells; beyond L_MAX the three envelopes are assumed non-increasing (hand argument in the note; every kernel is < 1e-9 there)
H_TAB = 0.01          # resolution of the tail tables
L_MIN = 0.4           # the tables start at l = L_MIN


def up(x):
    """A float that is >= every point of the ball x."""
    u = x.upper()
    return math.nextafter(float(u.mid()) + float(u.rad()), math.inf)


def lo(x):
    """A float that is <= every point of the ball x."""
    u = x.lower()
    return math.nextafter(float(u.mid()) - float(u.rad()), -math.inf)


def cell(a, b):
    """The ball [a, b] (a <= b floats)."""
    return arb((a + b) / 2.0, (b - a) / 2.0 * (1 + 1e-12) + 1e-300)


class Gram:
    """Base class.  A subclass defines `raw(order, x)`: Gamma_raw^{(order)}(x) for order 0, 1, 2 on an arb ball or float x >= 0 (the raw lengths), the constant `lam0` = -Gamma_raw''(0) and `n`.
    The "sigma-equivalent" Gram is Gamma(x) = Gamma_raw(s x) with s = 1/sqrt(2 lam0): then -Gamma''(0) = 1/2 as for the Gaussian of width sigma."""
    name = '?'
    n = 1
    lam0 = 0.5
    cusp = False

    def __init__(self):
        self.s = 1.0 / math.sqrt(2.0 * self.lam0)
        self._tab = None

    # ---- values in sigma-equivalent units (arb ball or float)
    def val(self, order, x):
        return self.raw(order, x * self.s if not isinstance(x, arb) else x * arb(self.s)) * (self.s ** order)

    def G(self, x): return self.val(0, x)
    def G1(self, x): return self.val(1, x)
    def G2(self, x): return self.val(2, x)

    # ---- numpy evaluation for the simulator (floats, arrays)
    def fG(self, x):
        return self.rawf(0, np.asarray(x, float) * self.s)

    def fdG(self, x):
        return self.rawf(1, np.asarray(x, float) * self.s) * self.s

    def fd2G(self, x):
        return self.rawf(2, np.asarray(x, float) * self.s) * self.s ** 2

    # ---- cells
    def _cells(self, a, b, n):
        if n is None:
            n = int(min(max(math.ceil((b - a) / 0.0015), 24), 1200))
        xs = np.linspace(a, b, n + 1)
        return [(float(xs[i]), float(xs[i + 1])) for i in range(n)]

    # ---- tails
    def _build_tab(self):
        n_ = int(round((L_MAX - L_MIN) / H_TAB))
        ls = [L_MIN + i * H_TAB for i in range(n_ + 1)]
        t0, t1, t2, th = [0.0] * (n_ + 1), [0.0] * (n_ + 1), [0.0] * (n_ + 1), [0.0] * (n_ + 1)
        # tail beyond L_MAX: the value at L_MAX (non-increasing envelopes assumed there)
        xe = arb(L_MAX)
        c0, c1, c2 = up(abs(self.G(xe))), up(abs(self.G1(xe))), up(abs(self.G2(xe)))
        ch = max(c2, c1 / L_MAX)
        t0[n_], t1[n_], t2[n_], th[n_] = c0, c1, c2, ch
        for i in range(n_ - 1, -1, -1):
            a, b = ls[i], ls[i + 1]
            v0, v1, v2, vh = self.cell_up(a, b)
            t0[i], t1[i], t2[i], th[i] = max(v0, t0[i + 1]), max(v1, t1[i + 1]), max(v2, t2[i + 1]), max(vh, th[i + 1])
        self._tab = (ls, t0, t1, t2, th)

    def cell_up(self, a, b):
        """Upper bounds over the cell [a, b] (a > 0, sigma-equivalent units) of |Gamma|, |Gamma'|, |Gamma''| and max(|Gamma''|, |Gamma'|/d)."""
        ball = cell(a, b)
        v0 = up(abs(self.G(ball)))
        v1 = up(abs(self.G1(ball)))
        v2 = up(abs(self.G2(ball)))
        return v0, v1, v2, max(v2, v1 / a)

    def tail(self, order, l, hess_n=None):
        """Upper bound of sup_{d >= l} |Gamma^{(order)}(d)| (order 0, 1, 2) or, order = 'H', of sup_{d>=l} ||Hess Gamma(|x|)|| for dimension hess_n (n = 1: |Gamma''|)."""
        if self._tab is None:
            self._build_tab()
        ls, t0, t1, t2, th = self._tab
        l = float(l)
        if l <= L_MIN:
            i = 0
        else:
            i = int(math.floor((l - L_MIN) / H_TAB + 1e-12))
        i = min(i, len(ls) - 1)
        if order == 0:
            return t0[i]
        if order == 1:
            return t1[i]
        if order == 2:
            return t2[i]
        nn = self.n if hess_n is None else hess_n
        return t2[i] if nn == 1 else th[i]

    # ---- own constants (directed bounds) -- generic cell versions; Bessel kernels override the pieces singular at 0
    def ptilde(self, x):
        """psi(x)/x = -Gamma'(x)/x = s^2 ptilde_raw(s x) (sigma-equivalent units) on a ball possibly containing 0; ptilde_raw is regular at 0 (subclasses)."""
        sx = x * arb(self.s) if isinstance(x, arb) else x * self.s
        return self.ptilde_raw(sx) * (self.s ** 2)

    def ptilde_raw(self, x):
        return -self.raw(1, x) / x

    def g_min(self, x, N=None):
        """lower bound of min_{[0,x]} Gamma.  Gamma is decreasing on [0, x] whenever lampsi(x) > 0 (psi/d > 0 on (0, x]), which `box` has checked before it calls this: the minimum is Gamma(x)."""
        return lo(self.G(arb(x)))

    def psi_at(self, x):
        return lo(-self.G1(arb(x)))

    def lampsi(self, x, N=None):
        """lower bound of inf_{(0,x]} psi(d)/d"""
        return min(lo(self.ptilde_ball(a, b)) for a, b in self._cells(0.0, x, N))

    def ptilde_ball(self, a, b):
        return self.ptilde(cell(a, b))

    def psihat(self, x, N=None):
        """upper bound of sup_{[0,x]} psi"""
        return max(up(-self.G1(cell(a, b))) for a, b in self._cells(0.0, x, N))

    def rho_min(self, x, N=None):
        """lower bound of inf_{[0,x]} (-Gamma'')"""
        return min(lo(-self.G2(cell(a, b))) for a, b in self._cells(0.0, x, N))

    # global sups
    def M1(self):
        return self.tail(1, L_MIN) if self.cusp else max(self.tail(1, L_MIN), self.psihat(L_MIN))

    def M2(self, hess_n=None):
        nn = self.n if hess_n is None else hess_n
        glob = max(abs(self.lam0 * self.s ** 2), self.tail(2, L_MIN))
        if nn > 1:
            glob = max(glob, self.tail('H', L_MIN, nn))
        # on [0, L_MIN]: |Gamma''| <= its sup over cells; |Gamma'|/d <= lampsi-type sup (ptilde is the Hessian's tangential eigenvalue)
        sup_small = max(up(abs(self.G2(cell(a, b)))) for a, b in self._cells(0.0, L_MIN, 16))
        sup_pt = max(up(abs(self.ptilde_ball(a, b))) for a, b in self._cells(0.0, L_MIN, 16)) if nn > 1 else 0.0
        return max(glob, sup_small, sup_pt)


# ------------------------------------------------------------------------------------------------------------------------------------------
class Gaussian(Gram):
    name = 'gauss'

    def __init__(self, n=1):
        self.n = n
        self.lam0 = 0.5
        super().__init__()

    def raw(self, order, x):
        if isinstance(x, arb):
            e = (-(x * x) / 4).exp()
            return (e, -x / 2 * e, (x * x / 4 - arb(1) / 2) * e)[order]
        e = math.exp(-x * x / 4)
        return (e, -x / 2 * e, (x * x / 4 - 0.5) * e)[order]

    def rawf(self, order, x):
        e = np.exp(-x * x / 4)
        return (e, -x / 2 * e, (x * x / 4 - 0.5) * e)[order]

    def ptilde_raw(self, x):
        return self.raw(0, x) / 2


class MaternHalf(Gram):
    """Gram = Matern nu = p + 1/2, spectral parameter 1: Gamma_raw(d) = P(d) exp(-d), P(0) = 1 (the Gram of the exponential basis function in R^1 (p = 1), R^3 (p = 2); of Matern basis ... see FAMILIES)."""

    def __init__(self, p, n=1):
        self.p = p
        self.nu = p + 0.5
        self.n = n
        self.name = 'matern%g' % self.nu
        self.P = [Fraction(0)] * (p + 1)
        f = math.factorial
        for i in range(p + 1):
            self.P[p - i] = Fraction(f(p), f(2 * p)) * Fraction(f(p + i), f(i) * f(p - i)) * (2 ** (p - i))
        assert self.P[0] == 1
        # derivatives: Gamma' = (P' - P) e^{-d}, Gamma'' = (P'' - 2P' + P) e^{-d}
        self.P1 = self._sub(self._der(self.P), self.P)
        self.P2 = self._sub(self._sub(self._der(self._der(self.P)), self._der(self.P)), self._sub(self._der(self.P), self.P))
        # -Gamma'(d)/d = Q(d) e^{-d}, Q = (P - P')/d exact
        num = self._sub(self.P, self._der(self.P))
        assert num[0] == 0
        self.Q = num[1:]
        self.lam0 = float(-(self.P2[0]))
        assert abs(self.lam0 - 1.0 / (2 * (self.nu - 1))) < 1e-12 and self.P1[0] == 0
        super().__init__()

    @staticmethod
    def _der(P):
        return [i * P[i] for i in range(1, len(P))] or [Fraction(0)]

    @staticmethod
    def _sub(A, B):
        m = max(len(A), len(B))
        return [(A[i] if i < len(A) else 0) - (B[i] if i < len(B) else 0) for i in range(m)]

    @staticmethod
    def _horner(C, x):
        if isinstance(x, arb):
            acc = arb(0)
            for c in reversed(C):
                acc = acc * x + arb(fmpq(c.numerator, c.denominator))
            return acc
        acc = 0.0
        for c in reversed(C):
            acc = acc * x + float(c)
        return acc

    def raw(self, order, x):
        P = (self.P, self.P1, self.P2)[order]
        e = (-x).exp() if isinstance(x, arb) else math.exp(-x)
        return self._horner(P, x) * e

    def rawf(self, order, x):
        P = (self.P, self.P1, self.P2)[order]
        return np.polyval([float(c) for c in reversed(P)], x) * np.exp(-x)

    def ptilde_raw(self, x):
        e = (-x).exp() if isinstance(x, arb) else math.exp(-x)
        return self._horner(self.Q, x) * e


class MaternBessel(Gram):
    """Gram = Matern nu, nu integer >= 2: Gamma_raw = c_nu d^nu K_nu(d), c_nu = 2^{1-nu}/Gamma(nu) (the Gram of a Matern/exponential basis function in the plane); with
    G_mu(d) := c_mu d^mu K_mu(d) (mu >= 1; decreasing from 1; d/dd G_mu = -d G_{mu-1}/(2(mu-1)) for mu >= 2):
        Gamma = G_nu,   -Gamma' = d G_{nu-1}/(2(nu-1)),   -Gamma'' = rho = G_{nu-1}/(2(nu-1)) - d^2 G_{nu-2}/(4(nu-1)(nu-2))  (nu > 2),   rho = G_1/2 - d^2 K_0(d)/2  (nu = 2).
    psi/d = G_{nu-1}/(2(nu-1)) is a Matern Gram: positive, decreasing.  arb's Bessel K on a ball is far too wide, so every directed bound below uses POINT evaluations of the decreasing
    functions G_mu and K_0 at cell endpoints (monotone brackets), never a ball."""

    def __init__(self, nu, n=2):
        assert nu >= 2
        self.nu = nu if int(nu) != nu else int(nu)
        self.n = n
        self.name = 'matern%g' % self.nu
        self.lam0 = 1.0 / (2 * (self.nu - 1))
        super().__init__()

    def c(self, mu):
        return 2.0 ** (1 - mu) / math.gamma(mu)

    @staticmethod
    def K(mu, x):
        """K_mu(x) for an arb point x > 0 (integer mu >= 0).  python-flint's K_0 is ill-conditioned for large x (a ball of radius 1e-3 at x = 26): K_0 = K_2 - 2 K_1/x there."""
        if mu == 0 and float(x.mid()) > 6.0:                            # (integer orders only: for half-integer orders the plain call is accurate)
            v = x.bessel_k(2) - 2 * x.bessel_k(1) / x
        else:
            v = x.bessel_k(mu)
        if not (float(v.rad()) <= 1e-7 * abs(float(v.mid()))):
            raise ArithmeticError('Bessel K ball too wide at %s: %s' % (x, v))
        return v

    def Gmu(self, mu, d):
        """G_mu(d) at a point d (float >= 0) as an arb (exact point evaluation)."""
        if d == 0:
            return arb(1)
        x = arb(d)
        return arb(self.c(mu)) * x ** mu * self.K(mu, x)

    def raw(self, order, x):
        nu = self.nu
        if isinstance(x, arb):
            assert x.rad() == 0 or float(x.rad()) < 1e-12, 'Bessel kernels: point evaluations only'
            if order == 0:
                return arb(self.c(nu)) * x ** nu * self.K(nu, x)
            if order == 1:
                return -arb(self.c(nu)) * x ** nu * self.K(nu - 1, x)
            return arb(self.c(nu)) * (x ** nu * self.K(nu - 2, x) - x ** (nu - 1) * self.K(nu - 1, x))
        return float(self.rawf(order, np.asarray(x, float)))

    def rawf(self, order, x):
        import scipy.special as sp
        nu = self.nu
        x = np.asarray(x, float)
        with np.errstate(all='ignore'):
            xs = np.where(x > 0, x, 1.0)
            if order == 0:
                return np.where(x > 0, self.c(nu) * xs ** nu * sp.kv(nu, xs), 1.0)
            if order == 1:
                return np.where(x > 0, -self.c(nu) * xs ** nu * sp.kv(nu - 1, xs), 0.0)
            v = self.c(nu) * (xs ** nu * sp.kv(nu - 2, xs) - xs ** (nu - 1) * sp.kv(nu - 1, xs))
            return np.where(x > 0, v, -self.lam0)

    # ---- raw-unit brackets on a cell [A, B] (0 <= A < B): u0 >= Gamma, u1 >= |Gamma'|, and P in [Plo, Pup], N in [Nlo, Nup] with rho = P - N (P = G_{nu-1}/(2(nu-1)) decreasing, N = d^2 G_{nu-2}/(4(nu-1)(nu-2)) (or d^2 K_0/2))
    def _cell_raw(self, A, B):
        nu = self.nu
        u0 = up(self.Gmu(nu, A))
        Pup = up(self.Gmu(nu - 1, A)) / (2 * (nu - 1))
        Plo = lo(self.Gmu(nu - 1, B)) / (2 * (nu - 1))
        u1 = B * Pup
        if nu > 2:
            Nup = B * B * up(self.Gmu(nu - 2, A)) / (4 * (nu - 1) * (nu - 2))          # G_{nu-2} decreasing, d^2 <= B^2
            Nlo = A * A * lo(self.Gmu(nu - 2, B)) / (4 * (nu - 1) * (nu - 2))
        else:
            if A == 0:
                assert B <= 1.4                                         # d^2 K_0(d) increases on (0, 1.4): <= B^2 K_0(B)
                Nup = B * B * up(self.K(0, arb(B))) / 2.0
                Nlo = 0.0
            else:
                Nup = B * B * up(self.K(0, arb(A))) / 2.0                # K_0 decreasing
                Nlo = A * A * lo(self.K(0, arb(B))) / 2.0
        return u0, u1, Pup, Nup, Plo, Nlo

    def cell_up(self, a, b):
        s = self.s
        u0, u1, Pup, Nup, Plo, Nlo = self._cell_raw(a * s, b * s)
        u2 = max(abs(Plo - Nup), abs(Pup - Nlo))
        return u0, u1 * s, u2 * s * s, max(u2, Pup) * s * s

    # ---- own constants (eq units)
    def g_min(self, x, N=None):
        return lo(self.Gmu(self.nu, x * self.s))

    def lampsi(self, x, N=None):
        return lo(self.Gmu(self.nu - 1, x * self.s)) / (2 * (self.nu - 1)) * self.s ** 2 * (1 - 1e-15)

    def psi_at(self, x):
        return lo(self.Gmu(self.nu - 1, x * self.s)) / (2 * (self.nu - 1)) * x * self.s * self.s * (1 - 1e-15)

    def psihat(self, x, N=None):
        s = self.s
        N = 64 if N is None else N
        xs = np.linspace(0.0, x, N + 1)
        best = max(self._cell_raw(float(xs[i]) * s, float(xs[i + 1]) * s)[1] for i in range(N))
        return best * s * (1 + 1e-15)

    def rho_min(self, x, N=None):
        nu, s = self.nu, self.s
        N = 64 if N is None else N
        xs = np.linspace(0.0, x, N + 1)
        worst = math.inf
        for i in range(N):
            A, B = float(xs[i]) * s, float(xs[i + 1]) * s
            _, _, _, Nup, Plo, _ = self._cell_raw(A, B)
            worst = min(worst, Plo - Nup)
        return worst * s * s

    def M1(self):
        return max(self.tail(1, L_MIN), self.psihat(L_MIN))

    def M2(self, hess_n=None):
        nn = self.n if hess_n is None else hess_n
        v = max(self.lam0 * self.s ** 2, self.tail(2, L_MIN))                 # on [0, L_MIN]: |Gamma''| <= max(P, N) <= max(lam0, N(L_MIN)) and the tangential eigenvalue psi/d <= lam0
        if nn > 1:
            v = max(v, self.tail('H', L_MIN, nn))
        return v * (1 + 1e-12)


class RationalQuadratic(Gram):
    """Gram Gamma(x) = q^{-beta}, q = 1 + x^2/(2 beta): a scale mixture of Gaussians (positive definite in every dimension) with ALGEBRAIC decay x^{-2 beta};
    Gamma' = -x q^{-beta-1}, Gamma'' = -q^{-beta-1} + ((beta+1)/beta) x^2 q^{-beta-2}, -Gamma''(0) = 1 (raw units)."""

    def __init__(self, beta, n=1):
        self.beta = float(beta)
        self.n = n
        self.name = 'ratquad%g' % beta
        self.lam0 = 1.0
        super().__init__()

    def raw(self, order, x):
        b = self.beta
        q = (arb(1) + x * x / (2 * b)) if isinstance(x, arb) else 1.0 + x * x / (2 * b)
        return (q ** (-b), -x * q ** (-b - 1), -q ** (-b - 1) + (b + 1) / b * x * x * q ** (-b - 2))[order]

    def rawf(self, order, x):
        return self.raw(order, np.asarray(x, float))

    def ptilde_raw(self, x):
        b = self.beta
        q = (arb(1) + x * x / (2 * b)) if isinstance(x, arb) else 1.0 + x * x / (2 * b)
        return q ** (-b - 1)


class Exp1(Gram):
    """Gamma(d) = exp(-d): a CUSP at 0 (Gamma'(0+) = -1, the Gram of phi with Fourier transform (1 + w^2)^{-1/2} in R^1, phi not in H^1): outside the smooth class of Proposition G;
    handled by `box_cusp` (Proposition G'')."""
    name = 'exp1'
    cusp = True

    def __init__(self, n=1):
        self.n = n
        self.lam0 = 0.5
        super().__init__()

    def raw(self, order, x):
        e = (-x).exp() if isinstance(x, arb) else math.exp(-x)
        return (e, -e, e)[order]

    def rawf(self, order, x):
        e = np.exp(-np.asarray(x, float))
        return (e, -e, e)[order]


# ------------------------------------------------------------------------------------------------------------------------------------------
def matern_gram(nu, n=1):
    """Gram = Matern-nu (nu half-integer >= 3/2 or integer >= 2)."""
    if abs(nu * 2 - round(nu * 2)) < 1e-12 and int(round(nu * 2)) % 2 == 1:
        return MaternHalf(int(round(nu - 0.5)), n)
    return MaternBessel(int(round(nu)), n)


def family(kind, n):
    """The Gram of a kernel family in R^n: 'gauss' (any n), 'exp' (phi = exp(-|x|): Gram Matern nu = 1 + n/2), 'mat32' (phi = Matern-3/2: nu_G = 3 + n/2), 'mat52' (nu_G = 5 + n/2).
    (A Matern-nu basis function in R^n has the Fourier transform ~ (1 + w^2)^{-(nu + n/2)}; the Gram's is the square, hence nu_G = 2 nu + n/2.)"""
    if kind == 'gauss':
        return Gaussian(n)
    nub = {'exp': 0.5, 'mat32': 1.5, 'mat52': 2.5}[kind]
    return matern_gram(2 * nub + n / 2.0, n)


# ------------------------------------------------------------------------------------------------------------------------------------------
def _sqrt(x):
    return x.nonnegative_part().sqrt()               # the radicands are sums of squares: a ball that dips below 0 by rounding is clipped


def box(kern, K, R, Delta, r, Omega=0.0, hc='norm2', need_N=True, sums=None, hess_n=None, cond_c=1.0):
    """Proposition G: the conditions (L), (C), (N) for K features of amplitudes in [1, R], separations >= Delta, box radius r (sigma-equivalent units), dormant budget Omega >= 0 (units
    of a_min).  `sums` = (S0, S1, S2) the coherence sums over the other features (default (K-1) m_j(Delta - 2r)); hc = 'norm2' (blockwise spectral norms, any n) or 'rowsum' (the
    maximum absolute row sum of the entries, the Gaussian note's H_c, n = 1 only).  cond_c >= 1: the condition number of the centre mobility matrix (anisotropic Gram with a common covariance:
    (C) becomes psi(r) > cond_c E_c).  Returns dict(ok, margins ...); `ok` = certainly true in ball arithmetic."""
    n = kern.n if hess_n is None else hess_n
    out = dict(ok=False)
    ell = Delta - 2 * r
    if ell <= L_MIN or r <= 0:
        return out
    if sums is None:
        S0 = arb(up(arb((K - 1) * kern.tail(0, ell))))
        S1 = arb(up(arb((K - 1) * kern.tail(1, ell))))
        S2 = arb(up(arb((K - 1) * kern.tail(2 if n == 1 else 'H', ell, n))))
        if K == 1:
            S0 = S1 = S2 = arb(0)
    else:
        S0, S1, S2 = [arb(float(v)) for v in sums]
    Rr, Om = arb(R), arb(Omega)
    one = arb(1)
    if not (S0 < 1):
        return out
    M1 = arb(kern.M1()) if Omega > 0 else arb(0)
    M2 = arb(kern.M2(n)) if Omega > 0 else arb(0)
    vbar = (Rr * (one + S0) + Om) / (one - S0)
    Ev = S0 * (Rr + vbar) + Om
    Ec = S1 * (Rr + vbar) + M1 * Om
    g_r = arb(kern.g_min(r))
    psi_r = arb(kern.psi_at(r))
    lps_r = arb(kern.lampsi(r))
    vlo = g_r - Ev
    Cm = psi_r - arb(cond_c) * Ec
    out.update(Ev=float(Ev.mid()), Ec=float(Ec.mid()), vbar=float(vbar.mid()), vlo=float(vlo.mid()), Cmargin=float(Cm.mid()), kappa=float(S0.mid()))
    if not (vlo > 0 and lps_r > 0 and Cm > 0):
        return out
    if not need_N:
        out['ok'] = True
        return out
    dt = Ec / lps_r                                         # |d| <= d-tilde at every critical point of the box
    out['dtilde'] = float(dt.mid())
    if not (dt < r):
        return out
    dtu = up(dt)
    gam_d = arb(kern.g_min(dtu))
    vlo_d = gam_d - Ev
    if not (vlo_d > 0):
        return out
    s_max = Rr * arb(kern.psihat(dtu))
    rho_d = arb(kern.rho_min(dtu))
    t_rad = vlo_d * rho_d
    lam2 = (one + t_rad - _sqrt((one - t_rad) * (one - t_rad) + 4 * s_max * s_max)) / 2
    if n > 1:
        t_perp = vlo_d * arb(kern.lampsi(dtu))
        lam = lam2 if lam2 < t_perp else t_perp                # min of the radial 2x2 block's least eigenvalue and the tangential stiffness
        if not (lam2 < t_perp) and not (t_perp < lam2):
            lam = arb(min(lo(lam2), lo(t_perp)))
    else:
        lam = lam2
    # the cross Hessian
    Rv = Rr + vbar
    if hc == 'rowsum':
        rowv = S0 + Rv * S1 + vbar * S1 + M1 * Om
        rowd = Rv * S1 + vbar * Rv * S2 + vbar * S1 + vbar * vbar * S2 + M1 * Om + vbar * M2 * Om
        Hc = rowv if float(rowv.mid()) > float(rowd.mid()) else rowd
    else:
        # D: [[0, Rv S1 + M1 Om], [., vbar Rv S2 + vbar M2 Om]],  B (sum over the other features): [[S0, vbar S1], [vbar S1, vbar^2 S2]]
        def n2(a, b, c):
            return (a + c + _sqrt((a - c) * (a - c) + 4 * b * b)) / 2
        Hc = n2(arb(0), Rv * S1 + M1 * Om, vbar * Rv * S2 + vbar * M2 * Om) + n2(S0, vbar * S1, vbar * vbar * S2)
    Nm = lam - Hc
    out.update(lam=float(lam.mid()), Hc=float(Hc.mid()), Nmargin=float(Nm.mid()), ok=bool(Nm > 0))
    return out


R_GRID_N = 36


def r_grid(d0, rmax=1.35, n=R_GRID_N):
    return [float(x) for x in np.linspace(max(d0, 1e-3) * 1.000001 + 1e-9, rmax, n)]


def thr_for_radius(kern, K, R, r, step=0.02, lo_=1.0, hi=30.0, **kw):
    """Smallest Delta on the grid lo_ + i step with box ok at radius r (bisection: the conditions are monotone in Delta); nan if none up to hi."""
    n_ = int(round((hi - lo_) / step))
    D = lambda i: round(lo_ + i * step, 10)
    if not box(kern, K, R, D(n_), r, **kw)['ok']:
        return float('nan')
    a, b = 0, n_
    if box(kern, K, R, D(0), r, **kw)['ok']:
        return D(0)
    while b - a > 1:
        m = (a + b) // 2
        if box(kern, K, R, D(m), r, **kw)['ok']:
            b = m
        else:
            a = m
    return D(b)


def entry_threshold(kern, K, R, d0, step=0.02, lo_=1.0, hi=30.0, rmax=1.35, n_r=R_GRID_N, return_r=False, **kw):
    """Delta_entry(K, R, d0; Gamma): the smallest Delta (grid `step`) at which a box of radius r > d0 exists."""
    best, br = float('nan'), None
    for r in r_grid(d0, rmax, n_r):
        t = thr_for_radius(kern, K, R, r, step, lo_, hi, **kw)
        if t == t and (best != best or t < best):
            best, br = t, r
    return (best, br) if return_r else best


def clean_index(kern, K, R, Delta, r, hess_n=None):
    """The coherence index  c(Delta; r) = (K-1) R^2 max(m0, m1, m2)(Delta - 2r)  of the clean corollary (NOTE.md, Corollary G1)."""
    n = kern.n if hess_n is None else hess_n
    ell = Delta - 2 * r
    m = max(kern.tail(0, ell), kern.tail(1, ell), kern.tail(2 if n == 1 else 'H', ell, n))
    return (K - 1) * R * R * m


# ------------------------------------------------------------------------------------------------------------------------------------------
# the cusp case (Gamma'(0+) < 0): Filippov / sliding-mode variant of the box conditions
# ------------------------------------------------------------------------------------------------------------------------------------------
def box_cusp(kern, K, R, Delta, r, sums=None):
    """Conditions of Proposition G'' for a Gram with a cusp, Gamma(d) = exp(-d) (kern = Exp1; lengths in correlation lengths): (L) Gamma(r) > E_v, (C') psi_min(r) = inf_{(0,r]} (-Gamma') = exp(-r) > E_c,
    kappa = (K-1) m0(Delta - 2r) < 1.  No Hessian condition: every critical point of the box has all offsets exactly 0 (the centres lock onto the features)."""
    out = dict(ok=False)
    ell = Delta - 2 * r
    if ell <= L_MIN or r <= 0:
        return out
    S0 = arb(up(arb((K - 1) * kern.tail(0, ell))))
    S1 = arb(up(arb((K - 1) * kern.tail(1, ell))))
    Rr, one = arb(R), arb(1)
    if not (S0 < 1):
        return out
    vbar = (Rr * (one + S0)) / (one - S0)
    Ev = S0 * (Rr + vbar)
    Ec = S1 * (Rr + vbar)
    g_r = arb(kern.g_min(r))
    psi_min = arb(lo((-kern.G1(arb(r)))))
    vlo, Cm = g_r - Ev, psi_min - Ec
    out.update(Ev=float(Ev.mid()), Ec=float(Ec.mid()), vlo=float(vlo.mid()), Cmargin=float(Cm.mid()), kappa=float(S0.mid()), ok=bool(vlo > 0 and Cm > 0))
    return out


def entry_threshold_cusp(kern, K, R, d0, step=0.02, lo_=1.0, hi=60.0, rmax=3.0, n_r=40):
    best = float('nan')
    for r in r_grid(d0, rmax, n_r):
        n_ = int(round((hi - lo_) / step))
        D = lambda i: round(lo_ + i * step, 10)
        if not box_cusp(kern, K, R, D(n_), r)['ok']:
            continue
        a, b = 0, n_
        if box_cusp(kern, K, R, D(0), r)['ok']:
            t = D(0)
        else:
            while b - a > 1:
                m = (a + b) // 2
                if box_cusp(kern, K, R, D(m), r)['ok']:
                    b = m
                else:
                    a = m
            t = D(b)
        if best != best or t < best:
            best = t
    return best


# ------------------------------------------------------------------------------------------------------------------------------------------
# the clean (coherence) form, Corollary G1
# ------------------------------------------------------------------------------------------------------------------------------------------
def clean_ok(kern, R, r, c1, hess_n=None):
    """The inequalities of Corollary G1 at the level c1 (floats with directed bounds): with c := (K-1) R^2 max(m0,m1,m2)(Delta - 2r) <= c1 <= 1/4 one has v_bar <= 1.7 R, E_v, E_c <= 2.7 c1, d-tilde <= d1 = 2.7 c1 / lam_r,
       ||H_cross|| <= 12 c1; need  2.7 c1 < Gamma(r),  2.7 c1 < psi(r),  d1 < r,  and  Lambda(d1) > 12 c1,  Lambda(d1) = the least eigenvalue floor of the one-unit block on |d| <= d1, v >= Gamma(d1) - 2.7 c1."""
    n = kern.n if hess_n is None else hess_n
    ev = 2.7 * c1
    g_r, psi_r, lam_r = kern.g_min(r), kern.psi_at(r), kern.lampsi(r)
    if not (c1 <= 0.25 and ev < g_r and ev < psi_r and lam_r > 0):
        return False
    d1 = ev / lam_r
    if not d1 < r:
        return False
    vlo = kern.g_min(d1) - ev
    if vlo <= 0:
        return False
    rho = kern.rho_min(d1)
    t = vlo * rho
    s_ = R * kern.psihat(d1)
    lam2 = 0.5 * (1 + t - math.sqrt((1 - t) ** 2 + 4 * s_ * s_))
    Lam = lam2 if n == 1 else min(lam2, vlo * kern.lampsi(d1))
    return Lam > 12.0 * c1


def clean_cstar(kern, R, r, hess_n=None):
    """c_*(r): the largest c1 (bisection) with clean_ok (constants 2.7 and 12 of Corollary G1).  If (K-1) R^2 max(m0, m1, m2)(Delta - 2r) <= c_*(r) for some r > d0 then (L), (C), (N) hold (Corollary G1)."""
    key = (float(R), float(r), hess_n)
    cache = kern.__dict__.setdefault('_cs_cache', {})
    if key in cache:
        return cache[key]
    if not clean_ok(kern, R, r, 1e-9, hess_n):
        cache[key] = 0.0
        return 0.0
    lo_, hi_ = 1e-9, 0.25
    if clean_ok(kern, R, r, hi_, hess_n):
        cache[key] = hi_
        return hi_
    for _ in range(30):
        mid = 0.5 * (lo_ + hi_)
        lo_, hi_ = (mid, hi_) if clean_ok(kern, R, r, mid, hess_n) else (lo_, mid)
    cache[key] = lo_
    return lo_


def clean_threshold(kern, K, R, d0, step=0.02, lo_=1.0, hi=80.0, rmax=1.35, n_r=R_GRID_N, hess_n=None):
    """Smallest Delta (grid step) such that for some r > d0:  (K-1) R^2 max(m0, m1, m2)(Delta - 2r) <= c_*(r)."""
    n = kern.n if hess_n is None else hess_n
    best = float('nan')
    for r in r_grid(d0, rmax, n_r):
        cs_ = clean_cstar(kern, R, r, n)
        if cs_ <= 0:
            continue
        # smallest ell on the table grid with (K-1) R^2 m(ell) <= cs_
        lo_l, hi_l = L_MIN, L_MAX
        f = lambda l: (K - 1) * R * R * max(kern.tail(0, l), kern.tail(1, l), kern.tail(2 if n == 1 else 'H', l, n)) <= cs_
        if not f(hi_l):
            continue
        a, b = 0, int(round((hi_l - lo_l) / H_TAB))
        while b - a > 1:
            m = (a + b) // 2
            if f(lo_l + m * H_TAB):
                b = m
            else:
                a = m
        ell = lo_l + b * H_TAB
        t = math.ceil((ell + 2 * r) / step - 1e-9) * step
        if best != best or t < best:
            best = t
    return best


def decay_modulus(kern, q, a, l0=3.0, hess_n=None):
    """C(a, q) = max over the table grid ell in [l0, L_MAX] of  m(ell) exp(a ell^q)  with m = max(m0, m1, m2) (monotone envelope: for ell in a table cell the bound m(ell_i) exp(a ell_{i+1}^q) is used):
    then m(ell) <= C exp(-a ell^q) for ell in [l0, L_MAX]."""
    n = kern.n if hess_n is None else hess_n
    kern.tail(0, l0)
    ls = kern._tab[0]
    best = 0.0
    for i in range(len(ls) - 1):
        if ls[i] < l0:
            continue
        m = max(kern._tab[1][i], kern._tab[2][i], kern._tab[3][i] if n == 1 else kern._tab[4][i])
        best = max(best, m * math.exp(a * ls[i + 1] ** q))
    return best
