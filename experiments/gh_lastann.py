"""Part VA of D‴.61 (`exp_gh_peano.py`), adopted from referee D3do (round 120; written by that referee, read in full),
verified independently by referee D3dp (round 121), who found one soundness defect — a NaN handoff interval dropping the
outer connector of the η-box at arg η ≈ 0 next to |η| = 2 — and its fix, applied here with that round's other changes
(`handoff_range`, `outer_box`, `handoff_check`, `soundness`, the asserted B0 check, exact `hi − lo`, `seg_bound`'s |σ| = 0
guard): a ball-arithmetic certificate of D‴.61 (iii)(b)'s valley condition on the last annulus
1.995 <= |η| < 2, 0 <= arg η <= π/2 (with η ↦ η̄, −η the whole annulus), along part VC's path family.

The bound (`cert_E`, exact): with σ = (η/2)(1 + z), E(z; η) = F_η(σ) − F_η(η/2) = log|1 − z²| + Re(c z²) + C(σ), c = η²/4,
C = 0 for |σ| ≤ 1 and h(|σ|) − (|Re σ| − 1)₊² beyond, h(r) = (r² − 1)/2 − log r. Write η = (2 − ε)e^{iθ}, ω = 1 − ε/2,
A = c − 1 (Re A < 0), κ = Im A/|Re A|, b = (1 + t)η/2 with t = ε/4 (z = t: part V's chord point, inside the unit disc).

The path for each η (z = X + iY):
  J  junction  X = t, Y = λκt, λ ∈ [0, 1]                        (from b straight up to the floor)
  F  floor     Y = κX, from X = t until κX = θ                     (the exact quadratic floor)
  C  cap       Y = θ, up to X = δ₀                                 (Im z = θ)
  O  outer     from δ₀ + i·min(κδ₀, θ): a vertical connector, then a descent polyline (per η-box) to |1 − z| ≤ 0.02
  E  end       straight to z = 1 (part V's end piece, at R2 = 2)
(J, F, C is part VC's family Im z = min(κ Re z, θ) started at b; the paper's 0-side and chord are part V's, exact.)

Two scaled forms, both exact identities, cover J, F and C for X ≤ δ₀ uniformly as ε → 0 and θ → 0:

  Form A, (u, q, ρ) = (ε/X, θ²/X, X), q ≤ QA:   E/X² = fA(u, q, ρ),
  Form B, (u, p, ϑ) = (ε/X, X/θ², θ²), p ≤ PB:   E/X² = fB(u, p, ϑ),

with every division by ε, θ or X done analytically (sinc-type series in θ², log(1 + v) − v and h(1 + δ) − δ² by series
with rigorous tails, |σ| − 1 = (ω²|1+z|² − 1)/(ω|1+z| + 1)). Both are C¹ in the expansion variable (ρ, resp. ϑ) on each
branch, and at ρ = 0 (ϑ = 0) they reduce to the weight-4 limit (Ψ̃ of part VC in (u, q) coordinates), given in closed form
(`fA0`, `fB0`). The certificate on a box is the mean value form

  f(ρ) ≤ f(0) + ρ·sup f′([0, ρ_b])      (sup f′ by forward differentiation in ball arithmetic)

which is negative when f(0) ≤ 0 and sup f′ < 0 (ρ > 0 on the path), or when f(0) + ρ_b·max(sup f′, 0) < 0. f(0) vanishes
exactly where the limits degenerate — at the junction's λ = 0 (b on the real axis), on q = 0 (charts φ and (ε/θ, θ/√ε)) and
at the origin u = q = 0 (the matching region, E ≈ −εθ² − θ⁴/2 − X³/3) — and there sup f′ < 0 carries the certificate.

Everything is in Arb (python-flint) at PREC bits; comparisons are rigorous (Arb's > on balls); floats are only used to
choose paths (the outer descent polyline) and boxes, and every float taken out of a ball that shapes the certificate (the
connector's ends) is checked finite (round 121: a NaN is not an empty interval). Deterministic: fixed grids, depth-first stacks, no randomness except
the seeded validation.

Usage: python experiments/gh_lastann.py validate | A-J | A-FC | B-J | B-F | outer | end | all (or a comma-separated list)
       [--eps0 0.005] [--delta0 0.05] [--json out.json]
(round 131, D3dz-4: `all` runs to the end — the outer piece's log no longer carries the complex joints — and an unknown name is an error.
Part VA of `exp_gh_peano.py` runs the same pieces at ε₀ = 0.005 and, as part VA011, at ε₀ = 0.11.)
"""
import json
import math
import sys
import time

from flint import arb, acb, ctx

PREC = 128
ctx.prec = PREC
EPS0 = 0.005          # ε ≤ EPS0 ⇔ |η| ≥ 1.995 (part V's last band ends there)
DELTA0 = 0.05         # the scaled forms cover X ≤ DELTA0; the outer piece starts at Re z = DELTA0
QA = 2.5              # form A for q = θ²/X ≤ QA
PB = 0.5              # form B for p = X/θ² ≤ PB (QA·PB > 1: the forms overlap)
XB = arb(1) / 4       # b at t = XB·ε, i.e. u = 1/XB = 4 at the junction
UJ = 4.0              # = 1/XB (exact)
TH_MAX2 = arb.pi() ** 2 / 4   # ϑ = θ² ≤ (π/2)²
END_R = 0.02          # the end piece's radius (part V's)


# ---------------------------------------------------------------------------------------------------------------- balls
def hull(a, b):
    return arb(a).union(arb(b))


def up(x):
    """Rigorous float upper bound of an arb."""
    return math.nextafter(float(x.upper()), math.inf)


def pos(x):
    """x₊ = max(x, 0) on a ball."""
    if x >= 0:
        return x
    if x <= 0:
        return arb(0)
    return hull(0, x.upper())


class I:
    """Endpoint interval arithmetic with arb endpoints (each operation rounded outward by Arb's own lower()/upper()):
    used for the closed-form limits, where midpoint-radius balls would lose the sign of products of intervals at 0."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = arb(lo)
        self.hi = arb(lo if hi is None else hi)

    @staticmethod
    def _c(o):
        return o if isinstance(o, I) else I(o)

    def __add__(s, o):
        o = I._c(o)
        return I((s.lo + o.lo).lower(), (s.hi + o.hi).upper())
    __radd__ = __add__

    def __neg__(s):
        return I(-s.hi, -s.lo)

    def __sub__(s, o):
        return s + (-I._c(o))

    def __rsub__(s, o):
        return I._c(o) + (-s)

    def __mul__(s, o):
        o = I._c(o)
        ps = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        lo, hi = ps[0].lower(), ps[0].upper()
        for p in ps[1:]:
            lo = p.lower() if p.lower() < lo else lo
            hi = p.upper() if p.upper() > hi else hi
        return I(lo, hi)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = I._c(o)
        if not (o.lo > 0 or o.hi < 0):
            raise ZeroDivisionError
        return s * I((1 / o.hi).lower(), (1 / o.lo).upper())

    def __rtruediv__(s, o):
        return I._c(o) / s

    def sq(s):
        if s.lo >= 0:
            return I((s.lo * s.lo).lower(), (s.hi * s.hi).upper())
        if s.hi <= 0:
            return I((s.hi * s.hi).lower(), (s.lo * s.lo).upper())
        m = s.hi if s.hi > -s.lo else -s.lo
        return I(0, (m * m).upper())

    def __pow__(s, n):
        r = I(1)
        for _ in range(n):
            r = r * s
        if n % 2 == 0:
            r = I(r.lo if r.lo > 0 else arb(0), r.hi)
        return r

    def pos(s):
        return I(s.lo if s.lo > 0 else arb(0), s.hi if s.hi > 0 else arb(0))


class D:
    """A dual number over arb: value v and derivative d with respect to the expansion variable (forward mode)."""
    __slots__ = ('v', 'd')

    def __init__(self, v, d=0):
        self.v = v if isinstance(v, arb) else arb(v)
        self.d = d if isinstance(d, arb) else arb(d)

    def __add__(s, o):
        return D(s.v + o.v, s.d + o.d) if isinstance(o, D) else D(s.v + o, s.d)
    __radd__ = __add__

    def __sub__(s, o):
        return D(s.v - o.v, s.d - o.d) if isinstance(o, D) else D(s.v - o, s.d)

    def __rsub__(s, o):
        return D(o - s.v, -s.d)

    def __neg__(s):
        return D(-s.v, -s.d)

    def __mul__(s, o):
        return D(s.v * o.v, s.d * o.v + s.v * o.d) if isinstance(o, D) else D(s.v * o, s.d * o)
    __rmul__ = __mul__

    def __truediv__(s, o):
        if isinstance(o, D):
            if not (o.v > 0 or o.v < 0):
                raise ZeroDivisionError
            return D(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))
        return D(s.v / o, s.d / o)

    def __rtruediv__(s, o):
        if not (s.v > 0 or s.v < 0):
            raise ZeroDivisionError
        return D(o / s.v, -o * s.d / (s.v * s.v))

    def sqrt(s):
        if not (s.v > 0):
            raise ValueError('sqrt of a non-positive ball')
        r = s.v.sqrt()
        return D(r, s.d / (2 * r))

    def pos(s):
        """x₊ with derivative 1{x > 0}x′ (the hull of 0 and x′ where x straddles 0; used only inside C¹ compositions)."""
        if s.v >= 0:
            return s
        if s.v <= 0:
            return D(0, 0)
        return D(hull(0, s.v.upper()), arb(0).union(s.d))


def val(x):
    return x.v if isinstance(x, D) else x


def mag(x):
    """A rigorous upper bound of |x| (float)."""
    return up(abs(val(x)))


# -------------------------------------------------------------------------------------- power series with rigorous tails
TAIL_FLOOR = 1e-300   # round 131 (D3dz-1): the least tail radius (the tails are evaluated in doubles)


def _series(x, coef, K, tail0, tail1):
    """Σ_k coef(k) x^k on a ball or dual x, truncated at K terms with the tail bounds tail0(r), tail1(r) (r ≥ |x|) of
    the value's and the derivative's remainder. Horner in Arb. The tails are evaluated in IEEE doubles with the factor
    1 + 10⁻⁷ (a handful of roundings, each below 2⁻⁵³ relative, whenever every intermediate is a normal double, which holds
    once the true tail exceeds 10⁻³⁰⁰); round 131 (D3dz-1): a double tail below TAIL_FLOOR — r^K underflowing to 0 or to a
    subnormal for r < 10⁻¹⁰ (Λ), which the shipped certificate meets — is raised to TAIL_FLOOR, so that the radius added
    is always at least the true remainder (a NaN tail stays NaN: the ball then decides nothing)."""
    r = mag(x)
    t0, t1 = tail0(r), tail1(r)
    t0, t1 = (TAIL_FLOOR if t0 < TAIL_FLOOR else t0), (TAIL_FLOOR if t1 < TAIL_FLOOR else t1)
    cs = [coef(k) for k in range(K)]
    xv = val(x)
    p = arb(0)
    for c in reversed(cs):
        p = p * xv + c
    p = p + arb(0, t0)
    if not isinstance(x, D):
        return p
    dp = arb(0)
    for k in range(K - 1, 0, -1):
        dp = dp * xv + k * cs[k]
    dp = dp + arb(0, t1)
    return D(p, dp * x.d)


_FACT = [arb(math.factorial(k)) for k in range(80)]


def _geo(first, ratio):
    if not ratio < 1:
        raise ValueError('series tail does not converge')
    return first / (1 - ratio)


def S(w, K=16):
    """sinc(√w) = Σ (−w)^k/(2k+1)! (entire; w = θ²)."""
    return _series(w, lambda k: (-1) ** k / _FACT[2 * k + 1], K,
                   lambda r: _geo(r ** K / math.factorial(2 * K + 1), r / ((2 * K + 2) * (2 * K + 3))) * 1.0000001,
                   lambda r: _geo(K * r ** (K - 1) / math.factorial(2 * K + 1), (K + 1) * r / (K * (2 * K + 2) * (2 * K + 3))) * 1.0000001)


def Cv(w, K=16):
    """(1 − cos √w)/w = Σ (−w)^k/(2k+2)!."""
    return _series(w, lambda k: (-1) ** k / _FACT[2 * k + 2], K,
                   lambda r: _geo(r ** K / math.factorial(2 * K + 2), r / ((2 * K + 3) * (2 * K + 4))) * 1.0000001,
                   lambda r: _geo(K * r ** (K - 1) / math.factorial(2 * K + 2), (K + 1) * r / (K * (2 * K + 3) * (2 * K + 4))) * 1.0000001)


def Lam(v, K=30):
    """(v − log(1 + v))/v² = Σ (−v)^k/(k + 2), |v| ≤ 1/2."""
    if mag(v) > 0.5:
        raise ValueError('|v| > 1/2')
    return _series(v, lambda k: arb((-1) ** k) / (k + 2), K,
                   lambda r: r ** K / ((K + 2) * (1 - r)) * 1.0000001,
                   lambda r: r ** (K - 1) / (1 - r) * 1.0000001)


def H1(d, K=26):
    """(h(1 + δ) − δ²)/δ³ = Σ (−1)^{k+1} δ^k/(k + 3), |δ| ≤ 1/2 (h(1 + δ) = δ + δ²/2 − log(1 + δ))."""
    if mag(d) > 0.5:
        raise ValueError('|δ| > 1/2')
    return _series(d, lambda k: arb((-1) ** (k + 1)) / (k + 3), K,
                   lambda r: r ** K / ((K + 3) * (1 - r)) * 1.0000001,
                   lambda r: r ** (K - 1) / (1 - r) * 1.0000001)


# -------------------------------------------------------------------------------------------- the bound, for validation
def E_centred(z, eta):
    """`cert_E` in Arb (z, η acb): log|1 − z²| + Re(c z²) + C(σ), the exact centred form of part V."""
    c = eta * eta / 4
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    C = arb(0)
    if r > 1:
        a = pos(abs(sg.real) - 1)
        C = (r * r - 1) / 2 - r.log() - a * a
    elif not (r <= 1):
        raise ValueError('|σ| straddles 1')
    return abs(1 - z * z).log() + (c * z * z).real + C


def E_def(z, eta):
    """The bound from F_η's definition (D61 (iii)(b)): F_η(σ) = log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + C(σ),
    C = Λ⁺(σ) = M(Re σ) − ½ − log|σ| − Re σ²/2 outside the unit disc; E = F_η(σ) − F_η(η/2)."""
    def F(s):
        t = eta - s
        v = abs(t).log() + (t * t).real / 2 + abs(s).log() + (s * s).real / 2
        r = abs(s)
        if r > 1:
            a = abs(s.real)
            M = a * a if a <= 1 else 2 * a - 1
            v += M - arb(1) / 2 - r.log() - (s * s).real / 2
        return v
    return F(eta / 2 * (1 + z)) - F(eta / 2)


def path_point(eps, th, kind, X=None, lam=None):
    """The path's point z (acb) for η = (2 − ε)e^{iθ} (eps, th arb): kind 'J' (X = t, λ), 'F' (Y = κX), 'C' (Y = θ)."""
    om = 1 - eps / 2
    A = om * om * acb(0, 2 * th).exp() - 1
    kap = A.imag / abs(A.real)
    if kind == 'J':
        t = eps * XB
        return acb(t, lam * kap * t)
    return acb(X, kap * X) if kind == 'F' else acb(X, th)


# --------------------------------------------------------------------------------------------------- form A: (u, q, ρ)
def fA(u, q, rho, br, lam=None):
    """E/X² at X = ρ, ε = uρ, θ² = qρ on the branch br ('J': Y = λκX at u = 1/XB, 'F': Y = κX, 'C': Y = θ); u, q, λ
    arb, ρ arb or dual. Exact (an identity) but for (|Re σ| − 1)₊ taken as (Re σ − 1)₊ (an upper bound; Re σ ≈ 1 here)."""
    w = q * rho                      # θ²
    eps = u * rho
    om = 1 - eps / 2
    om2 = om * om
    S1, S4, Cw = S(w), S(4 * w), Cv(w)
    alpha = u * (1 - eps / 4) + 2 * om2 * q * S1 * S1       # |Re A|/X
    bh = 2 * om2 * S4                                      # Im A/θ
    m = bh / alpha if br == 'F' else (lam * bh / alpha if br == 'J' else D(1) if isinstance(rho, D) else arb(1))
    Q = q * m * (alpha * m - 2 * bh) - alpha * rho          # Re(A z²)/X², Y = θm
    zz = rho + q * m * m                                    # |z|²/X
    vr = -2 * (rho - q * m * m) + rho * zz * zz             # v/X, v = |1 − z²|² − 1
    L = zz * zz / 2 - vr * vr * Lam(rho * vr) / 2           # (log|1 − z²| + Re z²)/X² = |z|⁴/2 − v²Λ(v)/2 over X²
    a1z = (1 + rho * (2 + zz)).sqrt()                       # |1 + z|
    Dh = (om2 * (2 + zz) - u * (1 - eps / 4)) / (om * a1z + 1)   # (|σ| − 1)/X
    Ph = -u / 2 + om * (1 - (1 + rho) * q * Cw - q * m * S1)    # (Re σ − 1)/X
    Dp = Dh.pos() if isinstance(Dh, D) else pos(Dh)
    Pp = Ph.pos() if isinstance(Ph, D) else pos(Ph)
    Cc = Dp * Dp * (1 + rho * Dp * H1(rho * Dp)) - Pp * Pp
    return Q + L + Cc


def _min_hi(a, b):
    """Two upper bounds of the same quantity: the smaller upper end (only .hi is used)."""
    return a if a.hi < b.hi else b


def fA0(u, q, br, lam=None):
    """fA at ρ = 0 in closed form (the weight-4 limit Ψ̃_ψ(x, path)/x² with u = 1/x, q = ψ²/x; Lemma VC's algebra), on
    endpoint intervals u, q, λ (`I`): an upper bound on the box. C₀ = D₀₊² − P₀₊² is bounded by D₀₊² − (P₀ rigorously ≥ 0)₊²
    and, by the convexity of x₊², by 2D₀₊(D₀ − P₀), which carries the factor q (D₀ − P₀ = 2q on the cap, q(2/s² + ½ + 2/s)
    on the floor); the smaller upper end is used."""
    s = u + 2 * q
    if br == 'C':
        P0 = 1 - u * arb(0.5) - q * arb(1.5)
        if P0.lo >= 0:
            return -(u * q) - q.sq() * arb(0.5)                     # D₀² − P₀² = 2q(2 − u − q): exact
        D0 = 1 - u * arb(0.5) + q * arb(0.5)
        b1 = q * (s - 4) - q.sq() * arb(0.5) + D0.pos().sq() - I(P0.pos().lo).sq()
        b2 = q * (s - 4 + 4 * D0.pos()) - q.sq() * arb(0.5)
        return _min_hi(b1, b2)
    if br == 'F':
        D0 = 1 - u * arb(0.5) + 2 * q / s.sq()
        base = -(4 * q / s) - 8 * q.sq() / s ** 4
        b1 = base + D0.pos().sq()
        b2 = -8 * q.sq() / s ** 4 + q * (-4 / s + 2 * D0.pos() * (2 / s.sq() + arb(0.5) + 2 / s))
        return _min_hi(b1, b2)
    D0 = 1 - u * arb(0.5) + 2 * q * lam.sq() / s.sq()
    return -((4 * q / s) * (lam * (2 - lam))) - 8 * lam ** 4 * q.sq() / s ** 4 + D0.pos().sq()


def branchesA(u, q, rho_b):
    """Which branches of the floor/cap switch (m = bh/α vs 1) occur on the box: compare α and bh over ρ ∈ [0, ρ_b]."""
    rho = hull(0, rho_b)
    w = q * rho
    om = 1 - u * rho / 2
    alpha = u * (1 - u * rho / 4) + 2 * om * om * q * S(w) ** 2
    bh = 2 * om * om * S(4 * w)
    d = alpha - bh
    if d > 0:
        return ('F',)
    if d < 0:
        return ('C',)
    return ('F', 'C')


# --------------------------------------------------------------------------------------------------- form B: (u, p, ϑ)
def fB(u, p, th2, br, lam=None):
    """E/X² at X = pϑ, ε = uX, θ² = ϑ (ϑ arb or dual) on 'J' (u = 1/XB) or 'F'; the (Re σ − 1)₊² term dropped (an upper
    bound; Re σ < 1 on these branches)."""
    X = p * th2
    eps = u * X
    om = 1 - eps / 2
    om2 = om * om
    S1, S4 = S(th2), S(4 * th2)
    alpha = u * p * (1 - eps / 4) + 2 * om2 * S1 * S1           # |Re A|/θ²
    bh = 2 * om2 * S4                                           # Im A/θ
    mt = bh / alpha if br == 'F' else lam * bh / alpha          # Y = θ p m̃ (Y/X = m̃/θ)
    Q = mt * (alpha * mt - 2 * bh) - th2 * alpha
    zz = p * (th2 + mt * mt)                                    # |z|²/X
    vr = -2 * p * (th2 - mt * mt) + X * zz * zz                 # v/X
    L = zz * zz / 2 - vr * vr * Lam(X * vr) / 2
    a1z = (1 + X * (2 + zz)).sqrt()
    Dh = (om2 * (2 + zz) - u * (1 - eps / 4)) / (om * a1z + 1)
    Dp = Dh.pos() if isinstance(Dh, D) else pos(Dh)
    Cc = Dp * Dp * (1 + X * Dp * H1(X * Dp))
    return Q + L + Cc


def fB0(u, p, br, lam=None):
    """fB at ϑ = 0 in closed form (fA0's floor/junction with q = 1/p), on endpoint intervals."""
    s = u * p + 2
    if br == 'F':
        return -(4 / s) - 8 * p.sq() / s ** 4 + (1 - u * arb(0.5) + 2 * p / s.sq()).pos().sq()
    return -((4 / s) * (lam * (2 - lam))) - 8 * lam ** 4 * p.sq() / s ** 4 + (1 - u * arb(0.5) + 2 * p * lam.sq() / s.sq()).pos().sq()


def Bfloor_ok(u, p, th2b):
    """Rigorous check that form B's boxes are on the floor: bh·p < α (κX < θ) over the box."""
    th2 = hull(0, th2b)
    X = p * th2
    om = 1 - u * X / 2
    alpha = u * p * (1 - u * X / 4) + 2 * om * om * S(th2) ** 2
    return bool(2 * om * om * S(4 * th2) * p < alpha)


# ------------------------------------------------------------------------------------------------------- the box bound
def mvbound(f, f0, ivs, lo, hi, branch):
    """Upper bound data for f on the expansion variable's interval [lo, hi] (lo ≥ 0), the other arguments the float
    intervals ivs (balls for f, endpoint intervals for the closed form f0). Returns (ok, bound, route, sup f′): with
    lo = 0, ok when f(0) ≤ 0 and sup f′ < 0 (route 'd': then f < 0 on (0, hi]) or f(0) + hi·max(sup f′, 0) < 0 (route
    'v'); with lo > 0, f(lo) is bounded both directly and by f(0) + lo·sup f′([0, lo]), and the box passes when
    f(lo) + (hi − lo)·max(sup f′([lo, hi]), 0) < 0."""
    bl = [hull(*iv) for iv in ivs]
    il = [I(*iv) for iv in ivs]
    a, b = bl[:2], bl[2:]
    ai, bi = il[:2], il[2:]
    try:
        der = f(*a, D(hull(lo, hi), 1), branch, *b).d.upper()
        dpos = der if der > 0 else (arb(0) if der <= 0 else None)      # None: a NaN derivative, never accepted
        if dpos is None:
            return False, math.inf, 'x', math.inf
        if lo == 0:
            v0 = f0(*ai, branch, *bi).hi
            if v0 <= 0 and der < 0:
                return True, float(v0), 'd', float(der)
            bd = v0 + hi * dpos
            return bool(bd < 0), float(bd), 'v', float(der)
        fl = f(*a, arb(lo), branch, *b).upper()
        d0 = f(*a, D(hull(0, lo), 1), branch, *b).d.upper()
        fl2 = f0(*ai, branch, *bi).hi + lo * d0
        fl = fl2 if fl2 < fl else fl                                # (a NaN fl2 is never taken)
        bd = fl + (arb(hi) - arb(lo)) * dpos                            # round 121 (D3dp-5): the width in balls
        return bool(bd < 0), float(bd), 'v', float(der)
    except (ValueError, ZeroDivisionError):
        return False, math.inf, 'x', math.inf


def bb(name, f, f0, roots, split, admissible, extra_check=None, max_boxes=2_000_000, min_w=1e-9, log=None, keep=None):
    """Branch and bound: boxes are dicts of intervals; `admissible(box)` prunes boxes outside the domain (returns the
    box with its expansion interval clipped, or None); a box passes when every branch it meets passes `mvbound`. With
    keep (a list), every accepted box's intervals and its branches' (branch, route, bound) are appended (round 121, for
    `soundness`)."""
    t0 = time.time()
    stack = list(reversed(roots))
    nok, nskip, nfail, worst, fails, nev = 0, 0, 0, -math.inf, [], 0
    worst_d, n_d = -math.inf, 0
    smallest = {}
    tick = 0
    while stack:
        tick += 1
        if tick % 50000 == 0:
            print(f"  [{name}] {tick} boxes popped, {nok} ok, {nfail} failed, {len(stack)} on the stack, {time.time() - t0:.0f}s",
                  file=sys.stderr, flush=True)
        box = stack.pop()
        box = admissible(box)
        if box is None:
            nskip += 1
            continue
        ok, w, brs, wd, nd, recs = True, -math.inf, box['branches'](), -math.inf, 0, []
        for br in brs:
            nev += 1
            o, v, route, der = box['bound'](br)
            if not o:
                ok = False
                break
            recs.append((br, route, v))
            if route == 'd':
                wd, nd = max(wd, der), nd + 1
            else:
                w = max(w, v)
        if ok and (extra_check is None or extra_check(box)):
            nok += 1
            worst = max(worst, w)
            worst_d = max(worst_d, wd)
            n_d += nd
            for k, (a, b) in box['iv'].items():
                smallest[k] = min(smallest.get(k, math.inf), b - a)
            if keep is not None:
                keep.append((dict(box['iv']), recs))
            continue
        kids = split(box)
        if kids is None:
            nfail += 1
            if len(fails) < 20:
                fails.append({k: v for k, v in box['iv'].items()})
            continue
        stack += list(reversed(kids))
        if nok + nfail + len(stack) > max_boxes:
            nfail += len(stack)
            stack = []
    out = dict(piece=name, boxes=nok, skipped=nskip, failed=nfail, worst=worst, worst_der=worst_d, n_der=n_d, evals=nev, smallest=smallest,
               fails=fails, seconds=round(time.time() - t0, 1))
    if log:
        log(json.dumps(out))
    return out


def _mk_box(iv, branches, bound):
    return dict(iv=iv, branches=branches, bound=bound)


# ------------------------------------------------------------------------------------------------ the pieces A and B
def pieceA(kind, eps0=EPS0, delta0=DELTA0, log=None, nu=8, nq=10, keep=None):
    """Form A over u ∈ [0, 4] (the floor and the cap) or u = 4 (the junction, λ ∈ [0, 1]), q ∈ [0, QA], ρ ∈ [0, ρmax]
    with ρ ≤ δ₀ and uρ ≤ ε₀ (ε ≤ ε₀). kind: 'FC' or 'J'."""
    def make(iv):
        u = hull(*iv['u'])
        q = hull(*iv['q'])
        lam = hull(*iv['l']) if 'l' in iv else None
        rb = iv['r'][1]
        if kind == 'J':
            brs = lambda: ('J',)
        else:
            brs = lambda: branchesA(u, q, rb)
        bound = lambda br: mvbound(fA, fA0, (iv['u'], iv['q']) + ((iv['l'],) if 'l' in iv else ()), iv['r'][0], rb, br)
        return _mk_box(iv, brs, bound)

    def adm(box):
        iv = dict(box['iv'])
        ua = iv['u'][0]
        rmax = delta0 if kind != 'J' else eps0 * float(XB)
        if ua > 0:
            rmax = min(rmax, eps0 / ua)
        rmax = rmax * (1 + 1e-12)                  # outward past the float divisions
        if iv['r'][0] > rmax:
            return None
        iv['r'] = (iv['r'][0], min(iv['r'][1], rmax))
        return make(iv)

    def split(box):
        iv = box['iv']
        cand = []
        for k, (a, b) in iv.items():
            if k == 'u' and kind == 'J':
                continue
            wdt = b - a
            if wdt < min_w(k):
                continue
            scale = {'u': 0.25, 'q': 0.25, 'l': 0.25, 'r': 4 * max(iv['r'][1], 1e-12)}[k]
            cand.append((wdt / scale, k))
        if not cand:
            return None
        _, k = max(cand)
        a, b = iv[k]
        m = (a + b) / 2
        return [make(dict(iv, **{k: (a, m)})), make(dict(iv, **{k: (m, b)}))]

    def min_w(k):
        return {'u': 1e-7, 'q': 1e-9, 'l': 1e-7, 'r': 1e-12}[k]

    roots = []
    if kind == 'J':
        for j in range(nq):
            for i in range(4):
                roots.append(make({'u': (UJ, UJ), 'q': (QA * j / nq, QA * (j + 1) / nq), 'l': (i / 4, (i + 1) / 4),
                                   'r': (0.0, eps0 * float(XB))}))
    else:
        for i in range(nu):
            for j in range(nq):
                roots.append(make({'u': (UJ * i / nu, UJ * (i + 1) / nu), 'q': (QA * j / nq, QA * (j + 1) / nq), 'r': (0.0, delta0)}))
    return bb('A-' + kind, fA, fA0, roots, split, adm, log=log, keep=keep)


def pieceB(kind, eps0=EPS0, delta0=DELTA0, log=None, nu=8, np_=8, nt=16, keep=None):
    """Form B over u ∈ [0, 4] (floor) or u = 4 (junction), p ∈ [0, PB], ϑ ∈ [0, (π/2)²] with X = pϑ ≤ δ₀, uX ≤ ε₀."""
    thm = up(TH_MAX2)

    def make(iv):
        u = hull(*iv['u'])
        p = hull(*iv['p'])
        lam = hull(*iv['l']) if 'l' in iv else None
        brs = (lambda: ('J',)) if kind == 'J' else (lambda: ('F',))
        bound = lambda br: mvbound(fB, fB0, (iv['u'], iv['p']) + ((iv['l'],) if 'l' in iv else ()), iv['t'][0], iv['t'][1], br)
        return _mk_box(iv, brs, bound)

    def adm(box):
        iv = dict(box['iv'])
        pa, ua = iv['p'][0], iv['u'][0]
        tmax = thm
        if pa > 0:
            tmax = min(tmax, delta0 / pa)
            if ua > 0:
                tmax = min(tmax, eps0 / (ua * pa))
        tmax = tmax * (1 + 1e-12)                  # outward past the float divisions
        if iv['t'][0] > tmax:
            return None
        iv['t'] = (iv['t'][0], min(iv['t'][1], tmax))
        return make(iv)

    def split(box):
        iv = box['iv']
        cand = []
        for k, (a, b) in iv.items():
            if k == 'u' and kind == 'J':
                continue
            if b - a < {'u': 1e-7, 'p': 1e-9, 'l': 1e-7, 't': 1e-12}[k]:
                continue
            scale = {'u': 0.25, 'p': 0.05, 'l': 0.25, 't': 0.05}[k]
            cand.append(((b - a) / scale, k))
        if not cand:
            return None
        _, k = max(cand)
        a, b = iv[k]
        m = (a + b) / 2
        return [make(dict(iv, **{k: (a, m)})), make(dict(iv, **{k: (m, b)}))]

    def chk(box):
        iv = box['iv']
        return Bfloor_ok(hull(*iv['u']), hull(*iv['p']), iv['t'][1])

    roots = []
    tg = [thm * (k / nt) ** 2 for k in range(nt + 1)]
    for j in range(np_):
        for k in range(nt):
            if kind == 'J':
                for i in range(4):
                    roots.append(make({'u': (UJ, UJ), 'p': (PB * j / np_, PB * (j + 1) / np_), 'l': (i / 4, (i + 1) / 4), 't': (tg[k], tg[k + 1])}))
            else:
                for i in range(nu):
                    roots.append(make({'u': (UJ * i / nu, UJ * (i + 1) / nu), 'p': (PB * j / np_, PB * (j + 1) / np_), 't': (tg[k], tg[k + 1])}))
    return bb('B-' + kind, fB, fB0, roots, split, adm, extra_check=chk, log=log, keep=keep)


# ------------------------------------------------------------------------------------------------------- the outer piece
def Ef(z, eta):
    """Float centred bound (path choice only)."""
    c = eta * eta / 4
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    C = 0.0 if r <= 1 else (r * r - 1) / 2 - math.log(r) - max(abs(sg.real) - 1, 0.0) ** 2
    return math.log(abs(1 - z * z)) + (c * z * z).real + C


def gradf(z, eta):
    c = eta * eta / 4
    v = (2 * z * (c - 1 / (1 - z * z))).conjugate()
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    if r > 1:
        a = sg.real
        gs = (r - 1 / r) * sg / r - (2 * (abs(a) - 1) * math.copysign(1, a) if abs(a) > 1 else 0.0)
        v += gs * (eta / 2).conjugate()
    return v


def flow(eta, z0, h=2e-3, every=8, stop=0.01, nmax=4000):
    """Normalised steepest descent of the float bound from z0 until |1 − z| < stop, decimated (None if it fails)."""
    z, pts = complex(z0), [complex(z0)]
    for _ in range(nmax):
        if abs(1 - z) < stop:
            out = pts[::every]
            if out[-1] != pts[-1]:
                out.append(pts[-1])
            return out
        g = gradf(z, eta)
        if abs(g) == 0:
            return None
        z = z - h * g / abs(g)
        pts.append(z)
    return None


def th_ball(f1, f2):
    """The angle interval [f1, f2] as a ball, the top edge taken as the true π/2 where it is the float π/2."""
    return hull(f1, arb.pi() / 2) if f2 == math.pi / 2 else hull(f1, f2)


def seg_bound(P, Q, R1, R2, f1, f2):
    """Rigorous upper bound of E over the segment [P, Q] (floats) × the η-box R1 ≤ |η| ≤ R2, f1 ≤ arg η ≤ f2: the mean
    value form in (s, R, θ) from the centre (z = (P + Q)/2, R_c, θ_c), every partial derivative enclosed on the box."""
    Pz, Qz = acb(P.real, P.imag), acb(Q.real, Q.imag)
    M, H = (Pz + Qz) / 2, (Qz - Pz) / 2
    Rc, Tc = (arb(R1) + arb(R2)) / 2, (arb(f1) + arb(f2)) / 2
    Tb = th_ball(f1, f2)
    rR = (arb(R2) - arb(R1)) / 2
    rT = (Tb.upper() - Tc).upper()
    rT2 = (Tc - arb(f1)).upper()
    rT = rT if rT > rT2 else rT2
    eta_c = Rc * acb(0, Tc).exp()
    E0 = E_centred_safe(M, eta_c)
    if E0 is None:
        return math.inf
    s = arb(0, 1)
    Rb = hull(R1, R2)
    z = M + s * H
    eta = Rb * acb(0, Tb).exp()
    c = eta * eta / 4
    w = 1 - z * z
    if not (abs(w) > 0):
        return math.inf
    gp = 2 * z * (c - 1 / w)
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    hp = pos(r - 1 / r)
    a = sg.real
    ta = pos(abs(a) - 1)
    if ta > 0 or not (ta <= 0):
        if a > 0:
            sgn = 1
        elif a < 0:
            sgn = -1
        else:
            return math.inf
        Gs = hp * sg / r - 2 * sgn * ta
    else:
        if not (r > 0):                                             # round 121 (D3dp-5): |σ|'s ball meets 0
            return math.inf
        Gs = hp * sg / r
    Gc = Gs.conjugate()
    dS = (gp * H).real + (Gc * eta / 2 * H).real
    e1 = acb(0, Tb).exp()
    dR = (z * z * Rb * e1 * e1 / 2).real + (Gc * e1 * (1 + z) / 2).real
    dT = (2 * acb(0, 1) * c * z * z).real + (Gc * acb(0, 1) * sg).real
    b = E0 + abs(dS) + abs(dR) * rR + abs(dT) * rT
    return up(b)


def E_centred_safe(z, eta):
    try:
        return E_centred(z, eta)
    except ValueError:
        c = eta * eta / 4
        sg = eta / 2 * (1 + z)
        r = abs(sg)
        h = (r * r - 1) / 2 - r.log()
        return abs(1 - z * z).log() + (c * z * z).real + hull(0, pos(h).upper())


def handoff_range(R1, R2, f1, f2, delta0):
    """Rigorous range [Y_lo, Y_hi] of the handoff's Im z = min(κδ₀, θ) over the η-box, and its value at the box centre."""
    Rb, Tb = hull(R1, R2), th_ball(f1, f2)
    om = Rb / 2
    A = om * om * acb(0, 2 * Tb).exp() - 1
    kap = A.imag / abs(A.real)
    kl, kh = float((kap * delta0).lower()), float((kap * delta0).upper())
    th_hi = up(Tb)
    if not (math.isfinite(kl) and math.isfinite(kh)):      # round 121 (D3dp-1): κ's ball NaN/∞ (Re A's ball contains 0)
        lo, hi = 0.0, th_hi                                 # 0 ≤ min(κδ₀, θ) ≤ θ ≤ θ_hi
    else:
        lo = min(math.nextafter(kl, -math.inf), f1)
        hi = max(min(math.nextafter(kh, math.inf), th_hi), 0.0)
        lo = max(min(lo, hi), 0.0)
    assert math.isfinite(lo) and math.isfinite(hi) and lo <= hi
    Rc, Tc = (R1 + R2) / 2, (f1 + f2) / 2
    Ac = (Rc / 2) ** 2 * cmath_exp(2j * Tc) - 1
    yc = min(Ac.imag / abs(Ac.real) * delta0, Tc)
    yc = min(max(yc, lo), hi)
    return lo, hi, yc


def cmath_exp(z):
    import cmath
    return cmath.exp(z)


def outer_box(R1, R2, f1, f2, delta0, min_seg=1e-6):
    """The outer piece on one η-box: the connector [δ₀ + iY_lo, δ₀ + iY_hi] and the centre's descent polyline from
    δ₀ + iY_c to |1 − z| ≤ END_R, every segment bisected until its rigorous bound is negative. Returns (ok, segments,
    worst bound, the polyline's end, the joint: lo, hi, y_c, the polyline's points and the certified pieces in the order
    certified); round 121 (D3dp-1): the connector is always certified (a point when lo = hi), never dropped; round 122
    (D3dq-1): the pieces actually certified are returned, so that `handoff_check` rebuilds the cover from them."""
    lo, hi, yc = handoff_range(R1, R2, f1, f2, delta0)
    if not (lo <= hi):
        return False, 0, math.inf, None, None
    eta_c = (R1 + R2) / 2 * cmath_exp(1j * (f1 + f2) / 2)
    pts = flow(eta_c, complex(delta0, yc))
    if pts is None:
        return False, 0, math.inf, None, None
    end = pts[-1]
    if not (abs(1 - acb(end.real, end.imag)) < END_R):
        return False, 0, math.inf, end, None
    segs = [(complex(delta0, lo), complex(delta0, hi))]
    segs += [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    n, worst, cert = 0, -math.inf, []
    stack = segs[::-1]
    while stack:
        P, Q = stack.pop()
        b = seg_bound(P, Q, R1, R2, f1, f2)
        if b < 0:
            n, worst = n + 1, max(worst, b)
            cert.append((P, Q))
            continue
        if abs(Q - P) < min_seg:
            return False, n, b, end, None
        m = (P + Q) / 2
        stack += [(m, Q), (P, m)]
    return True, n, worst, end, dict(lo=lo, hi=hi, yc=yc, pts=pts, cert=cert)


def piece_outer(eps0=EPS0, delta0=DELTA0, nT=64, depth=22, log=None):
    """The outer piece over 2 − ε₀ ≤ |η| ≤ 2, 0 ≤ arg η ≤ π/2 (the top edge's float π/2 rechecked below)."""
    t0 = time.time()
    R1, R2 = 2 - eps0, 2.0
    stack = [(R1, R2, (math.pi / 2) * j / nT, (math.pi / 2) * (j + 1) / nT, 0) for j in range(nT)][::-1]
    nb, ns, worst, fails, sides, acc, joints = 0, 0, -math.inf, [], [math.inf, math.inf], [], []
    while stack:
        a, b, f1, f2, dep = stack.pop()
        ok, n, w, end, conn = outer_box(a, b, f1, f2, delta0)
        if ok:
            nb, ns, worst = nb + 1, ns + n, max(worst, w)
            acc.append([a, b, f1, f2, conn['lo'], conn['hi']])
            joints.append(conn)
            sides = [min(sides[0], b - a), min(sides[1], f2 - f1)]
            continue
        if dep >= depth:
            fails.append((a, b, f1, f2))
            continue
        if (b - a) > 2 * (f2 - f1):
            m = (a + b) / 2
            stack += [(m, b, f1, f2, dep + 1), (a, m, f1, f2, dep + 1)]
        else:
            m = (f1 + f2) / 2
            stack += [(a, b, m, f2, dep + 1), (a, b, f1, m, dep + 1)]
    top = 'the top edge at the true π/2 (th_ball)'
    out = dict(piece='O', boxes=nb, segs=ns, worst=worst, failed=len(fails), fails=fails[:20], smallest_dR=sides[0],
               smallest_dth=sides[1], top=top, accepted=acc, R1=R1, joints=joints, seconds=round(time.time() - t0, 1))
    if log:                                     # round 131 (D3dz-4): the joints (complex points, in memory for `handoff_check`) not logged
        log(json.dumps({k: v for k, v in out.items() if k != 'joints'}))
    return out


def piece_end():
    """Part V's end piece at R2 = 2: from |1 − z_n| ≤ END_R straight to z = 1, log|1 − z²| ≤ log END_R + log(2 + END_R),
    Re(cz²) ≤ (R2²/4)(1 + END_R)², C ≤ h((R2/2)(2 + END_R))."""
    e, R = arb(END_R), arb(2)
    r = R / 2 * (2 + e)
    b = e.log() + (2 + e).log() + R * R / 4 * (1 + e) ** 2 + (r * r - 1) / 2 - r.log()
    return dict(piece='end', bound=up(b), ok=bool(b < 0))


# -------------------------------------------------------------------------------------------------------- validation
VALIDATE_TOL = ((0.005, 1e-20), (0.11, 1e-10))   # round 131 (D3dz-4): (ε₀, the bar on validate's deviations), ε₀ increasing


def validate_tol(eps0):
    """The bar on `validate`'s A_vs_centred, A_vs_def and B_vs_centred at ε₀ (round 131, D3dz-4) — a function of ε₀ through the
    series' tails. The deviations are of midpoints: the forms' balls are the truncated series (S, Cv to 16 terms, Λ to 30, H₁ to
    26) plus their rigorous tails as radius, so their midpoints miss the exact value by up to those tails at the arguments
    validate's points reach, and these grow with ε₀ (the seeded ε are 10^uniform(−9, log₁₀ε₀)): at ε₀ = 0.005 |δ| = |σ| − 1 ≤ 0.14
    and |v| ≤ 0.11 (H₁'s tail 2.9·10⁻²⁴; the deviations 7.6·10⁻²³), at ε₀ = 0.11 |δ| ≤ 0.42 and |v| ≤ 0.33 (H₁'s tail 8.4·10⁻¹²;
    1.1·10⁻¹²). The bar is VALIDATE_TOL's at the least tabulated ε₀′ ≥ ε₀ (the arguments, and so the tails, increase with ε₀),
    about two decades above the deviation measured there; None — no bar known, the seam check fails — beyond the table. The
    closed forms' deviations (A0, B0: ρ = 0, no ε) have the ε₀-free bar 10⁻⁴⁰."""
    for e, tol in VALIDATE_TOL:
        if eps0 <= e:
            return tol
    return None


def validate(log=print, n=200, seed=120, eps0=None):
    """The forms against `E_centred` (part V's centred form) and E_def (F_η's definition) at thin parameters, and fA0/fB0
    against fA/fB at ρ = 0 (ϑ = 0). Returns the largest relative deviations (their bar: `validate_tol`). The seeded ε range
    over [10⁻⁹, ε₀], ε₀ = EPS0 unless given (round 131: part VA011 passes 0.11)."""
    import random
    rnd = random.Random(seed)
    ctx.prec = 200
    eps0 = EPS0 if eps0 is None else eps0
    worst = dict(A_vs_centred=0.0, A_vs_def=0.0, B_vs_centred=0.0, centred_vs_def=0.0, A0=0.0, B0=0.0)
    cnt = dict(points=0, A=0, B=0, A0=0, B0=0)
    cases = []
    for _ in range(n):
        eps = arb(10.0 ** rnd.uniform(-9, math.log10(eps0)))
        reg = rnd.choice(('phi', 'psi', 'psi', 'th'))
        if reg == 'phi':
            th = eps * rnd.uniform(0, 20)
        elif reg == 'psi':
            th = eps.sqrt() * 10 ** rnd.uniform(-2, 1.3)
        else:
            th = arb(rnd.uniform(0.0, 1.5))
        th = arb(th)
        cases.append((eps, th))
    for eps, th in cases:
        om = 1 - eps / 2
        eta = (2 - eps) * acb(0, th).exp()
        A = om * om * acb(0, 2 * th).exp() - 1
        kap = A.imag / abs(A.real)
        t = eps * XB
        Xc = th / kap if th > 0 else arb(1)
        pts = [('J', None, arb(l)) for l in (0.1, 0.5, 1.0)]
        for f in (1.5, 3, 10, 100):
            X = t * f
            if X < Xc:
                pts.append(('F', X, None))
        for f in (1.2, 3, 30):
            X = max(Xc, t) * f
            if X < DELTA0:
                pts.append(('C', X, None))
        for br, X, lam in pts:
            z = path_point(eps, th, br, X, lam)
            Xr = z.real
            try:
                ec = E_centred(z, eta)
            except ValueError:
                continue
            ed = E_def(z, eta)
            cnt['points'] += 1
            u, q = eps / Xr, th * th / Xr
            rel = lambda a, b: float((abs(a - b) / abs(b)).mid()) if b != 0 else 0.0
            worst['centred_vs_def'] = max(worst['centred_vs_def'], rel(ec, ed))
            if float(q.mid()) <= QA:
                try:
                    fa = fA(u, q, Xr, br, lam) * Xr * Xr
                except ValueError:                                      # round 122 (D3dq-5): outside the forms' domain (ε₀ = 0.11)
                    cnt['skipped'] = cnt.get('skipped', 0) + 1
                    continue
                worst['A_vs_centred'] = max(worst['A_vs_centred'], rel(fa, ec))
                worst['A_vs_def'] = max(worst['A_vs_def'], rel(fa, ed))
                cnt['A'] += 1
            if th > 0 and br in ('J', 'F'):
                p = Xr / (th * th)
                if float(p.mid()) <= 1:
                    fb = fB(u, p, th * th, br, lam) * Xr * Xr
                    worst['B_vs_centred'] = max(worst['B_vs_centred'], rel(fb, ec))
                    cnt['B'] += 1
    for _ in range(n):
        u, q = arb(rnd.uniform(0, 4)), arb(rnd.uniform(0, QA))
        lam = arb(rnd.uniform(0, 1))
        for br in ('F', 'C', 'J'):
            uu = arb(UJ) if br == 'J' else u
            if br == 'F' and float((uu + 2 * q).mid()) < 2.05:
                continue
            if br == 'C' and float((uu + 2 * q).mid()) > 1.95:
                continue
            a = fA(uu, q, arb(0), br, lam)
            b = fA0(I(uu.mid()), I(q.mid()), br, I(lam.mid()))
            if not (a <= b.hi + arb(10) ** -50):
                raise AssertionError(('fA0 below fA(0)', br, a, b.hi))
            worst['A0'] = max(worst['A0'], float(abs(a - b.hi).mid()))
            cnt['A0'] += 1
        p = arb(rnd.uniform(0, 1))
        for br in ('F', 'J'):
            uu = arb(UJ) if br == 'J' else u
            a = fB(uu, p, arb(0), br, lam)
            b = fB0(I(uu.mid()), I(p.mid()), br, I(lam.mid()))
            if not (a <= b.hi + arb(10) ** -50):                       # round 121 (D3dp-4): asserted as A0 is
                raise AssertionError(('fB0 below fB(0)', br, a, b.hi))
            worst['B0'] = max(worst['B0'], float(abs(a - b.hi).mid()), float(abs(b.hi - b.lo).mid()))
            cnt['B0'] += 1
    ctx.prec = PREC
    out = dict(piece='validate', cases=len(cases), counts=cnt, **worst)
    log(json.dumps(out))
    return out


def _chain(pieces, A, B):
    """The certified pieces, in the order certified, form one chain from A to B (each piece's end the next one's start, exactly)."""
    return bool(pieces) and pieces[0][0] == A and pieces[-1][1] == B and all(p[1] == q[0] for p, q in zip(pieces, pieces[1:]))


def _y_enclosed(R1, R2, f1, f2, lo, hi, delta0, depth=40, budget=200_000, keep=None):
    """Rigorously, Y(η) = min(κ(η)δ₀, arg η) ∈ [lo, hi] for every η in the box (round 122, D3dq-1): on sub-boxes in balls, κ =
    ω²sin 2θ/((1 − ω²) + 2ω²sin²θ) — Im A/|Re A| without cancellation, ω = |η|/2 ≤ 1, |Re A| = 1 − ω²cos 2θ —, with sin 2θ ≥ 0
    and 1 − ω² ≥ 0 applied exactly (0 ≤ θ ≤ π/2, |η| ≤ 2), min of balls as [min of lowers, min of uppers], and where the
    denominator's ball meets 0 (only at η = 2) the enclosure 0 ≤ Y ≤ θ; the top edge at the true π/2; adaptive bisection, the
    corners where the connector's ends are attained needing the deepest boxes. Returns the number of sub-boxes used, or None if
    some sub-box stays undecided at the depth, the budget of sub-boxes is spent, or the box is not inside the sector 0 < |η| ≤ 2,
    0 ≤ arg η ≤ π/2 where the clamps hold (round 126, D3du-6). With keep (a list; round 135, I4: the certificate's export, not used
    otherwise), every sub-box proved (r1, r2, t1, t2) is appended."""
    if not (0 < R1 <= R2 <= 2.0 and 0.0 <= f1 <= f2 <= math.pi / 2):
        return None
    L, H = arb(lo), arb(hi)

    def ok(r1, r2, t1, t2):
        Rb = hull(r1, r2)
        Tb = hull(t1, arb.pi() / 2) if t2 == math.pi / 2 else hull(t1, t2)
        om2 = (Rb / 2) ** 2
        sn = Tb.sin()
        den = pos(1 - om2) + 2 * om2 * sn * sn
        t_lo, t_hi = arb(t1), (arb.pi() / 2).upper() if t2 == math.pi / 2 else arb(t2)   # θ's exact ends
        if not (den > 0):
            y_lo, y_hi = arb(0), t_hi                                   # 0 ≤ Y ≤ θ
        else:
            K = om2 * pos((2 * Tb).sin()) / den * delta0
            k_lo, k_hi = K.lower(), K.upper()
            k_lo = k_lo if k_lo > 0 else arb(0)                         # κ ≥ 0 exactly on 0 ≤ θ ≤ π/2 (a ball through 0 dips below)
            y_lo = k_lo if k_lo < t_lo else t_lo                        # Y = min(κδ₀, θ) ≥ min(inf κδ₀, θ₁)
            y_hi = k_hi if k_hi < t_hi else t_hi                        #                ≤ min(sup κδ₀, θ₂)
        return bool(y_lo >= L and y_hi <= H)

    stack, used = [(R1, R2, f1, f2, 0)], 0
    while stack:
        r1, r2, t1, t2, d = stack.pop()
        if ok(r1, r2, t1, t2):
            used += 1
            if keep is not None:
                keep.append((r1, r2, t1, t2))
            continue
        if d >= depth or used + len(stack) > budget:
            return None
        rm, tm = (r1 + r2) / 2, (t1 + t2) / 2
        stack += [(r1, rm, t1, tm, d + 1), (r1, rm, tm, t2, d + 1), (rm, r2, t1, tm, d + 1), (rm, r2, tm, t2, d + 1)]
    return used


def handoff_check(O, delta0=DELTA0, n=6, seed=121):
    """Rounds 121 (D3dp-2) and 122 (D3dq-1, -2): the joint between the inner pieces and the outer one, from what was certified.
    For every accepted η-box of the outer piece O: the connector [lo, hi] finite; Y(η) = min(κ(η)δ₀, arg η) ∈ [lo, hi] for every
    η in the box, proved in balls on sub-boxes (`_y_enclosed`), and checked at its corners, edges' midpoints, arg η = 0 and π/2
    and n seeded random η; and — when O carries the joints (`outer_box`'s certified pieces, in memory) — the pieces certified
    chaining the connector from δ₀ + i·lo to δ₀ + i·hi and the polyline from its start to its end through each of its points
    (a chain of certified segments: a path), y_c ∈ [lo, hi], the polyline starting at δ₀ + i·y_c and ending inside |1 − z| <
    END_R, where the end piece takes over (round 126, D3du-2). The boxes must tile 2 − ε₀ ≤ |η| ≤ 2, 0 ≤ arg η ≤ π/2: inside
    the sector, pairwise disjoint interiors, and the areas summing to the sector's exactly, in rationals (round 126, D3du-1:
    with containment and disjointness, equal measure of closed boxes is a tiling). Returns the counts and the failures."""
    import random
    rnd = random.Random(seed)
    acc = O['accepted']
    joints = O.get('joints')
    pts, bad, marg, subs = 0, [], math.inf, 0
    area = 0.0
    R1s = O.get('R1', 2 - EPS0)
    for k, (a, b, f1, f2, lo, hi) in enumerate(acc):
        area += (b - a) * (f2 - f1)
        if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
            bad.append((a, b, f1, f2, 'connector'))
            continue
        u = _y_enclosed(a, b, f1, f2, lo, hi, delta0)
        if u is None:
            bad.append((a, b, f1, f2, 'Y not enclosed'))
        else:
            subs += u
        if joints is not None:
            J = joints[k]
            pc = J['cert']
            con = [x for x in pc if x[0].real == delta0 and x[1].real == delta0 and lo <= x[0].imag <= hi and lo <= x[1].imag <= hi]
            ncon = 0
            while ncon < len(pc) and pc[ncon] in con and not (ncon > 0 and pc[ncon - 1][1] == complex(delta0, hi)):
                ncon += 1
            if not (lo <= J['yc'] <= hi):
                bad.append((a, b, f1, f2, 'y_c'))
            if not (J['pts'] and J['pts'][0] == complex(delta0, J['yc'])):
                bad.append((a, b, f1, f2, 'polyline start'))
            if not _chain(pc[:ncon], complex(delta0, lo), complex(delta0, hi)):
                bad.append((a, b, f1, f2, 'connector cover'))
            pl = pc[ncon:]
            verts = {x[0] for x in pl} | {x[1] for x in pl}
            if not (_chain(pl, J['pts'][0], J['pts'][-1]) and all(p_ in verts for p_ in J['pts'])):
                bad.append((a, b, f1, f2, 'polyline cover'))
            e_ = J['pts'][-1]
            if not (abs(1 - acb(e_.real, e_.imag)) < END_R):           # the last link: the end piece's start
                bad.append((a, b, f1, f2, 'polyline end'))
        top = arb.pi() / 2 if f2 == math.pi / 2 else arb(f2)
        thl = [arb(f1), top, (arb(f1) + top) / 2] + [arb(rnd.uniform(f1, f2)) for _ in range(n)]
        exact0 = lambda th: th == 0 or (f2 == math.pi / 2 and th is top)   # Im A = ω² sin 2θ = 0 exactly: κ = 0
        Rl = [arb(a), arb(b), (arb(a) + arb(b)) / 2] + [arb(rnd.uniform(a, b)) for _ in range(n)]
        for R, th in [(R, th) for R in Rl[:3] for th in thl[:3]] + list(zip(Rl[3:], thl[3:])):
            if exact0(th):
                Y = arb(0)                                              # (at θ = 0 and η = 2 the limit)
            else:
                om = R / 2
                A = om * om * acb(0, 2 * th).exp() - 1
                if not (A.real < 0):
                    bad.append((a, b, f1, f2, 'Re A'))
                    continue
                Y = A.imag / abs(A.real) * delta0
                if not (Y < th):                                        # min(κδ₀, θ) on balls: [min of lowers, min of uppers]
                    Y = th if th < Y else hull(Y.lower() if Y.lower() < th.lower() else th.lower(),
                                               Y.upper() if Y.upper() < th.upper() else th.upper())
            pts += 1
            m = min(float((Y - lo).lower()), float((hi - Y).lower()))
            if not (Y >= lo and Y <= hi):
                bad.append((a, b, f1, f2, float(R.mid()), float(th.mid())))
            marg = min(marg, m)
    eps0 = 2 - R1s
    rel = abs(area / (eps0 * math.pi / 2) - 1)
    from fractions import Fraction as Fr
    exact = (sum((Fr(x[1]) - Fr(x[0])) * (Fr(x[3]) - Fr(x[2])) for x in acc) == (Fr(2.0) - Fr(R1s)) * Fr(math.pi / 2))
    inside = all(R1s <= x[0] < x[1] <= 2.0 and 0.0 <= x[2] < x[3] <= math.pi / 2 for x in acc)
    disjoint = all(min(p[1], q[1]) <= max(p[0], q[0]) or min(p[3], q[3]) <= max(p[2], q[2])
                   for i, p in enumerate(acc) for q in acc[i + 1:])
    return dict(piece='handoff', boxes=len(acc), points=pts, sub_boxes=subs, joints=joints is not None, failed=len(bad), fails=bad[:10],
                min_margin=marg, area_rel=rel, area_exact=exact, inside=inside, disjoint=disjoint, tiled=bool(exact and inside and disjoint))


def soundness(name, f, keep, per=2, seed=121):
    """Round 121 (D3dp-4): the forms' value against the certified bound inside every accepted box of a piece (`bb`'s keep):
    per branch `per` seeded points and (rounds 122, 126) the corners over the box's non-degenerate coordinates at both ends of
    the expansion variable (its lower end when > 0, as on the path), each distinct point once; on route 'v' the form at the
    point (a ball at PREC bits) must lie below the box's bound, on route 'd' below 0. Returns the number of points, the
    largest excess (≤ 0 passes) and the failures."""
    import random
    rnd = random.Random(seed)
    ev = {'A': 'r', 'B': 't'}[name[0]]
    sec = {'A': 'q', 'B': 'p'}[name[0]]
    worst, bad, n, npts = -math.inf, 0, 0, 0
    import itertools
    for iv, recs in keep:
        e0, e1 = iv[ev]
        free = [k for k in ('u', sec, 'l') if k in iv and iv[k][1] > iv[k][0]]   # round 126 (D3du-5): the non-degenerate ones
        pts = set()
        for _ in range(per):
            pts.add(tuple((k, rnd.uniform(*iv[k]) if iv[k][1] > iv[k][0] else iv[k][0]) for k in ('u', sec, 'l') if k in iv)
                    + (('e', rnd.uniform(e0, e1) if e1 > e0 else e1),))
        for bits in itertools.product((0, 1), repeat=len(free)):     # corners, at both ends of the expansion variable (> 0)
            c = {k: iv[k][b_] for k, b_ in zip(free, bits)}
            for e_ in ([e0] if e0 > 0 else []) + [e1]:
                pts.add(tuple((k, c.get(k, iv[k][0])) for k in ('u', sec, 'l') if k in iv) + (('e', e_),))
        npts += sum(1 for pt in pts if dict(pt)['e'] > 0)          # round 127 (D3dv-4): distinct points, each branch evaluated
        for br, route, v in recs:
            for pt in sorted(pts):
                d_ = dict(pt)
                e = arb(d_['e'])
                if not (e > 0):
                    continue
                try:
                    y = f(arb(d_['u']), arb(d_[sec]), e, br, arb(d_['l']) if 'l' in d_ else None)
                except (ValueError, ZeroDivisionError):
                    bad += 1
                    continue
                n += 1
                ex = float(y.upper()) - (0.0 if route == 'd' else v)
                worst = max(worst, ex)
                if not (ex <= (0.0 if route == 'd' else 1e-30)) or (route == 'd' and not (y < 0)):
                    bad += 1
    return dict(piece='sound ' + name, points=npts, evals=n, failed=bad, worst=worst)


def main(argv):
    global EPS0, DELTA0
    what = argv[1] if len(argv) > 1 else 'all'
    if '--eps0' in argv:
        EPS0 = float(argv[argv.index('--eps0') + 1])
    if '--delta0' in argv:
        DELTA0 = float(argv[argv.index('--delta0') + 1])
    out = argv[argv.index('--json') + 1] if '--json' in argv else None
    res = []
    log = lambda s: print(s, flush=True)
    names = ['validate', 'A-J', 'A-FC', 'B-J', 'B-F', 'outer', 'end']
    todo = names if what == 'all' else what.split(',')
    if any(w not in names for w in todo):                               # round 131: an unknown piece is an error, not a no-op
        raise SystemExit(f"unknown piece(s) {[w for w in todo if w not in names]}; choose from {names} or all")
    keeps = {}
    for w in todo:
        if w == 'validate':
            res.append(validate(log, eps0=EPS0))
        elif w in ('A-J', 'A-FC'):
            keeps[w] = []
            res.append(pieceA(w[2:], eps0=EPS0, delta0=DELTA0, log=log, keep=keeps[w]))
            res.append(soundness(w, fA, keeps[w]))
            log(json.dumps(res[-1]))
        elif w in ('B-J', 'B-F'):
            keeps[w] = []
            res.append(pieceB(w[2:], eps0=EPS0, delta0=DELTA0, log=log, keep=keeps[w]))
            res.append(soundness(w, fB, keeps[w]))
            log(json.dumps(res[-1]))
        elif w == 'outer':
            res.append(piece_outer(eps0=EPS0, delta0=DELTA0, log=log))
            res.append(handoff_check(res[-1], delta0=DELTA0))
            log(json.dumps(res[-1]))
        elif w == 'end':
            r = piece_end()
            log(json.dumps(r))
            res.append(r)
    if out:
        json.dump(dict(prec=PREC, eps0=EPS0, delta0=DELTA0, QA=QA, PB=PB, xb=float(XB), end_r=END_R, results=res),
                  open(out, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main(sys.argv)
