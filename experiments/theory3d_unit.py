"""
Round theory3, direction D: the exact scalar unit of a factorised readout of depth p, in closed form (library; no experiment).

A unit has the coefficient w = W(parameters) and is trained by the plain gradient flow of F = 1/2 <r^2> with the learning rate ell for both parameters, from u(0) = v(0) = s,
s^2 = eps*b ('mass' b > 0). Whatever the other units and the unit's own centre do, its DRIVE g(t) = <r phi(.; c(t))> enters only through the dual eta(t) = int_0^t g, and
(Theorem 1 of the note)

    depth 1   w = u,                u' = ell g                                  w = ell*eta                                    (no wall)
    depth 2   w = u^2 - v^2,        u' = 2 ell u g, v' = -2 ell v g             w = 2 eps b sinh(4 ell eta)                    (wall ln(1/(eps b))/(4 ell): logarithmic)
    depth p>=3 w = u^p - v^p,       u' = ell p u^(p-1) g, v' = -ell p v^(p-1) g w = s^p [(1-eta/lam)^(-q) - (1+eta/lam)^(-q)]    (wall lam = s^(-(p-2))/(ell p (p-2)), q = p/(p-2))

(the code of `src/core/factorized_rbf.py`: Linear, Product2, Cubic3, PowerN). The distance form (d = lam - |eta| > 0, x = d/lam):

    depth p>=3  |w| = (ell p (p-2) d)^(-q) [1 - (x/(2-x))^q]        depth 2  |w| = exp(-4 ell d) - (eps b)^2 exp(4 ell d)

so that, in the distance to the wall, the coefficient does not depend on eps except through the factor in brackets (exponentially small in the former case, polynomially in the latter).
Everything here is a function of floats; the ball-arithmetic checks are in `ball_checks` (python-flint).
"""
import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# depth p >= 3
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def q_of(p):
    """the exponent p/(p-2) of the dormant coefficient in the margin"""
    return p / (p - 2.0)


def lam_p(p, eps, b=1.0, ell=1.0):
    """the wall of the dual: s^(-(p-2))/(ell p (p-2)) with s^2 = eps b; the code's `PowerN.scale` (p = 3: `Cubic3.scale = 1/(3 ell sqrt(eps))`)"""
    return (eps * b) ** (-(p - 2) / 2.0) / (ell * p * (p - 2.0))


def W_p(eta, p, eps, b=1.0, ell=1.0):
    """w as a function of the dual, |eta| < lam: s^p [(1 - eta/lam)^(-q) - (1 + eta/lam)^(-q)] (odd, increasing, onto R); evaluated as s^p (1+sg)^(-q) expm1(2 q atanh(sg)) (no cancellation)"""
    lam, q, s = lam_p(p, eps, b, ell), q_of(p), np.sqrt(eps * b)
    sg = np.asarray(eta, float) / lam
    return s ** p * (1.0 + sg) ** (-q) * np.expm1(2.0 * q * np.arctanh(sg))


def Psi_p(eta, p, eps, b=1.0, ell=1.0):
    """the mass identity: int_0^eta W = s^p lam/(q-1) [(1 - sg)^(-(q-1)) + (1 + sg)^(-(q-1)) - 2], sg = eta/lam (q - 1 = 2/(p-2))"""
    lam, q, s = lam_p(p, eps, b, ell), q_of(p), np.sqrt(eps * b)
    sg = np.asarray(eta, float) / lam
    return s ** p * lam / (q - 1.0) * ((1.0 - sg) ** (-(q - 1.0)) + (1.0 + sg) ** (-(q - 1.0)) - 2.0)


def uv_p(eta, p, eps, b=1.0, ell=1.0):
    """the parameters u = s (1 - eta/lam)^(-1/(p-2)), v = s (1 + eta/lam)^(-1/(p-2))"""
    lam, s = lam_p(p, eps, b, ell), np.sqrt(eps * b)
    sg = np.asarray(eta, float) / lam
    return s * (1.0 - sg) ** (-1.0 / (p - 2.0)), s * (1.0 + sg) ** (-1.0 / (p - 2.0))


def w_of_distance_p(d, p, eps, b=1.0, ell=1.0):
    """|w| at the distance d = lam - |eta| from the wall: (ell p (p-2) d)^(-q) [1 - (x/(2-x))^q], x = d/lam in (0, 1]"""
    lam, q = lam_p(p, eps, b, ell), q_of(p)
    d = np.asarray(d, float)
    x = d / lam
    return (ell * p * (p - 2.0) * d) ** (-q) * (1.0 - (x / (2.0 - x)) ** q)


def dormant_bound_p(mu, p, eps, b=1.0):
    """sup |w| over |eta| <= (1 - mu) lam, exactly s^p [mu^(-q) - (2-mu)^(-q)] (attained at |eta| = (1 - mu) lam); at most (eps b)^(p/2) mu^(-p/(p-2))"""
    q, s = q_of(p), np.sqrt(eps * b)
    mu = np.asarray(mu, float)
    return s ** p * (mu ** (-q) - (2.0 - mu) ** (-q))


def dormant_bound_simple_p(mu, p, eps, b=1.0):
    """the form of the theorem: (eps b)^(p/2) mu^(-p/(p-2)); p = 3: (eps b)^(3/2) mu^(-3)"""
    return (eps * b) ** (p / 2.0) * np.asarray(mu, float) ** (-q_of(p))


def mass_bound_p(mu, p, eps, b=1.0, ell=1.0):
    """int_0^eta W d eta at |eta| = (1 - mu) lam, at most (eps b)/(2 ell p) mu^(-2/(p-2)) (the centre drift of a dormant unit is this over the smallest drive)"""
    return Psi_p((1.0 - np.asarray(mu, float)) * lam_p(p, eps, b, ell), p, eps, b, ell)


def d_theta_p(theta, p, eps, b=1.0, ell=1.0):
    """the distance from the wall at which |w| = theta (the unique root of w_of_distance_p(d) = theta); the coefficient is formed iff d < d_theta"""
    lam = lam_p(p, eps, b, ell)
    f = lambda ld: float(w_of_distance_p(np.exp(ld), p, eps, b, ell)) - theta
    return float(np.exp(brentq(f, np.log(lam) - 50.0, np.log(lam), xtol=1e-14, rtol=1e-14)))


def d_theta_p_leading(theta, p, ell=1.0):
    """the eps-free leading value theta^(-(p-2)/p)/(ell p (p-2)) (an upper bound of d_theta; d_theta is below it by a factor 1 - O(eps^(p/2)...))"""
    return theta ** (-(p - 2.0) / p) / (ell * p * (p - 2.0))


def mu_theta_p(theta, p, eps, b=1.0, ell=1.0):
    """the margin below which a dormant unit is formed: d_theta/lam, about theta^(-(p-2)/p) (eps b)^((p-2)/2): p = 3: theta^(-1/3) sqrt(eps b)"""
    return d_theta_p(theta, p, eps, b, ell) / lam_p(p, eps, b, ell)


def eta_of_w_p(w, p, eps, b=1.0, ell=1.0):
    """the inverse of W_p (monotone): the dual that carries the coefficient w (vectorised; sign kept)"""
    w = np.atleast_1d(np.asarray(w, float))
    lam = lam_p(p, eps, b, ell)
    out = np.zeros_like(w)
    for i, wi in enumerate(w):
        if wi == 0.0:
            continue
        a = abs(wi)
        if a > 1.0:                                   # near the wall: solve in the distance (no cancellation)
            ld = brentq(lambda t: float(w_of_distance_p(np.exp(t), p, eps, b, ell)) - a, np.log(lam) - 80.0, np.log(lam), xtol=1e-15, rtol=1e-15)
            out[i] = np.sign(wi) * (lam - np.exp(ld))
        else:                                         # solve in log(sigma), sigma = eta/lam
            t = brentq(lambda t: float(W_p(np.exp(t) * lam, p, eps, b, ell)) / a - 1.0, -745.0 + np.log(lam) * 0 - 1.0, np.log1p(-1e-15), xtol=1e-14, rtol=1e-14, maxiter=500)
            out[i] = np.sign(wi) * np.exp(t) * lam
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# depth 2 (the logarithmic scale)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def lam_2(eps, b=1.0, ell=1.0):
    """the wall ln(1/(eps b))/(4 ell) of the code's `Product2.scale` (the dual at which |w| = 1 up to (eps b)^2)"""
    return np.log(1.0 / (eps * b)) / (4.0 * ell)


def W_2(eta, eps, b=1.0, ell=1.0):
    """w = u^2 - v^2 = 2 eps b sinh(4 ell eta) (u = s exp(2 ell eta), v = s exp(-2 ell eta))"""
    return 2.0 * eps * b * np.sinh(4.0 * ell * np.asarray(eta, float))


def Psi_2(eta, eps, b=1.0, ell=1.0):
    """int_0^eta W_2 = eps b (cosh(4 ell eta) - 1)/(2 ell)"""
    return eps * b * (np.cosh(4.0 * ell * np.asarray(eta, float)) - 1.0) / (2.0 * ell)


def eta_of_w_2(w, eps, b=1.0, ell=1.0):
    return np.arcsinh(np.asarray(w, float) / (2.0 * eps * b)) / (4.0 * ell)


def w_of_distance_2(d, eps, b=1.0, ell=1.0):
    """|w| at the distance d from the wall lam_2 with e^(4 ell lam_2) = 1/(eps b): exp(-4 ell d) - (eps b)^2 exp(4 ell d)"""
    d = np.asarray(d, float)
    return np.exp(-4.0 * ell * d) - (eps * b) ** 2 * np.exp(4.0 * ell * d)


def dormant_bound_2(mu, eps, b=1.0):
    """sup |w| over |eta| <= (1 - mu) lam_2: (eps b)^mu - (eps b)^(2-mu) (exact, attained)"""
    e = eps * b
    return e ** np.asarray(mu, float) - e ** (2.0 - np.asarray(mu, float))


def d_theta_2(theta, eps, b=1.0, ell=1.0):
    """exactly: |w| >= theta iff d <= ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/(4 ell)"""
    e = eps * b
    return np.log(2.0 / (theta + np.sqrt(theta ** 2 + 4.0 * e ** 2))) / (4.0 * ell)


def mu_theta_2(theta, eps, b=1.0, ell=1.0):
    """the margin below which a dormant depth-2 unit is formed: d_theta/lam_2 = ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/ln(1/(eps b)); within 2 (eps b/theta)^2/ln(1/eps b)
    of the paper's rule ln(1/theta)/ln(1/eps)"""
    return d_theta_2(theta, eps, b, ell) / lam_2(eps, b, ell)


def rule_depth2(theta, eps):
    """the paper's rule ln(1/theta)/ln(1/eps)"""
    return np.log(1.0 / theta) / np.log(1.0 / eps)


def wall_shift_2(A, eps, b=1.0, ell=1.0):
    """the dual at which a depth-2 unit reaches the amplitude A, in units of lam_2: 1 + ln(A/b)/ln(1/eps) (exact up to (eps b/A)^2): the wall depends on the amplitude, at relative order
    1/ln(1/eps)"""
    return 1.0 + np.log(A / b) / np.log(1.0 / eps)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# depth 1
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def d1_fixed_centres(H, b0, ell, t):
    """depth 1, fixed centres: with the Gram H = Phi^T Phi/N and b0 = Phi^T y/N, w(t) = H^+ (I - exp(-ell H t)) b0 (w(0) = 0), drives g(t) = exp(-ell H t) b0; returns (w, g)"""
    ev, V = np.linalg.eigh(H)
    c = V.T @ b0
    ex = np.exp(-ell * np.outer(np.atleast_1d(t), ev))
    inv = np.where(ev > 1e-14 * ev.max(), 1.0 / np.maximum(ev, 1e-300), 0.0)
    gt = (ex * c[None, :]) @ V.T
    wt = ((1.0 - ex) * (inv * c)[None, :]) @ V.T
    return wt, gt


def d1_formation_certificate(g0, normH, ell, theta):
    """depth 1, fixed centres (Theorem 4 of the note): for a unit whose initial drive g0_j satisfies g0_j^2 > 2 theta ||H|| ||g0||, |w_j(t)| >= ell t |g0_j|/2 for t <= T_j =
    |g0_j|/(ell ||H|| ||g0||); the unit is formed (|w| > theta) by t_j = 2 theta/(ell |g0_j|) <= T_j. Returns (certified mask, t_j, T_j)."""
    g0 = np.asarray(g0, float)
    ng = float(np.linalg.norm(g0))
    T = np.abs(g0) / (ell * normH * ng)
    tj = 2.0 * theta / (ell * np.maximum(np.abs(g0), 1e-300))
    return (g0 ** 2 > 2.0 * theta * normH * ng), tj, T


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the paper's unit (mass x saturating amplitude): Lemma 2.1 of the T1 note, closed forms
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)


class PaperUnit:
    """w = (a/m) tanh x (eps + 2 v_s sinh^2 x), eta = (a/(2 G v_s)) (x + kap arctan(kap tanh x)), x = alpha/a, kap = sqrt((2 v_s - eps)/eps), wall lam = pi a kap/(4 G v_s),
    distance s(x) = lam - eta = (a kap/(2 G v_s)) (arctan(1/(kap tanh x)) - x/kap)."""

    def __init__(self, eps, law=LAW):
        self.m, self.G, self.gam, self.a = law['m'], law['G'], law['gamma'], law['a_max']
        self.vs = self.gam * self.a ** 2 / (4.0 * self.G)
        self.eps = float(eps)
        self.kap = np.sqrt((2.0 * self.vs - eps) / eps)
        self.lam = np.pi * self.a * self.kap / (4.0 * self.G * self.vs)

    def w_of_x(self, x):
        return (self.a / self.m) * np.tanh(x) * (self.eps + 2.0 * self.vs * np.sinh(x) ** 2)

    def eta_of_x(self, x):
        return self.a / (2.0 * self.G * self.vs) * (x + self.kap * np.arctan(self.kap * np.tanh(x)))

    def dist_of_x(self, x):
        th = np.tanh(x)
        return self.a * self.kap / (2.0 * self.G * self.vs) * (np.arctan(1.0 / (self.kap * np.maximum(th, 1e-300))) - x / self.kap)

    def x_of_w(self, w):
        w = abs(float(w))
        hi = 1.0
        while self.w_of_x(hi) < w:
            hi *= 2.0
        return brentq(lambda x: self.w_of_x(x) - w, 0.0, hi, xtol=1e-300, rtol=1e-15)

    def eta_of_w(self, w):
        return np.sign(w) * self.eta_of_x(self.x_of_w(w)) if w != 0 else 0.0

    def dist_of_w(self, w):
        """the distance to the wall carrying the coefficient w, computed from the distance form (no cancellation)"""
        return self.dist_of_x(self.x_of_w(w))

    def d_theta(self, theta):
        return self.dist_of_w(theta)

    def mu_theta(self, theta):
        return self.d_theta(theta) / self.lam

    def dormant_bound(self, mu):
        """sup |w| over |eta| <= (1-mu) lam: w at the x with dist = mu lam"""
        target = mu * self.lam
        hi = 1.0
        while self.dist_of_x(hi) > target:
            hi *= 2.0
        x = brentq(lambda x: self.dist_of_x(x) - target, 1e-12, hi, xtol=1e-300, rtol=1e-15)
        return float(self.w_of_x(x))

    def dormant_bound_T1(self, mu):
        """Lemma 2.2 of THEORY_T1: rho^(3/2)/(1 - s_eps^2) a/(m sqrt(2 v_s)) eps^(3/2) cos(th)/sin(th)^3 <= (...) mu^(-3), th = pi mu/2, s_eps = cot(th)/kap (needs s_eps < 1)"""
        th = np.pi * mu / 2.0
        s_eps = 1.0 / (np.tan(th) * self.kap)
        rho = 2.0 * self.vs / (2.0 * self.vs - self.eps)
        pre = rho ** 1.5 / (1.0 - s_eps ** 2) * self.a / (self.m * np.sqrt(2.0 * self.vs)) * self.eps ** 1.5
        return pre * np.cos(th) / np.sin(th) ** 3, pre * mu ** (-3)


def paper_B_const(law=LAW):
    """B in |w| ~ B d^(-3) at the distance d from the wall (eps -> 0): 4/(G gamma^2 m) = 33.33 (Corollary 2.3 of the T1 note)"""
    return 4.0 / (law['G'] * law['gamma'] ** 2 * law['m'])


def paper_c_u(law=LAW):
    """lam sqrt(eps) -> pi a_max sqrt(2 v_s)/(4 G v_s) = 7.0248 (the wall of the paper's unit, c_u eps^(-1/2))"""
    vs = law['gamma'] * law['a_max'] ** 2 / (4.0 * law['G'])
    return float(np.pi * law['a_max'] * np.sqrt(2.0 * vs) / (4.0 * law['G'] * vs))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the sliding unit: the first integral, the drift lead and the competitor's overtaking margin
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def lead_const(readout, ell=1.0, law=LAW):
    """a_ro in Psi*(eta) ~ a_ro eps x^(-2) at the margin x (x = d/lam, small): the mass identity near the wall. depth 3: 1/(6 ell); the paper's unit: 2/(G gamma^2 m c_u^2) = 0.3377"""
    if readout == 'D3':
        return 1.0 / (6.0 * ell)
    if readout == 'P0':
        return 2.0 / (law['G'] * law['gamma'] ** 2 * law['m'] * paper_c_u(law) ** 2)
    raise ValueError(readout)


def drift_lead(kappa, g0, eps, x, x_a, a_ro):
    """the lead of a dormant unit's dual (in units of its wall) accrued while its margin goes from x_a to x, at the first order in kappa = mu_c h^2 (mu_c the centres' rate, h = <r d_c phi>), the
    landscape g = g0 + h dc being the same since the start of the run: kappa a_ro eps/g0^2 (1/x - 1/x_a), a_ro the mass-identity constant (Psi*(eta) ~ a_ro eps/x^2 near the wall).
    Derivation: first integral g^2 = g0^2 + 2 kappa Psi*(eta), so g - g0 ~ kappa Psi*/g0 and the dual gains int (g - g0) dt = kappa/g0^2 int Psi* d eta."""
    return kappa * a_ro * eps / g0 ** 2 * (1.0 / x - 1.0 / x_a)


def lead_first_order_D3(kappa, g0, eps, sg, ell=1.0):
    """the exact first-order (in kappa) lead of the depth-3 unit's dual from the start (eta = 0, c = 0), in units of the wall, when the frozen dual is at sg = eta_0/lam: kappa eps/(6 ell g0^2)
    [1/(1-sg) - 1/(1+sg) - 2 sg] (int_0^sg Psi* d sg', Psi* from Psi_p)"""
    return kappa * eps / (6.0 * ell * g0 ** 2) * (1.0 / (1.0 - sg) - 1.0 / (1.0 + sg) - 2.0 * sg)


def overtake_margin(Kc, x_a):
    """the margin x at which a competitor overtakes its entrant: the root of x = Kc (1/x - 1/x_a) (Kc = kappa a_ro eps/g0^2): x^2 = Kc (1 - x/x_a); about sqrt(Kc) = sqrt(kappa a_ro eps)/|g0|"""
    b = Kc / x_a
    return (-b + np.sqrt(b * b + 4.0 * Kc)) / 2.0


def integrate_sliding_unit(G, dG, eta0, c0, mu_c, W, T, lam, c_stop=None, rtol=1e-10, atol=1e-13, wall_frac=0.9999):
    """the frozen-environment unit: eta' = G(c), c' = mu_c W(eta) G'(c) up to time T, stopped when eta reaches wall_frac*lam or, if c_stop = (c*, r), when |c - c*| <= r;
    returns (t, eta, c) at the solver's steps"""
    from scipy.integrate import solve_ivp
    f = lambda t, y: [G(y[1]), mu_c * float(W(y[0])) * dG(y[1])]
    ev = [lambda t, y: wall_frac * lam - y[0]]
    ev[0].terminal = True
    if c_stop is not None:
        ev.append(lambda t, y: abs(y[1] - c_stop[0]) - c_stop[1])
        ev[1].terminal = True
    sol = solve_ivp(f, (0.0, T), [eta0, c0], method='DOP853', rtol=rtol, atol=atol, events=ev)
    return sol.t, sol.y[0], sol.y[1]


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ball arithmetic
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def ball_checks(ps=(3, 4, 5, 6), epss=('1e-4', '1e-6', '1e-8', '1e-12'), bs=('1', '3/2'), ells=('1', '1/4'), mus=("1/2", "1/10", "1/100", "1/1000"), thetas=("1/100", "1/1000"), prec=400):
    """Interval checks (python-flint, precision `prec` bits) of the closed-form statements of the note over the stated grids (exact rationals as inputs):
      (B1) |W_p((1-mu) lam)| = s^p [mu^-q - (2-mu)^-q]  (the exact dormant supremum) and it is <= s^p mu^-q;
      (B2) the distance form W_p(lam - d) = (ell p (p-2) d)^-q [1 - (x/(2-x))^q];
      (B3) the derivative identities d u/d eta = ell p u^(p-1), d v/d eta = -ell p v^(p-1) at eta = (1-mu) lam (arb_series), hence u' = ell p u^(p-1) g, v' = -ell p v^(p-1) g;
      (B4) the mass identity int_0^eta W_p = Psi_p(eta) (rigorous quadrature, acb.integral) and Psi_p((1-mu) lam) <= eps b mu^(-2/(p-2))/(2 ell p);
      (B5) the formation distance: d_lo < d_theta < d_hi with W_p(lam - d_lo) > theta > W_p(lam - d_hi), d_hi = d^0 = theta^(-1/q)/(ell p (p-2)), d_lo = d^0 (1 - (x0/(2-x0))^q)^(1/q);
      (B6) depth 2: 2 eps b sinh(4 ell lam_2(1-mu)) = (eps b)^mu - (eps b)^(2-mu) for lam_2 with e^(4 ell lam_2) = 1/(eps b); the exact distance form; the bracket
           mu_0 - delta <= mu_theta <= mu_0 for mu_0 = ln(1/theta)/ln(1/(eps b)), delta = -ln(1 - (eps b)^(2-2 mu_0))/ln(1/(eps b)) (monotone f(mu) = e^mu - e^(2-mu)).
    Returns dict(name -> (n checked, n failed))."""
    from flint import arb, ctx, fmpq, acb, arb_series
    ctx.prec = prec
    res = {k: [0, 0] for k in ('B1', 'B2', 'B3', 'B4', 'B5', 'B6')}

    def Q(sx):
        return arb(fmpq(*[int(t) for t in sx.split('/')])) if '/' in sx else (arb(sx) if 'e' in sx else arb(int(sx)))

    def ok(k, cond):
        res[k][0] += 1
        res[k][1] += 0 if cond else 1

    for p in ps:
        qa = arb(p) / (p - 2)
        for es in epss:
            for bs_ in bs:
                eb = Q(es) * Q(bs_)
                for ls in ells:
                    ell = Q(ls)
                    s = eb.sqrt()
                    lam = s ** (-(p - 2)) / (ell * p * (p - 2))
                    for ms in mus:
                        mu = Q(ms)
                        eta = (1 - mu) * lam
                        sg = eta / lam
                        Wv = s ** p * ((1 - sg) ** (-qa) - (1 + sg) ** (-qa))
                        exact = s ** p * (mu ** (-qa) - (2 - mu) ** (-qa))
                        simple = s ** p * mu ** (-qa)
                        ok('B1', (Wv - exact).abs_upper() < arb('1e-40') * exact.abs_lower() and (simple - exact).lower() > 0)
                        d = mu * lam
                        x = d / lam
                        dform = (ell * p * (p - 2) * d) ** (-qa) * (1 - (x / (2 - x)) ** qa)
                        ok('B2', (Wv - dform).abs_upper() < arb('1e-40') * exact.abs_lower())
                        e0 = arb_series([eta, 1], prec=3)
                        us = s * (1 - e0 / lam) ** (-arb(1) / (p - 2))
                        vs = s * (1 + e0 / lam) ** (-arb(1) / (p - 2))
                        u0, v0 = us[0], vs[0]
                        du, dv = us[1], vs[1]
                        ok('B3', (du - ell * p * u0 ** (p - 1)).abs_upper() < arb('1e-40') * abs(du).abs_lower() and (dv + ell * p * v0 ** (p - 1)).abs_upper() < arb('1e-40') * abs(dv).abs_lower())
                        if es in ('1e-4', '1e-8') and ms in ('1/2', '1/100') and bs_ == '1':
                            ps_ = s ** p * lam / (qa - 1) * ((1 - sg) ** (-(qa - 1)) + (1 + sg) ** (-(qa - 1)) - 2)
                            qq = arb(qa)

                            def f(z, _):
                                zz = z / lam
                                return s ** p * ((1 - zz) ** (-qq) - (1 + zz) ** (-qq))
                            I = acb.integral(lambda z, _: f(z, _), 0, acb(eta))
                            ok('B4', (I.real - ps_).abs_upper() < arb('1e-30') * ps_ and (eb / (2 * ell * p) * mu ** (-arb(2) / (p - 2)) - ps_).lower() > 0)
                        else:
                            ps_ = s ** p * lam / (qa - 1) * ((1 - sg) ** (-(qa - 1)) + (1 + sg) ** (-(qa - 1)) - 2)
                            ok('B4', (eb / (2 * ell * p) * mu ** (-arb(2) / (p - 2)) - ps_).lower() > 0)
                    for ts in thetas:
                        th = Q(ts)
                        d0 = th ** (-1 / qa) / (ell * p * (p - 2))
                        x0 = d0 / lam
                        d_lo = d0 * (1 - (x0 / (2 - x0)) ** qa) ** (1 / qa)

                        def Wd(dd):
                            xx = dd / lam
                            return (ell * p * (p - 2) * dd) ** (-qa) * (1 - (xx / (2 - xx)) ** qa)
                        ok('B5', (Wd(d_lo) - th).lower() >= 0 and (th - Wd(d0)).lower() > 0 and d_lo < d0)
    for es in epss:
        for bs_ in bs:
            eb = Q(es) * Q(bs_)
            for ls in ells:
                ell = Q(ls)
                lam2 = (1 / eb).log() / (4 * ell)
                for ms in mus:
                    mu = Q(ms)
                    val = 2 * eb * ((4 * ell) * lam2 * (1 - mu)).sinh()
                    ok('B6', (val - (eb ** mu - eb ** (2 - mu))).abs_upper() < arb('1e-40') * (eb ** mu).abs_lower())
                    d = mu * lam2
                    ok('B6', (val - ((-4 * ell * d).exp() - eb ** 2 * (4 * ell * d).exp())).abs_upper() < arb('1e-40') * (eb ** mu).abs_lower())
                for ts in thetas:
                    th = Q(ts)
                    mu0 = (1 / th).log() / (1 / eb).log()
                    delta = -((1 - eb ** (2 - 2 * mu0)).log()) / (1 / eb).log()
                    f = lambda m: eb ** m - eb ** (2 - m)
                    ok('B6', (f(mu0) - th).upper() < 0 and (f(mu0 - delta) - th).lower() >= 0 and mu0 < 1)
                    # the exact formation distance: w(d_theta) = theta
                    dth = (2 / (th + (th ** 2 + 4 * eb ** 2).sqrt())).log() / (4 * ell)
                    ok('B6', ((-4 * ell * dth).exp() - eb ** 2 * (4 * ell * dth).exp() - th).abs_upper() < arb('1e-40'))
    return {k: tuple(v) for k, v in res.items()}
