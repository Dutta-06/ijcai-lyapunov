"""
Two ML-facing statements about the escape from the dormant state: normalised optimizers (Adam, sign descent) erase the ordering, and a pure
small start in two layers does not leave its saddle
======================================================================================================================================

F1.  Adam and sign descent.  A unit w = W(p) of the factorised readouts of `factorized_rbf` (D3: w = u³ − v³ from (s, s); the paper's unit has
the same escape phase) is dormant while its drive g_i = ⟨r φ_i⟩ is a number that does not change sign.  Plain gradient descent moves it at the
speed lr·g_i·∂W/∂p (the escape takes about 1/(3 lr g_i s): ordered by g_i, unbounded as s → 0).  A per-coordinate normalised update moves
every parameter by lr·R_t, R_t = m̂_t/(√v̂_t + ε_adam), a number that does not depend on the SCALE of the gradient (ε_adam → 0): the escape takes
(u_θ − s)/lr steps for EVERY unit of nonzero drive.  Statements (proved in the note, checked here):
  F1a  sign descent: all units whose drive keeps its sign up to the threshold cross it at the same step (exactly, for equal starts), whatever
       the ratio of their drives; the escape time does not depend on the start s once lr > s.
  F1b  Adam (ε_adam → 0): a frozen drive gives the same trajectory for every drive (exact scale invariance); in general the escape steps of
       any two units with positive nondecreasing gradients are within the factor C/c, c = 1 − β1, C = sup_t C_t(β1, β2) (7.27 for the defaults).
  F1c  the coupled pair and duplicates: all units of positive drive rise together; a unit that is not in the sliding path reaches the
       coefficient ρ/(1 + ρ) (D3, sign descent) before its drive turns; duplicate atoms share the amplitude (ratio of the coefficients → 1).
  F1d  what keeps the mechanism: ε_adam large against the gradients up to the decision (Adam is then SGD with rate lr/ε_adam), a SHARED
       second moment, normalised gradient descent with the global norm (GD on a clock).
F2.  Two layers.  The balanced product of D small factors: z' = lr z^(D−1)(g − z^D), whose time from s to z_θ is the quadrature
       T = (1/(lr g))·[(s^(−(D−2)) − z_θ^(−(D−2)))/(D−2) + ∫_s^{z_θ} z dz/(g − z^D)]  (D ≥ 3),  T = (1/(lr g))·[ln(z_θ/s) + ∫ z dz/(g − z²)]  (D = 2):
a logarithmic escape for the diagonal two-layer network y = Σ u_j v_j x_j, the power s^(−(D−2)) for D ≥ 3.  The paper's start (w = 0 in both
layers) is an exact equilibrium of the two-layer product.  The hybrid (first layer factorised and small, second layer O(1)) obeys the exact
conservation law b² − κ²Σ_l u_l² = const and reduces to the one-layer escape with the effective drive b·g.

Usage: python experiments/theory3f_escape.py [--run [--parts A,B,C,D,D2,E,F,G,H,K]] [--from-raw] → results/synthetic/theory3f.md (+ _raw.npz).
Environment: T3F_OUT (output directory).  Deterministic: every draw is seeded; the report is regenerated from the raw alone.
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
from scipy.integrate import quad, solve_ivp

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.environ.get('T3F_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3f_raw.npz'), os.path.join(OUT, 'theory3f.md')

# ---- registries ------------------------------------------------------------------------------------------------------------------
BETAS = (0.9, 0.999)                                   # Adam's defaults
EA_TINY = 1e-30                                         # "ε_adam → 0"
# A: sign descent / Adam on D3 units with orthonormal atoms (drives g_i), threshold on the coefficient
G_A = (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
S_A, LR_A, THW_A, NSTEP_A = (1.3e-3, 1e-5, 1e-7), 1e-4, 1e-6, 400
H_GD_A, S_GD_A = 0.01, 0.05                             # the plain-GD check (explicit Euler, step h) against the closed form
# B: Adam on one coordinate with a frozen drive
SCALES_B, S_B, LR_B, UC_B = tuple(10.0 ** k for k in range(-8, 9, 4)), 1e-3, 1e-4, 0.1
BETAS_B = ((0.9, 0.999), (0.9, 0.99), (0.9, 0.95), (0.5, 0.999), (0.0, 0.999), (0.99, 0.9999))
EA_B = (EA_TINY, 1e-8, 1e-3, 1.0, 1e3)
# C: the coupled pair (Gram [[1, ρ], [ρ, 1]], target = atom 1) and duplicate atoms
RHO_C, S_C, THW_C = (0.2, 0.5, 0.8), (1.31e-2, 1.31e-4), 1e-3
OPTS_C = (('gd', dict(lr=0.05)), ('sign', dict(lr=1e-3)), ('adam', dict(lr=1e-3, eps_adam=EA_TINY)), ('adam_g', dict(lr=1e-3)), ('ngd', dict(lr=1e-3)))
S_DUP, DUP_RATIO = 1.31e-2, 0.7                         # duplicates: s_1 = 0.0131, s_2 = 0.7 s_1 (starts that are not multiples of lr: see the note)
# D: the scan of ε_adam on the pair
EA_D, S_D, RHO_D = (1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 1.0, 10.0), (1.31e-2, 1.31e-3), 0.5
# E, F: the repo's design (21 kernels, exp_parametrization), Adam / sign, then the modifications
EPS_E, STEPS_E, LRS_E = 1e-6, 20000, (1e-3, 1e-4)
STEPS_F, H_F, EPS_F = 200000, 0.05, 1e-6
EA_F = (1e-8, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 1e-1, 1.0, 10.0)
EA_H, STEPS_H = (1e-2, 1.0), 3000000                      # H: the paper's unit with a large eps_adam needs the long horizon of its wall (7000 / drive)
RHO_D2, S_D2, GRID_D2 = (0.2, 0.5, 0.8), (1.31e-2, 1.31e-3), tuple(np.geomspace(0.02, 8.0, 12))     # D2: eps_adam in units of s
# G: two layers
D_G, S_G = (2, 3, 4, 5, 6), (0.4, 0.3, 0.2, 0.1, 0.05)
S_DIAG = (1e-1, 1e-2, 1e-4, 1e-8)
S_PAIR, A0_PAIR = (0.5, 0.4, 0.3, 0.25), (1e-3, 1e-6)
S_H, B0_H, U_H = (0.1, 0.05, 0.02), (1.0, 0.3, 0.1, None), 0.2     # hybrid: b0 = None: the balanced pure start b0 = κ s
N_X, N_L, SIG_X, KH_H = 81, 21, 0.12, 2
EPS_KAN, STEPS_KAN = 1e-3, 3000
# K: the additive KAN [5, 1] of the spike (numpy, gradient descent against Adam), the development targets
TARGETS_K, NOISES_K, SEED_K, EPS_K, M_K, STEPS_K = (101, 102, 103), (0.1, 0.4), 0, 1e-3, 10, 12000
OPTS_K = (('GD (rate 1)', 'gd', 0.0), ('Adam (lr 0.01)', 'adam', 1e-8), ('Adam, eps_adam = 1 (GD rates)', 'adam_r', 1.0), ('Adam, eps_adam = 0.01 (GD rates)', 'adam_r', 0.01),
          ('Adam, shared v (lr 0.003)', 'adam_g', 1e-8))


# ===============================================================================================================================
# optimizers
# ===============================================================================================================================
class Opt:
    """One optimizer on the vector of ALL parameters (the callers pass −gradient).  kinds:
      'gd'      p += lr·pre·gneg                                      (plain gradient descent; `pre`: the group learning rates, default 1)
      'sign'    p += lr·sign(gneg)                                    (Adam with β1 = β2 = 0, ε_adam → 0)
      'adam'    p += lr·m̂/(√v̂ + ε_adam)                              (per-coordinate second moment, bias-corrected; the PyTorch rule)
      'adam_g'  Adam with ONE scalar second moment, v = EMA of mean((pre·gneg)²): keeps the relative gradients of the coordinates
      'ngd'     p += lr·(pre·gneg)/‖pre·gneg‖                         (normalised gradient descent, global norm: gradient flow on a clock)
    `lr` may be an array (per-parameter step sizes).  For 'adam_g' and 'ngd' the gradient is first multiplied by `pre`."""

    def __init__(self, kind, lr, beta1=0.9, beta2=0.999, eps_adam=1e-8, pre=None):
        assert kind in ('gd', 'sign', 'adam', 'adam_g', 'ngd')
        self.kind, self.lr, self.b1, self.b2, self.ea, self.pre = kind, lr, float(beta1), float(beta2), float(eps_adam), pre
        self.m = self.v = None
        self.t = 0

    def step(self, p, gneg):
        self.t += 1
        k = self.kind
        if self.pre is not None and k in ('gd', 'adam_g', 'ngd'):
            gneg = self.pre * gneg
        if k == 'gd':
            return p + self.lr * gneg
        if k == 'sign':
            return p + self.lr * np.sign(gneg)
        if k == 'ngd':
            n = np.linalg.norm(gneg)
            return p if n == 0.0 else p + self.lr * gneg / n
        if self.m is None:
            self.m, self.v = np.zeros_like(p), np.zeros_like(p)
        self.m = self.b1 * self.m + (1 - self.b1) * gneg
        gg = gneg ** 2 if k == 'adam' else np.mean(gneg ** 2) * np.ones_like(gneg)
        self.v = self.b2 * self.v + (1 - self.b2) * gg
        mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)
        den = np.sqrt(vh) + self.ea
        return p + self.lr * np.divide(mh, den, out=np.zeros_like(den), where=den > 0)


def adam_ratio_bounds(b1=0.9, b2=0.999, tmax=10 ** 6):
    """The bounds on R_t = m̂_t/√v̂_t (ε_adam = 0) used by the window theorem.
    Upper (Cauchy–Schwarz, any real gradient sequence): R_t² ≤ Σ_k a_k²/b_k with a, b the bias-corrected EMA weights, i.e.
        C_t² = (1−b1)²(1−b2^t)/((1−b2)(1−b1^t)²) · Σ_{j<t}(b1²/b2)^j ,   C = sup_t C_t.
    Lower (positive NONDECREASING gradients): m̂_t ≥ (1−b1)/(1−b1^t)·G_t and √v̂_t ≤ G_t, so R_t ≥ 1 − b1.   Returns (c, C, t at the sup)."""
    gam = b1 * b1 / b2
    t = np.unique(np.r_[np.arange(1, 4000), np.geomspace(4000, tmax, 400).astype(np.int64)]).astype(float)
    geo = (1 - gam ** t) / (1 - gam) if gam != 1 else t
    C2 = (1 - b1) ** 2 * (1 - b2 ** t) / ((1 - b2) * (1 - b1 ** t) ** 2) * geo
    k = int(np.argmax(C2))
    return 1.0 - b1, float(np.sqrt(C2[k])), int(t[k])


# ===============================================================================================================================
# F1: units, the diagonal model, the pair, duplicates
# ===============================================================================================================================
def d3_w(P):
    """P (d, 2) = (u, v): w = u³ − v³."""
    return P[:, 0] ** 3 - P[:, 1] ** 3


def d3_dir(P, g):
    """−dL/dP of the D3 units with drives g (d,): (3u²g, −3v²g)."""
    return np.stack([3 * P[:, 0] ** 2 * g, -3 * P[:, 1] ** 2 * g], 1)


def run_units(opt, G, c, s, n_steps, thr):
    """d D3 units (w_i = u_i³ − v_i³, (u, v)(0) = (s_i, s_i)) on F = ½(w − w*)ᵀG(w − w*), c = G w* (drive of unit i: g_i = c_i − (Gw)_i), by `opt`,
    n_steps steps.  s: scalar or (d,).  Returns W (n_steps + 1, d) and T (d,): the first step with |w_i| ≥ thr (nan: never)."""
    d = len(c)
    P = np.stack([np.broadcast_to(np.asarray(s, float), (d,))] * 2, 1).copy()
    W = np.empty((n_steps + 1, d))
    W[0] = d3_w(P)
    T = np.full(d, np.nan)
    with np.errstate(all='ignore'):
        for t in range(1, n_steps + 1):
            g = c - G @ d3_w(P)
            P = opt.step(P.ravel(), d3_dir(P, g).ravel()).reshape(d, 2)
            W[t] = d3_w(P)
            new = np.isnan(T) & (np.abs(W[t]) >= thr)
            T[new] = t
            if not np.all(np.isfinite(P)):
                W[t:] = np.nan
                break
    return W, T


def sign_escape_step(s, thw, lr):
    """The step at which a D3 unit moved by sign descent from (s, s) reaches |w| = thw: u_t = s + lr t, v_t = s − lr t, so w_t = 2(3 s² lr t + (lr t)³)."""
    t = 1
    while 2.0 * (3 * s * s * lr * t + (lr * t) ** 3) < thw:
        t += 1
    return t


def gd_escape_time(g, s, u_c):
    """Plain gradient FLOW of the one-sided unit w = u³, frozen drive g, unit rate: u' = 3 g u² ⇒ time to u_c from s: (1/s − 1/u_c)/(3g)."""
    return (1.0 / s - 1.0 / u_c) / (3.0 * g)


def frozen_adam_run(s, u_c, lr, beta1, beta2, ea, scale, max_steps=10 ** 6):
    """Adam on ONE coordinate with the gradient scale·3u² (a frozen drive `scale`): the steps to reach u = u_c from u = s, and the extreme values
    of R_t = m̂_t/(√v̂_t + ε_adam) on the way.  With ε_adam = 0 the steps do not depend on `scale` (invariance under a positive rescaling)."""
    u, m, v = float(s), 0.0, 0.0
    rmin, rmax = np.inf, 0.0
    for t in range(1, max_steps + 1):
        gr = scale * 3.0 * u * u
        m = beta1 * m + (1 - beta1) * gr
        v = beta2 * v + (1 - beta2) * gr * gr
        r = (m / (1 - beta1 ** t)) / (np.sqrt(v / (1 - beta2 ** t)) + ea)
        rmin, rmax = min(rmin, r), max(rmax, r)
        u += lr * r
        if u >= u_c:
            return t, rmin, rmax
    return -1, rmin, rmax


def pair_gram(rho):
    return np.array([[1.0, rho], [rho, 1.0]]), np.array([1.0, rho])      # G and c = G w*, w* = (1, 0): the target is atom 1


def run_pair(kind, kw, rho, s, n_steps, thw):
    """The pair (atom 1 = the target, atom 2 at overlap ρ) from s: returns (T1, T2, peak |w_2|, final w)."""
    G, c = pair_gram(rho)
    W, T = run_units(Opt(kind, **kw), G, c, s, n_steps, thw)
    return T[0], T[1], float(np.nanmax(np.abs(W[:, 1]))), W[-1]


def run_duplicates(kind, kw, s1, ratio, n_steps):
    """Duplicate atoms (φ_1 = φ_2, G = [[1, 1], [1, 1]], target amplitude 1) with the starts s_1, s_2 = ratio·s_1: the final coefficients."""
    G, c = np.ones((2, 2)), np.ones(2)
    W, T = run_units(Opt(kind, **kw), G, c, np.array([s1, ratio * s1]), n_steps, 1e-3)
    return W[-1]


# ===============================================================================================================================
# the network on the repo's design (E, F)
# ===============================================================================================================================
def train_opt(FR, readout, K, fstar, c0, eps, n_steps, opt, n_t=2001):
    """The repo's Adam loop (`factorized_rbf.train_adam`) with an arbitrary `Opt` acting on z = (parameters, centres) with the gradient of F in the raw
    coordinates (unit rates); recorded on n_t steps, time = step/n_steps.  `opt.kind == 'adam'` with the same betas and ε_adam reproduces train_adam."""
    c0 = FR.pts(c0)
    d, n = c0.shape
    p0 = readout.init(d, eps, None, 0.0)
    z = np.r_[p0.ravel(), c0.ravel()]
    base = FR.make_rhs(FR.Unit(readout, np.ones_like(readout.lr)), K, fstar, 1.0, d, n)
    rec = np.unique(np.linspace(0, n_steps, n_t).round().astype(np.int64))
    Z = np.empty((rec.size, z.size))
    Z[0] = z
    k_rec, status = 1, 'ok'
    with np.errstate(all='ignore'):
        for step in range(1, n_steps + 1):
            z = opt.step(z, base(0.0, z))
            if step == rec[k_rec]:
                Z[k_rec] = z
                k_rec += 1
                if not np.all(np.isfinite(z)) or np.abs(z).max() > 1e9:
                    status = 'diverged'
                    break
            if k_rec >= rec.size:
                break
    if k_rec < rec.size:
        Z[k_rec:] = Z[k_rec - 1]
    W, C = FR._unpack(readout, Z, d, n)
    taus = rec / n_steps
    return dict(taus=taus, W=W, C=C, status=status, nfev=n_steps, secs=0.0, scale=float(n_steps), eps=eps, readout=readout.name)


def formation_steps(FR, res, thr, n_steps):
    """The step at which each unit first has |w| ≥ thr (nan: never)."""
    return FR.formation_times(res['taus'], res['W'], thr) * n_steps


def design():
    """The repo's design (exp_parametrization): the dictionary, the target and its drives at the start."""
    import exp_parametrization as P
    from src.core import factorized_rbf as FR
    warnings.filterwarnings('ignore')
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    return P, FR, L, K, c0, cs, As, f


# ===============================================================================================================================
# F2: chains, the diagonal two-layer network, the paper-unit pair, the hybrid, the KAN toy
# ===============================================================================================================================
def chain_time_quad(D, g, s, z_t):
    """The time of the balanced product of D factors (rate 1 each) from s to z_t: ∫_s^{z_t} dz/(z^{D−1}(g − z^D))."""
    return quad(lambda z: 1.0 / (z ** (D - 1) * (g - z ** D)), s, z_t, epsabs=0, epsrel=1e-12, limit=400, points=None)[0]


def chain_time_leading(D, g, s):
    """Leading term: s^{−(D−2)}/((D−2) g) for D ≥ 3, ln(1/s)/g for D = 2."""
    return np.log(1.0 / s) / g if D == 2 else s ** (-(D - 2)) / ((D - 2) * g)


def chain_time_ode(D, g, s, z_t):
    ev = lambda t, z: z[0] - z_t
    ev.terminal = True
    sol = solve_ivp(lambda t, z: [z[0] ** (D - 1) * (g - z[0] ** D)], [0, 1e12], [s], events=ev, rtol=1e-11, atol=1e-16, method='LSODA')
    return float(sol.t_events[0][0])


def diag_two_layer_w(t, g, w0):
    """y = uvx, u = v balanced, w = uv, target g, rate 1 for u and v: w' = 2w(g − w), w(t) = g/(1 + (g/w0 − 1)e^{−2gt})."""
    return g / (1.0 + (g / w0 - 1.0) * np.exp(-2.0 * g * np.asarray(t, float)))


def diag_two_layer_time(g, w0, frac=0.5):
    """The time to frac·g: ln((g/w0 − 1)·frac/(1 − frac))/(2g)."""
    return float(np.log((g / w0 - 1.0) * frac / (1.0 - frac)) / (2.0 * g))


def diag_gd_check(g, s, h, t_end):
    """Explicit Euler (step h) of u' = v(g − uv), v' = u(g − uv) from (s, s): w(t_end)."""
    u = v = s
    for _ in range(int(round(t_end / h))):
        r = g - u * v
        u, v = u + h * v * r, v + h * u * r
    return u * v


def pair_rhs(t, z, cstar):
    """Two layers of the paper's unit (unsaturated, m = 1), rates 1: w_k = ν_k²α_k, c = w_1w_2, F = ½(c − c*)²; z = (α_1, ν_1, α_2, ν_2)."""
    a1, n1, a2, n2 = z
    w1, w2 = n1 ** 2 * a1, n2 ** 2 * a2
    G = cstar - w1 * w2
    return [n1 ** 2 * w2 * G, 2 * n1 * a1 * w2 * G, n2 ** 2 * w1 * G, 2 * n2 * a2 * w1 * G]


def pair_time(s, a0, cstar=1.0):
    """The time at which c = w_1w_2 reaches c*/2 from (α_1, ν_1, α_2, ν_2) = (0, s, a0, s)."""
    ev = lambda t, z, cs: z[0] * z[1] ** 2 * z[2] * z[3] ** 2 - 0.5 * cs
    ev.terminal = True
    sol = solve_ivp(pair_rhs, [0, 1e12], [0.0, s, a0, s], args=(cstar,), events=ev, rtol=1e-11, atol=1e-22, method='LSODA')
    return float(sol.t_events[0][0]) if len(sol.t_events[0]) else float('nan')


def pair_time_sym(s, a0, cstar=1.0):
    """The symmetric start α_1(0) = α_2(0) = a0 (both layers alike): α' = (s² + 2α²)²α(c* − (s² + 2α²)²α²) (ν² − 2α² is conserved), the time by quadrature
    to c = c*/2."""
    f = lambda a: 1.0 / ((s * s + 2 * a * a) ** 2 * a * (cstar - (s * s + 2 * a * a) ** 2 * a * a))
    from scipy.optimize import brentq
    af = brentq(lambda a: (s * s + 2 * a * a) ** 2 * a * a - 0.5 * cstar, a0, 10.0)
    return float(quad(f, a0, af, epsabs=0, epsrel=1e-11, limit=500)[0])


def rbf_dictionary():
    xs = np.linspace(-1, 1, N_X)
    ctr = np.linspace(-1, 1, N_L)
    Phi = np.exp(-(xs[:, None] - ctr[None]) ** 2 / (2 * SIG_X ** 2))
    fstar = 0.9 * np.exp(-(xs - 0.3) ** 2 / (2 * SIG_X ** 2))
    return xs, Phi, fstar


def hybrid_run(s, b0, U, Kh=KH_H):
    """f = Σ_k b_k Σ_l u_kl³ φ_l (first layer: one-sided cubic units from s; second layer b_k(0) = b0, rates 1 for u and for b), a feature at 0.3 (between kernels),
    gradient flow.  Stops when max u = U.  Returns the escape time, the conservation residual max|b² − κ²Σ_l(u_l² − s²) − b0²| (κ² = 1/3), the relative
    variation of the leading drive up to the stop, and b at the stop."""
    xs, Phi, fstar = rbf_dictionary()
    N, Lk = Phi.shape
    kappa2 = 1.0 / 3.0

    def rhs(t, z):
        u, b = z[:Kh * Lk].reshape(Kh, Lk), z[Kh * Lk:]
        w = u ** 3
        r = fstar - b @ (w @ Phi.T)
        gam = Phi.T @ r / N
        return np.r_[(3 * u ** 2 * b[:, None] * gam[None, :]).ravel(), w @ gam]

    z0 = np.r_[np.full(Kh * Lk, s), np.full(Kh, b0)]
    ev = lambda t, z: z[:Kh * Lk].max() - U
    ev.terminal = True
    sol = solve_ivp(rhs, [0, 1e9], z0, events=ev, rtol=1e-10, atol=1e-16, method='LSODA', dense_output=False)
    z = sol.y[:, -1]
    u, b = z[:Kh * Lk].reshape(Kh, Lk), z[Kh * Lk:]
    cons = float(np.abs(b ** 2 - kappa2 * (u ** 2 - s ** 2).sum(1) - b0 ** 2).max())
    lead = int(np.argmax((Phi.T @ fstar / N)))
    drive = []
    for k in range(sol.y.shape[1]):
        zz = sol.y[:, k]
        uu, bb = zz[:Kh * Lk].reshape(Kh, Lk), zz[Kh * Lk:]
        r = fstar - bb @ ((uu ** 3) @ Phi.T)
        drive.append((Phi.T @ r / N)[lead])
    drive = np.array(drive)
    delta = float(np.max(np.abs(drive / drive[0] - 1.0)))
    return float(sol.t_events[0][0]), cons, delta, b.copy(), float((Phi.T @ fstar / N).max()), int((Phi.T @ fstar / N > 0).sum())


def hybrid_bounds(s, b0, U, g0, delta, n_pos):
    """The theorem's bounds on the escape time of the leading unit (drive within (1 ± δ) of g0, positive-drive units of the node n_pos):
        lower  ∫_s^U du/(3(1+δ)g0 u² b_max(u)),  b_max(u)² = b0² + κ² n_pos (u² − s²);   upper  (1/s − 1/U)/(3(1−δ)g0 b0)."""
    kappa2 = 1.0 / 3.0
    lo = quad(lambda u: 1.0 / (3 * (1 + delta) * g0 * u ** 2 * np.sqrt(b0 ** 2 + kappa2 * n_pos * (u ** 2 - s ** 2))), s, U, epsabs=0, epsrel=1e-11, limit=400)[0]
    up = (1.0 / s - 1.0 / U) / (3 * (1 - delta) * g0 * b0)
    return float(lo), float(up)


# ---- the numpy [1, K, 1] KAN with RBF edges in both layers (the architecture of the KAN spike, without centre / width learning) ----
class KanToy:
    """x ∈ [−1, 1] → K hidden nodes h_k(x) = Σ_l W1[k,l] φ_l(x), standardised (detached mean / std, std clamped below at 0.05) → Σ_k Σ_m W2[k,m] ψ_m(h_k) + bias.
    Units: the paper's (w = ν² tanh α, m = a_max = 1) in the first layer; second layer 'fact' (the same unit) or 'free' (plain weights W2 = a). Rates: α 1, ν 0.05, plain 0.3."""

    def __init__(self, second, seed=0, K=3, n=200):
        rng = np.random.default_rng(seed)
        self.X = np.sort(rng.uniform(-1, 1, n))
        y = np.sin(3 * self.X) + 0.5 * np.cos(7 * self.X)
        self.y = (y - y.mean()) / y.std()
        self.c1, self.s1, self.c2, self.s2, self.K, self.second = np.linspace(-1, 1, 9), 0.25, np.linspace(-2.5, 2.5, 9), 0.6, K, second
        self.Ph1 = np.exp(-0.5 * ((self.X[:, None] - self.c1) / self.s1) ** 2)
        self.rng = rng

    def init(self, eps):
        K, M, M2 = self.K, 9, 9
        P = dict(a1=np.zeros((K, M)), n1=np.full((K, M), np.sqrt(eps)))
        if self.second == 'fact':
            P.update(a2=np.zeros((K, M2)), n2=np.full((K, M2), np.sqrt(eps)))
        else:
            P.update(w2=0.1 * self.rng.standard_normal((K, M2)))
        return P, 0.0

    def _ph2(self, h):
        return np.exp(-0.5 * ((h[..., None] - self.c2) / self.s2) ** 2)

    def stats(self, P):
        """The (detached) mean and clamped standard deviation of the hidden pre-activations."""
        h = (P['n1'] ** 2 * np.tanh(P['a1'])) @ self.Ph1.T
        return h.mean(1, keepdims=True), np.maximum(h.std(1, keepdims=True), 0.05)

    def grads(self, P, bias, stats=None):
        N = len(self.X)
        W1 = P['n1'] ** 2 * np.tanh(P['a1'])
        W2 = P['n2'] ** 2 * np.tanh(P['a2']) if self.second == 'fact' else P['w2']
        h = W1 @ self.Ph1.T
        mu, sd = self.stats(P) if stats is None else stats
        hs = (h - mu) / sd
        Ph2 = self._ph2(hs)
        out = np.einsum('knm,km->n', Ph2, W2) + bias
        r = self.y - out
        dPh2 = Ph2 * (-(hs[..., None] - self.c2) / self.s2 ** 2)
        gW2 = np.einsum('n,knm->km', r, Ph2) / N
        gW1 = np.einsum('n,kn,nm->km', r, np.einsum('knm,km->kn', dPh2, W2) / sd, self.Ph1) / N
        G = dict(a1=gW1 * P['n1'] ** 2 / np.cosh(P['a1']) ** 2, n1=gW1 * 2 * P['n1'] * np.tanh(P['a1']))
        if self.second == 'fact':
            G.update(a2=gW2 * P['n2'] ** 2 / np.cosh(P['a2']) ** 2, n2=gW2 * 2 * P['n2'] * np.tanh(P['a2']))
        else:
            G.update(w2=gW2)
        return G, float(r.mean()), float(np.sqrt(np.mean(r ** 2))), W1, W2

    def train(self, eps, steps, kind='gd', lr=1.0, nu_ratio=0.05, l2=0.3, ea=1e-8):
        P, bias = self.init(eps)
        rate = dict(a1=lr, n1=lr * nu_ratio, a2=lr * l2, n2=lr * l2 * nu_ratio, w2=lr * l2)
        m = {k: np.zeros_like(v) for k, v in P.items()}
        v2 = {k: np.zeros_like(v) for k, v in P.items()}
        first = None
        for t in range(1, steps + 1):
            G, gb, rms, W1, W2 = self.grads(P, bias)
            if first is None:
                first = max(float(np.abs(g).max()) for g in G.values())
            for k in P:
                if kind == 'gd':
                    P[k] = P[k] + rate[k] * G[k]
                else:
                    m[k] = 0.9 * m[k] + 0.1 * G[k]
                    v2[k] = 0.999 * v2[k] + 0.001 * G[k] ** 2
                    P[k] = P[k] + rate[k] * (m[k] / (1 - 0.9 ** t)) / (np.sqrt(v2[k] / (1 - 0.999 ** t)) + ea)
            bias += 0.1 * lr * gb
        G, gb, rms, W1, W2 = self.grads(P, bias)
        return dict(first_grad=first, rms=rms, w1max=float(np.abs(W1).max()), n_units1=int((np.abs(W1) > 0.01).sum()), w2max=float(np.abs(W2).max()),
                    amax=float(max(np.abs(P['a1']).max(), np.abs(P.get('a2', np.zeros(1))).max())))


# ===============================================================================================================================
# the parts
# ===============================================================================================================================
def part_A():
    """Sign descent and Adam on D3 units with orthonormal atoms and drives G_A: the step at which |w| reaches THW_A, for the starts S_A; GD flow
    time from the closed form (and an Euler check)."""
    G = np.eye(len(G_A))
    c = np.array(G_A)
    out = dict(A_g=c, A_s=np.array(S_A), A_thw=np.array([THW_A]), A_lr=np.array([LR_A]))
    Ts, Ta, Tg = [], [], []
    for s in S_A:
        Ts.append(run_units(Opt('sign', LR_A), G, c, s, NSTEP_A, THW_A)[1])
        Ta.append(run_units(Opt('adam', LR_A, eps_adam=EA_TINY), G, c, s, NSTEP_A, THW_A)[1])
        u_c = (THW_A / 2.0) ** (1.0 / 3.0)                       # w = u³ − v³ ≈ 2u³ − ... : the one-sided flow reference threshold
        Tg.append([gd_escape_time(g, s, u_c) for g in G_A])
    out.update(A_T_sign=np.array(Ts), A_T_adam=np.array(Ta), A_T_gd=np.array(Tg))
    # Euler check of the closed form of the GD flow (one-sided unit, drive 1 and 0.1, s = S_GD_A, h = H_GD_A, threshold u_c = 0.5)
    chk = []
    for g in (1.0, 0.1):
        u, t = S_GD_A, 0
        while u < 0.5:
            u += H_GD_A * 3 * g * u * u
            t += 1
        chk.append((g, t * H_GD_A, gd_escape_time(g, S_GD_A, 0.5)))
    out['A_euler'] = np.array(chk)
    return out


def part_B():
    """Adam on one coordinate with a frozen drive: steps to u_c for the drive scales SCALES_B, the ε_adam values EA_B, and the betas; the extreme R_t."""
    steps = np.zeros((len(EA_B), len(SCALES_B)))
    for i, ea in enumerate(EA_B):
        for j, sc in enumerate(SCALES_B):
            steps[i, j] = frozen_adam_run(S_B, UC_B, LR_B, BETAS[0], BETAS[1], ea, sc)[0]
    rows = []
    for b1, b2 in BETAS_B:
        c, C, tstar = adam_ratio_bounds(b1, b2)
        t, rmin, rmax = frozen_adam_run(S_B, UC_B, LR_B, b1, b2, EA_TINY, 1.0)
        rows.append((b1, b2, c, C, tstar, t, rmin, rmax))
    rng = np.random.default_rng(31)
    worst = 0.0                                                   # random nondecreasing positive sequences: max R_t / C_t, min R_t / (1 − b1)
    low = np.inf
    c, C, _ = adam_ratio_bounds(*BETAS)
    for _ in range(40):
        G = np.cumsum(rng.exponential(1.0, 600) * rng.choice([0.0, 1.0, 50.0], 600, p=[0.5, 0.45, 0.05])) + 1e-3
        G = G ** rng.uniform(0.5, 4.0)
        m = v = 0.0
        for t in range(1, len(G) + 1):
            m = BETAS[0] * m + (1 - BETAS[0]) * G[t - 1]
            v = BETAS[1] * v + (1 - BETAS[1]) * G[t - 1] ** 2
            R = (m / (1 - BETAS[0] ** t)) / np.sqrt(v / (1 - BETAS[1] ** t))
            worst, low = max(worst, R / C), min(low, R / c)
    return dict(B_steps=steps, B_scales=np.array(SCALES_B), B_ea=np.array(EA_B), B_rows=np.array(rows, float), B_rand=np.array([worst, low]))


def part_C():
    """The coupled pair for ρ in RHO_C, s in S_C and the optimizers of OPTS_C; the duplicates."""
    res = []
    for rho in RHO_C:
        for s in S_C:
            for kind, kw in OPTS_C:
                n_steps = int(min(400000, 4.0 / (3 * s * kw['lr'])) if kind == 'gd' else 4000)
                T1, T2, peak, wf = run_pair(kind, kw, rho, s, n_steps, THW_C)
                res.append((rho, s, [o[0] for o in OPTS_C].index(kind), T1, T2, peak, wf[0], wf[1]))
    dup = []
    for kind, kw in OPTS_C[:3]:
        n_steps = int(min(400000, 4.0 / (3 * S_DUP * kw['lr'])) if kind == 'gd' else 4000)
        w = run_duplicates(kind, kw, S_DUP, DUP_RATIO, n_steps)
        dup.append((w[0], w[1]))
    return dict(C_res=np.array(res, float), C_dup=np.array(dup, float))


def part_D():
    """ε_adam scan on the pair (ρ = RHO_D): lr = 5e-4 for ε_adam ≤ 1e-2, else 0.05 ε_adam (the SGD rate lr/ε_adam = 0.05)."""
    rows = []
    for s in S_D:
        for ea in EA_D:
            lr = 5e-4 if ea <= 1e-2 else 0.05 * ea
            n_steps = int(min(1.5e6, max(6000, 80.0 / (3 * s * lr / ea))))
            T1, T2, peak, wf = run_pair('adam', dict(lr=lr, eps_adam=ea), RHO_D, s, n_steps, THW_C)
            rows.append((s, ea, lr, T1, T2, peak, wf[0], wf[1]))
    return dict(D_rows=np.array(rows, float))


def _cell(P, FR, K, f, c0, cs, As, ref, res):
    a = FR.analyse(res, K, f, c0, cs, As, ref=ref)
    return [a['n_formed'], float(a['exact']), float(a['set_equal']), a['fit_rms'], a['gap'] if np.isfinite(a['gap']) else np.nan]


def part_E():
    """The design: Adam and sign descent from the paper's start, formation steps of the 21 units, against GD's drives."""
    P, FR, L, K, c0, cs, As, f = design()
    thr = FR.formed_threshold(As)
    g0 = K.corr(FR.pts(c0), f)[0]
    ref = P.reference(L, K, c0, f)
    out = dict(E_g0=g0, E_thr=np.array([thr]), E_ref_set=np.array(ref['set']), E_ref_n=np.array([len(ref['order'])]))
    names, T = [], []
    for name in ('D3', 'P0'):
        ro = P.make_readout(name)
        for lr in LRS_E:
            for lab, kw in (('adam', dict(beta1=0.9, beta2=0.999, eps_adam=1e-8)), ('sign', dict(beta1=0.0, beta2=0.0, eps_adam=1e-40))):
                res = FR.train_adam(ro, K, f, c0, EPS_E, STEPS_E, lr, n_t=2001, betas=(kw['beta1'], kw['beta2']), eps_adam=kw['eps_adam'])
                T.append(formation_steps(FR, res, thr, STEPS_E))
                names.append((0 if name == 'D3' else 1, lr, 0 if lab == 'adam' else 1))
    out.update(E_T=np.array(T), E_names=np.array(names, float))
    return out


def part_F():
    """The modifications on the design: final state after STEPS_F steps (formed units, exactness, the sliding path's set)."""
    P, FR, L, K, c0, cs, As, f = design()
    ref = P.reference(L, K, c0, f)
    rows, labels = [], []
    for name in ('D3', 'P0'):
        ro = P.make_readout(name)
        rates = np.r_[np.tile(ro.lr, len(c0)), np.full(len(c0), P.LR_C)] if name == 'P0' else np.r_[np.ones(2 * len(c0)), np.ones(len(c0))]
        # the paper's group rates as a per-parameter vector (centres last) for the GD-equivalents
        pre = rates
        specs = [('sign', Opt('sign', 1e-4), name)]
        for ea in EA_F:
            lr = 1e-4 if ea <= 1e-4 else H_F * ea
            specs.append(('adam %g' % ea, Opt('adam', lr * (rates if ea > 1e-4 else 1.0), eps_adam=ea), name))
        specs.append(('adam shared v', Opt('adam_g', 1e-4, pre=pre), name))
        specs.append(('norm. GD', Opt('ngd', 3e-4, pre=pre), name))
        for lab, opt, nm in specs:
            res = train_opt(FR, ro, K, f, c0, EPS_F, STEPS_F, opt, n_t=401)
            rows.append(_cell(P, FR, K, f, c0, cs, As, ref, res))
            labels.append(lab)
    return dict(F_rows=np.array(rows, float), F_labels=np.array(labels), F_ref_n=np.array([len(ref['order'])]))


def part_D2():
    """The crossover in eps_adam on the pair: for rho in RHO_D2 and s in S_D2, the peak |w_2| over the grid eps_adam = m s (m in GRID_D2), lr = 5e-4."""
    peaks = np.zeros((len(RHO_D2) * len(S_D2), len(GRID_D2)))
    cells = []
    for i, rho in enumerate(RHO_D2):
        for j, s in enumerate(S_D2):
            cells.append((rho, s))
            for k, m in enumerate(GRID_D2):
                ea, lr = m * s, 5e-4
                n_steps = int(min(1.5e6, max(6000, 80.0 / (3 * s * lr / ea))))
                peaks[i * len(S_D2) + j, k] = run_pair('adam', dict(lr=lr, eps_adam=ea), rho, s, n_steps, THW_C)[2]
    return dict(D2_cells=np.array(cells), D2_grid=np.array(GRID_D2), D2_peaks=peaks)


def part_H():
    """The paper's unit (P0) with a large eps_adam on the design: the horizon of its wall (3e6 steps at the GD step 0.05 of its group rates)."""
    P, FR, L, K, c0, cs, As, f = design()
    ref = P.reference(L, K, c0, f)
    ro = P.make_readout('P0')
    rates = np.r_[np.tile(ro.lr, len(c0)), np.full(len(c0), P.LR_C)]
    rows = []
    for ea in EA_H:
        res = train_opt(FR, ro, K, f, c0, EPS_F, STEPS_H, Opt('adam', H_F * ea * rates, eps_adam=ea), n_t=401)
        rows.append(_cell(P, FR, K, f, c0, cs, As, ref, res))
    return dict(H_rows=np.array(rows, float), H_ea=np.array(EA_H))


class AddKan:
    """The additive KAN [d, 1] of the spike in numpy: f(x) = bias + Σ_jl w_jl exp(−½((x_j − c_jl)/σ_jl)²), w = ν² A tanh(α/A) (A = 3), ν(0)² = eps, α(0) = 0, centres
    and log-widths trained, log σ clamped to [log 0.05, log 3]; the loss is ½ mean r²."""

    def __init__(self, d, M, eps):
        from exp_kan_spike import X_LO, X_HI
        grid = np.linspace(X_LO, X_HI, M)
        self.keys = ('a', 'n', 'c', 'ls', 'b')
        self.p = dict(a=np.zeros((d, M)), n=np.full((d, M), np.sqrt(eps)), c=np.tile(grid, (d, 1)), ls=np.full((d, M), np.log(grid[1] - grid[0])), b=np.zeros(1))
        self.A = 3.0

    def fwd(self, X):
        p = self.p
        s = np.exp(np.clip(p['ls'], np.log(0.05), np.log(3.0)))
        u = (X[:, :, None] - p['c'][None]) / s[None]
        phi = np.exp(-0.5 * u * u)
        w = p['n'] ** 2 * self.A * np.tanh(p['a'] / self.A)
        return (phi * w[None]).sum((1, 2)) + p['b'][0], phi, u, s, w

    def grad(self, X, y):
        p, A = self.p, self.A
        f, phi, u, s, w = self.fwd(X)
        r = y - f
        gw = (r[:, None, None] * phi).mean(0)
        inside = (p['ls'] > np.log(0.05)) & (p['ls'] < np.log(3.0))
        th = np.tanh(p['a'] / A)
        return dict(a=gw * p['n'] ** 2 * (1 - th ** 2), n=gw * 2 * p['n'] * A * th, c=(r[:, None, None] * phi * u / s[None]).mean(0) * w,
                    ls=(r[:, None, None] * phi * u * u).mean(0) * w * inside, b=np.array([r.mean()]))

    def vec(self, d):
        return np.concatenate([d[k].ravel() for k in self.keys])

    def setvec(self, z):
        i = 0
        for k in self.keys:
            n = self.p[k].size
            self.p[k] = z[i:i + n].reshape(self.p[k].shape)
            i += n

    def counts(self, X):
        from exp_kan_spike import layer_activity, edge_counts
        p = self.p
        w = p['n'] ** 2 * self.A * np.tanh(p['a'] / self.A)
        s = np.exp(np.clip(p['ls'], np.log(0.05), np.log(3.0)))
        return edge_counts(layer_activity(w[None], p['c'][None], s[None], X)[2])


def fit_add(model, data, opt, steps):
    """The spike's `fit` for the additive model: full batch, validation every 25 steps, best iterate kept, patience 2000 steps counted from the first time the
    validation MSE is below 0.9 of its initial value."""
    Xtr, ytr, Xva, yva = data['Xtr'], data['ytr'], data['Xva'], data['yva']
    best, best_p, v0, armed, since = np.inf, None, None, False, 0
    for t in range(steps + 1):
        if t % 25 == 0:
            v = float(np.mean((model.fwd(Xva)[0] - yva) ** 2))
            if v0 is None:
                v0 = v
            if v < best:
                best, since, best_p = v, 0, {k: x.copy() for k, x in model.p.items()}
            else:
                since += 25
            if v < 0.9 * v0:
                armed = True
            if armed and since >= 2000:
                break
        if t == steps:
            break
        model.setvec(opt.step(model.vec(model.p), model.vec(model.grad(Xtr, ytr))))
    model.p = best_p
    return t


def part_K():
    """GD against Adam (and the modifications) on the additive KAN, the spike's development targets (numpy re-implementation; centres and widths trained)."""
    import exp_kan_spike as KS
    counts, truth, steps = [], [], []
    for tid in TARGETS_K:
        tg = KS.make_target(tid)
        for noise in NOISES_K:
            data = KS.sample_synth(tg, noise, SEED_K)
            truth.append(tg['K'])
            c_o, s_o = [], []
            for lab, kind, ea in OPTS_K:
                m = AddKan(KS.D_SYN, M_K, EPS_K)
                rates = {k: 1.0 for k in m.keys}
                rates.update(n=0.05, b=0.1)
                rv = m.vec({k: np.full_like(m.p[k], rates[k]) for k in m.keys})
                opt = {'gd': Opt('gd', rv), 'adam': Opt('adam', 0.01, eps_adam=ea), 'adam_r': Opt('adam', ea * rv, eps_adam=ea), 'adam_g': Opt('adam_g', 0.003, eps_adam=ea)}[kind]
                s_o.append(fit_add(m, data, opt, STEPS_K))
                c_o.append(m.counts(data['Xtr']))
            counts.append(c_o)
            steps.append(s_o)
    return dict(K_counts=np.array(counts), K_truth=np.array(truth), K_steps=np.array(steps),
                K_cells=np.array([(t, n) for t in TARGETS_K for n in NOISES_K]))


def part_G():
    out = {}
    # chains: quadrature against the ODE, and the leading term
    g = 1.0
    Tq = np.full((len(D_G), len(S_G)), np.nan)
    To = np.full_like(Tq, np.nan)
    Tl = np.full_like(Tq, np.nan)
    for i, D in enumerate(D_G):
        for j, s in enumerate(S_G):
            if (D == 6 and s < 0.1) or (D == 5 and s < 0.05):
                continue
            z_t = (0.5 * g) ** (1.0 / D)
            Tq[i, j] = chain_time_quad(D, g, s, z_t)
            To[i, j] = chain_time_ode(D, g, s, z_t)
            Tl[i, j] = chain_time_leading(D, g, s)
    out.update(G_Tq=Tq, G_To=To, G_Tl=Tl)
    # the diagonal two-layer network: closed form against Euler, the times for small starts, and the pair of drives (1, 1/2)
    eul = []
    for gg, s in ((1.0, 0.1), (0.5, 0.1), (1.0, 0.01)):
        t_end = diag_two_layer_time(gg, s * s) + 2.0
        eul.append((gg, s, t_end, diag_gd_check(gg, s, 1e-3, t_end), float(diag_two_layer_w(t_end, gg, s * s))))
    out['G_diag_euler'] = np.array(eul)
    rows = []
    for s in S_DIAG:
        t1, t2 = diag_two_layer_time(1.0, s * s), diag_two_layer_time(0.5, s * s)
        # the D3 unit, same drives: the flow time to u_c = 0.5 (exact), for the ratio gap / transit
        d1, d2 = gd_escape_time(1.0, s, 0.5), gd_escape_time(0.5, s, 0.5)
        rows.append((s, t1, t2, d1, d2))
    out['G_diag_rows'] = np.array(rows)
    # the paper-unit pair: stationarity, the time law
    out['G_pair_stat'] = np.array(pair_rhs(0.0, [0.0, 0.3, 0.0, 0.3], 1.0))
    pr = []
    for s in S_PAIR:
        for a0 in A0_PAIR:
            pr.append((s, a0, pair_time(s, a0), pair_time_sym(s, a0 / 2.0), pair_time_sym(s, a0)))
    out['G_pair'] = np.array(pr)
    # the hybrid
    hy = []
    for s in S_H:
        for b0 in B0_H:
            b0v = np.sqrt(1.0 / 3.0) * s if b0 is None else b0
            T, cons, delta, b, g0, npos = hybrid_run(s, b0v, U_H)
            lo, up = hybrid_bounds(s, b0v, U_H, g0, delta, npos)
            hy.append((s, b0v, T, lo, up, delta, cons, g0, npos))
    out['G_hybrid'] = np.array(hy)
    # the KAN toy
    kan = []
    for second, kind in (('fact', 'gd'), ('fact', 'adam'), ('free', 'gd'), ('free', 'adam')):
        r = KanToy(second).train(EPS_KAN, STEPS_KAN, kind=kind, lr=(1.0 if kind == 'gd' else 0.01))
        kan.append((0 if second == 'fact' else 1, 0 if kind == 'gd' else 1, r['first_grad'], r['rms'], r['w1max'], r['n_units1'], r['w2max'], r['amax']))
    out['G_kan'] = np.array(kan)
    return out


PARTS = dict(A=part_A, B=part_B, C=part_C, D=part_D, D2=part_D2, E=part_E, F=part_F, G=part_G, H=part_H, K=part_K)


# ===============================================================================================================================
# the report
# ===============================================================================================================================
def _f(x, n=3):
    return 'nan' if x is None or (isinstance(x, float) and np.isnan(x)) else ('%.*g' % (n, x))


def _tab(head, rows):
    L = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    return L + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]


def report(R):
    L = ['# Normalised optimizers erase the ordering; the two-layer small start does not leave its saddle', '',
         'Generated by `experiments/theory3f_escape.py --from-raw`; proofs in the round note (outside the repository); tests `tests/test_theory3f_escape.py`.', '',
         '## Statements and their status', '',
         '* **F1a (sign descent; PROVED, checked below).** Units `w_i = W(p_i)` with a common start and common step sizes, trained by `p <- p + lr*sign(-dL/dp)`; let `T` be the first step at which '
         '`|W(p(t))| >= theta` along the path of a unit whose drive `g_i(t) = <r phi_i>` keeps its sign. Then every unit whose drive keeps its sign up to `T` has the SAME parameter path, hence crosses '
         'at the same step, whatever the ratio of the drives; for D3 (`w = u^3 - v^3`, start `(s, s)`) `T = min{t: 2(3 s^2 lr t + (lr t)^3) >= theta}`, which tends to `ceil((theta/2)^(1/3)/lr)` as `s -> 0`: '
         'the small start does not set the escape time. Hypotheses: `|g_i(0)| > theta*sum_j |G_ij|` (the drive cannot turn before `theta`); a generic start (if `v` hits exactly 0, `dw/dv = 0` and the sign update stops that coordinate: the runs use starts that are not multiples of the step). With starts `s_i` the steps differ by at most `eps |b_i - b_j|/(a_0 lr) + 1` (`s_i^2 = eps b_i`, `a_0 = (theta/2)^(1/3)`).',
         '* **F1b (Adam, `eps_adam = 0`; PROVED).** The Adam update is invariant under a positive rescaling of the gradient history, so units with PROPORTIONAL gradient histories (a frozen drive) have identical '
         'trajectories (exact). For positive nondecreasing gradient sequences `1 - beta1 <= R_t <= C(beta1, beta2)` (`C = 7.27` for `(0.9, 0.999)`), so the escape steps of any two such units are within the factor '
         '`C/(1 - beta1)`, independent of the drives. With `eps_adam > 0` the same holds while `eps_adam <= kappa*sqrt(vhat_t)` (up to `1/(1 + kappa)`); the dormant gradient is `3 g s^2` for D3.',
         '* **F1c (the pair, duplicates; sign descent).** Atom 1 is the target, atom 2 has overlap `rho`. PROVED: both units follow the same path (`w_1 = w_2`) until the drive of atom 2 turns, so '
         '`peak |w_2| >= rho/(1 + rho)` (up to one step\'s increment); gradient descent keeps atom 2 at the dormant scale. NUMERICAL: afterwards atom 2 slides along `w_2 = rho(1 - w_1)` and the end '
         'state is the sparse fit up to the chatter of size `lr` (the path, not the limit point, is wrong). PROVED for duplicate atoms (equal drives at every step): `u_1 - u_2 = s_1 - s_2` for all time, '
         '`w_2/w_1 = (3 s_2^2 + a^2)/(3 s_1^2 + a^2)` with a common `a`, so the amplitude is shared (ratio `1 - O(s^2)`); gradient descent: the larger start takes all.',
         '* **F1d (what keeps the mechanism).** PROVED: (0) any update of the form `p_i += lr*Phi(history of the gradient of unit i)` with `Phi` invariant under positive rescaling of its argument (sign descent, '
         'Adam with `eps_adam = 0`, RMSprop, per-coordinate normalisation in ANY coordinates, in particular in the mirror coordinates `eta` where the flow is `eta_i\' = g_i`; per-group step sizes only rescale `lr`) '
         'gives identical trajectories to units with proportional gradient histories: the ordering by the drive is lost; (i) GD: escape time `(1/s - 1/u_c)/(3 g)`, ordered by `g`; '
         '(ii) Adam with `eps_adam >= K max_t sqrt(vhat_t)` is momentum SGD with the rate `lr/eps_adam` up to the factor `1 + 1/K` (inequality on the update); (iii) normalised gradient descent with the GLOBAL norm '
         'is gradient flow run on the clock `d tau = dt/||grad||` (same orbits); (iv) one SHARED scalar second moment (`adam_g`, `beta1 = 0`) is gradient descent with a positive variable step: same orbits '
         '(a scale shared by all units keeps the relative speeds `g_i dW/dp`). NUMERICAL (D, D2, F, H, K): a large `eps_adam` restores the sliding path on the pair (crossover `eps* ~ c rho s/(1 - rho)`), on the design (D3: the number of formed units falls 21, 10, 8, 5, 3 as '
         '`eps_adam` rises from 1e-8 to 3e-3, then exactly the path\'s set; the paper\'s unit P0 at the long horizon of its wall) and on the additive KAN (exact-edge share 0.60 at `eps_adam = 1`, 0.43 at 0.01, '
         'against 0.67 for gradient descent and 0.27 for Adam); the shared second moment restores it on the pair (for rho = 0.8 and s = 0.013 the peak of atom 2 is 1.0e-3, at the formed threshold, against 2.5e-4 for gradient descent), on the design (D3 and P0) and on the additive KAN (0.70).',
         '* **F2a (the product of D small factors; PROVED, quadrature).** Balanced start `z = s`, `w = z^D`, drive `g`: time to `z_theta` is the quadrature in the module docstring: `ln(1/s)/g + O(1)` for `D = 2` '
         '(the diagonal two-layer network, `w` logistic: `w(t) = g/(1 + (g/w_0 - 1) exp(-2 g t))`) and `s^-(D-2)/((D-2) g) + O(1)` for `D >= 3`. With `s^2 = eps b` (the repo): `eps^(-(D-2)/2)`; `D = 3` is the repo\'s '
         'unit (`eps^(-1/2)`), `D = 4` (a cubic unit times a linear second layer) `eps^(-1)`, `D = 6` (two cubic layers) `eps^(-2)`.',
         '* **F2b (the paper\'s start in two layers; PROVED).** If the coefficients of both layers are zero (`w = 0`: alpha = 0 for the paper\'s unit, `(s, s)` for D3) and the target is centred, the vector field vanishes '
         '(each layer\'s gradient is proportional to the other layer\'s coefficient or to the correlation of the target with the other layer\'s output): the whole set is stationary, for every optimizer with `eps_adam > 0`. '
         'A perturbation `a_0` of one layer gives `T = (ln(s/a_0) + O(1))/(lr c* s^4)` for the unsaturated pair (exact quadrature for the symmetric start; the pair conserves `alpha_2^2 - alpha_1^2`, so the start `(0, a_0)` is the symmetric start `a_0/2` up to `O(a_0^2)`: equal to 9 digits at `a_0 = 1e-6`).',
         '* **F2c (the hybrid; PROVED).** First layer cubic units `u_kl^3` from `s`, second layer `b_k(0) = b_0 > 0`, rates 1: `b_k^2 - (1/3) sum_l u_kl^2 = b_0^2 - (1/3) n s^2` for all times, and the escape of the leading '
         'unit satisfies `int_s^U du/(3(1+delta) g u^2 b_max(u)) <= T <= (1/s - 1/U)/(3(1-delta) g b_0)`, `delta` = the relative variation of its drive, `b_max(u)^2 = b_0^2 + (n_+/3)(u^2 - s^2)`: the one-layer escape with '
         'the effective drive `b_0 g` (exponent `eps^(-1/2)` kept) when `b_0^2 >> U^2 n_+/3`.', '']
    # ---- A
    A = R.get('A_g')
    if A is not None:
        L += ['## A. Sign descent and Adam on the diagonal model (orthonormal atoms)', '',
              'Six D3 units, drives `g = %s`, coefficient threshold %s, step size %s, %d steps. Entries: the step at which `|w_i|` reaches the threshold (sign descent, Adam with `eps_adam = 1e-30`, default betas), and the plain gradient flow time of the closed form `(1/s - 1/u_c)/(3g)`.' %
              (', '.join(_f(float(x)) for x in R['A_g']), _f(float(R['A_thw'][0])), _f(float(R['A_lr'][0])), NSTEP_A), '']
        rows = []
        for i, s in enumerate(R['A_s']):
            for j, g in enumerate(R['A_g']):
                rows.append((_f(float(s)), _f(float(g)), _f(float(R['A_T_sign'][i, j])), _f(float(R['A_T_adam'][i, j])), _f(float(R['A_T_gd'][i, j]))))
        L += _tab(['s', 'drive g_i', 'sign: step', 'Adam: step', 'GD: flow time'], rows) + ['']
        Ts, Ta, Tg = R['A_T_sign'], R['A_T_adam'], R['A_T_gd']
        L += ['Window (max / min over the six units): sign descent %s, Adam %s, gradient descent %s (= g_max/g_min = %s). Sign-descent steps identical across the starts s: %s.' %
              (', '.join(_f(float(Ts[i].max() / Ts[i].min())) for i in range(len(Ts))), ', '.join(_f(float(Ta[i].max() / Ta[i].min())) for i in range(len(Ta))),
               ', '.join(_f(float(Tg[i].max() / Tg[i].min())) for i in range(len(Tg))), _f(float(max(R['A_g']) / min(R['A_g']))),
               'yes' if all(np.array_equal(Ts[0], Ts[i]) for i in range(1, len(Ts))) else 'no (the start enters through s^2 lr t)'),
              'Euler check of the closed form (drive, steps*h, formula): %s.' % '; '.join('g = %s: %s vs %s' % (_f(float(a)), _f(float(b)), _f(float(c))) for a, b, c in R['A_euler']), '']
    B = R.get('B_steps')
    if B is not None:
        L += ['## B. Adam on one coordinate with a frozen drive', '',
              'Gradient `scale*3u^2`, start `s = %s`, `lr = %s`, steps to reach `u = %s` (default betas). Rows: `eps_adam`; columns: the drive scale.' % (_f(S_B), _f(LR_B), _f(UC_B)), '']
        L += _tab(['eps_adam \\ scale'] + [_f(float(x)) for x in R['B_scales']], [[_f(float(e))] + [_f(float(x)) for x in B[i]] for i, e in enumerate(R['B_ea'])]) + ['']
        L += ['The window theorem\'s constants and the extreme values of `R_t` on the frozen-drive trajectory:', '']
        L += _tab(['beta1', 'beta2', 'c = 1 - beta1', 'C (bound)', 'argmax t', 'steps', 'min R_t', 'max R_t'],
                  [[_f(float(x), 5) for x in r] for r in R['B_rows']]) + ['']
        L += ['Random positive nondecreasing sequences (40, 600 steps, default betas): max `R_t / C_t` = %s (must be <= 1), min `R_t / c` = %s (must be >= 1).' % (_f(float(R['B_rand'][0])), _f(float(R['B_rand'][1]))), '']
    C = R.get('C_res')
    if C is not None:
        L += ['## C. The coupled pair and duplicate atoms', '',
              'Gram `[[1, rho], [rho, 1]]`, target = atom 1 (drives `1`, `rho`); D3 units; coefficient threshold %s. `T1`, `T2`: first step with `|w| >= threshold` (nan: never); `peak`: max `|w_2|`; final `w`. GD: step 0.05.' % _f(THW_C), '']
        rows = []
        for r in C:
            rows.append((_f(float(r[0])), _f(float(r[1])), OPTS_C[int(r[2])][0], _f(float(r[3])), _f(float(r[4])), _f(float(r[5])), _f(float(r[6])), _f(float(r[7]))))
        L += _tab(['rho', 's', 'optimizer', 'T1', 'T2', 'peak w2', 'final w1', 'final w2'], rows) + ['']
        L += ['Prediction for sign descent: peak `w_2 = rho/(1 + rho)` = %s (rho = 0.2, 0.5, 0.8).' % ', '.join(_f(r / (1 + r)) for r in RHO_C), '']
        Dd = R['C_dup']
        L += ['Duplicate atoms (`phi_1 = phi_2`, target amplitude 1, starts `s_1 = %s`, `s_2 = %s s_1`): final `(w_1, w_2)` and `w_2/w_1`:' % (_f(S_DUP), _f(DUP_RATIO)), '']
        L += _tab(['optimizer', 'w1', 'w2', 'w2/w1'], [(OPTS_C[i][0], _f(float(a), 4), _f(float(b), 4), _f(float(b / a), 4)) for i, (a, b) in enumerate(Dd)]) + ['']
    D = R.get('D_rows')
    if D is not None:
        L += ['## D. Does a large `eps_adam` restore the selection? (the pair, rho = %s)' % _f(RHO_D), '',
              'D3 units, Adam, `lr = 5e-4` for `eps_adam <= 1e-2` and `lr = 0.05 eps_adam` above (the SGD rate `lr/eps_adam` = 0.05). Selection restored = `peak w2` below the threshold %s (atom 2 never forms).' % _f(THW_C), '']
        rows = [(_f(float(r[0])), _f(float(r[1])), _f(float(r[2])), _f(float(r[3])), _f(float(r[4])), _f(float(r[5])), _f(float(r[6]), 4), _f(float(r[7]), 4)) for r in D]
        L += _tab(['s', 'eps_adam', 'lr', 'T1', 'T2', 'peak w2', 'final w1', 'final w2'], rows) + ['']
    D2 = R.get('D2_peaks')
    if D2 is not None:
        thr_ = THW_C
        L += ['### The crossover in `eps_adam` (pair, `lr = 5e-4`, D3)', '',
              '`eps*`: the smallest grid value above which the peak of `|w_2|` stays below the threshold %s (atom 2 never forms), `nan` if none; grid `eps_adam = m s`, `m` = %s.' % (_f(thr_), ', '.join(_f(float(m)) for m in R['D2_grid'])), '']
        rows = []
        for (rho, s), pk in zip(R['D2_cells'], D2):
            bad = np.flatnonzero(pk >= thr_)
            idx = int(bad.max()) + 1 if bad.size else 0
            m_star = float(R['D2_grid'][idx]) if idx < len(pk) else float('nan')
            rows.append((_f(float(rho)), _f(float(s)), _f(m_star), _f(m_star * float(s)), _f(float(pk[0])), _f(float(pk[-1]))))
        L += _tab(['rho', 's', 'eps*/s', 'eps*', 'peak w2 at the smallest eps_adam', 'peak w2 at the largest'], rows) + ['']
    E = R.get('E_T')
    if E is not None:
        g0 = R['E_g0']
        L += ['## E. The repo\'s design: Adam and sign descent form all 21 units together', '',
              'Dictionary of 21 Gaussian kernels (width 0.12), features at (-0.52, 0.07, 0.61) with amplitudes (0.9, -0.6, 0.75), `eps = %s`, %d steps, one step size for all parameters and centres '
              '(`factorized_rbf.train_adam`). The sliding path (the limit of gradient descent) forms %d units. Drives at the start (`<f* phi_i>`): min `|g|` %s, max `|g|` %s (ratio %s); gradient descent would '
              'separate the escape times of the extreme units by this ratio. Formation = first step with `|w| >= %s`.' % (_f(EPS_E), STEPS_E, int(R['E_ref_n'][0]), _f(float(np.abs(g0).min())), _f(float(np.abs(g0).max())), _f(float(np.abs(g0).max() / np.abs(g0).min())), _f(float(R['E_thr'][0]))), '']
        rows = []
        for nm, T in zip(R['E_names'], E):
            ok = ~np.isnan(T)
            rows.append(('D3' if nm[0] == 0 else 'P0', _f(float(nm[1])), 'Adam' if nm[2] == 0 else 'sign', int(ok.sum()), _f(float(np.nanmin(T))), _f(float(np.nanmedian(T))), _f(float(np.nanmax(T))),
                         _f(float(np.nanmax(T) / np.nanmin(T)))))
        L += _tab(['readout', 'lr', 'optimizer', 'units formed (of 21)', 'first step', 'median step', 'last step', 'last/first'], rows) + ['']
    F = R.get('F_rows')
    if F is not None:
        labs = [str(x) for x in R['F_labels']]
        L += ['## F. Modifications on the design (%d steps, eps = %s)' % (STEPS_F, _f(EPS_F)), '',
              'The sliding path (gradient descent in the limit) forms %d units, exactly. Adam with `eps_adam > 1e-4` uses `lr = %s eps_adam` times the readout\'s group rates (so `lr/eps_adam` is the GD step %s of the paper\'s rates); '
              '`adam shared v` (step 1e-4) and `norm. GD` (step 3e-4) act on the gradient multiplied by the group rates.' % (int(R['F_ref_n'][0]), _f(H_F), _f(H_F)), '']
        n = len(labs) // 2
        rows = []
        for i, lab in enumerate(labs):
            r = F[i]
            rows.append(('D3' if i < n else 'P0', lab, int(r[0]), 'yes' if r[1] else 'no', 'yes' if r[2] else 'no', _f(float(r[3]))))
        L += _tab(['readout', 'optimizer', 'units formed', 'exact', 'set = path', 'rms residual'], rows) + ['']
    Hh = R.get('H_rows')
    if Hh is not None:
        L += ['The paper\'s unit (P0) needs the horizon of its wall (about 7000/drive time units, 1.5e6 steps at the GD step %s): %d steps, `eps = %s`:' % (_f(H_F), STEPS_H, _f(EPS_F)), '']
        L += _tab(['readout', 'eps_adam', 'units formed', 'exact', 'set = path', 'rms residual'],
                  [('P0', _f(float(e)), int(r[0]), 'yes' if r[1] else 'no', 'yes' if r[2] else 'no', _f(float(r[3]))) for e, r in zip(R['H_ea'], Hh)]) + ['']
    G = R.get('G_Tq')
    if G is not None:
        L += ['## G. Two layers', '', '### Products of D small factors (drive 1; time from s to z_theta = (1/2)^(1/D); quadrature = ODE)', '']
        rows = []
        for i, D in enumerate(D_G):
            for j, s in enumerate(S_G):
                if not np.isnan(G[i, j]):
                    rows.append((D, _f(s), _f(float(G[i, j]), 6), _f(float(R['G_To'][i, j]), 6), _f(float(R['G_Tl'][i, j])), _f(float(G[i, j] / R['G_Tl'][i, j]))))
        L += _tab(['D', 's', 'quadrature', 'ODE', 'leading term', 'ratio'], rows) + ['']
        slopes = []
        for i, D in enumerate(D_G):
            ok = ~np.isnan(R['G_Tq'][i])
            ss = np.array(S_G)[ok]
            tt = R['G_Tq'][i][ok]
            slopes.append((D, _f(float(np.polyfit(np.log(ss[-3:]), np.log(tt[-3:]), 1)[0])), _f(-(D - 2.0)) if D > 2 else 'log'))
        L += ['Fitted exponents of `T` in `s` (last three starts) against `-(D - 2)`:', ''] + _tab(['D', 'fitted', 'predicted'], slopes) + ['']
        e = R['G_diag_euler']
        L += ['Diagonal two-layer network, closed form against Euler (step 1e-3) at `t_end`:', '']
        L += _tab(['g', 's', 't_end', 'Euler w', 'closed form w'], [[_f(float(x), 6) for x in r] for r in e]) + ['']
        L += ['Escape times of two units of drives 1 and 1/2 (to `w = g/2` for the two-layer network, to `u = 1/2` for the cubic unit); the gap between the two units against the duration of one transit:', '']
        rows = [(_f(float(r[0])), _f(float(r[1])), _f(float(r[2])), _f(float(r[2] - r[1])), _f(float(r[3])), _f(float(r[4])), _f(float(r[4] - r[3]))) for r in R['G_diag_rows']]
        L += _tab(['s', 'T (g=1)', 'T (g=1/2)', 'gap', 'cubic T (g=1)', 'cubic T (g=1/2)', 'cubic gap'], rows) + ['']
        L += ['### The paper\'s unit in two layers (unsaturated, rates 1, c* = 1)', '', 'Vector field at `alpha_1 = alpha_2 = 0, nu = 0.3`: %s (exactly zero).' % ', '.join(_f(float(x)) for x in R['G_pair_stat']), '']
        rows = [(_f(float(r[0])), _f(float(r[1])), _f(float(r[2])), _f(float(r[2] * r[0] ** 4)), _f(float(r[3] * r[0] ** 4)), _f(float(r[4] * r[0] ** 4)), _f(float(np.log(r[0] / r[1])))) for r in R['G_pair']]
        L += _tab(['s', 'a_0', 'T (alpha_2 = a_0)', 'T s^4', 'sym. quadrature (a_0/2) T s^4', 'sym. quadrature (a_0) T s^4', 'ln(s/a_0)'], rows) + ['']
        H = R['G_hybrid']
        L += ['### The hybrid: cubic first layer, second layer b_0 (RBF net with 21 kernels, 2 hidden nodes, feature at 0.3, escape of the leading unit to u = %s)' % _f(U_H), '']
        rows = [(_f(float(r[0])), _f(float(r[1])), _f(float(r[2])), _f(float(r[3])), _f(float(r[4])), 'yes' if r[3] <= r[2] <= r[4] else 'no', _f(float(r[5])), _f(float(r[6]))) for r in H]
        L += _tab(['s', 'b_0', 'T (flow)', 'lower bound', 'upper bound', 'inside', 'delta', 'conservation residual'], rows) + ['']
        K = R['G_kan']
        L += ['### A numpy [1, 3, 1] KAN with RBF edges in both layers (eps = %s, %d steps)' % (_f(EPS_KAN), STEPS_KAN), '']
        rows = [('factorised (pure small start)' if r[0] == 0 else 'plain weights (hybrid)', 'GD' if r[1] == 0 else 'Adam', _f(float(r[2])), _f(float(r[3])), _f(float(r[4])), int(r[5]), _f(float(r[6])), _f(float(r[7]))) for r in K]
        L += _tab(['second layer', 'optimizer', 'max grad at start', 'rms residual', 'max |w1|', 'first-layer units > 0.01', 'max |w2|', 'max |alpha|'], rows) + ['']
    Kc = R.get('K_counts')
    if Kc is not None:
        truth = R['K_truth']
        L += ['## K. The additive KAN [5, 1] of the spike in numpy: gradient descent against Adam', '',
              'The spike\'s development targets %s, noise sd %s, seed %d, %d units per edge, `eps = %s`, the spike\'s stopping rule (validation every 25 steps, patience 2000, at most %d steps), '
              'the paper\'s unit with trained centres and widths. The rule of exactness: the number of active units of every edge equals the number of true bumps (threshold 0.02 of the output standard deviation). '
              'Adam with `eps_adam = 1` and `eps_adam = 0.01` use `lr = eps_adam` times the group rates of gradient descent (so `lr/eps_adam` = the GD rates). The spike\'s torch records: exact-edge share '
              '0.70 (GD, `kan_dev.md`) against 0.23-0.28 (four Adam configurations).' % (', '.join(str(int(t)) for t in TARGETS_K), ', '.join(_f(x) for x in NOISES_K), SEED_K, M_K, _f(EPS_K), STEPS_K), '']
        rows = []
        for i, (lab, kind, ea) in enumerate(OPTS_K):
            cnt = Kc[:, i, :]
            exact = (cnt == truth)
            rows.append((lab, _f(float(exact.mean())), _f(float(np.maximum(cnt - truth, 0).sum(1).mean())), _f(float(np.maximum(truth - cnt, 0).sum(1).mean())),
                         _f(float(cnt.sum(1).mean())), _f(float(R['K_steps'][:, i].mean()))))
        L += _tab(['optimizer', 'exact-edge share', 'over-counted units / run', 'under-counted units / run', 'active units / run', 'mean steps'], rows) + ['']
        L += ['Counts per edge, noise 0.1 (true counts first):', '']
        rows = []
        for j, (tid, nz) in enumerate(R['K_cells']):
            if abs(nz - NOISES_K[0]) < 1e-12:
                rows.append([str(int(tid)), str(truth[j].tolist())] + [str(Kc[j, i].tolist()) for i in range(len(OPTS_K))])
        L += _tab(['target', 'true'] + [o[0] for o in OPTS_K], rows) + ['']
    return '\n'.join(L) + '\n'


# ===============================================================================================================================
def load_raw():
    with np.load(RAW) as z:
        return {k: z[k] for k in z.files}


def save_raw(raw):
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(RAW, **raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--parts', default=','.join(PARTS))
    a = ap.parse_args()
    if a.run:
        warnings.filterwarnings('ignore')
        raw = load_raw() if os.path.exists(RAW) else {}
        for p in a.parts.split(','):
            t0 = time.time()
            raw.update(PARTS[p]())
            save_raw(raw)
            print('part', p, 'done in %.0f s' % (time.time() - t0), flush=True)
    if a.run or a.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
