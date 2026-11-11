"""
Which factorised readouts of a Gaussian-kernel (RBF) network, trained by plain gradient descent from a small start, choose their own units?
====================================================================================================================================

The free mass–amplitude law (Paper 17, Lemma 2.1(o); `mirror_free.FreeLaw`, `rich_transit`) is the gradient flow of the fit
F = ½⟨(f* − f̂)²⟩ for the metric (Gm)⁻¹dα² + (γmβ)⁻¹dβ² + (Gm)⁻¹|dc|². With ν = √β the middle term is (4/(γm))dν², so in (α, ν, c) it is
the EUCLIDEAN gradient flow with the group learning rates (Gm, γm/4, Gm) of the network

    f̂(x) = Σ_i w_i φ(x; c_i),   w_i = ν_i² a_max tanh(α_i/a_max) / m,   α_i(0) = 0,  ν_i(0)² = ε b_i,

a standard RBF network with a factorised readout (`PaperLaw`; checked against `rich_transit.moving_rhs` in tests/test_factorized_rbf.py). The
ratio of the learning rates fixes the law's one dimensionless constant, v_s/a_max² = lr_ν/lr_α = γ/(4G).

This module trains the OTHER readouts the same way — every parameter p_i of unit i by ṗ_i = −lr·∂F/∂p_i, the centre by ċ_i = −lr_c·∂F/∂c_i
= lr_c w_i⟨r∇_cφ_i⟩ — to ask which part of the transit (units enter one at a time, an entrant slides onto its feature, the others stay
put, the formed set is the sliding inverse scale space's) belongs to the paper's unit and which to any factorised readout:

    D1  `Linear`       w = p                         depth 1, no factorisation (the control)
    D2  `Product2`     w = u² − v², u(0) = v(0) = s   depth 2 (the signed diagonal linear network)
    D3  `Cubic3`       w = u³ − v³, u(0) = v(0) = s   depth 3, no saturation, no free mass
    P0  `PaperLaw`     w = ν²a_max tanh(α/a_max)/m     the paper's unit (an effective depth-3 unit at small mass)
    P0u `Unsaturated`  w = ν²α/m                      the paper's unit without the amplitude's saturation
    P0p `PaperPairs`   P0 started on branches            two units per atom, α(0) = +a₀ and −a₀ (Paper 17, Proposition 'starts on branches')

with s² = ε b_i throughout, so one ε labels the start of all of them. Each factorisation conserves a quantity along the flow (P0: Q = β −
v_s cosh(2α/a_max); D2: uv = s²; D3: 1/u + 1/v = 2/s), which makes it a mirror flow η̇_i = ⟨rφ_i⟩, w_i = W(η_i) with its own mirror map:
D2's wall is at η ≈ ln(1/ε)/(4 lr) with the dormant coefficient ε^μ at the margin μ of the wall (a logarithmic scale), D3's and P0's at
∝ ε^{−1/2} with the dormant coefficient ∝ ε^{3/2}μ^{−3} (a power scale). `Readout.scale(ε)` is that wall's height, so τ = t/scale is the
sliding path's time in every case.

Also here: fixed-step gradient descent (`train_gd`), Adam (`train_adam`), and the measures that compare a trained network with the
sliding path (`analyse`): the formed set, the order and the times of formation, the gap between formed and dormant coefficients, the
dormant units' displacement, exactness.
"""
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from src.core.rich_transit import transit_scale

THETA = 0.01            # a unit has *formed* when |w| > THETA·max|A*| (Paper 17, Section 8)
STIFF_BUDGET = 150000   # right-hand-side evaluations of one LSODA piece before it is redone by Radau (`rich_transit.relax_chunk`: LSODA's stiffness
                        # detection can fail on long horizons with fast modes — 6.7·10⁶ evaluations against Radau's 100)


class _OverBudget(Exception):
    pass


def _counted(f, budget):
    n = [0]

    def g(t, y):
        n[0] += 1
        if n[0] > budget:
            raise _OverBudget
        return f(t, y)
    return g


def pts(x):
    """Points as an (m, n) array: a vector (one dimension) becomes a column."""
    x = np.asarray(x, float)
    return x[:, None] if x.ndim == 1 else x


# --------------------------------------------------------------------------------------------------------------------------------
# readouts: w_i = W(p_i), p_i ∈ ℝ^q, and the group learning rates
# --------------------------------------------------------------------------------------------------------------------------------
class Readout:
    """A unit-wise factorisation w = W(p), p ∈ ℝ^q, with the learning rate `lr[l]` of the l-th parameter. `w(p)`: (d,); `dw(p)`: ∂w/∂p
    (d, q); `init(d, eps, rng, rand)`: p(0) (d, q) with s² = ε b_i; `scale(eps)`: the wall's height in η (None: the readout has no wall)."""
    name = '?'
    q = 1

    def w(self, p):
        raise NotImplementedError

    def dw(self, p):
        raise NotImplementedError

    def masses(self, d, eps, rng=None, rand=0.0):
        """ε b_i with b_i = exp(rand·N(0, 1)) (rand = 0: equal masses)."""
        b = np.ones(d) if (rng is None or rand == 0.0) else np.exp(rand * rng.standard_normal(d))
        return eps * b

    def init(self, d, eps, rng=None, rand=0.0):
        raise NotImplementedError

    def scale(self, eps):
        raise NotImplementedError


class Linear(Readout):
    """D1, depth 1: w = p, started at w = 0 (or at ε b_i·sign with rand > 0). No factorisation, no wall: the control."""
    name = 'D1'
    q = 1

    def __init__(self, lr=1.0):
        self.lr = np.array([float(lr)])

    def w(self, p):
        return p[:, 0]

    def dw(self, p):
        return np.ones((p.shape[0], 1))

    def init(self, d, eps, rng=None, rand=0.0):
        if rng is None or rand == 0.0:
            return np.zeros((d, 1))
        return (self.masses(d, eps, rng, rand) * rng.choice([-1.0, 1.0], d))[:, None]

    def scale(self, eps):
        return None


class Product2(Readout):
    """D2, depth 2: w = u² − v², (u, v)(0) = (s, s) with s² = ε b_i (so w(0) = 0 and uv = s² is conserved: w = 2s² sinh(4·lr·η),
    η̇ = ⟨rφ⟩). The wall (|w| ≈ 1) is at η = ln(1/ε)/(4 lr): a LOGARITHMIC scale; an off-support coefficient at the margin μ of the wall
    is ε^μ, not ε^{3/2}μ^{−3}."""
    name = 'D2'
    q = 2

    def __init__(self, lr=1.0):
        self.lr = np.array([float(lr), float(lr)])

    def w(self, p):
        return p[:, 0] ** 2 - p[:, 1] ** 2

    def dw(self, p):
        return np.stack([2 * p[:, 0], -2 * p[:, 1]], 1)

    def init(self, d, eps, rng=None, rand=0.0):
        s = np.sqrt(self.masses(d, eps, rng, rand))
        return np.stack([s, s], 1)

    def scale(self, eps):
        return np.log(1.0 / eps) / (4 * self.lr[0])


class Cubic3(Readout):
    """D3, depth 3: w = u³ − v³, (u, v)(0) = (s, s) with s² = ε b_i (1/u + 1/v = 2/s conserved: w = (1/s − 3·lr·η)⁻³ − (1/s + 3·lr·η)⁻³,
    η̇ = ⟨rφ⟩). The wall is at η = 1/(3 lr s) ∝ ε^{−1/2}: a POWER scale; the coefficient at the margin μ is ε^{3/2}μ^{−3}/... — the
    paper's unit's structure (Paper 17, Proposition 3.1: an effective depth-3 unit w = u³ on the critical curve)."""
    name = 'D3'
    q = 2

    def __init__(self, lr=1.0):
        self.lr = np.array([float(lr), float(lr)])

    def w(self, p):
        return p[:, 0] ** 3 - p[:, 1] ** 3

    def dw(self, p):
        return np.stack([3 * p[:, 0] ** 2, -3 * p[:, 1] ** 2], 1)

    def init(self, d, eps, rng=None, rand=0.0):
        s = np.sqrt(self.masses(d, eps, rng, rand))
        return np.stack([s, s], 1)

    def scale(self, eps):
        return 1.0 / (3 * self.lr[0] * np.sqrt(eps))


class PowerN(Readout):
    """D_N, depth N ≥ 3: w = u^N − v^N, (u, v)(0) = (s, s), s² = ε b_i. The conserved quantity is u^{−(N−2)} + v^{−(N−2)} = 2s^{−(N−2)} and the
    wall is at η = s^{−(N−2)}/(N(N−2)·lr) ∝ ε^{−(N−2)/2}: a power scale, steeper the deeper the readout (N = 3 is `Cubic3`)."""

    def __init__(self, N, lr=1.0):
        assert N >= 3
        self.N, self.q, self.name = int(N), 2, f'D{int(N)}'
        self.lr = np.array([float(lr), float(lr)])

    def w(self, p):
        return p[:, 0] ** self.N - p[:, 1] ** self.N

    def dw(self, p):
        return np.stack([self.N * p[:, 0] ** (self.N - 1), -self.N * p[:, 1] ** (self.N - 1)], 1)

    def init(self, d, eps, rng=None, rand=0.0):
        s = np.sqrt(self.masses(d, eps, rng, rand))
        return np.stack([s, s], 1)

    def scale(self, eps):
        return eps ** (-(self.N - 2) / 2.0) / (self.N * (self.N - 2) * self.lr[0])


class PaperLaw(Readout):
    """P0, the paper's free mass–amplitude unit in raw coordinates: p = (α, ν), w = ν² a_max tanh(α/a_max)/m, the learning rates
    (Gm, γm/4) of `FreeLaw` L (and Gm for the centre), α(0) = amp0·N(0, 1) (0 = the paper's start), ν(0)² = ε b_i."""
    name = 'P0'
    q = 2

    def __init__(self, L, amp0=0.0):
        self.L, self.amp0 = L, float(amp0)
        self.lr = np.array([L.G * L.m, L.gamma * L.m / 4])

    def w(self, p):
        L = self.L
        return p[:, 1] ** 2 * L.A * np.tanh(p[:, 0] / L.A) / L.m

    def dw(self, p):
        L = self.L
        a = L.A * np.tanh(p[:, 0] / L.A)
        return np.stack([p[:, 1] ** 2 * (1.0 - (a / L.A) ** 2) / L.m, 2 * p[:, 1] * a / L.m], 1)

    def init(self, d, eps, rng=None, rand=0.0):
        al = np.zeros(d) if (rng is None or self.amp0 == 0.0) else self.amp0 * rng.standard_normal(d)
        return np.stack([al, np.sqrt(self.masses(d, eps, rng, rand))], 1)

    def scale(self, eps):
        return transit_scale(self.L, eps)


class Unsaturated(PaperLaw):
    """P0u: the paper's unit without the amplitude's saturation, w = ν²α/m (the rates and the start as `PaperLaw`). Its conserved
    quantity is ν² − (γ/(2G))α², and at small mass it has the same core and the same wall as P0."""
    name = 'P0u'

    def w(self, p):
        return p[:, 1] ** 2 * p[:, 0] / self.L.m

    def dw(self, p):
        return np.stack([p[:, 1] ** 2 / self.L.m, 2 * p[:, 1] * p[:, 0] / self.L.m], 1)


class PaperPairs(PaperLaw):
    """P0p, the paper's unit started on branches (Paper 17, Proposition 'starts on branches'): the units come in PAIRS on the same atom — the
    first half of the units has α(0) = +a₀, the second half α(0) = −a₀ (so d = 2·(atoms) and the units j and j + d/2, started at the same
    centre, are the pair of atom j; `atoms(d)` maps units to atoms), all with ν(0)² = ε b_i. A kernel started at α(0) ≠ 0 sits on a branch
    of its curve (β small, a ≈ a(0)), where its mass grows as β = ε·exp(γ a(0)(η − η(0))): the wall is reached after ln(1/ε)/(γ|a(0)|), the paper's
    LOGARITHMIC scale L_ε = ln(1/ε), weight ρ = 1/(γ|a(0)|) and τ = t/L_ε. A pair of equal magnitude gives the two-sided box |u| ≤ ρ, so the
    limit is the path of equal weights, and `scale(ε)` = L_ε ρ = ln(1/ε)/(γ|a(0)|), a(0) = a_max tanh(a₀/a_max), is the time at which that
    path (unit weights) has τ = 1. The formed threshold is crossed before the wall: the path's FIRST entrant is early by exactly the factor
    1 − c/L_ε of `form_ratio`, to leading order in 1/L_ε; a later atom is early by it only while its dual is straight until its wall (separated
    features), and by more or less when an earlier entry has bent the dual (part K's coupled targets). The proposition is for fixed centres; the
    sliding extension that part K tests is not proved."""
    name = 'P0p'

    def __init__(self, L, a0):
        assert a0 > 0
        super().__init__(L)
        self.a0 = float(a0)

    @staticmethod
    def atoms(d):
        """The atom of each of the d (even) units: unit j and unit j + d/2 are the pair of atom j."""
        assert d % 2 == 0
        return np.arange(d) % (d // 2)

    def init(self, d, eps, rng=None, rand=0.0):
        assert d % 2 == 0
        al = np.where(np.arange(d) < d // 2, self.a0, -self.a0)
        return np.stack([al, np.sqrt(self.masses(d, eps, rng, rand))], 1)

    def scale(self, eps):
        return np.log(1.0 / eps) / (self.L.gamma * self.L.A * np.tanh(self.a0 / self.L.A))

    def form_ratio(self, eps, thr):
        """The predicted ratio of the time at which a kernel of the pair reaches |w| = thr to the wall's time ρL_ε, for a kernel whose drive is
        constant until then: the path's FIRST entrant. Along the kernel's curve dη = dβ/(γaβ) (Lemma 2.1(v)) with β = ε + 2v_s(sinh²(α/a_max) −
        sinh²(a₀/a_max)), so the mirror coordinate has grown by (1/γ)∫ dβ/(aβ) from β = ε to the β_θ at which β_θ a/m = thr, against the wall's
        η = L_ε/(γ|a(0)|); the ratio is a(0)·∫ d ln β / a over L_ε = 1 − c/L_ε. The relation between w and η is the kernel's own, whatever the
        drive (Lemma 2.1(v)), so the formed level w = thr is crossed when the kernel's dual (η over its wall) reaches this ratio κ = 1 − c/L_ε, not
        at the wall (dual 1). The first entrant's dual is a straight line from τ = 0 (the residual is f*, the dormant centres do not move), so it
        reaches κ at the fraction κ of its entry time: it is early by exactly c/L_ε, to leading order in 1/L_ε (the entered kernel's growth beyond
        the threshold, to its O(1) amplitude, is the next order).

        This is NOT the earliness of every entry. "A path whose walls are all scaled by one factor is the same path in a scaled time" scales the entry
        levels (the walls), but the formed level is not a wall: a later atom is early by this factor only while its dual is straight until its wall
        (separated features). An earlier entry changes the dual's slope, κ is then reached on another part of the dual, and the earliness is
        slope-dependent at the same order 1/L_ε. Part K's coupled targets (16 targets of the benchmark with two features 1.5σ or 2σ apart, ε = 10⁻⁶ … 10⁻⁴⁸):
        the first entrant is early by 1 − ratio on every one (to 0.004), the second atom of an opposite-sign pair too (to 0.006: its dual is straight
        to its wall), the second atom of a same-sign pair by 1.4 to 6 times that (0.62–0.83 against 0.448 at ε = 10⁻⁶: its dual rises on its initial
        slope to about 0.9 of its wall by the first entry, then slows down, and the path enters it 1.4 to 3.4 times later than that straight line
        would have: the formed level is crossed on the early part of the dual; formed is not entered)."""
        L, s0 = self.L, np.sinh(self.a0 / self.L.A) ** 2

        def a_of(beta):
            s2 = s0 + (beta - eps) / (2 * L.c)
            return L.A * np.sqrt(s2 / (1.0 + s2))
        beta_th = brentq(lambda b: b * a_of(b) / L.m - thr, eps, 1e8)
        integral = quad(lambda lb: 1.0 / a_of(np.exp(lb)), np.log(eps), np.log(beta_th), limit=400)[0]
        return float(L.A * np.tanh(self.a0 / L.A) * integral / np.log(1.0 / eps))


# --------------------------------------------------------------------------------------------------------------------------------
# the gradient flow
# --------------------------------------------------------------------------------------------------------------------------------
def make_rhs(readout, K, fstar, lr_c, d, n):
    """The plain gradient flow z = (p_1..p_d, c_1..c_d), ṗ = lr·g_i·∂w/∂p, ċ = lr_c·w_i·h_i with g_i = ⟨rφ_i⟩, h_i = ⟨r∇_cφ_i⟩,
    r = f* − Σ w_iφ_i: the negative gradient of F = ½⟨r²⟩ times the group learning rates. K: a dictionary with `synth(c, w)` and
    `corr(c, r)` (`rich_transit_nd.BumpsND`, `BumpsGrid2`); centres (d, n)."""
    q = readout.q
    lr = readout.lr[None, :]

    def f(t, z):
        p = z[:d * q].reshape(d, q)
        c = z[d * q:].reshape(d, n)
        w = readout.w(p)
        r = fstar - K.synth(c, w)
        g, h = K.corr(c, r)
        return np.r_[(lr * g[:, None] * readout.dw(p)).ravel(), (lr_c * w[:, None] * h).ravel()]
    return f


def fit_value(K, fstar, w, c):
    """F = ½⟨(f* − Σ_i w_iφ(·; c_i))²⟩ and the rms residual."""
    r = fstar - K.synth(pts(c), w)
    return 0.5 * float(np.mean(r ** 2)), float(np.sqrt(np.mean(r ** 2)))


def _unpack(readout, Z, d, n):
    """The record Z (T, d·q + d·n) → W (T, d), C (T, d, n)."""
    q = readout.q
    P = Z[:, :d * q].reshape(-1, d, q)
    C = Z[:, d * q:].reshape(-1, d, n)
    W = np.stack([readout.w(P[k]) for k in range(P.shape[0])])
    return W, C


def train_flow(readout, K, fstar, c0, eps, tau_end, lr_c, n_t=1201, rtol=1e-9, atol_p=1e-10, atol_c=1e-12, method='LSODA',
               rng=None, rand=0.0, t_scale=None, n_seg=8, max_secs=None):
    """Integrate the gradient flow of `make_rhs` from the readout's start over τ ∈ [0, τ_end], τ = t/scale (scale: the readout's wall,
    `t_scale` overrides), recording on n_t uniform τ. Done in n_seg pieces so that a wall-clock budget `max_secs` can end it between
    pieces (status 'timeout'; the record then holds the pieces done, padded with the last state). Per-component absolute tolerances:
    atol_p·max(|p(0)|, √ε), atol_c for the centres. Returns dict(taus, W (n_t, d), C (n_t, d, n), status, nfev, secs,
    scale, p_end, radau_pieces, aborted_nfev): the pieces redone by Radau and the evaluations of the abandoned LSODA attempts (nfev counts the
    accepted solves only)."""
    c0 = pts(c0)
    d, n = c0.shape
    q = readout.q
    scale = readout.scale(eps) if t_scale is None else float(t_scale)
    if scale is None:
        raise ValueError('this readout has no wall: pass t_scale')
    p0 = readout.init(d, eps, rng, rand)
    z = np.r_[p0.ravel(), c0.ravel()]
    # absolute tolerances: atol_p times the parameter's size, but at least times √ε (α(0) = 0 would give a vanishing tolerance)
    atol = np.r_[atol_p * np.maximum(np.abs(p0).ravel(), np.sqrt(eps)), np.full(d * n, atol_c)]
    f = make_rhs(readout, K, fstar, lr_c, d, n)
    taus = np.linspace(0.0, tau_end, n_t)
    Z = np.empty((n_t, z.size))
    Z[0] = z
    cuts = np.unique(np.linspace(0, n_t - 1, n_seg + 1).astype(int))
    t0, nfev, status = time.time(), 0, 'ok'
    radau, aborted = 0, 0                    # pieces redone by Radau; right-hand-side evaluations of the LSODA attempts that were abandoned
    for a, b in zip(cuts[:-1], cuts[1:]):
        if max_secs is not None and time.time() - t0 > max_secs:
            Z[a + 1:] = Z[a]
            status = 'timeout'
            break
        for meth in ((method, 'Radau') if method == 'LSODA' else (method,)):
            try:
                sol = solve_ivp(_counted(f, STIFF_BUDGET) if meth == 'LSODA' else f, (scale * taus[a], scale * taus[b]), Z[a], method=meth,
                                rtol=rtol, atol=atol, t_eval=scale * taus[a + 1:b + 1])
            except _OverBudget:                        # LSODA stalled (stiffness not detected): the piece is redone by Radau
                aborted += STIFF_BUDGET
                radau += 1
                continue
            Y = np.asarray(sol.y)
            if sol.status >= 0 and Y.ndim == 2 and Y.shape[1] == b - a:
                nfev += sol.nfev
                break                                  # (also LSODA's "unexpected istate" on a relaxed state: redone by Radau)
            if meth == 'LSODA' and method == 'LSODA':
                aborted += sol.nfev
                radau += 1
            else:
                nfev += sol.nfev
        if sol.status < 0 or Y.ndim != 2 or Y.shape[1] != b - a:
            k = Y.shape[1] if Y.ndim == 2 else 0
            if k > 0:                                  # (an empty result — LSODA's "Unexpected istate" — has nothing to store)
                Z[a + 1:a + 1 + k] = Y.T[:k]
            Z[a + 1 + k:] = Z[a + k]
            status = 'failed: ' + str(sol.message)[:60]
            break
        Z[a + 1:b + 1] = Y.T
    W, C = _unpack(readout, Z, d, n)
    return dict(taus=taus, W=W, C=C, status=status, nfev=nfev, secs=time.time() - t0, scale=scale, eps=eps,
                p_end=Z[-1, :d * q].reshape(d, q), readout=readout.name, radau_pieces=radau, aborted_nfev=aborted)


def train_gd(readout, K, fstar, c0, eps, tau_end, lr_c, h, n_t=1201, rng=None, rand=0.0, t_scale=None, max_secs=None, blow=1e9):
    """Plain gradient DESCENT with the fixed step h (explicit Euler on `make_rhs`: every parameter moves by h·lr·(−∂F/∂p)) over the
    time τ_end·scale, recorded on n_t uniform τ (the nearest steps). status 'diverged' if a parameter exceeds `blow` or is not finite;
    'timeout' on the wall-clock budget. The number of steps is ⌈τ_end·scale/h⌉."""
    c0 = pts(c0)
    d, n = c0.shape
    q = readout.q
    scale = readout.scale(eps) if t_scale is None else float(t_scale)
    p0 = readout.init(d, eps, rng, rand)
    z = np.r_[p0.ravel(), c0.ravel()]
    f = make_rhs(readout, K, fstar, lr_c, d, n)
    n_steps = int(np.ceil(tau_end * scale / h))
    rec = np.unique(np.linspace(0, n_steps, n_t).round().astype(np.int64))
    Z = np.empty((rec.size, z.size))
    Z[0] = z
    t0, status, k_rec = time.time(), 'ok', 1
    for step in range(1, n_steps + 1):
        z = z + h * f(0.0, z)
        if step == rec[k_rec]:
            Z[k_rec] = z
            k_rec += 1
            if not np.all(np.isfinite(z)) or np.abs(z).max() > blow:
                status = 'diverged'
                break
            if max_secs is not None and time.time() - t0 > max_secs:
                status = 'timeout'
                break
        if k_rec >= rec.size:
            break
    if k_rec < rec.size:
        Z[k_rec:] = Z[k_rec - 1]
    taus_rec = rec / n_steps * tau_end
    W, C = _unpack(readout, Z, d, n)
    taus = np.linspace(0.0, tau_end, n_t)
    Wg = np.stack([np.interp(taus, taus_rec, W[:, i]) for i in range(d)], 1)
    Cg = np.stack([np.stack([np.interp(taus, taus_rec, C[:, i, k]) for k in range(n)], -1) for i in range(d)], 1)
    return dict(taus=taus, W=Wg, C=Cg, status=status, nfev=n_steps, secs=time.time() - t0, scale=scale, eps=eps, readout=readout.name,
                h=h, steps=n_steps)


def train_adam(readout, K, fstar, c0, eps, n_steps, lr_adam, n_t=1201, betas=(0.9, 0.999), eps_adam=1e-8, rng=None, rand=0.0,
               max_secs=None, blow=1e9):
    """Adam (one step size `lr_adam` for every parameter and centre; the gradient is that of F) from the readout's start, n_steps steps,
    recorded on n_t uniform steps. Time is the step count: τ in the record is step/n_steps. Reported, not claimed: Adam's per-parameter
    normalisation rescales the dormant parameters' gradients."""
    c0 = pts(c0)
    d, n = c0.shape
    q = readout.q
    p0 = readout.init(d, eps, rng, rand)
    z = np.r_[p0.ravel(), c0.ravel()]
    ones = np.ones_like(readout.lr)
    base = make_rhs(Unit(readout, ones), K, fstar, 1.0, d, n)         # the direction −∂F/∂θ (unit learning rates)
    m, v = np.zeros_like(z), np.zeros_like(z)
    rec = np.unique(np.linspace(0, n_steps, n_t).round().astype(np.int64))
    Z = np.empty((rec.size, z.size))
    Z[0] = z
    t0, status, k_rec = time.time(), 'ok', 1
    for step in range(1, n_steps + 1):
        gneg = base(0.0, z)
        m = betas[0] * m + (1 - betas[0]) * gneg
        v = betas[1] * v + (1 - betas[1]) * gneg ** 2
        z = z + lr_adam * (m / (1 - betas[0] ** step)) / (np.sqrt(v / (1 - betas[1] ** step)) + eps_adam)
        if step == rec[k_rec]:
            Z[k_rec] = z
            k_rec += 1
            if not np.all(np.isfinite(z)) or np.abs(z).max() > blow:
                status = 'diverged'
                break
            if max_secs is not None and time.time() - t0 > max_secs:
                status = 'timeout'
                break
        if k_rec >= rec.size:
            break
    if k_rec < rec.size:
        Z[k_rec:] = Z[k_rec - 1]
    W, C = _unpack(readout, Z, d, n)
    taus_rec = rec / n_steps
    taus = np.linspace(0.0, 1.0, n_t)
    Wg = np.stack([np.interp(taus, taus_rec, W[:, i]) for i in range(d)], 1)
    Cg = np.stack([np.stack([np.interp(taus, taus_rec, C[:, i, k]) for k in range(n)], -1) for i in range(d)], 1)
    return dict(taus=taus, W=Wg, C=Cg, status=status, nfev=n_steps, secs=time.time() - t0, scale=float(n_steps), eps=eps,
                readout=readout.name, lr_adam=lr_adam)


class Unit(Readout):
    """The readout `r` with its learning rates replaced (Adam's unit-rate direction)."""

    def __init__(self, r, lr):
        self._r, self.lr, self.q, self.name = r, np.asarray(lr, float), r.q, r.name

    def w(self, p):
        return self._r.w(p)

    def dw(self, p):
        return self._r.dw(p)


# --------------------------------------------------------------------------------------------------------------------------------
# comparing a trained network with the sliding path
# --------------------------------------------------------------------------------------------------------------------------------
def formation_times(taus, W, thr):
    """For every unit the first τ at which |w| ≥ thr (log-linear interpolation between the record's points; nan if never)."""
    A = np.abs(W)
    out = np.full(A.shape[1], np.nan)
    for i in range(A.shape[1]):
        above = np.flatnonzero(A[:, i] >= thr)
        if above.size == 0:
            continue
        k = above[0]
        if k == 0:
            out[i] = taus[0]
            continue
        a0, a1 = max(A[k - 1, i], 1e-300), A[k, i]
        fr = (np.log(thr) - np.log(a0)) / (np.log(a1) - np.log(a0)) if a1 > a0 else 1.0
        out[i] = taus[k - 1] + fr * (taus[k] - taus[k - 1])
    return out


def match_features(c_units, w_units, cstar, astar):
    """Each feature k to the nearest of the given units (centres (m, n), coefficients (m,)): (distinct, dist (K,), rel amp err (K,)).
    distinct: whether the nearest units are K different ones."""
    cstar, c_units = pts(cstar), pts(c_units)
    nearest, dist, aerr = [], [], []
    for k in range(cstar.shape[0]):
        dd = np.linalg.norm(c_units - cstar[k][None, :], axis=1)
        j = int(np.argmin(dd))
        nearest.append(j)
        dist.append(float(dd[j]))
        aerr.append(float(abs(w_units[j] - astar[k]) / abs(astar[k])))
    return len(set(nearest)) == len(nearest), np.array(dist), np.array(aerr)


def analyse(res, K, fstar, c0, cstar, astar, ref=None, thr_rel=THETA, tol_pos=0.05, tol_amp=0.05):
    """The trained network's record `res` against the target (features cstar, astar) and, with `ref`, against the sliding path
    (dict(order = the entrants in order, tau = their times, set = the final active set)). Returns dict:
      n_formed        number of units with |w| > thr = min(thr_rel·max|A*|, 0.1·min|A*|) at the end; also swept in `n_formed_sweep`
                      (1e-4 … 1e-1 times max|A*|);
      gap             log10(smallest formed |w| / largest unformed |w|) at the end (nan if none formed or none unformed);
      dormant_disp    largest displacement |c_i(end) − c_i(0)| over the units not formed at the end, in units of the kernel width;
      dormant_wmax    the largest |w_i| over all times of those units;
      gap_k, exact_topk  threshold-free: log₁₀ of the K-th over the (K+1)-th largest |w|, and whether the K largest units are exact;
      fit_rms         the rms residual at the end;
      exact           K units formed, one per feature, each within tol_pos·σ and tol_amp (relative) of its feature (pos_err, amp_err);
      set_equal       the formed set is the path's; order_equal: entrants in the path's order; tau_rel: the largest relative difference
                      of the formation times from the path's entry times (the units of the path's set)."""
    # a unit has formed when |w| > thr_rel·max|A*| (Paper 17's convention) — and at most 0.1·min|A*|, so that a feature 10³ times smaller than
    # the largest still counts (the two agree whenever min|A*| ≥ 0.1·thr_rel·max|A*|·10, i.e. on every target of amplitudes within a factor 10)
    thr = min(thr_rel * float(np.max(np.abs(astar))), 0.1 * float(np.min(np.abs(astar))))
    W, C, taus = res['W'], res['C'], res['taus']
    d = W.shape[1]
    c0 = pts(c0)
    formed = np.abs(W[-1]) > thr
    out = dict(n_formed=int(formed.sum()),
               n_formed_sweep={s: int((np.abs(W[-1]) > s * float(np.max(np.abs(astar)))).sum()) for s in (1e-4, 1e-3, 1e-2, 1e-1)})
    aw = np.abs(W[-1])
    out['gap'] = float(np.log10(aw[formed].min() / max(aw[~formed].max(), 1e-300))) if formed.any() and (~formed).any() else float('nan')
    sig = float(K.sig)
    disp = np.linalg.norm(C[-1] - c0, axis=1) / sig
    out['dormant_disp'] = float(disp[~formed].max()) if (~formed).any() else 0.0
    out['dormant_wmax'] = float(np.abs(W[:, ~formed]).max()) if (~formed).any() else 0.0
    out['fit_rms'] = fit_value(K, fstar, W[-1], C[-1])[1]
    # threshold-free: the K-th over the (K+1)-th largest |w| (decades), and whether the K largest units are exact (distinct, within the tolerances)
    srt = np.sort(aw)[::-1]
    K_ = len(astar)
    out['gap_k'] = float(np.log10(srt[K_ - 1] / max(srt[K_], 1e-300))) if d > K_ else float('nan')
    top = np.argsort(aw)[::-1][:K_]
    dist_ok, dist_t, aerr_t = match_features(C[-1][top], W[-1][top], cstar, astar)
    out['exact_topk'] = bool(dist_ok and dist_t.max() / sig <= tol_pos and aerr_t.max() <= tol_amp)
    ok = False
    pos_err = amp_err = float('nan')
    if formed.sum() == len(astar):
        idx = np.flatnonzero(formed)
        distinct, dist, aerr = match_features(C[-1][idx], W[-1][idx], cstar, astar)
        pos_err, amp_err = float(dist.max() / sig), float(aerr.max())
        ok = bool(distinct and pos_err <= tol_pos and amp_err <= tol_amp)
    out.update(exact=ok, pos_err=pos_err, amp_err=amp_err)
    ft = formation_times(taus, W, thr)
    out['form_order'] = [int(j) for j in np.argsort(np.where(np.isnan(ft), np.inf, ft)) if not np.isnan(ft[j])]
    if ref is not None:
        rset = sorted(int(j) for j in ref['set'])
        out['set_equal'] = sorted(int(j) for j in np.flatnonzero(formed)) == rset
        mine = [j for j in out['form_order'] if j in rset]
        out['order_equal'] = mine == [int(j) for j in ref['order']]
        rt = np.array([ref['tau'][list(ref['order']).index(j)] for j in rset])
        mt = ft[rset]
        out['tau_rel'] = float(np.max(np.abs(mt - rt) / rt)) if not np.any(np.isnan(mt)) else float('nan')
    return out


# --------------------------------------------------------------------------------------------------------------------------------
# units grouped by atom (the branch pairs of `PaperPairs`)
# --------------------------------------------------------------------------------------------------------------------------------
def formed_threshold(astar, thr_rel=THETA):
    """The formed threshold of `analyse`: min(thr_rel·max|A*|, 0.1·min|A*|)."""
    a = np.abs(np.asarray(astar, float))
    return min(thr_rel * float(a.max()), 0.1 * float(a.min()))


def atom_summary(res, astar, atom_of, ref=None, thr_rel=THETA):
    """The record `res` of a network whose units sit on atoms (`atom_of[i]`: the atom of unit i; `PaperPairs.atoms`), against the sliding path
    `ref` of the ATOMS' dictionary (dict(order, tau, set) as in `analyse`; its indices are the atoms'). An atom has *formed* when one of its
    units has |w| > the formed threshold at the end (`formed_threshold`), and its formation time is the earlier of its units'. Returns
      n_atoms, atoms   the formed atoms (their number and sorted indices);
      one_each         no atom has two formed units (with n_atoms = K, the K formed units are one per atom);
      form_order       the atoms in the order of their formation times (the atoms that ever crossed the threshold);
    and with ref
      set_equal        the formed atoms are the path's set;  order_equal: the path's atoms form in the path's order (other atoms ignored);
      tau_rel          the largest relative difference of the atoms' formation times from the path's entry times (nan if one of them
                       never formed);
      tau_ratio        the same for the ratios t_j/t_first to the path's τ_j/τ_first (first: the path's first entrant): the entry times
                       up to a common factor (nan if one never formed);
      form_times, earliness
                       the atoms of the path's set in the order of the path's entry times: their formation times t_j and
                       1 − t_j/τ_j, how early each forms against its entry (nan for an atom that never formed; one atom's nan does not
                       affect the others')."""
    thr = formed_threshold(astar, thr_rel)
    W, taus = res['W'], res['taus']
    atom_of = np.asarray(atom_of, int)
    na = int(atom_of.max()) + 1
    formed = np.abs(W[-1]) > thr
    per_atom = np.bincount(atom_of[formed], minlength=na)
    out = dict(n_atoms=int((per_atom > 0).sum()), atoms=[int(j) for j in np.flatnonzero(per_atom > 0)], one_each=bool((per_atom <= 1).all()))
    ft_u = formation_times(taus, W, thr)
    ft = np.full(na, np.nan)
    for j in range(na):
        v = ft_u[atom_of == j]
        v = v[~np.isnan(v)]
        if v.size:
            ft[j] = v.min()
    out['form_order'] = [int(j) for j in np.argsort(np.where(np.isnan(ft), np.inf, ft)) if not np.isnan(ft[j])]
    if ref is not None:
        rset = sorted(int(j) for j in ref['set'])
        out['set_equal'] = out['atoms'] == rset
        out['order_equal'] = [j for j in out['form_order'] if j in rset] == [int(j) for j in ref['order']]
        rt = np.array([ref['tau'][list(ref['order']).index(j)] for j in rset])
        mt = ft[rset]
        if np.any(np.isnan(mt)):
            out['tau_rel'] = out['tau_ratio'] = float('nan')
        else:
            out['tau_rel'] = float(np.max(np.abs(mt - rt) / rt))
            k0 = int(np.argmin(rt))
            out['tau_ratio'] = float(np.max(np.abs((mt / mt[k0]) / (rt / rt[k0]) - 1.0)))
        by_entry = np.argsort(rt, kind='stable')                                # the path's first entrant first
        out['form_times'] = [float(x) for x in mt[by_entry]]
        out['earliness'] = [float(1.0 - x) for x in (mt / rt)[by_entry]]
    return out
