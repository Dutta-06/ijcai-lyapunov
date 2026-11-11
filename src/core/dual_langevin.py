"""Proposition D‴.50: the free law under label noise — the dual of the loop is a Langevin diffusion.

On a dictionary (Proposition D‴.45: Φ ∈ ℝ^{N×d} of full row rank, η = η(0) + Φᵀζ) an isotropic white label noise of intensity σ, read
in Stratonovich's sense or run through the mirror update, gives

    dζ = (y − Φw)dt/N + (σ/N)dW = −∇F(ζ)dt/N + (2T/N)^{1/2}dW,   F(ζ) = Ψ*(η(0) + Φᵀζ) − yᵀζ,   T = σ²/(2N),

F the objective of the Bregman projection's dual (∇F = Φw − y), so the law of ζ converges to π ∝ exp(−F/T).

- `FastFree`: η, dα/dη, the inverse (warm-started Newton), w, dw/dη and ψ* on free curves (q > −c), vectorised;
- `mirror_update`: the loop's dual Euler–Maruyama under label noise (the mirror update; the unadjusted Langevin algorithm on F);
- `primal_sde`: the paradigm's law on (α, log β) under label noise, Stratonovich (stochastic Heun) or Itô (Euler–Maruyama);
- `gibbs_rwm`: a random-walk Metropolis sampler of π — an independent check of the Gibbs law;
- `fd_covariance`: the low-noise covariance of w, T·D^{1/2}PD^{1/2} (D = diag dw/dη, P the projector onto range(D^{1/2}Φᵀ));
- `ito_q_drift`: the drift of the conserved quantities when the noise is read in Itô's sense;
- `tangent_projectors`, `linear_response`: the transit's noise to first order along an inverse scale space path (the reflected
  dynamics linearised: δν ← P_T(δν + κdW)).
"""
import numpy as np


class FastFree:
    """The free law's mirror coordinate on the curves Q_i = q_i > −c (from α = 0, η odd in α), vectorised over the particles."""

    def __init__(self, law, q):
        q = np.asarray(q, float)
        if np.any(q <= -law.c):
            raise ValueError('FastFree: every curve must be free (q > −c)')
        self.L, self.q = law, q
        self.low = q <= law.c                                                  # κ′² = (c − q)/(q + c) ≥ 0
        self.kp = np.sqrt(np.where(self.low, (law.c - q) / (q + law.c), 0.0))
        self.k = np.sqrt(np.where(self.low, 0.0, (q - law.c) / (q + law.c)))  # κ² = (q − c)/(q + c) ∈ (0, 1)
        self.dl = np.where(self.low, 0.0, (2 * law.c / (q + law.c)) / (1 + self.k))   # δ = 1 − κ, without cancellation
        with np.errstate(divide='ignore'):                                     # η ≥ (α − a_maxκ artanh κ)/(2Gc) where q > c
            self.xcap = np.where(self.low, np.inf, law.A * self.k * 0.5 * np.log((1 + self.k) / np.where(self.low, 1.0, self.dl)))
        self.hi = np.flatnonzero(~self.low)

    def eta(self, x):
        """η at α = x ≥ 0 (FreeLaw.eta's closed forms): (x + a_maxκ′arctan(κ′t))/(2Gc) where q ≤ c; where q > c,
        (a_max/2)[log1p(δt/(1 + κt)) + log1p(δt/(1 − t)) + δ log((1 + κt)/(1 − κt))]/(2Gc), t = tanh(x/a_max), every term without
        cancellation."""
        L = self.L
        x = np.asarray(x, float)
        t = np.tanh(x / L.A)
        out = (x + L.A * self.kp * np.arctan(self.kp * t)) / (2 * L.G * L.c)
        if self.hi.size:                                                      # the particles index the last axis
            h = self.hi
            k, dl, th = self.k[h], self.dl[h], t[..., h]
            ex = np.exp(-2 * x[..., h] / L.A)
            u = np.maximum(2 * ex / (1 + ex), 1e-300)                         # 1 − t
            l1m = np.where(th < 0.5, np.log1p(-k * np.minimum(th, 0.5)), np.log(dl + k * u))
            v = np.log1p(dl * th / (1 + k * th)) + np.log1p(dl * th / u) + dl * (np.log1p(k * th) - l1m)
            out[..., h] = 0.5 * L.A * v / (2 * L.G * L.c)
        return out

    def dal(self, x):
        """dα/dη = Gβa′ = G(2c + (q − c)sech²(α/a_max)) — no underflow at large α."""
        L = self.L
        return L.G * (2 * L.c + (self.q - L.c) / np.cosh(np.minimum(np.abs(x) / L.A, 350.0)) ** 2)

    def inv(self, e, x0=None, it_max=80, tol=1e-14):
        """α with η(α) = e from the warm start x0: Newton on |α|. Where q < c η is concave on α ≥ 0: a step past 0 is clamped to 0, from
        where Newton is monotone; where q > c it is convex: the iterates are clamped below the root's bound 2Gc|e| + a_maxκ artanh κ, from
        where Newton is monotone. Stopped at the rounding level of η."""
        e = np.asarray(e, float)
        t = np.abs(e)
        up = 2 * self.L.G * self.L.c * t + self.xcap
        x = np.zeros_like(t) if x0 is None else np.abs(np.asarray(x0, float)).copy()
        x = np.minimum(x, up)
        floor = 4 * np.spacing(np.maximum(t, 1.0))
        for _ in range(it_max):
            g = self.dal(x)
            xn = np.clip(x - (self.eta(x) - t) * g, 0.0, up)
            if np.all(np.abs(xn - x) <= tol * (1 + xn) + floor * g):
                x = xn
                break
            x = xn
        else:
            if np.any(np.abs(self.eta(x) - t) > 1e-10 * (1 + t)):
                raise RuntimeError('FastFree.inv: no convergence')
        return np.copysign(x, e)

    def w(self, al):
        return self.L.w(al, self.q)

    def dw(self, al):
        """dw/dη = (β/m)(Gβa′² + γa²) = m·relevance."""
        return self.L.m * self.L.relevance(al, self.q)

    def psi_star(self, al):
        """ψ*(η) = 2c sinh²(α/a_max)/(γm) (zero at α = 0; dψ*/dη = w)."""
        L = self.L
        return 2 * L.c * np.sinh(np.asarray(al, float) / L.A) ** 2 / (L.gamma * L.m)


def potential(F, Phi, y, zeta, eta0=None, al0=None):
    """(F(ζ), ∇F(ζ) = Φw − y, α) for the dual objective F(ζ) = Ψ*(η(0) + Φᵀζ) − yᵀζ."""
    eta = (0.0 if eta0 is None else eta0) + Phi.T @ zeta
    al = F.inv(eta, al0)
    w = F.w(al)
    return float(F.psi_star(al).sum() - y @ zeta), Phi @ w - y, al


def mirror_update(F, Phi, y, sig, h, t_end, seed=0, rec_every=1, zeta0=None, al0=None, eta0=None, noise=None, rec_al=False):
    """ζ ← ζ + h(y − Φw)/N + (σ/N)√h ξ, w = w(η(0) + Φᵀζ): the mirror update under isotropic label noise (each particle exactly on its
    curve) — the unadjusted Langevin algorithm on F, whose law is π + O(h) over practical horizons (the potential's gradient grows
    exponentially, so the scheme is not globally stable in the strict sense). Returns (t, ζ, w) at every rec_every-th step and the last α
    (with rec_al, the recorded α instead — for dw/dη along the run). `noise(k)` may supply the standard normals ξ_k (else a generator
    seeded by `seed`)."""
    N, d = Phi.shape
    rng = np.random.default_rng(seed)
    z = np.zeros(N) if zeta0 is None else np.array(zeta0, float)
    e0 = np.zeros(d) if eta0 is None else eta0
    al = np.zeros(d) if al0 is None else np.array(al0, float)
    nst = int(round(t_end / h))
    sq = (sig / N) * np.sqrt(h)
    ts, Z, W, AL = [], [], [], []
    for k in range(nst + 1):
        al = F.inv(e0 + Phi.T @ z, al)
        w = F.w(al)
        if k % rec_every == 0 or k == nst:
            ts.append(k * h); Z.append(z.copy()); W.append(w)
            if rec_al:
                AL.append(al.copy())
        if k < nst:
            xi = rng.standard_normal(N) if noise is None else noise(k)
            z = z + (y - Phi @ w) * (h / N) + sq * xi
    return np.array(ts), np.array(Z), np.array(W), (np.array(AL) if rec_al else al)


def primal_sde(L, Phi, y, al0, lb0, sig, h, t_end, seed=0, rec_every=1, scheme='stratonovich'):
    """The paradigm's law under label noise: dα = Gβa′dξ, d log β = γa dξ, dξ = Φᵀ((y − Φw)dt + σdW)/N. 'stratonovich': read in
    Stratonovich's sense, which keeps every Q_i (Proposition D‴.45(iv)(a)) — so α alone, dα = G(2c + (q − c)sech²(α/a_max))∘dξ with β on
    its curve, by the stochastic Heun scheme (strong order 1: the noise is commutative, every sample's field being Σ_iΦ_ki V_i(α_i));
    'ito': the update on (α, log β) read in Itô's sense (Euler–Maruyama), which drifts Q (`ito_q_drift`). Returns (t, α, log β) at
    every rec_every-th step."""
    N, d = Phi.shape
    rng = np.random.default_rng(seed)
    al, lb = np.array(al0, float), np.array(lb0, float)
    q = L.Q(al, np.exp(lb))
    if scheme == 'stratonovich':
        F = FastFree(L, q)

    def drive(al, b, dW):
        return Phi.T @ ((y - Phi @ (b * L.amp(al)[0] / L.m)) * h + sig * dW) / N

    nst = int(round(t_end / h))
    ts, A, B = [], [], []
    for k in range(nst + 1):
        if k % rec_every == 0 or k == nst:
            ts.append(k * h); A.append(al.copy()); B.append(lb.copy())
        if k == nst:
            break
        dW = np.sqrt(h) * rng.standard_normal(N)
        if scheme == 'stratonovich':
            v1 = F.dal(al) * drive(al, L.beta(al, q), dW)
            a1 = al + v1
            al = al + 0.5 * (v1 + F.dal(a1) * drive(a1, L.beta(a1, q), dW))
            lb = np.log(L.beta(al, q))
        elif scheme == 'ito':
            b = np.exp(lb)
            a, da = L.amp(al)
            xi = drive(al, b, dW)
            al, lb = al + L.G * b * da * xi, lb + L.gamma * a * xi
        else:
            raise ValueError(scheme)
    return np.array(ts), np.array(A), np.array(B)


def gibbs_rwm(F, Phi, y, T, zeta0, cov, n, seed=0, thin=1, eta0=None):
    """Random-walk Metropolis for π ∝ exp(−F/T) in ζ with the Gaussian proposal N(0, cov). Returns (samples of ζ and of w every thin-th
    step, the acceptance rate)."""
    rng = np.random.default_rng(seed)
    C = np.linalg.cholesky(cov)
    z = np.array(zeta0, float)
    f, g, al = potential(F, Phi, y, z, eta0)
    Zs, Ws, acc = [], [], 0
    for k in range(n):
        zp = z + C @ rng.standard_normal(z.size)
        try:
            fp, gp, alp = potential(F, Phi, y, zp, eta0, al)
        except (RuntimeError, FloatingPointError, OverflowError):
            fp = np.inf
        if np.log(rng.uniform()) < (f - fp) / T:
            z, f, al = zp, fp, alp
            acc += 1
        if k % thin == 0:
            Zs.append(z.copy()); Ws.append(F.w(al))
    return np.array(Zs), np.array(Ws), acc / n


def fd_covariance(Phi, D, T):
    """The low-noise covariance of w under π: T·D^{1/2}PD^{1/2}, P the orthogonal projector onto range(D^{1/2}Φᵀ) — so Cov(Φw) = TΦDΦᵀ,
    and on a support S of full column rank with D negligible off S, Cov(w_S) = T·D_S."""
    s = np.sqrt(np.asarray(D, float))
    M = s[:, None] * Phi.T                                                     # D^{1/2}Φᵀ (d × N)
    Q, R = np.linalg.qr(M)
    r = np.abs(np.diag(R)) > 1e-12 * np.abs(np.diag(R)).max()
    Q = Q[:, r]
    return T * (s[:, None] * (Q @ Q.T) * s[None, :])


def ito_q_drift(L, al, beta, Phi, sig):
    """dQ_i/dt when the noise is read in Itô's sense on (α, log β) (Proposition D‴.45(iv)(b) per unit time):
    (γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²))·σ²|φ_i|²/N²."""
    a, da = L.amp(al)
    N = Phi.shape[0]
    return 0.5 * L.gamma * beta * (L.gamma * a ** 2 - L.G * beta * da * (1 + (a / L.A) ** 2)) * sig ** 2 * (Phi ** 2).sum(0) / N ** 2


def tangent_projectors(path, Phi, tol=0.0):
    """Per segment of an `iss_path`, the projector onto ker Φ_Aᵀ (the tangent space of the active face, A = supp W_k)."""
    N = Phi.shape[0]
    out = []
    for w in path['W']:
        A = np.flatnonzero(np.abs(w) > tol)
        if A.size == 0:
            out.append(np.eye(N))
        else:
            Q, _ = np.linalg.qr(Phi[:, A])
            out.append(np.eye(N) - Q @ Q.T)
    return out


def linear_response(path, Phi, taus, xi_tau, kappa, return_nu=False):
    """The transit's noise to first order along an `iss_path` of entries only: δν ← P_{T(τ)}(δν + κΔW) on the grid `taus` (increments
    xi_tau ~ N(0, Δτ)) — the active face's normal part absorbed by the active kernels at every step, the accumulated deviation projected
    at each entry, so that δν(τ) = κP_{T(τ)}W(τ). An exit adds a first-order jump this response does not carry: a path with exits is
    refused. Returns δu = Φᵀδν on the grid (with return_nu, δν)."""
    if any(kind != 'enter' for t, i, kind in path['events']):
        raise ValueError('linear_response: the path has exits, whose first-order jump is not included')
    P = tangent_projectors(path, Phi)
    seg = np.clip(np.searchsorted(path['taus'], taus[1:], side='right') - 1, 0, len(P) - 1)
    dnu = np.zeros((taus.size, Phi.shape[0]))
    for k in range(taus.size - 1):
        dnu[k + 1] = P[seg[k]] @ (dnu[k] + kappa * xi_tau[k])
    return dnu if return_nu else dnu @ Phi


def normal_ou_var(Phi, D, act, kappa, Lam):
    """The finite-ε residual of the first-order response: the active block absorbs the noise's normal part on the time scale of its
    relaxation, not at once, so the dual keeps an Ornstein–Uhlenbeck component in range(Φ_A) — dν_n = −Λ(Φ_AD_AΦ_Aᵀ/N)ν_ndτ + κP_ndW —
    with the stationary variance tr((Φ_AD_AΦ_Aᵀ/N)⁺)κ²/(2Λ) = κ²Σ_i 1/(2Λλ_i) (λ_i the nonzero eigenvalues). Returns that trace (the
    predicted mean square of |δν − δν_lin|)."""
    A = np.flatnonzero(act)
    if A.size == 0:
        return 0.0
    N = Phi.shape[0]
    M = Phi[:, A] @ (np.asarray(D, float)[A][:, None] * Phi[:, A].T) / N
    lam = np.linalg.eigvalsh(M)
    lam = lam[lam > 1e-12 * lam.max()]
    return float(kappa ** 2 * np.sum(1.0 / (2 * Lam * lam)))
