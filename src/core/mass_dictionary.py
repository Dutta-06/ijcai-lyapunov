"""
Proposition D‴.44 (PARADIGM.md §5d): the replicator's own masses on a kernel dictionary under measurement noise — the conserved mass
inflates the fit; a null kernel restores the colour's selection
=======================================================================================================================================

The loop of Proposition D‴.34 at fixed positions on a dictionary φ_1, …, φ_d sampled at N < d nodes (the design Φ, N × d): masses β_i > 0
(Σβ = m — the replicator (𝔉) — or free, the unbalanced law), amplitudes a_i = a_max tanh(α_i/a_max), coefficients w = β⊙a/m,
J = (m/2N)|Φw − y|², the laws

    dα_i/ds = Gβ_ia′_i⟨rφ_i⟩,    d log β_i/ds = γ(u_i − λ),    u_i = a_i⟨rφ_i⟩,  λ = Σβ_ju_j/Σβ_j  (λ ≡ 0 for free masses),

driven by the measured error r + ξ (label noise, as in Proposition D‴.43). In the flat coordinates θ = (α, ϱ), ϱ = √β, the flow is the
gradient flow of J (of its 0-homogeneous extension in ϱ when the mass is conserved) in the metric |dα|²/G + (4/γ)|dϱ|², and label noise
selects along the exact fits M by the gradient flow of R = tr Ψ(K) (Proposition D‴.41(iv)),

    K = (m/N)Φ diag(D)Φᵀ,   D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²)   (the atom's relevance: its mass times its two sensitivities),

in s = (m v/4)t — for the flat update, an Euler step on (α, √β). The paradigm's own update — an Euler step on (α, log β), renormalised —
adds the replicator drift of the fitness (Σ_f)_{ℓ_iℓ_i}/4 (the increment's second moment in ℓ = log β), γa_i²/(2Gβ_ia′_i² + γa_i²) of an
atom's mass drift — O(1) where the amplitudes saturate —, and is not a gradient on M (referee round 34). Free masses, the flat update: R
increasing in every D_i, the gauge settles each atom at its minimal relevance D*(|w_i|) ≈ (γa_max/m)|w_i| above γa_max³/(Gm). Free
masses, the paradigm's update: the selection is each atom's Euler lag behind its own exact flow; it moves the free law's conserved
Q_i = β_i − c cosh(2α_i/a_max) (Proposition D‴.45) at (ω_iγβ_i/(Nm))(2γa_i² − 2Gβ_i(1 − t_i⁴)) per unit s, for every colour, so the gauge
rests at Gβ_i(1 − t_i⁴) = γa_i² (D*'s split is 2Gβ_i(1 − t_i⁴) = γa_i²), where the mirror coordinates drift by ω_iζ(w_i),
ζ = −(2γa_max/(Nm))t(1 − t²)/(1 + t²): a logarithmic penalty's gradient above the threshold. Conserved mass: the excess
E = m(1 − ‖w‖₁/a_max) stays unsaturated, Σβ_ia′_i ≥ E, and tr K ≥ (G/(Nm))E²/Σ|φ_i|⁻²; the paradigm's centring adds the cross fitness
−(2γ²/(Nm))a_iφ_iᵀf(K)y. A null kernel (φ_0 ≡ 0 on the nodes) parks the excess at no cost.

This module: the relevance, K, the weights and the potential, the selection's drift on M (flat, and with the update's Itô term), its
integration with a retraction onto M, the paradigm's noisy scheme (paths in parallel) and its noiseless flow (the limit map onto M), the
gauge-minimal relevance, the free law's gauge invariant, its drift and the gauge equilibria, the conserved law's cross term, the excess
mass and its bound, and the general Katzenberger drift of the scheme for checks (`gauge_drift.limit_drift` on a chart of M).
"""
import time

import numpy as np

from . import gauge_drift as GDR
from . import label_noise as LN


class FlowStopped(Exception):
    """`MassDictionary.flow` stopped before its end; `s` and `th`: the last state reached (the end of the last finished piece, on M)."""

    def __init__(self, msg, s=None, th=None):
        super().__init__(msg)
        self.s, self.th = s, th


class FlowTimeout(FlowStopped):
    """`MassDictionary.flow` over its wall-clock limit (max_secs)."""


class FlowFailed(FlowStopped):
    """A piece of `MassDictionary.flow` that neither LSODA nor Radau finished (solve_ivp's status; referee round 34)."""


class MassDictionary:
    """The replicator–amplitude law on the design Phi (N × d) with target y; `balanced` (Σβ = m) or free masses; θ = (α, ϱ)."""

    def __init__(self, Phi, y, m, a_max=1.0, G=1.0, gamma=0.1, balanced=True):
        self.Phi, self.y = np.asarray(Phi, float), np.asarray(y, float)
        self.m, self.a_max, self.G, self.gamma, self.balanced = float(m), float(a_max), float(G), float(gamma), balanced
        self.N, self.d = self.Phi.shape
        self.Pm = np.r_[np.full(self.d, self.G), np.full(self.d, self.gamma / 4)]        # the mobility in θ (the metric's inverse)

    # --- the state ---------------------------------------------------------------------------------------------------------------
    def parts(self, th):
        """(α, ϱ, β, a, a′, a″) at θ."""
        d, am = self.d, self.a_max
        al, rho = th[:d], th[d:]
        a = am * np.tanh(al / am)
        da = 1.0 - (a / am) ** 2
        return al, rho, rho * rho, a, da, -2.0 * a * da / am ** 2

    def w(self, th):
        _, _, b, a, _, _ = self.parts(th)
        return b * a / self.m

    def theta(self, al, beta):
        return np.r_[np.asarray(al, float), np.sqrt(np.asarray(beta, float))]

    # --- the selection ------------------------------------------------------------------------------------------------------------
    def relevance(self, th):
        """D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²): the atom's column in the Gauss–Newton matrix of the flat coordinates."""
        _, _, b, a, da, _ = self.parts(th)
        return (b / self.m ** 2) * (self.G * b * da ** 2 + self.gamma * a ** 2)

    def gram(self, th):
        """K = (m/N)Φ diag(D)Φᵀ, whose nonzero spectrum is the normal rates of the flow at a point of M; with the conserved mass
        K − (γ/N)yyᵀ: the uniform scaling of the masses, which moves f̂ along y, is not a direction of the flow."""
        K = (self.m / self.N) * (self.Phi * self.relevance(th)[None, :]) @ self.Phi.T
        return K - (self.gamma / self.N) * np.outer(self.y, self.y) if self.balanced else K

    def weights(self, th, filt):
        """ω_i = φ_iᵀ f(K) φ_i for the filter `filt` (a function of K's eigenvalues)."""
        kap, Q = np.linalg.eigh(self.gram(th))
        return ((Q.T @ self.Phi) ** 2 * filt(np.clip(kap, 0.0, None))[:, None]).sum(axis=0)

    def potential(self, th, Psi):
        """R = tr Ψ(K) for the function Ψ (e.g. λ for white draws, log(1 + τλ)/τ for a slow Markov noise)."""
        kap = np.clip(np.linalg.eigvalsh(self.gram(th)), 0.0, None)
        return float(np.sum(Psi(kap)))

    def grad(self, th, filt):
        """∇_θ R (R = tr Ψ(K), Ψ′ = filt): (m/N)ω_i ∂D_i/∂θ."""
        al, rho, b, a, da, dda = self.parts(th)
        om = (self.m / self.N) * self.weights(th, filt)
        m2 = self.m ** 2
        dDa = (b / m2) * (2 * self.G * b * da * dda + 2 * self.gamma * a * da)
        dDb = (2 * self.G * b * da ** 2 + self.gamma * a ** 2) / m2
        return np.r_[om * dDa, om * dDb * 2 * rho]

    # --- the manifold of exact fits -------------------------------------------------------------------------------------------------
    def constraints(self, th):
        c = self.Phi @ self.w(th) - self.y
        return np.r_[c, float(th[self.d:] @ th[self.d:]) - self.m] if self.balanced else c

    def cjac(self, th):
        _, rho, b, a, da, _ = self.parts(th)
        J = np.c_[self.Phi * (b * da / self.m)[None, :], self.Phi * (2 * rho * a / self.m)[None, :]]
        return np.r_[J, np.r_[np.zeros(self.d), 2 * rho][None, :]] if self.balanced else J

    RCOND = 1e-10                                                                 # the projection's relative cut-off (see `project`)

    def project(self, th, q):
        """The flow-metric orthogonal projection of a velocity q onto the tangent space of M: q − P Jᵀ(J P Jᵀ)⁻¹ J q (P the mobility),
        computed as q − P^{1/2}A·A⁺P^{−1/2}q with A = P^{1/2}Jᵀ by least squares on A itself — not on the normal equations J P Jᵀ, whose
        condition number is A's squared: atoms with masses near e^{−30} make J P Jᵀ singular to rounding, the right-hand side noisy and the
        integrator's step collapse. Directions with singular values below RCOND·max (atoms too light to carry a constraint) are dropped."""
        s = np.sqrt(self.Pm)
        A = s[:, None] * self.cjac(th).T
        return q - s * (A @ np.linalg.lstsq(A, q / s, rcond=self.RCOND)[0])

    def retract(self, th, tol=1e-13, it=30):
        """Back onto M by Gauss–Newton steps in the flow's metric (the minimal-norm correction P Jᵀ(J P Jᵀ)⁻¹c = P^{1/2}(Aᵀ)⁺c)."""
        s = np.sqrt(self.Pm)
        for _ in range(it):
            c = self.constraints(th)
            if np.abs(c).max() < tol:
                break
            A = s[:, None] * self.cjac(th).T
            th = th - s * np.linalg.lstsq(A.T, c, rcond=self.RCOND)[0]
        return th

    def ito_term(self, th, Sig_ll=None, filt=None):
        """The paradigm's update (Euler on log β, renormalised) against the flat one: the drift (in θ, per unit v) whose part along M it
        adds — in ℓ = log β the fitness diag(Σ_f,ℓℓ)/4, Σ_f the increment's filtered second moment. Label noise is ∝ the Gauss–Newton
        matrix, diagonal in its eigenbasis, so the filter weighs each normal mode by f(κ): the noise fields in ℓ being U_ℓ,i(x) =
        γ(a_iφ_i(x) − f̂_λ(x)) (f̂_λ = Σβ_ja_jφ_j/Σβ_j with the conserved mass, 0 without), Σ_f,ℓℓ = (γ²/N)·X Φᵀf(K)Φ Xᵀ, X = diag(a) − (1/m)
        1(β⊙a)ᵀ (K the conserved mass's K̃); `filt` the filter (None: fresh draws as h → 0, f ≡ 1), or `Sig_ll` given directly."""
        al, rho, b, a, da, dda = self.parts(th)
        if Sig_ll is None:
            f = (lambda k: np.ones_like(k)) if filt is None else filt
            kap, Q = np.linalg.eigh(self.gram(th))
            M = (self.Phi.T @ Q) @ np.diag(f(np.clip(kap, 0.0, None))) @ (Q.T @ self.Phi)       # Φᵀ f(K) Φ
            X = np.diag(a)
            if self.balanced:
                X = X - np.outer(np.ones(self.d), b * a) / b.sum()
            Sig_ll = (self.gamma ** 2 / self.N) * np.einsum('ij,jk,ik->i', X, M, X)
        return np.r_[np.zeros(self.d), 0.5 * rho * np.asarray(Sig_ll, float) / 4.0]

    def drift(self, th, filt, ito=False):
        """The selection's drift along M per unit s = (m v/4)t: −P_g(P∇R) (the flat scheme), plus, with `ito`, the paradigm's update's
        Itô term at the same filter (`ito_term`)."""
        q = -self.project(th, self.Pm * self.grad(th, filt))
        if ito:
            q = q + self.project(th, self.ito_term(th, filt=filt) * (4.0 / self.m))
        return q

    def flow(self, th, S, filt, ito=False, t_eval=None, chunks=40, rtol=1e-8, atol=1e-10, method='LSODA', max_secs=None):
        """The selection flow over s ∈ [0, S] from θ on M → (s, Θ (len(s), 2d)), integrated in z = (α, log β) (an atom's mass decays
        exponentially once the selection drops it: linear in log β), piecewise, retracted onto M at the pieces' ends. max_secs: a
        wall-clock limit, checked at every evaluation (FlowTimeout): on one rbf instance of D‴.44's Part D a piece ran for hours, LSODA's
        steps collapsing. A piece the method does not finish is redone with Radau, and if that fails too the flow stops (FlowFailed);
        both carry the last state reached."""
        from scipy.integrate import solve_ivp
        d, lfloor = self.d, -1400.0                                              # a dead atom's mass read as ≥ e^{−1400}: its log-mass
        to_z = lambda x: np.r_[x[:d], 2 * np.log(x[d:])]                         # keeps decaying linearly, without underflow
        to_th = lambda z: np.r_[z[:d], np.exp(0.5 * np.maximum(z[d:], lfloor))]

        t_end = np.inf if max_secs is None else time.time() + max_secs
        last = [0.0, None]                                                       # the end of the last finished piece

        def rhs(s, z):
            if time.time() > t_end:
                raise FlowTimeout(f'flow over {max_secs} s at s = {s:g}', s=last[0], th=last[1])
            x = to_th(z)
            if not np.all(np.isfinite(x)):
                return np.full_like(z, np.nan)                                  # an overshooting trial stage: rejected
            try:
                q = self.drift(x, filt, ito)
            except np.linalg.LinAlgError:                                       # a trial stage far off M: rejected
                return np.full_like(z, np.nan)
            return np.r_[q[:d], 2 * q[d:] / x[d:]]
        edges = np.linspace(0.0, S, chunks + 1)
        t_eval = edges if t_eval is None else np.asarray(t_eval, float)
        th = self.retract(np.asarray(th, float))
        last[1] = th
        ts, ys = [], []
        for j, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            final = j == chunks - 1
            want = t_eval[(t_eval >= lo) & ((t_eval <= hi) if final else (t_eval < hi))]
            te = np.unique(np.r_[want, hi])
            sol = solve_ivp(rhs, (lo, hi), to_z(th), method=method, rtol=rtol, atol=atol, t_eval=te)
            if not sol.success:
                sol = solve_ivp(rhs, (lo, hi), to_z(th), method='Radau', rtol=rtol, atol=atol, t_eval=te)
                if not sol.success:
                    raise FlowFailed(f'flow: the piece [{lo:g}, {hi:g}] failed ({sol.message})', s=last[0], th=last[1])
            th = self.retract(to_th(sol.y[:, -1]))
            # the retraction's correction is additive in β (the metric's mass block is constant), so on an atom the selection has
            # dropped (β ≈ e^{−700}) it can cross zero; such an atom is kept at the floor, a change far below the constraint's tolerance
            th[d:] = np.maximum(th[d:], np.exp(0.5 * lfloor))
            if not np.all(np.isfinite(th)):
                raise FloatingPointError(f'flow: non-finite state at s = {hi:g}')
            last[0], last[1] = hi, th
            for s_, z_ in zip(sol.t, sol.y.T):
                if np.any(np.isclose(s_, want)):
                    ts.append(s_); ys.append(self.retract(to_th(z_)))
        return np.array(ts), np.array(ys)

    # --- the gauge, the excess, the bound -----------------------------------------------------------------------------------------------
    def excess(self, th):
        """E = m(1 − ‖w‖₁/a_max): the mass the fit leaves unused at saturation (conserved mass)."""
        return self.m * (1.0 - np.abs(self.w(th)).sum() / self.a_max)

    def white_bound(self, E):
        """R_white ≥ (G/(Nm))E²/Σ_i|φ_i|⁻² at every exact fit with excess E (Cauchy–Schwarz on Σβ_ia′_i ≥ E)."""
        return self.G / (self.N * self.m) * E ** 2 / float(np.sum(1.0 / (self.Phi ** 2).sum(axis=0)))

    def relevance_min(self, wabs):
        """D*(|w|) = min over the gauge (β = m|w|/|a|) of the relevance, free masses: min over t = |a|/a_max ∈ (0, 1) of
        G w²(1 − t²)²/(a_max²t²) + γ|w|a_max t/m (bounded scalar minimisation)."""
        from scipy.optimize import minimize_scalar
        G, g, am, m = self.G, self.gamma, self.a_max, self.m
        out = []
        for x in np.atleast_1d(np.asarray(wabs, float)):
            if x == 0:
                out.append(0.0); continue
            f = lambda t: G * x * x * (1 - t * t) ** 2 / (am * am * t * t) + g * x * am * t / m
            out.append(float(minimize_scalar(f, bounds=(1e-9, 1.0), method='bounded', options=dict(xatol=1e-13)).fun))
        return np.array(out)

    # --- free masses: the gauge invariant, its drift, the gauge equilibria (referee round 34) ----------------------------------------
    def gauge_invariant(self, th):
        """Q_i = β_i − c cosh(2α_i/a_max), c = γa_max²/(4G): constant along the noiseless free law for every drive (Proposition D‴.45(i)),
        with zero quadratic variation under label noise (the noise moves each atom along its own flow line)."""
        al, _, b, _, _, _ = self.parts(th)
        return b - self.gamma * self.a_max ** 2 / (4 * self.G) * np.cosh(2 * al / self.a_max)

    def gauge_ratio(self, th):
        """Gβ_i(1 − t_i⁴)/(γa_i²), t = a/a_max: Q_i's drift under the free masses' selection is (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴))
        per unit s (`q_drift`), so an alive atom rests at the ratio 1/2 under the flat update (κ = 1: D*'s split) and 1 under the
        paradigm's (κ = 2)."""
        _, _, b, a, _, _ = self.parts(th)
        t = a / self.a_max
        return self.G * b * (1 - t ** 4) / (self.gamma * a ** 2)

    def gauge_t(self, wabs, paradigm=True):
        """The gauge equilibrium's t = |a|/a_max at |w|: kGm|w|(1 − t⁴) = γa_max³t³ with k = 1 (the paradigm's update) or 2 (the flat one,
        D*'s minimiser) — the zero of `q_drift` on the orbit β = m|w|/(a_max t)."""
        from scipy.optimize import brentq
        k = 1.0 if paradigm else 2.0
        G, g, am, m = self.G, self.gamma, self.a_max, self.m
        return np.array([0.0 if x == 0 else brentq(lambda t: k * G * m * x * (1 - t ** 4) - g * am ** 3 * t ** 3, 0.0, 1.0, xtol=1e-15, rtol=1e-15)
                         for x in np.atleast_1d(np.asarray(wabs, float))])

    def q_drift(self, th, filt, ito=False):
        """Q's drift (`gauge_invariant`) along the free masses' selection per unit s, closed form: (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)),
        κ = 2 with the paradigm's Itô term, 1 without; ω_i = φ_iᵀf(K)φ_i — every colour. Q is conserved by the noiseless flow for every
        target, so the selection moves it only through the second-order term (1/2)∇²Q:Σ_f, whose block on atom i is the atom's own."""
        assert not self.balanced
        _, _, b, a, _, _ = self.parts(th)
        t = a / self.a_max
        return self.weights(th, filt) * self.gamma * b / (self.N * self.m) * ((2.0 if ito else 1.0) * self.gamma * a ** 2
                                                                              - 2 * self.G * b * (1 - t ** 4))

    def euler_lag(self, th, filt):
        """Free masses, the paradigm's update: its selection's drift before the projection onto M, per unit s, in θ — −(2ω_i/(Nm))·dv_i/dζ,
        v_i = (Gβ_ia′_i, γa_i) atom i's flow direction in (α, log β) per unit drive ζ: the Euler step's second-order lag behind the exact
        flow along the atom's own curve (which keeps Q and moves the mirror coordinate by exactly the drive)."""
        _, rho, b, a, da, dda = self.parts(th)
        om = self.weights(th, filt)
        lag_a = self.G * b * da * (self.gamma * a + self.G * b * dda)                 # dv/dζ, α-component
        lag_l = self.gamma * self.G * b * da ** 2                                     # dv/dζ, log β-component
        k = -2.0 * om / (self.N * self.m)
        return np.r_[k * lag_a, k * lag_l * 0.5 * rho]                                # dϱ = (ϱ/2)d log β

    def mirror_rate_eq(self, wabs):
        """ζ(|w|) at the paradigm's gauge equilibrium: the drift of the atom's mirror coordinate η = ∫dα/(Gβa′) (Proposition D‴.45(ii)) per
        unit ω, −(2γa_max/(Nm))t(1 − t²)/(1 + t²), t = `gauge_t`(|w|) (Q's drift vanishes there, so ζ does not depend on the mirror
        coordinate's normalisation); ≈ −γ²a_max⁴/(2NGm²|w|) above the threshold γa_max³/(Gm), −(2γa_max/(Nm))t below."""
        t = self.gauge_t(wabs, paradigm=True)
        return -(2 * self.gamma * self.a_max / (self.N * self.m)) * t * (1 - t * t) / (1 + t * t)

    def cross_term(self, th, filt):
        """The conserved law's leftover after the masses' own selection cancels: the centring's cross term in the Itô fitness, in θ per
        unit v (as `ito_term`): ϱ_i(γ²/N)(−2a_iφ_iᵀf(K)f̂)/8, f̂ = Φw (= y on M) — it favours atoms whose contribution a_iφ_i opposes the
        target in f(K)'s inner product."""
        _, rho, b, a, _, _ = self.parts(th)
        kap, Q = np.linalg.eigh(self.gram(th))
        fy = Q @ (filt(np.clip(kap, 0.0, None)) * (Q.T @ (self.Phi @ self.w(th))))
        return np.r_[np.zeros(self.d), 0.5 * rho * (self.gamma ** 2 / self.N) * (-2 * a * (self.Phi.T @ fy)) / 4.0]

    def settle(self, al, beta, T=2e5, rtol=1e-11, atol=1e-13):
        """The noiseless flow (the loop's own law, in (α, log β)) from (α, β) over t ∈ [0, T], retracted onto M: Katzenberger's limit map
        Φ when T is long against the slowest normal rate, the start of the selection's flows (Parts B–E)."""
        from scipy.integrate import solve_ivp
        d = self.d

        def F(z):
            x = np.r_[z[:d], np.exp(z[d:] / 2)]
            _, _, b, a, da, _ = self.parts(x)
            g = self.Phi.T @ (self.y - self.Phi @ (b * a / self.m)) / self.N
            u = a * g
            lam = b @ u / b.sum() if self.balanced else 0.0
            return np.r_[self.G * b * da * g, self.gamma * (u - lam)]
        z = solve_ivp(lambda t, z: F(z), (0.0, T), np.r_[np.asarray(al, float), np.log(np.asarray(beta, float))], method='LSODA',
                      rtol=rtol, atol=atol).y[:, -1]
        return self.retract(np.r_[z[:d], np.exp(z[d:] / 2)])

    # --- the noisy scheme -----------------------------------------------------------------------------------------------------------
    def noisy(self, al0, beta0, h, sigma, rho, n_steps, rng, paths, rec_every, keep=()):
        """The paradigm's scheme with label noise, `paths` in parallel: α += hGβa′⟨(r + ξ)φ⟩, log β += hγ(u − λ) renormalised (conserved
        mass), ξ AR(1) per node with per-step variance σ² and step correlation ρ (ρ = 0: drawn afresh), started stationary → the records
        (n_rec, paths, d) of w and (α, β) at the end; with `keep` (step numbers) also the list of the states (α, β) after those steps."""
        N, d, m, am = self.N, self.d, self.m, self.a_max
        al = np.tile(np.asarray(al0, float), (paths, 1)); lb = np.tile(np.log(np.asarray(beta0, float)), (paths, 1))
        xi = sigma * rng.standard_normal((paths, N)); s_in = np.sqrt(1 - rho * rho) * sigma
        rec, kept, keep = [], [], set(keep)
        for k in range(n_steps + 1):
            b = np.exp(lb); a = am * np.tanh(al / am)
            if k % rec_every == 0:
                rec.append(b * a / m)
            if k in keep:
                kept.append((al.copy(), b.copy()))
            if k == n_steps:
                break
            g = (self.y[None, :] - (b * a / m) @ self.Phi.T + xi) @ self.Phi / N
            u = a * g
            lam = ((b * u).sum(axis=1, keepdims=True) / b.sum(axis=1, keepdims=True)) if self.balanced else 0.0
            al = al + h * self.G * b * (1 - (a / am) ** 2) * g
            lb = lb + h * self.gamma * (u - lam)
            if self.balanced:
                lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / m)
            xi = rho * xi + s_in * rng.standard_normal((paths, N))
        return (np.array(rec), al, np.exp(lb)) + ((kept,) if keep else ())

    # --- checks: the general drift of the scheme --------------------------------------------------------------------------------------
    def scheme_drift(self, th, basis, h=0.0, corr=('fresh',), coords='log', hh=1e-5, hh2=1e-3):
        """Katzenberger's drift (per unit v, in θ) of the Euler scheme in z = (α, log β) (coords='log', the paradigm's, renormalised) or in
        θ (coords='flat', renormalised radially), for draws `corr` at the step h (h = 0: the fresh limit), by `gauge_drift.limit_drift` on
        the chart of M that solves the `basis` atoms' amplitudes."""
        d, N, m, am = self.d, self.N, self.m, self.a_max
        B = np.asarray(basis); Fr = np.setdiff1d(np.arange(d), B)
        to_z = (lambda x: np.r_[x[:d], 2 * np.log(x[d:])]) if coords == 'log' else (lambda x: x)
        from_z = (lambda z: np.r_[z[:d], np.exp(z[d:] / 2)]) if coords == 'log' else (lambda z: z)

        def F(z):
            al, rho, b, a, da, _ = self.parts(from_z(z))
            g = self.Phi.T @ (self.y - self.Phi @ (b * a / m)) / N
            u = a * g
            lam = b @ u / b.sum() if self.balanced else 0.0
            dl = self.gamma * (u - lam)
            return np.r_[self.G * b * da * g, dl if coords == 'log' else 0.5 * rho * dl]

        def chart(c):
            l, aF = c[:d], c[d:]
            b = np.exp(l)
            al = np.empty(d); al[Fr] = aF
            rhs = m * self.y - self.Phi[:, Fr] @ (b[Fr] * am * np.tanh(aF / am))
            aB = np.linalg.solve(self.Phi[:, B] * b[B][None, :], rhs)
            al[B] = am * np.arctanh(np.clip(aB / am, -1 + 1e-15, 1 - 1e-15))
            return to_z(np.r_[al, np.sqrt(b)])

        def coords_of(z):
            x = from_z(z)
            return np.r_[2 * np.log(x[d:]), x[:d][Fr]]

        z0 = to_z(np.asarray(th, float))
        c0 = coords_of(z0)
        al, rho, b, a, da, _ = self.parts(th)
        U = np.c_[self.G * (b * da)[None, :] * self.Phi, self.gamma * a[None, :] * self.Phi]
        if self.balanced:
            U[:, d:] -= (U[:, d:] @ b / b.sum())[:, None]
        if coords != 'log':
            U[:, d:] *= 0.5 * rho[None, :]
        Sig = U.T @ U / N
        A = GDR._jac(F, z0, hh)
        Jg = GDR._jac(chart, c0, hh)
        P = GDR.spectral_projection(A, Jg)[0]
        Sf = Sig if h == 0 else LN.filtered_covariance(A, P, Sig, h, corr)
        F1 = None
        if self.balanced:
            if coords == 'log':
                sh = -0.5 * float(b @ np.diag(Sf)[d:]) / b.sum()
                F1 = lambda x: np.r_[np.zeros(d), np.full(d, sh)]
            else:
                F1 = lambda x: np.r_[np.zeros(d), -x[d:] * np.trace(Sf[d:, d:]) / (2 * m)]
        dc = GDR.limit_drift(F, chart, coords_of, c0, Sf, h=hh, h2=hh2, F1=F1)[0]
        dz = GDR._jac(chart, c0, 1e-6) @ dc
        return np.r_[dz[:d], 0.5 * rho * dz[d:]] if coords == 'log' else dz, Sf
