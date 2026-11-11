"""
Propositions D‴.45–D‴.48 (PARADIGM.md §5d): the free mass–amplitude law is a mirror flow — a conservation law per particle, the implicit
bias as a Bregman projection, the noise's invariance in continuous time (noise with smooth paths: Stratonovich), the sampled-data drift,
the real loop's Lyapunov function; a leak in the mirror coordinates; the law's branches; the rich transit
====================================================================================================================================

The replicator–amplitude law of Proposition D‴.34 without the mass constraint (free masses: λ ≡ 0), one particle per class beyond the
cohesion's range (the cohesion's α-component vanishing), each particle driven by one scalar signal g_i — ⟨(r + ξ)φ_i⟩ on a dictionary,
eφ_i(x) in the real loop:

    dα_i = Gβ_ia′_ig_i dt,    dβ_i = γβ_ia_ig_i dt,    a_i = a_max tanh(α_i/a_max),   w_i = β_ia_i/m.

D‴.45:
(i)   Q_i := β_i − c cosh(2α_i/a_max), c := γa_max²/(4G), is constant for every drive (a noise path read in Stratonovich's sense).
(ii)  On a level set, η_i := ∫dα_i/(Gβ_ia′_i) obeys dη_i = g_i dt and dw_i/dη_i = mD_i, D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²): w follows the
      mirror flow of the separable potential Ψ_q(w) = Σ_iψ_{q_i}(w_i), ψ′ = η, whose conjugate is the mass gained, ψ*(η) = (β − β⁰)/(γm)
      = 2c sinh²(α/a_max)/(γm) (β⁰ = q + c the mass at α = 0), and ψ(w) = ηw − ψ*.
(iii) On a dictionary the noiseless loop's limit is the Bregman projection argmin{Ψ_q(w) − ⟨η(0), w⟩ : Φw = y} (every Q_i > −c);
      noise with smooth paths moves η only within η(0) + range(Φᵀ), so the projection of the current state never moves.
(iv)  An Euler step on (α, log β) changes Q_i by (h²g_i²γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) + O(h³) (`dQ_euler`; `dQ_exact` the
      whole change) — the Itô correction of white noise per unit time; the mirror update η⁺ = η + hg on each particle's curve
      (`mirror_step`) keeps Q and the projection at every step.
(v)   In the real loop ė = −k_ee + f(x) − f̂(x), f = Σ_iw*_iφ_i, V = ½e² + Σ_i[w*_i(η*_i − η_i) − (β*_i − β_i)/(γm)] (the Bregman divergence,
      (η*, β*) on the particle's own curve at w*_i) has V̇ = −k_ee².
(vi)  The conserved law (Σβ = m) shifts every level: dQ_i = −γβ_iλ dt, λ = Σ_jw_jg_j.
D‴.46: the leak −σ(η_i − η_p,i) on each particle's drive keeps Q; with a disturbance d, V̇ = −k_ee² + ed − σΣ[B(w*, w) + B(w, w_p)
      − B(w*, w_p)]; the averaged loop's equilibrium solves the monotone equation (H_L(w − w*) − d̄_L)/k_e + σ(η − η_p) = 0 with the
      error filter's Gram H_L = k_e⟨φL[φ]ᵀ⟩ and d̄_L = k_e⟨φL[d]⟩ (L = 1/(s + k_e), `filtered_gram`; `leak_equilibrium`).
D‴.47: the branches Q_i < −c (a start with β(0) < 2c sinh²(α(0)/a_max)): the sign frozen, an entropy at the branch point with the
      temperature 1/(γa₀); the limit for every start is the Bregman projection onto the sign-constrained least-squares fits, the
      coordinates vanishing there going extinct (`ParticleSet`, `sign_ls_set`, `projection`); the rich limits' leading order from the
      vertex's certificate (`certificate_slacks`, `rich_leading_order`, `FreeLaw.core_tail`).
D‴.48: the rich transit — as the starting masses vanish the rescaled loop follows the inverse scale space flow (`iss_path`).

This module: the curves, the mirror coordinates and their inverse (closed forms without cancellation: q > c, |q| ≤ c, the critical curve
q = −c and the branches q < −c; the inverse raises instead of clipping), the potential and its conjugate, the Bregman projection on a
dictionary (Newton on the dual with Armijo backtracking and an adaptive radius; a decimal-arithmetic step where float curvatures
underflow), the loop's vector fields, the Euler step's change of Q, the mirror update, the real loop and its Lyapunov function (from the
state's own mass), the leak and its equilibrium, the error filter's Gram, and the branch-aware particle set with the sign-constrained
projection, weighted basis pursuit, the nonnegative certificate and the irrepresentable condition.
"""
import numpy as np
from scipy.optimize import brentq, lsq_linear, linprog

ALPHA_MAX = 340.0            # |α| beyond which cosh(2α/a_max) overflows (a_max = 1): the inverse raises past it


def logsinh(x):
    """log sinh x for x > 0, without underflow near 0 or overflow at large x."""
    x = np.asarray(x, float)
    small = x < 1e-4
    xs = np.where(small, np.maximum(x, 1e-300), 1.0)
    xl = np.where(small, 1.0, x)
    big = np.log(-np.expm1(-2 * xl)) + xl - np.log(2.0)
    return np.where(small, np.log(xs) + xs * xs / 6, big)


class FreeLaw:
    """The free mass–amplitude law's geometry for constants (m, G, γ, a_max); c = γa_max²/(4G)."""

    def __init__(self, m=12.0, G=1.0, gamma=0.1, a_max=1.0):
        self.m, self.G, self.gamma, self.A = float(m), float(G), float(gamma), float(a_max)
        self.c = self.gamma * self.A ** 2 / (4 * self.G)

    # --- the state ---------------------------------------------------------------------------------------------------------------
    def amp(self, al):
        """(a, a′) at α."""
        a = self.A * np.tanh(np.asarray(al, float) / self.A)
        return a, 1.0 - (a / self.A) ** 2

    def Q(self, al, beta):
        """The conserved quantity β − c cosh(2α/a_max)."""
        return np.asarray(beta, float) - self.c * np.cosh(2 * np.asarray(al, float) / self.A)

    def beta(self, al, q):
        """The mass on the level set Q = q at α."""
        return np.asarray(q, float) + self.c * np.cosh(2 * np.asarray(al, float) / self.A)

    def w(self, al, q):
        return self.beta(al, q) * self.amp(al)[0] / self.m

    def relevance(self, al, q):
        """D = (β/m²)(Gβa′² + γa²) on the curve: dw/dη = mD."""
        b = self.beta(al, q)
        a, da = self.amp(al)
        return (b / self.m ** 2) * (self.G * b * da ** 2 + self.gamma * a ** 2)

    def core_tail(self, b0, mu):
        """From α = 0 at the mass b0, the coefficient at η = (1 − μ)λ₁ — an off-support atom of the rich limit with the certificate's slack
        μ (Proposition D‴.45(iii)): a_max b0^{3/2}cos(πμ/2)/(m(2c)^{1/2}sin³(πμ/2)), to leading order as b0 → 0."""
        b0, mu = np.asarray(b0, float), np.asarray(mu, float)
        with np.errstate(divide='ignore'):                                     # μ = 0 (the wall): infinite
            return self.A * b0 ** 1.5 * np.cos(np.pi * mu / 2) / (self.m * np.sqrt(2 * self.c) * np.sin(np.pi * mu / 2) ** 3)

    def branch_point(self, q):
        """α₀ ≥ 0 with β(±α₀) = 0 when q ≤ −c (the sign of w frozen beyond it; α₀ = 0 on the critical curve q = −c), else 0."""
        q = np.asarray(q, float)
        return np.where(q + self.c < 0, 0.5 * self.A * np.arccosh(np.maximum(-q / self.c, 1.0)), 0.0)

    # --- the mirror coordinates ----------------------------------------------------------------------------------------------------
    def _eta_branch(self, dl, a0):
        """η on the positive side of a branch at distance δ = α − α₀ from its branch point α₀ > 0 (cancellation-free):
        (α₀ + δ − (a_max/2)coth(α₀/a_max)[log sinh((2α₀ + δ)/a_max) − log sinh(δ/a_max)])/(2Gc); α₀ = 0: the critical curve,
        (α − a_max coth(α/a_max))/(2Gc)."""
        A, G, c = self.A, self.G, self.c
        dl = np.asarray(dl, float); a0 = np.asarray(a0, float)
        crit = a0 == 0
        a0s = np.where(crit, 1.0, a0)
        v = a0s + dl - 0.5 * A / np.tanh(a0s / A) * (logsinh((2 * a0s + dl) / A) - logsinh(dl / A))
        dls = np.where(crit & (dl > 0), dl, 1.0)
        vc = dl - A / np.tanh(dls / A)
        return np.where(crit, vc, v) / (2 * G * c)

    def _eta_branch_log(self, u, a0):
        """`_eta_branch` at δ = e^u from u itself (log sinh(δ/a_max) = u − log a_max + O(δ²) when δ is tiny), so that a particle at
        any depth below its branch point (δ underflowing) keeps a finite η (a₀ > 0)."""
        A, G, c = self.A, self.G, self.c
        u = np.asarray(u, float); a0 = np.asarray(a0, float)
        dl = np.exp(np.minimum(u, np.log(ALPHA_MAX)))
        tiny = u < np.log(1e-4 * A)
        lsd = np.where(tiny, u - np.log(A) + dl * dl / (6 * A * A), logsinh(np.where(tiny, 1.0, dl) / A))
        v = a0 + dl - 0.5 * A / np.tanh(a0 / A) * (logsinh((2 * a0 + dl) / A) - lsd)
        return v / (2 * G * c)

    def eta(self, al, q):
        """η = ∫dα/(Gβa′) along the curve: from α = 0 when q > −c (odd in α); on a branch (q < −c) and on the critical curve (q = −c)
        a primitive that → −∞ at the branch point on the positive side and → +∞ on the negative side (η odd in α). Closed forms
        without cancellation: η = (α − a_maxκ artanh(κt))/(2Gc), t = tanh(α/a_max), κ² = (q − c)/(q + c) ∈ (0, 1) (written with
        δ = 1 − κ and 1 − t = 2e^{−2|α|/a_max}/(1 + e^{−2|α|/a_max})); (α + a_maxκ′arctan(κ′t))/(2Gc) when κ′² = −κ² ≥ 0; on a branch
        the form of `_eta_branch`. NaN inside a branch's gap |α| < α₀."""
        A, G, c = self.A, self.G, self.c
        al = np.array(al, float, ndmin=1)
        q = np.broadcast_to(np.asarray(q, float), al.shape)
        sa = np.where(al < 0, -1.0, 1.0)
        x = np.abs(al)
        t = np.tanh(x / A)
        ex = np.exp(-2 * x / A)
        u = 2 * ex / (1 + ex)                                                     # 1 − t
        out = np.full(al.shape, np.nan)
        qc = q + c
        free = qc > 0
        with np.errstate(divide='ignore', invalid='ignore'):
            k2 = np.where(free, (q - c) / np.where(free, qc, 1.0), 0.0)
        mid = free & (k2 > 0)
        if np.any(mid):
            km, tm, um = np.sqrt(k2[mid]), t[mid], u[mid]
            dl = (2 * c / qc[mid]) / (1 + km)                                     # δ = 1 − κ, without cancellation
            small = tm < 0.5
            # log(1 − κt): log1p(−κt) for small t, log(δ + κ(1 − t)) near |t| = 1
            l1m = np.where(small, np.log1p(-km * np.where(small, tm, 0.0)), np.log(dl + km * np.where(small, 0.5, um)))
            # α − a_maxκ artanh(κt) = (a_max/2)[log((1 + t)/(1 + κt)) − log((1 − t)/(1 − κt)) + δ log((1 + κt)/(1 − κt))], with
            # (1 + t)/(1 + κt) = 1 + δt/(1 + κt) and (1 − κt)/(1 − t) = 1 + δt/(1 − t): every term without cancellation
            v = np.log1p(dl * tm / (1 + km * tm)) + np.log1p(dl * tm / um) + dl * (np.log1p(km * tm) - l1m)
            out[mid] = sa[mid] * 0.5 * A * v / (2 * G * c)
        low = free & (k2 <= 0)
        if np.any(low):
            kp = np.sqrt(-k2[low])
            out[low] = sa[low] * (x[low] + A * kp * np.arctan(kp * t[low])) / (2 * G * c)
        nb = ~free
        if np.any(nb):
            a0 = self.branch_point(q[nb])
            inside = x[nb] > a0
            v = np.full(a0.shape, np.nan)
            if np.any(inside):
                v[inside] = self._eta_branch(x[nb][inside] - a0[inside], a0[inside])
            out[nb] = sa[nb] * v
        return out

    def alpha_of_eta(self, e, q, sgn=None):
        """The inverse of η on each particle's curve (the side of a branch or of the critical curve from sgn, default +): brentq in α
        (q > −c) or in u = log(|α| − α₀) (q ≤ −c, never crossing the branch point). Raises ValueError when the target lies beyond
        |α| = ALPHA_MAX (or, on a branch, beyond the double-precision distance to the branch point) — no silent clipping."""
        e = np.array(e, float, ndmin=1)
        shp = e.shape
        e = e.ravel()
        q = np.broadcast_to(np.asarray(q, float), shp).ravel()
        sgn = np.ones_like(e) if sgn is None else np.broadcast_to(np.asarray(sgn, float), shp).ravel()
        a0 = self.branch_point(q)
        out = np.empty_like(e)
        for j in range(e.size):
            if q[j] + self.c > 0:                                                 # η odd and increasing: solve for |α|
                tgt = abs(e[j])
                if tgt == 0:
                    out[j] = 0.0
                    continue
                f = lambda x: float(self.eta(np.array([x]), q[j])[0]) - tgt
                hi = 1.0
                while f(hi) < 0:
                    if hi >= ALPHA_MAX:
                        raise ValueError(f"alpha_of_eta: η = {e[j]:.6g} beyond |α| = {ALPHA_MAX} on the curve q = {q[j]:.6g}")
                    hi = min(2 * hi, ALPHA_MAX)
                out[j] = np.copysign(brentq(f, 0.0, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=500), e[j])
            else:
                s = 1.0 if sgn[j] >= 0 else -1.0
                tgt = s * e[j]                                                    # η₊(δ) = s·η
                f = lambda u: float(self._eta_branch(np.exp(u), a0[j])) - tgt
                lo, hi = -5.0, 0.0
                while f(lo) > 0:
                    if lo <= -740:
                        raise ValueError(f"alpha_of_eta: η = {e[j]:.6g} closer to the branch point than double precision (q = {q[j]:.6g})")
                    lo = max(2 * lo, -740.0)
                while f(hi) < 0:
                    if hi >= np.log(ALPHA_MAX):
                        raise ValueError(f"alpha_of_eta: η = {e[j]:.6g} beyond |α| = {ALPHA_MAX} on the branch q = {q[j]:.6g}")
                    hi = min(hi + 1.0, np.log(ALPHA_MAX))
                u = brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)
                out[j] = s * (a0[j] + np.exp(u))
        return out.reshape(shp)

    def psi(self, al, q):
        """ψ_q(w(α)) = ηw − (β − β⁰)/(γm), β − β⁰ = 2c sinh²(α/a_max) (q > −c; on a branch up to the constant of η)."""
        return self.eta(al, q) * self.w(al, q) - self.psi_star(al, q)

    def psi_star(self, al, q):
        """ψ*_q(η(α)) = (β − β⁰)/(γm) = 2c sinh²(α/a_max)/(γm): the conjugate potential is the mass gained."""
        return 2 * self.c * np.sinh(np.asarray(al, float) / self.A) ** 2 / (self.gamma * self.m)

    def bregman_div(self, al_star, al, q):
        """B(w*, w) = w*(η* − η) − (β* − β)/(γm) on each curve (η* = η(α*), β* = β(α*)), summed."""
        return float(np.sum(self.w(al_star, q) * (self.eta(al_star, q) - self.eta(al, q))
                            - (self.beta(al_star, q) - self.beta(al, q)) / (self.gamma * self.m)))

    def alpha_of_w(self, w, q, sgn=None):
        """α on the curve with w(α) = w (w is increasing in α on each side; on a branch the side is w's sign)."""
        w = np.array(w, float, ndmin=1)
        shp = w.shape
        w = w.ravel()
        q = np.broadcast_to(np.asarray(q, float), shp).ravel()
        a0 = self.branch_point(q)
        out = np.empty_like(w)
        for j in range(w.size):
            if q[j] + self.c > 0:
                f = lambda x: float(self.w(x, q[j])) - abs(w[j])
                hi = 1.0
                while f(hi) < 0:
                    if hi >= ALPHA_MAX:
                        raise ValueError(f"alpha_of_w: w = {w[j]:.6g} beyond |α| = {ALPHA_MAX}")
                    hi = min(2 * hi, ALPHA_MAX)
                out[j] = np.copysign(brentq(f, 0.0, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=500), w[j]) if w[j] != 0 else 0.0
            else:
                if w[j] == 0:
                    raise ValueError("alpha_of_w: w = 0 is a branch's end point (β = 0), not a point of the curve")
                s = np.sign(w[j])
                f = lambda u: float(self.w(a0[j] + np.exp(u), q[j])) - abs(w[j])
                lo, hi = -5.0, 0.0
                while f(lo) > 0:
                    if lo <= -740:
                        raise ValueError("alpha_of_w: w closer to the branch point than double precision")
                    lo = max(2 * lo, -740.0)
                while f(hi) < 0:
                    if hi >= np.log(ALPHA_MAX):
                        raise ValueError(f"alpha_of_w: w = {w[j]:.6g} beyond |α| = {ALPHA_MAX}")
                    hi = min(hi + 1.0, np.log(ALPHA_MAX))
                out[j] = s * (a0[j] + np.exp(brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)))
        return out.reshape(shp)

    # --- the loop on a dictionary ------------------------------------------------------------------------------------------------
    def rhs(self, Phi, y, xi=None, conserved=False):
        """The loop on (α, log β) on the design Phi (N × d), measured error r + ξ(t); `conserved` adds the replicator's −λ."""
        N, d = Phi.shape
        G, gam, m = self.G, self.gamma, self.m

        def f(t, z):
            al, lb = z[:d], z[d:]
            b = np.exp(lb); a, da = self.amp(al)
            r = y - Phi @ (b * a / m)
            if xi is not None:
                r = r + xi(t)
            g = Phi.T @ r / N
            u = a * g
            lam = (b @ u) / b.sum() if conserved else 0.0
            return np.r_[G * b * da * g, gam * (u - lam)]
        return f

    def rhs_curve(self, Phi, y, b0):
        """The loop from α = 0 with the masses b0 on the design Phi, the state α alone: β = b0 + 2c sinh²(α/a_max) keeps Q exactly
        (Proposition D‴.48, round 37 — integrating (α, log β) drifts off the curves at small masses), α̇ = G(b0 sech² + 2c tanh²)(α/a_max)
        ·Φᵀr/N."""
        N, d = Phi.shape
        b0 = np.asarray(b0, float)

        def f(t, al):
            x = np.minimum(np.abs(al) / self.A, 350.0)
            b = b0 + 2 * self.c * np.sinh(al / self.A) ** 2
            r = y - Phi @ (b * self.A * np.tanh(al / self.A) / self.m)
            return self.G * (b0 / np.cosh(x) ** 2 + 2 * self.c * np.tanh(x) ** 2) * (Phi.T @ r) / N
        return f

    def bregman_projection(self, Phi, y, q, eta0, sgn=None, tol=1e-12, it_max=4000):
        """argmin{Σψ_{q_i}(w_i) − ⟨η0, w⟩ : Φw = y} (each particle within its curve's range): w_i = w(α(η0_i + (Φᵀz)_i)) at the
        minimiser of the convex dual 𝒟(z) = Σψ*_{q_i}(η0_i + (Φᵀz)_i) − yᵀz, by `dual_newton`. Returns (w, α, iterations, the residual's
        norm |Φw − y| of the returned point — compare it with the tolerance: no silent failure)."""
        sgn = None if sgn is None else np.asarray(sgn, float)
        at = lambda e: self.alpha_of_eta(e, q, sgn)
        w, al, it, res = dual_newton(at, lambda al_: self.w(al_, q), lambda al_: self.m * self.relevance(al_, q),
                                     lambda al_: self.psi_star(al_, q), Phi, y, np.asarray(eta0, float),
                                     radius=2.0 / (self.gamma * self.A), tol=tol, it_max=it_max)
        return w, al, it, res

    # --- sampled data -------------------------------------------------------------------------------------------------------------
    def dQ_exact(self, al, beta, g, h):
        """Q's change under one Euler step on (α, log β) with the drive g."""
        a, da = self.amp(al)
        al1 = al + h * self.G * beta * da * g
        b1 = beta * np.exp(h * self.gamma * a * g)
        return self.Q(al1, b1) - self.Q(al, beta)

    def dQ_euler(self, al, beta, g2, h):
        """Its second-order part (h²g²γβ/2)(γa² − Gβa′(1 + a²/a_max²)) (g2 = g², or its mean)."""
        a, da = self.amp(al)
        return 0.5 * h ** 2 * g2 * self.gamma * beta * (self.gamma * a ** 2 - self.G * beta * da * (1 + (a / self.A) ** 2))

    def mirror_step(self, al, q, g, h, tol=1e-12, it_max=40):
        """The mirror update: η⁺ = η + hg, each particle placed on its own curve at η⁺ (Newton in α from the current α, the step
        capped at a_max and kept on the particle's side of a branch point; entries not converged to |η(α⁺) − η⁺| ≤ tol(1 + |η⁺|) are
        solved by `alpha_of_eta`); the mass is β(α⁺, q) — Q is kept exactly, and η moves by hg."""
        al = np.asarray(al, float)
        q = np.broadcast_to(np.asarray(q, float), al.shape)
        e1 = self.eta(al, q) + h * np.asarray(g, float)
        a0 = self.branch_point(q)
        s = np.where(al < 0, -1.0, 1.0)
        a = np.array(al, float)
        for _ in range(it_max):
            r = self.eta(a, q) - e1
            if np.all(np.abs(r) <= tol * (1 + np.abs(e1))):
                break
            step = np.clip(r * self.G * self.beta(a, q) * self.amp(a)[1], -self.A, self.A)
            a1 = a - step
            cross = (q + self.c <= 0) & (s * a1 <= a0)                         # never past the branch point
            a = np.where(cross, s * (a0 + 0.5 * (s * a - a0)), a1)
        bad = ~(np.abs(self.eta(a, q) - e1) <= tol * (1 + np.abs(e1)))
        if np.any(bad):
            a[bad] = self.alpha_of_eta(e1[bad], q[bad], s[bad])
        return a

    # --- the real loop --------------------------------------------------------------------------------------------------------------
    def real_rhs(self, centres, widths, wstar, x_ref, x_ref_dot, ke, conserved=False):
        """The real loop with Gaussian atoms at fixed centres and widths: z = (x, α, log β), e = x − x_m(t),
        ẋ = f(x) + ẋ_m − k_ee − f̂(x), f = Σw*_iφ_i, f̂ = Σw_iφ_i; α̇_i = Gβ_ia′_ieφ_i(x), (log β_i)˙ = γ(a_ieφ_i(x) − λ)."""
        centres, widths, wstar = (np.asarray(v, float) for v in (centres, widths, wstar))
        d = centres.size
        G, gam, m = self.G, self.gamma, self.m

        def f(t, z):
            x, al, lb = z[0], z[1:1 + d], z[1 + d:]
            b = np.exp(lb); a, da = self.amp(al)
            phi = np.exp(-(x - centres) ** 2 / (2 * widths ** 2))
            e = x - x_ref(t)
            fh = phi @ (b * a / m)
            u = a * e * phi
            lam = (b @ u) / b.sum() if conserved else 0.0
            return np.r_[phi @ wstar + x_ref_dot(t) - ke * e - fh, G * b * da * e * phi, gam * (u - lam)]
        return f

    def real_lyapunov(self, z, t, x_ref, al_star, q):
        """V = ½e² + Σ_i[w*_i(η*_i − η_i) − (β*_i − β_i)/(γm)] at the state z = (x, α, log β): β_i the state's own mass and η_i on the
        level Q(α_i, β_i) it defines (α*, η*, β* on the start's level q at w*), so V sees the mass law, not only the amplitudes."""
        d = q.size
        al, b = z[1:1 + d], np.exp(z[1 + d:])
        qs = self.Q(al, b)
        ws_ = self.w(al_star, q)
        B = ws_ * (self.eta(al_star, q) - self.eta(al, qs)) - (self.beta(al_star, q) - b) / (self.gamma * self.m)
        return 0.5 * (z[0] - x_ref(t)) ** 2 + float(np.sum(B))

    # --- the leak in the mirror coordinates (Proposition D‴.46) -----------------------------------------------------------------------
    def real_rhs_leak(self, centres, widths, f_true, x_ref, x_ref_dot, ke, q, eta_p, sigma, dist=None, noise=None):
        """The real loop with the leak −σ(η_i − η_p,i) added to each particle's scalar drive: z = (x, α, log β),
        ẋ = f(x) + ẋ_m − k_ee − f̂(x) + d(t); the adaptation sees e + n(t); α̇_i = Gβ_ia′_iu_i, (log β_i)˙ = γa_iu_i,
        u_i = (e + n)φ_i(x) − σ(η_i − η_p,i), η_i = η(α_i; Q(α_i, β_i)) with Q computed from the state (the law on (α, β) as implemented;
        q is the start's level, kept by the law). f_true is a function of x (realisable or not)."""
        centres, widths, eta_p = (np.asarray(v, float) for v in (centres, widths, eta_p))
        d = centres.size
        G, gam, m = self.G, self.gamma, self.m

        def f(t, z):
            x, al, lb = z[0], z[1:1 + d], z[1 + d:]
            b = np.exp(lb); a, da = self.amp(al)
            phi = np.exp(-(x - centres) ** 2 / (2 * widths ** 2))
            e = x - x_ref(t)
            em = e + (0.0 if noise is None else noise(t))
            u = em * phi - (sigma * (self.eta(al, self.Q(al, b)) - eta_p) if sigma else 0.0)
            xdot = f_true(x) + x_ref_dot(t) - ke * e - phi @ (b * a / m) + (0.0 if dist is None else dist(t))
            return np.r_[xdot, G * b * da * u, gam * a * u]
        return f

    def leak_equilibrium(self, H, wstar, q, eta_p, sigma, ke, dphi=None, tol=1e-12, it_max=500):
        """The averaged loop's equilibrium with the leak: the root of R(η) = (H(w(η) − w*) − d̄)/k_e + σ(η − η_p) — for symmetric H the
        minimiser of (1/(2k_e))(w − w*)ᵀH(w − w*) − ⟨d̄, w⟩/k_e + σΣB(w_i, w_p,i); H may be the error filter's H_L (`filtered_gram`), not
        symmetric in general, the equation monotone when H + Hᵀ ⪰ 0. Newton in η (Jacobian H diag(mD)/k_e + σI) with Armijo backtracking
        on |R|², trial steps capped in η by an adaptive radius, trial points beyond the curves rejected, the best iterate kept.
        Returns (w, α, iterations, |R| at the returned point)."""
        H, wstar, q, eta_p = (np.asarray(v, float) for v in (H, wstar, q, eta_p))
        dphi = np.zeros_like(wstar) if dphi is None else np.asarray(dphi, float)
        n = len(wstar)

        def resid(eta):
            al = self.alpha_of_eta(eta, q)
            w = self.w(al, q)
            return (H @ (w - wstar) - dphi) / ke + sigma * (eta - eta_p), al, w
        eta = eta_p.copy()
        R, al, w = resid(eta)
        rad = 2.0 / (self.gamma * self.A)
        it = 0
        for it in range(it_max):
            nR = np.linalg.norm(R)
            if nR < tol:
                break
            Jm = H * (self.m * self.relevance(al, q))[None, :] / ke + sigma * np.eye(n)
            step = np.linalg.lstsq(Jm, R, rcond=None)[0]
            de = np.abs(step).max()
            lam0 = min(1.0, rad / de) if de > 0 else 1.0
            lam, acc = lam0, None
            while lam > 1e-18 * lam0:
                try:
                    R2, al2, w2 = resid(eta - lam * step)
                except ValueError:
                    R2 = None
                if R2 is not None and np.all(np.isfinite(R2)) and np.linalg.norm(R2) ** 2 <= (1 - 1e-4 * lam) * nR ** 2:
                    acc = (eta - lam * step, R2, al2, w2)
                    break
                lam /= 2
            if acc is None:
                break
            rad = min(4 * rad, 1e15) if lam == lam0 else max(lam * de, 1e-12)
            eta, R, al, w = acc
        return w, al, it, float(np.linalg.norm(R))


# --- the shared dual solver ---------------------------------------------------------------------------------------------------------
def _exact_newton_step(Phi, lmD, gr, radius, cap_digits=3000):
    """Newton's step for the Hessian Φ diag(e^lmD)Φᵀ solved in decimal arithmetic (precision 40 digits beyond lmD's range, at most
    cap_digits; exponents unbounded), then scaled so that its largest change of η is at most radius: the curvatures of particles near
    extinction, far below the float Hessian's rounding, still steer the step. Raises ZeroDivisionError if the decimal system is singular."""
    from decimal import Decimal, localcontext
    N, d = Phi.shape
    top = float(np.max(lmD))
    with localcontext() as ctx:
        ctx.prec = int(min(cap_digits, 40 + (top - float(np.min(lmD))) / np.log(10)))
        ctx.Emin, ctx.Emax = -10 ** 9, 10 ** 9
        wts = [(Decimal(float(v)) - Decimal(top)).exp() for v in lmD]                # e^{lmD − top}
        Pd = [[Decimal(float(v)) for v in row] for row in Phi]
        M = [[sum((wts[k] * Pd[i][k] * Pd[j][k] for k in range(d)), Decimal(0)) for j in range(N)] + [-Decimal(float(gr[i]))]
             for i in range(N)]
        for col in range(N):                                                    # Gaussian elimination, partial pivoting
            piv = max(range(col, N), key=lambda r_: abs(M[r_][col]))
            if M[piv][col] == 0:
                raise ZeroDivisionError("singular Newton system")
            M[col], M[piv] = M[piv], M[col]
            for r_ in range(col + 1, N):
                f = M[r_][col] / M[col][col]
                for cc in range(col, N + 1):
                    M[r_][cc] -= f * M[col][cc]
        st = [Decimal(0)] * N
        for i in range(N - 1, -1, -1):
            st[i] = (M[i][N] - sum((M[i][j] * st[j] for j in range(i + 1, N)), Decimal(0))) / M[i][i]
        sc = (-Decimal(top)).exp()                                              # undo the scaling e^{−top}
        st = [v * sc for v in st]
        de = max(abs(sum((Pd[i][k] * st[i] for i in range(N)), Decimal(0))) for k in range(d))
        fac = min(Decimal(1), Decimal(float(radius)) / de) if de > 0 else Decimal(1)
        return np.array([float(v * fac) for v in st])


def dual_newton(inv, w_of, mD_of, psistar_of, Phi, y, eta0, radius=20.0, tol=1e-12, it_max=4000, logmD_of=None, exact_steps=300,
                exact_seconds=600.0):
    """Minimise the convex dual 𝒟(z) = Σψ*_i(η0_i + (Φᵀz)_i) − yᵀz (the Bregman projection's): Newton's direction (Hessian
    Φ diag(mD)Φᵀ), Armijo backtracking on 𝒟 — with, near the solution, the residual's halving accepted instead (𝒟's rounding) —, an
    adaptive radius on the step's largest change of η (grown fourfold after a full step, shrunk to the accepted step otherwise), trial
    points the inverse rejects (ValueError) or that are not finite skipped; the best iterate (smallest residual) is returned.
    inv(η) → state, w_of(state), mD_of(state) = dw/dη, psistar_of(state). With logmD_of(state) = log(dw/dη) given, a float step that
    stalls (accepted without measurable progress, or none accepted) switches the remaining steps to `_exact_newton_step` — needed when
    some particles are near extinction and others not, their curvatures below the float Hessian's rounding —, and a stalled exact step
    ends the iteration (the rounding floor); without it the iteration is the float one. The exact steps are budgeted (exact_steps, and
    exact_seconds of wall time — a step at thousands of digits takes seconds; round 36): past the budget the best iterate is returned, its
    residual telling the caller. Returns (w, state, iterations, |Φw − y|)."""
    z = np.zeros(Phi.shape[0])
    tol_abs = tol * (1 + np.linalg.norm(y))

    def at(z):
        with np.errstate(over='ignore', invalid='ignore'):
            st = inv(eta0 + Phi.T @ z)
            w = w_of(st)
            return st, w, float(np.sum(psistar_of(st)) - y @ z)
    st, w, Dz = at(z)
    best = (w, st, np.linalg.norm(Phi @ w - y))
    R, it, exact = radius, 0, False
    import time as _time
    n_exact, t_exact = 0, None
    for it in range(it_max):
        gr = Phi @ w - y
        ng = np.linalg.norm(gr)
        if ng < tol_abs:
            break
        step = None
        if exact:
            t_exact = _time.time() if t_exact is None else t_exact
            n_exact += 1
            if n_exact > exact_steps or _time.time() - t_exact > exact_seconds:
                break
            try:
                step = _exact_newton_step(Phi, logmD_of(st), gr, R)
            except ZeroDivisionError:
                step = None
        if step is None:
            Hs = Phi @ (mD_of(st)[:, None] * Phi.T)
            step = -np.linalg.lstsq(Hs, gr, rcond=None)[0]
        de = np.abs(Phi.T @ step).max()
        lam0 = min(1.0, R / de) if de > 0 else 1.0
        lam, slope, acc = lam0, gr @ step, None
        while lam > 1e-18 * lam0:
            try:
                st2, w2, D2 = at(z + lam * step)
            except ValueError:
                st2 = None
            if st2 is not None and np.all(np.isfinite(w2)) and np.isfinite(D2):
                with np.errstate(over='ignore', invalid='ignore'):
                    ng2 = np.linalg.norm(Phi @ w2 - y)
                if D2 <= Dz + 1e-4 * lam * slope or ng2 < 0.5 * ng:
                    acc = (z + lam * step, st2, w2, D2, ng2)
                    break
            lam /= 2
        if acc is None:
            if logmD_of is not None and not exact:
                exact = True
                continue
            break
        stalled = logmD_of is not None and not (acc[4] < 0.99 * ng or acc[3] < Dz - 1e-13 * max(1.0, abs(Dz)))
        R = min(4 * R, 1e15) if lam == lam0 else max(lam * de, 1e-12)
        z, st, w, Dz, ng2 = acc
        if ng2 < best[2]:
            best = (w, st, ng2)
        if stalled:
            if exact:                                                         # the exact step stalls too: the rounding floor
                break
            exact = True
    return best[0], best[1], it, float(best[2])


# --- the error filter's Gram (Proposition D‴.46(ii)) ----------------------------------------------------------------------------------
def filtered_gram(Pt, T, ke, Dt=None):
    """H_L = k_e⟨φL[φ]ᵀ⟩ and d̄_L = k_e⟨φL[d]⟩ over one period T of a periodic regressor, L = 1/(s + k_e) taken at its periodic steady
    state (by FFT: the k-th coefficient divided by iω_k + k_e). Pt: n × d samples of φ on a uniform grid of the period (endpoint
    excluded), Dt: n samples of the disturbance. As k_e → ∞, H_L → H = ⟨φφᵀ⟩ and d̄_L → ⟨φd⟩."""
    Pt = np.asarray(Pt, float)
    n = Pt.shape[0]
    om = 2 * np.pi * np.fft.fftfreq(n, d=T / n)
    filt = lambda X: np.real(np.fft.ifft(np.fft.fft(X, axis=0) / (1j * om + ke)[:, None], axis=0))
    LP = filt(Pt)
    H = ke * Pt.T @ LP / n
    if Dt is None:
        return H
    LD = filt(np.asarray(Dt, float)[:, None])[:, 0]
    return H, ke * Pt.T @ LD / n


# --- sparse references ---------------------------------------------------------------------------------------------------------------
def lasso(H, wstar, mu, ke, it=20000, tol=1e-14):
    """min (1/(2k_e))(w − w*)ᵀH(w − w*) + Σ_iμ_i|w_i| by cyclic coordinate descent (H positive definite)."""
    H, wstar, mu = (np.asarray(v, float) for v in (H, wstar, mu))
    w = np.zeros_like(wstar)
    b = H @ wstar
    for _ in range(it):
        w_old = w.copy()
        for i in range(len(w)):
            r = b[i] - H[i] @ w + H[i, i] * w[i]
            w[i] = np.sign(r) * max(abs(r) - ke * mu[i], 0.0) / H[i, i]
        if np.abs(w - w_old).max() < tol:
            break
    return w


def irrepresentable(H, supp, signs, weights=None):
    """The LASSO's irrepresentable quantity max_{j∉S}|(H_{S^cS}H_SS⁻¹(λ_S⊙sign(w*_S)))_j|/λ_j for the penalty weights λ (default equal:
    ‖H_{S^cS}H_SS⁻¹sign(w*_S)‖_∞): < 1 suffices for the support to be recovered at small penalties, ≤ 1 is necessary (round 36: weighted
    when the starting masses differ)."""
    H = np.asarray(H, float); supp = np.asarray(supp, int)
    lam = np.ones(H.shape[0]) if weights is None else np.asarray(weights, float)
    off = np.setdiff1d(np.arange(H.shape[0]), supp)
    if off.size == 0:
        return 0.0
    v = H[np.ix_(off, supp)] @ np.linalg.solve(H[np.ix_(supp, supp)], lam[supp] * np.asarray(signs, float))
    return float(np.abs(v / lam[off]).max())


def weighted_bp(Phi, y, weights, signs=None):
    """min Σ_iweights_i|w_i| s.t. Φw = y (and s_iw_i ≥ 0 where signs_i ≠ 0), by linear programming (w = p − n); NaN if infeasible."""
    d = Phi.shape[1]
    weights = np.asarray(weights, float)
    signs = np.zeros(d) if signs is None else np.asarray(signs, float)
    bnd = [(0, None if s_ >= 0 else 0.0) for s_ in signs] + [(0, None if s_ <= 0 else 0.0) for s_ in signs]
    r = linprog(np.concatenate([weights, weights]), A_eq=np.hstack([Phi, -Phi]), b_eq=y, bounds=bnd, method='highs')
    return (r.x[:d] - r.x[d:]) if r.status == 0 else np.full(d, np.nan)


def nonneg_certificate(Phi, supp):
    """The largest margin μ of a certificate h = Φᵀv (|v|_∞ ≤ 1), h = 0 on supp and h ≥ μ off it: μ > 0 makes the nonnegative fit
    supported on supp the only nonnegative fit. Returns (μ, v)."""
    N, d = Phi.shape
    supp = np.asarray(supp, int)
    off = np.setdiff1d(np.arange(d), supp)
    cobj = np.zeros(N + 1); cobj[-1] = -1
    A_ub = np.hstack([-Phi[:, off].T, np.ones((off.size, 1))])
    A_eq = np.hstack([Phi[:, supp].T, np.zeros((supp.size, 1))])
    r = linprog(cobj, A_ub=A_ub, b_ub=np.zeros(off.size), A_eq=A_eq, b_eq=np.zeros(supp.size),
                bounds=[(-1, 1)] * N + [(None, 1)], method='highs')
    return (float(-r.fun), r.x[:N]) if r.status == 0 else (float('nan'), None)


# --- the branches (Proposition D‴.47) -------------------------------------------------------------------------------------------------
def certificate_slacks(Phi, w_lim, S, rho, frozen=None):
    """A nondegenerate vertex w_lim of weighted basis pursuit (support S, |S| = N, Φ_S invertible): its dual certificate ξ, Φ_Sᵀξ =
    ρ_S·sign(w_lim,S), and the relative slacks μ_j = 1 − s_j(Φᵀξ)_j/ρ_j (s_j the frozen sign where `frozen` gives one, else
    sign((Φᵀξ)_j)) — the exponents of the rich limits' approach (Propositions D‴.45(iii), D‴.47(iii)). Returns (μ, ξ)."""
    S = np.asarray(S, int)
    rho = np.asarray(rho, float)
    xi = np.linalg.solve(Phi[:, S].T, rho[S] * np.sign(w_lim[S]))
    h = Phi.T @ xi
    sj = np.sign(h) if frozen is None else np.where(np.asarray(frozen) != 0, frozen, np.sign(h))
    return 1.0 - sj * h / rho, xi


def rich_leading_order(P, Phi, y, w_lim, S):
    """The rich limit's leading order at a nondegenerate vertex w_lim (support S, |S| = N; Propositions D‴.45(iii), D‴.47(iii)) for the
    particles P (a ParticleSet) started at their x0: the dual z from the support's mirror coordinates at the vertex, Φ_Sᵀz = η_S(w_lim) −
    η_S(0); the off-support coefficients through each particle's own mirror map, w_j = w(η_j(0) + (Φᵀz)_j); the support's from the fit
    Φ_S w_S = y − Φ_{S^c}w_{S^c}. Returns w₁."""
    S = np.asarray(S, int)
    e0 = P.eta(P.x0)
    z = np.linalg.solve(Phi[:, S].T, P.eta(P.state_of_w(w_lim, S))[S] - e0[S])
    w1 = P.w(P.x_of_eta(e0 + Phi.T @ z))
    w1[S] = 0.0
    w1[S] = np.linalg.solve(Phi[:, S], y - Phi @ w1)
    return w1


def iss_path(Phi, y, rho, signs=None, tol=1e-10, max_events=1000):
    """The inverse scale space flow of R(w) = Σ_iρ_i|w_i| (frozen signs s_i: + the indicator of s_iw_i ≥ 0) — the rich transit's limit
    (Proposition D‴.48): du/dτ = Φᵀ(y − Φw)/N, u(0) = 0, w in the normal cone of the box |u_i| ≤ ρ_i (one-sided: s_iu_i ≤ ρ_i), the
    least-squares element of it (the minimal-norm section of the dual's subdifferential): w piecewise constant — on the wall's atoms the
    sign-constrained least-squares fit (NNLS after the signs) —, u piecewise linear, an atom entering when its u reaches the wall and
    leaving when its coefficient vanishes with an inward velocity (the adaptive inverse scale space method of Burger–Möller–Benning–Osher).
    Returns dict(taus = the event times, U = u there, W = w on each segment, events = [(τ, i, 'enter'/'exit')], w_end)."""
    from scipy.optimize import nnls
    N, d = Phi.shape
    rho = np.asarray(rho, float)
    one = signs is not None
    s = np.asarray(signs, float) if one else None
    u, tau = np.zeros(d), 0.0
    on = np.zeros(d, bool)
    taus, Us, Ws, ev = [0.0], [u.copy()], [], []
    w = np.zeros(d)
    for _ in range(max_events):
        idx = np.flatnonzero(on)
        w = np.zeros(d)
        if idx.size:
            sg = s[idx] if one else np.sign(u[idx])
            w[idx] = sg * nnls(Phi[:, idx] * sg, y, maxiter=50 * d)[0]
        r = y - Phi @ w
        v = Phi.T @ r / N
        Ws.append(w.copy())
        sg_all = s if one else np.sign(u)
        leave = on & (w == 0) & (sg_all * v < -tol * max(1.0, np.abs(v).max()))
        ev += [(tau, int(i), 'exit') for i in np.flatnonzero(leave)]
        on &= ~leave
        if np.linalg.norm(r) <= tol * (1 + np.linalg.norm(y)):
            break
        dt = np.full(d, np.inf)
        off = ~on
        if one:
            up = off & (s * v > 0)
            dt[up] = (rho[up] - s[up] * u[up]) / (s[up] * v[up])
        else:
            pos, neg = off & (v > 0), off & (v < 0)
            dt[pos] = (rho[pos] - u[pos]) / v[pos]
            dt[neg] = (-rho[neg] - u[neg]) / v[neg]
        j = int(np.argmin(dt))
        if not np.isfinite(dt[j]):
            break
        u = u + max(dt[j], 0.0) * v
        tau += max(dt[j], 0.0)
        u[j] = (s[j] if one else np.sign(v[j])) * rho[j]
        on[j] = True
        ev.append((tau, j, 'enter'))
        taus.append(tau); Us.append(u.copy())
    return dict(taus=np.array(taus), U=np.array(Us), W=np.array(Ws), events=ev, w_end=w)


def iss_eval(path, tau, Phi, y):
    """u(τ) and w(τ) on an `iss_path` (w of the segment containing τ)."""
    k = min(int(np.searchsorted(path['taus'], tau, side='right') - 1), len(path['W']) - 1)
    w = path['W'][k]
    return path['U'][k] + (tau - path['taus'][k]) * (Phi.T @ (y - Phi @ w)) / Phi.shape[0], w


class ParticleSet:
    """Particles of one free law from starts (α(0), β(0)), each on its curve Q_i: a free sign (Q_i > −c, the state α_i) or a branch
    (β(0) < 2c sinh²(α(0)/a_max), i.e. Q_i < −c; the state u_i = log δ_i, δ_i = |α_i| − α₀,i, the sign s_i frozen). On a branch every
    quantity is computed from δ without cancellation — β = 2c sinh((2α₀ + δ)/a_max)sinh(δ/a_max), η from u (`FreeLaw._eta_branch_log`)
    —, so a particle near extinction (β → 0, s·η → −∞) is represented to any depth (its β underflowing to 0 while u keeps its value); the start's δ(0) is solved from β(0) directly (no subtraction
    Q = β − c cosh(2α/a_max)), so starting masses far below c's precision are kept."""

    def __init__(self, law, al_init, b_init):
        L = self.law = law
        A, c = L.A, L.c
        al_init = np.array(al_init, float, ndmin=1)
        b_init = np.broadcast_to(np.asarray(b_init, float), al_init.shape).copy()
        self.d = al_init.size
        self.br = b_init < 2 * c * np.sinh(al_init / A) ** 2
        self.s = np.where(self.br, np.sign(al_init), 0.0)
        self.al0 = np.zeros(self.d)
        x0 = al_init.copy()
        q = np.empty(self.d)
        for i in range(self.d):
            if not self.br[i]:
                q[i] = b_init[i] - c * np.cosh(2 * al_init[i] / A)
                continue
            ab = abs(al_init[i])
            f = lambda dl: np.log(2 * c * np.sinh((2 * ab - dl) / A)) + logsinh(dl / A) - np.log(b_init[i])
            hi = ab * (1 - 1e-12)
            lo = min(b_init[i] * A / (2 * c * np.sinh(2 * ab / A)), 0.5 * hi)
            while f(lo) > 0:
                lo *= 0.5
            dl = brentq(f, lo, hi, xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=500) if f(hi) > 0 else hi
            self.al0[i] = ab - dl
            x0[i] = np.log(dl)
            q[i] = -c * np.cosh(2 * self.al0[i] / A)
        self.q, self.x0 = q, x0

    def _dl(self, x):
        return np.exp(np.where(self.br, x, 0.0))

    def alpha(self, x):
        return np.where(self.br, self.s * (self.al0 + self._dl(x)), x)

    def beta(self, x):
        L = self.law
        dl = self._dl(x)
        bb = 2 * L.c * np.sinh((2 * self.al0 + dl) / L.A) * np.sinh(dl / L.A)
        return np.where(self.br, bb, self.q + L.c * np.cosh(2 * self.alpha(x) / L.A))

    def w(self, x):
        L = self.law
        return self.beta(x) * L.A * np.tanh(self.alpha(x) / L.A) / L.m

    def mD(self, x):
        """dw/dη."""
        L = self.law
        b = self.beta(x); a, da = L.amp(self.alpha(x))
        return (b / L.m) * (L.G * b * da ** 2 + L.gamma * a ** 2)

    def log_mD(self, x):
        """log(dw/dη) without underflow: on a branch log β from u = log δ (β = 2c sinh((2α₀ + δ)/a_max)sinh(δ/a_max))."""
        L = self.law
        a, da = L.amp(self.alpha(x))
        dl = self._dl(x)
        xb = np.where(self.br, x, 0.0)
        tiny = dl / L.A < 1e-8
        lsh = np.where(tiny, xb - np.log(L.A) + dl ** 2 / (6 * L.A ** 2), logsinh(np.where(tiny, 1.0, dl) / L.A))
        with np.errstate(divide='ignore'):
            lb = np.where(self.br, np.log(2 * L.c) + logsinh((2 * self.al0 + dl) / L.A) + lsh, np.log(np.abs(self.beta(x))))
        return lb - np.log(L.m) + np.log(L.G * np.exp(lb) * da ** 2 + L.gamma * a ** 2)

    def temperature(self, x):
        """dη/dlog|w| = a/(γa² + Gβa′²): 1/(γa₀) at a branch point, 1/(γa_max) far out."""
        L = self.law
        b = self.beta(x); a, da = L.amp(self.alpha(x))
        return a / (L.gamma * a ** 2 + L.G * b * da ** 2)

    def eta(self, x):
        L = self.law
        out = np.empty(self.d)
        fr = ~self.br
        if np.any(fr):
            out[fr] = L.eta(x[fr], self.q[fr])
        if np.any(self.br):
            out[self.br] = self.s[self.br] * L._eta_branch_log(x[self.br], self.al0[self.br])
        return out

    def psistar(self, x):
        """ψ*(η) = 2c sinh²(α/a_max)/(γm) (up to its constant on a branch)."""
        L = self.law
        return 2 * L.c * np.sinh(self.alpha(x) / L.A) ** 2 / (L.gamma * L.m)

    def bdiv(self, wbar, xbar, x):
        """B(w̄, w) = Σ_i[w̄_i(η̄_i − η_i) − (β̄_i − β_i)/(γm)], the terms with w̄_i = 0 on a branch read as β_i/(γm) (xbar unused there)."""
        L = self.law
        wbar = np.asarray(wbar, float)
        eb, e = self.eta(xbar), self.eta(x)
        bb, b = self.beta(xbar), self.beta(x)
        zero = self.br & (wbar == 0)
        terms = np.where(zero, b / (L.gamma * L.m), wbar * (eb - e) - (bb - b) / (L.gamma * L.m))
        return float(np.sum(terms))

    def x_of_eta(self, e):
        """The inverse of η: `FreeLaw.alpha_of_eta` for free signs, brentq in u = log δ on a branch (raises beyond the range)."""
        L = self.law
        e = np.asarray(e, float)
        out = np.empty(self.d)
        fr = ~self.br
        if np.any(fr):
            out[fr] = L.alpha_of_eta(e[fr], self.q[fr])
        for j in np.flatnonzero(self.br):
            tgt = self.s[j] * e[j]
            f = lambda u: float(L._eta_branch_log(u, self.al0[j])) - tgt
            lo, hi = -5.0, 0.0
            while f(lo) > 0:
                if lo <= -1e300:
                    raise ValueError("x_of_eta: η beyond any depth")
                lo *= 2
            while f(hi) < 0:
                if hi >= np.log(ALPHA_MAX):
                    raise ValueError("x_of_eta: beyond |α| = ALPHA_MAX")
                hi = min(hi + 0.5, np.log(ALPHA_MAX))
            out[j] = brentq(f, lo, hi, xtol=1e-12, rtol=4 * np.finfo(float).eps, maxiter=500)
        return out

    def state_of_w(self, wv, idx):
        """The state with w_i = wv_i for i in idx (the others at the start), by brentq per coordinate."""
        xs_ = self.x0.copy()
        for i in idx:
            f = lambda xx: float(self.w(np.where(np.arange(self.d) == i, xx, self.x0))[i]) - wv[i]
            lo, hi = (-740.0, np.log(ALPHA_MAX)) if self.br[i] else (-ALPHA_MAX, ALPHA_MAX)
            xs_[i] = brentq(f, lo, hi, xtol=1e-14, rtol=4 * np.finfo(float).eps, maxiter=500)
        return xs_

    def _dstate(self, x, g):
        """ẋ for the drive g: α̇ = Gβa′g (free), u̇ = sGβa′g/δ with β/δ = 2c sinh((2α₀ + δ)/a_max)·sinh(δ/a_max)/δ (branch)."""
        L = self.law
        a, da = L.amp(self.alpha(x))
        dl = self._dl(x)
        shr = np.where(dl / L.A < 1e-8, 1.0 / L.A, np.sinh(dl / L.A) / np.where(dl > 0, dl, 1.0))
        bod = 2 * L.c * np.sinh((2 * self.al0 + dl) / L.A) * shr
        return np.where(self.br, self.s * L.G * bod * da * g, L.G * self.beta(x) * da * g)

    def rhs(self, Phi, y):
        """The noiseless loop on the design: dη = −∇J/m dt, J = (m/2N)|Φw − y|², i.e. g = Φᵀ(y − Φw)/N."""
        N = Phi.shape[0]
        return lambda t, x: self._dstate(x, -(Phi.T @ (Phi @ self.w(x) - y)) / N)

    def real_rhs(self, centres, widths, wstar, x_ref, x_ref_dot, ke):
        """The real loop (Proposition D‴.45(v)): z = (x, state), e = x − x_m(t), ẋ = Σw*_iφ_i(x) + ẋ_m − k_ee − Σw_iφ_i(x), g_i = eφ_i(x)."""
        centres, widths, wstar = (np.asarray(v, float) for v in (centres, widths, wstar))

        def f(t, z):
            xp, st = z[0], z[1:]
            phi = np.exp(-(xp - centres) ** 2 / (2 * widths ** 2))
            e = xp - x_ref(t)
            return np.r_[phi @ wstar + x_ref_dot(t) - ke * e - phi @ self.w(st), self._dstate(st, e * phi)]
        return f


def sign_ls_set(Phi, y, s, tol=1e-9):
    """The sign-constrained least-squares fits S = argmin{|Φw − y|² : s_iw_i ≥ 0 where s_i ≠ 0}: ŷ (Φw on S, unique), Z (the constrained
    coordinates vanishing on all of S, by a linear program per coordinate) and one point of S. The point by Lawson–Hanson's active set (nnls
    on the columns s_iφ_i and ±φ_j for the free coordinates) — exact whatever the rank —, its KKT conditions then checked (round 36:
    lsq_linear's status was never read, and it failed silently on rank-deficient and duplicated designs); raises when they fail or a
    linear program does not solve."""
    from scipy.optimize import nnls
    N, d = Phi.shape
    s = np.asarray(s, float)
    con, fre = np.flatnonzero(s != 0), np.flatnonzero(s == 0)
    B = np.hstack([Phi[:, con] * s[con], Phi[:, fre], -Phi[:, fre]])
    x = nnls(B, y, maxiter=100 * max(B.shape[1], 1))[0]
    w = np.zeros(d)
    w[con] = s[con] * x[:con.size]
    w[fre] = x[con.size:con.size + fre.size] - x[con.size + fre.size:]
    yh = Phi @ w
    g = Phi.T @ (yh - y)                                                          # ∇½|Φw − y|²
    scale = 1.0 + np.abs(Phi).max() * (np.linalg.norm(y) + np.linalg.norm(yh))
    on = np.abs(w) > tol * (1.0 + np.abs(w).max())
    viol = max(np.abs(g[fre]).max(initial=0.0), np.abs(g[con[on[con]]]).max(initial=0.0),
               (-s[con] * g[con])[~on[con]].max(initial=0.0))
    if viol > 1e3 * tol * scale:
        raise RuntimeError(f'sign_ls_set: KKT violated by {viol:.3g}')
    bounds = [(0, None) if s[j] > 0 else ((None, 0) if s[j] < 0 else (None, None)) for j in range(d)]
    Z = []
    for i in con:
        cc = np.zeros(d); cc[i] = -s[i]
        res = linprog(cc, A_eq=Phi, b_eq=yh, bounds=bounds, method='highs')
        if res.status == 3:                                                       # s_iw_i unbounded on S: not in Z
            continue
        if res.status != 0:
            raise RuntimeError(f'sign_ls_set: the linear program for coordinate {i} failed ({res.message})')
        if -res.fun < tol * (1 + np.abs(w).max()):
            Z.append(int(i))
    return yh, np.array(Z, int), w


def projection(P, Phi, yh, eta0, keep=None, tol=1e-12, it_max=4000, exact=False):
    """The Bregman projection argmin{Σ_{i∈keep}[ψ_i(w_i) − η0_iw_i] : Φ_keep w_keep = ŷ} for a ParticleSet (the others 0, each
    particle's range — a branch's sign — implied), by `dual_newton` on the kept coordinates (exact=True: with its decimal fallback,
    from `ParticleSet.log_mD`, for masses far below the float Hessian's range). Returns (w, iterations, the residual)."""
    keep = np.arange(P.d) if keep is None else np.asarray(keep, int)
    Pk = Phi[:, keep]
    eta0 = np.asarray(eta0, float)

    def inv(e_k):
        e = eta0.copy(); e[keep] = e_k
        return P.x_of_eta(e)
    wk, st, it, res = dual_newton(inv, lambda st_: P.w(st_)[keep], lambda st_: P.mD(st_)[keep],
                                  lambda st_: P.psistar(st_)[keep], Pk, yh, eta0[keep],
                                  radius=2.0 / (P.law.gamma * P.law.A), tol=tol, it_max=it_max,
                                  logmD_of=(lambda st_: P.log_mD(st_)[keep]) if exact else None)
    out = np.zeros(P.d); out[keep] = wk
    return out, it, res
