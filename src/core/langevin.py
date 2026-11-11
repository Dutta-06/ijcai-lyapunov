"""
Langevin adaptive control (PARADIGM.md §6c, Theorem H)
======================================================

Saturated lifted features + σ-modification (confinement) + parameter dither (entropy).

A particle is an UNCONSTRAINED vector ϑ = (α, γ, κ) ∈ ℝⁿ × ℝⁿ × ℝ, mapped smoothly into the box
by per-coordinate saturations

    a = a_max tanh(α / a_max),   c = c_max tanh(γ / c_max),   s = s_min + (s_max − s_min) sig(κ),

and the feature is Φ̃(x; ϑ) = a φ(x; c, s), so ‖Φ̃‖ ≤ K_Φ := a_max √n and every derivative of Φ̃ in ϑ
is bounded uniformly in x (s ≥ s_min > 0).  No projection is needed: the feature is bounded by
construction, and the parameters are confined by the leak.

Adaptive law (Euler–Maruyama in simulation):

    dϑ_i = [ G ∇_ϑ g_t(ϑ_i) − σ (ϑ_i − ϑ⁰) ] dt + √(2τ) G^{1/2} dB_i,    g_t(ϑ) = eᵀ Φ̃(x; ϑ),

i.e. the classical gradient law, σ-modification written in the G⁻¹ metric (velocity −σ(ϑ − ϑ⁰)),
and dither of temperature τ.  In the averaged mean-field limit this is the mean-field Langevin dynamics
(in the G⁻¹ metric) of

    F_τ(μ) = (1/k_e) L̄(μ) + (σ/2) ∫ ‖ϑ − ϑ⁰‖²_{G⁻¹} dμ + τ H(μ),     L̄(μ) = ½ ∫ ‖f − ∫Φ̃ dμ‖² dρ,

with ρ the occupation measure of the reference (Theorem H).
"""
import numpy as np
from scipy.special import expit as _sig        # numerically stable logistic function


# NOTE: occupation_loss and first_variation below assume an AUTONOMOUS plant f = f(x) (they evaluate
# plant.f(x, 0.0)); Theorem H is stated for autonomous plants (PARADIGM.md §6c).


class SaturatedLiftedNetwork:
    def __init__(self, n, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0):
        self.n = int(n)
        self.a_max, self.c_max, self.s_min, self.s_max = map(float, (a_max, c_max, s_min, s_max))

    # ---------------------------------------------------------------- maps
    def physical(self, Th):
        """Θ (m, p) unconstrained → (a (m,n), c (m,n), s (m,))."""
        n = self.n
        a = self.a_max * np.tanh(Th[:, :n] / self.a_max)
        c = self.c_max * np.tanh(Th[:, n:2 * n] / self.c_max)
        s = self.s_min + (self.s_max - self.s_min) * _sig(Th[:, 2 * n])
        return a, c, s

    def unconstrained(self, a, c, s):
        """Inverse map (for initialisation)."""
        al = self.a_max * np.arctanh(np.clip(a / self.a_max, -1 + 1e-12, 1 - 1e-12))
        ga = self.c_max * np.arctanh(np.clip(c / self.c_max, -1 + 1e-12, 1 - 1e-12))
        u = np.clip((s - self.s_min) / (self.s_max - self.s_min), 1e-12, 1 - 1e-12)
        return np.hstack([al, ga, np.log(u / (1 - u))[:, None]])

    @property
    def K_Phi(self):
        return self.a_max * np.sqrt(self.n)

    # ---------------------------------------------------------------- features
    def f_hat(self, x, Th):
        a, c, s = self.physical(Th)
        phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
        return (phi[:, None] * a).mean(axis=0)

    def f_hat_batch(self, X, Th):
        """X (N, n) → (N, n)."""
        a, c, s = self.physical(Th)
        d2 = np.sum((X[:, None, :] - c[None, :, :]) ** 2, axis=2)          # (N, m)
        phi = np.exp(-d2 / (2 * s[None, :] ** 2))
        return phi @ a / Th.shape[0]

    def grad_g(self, x, e, Th):
        """∇_ϑ (eᵀ Φ̃(x; ϑ_i)) for every particle, (m, p)."""
        n = self.n
        a, c, s = self.physical(Th)
        diff = x[None, :] - c
        sq = np.sum(diff ** 2, axis=1)
        phi = np.exp(-sq / (2 * s ** 2))
        ae = a @ e
        da = 1.0 - (a / self.a_max) ** 2                                   # da/dα
        dc = 1.0 - (c / self.c_max) ** 2                                   # dc/dγ
        sk = _sig(Th[:, 2 * n]); ds = (self.s_max - self.s_min) * sk * (1 - sk)   # ds/dκ
        out = np.empty_like(Th)
        out[:, :n] = phi[:, None] * e[None, :] * da
        out[:, n:2 * n] = (ae * phi / s ** 2)[:, None] * diff * dc
        out[:, 2 * n] = ae * phi * sq / s ** 3 * ds
        return out


class LangevinAdaptiveController:
    """
    Parameters: n, ke, G_diag (p,), sigma (leak gain), theta0 (p,) leak target, tau (dither temperature),
    net: SaturatedLiftedNetwork.  tau_schedule(t) optionally overrides tau (annealing).
    """

    def __init__(self, net, ke=10.0, G_diag=None, sigma=0.05, theta0=None, tau=0.0, tau_schedule=None):
        self.net, self.ke = net, float(ke)
        p = 2 * net.n + 1
        self.G = np.ones(p) if G_diag is None else np.asarray(G_diag, float)
        self.sigma = float(sigma)
        self.theta0 = np.zeros(p) if theta0 is None else np.asarray(theta0, float)
        self.tau, self.tau_schedule = float(tau), tau_schedule

    def control(self, x, x_ref, x_ref_dot, Th):
        return x_ref_dot - self.ke * (x - x_ref) - self.net.f_hat(x, Th)

    def drift(self, x, e, Th):
        return self.G[None, :] * self.net.grad_g(x, e, Th) - self.sigma * (Th - self.theta0[None, :])

    def temperature(self, t):
        return self.tau if self.tau_schedule is None else float(self.tau_schedule(t))

    def lsi_constant_lower_bound(self, F_max, tau=None):
        """Holley–Stroock ∘ Bakry–Émery (whitened coordinates): ρ_τ ≥ (σ/τ) exp(−2K_Φ(F_max + K_Φ)/(k_e τ))."""
        tau = self.tau if tau is None else tau
        if tau <= 0 or self.sigma <= 0:
            raise ValueError("the LSI bound needs σ > 0 and τ > 0")
        K = self.net.K_Phi
        return (self.sigma / tau) * np.exp(-2 * K * (F_max + K) / (self.ke * tau))

    def log_rate_lower_bound(self, F_max, tau=None):
        """log of the guaranteed rate 2τρ_τ (avoids underflow at small τ)."""
        tau = self.tau if tau is None else tau
        if tau <= 0 or self.sigma <= 0:
            raise ValueError("needs σ > 0 and τ > 0")
        K = self.net.K_Phi
        return float(np.log(2 * self.sigma) - 2 * K * (F_max + K) / (self.ke * tau))


def simulate(plant, ctrl: LangevinAdaptiveController, x0, Th0, T, dt, rng, record_every=50):
    """
    Closed loop: plant by RK4 with u held over the step; parameters by Euler–Maruyama.
    Returns dict(t, e_norm, Th (snapshots), x).
    """
    x, Th, t = np.array(x0, float), np.array(Th0, float), 0.0
    m, p = Th.shape
    sqG = np.sqrt(ctrl.G)
    n_steps = int(round(T / dt))
    ts, es, Ths, xs = [], [], [], []
    for k in range(n_steps):
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        e = x - xr
        if k % record_every == 0:
            ts.append(t); es.append(float(np.linalg.norm(e))); Ths.append(Th.copy()); xs.append(x.copy())
        u = ctrl.control(x, xr, xrd, Th)
        k1 = plant.f(x, t) + u
        k2 = plant.f(x + 0.5 * dt * k1, t + 0.5 * dt) + u
        k3 = plant.f(x + 0.5 * dt * k2, t + 0.5 * dt) + u
        k4 = plant.f(x + dt * k3, t + dt) + u
        dTh = ctrl.drift(x, e, Th) * dt
        tau = ctrl.temperature(t)
        if tau > 0:
            dTh += np.sqrt(2 * tau * dt) * sqG[None, :] * rng.standard_normal((m, p))
        x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Th = Th + dTh
        t = (k + 1) * dt
    return dict(t=np.array(ts), e_norm=np.array(es), Th=np.array(Ths), x=np.array(xs))


def occupation_loss(plant, net, Th, period, n_quad=400):
    """L̄(μ) = ½⟨‖f(x) − f̂(x)‖²⟩_ρ over one period of the reference (ρ = occupation measure)."""
    ts = np.arange(n_quad) * period / n_quad
    X = np.array([plant.x_ref(t) for t in ts])
    F = np.array([plant.f(x, 0.0) for x in X])
    R = F - net.f_hat_batch(X, Th)
    return 0.5 * float(np.mean(np.sum(R ** 2, axis=1)))


def first_variation(plant, net, Th, theta_eval, period, n_quad=400):
    """V_L[μ](ϑ) = −⟨r_μ(x)ᵀ Φ̃(x; ϑ)⟩_ρ for each row of theta_eval."""
    ts = np.arange(n_quad) * period / n_quad
    X = np.array([plant.x_ref(t) for t in ts])
    F = np.array([plant.f(x, 0.0) for x in X])
    R = F - net.f_hat_batch(X, Th)                                        # (N, n)
    a, c, s = net.physical(theta_eval)
    d2 = np.sum((X[:, None, :] - c[None, :, :]) ** 2, axis=2)
    phi = np.exp(-d2 / (2 * s[None, :] ** 2))                            # (N, k)
    return -np.mean((R @ a.T) * phi, axis=0)                             # (k,)
