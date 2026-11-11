"""
Oracle diagnostics: exact evaluation of the Lyapunov argument
=============================================================

Everything in the manuscript's stability argument involves the *unknown* ideal
parameters θ* = (W*, c*, σ*).  Those quantities can therefore never be computed
on a real plant — the LyapunovMonitor is an oracle diagnostic, not an
observable.  To check the argument numerically we use a plant that is EXACTLY
representable by the network with a known, constant θ*:

    f(x) = W*ᵀ φ(x; c*, σ*)        (ε ≡ 0,  Ẇ* = ċ* = σ̇* = 0  ⇒  C_Ω = 0)

Along a simulated trajectory we then evaluate exactly (THEORY_REVISED.md):

    V      = ½ eᵀe + ½ vec(W̃)ᵀΓ_W⁻¹vec(W̃) + ½ c̃ᵀΓ_c⁻¹c̃ + ½ σ̃ᵀΓ_σ⁻¹σ̃
    V̇      = eᵀė + vec(W̃)ᵀΓ_W⁻¹vec(Ẇ) + c̃ᵀΓ_c⁻¹ċ + σ̃ᵀΓ_σ⁻¹σ̇
    d(t)   = f − f̂ + W̃ᵀφ + Wᵀφ'_c c̃ + Wᵀφ'_σ σ̃            (exact lumped disturbance)
    cross  = ⟨W̃, Γ_W⁻¹Ẇ − φeᵀ⟩ + c̃ᵀ(Γ_c⁻¹ċ − φ'_cᵀWe) + σ̃ᵀ(Γ_σ⁻¹σ̇ − φ'_σᵀWe)

with the two facts the proof rests on:

    identity:    V̇ = −kₑ‖e‖² + eᵀd + cross           (algebra of Sections 2–4)
    property P:  cross ≤ 0                              (projection, Section 5)

and the consequence that is actually testable without any unknown constant:

    RMS bound:   RMS_{[T₁,T₂]}(‖e‖) ≤ D/(2kₑ) + √(D²/(4kₑ²) + V(T₁)/(kₑ (T₂−T₁)))
                 with D = sup_{[T₁,T₂]} ‖d(t)‖  (measured).
"""

import numpy as np
from .plant_models import PlantModel
from .rbf_kernel import RBFKernel
from .integrator import AdaptiveRK4Integrator, SimulationState
from .projection import min_center_separation


class RepresentablePlant(PlantModel):
    """f(x) = W*ᵀφ(x; c*, σ*) with a smooth multi-sine reference."""

    def __init__(self, W_star, c_star, sigma_star, amplitude=1.0, omega=1.0):
        W_star = np.asarray(W_star, float)
        super().__init__(n=W_star.shape[1], name="Representable")
        self.W_star, self.c_star, self.sigma_star = W_star, np.asarray(c_star, float), np.asarray(sigma_star, float)
        self.kernel = RBFKernel(W_star.shape[0], self.n)
        self.A, self.w = float(amplitude), float(omega)
        self._freq = self.w * (1.0 + 0.3 * np.arange(self.n))
        self._phase = 0.5 * np.arange(self.n)

    def f(self, x, t):
        return self.kernel.f_hat(x, self.W_star, self.c_star, self.sigma_star)

    def x_ref(self, t):
        return self.A * np.sin(self._freq * t + self._phase)

    def x_ref_dot(self, t):
        return self.A * self._freq * np.cos(self._freq * t + self._phase)

    def F_max(self):
        """sup_x ‖f(x)‖ ≤ ‖W*‖₂ ‖φ‖ ≤ ‖W*‖₂ √m  (exact for this plant class)."""
        return np.linalg.norm(self.W_star, 2) * np.sqrt(self.W_star.shape[0])


def sample_admissible_ideal(m, n, W_max, c_max, sigma_min, sigma_max, epsilon,
                            rng, center_spread=1.2, min_sep=0.6, w_scale=0.8):
    """Random θ* strictly inside the inner projection sets (Assumption A3)."""
    c_star = rng.uniform(-center_spread, center_spread, size=(m, n))
    for _ in range(500):
        moved = False
        for i in range(m):
            for j in range(i + 1, m):
                dv = c_star[i] - c_star[j]
                dist = np.linalg.norm(dv)
                if dist < min_sep:
                    moved = True
                    step = 0.05 * dv / (dist + 1e-9) if dist > 1e-9 else 0.05 * rng.normal(size=n)
                    c_star[i] += step
                    c_star[j] -= step
        if not moved:
            break
    sigma_star = rng.uniform(sigma_min + 2 * epsilon, min(sigma_max - 2 * epsilon, 1.6), size=m)
    W_star = rng.normal(0, w_scale, size=(m, n))
    nW = np.linalg.norm(W_star, 'fro')
    if nW > 0.8 * (W_max - epsilon):
        W_star *= 0.8 * (W_max - epsilon) / nW
    nc = np.linalg.norm(c_star)
    if nc > 0.8 * (c_max - epsilon):
        c_star *= 0.8 * (c_max - epsilon) / nc
    return W_star, c_star, sigma_star


class ExactLyapunov:
    """Evaluates V, V̇, d, cross exactly for a controller and known constant θ*."""

    def __init__(self, ctrl, plant: RepresentablePlant):
        self.ctrl, self.plant = ctrl, plant
        self.W_star, self.c_star, self.s_star = plant.W_star, plant.c_star, plant.sigma_star

    def evaluate(self, state: SimulationState, t: float) -> dict:
        ctrl, plant, ker = self.ctrl, self.plant, self.ctrl.kernel
        x, W, c, s = state.x, state.W, state.c, state.sigma
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        e = x - xr
        u = ctrl.control_law(x, xr, xrd, W, c, s)
        edot = plant.f(x, t) + u - xrd
        Wd, cd, sd = ctrl.adaptation_laws(x, e, W, c, s)
        phi = ker.phi(x, c, s)
        phi_c = ker.phi_jacobian_c(x, c, s)
        phi_s_diag = ker.phi_jacobian_sigma_diag(x, c, s)
        Wt = W - self.W_star
        ct = (c - self.c_star).reshape(-1)
        st = s - self.s_star
        We = W @ e
        GWi, Gci, Gsi = ctrl.gamma_W_inv, ctrl.gamma_c_inv, ctrl.gamma_sigma_inv
        d = (plant.f(x, t) - ker.f_hat(x, W, c, s) + Wt.T @ phi
             + W.T @ (phi_c @ ct) + W.T @ (phi_s_diag * st))
        GWi_Wd = (GWi @ Wd.reshape(-1)).reshape(W.shape)
        cross = (float(np.sum(Wt * (GWi_Wd - np.outer(phi, e))))
                 + float(ct @ (Gci @ cd - phi_c.T @ We))
                 + float(st @ (Gsi @ sd - phi_s_diag * We)))
        Vdot = (float(e @ edot) + float(np.sum(Wt * GWi_Wd))
                + float(ct @ (Gci @ cd)) + float(st @ (Gsi @ sd)))
        V = (0.5 * float(e @ e) + 0.5 * float(Wt.reshape(-1) @ GWi @ Wt.reshape(-1))
             + 0.5 * float(ct @ Gci @ ct) + 0.5 * float(st @ Gsi @ st))
        return dict(t=t, V=V, Vdot=Vdot, rhs=-ctrl.ke * float(e @ e) + float(e @ d) + cross,
                    cross=cross, e_norm=float(np.linalg.norm(e)), d_norm=float(np.linalg.norm(d)),
                    min_sep=min_center_separation(c),
                    inside=ctrl.params_inside_sets(W, c, s))


def run_oracle(ctrl, plant, state0: SimulationState, T: float, dt: float = 1e-3, every: int = 5):
    """Simulate and return a dict of arrays with the exact Lyapunov quantities."""
    integ = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
    ex = ExactLyapunov(ctrl, plant)
    st = state0.copy()
    rows, t = [], 0.0
    n_steps = int(round(T / dt))
    for k in range(n_steps):
        if k % every == 0:
            rows.append(ex.evaluate(st, t))
        st = integ.step(st, t)
        t = (k + 1) * dt
    keys = rows[0].keys()
    return {k: np.array([r[k] for r in rows]) for k in keys}


def rms_bound(D: float, ke: float, V_T1: float, window: float) -> float:
    """RMS tracking-error bound over a window of length `window` starting at V(T₁)."""
    return D / (2 * ke) + np.sqrt(D ** 2 / (4 * ke ** 2) + V_T1 / (ke * window))


def rms_over_window(res: dict, T1: float):
    mask = res['t'] >= T1
    return float(np.sqrt(np.mean(res['e_norm'][mask] ** 2))), mask
