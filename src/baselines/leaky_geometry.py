"""
Leaky-modification baselines on the SAME adaptive architecture
===============================================================

σ-modification and e-modification (Ioannou & Sun 1996) applied to the identical
network (adaptive W, c, σ with the same Jacobian-based gradient directions as
the proposed controller), with leakage toward a nominal parameter θ₀ instead of
projection:

    Ẇ = Γ_W (φ eᵀ       − κ(e) (W − W₀))
    ċ = Γ_c (φ'_cᵀ W e  − κ(e) (c − c₀))
    σ̇ = Γ_σ (φ'_σᵀ W e  − κ(e) (σ − σ₀))

    κ(e) = σ_mod           (σ-modification)
    κ(e) = σ_mod ‖e‖       (e-modification)

Setting Γ_c = Γ_σ = 0 recovers the classical fixed-grid baselines
(src/baselines/sigma_mod.py, e_mod.py) exactly.

The only hard constraint kept is the numerical guard σᵢ ≥ σ_min (a Gaussian
with σ → 0 is undefined); it is applied identically to every method in the
comparison, so the comparison isolates "leakage vs projection" on W and c.

These baselines are what a fair test of the manuscript's central claim
("projection removes the leakage bias") needs: same geometry, same gradient
directions, only the robustness modification differs.
"""

import numpy as np
from ..core.rbf_kernel import RBFKernel


class _FloorOnly:
    def __init__(self, sigma_min):
        self.sigma_min = sigma_min

    def hard_clip(self, v, n=None):
        return np.maximum(v, self.sigma_min)


class _NoClip:
    def hard_clip_matrix(self, W):
        return W

    def hard_clip(self, v, n=None):
        return v


class LeakyGeometryController:
    """
    Parameters
    ----------
    m, n, ke : as for AdaptiveController.
    gamma_W, gamma_c, gamma_sigma : scalar learning rates (0 freezes a group).
    sigma_mod : leakage gain.
    mode : 'sigma' (σ-modification) or 'e' (e-modification).
    W0, c0, sigma0 : nominal parameters the leakage pulls toward
                     (default: zeros for W, the initial geometry for c, σ —
                     i.e. the a-priori estimate, the standard choice).
    sigma_min : numerical floor on bandwidths.
    """

    def __init__(self, m, n, ke=5.0, gamma_W=1.0, gamma_c=0.0, gamma_sigma=0.0,
                 sigma_mod=0.05, mode='sigma', W0=None, c0=None, sigma0=None,
                 sigma_min=0.1):
        if mode not in ('sigma', 'e'):
            raise ValueError("mode must be 'sigma' or 'e'")
        self.m, self.n, self.ke = int(m), int(n), float(ke)
        self.gamma_W, self.gamma_c, self.gamma_sigma = float(gamma_W), float(gamma_c), float(gamma_sigma)
        self.sigma_mod, self.mode = float(sigma_mod), mode
        self.sigma_min = float(sigma_min)
        self.kernel = RBFKernel(m, n, sigma_min=sigma_min)
        self.W0 = np.zeros((m, n)) if W0 is None else np.array(W0, float)
        self.c0 = None if c0 is None else np.array(c0, float).reshape(-1)
        self.sigma0 = None if sigma0 is None else np.array(sigma0, float)
        self.proj_W = _NoClip()
        self.proj_c = _NoClip()
        self.proj_sigma = _FloorOnly(sigma_min)

    def control_law(self, x, x_ref, x_ref_dot, W, c, sigma):
        e = x - x_ref
        return x_ref_dot - self.ke * e - self.kernel.f_hat(x, W, c, sigma)

    def adaptation_laws(self, x, e, W, c, sigma):
        if self.c0 is None:
            self.c0 = np.array(c, float).reshape(-1)       # first call: nominal = initial
        if self.sigma0 is None:
            self.sigma0 = np.array(sigma, float)
        kappa = self.sigma_mod if self.mode == 'sigma' else self.sigma_mod * np.linalg.norm(e)
        We = W @ e
        diff = x[np.newaxis, :] - c
        sq_dist = np.sum(diff ** 2, axis=1)
        phi = np.exp(-sq_dist / (2.0 * sigma ** 2))
        W_dot = self.gamma_W * (np.outer(phi, e) - kappa * (W - self.W0))
        if self.gamma_c > 0:
            grad_c = ((phi * We) / sigma ** 2)[:, np.newaxis] * diff
            c_dot = self.gamma_c * (grad_c.reshape(-1) - kappa * (c.reshape(-1) - self.c0))
        else:
            c_dot = np.zeros(self.m * self.n)
        if self.gamma_sigma > 0:
            grad_s = phi * sq_dist / sigma ** 3 * We
            s_dot = self.gamma_sigma * (grad_s - kappa * (sigma - self.sigma0))
        else:
            s_dot = np.zeros(self.m)
        return W_dot, c_dot, s_dot
