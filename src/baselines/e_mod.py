"""
Classical e-modification Baseline Controller
=============================================

Implements the e-modification adaptive controller.
Mathematical reference: Ioannou & Sun (1996).

Similar to σ-modification, but the leakage term is proportional
to the error norm: -σ_mod * ‖e‖ * Γ_W * W. This reduces the
persistent bias compared to pure σ-modification but does not
eliminate it entirely like the proposed projection operator.
"""

import numpy as np
from ..core.rbf_kernel import RBFKernel


class EModController:
    """
    Classical adaptive controller with e-modification.

    Parameters
    ----------
    m : int
        Number of RBF kernels.
    n : int
        State-space dimension.
    ke : float
        Error feedback gain scalar.
    gamma_W : float or np.ndarray
        Weight learning rate.
    sigma_mod : float
        Leaky modification gain (σ_mod > 0).
    """

    def __init__(self, m: int, n: int, ke: float = 5.0,
                 gamma_W: float = 1.0, sigma_mod: float = 0.1):
        self.m = m
        self.n = n
        self.ke = ke
        self.sigma_mod = sigma_mod

        if np.isscalar(gamma_W):
            self.gamma_W = gamma_W * np.eye(m * n)
        else:
            self.gamma_W = np.asarray(gamma_W)

        self.kernel = RBFKernel(m, n)
        
    def control_law(self, x: np.ndarray, x_ref: np.ndarray,
                    x_ref_dot: np.ndarray,
                    W: np.ndarray, c: np.ndarray,
                    sigma: np.ndarray) -> np.ndarray:
        """
        u(t) = ẋₘ(t) - kₑ·e(t) - f̂(x, t)
        """
        e = x - x_ref
        f_hat = self.kernel.f_hat(x, W, c, sigma)
        return x_ref_dot - self.ke * e - f_hat

    def adaptation_laws(self, x: np.ndarray, e: np.ndarray,
                        W: np.ndarray, c: np.ndarray,
                        sigma: np.ndarray):
        """
        Ẇ = Γ_W · (φ · eᵀ - σ_mod · ‖e‖ · W)
        ċ = 0
        σ̇ = 0
        """
        phi = self.kernel.phi(x, c, sigma)  # (m,)
        e_norm = np.linalg.norm(e)

        # Raw gradient: φ · eᵀ
        phi_e_outer = np.outer(phi, e)  # (m, n)

        # Apply learning rate to both gradient and leakage
        # vectorized: Ẇ = Γ_W @ vec(φ·eᵀ - σ_mod·‖e‖·W)
        grad_flat = phi_e_outer.flatten()
        leak_flat = self.sigma_mod * e_norm * W.flatten()

        W_dot_flat = self.gamma_W @ (grad_flat - leak_flat)
        W_dot = W_dot_flat.reshape(self.m, self.n)

        # Fixed centers and bandwidths
        c_dot = np.zeros_like(c.flatten())
        sigma_dot = np.zeros_like(sigma)

        return W_dot, c_dot, sigma_dot

    # Dummy projection clips for integrator compatibility
    class _DummyProj:
        def hard_clip_matrix(self, W): return W
        def hard_clip(self, v, n=None): return v

    @property
    def proj_W(self): return self._DummyProj()
    @property
    def proj_c(self): return self._DummyProj()
    @property
    def proj_sigma(self): return self._DummyProj()
