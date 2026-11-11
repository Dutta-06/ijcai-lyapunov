"""
Self-Adaptive RBF Kernel Feature Map
=====================================

Gaussian RBF feature map with analytical Jacobians w.r.t. centers and
bandwidths (manuscript Section 1.1):

    φᵢ(x, cᵢ, σᵢ) = exp(−‖x − cᵢ‖² / (2σᵢ²)),      f̂(x) = Wᵀ φ(x, c, σ)

    ∂φᵢ/∂cᵢ = φᵢ (x − cᵢ)/σᵢ²            (φ'_c ∈ ℝ^{m×mn}, block-"diagonal")
    ∂φᵢ/∂σᵢ = φᵢ ‖x − cᵢ‖²/σᵢ³           (φ'_σ ∈ ℝ^{m×m}, diagonal)

Supremum norms used by the disturbance bound D_max (Section 2.1), valid for
all x and all σᵢ ≥ σ_min (r = ‖x − cᵢ‖):

    sup_r  r e^{−r²/2σ²}/σ²  = e^{−1/2}/σ       (at r = σ)
    sup_r  r² e^{−r²/2σ²}/σ³ = 2e^{−1}/σ        (at r = √2 σ)

Because rows of φ'_c have disjoint supports, ‖φ'_c‖₂ = maxᵢ ‖∂φᵢ/∂cᵢ‖ ≤ e^{−1/2}/σ_min,
and ‖φ'_σ‖₂ ≤ 2e^{−1}/σ_min.  Also 0 < φᵢ ≤ 1, so ‖φ‖ ≤ √m.

Note (terminology): Gaussian RBFs are *localized* but have global support;
nothing here relies on compact support.
"""

import numpy as np


class RBFKernel:
    """
    Parameters
    ----------
    m : number of kernels.
    n : state dimension.
    sigma_min : lower bandwidth bound (> 0) used for the Jacobian sup-norms.
    delta_min : nominal center-separation value (diagnostic only; see
                projection.py — separation is not enforced by the proved
                algorithm).
    """

    def __init__(self, m: int, n: int, sigma_min: float = 0.1,
                 delta_min: float = 0.1):
        if sigma_min <= 0:
            raise ValueError(
                f"sigma_min must be strictly positive, got {sigma_min} "
                "(bounded Jacobians require σ ≥ σ_min > 0).")
        if delta_min <= 0:
            raise ValueError(f"delta_min must be strictly positive, got {delta_min}.")
        self.m = int(m)
        self.n = int(n)
        self.sigma_min = float(sigma_min)
        self.delta_min = float(delta_min)

    # ----------------------------------------------------------- feature map
    def phi(self, x: np.ndarray, c: np.ndarray, sigma: np.ndarray) -> np.ndarray:
        diff = x[np.newaxis, :] - c
        sq_dist = np.sum(diff ** 2, axis=1)
        return np.exp(-sq_dist / (2.0 * sigma ** 2))

    def f_hat(self, x: np.ndarray, W: np.ndarray, c: np.ndarray,
              sigma: np.ndarray) -> np.ndarray:
        return W.T @ self.phi(x, c, sigma)

    # -------------------------------------------------------------- Jacobians
    def phi_jacobian_c_blocks(self, x, c, sigma):
        """(m, n) array whose i-th row is ∂φᵢ/∂cᵢ (the nonzero block of φ'_c)."""
        diff = x[np.newaxis, :] - c
        sq_dist = np.sum(diff ** 2, axis=1)
        phi_vec = np.exp(-sq_dist / (2.0 * sigma ** 2))
        return (phi_vec / sigma ** 2)[:, np.newaxis] * diff

    def phi_jacobian_c(self, x: np.ndarray, c: np.ndarray,
                       sigma: np.ndarray) -> np.ndarray:
        """Dense φ'_c ∈ ℝ^{m×mn} (block structure made explicit)."""
        m, n = self.m, self.n
        blocks = self.phi_jacobian_c_blocks(x, c, sigma)
        jac = np.zeros((m, m * n))
        rows = np.repeat(np.arange(m), n)
        cols = np.arange(m * n)
        jac[rows, cols] = blocks.reshape(-1)
        return jac

    def phi_jacobian_sigma_diag(self, x, c, sigma):
        """Diagonal of φ'_σ: ∂φᵢ/∂σᵢ = φᵢ ‖x − cᵢ‖²/σᵢ³."""
        diff = x[np.newaxis, :] - c
        sq_dist = np.sum(diff ** 2, axis=1)
        phi_vec = np.exp(-sq_dist / (2.0 * sigma ** 2))
        return phi_vec * sq_dist / sigma ** 3

    def phi_jacobian_sigma(self, x: np.ndarray, c: np.ndarray,
                           sigma: np.ndarray) -> np.ndarray:
        return np.diag(self.phi_jacobian_sigma_diag(x, c, sigma))

    # --------------------------------------------------------- sup-norm bounds
    def jacobian_sup_norms(self, sigma_min: float = None):
        """
        Analytic suprema (over all x, c and σ ≥ σ_min) of the spectral norms
        ‖φ'_c‖₂ and ‖φ'_σ‖₂.  Returns (J_c, J_sigma).
        """
        s = self.sigma_min if sigma_min is None else float(sigma_min)
        return np.exp(-0.5) / s, 2.0 * np.exp(-1.0) / s

    # ------------------------------------------------------------ diagnostics
    def check_center_separation(self, c: np.ndarray) -> bool:
        from .projection import min_center_separation
        return min_center_separation(c) >= self.delta_min

    def check_bandwidth_floor(self, sigma: np.ndarray) -> bool:
        return bool(np.all(sigma >= self.sigma_min))

    def gram_matrix(self, c: np.ndarray, sigma: np.ndarray) -> np.ndarray:
        """Φᵢⱼ = exp(−‖cᵢ − cⱼ‖² / (σᵢ² + σⱼ²)) — conditioning diagnostic."""
        sq = np.sum((c[:, None, :] - c[None, :, :]) ** 2, axis=2)
        return np.exp(-sq / (sigma[:, None] ** 2 + sigma[None, :] ** 2))
