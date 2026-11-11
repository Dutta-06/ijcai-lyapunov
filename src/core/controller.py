"""
Certainty-Equivalence Adaptive Controller
==========================================

Implements exactly the algorithm analysed in the manuscript (Sections 1, 5.1)
and in THEORY_REVISED.md.

Control law:
    u = ẋₘ − kₑ e − f̂(x),        f̂ = Wᵀ φ(x, c, σ),   e = x − xₘ

Adaptation laws (projected gradient flows):
    Ẇ  = Proj_W(W, Γ_W φ eᵀ)
    ċ  = Proj_c(c, Γ_c φ'_cᵀ W e)
    σ̇  = Proj_σ(σ, Γ_σ φ'_σᵀ W e)

Γ_W and Γ_c act on the row-major vectorisations vec(W) ∈ ℝ^{mn}, vec(c) ∈ ℝ^{mn};
Γ_σ acts on σ ∈ ℝ^m.  Scalars mean Γ = γI (fast path, O(mn) per step).  Matrix
Γ_W, Γ_c may be any symmetric positive-definite matrix (the Γ-weighted
projection is then used, as the proof requires); Γ_σ must be diagonal.

What the proof gives (THEORY_REVISED.md, Theorem 1): with θ* in the inner
projection sets, V̇ ≤ −kₑ‖e‖² + ‖e‖ D_max + C_Ω along every trajectory, hence
(i) a pointwise ultimate bound and (ii) an RMS bound  limsup RMS(‖e‖) ≤ R_UUB.
The theoretical constants D_max, C_Ω are geometry-sized; `theoretical_D_max`
computes the honest value from the projection geometry so that its looseness can
be reported rather than hidden.
"""

import numpy as np
from .rbf_kernel import RBFKernel
from .projection import ProjectionOperator, BandwidthProjection, CenterProjection


def _is_diag(M: np.ndarray, tol: float = 0.0) -> bool:
    return np.allclose(M, np.diag(np.diagonal(M)), atol=tol)


class AdaptiveController:
    """
    Parameters
    ----------
    m, n : kernels, state dimension.
    ke : error feedback gain (> 0).
    gamma_W, gamma_c : scalar or (m·n, m·n) SPD learning-rate matrices.
    gamma_sigma : scalar or (m, m) diagonal positive matrix.
    W_max, c_max : radii of the Frobenius / Euclidean norm balls for W and vec(c).
    sigma_min, sigma_max : bandwidth box.
    delta_min : nominal center separation (diagnostic only unless
                enforce_separation=True, which is outside the proof).
    proj_epsilon : boundary-layer width of every projection.
    enforce_separation : legacy option; see projection.CenterProjection.
    """

    def __init__(self, m: int, n: int, ke: float = 5.0,
                 gamma_W=1.0, gamma_c=0.5, gamma_sigma=0.3,
                 W_max: float = 10.0, c_max: float = 10.0,
                 sigma_min: float = 0.1, sigma_max: float = 5.0,
                 delta_min: float = 0.1, proj_epsilon: float = 0.1,
                 enforce_separation: bool = False):
        if ke <= 0:
            raise ValueError("ke must be positive")
        self.m, self.n, self.ke = int(m), int(n), float(ke)

        self._gamma_W_is_scalar = np.isscalar(gamma_W)
        self._gamma_c_is_scalar = np.isscalar(gamma_c)
        self._gamma_sigma_is_scalar = np.isscalar(gamma_sigma)

        self.gamma_W = (float(gamma_W) * np.eye(m * n) if self._gamma_W_is_scalar
                        else np.array(gamma_W, dtype=float))
        self.gamma_c = (float(gamma_c) * np.eye(m * n) if self._gamma_c_is_scalar
                        else np.array(gamma_c, dtype=float))
        self.gamma_sigma = (float(gamma_sigma) * np.eye(m) if self._gamma_sigma_is_scalar
                            else np.array(gamma_sigma, dtype=float))

        for name, G, d in (("gamma_W", self.gamma_W, m * n),
                           ("gamma_c", self.gamma_c, m * n),
                           ("gamma_sigma", self.gamma_sigma, m)):
            if G.shape != (d, d):
                raise ValueError(f"{name} must be scalar or {(d, d)}, got {G.shape}")
            if not np.allclose(G, G.T):
                raise ValueError(f"{name} must be symmetric")
            if np.min(np.linalg.eigvalsh(G)) < -1e-12:
                raise ValueError(f"{name} must be positive semidefinite "
                                 "(zero = that parameter group is frozen)")
        if not _is_diag(self.gamma_sigma):
            raise ValueError("gamma_sigma must be diagonal (box projection is "
                             "component-wise; see projection.py)")

        self._gamma_W_val = float(gamma_W) if self._gamma_W_is_scalar else None
        self._gamma_c_val = float(gamma_c) if self._gamma_c_is_scalar else None
        self._gamma_sigma_val = float(gamma_sigma) if self._gamma_sigma_is_scalar else None
        self._gamma_sigma_diag = np.diagonal(self.gamma_sigma).copy()

        # (Pseudo-)inverses for the Lyapunov function — oracle diagnostics only.
        # A zero learning rate freezes a parameter group; its error then enters
        # the disturbance d(t) instead of V (see THEORY_REVISED.md, Remark 3).
        self.gamma_W_inv = np.linalg.pinv(self.gamma_W)
        self.gamma_c_inv = np.linalg.pinv(self.gamma_c)
        with np.errstate(divide='ignore'):
            inv_diag = np.where(self._gamma_sigma_diag > 0,
                                1.0 / np.where(self._gamma_sigma_diag > 0,
                                               self._gamma_sigma_diag, 1.0), 0.0)
        self.gamma_sigma_inv = np.diag(inv_diag)

        self.kernel = RBFKernel(m, n, sigma_min=sigma_min, delta_min=delta_min)
        self.proj_W = ProjectionOperator(W_max, epsilon=proj_epsilon)
        self.proj_c = CenterProjection(c_max, delta_min, epsilon=proj_epsilon,
                                       enforce_separation=enforce_separation)
        self.proj_sigma = BandwidthProjection(sigma_min, sigma_max, epsilon=proj_epsilon)

        self.W_max, self.c_max = float(W_max), float(c_max)
        self.sigma_min, self.sigma_max = float(sigma_min), float(sigma_max)
        self.proj_epsilon = float(proj_epsilon)
        self.enforce_separation = bool(enforce_separation)

    # ------------------------------------------------------------ control law
    def control_law(self, x, x_ref, x_ref_dot, W, c, sigma):
        e = x - x_ref
        return x_ref_dot - self.ke * e - self.kernel.f_hat(x, W, c, sigma)

    # -------------------------------------------------------- adaptation laws
    def adaptation_laws(self, x, e, W, c, sigma):
        """
        Returns (Ẇ (m,n), ċ (m·n,), σ̇ (m,)) — the projected update rates.
        Both code paths compute the same quantities; the dense path exists for
        matrix-valued Γ and is checked against the fast path in the tests.
        """
        We = W @ e                                             # (m,)
        diff = x[np.newaxis, :] - c                            # (m, n)
        sq_dist = np.sum(diff ** 2, axis=1)                    # (m,)
        phi = np.exp(-sq_dist / (2.0 * sigma ** 2))            # (m,)
        grad_W = np.outer(phi, e)                              # φ eᵀ          (m, n)
        grad_c = ((phi * We) / sigma ** 2)[:, np.newaxis] * diff   # φ'_cᵀ W e  (m, n)
        grad_sigma = phi * sq_dist / sigma ** 3 * We           # φ'_σᵀ W e     (m,)

        # Weights
        if self._gamma_W_is_scalar:
            W_dot = self.proj_W.project_matrix(W, self._gamma_W_val * grad_W)
        else:
            raw_W = (self.gamma_W @ grad_W.reshape(-1)).reshape(self.m, self.n)
            W_dot = self.proj_W.project_matrix(W, raw_W, self.gamma_W)

        # Centers
        c_flat = c.reshape(-1)
        if self._gamma_c_is_scalar:
            c_dot = self.proj_c.project_vector(c_flat, self._gamma_c_val * grad_c.reshape(-1), self.n)
        else:
            raw_c = self.gamma_c @ grad_c.reshape(-1)
            c_dot = self.proj_c.project_vector(c_flat, raw_c, self.n, self.gamma_c)

        # Bandwidths (diagonal Γ_σ)
        raw_sigma = self._gamma_sigma_diag * grad_sigma
        sigma_dot = self.proj_sigma.project_vector(sigma, raw_sigma)
        return W_dot, c_dot, sigma_dot

    # ------------------------------------------------------- theoretical bounds
    def compute_theoretical_ruub(self, D_max: float, C_omega: float) -> float:
        """R_UUB = D/(2kₑ) + √(D²/(4kₑ²) + C_Ω/kₑ)  (Section 7)."""
        ke = self.ke
        return D_max / (2 * ke) + np.sqrt(D_max ** 2 / (4 * ke ** 2) + C_omega / ke)

    def theoretical_D_max(self, F_max: float) -> float:
        """
        Honest, computable bound on the lumped disturbance (THEORY_REVISED.md,
        Lemma 2).  With d = f − f̂ + W̃ᵀφ + Wᵀφ'_c c̃ + Wᵀφ'_σ σ̃ and all
        parameters inside their sets:

            ‖d‖ ≤ F_max + 3 W_max √m + W_max (J_c · 2 c_max + J_σ · √m (σ_max − σ_min))

        where F_max = sup_{x ∈ X, t} ‖f(x, t)‖ on the operating region and
        (J_c, J_σ) are the analytic Jacobian sup-norms.  This bound is
        geometry-sized by construction; it is reported, not tuned.
        """
        J_c, J_s = self.kernel.jacobian_sup_norms()
        m = self.m
        return (F_max + 3.0 * self.W_max * np.sqrt(m)
                + self.W_max * (J_c * 2.0 * self.c_max
                                + J_s * np.sqrt(m) * (self.sigma_max - self.sigma_min)))

    def estimate_C_omega(self, V_W: float, V_c: float, V_sigma: float) -> float:
        """
        C_Ω = 2W_max ‖Γ_W⁻¹‖₂ V_W + 2c_max ‖Γ_c⁻¹‖₂ V_c + √m(σ_max−σ_min) ‖Γ_σ⁻¹‖₂ V_σ
        (Section 6, with the tighter bandwidth-error bound).  V_W, V_c, V_σ are
        bounds on ‖Ẇ*‖_F, ‖ċ*‖, ‖σ̇*‖ (Assumption A4 of THEORY_REVISED.md).
        """
        return (2 * self.W_max * np.linalg.norm(self.gamma_W_inv, 2) * V_W
                + 2 * self.c_max * np.linalg.norm(self.gamma_c_inv, 2) * V_c
                + np.sqrt(self.m) * (self.sigma_max - self.sigma_min)
                * np.linalg.norm(self.gamma_sigma_inv, 2) * V_sigma)

    def parameter_energy_cap(self) -> float:
        """V_max,Θ = sup over the projection sets of the parameter part of V."""
        return (0.5 * (2 * self.W_max) ** 2 * np.linalg.norm(self.gamma_W_inv, 2)
                + 0.5 * (2 * self.c_max) ** 2 * np.linalg.norm(self.gamma_c_inv, 2)
                + 0.5 * self.m * (self.sigma_max - self.sigma_min) ** 2
                * np.linalg.norm(self.gamma_sigma_inv, 2))

    # --------------------------------------------------------------- checks
    def params_inside_sets(self, W, c, sigma, tol: float = 1e-9) -> bool:
        """Hard-bound invariant: true iff (W, c, σ) ∈ Ω_W × Ω_c × Ω_σ."""
        return (np.linalg.norm(W, 'fro') <= self.W_max + tol
                and np.linalg.norm(np.asarray(c).reshape(-1)) <= self.c_max + tol
                and np.all(sigma >= self.sigma_min - tol)
                and np.all(sigma <= self.sigma_max + tol))

    def ideal_params_admissible(self, W_star, c_star, sigma_star) -> bool:
        """Assumption A3: θ* strictly inside the boundary layers."""
        eps = self.proj_epsilon
        return (np.linalg.norm(W_star, 'fro') <= self.W_max - eps
                and np.linalg.norm(np.asarray(c_star).reshape(-1)) <= self.c_max - eps
                and np.all(np.asarray(sigma_star) >= self.sigma_min + eps)
                and np.all(np.asarray(sigma_star) <= self.sigma_max - eps))
