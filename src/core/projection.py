"""
Lipschitz ("smooth") Convex Projection Operator
================================================

Implements the projection operator Proj_Θ(θ, Γy) — called "smooth" in the Pomet–Praly
literature, but only LOCALLY LIPSCHITZ, not C¹ (see the regularity note below) — that keeps adaptive
parameters inside convex projection sets while preserving the inequality the
Lyapunov argument needs (manuscript Section 5):

    θ̃ᵀ Γ⁻¹ Proj_Θ(θ, Γy)  ≤  θ̃ᵀ y          with θ̃ = θ − θ*,           (P)

which holds for every θ in the set and every y, PROVIDED the ideal parameter θ*
lies in the *inner* set {θ : ‖θ‖ ≤ θ_max − ε} (strictly inside the boundary
layer).  (P) is the only property of the operator that the proof uses.

Operator (norm ball Ω = {θ : ‖θ‖ ≤ θ_max}, boundary layer of width ε,
Pomet–Praly / Ioannou–Sun / Lavretsky–Wise form):

    Proj(θ, Γy) = Γy − h(θ) · Γθ · (θᵀΓy) / (θᵀΓθ)   if h(θ) > 0 and θᵀΓy > 0
                = Γy                                  otherwise

    h(θ) = 0                       for ‖θ‖ ≤ θ_max − ε
         = ((‖θ‖ − θ_max + ε)/ε)²  for θ_max − ε < ‖θ‖ < θ_max       (C¹ ramp)
         = 1                       for ‖θ‖ ≥ θ_max

Regularity: h is C¹, but the operator equals Γy − h·Γθ·max(θᵀΓy, 0)/(θᵀΓθ), and max(·, 0) has a kink,
so Proj is locally Lipschitz in (θ, y) and not C¹.  Local Lipschitz continuity is all that existence
and uniqueness of classical solutions of the closed loop need.

Proof of (P):  θ̃ᵀΓ⁻¹(Proj − Γy) = −h (θᵀΓy)/(θᵀΓθ) · θ̃ᵀθ, and
θ̃ᵀθ = ‖θ‖² − θ*ᵀθ ≥ ‖θ‖(‖θ‖ − ‖θ*‖) ≥ 0 whenever h > 0 (‖θ‖ ≥ θ_max − ε ≥ ‖θ*‖).

For isotropic Γ = γI the operator reduces to  γ·[y − h (θ̂ᵀy) θ̂]  (the fast path).
For general symmetric positive-definite Γ the direction removed is Γθ, not θ;
using the isotropic formula with a non-isotropic Γ silently breaks (P).  This
module therefore takes Γ explicitly whenever it is not a scalar multiple of I.

Box sets (bandwidths): Ω_σ = {σ : σ_min ≤ σᵢ ≤ σ_max}.  Component-wise LINEAR ramps
(piecewise linear in σ, hence Lipschitz, not C¹) give (P) component-wise for DIAGONAL Γ_σ,
provided σ*ᵢ ∈ [σ_min+ε, σ_max−ε].

Center separation ‖cᵢ − cⱼ‖ ≥ δ is NOT a convex constraint.  The earlier
"repulsive correction" that tried to enforce it injected an update component
that is not covered by (P) (it was measured to violate (P) on ~30% of samples in
tests/test_theorem_to_code.py) and did not even keep the separation.  It is kept
only as an explicitly opt-in legacy option (`enforce_separation=True`) that is
OUTSIDE the proof; the default algorithm is exactly the proved one.
"""

import numpy as np


def _smooth_indicator(norm: float, theta_max: float, epsilon: float) -> float:
    """C¹ ramp h(‖θ‖): 0 in the interior, 1 on/outside the boundary."""
    if norm <= theta_max - epsilon:
        return 0.0
    if norm >= theta_max:
        return 1.0
    s = (norm - (theta_max - epsilon)) / epsilon
    return s * s


class ProjectionOperator:
    """
    Smooth projection onto a norm ball Ω = {θ : ‖θ‖ ≤ θ_max}.

    Parameters
    ----------
    theta_max : float
        Radius of the norm ball.
    epsilon : float
        Boundary-layer thickness, 0 < ε < θ_max.  The ideal parameter must
        satisfy ‖θ*‖ ≤ θ_max − ε for property (P).
    """

    def __init__(self, theta_max: float, epsilon: float = 0.1):
        if theta_max <= 0:
            raise ValueError(f"theta_max must be positive, got {theta_max}")
        if epsilon <= 0 or epsilon >= theta_max:
            raise ValueError(f"epsilon must be in (0, theta_max), got {epsilon}")
        self.theta_max = float(theta_max)
        self.epsilon = float(epsilon)

    # ------------------------------------------------------------------ core
    def indicator(self, theta_norm: float) -> float:
        return _smooth_indicator(theta_norm, self.theta_max, self.epsilon)

    def _smooth_indicator(self, theta_norm: float) -> float:  # backward compat
        return self.indicator(theta_norm)

    def project_vector(self, theta: np.ndarray, y: np.ndarray,
                       gamma: np.ndarray = None) -> np.ndarray:
        """
        Project the (already Γ-scaled) raw update `y = Γ·grad`.

        Parameters
        ----------
        theta : (d,) current parameter.
        y     : (d,) raw update Γ·grad.
        gamma : (d, d) symmetric positive-definite Γ, or None for Γ ∝ I.

        Returns
        -------
        (d,) projected update Proj(θ, Γ·grad).
        """
        theta = np.asarray(theta, dtype=float)
        y = np.asarray(y, dtype=float)
        theta_norm = np.linalg.norm(theta)
        if theta_norm < 1e-12:
            return y.copy()
        h = self.indicator(theta_norm)
        if h <= 0.0:
            return y.copy()
        outward = float(theta @ y)          # ∇gᵀ(Γ grad) up to the factor 2
        if outward <= 0.0:
            return y.copy()
        if gamma is None:
            return y - h * outward / (theta_norm ** 2) * theta
        g_theta = gamma @ theta
        return y - h * outward / float(theta @ g_theta) * g_theta

    def project_matrix(self, W: np.ndarray, Y: np.ndarray,
                       gamma: np.ndarray = None) -> np.ndarray:
        """
        Matrix version with the Frobenius inner product.  `gamma`, if given,
        is (m·n, m·n) and acts on vec(W) (row-major flatten, as everywhere in
        this code base).
        """
        W = np.asarray(W, dtype=float)
        Y = np.asarray(Y, dtype=float)
        if gamma is None:
            W_norm = np.linalg.norm(W, 'fro')
            if W_norm < 1e-12:
                return Y.copy()
            h = self.indicator(W_norm)
            if h <= 0.0:
                return Y.copy()
            outward = float(np.sum(W * Y))
            if outward <= 0.0:
                return Y.copy()
            return Y - h * outward / (W_norm ** 2) * W
        out = self.project_vector(W.flatten(), Y.flatten(), gamma)
        return out.reshape(W.shape)

    # --------------------------------------------------------- hard clipping
    def hard_clip_vector(self, theta: np.ndarray) -> np.ndarray:
        """Radial clip to the ball (discrete-time safety net after a step)."""
        theta_norm = np.linalg.norm(theta)
        if theta_norm > self.theta_max:
            return theta * (self.theta_max / theta_norm)
        return theta.copy()

    def hard_clip_matrix(self, W: np.ndarray) -> np.ndarray:
        W_norm = np.linalg.norm(W, 'fro')
        if W_norm > self.theta_max:
            return W * (self.theta_max / W_norm)
        return W.copy()


class BandwidthProjection:
    """
    Lipschitz box projection for bandwidths, Ω_σ = {σ : σ_min ≤ σᵢ ≤ σ_max}.

    Component-wise linear ramps of width ε at both faces (Lipschitz, not C¹).  Property (P) holds
    component-wise for diagonal Γ_σ provided σ_min + ε ≤ σ*ᵢ ≤ σ_max − ε.
    """

    def __init__(self, sigma_min: float, sigma_max: float, epsilon: float = 0.05):
        if sigma_min <= 0:
            raise ValueError(
                f"sigma_min must be strictly > 0, got {sigma_min}. "
                "Section 1.2 / Section 2 require this (bounded Jacobians).")
        if sigma_max <= sigma_min:
            raise ValueError(f"sigma_max ({sigma_max}) must exceed sigma_min ({sigma_min})")
        if epsilon <= 0 or 2 * epsilon >= sigma_max - sigma_min:
            raise ValueError("epsilon must be positive and < (sigma_max - sigma_min)/2")
        self.sigma_min = float(sigma_min)
        self.sigma_max = float(sigma_max)
        self.epsilon = float(epsilon)

    def project_vector(self, sigma: np.ndarray, y: np.ndarray) -> np.ndarray:
        sigma = np.asarray(sigma, dtype=float)
        y = np.asarray(y, dtype=float)
        result = y.copy()
        near_lower = (sigma < self.sigma_min + self.epsilon) & (y < 0)
        if np.any(near_lower):
            s = np.clip((sigma[near_lower] - self.sigma_min) / self.epsilon, 0.0, 1.0)
            result[near_lower] = y[near_lower] * s
        near_upper = (sigma > self.sigma_max - self.epsilon) & (y > 0)
        if np.any(near_upper):
            s = np.clip((self.sigma_max - sigma[near_upper]) / self.epsilon, 0.0, 1.0)
            result[near_upper] = y[near_upper] * s
        return result

    def hard_clip(self, sigma: np.ndarray) -> np.ndarray:
        return np.clip(sigma, self.sigma_min, self.sigma_max)


class CenterProjection:
    """
    Projection for the stacked center vector c ∈ ℝ^{m·n} onto the norm ball
    ‖c‖ ≤ c_max.  This is the operator used in the proof.

    `enforce_separation=True` additionally applies the legacy pairwise
    repulsion that tried to keep ‖cᵢ − cⱼ‖ ≥ δ_min.  That term is NOT covered
    by property (P) (the separation set is non-convex); it is retained only so
    that the pre-revision experiments can be reproduced, and is off by default.
    """

    def __init__(self, c_max: float, delta_min: float, epsilon: float = 0.1,
                 enforce_separation: bool = False):
        self.c_max = float(c_max)
        self.delta_min = float(delta_min)
        self.epsilon = float(epsilon)
        self.enforce_separation = bool(enforce_separation)
        self._norm_proj = ProjectionOperator(c_max, epsilon)

    def indicator(self, c_norm: float) -> float:
        return self._norm_proj.indicator(c_norm)

    def project_vector(self, c_flat: np.ndarray, y_flat: np.ndarray,
                       n_dim: int, gamma: np.ndarray = None) -> np.ndarray:
        projected = self._norm_proj.project_vector(c_flat, y_flat, gamma)
        if not self.enforce_separation:
            return projected
        return self._legacy_separation(c_flat, projected, n_dim)

    def _legacy_separation(self, c_flat, projected, n_dim):
        """Legacy repulsion (outside the proof).  O(m² n)."""
        m = len(c_flat) // n_dim
        if m <= 1:
            return projected
        c = c_flat.reshape(m, n_dim)
        proj = projected.reshape(m, n_dim).copy()
        diff_pairs = c[:, np.newaxis, :] - c[np.newaxis, :, :]
        dist_matrix = np.sqrt(np.sum(diff_pairs ** 2, axis=2))
        threshold = self.delta_min + self.epsilon
        close_i, close_j = np.where(np.triu(dist_matrix < threshold, k=1))
        for idx in range(len(close_i)):
            i, j = close_i[idx], close_j[idx]
            dist = dist_matrix[i, j]
            if dist > 1e-12:
                direction = diff_pairs[i, j] / dist
            else:
                direction = np.zeros(n_dim)
                direction[0] = 1.0
            relative_update = proj[i] - proj[j]
            compression = -np.dot(relative_update, direction)
            if compression > 0:
                s = max(0.0, (dist - self.delta_min) / self.epsilon)
                correction = (1 - s) * compression * 0.5
                proj[i] += correction * direction
                proj[j] -= correction * direction
        return proj.flatten()

    def hard_clip(self, c_flat: np.ndarray, n_dim: int = None) -> np.ndarray:
        return self._norm_proj.hard_clip_vector(c_flat)


def min_center_separation(c: np.ndarray) -> float:
    """Diagnostic: smallest pairwise distance between centers (m, n)."""
    c = np.asarray(c, dtype=float)
    if c.shape[0] < 2:
        return np.inf
    d = np.linalg.norm(c[:, None, :] - c[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(d.min())
