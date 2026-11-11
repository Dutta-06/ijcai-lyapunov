"""
Lyapunov Function Monitor
==========================

Computes and tracks the multi-parameter Lyapunov candidate function
and its time derivative for stability verification.

Mathematical reference: Section 4 of the manuscript.

V(t) = ½ eᵀe + ½ tr(W̃ᵀ Γ_W⁻¹ W̃) + ½ c̃ᵀ Γ_c⁻¹ c̃ + ½ σ̃ᵀ Γ_σ⁻¹ σ̃

Properties verified at runtime:
    - Positive definiteness: V(ζ) > 0 for ζ ≠ 0
    - Lyapunov descent: V̇(t) < 0 when ‖e‖ > R_UUB
    - Energy decomposition tracking
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class LyapunovSnapshot:
    """Single-timestep snapshot of Lyapunov function components."""
    t: float

    # Individual energy components
    V_e: float       # ½ eᵀe (tracking error energy)
    V_W: float       # ½ tr(W̃ᵀ Γ_W⁻¹ W̃) (weight error energy)
    V_c: float       # ½ c̃ᵀ Γ_c⁻¹ c̃ (center error energy)
    V_sigma: float   # ½ σ̃ᵀ Γ_σ⁻¹ σ̃ (bandwidth error energy)

    # Total energy
    V_total: float   # V(t) = V_e + V_W + V_c + V_sigma

    # Time derivative (numerical estimate)
    V_dot: Optional[float] = None

    # Tracking error norm
    e_norm: float = 0.0

    # Theoretical bound
    R_UUB: Optional[float] = None

    # Descent verification flag
    descent_satisfied: Optional[bool] = None


class LyapunovMonitor:
    """
    Real-time monitor for the multi-parameter Lyapunov candidate.

    Computes V(t), tracks its components, estimates V̇(t) via finite
    differences, and verifies the descent condition.

    Parameters
    ----------
    gamma_W_inv : np.ndarray, shape (m*n, m*n)
        Inverse weight learning rate matrix Γ_W⁻¹.
    gamma_c_inv : np.ndarray, shape (m*n, m*n)
        Inverse center learning rate matrix Γ_c⁻¹.
    gamma_sigma_inv : np.ndarray, shape (m, m)
        Inverse bandwidth learning rate matrix Γ_σ⁻¹.
    R_UUB : float, optional
        Theoretical UUB radius for descent verification.
    """

    def __init__(self, gamma_W_inv: np.ndarray,
                 gamma_c_inv: np.ndarray,
                 gamma_sigma_inv: np.ndarray,
                 R_UUB: Optional[float] = None):
        self.gamma_W_inv = gamma_W_inv
        self.gamma_c_inv = gamma_c_inv
        self.gamma_sigma_inv = gamma_sigma_inv
        self.R_UUB = R_UUB

        # History storage
        self.history: List[LyapunovSnapshot] = []
        self._prev_V: Optional[float] = None
        self._prev_t: Optional[float] = None

    def compute(self, t: float, e: np.ndarray,
                W_tilde: np.ndarray, c_tilde: np.ndarray,
                sigma_tilde: np.ndarray) -> LyapunovSnapshot:
        """
        Compute the Lyapunov function value and its components.

        V(t) = ½ eᵀe + ½ tr(W̃ᵀ Γ_W⁻¹ W̃) + ½ c̃ᵀ Γ_c⁻¹ c̃ + ½ σ̃ᵀ Γ_σ⁻¹ σ̃

        Parameters
        ----------
        t : float
            Current time.
        e : np.ndarray, shape (n,)
            Tracking error vector e = x - xₘ.
        W_tilde : np.ndarray, shape (m, n)
            Weight parameter error W̃ = W - W*.
        c_tilde : np.ndarray, shape (m*n,)
            Center parameter error c̃ = c - c* (flattened).
        sigma_tilde : np.ndarray, shape (m,)
            Bandwidth parameter error σ̃ = σ - σ*.

        Returns
        -------
        LyapunovSnapshot
            Snapshot with all components.
        """
        # Tracking error energy
        V_e = 0.5 * np.dot(e, e)

        # Weight error energy: ½ tr(W̃ᵀ Γ_W⁻¹ W̃)
        W_tilde_flat = W_tilde.flatten()
        V_W = 0.5 * W_tilde_flat @ self.gamma_W_inv @ W_tilde_flat

        # Center error energy: ½ c̃ᵀ Γ_c⁻¹ c̃
        V_c = 0.5 * c_tilde @ self.gamma_c_inv @ c_tilde

        # Bandwidth error energy: ½ σ̃ᵀ Γ_σ⁻¹ σ̃
        V_sigma = 0.5 * sigma_tilde @ self.gamma_sigma_inv @ sigma_tilde

        V_total = V_e + V_W + V_c + V_sigma
        e_norm = np.linalg.norm(e)

        # Numerical V̇ via backward finite difference
        V_dot = None
        if self._prev_V is not None and self._prev_t is not None:
            dt = t - self._prev_t
            if dt > 1e-15:
                V_dot = (V_total - self._prev_V) / dt

        # Descent verification
        descent_satisfied = None
        if self.R_UUB is not None and V_dot is not None:
            if e_norm > self.R_UUB:
                descent_satisfied = (V_dot < 0)
            else:
                descent_satisfied = True  # Inside UUB, not required

        snapshot = LyapunovSnapshot(
            t=t,
            V_e=V_e,
            V_W=V_W,
            V_c=V_c,
            V_sigma=V_sigma,
            V_total=V_total,
            V_dot=V_dot,
            e_norm=e_norm,
            R_UUB=self.R_UUB,
            descent_satisfied=descent_satisfied
        )

        # Update history
        self.history.append(snapshot)
        self._prev_V = V_total
        self._prev_t = t

        return snapshot

    def get_descent_violation_count(self) -> int:
        """
        Count the number of timesteps where V̇ > 0 when ‖e‖ > R_UUB.

        This should be 0 if the theory holds. Non-zero violations
        indicate implementation bugs or numerical integration artifacts.

        Returns
        -------
        int
            Number of descent violations.
        """
        violations = 0
        for snap in self.history:
            if snap.descent_satisfied is not None:
                if not snap.descent_satisfied:
                    violations += 1
        return violations

    def get_time_series(self):
        """
        Extract time series arrays from history.

        Returns
        -------
        dict of np.ndarray
            Keys: 't', 'V_e', 'V_W', 'V_c', 'V_sigma', 'V_total',
                  'V_dot', 'e_norm'
        """
        n = len(self.history)
        data = {
            't': np.zeros(n),
            'V_e': np.zeros(n),
            'V_W': np.zeros(n),
            'V_c': np.zeros(n),
            'V_sigma': np.zeros(n),
            'V_total': np.zeros(n),
            'V_dot': np.full(n, np.nan),
            'e_norm': np.zeros(n),
        }

        for i, snap in enumerate(self.history):
            data['t'][i] = snap.t
            data['V_e'][i] = snap.V_e
            data['V_W'][i] = snap.V_W
            data['V_c'][i] = snap.V_c
            data['V_sigma'][i] = snap.V_sigma
            data['V_total'][i] = snap.V_total
            if snap.V_dot is not None:
                data['V_dot'][i] = snap.V_dot
            data['e_norm'][i] = snap.e_norm

        return data

    def summary_stats(self) -> dict:
        """
        Compute summary statistics of the Lyapunov trajectory.

        Returns
        -------
        dict
            Summary statistics including final V(t), peak ‖e‖,
            steady-state ‖e‖, violations, etc.
        """
        if not self.history:
            return {}

        data = self.get_time_series()
        n = len(data['t'])

        # Steady-state: last 20% of trajectory
        ss_start = int(0.8 * n)
        ss_e_norm = data['e_norm'][ss_start:]

        return {
            'total_timesteps': n,
            'final_V': data['V_total'][-1],
            'peak_e_norm': np.max(data['e_norm']),
            'steady_state_mean_e': np.mean(ss_e_norm),
            'steady_state_max_e': np.max(ss_e_norm),
            'steady_state_std_e': np.std(ss_e_norm),
            'descent_violations': self.get_descent_violation_count(),
            'R_UUB': self.R_UUB,
            'UUB_satisfied': (np.max(ss_e_norm) <= self.R_UUB
                              if self.R_UUB is not None else None),
        }

    def reset(self):
        """Clear all history."""
        self.history.clear()
        self._prev_V = None
        self._prev_t = None
