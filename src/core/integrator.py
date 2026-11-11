"""
Fixed-Step RK4 Integrator with Projection Enforcement
======================================================

Integrates the coupled plant + controller + adaptation system
using a fixed-step 4th-order Runge-Kutta method.

The projection operator is enforced at every integration step
to maintain hard parameter boundedness (Section 5).

Design decision: Fixed-step RK4 over adaptive solvers for:
    1. Predictable projection enforcement at every substep
    2. Deterministic behavior for Monte Carlo reproducibility
    3. Simpler debugging of stability properties
"""

import numpy as np
from typing import Optional, Callable, Dict, List
from .controller import AdaptiveController
from .plant_models import PlantModel
from .lyapunov_monitor import LyapunovMonitor


def _validate_ideal_params(
    W_star: np.ndarray,
    c_star: np.ndarray,
    sigma_star: np.ndarray,
    W_max: float,
    c_max: float,
    sigma_min: float,
    sigma_max: float,
    epsilon: float,
    t: float = 0.0,
) -> None:
    """
    Verify that ideal (desired) parameters satisfy the interior-margin
    assumption required for the projection inequality to be a valid upper bound.

    The projection inequality
        tr(W̃ᵀ Γ_W⁻¹ Proj_W(W, y)) ≤ tr(W̃ᵀ y)
    holds only when W* lies strictly inside the boundary layer, i.e.
        ‖W*‖_F  ≤  W_max − ε
    (and analogously for c* and σ*). If this is violated, D_max and C_Ω
    are no longer valid sup-bounds and the UUB proof silently breaks.

    Parameters
    ----------
    W_star, c_star, sigma_star : np.ndarray
        Ideal parameter values at the query time.
    W_max, c_max, sigma_min, sigma_max, epsilon : float
        Projection-set parameters matching those used in projection.py.
    t : float
        Current time (used only for the error message).

    Raises
    ------
    ValueError
        If any ideal parameter violates its interior-margin constraint.
    """
    W_norm = np.linalg.norm(W_star, 'fro')
    if W_norm > W_max - epsilon:
        raise ValueError(
            f"[t={t:.4f}] Ideal weight ‖W*‖_F = {W_norm:.4f} violates "
            f"interior-margin constraint W_max - ε = {W_max - epsilon:.4f}. "
            "The projection inequality is no longer a valid upper bound. "
            "Increase W_max or decrease epsilon."
        )

    c_norm = np.linalg.norm(c_star)
    if c_norm > c_max - epsilon:
        raise ValueError(
            f"[t={t:.4f}] Ideal center ‖c*‖ = {c_norm:.4f} violates "
            f"interior-margin constraint c_max - ε = {c_max - epsilon:.4f}. "
            "Increase c_max or decrease epsilon."
        )

    sigma_star = np.atleast_1d(sigma_star)
    if not (np.all(sigma_star >= sigma_min + epsilon) and
            np.all(sigma_star <= sigma_max - epsilon)):
        bad = sigma_star[(sigma_star < sigma_min + epsilon) |
                         (sigma_star > sigma_max - epsilon)]
        raise ValueError(
            f"[t={t:.4f}] Ideal bandwidth σ* has {len(bad)} component(s) "
            f"violating interior-margin [{sigma_min + epsilon:.4f}, "
            f"{sigma_max - epsilon:.4f}]. Values: {bad}. "
            "Adjust sigma_min, sigma_max, or epsilon."
        )


class SimulationState:
    """
    Complete state of the coupled plant-controller system at one timestep.

    Attributes
    ----------
    x : np.ndarray, shape (n,)
        Plant state.
    W : np.ndarray, shape (m, n)
        Adaptive weight matrix.
    c : np.ndarray, shape (m, n)
        Adaptive kernel centers.
    sigma : np.ndarray, shape (m,)
        Adaptive bandwidths.
    """

    def __init__(self, x: np.ndarray, W: np.ndarray,
                 c: np.ndarray, sigma: np.ndarray):
        self.x = x.copy()
        self.W = W.copy()
        self.c = c.copy()
        self.sigma = sigma.copy()

    def copy(self):
        return SimulationState(self.x, self.W, self.c, self.sigma)


class SimulationResult:
    """
    Complete results from a simulation run.

    Stores full time series of all relevant quantities.
    """

    def __init__(self):
        self.times: List[float] = []
        self.states: List[np.ndarray] = []       # x(t)
        self.errors: List[np.ndarray] = []       # e(t)
        self.controls: List[np.ndarray] = []     # u(t)
        self.weights: List[np.ndarray] = []      # W(t)
        self.centers: List[np.ndarray] = []      # c(t)
        self.bandwidths: List[np.ndarray] = []   # σ(t)
        self.e_norms: List[float] = []           # ‖e(t)‖
        self.V_values: List[float] = []          # V(t)
        self.V_dot_values: List[float] = []      # V̇(t)

    def to_arrays(self) -> Dict[str, np.ndarray]:
        """Convert all lists to numpy arrays for analysis."""
        return {
            't': np.array(self.times),
            'x': np.array(self.states),
            'e': np.array(self.errors),
            'u': np.array(self.controls),
            'W': np.array(self.weights),
            'c': np.array(self.centers),
            'sigma': np.array(self.bandwidths),
            'e_norm': np.array(self.e_norms),
            'V': np.array(self.V_values),
            'V_dot': np.array(self.V_dot_values),
        }


class AdaptiveRK4Integrator:
    """
    Fixed-step RK4 integrator for the coupled adaptive control system.

    Integrates:
        ẋ  = f(x, t) + u(t)     (plant dynamics)
        Ẇ  = Proj_W(...)         (weight adaptation)
        ċ  = Proj_c(...)         (center adaptation)
        σ̇  = Proj_σ(...)         (bandwidth adaptation)

    Projection enforcement at every RK4 substep.

    Parameters
    ----------
    plant : PlantModel
        The nonlinear plant model.
    controller : AdaptiveController
        The projection-based adaptive controller.
    dt : float
        Fixed integration timestep. Should be small enough (≈1e-4)
        for stability of the coupled system.
    """

    def __init__(self, plant: PlantModel, controller: AdaptiveController,
                 dt: float = 1e-3):
        self.plant = plant
        self.controller = controller
        self.dt = dt

    def _compute_derivatives(self, state: SimulationState,
                             t: float):
        """
        Compute all state derivatives at a single point.

        Returns
        -------
        x_dot : np.ndarray, shape (n,)
        W_dot : np.ndarray, shape (m, n)
        c_dot : np.ndarray, shape (m*n,)
        sigma_dot : np.ndarray, shape (m,)
        u : np.ndarray, shape (n,)
            The control input (stored for logging).
        """
        x = state.x
        W = state.W
        c = state.c
        sigma = state.sigma

        # Reference trajectory
        x_ref = self.plant.x_ref(t)
        x_ref_dot = self.plant.x_ref_dot(t)

        # Tracking error
        e = x - x_ref

        # Control input
        u = self.controller.control_law(x, x_ref, x_ref_dot, W, c, sigma)

        # Plant dynamics: ẋ = f(x, t) + u
        x_dot = self.plant.f(x, t) + u

        # Adaptation laws
        W_dot, c_dot, sigma_dot = self.controller.adaptation_laws(
            x, e, W, c, sigma
        )

        return x_dot, W_dot, c_dot, sigma_dot, u

    def _apply_projection_clips(self, state: SimulationState):
        """
        Apply hard projection clips as safety net after each RK4 step.
        This should rarely activate if dt is small enough.
        """
        ctrl = self.controller

        # Weight clipping
        state.W = ctrl.proj_W.hard_clip_matrix(state.W)

        # Center clipping
        c_flat = ctrl.proj_c.hard_clip(state.c.flatten(), ctrl.n)
        state.c = c_flat.reshape(ctrl.m, ctrl.n)

        # Bandwidth clipping
        state.sigma = ctrl.proj_sigma.hard_clip(state.sigma)

    def step(self, state: SimulationState, t: float) -> SimulationState:
        """
        Perform one RK4 integration step.

        Parameters
        ----------
        state : SimulationState
            Current system state.
        t : float
            Current time.

        Returns
        -------
        SimulationState
            Updated system state at t + dt.
        """
        dt = self.dt
        m, n = self.controller.m, self.controller.n

        # --- k1 ---
        x1, W1, c1, s1, _ = self._compute_derivatives(state, t)

        # --- k2 ---
        state2 = SimulationState(
            x=state.x + 0.5*dt*x1,
            W=state.W + 0.5*dt*W1,
            c=state.c + (0.5*dt*c1).reshape(m, n),
            sigma=state.sigma + 0.5*dt*s1,
        )
        self._apply_projection_clips(state2)
        x2, W2, c2, s2, _ = self._compute_derivatives(state2, t + 0.5*dt)

        # --- k3 ---
        state3 = SimulationState(
            x=state.x + 0.5*dt*x2,
            W=state.W + 0.5*dt*W2,
            c=state.c + (0.5*dt*c2).reshape(m, n),
            sigma=state.sigma + 0.5*dt*s2,
        )
        self._apply_projection_clips(state3)
        x3, W3, c3, s3, _ = self._compute_derivatives(state3, t + 0.5*dt)

        # --- k4 ---
        state4 = SimulationState(
            x=state.x + dt*x3,
            W=state.W + dt*W3,
            c=state.c + (dt*c3).reshape(m, n),
            sigma=state.sigma + dt*s3,
        )
        self._apply_projection_clips(state4)
        x4, W4, c4, s4, _ = self._compute_derivatives(state4, t + dt)

        # --- RK4 combination ---
        new_state = SimulationState(
            x=state.x + (dt/6.0)*(x1 + 2*x2 + 2*x3 + x4),
            W=state.W + (dt/6.0)*(W1 + 2*W2 + 2*W3 + W4),
            c=state.c + ((dt/6.0)*(c1 + 2*c2 + 2*c3 + c4)).reshape(m, n),
            sigma=state.sigma + (dt/6.0)*(s1 + 2*s2 + 2*s3 + s4),
        )

        # Final projection enforcement
        self._apply_projection_clips(new_state)

        return new_state

    def simulate(self, initial_state: SimulationState,
                 t_start: float, t_end: float,
                 lyapunov_monitor: Optional[LyapunovMonitor] = None,
                 W_star_fn: Optional[Callable] = None,
                 c_star_fn: Optional[Callable] = None,
                 sigma_star_fn: Optional[Callable] = None,
                 record_interval: int = 10,
                 progress_callback: Optional[Callable] = None,
                 use_tqdm: bool = False
                 ) -> SimulationResult:
        """
        Run a complete simulation.

        Parameters
        ----------
        initial_state : SimulationState
            Initial conditions for all state variables.
        t_start : float
            Start time.
        t_end : float
            End time.
        lyapunov_monitor : LyapunovMonitor, optional
            If provided, V(t) is computed and verified at each recorded step.
        W_star_fn, c_star_fn, sigma_star_fn : callable, optional
            Functions returning ideal parameter values at time t.
            Required for Lyapunov monitoring (to compute parameter errors).
            Signature: f(t) -> np.ndarray

            IMPORTANT: The values returned by these functions must satisfy
            the interior-margin constraint at every t:
                ‖W*(t)‖_F  ≤  W_max − ε
                ‖c*(t)‖    ≤  c_max − ε
                σ_min + ε  ≤  σ*(t)  ≤  σ_max − ε
            This is checked at t_start and t_end. If these bounds are
            violated, D_max and C_Ω cease to be valid supremum bounds
            and the UUB radius formula is no longer guaranteed. See
            _validate_ideal_params for details.
        record_interval : int
            Record state every N integration steps (to manage memory).
        progress_callback : callable, optional
            Called with (current_step, total_steps) for progress reporting.

        Returns
        -------
        SimulationResult
            Complete simulation results.
        """
        n_steps = int((t_end - t_start) / self.dt)
        result = SimulationResult()
        state = initial_state.copy()
        t = t_start

        # --- Validate ideal parameter interior-margin at both endpoints ---
        # This enforces the geometric assumption underlying the projection
        # inequality. Without this, D_max/C_Ω may not be valid sup-bounds.
        if W_star_fn is not None or c_star_fn is not None or sigma_star_fn is not None:
            ctrl = self.controller
            _W_max    = ctrl.proj_W.theta_max
            _c_max    = ctrl.proj_c.c_max
            _sig_min  = ctrl.proj_sigma.sigma_min
            _sig_max  = ctrl.proj_sigma.sigma_max
            _eps      = ctrl.proj_W.epsilon  # shared epsilon convention

            for _t_check in [t_start, t_end]:
                _Ws = W_star_fn(_t_check)    if W_star_fn    else np.zeros_like(state.W)
                _cs = c_star_fn(_t_check)    if c_star_fn    else np.zeros_like(state.c.flatten())
                _ss = sigma_star_fn(_t_check) if sigma_star_fn else np.ones_like(state.sigma)
                _validate_ideal_params(
                    _Ws, _cs, _ss,
                    _W_max, _c_max, _sig_min, _sig_max, _eps,
                    t=_t_check
                )

        step_iterator = range(n_steps)
        if use_tqdm:
            try:
                from tqdm import tqdm
                step_iterator = tqdm(step_iterator, desc="Simulating")
            except ImportError:
                pass

        for step_idx in step_iterator:
            # Record state
            if step_idx % record_interval == 0:
                x_ref = self.plant.x_ref(t)
                e = state.x - x_ref

                # Get control for logging
                x_ref_dot = self.plant.x_ref_dot(t)
                u = self.controller.control_law(
                    state.x, x_ref, x_ref_dot,
                    state.W, state.c, state.sigma
                )

                result.times.append(t)
                result.states.append(state.x.copy())
                result.errors.append(e.copy())
                result.controls.append(u.copy())
                result.weights.append(state.W.copy())
                result.centers.append(state.c.copy())
                result.bandwidths.append(state.sigma.copy())
                result.e_norms.append(np.linalg.norm(e))

                # Lyapunov monitoring
                if lyapunov_monitor is not None:
                    W_star = W_star_fn(t) if W_star_fn else np.zeros_like(state.W)
                    c_star = c_star_fn(t) if c_star_fn else np.zeros_like(state.c.flatten())
                    sigma_star = sigma_star_fn(t) if sigma_star_fn else np.ones_like(state.sigma)

                    W_tilde = state.W - W_star
                    c_tilde = state.c.flatten() - c_star.flatten()
                    sigma_tilde = state.sigma - sigma_star

                    snap = lyapunov_monitor.compute(
                        t, e, W_tilde, c_tilde, sigma_tilde
                    )
                    result.V_values.append(snap.V_total)
                    result.V_dot_values.append(
                        snap.V_dot if snap.V_dot is not None else np.nan
                    )

            # Integration step
            state = self.step(state, t)
            t = t_start + (step_idx + 1) * self.dt

            # Progress reporting
            if progress_callback and step_idx % (n_steps // 20 + 1) == 0:
                progress_callback(step_idx, n_steps)

        return result
