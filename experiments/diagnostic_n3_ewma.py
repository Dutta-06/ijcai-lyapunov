"""
Diagnostic: n=3 Proposed Method with EWMA-filtered error for center adaptation
===============================================================================

Modification: Replace raw e in c-dot with filtered e_bar:
    e_bar_dot = -lambda * (e_bar - e)
    c_dot = Proj_c(c, Gamma_c * [phi'_c]^T * W * e_bar)

This low-pass filters the chaotic gradient signal while preserving the
Lyapunov structure — the projection inequality still holds because the
filtered signal is a causal transformation of e.
"""
import os, sys, time
import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core.controller import AdaptiveController
from utils_hd import f_true_lorenz3d, generate_reference, lorenz3d, rk4_step


def run_diagnostic_ewma(T=30, dt=0.001, m=125, c_max=50.0, lam=20.0):
    """
    Run n=3 diagnostic with EWMA-filtered error for center adaptation.

    Parameters
    ----------
    lam : float
        Filter bandwidth. Higher = faster tracking of e (less smoothing).
        lam=20 gives a time constant of 1/20 = 50ms.
    """
    n = 3
    t_eval, xm = generate_reference(lorenz3d, [1.0, 1.0, 1.0], T, dt)
    N = t_eval.shape[0]

    ctrl = AdaptiveController(
        m=m, n=n, ke=5.0,
        gamma_W=50.0, gamma_c=20.0, gamma_sigma=5.0,
        W_max=100.0, c_max=c_max,
        sigma_min=0.5, sigma_max=15.0, delta_min=0.1
    )

    x = np.array([1.0, 1.0, 1.0])
    W = np.zeros((m, n))
    c = np.random.randn(m, n) * 10.0
    sigma = np.ones(m) * 5.0

    # Filtered error state — initialized to zero (no prior knowledge)
    e_bar = np.zeros(n)

    e_hist, t_hist = [], []
    proj_c_active_count = 0

    for i in tqdm(range(N - 1), desc=f"EWMA n=3 (lam={lam})"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i + 1] - xm[:, i]) / dt if i < N - 2 else np.zeros(n)
        e = x - x_ref

        # ── Update filtered error: e_bar_dot = -lambda*(e_bar - e) ────────
        e_bar = e_bar + dt * (-lam * (e_bar - e))

        # ── Control law uses raw e (unchanged) ───────────────────────────
        u_vec = ctrl.control_law(x, x_ref, xm_dot, W, c, sigma)

        # ── Adaptation laws: W and sigma use raw e, c uses e_bar ─────────
        # We manually replicate the fast-path logic from the controller
        # but substitute e_bar into the center gradient computation.

        We_raw = W @ e          # for W-dot, sigma-dot (raw e)
        We_bar = W @ e_bar      # for c-dot (filtered e)

        diff = x[np.newaxis, :] - c                             # (m, n)
        sq_dist = np.sum(diff ** 2, axis=1)                     # (m,)
        phi = np.exp(-sq_dist / (2.0 * sigma ** 2))             # (m,)

        # W-dot: standard (raw e)
        raw_W = 50.0 * np.outer(phi, e)
        W_dot = ctrl.proj_W.project_matrix(W, raw_W)

        # c-dot: uses e_bar instead of e
        jvp_c = ((phi * We_bar)[:, np.newaxis] * diff
                 / (sigma[:, np.newaxis] ** 2))                  # (m, n)
        raw_c = 20.0 * jvp_c.flatten()
        c_dot = ctrl.proj_c.project_vector(c.flatten(), raw_c, n)

        # sigma-dot: standard (raw e)
        jvp_sigma = phi * sq_dist / (sigma ** 3) * We_raw
        raw_sigma = 5.0 * jvp_sigma
        sigma_dot = ctrl.proj_sigma.project_vector(sigma, raw_sigma)

        # ── Euler integration + projection ────────────────────────────────
        W = W + dt * W_dot
        c_pre = c + dt * c_dot.reshape(m, n)
        sigma = sigma + dt * sigma_dot

        W = ctrl.proj_W.hard_clip_matrix(W)

        c_flat = ctrl.proj_c.hard_clip(c_pre.flatten(), n)
        if not np.allclose(c_flat, c_pre.flatten()):
            proj_c_active_count += 1
        c = c_flat.reshape(m, n)

        sigma = ctrl.proj_sigma.hard_clip(sigma)

        # ── Plant step ────────────────────────────────────────────────────
        x = rk4_step(f_true_lorenz3d, x, u_vec, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i + 1]))
        t_hist.append(t)

    e_hist = np.array(e_hist)
    ss_mean = e_hist[15000:].mean()
    print(f"\n--- EWMA Diagnostic Results (m={m}, c_max={c_max}, lambda={lam}) ---")
    print(f"Mean SS |e|: {ss_mean:.4f}")
    print(f"Max |e|: {e_hist.max():.4f}")
    print(f"Projection active on c: {proj_c_active_count}/{(N - 1)} "
          f"timesteps ({(proj_c_active_count / (N - 1)) * 100:.2f}%)")
    return ss_mean


if __name__ == "__main__":
    # Run with the same m=125, c_max=50 as the unfiltered diagnostic,
    # so results are directly comparable.
    # Baseline (unfiltered) was: Mean SS |e| = 8.5907, proj active 6.14%
    print("=" * 60)
    print("  EWMA Error Filter Diagnostic for n=3 Lorenz")
    print("  Baseline (raw e): Mean SS |e| = 8.5907")
    print("=" * 60)
    run_diagnostic_ewma(m=125, c_max=50.0, lam=20.0)
