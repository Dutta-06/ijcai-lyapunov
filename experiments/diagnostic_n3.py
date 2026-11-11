"""Diagnostic script for n=3 Proposed Method"""
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core.controller import AdaptiveController
from utils_hd import f_true_lorenz3d, generate_reference, lorenz3d, rk4_step

def run_diagnostic(T=30, dt=0.001, m=125, c_max=50.0):
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
    
    # Init c randomly since m=125 doesn't fit a nice cube easily
    c = np.random.randn(m, n) * 10.0
    sigma = np.ones(m) * 5.0

    e_hist, t_hist = [], []
    proj_c_active_count = 0

    for i in tqdm(range(N-1), desc="Proposed n=3 Diagnostic"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        e = x - x_ref

        u_vec = ctrl.control_law(x, x_ref, xm_dot, W, c, sigma)

        W_dot, c_dot, sigma_dot = ctrl.adaptation_laws(x, e, W, c, sigma)
        W     = W     + dt * W_dot
        c_pre = c     + dt * c_dot.reshape(m, n)
        sigma = sigma + dt * sigma_dot
        W     = ctrl.proj_W.hard_clip_matrix(W)
        
        # Check if c projection activates
        c_flat = ctrl.proj_c.hard_clip(c_pre.flatten(), n)
        if not np.allclose(c_flat, c_pre.flatten()):
            proj_c_active_count += 1
            
        c     = c_flat.reshape(m, n)
        sigma = ctrl.proj_sigma.hard_clip(sigma)

        x = rk4_step(f_true_lorenz3d, x, u_vec, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)

    e_hist = np.array(e_hist)
    print(f"\n--- Diagnostic Results (m={m}, c_max={c_max}) ---")
    print(f"Mean SS |e|: {e_hist[15000:].mean():.4f}")
    print(f"Max |e|: {e_hist.max():.4f}")
    print(f"Projection active on c: {proj_c_active_count}/{(N-1)} timesteps ({(proj_c_active_count/(N-1))*100:.2f}%)")

if __name__ == "__main__":
    run_diagnostic(m=125, c_max=50.0)
