"""Proposed architecture — Van der Pol (n=2)"""
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core.controller import AdaptiveController
from src.core.projection import ProjectionOperator, BandwidthProjection, CenterProjection
from utils_hd import f_true_vdp, generate_reference, vanderpol, rk4_step

def run(T=30, dt=0.001, m=27):
    n = 2
    t_eval, xm = generate_reference(vanderpol, [2.0, 0.0], T, dt)
    N = t_eval.shape[0]

    ctrl = AdaptiveController(
        m=m, n=n, ke=5.0,
        gamma_W=50.0, gamma_c=20.0, gamma_sigma=5.0,
        W_max=50.0, c_max=5.0,
        sigma_min=0.1, sigma_max=5.0, delta_min=0.05
    )

    x = np.array([2.0, 0.0])
    W = np.zeros((m, n))
    # Uniform grid over 2D state space
    g = int(np.ceil(np.sqrt(m)))
    cx = np.linspace(-3, 3, g)
    cy = np.linspace(-3, 3, g)
    gx, gy = np.meshgrid(cx, cy)
    c = np.column_stack([gx.ravel()[:m], gy.ravel()[:m]])
    sigma = np.ones(m) * 1.5

    e_hist, t_hist, ct_hist = [], [], []

    for i in tqdm(range(N-1), desc="Proposed n=2"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        e = x - x_ref

        u_vec = ctrl.control_law(x, x_ref, xm_dot, W, c, sigma)

        t0 = time.perf_counter()
        W_dot, c_dot, sigma_dot = ctrl.adaptation_laws(x, e, W, c, sigma)
        W     = W     + dt * W_dot
        c     = c     + dt * c_dot.reshape(m, n)
        sigma = sigma + dt * sigma_dot
        W     = ctrl.proj_W.hard_clip_matrix(W)
        c_flat = ctrl.proj_c.hard_clip(c.flatten(), n)
        c     = c_flat.reshape(m, n)
        sigma = ctrl.proj_sigma.hard_clip(sigma)
        ct_hist.append(time.perf_counter() - t0)

        x = rk4_step(f_true_vdp, x, u_vec, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n2', exist_ok=True)
    t, e, ct = run()
    np.save('results/n2/proposed_error.npy', e)
    np.save('results/n2/proposed_time.npy', t)
    np.save('results/n2/proposed_comptime.npy', ct)
    print(f"n=2 Proposed | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
