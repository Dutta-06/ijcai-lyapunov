"""Proposed architecture — Coupled Duffing (n=6)"""
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core.controller import AdaptiveController
from utils_hd import f_true_duffing6d, generate_reference, coupled_duffing, rk4_step

def run(T=30, dt=0.001, m=125):
    n = 6
    x0 = [1.0, 0.0, -1.0, 0.5, 0.5, -0.5]
    t_eval, xm = generate_reference(coupled_duffing, x0, T, dt)
    N = t_eval.shape[0]

    ctrl = AdaptiveController(
        m=m, n=n, ke=5.0,
        gamma_W=30.0, gamma_c=10.0, gamma_sigma=3.0,
        W_max=200.0, c_max=10.0,
        sigma_min=0.3, sigma_max=8.0, delta_min=0.1
    )

    x = np.array(x0)
    W = np.zeros((m, n))
    c = np.random.randn(m, n) * 2.0  # Random init for 6D
    sigma = np.ones(m) * 3.0

    e_hist, t_hist, ct_hist = [], [], []

    for i in tqdm(range(N-1), desc="Proposed n=6"):
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

        x = rk4_step(f_true_duffing6d, x, u_vec, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n6', exist_ok=True)
    t, e, ct = run()
    np.save('results/n6/proposed_error.npy', e)
    np.save('results/n6/proposed_time.npy', t)
    np.save('results/n6/proposed_comptime.npy', ct)
    print(f"n=6 Proposed | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
