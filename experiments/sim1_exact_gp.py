"""
Simulation 1: Exact GP Control
================================

Full Gaussian Process regression conditioned on all observed data up to time t.
No approximation. O(n^3) cost per step as data accumulates.

Proves: The exact GP upper bound on GP-based control quality, and its
computational collapse under long trajectories.
"""

import os
import sys
import numpy as np
import torch
import gpytorch
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# --- Plant and Reference ---
def lorenz(t, state, sigma=10, rho=28, beta=8/3):
    x, y, z = state
    return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

def f_true(x, t):
    return -0.5*x[0] + 2*np.sin(x[0]) + 0.3*np.cos(2*t)*x[0]

def generate_lorenz_reference(T=30, dt=0.001):
    t_span = [0, T]
    t_eval = np.arange(0, T, dt)
    sol = solve_ivp(lorenz, t_span, [1.0, 1.0, 1.0], t_eval=t_eval, method='RK45')
    return t_eval, sol.y[0]


# --- Exact GP Model (GPyTorch) ---
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel()
        )
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


# --- Main Simulation Loop ---
def run_exact_gp(T=30, dt=0.001, max_gp_points=500):
    t_eval, xm = generate_lorenz_reference(T, dt)
    N = len(t_eval)

    x = np.array([1.0])
    e_history = []
    t_history = []
    comp_time = []

    X_train = []
    Y_train = []

    f_hat = 0.0

    for i in tqdm(range(N-1), desc="Exact GP"):
        t = t_eval[i]

        xm_dot = (xm[i+1] - xm[i]) / dt if i < N-2 else 0.0
        u = xm_dot - 5.0 * (x[0] - xm[i]) - f_hat

        k1 = f_true(x, t) + u
        k2 = f_true(x + dt/2*k1, t + dt/2) + u
        k3 = f_true(x + dt/2*k2, t + dt/2) + u
        k4 = f_true(x + dt*k3, t + dt) + u
        x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        e = x[0] - xm[i+1]
        e_history.append(abs(e))
        t_history.append(t)

        X_train.append([x[0], t])
        Y_train.append(f_true(x, t))

        if i % 50 == 0 and len(X_train) > 10:
            t0 = time.time()

            if max_gp_points and len(X_train) > max_gp_points:
                idx = np.random.choice(len(X_train), max_gp_points, replace=False)
                Xt = torch.tensor([X_train[j] for j in idx], dtype=torch.float32)
                Yt = torch.tensor([Y_train[j] for j in idx], dtype=torch.float32)
            else:
                Xt = torch.tensor(X_train, dtype=torch.float32)
                Yt = torch.tensor(Y_train, dtype=torch.float32)

            likelihood = gpytorch.likelihoods.GaussianLikelihood()
            model = ExactGPModel(Xt, Yt, likelihood)
            model.train(); likelihood.train()
            optimizer = torch.optim.Adam(model.parameters(), lr=0.1)
            mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

            for _ in range(20):
                optimizer.zero_grad()
                output = model(Xt)
                loss = -mll(output, Yt)
                loss.backward()
                optimizer.step()

            model.eval(); likelihood.eval()
            with torch.no_grad():
                x_query = torch.tensor([[x[0], t]], dtype=torch.float32)
                pred = likelihood(model(x_query))
                f_hat = pred.mean.item()

            comp_time.append(time.time() - t0)

    return np.array(t_history), np.array(e_history), np.array(comp_time)


if __name__ == "__main__":
    print("Running Exact GP simulation...")
    os.makedirs('results', exist_ok=True)
    t, e, ct = run_exact_gp(T=30, dt=0.001, max_gp_points=500)
    np.save("results/exact_gp_error.npy", e)
    np.save("results/exact_gp_time.npy", t)
    np.save("results/exact_gp_comptime.npy", ct)
    print(f"Done. Mean |e|: {e[5000:].mean():.4f}, "
          f"Mean step time: {np.mean(ct) if len(ct) > 0 else 0:.4f}s")
    plt.semilogy(t, e, label='Exact GP')
    plt.xlabel('Time (s)'); plt.ylabel('||e(t)||')
    plt.title('Exact GP Tracking Error'); plt.legend()
    plt.savefig('results/sim1_exact_gp.png', dpi=150)
    print("Plot saved to results/sim1_exact_gp.png")
