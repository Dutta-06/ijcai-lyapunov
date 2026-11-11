"""
Simulation 2: Sparse GP with Fixed Inducing Points
====================================================

GP regression using m=27 inducing points fixed at a uniform grid in state
space before the simulation starts. O(nm^2) cost per update.
Inducing points never move.

Proves: The best that offline-selected static kernel geometry can achieve.
The direct structural comparison to adaptive centers.
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


# --- Sparse GP (SVGP) with fixed inducing points ---
class SparseGPModel(gpytorch.models.ApproximateGP):
    def __init__(self, inducing_points):
        variational_distribution = gpytorch.variational.CholeskyVariationalDistribution(
            inducing_points.size(0)
        )
        variational_strategy = gpytorch.variational.VariationalStrategy(
            self, inducing_points, variational_distribution,
            learn_inducing_locations=False  # FIXED -- this is the key constraint
        )
        super().__init__(variational_strategy)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel()
        )
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def run_sparse_gp_fixed(T=30, dt=0.001, m=27):
    t_eval, xm = generate_lorenz_reference(T, dt)
    N = len(t_eval)

    # Fixed inducing points on uniform grid
    inducing_x = torch.linspace(-3, 3, m).unsqueeze(1)
    inducing_t = torch.linspace(0, T, m).unsqueeze(1)
    inducing_points = torch.cat([inducing_x, inducing_t], dim=1).float()

    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = SparseGPModel(inducing_points)

    optimizer = torch.optim.Adam([
        {'params': model.parameters()},
        {'params': likelihood.parameters()}
    ], lr=0.01)
    mll = gpytorch.mlls.VariationalELBO(likelihood, model, num_data=1000)

    x = np.array([1.0])
    e_history = []
    t_history = []
    comp_time = []
    f_hat = 0.0

    X_buffer = []
    Y_buffer = []

    for i in tqdm(range(N-1), desc="Sparse GP (Fixed)"):
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

        X_buffer.append([x[0], t])
        Y_buffer.append(f_true(x, t))

        # Update every 50 steps
        if i % 50 == 0 and len(X_buffer) > m:
            t0 = time.time()

            window = min(500, len(X_buffer))
            Xt = torch.tensor(X_buffer[-window:], dtype=torch.float32)
            Yt = torch.tensor(Y_buffer[-window:], dtype=torch.float32)

            model.train(); likelihood.train()
            for _ in range(10):
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
    print("Running Sparse GP (Fixed Inducing Points) simulation...")
    os.makedirs('results', exist_ok=True)
    t, e, ct = run_sparse_gp_fixed(T=30, dt=0.001, m=27)
    np.save("results/sparse_gp_fixed_error.npy", e)
    np.save("results/sparse_gp_fixed_time.npy", t)
    np.save("results/sparse_gp_fixed_comptime.npy", ct)
    print(f"Done. Mean |e|: {e[5000:].mean():.4f}")
    plt.semilogy(t, e, label='Sparse GP (Fixed)', color='orange')
    plt.xlabel('Time (s)'); plt.ylabel('||e(t)||')
    plt.title('Sparse GP Fixed Inducing Points -- Tracking Error')
    plt.legend(); plt.savefig('results/sim2_sparse_gp_fixed.png', dpi=150)
    print("Plot saved to results/sim2_sparse_gp_fixed.png")
