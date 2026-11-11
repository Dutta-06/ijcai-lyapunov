"""
Simulation 3: Deep Gaussian Process (2-Layer)
==============================================

A 2-layer Deep GP where the first GP layer learns a latent representation
of the state and the second GP layer maps that representation to the
dynamics estimate.

Proves: Whether learned hierarchical kernel composition can substitute
for online geometric adaptation. Deep GP is the most sophisticated
static kernel baseline.
"""

import os
import sys
import numpy as np
import torch
import gpytorch
from gpytorch.models.deep_gps import DeepGP, DeepGPLayer
from gpytorch.means import ConstantMean, LinearMean
from gpytorch.kernels import RBFKernel, ScaleKernel
from gpytorch.variational import VariationalStrategy, CholeskyVariationalDistribution
from gpytorch.distributions import MultivariateNormal
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


# --- Deep GP Layer ---
class DGPHiddenLayer(DeepGPLayer):
    def __init__(self, input_dims, output_dims=None, num_inducing=27, mean_type='constant'):
        if output_dims is None:
            inducing_points = torch.randn(num_inducing, input_dims)
            batch_shape = torch.Size([])
        else:
            inducing_points = torch.randn(output_dims, num_inducing, input_dims)
            batch_shape = torch.Size([output_dims])
            
        variational_distribution = CholeskyVariationalDistribution(
            num_inducing_points=num_inducing,
            batch_shape=batch_shape
        )
        variational_strategy = VariationalStrategy(
            self, inducing_points, variational_distribution,
            learn_inducing_locations=True
        )
        super().__init__(variational_strategy, input_dims, output_dims)
        self.mean_module = ConstantMean(batch_shape=batch_shape) \
            if mean_type == 'constant' else LinearMean(input_dims, batch_shape=batch_shape)
        self.covar_module = ScaleKernel(
            RBFKernel(batch_shape=batch_shape, ard_num_dims=input_dims),
            batch_shape=batch_shape
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return MultivariateNormal(mean_x, covar_x)

    def __call__(self, x, *other_inputs, **kwargs):
        if other_inputs:
            processed = [inp.unsqueeze(0).expand(
                self.num_outputs, *inp.shape) for inp in other_inputs]
            x = torch.cat([x] + processed, dim=-1)
        return super().__call__(x, **kwargs)


# --- 2-Layer Deep GP ---
class DeepGPModel(DeepGP):
    def __init__(self, train_x_shape):
        hidden_layer = DGPHiddenLayer(
            input_dims=train_x_shape[-1],
            output_dims=3,
            num_inducing=27,
            mean_type='linear'
        )
        last_layer = DGPHiddenLayer(
            input_dims=hidden_layer.output_dims,
            output_dims=None,
            num_inducing=27,
            mean_type='constant'
        )
        super().__init__()
        self.hidden_layer = hidden_layer
        self.last_layer = last_layer
        self.likelihood = gpytorch.likelihoods.GaussianLikelihood()

    def forward(self, inputs):
        hidden_rep = self.hidden_layer(inputs)
        output = self.last_layer(hidden_rep)
        return output

    def predict(self, test_loader):
        with torch.no_grad():
            means = torch.tensor([0.])
            for x_batch, _ in test_loader:
                preds = self.likelihood(self(x_batch))
                means = torch.cat([means, preds.mean.mean(0)])
        return means[1:]


def run_deep_gp(T=30, dt=0.001):
    t_eval, xm = generate_lorenz_reference(T, dt)
    N = len(t_eval)

    # Pre-collect initialization data
    print("  Collecting initialization data for Deep GP...")
    init_X, init_Y = [], []
    x_init = np.array([1.0])
    for i in tqdm(range(min(1000, N-1)), desc="  Init data"):
        t = t_eval[i]
        f_val = f_true(x_init, t)
        init_X.append([x_init[0], t])
        init_Y.append(f_val)
        u = -5.0 * (x_init[0] - xm[i])
        x_init = x_init + dt * (f_true(x_init, t) + u)

    Xt_init = torch.tensor(init_X, dtype=torch.float32)
    Yt_init = torch.tensor(init_Y, dtype=torch.float32)

    model = DeepGPModel(Xt_init.shape)
    optimizer = torch.optim.Adam([
        {'params': model.parameters()}
    ], lr=0.01)
    mll = gpytorch.mlls.DeepApproximateMLL(
        gpytorch.mlls.VariationalELBO(model.likelihood, model, num_data=len(Xt_init))
    )

    # Pre-train
    print("  Pre-training Deep GP...")
    model.train()
    for epoch in tqdm(range(100), desc="  Pre-train"):
        optimizer.zero_grad()
        with gpytorch.settings.num_likelihood_samples(10):
            output = model(Xt_init)
            loss = -mll(output, Yt_init)
            loss.backward()
        optimizer.step()

    # Control loop
    x = np.array([1.0])
    e_history, t_history, comp_time = [], [], []
    f_hat = 0.0
    X_buffer, Y_buffer = [], []

    print("  Running control loop...")
    for i in tqdm(range(N-1), desc="  Deep GP"):
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

        # Fine-tune every 100 steps on recent data
        if i % 100 == 0 and len(X_buffer) > 50:
            t0 = time.time()
            window = min(200, len(X_buffer))
            Xt = torch.tensor(X_buffer[-window:], dtype=torch.float32)
            Yt = torch.tensor(Y_buffer[-window:], dtype=torch.float32)

            model.train()
            for _ in range(5):
                optimizer.zero_grad()
                with gpytorch.settings.num_likelihood_samples(5):
                    output = model(Xt)
                    loss = -mll(output, Yt)
                    loss.backward()
                optimizer.step()

            model.eval()
            with torch.no_grad(), gpytorch.settings.num_likelihood_samples(10):
                x_query = torch.tensor([[x[0], t]], dtype=torch.float32)
                pred = model.likelihood(model(x_query))
                f_hat = pred.mean.mean().item()

            comp_time.append(time.time() - t0)

    return np.array(t_history), np.array(e_history), np.array(comp_time)


if __name__ == "__main__":
    print("Running Deep GP simulation...")
    os.makedirs('results', exist_ok=True)
    t, e, ct = run_deep_gp(T=30, dt=0.001)
    np.save("results/deep_gp_error.npy", e)
    np.save("results/deep_gp_time.npy", t)
    np.save("results/deep_gp_comptime.npy", ct)
    print(f"Done. Mean |e|: {e[5000:].mean():.4f}")
    plt.semilogy(t, e, label='Deep GP (2-layer)', color='purple')
    plt.xlabel('Time (s)'); plt.ylabel('||e(t)||')
    plt.title('Deep GP -- Tracking Error')
    plt.legend(); plt.savefig('results/sim3_deep_gp.png', dpi=150)
    print("Plot saved to results/sim3_deep_gp.png")
