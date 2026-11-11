"""Deep GP — Coupled Duffing (n=6)"""
import os, sys, time
import numpy as np
import torch, gpytorch
from tqdm import tqdm
from utils_hd import f_true_duffing6d, generate_reference, coupled_duffing, rk4_step

class DeepGPLayer(gpytorch.models.deep_gps.DeepGPLayer):
    def __init__(self, input_dims, output_dims, num_inducing=125, mean_type='constant'):
        if output_dims is None:
            inducing_points = torch.randn(num_inducing, input_dims)
            batch_shape = torch.Size([])
        else:
            inducing_points = torch.randn(output_dims, num_inducing, input_dims)
            batch_shape = torch.Size([output_dims])

        vd = gpytorch.variational.CholeskyVariationalDistribution(num_inducing, batch_shape=batch_shape)
        vs = gpytorch.variational.VariationalStrategy(self, inducing_points, vd, learn_inducing_locations=True)
        super().__init__(vs, input_dims, output_dims)

        if mean_type == 'constant':
            self.mean_module = gpytorch.means.ConstantMean(batch_shape=batch_shape)
        else:
            self.mean_module = gpytorch.means.LinearMean(input_dims)
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(batch_shape=batch_shape, ard_num_dims=input_dims),
            batch_shape=batch_shape, ard_num_dims=None
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


class DeepGPModel(gpytorch.models.deep_gps.DeepGP):
    def __init__(self, train_x_shape):
        hidden_layer = DeepGPLayer(
            input_dims=train_x_shape[-1],
            output_dims=6,
            mean_type='linear',
        )
        last_layer = DeepGPLayer(
            input_dims=hidden_layer.output_dims,
            output_dims=None,
            mean_type='constant',
        )
        super().__init__()
        self.hidden_layer = hidden_layer
        self.last_layer = last_layer
        self.likelihood = gpytorch.likelihoods.GaussianLikelihood()

    def forward(self, inputs):
        hidden_rep1 = self.hidden_layer(inputs)
        output = self.last_layer(hidden_rep1)
        return output

def run(T=30, dt=0.001):
    n = 6
    x0 = [1.0, 0.0, -1.0, 0.5, 0.5, -0.5]
    t_eval, xm = generate_reference(coupled_duffing, x0, T, dt)
    N = t_eval.shape[0]

    x = np.array(x0)
    e_hist, t_hist, ct_hist = [], [], []
    X_buf, Y_buf = [], []
    f_hat = np.zeros(n)

    models, opts, mlls = [], [], []
    for d in range(n):
        mdl = DeepGPModel(train_x_shape=torch.Size([1000, 7]))
        opt = torch.optim.Adam(mdl.parameters(), lr=0.01)
        mll = gpytorch.mlls.DeepApproximateMLL(
            gpytorch.mlls.VariationalELBO(mdl.likelihood, mdl, num_data=1000)
        )
        models.append(mdl); opts.append(opt); mlls.append(mll)

    for i in tqdm(range(N-1), desc="Deep GP n=6"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        u = xm_dot - 5.0*(x - x_ref) - f_hat

        x = rk4_step(f_true_duffing6d, x, u, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)
        X_buf.append(x.tolist() + [t])
        Y_buf.append(f_true_duffing6d(x, t).tolist())

        if i % 100 == 0 and len(X_buf) > 100:
            t0 = time.perf_counter()
            w = min(1000, len(X_buf))
            Xt = torch.tensor(X_buf[-w:], dtype=torch.float32)
            f_hat_new = np.zeros(n)
            for d in range(n):
                Yt = torch.tensor([Y_buf[j][d] for j in range(len(X_buf)-w, len(X_buf))], dtype=torch.float32)
                models[d].train()
                for _ in range(5):
                    opts[d].zero_grad()
                    output = models[d](Xt)
                    loss = -mlls[d](output, Yt)
                    loss.backward(); opts[d].step()
                models[d].eval()
                with torch.no_grad():
                    q = torch.tensor([x.tolist() + [t]], dtype=torch.float32)
                    preds = models[d].likelihood(models[d](q))
                    f_hat_new[d] = preds.mean.mean(0).item()
            f_hat = f_hat_new
            ct_hist.append(time.perf_counter() - t0)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n6', exist_ok=True)
    t, e, ct = run()
    np.save('results/n6/deep_gp_error.npy', e)
    np.save('results/n6/deep_gp_time.npy', t)
    np.save('results/n6/deep_gp_comptime.npy', ct)
    print(f"n=6 Deep GP | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
