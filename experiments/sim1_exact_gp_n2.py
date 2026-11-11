"""Exact GP — Van der Pol (n=2)"""
import os, sys, time
import numpy as np
import torch, gpytorch
import matplotlib.pyplot as plt
from tqdm import tqdm
from utils_hd import f_true_vdp, generate_reference, vanderpol, rk4_step

class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        # Multi-output: one GP per output dimension
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(ard_num_dims=3)  # [x1, x2, t]
        )
    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x))

def run(T=30, dt=0.001, max_pts=500):
    n = 2
    t_eval, xm = generate_reference(vanderpol, [2.0, 0.0], T, dt)
    N = t_eval.shape[0]

    x = np.array([2.0, 0.0])
    e_hist, t_hist, ct_hist = [], [], []
    X_buf, Y_buf = [], []
    f_hat = np.zeros(n)

    for i in tqdm(range(N-1), desc="Exact GP n=2"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        u = xm_dot - 5.0*(x - x_ref) - f_hat

        x = rk4_step(f_true_vdp, x, u, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)

        X_buf.append([x[0], x[1], t])
        Y_buf.append(f_true_vdp(x, t).tolist())

        if i % 50 == 0 and len(X_buf) > 20:
            t0 = time.perf_counter()
            idx = np.random.choice(len(X_buf),
                                   min(max_pts, len(X_buf)), replace=False)
            Xt = torch.tensor([X_buf[j] for j in idx], dtype=torch.float32)
            f_hat_new = np.zeros(n)
            for d in range(n):
                Yt = torch.tensor([Y_buf[j][d] for j in idx], dtype=torch.float32)
                lik = gpytorch.likelihoods.GaussianLikelihood()
                mdl = ExactGPModel(Xt, Yt, lik)
                mdl.train(); lik.train()
                opt = torch.optim.Adam(mdl.parameters(), lr=0.1)
                mll = gpytorch.mlls.ExactMarginalLogLikelihood(lik, mdl)
                for _ in range(20):
                    opt.zero_grad()
                    loss = -mll(mdl(Xt), Yt)
                    loss.backward(); opt.step()
                mdl.eval(); lik.eval()
                with torch.no_grad():
                    q = torch.tensor([[x[0], x[1], t]], dtype=torch.float32)
                    f_hat_new[d] = lik(mdl(q)).mean.item()
            f_hat = f_hat_new
            ct_hist.append(time.perf_counter() - t0)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n2', exist_ok=True)
    t, e, ct = run()
    np.save('results/n2/exact_gp_error.npy', e)
    np.save('results/n2/exact_gp_time.npy', t)
    np.save('results/n2/exact_gp_comptime.npy', ct)
    print(f"n=2 Exact GP | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
