"""Sparse GP Fixed — Van der Pol (n=2)"""
import os, sys, time
import numpy as np
import torch, gpytorch
from tqdm import tqdm
from utils_hd import f_true_vdp, generate_reference, vanderpol, rk4_step

class SparseGPModel(gpytorch.models.ApproximateGP):
    def __init__(self, inducing_points):
        vd = gpytorch.variational.CholeskyVariationalDistribution(
            inducing_points.size(0))
        vs = gpytorch.variational.VariationalStrategy(
            self, inducing_points, vd, learn_inducing_locations=False)
        super().__init__(vs)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(ard_num_dims=3))
    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x))

def run(T=30, dt=0.001, m=27):
    n = 2
    t_eval, xm = generate_reference(vanderpol, [2.0, 0.0], T, dt)
    N = t_eval.shape[0]

    # Fixed inducing points: uniform grid in [x1, x2, t] space
    g = int(np.ceil(m**(1/3))) + 1
    pts = []
    for xi in np.linspace(-3, 3, g):
        for xj in np.linspace(-3, 3, g):
            for ti in np.linspace(0, T, g):
                pts.append([xi, xj, ti])
                if len(pts) == m: break
            if len(pts) == m: break
        if len(pts) == m: break
    ip = torch.tensor(pts, dtype=torch.float32)

    x = np.array([2.0, 0.0])
    e_hist, t_hist, ct_hist = [], [], []
    X_buf, Y_buf = [], []
    f_hat = np.zeros(n)

    models, liks, opts, mlls = [], [], [], []
    for d in range(n):
        lik = gpytorch.likelihoods.GaussianLikelihood()
        mdl = SparseGPModel(ip.clone())
        opt = torch.optim.Adam(list(mdl.parameters()) +
                               list(lik.parameters()), lr=0.01)
        mll = gpytorch.mlls.VariationalELBO(lik, mdl, num_data=1000)
        models.append(mdl); liks.append(lik)
        opts.append(opt); mlls.append(mll)

    for i in tqdm(range(N-1), desc="Sparse GP n=2"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        u = xm_dot - 5.0*(x - x_ref) - f_hat

        x = rk4_step(f_true_vdp, x, u, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)
        X_buf.append([x[0], x[1], t])
        Y_buf.append(f_true_vdp(x, t).tolist())

        if i % 50 == 0 and len(X_buf) > m:
            t0 = time.perf_counter()
            w = min(500, len(X_buf))
            Xt = torch.tensor(X_buf[-w:], dtype=torch.float32)
            f_hat_new = np.zeros(n)
            for d in range(n):
                Yt = torch.tensor([Y_buf[j][d] for j in
                                   range(len(X_buf)-w, len(X_buf))],
                                  dtype=torch.float32)
                models[d].train(); liks[d].train()
                for _ in range(10):
                    opts[d].zero_grad()
                    loss = -mlls[d](models[d](Xt), Yt)
                    loss.backward(); opts[d].step()
                models[d].eval(); liks[d].eval()
                with torch.no_grad():
                    q = torch.tensor([[x[0], x[1], t]], dtype=torch.float32)
                    f_hat_new[d] = liks[d](models[d](q)).mean.item()
            f_hat = f_hat_new
            ct_hist.append(time.perf_counter() - t0)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n2', exist_ok=True)
    t, e, ct = run()
    np.save('results/n2/sparse_gp_error.npy', e)
    np.save('results/n2/sparse_gp_time.npy', t)
    np.save('results/n2/sparse_gp_comptime.npy', ct)
    print(f"n=2 Sparse GP | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
