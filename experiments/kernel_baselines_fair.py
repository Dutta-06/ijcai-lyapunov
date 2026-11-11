"""
Kernel-method comparison on the 1-D benchmark — fair protocol
==============================================================

Same task as sim0–sim4 (plant ẋ = −0.5x + 2 sin x + 0.3 cos(2t) x + u, Lorenz
x-component reference, kₑ = 5, dt = 1 ms, T = 30 s), with three protocol
defects of the original scripts removed:

  1. Stale predictions.  The GP/NTK baselines refreshed f̂ only every 50–100
     steps and held it constant in between while x moved; the proposed method
     re-evaluated f̂ every step.  Here EVERY method evaluates its current model
     at every step (prediction is cheap); only training is periodic.
  2. Inducing-point geometry.  The sparse GP placed its 27 inducing points on a
     diagonal line x ∈ [−3, 3] ↔ t ∈ [0, T] although the reference spans
     x ∈ [−20, 20].  Here they form a 9 × 3 grid covering the operating region.
  3. Statistics.  One trajectory and a Wilcoxon test over ~12 000 autocorrelated
     time samples.  Here K independent seeds (different reference and plant
     initial conditions) give seed-level paired statistics.

Timing: the full online step (prediction + amortised training) is measured for
every method after warm-up, on the same hardware, in the same process type.

Usage:  python kernel_baselines_fair.py <method> <seed>    (writes results/fair_n1/<method>_seed<seed>.npz)
        methods: proposed exact_gp sparse_gp deep_gp ntk
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import time
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'fair_n1')
T, DT, KE, M = 30.0, 1e-3, 5.0, 27


def lorenz(t, s, sigma=10, rho=28, beta=8 / 3):
    x, y, z = s
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]


def f_true(x, t):
    return -0.5 * x[0] + 2 * np.sin(x[0]) + 0.3 * np.cos(2 * t) * x[0]


def reference(seed):
    rng = np.random.default_rng(seed)
    x0 = np.array([1.0, 1.0, 1.0]) + rng.normal(0, 0.5, size=3)
    t_eval = np.arange(0, T, DT)
    sol = solve_ivp(lorenz, [0, T], x0, t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
    xm = sol.y[0]
    xm_dot = np.gradient(xm, DT)
    return t_eval, xm, xm_dot, rng


def rk4(x, u, t):
    k1 = f_true(x, t) + u
    k2 = f_true(x + DT / 2 * k1, t + DT / 2) + u
    k3 = f_true(x + DT / 2 * k2, t + DT / 2) + u
    k4 = f_true(x + DT * k3, t + DT) + u
    return x + (DT / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def loop(seed, predict, learn, warm=500):
    """Generic closed loop. predict(x, t) -> f̂ ; learn(i, x, t, y) trains periodically (returns seconds or 0)."""
    t_eval, xm, xm_dot, rng = reference(seed)
    N = len(t_eval)
    x = np.array([xm[0] + rng.normal(0, 1.0)])
    e_hist = np.zeros(N - 1); step_time = np.zeros(N - 1); train_time = np.zeros(N - 1)
    for i in range(N - 1):
        t = t_eval[i]
        t0 = time.perf_counter()
        f_hat = predict(x, t)
        u = xm_dot[i] - KE * (x[0] - xm[i]) - f_hat
        tp = time.perf_counter() - t0
        x = rk4(x, u, t)
        y = f_true(x, t + DT)                      # oracle label, as in the original scripts
        t1 = time.perf_counter()
        learn(i, x, t + DT, y)
        tl = time.perf_counter() - t1
        e_hist[i] = abs(x[0] - xm[i + 1]); step_time[i] = tp + tl; train_time[i] = tl
    return dict(t=t_eval[:-1], e=e_hist, step_time=step_time[warm:], train_time=train_time[warm:])


# ----------------------------------------------------------------- proposed
def run_proposed(seed):
    from src.core.controller import AdaptiveController
    ctrl = AdaptiveController(m=M, n=1, ke=KE, gamma_W=50.0, gamma_c=20.0, gamma_sigma=5.0,
                              W_max=50.0, c_max=20.0, sigma_min=0.2, sigma_max=10.0, delta_min=0.05)
    S = dict(W=np.zeros((M, 1)), c=np.linspace(-3.0, 3.0, M).reshape(M, 1), s=np.ones(M) * 1.5)
    _, xm, _, _ = reference(seed)

    def predict(x, t):
        return ctrl.kernel.f_hat(x, S['W'], S['c'], S['s'])[0]

    def learn(i, x, t, y):
        # the controller adapts from the tracking error only (no oracle label is used)
        e = x - np.array([xm[min(i + 1, len(xm) - 1)]])
        Wd, cd, sd = ctrl.adaptation_laws(x, e, S['W'], S['c'], S['s'])
        S['W'] = ctrl.proj_W.hard_clip_matrix(S['W'] + DT * Wd)
        S['c'] = ctrl.proj_c.hard_clip(S['c'].reshape(-1) + DT * cd, 1).reshape(M, 1)
        S['s'] = ctrl.proj_sigma.hard_clip(S['s'] + DT * sd)
    return loop(seed, predict, learn)


# ------------------------------------------------------------------- GP tools
def _torch():
    import torch, gpytorch
    torch.set_num_threads(2)
    torch.manual_seed(0)
    return torch, gpytorch


def run_exact_gp(seed, max_pts=500, every=50, iters=20):
    torch, gpytorch = _torch()

    class ExactGP(gpytorch.models.ExactGP):
        def __init__(self, X, Y, lik):
            super().__init__(X, Y, lik)
            self.mean_module = gpytorch.means.ZeroMean()
            self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=2))

        def forward(self, x):
            return gpytorch.distributions.MultivariateNormal(self.mean_module(x), self.covar_module(x))

    S = dict(model=None, lik=None, X=[], Y=[])
    rng = np.random.default_rng(10_000 + seed)

    # double precision + explicit jitter: in float32 the 500-point kernel matrix of
    # near-duplicate consecutive samples is numerically indefinite (NotPSDError).
    jitter = gpytorch.settings.cholesky_jitter(float=1e-4, double=1e-6) if 'double' in \
        gpytorch.settings.cholesky_jitter.__init__.__code__.co_varnames else gpytorch.settings.cholesky_jitter(1e-6)

    def predict(x, t):
        if S['model'] is None:
            return 0.0
        with torch.no_grad(), jitter:
            q = torch.tensor([[x[0], t]], dtype=torch.float64)
            return S['lik'](S['model'](q)).mean.item()

    def learn(i, x, t, y):
        S['X'].append([x[0], t]); S['Y'].append(y)
        if i % every == 0 and len(S['X']) > 10:
            idx = np.arange(len(S['X']))
            if len(idx) > max_pts:
                idx = rng.choice(idx, max_pts, replace=False)
            Xt = torch.tensor(np.array(S['X'])[idx], dtype=torch.float64)
            Yt = torch.tensor(np.array(S['Y'])[idx], dtype=torch.float64)
            lik = gpytorch.likelihoods.GaussianLikelihood().double()
            model = ExactGP(Xt, Yt, lik).double()
            if S['model'] is not None:                         # warm start hyper-parameters
                model.load_state_dict(S['model'].state_dict(), strict=False)
            model.train(); lik.train()
            opt = torch.optim.Adam(model.parameters(), lr=0.1)
            mll = gpytorch.mlls.ExactMarginalLogLikelihood(lik, model)
            with jitter:
                for _ in range(iters):
                    opt.zero_grad(); loss = -mll(model(Xt), Yt); loss.backward(); opt.step()
            model.eval(); lik.eval()
            S['model'], S['lik'] = model, lik
    return loop(seed, predict, learn)


def run_sparse_gp(seed, every=50, iters=10, window=500):
    torch, gpytorch = _torch()

    class SparseGP(gpytorch.models.ApproximateGP):
        def __init__(self, ip):
            vd = gpytorch.variational.CholeskyVariationalDistribution(ip.size(0))
            vs = gpytorch.variational.VariationalStrategy(self, ip, vd, learn_inducing_locations=False)
            super().__init__(vs)
            self.mean_module = gpytorch.means.ZeroMean()
            self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=2))

        def forward(self, x):
            return gpytorch.distributions.MultivariateNormal(self.mean_module(x), self.covar_module(x))

    gx, gt = np.meshgrid(np.linspace(-20, 20, 9), np.linspace(0, T, 3), indexing='ij')   # 27 points covering (x, t)
    ip = torch.tensor(np.stack([gx.ravel(), gt.ravel()], 1), dtype=torch.float32)
    lik = gpytorch.likelihoods.GaussianLikelihood()
    model = SparseGP(ip)
    opt = torch.optim.Adam(list(model.parameters()) + list(lik.parameters()), lr=0.01)
    mll = gpytorch.mlls.VariationalELBO(lik, model, num_data=1000)
    S = dict(ready=False, X=[], Y=[])

    def predict(x, t):
        if not S['ready']:
            return 0.0
        with torch.no_grad():
            return lik(model(torch.tensor([[x[0], t]], dtype=torch.float32))).mean.item()

    def learn(i, x, t, y):
        S['X'].append([x[0], t]); S['Y'].append(y)
        if i % every == 0 and len(S['X']) > M:
            Xt = torch.tensor(S['X'][-window:], dtype=torch.float32)
            Yt = torch.tensor(S['Y'][-window:], dtype=torch.float32)
            model.train(); lik.train()
            for _ in range(iters):
                opt.zero_grad(); loss = -mll(model(Xt), Yt); loss.backward(); opt.step()
            model.eval(); lik.eval(); S['ready'] = True
    return loop(seed, predict, learn)


def run_deep_gp(seed, every=100, iters=5, window=200):
    torch, gpytorch = _torch()
    from gpytorch.models.deep_gps import DeepGP, DeepGPLayer
    from gpytorch.means import ConstantMean, LinearMean
    from gpytorch.kernels import RBFKernel, ScaleKernel
    from gpytorch.variational import VariationalStrategy, CholeskyVariationalDistribution
    from gpytorch.distributions import MultivariateNormal

    class Layer(DeepGPLayer):
        def __init__(self, input_dims, output_dims=None, num_inducing=27, mean_type='constant'):
            if output_dims is None:
                ip = torch.randn(num_inducing, input_dims); bs = torch.Size([])
            else:
                ip = torch.randn(output_dims, num_inducing, input_dims); bs = torch.Size([output_dims])
            vd = CholeskyVariationalDistribution(num_inducing_points=num_inducing, batch_shape=bs)
            vs = VariationalStrategy(self, ip, vd, learn_inducing_locations=True)
            super().__init__(vs, input_dims, output_dims)
            self.mean_module = ConstantMean(batch_shape=bs) if mean_type == 'constant' else LinearMean(input_dims, batch_shape=bs)
            self.covar_module = ScaleKernel(RBFKernel(batch_shape=bs, ard_num_dims=input_dims), batch_shape=bs)

        def forward(self, x):
            return MultivariateNormal(self.mean_module(x), self.covar_module(x))

    class DGP(DeepGP):
        def __init__(self):
            super().__init__()
            self.hidden = Layer(2, 3, 27, 'linear'); self.last = Layer(3, None, 27, 'constant')
            self.likelihood = gpytorch.likelihoods.GaussianLikelihood()

        def forward(self, x):
            return self.last(self.hidden(x))

    # pre-training on a short open-loop-ish burn-in, as in the original script
    t_eval, xm, _, _ = reference(seed)
    xi = np.array([1.0]); X0, Y0 = [], []
    for i in range(min(1000, len(t_eval) - 1)):
        t = t_eval[i]; X0.append([xi[0], t]); Y0.append(f_true(xi, t))
        xi = xi + DT * (f_true(xi, t) - KE * (xi[0] - xm[i]))
    X0 = torch.tensor(X0, dtype=torch.float32); Y0 = torch.tensor(Y0, dtype=torch.float32)
    model = DGP()
    opt = torch.optim.Adam(model.parameters(), lr=0.01)
    mll = gpytorch.mlls.DeepApproximateMLL(gpytorch.mlls.VariationalELBO(model.likelihood, model, num_data=len(X0)))
    model.train()
    for _ in range(100):
        opt.zero_grad()
        with gpytorch.settings.num_likelihood_samples(10):
            loss = -mll(model(X0), Y0); loss.backward()
        opt.step()
    model.eval()
    S = dict(X=[], Y=[])

    def predict(x, t):
        with torch.no_grad(), gpytorch.settings.num_likelihood_samples(10):
            return model.likelihood(model(torch.tensor([[x[0], t]], dtype=torch.float32))).mean.mean().item()

    def learn(i, x, t, y):
        S['X'].append([x[0], t]); S['Y'].append(y)
        if i % every == 0 and len(S['X']) > 50:
            Xt = torch.tensor(S['X'][-window:], dtype=torch.float32); Yt = torch.tensor(S['Y'][-window:], dtype=torch.float32)
            model.train()
            for _ in range(iters):
                opt.zero_grad()
                with gpytorch.settings.num_likelihood_samples(5):
                    loss = -mll(model(Xt), Yt); loss.backward()
                opt.step()
            model.eval()
    return loop(seed, predict, learn)


def run_ntk(seed, every=50, iters=5, window=100, width=256):
    torch, _ = _torch()
    import torch.nn as nn
    net = nn.Sequential(nn.Linear(2, width), nn.Tanh(), nn.Linear(width, width), nn.Tanh(), nn.Linear(width, 1))
    for layer in net:
        if isinstance(layer, nn.Linear):
            nn.init.normal_(layer.weight, 0, 1.0 / np.sqrt(layer.in_features)); nn.init.zeros_(layer.bias)
    opt = torch.optim.SGD(net.parameters(), lr=1e-4)
    S = dict(X=[], Y=[])

    def predict(x, t):
        with torch.no_grad():
            return net(torch.tensor([[x[0], t]], dtype=torch.float32)).item()

    def learn(i, x, t, y):
        S['X'].append([x[0], t]); S['Y'].append(y)
        if i % every == 0 and len(S['X']) > 20:
            Xt = torch.tensor(S['X'][-window:], dtype=torch.float32); Yt = torch.tensor(S['Y'][-window:], dtype=torch.float32)
            for _ in range(iters):
                opt.zero_grad(); loss = nn.MSELoss()(net(Xt).squeeze(-1), Yt); loss.backward()
                torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step()
    return loop(seed, predict, learn)


RUNNERS = dict(proposed=run_proposed, exact_gp=run_exact_gp, sparse_gp=run_sparse_gp, deep_gp=run_deep_gp, ntk=run_ntk)

if __name__ == "__main__":
    method, seed = sys.argv[1], int(sys.argv[2])
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    r = RUNNERS[method](seed)
    np.savez(os.path.join(OUT, f"{method}_seed{seed}.npz"), **r)
    k = int(0.6 * len(r['e']))
    print(f"{method} seed {seed}: mean SS |e| = {r['e'][k:].mean():.4f}, max SS = {r['e'][k:].max():.3f}, "
          f"median step {np.median(r['step_time']) * 1e3:.3f} ms, mean step {r['step_time'].mean() * 1e3:.3f} ms, wall {time.time() - t0:.0f}s")
