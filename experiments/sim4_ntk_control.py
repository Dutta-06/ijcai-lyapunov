"""
Simulation 4: NTK-Based Control (Linearized Neural Network)
=============================================================

A finite-width neural network controller whose dynamics are approximated
via its Neural Tangent Kernel. The network is updated via gradient descent,
which is equivalent to kernel regression with the NTK.
No projection operator -- standard gradient clipping is used instead.

Proves: Whether NTK-based implicit kernel adaptation can match explicit
Lyapunov-guided kernel migration. The key comparison: NTK has no hard
parameter bounds, the proposed architecture does.
"""

import os
import sys
import numpy as np
import torch
import torch.nn as nn
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


# --- NTK Network: wide enough for NTK regime ---
class NTKNetwork(nn.Module):
    """
    Two hidden layers, width=256 (NTK regime approximation).
    Input: [x, t]. Output: f_hat scalar.
    Standard initialization -- NTK is defined at this initialization.
    """
    def __init__(self, width=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(2, width),
            nn.Tanh(),
            nn.Linear(width, width),
            nn.Tanh(),
            nn.Linear(width, 1)
        )
        # NTK-standard initialization (NTK parametrization)
        for layer in self.net:
            if isinstance(layer, nn.Linear):
                nn.init.normal_(layer.weight, 0, 1.0/np.sqrt(layer.in_features))
                nn.init.zeros_(layer.bias)

    def forward(self, x):
        return self.net(x).squeeze(-1)


def run_ntk_control(T=30, dt=0.001, width=256):
    t_eval, xm = generate_lorenz_reference(T, dt)
    N = len(t_eval)

    model = NTKNetwork(width=width)
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-4)  # SGD = NTK regime

    x = np.array([1.0])
    e_history, t_history, comp_time = [], [], []
    f_hat = 0.0
    X_buffer, Y_buffer = [], []

    for i in tqdm(range(N-1), desc="NTK Control"):
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

        # Gradient step on recent data (NTK regime: small LR, many steps)
        if i % 50 == 0 and len(X_buffer) > 20:
            t0 = time.time()
            window = min(100, len(X_buffer))
            Xt = torch.tensor(X_buffer[-window:], dtype=torch.float32)
            Yt = torch.tensor(Y_buffer[-window:], dtype=torch.float32)

            model.train()
            for _ in range(5):
                optimizer.zero_grad()
                preds = model(Xt)
                loss = nn.MSELoss()(preds, Yt)
                loss.backward()
                # Gradient clipping -- NTK has no hard bounds, this is the substitute
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

            model.eval()
            with torch.no_grad():
                x_query = torch.tensor([[x[0], t]], dtype=torch.float32)
                f_hat = model(x_query).item()

            comp_time.append(time.time() - t0)

    return np.array(t_history), np.array(e_history), np.array(comp_time)


if __name__ == "__main__":
    print("Running NTK-based control simulation...")
    os.makedirs('results', exist_ok=True)
    t, e, ct = run_ntk_control(T=30, dt=0.001, width=256)
    np.save("results/ntk_error.npy", e)
    np.save("results/ntk_time.npy", t)
    np.save("results/ntk_comptime.npy", ct)
    print(f"Done. Mean |e|: {e[5000:].mean():.4f}")
    plt.semilogy(t, e, label='NTK Control', color='red')
    plt.xlabel('Time (s)'); plt.ylabel('||e(t)||')
    plt.title('NTK-Based Control -- Tracking Error')
    plt.legend(); plt.savefig('results/sim4_ntk.png', dpi=150)
    print("Plot saved to results/sim4_ntk.png")
