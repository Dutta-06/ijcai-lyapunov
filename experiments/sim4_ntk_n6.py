"""NTK-Based Control — Coupled Duffing (n=6)"""
import os, sys, time
import numpy as np
import torch, torch.nn as nn
from tqdm import tqdm
from utils_hd import f_true_duffing6d, generate_reference, coupled_duffing, rk4_step

class WideNN(nn.Module):
    def __init__(self, in_dim=7, width=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, width),
            nn.ReLU(),
            nn.Linear(width, width),
            nn.ReLU(),
            nn.Linear(width, 1)
        )
        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, std=np.sqrt(2.0/m.in_features))
                nn.init.zeros_(m.bias)
    def forward(self, x):
        return self.net(x)

def run(T=30, dt=0.001, width=12):
    n = 6
    x0 = [1.0, 0.0, -1.0, 0.5, 0.5, -0.5]
    t_eval, xm = generate_reference(coupled_duffing, x0, T, dt)
    N = t_eval.shape[0]

    x = np.array(x0)
    e_hist, t_hist, ct_hist = [], [], []

    models = [WideNN(in_dim=7, width=width) for _ in range(n)]
    opts = [torch.optim.SGD(m.parameters(), lr=0.1) for m in models]
    criterion = nn.MSELoss()
    f_hat = np.zeros(n)

    for i in tqdm(range(N-1), desc="NTK n=6"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        u = xm_dot - 5.0*(x - x_ref) - f_hat

        x = rk4_step(f_true_duffing6d, x, u, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))
        t_hist.append(t)

        t0 = time.perf_counter()
        f_hat_new = np.zeros(n)
        for d in range(n):
            in_t = torch.tensor(x.tolist() + [t], dtype=torch.float32)
            y_t = torch.tensor([f_true_duffing6d(x, t)[d]], dtype=torch.float32)
            
            opts[d].zero_grad()
            out = models[d](in_t)
            loss = criterion(out, y_t)
            loss.backward()
            nn.utils.clip_grad_norm_(models[d].parameters(), 1.0)
            opts[d].step()

            with torch.no_grad():
                f_hat_new[d] = models[d](in_t).item()
        
        f_hat = f_hat_new
        ct_hist.append(time.perf_counter() - t0)

    return np.array(t_hist), np.array(e_hist), np.array(ct_hist)

if __name__ == "__main__":
    os.makedirs('results/n6', exist_ok=True)
    t, e, ct = run()
    np.save('results/n6/ntk_error.npy', e)
    np.save('results/n6/ntk_time.npy', t)
    np.save('results/n6/ntk_comptime.npy', ct)
    print(f"n=6 NTK | Mean SS |e|: {e[15000:].mean():.4f} | "
          f"Mean step: {ct.mean()*1000:.3f} ms")
