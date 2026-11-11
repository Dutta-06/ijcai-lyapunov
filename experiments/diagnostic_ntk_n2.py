"""Diagnostic script for n=2 NTK Control"""
import os, sys, time
import numpy as np
import torch, torch.nn as nn
from tqdm import tqdm
from utils_hd import f_true_vdp, generate_reference, vanderpol, rk4_step

class WideNN(nn.Module):
    def __init__(self, in_dim=3, width=256):
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

def get_params_vector(model):
    return torch.cat([p.view(-1) for p in model.parameters()])

def run_ntk_diagnostic(T=30, dt=0.001, width=256):
    n = 2
    t_eval, xm = generate_reference(vanderpol, [2.0, 0.0], T, dt)
    N = t_eval.shape[0]

    x = np.array([2.0, 0.0])
    e_hist = []

    models = [WideNN(in_dim=3, width=width) for _ in range(n)]
    
    # Store initial parameters to measure movement
    initial_params = [get_params_vector(m).clone() for m in models]
    
    opts = [torch.optim.SGD(m.parameters(), lr=0.1) for m in models]
    criterion = nn.MSELoss()
    f_hat = np.zeros(n)

    for i in tqdm(range(N-1), desc="NTK n=2 Diagnostic"):
        t = t_eval[i]
        x_ref = xm[:, i]
        xm_dot = (xm[:, i+1] - xm[:, i]) / dt if i < N-2 else np.zeros(n)
        u = xm_dot - 5.0*(x - x_ref) - f_hat

        x = rk4_step(f_true_vdp, x, u, t, dt)
        e_hist.append(np.linalg.norm(x - xm[:, i+1]))

        f_hat_new = np.zeros(n)
        for d in range(n):
            in_t = torch.tensor([x[0], x[1], t], dtype=torch.float32)
            y_t = torch.tensor([f_true_vdp(x, t)[d]], dtype=torch.float32)
            
            opts[d].zero_grad()
            out = models[d](in_t)
            loss = criterion(out, y_t)
            loss.backward()
            nn.utils.clip_grad_norm_(models[d].parameters(), 1.0)
            opts[d].step()

            with torch.no_grad():
                f_hat_new[d] = models[d](in_t).item()
        
        f_hat = f_hat_new

    e_hist = np.array(e_hist)
    
    print(f"\n--- NTK Diagnostic Results ---")
    print(f"Mean SS |e|: {e_hist[15000:].mean():.4f}")
    
    for d in range(n):
        final_params = get_params_vector(models[d])
        param_diff_norm = torch.norm(final_params - initial_params[d]).item()
        initial_norm = torch.norm(initial_params[d]).item()
        relative_movement = (param_diff_norm / initial_norm) * 100
        print(f"Model {d}: Parameter movement ||W_T - W_0|| = {param_diff_norm:.4f} (Relative: {relative_movement:.2f}%)")
        
if __name__ == "__main__":
    run_ntk_diagnostic()
