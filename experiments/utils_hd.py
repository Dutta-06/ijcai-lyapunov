# utils_hd.py
import numpy as np
from scipy.integrate import solve_ivp

# ── Plant Models ──────────────────────────────────────────────────────────────

def vanderpol(t, state, mu=1.5):
    """n=2. Van der Pol oscillator with parameter mu."""
    x1, x2 = state
    return [x2, mu*(1 - x1**2)*x2 - x1]

def lorenz3d(t, state, sigma=10, rho=28, beta=8/3):
    """n=3. Full Lorenz system as plant dynamics."""
    x, y, z = state
    return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

def coupled_duffing(t, state, alpha=1.0, beta_d=0.25, delta=0.3,
                    gamma=0.3, omega=1.2, k_couple=0.1):
    """
    n=6. Three coupled Duffing oscillators.
    State: [x1, v1, x2, v2, x3, v3]
    Each: ẍᵢ + delta·ẋᵢ - alphaᵢ·xᵢ + beta_d·xᵢ³ = gamma·cos(omega·t) + coupling
    """
    x1, v1, x2, v2, x3, v3 = state
    alpha_t = alpha + 0.3*np.sin(0.5*t)   # time-varying stiffness
    dx1 = v1
    dv1 = -delta*v1 + alpha_t*x1 - beta_d*x1**3 + gamma*np.cos(omega*t) + k_couple*(x2-x1)
    dx2 = v2
    dv2 = -delta*v2 + alpha_t*x2 - beta_d*x2**3 + gamma*np.cos(omega*t) + k_couple*(x1-x2+x3-x2)
    dx3 = v3
    dv3 = -delta*v3 + alpha_t*x3 - beta_d*x3**3 + gamma*np.cos(omega*t) + k_couple*(x2-x3)
    return [dx1, dv1, dx2, dv2, dx3, dv3]

# ── f_true per benchmark ──────────────────────────────────────────────────────

def f_true_vdp(x, t, mu=1.5):
    """True nonlinear dynamics for Van der Pol (n=2)."""
    x1, x2 = x[0], x[1]
    f1 = x2
    f2 = mu*(1 - x1**2)*x2 - x1
    return np.array([f1, f2])

def f_true_lorenz3d(x, t, sigma=10, rho=28, beta=8/3):
    """True nonlinear dynamics for 3D Lorenz plant."""
    return np.array([
        sigma*(x[1]-x[0]),
        x[0]*(rho-x[2])-x[1],
        x[0]*x[1]-beta*x[2]
    ])

def f_true_duffing6d(x, t, alpha=1.0, beta_d=0.25, delta=0.3,
                     gamma=0.3, omega=1.2, k_couple=0.1):
    """True nonlinear dynamics for coupled Duffing (n=6)."""
    alpha_t = alpha + 0.3*np.sin(0.5*t)
    x1,v1,x2,v2,x3,v3 = x
    f = np.array([
        v1,
        -delta*v1 + alpha_t*x1 - beta_d*x1**3 + gamma*np.cos(omega*t) + k_couple*(x2-x1),
        v2,
        -delta*v2 + alpha_t*x2 - beta_d*x2**3 + gamma*np.cos(omega*t) + k_couple*(x1-x2+x3-x2),
        v3,
        -delta*v3 + alpha_t*x3 - beta_d*x3**3 + gamma*np.cos(omega*t) + k_couple*(x2-x3)
    ])
    return f

# ── Reference Trajectory Generator ───────────────────────────────────────────

def generate_reference(plant_fn, x0, T=30, dt=0.001):
    """
    Generates a smooth reference trajectory by integrating the plant
    at slightly different initial conditions (offset reference).
    Returns t_eval, xm (n x N array).
    """
    t_span = [0, T]
    t_eval = np.arange(0, T, dt)
    x0_ref = np.array(x0) * 0.9   # slight offset to create non-trivial tracking
    sol = solve_ivp(plant_fn, t_span, x0_ref.tolist(),
                    t_eval=t_eval, method='RK45', max_step=dt)
    return t_eval, sol.y   # shape: (n, N)

# ── RK4 Step ──────────────────────────────────────────────────────────────────

def rk4_step(f_dyn, x, u, t, dt):
    """Single RK4 integration step. f_dyn(x,t) returns ndarray."""
    k1 = f_dyn(x, t) + u
    k2 = f_dyn(x + dt/2*k1, t + dt/2) + u
    k3 = f_dyn(x + dt/2*k2, t + dt/2) + u
    k4 = f_dyn(x + dt*k3,   t + dt  ) + u
    return x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

# ── Wilcoxon Test ─────────────────────────────────────────────────────────────

def wilcoxon_test(e_proposed, e_baseline, ss_frac=0.6):
    from scipy.stats import wilcoxon
    ss = int(ss_frac * len(e_proposed))
    a = e_proposed[ss:]
    b = e_baseline[ss:len(e_proposed)]
    stat, p = wilcoxon(a[:len(b)], b[:len(a)], alternative='two-sided')
    if p > 0.05:
        return f"p > 0.05 (n.s.)"
    elif p > 0.001:
        return f"p < 0.05"
    else:
        return f"p < 0.001"
