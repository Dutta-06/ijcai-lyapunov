"""
Proposition G(b) needs slow-manifold flatness (H2c), not only the rate condition (H2b)
=====================================================================================

Counterexample found in adversarial review round 2 (judge of target G), made reproducible here.

Setting: m = 1, n = 1, f(x) = 2 exp(−(x−1)²/0.5), reference x_m = 0.2 sin t (ω_ref = 1), law (L) with
k_e = 20, gains (g_a, δ, δ) (Γ = I), error started at its quasi-steady value.  The projection never activates.

Round 3: the closed loop and the joint averaged flow (a) are started ON the slow manifold of (b′),

    a(0) = a*(ξ₀) + (δ/g_a) M⁻¹ Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹ ∇L*(ξ₀),     ξ₀ = (c₀, s₀) = (1.75, 0.9),

so that no weight boundary layer enters the comparison.  A row with all gains divided by 4 (R and S unchanged,
averaging error ε_avg divided by 4, horizon multiplied by 4) shows the remaining closed-loop/averaged gap shrinking.

For each case we compare the displacement of the period-mean geometry ξ = (c, s) of the CLOSED LOOP (project code:
TransportController + rk4_closed_loop) with three averaged predictions, all integrated from the same start and
evaluated at the centres of the first and last periods:

  (a)   joint averaged flow       ϑ' = −(m/k_e) G ∇L̄(ϑ)                         (Proposition G(a))
  (b)   variable-projection flow  ξ' = −(mδ/k_e) Γ ∇L*(ξ)                        (Proposition G(b))
  (b′)  induced-metric flow       ξ' = −(mδ/k_e) (Γ⁻¹ + (δ/g_a) Da*ᵀDa*)⁻¹ ∇L*(ξ)

and report the hypotheses along the averaged trajectory: the (H2b) ratio, the slope S = (δ/g_a)‖Da*Γ^{1/2}‖² of
(H2c), and ω_ref over the joint fastest rate (H2a).

Outputs: results/variable_projection/table_slope.md
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import transport as T
from src.core.plant_models import PlantModel

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'variable_projection')
AMP, KE, M, PERIOD = 0.2, 20.0, 1, 2 * np.pi
C0, S0 = 1.75, 0.9
CASES = [  # (δ, g_a, horizon in s, label)
    (0.47, 20.0, 232.0, 'δ = 0.47'),
    (0.47 / 4, 5.0, 928.0, 'δ = 0.47, all gains / 4'),
    (0.047, 20.0, 1238.0, 'δ = 0.047'),
    (0.0047, 20.0, 3000.0, 'δ = 0.0047'),
]

ts_q = np.arange(4000) * PERIOD / 4000
xs_q = AMP * np.sin(ts_q)
fx_q = 2 * np.exp(-(xs_q - 1) ** 2 / 0.5)


class Plant(PlantModel):
    def __init__(self):
        super().__init__(n=1, name='slope')

    def f(self, x, t):
        return np.array([2 * np.exp(-(x[0] - 1) ** 2 / 0.5)])

    def x_ref(self, t):
        return np.array([AMP * np.sin(t)])

    def x_ref_dot(self, t):
        return np.array([AMP * np.cos(t)])


def _phis(c, s):
    d = xs_q - c
    ph = np.exp(-d ** 2 / (2 * s ** 2))
    return ph, ph * d / s ** 2, ph * d ** 2 / s ** 3                 # φ, ∂_cφ, ∂_sφ


def lbar_grad(a, c, s):
    ph, pc, ps = _phis(c, s)
    r = fx_q - a * ph / M
    return 0.5 * np.mean(r ** 2), np.array([-np.mean(r * ph) / M, -np.mean(r * a * pc) / M, -np.mean(r * a * ps) / M])


def a_star(c, s):
    ph = _phis(c, s)[0]
    return M * np.mean(fx_q * ph) / np.mean(ph ** 2)


def lstar_grad(c, s):
    L, g = lbar_grad(a_star(c, s), c, s)
    return L, g[1:]


def da_star(c, s, h=1e-6):
    return np.array([(a_star(c + h, s) - a_star(c - h, s)) / (2 * h), (a_star(c, s + h) - a_star(c, s - h)) / (2 * h)])


def hess(fun, z, h=1e-4):
    d = len(z)
    H = np.zeros((d, d))
    for k in range(d):
        e = np.zeros(d); e[k] = h
        H[:, k] = (fun(z + e) - fun(z - e)) / (2 * h)
    return 0.5 * (H + H.T)


def m0(c, s):
    return np.mean(_phis(c, s)[0] ** 2) / M ** 2


def slow_manifold_weight(delta, ga, c, s):
    """a*(ξ) + (δ/g_a) M⁻¹ Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹ ∇L*(ξ)  (m = 1, Γ = I): the corrected (b′) offset."""
    D = da_star(c, s)
    v = np.linalg.solve(np.eye(2) + (delta / ga) * np.outer(D, D), lstar_grad(c, s)[1])
    return a_star(c, s) + (delta / ga) * float(D @ v) / m0(c, s)


def averaged_flows(delta, ga, T_end, a0):
    Gd = np.array([ga, delta, delta])
    joint = lambda t, z: -(M / KE) * Gd * lbar_grad(*z)[1]
    red = lambda t, z: -(M * delta / KE) * lstar_grad(*z)[1]

    def ind(t, z):
        D = da_star(*z)
        return -(M * delta / KE) * np.linalg.solve(np.eye(2) + (delta / ga) * np.outer(D, D), lstar_grad(*z)[1])
    kw = dict(rtol=1e-10, atol=1e-12, method='DOP853', dense_output=True)
    return (solve_ivp(joint, (0, T_end), [a0, C0, S0], **kw), solve_ivp(red, (0, T_end), [C0, S0], **kw),
            solve_ivp(ind, (0, T_end), [C0, S0], **kw))


def closed_loop(args):
    delta, ga, T_end, label = args
    plant = Plant()
    ctrl = T.TransportController(1, 1, ke=KE, g_a=ga, g_c=delta, g_s=delta,
                                 a_max=20.0, c_max=6.0, s_min=0.25, s_max=4.0, eps=0.05)
    a0 = slow_manifold_weight(delta, ga, C0, S0)
    P0 = np.array([[a0, C0, S0]])
    x0 = plant.x_ref(0.0) + (plant.f(plant.x_ref(0.0), 0) - T.f_hat(plant.x_ref(0.0), P0)) / KE   # quasi-steady error
    dt = PERIOD / 4000
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=T_end, dt=dt, record_every=200)            # 20 samples per period
    P = res['P'][:, 0, :]
    per = 20
    n_per = len(P) // per
    first, last = P[:per].mean(0), P[(n_per - 1) * per:n_per * per].mean(0)
    return dict(delta=delta, ga=ga, T=T_end, label=label, a0=a0, first=first, last=last,
                t_first=0.5 * PERIOD, t_last=(n_per - 0.5) * PERIOD, max_abs_a=float(np.max(np.abs(P[:, 0]))))


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    with mp.Pool(len(CASES)) as pool:
        closed = pool.map(closed_loop, CASES)
    cos = lambda u, v: float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
    nr = lambda u, v: float(np.linalg.norm(u) / np.linalg.norm(v))
    fmt = lambda v: f"({v[0]:+.4f}, {v[1]:+.4f})"
    lines = ["# Proposition G(b) needs slow-manifold flatness: closed loop vs averaged predictions", "",
             f"m = 1, n = 1, f(x) = 2exp(−(x−1)²/0.5), x_m = 0.2 sin t, k_e = {KE:g}, gains (g_a, δ, δ).  The closed loop and the "
             f"joint averaged flow start ON the slow manifold of (b′) at ξ₀ = ({C0}, {S0}) (a*(ξ₀) = {a_star(C0, S0):.3f}).  "
             "Displacements of the period-mean geometry between the centres of the first and last periods.", "",
             "| case | horizon (s) | (H2b) ratio | slope S (H2c) | ω_ref / joint fastest rate (H2a) | closed-loop Δξ = (Δc, Δs) | "
             "(a) joint averaged Δξ | (b) VarPro Δξ | (b′) induced-metric Δξ |", "|---|---|---|---|---|---|---|---|---|"]
    lines2 = ["", "| case | ‖closed‖/‖(b)‖ | cos(closed, (b)) | ‖closed‖/‖(b′)‖ | cos(closed, (b′)) | ‖closed‖/‖(a)‖ | "
                  "cos(closed, (a)) | ‖(a)‖/‖(b′)‖ | max |a| along the run (box 20) |", "|---|---|---|---|---|---|---|---|---|"]
    for cl in closed:
        delta, ga, T_end = cl['delta'], cl['ga'], cl['T']
        sj, sr, si = averaged_flows(delta, ga, T_end, cl['a0'])
        t1, t2 = cl['t_first'], cl['t_last']
        d_closed = cl['last'][1:] - cl['first'][1:]
        d_joint = sj.sol(t2)[1:] - sj.sol(t1)[1:]
        d_red = sr.sol(t2) - sr.sol(t1)
        d_ind = si.sol(t2) - si.sol(t1)
        rat, S, h2a = [], [], []
        Gs = np.sqrt(np.array([ga, delta, delta]))
        for tt in np.linspace(0, T_end, 30):
            a, c, s = sj.sol(tt)
            Hs = hess(lambda z: lstar_grad(*z)[1], np.array([c, s]))
            rat.append(ga * m0(c, s) / (delta * np.max(np.abs(np.linalg.eigvalsh(Hs)))))
            D = da_star(c, s)
            S.append(delta / ga * float(D @ D))
            Hj = hess(lambda z: lbar_grad(*z)[1], np.array([a, c, s]))
            h2a.append(1.0 / ((M / KE) * float(np.max(np.linalg.eigvalsh(Gs[:, None] * Hj * Gs[None, :])))))
        lines.append(f"| {cl['label']} | {T_end:.0f} | {min(rat):.0f}–{max(rat):.0f} | {min(S):.3g}–{max(S):.3g} | "
                     f"{min(h2a):.1f}–{max(h2a):.1f} | {fmt(d_closed)} | {fmt(d_joint)} | {fmt(d_red)} | {fmt(d_ind)} |")
        lines2.append(f"| {cl['label']} | {nr(d_closed, d_red):.3f} | {cos(d_closed, d_red):.4f} | {nr(d_closed, d_ind):.3f} | "
                      f"{cos(d_closed, d_ind):.4f} | {nr(d_closed, d_joint):.3f} | {cos(d_closed, d_joint):.4f} | "
                      f"{nr(d_joint, d_ind):.3f} | {cl['max_abs_a']:.2f} |")
    lines += lines2
    lines += ["", "Reading: compare the cosine and size ratio of the closed loop against (b) and against (b′), row by row, "
                  "together with S.  (H2b) is satisfied in every row.  ‖(a)‖/‖(b′)‖ isolates the O(1/R) error of (b′); "
                  "‖closed‖/‖(a)‖ isolates the averaging and error-channel errors (compare the 'gains / 4' row)."]
    with open(os.path.join(OUT, 'table_slope.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
