"""
Proposition G in numbers: does the slow kernel geometry perform variable-projection descent?
=============================================================================================

Claim (PARADIGM.md §6b): in the three-scale regime

    k_e  ≫  ω_ref  ≫  (weight rates, ALL of them)  ≫  (geometry rate),

the period-averaged geometry velocity is  −(m/k_e)·diag(g_c, g_s)·∇L*(c, s), where
L*(c, s) = min_a ½⟨‖f − (1/m)Σ a_i φ_i‖²⟩_ρ is the profile (variable-projection) loss.

The weight rates are the eigenvalues of (g_a/(k_e m))⟨φφᵀ⟩_ρ; their SMALLEST must exceed the geometry rate.
Each configuration reports these measured ratios, so whether (H2) holds is visible next to the result.

Configurations
  'violating'  — m = 4 overlapping kernels: the Gram matrix is ill-conditioned, the slowest weight mode is
                 slower than the geometry, (H2) fails.  Kept as a documented negative.
  'satisfying' — m = 2 well-separated kernels (well-conditioned Gram), weights 10× faster than the
                 reference is slow, geometry ≥ 50× slower than the slowest weight mode.

Checks per configuration:
  G(a)  cosine (G⁻¹ metric) and magnitude ratio between the period-averaged JOINT velocity of (a, c, s) and
        −(m/k_e)G∇L̄ at the period-mean particles, plus the same split into the weight and geometry blocks;
  G(b)  cosine and magnitude ratio between the period-averaged geometry velocity and −(mδ/k_e)Γ∇L*;
  time scales along the run, including the slope S = (δ/g_a)‖Da*‖² of hypothesis (H2c);
  gap of the weights to a*(c, s); monotonicity of L* and of L̄ sampled once per period.
Outputs: results/variable_projection/{table_<config>.md, figure_<config>.png, raw_sim_<config>.npz}
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import transport as T
from src.core.plant_models import PlantModel

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'variable_projection')
N = 1
PERIOD = 2 * np.pi

CONFIGS = {
    'violating': dict(M=4, KE=20.0, GA=2.0, DELTA=0.1, T=800.0, DT=5e-4,
                      c0=[-1.6, -0.4, 0.4, 1.6], s0=0.7, win=5),
    'satisfying': dict(M=2, KE=20.0, GA=13.0, DELTA=0.004, T=3000.0, DT=2e-3,
                       c0=[-0.5, 0.5], s0=0.6, win=20),
}


class ScalarPlant(PlantModel):
    """ẋ = f(x) + u with f(x) = sin(2x) + 0.3x; periodic reference with fundamental ω_ref = 1."""

    def __init__(self):
        super().__init__(n=1, name="scalar")

    def f(self, x, t):
        return np.array([np.sin(2 * x[0]) + 0.3 * x[0]])

    def x_ref(self, t):
        return np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])

    def x_ref_dot(self, t):
        return np.array([1.5 * np.cos(t) + 1.0 * np.cos(2 * t)])


def features(xs, c, s, M):
    return np.exp(-(xs[:, None] - c[None, :]) ** 2 / (2 * s[None, :] ** 2)) / M      # (N, M): f̂ = F a


def profile_loss(cs, plant, xs, M):
    c, s = cs[:M], cs[M:]
    F = features(xs, c, s, M)
    y = np.array([plant.f(np.array([x]), 0.0)[0] for x in xs])
    a, *_ = np.linalg.lstsq(F, y, rcond=None)
    r = y - F @ a
    return 0.5 * np.mean(r ** 2), a


def grad_profile(cs, plant, xs, M, h=1e-5):
    g = np.zeros_like(cs)
    for k in range(len(cs)):
        e = np.zeros_like(cs); e[k] = h
        g[k] = (profile_loss(cs + e, plant, xs, M)[0] - profile_loss(cs - e, plant, xs, M)[0]) / (2 * h)
    return g


def hess_profile_norm(cs, plant, xs, M, h=1e-3):
    d = len(cs); H = np.zeros((d, d))
    for k in range(d):
        e = np.zeros(d); e[k] = h
        H[:, k] = (grad_profile(cs + e, plant, xs, M) - grad_profile(cs - e, plant, xs, M)) / (2 * h)
    return float(np.max(np.abs(np.linalg.eigvalsh(0.5 * (H + H.T)))))


def time_scales(cfg, plant, xs, cs):
    M, KE, GA, DELTA = cfg['M'], cfg['KE'], cfg['GA'], cfg['DELTA']
    c, s = cs[:M], cs[M:]
    Phi = np.exp(-(xs[:, None] - c[None, :]) ** 2 / (2 * s[None, :] ** 2))
    ev = np.linalg.eigvalsh(Phi.T @ Phi / len(xs))
    w_min, w_max = GA / (KE * M) * ev.min(), GA / (KE * M) * ev.max()
    geo = M * DELTA / KE * hess_profile_norm(cs, plant, xs, M)
    return dict(omega_ref=1.0, w_min=w_min, w_max=w_max, geo=geo, k_e=KE)


def run(name):
    cfg = CONFIGS[name]
    M, KE, GA, DELTA, T_end, dt = cfg['M'], cfg['KE'], cfg['GA'], cfg['DELTA'], cfg['T'], cfg['DT']
    os.makedirs(OUT, exist_ok=True)
    plant = ScalarPlant()
    ctrl = T.TransportController(M, N, ke=KE, g_a=GA, g_c=DELTA, g_s=DELTA,
                                 a_max=20.0, c_max=6.0, s_min=0.25, s_max=4.0, eps=0.05)
    P0 = np.zeros((M, 2 * N + 1)); P0[:, 1] = np.array(cfg['c0']); P0[:, 2] = cfg['s0']
    raw_path = os.path.join(OUT, f'raw_sim_{name}.npz')
    if os.path.exists(raw_path) and os.environ.get('VP_REUSE', '1') == '1':
        z = np.load(raw_path); res = dict(t=z['t'], P=z['P'], e_norm=z['e_norm'], x=z['x'])
    else:
        res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P0, T=T_end, dt=dt, record_every=int(round(PERIOD / dt / 20)))
        np.savez(raw_path, **res)
    ts = np.arange(0, PERIOD, PERIOD / 400); xs = np.array([plant.x_ref(t)[0] for t in ts])
    y = np.array([plant.f(np.array([x]), 0.0)[0] for x in xs])
    t, P = res['t'], res['P']
    per = int(round(PERIOD / (t[1] - t[0])))
    win = cfg.get('win', 5)
    G_diag = np.array([GA, DELTA, DELTA])                                   # n = 1: (a, c, s) gains

    def Lbar(a, c, s):
        return 0.5 * float(np.mean((y - features(xs, c, s, M) @ a) ** 2))

    def grad_Lbar(a, c, s, h=1e-6):
        """∇_{(a, c, s)} L̄, shape (M, 3)."""
        g = np.zeros((M, 3)); base = np.stack([a, c, s], 1)
        for i in range(M):
            for k in range(3):
                Pp, Pm = base.copy(), base.copy(); Pp[i, k] += h; Pm[i, k] -= h
                g[i, k] = (Lbar(*Pp.T) - Lbar(*Pm.T)) / (2 * h)
        return g

    def slope_S(cs, h=1e-5):
        """S := (δ/g_a)·‖D_{(c,s)} a*‖² (Γ = I here since g_c = g_s = δ)."""
        D = np.zeros((M, 2 * M))
        for k in range(2 * M):
            e = np.zeros(2 * M); e[k] = h
            D[:, k] = (profile_loss(cs + e, plant, xs, M)[1] - profile_loss(cs - e, plant, xs, M)[1]) / (2 * h)
        return DELTA / GA * float(np.linalg.norm(D, 2) ** 2)

    def pmean(k):                                                           # period-averaged particles at sample k
        return P[k:k + per].mean(0)

    rows = []
    for k in np.arange(win * per, len(t) - (win + 1) * per, per):
        Pk, Pb, Pf = pmean(k), pmean(k - win * per), pmean(k + win * per)
        geo_now = np.concatenate([Pk[:, 1], Pk[:, 2]])
        v_geo = (np.concatenate([Pf[:, 1], Pf[:, 2]]) - np.concatenate([Pb[:, 1], Pb[:, 2]])) / (2 * win * PERIOD)   # centred
        Lstar, a_star = profile_loss(geo_now, plant, xs, M)
        g = grad_profile(geo_now, plant, xs, M)
        v_pred = -(M / KE) * DELTA * g                                       # −(mδ/k_e)Γ∇L*, Γ = I
        cos = float(v_geo @ v_pred / (np.linalg.norm(v_geo) * np.linalg.norm(v_pred) + 1e-18))
        ratio = float(np.linalg.norm(v_geo) / (np.linalg.norm(v_pred) + 1e-18))
        # direct check of G(a): joint period-averaged velocity vs −(m/k_e) G ∇_ϑ L̄ at the current (period-mean) weights
        v_joint = ((Pf - Pb) / (2 * win * PERIOD)).ravel()
        vj_pred = (-(M / KE) * G_diag[None, :] * grad_Lbar(Pk[:, 0], Pk[:, 1], Pk[:, 2])).ravel()
        Ginv = np.tile(1 / G_diag, M)
        cos_j = float(np.sum(Ginv * v_joint * vj_pred) / np.sqrt(np.sum(Ginv * v_joint ** 2) * np.sum(Ginv * vj_pred ** 2) + 1e-30))
        ratio_j = float(np.sqrt(np.sum(Ginv * v_joint ** 2) / (np.sum(Ginv * vj_pred ** 2) + 1e-30)))
        wb = np.arange(M * 3) % 3 == 0                                       # weight block; ~wb = geometry block
        cosb = lambda u, v: float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-30))
        res_ = Ginv * (v_joint - vj_pred) ** 2
        rows.append(dict(t=t[k], Lstar=Lstar, Lbar=Lbar(Pk[:, 0], Pk[:, 1], Pk[:, 2]), gnorm=float(np.linalg.norm(g)),
                         cos=cos, ratio=ratio, cos_j=cos_j, ratio_j=ratio_j,
                         cos_jw=cosb(v_joint[wb], vj_pred[wb]), cos_jg=cosb(v_joint[~wb], vj_pred[~wb]),
                         res_w=float(res_[wb].sum() / (res_.sum() + 1e-30)),
                         pred_w=float(np.sqrt(np.sum((Ginv * vj_pred ** 2)[wb]) / (np.sum((Ginv * vj_pred ** 2)[~wb]) + 1e-30))),
                         a_gap=float(np.linalg.norm(Pk[:, 0] - a_star) / (np.linalg.norm(a_star) + 1e-12)),
                         e_rms=float(np.sqrt(np.mean(res['e_norm'][k:k + per] ** 2)))))
    R = {key: np.array([r[key] for r in rows]) for key in rows[0]}
    # time scales along the run (5 samples), including the slope S and the joint fastest rate
    samp = np.linspace(win * per, len(t) - (win + 1) * per - 1, 5).astype(int)
    scales = []
    for k in samp:
        Pk = pmean(k); cs = np.concatenate([Pk[:, 1], Pk[:, 2]])
        ts_k = time_scales(cfg, plant, xs, cs)
        a_st = profile_loss(cs, plant, xs, M)[1]
        h = 1e-4; base = np.stack([a_st, Pk[:, 1], Pk[:, 2]], 1).ravel(); d = len(base); H = np.zeros((d, d))
        f0 = lambda v: Lbar(*v.reshape(M, 3).T)
        for i in range(d):
            for j in range(i, d):
                ei = np.zeros(d); ej = np.zeros(d); ei[i] = h; ej[j] = h
                H[i, j] = H[j, i] = (f0(base + ei + ej) - f0(base + ei - ej) - f0(base - ei + ej) + f0(base - ei - ej)) / (4 * h * h)
        Gs = np.sqrt(np.tile(G_diag, M))
        joint_fast = M / KE * float(np.max(np.linalg.eigvalsh(Gs[:, None] * H * Gs[None, :])))
        scales.append(dict(ts_k, S=slope_S(cs), joint=joint_fast))
    S_ = {key: np.array([s[key] for s in scales]) for key in scales[0]}
    rng_ = lambda v: f"{v.min():.3g}–{v.max():.3g}"
    q = max(1, len(R['t']) // 4)
    lines = [f"# Variable-projection check — configuration '{name}' (m={M}, k_e={KE:g}, g_a={GA:g}, δ={DELTA:g}, T={T_end:.0f} s)", "",
             "Time scales along the run (min–max over 5 samples; rates in 1/s):", "",
             "| k_e | ω_ref | joint fastest rate (H2a) | fastest weight | slowest weight | geometry (mδ/k_e)‖Γ^{1/2}∇²L*Γ^{1/2}‖ | slope S (H2c) |",
             "|---|---|---|---|---|---|---|",
             f"| {KE:g} | 1 | {rng_(S_['joint'])} | {rng_(S_['w_max'])} | {rng_(S_['w_min'])} | {rng_(S_['geo'])} | {rng_(S_['S'])} |", "",
             f"(H2a) ω_ref / joint fastest: {rng_(1.0 / S_['joint'])}; (H2b) slowest weight / geometry: "
             f"{rng_(S_['w_min'] / np.maximum(S_['geo'], 1e-30))}; (H2c) S ≪ 1: S = {rng_(S_['S'])}.", "",
             "| quantity | first quarter (median) | last quarter (median) |", "|---|---|---|",
             f"| G(a): cosine(joint velocity, −(m/k_e)G∇L̄), G⁻¹ metric | {np.median(R['cos_j'][:q]):.4f} | {np.median(R['cos_j'][-q:]):.4f} |",
             f"| G(a): magnitude ratio | {np.median(R['ratio_j'][:q]):.3f} | {np.median(R['ratio_j'][-q:]):.3f} |",
             f"| G(a), weight block only: cosine | {np.median(R['cos_jw'][:q]):.4f} | {np.median(R['cos_jw'][-q:]):.4f} |",
             f"| G(a), geometry block only: cosine | {np.median(R['cos_jg'][:q]):.4f} | {np.median(R['cos_jg'][-q:]):.4f} |",
             f"| G(a): share of the residual ‖v − v_pred‖²_(G⁻¹) carried by the weight block | {np.median(R['res_w'][:q]):.3f} | {np.median(R['res_w'][-q:]):.3f} |",
             f"| G(a): predicted weight-block / geometry-block speed, G⁻¹ norm | {np.median(R['pred_w'][:q]):.3g} | {np.median(R['pred_w'][-q:]):.3g} |",
             f"| G(b): cosine(geometry velocity, −(mδ/k_e)Γ∇L*) | {np.median(R['cos'][:q]):.3f} | {np.median(R['cos'][-q:]):.3f} |",
             f"| G(b): magnitude ratio | {np.median(R['ratio'][:q]):.3f} | {np.median(R['ratio'][-q:]):.3f} |",
             f"| relative gap of the weights to a*(c,s) | {np.median(R['a_gap'][:q]):.3f} | {np.median(R['a_gap'][-q:]):.3f} |", "",
             "| quantity | first sample | last sample |", "|---|---|---|",
             f"| profile loss L* = min_a L̄ | {R['Lstar'][0]:.4e} | {R['Lstar'][-1]:.4e} |",
             f"| L̄ at the period-mean weights | {R['Lbar'][0]:.4e} | {R['Lbar'][-1]:.4e} |",
             f"| ‖∇L*‖ | {R['gnorm'][0]:.3e} | {R['gnorm'][-1]:.3e} |",
             f"| tracking RMS over a period | {R['e_rms'][0]:.4f} | {R['e_rms'][-1]:.4f} |",
             "", f"Sampled once per period (centred differences over ±{win} periods): L* non-increasing: "
                 f"**{bool(np.all(np.diff(R['Lstar']) <= 1e-6))}**; L̄ (period-mean weights) non-increasing: "
                 f"**{bool(np.all(np.diff(R['Lbar']) <= 1e-6))}** (fraction of decreasing steps {np.mean(np.diff(R['Lbar']) <= 1e-6):.2f}); "
                 f"fraction of samples with G(b) cosine > 0.9: {np.mean(R['cos'] > 0.9):.2f}; with G(a) cosine > 0.99: {np.mean(R['cos_j'] > 0.99):.2f}"]
    with open(os.path.join(OUT, f'table_{name}.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.5))
    axs[0].plot(R['t'], R['cos'], 'b.', ms=3); axs[0].axhline(1, color='k', lw=.5); axs[0].set_ylim(-1.05, 1.05)
    axs[0].set_xlabel('t'); axs[0].set_title(f'[{name}] cosine(actual, −(m/k_e)Γ∇L*)'); axs[0].grid(alpha=.3)
    axs[1].semilogy(R['t'], R['Lstar'], 'k', label='profile loss L*'); axs[1].semilogy(R['t'], R['gnorm'], 'r--', label='‖∇L*‖')
    axs[1].legend(); axs[1].set_xlabel('t'); axs[1].grid(alpha=.3)
    for i in range(M):
        axs[2].plot(t, P[:, i, 1], lw=1)
    axs[2].set_xlabel('t'); axs[2].set_ylabel('centres c_i(t)'); axs[2].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, f'figure_{name}.png'), dpi=140); print('saved figure')


if __name__ == "__main__":
    for name in (sys.argv[1:] or list(CONFIGS)):
        run(name)
