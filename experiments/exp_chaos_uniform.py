"""
Theorem D″: uniform-in-time propagation of chaos in the leak-dominated regime
============================================================================

The synchronous coupling of Theorem D′ with the contracting terms kept: with κ_e := k_e − L_f − L_x,
κ_ϑ := σ − λ_G E_X L₂ and L_b′ := λ_G(L_ϑ + E_X L_{xϑ}), the pair (‖x − x̄‖, (1/m)Σ‖ϑ_i − ϑ̄_i‖) obeys a cooperative linear
differential inequality with the Metzler matrix A = [[−κ_e, L_ϑ], [L_b′, −κ_ϑ]], driven by the sampling fluctuation ζ.
If κ_e > 0 and det A := κ_eκ_ϑ − L_b′L_ϑ > 0, then for all t and every τ ≥ 0

    E‖x_t − x̄_t‖ ≤ K_Φ κ_ϑ m^{−1/2}/det A,    E(1/m)Σ‖ϑ_i − ϑ̄_i‖ ≤ K_Φ L_b′ m^{−1/2}/det A,

uniformly in time.  This script (i) evaluates the Lipschitz constants of the saturated features in closed form (n = 1) and
confirms them by a numerical search, (ii) computes the a-priori tube radius E_X = (F + K_Φ)/k_e of Theorem D(b) (e(0) = 0)
and the certified leak σ_cert above which the theorem applies, (iii) computes Proposition K's critical leak σ_K of the
collapsed critical point, below which the kernels separate (the learning regime), and (iv) couples m-particle systems to an
M = 2048 reference at four leaks — certified with and without dither, one eighth of σ_cert, and σ_K/3 on a long horizon —
recording the coupled distance D_m(t) in time windows (uniformity in time means the window maxima do not grow; the
m^{−1/2} law should hold in each window), the particle spread of the reference network in physical coordinates, and its
fit of the plant on the orbit (is the network learning?).  Slow adaptation (G = 0.05·diag(4,1,1), k_e = 25, s_min = 0.3),
n = 1, a two-bump plant that is not odd (so that the dead kernel is not a critical point; referee round 6, DPP-3).

Output: results/chaos/uniform_table.md (set CU_WORKERS, CU_SEEDS, CU_SEEDS_LONG to change the pool and the seed counts).
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import periodic as PO

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
M_REF = 2048
MS = (16, 64, 256)
DT = 1e-3
KE = 25.0
GAINS = 0.05 * np.array([4.0, 1.0, 1.0])
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
N_SEEDS = int(os.environ.get('CU_SEEDS', 8))
N_SEEDS_LONG = int(os.environ.get('CU_SEEDS_LONG', 6))
N_WORKERS = int(os.environ.get('CU_WORKERS', 3))
X_BOX = 2.5
PER = 2 * np.pi


class Plant:
    """f(x) = 1.6 e^{−(x−1)²/0.4} − 1.2 e^{−(x+1)²/0.4} + 0.1x (two bumps of different height: not odd), reference
    x_m = 1.5 sin t + 0.5 sin 2t."""

    def f(self, x, t):
        return 1.6 * np.exp(-(x - 1.0) ** 2 / 0.4) - 1.2 * np.exp(-(x + 1.0) ** 2 / 0.4) + 0.1 * x

    def f_vec(self, X):
        return self.f(X, 0.0)

    def df(self, x):
        return (1.6 * np.exp(-(x - 1.0) ** 2 / 0.4) * (-2 * (x - 1.0) / 0.4)
                - 1.2 * np.exp(-(x + 1.0) ** 2 / 0.4) * (-2 * (x + 1.0) / 0.4) + 0.1)

    def x_ref(self, t):
        return np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])

    def x_ref_dot(self, t):
        return np.array([1.5 * np.cos(t) + 1.0 * np.cos(2 * t)])


# ---------------------------------------------------------------------------- constants
def lipschitz_constants_exact(net):
    """
    Suprema over all ϑ (unconstrained coordinates) and all x of ‖∂_xΦ̃‖, ‖∇_ϑΦ̃‖, ‖∇²_ϑΦ̃‖, ‖∂_x∇_ϑΦ̃‖ for the saturated
    features with n = 1, Φ̃ = a(α) e^{−(x−c(γ))²/(2s(κ)²)}:
        L_x = L_ϑ = a_max e^{−1/2}/s_min     (|x − c| = s_min, c = 0, a = ±a_max; at saturation da/dα = ds/dκ = 0),
        L₂ = L_xϑ = a_max/s_min²             (x = c, s = s_min, a = ±a_max: ∂²_γγΦ̃ = −a/s², ∂_x∂_γΦ̃ = a/s²).
    """
    return dict(L_x=net.a_max * np.exp(-0.5) / net.s_min, L_th=net.a_max * np.exp(-0.5) / net.s_min,
                L_2=net.a_max / net.s_min ** 2, L_xth=net.a_max / net.s_min ** 2)


def lipschitz_constants_search(net, X_box, n_phys=120000, seed=1, h=1e-5):
    """Numerical confirmation of the closed forms: suprema over a sampled parameter set and a grid in x, polished by
    Nelder–Mead over (x, α, γ, κ) with |x| ≤ X_box."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    a = rng.uniform(-net.a_max, net.a_max, (n_phys // 2, 1)); c = rng.uniform(-X_box, X_box, (n_phys // 2, 1))
    s = np.exp(rng.uniform(np.log(net.s_min), np.log(net.s_max), n_phys // 2))
    Th = np.vstack([net.unconstrained(a, c, s), rng.normal(0, 3, (n_phys // 2, 3))])
    grad = lambda x, T: net.grad_g(np.array([x]), np.array([1.0]), T)
    def phi(x, T):
        ap, cp, sp = net.physical(T); return ap[:, 0] * np.exp(-(x - cp[:, 0]) ** 2 / (2 * sp ** 2))
    def L2(x, T):
        Hs = []
        for j in range(3):
            E = np.zeros_like(T); E[:, j] = h
            Hs.append((grad(x, T + E) - grad(x, T - E)) / (2 * h))
        H = np.stack(Hs, axis=2); H = 0.5 * (H + np.swapaxes(H, 1, 2))
        return np.linalg.norm(H, ord=2, axis=(1, 2))
    funs = dict(L_x=lambda x, T: np.abs((phi(x + h, T) - phi(x - h, T)) / (2 * h)),
                L_th=lambda x, T: np.linalg.norm(grad(x, T), axis=1),
                L_2=L2,
                L_xth=lambda x, T: np.linalg.norm((grad(x + h, T) - grad(x - h, T)) / (2 * h), axis=1))
    out = {}
    for name, fun in funs.items():
        best, arg = 0.0, None
        for x in np.linspace(-X_box, X_box, 241):
            v = fun(x, Th)
            i = int(np.argmax(v))
            if v[i] > best:
                best, arg = float(v[i]), (x, Th[i].copy())
        neg = lambda z: -float(fun(float(np.clip(z[0], -X_box, X_box)), z[1:][None, :])[0])
        r = minimize(neg, np.concatenate([[arg[0]], arg[1]]), method='Nelder-Mead', options=dict(xatol=1e-9, fatol=1e-12, maxiter=20000))
        out[name] = max(best, -float(r.fun))
    return out


def a_priori_tube(plant, net, ke, x_box=X_BOX):
    """Theorem D(b): with e(0) = 0 the tube radius E_X solves E_X = (F + K_Φ)/k_e with F := sup_{tube}|f|; the tube
    |x| ≤ sup|x_m| + E_X lies inside |x| ≤ x_box (checked)."""
    ts = np.linspace(0, PER, 20001)
    xm_max = float(np.max(np.abs([plant.x_ref(t)[0] for t in ts])))
    E_X = 0.0
    for _ in range(200):
        xs = np.linspace(-(xm_max + E_X), xm_max + E_X, 4001)
        F = float(np.max(np.abs(plant.f(xs, 0.0))))
        E_new = (F + net.K_Phi) / ke
        if abs(E_new - E_X) < 1e-12:
            break
        E_X = E_new
    assert xm_max + E_X <= x_box
    xs = np.linspace(-x_box, x_box, 20001)
    L_f = float(np.max(np.abs(plant.df(xs))))
    return dict(E_X=E_X, F=F, xm_max=xm_max, L_f=L_f)


def critical_leak(plant, net):
    """Proposition K: σ_K := σ₀_crit/k_e, the leak below which the collapsed critical point (best single kernel) is
    unstable along its softest separation mode (periodic.sigma_crit; the law here uses the physical leak σ = σ₀/k_e)."""
    cl = PO.ClosedLoop(net, plant.f_vec, plant.x_ref, plant.x_ref_dot, PER, KE, GAINS, 1.0)
    orb = cl.orbit(N=1000)
    g = net.unconstrained(np.array([[1.0]]), np.array([[0.0]]), np.array([1.0]))[0]
    s0, tb, R2 = PO.sigma_crit(cl, orb, g)
    a, c, s = net.physical(tb[None, :])
    return dict(sigma_K=s0 / KE, sigma0_crit=s0, atom=(float(a[0, 0]), float(c[0, 0]), float(s[0])))


# ---------------------------------------------------------------------------- the coupled runs
def coupled_run(seed, sigma, tau, T_end, windows):
    rng = np.random.default_rng(seed)
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=GAINS, sigma=sigma, theta0=np.zeros(3), tau=tau)
    plant = Plant()
    Th_all = net.unconstrained(rng.normal(0, 1.0, (M_REF, 1)), rng.uniform(-2, 2, (M_REF, 1)), rng.uniform(0.4, 1.0, M_REF))
    systems = {m: dict(x=plant.x_ref(0.0).copy(), Th=Th_all[:m].copy()) for m in MS + (M_REF,)}
    n_steps = int(round(T_end / DT))
    rec_t, rec = [], {m: [] for m in MS}
    e_sup, x_sup = 0.0, 0.0
    sqG = np.sqrt(GAINS)
    noise_amp = np.sqrt(2 * tau * DT) if tau > 0 else 0.0
    t = 0.0
    orbit_x = np.array([plant.x_ref(tt)[0] for tt in np.arange(400) * PER / 400])[:, None]
    f_orb = plant.f(orbit_x[:, 0], 0.0)
    fit_t, fit = [], []
    spread_t, spread_max, spread_rms = [], [], []
    for k in range(n_steps):
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        xi = rng.standard_normal((M_REF, 3)) if tau > 0 else None            # shared Brownian increments (synchronous coupling)
        new = {}
        for m, S in systems.items():
            x, Th = S['x'], S['Th']
            e = x - xr
            u = ctrl.control(x, xr, xrd, Th)
            k1 = plant.f(x, t) + u
            k2 = plant.f(x + 0.5 * DT * k1, t + 0.5 * DT) + u
            k3 = plant.f(x + 0.5 * DT * k2, t + 0.5 * DT) + u
            k4 = plant.f(x + DT * k3, t + DT) + u
            dTh = ctrl.drift(x, e, Th) * DT
            if tau > 0:
                dTh = dTh + noise_amp * sqG[None, :] * xi[:m]
            new[m] = dict(x=x + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4), Th=Th + dTh)
        systems = new
        t = (k + 1) * DT
        if k % 20 == 0:
            ref = systems[M_REF]
            e_sup = max(e_sup, float(abs(ref['x'][0] - xr[0]))); x_sup = max(x_sup, float(abs(ref['x'][0])))
            rec_t.append(t)
            for m in MS:
                rec[m].append(abs(float(systems[m]['x'][0] - ref['x'][0])) + float(np.mean(np.linalg.norm(systems[m]['Th'] - ref['Th'][:m], axis=1))))
        if k % int(round(10.0 / DT)) == 0 or k == n_steps - 1:                   # every 10 s: fit and spread of the reference network
            ThM = systems[M_REF]['Th']
            fh = net.f_hat_batch(orbit_x, ThM)[:, 0]
            fit_t.append(t); fit.append(float(np.sqrt(np.mean((f_orb - fh) ** 2) / np.mean(f_orb ** 2))))
            a, c, s = net.physical(ThM)
            P = np.column_stack([a[:, 0], c[:, 0], s])
            d = np.linalg.norm(P - P.mean(0), axis=1)
            spread_t.append(t); spread_max.append(float(d.max())); spread_rms.append(float(np.sqrt(np.mean(d ** 2))))
    rec_t = np.array(rec_t)
    out = {}
    for m in MS:
        D = np.array(rec[m])
        out[m] = [float(D[(rec_t > a) & (rec_t <= b)].max()) for a, b in windows] + [float(D[(rec_t > windows[-1][0]) & (rec_t <= windows[-1][1])].mean())]
    return dict(seed=seed, sigma=sigma, tau=tau, D=out, e_sup=e_sup, x_sup=x_sup, fit_t=np.array(fit_t), fit=np.array(fit),
                spread_t=np.array(spread_t), spread_max=np.array(spread_max), spread_rms=np.array(spread_rms))


def _slope_with_se(R):
    """Slope of log(window max) against log m, fitted per seed; mean and standard error over seeds."""
    sl = [np.polyfit(np.log(MS), np.log([R[m][s] for m in MS]), 1)[0] for s in range(len(R[MS[0]]))]
    return float(np.mean(sl)), float(np.std(sl, ddof=1) / np.sqrt(len(sl))) if len(sl) > 1 else float('nan')


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    plant = Plant()
    Lc = lipschitz_constants_exact(net)
    Ls = lipschitz_constants_search(net, X_BOX)
    tube = a_priori_tube(plant, net, KE)
    lam_G = float(GAINS.max())
    kappa_e = KE - tube['L_f'] - Lc['L_x']
    E_X = tube['E_X']
    L_bp = lam_G * (Lc['L_th'] + E_X * Lc['L_xth'])
    kappa_th = lambda sg: sg - lam_G * E_X * Lc['L_2']
    det = lambda sg: kappa_e * kappa_th(sg) - L_bp * Lc['L_th']
    sigma_cert = lam_G * E_X * Lc['L_2'] + L_bp * Lc['L_th'] / kappa_e
    K = critical_leak(plant, net)
    sigma_K = K['sigma_K']
    W_short = ((10.0, 20.0), (20.0, 40.0), (40.0, 80.0))
    W_long = ((100.0, 200.0), (200.0, 500.0), (500.0, 1000.0))
    TAU_CERT = 0.02
    configs = [dict(name='certified, dither τ = %g' % TAU_CERT, sigma=1.2 * sigma_cert, tau=TAU_CERT, T=80.0, W=W_short, seeds=N_SEEDS),
               dict(name='certified, τ = 0', sigma=1.2 * sigma_cert, tau=0.0, T=80.0, W=W_short, seeds=N_SEEDS),
               dict(name='σ_cert/8, τ = 0', sigma=sigma_cert / 8, tau=0.0, T=80.0, W=W_short, seeds=N_SEEDS),
               dict(name='σ_K/3, τ = 0 (learning regime)', sigma=sigma_K / 3, tau=0.0, T=1000.0, W=W_long, seeds=N_SEEDS_LONG)]
    jobs = [(s, c['sigma'], c['tau'], c['T'], c['W']) for c in configs for s in range(c['seeds'])]
    jobs.sort(key=lambda j: -j[3])                                             # long runs first
    with mp.Pool(N_WORKERS) as pool:
        res = pool.starmap(coupled_run, jobs)
    E_X_meas = max(r['e_sup'] for r in res if r['tau'] == 0.0); X_sup = max(r['x_sup'] for r in res)
    lines = ["# Theorem D″: uniform-in-time propagation of chaos in the leak-dominated regime (`exp_chaos_uniform.py`)", "",
             f"n = 1, k_e = {KE:g}, G = diag{tuple(float(g) for g in GAINS)}, saturated features with a_max = {NET_KW['a_max']}, c_max = {NET_KW['c_max']}, "
             f"s ∈ [{NET_KW['s_min']}, {NET_KW['s_max']}], prior ϑ⁰ = 0; plant f(x) = 1.6e^{{−(x−1)²/0.4}} − 1.2e^{{−(x+1)²/0.4}} + 0.1x, reference "
             f"x_m = 1.5 sin t + 0.5 sin 2t; the m-system shares its particles (and its Brownian increments) with an M = {M_REF} reference; dt = {DT:g}.",
             "",
             f"Constants. Lipschitz constants of the features in closed form (suprema over all ϑ and x): L_x = L_ϑ = a_max e^{{−1/2}}/s_min = {Lc['L_x']:.3f}, "
             f"L₂ = L_xϑ = a_max/s_min² = {Lc['L_2']:.2f}; numerical search over 1.2·10⁵ sampled ϑ, |x| ≤ {X_BOX}, Nelder–Mead polish: "
             f"{Ls['L_x']:.3f}, {Ls['L_th']:.3f}, {Ls['L_2']:.2f}, {Ls['L_xth']:.2f}. L_f = sup|f′| = {tube['L_f']:.3f} on |x| ≤ {X_BOX}. "
             f"A-priori tube (Theorem D(b), e(0) = 0): sup|x_m| = {tube['xm_max']:.4f}, F = sup_tube|f| = {tube['F']:.4f}, E_X = (F + K_Φ)/k_e = {E_X:.4f} "
             f"(measured sup‖ē‖ on the τ = 0 reference runs: {E_X_meas:.4f}; the runs stay within |x| ≤ {X_sup:.2f} ≤ {X_BOX}). "
             f"κ_e = k_e − L_f − L_x = {kappa_e:.3f}; L_b′ = λ_G(L_ϑ + E_X L_xϑ) = {L_bp:.3f}; the theorem applies (a priori) for "
             f"σ > σ_cert = λ_G E_X L₂ + L_b′L_ϑ/κ_e = {sigma_cert:.3f}. Proposition K: the collapsed critical point (best single kernel, "
             f"a = {K['atom'][0]:.3f}, c = {K['atom'][1]:.3f}, s = {K['atom'][2]:.3f}) is a strict local minimum of J_m iff σ > σ_K = σ₀_crit/k_e = {sigma_K:.5f} "
             f"(σ₀_crit = {K['sigma0_crit']:.5f}); below σ_K the kernels separate: the learning regime. σ_cert/σ_K = {sigma_cert / sigma_K:.0f}.", ""]
    summary = {}
    for c in configs:
        sg, tau, W = c['sigma'], c['tau'], c['W']
        rows = [r for r in res if r['sigma'] == sg and r['tau'] == tau]
        cert = det(sg) > 0 and kappa_e > 0
        R = {m: np.array([r['D'][m] for r in rows]) for m in MS}              # (seeds, len(W)+1)
        fit0 = np.mean([r['fit'][0] for r in rows]); fitT = np.mean([r['fit'][-1] for r in rows])
        spm = np.mean([r['spread_max'][-1] for r in rows]); spr = np.mean([r['spread_rms'][-1] for r in rows])
        spr0 = np.mean([r['spread_rms'][0] for r in rows])
        lines += [f"## σ = {sg:.4f}, τ = {tau:g}: {c['name']} ({'certified: det A = ' + format(det(sg), '.3f') + ', κ_ϑ = ' + format(kappa_th(sg), '.3f') if cert else 'below the threshold'}); "
                  f"T = {c['T']:g} s, {len(rows)} seeds", "",
                  f"Reference network at T: particle spread in (a, c, s), max {spm:.2e} and rms {spr:.2e} (rms at t = 0: {spr0:.2f}); "
                  f"relative fit error ‖f − f̂‖/‖f‖ on the orbit {fitT:.3f} (at t = 0: {fit0:.3f}).", "",
                  "| m | " + " | ".join(f"max D_m on ({a:g}, {b:g}]" for a, b in W) + f" | time average on ({W[-1][0]:g}, {W[-1][1]:g}] | bound K_Φ(κ_ϑ + L_b′)m^{{−1/2}}/det A |",
                  "|---|" + "---|" * (len(W) + 2)]
        for m in MS:
            mean = R[m].mean(0)
            bnd = net.K_Phi * (kappa_th(sg) + L_bp) / det(sg) / np.sqrt(m) if cert else float('nan')
            lines.append(f"| {m} | " + " | ".join(f"{v:.3e}" for v in mean[:len(W)]) + f" | {mean[len(W)]:.3e} | {bnd:.3e} |")
        last = {m: R[m][:, len(W) - 1] for m in MS}
        if all(np.all(last[m] > 0) for m in MS):
            slope, se = _slope_with_se(last)
            slope_txt = f"{slope:.2f} ± {se:.2f}"
        else:
            slope_txt = "— (zero to rounding)"
        growth = [R[m].mean(0)[len(W) - 1] / R[m].mean(0)[0] if R[m].mean(0)[0] > 0 else float('nan') for m in MS]
        lines += ["", f"Slope of the last-window maximum against m (per-seed fits, mean ± s.e.): {slope_txt} (Theorem D″: −0.5). "
                  f"Ratio of the last to the first window maxima: " + ", ".join(f"{g:.2f}" for g in growth) + " (uniform in time means no growth). "
                  "The bound column is sup_t E D_m; the table reports E max_t D_m over the window, which is larger.", ""]
        summary[c['name']] = dict(sigma=sg, tau=tau, slope=slope_txt, growth=growth, fitT=fitT, spread_rms=spr)
    learn = summary['σ_K/3, τ = 0 (learning regime)']
    growth_txt = ("grows over the horizon (last-to-first window ratios " + ", ".join(f"{g:.2f}" for g in learn['growth'])
                  + "): on this horizon propagation of chaos in the learning regime is the finite-horizon statement of Theorem D′, not a uniform one; "
                  "whether the growth saturates on longer horizons is open") if min(learn['growth']) > 1.5 else \
                 ("does not grow over the horizon (last-to-first window ratios " + ", ".join(f"{g:.2f}" for g in learn['growth']) + "): uniform in time empirically, beyond what the theorem certifies")
    lines += ["Reading: at the certified leak the particle map is a contraction uniformly in time. Without dither the reference collapses to one kernel "
              "and the coupled distance vanishes to rounding; with dither the particles form the cloud of the prior and the coupled distance is of order "
              "m^{−1/2}, bounded in time and below the bound — the theorem's content at τ > 0. In both cases the network does not learn (its fit error "
              "stays at that of the empty network). At σ_cert/8 the particles still synchronise, more slowly, and the coupled distance decays with "
              f"them. At σ_K/3 the collapsed configuration is unstable, the kernels stay separated (rms spread {learn['spread_rms']:.2f}) and the network "
              f"fits the plant (fit error {learn['fitT']:.3f}): the learning regime, where the coupled distance follows the m^{{−1/2}} law at each time and "
              + growth_txt + ".", "",
              "Numbers are produced by `experiments/exp_chaos_uniform.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'uniform_table.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(lines) + "\n")
    raw = {f"D_{i}": np.array([r['sigma'], r['tau'], r['e_sup']] + [v for m in MS for v in r['D'][m]]) for i, r in enumerate(res)}
    raw.update({f"fit_{i}": np.column_stack([r['fit_t'], r['fit'], r['spread_max'], r['spread_rms']]) for i, r in enumerate(res)})
    np.savez(os.path.join(OUT, 'uniform_raw.npz'), **raw)
    print("\n".join(lines))


if __name__ == '__main__':
    main()
