"""
Proposition D‴.51: the transit at a finite mass — near ties share a feature, the entrant's overshoot, the creep, the predicted mass
====================================================================================================================================

The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025; Gaussian kernels of width 0.12
on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, masses ε) at a finite mass, against the ε-free inner problem of a
near tie and the predictions made from the sliding path alone (`rich_transit`).

A. The near-tie law: one bump between the kernels at 0.0 and 0.1, δ from the midpoint (amplitudes 0.3, 0.8, 2.0; δ = 10⁻⁴, 10⁻³, 10⁻²;
   ε = 10⁻⁵ … 10⁻¹⁰): the runner-up's share of the bump against x = μ*/ε^{1/2} (μ* its margin when the entrant enters), against the
   inner problem's share at the same x (`near_tie_inner`, at the exact tie: the limit ε → 0 at fixed x), the Gibbs split 1/(1 + e^{κx}) (κ = πρ(2c)^{1/2}) and the
   frozen coefficient K'³ρ⁻³/(|A|(x − Δ)³) (Δ the inner problem's pull, B). A2: the inner problems not stationary within 8 doublings of
   s, two doublings later (their drift).
B. The entrant's overshoot η_c(α_A)/(Λ₀ρ) (`overshoot`) against the pull Δ(A) the inner problem exerts on a frozen runner-up (x = 200;
   from the starts XI0_B behind the wall, extrapolated) and the loop's (the runner-up's final inner coordinate read from its coefficient,
   δ = 10⁻³, ε = EPS_BL), free and with the runner-up's centre held fixed (its approach creep removed); I∞ = overshoot − Δ.
C. The creep: target 10 of Proposition D‴.48's random targets (C4, C6), whose path leaves a kernel 0.0065 inside its wall beside an
   unfitted bump — its predicted creep threshold (`recovery_threshold`) bracketed by the loop.
D. Prediction from the path, out of sample: 48 random three-bump targets (D‴.48's recipe, a fresh seed); the predicted ε_c against the loop
   at ε = 10⁻⁵, 10⁻⁶, 10⁻⁷, 10⁻⁸ (formed: |w| > 10⁻²max|A*|; recovered: exactly the path's kernels formed).
E. The predicted ε_c over 400 random targets: the fraction predicted not to recover at ε, entry- and creep-limited, and their slopes in
   log ε (local, and fitted where between 10 targets and half of them fail).
F. The predicted tail checked: every target of E with 10⁻¹² ≤ ε_c ≤ 10⁻⁷, the loop at 3ε_c (predicted to fail) and ε_c/3 (to recover).
G. The loop's own threshold (referee round 39): on D's 48 targets the largest mass at which the loop still forms exactly the path's
   kernels, bracketed from ε_c by factors G_FAC and bisected G_STEPS times in log ε, against the predicted ε_c — the distribution of
   ε_loop/ε_c, for entry- and creep-bound targets.
G2. F's targets not as predicted (E's seed): their loop thresholds, bracketed and bisected as in G.
HI. The inner problem's own formation threshold x_θ at the exact tie (A = 0.3, 0.8, 2), bisected in x — H(1)'s reference.
HP. The predictor's own threshold for H(1)'s runner-up (`recovery_threshold`), with its parts.
H. Near-tie controls: (1) the single bump (A = 0.8, δ = 10⁻³, 10⁻², 3·10⁻²): the loop's formation threshold x_θ = μ*/ε^{1/2} of the
   runner-up, bisected, free and with the runner-up's centre or both centres held fixed until the entrant's arrival τ_e (the partners'
   approach creep removed); (2) exact ties made by unequal masses (the farther kernel heavier, b from the samples' rates): the farther
   kernel's share at ε = 10⁻⁸, 10⁻¹⁰ against ½; (3) D's target 2 with its middle bump moved (MIG_C2): the first entry's runner-up margin
   μ*, the loop's threshold (bisected) against ε_c, μ*² and μ*³, and where the extra kernel ends; (4) E's target 117, a runner-up on the
   entrant's way to its feature: its largest coefficient and its closest approach to the entrant at ε = 10⁻⁶, 10⁻⁵, 10⁻⁴.
The sliding paths are relaxed with the block's least-squares polish (`sliding_iss(polish=True)`: the same paths, without the fast flow's
crawl along the block's slow modes).

Output: results/chaos/transit_mass.md and transit_mass_raw.npz (`--from-raw` regenerates the md; `--only A,A2,B,C,D,E,F,G,G2,H,HI,HP`).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('TM_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TM_WORKERS', '3'))
AMPS_A, OFFS_A = (0.3, 0.8, 2.0), (1e-4, 1e-3, 1e-2)
EPS_A = (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
TAU_A = 1.32                                                          # the run to 1.32 τ_entry (the relaxation after it)
AMPS_B = (0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 4.0)
X_B, OFF_B, EPS_B = 200.0, 1e-3, 1e-10
EPS_BL = (1e-8, 1e-9, 1e-10)                                          # the loop's pull at x = μ*/ε^{1/2} ≈ 35, 110, 347
XI0_B = (400.0, 800.0, 1600.0, 3200.0)                                # B: the inner problem's starts behind the wall
CREEP_TARGET, EPS_C = 10, (1e-9, 2e-9, 3.5e-9, 7e-9, 1.4e-8, 3e-8)
N_D, SEED_D, EPS_D = 48, 51, (1e-5, 1e-6, 1e-7, 1e-8)
N_E, SEED_E = 400, 52
EPS_GRID_E = tuple(10.0 ** -k for k in np.arange(4, 12.5, 0.5))
F_RANGE, F_FAC = (1e-12, 1e-7), 3.0
G_FAC, G_STEPS = 4.0, 4                                               # G: the bracket's factor; bisections (the bracket ×1.09)
A_H, OFFS_H, H_STEPS = 0.8, (1e-3, 1e-2, 3e-2), 12                     # H(1): the single bump's near ties
EPS_H2 = (1e-8, 1e-10)                                                 # H(2): unequal-mass exact ties
MIG_N, MIG_C2 = 2, (0.030, 0.035, 0.040, 0.043, 0.046, 0.048)          # H(3): D's target 2, its middle bump moved
E117, EPS_H4 = 117, (1e-6, 1e-5, 1e-4)                                 # H(4): a runner-up on the entrant's way


def setup():
    warnings.filterwarnings('ignore')
    L = X.law()
    K = RT.Bumps(np.linspace(-1, 1, X.NS_C), X.SIG_C)
    c0 = np.linspace(-1, 1, X.D_C)
    return L, K, c0


def random_targets(seed, n):
    """D‴.48's recipe: three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
        out.append((cs, As))
    return out


def near_tie_case(A, off):
    L, K, c0 = setup()
    cs = 0.05 + off
    f = X.bumps_target(K, (cs,), (A,))
    g = K.phi(c0).T @ f / K.xs.size
    order = np.argsort(-np.abs(g))
    w_, n_ = int(order[0]), int(order[1])
    return L, K, c0, cs, f, g, w_, n_, 1 - abs(g[n_]) / abs(g[w_])


def a_run(args):
    A, off, eps = args
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A, off)
    t0 = time.time()
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, TAU_A / abs(g[w_]), 'paradigm', n_t=331)
    x = mu / np.sqrt(eps)
    inner = RT.near_tie_inner(L, K, A, 0.5 * (c0[w_] + c0[n_]), c0[w_], c0[n_], x)     # the limit ε → 0 at fixed x: the exact tie
    wn, ww = abs(W[-1][n_]), abs(W[-1][w_])
    return dict(A=A, off=off, eps=eps, mu=float(mu), x=float(x), entrant=w_, runner=n_, share=float(wn / (ww + wn)), w_r=float(wn),
                c_r=float(cc[-1][n_]), c_e=float(cc[-1][w_]), share_inner=float(inner['share']), inner_ok=bool(inner['stationary']),
                secs=time.time() - t0)


def part_A():
    from multiprocessing import Pool
    jobs = [(A, off, eps) for A in AMPS_A for off in OFFS_A for eps in EPS_A]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(a_run, jobs)
    return rows


def part_A2(rows):
    """The inner problems not stationary within 8 doublings of s: the share two doublings later (s four times longer) — their drift."""
    L, K, c0 = setup()
    out = []
    for r in rows:
        if r.get('inner_ok', True):
            continue
        w_, n_ = r['entrant'], r['runner']
        inner = RT.near_tie_inner(L, K, r['A'], 0.5 * (c0[w_] + c0[n_]), c0[w_], c0[n_], r['x'], max_doublings=10)
        out.append(dict(A=r['A'], off=r['off'], eps=r['eps'], x=r['x'], share8=r['share_inner'], share10=float(inner['share']),
                        s10=float(inner['s']), ok10=bool(inner['stationary'])))
    return out


def b_run(A):
    """The pull on a frozen runner-up in the inner problem — matched: its outer margin at the entrant's arrival less its final ξ — at the
    exact tie and at the offset δ, from the starts XI0_B behind the wall (the entrant's far-field creep makes the pull converge like
    a/xi0 + b/xi0²: extrapolated to xi0 → ∞ by second-order Richardson from the last three, referee round 39), against the loop's
    (x − ξ_r(∞), x the runner-up's margin at the entry over ε^{1/2}), free and with the runner-up's centre held fixed (its creep removed)."""
    L, K, c0 = setup()
    ov = float(RT.overshoot(L, A))
    tie = [RT.near_tie_inner(L, K, A, 0.05, 0.1, 0.0, X_B, both=False, xi0=x0) for x0 in XI0_B]
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A, OFF_B)
    xs, dl, dlf = [], [], []
    d = c0.size
    for eps in EPS_BL:
        taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, TAU_A / abs(g[w_]), 'paradigm', n_t=331)
        xs.append(float(mu / np.sqrt(eps)))
        dl.append(float(xs[-1] - RT.inner_xi(L, abs(W[-1][n_]))))
        b0 = np.full(d, eps)                                           # the runner-up's centre held fixed: its creep removed
        base = RT.critical_rhs(L, K, f, 'paradigm', b0=b0)

        def rhs(t, z, base=base):
            v = base(t, z)
            v[d + n_] = 0.0
            return v
        sol = solve_ivp(rhs, (0, RT.transit_scale(L, eps) * TAU_A / abs(g[w_])), np.r_[np.zeros(d), c0], method='LSODA', rtol=1e-9,
                        atol=1e-12)
        dlf.append(float(xs[-1] - RT.inner_xi(L, abs(RT.critical_w(L, sol.y[:d, -1], b0)[n_]))))
    # the runner-up's own creep while it waits at its margin (D‴.48(iii)(a)): ∫ of its creep times ⟨f*∂_cφ⟩ over its approach,
    # (a²b/(π²c))(g′_r/g_r)²/x in units of ε^{1/2} at leading order (b, not b^{3/2}: dμ/dτ = −g/ρ carries a ρ; b = 1 here)
    gr = float(K.phi(np.array([c0[n_]]))[:, 0] @ f) / K.xs.size
    gp = float(K.dphi(np.array([c0[n_]]))[:, 0] @ f) / K.xs.size
    creep = [L.A ** 2 / (np.pi ** 2 * L.c) * (gp / gr) ** 2 / x for x in xs]
    off = [RT.near_tie_inner(L, K, A, cs, c0[w_], c0[n_], X_B, both=False, xi0=x0) for x0 in XI0_B]
    ext = lambda v: (v[-3] - 6 * v[-2] + 8 * v[-1]) / 3                # second-order Richardson at xi0, 2xi0, 4xi0
    pt, po = [r['pull'] for r in tie], [r['pull'] for r in off]
    return dict(A=A, overshoot=ov, pull_tie=pt, pull_off=po, delta_inner=float(ext(pt)), delta_inner_off=float(ext(po)),
                xi_e=tie[-1]['xi_e'], x_loop=xs[-1], delta_loop=dl[-1], x_loops=xs, delta_loops=dl, delta_frozen=dlf, creep=creep,
                inner_ok=bool(all(r['stationary'] for r in tie + off)))


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_run, list(AMPS_B))


def c_run(args):
    eps, cs, As = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    te = 1.3 * R['events'][-1][0]
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, te, 'paradigm')
    sm = X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())
    path = X.net_entries(R)
    return dict(eps=eps, n_act=sm['n_act'], order=sm['order'], ok=bool(sm['n_act'] == len(path) and set(sm['order']) == set(path)),
                disp=float(np.abs(cc[-1] - c0).max()))


def part_C():
    from multiprocessing import Pool
    L, K, c0 = setup()
    cs, As = random_targets(21, X.N_RAND)[CREEP_TARGET]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    pred = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    with Pool(N_WORKERS) as pool:
        rows = pool.map(c_run, [(e, cs, As) for e in EPS_C])
    return dict(cs=cs.tolist(), As=As.tolist(), path=X.net_entries(R), pred=pred, rows=rows)


def d_pred(args):
    n, cs, As = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    pred = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    return dict(n=n, cs=cs.tolist(), As=As.tolist(), path=X.net_entries(R), slide_n=int(R['act'].sum()),
                slide_rms=float(R['hist'][-1][3]), te=float(1.3 * R['events'][-1][0]), pred=pred)


def d_loop(args):
    n, cs, As, eps, te, path = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, te, 'paradigm')
    sm = X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())
    return dict(n=n, eps=eps, n_act=sm['n_act'], order=sm['order'], ok=bool(sm['n_act'] == len(path) and set(sm['order']) == set(path)))


def part_D():
    from multiprocessing import Pool
    tg = random_targets(SEED_D, N_D)
    with Pool(N_WORKERS) as pool:
        preds = pool.map(d_pred, [(n, cs, As) for n, (cs, As) in enumerate(tg)])
        jobs = [(p['n'], np.array(p['cs']), np.array(p['As']), eps, p['te'], p['path']) for p in preds for eps in EPS_D]
        loops = pool.map(d_loop, jobs)
    return dict(preds=preds, loops=loops)


def f_run(args):
    n, eps = args
    L, K, c0 = setup()
    cs, As = random_targets(SEED_E, N_E)[n]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    path = X.net_entries(R)
    r = d_loop((n, cs, As, eps, 1.3 * R['events'][-1][0], path))
    return dict(n=n, eps=eps, ok=r['ok'], n_act=r['n_act'], order=r['order'], path=path)


def part_F(Ee):
    """The loop on E's predicted tail: each target with F_RANGE[0] ≤ ε_c ≤ F_RANGE[1] at F_FAC·ε_c and ε_c/F_FAC."""
    from multiprocessing import Pool
    pe = Ee['preds']
    ec = {p['n']: min(p['eps_entry'], p['eps_creep']) for p in pe}
    sel = sorted(n for n, v in ec.items() if F_RANGE[0] <= v <= F_RANGE[1])
    jobs = [(n, F_FAC * ec[n]) for n in sel] + [(n, ec[n] / F_FAC) for n in sel]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(f_run, jobs)
    return dict(rows=rows, sel=sel, below=int(sum(1 for v in ec.values() if v < F_RANGE[0])))


def part_E():
    from multiprocessing import Pool
    tg = random_targets(SEED_E, N_E)
    with Pool(N_WORKERS) as pool:
        preds = pool.map(d_pred, [(n, cs, As) for n, (cs, As) in enumerate(tg)])
    return dict(preds=[dict(n=p['n'], slide_n=p['slide_n'], eps_entry=p['pred']['eps_entry'], eps_creep=p['pred']['eps_creep'])
                       for p in preds])


def bisect_loop(ok, lo, hi, steps):
    """The transition of ok(ε) (True below, False above) in [lo, hi], bisected in log ε `steps` times: (lo, hi)."""
    a, b = np.log(lo), np.log(hi)
    for _ in range(steps):
        m = 0.5 * (a + b)
        if ok(np.exp(m)):
            a = m
        else:
            b = m
    return float(np.exp(a)), float(np.exp(b))


def g_job(args):
    """The loop's own threshold on one of D's targets: bracketed from the prediction ε_c by factors G_FAC (up to 10⁻³, down to 10⁻¹³), then
    bisected; (lo, hi) with the loop recovering at lo and not at hi; hi = inf if it recovers up to 10⁻³."""
    n, cs, As, ec, te, path = args
    ok = lambda e: d_loop((n, cs, As, e, te, path))['ok']
    lo, hi = ec / G_FAC, ec * G_FAC
    ok_hi = ok(hi)
    while ok_hi and hi < 1e-3:
        lo, hi = hi, hi * G_FAC
        ok_hi = ok(hi)
    if ok_hi:
        return dict(n=n, ec=ec, lo=float(hi), hi=np.inf)
    ok_lo = ok(lo)
    while not ok_lo and lo > 1e-13:
        hi, lo = lo, lo / G_FAC
        ok_lo = ok(lo)
    if not ok_lo:
        return dict(n=n, ec=ec, lo=0.0, hi=float(lo))
    a, b = bisect_loop(ok, lo, hi, G_STEPS)
    return dict(n=n, ec=ec, lo=a, hi=b)


def part_G(D):
    from multiprocessing import Pool
    jobs = [(p['n'], np.array(p['cs']), np.array(p['As']), min(p['pred']['eps_entry'], p['pred']['eps_creep']), p['te'], p['path'])
            for p in D['preds']]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(g_job, jobs, chunksize=1)
    return dict(rows=rows, bind={p['n']: 'entry' if p['pred']['eps_entry'] <= p['pred']['eps_creep'] else 'creep' for p in D['preds']})


def g2_job(n):
    """G's bisection on one of E's targets: its path, the prediction and the loop's threshold."""
    cs, As = random_targets(SEED_E, N_E)[n]
    p = d_pred((n, cs, As))
    r = g_job((n, np.array(cs), np.array(As), min(p['pred']['eps_entry'], p['pred']['eps_creep']), p['te'], p['path']))
    r['bind'] = 'entry' if p['pred']['eps_entry'] <= p['pred']['eps_creep'] else 'creep'
    return r


def part_G2(Ee, F):
    """F's targets not as predicted (recovered at F_FAC·ε_c, or not at ε_c/F_FAC): their loop thresholds as in G."""
    from multiprocessing import Pool
    pe = {p['n']: p for p in Ee['preds']}
    ec = lambda n: min(pe[n]['eps_entry'], pe[n]['eps_creep'])
    bad = sorted({r['n'] for r in F['rows'] if (r['eps'] > ec(r['n'])) == r['ok']})
    with Pool(N_WORKERS) as pool:
        return dict(rows=pool.map(g2_job, bad, chunksize=1))


def h_formed(off, frz, eps, A=A_H):
    """H(1): whether the single bump's runner-up (A, offset δ = off) forms in the loop at the mass ε, the centres in `frz` ('free',
    'runner-up', 'both') held fixed until the entrant's arrival on the path, τ_e = 1/|⟨fφ_e⟩|."""
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A, off)
    d = c0.size
    fix = {'free': [], 'runner-up': [n_], 'both': [w_, n_]}[frz]
    te = 1.0 / abs(g[w_])
    b0 = np.full(d, eps)
    Lam = RT.transit_scale(L, eps)
    base = RT.critical_rhs(L, K, f, 'paradigm', b0=b0)

    def rhs(t, z):
        v = base(t, z)
        if t < Lam * te:
            v[d + np.array(fix, int)] = 0.0
        return v
    sol = solve_ivp(rhs, (0, Lam * TAU_A * te), np.r_[np.zeros(d), c0], method='LSODA', rtol=1e-9, atol=1e-12)
    return abs(RT.critical_w(L, sol.y[:d, -1], b0)[n_]) > X.TH_ACT * abs(A)


def h_single(args):
    """H(1): the single bump's near tie (A_H, offset δ): the runner-up's formation threshold (`h_formed`, the centres in `frz` held fixed
    until the entrant's arrival), x_θ = μ*/ε^{1/2} at the bisected threshold."""
    off, frz = args
    mu = near_tie_case(A_H, off)[-1]
    a, b = bisect_loop(lambda e: not h_formed(off, frz, e), (mu / 12) ** 2, (mu / 4) ** 2, H_STEPS)
    return dict(off=off, frz=frz, mu=float(mu), lo=a, hi=b, x_theta=float(mu / np.sqrt(a * b) ** 0.5))


def h_inner(A):
    """H(1): the inner problem's own formation threshold at the exact tie (the limit ε → 0 at fixed x; the bump at the kernels' midpoint):
    the gap x at which the runner-up ends at |w_r| = θ = TH_ACT·|A|, bisected H_STEPS times in [4, 9]."""
    L, K, c0 = setup()
    lo, hi = 4.0, 9.0
    for _ in range(H_STEPS):
        x = 0.5 * (lo + hi)
        if RT.near_tie_inner(L, K, A, 0.05, 0.1, 0.0, x, max_doublings=10)['w_r'] > X.TH_ACT * abs(A):
            lo = x
        else:
            hi = x
    return dict(A=A, x_theta=0.5 * (lo + hi), gibbs=float(np.log(1 / X.TH_ACT - 1) / (np.pi * np.sqrt(2 * L.c))))


def part_HI():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return dict(rows=pool.map(h_inner, AMPS_A, chunksize=1))


def h_pred(off):
    """H(1): the predictor's own near-tie threshold on the single bump (`recovery_threshold`): x = μ*/ε_entry^{1/2} with its parts, the
    frozen runner-up's pull Δ + ξ_θ and the near-tie estimate x_θκ_ij."""
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A_H, off)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    p = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * A_H)
    i, j, m_, kij, xth, pull = p['entry']
    return dict(off=off, mu=float(mu), runner=int(n_), j=int(j), eps_entry=float(p['eps_entry']), x=float(m_ / np.sqrt(p['eps_entry'])),
                kij=float(kij), xth=float(xth), pull_x=float(pull + RT.inner_xi(L, X.TH_ACT * A_H)))


def part_HP():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return dict(rows=pool.map(h_pred, OFFS_H, chunksize=1))


def h_masses(args):
    """H(2): an exact tie made by unequal masses: the bump at c*, the nearer kernel jn (b = 1) and the farther jf, its mass b_f set so
    that both reach their walls together, v_f/ρ_f = v_n (ρ = b^{−1/2}); the loop at the masses εb: the farther kernel's share."""
    cstar, jn, jf, eps = args
    L, K, c0 = setup()
    f = X.bumps_target(K, (cstar,), (A_H,))
    g = K.phi(c0).T @ f / K.xs.size
    b = np.ones(c0.size)
    b[jf] = (g[jn] / g[jf]) ** 2
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, TAU_A / abs(g[jn]), 'paradigm', b=b, n_t=331)
    wn, wf = abs(W[-1][jn]), abs(W[-1][jf])
    R = RT.sliding_iss(L, K, f, c0, b ** -0.5, polish=True)
    return dict(cstar=cstar, jn=jn, jf=jf, b=float(b[jf]), eps=eps, share=float(wf / (wn + wf)), path=X.net_entries(R),
                formed=[int(j) for j in np.flatnonzero(np.abs(W[-1]) > X.TH_ACT * A_H)], c_f=float(cc[-1][jf]), c_n=float(cc[-1][jn]))


def h_migration(c2):
    """H(3): D's target MIG_N with its middle bump at c2: the first entry's runner-up and its margin μ*, the predicted ε_c, the loop's
    threshold (bracketed and bisected as in G) and the kernels the loop forms just above it, with where the extra ones end."""
    L, K, c0 = setup()
    cs, As = random_targets(SEED_D, N_D)[MIG_N]
    cs = np.array(cs, float); cs[1] = c2
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    P = RT.path_duals(K, f, R)
    path, te = X.net_entries(R), 1.3 * R['events'][-1][0]
    k1 = int(np.argmin(np.abs(np.array(P['taus']) - R['events'][0][0])))
    mus = {j: 1 - abs(P['U'][k1][j]) for j in range(c0.size) if j not in P['ever']}
    jr = min(mus, key=mus.get)
    pred = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    ec = min(pred['eps_entry'], pred['eps_creep'])
    g = g_job((MIG_N, cs, As, ec, te, path))
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, g['hi'], te, 'paradigm')
    sm = X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())
    extra = [j for j in sm['order'] if j not in path]
    return dict(c2=c2, path=path, runner=int(jr), mu=float(mus[jr]), eps_entry=float(pred['eps_entry']), eps_creep=float(pred['eps_creep']),
                lo=g['lo'], hi=g['hi'], order=sm['order'], extra={int(j): float(cc[-1][j]) for j in extra})


def h_near_side(eps):
    """H(4): E's target 117 — kernel 14 between the third entrant 13 and its bump, nearer to the bump, just inside its wall: its largest
    coefficient over the run and its closest approach to the entrant."""
    L, K, c0 = setup()
    cs, As = random_targets(SEED_E, N_E)[E117]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C), polish=True)
    P = RT.path_duals(K, f, R)
    k = int(np.argmin(np.abs(np.array(P['taus']) - R['events'][-1][0])))
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, 1.3 * R['events'][-1][0], 'paradigm', n_t=881)
    mu = 1 - abs(P['U'][k][14])
    return dict(eps=eps, mu=float(mu), x=float(mu / np.sqrt(eps)), w14=float(np.abs(W[:, 14]).max()), gap=float(np.abs(cc[:, 13] - cc[:, 14]).min()),
                c14=float(cc[-1][14]), order=X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())['order'], path=X.net_entries(R),
                theta=float(X.TH_ACT * np.abs(As).max()))


def part_H():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        h1 = pool.map_async(h_single, [(off, frz) for off in OFFS_H for frz in ('free', 'runner-up', 'both')], chunksize=1)
        h2 = pool.map_async(h_masses, [(cstar, jn, jf, eps) for cstar, jn, jf in ((0.07625, 11, 10), (0.08, 11, 12)) for eps in EPS_H2],
                            chunksize=1)
        h3 = pool.map_async(h_migration, MIG_C2, chunksize=1)
        h4 = pool.map_async(h_near_side, EPS_H4, chunksize=1)
        return dict(single=h1.get(), masses=h2.get(), migration=h3.get(), near_side=h4.get())


def report(A=None, B=None, C=None, D=None, Ee=None, A2=None, F=None, G=None, H=None, G2=None, HI=None, HP=None):
    fmt = E.fmt
    L, K, c0 = setup()
    kap = np.pi * np.sqrt(2 * L.c)
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    L_ = ['# Proposition D‴.51: the transit at a finite mass (`exp_transit_mass.py`)', '',
          "The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width "
          "0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, masses ε; the inner problem, the overshoot and the "
          "predictions from `rich_transit`.", '']
    dl = {r['A']: r['delta_inner'] for r in B} if B is not None else {}
    if A is not None:
        L_ += ['## A. The near-tie law', '',
               "One bump between the kernels at 0.0 and 0.1, δ from their midpoint; μ* the runner-up's margin when the entrant enters, "
               "x = μ*/ε^{1/2}; the runner-up's share |w_r|/(|w_e| + |w_r|) of the bump in the loop (at 1.32 τ_entry) against the ε-free inner "
               "problem's at the same x (its limit ε → 0 at fixed x: the bump at the midpoint, the runner-up x behind; the loop's ratio to it → 1 "
               f"as δ ∝ ε^{{1/2}}x → 0), the Gibbs split 1/(1 + e^{{κx}}) (κ = πρ(2c)^{{1/2}} = {kap:.4f}) and the frozen "
               f"coefficient K'³ρ⁻³/(|A|(x − Δ)³) (K'³ = {Kp3:.4f}; Δ from B):", '',
               "| A | δ | ε | x | share: loop | inner problem | loop/inner | Gibbs split | frozen | the runner-up's centre |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(A, key=lambda r: (r['A'], r['x'])):
            gib = 1 / (1 + np.exp(kap * r['x']))
            dd = dl.get(r['A'])
            fro = Kp3 / (r['A'] * (r['x'] - dd) ** 3) if dd is not None and r['x'] > dd + 1 else float('nan')
            L_.append(f"| {r['A']:g} | {r['off']:g} | {r['eps']:g} | {r['x']:.4g} | {fmt(r['share'])} | {fmt(r['share_inner'])} | "
                      f"{r['share'] / r['share_inner']:.4f} | {fmt(gib)} | {fmt(fro)} | {r['c_r']:.4f} |")
        by = []
        for off in sorted({r['off'] for r in A}):
            rat = np.array([r['share'] / r['share_inner'] for r in A if r['off'] == off])
            by.append(f"{rat.min():.4f}–{rat.max():.4f} at δ = {off:g}")
        nst = sum(1 for r in A if r.get('inner_ok', True))
        dr = ''
        if A2 is not None and len(A2):
            ch = max(abs(q['share10'] / q['share8'] - 1) for q in A2)
            dr = (f"; the other {len(A2)} (x ≈ {min(q['x'] for q in A2):.3g}–{max(q['x'] for q in A2):.3g}) still drifting, their share "
                  f"changing by at most {ch:.1e} (relative) from s = {min(q['s10'] for q in A2) / 4:.3g} to {min(q['s10'] for q in A2):.3g}")
        L_ += ['', f"The loop's share over the inner problem's: {'; '.join(by)} (x from {min(r['x'] for r in A):.3g} to "
               f"{max(r['x'] for r in A):.4g}; the inner problem stationary — every rate below 10⁻⁹|A| — on {nst} of {len(A)}{dr}).", '']
    if B is not None:
        L_ += ["## B. The entrant's overshoot and the runner-up's pull", '',
               "The overshoot η_c(α_A)/(Λ₀ρ) (an entered kernel's dual beyond its wall, in units of ε^{1/2}ρ) against the pull Δ the inner "
               f"problem exerts on a frozen runner-up (x = {X_B:g}; matched: the runner-up's outer margin at the entrant's arrival less its "
               f"final ξ) at the exact tie and at δ = {OFF_B:g}, from the starts xi0 = {', '.join(f'{v:g}' for v in XI0_B)} behind the wall "
               "(the entrant's far-field creep: the pull converges like a/xi0 + b/xi0²; ∞ the second-order Richardson extrapolation from the "
               f"last three), and the loop's (δ = {OFF_B:g}, ε = {', '.join(f'{e:g}' for e in EPS_BL)}: Δ = x − ξ_r(∞), ξ_r read from the "
               "runner-up's coefficient), free and with the runner-up's centre held fixed: the difference is the runner-up's own creep while "
               "it waits at its margin, (a²b/(π²c))(g′_r/g_r)²/x in units of ε^{1/2} at leading order (predicted; their ratio in "
               "parentheses), and the fixed loop's excess over the inner problem at that offset the rest; I∞ = overshoot − Δ (tie, ∞):", '',
               "| A | overshoot | Δ: inner, tie (xi0 = " + ', '.join(f'{v:g}' for v in XI0_B) + ", ∞) | I∞ | Δ: inner, δ (∞) | Δ: loop (x) | "
               "Δ: loop, the runner-up's centre fixed | free − fixed: predicted (ratio) | fixed − inner, δ |", "|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            xs_, dls_, cr_ = r.get('x_loops', [r['x_loop']]), r.get('delta_loops', [r['delta_loop']]), r.get('creep', [np.nan])
            dlf_ = r.get('delta_frozen', [np.nan] * len(dls_))
            loops = ', '.join(f'{d:.4f} ({x:.4g})' for d, x in zip(dls_, xs_))
            fixed = ', '.join(f'{d:.4f}' for d in dlf_)
            creep_ = ', '.join(f"{a - b:.4f}: {c:.4f} ({(a - b) / c:.3f})" for a, b, c in zip(dls_, dlf_, cr_))
            rest = ', '.join(f"{b - r['delta_inner_off']:+.4f}" for b in dlf_)
            L_.append(f"| {r['A']:g} | {r['overshoot']:.4f} | {', '.join(f'{v:.4f}' for v in r['pull_tie'])}, {r['delta_inner']:.4f}"
                      f"{'' if r.get('inner_ok', True) else ' (not stationary)'} | {r['overshoot'] - r['delta_inner']:.4f} | "
                      f"{r['delta_inner_off']:.4f} | {loops} | {fixed} | {creep_} | {rest} |")
        lam0 = RT.transit_scale(L, 1.0)
        Bs = sorted(B, key=lambda r: r['A'])
        inc = [(b['overshoot'] - a['overshoot']) / np.log(b['A'] / a['A']) for a, b in zip(Bs[:-1], Bs[1:])]
        xis = (15.0, 10.0, 8.3, 5.0, 3.0)
        prof = [float(RT.inner_coef(L, x) / (Kp3 * x ** -3)) - 1 for x in xis]
        L_ += ['', "The overshoot's growth per unit log A between consecutive amplitudes: " + ', '.join(f'{v:.4f}' for v in inc) +
               f"; the exponential branch's slope a/(4GcΛ₀ρ) = {L.A / (4 * L.G * L.c * lam0):.4f} (Λ₀ = {lam0:.4f}). The inner profile against "
               "its core power law, 𝒲(ξ)/(K'³ξ⁻³) − 1 at ξ = " + ', '.join(f'{x:g}' for x in xis) + ': ' +
               ', '.join(f'{100 * v:+.1f}%' for v in prof) + '.', '']
    if C is not None:
        p = C['pred']
        L_ += ['## C. The creep', '',
               f"Target {CREEP_TARGET} of D‴.48's random targets (centres {[round(v, 3) for v in C['cs']]}, amplitudes "
               f"{[round(v, 3) for v in C['As']]}; the path's kernels {C['path']}). Predicted from the path: the entry threshold "
               f"{fmt(p['eps_entry'])} (kernel {p['entry'][1] if p['entry'] else '—'} at the entry of {p['entry'][0] if p['entry'] else '—'}, "
               f"margin {fmt(p['entry'][2]) if p['entry'] else '—'}), the creep thresholds "
               f"{', '.join(f'{j}: {fmt(v)}' for j, v in p['creep'].items())}: ε_c = {fmt(p['eps_c'])}. The loop:", '',
               "| ε | ε/ε_c | formed | entrants | the path's kernels | the inactive kernels' largest displacement |", "|---|---|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {r['eps'] / p['eps_c']:.3g} | {r['n_act']} | {r['order']} | {'yes' if r['ok'] else 'no'} | {fmt(r['disp'])} |"
               for r in C['rows']]
        L_.append('')
    if D is not None:
        preds = {p['n']: p for p in D['preds']}
        L_ += ['## D. Prediction from the path (out of sample)', '',
               f"{len(preds)} random three-bump targets (seed {SEED_D}); ε_c = min(ε_entry, ε_creep) predicted from the sliding path; the loop "
               "recovers when it forms exactly the path's kernels. Predicted (ε < ε_c) against observed:", '',
               "| ε | recovered | predicted to recover | agree | recovered though not predicted | not recovered though predicted |",
               "|---|---|---|---|---|---|"]
        miss = []
        for eps in EPS_D:
            rr = [r for r in D['loops'] if r['eps'] == eps]
            pr = [eps < preds[r['n']]['pred']['eps_c'] for r in rr]
            ob = [r['ok'] for r in rr]
            L_.append(f"| {eps:g} | {sum(ob)} of {len(ob)} | {sum(pr)} | {sum(p_ == o for p_, o in zip(pr, ob))} | "
                      f"{sum((not p_) and o for p_, o in zip(pr, ob))} | {sum(p_ and (not o) for p_, o in zip(pr, ob))} |")
            miss += [(r['n'], eps) for r, p_, o in zip(rr, pr, ob) if p_ != o]
        dis = ', '.join(f"({n}, {e:g}, {preds[n]['pred']['eps_c']:.3g})" for n, e in miss) or 'none'
        L_ += ['', f"The disagreements (target, ε, ε_c): {dis}; the "
               f"path itself exact (three kernels, rms residual < 10⁻⁸) on {sum(1 for p in preds.values() if p['slide_n'] == 3 and p['slide_rms'] < 1e-8)} "
               f"of {len(preds)}.", '']
    if Ee is not None:
        pe = Ee['preds']
        ee = np.array([p['eps_entry'] for p in pe]); ec = np.array([p['eps_creep'] for p in pe]); eps_c = np.minimum(ee, ec)
        L_ += ['## E. The predicted mass over many targets', '',
               f"{len(pe)} random three-bump targets (seed {SEED_E}): the fraction predicted not to recover at ε (ε ≥ ε_c), and the fractions "
               "with ε ≥ ε_entry (an entry's near tie alone would fail them) and with ε ≥ ε_creep (the creep alone) — not exclusive; the "
               "local slope of each in log ε. These are the predictor's fractions; the loop's own were not measured:", '',
               "| ε | not recovered | by a near tie | by the creep | slope: all | near tie | creep |", "|---|---|---|---|---|---|---|"]
        G_ = np.array(EPS_GRID_E)
        fr = {k: np.array([np.mean(v <= e) for e in G_]) for k, v in (('all', eps_c), ('tie', ee), ('creep', ec))}

        def loc(f, k):
            return np.log(f[k - 1] / f[k]) / np.log(G_[k - 1] / G_[k]) if k > 0 and f[k] > 0 and f[k - 1] > 0 else float('nan')
        for k, e in enumerate(G_):
            L_.append(f"| {e:.3g} | {fr['all'][k]:.4f} | {fr['tie'][k]:.4f} | {fr['creep'][k]:.4f} | {fmt(loc(fr['all'], k))} | "
                      f"{fmt(loc(fr['tie'], k))} | {fmt(loc(fr['creep'], k))} |")
        fits = []
        for k, name in (('all', 'all'), ('tie', 'a near tie'), ('creep', 'the creep')):
            m = (fr[k] * len(pe) >= 10) & (fr[k] <= 0.5)
            sl = np.polyfit(np.log(G_[m]), np.log(fr[k][m]), 1)[0] if m.sum() >= 2 else float('nan')
            fits.append(f"{name} {sl:.3f} (ε = {G_[m].max():.3g} … {G_[m].min():.3g})")
        L_ += ['', "Least-squares slopes in log ε where between 10 targets and half of them fail: " + '; '.join(fits) + '.', '']
    if F is not None and Ee is not None:
        pe = {p['n']: p for p in Ee['preds']}
        rows = F['rows']
        hi = [r for r in rows if r['eps'] > min(pe[r['n']]['eps_entry'], pe[r['n']]['eps_creep'])]
        lo = [r for r in rows if r not in hi]
        bind = lambda n: 'entry' if pe[n]['eps_entry'] <= pe[n]['eps_creep'] else 'creep'
        L_ += ['## F. The predicted tail against the loop', '',
               f"Every target of E with {F_RANGE[0]:g} ≤ ε_c ≤ {F_RANGE[1]:g} ({len(F['sel'])} targets; {F['below']} more below "
               f"{F_RANGE[0]:g}, not run), the loop at {F_FAC:g}ε_c (predicted not to recover) and at ε_c/{F_FAC:g} (predicted to recover):", '',
               "| | runs | as predicted | binding: entry | binding: creep |", "|---|---|---|---|---|",
               f"| {F_FAC:g}ε_c | {len(hi)} | {sum(1 for r in hi if not r['ok'])} | {sum(1 for r in hi if bind(r['n']) == 'entry')} | "
               f"{sum(1 for r in hi if bind(r['n']) == 'creep')} |",
               f"| ε_c/{F_FAC:g} | {len(lo)} | {sum(1 for r in lo if r['ok'])} | {sum(1 for r in lo if bind(r['n']) == 'entry')} | "
               f"{sum(1 for r in lo if bind(r['n']) == 'creep')} |", '']
        bad = [r for r in hi if r['ok']] + [r for r in lo if not r['ok']]
        L_ += ["Not as predicted (target, ε/ε_c, the binding threshold, formed against the path's): " +
               ('; '.join(f"({r['n']}, {r['eps'] / min(pe[r['n']]['eps_entry'], pe[r['n']]['eps_creep']):.3g}, {bind(r['n'])}, "
                          f"{r['order']} vs {r['path']})" for r in bad) or 'none') + '.', '']
    if G is not None:
        rows = G['rows']
        r_ = {r['n']: (np.sqrt(r['lo'] * r['hi']) / r['ec'] if np.isfinite(r['hi']) and r['lo'] > 0 else (np.inf if not np.isfinite(r['hi'])
                                                                                                          else 0.0)) for r in rows}
        L_ += ['## G. The loop\'s own threshold against the prediction', '',
               f"D's {len(rows)} targets: the largest mass at which the loop still forms exactly the path's kernels, bracketed from ε_c by "
               f"factors {G_FAC:g} and bisected {G_STEPS} times in log ε (the bracket ×{G_FAC ** (0.5 ** G_STEPS):.3g}); ε_loop its geometric "
               "mid-point.", '',
               '| binding | targets | ε_loop/ε_c: median [quartiles] | range | below 1 (the loop fails below ε_c) | within ×1.5 | within ×2 |',
               '|---|---|---|---|---|---|---|']
        for lab in ('entry', 'creep', 'all'):
            v = np.array([r_[n] for n in r_ if lab == 'all' or G['bind'][n] == lab])
            if not v.size:
                continue
            L_.append(f"| {lab} | {v.size} | {np.median(v):.3g} [{np.percentile(v, 25):.3g}, {np.percentile(v, 75):.3g}] | {v.min():.3g}–{v.max():.3g} | "
                      f"{int(np.sum(v < 1))} | {int(np.sum((v >= 1 / 1.5) & (v <= 1.5)))} | {int(np.sum((v >= 0.5) & (v <= 2)))} |")
        # the entry-bound ratio against the predicted mass (ε_c ∝ the binding margin squared)
        bins = ((0.0, 2e-6), (2e-6, 1e-4), (1e-4, np.inf))
        parts_ = []
        for lo_, hi_ in bins:
            v = np.array([r_[r['n']] for r in rows if G['bind'][r['n']] == 'entry' and lo_ < r['ec'] <= hi_])
            if v.size:
                parts_.append(f"ε_c in ({lo_:g}, {hi_:g}]: {v.size} targets, median {np.median(v):.3g}, range {v.min():.3g}–{v.max():.3g}")
        L_ += ['', 'Entry-bound targets by the predicted mass (∝ the binding margin squared): ' + '; '.join(parts_) + '.']
        L_ += ['', 'Per target (ε_c, the binding threshold, ε_loop/ε_c): ' + '; '.join(
            f"{r['n']}: {r['ec']:.3g} ({G['bind'][r['n']]}) {r_[r['n']]:.3g}" for r in rows) + '.', '']
    if G2 is not None:
        L_ += ["F's targets not as predicted, their loop thresholds bracketed and bisected as in G:", '',
               '| target (E) | binding | ε_c | ε_loop | ε_loop/ε_c |', '|---|---|---|---|---|']
        for r in G2['rows']:
            el = np.sqrt(r['lo'] * r['hi']) if np.isfinite(r['hi']) and r['lo'] > 0 else (np.inf if not np.isfinite(r['hi']) else 0.0)
            L_.append(f"| {r['n']} | {r['bind']} | {r['ec']:.3g} | {el:.3g} | {el / r['ec']:.3g} |")
        L_.append('')
    if H is not None:
        L_ += ['## H. Near-tie controls', '',
               f"(1) The single bump (A = {A_H:g}) at δ from the midpoint of the kernels at 0.0 and 0.1: the runner-up's formation threshold, "
               f"x_θ = μ*/ε^{{1/2}} at the loop's bisected threshold ({H_STEPS} bisections), free and with centres held fixed until the "
               "entrant's arrival τ_e (the inner problem's x_θ and the Gibbs value in A):", '',
               "| δ | μ* | free | the runner-up's centre fixed before τ_e | both centres fixed before τ_e |", '|---|---|---|---|---|']
        for off in OFFS_H:
            rr = {r['frz']: r for r in H['single'] if r['off'] == off}
            L_.append(f"| {off:g} | {rr['free']['mu']:.4f} | {rr['free']['x_theta']:.3f} | {rr['runner-up']['x_theta']:.3f} | {rr['both']['x_theta']:.3f} |")
        if HI is not None:
            xi_ = {r['A']: r['x_theta'] for r in HI['rows']}
            L_ += ['', "The inner problem's own threshold at the exact tie (the limit ε → 0 at fixed x), bisected in x: " +
                   ', '.join(f"A = {r['A']:g}: {r['x_theta']:.3f}" for r in HI['rows']) +
                   f" (the Gibbs value ln(|A|/θ − 1)/κ = {HI['rows'][0]['gibbs']:.3f}). The free loop's shift from it against the partners' "
                   "approach-creep difference at the actual geometry, (a²b/(π²c))((g′_r/g_r)² − (g′_e/g_e)²)/x_θ (g = ⟨fφ⟩, g′ = ⟨f∂_cφ⟩ "
                   "before the entry):", '', '| δ | μ* | x_θ (free) − inner | creep difference / x_θ | ratio |', '|---|---|---|---|---|']
            if A_H in xi_:
                Lh, Kh, c0h = setup()
                for off in OFFS_H:
                    fr = [r for r in H['single'] if r['off'] == off and r['frz'] == 'free'][0]
                    _, _, _, _, fh, gh, we_, nr_, _ = near_tie_case(A_H, off)
                    gp = Kh.dphi(c0h).T @ fh / Kh.xs.size
                    dd = (gp[nr_] / gh[nr_]) ** 2 - (gp[we_] / gh[we_]) ** 2
                    pr = Lh.A ** 2 / (np.pi ** 2 * Lh.c) * dd / fr["x_theta"]
                    L_.append(f"| {off:g} | {fr['mu']:.4f} | {fr['x_theta'] - xi_[A_H]:.3f} | {pr:.3f} | {(fr['x_theta'] - xi_[A_H]) / pr:.3f} |")
        if HP is not None:
            L_ += ['', "The predictor's own threshold for the runner-up (`recovery_threshold`): x = μ*/ε_entry^{1/2} = max(Δ + ξ_θ, x_θκ_ij) — "
                   "the frozen runner-up's pull plus the formation level's inner coordinate, and the near-tie estimate: " +
                   '; '.join(f"δ = {r['off']:g}: {r['x']:.3f} (kernel {r['j']}{'' if r['j'] == r['runner'] else ', not the runner-up'}; "
                             f"Δ + ξ_θ = {r['pull_x']:.3f}, x_θκ_ij = {r['xth']:.3f}·{r['kij']:.3f} = {r['xth'] * r['kij']:.3f})" for r in HP['rows']) + '.']
        L_ += ['', "(2) Exact ties made by unequal masses (the farther kernel's mass b set so that both reach their walls together):", '',
               "| bump | nearer kernel | farther kernel (b) | ε | the farther kernel's share | formed | the path |", '|---|---|---|---|---|---|---|']
        for r in H['masses']:
            L_.append(f"| {r['cstar']:g} | {r['jn']} | {r['jf']} ({r['b']:.4f}) | {r['eps']:g} | {r['share']:.4f} | {r['formed']} | {r['path']} |")
        L_ += ['', f"(3) D's target {MIG_N} with its middle bump moved to c₂: the first entry's runner-up and its margin μ*, the predicted ε_c "
                   "(entry, creep), the loop's threshold (bracketed and bisected as in G), and the kernels it forms just above it (an extra "
                   "kernel's final centre):", '',
               '| c₂ | runner-up (μ*) | ε_c: entry, creep | ε_loop | ε_loop/ε_c | ε_loop/μ*² | ε_loop/μ*³ | formed just above (extra: final centre) |',
               '|---|---|---|---|---|---|---|---|']
        for r in H['migration']:
            el = np.sqrt(r['lo'] * r['hi'])
            ec = min(r['eps_entry'], r['eps_creep'])
            L_.append(f"| {r['c2']:g} | {r['runner']} ({r['mu']:.5f}) | {r['eps_entry']:.3g}, {r['eps_creep']:.3g} | {el:.3g} | {el / ec:.3g} | "
                      f"{el / r['mu'] ** 2:.3g} | {el / r['mu'] ** 3:.3g} | {r['order']} ({', '.join(f'{j}: {c:.3f}' for j, c in r['extra'].items()) or '—'}) |")
        mm = [(r['mu'], np.sqrt(r['lo'] * r['hi'])) for r in H['migration'] if r['extra'] and r['runner'] in r['extra']]
        if len(mm) >= 2:
            L_ += ['', f"Where the extra kernel is the first entry's runner-up, the loop's threshold against μ*: slope "
                       f"{np.polyfit(np.log([m for m, e in mm]), np.log([e for m, e in mm]), 1)[0]:.2f} in log–log.", '']
        L_ += ["(4) E's target 117 — kernel 14 between the third entrant (13) and its bump: its largest coefficient over the run (the "
               "formation level θ = 10⁻²max|A*|) and its closest approach to the entrant:", '',
               '| ε | x = μ/ε^{1/2} | max\\|w₁₄\\| | θ | min \\|c₁₃ − c₁₄\\| | formed | the path | the Gibbs split at x |', '|---|---|---|---|---|---|---|---|']
        kap_ = np.pi * np.sqrt(2 * setup()[0].c)
        for r in H['near_side']:
            L_.append(f"| {r['eps']:g} | {r['x']:.3f} | {r['w14']:.3g} | {r['theta']:.3g} | {r['gap']:.4f} | {r['order']} | {r['path']} | "
                      f"{1 / (1 + np.exp(kap_ * r['x'])):.3f} |")
        L_.append('')
    L_.append('Numbers are produced by `experiments/exp_transit_mass.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_mass_raw.npz')
    md = os.path.join(OUT, 'transit_mass.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    keys = ('A', 'B', 'C', 'D', 'E', 'A2', 'F', 'G', 'H', 'G2', 'HI', 'HP')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('B', part_B), ('A', part_A), ('A2', lambda: part_A2(res['A'])), ('C', part_C), ('D', part_D), ('E', part_E),
                  ('F', lambda: part_F(res['E'])), ('G', lambda: part_G(res['D'])), ('H', part_H),
                  ('G2', lambda: part_G2(res['E'], res['F'])), ('HI', part_HI),
                  ('HP', part_HP)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
