"""
Proposition D‴.49: the rich transit in the closed loop — from small masses the adaptive controller forms kernels one at a time
====================================================================================================================================

The scalar plant ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), e = x − x_m, x_m = x₀ + a sin(ωt); f̂ = Σ_iw_iφ(·; c_i), Gaussian kernels of width 0.12
started at α = 0 with the masses ε on a grid of 21 centres in [−1, 1]; the paradigm's free law with moving centres driven by the error,
dα_i = Gβ_ia′_ieφ_i(x), dβ_i = γβ_ia_ieφ_i(x), ċ_i = Gβ_ia_ie∂_cφ_i(x) (Proposition D‴.48(iii); `rich_transit`). f = three bumps of the kernels'
width at −0.52, 0.07, 0.61 with amplitudes 0.9, −0.6, 0.75 (max|f′| = 4.55). The real loop is integrated in C (`closed_loop`, the kernels
kept on their curves); the averaged loop along the frozen periodic orbit (`rich_transit.FrozenOrbit`, exact for k_e > Lip(f − f̂)) and
its leading order in 1/k_e, the reference's own orbit seen through the error filter (`OrbitBumps`). Entry times (a kernel formed when
|w| > 10⁻²max|A*|) are interpolated in log|w| between output points; τ = t/Λ(ε).

A. The base point (k_e = 10, ω = 0.2, x_m = sin ωt): the exact averaged loop and its leading order for ε = 10⁻⁴ … 10⁻⁷, the real loop at
   ε = 10⁻⁴, 10⁻⁵, 10⁻⁶ to the same τ = TAU_MAX — entries, kernels, coefficients, the inactive kernels' displacement, the tracking error.
B. The controller's knobs at ε = 10⁻⁵: k_e = 2 … 40, ω = 0.05 and 1, a = 0.6 and 0.8, a quasi-periodic reference — the real loop against
   the exact averaged loop (where the frozen orbit exists) to the same τ (1.3 times the averaged loop's last entry): the entries, kernels,
   coefficients and inactive displacements of both, the kernels pushed out of [−1.2, 1.2], the state's range, the formed centres' largest
   swing within a period during the transit (N_WIN windows) and at the end, the averaged active block's rates at the end (the real parts
   of the eigenvalues of its linearisation in t, against ω and k_e).
B2. Low gain: k_e = 2 against the mass (ε = 10⁻⁴ … 10⁻⁷, to τ = TAU_B2); f scaled by 2 at k_e = 4, 6, 8 and by ½ at k_e = 2, 3 (ε = 10⁻⁵);
   one ejection traced (k_e = 2, ε = 10⁻⁴, a dense record): each expelled kernel's formation, the period in which it left the range, the
   state's range in that period and the other formed kernels' swing.
C. A feature near the reference's range (x_m = 0.45 + 0.5 sin ωt, the range [−0.05, 0.95]): the first bump moved to 0.11, 0.19, 0.35 below
   it (0.9σ, 1.6σ, 2.9σ), the exact averaged loop at ε = 10⁻⁶ to τ = TAU_C, with the bump at 0.07 as given and flipped to +0.6: whether and
   when the bump is formed, against the averaged drive's size at its nearest kernel before any entry, the orbit's error at CHECK_C; and
   Part C's original bump at −0.52.
D. The positions' metric: the exact averaged loop under the Wasserstein–Fisher–Rao metric and with fixed centres, ε = 10⁻⁶ and 10⁻⁸, to
   τ = TAU_D; for fixed centres the averaged drive's zero on the formed kernels and on the whole grid (the orbit's error of each).
E. Measurement noise (ε = 10⁻⁴; 16 tones in 10–30 rad/s of rms s in the measured x, seen by the error, the control and the regressor):
   s = 0.005, 0.01, 0.02 with the free law, and with Proposition D‴.46's leak −λη toward the start (λΛ(ε) = LEAK_E in the time τ — the
   leak's weight on the scaled dual u = η/Λ, a LASSO weight λΛρ on the averaged fit); E2: weaker leaks; E3: the desired regressor
   (D‴.31) at s = 0, 0.01, 0.02 and stronger leaks at s = 0.02. Each row: the measured and the true error, the in-range fit's rms error
   and the largest coefficient left below the formation threshold.

Output: results/chaos/closed_transit.md and closed_transit_raw.npz (`--from-raw` regenerates the md; `--only A,B,B2,C,D,E,E2,E3`; progress
to stdout).
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
from src.core import mirror_free as MF
from src.core import rich_transit as RT
from src.core import closed_loop as CL
import exp_cohesion_leak as E

OUT = os.environ.get('CT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CT_WORKERS', '3'))
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                      # c = 0.025
KE, OM, SIG = 10.0, 0.2, 0.12
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)
M_ORBIT, TAU_MAX, TH = 512, 600.0, 1e-2                              # orbit samples; horizon in τ; formed: |w| > TH·max|A*|
EPS_A_AVG, EPS_A_REAL = (1e-4, 1e-5, 1e-6, 1e-7), (1e-4, 1e-5, 1e-6)
EPS_B = 1e-5
KE_B = (2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 40.0)
REF_B = (('ω = 0.05', dict(om=0.05)), ('ω = 1', dict(om=1.0)), ('a = 0.6', dict(a=0.6)), ('a = 0.8', dict(a=0.8)),
         ('quasi-periodic 0.6 sin ωt + 0.4 sin √2ωt', dict(a=0.6, a2=0.4)))
EPS_B2, KE_B2, TAU_B2 = (1e-4, 1e-5, 1e-6, 1e-7), 2.0, 150.0          # B2: low gain against the mass
AMP_B2 = ((2.0, (4.0, 6.0, 8.0)), (0.5, (2.0, 3.0)))                  # B2: f scaled, at these gains (ε = EPS_B)
TRACE_B2 = (2.0, 1e-4, 100.0)                                         # B2: the traced ejection (k_e, ε, τ)
N_WIN = 24                                                            # windows of one period for the transit-time swing
REF_HALF, D_C, EPS_C = (0.45, 0.5), (0.11, 0.19, 0.35), 1e-6
TAU_C, CHECK_C, FLIP_C = 20000.0, (1000.0, 5000.0, 10000.0, 20000.0), 0.6
EPS_D, TAU_D = (1e-6, 1e-8), 30000.0
EPS_E, SIG_E, LEAK_E, TAU_E = 1e-4, (0.005, 0.01, 0.02), 3e-3, 600.0
LEAKS_E2 = (3e-4, 1e-3)                                               # E2: weaker leaks, without noise and at s = 0.01
E3 = ((0.0, 0.0, True), (0.01, 0.0, True), (0.02, 0.0, True), (0.02, 5e-3, False), (0.02, 1e-2, False))   # E3: (s, λΛ, regx)
OUTSIDE = 1.2                                                         # a kernel pushed out: |c| > OUTSIDE


def law():
    return MF.FreeLaw(**LAW)


def phi(x, c):
    return np.exp(-(x - c) ** 2 / (2 * SIG ** 2))


def make_f(cstar=CSTAR, astar=ASTAR):
    cstar, astar = np.asarray(cstar, float), np.asarray(astar, float)
    return lambda x: np.exp(-(np.asarray(x, float)[..., None] - cstar) ** 2 / (2 * SIG ** 2)) @ astar


def th_abs(astar=ASTAR):
    return TH * np.abs(astar).max()


def crossings(W, taus, th):
    """(τ, i) of each kernel's first |w| > th, interpolated in log|w| between output points."""
    out = []
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            k = k[0]
            if k == 0:
                out.append((float(taus[0]), i))
                continue
            a, b = np.log(max(abs(W[k - 1, i]), 1e-300)), np.log(abs(W[k, i]))
            out.append((float(taus[k - 1] + (taus[k] - taus[k - 1]) * (np.log(th) - a) / (b - a)), i))
    return sorted(out)


def averaged(L, eps, ke=KE, om=OM, x0=0.0, a=1.0, cstar=CSTAR, astar=ASTAR, metric='paradigm', exact=True, tau_max=TAU_MAX, n_t=2401,
             checkpoints=()):
    """The averaged loop from α = 0, masses ε, to τ_max: exact (the frozen orbit) or its leading order (the reference's orbit through the
    filter). Returns None for the exact loop when the frozen orbit does not converge (k_e too small against Lip(f − f̂)). checkpoints: the
    τ's at which the state (c, w) is also returned ('checks')."""
    d = C0.size
    Lam = RT.transit_scale(L, eps)
    taus = np.unique(np.r_[np.linspace(0, tau_max, n_t), [t for t in checkpoints if t <= tau_max]])
    f = make_f(cstar, astar)
    xm = lambda t: x0 + a * np.sin(om * t)
    if exact:
        O = RT.FrozenOrbit(xm, 2 * np.pi / om, M_ORBIT, SIG, ke, f)
        rhs = RT.frozen_orbit_rhs(L, O, np.full(d, eps), metric)
    else:
        K = RT.OrbitBumps(xm, 2 * np.pi / om, M_ORBIT, SIG, ke)
        rhs = RT.critical_rhs(L, K, f(K.xs), 'paradigm' if metric == 'fixed' else metric, b0=np.full(d, eps))
        if metric == 'fixed':
            base = rhs
            rhs = lambda t, z: np.r_[base(t, z)[:d], np.zeros(d)]
    try:
        sol = solve_ivp(rhs, (0, Lam * tau_max), np.r_[np.zeros(d), C0], method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    except RuntimeError:
        return None
    al, cc = sol.y[:d].T, sol.y[d:].T
    W = (eps + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    th = th_abs(astar)
    act = np.abs(W[-1]) > th
    ent = crossings(W, taus, th)
    cl = np.sort(cc[-1][act])
    if exact:
        e = O.orbit(cc[-1], W[-1])
        rms = float(np.sqrt(np.mean(e ** 2)))
    else:
        rms = float(np.sqrt(np.mean((f(K.xs) - K.phi(cc[-1]) @ W[-1]) ** 2)))
    checks = [(float(t), cc[int(np.argmin(np.abs(taus - t)))].tolist(), W[int(np.argmin(np.abs(taus - t)))].tolist()) for t in checkpoints
              if t <= tau_max]
    return dict(eps=eps, Lam=float(Lam), tau_max=float(tau_max), entries=ent, n_act=int(act.sum()), kernels=np.flatnonzero(act).tolist(),
                centres=cc[-1][act].tolist(), w=W[-1][act].tolist(), clusters=int(np.sum(np.diff(cl) > 1e-3) + 1) if cl.size else 0,
                rms=rms, disp_in=float(np.abs(cc[-1] - C0)[~act].max()) if (~act).any() else 0.0,
                sub=float(np.abs(W[-1][~act]).max()) if (~act).any() else 0.0, al_end=al[-1].tolist(), c_end=cc[-1].tolist(),
                w_end=W[-1].tolist(), checks=checks, t_last=float(Lam * max(t for t, i in ent)) if ent else np.nan, done=sol.status == 0)


def active_rates(L, eps, av, ke=KE, om=OM, x0=0.0, a=1.0, cstar=CSTAR, astar=ASTAR, metric='paradigm'):
    """The exact averaged loop's linearisation at its end, restricted to the formed kernels' (α, c): the real parts of its eigenvalues
    in the time t (central differences of `frozen_orbit_rhs`)."""
    d = C0.size
    O = RT.FrozenOrbit(lambda t: x0 + a * np.sin(om * t), 2 * np.pi / om, M_ORBIT, SIG, ke, make_f(cstar, astar))
    rhs = RT.frozen_orbit_rhs(L, O, np.full(d, eps), metric)
    z = np.r_[av['al_end'], av['c_end']]
    idx = np.r_[av['kernels'], d + np.asarray(av['kernels'], int)].astype(int)
    if idx.size == 0:
        return []
    J = np.zeros((idx.size, idx.size))
    for k, i in enumerate(idx):
        h = 1e-6 * max(1.0, abs(z[i]))
        zp, zm = z.copy(), z.copy()
        zp[i] += h
        zm[i] -= h
        J[:, k] = (rhs(0.0, zp)[idx] - rhs(0.0, zm)[idx]) / (2 * h)
    return sorted(float(v) for v in np.linalg.eigvals(J).real)


def real(L, eps, T, ke=KE, om=OM, x0=0.0, a=1.0, a2=0.0, cstar=CSTAR, astar=ASTAR, noise=None, lam=0.0, regx=False, n_end=6, n_out=401,
         n_win=N_WIN):
    """The real loop (C) from x(0) = x_m(0), α = 0, masses ε, to the time T; a dense record of the last n_end periods for the end's swing
    and n_win windows of one period over [0.05T, T] for the transit-time swing (the formed centres' range within a period). Also the true
    error x − x_m over the last period, the state's range, the in-range fit's rms error (over [x₀ − a − a₂, x₀ + a + a₂]) and the largest
    coefficient below the formation threshold."""
    per = 2 * np.pi / om
    wins = [np.linspace(t0, t0 + per, 25) for t0 in np.linspace(0.05 * T, T - per, n_win)] if n_win else []
    ts = np.unique(np.concatenate([np.linspace(0, T, n_out), np.linspace(max(T - n_end * per, 0), T, 60 * n_end)] + wins))
    ref = (x0, a, om, a2, np.sqrt(2) * om)
    th = th_abs(astar)
    t0 = time.time()
    r = CL.simulate(L, eps, C0, SIG, (cstar, astar, SIG), ke, ref, ts, noise=noise, lam=lam, th=th, regx=regx)
    Lam = RT.transit_scale(L, eps)
    act = np.abs(r['w'][-1]) > th
    ent = sorted((float(t / Lam), i) for i, t in enumerate(r['entries']) if np.isfinite(t))
    last = ts >= T - per
    swing = float(np.ptp(r['c'][last][:, act], axis=0).max()) if act.any() else 0.0
    sw = []
    for wv in wins:
        idx = np.searchsorted(ts, wv).clip(0, ts.size - 1)
        k = np.abs(r['w'][idx[-1]]) > th
        sw.append((float(wv[0] / Lam), float(np.ptp(r['c'][idx][:, k], axis=0).max()) if k.any() else 0.0))
    xs = np.linspace(x0 - a - a2, x0 + a + a2, 401)
    fit = make_f(cstar, astar)(xs) - phi(xs[:, None], r['c'][-1][None, :]) @ r['w'][-1]
    return dict(eps=eps, T=float(T), entries=ent, n_act=int(act.sum()), kernels=np.flatnonzero(act).tolist(),
                centres=r['c'][-1][act].tolist(), w=r['w'][-1][act].tolist(), swing=swing, swing_transit=max((v for t, v in sw), default=0.0),
                swings=sw, outside=int(np.sum(np.abs(r['c'][-1][act]) > OUTSIDE)), wmax=float(np.abs(r['w'][-1]).max()),
                disp_in=float(np.abs(r['c'][-1] - C0)[~act].max()) if (~act).any() else 0.0,
                sub=float(np.abs(r['w'][-1][~act]).max()) if (~act).any() else 0.0, fit_rms=float(np.sqrt(np.mean(fit ** 2))),
                e_first=r['e_first'], e_last=r['e_last'], e_true=r['e_true_last'], e_true_rms=r['e_true_rms_last'],
                x_range=list(r['x_range']), secs=time.time() - t0, steps=r['steps'])


def _pool():
    from multiprocessing import Pool
    return Pool(N_WORKERS)


def a_avg(args):
    eps, exact = args
    warnings.filterwarnings('ignore')
    return averaged(law(), eps, exact=exact)


def a_real(eps):
    warnings.filterwarnings('ignore')
    L = law()
    return real(L, eps, RT.transit_scale(L, eps) * TAU_MAX)


def part_A():
    with _pool() as pool:
        ex = pool.map(a_avg, [(e, True) for e in EPS_A_AVG])
        lo = pool.map(a_avg, [(e, False) for e in EPS_A_AVG])
        rl = pool.map(a_real, list(EPS_A_REAL))
    return dict(exact=ex, lead=lo, real=rl)


def b_job(args):
    """One knob (or B2's case): the exact averaged loop to TAU_MAX for its last entry, then both loops to the same τ_h = 1.3 times it (or
    390·k_e/10 where there is no averaged loop, eps the mass, amp the factor on f)."""
    label, kw = args
    warnings.filterwarnings('ignore')
    L = law()
    ke, om, a, a2 = kw.get('ke', KE), kw.get('om', OM), kw.get('a', 1.0), kw.get('a2', 0.0)
    eps, astar = kw.get('eps', EPS_B), ASTAR * kw.get('amp', 1.0)
    Lam = RT.transit_scale(L, eps)
    av0 = averaged(L, eps, ke=ke, om=om, a=a, astar=astar) if (a2 == 0.0 and 'tau' not in kw) else None
    if 'tau' in kw:
        tau_h = kw['tau']
    elif av0 is not None and np.isfinite(av0['t_last']):
        tau_h = 1.3 * av0['t_last'] / Lam
    else:
        tau_h = 390.0 * ke / KE
    av = (averaged(L, eps, ke=ke, om=om, a=a, astar=astar, tau_max=tau_h, n_t=int(40 * tau_h) + 1)
          if (a2 == 0.0 and ('tau' in kw or av0 is not None)) else None)
    rates = active_rates(L, eps, av, ke=ke, om=om, a=a, astar=astar) if av is not None else None
    rl = real(L, eps, Lam * tau_h, ke=ke, om=om, a=a, a2=a2, astar=astar)
    print(f"  B {label}: averaged {av['n_act'] if av else '—'}, real {rl['n_act']} ({rl['secs']:.0f}s)", flush=True)
    return dict(label=label, kw=kw, tau=float(tau_h), avg=av, real=rl, rates=rates)


def part_B():
    jobs = [(f'k_e = {k:g}', dict(ke=k)) for k in KE_B] + list(REF_B)
    with _pool() as pool:
        return pool.map(b_job, jobs)


def trace_job(args):
    """B2's traced ejection: the real loop with a record of 30 points per period; for each kernel that formed and then left
    [−OUTSIDE, OUTSIDE], its formation, the first output at which it is out, the period before it — its centre and coefficient at the
    period's start, the state's range over the period and the other formed kernels' largest swing in it."""
    ke, eps, tau = args
    warnings.filterwarnings('ignore')
    L = law()
    Lam = RT.transit_scale(L, eps)
    per = 2 * np.pi / OM
    ts = np.arange(0.0, Lam * tau + per / 60, per / 30)
    th = th_abs()
    r = CL.simulate(L, eps, C0, SIG, (CSTAR, ASTAR, SIG), ke, (0.0, 1.0, OM, 0.0, np.sqrt(2) * OM), ts, th=th)
    W, Cc, X = r['w'], r['c'], r['x']
    ev = []
    for j in range(C0.size):
        fm = np.flatnonzero(np.abs(W[:, j]) > th)
        if not fm.size:
            continue
        out = np.flatnonzero((np.abs(Cc[:, j]) > OUTSIDE) & (np.arange(ts.size) > fm[0]))
        if not out.size:
            continue
        ix = int(out[0])
        i0 = max(ix - 30, 0)
        oth = [k for k in range(C0.size) if k != j and np.abs(W[i0, k]) > th]
        ev.append(dict(kernel=j, tau_formed=float(ts[fm[0]] / Lam), c_formed=float(Cc[fm[0], j]), tau_out=float(ts[ix] / Lam),
                       c_from=float(Cc[i0, j]), c_to=float(Cc[ix, j]), w=float(W[i0, j]), x_lo=float(X[i0:ix + 1].min()),
                       x_hi=float(X[i0:ix + 1].max()), swing_others=float(np.ptp(Cc[i0:ix + 1][:, oth], axis=0).max()) if oth else 0.0))
    return dict(ke=ke, eps=eps, tau=tau, x_range=list(r['x_range']), n_act=int((np.abs(W[-1]) > th).sum()), events=ev)


def part_B2():
    jobs = ([(f'k_e = {KE_B2:g}, ε = {e:g}', dict(ke=KE_B2, eps=e, tau=TAU_B2)) for e in EPS_B2]
            + [(f'f × {m:g}, k_e = {k:g}', dict(ke=k, amp=m)) for m, ks in AMP_B2 for k in ks])
    with _pool() as pool:
        tr = pool.apply_async(trace_job, (TRACE_B2,))
        rows = pool.map(b_job, jobs)
        return dict(rows=rows, trace=tr.get())


def c_job(args):
    dist, flip = args
    warnings.filterwarnings('ignore')
    L = law()
    lo = REF_HALF[0] - REF_HALF[1]
    cstar, astar = CSTAR.copy(), ASTAR.copy()
    cstar[0] = lo - dist if dist is not None else CSTAR[0]
    if flip:
        astar[1] = FLIP_C
    av = averaged(L, EPS_C, x0=REF_HALF[0], a=REF_HALF[1], cstar=cstar, astar=astar, tau_max=TAU_C, n_t=int(2 * TAU_C) + 1,
                  checkpoints=CHECK_C)
    O = RT.FrozenOrbit(lambda t: REF_HALF[0] + REF_HALF[1] * np.sin(OM * t), 2 * np.pi / OM, M_ORBIT, SIG, KE, make_f(cstar, astar))
    g, h, e = O.drive(C0, np.zeros(C0.size))
    j = int(np.argmin(np.abs(C0 - cstar[0])))
    first = [t for t, i in av['entries'] if abs(C0[i] - cstar[0]) <= 0.15]
    hist = [(t, float(np.sqrt(np.mean(O.orbit(np.array(c), np.array(w)) ** 2)))) for t, c, w in av['checks']]
    return dict(dist=dist, flip=bool(flip), bump=float(cstar[0]), avg=av, drive_near=float(abs(g[j])), drive_max=float(np.abs(g).max()),
                kernel=j, t_first=float(min(first)) if first else np.nan, err_hist=hist,
                fitted=bool(any(abs(c - cstar[0]) < 0.01 for c in av['centres'])))


def part_C():
    with _pool() as pool:
        return pool.map(c_job, [(d, f) for f in (False, True) for d in list(D_C) + [None]])


def drive_zero(L, S, w0):
    """Fixed centres: the exact averaged drive's zero on the grid kernels S (the others at 0), from w0 — the closed loop's least-squares
    fit on S; returns (w, the orbit's error rms, the drive's largest remaining component)."""
    from scipy.optimize import least_squares
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), 2 * np.pi / OM, M_ORBIT, SIG, KE, make_f())
    S = np.asarray(S, int)
    if S.size == 0:
        return np.zeros(C0.size).tolist(), float(np.sqrt(np.mean(O.orbit(C0, np.zeros(C0.size)) ** 2))), 0.0

    def g(v):
        w = np.zeros(C0.size)
        w[S] = v
        return O.drive(C0, w)[0][S]
    sol = least_squares(g, np.asarray(w0, float)[S], xtol=1e-15, ftol=1e-15, gtol=1e-15)
    w = np.zeros(C0.size)
    w[S] = sol.x
    return w.tolist(), float(np.sqrt(np.mean(O.orbit(C0, w) ** 2))), float(np.abs(sol.fun).max())


def d_job(args):
    metric, eps = args
    warnings.filterwarnings('ignore')
    L = law()
    av = averaged(L, eps, metric=metric, tau_max=TAU_D, n_t=int(TAU_D / 2) + 1)
    if metric == 'fixed':
        av['zero_formed'] = drive_zero(L, av['kernels'], av['w_end'])
        av['zero_grid'] = drive_zero(L, np.arange(C0.size), av['w_end'])
    return av


def part_D():
    with _pool() as pool:
        rows = pool.map(d_job, [(m, e) for m in ('wfr', 'fixed') for e in EPS_D])
    return [dict(metric=m, **r) for (m, e), r in zip([(m, e) for m in ('wfr', 'fixed') for e in EPS_D], rows)]


def e_job(args):
    s, lam = args[:2]
    regx = args[2] if len(args) > 2 else False
    warnings.filterwarnings('ignore')
    L = law()
    Lam = RT.transit_scale(L, EPS_E)
    r = real(L, EPS_E, Lam * TAU_E, noise=CL.tones(s, 10.0, 30.0, seed=7) if s > 0 else None, lam=lam / Lam, regx=regx)
    print(f"  E s {s} λ {lam} regx {regx}: {r['n_act']} kernels, {r['outside']} outside, max|w| {r['wmax']:.3g} ({r['secs']:.0f}s)",
          flush=True)
    return dict(s=s, lam=lam, regx=bool(regx), **r)


def part_E():
    jobs = [(0.0, 0.0), (0.0, LEAK_E)] + [(s, lam) for lam in (0.0, LEAK_E) for s in SIG_E]
    with _pool() as pool:
        return pool.map(e_job, jobs)


def part_E2():
    with _pool() as pool:
        return pool.map(e_job, [(s, lam) for lam in LEAKS_E2 for s in (0.0, 0.01)])


def part_E3():
    with _pool() as pool:
        return pool.map(e_job, list(E3))


def lip_f(cstar=CSTAR, astar=ASTAR):
    """max|f′| on a fine grid (the hypothesis k_e > Lip(r) of (i) at r = f, the start)."""
    x = np.linspace(-3, 3, 60001)
    return float(np.abs(np.gradient(make_f(cstar, astar)(x), x)).max())


def report(A=None, B=None, C=None, D=None, E_=None, E2=None, B2=None, E3=None):
    fmt = E.fmt
    L = law()
    ev = lambda ent: ', '.join(f"{i}: {t:.4g}" for t, i in ent)
    rnd = lambda v: [round(x, 4) for x in v]
    lip = lip_f()
    L_ = ["# Proposition D‴.49: the rich transit in the closed loop (`exp_closed_transit.py`)", "",
          f"ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), x_m = x₀ + a sin(ωt) (base point k_e = {KE:g}, ω = {OM:g}, x₀ = 0, a = 1); f three bumps of "
          f"width {SIG} at {CSTAR.tolist()} with amplitudes {ASTAR.tolist()} (max|f′| = {lip:.3g}); "
          f"21 kernels of width {SIG} started on a grid of [−1, 1] at α = 0 with the masses ε; the free law (m = {LAW['m']:g}, G = "
          f"{LAW['G']:g}, γ = {LAW['gamma']:g}, a_max = {LAW['a_max']:g}) with moving centres; τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{{1/2}}. "
          f"A kernel has formed when |w| > {TH:g}·max|A*|; entries listed as kernel: τ (kernel i started at −1 + 0.1i), interpolated. The "
          "real loop in C (`closed_loop`); the exact averaged loop along the frozen periodic orbit (`FrozenOrbit`; (i)'s contraction needs "
          f"k_e > Lip(f − f̂), at the start k_e > {lip:.3g} — below it the fixed-point iteration may still converge, outside (i)); its "
          "leading order in 1/k_e along the reference's orbit (`OrbitBumps`). The errors: max|e| the measured error x̃ − x_m (= x − x_m "
          "without noise), true the plant's x − x_m, over the last period.", ""]
    if A is not None:
        L_ += [f"## A. The base point: the exact averaged loop, its leading order and the real loop (all to τ = {TAU_MAX:g})", "",
               "| ε | loop | entries (kernel: τ) | formed | centres | coefficients | inactive: max \\|Δc\\| | the largest sub-threshold \\|w\\| | "
               "max \\|e\\|, last period |", "|---|---|---|---|---|---|---|---|---|"]
        reals = {r['eps']: r for r in A['real']}
        for ex, lo in zip(A['exact'], A['lead']):
            L_.append(f"| {ex['eps']:g} | averaged, exact | {ev(ex['entries'])} | {ex['n_act']} | {rnd(ex['centres'])} | {rnd(ex['w'])} | "
                      f"{fmt(ex['disp_in'])} | {fmt(ex.get('sub'))} | — |")
            L_.append(f"| {lo['eps']:g} | averaged, leading order | {ev(lo['entries'])} | {lo['n_act']} | {rnd(lo['centres'])} | "
                      f"{rnd(lo['w'])} | {fmt(lo['disp_in'])} | {fmt(lo.get('sub'))} | — |")
            if ex['eps'] in reals:
                r = reals[ex['eps']]
                L_.append(f"| {r['eps']:g} | real (τ = {r['T'] / RT.transit_scale(L, r['eps']):.4g}, {r['secs']:.0f} s) | {ev(r['entries'])} | "
                          f"{r['n_act']} | {rnd(r['centres'])} | {rnd(r['w'])} | {fmt(r['disp_in'])} | {fmt(r.get('sub'))} | {fmt(r['e_last'])} |")
        L_.append("")
    if B is not None:
        L_ += [f"## B. The controller's knobs (ε = {EPS_B:g})", "",
               "Both loops to the same τ (1.3 times the exact averaged loop's last entry; 390·k_e/10 without one); — where the frozen orbit "
               f"does not converge or the reference is not periodic; * below (i)'s hypothesis (k_e ≤ max|f′| = {lip:.3g}). (a) The two loops:", "",
               "| knob | τ | averaged: entries (kernel: τ) | averaged: centres, coefficients | real: entries (kernel: τ) | real: centres, "
               "coefficients | inactive max \\|Δc\\|: averaged, real |", "|---|---|---|---|---|---|---|"]
        for b in B:
            av, r = b['avg'], b['real']
            star = '*' if b['kw'].get('ke', KE) <= lip else ''
            L_.append(f"| {b['label']}{star} | {b['tau']:.4g} | {ev(av['entries']) if av else '—'} | "
                      f"{(rnd(av['centres']), rnd(av['w'])) if av else '—'} | {ev(r['entries'])} | {rnd(r['centres'])}, {rnd(r['w'])} | "
                      f"{fmt(av['disp_in']) if av else '—'}, {fmt(r['disp_in'])} |")
        L_ += ["", "(b) The real loop's state and active block: the kernels pushed out of "
               f"[−{OUTSIDE:g}, {OUTSIDE:g}], the state's range over the run, the formed centres' largest swing within one period during the "
               f"transit ({N_WIN} windows) and at the end, the averaged active block's rates at the end (the real parts of its "
               "linearisation's eigenvalues in t: the slowest and the fastest of those above 10⁻² of the largest, against ω and k_e, and the "
               "number below — the gauge of coinciding kernels), the errors over the last period:", "",
               "| knob | pushed out | the state's range | swing: transit, end | active rates: slowest, fastest (/ω; /k_e) | max \\|e\\|, true "
               "max \\|e\\| | the largest sub-threshold \\|w\\| |", "|---|---|---|---|---|---|---|"]
        for b in B:
            r, rt = b['real'], b.get('rates')
            om, ke = b['kw'].get('om', OM), b['kw'].get('ke', KE)
            rr = '—'
            if rt is not None and len(rt):
                rt = np.asarray(rt, float)
                nz_ = rt[np.abs(rt) > 1e-2 * np.abs(rt).max()]
                rr = (f"{nz_.max():.3g}, {nz_.min():.3g} ({abs(nz_.max()) / om:.3g}–{abs(nz_.min()) / om:.3g}ω; {abs(nz_.min()) / ke:.3g}k_e); "
                      f"{rt.size - nz_.size} null")
            L_.append(f"| {b['label']} | {r['outside']} | [{r['x_range'][0]:.4f}, {r['x_range'][1]:.4f}] | {fmt(r['swing_transit'])}, "
                      f"{fmt(r['swing'])} | {rr} | {fmt(r['e_last'])}, {fmt(r['e_true'])} | {fmt(r['sub'])} |")
        L_.append("")
    if B2 is not None:
        L_ += [f"## B2. Low gain", "",
               f"(a) k_e = {KE_B2:g} against the mass (to τ = {TAU_B2:g}) and f scaled (ε = {EPS_B:g}, both loops to 1.3 times the averaged "
               "loop's last entry): the kernels formed and pushed out, the state's range, the transit-time swing, the averaged loop's "
               "kernels:", "",
               "| case | τ | real: entries (kernel: τ) | real: kernels (pushed out) | final centres | the state's range | swing: transit | "
               "averaged: kernels (pushed out) |", "|---|---|---|---|---|---|---|---|"]
        for b in B2['rows']:
            av, r = b['avg'], b['real']
            avs = (f"{av['n_act']} ({int(np.sum(np.abs(np.array(av['centres'])) > OUTSIDE))})" if av else '—')
            L_.append(f"| {b['label']} | {b['tau']:.4g} | {ev(r['entries'])} | {r['n_act']} ({r['outside']}) | {rnd(r['centres'])} | "
                      f"[{r['x_range'][0]:.4f}, {r['x_range'][1]:.4f}] | {fmt(r['swing_transit'])} | {avs} |")
        tr = B2['trace']
        L_ += ["", f"(b) The traced ejections (k_e = {tr['ke']:g}, ε = {tr['eps']:g}, to τ = {tr['tau']:g}; the state's range over the run "
               f"[{tr['x_range'][0]:.4f}, {tr['x_range'][1]:.4f}], {tr['n_act']} kernels formed at the end): each expelled kernel's formation "
               f"(τ, centre), the first record beyond {OUTSIDE:g} (τ), one period earlier its centre and coefficient, the state's range and "
               "the other formed kernels' largest swing over that period:", "",
               "| kernel | formed at τ (centre) | out at τ | one period earlier: centre → | its \\|w\\| | the state's range in the period | "
               "the others' swing |", "|---|---|---|---|---|---|---|"]
        for e in tr['events']:
            L_.append(f"| {e['kernel']} | {e['tau_formed']:.4g} ({e['c_formed']:.3f}) | {e['tau_out']:.4g} | {e['c_from']:.3f} → {e['c_to']:.3f} | "
                      f"{abs(e['w']):.3g} | [{e['x_lo']:.4f}, {e['x_hi']:.4f}] | {fmt(e['swing_others'])} |")
        L_.append("")
    if C is not None:
        lo = REF_HALF[0] - REF_HALF[1]
        L_ += [f"## C. A feature near the range (x_m = {REF_HALF[0]} + {REF_HALF[1]} sin ωt, the range [{lo:g}, {REF_HALF[0] + REF_HALF[1]:g}]; "
               f"the exact averaged loop, ε = {EPS_C:g}, to τ = {TAU_C:g}; the bump at 0.07 as given (−0.6) and flipped (+{FLIP_C:g}))", "",
               "| the first bump | below the range | the bump at 0.07 | the drive at its nearest kernel before any entry (of the largest) | "
               "entry near it (τ) | fitted | kernels | centres | coefficients | the largest sub-threshold \\|w\\| | the orbit's error rms at τ = "
               + ', '.join(f'{t:g}' for t in CHECK_C) + " |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for c in C:
            av = c['avg']
            below = f"{c['dist']:g} ({c['dist'] / SIG:.2g}σ)" if c['dist'] is not None else f"{lo - c['bump']:g} ({(lo - c['bump']) / SIG:.2g}σ)"
            L_.append(f"| {c['bump']:g} | {below} | {('+' + format(FLIP_C, 'g')) if c.get('flip') else ASTAR[1]:} | "
                      f"{c['drive_near']:.3g} ({c['drive_near'] / c['drive_max']:.3g}) | "
                      f"{fmt(c['t_first']) if np.isfinite(c['t_first']) else 'none'} | {'yes' if c['fitted'] else 'no'} | {av['kernels']} | "
                      f"{rnd(av['centres'])} | {rnd(av['w'])} | {fmt(av.get('sub'))} | {', '.join(fmt(v) for t, v in c.get('err_hist', []))} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. The positions' metric (the exact averaged loop, to τ = {TAU_D:g})", "",
               "| metric | ε | formed | clusters | centres | the largest sub-threshold \\|w\\| | rms of the orbit's error |", "|---|---|---|---|---|---|---|"]
        for r in D:
            L_.append(f"| {r['metric']} | {r['eps']:g} | {r['n_act']} | {r['clusters']} | {sorted(set(round(v, 3) for v in r['centres']))} | "
                      f"{fmt(r.get('sub'))} | {fmt(r['rms'])} |")
        zs = [r for r in D if 'zero_grid' in r]
        if zs:
            L_ += ["", "Fixed centres: the exact averaged drive's zero (the closed loop's least-squares fit along the orbit) on the formed kernels "
                   "and on the whole grid — the orbit's error rms (the drive's largest remaining component): " + '; '.join(
                       f"ε = {r['eps']:g}: formed {fmt(r['zero_formed'][1])} ({fmt(r['zero_formed'][2])}), grid {fmt(r['zero_grid'][1])} "
                       f"({fmt(r['zero_grid'][2])})" for r in zs) + "."]
        L_.append("")
    if E_ is not None:
        L_ += [f"## E. Measurement noise (ε = {EPS_E:g}, τ ≤ {TAU_E:g}; 16 tones in 10–30 rad/s of rms s in the measured state)", "",
               "| s | leak λΛ | regressor | kernels | pushed out | max \\|w\\| | max \\|e\\| (measured), last period | true: max, rms | the in-range "
               "fit's rms error | the largest sub-threshold \\|w\\| | entries (kernel: τ) |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in list(E_) + ([] if E2 is None else list(E2)) + ([] if E3 is None else list(E3)):
            L_.append(f"| {r['s']:g} | {r['lam']:g} | {'x_m (desired)' if r.get('regx') else 'x̃'} | {r['kernels']} | {r['outside']} | "
                      f"{r['wmax']:.3g} | {fmt(r['e_last'])} | {fmt(r.get('e_true'))}, {fmt(r.get('e_true_rms'))} | {fmt(r.get('fit_rms'))} | "
                      f"{fmt(r.get('sub'))} | {ev(r['entries'])} |")
        if E2 is not None and len(E2):
            L_ += ["", f"The rows with λΛ = {', '.join(f'{v:g}' for v in LEAKS_E2)} (E2): the leak's strength between none and {LEAK_E:g}, "
                   "without noise and at s = 0.01 — the mirror leak's bias (the amplitudes shrink toward the leak's centre) against its "
                   "hold on the kernels' positions."]
        if E3 is not None and len(E3):
            L_ += ["", "E3: the regressor at the reference (the desired regressor of D‴.31) at s = 0, 0.01, 0.02, and the mirror leak at "
                   + ', '.join(f'{lam:g}' for lam in sorted(set(r['lam'] for r in E3 if r['lam'] > 0))) + " at s = 0.02. The dead zones are "
                   "Proposition D‴.57's "
                   "(`exp_dead_zone.py`, Part C)."]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_closed_transit.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'closed_transit_raw.npz')
    md = os.path.join(OUT, 'closed_transit.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'E2', 'B2', 'E3')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = {k: v for k, v in old.items() if k in keys and k not in parts}
    t0 = time.time()
    for k, fn in (('A', part_A), ('C', part_C), ('D', part_D), ('E', part_E), ('E2', part_E2), ('E3', part_E3), ('B', part_B),
                  ('B2', part_B2)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
