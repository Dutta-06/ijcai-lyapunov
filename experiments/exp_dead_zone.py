"""
Proposition D‴.57: the transit's dead zone — gated by the residual's energy the free law stops at the discrepancy point; gating every
particle freezes the last block before it has relaxed, gating only the unformed ones gives the oracle's fit
====================================================================================================================================

D‴.56's setting (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least 0.45 apart
with m > 0.8q; data y = f* + ξ, ξ i.i.d. N(0, s²)), the gate open while the residual's rms exceeds (1 + δ)s:

A. The sliding limit: the noisy path; the dead zone stops it at the first entry k whose relaxed residual's rms is ≤ (1 + δ)s (the
   discrepancy rule, D‴.56 D). With every particle gated (the energy dead zone) the k-th block's fast flow from its entrant's seed
   (`critical_rhs` on the block) stops where the rms falls to (1 + δ)s: the frozen fit's error against the oracle's (the k-th block
   relaxed, the k-kernel least-squares fit of y) and against s(4K/N + 2δ + δ²)^{1/2}; the frozen residual's tangent part u = r_c − r_K,
   |u|² against (1 + δ)²s²N − |r_K|², the cross term ⟨f̂_K − f*, u⟩ and u's direction across a target's draws with the same last
   entrant, the last entrant's amplitude deficit δw (at first order u = δw·P⊥φ_last) — every statistic split by whether the noisy path
   keeps the noiseless entry order (D‴.56's classification; referee round 50: a swapped order correlates the oracle's error with u).
   With only the unformed particles gated (the formation dead zone) the block relaxes: the oracle's fit. A3: the first NT_A2 targets on
   N_A2 samples at s = S_A2 — the excess over the oracle grows like N^{1/2}, the amplitude deficit does not shrink.
B. The loop at the mass EPS_B from α = 0 on N_B targets × N_DRAW_B draws at LEVELS_B: ungated to τ_K + 3W (W the path's window to its
   first noise kernel), with the energy dead zone and with the formation dead zone (|w_i| ≤ THETA_B gated): the kernels formed (|w| >
   THETA_B), the fits' errors and their distance from the sliding limit's fits; the noise drawn afresh at each level.

C. The closed loop of Proposition D‴.49 under its measurement noise (Part E there: 16 tones in 10–30 rad/s of rms s, ε = EPS_C, to
   τ = TAU_C), the regressor at the measured state or at the reference (D‴.31's desired regressor): ungated, with D‴.49's mirror leak
   (λΛ = LEAK_C), and with the dead zone on the error's energy filtered over TF_C, gated between the rms GATE_C·s (every kernel, or
   only the unformed ones): the kernels formed, those pushed beyond the range, the coefficients and the tracking error at the end;
   and the noise-free loop with either regressor as the reference. Dense outputs (N_OUT_C times): the gate's state over time (the
   fractions shut and fully open, the first closing, the last reopening), each kernel's formation and its ejection beyond the range,
   and the fit's rms error over the last window (median, 5–95%). At s = S_EXTRA_C with the desired regressor, EXTRA_C's tests of the
   mechanism (referee round 50): the adaptation driven by the true error, the centres bounded, other gate thresholds and filters, another
   noise seed. C2, ungated at x_m: the noise's structure — one tone (no difference frequencies), two tones far apart (a fast beat), two
   close ones (a slow beat), the frequencies incommensurate with the reference's —, and the centres' metric: Wasserstein–Fisher–Rao
   (mobility G·a, bounded) against the paradigm's (G·β·a = Gm|w|, growing with the coefficient), noise-free and noisy.
D. From the fitted state (each feature's kernel at its centre with its amplitude): the closed loop at x_m, ungated, under D‴.49's noise
   for D_PERIODS periods — which formed kernels leave their feature, at which level, against the amplitudes (λA*) and the metric.

Usage: python experiments/exp_dead_zone.py [--only A,A2,B,C,C2,D] [--from-raw] → results/chaos/dead_zone.md (+ _raw.npz).
"""
import os
import sys
import time
import warnings

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X
import exp_noisy_certificate as NC
import exp_transit_mass as TM

OUT = os.environ.get('DZ_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('DZ_WORKERS', '2'))
DELTA = 0.1
LEVELS_A, N_DRAW_A = (0.01, 0.03, 0.1), 10
N_A2, S_A2, NT_A2, ND_A2 = (61, 181, 541, 1621), 0.03, 20, 20          # A3: sample counts (61 the setting's), level, targets, draws
EPS_B, N_B, N_DRAW_B, LEVELS_B, THETA_B = 1e-8, 16, 3, (0.01, 0.03, 0.1), 1e-3
T_RELAX = 1e8                                                          # the block's fast flow when the gate does not close
EPS_C, TAU_C, SIG_C, LEAK_C, TF_C, GATE_C = 1e-4, 600.0, (0.01, 0.02), 3e-3, 5.0, (1.2, 1.5)
N_OUT_C, WIN_C = 1201, 0.5                                             # C: output times on [0, T]; the last window, a share of T
S_EXTRA_C = 0.02
EXTRA_C = (('true error', 'none', dict(etrue=True)), ('box 1.2', 'none', dict(cbox=1.2)), ('box 1.2', 'energy', dict(cbox=1.2)),
           ('thresholds 2.0s, 2.5s', 'energy', dict(gate=(2.0, 2.5))), ('thresholds 2.0s, 2.5s', 'formation', dict(gate=(2.0, 2.5))),
           ('T_f = 40', 'energy', dict(Tf=40.0)), ('T_f = 40', 'formation', dict(Tf=40.0)),
           ('noise seed 8', 'none', dict(seed=8)), ('noise seed 8', 'energy', dict(seed=8)))
C2_RUNS = (('one tone, 21.37 rad/s', 0.02, dict(freqs=(21.37,))), ('one tone, 21.37 rad/s', 0.04, dict(freqs=(21.37,))),
           ('two tones, 13.1 and 27.9 rad/s', 0.02, dict(freqs=(13.1, 27.9))), ('two tones, 21.37 and 21.44 rad/s', 0.02, dict(freqs=(21.37, 21.44))),
           ('WFR metric', 0.0, dict(metric='wfr')), ('WFR metric', 0.02, dict(metric='wfr')), ('WFR metric', 0.03, dict(metric='wfr')))
D_LAM, D_PERIODS, D_LEFT = (0.5, 0.75, 1.0, 1.5), 300, 0.3                # D: amplitude scales, periods of the reference, the distance left
D_LEVELS = (0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.04)
D_WFR_LEVELS = (0.02, 0.04, 0.06, 0.08, 0.12, 0.16)


def stop_law(N, nK, delta=DELTA):
    """The stop at the K-th entry at first order: |r_K|²/s² ~ χ²_{N−2K} and the stop is the event |r_K|² ≤ cs², c = (1 + δ)²N. Returns
    c, its probability p, m = E[χ²_{N−2K} | ≤ c] (from E[X1{X ≤ c}] = kP(χ²_{k+2} ≤ c), X ~ χ²_k), and given the stop E|u|²/(s²N) =
    (c − m)/N and E err_dz²/s² = (2K + c − m)/N."""
    from scipy.stats import chi2
    c = (1 + delta) ** 2 * N
    p = chi2.cdf(c, N - 2 * nK)
    m = (N - 2 * nK) * chi2.cdf(c, N - 2 * nK + 2) / p
    return dict(c=c, p=p, m=m, u=(c - m) / N, dz=(2 * nK + c - m) / N)


def gated_relax(L, K, y, c_blk, w_blk, thr):
    """The block's fast flow on the critical curve from (c_blk, w_blk) until the residual's rms falls to thr (event), or T_RELAX.
    Returns (crossed, c, w)."""
    n = len(c_blk)
    al = np.sign(w_blk) * RT._alpha_of(L, np.abs(w_blk), RT.coef_crit)
    rhs = RT.critical_rhs(L, K, y)

    def ev(t, z):
        return float(np.sqrt(np.mean((y - K.phi(z[n:]) @ RT.critical_w(L, z[:n])) ** 2))) - thr
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, T_RELAX), np.r_[al, c_blk], method='LSODA', rtol=1e-10, atol=1e-13, events=ev)
    z = sol.y[:, -1]
    return bool(sol.t_events[0].size), z[n:], RT.critical_w(L, z[:n])


def a_draw(L, K, c0, f, order, xi, s, nK):
    """One draw: the entries classified against the noiseless order (D‴.56's `classify`: the set, the order), the stop, the energy dead
    zone's frozen fit and the oracle's, the tangent part u, the last entrant and its amplitude deficit δw = sign(w_j)(w_j^oracle −
    w_j^frozen), with P⊥φ_j its kernel off the rest of the block's tangent at the oracle (the other kernels' coefficients and centres, its
    own centre): at first order u = δw·P⊥φ_j. The fits f̂_dz, f̂_or and u are returned as arrays (`a_job` drops them)."""
    y = f + xi
    thr = (1 + DELTA) * s
    # to the (K + 1)-th entry, not relaxed (a noise kernel's relaxation is slow and not needed: a stop after the K-th counts as late)
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)
    ent, exits, exact, exact_order = NC.classify(R, order)
    d = dict(exact=exact, exact_order=exact_order, exits=exits)
    if not exact:
        return d
    rms = [R['hist'][NC.fit_after(ent, i)][3] for i in range(1, nK + 1)]   # relaxed after entries 1, …, K
    k = next((i + 1 for i, r_ in enumerate(rms) if r_ <= thr), -1)
    d.update(k=k, last=int(ent[nK - 1][1]))
    if k < 1:
        return d
    # the k-th block before it relaxes: the state before the k-th entry's iteration (the (k − 1)-th block relaxed) and the entrant
    tau, c_prev, w_prev, _ = R['hist'][ent[k - 1][2]]
    j = ent[k - 1][1]
    r_prev = y - K.phi(c_prev) @ w_prev
    sg = np.sign(K.phi(c_prev)[:, j] @ r_prev)
    act = np.flatnonzero(w_prev != 0)
    blk = np.r_[act, j]
    w0 = np.r_[w_prev[act], sg * 1e-7]
    crossed, cc, wc = gated_relax(L, K, y, c_prev[blk], w0, thr)
    _, cK, wK, rmsK = R['hist'][NC.fit_after(ent, k)]                 # the k-th block relaxed (polished): the oracle at k = K
    fK, fc = K.phi(cK) @ wK, K.phi(cc) @ wc
    u = fK - fc                                                        # r_c − r_K
    N = y.size
    blkK = np.flatnonzero(wK != 0)
    prev = blkK[blkK != j]
    Ph, dPh = K.phi(cK), K.dphi(cK)
    T = np.c_[Ph[:, prev], dPh[:, prev] * wK[prev], dPh[:, j] * wK[j]]
    pj = Ph[:, j] - T @ np.linalg.lstsq(T, Ph[:, j], rcond=None)[0]
    d.update(crossed=crossed, err_dz=float(np.sqrt(np.mean((fc - f) ** 2))), err_or=float(np.sqrt(np.mean((fK - f) ** 2))),
             u2=float(u @ u), u2_pred=float(thr ** 2 * N - np.sum((y - fK) ** 2)),
             cross=float((fK - f) @ u / (np.linalg.norm(fK - f) * np.linalg.norm(u) + 1e-300)),
             dw=float(np.sign(wK[j]) * (wK[j] - wc[-1])), wj=float(abs(wK[j])), pnorm=float(np.linalg.norm(pj)),
             cos_last=float(abs(u @ pj) / (np.linalg.norm(u) * np.linalg.norm(pj) + 1e-300)), u=u, fc=fc, fK=fK)
    return d


def a_job(args):
    q, cs, As = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    order = [e[1] for e in R0['events'] if e[2] == 'enter'][:nK]
    rng = np.random.default_rng(5700 + q)
    out = dict(q=q, levels=[])
    for s in LEVELS_A:
        draws = [a_draw(L, K, c0, f, order, s * rng.standard_normal(nx), s, nK) for _ in range(N_DRAW_A)]
        # u's direction across the draws with the same last entrant: |cos| against the first such draw's (the reference itself left out)
        ref = {}
        for d in draws:
            if d.get('k') == nK and 'u' in d:
                if d['last'] not in ref:
                    ref[d['last']] = d['u']
                    d['cos'] = np.nan
                else:
                    u0 = ref[d['last']]
                    d['cos'] = float(abs(d['u'] @ u0) / (np.linalg.norm(d['u']) * np.linalg.norm(u0) + 1e-300))
        for d in draws:
            for key in ('u', 'fc', 'fK'):
                d.pop(key, None)
        out['levels'].append(dict(s=s, draws=draws))
    out['secs'] = time.time() - t0
    print(f"  A target {q}: " + ', '.join(f"{sum(d.get('k') == nK for d in lv['draws'])}/{len(lv['draws'])}" for lv in out['levels'])
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_A():
    from multiprocessing import Pool
    jobs = [(q, cs, As) for q, (cs, As) in enumerate(NC.targets())]
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, jobs, chunksize=1)


def a2_job(args):
    """A3: target q on N equispaced samples of [−1, 1] (the kernels and the feature unchanged): the noiseless path (exact if its K-th
    relaxed fit has rms ≤ 10⁻⁸max|f*| with no exit), then ND_A2 draws at S_A2 as in Part A."""
    q, cs, As, N = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, _, c0 = TM.setup()
    K = RT.Bumps(np.linspace(-1, 1, N), X.SIG_C)
    f = X.bumps_target(K, cs, As)
    nK = len(cs)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)
    order = [e[1] for e in R0['events'] if e[2] == 'enter'][:nK]
    ent0, exits0, _, _ = NC.classify(R0, order)
    exact0 = len(order) == nK and exits0 == 0 and R0['hist'][NC.fit_after(ent0, nK)][3] <= 1e-8 * np.abs(f).max()
    out = dict(q=q, N=N, exact0=bool(exact0), draws=[])
    if exact0:
        rng = np.random.default_rng([5900, q, N])
        for _ in range(ND_A2):
            d = a_draw(L, K, c0, f, order, S_A2 * rng.standard_normal(N), S_A2, nK)
            for key in ('u', 'fc', 'fK'):
                d.pop(key, None)
            out['draws'].append(d)
    out['secs'] = time.time() - t0
    print(f"  A2 target {q} N {N}: exact {exact0}, stop at K {sum(d.get('k') == nK for d in out['draws'])}/{len(out['draws'])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def part_A2():
    from multiprocessing import Pool
    jobs = [(q, cs, As, N) for N in N_A2[::-1] for q, (cs, As) in enumerate(NC.targets()[:NT_A2])]
    with Pool(N_WORKERS) as pool:
        return pool.map(a2_job, jobs, chunksize=1)


def loop_gated(L, K, y, c0, eps, tau_max, thr, mode):
    """The loop from α = 0 at the masses ε to τ_max (τ = t/Λ(ε)): mode 'none' ungated; 'energy' every particle's drive gated by the
    residual's rms (frozen once it falls to thr); 'formation' the unformed particles (|w| ≤ THETA_B) frozen at that moment, the formed
    block relaxing to τ_max. Returns (w, c, τ at which the gate closed or nan)."""
    d = len(c0)
    b0 = eps * np.ones(d)
    lam = RT.transit_scale(L, eps)
    rhs = RT.critical_rhs(L, K, y, b0=b0)
    z0 = np.r_[np.zeros(d), c0]
    wz = lambda z: RT.critical_w(L, z[:d], b0)
    if mode == 'none':
        sol = solve_ivp(rhs, (0, lam * tau_max), z0, method='LSODA', rtol=1e-9, atol=1e-12)
        z = sol.y[:, -1]
        return wz(z), z[d:], np.nan

    def ev(t, z):
        return float(np.sqrt(np.mean((y - K.phi(z[d:]) @ wz(z)) ** 2))) - thr
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, lam * tau_max), z0, method='LSODA', rtol=1e-9, atol=1e-12, events=ev)
    z = sol.y[:, -1]
    if not sol.t_events[0].size:
        return wz(z), z[d:], np.nan
    t_c = float(sol.t_events[0][0])
    if mode == 'formation':
        frozen = np.abs(wz(z)) <= THETA_B
        mask = np.r_[~frozen, ~frozen].astype(float)
        sol2 = solve_ivp(lambda t, zz: mask * rhs(t, zz), (t_c, max(lam * tau_max, t_c * 1.5)), z, method='LSODA', rtol=1e-9, atol=1e-12)
        z = sol2.y[:, -1]
    return wz(z), z[d:], t_c / lam


def b_noise(q, k_, s, nx):
    """Part B's draw k_ on target q at the level s: its own stream per level (referee round 50: the levels had shared one vector)."""
    return s * np.random.default_rng(5800 + 97 * q + k_ + 10_000 * LEVELS_B.index(s)).standard_normal(nx)


def b_job(args):
    q, cs, As, s, k_ = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    xi = b_noise(q, k_, s, nx)
    y = f + xi
    thr = (1 + DELTA) * s
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=nK + 4, tol=0.0, max_entries=nK + 1)
    ent = [e for e in R['events'] if e[2] == 'enter']
    out = dict(q=q, s=s, k=k_, exact=sorted(e[1] for e in ent[:nK]) == sorted(cert['e']) and not any(e[2] == 'exit' for e in R['events']))
    if not out['exact'] or len(ent) <= nK:
        return out
    tau_K, W = ent[nK - 1][0], ent[nK][0] - ent[nK - 1][0]
    tau_max = tau_K + 3 * W
    order = [e[1] for e in R0['events'] if e[2] == 'enter'][:nK]
    lim = a_draw(L, K, c0, f, order, xi, s, nK)                        # the sliding limit's stop and fits on the same draw
    out.update(tau_K=tau_K, W=W, k_lim=lim.get('k'), err_dz_lim=lim.get('err_dz', np.nan), err_or_lim=lim.get('err_or', np.nan))
    for mode in ('none', 'energy', 'formation'):
        w, c, tc = loop_gated(L, K, y, c0, EPS_B, tau_max, thr, mode)
        fit = K.phi(c) @ w
        out[mode] = dict(formed=int((np.abs(w) > THETA_B).sum()), err=float(np.sqrt(np.mean((fit - f) ** 2))), tau_gate=tc,
                         rms=float(np.sqrt(np.mean((y - fit) ** 2))),
                         d_dz=float(np.sqrt(np.mean((fit - lim['fc']) ** 2))) if 'fc' in lim else np.nan,
                         d_or=float(np.sqrt(np.mean((fit - lim['fK']) ** 2))) if 'fK' in lim else np.nan)
    out['secs'] = time.time() - t0
    print(f"  B target {q} s {s:g} draw {k_}: formed " + ', '.join(f"{m} {out[m]['formed']}" for m in ('none', 'energy', 'formation'))
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_B():
    from multiprocessing import Pool
    tg = NC.targets()[:N_B]
    jobs = [(q, cs, As, s, k_) for q, (cs, As) in enumerate(tg) for s in LEVELS_B for k_ in range(N_DRAW_B)]
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, jobs, chunksize=1)


def gate_stats(z, ts, lo, hi, Lam):
    """The gate's state along the dense outputs: the shares of output times with the gate shut (z ≤ lo²) and fully open (z ≥ hi²), the
    first closing and the last reopening (τ = t/Λ), and whether it is shut at the end."""
    shut, full = z <= lo ** 2, z >= hi ** 2
    up = np.flatnonzero(shut[:-1] & ~shut[1:])                          # shut at one output time, not at the next
    return dict(shut=float(shut.mean()), open=float(full.mean()), first_close=float(ts[np.argmax(shut)] / Lam) if shut.any() else np.nan,
                last_reopen=float(ts[up[-1] + 1] / Lam) if up.size else np.nan, shut_end=bool(shut[-1]))


def kernel_events(entries, ts, C, Lam, cstar, outside):
    """Each kernel's formation (its first |w| > th, located by the integrator) and its ejection — the first output time after the
    formation with |c| > outside —, with the feature nearest its centre at the formation; τ units, in the order of formation."""
    out = []
    for i in np.flatnonzero(np.isfinite(entries)):
        k0 = min(int(np.searchsorted(ts, entries[i])), ts.size - 1)
        beyond = np.flatnonzero(np.abs(C[k0:, i]) > outside)
        out.append(dict(i=int(i), formed=float(entries[i] / Lam), feature=int(np.argmin(np.abs(np.asarray(cstar) - C[k0, i]))),
                        ejected=float(ts[k0 + beyond[0]] / Lam) if beyond.size else np.nan))
    return sorted(out, key=lambda e: e['formed'])


def noise_of(s, opts):
    """C's measurement noise of rms s: 16 tones in 10–30 rad/s (seed opts['seed'], default 7), or the tones at opts['freqs'] with equal
    amplitudes and random phases; None at s = 0."""
    from src.core import closed_loop as CL
    if s <= 0:
        return None
    if 'freqs' in opts:
        fr = np.asarray(opts['freqs'], float)
        return dict(w=fr, ph=np.random.default_rng(opts.get('seed', 7)).uniform(0, 2 * np.pi, fr.size), amp=np.full(fr.size, s * np.sqrt(2.0 / fr.size)))
    return CL.tones(s, 10.0, 30.0, seed=opts.get('seed', 7))


def c_job(args):
    """C: one closed loop (D‴.49's setting, the C integrator) — mode 'none', 'leak', 'energy' or 'formation', the regressor at x̃ or x_m,
    a variant `name` with its options (EXTRA_C: the adaptation driven by the true error, the centres bounded, the gate's thresholds or
    filter, the noise's seed; C2: the noise's frequencies `freqs`, equal amplitudes of total rms s, or the centres' metric), dense outputs."""
    s, regx, mode, name, opts = args
    import exp_closed_transit as CT
    from src.core import closed_loop as CL
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = CT.law()
    Lam = RT.transit_scale(L, EPS_C)
    T = Lam * TAU_C
    per = 2 * np.pi / CT.OM
    ts = np.unique(np.r_[np.linspace(0, T, N_OUT_C), np.linspace(T - per, T, 60)])
    th = CT.th_abs()
    lo_f, hi_f = opts.get('gate', GATE_C)
    gate = (dict(mode=mode, Tf=opts.get('Tf', TF_C), lo=lo_f * s, hi=hi_f * s, theta=th) if mode in ('energy', 'formation') and s > 0
            else None)
    r = CL.simulate(L, EPS_C, CT.C0, CT.SIG, (CT.CSTAR, CT.ASTAR, CT.SIG), CT.KE, (0.0, 1.0, CT.OM, 0.0, np.sqrt(2) * CT.OM), ts,
                    noise=noise_of(s, opts), th=th, regx=regx, metric=opts.get('metric', 'paradigm'),
                    lam=(LEAK_C / Lam if mode == 'leak' else 0.0), gate=gate, etrue=opts.get('etrue', False),
                    cbox=opts.get('cbox', np.inf))
    act = np.abs(r['w'][-1]) > th
    gs = gate_stats(r['z'], ts, lo_f * s, hi_f * s, Lam) if gate is not None else None
    xr = np.linspace(-1, 1, 201)                                       # the reference's range
    phi = lambda c: np.exp(-(xr[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * CT.SIG ** 2))
    truth = phi(CT.CSTAR) @ CT.ASTAR
    fit = lambda k: float(np.sqrt(np.mean((phi(r['c'][k]) @ r['w'][k] - truth) ** 2)))
    win = np.flatnonzero(ts >= T * (1 - WIN_C))
    fw = np.array([fit(k) for k in win])
    inw = np.array([int(np.sum((np.abs(r['w'][k]) > th) & (np.abs(r['c'][k]) <= CT.OUTSIDE))) for k in win])
    ev = kernel_events(r['entries'], ts, r['c'], Lam, CT.CSTAR, CT.OUTSIDE)
    out = dict(s=s, regx=bool(regx), mode=mode, name=name, n_act=int(act.sum()), outside=int(np.sum(np.abs(r['c'][-1][act]) > CT.OUTSIDE)),
               w=r['w'][-1][act].tolist(), c=r['c'][-1][act].tolist(), wmax=float(np.abs(r['w'][-1]).max()), e_last=r['e_last'],
               e_true=r['e_true_rms_last'], e_true_max=r['e_true_last'], fit_rms=fit(-1), w_all=r['w'][-1].tolist(),
               c_all=r['c'][-1].tolist(), tau_closed=gs['first_close'] if gs else np.nan, gate=gs, events=ev,
               fit_win=(float(np.median(fw)), float(np.quantile(fw, 0.05)), float(np.quantile(fw, 0.95))), formed_win=float(np.median(inw)),
               cmax=float(np.abs(r['c']).max()), secs=time.time() - t0)
    print(f"  C s {s:g} regx {regx} {mode} {name}: {out['n_act']} kernels, {out['outside']} outside, formed {len(ev)}, ejected "
          f"{sum(np.isfinite(e['ejected']) for e in ev)}, fit (window) {out['fit_win'][0]:.4f} [{out['secs']:.0f}s]", flush=True)
    return out


def part_C():
    from multiprocessing import Pool
    jobs = ([(0.0, regx, 'none', 'base', {}) for regx in (False, True)]                # the noise-free references
            + [(s, regx, mode, 'base', {}) for regx in (False, True) for s in SIG_C for mode in ('none', 'leak', 'energy', 'formation')]
            + [(S_EXTRA_C, True, mode, name, opts) for name, mode, opts in EXTRA_C])
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=1)


def part_C2():
    from multiprocessing import Pool
    jobs = [(s, True, 'none', f"{name}, s = {s:g}", opts) for name, s, opts in C2_RUNS]
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=1)


def fitted_state(L, eps, cstar, astar):
    """D's start: for each feature the kernel started nearest it, moved to the feature's centre with its coefficient equal to the
    feature's amplitude (α by root-finding on its curve from the masses eps); the other kernels at α = 0 on their grid points."""
    from scipy.optimize import brentq
    import exp_closed_transit as CT
    d = CT.C0.size
    wfun = lambda al: (eps + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    al0, c0, idx = np.zeros(d), CT.C0.copy(), []
    for cs, As in zip(cstar, astar):
        i = int(np.argmin(np.abs(CT.C0 - cs)))
        al0[i] = np.sign(As) * brentq(lambda a: wfun(a) - abs(As), 1e-12, 50 * L.A)
        c0[i] = cs
        idx.append(i)
    return al0, c0, idx


def d_job(args):
    """D: C's closed loop (x_m, ungated, D‴.49's 16 tones of rms s) from the fitted state with the amplitudes λA*, the centres' metric
    given, for D_PERIODS periods of the reference: per formed kernel the first period whose period-averaged centre is D_LEFT from its
    feature (−1: never), the largest such distance and the coefficient at the end; the other kernels' largest |w| at the end."""
    metric, lam, s = args
    import exp_closed_transit as CT
    from src.core import closed_loop as CL
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = CT.law()
    As_ = lam * CT.ASTAR
    al0, c0, idx = fitted_state(L, EPS_C, CT.CSTAR, As_)
    k = 8                                                              # output times per period
    ts = np.linspace(0, D_PERIODS * 2 * np.pi / CT.OM, D_PERIODS * k + 1)
    r = CL.simulate(L, EPS_C, c0, CT.SIG, (CT.CSTAR, As_, CT.SIG), CT.KE, (0.0, 1.0, CT.OM, 0.0, np.sqrt(2) * CT.OM), ts,
                    noise=noise_of(s, {}), th=CT.th_abs(As_), regx=True, al0=al0, metric=metric)
    dev = np.abs(r['c'][:-1].reshape(D_PERIODS, k, c0.size).mean(1)[:, idx] - CT.CSTAR)
    left = [int(np.argmax(dev[:, q] > D_LEFT)) + 1 if (dev[:, q] > D_LEFT).any() else -1 for q in range(len(idx))]
    others = [j for j in range(c0.size) if j not in idx]
    out = dict(metric=metric, lam=lam, s=s, A=As_.tolist(), left=left, dev_max=dev.max(0).tolist(),
               w_end=[float(r['w'][-1, i]) for i in idx], others=float(np.abs(r['w'][-1, others]).max()), secs=time.time() - t0)
    print(f"  D {metric} lam {lam} s {s:g}: left {left}, w end {', '.join(f'{v:+.3f}' for v in out['w_end'])} [{out['secs']:.0f}s]",
          flush=True)
    return out


def part_D():
    from multiprocessing import Pool
    jobs = [('paradigm', lam, s) for lam in D_LAM for s in D_LEVELS] + [('wfr', 1.0, s) for s in D_WFR_LEVELS]
    with Pool(N_WORKERS) as pool:
        return pool.map(d_job, jobs, chunksize=1)


def _mq(v, fmt='.3f'):
    """median [5%, 95%] of v"""
    v = np.asarray(v, float)
    return f"{np.median(v):{fmt}} [{np.quantile(v, 0.05):{fmt}}, {np.quantile(v, 0.95):{fmt}}]"


def report(A=None, B=None, C=None, A2=None, C2=None, D=None):
    nK, N = 3, 61
    L_ = ["# Proposition D‴.57: the transit's dead zone (`exp_dead_zone.py`)", "",
          "D‴.56's setting — 21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least "
          "0.45 apart with m > 0.8q; y = f* + ξ, ξ i.i.d. N(0, s²) —, the gate open while the residual's rms exceeds (1 + δ)s, "
          f"δ = {DELTA:g}.", ""]
    se = lambda v: float(np.std(v) / np.sqrt(len(v))) if len(v) > 1 else np.nan
    if A is not None:
        law = stop_law(N, nK)
        L_ += ["## A. The sliding limit", "",
               f"Per level s: {N_DRAW_A} draws on each of the {len(A)} targets; the draws whose noisy path enters E first (no exit before "
               "the (K + 1)-th entry), split by whether it keeps the noiseless entry order (D‴.56's classification): the dead zone's stop "
               "(the first entry whose relaxed rms is ≤ (1 + δ)s) at the K-th entry, later (none of the first K) or earlier; at the K-th, "
               "the energy dead zone's frozen fit against the oracle's (the K-th block relaxed): the errors' ratio, the mean square over "
               "its prediction given the stop (± its standard error) and over s²(4K/N + 2δ + δ²) (the unconditional value), the oracle's "
               "over s²·2K/N, and the correlation of the oracle's error with the frozen part u = r_c − r_K:", "",
               "| s | entry order | draws | stop at K, later, earlier | err_dz/err_or: median [5–95%] | mean err_dz²: / prediction given "
               "the stop; / s²(4K/N + 2δ + δ²) | mean err_or² / (s²2K/N) | corr(err_or, u): mean ± s.e. |",
               "|---|---|---|---|---|---|---|---|"]
        T2 = []
        for i, s in enumerate(LEVELS_A):
            for kept in (True, False):
                ds = [d for t in A for d in t['levels'][i]['draws'] if d['exact'] and d.get('exact_order') == kept]
                atK = [d for d in ds if d['k'] == nK and d.get('crossed')]
                later = sum(1 for d in ds if d['k'] == -1 or d['k'] > nK)
                earlier = sum(1 for d in ds if 0 < d['k'] < nK)
                lab_ = 'noiseless' if kept else 'swapped'
                head = f"| {s:g} | {lab_} | {len(ds)} | {len(atK)}, {later}, {earlier} | "
                if not atK:
                    L_.append(head + "— | — | — | — |")
                    T2.append(f"| {s:g} | {lab_} | 0 | — | — | — | — | — | — |")
                    continue
                e2 = np.array([d['err_dz'] ** 2 for d in atK]) / s ** 2
                eo = np.mean([d['err_or'] ** 2 for d in atK]) / (s ** 2 * 2 * nK / N)
                cr = [d['cross'] for d in atK]
                L_.append(head + f"{_mq([d['err_dz'] / d['err_or'] for d in atK])} | {e2.mean() / law['dz']:.3f} ± "
                          f"{se(e2) / law['dz']:.3f}; {e2.mean() / (4 * nK / N + 2 * DELTA + DELTA ** 2):.3f} | {eo:.3f} | "
                          f"{np.mean(cr):+.3f} ± {se(cr):.3f} |")
                ur = np.array([d['u2'] / d['u2_pred'] for d in atK if d['u2_pred'] > 0])
                cs_ = np.array([d['cos'] for d in atK if np.isfinite(d.get('cos', np.nan))])
                al = [d['dw'] * d['pnorm'] / np.sqrt(d['u2']) for d in atK]
                T2.append(f"| {s:g} | {lab_} | {len(atK)} | {_mq(ur, '.4f')}; {np.abs(ur - 1).max():.1e} | "
                          f"{np.sqrt(np.mean([d['u2'] for d in atK]) / N) / s:.3f} ({np.sqrt(law['u']):.3f}) | "
                          + (f"{np.median(cs_):.4f}, {np.quantile(cs_, 0.05):.4f}, {cs_.min():.4f} ({cs_.size})" if cs_.size else "— (0)")
                          + f" | {_mq([d['dw'] / s for d in atK], '.2f')} | {_mq([d['dw'] / d['wj'] for d in atK])} | "
                          f"{np.median(al):.4f} [{min(al):.4f}, {max(al):.4f}]; {min(d['cos_last'] for d in atK):.4f} |")
        L_ += ["", "The frozen part u at the K-th stop: the identity |u|² = (1 + δ)²s²N − |r_K|² (first order) with its largest deviation, "
                   "the rms of |u|/(s√N) against its prediction given the stop, |cos| of u against the target's first draw with the same "
                   "last entrant (that draw left out), and the last entrant's amplitude deficit δw = sign(w_last)(w_last of the oracle − "
                   "w_last frozen), P⊥φ_last its kernel off the rest of the block's tangent at the oracle (the other kernels' coefficients "
                   "and centres, its own centre) — at first order u = δw·P⊥φ_last:", "",
               "| s | entry order | at K | \\|u\\|² / ((1 + δ)²s²N − \\|r_K\\|²): median [5–95%]; max\\|· − 1\\| | rms \\|u\\|/(s√N) "
               "(prediction given the stop) | \\|cos\\| of u, same last entrant: median, 5%, min (draws) | δw/s: median [5–95%] | "
               "δw/\\|w_last\\|: median [5–95%] | δw\\|P⊥φ_last\\|/\\|u\\|: median [range]; min \\|cos(u, P⊥φ_last)\\| |",
               "|---|---|---|---|---|---|---|---|---|"] + T2
        L_ += ["", f"At first order |u|²/N = (1 + δ)²s² − |r_K|²/N, of mean s²(2δ + δ² + 2K/N); without the cross term the errors' rms "
                   f"ratio would be ((4K/N + 2δ + δ²)/(2K/N))^{{1/2}} = {np.sqrt((4 * nK / N + 2 * DELTA + DELTA ** 2) / (2 * nK / N)):.3f}. "
                   f"The stop at the K-th entry is the event |r_K|² ≤ (1 + δ)²s²N, |r_K|²/s² ~ χ²_{{N−2K}} (probability {law['p']:.3f}); "
                   f"given it E|r_K|² = s²(N − 2K)P(χ²_{{N−2K+2}} ≤ c)/P(χ²_{{N−2K}} ≤ c), c = (1 + δ)²N, = {law['m']:.2f}s² against "
                   f"{N - 2 * nK}s², so E|u|²/N = {law['u']:.4f}s² (rms |u|/(s√N) = {np.sqrt(law['u']):.3f}) and E err_dz² = s²(2K + c − "
                   f"E[χ²_{{N−2K}} | ≤ c])/N = {law['dz']:.4f}s² (×{law['dz'] / (4 * nK / N + 2 * DELTA + DELTA ** 2):.3f} the "
                   f"unconditional value; the rms ratio {np.sqrt(law['dz'] / (2 * nK / N)):.3f}). With u = δw·P⊥φ_last, E(δw|P⊥φ_last|)² "
                   f"= s²(c − E[χ²_{{N−2K}} | ≤ c]) = {law['c'] - law['m']:.2f}s²: both |P⊥φ_last|² and c − E[χ²_{{N−2K}} | ≤ c] grow like "
                   "N, so the deficit δw stays O(s) as N grows while the oracle's error falls like s(2K/N)^{1/2}: the rms ratio grows like "
                   "(2 + N(2δ + δ²)/(2K))^{1/2} (A3).", ""]
    if A2 is not None:
        L_ += ["### A3. The number of samples", "",
               f"The first {NT_A2} targets on N equispaced samples of [−1, 1] (the kernels and the features unchanged), s = {S_A2:g}, "
               f"{ND_A2} draws per target (the targets whose noiseless path recovers the features exactly); the draws in the noiseless "
               "entry order (the swapped ones counted), at the stop at the K-th entry: the mean squared errors against their predictions "
               "given the stop (stop_law(N)), the errors' rms ratio, and the last entrant's amplitude deficit:", "",
               "| N | targets: exact noiseless path | draws in the noiseless order: stop at K, later, earlier; swapped | mean err_dz²/s²: "
               "measured (prediction) | mean err_or² / (s²2K/N) | rms err_dz / rms err_or: measured (prediction) | δw/s: median "
               "[5–95%] | rms δw\\|P⊥φ_last\\|/s: measured (prediction (c − E[χ²_{N−2K} \\| ≤ c])^{1/2}) |",
               "|---|---|---|---|---|---|---|---|"]
        for N_ in N_A2:
            rs = [r for r in A2 if r['N'] == N_]
            if not rs:
                continue
            ok = [r for r in rs if r['exact0']]
            dd = [d for r in ok for d in r['draws'] if d['exact']]
            ds = [d for d in dd if d.get('exact_order')]
            atK = [d for d in ds if d['k'] == nK and d.get('crossed')]
            head = (f"| {N_} | {len(ok)}/{len(rs)} | {len(atK)}, {sum(1 for d in ds if d['k'] == -1 or d['k'] > nK)}, "
                    f"{sum(1 for d in ds if 0 < d['k'] < nK)}; {len(dd) - len(ds)} | ")
            if not atK:
                L_.append(head + "— | — | — | — | — |")
                continue
            lw = stop_law(N_, nK)
            e2 = np.mean([d['err_dz'] ** 2 for d in atK]) / S_A2 ** 2
            eo = np.mean([d['err_or'] ** 2 for d in atK]) / S_A2 ** 2
            amp = np.sqrt(np.mean([(d['dw'] * d['pnorm']) ** 2 for d in atK])) / S_A2
            L_.append(head + f"{e2:.4f} ({lw['dz']:.4f}) | {eo / (2 * nK / N_):.3f} | {np.sqrt(e2 / eo):.3f} "
                      f"({np.sqrt(lw['dz'] / (2 * nK / N_)):.3f}) | {_mq([d['dw'] / S_A2 for d in atK], '.2f')} | {amp:.3f} "
                      f"({np.sqrt(lw['c'] - lw['m']):.3f}) |")
        L_.append("")
    if B is not None:
        rs = [r for r in B if r.get('exact') and 'none' in r]
        L_ += ["## B. The loop", "",
               f"The loop at ε = {EPS_B:g} from α = 0 on {N_B} targets × {N_DRAW_B} draws at s = {', '.join(f'{s:g}' for s in LEVELS_B)} "
               "(the draws whose noisy path is exactly E; the noise drawn afresh at each level), to τ_K + 3W: on the draws whose sliding "
               f"limit stops at the K-th entry, the kernels formed (|w| > {THETA_B:g}) and the fits' rms errors over the oracle's (the "
               "limit's K-th block relaxed), the energy dead zone's against the ungated loop's, and the gated loops' fits against the "
               "sliding limit's (the energy dead zone's frozen fit f̂_dz, the oracle's f̂_or: the rms of the difference); on the late draws "
               "(the stop after the K-th entry) the kernels formed:", "",
               "| s | draws: stop at K, late | formed at K: ungated, energy, formation (median [range]) | err / err_or (medians): ungated, "
               "energy, formation | energy against the ungated (median) | the loop against the limit: max rms(f̂ − f̂_dz) (energy), "
               "max rms(f̂ − f̂_or) (formation) | the gate's τ / τ_K (median) | formed, late: ungated, energy, formation |",
               "|---|---|---|---|---|---|---|---|"]
        for s in LEVELS_B:
            ss = [r for r in rs if r['s'] == s]
            atK = [r for r in ss if r.get('k_lim') == nK]
            late = [r for r in ss if r.get('k_lim') != nK]
            if not atK:
                continue
            fm = {m: [r[m]['formed'] for r in atK] for m in ('none', 'energy', 'formation')}
            med = lambda v: f"{np.median(v):.3f}"
            L_.append(f"| {s:g} | {len(atK)}, {len(late)} | " + '; '.join(f"{int(np.median(fm[m]))} [{min(fm[m])}, {max(fm[m])}]" for m in fm)
                      + " | " + ', '.join(med([r[m]['err'] / r['err_or_lim'] for r in atK]) for m in ('none', 'energy', 'formation')) + " | "
                      + med([r['energy']['err'] / r['none']['err'] for r in atK]) + " | "
                      + f"{max(r['energy']['d_dz'] for r in atK):.1e}, {max(r['formation']['d_or'] for r in atK):.1e} | "
                      + f"{np.median([r['energy']['tau_gate'] / r['tau_K'] for r in atK if np.isfinite(r['energy']['tau_gate'])]):.4f} | "
                      + ('; '.join(', '.join(str(r[m]['formed']) for m in ('none', 'energy', 'formation')) for r in late) or '—') + " |")
        L_.append("")
    if C is not None:
        import exp_closed_transit as CT
        from src.core import closed_loop as CL
        lab = dict(none='ungated', leak=f'the mirror leak λΛ = {LEAK_C:g}', energy='the energy dead zone', formation='the formation dead zone')
        modes = ['none', 'leak', 'energy', 'formation']
        base = sorted((r for r in C if r.get('name', 'base') == 'base'), key=lambda r: (r['regx'], r['s'], modes.index(r['mode'])))
        extra = [r for r in C if r.get("name", "base") != "base"] + (list(C2) if C2 is not None else [])
        reg = lambda r: 'x_m' if r['regx'] else 'x̃'

        def crow(r):
            ev, g, fw = r['events'], r['gate'], r['fit_win']
            gs = ('— | —' if g is None else f"{g['shut']:.2f}, {g['open']:.2f} | "
                  + ', '.join('—' if not np.isfinite(v) else f"{v:.1f}" for v in (g['first_close'], g['last_reopen'])))
            return (f"{len(ev)}, {sum(np.isfinite(e['ejected']) for e in ev)} | {r['n_act']}, {r['outside']} | {fw[0]:.4f} [{fw[1]:.4f}, "
                    f"{fw[2]:.4f}] | {r['formed_win']:g} | {r['e_true']:.4f} | {r['e_last']:.4f} | {gs} |")
        cols = ("formed, ejected | at the end: formed, outside | fit rms on [−1, 1], last half: median [5–95%] | formed inside, last "
                "half (median) | true error rms, last period | max\\|ẽ\\|, last period | the gate: shut, fully open (shares of time) | "
                "the gate: first closing, last reopening (τ) |")
        L_ += ["## C. The closed loop under measurement noise", "",
               f"D‴.49's Part E (k_e = 10, ω = 0.2, three bumps of width 0.12; 21 kernels from α = 0 at ε = {EPS_C:g} to τ = {TAU_C:g}; the "
               "noise 16 tones in 10–30 rad/s of rms s, seed 7); the dead zone on the error's energy filtered over "
               f"T_f = {TF_C:g} time units, the drive gated smoothly between the rms {GATE_C[0]:g}s and {GATE_C[1]:g}s (every kernel, or "
               "the kernels with |w| below about the threshold); formed: |w| > 1% of max|A*| (a kernel's first crossing), ejected: |c| > "
               f"{CT.OUTSIDE:g} at an output time after its formation, outside: |c| > {CT.OUTSIDE:g} at the end; {N_OUT_C} output times on "
               f"[0, T]: the fit's rms error on the reference's range [−1, 1] (every kernel) and the formed kernels inside over the last "
               f"{WIN_C:.0%} of the run, the true error x − x_m's rms and the measured error's largest value over the last period, the "
               "gate's state over time (shut: the filtered energy at or below the lower threshold; fully open: at or above the upper):", "",
               "| regressor | s | adaptation | " + cols, "|---|---|---|---|---|---|---|---|---|---|---|"]
        L_ += [f"| {reg(r)} | {r['s']:g} | {lab[r['mode']]} | " + crow(r) for r in base]
        L_ += ["", "The formed kernels at the end (snapshots of a moving state; the window's statistics above):", "",
               "| regressor | s | adaptation | coefficients | centres |", "|---|---|---|---|---|"]
        L_ += [f"| {reg(r)} | {r['s']:g} | {lab[r['mode']]} | {', '.join(f'{v:.3f}' for v in r['w'])} | "
               f"{', '.join(f'{v:.3f}' for v in r['c'])} |" for r in base]
        if extra:
            L_ += ["", f"Tests of the mechanism at x_m, s = {S_EXTRA_C:g} unless stated (referee round 50): the adaptation driven by the true "
                       "error x − x_m (the control and the gate by the measured one), the centres bounded to [−1.2, 1.2], other gate "
                       "thresholds and filters, another noise seed; C2: the noise's structure — one tone, two tones far apart, two close "
                       "ones (the same rms) — and the centres' Wasserstein–Fisher–Rao metric (mobility G·a, bounded, against the "
                       "paradigm's Gm|w|):", "",
                   "| variant | adaptation | " + cols, "|---|---|---|---|---|---|---|---|---|---|"]
            L_ += [f"| {r['name']} | {lab[r['mode']]} | " + crow(r) for r in extra]
        tl = [r for r in base if r['regx'] and r['s'] == S_EXTRA_C] + extra
        if tl:
            L_ += ["", f"The kernels' timelines at x_m, s = {S_EXTRA_C:g}: kernel: formed at τ (the feature nearest its centre then) → "
                       f"ejected beyond |c| = {CT.OUTSIDE:g} at τ:", "", "| run | timeline |", "|---|---|"]
            for r in tl:
                L_.append(f"| {r.get('name', 'base')}, {lab[r['mode']]} | "
                          + ('; '.join(f"{e['i']}: {e['formed']:.1f} ({CT.CSTAR[e['feature']]:+.2f})"
                                       + (f" → {e['ejected']:.1f}" if np.isfinite(e['ejected']) else '') for e in r['events']) or '—') + " |")
        nz = CL.tones(1.0, 10.0, 30.0, seed=7)
        fl = float(np.sqrt(np.sum(nz['amp'] ** 2 / 2 * nz['w'] ** 2 / (nz['w'] ** 2 + CT.KE ** 2))))
        L_ += ["", f"Near the fit the loop high-passes the noise in the measured error (ẽ = e + n, e's response to n through −k_eẽ: "
                   f"|iω/(iω + k_e)|), a floor of rms {fl:.3f}s with the desired regressor (at x̃ the regressor's noise adds f̂′(x)n): the "
                   f"energy dead zone starts to close when the slow error's rms falls to {np.sqrt(GATE_C[1] ** 2 - fl ** 2):.2f}s and is "
                   f"shut below {np.sqrt(GATE_C[0] ** 2 - fl ** 2):.2f}s."]
        ref = {r['regx']: r for r in base if r['s'] == 0 and r['mode'] == 'none'}
        dev = []
        for s_ in SIG_C:
            for regx in (False, True):
                r = next((r for r in base if r['s'] == s_ and r['regx'] == regx and r['mode'] == 'none'), None)
                if r is None or regx not in ref:
                    continue
                r0 = ref[regx]
                if r['n_act'] != r0['n_act'] or r['outside'] or r0['outside']:
                    dev.append(f"{reg(r)}, s = {s_:g}: {r['n_act']} kernels against {r0['n_act']}")
                    continue
                o, o0 = np.argsort(r['c']), np.argsort(r0['c'])
                dev.append(f"{reg(r)}, s = {s_:g}: max|Δw| {np.abs(np.array(r['w'])[o] - np.array(r0['w'])[o0]).max():.2e}, "
                           f"max|Δc| {np.abs(np.array(r['c'])[o] - np.array(r0['c'])[o0]).max():.2e}")
        if dev:
            L_.append("The ungated loops at the end against the noise-free loop with the same regressor (the formed kernels in the order of "
                      "their centres): " + '; '.join(dev) + ".")
        L_.append("")
    if D is not None:
        L_ += ["## D. From the fitted state: the formed kernels under measurement noise", "",
               f"C's closed loop (the regressor at x_m, ungated, the 16 tones of rms s) started with each feature fitted by one kernel — the "
               "kernel started nearest it, at the feature's centre with the coefficient λA*_k —, for "
               f"{D_PERIODS} periods of the reference: per kernel the first period at which its period-averaged centre is {D_LEFT:g} from "
               "its feature (— if never) and the coefficient's growth w/(λA*) − 1 at the end; the paradigm's metric (the centres' "
               "mobility Gβa = Gm|w|) at λ = " + ', '.join(f'{v:g}' for v in D_LAM) + ", the Wasserstein–Fisher–Rao metric (G·a, at most "
               "G·a_max) at λ = 1:", "",
               "| metric | λ | s | " + ' | '.join(f"feature {q + 1}: left at period, w/(λA*) − 1" for q in range(3)) + " | the others' max\\|w\\| |",
               "|---|---|---|---|---|---|---|"]
        Ds = sorted(D, key=lambda r: (r['metric'] != 'paradigm', r['lam'], r['s']))
        for r in Ds:
            cells = [f"{'—' if r['left'][q] < 0 else r['left'][q]}, {r['w_end'][q] / r['A'][q] - 1:+.3f}" for q in range(3)]
            L_.append(f"| {r['metric']} | {r['lam']:g} | {r['s']:g} | " + ' | '.join(cells) + f" | {r['others']:.1e} |")
        th = []
        for metric in ('paradigm', 'wfr'):
            for lam in sorted({r['lam'] for r in Ds if r['metric'] == metric}):
                rs = [r for r in Ds if r['metric'] == metric and r['lam'] == lam]
                for q in range(3):
                    A = abs(rs[0]['A'][q])
                    s_ = next((r['s'] for r in rs if r['left'][q] >= 0), None)
                    th.append(f"{metric}, |A| = {A:.3g}: " + (f"{s_:g} (s·|A| = {s_ * A:.4f})" if s_ is not None
                                                             else f"kept to s = {rs[-1]['s']:g}"))
        bl = [(r['w_end'][q] / r['A'][q] - 1) / (r['A'][q] * r['s']) ** 2 for r in Ds if r['metric'] == 'paradigm' for q in range(3)
              if r['left'][q] < 0]
        L_ += ["", "The smallest level at which each kernel left its feature: " + '; '.join(th) + "."]
        if bl:
            L_ += ["", f"The kept kernels' coefficient growth over (|A|s)² under the paradigm's metric: median {np.median(bl):.0f} [5–95% "
                       f"{np.quantile(bl, 0.05):.0f}, {np.quantile(bl, 0.95):.0f}] ({len(bl)} kernels and levels)."]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_dead_zone.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'dead_zone_raw.npz')
    md = os.path.join(OUT, 'dead_zone.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('C'), old.get('A2'), old.get('C2'), old.get('D')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'A2', 'B', 'C', 'C2', 'D'}
    res = {k: v for k, v in old.items() if k not in parts}             # a part being rerun is not reported from the old raw
    t0 = time.time()
    for k, fn in (('A', part_A), ('A2', part_A2), ('B', part_B), ('C', part_C), ('C2', part_C2), ('D', part_D)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('C'), res.get('A2'), res.get('C2'), res.get('D')))


if __name__ == '__main__':
    main()
