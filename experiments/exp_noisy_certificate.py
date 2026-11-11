"""
Proposition D‴.56: the transit under measurement noise — the path's linear response decides the order and the support, the window and the stop
====================================================================================================================================

The sliding path (`sliding_iss`, polished) on noisy data y = f* + ξ against its first-order response along the noiseless path
(`rich_transit.path_response`): the duals at the entries u⁰ + Uξ, the entry times τ⁰ + Tξ and, after the K-th entry, the residual R_Kξ
— the projection of ξ off the K-kernel fit's tangent space —, which sets the first noise kernel's entry. The noisy paths keep their
centres in the samples' window ± BOX·σ (a relaxation that carries one out stops the path: on noisy data the fit can have no minimiser,
a kernel leaving the window fitting the noise at the edge with a coefficient growing without bound) and judge exits and stationarity on
the data's scale max|y| (referee round 45). Targets: D‴.55's recipe, three features of the kernels' width at least DMIN apart, amplitudes
±U(0.4, 1), whose certificate has m > 0.8q (so that without noise the path enters exactly the nearest kernels E).

The first-order criteria (`crit_rows`): at every entry k ≤ K the duals |u⁰_j(τ_k) + U_{k,j}ξ| of the kernels j not yet entered — 'set':
those off the path; 'ord': all of them, the path's later entrants included (a later entrant reaching its wall out of turn reorders the
entries, and the linearisation along the noiseless order no longer applies). Their deterministic radius r₁ = min (1 − |u⁰|)/|U| and the
union bound Σ 2Φ̄((1 − |u⁰|)/(s|U|)) for i.i.d. noise of variance s² per sample.

A. The response is the first order: per target, one noise direction ξ at the levels h: the noisy path's entry times and its duals at the
   K-th entry against τ⁰ + hTξ and u⁰ + hUξ — their errors ∝ h².
B. The order and the support: per target and noise level s, N_DRAW draws of i.i.d. Gaussian noise: whether the path enters exactly E first
   (as a set, and in the noiseless order), against the first-order criteria per draw, the linear model's probabilities (Monte Carlo over
   its Gaussian law) and the union bound; the linear duals against the path's own (`path_duals`) on the draws that keep the order.
B2. The same on targets of small radius r₁ (in R1_B2), where the order or the support fails at the small s where the first order holds.
C. The window: on the draws that keep the order, the first noise kernel's entry τ_{K+1} − τ_K against W = min_j (1 − sign(v_j)(u⁰_j +
   U_jξ))/|v_j|, v_j = ⟨R_Kξ, φ_j⟩, and the kernel that enters; W·s against its limit law as s → 0 (`law_W`; the probability integral
   transform of each draw's W·s under its target's law, Kolmogorov–Smirnov against the uniform).
D. Stopping and the oracle: the discrepancy rule (stop at the first entry after which the relaxed residual's rms is ≤ (1 + DELTA)s)
   against its first-order prediction (at the K-th iff rms(R_Kξ) ≤ (1 + DELTA)s; late with the probability P(χ²_{N−2K} > (1 + DELTA)²N))
   and the largest relative gap (τ_{k+1} − τ_k)/τ_k over the first K + 1 entries; the K-th fit — the K-kernel least-squares fit of y with
   free centres — against f*: its error against the projection P_Tξ of the noise on the tangent space and against s(2K/N)^{1/2}.
E. Past the first noise kernel, on four families — 'base' (N_E of B's targets), 'wide' (amplitudes WIDE_A), 'K5' (five features) and
   'long' (the horizon K + 8) —, the path to K + more entries under a budget of MAX_RHS right-hand-side evaluations: its health (runaway
   kernels, unrelaxed blocks), the largest relative gap over horizons of K + 1 … K + more entries (`gap_stop`), the discrepancy rule over
   the relaxed fits (`disc_stop`), and the fits the rules stop at against the oracle's — the K-kernel least-squares fit of y on the true
   support (`polish_block` from the truth) — on every finished path.

Usage: python experiments/exp_noisy_certificate.py [--from-raw] [--only A,B,B2,E] → results/chaos/noisy_certificate.md (+ _raw.npz).
NC_WORKERS (default 2).
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X
import exp_transit_certificate as TC
import exp_transit_mass as TM

OUT = os.environ.get('NC_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('NC_WORKERS', '2'))
SEED, N_T, DMIN = 561, 40, 0.45
N_A, H_A = 10, (1e-3, 1e-4, 1e-5)
LEVELS, N_DRAW, N_MC, DELTA, N_LAW = (0.003, 0.01, 0.03, 0.1, 0.3), 10, 20000, 0.1, 20000
BOX = 2.0                                   # the noisy paths' centres in the samples' window ± BOX·σ (referee round 45)
SEED_B2, N_B2, R1_B2, LEVELS_B2, N_DRAW_B2 = 562, 20, (0.002, 0.02), (0.001, 0.003, 0.01), 20
MAX_RHS = 3_000_000                         # E: a path's budget of right-hand-side evaluations (deterministic, unlike a clock)
FAMILIES = dict(                            # E: (targets, levels, draws per level, entries past the K-th, noise seed)
    base=dict(n=20, levels=(0.01, 0.03, 0.1, 0.3), draws=5, more=4, seed=3000),
    wide=dict(n=5, levels=(0.01, 0.03, 0.1), draws=4, more=4, seed=3100),
    K5=dict(n=4, levels=(0.01, 0.03), draws=3, more=4, seed=3200),
    long=dict(n=5, levels=(0.03, 0.1), draws=3, more=8, seed=3300))
SEED_WIDE, WIDE_A = 563, (1.0, 0.3, 0.1)    # 'wide': three features at D‴.55's positions with these amplitudes (random order, signs)
SEED_K5, K5, DMIN5 = 564, 5, 0.3            # 'K5': five features in [−0.8, 0.8] at least DMIN5 apart, amplitudes ±U(0.4, 1)


def box_of(K):
    return float(K.xs.min()) - BOX * K.sig, float(K.xs.max()) + BOX * K.sig


def targets():
    rng = np.random.default_rng(SEED)
    L, K, c0 = TM.setup()
    out = []
    while len(out) < N_T:
        (cs, sg, a, u), = TC.targets(rng, 1, DMIN)
        As = sg * a
        c = TC.certificate(K, c0, cs, As)
        if c['unique'] and c['m'] > 0.8 * c['q']:
            out.append((cs, As))
    return out


def noiseless(cs, As):
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, cs, As)
    cert = TC.certificate(K, c0, cs, As)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    return L, K, c0, f, cert, R0, RT.path_response(K, f, R0, c0)


def noisy(L, K, c0, f, xi, nK, more=1, max_rhs=None):
    """The path on y = f* + ξ up to its (K + more)-th entry, the centres in `box_of`, the thresholds on the scale max|y|, the exits by
    the rule that keeps a kernel whose removal would turn its correlation outward (`exit_kkt`: no cycles of exits and re-entries at one
    τ, referee round 48)."""
    y = f + xi
    return RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=2 * (nK + more) + 2, tol=0.0, max_entries=nK + more,
                          box=box_of(K), w_scale=float(np.abs(y).max()), max_rhs=max_rhs, exit_kkt=True)


def noise(seed, q, i, k, s, nx):
    """Draw k at level i of target q: i.i.d. N(0, s²) on the samples, its own stream (reproducible draw by draw)."""
    return s * np.random.default_rng([seed, q, i, k]).standard_normal(nx)


def crit_rows(P, d):
    """The criteria's rows (k, j), k the entry (0-based), j a kernel not yet entered there: ('set', off the path) and ('ord', all)."""
    ent = list(P['entrants'])
    rs, ro = [], []
    for k in range(len(ent)):
        for j in range(d):
            if j not in ent[:k + 1]:
                ro.append((k, j))
                if j not in ent:
                    rs.append((k, j))
    return rs, ro


def crit_arrays(P, rows):
    return np.array([P['u_all'][k][j] for k, j in rows]), np.array([P['U_all'][k][j] for k, j in rows])


def radius(u, U):
    """r₁ = min (1 − |u⁰|)/|U|: every ξ with |ξ| < r₁ keeps every row below the wall at first order."""
    return float(np.min((1 - np.abs(u)) / np.linalg.norm(U, axis=1)))


def union_bound(u, U, s):
    """Σ 2Φ̄((1 − |u⁰|)/(s|U|)) over the rows: the first-order probability of a failure is at most this."""
    from scipy.stats import norm
    return float(np.sum(2 * norm.sf((1 - np.abs(u)) / (s * np.linalg.norm(U, axis=1)))))


def law_W(P, K, c0, off, n=N_LAW, seed=0):
    """Samples of sW's limit law as s → 0: min_j (1 − sign(Z_j)u⁰_j(τ_K))/|Z_j| over the kernels off E, Z = Φ_offᵀR_KG/N with G a
    standard normal vector on the samples (Z ~ N(0, Σ⊥), Σ⊥ = Φ_offᵀR_KΦ_off/N²)."""
    nx = K.xs.size
    G = np.random.default_rng(seed).standard_normal((nx, n))
    Z = K.phi(np.asarray(c0, float)[off]).T @ (P['Rlast'] @ G) / nx
    return np.min((1 - np.sign(Z) * P['u'][off][:, None]) / np.abs(Z), axis=0)


def classify(R, order):
    """The path's entries (τ, kernel, iteration) and its exits before the (K + 1)-th entry; exact_set: its first K entries are the
    noiseless path's kernels with no exit before the (K + 1)-th; exact_order: and in the noiseless order."""
    nK = len(order)
    ev, it = R['events'], R['ev_iter']
    ent = [(e[0], e[1], i) for e, i in zip(ev, it) if e[2] == 'enter']
    pos = [p for p, e in enumerate(ev) if e[2] == 'enter']
    i1 = pos[nK] if len(pos) > nK else len(ev)
    exits = sum(e[2] == 'exit' for e in ev[:i1])
    first = [e[1] for e in ent[:nK]]
    exact_set = len(ent) >= nK and sorted(first) == sorted(order) and exits == 0
    return ent, exits, exact_set, exact_set and first == list(order)


def fit_after(ent, k):
    """The hist index of the fit relaxed after the k-th entry (1-based): the iteration of that entry plus one."""
    return ent[k - 1][2] + 1


def disc_stop(fit_rms, s, delta=DELTA):
    """The discrepancy rule: the first relaxed fit (hist index ≥ 1; fit_rms[i − 1] its rms) with rms ≤ (1 + δ)s, None if none."""
    return next((i + 1 for i, r in enumerate(fit_rms) if r <= (1 + delta) * s), None)


def gap_stop(taus, h):
    """The largest-relative-gap rule over the first h entries: the entry k* (1-based) after which (τ_{k+1} − τ_k)/τ_k is largest."""
    g = [(taus[i + 1] - taus[i]) / taus[i] for i in range(min(h, len(taus)) - 1)]
    return int(np.argmax(g)) + 1 if g else None


def a_job(args):
    q, cs, As = args
    L, K, c0, f, cert, R0, P = noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    xi = np.random.default_rng(1000 + q).standard_normal(nx)
    off = [j for j in range(c0.size) if j not in P['entrants']]
    rows = []
    for h in H_A:
        R = noisy(L, K, c0, f, h * xi, nK)
        ent = [e for e in R['events'] if e[2] == 'enter']
        t_obs = np.array([e[0] for e in ent[:nK]])
        u_obs = RT.path_duals(K, f + h * xi, R)['U'][nK]
        dt = t_obs - (P['taus'] + h * P['T'] @ xi)
        du = (u_obs - (P['u'] + h * P['U'] @ xi))[off]
        rows.append(dict(h=h, same=[e[1] for e in ent[:nK]] == P['entrants'], dt=float(np.abs(dt).max()),
                         dt1=float(h * np.abs(P['T'] @ xi).max()), du=float(np.abs(du).max()), du1=float(h * np.abs((P['U'] @ xi)[off]).max())))
    return dict(q=q, rows=rows)


def b_draw(L, K, c0, f, P, order, off, rows, s, xi):
    """One draw of B/B2: the classification, the criteria, the path's own duals and, on the exact draws, the window, the stop, the fits."""
    nK, nx = len(order), K.xs.size
    y = f + xi
    R = noisy(L, K, c0, f, xi, nK)
    ent, exits, exact_set, exact_order = classify(R, order)
    (us, Us), (uo, Uo) = rows
    d = dict(exact=exact_set, exact_order=exact_order, exits=exits, n_ent=len(ent), runaway=R['runaway'],
             crit=float(np.abs(us + Us @ xi).max()), crit_ord=float(np.abs(uo + Uo @ xi).max()),
             rms_RK=float(np.sqrt(np.mean((P['Rlast'] @ xi) ** 2))))
    if exact_set and fit_after(ent, nK) >= len(R['hist']):                 # the K-th block ran away: no K-th fit
        d.update(exact=False, exact_order=False, runaway_K=True)
        return d
    if not exact_set:
        return d
    Dp = RT.path_duals(K, y, R)['U']
    up = np.array([Dp[fit_after(ent, k)] for k in range(1, nK + 1)])        # the path's duals at its first K entries
    d['crit_path'] = float(np.abs(up[:, off]).max())
    if exact_order:
        ul = np.array([P['u_all'][k] + P['U_all'][k] @ xi for k in range(nK)])
        d['lin_err'] = float(np.abs((ul - up)[:, off]).max())
    v = K.phi(c0).T @ (P['Rlast'] @ xi) / nx
    u = P['u'] + P['U'] @ xi
    W = {j: (1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off if v[j] != 0}
    jn = min(W, key=W.get)
    Ws = sorted(W.values())
    taus = [e[0] for e in ent]
    d.update(tauK=taus[nK - 1], W_pred=float(W[jn]), j_pred=int(jn), W_gap=float((Ws[1] - Ws[0]) / abs(Ws[0])),
             W_obs=float(taus[nK] - taus[nK - 1]) if len(ent) > nK else np.nan, j_obs=int(ent[nK][1]) if len(ent) > nK else -1,
             gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)])
    fits = [R['hist'][fit_after(ent, k)] for k in range(1, nK + 1)]            # the relaxed fits after entries 1 … K
    rms = [h[3] for h in fits]
    d['k_disc'] = next((i + 1 for i, r_ in enumerate(rms) if r_ <= (1 + DELTA) * s), -1)
    d['rms_K'] = float(rms[-1])
    d['relaxed_K'] = bool(all(R['relaxed'][:ent[nK - 1][2] + 1]))
    fhat = K.phi(fits[-1][1]) @ fits[-1][2]
    d['err'] = float(np.sqrt(np.mean((fhat - f) ** 2)))
    d['err_T'] = float(np.sqrt(np.mean(((np.eye(nx) - P['Rlast']) @ xi) ** 2)))
    return d


def b_target(tag, q, cs, As, levels, n_draw, seed):
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = noiseless(cs, As)
    nx, d_ = K.xs.size, c0.size
    order = list(P['entrants'])
    off = [j for j in range(d_) if j not in order]
    rs, ro = crit_rows(P, d_)
    rows = (crit_arrays(P, rs), crit_arrays(P, ro))
    (us, Us), (uo, Uo) = rows
    law = np.sort(law_W(P, K, c0, off, seed=seed + q))
    taus0 = P['taus']
    out = dict(tag=tag, q=q, m=cert['m'], q_ct=cert['q'], path_exact=sorted(order) == sorted(cert['e']), r1=radius(us, Us),
               r1_ord=radius(uo, Uo), rel_gap0=float(np.min(np.diff(taus0) / taus0[:-1])), law=law[::N_LAW // 1000], levels=[])
    for i, s in enumerate(levels):
        Z = s * np.random.default_rng([seed + 1, q, i]).standard_normal((nx, N_MC))
        lv = dict(s=s, p_mc=float(np.mean(np.abs(us[:, None] + Us @ Z).max(0) < 1)),
                  p_ord=float(np.mean(np.abs(uo[:, None] + Uo @ Z).max(0) < 1)),
                  ub=union_bound(us, Us, s), ub_ord=union_bound(uo, Uo, s), draws=[])
        for k in range(n_draw):
            d = b_draw(L, K, c0, f, P, order, off, rows, s, noise(seed, q, i, k, s, nx))
            if d['exact_order'] and np.isfinite(d.get('W_obs', np.nan)):
                d['pit'] = float(np.searchsorted(law, d['W_obs'] * s) / law.size)
            lv['draws'].append(d)
        out['levels'].append(lv)
    out['secs'] = time.time() - t0
    print(f"  {tag} target {q}: r1 {out['r1_ord']:.4f}, exact " +
          ', '.join(f"{sum(d['exact'] for d in lv['draws'])}/{sum(d['exact_order'] for d in lv['draws'])}" for lv in out['levels']) +
          f" [{out['secs']:.0f}s]", flush=True)
    return out


def b_job(args):
    q, cs, As = args
    return b_target('B', q, cs, As, LEVELS, N_DRAW, 2000)


def b2_job(args):
    q, cs, As = args
    return b_target('B2', q, cs, As, LEVELS_B2, N_DRAW_B2, 2500)


def b2_candidate(args):
    """B2's candidate: (cs, As) → its radius r₁ (the order criterion's) if the noiseless path enters three kernels without exits at
    distinct times, else None."""
    cs, As = args
    warnings.filterwarnings('ignore')
    try:
        L, K, c0, f, cert, R0, P = noiseless(cs, As)
    except ValueError:
        return None
    if len(R0['events']) != len(cs) or R0['hist'][-1][3] > 1e-8:
        return None
    return radius(*crit_arrays(P, crit_rows(P, c0.size)[1]))


def targets_B2(pool=None, batch=40):
    """B2: D‴.55's recipe (no certificate condition) for targets whose noiseless path is exact with the radius r₁ in R1_B2 — the order
    or the support fails at noise levels where the first order holds. The candidates in the stream's order, evaluated in batches."""
    rng = np.random.default_rng(SEED_B2)
    out = []
    while len(out) < N_B2:
        cand = []
        for _ in range(batch):
            (cs, sg, a, u), = TC.targets(rng, 1, DMIN)
            cand.append((cs, sg * a))
        rr = pool.map(b2_candidate, cand, chunksize=1) if pool is not None else [b2_candidate(c) for c in cand]
        out += [c for c, r in zip(cand, rr) if r is not None and R1_B2[0] <= r <= R1_B2[1]]
    return out[:N_B2]


def path_target(cs, As):
    """A noiseless path without the certificate: (L, K, c0, f, the entrants in order) if it enters len(cs) kernels without exits and
    fits exactly, else None."""
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, cs, As)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    if any(e[2] != 'enter' for e in R0['events']) or len(R0['events']) != len(cs) or R0['hist'][-1][3] > 1e-8:
        return None
    return L, K, c0, f, [e[1] for e in R0['events']]


def targets_wide():
    rng = np.random.default_rng(SEED_WIDE)
    out = []
    while len(out) < FAMILIES['wide']['n']:
        (cs, sg, a, u), = TC.targets(rng, 1, DMIN)
        As = sg * rng.permutation(WIDE_A)
        if path_target(cs, As) is not None:
            out.append((cs, As))
    return out


def targets_K5():
    rng = np.random.default_rng(SEED_K5)
    out = []
    while len(out) < FAMILIES['K5']['n']:
        cs = np.sort(rng.uniform(-0.8, 0.8, K5))
        if np.diff(cs).min() < DMIN5:
            continue
        As = rng.choice([-1.0, 1.0], K5) * rng.uniform(0.4, 1.0, K5)
        if path_target(cs, As) is not None:
            out.append((cs, As))
    return out


def e_targets(fam):
    return dict(base=lambda: targets()[:FAMILIES['base']['n']], wide=targets_wide, K5=targets_K5,
                long=lambda: targets()[:FAMILIES['long']['n']])[fam]()


def e_path(L, K, c0, f, cs, As, order, s, xi, more):
    """E: one path to its (K + more)-th entry under MAX_RHS: its classification, entries and fits, their health, and the oracle's error."""
    nK, nx = len(order), K.xs.size
    lo, hi = float(K.xs.min()), float(K.xs.max())
    y = f + xi
    R = noisy(L, K, c0, f, xi, nK, more, max_rhs=MAX_RHS)
    ent, exits, exact_set, exact_order = classify(R, order)
    hs = R['hist'][1:]
    out_ = lambda h: float(max([0.0] + [max(c - hi, lo - c) for c, w in zip(h[1], h[2]) if w != 0]))
    wo, co = RT.polish_block(K, y, np.asarray(As, float), np.asarray(cs, float), np.sign(As), cbox=box_of(K))
    return dict(s=s, exact=exact_set, exact_order=exact_order, exits_before=exits, exits=sum(e[2] == 'exit' for e in R['events']),
                runaway=R['runaway'], budget=R['out_of_budget'], nfev=R['nfev'], stiff=R['stiff_redone'], n_ent=len(ent),
                taus=[e[0] for e in ent], ent_j=[e[1] for e in ent], ent_iter=[e[2] for e in ent], relaxed=list(R['relaxed']),
                fit_rms=[float(h[3]) for h in hs], fit_err=[float(np.sqrt(np.mean((K.phi(h[1]) @ h[2] - f) ** 2))) for h in hs],
                fit_out=[out_(h) for h in hs], fit_maxw=[float(np.abs(h[2]).max()) if np.any(h[2]) else 0.0 for h in hs],
                oracle=float(np.sqrt(np.mean((K.phi(co) @ wo - f) ** 2))), x_end=float(np.abs(R['c'][R['w'] != 0]).max()),
                maxw_end=float(np.abs(R['w']).max()))


def e_job(args):
    fam, q, cs, As = args
    warnings.filterwarnings('ignore')
    F = FAMILIES[fam]
    t0 = time.time()
    L, K, c0, f, order = path_target(cs, As)
    nx = K.xs.size
    out = dict(fam=fam, q=q, K=len(cs), more=F['more'], order=order, amps=list(np.abs(As)), levels=[])
    for i, s in enumerate(F['levels']):
        out['levels'].append(dict(s=s, draws=[e_path(L, K, c0, f, cs, As, order, s, noise(F['seed'], q, i, k, s, nx), F['more'])
                                              for k in range(F['draws'])]))
    out['secs'] = time.time() - t0
    print(f"  E {fam} target {q}: " + ', '.join(f"{sum(d['exact'] for d in lv['draws'])}/{len(lv['draws'])}"
                                               f" (runaway {sum(d['runaway'] for d in lv['draws'])})" for lv in out['levels']) +
          f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_A():
    from multiprocessing import Pool
    tg = targets()[:N_A]
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def part_B():
    from multiprocessing import Pool
    tg = targets()
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def part_B2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        tg = targets_B2(pool)
        return pool.map(b2_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def part_E():
    from multiprocessing import Pool
    jobs = [(fam, q, cs, As) for fam in FAMILIES for q, (cs, As) in enumerate(e_targets(fam))]
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, jobs, chunksize=1)


def _z(obs, p):
    """The observed count against Σp (independent Bernoulli draws): (expected, sd, z)."""
    p = np.asarray(p, float)
    e, sd = p.sum(), np.sqrt(np.sum(p * (1 - p)))
    return e, sd, (obs - e) / sd if sd > 0 else 0.0


def _b_tables(B, levels, n_draw, title):
    L_ = [f"## {title}", ""]
    r1 = np.array([b['r1_ord'] for b in B])
    L_ += [f"{len(B)} targets ({sum(b['path_exact'] for b in B)} with the noiseless path exactly E), {n_draw} draws per level. The "
           f"radius r₁ = min (1 − |u⁰|)/|U| over the order criterion's rows: {r1.min():.4f} … {r1.max():.4f} (median {np.median(r1):.4f}); "
           f"over the support criterion's alone: {min(b['r1'] for b in B):.4f} … {max(b['r1'] for b in B):.4f}; the targets whose consecutive "
           f"entries lie within 5% (τ_(k+1)/τ_k < 1.05): {sum(b['rel_gap0'] < 0.05 for b in B)}; draws whose K-th block ran away (not "
           f"exact): {sum(d.get('runaway_K', False) for b in B for lv in b['levels'] for d in lv['draws'])}.", "",
           "Observed: the path enters exactly E first (as a set; in the noiseless order — the swaps the difference); predicted per draw: the "
           "order criterion max|u⁰_j + U_jξ| < 1 over the entries and the kernels not yet entered (E's later entrants included), and the "
           "support criterion over the kernels off E alone; the linear model's expected count of the order kept (Monte Carlo) ± its sd and "
           "the observed count's z; the union bound on the expected failures (Σ over targets of min(1, Σ_rows 2Φ̄)·draws) against the "
           "observed; disagreements of the order criterion, and of these the ones within the level's largest linear-dual error of 1.", "",
           "| s | exact: set, order (swaps) | criterion < 1: order, support | agree per draw: order, support | expected order (linear "
           "model) ± sd, z | union bound: order [support] vs failures | disagreements (within the linear error of 1) |",
           "|---|---|---|---|---|---|---|"]
    for i, s in enumerate(levels):
        ds = [d for b in B for d in b['levels'][i]['draws']]
        es = np.array([d['exact'] for d in ds]); eo = np.array([d['exact_order'] for d in ds])
        co = np.array([d['crit_ord'] < 1 for d in ds]); cs_ = np.array([d['crit'] < 1 for d in ds])
        e_, sd, z = _z(eo.sum(), [b['levels'][i]['p_ord'] for b in B for _ in range(n_draw)])
        ub = sum(min(1.0, b['levels'][i]['ub_ord']) * n_draw for b in B)
        ubs = sum(min(1.0, b['levels'][i]['ub']) * n_draw for b in B)
        le = max([d['lin_err'] for d in ds if 'lin_err' in d] or [0.0])
        dis = [d for d in ds if d['exact_order'] != (d['crit_ord'] < 1)]
        L_.append(f"| {s:g} | {es.sum()}, {eo.sum()} ({es.sum() - eo.sum()}) | {co.sum()}, {cs_.sum()} | {np.sum(eo == co)}, "
                  f"{np.sum(es == cs_)} | {e_:.1f} ± {sd:.1f}, {z:+.1f} | {ub:.1f} [{ubs:.1f}] vs {len(ds) - eo.sum()} [{len(ds) - es.sum()}] | "
                  f"{len(dis)} ({sum(abs(d['crit_ord'] - 1) <= le for d in dis)}) |")
    L_ += ["", "The linear duals against the path's own (`path_duals`) at its first K entries, the kernels off E: their largest difference "
           "on the draws that keep the order (and over s²); on the swaps, the support criterion against the path's own largest dual off E.", "",
           "| s | order kept: max error | / s² | swaps: max |crit − path's| |", "|---|---|---|---|"]
    for i, s in enumerate(levels):
        ds = [d for b in B for d in b['levels'][i]['draws']]
        le = [d['lin_err'] for d in ds if 'lin_err' in d]
        sw = [abs(d['crit'] - d['crit_path']) for d in ds if d['exact'] and not d['exact_order']]
        L_.append(f"| {s:g} | {max(le):.2e} | {max(le) / s ** 2:.2f} | " + (f"{max(sw):.3f} |" if sw else "— |") if le else f"| {s:g} | — | — | — |")
    return L_


def report(A=None, B=None, B2=None, E_=None):
    from scipy.stats import kstest, binomtest, chi2
    L_ = ['# Proposition D‴.56: the transit under measurement noise (`exp_noisy_certificate.py`)', '',
          "The sliding path on y = f* + ξ (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]) against its "
          f"first-order response along the noiseless path (`path_response`); targets: three features at least {DMIN} apart (D‴.55's "
          f"recipe), certificate margin m > 0.8q. The noisy paths keep their centres in the samples' window ± {BOX:g}σ (a relaxation that "
          "carries one out stops the path: runaway) and judge exits and stationarity on the data's scale max|y|; a kernel exits only if "
          "its correlation stays inward with its own contribution removed.", '']
    if A is not None:
        L_ += ['## A. The response is the first order', '',
               "Per target (one noise direction ξ, unit variance per sample) and level h: the noisy path's entry times and its duals at the "
               "K-th entry against the first-order prediction — the largest error, and the first-order term's own size.", '',
               '| target | h | same entries | entry times: error (first-order term) | duals: error (first-order term) |', '|---|---|---|---|---|']
        for a in A:
            for r in a['rows']:
                L_.append(f"| {a['q']} | {r['h']:g} | {r['same']} | {r['dt']:.2e} ({r['dt1']:.2e}) | {r['du']:.2e} ({r['du1']:.2e}) |")
        rat = [a['rows'][i]['du'] / a['rows'][i + 1]['du'] for a in A for i in range(len(a['rows']) - 1) if a['rows'][i + 1]['du'] > 0]
        L_ += ['', f"The duals' error falls by {min(rat):.1f}–{max(rat):.1f} per decade of h (second order: 100).", '']
    if B is not None:
        L_ += _b_tables(B, LEVELS, N_DRAW, 'B. The order and the support')
        if B2 is not None:
            L_ += [''] + _b_tables(B2, LEVELS_B2, N_DRAW_B2, f'B2. Small radius (r₁ in [{R1_B2[0]:g}, {R1_B2[1]:g}])')
        pooled = np.concatenate([b['law'] for b in B])
        pq = np.percentile(pooled, [25, 50, 75])
        L_ += ['', '## C. The window', '',
               "On the draws that keep the order with the order criterion below 1 (the first order's regime): the first noise kernel's entry "
               "after the K-th, τ_{K+1} − τ_K, against the first-order W — "
               "the ratio's median, range and 5–95% quantiles, and their half-width over s (constant if the relative error is O(s)) —; the "
               "kernel as predicted, and on the draws where it is not, the relative gap between the two smallest predicted windows "
               "(W₍₂₎ − W₍₁₎)/W₍₁₎ (median; all draws' median in brackets); W·s (median, quartiles) against its limit law (`law_W`, "
               f"{N_LAW} samples per target; pooled over the targets: median {pq[1]:.3g}, quartiles [{pq[0]:.3g}, {pq[2]:.3g}]) — the "
               "Kolmogorov–Smirnov p-value of the draws' probability integral transforms under their targets' laws —; outside the regime: "
               "the draws that keep the order with the criterion ≥ 1, and those with the support exact but the order swapped (of each, "
               "those with W_pred < 0: the first-order dual already past its wall).", '',
               '| s | draws | W_obs/W_pred: median [min, max] | 5–95% | half-width / s | next kernel as predicted | '
               'missed: (W₍₂₎ − W₍₁₎)/W₍₁₎ [all] | W·s: median [quartiles] | KS p | order kept, criterion ≥ 1 (W_pred < 0) | '
               'swaps (W_pred < 0) |', '|---|---|---|---|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            al = [d for b in B for d in b['levels'][i]['draws']]
            ds = [d for d in al if d['exact_order'] and d['crit_ord'] < 1 and np.isfinite(d['W_obs'])]
            out = [d for d in al if d['exact_order'] and d['crit_ord'] >= 1]
            sw = [d for d in al if d['exact'] and not d['exact_order']]
            oo = f"{len(out)} ({sum(d['W_pred'] < 0 for d in out)}) | {len(sw)} ({sum(d['W_pred'] < 0 for d in sw)})"
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — | — | — | — | {oo} |")
                continue
            r = np.array([d['W_obs'] / d['W_pred'] for d in ds])
            q5, q95 = np.percentile(r, [5, 95])
            ws = np.array([d['W_obs'] * s for d in ds])
            miss = [d['W_gap'] for d in ds if d['j_obs'] != d['j_pred']]
            mg = f"{np.median(miss):.3g}" if miss else '—'
            p_ks = kstest([d['pit'] for d in ds], 'uniform').pvalue
            L_.append(f"| {s:g} | {len(ds)} | {np.median(r):.4f} [{r.min():.3f}, {r.max():.3f}] | [{q5:.4f}, {q95:.4f}] | "
                      f"{(q95 - q5) / 2 / s:.2f} | {sum(d['j_obs'] == d['j_pred'] for d in ds)}/{len(ds)} | {mg} "
                      f"[{np.median([d['W_gap'] for d in ds]):.3g}] | {np.median(ws):.3g} [{np.percentile(ws, 25):.3g}, "
                      f"{np.percentile(ws, 75):.3g}] | {p_ks:.2g} | {oo} |")
        nK, nx = 3, 61
        p_late = float(chi2.sf((1 + DELTA) ** 2 * nx, nx - 2 * nK))
        L_ += ['', '## D. Stopping and the oracle', '',
               f"The discrepancy rule: the first entry after which the relaxed residual's rms is ≤ {1 + DELTA:g}s, on the draws with the "
               f"support exact: at the K-th, not by the K-th (late), earlier; its first-order prediction — at the K-th iff rms(R_Kξ) ≤ "
               f"{1 + DELTA:g}s, R_Kξ the first-order residual (|R_Kξ|²/s² ~ χ²_(N−2K)) — and its agreement per draw; the path's late stops "
               f"against the first order's binomial law, P(χ²_(N−2K) > {1 + DELTA:g}²N) = {p_late:.4f} per draw (expected ± sd, the binomial "
               "test's p); the relaxed residual after the K-th entry against rms(R_Kξ).", '',
               '| s | draws | at K | late | early | first order: at K | agree per draw | late: expected ± sd, p | '
               'rms_K / rms(R_Kξ): median [5–95%] |', '|---|---|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact']]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — | — | — | — |")
                continue
            kd = np.array([d['k_disc'] for d in ds])
            fo = np.array([d['rms_RK'] <= (1 + DELTA) * s for d in ds])
            rr = np.array([d['rms_K'] / d['rms_RK'] for d in ds])
            late = int(np.sum((kd > nK) | (kd < 0)))
            L_.append(f"| {s:g} | {len(ds)} | {int(np.sum(kd == nK))} | {late} | {int(np.sum((kd > 0) & (kd < nK)))} | {int(fo.sum())} | "
                      f"{int(np.sum(fo == (kd == nK)))} | {p_late * len(ds):.1f} ± {np.sqrt(len(ds) * p_late * (1 - p_late)):.1f}, "
                      f"{binomtest(late, len(ds), p_late).pvalue:.2g} | {np.median(rr):.4f} [{np.percentile(rr, 5):.4f}, "
                      f"{np.percentile(rr, 95):.4f}] |")
        L_ += ['', "The largest relative gap (τ_{k+1} − τ_k)/τ_k over the first K + 1 entries falling after the K-th; the K-th entry's fit (the "
               "K-kernel least-squares fit of y): its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against "
               "s(2K/N)^{1/2}.", '',
               '| s | draws | largest gap at K | error / rms(P_Tξ): median [range] | 5–95% | error / (s(2K/N)^{1/2}): median |',
               '|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact']]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — |")
                continue
            e1 = np.array([d['err'] / d['err_T'] for d in ds])
            e2 = np.array([d['err'] / (s * np.sqrt(2 * nK / nx)) for d in ds])
            g1 = [int(np.argmax(d['gaps'][:nK])) == nK - 1 for d in ds if len(d['gaps']) >= nK]
            L_.append(f"| {s:g} | {len(ds)} | {sum(g1)}/{len(g1)} | {np.median(e1):.4f} [{e1.min():.3f}, {e1.max():.3f}] | "
                      f"[{np.percentile(e1, 5):.4f}, {np.percentile(e1, 95):.4f}] | {np.median(e2):.3f} |")
    if E_ is not None:
        L_ += ['', '## E. Past the first noise kernel', '',
               "Four families: 'base' (B's first targets), 'wide' (amplitudes " + ', '.join(f'{a:g}' for a in WIDE_A) + " in random order "
               "and signs), 'K5' (five features), 'long' (B's first targets, the horizon K + 8); each path to its (K + more)-th entry under "
               f"{MAX_RHS:.0e} right-hand-side evaluations. Paths: finished (out of the budget), stopped by a runaway kernel before the "
               "horizon; the first K entries E as a set; the largest relative gap over the first K + 1, K + 4 (and K + 8) entries at the K-th "
               "(on the exact paths); its condition, W/τ_K above the features' own largest relative gap and above the later ones (on the "
               "exact paths that reach the horizon); the discrepancy rule over the relaxed fits (at the K-th fit, later, earlier, never); the "
               "fits the two rules stop at — the gap rule over the longest horizon — against the oracle (the K-kernel least-squares fit of y "
               "on the true support) on every finished path; the stopped fits unrelaxed or with a centre outside the samples' window. "
               "Later exits: paths with an exit after the (K + 1)-th entry.", '',
               '| family | s | paths (budget) | runaway, later exits | exact | largest gap at K: K+1, K+4[, K+8] | W/τ_K > features\' '
               'gaps, > later gaps | discrepancy: at K, later, early, never | gap-stopped error / oracle: median [max] | discrepancy-stopped: '
               'median [max] | stopped fits unrelaxed or outside: gap, discrepancy |', '|---|---|---|---|---|---|---|---|---|---|---|']
        for fam in dict.fromkeys(e['fam'] for e in E_):
            es = [e for e in E_ if e['fam'] == fam]
            if not es:
                continue
            nK, more = es[0]['K'], es[0]['more']
            hz = [nK + 1, nK + 4] + ([nK + 8] if more >= 8 else [])
            for i, s in enumerate(lv['s'] for lv in es[0]['levels']):
                al = [d for e in es for d in e['levels'][i]['draws']]
                fin = [d for d in al if not d['budget']]
                ex = [d for d in fin if d['exact']]
                g_at = [sum(gap_stop(d['taus'], h) == nK for d in ex if d['n_ent'] >= h) for h in hz]
                g_n = [sum(d['n_ent'] >= h for d in ex) for h in hz]
                full = [d for d in ex if d['n_ent'] >= nK + more]
                gp = lambda d: [(d['taus'][k + 1] - d['taus'][k]) / d['taus'][k] for k in range(len(d['taus']) - 1)]
                c1 = sum(gp(d)[nK - 1] > max(gp(d)[:nK - 1]) for d in full)
                c2 = sum(gp(d)[nK - 1] > max(gp(d)[nK:]) for d in full)
                kd = []
                rg, rd, bad_d, bad_g = [], [], 0, 0
                for d in fin:
                    iK = d['ent_iter'][nK - 1] + 1 if d['n_ent'] >= nK else None
                    j = disc_stop(d['fit_rms'], s)
                    kd.append('never' if j is None else 'K' if j == iK else 'later' if iK is not None and j > iK else 'early')
                    if j is not None:
                        rd.append(d['fit_err'][j - 1] / d['oracle'])
                        bad_d += (not d['relaxed'][j - 1]) or d['fit_out'][j - 1] > 0
                    k = gap_stop(d['taus'], len(d['taus']))
                    if k is not None and d['ent_iter'][k - 1] + 1 <= len(d['fit_err']):
                        jg = d['ent_iter'][k - 1] + 1
                        rg.append(d['fit_err'][jg - 1] / d['oracle'])
                        bad_g += (not d['relaxed'][jg - 1]) or d['fit_out'][jg - 1] > 0
                fm = lambda x: f"{np.median(x):.3f} [{np.max(x):.3g}]" if x else '—'
                L_.append(f"| {fam} | {s:g} | {len(fin)} ({len(al) - len(fin)}) | {sum(d['runaway'] for d in fin)}, "
                          f"{sum(d['exits'] > d['exits_before'] for d in fin)} | {len(ex)} | "
                          + ', '.join(f"{a}/{n}" for a, n in zip(g_at, g_n)) + f" | {c1}, {c2} of {len(full)} | "
                          f"{kd.count('K')}, {kd.count('later')}, {kd.count('early')}, {kd.count('never')} | {fm(rg)} | {fm(rd)} | {bad_g}, {bad_d} |")
        L_ += ['', "The features' own relative gaps (τ_{k+1} − τ_k)/τ_k, k < K, and the later ones (k > K), on the exact paths that reach "
               "the horizon: their largest (median [max]), against W/τ_K (median [min]).", '',
               "| family | s | features' largest | later largest | W/τ_K |", "|---|---|---|---|---|"]
        for fam in dict.fromkeys(e['fam'] for e in E_):
            es = [e for e in E_ if e['fam'] == fam]
            if not es:
                continue
            nK, more = es[0]['K'], es[0]['more']
            for i, s in enumerate(lv['s'] for lv in es[0]['levels']):
                full = [d for e in es for d in e['levels'][i]['draws'] if not d['budget'] and d['exact'] and d['n_ent'] >= nK + more]
                if not full:
                    L_.append(f"| {fam} | {s:g} | — | — | — |")
                    continue
                gp = [[(d['taus'][k + 1] - d['taus'][k]) / d['taus'][k] for k in range(len(d['taus']) - 1)] for d in full]
                a = [max(g[:nK - 1]) for g in gp]; b = [max(g[nK:]) for g in gp]; w = [g[nK - 1] for g in gp]
                L_.append(f"| {fam} | {s:g} | {np.median(a):.3g} [{max(a):.3g}] | {np.median(b):.3g} [{max(b):.3g}] | "
                          f"{np.median(w):.3g} [{min(w):.3g}] |")
    L_ += ['', 'Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'noisy_certificate_raw.npz')
    md = os.path.join(OUT, 'noisy_certificate.md')
    old = E.load(raw)
    keys = ('A', 'B', 'B2', 'E')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('B2', part_B2), ('E', part_E)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
