"""
Proposition D‴.48: the rich transit — from small masses the kernels enter one at a time
====================================================================================================================================

The free mass–amplitude law (Proposition D‴.45; `mirror_free`, `rich_transit`) from small masses, in the time rescaled by the transit's
scale.

A. From α = 0 (masses εb_i, three 6 × 14 Gaussian designs, unequal b_i): the inverse scale space path of Σ|w_i|/√b_i (`iss_path`) against
   the loop in τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{1/2}, u = η/Λ, for ε = 10⁻² … 10⁻⁸: sup_τ|u − u_ISS| (with its log-log slope),
   the residual's rms difference over τ, the entry times, the end's distance to weighted basis pursuit.
B. On branches (amplitudes s_i·U(0.2, 0.6), masses ε): the one-sided path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in
   τ = t/log(1/ε), u = (η − η(0))/log(1/ε), for ε = 10⁻⁴ … 10⁻⁶⁴; a ± pair per atom (amplitudes ±0.4, masses with w(0) = u) against the
   two-sided path (weights 1/(γa₀)).
C. Moving centres (Gaussian kernels of width 0.12 on 61 samples of [−1, 1]; 21 kernels on a grid; the target three bumps of the same
   width off the grid): (1) the paradigm's law (ċ = Gβa⟨r∂_cφ⟩) against the sliding inverse scale space (`sliding_iss`) for ε = 10⁻³ …
   10⁻⁷ — entries, the formed kernels, the inactive kernels' displacement and mass; (2) the Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩):
   the inactive kernels' centres against their limit (`wfr_inactive_rhs`) before the first entry — condensation —, the formed kernels;
   (3) the centres held fixed (the grid dictionary's transit); (4) random three-bump targets: the formed kernels under each metric and with
   fixed centres; (5) two bumps closer than the kernels' width, and a target of wider bumps than the kernels; (6) a lost near-tie's creep.
D. Moving widths too (the paradigm's full kernel law, `moving_rhs_widths`): targets of bumps wider (0.2) and narrower (0.08) than the
   kernels' starting width 0.12, ε = 10⁻⁵, 10⁻⁶, against fixed widths.

Output: results/chaos/mirror_transit.md and mirror_transit_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`).
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
import exp_cohesion_leak as E

OUT = os.environ.get('MT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                       # c = 0.025
NA, DA, SEEDS_A = 6, 14, (3, 4, 5)
EPS_A = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 1e-12)
NT_A = 2001                                                           # τ-grid points on [0, 1.5τ*]
EPS_B, U_B = (1e-4, 1e-8, 1e-16, 1e-32, 1e-64), (1e-4, 1e-8, 1e-16, 1e-32, 1e-64)
NS_C, SIG_C, D_C = 61, 0.12, 21
CSTAR, ASTAR = (-0.52, 0.07, 0.61), (0.9, -0.6, 0.75)
EPS_C1, EPS_C2 = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), (1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
TAU_PRE = 10.0                                                        # before the first entry (C2)
N_RAND, EPS_RAND, EPS_RAND2, EPS_RAND3 = 16, 1e-5, 1e-6, 1e-8                # C4: the paradigm at all three, fixed centres at the second
EPS_WFR_RAND = (1e-6, 1e-8, 1e-10)                                          # C4: the WFR metric's count as ε → 0 (round 37)
EPS_C2B = 1e-6                                                              # C2: the WFR limit with unequal masses
EPS_C5 = (1e-6, 1e-8, 1e-10)                                                # C5: ties, near ties, close off-grid bumps
CREEP_TARGET, EPS_CREEP = 10, (1e-8, 1e-10)                                  # C6: a lost near-tie's creep
TH_ACT = 1e-2                                                         # a formed kernel: |w| > TH_ACT·max|A*|
WIDTHS_D, EPS_D, TAU_D = (0.2, 0.08), (1e-5, 1e-6), 60.0                    # D: the targets' widths, the masses, the horizon


def law():
    return MF.FreeLaw(**LAW)


def design(seed):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((NA, DA))
    ws = np.zeros(DA); idx = rng.choice(DA, 3, replace=False); ws[idx] = rng.choice([-1, 1], 3) * rng.uniform(0.3, 0.9, 3)
    b = np.exp(rng.uniform(np.log(0.25), np.log(4.0), DA))
    return Phi, ws, b


def entries_of(W, taus, th):
    """The first τ at which |w_i| exceeds th, per coordinate (inf if never)."""
    out = np.full(W.shape[1], np.inf)
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            out[i] = taus[k[0]]
    return out


def part_A():
    """The inverse scale space limit from α = 0."""
    warnings.filterwarnings('ignore')
    L = law()
    res = []
    for seed in SEEDS_A:
        Phi, ws, b = design(seed); y = Phi @ ws
        rho = b ** -0.5
        path = MF.iss_path(Phi, y, rho)
        wbp = MF.weighted_bp(Phi, y, rho)
        tau_end = 1.5 * path['taus'][-1]
        taus = np.linspace(0, tau_end, NT_A)
        Uref = np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])
        Wref = np.array([MF.iss_eval(path, t, Phi, y)[1] for t in taus])
        th = 1e-3 * np.abs(wbp).max()
        ent_ref = entries_of(Wref, taus, th)                                    # the path's own crossings of the threshold
        rows = []
        for eps in EPS_A:
            t0 = time.time()
            Lam = RT.transit_scale(L, eps)
            b0 = eps * b; q = L.Q(np.zeros(DA), b0)
            sol = solve_ivp(L.rhs_curve(Phi, y, b0), (0, Lam * tau_end), np.zeros(DA), method='LSODA', rtol=1e-10, atol=1e-13,
                            t_eval=Lam * taus)                                  # α alone, β on its curve: Q kept exactly (round 37)
            al = sol.y.T
            U = np.array([L.eta(al[k], q) for k in range(len(taus))]) / Lam
            W = np.array([L.w(al[k], q) for k in range(len(taus))])
            ent = entries_of(W, taus, th)
            both = np.isfinite(ent) & np.isfinite(ent_ref)
            dr = np.sqrt(np.sum((W @ Phi.T - Wref @ Phi.T) ** 2, axis=1))         # |Φ(w − w_ISS)| at each sample
            rows.append(dict(eps=eps, Lam=float(Lam), du=float(np.abs(U - Uref).max()),
                             dres=float(np.sqrt(np.mean(dr ** 2))), dres_med=float(np.median(dr)), dres_n=int(np.sum(dr > 1e-3)),
                             dent=float(np.abs(ent - ent_ref)[both].max()) if both.any() else np.nan, n_only=int(np.sum(np.isfinite(ent) != np.isfinite(ent_ref))),
                             dend=float(np.abs(W[-1] - wbp).max()), secs=time.time() - t0))
        res.append(dict(seed=seed, events=[(float(t), int(i), k) for t, i, k in path['events']], tau_end=float(tau_end),
                        end_is_bp=float(np.abs(path['w_end'] - wbp).max()), rows=rows, dtau=float(taus[1] - taus[0])))
    return res


def part_B():
    """Branches (one-sided path, the scale log(1/ε)) and ± pairs (two-sided)."""
    warnings.filterwarnings('ignore')
    L = law()
    Phi, ws, b = design(SEEDS_A[0]); y = Phi @ ws
    rng = np.random.default_rng(11)
    s = rng.choice([-1.0, 1.0], DA); s[ws != 0] = np.sign(ws[ws != 0])
    amp = rng.uniform(0.2, 0.6, DA)
    a_init = L.A * np.tanh(amp / L.A)
    rho = 1 / (L.gamma * a_init)
    path = MF.iss_path(Phi, y, rho, signs=s)
    wlim = MF.weighted_bp(Phi, y, rho, signs=s)
    out = dict(events=[(float(t), int(i), k) for t, i, k in path['events']], end_is_lim=float(np.abs(path['w_end'] - wlim).max()))
    tau_end = 1.5 * path['taus'][-1]
    taus = np.linspace(0, tau_end, NT_A)
    Uref = np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])
    Wref = np.array([MF.iss_eval(path, t, Phi, y)[1] for t in taus])
    th = 1e-3 * np.abs(wlim).max()
    ent_ref = entries_of(Wref, taus, th)
    rows = []
    for eps in EPS_B:
        t0 = time.time()
        Lg = np.log(1 / eps)
        P = MF.ParticleSet(L, s * amp, np.full(DA, eps))
        e0 = P.eta(P.x0)
        sol = solve_ivp(P.rhs(Phi, y), (0, Lg * tau_end), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=Lg * taus,
                        max_step=Lg * tau_end / NT_A)                    # a capped step: no trial state beyond the curves' range
        U = np.array([(P.eta(sol.y[:, k]) - e0) / Lg for k in range(len(taus))])
        W = np.array([P.w(sol.y[:, k]) for k in range(len(taus))])
        ent = entries_of(W, taus, th)
        both = np.isfinite(ent) & np.isfinite(ent_ref)
        rows.append(dict(eps=eps, all_branch=bool(P.br.all()), du=float(np.abs(U - Uref).max()),
                         dres=float(np.sqrt(np.mean(np.sum((W @ Phi.T - Wref @ Phi.T) ** 2, axis=1)))),
                         dent=float(np.abs(ent - ent_ref)[both].max()) if both.any() else np.nan, n_only=int(np.sum(np.isfinite(ent) != np.isfinite(ent_ref))),
                         dend=float(np.abs(W[-1] - wlim).max()), secs=time.time() - t0))
    out['rows'] = rows
    # a ± pair per atom: the two-sided path with the weights 1/(γa₀)
    al = 0.4
    Phi2 = np.hstack([Phi, Phi])
    rows2 = []
    for u0 in U_B:
        t0 = time.time()
        bb = u0 * L.m / (L.A * np.tanh(al / L.A))
        P = MF.ParticleSet(L, np.r_[np.full(DA, al), np.full(DA, -al)], np.full(2 * DA, bb))
        a0 = L.A * np.tanh(P.al0[0] / L.A)
        rho2 = np.full(DA, 1 / (L.gamma * a0))
        path2 = MF.iss_path(Phi, y, rho2)
        te = 1.5 * path2['taus'][-1]
        tg = np.linspace(0, te, NT_A)
        Lg = np.log(1 / u0)
        e0 = P.eta(P.x0)
        sol = solve_ivp(P.rhs(Phi2, y), (0, Lg * te), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=Lg * tg,
                        max_step=Lg * te / NT_A)
        Wp = np.array([P.w(sol.y[:, k]) for k in range(len(tg))])
        W = Wp[:, :DA] + Wp[:, DA:]
        Uplus = np.array([(P.eta(sol.y[:, k]) - e0)[:DA] / Lg for k in range(len(tg))])   # the atom's common increment / log(1/u)
        Uref2 = np.array([MF.iss_eval(path2, t, Phi, y)[0] for t in tg])
        Wref2 = np.array([MF.iss_eval(path2, t, Phi, y)[1] for t in tg])
        rows2.append(dict(u=u0, du=float(np.abs(Uplus - Uref2).max()),
                          dres=float(np.sqrt(np.mean(np.sum((W @ Phi.T - Wref2 @ Phi.T) ** 2, axis=1)))),
                          dend=float(np.abs(W[-1] - path2['w_end']).max()), n_ent=sum(1 for e in path2['events'] if e[2] == 'enter'),
                          secs=time.time() - t0))
    out['pairs'] = rows2
    return out


def bumps_target(K, cstar, astar):
    return K.phi(np.asarray(cstar, float)) @ np.asarray(astar, float)


def run_moving(L, K, fstar, c0, eps, tau_end, metric, n_t=441, fixed=False, b=None, full=False, taus=None):
    """The free law with moving (or fixed) centres from α = 0, masses εb (b = 1 by default), to τ_end in τ = t/Λ(ε): the state (α, c), each
    kernel kept on its curve β = εb + 2c sinh²(α/a_max) (`critical_rhs` with b0 = εb; round 37). The record on n_t uniform times, or on
    the given grid `taus` (ending at τ_end)."""
    d = len(c0)
    b0 = eps * (np.ones(d) if b is None else np.asarray(b, float))
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_end, n_t) if taus is None else np.asarray(taus, float)
    f = RT.critical_rhs(L, K, fstar, metric, b0=b0)
    rhs = (lambda t, z: np.r_[f(t, z)[:d], np.zeros(d)]) if fixed else f
    sol = solve_ivp(rhs, (0, Lam * taus[-1]), np.r_[np.zeros(d), c0], method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    al, cc = sol.y[:d].T, sol.y[d:].T
    if full:
        return taus, RT.critical_w(L, al, b0), cc, Lam, al
    return taus, RT.critical_w(L, al, b0), cc, Lam


def summarize(K, fstar, c0, taus, W, cc, th):
    act = np.abs(W[-1]) > th
    r = fstar - K.phi(cc[-1]) @ W[-1]
    ent = entries_of(W, taus, th)
    order = [int(i) for i in np.argsort(ent) if np.isfinite(ent[i])]
    cl = np.sort(cc[-1][act])
    clusters = int(np.sum(np.diff(cl) > 1e-3) + 1) if cl.size else 0
    return dict(n_act=int(act.sum()), centres=cc[-1][act].tolist(), w=W[-1][act].tolist(), order=order, ent=[float(ent[i]) for i in order],
                rms=float(np.sqrt(np.mean(r ** 2))), disp_in=float(np.abs(cc[-1] - c0)[~act].max()) if (~act).any() else 0.0,
                w_in=float(np.abs(W[-1][~act]).max()) if (~act).any() else 0.0, clusters=clusters)


def net_entries(R):
    """A sliding path's formed kernels in the order of their (last) entry — exit/re-entry pairs removed (round 37)."""
    last = {}
    for t, j, k, c in R['events']:
        if k == 'enter':
            last[j] = t
    return [int(j) for j, t in sorted(last.items(), key=lambda kv: kv[1]) if R['act'][j]]


def creep_case(L, K, c0, cs, As):
    # (6) a lost near-tie: the target C4's CREEP_TARGET; the kernel left just inside its wall creeps — its instantaneous speed against
    # (πa_max²/(4c))εb^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ (the core's profile; round 37) and against the bound's μ^{−3} form
    f_c = bumps_target(K, cs, As)
    Rc = RT.sliding_iss(L, K, f_c, c0, np.ones(D_C))
    te = 1.3 * Rc['events'][-1][0]
    creep = dict(slide_order=net_entries(Rc), slide_tau=[float(e[0]) for e in Rc['events']], rows=[])
    t2 = Rc['events'][1][0]; t3 = Rc['events'][2][0]
    kc = np.pi * L.A ** 2 / (4 * L.c)
    for eps in EPS_CREEP:
        b0 = np.full(D_C, eps)
        taus, W, cc, Lam, AL = run_moving(L, K, f_c, c0, eps, te, 'paradigm', n_t=881, full=True)
        q = L.Q(np.zeros(D_C), b0)
        lam1 = np.pi * L.A * np.sqrt((L.c - q) / (q + L.c)) / (4 * L.G * L.c)
        ent = entries_of(W, taus, TH_ACT * np.abs(As).max())
        order = [int(i) for i in np.argsort(ent) if np.isfinite(ent[i])]
        k2 = int(np.searchsorted(taus, t2 + 0.2)); kb = int(np.searchsorted(taus, min(t3, ent[order[2]]) - 0.1))
        f = RT.critical_rhs(L, K, f_c, 'paradigm', b0=b0)
        al2 = AL[k2]
        mu = 1 - np.abs(L.eta(al2, q)) / lam1
        near = [int(j) for j in np.flatnonzero((np.abs(W[k2]) < 1e-3) & (mu < 0.01))]   # the inactive kernels lingering at their walls
        rows_near = []
        for j in near:
            pts = []
            for kk in np.linspace(k2, kb, 4).astype(int):
                alk = AL[kk]
                dz = f(0.0, np.r_[alk, cc[kk]])
                rk = f_c - K.phi(cc[kk]) @ W[kk]
                hk = float((K.dphi(cc[kk][j:j + 1]).T @ rk)[0] / K.xs.size)
                muk = float(1 - abs(L.eta(alk[j:j + 1], q[j:j + 1])[0]) / lam1[j])
                sk = float(np.sign(W[kk, j]))
                pts.append(dict(tau=float(taus[kk]), mu=muk, h=hk, ode=float(Lam * dz[D_C + j]),
                                profile=float(kc * eps * np.cos(np.pi * muk / 2) / np.sin(np.pi * muk / 2) ** 3 * sk * hk),
                                old=float(kc * eps * muk ** -3 * abs(hk))))
            rows_near.append(dict(j=j, mu=float(mu[j]), pts=pts, disp_window=float(cc[kb, j] - cc[k2, j]), disp_c0=float(cc[kb, j] - c0[j])))
        creep['rows'].append(dict(eps=eps, order=order, centres=np.sort(cc[-1][np.abs(W[-1]) > TH_ACT * np.abs(As).max()]).tolist(),
                                  near=rows_near))
    return creep


C5_CASES = (('one bump midway between two kernels (a tie)', (0.05,), (0.8,), SIG_C, EPS_C5),
            ('one bump 0.001 off the midpoint (a near tie)', (0.051,), (0.8,), SIG_C, EPS_C5),
            ('three bumps, one at a midpoint', (-0.52, 0.05, 0.61), ASTAR, SIG_C, EPS_C5),
            ('two bumps 0.15 apart, off the grid', (-0.03, 0.12), (0.8, 0.6), SIG_C, EPS_C5),
            ('two bumps 0.25 apart, off the grid', (-0.07, 0.18), (0.8, 0.6), SIG_C, EPS_C5),
            ('three wider bumps (width 0.2)', CSTAR, ASTAR, 0.2, (EPS_RAND,)))


def part_C5(L, K, c0):
    """C (5): ties and near ties, close bumps off the grid, a target of bumps wider than the kernels (not representable by three of them).
    The sliding path's kernels counted, like the loop's, by |w| > TH_ACT·max|A*| (its active set given beside them); a path that does not
    terminate (its rms residual still above 10⁻⁹ when no kernel is left to enter) is followed on a geometric τ grid (round 38), with its
    rms residual at each event and whether each relaxation became stationary."""
    hard = []
    for lab, cs, As, wid, eps_h in C5_CASES:
        Kt = RT.Bumps(K.xs, wid)
        f_h = bumps_target(Kt, cs, As)
        th_h = TH_ACT * np.abs(As).max()
        Rh = RT.sliding_iss(L, K, f_h, c0, np.ones(D_C), max_events=30)
        done = bool(Rh['hist'][-1][3] < 1e-9)
        te = 1.3 * Rh['events'][-1][0]
        formed = np.abs(Rh['w']) > th_h
        rec = dict(lab=lab, slide_n=int(formed.sum()), slide_act=int(Rh['act'].sum()), slide_c=Rh['c'][formed].tolist(),
                   slide_w=Rh['w'][formed].tolist(), slide_rms=Rh['hist'][-1][3], slide_events=len(Rh['events']), slide_done=done,
                   slide_entries=[(float(t), int(jj)) for t, jj, k, c in Rh['events'] if k == 'enter'],
                   slide_hist=[(float(h[0]), float(h[3])) for h in Rh['hist']], slide_relaxed=[bool(v) for v in Rh['relaxed']],
                   te=float(te), loops=[])
        grid = None if done else np.r_[0.0, np.geomspace(1.0, te, 1760)]
        for eps in eps_h:
            taus, W, cc, Lam = run_moving(L, K, f_h, c0, eps, te, 'paradigm', n_t=1761, taus=grid)
            sm = summarize(K, f_h, c0, taus, W, cc, th_h)
            rec['loops'].append(dict(eps=eps, grid='uniform' if grid is None else 'geometric',
                                     **{k: sm[k] for k in ('n_act', 'centres', 'w', 'rms', 'clusters', 'order', 'ent')}))
        hard.append(rec)
        print(f"  C5 {lab}: slide {rec['slide_n']} (active {rec['slide_act']}), loop {[r['n_act'] for r in rec['loops']]}", flush=True)
    return hard


def part_C():
    warnings.filterwarnings('ignore')
    L = law()
    K = RT.Bumps(np.linspace(-1, 1, NS_C), SIG_C)
    c0 = np.linspace(-1, 1, D_C)
    fstar = bumps_target(K, CSTAR, ASTAR)
    th = TH_ACT * np.abs(ASTAR).max()
    out = {}
    # (1) the paradigm's metric against the sliding inverse scale space
    R = RT.sliding_iss(L, K, fstar, c0, np.ones(D_C))
    out['slide'] = dict(events=[(float(t), int(j), k, float(c)) for t, j, k, c in R['events']], centres=R['c'][R['act']].tolist(),
                        w=R['w'][R['act']].tolist(), rms=R['hist'][-1][3])
    tau_end = 1.3 * R['events'][-1][0]
    rows = []
    for eps in EPS_C1:
        t0 = time.time()
        taus, W, cc, Lam = run_moving(L, K, fstar, c0, eps, tau_end, 'paradigm')
        rows.append(dict(eps=eps, Lam=float(Lam), secs=time.time() - t0, **summarize(K, fstar, c0, taus, W, cc, th)))
    out['c1'] = rows
    out['dtau1'] = float(tau_end / 440)
    # (2) the Wasserstein–Fisher–Rao metric: the inactive kernels against their limit along the trajectory before the first entry (equal
    # and unequal masses), and the formed kernels as ε → 0 (round 37: the count falls to one per feature)
    def wfr_case(b, eps_list, tau_end_w):
        rho = b ** -0.5
        lim = solve_ivp(RT.wfr_inactive_rhs(L, K, fstar, rho), (0, tau_end_w), np.r_[np.zeros(D_C), c0], method='LSODA', rtol=1e-10,
                        atol=1e-12, dense_output=True)
        lim1 = solve_ivp(RT.wfr_inactive_rhs(L, K, fstar, np.ones(D_C)), (0, tau_end_w), np.r_[np.zeros(D_C), c0], method='LSODA',
                         rtol=1e-10, atol=1e-12, dense_output=True)          # the masses ignored (ρ = 1): a wrong limit, for contrast
        rows = []
        for eps in eps_list:
            t0 = time.time()
            taus, W, cc, Lam = run_moving(L, K, fstar, c0, eps, tau_end_w, 'wfr', n_t=881, b=b)
            first = np.flatnonzero((np.abs(W) > th).any(axis=1))
            t_first = taus[first[0]] if first.size else taus[-1]
            pre = taus < 0.9 * t_first
            cl = lim.sol(taus[pre])[D_C:].T
            rows.append(dict(eps=eps, t_first=float(t_first), d_pre=float(np.abs(cc[pre] - cl).max()),
                             d_pre_rho1=float(np.abs(cc[pre] - lim1.sol(taus[pre])[D_C:].T).max()), w_pre=float(np.abs(W[pre]).max()),
                             secs=time.time() - t0, **summarize(K, fstar, c0, taus, W, cc, th)))
        return rows, lim
    rows, lim = wfr_case(np.ones(D_C), EPS_C2, tau_end)
    c_lim = lim.sol(TAU_PRE)[D_C:]
    cl = np.sort(c_lim)
    b_un = np.exp(np.random.default_rng(31).uniform(np.log(0.25), np.log(4.0), D_C))
    rows_b, lim_b = wfr_case(b_un, (EPS_C2B,), tau_end)
    out['c2'] = dict(rows=rows, c_lim=c_lim.tolist(), lim_clusters=int(np.sum(np.diff(cl) > 1e-3) + 1), unequal=rows_b,
                     b_range=[float(b_un.min()), float(b_un.max())])
    # (3) fixed centres (the grid dictionary)
    taus, W, cc, Lam = run_moving(L, K, fstar, c0, EPS_RAND2, tau_end * 3, 'paradigm', fixed=True)
    Phi = K.phi(c0)
    wls = np.linalg.lstsq(Phi, fstar, rcond=None)[0]
    out['c3'] = dict(eps=EPS_RAND2, tau=float(taus[-1]), **summarize(K, fstar, c0, taus, W, cc, th),
                     rms_ls=float(np.sqrt(np.mean((fstar - Phi @ wls) ** 2))))
    # (4) random three-bump targets
    rng = np.random.default_rng(21)
    rnd = []
    for n in range(N_RAND):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
        f_n = bumps_target(K, cs, As)
        th_n = TH_ACT * np.abs(As).max()
        Rn = RT.sliding_iss(L, K, f_n, c0, np.ones(D_C))
        te = 1.3 * Rn['events'][-1][0]
        rec = dict(cs=cs.tolist(), As=As.tolist(), slide_n=int(Rn['act'].sum()), slide_rms=Rn['hist'][-1][3],
                   slide_err=float(np.abs(np.sort(Rn['c'][Rn['act']]) - cs).max()) if Rn['act'].sum() == 3 else np.inf)
        rec['slide_order'] = net_entries(Rn)
        rec['slide_exits'] = sum(1 for e in Rn['events'] if e[2] == 'exit')
        runs = [('paradigm', 'paradigm', EPS_RAND), ('paradigm2', 'paradigm', EPS_RAND2), ('paradigm3', 'paradigm', EPS_RAND3)]
        runs += [(f'wfr{k}', 'wfr', e) for k, e in enumerate(EPS_WFR_RAND)]
        for key, metric, eps_n in runs:
            taus, W, cc, Lam = run_moving(L, K, f_n, c0, eps_n, te, metric)
            sm = summarize(K, f_n, c0, taus, W, cc, th_n)
            rec[key] = dict(n_act=sm['n_act'], clusters=sm['clusters'], rms=sm['rms'], disp_in=sm['disp_in'], order=sm['order'],
                            err=float(np.abs(np.sort(sm['centres']) - cs).max()) if sm['n_act'] == 3 else np.inf)
        taus, W, cc, Lam = run_moving(L, K, f_n, c0, EPS_RAND2, 3 * te, 'paradigm', fixed=True)
        sm = summarize(K, f_n, c0, taus, W, cc, th_n)
        rec['fixed'] = dict(n_act=sm['n_act'], rms=sm['rms'])
        rnd.append(rec)
        print(f"  C4 target {n}: slide {rec['slide_n']} paradigm {rec['paradigm']['n_act']}, {rec['paradigm2']['n_act']}, "
              f"{rec['paradigm3']['n_act']} wfr {[rec[f'wfr{k}']['n_act'] for k in range(len(EPS_WFR_RAND))]} fixed {rec['fixed']['n_act']}",
              flush=True)
    out['c4'] = rnd
    # (5) ties and near ties (round 37), close bumps off the grid, a target wider than the kernels
    hard = part_C5(L, K, c0)
    out['c5'] = hard
    rc = rnd[CREEP_TARGET]
    creep = creep_case(L, K, c0, np.array(rc['cs']), np.array(rc['As']))
    out['c6'] = creep
    return out


def part_D():
    """Moving widths: the free law with moving centres and widths against fixed widths, on bumps wider and narrower than the kernels."""
    warnings.filterwarnings('ignore')
    L = law()
    K = RT.Bumps(np.linspace(-1, 1, NS_C), SIG_C)
    c0 = np.linspace(-1, 1, D_C)
    th = TH_ACT * np.abs(ASTAR).max()
    rows = []
    for wid in WIDTHS_D:
        fs = RT.Bumps(K.xs, wid).phi(np.asarray(CSTAR)) @ np.asarray(ASTAR)
        for eps in EPS_D:
            for moving in (True, False):
                Lam = RT.transit_scale(L, eps)
                taus = np.linspace(0, TAU_D, 601)
                z0 = np.r_[np.zeros(D_C), np.log(np.full(D_C, eps)), c0, np.full(D_C, SIG_C)]
                f = RT.moving_rhs_widths(L, K, fs)
                rhs = f if moving else (lambda t, z, f=f: np.r_[f(t, z)[:3 * D_C], np.zeros(D_C)])
                sol = solve_ivp(rhs, (0, Lam * TAU_D), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
                W = np.exp(sol.y[D_C:2 * D_C].T) * L.amp(sol.y[:D_C].T)[0] / L.m
                cc, sg = sol.y[2 * D_C:3 * D_C].T, sol.y[3 * D_C:].T
                act = np.abs(W[-1]) > th
                Ph = np.exp(-(K.xs[:, None] - cc[-1][None, :]) ** 2 / (2 * sg[-1] ** 2))
                o = np.argsort(cc[-1][act])
                rows.append(dict(wid=wid, eps=eps, moving=moving, n_act=int(act.sum()), centres=cc[-1][act][o].tolist(),
                                 widths=sg[-1][act][o].tolist(), w=W[-1][act][o].tolist(),
                                 rms=float(np.sqrt(np.mean((fs - Ph @ W[-1]) ** 2))),
                                 order=[i for t, i in sorted((taus[np.flatnonzero(np.abs(W[:, i]) > th)[0]], i) for i in range(D_C)
                                                             if (np.abs(W[:, i]) > th).any())]))
    return rows


def report(A=None, B=None, C=None, D=None):
    fmt = E.fmt
    L = law()
    L_ = ["# Proposition D‴.48: the rich transit (`exp_mirror_transit.py`)", "",
          f"The free mass–amplitude law (Proposition D‴.45; m = {LAW['m']:g}, G = {LAW['G']:g}, γ = {LAW['gamma']:g}, a_max = {LAW['a_max']:g}, "
          f"c = {L.c:g}) from small masses, in the time rescaled by the transit's scale (`mirror_free.iss_path`, `rich_transit`).", ""]
    if A is not None:
        L_ += ["## A. From α = 0: the inverse scale space flow on the scale Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}", "",
               f"Three {NA} × {DA} Gaussian designs, 3-sparse truths, masses εb_i (b_i log-uniform in [1/4, 4]); the path of Σ|w_i|/√b_i against "
               "the loop in τ = t/Λ, u = η/Λ, on 2001 times in [0, 1.5τ*] (τ* the last event): sup_τ|u − u_ISS|, the residual's difference "
               "|Φ(w − w_ISS)| over the samples — its rms, its median and the number of samples where it exceeds 10⁻³ (at an event the "
               "difference is O(1) over a τ-width ∝ ε^{1/2}, so its L² norm over τ is ∝ ε^{1/4}, while away from the events it is O(ε^{3/2}), the "
               "inactive atoms' coefficients; once the width falls below the grid's step the rms counts only the samples that fall within an "
               "entry), the largest difference of the entry times (the first sample where |w_i| > 10⁻³ max|w_bp|, on the loop and on the path; "
               "in brackets the atoms crossing on one of them only; the grid's step given), the end's distance to weighted basis pursuit.", ""]
        for a in A:
            ev = '; '.join(f"{k} {i} at {t:.4g}" for t, i, k in a['events'])
            L_ += [f"Design {a['seed']}: the path's events — {ev}; its end against weighted basis pursuit {fmt(a['end_is_bp'])}; the grid's step "
                   f"{fmt(a['dtau'])}.", "",
                   "| ε | Λ | sup|u − u_ISS| | residual's difference: rms, median, samples > 10⁻³ | entry times: largest difference | "
                   "|w(1.5τ*) − w_bp|∞ |", "|---|---|---|---|---|---|"]
            L_ += [f"| {r['eps']:g} | {r['Lam']:.4g} | {fmt(r['du'])} | {fmt(r['dres'])}, {fmt(r['dres_med'])}, {r['dres_n']} | "
                   f"{fmt(r['dent'])} ({r['n_only']}) | {fmt(r['dend'])} |" for r in a['rows']]
            du = {r['eps']: r['du'] for r in a['rows']}
            sl = f"the slope {np.log(du[1e-6] / du[1e-8]) / np.log(100):.3f} between ε = 10⁻⁶ and 10⁻⁸"
            if 1e-12 in du:
                sl += f", {np.log(du[1e-8] / du[1e-12]) / np.log(1e4):.3f} between 10⁻⁸ and 10⁻¹²"
            dm = {r['eps']: r['dres_med'] for r in a['rows']}
            dr = {r['eps']: r['dres'] for r in a['rows']}
            if 1e-12 in dm:
                sl += (f"; the residual's difference — rms: the slope {np.log(dr[1e-6] / dr[1e-8]) / np.log(100):.3f} between 10⁻⁶ and 10⁻⁸ "
                       f"(the L² norm's ε^{{1/4}}: 0.25); median: {np.log(dm[1e-6] / dm[1e-10]) / np.log(1e4):.3f} between 10⁻⁶ and 10⁻¹⁰ "
                       "(ε^{3/2}: 1.5; at 10⁻¹² it nears the rounding floor)")
            L_ += ["", f"sup|u − u_ISS| in log ε: {sl} (1/Λ ∝ ε^{{1/2}}: 0.5). Each atom is integrated on its own curve (α alone, "
                   "β = εb + 2c sinh²(α/a_max)), so Q is kept exactly.", ""]
    if B is not None:
        L_ += ["## B. On branches: the logarithmic scale", "",
               f"Design {SEEDS_A[0]} of A; amplitudes s_i·U(0.2, 0.6) (the truth's signs on its atoms, random elsewhere), masses ε; the one-sided "
               "path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in τ = t/log(1/ε), u = (η − η(0))/log(1/ε). The path's events: "
               + '; '.join(f"{k} {i} at {t:.4g}" for t, i, k in B['events']) + f"; its end against the sign-constrained weighted basis pursuit "
               f"{fmt(B['end_is_lim'])}.", "",
               "| ε | all on branches | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_lim|∞ |",
               "|---|---|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {'yes' if r['all_branch'] else 'no'} | {fmt(r['du'])} | {fmt(r['dres'])} | {fmt(r['dent'])} ({r['n_only']}) | "
               f"{fmt(r['dend'])} |" for r in B['rows']]
        L_ += ["", "A ± pair per atom (amplitudes ±0.4, equal masses, w(0) = u): the atom's common increment over log(1/u) against the two-sided "
               "path with the weights 1/(γa₀):", "", "| u | sup|u − u_ISS| | residual: rms difference | |W(1.5τ*) − w_end|∞ | entries on the path |",
               "|---|---|---|---|---|"]
        L_ += [f"| {r['u']:g} | {fmt(r['du'])} | {fmt(r['dres'])} | {fmt(r['dend'])} | {r['n_ent']} |" for r in B['pairs']]
        L_.append("")
    if C is not None:
        sl = C['slide']
        L_ += ["## C. Moving centres", "",
               f"Gaussian kernels of width {SIG_C} on {NS_C} samples of [−1, 1]; {D_C} kernels started on a grid of [−1, 1] (spacing 0.1), α = 0, "
               f"masses ε; the target {len(CSTAR)} bumps of the same width at {list(CSTAR)} with amplitudes {list(ASTAR)}. A kernel has formed when "
               f"|w| > {TH_ACT:g}·max|A*|.", "",
               "(1) The paradigm's metric (ċ = Gβa⟨r∂_cφ⟩). The sliding inverse scale space: " +
               '; '.join(f"{k} kernel {j} (start {c:+.1f}) at τ = {t:.4g}" for t, j, k, c in sl['events']) +
               f"; it ends with {len(sl['centres'])} kernels at {[round(v, 5) for v in sl['centres']]} with coefficients {[round(v, 5) for v in sl['w']]} "
               f"(rms residual {fmt(sl['rms'])}). The loop (the entries at the first sample — step {C['dtau1']:.3g} in τ — where |w| > {TH_ACT:g}·max|A*|):", "",
               "| ε | Λ | entries (kernel: τ) | formed | their centres | their coefficients | rms residual | inactive: max |Δc|, max |w| |",
               "|---|---|---|---|---|---|---|---|"]
        for r in C['c1']:
            ent = ', '.join(f"{i}: {t:.4g}" for i, t in zip(r['order'], r['ent']))
            L_.append(f"| {r['eps']:g} | {r['Lam']:.4g} | {ent} | {r['n_act']} | {[round(v, 4) for v in r['centres']]} | {[round(v, 4) for v in r['w']]} | "
                      f"{fmt(r['rms'])} | {fmt(r['disp_in'])}, {fmt(r['w_in'])} |")
        c2 = C['c2']
        L_ += ["", f"(2) The Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩). The inactive kernels' limit (du/dτ = ⟨f*φ⟩, dc/dτ = "
               f"(πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩) gathers the {D_C} kernels at {c2['lim_clusters']} positions by τ = {TAU_PRE:g} "
               f"{sorted(set(round(v, 3) for v in c2['c_lim']))}. The loop against that limit along the trajectory (the largest |c − c_lim| over "
               "τ < 0.9τ_first, τ_first the first time a kernel's |w| exceeds the threshold), and its end:", "",
               "| ε | τ_first | max |c − c_lim| before | max |w| before | formed | their positions (clusters) | rms residual |", "|---|---|---|---|---|---|---|"]
        for r in c2['rows']:
            L_.append(f"| {r['eps']:g} | {r['t_first']:.4g} | {fmt(r['d_pre'])} | {fmt(r['w_pre'])} | {r['n_act']} | "
                      f"{sorted(set(round(v, 3) for v in r['centres']))} ({r['clusters']}) | {fmt(r['rms'])} |")
        for r in c2['unequal']:
            L_ += ["", f"Unequal masses (εb_i, b_i log-uniform in [{c2['b_range'][0]:.3g}, {c2['b_range'][1]:.3g}], ε = {r['eps']:g}): the "
                   f"limit with ρ_i = b_i^{{−1/2}} followed to {fmt(r['d_pre'])} before τ_first = {r['t_first']:.4g} (the limit with the masses "
                   f"ignored, ρ = 1: {fmt(r['d_pre_rho1'])}); formed {r['n_act']} kernels in {r['clusters']} clusters."]
        c3 = C['c3']
        L_ += ["", f"(3) The centres held on the grid (ε = {c3['eps']:g}, to τ = {c3['tau']:.3g}): {c3['n_act']} kernels formed, rms residual "
               f"{fmt(c3['rms'])} (the grid's least-squares fit: {fmt(c3['rms_ls'])}).", ""]
        c4 = C['c4']
        n = len(c4)
        ok_s = sum(1 for r in c4 if r['slide_n'] == 3 and r['slide_err'] < 1e-3)
        ok = {k: sum(1 for r in c4 if r[k]['n_act'] == 3 and r[k]['err'] < 1e-3) for k in ('paradigm', 'paradigm2', 'paradigm3')}
        wfr = '; '.join(f"ε = {e:g}: {sorted(r[f'wfr{k}']['n_act'] for r in c4)} kernels, exactly three at the bumps on "
                        f"{sum(1 for r in c4 if r[f'wfr{k}']['n_act'] == 3 and r[f'wfr{k}']['err'] < 1e-3)} of {n}"
                        for k, e in enumerate(EPS_WFR_RAND))
        L_ += [f"(4) {n} random targets (three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)): "
               f"the sliding inverse scale space ends with exactly three kernels within 10⁻³ of the bumps' centres on {ok_s} of {n} "
               f"({sum(r['slide_exits'] for r in c4)} exits in all); the loop with the paradigm's metric on {ok['paradigm']}, {ok['paradigm2']} and "
               f"{ok['paradigm3']} of {n} at ε = {EPS_RAND:g}, {EPS_RAND2:g}, {EPS_RAND3:g} (formed kernels {sorted(r['paradigm']['n_act'] for r in c4)}, "
               f"{sorted(r['paradigm2']['n_act'] for r in c4)}, {sorted(r['paradigm3']['n_act'] for r in c4)}; the inactive kernels' largest "
               f"displacement {fmt(max(r['paradigm']['disp_in'] for r in c4))}, {fmt(max(r['paradigm2']['disp_in'] for r in c4))}, "
               f"{fmt(max(r['paradigm3']['disp_in'] for r in c4))}), the formed kernels those of the sliding path, in its order, on "
               f"{sum(1 for r in c4 if r['paradigm3']['order'] == r['slide_order'])} of {n} at ε = {EPS_RAND3:g}; with the Wasserstein–Fisher–Rao "
               f"metric {wfr}; with the centres held on the grid {sorted(r['fixed']['n_act'] for r in c4)} kernels (rms residual "
               f"{fmt(min(r['fixed']['rms'] for r in c4))}–{fmt(max(r['fixed']['rms'] for r in c4))}).", ""]
        L_ += ["(5) Ties, near ties, close bumps off the grid, wider bumps (the paradigm's metric): the sliding inverse scale space (its "
               "entries; a tie enters together), and the loop at each ε:", "",
               "| target | sliding ISS: entries (τ, kernel), kernels, centres, coefficients, rms residual | ε | loop: entries (kernel: τ), formed, "
               "centres, coefficients (clusters), rms residual |", "|---|---|---|---|"]
        for r in C['c5']:
            if r['slide_done']:
                sl = (f"{[(round(t, 3), jj) for t, jj in r['slide_entries']]}, {r['slide_n']}, {[round(v, 3) for v in r['slide_c']]}, "
                      f"{[round(v, 3) for v in r['slide_w']]}, {fmt(r['slide_rms'])}")
            else:
                sl = (f"{len(r['slide_entries'])} entries, at τ = {', '.join(f'{t:.3g}' for t, jj in r['slide_entries'])} (not terminated: "
                      f"{r['slide_act']} kernels active, {r['slide_n']} with |w| above the threshold, rms residual {fmt(r['slide_rms'])})")
            for k, lp in enumerate(r['loops']):
                ent = ', '.join(f"{i}: {t:.4g}" for i, t in zip(lp['order'], lp['ent']))
                cw = (f"{[round(v, 3) for v in lp['centres']]}, {[round(v, 3) for v in lp['w']]} " if r['slide_done'] else "")
                L_.append(f"| {r['lab'] if k == 0 else ''} | {sl if k == 0 else ''} | {lp['eps']:g} | {ent}; {lp['n_act']}, "
                          f"{cw}({lp['clusters']}), {fmt(lp['rms'])} |")
        L_.append("")
        for r in C['c5']:
            if not r['slide_done'] and 'slide_hist' in r:
                hh = r['slide_hist']
                L_ += [f"The target '{r['lab']}' is not three kernels of width {SIG_C}: the sliding path enters a kernel at a time, at times "
                       f"growing geometrically, its rms residual at the events {', '.join(fmt(v) for t, v in hh)} (at τ = "
                       f"{', '.join(f'{t:.3g}' for t, v in hh)}); {sum(r['slide_relaxed'])} of its {len(r['slide_relaxed'])} relaxations "
                       f"became stationary within the budget."]
                path_order = [jj for t, jj in r['slide_entries']]
                for lp in r['loops']:
                    n_ = min(len(lp['order']), len(path_order))
                    k = 0
                    while k < n_ and lp['order'][k] == path_order[k]:
                        k += 1
                    ks = max([j for j in range(1, n_ + 1) if set(lp['order'][:j]) == set(path_order[:j])], default=0)
                    L_[-1] += (f" The loop at ε = {lp['eps']:g} (on a geometric τ grid to τ = {r['te']:.3g}) enters {len(lp['order'])} kernels, "
                               f"at τ = {lp['ent'][0]:.3g} to {lp['ent'][-1]:.3g}: the first {k} of its entrants the path's in its order, the first "
                               f"{ks} the path's first {ks} (its entries at τ = {', '.join(f'{t:.3g}' for t in lp['ent'][:ks])}; the path's "
                               f"{', '.join(f'{t:.3g}' for t, jj in r['slide_entries'][:ks])}); it ends with {lp['n_act']} kernels above the "
                               f"threshold at the rms residual {fmt(lp['rms'])} (the path {fmt(r['slide_rms'])}).")
                L_.append("")
    if C is not None and 'c6' in C:
        c6 = C['c6']
        L_ += ["", f"(6) A lost near-tie (target {CREEP_TARGET} of (4)): the sliding path enters kernels {c6['slide_order']} at τ = "
               f"{[round(t, 3) for t in c6['slide_tau']]}. After the second entry an inactive kernel left near its wall lingers at a small margin μ "
               "and creeps; its instantaneous speed dc/dτ (the loop's own right-hand side) against the core's profile "
               "(πa_max²/(4c))εb^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ and against the bound's μ^{−3} form (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩|:", "",
               "| ε | kernel | τ | μ | ⟨r∂_cφ⟩ | dc/dτ | the profile's speed (ratio) | the μ^{−3} form (ratio) |", "|---|---|---|---|---|---|---|---|"]
        for r in c6['rows']:
            for nn in r['near']:
                for p in nn['pts']:
                    L_.append(f"| {r['eps']:g} | {nn['j']} | {p['tau']:.4g} | {fmt(p['mu'])} | {fmt(p['h'])} | {fmt(p['ode'])} | {fmt(p['profile'])} "
                              f"({p['ode'] / p['profile']:.4f}) | {fmt(p['old'])} ({abs(p['ode']) / p['old']:.4f}) |")
        L_ += ["", "| ε | kernel: displacement over the window (since τ_2 + 0.2), since its start | the loop's entrants | the formed kernels' centres |",
               "|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | " + '; '.join(f"{nn['j']}: {fmt(nn['disp_window'])}, {fmt(nn['disp_c0'])}" for nn in r['near'])
               + f" | {r['order']} | {[round(v, 4) for v in r['centres']]} |" for r in c6['rows']]
        L_.append("")
    if D is not None:
        L_ += ["## D. Moving widths (σ̇ = Gβa⟨r∂_σφ⟩ beside ċ = Gβa⟨r∂_cφ⟩; kernels started at width 0.12)", "",
               f"The target: the three bumps of C at widths {list(WIDTHS_D)}; τ ≤ {TAU_D:g}.", "",
               "| the bumps' width | ε | widths | entrants | formed | centres | widths | coefficients | rms residual |", "|---|---|---|---|---|---|---|---|---|"]
        L_ += [f"| {r['wid']:g} | {r['eps']:g} | {'moving' if r['moving'] else 'fixed'} | {r['order']} | {r['n_act']} | "
               f"{[round(v, 4) for v in r['centres']]} | {[round(v, 4) for v in r['widths']]} | {[round(v, 4) for v in r['w']]} | {fmt(r['rms'])} |"
               for r in D]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_transit.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_transit_raw.npz')
    md = os.path.join(OUT, 'mirror_transit.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    if 'C5' in parts and 'C' not in parts:                            # C (5) alone, into the stored C (round 38)
        warnings.filterwarnings('ignore')
        res['C'] = dict(res['C'], c5=part_C5(law(), RT.Bumps(np.linspace(-1, 1, NS_C), SIG_C), np.linspace(-1, 1, D_C)))
        E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] C5", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
