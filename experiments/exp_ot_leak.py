"""
Proposition D‴.12: the barycentric OT leak toward a prior measure — the mass constraint against the weight fluctuation
=====================================================================================================================

The averaged slow loop of Theorem D‴ on the three-bump instance of `exp_chaos_learning.py` at σ₀ = 0.15σ_crit, with the matching prior
of `exp_prior_design.py` (three atoms at the exact representation, masses q = (1/3, 1/3, 1/3) = the plant's representation weights), under
two leaks on the same clouds: the nearest-atom leak (𝔄_m^ν) of Proposition D‴.9, which conserves the cluster masses, and the barycentric
OT leak (𝔄_m^{OT}) of Theorem C′, which constrains them to q.

1. Commensurate sizes (mq ∈ ℕ³): m = 15, 60, 240 against M = 510, warm start = the jittered mixture around the exact representation,
   600 slow units. Records: the cluster weights (conserved against driven to q; the equalisation time), the mass in transit, the fit
   errors, and — against the reference — the synchronous coupled distance D, the function-level distance ‖f̂_m − f̂_M‖²_{L²(ρ)} and the exact
   W²_{G⁻¹}(μ_m, μ_M); the last-window values of the three distances under the two leaks, and the nearest-atom leak's predicted function-level
   floor (1/m − 1/M)tr(FC̃Fᵀ)/N.
2. Incommensurate sizes: m = 16, 64, 256 against M = 512 (one or two split particles): the same records, the split particles at the end
   (their plan rows and positions), and the m-laws of the stationary floors under the two leaks.
3. The masses as a design input (no dynamics): for q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) the weighted stationary point Θ*(q_δ), its fit
   error (the bias of the mass constraint), κ_R at the atoms and over N_{ϱ′}, the margin, the splitting margin (admissibility of the
   designed masses) and the half-separation check ε < Δ_ν/2 of the diagonal plan.

Output: results/chaos/ot_leak.md and ot_leak_raw.npz (CLW_WORKERS, OTL_SEEDS, OTL_S_END to change the pool).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import prior_design as PD
from src.core import ot_leak as OT

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846                                       # Proposition K on this instance (asserted in the tests)
SIGMA0 = 0.15 * SIGMA_CRIT
SETS = {'commensurate': (510, (15, 60, 240)), 'incommensurate': (512, (16, 64, 256))}
S_END, DT_S, REC, PLAN_EVERY = float(os.environ.get('OTL_S_END', 600.0)), 0.02, 50, 10
WINDOWS = tuple((a * S_END / 600.0, b * S_END / 600.0) for a, b in ((0.0, 75.0), (75.0, 150.0), (150.0, 300.0), (300.0, 450.0), (450.0, 600.0)))
RHO_PRIME = 0.1
JITTER = 0.02
N_SEEDS = int(os.environ.get('OTL_SEEDS', 4))
N_WORKERS = int(os.environ.get('CLW_WORKERS', 3))
Q = np.ones(3) / 3
DELTAS = (0.0, 0.01, 0.02, 0.03, 0.05, 0.1)


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def make_loop(net, q=None):
    return OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact_representation(net), q=q, n_quad=N_Q)


# ---------------------------------------------------------------------------- coupled runs
def coupled_job(args):
    seed, set_name, leak = args
    M, ms = SETS[set_name]
    net = make_net()
    exact = exact_representation(net)
    loop = make_loop(net)
    rng = np.random.default_rng(seed)
    Th_all = CL.mixture_start(exact, M, rng, JITTER, G)
    rec = loop.coupled_run_ot(Th_all, ms, S_END, DT_S, rec_every=REC, plan_every=PLAN_EVERY, leak=leak, label_atoms=exact, w2_every=1)
    out = dict(seed=seed, set=set_name, leak=leak, M=M, ms=ms, s=rec['s'])
    for m in list(ms) + [M]:
        out[f'w{m}'] = rec['w'][m]; out[f'transit{m}'] = rec['transit'][m]; out[f'nsplit{m}'] = rec['n_split'][m]; out[f'fit{m}'] = rec['fit'][m]; out[f'W2nu{m}'] = rec['W2nu'][m]
        out[f'final{m}'] = rec['final'][m]
        if leak == 'ot':
            out[f'bridge{m}'] = OT.bridge_report(loop, rec['final'][m], rec['final_plan'][m])
    for m in ms:
        out[f'D{m}'] = rec['D'][m]; out[f'h2{m}'] = rec['hnorm2'][m]; out[f'W2{m}'] = rec['W2'][m]
    return out


# ---------------------------------------------------------------------------- the masses as a design input
def mass_design_row(delta):
    net = make_net()
    exact = exact_representation(net)
    q = Q + delta * np.array([1.0, -0.5, -0.5])
    loop = make_loop(net, q=q)
    atoms, res, (eps, half) = loop.stationary_atoms(exact, w0=Q, n_cont=4)
    a, c, s, diff, phi = loop.features(atoms)
    r = loop.F - phi @ (q * a)                                                                     # residual of the weighted measure Σ_k q_kδ_{ϑ*_k}
    k_at = loop.kappa_R(r, atoms)
    k_nb = loop.kappa_R_neighbourhood(r, atoms, RHO_PRIME, n_dirs=200, radii=(1.0, 0.5))
    fit = float(np.sqrt(np.mean(r ** 2) / np.mean(loop.F ** 2)))
    a, c, s = net.physical(atoms)
    return dict(delta=delta, q=q, atoms=atoms, res=res, eps=eps, half=half, fit=fit, k_at=k_at, k_nb=k_nb, margin=SIGMA0 - k_nb,
                split_margin=SIGMA0 - k_at.max(), a=a[:, 0], c=c[:, 0], s=s)


# ---------------------------------------------------------------------------- statistics
def window_mean(s, v, a, b):
    msk = (s > a) & (s <= b)
    return float(np.nanmean(v[msk]))


def transit_end_time(s, transit, nsplit, m):
    """The first record time from which the mass in transit stays at or below (split particles)/m — no unsplit particle is assigned away
    from its nearest atom any more; nan if never.  For commensurate m this is the equalisation time."""
    ok = np.asarray(transit) <= np.asarray(nsplit) / m + 1e-12
    if not ok[-1]:
        return float('nan')
    idx = len(ok) - 1
    while idx > 0 and ok[idx - 1]:
        idx -= 1
    return float(s[idx])


def mean_se(v):
    v = np.asarray(v, float)
    return float(np.mean(v)), (float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float('nan'))


def fmt(x, se=None, p=3):
    if not np.isfinite(x):
        return "—"
    return f"{x:.{p}e}" if se is None or not np.isfinite(se) else f"{x:.{p}e} ± {se:.1e}"


def slope_se(vals_by_m):
    ms = sorted(vals_by_m); n = len(vals_by_m[ms[0]])
    sl = []
    for i in range(n):
        y = np.array([vals_by_m[m][i] for m in ms])
        if np.all(y > 0):
            sl.append(np.polyfit(np.log(ms), np.log(y), 1)[0])
    if not sl:
        return float('nan'), float('nan')
    return float(np.mean(sl)), (float(np.std(sl, ddof=1) / np.sqrt(len(sl))) if len(sl) > 1 else float('nan'))


def start_clouds(seed, M):
    """The warm clouds of `coupled_job` regenerated from the seed (the systems are autonomous: the first m rows are the m-system's start)."""
    net = make_net()
    exact = exact_representation(net)
    rng = np.random.default_rng(seed)
    return CL.mixture_start(exact, M, rng, JITTER, G)


def w2_weighted_atoms(w_a, atoms_a, w_b, atoms_b):
    """W²_{G⁻¹} between the two weighted K-atom measures (the K × K transportation LP)."""
    C = OT._cost_matrix(atoms_a, atoms_b, G)
    gamma = OT.transport_plan(C, np.asarray(w_a, float), np.asarray(w_b, float))
    return float(np.sum(gamma * C))


def report(runs, design, t0=None):
    net = make_net()
    exact = exact_representation(net)
    loop = make_loop(net)
    sqG = np.sqrt(G)
    inter = OT._cost_matrix(exact, exact, G)                                                   # Δ²_{kk'} between the prior atoms
    # the nearest-atom leak's predicted floors at the uniform weights (Proposition D‴.9(iii))
    atoms, res = PD.weighted_stationary_atoms_prior(loop.nearest, Q, exact)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, atoms)
    F = PD.function_sensitivity(loop.nearest, Q, atoms, S)
    Delta = loop.Delta_nu
    lines = ["# Proposition D‴.12: the barycentric OT leak toward a prior measure — the mass constraint against the weight fluctuation (`exp_ot_leak.py`)", "",
             f"The three-bump instance of `exp_chaos_learning.py` (n = 1, a_max = {NET_KW['a_max']}, c_max = {NET_KW['c_max']}, s ∈ [{NET_KW['s_min']}, {NET_KW['s_max']}], G = diag{tuple(float(g) for g in G)}, "
             f"(A, C, S) = {BUMPS}, x_m = 2 sin t, ρ on {N_Q} times) at σ₀ = 0.15σ_crit = {SIGMA0:.5f}, with the matching prior ν = Σ_j q_jδ_{{ϑ⁰_j}} of `exp_prior_design.py` (three atoms at the exact "
             f"representation (3A_k, C_k, S_k), masses q = (1/3, 1/3, 1/3), minimal inter-atom distance Δ_ν = {Delta:.3f} in the G⁻¹ metric). Two leaks on the same clouds: the nearest-atom leak (𝔄_m^ν) and the "
             f"barycentric OT leak (𝔄_m^{{OT}}) with the plan recomputed every {PLAN_EVERY} steps of {DT_S} slow units; warm start = the mixture of Gaussians of width {JITTER:g} around the exact representation; "
             f"{S_END:g} slow units ({S_END * SIGMA0:.1f}/σ₀); {N_SEEDS} seeds; records every {REC * DT_S:g} slow units. The reference M stands in for the mean-field flow; the m-systems share their particles with it.", ""]
    for set_name, (M, ms) in SETS.items():
        R = {leak: [r for r in runs if r['set'] == set_name and r['leak'] == leak] for leak in ('ot', 'nearest')}
        lines += [f"## {'1' if set_name == 'commensurate' else '2'}. {set_name.capitalize()} sizes: m = {', '.join(str(m) for m in ms)} against M = {M}" + (" (mq ∈ ℕ³: every plan integral)" if set_name == 'commensurate' else " (mq ∉ ℕ³: one or two split particles)"), ""]
        # weights, transit, equalisation
        lines += ["### Cluster masses: conserved by the nearest-atom leak, driven to q by the OT leak", "",
                  "| m | mean |w(0) − q|₁ over seeds (both leaks start from the same clouds) | ½|w^m(0) − w^M(0)|₁ (mean; the mass the two systems' clusters disagree on) | nearest-atom leak: max_s |w(s) − w(0)|_∞ | OT leak: |w(end) − q|_∞ (by nearest atom) | OT leak: equalisation time s_eq (first record from which w = q to within half a particle; per seed) | OT leak: end of the transit s_tr (first record from which no unsplit particle is assigned away from its nearest atom; per seed; recorded series — on the incommensurate sizes the recorded mask counted split rows, so only the commensurate column is exact) | OT leak: unsplit particles in transit at s = 0 → at the end (per seed; recomputed from the clouds) | split particles (OT, end) |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for m in list(ms) + [M]:
            w0 = [float(np.abs(r[f'w{m}'][0] - Q).sum()) for r in R['ot']]
            dwM = [0.5 * float(np.abs(r[f'w{m}'][0] - r[f'w{M}'][0]).sum()) for r in R['ot']]
            cons = [float(np.abs(r[f'w{m}'] - r[f'w{m}'][0][None, :]).max()) for r in R['nearest']]
            wend = [float(np.abs(r[f'w{m}'][-1] - Q).max()) for r in R['ot']]
            seq = [OT.equalisation_time(r['s'], r[f'w{m}'], Q, m) for r in R['ot']]
            str_ = [transit_end_time(r['s'], r[f'transit{m}'], r[f'nsplit{m}'], m) for r in R['ot']] if set_name == 'commensurate' else []
            tr0 = [int(loop.transit_mask(start_clouds(r['seed'], M)[:m]).sum()) for r in R['ot']]
            tr1 = [int(loop.transit_mask(r[f'final{m}']).sum()) for r in R['ot']]
            nsp = [int(r[f'nsplit{m}'][-1]) for r in R['ot']]
            lines.append(f"| {m}{' (reference)' if m == M else ''} | {np.mean(w0):.4f} | {np.mean(dwM):.4f} | {max(cons):.1e} | {max(wend):.1e} | {', '.join(f'{x:.0f}' if np.isfinite(x) else '—' for x in seq)} | "
                         f"{', '.join(f'{x:.0f}' if np.isfinite(x) else '—' for x in str_) if str_ else 'n/a (split rows)'} | {', '.join(str(n) for n in tr0)} → {', '.join(str(n) for n in tr1)} | {', '.join(str(n) for n in nsp)} |")
        ln2 = np.log(2) / SIGMA0
        lines += ["", f"The pure leak carries a migrating particle across the Voronoi midpoint at s = ln 2/σ₀ = {ln2:.0f} slow units; the equalisation times above are where the last migrating particle's nearest atom flips.", ""]
        # fit errors
        lines += ["### Fit errors ‖f − f̂‖/‖f‖ at the end of the run", "", "| m | nearest-atom leak (mean over seeds; per seed) | OT leak (mean over seeds; per seed) | max whitened drift of the final configuration: nearest / OT (per seed) | migrated particles (nearest atom changed): count per seed; their distance to the new atom (G⁻¹; min–max over all) |", "|---|---|---|---|---|"]
        for m in list(ms) + [M]:
            fn = mean_se([float(r[f'fit{m}'][-1]) for r in R['nearest']]); fo = mean_se([float(r[f'fit{m}'][-1]) for r in R['ot']])
            dn = [float(np.abs(loop.drift(r[f'final{m}'], loop.nearest.leak_target(r[f'final{m}'])) / sqG[None, :]).max()) for r in R['nearest']]
            do = [float(np.abs(loop.drift(r[f'final{m}']) / sqG[None, :]).max()) for r in R['ot']]
            mig_n, mig_d = [], []
            for r in R['ot']:
                Th0 = start_clouds(r['seed'], M)[:m]; Th1 = r[f'final{m}']
                l0 = loop.nearest.assignment(Th0); l1 = loop.nearest.assignment(Th1)
                mig = np.flatnonzero(l0 != l1); mig_n.append(len(mig))
                if len(mig):
                    mig_d += list(np.sqrt(np.sum((Th1[mig] - exact[l1[mig]]) ** 2 / G[None, :], axis=1)))
            lines.append(f"| {m}{' (reference)' if m == M else ''} | {fmt(*fn, p=4)}; {', '.join(f'{float(r[f'fit{m}'][-1]):.4f}' for r in R['nearest'])} | {fmt(*fo, p=4)}; {', '.join(f'{float(r[f'fit{m}'][-1]):.1e}' for r in R['ot'])} | "
                         f"{', '.join(f'{x:.1e}' for x in dn)} / {', '.join(f'{x:.1e}' for x in do)} | {', '.join(str(n) for n in mig_n)}; {min(mig_d):.3f}–{max(mig_d):.3f} |" if mig_d else
                         f"| {m}{' (reference)' if m == M else ''} | {fmt(*fn, p=4)}; {', '.join(f'{float(r[f'fit{m}'][-1]):.4f}' for r in R['nearest'])} | {fmt(*fo, p=4)}; {', '.join(f'{float(r[f'fit{m}'][-1]):.1e}' for r in R['ot'])} | "
                         f"{', '.join(f'{x:.1e}' for x in dn)} / {', '.join(f'{x:.1e}' for x in do)} | {', '.join(str(n) for n in mig_n)}; — |")
        lines += ["", f"The pure leak's remaining approach of a migrated particle at s = {S_END:g}: Δ_{{kk'}}e^{{−σ₀·{S_END:g}}} = {', '.join(f'{np.sqrt(inter[k, kp]) * np.exp(-SIGMA0 * S_END):.3f}' for k, kp in ((0, 1), (1, 2), (0, 2)))} for the three atom pairs."]
        # the three distances on the last window
        a_w, b_w = WINDOWS[-1]
        lines += ["", f"### The distances to the reference on the last window ({a_w:g}, {b_w:g}] (means over the window and over seeds)", "",
                  "| m | ‖f̂_m − f̂_M‖²_{L²(ρ)}: nearest | OT | predicted (1/m − 1/M)tr(FC̃Fᵀ)/N (nearest, Proposition D‴.9(iii)) | W²_{G⁻¹}(μ_m, μ_M): nearest | OT | ½|w^m − w^M|₁·Δ_ν² (nearest: the cross-cluster transport of the weight mismatch, mean over seeds) | synchronous D: nearest | OT |",
                  "|---|---|---|---|---|---|---|---|---|"]
        H = {leak: {m: [window_mean(r['s'], r[f'h2{m}'], a_w, b_w) for r in R[leak]] for m in ms} for leak in R}
        W = {leak: {m: [window_mean(r['s'], r[f'W2{m}'], a_w, b_w) for r in R[leak]] for m in ms} for leak in R}
        D = {leak: {m: [window_mean(r['s'], r[f'D{m}'], a_w, b_w) for r in R[leak]] for m in ms} for leak in R}
        a_p, b_p = WINDOWS[-2]
        Hp = {leak: {m: [window_mean(r['s'], r[f'h2{m}'], a_p, b_p) for r in R[leak]] for m in ms} for leak in R}
        Wp = {leak: {m: [window_mean(r['s'], r[f'W2{m}'], a_p, b_p) for r in R[leak]] for m in ms} for leak in R}
        Dp = {leak: {m: [window_mean(r['s'], r[f'D{m}'], a_p, b_p) for r in R[leak]] for m in ms} for leak in R}
        cross_ratio, exact_w2, exact_h2, saddle_pairs, plateau_ratio = [], [], [], [], []
        for m in ms:
            cross = np.mean([0.5 * float(np.abs(r[f'w{m}'][-1] - r[f'w{M}'][-1]).sum()) * Delta ** 2 for r in R['nearest']])
            for r in R['nearest']:
                cross_ratio.append(window_mean(r['s'], r[f'W2{m}'], a_w, b_w) / (0.5 * float(np.abs(r[f'w{m}'][-1] - r[f'w{M}'][-1]).sum()) * Delta ** 2))
                wm, wM = r[f'w{m}'][-1], r[f'w{M}'][-1]
                Tm, _ = PD.weighted_stationary_atoms_prior(loop.nearest, wm, exact, w0=Q, n_cont=4); TM, _ = PD.weighted_stationary_atoms_prior(loop.nearest, wM, exact, w0=Q, n_cont=4)
                exact_w2.append(window_mean(r['s'], r[f'W2{m}'], a_w, b_w) / w2_weighted_atoms(wm, Tm, wM, TM))
                am, cm, sm, dm, pm = loop.features(Tm); aM, cM, sM, dM, pM = loop.features(TM)
                exact_h2.append(window_mean(r['s'], r[f'h2{m}'], a_w, b_w) / float(np.mean((pm @ (wm * am) - pM @ (wM * aM)) ** 2)))
                saddle_pairs.append(PD.splitting_margin_prior(loop.nearest, wm, exact)[0] <= 0)
            for r in R['ot']:                                                                  # the synchronous D's plateau: the pairs reassigned in one system and not in the other
                lm = loop.nearest.assignment(r[f'final{m}']); lM = loop.nearest.assignment(r[f'final{M}'][:m])
                plat = float(np.sum(inter[lm, lM]) / m)
                plateau_ratio.append(window_mean(r['s'], r[f'D{m}'], a_w, b_w) / plat if plat > 0 else float('nan'))
            lines.append(f"| {m} | {fmt(*mean_se(H['nearest'][m]))} | {fmt(*mean_se(H['ot'][m]))} | {PD.predicted_stationary_function_D(F, Q, m, M):.3e} | {fmt(*mean_se(W['nearest'][m]))} | {fmt(*mean_se(W['ot'][m]))} | {cross:.3e} | {fmt(*mean_se(D['nearest'][m]))} | {fmt(*mean_se(D['ot'][m]))} |")
        ratios_h = [np.mean(H['ot'][m]) / np.mean(H['nearest'][m]) for m in ms]; ratios_w = [np.mean(W['ot'][m]) / np.mean(W['nearest'][m]) for m in ms]
        lp = {leak: ([np.mean(H[leak][m]) / np.mean(Hp[leak][m]) for m in ms], [np.mean(W[leak][m]) / np.mean(Wp[leak][m]) for m in ms], [np.mean(D[leak][m]) / np.mean(Dp[leak][m]) for m in ms]) for leak in R}
        reg = [x for x, s_ in zip(exact_w2, saddle_pairs) if not s_]; sad = [x for x, s_ in zip(exact_w2, saddle_pairs) if s_]
        regh = [x for x, s_ in zip(exact_h2, saddle_pairs) if not s_]; sadh = [x for x, s_ in zip(exact_h2, saddle_pairs) if s_]
        lines += ["", f"Nearest-atom leak, exact laws at the realised weights (last-window means over the predictions from the weighted stationary points Θ*(w^m), Θ*(w^M)): W² against the K × K transport between the two weighted stationary measures, "
                  f"measured/exact {min(reg):.3f}–{max(reg):.3f} in the {len(reg)} pairs whose realised weights are admissible" + (f" and {', '.join(f'{x:.3f}' for x in sorted(sad))} at the {len(sad)} saddle weights (where the m-system is not at the replicated configuration)" if sad else "")
                  + f"; the formula ½|w^m − w^M|₁Δ_ν² is within {100 * (min(cross_ratio) - 1):+.0f}% to {100 * (max(cross_ratio) - 1):+.0f}% of the measured W² (measured/formula {min(cross_ratio):.3f}–{max(cross_ratio):.3f}). Function level against ‖f̂_{{Θ*(w^m)}} − f̂_{{Θ*(w^M)}}‖²_{{L²(ρ)}}: "
                  f"measured/exact {min(regh):.3f}–{max(regh):.3f} (admissible pairs)" + (f", {', '.join(f'{x:.3f}' for x in sorted(sadh))} at the saddles" if sadh else "") + f"; the w̄-law of the table is the leading term at the mean weights."]
        if plateau_ratio:
            lines += ["", f"OT leak, the synchronous D: its plateau is the mass reassigned in one system and not in the other, Σ_{{i: label_m(i) ≠ label_M(i)}}Δ²_{{kk'}}/m from the end labels; measured last-window D/plateau = {np.nanmin(plateau_ratio):.3f}–{np.nanmax(plateau_ratio):.3f} over the {int(np.sum(np.isfinite(plateau_ratio)))} pairs with mismatched particles — " + ("D converges to the plateau from below (the migrated particles still approaching their atoms) and does not vanish." if np.nanmax(plateau_ratio) < 1 else "D does not vanish; on these sizes the split (bridge) particles of the two systems add to and subtract from the plateau, so the ratio is not one-sided.")]
        sl = {leak: {k: slope_se(X[leak]) for k, X in (('h2', H), ('W2', W), ('D', D))} for leak in R}
        lines += ["", f"Ratios OT/nearest of the function-level distance: {', '.join(f'{x:.1e}' for x in ratios_h)}; of W²: {', '.join(f'{x:.1e}' for x in ratios_w)}. "
                  f"Ratios of the last window to the one before (({a_p:g}, {b_p:g}]; 1 = stationary): function level nearest {', '.join(f'{x:.2f}' for x in lp['nearest'][0])}, OT {', '.join(f'{x:.3f}' for x in lp['ot'][0])}; "
                  f"W² nearest {', '.join(f'{x:.2f}' for x in lp['nearest'][1])}, OT {', '.join(f'{x:.3f}' for x in lp['ot'][1])}; synchronous D nearest {', '.join(f'{x:.2f}' for x in lp['nearest'][2])}, OT {', '.join(f'{x:.3f}' for x in lp['ot'][2])}. "
                  f"Slopes in m (log–log, per seed, mean ± s.e.): function level nearest {sl['nearest']['h2'][0]:.2f} ± {sl['nearest']['h2'][1]:.2f}, OT {sl['ot']['h2'][0]:.2f} ± {sl['ot']['h2'][1]:.2f}; "
                  f"W² nearest {sl['nearest']['W2'][0]:.2f} ± {sl['nearest']['W2'][1]:.2f}, OT {sl['ot']['W2'][0]:.2f} ± {sl['ot']['W2'][1]:.2f}; synchronous D nearest {sl['nearest']['D'][0]:.2f} ± {sl['nearest']['D'][1]:.2f}, OT {sl['ot']['D'][0]:.2f} ± {sl['ot']['D'][1]:.2f}."]
        # time course of the OT distances: window means
        lines += ["", "### Time course under the OT leak: window means of the function-level distance² and of W² (the transit, then the stationary state)", "",
                  "| m | " + " | ".join(f"‖f̂_m − f̂_M‖² on ({a:g}, {b:g}]" for a, b in WINDOWS) + " | " + " | ".join(f"W² on ({a:g}, {b:g}]" for a, b in WINDOWS) + " |", "|---|" + "---|" * (2 * len(WINDOWS))]
        for m in ms:
            hw = [np.mean([window_mean(r['s'], r[f'h2{m}'], a, b) for r in R['ot']]) for a, b in WINDOWS]
            ww = [np.mean([window_mean(r['s'], r[f'W2{m}'], a, b) for r in R['ot']]) for a, b in WINDOWS]
            lines.append(f"| {m} | " + " | ".join(f"{x:.3e}" for x in hw) + " | " + " | ".join(f"{x:.3e}" for x in ww) + " |")
        if set_name == 'incommensurate':
            lines += ["", "### The quantisation laws of (iv) against the mean-field stationary measure μ*(q) = ν (OT leak; the reference carries its own bridge particles)", "",
                      "| m | W²(μ_m, ν) on the last window (mean over seeds) | m·W²(μ_m, ν) | fit error ‖f̂_m − f‖/‖f‖ at the end (mean) | m·fit |", "|---|---|---|---|---|"]
            Wn = {m: [window_mean(r['s'], r[f'W2nu{m}'], a_w, b_w) for r in R['ot']] for m in list(ms) + [M]}
            Fn = {m: [float(r[f'fit{m}'][-1]) for r in R['ot']] for m in list(ms) + [M]}
            for m in list(ms) + [M]:
                lines.append(f"| {m}{' (reference)' if m == M else ''} | {fmt(*mean_se(Wn[m]))} | {m * np.mean(Wn[m]):.1f} | {fmt(*mean_se(Fn[m]), p=4)} | {m * np.mean(Fn[m]):.2f} |")
            slW = slope_se({m: Wn[m] for m in ms}); slWM = slope_se(Wn); slF = slope_se({m: [x ** 2 for x in Fn[m]] for m in ms}); slFM = slope_se({m: [x ** 2 for x in Fn[m]] for m in Fn})
            lines += ["", f"Slopes in m: W²(μ_m, ν) {slW[0]:.2f} ± {slW[1]:.2f} over m = {', '.join(str(m) for m in ms)} ({slWM[0]:.2f} ± {slWM[1]:.2f} with M = {M}); fit error squared {slF[0]:.2f} ± {slF[1]:.2f} ({slFM[0]:.2f} ± {slFM[1]:.2f} with M): the order 1/m of (iv) in W² and 1/m² at function level, with the constants m·W² ≈ {np.mean([m * np.mean(Wn[m]) for m in Wn]):.0f} = {np.mean([m * np.mean(Wn[m]) for m in Wn]) / Delta ** 2:.2f}Δ_ν² and m·fit ≈ {np.mean([m * np.mean(Fn[m]) for m in Fn]):.2f}."]
            lines += ["", "### The split particles at the end of the OT runs (plan rows m·γ_i; distance to the nearest prior atom in units of Δ_ν; distance to the barycentric target T̄_i)", ""]
            for r in R['ot']:
                for m in list(ms) + [M]:
                    for b in r[f'bridge{m}']:
                        lines.append(f"- seed {r['seed']}, m = {m}: particle {b['i']}, weights ({', '.join(f'{x:.3f}' for x in b['weights'])}), dist/Δ_ν = {b['dist_nearest_over_Delta']:.3f}, ‖ϑ − T̄‖ = {b['dist_target']:.3f}")
            bound = [(loop.K0 - 1) / m * Delta ** 2 for m in ms]
            lines += ["", f"Forest bound on the split mass: (K₀ − 1)/m = {', '.join(f'{2 / m:.4f}' for m in ms)} for m = {', '.join(str(m) for m in ms)}; its worst-case W² contribution (K₀ − 1)Δ_ν²/m = {', '.join(f'{x:.2f}' for x in bound)}."]
        lines.append("")
    # ---- the masses as a design input
    lines += ["## 3. The masses as a design input: the weighted stationary point Θ*(q_δ) for q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) (no dynamics)", "",
              "| δ | q_δ | fit error of Θ*(q_δ) (the bias of the mass constraint) | atoms (a, c, s) | κ_R at the atoms | κ_R over N_ϱ′ | margin κ = σ₀ − κ_R(N_ϱ′) | splitting margin σ₀ − max κ_R(atoms) (admissibility of q_δ) | ε = max‖ϑ*_k − ϑ⁰_k‖ against Δ_ν/2 (diagonal plan) |",
              "|---|---|---|---|---|---|---|---|---|"]
    for d in design:
        lines.append(f"| {d['delta']:g} | ({', '.join(f'{x:.4f}' for x in d['q'])}) | {d['fit']:.4f}" + (f" ({d['fit'] / d['delta']:.2f}δ)" if d['delta'] > 0 else "") + " | " + ", ".join(f"({d['a'][k]:.3f}, {d['c'][k]:.3f}, {d['s'][k]:.3f})" for k in range(3))
                     + f" | {', '.join(f'{x:+.5f}' for x in d['k_at'])} | {d['k_nb']:+.5f} | {d['margin']:+.5f} | {d['split_margin']:+.5f} | {d['eps']:.3f} < {d['half']:.3f}: {'yes' if d['eps'] < d['half'] else 'no'} |")
    bad = [i for i, d in enumerate(design) if d['split_margin'] <= 0]
    where = (f"between δ = {design[bad[0] - 1]['delta']:g} and {design[bad[0]]['delta']:g}" if bad and bad[0] > 0 else (f"below δ = {design[0]['delta']:g}" if bad else f"beyond δ = {design[-1]['delta']:g}"))
    lines += ["", f"Reading: the stationary measure of the OT leak carries the prior's masses, so a mass mismatch is a bias (fit error {design[0]['fit']:.1e} at δ = 0, {design[-1]['fit']:.3f} at δ = {design[-1]['delta']:g}) "
              f"and, from a threshold {where}, an instability: the over-weighted first cluster overshoots the plant, its splitting margin turns negative and the equal-mass "
              f"clusters split (Remark D‴.11); the designed masses must be admissible — the same asymmetry as the matching prior's admissible set in `prior_design.md`. The diagonal plan is the unique optimal plan at every row (ε < Δ_ν/2).", ""]
    lines += ["Numbers are produced by `experiments/exp_ot_leak.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'ot_leak.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'ot_leak_raw.npz'), allow_pickle=True)
        report(list(z['runs']), list(z['design']), t0)
        return
    design = [mass_design_row(d) for d in DELTAS]
    jobs = [(seed, set_name, leak) for seed in range(N_SEEDS) for set_name in SETS for leak in ('ot', 'nearest')]
    with mp.Pool(N_WORKERS) as pool:
        runs = pool.map(coupled_job, jobs, chunksize=1)
    np.savez(os.path.join(OUT, 'ot_leak_raw.npz'), runs=np.array(runs, dtype=object), design=np.array(design, dtype=object), allow_pickle=True)
    report(runs, design, t0)


if __name__ == '__main__':
    main()
