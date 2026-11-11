# (id, what, old, new, k): single-point mutants of the round-150 revision (referee D3es): exp_gauge_fixman.py Part K (`cut_K`, `local_K`,
# `kat_density`, `law_local_K`, `part_K`, `report_K`: the seed-level jackknife of Part J's pooled exponents, the κ_F = 2 − 1 difference on
# the paired seeds, the limit equation's own window exponent at γ = 0.1, the cap over the threshold) — ids K..; and exp_gauge_plateau.py
# Part H (`rows_G`, `fixed_barrier`, `extrap`, `boundary_minima`, `exit_points`, `exits_by_tau`, `part_H`, `report_H`) — ids H.. .
# Run each half with --only (an id whose anchor is in the other script is not attempted):
#   tools/mutants.py --exp exp_gauge_fixman --tests tests/test_gauge_fixman.py --raw results/chaos/gauge_fixman_raw.npz
#     --md results/chaos/gauge_fixman.md --out-env GF_OUT --workers-env GF_WORKERS --extra results/chaos/gauge_drift_raw.npz
#     --spec tools/mutants/rev150_gauge_own.py --work <scratch> --slot k --only K00,K01,…
#   tools/mutants.py --exp exp_gauge_plateau --tests tests/test_gauge_plateau.py --raw results/chaos/gauge_plateau_raw.npz
#     --md results/chaos/gauge_plateau.md --out-env GP_OUT --workers-env GP_WORKERS --spec tools/mutants/rev150_gauge_own.py
#     --work <scratch> --slot h --only H00,H01,…
# Parts K and H are stored in the raw, so the md's regeneration sees only report_K and report_H; the tests recompute them.
# Result (on the final scripts): 25 mutants and two controls (K00, H00 pass the tests and leave the md unchanged). All 25 killed by the
# tests, K12 also by the md. H09 (the minima clustered by energy only) survived a first run, being equivalent on the data (the two minima
# differ in energy); it is killed by test_cluster_minima_on_plants, added then (clustering factored out as cluster_minima). No survivor.
K_IDS = [f'K{i:02d}' for i in range(13)]
H_IDS = [f'H{i:02d}' for i in range(14)]
M = [
 ('K00', 'control', "def cut_K(k, dl):", "def cut_K(k, dl):", 0),
 ('K01', "the cut's other masses not equal", "Wk = W_MIN[k] + dl; W = np.full(3, (1 - Wk) / 2); W[k] = Wk",
  "Wk = W_MIN[k] + dl; W = np.array([0.4, 0.35, 0.25]) * (1 - Wk) / 0.75; W[k] = Wk", 0),
 ('K02', "the cut's direction moves one other mass", "v = np.full(3, -0.5); v[k] = 1.0", "v = np.zeros(3); v[(k + 1) % 3] = -1.0; v[k] = 1.0", 0),
 ('K03', "the zero-flux exponent without the factor 2", "return float(dl * 2 * (a @ b) / (a @ D @ a)", "return float(dl * (a @ b) / (a @ D @ a)", 0),
 ('K04', "d log D with the wrong sign", "- dl * (np.log(a @ Dp @ a) - np.log(a @ Dm @ a)) / (2 * hh))",
  "+ dl * (np.log(a @ Dp @ a) - np.log(a @ Dm @ a)) / (2 * hh))", 0),
 ('K05', "the density integrated in δ, not log δ", "lp = np.r_[0.0, np.cumsum(0.5 * (ex[1:] + ex[:-1]) * np.diff(u))]",
  "lp = np.r_[0.0, np.cumsum(0.5 * (ex[1:] + ex[:-1]) * np.diff(grid))]", 0),
 ('K06', "the bin average without the Jacobian δ", "return float(np.sum(np.exp(np.interp(vm, u, lp) + vm)) * (v[1] - v[0]) / (hi - lo))",
  "return float(np.sum(np.exp(np.interp(vm, u, lp))) * (v[1] - v[0]) / (hi - lo))", 0),
 ('K07', "the layer law's local exponent without its middle factor", "return (kF - 1) + (kF - 1) * dl / (dl + 2 * W_MIN[k]) - kF * dl / (W_MIN[k] + dl)",
  "return (kF - 1) - kF * dl / (W_MIN[k] + dl)", 0),
 ('K08', "the jackknife without (S − 1)", "jse = lambda v: float(np.sqrt((S - 1) / S * np.sum((v - v.mean()) ** 2)))",
  "jse = lambda v: float(np.sqrt(1 / S * np.sum((v - v.mean()) ** 2)))", 0),
 ('K09', "the pooled weights equal", "w = np.array([1 / f['se'] ** 2 for f in ff]); w = w / w.sum()", "w = np.ones(len(ff)); w = w / w.sum()", 0),
 ('K10', "the difference's jackknife unpaired (sum)", "diff = dict(d=e[2.0] - e[1.0], se=jse(jj[2.0] - jj[1.0])", "diff = dict(d=e[2.0] - e[1.0], se=jse(jj[2.0] + jj[1.0])", 0),
 ('K11', "the ratio threshold over cap", "ratio=np.asarray(r['cap']) / np.asarray(r['thr'])", "ratio=np.asarray(r['thr']) / np.asarray(r['cap'])", 0),
 ('K12', "the md's jackknife s.e. with two decimals", "| {p['se_ind']:.3f} | {p['se']:.3f} |", "| {p['se_ind']:.3f} | {p['se']:.2f} |", 0),
 ('H00', 'control', "def rows_G(runs, t0=0.0):", "def rows_G(runs, t0=0.0):", 0),
 ('H01', "the lag keeps the exits before t₀", "n = int(np.sum(ex & (t > t0)))", "n = int(np.sum(ex))", 0),
 ('H02', "the lag does not shorten the runs", "T = float(np.sum(np.maximum(t - t0, 0.0)))", "T = float(np.sum(t))", 0),
 ('H03', "the fixed barrier added, not removed", "y = np.log([r['rate'] for r in ok]) + dJ / tau", "y = np.log([r['rate'] for r in ok]) - dJ / tau", 0),
 ('H04', "the fit unweighted", "tau = np.array([r['tau'] for r in ok]); lt = np.log(tau); w = np.array([r['n'] for r in ok], float)",
  "tau = np.array([r['tau'] for r in ok]); lt = np.log(tau); w = np.ones(len(ok))", 0),
 ('H05', "Δχ² = 4 at three s.e.", "int2=(c - 2 * se, c + 2 * se)", "int2=(c - 3 * se, c + 3 * se)", 0),
 ('H06', "the fixed prefactors' χ² on n − 2 dof", "fixed = [dict(c=cc, chi2=chi2(cc), dof=len(ok) - 1, p=float(C2.sf(chi2(cc), len(ok) - 1))",
  "fixed = [dict(c=cc, chi2=chi2(cc), dof=len(ok) - 1, p=float(C2.sf(chi2(cc), len(ok) - 2))", 0),
 ('H07', "the extrapolation at 1/τ = 10³", "cov = np.linalg.inv(X.T @ (w[:, None] * X)); x4 = np.r_[1.0, -1e4]", "cov = np.linalg.inv(X.T @ (w[:, None] * X)); x4 = np.r_[1.0, -1e3]", 0),
 ('H08', "the barrier's part without 10⁴", "lseB=float(1e4 * f['B_se'])", "lseB=float(f['B_se'])", 0),
 ('H09', "the minima clustered by energy only", "if abs(vals[i] - c['dJ']) < CLUSTER_H[0] and np.abs(U[i] - c['u']).max() < CLUSTER_H[1]:",
  "if abs(vals[i] - c['dJ']) < CLUSTER_H[0]:", 0),
 ('H10', "the flow's largest φ includes the start", "T = np.array([c['Th'] for c in cl]); lb = np.zeros((len(cl), 3)); phi = np.zeros(len(cl))",
  "T = np.array([c['Th'] for c in cl]); lb = np.zeros((len(cl), 3)); phi = dom_phi(T, Th0)", 0),
 ('H11', "the exit point not checked against the stored time", "assert (None if s['exit'] is None else s['exit'] * DT) == r['t'], (r['tau'], r['seed'])", "pass", 0),
 ('H12', "broadened and narrowed swapped", "broad=int(np.sum(w2 & (U[:, 1, 2] > 0)))", "broad=int(np.sum(w2 & (U[:, 1, 2] < 0)))", 0),
 ('H13', "the exits before the lag's ½, not 1", "and x['t'] < LAGS_H[1]))) for r in rr]", "and x['t'] < LAGS_H[0]))) for r in rr]", 0),
]
