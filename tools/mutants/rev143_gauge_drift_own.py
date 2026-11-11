# (id, what, old, new, k): single-point mutants of exp_gauge_drift.py — round 143's revision (referee D3el-1, -4, -5), never mutation-tested
# when applied: Part I (`normal_rates`: the slowest normal rates and the slowest mode's shares of its g-norm; `part_I`: the leak's stiffness
# and spread against the distance to the faces, the rates toward the first face, the Péclet number and the free spread, F's ensemble against
# the limit SDE's; `report_I`), with `kappa` and `sde_paths`, which Part I calls. Part I's numbers are stored in the raw, so the md's
# regeneration sees only `report_I`.
# Run: tools/mutants.py --exp exp_gauge_drift --tests tests/test_gauge_drift.py --raw results/chaos/gauge_drift_raw.npz
# --md results/chaos/gauge_drift.md --out-env GDR_OUT --workers-env GDR_WORKERS (one worker; the whole test file, ~2–8 s per mutant).
# Result: 25 mutants and a control. Against the tests as committed: 11 killed by the tests (D01–D03, D06, D07, D13, D14, D17–D19, D23), 2 by
# the md (D20: the SDE's seed is printed in the header; D25); the 12 others (D04, D05, D08–D12, D15, D16, D21, D22, D24) killed by the tests
# added then: test_part_I_is_recomputed_from_the_raw (every stored number of Part I recomputed from the stored field E and ensemble F, to 1e-9),
# test_normal_rates_shares_on_a_planted_linearisation, and a planted stiffness in test_part_I_on_the_stored_field_and_ensemble (D09: on the
# instance numpy returns the 2×2 stiffness's eigenvalues already in ascending order, so only a plant whose eigvals come unsorted sees the
# missing sort). No survivor; the control passes.
M = [
 ('D00', 'control', "def normal_rates(w, gamma, n=4):", "def normal_rates(w, gamma, n=4):", 0),
 ('D01', "the slowest normal mode one further", "v = V[:, 2].real", "v = V[:, 3].real", 0),
 ('D02', "the rates from the tangent's second eigenvalue", "rates=-lam[2:2 + n].real", "rates=-lam[1:1 + n].real", 0),
 ('D03', "the eigenvalues sorted by their real part", "o = np.argsort(np.abs(lam))", "o = np.argsort(lam.real)", 0),
 ('D04', "δβ without the projection onto Σδβ = 0", "db = beta * (dz - W @ dz)", "db = beta * dz", 0),
 ('D05', "the masses' g-norm without γ", "mas = np.sum(db ** 2 / (gamma * beta))", "mas = np.sum(db ** 2 / beta)", 0),
 ('D06', "the positions' g-norm with G, not G⁻¹", "(v[:9].reshape(3, 3) ** 2 / E.G[None, :]).sum(axis=0)", "(v[:9].reshape(3, 3) ** 2 * E.G[None, :]).sum(axis=0)", 0),
 ('D07', "the distance to the faces W − W_min²", "w = np.ones(3) / 3; dist = w - W_MIN", "w = np.ones(3) / 3; dist = w - W_MIN ** 2", 0),
 ('D08', "the leak's spread over W, not over the distance", "rel=sp / dist))", "rel=sp / w))", 0),
 ('D09', "the stiffness's eigenvalues unsorted", "eig=np.sort(np.linalg.eigvals(kappa(w, eps)).real)", "eig=np.linalg.eigvals(kappa(w, eps)).real", 0),
 ('D10', "the line toward the face W_2 = 2W_3", "wf = np.array([W1, (1 - W1) / 2])", "wf = np.array([W1, (1 - W1) / 3])", 0),
 ('D11', "α_1 read as the first centre", "r['alpha1'] = float(instance(g)[1](wf)[0])", "r['alpha1'] = float(instance(g)[1](wf)[1])", 0),
 ('D12', "the face rates at the instance's γ", "r = normal_rates(wf, g); r['delta'] = d", "r = normal_rates(wf, GAMMA); r['delta'] = d", 0),
 ('D13', "Pe over D_22", "pe=abs(b[0]) * dist[0] / D[0, 0]", "pe=abs(b[0]) * dist[0] / D[1, 1]", 0),
 ('D14', "the free spread at the time 1/τ", "spread=np.sqrt(D[0, 0] / (g * TAU))", "spread=np.sqrt(D[0, 0] / TAU)", 0),
 ('D15', "Pe's drift at the instance's γ", "b, D, _ = drift_at(w[:2], g)", "b, D, _ = drift_at(w[:2])", 0),
 ('D16', "the SDE's record farthest from s", "t = min(T_I, key=lambda u: abs(u - s_))", "t = max(T_I, key=lambda u: abs(u - s_))", 0),
 ('D17', "the SDE's s.e. over n, not √n", "ys = Y.std(axis=0, ddof=1) / np.sqrt(len(Y))", "ys = Y.std(axis=0, ddof=1) / len(Y)", 0),
 ('D18', "the combined distance with se + ys", "z=(m_ - Y.mean(axis=0)) / np.hypot(se, ys)", "z=(m_ - Y.mean(axis=0)) / (se + ys)", 0),
 ('D19', "the loop's mean drift over sτ", "rate_loop=(m_ - 1 / 3) / (s_ * GAMMA * TAU)", "rate_loop=(m_ - 1 / 3) / (s_ * TAU)", 0),
 ('D20', "the SDE's own seed 2", "T_I, SEED_I = (1e4, 2e4, 3e4, 5e4), 1", "T_I, SEED_I = (1e4, 2e4, 3e4, 5e4), 2", 0),
 ('D21', "Katzenberger's drift at 1/3 over τ", "out['b0'] = np.r_[b0, -b0.sum()] / (GAMMA * TAU)", "out['b0'] = np.r_[b0, -b0.sum()] / TAU", 0),
 ('D22', "the face distances to 10⁻³", "(0.08, 0.02, 0.005, 0.002, 0.0005)", "(0.08, 0.02, 0.005, 0.002, 0.001)", 0),
 ('D23', "the leak's metric with W, not 1/W", "g = np.diag(np.r_[np.full(K, 1 / Ga), (m / GAMMA) / W])", "g = np.diag(np.r_[np.full(K, 1 / Ga), (m / GAMMA) * W])", 0),
 ('D24', "the SDE's increments with ds, not √ds", "z = rng.standard_normal((N_PATHS, 2)) * np.sqrt(DS_SDE)", "z = rng.standard_normal((N_PATHS, 2)) * DS_SDE", 0),
 ('D25', "the md's face rates two, not three", "{', '.join(fmt(v) for v in r['rates'][:3])} |\")", "{', '.join(fmt(v) for v in r['rates'][:2])} |\")", 0),
]
