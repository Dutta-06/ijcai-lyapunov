# (id, what, old, new, k): single-point mutants of exp_mirror_leak.py — round 141's revision (referee D3ej), never mutation-tested when
# applied: Part E, the 1760 patterns of D(2) against the exact weighted LASSO (`lasso_exact`: coordinate descent, the active-set
# equations, the sign flips, the largest KKT violator, the KKT residual and slacks; `tangent_profile`; `pattern_record`; `part_E`) and
# `report_E`.
# Result: 24 mutants and a control. Against the tests as committed: 4 killed by the tests (E02, E04, E09, E11), 4 by the tests and the md
# (E21–E24); 14 killed by the tests added then (test_lasso_exact_corrects_any_warm_start, the kkt ≥ 0 check, test_part_E_replays_the_stored_
# patterns, test_pattern_record_prediction_and_threshold_support_on_planted_inputs — E10, E12 change no stored pattern). Two survivors,
# equivalent: E01 (coordinate descent's warm start is only a start: the active-set loop ends on the unique minimiser's support and signs and
# solves the same equations there, bit for bit, as the three planted starts show); E08 (the KKT residual's max with 0 off S is absorbed by
# the outer max with the support term, which is ≥ 0, or 0.0 when the support is empty).
M = [
 ('E00', 'control', "def tangent_profile(L, b0, slack, sgn):", "def tangent_profile(L, b0, slack, sgn):", 0),
 ('E01', "coordinate descent's warm start at the weight 2μ", "w = MF.lasso(H, ws, np.full(n, mu), ke)", "w = MF.lasso(H, ws, np.full(n, 2 * mu), ke)", 0),
 ('E02', "the active-set equations with +k_eμs", "(H @ ws)[Sl] - ke * mu * sg[Sl])", "(H @ ws)[Sl] + ke * mu * sg[Sl])", 0),
 ('E03', "no sign flips checked", "flip = [j for j in Sl if np.sign(w[j]) != sg[j]]", "flip = []", 0),
 ('E04', "the KKT multipliers over μ (k_e dropped)", "c = H @ (ws - w) / (ke * mu)                                              # the KKT",
  "c = H @ (ws - w) / mu                                              # the KKT", 0),
 ('E05', "a violator enters only beyond 1 + 10⁻³", "np.abs(c[off]).max() > 1 + 1e-12", "np.abs(c[off]).max() > 1 + 1e-3", 0),
 ('E06', "the smallest violator enters", "j = int(off[np.argmax(np.abs(c[off]))])", "j = int(off[np.argmin(np.abs(c[off]))])", 0),
 ('E07', "the entrant's sign reversed", "S.add(j); sg[j] = np.sign(c[j])", "S.add(j); sg[j] = -np.sign(c[j])", 0),
 ('E08', "the KKT residual off S signed (no max with 0)", "float(max(np.abs(c[w == 0]).max() - 1, 0.0))", "float(np.abs(c[w == 0]).max() - 1)", 0),
 ('E09', "the tangent profile with sin²", "np.cos(x) / np.sin(x) ** 3", "np.cos(x) / np.sin(x) ** 2", 0),
 ('E10', "the prediction's shift with +μ", "wS = wt[S_] - mu * np.linalg.solve(H[np.ix_(S_, S_)], sg_)", "wS = wt[S_] + mu * np.linalg.solve(H[np.ix_(S_, S_)], sg_)", 0),
 ('E11', "the prediction without the irrepresentable condition", "pred = bool(ic < 1 and np.all(np.sign(wS) == sg_))", "pred = bool(np.all(np.sign(wS) == sg_))", 1),
 ('E12', "the LASSO's threshold relative to the equilibrium's largest", "np.abs(wl) > thr * np.abs(wl).max()", "np.abs(wl) > thr * np.abs(w).max()", 0),
 ('E13', "the tangent profile's sign reversed", "tp = tangent_profile(L, b0, slack[inact], np.sign(c[inact]))", "tp = tangent_profile(L, b0, slack[inact], -np.sign(c[inact]))", 0),
 ('E14', "the profile's error relative to the equilibrium", "tp_err = np.abs(w[inact] - tp) / np.abs(tp)", "tp_err = np.abs(w[inact] - tp) / np.abs(w[inact])", 0),
 ('E15', "atoms beyond the LASSO counted the other way", "eq_extra_vs_l0=int(len(sup_eq - sup_l0))", "eq_extra_vs_l0=int(len(sup_l0 - sup_eq))", 0),
 ('E16', "the smallest LASSO coefficient absolute", "l_small=float(np.abs(wl[wl != 0]).min() / np.abs(wl).max())", "l_small=float(np.abs(wl[wl != 0]).min())", 0),
 ('E17', "E's λ₁ without the square root", "lam1 = np.pi * L.A * np.sqrt(abs(2 * L.c - b0) / b0) / (4 * L.G * L.c)",
  "lam1 = np.pi * L.A * abs(2 * L.c - b0) / b0 / (4 * L.G * L.c)", 0),
 ('E18', "E's magnitudes from seed 24", "rng = np.random.default_rng(23)", "rng = np.random.default_rng(24)", 1),
 ('E19', "E's magnitudes in [0.1, 0.3]", "mag = rng.uniform(0.15, 0.3, 3)", "mag = rng.uniform(0.1, 0.3, 3)", 1),
 ('E20', "the core's prefactor with β⁰^{1/2}", "core=float(L.A * b0 ** 1.5 / (L.m * np.sqrt(2 * L.c)))", "core=float(L.A * b0 ** 0.5 / (L.m * np.sqrt(2 * L.c)))", 0),
 ('E21', "the disagreements against the LASSO's support", "mm = [p for p in P if p['pred'] != p['eq_exact']]", "mm = [p for p in P if p['pred'] != p['lt_exact']]", 0),
 ('E22', "the threshold support's agreement counted with the exact one", "cnt(lambda p: p['lt_exact'] == p['pred'])", "cnt(lambda p: p['l0_exact'] == p['pred'])", 0),
 ('E23', "'atoms it drops' counted as kept", "LASSO atoms it drops {sum(p['eq_missing_vs_l0'] for p in eqd)}", "LASSO atoms it drops {sum(p['eq_extra_vs_l0'] for p in eqd)}", 0),
 ('E24', "the tangent profile's largest error over the medians", "largest {fmt(max(p['tp_err_max'] for p in P if np.isfinite(p['tp_err_max'])))}",
  "largest {fmt(max(p['tp_err_med'] for p in P if np.isfinite(p['tp_err_med'])))}", 0),
]
