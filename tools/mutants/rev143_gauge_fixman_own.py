# (id, what, old, new, k): single-point mutants of exp_gauge_fixman.py — round 143's revision (referee D3el-3, -10) and the confinement runs
# after it, never mutation-tested: Part I (`gamma_kk`, `part_I`: the remainder at smaller γ, Γ_kk, the floor's cap, the finite-τ threshold,
# the stored records' fractions; `report_I`) and Part J (`floor_J`, `thr_J`, `layer_J`, `record_J`, `j_runs`, `window_bins`, `expo_fit`,
# `layer_density`, `expo_jack`, `verdict_J`, `pool_J`, `confined_J`, `part_J`, `report_J`). Parts I and J are stored in the raw (J's runs
# as their accumulators), so the md's regeneration sees only the reports.
# Run: tools/mutants.py --exp exp_gauge_fixman --tests tests/test_gauge_fixman.py --raw results/chaos/gauge_fixman_raw.npz
# --md results/chaos/gauge_fixman.md --out-env GF_OUT --workers-env GF_WORKERS --extra results/chaos/gauge_drift_raw.npz (one worker; ~4–10 s).
# Result: 25 mutants and a control. Against the tests as committed: 11 killed by the tests (F02, F03, F05–F09, F15, F19–F21), F25 by the tests
# and the md; the 13 others (F01, F04, F10–F14, F16–F18, F22–F24) killed by the tests added then: test_part_I_is_recomputed_from_the_raw,
# test_part_J_is_recomputed_in_full_from_the_stored_accumulators (every number and label, and the stored accumulators' set-up),
# test_j_runs_set_up_and_a_short_run_against_run_dithered (j_runs' own batch for 2000 time units against run_dithered and the stored first
# rows), the in_cap assertion in test_part_J_record_accumulator (F10), and test_j_runs_hands_the_batch_each_runs_parameters (F12: the stored
# runs never come within the cap, so a fixed floor 0.006 at τ = 5·10⁻⁵ changes no record; only the batch's parameters show it). No survivor.
M = [
 ('F00', 'control', "def floor_J(tau):", "def floor_J(tau):", 0),
 # Part I
 ('F01', "I's remainder at γ = 4·10⁻⁵", "GAMMAS_I, GAINS_I, FLOOR_PLUG = (5e-5,", "GAMMAS_I, GAINS_I, FLOOR_PLUG = (4e-5,", 0),
 ('F02', "Γ_kk = ⟨φ_k⟩², not ⟨φ_k²⟩", "return np.mean(phi ** 2, axis=0)", "return np.mean(phi, axis=0) ** 2", 0),
 ('F03', "the threshold with m = 2", "thr = {k: np.sqrt(k * TAU / (3 * A_MAX ** 2 * Gam)) for k in GAINS_I}",
  "thr = {k: np.sqrt(k * TAU / (2 * A_MAX ** 2 * Gam)) for k in GAINS_I}", 0),
 ('F04', "a′ = 1 − tanh(α/a_max)", "apr = 1 - np.tanh(R[:, :, 1:4] / A_MAX) ** 2", "apr = 1 - np.tanh(R[:, :, 1:4] / A_MAX)", 0),
 # Part J's constants and closed forms
 ('F05', "the floor linear in τ (not τ^{1/2})", "return FLOOR_J * np.sqrt(tau / 1e-4)", "return FLOOR_J * (tau / 1e-4)", 0),
 ('F06', "the threshold at κ_F = 0 zero (not κ_F = 1's)", "return np.sqrt(max(kF, 1.0) * tau / (3 * A_MAX ** 2 * np.asarray(Gam)))",
  "return np.sqrt(kF * tau / (3 * A_MAX ** 2 * np.asarray(Gam)) + 1e-300)", 0),
 ('F07', "the layer without the ½", "return 0.5 * np.sqrt(gamma * A_MAX * np.abs(E.A_BUMP) / (3 * E.G[0]))",
  "return np.sqrt(gamma * A_MAX * np.abs(E.A_BUMP) / (3 * E.G[0]))", 0),
 # record_J and j_runs
 ('F08', "the time segment from k/(n + 1)", "g = min(k * N_SEG // acc['n'], N_SEG - 1)", "g = min(k * N_SEG // (acc['n'] + 4 * REC_J), N_SEG - 1)", 0),
 ('F09', "the masses W = β, not β/3", "W = np.exp(lb) / 3; d = W - W_MIN", "W = np.exp(lb); d = W - W_MIN", 0),
 ('F10', "inside δ_cap with ≤ the threshold", "acc['in_cap'][g] += d < acc['cap']", "acc['in_cap'][g] += d < acc['thr']", 0),
 ('F11', "the run's noise √(τ dt G)", "sq=np.sqrt(2 * tau * DT * E.G), Th0=Th", "sq=np.sqrt(tau * DT * E.G), Th0=Th", 0),
 ('F12', "the batch's floor the fixed 0.006", "par.append((kF, tau, floor_J(tau)))", "par.append((kF, tau, FLOOR_J))", 0),
 ('F13', "the runs' cap at the floor 0.01", "cap = np.array([delta_cap(k, floor_J(tau)) for k in range(3)])", "cap = np.array([delta_cap(k) for k in range(3)])", 0),
 ('F14', "the batch's Â the plug-in form", "fixman=1.0, fixman_tau=1e-4, fixman_A=A0)", "fixman=1.0, fixman_tau=1e-4, fixman_A=None)", 0),
 # the analysis
 ('F15', "the window's first bin one lower", "i0 = int(np.searchsorted(EDGES_J, lo, side='left')) + 1", "i0 = int(np.searchsorted(EDGES_J, lo, side='left'))", 0),
 ('F16', "the fit unweighted", "x = np.log(np.sqrt(lo * hi)); w = c", "x = np.log(np.sqrt(lo * hi)); w = np.ones_like(c)", 0),
 ('F17', "the layer law's W_k^{−κ_F} at 2W_min", "* (W_MIN[k] + d) ** (-kF)", "* (2 * W_MIN[k] + d) ** (-kF)", 0),
 ('F18', "the jackknife without (S − 1)", "se = float(np.sqrt((S - 1) / S * np.sum((ej - ej.mean()) ** 2)))", "se = float(np.sqrt(1 / S * np.sum((ej - ej.mean()) ** 2)))", 0),
 ('F19', "the layer law's verdict at |z| ≤ 1", "law = 'consistent with the layer law' if abs(z) <= 2 else", "law = 'consistent with the layer law' if abs(z) <= 1 else", 0),
 ('F20', "the pooled s.e. 1/Σw", "se = float(1 / np.sqrt(w.sum()))", "se = float(1 / w.sum())", 0),
 ('F21', "confined ignores δ_cap", "ok = beyond == 0 and incap == 0", "ok = beyond == 0", 0),
 ('F22', "the window from the threshold itself", "groups = window_bins(2 * thr[k], ds[k] / 3)", "groups = window_bins(thr[k], ds[k] / 3)", 0),
 ('F23', "no burn-in in the exponent", "H = np.array([a['hist'][BURN_J:, k].sum(axis=0) for a in rr])", "H = np.array([a['hist'][:, k].sum(axis=0) for a in rr])", 0),
 ('F24', "the second half's mean over the whole run", "half = np.array([a['dec'][len(a['dec']) // 2:, 4:7].mean(axis=0) for a in rr])",
  "half = np.array([a['dec'][:, 4:7].mean(axis=0) for a in rr])", 0),
 ('F25', "the md's z with one decimal", "{c['v']['z']:+.2f}; {c['v']['law']}", "{c['v']['z']:+.1f}; {c['v']['law']}", 0),
]
