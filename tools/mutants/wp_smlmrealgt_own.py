# (id, what, old, new, k): single-point mutants of exp_smlm_real_gt.py (the pre-registered real-data test on nuclear pores): the pore finder, circle fit and 8-fold phase,
# the pore model (convention shift, assignment, areas), the radial-profile fit's J_s / E_s, the coordinate ascent, Picasso's post-filter semantics, the ISS configuration
# grid and the operating-point choice, the equal-yield thinning, the intervals, the C8 harmonic, the reading of the outcomes (margin, decision rule) and the frame offsets.
# Run: tools/mutants.py --exp exp_smlm_real_gt --tests tests/test_exp_smlm_real_gt.py --raw results/smlm/real_gt_raw.npz --md results/smlm/real_gt.md --out-env REALGT_OUT
#      --spec tools/mutants/wp_smlmrealgt_own.py --work <scratch dir> (one slot; the tests take about 30 s, the --from-raw regeneration of the md the rest).
# Result: 28 mutants and a control. Against the first version of tests/test_exp_smlm_real_gt.py (26 tests) 14 survived (R01-R03, R05, R15-R21, R23, R27, R28; the tests had not looked at the
# pore finder's radius rule and merge distance, the phase's range and window, the inclusive width bound, the equal-yield thinning's keep probability and number, the interval's level, the
# block size, the table's frame base, the margin, the rule's upper limit, the order forced into the path and the kappa of every block); the tests added then
# (test_find_pores_rejects..., test_pore_phase_is_reported..., test_post_width_keeps..., test_equal_yield_thinning_thins..., test_percentile_ci..., test_prep_row_assigns...,
# test_table_loc_gives..., test_verdict_uses_the_three_point_margin..., test_the_iss_calls_the_guarded_path...) kill all of them. No survivor; the control passes.
# The criterion is the tests file alone (pytest -x, one private copy of src/, experiments/ and tests/, run by a driver equivalent to tools/mutants.py without its --from-raw step: the
# regeneration of the md takes 12 minutes per mutant, so the md column was NOT run).
M = [
 ('R00', 'control', "def find_pores(x, y, region, sp):", "def find_pores(x, y, region, sp):", 0),
 ('R01', "accepted radius within +-50 % of the median", "ok = np.abs(fits[:, 2] - med) <= RADIUS_TOL * med", "ok = np.abs(fits[:, 2] - med) <= 2 * RADIUS_TOL * med", 0),
 ('R02', "merge distance 3 R_start", "dmin = 1.5 * sp['r_start_px'] * pix", "dmin = 3.0 * sp['r_start_px'] * pix", 0),
 ('R03', "the 8-fold phase modulo 90 degrees", "np.mod(np.angle(Z) / 8.0, np.pi / 4)", "np.mod(np.angle(Z) / 8.0, np.pi / 2)", 0),
 ('R04', "the harmonic of order 4", "Z = np.exp(8j * th).sum() if n else 0j", "Z = np.exp(4j * th).sum() if n else 0j", 0),
 ('R05', "the ring window of the phases without the half width", "m = np.abs(r - R) <= sp['ring_half']\n        th", "m = np.abs(r - R) <= 2 * sp['ring_half']\n        th", 0),
 ('R06', "the convention rounded up to a whole pixel", "return np.round(off / half) * half", "return np.round(off / (2 * half)) * (2 * half)", 0),
 ('R07', "the cell area with the pitch, not its square", "self.area = np.bincount(b, minlength=self.nb) * pitch * pitch", "self.area = np.bincount(b, minlength=self.nb) * pitch", 0),
 ('R08', "the per-pore area with the pitch only", "np.add.at(ap, (pos[j[ok]], b), pitch * pitch)", "np.add.at(ap, (pos[j[ok]], b), pitch)", 0),
 ('R09', "every pore is scored", "keep[ok] &= self.scored[j[ok]]", "keep[ok] &= True", 0),
 ('R10', "the missed yield without the clip at 0", "miss = max(0.0, n_exp - fit['N'])", "miss = n_exp - fit['N']", 0),
 ('R11', "E_s with the width in nm / 2", "(alpha * fit['w']) ** 2)", "(alpha * fit['w'] / 2) ** 2)", 0),
 ('R12', "the ascent moves on ties", "if x > best + 1e-12:", "if x >= best:", 0),
 ('R13', "the net gradient filter inclusive", "m &= a[:, C.NG_COL] > f['min_net_gradient']", "m &= a[:, C.NG_COL] >= f['min_net_gradient']", 0),
 ('R14', "the r_ex variant at 2 s", "cfg['r_ex'] = 2.5 * s", "cfg['r_ex'] = 2.0 * s", 0),
 ('R15', "the width post-filter exclusive", "d[d[:, 3] <= smax + 1e-9]", "d[d[:, 3] < smax - 1e-9]", 0),
 ('R16', "the equal-yield keep probability inverted", "q = min(small['N'] / max(big['N'], 1e-9), 1.0)", "q = min(big['N'] / max(small['N'], 1e-9), 1.0)", 0),
 ('R17', "one thinning only", "    for _ in range(N_THIN):\n        Ht = rng.binomial", "    for _ in range(1):\n        Ht = rng.binomial", 0),
 ('R18', "the upper percentile at the 90 %", "float(np.percentile(v, 100 * (1 + level) / 2))", "float(np.percentile(v, 100 * (1 + 0.9) / 2))", 0),
 ('R19', "blocks of 1 frame", "blk = np.minimum(loc.frame // BLOCK, nblocks - 1)", "blk = np.minimum(loc.frame // 1, nblocks - 1)", 0),
 ('R20', "the table's frames 1-based in the set", "return Loc(tab['frame'] - a, np.c_[tab['x'], tab['y']] + conv)", "return Loc(tab['frame'] - a + 1, np.c_[tab['x'], tab['y']] + conv)", 0),
 ('R21', "the equivalence margin 5", "MARGIN_E = 3.0 ", "MARGIN_E = 5.0 ", 0),
 ('R22', "ISS better needs only the E interval", "better = dE[0] > 0 and dw[1] < 0", "better = dE[0] > 0", 0),
 ('R23', "Picasso better on the upper limit above 0", "worse = dE[1] < 0", "worse = dE[1] < 1", 0),
 ('R24', "the pore sensitivity ignored", "better_pores = pr['dE_ci_pores'][0] > 0 and pr['dw_ci_pores'][1] < 0", "better_pores = True", 0),
 ('R25', "the odd blocks as even", "return {k: v[inside & (blk == parity)] for k, v in tab.items()}", "return {k: v[inside & (blk != parity)] for k, v in tab.items()}", 0),
 ('R26', "the disc margin ignored", "d = R + margin\n    return", "d = R\n    return", 0),
 ('R27', "the greedy order not forced", "fit_radius=cfg['fit_radius'], order=order, snr_min", "fit_radius=cfg['fit_radius'], order='dual', snr_min", 0),
 ('R28', "kappa blocks all taken at 4.5", "i = IM.stop_snr(R['hist'], R['taus'], R['snr'], kap)[0]", "i = IM.stop_snr(R['hist'], R['taus'], R['snr'], 4.5)[0]", 0),
]
