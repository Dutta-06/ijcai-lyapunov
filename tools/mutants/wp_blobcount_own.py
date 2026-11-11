"""Single-point mutants of experiments/exp_blobcount.py (nucleus counting spike), for tools/mutants.py:

    PYTHONPATH=<dir with scikit-image> venv/bin/python tools/mutants.py --exp exp_blobcount --tests tests/test_exp_blobcount.py --raw results/ml/blobcount_raw.npz \\
        --md results/ml/blobcount.md --extra results/ml/blobcount.md --extra results/ml/blobcount_dev_raw.npz --extra results/ml/blobcount_dev.md --out-env BLOBCOUNT_OUT \\
        --spec tools/mutants/wp_blobcount_own.py --work <scratch dir>

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence (0-based). The kill criterion is the tests file alone. The run parts (the notebook driver, `stage_*`)
need the dataset and are not mutated here; the data rules, front end, scoring, statistics, decision rule and report are.

Result (tests/test_exp_blobcount.py, 47 tests, OMP_NUM_THREADS=1; 58 mutants and the control B00, which passes): 57 killed, 1 survivor. First version of the tests (34 tests): 13 survived (B03, B04, B10, B15, B18-B23, B26-B28, B38) and B30's anchor was stale;
the tests added then (third-channel-only colour test, colour-branch median limit, the y range of the patch membership, an isolated multi-pixel mask for the touching test, a tighter chi-square-median check of the noise model, the floor and the clamped slope of the
line fit, the variance-map formula and floor, a scale estimate that must not follow the strongest blob only, fake-path tests of the BIC calibration (floor at 1, scaling above it, off) and of the clean-up threshold sqrt(lam ln N)) kill all but one.
Equivalent mutant: B38 (the columns of a matrix of independent uniform draws reversed: the resampling distribution is the same; the interval changes only through the generator's stream, and with 10 000 draws the percentile moves by
less than the md's rounding, so the md is byte-identical too).
"""
M = [
 ('B00', 'control', 'def is_fluorescence(rgb):', 'def is_fluorescence(rgb):', 0),
 ('B01', 'test images shifted by one', 'return order[:N_DEV], order[N_DEV:N_DEV + N_TEST]', 'return order[:N_DEV], order[N_DEV + 1:N_DEV + N_TEST + 1]', 0),
 ('B02', 'split not seeded by the plan seed', 'default_rng(SEED_SPLIT).permutation(len(ids))', 'default_rng(SEED_SPLIT + 1).permutation(len(ids))', 0),
 ('B03', 'colour channels need not be identical', 'np.all(c[..., 0] == c[..., 1]) and np.all(c[..., 1] == c[..., 2])', 'np.all(c[..., 0] == c[..., 1])', 0),
 ('B04', 'median limit inclusive (colour branch)', 'np.median(c[..., 0]) < MED_MAX', 'np.median(c[..., 0]) <= MED_MAX', 0),
 ('B05', 'median limit inclusive (grey branch)', 'return bool(np.median(rgb) < MED_MAX)', 'return bool(np.median(rgb) <= MED_MAX)', 0),
 ('B06', 'debris limit exclusive', 'if a < MIN_AREA:', 'if a <= MIN_AREA:', 0),
 ('B07', 'centroid from rows and columns swapped', 'xy.append((c.mean(), r.mean()))', 'xy.append((r.mean(), c.mean()))', 0),
 ('B08', 'patch membership lower edge exclusive', 'return (x >= -0.5) & (x < PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)', 'return (x > -0.5) & (x < PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)', 0),
 ('B09', 'patch membership upper edge inclusive', 'return (x >= -0.5) & (x < PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)', 'return (x >= -0.5) & (x <= PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)', 0),
 ('B10', 'y edge of the membership ignored', 'return (x >= -0.5) & (x < PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)', 'return (x >= -0.5) & (x < PATCH - 0.5) & (y >= -0.5)', 0),
 ('B11', 'one base tile too many', 'base = [int(i) for i in perm[:PER_IMAGE]]', 'base = [int(i) for i in perm[:PER_IMAGE + 1]]', 0),
 ('B12', 'boost counts exclusive', 'if counts[i] >= CROWD][:BOOST]', 'if counts[i] > CROWD][:BOOST]', 0),
 ('B13', 'boost may repeat the base tiles', 'for i in perm[PER_IMAGE:] if counts[i] >= CROWD]', 'for i in perm if counts[i] >= CROWD]', 0),
 ('B14', 'boost not limited', 'if counts[i] >= CROWD][:BOOST]', 'if counts[i] >= CROWD][:BOOST + 1]', 0),
 ('B15', 'touching without the different-label test', 'm = (a > 0) & (b > 0) & (a != b)', 'm = (a > 0) & (b > 0)', 0),
 ('B16', 'boost patches inside the base sample', "base = np.flatnonzero(sel & (bank['boost'] == 0))", "base = np.flatnonzero(sel)", 0),
 ('B17', 'crowded stratum exclusive', "crowd = np.flatnonzero(sel & (bank['n_gt'] >= CROWD))", "crowd = np.flatnonzero(sel & (bank['n_gt'] > CROWD))", 0),
 ('B18', 'median of chi-square wrong', 'vs.append(float(np.median(d2[s])) / 0.4549)', 'vs.append(float(np.median(d2[s])) / 0.5)', 0),
 ('B19', 'noise floor of the line fit', 'return float(a), float(max(b, 0.25))', 'return float(a), float(max(b, 0.5))', 0),
 ('B20', 'negative slope allowed', 'if a < 0:\n        a, b = 0.0', 'if a < -1e9:\n        a, b = 0.0', 0),
 ('B21', 'variance floor', 'var = np.maximum(a * ndi.gaussian_filter(g, 1.5) + b, 0.25)', 'var = np.maximum(a * ndi.gaussian_filter(g, 1.5) + b, 0.5)', 0),
 ('B22', 'variance map unsmoothed', 'var = np.maximum(a * ndi.gaussian_filter(g, 1.5) + b, 0.25)', 'var = np.maximum(a * g + b, 0.25)', 0),
 ('B23', 'scale from the median of too few maxima', 'def scale_estimate(p, sigmas=None, top=12, default=5.0):', 'def scale_estimate(p, sigmas=None, top=1, default=5.0):', 0),
 ('B24', 'scale default', 'def scale_estimate(p, sigmas=None, top=12, default=5.0):', 'def scale_estimate(p, sigmas=None, top=12, default=4.0):', 0),
 ('B25', 'LoG not scale-normalised', 'R = np.stack([-(s * s) * ndi.gaussian_laplace(np.asarray(p, float), s, mode=\'reflect\') for s in sigmas])', 'R = np.stack([-s * ndi.gaussian_laplace(np.asarray(p, float), s, mode=\'reflect\') for s in sigmas])', 0),
 ('B26', 'calibration without the floor at 1', "max(hist[-1]['chi2'], 1.0)", "max(hist[-1]['chi2'], 0.0)", 0),
 ('B27', 'calibration switched off', "if stop.get('cal', 1) else 1.0)", "if False else 1.0)", 0),
 ('B28', 'clean-up threshold without the log', "np.sqrt(lam * np.log(N))", "np.sqrt(lam)", 0),
 ('B29', 'penalty without lam', 'lam = stop[\'lam\'] * (', 'lam = 1.0 * (', 0),
 ('B30', 'orders swapped', "return 'dual' if method in ('ours', 'ours_nc') else 'greedy'", "return 'greedy' if method in ('ours', 'ours_nc') else 'dual'", 0),
 ('B31', 'match radius', 'R_MATCH = 6.0 ', 'R_MATCH = 7.0 ', 0),
 ('B32', 'match at the radius inclusive -> exclusive', 'cost = np.where(d <= r, d, 1e6)', 'cost = np.where(d < r, d, 1e6)', 0),
 ('B33', 'matched pairs counted at the radius exclusive', 'return int(np.sum(cost[i, j] <= r))', 'return int(np.sum(cost[i, j] < r))', 0),
 ('B34', 'greedy nearest instead of the optimal assignment', 'i, j = linear_sum_assignment(cost)', 'i, j = np.arange(min(cost.shape)), cost[:min(cost.shape)].argmin(1)[:min(cost.shape)]', 0),
 ('B35', 'detections outside the patch counted', 'd = inside(np.asarray(dets, float)[:, :2]) if len(dets) else np.zeros((0, 2))', 'd = np.asarray(dets, float)[:, :2] if len(dets) else np.zeros((0, 2))', 0),
 ('B36', 'lower bootstrap percentile', 'float(np.percentile(est, 2.5))', 'float(np.percentile(est, 5.0))', 0),
 ('B37', 'upper bootstrap percentile', 'float(np.percentile(est, 97.5))', 'float(np.percentile(est, 95.0))', 0),
 ('B38', 'bootstrap resamples patches not images', 'k = rng.integers(0, len(ids), size=(n_boot, len(ids)))', 'k = rng.integers(0, len(ids), size=(n_boot, len(ids)))[:, ::-1]', 0),
 ('B39', 'mean of the image means', 'est = s[k].sum(1) / n[k].sum(1)', 'est = (s[k] / n[k]).mean(1)', 0),
 ('B40', 'absolute value dropped', 'd = np.abs(ea[idx]) - np.abs(eb[idx])', 'd = ea[idx] - eb[idx]', 0),
 ('B41', 'C1 margin', 'M_REL = 0.05 ', 'M_REL = 0.10 ', 0),
 ('B42', 'C1 interval may touch 0', "table['diff_c'][1] > 0.0", "table['diff_c'][1] >= 0.0", 0),
 ('B43', 'C1 uses the upper end', "table['diff_c'][1] > 0.0", "table['diff_c'][2] > 0.0", 0),
 ('B44', 'C2 strict', "c2 = table['mae_o'] <= table['mae_g'] + M_GREEDY", "c2 = table['mae_o'] < table['mae_g'] + M_GREEDY", 0),
 ('B45', 'C2 margin', 'M_GREEDY = 0.5 ', 'M_GREEDY = 1.0 ', 0),
 ('B46', 'C3 strict', "c3 = table['diff_all'][2] <= M_NONINF * table['mae_all_c']", "c3 = table['diff_all'][2] < M_NONINF * table['mae_all_c']", 0),
 ('B47', 'C3 margin', 'M_NONINF = 0.10 ', 'M_NONINF = 0.20 ', 0),
 ('B48', 'C3 uses the lower end', "c3 = table['diff_all'][2] <=", "c3 = table['diff_all'][1] <=", 0),
 ('B49', 'GO needs only C1 and C3', 'GO=bool(c1 and c2 and c3)', 'GO=bool(c1 and c3)', 0),
 ('B50', 'dual contributes on the point estimate', "dual_contributes=bool(table['diff_og'][2] < 0.0)", "dual_contributes=bool(table['diff_og'][0] < 0.0)", 0),
 ('B51', 'objective weights', 'return 0.5 * (np.abs(err[:, base]).mean(1) + np.abs(err[:, crowd]).mean(1))', 'return np.abs(err[:, base]).mean(1)', 0),
 ('B52', 'F1 precision and recall swapped roles (recall over detections)', 'r = tp.sum() / max(ngt.sum(), 1)', 'r = tp.sum() / max(ndet.sum(), 1)', 0),
 ('B53', 'relative error over all patches', 'pos = ngt[idx] >= 1', 'pos = ngt[idx] >= 0', 0),
 ('B54', 'best classical chosen by the highest MAE', 'best[st] = min(cl, key=lambda m: pooled(err[m], idx[st]))', 'best[st] = max(cl, key=lambda m: pooled(err[m], idx[st]))', 0),
 ('B55', 'best classical includes greedy', 'CLASSICAL = (\'log\', \'dog\', \'plm\', \'ws\')', 'CLASSICAL = (\'log\', \'dog\', \'plm\', \'ws\', \'greedy\')', 0),
 ('B56', 'verdict taken against the all-patch classical in the crowded stratum', "ea, eo = err[best['crowded']], err['ours']", "ea, eo = err[best['all']], err['ours']", 0),
 ('B57', 'freeze picks the maximum', 'k = int(np.argmin(obj))\n        res[m]', 'k = int(np.argmax(obj))\n        res[m]', 0),
 ('B58', 'bias sign', "e = (nd[idx] - ngt[idx]).astype(float)", "e = (ngt[idx] - nd[idx]).astype(float)", 0),
]
