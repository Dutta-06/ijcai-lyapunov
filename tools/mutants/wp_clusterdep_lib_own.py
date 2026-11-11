"""Single-point mutants of experiments/astro_clusterdep.py (the coherence table, the merge of units the kernel cannot tell apart, the residual calibration and the weighted path), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp astro_clusterdep --tests tests/test_exp_astro_clusterdep.py --raw results/astro/clusterdep_raw.npz --md results/astro/clusterdep.md --out-env CLUSTERDEP_OUT \\
        --spec tools/mutants/wp_clusterdep_lib_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/astro_clusterdep.py. The kill criterion is the tests file alone (planted worlds; no real patch is read).

Final result (a first pass, then the survivors against the sharpened tests): 68 mutants, 66 killed, 2 survivors, both equivalent on every input:
  G09  `fi >= 0 and fj >= 0` for `fi > 0 and fj > 0` in the merge: a unit of exactly zero flux weighs nothing in the centroid, so the merged position is the other unit's either way.
  G18  `j != i` for `j > i` in the pair list: every pair is listed twice, the two keys are removed together after a merge, and the tie-break prefers the smaller first index, so the survivors are the same.
G10 (the absorbed unit stays alive) makes the merge loop forever; its pytest process was ended by hand after 5 minutes and counts as killed (the tool records it as a failure). A first pass left 22 survivors (the calibration
weights, the weighted path's clean-up, the anisotropic table, tie-breaks, the recomputation of a merged unit's pairs); each was killed by a test added for it.
"""
M = [
    # ---- weights and kernels
    ('W01', 'nominal_weights: the gain dropped', "var0 = expected / patch['gain'] + patch['dark']", "var0 = expected + patch['dark']"),
    ('W02', 'nominal_weights: the dark variance dropped', "var0 = expected / patch['gain'] + patch['dark']", "var0 = expected / patch['gain']"),
    ('W03', 'nominal_weights: the smoothing width', "expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))", "expected = np.maximum(gaussian_filter(dn, 2.5, mode='reflect'), np.percentile(dn, 1))"),
    ('W04', 'nominal_weights: the floor percentile', "np.percentile(dn, 1))\n    var0", "np.percentile(dn, 50))\n    var0"),
    ('W05', 'nominal_weights: sw is the variance', "return 1.0 / np.sqrt(var0), var0, expected", "return 1.0 / var0, var0, expected"),
    ('W06', 'make_kernel: the default fit is adequate-or-general', "fit = PG.fit_adequate(patch['psf'], rms_goal=np.inf)", "fit = PG.fit_adequate(patch['psf'])"),
    ('W07', 'make_kernel: the high-pass dropped', "return PG.make_psf(fit, ny, nx, hp)", "return PG.make_psf(fit, ny, nx, None)"),
    ('W08', 'robust_var: not squared', "return float((1.4826 * np.median(np.abs(v - np.median(v)))) ** 2)", "return float((1.4826 * np.median(np.abs(v - np.median(v)))))"),
    ('W09', 'robust_var: the constant of the MAD', "return float((1.4826 * np.median(np.abs(v - np.median(v)))) ** 2)", "return float((1.0 * np.median(np.abs(v - np.median(v)))) ** 2)"),
    ('W10', 'model_error_image: the signed model', "return gaussian_filter(np.abs(M), sigma, mode='reflect'), M", "return gaussian_filter(M, sigma, mode='reflect'), M"),
    ('W11', 'model_error_image: no smoothing', "return gaussian_filter(np.abs(M), sigma, mode='reflect'), M", "return np.abs(M), M"),
    # ---- the calibration
    ('K01', 'calibrate: the residual against the unfiltered frame', "y = P.highpass(patch['dn'] - np.median(patch['dn']))\n    r = y - M", "y = patch['dn'] - np.median(patch['dn'])\n    r = y - M"),
    ('K02', 'calibrate: the residual sign', "r = y - M", "r = y + M"),
    ('K03', 'calibrate: the border not removed', "sl = (slice(border, -border), slice(border, -border))", "sl = (slice(0, None), slice(0, None))"),
    ('K04', 'calibrate: bins by the variance', "order = np.argsort(q, kind='stable')", "order = np.argsort(v0, kind='stable')"),
    ('K05', 'calibrate: the squared model error not squared', "A1 = np.array([(q[b] ** 2).mean() for b in bins])", "A1 = np.array([(q[b]).mean() for b in bins])"),
    ('K06', 'calibrate: unweighted bins', "w = np.array([len(b) for b in bins], float)", "w = np.ones(len(bins))"),
    ('K07', 'calibrate: s0 below one allowed', "T = V - S0_MIN ** 2 * A0", "T = V - 0.0 * A0"),
    ('K08', 'calibrate: the unconstrained solution not offered', "cand.append(sol)\n", "pass\n"),
    ('K09', 'calibrate: negative coefficients accepted', "if s[0] < 0 or s[1] < 0:\n            continue", "if False:\n            continue"),
    ('K10', 'calibrate: the eta-only candidate', "cand.append(np.array([0.0, max(float(X2[:, 1] @ yv) / max(float(X2[:, 1] @ X2[:, 1]), 1e-300), 0.0)]))", "cand.append(np.array([0.0, 0.0]))"),
    ('K11', 'calibrate: the s0-only candidate', "cand.append(np.array([max(float(X2[:, 0] @ yv) / max(float(X2[:, 0] @ X2[:, 0]), 1e-300), 0.0), 0.0]))", "cand.append(np.array([0.0, 0.0]))"),
    ('K12', 'calibrate: the best candidate by the largest cost', "if best is None or c < best[0]:", "if best is None or c > best[0]:"),
    ('K13', 'calibrate: s0 without the unit floor', "s0=float(np.sqrt(S0_MIN ** 2 + d0))", "s0=float(np.sqrt(d0))"),
    ('K14', 'calibrate: eta is its square', "eta=float(np.sqrt(eta2))", "eta=float(eta2)"),
    ('K15', 'calibrated_weights: eta not applied', "sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'] + c['eta'] ** 2 * c['qd'] ** 2)", "sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'])"),
    ('K16', 'calibrated_weights: s0 not squared', "sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'] + c['eta'] ** 2 * c['qd'] ** 2)", "sw = 1.0 / np.sqrt(c['s0'] * c['var0'] + c['eta'] ** 2 * c['qd'] ** 2)"),
    ('K17', 'calibrated_weights: eta not squared', "sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'] + c['eta'] ** 2 * c['qd'] ** 2)", "sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'] + c['eta'] * c['qd'] ** 2)"),
    # ---- the weighted path
    ('P01', 'run_path_weighted: the kernel without the high-pass of the data', "y = P.highpass(patch['dn'] - np.median(patch['dn']))\n    reach", "y = patch['dn'] - np.median(patch['dn'])\n    reach"),
    ('P02', 'run_path_weighted: the noise scale from the nominal weights', "scale = IM.noise_scale(L, sw ** 2)\n    t0", "scale = IM.noise_scale(L, C_nominal(patch) ** 2)\n    t0"),
    ('P03', 'run_path_weighted: the stop time ignores kappa', "i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (k * scale))", "i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (2.0 * scale))"),
    ('P04', 'run_path_weighted: no clean-up', "keep = IM.unit_significance(P, sw, h['c'], h['w']) >= k\n        out[k]", "keep = np.ones(len(h['w']), bool)\n        out[k]"),
    ('P05', 'run_path_weighted: the clean-up with the wrong threshold', "keep = IM.unit_significance(P, sw, h['c'], h['w']) >= k\n        out[k]", "keep = IM.unit_significance(P, sw, h['c'], h['w']) >= 2.0 * k\n        out[k]"),
    ('P06', 'run_path_weighted: the exclusion radius ignored', "r_ex=rex, fit_radius", "r_ex=0.0, fit_radius"),
    ('P07', 'run_path_weighted: the BIC clean-up threshold', "kb = float(np.sqrt(3.0 * np.log(ny * nx)))", "kb = float(np.sqrt(1.0 * np.log(ny * nx)))"),
    ('P08', 'run_path_weighted: the entries of the path counted from the blocks', "return out, secs, len(R['taus']), dict(bic_index=int(i), n_blocks=len(R['hist']))", "return out, secs, len(R['hist']), dict(bic_index=int(i), n_blocks=len(R['hist']))"),
    # ---- the coherence table
    ('T01', 'coherence_table: the range of offsets', "off = np.arange(-COH_GRID, COH_GRID + 1e-9, COH_STEP)", "off = np.arange(-COH_GRID, COH_GRID - 1e-9, COH_STEP)"),
    ('T02', 'coherence_table: the normalisation', "tab[iy] = (ph0 @ ph) / (ph0 @ ph0)", "tab[iy] = (ph0 @ ph)"),
    ('T03', 'coherence_table: the rows are the x offsets', "ph = P.phi(np.c_[c0 + off, np.full(len(off), c0 + dy)])", "ph = P.phi(np.c_[np.full(len(off), c0 + dy), c0 + off])"),
    ('T04', 'coherence_table: the high-pass kept', "psf=patch['psf']), None, fit)", "psf=patch['psf']), 3.0, fit)"),
    ('T05', 'coherence_table: the frame too small for the offsets', "COH_GRID, COH_STEP, COH_FRAME = 4.0, 0.25, 41", "COH_GRID, COH_STEP, COH_FRAME = 4.0, 0.25, 7"),
    ('T06', 'coherence_at: the step of the grid', "fx = np.atleast_1d((np.asarray(dx, float) - off[0]) / (off[1] - off[0]))", "fx = np.atleast_1d((np.asarray(dx, float) - off[0]) / (off[1] - off[0]) * 0.5)"),
    ('T07', 'coherence_at: x and y exchanged', "out = map_coordinates(tab, [fy, fx], order=1, mode='constant', cval=0.0)", "out = map_coordinates(tab, [fx, fy], order=1, mode='constant', cval=0.0)"),
    ('T08', 'coherence_at: nearest-neighbour', "out = map_coordinates(tab, [fy, fx], order=1, mode='constant', cval=0.0)", "out = map_coordinates(tab, [fy, fx], order=0, mode='constant', cval=0.0)"),
    ('T09', 'coherence_at: the outside held at the edge', "out = map_coordinates(tab, [fy, fx], order=1, mode='constant', cval=0.0)", "out = map_coordinates(tab, [fy, fx], order=1, mode='nearest')"),
    ('T10', 'coherence_at: the scalar not returned as a float', "return out.reshape(shape) if shape else float(out[0])", "return out.reshape(shape) if shape else out"),
    # ---- the merge
    ('G01', 'merge_coherent: the threshold strict', "                if r >= rho:\n                    pairs[(i, int(j))] = float(r)", "                if r > rho:\n                    pairs[(i, int(j))] = float(r)"),
    ('G02', 'merge_coherent: the threshold strict for the neighbours of a merged unit', "            if rr >= rho:\n                pairs[(min(i, int(k)), max(i, int(k)))] = float(rr)", "            if rr > rho:\n                pairs[(min(i, int(k)), max(i, int(k)))] = float(rr)"),
    ('G03', 'merge_coherent: the least coherent pair first', "(i, j), r = max(pairs.items(), key=lambda kv: (kv[1], -kv[0][0], -kv[0][1]))", "(i, j), r = min(pairs.items(), key=lambda kv: (kv[1], -kv[0][0], -kv[0][1]))"),
    ('G04', 'merge_coherent: ties to the larger indices', "(i, j), r = max(pairs.items(), key=lambda kv: (kv[1], -kv[0][0], -kv[0][1]))", "(i, j), r = max(pairs.items(), key=lambda kv: (kv[1], kv[0][0], kv[0][1]))"),
    ('G05', 'merge_coherent: the unweighted centroid', "det[i, 0] = (det[i, 0] * fi + det[j, 0] * fj) / (fi + fj)", "det[i, 0] = 0.5 * (det[i, 0] + det[j, 0])"),
    ('G06', 'merge_coherent: y not updated', "det[i, 1] = (det[i, 1] * fi + det[j, 1] * fj) / (fi + fj)", "det[i, 1] = det[i, 1]"),
    ('G07', 'merge_coherent: the flux not summed', "det[i, 2] = fi + fj", "det[i, 2] = max(fi, fj)"),
    ('G08', 'merge_coherent: the non-positive unit moves the position', "elif fj > fi:\n            det[i, 0], det[i, 1] = det[j, 0], det[j, 1]", "elif fj < fi:\n            det[i, 0], det[i, 1] = det[j, 0], det[j, 1]"),
    ('G09', 'merge_coherent: the positive-flux test', "if fi > 0 and fj > 0:", "if fi >= 0 and fj >= 0:"),
    ('G10', 'merge_coherent: the absorbed unit stays alive', "alive[j] = False\n        merges += 1", "merges += 1"),
    ('G11', 'merge_coherent: the merge count', "alive[j] = False\n        merges += 1", "alive[j] = False\n        merges += 2"),
    ('G12', 'merge_coherent: the stale pairs of the merged unit kept', "pairs = {k: v for k, v in pairs.items() if i not in k and j not in k}", "pairs = {k: v for k, v in pairs.items() if j not in k}"),
    ('G13', 'merge_coherent: the pairs of the absorbed unit kept', "pairs = {k: v for k, v in pairs.items() if i not in k and j not in k}", "pairs = {k: v for k, v in pairs.items() if i not in k}"),
    ('G14', 'merge_coherent: the new pairs of the merged unit not searched', "        for k, rr in zip(near, pair_r(i, near)):", "        for k, rr in zip(near[:0], pair_r(i, near[:0])):"),
    ('G15', 'merge_coherent: the merged unit compared with itself', "others = others[others != i]", "others = others"),
    ('G16', 'merge_coherent: the search window around the merged unit shrunk', "near = others[(np.abs(det[others, 0] - det[i, 0]) <= reach) & (np.abs(det[others, 1] - det[i, 1]) <= reach)]", "near = others[(np.abs(det[others, 0] - det[i, 0]) <= 1.0) & (np.abs(det[others, 1] - det[i, 1]) <= 1.0)]"),
    ('G17', 'merge_coherent: the reach', "reach = COH_GRID - 1e-6", "reach = COH_GRID - 2.0"),
    ('G18', 'merge_coherent: pairs counted twice (j > i dropped)', "js = [j for j in js if j > i]", "js = [j for j in js if j != i]"),
    ('G19', 'merge_coherent: the input modified', "det = np.asarray(det, float).reshape(-1, 3).copy()", "det = np.asarray(det, float).reshape(-1, 3)"),
    ('G20', 'merge_coherent: the survivors returned with the absorbed', "return det[alive], merges", "return det, merges"),
    ('G21', 'merge_coherent: the query radius', "r=np.sqrt(2.0) * reach)", "r=1.0)"),
    ('G22', 'merge_coherent: the offset sign in pair_r', "dx, dy = det[i, 0] - det[js, 0], det[i, 1] - det[js, 1]", "dx, dy = det[i, 0] - det[js, 0], det[i, 1] + det[js, 1]"),
]
