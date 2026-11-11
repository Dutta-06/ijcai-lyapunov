"""Single-point mutants of experiments/exp_gam_spike.py (the pure parts: the dataset rules, the metrics and counts, the dev -> frozen pick, the
criteria and the verdict, the report), run with

    venv/bin/python tools/mutants.py --exp exp_gam_spike --tests tests/test_gam_spike.py --raw results/ml/gam_raw.npz --md results/ml/gam.md \\
        --out-env GAM_OUT --spec tools/mutants/gam_spike_own.py --work <scratch>/mut --slot a \\
        --extra results/ml/gam_probe.json --extra results/ml/gam_frozen.json --extra results/ml/gam_dev_all.jsonl

(id, what, old, new, k): `old` replaced by `new` at its k-th occurrence. G00 is the control (no change). The torch / interpret / xgboost parts
(training loops, baselines) are not mutated: they cannot run locally and their numbers are the experiment's output, not a tested rule.

Run (80 entries, the control and 79 mutants, 7 min): 76 of the 79 killed by the tests (each also leaves `--from-raw` working; the control passes). Three survivors, all equivalent:
G10 (`n <= cap` -> `n < cap`: at n == cap the sorted first `cap` of a permutation of n = cap rows is `arange(cap)`, the same as the early return);
G22 (the AUC label threshold 0.5 -> 0.0: the labels are 0 / 1, so `y > 0` and `y > 0.5` coincide);
G74 (the post hoc "not worse than the control" count with `<=` -> `<`: the comparison already carries the 1e-9 slack, so the two differ only when the gap equals the tie band plus 1e-9 exactly)."""
M = [
    ('G00', 'control', '', '', 0),
    # data rules
    ('G01', 'rule: n >= N_MIN -> n > N_MIN', "if r['n'] >= N_MIN and", "if r['n'] > N_MIN and", 0),
    ('G02', 'rule: d_enc <= D_MAX -> <', "3 <= r['d_enc'] <= D_MAX", "3 <= r['d_enc'] < D_MAX", 0),
    ('G03', 'rule: d_enc >= 3 -> >= 2', "3 <= r['d_enc'] <= D_MAX", "2 <= r['d_enc'] <= D_MAX", 0),
    ('G04', 'rule: the kind must match', "r is None or r['kind'] != kind or", "r is None or", 0),
    ('G05', 'rule: take the first n (>= instead of >)', "len(got[kind]) >= n_names[kind]", "len(got[kind]) > n_names[kind]", 0),
    ('G06', 'roles: dev if rank < n_dev -> <=', "'dev' if rank < n_dev[kind] else 'test'", "'dev' if rank <= n_dev[kind] else 'test'", 0),
    ('G07', 'roles: sorted names -> given order', "ns = sorted(n for n in names if kinds[n] == kind)", "ns = [n for n in names if kinds[n] == kind]", 0),
    ('G08', 'roles: one permutation per kind from the same stream -> reseeded', "    rng = np.random.RandomState(seed)\n    roles = {}", "    roles = {}", 0),
    ('G09', 'subsample: unsorted', "return np.sort(np.random.RandomState(seed).permutation(n)[:cap])", "return np.random.RandomState(seed).permutation(n)[:cap]", 0),
    ('G10', 'subsample: n <= cap -> n < cap', "if n <= cap:", "if n < cap:", 0),
    ('G11', 'split: seed offset 1000 -> 0', "RandomState(1000 + seed).permutation(n)", "RandomState(seed).permutation(n)", 0),
    ('G12', 'split: test = last part', "return p[:a], p[a:b], p[b:]", "return p[:a], p[b:], p[a:b]", 0),
    ('G13', 'encode: two-level -> the earlier level', "use = levels[1:] if len(levels) == 2 else levels", "use = levels[:1] if len(levels) == 2 else levels", 0),
    ('G14', 'encode: max_levels inclusive -> exclusive', "len(levels) > max_levels", "len(levels) >= max_levels", 0),
    ('G15', 'encode: missing numeric -> mean', "np.where(np.isnan(a), np.nanmedian(a), a)", "np.where(np.isnan(a), np.nanmean(a), a)", 0),
    ('G16', 'encode: keep constant numeric columns', "            if a.max() == a.min():\n                continue\n", "", 0),
    ('G17', 'quantizer: binary threshold >', "(X[:, j] > u[0]).astype(float)", "(X[:, j] >= u[0]).astype(float)", 0),
    ('G18', 'quantizer: binary cut at <= 2 -> <= 3', "if len(u) <= 2:", "if len(u) <= 3:", 0),
    ('G19', 'quantizer: ties not averaged', "pr /= np.bincount(inv)", "pr /= 1", 0),
    ('G20', 'quantizer: 513 -> 257 nodes in the CDF (fit only)', "q = np.quantile(X[:, j], np.linspace(0, 1, 513))", "q = np.quantile(X[:, j], np.linspace(0, 1, 5))", 0),
    # metrics and counts
    ('G21', 'auc: sign of the rank sum', "(r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)", "(r[y].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)", 0),
    ('G22', 'auc: threshold on y', "y = np.asarray(y) > 0.5", "y = np.asarray(y) > 0.0", 0),
    ('G23', 'rmse: no square root', "return float(np.sqrt(np.mean((np.asarray(yz) - np.asarray(pz)) ** 2)))", "return float(np.mean((np.asarray(yz) - np.asarray(pz)) ** 2))", 0),
    ('G24', 'score: 1 - AUC -> AUC', "return metric if kind == 'reg' else 1.0 - metric", "return metric if kind == 'reg' else metric", 0),
    ('G25', 'active: |w| > theta -> >=', "return np.abs(np.asarray(w)) > theta", "return np.abs(np.asarray(w)) >= theta", 0),
    ('G26', 'active: THETA 0.01 -> 0.02', "THETA = 0.01 ", "THETA = 0.02 ", 0),
    ('G27', 'params: 2 per bump -> 3', "return int(2 * np.sum(mask))", "return int(3 * np.sum(mask))", 0),
    ('G28', 'segments: tolerance', "np.abs(np.diff(v)) > rtol * scale", "np.abs(np.diff(v)) > 0", 0),
    ('G29', 'ebm params: 2s - 1 -> 2s', "2 * segments(np.asarray(t)[1:-1]) - 1", "2 * segments(np.asarray(t)[1:-1])", 0),
    ('G30', 'ebm params: keep the missing / unknown bins', "segments(np.asarray(t)[1:-1])", "segments(np.asarray(t))", 0),
    ('G31', 'nam params: h^2 + 4h + 1 -> h^2 + 3h + 1', "return int(d * (h * h + 4 * h + 1) + 1)", "return int(d * (h * h + 3 * h + 1) + 1)", 0),
    ('G32', 'nam1 params: 3h + 1 -> 3h', "return int(d * (3 * h + 1) + 1)", "return int(d * (3 * h) + 1)", 0),
    ('G33', 'matched width: >= target -> >', "if nam1_params(d, h) >= target:", "if nam1_params(d, h) > target:", 0),
    ('G34', 'extrema: reversal tolerance', "tol = rel_tol * rng", "tol = 0.0", 0),
    ('G35', 'extrema: count only maxima', "                ext, direction, ref = ext + 1, -1, x", "                ext, direction, ref = ext, -1, x", 0),
    ('G36', 'stability: one flat shape counts 0 -> skipped', "                elif fa != fb:\n                    cs.append(0.0)\n", "", 0),
    ('G37', 'stability: unweighted', "wgt = float(sd.mean())", "wgt = 1.0", 0),
    ('G38', 'stability: all pairs -> neighbours', "for b in range(a + 1, S):", "for b in range(a + 1, min(S, a + 2)):", 0),
    # dev -> frozen, tests and criteria
    ('G39', 'pick: ratio to the best -> raw mean', "rs.append(mean[(c, d)] / max(mn, 1e-12))", "rs.append(mean[(c, d)])", 0),
    ('G40', 'pick: ties go to the highest id', "best[m] = min(cfgs, key=lambda c: (ratio[c], c))", "best[m] = min(cfgs, key=lambda c: (ratio[c], -c))", 0),
    ('G41', 'pick: mean over seeds -> first seed', "np.mean([r['score'] for r in rm if r['cfg_id'] == c and r['dataset'] == d])", "[r['score'] for r in rm if r['cfg_id'] == c and r['dataset'] == d][0]", 0),
    ('G42', 'paired p: one-sided -> two-sided', "return float(tdist.sf(d.mean() / (sd / np.sqrt(n)), n - 1))", "return float(2 * tdist.sf(abs(d.mean() / (sd / np.sqrt(n))), n - 1))", 0),
    ('G43', 'paired p: df n - 1 -> n', "np.sqrt(n)), n - 1))", "np.sqrt(n)), n))", 0),
    ('G44', 'paired p: ddof 1 -> 0', "sd = d.std(ddof=1)", "sd = d.std(ddof=0)", 0),
    ('G45', 'gap: relative -> absolute for regression', "return (s_ours - s_other) / s_other if kind == 'reg' else s_ours - s_other", "return s_ours - s_other", 0),
    ('G46', 'criteria A: margin 2 % -> 5 %', "MARGIN = {'reg': 0.02, 'clf': 0.005}", "MARGIN = {'reg': 0.05, 'clf': 0.005}", 0),
    ('G47', 'criteria A: AUC margin', "MARGIN = {'reg': 0.02, 'clf': 0.005}", "MARGIN = {'reg': 0.02, 'clf': 0.01}", 0),
    ('G48', 'criteria A: parameters 0.7 -> 0.8', "PARAM_RATIO = 0.7 ", "PARAM_RATIO = 0.8 ", 0),
    ('G49', 'criteria A: best of NAM / EBM -> NAM only', "best = min(('nam', 'ebm'), key=lambda m: mean[m])", "best = 'nam'", 0),
    ('G50', 'criteria A: parameters condition dropped', "A = bool(ni and (fewer or not a_params))", "A = bool(ni)", 0),
    ('G51', 'criteria B: p 0.05 -> 0.10', "P_ALPHA = 0.05", "P_ALPHA = 0.10", 0),
    ('G52', 'criteria B: the sign of the paired difference', "diff = np.asarray(sc['nam_matched']) - ours", "diff = ours - np.asarray(sc['nam_matched'])", 0),
    ('G53', 'criteria B: matched needs >= ours parameters', "and par['nam_matched'] >= par['ours_gd'])", ")", 0),
    ('G54', 'criteria C: tie band', "TIE = {'reg': 0.005, 'clf': 0.002}", "TIE = {'reg': 0.01, 'clf': 0.002}", 0),
    ('G55', 'criteria C: AUC tie band', "TIE = {'reg': 0.005, 'clf': 0.002}", "TIE = {'reg': 0.005, 'clf': 0.004}", 0),
    ('G56', 'criteria C: better control by L1 only', "cb = min(ctrls, key=lambda m: mean[m])", "cb = 'ctrl_l1'", 0),
    ('G57', 'criteria C: parsimony against the L2 control', "par_win = par['ours_gd'] <= PARAM_RATIO * par['ctrl_l1'] + EPS", "par_win = par['ours_gd'] <= PARAM_RATIO * par['ctrl_l2'] + EPS", 0),
    ('G58', 'criteria C: a tie with fewer parameters is not a win', "(acc == 'tie' and par_win)", "(acc == 'tie' and False)", 0),
    ('G59', 'criteria C: loss band', "'loss' if g > TIE[kind] * tie_scale + EPS", "'loss' if g > 2 * TIE[kind] * tie_scale + EPS", 0),
    ('G60', 'verdict: A needs 8 -> 7', "N_NEED_A, N_NEED_B, N_NEED_C = 8, 8, 7", "N_NEED_A, N_NEED_B, N_NEED_C = 7, 8, 7", 0),
    ('G61', 'verdict: B needs 8 -> 7', "N_NEED_A, N_NEED_B, N_NEED_C = 8, 8, 7", "N_NEED_A, N_NEED_B, N_NEED_C = 8, 7, 7", 0),
    ('G62', 'verdict: C needs 7 -> 6', "N_NEED_A, N_NEED_B, N_NEED_C = 8, 8, 7", "N_NEED_A, N_NEED_B, N_NEED_C = 8, 8, 6", 0),
    ('G63', 'verdict: (A or B) and C -> A and B and C', "go = (nA >= N_NEED_A or nB >= N_NEED_B) and nC >= N_NEED_C", "go = (nA >= N_NEED_A and nB >= N_NEED_B) and nC >= N_NEED_C", 0),
    ('G64', 'verdict: C or the rest', "go = (nA >= N_NEED_A or nB >= N_NEED_B) and nC >= N_NEED_C", "go = (nA >= N_NEED_A or nB >= N_NEED_B) or nC >= N_NEED_C", 0),
    ('G65', 'verdict: C counts ties', "nC = sum(c['C'] == 'win' for c in crit.values())", "nC = sum(c['C'] != 'loss' for c in crit.values())", 0),
    ('G66', 'boundary slack off (EPS)', "EPS = 1e-9 ", "EPS = 0.0 ", 0),
    # report
    ('G67', 'report: Errors count', "f'Errors: {sum(\"error\" in r for r in recs)}.'", "f'Errors: 0.'", 0),
    ('G68', 'report: criteria computed only on complete datasets -> all', "if all(sco(ds, m) is not None and len(sco(ds, m)) == len(seeds) for m in need):", "if True:", 0),
    ('G69', 'report: verdict label', "{\"GO\" if v[\"go\"] else \"KILL\"}", "{\"GO\" if not v[\"go\"] else \"KILL\"}", 0),
    # post hoc section
    ('G70', 'post hoc: tie band x2 -> x1', "dict(tie_scale=2.0)", "dict(tie_scale=1.0)", 0),
    ('G71', 'post hoc: ctrl_l1 only -> both', "dict(ctrls=('ctrl_l1',))", "dict(ctrls=('ctrl_l1', 'ctrl_l2'))", 0),
    ('G72', 'post hoc: margin x1.5 -> x1', "dict(margin_scale=1.5, a_params=False)", "dict(margin_scale=1.0, a_params=False)", 0),
    ('G73', 'post hoc: p < 0.10 -> 0.05', "dict(p_alpha=0.10)", "dict(p_alpha=0.05)", 0),
    ('G74', 'post hoc: not worse than the control: <= -> <', "<= TIE[kinds[ds]] + EPS) for m in ('ctrl_l1', 'ctrl_l2')}", "< TIE[kinds[ds]] + EPS) for m in ('ctrl_l1', 'ctrl_l2')}", 0),
    ('G75', 'post hoc: parameter ratio inverted', 'np.mean([pars(ds, "ours_gd") for ds in dss]) / np.mean([pars(ds, "ctrl_l1") for ds in dss])', 'np.mean([pars(ds, "ctrl_l1") for ds in dss]) / np.mean([pars(ds, "ours_gd") for ds in dss])', 0),
    ('G76', 'dev summary: the chosen configuration marked wrongly', '{"*" if pick[m] == c else ""}', '{"*" if pick[m] != c else ""}', 0),
    ('G77', 'dev summary: ratio to the best -> raw', "ratios = [mean[(c, d)] / max(min(mean[(c2, d)] for c2 in cfgs if (c2, d) in mean), 1e-12) for d in dss if (c, d) in mean]", "ratios = [mean[(c, d)] for d in dss if (c, d) in mean]", 0),
    ('G78', 'dev summary: mean params over all configurations', "par = float(np.mean([r['params'] for r in rm if r['cfg_id'] == c]))", "par = float(np.mean([r['params'] for r in rm]))", 0),
    ('G79', 'dev summary: first-round controls dropped', "    for m in sorted(set(r['method'] for r in dev if r['method'] not in GRIDS)):", "    for m in []:", 0),
]
