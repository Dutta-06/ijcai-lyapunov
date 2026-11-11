"""Single-point mutants of experiments/exp_kan_merge.py (registries, selection of the control row, merging rule, strata, duplicate statistics, matching, bootstrap, oracle
selection, report), run with

    venv/bin/python tools/mutants.py --exp exp_kan_merge --tests tests/test_exp_kan_merge.py --raw results/ml/kan_raw.npz --md results/ml/kan_merge.md \\
        --out-env KAN_OUT --spec tools/mutants/kan_merge_own.py --work <scratch>/mut --slot a --extra results/ml/kan.md --extra results/ml/kan_merge.md

(id, what, old, new, k): `old` replaced by `new` at its k-th occurrence. M00 is the control (no change). The kill criterion is the tests file alone.

Result (58 entries: the control and 57 mutants): 53 killed by tests/test_exp_kan_merge.py (16 tests, about 1.5 s), the control passes. Four survivors, all equivalent on the stored data:
  S06  `o = o[0]` -> `o[-1]`: every complete run holds exactly one 'ours' record, so first and last are the same record.
  P02  `<= 0.0` -> `< 0.0` in the "point estimate first not positive" delta: no gain of the grid is exactly 0 (the gain is a mean of 10 target gains in 1/30 steps), so the two rules give the same delta.
  P03  `>=` -> `>` against the plan's 0.10 floor of the gain: no gain on the grid equals 0.10 exactly (0.187, 0.097, ...).
  P07  the "match: yes" check of the headline reproduction replaced by True: the check is always true on the stored raw (it is a printed self-check, asserted by the tests through the headline numbers).
"""
M = [
    ('M00', 'control (no change)', '', ''),
    # ---- registries
    ('R01', 'delta grid', "DELTAS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)", "DELTAS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.6)"),
    ('R02', 'gap thresholds', "GAP_THRESH = (0.1, 0.3)", "GAP_THRESH = (0.1, 0.4)"),
    ('R03', 'match tolerances', "MATCH_TOLS = (0.2, 0.5)", "MATCH_TOLS = (0.2, 0.4)"),
    ('R04', 'bootstrap resamples', "N_BOOT, BOOT_SEED = 20000, 0", "N_BOOT, BOOT_SEED = 10000, 0"),
    ('R05', 'bootstrap seed', "N_BOOT, BOOT_SEED = 20000, 0", "N_BOOT, BOOT_SEED = 20000, 1"),
    # ---- selection (the spike module's rules)
    ('S01', 'control by lowest validation MSE instead of one-SE', "cv = K.one_se_pick(cands, ctrl[0]['n_val'])['r']", "cv = min(ctrl, key=lambda r: r['val_mse'])"),
    ('S02', 'oracle tie rule dropped', "co = max(ctrl, key=lambda r: (K.exact_share(r, tau), -r['P_%g' % tau]))", "co = max(ctrl, key=lambda r: K.exact_share(r, tau))"),
    ('S03', 'oracle tie by more parameters', "co = max(ctrl, key=lambda r: (K.exact_share(r, tau), -r['P_%g' % tau]))", "co = max(ctrl, key=lambda r: (K.exact_share(r, tau), r['P_%g' % tau]))"),
    ('S04', 'incomplete runs scored', "if not o or not ctrl or not any(r['method'] == 'fixed' for r in rs):", "if not o or not ctrl:"),
    ('S05', 'only synthetic records dropped', "if r['kind'] == 'synth':\n            runs.setdefault", "if True:\n            runs.setdefault"),
    ('S06', 'ours is the last record', "o = o[0]\n        cands", "o = o[-1]\n        cands"),
    ('S07', 'ctrl_all truncated', "ctrl_all=ctrl))", "ctrl_all=ctrl[:5]))"),
    # ---- merging
    ('G01', 'gaps not sorted', "s = sorted_centres(cs)\n    return [round(b - a, 6)", "s = [float(c) for c in cs]\n    return [round(abs(b - a), 6)"),
    ('G02', 'gap rounding', "round(b - a, 6)", "round(b - a, 1)"),
    ('G03', 'merge when gap <= delta', "return 1 + sum(1 for g in gaps(cs) if g >= delta)", "return 1 + sum(1 for g in gaps(cs) if g > delta)"),
    ('G04', 'merge only at twice the radius', "return 1 + sum(1 for g in gaps(cs) if g >= delta)", "return 1 + sum(1 for g in gaps(cs) if g >= 2 * delta)"),
    ('G05', 'count of gaps not clusters', "return 1 + sum(1 for g in gaps(cs) if g >= delta)", "return sum(1 for g in gaps(cs) if g >= delta)"),
    ('G06', 'empty edge counted 1', "    if n == 0:\n        return 0\n", "    if n == 0:\n        return 1\n"),
    ('G07', 'delta 0 merges identical centres', "if g >= delta)", "if g > delta)"),
    ('G08', 'complete linkage proxy: span', "return 1 + sum(1 for g in gaps(cs) if g >= delta)", "return 1 + int(max(sorted_centres(cs)[-1] - sorted_centres(cs)[0], 0) // max(delta, 1e-9)) if delta > 0 else n"),
    ('G09', 'counts matrix of oracle uses ours', "return np.array([[merged_count(cs, delta) for cs in r[which][CKEY]] for r in runs])", "return np.array([[merged_count(cs, delta) for cs in r['ours'][CKEY]] for r in runs])"),
    ('G10', 'consistency check weakened', "assert [len(c) for c in rec[CKEY]] == list(rec[NKEY]), (r['key'], rec['conf'])", "assert True"),
    ('G11', 'consistency counts only ours', "for rec in [r['ours']] + r['ctrl_all']:", "for rec in [r['ours']]:"),
    # ---- oracle with merge
    ('O01', 'oracle picks the worst', "best = max(r['ctrl_all'], key=lambda rec: (float(np.mean(counts_matrix_rec(rec, delta) == t)), -rec['P_%g' % K.TAU]))",
     "best = min(r['ctrl_all'], key=lambda rec: (float(np.mean(counts_matrix_rec(rec, delta) == t)), -rec['P_%g' % K.TAU]))"),
    ('O02', 'oracle merged at delta 0', "best = max(r['ctrl_all'], key=lambda rec: (float(np.mean(counts_matrix_rec(rec, delta) == t))", "best = max(r['ctrl_all'], key=lambda rec: (float(np.mean(counts_matrix_rec(rec, 0.0) == t))"),
    ('O03', 'oracle tie by more parameters', "-rec['P_%g' % K.TAU]))\n        out.append", "rec['P_%g' % K.TAU]))\n        out.append"),
    ('O04', 'oracle returns unmerged counts', "out.append(counts_matrix_rec(best, delta))", "out.append(counts_matrix_rec(best, 0.0))"),
    # ---- bootstrap
    ('B01', 'bootstrap over runs not targets', "tids = np.array([r['key'][0] for r in runs])", "tids = np.array([hash(r['key']) for r in runs])"),
    ('B02', 'gain sign', "d = np.asarray(share_a, float) - np.asarray(share_b, float)", "d = np.asarray(share_b, float) - np.asarray(share_a, float)"),
    ('B03', 'gain not averaged within target', "by = [float(np.mean(d[tids == t])) for t in sorted(set(tids.tolist()))]", "by = [float(d[tids == t][0]) for t in sorted(set(tids.tolist()))]"),
    ('B04', 'boot called with default settings', "return K.boot_ci(by, n_boot=N_BOOT, seed=BOOT_SEED)", "return K.boot_ci(by)"),
    ('B05', 'gain of ctrl_best vs ctrl', "gain={nm: gain_ci(runs, sh['ours'], sh[nm]) for nm in ('ctrl', 'ctrl_or', 'ctrl_best')}", "gain={nm: gain_ci(runs, sh['ours'], sh['ctrl']) for nm in ('ctrl', 'ctrl_or', 'ctrl_best')}"),
    # ---- strata
    ('T01', 'empty mask nonzero', "(('empty', T == 0), ('nonempty', T > 0))", "(('empty', T == 0), ('nonempty', T >= 0))"),
    ('T02', 'over counted in under', "over = np.maximum(C - T, 0) * mask\n        under = np.maximum(T - C, 0) * mask", "over = np.maximum(T - C, 0) * mask\n        under = np.maximum(C - T, 0) * mask"),
    ('T03', 'surplus per edge divided by runs', "over_edge=float(over.sum() / n) if n else float('nan')", "over_edge=float(over.sum() / T.shape[0]) if n else float('nan')"),
    ('T04', 'share over all edges', "share=float(ex[mask].mean()) if n else float('nan')", "share=float(ex.mean()) if n else float('nan')"),
    ('T05', 'surplus not masked', "over = np.maximum(C - T, 0) * mask", "over = np.maximum(C - T, 0)"),
    ('T06', 'all5 uses any', "all5={nm: float(np.mean((C[nm] == T).all(1))) for nm in C}", "all5={nm: float(np.mean((C[nm] == T).any(1))) for nm in C}"),
    ('T07', 'nonempty share over all edges', "nonempty={nm: float(((C[nm] == T)[T > 0]).mean()) for nm in C}", "nonempty={nm: float(((C[nm] == T)).mean()) for nm in C}"),
    ('T08', 'units per run summed over runs', "units={nm: float(C[nm].sum(1).mean()) for nm in C}", "units={nm: float(C[nm].sum(1).sum()) for nm in C}"),
    # ---- duplicates and matching
    ('D01', 'gap threshold strict', "return len(g), [float(np.mean(g < t)) for t in GAP_THRESH]", "return len(g), [float(np.mean(g <= t)) for t in GAP_THRESH]"),
    ('D02', 'gaps of ours always', "g = [x for r in runs for cs in r[which][CKEY] for x in gaps(cs)]", "g = [x for r in runs for cs in r['ours'][CKEY] for x in gaps(cs)]"),
    ('D03', 'match recall and precision swapped', "return float(m[:, 0].mean()), float(m[:, 1].mean())", "return float(m[:, 1].mean()), float(m[:, 0].mean())"),
    ('D04', 'match tolerance ignored', "K.centre_match(dict(r[which], truth_c=r['truth_c']), tol)", "K.centre_match(dict(r[which], truth_c=r['truth_c']))"),
    ('D05', 'seed sd uses pooled runs', "return float(np.mean([np.std(v) for v in cells.values()]))", "return float(np.std([s for v in cells.values() for s in v]))"),
    # ---- report
    ('P01', 'first delta uses upper bound', "d_ci = first_delta(sw, lambda s: s['gain']['ctrl'][1] <= 0.0)", "d_ci = first_delta(sw, lambda s: s['gain']['ctrl'][2] <= 0.0)"),
    ('P02', 'first delta point estimate', "d_pt = first_delta(sw, lambda s: s['gain']['ctrl'][0] <= 0.0)", "d_pt = first_delta(sw, lambda s: s['gain']['ctrl'][0] < 0.0)"),
    ('P03', 'plan criterion floor', "d_s = first_delta(sw, lambda s: not (s['gain']['ctrl'][0] >= K.GO_SHARE_GAIN and s['gain']['ctrl'][1] > 0.0))", "d_s = first_delta(sw, lambda s: not (s['gain']['ctrl'][0] > K.GO_SHARE_GAIN and s['gain']['ctrl'][1] > 0.0))"),
    ('P04', 'worse-than share inverted', "worse = float(np.mean([r['ours']['test_rmse'] > r['ctrl']['test_rmse'] for r in runs]))", "worse = float(np.mean([r['ours']['test_rmse'] < r['ctrl']['test_rmse'] for r in runs]))"),
    ('P05', 'cap counted with >', "cap = int(sum(r['ours']['steps'] >= K.STEPS for r in runs))", "cap = int(sum(r['ours']['steps'] > K.STEPS for r in runs))"),
    ('P06', 'rmse of oracle for ctrl', "rm = {nm: float(np.mean([r[nm]['test_rmse'] for r in runs])) for nm, _ in METHODS}", "rm = {nm: float(np.mean([r['ctrl_or']['test_rmse'] for r in runs])) for nm, _ in METHODS}"),
    ('P07', 'reproduction check loosened', "all(abs(shs[nm] - pub[nm]) < 1e-12 for nm in pub)", "True"),
    ('P08', 'ties count loosened', "abs(x - d) < 1e-9)", "abs(x - d) < 1e-1)"),
    ('P09', 'per-K column at the wrong delta', "*(fmt(float((sw[0.3]['C'][nm] == T)[m].mean())) for nm in ('ours', 'ctrl'))", "*(fmt(float((sw[0.2]['C'][nm] == T)[m].mean())) for nm in ('ours', 'ctrl'))"),
    ('P10', 'default output name', "os.path.join(outd, 'kan_merge.md')", "os.path.join(outd, 'kan.md')"),
    ('P11', 'default raw name', "os.path.join(outd, 'kan_raw.npz'))\n    txt", "os.path.join(outd, 'kan_dev_raw.npz'))\n    txt"),
    ('P12', 'verdict from arch A', "ver = K.verdict(rows0, K.real_summary(recs, K.PRIMARY))", "ver = K.verdict(rows0, K.real_summary(recs, 'A'))"),
]
