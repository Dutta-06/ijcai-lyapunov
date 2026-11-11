# (id, what, old, new, k): single-point mutants of experiments/theory3y_run.py (the layer list, the loops over eps and v_0, the aligned windows, the report tables).  Run:
#   venv/bin/python tools/mutants.py --exp theory3y_run --tests tests/test_theory3y.py --raw results/synthetic/theory3y_raw.npz --md results/synthetic/theory3y.md --out-env T3Y_OUT \
#       --spec tools/mutants/theory3y_run_own.py --work <scratch>/mut --slot c
# (kill criterion: the tests file alone; the `md` column of the result is the --from-raw comparison and is informative here).
# Result (21 entries, the control and 20 mutants; first pass 9 killed, then tests added (a layer outside the regime is not integrated; the report's tables on tampered raw data; the quick report's rows)
#   and re-run with --only: 18 of 20 killed; control passes).  Two survivors, equivalent: R09 (the window of part D: 0.9 or 0.5 of the smallest t_# gives the same maximum, attained after the crossing);
#   R17 (adding a constant to log10 y does not change a fitted slope).
M = [
 ('R00', 'control (unmutated)', '', '', 0),
 ('R01', 'quick: the layer list', "return [(Y.layers_of(Y.D_TARGETS[0])[0], 'D'), (Y.stress_layers()[0], 'S')]", "return [(Y.layers_of(Y.D_TARGETS[0])[0], 'D'), (Y.stress_layers()[1], 'S')]", 0),
 ('R02', 'quick: eps list', "eps_l = [1e-6, 1e-10] if quick else Y.EPS_LIST", "eps_l = [1e-6, 1e-8] if quick else Y.EPS_LIST", 0),
 ('R03', 'quick: v0 list', "v0_l = [1e-3, 1e-6] if quick else Y.V0_LIST", "v0_l = [1e-3, 1e-5] if quick else Y.V0_LIST", 0),
 ('R04', 'regime flag', "regime=float(b is not None)", "regime=1.0", 0),
 ('R05', 'reference flow of the aligned distances', "ref = Y.integrate_layer(L, 0.0, v0=Y.V0_REF, w_hash=bd0['w_hash'])", "ref = Y.integrate_layer(L, 0.0, v0=1e-12, w_hash=bd0['w_hash'])", 0),
 ('R06', 'part C window', "lo = 0.9 * tr['t_hash']\n                C_win[i, j] = lo", "lo = 0.5 * tr['t_hash']\n                C_win[i, j] = lo", 0),
 ('R07', 'part C: bound of eps = 0 replaced by the layer\'s eps', "ok, vec = Y.check_layer(L, b, bd0, tr, 0.0)\n                C_ok", "ok, vec = Y.check_layer(L, b, bd0, tr, 1e-3)\n                C_ok", 0),
 ('R08', 'part B: bounds of eps = 0', "bd = Y.bounds(L, b, eps)\n                tr = Y.integrate_layer(L, eps, w_hash=bd['w_hash'])", "bd = Y.bounds(L, b, 0.0)\n                tr = Y.integrate_layer(L, eps, w_hash=bd['w_hash'])", 0),
 ('R09', 'part D window', "lo = 0.9 * min(trB[j]['t_hash'] for j in trB)", "lo = 0.5 * min(trB[j]['t_hash'] for j in trB)", 0),
 ('R10', 'part D: orbit start level', "max(vs_, 1e-12)", "max(vs_, 1e-3)", 0),
 ('R11', 'report: violation count', "int(sum(len(ii) - int(B_ok[ii, j].all(axis=1).sum()) for j in range(len(eps_l))))", "0", 0),
 ('R12', 'report: all-checks column', "int(okj.sum()), fe(fmax(rh))", "len(ii), fe(fmax(rh))", 0),
 ('R13', 'report: ratio of the hitting time', "rh = vv[:, vi['ratio_hat']]", "rh = vv[:, vi['ratio_decay']]", 0),
 ('R14', 'report: min t_# column', "fe(th.min()), fe(th.max()),", "fe(th.max()), fe(th.max()),", 0),
 ('R15', 'report: aligned vs orbit column of part C', "fe(fmax(dd)), fe(fmed(dd)), fe(fmax(oo)), fe(fmed(oo)),", "fe(fmax(oo)), fe(fmed(oo)), fe(fmax(dd)), fe(fmed(dd)),", 0),
 ('R16', 'report: regime shares by set', "int(reg[st == 'S'].sum()), int((st == 'S').sum()),", "int(reg[st == 'R'].sum()), int((st == 'R').sum()),", 0),
 ('R17', 'report: slope fit', "np.log10(x[ok]), np.log10(y[ok])", "np.log10(x[ok]), np.log10(y[ok]) + 1.0", 0),
 ('R18', 'run: layers outside the regime integrated', "        if b is None:\n            continue\n        bd0", "        if False:\n            continue\n        bd0", 0),
 ('R19', 'raw: checks array shape', "B_ok = np.zeros((nL, ne, len(CHK)), bool)", "B_ok = np.ones((nL, ne, len(CHK)), bool)", 0),
 ('R20', 'label format', "return '(%d,%g,%g,%d)/%d' % (t[0], t[1], t[2], t[3], L['k'] + 1)", "return '(%d,%g,%g,%d)/%d' % (t[0], t[1], t[2], t[3], L['k'])", 0),
]
