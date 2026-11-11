# (id, what, old, new, k): single-point mutants of experiments/theory3y_sep.py (the closed-form sufficient condition of Theorem S6' for two features of equal amplitude, its box, Lemma phi, the exact-recursion
# oracle, the part-S run and its report).  Run:
#   venv/bin/python tools/mutants.py --exp theory3y_sep --tests tests/test_theory3y.py --raw results/synthetic/theory3y_raw.npz --md results/synthetic/theory3y.md --out-env T3Y_OUT \
#       --spec tools/mutants/theory3y_sep_own.py --work <scratch>/mut --slot b
# (kill criterion: the tests file alone; `--from-raw` is a mode of theory3y_run.py, not of this library: fr_rc is not informative).
# Result (41 entries, the control and 40 mutants; first pass 24 killed, then tests added for the box threshold in Delta (S02), the target geometry (S29, S30) and the benchmark rows (S39) and re-run with
#   --only: 28 of 40 killed; control passes).  Twelve survivors, equivalent on the data or redundant by construction: S03 ((C) fails => box_t1 returns before computing nu, so `'nu' in b` already excludes it);
#   S06 (Delta < 2 sqrt 2 never passes the box: it needs Delta - 2r >= sqrt 6); S10 (a_1 = a_2 or a_2 = h/2 make (S1) fail by themselves: P_- < phi(a_2) = phi(a_1) resp. = phi(h - a_2));
#   S14, S15 (the terms phi(h), phi(h - a_2) of L1 are never the maximum: phi(a_2) dominates them for a_1 < a_2 < h/2, so they are redundant for the decisions; they are kept in the theorem so that its proof needs
#   no monotonicity of phi below h); S20, S21, S22 (the competitor classes C2, C3, C4 are never the minimum: C1, the second-nearest kernel of the second feature, binds on all of the sampled region, which
#   `test_sep_benchmark_rows_and_the_competitor_classes` asserts; they are needed for the proof, not for the numbers); S26 (q_2 > 0 always: g(a_2) - rho(Delta - a_2) is above 0.9); S31 (the limit path of K = 2 has exactly two
#   events); S34 (sign_ok holds on all drawn targets); S35 (with two entries the final block is the target, Prop. A(iv)).
M = [
 ('S00', 'control (unmutated)', '', '', 0),
 ('S01', 'box: radius h/2 (not above the offsets)', "r = h / 2 + 0.02\n    b = T.box_t1(2, 1.0, Delta, r, 0.0, 0.05, 0.01)", "r = h / 2\n    b = T.box_t1(2, 1.0, Delta, r, 0.0, 0.05, 0.01)", 0),
 ('S02', 'box: no (N)', "if 'nu' in b and b['nu'] > 0 and b['vlo'] > 0 and b['Cmargin'] > 0 and b['vD_lo'] > 0:", "if 'nu' in b and b['vlo'] > 0 and b['Cmargin'] > 0 and b['vD_lo'] > 0:", 0),
 ('S03', 'box: no (C)', "and b['Cmargin'] > 0 and b['vD_lo'] > 0:", "and b['vD_lo'] > 0:", 0),
 ('S04', 'box: R = 2', "T.box_t1(2, 1.0, Delta, r, 0.0, 0.05, 0.01)", "T.box_t1(2, 2.0, Delta, r, 0.0, 0.05, 0.01)", 0),
 ('S05', 'lemma phi: threshold on Delta - sqrt 2', "T.psi(h) > T.psi(Delta - np.sqrt(2.0)) + 1e-12", "T.psi(h) > T.psi(Delta) + 1e-12", 0),
 ('S06', 'lemma phi: Delta below 2 sqrt 2 allowed', "if Delta < 2 * np.sqrt(2.0) or h >= np.sqrt(2.0):", "if h >= np.sqrt(2.0):", 0),
 ('S07', 'phi: no second term', "return g(D) + g(Delta - D)", "return g(D)", 0),
 ('S08', 'rho: no psi term', "return dt * T.psi_sup(xx) + Q * g(xx)", "return Q * g(xx)", 0),
 ('S09', 'rho: no g term', "return dt * T.psi_sup(xx) + Q * g(xx)", "return dt * T.psi_sup(xx)", 0),
 ('S10', 'conditions: offsets not required ordered', "or not (0 <= a1 < a2 < h / 2):", "or not (0 <= a1 <= a2 <= h / 2):", 0),
 ('S11', 'conditions: Q without E_v', "Q = 1 - float(g(dt)) + Ev", "Q = 1 - float(g(dt))", 0),
 ('S12', 'conditions: P_- with +', "Pm, Pp = float(g(a1) - g(Delta - a1)), float(g(a1) + g(Delta - a1))", "Pm, Pp = float(g(a1) + g(Delta - a1)), float(g(a1) + g(Delta - a1))", 0),
 ('S13', 'conditions: L1 without phi(a2)', "L1 = float(max(phi(a2, Delta), phi(h - a1, Delta), phi(h - a2, Delta), phi(h, Delta)))", "L1 = float(max(phi(h - a1, Delta), phi(h - a2, Delta), phi(h, Delta)))", 0),
 ('S14', 'conditions: L1 without phi(h)', "phi(h - a2, Delta), phi(h, Delta)))\n    out.update", "phi(h - a2, Delta)))\n    out.update", 0),
 ('S15', 'conditions: L1 without phi(h - a2)', "phi(h - a1, Delta), phi(h - a2, Delta), phi(h, Delta)))\n    out.update", "phi(h - a1, Delta), phi(h, Delta)))\n    out.update", 0),
 ('S16', 'conditions: S2 against rho(Delta)', "q2 = float(g(a2) - rho(Delta - a2))", "q2 = float(g(a2) - rho(Delta))", 0),
 ('S17', 'conditions: T_* without P_+', "(1 - (g(a2) - g(Delta - a2)) / Pp) / q2", "(1 - (g(a2) - g(Delta - a2)) / Pm) / q2", 0),
 ('S18', 'conditions: T_* numerator without the cross term', "(1 - (g(a2) - g(Delta - a2)) / Pp) / q2", "(1 - g(a2) / Pp) / q2", 0),
 ('S19', 'conditions: competitor C1 rate', "float(g(h - a2)) + rho(Delta - h + a2)),", "float(g(h - a2)) + rho(Delta)),", 0),
 ('S20', 'conditions: competitor C2 dropped', "             (float(phi(h - a1, Delta)) / Pm, float(g(Delta - h + a1)) + rho(h - a1)),\n", "", 0),
 ('S21', 'conditions: competitor C3 dropped', "             (float(phi(h, Delta)) / Pm, float(g(h)) + rho(Delta / 2)),\n", "", 0),
 ('S22', 'conditions: competitor C4 dropped', "             (float(phi(h, Delta)) / Pm, float(g(Delta / 2)) + rho(h))]", "             ]", 0),
 ('S23', 'conditions: strict versus margin of C1 (u_bar of the wrong kernel)', "(phi(h - a2, Delta)) / Pm,", "(phi(h - a1, Delta)) / Pm,", 0),
 ('S24', 'conditions: S3 inequality reversed', "out['s3'] = bool(Tstar < Tc)", "out['s3'] = bool(Tstar > Tc)", 0),
 ('S25', 'conditions: S1 not strict against L1', "out['s1'] = bool(Pm > L1)", "out['s1'] = bool(Pm > 0.5 * L1)", 0),
 ('S26', 'conditions: S2 not required', "out['s2'] = bool(q2 > 0)", "out['s2'] = True", 0),
 ('S27', 'conditions: ok without S3', "out['ok'] = bool(out['s1'] and out['s2'] and out['s3'])", "out['ok'] = bool(out['s1'] and out['s2'])", 0),
 ('S28', 'build_target: cells apart', "N = int(round(Delta / h))", "N = int(Delta / h)", 0),
 ('S29', 'build_target: amplitude sign of the second feature lost', "As = AMP * np.array([1.0, sgn_sep])", "As = AMP * np.array([1.0, 1.0])", 0),
 ('S30', 'build_target: flip ignored', "left_off, right_off = (sg1 * a1, sg2 * a2) if not flip else (sg2 * a2, sg1 * a1)", "left_off, right_off = (sg1 * a1, sg2 * a2)", 0),
 ('S31', 'oracle: entry count', "if len(ev) != 2:\n        return False", "if len(ev) < 2:\n        return False", 0),
 ('S32', 'oracle: order of the entries', "first = int(np.argmin(off))", "first = int(np.argmax(off))", 0),
 ('S33', 'oracle: uniqueness not required', "all(e['gap'] > 1e-9 for e in ev)", "all(e['gap'] > -1.0 for e in ev)", 0),
 ('S34', 'oracle: signs not required', "if not all(e['sign_ok'] for e in ev) or", "if False or", 0),
 ('S35', 'oracle: final block not checked', "return bool(np.allclose(np.sort(fin['c']), np.sort(cs), atol=1e-8) and np.allclose(np.sort(np.abs(fin['w'])), np.sort(np.abs(As)), atol=1e-8))", "return True", 0),
 ('S36', 'run_part: region counts only the diagonal side', "for j in range(i + 1, ng):", "for j in range(i, ng):", 0),
 ('S37', 'run_part: offsets sorted descending', "a1, a2 = np.sort(rng.uniform(0.0, h / 2, 2))", "a2, a1 = np.sort(rng.uniform(0.0, h / 2, 2))", 0),
 ('S38', 'run_part: violation counted without the condition', "n_viol += int(c['ok'] and not ok)", "n_viol += int(not ok)", 0),
 ('S39', 'run_part: benchmark offsets from the first kernels', "a = [float(np.min(np.abs(c0 - c))) for c in cs]", "a = [float(abs(c0[0] - c)) for c in cs]", 0),
 ('S40', 'report: violation column', "L.append('| %g | %g | %d | %d | %d | %d | **%d** |' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))", "L.append('| %g | %g | %d | %d | %d | %d | **%d** |' % (r[0], r[1], r[2], r[3], r[4], r[5], r[5]))", 0),
]
