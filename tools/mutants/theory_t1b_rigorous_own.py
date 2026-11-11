# (id, what, old, new, k): single-point mutants of the changed lines of experiments/theory_t1b_rigorous.py (referee round T1b, findings 2, 5, 6): the hypotheses asserted inside
# `verify` (r > d_0, the block critical point in D-hat, theta < v_#), the strict bootstrap test E < D/1.0001, the upper-bound guards for undecided balls.
# Run: tools/mutants.py --exp theory_t1b_rigorous --tests tests/theory_t1b_test.py --raw results/synthetic/theory_t1b_raw.npz --md results/synthetic/theory_t1b.md --out-env TT1B_OUT
#      --extra results/synthetic/theory_t1b.md --spec tools/mutants/theory_t1b_rigorous_own.py --work <scratch dir> -k "not regenerated"
# (the md and raw are those of theory_t1b_wide.py, which does not import the mutated lines for its report; the criterion is the tests file alone)
# Result (33 mutants, criterion: the tests file alone, -k selects the verifier tests): 29 killed.  Survivors, equivalent on the stored data:
#   H00 control (unmutated);  G06 (T_2 as a plain log: the ratio is decided and > 1 on every tested target, d_0 ~ 0.4 against d-hat ~ 0.04, so log and the clamped log coincide; the straddling case is
#   covered by the unit test of `logpos_up`);  G08 (Y_0 by its lower end: the ball Y_0 is point-like, width 1e-30, so lower and upper ends coincide);  G10 (guard on the lower end of vbar - vup:
#   point-like ball, same);  G11 (T_3 of the high side dropped: T_3 is unchanged to 1e-5 in the tests, i.e. the high-side time never exceeds the low-side time on K = 1, 2, 3 targets).
M = [
 ('H00', 'control', "def pos_up(x):", "def pos_up(x):", 0),
 ('H01', "r > d_0 non-strict", "h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))", "h['r>d0'] = bool(float(lo(r)) >= max(prep['d0']))", 0),
 ('H02', "r > d_0 against the smallest offset", "h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))", "h['r>d0'] = bool(float(lo(r)) > min(prep['d0']))", 0),
 ('H03', "r > d_0 against the first offset", "h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))", "h['r>d0'] = bool(float(lo(r)) > prep['d0'][0])", 0),
 ('H04', "r > d_0 against the upper end of r", "h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))", "h['r>d0'] = bool(float(up(r)) > max(prep['d0']))", 0),
 ('H05', "r > d_0 always true", "h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))", "h['r>d0'] = True", 0),
 ('H06', "D-hat membership: offset bound reversed", "off <= float(lo(dh))", "off >= float(lo(dh))", 0),
 ('H07', "D-hat membership: offset against the upper end of d-hat", "off <= float(lo(dh))", "off <= float(up(dh))", 0),
 ('H08', "D-hat membership: lower coefficient bound dropped", "wv_ >= float(up(amin_f * Bx['vD_lo'])) and ", "", 0),
 ('H09', "D-hat membership: upper coefficient bound dropped", " and wv_ <= float(lo(amin_f * Bx['vD_hi']))", "", 0),
 ('H10', "D-hat membership: or instead of and", "inD = inD and off <= float(lo(dh))", "inD = inD or off <= float(lo(dh))", 0),
 ('H11', "D-hat membership: only the first unit", "for i, e in enumerate(ev):\n        off", "for i, e in enumerate(ev[:1]):\n        off", 0),
 ('H12', "D-hat membership: the offset against the wrong feature (the unit's own centre index)", "float(prep['cs'][e['feature']])", "float(prep['cs'][i])", 0),
 ('H13', "D-hat membership: coefficient bounds with the lower end of v_-", "wv_ >= float(up(amin_f * Bx['vD_lo']))", "wv_ >= float(lo(amin_f * Bx['vD_lo']))", 0),
 ('H14', "theta < v_#: non-certain end of the ball", "h['theta<v#'] = bool(theta < float(lo(vsharp)))", "h['theta<v#'] = bool(theta < float(up(vsharp)))", 0),
 ('H15', "theta < v_#: reversed", "h['theta<v#'] = bool(theta < float(lo(vsharp)))", "h['theta<v#'] = bool(theta > float(lo(vsharp)))", 0),
 ('H16', "theta = the larger of the two caps", "theta = min(T.THETA_REL * float(prep['a'].max()), 0.1 * amin_f)\n    h[", "theta = max(T.THETA_REL * float(prep['a'].max()), 0.1 * amin_f)\n    h[", 0),
 ('H17', "theta < v_# always true", "h['theta<v#'] = bool(theta < float(lo(vsharp)))", "h['theta<v#'] = True", 0),
 ('H18', "the hypotheses not entered in the conditions", "    cond.update(hypotheses(prep, par, Bx, r, dh, vsharp))\n", "", 0),
 ('H19', "the hypotheses called with the wrong radius (the box offset)", "cond.update(hypotheses(prep, par, Bx, r, dh, vsharp))", "cond.update(hypotheses(prep, par, Bx, dh, dh, vsharp))", 0),
 ('S01', "bootstrap non-strict", "if not (float(Etot) < float(Dv[j]) / 1.0001):", "if not (float(Etot) <= float(Dv[j])):", 0),
 ('S02', "bootstrap strict without the margin", "if not (float(Etot) < float(Dv[j]) / 1.0001):", "if not (float(Etot) < float(Dv[j])):", 0),
 ('S03', "bootstrap margin 1.00001", "float(Dv[j]) / 1.0001)", "float(Dv[j]) / 1.00001)", 0),
 ('S04', "bootstrap margin 1.001", "float(Dv[j]) / 1.0001)", "float(Dv[j]) / 1.001)", 0),
 ('S05', "bootstrap margin 1.0002", "float(Dv[j]) / 1.0001)", "float(Dv[j]) / 1.0002)", 0),
 ('G01', "pos_up takes the lower end", "return arb(max(float(up(x)), 0.0))", "return arb(max(float(lo(x)), 0.0))", 0),
 ('G02', "pos_up returns 0 for an undecided ball", "return arb(max(float(up(x)), 0.0))", "return arb(max(float(up(x)), 0.0)) if float(lo(x)) > 0 else arb(0)", 0),
 ('G03', "logpos_up takes the log of the lower end", "return arb(0).union(up(x).log())", "return arb(0).union(lo(x).log())", 0),
 ('G04', "logpos_up returns 0 for an undecided ball", "    u_ = float(up(x))\n    if u_ <= 1.0:", "    u_ = float(lo(x))\n    if u_ <= 1.0:", 0),
 ('G05', "logpos_up without the clamp at 1", "    if u_ <= 1.0:\n        return arb(0)\n    return arb(0).union(up(x).log())", "    return arb(0).union(up(x).log()) if u_ > 0 else arb(0)", 0),
 ('G06', "T_2 as a plain log (no clamp)", "T2 = logpos_up(ratio) / kd", "T2 = ratio.log() / kd", 0),
 ('G07', "T_2 dropped", "T2 = logpos_up(ratio) / kd", "T2 = arb(0) / kd", 0),
 ('G08', "Y_0 the lower end", "Y0 = pos_up(Y0) ", "Y0 = arb(max(float(lo(Y0)), 0.0)) ", 0),
 ('G09', "T_3 (low side) dropped", "T3lo = logpos_up(q1) / mV", "T3lo = arb(0) / mV", 0),
 ('G10', "the high-side relaxation guarded by the lower end", "if float(up(vbar - vup_)) > 0:", "if float(lo(vbar - vup_)) > 0:", 0),
 ('G11', "T_3 (high side) dropped", "T3hi = logpos_up(q2) / mV2", "T3hi = arb(0) / mV2", 0),
]
