# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 96's new code (up_e, dn_f, up_f, the split by frame ×
# branch × nearness, the inward thresholds, the labels from D61_RSTAR, the Mercer frame's proved label).
M = [
 ('X00', 'control', 'D61_RSTAR = 0.8052 ', 'D61_RSTAR = 0.8052 ', 0),
 ('L01', 'up_e: round to nearest', "m = math.ceil(v / 10.0 ** e * 10 ** d * (1 - 1e-12)) / 10 ** d", "m = round(v / 10.0 ** e * 10 ** d) / 10 ** d", 0),
 ('L02', 'up_e: no carry at 10', "    if m >= 10:\n        m, e = m / 10, e + 1\n", "", 0),
 ('L03', 'up_e: the guard 1e-6', "(1 - 1e-12)) / 10 ** d", "(1 - 1e-6)) / 10 ** d", 0),
 ('L04', 'dn_f: ceil', "return f\"{math.floor(v * 10 ** d) / 10 ** d:.{d}f}\"", "return f\"{math.ceil(v * 10 ** d) / 10 ** d:.{d}f}\"", 0),
 ('L05', 'up_f: floor', "return f\"{math.ceil(v * 10 ** d) / 10 ** d:.{d}f}\"", "return f\"{math.floor(v * 10 ** d) / 10 ** d:.{d}f}\"", 0),
 ('L06', 'cells: nearness inverted', "(abs(np.log(z[0]['v_hs'])) >= A6_NEAR) == nf]", "(abs(np.log(z[0]['v_hs'])) < A6_NEAR) == nf]", 0),
 ('L07', 'cells: the frame flag ignored', "zs = [z for z in dev if fr_(z) == fz and z[0]['arg'] == a_", "zs = [z for z in dev if z[0]['arg'] == a_", 0),
 ('L08', 'cells: the bracket by the first fit alone', "bracketing at {sum(1 for z in zs if z[1] * z[4] <= 0)}\")", "bracketing at {sum(1 for z in zs if z[1] <= 0)}\")", 0),
 ('L09', 'cells: the first fit for the second', "f\"{up_e(max(abs(z[4]) for z in zs))}, bracketing at", "f\"{up_e(max(abs(z[1]) for z in zs))}, bracketing at", 0),
 ('L10', 'cells: the branch label swapped', "for a_, bl in ((1, 'τ²/v'), (2, 'vσ*²/τ²')):\n                for nl", "for a_, bl in ((2, 'τ²/v'), (1, 'vσ*²/τ²')):\n                for nl", 0),
 ('L11', 'thresholds: the Fock σ*/σ to nearest', "σ*/σ < {dn_f(math.sqrt(1 + rf * rf))}", "σ*/σ < {math.sqrt(1 + rf * rf):.3f}", 0),
 ('L12', 'thresholds: the end to nearest', "σ*/σ > {up_f(math.sqrt((1 + D61_RSTAR) / (1 - D61_RSTAR)))}", "σ*/σ > {math.sqrt((1 + D61_RSTAR) / (1 - D61_RSTAR)):.3f}", 0),
 ('L13', 'the Mercer label: t(v) as conj', "most {max(1 / c for _, _, c, _ in lims):.4f} < r*", "most {max(c for _, _, c, _ in lims):.4f} < r*", 0),
 ('L14', 'the cusp bound to nearest', "(f\"; at v_hs within {up_e(max(abs(z[1]) for z in cusp))} of 1/t\"", "(f\"; at v_hs within {max(abs(z[1]) for z in cusp):.1e} of 1/t\"", 0),
 ('L15', "the frame rules' second bound to nearest", "the second \"\n                            f\"within {up_e(max(abs(z[4]) for z in zs))}", "the second \"\n                            f\"within {max(abs(z[4]) for z in zs):.1e}", 0),
 ('L16', 'the A6 heading literal', "where a rule has a frame and t(v) < ' + f\"{D61_RSTAR:g}\" + ', '", "where a rule has a frame and t(v) < 0.8053, '", 0),
]
