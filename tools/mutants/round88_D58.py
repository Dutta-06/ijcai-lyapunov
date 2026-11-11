# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 88's new or changed code (git diff 78a3a04: a6_lines'
# second fit and its bracketing count, the cusp's label, b5_lines' censored cells) and the referee's round-88 survivors that change
# results on a rerun (D3ci-4: P01, P02, R01, its ids kept). (eq?) marks a mutant expected to change no stored result and no printed line.
M = [
 ('X00', 'control', 'GN_PREC = 80 ', 'GN_PREC = 80 ', 0),
 # the referee's survivors (D3ci-4)
 ('P01', "a6_error: gh_nodes' default 110 digits", "zs, ws = gh_nodes(K, prec)", "zs, ws = gh_nodes(K)", 1),
 ('P02', "gen_ratio: gh_nodes' default 110 digits", "zs, ws = gh_nodes(K, prec)", "zs, ws = gh_nodes(K)", 2),
 ('R01', 'rate_conj: the tie to the second argument', "return (a, 1) if a >= b else (b, 2)", "return (a, 1) if a > b else (b, 2)", 0),
 # a6_lines: the second fit
 ('F21', 'a6_lines: the K⁻¹ fit through K⁻¹, K⁻³', "[[1.0, K ** -1.0, K ** -2.0] for K in Ks[-3:]]", "[[1.0, K ** -1.0, K ** -3.0] for K in Ks[-3:]]", 0),
 ('F22', 'a6_lines: the K⁻¹ fit through the first three K', "[[1.0, K ** -1.0, K ** -2.0] for K in Ks[-3:]]), rr[-3:])",
  "[[1.0, K ** -1.0, K ** -2.0] for K in Ks[:3]]), rr[:3])", 0),
 ('F23', 'a6_lines: the K⁻¹ limit unscaled', "(x, L * x['rate'] - 1, y1, ex, L2 * x['rate'] - 1)", "(x, L * x['rate'] - 1, y1, ex, L2 - 1)", 0),
 ('F24', 'a6_lines: the farthest by the first fit', "(max(dev, key=lambda z: abs(z[4])))]", "(max(dev, key=lambda z: abs(z[1])))]", 0),
 ('F25', 'a6_lines: the bracket strict', "if z[1] * z[4] <= 0)", "if z[1] * z[4] < 0)", 0),
 ('F26', 'a6_lines: the bracket over every rule', "of the {len(far(a_))}\"", "of the {len(dev)}\"", 0),
 ('F27', 'a6_lines: the K⁻¹ range by the first fit', "+ \"; with a K⁻¹ term \" + rng([z[4] for z in far(a_)]",
  "+ \"; with a K⁻¹ term \" + rng([z[1] for z in far(a_)]", 0),
 ('F28', "a6_lines: the bracket counts the first fit's sign", "if z[1] * z[4] <= 0)", "if z[1] <= 0)", 0),
 # the cusp's label
 ('C01', "a6_lines: the cusp labelled by its argument", "lab = 'both arguments, the cusp' if at else",
  "lab = ('τ²/v' if x['arg'] == 1 else 'vσ*²/τ²') + ', the cusp' if at else", 0),
 # b5_lines' censored cells
 ('B01', 'b5_lines: censored at one accepted step', "('*' if y['acc'] <= 2 else '')", "('*' if y['acc'] <= 1 else '')", 0),
 ('B02', 'b5_lines: censored at three', "('*' if y['acc'] <= 2 else '')", "('*' if y['acc'] <= 3 else '')", 0),
 ('B03', 'b5_lines: the footnote never', "if any(y['acc'] <= 2 for x in dirs for y in x['rows']):", "if False:", 0),
]
