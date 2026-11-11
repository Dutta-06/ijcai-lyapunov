# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 114's revision (a pass must pass the alarm's period;
# late leaves named; the samples checked against the scan; offset 0's distance; the stored sample offsets; the early starts' guards).
M = [
 ('K00', 'control', "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", 0),
 ('K01', 'early leaves counted as passes', "        cold = [(r, x) for r, x in rows if (r, x) not in warm and not early(x)]", "        cold = [(r, x) for r, x in rows if (r, x) not in warm]", 0),
 ('K02', 'early: a period too strict', "if 0 <= y[0] < x['reduced'] + k - 1]", "if 0 <= y[0] < x['reduced'] + k - 2]", 0),
 ('K03', 'late: a period too loose', "if y[0] > x['reduced'] + k + 1]", "if y[0] > x['reduced'] + k + 2]", 0),
 ('K04', 'coord without the samples', "and zero is not None and zero[1] == x['reduced'] and agree(x))", "and zero is not None and zero[1] == x['reduced'])", 0),
 ('K05', 'agree: the keep side unchecked', " or (keep and x[o_] >= min(keep) and x[tag][0] < 0)", "", 0),
 ('K06', "starts without 'warm_p'", " + ([(1, x['warm_p'])] if 'warm_p' in x else [])", "", 0),
 ('K07', 'neark exact', "        neark = lambda k, w, x: w[0] >= 0 and abs(w[0] - (x['reduced'] + k)) <= 1", "        neark = lambda k, w, x: w[0] >= 0 and abs(w[0] - (x['reduced'] + k)) <= 0", 0),
 ('K08', "the label's distance the scan's largest", "zd = lambda x: next((y[2] for y in x['scan'] if y[0] == 0.0), float('nan'))", "zd = lambda x: max(y[2] for y in x['scan'])", 0),
 ('K09', 'the warm guard p ≥ 1', "        if p >= 2:                                                      # round 114 (D3di-6)", "        if p >= 1:                                                      # round 114 (D3di-6)", 0),
 ('K10', 'the samples swapped', "        row['c1i'], row['c1p'] = float(Ci[0] - CT.CSTAR[0]), float(Cq[0] - CT.CSTAR[0])", "        row['c1i'], row['c1p'] = float(Cq[0] - CT.CSTAR[0]), float(Ci[0] - CT.CSTAR[0])", 0),
]
