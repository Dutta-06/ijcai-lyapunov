# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 111's revision (the earlier starts reaching the
# alarm, the coordinate's and the cold restart's conditions) and the referee's rerun-path survivors R01, R02.
M = [
 ('X00', 'control', "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", 0),
 ('R01', 'rerun path: the means one period late', "        W, C, src = Wp.tolist(), Cp.tolist(), 'rerun'", "        W, C, src = np.roll(Wp, -1, 0).tolist(), np.roll(Cp, -1, 0).tolist(), 'rerun'", 0),
 ('R02', "rerun path: the first period's centres", "        W, C, src = Wp.tolist(), Cp.tolist(), 'rerun'", "        W, C, src = Wp.tolist(), [Cp[0].tolist()] * len(Cp), 'rerun'", 0),
 ('T01', 'warm: not lengthened', "        row['warm'] = list(win_run(m, ke, mt, s, opts, 0.0, p - 1, NP=NP_C7 + 1,", "        row['warm'] = list(win_run(m, ke, mt, s, opts, 0.0, p - 1, NP=NP_C7,", 0),
 ('T02', 'warm_pk dropped', "        for k in (2, 3):\n            if p - k - 1 >= 0:", "        for k in ():\n            if p - k - 1 >= 0:", 0),
 ('T03', 'warm_pk with δ(p − k)', "                Wk[p - k - 1] = list(W[p - 1])", "                Wk[p - k - 1] = list(W[p - k - 1])", 0),
 ('T04', 'coord: both outcomes unchecked', " and any(seq) and not all(seq)\n", "\n", 0),
 ('T05', 'coord: offset 0 unchecked', "                    and zero is not None and zero[1] == x['reduced'])", "                    and True)", 0),
 ('T06', 'near: the exact period', "        near = lambda w, x: w[0] >= 0 and abs(w[0] - (x['reduced'] + 1)) <= 1", "        near = lambda w, x: w[0] == x['reduced'] + 1", 0),
 ('T07', 'warm_pk ignored', "                or any(neark(y[0], y[1:], x) for y in x.get('warm_pk', []))]", "]", 0),
 ('T08', 'one coordinate whatever the cold windows', "                       if cold and all((r, x) in nco for r, x in cold) else '')", "                       if cold else '')", 0),
]
