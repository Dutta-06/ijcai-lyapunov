# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 105's new code (kernel 1's centre alone, the scan, the
# warm start, the probe's 8 subsets, the stricter conditions, the C8R aggregate).
M = [
 ('X00', 'control', "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", "C8A_SCAN = (-0.06, -0.04, -0.03, -0.02, 0.0, 0.04)", 0),
 ('Q01', 'centre alone: from the fitted coefficient', "        c1only = (('mean', 'mean'), ('fit', 'fit'), ('fit', 'fit'))", "        c1only = (('fit', 'mean'), ('fit', 'fit'), ('fit', 'fit'))", 0),
 ('Q02', '(p − 1)T: the state at pT', "        Cq = [float(v) for v in r['c'][(p - 1) * kk, idx]]", "        Cq = [float(v) for v in r['c'][p * kk, idx]]", 0),
 ('Q03', 'scan: offset from the fitted centre', "            Cm[p - 1][0] = float(CT.CSTAR[0] + off)", "            Cm[p - 1][0] = float(C[p - 1][0] + off)", 0),
 ('Q04', 'warm: from the fitted state', "start=offsets_start(W, C, p - 1, False, (0,))))", "start=None))", 0),
 ('Q05', 'warm: kept at the reduced period', "warm = sum(1 for _, x in rows if x['warm'][0] == x['reduced'] + 1)", "warm = sum(1 for _, x in rows if x['warm'][0] == x['reduced'])", 0),
 ('Q06', 'monotone unchecked', "        nmono = sum(mono(x) for _, x in rows)", "        nmono = len(rows)", 0),
 ('Q07', "centre clause: classes", "lpv('w+c1') == lpv('wc') and lpv('w+c23') == lpv('w') and oc['wc'] != oc['w']", "oc['w+c1'] == oc['wc'] and oc['w+c23'] == oc['w'] and oc['wc'] != oc['w']", 0),
 ('Q08', "centre clause: without wc ≠ w", " and oc['wc'] != oc['w'])", ")", 0),
 ('Q09', "needs: subsets unchecked", " and {(), (0,), (1,), (0, 1)} <= {tuple(o) for o, _ in has3}", "", 0),
 ('Q10', 'C8R aggregate unchecked', " and C8R is not None and not any(\n                                   y['left'] >= 0 and y['left'] - r_['left'] < -1 for r_ in C8R for y in r_['rows'])", "", 0),
 ('Q11', 'probe: 5 subsets', "((), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)))", "((2,), (0,), (0, 1), (0, 2), (0, 1, 2)))", 0),
 ('Q12', 'cold restart unconditional', '(" — the alarms are the cold restart\'s" if warm == 0 else \'\')', '" — the alarms are the cold restart\'s"', 0),
]
