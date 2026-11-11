# (id, what, old, new, k): single-point mutants of exp_moving_width.py — round 109's new code (the exits with exact and reversed sums,
# the conditions of the exit line, the 50-digit offset, the resolved odd part, the lemma's roots, the sup-norm half, the held criterion).
M = [
 ('X00', 'control', 'class _FsumBumps(MW.BumpsCW):', 'class _FsumBumps(MW.BumpsCW):', 0),
 ('S01', 'fsum: numpy sums for h', "math.fsum((rP * a)[:, j]), math.fsum((rP * b)[:, j])", "float(np.sum((rP * a)[:, j])), float(np.sum((rP * b)[:, j]))", 0),
 ('S02', "rev: not reversed", "    xs = xs[::-1] if summ == 'rev' else xs\n", "    xs = xs\n", 0),
 ('S03', 'offset from the double zero', "{r['wmul']: a2_ulps_exact(r['ratio'], x['w'], x['sig']) for r, x in held}", "{r['wmul']: a2_ulps(r['ratio'], x['w'], x['sig'])[1] for r, x in held}", 0),
 ('S04', 'data: times unchecked', "and all(abs(t[0] / r['a'][0][0] - 1) < 0.02 for t in", "and all(True or abs(t[0] / r['a'][0][0] - 1) < 0.02 for t in", 0),
 ('S05', 'arith: the exact sums unchecked', "arith = all(r['summ']['fsum'][1][0] is None and", "arith = all(True and", 0),
 ('S06', 'arith: the flip unchecked', "\n                    and r['summ']['rev'][1][1] == -r['sym'][0][1] for r in ex.values())", " for r in ex.values())", 0),
 ('S07', 'odd part: not resolved', "    if xs is None:                                                      # round 109", "    if False:                                                         # round 109", 0),
 ('S08', 'lemma: grid floors', "            if q is None or nxt is None:", "            if True:", 0),
 ('S09', 'mirror: the subspace distance not halved', "({bn / 2:.2g} from the subspace", "({bn:.2g} from the subspace", 0),
 ('S10', "held: the outputs' maximum", "if abs(r['rows'][-1]['dc']) < 1e-12 else \" — leaving", "if mx < 1e-12 else \" — leaving", 0),
]
