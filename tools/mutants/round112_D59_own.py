# (id, what, old, new, k): single-point mutants of exp_moving_width.py — round 112's revision (the long-double exit with the offset as
# the state, the route's seed, the exit line's conditions decided in long double, the arithmetic's clause, the mirror evidence by width).
M = [
 ('V00', 'control', "def a2_exit_ld(ratio, wmul, sym=False, level=1e-2,", "def a2_exit_ld(ratio, wmul, sym=False, level=1e-2,", 0),
 ('V01', '_LDBumps in double', "        d, s = th[:, 0].astype(LD), th[:, 1].astype(LD)", "        d, s = th[:, 0].astype(float), th[:, 1].astype(float)", 0),
 ('V02', 'the seed without the linearisation', "    g = np.exp(-Lam) * dd", "    g = dd", 0),
 ('V03', 'the seeded run to 10%', "and rel(r['ld']['symS'], r['ld']['a'][1]) < 1e-4", "and rel(r['ld']['symS'], r['ld']['a'][1]) < 1e-1", 0),
 ('V04', "the seeded run's direction unchecked", "and ok_(r['ld']['symS']) and r['ld']['symS'][1] == r['ld']['a'][1][1] and", "and ok_(r['ld']['symS']) and", 0),
 ('V05', 'the arithmetic clause without the precision', "and len(late) == 2 * len(ex) and min(late) > 0.05", "and len(late) == 2 * len(ex) and min(late) > -1", 0),
 ('V06', 'the 2% against the double run', "and all(abs(t[0] / r['ld']['a'][1][0] - 1) < 0.02 for t in", "and all(abs(t[0] / r['a'][0][0] - 1) < 0.02 for t in", 0),
 ('V07', 'the flip in long double only', "flip = lambda r: r['sym'][0][1] != r['sym'][1][1] or r['ld']['sym'][0][1] != r['ld']['sym'][1][1]", "flip = lambda r: r['ld']['sym'][0][1] != r['ld']['sym'][1][1]", 0),
 ('V08', 'the long-double tolerances unchecked', "                   and rel(r['ld']['a'][0], r['ld']['a'][1]) < 1e-4\n", "\n", 0),
 ('V09', 'the width split at 2h', "ac & (np.abs(th_[:, 1]) >= XS[1] - XS[0]), ac & (np.abs(th_[:, 1]) < XS[1] - XS[0])", "ac & (np.abs(th_[:, 1]) >= 2 * (XS[1] - XS[0])), ac & (np.abs(th_[:, 1]) < 2 * (XS[1] - XS[0]))", 0),
 ('V10', 'the 1e-10 long-double run dropped', "dict(a=[a2_exit_ld(ratio, wm, rtol=1e-10), exl]", "dict(a=[exl, exl]", 0),
 ('V11', "A's exact and reversed exits' directions dropped", "    A_ = lambda r: r['a'] + [m[1] for m in r['imp'] if m[1][0] is not None] + [r['summ']['fsum'][0], r['summ']['rev'][0]]", "    A_ = lambda r: r['a'] + [m[1] for m in r['imp'] if m[1][0] is not None]", 0),
 ('V12', 'the seed by a one-sided difference', "    lam = np.array([(rhs(t, np.r_[y[0], hd, y[2]])[1] - rhs(t, np.r_[y[0], -hd, y[2]])[1]) / (2 * hd)", "    lam = np.array([(rhs(t, np.r_[y[0], hd, y[2]])[1] - rhs(t, np.r_[y[0], 0.0, y[2]])[1]) / hd", 0),
]
