# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 113's revision of `a6x_lines` (successive orders from
# the 25-digit errors with t(v) exact, the second branch's transient, the window geometry, the verdict's conditions, the coverage count).
M = [
 ('U00', 'control', "    pw = lambda e: (lambda K: Fr(1, int(K) ** e))", "    pw = lambda e: (lambda K: Fr(1, int(K) ** e))", 0),
 ('U01', 't(v) from the double', "            tx = Fr(rate_conj(Decimal(repr(x['r'])), vv)[0]) if vv is not None else Fr(x['rate'])", "            tx = Fr(x['rate'])", 0),
 ('U02', 'the orders with a K⁻¹ term', "[one] + [pw(e) for e in range(2, n + 2)])[0] * tx - 1) for n in (3, 4, 5)]", "[one] + [pw(e) for e in range(1, n + 1)])[0] * tx - 1) for n in (3, 4, 5)]", 0),
 ('U03', 'the transient only above q = 0.9', "            if q > 0.5:", "            if q > 0.9:", 0),
 ('U04', 'q not squared', "            q = float((ho / tx) ** 2) if ho is not None else 0.0", "            q = float(ho / tx) if ho is not None else 0.0", 0),
 ('U05', 'convergence by 2 times', "abs(z['orders'][1]) < abs(z['orders'][0]) / 10 and abs(z['orders'][2]) < abs(z['orders'][1]) / 10", "abs(z['orders'][1]) < abs(z['orders'][0]) / 2 and abs(z['orders'][2]) < abs(z['orders'][1]) / 2", 0),
 ('U06', 'the worsening unchecked', "                and abs(z['orders'][0]) < abs(z['orders'][1]) < abs(z['orders'][2]) for z in trs)", "                for z in trs)", 0),
 ('U07', 'the decay against q unchecked', "all(abs(v / z['tr']['q'] - 1) < 2e-3 for v in z['tr']['dec'][1:])", "True", 0),
 ('U08', 'window geometry with K⁻⁵', "    g = [geom(Kw(z0['i0']), Kw(z0['iL']), 4, (2, 3)),", "    g = [geom(Kw(z0['i0']), Kw(z0['iL']), 5, (2, 3)),", 0),
 ('U09', 'nfast over every fast-branch rule', "    nfast = sum(1 for y in A6 if y['arg'] == 2 and not (abs(1 + 2 / y['r'] ** 2 - 2 / float(y['v'])) < 1 - 1e-12 and y['rate'] < D61_RSTAR_V))", "    nfast = sum(1 for y in A6 if y['arg'] == 2)", 0),
 ('U10', 'the residual against the five-term fit', "                c4 = efit(R, Ks[-4:], [one, pw(2), pw(3), pw(4)])", "                c4 = efit(R, Ks[-5:], [one, pw(2), pw(3), pw(4), pw(5)])", 0),
 ('U11', 'ratio_hp inverted', "    ratio_hp = lambda y: Fr(y['e0']) / Fr(y['e1']) if y.get('e0') and y.get('e1') else None", "    ratio_hp = lambda y: Fr(y['e1']) / Fr(y['e0']) if y.get('e0') and y.get('e1') else None", 0),
 ('U12', 'the transient through K = 160', "                ct = efit(R, Ks[-5:], [one, pw(2), pw(3), pw(4), lambda K: qf ** int(K)])", "                ct = efit(R, Ks[-4:], [one, pw(2), pw(3), lambda K: qf ** int(K)])", 0),
]
