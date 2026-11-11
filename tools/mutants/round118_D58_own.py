# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 118's revision of A6X's verdict (the orders through eight
# terms decreasing, the limit against the last change, the transient's resolvability, B_K's flatness, its next order).
M = [
 ('T00', 'control', "    return q > 0 and q ** K > 10.0 ** (6 - digits)", "    return q > 0 and q ** K > 10.0 ** (6 - digits)", 0),
 ('T01', 'the old trigger q > ½', "    return q > 0 and q ** K > 10.0 ** (6 - digits)", "    return q > 0.5", 0),
 ('T02', 'decreasing: the first pair only', "all(abs(b_) < abs(a_) for a_, b_ in zip(z['orders'], z['orders'][1:]))", "all(abs(b_) < abs(a_) for a_, b_ in zip(z['orders'][:2], z['orders'][1:2]))", 0),
 ('T03', 'the last change unchecked', "    cons = all(abs(z['orders'][-1]) <= 10 * abs(z['orders'][-1] - z['orders'][-2]) for z in nofit)", "    cons = True", 0),
 ('T04', "B_K's flatness unchecked", " and z['tr']['flat'] < 1e-3\n", "\n", 0),
 ('T05', 'B_K over K = 120, 160 only', "for K in (60, 80, 120, 160) if K in R]", "for K in (120, 160) if K in R]", 0),
 ('T06', 'orders through seven terms', "for n in (3, 4, 5, 6, 7)]", "for n in (3, 4, 5, 6)]", 0),
]
