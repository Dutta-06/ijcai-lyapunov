# (id, what, old, new, k): single-point mutants of gh_lastann.py — round 126's fixes (the end link, the exact tiling, `_y_enclosed`'s
# preconditions, budget and true π/2, the chain's ends, the soundness corners). Run with --exp gh_lastann against tests/test_gh_peano.py
# -k part_va (gh_lastann writes no md: only the tests decide).
M = [
 ('R00', 'control', "def _chain(pieces, A, B):", "def _chain(pieces, A, B):", 0),
 ('R01', "the chain's end unchecked", "pieces[0][0] == A and pieces[-1][1] == B and all(", "pieces[0][0] == A and all(", 0),
 ('R02', 'an empty chain passing', "    return bool(pieces) and pieces[0][0] == A", "    return not pieces or pieces[0][0] == A", 0),
 ('R03', 'the end link unchecked', "            if not (abs(1 - acb(e_.real, e_.imag)) < END_R):", "            if False:", 0),
 ('R04', 'the tiling to 1e-12 again', "tiled=bool(exact and inside and disjoint))", "tiled=bool(rel < 1e-12 and inside and disjoint))", 0),
 ('R05', "the top edge's ball at the float π/2", "        Tb = hull(t1, arb.pi() / 2) if t2 == math.pi / 2 else hull(t1, t2)\n        om2", "        Tb = hull(t1, t2)\n        om2", 0),
 ('R06', "θ's upper end at the float π/2", "t_lo, t_hi = arb(t1), (arb.pi() / 2).upper() if t2 == math.pi / 2 else arb(t2)", "t_lo, t_hi = arb(t1), arb(t2)", 0),
 ('R07', 'the preconditions dropped', "    if not (0 < R1 <= R2 <= 2.0 and 0.0 <= f1 <= f2 <= math.pi / 2):\n        return None\n", "", 0),
 ('R08', "the fallback's lower end θ₁", "            y_lo, y_hi = arb(0), t_hi ", "            y_lo, y_hi = t_lo, t_hi ", 0),
 ('R09', "the expansion variable's lower end dropped", "            for e_ in ([e0] if e0 > 0 else []) + [e1]:", "            for e_ in [e1]:", 0),
 ('R10', "y_lo without θ₁", "            y_lo = k_lo if k_lo < t_lo else t_lo ", "            y_lo = k_lo ", 0),
]
