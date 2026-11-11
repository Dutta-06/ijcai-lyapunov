# (id, what, old, new, k): single-point mutants of gh_lastann.py — round 121's fixes (the NaN-safe handoff, the connector always
# certified, the joint's check, the soundness sampling, the tiling, the top edge's exact Im A = 0). Run with --exp gh_lastann against
# tests/test_gh_peano.py (gh_lastann has no md of its own: only the tests decide).
M = [
 ('L00', 'control', "def handoff_check(O, delta0=DELTA0, n=6, seed=121):", "def handoff_check(O, delta0=DELTA0, n=6, seed=121):", 0),
 ('L01', 'the NaN branch removed', "    if not (math.isfinite(kl) and math.isfinite(kh)):", "    if False:", 0),
 ('L02', 'the connector dropped', "    segs = [(complex(delta0, lo), complex(delta0, hi))]\n", "    segs = []\n", 0),
 ('L03', "the joint's containment unchecked", "            if not (Y >= lo and Y <= hi):", "            if False:", 0),
 ('L04', 'the soundness tolerance at 10³⁰', "if not (ex <= (0.0 if route == 'd' else 1e-30))", "if not (ex <= (0.0 if route == 'd' else 1e30))", 0),
 ('L05', 'the tiling unchecked', "tiled=bool(rel < 1e-12))", "tiled=True)", 0),
 ('L06', "the top edge's exact Im A = 0 dropped", "    exact0 = lambda th: th == 0 or (f2 == math.pi / 2 and th is top)", "    exact0 = lambda th: th == 0", 0),
 ('L07', "the NaN branch's connector a point", "        lo, hi = 0.0, th_hi ", "        lo, hi = 0.0, 0.0 ", 0),
 ('L08', 'the soundness sampling only at the upper end', "                e = arb(rnd.uniform(e0, e1) if e1 > e0 else e1)", "                e = arb(e1)", 0),
 ('L09', "B0's assertion removed", "            if not (a <= b.hi + arb(10) ** -50):                       # round 121", "            if False:                       # round 121", 0),
]
