# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — round 127's seam check (`va_seams`) and the distinct soundness counts'
# printing. Run with --exp exp_gh_peano against tests/test_gh_peano.py -k "part_va or md_is_the_report".
M = [
 ('S00', 'control', "def va_seams(VA, V=None, D=None):", "def va_seams(VA, V=None, D=None):", 0),
 ('S01', 'the seam R_v == R1 dropped', "bool(R_bb > 1.0 and R1 is not None and R_v == R1)", "bool(R_bb > 1.0 and R1 is not None)", 0),
 ('S02', "D's resolution unchecked", "bool(R_bb > 1.0 and R1 is not None and R_v == R1)", "bool(R1 is not None and R_v == R1)", 0),
 ('S03', "V's failed bands counted", "        if b['fails'] == 0 and b['end'] < 0 and b.get('ball_worst', 0.0) < 0", "        if b['end'] < 0 and b.get('ball_worst', 0.0) < 0", 0),
 ('S04', "D's open boxes counted", "if x.get('open', 1) == 0 and x.get('ball_worst', 0.0) < 0]", "if x.get('ball_worst', 0.0) < 0]", 0),
 ('S05', "the annulus from a literal", "So the valley condition holds on {R1_va:g} ≤ |η| < 2", "So the valley condition holds on 1.995 ≤ |η| < 2", 0),
 ('S06', "V and D not passed to va_lines", "+ va_lines(VA, V, D) +", "+ va_lines(VA) +", 0),
]
