# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — D‴.58's tiers after D‴.61's part VA (round 121): proved iff a
# Mehler frame, the three tiers named, the frameless fast-branch count.
M = [
 ('V00', 'control', "def d61_step(t):", "def d61_step(t):", 0),
 ('V01', 'proved only below V again', "        fr = frame                                                      # proved iff a Mehler frame", "        fr = frame and x['rate'] < D61_RSTAR_V  # proved iff a Mehler frame", 0),
 ('V02', 'the tiers swapped', "if t < D61_RSTAR else ('V', D61_RSTAR_V) if t < D61_RSTAR_V else", "if t < D61_RSTAR_V else ('V', D61_RSTAR_V) if t < D61_RSTAR else", 0),
 ('V03', 'frameless count over every fast-branch rule', "    nfast = sum(1 for y in A6 if y['arg'] == 2 and not abs(1 + 2 / y['r'] ** 2 - 2 / float(y['v'])) < 1 - 1e-12)", "    nfast = sum(1 for y in A6 if y['arg'] == 2)", 0),
 ('V04', "VA's tier printed as V's", "frame, r* = {D61_RSTAR_VA:g};", "frame, r* = {D61_RSTAR_V:g};", 0),
]
