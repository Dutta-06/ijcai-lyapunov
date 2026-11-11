# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 119's revision (the earlier starts with their own δ;
# the δ(p) starts as their control's drift; the samples merged into the scan; the late leaves labelled by their δ).
M = [
 ('Y00', 'control', "        def merged(x):", "        def merged(x):", 0),
 ('Y01', 'warm over every start', "        warm = [(r, x) for r, x in rows if any(neark(k, y, x) for k, y in starts(x, 'own'))]", "        warm = [(r, x) for r, x in rows if any(neark(k, y, x) for k, y in starts(x))]", 0),
 ('Y02', "the δ(p) reproductions dropped", "        dpo = [(r, x) for r, x in rows if (r, x) not in warm and any(neark(k, y, x) for k, y in starts(x, 'p'))]", "        dpo = []", 0),
 ('Y03', 'cold over every start', "not [1 for k, y in starts(x, 'own') if 0 <= y[0] < x['reduced'] + k - 1]]", "not [1 for k, y in starts(x) if 0 <= y[0] < x['reduced'] + k - 1]]", 0),
 ('Y04', "the δ(p) label unconditional", "(\" — the δ(p) starts' own drift (referee rounds 114, 119)\" if all(", "(\" — the δ(p) starts' own drift (referee rounds 114, 119)\" if True or all(", 0),
 ('Y05', 'merged without c1p', "(('c1', 'r+c1m'), ('c1i', 'r+c1i'), ('c1p', 'r+c1p')) if o_ in x and tag in x])", "(('c1', 'r+c1m'), ('c1i', 'r+c1i')) if o_ in x and tag in x])", 0),
 ('Y06', 'the keeping set a prefix', "            seq = [b_ for _, b_ in merged(x)]", "            seq = [not b_ for _, b_ in merged(x)]", 0),
 ('Y07', 'inside the bracket non-strict', "max(rem) < x[o_] < min(keep))", "max(rem) <= x[o_] <= min(keep))", 0),
 ('Y08', "warm_k started with δ(p)", "start=offsets_start(W, C, p - k, False, (0,)))))", "start=offsets_start(Wk, C, p - k, False, (0,)))))", 0),
 ('Y09', "warm_p's late leave labelled its own", "(1, x.get('warm_p'), 'δ(p)')", "(1, x.get('warm_p'), 'its own δ')", 0),
 ('Y10', "warm_k's late leave labelled δ(p)", "(z[0], z[1:], 'its own δ') for z in x.get('warm_k', [])]", "(z[0], z[1:], 'δ(p)') for z in x.get('warm_k', [])]", 0),
]
