# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 94's new code (c7s_lines' least raise and overlap,
# the repeat at a ± 1, the C8R misses' names, the C7L header and superset sentence).
M = [
 ('X00', 'control', 'C7S_D = (-1, 1) ', 'C7S_D = (-1, 1) ', 0),
 ('B01', 'repeat by then again', "abs(y['p'] + y['left'] - a) <= 1]", "y['p'] + y['left'] <= a]", 0),
 ('B02', 'repeat within two periods', "abs(y['p'] + y['left'] - a) <= 1]", "abs(y['p'] + y['left'] - a) <= 2]", 0),
 ('B03', 'Jaccard by intersection over the base', "jac = lambda a, b: len(a & b) / len(a | b) if a | b else 1.0", "jac = lambda a, b: len(a & b) / len(a) if a else 1.0", 0),
 ('B04', "the least raise's change by the base's over the neighbour's", "lr.append((max(least(b), least(q)) / min(least(b), least(q)), k_, q['opts'].get('shift', 0)))",
  "lr.append((least(b) / least(q), k_, q['opts'].get('shift', 0)))", 0),
 ('B05', 'far pairs 10 periods apart', "or kb[5] - ka[5] < 25:", "or kb[5] - ka[5] < 10:", 0),
 ('B06', "C8R's misses without the never-left", "miss = [r for r, y in xs if y['left'] < 0 or abs(y['left'] - r['left']) > 1]", "miss = [r for r, y in xs if abs(y['left'] - r['left']) > 1 and y['left'] >= 0]", 0),
 ('B07', 'the median overlap by the mean', "median {float(np.median(nbj)):.3f}", "median {float(np.mean(nbj)):.3f}", 0),
]
