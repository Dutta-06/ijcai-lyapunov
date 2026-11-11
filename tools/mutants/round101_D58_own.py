# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 101's new code (out_e and the outward ranges, the
# distinct counts, the fast branch's outlier, the slow branch by t(v), a5's computed comparison, part V's tier of thresholds, the Mercer
# constant rounded up).
M = [
 ('X00', 'control', 'D61_RSTAR_V = 0.9900 ', 'D61_RSTAR_V = 0.9900 ', 0),
 ('N01', 'out_e: the direction by up alone', "(math.ceil(m * (1 - 1e-12)) if (v > 0) == up else math.floor(m * (1 + 1e-12)))", "(math.ceil(m * (1 - 1e-12)) if up else math.floor(m * (1 + 1e-12)))", 0),
 ('N02', 'out_e: no carry at 10', "if m >= 10:\n        m, e = m / 10, e + 1\n    return f\"{'-'", "if False:\n        m, e = m / 10, e + 1\n    return f\"{'-'", 0),
 ('N03', 'rgo: the lower end rounded up', "rgo = lambda vs: f\"{out_e(min(vs), False)}", "rgo = lambda vs: f\"{out_e(min(vs), True)}", 0),
 ('N04', 'dist: by name', "len({(z[0]['r'], round(float(z[0]['v']), 12)) for z in zs})", "len({(z[0]['r'], z[0]['name']) for z in zs})", 0),
 ('N05', 'dist: printed whenever', "if dist(zs) < len(zs) else ''", "if dist(zs) <= len(zs) else ''", 0),
 ('N06', "outlier: the second ratio by the first fit", "for w, j in ((w1, 1), (w2, 4)))", "for w, j in ((w1, 1), (w2, 1)))", 0),
 ('N07', 'outlier: the first ratio rounded up', "{dn_f(q1, 1)} and {dn_f(q2, 1)} times", "{up_f(q1, 1)} and {dn_f(q2, 1)} times", 0),
 ('N08', 'outlier: near rules too', "fz = lambda fr: [z for z in dev if fr_(z) == fr and z[0]['arg'] == 2 and abs(np.log(z[0]['v_hs'])) >= A6_NEAR]", "fz = lambda fr: [z for z in dev if fr_(z) == fr and z[0]['arg'] == 2]", 0),
 ('N09', 'slow list ascending', "key=lambda z: 1 / z[0]['rate'])", "key=lambda z: z[0]['rate'])", 0),
 ('N10', 'slow list: near rules too', "abs(np.log(z[0]['v_hs'])) >= A6_NEAR and 't' in z[0]], key", "'t' in z[0]], key", 0),
 ('N11', 'slow list: frame label inverted', "({'a frame' if fr_(z) else 'no frame'})", "({'a frame' if not fr_(z) else 'no frame'})", 0),
 ('N12', 'a5: compared with 1', "{'<' if max(1 / c for _, _, c, _ in lims) < D61_RSTAR else '≥'}", "{'<' if max(1 / c for _, _, c, _ in lims) < 1 else '≥'}", 0),
 ('N13', 'a5: t(v) to nearest', "{up_f(max(1 / c for _, _, c, _ in lims), 4)}", "{max(1 / c for _, _, c, _ in lims):.4f}", 0),
 ('N14', "tier 2: Mercer-zero at the first r*", "cross(mz, D61_RSTAR_V, hi=1e3)", "cross(mz, D61_RSTAR, hi=1e3)", 0),
 ('N15', 'tier 2: the bisection to 50', "rfv, rmv = cross(fock, D61_RSTAR_V, hi=1e3), cross(mz, D61_RSTAR_V, hi=1e3)", "rfv, rmv = cross(fock, D61_RSTAR_V), cross(mz, D61_RSTAR_V)", 0),
 ('N16', 'tier 2: the frame end rounded down', "{up_f(math.sqrt((1 + D61_RSTAR_V) / (1 - D61_RSTAR_V)))}", "{dn_f(math.sqrt((1 + D61_RSTAR_V) / (1 - D61_RSTAR_V)))}", 0),
 ('N17', 'the Mercer constant to nearest', "{up_e(c)} at τ/σ", "{c:.1e} at τ/σ", 0),
 ('N18', "part V's r* 0.9905", 'D61_RSTAR_V = 0.9900 ', 'D61_RSTAR_V = 0.9905 ', 0),
 ('N19', "tier 2: Fock's σ*/σ rounded up", "{dn_f(math.sqrt(1 + rfv * rfv))}", "{up_f(math.sqrt(1 + rfv * rfv))}", 0),
]
