# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — part VC (the corner's blow-up limit).
M = [
 ('X00', 'control', 'def corner_psi(x, y, phi):', 'def corner_psi(x, y, phi):', 0),
 ('V01', 'Ψ: the cubic over 2', "    return -(x * x - y * y) - 4 * phi * x * y + u * ((y + phi) ** 2 - u * u / 3)", "    return -(x * x - y * y) - 4 * phi * x * y + u * ((y + phi) ** 2 - u * u / 2)", 0),
 ('V02', 'floor: slope φ', "    return min(2 * phi * x, phi)", "    return min(phi * x, phi)", 0),
 ('V03', 'φ_dom: the cubic dropped', "        return float(np.max(-(us + 0.5) ** 2 + phi * phi * us - us ** 3 / 3))", "        return float(np.max(-(us + 0.5) ** 2 + phi * phi * us))", 0),
 ('V04', 'sector over |η|', "f\"{x['lo'] / (2 - x['R']):.4f}\"", "f\"{x['lo'] / x['R']:.4f}\"", 0),
 ('V05', 'the negativity unchecked', "    neg = all(p_['E'] < 0 for r in VC['rows'] for p_ in r['pts'])", "    neg = True", 0),
]
