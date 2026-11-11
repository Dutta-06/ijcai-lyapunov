# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — round 99's new code (chord_q, the seed pixel's guard, path_ok's
# bool, dom_sector, part_dom's sect, d_lines' sector and chord lines).
M = [
 ('X00', 'control', 'DOM_SLACK = 1e-12 ', 'DOM_SLACK = 1e-12 ', 0),
 ('T01', 'chord_q: (|η|/2 − D)³', "return D * D / ((abs(eta) / 2 - D) ** 4 * abs(2 - 8 / eta ** 2))", "return D * D / ((abs(eta) / 2 - D) ** 3 * abs(2 - 8 / eta ** 2))", 0),
 ('T02', 'chord_q: |φ″| from 2 − 4/η²', "return D * D / ((abs(eta) / 2 - D) ** 4 * abs(2 - 8 / eta ** 2))", "return D * D / ((abs(eta) / 2 - D) ** 4 * abs(2 - 4 / eta ** 2))", 0),
 ('T03', 'chord_q: D for D²', "return D * D / ((abs(eta) / 2 - D) ** 4", "return D / ((abs(eta) / 2 - D) ** 4", 0),
 ('T04', 'path_check: q ≤ 1', "and l0 != le and q < 1), fv(a_)", "and l0 != le and q <= 1), fv(a_)", 0),
 ('T05', 'path_check: the chord not required', "and l0 != le and q < 1), fv(a_)", "and l0 != le), fv(a_)", 0),
 ('T06', 'path_check: D the ball radius', "q = chord_q(eta, abs(b_ - eta / 2))", "q = chord_q(eta, rb)", 0),
 ('T07', 'the pixel guard: ≥ rb', "if not abs(complex(X[i, j], Y[i, j]) - eta / 2) > rb:", "if not abs(complex(X[i, j], Y[i, j]) - eta / 2) > 0.5 * rb:", 0),
 ('T08', 'the pixel guard on a only', "    for z in (a_, b_):                                                  # round 99", "    for z in (a_,):                                                  # round 99", 0),
 ('T09', 'path_ok: the margin', "    return bool(path_check(eta, h)[0])", "    return path_check(eta, h)[1]", 0),
 ('T10', 'dom_sector: 20 bisections', "        for _ in range(80):\n            m = (a + b) / 2\n            a, b = (m, b) if (gap_closed(R, m) > 0) == sa", "        for _ in range(20):\n            m = (a + b) / 2\n            a, b = (m, b) if (gap_closed(R, m) > 0) == sa", 0),
 ('T11', 'dom_sector: the ends left out of Re', "re = [gap_arg(R, f)[2] for f in [lo] + [float(ph[i]) for i in idx] + [hi]]", "re = [gap_arg(R, f)[2] for f in [float(ph[i]) for i in idx]]", 0),
 ('T12', 'dom_sector: the grid ends for the crossings', "lo = cross(float(ph[idx[0] - 1]), float(ph[idx[0]])) if idx[0] > 0 else 0.0", "lo = float(ph[idx[0]])", 0),
 ('T13', 'part_dom: no sector', "margin=margin, sect=dom_sector())", "margin=margin)", 0),
 ('T14', 'd_lines: the sector rounded to nearest', "f\"({fl(sc['lo'], 5):.5f}, {cl(sc['hi'], 5):.5f}), and there", "f\"({sc['lo']:.5f}, {sc['hi']:.5f}), and there", 0),
 ('T15', 'd_lines: Re rounded to nearest', "f\"{fl(sc['re_min'], 4):.4f} and {cl(sc['re_max'], 4):.4f} (", "f\"{sc['re_min']:.4f} and {sc['re_max']:.4f} (", 0),
 ('T16', 'd_lines: the chord ratio min', "(f\", the chord's ratio q at most {max(x['q']):.1e}\"", "(f\", the chord's ratio q at most {min(x['q']):.1e}\"", 0),
]
