# Anchors of QL2, QL3, QR1, QR2 refreshed after round 91's part D moved report and q_lines' neighbour.
# Referee round 91 (D3cl): single-point mutants of part Q (new in a02bbba) and of round-90 code the round's specs leave out
M = [
 ('Z00', 'control', 'Q_K, Q_MU, Q_TH = ', 'Q_K, Q_MU, Q_TH = ', 0),
 # corner_ratio
 ('QA1', 'corner_ratio default he: He_{K-2}', "he = [WM.he_all(z, K)[K - 1] for z in zs] if he is None else he", "he = [WM.he_all(z, K)[K - 2] for z in zs] if he is None else he", 0),
 ('QA2', 'corner_ratio: |u|^K', "mag = r ** (K - 1) * ((x * x - y * y) / 2).exp()", "mag = r ** K * ((x * x - y * y) / 2).exp()", 0),
 ('QA3', 'corner_ratio: e^{|u|^2/2}', "mag = r ** (K - 1) * ((x * x - y * y) / 2).exp()", "mag = r ** (K - 1) * ((x * x + y * y) / 2).exp()", 0),
 ('QA4', 'corner_ratio: phase -xy', "ph = (K - 1) * th + x * y ", "ph = (K - 1) * th - x * y ", 0),
 ('QA5', 'corner_ratio: phase K th', "ph = (K - 1) * th + x * y ", "ph = K * th + x * y ", 0),
 ('QA6', 'corner_ratio: Im part sign', "qi = sum((w * h * (x * z).exp() * WM.dcos(y * z - pi / 2, pi)", "qi = sum((w * h * (x * z).exp() * WM.dcos(y * z + pi / 2, pi)", 0),
 ('QA7', 'corner_ratio: no conjugate in the real part', "return float((qr * dr + qi * di) / den), float((qi * dr - qr * di) / den)", "return float((qr * dr - qi * di) / den), float((qi * dr - qr * di) / den)", 0),
 ('QA8', 'corner_ratio: theta = j pi/4', "th = pi * j8 / 8", "th = pi * j8 / 4", 1),
 ('QA9', 'corner_ratio: u = mu K (not sqrt K)', "r = Decimal(repr(mu)) * Decimal(K).sqrt()", "r = Decimal(repr(mu)) * Decimal(K)", 0),
 ('QA10', 'corner_ratio: x, y swapped', "x, y = r * WM.dcos(th, pi), r * WM.dcos(th - pi / 2, pi)", "y, x = r * WM.dcos(th, pi), r * WM.dcos(th - pi / 2, pi)", 0),
 ('QA11', 'corner_ratio: Im part from the real weights only (no h)', "qi = sum((w * h * (x * z).exp() * WM.dcos(y * z - pi / 2, pi) for w, z, h in zip(ws, zs, he))", "qi = sum((w * (x * z).exp() * WM.dcos(y * z - pi / 2, pi) for w, z, h in zip(ws, zs, he))", 0),
 # q_job, part_Q
 ('QJ1', 'q_job: prec 3K + 20', "prec = 3 * K + 80", "prec = 3 * K + 20", 0),
 ('QJ2', 'q_job: limit 1 + mu^2', "dev.append(float(abs(c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))", "dev.append(float(abs(c / (1 + (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))", 0),
 ('QJ3', 'q_job: limit at conj(mu)', "dev.append(float(abs(c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))", "dev.append(float(abs(c / (1 - (m * cmath.exp(-1j * j8 * np.pi / 8)) ** 2) - 1)))", 0),
 ('QJ4', 'q_job: He_K instead of He_{K-1}', "he = [WM.he_all(z, K)[K - 1] for z in nodes[0]]", "he = [WM.he_all(z, K)[K] for z in nodes[0]]", 0),
 ('QJ5', 'q_job: the real part of the deviation', "dev.append(float(abs(c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))", "dev.append(float(abs((c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1).real)))", 0),
 ('QJ6', 'q_job: the limit 1 - |mu|^2', "dev.append(float(abs(c / (1 - (m * cmath.exp(1j * j8 * np.pi / 8)) ** 2) - 1)))", "dev.append(float(abs(c / (1 - m ** 2) - 1)))", 0),
 ('QP1', 'part_Q: angles of LC', "pool.map(q_job, [(K, j8) for K in Q_K for j8 in Q_TH]", "pool.map(q_job, [(K, j8) for K in Q_K for j8 in LC_TH]", 0),
 ('QP2', 'part_Q: the K of LC (same values)', "pool.map(q_job, [(K, j8) for K in Q_K for j8 in Q_TH]", "pool.map(q_job, [(K, j8) for K in LC_K for j8 in Q_TH]", 0),
 ('QM1', 'main: Q runs part_LC', "('LC', part_LC), ('Q', part_Q)", "('LC', part_LC), ('Q', part_LC)", 0),
 # q_lines, report
 ('QL1', 'q_lines: the smallest deviation summarised', "by.setdefault(x['j8'], []).append(max(x['dev']))\n        out.append(f\"| {x['K']} | {x['j8']}π/8 | \" + ' | '.join(f\"{d:.2e}\"", "by.setdefault(x['j8'], []).append(min(x['dev']))\n        out.append(f\"| {x['K']} | {x['j8']}π/8 | \" + ' | '.join(f\"{d:.2e}\"", 0),
 ('QL2', 'q_lines: ratios inverted', "f\"{v[i + 1] / v[i]:.3f}\" for i in range(len(v) - 1)) + ')'\n        for j, v in sorted(by.items())) + '.']\n    return out\n\n\ndef alias_rate", "f\"{v[i] / v[i + 1]:.3f}\" for i in range(len(v) - 1)) + ')'\n        for j, v in sorted(by.items())) + '.']\n    return out\n\n\ndef alias_rate", 0),
 ('QL3', 'q_lines: sorted by angle reversed', "for j, v in sorted(by.items())) + '.']\n    return out\n\n\ndef alias_rate", "for j, v in sorted(by.items(), reverse=True)) + '.']\n    return out\n\n\ndef alias_rate", 0),
 ('QR1', 'report: Q not reported', "+ lc_lines(LC) + q_lines(Q) + d_lines(D, S, R)", "+ lc_lines(LC) + d_lines(D, S, R)", 0),
 ('QR2', 'report: Q after R', "+ lc_lines(LC) + q_lines(Q) + d_lines(D, S, R) + s_lines(S) + r_lines(R)", "+ lc_lines(LC) + d_lines(D, S, R) + s_lines(S) + r_lines(R) + q_lines(Q)", 0),
 # round-90 code the round's specs leave out
 ('FA8', 'frame_series_any: stop after 2 terms', "if len(terms) > 10 and 0 < max(terms[-2:]) < tol * S:", "if len(terms) > 2 and 0 < max(terms[-2:]) < tol * S:", 0),
 ('FA9', 'frame_series_any: absolute low defect', "low = max(low, abs(D) / (abs(E_) + abs(Q) + 1))", "low = max(low, abs(D))", 0),
 ('FA10', 'frame_series_any: zero frame as He_n', "h_prev, h_cur = h_cur, [x_ * hc for x_, hc in zip(xs, h_cur)]", "h_prev, h_cur = h_cur, [x_ * hc + n * hp for x_, hc, hp in zip(xs, h_cur, h_prev)]", 0),
 ('FA11', 'frame_series_any: tol 1e-20', "def frame_series_any(K, rho, a, prec, tol_exp=40, zero=None):", "def frame_series_any(K, rho, a, prec, tol_exp=20, zero=None):", 0),
 ('SJ1', 's_job: check difference at the lower precision', "ctx.prec = rows[-1]['prec'] + S_CHECK\n        d = abs(S2 / S_ - 1)", "ctx.prec = rows[-1]['prec']\n        d = abs(S2 / S_ - 1)", 0),
 ('LJ2', 'lc_job: precision with cos', "y2 = (max(L_S) * math.sqrt(K) * math.sin(j8 * math.pi / 8)) ** 2", "y2 = (max(L_S) * math.sqrt(K) * math.cos(j8 * math.pi / 8)) ** 2", 0),
 ('LJ3', 'lc_job: the arg of the ratio dropped (abs)', "arg=[float(np.angle(v)) for v in q]", "arg=[float(abs(np.angle(v))) for v in q]", 0),
 ('RL4', 'r_lines: the none-list sign format', "f\"τ/σ = {x['r']:g}, {x['name']} ({x['rho']:+.3f})\" for x in none", "f\"τ/σ = {x['r']:g}, {x['name']} ({abs(x['rho']):.3f})\" for x in none", 0),
 ('RR2', 'r_rules: convergent series by b < 1 + a', "if a is not None and 1 - b > 0:", "if a is not None and 1 - b > -a:", 0),
 ('RF4', 'rule_frame: C with (1 - rho)^{1/2}', "return rho, a, a - rho, np.sqrt(1 + tau2) * np.sqrt(1 - rho ** 2) * v / tau2", "return rho, a, a - rho, np.sqrt(1 + tau2) * np.sqrt(1 - rho) * v / tau2", 0),
]
