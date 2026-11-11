"""Single-point mutants of experiments/exp_frameless.py (Proposition D‴.62's experiment), run with

    venv/bin/python tools/mutants.py --exp exp_frameless --tests tests/test_frameless.py --raw results/chaos/frameless_raw.npz \\
        --md results/chaos/frameless.md --out-env FL_OUT --spec tools/mutants/d62_initial_own.py --work <scratch>/mut --slot a \\
        --extra results/chaos/width_mismatch_raw.npz

(id, what, old, new, k): `old` replaced by `new` at its k-th occurrence. D00 is the control (no change).

First run (62 mutants and the control, 8 min): 60 killed (by the tests; 36 of them also change the md), two survivors. D33 (the recheck
with fewer bits) is killed since by test_s_rows_recomputed_and_the_precision's check of the stored recheck's precision. D40 (r_2d with
the two variances swapped) is equivalent: φ_K is even (Z and the rule are symmetric), so E[φ_K(a′ + b′)φ_K(a′ − b′)] with the variances
exchanged is E[φ_K(b + a)φ_K(b − a)] = E[φ_K(a + b)φ_K(a − b)], and the product rule uses the same nodes on both axes."""
M = [
    ('D00', 'control', '', '', 0),
    # the closed forms
    ('D01', 'C(v): 1 − s_d → 1 + s_d', "ss * v * (1 - P['s_d'])", "ss * v * (1 + P['s_d'])", 0),
    ('D02', 'C(v): 1 + ϖ → 1 − ϖ', "(t2 + v * ss ** 2) * (1 + P['varpi'])).sqrt()", "(t2 + v * ss ** 2) * (1 - P['varpi'])).sqrt()", 0),
    ('D03', 'C(v): √2 → 2^{1/4}', "return 2 ** A_('0.5') * ss * v * (1 - P['s_d'])", "return 2 ** A_('0.25') * ss * v * (1 - P['s_d'])", 0),
    ('D04', 'C_o: 1 + 1/ϖ → 1 + ϖ', "(1 + 1 / P['varpi'])).sqrt()", "(1 + P['varpi'])).sqrt()", 0),
    ('D05', 'C_o: also at ρ = −1', "if P['kind'] not in ('rho>1', 'rho<-1'):", "if P['kind'] not in ('rho>1', 'rho<-1', 'MM'):", 0),
    ('D06', 's₁: 1/16 → 1/8', "return one / 16 + 2 * s * (1 + s)", "return one / 8 + 2 * s * (1 + s)", 0),
    ('D07', 's₁: 2s(1 + s) → 2s(1 − s)', "one / 16 + 2 * s * (1 + s) / (1 - s) ** 3", "one / 16 + 2 * s * (1 - s) / (1 - s) ** 3", 0),
    ('D08', 's₁: the ϖ term (1 + s)² → (1 − s)²', "- w * (1 + s) ** 2 / (8 * (1 + w)", "- w * (1 - s) ** 2 / (8 * (1 + w)", 0),
    ('D09', 's₁: s/(4(1 − s)²) → s/(4(1 − s))', "+ s / (4 * (1 - s) ** 2) - w", "+ s / (4 * (1 - s)) - w", 0),
    ('D10', 'params: s_b = c_b/(1 + c_b)', "sa, sb = c / (1 - c), cb / (1 - cb)", "sa, sb = c / (1 - c), cb / (1 + cb)", 0),
    ('D11', 'params: c_b = c + v/2', "    cb = c - v / 2\n    sa, sb", "    cb = c + v / 2\n    sa, sb", 0),
    ('D12', 'params: the dominant by the smaller modulus', "if float(abs(sa).mid()) >= float(abs(sb).mid()):", "if float(abs(sa).mid()) <= float(abs(sb).mid()):", 0),
    ('D13', 'params: ρ = 1 + 2/τ² + 2/v', "rho = 1 + 2 / t2 - 2 / v\n    if float", "rho = 1 + 2 / t2 + 2 / v\n    if float", 0),
    ('D14', 'params: the cusp test loose', "kind = 'cusp' if abs(float((1 + w).mid())) < 1e-12", "kind = 'cusp' if abs(float((1 + w).mid())) < 1e-1", 0),
    ('D15', 'params: frames include ρ = 1', "if abs(rf) < 1 - 1e-12:", "if abs(rf) <= 1 + 1e-12:", 0),
    ('D16', 'B₀: 1 − q → 1 + q', "return (Co / C_const(r, v, P) * (1 - q) / (2 * P['t'] * q),", "return (Co / C_const(r, v, P) * (1 + q) / (2 * P['t'] * q),", 0),
    ('D17', 'B₁/B₀: s₁′ with ϖ instead of 1/ϖ', "s1_formal(P['s_o'], 1 / P['varpi']) - s1_formal", "s1_formal(P['s_o'], P['varpi']) - s1_formal", 0),
    ('D18', 'mm_s1: 2λ(1 − λ) → 2λ(1 + λ)', "- 2 * lam * (1 - lam) / (1 + lam) ** 3", "- 2 * lam * (1 + lam) / (1 + lam) ** 3", 0),
    ('D19', 'mm_C2: (1 + λ)³/4', "return (1 + lam) ** 3 / 2", "return (1 + lam) ** 3 / 4", 0),
    ('D20', 'D‴.61 constant: 1 − r → 1 + r on the real axis', "rr, S = 1 / (2 * gx), (1 - 1 / (2 * gx))", "rr, S = 1 / (2 * gx), (1 + 1 / (2 * gx))", 0),
    ('D21', 'D‴.61 constant: γy = B + ½', "v * ss2 / (v * ss2 - t2) - A_('0.5')", "v * ss2 / (v * ss2 - t2) + A_('0.5')", 0),
    ('D22', "rule_v: the Mercer frame's 2τ²/(b − 1)", "return 2 * t2 / (b + 1)", "return 2 * t2 / (b - 1)", 0),
    ('D23', 'rule_v: x·τ²/σ*² for x·v_hs', "return A_(w) * t2 / (1 + t2).sqrt()", "return A_(w) * t2 / (1 + t2)", 0),
    ('D24', 's1_r129: the √q term with a minus', "+ sq * (1 - rho_) ** 2 / (8 * (1 - sq)", "- sq * (1 - rho_) ** 2 / (8 * (1 - sq)", 0),
    ('D25', 's1_unsimplified: + ϖω/8', "return (one / 16 - w * om / 8 + s * f1 / 2", "return (one / 16 + w * om / 8 + s * f1 / 2", 0),
    ('D26', 'F₁: y/(4(1 − y)²)', "y / (8 * (1 - y) ** 2) + y / (8 * (1 - y))", "y / (4 * (1 - y) ** 2) + y / (8 * (1 - y))", 0),
    ('D27', 'phi_rate: b²/|s_b| without the ½', "- b * b / (2 * abs(sb))", "- b * b / abs(sb)", 0),
    # S: the exact errors
    ('D28', 'gh_nodes: Newton in the wrong direction', "x = arb((x - dx).mid())", "x = arb((x + dx).mid())", 0),
    ('D29', "gh_nodes: weights K!/(K He_{K−1})", "ws.append(fK / (K * hm) ** 2)", "ws.append(fK / (K * hm))", 0),
    ('D30', 'err2: the Gram sum without the factor 2', "gram += 2 * a[j] * s", "gram += a[j] * s", 0),
    ('D31', 'err2: κ = e^{−(y − y′)²/2}', "s += a[k] * (-(d * d) / 4).exp()", "s += a[k] * (-(d * d) / 2).exp()", 0),
    ('D32', 'err2: the reweighting e^{cz²} with c of the wrong sign', "g = (1 - v / tau2) / 2\n    a = [", "g = (v / tau2 - 1) / 2\n    a = [", 0),
    ('D33', 's_bits: the recheck with fewer bits', "+ S_MARGIN + extra)", "+ S_MARGIN - extra)", 0),
    ('D34', 'e_row: one digit fewer', "return (K, e.str(S_DIGITS, radius=False)", "return (K, e.str(S_DIGITS - 1, radius=False)", 0),
    ('D35', 't_float: min instead of max', "return max(abs(float(h_(t2 / v).mid())), abs(float(h_(v * (1 + t2) / t2).mid())))",
     "return min(abs(float(h_(t2 / v).mid())), abs(float(h_(v * (1 + t2) / t2).mid())))", 0),
    ('D36', 'part_S: resume key off by one', "r = have.get((tau, w, kmax))", "r = have.get((tau, w, kmax + 1))", 0),
    # R: the representation
    ('D37', 'r_mm: B ~ N(0, τ²) instead of N(0, τ²/2)', "return p * p * (-(W * W) / t2).exp()", "return p * p * (-(W * W) / (2 * t2)).exp()", 0),
    ('D38', 'r_mm: the truncation bound 2e^{−L²/τ²}', "tail = (4 * (-(L * L) / t2).exp()).upper()", "tail = (2 * (-(L * L) / t2).exp()).upper()", 0),
    ('D39', 'r_2d: φ(a + b)φ(a + b)', "tot += wi * wj * phi(a + b) * phi(a - b)", "tot += wi * wj * phi(a + b) * phi(a + b)", 0),
    ('D40', 'r_2d: the variances swapped', "ka, kb = abs(P['c']).sqrt(), abs(P['c_b']).sqrt()", "ka, kb = abs(P['c_b']).sqrt(), abs(P['c']).sqrt()", 0),
    ('D41', 'phi_r: e^{x²/4}', "return (x * x / 2).exp() - sum", "return (x * x / 4).exp() - sum", 0),
    ('D42', 'phi_i: cos → sin', "return (-(w * w) / 2).exp() - sum((wk * (w * z).cos()", "return (-(w * w) / 2).exp() - sum((wk * (w * z).sin()", 0),
    # A
    ('D43', 'part_A: e0 from e1', "e0=str(y['e0'])", "e0=str(y['e1'])", 0),
    ('D44', "a6_w: 'moment matching' unmapped", "if name == 'moment matching':\n        return 'mm'", "if name == 'moment-matching':\n        return 'mm'", 0),
    # T
    ('D45', 't_reduced: k = |s|/(1 + |s|)', "ka, kb = abs(sa) / (1 - abs(sa)), abs(sb) / (1 - abs(sb))\n    y = np",
     "ka, kb = abs(sa) / (1 + abs(sa)), abs(sb) / (1 - abs(sb))\n    y = np", 0),
    ('D46', 't_reduced: Q(2, y₂) with k_a, k_b swapped', "e = L_rate(y2) - ((2 + y2) ** 2 / (8 * ka) + (2 - y2) ** 2 / (8 * kb))",
     "e = L_rate(y2) - ((2 + y2) ** 2 / (8 * kb) + (2 - y2) ** 2 / (8 * ka))", 0),
    ('D47', 't_grid: the real axis unsubtracted on the wrong branch', "real = c > 0", "real = c < 0", 0),
    ('D48', 't_grid: the tail from 2.5', "tail = np.maximum(np.abs(Y1), np.abs(Y2)) >= 2", "tail = np.maximum(np.abs(Y1), np.abs(Y2)) >= 2.5", 0),
    ('D49', 'L_rate: + y²/4', "return 1 + np.log(y * y / 4) - y * y / 4", "return 1 + np.log(y * y / 4) + y * y / 4", 0),
    ('D50', "t_margins: the o-sector against the d-sector's sup", "2 * math.log(abs(so)) - max(red[o]))", "2 * math.log(abs(so)) - max(red[d]))", 0),
    # F
    ('D51', 'F: the third window from 40', "kmin = max(60, kmax // j) if j == 4", "kmin = max(40, kmax // j) if j == 4", 0),
    ('D52', 'F: q^K terms where q^{K_min} is below the resolution', "m = F_M if (two and float(q.mid()) ** kmin > F_QRES)", "m = F_M if (two and float(q.mid()) ** kmin < F_QRES)", 0),
    ('D53', "F: g₁/g₀ − s₁ with the formal s₁'s sign", "b1=g1 / g0 - x[0]", "b1=g1 / g0 + x[0]", 0),
    ('D54', 'F: the verdict on the worse window', "best = min(float(r['g0'].mid()) for r in rid)", "best = max(float(r['g0'].mid()) for r in rid)", 0),
    ('D55', "F: the overlap's K − 1 paired with e_K", "for K, s in ((row['K'] - 1, row['e0']), (row['K'], row['e1'])):", "for K, s in ((row['K'] - 1, row['e1']), (row['K'], row['e0'])):", 0),
    ('D56', 'F: A6 at K = 160 with e_{K−1}', "R = arb(row['e1']) ** 2 / (C * P['t'] ** (2 * F_A6K))", "R = arb(row['e0']) ** 2 / (C * P['t'] ** (2 * F_A6K))", 0),
    ('D57', 'F: the data y_K = e_K/(C t^{2K})', "y = {K: e[K] ** 2 / (C * t ** (2 * K)) for K in e}", "y = {K: e[K] / (C * t ** (2 * K)) for K in e}", 0),
    ('D58', 'F: k₀ with 20 digits', "F_K0DIG = 30 ", "F_K0DIG = 20 ", 0),
    # the report and the raw
    ('D59', 'report: B₀ to 9 digits', "fs(r['B0'], 10)", "fs(r['B0'], 9)", 0),
    ('D60', 'fs: the exponent as arb prints it', "s = f\"{m}e{int(e):+03d}\"", "s = f\"{m}e{int(e)}\"", 0),
    ('D61', 'load: no list conversion', "return {k: (list(v) if isinstance(v, np.ndarray) else v) for k, v in out.items()}", "return out", 0),
    ('D62', "main: --only A runs T's part", "res[k] = {'R': part_R, 'A': part_A, 'T': part_T}[k]()", "res[k] = {'R': part_R, 'A': part_T, 'T': part_T}[k]()", 0),
]
