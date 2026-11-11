"""Single-point mutants of experiments/exp_frameless.py for referee round 132's changes to Proposition D‴.62's experiment (ρ = −1 in R
and S, moment matching at Kuo–Woźniakowski's threshold 2αγ² = 1 and beyond (S_NEW, F5), the pinned small product rules R_PIN, the
sweep's rules `t_sweep_rules`, assumption (F) on the axes (G)), and the referee's survivors re-run against the new tests. Run with

    venv/bin/python tools/mutants.py --exp exp_frameless --tests tests/test_frameless.py --raw results/chaos/frameless_raw.npz \\
        --md results/chaos/frameless.md --out-env FL_OUT --spec tools/mutants/round132_D62_own.py --work <scratch>/mut --slot a \\
        --extra results/chaos/width_mismatch_raw.npz

(id, what, old, new, k): `old` replaced by `new` at its k-th occurrence. B00 is the control (no change). E20 … E40 are referee round
132's mutants (its mutants_ref132.py), whose anchors still apply.

Run (42 mutants and the control, 10 min): all 33 B mutants killed (by the tests; 16 of them also change the md); the referee's E20, E21,
E32, E33 and E40 killed by the new tests (the two-dimensional R rows and the pins recomputed, a ρ > 1 rule of T, the sweep's rules
and one sweep grid recomputed); E30, E35, E36 and E38 survive, equivalent as below.

Equivalent on the data (recorded, not killed): E30 (`f_solve` without taking the data's midpoints exact: the printed midpoints of the
fits are unchanged at F_PREC bits), E35 (`t_margins`' supremum without sector b's R1 corner: R1 < R2 in sector b at every stored rule,
so the maximum is R2's either way), E38 (the middle node of an odd K not zeroed: Newton from LAPACK's 0 stays at 0) and E36 (`t_grid`'s
tail_b with y₁y₂ < 0 instead of ≤ 0: on y₁y₂ = 0 with max|y_i| ≥ 2 one coordinate is 0, where L⁺ = −∞, so the maximum is unchanged)."""
M = [
    ('B00', 'control', '', '', 0),
    ('B01', "A_: 'a√b' as a·b", "return (arb(a) if a else arb(1)) * arb(b).sqrt()", "return (arb(a) if a else arb(1)) * arb(b)", 0),
    ('B02', 'tau2_frac: a√b → ab', "return Fraction(a or 1) ** 2 * Fraction(b)", "return Fraction(a or 1) * Fraction(b)", 0),
    ('B03', "rule_v: 'rho=-1' as τ²/(2 + τ²)", "        return t2 / (1 + t2)\n", "        return t2 / (2 + t2)\n", 0),
    ('B04', 'rm1_C2: (1 − λ) for (1 − λ)²', "return (1 - lam) ** 2 / (2 + tau2)", "return (1 - lam) / (2 + tau2)", 0),
    ('B05', 'rm1_s1: s = −λ', "return s1_formal(lam, Fraction(0))", "return s1_formal(-lam, Fraction(0))", 0),
    ('B06', 'r_rm1: c = τ²/σ*²', "    c = t2 / (2 * ss2)\n", "    c = t2 / ss2\n", 0),
    ('B07', "r_rm1: the density's normalisation (πc)^{1/2}", "/ (2 * arb.pi() * c).sqrt() + arb(0, tail)", "/ (arb.pi() * c).sqrt() + arb(0, tail)", 0),
    ('B08', 'r_rm1: without the factor 1/σ*', "arb(0, tail)) / ss2.sqrt()", "arb(0, tail))", 0),
    ('B09', "r_rm1: the tail bound e^{−2βL²}", "(-(beta * L * L)).exp()).upper()", "(-(2 * beta * L * L)).exp()).upper()", 0),
    ('B10', 'r_rm1: L without the + 5', "L = (R_LCUT / beta).sqrt() + 5", "L = (R_LCUT / beta).sqrt()", 0),
    ('B11', "r_rm1: φ_K(x) with e^{x²/4}", "p = (x * x / 2).exp() - sum((wk * (x * z).exp()", "p = (x * x / 4).exp() - sum((wk * (x * z).exp()", 0),
    ('B12', 'R_RM1: K = 5', "R_RM1 = (('1', 6),)", "R_RM1 = (('1', 5),)", 0),
    ('B13', 'part_R: the pins at 50 digits', "pin=[dict(r_2d(tau, x, K, Ms, 60)", "pin=[dict(r_2d(tau, x, K, Ms, 50)", 0),
    ('B14', 'r_2d: Ms ignored', "for M in (Ms or R_M):", "for M in R_M:", 0),
    ('B15', "t_sweep_rules: v/v_hs perturbed by 10⁻¹⁵", "out.append((repr(float(tau)), repr(float(w))))", "out.append((repr(float(tau)), repr(float(w) * (1 + 1e-15))))", 0),
    ('B16', "g_D: the real axis's +½log(1 − η²/4)", "- (1 - e * e / 4).log() / 2)", "+ (1 - e * e / 4).log() / 2)", 0),
    ('B17', "g_D: the imaginary axis's ½log(1 − y²/4)", "- (1 + e * e / 4).log() / 2)", "- (1 - e * e / 4).log() / 2)", 0),
    ('B18', 'g_D: log(K!/(2K)!)', "lf = arb(math.factorial(2 * K)).log() - arb(math.factorial(K)).log()", "lf = arb(math.factorial(K)).log() - arb(math.factorial(2 * K)).log()", 0),
    ('B19', 'g_D: the real axis −Kη²/2', "phi_r(xi).log() - K * e * e / 4", "phi_r(xi).log() - K * e * e / 2", 0),
    ('B20', 'part_G: 1.5K bits', "bits = int(2.5 * K + G_MARGIN)", "bits = int(1.5 * K + G_MARGIN)", 0),
    ('B21', 'G_DIGITS 40', "G_DIGITS = 50 ", "G_DIGITS = 40 ", 0),
    ('B22', 'g_lines: the last three points', "g_extrap(Ks[-4:], D[-4:])", "g_extrap(Ks[-3:], D[-3:])", 0),
    ('B23', 'g_extrap: the powers from 1/K', "M = arb_mat([[arb(1) / arb(K) ** j for j in range(n)] for K in Ks])\n    return M.solve",
     "M = arb_mat([[arb(1) / arb(K) ** (j + 1) for j in range(n)] for K in Ks])\n    return M.solve", 0),
    ('B24', 'f_new: y_K without the − 1', "e[K] ** 2 / (out['C'] * out['t'] ** (2 * K)) - 1) for K in ks]", "e[K] ** 2 / (out['C'] * out['t'] ** (2 * K))) for K in ks]", 0),
    ('B25', 'f_new: the cheap bound (τ²/2)^K', "(P['t2'] / 2) ** (2 * K)", "(P['t2'] / 2) ** K", 0),
    ('B26', 'f_new: (4K − 1)!! with 2^K', "(arb(2) ** (2 * K) * arb(math.factorial(2 * K)))", "(arb(2) ** K * arb(math.factorial(2 * K)))", 0),
    ('B27', "s_is_a6: S_RULES for S_A6", "in {(t, w) for t, w, _ in S_A6}", "in {(t, w) for t, w, _ in S_RULES}", 0),
    ('B28', 'part_F: F5 without the last new rule', "new=[f_new(r) for r in S if not s_is_a6(r)]", "new=[f_new(r) for r in S[:-1] if not s_is_a6(r)]", 0),
    ('B29', "s_lines: S_NEW's time from A6's", "{sum(r['secs'] for r in SN):.0f} s.", "{sum(r['secs'] for r in S6):.0f} s.", 0),
    ('B30', 'new_lines: 2αγ² = τ²', "| {t2 / 2} |", "| {t2} |", 0),
    ('B31', 'S_NEW: ρ = −1 to K = 280', "S_NEW = (('1', 'rho=-1', 300)", "S_NEW = (('1', 'rho=-1', 280)", 0),
    ('B32', "rule_name: ρ = −1 named as a multiple of v_hs", "else 'v = τ²/σ*² (ρ = −1)' if w == 'rho=-1' else", "else", 0),
    ('B33', "params: ρ = −1's tolerance 10⁻¹⁶", "elif abs(rf + 1) <= 1e-12:", "elif abs(rf + 1) <= 1e-16:", 0),
    # referee round 132's survivors and equivalents, re-run
    ('E20', 'r_2d: b perturbed by 2^{−100} (relative)', "                b = kb * yj\n", "                b = kb * yj * (1 + arb(2) ** -100)\n", 0),
    ('E21', 'r_2d: dM at the midpoint', "dM=float(abs((vals[-1] - vals[0]) / ex).upper())", "dM=float(abs((vals[-1] - vals[0]) / ex).mid())", 0),
    ('E30', 'F: midpoints not taken exact', "b = arb_mat([[arb(y.mid())] for y in ys])", "b = arb_mat([[y] for y in ys])", 0),
    ('E32', 'T: t_grid on the imaginary branch, L⁺ + 10⁻³', "return out + (y * y / 2 if real else 0)", "return out + (y * y / 2 if real else 1e-3)", 0),
    ('E33', 'T: the sweep upper branch from 1.01σ*', "list(np.geomspace(1.001 * ss, 50, T_WN))", "list(np.geomspace(1.01 * ss, 50, T_WN))", 0),
    ('E35', 'T: t_margins sup without the R1 corner of sector b', "sup = max(max(red['a']), max(red['b']))", "sup = max(max(red['a']), red['b'][1])", 0),
    ('E36', "T: t_grid's tail_b sector y₁y₂ < 0", "tail_b=float(F[tail & (Y1 * Y2 <= 0)].max())", "tail_b=float(F[tail & (Y1 * Y2 < 0)].max())", 0),
    ('E38', 'gh_nodes: the middle node not zeroed (odd K)', "        zs[-1] = arb(0)\n", "        pass\n", 0),
    ('E40', 'T: the sweep v/v_hs lower end 0.03', "list(np.geomspace(0.02, 0.999 / ss, T_WN))", "list(np.geomspace(0.03, 0.999 / ss, T_WN))", 0),
]
