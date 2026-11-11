# (id, what, old, new, k): single-point mutants of round 90's new or changed code in exp_gh_peano.py (k: the k-th occurrence)
M = [
 ('X00', 'control', 'R_PREC = 80 ', 'R_PREC = 80 ', 0),
 # gammas_ab, the (a, b) forms
 ('GA1', 'gammas_ab: A from 1 + b', "return (1 - b) / a + 0.5, (1 + b) / a - 0.5", "return (1 + b) / a + 0.5, (1 + b) / a - 0.5", 0),
 ('GA2', 'gammas: b = ρv′', "return gammas_ab(rho * vp, rho * (vp - 1))", "return gammas_ab(rho * vp, rho * vp)", 0),
 ('GA3', 'log_s_pred_ab: prefactor π·b', "return -log(pi * a) + 2 * lgamma(K + 1)", "return -log(pi * b) + 2 * lgamma(K + 1)", 0),
 ('GA4', 'log_s_simple_ab: cusp 2r/b', "np.log(hyp_half(r * r)) + np.log(2 * r / a)", "np.log(hyp_half(r * r)) + np.log(2 * r / b)", 0),
 ('GA5', 'log_s_simple_ab: off-cusp a²', "- np.log(a * np.sqrt(2 * g * (g2 - g)))", "- np.log(a * a * np.sqrt(2 * g * (g2 - g)))", 0),
 # frame_series_any
 ('FA1', 'frame_series_any: zero threshold 1e-10', "zero = abs(rho) < Decimal(10) ** -30 if zero is None else zero", "zero = abs(rho) < Decimal(10) ** -10 if zero is None else zero", 0),
 ('FA2', 'frame_series_any: ρ < 0 with He_n (no sign flip in the recurrence)', "[x_ * hc + n * hp for x_, hc, hp in zip(xs, h_cur, h_prev)]", "[x_ * hc - n * hp for x_, hc, hp in zip(xs, h_cur, h_prev)]", 0),
 ('FA3', 'frame_series_any: E H̃_n with (u − 1)', "E_ = dfact if zero else (u + 1) ** (n // 2) * dfact", "E_ = dfact if zero else (u - 1) ** (n // 2) * dfact", 0),
 ('FA4', 'frame_series_any: ρ < 0 weights ρ (signed)', "w_ = a if zero else -rho", "w_ = a if zero else rho", 0),
 ('FA5', 'frame_series_any: ρ = 0 weights 1', "w_ = a if zero else -rho", "w_ = 1 if zero else -rho", 0),
 ('FA6', 'frame_series_any: u = a/ρ', "u = Decimal(0) if zero else -a / rho", "u = Decimal(0) if zero else a / rho", 0),
 ('FA7', 'frame_series_any: ρ > 0 at v′ = a', "return WM.frame_series(K, rho, a / rho, prec=prec, tol_exp=tol_exp)", "return WM.frame_series(K, rho, a, prec=prec, tol_exp=tol_exp)", 0),
 # rule_frame, r_rules, r_job
 ('RF1', 'rule_frame: a = v(1 − ρ)/2', "a = v * (1 - rho ** 2) / 2", "a = v * (1 - rho) / 2", 0),
 ('RF2', 'rule_frame: b = a + ρ', "return rho, a, a - rho, np.sqrt(1 + tau2)", "return rho, a, a + rho, np.sqrt(1 + tau2)", 0),
 ('RF3', 'rule_frame: |ρ| ≤ 1 admitted', "if not -1 < rho < 1:\n        return rho, None, None, None", "if not -1 <= rho <= 1:\n        return rho, None, None, None", 0),
 ('RR1', 'r_rules: A > 0 not required', "if a is not None and 1 - b > 0:", "if a is not None:", 0),
 ('RJ1', 'r_job: the constant without (1 − ρ²)^{1/2}', "Cd = (1 + tau2).sqrt() * (1 - rh * rh).sqrt() * v / tau2", "Cd = (1 + tau2).sqrt() * v / tau2", 0),
 ('RJ2', 'r_job: a in the check from v/2', "ad = v * (1 - rh * rh) / 2", "ad = v / 2", 0),
 ('RJ3', 'r_job: closed from the per-angle prediction', "closed=float(np.exp(le(y) - log_s_simple_ab(y['K'], a, b)))", "closed=float(np.exp(le(y) - log_s_pred_ab(y['K'], a, b)))", 0),
 ('RJ4', 'r_job: e_K, not e_K²', "le = lambda y: 2 * float(Decimal(y['e1']).ln()) - np.log(C)", "le = lambda y: float(Decimal(y['e1']).ln()) - np.log(C)", 0),
 # gen_ratio_c, lc_limit, lc_job
 ('GC1', 'gen_ratio_c: e^{−ξ²/2} phase +xy', "er, ei = m * WM.dcos(x * y, pi), -m * WM.dcos(x * y - pi / 2, pi)", "er, ei = m * WM.dcos(x * y, pi), m * WM.dcos(x * y - pi / 2, pi)", 0),
 ('GC2', 'gen_ratio_c: the denominator phase without −xy/2', "ph = 2 * K * th - x * y / 2", "ph = 2 * K * th", 0),
 ('GC3', 'gen_ratio_c: θ = jπ/4', "th = pi * j8 / 8", "th = pi * j8 / 4", 0),
 ('GC4', 'gen_ratio_c: the division by the conjugate dropped', "return float((nr * dr + ni * di) / den), float((ni * dr - nr * di) / den)", "return float((nr * dr - ni * di) / den), float((ni * dr + nr * di) / den)", 0),
 ('LL1', 'lc_limit: the conjugate root', "return cmath.sqrt(1 - s * s / 4 * cmath.exp(2j * j8 * np.pi / 8))", "return cmath.sqrt(1 - s * s / 4 * cmath.exp(-2j * j8 * np.pi / 8))", 0),
 ('LJ1', 'lc_job: precision without e^{y²/2}', "prec = 3 * K + 60 + int(math.ceil(y2 / (2 * math.log(10))))", "prec = 3 * K + 60", 0),
 # the report: s_lines' digits and summaries, r_lines' census, lc_lines
 ('SL1', 's_lines: the digits from the first K', "chk = float(-d.log10()) if d != 0 else float(rows[-1]['prec'])", "chk = float(-d.log10()) if d != 0 else float(rows[0]['prec'])", 0),
 ('SL2', 's_lines: the closed-form range by min only', "over the closed form {min(c for _, _, c in v):.5f}–{max(c for _, _, c in v):.5f}", "over the closed form {min(c for _, _, c in v):.5f}–{min(c for _, _, c in v):.5f}", 0),
 ('RL1', "r_lines: ρ = 0 counted as ρ > 0", "sgn = lambda x: 'ρ > 0' if x['rho'] > 1e-12 else", "sgn = lambda x: 'ρ > 0' if x['rho'] >= 0 else", 0),
 ('RL2', 'r_lines: within R_TOL by <', "near = [v for _, v in last if abs(v - 1) <= R_TOL]", "near = [v for _, v in last if abs(v - 1) < R_TOL / 2]", 0),
 ('RL3', 'r_lines: the closed summary from the ratios', "for key, lab in (('ratio', 'the per-angle prediction'), ('closed', 'the closed forms')):", "for key, lab in (('ratio', 'the per-angle prediction'), ('ratio', 'the closed forms')):", 0),
 ('LC1', 'lc_lines: the largest |arg| over the first job only', "the largest |arg| {max(abs(a) for x in LC for a in x['arg']):.1e}", "the largest |arg| {max(abs(a) for a in LC[0]['arg']):.1e}", 0),
 ('RP1', 'report: LC not reported', "body = p_lines(P) + m_lines(M) + l_lines(L) + lc_lines(LC) + s_lines(S) + r_lines(R)", "body = p_lines(P) + m_lines(M) + l_lines(L) + s_lines(S) + r_lines(R)", 0),
 ('KY1', 'KEYS without LC', "KEYS = ('P', 'M', 'L', 'S', 'R', 'LC')", "KEYS = ('P', 'M', 'L', 'S', 'R')", 0),
]
