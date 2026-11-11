# (id, what, old, new, k): single-point mutants of experiments/theory3a_tail.py (round theory3 A): the h = 1/kappa unit, the nonnegative squares, the ball series of atan(y)/y, the layer constants over the
# ball [0, h_hi], the envelope-based dormant terms (mass, peak, drift with the stage-resolved displacement integral and the window gating), the conditions, the candidate iteration.
# Run: tools/mutants.py --exp theory3a_tail --tests tests/theory3a_test.py --raw results/synthetic/theory3a_raw.npz --md results/synthetic/theory3a.md --out-env TH3A_OUT
#      --extra results/synthetic/theory_t1b_raw.npz --spec tools/mutants/theory3a_tail_own.py --work <scratch dir>
# Result (51 mutants and the control; the criterion is the tests file alone, 16 tests, about 15 s per mutant): 45 of the 51 mutants killed.  The 6 survivors are equivalent on the stored data (and the control survives, as it must):
#   V00 control (unmutated);  V07 (h_hi rounded down instead of up: the two differ by 1e-30 relative, the ball [0, h_hi] reaches every eps of the grid either way);  V28 (the part x < 1e-9 of the pre-layer
#   drift integral: contributes 1e-6 of 11, below every tolerance of the tests);  V29 (no cap 2r on the pre-layer displacement: the oscillation cap is the active one wherever 2r would be, pre_drift
#   unchanged to 1e-6);  V41 (H of a dormant kernel without the 0.43 a_min omega_pk term: omega_pk ~ 1e-9);  V44 (the dc condition on the lower end of dcmax: dcmax is a point ball, lower end = upper end);
#   V49 (theta condition with 2 theta: the dormant peaks are ~1e-9 against theta ~ 7e-3).  V21, V20, V35 and V10, V33, V34, V43, V46, V51 were survivors of the first test file and are killed by tests added for them
#   (kernel margin with a last-stage rate that reduces |u|, golden values of drift/dcmax/mass/pre-drift, window order with a gap between the correct and the shortened requirement, strictness factor, atan point ball, solve_D iteration).
# (the criterion is the tests file alone; `--from-raw` of theory3a_tail.py does not exist, its exit code is not part of the criterion; the tests read theory_t1b_raw.npz through TH3A_Q's default)
M = [
 ('V00', 'control (unmutated)', "def sqn(z):", "def sqn(z):", 0),
 ('V01', 'squares not clamped to >= 0', "return (z * z).nonnegative_part()", "return z * z", 0),
 ('V02', 'NaN guard of bmax_up removed', "    if a.is_nan() or b.is_nan():\n        raise ValueError('nan in bmax_up')\n", "", 0),
 ('V03', 'bmax_up returns the smaller upper end', "return a if float(a.upper()) >= float(b.upper()) else b", "return a if float(a.upper()) <= float(b.upper()) else b", 0),
 ('V04', 'sigma_0 = 19', "SIG0 = arb(20)", "SIG0 = arb(19)", 0),
 ('V05', 'envelope constant C = 25', "C_ENV = arb(100) / 3", "C_ENV = arb(100) / 4", 0),
 ('V06', 'envelope used from s >= 10', "S_MIN = arb(25)", "S_MIN = arb(10)", 0),
 ('V07', 'h_hi rounded down', "return arb((e / (2 * VS - e)).sqrt().upper())", "return arb((e / (2 * VS - e)).sqrt().lower())", 0),
 ('V08', 'eps(h) without the 1 + h^2', "return 2 * VS * h2 / (1 + h2)\n\n\ndef g_atan", "return 2 * VS * h2\n\n\ndef g_atan", 0),
 ('V09', 'atan series without alternating signs', "s_ = s_ + (pw / (2 * n + 1) if n % 2 == 0 else -pw / (2 * n + 1))", "s_ = s_ + pw / (2 * n + 1)", 0),
 ('V10', 'atan series without its error bar', "return s_ + (-err).union(err)", "return s_", 0),
 ('V11', 'atan series used up to y < 2', "if yu < 0.5:", "if yu < 2.0:", 0),
 ('V12', 'atan(y)/y lower end taken at the lower end of y', "lo_ = b_.atan() / b_", "lo_ = a_.atan() / a_ if yl > 0 else arb(1)", 0),
 ('V13', 'L = 1/lambda with 2 instead of 4', "self.L = 4 * self.G * self.vs * self.h / (arb.pi() * self.amax)", "self.L = 2 * self.G * self.vs * self.h / (arb.pi() * self.amax)", 0),
 ('V14', 'q without the 1 + h^2', "return 2 * self.vs * (self.h2 + sqn(x.tanh())) / (1 + self.h2)", "return 2 * self.vs * (self.h2 + sqn(x.tanh()))", 0),
 ('V15', 'w without the eps term', "* (self.eps + 2 * self.vs * sqn(x.sinh()))", "* (2 * self.vs * sqn(x.sinh()))", 0),
 ('V16', 'y = h tanh x instead of h/tanh x', "y = self.h / t", "y = self.h * t", 0),
 ('V17', 'd eta/dx bound too small (factor 2)', "return self.amax / (self.G * 2 * self.vs * sqn(x.tanh()))", "return self.amax / (self.G * 4 * self.vs * sqn(x.tanh()))", 0),
 ('V18', 'M with 2 v_s instead of 4 v_s', "(q + 4 * self.vs * sqn(x.sinh()))", "(q + 2 * self.vs * sqn(x.sinh()))", 0),
 ('V19', 'Psi without the factor 2', "return self.amax ** 2 * sqn(x.sinh()) / (2 * self.G * self.m)", "return self.amax ** 2 * sqn(x.sinh()) / (self.G * self.m)", 0),
 ('V20', 'kernel margin without the last breakpoint', "vals = [abs(float(Lm['U_br'][i][j])) for i in range(K + 1)]", "vals = [abs(float(Lm['U_br'][i][j])) for i in range(K)]", 0),
 ('V21', 'kernel margin: horizon end without the rate', "uend = float(Lm['U_br'][K][j]) + float(Lm['P'][K][j]) * (tau_end - float(taus[K]))", "uend = float(Lm['U_br'][K][j])", 0),
 ('V22', 'T1 sum with the wrong weight', "T1 = T1 + up(f) * wd_", "T1 = T1 + up(f) * wd_ * arb(9) / 10", 0),
 ('V23', 'Atr without the w factor', "Atr = Atr + up(U.w(xi) * f) * wd_", "Atr = Atr + up(f) * wd_", 0),
 ('V24', 'T2 with half the rate', "T2 = logpos_up(ratio) / kd", "T2 = logpos_up(ratio) / (2 * kd)", 0),
 ('V25', 'mV at the upper end of x_#', "mV = U.M(arb(float(lo(xs)))).lower()", "mV = U.M(arb(float(up(xs)))).upper()", 0),
 ('V26', 'pre with v_lo instead of g_-', "pre = U.Psi(x1) / gminus", "pre = U.Psi(x1) / vlo", 0),
 ('V27', 'kappa the larger rate', "kappa = k1 if float(up(k1)) < float(up(k2)) else k2", "kappa = k1 if float(up(k1)) > float(up(k2)) else k2", 0),
 ('V28', 'pre-layer drift: the part x < 1e-9 dropped', "acc = up(Hs * muc * He / gminus ** 2 * c_in) * arb(1e-9)", "acc = arb(0)", 0),
 ('V29', 'pre-layer drift: no cap by 2r', "dcx = dcx if float(up(dcx)) < float(up(2 * r)) else 2 * r", "dcx = dcx", 0),
 ('V30', 'envelope kappa_j without the cube', "kj[j] = up(cmin[j] ** 3 / den)", "kj[j] = up(cmin[j] / den)", 0),
 ('V31', 'distance test against 0 instead of 26 L', "if not float(lo(aj[j])) > float(up((S_MIN + 1) * Lb)):", "if not float(lo(aj[j])) > 0:", 0),
 ('V32', 'mass without kappa_j', "massL = {j: up(C_ENV * sqn(Lb) * kj[j] * Qj[j]) for j in non}", "massL = {j: up(C_ENV * sqn(Lb) * Qj[j]) for j in non}", 0),
 ('V33', 'peak with L^2', "pk = {j: up(Lb ** 3 * ej[j]) for j in non}", "pk = {j: up(Lb ** 2 * ej[j]) for j in non}", 0),
 ('V34', 'omega_pk x Tfin with the mass sum', "omega_T = up(sum((TE * sqn(Lb) * ej[j] for j in non), arb(0)) / amin)", "omega_T = up(sum((sqn(Lb) * ej[j] for j in non), arb(0)) / amin)", 0),
 ('V35', 'windows gate: all windows counted from the start of the piece', "hs = up(sum((hd[i] for i in range(K) if gate[i] <= q0), arb(0)) + arb(43) / 100 * mass_tot)", "hs = up(sum((hd[i] for i in range(K)), arb(0)) + arb(43) / 100 * mass_tot)", 0),
 ('V36', 'gate: windows counted never', "hs = up(sum((hd[i] for i in range(K) if gate[i] <= q0), arb(0)) + arb(43) / 100 * mass_tot)", "hs = up(sum((hd[i] for i in range(K) if gate[i] > q0), arb(0)) + arb(43) / 100 * mass_tot)", 0),
 ('V37', 'displacement integral with the current stage H only (the legacy accounting)', "DM = up(Wq + Hst[st] * Qs)", "DM = up(Hst[st] * (Wq + Qs))", 0),
 ('V38', 'displacement integral not accumulated', "                Wq = DM\n", "                Wq = arb(0)\n", 0),
 ('V39', 'drift without the stage factor H', "dsum = dsum + up((f1 - f0) * dtp * Hst[st] * muc * (part1 + part2))", "dsum = dsum + up((f1 - f0) * dtp * muc * (part1 + part2))", 0),
 ('V40', 'drift part2 with L instead of L^2', "part2 = sqn(Lb) * ej[j] * hs", "part2 = Lb * ej[j] * hs", 0),
 ('V41', 'H of the non-entrant: without the omega term', "Hst = [up(h_[0] + arb(43) / 100 * amin * omega_pk) for h_ in HO]", "Hst = [up(h_[0]) for h_ in HO]", 0),
 ('V42', 'first order condition without L', "ordok = float(lo(arb(float(taus[1])))) >= float(up(bmax_up((d1 + De[0]) / p[0], arb(0)) * Lb))", "ordok = float(lo(arb(float(taus[1])))) >= float(up(bmax_up((d1 + De[0]) / p[0], arb(0))))", 0),
 ('V43', 'window order condition: the next window against the layer length only', "hi_ = bmax_up((De[k] - d1) / p[k], arb(0)) + Tlay", "hi_ = Tlay", 0),
 ('V44', 'dc condition on the lower end', "cond['dc'] = float(up(dcmax)) <= float(dc_f)", "cond['dc'] = float(lo(dcmax)) <= float(dc_f)", 0),
 ('V45', 'Omega condition dropped', "cond['Omega'] = float(omega_pk) <= float(Om_f)", "cond['Omega'] = True", 0),
 ('V46', 'strictness factor of the fixed point dropped', "if not (float(Etot) < float(Dv[j]) / 1.0001):", "if not (float(Etot) < float(Dv[j]) * 1.0001):", 0),
 ('V47', 'E_j without the drift', "Etot = up(cum[j] + mass_tot + drift[j])", "Etot = up(cum[j] + mass_tot)", 0),
 ('V48', 'E_j without the mass', "Etot = up(cum[j] + mass_tot + drift[j])", "Etot = up(cum[j] + drift[j])", 0),
 ('V49', 'theta condition: peak against 2 theta', "cond['theta'] = bool(max(float(up(pk[j])) for j in non) < theta and", "cond['theta'] = bool(max(float(up(pk[j])) for j in non) < 2 * theta and", 0),
 ('V50', 'tail ball [0, h_hi] replaced by the point h_hi', "hb = arb(0).union(A(h_hi))", "hb = A(h_hi)", 0),
 ('V51', 'candidate iteration: strictness factor dropped', "new = [Dv[j] if j in ent else float(fac * info['Etot'][j]) for j in range(J)]", "new = [Dv[j] if j in ent else float(0.5 * info['Etot'][j]) for j in range(J)]", 0),
]
