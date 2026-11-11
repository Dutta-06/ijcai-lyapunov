# (id, what, old, new, k): single-point mutants of experiments/theory3e_kernels.py (the box conditions of Proposition G, the tail envelopes, the own constants, the cusp conditions, the clean corollary).
# Run: venv/bin/python tools/mutants.py --exp theory3e_kernels --tests tests/test_theory3e_core.py --raw results/synthetic/theory3e_raw.npz --md results/synthetic/theory3e.md --out-env T3E_OUT \
#        --extra results/synthetic/theory3e.md --spec tools/mutants/theory3e_kernels_own.py --work <scratch dir>
# (the raw and md are those of theory3e_run.py, which the mutated module does not regenerate; the criterion is the tests file alone.)
# Result (57 mutants + the control K00, criterion: tests/test_theory3e_core.py alone): 51 killed. 4 first survived and were killed by tests added for them: K09 (d-tilde < r as a condition of its own, stub with
#   lambda_psi(r) < psi(r)/r), K22 (the kappa < 1 guard), T03 (the tail value at L_MAX), C03 (d_1 < r in the clean corollary). Survivors, equivalent on the data:
#   T05 (Hessian-norm envelope |Gamma'|/d with the cell's right end: compensated by the looseness of the ball enclosure of the product, difference < 1e-3);  T07 (the 'lo' helper returning the midpoint: the
#   balls have radius ~1e-16);  H05 (cusp E_c from m_0 instead of m_1: for Gamma = exp(-d) they are the same function);  C06 (cap c_1 <= 1/4 -> 1: c_* ~ 0.035 never reaches the cap);  C07 (decay modulus
#   with the cell's left end: changes C by <= 1.5%, below the tests' 10% window);  C08 (no rounding of the clean threshold up to the 0.02 grid: within the grid resolution).
M = [
 ('K00', 'control (unmutated)', 'def up(x):', 'def up(x):', 0),
 # ---- box(): formulas
 ('K01', 'vbar: denominator 1 + kappa', 'vbar = (Rr * (one + S0) + Om) / (one - S0)', 'vbar = (Rr * (one + S0) + Om) / (one + S0)', 0),
 ('K02', 'E_v without the R + vbar factor', 'Ev = S0 * (Rr + vbar) + Om', 'Ev = S0 * Rr + Om', 0),
 ('K03', 'E_v without the dormant budget', 'Ev = S0 * (Rr + vbar) + Om', 'Ev = S0 * (Rr + vbar)', 0),
 ('K04', 'E_c without R + vbar', 'Ec = S1 * (Rr + vbar) + M1 * Om', 'Ec = S1 * Rr + M1 * Om', 0),
 ('K05', 'E_c without the dormant term', 'Ec = S1 * (Rr + vbar) + M1 * Om', 'Ec = S1 * (Rr + vbar)', 0),
 ('K06', 'v_lo without E_v', 'vlo = g_r - Ev', 'vlo = g_r', 0),
 ('K07', '(C) without the condition number', 'Cm = psi_r - arb(cond_c) * Ec', 'Cm = psi_r - Ec', 0),
 ('K08', 'd-tilde divided by psi(r) instead of lambda_psi(r)', 'dt = Ec / lps_r', 'dt = Ec / psi_r', 0),
 ('K09', 'd-tilde < 2 r', 'if not (dt < r):', 'if not (dt < 2 * r):', 0),
 ('K10', 's without the amplitude ratio', 's_max = Rr * arb(kern.psihat(dtu))', 's_max = arb(kern.psihat(dtu))', 0),
 ('K11', 't without the coefficient floor', 't_rad = vlo_d * rho_d', 't_rad = rho_d', 0),
 ('K12', 'least eigenvalue with +sqrt', 'lam2 = (one + t_rad - _sqrt((one - t_rad) * (one - t_rad) + 4 * s_max * s_max)) / 2', 'lam2 = (one + t_rad + _sqrt((one - t_rad) * (one - t_rad) + 4 * s_max * s_max)) / 2', 0),
 ('K13', 'least eigenvalue without the off-diagonal s', 'lam2 = (one + t_rad - _sqrt((one - t_rad) * (one - t_rad) + 4 * s_max * s_max)) / 2', 'lam2 = (one + t_rad - _sqrt((one - t_rad) * (one - t_rad))) / 2', 0),
 ('K14', 'tangential stiffness dropped', 'if n > 1:\n        t_perp', 'if n > 99:\n        t_perp', 0),
 ('K15', 'tangential stiffness with the radial lower bound', 't_perp = vlo_d * arb(kern.lampsi(dtu))', 't_perp = vlo_d * rho_d', 0),
 ('K16', 'block-norm H_c without the diagonal block D', 'Hc = n2(arb(0), Rv * S1 + M1 * Om, vbar * Rv * S2 + vbar * M2 * Om) + n2(S0, vbar * S1, vbar * vbar * S2)', 'Hc = n2(S0, vbar * S1, vbar * vbar * S2)', 0),
 ('K17', 'block-norm H_c without the pair block B', 'Hc = n2(arb(0), Rv * S1 + M1 * Om, vbar * Rv * S2 + vbar * M2 * Om) + n2(S0, vbar * S1, vbar * vbar * S2)', 'Hc = n2(arb(0), Rv * S1 + M1 * Om, vbar * Rv * S2 + vbar * M2 * Om)', 0),
 ('K18', 'norm of the 2x2 block: max replaced by sum of the diagonal only', 'return (a + c + _sqrt((a - c) * (a - c) + 4 * b * b)) / 2', 'return (a + c + _sqrt((a - c) * (a - c))) / 2', 0),
 ('K19', 'row-sum H_c: the second row dropped', 'Hc = rowv if float(rowv.mid()) > float(rowd.mid()) else rowd', 'Hc = rowv', 0),
 ('K20', 'row-sum row d without the vbar S1 term', 'rowd = Rv * S1 + vbar * Rv * S2 + vbar * S1 + vbar * vbar * S2 + M1 * Om + vbar * M2 * Om', 'rowd = Rv * S1 + vbar * Rv * S2 + vbar * vbar * S2 + M1 * Om + vbar * M2 * Om', 0),
 ('K21', 'dormant M_2 term dropped (norm)', 'vbar * Rv * S2 + vbar * M2 * Om)', 'vbar * Rv * S2)', 0),
 ('K22', 'kappa < 1 dropped', 'if not (S0 < 1):\n        return out\n    M1', 'M1', 0),
 ('K23', 'n >= 2 uses the n = 1 Hessian envelope', "S2 = arb(up(arb((K - 1) * kern.tail(2 if n == 1 else 'H', ell, n))))", 'S2 = arb(up(arb((K - 1) * kern.tail(2, ell))))', 0),
 ('K24', 'sums= ignored', 'S0, S1, S2 = [arb(float(v)) for v in sums]', 'S0, S1, S2 = arb(0.001), arb(0.001), arb(0.001)', 0),
 ('K25', 'ell threshold L_MIN dropped', 'if ell <= L_MIN or r <= 0:\n        return out\n    if sums is None:', 'if r <= 0:\n        return out\n    if sums is None:', 0),
 # ---- tails
 ('T01', 'tail table: no running maximum', 't0[i], t1[i], t2[i], th[i] = max(v0, t0[i + 1]), max(v1, t1[i + 1]), max(v2, t2[i + 1]), max(vh, th[i + 1])', 't0[i], t1[i], t2[i], th[i] = v0, v1, v2, vh', 0),
 ('T02', 'tail lookup rounds up', 'i = int(math.floor((l - L_MIN) / H_TAB + 1e-12))', 'i = int(math.ceil((l - L_MIN) / H_TAB + 1e-12))', 0),
 ('T03', 'tail beyond L_MAX: zero', 'c0, c1, c2 = up(abs(self.G(xe))), up(abs(self.G1(xe))), up(abs(self.G2(xe)))', 'c0, c1, c2 = 0.0, 0.0, 0.0', 0),
 ('T05', 'cell_up: |Gamma\'|/d with the right end of the cell', 'return v0, v1, v2, max(v2, v1 / a)', 'return v0, v1, v2, max(v2, v1 / b)', 0),
 ('T06', 'cell: radius without the outward margin', 'return arb((a + b) / 2.0, (b - a) / 2.0 * (1 + 1e-12) + 1e-300)', 'return arb((a + b) / 2.0, (b - a) / 2.0 * 0.5)', 0),
 ('T07', 'lower bound helper returns the midpoint', 'return math.nextafter(float(u.mid()) - float(u.rad()), -math.inf)', 'return float(u.mid())', 0),
 # ---- own constants
 ('O01', 'g_min at 0', 'return lo(self.G(arb(x)))', 'return 1.0', 0),
 ('O02', 'psi_at sign', 'return lo(-self.G1(arb(x)))', 'return lo(self.G1(arb(x)))', 0),
 ('O03', 'cells: fixed 1 cell', 'n = int(min(max(math.ceil((b - a) / 0.0015), 24), 1200))', 'n = 1', 0),
 ('O04', 'lampsi: sup instead of inf', 'return min(lo(self.ptilde_ball(a, b)) for a, b in self._cells(0.0, x, N))', 'return max(lo(self.ptilde_ball(a, b)) for a, b in self._cells(0.0, x, N))', 0),
 ('O05', 'psihat lower end', 'return max(up(-self.G1(cell(a, b))) for a, b in self._cells(0.0, x, N))', 'return max(lo(-self.G1(cell(a, b))) for a, b in self._cells(0.0, x, N))', 0),
 ('O06', 'ptilde scaling dropped (sigma-equivalent units)', 'return self.ptilde_raw(sx) * (self.s ** 2)', 'return self.ptilde_raw(sx)', 0),
 ('O07', 'sigma-equivalent scale: s = 1/sqrt(lam0)', 'self.s = 1.0 / math.sqrt(2.0 * self.lam0)', 'self.s = 1.0 / math.sqrt(self.lam0)', 0),
 ('O08', 'Matern P coefficient (binomial factor)', 'self.P[p - i] = Fraction(f(p), f(2 * p)) * Fraction(f(p + i), f(i) * f(p - i)) * (2 ** (p - i))', 'self.P[p - i] = Fraction(f(p), f(2 * p)) * Fraction(f(p + i), f(i) * f(p - i)) * (2 ** (p - i + 1))', 0),
 ('O09', 'Bessel Gram: c_nu', 'return 2.0 ** (1 - mu) / math.gamma(mu)', 'return 2.0 ** (1 - mu) / math.gamma(mu + 1)', 0),
 ('O10', 'Bessel kernel rho_min: P lower bound at the left end', 'worst = min(worst, Plo - Nup)', 'worst = min(worst, Plo + Nup)', 0),
 ('O11', 'Bessel psihat: cell upper bound without B', 'u1 = B * Pup', 'u1 = Pup', 0),
 ('O12', 'family(): Gram index 2 nu + n/2 -> nu + n/2', 'return matern_gram(2 * nub + n / 2.0, n)', 'return matern_gram(nub + n / 2.0 + 1.0, n)', 0),
 # ---- thresholds and the cusp
 ('H01', 'threshold search keeps the largest', 'if t == t and (best != best or t < best):\n            best, br = t, r', 'if t == t and (best != best or t > best):\n            best, br = t, r', 0),
 ('H02', 'r grid starts at 0', 'return [float(x) for x in np.linspace(max(d0, 1e-3) * 1.000001 + 1e-9, rmax, n)]', 'return [float(x) for x in np.linspace(max(d0, 1e-3) * 0.5, rmax, n)]', 0),
 ('H03', 'bisection returns the last failing grid point', '    return D(b)\n\n\ndef entry_threshold(', '    return D(a)\n\n\ndef entry_threshold(', 0),
 ('H04', 'cusp: (C\') with psi at 0', 'psi_min = arb(lo((-kern.G1(arb(r)))))', 'psi_min = arb(1.0)', 0),
 ('H05', 'cusp: E_c uses m_0', 'S1 = arb(up(arb((K - 1) * kern.tail(1, ell))))\n    Rr, one = arb(R), arb(1)', 'S1 = arb(up(arb((K - 1) * kern.tail(0, ell))))\n    Rr, one = arb(R), arb(1)', 0),
 ('H06', 'cusp: (L) without E_v', 'vlo, Cm = g_r - Ev, psi_min - Ec', 'vlo, Cm = g_r, psi_min - Ec', 0),
 # ---- the clean corollary
 ('C01', 'clean: constant 12 -> 6', 'return Lam > 12.0 * c1', 'return Lam > 6.0 * c1', 0),
 ('C02', 'clean: e = 2.7 c1 -> 1.35 c1', 'ev = 2.7 * c1', 'ev = 1.35 * c1', 0),
 ('C03', 'clean: d_1 < r dropped', 'if not d1 < r:\n        return False', 'if False:\n        return False', 0),
 ('C04', 'clean: index without R^2', 'return (K - 1) * R * R * m', 'return (K - 1) * R * m', 0),
 ('C05', 'clean: index without (K-1)', 'return (K - 1) * R * R * m', 'return R * R * m', 0),
 ('C06', 'clean: c_* cap 1/4 -> 1', 'if not (c1 <= 0.25 and', 'if not (c1 <= 1.0 and', 0),
 ('C07', 'decay modulus uses the cell left end', 'best = max(best, m * math.exp(a * ls[i + 1] ** q))', 'best = max(best, m * math.exp(a * ls[i] ** q))', 0),
 ('C08', 'clean threshold: no rounding up to the grid', 't = math.ceil((ell + 2 * r) / step - 1e-9) * step', 't = (ell + 2 * r)', 0),
]
