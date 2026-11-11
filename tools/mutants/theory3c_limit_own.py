# (id, what, old, new, k): single-point mutants of experiments/theory3c_limit.py (rigorous enclosure of the limit object of Theorem T1; round theory3, direction C):
# the Krawczyk operator (interior test, Jacobian entries, preconditioner, radius term, contraction), the strict-argmin decisions, the nearest-feature and sign tests.
# Run: tools/mutants.py --exp theory3c_limit --tests tests/theory3c_test.py --raw results/synthetic/theory3c_raw.npz --md results/synthetic/theory3c.md --out-env T3C_OUT
#      --extra results/synthetic/theory3c.md --spec tools/mutants/theory3c_limit_own.py --work <scratch dir>
# (the criterion is the tests file alone; `--from-raw` of a library module regenerates nothing, the md column is not informative)
# Result (44 mutants, criterion: the tests file alone): 39 killed.  Survivors, equivalent on the data (dead or implied branches):
#   D04 (strictness against the lower bounds of undecided-rate kernels: those bounds are about 10^60 times the entrant's time at 220 bits);  D09 (the wall margin |u_j| < 1: implied by the strict
#   argmin, since a dormant kernel with |u_j| >= 1 would have been the entrant);  D16 (undecided-rate kernels ignored: same reason as D04);  D19 (ok without len(events) == K: a break before
#   K events always sets a failed check);  D20 (final-stage rates replaced by 0: the final block is the target, its residual rates are 0 to 10^-60).
M = [
 ('H00', 'control (unmutated)', '', '', 0),
 ('K01', 'interior test: lower end dropped', "Kl[i].lower() > X[i].lower() and Kl[i].upper() < X[i].upper()", "Kl[i].upper() < X[i].upper()", 0),
 ('K02', 'interior test: upper end dropped', "Kl[i].lower() > X[i].lower() and Kl[i].upper() < X[i].upper()", "Kl[i].lower() > X[i].lower()", 0),
 ('K03', 'interior test: any coordinate instead of all', "inside = all(Kl[i]", "inside = any(Kl[i]", 0),
 ('K04', 'interior test: reversed', "Kl[i].lower() > X[i].lower() and Kl[i].upper() < X[i].upper()", "Kl[i].lower() < X[i].lower() and Kl[i].upper() > X[i].upper()", 0),
 ('K05', 'radius term dropped (Newton image only)', "d = _mat_of_col([X[i] - z0[i] for i in range(n)])", "d = _mat_of_col([arb(0) for i in range(n)])", 0),
 ('K06', 'J(X) replaced by J(z0)', "_, JX = FJ(X[:k], X[k:], cs, As)", "JX = J0", 0),
 ('K07', 'no preconditioner (C = I)', "Cm = J0.mid().inv().mid()", "Cm = arb_mat(n, n, [arb(1) if i == j else arb(0) for i in range(n) for j in range(n)])", 0),
 ('K08', 'Newton term sign', "Kv = _mat_of_col(z0) - Cm * _mat_of_col(F0) + Mx * d", "Kv = _mat_of_col(z0) + Cm * _mat_of_col(F0) + Mx * d", 0),
 ('K09', 'contraction dropped', "inter = a_.intersection(b_)", "inter = b_", 0),
 ('K10', 'box placement without polishing (equivalent: the boxes grow)', "z0 = newton_polish(z0, cs, As)", "pass", 0),
 ('J01', 'Jacobian dg_i/dw_i', "ent[i][i] = arb(-1)", "ent[i][i] = arb(1)", 0),
 ('J02', 'Jacobian dg_i/dc_i', "ent[i][k + i] = h[i]", "ent[i][k + i] = arb(0)", 0),
 ('J03', 'Jacobian dg_i/dc_j sign', "ent[i][k + j] = w[j] * dGam(x)", "ent[i][k + j] = -w[j] * dGam(x)", 0),
 ('J04', 'Jacobian dh_i/dw_j sign', "ent[k + i][j] = -dGam(x)", "ent[k + i][j] = dGam(x)", 0),
 ('J05', 'Jacobian dh_i/dc_j sign', "ent[k + i][k + j] = w[j] * d2Gam(x)", "ent[k + i][k + j] = -w[j] * d2Gam(x)", 0),
 ('J06', 'Jacobian dh_i/dc_i cross term sign', "h2 = h2 - w[j] * d2Gam(x)", "h2 = h2 + w[j] * d2Gam(x)", 0),
 ('J07', 'Jacobian dg_i/dw_j sign', "ent[i][j] = -Gam(x)", "ent[i][j] = Gam(x)", 0),
 ('F01', 'second derivative of the Gaussian', "return (sqb(x) / 4 - arb(1) / 2) * Gam(x)", "return (arb(1) / 2 - sqb(x) / 4) * Gam(x)", 0),
 ('F02', 'Gaussian width', "return (-sqb(x) / 4).exp()", "return (-sqb(x) / 2).exp()", 0),
 ('F03', 'h without the dictionary term sign', "hi = hi + As[m] * dGam(x)", "hi = hi - As[m] * dGam(x)", 0),
 ('F04', 'g: the unit\'s own term dropped', "gi = gi - w[i]                     # G(0) = 1", "gi = gi", 0),
 ('F05', 'sq of a ball containing 0: negative lower end', ".nonnegative_part()", "", 0),
 ('D01', 'entrant: the largest hitting time instead of the smallest', "j1 = min(times, key=lambda j: float(times[j][0].mid()))", "j1 = max(times, key=lambda j: float(times[j][0].mid()))", 0),
 ('D02', 'strictness against the other decided times: tie accepted', "if j != j1 and not (dt1.upper() < dt.lower()):", "if j != j1 and not (dt1.upper() <= dt.upper()):", 0),
 ('D03', 'strictness against the other decided times dropped', "if j != j1 and not (dt1.upper() < dt.lower()):", "if False:", 0),
 ('D04', 'strictness against the undecided-rate lower bounds dropped', "if not (dt1.upper() < lb):", "if False:", 0),
 ('D05', 'hitting time without the dual: (1 - u)', "times[j] = ((1 - s * u[j]) / (s * p[j]), s)", "times[j] = ((1 - s) / (s * p[j]), s)", 0),
 ('D06', 'entry time not accumulated', "tau = tau + dt1", "tau = dt1", 0),
 ('D07', 'entrant dual not set to the wall', "un[j1] = arb(s1)", "un[j1] = arb(0)", 0),
 ('D08', 'duals not advanced', "un = [u[j] + dt1 * p[j] for j in range(J)]", "un = [u[j] for j in range(J)]", 0),
 ('D09', 'wall margin test dropped (implied by the strict argmin: equivalent on data)', "if not (abs(un[j]).upper() < 1):", "if False:", 0),
 ('D10', 'nearest feature: the farthest', "mbest = min(dist, key=lambda m: float(dist[m].mid()))", "mbest = max(dist, key=lambda m: float(dist[m].mid()))", 0),
 ('D11', 'nearest feature: tie accepted', "if m != mbest and not (dist[mbest].upper() < dist[m].lower()):", "if False:", 0),
 ('D12', 'sign test dropped', "sgn_ok = (s1 == (1 if float(As_f[f]) > 0 else -1))", "sgn_ok = True", 0),
 ('D13', 'sign test reversed', "sgn_ok = (s1 == (1 if float(As_f[f]) > 0 else -1))", "sgn_ok = (s1 != (1 if float(As_f[f]) > 0 else -1))", 0),
 ('D14', 'rate of the active kernels not zeroed', "                p.append(arb(0))\n                continue", "                pass", 0),
 ('D15', 'rate without the block term', "pj = pj - w[i] * Gam(c0[j] - c[i])", "pj = pj", 0),
 ('D16', 'undecided-rate kernels ignored (their lower bound is astronomically large: equivalent on data)', "lows[j] = den / pu if pu > 0 else arb.pos_inf()", "lows[j] = arb.pos_inf()", 0),
 ('D17', 'the block of the wrong feature set', "wf, cf = crit(np.asarray(cs_f, float), np.asarray(As_f, float), assigned)", "wf, cf = crit(np.asarray(cs_f, float), np.asarray(As_f, float), assigned[:1])", 0),
 ('D18', 'limit declared ok without the checks', "ok = (len(events) == K) and all(chk.values())", "ok = (len(events) == K)", 0),
 ('D19', 'limit declared ok with fewer events', "ok = (len(events) == K) and all(chk.values())", "ok = all(chk.values())", 0),
 ('D20', 'final-stage rates dropped', "        P.append(rates(w, c, active))", "        P.append([arb(0)] * J)", 0),
 ('D21', 'offset d_0 of the wrong kernel', "d0=abs(c0[j1] - cs[f])", "d0=abs(c0[j1] - cs[0])", 0),
 ('D22', 'p_e without the absolute value', "p_e=abs(p[j1])", "p_e=p[j1]", 0),
]
