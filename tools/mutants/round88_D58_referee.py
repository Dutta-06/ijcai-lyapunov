# Referee round 88 (D3ci): single-point mutants of round 84's new or changed code in exp_width_mismatch.py, beyond the round's own lists.
M = [
 ('Z00', 'control', 'A6_MARGIN = 60 ', 'A6_MARGIN = 60 ', 0),
 # precision bookkeeping of the new exact computations (the real-job tests run only at small K)
 ('P01', 'a6_error: nodes at 110 digits', "        zs, ws = gh_nodes(K, prec)\n        sv = v.sqrt()", "        zs, ws = gh_nodes(K)\n        sv = v.sqrt()", 0),
 ('P02', 'gen_ratio: nodes at 110 digits', "        zs, ws = gh_nodes(K, prec)\n        x = Decimal(repr(s))", "        zs, ws = gh_nodes(K)\n        x = Decimal(repr(s))", 0),
 ('P03', 'a6_job: precision from t(v), not t', "prec = int(math.ceil(2 * K * math.log10(1 / t))) + A6_MARGIN", "prec = int(math.ceil(2 * K * math.log10(1 / rate_conj(r_, float(a6_scalings(r_)[i][1]))[0]))) + A6_MARGIN", 0),
 ('P04', 'a6_error: the context at half the precision', "        ctx.prec = prec\n        v = +v", "        ctx.prec = prec // 2\n        v = +v", 0),
 ('P05', 'gh_nodes: 6 Newton steps', "for _ in range(8):", "for _ in range(6):", 0),
 ('P06', 'part_A6E: 2K + 60 digits for K = 80 only', "ratio=[gen_ratio(K, s, ax, 3 * K + 60) for s in A6E_S]", "ratio=[gen_ratio(K, s, ax, 3 * K + 60 if K < 80 else 2 * K + 60) for s in A6E_S]", 0),
 ('P07', 'a5_far_row: S0 at prec - 40', "S0 = frame_series(K - 1, q, 1 + q, prec=prec)", "S0 = frame_series(K - 1, q, 1 + q, prec=prec - 40)", 0),
 # a6_lines
 ('L01', 'a6_lines: exponent from K = 80', "ex = -float(np.log(y1 / y0) / np.log(Ks[-1] / Ks[-2])) if y1 * y0 > 0 else None", "ex = -float(np.log(y1 / y0) / np.log(Ks[-1] / Ks[-3])) if y1 * y0 > 0 else None", 0),
 ('L02', 'a6_lines: near by <=', "near = [z for z in dev if abs(np.log(z[0]['v_hs'])) < A6_NEAR]", "near = [z for z in dev if abs(np.log(z[0]['v_hs'])) <= A6_NEAR]", 0),
 ('L03', 'a6_lines: the farthest by the least', "w = max(dev, key=lambda z: abs(z[1]))", "w = min(dev, key=lambda z: abs(z[1]))", 0),
 ('L04', 'a6_lines: cusp residual 1/(2K)', "res = [y['ratio'] * x['rate'] / (1 - 1 / (4 * K)) - 1 for y, K in zip(x['rows'], Ks)]", "res = [y['ratio'] * x['rate'] / (1 - 1 / (2 * K)) - 1 for y, K in zip(x['rows'], Ks)]", 0),
 ('L05', 'a6_lines: the cusp within, by min', "f\"; at v_hs within {max(abs(z[1]) for z in cusp):.1e} of 1/t\"", "f\"; at v_hs within {min(abs(z[1]) for z in cusp):.1e} of 1/t\"", 0),
 ('L06', 'a6_lines: Mercer check squared error unsquared', "e2[K_] = float(e_) ** 2 / sr[x['r']][K_]", "e2[K_] = float(e_) / sr[x['r']][K_]", 0),
 ('L07', 'a6_lines: excess at K = 120 printed', "(cusp if at else dev).append((x, L * x['rate'] - 1, y1, ex))", "(cusp if at else dev).append((x, L * x['rate'] - 1, y0, ex))", 0),
 ('L08', 'a6_lines: the limit fit through K = 60, 80, 120', "L = float(np.linalg.solve(np.array([[1.0, K ** p, K ** q] for K in Ks[-3:]]), rr[-3:])[0])", "L = float(np.linalg.solve(np.array([[1.0, K ** p, K ** q] for K in Ks[-4:-1]]), rr[-4:-1])[0])", 0),
 ('L09', 'a6_lines: far filter by arg ignored', "far = lambda a_: [z for z in dev if z[0]['arg'] == a_ and abs(np.log(z[0]['v_hs'])) >= A6_NEAR]", "far = lambda a_: [z for z in dev if abs(np.log(z[0]['v_hs'])) >= A6_NEAR]", 0),
 # rate_conj, a6_scalings
 ('R01', 'rate_conj: ties to the second argument', "return (a, 1) if a >= b else (b, 2)", "return (a, 1) if a > b else (b, 2)", 0),
 ('R02', 'a6_scalings: Mercer-zero over b²', "('Mercer-zero', tau2 / b)", "('Mercer-zero', tau2 / (b * b))", 0),
 ('R03', 'a6_scalings: v_hs = τ²σ*', "    vhs = tau2 / ss\n    return", "    vhs = tau2 * ss\n    return", 0),
 # a6_job
 ('J01', 'a6_job: e0 to 20 digits', "e0=format(E0.sqrt(), '.25e')", "e0=format(E0.sqrt(), '.20e')", 0),
 ('J02', 'a6_job: v_hs ratio without the square root', "vv, vr = float(v), float(v / (tau2 / (1 + tau2).sqrt()))", "vv, vr = float(v), float(v / (tau2 / (1 + tau2)))", 0),
 ('J03', 'a6_job: ok by E1 only', "ok = E0 > 0 and E1 > 0", "ok = E1 > 0", 0),
 # gen_ratio, dpi, dcos, a6e_lines
 ('G01', 'gen_ratio: the sign (−1)^K dropped', "return float(e / ((-1) ** K * a * ((x * x) / 4).exp()))", "return float(e / (a * ((x * x) / 4).exp()))", 0),
 ('G02', 'dcos: no reduction modulo 2π', "    t = t % (2 * pi)\n", "", 0),
 ('G03', 'gen_ratio: real axis e^{-ξ²/2} for e^{-ξ²/4}', "return float(e / (a * (-(x * x) / 4).exp()))", "return float(e / (a * (-(x * x) / 2).exp()))", 0),
 ('G04', 'a6e_lines: values at s = 1.9 for 2', "for s in (0.5, 1.0, 1.5, 2.0)) + ')' for x in A6E)", "for s in (0.5, 1.0, 1.5, 1.9)) + ')' for x in A6E)", 0),
 ('G05', 'A6E_S from 0.2', "A6E_K, A6E_S = (20, 40, 80), tuple(i / 10 for i in range(1, 21))", "A6E_K, A6E_S = (20, 40, 80), tuple(i / 10 for i in range(2, 21))", 0),
 # b5_lines, b5_jobs, b5_job
 ('F01', 'b5_lines: per-start spread from the least', "v, c_, k_ = max(sp_all, key=lambda z: z[0])", "v, c_, k_ = min(sp_all, key=lambda z: z[0])", 0),
 ('F02', 'b5_lines: anchored r0 from the other start', "qb = next((q for q in B4 if (q['r'], q['k']) == c_ and q['start'] == best[c_]['start']), None)", "qb = next((q for q in B4 if (q['r'], q['k']) == c_ and q['start'] != best[c_]['start']), None)", 0),
 ('F03', 'b5_lines: flat count by not stall', "f\"{fr(k_)}: {sum(1 for x in anc for y in x['rows'] if y['kappa'] == k_ and y['stop'] == 'flat')} of \"", "f\"{fr(k_)}: {sum(1 for x in anc for y in x['rows'] if y['kappa'] == k_ and y['stop'] != 'stall')} of \"", 0),
 ('F04', 'b5_lines: the spread over κ ≥ 1/16 only', "sp_all = [(v, c_, k_) for c_, n, g, a, b, sp in gm for v, k_ in zip(sp, B5_KAPPA)]", "sp_all = [(v, c_, k_) for c_, n, g, a, b, sp in gm for v, k_ in zip(sp, B5_KAPPA) if k_ > 0]", 0),
 ('F05', 'b5_job: e_start from lm_rms', "e_start=exact_rms(f, np.r_[w0, c0], k) / rf / e_opt, rows=rows,", "e_start=lm_rms(f, np.r_[w0, c0], k) / rf / e_opt, rows=rows,", 0),
 ('F06', 'b5_lines: the κ = 1 check skips the block', "st = {0.0: 'optimum', 1.0: 'block'}.get(x['theta'])", "st = {0.0: 'optimum'}.get(x['theta'])", 0),
 ('F07', 'kappa_blend: κ in double', "return (rl + LD(kappa) * (np.asarray(rd).astype(LD) - rl)).astype(np.float64)", "return (rl + LD(float(np.float32(kappa))) * (np.asarray(rd).astype(LD) - rl)).astype(np.float64)", 0),
 ('F08', 'lm_fit: kappa 0 treated as None', "mix = (lambda rd, p_: rd) if kappa is None else", "mix = (lambda rd, p_: rd) if not kappa else", 0),
 # b4_lines' new summaries, b4_floors, anchor_info, exact_resid's decimal
 ('Q01', 'b4_lines: the gap-ratio threshold 10² for 10³', "eg = [ex(q, k_) / q[k_]['gap'] for q in B4 for k_, _ in B4_ENDS if ex(q, k_) >= 1e3 * max(gaps)]", "eg = [ex(q, k_) / q[k_]['gap'] for q in B4 for k_, _ in B4_ENDS if ex(q, k_) >= 1e2 * max(gaps)]", 0),
 ('Q02', 'b4_lines: last-anchor kinds miss start', "for k_ in ('stall', 'period', 'start') if any(x['kind'] == k_ for x in an))", "for k_ in ('stall', 'period') if any(x['kind'] == k_ for x in an))", 0),
 ('Q03', 'b4_lines: double over long double from the blocks only', "dl = [(q, ex(q, 'double'), ex(q, 'long double')) for q in B4]", "dl = [(q, ex(q, 'double'), ex(q, 'long double')) for q in blk]", 0),
 ('Q04', 'b4_floors: e_end against the last anchor exact residual', "e_end = [Decimal(float(a)) - b for a, b in zip(direct_double(f, p, k), r)]", "e_end = [Decimal(float(a)) - b for a, b in zip(direct_double(f, p, k), ra)]", 0),
 ('Q05', 'exact_parts: sg2 in double', "        sg2 = Decimal(SIG) ** 2\n        s2 = 2 * sg2\n        r, J = [], []", "        sg2 = Decimal(float(np.float64(SIG) ** 2))\n        s2 = 2 * sg2\n        r, J = [], []", 0),
 ('Q06', 'b4_job: rms_d at the long double end itself', "rms_d=str(exact_resid(f, pl.astype(float), k, dec=True)[2])", "rms_d=str(exact_resid(f, pl, k, dec=True)[2])", 0),
 ('Q07', 'lm_fit: the anchor log copies r0 by reference', "anchor_log.update(p=p_.copy(), r0=anc['r'].copy(), it=it_, kind=kind)", "anchor_log.update(p=p_, r0=anc['r'], it=it_, kind=kind)", 0),
 ('Q08', 'b4_lines: d2 over 2σ² in double', "d2 = float((Decimal(float(2 * np.float64(SIG) ** 2)) - 2 * Decimal(SIG) ** 2) / (2 * Decimal(SIG) ** 2))", "d2 = float((Decimal(float(2 * np.float64(SIG) ** 2)) - 2 * Decimal(SIG) ** 2) / Decimal(float(2 * np.float64(SIG) ** 2)))", 0),
 # a5 far line
 ('A01', 'a5 far line: check compared to 0 by <', "f\"; the last K with {A5_FAR_CHECK} more digits the same in double\" if max(ck) == 0 else", "f\"; the last K with {A5_FAR_CHECK} more digits the same in double\" if max(ck) < 1e-12 else", 0),
 ('A02', 'a5 far line: within, by min', "+ f\", within {max(abs(L / c - 1) for _, L, c, _ in lims):.1e} of the conjecture's.\"]", "+ f\", within {min(abs(L / c - 1) for _, L, c, _ in lims):.1e} of the conjecture's.\"]", 0),
 ('A03', 'frame_series: the guard 10 terms → 3', "if len(terms) > 10 and 0 < max(terms[-2:]) < tol * S:", "if len(terms) > 3 and 0 < max(terms[-2:]) < tol * S:", 0),
]
