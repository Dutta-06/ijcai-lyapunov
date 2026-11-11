# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 129's revision (referee round 129, D3dx-1 … -6): A6C
# (the exact errors in Arb at every K = 159 … 400: `a6c_v`, `a6c_nodes`, `a6c_error`, `a6c_str`, `a6c_job`, `part_A6C`, main's part),
# its analysis (`a6c_win`, `a6c_analysis`: the orders on K ≥ 160, the residual floor, the coefficients and their growth, c₂ fitted, the
# outlier with q^K(B₀ + B₁/K [+ B₂/K²])) and its line (`a6c_lines`: the cross-check); in `a6x_lines` the window verdict (D3dx-1), the
# transient carried only above the order's own remainder and on the order's own K (`a6x_carry`, D3dx-2), the next term c₂ (`a6x_c2`,
# `a6x_branches`; its fit and the check with c₂ fixed, D3dx-3), the outlier's criterion and carried orders (D3dx-4), D61_STATUS (D3dx-5);
# the referee's Y01, Y07 paths (D3dx-6). N00 is the control.
M = [
 ('N00', 'control', "def a6x_carry(q, K, rem):", "def a6x_carry(q, K, rem):", 0),
 # A6C: the exact errors (D3dx-1)
 ('N01', "a6c_v: the Mercer frame's variance 2τ²/(b + 2)", "2 * t2 / (b + 1)][i] if i < 4", "2 * t2 / (b + 2)][i] if i < 4", 0),
 ('N02', "a6c_nodes: an odd rule's middle node LAPACK's, unrefined", "            x = arb(0)\n        else:", "            x = arb(float(x0))\n        else:", 0),
 ('N03', "a6c_nodes: the weights by He_K", "ws.append(fK / (K * he(x)[0]) ** 2)", "ws.append(fK / (K * he(x)[1] + K) ** 2)", 0),
 ('N04', "a6c_error: the amplitudes' exponent with + v/τ²", "g = ss2.sqrt(), v.sqrt(), (1 - v / t2) / 2", "g = ss2.sqrt(), v.sqrt(), (1 + v / t2) / 2", 0),
 ('N05', "a6c_str: one digit fewer", "        c.prec = digits\n        d = Decimal(num) / Decimal(den)", "        c.prec = digits - 1\n        d = Decimal(num) / Decimal(den)", 0),
 ('N06', "a6c_job: half the margin bits", "bits.append(int(math.ceil(2 * K * math.log2(1 / rate))) + A6C_BITS)", "bits.append(int(math.ceil(2 * K * math.log2(1 / rate))) + A6C_BITS // 2)", 0),
 ('N07', "a6c_job: flint's precision not restored", "    finally:\n        flint.ctx.prec = p0", "    finally:\n        pass", 0),
 ('N08', "a6c_job: the range from A6C_K[0] (no K₀ − 1)", "Ks, es, bits, rad = list(range(A6C_K[0] - 1, A6C_K[1] + 1))", "Ks, es, bits, rad = list(range(A6C_K[0], A6C_K[1] + 1))", 0),
 ('N09', "main: A6C not a part", "                  ('A6C', part_A6C)):", "                  ('A6C', lambda: None)):", 0),
 ('N10', "main --from-raw: A6C not read", "    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X', 'A6C')", "    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X')", 0),
 # A6C's analysis
 ('N11', "a6c_win: uniform in K, not in 1/K", "{int(round(1 / u)) for u in np.linspace(1 / K0, 1 / K1, N)}", "{int(round(u)) for u in np.linspace(K0, K1, N)}", 0),
 ('N12', "the residual floor over the nodes only", "                for K in R:\n                    u, p = Decimal(1) / K, Decimal(0)", "                for K in Ks:\n                    u, p = Decimal(1) / K, Decimal(0)", 0),
 ('N13', "the residual floor one power short", "mx = max(mx, abs(ed[K] - p * u * u))", "mx = max(mx, abs(ed[K] - p * u))", 0),
 ('N14', "the coefficients' spread threshold ×10", "            if sp > A6C_SPREAD:", "            if sp > 10 * A6C_SPREAD:", 0),
 ('N15', "the growth over the last three ratios", "g4=(min(rat[-4:]), max(rat[-4:])) if len(rat) >= 4", "g4=(min(rat[-3:]), max(rat[-3:])) if len(rat) >= 4", 0),
 ('N16', "c₂ fitted: the smallest fit, not the median", "c2q = (c2s[(len(c2s) - 1) // 2] + c2s[len(c2s) // 2]) / 2", "c2q = c2s[0]", 0),
 ('N17', "the transient's rule analysed as a power rule", "    if not a6x_resolvable(float(qd), K0):", "    if True:", 0),
 ('N18', "B₁/B₀ printed as B₁", "b10=bj[-1][1] / bj[-1][0])", "b10=bj[-1][1])", 0),
 ('N19', "the stored digits' rounding ten times larger", "dE = {K: Fraction(1, 2 * 10 ** (A6C_DIGITS - 1))", "dE = {K: Fraction(1, 2 * 10 ** (A6C_DIGITS - 2))", 0),
 ('N20', "the memo keyed on the rule alone", "    key = (x['r'], x['name'], x['K0'], tuple(x['e']), x.get('rad'))", "    key = (x['r'], x['name'])", 0),
 ('N21', "a6c_lines: the cross-check's K uncollapsed", "    kr = sorted(K for K in kx if K - 1 in kx)", "    kr = sorted(kx)", 0),
 ('N22', "a6c_lines: the 26-digit rounding without the half", "rnd = max(rnd, 0.5 * 10.0 ** (1 - A6X_DIGITS)", "rnd = max(rnd, 10.0 ** (1 - A6X_DIGITS)", 0),
 # the window (D3dx-1)
 ('N23', "the window wherever A6C exists", "            w_ = [z for z in zs if conv(z) and z['a6c'].get('g4') is not None]", "            w_ = [z for z in zs if z.get('a6c') is not None and z['a6c'].get('g4') is not None]", 0),
 ('N24', "the controls' growth over every rule", "        oth = [z for z in ok_ if z.get('a6c') is not None", "        oth = [z for z in nofit if z.get('a6c') is not None", 0),
 ('N48', "the shrink factor at the seven-term order's K", "        k8 = (list(A6_K) + list(A6X_K))[-(max(A6X_TERMS) + 1)]", "        k8 = (list(A6_K) + list(A6X_K))[-max(A6X_TERMS)]", 0),
 ('N25', "q^160 printed wherever A6C exists", "q^{A6C_K[0]} = {z['a6c']['qK0']:.1e}\" if conv(z) else ''))", "q^{A6C_K[0]} = {z['a6c']['qK0']:.1e}\" if z.get('a6c') is not None else ''))", 0),
 # the transient's carry (D3dx-2)
 ('N26', "a6x_carry: a tenth of the remainder", "    return a6x_resolvable(q, K) and rem is not None and q ** K > rem", "    return a6x_resolvable(q, K) and rem is not None and q ** K > 0.1 * rem", 0),
 ('N27', "a6x_carry: resolvability dropped", "    return a6x_resolvable(q, K) and rem is not None and q ** K > rem", "    return rem is not None and q ** K > rem", 0),
 ('N28', "the remainder one power lower (K^{−n−2})", "/ k0s[-1] ** (n + 1)) if n + 1 in cc else None)", "/ k0s[-1] ** (n + 2)) if n + 1 in cc else None)", 0),
 ('N29', "the remainder from the next order's B", "rem.append(float(abs(cc[n + 1][0][-2] * tx)", "rem.append(float(abs(cc[n + 1][0][-1] * tx)", 0),
 ('N30', "the carried order keeps the top power and adds B (non-square)", "[one] + [pw(e) for e in range(2, n + 1)] + [qK])", "[one] + [pw(e) for e in range(2, n + 2)] + [qK])", 0),
 ('N31', "the transient's rule gated at the five-term order's K (Y07)", "    trr = a6x_resolvable(q, Ks[-(min(A6X_TERMS) + 1)])", "    trr = a6x_resolvable(q, Ks[-(min(A6X_TERMS) + 2)])", 0),
 ('N32', "the printed B: the dropped power's slot", "Bw=[float(c[-1]) if r_ else None", "Bw=[float(c[-2]) if r_ else None", 0),
 ('N33', "the negative result at every rule not carried", " and any(a6x_resolvable(z['q'], k0) for k0 in z['k0s'])]", "]", 0),
 ('N34', "'not carried' printed though unresolvable", " if not r_ and a6x_resolvable(z['q'], k0)]", " if not r_]", 0),
 # the outlier (D3dx-4)
 ('N35', "the outlier's fit the four-term carried order", "                ct = cw[1]  ", "                ct = cw[0]  ", 0),
 ('N36', "the outlier's deviation from the four-term carried order", "tr = dict(q=q, B=float(ct[4]), dev=ow[1], dev5=ow[2], bnd=ob[1],", "tr = dict(q=q, B=float(ct[4]), dev=ow[0], dev5=ow[1], bnd=ob[0],", 0),
 ('N37', "the outlier's carried orders called decreasing", "flat=max(abs(b_ / float(ct[4]) - 1) for b_ in Bk), mono=a6x_decreasing(ow, ob))", "flat=max(abs(b_ / float(ct[4]) - 1) for b_ in Bk), mono=True)", 0),
 # the next term (D3dx-3)
 ('N38', "a6x_c2: ϱ's sign by the other branch", "    g = t if fast else -t", "    g = -t if fast else t", 0),
 ('N39', "a6x_c2: the √q term over 4", "sq * (1 - g) ** 2 / (8 * (1 - sq)", "sq * (1 - g) ** 2 / (4 * (1 - sq)", 0),
 ('N40', "a6x_c2: 1/12 for 1/16", "    return (one / 16 - g", "    return (one / 12 - g", 0),
 ('N41', "a6x_c2: the cusp guard at 10⁻⁸⁰", "    if abs(one - sq) < one / 10 ** 20:", "    if abs(one - sq) < one / 10 ** 80:", 0),
 ('N42', "c₂ from the stored errors: six powers", "                m = max(A6X_TERMS)\n", "                m = max(A6X_TERMS) - 1\n", 0),
 ('N43', "the check with c₂ fixed: the shift as K⁻³", "Rm = {K: R[K] - c2 / (tx * K * K) for K in Ks}", "Rm = {K: R[K] - c2 / (tx * K * K * K) for K in Ks}", 0),
 ('N44', "the check with c₂ fixed: cons always", "                z['cons2'] = abs(o2[-1]) <= z['res2']", "                z['cons2'] = True", 0),
 ('N45', "c₂'s fit never A6C's", "                if z['a6c'] is not None and z['a6c'].get('c2rel') is not None:", "                if False:", 0),
 ('N46', "the known cases with c₂ fixed: every rule", "            kp2 = [z for z in cz if z['fr']]", "            kp2 = cz", 0),
 # D61_STATUS (D3dx-5)
 ('N47', "D61_STATUS without Arb's correctness", " — all rigorous in ball arithmetic, Arb's correctness assumed\")", " — all rigorous in ball arithmetic\")", 0),
]
