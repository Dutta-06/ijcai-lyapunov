# (id, what, old, new, k): single-point mutants of exp_transit_certificate.py — round 140's revision (referee D3ei), never
# mutation-tested when applied: Part G (`job_2d_flow`, `part_G`: the unpolished flow on B, C, D(a), D(c)) and Part H (`eta_v`: Duval–Peyré's
# vanishing-derivative pre-certificate, its margin off E, its nondegeneracy beyond EX_R·σ; `bp_certificate`: basis pursuit's least
# certificate on the grid by a linear program; `omp`: continuous OMP without and with the slide; `h_job`, `part_H`) and the G and H
# sections of `report`.
# Result: 26 mutants and a control. Against the tests as committed: 3 killed by the tests (G07, G10, G17), 2 by the md (G25, G26); the 21
# others killed by the tests added then (test_parts_G_and_H_job_lists_are_the_stored_rows, test_h_job_reproduces_stored_baselines,
# test_eta_v_nondegenerate_needs_the_hessian, test_job_2d_flow_reports_the_unpolished_flow, test_job_d_with_the_flow_reproduces_stored_G_
# rows, test_report_G_and_H_count_planted_rows, test_omp_selects_by_the_normalised_correlation). G15 changes the stored OMP results by at
# most 5·10⁻⁵ (relative) on 4 of the 1154 targets; the planted edge feature kills it. No survivor.
M = [
 ('G00', 'control', "def _k1(K):", "def _k1(K):", 0),
 ('G01', "G's plane flow polished", "Rf = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)), polish=False)",
  "Rf = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)), polish=True)", 0),
 ('G02', "G's plane flow's outcome the path's", "out_flow=outcome(Rf, cert['e'], len(cs))", "out_flow=outcome(R, cert['e'], len(cs))", 1),
 ('G03', "G's blocks relaxed if any is", "flow_relaxed=bool(all(Rf['relaxed']))", "flow_relaxed=bool(any(Rf['relaxed']))", 1),
 ('G04', "G's jobs include D(b)", "if r['tag'][0] in 'ac' and not r['skip']]", "if r['tag'][0] in 'abc' and not r['skip']]", 0),
 ('G05', "G's B jobs without the flow", "for i, r in enumerate(B) if not r['skip']]\n    j1 +=",
  "for i, r in enumerate(B) if not r['skip']][:0]\n    j1 +=", 0),
 ('G06', "η_V interpolating the amplitudes, not their signs", "np.r_[np.sign(As), np.zeros(len(As) * n)]", "np.r_[As, np.zeros(len(As) * n)]", 0),
 ('G07', "ν_V's margin over every starting point", "off = [j for j in range(len(c0)) if j not in e]", "off = [j for j in range(len(c0))]", 1),
 ('G08', "η_V's sup beyond 0.5σ", "axis=1) > EX_R * K.sig", "axis=1) > 0.5 * K.sig", 0),
 ('G09', "η_V's grid at the plane's spacing", "g = _grid(n, H_ETA[n - 1])", "g = _grid(n, H_ETA[-1])", 0),
 ('G10', "the Hessian without the feature's sign", "np.linalg.eigvalsh(np.sign(As[k]) * Q[k])", "np.linalg.eigvalsh(Q[k])", 0),
 ('G11', "nondegenerate without the Hessian", "nondeg=bool(sup < 1 and hmax < 0)", "nondeg=bool(sup < 1)", 0),
 ('G12', "basis pursuit drops the grid within 1.5 spacings", "axis=1) > 0.5 * H_BP[n - 1]]", "axis=1) > 1.5 * H_BP[n - 1]]", 0),
 ('G13', "basis pursuit's interpolation of the amplitudes", "b_eq=np.sign(As), bounds", "b_eq=As, bounds", 0),
 ('G14', "basis pursuit's basis truncated at 10⁻⁸", "Ub = U[:, sv > sv[0] * 1e-13]", "Ub = U[:, sv > sv[0] * 1e-8]", 0),
 ('G15', "OMP's selection unnormalised", "cc = np.abs(_corr(K, g, r)) / nrm", "cc = np.abs(_corr(K, g, r))", 0),
 ('G16', "OMP's refinement unnormalised",
  "obj = lambda c: -abs(float(K.phi(c[None, :])[:, 0] @ r)) / float(np.sqrt((K.phi(c[None, :])[:, 0] ** 2).sum()))",
  "obj = lambda c: -abs(float(K.phi(c[None, :])[:, 0] @ r))", 0),
 ('G17', "OMP's residual on the positions before the slide", "r = f - K.phi(sel) @ w\n", "r = f - Ph @ w\n", 0),
 ('G18', "the slide's signs all +", "RN.polish_block_nd(K, f, w, sel, np.sign(w))", "RN.polish_block_nd(K, f, w, sel, np.ones_like(w))", 0),
 ('G19', "OMP's distance in units of the samples", "axis=2), axis=1))) / K.sig", "axis=2), axis=1)))", 0),
 ('G20', "the plane's H jobs in one dimension", "(('C', i), r['cs'], r['As'], 2) for i, r in enumerate(C)", "(('C', i), r['cs'], r['As'], 1) for i, r in enumerate(C)", 0),
 ('G21', "G's 'same' from the recomputed path", "same = sum(x['out_flow'] == src[k[0]][k[1]]['out'] for k, x in rs)",
  "same = sum(x['out'] == src[k[0]][k[1]]['out'] for k, x in rs)", 0),
 ('G22', "G's total relaxed from 'same'", "tot = [tot[0] + len(rs), tot[1] + same, tot[2] + rel]", "tot = [tot[0] + len(rs), tot[1] + same, tot[2] + same]", 0),
 ('G23', "G's rejections with the timeouts' −1", "rej = [max(src[k[0]][k[1]].get('rej', 0), 0) for k in keys]", "rej = [src[k[0]][k[1]].get('rej', 0) for k in keys]", 0),
 ('G24', "H's rows with the timeouts", "for k in keys if k in hb and src[k[0]][k[1]]['out'] != 'timeout']", "for k in keys if k in hb]", 0),
 ('G25', "H's basis pursuit out of all targets", "n1 += row[0] if fam != 'C' else 0", "n1 += row[0]", 0),
 ('G26', "OMP's distance by its mean", "f\"{np.median(dm):.3g}σ, range", "f\"{np.mean(dm):.3g}σ, range", 0),
]
