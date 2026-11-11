# (id, what, old, new, k): single-point mutants of exp_transit_certificate.py — the runs after round 140 (D3ei-5): Part I (`bp_grid`,
# `lasso_fs`: the LASSO by feature-sign search with its KKT check; `blasso`: the BLASSO on Part H's basis-pursuit grid at λ = ℓλ_max;
# `blasso_verdict`; `i_job`, `i_jobs`), Part J (`jets`, `bp_jet_certificate`: basis pursuit certified by a p in the features' second-order
# jets; `j_job`, `j_jobs`), the resumable `run_chunks`, and the I and J sections of `report`.
# Run with -k "lasso or blasso or bp_grid or parts_I or run_chunks or i_job or part_I or jets or part_J or report_I" and
# --extra results/chaos/transit_certificate.md (the tests recount printed rows from the raw).
# Result: 23 mutants and a control (the control survives: 10 tests pass, the md regenerates). 19 killed by the tests, 2 by the md alone
# (P11, P21: counts printed from the raw); P05 and P07 survived the first run (the tests' stored targets had no full step that flips a
# sign, and compared no iteration counts) and are killed by test_i_job_reproduces_stored_blasso as extended (B 13, B 83, D 21; the
# iterations compared). P13 leaves the jets' span unchanged (φ is in it), so t_J is the same; test_jets_and_the_jet_certificate kills it on
# the columns. No survivor.
M = [
 ('P00', 'control', "def bp_grid(K, cs, h):", "def bp_grid(K, cs, h):", 0),
 ('P01', "the grid's neighbours dropped within 1.5 spacings", "axis=1) > 0.5 * h]", "axis=1) > 1.5 * h]", 0),
 ('P02', "feature-sign stops 1% short of λ", "if abs(gr[i]) <= lam * (1 + tol):", "if abs(gr[i]) <= 1.01 * lam:", 0),
 ('P03', "the added coordinate's sign reversed", "th[i] = -np.sign(gr[i])", "th[i] = np.sign(gr[i])", 0),
 ('P04', "the line search without the zero crossings", "ts = [1.0] + [t for t in tc if 0 < t < 1]", "ts = [1.0]", 0),
 ('P05', "a full step trusted whatever its signs", "full = bt == 1 and bool(np.all(np.sign(x[S]) == th[S]))", "full = bt == 1", 0),
 ('P06', "the support's KKT residual not relative to λ", "float(np.abs(c[S] - lam * np.sign(x[S])).max() / lam) if len(S) else 0.0)",
  "float(np.abs(c[S] - lam * np.sign(x[S])).max()) if len(S) else 0.0)", 0),
 ('P07', "no warm start", "x, kkt, it = lasso_fs(G, b, l * lmax, x)", "x, kkt, it = lasso_fs(G, b, l * lmax)", 0),
 ('P08', "the correlations not averaged over the samples", "G, b = P.T @ P / K.N, P.T @ f / K.N", "G, b = P.T @ P / K.N, P.T @ f", 0),
 ('P09', "recovery's tolerance doubled for the far atoms", "far = int(np.sum(d.min(1) > tol))", "far = int(np.sum(d.min(1) > 2 * tol))", 0),
 ('P10', "the clusters' signs ignored", "signs = all(c > 0 and np.sign(a) == np.sign(A) for c, a, A in zip(cnt, amp, As))",
  "signs = all(c > 0 for c in cnt)", 0),
 ('P11', "the support exact within the tolerance", "exact = len(x) == len(As) and bool(np.all(d.min(1) < 1e-12))",
  "exact = len(x) == len(As) and bool(np.all(d.min(1) <= tol))", 0),
 ('P12', "the amplitude error absolute", "err = float(max(abs(a - A) / abs(A) for a, A in zip(amp, As)))",
  "err = float(max(abs(a - A) for a, A in zip(amp, As)))", 0),
 ('P13', "the second derivatives without −δ/σ² (the same span with φ)", "cols += [P * (dx[:, :, a] * dx[:, :, b] - (a == b) / K.sig ** 2)",
  "cols += [P * (dx[:, :, a] * dx[:, :, b])", 0),
 ('P14', "first-order jets only", "    return np.concatenate(cols, 1)\n", "    return np.concatenate(cols[:1 + n], 1)\n", 0),
 ('P15', "the jets' basis truncated at 10⁻²", "Ub = U[:, sv > sv[0] * 1e-13]\n    Qg, Qs", "Ub = U[:, sv > sv[0] * 1e-2]\n    Qg, Qs", 0),
 ('P16', "the jet certificate interpolates the amplitudes", "A_eq=np.c_[Qs, np.zeros(len(As))], b_eq=np.sign(As), bounds",
  "A_eq=np.c_[Qs, np.zeros(len(As))], b_eq=As, bounds", 0),
 ('P17', "η_V's sup without the modulus", "t_etav=float(np.abs(Pg.T @ q).max())", "t_etav=float((Pg.T @ q).max())", 0),
 ('P18', "the plane's J jobs in one dimension", "(('C', i), r['cs'], r['As'], 2) for i, r in enumerate(C) if not r['skip']]\n    return jobs",
  "(('C', i), r['cs'], r['As'], 1) for i, r in enumerate(C) if not r['skip']]\n    return jobs", 0),
 ('P19', "no save after a chunk", "        save()\n    order", "        pass\n    order", 0),
 ('P20', "the done rows recomputed", "todo = [j for j in jobs if tuple(j[0]) not in have]", "todo = list(jobs)", 0),
 ('P21', "η_V's split over every target", "nd = [k for k in rows1 if k in hb and hb[k]['nondeg']]", "nd = [k for k in rows1 if k in hb]", 0),
 ('P22', "J's agreement as the count below 1", "{sum((x['t_jet'] < 1) == (hb[k]['t_bp'] < 1) for k, x in one)}",
  "{sum(x['t_jet'] < 1 for k, x in one)}", 0),
 ('P23', "I's recovery split by the path 'not' for 'other'",
  "{sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'other')} \"", "{sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'not')} \"", 0),
]
