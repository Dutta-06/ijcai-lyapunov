# (id, what, old, new, k): single-point mutants of exp_gauge_noise.py — Part N (D‴.40 after Part M, the finite-τ question):
# `part_N`, `n_pool`, `n_fit`, `n_ratios`, `report_N`.
# Result (2026-11-05, -k part_N, one worker): 10 mutants and the control; all 10 killed by the tests (N01–N04, N07–N10 also change the md;
# N04 also breaks --from-raw). See review/REVIEW_LOG.md, "D‴.40: the finite-τ question (Part N)".
M = [
 ('N00', 'control', "def n_fit(tau, y, e):", "def n_fit(tau, y, e):", 0),
 ('N01', "n_fit: the slope's numerator with the wrong sign", "b = (S * Sxy - Sx * Sy) / D", "b = (S * Sxy + Sx * Sy) / D", 0),
 ('N02', "n_fit: the intercept's s.e. from S", "ec=np.sqrt(Sxx / D)", "ec=np.sqrt(S / D)", 0),
 ('N03', "n_fit: the line's degrees of freedom", "dof=len(x) - 2", "dof=len(x) - 1", 0),
 ('N04', "n_pool: Part M(d)'s seeds left out", "ch = {(0, c): v for c, v in ((M.get('fixed') or {}).get(t) or {}).items()}", "ch = {}", 0),
 ('N05', "part_N: a chunk short of its first seed", "fx[t][c0] = m_fixed((t, sd[c0:c0 + CHUNK_N], S_N[t], BURN_MF))",
  "fx[t][c0] = m_fixed((t, sd[c0 + 1:c0 + CHUNK_N], S_N[t], BURN_MF))", 0),
 ('N06', "part_N: one run length for every τ", "fx[t][c0] = m_fixed((t, sd[c0:c0 + CHUNK_N], S_N[t], BURN_MF))",
  "fx[t][c0] = m_fixed((t, sd[c0:c0 + CHUNK_N], 1e5, BURN_MF))", 0),
 ('N07', "n_ratios: the Δt chain in place of Richardson", "dict(r=v['rich'] / d", "dict(r=v['full'] / d", 0),
 ('N08', "report_N: the old-new difference's errors added, not in quadrature", "dz = z1((m2 - m1) / np.sqrt(e1 ** 2 + e2 ** 2))",
  "dz = z1((m2 - m1) / (e1 + e2))", 0),
 ('N09', "report_N: the seeds' spread with ddof 0", "q['pool']['r'].std(axis=0, ddof=1)", "q['pool']['r'].std(axis=0)", 0),
 ('N10', "report_N: the prediction's error without the covariance", "(tx * f['eb']) ** 2 + 2 * tx * f['cov'])", "(tx * f['eb']) ** 2)", 0),
]
