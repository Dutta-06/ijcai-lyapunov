# (id, what, old, new, k): single-point mutants of exp_colour_selection.py — round 142's revision (referee D3ek), never mutation-tested
# when applied: Part K — K1 (`clopper_pearson`, `mcnemar_exact`, `k_stats`: the exact intervals, the paired discordant counts and the
# exact McNemar p-value), K2 (`vertex_limit`, `k_job`: the flow run to convergence on Part C's failures, the end's support, basic or not,
# the certificate for each τ_c and its τ_c → ∞ limit, R(end) − R(truth)) and the K section of `report`.
# Result: 25 mutants and a control. Against the tests as committed: 9 killed by the tests (K01–K05, K07–K10), 3 by the md (K22, K24, K25);
# 12 killed by the tests added then (test_part_K1_is_recomputed_from_parts_C_and_G, test_part_K2_k_job_reproduces_stored_ends,
# test_report_K_counts_planted_ends, test_k_job_on_a_planted_end — K17, K20 change no stored end: P_max ≤ 3·10⁻¹¹ and no end coefficient
# lies between 10⁻⁹ and 10⁻⁶ of its largest). One survivor, equivalent: K06 (Binomial(n, 1/2) is symmetric, so the two-sided exact test
# of max(b, c) has the p-value of min(b, c)).
M = [
 ('K00', 'control', "def mcnemar_exact(b, c):", "def mcnemar_exact(b, c):", 0),
 ('K01', "the interval's level 90%", "def clopper_pearson(x, n, level=0.95):", "def clopper_pearson(x, n, level=0.9):", 0),
 ('K02', "the lower end one-sided", "float(beta.ppf(a / 2, x, n - x + 1))", "float(beta.ppf(a, x, n - x + 1))", 0),
 ('K03', "the upper end's beta with n − x + 1", "float(beta.ppf(1 - a / 2, x + 1, n - x))", "float(beta.ppf(1 - a / 2, x + 1, n - x + 1))", 0),
 ('K04', "the lower end's beta with n − x", "float(beta.ppf(a / 2, x, n - x + 1))", "float(beta.ppf(a / 2, x, n - x))", 0),
 ('K05', "McNemar one-sided", "binomtest(min(b, c), b + c, 0.5).pvalue", "binomtest(min(b, c), b + c, 0.5, alternative='less').pvalue", 0),
 ('K06', "McNemar on max(b, c)", "binomtest(min(b, c), b + c, 0.5).pvalue", "binomtest(max(b, c), b + c, 0.5).pvalue", 0),
 ('K07', "McNemar's p 0 without discordant pairs", "return 1.0 if b + c == 0 else", "return 0.0 if b + c == 0 else", 0),
 ('K08', "the vertex limit's numerator Σα_ij·w_j", "np.abs((al / w[S][:, None]).sum(axis=0))", "np.abs((al * w[S][:, None]).sum(axis=0))", 0),
 ('K09', "the vertex limit's denominator over w_j (signed)", "(al ** 2 / np.abs(w[S])[:, None]).sum(axis=0)", "(al ** 2 / w[S][:, None]).sum(axis=0)", 0),
 ('K10', "the vertex limit's min over the atoms", "return float(np.max(np.abs((al / w[S][:, None])", "return float(np.min(np.abs((al / w[S][:, None])", 0),
 ('K11', "K1's pairs without the slow flow at τ_c = 1000", "for m in ('bp', 'sparse Bayesian learning (EM)', 'slow 1000'):",
  "for m in ('bp', 'sparse Bayesian learning (EM)'):", 0),
 ('K12', "K1's discordant 'slow only' as concordant", "nb, nc = int(np.sum(a & ~b_)), int(np.sum(~a & b_))", "nb, nc = int(np.sum(a & b_)), int(np.sum(~a & b_))", 0),
 ('K13', "K1's pairs against the slow flow at τ_c = 10", "a, b_ = np.asarray(out['slow 100']), np.asarray(out[m])", "a, b_ = np.asarray(out['slow 10']), np.asarray(out[m])", 0),
 ('K14', "K1's intervals at 90%", "ci=clopper_pearson(x, len(v))", "ci=clopper_pearson(x, len(v), 0.9)", 0),
 ('K15', "K1's SBL column with a G instance missing", "if (N, d) == (10, 40) and all((k, i) in gs for i in ii):", "if (N, d) == (10, 40) and any((k, i) in gs for i in ii):", 0),
 ('K16', "K2's basic ends below N atoms only", "basic = bool(len(S) <= Phi.shape[0] and", "basic = bool(len(S) < Phi.shape[0] and", 0),
 ('K17', "K2's relevances without the 4P²", "De = np.sqrt(w ** 2 + 4 * P ** 2)", "De = np.sqrt(w ** 2 + P ** 2)", 0),
 ('K18', "K2's R(truth) at τ_c = 100 for every τ_c", "CS.potential(Phi, np.abs(ws), 0, ('ou', tc)) for tc in TAUS_KR}",
  "CS.potential(Phi, np.abs(ws), 0, ('ou', 100.0)) for tc in TAUS_KR}", 0),
 ('K19', "K2's flow at τ_c = 1000", "S_SLOW(100.0), 0, ('ou', 100.0), S_max=S_MAX)", "S_SLOW(100.0), 0, ('ou', 1000.0), S_max=S_MAX)", 1),
 ('K20', "K2's end support at 10⁻⁶", "we = np.where(np.abs(w) > 1e-9 * np.abs(w).max(), w, 0.0)", "we = np.where(np.abs(w) > 1e-6 * np.abs(w).max(), w, 0.0)", 0),
 ('K21', "K2's N-atom ends any basic end", "nb = [r for r in rs if r['basic'] and r['nnz'] == 10]", "nb = [r for r in rs if r['basic']]", 0),
 ('K22', "K2's limit by its mean", "{np.median(lim):.3g} |", "{np.mean(lim):.3g} |", 0),
 ('K23', "K2's 'converged, failing' without convergence", "{sum(r['conv'] and not r['recovers'] for r in rs)} | {int(np.median(",
  "{sum(not r['recovers'] for r in rs)} | {int(np.median(", 0),
 ('K24', "R's comparison the other way", "np.mean([r['dR'][tc] < 0 for r in rs])", "np.mean([r['dR'][tc] > 0 for r in rs])", 0),
 ('K25', "R's median unscaled by τ_c", "np.median([tc * r['dR'][tc] for r in rs])", "np.median([r['dR'][tc] for r in rs])", 0),
]
