# (id, what, old, new, k): single-point mutants of exp_transit_certificate.py — Parts L and M after round 149 (the plane's two remaining
# computations): `_kron_basis` and `bp_certificate_kron` (basis pursuit's least certificate in the plane by an interior-point method with
# certified bounds), `KronGram` and `blasso_kron` (the BLASSO in the plane), and the L and M sections of `report`.
# Run with -k "kron or dual_lower or part_L or part_M or report_L" and --extra results/chaos/transit_certificate.md (the tests recount printed rows).
# Result: 19 mutants and a control (the control survives: 11 tests pass, the md regenerates). First run: 16 killed by the tests, R12 and
# R18 only by the md, R09 (the dual correction of the wrong sign) survived: the stationarity residual is at rounding level on the data, so
# the bound moved below the tests' 1e-12. The correction was made a function (`_dual_lower`, bit-identical results) and tested on an
# explicit program with a planted large residual; the report test was extended (the monotonicity count, the tolerance and radius). Second
# and third runs: all 19 killed by the tests (R12, R17, R18 also by the md). No survivor.
M = [
 ('R00', 'control', "def bp_certificate_kron(K, cs, As, h=None, tol=L_TOL, it_max=L_IT, cut=1e-13):",
  "def bp_certificate_kron(K, cs, As, h=None, tol=L_TOL, it_max=L_IT, cut=1e-13):", 0),
 ('R01', "the truncation on s_is_j against s_0, not s_0²", "np.outer(s1, s1) > s1[0] ** 2 * cut", "np.outer(s1, s1) > s1[0] * cut", 0),
 ('R02', "the features' grid neighbours kept within a quarter spacing only",
  "    keep = np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.5 * h\n    km = keep",
  "    keep = np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.25 * h\n    km = keep", 0),
 ('R03', "the lower bound without the mass's renormalisation", "return float(-(s @ nu) / (a.sum() + b.sum())), a, b",
  "return float(-(s @ nu)), a, b", 0),
 ('R04', "the upper bound over the dropped grid points too", "up = float(np.abs(eta[keep]).max())", "up = float(np.abs(eta).max())", 0),
 ('R05', "the corrector without the second-order term", "step(l1 * u1 + du1 * dl1 - sig * mu * km,", "step(l1 * u1 - sig * mu * km,", 0),
 ('R06', "the centring exponent 2", "/ (2 * n) / mu) ** 3", "/ (2 * n) / mu) ** 2", 0),
 ('R07', "the primal step to the boundary", "ap = 0.995 * min(amax(u1[keep], du1[keep])", "ap = 1.0 * min(amax(u1[keep], du1[keep])", 0),
 ('R08', "no stall test", "it - best['it'] >= L_STALL", "it - best['it'] >= L_STALL + 1000", 0),
 ('R09', "the dual correction of the wrong sign", "dlt = Bw(solve(-r)) * keep", "dlt = Bw(solve(r)) * keep", 0),
 ('R10', "the equalities' rows without the singular values", "np.exp(-(K.x1 - x) ** 2 / (2 * K.sig ** 2)) @ U1 / s1",
  "np.exp(-(K.x1 - x) ** 2 / (2 * K.sig ** 2)) @ U1", 0),
 ('R11', "the direct check over the dropped grid points too", "t_dir=float(np.abs(_corr(K, g[keep], p)).max())",
  "t_dir=float(np.abs(_corr(K, g, p)).max())", 0),
 ('R12', "the report's monotonicity reversed", "all(r['cuts'][i]['t'] <= r['cuts'][i + 1]['t'] + 1e-9", "all(r['cuts'][i]['t'] >= r['cuts'][i + 1]['t'] - 1e-9", 0),
 ('R13', "KronGram's product with the first factor twice", "np.einsum('am,ab,bm->m', self.F1, R, self.F2, optimize=True)",
  "np.einsum('am,ab,bm->m', self.F1, R, self.F1, optimize=True)", 0),
 ('R14', "KronGram's blocks without the mean", "(self.F2[:, I].T @ self.F2[:, J]) / self.N", "(self.F2[:, I].T @ self.F2[:, J])", 0),
 ('R15', "KronGram from the left with |z|", "        return self @ z\n", "        return self @ np.abs(z)\n", 0),
 ('R16', "the plane's BLASSO without warm starts", "x, kkt, it = lasso_fs(G, b, l * lmax, x)", "x, kkt, it = lasso_fs(G, b, l * lmax, None)", 1),
 ('R17', "the report's tolerance in the line's spacing", "tol, rd = I_TOL * H_BP[1], K_RD * K2.sig", "tol, rd = I_TOL * H_BP[0], K_RD * K2.sig", 0),
 ('R18', "the report's radius with the line's σ", "tol, rd = I_TOL * H_BP[1], K_RD * K2.sig", "tol, rd = I_TOL * H_BP[1], K_RD * 0.12", 0),
 ('R19', "λ_max from the features' columns only", "lmax = float(np.abs(b).max())", "lmax = float(np.abs(b[-len(As):]).max())", 1),
]
