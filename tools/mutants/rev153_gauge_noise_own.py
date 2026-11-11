# (id, what, old, new, k): single-point mutants of exp_gauge_noise.py — round 153's revision (referee D3ev): Part J (`gauge_terms`,
# `j_runs`, `ratio_se`, `report_J`), Part K (`k_runs`, `k_half`, `k_identity`, `k_offset`, `k_response`, `k_forcing`, `k_effects`,
# `report_K`) and the report section L (`report_L`). R08, R10 mutate k_runs (the first occurrence), R23, R24 the same lines in k_half.
# Result: 24 mutants and a control (tools/mutants.py, one worker, -k 'part_J or part_K or part_L or regenerate'): the control passes; 22
# killed by the tests at the first pass (10 also by the md); R24 survived and was killed by a test comparing the half run's accumulators
# term by term with the raw; R09 is equivalent on the data: the box's term is identically zero in these runs (|α| ≤ 4.09 < L_α = 5).
M = [
 ('R00', 'control', "def ratio_se(x, y):", "def ratio_se(x, y):", 0),
 ('R01', "the half steps at the full step's noise", "rng = np.random.default_rng(seed); sq = np.sqrt(h)", "rng = np.random.default_rng(seed); sq = np.sqrt(DT)", 0),
 ('R02', "the full step driven by z₁ alone", "Th[F_] = Th[F_] + DT * dT[F_] + noise(z1 + z2)", "Th[F_] = Th[F_] + DT * dT[F_] + noise(z1)", 0),
 ('R03', "the Gaussian factor without its ½", "np.exp(0.5 * kk * kk * var2[sl] * hh)", "np.exp(kk * kk * var2[sl] * hh)", 0),
 ('R04', "the other step the own step", "alt2 = np.r_[np.full(B, h), np.full(B, DT)][:, None]", "alt2 = np.r_[np.full(B, DT), np.full(B, h)][:, None]", 0),
 ('R05', "the half chain's states counted once", "cnt = np.r_[np.full(B, nb), np.full(B, 2 * nb)][:, None]", "cnt = np.r_[np.full(B, nb), np.full(B, nb)][:, None]", 0),
 ('R06', "Lψ without the Itô term in Part J", "lpsi = u - g * nonfit + tk2[sl] * gp", "lpsi = u - g * nonfit", 0),
 ('R07', "the ratio's s.e. without the ratio", "return r, (x - r * y).std(axis=0, ddof=1)", "return r, (x - y).std(axis=0, ddof=1)", 0),
 ('R08', "the barrier given to the frozen variant", "fixman=np.array([{'nonlinear only': 0.0, 'nonlinear only + barrier': 1.0}[v]", "fixman=np.array([{'nonlinear only': 1.0, 'nonlinear only + barrier': 1.0}[v]", 0),
 ('R09', "the chain's term keeps the box's", "acc[:, 3] += u * DT - dphi - bx - co", "acc[:, 3] += u * DT - dphi - co", 0),
 ('R10', "the boundary term over the new masses", "dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / b;", "dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / np.exp(lb);", 0),
 ('R11', "the offset's W_k taken as 1/3", "return GAMMA_E / 3 * (d - W * d.sum(axis=-1, keepdims=True))", "return GAMMA_E / 3 * (d - d.sum(axis=-1, keepdims=True) / 3)", 0),
 ('R12', "the response's projection without −W·du", "GAMMA_E * (W * (d - (d @ W)[:, None]))[:, :2]", "GAMMA_E * (W * d)[:, :2]", 0),
 ('R13', "the response's Jacobian reversed", "Jc = np.array([(f(W + eps * e) - f(W - eps * e))[:2]", "Jc = np.array([(f(W - eps * e) - f(W + eps * e))[:2]", 0),
 ('R14', "the effects over the whole run", "r = K['var'][v]; s = K['s']; R = np.asarray(E_rec[v]); hi = R.shape[1] // 2\n    du", "r = K['var'][v]; s = K['s']; R = np.asarray(E_rec[v]); hi = 0\n    du", 0),
 ('R15', "the delta method without the white run's error", "** 2 + (xw.std(axis=0, ddof=1) / np.sqrt(n) / xw.mean(axis=0)) ** 2)", "** 2)", 0),
 ('R16', "b's estimation without n", "e2 = e * np.sqrt(1 + nn * m.mean() ** 2", "e2 = e * np.sqrt(1 + m.mean() ** 2", 0),
 ('R17', "the zero shift at γ = 10⁻³", "if z['variant'] == 'correlated E + class levels' and z['gamma'] == GAMMA_E][0]", "if z['variant'] == 'correlated E + class levels' and z['gamma'] == 1e-3][0]", 0),
 ('R18', "the half steps without the bridge", "for dz in (z + w, z - w):", "for dz in (z + w, z + w):", 0),
 ('R19', "the half steps' increments not halved", "Th = Th + h * dT + 0.5 * sq[None, None, :] * dz", "Th = Th + h * dT + sq[None, None, :] * dz", 0),
 ('R20', "the bridge's generator another seed", "gw = [np.random.default_rng(s + W_OFF_KH)", "gw = [np.random.default_rng(s + 1 + W_OFF_KH)", 0),
 ('R21', "the step's forcing not doubled", "f = 2 * (fv['acc'][:, :n, :4] - hv['acc'][:, :, :4]).sum(axis=2)", "f = (fv['acc'][:, :n, :4] - hv['acc'][:, :, :4]).sum(axis=2)", 0),
 ('R22', "the forcing with the masses form's change", "f = 2 * (fv['acc'][:, :n, :4] - hv['acc'][:, :, :4]).sum(axis=2)", "f = 2 * (fv['acc'][:, :n, :5] - hv['acc'][:, :, :5]).sum(axis=2)", 0),
 ('R23', "k_half: the barrier given to the frozen variant", "fixman=np.array([{'nonlinear only': 0.0, 'nonlinear only + barrier': 1.0}[v]", "fixman=np.array([{'nonlinear only': 1.0, 'nonlinear only + barrier': 1.0}[v]", 1),
 ('R24', "k_half: the boundary term over the new masses", "dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / b;", "dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / np.exp(lb);", 1),
]
