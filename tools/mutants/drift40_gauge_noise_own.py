# (id, what, old, new, k): single-point mutants of exp_gauge_noise.py — Part M (D‴.40 after round 153, the missing O(γ²τ) drift):
# `m_limit`, `m_offset_lin`, `m_lim_K`, `m_loop`, `m_fixed`, `m_anti`, `m_ou`, `m_fd`, `m_loop_rates`, `report_M`
# (N18–N22: its comparisons and part (e)). N06, N07 mutate m_loop (the first occurrence), N08 the same line in m_fixed.
# Result (2026-11-05, -k part_M, one worker): 22 mutants and the control; all 22 killed by the tests (N04, N05, N12, N15, N17–N22 also
# change the md). N14 survived the first pass (the raw stores m_fd's Richardson value, and no test recomputed m_fd) and was killed by
# test_part_M_the_second_differences_and_richardson_on_a_known_map. See review/REVIEW_LOG.md, "D‴.40: the missing O(γ²τ) drift".
M = [
 ('N00', 'control', "def m_anti(a, tau):", "def m_anti(a, tau):", 0),
 ('N01', "the limit map at the step 4 throughout", "LIM_STEPS = ((1.0, 4000), (8.0, 4500))", "LIM_STEPS = ((4.0, 1000), (8.0, 4500))", 0),
 ('N02', "RK4's second stage at the full step", "k2 = f(Th + ds / 2 * k1[0], lr + ds / 2 * k1[1])", "k2 = f(Th + ds * k1[0], lr + ds / 2 * k1[1])", 0),
 ('N03', "the limit map's log-ratios without β_3", "lr = np.asarray(lb, float)[:, :2] - np.asarray(lb, float)[:, 2:3]", "lr = np.asarray(lb, float)[:, :2]", 0),
 ('N04', "the fit's motion with the wrong sign", "D = -h ** 2 / (E.G[0] * W)", "D = h ** 2 / (E.G[0] * W)", 0),
 ('N05', "the fit's motion without −W1ᵀ", "Mm = (gamma / 3) * (np.eye(3) - W[..., :, None]) * D[..., None, :]", "Mm = (gamma / 3) * np.eye(3) * D[..., None, :]", 0),
 ('N06', "m_loop: the half chain's second half step at z + w", "Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z - w)", "Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z + w)", 0),
 ('N07', "m_loop: the bridge's generator another seed", "gw = [np.random.default_rng(s + W_OFF_KH) for s in seeds]", "gw = [np.random.default_rng(s + 1 + W_OFF_KH) for s in seeds]", 0),
 ('N08', "m_fixed: the half chain's second half step at z + w", "Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z - w)", "Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z + w)", 1),
 ('N09', "m_loop: the half chain's second half step not renormalised", "lb[H_] = renorm(lb[H_] + h * db / b)", "lb[H_] = lb[H_] + h * db / b", 0),
 ('N10', "m_fixed: the burn-in counted in the half chain", "        if k >= nb:\n            add(Th[H_], dT, H_)", "        if True:\n            add(Th[H_], dT, H_)", 0),
 ('N11', "m_fixed: φ̇ without the chain rule's 2/a_max", "pd = cphi * kk * np.sinh(kk * T[:, :, 0]) * dT[:, :, 0]", "pd = cphi * np.sinh(kk * T[:, :, 0]) * dT[:, :, 0]", 0),
 ('N12', "m_anti summed over k", "return 0.5 * (c - c.swapaxes(1, 2)).sum(axis=1) / (9 * tau)", "return 0.5 * (c - c.swapaxes(1, 2)).sum(axis=2) / (9 * tau)", 0),
 ('N13', "m_ou: BC for CBᵀ", "X = (dphi[:, None] * dphi[None, :]) * (C @ Bm.T)[0::3, 0::3]", "X = (dphi[:, None] * dphi[None, :]) * (Bm @ C)[0::3, 0::3]", 0),
 ('N14', "m_fd: Richardson's weights", "rich=(4 * dr[1] - dr[0]) / 3", "rich=(2 * dr[1] - dr[0])", 0),
 ('N15', "m_loop_rates: no extrapolation", "rich=2 * rt[1] - rt[0]", "rich=rt[1]", 0),
 ('N16', "m_lim_K: the records from the next one", "idx = np.arange(n // 2, n, EVERY_M)", "idx = np.arange(n // 2 + 1, n, EVERY_M)", 0),
 ('N17', "report_M: z over twice the error", "{', '.join(f'{v:+.1f}' for v in (m_ - kz) / e_)}", "{', '.join(f'{v:+.1f}' for v in (m_ - kz) / e_ / 2)}", 0),
 ('N18', "report_M(e): the step effect not carried to the half's mean time", "ec = (q['full'] - q['rich']) * hs * GAMMA_E ** 2 * TAU", "ec = (q['full'] - q['rich']) * GAMMA_E ** 2 * TAU", 0),
 ('N19', "report_M: the z's errors added, not in quadrature", "(m1 - m2) / np.sqrt(e1 ** 2 + e2 ** 2)", "(m1 - m2) / (e1 + e2)", 0),
 ('N20', "report_M(c): the first transient without − 1/3", "2 * rr[g]['y0'][1] - rr[g]['y0'][0] - 1 / 3", "2 * rr[g]['y0'][1] - rr[g]['y0'][0]", 0),
 ('N21', "report_M(e): the rest's error without (c)'s", "kp - y2.mean(axis=0), np.sqrt(y2e ** 2 + ece ** 2))", "kp - y2.mean(axis=0), y2e)", 0),
 ('N22', "report_M: the γ dependence against the same γ", "zip(*ratio_se(rr[GAMMAS_M[1]]['rich'], rr[GAMMAS_M[0]]['rich']))", "zip(*ratio_se(rr[GAMMAS_M[1]]['rich'], rr[GAMMAS_M[1]]['rich']))", 0),
]
