# (id, what, old, new, k): single-point mutants of exp_gauge_noise.py — Part I, the early phase's gap of Proposition D‴.40 (round 147's
# open item): the batched loop with the gauge identity's decomposition (`i_runs`), the fixed-mass box term (`i_fixed`), its interpolation
# (`i_cbox`), the extended field and the box's projection (`i_field_row`, `grid_ext`), the SDE variants (`i_sde`), the start's transient
# (`i_transient`), the statistics (`_p2`, `holm`, `i_regressed`, `i_contrib`) and `report_I`.
# Run: venv/bin/python tools/mutants.py --exp exp_gauge_noise --tests tests/test_gauge_noise.py --raw results/chaos/gauge_noise_raw.npz
#      --md results/chaos/gauge_noise.md --out-env GN_OUT --workers-env GN_WORKERS --spec tools/mutants/gap40_gauge_noise_own.py
#      --work <scratch> -k "part_I or md"
# Result (one worker): 16 mutants and the control; all 16 killed by the tests (Q09–Q12 and Q16 also by the md). Q14 survived the first
# run (over 2·10³ slow units no amplitude reaches the box's edge, so the box changes nothing) and is killed by
# test_part_I_the_loop_without_the_box_leaves_it (seed 507 passes the edge by 5·10³).
M = [
 ('Q00', 'control', "def i_bd(box, gamma=GAMMA_E):", "def i_bd(box, gamma=GAMMA_E):", 0),
 ('Q01', "E_nφ without the Gaussian increment's factor", "Ephi = cphi * np.cosh(kk * (Th[:, :, 0] + DT * dT[:, :, 0])) * np.exp(0.5 * kk * kk * s2)",
  "Ephi = cphi * np.cosh(kk * (Th[:, :, 0] + DT * dT[:, :, 0]))", 0),
 ('Q02', "the burn-in on the seed's own generator", "g0 = [np.random.default_rng(s + BURN_OFF_I) for s in seeds]",
  "g0 = [np.random.default_rng(s) for s in seeds]", 0),
 ('Q03', "the box's term at the centres' half-width", "np.clip(Th[:, :, 0], -P.BOX[0], P.BOX[0])) / E.G[0])",
  "np.clip(Th[:, :, 0], -P.BOX[1], P.BOX[1])) / E.G[0])", 0),
 ('Q04', "the fixed-mass average started a step late", "        if k >= nb:", "        if k > nb:", 0),
 ('Q05', "the box term extrapolated upward below W_1's grid", "np.where(W1 < x[0], (W1 - x[0]) * lo, 0.0)",
  "np.where(W1 < x[0], (x[0] - W1) * lo, 0.0)", 0),
 ('Q06', "the box's projection with the force reversed", "f[3 * j] = -E.G[0] / gk[j]", "f[3 * j] = E.G[0] / gk[j]", 0),
 ('Q07', "the SDE's paths put on the wall, not reflected", "        for f in (2.0, 1.0):", "        for f in (1.0,):", 0),
 ('Q08', "the box's term given to class 2", "v[:, 5:11].reshape(-1, 2, 3)[:, :, 0]", "v[:, 5:11].reshape(-1, 2, 3)[:, :, 1]", 0),
 ('Q09', "Bonferroni for Holm", "run = max(run, min(1.0, (m - r) * p[i]))", "run = max(run, min(1.0, m * p[i]))", 0),
 ('Q10', "the p-value without the √2", "math.erfc(abs(z) / math.sqrt(2))", "math.erfc(abs(z))", 0),
 ('Q11', "the regressed mean's s.e. with n − 1", "r.std(ddof=2) / np.sqrt(n)", "r.std(ddof=1) / np.sqrt(n)", 0),
 ('Q12', "z against the sample's s.e. only", "z = [(m[0] - S[k]['rec'][t][:, 0].mean()) / zden(e[0], S[k]['rec'][t][:, 0].std(ddof=1), n) for k in (ref, soft)]",
  "z = [(m[0] - S[k]['rec'][t][:, 0].mean()) / e[0] for k in (ref, soft)]", 0),
 ('Q13', "the transient without the rate γ", "return JAC13 @ (GAMMA_E * dphi / beta)", "return JAC13 @ (dphi / beta)", 0),
 ('Q14', "the no-box runs with the box", "net, bd = i_bd(kind != 'no box')", "net, bd = i_bd(True)", 0),
 ('Q15', "the extended grid 0.010 above the saturations", "L_EXT_I = np.abs(E.A_BUMP) / A_MAX + 0.012", "L_EXT_I = np.abs(E.A_BUMP) / A_MAX + 0.010", 0),
 ('Q16', "the martingale's contribution with the wrong sign", "sg = np.array([1, -1, 1, 1, 1, 1])[None, :]", "sg = np.array([1, 1, 1, 1, 1, 1])[None, :]", 0),
]
