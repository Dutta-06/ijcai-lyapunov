"""Single-point mutants of experiments/theory3f_escape.py (the optimizers, the closed forms and quadratures, the pair / duplicate / design runs, the hybrid,
the numpy KAN, the registries and the report), run with

    venv/bin/python tools/mutants.py --exp theory3f_escape --tests tests/test_theory3f_escape.py --raw results/synthetic/theory3f_raw.npz --md results/synthetic/theory3f.md \\
        --out-env T3F_OUT --spec tools/mutants/theory3f_escape_own.py --work <scratch>/mut --slot a

(id, what, old, new, k): `old` replaced by `new` at its k-th occurrence. F00 is the control (no change). The kill criterion is the tests file alone.

Run (105 entries, the control and 104 mutants, about 14 min): the control passes; 99 of the 104 mutants are killed by the tests (each also leaves `--from-raw` working). Five survivors, all equivalent:
F04 (ngd: a zero gradient gives 0/max(n, 1e-300) = 0, the same step); F14 (adam: a zero denominator needs v = 0, hence every past gradient 0 and m = 0: the step is 0 either way);
F19 (the bound's bias correction of m: the supremum C over t is attained at t = 37640, where 1 - b1^t = 1, so the value of C is unchanged); F56 (diag Euler: u = v along the balanced flow, so swapping the two updates is the same map);
F85 (report: the upper-bound half of the `inside` test: all 12 hybrid cells are inside the bounds, so the md is unchanged). The training loops of the design runs (`part_*`) are not mutated beyond the checks of
`train_opt` against `factorized_rbf.train_adam` and the regenerated `part_A`, `part_G`: their numbers are the experiment's output.
"""
M = [
    ('F00', 'control (no change)', '', ''),
    # ---- optimizers
    ('F01', 'gd: no learning rate', "            return p + self.lr * gneg\n        if k == 'sign':", "            return p + gneg\n        if k == 'sign':"),
    ('F02', 'sign: sign of the wrong gradient', "return p + self.lr * np.sign(gneg)", "return p - self.lr * np.sign(gneg)"),
    ('F03', 'ngd: L1 norm', "n = np.linalg.norm(gneg)", "n = np.abs(gneg).sum()"),
    ('F04', 'ngd: zero gradient moves', "return p if n == 0.0 else p + self.lr * gneg / n", "return p + self.lr * gneg / max(n, 1e-300)"),
    ('F05', 'pre applies to adam (per-coordinate)', "k in ('gd', 'adam_g', 'ngd')", "k in ('gd', 'adam', 'adam_g', 'ngd')"),
    ('F06', 'pre not applied to ngd', "k in ('gd', 'adam_g', 'ngd')", "k in ('gd', 'adam_g')"),
    ('F07', 'adam: no bias correction of m', "mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)", "mh, vh = self.m, self.v / (1 - self.b2 ** self.t)"),
    ('F08', 'adam: no bias correction of v', "mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)", "mh, vh = self.m / (1 - self.b1 ** self.t), self.v"),
    ('F09', 'adam: eps inside the root', "den = np.sqrt(vh) + self.ea", "den = np.sqrt(vh + self.ea)"),
    ('F10', 'adam: m uses b2', "self.m = self.b1 * self.m + (1 - self.b1) * gneg", "self.m = self.b2 * self.m + (1 - self.b2) * gneg"),
    ('F11', 'adam: v without the square', "gg = gneg ** 2 if k == 'adam' else", "gg = np.abs(gneg) if k == 'adam' else"),
    ('F12', 'adam_g: per-coordinate v', "gg = gneg ** 2 if k == 'adam' else np.mean(gneg ** 2) * np.ones_like(gneg)", "gg = gneg ** 2"),
    ('F13', 'adam_g: sum instead of mean', "np.mean(gneg ** 2) * np.ones_like(gneg)", "np.sum(gneg ** 2) * np.ones_like(gneg)"),
    ('F14', 'adam: zero denominator moves', "np.divide(mh, den, out=np.zeros_like(den), where=den > 0)", "mh / np.maximum(den, 1e-300)"),
    ('F15', 'default eps_adam', "beta1=0.9, beta2=0.999, eps_adam=1e-8, pre=None):\n        assert", "beta1=0.9, beta2=0.999, eps_adam=1e-7, pre=None):\n        assert"),
    # ---- the ratio bounds
    ('F16', 'bounds: c = 1 - b2', "return 1.0 - b1, float(np.sqrt(C2[k])), int(t[k])", "return 1.0 - b2, float(np.sqrt(C2[k])), int(t[k])"),
    ('F17', 'bounds: geometric ratio b1/b2', "gam = b1 * b1 / b2", "gam = b1 / b2"),
    ('F18', 'bounds: (1 - b1) factor', "C2 = (1 - b1) ** 2 * (1 - b2 ** t)", "C2 = (1 - b1) * (1 - b2 ** t)"),
    ('F19', 'bounds: bias correction of m', "((1 - b2) * (1 - b1 ** t) ** 2)", "((1 - b2) * (1 - b1 ** t))"),
    ('F20', 'bounds: tmax', "def adam_ratio_bounds(b1=0.9, b2=0.999, tmax=10 ** 6):", "def adam_ratio_bounds(b1=0.9, b2=0.999, tmax=10 ** 3):"),
    # ---- the D3 units and the runs
    ('F21', 'D3: w = u^3 + v^3', "return P[:, 0] ** 3 - P[:, 1] ** 3", "return P[:, 0] ** 3 + P[:, 1] ** 3"),
    ('F22', 'D3: direction of v', "-3 * P[:, 1] ** 2 * g], 1)", "3 * P[:, 1] ** 2 * g], 1)"),
    ('F23', 'D3: slope 2u^2', "np.stack([3 * P[:, 0] ** 2 * g,", "np.stack([2 * P[:, 0] ** 2 * g,"),
    ('F24', 'run_units: drive without the Gram', "g = c - G @ d3_w(P)", "g = c"),
    ('F25', 'run_units: drive sign', "g = c - G @ d3_w(P)", "g = (G @ d3_w(P)) - c"),
    ('F26', 'run_units: threshold strict', "new = np.isnan(T) & (np.abs(W[t]) >= thr)", "new = np.isnan(T) & (np.abs(W[t]) > 2 * thr)"),
    ('F27', 'run_units: threshold on w, not |w|', "new = np.isnan(T) & (np.abs(W[t]) >= thr)", "new = np.isnan(T) & (W[t] >= thr)"),
    ('F28', 'run_units: start v = 0', "* 2, 1).copy()", "* 2, 1).copy() * np.array([1.0, 0.0])"),
    ('F29', 'sign_escape_step: a plain cube', "while 2.0 * (3 * s * s * lr * t + (lr * t) ** 3) < thw:", "while 2.0 * ((lr * t) ** 3) < thw:"),
    ('F30', 'sign_escape_step: factor 2 dropped', "while 2.0 * (3 * s * s * lr * t + (lr * t) ** 3) < thw:", "while (3 * s * s * lr * t + (lr * t) ** 3) < thw:"),
    ('F31', 'gd_escape_time: no 3', "return (1.0 / s - 1.0 / u_c) / (3.0 * g)", "return (1.0 / s - 1.0 / u_c) / g"),
    ('F32', 'gd_escape_time: 1/s only', "return (1.0 / s - 1.0 / u_c) / (3.0 * g)", "return (1.0 / s) / (3.0 * g)"),
    ('F33', 'frozen run: gradient not 3u^2', "gr = scale * 3.0 * u * u", "gr = scale * 3.0 * u"),
    ('F34', 'frozen run: no v square', "v = beta2 * v + (1 - beta2) * gr * gr", "v = beta2 * v + (1 - beta2) * gr"),
    ('F35', 'frozen run: eps ignored', "np.sqrt(v / (1 - beta2 ** t)) + ea)", "np.sqrt(v / (1 - beta2 ** t)) + 0.0)"),
    ('F36', 'frozen run: threshold >', "if u >= u_c:", "if u > 2 * u_c:"),
    ('F37', 'pair gram: c = (1, 1)', "np.array([1.0, rho])      # G", "np.array([1.0, 1.0])      # G"),
    ('F38', 'pair gram: off-diagonal 0', "return np.array([[1.0, rho], [rho, 1.0]]),", "return np.array([[1.0, 0.0], [0.0, 1.0]]),"),
    ('F39', 'pair: peak of w_1', "float(np.nanmax(np.abs(W[:, 1])))", "float(np.nanmax(np.abs(W[:, 0])))"),
    ('F40', 'pair: T1 and T2 swapped', "return T[0], T[1], float", "return T[1], T[0], float"),
    ('F41', 'duplicates: Gram identity', "G, c = np.ones((2, 2)), np.ones(2)", "G, c = np.eye(2), np.ones(2)"),
    ('F42', 'duplicates: same start', "np.array([s1, ratio * s1]), n_steps, 1e-3)", "np.array([s1, s1]), n_steps, 1e-3)"),
    # ---- the design runs
    ('F43', 'train_opt: not the unit rate', "FR.Unit(readout, np.ones_like(readout.lr))", "readout"),
    ('F44', 'train_opt: centres not moved', "z = opt.step(z, base(0.0, z))", "z = opt.step(z, base(0.0, z) * np.r_[np.ones(readout.q * d), np.zeros(d * n)])"),
    ('F45', 'train_opt: clock', "taus = rec / n_steps", "taus = rec / (2 * n_steps)"),
    ('F46', 'formation_steps: clock', "return FR.formation_times(res['taus'], res['W'], thr) * n_steps", "return FR.formation_times(res['taus'], res['W'], thr)"),
    # ---- chains and the diagonal network
    ('F47', 'chain quad: z^D', "1.0 / (z ** (D - 1) * (g - z ** D))", "1.0 / (z ** D * (g - z ** D))"),
    ('F48', 'chain quad: without g - z^D', "1.0 / (z ** (D - 1) * (g - z ** D))", "1.0 / (z ** (D - 1) * g)"),
    ('F49', 'chain leading: D - 1', "s ** (-(D - 2)) / ((D - 2) * g)", "s ** (-(D - 2)) / ((D - 1) * g)"),
    ('F50', 'chain leading: log branch', "return np.log(1.0 / s) / g if D == 2 else", "return np.log(1.0 / s) / g if D == 3 else"),
    ('F51', 'chain ODE: wrong power', "[z[0] ** (D - 1) * (g - z[0] ** D)], [0, 1e12]", "[z[0] ** D * (g - z[0] ** D)], [0, 1e12]"),
    ('F52', 'diag closed form: rate 1', "np.exp(-2.0 * g * np.asarray(t, float))", "np.exp(-1.0 * g * np.asarray(t, float))"),
    ('F53', 'diag closed form: g/w0', "(g / w0 - 1.0) * np.exp", "(g / w0 + 1.0) * np.exp"),
    ('F54', 'diag time: frac odds', "(g / w0 - 1.0) * frac / (1.0 - frac)", "(g / w0 - 1.0) * frac"),
    ('F55', 'diag time: 2g', "/ (2.0 * g))\n\n\ndef diag_gd_check", "/ g)\n\n\ndef diag_gd_check"),
    ('F56', 'diag Euler: gradient of u', "u, v = u + h * v * r, v + h * u * r", "u, v = u + h * u * r, v + h * v * r"),
    # ---- the paper-unit pair
    ('F57', 'pair rhs: alpha gradient', "return [n1 ** 2 * w2 * G, 2 * n1 * a1 * w2 * G, n2 ** 2 * w1 * G, 2 * n2 * a2 * w1 * G]", "return [n1 ** 2 * w2 * G, 2 * n1 * a1 * w2 * G, n2 ** 2 * w1 * G, n2 * a2 * w1 * G]"),
    ('F58', 'pair rhs: drive', "G = cstar - w1 * w2", "G = cstar - w1"),
    ('F59', 'pair time: half', "z[2] * z[3] ** 2 - 0.5 * cs", "z[2] * z[3] ** 2 - 0.25 * cs"),
    ('F60', 'pair sym: conserved nu', "f = lambda a: 1.0 / ((s * s + 2 * a * a) ** 2 * a *", "f = lambda a: 1.0 / ((s * s + a * a) ** 2 * a *"),
    ('F61', 'pair sym: end point', "(s * s + 2 * a * a) ** 2 * a * a - 0.5 * cstar, a0, 10.0)", "(s * s + 2 * a * a) ** 2 * a * a - 0.25 * cstar, a0, 10.0)"),
    # ---- the hybrid
    ('F62', 'hybrid: conservation constant', "kappa2 = 1.0 / 3.0\n\n    def rhs(t, z):", "kappa2 = 1.0 / 2.0\n\n    def rhs(t, z):"),
    ('F63', 'hybrid: second layer drive', "return np.r_[(3 * u ** 2 * b[:, None] * gam[None, :]).ravel(), w @ gam]", "return np.r_[(3 * u ** 2 * b[:, None] * gam[None, :]).ravel(), 2 * (w @ gam)]"),
    ('F64', 'hybrid: first layer slope', "(3 * u ** 2 * b[:, None] * gam[None, :]).ravel()", "(2 * u ** 2 * b[:, None] * gam[None, :]).ravel()"),
    ('F65', 'hybrid: stop on the mean', "ev = lambda t, z: z[:Kh * Lk].max() - U", "ev = lambda t, z: z[:Kh * Lk].mean() - U"),
    ('F66', 'hybrid: leader', "lead = int(np.argmax((Phi.T @ fstar / N)))", "lead = int(np.argmin((Phi.T @ fstar / N)))"),
    ('F67', 'hybrid bounds: lower b', "lo = quad(lambda u: 1.0 / (3 * (1 + delta) * g0 * u ** 2 * np.sqrt(b0 ** 2 + kappa2 * n_pos * (u ** 2 - s ** 2))), s, U", "lo = quad(lambda u: 1.0 / (3 * (1 + delta) * g0 * u ** 2 * np.sqrt(b0 ** 2)), s, U"),
    ('F68', 'hybrid bounds: delta side', "up = (1.0 / s - 1.0 / U) / (3 * (1 - delta) * g0 * b0)", "up = (1.0 / s - 1.0 / U) / (3 * (1 + delta) * g0 * b0)"),
    ('F69', 'hybrid: dictionary width', "N_X, N_L, SIG_X, KH_H = 81, 21, 0.12, 2", "N_X, N_L, SIG_X, KH_H = 81, 21, 0.10, 2"),
    # ---- the numpy KAN
    ('F70', 'kan: no standardisation floor', "np.maximum(h.std(1, keepdims=True), 0.05)", "np.maximum(h.std(1, keepdims=True), 1e-9)"),
    ('F71', 'kan: dtanh', "gW1 * P['n1'] ** 2 / np.cosh(P['a1']) ** 2", "gW1 * P['n1'] ** 2"),
    ('F72', 'kan: dnu', "n1=gW1 * 2 * P['n1'] * np.tanh(P['a1']))", "n1=gW1 * P['n1'] * np.tanh(P['a1']))"),
    ('F73', 'kan: hybrid second-layer rate', "l2=0.3, ea=1e-8", "l2=0.0, ea=1e-8"),
    ('F74', 'kan: bias rate', "bias += 0.1 * lr * gb", "bias += 0.0 * lr * gb"),
    ('F75', 'kan: second layer starts at ones', "P.update(w2=0.1 * self.rng.standard_normal((K, M2)))", "P.update(w2=np.zeros((K, M2)))"),
    ('F76', 'kan: nu ratio', "nu_ratio=0.05, l2=0.3", "nu_ratio=0.5, l2=0.3"),
    # ---- registries that reach the report
    ('F77', 'registry: drives of A', "G_A = (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)", "G_A = (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 2e-5)"),
    ('F78', 'registry: D set', "D_G, S_G = (2, 3, 4, 5, 6), (0.4, 0.3, 0.2, 0.1, 0.05)", "D_G, S_G = (2, 3, 4, 5, 7), (0.4, 0.3, 0.2, 0.1, 0.05)"),
    ('F79', 'registry: rho set', "RHO_C, S_C, THW_C = (0.2, 0.5, 0.8)", "RHO_C, S_C, THW_C = (0.2, 0.5, 0.9)"),
    ('F80', 'registry: duplicate ratio (report text)', "DUP_RATIO = 1.31e-2, 0.7 ", "DUP_RATIO = 1.31e-2, 0.8 "),
    ('F81', 'registry: step budget (report text)', "NSTEP_A = (1.3e-3, 1e-5, 1e-7), 1e-4, 1e-6, 400", "NSTEP_A = (1.3e-3, 1e-5, 1e-7), 1e-4, 1e-6, 500"),
    ('F82', 'registry: STEPS_F in the text', "STEPS_F, H_F, EPS_F = 200000, 0.05, 1e-6", "STEPS_F, H_F, EPS_F = 100000, 0.05, 1e-6"),
    # ---- report
    ('F83', 'report: peak prediction', "_f(r / (1 + r)) for r in RHO_C", "_f(r / (1 + 2 * r)) for r in RHO_C"),
    ('F84', 'report: window of sign descent', "_f(float(Ts[i].max() / Ts[i].min()))", "_f(float(Ts[i].max() - Ts[i].min()))"),
    ('F85', 'report: inside test', "'yes' if r[3] <= r[2] <= r[4] else 'no'", "'yes' if r[3] <= r[2] else 'no'"),
    ('F86', 'report: fitted exponent window', "np.polyfit(np.log(ss[-3:]), np.log(tt[-3:]), 1)[0]", "np.polyfit(np.log(ss[:3]), np.log(tt[:3]), 1)[0]"),
    ('F87', 'report: formed count', "int(ok.sum())", "int(ok.sum()) + 1"),
    # ---- the additive KAN and the report of K, D2, H
    ('F88', 'addkan: width clamp', "s = np.exp(np.clip(p['ls'], np.log(0.05), np.log(3.0)))\n        u =", "s = np.exp(np.clip(p['ls'], np.log(0.05), np.log(2.0)))\n        u ="),
    ('F89', 'addkan: clamp mask', "inside = (p['ls'] > np.log(0.05)) & (p['ls'] < np.log(3.0))", "inside = True"),
    ('F90', 'addkan: centre gradient without w', "c=(r[:, None, None] * phi * u / s[None]).mean(0) * w,", "c=(r[:, None, None] * phi * u / s[None]).mean(0),"),
    ('F91', 'addkan: amplitude cap', "self.A = 3.0", "self.A = 2.0"),
    ('F92', 'addkan: bias gradient', "b=np.array([r.mean()]))", "b=np.array([0.0]))"),
    ('F93', 'fit_add: arming level', "if v < 0.9 * v0:\n                armed = True\n            if armed and since >= 2000:", "if v < 0.99 * v0:\n                armed = True\n            if armed and since >= 2000:"),
    ('F94', 'fit_add: last iterate instead of the best', "    model.p = best_p\n    return t", "    return t"),
    ('F95', 'fit_add: validation every 25 -> 50', "if t % 25 == 0:\n            v = float(np.mean((model.fwd(Xva)[0] - yva) ** 2))", "if t % 50 == 0:\n            v = float(np.mean((model.fwd(Xva)[0] - yva) ** 2))"),
    ('F96', 'report K: over-count', "np.maximum(cnt - truth, 0).sum(1).mean()", "(cnt - truth).sum(1).mean()"),
    ('F97', 'report K: exact share', "exact = (cnt == truth)", "exact = (cnt >= truth)"),
    ('F98', 'report D2: crossing index', "idx = int(bad.max()) + 1 if bad.size else 0", "idx = int(bad.max()) if bad.size else 0"),
    ('F99', 'report H: formed count', "int(r[0]), 'yes' if r[1] else 'no', 'yes' if r[2] else 'no', _f(float(r[3]))) for e, r in zip(R['H_ea'], Hh)", "int(r[0]) + 1, 'yes' if r[1] else 'no', 'yes' if r[2] else 'no', _f(float(r[3]))) for e, r in zip(R['H_ea'], Hh)"),
    ('F100', 'registry: grid of D2', "tuple(np.geomspace(0.02, 8.0, 12))", "tuple(np.geomspace(0.02, 9.0, 12))"),
    ('F101', 'registry: EA_F', "EA_F = (1e-8, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 1e-1, 1.0, 10.0)", "EA_F = (1e-8, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 1e-1, 1.0, 20.0)"),
    ('F102', 'registry: targets of K', "TARGETS_K, NOISES_K, SEED_K, EPS_K, M_K, STEPS_K = (101, 102, 103)", "TARGETS_K, NOISES_K, SEED_K, EPS_K, M_K, STEPS_K = (101, 102, 104)"),
    ('F103', 'report E: last/first', "_f(float(np.nanmax(T) / np.nanmin(T)))", "_f(float(np.nanmax(T) - np.nanmin(T)))"),
    ('F104', 'report: Euler check', "'Euler check of the closed form (drive, steps*h, formula): %s.'", "'Euler check of the closed form (drive, steps*h, formula) %s.'"),
]
