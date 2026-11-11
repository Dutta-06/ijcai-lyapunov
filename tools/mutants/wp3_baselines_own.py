"""Single-point mutants of src/core/baselines.py (WP3), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/baselines --tests tests/test_baselines.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp3_baselines_own.py --work <scratch dir> [--slot a --i0 0 --step 2] --extra experiments/exp_parametrization.py

(the tool's md / `--from-raw` columns mean nothing for these mutants: only `tests` does — True: tests/test_baselines.py kills the mutant). A tuple is (id, what, old, new[, k]): `old`
replaced by `new` at its k-th occurrence in src/core/baselines.py.

Result of the run (2026-11-08, one slot, 60 mutants, about 17 s each for the tests): 57 killed, 3 survivors, all in `iss_path`'s exit branch. I03 (the duals of active units also advance) changes a quantity that is
never read again once a unit has entered: equivalent. I05 (exit without the inward-correlation condition) and I06 (an exited kernel keeps its centre) change the behaviour only when a unit leaves the block, and
none of the planted targets of tests/test_baselines.py produces an exit: unreachable by the tests, not covered (the exit branch of the imaging path is covered by tests/test_imaging.py).
"""
M = [
    # refit_block / coef_lsq / entry_significance
    ('R01', 'refit_block: Jacobian without the sqrt(N) scale', "return -np.c_[K.phi(cc), (K.dphi(cc) * z[:k][None, :, None]).reshape(K.N, k * n)] / np.sqrt(K.N)", "return -np.c_[K.phi(cc), (K.dphi(cc) * z[:k][None, :, None]).reshape(K.N, k * n)]"),
    ('R02', 'refit_block: centre columns without the coefficient', "(K.dphi(cc) * z[:k][None, :, None]).reshape(K.N, k * n)", "(K.dphi(cc)).reshape(K.N, k * n)"),
    ('R03', 'entry_significance: noise level not applied', "return np.abs(g) * K.N / (s * np.linalg.norm(K.phi(c), axis=0))", "return np.abs(g) * K.N / (np.linalg.norm(K.phi(c), axis=0))"),
    # omp_slide
    ('O01', 'omp: signed correlation instead of absolute', "g = np.abs(K.corr(c_grid, r)[0])\n        g[used] = -1.0", "g = K.corr(c_grid, r)[0]\n        g[used] = -1.0"),
    ('O02', 'omp: a used kernel may enter again', "g[used] = -1.0\n        j = int(np.argmax(g))\n        if g[j] <= 0:", "j = int(np.argmax(g))\n        if g[j] <= 0:"),
    ('O03', 'omp: no slide', "w, c = refit_block(K, y, w, c)\n        r = y - K.synth(c, w)\n        blocks.append(dict(w=w, c=c, rss=float(r @ r)))\n    return dict(blocks=blocks, order=order, sig=np.array(sig))", "pass\n        r = y - K.synth(c, w)\n        blocks.append(dict(w=w, c=c, rss=float(r @ r)))\n    return dict(blocks=blocks, order=order, sig=np.array(sig))"),
    ('O04', 'omp: significance of the entrant after it entered', "sig.append(float(entry_significance(K, r, c_grid[j:j + 1], s)[0]))", "sig.append(float(entry_significance(K, y, c_grid[j:j + 1], s)[0]))"),
    ('O05', 'omp: the empty block has rss of the zero fit', "blocks = [dict(w=w, c=c, rss=float(r @ r))]\n    order, sig, used", "blocks = [dict(w=w, c=c, rss=0.0)]\n    order, sig, used"),
    # stepwise_bic
    ('T01', 'stepwise: best by largest rss', "if best is None or rss < best[0]:", "if best is None or rss > best[0]:"),
    ('T02', 'stepwise: chosen kernels may repeat', "if j in chosen:\n                continue", "if False:\n                continue"),
    # iss_path
    ('I01', 'iss: dual time with the wrong sign', "dt[ina] = (rho - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])", "dt[ina] = (rho + np.sign(v[ina]) * u[ina]) / np.abs(v[ina])"),
    ('I02', 'iss: ties do not enter together', "ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + 1e-300)", "ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + 1e-300)[:1]"),
    ('I03', 'iss: active duals also grow', "u[~act] += dmin * v[~act]", "u += dmin * v"),
    ('I04', 'iss: entrant sign from the coefficient, not the dual', "sg[j] = np.sign(v[j])", "sg[j] = 1.0"),
    ('I05', 'iss: exit without the inward-correlation condition', "out = (np.abs(wp) <= 1e-9 * max(np.abs(wp).max(), 1e-300)) & (sg[idx] * g <= vtol * vscale)", "out = (np.abs(wp) <= 1e-9 * max(np.abs(wp).max(), 1e-300))"),
    ('I06', 'iss: an exited kernel keeps its centre', "c[i] = c_grid[i]\n            u[i] = sg[i] * rho", "u[i] = sg[i] * rho"),
    ('I07', 'iss: entry times not accumulated', "tau += dmin\n", "tau = dmin\n"),
    # lasso
    ('A01', 'lasso_cd: soft threshold at lam / 2', "wj = np.sign(rho) * max(abs(rho) - lam, 0.0) / col2[j]", "wj = np.sign(rho) * max(abs(rho) - lam / 2, 0.0) / col2[j]"),
    ('A02', 'lasso_cd: no residual update', "r -= Phi[:, j] * (wj - w[j])", "pass"),
    ('A03', 'lasso_grid: lam_max below the first activation', "lam_max = float(np.abs(Phi.T @ y).max()) * (1.0 + 1e-9)", "lam_max = float(np.abs(Phi.T @ y).max()) * 0.5"),
    ('A04', 'merge_adjacent: every nonzero atom its own spike', "return int(1 + np.sum(np.diff(pos) > width))", "return int(idx.size)"),
    # nets
    ('N01', 'train_net: centres move against the gradient', "c = c + lc * (w[:, None] * g[1])", "c = c - lc * (w[:, None] * g[1])"),
    ('N02', 'train_net: no prox step', "if lam > 0 and prox:", "if False:"),
    ('N03', 'train_net: prox threshold without the step', "w = np.sign(w) * np.maximum(np.abs(w) - lw * lam, 0.0)", "w = np.sign(w) * np.maximum(np.abs(w) - lam, 0.0)"),
    ('N04', 'prune_net: keeps the small units', "kept = np.abs(w) >= thr * np.abs(w).max()", "kept = np.abs(w) <= thr * np.abs(w).max()"),
    # sfw
    ('F01', 'sfw: stops when the certificate is below 1/2', "if g[j] <= lam * (1.0 + tol):", "if g[j] <= lam * 0.5:"),
    ('F02', 'sfw: tolerance on the wrong side', "if g[j] <= lam * (1.0 + tol):", "if g[j] <= lam * (1.0 - tol):"),
    ('F03', 'sfw: lasso with the penalty not scaled by N', "w = lasso_cd(K.phi(c), y, lam * K.N, np.r_[w, 0.0])", "w = lasso_cd(K.phi(c), y, lam, np.r_[w, 0.0])"),
    ('F04', 'sfw: unit cap ignored', "if len(w) >= max_units:\n            break", "if False:\n            break"),
    ('F05', 'sfw: no slide', "w, c = _slide_l1(K, y, w[keep], c[keep], lam, lo, hi)\n            keep = np.abs(w)", "w, c = w[keep], c[keep]\n            keep = np.abs(w)"),
    ('F06', 'sfw: zero-weight units kept', "keep = np.abs(w) > 1e-9 * max(np.abs(w).max(), 1e-300)", "keep = np.abs(w) >= 0"),
    ('F07', 'sfw slide: penalty gradient with the wrong sign', "np.r_[lam - s * g0, (-(s * b)[:, None] * g1).ravel()]", "np.r_[-lam - s * g0, (-(s * b)[:, None] * g1).ravel()]"),
    ('F08', 'sfw slide: centre gradient with the wrong sign', "(-(s * b)[:, None] * g1).ravel()", "((s * b)[:, None] * g1).ravel()"),
    ('F09', 'sfw slide: coefficients may change sign', "bounds = [(0.0, None)] * k", "bounds = [(None, None)] * k"),
    ('F10', 'sfw slide: centres unbounded', "[(float(lo[j]), float(hi[j])) for _ in range(k) for j in range(n)]", "[(None, None) for _ in range(k) for j in range(n)]"),
    ('F11', 'sfw: the penalty of the objective without the absolute value', "lam * float(np.abs(w).sum()), cert=cert", "lam * float(w.sum()), cert=cert"),
    ('F12', 'sfw: certificate not divided by lam', "cert = float(np.abs(K.corr(c_grid, r)[0]).max() / lam)", "cert = float(np.abs(K.corr(c_grid, r)[0]).max())"),
    # cpgd
    ('C01', 'cpgd: potential with the wrong sign of the correlation', "V = lam - s * g0", "V = lam + s * g0"),
    ('C02', 'cpgd: multiplicative update against the potential', "w_new = np.maximum(w * np.exp(np.clip(-eta * V, -30.0, 30.0)), 1e-300)", "w_new = np.maximum(w * np.exp(np.clip(eta * V, -30.0, 30.0)), 1e-300)"),
    ('C03', 'cpgd: positions move against the gradient', "x_new = np.clip(x + eta * ratio_x * (s[:, None] * g1), lo, hi)", "x_new = np.clip(x - eta * ratio_x * (s[:, None] * g1), lo, hi)"),
    ('C04', 'cpgd: position gradient without the sign of the particle', "x_new = np.clip(x + eta * ratio_x * (s[:, None] * g1), lo, hi)", "x_new = np.clip(x + eta * ratio_x * g1, lo, hi)"),
    ('C05', 'cpgd: no backtracking acceptance test', "if obj_new <= obj:\n                break", "if True:\n                break"),
    ('C06', 'cpgd: step never recovers', "eta = min(eta * 1.05, 0.5 / G)", "eta = eta"),
    ('C07', 'cpgd: no negative particles', "s = np.r_[np.ones(d), -np.ones(d)]", "s = np.r_[np.ones(d), np.ones(d)]"),
    ('C08', 'cpgd: particles that did not grow count as spikes', "live = w > 10.0 * w_init", "live = w > 0.0"),
    ('C15', 'cpgd: dormant particles take part in the clustering', "lab = _cluster(xl, sl, link) if live.any() else np.zeros(0, int)", "lab = np.zeros(len(xl), int)"),
    ('C09', 'cpgd: clusters ignore the sign', "adj = (d <= radius) & (s[:, None] == s[None, :])", "adj = (d <= radius)"),
    ('C10', 'cpgd: spike mass is the largest particle', "W.append(sl[m][0] * wl[m].sum())", "W.append(sl[m][0] * wl[m].max())"),
    ('C11', 'cpgd: spike position is the plain mean', "C.append((xl[m] * wl[m][:, None]).sum(0) / wl[m].sum())", "C.append(xl[m].mean(0))"),
    ('C12', 'cpgd: weak spikes judged against the smallest', "big = np.abs(W) >= thr * np.abs(W).max()", "big = np.abs(W) >= thr * np.abs(W).min()"),
    ('C13', 'cpgd: link radius too small', "link = 1.5 * float(dd[dd > 0].min()) if (dd > 0).any() else 1.0", "link = 0.2 * float(dd[dd > 0].min()) if (dd > 0).any() else 1.0"),
    ('C14', 'cpgd: objective of the clustered measure without the penalty', "obj=0.5 * float(rr @ rr) / K.N + lam * float(np.abs(W).sum()), w_part=w", "obj=0.5 * float(rr @ rr) / K.N, w_part=w"),
    # rules
    ('U01', 'bic: N ln(rss) without the division by N', "return N * np.log(max(rss, 1e-300) / N) + n_params * np.log(N)", "return N * np.log(max(rss, 1e-300)) + n_params * np.log(N)"),
    ('U02', 'bic: no penalty', "return N * np.log(max(rss, 1e-300) / N) + n_params * np.log(N)", "return N * np.log(max(rss, 1e-300) / N)"),
    ('U03', 'select_bic: the largest value', "return int(np.argmin(vals))", "return int(np.argmax(vals))"),
    ('U04', 'discrepancy: the last block instead of the first within the noise', "if b['rss'] / N <= thr:\n            return k", "if b['rss'] / N <= thr:\n            return len(blocks) - 1"),
    ('U05', 'discrepancy: threshold without the factor', "thr = s * s * (1.0 + z * np.sqrt(2.0 / N))", "thr = s * s"),
    ('U06', 'significance: counts all entrants', "return int(bad[0]) if bad.size else len(sig)", "return len(sig)"),
    ('U07', 'gap: smallest gap', "return int(np.argmax(taus[k_min:] / np.maximum(taus[k_min - 1:-1], 1e-300))) + k_min", "return int(np.argmin(taus[k_min:] / np.maximum(taus[k_min - 1:-1], 1e-300))) + k_min"),
    ('U08', 'match_units: tolerance ignored', "rr, cc = linear_sum_assignment(np.where(d <= tol, d, big))\n    ok = d[rr, cc] <= tol", "rr, cc = linear_sum_assignment(np.where(d <= tol, d, big))\n    ok = d[rr, cc] <= 1e9"),
]
