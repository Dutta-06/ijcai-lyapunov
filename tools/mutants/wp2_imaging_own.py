"""Single-point mutants of src/core/imaging.py (WP2), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/imaging --tests tests/test_imaging.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp2_imaging_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(the tool's md / `--from-raw` columns mean nothing for these mutants: only `tests` does — True: tests/test_imaging.py kills the mutant). A tuple is (id, what, old, new[, k]): `old`
replaced by `new` at its k-th occurrence in src/core/imaging.py.

Result (tests/test_imaging.py, 61 tests): 52 of the 58 mutants are killed. The first run left 18 survivors; tests were added for the ones that change the result (P01, P10, B02, S03, S14, T02: smoother support,
reach with a high-pass, start clipped to the box, ties entering together, the record ratio of the gap-jump stop, the degrees of freedom of the discrepancy rule; S17, S19, S20, T11, T12, F02: the fit radius, the
local window against the global fit in noise [noise-free data cannot tell a window from a larger one: the model fits exactly in any window], the window reach of the widest kernel, the refine box and its
reweighting, the stopping-time rule of `fit_frame` against the significance rule). The six that remain do not change any output the code can produce:
  P12  equivalent: in the branch where the window is cut from a root PSF, `root` is `self`, so `self.y0 - root.y0` is 0 and the offset that the mutant drops is zero;
  S04  equivalent: `v` is already 0 where `blocked` is true, so `U += dmin * v` and `U[~blocked] += dmin * v[~blocked]` are the same update;
  S05  equivalent: the dual of a kernel that has entered is never read while the kernel is taken (it is blocked), and it is reset when the kernel exits (S08's line);
  S06  redundant under the optimality conditions of the non-negative fit: a unit at zero amplitude after a converged block fit has a non-positive correlation, so the second condition never decides;
  S07, S08  concern the re-entry of a unit that has exited: in 345 random frames (2 to 3 emitters, noise-free and noisy) with 25 exits no exited unit entered again, so no frame is known that separates
        the mutants from the code; the behaviour (the exited kernel stays at the wall and may re-enter) is the algorithm's definition, documented in `sliding_path`.
"""
M = [
    ('P01', 'smoother truncated at 2 sb', "(np.abs(d[:, None] - d[None, :]) <= 4 * sb)", "(np.abs(d[:, None] - d[None, :]) <= 2 * sb)"),
    ('P02', 'smoother normalised by columns', "return S / S.sum(1, keepdims=True)", "return S / S.sum(0, keepdims=True)"),
    ('P03', 'pixel integral with a sum of erfs', "0.5 * (erf(up) - erf(um))", "0.5 * (erf(up) + erf(um))"),
    ('P04', 'sign of the centre derivative', "-(eu - em) / (_SQ2PI * sg)", "(eu - em) / (_SQ2PI * sg)"),
    ('P05', 'sign of the width derivative', "out += (-(up * eu - um * em) / (np.sqrt(np.pi) * sg),)", "out += ((up * eu - um * em) / (np.sqrt(np.pi) * sg),)"),
    ('P06', 'upper pixel edge at +1/4', "up = (p[:, None] + 0.5 - c[None, :])", "up = (p[:, None] + 0.25 - c[None, :])"),
    ('P07', 'high-pass partner with the wrong sign', "terms.append((-a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))", "terms.append((a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))"),
    ('P08', 'part widths not scaled by kappa (x)', "ax = self.axis(px, c[:, 0], kap * s0, deriv_s)", "ax = self.axis(px, c[:, 0], s0, deriv_s)"),
    ('P09', 'width derivative without kappa', "fx = [ax[0], ax[1]] + ([kap * ax[2]] if deriv_s else [])", "fx = [ax[0], ax[1]] + ([ax[2]] if deriv_s else [])"),
    ('P10', 'reach: max instead of hypot', "return float(np.hypot(r, self.hp) if self.hp else r)", "return float(max(r, self.hp) if self.hp else r)"),
    ('P11', 'high-pass without the x smoother', "return y if not self.hp else y - root.Sy @ y @ root.Sx.T", "return y if not self.hp else y - root.Sy @ y"),
    ('P12', 'window crop offset of the y axis lost', "w._crop = (ya + (self.y0 - root.y0), yb + (self.y0 - root.y0),", "w._crop = (ya, yb,"),
    ('P13', 'synth without the part coefficient', "return sum(a * ((gy * w) @ gx.T) for a, gy, gx, *_ in self.comps(c))", "return sum(((gy * w) @ gx.T) for a, gy, gx, *_ in self.comps(c))"),
    ('P14', 'corr: x gradient with the kernel itself', "vx + a * (A * dgx.T).sum(1)", "vx + a * (A * gx.T).sum(1)"),
    ('P15', 'jac: x column without the coefficient', "cols = [cw, cx * W, cy * W]", "cols = [cw, cx, cy * W]"),
    ('P16', 'jac: width column without the coefficient', "([cs * W] if wid else [])", "([cs] if wid else [])"),
    ('P17', 'phi without the part coefficient', "return sum(a * (gy[:, None, :] * gx[None, :, :]).reshape(self.N, -1) for a, gy, gx, *_ in self.comps(c))", "return sum((gy[:, None, :] * gx[None, :, :]).reshape(self.N, -1) for a, gy, gx, *_ in self.comps(c))"),
    ('L01', 'lattice margin shifted', "self.xs = np.arange(-margin,", "self.xs = np.arange(-margin + 0.1,"),
    ('L02', 'lattice correlation of the first term only', "return sum(a * (gy.T @ Rw @ gx) for a, gy, gx in self.terms)", "return sum(a * (gy.T @ Rw @ gx) for a, gy, gx in self.terms[:1])"),
    ('L03', 'near mask stamp too small', "reach = int(np.ceil(r / self.h)) + 1", "reach = int(np.ceil(r / self.h)) - 1"),
    ('L04', 'norm2 without cross terms', "out = out + a * a2 * ((gy * gy2).T @ W @ (gx * gx2))", "out = out + a * a * ((gy * gy).T @ W @ (gx * gx))"),
    ('B01', 'block fit allows negative amplitudes', "lo = np.r_[np.zeros(k), c_lo.T.ravel()]", "lo = np.r_[np.full(k, -np.inf), c_lo.T.ravel()]"),
    ('B02', 'start not clipped to the box', "z0 = np.r_[w1, np.clip(c, c_lo, c_hi).T.ravel()]", "z0 = np.r_[w1, c.T.ravel()]"),
    ('B03', 'only the coefficients are fitted (stage 1) when joint', "if not joint:\n        return c, w1, dict(nfev=0, status=0)", "if True:\n        return c, w1, dict(nfev=0, status=0)"),
    ('S01', 'negative correlations grow too', "grow = v > vtol * vscale", "grow = np.abs(v) > vtol * vscale"),
    ('S02', 'dual time with the wrong sign', "dt[grow] = (rho - U[grow]) / v[grow]", "dt[grow] = (rho + U[grow]) / v[grow]"),
    ('S03', 'ties do not enter together', "ent = np.argwhere(dt <= dmin * (1.0 + tie_tol) + tie_tol * 1e-300)", "ent = np.argwhere(dt == dmin)[:1]"),
    ('S04', 'blocked duals also grow', "U[~blocked] += dmin * v[~blocked]", "U += dmin * v"),
    ('S05', 'entrant dual not put on the wall', "            U[b, a] = rho\n            taken[b, a] = True", "            U[b, a] = 0.0\n            taken[b, a] = True"),
    ('S06', 'exit needs only the coefficient at zero', "gone = (w <= 1e-12 * max(w.max(), 1e-300)) & (g <= vtol * vscale)", "gone = (w <= 1e-12 * max(w.max(), 1e-300))"),
    ('S07', 'an exited kernel stays taken', "            taken[ids[i]] = False\n", "            pass\n"),
    ('S08', 'an exited kernel is reset to zero dual', "            U[ids[i]] = rho\n            del ids[i]", "            U[ids[i]] = 0.0\n            del ids[i]"),
    ('S09', 'tau_max compares the time before the step', "if tau_max is not None and tau + dmin > tau_max:", "if tau_max is not None and tau > tau_max:"),
    ('S10', 'greedy step count', "            tau += 1.0\n", "            tau += 2.0\n"),
    ('S11', 'the resolution limit is ignored', "blocked = taken | lat.near(c[:, :2], r_ex)", "blocked = taken"),
    ('S12', 'significance with unsquared weights', "snrs.append(float(v[b, a] / np.sqrt(np.sum(swv.ravel() ** 2 * ph_ ** 2))))", "snrs.append(float(v[b, a] / np.sqrt(np.sum(swv.ravel() * ph_ ** 2))))"),
    ('S13', 'snr early stop one entry early', "stop_snr_at = n_ent + 1 + snr_extra", "stop_snr_at = n_ent + snr_extra"),
    ('S14', 'record ratio inverted', "if taus and tau / taus[-1] > best[0]:", "if taus and tau / taus[-1] < best[0]:"),
    ('S15', 'only the entrant is refitted', "            S = _cluster(c, seeds, dint)", "            S = np.array(seeds)"),
    ('S16', 'the others are not subtracted in the window', "oth = np.setdiff1d(np.arange(len(w)), S)", "oth = np.array([], int)"),
    ('S17', 'fit radius halved', "< fit_radius)", "< 0.5 * fit_radius)"),
    ('S18', 'chi2 is a sum', "return float(np.mean((swv * (y - model)) ** 2))", "return float(np.sum((swv * (y - model)) ** 2))"),
    ('S19', 'window radius halved', "                r = rtrunc * reach\n", "                r = 0.5 * rtrunc * reach\n"),
    ('S20', 'reach ignores the width bound', "reach = psf.reach(psf.sigma if sig_bounds is None else sig_bounds[1])", "reach = psf.reach(psf.sigma)"),
    ('T01', 'discrepancy band below', "(1.0 + z * np.sqrt(2.0 / dof))", "(1.0 - z * np.sqrt(2.0 / dof))"),
    ('T02', 'discrepancy degrees of freedom', "dof = max(N - 3 * k, 1)", "dof = max(N - 2 * k, 1)"),
    ('T03', 'gap rule on differences', "ratio = taus[k_min:] / np.maximum(taus[k_min - 1:-1], 1e-300)", "ratio = taus[k_min:] - taus[k_min - 1:-1]"),
    ('T04', 'gate inverted', "np.asarray(snr)[0] < gate", "np.asarray(snr)[0] > gate"),
    ('T05', 'significance rule inverted', "bad = np.flatnonzero(snr < kappa)", "bad = np.flatnonzero(snr > kappa)"),
    ('T06', 'stopping time keeps the late entrants', "K = int(np.sum(np.asarray(taus, float) <= tau_stop))", "K = int(np.sum(np.asarray(taus, float) >= tau_stop))"),
    ('T07', 'noise scale by the mean', "return float(np.sqrt(np.median(lat.norm2(W))))", "return float(np.sqrt(np.mean(lat.norm2(W))))"),
    ('T08', 'PSF stamp not normalised', "    stamp = np.asarray(stamp, float)\n    stamp = stamp / stamp.sum()\n    n = stamp.shape[0]", "    stamp = np.asarray(stamp, float)\n    n = stamp.shape[0]"),
    ('T09', 'EM excess factor 1', "return (2.0 / cam['qe']) * np.asarray(mean_photons, float)", "return (1.0 / cam['qe']) * np.asarray(mean_photons, float)"),
    ('T10', 'deterministic EM gain', "out = np.where(n > 0, rng.gamma(np.maximum(n, 1), cam['em_gain']), 0.0)", "out = np.where(n > 0, n * cam['em_gain'], 0.0)"),
    ('T11', 'refine box doubled', "lo, hi = c - box, c + box\n        if c.shape[1] == 3:", "lo, hi = c - 2 * box, c + 2 * box\n        if c.shape[1] == 3:"),
    ('T12', 'refine without reweighting', "        if sw_fn is not None:\n            swv = np.asarray(sw_fn(psf.synth(c, w)) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((psf.ny, psf.nx))\n    return c, w, swv", "    return c, w, swv"),
    ('F01', 'fit_frame ignores the high-pass of the data', "    y = psf.highpass(y0 - np.median(y0))\n    sw0 = sw_fn(np.zeros((0, nc)), np.zeros(0))", "    y = y0 - np.median(y0)\n    sw0 = sw_fn(np.zeros((0, nc)), np.zeros(0))"),
    ('F02', 'fit_frame time rule uses the significance', "        i = stop_time(R['hist'], R['taus'], 1.0 / (par * noise_scale(lat, sw0 ** 2)))[0]", "        i = stop_snr(R['hist'], R['taus'], R['snr'], par)[0]"),
]
