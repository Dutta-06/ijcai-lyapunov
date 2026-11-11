"""Single-point mutants of src/core/imaging_fast.py (the fast evaluation of the sliding path), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/imaging_fast --tests tests/test_imaging_fast.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp_imagingfast_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(the tool's md / `--from-raw` columns mean nothing for these mutants: only `tests` does, True: tests/test_imaging_fast.py kills the mutant). A tuple is (id, what, old, new[, k]): `old`
replaced by `new` at its k-th occurrence (0-based) in src/core/imaging_fast.py. The mutants are the kernel factors (edges, signs, widths, the high-pass partners, the crop), the Jacobian's blocks, the Gram
iteration's lines (scalings, trust-region quantities, rank test, termination), the block fit, the lattice correlation, the early exit, the incremental model and the exits.

Result (tests/test_imaging_fast.py, 74 tests, 113 mutants; run with OMP_NUM_THREADS=1): 97 are killed, 16 survive. A first run killed 83 and left 30; the tests of the
non-square frame, the step-for-step agreement with scipy, the zero column, the unit initial radius, the step tolerance, the evaluation counts, the box, the second pass and the tie were added and
killed 14 of them. The 16 survivors, none of which any test can distinguish from the code:
  EQUIVALENT (the mutated line is redundant or its effect is below the rounding):
  T11  |lambda| for max(lambda, 0): only eigenvalues at the rounding noise of a positive semi-definite matrix are negative;
  T19  `actual_reduction >= 0`: a step with exactly zero reduction;  T24 `g_norm <= gtol`: exact equality;  P03 `snrs[-1] <= stop_kappa`, P25 `tau + dmin >= tau_max`, P28 `snrs[-1] <= snr_min`: the
       same, equality of two floating-point numbers that are computed independently;
  B03  no clip of the start into the box: `trf_gram` makes the start feasible itself (scipy's `least_squares` would raise on a start outside the box, the clip exists in `imaging.block_fit` for that);
  B08  the cached factors used for any point: `jacT` is only ever called at the point of the last `fun` call (the accepted step), so the guard is a safety check;
  P11  `U[ids[i]] = rho` at an exit: the dual of an entered point is already rho and blocked points do not grow, so the line is redundant (it is in `imaging.sliding_path` too);
  P13  the incremental model also when the fit is global: `model - pre + post` is `post` then;
  P19  an exited unit's image is not subtracted from the carried model: an exit needs w <= 1e-12 of the largest coefficient, the image is below 1e-12 of the model's peak (the test worlds' planted exits have w = 0);
  P24  no refresh of the carried model from the units (every 64 steps): it only removes the drift of the additions (about 1e-13 of the peak after 1300 steps);
  P27  `hist` keeps the live arrays: the arrays are replaced (vstack/delete) at each entry, so the stored blocks are never modified afterwards.
  NOT DISTINGUISHED, effect on the iteration count not shown by any tested problem:
  T10  eigenvalues in ascending order: only the rank test reads which end is largest; with the rank-deficient (zero-column) test the steps agree for 40 evaluations either way;
  T17  the Levenberg-Marquardt parameter is reset instead of carried: it only seeds the root search of the trust-region sub-problem, which stops at 1 % of the radius;
  T27  theta fixed at 0.995: differs only for gradients below 0.005 near a bound.
"""
M = [
    # ---- the kernel factors
    ('A01', 'lower edge of the first pixel shifted', "edges[:-1] = p - 0.5", "edges[:-1] = p + 0.5"),
    ('A02', 'upper edge of the last pixel', "edges[-1] = p[-1] + 0.5", "edges[-1] = p[-1] - 0.5"),
    ('A03', 'the argument is not divided by sqrt 2 sigma', "u = (edges[:, None] - c[None, :]) / (_SQ2 * sg[None, :])", "u = (edges[:, None] - c[None, :]) / sg[None, :]"),
    ('A04', 'exp(-u) for exp(-u^2)', "e = np.exp(-u ** 2)", "e = np.exp(-u)"),
    ('A05', 'the pixel integral without the one half', "out = (0.5 * (E[1:] - E[:-1]),", "out = ((E[1:] - E[:-1]),"),
    ('A06', 'the derivative without its minus sign', "-(e[1:] - e[:-1]) / (_SQ2PI * sg[None, :])", "(e[1:] - e[:-1]) / (_SQ2PI * sg[None, :])"),
    ('A07', 'the width derivative without the factor u', "ue = u * e", "ue = e"),
    ('A08', 'the width derivative not divided by sigma', "(_SQPI * sg[None, :])", "_SQPI"),
    ('K01', 'the partner terms have the same sign', "(np.array([1.0, -1.0] * len(psf.parts)) if psf.hp else 1.0)", "(np.array([1.0, 1.0] * len(psf.parts)) if psf.hp else 1.0)"),
    ('K02', 'no partner terms in the coefficients', "for _ in range(2 if psf.hp else 1)]", "for _ in range(1)]"),
    ('K03', 'a kernel width is the PSF width even when the kernel has its own', "s0 = c[:, 2] if c.shape[1] > 2 else np.full(k, psf.sigma)", "s0 = np.full(k, psf.sigma)"),
    ('K04', 'part width multiplier dropped (x)', "ax = _axis(self.px, c[:, 0], kap * s0, deriv_s, deriv)", "ax = _axis(self.px, c[:, 0], s0, deriv_s, deriv)"),
    ('K05', 'part width multiplier dropped (y)', "ay = _axis(self.py, c[:, 1], kap * s0, deriv_s, deriv)", "ay = _axis(self.py, c[:, 1], s0, deriv_s, deriv)"),
    ('K06', 'the width derivative is not scaled by the part multiplier', "f = kap if j == 2 else 1.0", "f = 1.0"),
    ('K07', 'no width derivative', "nf = (3 if deriv_s else 2) if deriv else 1", "nf = 2 if deriv else 1"),
    ('K08', 'the high-pass of x applied to the rows', "sx = (root.Sx @ np.concatenate(ax, axis=1))[xa:xb]", "sx = (root.Sy @ np.concatenate(ax, axis=1))[xa:xb]"),
    ('K09', 'the y crop applied to x', "sx = (root.Sx @ np.concatenate(ax, axis=1))[xa:xb]", "sx = (root.Sx @ np.concatenate(ax, axis=1))[ya:yb]"),
    ('K10', 'the partner of y not smoothed by the y smoother', "sy = (root.Sy @ np.concatenate(ay, axis=1))[ya:yb]", "sy = (root.Sx @ np.concatenate(ay, axis=1))[ya:yb]"),
    ('K11', 'the crop of a window is ignored', "self.crop = psf._crop if (psf.hp and psf._parent is not None) else None", "self.crop = None"),
    ('K12', 'the pixel axes of a high-passed window are the window\'s own', "self.py, self.px = (root.py, root.px) if psf.hp else (psf.py, psf.px)", "self.py, self.px = (psf.py, psf.px)"),
    ('K13', 'the products do not carry the coefficients', "A = (self.coef[:, None, None] * FY).transpose(2, 1, 0)", "A = FY.transpose(2, 1, 0)"),
    ('K14', 'the image is not weighted by the coefficients w', "* np.asarray(w, float)[None, None, :]).transpose(1, 0, 2)", "* 1.0).transpose(1, 0, 2)"),
    ('J01', 'the x derivative column uses the y derivative factor', "JT[k:2 * k] = self._mm(GY, DGX)", "JT[k:2 * k] = self._mm(DGY, GX)"),
    ('J02', 'the y derivative column uses the x derivative factor', "JT[2 * k:3 * k] = self._mm(DGY, GX)", "JT[2 * k:3 * k] = self._mm(GY, DGX)"),
    ('J03', 'the width column without its y term', "JT[3 * k:] = self._mm(f[4], GX) + self._mm(GY, f[5])", "JT[3 * k:] = self._mm(GY, f[5])"),
    ('J04', 'the width column without its x term', "JT[3 * k:] = self._mm(f[4], GX) + self._mm(GY, f[5])", "JT[3 * k:] = self._mm(f[4], GX)"),
    ('J05', 'the derivative columns not multiplied by the amplitudes', "JT[k:] *= np.tile(w, 3 if wid else 2)[:, None]", "JT[k:] *= 1.0"),
    ('J06', 'the Jacobian has the sign of the model, not of the residual', "JT *= -np.asarray(sw, float).reshape(1, -1)", "JT *= np.asarray(sw, float).reshape(1, -1)"),
    ('J07', 'the Jacobian is not weighted', "JT *= -np.asarray(sw, float).reshape(1, -1)", "JT *= -1.0"),
    ('J08', 'the amplitudes are tiled twice with a width', "np.tile(w, 3 if wid else 2)", "np.tile(w, 2)"),
    ('C01', 'the correlation without the part coefficients', "a = self.coef[t]", "a = 1.0"),
    ('C02', 'the x gradient of the correlation uses the value factor', "vx + a * (A * DGX[t].T).sum(1)", "vx + a * (A * GX[t].T).sum(1)"),
    ('C03', 'the y gradient of the correlation uses the value derivative of x', "vy + a * (B * GX[t].T).sum(1)", "vy + a * (B * DGX[t].T).sum(1)"),
    # ---- the lattice correlation
    ('L01', 'the transposed high-pass is the high-pass', "Rh = Rw - root.Sy.T @ Rw @ root.Sx if self.psf.hp else Rw", "Rh = Rw - root.Sy @ Rw @ root.Sx.T if self.psf.hp else Rw"),
    ('L02', 'no high-pass on the residual', "Rh = Rw - root.Sy.T @ Rw @ root.Sx if self.psf.hp else Rw", "Rh = Rw"),
    ('L03', 'the parts are indexed as if there were no partners', "step = 2 if psf.hp else 1", "step = 1"),
    ('L04', 'a window\'s lattice uses the fast formula', "self._ok = (psf._parent is None)", "self._ok = True"),
    ('L05', 'the lattice correlation without the part coefficients', "out = out + a * (B[m * self._nly:(m + 1) * self._nly] @ gx)", "out = out + (B[m * self._nly:(m + 1) * self._nly] @ gx)"),
    ('L06', 'the stacked parts are cut at the wrong row', "B[m * self._nly:(m + 1) * self._nly]", "B[m * self._nly:(m + 1) * self._nly + 1]"),
    # ---- the iteration
    ('T01', 'a zero column keeps its zero scale', "scale_inv[scale_inv == 0] = 1", "scale_inv[scale_inv == 0] = 0"),
    ('T02', 'the initial radius without the Coleman-Li scaling', "Delta = np.linalg.norm(x * scale_inv / v ** 0.5)", "Delta = np.linalg.norm(x * scale_inv)"),
    ('T03', 'a zero initial radius becomes 5', "    if Delta == 0:\n        Delta = 1.0", "    if Delta == 0:\n        Delta = 5.0"),
    ('T04', 'the default evaluation cap', "max_nfev = x.size * 100", "max_nfev = x.size"),
    ('T05', 'the evaluation cap is not enforced exactly', "if status is not None or nfev == max_nfev:", "if status is not None or nfev > max_nfev:"),
    ('T06', 'the scaling vector is not rescaled in the loop', "        v[dv != 0] *= scale_inv[dv != 0]\n        d = v ** 0.5 * scale", "        d = v ** 0.5 * scale"),
    ('T07', 'd without the square root of the Coleman-Li vector', "d = v ** 0.5 * scale", "d = v * scale"),
    ('T08', 'the bound regularisation without the scale', "diag_h = g * dv * scale", "diag_h = g * dv"),
    ('T09', 'the bound regularisation not added to the Gram matrix', "        Ah[np.diag_indices(n)] += diag_h\n", "        pass\n"),
    ('T10', 'the eigenvalues in ascending order', "lam, V = lam[::-1], V[:, ::-1]", "lam, V = lam, V"),
    ('T11', 'negative eigenvalues made positive', "s = np.sqrt(np.maximum(lam, 0.0))", "s = np.sqrt(np.abs(lam))"),
    ('T12', 'the projected gradient uses the unscaled gradient', "suf = V.T @ g_h", "suf = V.T @ g"),
    ('T13', 'the rank test accepts everything', "m_rank = RANK_RATIO / _EPS", "m_rank = 1.0"),
    ('T14', 'the step is not made strictly feasible', "x_new = make_strictly_feasible(x + step, lb, ub, rstep=0)", "x_new = x + step"),
    ('T15', 'the cost reduction with the wrong sign', "actual_reduction = cost - cost_new", "actual_reduction = cost_new - cost"),
    ('T16', 'the bound-hit test of the radius update', "step_h_norm > 0.95 * Delta", "step_h_norm > 0.5 * Delta"),
    ('T17', 'the Levenberg-Marquardt parameter is not carried', "            alpha *= Delta / Delta_new\n", "            alpha = 0.0\n"),
    ('T18', 'the column scales may shrink', "scale_inv = np.maximum(np.sqrt(np.diag(A)), scale_inv)", "scale_inv = np.sqrt(np.diag(A))"),
    ('T19', 'a step with no reduction is accepted', "        if actual_reduction > 0:\n            x, f, cost = x_new, f_new, cost_new", "        if actual_reduction >= 0:\n            x, f, cost = x_new, f_new, cost_new"),
    ('T20', 'the status is lost', "return x, nfev, (0 if status is None else status)", "return x, nfev, 0"),
    ('T21', 'the initial cost without the half', "    cost = 0.5 * np.dot(f, f)\n    g = JT @ f\n    A = JT @ JT.T\n    scale_inv", "    cost = np.dot(f, f)\n    g = JT @ f\n    A = JT @ JT.T\n    scale_inv"),
    ('T22', 'the Gram matrix is not refreshed after an accepted step', "            g = JT @ f\n            A = JT @ JT.T\n            scale_inv", "            g = JT @ f\n            scale_inv"),
    ('T23', 'the gradient is not refreshed after an accepted step', "            JT = jacT(x)\n            g = JT @ f\n", "            JT = jacT(x)\n"),
    ('T24', 'gradient tolerance test is not strict', "if g_norm < gtol:", "if g_norm <= gtol:"),
    ('T25', 'the step quality uses the wrong scaled norm', "np.linalg.norm(x), ratio, ftol, xtol)", "0.0, ratio, ftol, xtol)"),
    ('T26', 'the non-finite residual is accepted', "if not np.all(np.isfinite(f_new)):", "if False:"),
    ('T27', 'theta fixed', "theta = max(0.995, 1 - g_norm)", "theta = 0.995"),
    # ---- the block fit
    ('B01', 'the NNLS start is unweighted', "w1, _ = nnls((K.phiT(c) * sv).T, yv, maxiter=200 * k)", "w1, _ = nnls(K.phiT(c).T, yv, maxiter=200 * k)"),
    ('B02', 'the joint fit is skipped when asked for', "    if not joint:\n        return c, w1, dict(nfev=0, status=0)\n    lo = ", "    if joint:\n        return c, w1, dict(nfev=0, status=0)\n    lo = "),
    ('B03', 'the start is not clipped to the box', "z0 = np.r_[w1, np.clip(c, c_lo, c_hi).T.ravel()]", "z0 = np.r_[w1, c.T.ravel()]"),
    ('B04', 'the amplitudes may be negative', "lo = np.r_[np.zeros(k), c_lo.T.ravel()]", "lo = np.r_[np.full(k, -np.inf), c_lo.T.ravel()]"),
    ('B05', 'the amplitudes are bounded above', "hi = np.r_[np.full(k, np.inf), c_hi.T.ravel()]", "hi = np.r_[np.full(k, 1e3), c_hi.T.ravel()]"),
    ('B06', 'the fall-back catches the wrong error', "except np.linalg.LinAlgError:", "except KeyError:"),
    ('B07', 'the fall-back ignores the options', "return IM.block_fit(psf, y, sw, c, w, c_lo, c_hi, joint=joint, max_nfev=max_nfev, tol=tol)", "return IM.block_fit(psf, y, sw, c, w, c_lo, c_hi)"),
    ('B08', 'the cached factors are used for any point', "f = memo['f'] if ('z' in memo and np.array_equal(memo['z'], z)) else None", "f = memo['f'] if 'z' in memo else None"),
    ('B09', 'the residual has the wrong sign', "return (swv * (y - K.synth(cc, z[:k], f))).ravel()", "return (swv * (K.synth(cc, z[:k], f) - y)).ravel()"),
    ('B10', 'the evaluation cap is not passed on', "trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, max_nfev)", "trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, None)"),
    ('B11', 'the tolerance is not passed on (ftol)', "trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, max_nfev)", "trf_gram(fun, jacT, z0, lo, hi, 1e-3, tol, tol, max_nfev)"),
    ('B12', 'the tolerance is not passed on (gtol)', "trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, max_nfev)", "trf_gram(fun, jacT, z0, lo, hi, tol, tol, 1e-3, max_nfev)"),
    ('B13', 'an empty block is fitted', "    if k == 0:\n        return c, w, dict(nfev=0, status=0)\n    K = kern", "    if k == -1:\n        return c, w, dict(nfev=0, status=0)\n    K = kern"),
    ('R01', 'refine\'s box', "lo, hi = c - box, c + box", "lo, hi = c - 2 * box, c + 2 * box"),
    ('R02', 'refine runs one round', "rounds=2, tol=1e-10", "rounds=1, tol=1e-10"),
    ('R03', 'refine does not reweight', "        if sw_fn is not None:\n            swv = np.asarray(sw_fn(K.synth(c, w)) if sw_on == 'model'", "        if sw_fn is None:\n            swv = np.asarray(sw_fn(K.synth(c, w)) if sw_on == 'model'"),
    ('R04', 'refine reweights with the model when asked for the parameters', "sw_fn(K.synth(c, w)) if sw_on == 'model' else sw_fn(c, w)", "sw_fn(K.synth(c, w)) if sw_on != 'model' else sw_fn(c, w)"),
    # ---- the path
    ('P01', 'the early exit needs a strictly larger threshold', "stop_kappa >= snr_min", "stop_kappa > snr_min"),
    ('P02', 'the early exit works without snr_min', "early = stop_kappa is not None and snr_min is not None and stop_kappa >= snr_min", "early = stop_kappa is not None"),
    ('P03', 'the early exit on the non-strict comparison', "if early and snrs[-1] < stop_kappa:", "if early and snrs[-1] <= stop_kappa:"),
    ('P04', 'the early exit never fires', "            if early and snrs[-1] < stop_kappa:\n                done = True", "            if early and snrs[-1] < stop_kappa:\n                done = False"),
    ('P05', 'the early exit does not end the path', "        if done:\n            break\n", "        if done:\n            pass\n"),
    ('P06', 'the guard bar is not applied', "blocked = taken | barred if cycle_guard else taken", "blocked = taken"),
    ('P07', 'the guard\'s two-strike rule is one strike', "if abs(tau - last_exit.get(ids[i], -np.inf)) <= 1e-12 * max(1.0, abs(tau)):", "if True:"),
    ('P08', 'the last exit is not recorded', "                last_exit[ids[i]] = tau\n", "                pass\n"),
    ('P09', 'the bar is not recorded', "                    barred[ids[i]] = True\n", "                    pass\n"),
    ('P10', 'the guard acts with the switch off', "            if cycle_guard:\n                if abs(", "            if True:\n                if abs("),
    ('P11', 'an exit does not put the dual back', "            U[ids[i]] = rho\n            if cycle_guard", "            if cycle_guard"),
    ('P12', 'an exit does not release the point', "            taken[ids[i]] = False\n", "            pass\n"),
    ('P13', 'the incremental model also when the fit is global', "inc = incremental and local and len(S) < len(w)", "inc = incremental and local"),
    ('P14', 'the other units\' image is the whole model', "m_oth = (model - pre)[ya:yb, xa:xb]", "m_oth = model[ya:yb, xa:xb]"),
    ('P15', 'the carried model forgets the old image', "model = model - pre + post", "model = model + post"),
    ('P16', 'the carried model is replaced by the refit units', "model = model - pre + post", "model = post"),
    ('P17', 'a second pass uses the stale old image', "                    pre = post\n", "                    pass\n"),
    ('P18', 'the other units are the cluster itself', "oth = np.setdiff1d(np.arange(len(w)), S)", "oth = S"),
    ('P19', 'an exit does not subtract the unit\'s image', "model = model - gone_img", "model = model"),
    ('P20', 'the exit test reads a coefficient above the cut', "cand = w <= 1e-12 * max(w.max(), 1e-300)", "cand = w >= 1e-12 * max(w.max(), 1e-300)"),
    ('P21', 'the exit test reads the wrong sign of the correlation', "gone[ic] = kern_full.corr(c[ic], Rw)[0] <= vtol * vscale", "gone[ic] = kern_full.corr(c[ic], Rw)[0] >= vtol * vscale"),
    ('P22', 'the exit test uses the correlation of the first candidate only', "kern_full.corr(c[ic], Rw)[0] <= vtol * vscale", "kern_full.corr(c[ic[:1]], Rw)[0][0] <= vtol * vscale"),
    ('P23', 'the cluster is the fit radius ball of the first seed only', "S = np.flatnonzero(np.linalg.norm(c[:, None, :2] - c[seeds][None, :, :2], axis=2).min(1) < fit_radius)", "S = np.flatnonzero(np.linalg.norm(c[:, None, :2] - c[seeds][None, :1, :2], axis=2).min(1) < fit_radius)"),
    ('P24', 'the model is not refreshed from scratch at the refresh step (it must agree anyway)', "if incremental and refresh and len(hist) % refresh == 0:", "if False:"),
    ('P25', 'the time stop is not strict', "if tau_max is not None and tau + dmin > tau_max:", "if tau_max is not None and tau + dmin >= tau_max:"),
    ('P26', 'the time stop never acts', "if tau_max is not None and tau + dmin > tau_max:", "if False:"),
    ('P27', 'the block history keeps the live arrays', "hist.append(dict(tau=tau, n_ent=n_ent, c=c.copy(), w=w.copy(), chi2=chi2(), ids=list(ids), nfev=nfev))", "hist.append(dict(tau=tau, n_ent=n_ent, c=c, w=w, chi2=chi2(), ids=list(ids), nfev=nfev))"),
    ('P28', 'the significance threshold of the path\'s own stop', "if snr_min is not None and snrs[-1] < snr_min and stop_snr_at is None:", "if snr_min is not None and snrs[-1] <= snr_min and stop_snr_at is None:"),
    ('P29', 'the model is the carried one before the first fit of a global path', "                c, w, info = fit_block_fast(psf, y, swv, c, w, lo, hi, joint=joint, max_nfev=max_nfev, tol=tol, kern=kern_full)\n                model = kern_full.synth(c, w)", "                c, w, info = fit_block_fast(psf, y, swv, c, w, lo, hi, joint=joint, max_nfev=max_nfev, tol=tol, kern=kern_full)\n                model = model"),
    ('P30', 'the weights are not updated after a fit', "                swv = np.asarray(sw_fn(model) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((ny, nx))\n        Rw", "                swv = swv\n        Rw"),
]
