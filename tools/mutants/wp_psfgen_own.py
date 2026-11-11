"""Single-point mutants of src/core/psf_general.py, for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/psf_general --tests tests/test_psf_general.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp_psfgen_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(the tool's md / `--from-raw` columns mean nothing for these mutants: only `tests` does — True: tests/test_psf_general.py kills the mutant). A tuple is (id, what, old, new[, k]): `old`
replaced by `new` at its k-th occurrence in src/core/psf_general.py.

Result: 34 of the 40 mutants are killed by tests/test_psf_general.py. The six that survive are equivalent or insensitive on every test stamp: S04 (the weights' bounds already keep them
non-negative) and S08 (the goal met at exactly equality, a measure-zero boundary) are equivalent; S09, S10, S13 and S17 change only where or from which weight/width the greedy step STARTS
a new part, and the least-squares refit of all parameters reaches the same fit from the other starts on the stamps tested (the initialisation is a heuristic, the fit's objective is what
the tests pin).
"""
M = [
    ('G01', 'the base class keeps its isotropic parts (a base method would return another kernel)', "self.parts = None", "self.parts = ((1.0, 1.0),)"),
    ('G02', 'window: the origin of the sub-frame is lost', "self.y0 + ya, self.x0 + xa, _parent=root)", "self.y0, self.x0, _parent=root)"),
    ('G03', 'window of a window: the crop is not composed', "else (ya + self._crop[0], yb + self._crop[0], xa + self._crop[2], xb + self._crop[2])", "else (ya, yb, xa, xb)"),
    ('G04', 'reach: the offsets are ignored', "r = max(max(kx, ky) + max(abs(dx), abs(dy)) for a, kx, ky, dx, dy in self.parts5)", "r = max(max(kx, ky) for a, kx, ky, dx, dy in self.parts5)"),
    ('G05', 'reach: the smaller offset', "max(kx, ky) + max(abs(dx), abs(dy))", "max(kx, ky) + min(abs(dx), abs(dy))"),
    ('G06', 'reach: the smaller width', "max(kx, ky) + max(abs(dx), abs(dy))", "min(kx, ky) + max(abs(dx), abs(dy))"),
    ('G07', 'reach: the high-pass added linearly', "np.hypot(r, self.hp) if self.hp else r", "(r + self.hp) if self.hp else r"),
    ('G08', 'comps: the x offset has the wrong sign', "c[:, 0] + dx if dx else c[:, 0]", "c[:, 0] - dx if dx else c[:, 0]"),
    ('G09', 'comps: the y offset is the x offset', "c[:, 1] + dy if dy else c[:, 1]", "c[:, 1] + dx if dy else c[:, 1]"),
    ('G10', 'comps: the x width is the y width', "np.full(len(c), kx), False)", "np.full(len(c), ky), False)"),
    ('G11', 'comps: the y width is the x width', "np.full(len(c), ky), False)", "np.full(len(c), kx), False)"),
    ('G12', 'comps: the high-pass partner has the sign of the part', "terms.append((-a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))", "terms.append((a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))"),
    ('G13', 'comps: the partner smooths the axes the wrong way', "terms.append((-a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))", "terms.append((-a, [root.Sx @ f for f in fy], [root.Sy @ f for f in fx]))"),
    ('G14', 'comps: the derivative columns are exchanged', "fx, fy = [ax[0], ax[1]], [ay[0], ay[1]]", "fx, fy = [ax[0], ay[1]], [ay[0], ax[1]]"),
    ('G15', 'comps: a per-kernel width derivative is accepted', "        assert not deriv_s, 'no per-kernel width'\n", ""),
    ('G16', 'comps: a kernel with a width of its own is accepted', "        assert c.shape[1] == 2\n", ""),
    ('G17', 'reach: a width is accepted', "assert s is None or s == self.sigma, 'no per-kernel width'", "pass"),
    ('S01', 'isotropic parts: the widths are not times sigma', "float(k) * float(sigma), float(k) * float(sigma), 0.0, 0.0", "float(k), float(k), 0.0, 0.0"),
    ('S02', 'stamp model: the kernel is centred at n / 2', "np.array([[(n - 1) / 2.0, (n - 1) / 2.0]])", "np.array([[n / 2.0, n / 2.0]])"),
    ('S03', 'normalisation of the weights dropped', "P[:, 0] = w / w.sum()", "P[:, 0] = w"),
    ('S04', 'the weights are not made positive (bounds keep them so: equivalent)', "w = np.abs(P[:, 0])", "w = P[:, 0]"),
    ('S05', 'the lower bound of the widths is 0', "np.tile([0.0, W_LO, W_LO, -D_MAX, -D_MAX], K)", "np.tile([0.0, 0.0, 0.0, -D_MAX, -D_MAX], K)"),
    ('S06', 'the offsets are bounded by the widths bound', "np.tile([1.0, W_HI, W_HI, D_MAX, D_MAX], K)", "np.tile([1.0, W_HI, W_HI, W_HI, W_HI], K)"),
    ('S07', 'one part more than k_max', "while rms > rms_goal and len(P) < k_max:", "while rms > rms_goal and len(P) <= k_max:"),
    ('S08', 'the goal is met only strictly beyond', "while rms > rms_goal and len(P) < k_max:", "while rms >= rms_goal and len(P) < k_max:"),
    ('S09', 'the part goes where the signed residual is largest', "np.argmax(np.abs(R))", "np.argmax(R)"),
    ('S10', 'the first part lands at the stamp corner', "c = (n - 1) / 2.0\n", "c = 0.0\n"),
    ('S11', 'a part is accepted that does not lower the rms by 1e-6', "if best[0] >= rms - 1e-6:", "if best[0] >= rms + 1.0:"),
    ('S12', 'a part is kept that makes the fit worse', "if best[0] >= rms - 1e-6:\n            break\n        rms, P = best", "rms, P = best"),
    ('S13', 'only the weakest start weight is tried', "for w0 in (0.01, 0.03, 0.08):", "for w0 in (0.01,):"),
    ('S14', 'the adequacy rule reversed', "if rms0 <= rms_goal:", "if rms0 > rms_goal:"),
    ('S15', 'the adequacy goal ignored in the general fit', "isotropic_as_general(sigma, parts), rms_goal, k_max)", "isotropic_as_general(sigma, parts), RMS_GOAL, k_max)"),
    ('S16', 'the part bound ignored in the general fit', "isotropic_as_general(sigma, parts), rms_goal, k_max)", "isotropic_as_general(sigma, parts), rms_goal, K_MAX)"),
    ('S17', 'the general fit starts from nothing but the isotropic start of another width', "fit_general_parts(stamp, isotropic_as_general(sigma, parts), rms_goal, k_max)", "fit_general_parts(stamp, isotropic_as_general(sigma * 1.5, parts), rms_goal, k_max)"),
    ('S18', 'the kind of the isotropic fit is called general', "return dict(kind='isotropic', sigma=float(sigma)", "return dict(kind='general', sigma=float(sigma)"),
    ('S19', 'make_psf: the kinds exchanged', "if fit['kind'] == 'isotropic':", "if fit['kind'] != 'isotropic':"),
    ('S20', 'make_psf: the high-pass dropped for the isotropic kernel', "parts=fit['parts'], hp=hp or None)", "parts=fit['parts'], hp=None)"),
    ('S21', 'make_psf: the high-pass dropped for the general kernel', "GeneralPixelPSF(ny, nx, fit['parts5'], hp=hp or None)", "GeneralPixelPSF(ny, nx, fit['parts5'], hp=None)"),
    ('S22', 'the constants: goal 5e-3', "RMS_GOAL = 2e-3 ", "RMS_GOAL = 5e-3 "),
    ('S23', 'the constants: at most 6 parts', "K_MAX = 10 ", "K_MAX = 6 "),
]
