"""Single-point mutants of experiments/mechreal_core.py (the numerics of the real-window mechanism test: readouts, the network's flow and its analytic Jacobian, the discrete
optimizers, the job files), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp mechreal_core --tests tests/test_mechreal.py --raw results/astro/mechreal_raw.npz --md results/astro/mechreal.md --out-env MECHREAL_OUT \\
        --spec tools/mutants/wp_mechreal_core_own.py --extra results/synthetic/network_imaging_raw.npz --work <scratch dir>

A tuple is (id, what, old, new[, k]). The kill criterion is the tests file alone (the md column is meaningless here: this module has no --from-raw). The harness's stored raw
(`network_imaging_raw.npz`) has to be in the slot for the trajectory test (--extra); without it that test is skipped and the mutants it alone kills survive.
Result (final code, run of 2026-11-10): 43 mutants, 43 killed by the tests file alone (recorded also in the docstring of `wp_mechreal_own.py`).
"""
M = [
    ('L01', 'lam_of: v_s without the factor 4', "vs = law['gamma'] * law['a_max'] ** 2 / (4.0 * law['G'])", "vs = law['gamma'] * law['a_max'] ** 2 / law['G']"),
    ('L02', 'lam_of: pi / 2', "return float(np.pi * law['a_max'] / (4.0 * law['G']) * np.sqrt(2.0 / (vs * eps)))", "return float(np.pi * law['a_max'] / (2.0 * law['G']) * np.sqrt(2.0 / (vs * eps)))"),
    ('P01', 'Paper: the rate of nu', "self.lr = (self.G * self.m, self.gam * self.m / 4.0)", "self.lr = (self.G * self.m, self.gam * self.m / 2.0)"),
    ('P02', 'Paper: start nu = eps', "return [np.zeros(d), np.full(d, np.sqrt(eps))]", "return [np.zeros(d), np.full(d, eps)]"),
    ('P03', 'Paper: w without the saturation', "return p[1] ** 2 * self.A * np.tanh(p[0] / self.A) / self.m", "return p[1] ** 2 * p[0] / self.m"),
    ('P04', 'Paper: dw/dnu halved', "return [p[1] ** 2 * ap / self.m, 2.0 * p[1] * a / self.m]", "return [p[1] ** 2 * ap / self.m, p[1] * a / self.m]"),
    ('P05', 'Paper: second derivative sign', "app = -2.0 * a * ap / self.A ** 2", "app = 2.0 * a * ap / self.A ** 2"),
    ('P06', 'Paper: active threshold', "return p[1] ** 2 > ACT_FACTOR * eps", "return p[1] ** 2 > eps"),
    ('P07', 'Paper: atol without sqrt(eps)', "return atol_p * np.sqrt(eps)                                    # the harness's rule", "return atol_p                                    # the harness's rule"),
    ('D01', 'Depth2: w = u^2 + v^2', "return p[0] ** 2 - p[1] ** 2", "return p[0] ** 2 + p[1] ** 2"),
    ('D02', 'Depth2: dw sign', "return [2.0 * p[0], -2.0 * p[1]]", "return [2.0 * p[0], 2.0 * p[1]]"),
    ('D03', 'Depth2: clock', "return float(np.log(1.0 / eps) / (4.0 * self.lr[0]))", "return float(np.log(1.0 / eps) / (2.0 * self.lr[0]))"),
    ('D04', 'Depth2: ddw', "return [[2.0 + z, z], [z, -2.0 + z]]", "return [[2.0 + z, z], [z, 2.0 + z]]"),
    ('D05', 'Depth1: rate 1', "self.lr = (law['G'] * law['m'],)\n        self.lr_c = law['G'] * law['m']\n        self._law = law", "self.lr = (1.0,)\n        self.lr_c = law['G'] * law['m']\n        self._law = law"),
    ('D06', 'Depth1: atol with sqrt(eps)', "return atol_p                                                   # the amplitudes are O(1)", "return atol_p * np.sqrt(eps)                                                   # the amplitudes are O(1)"),
    ('E01', 'fast axis: the edges shifted', "e = np.r_[p - 0.5, p[-1] + 0.5]", "e = np.r_[p, p[-1] + 1.0]"),
    ('E02', 'fast axis: derivative sign', "return 0.5 * (er[1:] - er[:-1]), -(ex[1:] - ex[:-1]) / (_SQ2PI * s)", "return 0.5 * (er[1:] - er[:-1]), (ex[1:] - ex[:-1]) / (_SQ2PI * s)"),
    ('E03', 'fast comps: the width of a part', "s = kap * P.sigma\n        fx, fy", "s = P.sigma\n        fx, fy"),
    ('E04', 'fast comps: the high-pass partner with the wrong sign', "out.append((-a, P.Sy @ fy[0], P.Sx @ fx[0], P.Sy @ fy[1], P.Sx @ fx[1]))", "out.append((a, P.Sy @ fy[0], P.Sx @ fx[0], P.Sy @ fy[1], P.Sx @ fx[1]))"),
    ('E05', 'fast comps: x and y exchanged', "fx, fy = _axis_edges(P.px, c[:, 0], s), _axis_edges(P.py, c[:, 1], s)", "fx, fy = _axis_edges(P.px, c[:, 1], s), _axis_edges(P.py, c[:, 0], s)"),
    ('N01', 'pieces: residual without the weights', "Rw = self.W * (self.y - model)\n        v = vx = vy = 0.0", "Rw = (self.y - model)\n        v = vx = vy = 0.0"),
    ('N02', 'pieces: vx and vy exchanged', "v, vx, vy = v + co * (A_ * gx.T).sum(1), vx + co * (A_ * dgx.T).sum(1), vy + co * (B_ * gx.T).sum(1)\n        return comps, p, w, Rw, v, vx, vy",
     "v, vy, vx = v + co * (A_ * gx.T).sum(1), vx + co * (A_ * dgx.T).sum(1), vy + co * (B_ * gx.T).sum(1)\n        return comps, p, w, Rw, v, vx, vy"),
    ('N03', 'neg_grad: the centres without w', "return np.concatenate([dwl * v for dwl in self.ro.dw(p)] + [w * vx, w * vy])", "return np.concatenate([dwl * v for dwl in self.ro.dw(p)] + [vx, vy])"),
    ('N04', 'rates: the centres at rate 1', "return np.concatenate([np.full(d, r) for r in self.ro.lr] + [np.full(2 * d, self.ro.lr_c)])", "return np.concatenate([np.full(d, r) for r in self.ro.lr] + [np.full(2 * d, 1.0)])"),
    ('N05', 'init: the centres exchanged', "return np.concatenate(self.ro.init(self.d, self.eps) + [self.c0[:, 0], self.c0[:, 1]])", "return np.concatenate(self.ro.init(self.d, self.eps) + [self.c0[:, 1], self.c0[:, 0]])"),
    ('J01', 'jac: the Gram term of the parameter columns dropped', "blk = -lr[l] * dw[l][:, None] * G00 * dwk[None, :]", "blk = 0.0 * dw[l][:, None] * G00 * dwk[None, :]"),
    ('J02', 'jac: own centre term of the x column', "bx[act, ar] += lc * wa * own_x[act]", "bx[act, ar] += lc * wa * own_y[act]"),
    ('J03', 'jac: own block (x, p) dropped', "vals += [lc * dw[k] * vx, lc * dw[k] * vy]", "vals += [0.0 * dw[k] * vx, lc * dw[k] * vy]"),
    ('J04', 'jac: own block (p, p) transposed rates', "vals.append(lr[l] * ddw[l][k] * v)", "vals.append(lr[k] * ddw[l][k] * v)"),
    ('F01', 'run_flow: the snapshot time without the clock', "z, t, status, c_ = advance(net, z, t, clock * tau, rtol=rtol)", "z, t, status, c_ = advance(net, z, t, tau, rtol=rtol)"),
    ('F02', 'run_flow: a failed run counts as done', "ok = status == 'ok' and bool(np.all(np.isfinite(z)))", "ok = True"),
    ('F03', 'advance: no budget', "        except _Budget:\n            return z, t, 'timeout', cnt", "        except _Budget:\n            return z, t, 'ok', cnt"),
    ('O01', 'Opt adam: no bias correction of the second moment', "mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)", "mh, vh = self.m / (1 - self.b1 ** self.t), self.v"),
    ('O02', 'Opt adam_g: per-coordinate second moment', "gg = gneg ** 2 if self.kind == 'adam' else np.mean(gneg ** 2) * np.ones_like(gneg)", "gg = gneg ** 2"),
    ('O03', 'Opt adam_g: the rates not applied', "if self.pre is not None and self.kind == 'adam_g':\n            gneg = self.pre * gneg", "if self.pre is not None and self.kind == 'adam':\n            gneg = self.pre * gneg"),
    ('O04', 'run_discrete: one step short', "while status == 'ok' and k < target:", "while status == 'ok' and k < target - 1:"),
    ('O05', 'run_discrete: adam_g for both', "opt = Opt(kind, lr, pre=net.rates() if kind == 'adam_g' else None)", "opt = Opt('adam_g', lr, pre=net.rates() if kind == 'adam_g' else None)"),
    ('K01', 'pack_result: keeps |w| >= 0', "f = np.flatnonzero(np.isfinite(Wk) & (np.abs(Wk) > thr))", "f = np.flatnonzero(np.isfinite(Wk))"),
    ('K02', 'pack_result: arm A in single precision', "dt = np.float64 if arm == 'A' else np.float32", "dt = np.float32"),
    ('K03', 'unpack_result: snapshots shifted', "snaps = [(np.asarray(raw[pre + 'idx'][e - k:e], int)", "snaps = [(np.asarray(raw[pre + 'idx'][e - k:e][::-1], int)"),
    ('K04', 'work: reruns finished jobs', "todo = [s for s in specs if 'R|%s|run' % s['key'] not in done]", "todo = list(specs)"),
    ('K05', 'make_problem: weights not squared', "return dict(P=P, lat=lat, c0=np.c_[X.ravel(), Y.ravel()], y=y, sw=swv, W=swv ** 2)", "return dict(P=P, lat=lat, c0=np.c_[X.ravel(), Y.ravel()], y=y, sw=swv, W=swv)"),
    ('K06', 'make_problem: lattice order', "return dict(P=P, lat=lat, c0=np.c_[X.ravel(), Y.ravel()], y=y, sw=swv, W=swv ** 2)", "return dict(P=P, lat=lat, c0=np.c_[Y.ravel(), X.ravel()], y=y, sw=swv, W=swv ** 2)"),
]
