# (id, what, old, new, k): single-point mutants of experiments/theory3n_system.py (the corrected certificate: anchored source repairs (f) stage-resolved displacement, (g) pre-layer term of windows not
# started, (6) stage maximum of the old units' coefficient, kappa by the lower end, certain `inv_of`; the float / planar / ball builders; the chain cover).  The `old` strings are the text of the
# replacement strings of the repairs, i.e. of the code that is compiled into the corrected functions.
# Run: tools/mutants.py --exp theory3n_system --tests tests/theory3n_test.py --raw results/synthetic/theory3n_raw.npz --md results/synthetic/theory3n.md --out-env T3N_OUT
#      --extra results/synthetic/theory3n.md --extra tests/theory3c_test.py --spec tools/mutants/theory3n_system_own.py --work <scratch dir>
# Result (32 mutants, criterion: the tests file alone, tests/theory3n_test.py): 28 killed (C02, a failed ball width that is never halved, by NON-TERMINATION: the run was stopped by hand after 400 s and counts as
#   killed by timeout).  Survivors, equivalent on the data: V04 (stage maximum of the balls compared by the lower instead of the upper ends: the balls are points to 1e-60);  K01 (kappa = the first of the two
#   stability constants: k1 < k2 by a factor 200 on all 8 targets, so the selection never matters; the selection rule itself is unit-tested through `min_lower`);  C03 (the candidate bounds D_j of the cover taken
#   one decade lower: a candidate only proposes, the ball verifier decides, and the cover still closes);  C04 (cover balls built from the rounded-to-nearest powers: the ball of two doubles is padded outward by
#   about 1e-9 of its width, as for O03 of theory3c).
M = [
 ('H00', 'control (unmutated)', '', '', 0),
 ('F01', 'float (f): the earlier stages\' mass dropped (H^(st) on the cumulative mass, the legacy form)', "dcb = muc * (Wprev + Hst[st] * (Ipts - Icum) + hinc)", "dcb = muc * (Hst[st] * Ipts + hinc)", 0),
 ('F02', 'float (f): the stage-resolved sum is not accumulated', "Wprev = Wprev + Hst[st] * (float(Ipts[-1]) - Icum)", "Wprev = Wprev", 0),
 ('F03', 'float (f): the accumulation without the stage weight', "Wprev = Wprev + Hst[st] * (float(Ipts[-1]) - Icum)", "Wprev = Wprev + (float(Ipts[-1]) - Icum)", 0),
 ('G01', 'float (g): the pre-layer term of unstarted windows dropped', " + float(sum(psi_sup(Xj)[i] * pre_mass for i in range(K) if not (wstart[i] <= q0 + 1e-9)))", "", 0),
 ('G02', 'float (g): the term over the STARTED windows', "pre_mass for i in range(K) if not (wstart[i] <= q0 + 1e-9)))", "pre_mass for i in range(K) if (wstart[i] <= q0 + 1e-9)))", 0),
 ('G03', 'float (g): psi_sup replaced by Gamma', "float(sum(psi_sup(Xj)[i] * pre_mass", "float(sum(Gam(Xj)[i] * pre_mass", 0),
 ('V01', 'float (6): the final block\'s coefficient again', "    vhat_fin = _stage_max_vhat(ev, K)\\n", "    vhat_fin = np.abs(ev[-1]['w'])\\n", 0),
 ('V02', 'stage maximum (float) replaced by the minimum', "max(abs(float(e['w'][i])) for e in ev if len(e['w']) > i)", "min(abs(float(e['w'][i])) for e in ev if len(e['w']) > i)", 0),
 ('V03', 'stage maximum (balls): the last stage only', "m = v if m is None or R3.up(v) >= R3.up(m) else m", "m = v", 0),
 ('V04', 'stage maximum (balls): by the lower ends', "m = v if m is None or R3.up(v) >= R3.up(m) else m", "m = v if m is None or R3.lo(v) >= R3.lo(m) else m", 0),
 ('R01', 'rig (f): the legacy form', "dcb = muc * (Wprev + Hst[st] * (Ig - Icum) + hinc)", "dcb = muc * (Hst[st] * Ig + hinc)", 0),
 ('R02', 'rig (f): the sum is not accumulated', "Wprev = up(Wprev + Hst[st] * (up(last) - Icum))", "Wprev = Wprev", 0),
 ('R03', 'rig (g): the term dropped', "+ sum((hpre[i] for i in range(nw, K)), arb(0))", "", 0),
 ('R04', 'rig (g): pre-layer terms are zero', "hpre.append(psi_sup_r(Xl) * pre)", "hpre.append(arb(0))", 0),
 ('R05', 'rig (g): one window too few', "for i in range(nw, K)), arb(0)) + arb(43)", "for i in range(nw + 1, K)), arb(0)) + arb(43)", 0),
 ('R06', 'rig (g): psi_sup replaced by Gamma', "hpre.append(psi_sup_r(Xl) * pre)", "hpre.append(Gam_r(Xl) * pre)", 0),
 ('K01', 'rig: kappa by the first', "kappa = _min_lower(k1, k2)", "kappa = k1", 0),
 ('I01', 'inv_of: only the positive side is certain', "if R3.lo(d_) > tol or R3.up(d_) < -tol:", "if R3.lo(d_) > tol:", 0),
 ('I02', 'inv_of: the tolerance multiplied instead of divided', "tol = R3.up(abs(x) + abs(y)) / 1000000", "tol = R3.up(abs(x) + abs(y)) * 1000000", 0),
 ('I03', 'inv_of: no tolerance', "tol = R3.up(abs(x) + abs(y)) / 1000000", "tol = arb(0)", 0),
 ('B01', 'best_cert_n uses the legacy certificate', "    return _best(cert_n, prep, eps, pars)", "    return _best(T.cert, prep, eps, pars)", 0),
 ('B02', 'the planar certificate without the repairs', "cert2d_n = build_cert(planar=True)", "cert2d_n = build_cert(True, True, True, planar=True)", 0),
 ('B03', 'the ball verifier without the repairs', "verify_rig_n = build_rig()", "verify_rig_n = build_rig(True, True, True, True, True)", 0),
 ('B04', 'build_rig applies the repairs that are switched off', "        if not flag:\n            src = _apply(src, ed)", "        if flag:\n            src = _apply(src, ed)", 0),
 ('B05', 'anchor test accepts a missing anchor', "    if src.count(a) != 1:\n            raise RuntimeError", "    if src.count(a) < 0:\n            raise RuntimeError", 0),
 ('B06', 'the cert builder applies (g) when switched off', "    if not legacy_g:\n        src = _apply(src, G_FLOAT)", "    if legacy_g:\n        src = _apply(src, G_FLOAT)", 0),
 ('P01', 'planar: the norm of the centre differences replaced', "_nrm(c0[:, None, :] - cs[featarr][None, :, :])", "np.abs(c0[:, None, 0] - cs[featarr][None, :, 0])", 0),
 ('C01', 'cover: the lowest ball is not clipped at the bottom', "lo_ = max(cur - w, e_bot)", "lo_ = cur - w", 0),
 ('C02', 'cover: a failed width is not halved', "            w /= 2\n", "            w /= 1\n", 0),
 ('C03', 'cover: candidate bounds from the lower end', "o, pi = best_cert_n(prep, 10.0 ** cur)", "o, pi = best_cert_n(prep, 10.0 ** (cur - 1.0))", 0),
 ('C04', 'cover: the ball is not outward rounded', "ball = R3.eps_ball(lo_, cur)", "ball = arb(10.0 ** lo_).union(arb(10.0 ** cur))", 0),
 ('C05', 'cover: bounds without the factor', "Dv = [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]", "Dv = [1.0 if j in ent else max(o['Etot'][j], 1.0) for j in range(len(prep['c0']))]", 0),
]
