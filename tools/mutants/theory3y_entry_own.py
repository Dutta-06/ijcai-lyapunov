# (id, what, old, new, k): single-point mutants of experiments/theory3y_entry.py (the unit on the critical curve, the layer data and the stress layers, the box search, the ball-arithmetic bounds of Theorem E3, the layer
# integrator, the checks of Theorems E0-E4, the aligned / orbit distances, the tail constant).  Run:
#   venv/bin/python tools/mutants.py --exp theory3y_entry --tests tests/test_theory3y.py --raw results/synthetic/theory3y_raw.npz --md results/synthetic/theory3y.md --out-env T3Y_OUT \
#       --spec tools/mutants/theory3y_entry_own.py --work <scratch>/mut --slot a
# (the kill criterion is the tests file alone; `--from-raw` is not a mode of the library file: fr_rc is not informative for this spec).
# Result (66 entries, the control and 65 mutants; first pass 53 killed, then tests added for Y11 (signs of the layer's units: target with a non-identity feature order and mixed signs), Y74 and Y76 (tampered
#   trajectories) and re-run with --only: 56 of 65 killed; control passes).  Nine survivors, all equivalent on the data: Y23 (the separation is not used when K = 1); Y31 (upper instead of lower end of the ball
#   M_lower: the ball's radius is about 1e-17 relative); Y37 (u = max(d0, d-tilde) against d0: ln+ of a ratio below 1 is 0 either way, guarded by `u > dt`); Y38 (the centre mobility mu_c a_min^2 = 4e2 never
#   is the smaller of the two mobilities, M about 0.07); Y50 (a_max = 1 in the repository's law: G/a_max = G); Y53 (the diagonal of Gamma' is 0 already); Y77 (the tail check with c_M w_#|s| instead of
#   c_M w_#|s|/3: the actual tail is 28 times better than the lemma's constant, so the stronger inequality also holds); Y92 (levels above the end value: both flows are at the end value there, distance 0);
#   Y93 (the entrant's coefficient is monotone on the sampled flows, the running maximum changes nothing).
M = [
 ('Y00', 'control (unmutated)', '', '', 0),
 # ---- the unit
 ('Y01', 'W: no eps term', "np.tanh(x) * (eps + 2 * VS * np.sinh(x) ** 2)", "np.tanh(x) * (2 * VS * np.sinh(x) ** 2)", 0),
 ('Y02', 'W: x instead of tanh x', "return (LAW['a_max'] / LAW['m']) * np.tanh(x) * (eps", "return (LAW['a_max'] / LAW['m']) * x * (eps", 0),
 ('Y03', 'x_of_w: bisection bracket too small', "lo, hi = 0.0, 30.0", "lo, hi = 0.0, 3.0", 0),
 ('Y04', 'x_of_w: sign lost', "s = 1.0 if w >= 0 else -1.0", "s = 1.0", 0),
 ('Y05', 'M: 2 v_s instead of 4 v_s', "(q + 4 * VS * np.sinh(x) ** 2)", "(q + 2 * VS * np.sinh(x) ** 2)", 0),
 ('Y06', 'M: q without the -eps', "q = eps + (2 * VS - eps) * np.tanh(x) ** 2\n    return (LAW['G'] / LAW['m'])", "q = eps + (2 * VS) * np.tanh(x) ** 2\n    return (LAW['G'] / LAW['m'])", 0),
 ('Y07', 'M: G/m inverted', "return (LAW['G'] / LAW['m']) * q * (q +", "return (LAW['m'] / LAW['G']) * q * (q +", 0),
 # ---- layers
 ('Y10', 'layers: block features off by one', "feats = [e['feature'] for e in ev[:k + 1]]", "feats = [e['feature'] for e in ev[:k]]", 0),
 ('Y11', 'layers: signs of the target not of the units', "sgn=np.sign(As[feats])", "sgn=np.sign(As[:k + 1])", 0),
 ('Y12', 'layers: entrant offset of the previous event', "d0=float(ev[k]['d0'])", "d0=float(ev[k - 1]['d0'])", 0),
 ('Y13', 'layers: block coefficients of the current event', "wb=np.asarray(ev[k - 1]['w']) if k > 0 else np.zeros(0)", "wb=np.asarray(ev[k]['w'][:k]) if k > 0 else np.zeros(0)", 0),
 ('Y14', 'layers: entrant kernel centre from the feature', "ce=float(c0[ev[k]['kernel']])", "ce=float(cs[ev[k]['feature']])", 0),
 ('Y15', 'layers: end point of the previous layer', "what=np.asarray(ev[k]['w']), chat=np.asarray(ev[k]['c'])", "what=np.asarray(ev[k - 1]['w']), chat=np.asarray(ev[k - 1]['c'])", 0),
 ('Y16', 'stress: entrant away from the block', "ce=float(cs[1] - STRESS_D0)", "ce=float(cs[1] + STRESS_D0)", 0),
 ('Y17', 'stress: no opposite signs', "As = a * np.array([1.0, -1.0 if opp else 1.0])", "As = a * np.array([1.0, 1.0])", 0),
 ('Y18', 'stress: block fits the second feature', "wb, cb = T.crit_block(cs, As, [0])", "wb, cb = T.crit_block(cs, As, [1])", 0),
 ('Y19', 'stress: offset 0.30', "STRESS_D0 = 0.40", "STRESS_D0 = 0.30", 0),
 # ---- the box search
 ('Y20', 'layer_box: radius not above d0', "if r <= L['d0'] + 1e-9:\n            continue", "if r <= 0.0:\n            continue", 0),
 ('Y21', 'layer_box: no (N-hat)', "if 'nu' in b and b['nu'] > 0 and b['vD_lo'] > 0:", "if 'nu' in b and b['vD_lo'] > 0:", 0),
 ('Y22', 'layer_box: dormant budget 0.05', "T.box_t1(K, L['R'], Delta, r, 0.0, DHAT_FRAC, VSIG)", "T.box_t1(K, L['R'], Delta, r, 0.05, DHAT_FRAC, VSIG)", 0),
 ('Y23', 'layer_box: finite Delta for K = 1 replaced', "Delta = L['Delta'] if np.isfinite(L['Delta']) else 0.0", "Delta = L['Delta'] if np.isfinite(L['Delta']) else 6.0", 0),
 # ---- the bounds
 ('Y30', 'M_lower: wrong direction of the certain test', "if bool(_arb_W(xl, eps) < arb(w)):", "if bool(_arb_W(xl, eps) > arb(w)):", 0),
 ('Y31', 'M_lower: upper end', "return float(_arb_M(xl, eps).lower())", "return float(_arb_M(xl, eps).upper())", 0),
 ('Y32', 'bounds: level 0.4 a_min v_lo', "w_hash = 0.5 * amin * vlo", "w_hash = 0.4 * amin * vlo", 0),
 ('Y33', 'bounds: kd without the 1/2', "(arb(-r * r / 4)).exp() / 2", "(arb(-r * r / 4)).exp()", 0),
 ('Y34', 'bounds: T2 against r instead of d-hat', "(u - dt) / (dh - dt)", "(u - dt) / (r - dt)", 0),
 ('Y35', 'bounds: T3 with the upper-face mobility', "np.log(max((amin * Vt - w_hash) / (amin * VSIG), 1.0)) / mV if", "np.log(max((amin * Vt - w_hash) / (amin * VSIG), 1.0)) / mVp if", 0),
 ('Y36', 'bounds: V_t without E_v', "Vt = G2(dh) - Ev", "Vt = G2(dh)", 0),
 ('Y37', 'bounds: T2 uses d0 only', "u = max(L['d0'], dt)", "u = L['d0']", 0),
 ('Y38', 'bounds: centre mobility not a_min squared', "m_lo = min(M_lower(amin * box['vD_lo'], eps), MC * amin ** 2)", "m_lo = min(M_lower(amin * box['vD_lo'], eps), MC)", 0),
 ('Y39', 'bounds: rate without nu', "rho = m_lo * nu", "rho = m_lo", 0),
 ('Y40', 'bounds: prefactor sqrt(Hup/nu) -> Hup/nu', "pref=float(np.sqrt(Hup / nu))", "pref=float(Hup / nu)", 0),
 ('Y41', 'bounds: lower mobility from the end of the range', "M_lower(amin * box['vD_lo'], eps)", "M_lower(amin * box['vD_lo'] * 2, eps)", 0),
 # ---- the integrator
 ('Y50', 'rhs: no a_max in the x-flow', "return np.r_[(G / amax) * (eps + (2 * VS - eps) * np.tanh(x) ** 2) * gg, MC * w * h]", "return np.r_[G * (eps + (2 * VS - eps) * np.tanh(x) ** 2) * gg, MC * w * h]", 0),
 ('Y51', 'rhs: centre mobility G m', "MC * w * h]", "w * h]", 0),
 ('Y52', 'rhs: drive without the block', "gg = Gam(c[:, None] - cs[None, :]) @ As - Gam(c[:, None] - c[None, :]) @ w", "gg = Gam(c[:, None] - cs[None, :]) @ As", 0),
 ('Y53', 'rhs: self term in the centre force', "np.fill_diagonal(M, 0.0)", "pass", 0),
 ('Y54', 'initial: entrant not dormant', "        xe = 0.0\n", "        xe = 1e-3\n", 0),
 ('Y55', 'initial: v0 sign lost', "xe = x_of_w(L['sgn'][-1] * v0 * L['amin'], 0.0)", "xe = x_of_w(v0 * L['amin'], 0.0)", 0),
 ('Y56', 'integrate: event direction', "lev.direction = 1", "lev.direction = -1", 0),
 ('Y57', 'integrate: level of the event', "return abs(float(W_of_x(z[n - 1], eps))) - w_hash", "return abs(float(W_of_x(z[n - 1], eps))) - 0.5 * w_hash", 0),
 ('Y58', 'integrate: coarse pre-crossing sampling', "np.geomspace(1e-3, th, 4000)", "np.geomspace(1e-3, th, 40)", 0),
 ('Y59', 'derived: no normalisation by a_min', "v = (L['sgn'][:, None] * w) / L['amin']", "v = (L['sgn'][:, None] * w)", 0),
 ('Y60', 'derived: offsets from the first features', "d = tr['c'] - L['cs'][L['feats']][:, None]", "d = tr['c'] - L['cs'][:n][:, None]", 0),
 # ---- the checks
 ('Y70', 'check: box radius doubled', "ok['box_d'] = float(np.abs(d).max()) < r", "ok['box_d'] = float(np.abs(d).max()) < 2 * r", 0),
 ('Y71', 'check: no tolerance on v_bar', "float(v.max()) <= vbar * (1 + 1e-9)", "float(v.max()) <= vbar", 0),
 ('Y72', 'check: block lower face from the entrant', "bool(n == 1 or v[:-1].min() >= vlo)", "bool(n == 1 or v[:-1].min() >= 0.0)", 0),
 ('Y73', 'check: monotonicity tolerance', "np.all(np.diff(ve[:i_lo + 1]) >= -1e-13)", "np.all(np.diff(ve[:i_lo + 1]) >= -1.0)", 0),
 ('Y74', 'check: entrant above v_lo afterwards', "ve[i_lo:].min() >= vlo - 1e-12", "ve[i_lo:].min() >= 0.0", 0),
 ('Y75', 'check: relaxation set without the margin', "(v >= np.exp(-dh * dh / 4) - Ev - vsig)", "(v >= np.exp(-dh * dh / 4) - Ev)", 0),
 ('Y76', 'check: block in D-hat not required', "ok['block_in_hat'] = bool(in_hat[:-1].all()) if n > 1 else True", "ok['block_in_hat'] = True", 0),
 ('Y77', 'check: tail constant /1', "c_tail(wh) * wh / 3.0 * np.abs(s[pre_])", "c_tail(wh) * wh * np.abs(s[pre_])", 0),
 ('Y78', 'check: decay slack', "viol = err[i_hat:] - env - 1e-9", "viol = err[i_hat:] - env - 1.0", 0),
 ('Y79', 'check: decay from the crossing', "env = bd['pref'] * err[i_hat] * np.exp(-bd['rho'] * (s[i_hat:] - s[i_hat]))", "env = bd['pref'] * err[i_hat] * np.exp(-bd['rho'] * (s[i_hat:]))", 0),
 ('Y80', 'check: D* without E_v', "in_star = (np.abs(d) <= dt) & (v >= np.exp(-dt * dt / 4) - Ev) & (v <= R + Ev)", "in_star = (np.abs(d) <= dt) & (v >= np.exp(-dt * dt / 4)) & (v <= R + Ev)", 0),
 ('Y81', 'check: bound of the D* time without the diameter', "diam = float(np.sqrt(n * (wid ** 2 + (2 * dh) ** 2)))", "diam = float(wid)", 0),
 ('Y82', 'check: margin of D* from the lower face only', "marg = float(min(np.min(dt - np.abs(dhat_c)), np.min(vh - (np.exp(-dt * dt / 4) - Ev)), np.min(R + Ev - vh)))", "marg = float(np.min(vh - (np.exp(-dt * dt / 4) - Ev)))", 0),
 # ---- the distances
 ('Y90', 'aligned distance: mean instead of sup', "return float(np.sqrt(((A - B) ** 2).sum(axis=0)).max())\n\n\ndef orbit_distance", "return float(np.sqrt(((A - B) ** 2).sum(axis=0)).mean())\n\n\ndef orbit_distance", 0),
 ('Y91', 'aligned distance: block coordinates dropped', "n_u = trA['n']\n    VA, DA = derived(LA, trA, epsA, None)\n    VB, DB = derived(LA, trB, epsB, None)\n    A =", "n_u = 1\n    VA, DA = derived(LA, trA, epsA, None)\n    VB, DB = derived(LA, trB, epsB, None)\n    A =", 0),
 ('Y92', 'orbit distance: levels up to the end value', "top = frac * float(np.abs(LA['what'][-1]) / LA['amin'])", "top = 2.0 * float(np.abs(LA['what'][-1]) / LA['amin'])", 0),
 ('Y93', 'orbit distance: no monotone envelope', "ve = np.maximum.accumulate(V[-1])", "ve = V[-1]", 0),
 # ---- the tail constant, the random targets
 ('Y95', 'c_tail: cosh squared', "np.cosh(xh) ** 4", "np.cosh(xh) ** 2", 0),
 ('Y96', 'c_tail: exponent 1/2', "** (2.0 / 3.0)\n", "** 0.5\n", 0),
 ('Y97', 'random targets: seed', "R_SEED = 20260111", "R_SEED = 20260112", 0),
 ('Y98', 'random targets: separation range', "rng.uniform(5.5, 9.5)", "rng.uniform(5.5, 9.0)", 0),
]
