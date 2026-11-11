# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — round 115's revision (part VC recomputed in the tests, its first
# correction, the floor far out, the θ ~ √ε limit and its path, φ_dom, the sector's bisection, dom_bb's acceptance, the guards' plants).
M = [
 ('W00', 'control', "VC_PREC = 320 ", "VC_PREC = 320 ", 0),
 ('W01', 'part_VC: θ of the other sign', "            eta = (2 - eps) * cmath.exp(1j * eps * phi)\n            pts = []", "            eta = (2 - eps) * cmath.exp(-1j * eps * phi)\n            pts = []", 0),
 ('W02', 'part_VC: the floor at y = 0', "                y = corner_floor(x, phi)\n                E_ = cert_E", "                y = 0.0\n                E_ = cert_E", 0),
 ('W03', 'corner_phi_dom: the cubic over 3.1', "- us ** 3 / 3))", "- us ** 3 / 3.1))", 0),
 ('W04', 'corner_sector: the bisection reversed', "        a, b = (a, m) if gap_closed(R, m) > 0 else (m, b)", "        a, b = (m, b) if gap_closed(R, m) > 0 else (a, m)", 0),
 ('W05', "corner_psi1: the 2iφ term's sign", "(0.25 - 2 * phi * phi - 2j * phi)", "(0.25 - 2 * phi * phi + 2j * phi)", 0),
 ('W06', 'corner_psi_mid: D with y²', "    D, P = x - 0.5 + y * y / 2, x - 0.5 - psi * psi / 2 - psi * y", "    D, P = x - 0.5 + y * y, x - 0.5 - psi * psi / 2 - psi * y", 0),
 ('W07', 'corner_mid_path: slope ψ', "np.minimum(2 * psi / (1 + 2 * psi * psi) * xs, psi)", "np.minimum(psi * xs, psi)", 0),
 ('W08', "far: the least over y up to the floor only", "np.linspace(0.0, 3 * yf, 3001)", "np.linspace(0.0, yf, 3001)", 0),
 ('W09', 'VC_PREC 120', "VC_PREC = 320 ", "VC_PREC = 120 ", 0),
 ('W10', 'mid: over ε, not ε²', "ee.sqrt() * psi)) / ee ** 2).mid()", "ee.sqrt() * psi)) / ee).mid()", 0),
 ('W11', "dom_bb: a ball bound below 10⁻⁹ accepted", "if not ball or (u is not None and u < 0):", "if not ball or (u is not None and u < 1e-9):", 0),
 ('W12', "d_lines: 'inside' against 2", "{'inside' if c_need_ball(rs_) < R_ver else 'NOT inside'}", "{'inside' if c_need_ball(rs_) < 2 else 'NOT inside'}", 0),
 ('W13', "v_lines: 'inside' against 2", "{'inside' if c_need_ball(rv) < R_v else 'NOT inside'}", "{'inside' if c_need_ball(rv) < 2 else 'NOT inside'}", 0),
]
