# (id, what, old, new, k): single-point mutants of exp_transit_limit.py — round 140's revision (referee D3ei), never mutation-tested
# when applied: `entry_shift` (σ_k = −s_k·dev_{j_k}/(ρ_{j_k}v_k)); Part U (`u_masses`, `path_duals_rho`, `u_job`: the walls ρ_j = b_j^{−1/2},
# the margin over ρ_j, the relative speed, the entrant's ρ in the pulls, σ_k with and without 1/ρ, the recursion of dev); Part W (`w_job`,
# `part_W`: the window start τ_k + (C_W/λ_k)ε^{1/2}log(1/ε), the skipped windows, the block) and the W and U sections of `report`.
# Run with -k "entry_shift or unequal or path_duals_rho or part_W or report_W" (the round's tests; the older ones take ~1.5 min more per
# mutant).
# Result: 26 mutants and a control. Against the tests as committed (-k "entry_shift or unequal"): 11 killed by the tests (L01–L03,
# L06–L12, L15), 5 by the md (L21, L22, L24–L26); the 10 others (L04, L05, L13, L14, L16–L20, L23; L18 at its w_job occurrence, k = 1)
# killed by the tests added then: Part U's stored entries and row reproduced in test_unequal_masses_the_shift_needs_one_over_rho;
# test_path_duals_rho_pins_the_entrants_at_their_walls (L04, L05: `u_job` reads only the signs of the pinned duals);
# test_part_W_reproduces_the_stored_windows (with a planted λ_2 for L18: no stored window starts between τ_{k+1} − ε^{1/4} and
# τ_{k+1} − ε^{1/2}); test_report_W_takes_the_transient_over_the_intervals_present (L23: λ_k decreases on every stored target, so the
# last interval, always present, decides). No survivor.
M = [
 ('L00', 'control', "def u_masses(d):", "def u_masses(d):", 0),
 ('L01', "σ_k without 1/ρ_{j_k} (the theorem's old statement)", "return -s_k * dev_j / (rho_j * v_k)", "return -s_k * dev_j / v_k", 0),
 ('L02', "σ_k's sign", "return -s_k * dev_j / (rho_j * v_k)", "return s_k * dev_j / (rho_j * v_k)", 0),
 ('L03', "the masses uniform, not log-uniform", "return np.exp(rng.uniform(np.log(U_BRANGE[0]), np.log(U_BRANGE[1]), d))",
  "return rng.uniform(U_BRANGE[0], U_BRANGE[1], d)", 0),
 ('L04', "the walls at ±1 (ρ ignored in the path's duals)", "U[k][act] = np.sign(U[k][act]) * rho[act]", "U[k][act] = np.sign(U[k][act])", 0),
 ('L05', "the first entrant not pinned at its wall", "for k in range(1, len(T)):", "for k in range(2, len(T)):", 0),
 ('L06', "the walls ρ_j = b_j^{1/2}", "rho = b ** -0.5", "rho = b ** 0.5", 0),
 ('L07', "the margin 1 − |u_j| (not over ρ_j)", "1 - abs(P['U'][k][j]) / rho[j] for k", "1 - abs(P['U'][k][j]) for k", 0),
 ('L08', "v_k the absolute rate (not over ρ_{j_k})", "@ r_pre) / nx / rho[j])", "@ r_pre) / nx)", 0),
 ('L09', "the pulls with ρ = 1", "rho=float(rho[j]), xi0=x0", "rho=1.0, xi0=x0", 0),
 ('L10', "σ_k computed without 1/ρ", "sig = entry_shift(dev[j], s, v, float(rho[j]))", "sig = entry_shift(dev[j], s, v, 1.0)", 0),
 ('L11', "the 'without 1/ρ' column with it", "sig_w = entry_shift(dev[j], s, v, 1.0)", "sig_w = entry_shift(dev[j], s, v, float(rho[j]))", 0),
 ('L12', "the tracers' pulls in relative units", "dev[jj] += sj * float(rho[j]) * p_hi[jj] + sig * (a0 - a1)",
  "dev[jj] += sj * p_hi[jj] + sig * (a0 - a1)", 0),
 ('L13', "the recursion with the shift without 1/ρ", "dev[jj] += sj * float(rho[j]) * p_hi[jj] + sig * (a0 - a1)",
  "dev[jj] += sj * float(rho[j]) * p_hi[jj] + sig_w * (a0 - a1)", 0),
 ('L14', "the last midpoint a full interval on", "[tks[-1] + 0.5 * (tks[-1] - tks[-2])]", "[tks[-1] + (tks[-1] - tks[-2])]", 0),
 ('L15', "the path's dual drifted from the previous entry", "Up = P['U'][k] + (tm - T[k])", "Up = P['U'][k] + (tm - T[k - 1])", 0),
 ('L16', "W's start at C_W·λ_k", "C_W / C['lam'][q] * se", "C_W * C['lam'][q] * se", 0),
 ('L17', "W's start with the previous interval's λ", "C_W / C['lam'][q] * se", "C_W / C['lam'][q - 1] * se", 0),
 ('L18', "W's windows end at τ_{k+1} − ε^{1/2}", "b_ = tk[q + 1] - eps ** 0.25 if q + 1 < len(tk) else t_end",
  "b_ = tk[q + 1] - eps ** 0.5 if q + 1 < len(tk) else t_end", 1),
 ('L19', "W's block on the coefficients alone",
  "float(max(np.abs(lp['w'][i][a] - R['hist'][k][2][a]).max(), np.abs(lp['c'][i][a] - R['hist'][k][1][a]).max()))",
  "float(np.abs(lp['w'][i][a] - R['hist'][k][2][a]).max())", 0),
 ('L20', "part_W on the first seven targets", "[by[n] for n in chosen(B)]", "[by[n] for n in chosen(B)[:-1]]", 0),
 ('L21', "W's slopes down to ε = 10⁻¹⁰", "if r['eps'] >= 1e-9 and q in dict(r['wins'])", "if r['eps'] >= 1e-10 and q in dict(r['wins'])", 0),
 ('L22', "A's transient at the fastest interval present", "C_B * min(lam[q] for q in present)))", "C_B * max(lam[q] for q in present)))", 0),
 ('L23', "'present' with the end ε^{1/2}", "(tk[q + 1] - e ** 0.25 if q + 1 < nq", "(tk[q + 1] - e ** 0.5 if q + 1 < nq", 0),
 ('L24', "U's σ_k printed as −0", "{e['sigma'] + 0.0:.5g}", "{e['sigma']:.5g}", 0),
 ('L25', "U's crossings against s_θ − σ_k", "abs(c - e['scross'] - e['sigma'])", "abs(c - e['scross'] + e['sigma'])", 0),
 ('L26', "U's second column against σ_k", "abs(c - e['scross'] - e['sigma_wrong'])", "abs(c - e['scross'] - e['sigma'])", 0),
]
