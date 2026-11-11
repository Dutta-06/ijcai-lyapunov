# (id, what, old, new, k): single-point mutants of exp_gauge_plateau.py — Part G, the bounded-domain escape added with the runs declined after
# round 143 (`dom_phi`, `domain_job`: the barrier on the boundary of the rounded box and the deterministic flow back; `escape_runs`: the
# batched dithered runs and their exits; `arrhenius`, `part_G`, `report_G`) and its constants (the domain's L and quartic, the τ, the caps,
# the seeds). G's domain, runs and statistics are stored in the raw, so the md's regeneration sees only `report_G`.
# Run: tools/mutants.py --exp exp_gauge_plateau --tests tests/test_gauge_plateau.py --raw results/chaos/gauge_plateau_raw.npz
# --md results/chaos/gauge_plateau.md --out-env GP_OUT --workers-env GP_WORKERS (one worker; ~2–20 s per mutant).
# Result: 25 mutants and a control. Against the tests as committed: 6 killed by the tests (G03, G04, G18–G20, G22), 4 by the tests and the md
# (G01, G02: the md prints the domain's L and p; G24, G25); the 15 others (G05–G17, G21, G23) killed by the tests added then:
# test_part_G_is_recomputed_in_full_from_the_stored_runs, test_domain_job_reproduces_the_stored_domain (the barrier from the stored set-up,
# the other 198 L-BFGS starts stubbed so that the 63 random boundary points are the stored ones; the flow back over the full S_FLOW) and
# test_escape_runs_reproduce_the_stored_exits (seed 716 at the two largest τ, its stored exit times exactly). No survivor.
M = [
 ('G00', 'control', "def dom_phi(Th, Th0):", "def dom_phi(Th, Th0):", 0),
 ('G01', "the domain's L_α = 3", "DOM_L, DOM_P, N_DOM, N_FLOW, S_FLOW = (2.5, 0.8, 1.2), 4,", "DOM_L, DOM_P, N_DOM, N_FLOW, S_FLOW = (3.0, 0.8, 1.2), 4,", 0),
 ('G02', "the domain a quadratic (p = 2)", "DOM_L, DOM_P, N_DOM, N_FLOW, S_FLOW = (2.5, 0.8, 1.2), 4,", "DOM_L, DOM_P, N_DOM, N_FLOW, S_FLOW = (2.5, 0.8, 1.2), 2,", 0),
 ('G03', "the boundary function without the power", "/ np.asarray(DOM_L)[None, None, :]) ** DOM_P, axis=(1, 2))", "/ np.asarray(DOM_L)[None, None, :]) ** 2, axis=(1, 2))", 0),
 ('G04', "L applied by class, not by coordinate", "np.asarray(DOM_L)[None, None, :]) ** DOM_P", "np.asarray(DOM_L)[None, :, None]) ** DOM_P", 0),
 ('G05', "the boundary point's ℓ_p norm without the 1/p root", "u / np.sum(np.abs(u) ** DOM_P) ** (1 / DOM_P)).reshape(3, 3)", "u / np.sum(np.abs(u) ** DOM_P)).reshape(3, 3)", 0),
 ('G06', "the boundary point with L by class", "Lv = np.tile(np.asarray(DOM_L, float), 3)", "Lv = np.repeat(np.asarray(DOM_L, float), 3)", 0),
 ('G07', "the barrier's starts from seed 1701", "rng = np.random.default_rng(1700)", "rng = np.random.default_rng(1701)", 0),
 ('G08', "ΔJ_𝒪 in units of τ", "return dict(dJ=best.fun * TAU, at=Tb - Th0,", "return dict(dJ=best.fun, at=Tb - Th0,", 0),
 ('G09', "the flow's distance with G, not G⁻¹", "dist = np.sqrt(np.sum((T - Th0[None]) ** 2 / E.G[None, None, :], axis=(1, 2)))",
  "dist = np.sqrt(np.sum((T - Th0[None]) ** 2 * E.G[None, None, :], axis=(1, 2)))", 0),
 ('G10', "the flow half as long", "for _ in range(int(round(S_FLOW / DT))):", "for _ in range(int(round(S_FLOW / DT / 2))):", 0),
 ('G11', "the escapes' noise √(τ dt G)", "sq=np.sqrt(2 * tau * DT * E.G), Th0=Th0, lb0=np.zeros(3))", "sq=np.sqrt(tau * DT * E.G), Th0=Th0, lb0=np.zeros(3))", 0),
 ('G12', "an exit at Σ(·)^p ≥ 1.1", "inside=lambda T: dom_phi(T, Th0) < 1.0,", "inside=lambda T: dom_phi(T, Th0) < 1.1,", 0),
 ('G13', "the exit time one step late", "t=(None if st[r['key']]['exit'] is None else st[r['key']]['exit'] * DT))",
  "t=(None if st[r['key']]['exit'] is None else (st[r['key']]['exit'] + 1) * DT))", 0),
 ('G14', "the runs at γ = 10⁻²", "loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)\n    Th0 = E.exact_rep(net, np.ones(3) / 3)\n    runs",
  "loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)\n    Th0 = E.exact_rep(net, np.ones(3) / 3)\n    runs", 0),
 ('G15', "the rate over the exits' time only", "t = np.array([r['t'] if r['t'] is not None else r['cap'] for r in rr])",
  "t = np.array([r['t'] if r['t'] is not None else 0.0 for r in rr])", 0),
 ('G16', "the rate's s.e. over n", "rate_se=rate / np.sqrt(max(n, 1))", "rate_se=rate / max(n, 1)", 1),
 ('G17', "the scaled rate with e^{−ΔJ/τ}", "scaled=rate * np.exp(dom['dJ'] / tau)", "scaled=rate * np.exp(-dom['dJ'] / tau)", 0),
 ('G18', "the fits unweighted", "w = np.array([r['n'] for r in ok], float)", "w = np.ones(len(ok))", 0),
 ('G19', "the prefactor's sign reversed", "yy = y + c * lt if c is not None else y", "yy = y - c * lt if c is not None else y", 0),
 ('G20', "the free fit's log τ column +", "X = np.c_[np.ones_like(x), -x] if c is not None else np.c_[np.ones_like(x), -x, -lt]",
  "X = np.c_[np.ones_like(x), -x] if c is not None else np.c_[np.ones_like(x), -x, lt]", 0),
 ('G21', "χ² unweighted", "chi2=float(w @ res ** 2)", "chi2=float(np.sum(res ** 2))", 0),
 ('G22', "the extrapolation to τ = 10⁻³", "x4 = np.r_[1.0, -1e4] if c is not None else", "x4 = np.r_[1.0, -1e3] if c is not None else", 0),
 ('G23', "the fits' prefactors τ⁰, τ^{−1}, τ^{−2}", "fits = [arrhenius(rows, dom['dJ'], c) for c in (0.0, 0.5, 1.0)]", "fits = [arrhenius(rows, dom['dJ'], c) for c in (0.0, 1.0, 2.0)]", 0),
 ('G24', "the md's exit time over 1/τ", "{f['t4'] * GAMMA_REF * 1e-4:.2f} |", "{f['t4'] * 1e-4:.2f} |", 0),
 ('G25', "the md's ΔJ_𝒪/τ at one decimal", "ΔJ_𝒪/τ = {dom['dJ'] / TAU:.3f} at τ", "ΔJ_𝒪/τ = {dom['dJ'] / TAU:.1f} at τ", 0),
]
