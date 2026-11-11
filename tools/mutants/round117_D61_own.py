# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — round 117's revision and the referee's survivors X02, X09, X10
# (the (Re σ − 1)₊² branch, the path grid's range), the ψ grid's duplicate, Ψ_φ in balls, the stored ε lists.
M = [
 ('Q00', 'control', "VC_PREC = 320 ", "VC_PREC = 320 ", 0),
 ('Q01', "_cert_E_ball: (Re σ − 1)₊² dropped (X02)", "        C = (rho * rho - 1) / 2 - rho.log() - (a * a if a > 0 else arb(0))\n    elif", "        C = (rho * rho - 1) / 2 - rho.log()\n    elif", 0),
 ('Q02', 'the path grid from 10⁻² (X09)', "    xq = lambda q: np.r_[np.linspace(1e-4, 1.0, 2001),", "    xq = lambda q: np.r_[np.linspace(1e-2, 1.0, 2001),", 0),
 ('Q03', 'the path grid to max(20, …) (X10)', "np.linspace(1.0, max(200.0, 2 + 4 * q * q), 20001)]", "np.linspace(1.0, max(20.0, 2 + 4 * q * q), 20001)]", 0),
 ('Q04', 'ψ = 1 twice', "    pg = np.unique(np.r_[np.linspace(0.02, 1.0, 50), np.linspace(1.0, 20.0, 96)])", "    pg = np.r_[np.linspace(0.02, 1.0, 50), np.linspace(1.0, 20.0, 96)]", 0),
 ('Q05', 'Ψ_φ in balls with 3φ', "            ps = -(xa * xa - ya * ya) - 4 * pa * xa * ya", "            ps = -(xa * xa - ya * ya) - 3 * pa * xa * ya", 0),
 ('Q06', 'the Re σ > 1 point dropped', ", (3.0, 5.0, 3.0), (0.2, 3.0, 0.1))", ", (3.0, 5.0, 3.0))", 0),
 ('Q07', 'the worst rounded away from 0', "{math.trunc(pt['worst'] * 1e4) / 1e4:.4f}", "{math.floor(pt['worst'] * 1e4) / 1e4:.4f}", 0),
]
