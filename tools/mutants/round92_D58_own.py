# (id, what, old, new, k): single-point mutants of exp_width_mismatch.py — round 92's new code (a6_lines' split by frame, the two fits'
# difference, the Fock and Mercer-zero thresholds below D‴.61's r*).
M = [
 ('X00', 'control', 'D61_RSTAR = 0.8052 ', 'D61_RSTAR = 0.8052 ', 0),
 ('F01', 'the frame test by |ρ| ≤ 1', "< 1 - 1e-12 if 'v' in z[0] else False", "<= 1 + 1e-12 if 'v' in z[0] else False", 0),
 ('F02', 'the two fits differ by the first only', "{max(abs(z[1] - z[4]) for z in dev):.1e}", "{max(abs(z[1]) for z in dev):.1e}", 0),
 ('F03', 'the bracket count strict', "the two bracketing 1/t(v) at {sum(1 for z in zs if z[1] * z[4] <= 0)}", "the two bracketing 1/t(v) at {sum(1 for z in zs if z[1] * z[4] < 0)}", 0),
 ('F04', 'Fock with τ²/(1 + σ*²)', "fock = lambda r_: rate_conj(r_, 2 * r_ * r_ / (2 + r_ * r_))[0]", "fock = lambda r_: rate_conj(r_, r_ * r_ / (2 + r_ * r_))[0]", 0),
 ('F05', 'Mercer-zero with (1 + τ²)', "mz = lambda r_: rate_conj(r_, r_ * r_ / math.sqrt(1 + 2 * r_ * r_))[0]", "mz = lambda r_: rate_conj(r_, r_ * r_ / math.sqrt(1 + r_ * r_))[0]", 0),
 ('F06', "the ends' threshold without the root", "{math.sqrt((1 + D61_RSTAR) / (1 - D61_RSTAR)):.3f}", "{(1 + D61_RSTAR) / (1 - D61_RSTAR):.3f}", 0),
]
