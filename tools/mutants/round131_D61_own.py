# (id, what, old, new, k): single-point mutants of round 131's new code (referee D3dz) — exp_gh_peano.py's VA011 part, `va_complete`,
# `va_seams`' defaults and completeness, `va011_seam`, `va011_lines`, the eps0 threading of `part_VA`; gh_lastann.py's tail floor,
# `validate`'s eps0 and `validate_tol`, the outer piece's log and the CLI's piece names. A spec names no target file: `tools/mutants.py`
# mutates experiments/<--exp>.py, so this spec picks its list from --exp — run it twice (round 131's run: U01–U16, U18–U26 and U31–U39
# killed; U17 survives, equivalent: with va_complete's R1 = fl(2 − ε₀), 2 − ε₀ ≤ R_bb gives fl(2 − ε₀) ≤ fl(R_bb) = R_bb, rounding being
# monotone, so the R1 check is implied):
#   venv/bin/python tools/mutants.py --exp exp_gh_peano --tests tests/test_gh_peano.py --raw results/chaos/gh_peano_raw.npz \
#       --md results/chaos/gh_peano.md --out-env GP_OUT --spec tools/mutants/round131_D61_own.py --work <scratch>/mut --slot a
#   ... the same with --exp gh_lastann --slot b (gh_lastann writes no md: only the tests decide; its md column is meaningless).
import sys

_exp = sys.argv[sys.argv.index('--exp') + 1] if '--exp' in sys.argv else 'exp_gh_peano'

M_GP = [
 ('U00', 'control', "def va_complete(VA):", "def va_complete(VA):", 0),
 ('U01', "V's ball_end default back to round 127's −1.0", "b.get('ball_end', 0.0) < 0 and b['R1'] <= R_v:", "b.get('ball_end', -1.0) < 0 and b['R1'] <= R_v:", 0),
 ('U02', "va_seams without the completeness", "R_v == R1 and va_complete(VA)[0])", "R_v == R1)", 0),
 ('U03', "va_seams on the numpy array (VA or [])", "    VA = list(VA) if VA is not None else []                           # round 131",
  "    VA = list(VA or [])                           # round 131", 0),
 ('U04', "zero boxes accepted", "if not (g(k, 'boxes', 0) > 0 and g(k, 'failed', 1) == 0", "if not (g(k, 'boxes', 0) >= 0 and g(k, 'failed', 1) == 0", 0),
 ('U05', "the pieces' largest bound unchecked", "and g(k, 'worst', 0.0) < 0):", "and g(k, 'worst', 0.0) < 1):", 0),
 ('U06', "the outer piece's segments unchecked", "    if not g('O', 'segs', 0) > 0:", "    if False:", 0),
 ('U07', "the end piece's bound unchecked (ok only)", "if not (g('end', 'ok', False) and g('end', 'bound', 0.0) < 0):", "if not (g('end', 'ok', False)):", 0),
 ('U08', "the handoff's sub-boxes may be 0", "if not (h.get('sub_boxes', 0) > 0 and", "if not (h.get('sub_boxes', 0) >= 0 and", 0),
 ('U09', "one connector per outer η-box unchecked", "h.get('boxes') == g('O', 'boxes') and", "True and", 0),
 ('U10', "the soundness excess unchecked", "s.get('worst', 1.0) <= 1e-30):", "s.get('worst', 1.0) <= 1.0):", 0),
 ('U11', "R1 against meta's ε₀ unchecked", "    if e0 is None or g('O', 'R1') != 2 - e0:", "    if e0 is None:", 0),
 ('U12', "R1 against meta's ε₀ within 1e-6", "    if e0 is None or g('O', 'R1') != 2 - e0:", "    if e0 is None or abs(g('O', 'R1') - (2 - e0)) > 1e-6:", 0),
 ('U13', "the validation's bar ε₀-free (1e-10)", "for k in ('A_vs_centred', 'A_vs_def', 'B_vs_centred')) < tol", "for k in ('A_vs_centred', 'A_vs_def', 'B_vs_centred')) < 1e-10", 0),
 ('U14', "the closed forms' bar dropped", "            and max(v.get(k, math.inf) for k in ('A0', 'B0')) < 1e-40):", "            ):", 0),
 ('U15', "the validation's points unchecked", "if not (v.get('counts', {}).get('points', 0) > 0 and tol is not None", "if not (tol is not None", 0),
 ('U16', "VA011's seam on R1 only (2 − ε₀ unchecked)", "and Fraction(2) - Fraction(e0) <= Fraction(R_bb) and Fraction(R1) <= Fraction(R_bb))",
  "and Fraction(R1) <= Fraction(R_bb))", 0),
 ('U17', "VA011's seam on 2 − ε₀ only (R1 unchecked)", "and Fraction(2) - Fraction(e0) <= Fraction(R_bb) and Fraction(R1) <= Fraction(R_bb))",
  "and Fraction(2) - Fraction(e0) <= Fraction(R_bb))", 0),
 ('U18', "VA011's seam without the completeness", "    ok = (R_bb > 1.0 and R1 is not None and e0 is not None and va_complete(VA011)[0]",
  "    ok = (R_bb > 1.0 and R1 is not None and e0 is not None", 0),
 ('U19', "va011_lines ignores the seam", "    if seam:\n        out[-1] += (f\". VA011 certifies", "    if True:\n        out[-1] += (f\". VA011 certifies", 0),
 ('U20', "part_VA's validation at the default ε₀", "LA.validate(log=lambda s: None, eps0=e0)", "LA.validate(log=lambda s: None)", 0),
 ('U21', "part_VA's outer piece at the default ε₀", "LA.piece_outer(eps0=e0)", "LA.piece_outer()", 0),
 ('U22', "part_VA's meta ε₀ the default", "delta0=LA.DELTA0, eps0=e0))", "delta0=LA.DELTA0, eps0=LA.EPS0))", 0),
 ('U23', "VA011 at ε₀ = 0.1", "VA011_EPS0 = 0.11 ", "VA011_EPS0 = 0.1 ", 0),
 ('U24', "report without VA011's line", " + va011_lines(VA011, D, V) +", " +", 0),
 ('U25', "va_lines' conclusion without the completeness", "connector `outer_box` reported and an area sum; they now rebuild the cover from the pieces certified (above).\" if ok and complete",
  "connector `outer_box` reported and an area sum; they now rebuild the cover from the pieces certified (above).\" if ok", 0),
 ('U26', "D's ball_worst default lax in va011_seam", "x.get('open', 1) == 0 and x.get('ball_worst', 0.0) < 0]", "x.get('open', 1) == 0 and x.get('ball_worst', -1.0) < 0]", 1),
]

M_LA = [
 ('U30', 'control', "def validate_tol(eps0):", "def validate_tol(eps0):", 0),
 ('U31', "no tail floor", "TAIL_FLOOR = 1e-300 ", "TAIL_FLOOR = 0.0 ", 0),
 ('U32', "the derivative's tail not floored", "(TAIL_FLOOR if t1 < TAIL_FLOOR else t1)", "t1", 0),
 ('U33', "a NaN tail floored away", "(TAIL_FLOOR if t0 < TAIL_FLOOR else t0)", "(t0 if t0 >= TAIL_FLOOR else TAIL_FLOOR)", 0),
 ('U34', "the outer piece's log with the joints (round 131's crash)", "log(json.dumps({k: v for k, v in out.items() if k != 'joints'}))",
  "log(json.dumps(out))", 0),
 ('U35', "validate ignores its eps0", "    eps0 = EPS0 if eps0 is None else eps0\n", "    eps0 = EPS0\n", 0),
 ('U36', "validate_tol's step the wrong way", "        if eps0 <= e:\n", "        if eps0 >= e:\n", 0),
 ('U37', "validate_tol beyond the table: a bar of 1", "            return tol\n    return None", "            return tol\n    return 1.0", 0),
 ('U38', "ε₀ = 0.11's bar 1e-12 (below its 1.1e-12)", "(0.11, 1e-10))", "(0.11, 1e-12))", 0),
 ('U39', "an unknown CLI piece a no-op", "    if any(w not in names for w in todo):", "    if False:", 0),
]

M = M_GP if _exp == 'exp_gh_peano' else M_LA
