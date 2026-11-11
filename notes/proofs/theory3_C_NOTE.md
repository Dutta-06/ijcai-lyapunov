# NOTE (theory3, direction C): the interval verification of Theorem T1 without floating-point inputs

Status words: PROVED (argument written, each step checked), COMPUTER-ASSISTED (what is computed, in which arithmetic), OUTLINE, FAILED.
Code: `experiments/theory3c_limit.py` (limit object), `theory3c_rigorous.py` (verifier), `theory3c_run.py` (driver and report); tests `tests/theory3c_test.py`; mutants `tools/mutants/theory3c_{limit,rigorous}_own.py`; report `results/synthetic/theory3c.md` (+ `_raw.npz`, byte-identical regeneration with `--from-raw`).

## 0. What the referee asked (REVIEW_roundT1b.md findings 2 and 5) and what was done

Finding 2: the limit object's data (block critical points, entry times, duals at breakpoints, rates p_e, margins, entrant/feature assignment) were floating-point numbers (hybr) used as exact; the hypotheses r > d_0, "block critical point in D-hat", (a), (b) were checked in floating point.
Finding 5: guards replacing an undecided ball by 0; `np.linspace` edges not outward rounded; `float(ball.upper())` comparisons.

All of it is replaced by ball arithmetic (python-flint 0.9.0 / arb, 100 bits for S(eps), 220 bits for the limit object):

| item | before | now |
|---|---|---|
| target (centres, amplitudes, dictionary) | doubles | the exact binary numbers of the doubles (the same objects, no rounding) |
| block critical point (w, c), k = 1..K | hybr, tol 1e-14, used as exact | Krawczyk operator K(X) = z0 - C F(z0) + (I - C J(X))(X - z0) on the reduced critical system F = (g, h) = 0 (closed-form Jacobian, tested against finite differences); K(X) in int X proves existence AND uniqueness of a root in X; X is then contracted by X <- K(X) cap X. Final box radius about 1e-65 (220 bits). hybr + 6 Newton steps in arb only PLACE the box |
| uniqueness in the whole relaxation set | by hand (Gershgorin) | additionally: Krawczyk on D-hat itself (a box containing it, for every stage): K(D-hat) in int D-hat. One critical point of each stage block lies in D-hat, computed |
| entry times tau_k, duals u_j at the breakpoints, rates p_j, margins | float recursion | the recursion on the enclosed blocks in balls: dt_j = (1 - sgn(p_j) u_j)/\|p_j\| (sign of p_j decided; a rate whose enclosure contains 0 only gets the lower bound (1-\|u_j\|)/sup\|p_j\|) |
| entrant | argmin of floats | STRICT argmin: upper end of the entrant's time < lower end of every other time (or lower bound); otherwise `limit_rig` returns ok = False |
| nearest-unfitted-feature assignment | argmin of floats | strict (balls of the distances; ties refused) |
| sign test sgn p_e = sgn A_f | float, in `cert` | on the enclosures, inside `verify_rig` (`cond['sign']`) |
| margin \|u_j\| < 1 (Thm K (b)) | float | strict on the balls (and implied by the strict argmin) |
| r > d_0 | float, before `verify` | on the balls d_0 = \|c0_e - c*_f\|, certain comparison, in `verify_rig` |
| block critical point in D-hat | only the final block, floats | EVERY stage k = 1..K, every unit: offset <= d-hat, \|w\| in [v_D^-, v_D^+] a_min, sign = sgn A_f, certain comparisons on the enclosures |
| theta < v_# | in `verify` (float) | certain comparison against the ball theta = min(0.01 a_max, 0.1 a_min) (the double 0.01 and the rational 1/100 both inside) |
| r < sqrt 2 (3.0) | not asserted | asserted in `rig_box` |
| all final comparisons | `float(ball.upper())` | certain arb comparisons (`a < b` is True only if it holds for all points of both balls); no `float()` of a ball in any decision or edge |
| undecided guards | partly | `pos_up`, `logpos_up`, `bmax_up` (use the upper end), `clamp01` (corrected: an entirely negative ball clamps to 0, it used to give [-1, 0]), `sq` (nonnegative) |
| quadrature edges (`np.linspace`, `geomspace`) | float, not outward rounded | `fdown`/`fup` (next double below/above): the subintervals cover the integration range; the weights (b - a) and the partition of [0,1] are exact |
| bisection inversions (`x_of_s`, `x_of_w`) | `float(lo(.)) > ...` | certain comparisons |
| eps-ball ends | `10.0**lo` | ball between the double below 10^lo and the double above 10^hi |
| horizon tau_end = 1.05 tau_K + 1 | float | an upper double of the ball |
| bootstrap | E <= D, or E < D/1.0001 in floats | E * (10001/10000) < D, ball comparison |
| window-start gating of the drift terms | gate `up(ta) < lo(wsi) < lo(tb)` on balls, sorted by midpoints | structural (see section 3) |

## 1. The reduced critical system (PROVED, tested)

For k units (w_i, c_i) and all K features in the residual, F(y) = (1/2)||sum_m A_m phi*_m - sum_i w_i phi_i||^2 in the Gram G(x) = exp(-x^2/4):
g_i = sum_m A_m G(c_i - c*_m) - sum_j w_j G(c_i - c_j) = -dF/dw_i, and  w_i h_i = -dF/dc_i, with h_i = sum_m A_m G'(c_i - c*_m) - sum_{j != i} w_j G'(c_i - c_j).
The code's critical equations (g, w h) = 0 and the reduced system (g, h) = 0 have the same roots with all w_i != 0. Jacobian in (w, c):
dg_i/dw_j = -G(c_i - c_j) (j != i), -1 (j = i); dg_i/dc_i = h_i; dg_i/dc_j = w_j G'(c_i - c_j); dh_i/dw_j = -G'(c_i - c_j), dh_i/dw_i = 0; dh_i/dc_i = sum_m A_m G''(c_i - c*_m) - sum_{j != i} w_j G''(c_i - c_j); dh_i/dc_j = w_j G''(c_i - c_j).
(Test `test_closed_form_jacobian_matches_finite_differences`.)

## 2. Krawczyk (PROVED: standard theorem; COMPUTER-ASSISTED: the application)

Theorem (Krawczyk 1969, Moore): F in C^1(X), X a box, z0 in X, C any matrix, J(X) an interval extension of F' on X. If K(X) := z0 - C F(z0) + (I - C J(X))(X - z0) lies in the interior of X, then F has exactly one root in X, and any root in X lies in K(X) (so X <- K(X) cap X keeps the root).
The matrix C is the inverse of the midpoint Jacobian; its accuracy only affects whether the test passes, not its validity. Every operation is in arb. For each of the 8 targets and every stage the test passed at the first box tried (radius 1e-14 around the polished root; checked for all 14 blocks); the contraction reaches radius about 1e-65.
Negative controls (tests): a box of the same radius centred 1e-6 (or 3e-5) off the root in any coordinate is rejected; a guess 1e-3 off with unpolished tiny boxes raises `LimitFail`; a box containing three roots (one unit, two features) is not certified; a box without any root is not certified.

Relaxation-set uniqueness (COMPUTER-ASSISTED): X = D-hat-box = {c_i in c*_f(i) +- up(d-hat), w_i in sign * [lo(a_min v_D^-), up(a_min v_D^+)]}. K(X) in int X was verified for every stage of the 8 targets (`unique_in_Dhat` column of the report). Hence "the block critical point lies in D-hat" (hypothesis of Lemma 4.1) is not only checked for the enclosed point but also unique in D-hat, which is exactly what the hand proof of Lemma 4.1 uses (the target of the Lyapunov argument is THE critical point in D-hat).

## 3. The window-start gating (rewritten, PROVED given the check 'order')

Drift terms of a dormant kernel add the forcing of entrant i only after its window opens at t >= wstart_i = t_i - S^-_i (t_i = lambda tau_i). The old code decided "split inside the piece" by a gate on balls (`up(ta) < lo(wsi) < lo(tb)`), sorted split points by midpoint, and lost the dependency between wsi = tk[i] - Sm[i] and the piece end tk[i] (with Sm = 0 the gate failed to be decided). Now:
- by the condition 'order' (which `verify_rig` checks on every call: tk[0] >= S^-_1 and lambda(tau_{k+1} - tau_k) >= S^-_{k+1} + S^+_k + T_lay), window i opens inside the stage-i interval [t_{i-1}, t_i]; pieces of stages > i see it open, stages < i see it closed. So only the piece(s) of stage st can contain a window start, and only window st.
- inside a stage-st piece the start is placed by its distance back from the stage end: dback = Sm[st] - back (back = lambda * (stage end - piece end), an expression with no eps-dependency loss); fraction f = clamp01(1 - dback/(lambda dtau)). Cases: dback surely > piece length: window open throughout; dback surely < 0: closed throughout; else two sub-pieces (windows < st) on [0,f], (windows <= st) on [f,1].
- The sub-piece forcing counts the windows that may be open (a superset). Sub-pieces of length exactly 0 are skipped.
Comparison with the float-limit verifier on 4 (target, eps) cases of the tests: Etot within 6 percent (typically 3e-5 relative at the same Dv), dcmax between 1x and 4x (more conservative), all other ingredients (T_1, T_2, T_3, Atr, pre, kappa, mass, omega) identical to 1e-3 or better. Golden values stored in the tests.

## 4. Results (generated: results/synthetic/theory3c.md)

8 of 8 targets close: S(eps) verified in ball arithmetic with the enclosed limit object for EVERY real eps in [10^-24, 10^top], and every top is the SAME as with floating-point limit data:
(K=1 seeds 1, 2, 4): -7.5, -5.5, -6.5; (K=2: 6/1/11, 8/1/13, 8/1/22, 6/2.4/12): -7.5, -8.0, -7.0, -7.0; (K=3, 8/1/31): -8.5.
Nothing failed to enclose: all 14 block critical points (stage blocks) pass Krawczyk at the first box; all 14 entry decisions are strict (relative gap to the second hitting time 0.0155 ... 1.30); all signs agree; r - max d_0 >= 0.099; d-hat - max offset >= 0.042; \|w\| margins to [v_D^-, v_D^+] a_min >= 0.049; v_# - theta >= 0.30; Krawczyk on D-hat passes for every stage.
The two targets of Part Q with incomplete chains (rows 7 and 9) were not rerun: their outcome depends on the budget.

## 5. Gaps (stated, not hidden)

1. Theorem K's statement (a), that the block critical point is the global minimiser of the sign-constrained fit, is NOT enclosed. What T1 uses, and what is proved here, is that each stage block has exactly one critical point in D-hat and that it is the enclosed one; the limit object is the recursion with these critical points. (If one defines the limit object by global minimisation, the identification with the D-hat critical point is a separate global-optimisation statement; the paper's wording should define the limit object through the D-hat critical point or state (a) as an assumption.)
2. S(eps) => theorem is the hand proof of THEORY_T1.md (unchanged, found sound by the referee). The ball implementation is a transcription of that proof; its faithfulness is checked by tests against the float implementation (agreement to 1e-3 or the stated conservative differences) and by the referee's independent transcription (0.13 percent), not formally verified.
3. Candidate bounds D_j come from the float certificate; only the verification of E_j(D) < D_j/1.0001 counts, so this is harmless (as before).
4. Soundness is relative to python-flint/arb 0.9.0 (including `arb.__lt__` etc. returning True only if certain, `upper()`/`lower()` rounding, transcendental enclosures) and to the correctness of the Python code above (tests and mutants, section 6).
5. The model constants sigma = 0.12 (mu_c = G m / sigma^2 = 625/9 * 12), the law (12, 1, 1/10, 1) and the formed threshold 0.01 are parameters of the model, not enclosed quantities.
6. Scope unchanged: K <= 3, n = 1, equal masses, continuum Gram, noiseless; the planar case (double precision `cert2d`) and K = 4..6 (float certificate) are NOT in ball arithmetic. Monotonicity of S in eps is not proved: the chains are explicit balls covering [10^-24, 10^top].
7. The upper endpoints 10^top: the covers start at the old float eps_0 (minus a half-decade); no attempt was made to push the tops up beyond the old ones.

## 6. Tests and mutants

`tests/theory3c_test.py` (38 tests, about 40 s). Mutants (tools/mutants/theory3c_limit_own.py, theory3c_rigorous_own.py): limit object 44 mutants, 39 killed; verifier 69 mutants, 52 killed; the survivors are listed with reasons (dead or implied branches, one-ulp changes masked by the padding of arb unions) in the headers of the spec files.
