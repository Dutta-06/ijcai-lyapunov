# Referee report D3dc, round 108: Proposition D‴.61 as revised after round 104 (V in ball arithmetic)

## Scope

**Commit.** I reviewed d5d8438 on branch claude/vigilant-clarke-hf939b. The branch has since moved to 319e8a5, a D‴.59 commit from another session. That commit touches none of D61.md, exp_gh_peano.py, test_gh_peano.py, gh_peano.md or the raw.

**Setup**
- **Work directory.** W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref108.
- **Export.** I ran `git archive d5d8438 | tar -x -C W/repo` and symlinked the repository's venv (python-flint 0.9.0, numpy 2.5.4, scipy 1.18.1) as W/repo/venv. GP_OUT = W/gp_out holds a copy of the raw.
- **Run conditions.** Every job ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with a `timeout`. Pytest used `--basetemp` inside W (W/pyt; for tools/mutants.py through PYTEST_ADDOPTS, one basetemp per slot).
- **Processes.** At most two heavy processes of mine ran at a time. I stopped no process, mine or anyone else's.
- **Repository untouched.** I used only read-only git commands (archive, log, diff, status, rev-parse). The sha256 of W/repo's exp_gh_peano.py plus test_gh_peano.py equals that of d5d8438's.

**Read in full**
- D61.md.
- exp_gh_peano.py:
  - part V, from `cert_E` to `v_lines`;
  - `ratio_floor` and the CERT_* constants;
  - `d_lines`' ratio lines;
  - D's `gap_closed`, `quad_roots`, `bound_sup`, `box_bound` and `dom_bb`.
- gh_peano.md, parts D and V.
- tests/test_gh_peano.py: the part V tests and the round 104 test in full, the rest skimmed.
- tools/mutants/round104_D61_own.py and tools/mutants.py.
- The Round 104 entry of REVIEW_LOG.md and REVIEW_round104_D3cy.md.
- `git diff b98e652 d5d8438` of these files and of the requirements files.
- HANDOFF.md's and PAPERS.md's D‴.61 lines.

**Web.** Search snippets only. python-flint.readthedocs.io does not resolve from this session (ENOTFOUND); it would have to be allowed in the session's network settings. For python-flint's semantics I used the installed package's docstrings and my own probes.

## Summary

**Design.** Part V certifies the valley condition on the closed annulus 1.8927 ≤ |η| ≤ 1.995, 0 ≤ arg η ≤ π/2.
- **Coordinates.** In the saddle's coordinates σ = (η/2)(1 + z), E(z; η) = F_η(σ) − F_η(η/2) = log|1 − z²| + Re(η²z²/4) + C(σ).
- **The straight piece** from z = −1 to b = t = 0.9(2 − R2)/R2 is negative by hand.
- **The polyline.** One float polyline per η-box, the box centre's descent flow, is bounded segment by segment by E(m; η_c) + ℓL_z + δL_η. Since round 104 the bound uses |∇_σC| ≤ h′(|σ|).
- **The end piece** runs from |1 − z_n| ≤ 0.02 to z = 1 and has a closed-form bound.
- **The ball recheck.** Every accepted segment is now rechecked in Arb at 120 bits (`cert_seg_ball`). So are each box's covering radius (top edge at the true π/2), b's place in the disc, the end check and the end piece. The results are stored in the raw and printed.
- **What the draft concludes.** V is rigorous. (iii)(b) holds on |η| ≤ 1.995, and (iv) holds for r < 0.9900 (σ*/σ < 199.00), "modulo D's floating-point branch and bound alone".

**Verdicts**
- **Reproducibility is exact.**
  - The md regenerates byte for byte (sha256 a7287107…2044).
  - All 35 tests pass (10.2 s).
  - All three bands rerun with the ball recheck are repr-identical to the raw (7.9 s, 31.2 s, 280.4 s on the shared machine).
  - Every number new in D61 is in the md. The one exception is "77 million", which cites round 104's report.
- **The gradient identity is exact.** C is C¹ across |σ| = 1 and |Re σ| = 1, so the mean-value argument is valid on segments that cross them. Every term of `cert_seg_parts` is a valid bound.
- **The ball code is sound as a proof** (details under What I checked, item 3).
  - Every `if` compares exact arbs, or fails in the safe direction. Every monotone term is taken at its rigorous end, and the float result is rounded up correctly.
  - The tiling is exact. The only remaining floating-point dependencies are trivial (D3dc-3).
- **My own second recheck agrees.** I re-bounded all 10,687,833 accepted segments with my own differently built Arb enclosure: 0 failures. I also redid the corners, b, the end, the chain and the exact rational tiling, with no fault.
- **D's branch and bound is correctly called the one remaining floating-point step.** Porting it to ball arithmetic is easy: my 60-line port verifies all 402,664 boxes rigorously in 18.5 s, with the largest rigorous bound −1.9432·10⁻¹¹ (D3dc-1).
- **Ω′ and the proofs.** Ω′ as restated is exactly what the proof uses. The proofs of (iii)(b) and (iv) are correct as written. "Certified in ball arithmetic" and "rigorous" are earned for V. "Modulo the branch and bound alone" is accurate.
- **Tests and mutants.**
  - Round 104's 21 mutants: 20 are killed. B02 survives, equivalent as logged.
  - Of my 27, 18 are killed, including every "cheaper" change to the float pass.
  - The unsound survivors Z18, Z19 and Z30 change the recomputed rigorous end-piece bound, and no test sees it. The threshold mutants Z21 and Z22 are unsound but equivalent on the data (D3dc-2).

No major issue.

**Findings**
- **D3dc-1 (minor, constructive). D's branch and bound in ball arithmetic is an 18-second computation.** My port of `box_bound` to Arb checks every one of the 402,664 boxes on 1 < |η| ≤ 1.8927 rigorously:
  - with R1 exact (no 1 + 10⁻¹⁵) and the top edge at the true π/2;
  - with the critical points of both pieces enclosed by ball roots;
  - result: 0 boxes unresolved, the largest rigorous bound −1.94320·10⁻¹¹, the tiling exact in rationals.

  Shipping it removes "modulo the branch and bound" everywhere.
- **D3dc-2 (minor). The rigorous end piece is not tested against its formula, and two ball thresholds are untested.**
  - Z18 (h dropped from `cert_end_ball`), Z19 ((1 + ε₀)² dropped) and Z30 (`ball_end` at Ra) are unsound and change the recomputed `ball_end`. All three pass the tests.
  - Z21 (accept a rigorous bound below 10⁻¹² instead of below 0) and Z22 (b inside the disc to 1.001) are unsound but equivalent on the data.
- **D3dc-3 (cosmetic). Floating-point residues besides D, all harmless on the data, none asserted:**
  - `v_lines` accepts a contiguity gap of 10⁻¹²;
  - the band grid's last edge Ra + (Rb − Ra)·nR/nR == Rb holds but is not asserted;
  - |η_c| ≤ R2, which L_η needs, is not asserted;
  - r* and c_need are float bisections. The margins are 2.3·10⁻⁵ at 1.995 and 1.8·10⁻⁶ at 1.8927; I verified both in Arb.
- **D3dc-4 (cosmetic). Wording.**
  - "|∇_σC| ≤ h′(|σ|), which is exact": the bound is attained only where |Re σ| ≤ 1 < |σ| and is strict beyond. Say "sharp".
  - New says "proved … on |η| < 1.8927" while (iii)(b) and the Status say "≤", and New omits the limit on |η| ≤ 1.995.
  - The proof's parenthetical "(strictly within 2rb of η/2)" belongs to condition (1). It also holds for (2), but the restated Ω′ needs only "in the open disc".
- **D3dc-5 (cosmetic). Stale text.**
  - PAPERS.md, Paper 18, still has:
    - "σ*/σ < 198.99";
    - "modulo both computer-assisted steps";
    - "Missing. The two computer-assisted steps in interval arithmetic";
    - "about sixfold per halving".
  - HANDOFF.md:
    - §7 (line 243) still says "about sixfold per halving";
    - the experiments table lists `exp_gh_peano.py` parts as "P, M, L, LC, S, R" (Q, D and V missing);
    - the test time is still "about 2 s" (now about 10 s).
  - `cert_seg_parts`' docstring still opens with the old bound h′ + 2(|Re σ| − 1)₊.

## What I checked and how

**(1) Reproducibility** [own][md][raw]
- **The md.** `--from-raw` into W/gp_out gives an md identical to the committed one (cmp; sha256 a728710726d5…2044 for both).
- **Tests.** `pytest -q -p no:cacheprovider --basetemp=W/pyt/t1 tests/test_gh_peano.py`: 35 passed in 10.21 s. The slowest is the round 104 test, at 4.1 s.
- **Part V against the raw.** W/own/rerun_band.py reran `cert_band(*CERT_BANDS[i], ball=True)` with the script's own functions. All three result dicts are `==` to the raw's, for every key (boxes, segments, worst, side, end, ball_worst, ball_excess, ball_short, ball_end):

  | band | boxes | segments | ball_worst | ball_excess | ball_short | time |
  |---|---|---|---|---|---|---|
  | 0 | 760 | 160111 | −2.4767054·10⁻¹² | 2.89·10⁻¹⁴ | 220 | 7.9 s |
  | 1 | 2011 | 946644 | −1.0258402·10⁻¹² | 2.63·10⁻¹⁴ | 495 | 31.2 s |
  | 2 | 10658 | 9581078 | −1.0004003·10⁻¹² | 3.15·10⁻¹⁴ | 2584 | 280.4 s |

- **D's branch and bound.** Its float recompute, (402664, −1.9432233600014115·10⁻¹¹, 0), equals the raw's `dom['bb']`.
- **Numbers.** Applying tools/numcheck.py's logic to `git diff b98e652 d5d8438 -- D61.md`, every added number is in the md except "77", round 104's 77 million.

**(2) The gradient identity and C's regularity** [own] (W/own/identity_grad.py, Arb at 200 bits)
- **By hand.** Let a = |Re σ| and ρ = |σ|.
  - For a > 1, ∇C = h′(ρ)σ/ρ − 2(a − 1)sgn(Re σ)e_x. So |∇C|² = h′² − 4(a − 1)h′a/ρ + 4(a − 1)².
  - With h′ = ρ − 1/ρ this is h′² − 4(a − 1)(1 − a/ρ²).
  - Since a ≤ ρ ≤ ρ² for ρ ≥ 1, the subtracted term is nonnegative.
  - For a ≤ 1 < ρ, |∇C| = h′(ρ) exactly. Inside the disc ∇C = 0, consistent with hp = 0.
- **At 200 bits.** Over 17,527 random σ with |σ| > 1:
  - |(∇C)²computed − identity| ≤ 1.7·10⁻³⁴;
  - max |∇C|/h′ = 1.0 (attained where a < 1).
- **Regularity.** One-sided directional derivatives agree at (0.6, 0.8) (on |σ| = 1), at (1, 0.5) and (1, 1.2) (on |Re σ| = 1), and at (1, 0) and (0, 1). So C is C¹, as h(1) = h′(1) = 0 and (a − 1)₊² is C¹.
  - Hence E(·; η) is C¹ off z = ±1, and E(z; ·) is C¹ for z ≠ ±1.
  - The mean-value inequality therefore holds on segments that cross |σ| = 1 or |Re σ| = 1 (W > 0 excludes z = ±1).
- **The terms of `cert_seg_parts`**, all correct:
  - g = log(1 − z²) + cz², g′ = 2z(c − 1/(1 − z²)) and g″ = 2c − 2(1 + z²)/(1 − z²)²;
  - |1 − z²| ≥ |1 − m²| − ℓ(2|m| + ℓ) and |z| ≤ |m| + ℓ;
  - |∇_zC| = (|η_c|/2)|∇_σC| ≤ (|η_c|/2)h′(|σ_m| + |η_c|ℓ/2), since hp = max(0, r − 1/r) is nondecreasing;
  - L_η = R2·Zmax²/2 + (|1 + m| + ℓ)/2 · h′((R2/2)(|1 + m| + ℓ)).
- **The two legs.** The leg (m, η_c) → (z, η_c) → (z, η) needs [η_c, η] ⊂ {|η′| ≤ R2}. That holds because |η_c| ≈ (R1 + R2)/2 < R2 by at least 2.4·10⁻⁵, which is not asserted (D3dc-3).

**(3) The ball code as a proof** [own] (W/own/flint_semantics.py, W/own/indep_ball.py, W/own/attack.py)
- **python-flint 0.9.0 semantics.**
  - **Conversions.** `arb(float)` is exact: radius 0 for 0.1, 1.8927, 1.995, the float π/2, 0.02, 10⁻³⁰⁰ and 5·10⁻³²⁴.
  - **Bounds.** `.upper()`/`.lower()` return exact arbs. The docstring says "an exact floating-point number that has been rounded up [down] to the current precision" [own, installed docstring]. `.mid()`/`.rad()` return exact arbs.
  - **Comparisons.** For an overlapping ball [0.999, 1.001], `<`, `>`, `<=`, `>=`, `==` and `!=` against 1 are all False, so `not (r > 1)` is True. For exact 1, `>` is False and `<=`/`>=` are True.
  - **NaN.** `log` of a ball containing 0 is nan, and comparisons with nan are False. `1/0` is nan.
  - **float().** `float(arb)` rounds the midpoint to nearest: float(0.1 − 2⁻⁸⁰) = 0.1, which is above. So `math.nextafter(float(u), inf)` is needed for a rigorous float upper end, and it is there. The final result is ≥ u in every case.
- **Every `if` in the ball code** decides on an exact arb or fails in the safe direction:
  - `hp`'s `r > 1`, on an exact upper end;
  - `W > 0`, on an exact lower end;
  - `up(rho) <= 1`, exact;
  - `rho.lower() > 1`, exact;
  - `al > 0`: a straddling al gives 0, the safe lower bound of (·)₊²;
  - `up(h) > 0`: C ≤ max(0, h), and h ≥ 0 anyway;
  - `v > best`, exact;
  - `db > arb(d)`, exact;
  - b's `< 1`, exact;
  - the end's `<= arb(CERT_EPS0)`, exact; arb(0.02) is the double that `cert_end_ball` uses;
  - `u < 0`, on a rounded-up float.

  A nan anywhere makes `u < 0` False, so the segment is rejected or bisected.
- **Monotone ends.** Each term is taken at the correct end:
  - ℓ, Zmax and rmax = |σ_m| + |η_c|ℓ/2 at their upper ends;
  - W at its lower end; Lg, L_z and L_η at their upper ends;
  - opz and rB = (R2/2)·opz at their upper ends;
  - h′ at the exact upper end of its argument;
  - C(σ_m) at h's upper end minus ((lower end of |Re σ_m|) − 1)₊²;
  - E(m) and the total at their upper ends;
  - δ the rigorous corner maximum, or the float δ when that is not smaller.
- **The midpoint ball.** m = (P + Q)/2 is the exact midpoint of the exact double endpoints, so the disc |z − m| ≤ ℓ covers the segment between those doubles.
- **Closing the gaps between the float computation and the claim**
  - **The tiling of the polyline.** Children split at the float midpoint m_f, so the accepted segments form a continuous chain of exact segments from pts[0] = b to pts[−1]. That is all the proof needs; the original polyline is irrelevant. My check found 0 chain breaks in 13,429 boxes.
  - **b.** It is the double t; (R2/2)(1 + t) < 1 is checked in Arb for the box's largest radius, and |σ_b| ≤ 0.99975.
  - **The end piece.** |1 − z_n| ≤ double(0.02) is checked in Arb. `cert_end_ball(Rb)` uses the same double, and cert_end increases in R2 ≥ every box's R2.
  - **The top edge.** It is extended to the true π/2 in δ. The float π/2 is 6.1·10⁻¹⁷ below the true π/2.
- **The box tiling** (W/own/indep_ball.py).
  - The leaves' Σ(R2 − R1)(f2 − f1) equals (Rb − Ra)·(float π/2) exactly in Fractions for each band, and every leaf lies inside its band. Splits share the float midpoint, so the leaves are an exact partition.
  - The grid's last edge equals Rb exactly in all three bands, and the bands share their edges exactly.
  - CERT_BANDS[0][0] == DOM_BB[−1] (the test asserts it).
  - The float 1.995 is decimal 1.995 + 1.07·10⁻¹⁶, and the float 1.8927 is decimal 1.8927 + 4.9·10⁻¹⁷, so the decimal claims are covered.
- **Symmetries, by hand.** E(z̄; η̄) = E(z; η) and E(z; −η) = E(z; η). Moreover ĝ_K is even and real on ℝ, so the limit transfers without them.
- **My own recheck (indep_ball.py, Arb at 128 bits).** It uses the script's polyline and box tree, but re-bounds each accepted segment with a differently built enclosure:
  - g′ by its natural interval extension on an acb rectangle containing the disc |z − m| ≤ ℓ, not |g′(m)| + ℓ·sup|g″|;
  - C's gradient through h′ at the rectangle's |σ| upper end;
  - my own corner enclosure via exp(if), with the top edge at the true π/2.

  Where my bound was not below 0, I bisected at exact midpoints.

  | band | script segments | my bounds | my extra splits | failures | my largest bound | chain, b, end | time |
  |---|---|---|---|---|---|---|---|
  | 0 | 160111 | 161123 | 1012 | 0 | −2.48·10⁻¹² | all OK | 7.8 s |
  | 1 | 946644 | 948742 | 2098 | 0 | −6.33·10⁻¹³ | all OK | 40.1 s |
  | 2 | 9581078 | 9588932 | 7854 | 0 | −1.30·10⁻¹³ | all OK | 401 s |

  My end pieces are −1.382677, −1.357139 and −1.344304. Sampled against E on 62 segments, my bound is above E everywhere, by at least 4.2·10⁻⁵.
- **An attack from F's definition** (attack.py, 128 bits).
  - Inputs: 336 η (300 random, plus 36 at the band edges × {0, 10⁻⁹, 0.0058, 0.08605, 0.345, 0.61756, π/4, π/2 − 10⁻¹², π/2}), giving 360 (η, box) pairs. Points on each path: the straight piece, 8 per polyline segment, 20 on the end piece.
  - Every η lies in a box.
  - The largest F − F(η/2) away from the saddle is −2.54·10⁻⁸, at η = 1.995 on the real axis. By band: −1.65·10⁻⁶, −2.04·10⁻⁷, −2.54·10⁻⁸.
- **Still only floating point.**
  - D's branch and bound (item 4).
  - The trivial residues of D3dc-3.
  - The flow and CERT_SLACK, which are not needed for soundness.

**(4) D's branch and bound** [own] (W/own/dom_bb_ball.py, W/own/dom_sanity.py)
- **The description is right.** `dom_bb` is the one floating-point step left in (iii)(b) on |η| ≤ 1.995.
  - Its box bound is analytically valid. R² cos 2φ is smallest at a corner, R cos φ largest at (R2, f1), and −log R largest at R1.
  - w² < ½ on (0, 1 − 1/R2] (w ≤ 0.4717), so the least p dominates. M is increasing and C¹, so the supremum lies at wm or at a root of one of the two pieces' quadratics: (p + 2B²)w² − 2B²w + 1 and pw² − 2Bw + 1.
- **The port.**
  - It follows the script's float traversal, and every accepted box is re-bounded in Arb.
  - The parameters are exact numbers: R1 exact (R1 = 1 kept), p at the rigorous lower end over the corners with the top edge at the true π/2, B and wm at their upper ends.
  - The roots are enclosed by the stable ball formula; where qa may vanish, |q/qa| > wm is proved. Candidates are wm, the kink and the root balls. h is taken at its upper end with M at the upper end of a. A root ball the code cannot place would mark the box for my own bisection.
- **Results.**
  - 1 < |η| ≤ 1.78: 4935 boxes; rigorous largest bound −3.2812573·10⁻⁷ (float −3.2812573·10⁻⁷).
  - 1 < |η| ≤ 1.8927: all 402,664 boxes rigorous, 0 unresolved; largest rigorous bound −1.9432018·10⁻¹¹ (float −1.9432234·10⁻¹¹); tiling exact in Fractions; 18.5 s.
  - Sanity: on 3000 random boxes (10% touching π/2), against 74,745 dense samples of `dom_gap`, max(sampled gap − ball bound) = −1.57·10⁻⁴.
- **Effort.** Rechecking it "the same way" takes about 60 lines and 18.5 s (D3dc-1).

**(5) Ω′, the proofs of (iii)(b) and (iv), and the wording** [md][own]
- **Ω′'s definition is the proof's hypothesis:** a chord through η/2 in the open disc along u with Re(φ″u²) < 0, below F(η/2) except at η/2, and paths 0 → a and b → η below F(η/2).
- **Ω ⊆ Ω′.** On the straight segment, F_η(ηs) − F_η(η/2) = g_η(s) + C(ηs) − g_η(½). Below 1/|η| this is negative except at ½; beyond, the domination condition makes it negative. The chord runs along u = η/|η|, with Re(φ″u²) = (2Re η² − 8)/|η|² < 0.
- **Condition (2)'s analytic part.** E(x) = log(1 − x²) + Re(c)x² < −x²(1 − Re c) < 0 for 0 < |x| < 1 and |σ| < 1. Its chord has half-length t|η|/2 < 2rb.
- **Condition (1).** It is stated with the descent line and Taylor's bound q < 1. It is unchanged and checked on a grid only, as said.
- **The proof's step near |σ| = 1.** "Λ⁺'s continuity across |σ| = 1 covering the unit circle as above" works because V's margin is strict: max rigorous bound −1.0004·10⁻¹². Then h(r₀) < 10⁻¹² for r₀ = 1 − 10⁻⁶, and (a) is uniform on |σ| ≤ r₀.
- **(iv)'s proof.** It now says "for r < 0.9900, with V" and "inside |η| ≤ 1.995". In Arb:
  - c_need(0.99) = 1.9949769 < 1.995;
  - c_need(0.8052) = 1.8926982 < 1.8927;
  - c_need(0.0966) = 0.9998257 < 1;
  - (1 + 0.99)/(1 − 0.99) = 199 exactly.
- **Wording.**
  - "Certified in ball arithmetic" and "computer-assisted, and rigorous" are now earned for V, by the standard of validated numerics (items 3 and 7; bit-identical rerun; my independent second recheck).
  - "Modulo the branch and bound alone" is accurate, apart from D3dc-3's trivia.
  - New's "|η| < 1.8927" is out of step (D3dc-4).
- **Box concentration** (W/logs/leaves*.pkl): within 5° of the real axis, 71.8%, 90.2% and 96.9% of the bands' boxes; within 1°, 33.6%, 66.4% and 88.7%. So "almost all of them near the real axis" in the Status is right.

**(6) Tests and mutants** [own] (tools/mutants.py, one slot at a time; survivors recomputed on band 0 and on the slice 1.99 ≤ |η| ≤ 1.9902 with ball=True by W/mut/recompute.py, against a recomputed control)
- **Round 104's spec** (`round104_D61_own.py`):
  - 20 of 21 are killed.
  - B02 (ℓ at its lower end) survives. It is equivalent: the change is about 10⁻³⁶ and the final rounding up absorbs it.
  - This matches the log.
- **My spec** (W/mut/spec108.py, 27 mutants).
  - **Killed (18):**

    | id | change |
    |---|---|
    | Z02 | Zmax without ℓ |
    | Z03 | W without ℓ² |
    | Z04 | sup\|g″\| without 2\|c\| |
    | Z05 | ball rmax without \|η\|ℓ/2 |
    | Z06 | ball rB with \|η_c\| |
    | Z10 | no final rounding up |
    | Z11 | the η term dropped |
    | Z13 | ℓ halved |
    | Z15 | the f1 corners only |
    | Z27 | E at the endpoint p |
    | Z33 | ball L_η's \|1 + z\|/4 |
    | Z34 | L_z without C |
    | Z35 | cheaper and sound: the float slack skipped where the ball decides |
    | Z36 | prec 53 |
    | D01 | box p at f1 only |
    | D02 | B at f2 |
    | D04 | DOM_SLACK 10⁻¹⁵, cheaper; tests and md |
    | D05 | wm from R1 |

  - **Unsound survivors that change the recomputed values (3):**

    | mutant | change | recomputed ball_end, slice | band 0 | control |
    |---|---|---|---|---|
    | Z18 | h dropped from `cert_end_ball` | −2.1787 | −2.1892 | −1.3566 and −1.3827 |
    | Z19 | (1 + ε₀)² dropped | −1.3966 | −1.4223 | same |
    | Z30 | `ball_end` at Ra | −1.3571 | −1.5980 | same |

    The tests check only `ball_end < 0` and `< −1.35`.
  - **Unsound but equivalent on the data (5).** Recomputes are identical for all five.
    - Z07, Z08 and Z09: |Re σ_m| at its upper end; C = 0 where |σ_m| straddles 1; the inside test on the lower end. The midpoint balls are exact, so these act below double resolution.
    - Z21: rigorous bounds below 10⁻¹² accepted.
    - Z22: b inside the disc to 1.001.

    The planted tests use u = 1.0 and |σ_b| = 1.0045, which are too far out to tell these thresholds apart.
  - **Sound and equivalent (1):** D03, the kink candidate dropped. The dom recompute is identical, (402664, −1.9432233600014115·10⁻¹¹, 0).
- **Do the pins catch a cheaper change?** Yes for any change to the float pass or the tree (Z35, Z36, D04, round 104's Y-series), through the pinned segment counts and values. Yes for the ball segment formula, through the 60 random segments with |u − ub| ≤ 10⁻¹³ that reach |σ| > 1. No for the end piece's ball twin or the two ball thresholds (D3dc-2).

**(7) Literature** (web search; [snippet] unless marked)
- **Arb.** Midpoint–radius ball arithmetic with automatic rigorous error bounds (Johansson, ARITH 24, 2017; arXiv:1611.02831). The rigorous integrator in ball arithmetic gets its error bounds from complex magnitudes on ellipses, without derivatives (Johansson, arXiv:1802.07942).
- **Mean-value forms.** The mean-value or centred form, f(m) plus an interval gradient times the radius, used on boxes with recursive subdivision, is the standard enclosure. Its overestimation decreases quadratically in the width:
  - Moore–Kearfott–Cloud, Introduction to Interval Analysis, SIAM 2009 [snippet: bibliographic only];
  - Caprani–Madsen 1980 and Krawczyk–Nickel 1982, as cited in search results [snippet];
  - arXiv:1409.6414 (mean-value form evaluated in interval arithmetic) [snippet].

  V's E(m) + ℓL_z + δL_η is exactly such a form.
- **Arb with adaptive bisection in a recent computer-assisted proof:** arXiv:2411.10885 (restricted three-body problem, 256 bits) [snippet].
- **The standard of proof.** Floating point without rigorous error analysis is not a proof; interval and ball arithmetic is. See the ICM 2018 panel (arXiv:1809.08062) and the Gómez-Serrano survey (arXiv:1810.00745) [snippets].
- **No direct precedent.** I found no validated (interval) certificate of a steepest-descent or Laplace step for a quadrature or Krylov error.
  - The nearest is a FAU thesis that combines the saddle-point method with interval arithmetic in MATLAB [snippet, record only].
  - O'Sullivan (arXiv:1702.03777) gives explicit, non-validated saddle-point error control [snippet].

  So the draft's "new" claim for the certified deformed path is not contradicted by what I saw.

## Verdicts per part

- **(i), (ii), (iii)(a), (iii)(c) and (v):** unchanged by this revision, and proved as before. P_j's polynomial form is still a sketch.
- **(iii)(b)**
  - **|η| < 1 and both axes:** proved.
  - **1 < |η| ≤ 1.8927:** proved, modulo D's floating-point branch and bound. My Arb port verifies that branch and bound rigorously (D3dc-1); it is not yet in the repository.
  - **The closed annulus 1.8927 ≤ |η| ≤ 1.995:** proved, by a computer-assisted proof in ball arithmetic.
    - The repository's recheck reran bit-identically.
    - My independent enclosure re-verified all 10.7 M segments, the corners, b, the end, the chain and the tiling.
    - Ω′ as restated matches the proof, and both sufficient conditions are correct.
  - **1.995 < |η| < 2:** a conjecture.
- **(iv):** correct. It holds for r < 0.8052 (σ*/σ < 9.26) and for r < 0.9900 (σ*/σ < 199.00, exact), with the proviso of D alone, which D3dc-1 removes. The constants are right in Arb.

## Status of the round-104 issues

| id | status |
|---|---|
| D3cy-1 (Ω′) | **Resolved.** Ω′ is the proof's hypothesis, with (1) and (2) as sufficient conditions and "the proof applies as soon as…". |
| D3cy-2 ("certified") | **Resolved.** `cert_seg_ball`, `cert_box_disc_ball` (true π/2), `cert_end_ball`, b and the end in Arb are shipped, stored in the raw and printed. python-flint 0.9.0 is in both requirements files. Rerun bit-identical; my independent recheck agrees. |
| D3cy-3 (tests) | **Largely resolved.** Pinned boxes and slice, float and ball. Round 104's spec: 20 of 21 killed (B02 equivalent). Remaining: the end piece's ball twin and the ball thresholds (D3dc-2). |
| D3cy-4 (redundant term) | **Resolved.** The term is dropped and the identity is in the docstring; I verified it to 1.7·10⁻³⁴. The docstring's first sentence is stale (D3dc-5). |
| D3cy-5 (198.99) | **Resolved** in D61, the md and `ratio_floor` (exact Fraction). Still stale in PAPERS.md (D3dc-5). |
| D3cy-6 (stale text) | **(iv)'s proof and HANDOFF's D‴.61 state and test count: fixed.** "Sixfold" remains in HANDOFF §7 and PAPERS.md, and PAPERS.md's Paper 18 predates round 104 (D3dc-5). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dc-1 | minor | D's branch and bound, the step "alone" left in floating point, ports to Arb in about 60 lines. With the script's float tree, R1 exact, p at its rigorous lower end over the corners (top edge at the true π/2), B and wm at their upper ends, and ball roots for both pieces' quadratics: all 402,664 boxes rigorous, 0 unresolved, largest bound −1.9432018·10⁻¹¹ (float −1.9432234·10⁻¹¹), tiling exact, 18.5 s; 4935 boxes at 1.78. Sanity: 74,745 samples on 3000 random boxes, max(gap − bound) = −1.57·10⁻⁴ [own: W/own/dom_bb_ball.py, W/logs/dom_bb_ball.log]. | exp_gh_peano.py: add `box_bound_ball` and run `dom_bb` with ball (every accepted box rechecked, bisected further unless below 0). Store `ball_worst`/`unresolved` in `dom['bb']` and print them in `d_lines` and `v_lines`. Test: pin the 1.78 run in ball, plant an R1 = 1 box and a top box. D61 (iii)(b), (iv), Status, New, Numerics: "verified in ball arithmetic", and drop "modulo the branch and bound". |
| D3dc-2 | minor | The rigorous end piece is untested against its formula. Z18 (h dropped), Z19 ((1 + ε₀)² dropped) and Z30 (`ball_end` at Ra) are unsound, change the recomputed `ball_end` (slice −2.1787, −1.3966, −1.3571 against −1.3566; band 0 −2.1892, −1.4223, −1.5980 against −1.3827) and pass. Z21 (u < 10⁻¹² accepted) and Z22 (b to 1.001) are unsound, equivalent on the data, and pass: the plants (u = 1.0, \|σ_b\| = 1.0045) are too far out. Z07–Z09 act below double resolution. Round 104's 21: 20 killed, B02 equivalent [own: W/mut/spec108.py, W/logs/mut108.log, W/logs/rc_*.log]. | tests/test_gh_peano.py: <br>• cert_end(R) ≤ cert_end_ball(R) ≤ cert_end(R) + 10⁻¹⁴ for R in (1.9, 1.98, 1.995); <br>• pin bb['ball_end'] == −1.3566268758611253; <br>• plant `cert_seg_ball` → 10⁻¹⁵ (a positive tiny bound must be rejected); <br>• plant cert_b with (R2/2)(1 + t) = 1.0005. <br>Record Z07–Z09 and D03 as equivalent. |
| D3dc-3 | cosmetic | Floating-point residues besides D, true on the data but not asserted: `v_lines` accepts a contiguity gap of 10⁻¹²; the grid's last edge Ra + (Rb − Ra)·nR/nR == Rb (it is exact in all three bands); \|η_c\| ≤ R2 (needed by L_η; margin ≥ 2.4·10⁻⁵); r*/c_need by float bisection. In Arb, c_need(0.99) = 1.9949769 < 1.995, c_need(0.8052) = 1.8926982 < 1.8927 (margin 1.8·10⁻⁶) and c_need(0.0966) = 0.9998257 < 1 [own]. | `v_lines`: require `b['R1'] == R_v` (exact). `cert_band`: set the last edge to Rb explicitly. `cert_box` with ball: assert abs(ec) < R2 in Arb. `v_lines`/`d_lines` or a test: c_need at the printed floors in Arb. Then "alone" is literally true. |
| D3dc-4 | cosmetic | Wording. "\|∇_σC\| ≤ h′(\|σ\|), which is exact": the bound is attained only where \|Re σ\| ≤ 1 < \|σ\| and is strict beyond. New: "proved … on \|η\| < 1.8927", where (iii)(b) and the Status say ≤, and the limit on \|η\| ≤ 1.995 is not stated there. In the Ω′ proof, "(strictly within 2rb of η/2)" is condition (1)'s wording; it holds for (2) too, but the hypothesis needs only "in the open disc". The step near \|σ\| = 1 rests on V's strict margin (−1.0004·10⁻¹²) beating h(r₀), which goes unsaid [md][own]. | D61 (iii)(b): "sharp (attained where \|Re σ\| ≤ 1 < \|σ\|)". New: "on \|η\| ≤ 1.8927 … and, by V, on \|η\| ≤ 1.995". Proof: "(in the open disc; within 2rb of η/2 in both sufficient conditions)" and one clause: "V's margin is uniform, so h(r₀) below it for r₀ near 1 covers the circle". |
| D3dc-5 | cosmetic | Stale text. PAPERS.md Paper 18: "σ*/σ < 198.99", "modulo both computer-assisted steps", "Missing. The two computer-assisted steps in interval arithmetic", "about sixfold per halving" (now 2011 to 10658 for 0.02→0.01 to 0.01→0.005, i.e. 5.3). HANDOFF: line 243 "about sixfold per halving"; line 62 lists `exp_gh_peano.py` parts as "P, M, L, LC, S, R" (Q, D, V missing); the test time is "about 2 s" (now about 10 s). `cert_seg_parts`' docstring opens with the old bound h′ + 2(\|Re σ\| − 1)₊ [own]. | PAPERS.md Paper 18: D‴.61's state as in HANDOFF §305–311, with 199.00 and "V in ball arithmetic", and Missing: D in ball arithmetic (or nothing, after D3dc-1) and the annulus 1.995 < \|η\| < 2. HANDOFF: the parts list and the timing. Docstring: state h′(\|σ\|) alone. |

## Sources

- [Johansson, Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic (arXiv:1611.02831; ARITH 24, 2017)](https://arxiv.org/abs/1611.02831v1) [snippet]
- [Johansson, Numerical integration in arbitrary-precision ball arithmetic (arXiv:1802.07942)](https://arxiv.org/pdf/1802.07942) [snippet]
- [Machine-Assisted Proofs, ICM 2018 panel (arXiv:1809.08062)](https://arxiv.org/pdf/1809.08062) [snippet]
- [Gómez-Serrano, Computer-assisted proofs in PDE: a survey (arXiv:1810.00745)](https://ar5iv.arxiv.org/html/1810.00745) [snippet]
- [Proof Generation from Delta-Decisions (arXiv:1409.6414): the mean-value form in interval arithmetic](https://arxiv.org/pdf/1409.6414) [snippet]
- [Contact geometry of the restricted three-body problem on S² (arXiv:2411.10885): Arb with adaptive bisection](https://arxiv.org/pdf/2411.10885) [snippet]
- [Moore–Kearfott–Cloud, Introduction to Interval Analysis (SIAM 2009)](https://woodslane.com.au/products/9780898716696) [snippet, bibliographic]
- [Interval methods in mathematical programming (centred form, overestimation)](https://www.reliable-computing.org/archive/preprints/interval-methods-in-mathematical-programming.pdf) [snippet]
- [Springer article citing Caprani–Madsen 1980 and Krawczyk–Nickel 1982 on the quadratic convergence of mean-value forms](https://link.springer.com/article/10.1007/BF02684417) [snippet, via search summary]
- [O'Sullivan, Revisiting the saddle-point method of Perron (arXiv:1702.03777)](https://arxiv.org/pdf/1702.03777) [snippet]
- [FAU dissertation: Arnold's problem, saddle-point method with interval arithmetic](https://fau.isle.flvc.org/node/57710) [snippet, record only]
- [Arb's `arb_overlaps` (flint-sys source)](https://docs.rs/crate/flint-sys/latest/source/flint-3.5.0/src/arb/overlaps.c) and [Sage's real_arb documentation (comparison semantics)](https://doc.sagemath.org/html/en/reference/rings_numerical/sage/rings/real_arb.html) [snippet]
- python-flint 0.9.0's installed docstrings for `arb`, `arb.upper`, `arb.lower`, `arb.mid`, `arb.rad` and `arb.log` [own]

## Lesson

1. **Port all the steps at once.** When one computer-assisted step is ported to ball arithmetic, port the others in the same round. The second port here (D's branch and bound) was 60 lines and 18 seconds; leaving it as "modulo the branch and bound alone" cost a round of qualified claims.
2. **Test a rigorous twin against its float twin, not only its sign.** Each ball function should be tested as an upper bound of its float version within a few ulps (u ≥ float and u − float ≤ ε). A test that the rigorous value is merely negative let three unsound mutants of the end piece pass.
3. **Plant threshold tests at the threshold.** Plants far from the threshold (u = 1.0, |σ_b| = 1.0045) do not tell 0 from 10⁻¹² or 1 from 1.001.
4. **Recheck independently.** A second ball enclosure built differently (here the natural interval extension of g′ on a rectangle) confirms soundness cheaply: 10.7 M bounds in about 7.5 minutes, 0 failures, 0.1% extra splits.

My files are all under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref108:
- own/: rerun_band.py, indep_ball.py, sanity_mybound.py, dom_bb_ball.py, dom_sanity.py, identity_grad.py, flint_semantics.py, attack.py, numcheck_own.py, d61.diff;
- mut/: spec108.py, recompute.py, r104/, r108/;
- logs/: band0–2.log, indep0–2.log, leaves0–2.pkl, dom_bb_ball.log, mut104.log, mut108.log, rc_*.log;
- gp_out/: the regenerated md and the raw copy.
