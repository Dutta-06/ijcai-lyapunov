# Referee report D3de, round 110: Proposition D‴.61 as revised after round 108 (both steps in ball arithmetic)

## Scope

**Commit.** I reviewed 88ef79f on claude/vigilant-clarke-hf939b. The branch has since moved to 8139e12 ("D‴.61 part VC: the corner at η = 2 as a blow-up limit"), via cd34fed. That commit changes D61.md, exp_gh_peano.py, the tests, the md and the raw. It is outside this round and I did not review it. The HANDOFF lines I flag below are still present at HEAD (checked with `git show HEAD:HANDOFF.md`).

**Setup**
- W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref110.
- I exported with `git archive 88ef79f | tar -x -C W/repo`. The sha256 of exp_gh_peano.py, test_gh_peano.py, D61.md, gh_peano.md and the raw all match 88ef79f's.
- W/repo/venv is a symlink to the repository's venv (python-flint 0.9.0, numpy 2.5.4, scipy 1.18.1).
- GP_OUT pointed at copies of the raw (W/out, W/outDB).
- Every job ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with a timeout. Pytest used `--basetemp` in W, and my copy of tools/mutants.py (W/tools_ref/mutants_bt.py) adds a per-slot basetemp and a timeout.
- At most two heavy processes of mine ran at a time. I stopped one of my own background tasks (a too-slow tight run) with TaskStop. I touched no other process.
- I used only read-only git commands. `git status` is clean.

**Read in full**
- review/drafts/D61.md.
- exp_gh_peano.py:
  - D: `gap_closed`, `quad_roots`, `bound_sup`, `box_bound`, `box_bound_ball`, `dom_bb`, `part_dom`, `part_D`, `part_DB`, `d_lines`;
  - V: the `cert_*` functions, `v_lines`, `ratio_floor`;
  - `report`, KEYS, `main`.
- gh_peano.md, parts D and V.
- tests/test_gh_peano.py: the D, V, round-104 and round-108 tests.
- tools/mutants/round108_D61_own.py and tools/mutants.py.
- The Round 108 entry of REVIEW_LOG.md and REVIEW_round108_D3dc.md.
- PAPERS.md Paper 18, and HANDOFF.md's D‴.61 lines (§1 table, §7, §8, §10).
- `git diff d5d8438 88ef79f` of these files.

## Summary

**Design.** D's branch and bound now has a ball twin.
- **`box_bound_ball`** takes these parameters in Arb at 120 bits:
  - R1 exact;
  - p, the rigorous least of R²cos 2φ over the corners, with the top edge at the true π/2;
  - B = R2 cos f1 and wm = 1 − 1/R2 at their upper ends, with wm² < ½ checked.

  It bounds h(w) = log w + (p/2)(w² − ½) + M(B(1 − w)) − ½ − log R1 + log 4 on (0, wm]. The candidates are wm, the kink and both pieces' roots, as balls from the stable formula. The result is a float rounded up.
- **`dom_bb(…, ball=True)`** accepts a box only if the float bound is below −10⁻¹² and the ball bound is below 0.
- **Part DB** stores the result, and `report` merges it in as D's branch and bound.
- **The claim.** On this, the draft claims (iii)(b) on |η| ≤ 1.995 and (iv) for r < 0.9900 (σ*/σ < 199.00), with both computer-assisted steps rigorous and no floating-point proviso.

**Verdicts**
- **Reproducibility is exact.**
  - The md regenerates byte for byte (sha256 e2a9052b4ddd8228…14c3d).
  - All 36 tests pass (16.2 s).
  - Part DB reran in 14.6 s and is repr-identical to the raw: (4935, −3.281257330289833·10⁻⁷, 0, −3.2812573349693787·10⁻⁷) and (402664, −1.9432233600014115·10⁻¹¹, 0, −1.9432017877494493·10⁻¹¹).
  - All three V bands with the ball recheck equal the raw dict for dict (7.0 s, 29.8 s, 275.1 s).
- **`box_bound_ball` is sound as a proof.** The analytic box bound, the candidate set, every enclosure and comparison, the root formula, the final rounding up and the tiling are all correct.
- **An independent check agrees.** I re-verified all 402,664 leaves with a differently built rigorous bound: Taylor-form bisection in w, no candidate set, no root formula. There were 0 failures.
- **An attack from the definition** found nothing: 3776 η, minimum slack 2.83·10⁻⁶.
- **"No floating-point proviso" is true in substance.** What remains in floating point is threshold arithmetic (r*, c_need), with margins of at least 1.8·10⁻⁶ against errors near 10⁻¹⁵. I rechecked it in Arb.
- **No surviving mutant changes a rigorous result.** I recomputed DB under all 20 survivors of the ball code.

No major issue.

**Findings**
- **D3de-1 (minor). Tests.**
  - The quadratic piece of `box_bound_ball` is never exercised. Y05 (first piece dropped) and Y06 (wrong quadratic) are unsound and survive. On 110 of my 20,000 random boxes they put the bound up to 1.69·10⁻² below the dense supremum of h. None of those boxes is on the data, and all have positive bounds.
  - The guards on the reporting lines are untested. Y26 and Y27 (an unresolved run counted), Y30 (the float worst stored as `ball_worst`) and Y31 (the "both … in ball arithmetic" condition reduced to a key-presence check) all survive.
  - Round 108's E10 is "killed" only by a TypeError crash, not by a semantic test.
- **D3de-2 (minor). HANDOFF still describes D‴.61 as floating point and qualified.**
  - The §1 status table (line 34): "modulo two floating-point computer-assisted steps".
  - §8: "r < 0.8053 (modulo a grid-checked inequality)".
  - §10: "by the floating-point branch and bound", which contradicts the same paragraph's last sentence; "h′(|σ|), exact"; "test_gh_peano.py 35".
  - All of these persist at HEAD.
- **D3de-3 (cosmetic). HANDOFF §7's probe numbers contradict its own blow-up limit.**
  - It says θ_dom ≈ 1.7(2 − |η|), with 0.087, 0.0153, 0.0087 and "<0.0022".
  - The exact values are 0.08606, 0.01484, 0.00731 and 0.00144. θ_dom/ε → 1.4397, the root of (φ − 1)²(2φ + 1) = ¾.
  - The rest of the blow-up checks out.
- **D3de-4 (cosmetic). "No floating-point proviso" is not literal yet.**
  - r* = 0.8052, 0.9900 and 0.0966 come from float bisections, checked by a float test with margin 10⁻⁷.
  - `d_lines` falls back to the grid radius R_dom as "verified" when no branch-and-bound run is complete.
- **D3de-5 (cosmetic). PAPERS.md Paper 18.**
  - The Claim is unqualified.
  - "Missing" lists the ρ < 0 rules as outside the frame picture. D61 puts ρ ∈ (−1, 0) inside it; the rules outside are those with |ρ| ≥ 1.
- **D3de-6 (cosmetic). Odds and ends.**
  - D61 cites no ball-arithmetic or validated-numerics source.
  - Only the smaller root 1/q of each quadratic can be a maximum, so the q/qa candidates, the qa-straddle branch and the kink are redundant. Record Y04, Y09 and E03 as sound, or simplify.
  - REVIEW_LOG says "`part_dom` itself now runs it"; it is `part_D`.
  - `report` lets a stored DB override a freshly run D.

## What I checked and how

**(1) Reproducibility** [own][md][raw]
- **The md.** `--from-raw` into W/out matches the committed md byte for byte (cmp; sha256 e2a9052b4ddd8228ca269c0ad6f85d19f9db86f6635ef3010faabba298414c3d for both).
- **Tests.** `pytest -q -p no:cacheprovider --basetemp=W/pt tests/test_gh_peano.py`: 36 passed in 16.24 s.
- **Part DB.** `--only DB` against a copy of the raw ran in 14.6 s. `DB` is identical to the raw's, and the regenerated md is still byte-identical.
- **Part V.** W/own/rerun_v.py reran `cert_band(*CERT_BANDS[i], ball=True)`. All keys equal the raw for bands 0, 1 and 2 (7.0 s, 29.8 s, 275.1 s). So round 108's added checks (last edge exactly Rb, |η_c| < R2 in Arb) change nothing on the data.
- **The raw's two copies of D's branch and bound.** D's own `dom['bb']` in the raw predates the change (no `ball_worst`). The printed ball values come from DB through `report`'s merge.

**(2) `box_bound_ball` as a proof** [own]
- **The analytic box bound** (by hand).
  - With w = 1 − s, the gap is log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4. This follows exactly from g_η and Λ⁺.
  - Over a box, −log R ≤ −log R1.
  - Since w² − ½ < 0 (wm ≤ 0.4717), the p-term decreases in p, so the least p bounds it. cos 2φ decreases on [0, π/2], and R² cos 2φ is monotone in R at fixed sign, so the least p is at a corner.
  - M is even and increasing in |a|, and 0 ≤ R cos φ ≤ R2 cos f1 = B.
  - The domain (0, 1 − 1/R] lies in (0, wm].
  - R1 = 1 is valid because the gap is −∞ (void) for R ≤ 1.
- **The candidate set is complete.**
  - h → −∞ as w → 0, so the supremum is attained.
  - M is C¹ (M(1) = 1 and M′(1) = 2 from both sides), so h is C¹ and the maximum is at wm or at a root of h′. On each piece w·h′ = Q(w) = qa w² + qb w + 1, with Q(0) = 1 > 0.
  - Hence only the smaller root 1/q can be a local maximum. With qa > 0, Q is positive, then negative, then positive; with qa ≤ 0 the larger root is negative or absent.
  - So the kink (E03), q/qa (my Y04) and the qa-straddle branch (my Y09) are redundant, and sound to drop.
- **Ball semantics** (my probes on python-flint 0.9.0).
  - Overlapping-ball comparisons are all False.
  - `arb(float)` is exact; `.upper()`/`.lower()` return exact arbs.
  - `float(arb)` rounds to nearest: −(0.1 − 2⁻⁸⁰) comes back as −0.1, below the true value. So the `nextafter` up is needed, and present.
  - Comparisons with log 0 are False.
- **Every `if` in `box_bound_ball`** compares exact arbs or fails safe:
  - `wm*wm < half`, `B > 1`, `qb < 0`;
  - `disc < 0` then `continue`; a straddling disc returns None;
  - `qa > 0 or qa < 0`, else |q/qa| > wm proved;
  - candidates that may lie in (0, wm] are kept; a candidate not certainly > 0 returns None;
  - M is taken at the upper end of a, with `v > best` exact.
- **The stable root formula is correct with balls.** The roots are q/qa and 1/q with q = (√disc − qb)/2 > 0, and the linear limit 1/q = −1/qb.
  - qa near 0: qa = p for the second piece is an exact point.
  - The first piece's qa = p + 2B² crosses 0 near arg η = π/3, but a straddle needs |qa| < 10⁻³⁶.
  - The discriminants vanish on curves, and a straddle likewise needs |disc| < 10⁻³⁶.
  - None of these branches fires on the data: 0 None results and 0 ball rejections over the whole tree.
- **The tiling** (W/own/leaves.py, with the script's own `box_bound` and `box_bound_ball`).

  | Rmax | leaves | ball rejections | area exact in Fractions | leaves at R1 = 1 | top leaves (f2 = float π/2) | bottom leaves | leaves at Rmax | smallest sides (R, φ) | ball − float |
  |---|---|---|---|---|---|---|---|---|---|
  | 1.78 | 4935 | 0 | yes | 8 | 9 | 50 | 360 | 1.52·10⁻³, 1.53·10⁻³ | – |
  | 1.8927 | 402664 | 0 | yes | 6 | 16 | 231 | 13640 | 1.70·10⁻⁶, 1.50·10⁻⁶ | ∈ [−9.1·10⁻¹⁶, 1.7·10⁻¹⁵] |

  - At 1.8927, every leaf lies inside [1, Rmax] × [0, float π/2].
  - The float π/2 is 6.12·10⁻¹⁷ below the true π/2. The top boxes take p at the true π/2, so the sliver is covered.
  - The float 1.8927 is the decimal + 4.90·10⁻¹⁷, and the float 1.995 is the decimal + 1.07·10⁻¹⁶. DOM_BB[−1] == CERT_BANDS[0][0] exactly, so the decimal disc |η| ≤ 1.995 is covered.
- **Independent rigorous recheck** (W/own/dom_indep.py, Arb at 128 bits).
  - The parameters come from monotonicity alone: cos 2φ is least at f2 (the true π/2 on the top edge), then R1 or R2 by its sign.
  - The supremum over w is bounded by adaptive bisection, with the second-order Taylor form h(c) + |h′(c)|r + ½max(0, sup h″)r². h′ is continuous and M″ ∈ {0, 2}. On (0, 2⁻¹²] a term-wise monotone bound is used.
  - Results: 1.8927, all 402,664 leaves below 0, 0 failures, 10,514,450 interval evaluations, 79.5 s; 1.78, all 4935 leaves, 0 failures.
  - The 12 worst leaves by branch and bound in w (W/own/dom_tight.py). At the worst, [1.892669351577759, 1.892672756958008] × [0.3426829190562185, 0.34268591511244484], my enclosure of sup h is [−1.94320183·10⁻¹¹, −1.94319956·10⁻¹¹]. The script's −1.9432017877·10⁻¹¹ lies 4·10⁻¹⁹ above my lower end. All 12 agree the same way.
- **Attack from the definition** (W/own/attack_dom.py). The gap sup_s[g_η(s) + Λ⁺(ηs)] − g_η(½) comes straight from the definitions, not the closed form: a 4001-point s-grid plus golden-section refinement, with the worst 8 redone at 200 bits.
  - Samples: 3776 η (3000 random; 360 at corners and centres of the 40 worst leaves; 400 on |η| = 1.8927; 16 at |η| = 1 + 10⁻³ … 1 + 10⁻¹²).
  - min(leaf bound − gap) = 2.83·10⁻⁶, with 0 violations. The largest gap is −3.118·10⁻⁶ (|η| = 1.8927, arg η = 0.3425).
  - The box bound's first-order overestimate (~3·10⁻⁶ on the smallest boxes) explains why the bound (−1.94·10⁻¹¹) sits so near 0.
- **Where the supremum lies.** On all 402,664 leaves the box bound's supremum is at wm (87,717) or at the linear piece's own critical point (314,947), never on the quadratic piece. Pointwise, the quadratic piece's interior maximum occurs only at |η| ≥ 1.957, arg η ∈ [0.637, 0.779] (584 of 640,000 grid points). This is why Y05 and Y06 are equivalent on the data (D3de-1).

**(3) What is still floating point or by hand** [own]
- **r* and c_need.** Float bisections; the test checks c_need(r) < R − 10⁻⁷ in float. In Arb at 200 bits:
  - c_need(0.99) = 1.99497694172, margin 2.3058·10⁻⁵ to 1.995;
  - c_need(0.8052) = 1.89269819252, margin 1.8075·10⁻⁶ to 1.8927;
  - c_need(0.0966) = 0.999825651279, margin 1.7435·10⁻⁴ to 1;
  - (1 + 0.99)/(1 − 0.99) = 199 exactly, and 1.8052/0.1948 = 9.26694, so "9.26" is a safe floor.
  - c_need increases with r, because ln u/(u − 1) decreases in u.
  - The frames' r (R ≤ 0.6345, S ≤ 0.5279) are floats with margins of at least 0.17.
- **The axes' proofs, by hand, checked.**
  - Real axis: the derivative is (w′ − 1)²/w′, with w′ = η(1 − s). The gap at s = 1/η and its first two derivatives vanish at η = 2. The third derivative, 2/(η − 1)³ − 4/η³, is positive iff η < 2^{1/3}/(2^{1/3} − 1) ≈ 4.85.
  - Imaginary axis: the gap is concave on (1, 3.41).
  - The closed forms equal the definition's gap at R = 1.1 … 1.999 (−2.502·10⁻¹⁰ ≈ −¼(2 − R)³; −5.0025·10⁻⁷ ≈ −½(2 − R)²).
- **The other analytic pieces.** Λ⁺'s max_θ(2a cos θ − ½cos 2θ) = M(a) + ½ is right.
  - The disc |η| < 1 is void.
  - The symmetries: the gap depends on R, Re η² and |Re ηs|.
  - V's straight piece: log(1 − x²) + x² < 0.
  - The closed-form algebra I rederived; the leaf bounds also agree with the definition-based attack above.
- **Condition (1) is not used by the rigorous claim.** It remains a grid check, said as such.
- **Verdict on (3).** No computer-assisted step of (iii)(b) on |η| ≤ 1.995, or of (iv) for r < 0.9900, is float-only. Only the threshold arithmetic is float, with margins about 10⁹ times the float error. "No floating-point proviso" is true in substance (D3de-4).

**(4) The draft's statements** [md]
- **These match the raw and are earned:**
  - (iii)(b): "verified by a branch and bound rechecked in ball arithmetic", "below −1.0004·10⁻¹²", "with the branch and bound, the limit holds on the disc |η| ≤ 1.995";
  - (iv): "for r < 0.9900 … both computer-assisted steps rigorous in ball arithmetic";
  - the Status: 402664 boxes, the largest rigorous bound −1.9432·10⁻¹¹;
  - New: "on |η| ≤ 1.8927 … and, by V, on |η| ≤ 1.995";
  - the Numerics: −3.2813·10⁻⁷ and −1.9432·10⁻¹¹.
- **Round 108's wording items are applied:** "sharp", "within 2rb in both sufficient conditions", and the circle step's strict margin.
- **Outside the draft (D3de-2, D3de-5).** PAPERS.md Paper 18 is consistent on the rigour, but its Claim and "Missing" bullets need fixing. HANDOFF is stale.

**(5) Tests and mutants** [own] (W/mut/log_r108.txt, log_r110.txt, recompute.log)
- **Round 108's spec.** 13 of 17 are killed (Z18, Z19, Z21, Z22, Z30, E02, E04, E05, E07, E08, E09, E10, E11). E01, E03, E06 and E12 survive, as logged. E10's kill is `all(True)` raising TypeError (fr_rc 1), not a semantic test.
- **My spec** (W/mut/spec110.py, 32 mutants): 11 killed (Y03, Y07, Y14, Y16, Y17, Y20, Y23, Y24, Y25, Y28, Y29); 21 survive.
- **The survivors in the D ball code**, plus E01, E03 and E06, re-run through `dom_bb(1.78, ball=True)` and `dom_bb(1.8927, ball=True)`: all 20 give results identical to the control.

  On 20,000 random boxes against the control (the same boxes for every mutant):

  | class | mutants | effect |
  |---|---|---|
  | sound by mathematics | Y04 (q/qa dropped), Y09 (qa-straddle check skipped), Y32 (kink for B ≥ 1), E03 (kink dropped) | identical; only 1/q can be a maximum |
  | below double resolution, or a branch the data never take | Y01, Y02, Y08, Y10, Y11, Y12, Y13, Y18, Y19, Y21, Y22, E01, E06 | identical |
  | unsound, harmless on the data | Y15 (R1 = 1 lifted to 1 + 10⁻¹⁵) | 430 boxes, −1.8·10⁻¹⁵; the gap there is ≤ −30 |
  | **unsound, untested branch** | **Y05 (first piece dropped), Y06 (qa = p + B²)** | **110 boxes lowered by up to 1.69·10⁻² and 9.6·10⁻³, below the dense supremum of h. Example: (1.8080676209435946, 1.821991331786887) × (0.9495496783755496, π/2), control 0.349033 = dense sup, Y05 0.332124. No sign flips (all changed values ≥ 0.047).** |
  | reporting guards untested | Y26, Y27 (`open` ignored in R_bb, R_ver), Y30 (part_DB's `ball_worst` = the float worst; identical at the printed .4e), Y31 ("both … in ball arithmetic" on key presence only) | equivalent on the data |

  The test's 200 random boxes (widths ≤ 5·10⁻³, R1 never exactly 1) and its 1.78 pin do not reach the quadratic piece.

**(6) Literature** [snippet] unless marked
- **Interval and ball branch and bound with Taylor bounds as proof.**
  - Hales' Kepler algorithms (arXiv:math/0205209).
  - Recent Arb box computations: the real Grothendieck constant (arXiv:2609.20074; "each accepted box … evaluated with 96-bit ball arithmetic"); a dependent Gabor system (arXiv:2608.08190; "no pointwise floating-point sample enters the proof").

  D now meets that standard: the float pass only shapes the tree, and every accepted leaf is ball-checked. That is the precise sense of "no floating-point proviso".
- **Background.** Arb (Johansson, IEEE Trans. Comput. 2017; arXiv:1611.02831) [memory; round 108's snippet]. Arb's `arb_overlaps` is a conservative test (flint source) [snippet]. Validated numerics: Tucker 2011 [memory].
- **The last annulus.**
  - HANDOFF §7's blow-up checks out [own]:
    - the formula for Ψ_φ, from h(1 + δ) = δ² − δ³/3 and the identity that the squares of |σ| − 1 and Re σ − 1 differ by ε³(x − ½)(y + φ)²;
    - the valley floor y = min(2φx, φ), with −(1 + 4φ²)x² and −x² − φ² − (x − ½)³/3 on it;
    - φ_dom = 1.4396926 from (φ − 1)²(2φ + 1) = ¾. Measured θ_dom/ε = 1.721, 1.484, 1.461, 1.444 and 1.440 at |η| = 1.95, 1.99, 1.995, 1.999 and 1.9999;
    - the gap at θ = 0 is −ε³/4;
    - on the floor, E/(ε³Ψ) ∈ [0.951, 1.011] at ε = 10⁻³ (and [0.51, 1.09] at 10⁻², where the inner scale is already strained).
  - The degeneracy is a saddle (η/2) meeting the curve |σ| = 1, where the bound changes form. That is the analogue of a stationary point approaching an endpoint or singularity: Chester–Friedman–Ursell 1957 and Bleistein 1966; uniform asymptotics as a stationary point approaches an endpoint (arXiv:1707.07927); for Hermite polynomials, Airy-type coalescence at the turning point (Gil–Segura–Temme, arXiv:2001.02097).
  - Only the valley inequality is needed, not uniform asymptotics, because the claim is on compact subsets of |η| < 2. A rigorous finish could certify the rescaled bound Ψ_ε(w, φ) = E/ε³ in ball arithmetic with ε ∈ [0, 0.005] as an interval variable, ε = 0 included.
  - That needs a cancellation-free form (a Taylor model with remainder): computer-assisted blow-up proofs replace f(x)/x at 0 by f′(0) (De Gregorio model, arXiv:1905.06387; Hou–Luo, arXiv:2106.05422) [snippets]. This is the desingularisation idea of geometric blow-up (Krupa–Szmolyan 2001) [memory].
  - Plus the fixed-θ regime at scale ε², the outer region at η = 2 (E ≤ −x³/3 on the real segment), and the overlaps, as HANDOFF lists.
  - I found no validated certificate of a degenerate-saddle Laplace or valley step in the literature.

## Verdicts per part

- **(i), (ii), (iii)(a), (iii)(c), (v):** unchanged and proved as before (P_j's polynomial form is still a sketch).
- **(iii)(b)**
  - **|η| < 1 and both axes:** proved; I rechecked the axes' argument.
  - **1 < |η| ≤ 1.8927:** proved, computer-assisted and rigorous in ball arithmetic. The tiling is exact, my independent bound verifies every leaf, and the attack from the definition finds nothing.
  - **1.8927 ≤ |η| ≤ 1.995:** proved by V, rerun identical.
  - **Union:** the closed disc |η| ≤ 1.995, decimal included.
  - **1.995 < |η| < 2:** conjecture P.
- **(iv):** proved for r < 0.8052 (σ*/σ < 9.26) and r < 0.9900 (σ*/σ < 199.00). The thresholds are right in Arb. "No floating-point proviso" is accurate except for the float threshold arithmetic (D3de-4).

## Status of the round-108 issues

| id | status |
|---|---|
| D3dc-1 (D in ball arithmetic) | **Resolved.** `box_bound_ball`, `dom_bb(ball=True)` and part DB are shipped; 402,664 boxes rigorous, 0 rejections; rerun identical; my independent recheck agrees. |
| D3dc-2 (end piece and thresholds untested) | **Resolved.** Z18, Z19, Z21, Z22 and Z30 are now killed. New gap: the quadratic piece and the reporting guards (D3de-1). |
| D3dc-3 (residues) | **Largely resolved.** Exact contiguity (E11 killed), the last edge = Rb, and \|η_c\| < R2 in Arb are in. c_need is checked only in float (D3de-4). |
| D3dc-4 (wording) | **Resolved in D61.** HANDOFF §10 still says "exact" (D3de-2). |
| D3dc-5 (stale text) | **Partly resolved.** PAPERS.md now says 199, "rigorous" and the new Missing; HANDOFF no longer says "sixfold" and has the parts list and the timing; the docstring is fixed. Still stale: HANDOFF lines 34, §8 and §10 (D3de-2); PAPERS.md's Claim and Missing (D3de-5). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3de-1 | minor | Untested branches. The unsound Y05 (quadratic piece dropped) and Y06 (wrong quadratic) survive: on 110 of 20,000 random boxes they lower `box_bound_ball` by up to 1.69·10⁻² below the dense sup of h (e.g. (1.8080676209435946, 1.821991331786887) × (0.9495496783755496, π/2): 0.349033 → 0.332124). None of these boxes is on the data, and DB recomputes identically. The reporting guards are untested: Y26 and Y27 (`open` ignored), Y30 (float worst stored as `ball_worst`), Y31 (the "both rigorous" condition without teeth) all survive. E10 is killed only by a TypeError [own: W/mut/res_r110.jsonl, recompute.log, rc/*_boxes.npy]. | tests/test_gh_peano.py: <br>• `box_bound_ball` ≥ the dense max of its own h (with the box parameters) on wide boxes, including the one above and boxes at \|η\| ≈ 1.96–1.99, arg η ≈ 0.64–0.78; <br>• plant a bb with open > 0, and one without `ball_worst`, into `d_lines`/`v_lines` and assert the radius and wording fall back; <br>• `part_DB` with DOM_BB = (1.78,) equals `dom_bb(1.78, ball=True)`'s fourth entry. <br>Record Y04, Y09 and E03 as sound (only 1/q can be a maximum: w·h′ = Q, Q(0) = 1 > 0). |
| D3de-2 | minor | HANDOFF is stale on D‴.61 (all still at HEAD): <br>• line 34: "modulo two floating-point computer-assisted steps" and "partly proved (modulo computer-assisted steps)"; <br>• §8: "r < 0.8053 (modulo a grid-checked inequality)" and "the deformed path may close the rest"; <br>• §10: "by the floating-point branch and bound" (contradicted two lines later), "h′(\|σ\|), exact", "test_gh_peano.py 35" (36 now). [md] | HANDOFF: the table row → "proved on \|η\| ≤ 1.995 (both computer-assisted steps rigorous in ball arithmetic)"; §8 → r < 0.9900, rigorous; §10 → "by the branch and bound (ball arithmetic)", "sharp", 36. |
| D3de-3 | cosmetic | HANDOFF §7's probe says θ_dom ≈ 1.7(2 − \|η\|) (0.087, 0.0153, 0.0087, <0.0022). By bisection on `gap_closed` the values are 0.08606, 0.01484, 0.00731, 0.00144, and θ_dom/ε → 1.4397 = φ_dom, the limit §7 itself derives. The rest of the blow-up (Ψ_φ, the floor, −ε³/4, floor ratio 0.951–1.011 at ε = 10⁻³) checks out [own]. | HANDOFF §7: the exact θ_dom and "θ_dom/ε → φ_dom = 1.4397". When written up (part VC), certify Ψ_ε with ε ∈ [0, ε₀] as an interval through a cancellation-free Taylor form. |
| D3de-4 | cosmetic | "No floating-point proviso" is not literal. r* = 0.8052, 0.9900 and 0.0966 are float bisections checked in float with margin 10⁻⁷. In Arb: c_need(0.99) = 1.99497694 (margin 2.31·10⁻⁵), c_need(0.8052) = 1.89269819 (1.81·10⁻⁶), c_need(0.0966) = 0.99982565 (1.74·10⁻⁴). `d_lines` falls back to R_dom (a grid radius) as "verified" if no run has open = 0 [own]. | test (or md): c_need at the printed floors in Arb, strictly below their radii. `d_lines`: no fallback to R_dom (print "unverified"). |
| D3de-5 | cosmetic | PAPERS.md Paper 18. The Claim states t(v) and the sharp constants for "every scaled Gauss–Hermite rule", but they are proved only for frames with r < 0.99 and conjectural beyond (D‴.58 (ii)). "Missing" lists "the rules below the Fock scaling (ρ < 0) … in the frame picture", yet D61 (v) and (R) include ρ ∈ (−1, 0) (5 of the 22); the ones outside have \|ρ\| ≥ 1 [md]. | Paper 18: qualify the Claim ("proved for every rule with a Mehler frame and t(v) < 0.99"); Missing: "the rules with \|ρ\| ≥ 1 (moment matching among them)". |
| D3de-6 | cosmetic | D61's sources cite nothing for ball arithmetic or validated numerics, though the main claims rest on it. REVIEW_LOG round 108: "`part_dom` itself now runs it", but it is `part_D`. `report` replaces even a fresh D's bb by a stored DB (stale if DOM_BB changes and only D is rerun) [md]. | D61 Classical/sources: Johansson's Arb (2017), and a validated-numerics reference or Hales' interval branch and bound as precedent. REVIEW_LOG: `part_D`. `report`: prefer D's bb when it carries `ball_worst`. |

## Sources

- [Hales, Some algorithms arising in the proof of the Kepler conjecture (arXiv:math/0205209)](https://arxiv.org/pdf/math/0205209) [snippet]
- [A computer-assisted upper bound of 1.7813 for the real Grothendieck constant (arXiv:2609.20074)](https://arxiv.org/pdf/2609.20074) [snippet]
- [The Minimum Cardinality of a Dependent Finite Gabor System Is Four (arXiv:2608.08190)](https://arxiv.org/pdf/2608.08190) [snippet]
- [On the Finite Time Blowup of the De Gregorio Model for the 3D Euler Equation (arXiv:1905.06387)](https://arxiv.org/pdf/1905.06387) [snippet]
- [Asymptotically self-similar blowup of the Hou–Luo model (arXiv:2106.05422)](https://arxiv.org/pdf/2106.05422) [snippet]
- [Numerical validation of blow-up solutions of ODEs (arXiv:1606.03039)](https://arxiv.org/pdf/1606.03039) [snippet]
- [Method of Chester–Friedman–Ursell (Proc. Camb. Phil. Soc. 53 (1957) 599–611)](https://en.wikipedia.org/wiki/Method_of_Chester%E2%80%93Friedman%E2%80%93Ursell) [snippet]
- Bleistein, Uniform asymptotic expansions of integrals with stationary point near algebraic singularity, CPAM 19 (1966) 353–370, via [arXiv:2204.09305](https://arxiv.org/pdf/2204.09305) and [DLMF §2.3](https://dlmf.nist.gov/2.3) [snippet]
- [Uniform asymptotics as a stationary point approaches an endpoint (arXiv:1707.07927)](https://arxiv.org/pdf/1707.07927) [snippet]
- [Gil–Segura–Temme, Numerical evaluation of Airy-type integrals arising in uniform asymptotic analysis (arXiv:2001.02097)](https://arxiv.org/pdf/2001.02097) [snippet]
- [Huybrechs–Kuijlaars–Lejon, Gaussian quadrature with coalescing saddle points (arXiv:1806.06549)](https://arxiv.org/abs/1806.06549v2) [snippet]
- [FLINT `arb_overlaps` source (flint-sys)](https://docs.rs/crate/flint-sys/latest/source/flint-3.5.0/src/arb/overlaps.c) [snippet]
- Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017), arXiv:1611.02831 [memory; round 108's snippet]
- Tucker, Validated Numerics (Princeton, 2011); Krupa–Szmolyan, Extending geometric singular perturbation theory to nonhyperbolic points (SIAM J. Math. Anal. 2001) [memory]
- python-flint 0.9.0, my own probes of comparisons, `float()`, `upper`/`lower` and conversions [own]

## Lesson

1. **Recheck a rigorous port with a different algorithm, not the same one in balls.** A Taylor-form bisection without the candidate set verified all 402,664 leaves in 80 s. That is independent of the analytic argument the port shares with its float twin.
2. **A mutant that survives on the data can still be wrong as a function.** Test a rigorous bound against a dense evaluation of what it bounds, on boxes where each branch is the maximiser. Here the quadratic piece only matters for |η| ≥ 1.957 or wide boxes, so the shipped tests never reach it.
3. **A crash is not a kill.** Check that a mutant's failing test fails on the semantics (E10 failed on a TypeError).
4. **When the main claim changes, update every status line in the same commit.** That includes the project's status table and the handoff's summaries: line 34 still says "floating point" one round after both steps became rigorous.

My files are all under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref110:
- own/: leaves.py, dom_indep.py, dom_tight.py, attack_dom.py, rerun_v.py;
- logs/: leaves.log, leaves_1.78.npy, leaves_1.8927.npy, indep_full.log, indep_*.txt/.npy, dom_tight.log, attack_dom.log, rerun_v.log;
- mut/: spec110.py, recompute.py, log_r108.txt, log_r110.txt, res_r108.jsonl, res_r110.jsonl, recompute.log, rc/;
- tools_ref/mutants_bt.py;
- out/ and outDB/: the regenerated md and the raw copies;
- pytest_log.txt.
