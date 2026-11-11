# Referee report D3cy, round 104: Proposition D‴.61 with part V (the valley certificate)

## Scope

I reviewed commit b98e652 on branch claude/vigilant-clarke-hf939b.

- **Setup**
  - **Work directory.** I exported the commit with `git archive b98e652 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref104. GP_OUT = W/gp_out holds a copy of the raw. The repository's venv is symlinked as W/repo/venv.
  - **Ball arithmetic and high precision.** These ran in a separate venv, W/ownenv, with mpmath 1.4.1 and python-flint 0.9.0 (Arb). Nothing was installed into the repository or its venv.
  - **Run conditions.** Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, each job under `timeout`. Pytest used `--basetemp` under W/pyt; for tools/mutants.py this was passed through PYTEST_ADDOPTS.
- **Processes.** At most two heavy processes of mine ran at a time. Other referees' jobs (ref105, ref106) shared the machine and I left them alone. The only process I stopped was my own exploration `own/beyond.py` (PID 19590), after checking /proc/19590/cmdline. It had passed its time limit's useful part.
- **Repository untouched.** I used only read-only git commands (archive, log, show, diff, rev-parse).
- **Read in full**
  - D61.md.
  - experiments/exp_gh_peano.py part V and the CERT_* constants. Also what part V relies on: `path_F`, `lam_plus`, `r_star`, `c_need`, `dom_bb`, `box_bound`, `bound_sup`, `v_lines`.
  - gh_peano.md, parts D and V.
  - tests/test_gh_peano.py.
  - tools/mutants/partV_D61_own.py and tools/mutants.py.
  - The Round 99 entry of REVIEW_LOG.md and REVIEW_round99_D3ct.md.
  - The diff 49d78c8..b98e652 of these files: round 99 applied in 38fcf58, part V in b98e652.
  - HANDOFF.md's D‴.61 lines.
- **Web.** Search snippets only. arxiv.org does not resolve from this session (ENOTFOUND), so no full text was opened. To open full texts, arxiv.org would have to be allowed in the session's network settings.

## Summary

**Design.** In the saddle's coordinates σ = (η/2)(1 + z), part V uses E(z; η) = F_η(σ) − F_η(η/2) = log|1 − z²| + Re(cz²) + C(σ), with c = η²/4.
- **The analytic piece.** The straight piece from z = −1 to b, at z = t = 0.9(2 − R2)/R2, is shown negative by hand, except at the saddle.
- **The verified piece.** From b, one polyline per η-box: the box centre's normalised descent flow, decimated. On each segment the bound E(m; η_c) + ℓL_z + δL_η is computed over the segment × the box's covering disc, and the segment is bisected until the bound is below −10⁻¹².
- **The end piece.** From |1 − z| ≤ 0.02 to z = 1, a closed-form bound.
- **The result.** With 2454 + 7338 + 41008 boxes, by the symmetries and with D's branch and bound, the limit (iii)(b) holds on |η| ≤ 1.995 and (iv) holds for r < 0.9900.

**Verdicts**
- **Reproducibility.**
  - `--from-raw` reproduces the md byte for byte.
  - All 34 tests pass.
  - All three bands rerun bit-identical to the raw (9.4 s, 54 s, 434 s).
  - Every number new in D61 is in the md.
- **(1) The identity is exact.**
- **(2) The analytic part and the deformed-path proof are correct.**
  - The proof needs less than Ω′'s definition states, and V's path does not always meet that definition (D3cy-1).
- **(3) `cert_seg_ub` is sound, term by term.**
  - Its gradient bound for C carries a redundant term (D3cy-4).
  - I re-verified all 77,355,579 segment bounds of the three bands in ball arithmetic: no failure. The rigorous bound exceeds the float bound by at most 3.4·10⁻¹⁴, against a slack of 10⁻¹². So V is a computer-assisted proof in the strict sense once such a recheck is part of the repository (D3cy-2).
- **(4) `cert_end` is sound.**
- **(5) Coverage is exact**, up to a sliver of width 1.2·10⁻¹⁶ rad at the imaginary axis, which the slack absorbs.
- **(7) No point found where the bound is ≥ 0.**
  - 2150 paths were sampled; the largest F − F(η/2) is −2.54·10⁻⁸.
  - 58858 random segment bounds were checked against dense sampling: 0 violations.
- **(8) Constants.** r*(1.995) = 0.990046 and the floor 0.9900 are right. The printed 198.99 is a floating-point artifact of 199 (D3cy-5).
- **(9) Mutants.**
  - Of my 33 part-V mutants, the tests kill 17.
  - 7 change the recomputed bands, and 6 of these are unsound. None is caught, because no test recomputes V against the raw (D3cy-3).
  - Of the round's own 21, W19 survives, contrary to "all killed"; it is equivalent on the data.
- **(10) Literature and wording.**
  - In validated numerics, floating point without rounding control is not a proof.
  - The draft flags "floating point, not interval arithmetic" every time, which is honest. But it uses "certified" and "certificate", which name the rigorous version.

No major issue.

**Findings**
- **D3cy-1 (minor). Ω′'s definition and what V certifies.** For 75 band-0 boxes, b lies in B(η/2, rb) for part of the box (d/rb ≥ 0.851). There "the condition holds as soon as…" is false as written. The proof never uses the ball.
- **D3cy-2 (minor). "Certified" for a floating-point computation.**
  - A ball-arithmetic recheck of the whole of V (about 29 min of CPU) finds no failure.
  - The float corner distances fall one ulp short of the rigorous ones in 24–27% of boxes.
  - The float π/2 leaves a sliver at the imaginary axis outside every box.
  - All of this is absorbed by the slack.
- **D3cy-3 (minor). Tests and mutants.**
  - 16 of my 33 mutants survive.
  - The unsound survivors Y02, Y03, Y07, Y16, Y21 and Y22 change the recomputed bands, and nothing catches them.
  - W19 of the round's own spec survives.
- **D3cy-4 (cosmetic, constructive). A redundant gradient term.**
  - |∇C| ≤ h′(|σ|) holds exactly. Dropping the (·)₊ terms (sound) cuts the bands to 760, 2011 and 10658 boxes.
  - The cost is concentrated in a sliver at the real axis: 98.3% of band 2's boxes lie within 5° of it.
  - Beyond 1.995 the method fails even with the tighter bound.
- **D3cy-5 (cosmetic). 198.99.** It is (1.99)/(0.01) = 199 in floating point.
- **D3cy-6 (cosmetic). Stale text.**
  - The proof of (iv) still says "by it for r < r*, conditionally on P beyond" and "inside |η| < 1.8927 when r < r*".
  - HANDOFF.md still states the round-99 position.
  - "Sixfold" rests on one ratio.

## What I checked and how

**Reproducibility** [own][md][raw]
- **The md.** `--from-raw` gives an identical md (sha256 a48ceb57…).
- **Tests.** `pytest tests/test_gh_peano.py`: 34 passed in 6.1 s.
- **The bands.** W/own/rerun_band.py reran each band with the script's own functions. The bisection was re-implemented only to record the segments. All three result dicts are `repr`-identical to the raw: 2454/892284/−2.266867430732503e-12, 7338/6276814/−1.014610251317985e-12 and 41008/70186481/−1.0002938226279955e-12. They took 9.4 s, 54.1 s and 434 s; the last includes recording.
- **The raw outside V.** Parts P…D are unchanged since 38fcf58.
- **Numbers.** Replicating tools/numcheck.py on the diff 38fcf58..b98e652 of D61.md: every added number is in the md.

**(1) The identity** [own] (W/own/identity.py)
- **By hand.** η − σ = (η/2)(1 − z).
  - The logs: log|η − σ| + log|σ| = 2log(|η|/2) + log|1 − z²|.
  - The quadratic terms: (η − σ)²/2 + σ²/2 = (η²/8)((1 − z)² + (1 + z)²) = (η²/4)(1 + z²).
  - At the saddle: F(η/2) = 2log(|η|/2) + Re η²/4.
  - Outside the disc, with a = Re σ: M(a) − ½ − log ρ − (a² − b²)/2 equals h(ρ) for |a| ≤ 1 and h(ρ) − (|a| − 1)² for |a| > 1.
- **At 50 digits.** Over 20000 random (η, z), with |η| ∈ [1, 2), all angles and z ∈ [−2, 2.5] × [−1.5, 1.5]:
  - my F from its definition against the z-form: max difference 3.2·10⁻⁵⁰;
  - Λ⁺'s two forms: 2.1·10⁻⁵⁰;
  - `cert_E` against the 50-digit value: 2.6·10⁻¹⁵ relative.

**(2) The proof route, re-derived** [own]
- **The rate.** At u = √Kσ the u-integrand's rate per K is F, with C replaced by (1/K)log|c_K(√Kσ)|.
- **Two valid upper rates.**
  - (a) gives rate 0 uniformly on |σ| ≤ r₀, for every r₀ < 1.
  - The contour |ζ| = ρ₀ gives Λ⁺(σ) + ε uniformly on compacts of σ ≠ 0, everywhere, inside the disc too, where Λ⁺ = h(|σ|) ≥ 0. This holds for ρ₀ near 1 and K large, and needs I_K → (1 − ζ²)^{−1/2} on |ζ| = ρ₀.
  - So C̃ = 0 on |σ| ≤ r₀ and Λ⁺ beyond is a valid rate, and C̃ − C ≤ h(r₀) on r₀ < |σ| < 1.
- **The margin carries over.** V's bounds are < −10⁻¹² on every segment of a box, and the end piece is ≤ −1.34.
  - Take r₀ between max|σ_b| = 0.99975 and 1 with h(r₀) < 5·10⁻¹³.
  - Then the true rate on b → η is uniformly below F(η/2).
  - The unit circle is therefore harmless, but only because V's margin is strict and uniform per box.
- **The other pieces.**
  - The piece 0 → a, at z ∈ [−1, −t], is ≤ −ε′ away from the saddle.
  - The chord [a, b] lies in |σ| ≤ 0.99975, a compact subset of the open disc, so (a) holds uniformly there.
  - Laplace's method along u = η/|η| applies, since Re(φ″u²) = (2Re η² − 8)/|η|² < 0, and gives the same main term as the steepest-descent direction.
  - At the endpoints, s^{K−1} and (1 − s)^K send F → −∞, and the bounds are uniform near η because |η| > 1.
  - Cauchy's theorem applies: the integrand is entire.
- **Uniformity.** The closed disc |η| ≤ 1.995 is covered by D's region and V's closed boxes, each with a strict margin, and the bounds are continuous in η. So the limit is uniform on |η| ≤ c for c ≤ 1.995, which is what (iv) needs.
- **Hypotheses used:**
  - (a) on a compact neighbourhood of the chord;
  - the uniform contour rate Λ⁺;
  - entireness;
  - strict F < F(η/2) on the compact path except at η/2;
  - a chord direction with Re(φ″u²) < 0, inside the open disc;
  - one path per η-box.
- **Not used:** the ball B(η/2, rb), d > rb, and distinctness of the components (which follows by Cauchy, as round 99 noted).
- **The straight piece, numerically.** At 50 digits, over 3000 η up to |η| = 1.9999, max E/x² = −0.0017. In the attack below, at 30 digits, max (F − F₀)/x² = −0.0050.
- **D3cy-1's evidence** (from the recorded boxes).
  - In band 0, 75 of 2454 boxes have t < (2 − R1)/(2R1), so b ∈ B̄(η/2, rb) near |η| = R1. These are the outer first-grid boxes, e.g. [1.958475, 1.98]; d/rb ≥ 0.8512.
  - In bands 1 and 2 no box has this.
  - The polylines come closest to the saddle at b in bands 0 and 1. In band 2, 10 come closer than b, but none closer than 1.1985·rb.

**(3) `cert_seg_ub`, term by term** [own]
- **The split.** E(z; η) − E(m; η_c) is the z-part at η_c (≤ ℓ·sup|∇_zE(·; η_c)|) plus the η-part at fixed z (≤ δ·sup|∇_ηE| along [η_c, η]).
  - [η_c, η] lies in disc(η_c, δ) ∩ {|η| ≤ R2}, which is convex and contains η_c (|η_c| = (R1 + R2)/2) and the box. So using R2 in L_η is valid.
  - E is C¹ off z = ±1, since h(1) = h′(1) = 0 and (|a| − 1)₊² is C¹. So the mean-value inequality applies.
- **g = log(1 − z²) + cz².**
  - g′ = 2z(c − 1/(1 − z²)) and g″ = 2c − 2(1 + z²)/(1 − z²)².
  - |1 − z²| ≥ |1 − m²| − ℓ(2|m| + ℓ) and |z| ≤ |m| + ℓ.
  - |g′(z)| ≤ |g′(m)| + ℓ·sup|g″|. All correct.
- **C's gradient.**
  - |∇_zC| = (|η_c|/2)|∇_σC|, by conformality.
  - h′(ρ) = ρ − 1/ρ is increasing; rmax and amax bound |σ| and |Re σ| on the segment.
  - L_η: |η||z|²/2 from Re(cz²); (|1 + z|/2)|∇_σC| at |σ| ≤ (R2/2)(|1 + m| + ℓ).
  - All correct.
- **A redundant term.** For |a| = |Re σ| > 1, |∇C|² = h′² − 4(|a| − 1)(1 − |a|/ρ²) ≤ h′². So |∇C| ≤ h′(|σ|) exactly, and the 2(·)₊ terms are redundant. On 10⁶ points the maximum of |∇C|/h′ was 1.0000000000000004 (D3cy-4).
- **The covering disc.** |Re^{if} − η_c|² is convex in R and largest at the extreme angle, so the farthest box point is a corner. In fact it is always an R2 corner: the difference is (R2 − R1)·2R̄(1 − cos Δ) ≥ 0.
- **Random segments** (W/own/seg_attack.py): 58858 finite bounds on random boxes of the annulus, with segments in five kinds of region: crossing |σ| = 1, at Re σ = 1, near z = 1, near the saddle, and the paths' region.
  - Each bound was compared with my F from its definition at 33 points × 121 values of η (a 9 × 9 grid plus 40 random interior points).
  - 0 violations; the least slack is 6.4·10⁻⁹, near the saddle, and the median 3.8·10⁻³.

**(3′) Floating point: a rigorous re-verification** [own] (W/own/arb_verify.py and W/own/arb_stream.py, Arb at 120 bits)
- **Method.** The same formula as `cert_seg_parts`, evaluated on the exact binary inputs.
  - Each monotone term is taken at its rigorous upper end, with h′ and (·)₊ evaluated there.
  - E is enclosed on the exact midpoint's ball.
  - δ is recomputed rigorously from the corners, with each top edge equal to the float π/2 extended to the true π/2.
  - It also checks that the verified segments tile the polyline exactly in floats (children share the float midpoint, so the verified path is a continuous chain from b), that b is inside the disc at R2, and that |1 − z_n| ≤ 0.02.

| band | segments | failures | largest rigorous bound | rigorous − float, at most | boxes whose float δ falls short | end piece | time |
|---|---|---|---|---|---|---|---|
| 0 | 892284 | 0 | −2.2669·10⁻¹² | 2.85·10⁻¹⁴ | 663 of 2454 | −1.382677 | 18 s |
| 1 | 6276814 | 0 | −1.0146·10⁻¹² | 3.41·10⁻¹⁴ | 1787 of 7338 | −1.357139 | 115 s |
| 2 | 70186481 | 0 | −1.0002·10⁻¹² | 3.30·10⁻¹⁴ | 10123 of 41008 | −1.344304 | ≈ 27 min |

- The shortfall in the corner distances is about one ulp. The time for band 2 is beyond the 434 s float rerun streamed alongside it.
- So the slack of 10⁻¹² exceeds the actual rounding by about 30× everywhere, and V's claim holds rigorously, modulo my 120-line arb code.

**(4) `cert_end`.**
- **The terms.** On the piece from z_n to 1:
  - |1 − z| ≤ ε₀ and |1 + z| ≤ 2 + ε₀;
  - Re(cz²) ≤ (R2²/4)(1 + ε₀)²;
  - C ≤ h((R2/2)(2 + ε₀)), since h ≥ 0 and is increasing on [1, ∞), and −(·)₊² ≤ 0.
- **Use per band.** The bound increases in R2, so the band's value at Rb covers every box, and `v_lines` requires it to be negative.
- **Values.** Analytically sound; rigorously −1.382677, −1.357139 and −1.344304.

**(5) Coverage** [own] (W/own/attack.py)
- **Exact tiling.** Each band's leaf areas equal (Rb − Ra)·π/2 exactly in floating point. Of 200000 random points, none is in no box and none is interior to two.
- **Edges.** Ra + (Rb − Ra)·nR/nR == Rb exactly (nR = 4, 2). Band edges coincide with each other, and band 0 starts at DOM_BB's 1.8927.
- **Decimal ends.** The float 1.995 is 1.995 + 1.1·10⁻¹⁶, so the decimal 1.995 is covered.
- **The π/2 sliver.** `math.pi/2` is 6.1·10⁻¹⁷ below π/2. The quadrant and its images therefore leave the sectors (π/2_fp, π − π/2_fp) and their mirror images, of width 1.2·10⁻¹⁶ rad around the imaginary axis, in no box. The axis itself is proved by hand; the rest is covered only by continuity and the slack. My arb recheck extends the top boxes to π/2.
- **Symmetries.** E(z̄; η̄) = E(z; η), since c̄ and C(σ̄) = C(σ). E(z; −η) = E(z; η), since σ ↦ −σ and C is even. ĝ_K is even and real on ℝ, so the limit transfers.
- **b inside the disc.** For every η of every box, |σ_b| ≤ R2(1 + t)/2 = 0.9 + 0.05R2 ≤ 0.999, 0.9995 and 0.99975 in the three bands.

**(7) Independent attack** [own]
- **Paths.** 1832 η: 1500 random, 32 at the band edges × {0, 10⁻⁹, 0.08605, 0.345, 0.61756, π/4, π/2 − 10⁻⁹, π/2}, and 300 box corners. For each, 2150 (η, box) pairs (corners in up to two boxes).
  - My F from its definition was evaluated at 30 digits, 12 points per polyline segment and 40 on the end piece.
  - The largest F − F(η/2) on the certified pieces is −2.54·10⁻⁸, at b for η = 1.995. Per band: −1.6·10⁻⁶, −2.0·10⁻⁷, −2.5·10⁻⁸.
- **Where the cost lies.**
  - The share of boxes within 5° of the real axis grows from band to band: 80.4%, 93.9%, 98.3%.
  - Band 2's smallest box sits at arg η ≈ 0.0058, |η| ≈ 1.99443. Its path's maximum is at b (−3.5·10⁻⁸): there 1 − Re c ≈ 2 − |η|, and b's margin is about −t²(1 − Re c).
  - The domination branch and bound does not close there either: 400000 boxes on [1.99, 1.995] × [0, 0.005] without finishing.
- **Beyond 1.995** (exploration, not a claim; with the tighter, sound gradient bound of D3cy-4). On 1.995 ≤ |η| ≤ 1.9975, 59275 boxes fail at depth 20 (side 9.8·10⁻⁶) in 1484 s. The method stops near 1.995, as the Status says.

**(8) Constants** [own]
- **r\*.** By my own 50-digit bisection, r*(1.995) = 0.990045694778; floored 0.9900; c_need(0.99) = 1.994977 < 1.995.
- **The ratio.** (1 + r)/(1 − r) is 199.0000 exactly at r = 0.99, and 199.918 at r*. In floating point (1 + 0.99)/(1 − 0.99) = 198.99999999999983, which floors to 198.99.
- **Earlier floors.** r*(1.8927) = 0.8052030 and r*(1.0) = 0.0966496: unchanged and right.

**(9) Tests and mutants** [own] (W/mut/spec104.py; tools/mutants.py, one slot; survivors recomputed on band 0 and on the slice 1.99 ≤ |η| ≤ 1.9902 by W/mut/recompute.py, against a recomputed control)
- **Killed (17):** Y01, Y06, Y11, Y12, Y13, Y14, Y17, Y19, Y20, Y23, Y25, Y26, Y27, Y28, Y30, Y32, Y33.
- **Survivors that change the recomputed bands (7):**

  | mutant | change | recomputed band 0 |
  |---|---|---|
  | Y02 (unsound) | rB with \|η_c\| for R2 | 2305 boxes |
  | Y03 (unsound) | rmax without + ℓ\|η\|/2 | 887700 segments |
  | Y07 (unsound) | h′ halved; unsound in L_z where \|Re σ\| < 1 < \|σ\| | 1503 boxes |
  | Y16 (unsound) | `cert_box` starts from `cert_b(R1)`; \|σ_b\| up to 1.0089, b outside the disc | 2796 boxes |
  | Y21 (unsound) | 10⁻⁹ gaps at the splits in \|η\| | 892282 segments |
  | Y22 (unsound) | 10⁻⁹ gaps at the splits in angle | 892266 segments |
  | Y31 (reporting only) | CERT_BANDS from 1.893 | `v_lines` then says \|η\| ≤ 1.8927, r < 0.8052 |

- **Unsound but equivalent on the data (3):**
  - Y05: W without ℓ².
  - Y18: the end checked at 2ε₀; the flow stops at ε₀/2.
  - Y24: `v_lines` accepts contiguity within 10⁻².
- **Sound or equivalent by argument (6):**
  - Y04, Y08 and Y09: the redundant terms. Y09 alone cuts band 0 to 760 boxes and 256435 segments.
  - Y10: |1 + z|/2 replaced by ½ in L_η. This is compensated by the redundant 2(rB − 1)₊ on |σ| ≤ 2.02 (checked on 10⁶ points).
  - Y15: the R2 corners only; they are always the farthest.
  - Y29: ρ < 1 in place of ρ ≤ 1; C is continuous at 1.
- **The round's own 21:** 20 are killed. W19 (round for floor in `v_lines`) survives. It is equivalent on the data: round and floor agree at 1.995 (0.990046), at 1.98 and at the test's 1.902 (0.820609). They differ at 1.99 (0.9801 against 0.9802). The commit message says "all killed".

**(10) Literature** (web search; [snippet] unless marked)
- **The standard of proof.** The field holds that floating point is not a proof without rigorous error analysis.
  - ICM 2018 panel: "floating-point arithmetic alone, in the absence of rigorous error analysis, cannot constitute a proof"; interval arithmetic "keeps track of rounding errors automatically"; analysing errors a priori by hand is an alternative, but harder.
  - Rump: INTLAB's results are "correct with mathematical certainty" through IEEE-754 directed rounding.
  - Tucker: validated numerics "account[s] for both rounding and discretization errors in computer-aided proofs".
  - Arb: ball arithmetic with automatic rigorous bounds.
- **Precedents for a verified saddle-point step.** I found no validated steepest-descent certificate of a Gauss quadrature or Krylov error. The nearest is a 2023 FAU dissertation that validates saddle-point asymptotics (Morse functions on S²) with interval arithmetic in MATLAB [snippet, record only].
- **Explicit (non-validated) bounds for steepest descents** exist: Olver, Boyd 1993, Nemes [snippet, from search summaries].
- **The draft's weight.** The draft is right to say "floating point, not interval arithmetic" every time. "Certified" and "the limit holds … modulo both computer-assisted steps" read as stronger than the computation was. After my recheck, V can honestly carry "certified"; D's branch and bound cannot yet.

## Verdicts per part

- **(i), (ii), (iii)(a), (iii)(c) and (v):** unchanged, proved as before. P_j's form is still a sound sketch.
- **(iii)(b)**
  - **On |η| < 1 and on both axes:** proved.
  - **On |η| ≤ 1.8927:** proved modulo D's floating-point branch and bound, which I did not re-verify rigorously.
  - **On the closed annulus 1.8927 ≤ |η| ≤ 1.995:**
    - The analytic piece is proved.
    - The b → η piece is verified, rigorously so by my ball-arithmetic recheck of all 77,355,579 segment bounds, with no failure and margins below −10⁻¹² throughout.
    - The deformed-path argument is correct, with the hypotheses stated above. Ω′'s definition is stronger than what V shows and than what the proof needs (D3cy-1).
  - **So:** the limit holds on |η| ≤ 1.995, modulo the branch and bound. It is a conjecture on 1.995 < |η| < 2, where the method fails.
- **(iv):** correct. r < 0.9900 and σ*/σ < 199 (198.99 as printed) follow, with the same proviso. The proof's text is stale (D3cy-6).

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cy-1 | minor | Ω′ is defined by components of {F < F₀} ∖ B(η/2, rb) with rb < d < 2rb. (iii)(b) says "the condition holds as soon as a path from b = (1 + t)η/2, some 0 < t < (2 − \|η\|)/\|η\|, to η keeps F < F(η/2)". For t ≤ (2 − \|η\|)/(2\|η\|), b ∈ B̄. In 75 of band 0's 2454 boxes this happens for part of the box (d/rb ≥ 0.851). The proof uses neither the ball nor d > rb [own]. | D61 (iii)(b) and the proof: define Ω′ by the proof's own hypothesis: a chord through η/2 in the open disc with Re(φ″u²) < 0, below F(η/2) except at η/2, and paths from its ends to 0 and to η below F(η/2). Name the ball version and V's version as two sufficient conditions, and replace "the condition holds" by "the proof applies". Alternatively, in `cert_box`, require t > (2 − R1)/(2R1) (split those 75 boxes) and assert that the polyline avoids B (it stays ≥ 1.197·rb away in bands 1–2). |
| D3cy-2 | minor | "Certified" and "certificate" are used for a floating-point computation with an unanalysed slack. The float corner distances fall one ulp short of the rigorous ones in 24–27% of boxes, and the float π/2 leaves a 1.2·10⁻¹⁶ rad sliver at the imaginary axis in no box. My Arb recheck of all 77.4 M segment bounds, plus the end pieces, b and the tiling, finds 0 failures and a rigorous excess of at most 3.4·10⁻¹⁴ [own]. | exp_gh_peano.py: ship the ball-arithmetic recheck (python-flint in tools/requirements-pinned.txt), store its result in the raw, and print it in V. D61 Status and New: then "certified in ball arithmetic", and "modulo the branch and bound" alone. Otherwise say "verified in floating point" and state that the edges and the corner distances rely on the slack, with an a priori rounding bound. |
| D3cy-3 | minor | 16 of my 33 part-V mutants survive. Y02, Y03, Y07, Y16, Y21 and Y22 are unsound and change the recomputed bands (band 0: 2305, 887700 segments, 1503, 2796, 892282, 892266). Y31 drops the claim to 1.8927 on recompute. Y05, Y18 and Y24 are unsound but equivalent on the data. No test recomputes V against the raw. W19 of the round's own spec survives (the commit says "all killed") [own]. | tests/test_gh_peano.py: <br>• pin `cert_box` on three fixed boxes (one per band, one in the violated sector, one at the real-axis sliver), and `cert_band` on a thin slice, to the committed values of (ok, n, worst); <br>• oracle tests at each term's binding configuration: \|Re σ\| < 1 < \|σ\| (Y07), a segment crossing \|σ\| = 1 (Y03), η at the box's outer corner (Y02), ℓ ~ \|1 − m\| near z = 1 (Y05); <br>• capture z₀ = `cert_b(R2)` in `cert_box` (Y16); <br>• a forced-split tiling test (Y21, Y22); <br>• CERT_BANDS contiguous with DOM_BB[-1], and a 10⁻⁶ gap rejected by `v_lines` (Y31, Y24); <br>• a path ending at 1.5ε₀ rejected (Y18); <br>• W19 with R_v = 1.99. <br>Record Y04, Y08, Y09, Y10, Y15 and Y29 as equivalent; they are sound. |
| D3cy-4 | cosmetic | C's gradient bound h′ + 2(\|Re σ\| − 1)₊ is loose: \|∇C\|² = h′² − 4(\|a\| − 1)(1 − \|a\|/ρ²) ≤ h′². Dropping both (·)₊ terms (sound) gives 760, 2011 and 10658 boxes (160111, 946644 and 9581078 segments; 2 s, 11 s, 93 s), against 2454, 7338 and 41008. The cost lies at the real axis (98.3% of band 2's boxes within 5°; b's margin ≈ −t²(1 − Re c) ≈ −3.5·10⁻⁸). Beyond 1.995 the method fails even with the tighter bound (59275 boxes failing on [1.995, 1.9975]) [own]. | `cert_seg_parts`: use h′(\|σ\|) alone, and give the identity in the docstring. D61 Status: say the cost and the failure sit in a sliver at the real axis, where the domination branch and bound does not close either. That supports "an analysis at η = 2"; "sixfold" is one ratio, 5.6. |
| D3cy-5 | cosmetic | σ*/σ < 198.99 is a floating-point artifact: (1 + 0.99)/(1 − 0.99) = 198.99999999999983 floors to 198.99. The exact value from the floored r is 199.00, and 199.92 from r* [own]. | `v_lines`, `d_lines`: compute the ratio from the floored r in Fraction or Decimal, or quote "σ*/σ < 199". D61 (iv), Status, New, Numerics. |
| D3cy-6 | cosmetic | The proof of (iv) still says "(by it for r < r*, conditionally on P beyond)" and "inside \|η\| < 1.8927 when r < r*". HANDOFF.md §7–8 predates V: the valley condition "Open: the components certified", and `test_gh_peano.py` 32 tests [own]. | D61 proof (iv): "for r < 0.9900 with V, conditionally on P beyond" and "inside \|η\| ≤ 1.995 when r < 0.9900". HANDOFF.md: D‴.61's state and the test count (34). |

## Sources

- [Machine-Assisted Proofs, ICM 2018 panel (arXiv:1809.08062)](https://arxiv.org/pdf/1809.08062) [snippet]
- [Rump, Verification methods (Acta Numerica 2010)](https://www.tuhh.de/ti3/paper/rump/Ru10.pdf) [snippet, record only]; [INTLAB paper](https://www.tuhh.de/ti3/paper/rump/INTLAB2026final.pdf) [snippet]
- [Tucker, Validated Numerics (Princeton 2011)](https://press.princeton.edu/titles/9488.html) and [Tucker, ECM 2004](https://www2.math.uu.se/~watuc425/main/papers/ECM04Tucker.pdf) [snippet]
- [Johansson, Arb (arXiv:1611.02831)](https://arxiv.org/abs/1611.02831v1); [arblib.org](https://arblib.org/) [snippet]
- [Gómez-Serrano, Computer-assisted proofs in PDE: a survey (arXiv:1810.00745)](https://www.arxiv.org/pdf/1810.00745) [snippet]
- [FAU dissertation on Arnold's problem: Morse functions on the 2-sphere, saddle point method with interval arithmetic](https://fau.isle.flvc.org/node/57710) [snippet, record only]
- Explicit error bounds for steepest descents and Laplace's method (Olver; Boyd 1993; Nemes), via search summaries: [arXiv:1310.0166](https://arxiv.org/pdf/1310.0166) and [error bounds for Hermite asymptotics, Proc. Roy. Soc. Edinburgh A](https://resolve.cambridge.org/core/journals/proceedings-of-the-royal-society-of-edinburgh-section-a-mathematics/article/error-bounds-for-the-asymptotic-expansions-of-the-hermite-polynomials/8FF6A73E60EA4B11B491F9AB64E3D2DD) [snippet, record only]
- [Wang–Zhang, Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/abs/2312.07940) [snippet; already cited in D61]

## Lesson

1. **A floating-point certificate whose terms are monotone is one step from a rigorous one.** Rerun the same formula in ball arithmetic. Here 77 M bounds took about half an hour, and the rigorous excess over the float bounds was 3·10⁻¹⁴.
2. **Prove the tightest bound you can, then test each term where it binds.** A redundant term in a Lipschitz bound multiplies the cost (by 4–6× here) and masks unsound mutants (it made Y10 equivalent).
3. **When a computation shows a weaker condition than a definition states, restate the definition as the proof uses it.** Do not claim the stronger one.
4. **A certificate's tests must recompute it against the stored result.** Otherwise an unsound change that only makes it cheaper passes unseen.

My files are all under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref104:
- own/: rerun_band.py, arb_verify.py, arb_stream.py, identity.py, attack.py, seg_attack.py, tight.py, beyond.py;
- mut/: spec104.py, recompute.py, run_rc.sh;
- logs/: band0–2.log, arb0.log, arb1.log, arb2_stream.log, identity.log, attack.log, seg_attack.log, tight.log, beyond.log, mut104_a.log, mut_own.log, rc104.log, pytest.log.
