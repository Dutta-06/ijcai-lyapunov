# Referee report D3dl, round 117: Proposition D‴.61 as revised after round 115

## Scope

**Commit and setup.**
- Commit eb22ec4 on claude/vigilant-clarke-hf939b, exported with `git archive eb22ec4 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref117.
- The sha256 of exp_gh_peano.py, test_gh_peano.py, D61.md and gh_peano_raw.npz all equal `git show eb22ec4:…`.
- W/repo/venv is a symlink to the repository's venv (python-flint 0.9.0, numpy 2.5.4, scipy 1.18.1).
- GP_OUT pointed at W/out, which holds copies of the raw and the md. Tests and mutant slots read width_mismatch_raw.npz from their own copy.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. Pytest always ran with `--basetemp` inside W, and at most two heavy processes of mine ran at a time.
- I used only read-only git commands; `git status` is clean. I touched no process of anyone else and killed none.

**Read in full.**
- review/drafts/D61.md, especially (iii)(b)'s Status, "What (b) still lacks" with Lemma VC and its proof, and the Numerics' (VC).
- experiments/exp_gh_peano.py:
  - every hunk since 5be79c1: `_ball(prec)`, `VC_PREC`, `_cert_E_ball`, `corner_psi1`, `corner_psi_mid`, `corner_mid_path`, `part_VC`'s keys coef/far/mid/path, `vc_lines`;
  - the functions they call: `cert_E`, `corner_psi`, `corner_floor`, `corner_phi_dom`, `corner_sector`, `gap_closed`/`bound_sup`, `dom_bb`/`box_bound_ball`.
- results/chaos/gh_peano.md.
- tests/test_gh_peano.py (the round-115 test line by line, and the part-VC and md tests).
- tools/mutants/round115_D61_own.py.
- review/REVIEW_round115_D3dj.md and the Round 115 entry of REVIEW_LOG.md.
- PAPERS.md Paper 18 and HANDOFF.md §§5, 7 (and §§0, 10).
- `git diff 5be79c1 eb22ec4` of these files.

## Summary

**Design.** Round 115 asked for a sharper account of what part VC leaves open. The revision adds:
- the first correction Ψ¹_φ in closed form, checked in Arb;
- the floor's failure far out (ε = 10⁻³, φ = 20, 25, 30);
- the θ = √εψ limit Ψ̃_ψ, checked in Arb, with an explicit path y = min(sx, ψ), s = 2ψ/(1 + 2ψ²), and *Lemma VC* proving the path negative;
- the matching limits Ψ_∞ and the fixed-θ limit with their floors;
- a statement that four regimes each have a proved negative path, and that "what remains is uniform".

Tests now recompute part VC. The other edits are 13 mutants and HANDOFF/PAPERS text. No proved claim changes.

**Verdicts.**
- **Reproducibility is exact.**
  - The md regenerates byte for byte from the raw (sha256 459bd64b…df25e1).
  - 39 tests pass (18.1 s).
  - `part_VC()` equals the raw (0.14 s).
- **All new mathematics is correct.** I re-derived it independently, recomputed it in Arb from the definition of the bound (not from `cert_E`), and verified it line by line:
  - Ψ¹_φ and its ratio on the floor;
  - Ψ̃_ψ, including the band where |σ| > 1 ≥ Re σ;
  - Lemma VC;
  - Ψ_∞ and the fixed-θ limit, with their floors.

  I could not break Lemma VC.
- **"Four regimes, each limit's negative path proved" is accurate for the inner limits.** The four regimes do cover [0, ε₀] × [0, π/2] once two of the charts take a reciprocal parameter as an interval including 0. The account of what remains has three problems:
  - it omits the junction with V's chord, where three of the four limits vanish identically;
  - it treats as four paths what is one explicit path family;
  - it calls the per-regime certificate wholly open, although the Ψ̃ chart's inner part closes in ball arithmetic in seconds (my prototype).
- **Mutants.** All 13 shipped mutants are killed and the control survives. Of my 19, 14 are killed. One non-equivalent survivor shows that `_cert_E_ball`'s (Re σ − 1)₊² branch is never exercised.
- **No major issue.**

**Findings.**
1. **D3dl-1 (minor). The junction from b is missing from "what remains".**
   - V's chord and 0-side lie on the real z-segment and are proved exactly.
   - But Ψ̃_ψ(x, 0) = Ψ_∞(x, 0) = 0 inside the disc, and the small-θ fixed-θ form is 0 at ŷ = 0. So b sits on the zero level at the regime's scale, and its margin is the next order.
   - Round 115 noted this in its body; it is in neither D61, the md nor HANDOFF.
2. **D3dl-2 (minor). The regimes' paths are one family, and the overlaps are by construction.**
   - The family is Im z = min(κ·Re z, θ), with κ = Im A/|Re A| ≈ 2θ/(ε + 2θ²) and A = η²/4 − 1. In Ψ̃ coordinates this is exactly y = min(sx, ψ); it tends to the Ψ_φ/Ψ_∞ floor as ψ → 0 and to the fixed-θ floor as ψ → ∞.
   - The Ψ_φ floor is valid precisely for θ < √(ε/2) (ψ < 1/√2), not merely for θ ≪ √ε.
   - "A path per regime plus the regimes' overlaps" should read "one path, four parameter charts".
3. **D3dl-3 (minor, constructive). The Ψ̃ chart's inner part is certifiable now.**
   - I wrote a cancellation-free form in Arb with ε ∈ [0, 0.005] as an interval including 0.
   - It certifies ψ ∈ [10⁻³, 5] and x up to 4–60 (depending on ψ) with no failures, in 0.1–11 s per run.
4. **D3dl-4 (minor). Test gaps.**
   - X02 survives: no Ψ̃ check point has Re σ > 1.
   - X09 and X10 survive: the path grid's x-range is not in the raw.
5. **D3dl-5 (cosmetic). md wording.** "0 < x", "146 ψ", "at most −0.0044", "D² < 2ψxy + y⁴/2", "by cases in D61, (iii)(b)", and "(Arb, 320 bits)" with Ψ in float. The ε labels are printed from code constants.
6. **D3dl-6 (cosmetic). Residue and an overstatement.**
   - The md and HANDOFF §5 say round 115's prototype "closes the φ regime … for φ ≲ 7". It enclosed 7 sample boxes.
   - HANDOFF §7 still has the "straight segment covers φ < φ₁ … plus the saddle's degeneration" plan.
   - D61's Numerics (VC) omits round 115's numbers.
   - The matching chart's E ≈ −X³/3 holds for X ≫ max(ε, θ) only.

## What I checked and how

**(1) Reproducibility** [own][md][raw]
- `--from-raw` into W/out: `cmp` identical to the committed md. Both sha256 are 459bd64bb0ab02da4d5b5f3f8dc132d3cdf65cfe3f979ed20ed7437d10df25e1.
- `pytest -q -p no:cacheprovider --basetemp=W/pt/base1 tests/test_gh_peano.py`: 39 passed in 18.10 s (W/pytest_gh.log). HANDOFF §§0 and 10 say 39 for this file.
- **Part VC recomputed** (W/logs/rerun_vc.log):
  - `GP.part_VC() == raw['VC']` is True, in 0.14 s; `vc_lines` of the recomputation equals that of the raw;
  - φ_dom = 1.4396926208478287, and the sector lo values are 0.014842…, 0.007307…, 0.001443….
  - `part_VC` leaves flint's global precision at 320 bits. This is harmless: every other ball user sets its own.

**(2) The first correction** [own] (by hand, then Arb from the definition F_η(σ) − F_η(η/2), 400 bits; W/own/vc_check.py, W/logs/vc_check.log)
- **Derivation.**
  - Inside the disc the bound is Re((c − 1)z²) + Re(log(1 − z²) + z²), with c − 1 = (1 − ε/2)²e^{2iεφ} − 1 = ε(−1 + 2iφ) + ε²(¼ − 2φ² − 2iφ) + O(ε³).
  - The log term is −Re z⁴/2 + ….
  - So, with z = εw, the bound is ε³Ψ_φ + ε⁴[Re((¼ − 2φ² − 2iφ)w²) − Re(w⁴)/2] + O(ε⁵), which is the stated Ψ¹_φ.
- **On the floor.** With w = x(1 + 2iφ):
  - Ψ_φ = −(1 + 4φ²)x²;
  - Re((¼ − 2φ² − 2iφ)(1 − 4φ² + 4iφ)) = ¼ + 5φ² + 8φ⁴;
  - so Ψ¹/Ψ = −(¼ + 5φ² + 8φ⁴ − x²Re((1 + 2iφ)⁴)/2)/(1 + 4φ²) ≈ −2φ²(1 − x²) for large φ. The md is correct.
- **In Arb**, with Ψ_φ exact and ε = 10⁻⁶, 10⁻⁸, 10⁻¹⁰, the coefficient equals the closed form to at most 4.1·10⁻¹¹ relative at ε = 10⁻¹⁰:

  | φ | x = 0.05 | x = 0.25 | x = 0.45 |
  |---|---|---|---|
  | 0 | −0.24875 | −0.21875 | −0.14875 |
  | 10 | −200.2575 | −188.4669 | −160.9555 |

  - The repository's ε = 10⁻⁸ column is off by 7.8·10⁻⁸ relative at φ = 0 (−0.248750019 against −0.24875). It uses Ψ_φ and y in float and keeps only midpoints. The printed 6 digits are unaffected.

**(3) The θ = √εψ limit** [own]
- **Inside part.** Write Re(c − 1) = −ε(1 + 2ψ²) + O(ε²) and Im(c − 1) = 2√εψ + O(ε^{3/2}). With z² = ε²x² − εy² + 2iε^{3/2}xy this gives ε²[(1 + 2ψ²)y² − 4ψxy], and −Re z⁴/2 gives −ε²y⁴/2.
- **Outside part.**
  - |σ| − 1 = ε(x − ½ + y²/2) + O(ε²) = εD + O(ε²), and Re σ − 1 = ε(x − ½ − ψ²/2 − ψy) + O(ε²) = εP + O(ε²).
  - C = h(|σ|)·1{|σ| > 1} − (Re σ − 1)₊², with h(1 + δ) = δ² − δ³/3 + ….
  - When |σ| > 1 ≥ Re σ, C = h(|σ|) = ε²D² + O(ε³) and P ≤ 0, so the limit's D₊² − P₊² is D². Re σ > 1 forces |σ| > 1 exactly, so P > 0 implies D > 0. Both positive parts are C¹, so the limit is uniform across both switches.
  - Since cert_E(z̄, η̄) = cert_E(z, η), the expansion runs in integer powers of ε, and the remainder is O(ε).
- **Arb check at seven (ψ, x, y)** (ε = 10⁻⁴ … 10⁻¹⁰). In every case (B/ε² − Ψ̃)/ε is constant to 4 digits. The points include:
  - D > 0 > P: (0.5, 0.6, 0.2), with −0.16640000 against −0.1664, and (1, 0.45, 0.4);
  - P > 0: (0.2, 3, 0.1), with −0.00582514 at 10⁻⁸ against −0.005825 and coefficient −14.47;
  - the md's three points, verbatim (for example −0.89068410, −0.89062559, −0.89062501).
- **Uniformity.** The limit is uniform with O(ε) error on compact (x, y, ψ)-sets with ψ ∈ [ψ₀, ψ₁] ⊂ (0, ∞). Along the path it is relatively uniform near the saddle, with relative error about ε/(4ψ²) = 1/(4φ²): B/(ε²Ψ̃) = 1.0024 at ψ = 0.1, ε = 10⁻⁴, x = 0.01. It is not uniform as x → ∞:

  | ψ | ε | x | B/(ε²Ψ̃) |
  |---|---|---|---|
  | 0.1 | 10⁻⁴ | 10 | 4.84 |
  | 0.1 | 10⁻⁴ | 100 | 3392 |
  | 10 | 10⁻⁴ | 100.5 | 0.9869 |

  - A ratio above 1 is favourable. The next order on y = ψ is ε[(ψ² − 1)x² − x³/3 + O(x)], which is unfavourable at intermediate x for ψ > 1.
  - So the inner chart must stop at a finite X and hand over to the matching chart.

**(4) Lemma VC, line by line** [own]
- Every step holds:
  - D − P = (y + ψ)²/2;
  - on y = sx the quadratic part is −2ψxy, and P has slope 1/(1 + 2ψ²) and value −ψ²/2 at x₁ = ½ + ψ²;
  - the case D ≤ 0;
  - the case D > 0, u ≤ 0, where 0 < D ≤ a gives D² ≤ y⁴/4 < y⁴/2;
  - the case 0 < u ≤ ψ², where D² − y⁴/2 = 2u² − (u − a)² ≤ 2u², and u² ≤ ψ²u < 2ψ²u + ψ²/2 is the needed inequality, strict;
  - on y = ψ: D ≥ 3ψ²/2, convexity on [x₁, x₂], the end values −ψ² − ψ⁴/4 and −ψ² − ψ⁴/2, and D² − P² = 2ψ²(2x − 1 − ψ²) giving −ψ² − ψ⁴/2 beyond.
- The pieces agree at x₁ and x₂. Strictness holds for x > 0 and ψ > 0.
- **Trying to break it.** In exact rationals I evaluated 126 adversarial points: ψ from 10⁻⁶ to 10⁴; x₁ ± 10⁻¹² and x₂ ± 10⁻¹²; the D = 0 crossing ± 10⁻⁹; x = 10⁹·x₂. All are negative.

**(5) Matching limits and floors** [own]
- **Ψ_∞ from Ψ_φ (exact identity).** Ψ_φ(x, φỹ)/φ² = Ψ_∞ − (x² + (x − ½)₊³/3)/φ², checked to 3.6·10⁻¹⁵. The correction is non-positive, so regime 1's floor inherits regime 2's margin.
- **Ψ_∞ from Ψ̃.** Ψ̃_ψ(x, ψỹ)/ψ² = Ψ_∞ + O(ψ²), uniformly, including across x = ½.
- **Ψ_∞'s floor.** Ψ_∞ is convex in ỹ, with minimiser ỹ = 2x (value −4x²) for x ≤ ½ and ỹ = 1 (value −1) for every x > ½.
- **The fixed-θ limit.**
  - It is Re((e^{2iθ} − 1)(x + iY)²) + (x − ½)₊², with floor Y = x·cot θ and value −2x² + (x − ½)₊².
  - On the true bound at θ = 0.3 and 1.0, ε = 10⁻⁹: −0.08, −0.5, −1.75, −11.75 at x = 0.2, 0.5, 1, 3, equal to the formula to 6 digits.
  - Its small-θ form is the md's 2ŷ² − 4xŷ + (x − ½)₊², with floor ŷ = x and values −2x², and −x² − x + ¼ beyond.

**(6) Coverage, overlaps, and where each path suffices** [own] (W/own/composite_path.py, W/logs/composite_path.log)
- **Coverage.** Four closed parameter boxes cover [0, ε₀] × [0, π/2]:
  1. φ = θ/ε ∈ [0, Φ];
  2. (u, ψ) = (ε/θ, θ/√ε) ∈ [0, 1/Φ] × [0, ψ₀], with ε = ψ²u² substituted and B/(εθ²) as the scaled bound;
  3. ψ ∈ [ψ₀, ψ₁];
  4. (v, θ) = (√ε/θ, θ) ∈ [0, 1/ψ₁] × [0, π/2].

  Every (ε, θ) falls in one of them. Because the path is chosen per η, the overlaps need no matching of paths, only closed boxes. The fourth chart is the rest of the annulus, not a corner chart.
- **Where each regime's path suffices.**
  - *The Ψ_φ floor* y = min(2φx, φ) is, in Ψ̃ coordinates, y = min(2ψx, ψ). There Ψ̃ = 4ψ²x²(2ψ² − 1) − 8ψ⁴x⁴ + D₊² inside the disc, which is negative iff ψ < 1/√2 (bisection: 0.70711). So the floor works asymptotically exactly for θ < √(ε/2).
  - This explains the md's far-out values. At ε = 10⁻³, φ = 20 gives ψ = 0.63 (floor negative) and φ = 25 gives ψ = 0.79 (floor positive). My Arb values match ε²Ψ̃ to 3 digits: +1.544·10⁻⁹ against +1.543·10⁻⁹, and −8.070·10⁻¹⁰ against −8.080·10⁻¹⁰.
  - Regime 1's own expansion controls that floor only while 2εφ²(1 − x²) < 1, which is the same threshold to first order.
  - *The Ψ̃ path* works for every ψ > 0 in the limit, and on compact ψ-ranges uniformly.
- **The paths are one family.** In z, the Ψ̃ path is Im z = min(2θ·Re z/(ε + 2θ²), θ). With the exact κ = Im A/|Re A| this is the exact quadratic floor (HANDOFF §7) capped at Im z = θ. It reduces to the Ψ_φ floor (cap Im z = εφ) and to the fixed-θ floor (κ → cot θ).
- **On the true bound** (Arb, 300 bits), this family is negative at 12 parameter points for each of ε = 10⁻⁴, 10⁻⁶, 10⁻⁸ — φ = 0, 0.5, 5; ψ = 0.3, 0.7, 1, 3, 10; θ = 0.05, 0.3, 1, 1.5 — on x ∈ [0.01, 20].
- **The junction (D3dl-1).** On the vertical segment from b = (x_b, 0) to the path, Ψ̃ = (1 + 2ψ²)y² − 4ψx_b·y − y⁴/2, and ∂_yΨ̃ = 2(1 + 2ψ²)y − 4ψx_b − 2y³ < 0 for 0 ≤ y < sx_b. The value at y = 0 is 0 at scale ε². The margin there is V's exact real-axis bound, of order ε³.

**(7) A uniform certificate per regime, and the Ψ̃ prototype** [own] (W/own/chart3_ball.py, chart3_validate.py; W/logs/chart3_*.log)
- **What each chart needs:**
  - a chart in which the scaled bound is continuous in its small parameters at 0, with every division by ε, √ε or θ done analytically (sinc-type factors, the log and h series with rigorous tails, and |σ| − Re σ = (Im σ)²/(|σ| + Re σ));
  - the path family above;
  - branch and bound over parameter × path boxes, with the small parameters as intervals including 0, and splitting in ε as well;
  - near the saddle, division by the path's vanishing order x²;
  - the junction from b;
  - a hand-over to a matching chart. On the cap the bound equals log(1 − X) + X + X²/2 ≈ −X³/3 to ratio 1.0000 for X ≫ max(ε, θ); at ψ = 10, X = 10⁻³ the ratio is 16.1. So that chart has small parameters ε/X and θ/X.
  - Then an outer chart and V's end piece.
- **The Ψ̃ chart's forms** (t = √ε, θ = tψ, z = εx + ity):
  - A_r = Re(c − 1)/ε = −1 + ε/4 − 2(1 − ε/2)²ψ²sinc²θ;
  - A_i = Im(c − 1)/t = 2ψ(1 − ε/2)²sinc 2θ;
  - Q = A_r(εx² − y²) − 2A_i·xy. On y = sx, substituting (1 + 2ψ²)s = 2ψ analytically: Q = x²[εA_r − s(2ψ + Rs + 2(A_i − 2ψ))], with R = A_r + 1 + 2ψ² = ε/4 + 2ψ²(1 − (1 − ε/2)²sinc²θ);
  - L = −Re(Z⁴G(z²)), with Z = tx + iy and G(v) = Σ_{k≥2}v^{k−2}/k plus a tail;
  - D̃ = (1 − ε/2)(2x + y² + εx²)/(|1 + z| + 1) − ½;
  - P̃ = x − ½ − εx/2 − (1 − ε/2)[(ψ²/2)sinc²(θ/2)(1 + εx) + ψ·sinc θ·y];
  - C/ε² = D̃₊²(1 + εD̃₊H₁(εD̃₊)) − P̃₊², with H₁(δ) = Σ_{k≥1}(−δ)^{k−1}·(−1)/(k + 2)… plus a tail;
  - on y = ψ with P̃ > 0: A_r(εx² − ψ²) + 2xψ²[K₁²/(|σ| + Re σ) − 2K₃] + 2J(D̃ − x) − J² + εD̃³H₁ + L, where J = ψ²K₁²/(|σ| + Re σ) = D̃ − P̃ exactly;
  - sinc and cos on θ-intervals are taken by monotonicity. Arb's sinc of a ball through 0 is loose: sinc(0.5 ± 0.6) = [± 1.26].
- **Validation.** At thin parameters (ε = 10⁻², 10⁻³, 10⁻⁵ × 24 points) the forms match the definition to at most 7.9·10⁻¹⁵ relative, within the radius.
- **Certified with ε ∈ [0, 0.005]** (the whole last annulus), none failed:

  | ψ | x range | boxes |
  |---|---|---|
  | [10⁻³, 10⁻²] | ≤ 4 | 2675 |
  | [0.01, 0.05] | ≤ 4 | 837 |
  | [0.05, 0.2] | ≤ 4 | 372 |
  | [0.2, 0.5] | ≤ 4 | 271 |
  | [0.05, 0.5] | ≤ 60 | 2332 |
  | [0.5, 2] | ≤ 8 | 2436 |
  | [0.5, 2] | ≤ 30 | 8941 |
  | [0.5, 2] | ≤ 60 | 11742 |
  | [2, 3] | ≤ 15 | 3995 |
  | [3, 5] | ≤ 40 | 28694 (11.3 s) |

  - With ε ∈ [0, 10⁻⁴], ψ ∈ [4, 8] and x ≤ 100 certify in 3030 boxes.
  - The smallest ε-box was 1.5·10⁻⁷, at ψ ≈ 10⁻³: the margin ψ² against O(ε) terms.
  - At x ≤ 100 with ε₀ = 0.005 the run leaves the series' domain (εD̃ ≥ 0.5). That is the matching chart's job.
- **Not covered by the prototype:** the junction, the matching and outer charts, ψ < 10⁻³ (charts 1 and 2), and ψ > 5 at ε ≤ 0.005 (chart 4).

**(8) Claims against the md and raw** [md][raw]
- **Verbatim and earned:** D61's −200.257 at φ = 10, "from φ = 25", ψ ∈ [0.02, 20], −ψ² − ψ⁴/2; the md's coefficients, far-out values, mid values and path summary (raw: worst −0.004365644, far 8.946·10⁻¹³); φ_dom = 1.4397; PAPERS' list.
- **Wording slips:** see D3dl-5 and D3dl-6.
  - Round 115's prototype (ref115/own/chart1_ball.py, read-only) evaluated 7 fixed boxes of half-width 0.01. It did not run a branch and bound over the φ regime.

**(9) Tests and mutants** [own] (W/mut/; W/tools_ref/mutants_bt.py is tools/mutants.py with REPO fixed to W/repo and a per-slot basetemp)
- **Shipped round115_D61_own.py.** W00 (control) survives; W01–W13 are all killed by the tests. W09, W11 and W12 also change the md on regeneration.
- **My 19 mutants and a control** (W/mut/spec117.py; logs log_own_0.txt and log_own_1.txt):
  - Killed (14): X01 (`_ball` ignores prec), X04 (h halved), X05 (w⁴/4), X06 (P₊³), X07 ((1 + ψ²)y²), X08 (coef at 10⁻⁵), X11 (floor positive taken from the least), X12 (negativity tested against −1), X13 (θ = εψ), X14 (least by max), X15 (min(x, 2)), X16 (x = 0.25 in the ball), X17 (VC_MID_EPS 10⁻³; also changes the md on --from-raw, because the labels come from the constant), X18 (xmax).
  - Survive, non-equivalent (3): X02 (the (Re σ − 1)₊² term dropped), X09 (path grid from 10⁻²) and X10 (grid to max(20, …)). These are D3dl-4.
  - Survive, equivalent on the data (2): X03 (|Re σ| → Re σ) and X19 (no straddle guard).
- **The round-115 test.** I read it line by line:
  - it recomputes `part_VC()`;
  - it pins φ_dom to 10⁻⁹ and checks the closed ratio against ball[1] to 10⁻⁶;
  - it checks the far floor, the mid values to 10⁻⁵, path.worst and path.far, the inside closed form at ψ = 0.3 and 4 (inside the disc only);
  - it checks Lemma VC in exact rationals at 3000 points, and the three matching limits with tolerance 10⁻³;
  - it plants a ball bound of 5·10⁻⁹/10 and requires `dom_bb` to accept none. It does exercise acceptance: W11 is killed.
  - An unused local `w` sits in the coef loop.

**(10) Literature** [snippet] (web search; only what I saw)
- **Precedents for D61's route, and a contrast.**
  - The precedent for the deformed path, through saddles for entire integrands: Chawla–Jain 1968 (AMS PDF).
  - Gauss–Hermite root-exponential convergence: arXiv:2312.07940.
  - Gaussian rules uniform under coalescing saddles (Huybrechs–Kuijlaars–Lejon, arXiv:1806.06549). These are uniform expansions, not validated valley inequalities.
- **The regime structure.** It matches the classical Kruskal–Newton diagram and distinguished limits of dominant balance (Wikipedia, "Method of dominant balance"; arXiv:1807.03107, where leading balances are faces of a Newton polytope). Here Re A ≈ −(ε + 2θ²) balances at θ ~ √ε, and the outside term balances at θ ~ ε.
  - Quasi-homogeneous multi-order expansions: arXiv:2211.06865.
- **The "factor ε out before the interval check" technique** for ε ∈ (0, ε₀]: arXiv:1909.06207. I found no validated certificate of a degenerate saddle meeting a branch point.

## Verdicts per part

- **(i), (ii), (iii)(a), (iii)(c), (v):** unchanged by this diff; proved as before (P_j's polynomial form still a sketch).
- **(iii)(b):**
  - **|η| ≤ 1.995:** unchanged; proved, computer-assisted, both steps in ball arithmetic.
  - **1.995 < |η| < 2:** conjecture P, correctly labelled.
  - **Part VC:** Ψ¹_φ and its ratio, Ψ̃_ψ, Lemma VC, Ψ_∞, the fixed-θ limit and all the floors are correct, and the printed numbers are verbatim. "Four regimes, each limit's negative path proved" is accurate for the inner limits. The list of what remains omits the junction (D3dl-1), overstates the per-regime work (one path family; overlaps by construction; D3dl-2), and is staler than the Ψ̃ chart's actual state (D3dl-3).
- **(iv):** unchanged: r < 0.8052 and r < 0.9900 (σ*/σ < 9.26 and 199.00).

## Status of the round-115 issues

| id | status |
|---|---|
| D3dj-1 (VC's remaining gap) | **Resolved.** Ψ¹_φ, the far floor, Ψ̃_ψ with its path and Lemma VC, the matching limits, "φ″ → 0 harmless" and "no overlap with the straight segment" are all in. Residue: the junction from b, which round 115 noted in its §3, is not carried over (D3dl-1); the prototype's reach is overstated (D3dl-6). |
| D3dj-2 (tests) | **Resolved.** `part_VC()` is recomputed, φ_dom is pinned to 10⁻⁹, the plants are 1.9999/1.8928/1.5, and all 13 shipped mutants are killed (rerun). Residue: X02, X09, X10 (D3dl-4). |
| D3dj-3 (the straight segment's maximum) | **Resolved** in `vc_lines`, D61's Numerics and `corner_phi_dom`'s docstring. |
| D3dj-4 (float's limit) | **Resolved.** `corner_sector`'s docstring and the md say "valid for ε ≥ 10⁻³", and the new checks run in Arb. The far-out floats at ε = 10⁻³ agree with my Arb values to 3 digits. |
| D3dj-5 (HANDOFF) | **Largely resolved.** §5 is rewritten, §10's counts and ids are right, and the re-review line is dropped. Residue in §7: the "straight segment covers φ < φ₁ … descent path covers φ > φ₀ … plus the saddle's degeneration" sentence (contradicting §5), and "Write it as D‴.61 (iii)(b) on the whole disc, add the connectivity check to part D" (superseded by V). |
| D3dj-6 (odds and ends) | **Resolved.** `vc_lines` takes x and ε from the raw's rows and skips Ψ = 0; the dense-test tolerance is 10⁻¹⁴; the planted ball bound is rejected. Residue: the coef and mid ε labels still come from code constants (X17). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dl-1 | minor | **The junction from b to the regime's path is missing from "what remains".** V's chord and 0-side are the real z-segment through b = (x_b, 0), x_b < ½, proved exactly. But Ψ̃_ψ(x, 0) = 0 and Ψ_∞(x, 0) = 0 inside the disc, and 2ŷ² − 4xŷ is 0 at ŷ = 0. So in regimes 2–4 the path starts on the limit's zero level, and the margin there is V's ε³-order real-axis bound. On the segment up to the path, ∂_yΨ̃ = 2(1 + 2ψ²)y − 4ψx_b − 2y³ < 0 for y < sx_b [own]. Round 115 noted it in §3; it is absent from D61, the md and HANDOFF [md]. | D61 "What (b) still lacks", `vc_lines`' last sentence, HANDOFF §5: add "the junction from b: B(x_b, 0) < 0 exactly (V), and ∂_yB < 0 near y = 0 in a cancellation-free form". |
| D3dl-2 | minor | **The regimes' paths are one family, and the overlaps are by construction.** Im z = min(κRe z, θ), κ = Im A/\|Re A\| ≈ 2θ/(ε + 2θ²), is exactly the Ψ̃ path. It tends to the Ψ_φ/Ψ_∞ floor as ψ → 0 and to the fixed-θ floor as ψ → ∞, and it is negative on the true bound at 36 (ε, θ) spanning all regimes, θ up to 1.5 [own]. The Ψ_φ floor fails precisely at ψ = 1/√2 (θ = √(ε/2)), which explains φ = 25 against φ = 20 at ε = 10⁻³ (Arb: +1.544·10⁻⁹ against ε²Ψ̃ +1.543·10⁻⁹) [own]. The four closed boxes listed in (6) cover [0, ε₀] × [0, π/2]; the path is per η, so no path matching is needed, and the fourth chart is the rest of the annulus. | D61 (Status, "What (b) still lacks"), `vc_lines`, PAPERS Missing, HANDOFF §5: "one path family (the exact quadratic floor capped at Im z = θ), four parameter charts — two with a reciprocal parameter as an interval including 0 —, closed overlapping boxes; the Ψ_φ floor valid for θ < √(ε/2)". |
| D3dl-3 | minor | **The Ψ̃ chart's inner part is certifiable now.** VC calls every per-regime certificate open. My Arb prototype (forms in (7); validated to 7.9·10⁻¹⁵) certifies ε ∈ [0, 0.005] for ψ ∈ [10⁻³, 5], x up to 4–60, in 271–28694 boxes per run, none failed, 0.1–11 s [own]. | Part VC: add the chart (a `corner_mid_cert` with the forms in (7)), its counts in the md and D61's Numerics; state what remains: ψ < 10⁻³ (charts 1 and 2), ψ > 5 (chart 4), the junction, the matching chart (small parameters ε/X and θ/X) and the outer chart. |
| D3dl-4 | minor | **Test gaps.** `_cert_E_ball`'s (Re σ − 1)₊² branch is never exercised: no VC_MID point has P > 0 and the coef points are inside the disc, so X02 survives. The far piece of the path is never checked against the true bound. X09 and X10 survive because the raw's `path` does not record the grid's x-range [own]. | exp_gh_peano.py VC_MID: add (0.2, 3, 0.1) (P = +2.46; Arb −0.00582514 at 10⁻⁸ against −0.005825, coefficient −14.47) and a far-piece point. Store xq's range in `path`. tests/test_gh_peano.py: test both. tools/mutants: add X02, X09, X10. |
| D3dl-5 | cosmetic | **md wording** [md][raw]. "0 < x" (the grid starts at 10⁻⁴). "146 ψ" (145 distinct; ψ = 1 twice). "At most −0.0044·ψ²min(x, 1)²" (raw −0.0043656, so "at most −0.0043"). "D² < 2ψxy + y⁴/2" (false near x = 0 where D < 0; should be D₊²). "By cases in D61, (iii)(b)" (the lemma is in the Status). "(Arb, 320 bits)" while Ψ_φ and y are float and only midpoints are kept (the 10⁻⁸ column is off by 7.8·10⁻⁸ relative). The coef and mid ε labels are printed from constants (X17). | `vc_lines` and `part_VC`: round toward 0, say "10⁻⁴ ≤ x" and "145", use D₊², point to the Status, compute Ψ_φ in Arb, and store the ε lists in the raw. |
| D3dl-6 | cosmetic | **Residue and an overstatement** [md]. The md and HANDOFF §5 say round 115's prototype "closes the φ regime in balls for φ ≲ 7 at ε ≤ 0.005"; it enclosed 7 fixed boxes (ref115/own/chart1_ball.py). HANDOFF §7 still has the φ₀ < φ₁ / "saddle's degeneration" plan and "add the connectivity check to part D". D61's Numerics (VC) omits round 115's numbers. The matching chart "E ≈ −X³/3" holds only for X ≫ max(ε, θ) [own]. | md and HANDOFF §5: "enclosed sample boxes of the φ chart below 0 for φ ≲ 7". HANDOFF §7: drop the stale sentences. D61 Numerics (VC): add the coefficients, the far values, the mid convergence and the grid. HANDOFF §5: give the matching chart's two small parameters. |

## Sources

- Chawla–Jain, Asymptotic error estimates for the Gauss quadrature formula, Math. Comp. 22 (1968): https://www.ams.org/mcom/1968-22-101/S0025-5718-1968-0223094-5/S0025-5718-1968-0223094-5.pdf [snippet]
- Gauss–Hermite root-exponential convergence, arXiv:2312.07940: https://arxiv.org/pdf/2312.07940v2 [snippet]
- Huybrechs–Kuijlaars–Lejon, Gaussian rules for coalescing saddle points, arXiv:1806.06549: https://arxiv.org/abs/1806.06549v2 [snippet]
- Rigorous numerics for a singular perturbation problem (ε factored out before interval verification), arXiv:1909.06207: https://arxiv.org/pdf/1909.06207 [snippet]
- Singular perturbations and scaling (dominant balances as Newton-polytope faces), arXiv:1807.03107: https://arxiv.org/pdf/1807.03107 [snippet]
- Method of dominant balance (Kruskal–Newton diagram, distinguished limits): https://en.wikipedia.org/wiki/Method_of_dominant_balance [snippet]
- Multi-order asymptotic expansion of blow-up solutions (quasi-homogeneous weights), arXiv:2211.06865: https://arxiv.org/pdf/2211.06865 [snippet]
- Johansson, Arb, IEEE Trans. Comput. 66 (2017), arXiv:1611.02831 [memory; as cited in D61]

## Lesson

1. **Write each regime's path back in the original coordinates before planning one certificate per regime.** Here the four paths are restrictions of one family: the exact quadratic floor Im z = (Im A/|Re A|)Re z, capped at Im z = θ. "Paths per regime and their overlaps" then becomes "one path, four remainder charts".
2. **A limit that vanishes on the chord is a junction to certify.** Ψ̃(x, 0) = 0 means the path's first step carries a lower-order margin, and only the limit's zero level shows it.
3. **Two things make a ball certificate with a small parameter as an interval including 0 work:** every division by ε, √ε or θ done analytically, and splitting in ε near 0. Arb's own sinc on a ball through 0 is too loose, so use monotonicity. With both, the Ψ̃ chart's inner part closed in seconds.

My files, all under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref117:
- own/: vc_check.py, chart3_ball.py, chart3_validate.py, composite_path.py, matching_form.py;
- logs/: vc_check.log, chart3_validate.log, chart3_run1.log … chart3_run5.log, composite_path.log, matching_form.log, rerun_vc.log;
- mut/: spec117.py, log_ship.txt, log_own_0.txt, log_own_1.txt, ship/res_a.jsonl, own/res_s0.jsonl, own/res_s1.jsonl;
- tools_ref/mutants_bt.py;
- pytest_gh.log;
- out/: the regenerated md and the raw copy.
