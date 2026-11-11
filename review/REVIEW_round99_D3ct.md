# Referee report D3ct, round 99: Proposition D‴.61 as revised after round 95

## Scope

I reviewed commit cace1fd on branch claude/vigilant-clarke-hf939b.

- **Setup**
  - **Work directory.** I exported the commit with `git archive cace1fd | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref99. GP_OUT = W/gp_out holds a copy of the raw.
  - **Run conditions.** Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` in the repository's venv, symlinked as W/repo/venv for tools/mutants.py. Pytest used `--basetemp` under W/pyt, one directory per mutant slot.
  - **Limits.** Every mutant and every heavy job ran under `ulimit -v` (4–7 GB) and `timeout`. Pytest had a 600 s timeout per mutant.
- **Processes.** I kept to at most two heavy processes at a time. One exception: a 0.7 s smoke test, and possibly a few seconds' overlap between two runs. The only process I stopped was my own `wait_log.sh` (PID 5032), after checking /proc/5032/cmdline.
- **Comparison commit.** The round-95 referee's commit is 5a171c7: the last commit before c72a571 that touched D‴.61's files. I diffed every file from there to cace1fd.
- **Read in full**
  - review/drafts/D61.md and results/chaos/gh_peano.md.
  - experiments/exp_gh_peano.py: all of part D, and the round's changes line by line (`quad_roots`, `bound_sup`, `gap_arg`, `box_bound`, `dom_bb`, `path_seeds`, `path_F`, `path_check`, `part_dpath`, `d_lines`).
  - tests/test_gh_peano.py.
  - review/REVIEW_round95_D3cp.md and the Round 95 entry of REVIEW_LOG.md.
  - tools/mutants/round95_D61_own.py and tools/mutants.py.
- **Repository untouched.**
  - My only git commands were archive, log, show, diff, status and rev-parse, plus tools/numcheck.py's own `git diff`.
  - During the review, HEAD moved to d0ec681 through three commits that are not mine, and the working tree gained edits to D‴.58/width_mismatch files that are not mine. No D‴.61 file changed.
- **Web.** I used search snippets only. www.ams.org and ir.iitd.ac.in do not resolve from this session (ENOTFOUND), so I stopped fetching. For that, ams.org would have to be allowed in the session's network settings.

## Summary

**Design.** The revision applies round 95.
- `quad_roots` computes the roots of `bound_sup`'s two quadratics by the stable formula. The branch and bound was rerun.
- The valley check's seeds now lie on the saddle's steepest-descent line (`path_seeds`), with one formula for F (`path_F`). `path_check` adds a coarse-grid guard and requires the chord [a, b] to stay below F(η/2) (481 samples, |t| ≥ 0.05), and the margins are reported.
- (iii)(b) restates Ω′: a and b on a line through η/2, within 2rb of it, in distinct components, with the chord below F(η/2) except at η/2.
- (iv)'s conditional hypothesis is now the annulus 1.8927 < |η| ≤ 1.98, with σ*/σ < 49.89.
- The alias continuation is labelled a heuristic, with the rate max(0, λ₁).
- `gap_arg` prints the violated sector and the maximiser's Re(ηs).
- Chawla–Jain 1968 is cited.
- There are 4 new tests and 30 own mutants.

**Verdicts**
- **Reproducibility.**
  - `--from-raw` regenerates the md byte for byte, and all 30 tests pass.
  - `dom_bb(1.78)`, `dom_bb(1.8927)` and `part_dpath` rerun bit-identically against the raw, as do `part_dom`'s angles, axes, violations and margins.
  - `part_dom`'s list of radii is stale: the raw still holds the round-93 roots' values (D3ct-5, cosmetic).
- **Numbers.** Every decimal of the draft is in the md verbatim. No overreach in substance.
- **Stable roots.** The D3cp-2 hole is closed. On 36 boxes chosen for their edges and 400 random dyadic boxes, no box bound falls below a dense supremum, and `bound_sup` matches an independent dense supremum of its own function to 4.4·10⁻¹⁶.
- **Steepest-descent seeds.** They are correct and correctly oriented.
- **Ω′ argument.** It is correct, and it needs essentially the hypotheses now stated. But "distinct components" now does real work: it is what puts a and b on opposite sides of η/2. And the proof uses a steepest-descent direction that the hypothesis does not ask for (D3ct-1).
- **The chord condition is provable, not just sampled.** A two-line Taylor bound gives it at all 372 annulus points I tried. That also justifies the 0.05 exclusion.
- **A certificate of the valley condition looks feasible.** Explicit monotone gradient-flow paths exist at 372 of 372 points. A first-order estimate needs about 1.5·10³ η-boxes, provided F − F₀ is enclosed as one centred form; naive intervals would need about 8·10⁹ (D3ct-2).
- **(iv)'s constants are right:** 0.9607, 49.89, 0.8052, 9.26, 0.0966.
- **Alias label.** It is correct where the alias is discussed, but one unlabelled clause remains. "Λ⁺ attained numerically" in the sector is not shown at K ≤ 340 (D3ct-3).
- **Mutants.** 19 of my 34 survive, and 3 of them change published results (D3ct-4). The round's own spec reruns as logged: 28 killed, 2 equivalent.
- **Round-95 issues.** All seven are resolved or acknowledged, with residues. No major issue.

**Findings**
- **D3ct-1 (minor). Ω′'s hypotheses against its proof.**
  - As now worded (two points on a line through η/2), "distinct components" is load-bearing: it excludes a and b on the same side of η/2.
  - The proof's "steepest-descent line" is not a hypothesis. Only Re(φ″u²) < 0 is needed, and the chord hypothesis gives that up to the degenerate case.
  - "Within 2rb" must be strict.
  - The chord condition holds by Taylor's theorem: F − F₀ ≤ −|φ″|t²/2 + sup|φ⁗|t⁴/24, since φ‴(η/2) = 0, with ratio ≤ 6.4·10⁻³. This should replace the 481-point sampling.
- **D3ct-2 (minor). The certificate is feasible, and the draft's recipe omits what decides its cost.**
  - Gradient-flow paths of F descend monotonically from the seeds to 0 and to η at all 57 points and at 372 of 372 points of a 12 × 31 grid of the annulus.
  - The η-boxes can be about 2.4·10⁻³ … 1.6·10⁻² wide only if F(σ; η) − F₀(η) is enclosed as one centred form around the moving saddle. Naive enclosures shrink δ to 10⁻⁷.
- **D3ct-3 (minor). Residue of D3cp-4.**
  - "Λ⁺ … is attained where |Re σ| > 1 …, the first alias's lower rate holding where |Re σ| < 1" is unlabelled.
  - It is wrong where λ₁ < 0 (the rate there is max(0, λ₁)), and λ₁ is undefined off the axis.
  - At the sector's maximiser, (1/K)log|c_K| grows at 0.0158 … 0.0204 for K = 60 … 340, against Λ⁺ = 0.0233, and the straight segment has not overtaken s = ½ by K = 340.
- **D3ct-4 (minor). Mutants.** 19 of 34 survive.
  - K25, K27 and K32 change published results.
  - K29 and K31 are unsound but equivalent on the data.
- **D3ct-5 (cosmetic). The raw's D.dom.radius predates `quad_roots`.** 113 of 361 radii differ from the committed code's, by up to 4.4·10⁻¹². The md is unaffected.
- **D3ct-6 (cosmetic). Rounding and storage.**
  - The maximiser's Re(ηs) has the lower end 1.0052: the grid minimum rounded to nearest. On the stated sector [0.0861, 0.6175] the infimum is 1.00514, and on the exact sector [0.08606, 0.61756] it is 1.00511.
  - The sector is computed at report time and not stored in the raw.
- **D3ct-7 (cosmetic). Wording and citation.**
  - "Ω′ is known only at those points" and (iv)'s "the condition holds at the 57 points" should say *numerically*.
  - Chawla–Jain treat Gauss–Legendre.

## What I checked and how

**Reproducibility** [own]
- **The md.** `--from-raw` gives an identical md (sha256 89396fa3…).
- **Tests.** `pytest tests/test_gh_peano.py`: 30 passed in 4.6 s.
- **Rerun** with the script's own functions (W/own/rerun.py), compared by `repr`:
  - **Bit-identical:**
    - `dom_bb(1.78)` = (4935, −3.281257330289833e-07, 0) and `dom_bb(1.8927)` = (402664, −1.9432233600014115e-11, 0), in 4 s;
    - `part_dpath()` (57 checks, 51 s);
    - `part_dom`'s th, axes, viol and margin;
    - `duhamel_job` at (12, 1.5, 2) and (12, 1.9, 1), and `qx_job(20)`.
  - **Not identical: `part_dom`'s radius.**
    - 113 of 361 radii differ, by up to 4.36·10⁻¹².
    - The least radius is unchanged: 1.8927345433723064 at 0.3447.
    - The 5a171c7 code reproduces all five of `part_dom`'s fields exactly. So the raw's `dom` was not recomputed after `quad_roots` changed `gap_closed`.
- **The raw from 5a171c7 to cace1fd** [raw]. The only change is D.dpath: each row gained ma, mb and ms, and `seeds` was added. D.dom.bb is identical, which is consistent with "rerun, unchanged" for the branch and bound.

**The draft against the md and the raw** [md][raw][own]
- **tools/numcheck.py from 5a171c7.** Every added number is in the md except 1968 and 97, the Chawla–Jain citation.
- **All 222 distinct decimals** (W/own/numall.py) are in the md. The exceptions are "Theorem 2.10", three arXiv ids, and 1.72/1.44/1.16, which the md prints as 1.72e+00 etc.
- **Semantics, recomputed by my own code:**
  - the seeds' margins −9.77·10⁻⁵ … −1.90·10⁻³, −1.17·10⁻⁵ … −4.62·10⁻⁴ and −7.32·10⁻⁷ … −7.27·10⁻⁵;
  - the chord's maximum at |t| = 0.05: mutant K27 moves it to |t| = 0.5 and multiplies it by exactly (0.5/0.05)² = 100;
  - 19 of 19 angles at each radius;
  - r*, 0.9607, 49.89, 9.26 and 0.0966 (below).

**The domination gap and the box bound** [own] (W/own/dom_own.py: a dense supremum over s from the definition g_η(s) + Λ⁺(ηs) − g_η(½), refined by bounded minimisation)
- **The closed form against the dense supremum.** At 300 random η and 16 special η (φ = 0, π/4, π/3, π/2 and 0.345; R up to 1.999), `gap_closed` agrees within ±8.9·10⁻¹⁶.
- **Boxes chosen for their edges (36):**
  - the old hole [1.5, 1.8927] × [0, π/4], now bounded by 0.935078;
  - boxes ending at or starting from π/4, of widths π/2² … π/2¹² (where p = R²cos 2φ ≈ 0);
  - boxes at π/3 (where the first quadratic's leading coefficient R²(4cos²φ − 1) vanishes);
  - R₁ = 1 and [1, 1 + 10⁻⁷];
  - the least radius from inside, and straddling it (dense supremum +6.4·10⁻⁶, bound +2.4·10⁻³);
  - B = R₂cos f₁ crossing 1;
  - f₂ = π/2;
  - a degenerate point box.

  On all of them the bound is never below a 17 × 17 lattice of dense suprema: the largest (dense − bound) is −6.0·10⁻¹¹, at the point box. `bound_sup` against my own dense supremum of the same function differs by at most 2.2·10⁻¹⁶.
- **400 random dyadic boxes** (angle edges j·(π/2)/2ᵏ, as `dom_bb` makes them): largest (dense − bound) −6.0·10⁻⁶; largest (own dense h − `bound_sup`) 4.4·10⁻¹⁶.
- **Analytically** (the soundness of `bound_sup`'s candidate set):
  - M is C¹ at 1, so the kink candidate is redundant (mutant K07 is equivalent).
  - A computed disc < 0 at a near-double root loses at most O(δ³), because h′ is ≥ 0 on both sides.
  - For qa < 0 there is exactly one positive root.
  - Below the 10⁻¹⁴ threshold the dropped root is ≈ −qb/qa, far outside (0, wm ≤ 0.48]; the threshold could be raised to about 1 (mutant K01, threshold 1.0, does change the md).

  No box bound is unsound.
- **Sector and maximiser at 1.95** (own bisection on the dense gap):
  - The exact sector is [0.086058, 0.617558]; the md's [0.0861, 0.6175] are the first and last of its 20001 grid angles.
  - The maximiser's Re(ηs*) is 1.00511 … 1.09995 on the exact sector and 1.00514 … 1.09995 on [0.0861, 0.6175].
  - The code's grid minimum is 1.005156, at φ = 0.617480, printed as 1.0052 (D3ct-6).
  - The minimum is at the sector's upper end. Re > 1 holds throughout, which is the substantive claim.
- **Constants** (own c_need and bisection):

  | radius | r* | floored | c_need(floored) | σ*/σ = (1 + r)/(1 − r) |
  |---|---|---|---|---|
  | 1.8927 | 0.8052030 | 0.8052 | 1.8926982 | 9.26710; 9.26694 from the floored r → 9.26 |
  | 1.98 | 0.9607245 | 0.9607 | 1.9799873 | 49.92237; 49.89059 from the floored r → 49.89 |
  | 1.0 | 0.0966496 | 0.0966 | 0.9998257 | — |

  - The least radius on 101 angles in [0.30, 0.40] is 1.8927327 at 0.343. That is consistent with the md's 361-angle grid value, 1.892735 at 0.345. The verified radius 1.8927 is below both.

**The deformed path, re-derived** [own]
- **The bound F.**
  - The exponent per K of the u-integrand is F with C = (1/K)log|c_K|.
  - Inside the disc, (a) gives 0 uniformly on compacts. In the band near |σ| = 1, Λ⁺ ≈ (1 − |σ|)² → 0, and F is continuous across the circle.
  - Outside the disc, the ρ₀-contour bound tends to Λ⁺ uniformly on compacts.
  - So F is a valid upper rate along any path, and exact on the chord.
- **Orientation.** For |η| < 2, arg(η² − 4) ∈ (π/2, 3π/2). So the steepest-descent direction u satisfies |arg u − arg η| < π/4: Re(uη̄) > 0 picks the branch unambiguously, and the seed "toward 0" really is a. My flows confirm a → 0 and b → η at every point.
- **Hypotheses.**
  - **Within 2rb (strictly).** This is sufficient to keep the chord in the open disc. It is not necessary: Laplace needs (a) only near η/2, and elsewhere the upper bound F < F₀ suffices.
  - **Chord below F(η/2) except at η/2.** This is necessary for Laplace's method on the chord.
  - **Distinct components.** As worded, this is necessary. If a and b lay on one ray from η/2 outside B, then [a, b] would lie outside B and, by the chord hypothesis, in {F < F₀}, so a and b would share a component. If instead a = η/2 − du and b = η/2 + du with rb < d < 2rb, distinctness follows from the rest: otherwise Cauchy's theorem would make the chord's nonzero Laplace term exponentially small.
  - **Non-degeneracy.** Re(φ″u²) < 0 is used ("along its steepest-descent line") but not assumed. The chord hypothesis gives ≤ 0, and the degenerate case is curable by a local deformation, but it should be stated.
- **The chord by Taylor.**
  - F − F₀ = Re[φ(η/2 + tu) − φ(η/2)] ≤ −|φ″|t²/2 + sup|φ⁗|t⁴/24, since φ‴(η/2) = 0.
  - On the ball |σ − η/2| ≤ D, |φ⁗| ≤ 12/(|η|/2 − D)⁴, because φ⁗ = −6/(η − σ)⁴ − 6/σ⁴.
  - Set q = D²·sup|φ⁗|/(12|φ″|). Then F − F₀ ≤ −(|φ″|t²/2)(1 − q(t/D)²) < 0 on the whole chord once q < 1.
  - q at D = 1.2rb (and at D = 2rb):

    | \|η\| | q at 1.2rb | q at 2rb |
    |---|---|---|
    | 1.9 | ≤ 5.8·10⁻³ | ≤ 1.8·10⁻² |
    | 1.95 | ≤ 2.6·10⁻³ | ≤ 7.4·10⁻³ |
    | 1.98 | ≤ 9.5·10⁻⁴ | ≤ 2.7·10⁻³ |

  - Over the 372-point annulus grid, q ≤ 6.4·10⁻³ at 1.2rb.
  - So the 0.05 exclusion is justified (indeed superfluous), and the chord condition is proved rather than sampled.

**The valley check, numerically** [own] (W/own/valley_own.py: my own F, seeds, labelling and flows)
- **Grid variants.** Every variant gives 19 of 19 at each of 1.9, 1.95 and 1.98:
  - seeds at 1.2rb, 1.5rb and 1.9rb on the descent line;
  - 1.2rb rotated by ±0.35 rad;
  - 4-connectivity;
  - pad 0.6 and 0.4;
  - h/2 at 1.95 (0.001) and at 1.98 (0.0005, also with 1.9rb seeds).

  The seeds' margins at 1.98 are −1.84·10⁻⁶ … −1.82·10⁻⁴ with 1.9rb seeds, and −5.6·10⁻⁷ … −5.6·10⁻⁵ when rotated.
- **Explicit paths.** Gradient flows of F from the 1.2rb seeds reach 0 and η at all 57 points.
  - On each path, max(F − F₀) equals the seed's own value: the paths descend monotonically.
  - They never come closer to η/2 than 1.2rb, are 0.92–1.14 long, and on the η-side reach C = Λ⁺ up to 0.78.
  - The same holds at 372 of 372 points of a 12 × 31 grid of 1.8927 ≤ |η| ≤ 1.98 × [0, π/2].
- **Feasibility of a certificate** (W/own/feas.py).
  - Along the flow paths, take the margin m = F₀ − F, the σ-gradient |∇_σF| and the η-sensitivity |∇_η(F − F₀)| at fixed σ.
  - A first-order (mean-value) enclosure needs |∇_σF|ℓ/2 + |∇_η(F − F₀)|δ ≤ m. Add the seed's connector, δ ≤ m_seed/|∇_σF|_seed, and the Taylor-certified chord.
  - **With a centred enclosure:**
    - δ = 2.4·10⁻³ … 1.6·10⁻². The least δ is on the real axis: 1.24·10⁻² at 1.8927, 7.1·10⁻³ at 1.94, 2.39·10⁻³ at 1.98.
    - That gives N_box ≈ 1.45·10³ boxes for the quadrant, with 35–50 segments per path: ≈ 5.8·10⁴ box-segment enclosures.
    - Even with a factor of 100 for overestimation, this is minutes with ball arithmetic.
  - **The decisive point.** Near the saddle, ∂_η(F − F₀) ≈ −(φ″/2)(σ − η/2) is small, because dF₀/dη = ∂_ηφ(η/2). An enclosure that bounds F and F₀ separately loses this cancellation: δ falls to 1.0·10⁻⁷ and N_box ≈ 7.9·10⁹. So the certificate must treat F − F₀ as one centred (Taylor-model) form around the moving saddle.
  - The band near |σ| = 1 needs 1 − r ≲ √m (about 10⁻³ at 1.98) and is harmless.

**Is Λ⁺ attained in the violated sector?** [own] (W/own/lam_att.py, the script's `corner_ratio`, 4K + 100 digits)
- **The point.** η = 1.95e^{iπ/8} has gap +7.07·10⁻³ and s* = 0.607504, with Re(ηs*) = 1.09446 and Λ⁺ = 0.023322.

  | K | 60 | 100 | 140 | 180 | 220 | 260 | 300 | 340 |
  |---|---|---|---|---|---|---|---|---|
  | (1/K)log\|c_K\| | −0.00559 | 0.00297 | 0.00675 | 0.00895 | 0.01050 | 0.01174 | 0.01279 | 0.01369 |

- **Slopes per 40 K:** 0.0158, 0.0162, 0.0167, 0.0175, 0.0185, 0.0196, 0.0204. They are accelerating, like a crossover from a lower-rate term to the edge term.
- **Fit λK + α log K + β on K = 220–340:** λ = 0.0261 with α = −1.82; with α = −3/2 fixed, λ = 0.0250. Both exceed Λ⁺, so this is still pre-asymptotic.
- **The straight segment.** At s*, g(s*) − g(½) = −0.016257. The integrand at s* is still below s = ½ at K = 340 (−0.0026 per K), though its local growth exceeds s = ½'s from K ≈ 200.
- **At |σ| = 1.25** (Re σ = 1.155): Λ⁺ = 0.034128; fit 0.0346; slopes 0.0293–0.0306.
- **Conclusion.** The data are consistent with Λ⁺ being attained and with eventual overtaking. Neither is observed at K ≤ 340.

**Mutants** [own] (W/mut/spec99.py: 34 single-point mutants of this round's new or changed code, plus a control; tools/mutants.py with two slots; survivors recomputed by W/mut/recompute.py against a recomputed control)
- **Killed (15):** K01, K04, K05, K06, K09, K10, K11, K13, K14, K15, K16, K18, K19, K30, K34.
- **Survivors that change published results on recompute (3):**

  | mutant | change | effect |
  |---|---|---|
  | K25 | seed pixel by `int()` (floor), not `round` | dpath 19, 19, 16: false negatives at 1.98 for θ = 0, 0.698, π/2, where the pixel falls into the ball; the guard h/√2 < 0.2rb is tight at 1.98 |
  | K27 | chord exclusion \|t\| ≥ 0.5 | the md's chord values become −2.4·10⁻⁵, −2.9·10⁻⁶, −1.8·10⁻⁷ |
  | K32 | `zip(DP_R, DP_H[::-1])` | `part_dpath` raises ValueError at 1.98 (h = 0.004), so part D cannot be recomputed |

- **Unsound but equivalent on the data (2):**
  - K29: the chord is accepted up to 10⁻³ above F(η/2); the test's bump is 0.1.
  - K31: `path_ok` returns the margin, a truthy float; the md names `path_ok`.
- **Equivalent (14):**
  - K02, K03, K07 (the kink) and K08: the recomputed dom and branch and bound are identical to the control's.
  - K12, K17 (Λ⁺ = 0 on the circle), K22 (F₀ via `path_F`: the raw changes in its last bits only) and K33 (round vs floor: still 0.9607).
  - K20 and K21: pad 0.3 suffices, and dpath is identical.
  - K23 (the ball halved), K24 (4-connectivity), K26 (probe at 0.2η) and K28 (inner half of the chord).
- **The round's own spec, rerun:** 28 killed; J12 and J13 survive (equivalent by F's symmetry about η/2 in the disc); control clean. This matches the log.

## Verdicts per part

- **(i): proved** (unchanged).
- **(ii): proved.** P_j's polynomial form remains a sound sketch.
- **(iii)(a): proved.**
- **(iii)(b)**
  - On |η| < 1 and on both axes: proved.
  - On |η| ≤ 1.8927: proved modulo the floating-point branch and bound. That is now free of the D3cp-2 hole, and I found no unsound box.
  - On Ω′: the theorem is correct, with the hypothesis wording of D3ct-1.
  - The chord step is now provable by Taylor.
  - Membership of Ω′ beyond 1.8927 is numerical (57 grid points; 372 explicit flow paths). It is certifiable at modest cost (D3ct-2).
  - Elsewhere near |η| = 2: conjecture.
- **(iii)(c): proved.**
- **The alias continuation:** a heuristic, checked numerically and labelled so in its own sentence. One clause is still unlabelled and partly wrong (D3ct-3).
- **(iv): correct.**
  - r* = 0.8052 and 9.26 are valid floors; r < 0.0966 needs no check.
  - The conditional r < 0.9607 (σ*/σ < 49.89) has the right hypothesis: the annulus 1.8927 < |η| ≤ 1.98, closed at 1.98, which makes the disc |η| ≤ 1.98 a compact subset of Ω′.
- **(v): proved** (unchanged).

## Status of the round-95 issues

| id | status |
|---|---|
| D3cp-1 | Resolved as asked for meanwhile. The seeds are on the descent line, the guard is asserted, the margins are printed, the chord is required, and the Status says "numerical, not certified" and what a certificate needs. The certificate itself was not done (acknowledged). Residues: the refinement is not reported in part D (mine passes at h/2); the chord can be proved by Taylor rather than sampled; the recipe omits centring (D3ct-1, D3ct-2). |
| D3cp-2 | Resolved. The stable roots give the hole box 0.935078, and no unsound bound is found. Residue: `part_dom`'s stored radii predate the fix (D3ct-5). |
| D3cp-3 | Resolved for the named survivors: the round's 30 mutants rerun as logged, and there are negative controls. Residue: new survivors (D3ct-4). |
| D3cp-4 | Largely resolved: the label is in (iii), Status and New, and max(0, λ₁) is stated on the axis. Residue: the clause "Λ⁺ … attained …, the first alias's lower rate holding where \|Re σ\| < 1", and "attained numerically" without supporting data (D3ct-3). |
| D3cp-5 | Resolved (49.89; Re(ηs) printed). Residue: 1.0052 should be floored to 1.0051 (D3ct-6). |
| D3cp-6 | Resolved: 2rb, the annulus, New's dangling clause, the branch-and-bound proviso. Residue: the hypothesis vs the proof (D3ct-1); "known only at those points" (D3ct-7). |
| D3cp-7 | Resolved: Chawla–Jain is cited [snippet]. Residue: it is Gauss–Legendre (D3ct-7). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ct-1 | minor | Ω′'s statement allows two points "on a line through η/2", so "distinct components" is what forces them onto opposite sides; this is now load-bearing, not redundant. The proof uses "its steepest-descent line", which is not assumed (only Re(φ″u²) < 0 is needed). "Within 2rb" must be strict. The chord condition, sampled at 481 points with \|t\| ≥ 0.05, follows from F − F₀ ≤ −\|φ″\|t²/2 + sup\|φ⁗\|t⁴/24 (φ‴(η/2) = 0, sup\|φ⁗\| ≤ 12/(\|η\|/2 − D)⁴), with ratio ≤ 6.4·10⁻³ on the annulus and ≤ 5.8·10⁻³ at the 57 points [own]. | State a = η/2 − du, b = η/2 + du, rb < d < 2rb, Re(φ″u²) < 0 (e.g. the descent line); distinctness then follows. Replace the chord sampling by the Taylor bound, which also justifies the 0.05 exclusion. D61 (iii)(b) and Proof; `path_check`; a test of the bound. |
| D3ct-2 | minor | A certificate is within reach, but the draft's recipe misses its decisive point. Monotone gradient-flow paths exist at 57/57 and 372/372 annulus points. With centred enclosures of F − F₀ the first-order estimate is δ = 2.4·10⁻³ … 1.6·10⁻², ≈ 1.5·10³ η-boxes, 35–50 segments per path, ≈ 6·10⁴ enclosures. Naive enclosures need δ ≈ 10⁻⁷ and ≈ 8·10⁹ boxes [own]. | Certify along the flow polylines with ball or interval arithmetic in centred (Taylor-model) form around the moving saddle, plus the Taylor chord; or at least record the estimate and "centred" in the Status. D61 Status, "What (b) still lacks"; part D. |
| D3ct-3 | minor | (iii) still says "Λ⁺ … is attained where \|Re σ\| > 1 (the edge saddles), the first alias's lower rate holding where \|Re σ\| < 1 (referee round 93)" without a label. It is wrong where λ₁ < 0 (the rate is max(0, λ₁), e.g. 0 at 1.05i), and λ₁ is undefined off the axis. At the sector's maximiser (Re ηs* = 1.0945) the local rate is 0.0158 … 0.0204 for K = 60 … 340, against Λ⁺ = 0.0233, and the segment has not overtaken s = ½ by K = 340 [own]. The md's "where Λ⁺ is attained" is not computed by `gap_arg`. | "Λ⁺ bounds the rate (attained, numerically, where \|Re σ\| > 1, as K → ∞: K ≤ 340 fits approach it); on the imaginary axis the rate is max(0, λ₁), a heuristic"; drop "attained" from the md's sentence, or compute it. D61 (iii); `d_lines`. |
| D3ct-4 | minor | 19 of 34 of my mutants survive. K25 (floor pixel: dpath 19, 19, 16), K27 (chord values 100×) and K32 (h pairing reversed: `part_dpath` raises) change published results. K29 (chord tolerance 10⁻³) and K31 (`path_ok` returns a float) are unsound but equivalent [own]. | Tests: pin the seed pixel's rounding (assert the used pixel lies outside the ball inside `path_check`); a chord bump of ~10·\|ms\|, not 0.1; `path_ok(…) is True`; `part_dpath` with two (R, h) pairs; pin the chord's maximum at \|t\| = 0.05, or test the Taylor bound. Log the 14 equivalents. tests/test_gh_peano.py. |
| D3ct-5 | cosmetic | The raw's D.dom.radius was not recomputed after `quad_roots`: 113 of 361 differ from the committed code by up to 4.36·10⁻¹², and the 5a171c7 code reproduces them exactly. The md is unaffected (6 digits). Only `dom_bb` was rerun [own]. | Rerun `part_dom` (with the branch and bound), store it, and say so in the log. raw; REVIEW_LOG. |
| D3ct-6 | cosmetic | "Re(ηs) between 1.0052 and 1.0999 on the sector [0.0861, 0.6175]": 1.0052 is the grid minimum (1.005156 at 0.617480) rounded to nearest. On the stated sector the infimum is 1.00514, and on the exact sector [0.08606, 0.61756] it is 1.00511. These numbers are computed at report time and are not stored in the raw [own]. | Floor (1.0051), or quote the exact sector. Store the sector in the raw. `d_lines`; D61 (iii), Numerics. |
| D3ct-7 | cosmetic | "Ω′ is known only at those points" and (iv)'s "the condition holds at the 57 points checked" overstate a numerical check. Chawla–Jain 1968 treat Gauss–Legendre: the Legendre function of the second kind, Barrett's ratio, steepest descents for entire f [snippet]. | "checked numerically only at"; "holds numerically"; "(for Gauss–Legendre)". D61 (iii)(b), (iv), Status. |

## Sources

Web searches only. Everything below is [snippet] unless marked otherwise.
- [Chawla–Jain, Asymptotic error estimates for the Gauss quadrature formula, Math. Comp. 22 (1968) 91–97](https://www.ams.org/mcom/1968-22-101/S0025-5718-1968-0223094-5/): for large n, according to the integrand's nature; the error via the Legendre function of the second kind and Barrett's ratio asymptotics; for entire f the contour is displaced freely and estimated by steepest descents through saddle points.
- [Chawla–Jain, Error estimates for Gauss quadrature formulas for analytic functions, Math. Comp. 22 (1968) 82–90](https://www.ams.org/journals/mcom/1968-22-101/S0025-5718-1968-0223093-3) [record only]. A BIT 1982 improvement of their Gauss–Legendre bound (Chen) [snippet, record only].
- [Convergence analysis of Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/pdf/2312.07940): contour representations of the Gauss–Hermite remainder (Donaldson–Elliott); a root-exponential bound.
- [Asymptotic error analysis of Gauss quadrature for nonsmooth functions (arXiv:2509.23532)](https://www.arxiv.org/abs/2509.23532) [record only].
- [Kazashi et al., Suboptimality of Gauss–Hermite quadrature (JNA 2023)](https://strathprints.strath.ac.uk/85830/7/Kazashi_etal_JNA_2023_Suboptimality_of_Gauss_Hermite_quadrature_and_optimality_of_the_trapezoidal.pdf): worst-case rates.
- [Computable upper error bounds for Krylov approximations to matrix exponentials (arXiv:1809.03369)](https://arxiv.org/pdf/1809.03369) and [Wang–Ye, Lanczos error bounds for the matrix exponential](https://www.ms.uky.edu/~qye/reports/M106373.pdf): non-asymptotic bounds.
- A search for certified (interval-arithmetic) saddle-point or contour-deformation bounds for Gauss–Hermite or Krylov errors found nothing relevant.
- [memory]: Olver, *Asymptotics and Special Functions* (Laplace's method for contour integrals); Saad 1992.

## Lesson

1. **Prove a step if it can be proved in two lines.** The chord's "below F(η/2) but at η/2" was sampled at 481 points with a hand-tuned exclusion. Taylor's theorem with φ‴(η/2) = 0 proves it with a margin of 99.4%.
2. **The cost of a certificate is set by how a difference is enclosed.** Near a moving saddle, F(σ; η) − F₀(η) varies like (σ − η/2)·φ″. Enclosed as one centred form, the annulus needs about 10³ boxes; enclosed term by term, about 10¹⁰.
3. **Rerun everything a changed function feeds.** `quad_roots` changed `gap_closed`, which feeds `dom_radius` as well as `dom_bb`. Only the latter was rerun, so the raw no longer matches the code bit for bit.
4. **A negative control must sit near the threshold.** A 0.1 bump cannot tell a strict test from one with a 10⁻³ tolerance, and a truthiness assert cannot tell True from a margin.

All my files are under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref99:
- own/: rerun.py, rawdiff.py, radius_diff.py, radius_old.py, numall.py, dom_own.py, gapfun.py, sector_chk.py, valley_own.py, feas.py, lam_att.py;
- mut/: spec99.py, recompute.py, rc/ and w99/, w95o/;
- logs/: rerun_bb.log, rerun_dpath.log, dom_own.log, taylor.log, flow.log, grid_var.log, feas2.log, lam_att.log, lam_att2.log, mut99_a.log, mut99_b.log, rc_a–d.log, mut95o.log.
