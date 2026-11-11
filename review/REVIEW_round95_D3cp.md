# Referee report D3cp, round 95: Proposition D‴.61 as revised after round 93

## Scope

I reviewed commit 5a171c7 on branch claude/vigilant-clarke-hf939b.

- **Setup**
  - **Work directory.** The commit was exported with `git archive 5a171c7 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref95.
  - **Outputs.** GP_OUT = W/gp_out and WM_OUT = W/wm_out each hold a copy of their raw.
  - **Run conditions.** Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. Pytest ran with `--basetemp` W/pyt.
  - **Interpreters.** Repository code ran in the repository's venv; for tools/mutants.py I symlinked it as W/repo/venv. My own code is in W/own and is all new; it ran in …/ref90/venv2 where it needed mpmath.
- **Heavy processes.** I kept to at most two of mine at a time, with two exceptions.
  - A third (my own dom_check.py) started by mistake. I stopped it after checking /proc/PID/cmdline, then reran it later.
  - My own mutant H11 made its pytest grow without bound (1.7 GB after 114 s). I stopped that pytest after checking its cmdline and its cwd, the mutant slot.
  - I touched no other process.
- **Read in full**
  - review/drafts/D61.md, experiments/exp_gh_peano.py, results/chaos/gh_peano.md and tests/test_gh_peano.py.
  - review/REVIEW_round93_D3cn.md and the Round 93 entry of review/REVIEW_LOG.md.
  - tools/mutants/round93_D61_own.py and tools/mutants.py.
  - `git diff 363855d 5a171c7` for these files.
- **Repository untouched.** My only git commands were archive, diff, log, status, rev-parse and show.
  - During the review HEAD moved to bc0c89b (D‴.58/59/60 commits, not mine). No D‴.61 file changed.
  - Commit 420a3ef updated D58.md's cross-references to D‴.61.
- **Web.** I used search snippets only. www.ams.org does not resolve from this session (ENOTFOUND).

## Summary

**Design.** The revision applies round 93.
- It corrects the alias rate to (5/8)y² − log y − ½ − (3/2)log 2, with Prony fits (`qxp_job`, `prony`).
- It puts the domination condition in closed form (`gap_closed`) and verifies it on |η| ≤ 1.8927 with a floating-point branch and bound (`box_bound`, `dom_bb`).
- It sets r* = 0.8052 and σ*/σ < 9.26, both rounded down.
- It adds a new argument. The u-integrand is entire, so Duhamel's path can be deformed: 0 → a → the segment through η/2 inside B(η/2, ½(1 − |η|/2)) → b → η, with a and b in the components of {F_η < F_η(η/2)} that contain 0 and η.
  - This proves the limit on an open set Ω′ ⊇ Ω.
  - Membership of Ω′ is grid-checked (`path_ok`) at 57 points.
  - (iv) is stated for r < 0.9607 (σ*/σ < 49.92), conditional on the valley condition holding on |η| ≤ 1.98.

**Verdicts**
- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte, and all 26 tests pass.
  - All 27 comparisons I made are bit-identical, including every new part-D job.
  - Every decimal in the draft is in the md verbatim.
- **The deformed-path argument is mathematically correct, as a theorem conditional on membership in Ω′.**
  - The bound F is a valid asymptotic upper bound everywhere the path goes, including near and on |σ| = 1 by the band argument, and at the endpoints.
  - The pieces outside the ball are exponentially smaller by compactness.
  - Laplace's method on the ball segment is legitimate.
- **The 57-point valley check is a numerical step, as labelled.** It is stable under grid refinement and connectivity changes, but carries no certified margin. Near |η| = 1.98 the seeds' margins are about 10⁻⁶ (D3cp-1).
- **The hypotheses of the conditional consequences are right.**
  - r_star(1.98) = 0.9607245, which gives r < 0.9607.
  - 49.92 is the floor of the unrounded value 49.9224. The floored r gives 49.89 (D3cp-5).
- **The closed form and its candidate set are exact, and `box_bound` is a valid bound in exact arithmetic.**
  - `dom_bb`'s 402,664 leaves all recompute to 8·10⁻¹⁶ at 40 digits.
  - But `bound_sup`'s quadratic formula drops a root when the leading coefficient is cos(π/2) ≈ 6·10⁻¹⁷. It then returns values up to 0.09 below the true supremum on boxes whose upper angle is π/4. The data are not affected (D3cp-2).
- **The alias rate is right.** My own exact c_K, Prony-fitted at K = 300–312, gives λ₁ to 6–9·10⁻⁶. But the continuation to y = 1.8589 and the alias rate are now stated without their heuristic label (D3cp-4).
- **Mutants**
  - Of my 32 mutants of the new code, 15 are killed and one diverges.
  - 16 survive; 9 of them change published results on recompute (D3cp-3).
  - Round 93's own spec reruns as logged: 12 killed, 6 equivalent.
- **Round-93 issues:** all 8 resolved, with residues noted below. No major issue.

**Findings**
- **D3cp-1 (minor). The valley check has no certified margin, and its seeds sit where the margin is tiny.**
  - The seeds a, b = η/2 ∓ 1.2·rb·η/|η| lie on the ray. Along the ray the margin is only F − F0 ≈ −(4 − Re η²)(1.2rb/|η|)², because near |η| = 2 that direction lies at the edge of the saddle's descent sector.
  - At |η| = 1.98 this is −7·10⁻⁷ … −5·10⁻⁵, comparable to F's change across one pixel (≈ 8·10⁻⁶ at h = 0.001) [own].
  - The verdicts are stable under h → h/2, 4- vs 8-connectivity, and seeds on the steepest-descent line. A dense scan of 222 points on |η| = 1.90 … 1.98 passes [own].
  - The grid breaks down beyond the claimed radii. At |η| ≥ 1.99 the rounded seed pixel can fall inside the ball, giving false negatives even at points of Ω, for example 1.995e^{0.785i}.
- **D3cp-2 (minor). A hole in `bound_sup` (where `dom_bb` can fail).**
  - The test `abs(qa) < 1e-300` misses qa = R²cos(π/2) ≈ 6·10⁻¹⁷. The root formula then cancels to w = 0, and the true critical point w ≈ 1/(2B) is lost.
  - Example: box [1.5, 1.8927] × [0, π/4] gives 0.8458, against the true supremum 0.9351.
  - No verified leaf is affected [own]. But the routine is not "sound modulo rounding" as written.
- **D3cp-3 (minor). The tests miss result-changing mutants of the new code.**
  - H02 (wrong quadratic) changes the violated sector at 1.999 to [0.004, 0.720].
  - H06 (−log R₂, an unsound bound) still verifies, with 201,762 boxes.
  - H10 (only [0, π/4] covered) changes only the box counts, while the md still says π/2.
  - H12 and H13 change the B&B's worst bound and box counts.
  - H15 and H21 change the dpath counts.
  - H23 (angles [0, π/4]) changes the raw while the md still says [0, π/2].
  - H26 adds y = 1.8 to the Prony fits.
  - H18 and H22 (C = 0 beyond the disc, i.e. no bound at all) are equivalent on the data, so the valley check never exercises Λ⁺. There is no negative control.
- **D3cp-4 (minor). Unqualified claims about the corner ratio beyond the disc.**
  - The draft says: "c_K(i√Ky) → 1 + y² beyond the disc up to y = 1.8589, and beyond it … grows at the first alias's rate". This appears in (iii), and in New as a new result.
  - It is not proved. (a) covers |μ| < 1 only. The alias is the saddle of the first term of an expansion of 1/I_K that does not converge on all of |ζ| = 1/√2.
  - This round dropped "as a heuristic" from New.
  - "The first alias's lower rate holding where |Re σ| < 1" misstates the rate, which on the axis is max(0, λ₁).
  - "Does overtake" is not observed at K ≤ 160 [own].
- **D3cp-5 (cosmetic). Rounding.**
  - "r < 0.9607 (σ*/σ < 49.92)" pairs a floored r with an unfloored ratio; the floored r gives 49.89.
  - The maximiser's Re(ηs*) on the stated sector is 1.0067–1.0999, not "1.009 and 1.099".
- **D3cp-6 (cosmetic). The wording of Ω′ and of New.**
  - The two points "beside the ball" are unspecified. They must lie within 2rb of η/2 for the segment a–b to lie in the unit disc, as the proof says.
  - The clause "those differ" is implied by the other conditions.
  - The hypothesis of (iv) really concerns the annulus 1.8927 < |η| ≤ 1.98.
  - New has a dangling "and those rules' sharp asymptotics", and its "proved … on |η| < 1.8927" lacks "modulo the branch and bound".
- **D3cp-7 (cosmetic). Literature.** Chawla–Jain (1968) deform the Gauss error contour through saddle points by steepest descents for entire integrands [snippet]. That is the classical precedent of the deformed path, and it is not cited.

## What I checked and how

**Reproducibility** [own]
- `--from-raw` gives an identical md (sha256 9a45399b…).
- tests/test_gh_peano.py: 26 passed in 2.0 s.
- With the script's own functions (W/own/rerun.py; `repr` equality, `secs` dropped) these are all bit-identical against the raw:
  - part D: `part_dpath`, `part_dom`, `dom_bb(1.78)` (4935 boxes) and `dom_bb(1.8927)` (402,664 boxes);
  - `qxp_job(100)` and `qxp_job(150)` (26 K × 4 y);
  - all 5 `duhamel_job`s, `qx_job(20)` and `qx_job(160)`, and `d_job` at (160, 0), (160, 1), (160, 2), (80, 4);
  - parts P and M; L (160, real), LC (160, π/8), Q (160, 0), S[10] (τ/σ = 2, v′ = 2.5), R[0], R[13], R[17] and R's 14 rules without a frame.

**The claims against the md and the raw** [md][raw][own]
- W/own/numcheck.py finds 211 distinct decimals in the draft. All are in the md except "Theorem 2.10" and three arXiv ids.
- I spot-checked the semantics of these (raw `D['dom']`, `D['dpath']`, `D['qxp']`):
  - 0.9607 and 49.92, 0.0966, 1.892735 at 0.345;
  - 402664 and 4935, and −1.9·10⁻¹¹;
  - the margins' angles (0.000, i.e. the real axis);
  - the 19/19 counts;
  - the Prony rates and frequencies.
- Title and statement: no overreach, except D3cp-4.

**The deformed path, re-derived** [own]
- **The bound F.**
  - The exponent per K of Duhamel's u-integrand is log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + (1/K)log|c_K(√Kσ)|.
  - Inside the disc, (a) gives (1/K)log|c_K| → 0 uniformly on |σ| ≤ r < 1.
  - On the band r < |σ| < 1 the ρ₀-contour gives Λ⁺ = |σ|²/2 − ½ − log|σ| ≥ 0. Λ⁺ vanishes on |σ| = 1, where C is continuous, so a path with margin ε survives once the band is thin.
  - Outside the disc the ρ₀-contour bound holds for every σ, uniformly on compact sets: the ellipse encloses the nodes and 1/I_K is bounded on |ζ| = ρ₀. Λ_{ρ₀} → Λ⁺ uniformly as ρ₀ → 1.
  - At the endpoints, F → −∞ (log|σ|, log|η − σ|), and the integrand is bounded there.
  - So F is a valid limsup bound wherever the path goes.
- **The ball.** It lies in |σ| ≤ ½ + |η|/4 < 1. η/2 is a saddle of φ(σ) = log(η − σ) + (η − σ)²/2 + log σ + σ²/2, with φ″ = 2 − 8/η². Along the ray, Re(φ″u²) = −2(4 − Re η²)/|η|² < 0. Laplace's method with complex phase applies, with c_K → 1 − σ² uniformly on the ball.
- **The pieces.** They are compact and lie in the open sublevel set, so they have a positive margin. The condition is open in η. Ω ⊂ Ω′, because on Ω the straight segment itself lies in the sublevel set.
- **The separation clause** ("those differ") is implied. Otherwise a path avoiding the saddle would make the integral exponentially smaller than its nonzero main term.

**The valley check, numerically** [own]
- **Code.** W/own/valley.py, scan.py and focus.py: my own F, bottleneck levels by bisection on `ndimage.label`, and variants.
- **Margins at the 57 points.** These are F(seed) − F0 and the bottleneck levels ta, tb (seed to 0 or η, ball excluded):

  | \|η\| | seeds' F − F0 | separation level tsep | straight segment's max where it fails |
  |---|---|---|---|
  | 1.9 | −9.8·10⁻⁵ … −1.9·10⁻³ | +8·10⁻⁵ … +1.6·10⁻³ | — |
  | 1.95 | −1.2·10⁻⁵ … −4.6·10⁻⁴ | +1·10⁻⁵ … +4·10⁻⁴ | +3.8·10⁻⁶ … +7.8·10⁻³ |
  | 1.98 | −7.3·10⁻⁷ … −5·10⁻⁵ | +5·10⁻⁷ … +5·10⁻⁵ | +2.1·10⁻⁴ … +1.4·10⁻² |

  - At the steepest-descent points 0.05 from η/2 the margins are −7·10⁻⁵ … −5·10⁻³. The bottleneck from there to η equals the point's own value, so the detour descends from there.
- **Dense scan.** 222 points on |η| = 1.90, 1.92, 1.94, 1.96, 1.97, 1.98 at 37 angles each.
  - `path_ok` holds at all of them, as do the 4-connected variant, steepest-descent seeds, and h/2 for R < 1.98.
  - At 1.98 with h = 0.0005 (pad 0.3), the verdicts and bottlenecks agree with h = 0.001: tb = −1.55·10⁻⁶ at both h for θ = 0.087.
- **Beyond the claim (|η| = 1.99, 1.995).** The verdicts flip with the grid's alignment. Every failure has a rounded seed pixel inside the ball (distance/rb = 0.75–0.97), so these are false negatives, including at 1.995e^{0.785i} ∈ Ω.
  - `path_ok` therefore needs h/√2 < 0.2·rb. The choices for 1.9, 1.95 and 1.98 satisfy it, the last only barely.
- **No false positive found.** Possible sources remain uncertified: 8-connectivity, the 10⁻¹² tolerance, and sampling between pixels.
- **What rigour needs.** Seeds on the steepest-descent line. The local piece by an explicit Taylor remainder of φ, with Λ⁺ ≲ (|σ| − 1)² near the circle. The rest along explicit polylines, with F bounded by interval arithmetic (it is elementary). η-boxes covering the annulus 1.8927 < |η| ≤ 1.98, as `dom_bb` does for Ω.

**Finite K along a detour** [own] (W/own/finiteK.py: exact c_K from my integer series)
- At η = 1.95e^{0.35i}, where the straight segment fails: a minimax path from b to η, with bottleneck −7.6·10⁻⁵ below F0.
- At 14 points along it, the exact (1/K)log|integrand| minus its value at η/2 is negative at K = 40, 80 and 160. At K = 160 it runs −0.0000, −0.0025, −0.0107 … −2.06.
- The bound F − F0 at the same points is −0.0001 … −1.99. The exact values stay below it and rise toward it as K grows.
- The straight segment's worst point, s = 0.600 (Re ηs = 1.099, bound +0.0059), gives −0.0230, −0.0135 and −0.0076 at K = 40, 80, 160. It has not yet overtaken s = ½.

**The domination condition and the branch and bound** [own] (W/own/dom_check.py, bb_edge.py)
- **The closed form.** I checked the candidate set analytically: the endpoint 1 − 1/R, and the roots of (p + 2B²)w² − 2B²w + 1 and pw² − 2Bw + 1. M is C¹ at 1, so the kink is harmless. Extra candidates from the wrong piece are harmless too.
  - Against a direct 200,001-point supremum at 300 random η, the closed form is never below it; the grid supremum exceeds it by at most 8.9·10⁻¹⁶.
- **The box bound.** −log R₁, p_lo at a corner (monotone at fixed sign) and B_hi = R₂cos φ₁ make it a valid bound, using w² < ½ (since wm < 0.4716) and M increasing.
- **The leaves.** My own reimplementation collects the same 4935 and 402,664 leaves. Recomputed at 40 digits with stable roots:
  - max |float − mp| is 5.0·10⁻¹⁶ and 8.1·10⁻¹⁶;
  - the worst bound is −3.2813·10⁻⁷ and −1.9432·10⁻¹¹;
  - a dense 20,001-point check on 3000 leaves, plus 9 lattice points each, never exceeds the box bound (excess ≤ 4.4·10⁻¹⁶).
- **The hole.** The smallest |leading coefficient| among leaves is 6.1·10⁻¹⁷, at φ = π/4. On boxes ending at π/4 with B > 1.5, `bound_sup` is too low:

  | box | float | true |
  |---|---|---|
  | [1.5, 1.8927] × [0, π/4] | 0.8458 | 0.9351 |
  | [1.7, 1.8927] × [0.5, π/4] | 0.4513 | 0.4771 |
  | [1.88, 1.8927] × [π/4 − π/16, π/4] | 0.2485 | 0.2559 |

  - Boxes from π/4 − π/32 on are unaffected (B < 1.5).
- **Radii and constants.** The least radius is 1.892732716 at arg 0.34320 (2001 angles plus bisection).
  - r_star(1.8927) = 0.8052030, which gives 0.8052; (1 + r)/(1 − r) = 9.2671, which gives 9.26.
  - r_star(1.98) = 0.9607245, with 49.9224 unfloored and 49.8906 from the floored r.
  - r_star(1.0) = 0.0966496.
  - The violated sector at 1.95 is [0.0861, 0.6175]. The maximiser's Re(ηs*) is 1.0067–1.0999 on [0.087, 0.615] and 1.0051–1.0999 on the exact sector.

**The alias rate** [own] (W/own/misc.py, alias_check.py)
- **The saddle.** At ζ = (iy ± √(8 − y²))/4: |ζ|² = ½ exactly, |f₁′| ≤ 1.3·10⁻¹⁵, and Re f₁ − Re f_σ(σ) = λ₁ to 12 digits.
  - Im mod 2π is 2.13177, 2.28540, 2.56931, 2.92170, the md's Prony frequencies.
  - The zero is at 1.8589069607.
- **Exact c_K.** I computed c_K(u) = (K − 1)!e^{−u²/2}Σ_m u^m h_m/(m + K − 1)!, with h_m the nodes' complete homogeneous polynomials from He_K's integer coefficients. It matches the raw's `qxp` values to 12 digits at K = 100–102.
- **Prony fits.** My own fits on the raw's rows reproduce the md's rates; dropping two K changes them by ≤ 10⁻⁵. At K = 300–312:

  | y | Prony rate | λ₁ | frequency |
  |---|---|---|---|
  | 1.9 | 0.0746838 | 0.0746753 | 2.13177 |
  | 2.0 | 0.2671395 | 0.2671320 | 2.28540 |
  | 2.2 | 0.6968281 | 0.6968219 | 2.56931 |

  - At y = 1.9 the error shrinks from 7·10⁻⁵ (K₀ = 100) to 3·10⁻⁵ (K₀ = 150) to 9·10⁻⁶ (K₀ = 300), about as K⁻².

**Mutants** [own]
- **Round 93's own spec,** rerun with tools/mutants.py: 12 of 18 killed. The survivors are G02, G06, G10, G11, G12, G17, as the log says. The control is clean.
- **My spec,** W/mut/spec95.py: 32 single-point mutants of round 93's new code.
  - Killed by the tests (15): H01, H03, H04, H05, H07, H08, H09, H17, H19, H24, H25, H27, H30, H31, H32.
  - Divergent (1): H11, which splits the shorter side, so its stack grows without bound. I stopped it.
  - **Survivors, recomputed** with W/mut/recompute.py (the part touched, then the md rebuilt). These 9 change published results:

    | mutant | change | effect on recompute |
    |---|---|---|
    | H02 | first quadratic with B² for 2B² | violated sector at 1.999 becomes [0.004, 0.720] |
    | H06 | −log R₂ in the box bound | B&B verifies with 1701 / 201,762 boxes |
    | H10 | quadrant halved to [0, π/4] | 4594 / 399,341 boxes; md still says π/2 |
    | H12 | worst = min | −inf |
    | H13 | slack 10⁻⁹ in the code | 402,724 boxes, worst −1.0·10⁻⁹ |
    | H15 | tolerance F < F0 + 10⁻⁴ | dpath counts 18, 15, 0 |
    | H21 | seeds at 3rb | dpath counts 19, 18, 17 |
    | H23 | angles on [0, π/4] | raw changes; md unchanged but now wrong |
    | H26 | y with λ₁ > −0.2 | y = 1.8 added to the Prony fits |

  - **Equivalent on the data:** H14, H16, H20, H28, H29, and two unsound ones, H18 (C = 0 up to |σ| = 1.05) and H22 (C = 0 everywhere).

**Literature** (web search)
- Chawla–Jain, *Asymptotic error estimates for the Gauss quadrature formula*, Math. Comp. 22 (1968) 91–97, §3.1 [snippet]. For an entire integrand they displace the error contour freely, avoiding the branch points, and estimate it by steepest descents through the saddle points.
- The companion paper (Math. Comp. 22, 82–90) [snippet: record only].
- Hanke: precise asymptotics of the Lanczos error via orthogonal polynomials [snippet].
- Ye (arXiv:1603.07358): Faber-polynomial bounds [snippet].
- arXiv:2509.23532: contour-integral asymptotics of Gauss errors [snippet].
- Pozza–Pranić–Strakoš (ETNA 50): Gauss quadrature and complex Jacobi matrices [snippet].
- I found no sharp asymptotics of the Gauss–Hermite or Lanczos error of e^{ξz} at ξ ~ √K, and no deformed-path argument on the Krylov (Duhamel) representation. Saad 1992 [memory].

## Verdicts per part

- **(i): proved** (unchanged).
- **(ii): proved.** P_j's polynomial form remains a sound sketch.
- **(iii)(a): proved.** The new wording is correct: the ellipse for every K, I_K zero-free, c_K > 0 from the Dirichlet form.
- **(iii)(b)**
  - On |η| < 1 and on both axes: proved.
  - On |η| ≤ 1.8927: proved modulo the floating-point B&B. That is sound on the data, confirmed at 40 digits for every leaf. The code has the hole of D3cp-2.
  - On Ω′: the deformed-path theorem is correct.
  - The 57 points' membership is numerical. It is robust to refinement but uncertified (D3cp-1).
  - The sector stays conjectural elsewhere near |η| = 2.
- **(iii)(c): proved.**
- **The alias rate:** the formula is right and is numerically confirmed to K = 312. It is not proved, and must be labelled so (D3cp-4).
- **(iv): correct.**
  - r* = 0.8052 and 9.26 are valid floors.
  - The conditional r < 0.9607 has the right hypothesis: the valley condition on the closed disc |η| ≤ 1.98, compact in Ω′, which in effect means on 1.8927 < |η| ≤ 1.98.
  - 49.92 is a valid floor but rounded inconsistently (D3cp-5).
- **(v): proved** (unchanged).

## Status of the round-93 issues

| id | status |
|---|---|
| D3cn-1 | Resolved. The rate, the zero 1.8589, the Prony fits and the test pin are all in; my exact c_K confirms them at K = 300. Residue: the claim is now unlabelled (D3cp-4). |
| D3cn-2 | Resolved. Residues: the maximiser range is 1.0067–1.0999, and the "lower rate" clause is imprecise (D3cp-4, D3cp-5). |
| D3cn-3 | Largely resolved: the closed form and the B&B. Interval arithmetic was not adopted (acknowledged). Residues: D3cp-2 and D3cp-3. |
| D3cn-4 | Resolved in D61. At 5a171c7, D58.md still said 0.8053, 9.27 and "grid-checked"; commit 420a3ef (outside the range) fixes this. |
| D3cn-5 | Resolved. |
| D3cn-6 | Resolved for r*, 9.26 and "at least". Residue: 49.92 (D3cp-5). |
| D3cn-7 | Resolved. |
| D3cn-8 | Resolved: qx's `prec`, DOM_TH, DOM_NR and DOM_MARGIN are pinned. |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cp-1 | minor | The valley check is uncertified and its ray seeds have margins of about (4 − Re η²)(1.2rb/\|η\|)²: −7·10⁻⁷ … −5·10⁻⁵ at 1.98, comparable to F's change per pixel. The 1e-12 tolerance, 8-connectivity and pixel sampling certify nothing. The verdicts are stable (h/2, 4-connectivity, steepest-descent seeds; 222 annulus points pass), but at \|η\| ≥ 1.99 seed pixels fall in the ball, giving false negatives even in Ω [own]. | Seeds on the steepest-descent line. The local piece by a Taylor remainder; the rest along explicit polylines with interval bounds on F; η-boxes over 1.8927 < \|η\| ≤ 1.98. Meanwhile assert h/√2 < 0.2rb in `path_ok`, and report the margins and the refinement in part D. `path_ok`, `part_dpath`; D61 (iii)(b), Status. |
| D3cp-2 | minor | `bound_sup` loses the critical point w ≈ 1/(2B) when qa = R²cos(π/2) ≈ 6·10⁻¹⁷, because the 1e-300 test misses it and the root cancels to 0. Box [1.5, 1.8927] × [0, π/4] gives 0.8458 against the true 0.9351. No leaf is affected: 402,664 leaves agree with 40-digit values to 8.1·10⁻¹⁶ [own]. | Stable roots (q = −(b + sgn b·√disc)/2; q/a and c/q) or a relative threshold. Add a test on a box ending at π/4. Optionally record a high-precision leaf recheck. `bound_sup`; tests. |
| D3cp-3 | minor | 16 of 32 new-code mutants survive. 9 change published results on recompute: H02, H06, H10, H12, H13, H15, H21, H23 (the md would misreport [0, π/2]) and H26. H18 and H22 (no bound beyond the disc) are equivalent on the data: there is no negative control [own]. | Tests: `gap_closed` against a dense supremum at random η; `bound_sup`'s R₁ against R₂; `dom_bb`'s coverage (area, or a planted violation at arg η > π/4); `worst` as the maximum; a negative control for `path_ok` (the segment's violated point excluded at 1.95e^{0.35i}); pins for part_dpath's range, the seed factor 1.2 and qxp's y on the real QX_Y. Log the equivalents. tests/test_gh_peano.py. |
| D3cp-4 | minor | (iii) and New state the continuation to 1.8589 and the alias rate as results; New dropped "as a heuristic". (a) covers \|μ\| < 1 only, and the alias is the saddle of a non-uniformly convergent expansion of 1/I_K. On the axis the rate is max(0, λ₁), e.g. 0 at 1.05i, where λ₁ < 0 < Λ⁺. Off the axis λ₁ is undefined. The overtaking is not seen at K ≤ 160 (−0.0076) [own]. | "Numerically (K ≤ 312, Prony), with a saddle-point heuristic"; restore the label in New; give the rate as max(0, λ₁); "overtakes asymptotically, Λ⁺ being attained there numerically". D61 (iii), Status, New. |
| D3cp-5 | cosmetic | "r < 0.9607 (σ*/σ < 49.92)": the floored r gives 49.89. "Re(ηs) between 1.009 and 1.099": it is 1.0067–1.0999 (1.0051 on the exact sector) [own]. | Floor consistently in `d_lines`' dpath line. 1.007 … 1.100. D61 (iii), (iv), Status, New; md. |
| D3cp-6 | cosmetic | The points "beside the ball" are unspecified (they need to be within 2rb). The clause "those differ" is implied. The "where it is unchecked" scope is loose. (iv)'s hypothesis is in effect the annulus. New has a dangling "and those rules' sharp asymptotics", and its "proved … on \|η\| < 1.8927" has no B&B proviso [md]. | One phrase each. D61 (iii)(b), (iv), New. |
| D3cp-7 | cosmetic | Chawla–Jain 1968 already deform the Gauss error contour through saddle points for entire integrands [snippet]. | Cite it as the precedent; the novelty is applying it to the Krylov (Duhamel) representation. D61 Status. |

## Sources

Web searches only. Everything below is [snippet] unless marked otherwise.
- [Chawla–Jain, Asymptotic error estimates for the Gauss quadrature formula, Math. Comp. 22 (1968) 91–97](https://www.ams.org/mcom/1968-22-101/S0025-5718-1968-0223094-5/): §3.1, steepest descents on the displaced error contour.
- [Chawla–Jain, Error estimates for Gauss quadrature formulas for analytic functions, Math. Comp. 22 (1968) 82–90](https://www.ams.org/journals/mcom/1968-22-101/S0025-5718-1968-0223093-3) [record only].
- [Hanke, Superlinear convergence rates for the Lanczos method applied to elliptic operators](https://kluedo.ub.rptu.de/files/4864/Hanke_Superlinear+convergence+rates+for+the+Lanczos+method+applied+to+elliptic+operators.pdf)
- [Ye, error bounds for Krylov approximation of matrix exponentials (arXiv:1603.07358)](https://arxiv.org/pdf/1603.07358)
- [Asymptotic error analysis of Gauss quadrature for nonsmooth functions (arXiv:2509.23532)](https://arxiv.org/abs/2509.23532)
- [Pozza–Pranić–Strakoš, ETNA 50 (2018)](https://etna.ricam.oeaw.ac.at/vol.50.2018/contents.pdf)
- [Convergence analysis of Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/pdf/2312.07940)
- [Method of steepest descent (Wikipedia)](https://en.wikipedia.org/wiki/Method_of_steepest_descent)

From memory [memory], not seen here: Saad, SIAM J. Numer. Anal. 29 (1992); Olver, *Asymptotics and Special Functions* (Laplace's method with complex phase).

## Lesson

1. **A grid check is only as good as its seeds.** On the ray, near the edge of the saddle's descent sector, the margin is O((4 − Re η²)rb²), about 10⁻⁶, which no pixel grid certifies. Seed on the steepest-descent line, prove the saddle's neighbourhood by hand, and certify the rest along explicit paths with interval bounds.
2. **A floating-point branch and bound is only as sound as its formulas.** An absolute "zero" test of 10⁻³⁰⁰ misses cos(π/2) = 6·10⁻¹⁷ and silently drops a root. Use stable roots, and recheck the leaves at high precision.
3. **A test suite needs negative controls.** A check that passes with no bound at all (C ≡ 0) is not testing the bound.

All my files are under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref95:
- own/: rerun.py, numcheck.py, dom_check.py, bb_edge.py, valley.py, tab_valley.py, scan.py, focus.py, misc.py, alias_check.py, finiteK.py;
- mut/: spec95.py, recompute.py, and the slots w95 and w93o;
- logs/: rerun_a–d.log, dom_check.log, valley_pts.log, scan.log, focus.log, mut95.log, mut93o.log, recompute.log, alias_finiteK2.log.
