## Referee D3at, round 47: Proposition D‴.54, revised after round 42

**Summary.** Most of the round-42 fixes hold in the plane. These are x_θ(t) and the face averages, the edge factors posed at t = H/2, the power diagram, the aliasing and window formulas, the hexagonal check and the softmax error stated per share. Part (v), new in this revision, fails on every count I tested:
- Its balance is inconsistent. The entrant's shift that makes the binding quartic also cancels the linear pull. The pull that survives is a curvature term, so the stated unbinding level is 1–2 orders of magnitude too low.
- Runner-ups that *form* are captured too, not only tails below θ. (iv)'s count therefore fails near a same-sign unfitted feature.
- "The count per feature by final positions does not differ" is false. A capture that arrives as the feature's own entrant enters forms two kernels at that feature. This happens on C2 itself and on D‴.51's own 1-D target.
- (v) is D‴.51(ii)(c)'s "migrating loser", and it is not cited.

The n-D product rule repeats round 42's wrong-geometry error for n ≥ 3. B4's 213 unsettled inner problems do not affect the areas, but they ran out of doublings, not wall-clock time. The tests kill 4 of 24 mutants.

### What I checked and how
- **Setup.**
  - Scratch copy at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review47/repo` (venv and .git excluded); the original venv's python throughout.
  - Scripts and logs in `.../review47/work` and `.../review47/work/v`; mutation sandbox in `.../review47/mut`.
  - Every run at nice 15 with OMP_NUM_THREADS=1, one heavy process at a time.
  - Nothing in the repository was touched; git was used read-only. I killed only my own PIDs (7816, 8449).
- **Read:** D54.md, the round-42 report and its log entry, `rich_transit_nd.py`, `exp_transit_nd.py`, the md and raw, and the tests. Context: D‴.51 and D‴.53 in PARADIGM.md and transit_mass.md.
- **Tests and regeneration.**
  - 9 of 9 tests pass in 6.4 s.
  - The md regenerates byte-identically from the raw (`--from-raw`).
  - `c2_job` rerun at ε = 10⁻⁹ and 10⁻¹⁰ reproduces the raw exactly.
  - Raw provenance: A's loops, D's loops and all of Part B (loops and inner problems) are bit-identical to round 42's raw (72c4fe5). They were not rerun; B's inner column comes from B4.
- **Numerics paragraph:** every number checked against the md and raw.
- **New computations:**
  - V1, a static test of (v)'s binding in 1-D at a fixed tail coefficient, for d = 4–7σ.
  - C2's near tie at 20 further masses, from 1.5·10⁻¹⁰ to 7·10⁻⁹. Key points rechecked at rtol 10⁻¹¹, 10⁻¹² and 10⁻¹³.
  - Two C2 variants with the unfitted feature at 5σ, at −60° and at 160°.
  - D‴.51's H(3) target in 1-D (c₂ = 0.043 and 0.048) over 30 masses and at rtol 10⁻⁹ to 10⁻¹³.
  - 105 of B4's 213 unsettled problems rerun.
  - Part D's rule with offset-blind thresholds.
  - Product-rule areas on grids of step 1 down to 0.02.
  - 24 mutants, plus the md regenerated under 6 surviving mutants.

### Verdicts
- **(i) Dimension-free.** The creep's vector law holds (C). The claim still includes D‴.51(iv) and D‴.53's ε^{1/2} rate and first-order correction in ℝⁿ. Neither is implemented or tested in n-D, and C measures only the end state.
- **(ii) Selection.**
  - (a) and (b) are correct; W2 checks the entrant.
  - (c)'s aliasing and window formulas are correct, and the W1 zone is reproduced. Cross-talk is stated as a condition but never tested.
- **(iii) The k-fold near tie.**
  - 𝒮(x; y) and x_θ(t) (T, T2) are right at A = 0.8.
  - (a) is right. (b), the profile split at the centre and its departure off-centre (W4), is right; the per-share softmax error (W5) is right.
  - (c) is fine.
  - (d): the multiplicativity is exact, and the plane's near-factorisation is verified (9.1·10⁻⁴ over all shares, 2.0·10⁻⁴ for the edges' runner-ups). But "the offset h/2" is wrong for n ≥ 3 (D3at-6).
  - Everything in (iii) is for equal masses only; this is unstated (D3at-7).
- **(iv) Formation rule.**
  - E[N − 1] with the face averages is correct, and 71.5, 105.8, 139.0, 171.3 ε^{1/2} reproduce.
  - The codimension ranges are stated correctly.
  - B4's areas reproduce: 33.5 and 15.5, so P(N ≥ 3) = 1112ε and P(N = 4) = 513ε.
  - The unsettled problems are irrelevant to the areas. The 0.25 grid carries about 1% error on N = 4: the product rule gives 15.08 there against 15.27 at step 0.02.
  - The rule fails in multi-feature targets (D3at-2), and Part D cannot distinguish x_θ(t) from x_θ(0) (D3at-11).
- **(v) The tail's capture.**
  - Right: the quartic binding and the speed ∝ w.
  - Wrong: the pull's order, the level, the scope and the count claim.
  - The C2 window exists but its upper edge sits 35–56× above the stated level.
  - Not new (D‴.51(ii)(c)).

**Round-42 issues**

| issue | status | note |
|---|---|---|
| D3ao-1 | fixed in the plane, partly in n-D | D3at-6, D3at-11 |
| D3ao-2 | fixed in 2-D, recurs for n ≥ 3 | D3at-6 |
| D3ao-3 | fixed for (ii) | unequal-mass split not addressed (D3at-7) |
| D3ao-4 | fixed | the generic-samples split is checked only by round 42 |
| D3ao-5 | partly | D3at-10 |
| D3ao-6 | fixed | |
| D3ao-7 | fixed | D‴.51's wording now says "at the feature's centre" |
| D3ao-8 | fixed | |
| D3ao-9 | fixed | |
| D3ao-10 | partly | D3at-8, D3at-9 |
| D3ao-11 | fixed | |
| D3ao-12 | partly | D3at-11 |
| D3ao-13 | fixed | labelled |
| D3ao-14 | mostly | new mismatches in D3at-13 |
| D3ao-15 | partly | 20 of 24 mutants survive |
| D3ao-16 | fixed for (ii)–(iv) | not for (v) (D3at-5) |
| D3ao-17 | fixed | |

### Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3at-1 | major | **(v)'s balance is wrong.**<br>- The entrant's shift −(w/A)δ that makes the binding O(w²\|δ\|⁴/σ⁴) also cancels the linear pull −wδ·⟨r′∇φ⟩. With the entrant slaved, the force on the tail is w⟨r∇²φ⟩Δ.<br>- The radial curvature is ⟨r∇²φ⟩ ≈ A′‖φ‖²e^{−d²/4σ²}d²/(4σ⁴). The transverse parts cancel against the coefficient's absorption.<br>- So the minimum is at Δ = σ[A′(d/σ)²e^{−d²/4σ²}/μ_r]^{1/2}, μ_r = w_ew/(w_e + w), not σ(σ\|⟨r′∇φ⟩\|/(w‖φ‖²))^{1/3}.<br>- V1 check (d = 7σ, A′ = 0.7): Δ/σ = 0.0298 and 0.0433 at w = 0.3 and 0.107. My law gives 0.0296 and 0.0421; the draft's gives 0.0339 and 0.0479. The local slope over two decades is −0.52, not −1/3.<br>- Unbinding (the fold) occurs at w ≈ 0.13, 0.023 and 0.0015 for d = 5, 6, 7σ. That is 19×, 44× and 64× the stated A′(d/σ)e^{−d²/4σ²}, so it is not "O(1)". | (v): statement and proof. Use the curvature pull and the square root. The level should be A′(d/σ)²e^{−d²/4σ²} times a measured c(d) (about 4–9 here). |
| D3at-2 | major | **Runner-ups that form are captured; (iv)'s count fails near same-sign unfitted features.** All results below agree at rtol 10⁻¹² and 10⁻¹³.<br>- C2 at 3·10⁻⁹: gap 5.48 < x_θ(0) = 6.30, share 1.8% after the entry. It is captured, so feature 0 gets N = 1 where the rule gives 2.<br>- C2's window edge lies between tails 0.0143 (captured) and 0.0231 (bound). That is 35–56× the stated level of 4.1·10⁻⁴.<br>- With feature 2 at 5σ, shares of 10% (10⁻⁸) and 21% (3·10⁻⁸) are captured, 26× the stated level.<br>- 1-D (D‴.51 H(3), c₂ = 0.048): captured at every ε from 4·10⁻¹⁰ to 10⁻⁶ (x down to 1.77, share about 20%); bound at 2·10⁻⁶.<br>- Part C's recipe (features ≥ 3σ apart in [−0.7, 0.7]²) makes neighbours at 4–5σ typical. | (iv): restrict the rule to features with no same-sign unfitted feature within capture range, or state the exclusion. (ii)(c): separation at the entries is not enough. (v): drop "below θ". |
| D3at-3 | major | **"The count per feature by final positions does not differ" is false: a capture tie.**<br>- C2 at ε = 3.13·10⁻¹⁰: feature 2 ends formed by 163 (0.055) and its own entrant 208 (0.645), rms 1.4·10⁻¹¹ (converged).<br>- 1-D (c₂ = 0.043), ε = 1.3–1.4·10⁻⁶: bump 2 at 0.632, which has no near tie of its own (gap ≈ 53), is formed by 11 and 16. The split varies continuously: 0.039/0.379, 0.177/0.241, 0.342/0.076, 0.408/0.010. D‴.51's own H(3) row already lists [10, 11, 16, 3].<br>- Just below the capture (3.0–3.12·10⁻¹⁰) the tail is stranded about 1σ from feature 0, and the fit at τ_end has rms up to 7·10⁻⁶. | (v) and (iv): add the arrival-time tie, where N = 2 occurs at the capturing feature. Add it to C2 as a designed case. |
| D3at-4 | minor | **(v) rests on one capture, at one mass, on one target, at one tolerance.**<br>- Capture outcomes sit on a saddle and flip with the integrator's tolerance.<br>- 1-D at the library default (rtol 10⁻⁹): spurious "bound" at 6 of 22 masses, which also change with atol and the output grid.<br>- 2-D 5σ variant: rtol 10⁻¹¹ says bound, while 10⁻⁹, 10⁻¹² and 10⁻¹³ say captured.<br>- C2's outcomes are robust (checked at 10⁻⁸, 4·10⁻⁹, 3·10⁻⁹, 10⁻⁹ and 3.13·10⁻¹⁰). | C2: scan ε, d and A′. Check tolerance at rtol ≤ 10⁻¹² for every capture case. |
| D3at-5 | minor | **(v) is not new.** D‴.51(ii)(c) already says "with a feature of its sign still unfitted within reach it migrates there and enters". D‴.51(iii) gives its migration scale, and D‴.51's numerics (G, target 2; H(3)) observe it. The revision does not cite this, and the round-42 log calls it "a new mechanism". | (v) and the log: cite D‴.51(ii)(c), (iii), H(3); reconcile with D‴.51(iii)'s criterion. |
| D3at-6 | major | **The n-D product rule and near-factorisation use the planar offset.**<br>- At a vertex of hℤⁿ each edge's transverse offset is h(n−1)^{1/2}/2. For n = 3 that is 0.707h, where x_θ ≈ 5.77 (T2) against 6.03 at h/2.<br>- (iii)(d)'s "the offset h/2" and (iv)'s "p(x; h/2)" hold only for n = 2. This is D3ao-2 again.<br>- There are no n ≥ 3 numerics. | (iii)(d) and (iv): use p(x_i; \|t_i\|), \|t_i\| the feature's distance from edge i's axis. Add a 3-D cube check. |
| D3at-7 | minor | **(iii)–(iv) hold for equal masses only; this is unstated.**<br>- W2's own loops: at the power face's exact tie the split is 0.498/0.502 (β = 1.001), 0.478/0.522 (1.01) and 0.021/0.979 (1.5).<br>- At gap 12.4 the heavier runner-up keeps 0.81% at β = 1.5 against 0.03% at β = 1.<br>- So the layers and x_θ depend on the mass ratio (D‴.51(ii)(b)). | Statement: say "equal masses" in (iii)–(iv), or add the ρ-dependence. |
| D3at-8 | minor | **B4's 213 unsettled problems do not affect the areas, but their cause is misstated.**<br>- All but (5, 5) have N ≤ 2 with third share ≤ 0.0059 (θ = 0.01). (5, 5) has third share 0.0209.<br>- I reran 105: every count is unchanged.<br>- 101 finish all 12 doublings in ≤ 0.5 s (s = 1.16·10⁸), with last drift 1.0–3.7·10⁻⁷ against share_tol 10⁻⁷. They are bit-identical to the raw.<br>- 4 take 139–189 s. Of these, (3.75, 7.25), (3.75, 7.75) and (5, 5) settle on rerun, so their raw flags were budget artefacts.<br>- The numerics' "648 settled within the 300 s budget", the log's "the counts of the rest at the budget", and B's label "(B4's, settled)" (17 of B's 66 points are unsettled) are all wrong. | Numerics (B4), log, report(): state the doubling cap. Quantify the grid and xi0 = 400 errors (about 1%). |
| D3at-9 | minor | **inner_k's "stationary" is a cap and clock flag, not convergence.**<br>- A drift ≤ 10⁻⁷ over one doubling does not bound the remaining change.<br>- max_doublings = 12 decides most cases, and the wall-clock budget makes the rest machine-dependent: (4.75, 8.75) differs from the raw by 6·10⁻⁷.<br>- No test covers any of this (mutants IK1–IK3 survive). | `inner_k`: report why it stopped; estimate the remaining drift from the doubling sequence; add tests. |
| D3at-10 | minor | **D3ao-5 is only partly fixed.** (i) still asserts D‴.51(iv) and D‴.53(c)–(d) in ℝⁿ, with no n-D `recovery_threshold`, `entry_pull` or `loop_with_duals`. D‴.51(iv)'s near-tie heuristic uses a t-independent Gibbs threshold, which (iii) contradicts, and it misses captures. | (i) and Status: restrict to D‴.51(i), (iii) and D‴.53(a)–(b). |
| D3at-11 | minor | **Part D does not resolve x_θ(t).**<br>- An offset-blind rule with x_θ(0) = 6.30 everywhere matches the loop at 596 and 600 of 600; the offset rule matches at 592 and 597. Every miss has the loop above the rule, a finite-ε shift.<br>- "Within 0.3 of x_θ but one" has two exceptions: +0.37 at gap 6.52, and the 4/3 miss decided by a diagonal share of 0.0094.<br>- P(N ≥ 3) still falls only 9 → 3 between the masses, against the law's factor 10. | Numerics (D): stratify near the faces, or give up resolving x_θ(t) from D. |
| D3at-12 | minor | **20 of 24 mutants survive.**<br>- Product rule: offsets swapped (md D 592 → 589), diagonal share replaced, edges independent, p_edge without \|t\| (D 592 → 594). The tests use x₁ = x₂ and t₁ = t₂ only.<br>- Report formulas: the probability prefactor; no B4 symmetrisation (area 33.5 → 17.5); the face-average norm; E[N−1]'s factor; x̄_θ's quadrature.<br>- x_θ(t): `x_theta_t` at 2θ; T2's bisection direction.<br>- All four C2 bookkeeping mutants survive: final positions read from the start, the captured test inverted, n0_fin by start, the unbinding level's factor.<br>- inner_k: share_tol ×1000, no budget, drift measured against the first chunk. Also `margins_diag`'s default ε.<br>- Killed: the half weights in the areas, the face cube [0, H], the edge offset along the axis, and profile_split with max in place of sum. | tests: product_count with x₁ ≠ x₂ and t₁ ≠ t₂; x_theta_t on a synthetic table; face_average against a quadratic x_θ; report() on a synthetic raw; c2_job's bookkeeping on the stored C2 raw; inner_k's stop on a slow tail. |
| D3at-13 | cosmetic | **Text and md mismatches.**<br>- C2 md's "runner-up's share: loop" at 10⁻⁹ is 0.4667: the captured kernel's coefficient at feature 2 over feature 0's kernels, not a share.<br>- C2 judges "formed" by the start feature.<br>- "2.5·10⁻⁸ to 2.4·10⁻¹⁴" are the edge values (the vertex runs 3.8·10⁻⁹ to 3.8·10⁻¹⁵).<br>- "α_j its mirror coordinate": the mirror coordinate is η.<br>- The amplitude A = 0.8 is not stated with x̄_θ and E[N−1].<br>- The (4, 6)/(5, 6) non-monotonicity rests on shares of 0.009801 and 0.010031, i.e. 0.3% from θ. | Numerics, md report(), statement. |

### Lesson
A mechanism found by balancing two terms must expand both under the same constraint. The compensation that made the binding quartic also cancelled the pull it was weighed against. The level came out 1–2 orders of magnitude too low, so the one observed capture looked "O(1)-consistent", and the formed runner-ups that the true level lets go were never looked for. One designed target at five masses cannot validate a threshold law. The checks that exposed the error were: scanning the parameters that move the threshold (ε, d, A′); converging the integrator where the dynamics sit on a saddle (1-D outcomes flipped between rtol 10⁻⁹ and 10⁻¹²); and reading the earlier proposition that had already recorded the phenomenon. D‴.51's H(3) table already held the counterexample to the count claim.
