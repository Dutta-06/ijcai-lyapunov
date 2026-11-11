# Referee D3al (round 39): Proposition D‴.51, the transit at a finite mass

## What I checked and how

- **What was reviewed.** Commit aab77ca, copied to /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review39/repo.
  - I wrote nothing inside the repository.
  - While I worked, HEAD moved to 8ddfaba, and other work left uncommitted edits to `near_tie_inner` and `entry_pull` in the live tree. I did not review those edits.
  - Every run used the venv, OMP_NUM_THREADS=1 and nice -n 10, one heavy process at a time. The container restarted once; I redid the interrupted run.
  - All scripts and logs are in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review39/work/ (a summary is in notes.md there).
- **Tests.** tests/test_transit_mass.py passes, 8 of 8, in 77 s.
- **Claim (i).** I re-derived it by hand and compared 𝒲 with both asymptotes numerically.
- **Numbers.** I checked every number in the NUMERICS paragraph against results/chaos/transit_mass.md and the raw npz. All match.
  - My harness reproduces the committed numbers exactly: Part D's 48 predictions are bit-identical, and Part B's loop Δ values agree to 4 decimals.
- **Loop thresholds by bisection.** I used Part D's own code path (`run_moving` to 1.3 τ_last, θ = 0.01·max|A*|, "recovered" = the path's set of kernels), with brackets of ×1.04–1.09. Targets bisected:
  - all 48 Part D targets;
  - 16 Part F targets: the five that "recover at 3ε_c" and 11 others;
  - single-bump near ties at δ = 10⁻³, 10⁻² and 3·10⁻²;
  - D‴.48's target 10;
  - D2's target with its middle bump moved to 6 positions.
  - Check on the harness: these thresholds give exactly Part D's 184/192 on its decade grid.
- **Diagnostics.** The internals of `entry_pull` on D27, and the trajectories of D2, E115 and E117.
- **Controlled loops.**
  - The runner-up's centre frozen (Part B).
  - Both centres, or one, held fixed until τ_e (the "approach creep").
  - Unequal masses: exact ties made by b = 1.2 (far side) and b = 1.626 (near side).
  - The frozen inner pull at xi0 = 1600 and 3200.
- **Statistics.** On the 400 Part E paths: whether the entrant is its feature's nearest kernel, and the runner-up margin distribution against 4σ²y/h². Bootstrap confidence intervals for Part E's slopes.
- **Mutations.** 34 mutants of the D‴.51 section:
  - eta_crit/inner_coef I1–I5, near_tie_inner N1–N7, path_duals P1–P4, creep_reaches C1–C6, recovery_threshold R1–R12.
  - Each was run against test_transit_mass.py.
  - Survivors were re-scored on Part D (agreement on its decade grid, and agreement with the bisected loop thresholds) or on Part A/B inner quantities.
- **Literature.** Web search only; arxiv.org and neurips.cc were unreachable. Statements marked [memory] are from memory.

## Verdicts

### (i) The inner profile: correct, with two misstatements

- **Confirmed by re-derivation:** η_ε(α) − Λρ → η_c = (α − a coth(α/a))/(2Gc); K'³ = (8/π³)a/(m(2c)^{1/2}) = 0.0962; κ = πρ(2c)^{1/2} = 0.7025; the exponential branch (ca/2m)e^{2α/a}; the overshoot slope 1/κ = 1.4235; and the overshoot values 3.65…8.90.
- **The rate is ε^{1/2}, not ε.** The proof's "λ₁ − a²/(2Gcα) + O(ε)" holds for Λρ, but λ₁ − Λρ = O(ε^{1/2}). The second integral also converges only at O(ε^{1/2}), because its ε-dependence sits at α ~ ε^{1/2}.
- **"Power law off by more than 20% below ξ ≈ 15" is false.** 𝒲/(K'³ξ⁻³) − 1 is:

  | ξ | 15 | 10 | ≈8.3 | 5 | 3 |
  |---|---|---|---|---|---|
  | 𝒲/(K'³ξ⁻³) − 1 | −6.7% | −14% | −20% | −40% | −65% |

  - Equivalently 𝒲 ≈ K'³ρ⁻³ξ⁻³(1 − 16.2/ξ²).
  - The exponential branch is within 2.7% from ξ = 0 downward.
  - Consequence: `creep_reaches` uses the core law down to ξ = 3, where it overstates the coefficient by 65%.

### (ii) The near-tie law: right only for symmetric ties; the loser is not frozen

**What holds.**
- The leading order at fixed x is fine for the symmetric single bump.
- x_θ from the inner problem is 6.41, 6.30 and 6.27 for A = 0.3, 0.8 and 2, against the Gibbs value 6.54.

**"The gap shrinking a little during the slide" has the wrong sign.**
- The inner share is below the Gibbs split at every x tested: 0.437 vs 0.439, 0.309 vs 0.316, 0.073 vs 0.080.
- At A = 0.8, x = 6.5 the locked gap is 6.79, i.e. it grew.

**The finite-margin correction is the partners' approach creep, and it is not small.**
- In the loop (single bump, A = 0.8) the formation threshold is x_θ = 6.39, 7.10 and 8.19 at μ* = 0.0035, 0.034 and 0.099.
- With both centres held fixed until τ_e it is 6.56, 6.55 and 6.52 — the Gibbs value at every margin.
- Fixing the approach also cuts Part A's worst excess (δ = 10⁻², x ≈ 11) from 1.39–1.55 to 1.16–1.22.
- So the law is exact only when both partners creep equally, which in practice means a symmetric geometry.

**𝒮(0) = ½ needs equal creep, not just a tie.**
- Exact ties made by unequal masses give 𝒮 = 0.885 (far side: b = 1.2 at distance 0.076 against b = 1 at 0.024) and 0.997 (near side: b = 1.626 at 0.12 against 0.02).
- Both values are the same at ε = 10⁻⁸ and 10⁻¹⁰.
- In the near-side case the loop forms one kernel where the path forms two.
- The draft does say ½ "at a symmetric tie". But it presents the Gibbs split ("exact at x = 0") as the near-tie law, and (iv) applies it to every partner with 0 < κ_ij ≤ 1.

**Near-side runner-up (E117).**
- Kernel 14 sits between the entrant 13 and bump 3, and is the bump's nearest kernel.
- It never forms at x = 5.2, 1.6 or 0.52: max|w₁₄| = 10⁻⁷, 4·10⁻⁶ and 6.5·10⁻⁴, against θ = 9·10⁻³.
- The Gibbs split predicts a 41% share at x = 0.52.
- What happens instead: the entrant slides through it (min|c₁₃ − c₁₄| = 0), and the two never lock.

**"Far behind, the runner-up stays out with the frozen coefficient" is false whenever features remain unfitted.**
- The loser's coefficient 𝒲(x − Δ) does not depend on ε, so on the τ scale its centre moves at ε^{−1/2}Λ₀Gm𝒲⟨r∂_cφ⟩.
- In D2 the loser of the first entry (x = 16.8–18.4, |w| ≈ 5·10⁻⁵) travels from 0.1 to 0.632 and takes the unfitted same-sign bump 3. The loop fails at 1.2·10⁻⁶ against ε_c = 8.85·10⁻⁶.
- D2 sweep: migration thresholds of 1.27·10⁻⁶, 9.6·10⁻⁸ and 3.3·10⁻¹⁰ at μ* = 0.0189, 0.0087 and 0.0018, i.e. ∝ μ*^{3.3–3.6}.
- These are ×7, ×20 and ×240 below the entry law's (μ*/6.3)².
- So "costs a kernel unless μ* ≫ ε^{1/2}, not ε ≪ μ*³" fails here; the binding condition is close to D‴.48's cubic one.

**Proof gaps.**
- The sublevel set contains {w_r = 0, w_e = A, c_e = c*}, where E = 0 and |𝒲′(ξ_r)| vanishes. So "away from the zero coefficients, where the metric is bounded" is not established.
- For x ≳ 34 the computed "end" is a plateau with the runner-up still at its start and still drifting (A2: −9·10⁻⁴ in share from s = 4.7·10⁶ to 1.9·10⁷). That is not the Łojasiewicz limit, which has c_r = c*.
- "Theorem D‴.53 Steps 3–4 with both kernels in the layer" does not cover the time after the entry: D‴.53's outer step assumes the non-block kernels are O(ε^{3/2}μ₀⁻³), while the loser is O(x⁻³).

**Part B's creep explanation is confirmed.**
- With the runner-up's centre frozen in the loop, Δ_loop − Δ_inner(δ) = +0.0014 to +0.0018 at x = 34.7, −0.0003 at 110 and −0.0018 at 347, for A = 0.1, 0.8 and 4.
- At x = 347, loop(free) − loop(frozen) = 0.0366–0.0371, which is 0.997–1.011 of the formula.
- The draft's "0.95–0.96" is an error in its reference value, not in the formula:
  - the Richardson extrapolation from xi0 = 400/800 overestimates Δ∞ by 0.0010 (converged values: 2.2776, not 2.2786; 7.2970, not 7.2980);
  - there is also an offset of about 0.0008 between loop and inner problem.
- The A-dependent excess at x ≈ 35 and x ≈ 110 is a genuine 1/x² term, and the draft does not explain it:
  - at x ≈ 35 the ratios are 1.023, 1.114 and 1.195; at x ≈ 110 they are 1.009, 1.037 and 1.060;
  - x·(ratio − 1) ≈ 0.9, 4 and 7 for A = 0.1, 0.8 and 4, at both x.
- The formula's mass factor is wrong: it should be a²b/(π²c), not a²b^{3/2}/(π²c), because dμ/dτ = −g/ρ contributes a factor ρ.

### (iii) The creep: the μ⁴ algebra is right; the scale is not predictive; the threshold criterion is flawed

- **Target 10.**
  - The first-order formula gives anywhere from 2.6·10⁻⁹ to 4.2·10⁻⁸, depending on which μ, g and T one picks.
  - `creep_reaches` gives 4.88·10⁻⁹.
  - The loop's threshold is 3.82–3.99·10⁻⁹, i.e. 0.80 ε_c. Part C's bracket (0.72–1.43 ε_c) hides that the prediction is 20% on the unsafe side.
- **(a) The criterion fires on kernels that merely touch the band at another kernel's entry.**
  - E115: the margin bottoms out at 6.2·10⁻⁴ at the third entry, giving ε_creep = 3.85·10⁻¹⁰. The loop recovers up to 4.4·10⁻⁵, a factor of 1.15·10⁵.
  - The other four "recover at 3ε_c" targets are the same kind: E117 ×2.3·10³, E8 ×9.1, E216 ×9.9, E267 ×225.
  - In all five, the binding kernel's minimum margin sits at an entry, with x = 9–32.
  - The predictor's own entry analysis says these kernels do not form.
- **(b) It ignores that a path kernel's own creep can fit the driving feature earlier.** D15 ×94.5, D28 ×25.8, D45 ×16.7, D46 ×19.9. In these the loop enters the path's kernels in a different order, and the predicted kernel never forms.
- **(c) It starts a near-tie loser from its unpulled margin.** D2, ratio 0.14.
- **Otherwise it is mildly unsafe:** 0.77–0.92 on 7 small-margin Part F targets, and 0.32–0.92 on 11 Part D targets.
- So "errs on the safe side" is wrong as a summary: the errors are bimodal.

### (iv) The predicted mass: false as stated and biased toward the unsafe side

**Bisection on all 48 Part D targets.**
- The loop fails strictly below ε_c on 32 of 48 targets, and already at ε_c itself on D0, D3, D19, D39 and D40.
- Entry-bound targets (32): median ratio ε_loop/ε_c = 0.65, range 0.13–1.09, plus one false alarm (D27, ×43).
- The bias grows with the margin:

  | ε_c range | ratio ε_loop/ε_c |
  |---|---|
  | ≤ 2·10⁻⁶ | ≈ 0.92 |
  | 10⁻⁵ to 10⁻⁴ | 0.55–0.92 |
  | above 10⁻⁴ | 0.27–0.65 |
  | D39, ε_c = 8·10⁻⁴ | 0.39 (the loop forms 6–11 kernels there) |

- Creep-bound targets (16): 11 at 0.32–0.92, and 5 at ×5.4 to ×94.5.
- Only 20 of 48 are within ×1.5, 36 within ×2, and 9 are off by more than ×3.

**Single near tie.**
- ε_loop/ε_c = 0.96, 0.74 and 0.52 at μ* = 0.0035, 0.034 and 0.099.
- The predictor's x for the partner is max(Δ + ξ_θ, x_θκ_ij) = max(5.05, 6.54κ_ij) = 6.25, 6.14 and 5.85.
- The loop needs 6.39, 7.10 and 8.19, so κ_ij moves x the wrong way.
- The frozen-tracer pull misses the partner's own creep and its share of the bump. The ad-hoc floor x_θκ_ij is what actually carries the prediction.

**The decade grid hides all of this.**
- My bisected thresholds reproduce Part D's 184/192 exactly.
- Ablations score almost the same on that grid:

  | Mutant | Decade-grid agreement | Other effect |
  |---|---|---|
  | R6: ε_c = ε_creep (entry threshold dropped) | 183/192 | 15/48 within ×2 of the loop (committed: 36/48) |
  | R5: x_θ × 1.5 | 177/192 | unsafe count falls from 32 to 11 |
  | C3: creep cut-off 1·ε^{1/2} | 186/192 | better than the committed predictor |

**D27, a false alarm.**
- The measured pull Δ = 70.5 on kernel 10 at the entry of 15 (with 171.7 on kernel 11 and −850 on kernel 17) is the block's slow relaxation, integrated to s = 1.19·10⁵.
- The next entry comes 0.27 later in τ, i.e. s ≈ 270 at ε = 10⁻⁶.
- Result: ε_c = 6.2·10⁻⁷ against the loop's 2.7·10⁻⁵.
- The earlier predictor in transit_mass_raw_v1.npz gave 1.84·10⁻⁵ and agreed. Part D fell from 186/192 to 184/192, and the draft does not report this.

**"Out of sample".**
- Part D's seed is fresh.
- But the measured pull was introduced after Part D and Part E had been run, and its regression test uses Part E's target 327.
- Part C is the tuning target.

**Unequal masses are not handled.** `path_duals` pins duals at ±1, and `creep_reaches` and `recovery_threshold` use 1 − |u| with no b. The proposition is stated for masses εb_i.

### (v) Statistics: these are the predictor's, not the loop's

- **The per-feature margin law checks out.**
  - P(runner-up margin < y) = 0.0125, 0.030, 0.051, 0.106, 0.263 at y = 0.002, 0.005, 0.01, 0.02, 0.05.
  - 4σ²y/h² gives 0.0115, 0.029, 0.058, 0.115, 0.288.
  - The h/σ = 0.83 objection does not matter at small μ: there μ is linear in h − 2|δ| whatever h/σ is.
- **But the premise "entrants are the nearest kernels" fails on 265 of 1200 entries.**
- **Small margins are not only runner-ups.** The smallest margin of any never-entering kernel is 2–3 times as frequent: 0.033, 0.072, 0.114, 0.24, 0.53 at the same y.
- **Slopes (bootstrap 95%):** all targets 0.294 [0.24, 0.37]; near ties 0.462 [0.39, 0.57]; creep 0.285 [0.24, 0.36].
  - These are consistent with ε^{1/4} and ε^{1/2}, but only for the predicted fractions.
  - The loop's failure rate over the ensemble was never measured; Part D's 16, 10, 6 and 1 out of 48 are too few.
- **"Mass ∝ h⁴" is untested.** No experiment varies h. If the creep that is said to dominate acts at margins ∝ h², its threshold would scale like h⁸.

### Mutations: 18 killed, 16 survive

**Killed.**
- I1–I4 by the inner-profile test.
- N3 and N4 by test_near_tie_inner_problem.
- N1/N2 (slide factor ×½ or ×2) and N7 (matched pull without s₀(d_e − d_r)) only by the internal-consistency test between `entry_pull` and `near_tie_inner`.
  - The physical loop test passes N1/N2: the shares move by ≤ 10⁻⁴, while the frozen pull goes from 5.12 to 5.47 or 4.65.
- P3 by test_path_duals.
- C1, C6 and R7 by the single-target test on target 10.
- R1, R2, R3, R8 and R12 only by the single-target test on target 327.

**Surviving with no effect on any D‴.51 output.**
- I5 (ρ dropped; ρ = 1 everywhere).
- P1, P2, P4 (pinning removed or set to 0.99, next interval's centres). The never-entering kernels' duals are identical, so this is dead code for D‴.51.

**Surviving with small effects.**
- N5/N6 (xi0 = 100 or 1600): shares change by 0.03–0.8%, the pull goes from 5.12 to 5.28 or 5.08.
- C5, R4, R9, R10: change 6, 28, 18 and 5 of the 48 predictions, with no change in aggregate scores.

**Surviving and changing Part D's results.**

| Mutant | Decade-grid agreement | Within ×2 of the loop |
|---|---|---|
| R6 (entry threshold dropped) | 183/192 | 15/48 |
| R5 (x_θ × 1.5) | 177/192 | — |
| C4 (cut-off 10·ε^{1/2}) | 172/192 | 13/48 |
| C3 (cut-off 1·ε^{1/2}) | 186/192 | 35/48 |
| C2 (creep × 1.5) | 184/192 | 38/48 |
| R11 (n_creep = 1) | 181/192 | 34/48 |

C3 and C2 come out at least as close to the loop on Part D as the committed constants (C3: 186 vs 184 on the grid; C2: 38 vs 36 within ×2). The constants are not pinned by any test.

### Novelty

- **Pesme–Flammarion and Berthier.** Asymptotic results, with no finite-initialisation near-tie analysis that I could find [memory; I could not re-read the appendices either]. The lit note is fair here.
- **Gissin–Shalev-Shwartz–Daniely (ICLR 2020)** [memory]. For depth ≥ 3, incremental learning needs the relative gap to be large compared with a power of the initialisation. At depth 3 with w₀ ∝ ε^{3/2} (D‴.48(o)) that is a margin ≫ ε^{1/2}, which is D‴.51's margin. Not cited.
- **Exact regularisation (Friedlander–Tseng 2007 [memory]; Yin, linearized Bregman, SIIMS 2010).** There is a finite smoothing threshold, set by the dual certificate, above which the exact sparse solution is recovered. This is the static analogue of a "mass for exact recovery read off the limit path". Not cited.
- **Osher–Ruan–Xiong–Yao–Yin (ACHA 2016).** The damped (finite-κ) linearized-Bregman ISS is a smoothing family of the inverse scale space, with sign-consistency results at finite κ. This contradicts lit_S76's "not a smoothing family".
- **AGF (Kunin et al.).** Flagged in round 41; should be positioned against.
- **Hedge/EG and Duval–Peyré.** Correctly acknowledged.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3al-1 | major | (iv) is false as stated. The loop fails below ε_c on 32/48 Part D targets. Entry-bound median ratio is 0.65, falling to 0.13–0.27 at larger margins; only 36/48 are within ×2. Part D's decade grid cannot see this: dropping the entry threshold still scores 183/192, x_θ×1.5 scores 177. | D51_draft (iv): state ε_c as an estimate with its measured error band. exp_transit_mass D/F: bisect the loop's threshold and report the distribution of ε_loop/ε_c. Drop the "safe side" narrative. |
| D3al-2 | major | The finite-margin correction to the near-tie law is the partners' differential approach creep. Freezing both centres until τ_e gives x_θ = 6.52–6.56 at every margin; free, it rises to 7.10 and 8.19. This creep is the cause of D3al-1's bias at larger margins. It is also why asymmetric exact ties do not split ½: 0.885 and 0.997 with unequal masses. | (ii): add the O(1/x) creep term with coefficient ∝ b_r(g′_r/g_r)² − b_e(g′_e/g_e)², and restrict "Gibbs, exact at x = 0" to equal creep. `near_tie_inner`/`entry_pull`: pose the inner problem at the actual geometry, with the pre-entry creep and the tracers' 𝒲. Part A: add asymmetric and unequal-mass cases. |
| D3al-3 | major | A near-tie loser is mobile at ε^{−1/2}Λ₀Gm𝒲⟨r∂φ⟩ and migrates to an unfitted same-sign feature. In D2 it travels 0.1→0.632 and the loop fails at 0.14 ε_c. Thresholds scale ∝ μ*^{3.3–3.6}, ×7–240 below (μ*/6.3)². "Stays out with the frozen coefficient" and "not ε ≪ μ*³" are false in these configurations. | (ii)/(iii) text. `creep_reaches`: after each entry restart the never-entering kernels from their pulled margin μ − ε^{1/2}Δ_ij. Add the D2 sweep as an experiment. |
| D3al-4 | major | `entry_pull` integrates to stationarity whatever the time to the next entry. In D27 a slow block relaxation becomes a "pull" of 70.5, giving ε_c = 6.2·10⁻⁷ against the loop's 2.7·10⁻⁵. The v1→current change made Part D worse (186→184) without being reported. | `entry_pull`: cap s at (τ_{k+1} − τ_k)/ε^{1/2}, which makes Δ depend on ε, or flag relaxations that are not separated from the next entry. Numerics: report the predictor's revision history. |
| D3al-5 | major | The `creep_reaches` criterion (3ε^{1/2} band before the last event, along the path's residuals) fires on kernels whose margin bottoms out at another kernel's entry. All five Part F "recover at 3ε_c" are of this kind, off by ×9 to ×1.15·10⁵. It also misses path kernels that creep in early (×17–95) and is mildly unsafe otherwise (0.32–0.92). | `creep_reaches`: test only the creep-induced deviation; require ξ < ξ_θ; stop at the entry that fits the driving feature (use that entry's inner problem). Part F: report the bisected factors. |
| D3al-6 | major | Near-side runner-up: E117's kernel 14 does not form at x = 0.52, where Gibbs predicts 41%. The entrant slides through it, so the lockstep premise fails. | (ii): state the far-side, coinciding-and-staying premise. Add E117 as a test. |
| D3al-7 | minor | "The gap shrinking during the slide" has the wrong sign: the inner share is below Gibbs everywhere; at x = 6.5 the gap grows to 6.79. | D51_draft (ii). |
| D3al-8 | minor | The runner-up creep formula's mass factor should be b, not b^{3/2}. ρ ≠ 1 is never exercised (I5 survives). | Draft (ii)/status; exp_transit_mass b_run comment; add a ρ ≠ 1 test. |
| D3al-9 | minor | (i): the power law is 6.7% off at ξ = 15, not >20% (20% at about 8.3). The proof's O(ε) is O(ε^{1/2}). | Draft (i) and proof; inner_coef docstring. |
| D3al-10 | minor | Part B: the "0.95–0.96 at x ≈ 347" is the reference's error (Richardson +0.0010 plus a 0.0008 offset); the loop with the runner-up's centre frozen gives 0.997–1.011. The A-dependent 1/x² excess at x ≈ 35 and 110 is left unexplained. | Numerics B: converge xi0 to 1600/3200 or fit three terms; add the frozen-centre control; state the 1/x² term. |
| D3al-11 | minor | (ii) proof: the sublevel set includes w_r = 0, where the metric degenerates. For large x the computed end is a plateau, not the limit (A2 is still drifting). The D‴.53 Steps 3–4 citation does not cover the loser after the entry. | Proof of (ii): define 𝒮 as the plateau at s ≲ ε^{−1/2} and prove it; drop or extend the D‴.53 citation. |
| D3al-12 | minor | (iii): the first-order scale is not predictive (2.6·10⁻⁹ to 4.2·10⁻⁸ on target 10, against the loop's 3.9·10⁻⁹). Part C's ε_c is 20% unsafe. `creep_reaches` uses the core law down to ξ = 3, where it is 65% too large. | (iii) wording; use 𝒲 inside ξ ≲ 30; Part C: bisect the threshold. |
| D3al-13 | minor | (v): Part E fractions are the predictor's, never the loop's. Entrant ≠ nearest kernel on 22% of entries. Other kernels supply half the small margins. h⁴ is untested and conflicts with creep dominance (h⁸). | (v) text; Part E: run the loop on a random subset at a few ε, or bisect it; drop "h⁴" or test it. |
| D3al-14 | minor | Out-of-sample history: the measured pull was designed on Part E's target 327 after Part D/E were run; Part C is the tuning target. | Numerics D/E: disclose this. |
| D3al-15 | minor | κ_ij is an underived heuristic that moves the partner's x the wrong way. The frozen-tracer Δ_ij + ξ_θ (5.0–5.2 for the single near tie) is far below the inner problem's 6.3, so the ad-hoc floor x_θκ_ij carries the prediction, and only for 0 < κ_ij ≤ 1. | `recovery_threshold`: derive κ_ij or replace it by the two-kernel inner problem with creep. |
| D3al-16 | minor | Tests are weak: 16 of 34 mutants survive, including dropping the entry threshold, x_θ×1.5, the cut-off, creep×1.5 and n_creep. A slide factor ×2 passes the physical tests. The pinning is dead code for D‴.51. There is no test against a bisected loop threshold and none with ρ ≠ 1. | tests/test_transit_mass.py: add bisected-threshold tests on 3–4 targets, asymmetric and unequal-mass ties, and tighter loop-vs-inner tolerances. |
| D3al-17 | minor | Novelty: missing Gissin et al. (the depth-3 condition is the ε^{1/2} margin), exact regularisation (Friedlander–Tseng; Yin), and LBISS at finite κ (Osher et al. 2016), which refutes "not a smoothing family"; AGF is not positioned. | lit_S76.md; the draft's Classical/New paragraph. |
| D3al-18 | minor | The code covers equal masses only, but the proposition is stated for εb_i. | `path_duals`/`creep_reaches`/`recovery_threshold`: take ρ, or state the restriction in (iv). |
| D3al-19 | cosmetic | (1) The lit note's temperature a²/(4Gc) should be a/(4Gc). (2) The `near_tie_inner` docstring says xi0⁻²; the data show 1/xi0. (3) The `recovery_threshold` docstring cites "(iii)" instead of (iv). (4) The code uses abs(gi) where the draft writes s_i⟨Δr_iφ_i⟩. (5) The exp_transit_mass B docstring says ε = 10⁻¹⁰, but EPS_BL has three values. (6) Part E's "creep binds on most": the columns count ε ≥ ε_creep, not where the creep binds. (7) D39's ε_c = 8·10⁻⁴ is outside any asymptotic regime. | Respective files. |

## Lesson

A threshold predictor checked at masses a decade apart can agree on 96% of runs while being biased by a factor of about 2 on most targets and wrong by 10–10⁵ on others. The same grid scores a predictor with half its mechanism removed at 183/192. The decisive check is cheap: bisect the loop's own threshold and report ε_loop/ε_c, alongside ablations that must score worse. Here that check exposes three structural gaps:
- the leading-order inner problem, posed at the symmetric tie, is applied at margins of 10⁻²–10⁻¹, where the partners' asymmetric creep shifts x_θ by up to 25% (by 30% at μ* = 0.1);
- a loser in the inner layer keeps an ε-independent coefficient, so it remains mobile on the outer time scale;
- entry and creep are computed separately and combined by a minimum, although on the loop they compound.

A leading-order law at fixed x should be validated where the predictions are actually used, with the next-order mechanism isolated by a controlled run (freezing centres before τ_e) before it is folded into a predictor.
