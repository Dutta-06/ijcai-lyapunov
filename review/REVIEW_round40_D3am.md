# Referee D3am (round 40): Proposition D‴.50, the rich transit under label noise

I modified nothing in the repository and ran only read-only git commands (log, show, diff --stat, status). All runs were in `SCRATCH/repo`, a `cp -r` copy whose HEAD is 0369345: the code is at 284bde6 and the experiment and results at 0cf2e9b. Every Python process ran with `OMP_NUM_THREADS=1` and `nice -n 10`, one heavy process at a time, each under 15 min. EXP_CACHE was not set. My scripts and logs are in the session scratchpad's `review40/work/`.

## What I checked and how

**1. Statement and proof.**
- I re-derived (i)–(iv) by hand: the sign and factor N of ∇F, the Itô/Stratonovich bookkeeping, the projector P, the critical-curve asymptotes, and the first-order jump at an entry.
- I compared the draft against D‴.45 and D‴.48 in PARADIGM.md.

**2. Tests.** `tests/test_dual_langevin.py` gives 7 passed in 121 s.

**3. Mutation testing** (`mutate.py` on a private copy of `src/`).
- 28 mutants covering `mirror_update`, `gibbs_rwm`, `FastFree`, `potential`, `fd_covariance`, `ito_q_drift`, Heun, `tangent_projectors` and `linear_response`.
- I evaluated the response and metric mutants on Part C's stored runs (`cmut.py`, `cblocks.py`).

**4. Part C reproduction and decomposition.** I re-ran Part C with the experiment's own noise arrays (`cdiag*.py`, `cprof.py`).
- The replica reproduces the committed numbers exactly. For example, ε = 10⁻⁶, σ = 0.3, path 0 gives 0.389 / 0.0263 / 0.062.
- I added the exact tangent-linear model (TLM) of the finite-ε scheme. This splits the reported error into a nonlinear part and a "finite-ε model" part.
- I split the TLM further into tangent and normal parts.
- I computed the error profile against the distance to the events at ε = 10⁻⁶ (4 paths) and 10⁻⁸ (2 paths), for σ = 0.1 and 1.
- I bootstrapped the share slope over the paths.

**5. Part B.**
- Design 4 over a long window (`b4long.py`): t = 4·10⁵, which is 18 relaxation times of its slowest mode, at h = 0.5, σ = 0.3 and 0.03. I compared against both T·D(w_BP) and the time-averaged T·E[D].
- Design 5 on 4 independent noise paths (`dseeds.py`). The experiment's seed 105 reproduces its row.

**6. Part D.** Design 3 at ε = 10⁻⁶, σ = 0.3 on 6 independent paths. The experiment's seed 7 reproduces its z-scores (2.01, 1.33, 1.05).

**7. Part A.**
- First-order Laplace prediction of E_π[w] − w* and of the sd (`laplaceA.py`).
- The exact integration-by-parts identity Cov_π(Φw) = T·Φ E_π[D] Φᵀ, from a mirror-update run (`fdexact.py`).
- Metropolis integrated autocorrelation times (`rwm_iat.py`).
- Growth law of the Heun conserved-direction drift, from the stored raw data.

**8. Further checks.**
- Eigenvalues of ∇²F at design 3's plateau, to judge whether the Laplace expansion is uniform.
- Finite-difference check of dw/dη against `FastFree.dw` (agrees to 3·10⁻⁷).
- Every number in `numerics_D50.md` against `results/chaos/transit_noise.md` and the raw npz.

**9. Literature.** Web search found the Hsieh et al. 2018 abstract (mirrored Langevin targets the dual pushforward). Krichene–Bartlett give a Gibbs invariant law only for Euclidean stochastic gradient dynamics. I found nothing on Bayesian readings of inverse scale space. arXiv full texts were unreachable, so everything else on literature below is from memory.

## Verdicts

**(i) The Gibbs law: correct; proof incomplete in one place; misses a stronger exact identity.**
- **Correct parts.**
  - ∇F = Φw − y holds on every free curve. ψ*′ = w is consistent with α(η; q), and the signs and the factor N are right.
  - The ζ-equation is exact for the Stratonovich reading on the curves.
  - The mirror update is the Euler–Maruyama scheme (ULA) of that equation.
  - The paradigm's Euler step on (α, log β) is the Itô reading. Q drifts under it, by −2.2 to −3.8 over t = 2000 in Part A, so it does not have this law.
  - The diffusion is overdamped Langevin at T = σ²/(2N).
  - F is convex and superlinear on free curves, since w is unbounded both ways.
  - E_π[Φw] = y and the pathwise identity are correct. The exact mean fit was already in D‴.45(iv)(a).
- **The convergence proof.** Bobkov's Poincaré inequality covers only initial densities in L²(π), not "every start". The claim is true by another route: inf_α D_i > 0, so ∇²F ⪰ δΦΦᵀ and Bakry–Émery or W₂ contraction applies. The rate is O(min D) ≈ O(ε²), and the real mixing time is Λ²N²/σ².
- **(b) is correct at fixed ε but not uniform.** At design 3's plateau the three face eigenvalues of ΦDΦᵀ are 2.2·10⁻¹², 3.9·10⁻¹¹ and 6.2·10⁻¹⁰.
  - The Laplace standard deviations along them are 350–5 800 at σ = 0.03 and up to 58 000 at σ = 0.3. The face size is Λ = 7 025 and the nearest wall slack is 0.14.
  - So (b) needs roughly T ≲ ε, while Parts B and D run at T/ε = 75–7 500.
- **(b) fails quantitatively at Part A's T = 0.125.** The first-order mean shift is off by ×1.4–3.7, the leading-order sd by −39% to +45%, and T·ΦD(w*)Φᵀ by up to 135%.
- **The exact identity the draft misses.** Integration by parts gives E_π[∇F∇Fᵀ] = T·E_π[∇²F], that is Cov_π(Φw) = T·Φ E_π[D] Φᵀ at every noise level. In Part A it holds to 4.6% (sampling).

**(ii) Fluctuation–dissipation on a support: the result holds for well-sized coefficients; the proof and a hypothesis are missing.**
- The algebra P_SS → I is right. The stated limit of range(D^{1/2}Φᵀ) is dimensionally wrong when |S| < N, but the conclusion survives.
- **The derivation interchanges T → 0 and D_off → 0.** In the regime the proposition claims (fixed T, rich limit) the Laplace premise fails along the face. The result instead follows exactly from the identity above: Cov_π(w_S) = T·E_π[D_S] + off-support terms. This gives exactly zero correlations when |S| = N.
- **Missing hypothesis: T·D_i ≪ w_i².** D varies on the scale |w_i|. Design 4 has w_BP,3 = 0.0018. At σ = 0.3:
  - predicted sd/|w₃| = 0.59;
  - E[D₃]/D₃(w_BP) = 0.83;
  - skewness 0.73, and |w₃| drops to 2.6·10⁻⁴;
  - Var/(T·D) = 0.44 over 18 relaxation times (0.65 at σ = 0.03 on the same path).
  A longer window fixes the well-sized coefficients (1.02, 1.14, 1.02, 0.97, 1.03) but not this one.
- **Design 4's ratio 1.63 is sampling error, not non-equilibration.** Starting from the mode in the linear regime, the expected sample covariance is at most the stationary one (Loewner order), so non-equilibration cannot push a ratio above 1.
- **Design 3 confirms (ii).** Over 6 independent paths the ratios are 0.987 ± 0.033, 1.010 ± 0.019 and 1.006 ± 0.059. |corr| is 0.004–0.042, which is sampling noise; the design's cosines are not a relevant null.
- **Design 5's "1.00–1.05 at every σ" is one path.** On 4 independent paths the coefficient with w = −0.022 gives 0.58, 0.62, 1.04 and 1.36, and max|corr| is 0.053–0.135.
- The asymptotes are verified: 4Gc/a_max = γa_max, and γ + 2Gc/a_max² = 6Gc/a_max².

**(iii) The rich limit: formally consistent, numerically untested.**
- The epi-limit −yᵀν + χ_P, the Laplace principle at temperature T/Λ → 0, the near-uniform law on the face and the crossing time Λ²N²/σ² are all consistent. My own estimate gives E_π[w_off] = O(T·ε^{1/2}) → 0.
- "An unbiased sparse estimate" is a misnomer.
- The O(T^{3/2}) penetration uses the small-α branch, which is not in force at σ = 1: β_j ≈ γmT = 0.1 gives α_j ≈ 1.15.
- **Part D does not sample π.** Its window covers 2·10⁻⁸ to 2·10⁻³ of the crossing time; at ε = 10⁻⁶, σ = 0.3 the crossing time is 2·10¹⁰.
  - The mean agreement is the scheme's exact pathwise identity: w̄_S − w_BP = Φ_S⁺(NΔζ − σΔW)/t.
  - All six rows are one Brownian draw scaled.
  - On 6 paths the z-score spread is about 1, so the batch standard errors are credible, but "within 2.3 se" is one sample.
  - The inactive coefficients ∝ ε^{3/2} only show that the boundary was not reached. They do not contradict O(T^{3/2}).

**(iv) The noisy transit.**
- κ = σ/(NΛ^{1/2}) is right.
- **The first-order response is confirmed sharply.** The TLM's tangent part equals κP_{T(τ)}W(τ) to 10⁻⁵ relative.
- **The residual outside the windows is finite-ε, as the paragraph says.** It is the active block's unabsorbed normal Ornstein–Uhlenbeck component:
  - rms 4.03·10⁻⁵ against the predicted κ(Σᵢ 1/(2Λλᵢ))^{1/2} = 4.09·10⁻⁵ at ε = 10⁻⁶;
  - 4.09·10⁻⁶ against 4.09·10⁻⁶ at 10⁻⁸;
  - its relative size scales like (Λλ)^{-1/2}·√log.
  The paragraph does not state this mechanism or rate.
- **The windows are misexplained.** At the tested ε the failure inside the windows is linear in the noise, not nonlinear.
  - The sup-error column at 10⁻⁶ is σ-independent: 0.694, 0.696, 0.699.
  - Actual vs TLM is 0.04–0.19 in the windows; TLM vs response is 0.31–0.69.
  - The window width is the entering kernel's finite-ε absorption time, ∝ 1/(Λ·D_j) ∝ ε^{1/2}: about 0.02–0.03 at 10⁻⁶ and 0.002–0.005 at 10⁻⁸ for σ = 0.1.
  - The noise-induced entry shift at σ = 0.1 is far smaller (7·10⁻⁴ to 3·10⁻³ at ε = 10⁻⁶, 5–9·10⁻⁴ at 10⁻⁸).
  - The noise window, about 0.01, dominates only at σ = 1 and ε = 10⁻⁸; at 10⁻⁶ and σ = 1 the two widths are comparable.
  - The proof's "rates O(1) in t" fails just after an entry, where D_j is small.
- **The rate σε^{1/4} is consistent.** The slope is 0.267 with bootstrap 95% CI 0.241–0.291, but from three different path sets; "at every σ" is the same 4 paths.

**Numerics A–D, by question.**
- **A — Metropolis independence.** The sampler is unbiased: the proposal is fixed and symmetric, and it starts at the mode. It is independent of the time step but shares the `FastFree` code.
- **A — batch-means standard errors.** Credible. The Metropolis autocorrelation times are 10–19 steps against 45 000-step batches; the mirror update's are about 7 time units against 900.
- **A — 2.4 combined standard errors at h = 0.005.** This is chance: it is the maximum over 8 coefficients, and it is larger than at h = 0.01, which an O(h) bias would not produce. The O(h) bias cannot be resolved by this comparison.
- **A — Heun drift.** A scheme artefact, |P_{ker Φ}η| ≈ (0.4–1.3)·h√t, with a ratio of 4.06 between the two steps at t = 2000. At fixed h there is no stationary law at all, rather than "another law".
- **A — Itô rates 0.54–1.09.** Not a test: one path over one time unit while the state moves. The unit test is the real test.
- **B, C and D** are covered above.

**Code and tests.**
- **Killed.** The integrator test pins `mirror_update` exactly: noise ×√N, ×1.1 and ×1.05, drift sign, missing N and ×1.1 are all killed. Also killed:
  - FastFree w, ψ* and `potential` mutants, by the gradient test;
  - dw mutants, by the library test;
  - `fd_covariance`, `ito_q_drift`, and Heun→Euler;
  - `gibbs_rwm` with T removed or 2T.
- **Survived.**
  - `gibbs_rwm` acceptance at 1.1T and 1.2T.
  - Newton tolerance 10⁻⁶ (equivalent).
  - Every `tangent_projectors` / `linear_response` mutant: no test calls these functions.
- **What the experiment catches instead.** Part C's windowed error moves from 0.054/0.026 to 0.64–1.73 for every response mutant except the interval-start segment lookup, which is equivalent outside the windows. These mutants are:
  - normal projector, identity, tol 0.5;
  - no projection of the accumulated deviation;
  - segment −2 or +0;
  - κ×2.
- **Metric sensitivity.** Shifting the noise blocks by one step changes the metric only from 0.0263 to 0.0269 and from 0.0473 to 0.0446. At its 2.6–5% floor the metric cannot see alignment bugs; the TLM tangent comparison can.

**Number check.** Every number in the paragraph matches the md, with three exceptions:
- "0.57–0.72" omits ε = 10⁻⁴, which has 0.84–0.85.
- "∝ ε^{3/2} and independent of σ" is contradicted at 10⁻⁴ by the paragraph's own figures. 3.0·10⁻⁴ against 1.3·10⁻³ is σ-dependent, and the ratio to the 10⁻⁶ values is about 4 100 and 16 000 rather than 1 000.
- Entry shifts at 10⁻⁴ are quantised by the record grid: Δτ = 7.1·10⁻³, and one path shows a shift of 0.

**Novelty.**
- (i) is a one-line corollary of D‴.45(iii)'s dual gradient flow ż = −∇𝒟(z)/N plus (iv)(a). Equivalently, it is preconditioned Langevin with the constant degenerate preconditioner ΦᵀΦ/N.
- A Gibbs law in dual coordinates under additive noise is the construction of mirrored Langevin (Hsieh–Kavis–Rolland–Cevher 2018), which is not cited.
- The contrast with mirror-Langevin diffusions should be stated: Hessian-shaped dual noise gives a primal Gibbs law (Zhang et al. 2020; Chewi et al. 2020; Ahn–Chewi 2021 — from memory).
- The label-noise works on diagonal networks are in the Itô reading, so the draft's framing of them is fair.
- "For every mirror map" fails on branches: ψ* is bounded on one side and π may not be normalisable.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3am-1 | major | (iv)'s window statement, the C paragraph, the 0cf2e9b commit message and the experiment docstring (lines 20–21) say the in-window failure is noise-nonlinear ("one run has the kernel in, the other not", width O(κ/\|u̇\|)). At the tested ε it is linear and finite-ε: the sup error is σ-independent (0.694/0.696/0.699); actual vs TLM is 0.04–0.19 against response error 0.31–0.69; the width is 0.02–0.03 at 10⁻⁶ and 0.002–0.005 at 10⁻⁸ against noise shifts of 7·10⁻⁴–3·10⁻³. The cause is the entering kernel's slow absorption, about 1/(Λ·D_j\|P_{T⁻}φ_j\|²/N). The proof's "rates O(1) in t" fails just after an entry. | D50_draft (iv) and its proof: windows of width O(κ/\|u̇_j\|) + O(1/(Λ·m_j)); numerics_D50 (C) and the exp_transit_noise.py docstring: report the error against distance to the events and its σ-scaling inside the windows |
| D3am-2 | major | (ii) is derived from (i)(b) by letting D_off → 0 at fixed T, but (i)(b)'s Laplace expansion is not uniform in ε when \|S\| < N. Design 3's face Hessian eigenvalues are 2·10⁻¹²–6·10⁻¹⁰, giving Laplace sd of 350–58 000 against Λ = 7 025, so it needs T ≲ ε while B/D have T/ε = 75–7 500. "(ii) proved" is overstated. | D50 (i)(b), (ii), proof, Status: state the exact identity Cov_π(Φw) = T·ΦE_π[D]Φᵀ (holds to 4.6% at Part A's T = 0.125, where T·ΦD(w*)Φᵀ is off by 135%), deduce Cov(w_S) = T·E_π[D_S] + off-support terms, and state (i)(b)'s validity range T ≪ λ_min(ΦDΦᵀ)·(distance to the walls)² |
| D3am-3 | major | (ii) lacks the hypothesis T·D_i ≪ w_i². Design 4 (w₃ = 0.0018) at σ = 0.3: sd_pred/\|w₃\| = 0.59, E[D₃]/D₃ = 0.83, skewness 0.73, ratio 0.44 even over 18 relaxation times. The paragraph blames everything on non-equilibration, but the ratio 1.63 > 1 cannot come from it (Loewner bound from the mode); it is sampling error. | D50 (ii): add the condition and use E_π[D]; numerics B: rerun design 4 for ≥ 30 slow relaxation times at h = 0.5 and compare with T·E[D]; report sd/\|w\| |
| D3am-4 | major | No test covers `tangent_projectors` / `linear_response`: all 8 mutants survive (normal projector, identity, tol, no projection of the accumulated deviation, segment ±1, κ×2). `gibbs_rwm` at 1.1T and 1.2T also survives, so the test pins T only to about ±20%. | tests/test_dual_langevin.py: (a) entries-only path, response = κP_{T(τ)}W(τ) in closed form; (b) response against the scheme's TLM tangent part at small ε (agrees to 10⁻⁵); (c) a Gibbs test using the exact identity Cov(Φw) = TΦE[D]Φᵀ, which is linear in T |
| D3am-5 | major | Part D cannot test (iii). Its window is 2·10⁻⁸–2·10⁻³ of the crossing time Λ²N²/σ². The mean agreement is the scheme's pathwise identity, and all six rows are one Brownian path. (iii) has no numerical support. | numerics D: relabel as the local quasi-equilibrium; to test (iii), sample π directly (Metropolis on F in ν with a face-aware proposal) or use ε with a reachable Λ²N²/σ², and measure E_π[w_off] and the boundary's O(T^{3/2}) |
| D3am-6 | major | Same-path replication overstates precision. B uses seed 100+design at every σ, D seed 7, C the same 4 paths. Independent paths for design 5 give 0.58–1.36 for one coefficient (the paragraph says "1.00–1.05") and max\|corr\| 0.053–0.135. The corr comparison with the design's cosines is not a null. | exp_transit_noise.py b_job/d_job: independent seeds per σ, ≥ 4 paths, report across-path spreads; compare correlations with their sampling error (exactly 0 is predicted when \|S\| = N) |
| D3am-7 | minor | Convergence "from every start" is not covered by Bobkov's Poincaré inequality (needs an L²(π) start). The rate degenerates (O(ε²); mixing Λ²N²/σ²) and is not stated. | D50 proof (i): use inf D_i > 0 ⇒ strong convexity ⇒ W₂/entropy contraction; state the rate |
| D3am-8 | minor | "Run through the mirror update": the scheme is ULA on a potential whose gradient grows exponentially. In the strict sense it is transient and its moments diverge (Roberts–Tweedie 1996; Hutzenthaler–Jentzen–Kloeden 2011 — from memory), so π_h = π + O(h) holds only on practical horizons. | D50 (i): state the law for the SDE and the scheme as an O(h), finite-horizon statement, or use a tamed or implicit step |
| D3am-9 | minor | Scope and novelty overclaims. "Every mirror map" fails on branches (π may be non-normalisable). The paradigm's Euler/Itô update is not covered (Q drifts −2.2 to −3.8). (i) is a corollary of D‴.45(iii)+(iv)(a), and the exact mean fit is already in D‴.45(iv)(a). Mirrored Langevin (Hsieh et al. 2018) is uncited, and the contrast with mirror Langevin is not drawn. | D50 title, "New", Status |
| D3am-10 | minor | Part A calls 0.150 → 0.166 "the O(T) shift of (i)(b)", but the first-order prediction is off by ×1.4–3.7 (coefficient 3: 0.048 vs 0.013) and the sd by −39% to +45%. (i)(b) is never tested in the nonlinear model; the unit test is linear algebra. | numerics A: drop the label; add a small-T run (σ ≤ 0.1) testing the mean shift and T·D^{1/2}PD^{1/2} with P ≠ I |
| D3am-11 | minor | A: the Heun drift is ∝ h√t (ratio 4.06 at t = 2000), a scheme artefact with no stationary law, not "another law". The Itô "0.54–1.09" is one path over one time unit and not a test, and the Itô standard deviations come from a non-stationary window. | numerics A: reword; cite the unit test |
| D3am-12 | minor | C: the outside-window residual's mechanism and rate are not stated; it is the normal OU component, rms κ(Σ1/(2Λλᵢ))^{1/2}, relative size ∝ ε^{1/4}√log. The metric's 2.6–5% floor masks noise-alignment bugs (one-step shift: 0.0263 → 0.0269). | numerics C: state the mechanism and rate; check the response against the TLM tangent part |
| D3am-13 | minor | C: slope 0.267 (CI 0.241–0.291) from three different Brownian path sets (the noise is indexed by step and the τ-grid scales with Λ). Shifts at 10⁻⁴ are grid-quantised (Δτ = 7.1·10⁻³; one path has shift 0). "Order kept" is trivial (gaps 0.65 against shifts ≤ 0.057). | c_job: drive every ε with one Brownian path on a common τ-grid, record finer, report CIs |
| D3am-14 | minor | Numbers: "0.57–0.72" omits 0.84–0.85 at 10⁻⁴; "∝ ε^{3/2}, independent of σ" is contradicted at 10⁻⁴; B's "relaxation times" column counts only the support's normal modes (design 3's face modes have rates O(10⁻¹²–10⁻¹⁰)). | numerics_D50 C, D; transit_noise.md B header |
| D3am-15 | minor | `linear_response` silently drops the first-order jump at exits (designs 4 and 5 have exits). C's windows sit on the ISS events, but at 10⁻⁴ the finite-ε entries are 0.04–0.13 earlier, outside ±0.05. | dual_langevin.py: assert entries-only or add the jump term; C: also report windows at the run's own entries |
| D3am-16 | minor | (iii) wording: "unbiased sparse estimate"; "at every noise level" is ε → 0 at fixed σ and not uniform; the O(T^{3/2}) uses the small-α branch, invalid at σ = 1. | D50 (iii) |
| D3am-17 | cosmetic | The range limit in (ii)'s proof is dimensionally wrong for \|S\| < N; "Cov_π(Φw) = T·ΦDΦᵀ" lacks "+O(T²)"; "entries keep their order" needs κ ≪ the gaps between entries. | D50 (i)(b), (ii) proof, (iv) |

## Lesson

A relative error that does not change with σ cannot be a nonlinearity in the noise. Before naming the cause of a gap between a first-order theory and a simulation, split it with the scheme's own tangent-linear model: actual vs TLM measures the noise nonlinearity, and TLM vs theory measures the finite-ε model error. Size the exclusion windows by each mechanism's predicted width, κ/|u̇| and 1/(Λ·D_j), not by a fixed ±0.05. Prefer exact identities, such as integration by parts E_π[∇F∇Fᵀ] = T·E_π[∇²F], to leading-order expansions whose uniformity fails in exactly the regime of interest. Finally, do not present one noise path scaled across σ, or a window a millionth of the mixing time, as replication or as evidence about a stationary law, and test new functions directly rather than relying on the experiment to catch their mutations.
