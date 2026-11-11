# Referee D3ay, round 52: Proposition D‴.50 (the rich transit under label noise, as revised after round 40)

## 1. What I checked and how

**Scope and setup**
- I reviewed `review/drafts/D50.md` as revised at bc14fdb, the code at 340a61e (`experiments/exp_transit_noise.py`, `src/core/dual_langevin.py` and the `mirror_free`, `rich_transit`, `exp_mirror_transit` and `exp_cohesion_leak` helpers it calls), `results/chaos/transit_noise.md` and its raw npz, `tests/test_dual_langevin.py` (the only test that imports this code; no test imports `exp_transit_noise`), and `review/REVIEW_round40_D3am.md` with its `REVIEW_LOG.md` entry.
- I modified nothing in the repository and ran only read-only git commands.
- During the review other sessions committed to the repo (HEAD moved from 95c3afe to 265a133), and an untracked `review/REVIEW_round53_D3az.md` appeared; neither is mine. `git diff 95c3afe HEAD` is empty for every D‴.50 file.
- All my work is in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref52/`:
  - `repo/` is a copy of `src`, `experiments`, `tests` and the two result files;
  - `work/` holds my scripts and logs;
  - `mut/` and `mutpost/` hold the mutant trees.
- Every Python run used `OMP_NUM_THREADS=1`, `nice -n 15`, `PYTHONDONTWRITEBYTECODE=1` and `-B`, with at most two of my processes at once. Each check stayed under 10 CPU-minutes.
- The harness stopped one background mutation runner at its 30-minute limit, after mutant D17. I reran D18 on its own.

**1. Analysis.** I re-derived by hand:
- Q's conservation and its Itô drift. My result, ½(σ²|φ_i|²/N²)βγ[γa² − Gβa′(1 + a²/a_max²)], equals `ito_q_drift`.
- ∇F = Φw − y, the Langevin form and T = σ²/(2N).
- Normalisability on free curves: ψ*(η) grows like exp(4Gc|η|/a_max) and Φᵀ is injective.
- inf D_i > 0, and the W₂ rate by synchronous coupling.
- The integration-by-parts identity.
- Laplace's Cov(w) = T·DΦᵀ(ΦDΦᵀ)⁻¹ΦD = T·D^{1/2}PD^{1/2}.
- The R_off algebra of (ii), term by term.
- Both critical-curve asymptotes.
- (iii)'s epi-limit, near-uniform face law and the O(T^{3/2}) penetration, which is 0.49·T^{3/2} on the instance.
- In (iv):
  - κ = σ/(NΛ^{1/2});
  - the inverse scale space as the normal Skorokhod problem (by Moreau, the NNLS residual is the projection onto the tangent cone);
  - the entry jump δν⁺ = δν⁻ + (v⁻ − v⁺)δτ = P_{T⁺}δν⁻;
  - the nested composition κP_{T(τ)}W(τ) and the exits' jump;
  - the OU variance κ²tr(M⁺)/(2Λ);
  - the absorption rate ΛD_j|P_{T⁻}φ_j|²/N.

**2. Numbers.**
- I ran `--from-raw` on a copy (`TN_OUT` set to scratch). The regenerated md is byte-identical to the repo's.
- I checked every number in the draft's Numerics against the md, and inspected the raw per path.

**3. Fixes.** I checked each of the 17 round-40 fixes in the code, the md and the draft.

**4. Tests.** The scratch copy gives 9 passed in 249 s.
- I built 25 single-replacement mutants: 19 in `dual_langevin.py` and 6 in `exp_transit_noise.py`.
- Each mutant ran the suite fast-first with `-x`; survivors ran the full suite.

**5. Targeted checks.**
- **(a) Part C replica at ε = 10⁻⁶.** I re-ran the noiseless run plus 3 σ × 4 paths with the repo's own `c_job`. The md's three ε = 10⁻⁶ rows and every per-path raw value are reproduced exactly.
- **(b) Part B replica** for design 4, σ = 0.3, path 0. Its ratios and E[D]/D(w_BP) equal the raw exactly.
- **(c) Window analysis.** On design 3's ISS path I computed the windows' coverage and the entry-shift standard deviation.
- **(d) Time-step bias of Part B.** I computed the linearised unadjusted-Langevin (ULA) stationary covariance at B's steps.
- **(e) Part A's mean fit**, reconstructed from the run's seed.
- **(f) Heun drift ensembles:** 12 seeds at h = 0.01 and 8 seeds at h = 0.02, each to t = 2000.
- **(g) Mutants' effect on the md**, using the replicas and the raw.

**6. Literature.** I used web search only. WebFetch could not reach arxiv.org, ar5iv, projecteuclid.org or sufe.edu.cn (hosts not resolvable in this session), so the citations below rest on search summaries or abstracts.

## 2. Verdicts

### (i) The Gibbs law and the exact identity
**Correct.**
- The ζ-equation is overdamped Langevin at T = σ²/(2N). The noise is additive in ζ, so the mirror update is a consistent Euler–Maruyama scheme of the Stratonovich (α, β) dynamics.
- π is normalisable on free curves: F grows superlinearly because each ψ*_i grows exponentially.
- The W₂ rate min_i inf D_i·λ_min(ΦΦᵀ)/N is right. inf D > 0 holds, although for β⁰ > 4c the minimum of D is not at α = 0.
- E_π[∇F∇Fᵀ] = T·E_π[∇²F] is exact; the boundary terms decay doubly exponentially.
- (c) is right as T → 0 asymptotics. Its validity range T ≪ λ_min·δ² is a heuristic.
- **Unproved:** "the scheme samples π up to O(h)" is not covered by the proof. Its bias is visible in B (D3ay-7). In A, the "O(h)" label on the mean fit is wrong (D3ay-5).
- **Measured in A:**
  - the identity's error is 0.0213 and 0.0218 over the whole run, against 0.086 and 0.088 per batch, which is √20 times larger: the sampling floor;
  - Laplace's value is off by 0.554 and 0.558;
  - Cov(w) agrees with T·D^{1/2}PD^{1/2} to 0.016 at T = 1.25·10⁻³.

### (ii) Fluctuation–dissipation on a support
**Correct as an exact consequence of (i)(b).** I checked the R_off formula and the asymptotes 4Gc|w|/a_max and 12Gc²α⁴/(ma_max⁴).
- B's Var/(T·E[D]) values (0.93–1.12) are consistent with 1, except where they show the scheme's predicted time-step bias. My linearised ULA factor Φ_S⁻¹(I − hA/2)⁻¹Φ_S predicts:
  - **design 4 (h = 0.5):** 1.0249 and 1.0185 for its two large coefficients; measured 1.024–1.026 and 1.019–1.020 at all three σ, which is 5–10 se from 1 and within 0.5 se of the prediction;
  - designs 3 and 5: at most 1.006 and 1.009.
- The predicted time-step correlations are at most 0.0066, so the measured |corr| of 0.034–0.146 are sampling noise.
- The draft's evidence that "D at w_BP misses" for design 4's third coefficient does not discriminate between the two predictions (D3ay-8).
- For |S| < N, R_off under π is untested (D3ay-10).

### (iii) The rich limit
**Formally consistent and honestly labelled untested.** Inside the face ψ* = O(ε), so the excess over the temperature T/Λ is O(ε/T) → 0, which gives a near-uniform law on the face. The small-α condition T ≪ 2c/(γm) = 0.042 indeed fails at σ = 1 (T = 0.083).

### (iv) The noisy transit
**The first-order response is correct and well confirmed.**
- I re-derived δν = κP_{T(τ)}W(τ) with the projection at each entry.
- The first-order entry shifts −φ_jᵀδν_lin(τ_e)/u̇_j match the measured ones. For example, at σ = 1, path 1: (−0.0039, +0.0104, −0.0260) predicted against (−0.0044, +0.0101, −0.0252) measured.
- The OU residual ratio is 1.005–1.010 at 10⁻⁶ and 10⁻⁸; my replica gives 0.972–1.038 per path.
- Share/κ is ε-independent per path to within 10% (for example 6.385 → 6.444 and 8.20 → 7.43), and the slope is 0.252 ± 0.006.

**Not right as stated:**
- the shift-window width O(κ/|u̇_j|), which is missing |P_{T⁻}φ_j|τ_e^{1/2} (D3ay-2);
- "σ-independent" inside the windows (D3ay-3);
- the 10⁻⁴ "outside" numbers, which test no part of the transit (D3ay-9).

### Numbers
- Every number in the Numerics paragraph matches the md, and the md matches the raw.
- The two hand-derived numbers are right: ε^{0.20} (ln(0.1675/0.026)/ln 10⁴ = 0.202) and "≥ 30 relaxation times".
- Four labels are wrong:
  - "(the scheme's O(h))" (D3ay-5);
  - "± batch errors" (D3ay-6);
  - "∝ h√t" (D3ay-4);
  - "the windows' widths as predicted", which quotes the formulas' own values (D3ay-2).
- One comparison is unlike-for-like: the sampler's 2.8·10⁻³ against the scheme's 10⁻² (D3ay-5).

### Round-40 fixes

| D3am | Status |
|---|---|
| 1 | Present, partly wrong: width formula (D3ay-2), σ-independence (D3ay-3), 10⁻⁴ coverage (D3ay-9) |
| 2 | Present and correct; its ± is mislabelled (D3ay-6) |
| 3 | Present; the numerical E[D]-versus-D(w_BP) claim is unresolved (D3ay-8) |
| 4 | Partial. The identity test now kills 1.05T and 1.1T (rel 0.064 and 0.093 against 0.04). The closed-form test is circular; `normal_ou_var` is untested; round 40's tangent-linear-model test is not implemented (D3ay-1) |
| 5 | Correct |
| 6 | Correct for B and D (independent seeds); D's ε rows share seeds (identical z columns) |
| 7 | Correct |
| 8 | Stated but unquantified, and listed as "proved" (D3ay-7) |
| 9 | Present; further prior art is missing (D3ay-11) |
| 10 | Correct; md and code still say "(i)(b)" (D3ay-12) |
| 11 | The Itô part is correct; "∝ h√t" is wrong (D3ay-4) |
| 12 | Correct |
| 13 | Correct: the coarse increments are sums of the finest's (r = 10, 100) |
| 14 | Correct |
| 15 | Correct (mutant D19 killed) |
| 16 | Correct |
| 17 | Correct |

### Tests and mutants
**15 of 25 killed** (15 of 19 in the library, 0 of 6 in the script).

**Killed:**
- D1 (noise ×√1.02) and D2 (drift h/(N+1)): the integrator test.
- D3 (dw without m): the FastFree test.
- D4 (ψ*/1.1) and D5 (+yᵀζ): the gradient test.
- D6 (1.05T), D7 (1.1T) and D8 (proposal recorded): the Gibbs-identity test.
- D9 (D for D^{1/2}): the `fd_covariance` test.
- D10 (Itô sign): the Itô test.
- D11 (Heun → Euler): the Heun test.
- D13 (normal projector), D15 (interval-start lookup), D16 (no projection of the accumulated deviation) and D19 (exit guard off): the response test.

**Survivors, all of which change the md:**
- D12 (identity tangent projector): sup 0.708 → 1.482, outside the windows 0.066 → 1.455, OU ratio 1.005 → 35.1.
- D14 (drop the lowest-index active atom): 1.276, 1.276 and 24.1.
- D17 (OU without the factor 2): 1.005 → 0.710.
- D18 (OU without /N): 1.005 → 2.461.
- E1 (shift width κ/|u̇|): widths 1.47·10⁻³ → 4.91·10⁻⁴.
- E2 (P_{T⁺} in the absorption width): widths 1.9·10³⁰, mechanism columns NaN.
- E3 (response blocks without the 1/√Λ factor): sup 84.
- E4 (D at the mean α): design 4's third coefficient 1.018 → 1.217, E[D]/D(w_BP) 0.7315 → 0.6117.
- E5 (OU ratio over the ×1 windows): 1.005 → 1.003.
- E6 (batch se divided by √(nb − 1)): 26 md lines change, for example D's z 1.70 → 1.66.

All six script mutants applied together still give 9 passed.

### Literature
- **Classical:**
  - the Gibbs law of preconditioned or mirrored Langevin, and the dual Gibbs law of stochastic mirror descent with additive dual noise when ∇f∘∇Φ* is a gradient (Borovykh–Kantas–Parpas–Pavliotis, Physica D 2021; model of Raginsky–Bouvrie 2012);
  - (i)(b), which is the configurational-temperature identity in matrix form (Rugh 1997; Butler–Ayton–Jepps–Evans), with Yaida's fluctuation–dissipation relations for SGD (ICLR 2019) as the discrete analogue;
  - (iv)'s first-order response, which is the directional derivative of the Skorokhod map in a convex polyhedron (Mandelbaum–Ramanan 2010; Lipshutz–Ramanan 2018, 2019) — the functional CLT of a reflected process around its fluid limit;
  - the Itô reading's drift of conserved quantities, which is the law of balance (Ziyin–Li–Ueda).
- **Closest contemporary:** Corlouer–Semler–Strang–Gietelink Oldenziel (arXiv 2604.06366, 2026). They find an approximately Boltzmann stationary law for label-noise SGD in the saddle-to-saddle regime, and that the noise does not alter the saddle-to-saddle path.
- **Nothing found on:** noisy inverse scale space flows; the coefficient-level identity Cov(w_S) = T·E[D_S] on a support; weighted basis pursuit as a stationary mean; the finite-ε absorption window; the rate σε^{1/4}. These appear new.

## 3. Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ay-1 | major | **The tests cannot see the projector or the OU code.** `test_linear_response…` builds its closed form from `DL.tangent_projectors` itself (`tests/test_dual_langevin.py:177–180`), so any nested family of projectors passes, including the identity (D12) and dropping the lowest-index atom (D14; the test path's entries are 5, 0, 4). `normal_ou_var` (`dual_langevin.py:240`) is never called by a test, so D17 and D18 pass; it produces the headline OU ratio. No test touches `exp_transit_noise`, so E1–E6 pass. All 10 survivors change the md (numbers in §2). | `tests/test_dual_langevin.py`: test P against independent properties (P² = P = Pᵀ, PΦ_A = 0, rank N − \|A\|) on a three-entry path whose entry order is not the index order (design 3: 8, 13, 7); test `normal_ou_var` against `solve_continuous_lyapunov` for the block OU; implement round 40's tangent-linear-model comparison; add a tiny smoke test of `windows_C` and the `part_C` metrics |
| D3ay-2 | major | **(iv)'s shift width O(κ/\|u̇_j\|) is dimensionally incomplete** (`D50.md:52`, code `exp_transit_noise.py:292`). δτ_j = −κφ_jᵀP_{T⁻}W(τ_e)/u̇_j has sd κ\|P_{T⁻}φ_j\|τ_e^{1/2}/\|u̇_j\|. The missing factor is 3.07, 2.06 and 2.74 for design 3's three entries, so 3κ/\|u̇_j\| is only about 1–1.5 sd. The md's measured maximum shifts are 1.49–1.77 times the ×1 width in all nine cells. The ×1 windows fail at 10⁻⁸, σ = 1 (0.151 and 0.130 on paths 1–2, against 0.026 at ×2). At 10⁻⁶, σ = 1 the real shifts reach 0.025, which is 0.7 of the absorption width, not "dominated" by it (`D50.md:131`). "The windows' widths as predicted" (`:132`) quotes the formulas' own values; the Status's "predicted and checked" (`:90`) overstates this. As ε → 0 the shift window (∝ σε^{1/4}) dominates the absorption window (∝ ε^{1/2}), so this is the regime that matters. | (iv), its proof and `windows_C`: width κ\|P_{T⁻}φ_j\|τ_e^{1/2}/\|u̇_j\|. md: report measured against predicted widths. Store the per-run C arrays in the raw so the windows can be recomputed without a rerun |
| D3ay-3 | minor | **"Inside the windows … σ-independent too"** (`D50.md:131`; md header `exp_transit_noise.py:467`; docstring line 25) is contradicted. Each Brownian path drives every σ, so a linear error is exactly σ-independent per path. Per path at σ = 0.1 → 1, the in-window error falls by up to 39% at 10⁻⁶ (0.560 → 0.340) and by 34% at 10⁻⁸ already at σ = 0.3 (0.517 → 0.339). The md maxima are 0.708/0.698/0.559 and 0.717/0.620/0.572; L² is 0.118 → 0.080. Removing the interval between the two runs' entries leaves the drop unchanged (0.397, 0.330, 0.559), so it is nonlinear absorption once the real shift is comparable to the absorption width. σ-independence holds at 10⁻⁴ (≤ 6%) and at 10⁻⁶ for σ ≤ 0.3 (≤ 5%). | Draft C, md header, docstring: "σ-independent where the shift is ≪ the absorption width"; report the per-path σ-ratio |
| D3ay-4 | minor | **The Heun drift along ker Φ is ballistic (∝ h·t), not ∝ h√t** (`D50.md:114`; md text `exp_transit_noise.py:420`; docstring lines 14–15). Over 12 seeds at h = 0.01, \|E P_{ker Φ}η(2000)\| = 0.778 ± 0.064; \|E\|/(ht) = 0.039 for t from 500 to 2000; the rms log-slope is 0.96. Over 8 seeds at h = 0.02: 1.419 ± 0.101, and \|E\|/(ht) = 0.036. Only the random part is ≈ 0.5·h√t. The draft's own single paths grow ×14 and ×17 from t = 50 to 2000, against ×6.3 for √t. "No stationary law at a fixed step" still holds. | Say: a weak-order-1 bias v̄·h·t (\|v̄\| ≈ 0.04 on A's dictionary) plus O(h√t) |
| D3ay-5 | minor | **"E[Φw] = y to … 1.0·10⁻², 8.8·10⁻³ (the scheme's O(h))"** (`D50.md:107`) is wrong. The values are the noise term σ√hΣξ_k/(t − t₀) of (i)'s own pathwise identity. Reconstructed from seed 2: (−0.0023, 0.0091, −0.0012, 0.0106) against the stored (−0.0024, 0.0096, −0.0013, 0.0104), and similarly at h = 0.005; the sd is 0.0075 per component. ULA has E_{π_h}[∇F] = 0 exactly, so there is no O(h) in the mean fit. The comparison with the sampler's 2.8·10⁻³ (10⁶ steps) is unlike-for-like. | Draft A: "the pathwise identity's σW(t)/t, 0.0075 per component" |
| D3ay-6 | minor | **"0.021 ± 0.003 and 0.022 ± 0.005 (batch errors)"** (`D50.md:109`; md; code `:137`, `:411`). The ± is the standard error of the mean of the 20 per-batch errors (0.086 ± 0.003, 0.088 ± 0.005), not of 0.021. Since 0.0213·√20 = 0.095 ≈ 0.086, 0.021 is the sampling floor: the identity holds to sampling accuracy. | Report it as "at the sampling floor (per batch 0.086, scaling as length^{−1/2})" |
| D3ay-7 | minor | **B carries the scheme's O(h) bias, unmentioned.** For design 4 at h = 0.5 the linearised factor diag Φ_S⁻¹(I − hA/2)⁻¹Φ_S = (1.0249, …, 1.0185, …) matches the measured 1.024–1.026 and 1.019–1.020: 12 of 12 path values lie in [1.018, 1.032] and [1.009, 1.027]. Separately, the Status's "(i) proved" (`D50.md:87`) includes two heuristics: the scheme's "π up to O(h)" (`:18`) and (c)'s validity range (`:26`). | (i) and Status: mark both as heuristic and give the linear-response factor; B: quote the predicted bias or use `H_MAX_B` ≤ 0.1 |
| D3ay-8 | minor | **Design 4's third coefficient (σ = 0.3): "D at w_BP misses" is unresolved** (`D50.md:119–120`; commit bc14fdb). Per path: E[D]/D(w_BP) = 0.73, 0.86, 0.80, 1.19 (0.894 ± 0.10 se); Var/(T·D(w_BP)) = 0.745, 0.726, 1.011, 1.100 (0.90 ± 0.10); Var/(T·E[D]) = 1.02, 0.85, 1.27, 0.92 (1.01 ± 0.09). Both predictions are within about 1 se. Thirty relaxation times give a per-path relative sd of about (2/30)^{1/2} = 0.26. Also, the batches on designs 4–5 are 1.5 slowest relaxation times long, not "a few" (`:121`). The lag-1 correlation of OU batch means is about 0.42, so the se are low by about 1.4×, which inflates the z of 4.1 and 4.8. | Draft B: "not resolved at this precision"; use ≥ 300 relaxation times or more paths; batch on the slowest mode, or quote paths only |
| D3ay-9 | minor | **At 10⁻⁴ the windows cover the whole transit.** The ×2 windows leave 26.7% of the τ-range outside but 0.0% of the transit (τ₁, τ*); the ×4 windows cover everything (NaN). So the 10⁻⁴ figures "0.167–0.168" and "OU 0.939" (per path 0.745–1.214) test only the stretch before the first entry, where there is no active block, and the tail after the transit. At 10⁻⁶ and 10⁻⁸, 36% and 43% of the range is transit outside ×2. The ε^{0.20} (`:128`) mixes these regions and departs from (Λλ)^{−1/2} ∝ ε^{1/4} through the √log of a sup. | md: report the fraction outside each window; draft: state the 10⁻⁴ caveat and the √log |
| D3ay-10 | minor | **R_off for \|S\| < N is untested under π.** The Status says "R_off's size numerical" (`:88`). Under π the dual covers the whole face, including its boundary, where D already finds w_off up to 5.4·10⁻² (10⁻⁴, σ = 1). B's design-3 windows (≈ 2·10⁻⁶ of the crossing time) cannot see this; only designs 4–5 (\|S\| = N) sample π. | Status (ii): measured for \|S\| = N; for \|S\| < N only locally |
| D3ay-11 | minor | **Prior art missing:** Skorokhod-map directional derivatives (Mandelbaum–Ramanan 2010; Lipshutz–Ramanan 2018, 2019) for (iv)'s κP_TW; stochastic mirror descent's dual Gibbs law (Borovykh et al. 2021; Raginsky–Bouvrie 2012); the configurational temperature (Rugh; Butler et al.; already cited in the repo's D‴.40) and Yaida's fluctuation–dissipation relations for SGD; the law of balance (Ziyin–Li–Ueda; cited in D‴.41); Corlouer et al. 2026. | D50 classical list and "New"; `LITERATURE_TRANSPORT.md` |
| D3ay-12 | cosmetic | **Labels and wording:** | |
| | | – md and code say "(i)(b) at small T" (`exp_transit_noise.py:16, 57, 139, 415`); it is now (i)(c). | |
| | | – "γa_max = 4Gc/a_max on the paradigm's instance" holds by the definition of c on every instance. | |
| | | – The three per-σ slopes (`:124–125`) are one measurement (common paths, linear response). | |
| | | – "One reaches 5.4·10⁻² … 70%" (`:138`) is one path of four; the others give 0.4%, 11% and 16%. | |
| | | – "Independent of σ" at 10⁻⁶ compares 7.4·10⁻⁸ with 1.1·10⁻⁷. | |
| | | – "∝ σt^{−1/2}" is asserted from one t. | |

## 4. Lesson

A closed-form test built from the function's own ingredient cannot catch an error in that ingredient. Here the response test accepts every nested family of projectors, the identity included. Test each building block against properties it must satisfy independently: idempotence, the annihilated columns, a Lyapunov solve. Dimensional analysis would have caught the window width: κ/|u̇| has units of τ^{1/2}, so something with those units was missing, and the md's own column of measured shifts — 1.5–1.8 times the "width" in every cell — said so. When the noise paths are common across σ, any σ-dependence of a relative error is nonlinearity and can be measured path by path; report it that way rather than calling a maximum over paths σ-independent. Name an error term by testing it: a mean fit reproduced exactly from the seed's noise is not O(h). And settle a growth law (√t against t) with an ensemble of a dozen seeds, which shows the mean drift; one or two paths cannot.
