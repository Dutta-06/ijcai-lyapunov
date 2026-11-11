# Referee report D3ep, round 147: Proposition D‴.40 as revised after round 31

## Scope
This round reviews Proposition D‴.40 ("the gauge under general noise"). The text is at PARADIGM.md ll. 4953–5014. Its summary row is l. 152, its numerics row l. 9182 and its roadmap entry ll. 9549–9557. The commit is 1510c79, in the state left after the round-31 revision.

The numerics come from `experiments/exp_gauge_noise.py`, which writes `results/chaos/gauge_noise.md` and its raw file `gauge_noise_raw.npz`. Supporting code is `src/core/gauge_noise.py`, plus the `noise_L`, `mask`, `fixman_*` and `box` paths of `src/core/mass_transfer.py`. The tests are `tests/test_gauge_noise.py` and `tests/test_runcache.py`.

D‴.40 is checked against the round-31 issues (D3ad-1…9 and 13) and against the dated correction blocks after rounds 143 (l. 4928) and 142 (l. 5273). I did not modify the repository. All work was done in an export at …/scratchpad/ref147/repo, with scripts in …/scratchpad/ref147/work.

## Summary
**Design.** One referee covered four things:
- **Re-derivation** of every closed form, by hand and with numeric scripts.
- **The md reproduced from its raw file**, plus the tests.
- **New computations:**
  - Katzenberger's drift at finite γ for Part E's noise.
  - The zero of the corrected drift.
  - The leading-order ODE set against Part E's time series.
  - The coloured-noise linear response, with a Monte Carlo on the linearisation.
  - A nonlinear Monte Carlo on the instance: fixed masses, the box, OU-coloured against white amplitude noise, 6·10⁵ slow units each.
- **Literature:** the Shalova–Schlichting–Peletier (SSP) and Li–Wang–Arora (LWA) abstracts read on arXiv, and SSP's theorems read in full.

**Verdict.** D‴.40 as revised is correct within its stated hypotheses, which are white noise with constant covariance in the unconstrained chart ϑ.
- (i)'s generator identity is exact, and the round-31 repairs hold.
- (ii) is solid. Its zero-drift and freezing statements extend to any noise without an amplitude component, coloured or not.
- (iii) is right at leading order but drops two hypotheses.

What is wrong is the scope and the statuses, not the algebra:
- **Scope.** The title and summary say "general noise", but the result is chart-bound and white-only.
- **Leading order in the boxed law.** It is a formal small-noise expansion of a non-Gibbs law, not "proved".
- **Slow-limit flow.** It inherits D‴.37(ii)'s sketched averaging step.
- **Part E.** Its uncorrected and κ_F = 2 runs are reported without the predictions they should be tested against.

Revision needed: minor-to-major.

**Findings**
- **D3ep-1 (major): chart and noise class.** The statement holds for white noise of constant covariance in ϑ. It fails beyond that:
  - **Coloured noise.** The leading order becomes −g′_k[(I+τ_cGH)⁻¹Σ]_{α_kα_k}/(2G_α). On the instance at τ_c = 100 this is 0.47, 0.37 and 0.28 of the white value. Correlations enter at O(τ_cλ): 8% at τ_c = 100 and 28% at τ_c = 1000 for 'correlated E'. A nonlinear Monte Carlo confirms the reduction: E[u_k] − box is (−0.56 ± 0.44, −0.73 ± 0.16, −0.15 ± 0.21)·10⁻³ coloured against (−2.53 ± 0.46, −1.87 ± 0.16, −1.72 ± 0.22)·10⁻³ white, same seed.
  - **Quadratic-variation reading.** Under OU noise α is C¹, so the level read from the quadratic variation in (iii) is wrong: it tends to 0 as Δt → 0.
  - **Stratonovich multiplicative noise.** It adds −g_k c_α/G_α at O(τ), where c_α involves the derivatives of every noise coefficient.
  - **Noise white in the physical amplitude a.** The same identity gives ū_k = −τ_k/β_k. The selection is then Σ_kτ_k log W_k, with no divergence at the saturation, plus a reflecting-boundary term. The statement does not mention the chart, although round 143 (2) corrects exactly this for D‴.37.
- **D3ep-2 (major): status of the leading order in the boxed law.** "The leading order proved" is formal there. For non-isotropic Σ the stationary law is not Gibbs. The plateau's weight is governed by a quasi-potential that is not computed. At the instance's τ, D‴.39(iii)(b) itself says the box's term is as large as the noise's selection.
- **D3ep-3 (major): status of the slow-limit flow.** "In the slow-replicator limit the masses follow the replicator flow of ū … whatever Σ" rests on D‴.37(ii)'s averaging, which is sketched (Paper 14 calls it a "Derived statement"). It also needs:
  - compact sets of M (round 143 (3));
  - times short against the exit from the well, whose margin is only about ×7 at τ = 10⁻⁴ (round 143 (5));
  - a large Péclet number (round 143 (7));
  - for degenerate Σ, unique ergodicity, which is not given by ellipticity.
- **D3ep-4 (major): Part E is not confronted with predictions.** The uncorrected 'correlated E' run:
  - The leading-order ODE sends W_1 to the box wall 0.2947 by 2.2·10⁴ slow units, then toward the vertex (0.295, 0.206, 0.499).
  - Katzenberger's drift at γ = 10⁻² is (−0.834, +0.318, +0.515)γτ, 65% of the slow limit.
  - The 12-seed mean of W_1 instead rises to 0.343 ± 0.007 at 2·10⁴ (about 3.9 s.e. above about 0.316) and never falls below 0.326.
  - W_2 falls (0.246 at 10⁶) where the leading order and Katzenberger's drift give a rise at 1/3.
  - The late mean drifts toward that vertex; seed spread 0.03–0.04, so the Péclet number is about 1.
  - The positions' spreading transient is only 2·10⁻⁴. The early discrepancy is therefore unexplained. My 36-seed rerun to settle it was stopped after one seed when the report was called for. That seed gave W = (0.300, 0.337, 0.363) at 2·10⁴, consistent with the prediction.

  The κ_F = 2 run:
  - Its second-half mean is 0.022 from W♯_τ (4.4–5.8 s.e.).
  - The zero of the corrected Katzenberger drift at γ = 10⁻² is (0.4508, 0.2821, 0.2671). That is already 0.010 from W♯_τ, a finite-γ shift like D‴.38(i)'s over-compensation, which the text does not mention.
  - At γ = 10⁻³ and 10⁻⁵ the zero is within 3·10⁻⁴ of W♯_τ.
- **D3ep-5 (minor/major): two dropped hypotheses in (iii).**
  - Every τ_k must be > 0. If one class has no amplitude noise, F_τ is flat in its mass and the κ_F > 1 maximiser sits on that class's saturation, not in the interior.
  - With the box, W♯_τ must lie in P_B. For levels (2, 1, 0.2) it is (0.4813, 0.2926, 0.2262), below W_B,3 = 0.2358; for (1, 1, 0.1) likewise. This is the analogue of D‴.39(iii)(d).
  - The helper `w_sharp_weighted` raises a brentq sign error at levels (2, 1, 0.2).
- **D3ep-6 (minor): unused sign statement.** g_k∂_{α_k}V ≥ 0 pointwise, because sign g_k = sign α_k = sign ∂_{α_k}V outside B. The box's term therefore always raises the class's fitness and opposes the noise's drive into the saturation. This is provable and matches all of Part D, but is not stated.
- **D3ep-7 (minor): (ii) can be strengthened.** Without amplitude noise, u_k = dψ_k/ds + g_k∂_{α_k}(J_coh + V) pathwise. The masses' increments therefore telescope, and the freezing holds for any noise process with no amplitude component (coloured, non-Gaussian, multiplicative), not only in Katzenberger's limit.
- **D3ep-8 (minor): proof gaps.**
  - The Gaussian tails of the non-Gibbs boxed law are asserted without an argument. A Lyapunov function exp(ε|ϑ|²) supplies one.
  - "½E[g‴_k]C" should read g‴_k at the fit.
  - "log-barrier" in (ii) for the floored correction (round 143 (6)).
  - "weighted log-determinant" in rows 152 and 9549 for F_τ; it is a log-determinant only for the α-chart and equal levels.
- **D3ep-9 (minor): literature.** Missing are noisy replicator dynamics (Fudenberg–Harris's Itô correction from own payoff noise, the exact analogue of "each class's own level"; Imhof; Hofbauer–Imhof) and SSP's result that noise structure fixes the regulariser (∇²L:C) and the time scale.

## What I checked and how
- **[read]** D‴.40's statement, proof, Status and Numerics, rows 152, 9182 and 9549. Also D‴.36–D‴.39 (statements and statuses), the round-143 and round-142 correction blocks, REVIEW_LOG round 31 and REVIEW_round31_D3ad.md (all D‴.40 parts). In the code: `gauge_noise.py`, `exp_gauge_noise.py`, `run_dithered`, `fixman_correction`, `limit_drift`. In Paper 14: README items on D‴.40, §9 (1), the dithered-chart remark, Derived statement der:slow.
- **[own] Closed forms.** I checked numerically that:
  - ψ′ = g/G_α and g′ = (1+q)/(βa′) = cosh(2α/a_max)/β, to 5·10⁻¹¹;
  - ū_k = −∂_{β_k}F_τ, to 10⁻⁸;
  - log(W − c/W) is strictly concave (h″ < 0; the numerator −W⁴ − 4cW² + c² < 0);
  - W♯_τ by direct maximisation is (0.46088, 0.28238, 0.25675), equal to the helper's value.

  Also [own] (hand derivations): the generator identity term by term, D = 0 from the block structure, the τ̂ bias, and the coloured, Stratonovich and a-chart variants of the Stein identity.
- **[own] Dynamics.**
  - The leading-order replicator ODE for Part E's four variants, with a wall event (`ode_E.py`).
  - Katzenberger's limit drift at γ = 10⁻², 10⁻³ for 'correlated E', its diagonal part, amplitude levels only, isotropic and amplitudes only. Correlations change the drift by ≤ 1% and D by ≤ 3% at finite γ.
  - The zero of the corrected Katzenberger drift (κ_F = 2, levels τ_k) at γ = 10⁻², 10⁻³, 10⁻⁵ (`zero_corrE.py`).
  - Part E's per-seed records (`peekE*.py`).
  - The positions' spreading transient (`transient.py`).
- **[own] Coloured noise.**
  - Linear response on the instance's GH (rates 3.06·10⁻⁴ … 0.195).
  - A Monte Carlo on the linearisation at τ_c = 100: E[δα η_α]/(Gτ) = 0.903, 0.394, 0.147 against the formula's 0.923, 0.397, 0.154.
  - A nonlinear Monte Carlo on the instance (`mc_colour.py`): fixed masses 1/3, the box, amplitudes-only noise at τ = 10⁻³, OU τ_c = 100 against white, same seed, 6·10⁵ slow units. The exact coloured Stein identity E[u] = −E[gη_α]/G_α + box holds to (−2 ± 27, −1.5 ± 4, −0.6 ± 4.5)·10⁻⁶.
- **[md]/[raw] Reproduction.** I copied the raw file into a directory set by GN_OUT and ran `--from-raw`. The regenerated md is byte-identical to the committed one (cmp). Tests: 9/9 pass (7 in test_gauge_noise, 2 in test_runcache).
- **Decimals.** Every decimal in D‴.40's text and rows is printed in the md: 1.3·10⁻⁸–1.0·10⁻⁶, 5.3·10⁻¹¹, 1.46–5.77, 1.1–8.8·10⁻⁴, 8·10⁻¹⁶, 10⁻³⁶, 3.2·10⁻⁴…8.0·10⁻⁷, 0.03–0.93, 1.1·10⁻²², 1.6–13, 58, 3.9/3.7 ± 0.3·10⁻², 10⁻⁷, 0.2–1.9%, χ² = 24.8, |z| 1.71, 0.03 errors, 7·10⁻¹⁹ (6.9e-19), 1.4%/0.25%, 2·10⁻⁴, 3–10·10⁻⁴, the W♯ triples, 0.326/0.246/0.428, 0.10, 0.03–0.04, 0.409/0.428/0.440, 0.022, 0.062, and 9·10⁻⁴ / 10⁻⁶ in the rows. One reading needs care: the τ_k/τ̂_k column is (with, without) pairs per class, and only that reading gives 0.25%. Recommend labelling the column.
- **[paper]** I read SSP arXiv:2404.12293 (theorems A and B, Assumption 3.5, §6.2 on correlated noise) and the LWA abstract (arXiv:2110.06914), both via pdftotext with python -I.
- **[memory]** Katzenberger 1991; Funaki 1995; Fudenberg–Harris 1992; Imhof 2005; Hofbauer–Imhof 2009; Rugh; Jepps–Ayton–Evans.

## Verdicts per part
- **(i), the generator identity.** Correct and exact for the diffusion in any stationary law. The Euler chain has an exactly measured O(Δt) bias. The well form is right up to exit terms, which are beyond every order. The leading-order formula and the F_τ reading are correct algebra. Its standing as the small-τ asymptotics of the boxed average is formal (D3ep-2). "Only own intensity" holds for white additive noise in ϑ, even at finite γ (≤ 1% in Katzenberger's drift). It fails for colour, Stratonovich multiplicative noise and other charts (D3ep-1).
- **Signs.** ū_k < 0 in every class, and the selection runs toward the saturation in the α-chart. The box term is ≥ 0. The replicator's sign depends on the differences only. All consistent with Parts A–D.
- **Orders.** The O(τ²) term is correctly attributed to all of Σ. Part B's 1.1–8.8·10⁻⁴ is the O(γ) correction (coefficient about 30–90, consistent with D‴.38(i)'s 57% at γ = 10⁻²), not an error.
- **(ii).** Correct; D = 0 holds at every γ by the block structure. Stronger than stated (D3ep-7). Part C confirms the residual drift is O(γ²τ).
- **(iii).** Correct at leading order with all τ_k > 0 and W♯_τ inside the representable set, or inside P_B with the box (D3ep-5). The plug-in cancellation is exact for the noise part in the slow limit (floor inactive). The quadratic-variation reading of τ_k is valid for white noise only.
- **Statuses.** The status of the leading order in the boxed law and of the slow-limit flow is overstated (D3ep-2, D3ep-3). The box-free well form is correctly "not quantified". (ii) is fine.
- **Numerics.**
  - A, B, C and D: sound.
  - E, 'nonlinear only' and '+ barrier': quantitatively confirmed. The ODE gives 0.3565 and 0.3715 at 10⁵ and 3·10⁵ against measured 0.3528 and 0.3692, the expected finite-γ slowdown, and the target is reached.
  - E, 'correlated E' and '+ class levels': not tests of (i) and (iii) as written (D3ep-4).

## Status of the round-31 issues
| Issue | Status |
|---|---|
| D3ad-1 | Resolved. The identity carries the box and cohesion terms; the well form is up to exit terms; Part D has separate columns. |
| D3ad-2 | Resolved: "at leading order", and the O(τ²) term depends on all of Σ (statement and row). The amplitudes-only row was added to Part D. |
| D3ad-3 | Resolved as to the box and to a Σ whose centre and width block is at the dither's level. W♯_τ is now 0.084 from W♯ (separable). But the comparison with predictions is missing (D3ep-4). |
| D3ad-5 | Resolved: (ii) qualified, "no block coupling". |
| D3ad-6 | Resolved: the noise increments Δα − Δtb, and the O(Δt) bias stated. |
| D3ad-7 | Resolved: the residual and its z, the chain-exact bias, martingale errors, the boundary term, and the shared-seed caveat. |
| D3ad-8 | Resolved: the numbers and the docstring now agree with (i). |
| D3ad-9 | Resolved: [S68] and the narrowed "New". |
| D3ad-13 | Resolved: the noise_L, level, per-class, pointwise-identity and generator tests exist and pass. I did not rerun mutants. |

## Consistency with the corrections after rounds 142 and 143
- **Round 143 (1), the leak window.** Not applicable; D‴.40 has no leak.
- **Round 143 (2), the chart of ℓ's divergence.** Not applied to D‴.40. F_τ's divergence and "the weighted log-determinant" are the α-chart's; a noise white in a gives Σ_kτ_k log W_k (D3ep-1).
- **Round 143 (3), non-compactness as the obstruction to Katzenberger.** D‴.40's Katzenberger and slow-limit statements need "on compact parts of M" (D3ep-3).
- **Round 143 (5), D‴.39 at fixed masses.** Consistent: (i) says "at fixed masses … none exists under the cohesion leak alone". The well form needs times short against the exit; the margin at τ = 10⁻⁴ is only about ×7 the slow time (D3ep-3).
- **Round 143 (6), the floored correction does not confine.** (ii)'s "log-barrier" should read "the masses form (floored)". Harmless there, since D = 0 and the noise drift is zero.
- **Round 143 (7), the Péclet number.** The deterministic reading of Part E at γ = 10⁻² fails for amplitude noise (spread 0.03–0.04, Pe about 1) (D3ep-4).
- **Round 142, D‴.41's statuses.** The Euler-chain identity in D‴.40 is exact and needs no fixed-step Katzenberger analogue. The full-loop Euler runs are read through the continuous limit under that assumed analogue, which is harmless at λ_maxΔt = 0.05 but should be said. "Independence over the nodes = per-sample noise" does not touch D‴.40. The roadmap's open item (noise reaching the amplitudes through the measurement) is consistent with D‴.41.

## Overall assessment
The mathematics of D‴.40 is right where it claims exactness: the gauge identity, the Stein/generator identity with the box and cohesion terms, D = 0, and the chain-exact Euler bias. The md reproduces byte for byte, and every decimal is traceable. The weaknesses are scope ("general"), status inflation for the boxed-law leading order and the slow-limit flow, two dropped hypotheses in (iii), and a Part E that illustrates rather than tests.

**Could it support a section of Paper 14?** Yes, after the fixes below, as "white noise of any covariance in the unconstrained chart". The section could claim:
- **(a) Proved:** the exact identity E[u_k] = −τ_kE[g′_k] + E[g_k∂(J_coh + V)] in the boxed law, with the box term ≥ 0.
- **(b) Proved algebra:** the leading-order fitness −∂_{β_k}F_τ, which depends only on each amplitude's own level, in the α-chart.
- **(c) Derived, like der:slow:** the slow-limit flow on compacts, for times short against exit, at large Péclet.
- **(d) Proved:** no amplitude noise means a frozen gauge. With D3ep-7's argument this holds for any noise.
- **(e) Leading order:** the class-level correction, with κ_F > 1 targeting W♯_τ when all τ_k > 0 and W♯_τ ∈ P_B.
- **Remarks:** colour (the filtered level, the quadratic-variation reading fails) and the a-chart (Σ_kτ_k log W_k).

## Issues
| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ep-1 | major | "General noise" means white noise with constant covariance in ϑ. Colour: on the instance the levels fall to 0.47, 0.37, 0.28 of white at τ_c = 100, correlations enter at 8–28%, and the nonlinear Monte Carlo halves or cuts by more than four the fitness (z ≈ 3–5 against white, same seed). The quadratic-variation reading of τ_k fails (α is C¹). Stratonovich multiplicative noise adds an O(τ) term depending on all noise coefficients. In the a-chart ū_k = −τ_k/β_k (round 143 (2)). | Title, l. 4953 statement and row 152: "white noise of constant covariance in the unconstrained chart". Add a remark with the coloured formula −g′_k[(I+τ_cGH)⁻¹Σ]_{αα}/(2G_α), the Stratonovich term −g_kc_α/G_α, and the a-chart result. (iii): the quadratic-variation estimate for white noise only. |
| D3ep-2 | major | "The leading order proved" in the boxed law: the law is non-Gibbs for general Σ, and the plateau weight is non-perturbative at the instance's τ (D‴.39(iii)(b)). | Status (i): "the leading-order fitness at the fit proved; as the small-τ asymptotics of the boxed average, formal (Gaussian/quasi-potential expansion, plateau weight e^{−V/τ} not controlled)". |
| D3ep-3 | major | The slow-limit flow "whatever Σ" rests on D‴.37(ii)'s sketched averaging. Missing: compactness (round 143 (3)), the exit-time window (about ×7 at τ = 10⁻⁴), Péclet (round 143 (7)), and unique ergodicity for degenerate Σ. | (i), the proof and Status: "derived (as D‴.37(ii)/der:slow) on compact parts of M, for times short against the exit, at large Péclet; for degenerate Σ, given unique ergodicity (hypoellipticity)". |
| D3ep-4 | major | Part E is unconfronted. 'Correlated E': the ODE puts W_1 at the wall by 2.2·10⁴; measured 0.343 ± 0.007 at 2·10⁴ against about 0.316 from finite-γ Katzenberger; W_2 falls where the prediction at 1/3 rises; Pe about 1. κ_F = 2: 0.022 from W♯_τ, of which 0.010 is the finite-γ zero (0.4508, 0.2821, 0.2671). | `exp_gauge_noise.py` E: tabulate the slow-limit ODE and finite-γ Katzenberger predictions next to each variant, plus Katzenberger's zero for the class-level run. Add 36+ seeds over 2·10⁴ for the early phase (my rerun was stopped after 1 seed). Numerics (E): replace "move toward" with the comparison. |
| D3ep-5 | minor/major | (iii) needs all τ_k > 0: with some τ_k = 0, κ_F > 1 drives that class to its saturation. With the box it needs W♯_τ ∈ P_B, which fails for levels (2, 1, 0.2) and (1, 1, 0.1). The helper `w_sharp_weighted` crashes at (2, 1, 0.2). | (iii) and the proof: add both hypotheses (cf. D‴.39(iii)(d)), or use the levels max(τ_k, τ_B). Fix the brentq bracket in `w_sharp_weighted`. |
| D3ep-6 | minor | g_k∂_{α_k}V ≥ 0 pointwise is not stated. | (i): add the sign. The box term only opposes the noise's drive into the saturation. |
| D3ep-7 | minor | (ii) is stronger than stated: u_k = dψ_k/ds + g_k∂(J_coh + V) pathwise without amplitude noise, so freezing holds for any noise process. | (ii) and the proof: add the pathwise telescoping. |
| D3ep-8 | minor | Gaussian tails of the boxed non-Gibbs law are asserted; "½E[g‴_k]C"; "log-barrier" for the floored form; "weighted log-determinant". | The proof (Lyapunov function exp(ε\|ϑ\|²)); g‴_k at the fit; (ii) "masses form (floored)"; rows 152 and 9549: "Σ_kτ_k log(W_ka′_k)". Label the τ̂ column's pairs. |
| D3ep-9 | minor | Precedents are missing: Fudenberg–Harris (own-noise Itô correction in the stochastic replicator), Imhof, Hofbauer–Imhof, and SSP (the noise structure fixes the regulariser ∇²L:C and the time scale; Assumption 3.5 is localised to compacts). | [S68], Status "Classical". |

## Sources
- **Repository** at 1510c79: PARADIGM.md ll. 152, 4495–4935, 4953–5014, 5273–5285, 9182, 9549; REVIEW_LOG.md l. 1334; review/REVIEW_round31_D3ad.md; results/chaos/gauge_noise.md and its raw file; experiments/exp_gauge_noise.py; src/core/{gauge_noise,mass_transfer,gauge_drift}.py; tests/test_gauge_noise.py; manuscript/paper14/{README.md, paper14.tex}.
- **Shalova, Schlichting, Peletier**, "Singular-limit analysis of gradient descent with noise injection", arXiv:2404.12293 [paper]. **Li, Wang, Arora**, arXiv:2110.06914 (abstract) [paper].
- **[memory]** Katzenberger, Ann. Probab. 19 (1991) 1587; Funaki, PTRF 102 (1995); Fudenberg–Harris, JET 57 (1992); Imhof, Ann. Appl. Probab. 15 (2005); Hofbauer–Imhof, Ann. Appl. Probab. 19 (2009); Rugh, PRL 78 (1997); Jepps–Ayton–Evans, PRE 62 (2000).
- **My scripts and outputs** are in …/scratchpad/ref147/work. The ones behind the numbers above are ode_E.py, katz_corrE.py, katz_finite.py, zero_corrE.py, colour_lr.py, mc_colour.py with mc_tc100.log and mc_white.log, wsharp_PB2.py, transient.py, peekE.py and peekE2.py. The early-phase rerun (e_early.py, e_early.log) was stopped after one seed.

## Lesson
An identity that is exact for one noise model can still be presented as "general" while leaning on that model three times over: white noise, additive noise, and a fixed chart. Each relaxation changes the leading order, here by factors of 2–10 on the instance for an ordinary OU colour. Name the noise class and the chart in the title. Write the identity for the next class out (coloured, Stratonovich, another chart) to see what survives; here the freezing of (ii) survives and the levels do not. And put every simulation next to the prediction it is supposed to test, including the finite-γ zero, before calling it confirmation.
