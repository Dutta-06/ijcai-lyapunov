# Referee D3aw, round 50: Proposition D‴.57 (the transit's dead zone)

I did not modify the repository. At the end I checked that D57.md, exp_dead_zone.py, dead_zone.md, test_dead_zone.py, closed_loop.c/.py
and exp_closed_transit.py are byte-identical to my copy. All my processes are finished, and I killed only my own PIDs (28005, 28006).

SCRATCH = the session's scratchpad `review50/`. Scripts are in SCRATCH/scripts and outputs in SCRATCH/out.

## 1. What I checked and how

- **Regeneration.** I copied the repository to SCRATCH/repo (HEAD 502c522; another agent's uncommitted edit to exp_noisy_certificate's
  report does not touch `targets`/`noiseless`). Running `exp_dead_zone.py --from-raw` reproduces dead_zone.md byte for byte.
- **Theory.**
  - I re-derived the identity, the orthogonality and the cross-term argument.
  - I checked every χ² number with mpmath, both by quadrature and by the χ²_{k+2} identity: P(late) = 0.046085, m_c = 53.82663,
    E|u|²/N = 0.327596 s² (rms 0.57236 s), E err² = 0.425957 s², ratio 2.08100.
- **Part A re-run.** I re-ran all 1200 draws (same seeds) with the entry order and the vectors u, f̂_K, f̂_c recorded
  (SCRATCH/scripts/order_cos.py). It reproduces 386/381/356 and the |cos| medians exactly.
  - I decomposed u on the fit's tangent basis.
  - I linearised the fast flow (`critical_rhs`) at each target's noiseless block (SCRATCH/out/slow.log).
- **Part B.** I took the two worst draws (target 11, s = 0.1, draw 1; target 6, s = 0.1, draw 2) and ran them at rtol 10⁻⁹ and 10⁻¹¹. I
  compared the loop's fits with the limit's fits as vectors, not only as rms errors.
- **Part C.** I re-ran the C loop with a dense record (1201 output times, plus the filtered energy z, the entry times and the
  trajectories) for these cases at x_m:
  - s = 0.02: formation gate; energy gate (seeds 7 and 8); ungated, seed 8; energy gate with thresholds 2.0s/2.5s.
  - s = 0.02, two C variants in copies outside the repo: f̂ evaluated at the true state in the control; adaptation driven by x − x_m
    instead of x̃ − x_m.
  - s = 0.01: ungated to τ = 800.
- **Mutation tests.** I made three mutants in copies (SCRATCH/mut/M1, M3, M8) and ran the whole test file on each.
- **Literature.** I used web search. arxiv/ar5iv and projecteuclid pages would not load, so items marked "memory" come from memory.

## 2. Verdicts

**(i) The formation dead zone is the discrepancy principle — correct in the sliding limit on fixed data, but close to tautological.**
- Freezing the unformed duals is, by construction, the inverse-scale-space path stopped by the discrepancy rule.
- Part A does not simulate this gate. It takes the polished K-th block as its outcome. Only Part B tests it (46, 46, 44 draws), and there
  it holds.
- A formed kernel whose |w| drops below ϑ would be frozen at about ϑ instead of exiting. Only the no-exit hypothesis excludes this.
- The implementations are not gate (b):
  - Part B applies a one-time mask at the closing. This is equivalent on fixed data.
  - The C gate is smooth between θ/2 and θ. With the gate shut, kernels in that band keep a drive s_i ∈ (0, 1).
- "No noise kernel forms at any later time" holds only for the sharp gate on fixed data.

**(ii) The energy dead zone freezes the last block.**
- **Identity.** |u|² = (1+δ)²s²N − |r_K|² is exact Pythagoras plus the fit's first-order condition Jᵀr_K = 0; the error is
  ⟨r_K, curvature⟩ = O(s³). Measured: 0.99946–1.00075 at worst.
- **Law.** The law is correct at leading order for a fixed target and a fixed entry order. The cross-term argument breaks when the entry
  order is random:
  - u lies along the last entrant's coefficient. Which kernel enters last is decided by components of ξ that include P_Tξ.
  - Conditional on a swapped order, ⟨P_Tξ, u⟩ has a large negative mean.
  - On draws whose last entrant differs from the noiseless path's (5/386, 5/381, 29/356): correlation −0.32 ± 0.08, −0.31 ± 0.15,
    −0.41 ± 0.06. Mean err_dz² over its prediction: 1.53 ± 0.27, 1.12 ± 0.33, 1.54 ± 0.11.
  - On the noiseless-order draws: 1.016, 0.978, 0.958 ± 0.028.
  - So the reported 1.006 at s = 0.1 is a mixture. The o(s²) is not uniform in the entry-order margins.
- **Slowest-mode step: not the mechanism.**
  - λ₂/λ₁ ranges from 1.0025 to 2.20 (≤ 1.1 on 11 of 40 targets), so s^{λ₂/λ₁−1} is not small at s = 0.01–0.1.
  - The alignment comes from the other modes' amplitudes being O(feature overlap), not O(1): the block before the entry is already
    relaxed.
  - u is the last entrant's amplitude deficit: ≥ 99.4% of its norm on targets 4 and 10. The deficit is δw ≈ 1.76 s (median; 5–95%
    [0.7, 2.4] s), which is 31% of that amplitude at s = 0.1.
  - It coincides with the slowest mode only because the last entrant has the smallest fitted |w| (true on 40 of 40 targets).
  - At near-degeneracy, "v₁ is the noiseless block's up to O(s)" is false. On target 20 (λ₂/λ₁ = 1.003) the slowest eigenvector swaps
    between draws.
  - The law does not need this step: at leading order u depends on ξ only through the crossing time, i.e. through |R_Kξ|.
- **Numbers.** 0.4260 s², 2.08, 0.57 and P(late) = 0.046 are correct.
- **"About twice".** This is a coincidence of N/K ≈ 20 and δ = 0.1. The ratio is (2 + N(2δ+δ²)/(2K))^{1/2}: 3.0 at N = 200 and 6.1 at
  N = 1000. The excess per sample does not vanish with N.

**(iii) "A closed gate stays closed" — correct for both sharp gates on fixed data.**
- The energy gate leaves the state stationary.
- Under the formation gate the formed block runs a gradient flow, so the energy cannot rise.
- "Theorem D‴.53 on y" assumes (H1)–(H4) for noisy data. (H3) was checked only on noiseless targets.

**(iv) Closed loop — the failure is real, but its description is wrong.**
- The negative result for x_m at s = 0.02 is reproduced and robust to the noise's phases (seed 8: ungated 12 kernels / 11 outside / fit
  0.386; energy 7 / 3 / 0.121).
- The gates chatter (D3aw-4).
- The "12 formed kernels" are an eject-and-replace cycle on the features at −0.52 and 0.61, not noise kernels (D3aw-5).
- Nothing supports the "fits the noise at the edge" explanation: the ejected kernels end at |c| ≈ 1.5, where φ(x_m) ≤ 3·10⁻⁵.
- **Mechanism tests.**
  - With f̂ evaluated at the true state in the control, nothing changes (12 / 11 / 0.387).
  - Driving the adaptation with x − x_m keeps the 0.61 feature (8 kernels / 6 outside / 0.308).
  - This points to the noise in the adaptation's error moving the formed kernels' centres, whose mobility is ∝ w. At s = 0.01 the centre
    of kernel 16 wanders over [0.550, 0.673].
- **Overgeneralised.**
  - At x_m, s = 0.01, the kernels stay in the window to τ = 800 (3 / 0 / 0.0225).
  - The thresholds decide the energy gate's outcome: at 2.0s/2.5s it ends with 1 kernel and none outside.
  - A bound on the centres is never run.
  - The x̃ rows are confounded: without noise the x̃ loop already has 5 kernels (duplicates), and noise enters its regressor.

**Numerics.**
- Every number in the draft matches the md, except "|u|/(s√N) 0.569, 0.566, 0.570". Those are hand-computed medians (not in the md),
  compared against an rms (0.572); the conditional median is 0.574.
- **Part A.** The energy gate is implemented as defined; the formation gate is not simulated.
- **Part B.** Converged:
  - Unchanged at rtol 10⁻¹¹.
  - The loop's fits match the limit's to 1.3·10⁻⁸ rms; the draft's 2.2·10⁻⁹ compares rms errors only.
  - Its noise vectors are the same at the three levels (D3aw-13).
- **Part C.** Robust to the noise's phases, but the fit, coefficient and centre values are single snapshots (D3aw-7).

**Tests.** All 5 pass, and all 5 still pass on each of these mutants:
- M1: the C formation gate replaced by the energy gate (`gi = g`).
- M3: m_c with χ²_{N−2K} in the numerator. The md changes (the "given the stop" ratios become 1.071, 1.026, 1.054).
- M8: the closing time taken on the upper threshold, and the late count dropping k = −1. The md's "later" column becomes 0, 0, 0.

**Novelty.** The mechanism is classical:
- A dead zone freezes the parameters where the error enters it, leaving the residual at its boundary: Egardt 1979; Peterson–Narendra 1982;
  Samson 1983; Ioannou–Sun 1996 §8.5 (memory).
- A continuous flow stopped exactly at the discrepancy level, and its excess over the oracle, is studied in the early-stopping
  literature: Blanchard–Hoffmann–Reiß 2018, who interpolate the estimator so that the residual equals κ (memory); Celisse–Wahl 2021
  (JMLR 22(76)); Blanchard–Mathé 2012 (memory).
- ISS stopped by the discrepancy principle is Burger–Gilboa–Osher–Xu 2006 and the aISS of Burger–Möller–Benning–Osher 2013.
- I know of no source stating "twice by freezing the last block", but "twice" is not a general fact. What is new is narrower: in the
  sliding limit the excess is an O(s), N-independent shrinkage of the last entrant's amplitude.

## 3. Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3aw-1 | major | "u's direction is one per target … falling with s as s^{λ₂/λ₁−1}" is false. At s = 0.01, 0.03, 0.1, 5, 5 and 46 draws (on 2, 2 and 18 targets) have \|cos\| < 0.5; at s = 0.1 the 10% quantile is 0.02. Every one of them has a different last entrant from the reference draw. Draws with the same last entrant have \|cos\| ≥ 0.9996, 0.9954, 0.9388. The falling median comes from order swaps. | D57 numerics A; `a_job`: compare u within one entry order and report quantiles. |
| D3aw-2 | major | The vanishing cross term assumes e is independent of P_Tξ. The entry order (hence e) depends on P_Tξ. On swapped draws err_dz² is 1.53, 1.12, 1.54 × the prediction; the aggregate is a mixture (0.958 on noiseless-order draws at s = 0.1). | (ii) and Status: add the hypothesis "noiseless order kept" (D‴.56's crit_ord < 1) and state the non-uniformity; Part A: use exact_order and report swaps separately. |
| D3aw-3 | major | "Costs about twice the oracle's error" (New, title gloss) holds only for N/K ≈ 20, δ = 0.1. The ratio grows like (Nδ/K)^{1/2}. The energy dead zone is inconsistent in N: it leaves a ≈1.4–1.8 s bias on the last amplitude. The statistical discrepancy principle uses τ² − 1 = O(N^{−1/2}). | Title and New: state the N-dependence and the amplitude bias. |
| D3aw-4 | major | Part C's "gate closed at τ" is the first sample of a 3-τ grid with z < lo². In the dense re-run (x_m, s = 0.02): formation gate: first closes at τ = 149.5, closed 11% and fully open 72% of output times; 11 of 12 kernels formed after its first closing (9 after the reported τ = 177); energy gate: shuts for good only after τ ≈ 400. "It holds only the unformed particles" is false. | (iv), numerics C; `c_job`: store entries and z; report the open fraction and the last reopening. |
| D3aw-5 | major | The mechanism is misdescribed. Each kernel formed on the features at −0.52 and 0.61 is ejected 20–90 τ later (4: 136.5 → 164.5; 5: 169 → 197, to +1.57; 3, 6, 2, 7, 16, 15, 17, 18, 14 likewise). Ejected kernels end where φ(x_m) ≤ 3·10⁻⁵ and fit nothing, so "fits the noise at the edge" is unsupported. The energy gate's "kept" features are replacements (kernels 3 and 15); kernel 11 is frozen at w = −0.253 of −0.6 and kernel 17 at w = 0.010. The same holds at seed 8. | (iv): describe the eject-and-replace cycle with timings. Test the mechanism (my variants implicate the noise in the adaptation's error, not f̂(x̃)). |
| D3aw-6 | major | The negative result is overgeneralised: at x_m, s = 0.01 all kernels stay in the window (to τ = 800); gate thresholds 2.0s/2.5s keep the single formed kernel in the window; "only … a bound on the centres" was never run; "no gate on the drive's energy prevents it" rests on one parameter set (Tf = 5, shorter than the period 31.4). | Title, (iv), New: restrict to the tested level, regressor and gate settings, or run the box and other gate settings. |
| D3aw-7 | minor | Part C's fit, coefficients and centres are single snapshots of a jittering state. At x_m, s = 0.01 (ungated), the fit rms over τ ∈ [300, 600] has median 0.0257, 5–95% [0.015, 0.070]; the reported 0.058 is a high sample. The energy gate's frozen 0.062 is about 2.4× the ungated median, not "0.062 vs 0.058". | Report time-medians and quantiles over the last window. |
| D3aw-8 | minor | The proof's slowest-mode step: "a_i = O(1)" is wrong (O(overlap)); s^{λ₂/λ₁−1} is vacuous (λ₂/λ₁ down to 1.0025); "v₁ is the noiseless block's up to O(s)" fails at near-degeneracy (target 20). | Replace with "u depends on ξ through the crossing time"; name the last entrant's amplitude. |
| D3aw-9 | minor | The gates as implemented are not (b): the C gate is smooth and partly drives kernels in (θ/2, θ); Part B uses a one-time mask; a formed kernel falling below ϑ is frozen, not exited; Part A does not simulate the formation gate. | (b), numerics: state the differences or align the code. |
| D3aw-10 | minor | Tests: mutants M1 (C formation gate made the energy gate), M3 (m_c), M8 (closing time, late count) all survive. | test_dead_zone.py: a closed gate with some kernels formed; a report test on a synthetic raw; a `c_job` closing test. |
| D3aw-11 | minor | Missing literature: the dead zone's freezing (Egardt, Peterson–Narendra, Samson, Ioannou–Sun §8.5); continuous discrepancy stopping and its excess (Blanchard–Hoffmann–Reiß 2018, Celisse–Wahl 2021, Blanchard–Mathé 2012); ISS with the discrepancy principle (Burger–Gilboa–Osher–Xu 2006; Burger–Möller–Benning–Osher 2013). | Status. |
| D3aw-12 | minor | Part A's exactness checks the entry set, not the order. `a_draw` indexes `hist[1:nK+1]`, which assumes one entry per iteration. The \|cos\| median includes each reference draw (cos = 1). | Use `NC.classify` and `fit_after`. |
| D3aw-13 | minor | Part B seeds the noise with 5800 + 97q + k, independent of s. The levels reuse the same noise vectors (the "2 late per level" are the same two draws). | Seed by level, or state it. |
| D3aw-14 | minor | The x̃-regressor rows are confounded: the noise-free loop already has 5 kernels, and the noise also enters the regressor. | (iv): say so. |
| D3aw-15 | cosmetic | "\|u\|/(s√N) 0.569, 0.566, 0.570 against 0.572" are hand-computed medians compared with an rms (conditional median 0.574, mean 0.553). The md's own 0.555 is the unconditional rms. | Have the report print the like-for-like statistic. |
| D3aw-16 | cosmetic | "Identity to 10⁻⁴" is true only for the 5–95% range; the extremes reach 7.5·10⁻⁴. "Fit to 2.2·10⁻⁹" compares rms errors; the fits themselves differ by 1.3·10⁻⁸. | Numerics A and B. |
| D3aw-17 | cosmetic | Wording and minor technical slips: "(1+δ)²s²N/N" in the proof of (i); "minimiser over a larger set": monotone decrease follows from the flow, not from minimality; the noise floor ignores the f̂′(x)n term (up to about 20% at 10 rad/s). | Proof; md. |

## 4. Lesson

The fixed-data algebra is sound; both the statistics and the narrative fail for lack of the right per-draw and per-time breakdowns.
Medians hid a bimodal direction set by the entry order. That same order also breaks the independence argument behind the cross term,
which the aggregate mean happened to mask. "First closing" hid a gate that is open most of the time. Single τ = 600 snapshots hid kernels
that jitter by half a kernel width and are ejected and replaced in a cycle. Before writing a mechanism into a proposition, condition on
the discrete events (the order, gate state, entries and ejections) and plot each quantity's time-course. Before writing a ratio into the
"New" paragraph, write down how it scales with N, K and δ.

Sources:
- [Celisse & Wahl, JMLR 22(76), 2021](https://jmlr.org/papers/v22/20-358.html)
- [Blanchard, Hoffmann & Reiß, Optimal adaptation for early stopping (arXiv 1606.07702)](https://arxiv.org/pdf/1606.07702)
- [Blanchard, Hoffmann & Reiß, EJS 2018](https://projecteuclid.org/euclid.ejs/1538121641)
- [Burger, Möller, Benning & Osher, Math. Comp. 82 (2013)](https://www.ams.org/mcom/2013-82-281/S0025-5718-2012-02599-3/)
- [Burger, Gilboa, Osher & Xu 2006](https://www.projecteuclid.org/euclid.cms/1145905942)
- The other references are from memory: arxiv and projecteuclid pages could not be fetched from this session.
