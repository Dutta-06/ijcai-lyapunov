# Referee report D3ak (round 38): Proposition D‴.49 (the rich transit in the closed loop), objects at b59721f

## What I checked and how
- **Setup.**
  - I worked only in a scratch copy, `SCRATCH/repo`, at HEAD 2773f33. The four D49 objects there are identical to b59721f.
  - The working repository's HEAD moved to 0369345 during the review. That commit only adds D‴.51 functions to `rich_transit.py`; `OrbitBumps`, `moving_rhs`, the experiment and the test are unchanged, so every finding below still applies at HEAD.
  - Nothing in the working repository was touched, and I ran only read-only git commands there.
  - My scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review38/work` and the logs in `…/review38/logs`.
- **Tests.**
  - `tests/test_rich_transit.py`: 11/11 pass (12.8 s).
  - The full suite reached 23% before the 15-minute cap; everything up to that point passed.
- **Reproduction.**
  - The repository's `averaged` reproduces Part A (ε = 10⁻⁵) and Part C (ε = 10⁻⁶) exactly.
  - `report()` run on the committed raw npz regenerates the committed md with 0 differing lines.
  - The real loop takes 4155 s in the repository's code, so I re-ran it with my own C integrator of the same ODE (DOPRI5, `rloop.c`/`rl.py`/`cex.py`), which keeps Q exactly.
  - Against the repository's `real_loop` (T = 1.2·10⁴, ε = 10⁻²) the full state agrees to 7·10⁻⁸ in x, 3.5·10⁻⁷ in c and 1.6·10⁻⁷ in w.
  - It reproduces every real-loop row of the md:
    - B at 10⁻⁵: entries 5: 132.85, 17: 157.67, 10: 254.13. The md's 133.2, 158, 254.3 are read on an output grid with Δτ = 0.41. Coefficients [0.9, −0.5993, 0.7499], inactive displacement 0.0346, max|e| 0.0894 then 3.8·10⁻⁸.
    - B at 10⁻⁴: five kernels, inactive displacement 0.119, final error 1.66·10⁻⁵.
    - C: 10: 76.28, 17: 92.19, 16: 95.08, coefficients [−0.5998, 0.6992, 0.0501].
- **Numbers.** I checked every number in numerics_D49.md against the md.
- **Independent numerics.**
  - The **exact averaged loop**: the drive taken along the frozen-parameter periodic orbit x_m + e_p, where e_p = L[r(x_m + e_p)] is found as an FFT fixed point (`exact_avg.py`).
  - First-entry drives under three definitions of the drive.
  - The sliding inverse-scale-space path with `OrbitBumps`.
  - The Jacobian of the averaged active block.
  - The real loop's motion within one period, and a trace of how formed kernels are expelled.
  - Real-loop sweeps: k_e ∈ {2, 3, 4, 5, 7, 10, 20, 40}, ω = 1, amplitude a = 0.8, the regressor at x_m (also with f̂ at x_m), measurement noise, a quasi-periodic reference, and a bump just outside the range.
  - Averaged sweeps: k_e, ω ∈ {0.05, 1}, a ∈ {0.6, 0.8}, the half-range reference, bump positions near the range, close bumps, and a briefly visited bump.
  - Part D at ε = 10⁻⁸ and 10⁻¹⁰, and out to τ = 10⁵.
- **Mutation testing.** 30 single-edit mutants (`mutate.py`). Each was run against the test file and against the experiment's outputs (`mut_effect.py`).
- **Literature.** Web search only. Details marked *(memory)* are not confirmed.

## Verdicts

**(i) The averaged loop.**
- **`OrbitBumps`.** It is the correct periodic steady state of the *linear* equation ė = −k_e e + r(x_m(t)).
  - With numpy's convention, E_k = R_k/(iω_k + k_e) and ω_k = 2πk/T. The grid excludes the endpoint, the mean gets the gain 1/k_e, and the Nyquist bin's error is negligible. The test pins the sign of iω.
  - Neither the filter's phase nor its dynamics matter here. Flipping the sign of iω changes no reported number, because the first-order term vanishes on a scalar orbit. Replacing the filter by the static gain 1/k_e moves the entries by about 1% (141.0, 155.5, 251.0).
- **The 1/k_e factor.** It is right asymptotically: k_e·max|g| = 0.0703 → 0.0709 for k_e = 10 → 100.
  - The entry times are 13.3, 12.2 and 15.5 times D‴.48's, not "13–16". The factor is k_e times the ratio of the arcsine to the uniform measure (0.0712 vs 0.0933), plus a 1–2% filter term.
  - The sliding path with `OrbitBumps` reproduces Part A (142.34, 156.54, 254.40), so "D‴.48's transit with this form" holds for the draft's averaged loop.
- **The truncation fails at the experiment's parameters**, where |f|/(k_eσ) = 0.75.
  - First entry: the draft's form gives 142.3, the frozen orbit 132.8, the real loop 132.85 (a 7% error). At k_e = 2 they are 33.9, 9.82 and 9.85 (a factor 3.4).
  - The error halves with each doubling of k_e (3.6% at 20, 1.8% at 40). The near-tie margins in the draft's form are only 6–7.5% (kernel 5 vs 4: 142.3 vs 151.0; 16 vs 17: 156.5 vs 168.2).
  - The exact averaged loop takes 13 s and reproduces the real loop at k_e = 10, ε = 10⁻⁵:

    | | entries (kernel 5, 17, 10) | coefficients | inactive displacement |
    |---|---|---|---|
    | exact averaged | 133.0, 157.75, 254.0 | [0.9, −0.5992, 0.7499] | 0.0334 |
    | real loop | 132.85, 157.67, 254.13 | [0.9, −0.5993, 0.7499] | 0.0346 |

  - It matches equally at k_e = 5 (64.5/80.25/128.75 vs 64.47/80.02/128.78), k_e = 20 (five kernels), k_e = 40 (554.25, 628.75 vs 554.24, 628.66), ω = 1 (four kernels) and a = 0.8.
  - The draft's form gets the kernel identity wrong at k_e = 2, 5, 10, ω = 0.05, ω = 1, a = 0.6, a = 0.8 and on the half-range reference.
  - It gets the count wrong at k_e = 2 (5 vs 3), 20 (3 vs 5), ω = 1 (5 vs 4) and a = 0.8 (5 vs 3). It is right only at k_e = 40.
  - So the status line "kernels' identity not predicted by the averaged loop" describes the truncation, not averaging.
- **Proof sketch.**
  - The periodic steady state solves an implicit, nonlinear equation. It needs k_e > Lip(f − f̂) ≈ 4.5, which is not stated.
  - The multiplicative error (1 + O(|f − f̂|/(k_eσ))) is wrong where the leading drive is small (near ties and tails). The correction is additive, O(|r|²/(k_e²σ)).
- **"Never formed"**: false; see (iii).

**(ii) The real loop.**
- **The explanation of the kernel identity is right but incomplete.**
  - At k_e = 40 the real loop recovers the draft's kernels 5 and 16.
  - With the adaptation regressor at x_m it forms the first bump with kernel 4 (τ = 136.3; ⟨e_pφ_i(x_m)⟩ predicts 136.3), then 16 and 11. Three natural drives pick three different kernel sets.
  - Averaging gives O(rate) closeness, but the entries are near-saddle events whose identity O(rate) perturbations decide. Count and identity therefore need margin conditions, as in D‴.48(iii).
- **The averaging hypothesis is false on the active block.**
  - At the fitted state the averaged Jacobian's centre eigenvalues are −2.70, −2.08 and −0.98 per unit time (mobility Gmw). That is 5–13.5 times 2π/T = 0.2 and only 4–10 times below k_e. The amplitude eigenvalues (−7·10⁻⁴) are slow.
  - Near a fitted kernel the (e, c) pair has roots −5 ± 9–15i, faster than k_e.
  - In the real loop the entered centre swings by 0.060 at τ = 132.86 and by 0.046 at τ = 140 within every period, for at least 100 periods. At k_e = 2 the swing is 0.24 (two kernel widths).
- **Low gain expels formed kernels** (ε = 10⁻⁵).
  - k_e = 2: six kernels. Three sit at the bumps; three were expelled to −1.482, 1.498 and 1.513 with w = −0.508, 0.773, 0.817.
  - The trace shows each was formed first. Kernel 17 sat at 0.633 with w = 0.812 until τ = 63.5, then jumped to 1.343. Kernel 10 left at τ = 65.5. Kernel 16 was expelled at τ = 72.7.
  - k_e = 3: six kernels, three expelled. k_e = 4, 5, 7 and 10: three kernels.
  - The proposer's own scratch log `d49/closed_1e-4.log` (k_e = 2) already shows a kernel at 1.4826 with w = 0.8151; I reproduce it as 1.4907, 0.8151. It is not reported.
  - "A smaller mass forms fewer kernels" fails here: five kernels at 10⁻⁴, six at 10⁻⁵.
- **The count at ε = 10⁻⁵ is not robust.**
  - k_e = 20: five kernels (pairs at 0.07 and 0.61).
  - ω = 1: four kernels.
  - Quasi-periodic x_m = 0.6 sin ωt + 0.4 sin √2ωt: four kernels.
  - Briefly visited bump (both averaged forms): four kernels.
- **"The others keep their starting positions."** Kernel 6 moves 0.119 (about one width) at 10⁻⁴, and 0.15–0.21 at k_e = 3–4.
- **Tracking error.**
  - It stays at the unadapted 0.089 for about 9400 periods, then at 0.074–0.059 for about 8600 more.
  - It then plateaus at 3.8·10⁻⁸. At 10⁻⁴ it is still 4.5·10⁻⁶ at τ = 1000, because the pair at 0.61 is slowly merging (0.6094/0.6256).
  - Part C's real loop ends at 9.4·10⁻⁵, which is not reported.
- **Measurement noise** (not discussed; the free law has no leak). ε = 10⁻⁴, noise in x at 10–30 rad/s, τ ≤ 600:

  | noise rms | outcome | final tracking error |
  |---|---|---|
  | 0.005 | five kernels with biased centres (0.014, 0.082, 0.755) | 0.013–0.021 (clean: ≤ 1.7·10⁻⁵) |
  | 0.01 | 14 kernels, 12 expelled to 1.392 with w = 0.55–0.87 each | ≥ 0.093 (unadapted: 0.089) |
  | 0.02 (two seeds) | 16 kernels, all expelled to about ±1.3 | 0.10–0.12 |
  | 0.05 | w up to 18.6 and −5.4 | 0.205 (identical entries at rtol 10⁻⁹) |
  | 0.1 | w = 25.6 and −10.4 | 0.43 |

  An averaged estimate adds s²|H|²⟨φ′⟩ to every kernel's dual coordinate; the edge kernels have ⟨φ′⟩ = ∓0.69.

**(iii) The price of parsimony.**
- **"Never formed" is false near the range.** These runs use the draft's own `averaged` with x_m = 0.45 + 0.5 sin ωt at ε = 10⁻⁶.
  - A bump at −0.16 (0.9σ outside the range) is formed exactly: kernel 9 enters at τ = 688.5, and the fit is centres [−0.16, 0.07, 0.61], coefficients [0.9, −0.6, 0.75], rms 2.2·10⁻⁷.
  - At 1.6σ, a kernel sits at −0.2345 with w = 0.844 by τ = 5000.
  - At 2.9σ, nothing has formed by τ = 5000.
  - For Part C's own bump at −0.52, the sliding path keeps entering kernels (τ = 1.7·10⁶, 6.2·10⁶). That late path is not relaxed, so only its early entries are reliable.
  - Meanwhile the *visited* neighbour is mis-formed: 0.158 with −0.384 at τ = 600 (rms 0.10). In the real loop (ε = 10⁻⁴) it sits at 0.084 with −0.49 and the tracking error is stuck at 0.020 from τ = 900 to 2000.
- **Transit time ∝ k_e only asymptotically.**
  - τ₁/k_e in the real loop is 4.9, 11.9, 12.7, 12.9, 13.1, 13.3, 13.6 and 13.9 for k_e = 2, 3, 4, 5, 7, 10, 20 and 40; the limit is about 14.0.
  - "Exactly ∝ ε^{−1/2}" restates t = Λ(ε)τ with τ read on a 0.5 grid. Between 10⁻⁴ and 10⁻⁵ the ratio is 2.94, because the last entries are different events.
- **Part D compares at one ε and one horizon.**
  - Wasserstein–Fisher–Rao forms three kernels at 10⁻⁸ and 10⁻¹⁰, with the same entrants as the paradigm's metric.
  - Fixed centres are still mid-transit at τ = 600. Kernels keep entering until τ = 8.4·10⁴, and the rms falls from 0.0127 to 0.0016 at τ = 10⁵. The grid's least-squares fit is 8.2·10⁻⁵.
- **Close bumps.** At 0.83σ apart a single merged kernel (w = 1.151, rms 0.0185) persists from τ = 130 to 835.5, then the bumps are resolved.

**Novelty.**
- "Which features form is decided by the reference" is partial persistent excitation along the orbit:
  - Wang–Hill deterministic learning (IEEE TNN 2006): neurons along a periodic orbit converge, the others are barely updated *(the latter from memory)*.
  - Kurdila–Narcowich–Ward (SIAM J. Control Optim. 1995) *(authorship from memory)*.
- "One unit at a time along the trajectory, compact" — all with explicit growth rules:
  - Platt's resource-allocating network (1991), MRAN (1997), GAP-RBF and GGAP-RBF (2004–05).
  - Fabri–Kadirkamanathan, dynamic-structure Gaussian networks that grow where the state travels (IEEE TNN 1996).
  - Sundararajan–Saratchandran–Li, fully tuned RBF control with MRAN.
  - Kingravi–Chowdhary–Vela–Johnson, online centre restructuring in model-reference adaptive control (TNNLS 2012).
- "No growth rule":
  - Sparsity-promoting implicit regularisation of mirror-descent adaptation laws, Boffi–Slotine (2021); not cited for this.
  - One-at-a-time activation from small initialisation: Pesme–Flammarion; Kunin et al.'s alternating gradient flows (2025); D‴.48 itself.
- So D49 is D‴.48 combined with leading-order classical averaging. Its genuinely closed-loop phenomena — near ties decided by the orbit, expulsion at low gain, expulsion under noise — are not in the statement.

**Code and tests.**
- The only D49 test is `test_orbit_filter_is_the_periodic_error`, which feeds a zero-mean sinusoid.
- 6 of 30 mutants are killed, all in `OrbitBumps`. The killed sign-of-iω mutant changes no reported number.
- 24 mutants survive. Each changes the experiment's outputs except the 64-sample mutant:
  - `OrbitBumps` dropping the mean (H(0) = 0): reorders the entries to 5: 184, 11: 191, 16: 213.
  - All 5 `averaged` mutants: no filter (τ ÷ 10), a one-sample lag, k_e doubled, 64 samples (no change), τ normalised by k_eΛ.
  - All 8 `real_loop` mutants: regressor at x_m for the amplitudes or the centres, the error's sign, a position law without the mass, 1/σ instead of 1/σ², the position law's sign, f evaluated at x_m, a wrong reference in the reported error.
  - All 5 entry and threshold mutants.
  - All 4 `report` mutants.
- The md's tables A and B do not render as tables: unescaped "|" makes the headers 10 and 12 cells against 8.

**Numerics paragraph.**
- Every transcribed number matches the md.
- Wrong or overstated:
  - "13–16" should be 12.2–15.5.
  - "Exactly ∝ ε^{−1/2}" is tautological.
  - "Never formed" holds only within the horizon.
  - The real-loop entry times are read on a 0.41 grid but quoted to 0.1.
  - Part C's averaged pair at 0.07 versus the real loop's pair at 0.61 is not discussed.

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ak-1 | major | (i)'s averaged loop takes the orbit to be the reference, which is only the leading-order truncation. It misplaces the first entry by 7% (a factor 3.4 at k_e = 2) and gets the identity or count wrong in 9 of 10 settings. The frozen-orbit averaged loop reproduces the real loop. | D49_draft (i), proof, status; `rich_transit`: add the exact frozen-orbit drive (requires k_e > Lip(f − f̂)) and use it in Part A |
| D3ak-2 | major | The averaging hypothesis fails on the active block. Centre rates of 0.98–2.70 are 5–13.5 times 2π/T; the (e, c) roots are −5 ± 9–15i; the centre swings 0.046–0.06 per period (0.24 at k_e = 2). | (ii) and proof: average the inactive kernels only and treat the active block quasi-statically, stating Gmw²⟨(∂_cφ)²⟩/k_e ≪ 2π/T, or label the agreement numerical; report the swing |
| D3ak-3 | major | At k_e = 2–3 formed kernels are expelled to \|x\| ≈ 1.5 with \|w\| = 0.51–0.82, giving six kernels. The proposer's own k_e = 2 scratch run shows this and it is not reported. "Fewer kernels at smaller mass" fails. | (ii), (iii): state the gain range or a condition; Part B: sweep k_e from 2 to 40 |
| D3ak-4 | major | At ε = 10⁻⁵ the count is 5 at k_e = 20 and 4 at ω = 1, with a quasi-periodic reference, and with a brief visit. One kernel per feature depends on where the near ties fall. | (ii): claim one per feature only for ε ≪ μ³, with margins; sweep k_e, ω and a in Part B |
| D3ak-5 | major | "Never formed" is false within about 1–2σ of the range: a bump at −0.16 is formed exactly at τ = 688.5. The visited neighbour is mis-formed meanwhile, and the real loop's error is stuck at 0.020. | (i), (iii), numerics C: define "visit" via the occupation measure convolved with the kernel; give the threshold; add an edge-bump sweep |
| D3ak-6 | major | Noise of 1% of the reference expels 12 kernels and leaves tracking worse than with no adaptation; 5% drives w to 18.6. No robustness statement. | Draft: make "noise-free" a hypothesis and add the counterexample; test with D‴.46's leak |
| D3ak-7 | major | One test covers D49 and 24 of 30 mutants survive, including the dropped mean and every mutant of `averaged`, `real_loop`, the entry detection and `report`. | `test_rich_transit`: filter test with a nonzero mean and several harmonics; sliding path against Part A; a short real-loop run against the exact averaged loop |
| D3ak-8 | minor | Transit time ∝ k_e only for large k_e (τ₁/k_e = 4.9 at k_e = 2); "exactly ∝ ε^{−1/2}" is the definition, and the 10⁻⁴ to 10⁻⁵ ratio is 2.94. | (iii); numerics A |
| D3ak-9 | minor | Part D: Wasserstein–Fisher–Rao forms 3 kernels at 10⁻⁸ and 10⁻¹⁰; fixed centres are still entering kernels at τ = 8.4·10⁴. | Part D: ε sweep and full horizon; (iii) |
| D3ak-10 | minor | The price is also in tracking: the unadapted error persists for about 9400 periods; the end is a plateau (4.5·10⁻⁶ at 10⁻⁴; Part C 9.4·10⁻⁵, unreported). | (ii), (iii) wording; report e_last for Part C and the error history |
| D3ak-11 | minor | "Keep their starting positions" fails: 0.119 at 10⁻⁴ and 0.15–0.21 at k_e = 3–4. | (ii) |
| D3ak-12 | minor | The near ties also flip with the regressor's evaluation point: at x_m the first bump goes to kernel 4. | (ii) wording |
| D3ak-13 | minor | The steady state is implicit and needs k_e > Lip(f − f̂); the multiplicative error bound is wrong; the O(k_e⁻³) symmetry remark concerns a phase that does not matter. | Proof sketch (i) |
| D3ak-14 | minor | Novelty is overstated; Wang–Hill, Boffi–Slotine, Fabri–Kadirkamanathan and the RAN/MRAN/GAP-RBF line are uncited. | Status / "New"; literature notes |
| D3ak-15 | minor | Close bumps at 0.83σ give a merged kernel with w = 1.15 for τ from 130 to 835. | (ii): quantify "well-separated" |
| D3ak-16 | cosmetic | md tables A and B do not render; entry times quoted beyond their grid; "13–16"; Λ vs k_eΛ is ambiguous (the τ-normalisation mutant survives). | `exp_closed_transit.report`; numerics |

## Lesson
The draft was tested at one point of the controller's own parameter space — k_e = 10, ω = 0.2, a = 1, no noise — where the leading-order picture happens to give the right features.
- **Compute the averaged system the theorem actually delivers** — the frozen nonlinear orbit, a 13-second job — not its truncation.
- **Sweep the loop's own knobs** (gain, frequency, amplitude, reference shape). Each one crossed a near tie or an instability.
- **Move the "unvisited" feature toward the range** instead of keeping it four widths away.
- **Add the smallest realistic perturbation** — 1% measurement noise — to a law with no leak.

For a controller proposition, the parameters to vary are the controller's and the plant's, not only ε.
