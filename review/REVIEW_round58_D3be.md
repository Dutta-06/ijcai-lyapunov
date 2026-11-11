# Referee report D3be, round 58: Proposition D‴.54 as revised after round 53

## Summary

**What checks out**
- **Reproducibility.**
  - `--from-raw` (TN_OUT pointed at scratch) reproduces `results/chaos/transit_nd.md` byte for byte (md5 fdf1e396…).
  - 33 of 33 tests pass.
  - Four reruns match the raw bitwise: R at x = 10, R's mirror at x = 7, the R-scan cell (4, 10⁻¹¹), and VF at 6σ.
- **The third-moment law of (v) is right.**
  - M₃, the projection (∂³_eφ)_⊥ = σ⁻³(u³ − (3/2)u)φ, the order counting behind the envelope step, and the sign all hold.
  - P, B, P₃ and M₂* agree with quadrature to 8 digits.
  - `e3_law` is the formula, and the "(√πσ)^{n−1}" factor is right.
  - The two sides' secant ratios behave as the law plus an even next order predicts: the unrounded mean at 7σ is 1.00242.
- **The reduction on the rests holds numerically.**
  - R's held runner-up sits at the static rest for its current coefficient to 0.1–0.4%, at all 25 marks I compared, on both sides.
  - dw/ds = −Λ₀w′(α)M(α)∂_wE* reproduces R's far-side growth times to 4.5% and 13%.
- **The far tail's rule s_θ < s₂ holds at its boundary.** Four new loop runs at s₂/s_θ = 0.80–1.24 give the counts the rule predicts.
- **Held coefficients reach the fold before leaving.** At 0.9–1σ they are 0.97–0.98 of FC's fold (dense output).

**Two statements of (v) are wrong, and both feed (iv)'s count**
- **The far side ends in a swap (major).**
  - The third-moment drift never vanishes for 0 < w < A. At w = w_e, E*₃ = 0 but ∂_wE*₃ = −1.00·10⁻⁷ (6σ).
  - So the far runner-up outgrows the entrant. The entrant, now the lighter kernel and displaced toward the unfitted feature, then decays and races.
  - In R's own mirrored problem continued past its horizon (x = 7), the coefficients cross at s = 8.67·10⁸. The entrant leaves at s = 1.40·10⁹ and fits feature 2.
  - The loop on the mirrored target at x = 7, ε = 10⁻¹⁷ ends [1, 1, 1], with the entrant captured by feature 2. (iv) and (v) predict [2, 1, 1].
- **For formed runner-ups, the side and s_d come from a rounding floor (major).**
  - The collapse keeps the layer's direction down to about 10⁻¹³σ. Then the pair pins at a floor of about 3·10⁻¹⁶σ/w, and the floor's direction is −y on C2's target and on its mirror alike.
  - At rtol 10⁻¹³ on the mirror, c_r − c_e flips at s ≈ 450 from the layer's (far) side to toward.
  - So "at 10⁻¹³ every run takes the side the layer keeps" is a coincidence on C2's target. 26 of the 74 rule comparisons involve floor-decided sides.
  - At x = 5.48 and 6.04, the floor-to-σ climb is 94% and 97% of s_d. When the floor is reached, the relative mode is still contracting at 0.135 per unit s.
  - Hence the formed part of s_d(x), its minimum at 5.48 and the capture set's lower edge are not established.

**Minor issues**
- The proof's handling of the linear pull.
- "Collapsed onto the entrant", said of far tails.
- R's s_d is quantised to its 3.8% output grid.
- Hand-computed or misquoted numbers.
- A missing literature neighbour.
- Mutation testing: 16 of 29 mutants survive the tests, 6 of them changing the md.
- Dead outputs.

## What I checked and how

- **Setup.**
  - Reviewed edb5dd3 via `git archive` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref58/D3be_work/repo`. Scripts are in `…/D3be_work/scripts`, outputs in `…/D3be_work/out`, and the mutant copy in `…/D3be_work/mut`.
  - I wrote nothing to the repository and killed no process.
  - HEAD has since moved to cfdf29c (D‴.52/58 sessions). All D54 files equal edb5dd3's except `src/core/rich_transit.py`, whose change only adds defaulted arguments that D54's code does not use.
  - Runs used the venv python with OMP_NUM_THREADS=1, nice -n 10 and PYTHONDONTWRITEBYTECODE=1, with at most two heavy processes of mine at once. The load average was 13–14 on 4 cores.
- **Regeneration and tests.** `cmp` identical; `pytest -p no:cacheprovider` gives 33 passed in 48.9 s.
- **Reruns against the raw** (`rerun_cells.py`):
  - `race_job(('near tie', 10))`;
  - `race_job(('near tie mirrored', 7))`;
  - `race_job(('near tie', 4, 1e-11, 3e6))`;
  - `v_job((1, 6, 0.7, −1))`.

  Every field is identical.
- **(v), static.**
  - Re-derived the moments, projections, P, B, P₃ and M₂*, including the bookkeeping at O(|δ|⁴) and O(|δ|⁵) with the first moment free.
  - Quadrature at 5, 6, 7σ (`e3check.py`).
  - The derivative of `e3_law` across w_e.
  - The ratios from the raw (`e3raw.py`).
  - V's far branch continued from w = 0.001 to 0.78 (`far_up.py`).
- **(v), dynamics.**
  - R's problem continued to s = 10¹¹ (`far_long.py`).
  - The loop at ε = 10⁻¹⁷ on the mirror (`loop_swap.py`).
  - Marks against static rests, and the reduced rate (`slave_check.py`, `rate_check.py`).
  - Dense output at the end of the hold (`hold_end.py`).
  - c4_job at the far rule's boundary (`c4m_boundary.py`).
  - R at x = 30.
- **The side.**
  - c_r − c_e traced in x and y through the collapse, at rtol 10⁻¹³ and 10⁻¹¹, on three targets: C2's target, its x-mirror (the repository's) and a y-reflection (`race_trace.py`, `race_xy.py`, `race_flip.py`).
  - The full variational equation along the computed trajectory, the Jacobian's eigenvalues and λ_u per gap (`var_full*.py`, `lam_u.py`).
- **Mutants.** 29 single replacements in `exp_transit_nd.py`, each run against `pytest -x` and against `--from-raw` with the md diffed (`mutate.py`).
- **Numerics paragraph.** Every (v) number checked against the md and raw; spot checks of A, T, B, B4, C, D, D2, H and W (all match).
- **Literature.** Web searches. arxiv.org was not reachable, so I used jmlr.org, proceedings, ar5iv and epubs snippets.

## Verdicts per part

**Parts (i)–(iii)** were not changed in round 53, and their quoted numbers match the md. T3's new frozen-tail column is printed, but the numerics overstate it (D3be-7).

**Part (iv)**
- The count and the codimension law stand.
- The rewritten exception sentence is grammatical but has two problems:
  - it misses the far side's swap, so its "+1" holds only in a window and "the reach of all three grows as ε → 0" is false (D3be-1);
  - its "−1" rests on sides the numerics do not resolve (D3be-2).

**Part (v)**

| claim | verdict |
|---|---|
| second-moment law; E* − E₀ = −P²/(2B), the same on both sides and independent of w | holds |
| M₃ = w(1 − w²/w_e²)δ³, P₃, envelope step, sign, ℝⁿ factor | holds (quadrature ratio 1.00000000 at 5, 6, 7σ) |
| two-side mean "1.0025", even next order cancels | the explanation is right (even terms enter the two ratios with opposite signs); the number is hand-computed (1.00242) |
| toward: E* increases, w decays to the fold | holds (leaves at 0.97–0.98 of FC's fold) |
| far: E* decreases, no fold | holds, up to w = 0.78 |
| far: "drift vanishing only at w = w_e", "a formed runner-up stays bound", the far +1 | wrong (D3be-1) |
| dw/ds = −(\|𝒲′\|/ρ)∂_wE* (slaving) | holds numerically; the proof omits the time-scale separation, which fails at the fold (D3be-4) |
| δ = 0 "a saddle … left along the direction it has" | a degenerate saddle (−kP along e, flat to second order across e), left along ±e with the sign of the e-component; wording only |
| side rule (c_r⁰ − c_e⁰)·e > 0 | established for tails only (x ≥ 7, closest approach ≥ 1.8·10⁻⁸σ: C2's target toward, mirror far at x = 7, 10, 17.3); untested for formed runner-ups (D3be-2) |
| race: captured iff s_d < s₂, s_d ε-free | holds for tails (48 of the 74 comparisons); the formed part is floor-limited (D3be-2) |
| capture set an interval of gaps widening on both sides | R's 23 gaps: s_d falls on 2.5–5.48 and is non-decreasing on 5.48–24.5 (3.26·10⁶ at x = 30); the lower branch is floor-limited |
| far: formed iff s_θ < s₂ | holds at the boundary (D3be-5), and only while s₂ < s_{d,e} |

**PARADIGM.md D‴.51(ii)(c)** (l. 6224; table l. 160) repeats the far-side error: "held … and its coefficient grows" (D3be-1).

**D3az issues**

| D3az | status | note |
|---|---|---|
| -1 (major, far side) | partly fixed | Both sides, the third-moment law, VF, FC, R's mirror and C4M added; all right within R's 3·10⁸ horizon. The end state is wrong (D3be-1). |
| -2 (major, rounding) | partly fixed | Sides and closest approaches are recorded and C4b is re-diagnosed correctly. But the scan is read as convergence although the mirror contradicts it, there is no relative-offset or seed study, and the λ⁻¹ln(rest/offset) caveat the log claims is not in D54.md (D3be-2). |
| -3 | fixed | c3_outcome uses the distance. C2's gap-30 tail is argued from R beyond its range, but confirmed: R at x = 30 gives s_d = 3.26·10⁶ > 9.7·10⁵, the tail 0.36σ from feature 0. |
| -4 | fixed | 3.135 (fresh start); test corrected |
| -5 | fixed | FC's 0.01668 |
| -6 | fixed | capture set stated; R at 23 gaps (quantised: D3be-6) |
| -7 | partly fixed | new hand-computed or misquoted numbers (D3be-7) |
| -8 | fixed in form | content incomplete (D3be-1) |
| -9 | partly fixed | 11 new tests; 16 of 29 mutants survive (D3be-8) |
| -10 | fixed | |
| -11 | fixed | characterisations acceptable; one more neighbour suggested (D3be-9) |
| -12 | fixed | new dead outputs (D3be-10) |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3be-1 | major | **The far side ends in a swap: the entrant races in its turn.**<br>- D54.md:108 "(the third-moment drift vanishing only at w = w_e)" is false. E*₃ ∝ ±(A − 2w)/(Aw(A − w))^{1/2}, which is strictly monotone on (0, A). At w = w_e, E*₃ = 0 but ∂_wE*₃ = −1.002·10⁻⁷ (6σ), −7.48·10⁻¹¹ (7σ).<br>- Static check: V's far branch continued upward has dE*/dw < 0 at 77 of 77 steps from w = 0.001 to 0.78, through w_e. At w ≈ 0.4 the secant is −1.025·10⁻⁷ against the law's −1.002·10⁻⁷. The entrant moves +0.049σ → +0.64σ toward A′ as w_e falls to 0.022.<br>- R's mirrored problem, continued past 3·10⁸ at x = 7: w = 0.153 (3·10⁸), 0.283 (6·10⁸), 0.467 (10⁹). The coefficients cross at s = 8.67·10⁸, and the entrant leaves feature 0 at 1.40·10⁹ and fits feature 2. Feature 0 keeps the former runner-up at w = 0.8.<br>- A formed far runner-up does the same: x = 4.9, C2's target, rtol 10⁻¹¹, far by rounding. The coefficients cross at 8.40·10⁸ and the entrant leaves at 1.37·10⁹.<br>- Loop, mirrored target, x = 7, ε = 10⁻¹⁷ (s₂ = 3.07·10⁹): counts [1, 1, 1], kernel 146 (the entrant) captured 0 → 2, rms 4·10⁻¹⁷.<br>- Consequences:<br>  • the far +1 (D54.md:117–118, 68) needs s_θ < s₂ < s_{d,e}; at x = 7 that is 4.8·10⁻¹⁷ < ε < 2.9·10⁻¹²;<br>  • "a formed runner-up stays bound" fails for ε < 5·10⁻¹⁷;<br>  • (iv) (l. 65–69) misses the entrant's capture (−1 at a formed runner-up's feature), and "the reach of all three grows as ε → 0" is false on the far side;<br>  • l. 5 ("grows on the other") and l. 190 ("adds a kernel on the far side") overstate;<br>  • PARADIGM.md:6224 and :160 repeat the error. | **(v):**<br>- State that the pair's member displaced toward the unfitted feature, whichever it is, decays and races; on the far side the runner-up grows through w_e and the entrant, once lighter, folds and leaves at s_{d,e}(x).<br>- Far +1 iff s_θ < s₂ < s_{d,e}.<br>- Drop the w_e parenthesis.<br>- Proof: ∂_wE*₃ ≠ 0 on (0, A).<br>**(iv):** add the exception (and "by one per runner-up"); replace "the reach … grows" with the windows.<br>**R:** continue the mirror to ≈ 3·10⁹ and print s_swap and s_{d,e}.<br>**C4M:** add the 10⁻¹⁷ run.<br>**PARADIGM D‴.51(ii)(c):** "…grows until it outgrows the entrant, which then migrates". |
| D3be-2 | major | **For formed runner-ups the side and s_d are a rounding floor's; the tolerance scan is not convergence.**<br>- D54.md:280 ("at 10⁻¹³ every run takes the side the layer keeps") is contradicted by l. 277–278: on the mirror, the formed x = 3 and 4.9 go toward at the same rtol (RTOL_R = 10⁻¹³).<br>- My traces at rtol 10⁻¹³, x = 4.9:<br>  • the collapse keeps the layer's direction to ≈10⁻¹³σ: e-share +0.48 → +0.50 on C2's target and −0.48 → −0.50 on the mirror, both at 1.3·10⁻¹⁰σ;<br>  • the full variational equation reproduces the computed gap (1.0·10⁻¹⁰ against 1.06·10⁻¹⁰σ at s = 400);<br>  • at s ≈ 450 the x-gap falls to 1 ulp (9·10⁻¹⁷σ) while the y-gap pins: at −1.54·10⁻¹⁴σ on C2's target, at −1.57·10⁻¹⁴σ on the mirror (where it was +3.0·10⁻¹³σ at s = 419, the layer's far side), and at +1.5·10⁻¹⁴σ on the y-reflection;<br>  • the e-share then jumps to +0.87–0.88 on all three targets, so all three go toward.<br>- The floor does not depend on rtol: at x = 4.9 the closest approach is 1.5·10⁻¹⁴, 6.1·10⁻¹⁵, 1.2·10⁻¹⁴, 1.5·10⁻¹⁴ and 1.5·10⁻¹⁴σ at rtol 10⁻⁹ … 10⁻¹³.<br>- It scales as w·floor = 2.6–3.3·10⁻¹⁶σ (1–2 ulps of the centres) for x ≤ 5.48. It is the rounding of two nearby kernels' forces, not "the integration's error" or "the centres' resolution" (l. 110).<br>- 26 of the 74 rule comparisons involve such sides (C3 x ≤ 6.04: 6; C4/C4b x = 3, 4.9: 20). In 9 of them the side decides the outcome.<br>- s_d of a formed runner-up is the climb from the floor. λ_u = 8.2·10⁻⁴, 5.4·10⁻⁴, 3.6·10⁻⁴ at x = 4.9, 5.48, 6.04, i.e. 2.8·10³, 4.3·10³, 6.4·10³ per decade. ln(σ/floor)/λ_u = 5.8·10⁴ and 8.4·10⁴, against s_d = 6.20·10⁴ and 8.68·10⁴ (94%, 97%).<br>- When the floor is hit, the relative mode is still contracting: Jacobian eigenvalues −0.135, −0.105, −0.026 at s = 450, 600, 800, about 34 e-folds more to s = 1000. So the ε-free problem's offset lies far below the floor. If it ends 15 decades lower, its s_d is longer by ~4·10⁴, 6·10⁴ and 10⁵ at x = 4.9, 5.48, 6.04.<br>- Hence "least at x = 5.48" (l. 268), the capture set's lower edge (C3 at x = 5.07: s₂ = 1.64·10⁵ vs R's 1.52·10⁵) and the "ε-free s_d" of formed runner-ups are not established.<br>- The s_d caveat that the round-53 log (D3az-2 row) says was added is absent from D54.md. | **Code:** integrate D = c_r − c_e and the block's centre as state variables, with a cancellation-free relative force (or extended-precision correlations); show the side on C2's target, the mirror and the y-reflection, and s_d against the offset.<br>**Until then:**<br>- (v) l. 109–110: the rule is established for tails only;<br>- l. 278–280: replace the convergence sentence with the mirror counter-example;<br>- R and the 73/74 count (l. 264–275): mark the formed rows (x ≤ 6.04) as floor-limited, with s_d a lower bound;<br>- state s_d's λ_u⁻¹ln(rest/offset) term. |
| D3be-3 | minor | **The linear pull** (D54.md:88–89, 142–148).<br>- "The block's centre of mass kept … an unfitted feature's linear pull alike [cancels]": slaving the entrant's centre moves the block by ⟨r′, ∂_eφ⟩/(A‖∂_eφ‖²). In one dimension that is A′σ(d/σ)e^{−d²/(4σ²)}/A, i.e. 6.5·10⁻⁴σ at 6σ.<br>- R's block sits (+4.10, −2.39)·10⁻⁴σ from feature 0, constant while W grows from 0.15 to 0.80.<br>- The proof's prescriptions m₁ = 0 and "projected off span{φ, ∇φ}" are inconsistent. E* survives, because the O(\|δ\|⁵) cross term ⟨m₁∂φ, M₃∂³φ⟩ is exactly what the projection removes.<br>- δ is then relative to the block's centre, whereas V measures from c* (0.8% of the 6σ, w = 0.3 rest). | Proof: slave c_e by the projection only, and define δ from the block's centre. Statement: "the linear pull moves the block". |
| D3be-4 | minor | **The reduction's justification; far tails.**<br>- (a) l. 157–158 assert the reduction without the separation that makes it valid: the rest's restoring rate (~2λ_u ≈ 10⁻³ per unit s) against the coefficient's relative drift (10⁻⁹–10⁻⁸). It is numerically fine (positions 0.1–0.4%, rates 4.5–13%) but fails at the fold, which is where the departure happens.<br>- (b) l. 86–87 say "collapsed onto the entrant … held by the second moment". This is false for far tails: R's closest approach is 0.056–0.33σ at x = 12.2–24.5 (0.36σ at 30), and the tail sits 0.11–0.36σ from feature 0 at s = 10⁵, below the fold, with no rest. It waits because its mobility ∝ w; "travels the slower" (l. 116) is mostly the slower exit (λ_u ∝ w). | Proof: state the separation and where it holds. (v): restrict "collapsed … held" to x ≲ 8.6 (closest approach ≤ 2.4·10⁻⁴σ). |
| D3be-5 | minor | **C4M tests the far rule only far from its boundary** (l. 261–264: s₂/s_θ = 5.4 and 0.17).<br>- The rule does hold at the boundary. My runs: x = 7 at 2·10⁻¹² (s₂/s_θ = 1.20) gives [2, 1, 1] with w(τ₂) = 0.00871, and at 4.5·10⁻¹² (0.80) gives [1, 1, 1] with 0.00741; x = 10 at 2.5·10⁻¹³ (1.24) gives [2, 1, 1] with 0.01046, and at 6·10⁻¹³ (0.80) gives [1, 1, 1] with 0.00655.<br>- c3_outcome labels the two unformed far tails differently ("not departed" and "one kernel"), by the residual cut. | Add the four runs to C4M; label held tails by their side. |
| D3be-6 | minor | **R's s_d is quantised.**<br>- It is the first output beyond 2σ on a geometric grid of ratio 1.038 (exp_transit_nd.py:1033).<br>- x = 5.07 and 10 share the grid point 1.520·10⁵; x = 17.26 and 17.3 share 5.837·10⁵.<br>- The draft quotes three digits (l. 266–270), and its capture-tie comparison (5.48·10⁵ vs 5.62·10⁵, 2.5%) lies within one grid step. | race_job: find the 2σ crossing by dense output or an event; print the bracket. |
| D3be-7 | minor | **Numbers.**<br>- l. 242 "mean 1.0025": not in the md (1.00242 unrounded).<br>- l. 286–287 "below the unformed corners by up to 6.5–50% per triple": (5.9/5.95/6, 9, 9) print −0.2 … +0.1%, and (1, 1, 1) has no such corners.<br>- l. 255 "captured from 10⁻¹⁰ on … 3.30–3.52·10⁵ … 0.97–1.03": at the paragraph's rtol 10⁻¹², 10⁻¹¹ is bound; 3.30·10⁵ and 1.03 are C4b's.<br>- l. 110 "the centres' floating-point resolution (2·10⁻¹⁵ … 6·10⁻¹⁴σ)": these are 6–150 ulps.<br>- l. 160: the proof cites a referee's trace as evidence.<br>- l. 227: C2's gap-30 tail is argued from R beyond R's range (true, see the D3az-3 row above). | Print the mean and R at x = 30 in the md; fix T3's range and C4's tolerance; record the e-share in R. |
| D3be-8 | minor | **Tests: 16 of 29 mutants survive** (13 killed).<br>- 6 change the md under `--from-raw` (lines changed): T3's tail column also at H/2 (9); `side_of` returning empty (31); τ₂ − τ_e by the mean (7); e3_ratios' displacement at the step's start (10); VF's cell tolerance 0.5 (1); c4m_lines printing d_held (2).<br>- 4 change only a rerun:<br>  • race_job's θ taken absolute (R's s_θ column);<br>  • its s_θ guard removed (s_θ ≈ 10³ for the formed mirror rows);<br>  • w_dep read at the departure output;<br>  • the departure taken at 1.5σ (s_d one grid step earlier at x = 16.96 and 20.56).<br>- 2 are dead (c4_job's e_t2 ×2).<br>- 4 are equivalent on the data.<br>- No test reaches s_θ. | Add a golden test, `report(*raw) == committed md` (4 s; kills all six). Add a mirrored race to 10⁷ at x = 7 asserting s_θ ≈ 5.7·10⁶, and None for a formed runner-up. Check w_dep at the last output before the departure and d_pre < 2σ ≤ d_dep on a dense grid. Add a c3_outcome case at d_held = 1.5σ. Drop or test e_t2. |
| D3be-9 | minor | **Literature.**<br>- The hold's moment structure is classical. M₂ = w(1 + w/w_e)δ² and M₃ = w(1 − w²/w_e²)δ³, odd and vanishing at equal weights, are the moments of an over-specified two-component mixture fitting one component: Chen's n^{−1/4}; Ho–Nguyen 2016; Dwivedi–Ho–Khamaru–Wainwright–Jordan–Yu, Ann. Statist. 2020, where unbalanced weights converge geometrically and balanced ones exponentially slower. It should be cited as the static neighbour.<br>- Pesme–Flammarion's coordinates "are successively activated"; "dormant" is Kunin et al.'s term.<br>- The other characterisations are acceptable: Chizat's conic particle GD [memory: WFR geometry]; Chizat–Bach's 6-particle failure; Boursier–Flammarion; Poon–Peyré's relative positions in d ≥ 2.<br>- "New" (l. 189–190) should drop "the side kept by the inner layer" for formed runner-ups and "adds a kernel on the far side". | Status, l. 176–190 |
| D3be-10 | cosmetic | **Dead outputs and stale fields.**<br>- c4_job's e_t2 and w_t2 are recorded but never reported or tested (the log says the displacement at τ₂ is recorded).<br>- C3's raw still carries w_unbind.<br>- R's main rows carry s_theta equal to s_d for tails, which is meaningless. | Report them or drop them. |

Sources:
- https://jmlr.org/papers/v26/24-1523.html
- https://arxiv.org/pdf/1907.10300
- https://arxiv.org/pdf/1805.09545
- https://proceedings.neurips.cc/paper_files/paper/2023/hash/17a9ab4190289f0e1504bbb98d1d111a-Abstract.html
- https://proceedings.nips.cc/paper_files/paper/2025/hash/06cbd2e81dfbd3bb4cb0abce95b32584-Abstract-Conference.html
- https://epubs.siam.org/doi/10.1137/17M1147822
- https://ar5iv.arxiv.org/html/1810.00828
- https://arxiv.org/pdf/2006.00704

## Lesson

A potential that changes sign at a point does not stop the flow there. The far side's E*₃ vanishes at w = w_e, but its derivative does not. Following the branch through the point where the two kernels' roles exchange — 77 static steps, and R continued past its own horizon — turns "held and grows" into a swap, after which the entrant races. Integrate past the horizon chosen for the outcome you expect.

Round 53's lesson also returns with a twist. A tolerance scan that converges says something about the dynamics only if the converged answer moves when the geometry that is supposed to decide it moves. At rtol 10⁻¹³ the formed runner-up goes toward on the target and on its mirror alike, because the pair is pinned at a rounding floor whose direction belongs to the samples, not to the inner layer; the tighter tolerance only made the artefact reproducible. Mirror the geometry at the converged tolerance before calling a side converged. Then follow the decisive small quantity below the floor with a formulation that resolves it: relative coordinates, extended precision, or the variational equation.
