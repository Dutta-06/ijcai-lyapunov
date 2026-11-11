# Referee report D3bi, round 62: Proposition D‴.54 as revised after round 58

## Summary

The new pin law is right about the pin, but the revised (v) draws three wrong conclusions from it:
- **The continuum.** Without a window edge, the near tie (the inner layer's offset) decides a formed runner-up's side.
- **The rate.** The law's relaxation and growth rate is about twice the true one.
- **s_d.** The window comparisons offered as evidence about s_d are mostly a rescaling of time.

**What checks out**
- **Reproducibility.**
  - `--from-raw` into my own TN_OUT reproduces `results/chaos/transit_nd.md` byte for byte (md5 9c8370b58b6a…).
  - The 38 tests pass (59.7 s).
  - Fresh reruns equal the raw bitwise, all fields except `secs`:
    - R at x = 4.9, 5.4772 and 10;
    - R's tolerance-scan cell (4.9, 10⁻¹¹, s_end 3·10⁶);
    - R's mirror at x = 7 to 3·10⁸, and to 3·10⁹ (the swap and s_{d,e});
    - all 14 RP rows;
    - the C4M loop (7, 10⁻¹⁰).
- **The relative-mode identity is exact.** I re-derived dD/ds = (2Gm ww_e/W)⟨r, ∇²φ(c̄)⟩D + (2(w − w_e)/W)dc̄/ds. It holds to O(|D|³) when dc̄/ds is the midpoint's actual velocity. Λ₀ is missing (D3bi-8).
- **The rest's shift with W is right.** ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩, and it vanishes in the continuum.
  - Independent check: a long-double Newton solve for the rest, differenced in W, gives (−1.2152·10⁻¹⁵, +1.14039·10⁻¹²)σ on Part C's window. The formula gives (−1.2156·10⁻¹⁵, +1.14039·10⁻¹²)σ.
  - With a right edge at 0.8 the two give 2.7203·10⁻¹² and 2.7203·10⁻¹²σ.
- **The ~10⁻¹⁴σ pin of formed runner-ups is the window edge's.** My own tests:
  - **Right edge.** I moved the right edge to 4.92σ from the block and the top edge to 8.33σ. The rest's shift turns to +x, 2.72·10⁻¹²σ. The pins turn to −x and the law matches them:

    | x | measured pin (σ) | law (σ) |
    |---|---|---|
    | 3 | −9.02·10⁻¹⁵ | −9.10·10⁻¹⁵ |
    | 4.9 | −3.69·10⁻¹⁴ | −3.68·10⁻¹⁴ |
  - **Two near edges** (right and top). The pin is their vector sum: (−3.68, −1.48)·10⁻¹⁴σ against the law's (−3.67, −1.50)·10⁻¹⁴σ.
  - **Sample spacing.** At spacings 0.04 and 0.0625 the rest's shift becomes 1.887·10⁻¹² and 5.81·10⁻¹³σ, and the pin follows it:

    | spacing | x | measured pin (σ) | law (σ) |
    |---|---|---|---|
    | 0.04 | 3 | −6.29·10⁻¹⁵ | −6.35·10⁻¹⁵ |
    | 0.04 | 4.9 | −2.59·10⁻¹⁴ | −2.57·10⁻¹⁴ |
    | 0.0625 | 3 | −2.04·10⁻¹⁵ | −1.95·10⁻¹⁵ |
    | 0.0625 | 4.9 | −7.77·10⁻¹⁵ | −7.91·10⁻¹⁵ |
  - **The law predicts, not merely tracks.** At the pin it is independent of its starting value: started from D(400) = 0 instead of the run's D(400), it gives the same values at x = 3 and 4.9.
- **(v)'s static laws and the far side hold.**
  - The second- and third-moment laws hold.
  - E*₃ is strictly monotone: E*₃ ∝ ∓(A − 2w)/(Aw(A − w))^{1/2}, and d/dw[(A − 2w)(w(A − w))^{−1/2}] = −A²/(2(w(A − w))^{3/2}) < 0.
  - The swap and the entrant's race hold.
  - The far rule s_θ < s₂ < s_{d,e} holds at both boundaries. My two new loops sit at s₂ = 0.87 and 1.15 times s_{d,e} (D3bi-9).
- **Numbers.** Every number I checked in (v), the Status, and the (C4), (C4M), (R), (RP), (V)/(VF), (V3), (FC) and (T3) paragraphs matches the md and raw, except those in D3bi-8. Spot checks of A, T and B4 also match.

**What is wrong: majors**
- **D3bi-1. In the continuum the near tie does decide a formed runner-up's side.** D54.md:122–123 says the opposite.
  - **Setup.** On a window whose edges are at least 11.7σ from the block, the rest's shift with W is ≤ 4·10⁻¹⁶σ, i.e. zero for this purpose.
  - **x = 6.04, resolved directly.** Float64 resolves the layer's offset: the minimum gap is 3.7·10⁻¹⁴σ. The runner-up goes toward on C2's target and far on the mirror, which is the tails' rule. The authors' own top-edge-1.5 row at 6.04 shows the same.
  - **x = 3, 4.9, 5.48: the float64 runs cannot decide.** The pairs merge to D = 0 exactly, after which the side belongs to the integrator. At x = 4.9 the same problem went toward with one output grid and far with another.
  - **x = 3, 4.9, 5.48: a linear continuation can.** I continued the relative mode in log-scaled form below the floor. Against the run it agrees to |log ratio| ≤ 0.073 and 0.007 rad over the stretch the run resolves:

    | x | min \|D\| (σ) | share along e at the minimum (C2's target) |
    |---|---|---|
    | 3 | 10^{−54.8} | +0.96 |
    | 4.9 | 10^{−27.3} | +0.79 |
    | 5.48 | 10^{−19.0} | +0.72 |

    The mirror gives the negatives, and both end on ±e. That is again the tails' rule.
  - **So the window's pin overrides a deterministic offset by 5 to 40 decades; it does not fill a void.**
  - **s_d in the continuum is longer.** The continuum's s_d is about 3.6·10⁵ at x = 4.9 (+12%) and about 8.4·10⁴ at 5.48 (+38%). The latter roughly equals the 6.04 value, so "the least s_d, at 5.48" is a property of the window.
- **D3bi-2. The coefficient in the identity is not the relative mode's rate.**
  - **Why.** The block's rest is that of its centroid c_b, not of the midpoint c̄. Whenever D changes, the midpoint moves by ((w_e − w)/(2W))dD/ds. So the relative mode relaxes, and later grows, at Λ₀Gm ww_eW|⟨r, ∇²φ⟩|/(w² + w_e²). That is λ_D·W²/(2(w² + w_e²)), i.e. 0.51–0.53 of what (v) calls the contraction rate.
  - **The full Jacobian confirms it.**

    | case | Jacobian eigenvalues | λ_D·eig H | with the centroid factor | RP's measured growth |
    |---|---|---|---|---|
    | x = 4.9, s = 450 | −0.1347, −0.1352 | −0.2600 | −0.1373 | — |
    | x = 4.9, late, unstable | +0.00082 | +0.00156 | +0.00082 | 0.000829 |
    | x = 7, late, unstable | +0.00018 | +0.00036 | +0.00018 | — |

  - **With the factor, RP's law matches the cases the draft calls poor.**

    | case | committed law / run | law with the factor / run |
    |---|---|---|
    | x = 6.04, D_y at s = 700 | 0.427 | 0.905 |
    | x = 6.04, D_y at s = 1000 | 0.867 | 0.967 |
    | tail x = 7 | off by 27× to 15 000× | within 1–13% |

  - So the poor agreement at 6.04 and for the tail is the law's error, not "the layer's offset still above the pin" (D54.md:318–319).
- **D3bi-3. The s_d conclusion rests on a time rescaling and does not hold for held runner-ups.**
  - **The rescaling.** ⟨·⟩ is a mean over the samples, so moving the top edge changes their number N and rescales s by N. At x = 3, s_d with the top edge at 0.95 and at 1.5 is 0.97557 and 1.24400 times its value at 1.0. The sample-count ratios are 40/41 = 0.97561 and 51/41 = 1.24390.
  - **So the "s_d falling" at the 0.95 edge (D54.md:322) is entirely the rescaling at x = 3.** Beyond the rescaling, s_d moves by −1.3% at 4.9 and −9.0% at 6.04.
  - **The hold dominates for held runner-ups.** The climb from the pin to 10⁻²σ takes 9.9·10³ of s_d = 2.35·10⁷ at x = 3 (0.04%) and 3.46·10⁴ of 3.19·10⁵ at 4.9 (11%). "s_d — the climb from the pin" (D54.md:125) is true only for x ≳ 5.

**Minor**
- **D3bi-4. Window direction decides formed runners' outcomes.** A right edge 4.9σ away turns C4's (4.9, 10⁻¹⁰) from captured into bound on the far side. So the formed rows behind R's rule and the capture set's lower edge belong to Part C's top edge.
- **D3bi-5. (iv)'s exception is incomplete.**
  - It omits the −1 of a formed far-side runner-up when s₂ > s_{d,e}.
  - It does not say what decides a formed runner-up's side.
- **D3bi-6. "Held by the block's second moment" overreaches.** It is said of every collapsed gap (x ≲ 8.6), but 5.48 ≤ x ≤ 8.59 lie below the fold and have no hold.
- **D3bi-7. The pins are only resolved to the state's quantisation.** That is 1.85·10⁻¹⁶σ, so x = 3's pin is 23 ulps. The long-double run changes only the sums.
- **D3bi-8. Numbers and text:**
  - "2.36·10⁷" is the old grid value;
  - Λ₀ is missing from the identity;
  - the proof ends with a garbled duplicate;
  - s_θ is still read on the output grid.
- **D3bi-9. C4M's ε = 10⁻¹⁷ row.** Its "closest approach 0.11σ" is not the inner layer's. The far rule's s_{d,e} boundary was untested, and my loops confirm it.
- **D3bi-10. Mutation testing.** 21 of 36 mutants survive, and 14 of the survivors change rerun results. The golden test protects the report, not the job functions that produce the raw.
- **D3bi-11. Literature.**
  - Conic/Wasserstein–Fisher–Rao particle flows move positions independently of the weight, so they have neither the trailing nor the pin.
  - Clustered-spike super-resolution is the hold's static neighbour.
- **D3bi-12 (cosmetic).** The C4 and C4b raw rows still carry `e_t2` and `w_t2`.

## What I checked and how

**Setup.**
- I reviewed 86fbc02 through `git archive` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref62/D3bi_work/repo`. Scripts are in `…/D3bi_work/scripts`, outputs in `…/D3bi_work/out`, mutant copies in `…/D3bi_work/mut/w1` and `…/mut/w2`.
- Every run used the venv python with OMP_NUM_THREADS=1, nice -n 10 and PYTHONDONTWRITEBYTECODE=1, with at most two heavy processes of mine at once.
- I wrote nothing to the repository and signalled no process.
- The live HEAD has since moved to e025861. `git diff 86fbc02 e025861` is empty for all D54 files. The live tree has uncommitted edits to `src/core/rich_transit.py` and to the close-features files; they are not mine.

**Regeneration, tests and reruns.**
- `cmp` reports the regenerated md identical.
- pytest: 38 passed.
- Bitwise reruns: `rerun_cells.py` for R and RP, `rerun_c4m.py` for C4M (7, 10⁻¹⁰).

**The identity, the rest and the rate.**
- I re-derived the identity and its centroid form: with c_b slaved,

  dD/ds = (W/(w² + w_e²))·[(w − w_e)dc*/ds + Λ₀Gm ww_e⟨r, ∇²φ⟩D].

  The pin is unchanged, D* = −(w − w_e)(dc*/ds)/(Λ₀Gm ww_e⟨r, ∇²φ⟩), and the rate is ×W²/(2(w² + w_e²)).
- `jac_rate.py` takes the eigenvalues of Λ₀·`critical_jac_nd` along R's race at x = 4.9 and 7.
- `rest_fd.py` solves for the rest by long-double Newton at W = 0.70–0.90 and differences in W.

**The pin law.**
- `pinlaw.py` replicates `rp_job`; its "committed" variant reproduces RP's law values exactly. It adds variants: started from D(400) = 0, and with rate and forcing multiplied by W²/(2(w² + w_e²)).
- It also takes configurable windows: a right edge at 0.8; edges at 0.8 and 1.0; spacing 0.04 and 0.0625; far edges x ∈ [−2, 2], y ∈ [−1.5, 2].
- Each window was run on C2's target and on the mirror.

**The continuum's side.**
- `continuum_side3/4.py` work on the far-edge window.
- From a resolved point (|D| ≈ 10⁻⁹σ) they continue D in log-scaled form with the centroid-factor law along the float64 base trajectory.
- They are run without the block's forcing and with it: the measured centroid velocity while resolved, then extrapolated exponentially at the fitted rate 0.028–0.049.
- I validated the continuation against the run over the stretch where the run still resolves D.

**Loops.**
- `c4_window.py`: C4 at (4.9, 10⁻¹⁰) on the window with the right edge at 0.8 and the top at 1.5.
- `c4m_trace.py`: C4M on the mirror at x = 7 and ε = 6.5·10⁻¹⁷, 3.7·10⁻¹⁷ and 10⁻¹⁷, recording the pair's gap, distances and coefficients after the entry. Each took about 110–145 s.

**Mutants.**
- `mutants.py` and `mutdriver.py`: 36 single replacements in `exp_transit_nd.py`. Each was run against `pytest -x` and against `--from-raw` with the md diffed.
- Survivors with no md change were also rerun on the affected cells (`cellcheck.py`): R at 4.9, R's long mirror run at 7, RP at 4.9 and 6.04, and c4_job.

**Literature.** WebSearch. arxiv.org and the ar5iv hosts (ar5iv.labs.arxiv.org, ar5iv.arxiv.org) did not resolve, so the mixture and WFR points rest on publisher pages, abstracts and search summaries.

## Verdicts per part

| claim | verdict |
|---|---|
| (i)–(iii); (iv)'s count and codimension law | unchanged since round 53; spot-checked numbers (A, T, B4, T3) match |
| (iv) exception: −1 toward, +1 capture tie, far +1 iff s_θ < s₂ < s_{d,e}, the entrant's −1 | holds. C4M covers the s_θ boundary; my loops cover s_{d,e}. Omits the formed far-side −1 (D3bi-5) |
| (v) M₂ and M₃ laws; E* − E₀ = −P²/(2B); E*₃ monotone on (0, A) | holds (re-derived) |
| fold c(d, A′/A); no fold on the far side; swap at s_swap; entrant leaves at s_{d,e} | holds numerically; R's long runs rerun bitwise |
| "held by the block's second moment" for x ≲ 8.6 | overreach: 5.48–8.59 collapse without a hold (D3bi-6) |
| relative-mode identity, D54.md:116 | holds exactly (Λ₀ missing) |
| "λ_D the contraction rate", D54.md:121 | wrong by W²/(2(w² + w_e²)) ≈ 0.51–0.53 (D3bi-2) |
| ∂c̄*/∂W formula; zero in the continuum | holds (independent finite difference agrees to 4 digits) |
| the formed runner-ups' pin is the window edge's | holds (three own discriminating tests) |
| pin law 0.902, 0.988, 1.005 at x = 3, 4.9, 5.48 | holds, and is a genuine prediction; the x = 3 value is within about 2 ulps of the state |
| "poor at 6.04 and for the tail: the layer's offset", D54.md:318–319 | wrong explanation; the law's rate (D3bi-2) |
| a tail's side is the layer's | holds |
| "in the continuum … nothing in the near tie decides a formed runner-up's side", D54.md:122–123 | wrong (D3bi-1) |
| a formed runner-up's side belongs to the sampled problem | holds on samples; a right edge 4.9σ away flips it (D3bi-4) |
| "s_d — the climb from the pin — belong to the sampled problem", D54.md:125, :322 | wrong for x ≤ 4.9; the window comparison is a rescaling of time (D3bi-3) |
| s_d least at 5.48; capture set's lower edge "the pin's" | true of Part C's window only; the continuum differs (D3bi-1, D3bi-4) |
| R's rule, 73 of 74 | reproduced from the raw; the 26 formed rows test the side the window imposes |
| Status: "the relative mode's law exact to first order and its pin checked"; "New: … a formed runner-up's [side] is set by whatever moves the block" | the pin yes; the law's rate no; the "New" claim holds only for perturbations above 10⁻¹⁹ to 10⁻⁵⁵σ |

**D3be issues**

| D3be | status | note |
|---|---|---|
| -1 (major, far swap) | fixed | R's long runs (bitwise), C4M at 10⁻¹⁷, text, PARADIGM. My loops at s₂/s_{d,e} = 0.87 and 1.15 give [2, 1, 1] and [1, 1, 1] as predicted. (iv) still omits the formed far-side −1 (D3bi-5) |
| -2 (major, the floor) | partly fixed | The floor is the window's pin, and my tests confirm the law. But the continuum's answer was asserted rather than computed, and it is wrong (D3bi-1). Also the rate (D3bi-2) and the s_d attribution (D3bi-3). The long double tests only the sums (D3bi-7) |
| -3 | fixed | |
| -4 | fixed | wording overreach (D3bi-6) |
| -5 | fixed | |
| -6 | mostly fixed | s_d, s_swap and s_{d,e} are bisected. s_θ is still on the grid, and the draft still quotes 2.36·10⁷ |
| -7 | partly fixed | new misquote and text errors (D3bi-8) |
| -8 | partly fixed | The golden test kills all 7 md-changing mutants. 21 of 36 survive, 14 of them changing reruns (D3bi-10) |
| -9 | fixed | further neighbours (D3bi-11) |
| -10 | mostly fixed | C4 and C4b raw still carry `e_t2` and `w_t2` |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bi-1 | major | **In the continuum a formed runner-up's side is the near tie's.**<br>- D54.md:122–123 says that in the continuum "nothing in the near tie decides a formed runner-up's side". Without forcing, the relative mode is linear and never reaches zero, so the layer's offset keeps its sign.<br>- **Far-edge window** (edges ≥ 11.7σ away; rest's shift with W ≤ 4·10⁻¹⁶σ): x = 6.04 is resolved in float64. The minimum gap is 3.7·10⁻¹⁴σ, with share along e +0.63 on C2's target and −0.64 on the mirror; the runner-up goes toward and far respectively. The authors' top-1.5 row at 6.04 also goes toward.<br>- At x = 3, 4.9 and 5.48 the float64 pairs merge to D = 0 exactly. After that the side is the integrator's: x = 4.9 went toward with a 1200-point output grid and far with a 6000-point one, and the mirror did the opposite.<br>- **Linear continuation below the floor**, validated to \|log ratio\| ≤ 0.073 and angle ≤ 0.007 rad. On C2's target:<br>  • x = 3: min\|D\| = 10^{−54.8}σ with the forcing extrapolated (10^{−102.7}σ homogeneous), share at the minimum +0.96;<br>  • x = 4.9: 10^{−27.3}σ, +0.79;<br>  • x = 5.48: 10^{−19.0}σ, +0.72.<br>  The mirror gives the negatives, and all end at ±0.99 along e. That is the tails' rule, toward iff (c_r⁰ − c_e⁰)·e > 0.<br>- So Part C's pin (10⁻¹⁵ to 10⁻¹⁴σ) overrides the near tie's offset by 5 to 40 decades.<br>- **Consequences for s_d** (estimates from the continuation and the growth rate):<br>  • x = 4.9: s_d ≈ 3.6·10⁵, i.e. +12%;<br>  • x = 5.48: \|D\| reaches 10⁻²σ at s ≈ 7.45·10⁴ instead of 5.13·10⁴, so s_d ≈ 8.4·10⁴ (+38%). That is about the 6.04 value (8.65·10⁴ in the far window, rescaled), so "the least at 5.48" (D54.md:299) belongs to the window.<br>  • C3's capture at x = 5.07 (s_d 1.50·10⁵ against s₂ 1.64·10⁵) likely becomes bound in the continuum. This is an estimate. | **(v), l. 122–125:** "In the continuum the forcing vanishes and the relative mode keeps the layer's sign, as a tail's does. For formed runner-ups the offset falls to 10^{−19} (x = 5.48) … 10^{−55}σ (x = 3) before it grows, so any shift of the block's rest with W above that decides instead: on samples the window's edge, in floating point without one the rounding."<br>**Status/New, l. 215–217:** restrict "set by whatever moves the block" accordingly.<br>**RP:** add a far-edge row and the relative-mode continuation in log-scaled form, which is the variational treatment D3be-2 asked for.<br>**R:** mark s_d at x ≥ 5 as the window's. |
| D3bi-2 | major | **The relative mode's rate is mislabelled, and RP's law integrates with it.**<br>- The identity (D54.md:116) is exact, but its dc̄/ds is the midpoint's velocity.<br>- The block's rest belongs to the centroid c_b, where F(c_b) = 0, not to the midpoint; the proof's "the block's rest solves F(c̄) = 0" (l. 180) conflates the two. With c_b slaved, the midpoint moves by ((w_e − w)/(2W))dD/ds.<br>- Hence the effective rate is λ_D·W²/(2(w² + w_e²)) = 0.508 to 0.530 of λ_D here, and "λ_D the contraction rate" (l. 121) is wrong.<br>- **Full Jacobian at x = 4.9:**<br>  • s = 450: −0.1347 and −0.1352, against λ_D·eig H = −0.2600, −0.2590 and the corrected −0.1373, −0.1368;<br>  • s = 600: −0.1052, −0.1061, against −0.2025 and the corrected −0.1069;<br>  • late: −0.00041 and +0.00082, against −0.00078 and +0.00156. RP's measured growth is 0.000829.<br>- At x = 7 the late unstable eigenvalue is +0.00018, against λ_D's +0.00036.<br>- **`rp_job` (exp_transit_nd.py:1217, 1223)** treats (∂c̄*/∂W)dW/ds as external. With both terms multiplied by the factor, the law tracks the run:<br>  • x = 6.04: D_y ratio 0.905 and 0.967 at s = 700 and 1000, against the committed 0.427 and 0.867; at s = 500, D = (6.82·10⁻⁹, −1.74·10⁻¹⁰) against the run's (7.06·10⁻⁹, −1.80·10⁻¹⁰), while the committed law gives 1.77·10⁻¹¹;<br>  • tail x = 7: D_x 1.009·10⁻⁵, 7.81·10⁻⁸ and 1.72·10⁻⁸ against the run's 1.019·10⁻⁵, 8.46·10⁻⁸ and 1.98·10⁻⁸. The committed law is off by 27×, 3200× and 15 000×;<br>  • x = 5.48 at s = 500: the D_y ratio is 0.934, against the committed 0.092.<br>- The pin's value is unaffected. Forcing and rate scale alike, and at the pin the two laws agree to ≤ 1%.<br>- So "at x = 6.04 the layer's offset is still above the pin" (l. 318–319) explains the run, not the law's failure. | **(v):** state the slaved-centroid form above, with rate Λ₀Gm ww_eW\|⟨r, ∇²φ⟩\|/(w² + w_e²) and the same pin. Add Λ₀ to l. 116.<br>**Proof l. 177–181:** derive the pin and the rate from it; distinguish c_b from c̄.<br>**`rp_job`:** multiply rate and forcing by W²/(2(w² + w_e²)). Report and quote the law at 6.04 and 7, and compare its rate with the Jacobian's eigenvalue.<br>**Test:** law/run at 6.04, s = 1000, within 10%. |
| D3bi-3 | major | **The s_d attribution: a time rescaling and the hold.**<br>- ⟨·⟩ is the mean over the samples, so the top-edge windows (41×40, 41×41, 41×51 samples) rescale every rate, and with it s, by N/1681.<br>- At x = 3, s_d at the 0.95 and 1.5 edges is 2.28912·10⁷ and 2.91898·10⁷ against 2.34645·10⁷. The ratios 0.97557 and 1.24400 match 40/41 = 0.97561 and 51/41 = 1.24390.<br>- The growth rate at x = 3 and top 1.5 is 0.002625/0.003257 = 0.806 of the top-1.0 value, against 41/51 = 0.804.<br>- So "s_d falling from 2.346·10⁷, 3.191·10⁵, 8.439·10⁴ to 2.289·10⁷, 3.073·10⁵, 7.489·10⁴" (l. 322) is pure rescaling at x = 3. Beyond the rescaling the change is −1.3% at 4.9 and −9.0% at 6.04. At top 1.5, x = 6.04 is +2.5% beyond 51/41.<br>- On Part C's window, \|D\| reaches 10⁻²σ at s = 9.9·10³, 3.46·10⁴, 5.13·10⁴ and 7.27·10⁴ for x = 3, 4.9, 5.48 and 6.04. Against s_d = 2.35·10⁷, 3.19·10⁵, 6.10·10⁴ and 8.44·10⁴, the climb is 0.04%, 11%, 84% and 86% of s_d.<br>- So "s_d — the climb from the pin" (l. 125) is false for held runner-ups (x ≤ 4.9), whose s_d is the hold's and pin-independent. | **(v), l. 124–125:** "the side belongs to the sampled problem; s_d does only where there is no hold (x ≳ 5), through the climb".<br>**RP:** compare s_d in rescaled time (s·N₀/N), or hold N fixed by moving the bottom edge with the top. State the rescaling in the md. |
| D3bi-4 | minor | **Which edge is nearest decides formed runner-ups' outcomes.**<br>- Window x ∈ [−1, 0.8], y ∈ [−1, 1.5]: the right edge 4.92σ from the block, the top 8.33σ.<br>  • The formed runner-ups go to the far side on C2's target and do not depart by 5·10⁷.<br>  • They go toward on the mirror, with s_d = 2.64·10⁷ and 3.58·10⁵. Divided by 1887/1681, the latter is 3.19·10⁵, Part C's value.<br>- C4's loop at (4.9, 10⁻¹⁰) on this window ends bound, on the far side, with closest approach 3.09·10⁻¹⁴σ, w_max 0.0221 and counts [2, 1, 1]. On Part C's window it is captured.<br>- So the formed rows behind R's rule (26 of 74: C3 at x ≤ 6.04, C4/C4b at x = 3 and 4.9), the captures in (iv)'s −1 and "its lower edge the pin's" (l. 131) are properties of Part C's top edge. The draft says so only implicitly. | **R/C3/C4 text:** state that formed rows test the side the window forces.<br>**RP:** add the right-edge row.<br>**C4:** add one loop on that window. |
| D3bi-5 | minor | **(iv)'s exception** (D54.md:65–70):<br>- It omits that a formed far-side runner-up loses its entrant when s₂ > s_{d,e}: −1, back to one. (v) l. 134 has it.<br>- It does not say how a formed runner-up's side is decided (D3bi-1, D3bi-4).<br>- "and closes" is vague: at fixed x the window in ε is the interval (9.7/s_{d,e})² < ε < (9.7/s_θ)², which is 4.9·10⁻¹⁷ to 2.9·10⁻¹² at x = 7. | Add the formed far-side −1. Name the side's decider for formed runner-ups. Replace "closes" with the interval. |
| D3bi-6 | minor | D54.md:88–89 says that wherever the layer collapses the pair (x ≲ 8.6) "it is held by the block's second moment". R shows 5.48 ≤ x ≤ 8.59 leaving "without a hold": below the fold there is no rest. | "collapsed for x ≲ 8.6; held (a rest above the fold) only for x ≲ 5.07". |
| D3bi-7 | minor | **Resolution of the pins.**<br>- D is the difference of float64 centres, so D_y is quantised in steps of ulp(c_y)/σ = 1.85·10⁻¹⁶σ.<br>- x = 3's pin (−4.26·10⁻¹⁵σ) is 23 steps, and the law's 0.902 is 2.3 steps off; RP's x components (±4.6·10⁻¹⁷ to 9·10⁻¹⁶σ) are 0.5 to 10 steps.<br>- The long-double rows change only the sums; the state is float64, so they do not test this. | State the quantisation. Quote x = 3's agreement as "within 2 ulps". |
| D3bi-8 | minor | **Numbers and text:**<br>- D54.md:281 "R's s_d = 2.36·10⁷" is the old grid value (s_d_grid 2.3618·10⁷); the bisected value is 2.346·10⁷.<br>- The identity in (v) and its proof (l. 116, 177) lack Λ₀. (iii) l. 42 and 146 have it, and the md's law has it.<br>- l. 181 ends with a garbled duplicate: "… (V, VF, V3, FC, R, RP). ∎race lie beyond the expansion (V, VF, V3, FC, R). ∎".<br>- s_θ is still read on the output grid (exp_transit_nd.py:1070; ratio 1.038 for the 3·10⁸ runs and 1.044 for the 3·10⁹ runs) but quoted to 3 digits. x = 10's s_θ is 1.566·10⁷ in the main mirror run and 1.552·10⁷ in the long run. | Fix the text; bisect s_θ or quote it to 2 digits. |
| D3bi-9 | minor | **C4M at 10⁻¹⁷, and the far rule's upper boundary.**<br>- The md's "closest approach 1.1e-01" for (7, 10⁻¹⁷) (transit_nd.md:356) is the swap-phase minimum: 0.109σ at s = 8.2·10⁸. The dense grid starts at s = 307, by which point the pair is already 0.207σ apart: the runner-up is 0.205σ from feature 0 with w 0.00494.<br>- R reaches its closest approach, 1.8·10⁻⁸σ, at s ≈ 1340, and is still at 4.7·10⁻⁴σ at 10⁴.<br>- At 6.5·10⁻¹⁷ the first output (s = 120) already shows the collapsed pair at 2.0·10⁻⁷σ, growing at 1.8·10⁻⁴ per unit s.<br>- So at ε ≲ 10⁻¹⁶ the loop's early phase precedes the grid or runs ahead of R's s. Both loops agree with R from s ≈ 10⁶ on; at 10⁷ the distance is 0.362 against 0.364 and w is 0.0107 against 0.0105.<br>- **Far rule, newly tested at s_{d,e}** (C4M had only s₂ = 2.2·s_{d,e}):<br>  • ε = 6.5·10⁻¹⁷ (s₂ = 1.203·10⁹ = 0.87·s_{d,e}): [2, 1, 1]. The swap falls between outputs at 8.9·10⁸; R's is 8.63·10⁸. The entrant stays, with w = 0.219 against the former runner-up's 0.581.<br>  • ε = 3.7·10⁻¹⁷ (1.15·s_{d,e}): [1, 1, 1]. The entrant is captured 0 → 2 and is beyond 2σ by s = 1.44·10⁹; R's s_{d,e} is 1.385·10⁹. | Start C4M's dense grid at the loop's own entry. Mark the 10⁻¹⁷ closest-approach cell as not the layer's. Add the two boundary loops to C4M. |
| D3bi-10 | minor | **Mutation testing: 36 single-point mutants; 15 killed, 21 survive.**<br>- All 7 md-changing mutants are killed by the golden test.<br>- 7 survivors are equivalent on the data (crossing returns its left end; entrant's departure at 1σ; side hard-coded to "far" in `c4m_outcome`; Q2's off-diagonals; least over all of `post`; two-side selection at 0.005; R's side at SIDE_R/2).<br>- **14 survivors change rerun results but never the md**, because `--from-raw` does not run the jobs:<br>  • crossing's bisection inverted: s_d −3.1%, s_swap +3.8%, s_{d,e} +4.0%. The test at tests/test_transit_nd.py:417 only brackets one grid step;<br>  • two bisection steps: +0.13%; s_d at 1.9σ: 2·10⁻⁵;<br>  • e_closest unnormalised (1.3·10⁻¹⁴ for 0.876) or read one output late;<br>  • swap at 0.9·w_e (s_swap −4.8%) or the swap crossing's sign (+3.8%);<br>  • e_dest by argmax: "ends at feature 1";<br>  • `rp_job`'s forcing ×1.03 (inside the 3% test at :547), its start s₀ ×1.2, or D(s₀) = 0 (both change only the 6.04 and 7 rows);<br>  • `rp_job`'s s_d at 1σ (−0.17%), growth window from 10³ (−2.1%), e_least unnormalised.<br>- Killed by other tests: s_θ's restriction and guard, `c4m_outcome`'s w_max, `PairSums`' kernel swap, the block's H and ⟨φ, ∇φ⟩ ×1.02, the rate halved, and c4_job's side at 2·SIDE_R. | Golden job tests (cheap):<br>- `race_job(('near tie', 4.9))`: s_d = 3.19103·10⁵ ± 10⁻⁶ relative, e_closest 0.876 ± 0.01 (2 s);<br>- the long mirror run at x = 7: s_swap and s_{d,e} to 10⁻⁶, e_dest = 2 (13 s);<br>- `rp_job` at 4.9 and 6.04: law at all marks to 1%;<br>- `c4m_outcome` with side "toward". |
| D3bi-11 | minor | **Literature.**<br>- Chizat's conic particle gradient descent is called the free law's "nearest dynamical relative" (l. 201–203). It couples a multiplicative weight update with an additive position update, and its mean-field limit is the Wasserstein–Fisher–Rao flow. There a particle's position velocity does not scale with its weight, so the forcing (2(w − w_e)/W)dc̄/ds is absent: no trailing and no pin, and a formed runner-up's side would always be the layer's. The same holds for Yan–Wang–Rigollet's WFR particle flow for Gaussian mixtures (Ann. Statist. 2024).<br>- So "the lighter kernel trails the block" (New, l. 215–216) is an elementary consequence of the paradigm's Euclidean metric on (w, c) (mobility ∝ w), not of particle flows in general.<br>- The window-edge pin rests on the truncated Gaussian's ⟨φ, ∇_cφ⟩ ≠ 0. I found no prior use of it as a symmetry breaker; it is an artefact of the sampled problem.<br>- The over-specified mixtures are correctly cited. Dwivedi et al., Ann. Statist. 48(6) 2020: unbalanced weights converge geometrically to O((d/n)^{1/2}); balanced ones exponentially slower to O((d/n)^{1/4}).<br>- Suggested static neighbours for the hold's moments about the cluster's barycentre: Denoyelle–Duval–Peyré (JFAA 2017), clustered positive spikes; Batenkov–Goldman–Yomdin (Inf. Inference), near-colliding sources.<br>- No prior statement of the far-side swap found.<br>- Caveat: full texts were unreachable (arxiv and ar5iv); the WFR position-mobility point rests on the summaries and on the standard form of the WFR flow. | Status, l. 201–217: state the metric difference, add the two super-resolution references, and recast the "New" item as a consequence of mobility ∝ w. |
| D3bi-12 | cosmetic | The C4 and C4b raw rows still carry `e_t2` and `w_t2`. The round-58 log reports them as "dropped"; only the code dropped them. | Strip them as C3's `w_unbind` was. |

Sources:
- https://projecteuclid.org/journals/annals-of-statistics/volume-52/issue-4/Learning-Gaussian-mixtures-using-the-WassersteinFisherRao-gradient-flow/10.1214/24-AOS2416.full
- https://www.doi.org/10.1214/19-AOS1924
- https://deepai.org/publication/singularity-misspecification-and-the-convergence-rate-of-em
- https://arxiv.org/pdf/1907.10300 (search-result summary only; arxiv.org unreachable)
- https://arxiv.org/pdf/2312.05993 (search-result summary of conic PGD's updates)
- https://proceedings.neurips.cc/paper/2018/hash/a1afc58c6ca9540d057299ec3016d726-Abstract.html
- https://hal.archives-ouvertes.fr/hal-01271269
- https://arxiv.org/pdf/1904.09186 (search-result summary)

## Lesson

A law that reproduces a pin says what the pin is, not what happens without it. Removing the forcing (a window whose edges are 12σ away) and continuing the same law below the floating-point floor, in log-scaled form, shows the near tie still decides, at 10⁻¹⁹ to 10⁻⁵⁵σ. The window only overrides it. Two habits would have caught the rest:
- **Check a rate against the Jacobian before integrating with it.** A coefficient read off an exact identity need not be the rate: here the midpoint in the forcing moves with D, which halves the rate.
- **Change one thing at a time.** Moving a window's edge also changed the number of samples that ⟨·⟩ averages over. That rescaled time and produced the whole s_d change at x = 3.
