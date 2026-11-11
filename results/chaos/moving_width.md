# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)

The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = 0.12 (`rich_transit_nd.BumpsND`); 21 kernels started at −1, −0.9, …, 1; 321 samples of [-1.6, 1.6]; exact: as many formed kernels (|w| above 1% of the smallest amplitude) as features, each feature's centre, width and amplitude within 1e-06, the residual's rms below 1e-08.

## A. One feature

A = 0.8 at c* = 0.03, the width σ* = ratio·σ₀; moving widths: the kernels, the largest error of (c, σ, w), the residual; fixed widths (up to 4 entries, 3e+05 evaluations — † the budget spent first): the kernels and the residual:

| σ*/σ₀ | moving: kernels | exact | max error (c, σ, w) | residual | fixed: kernels | residual |
|---|---|---|---|---|---|---|
| 0.25 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 6.8e-17 | 2† | 7.49e-02 |
| 0.5 | 1 | yes | 0.0e+00, 1.2e-16, 1.4e-15 | 2.3e-16 | 2† | 6.51e-02 |
| 0.75 | 1 | yes | 3.5e-18, 1.4e-17, 5.6e-16 | 1.1e-16 | 2† | 3.57e-02 |
| 1.031 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 1.5e-16 | 4 | 1.30e-06 |
| 1.118 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 1.5e-16 | 3† | 1.10e-03 |
| 1.414 | 1 | yes | 1.4e-17, 0.0e+00, 5.6e-16 | 1.7e-16 | 4 | 2.19e-03 |
| 2.236 | 1 | yes | 3.5e-18, 0.0e+00, 5.6e-16 | 2.1e-16 | 4 | 2.90e-02 |
| 3 | 1 | yes | 3.5e-17, 0.0e+00, 5.6e-16 | 2.4e-16 | 4 | 7.17e-02 |
| 4 | 1 | yes | 1.7e-17, 0.0e+00, 5.6e-16 | 2.8e-16 | 4 | 1.35e-01 |

### A2. One kernel from a large coefficient (referee rounds 74, 80, 83)

A's feature at σ*/σ₀ = 0.25, 1, 4: one kernel at (c*, σ*) with w = 2.5, 3, 5, 10, 20·A — E(start) (1 − w/A)²·½⟨f*²⟩, (i)'s condition E(start) < ½‖f*‖² holding for w < 2A only —, its block's flow on the critical curve to s = 1e+13: at the end w/A, σ/σ*, (c − c*)/σ* and E over ½⟨f*²⟩; the smallest σ/σ* on the way (at the outputs); the plateau — the outputs with E/½⟨f*²⟩ within 0.001 of E₂, the kernel interpolating the two samples nearest c* (below) —, its E over ½⟨f*²⟩, and the first output after it with σ/σ* back above 0.5 (the escape):

| σ*/σ₀ | w(start)/A | E(start)/½⟨f*²⟩ | end: w/A, σ/σ*, (c − c*)/σ*, E/½⟨f*²⟩ | smallest σ/σ* | on the plateau: s, E/½⟨f*²⟩ | escaped at s | status |
|---|---|---|---|---|---|---|---|
| 0.25 | 2.5 | 2.25 | 1.0000, 1, +1.2e-16, 2.51e-32 (s = 1e+13) | 0.304 | — | — | 0 |
| 0.25 | 3 | 4 | 1.0000, 1, +0.0e+00, 2.23e-32 (s = 1e+13) | 0.111 | 1e+00 … 3e+05, 0.644–0.644 | 1e+06 | 0 |
| 0.25 | 5 | 16 | 1.0000, 1, +0.0e+00, 2.23e-32 (s = 1e+13) | 0.0921 | 1e+00 … 3e+07, 0.644–0.644 | 1e+08 | 0 |
| 0.25 | 10 | 81 | 1.0000, 1, +0.0e+00, 2.23e-32 (s = 1e+13) | 0.0772 | 1e+00 … 3e+10, 0.644–0.644 | 1e+11 | 0 |
| 0.25 | 20 | 361 | 18.5905, 0.06861, +1.7e-01, 6.44e-01 (s = 1e+13) | 0.0678 | 1e+00 … 1e+13, 0.644–0.644 | not by the end | 0 |
| 1 | 2.5 | 2.25 | 1.0000, 1, -2.9e-17, 2.62e-32 (s = 1e+13) | 0.295 | — | — | 0 |
| 1 | 3 | 4 | 1.0000, 1, +5.8e-17, 2.76e-32 (s = 1e+13) | 0.0281 | 1e+00 … 3e+05, 0.906–0.906 | 1e+06 | 0 |
| 1 | 5 | 16 | 1.0000, 1, -2.9e-17, 2.62e-32 (s = 1e+13) | 0.0232 | 1e+00 … 3e+07, 0.906–0.906 | 1e+08 | 0 |
| 1 | 10 | 81 | 1.0000, 1, +2.9e-17, 2.49e-32 (s = 1e+13) | 0.0194 | 1e+00 … 1e+10, 0.906–0.906 | 3e+10 | 0 |
| 1 | 20 | 361 | 17.5957, 0.01739, +4.2e-02, 9.06e-01 (s = 1e+13) | 0.017 | 1e+00 … 1e+13, 0.906–0.906 | not by the end | 0 |
| 4 | 2.5 | 2.25 | 1.0000, 1, +9.4e-17, 2.63e-32 (s = 1e+13) | 0.302 | — | — | 0 |
| 4 | 3 | 4 | 1.0000, 1, +1.1e-16, 2.63e-32 (s = 1e+13) | 0.00704 | 1e+00 … 3e+05, 0.976–0.976 | 1e+06 | 0 |
| 4 | 5 | 16 | 1.0000, 1, +8.0e-17, 2.49e-32 (s = 1e+13) | 0.00581 | 1e+00 … 3e+07, 0.976–0.976 | 1e+08 | 0 |
| 4 | 10 | 81 | 1.0000, 1, +8.7e-17, 2.49e-32 (s = 1e+13) | 0.00485 | 1e+00 … 1e+10, 0.976–0.976 | 3e+10 | 0 |
| 4 | 20 | 361 | 17.4894, 0.004354, +1.0e-02, 9.76e-01 (s = 1e+13) | 0.00426 | 1e+00 … 1e+13, 0.976–0.976 | not by the end | 0 |

The plateau (referee round 83): a kernel narrower than the samples' spacing fits at most two neighbouring samples, so lim_{σ→0} min_c J_N interpolates the best pair, E₂ = 1 − max_k (f*(x_k)² + f*(x_{k+1})²)/Σf*² — here the two samples nearest c* —, and J_N's σ → 0 limit at a fixed c, on the best sample, is E₁ = 1 − max_k f*(x_k)²/Σf*²; over ½⟨f*²⟩ at σ*/σ₀ = 0.25, 1, 4: E₂ 0.643650, 0.906294, 0.976497, E₁ 0.8119, 0.9530, 0.9882.

The symmetric line (referee rounds 97, 100): on samples exactly symmetric about c* the line c = c* is invariant (Lemma D‴.59a: the flow commutes with the reflection c ↦ 2c* − c) — the flow's |ċ| on it (w = 3A, 5A; σ/σ* = 1, 0.5, 0.2, 0.1) at most 4.2e-15, 4.2e-15, 2.1e-15 there, roundoff —; A's samples are symmetric about c* but for 6 at the window's far end, whose f*/A is at most 10^{-602.3}, 10^{-37.6}, 10^{-2.4} (0.0e+00, 2.3e-38, 4.4e-03 in double) at σ*/σ₀ = 0.25, 1, 4, and in double the samples' and f*'s roundoff (the centre sample lies 8 ulps of c* above c*) — and where their f* is not negligible the unpaired samples — drift c on the line: |ċ| there at most 1.5e-12, 7.2e-14, 2.9e-04. The problem as posed is a perturbation of the symmetric one, and the symmetric state, for σ well below h, a saddle in c (Lemma D‴.59a): from c* exactly the flow leaves it too, after a time set by the asymmetry's logarithm over the unstable rate, and a solver can stay on it by damping the unstable direction (LSODA's steps; below). The samples XS are symmetric about 0 to 4.4e-16 in double (referee round 106). A2 starts off the line (c* + 1e-09σ*): the generic flow.

The early collapse (referee round 83): from 5A at σ*/σ₀ = 0.25 — c* = 0.03 is a sample — the block's flow (`a2_job`'s, LSODA) over the span (0, S), S = 10¹³ being A2's (a terminal event past s = 10 leaves its steps alone): (c − c*)/σ*, σ/σ* and E/½⟨f*²⟩ at s ≤ 10, and the first output on the plateau (E within 0.001 of E₂):

| S | s = 1e-04 | s = 1e-03 | s = 1e-02 | s = 1e-01 | s = 1e+00 | s = 1e+01 | on the plateau from s |
|---|---|---|---|---|---|---|---|
| 1e+01 | +1.0e-15, 0.8626, 13.327 | +6.9e-16, 0.1846, 3.480 | -1.2e-16, 0.1840, 3.480 | +9.3e-16, 0.1840, 3.480 | +9.3e-16, 0.1841, 3.473 | +9.3e-16, 0.1846, 3.406 | not by s = 10 |
| 1e+13 | +0.0e+00, 0.8626, 13.327 | -1.0e-14, 0.1846, 3.480 | -1.7e-01, 0.0921, 0.644 | -1.7e-01, 0.0921, 0.644 | -1.7e-01, 0.0921, 0.644 | -1.7e-01, 0.0921, 0.644 | 0.0056 |

On the symmetric state (\|c − c*\| below 10⁻¹²σ*) the flow's ∂ċ/∂c is -3.833e+02 … -2.789e+02 while the width falls (σ/σ* 0.3440 … 0.8626, 8 outputs: stable in c; the (α, σ) block's eigenvalues there -886 … -644 and +0.013 … +0.022, the width contracting faster than c (∂ċ/∂c above the smallest at every output): not normally hyperbolic, referee rounds 125, 130) and +1.050e+04 … +1.089e+04 at σ/σ* 0.1840 … 0.1846 (17 outputs): a saddle in c, an offset growing ×e per 9.2e-05 … 9.5e-05 of s, left once the seed — on A's samples the data's drift, on exactly symmetric ones the sums' rounding (below) — has grown.
LSODA's steps on the symmetric state (s in (0.001, 0.01], the offset below 10⁻⁶σ*; λ = 1.089e+04), S = 1e+01: Adams at hλ 0.14–0.42, BDF at hλ 2.4–50 from s = 0.00264; the offset's largest 2.2e-11σ*, at the last such step +9.3e-16σ*; over the Adams steps above 10⁻¹²σ* it grows at 1.101e+04 per unit s; its BDF order-1 steps at hλ 3.99–7.56 multiply the offset by -0.335 … -0.152 (1/(1 − hλ): -0.335 … -0.153); its BDF order-2 steps at hλ 2.41–2.41 multiply the offset by -2.54 … -1.67; 1 further BDF step left out, where the offset falls below 1e-14σ* — to the scale of the drift's fixed point −d/λ (below; referee round 100).
LSODA's steps on the symmetric state (s in (0.001, 0.01], the offset below 10⁻⁶σ*; λ = 1.089e+04), S = 1e+13: Adams at hλ 0.14–0.55; the offset's largest 8.9e-07σ*, at the last such step -8.9e-07σ*; over the Adams steps above 10⁻¹²σ* it grows at 1.092e+04 per unit s.
Lemma D‴.59a's split of ∂ċ/∂c = Gm·w·⟨r∂²_cφ − w(∂_cφ)²⟩ at its largest output (s = 0.0018, σ/σ* 0.1840, the samples' spacing 1.81σ): the centre sample +15703 — (Gm w/(Nσ²))(w − f*(c*)), the kernel overshooting the sample it sits on —, the samples at ±1 -4919.8, at ±2 +105.6, the rest +0.052; their sum 10888.72207 against the exact Jacobian's 10888.72207.
(BDF-k damps a positive real eigenvalue only for hλ beyond 2, 4, 6.667, 10.67, 17.07 (k = 1 … 5; referee round 87): below that an implicit step amplifies the unstable offset; on the Adams steps it grows at λ, the flow's own instability.)
From 3A at σ*/σ₀ = 0.25 with the centre at c* exactly (A2's start before referee round 87), over the span 1e+13: |c − c*|/σ* at most 6.3e-13 at the outputs to s = 10, σ/σ* down to 0.2345, not on the plateau by then — held on the symmetric state by LSODA's damping (A2 starts off it, at c* + 1e-09σ*; referee round 106: the clause from the data). Its steps to s = 10 (referee round 100): Adams to s = 0.0117, then BDF at hλ 6.96–9.44 (order 1), 0.696–111 (order 2), 295–1.24e+03 (order 3), λ the outputs' ∂ċ/∂c; the offset at most 9.4e-13σ*, at the last step +1.04e-15σ*.
The held runs at s = 10 (referee round 100): 5A over the span 1e+01: the line's drift ċ at c = c* -9.84e-12σ* per unit s, λ = 1.05e+04 (the (α, σ) block's eigenvalues -0.0012 and -1.446e+04: λ above both, normally hyperbolic, for every r as both are negative; referee rounds 125, 130), the implicit step's fixed point −d/λ = +9.372e-16σ* against the held offset +9.252e-16σ*; 3A over the span 1e+13: the line's drift ċ at c = c* -6.34e-13σ* per unit s, λ = 628.9 (the (α, σ) block's eigenvalues -0.00062 and -2875: λ above both, normally hyperbolic, for every r as both are negative; referee rounds 125, 130), the implicit step's fixed point −d/λ = +1.008e-15σ* against the held offset +1.041e-15σ* — the offset sits within an ulp of the zero of ċ (−d/λ with d summed in double carries the double sum's error in d), the computed problem's own near-symmetric state (referee rounds 103, 106, 109: 5A: held 8 ulps of c* above c*, the zero of ċ at 8.12 ulps in 50-digit arithmetic on the double data (8.1 in double), the centre sample 8 ulps above c*; 3A: held 9 ulps of c* above c* (its last 10 steps to s = 10 at 9, 10 ulps), the zero of ċ at 9.57 ulps in 50-digit arithmetic on the double data (9.9 in double), the centre sample 8 ulps above c*): a repelling state the flow leaves, which LSODA's steps (BDF beyond its bounds) make attracting — the location the problem's, the stability LSODA's step sequence's (referee round 106: scipy's Radau and BDF need not hold it, below).
The flow from c* exactly (referee rounds 103, 106, 109): |c − c*| reaches 10⁻²σ* (the exit's direction in brackets) at s = 3A: on A's samples 0.0483 (−), 0.0478 (−) (DOP853, rtol 1e-10, 1e-12), 0.0567 (−) (Radau), held to s = 1 (BDF), 0.0516 (−) (BDF at rtol 1e-12) (scipy, rtol 1e-10 unless stated), 0.0481 (−) with the samples' sums exactly rounded, 0.0480 (−) summed in reversed order; on exactly symmetric samples from their centre 0.0507 (+), 0.0514 (+) (DOP853), no exit by s = 1 exactly rounded, 0.0507 (−) reversed; on them started 9.57 ulps of c* below their centre — the posed problem's c* below its zero of ċ — 0.0495 (−); 5A: on A's samples 0.0035 (−), 0.0035 (−) (DOP853, rtol 1e-10, 1e-12), 0.0039 (−) (Radau), 0.0039 (−) (BDF), 0.0041 (−) (BDF at rtol 1e-12) (scipy, rtol 1e-10 unless stated), 0.0035 (−) with the samples' sums exactly rounded, 0.0035 (−) summed in reversed order; on exactly symmetric samples from their centre 0.0040 (+), 0.0039 (+) (DOP853), no exit by s = 1 exactly rounded, 0.0040 (−) reversed; on them started 8.12 ulps of c* below their centre — the posed problem's c* below its zero of ċ — 0.0035 (−). In long double with the offset c − c* as the state (`a2_exit_ld`; referee round 112): 3A: on A's samples 0.048297 (−) at rtol 1e-10 and 1e-12 (apart by 2.1e-06 relatively), on the symmetric samples 0.0609 (−), 0.0643 (−); the data's drift carried along the route (`a2_seed`: S = ∫e^{−Λ}ċ|_{c*}ds to the exit, Λ = ∫∂ċ/∂c|_{c*}) S = -7.762e-17 (-2.59e-15σ*), 2.34 times c*'s distance below the zero of ċ (9.57 ulps), e^{Λ(s_e)}|S| 1.004 times the exit's level (the linear response), and the symmetric problem started at c* + S leaving at 0.048297 (−), A's time within 3·10⁻⁵ (a 10⁻³ change of S; the floor of a few 10⁻⁶ is the long double's own rounding of ċ, which no tolerance moves — referee round 125); 5A: on A's samples 0.0034820 (−) at rtol 1e-10 and 1e-12 (apart by 8.0e-07 relatively), on the symmetric samples 0.00452 (−), 0.00457 (−); the data's drift carried along the route (`a2_seed`: S = ∫e^{−Λ}ċ|_{c*}ds to the exit, Λ = ∫∂ċ/∂c|_{c*}) S = -5.473e-17 (-1.82e-15σ*), 1.94 times c*'s distance below the zero of ċ (8.12 ulps), e^{Λ(s_e)}|S| 1.003 times the exit's level (the linear response), and the symmetric problem started at c* + S leaving at 0.0034820 (−), A's time within 3·10⁻⁵ (a 10⁻³ change of S; the floor of a few 10⁻⁶ is the long double's own rounding of ċ, which no tolerance moves — referee round 125) — on A's samples every exit downward, in double its time within 1.1% of the long-double one whatever the tolerance and the sums' order or rounding (DOP853), in long double the same at both tolerances to 10⁻⁴, and the symmetric problem started at the data's seed leaving with A's direction at A's time to 10⁻⁴: the exit is the data's; on exactly symmetric samples from their centre exact sums keep c at the centre (no exit by s = 1) and reversed sums mirror the run — identities of the symmetric data, not tests —, in the code's order of the sums the arithmetic flips it at both tolerances (3A: +, + in double and −, − in long double at rtol 1e-10, 1e-12; 5A: +, + in double and −, − in long double at rtol 1e-10, 1e-12), and in long double the exit comes 13–26% later — by Λ = ∫λ ds 5.91 … 8.59 along the route (`a2_sym_lag`, from the double exit to the long-double one at each tolerance: the two seeds' ratio e^Λ) —, each in the code's order one draw (A2P, below): its seed is the arithmetic's.
The symmetric samples' exit as a draw (A2P; referee rounds 125, 130, 133, 137): over summation orders of the same data — the exactly symmetric samples in 100 random permutations from each of the seeds 131 and 132 (drawn in round 130, after round 125's rule was set: round 125's stored draws, seed 125, were its referee's own) and in the code's order, one draw apart (its Λ within 7.8e-03 of `a2_sym_lag`'s, along the long-double run's own route) —, DOP853 from their centre in double and in long double at its own steps, each exit's effective seed |S| = 10⁻²σ*·e^{−Λ(0, s_e)} (Λ along the line c = c*, `a2_line_lam`) in units of uσ*, u the unit roundoff (2^{-53} in double, 2^{-64} in long double), and per order the paired Λ = ln(|S_d|/|S_ld|), from the double exit to the long-double one, and the exponent Λ/(11 ln 2) with its 99% interval (mean ± 2.576 sem): 3A at rtol 1e-10: Λ 8.33 ± 0.20 and 7.68 ± 0.20 (seeds 131 and 132, over 100 and 100: +3.6 and +0.3 sem from 11 ln 2), pooled 8.01 ± 0.14 (over 200, +2.7 sem; median 7.97, interquartile 6.83 … 9.28, extremes -0.19 … 13.58 of 200), the exponent 1.050 (99% interval 1.001 … 1.098), the seeds' geometric mean 1.43 (interquartile 0.66–3.7) u_dσ* in double and 0.98 (interquartile 0.44–2.9) u_ldσ* in long double, upward in 95 of 200 (double) and 104 of 200 (long double), the long double later by a median 23% (interquartile 19–27%, extremes -0.5% … 40.8% of 200), DOP853's 60 and 73 steps to the exit on average (mean step 8.71e-04 and 8.83e-04); 3A at rtol 1e-12: Λ 8.06 ± 0.24 and 7.45 ± 0.23 (seeds 131 and 132, over 100 and 100: +1.8 and -0.8 sem from 11 ln 2), pooled 7.76 ± 0.17 (over 200, +0.8 sem; median 7.69, interquartile 6.60 … 8.94, extremes 1.57 … 16.27 of 200), the exponent 1.017 (99% interval 0.961 … 1.073), the seeds' geometric mean 0.96 (interquartile 0.39–2.9) u_dσ* in double and 0.84 (interquartile 0.38–2.8) u_ldσ* in long double, upward in 104 of 200 (double) and 103 of 200 (long double), the long double later by a median 22% (interquartile 18–27%, extremes 3.8% … 47.3% of 200), DOP853's 79 and 100 steps to the exit on average (mean step 6.75e-04 and 6.51e-04); 5A at rtol 1e-10: Λ 7.57 ± 0.17 and 7.45 ± 0.19 (seeds 131 and 132, over 100 and 100: -0.3 and -0.9 sem from 11 ln 2), pooled 7.51 ± 0.13 (over 200, -0.9 sem; median 7.75, interquartile 6.33 … 8.70, extremes 0.57 … 12.58 of 200), the exponent 0.985 (99% interval 0.942 … 1.029), the seeds' geometric mean 0.60 (interquartile 0.21–1.9) u_dσ* in double and 0.67 (interquartile 0.34–1.5) u_ldσ* in long double, upward in 93 of 200 (double) and 92 of 200 (long double), the long double later by a median 19% (interquartile 15–22%, extremes 1.2% … 32.3% of 200), DOP853's 55 and 66 steps to the exit on average (mean step 6.85e-05 and 6.76e-05); 5A at rtol 1e-12: Λ 7.50 ± 0.17 and 7.74 ± 0.16 (seeds 131 and 132, over 100 and 100: -0.8 and +0.7 sem from 11 ln 2), pooled 7.62 ± 0.12 (over 200, -0.1 sem; median 7.72, interquartile 6.46 … 8.63, extremes 0.27 … 12.22 of 200), the exponent 0.999 (99% interval 0.960 … 1.038), the seeds' geometric mean 0.60 (interquartile 0.24–1.7) u_dσ* in double and 0.60 (interquartile 0.29–1.2) u_ldσ* in long double, upward in 105 of 200 (double) and 103 of 200 (long double), the long double later by a median 19% (interquartile 15–21%, extremes 0.6% … 30.6% of 200), DOP853's 74 and 93 steps to the exit on average (mean step 5.13e-05 and 4.81e-05) — against ln(u_d/u_ld) = 11 ln 2 = 7.62, the unit roundoffs' ratio; over all 800 draws Λ = 7.723 ± 0.070, the exponent Λ/(11 ln 2) = 1.013 ± 0.009 (± 0.0098 over the 200 orders, each order's cells averaged: the same orders recur in every cell), the exits upward in 397 of 800 (50%) in double and 402 of 800 (50%) in long double, the per-seed means' largest offset +3.6 sem (3A at rtol 1e-10, seed 131); the rule on all of each cell's draws (A2P with A2R's further seeds: 3A at rtol 1e-10 over 600, 1.015 … 1.073, 3A at rtol 1e-12 over 200, 0.961 … 1.073, 5A at rtol 1e-10 over 200, 0.942 … 1.029, 5A at rtol 1e-12 over 200, 0.960 … 1.038) with the band [0.75, 1.25] — halfway to the neighbouring exponents ½ and 3/2; at γ = 1 the two prefactors' ratio within e^{±1.91}, while each geometric mean within 0.05–20 uσ* asks each prefactor O(1); round 133's band [0.9, 1.1], set after the data (e^{±0.76}): A2P's 200 alone passed at 3A, rtol 1e-10 by 0.0017 at its upper end, against the interval's half-width 0.0484, marginal —: every cell's interval inside the band and every geometric mean within it, so the seed scales with the arithmetic's unit roundoff at DOP853's steps (A2H, below) — the two seeds' ratio the unit roundoffs' to within the band, the exponent and the prefactors' ratio not separable from two arithmetics (A2T, below: three) —, while the direction, Λ and the lateness of one order are draws; Λ = 11 ln 2 itself asks the two arithmetics' prefactors to be equal — no cell beyond 3 sem of it (the largest +2.7, 3A at rtol 1e-10: A2R, below).
The steps (A2H; referee rounds 130, 133, 137): on A2P's first 100 random orders of seed 131 at rtol 1e-10, DOP853 with its step at most 0.25·h₀ and 0.0625·h₀ in both arithmetics, h₀ the double's mean step over A2P's draws at that rtol (the line's Λ at its own steps), ln(S/uσ*) as mean ± sem and the shift of its mean from DOP853's own steps on the same orders (± the two sems in quadrature: the orders' draws are independent at different steps, their correlation -0.10 … +0.24), against Brouwer's law for rounding in the derivative, S ∝ h^{1/2} (a shift ½ ln f at the cap f·h₀): 3A (h₀ 8.71e-04): at DOP853's own steps mean step 8.701e-04 and 8.866e-04 (double, long double; 60 and 73 steps), ln(S/uσ*) +0.41 ± 0.13 and -0.29 ± 0.15, Λ 8.33 ± 0.20 (+3.6 sem from 11 ln 2, the exponent 1.093 ± 0.026); at most 0.25·h₀: mean step 2.161e-04 and 2.164e-04 (246 and 300 steps), ln(S/uσ*) -0.13 ± 0.19 and -0.21 ± 0.21, its mean moved by -0.54 ± 0.23 and +0.08 ± 0.26 (Brouwer's -0.69: +0.66 and +2.96 sem from it), Λ 7.71 ± 0.28 (+0.3 sem from 11 ln 2, the exponent 1.011 ± 0.037), upward in 47 and 49 of 100; at most 0.0625·h₀: mean step 5.443e-05 and 5.443e-05 (989 and 1195 steps), ln(S/uσ*) -0.52 ± 0.10 and -0.36 ± 0.07, its mean moved by -0.94 ± 0.16 and -0.07 ± 0.17 (Brouwer's -1.39: +2.78 and +7.93 sem from it), Λ 7.46 ± 0.12 (-1.4 sem from 11 ln 2, the exponent 0.979 ± 0.015), upward in 54 and 46 of 100; 5A (h₀ 6.85e-05): at DOP853's own steps mean step 6.846e-05 and 6.767e-05 (double, long double; 55 and 66 steps), ln(S/uσ*) -0.53 ± 0.15 and -0.48 ± 0.10, Λ 7.57 ± 0.17 (-0.3 sem from 11 ln 2, the exponent 0.993 ± 0.022); at most 0.25·h₀: mean step 1.708e-05 and 1.708e-05 (225 and 265 steps), ln(S/uσ*) -1.17 ± 0.08 and -1.01 ± 0.10, its mean moved by -0.64 ± 0.17 and -0.53 ± 0.15 (Brouwer's -0.69: +0.31 and +1.10 sem from it), Λ 7.46 ± 0.14 (-1.2 sem from 11 ln 2, the exponent 0.979 ± 0.018), upward in 51 and 57 of 100; at most 0.0625·h₀: mean step 4.276e-06 and 4.277e-06 (901 and 1073 steps), ln(S/uσ*) -1.24 ± 0.07 and -1.61 ± 0.04, its mean moved by -0.70 ± 0.16 and -1.13 ± 0.11 (Brouwer's -1.39: +4.15 and +2.26 sem from it), Λ 8.00 ± 0.08 (+4.7 sem from 11 ln 2, the exponent 1.049 ± 0.011), upward in 57 and 47 of 100 — the seed's prefactor in units of uσ* moves with the step sequence (the shifts of the means -1.13 … +0.08 in ln, 5 of 8 beyond 3 sem), the integrator sampling the rounding, 2 of 8 beyond 3 sem of Brouwer's law (3A at 0.0625·h₀ in long double, 5A at 0.0625·h₀ in double); at DOP853's own steps the two arithmetics' mean steps differ by at most 6% (A2P's cells), under a cap both step at it, and there Λ is 7.46 … 8.00 (the exponent 0.979 … 1.049, 1 of 4 beyond 3 sem of 11 ln 2: 5A at 0.0625·h₀, +4.7 — the two arithmetics' prefactors differ at equal steps too), the step sequence moving the two arithmetics' prefactors in part alike; over every configuration (4 cells at DOP853's own steps, A2P; 4 caps) the exponent Λ/(11 ln 2) spans 0.979 … 1.050 — a range of 8 estimates of sems 0.011–0.037; their heterogeneity Q = 27.8 on 7 dof, the between-configuration sd τ = 0.029 (DerSimonian–Laird), the random-effects exponent 1.009 ± 0.012 over configurations —, beside the pooled exponent's statistical ± 0.009 (A2P).
The replications (A2R; referee rounds 133, 137): at 3A, rtol 1e-10 — selected as A2P's cell farthest from 11 ln 2 (of 4, on its 200 draws: the selection's) — 100 random orders from each of the further seeds 133, 134, 135, 136: Λ 7.78 ± 0.22, 7.97 ± 0.19, 8.05 ± 0.21, 7.96 ± 0.23, over their 400 alone 7.941 ± 0.107 (+2.97 sem); with the selected 200, over all 600 Λ = 7.963 ± 0.086 (+3.95 sem) from 11 ln 2, the exponent 1.044 (99% interval 1.015 … 1.073), ln(S/uσ*) +0.26 ± 0.06 in double and -0.08 ± 0.06 in long double — their difference Λ − 11 ln 2 = +0.34 ± 0.09: on the further draws alone within 3 sem (+2.97), not resolved by this replication alone; with the double run at the long double's atol (3.0e-16 = 1e-14·σ*, against 1e-14) on the same 600 orders (seeds 131, 132, 133, 134, 135, 136) Λ = 7.928 ± 0.091 (+3.33 sem), on the further seeds' 400 alone 7.959 ± 0.106 (+3.16 sem), DOP853's double taking 67.67 steps against 60.28 at its own atol (the long double 73.08), its mean step 7.783e-04 against 8.730e-04 (the long double's 8.841e-04); Λ's paired shift from the double's own atol on the same 600 orders (the long double's runs the same) -0.035 ± 0.092: not the atol; A2H's cap at 0.0625·h₀ at 5A — selected as A2H's farthest from 11 ln 2 (of 4, +4.65 sem on seed 131: the selection's) — on seed 132's first 100 orders alone: Λ 7.825 ± 0.077 (+2.59 sem), ln(S/uσ*) -1.31 ± 0.07 in double and -1.51 ± 0.03 in long double; with the selected seed's, over both seeds' 200 Λ = 7.912 ± 0.056 (+5.13 sem), the exponent 1.038 (99% interval 1.019 … 1.057): on the further draws alone within 3 sem (+2.59), not replicated beyond 3 sem by itself; the 2 replications alone combined (Stouffer, Σt/√2) +3.93 sem: the two arithmetics' runs' prefactors differ, on the further draws alone.
The computed field at the scale of the rounding and three arithmetics at equal steps (A2T; referee round 137): the exactly symmetric samples with the centre's offset δ as the state and the same right-hand side in float32, double, long double (`_DTBumps`: at double and long double bitwise `BumpsCW`'s and `_LDBumps`') — (a) ċ(δ) at the long-double line flow's states, on 801 points of |δ| ≤ 4 ulp(x₁) (x₁ = 0.01, the ±1 samples' coordinate, whose x₁ − δ rounds at half an ulp; ulp(x₁) = 0.521 uσ* in every arithmetic), in each arithmetic on the code's order and seed 141's first two random orders, against λδ (λ from ±10⁻⁶σ*): 3A (s = 0.02, 0.04, 0.06: σ/σ* 0.2345, 0.2345, 0.2345; λ 661.7 … 704.5): attracting zeros in 27 of 27 scans, the nearest to 0 at -1.495 … -0.505 and +0.505 … +1.495 ulp, by position (in how many scans) -2.505 (9), -1.495 (20), -0.505 (24), +0.505 (24), +1.495 (21), +2.505 (11), the slope inside ±0.45 ulp 4.11 … 4.39λ, the rms of ċ − λδ 0.61 … 0.98 λuσ*; 5A (s = 0.001, 0.002, 0.004: σ/σ* 0.1846, 0.1840, 0.1840; λ 1.071e+04 … 1.093e+04): attracting zeros in 0 of 27 scans, the slope inside ±0.45 ulp 1.44 … 1.46λ, the rms of ċ − λδ 0.07 … 0.09 λuσ* — (b) classical RK4 at the fixed step f·h₀ (h₀ the double's mean DOP853 step at rtol 1e-10, A2P), the same step sequence in every arithmetic and the same orders in each, the exit at |δ| = 10⁻²σ* by s = 0.5 and its effective seed as A2P's (Λ along the line at rtol 1e-12), exits and upward ones per arithmetic (float32, double, long double): 3A (h₀ 8.713e-04): at 1·h₀ (8.713e-04, 30 orders of seed 142): exits 30, 30, 30 of 30 by s = 0.5; float32 s_e 0.0245, upward 12, ln(S/uσ*) -1.13 ± 0.22 (median -1.09); double s_e 0.0550, upward 19, ln(S/uσ*) -1.33 ± 0.24 (median -1.03); long double s_e 0.0663, upward 10, ln(S/uσ*) -1.22 ± 0.24 (median -0.90); over the 30 orders where all three leave the three-point slope γ = 1.0047 ± 0.0115 (99% interval 0.9751 … 1.0343, median 0.9998), pairwise float32/double 1.0103 ± 0.0131, float32/long double 1.0033 ± 0.0133, double/long double 0.9846 ± 0.0531, the prefactors' log-ratios ln(C_float32/C_long double) +0.09 ± 0.37, ln(C_double/C_long double) -0.12 ± 0.40; at 0.5·h₀ (4.356e-04, 30 orders of seed 142): exits 30, 30, 30 of 30 by s = 0.5; float32 s_e 0.0285, upward 17, ln(S/uσ*) -3.77 ± 0.73 (median -2.77); double s_e 0.0614, upward 14, ln(S/uσ*) -5.53 ± 0.67 (median -4.88); long double s_e 0.0710, upward 16, ln(S/uσ*) -4.33 ± 0.61 (median -3.91); over the 30 orders where all three leave the three-point slope γ = 1.0341 ± 0.0315 (99% interval 0.9530 … 1.1153, median 1.0531), pairwise float32/double 1.0877 ± 0.0472, float32/long double 1.0204 ± 0.0327, double/long double 0.8429 ± 0.1278, the prefactors' log-ratios ln(C_float32/C_long double) +0.57 ± 0.91, ln(C_double/C_long double) -1.20 ± 0.97; at 0.25·h₀ (2.178e-04, 30 orders of seed 142): exits 30, 0, 0 of 30 by s = 0.5; float32 s_e 0.1124, upward 14, ln(S/uσ*) -59.52 ± 8.21 (median -37.90); double —; long double —; 5A (h₀ 6.846e-05): at 0.25·h₀ (1.712e-05, 200 orders of seed 143): exits 200, 199, 200 of 200 by s = 0.5; float32 s_e 0.0021, upward 107, ln(S/uσ*) -2.58 ± 0.06 (median -2.30); double s_e 0.0040, upward 112, ln(S/uσ*) -2.39 ± 0.06 (median -2.06); long double s_e 0.0047, upward 92, ln(S/uσ*) -2.45 ± 0.06 (median -2.11); over the 199 orders where all three leave the three-point slope γ = 0.9939 ± 0.0028 (99% interval 0.9867 … 1.0010, median 0.9979), pairwise float32/double 0.9907 ± 0.0043, float32/long double 0.9947 ± 0.0028, double/long double 1.0050 ± 0.0110, the prefactors' log-ratios ln(C_float32/C_long double) -0.15 ± 0.08, ln(C_double/C_long double) +0.04 ± 0.08; at 0.0625·h₀ (4.279e-06, 100 orders of seed 144): exits 100, 100, 100 of 100 by s = 0.5; float32 s_e 0.0022, upward 59, ln(S/uσ*) -3.19 ± 0.13 (median -2.73); double s_e 0.0040, upward 47, ln(S/uσ*) -3.05 ± 0.12 (median -2.73); long double s_e 0.0047, upward 53, ln(S/uσ*) -2.99 ± 0.11 (median -2.61); over the 100 orders where all three leave the three-point slope γ = 0.9925 ± 0.0060 (99% interval 0.9770 … 1.0081, median 0.9994), pairwise float32/double 0.9928 ± 0.0081, float32/long double 0.9925 ± 0.0061, double/long double 0.9917 ± 0.0185, the prefactors' log-ratios ln(C_float32/C_long double) -0.21 ± 0.17, ln(C_double/C_long double) -0.06 ± 0.14 — at 3A the computed field's ulp-scale zeros decide: RK4 at 0.25·h₀ holds the line in double and long double to s = 0.5, at 1·h₀ and 0.5·h₀ it leaves: the exit's existence and its seed are the integrator's, S ∝ u only through ulp(x₁); where the scans find no attracting zero (5A) three arithmetics at equal steps put γ's 99% interval inside [0.75, 1.25] at every arm, within 0.023 of 1: the seed scales with the unit roundoff, its prefactors a separate O(1) scatter (log-ratios up to 0.21 in size).
Lemma D‴.59a's terms on exactly symmetric samples at σ*/σ₀ = 0.25 against σ/h (h the spacing; the centre, the ±1 samples, the total; referee round 103): w = 4A: σ/h 1: +2871, -2817, -213.7, σ/h 0.8: +4486, -4600, +193.7, σ/h 0.6: +7975, -3794, +4341, σ/h 0.45: +1.418e+04, +575.8, +1.476e+04, σ/h 0.3: +3.19e+04, +759.4, +3.266e+04, σ/h 0.2: +7.178e+04, +4.048, +7.178e+04; the least-squares w: σ/h 1: +109.6, -316.8, -74.77 (w/A 1.3415), σ/h 0.8: +184, -369.4, -20.96 (w/A 1.3616), σ/h 0.6: +274.7, +82.15, +410.9 (w/A 1.3145), σ/h 0.45: +194.4, +655.8, +852.3 (w/A 1.1439), σ/h 0.3: +19.5, +196.3, +215.8 (w/A 1.0073), σ/h 0.2: +0.04217, +1.012, +1.054 (w/A 1.0000). On a grid of σ/h step 0.01 from 1.2 down to 0.15 (`a2_lemma_signs`; referee round 106) the σ/h below which the total, the ±1 term, and the ±1 term above the centre's hold positive (each a root, by bisection between the grid's points; referee round 109): w = 4A: 0.8437, 0.4696, —; the least-squares w: 0.7746, 0.6195, 0.5583 — the scan's states, not A2E's.

## B. The one-kernel reduced energy (the continuum)

J(c, σ) = ½‖f*‖² − ½⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² on a 241 × 241 grid (c within ±6σ*, σ from σ*/8 to 8σ*, A = 1, c* = 0): its grid-local minima and the closed-form gradient's largest relative error against central differences (20 random points):

| σ*/σ₀ | local minima (c, σ/σ*) | min J | gradient error |
|---|---|---|---|
| 0.25 | (+0.0e+00, 1.000) | 0.0e+00 | 4.6e-10 |
| 0.5 | (+0.0e+00, 1.000) | 6.9e-18 | 2.2e-09 |
| 0.75 | (+0.0e+00, 1.000) | 0.0e+00 | 1.3e-09 |
| 1.031 | (+0.0e+00, 1.000) | 0.0e+00 | 6.9e-10 |
| 1.118 | (+0.0e+00, 1.000) | 0.0e+00 | 8.9e-10 |
| 1.414 | (+0.0e+00, 1.000) | -2.8e-17 | 1.3e-09 |
| 2.236 | (+0.0e+00, 1.000) | 5.6e-17 | 1.4e-09 |
| 3 | (+0.0e+00, 1.000) | 0.0e+00 | 1.0e-09 |
| 4 | (+0.0e+00, 1.000) | 0.0e+00 | 7.9e-10 |

## C. Random three-feature targets

40 targets per minimal separation (seed 590; widths log-uniform in [0.05, 0.3], amplitudes ±U(0.4, 1); a moving-width path stopped by its budget counted at its stalled block — G — not at the block before it, referee round 68):

| Δ ≥ | moving widths: exact | formed kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median [range]) | residual (median [max]) |
|---|---|---|---|---|---|
| 0.4 | 36 of 40 | 3 [2, 4] | 3.7e-16 [1.2e-01] | 4 [4, 4] (0† stopped by the budget) | 8.9e-02 [1.9e-01] |
| 0.6 | 35 of 40 | 3 [2, 4] | 3.8e-16 [1.4e-01] | 4 [4, 4] (0† stopped by the budget) | 1.1e-01 [2.1e-01] |

The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; exact with a null cluster where each feature's nearest kernel is exact and the others sum to nothing; a coinciding pair within 0.01 in |Δc|/σ + |Δlog σ| whose net coefficient and centre of mass fit a feature within 0.01, a coinciding pair off the features where it fits none (referee round 106); a budget stop a dipole run-off if the stalled block has an opposite-sign pair under 0.05 with |w| beyond 10 times the largest amplitude): Δ ≥ 0.4: 36 exact, 2 coinciding pair, 2 dipole run-off; Δ ≥ 0.6: 35 exact, 1 coinciding pair, 4 dipole run-off.

Not exact with moving widths (the end's state; for a budget stop the state before the stalled block, and the stalled block): Δ ≥ 0.4, target 1 (widths 0.064, 0.273, 0.079; centres -0.06, +0.35, -0.61; amplitudes -0.85, +0.78, +0.64): coinciding pair, 4 kernels (the pair: +0.334/+0.444, Δc 9.57e-04, the width split +7.9453e-08 (predicted +7.9449e-08)), residual 2.5e-13; Δ ≥ 0.4, target 7 (widths 0.064, 0.234, 0.077; centres +0.27, -0.68, -0.23; amplitudes +0.58, -0.91, -0.43): coinciding pair, 4 kernels (the pair: -0.786/+0.357, Δc 1.53e-04, the width split -1.3404e-07 (predicted -1.3405e-07)), residual 1.5e-13; Δ ≥ 0.4, target 9 (widths 0.141, 0.176, 0.144; centres -0.61, +0.18, +0.65; amplitudes +0.79, -0.85, +0.56): dipole run-off, 2 kernels, residual 2.7e-01; the stalled block: 2 kernels, max|w| 3910 (the accepted polish's, solver-dependent), its closest opposite pair 3.12e-04 apart; Δ ≥ 0.4, target 18 (widths 0.102, 0.280, 0.185; centres -0.58, -0.02, +0.54; amplitudes +0.72, -0.51, +0.60): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 3964 (the accepted polish's, solver-dependent), its closest opposite pair 2.30e-04 apart; Δ ≥ 0.6, target 2 (widths 0.062, 0.295, 0.071; centres +0.69, +0.07, -0.70; amplitudes +0.43, +0.82, +0.47): coinciding pair, 4 kernels (the pair: +0.231/+0.587, Δc 2.78e-08, the width split +3.2834e-08 (predicted +1.8958e-16)), residual 1.0e-15; Δ ≥ 0.6, target 12 (widths 0.219, 0.274, 0.133; centres +0.04, -0.67, +0.68; amplitudes +0.83, -0.49, -0.46): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 1821 (the accepted polish's, solver-dependent), its closest opposite pair 6.13e-04 apart; Δ ≥ 0.6, target 20 (widths 0.080, 0.123, 0.272; centres -0.63, +0.01, +0.66; amplitudes +0.72, +0.83, -0.51): dipole run-off, 3 kernels, residual 1.6e-01; the stalled block: 3 kernels, max|w| 1431 (the accepted polish's, solver-dependent), its closest opposite pair 4.08e-04 apart; Δ ≥ 0.6, target 27 (widths 0.096, 0.124, 0.142; centres -0.68, +0.63, -0.05; amplitudes +0.48, +0.48, -0.56): dipole run-off, 2 kernels, residual 1.7e-01; the stalled block: 2 kernels, max|w| 3002 (the accepted polish's, solver-dependent), its closest opposite pair 2.30e-04 apart; Δ ≥ 0.6, target 29 (widths 0.197, 0.054, 0.088; centres +0.04, +0.69, -0.64; amplitudes -0.59, +0.85, +0.72): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 3956 (the accepted polish's, solver-dependent), its closest opposite pair 1.92e-04 apart. The width split (σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), the pair's third moment cancelled) against the observed.

## D. A narrow feature on a wide pedestal

The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):

| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |
|---|---|---|---|---|---|---|
| 0 | 3 | coinciding pair (+0.250/+0.250, Δc 2.94e-04, the width split -1.6515e-15 (predicted -4.0779e-18)) | 1.5e-04, 1.8e-07, 2.5e-01 | 2.6e-13 | 5 | 3.53e-02 |
| 0.1 | 2 | exact | 7.8e-18, 6.9e-18, 5.6e-16 | 1.8e-16 | 4 | 4.35e-02 |
| 0.2 | 2 | exact | 1.5e-17, 6.9e-18, 5.6e-16 | 2.1e-16 | 4 | 4.52e-02 |
| 0.3 | 2 | exact | 1.0e-15, 3.3e-16, 7.2e-15 | 1.0e-15 | 4 | 3.68e-02 |
| 0.45 | 2 | exact | 5.3e-18, 0.0e+00, 5.6e-16 | 2.2e-16 | 4 | 5.05e-02 |
| 0.6 | 3 | coinciding pair (+0.606/+0.194, Δc 1.03e-03, the width split +3.0498e-07 (predicted +3.0516e-07)) | 2.5e-04, 2.5e-07, 1.9e-01 | 1.9e-13 | 4 | 7.25e-02 |

## E. Under noise: the discrepancy stop

C's first 20 targets at Δ ≥ 0.6, one draw each, the path with moving widths to 5 entries (past the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against the data) is within (1 + 0.1)s; its kernels and its error against the truth (rms on the samples), median [range] and the rms over the draws, against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3, the rms of the errors), and per draw against the oracle fit (the K-kernel least-squares fit of the data started at the truth):

| s | stopped | kernels | error: median [range]; rms | the oracle | the stopped error over the oracle fit's: median [range] (within 10⁻⁴ of it: how many, the largest deviation among them) | not stopped: the last block's rms |
|---|---|---|---|---|---|---|
| 0.001 | 19 of 20 | 3 [3, 4] | 1.42e-04 [7.68e-05, 2.48e-04]; 1.648e-04 | 1.674e-04 | 1.0000 [1.0000, 1.0486] (18 of 19, 7e-08) | 2.16e-01 (budget) |
| 0.01 | 19 of 20 | 3 [3, 4] | 1.77e-03 [6.58e-04, 2.17e-03]; 1.750e-03 | 1.674e-03 | 1.0000 [1.0000, 1.2670] (18 of 19, 7e-08) | 2.16e-01 (budget) |

## F. The real loop at a small mass

The free law with moving centres and widths from α = 0, masses ε = 1e-06, 1e-08 (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to 1.3, 6, 30, 150 times the path's last entry, on C's first 4 targets at Δ ≥ 0.6: its formed kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest |w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair and residual:

| target | the path's kernels (closest pair; residual) | ε | times | formed by the loop | max \|Δθ\| | max \|Δw\| | the others' largest \|w\| | the loop's closest pair: distance (\|Δc\|, \|Δσ\|; coefficients) | residual |
|---|---|---|---|---|---|---|---|---|---|
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 1.3 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 6 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 30 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 150 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 1.3 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 6 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 30 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 150 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 1.3 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 6 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 30 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 150 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 8.9e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 1.3 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 6 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 30 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 150 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 1.3 | [3, 10, 11, 17] | 9.0e-04 | 9.8e-02 | 5.4e-06 | 7.0e-03 (1.5e-03, 5.5e-04; +0.329/+0.489) | 5.2e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 6 | [3, 10, 11, 17] | 1.3e-03 | 9.8e-02 | 5.4e-06 | 8.8e-03 (2.2e-03, 3.9e-04; +0.329/+0.489) | 5.0e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 30 | [3, 10, 11, 17] | 1.9e-03 | 9.8e-02 | 5.3e-06 | 1.2e-02 (3.2e-03, 2.8e-04; +0.329/+0.489) | 4.8e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 150 | [3, 10, 11, 17] | 2.7e-03 | 9.8e-02 | 5.0e-06 | 1.6e-02 (4.5e-03, 1.9e-04; +0.329/+0.489) | 4.5e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 1.3 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 1.2e-06; +0.329/+0.489) | 2.6e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 6 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 30 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 150 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 1.3 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 6 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 30 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 150 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 1.3 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 6 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 30 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 150 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |

## G. The blocks that never settle: a dipole

C's paths stopped by their budget. The stalled block (C's path again; its state when the budget ran out): its kernels, the largest |w|, the closest opposite-sign pair's distance (|Δc|/σ + |Δlog σ|), its residual against the optimal dipole fit's (below), its polishes (status/evaluations/accepted, the largest coefficient over the block's). The block from its entry, flowed on the critical curve to s = 1e+11: the largest |w| and the pair's distance at s = 1000, 100000, 1e+07, 1e+09, 3.2e+10, 1e+11; over the last factor in s the exponents of the largest |w|, the distance, the pair's |Δc| and |Δσ| and the residual's excess over the dipole infimum; the dipole law's constant K = v²(E − E_inf) at the last two outputs (v the pair's half-difference |w_a − w_b|/2, E_inf the optimal dipole fit's — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄) by least squares), the growth of v³ over the last interval against 3μ₀K (μ₀ = 4Gc/A: the coefficient's mobility μ₀|w| on the exponential branch), the moment's σ-share |p_σ|/|p|:

| Δ ≥ | target | entries | the stalled block: kernels, max\|w\|, distance, residual (infimum), polishes | the flow: largest \|w\| | distance | exponents: \|w\|, distance, \|Δc\|, \|Δσ\|, the residual's excess | K | Δv³/Δs, 3μ₀K | \|p_σ\|/\|p\| |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 2 | 2, 3910, 3.12e-04, 1.229967e-01 (1.229967e-01), 3/3358/a/846 | 0.895, 4.62, 20.2, 92, 291, 425 | 5.95e+00, 2.91e-01, 6.18e-02, 1.33e-02, 4.20e-03, 2.87e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 2.55977e-03, 2.55976e-03 | 7.67931e-04, 7.67927e-04 | 0.962 |
| 0.4 | 18 | 2 | 2, 3964, 2.30e-04, 9.689918e-02 (9.689918e-02), 3/2904/a/1.01e+03 | 0.572, 3.92, 17.8, 81.4, 258, 377 | 4.88e+00, 2.53e-01, 5.20e-02, 1.12e-02, 3.53e-03, 2.42e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 1.77965e-03, 1.77964e-03 | 5.33894e-04, 5.33893e-04 | 0.987 |
| 0.6 | 12 | 2 | 2, 1821, 6.13e-04, 5.423548e-02 (5.423547e-02), 3/4572/a/592 | 0.867, 3.07, 9.64, 40.4, 127, 185 | 5.55e+00, 4.29e-01, 1.21e-01, 2.79e-02, 8.82e-03, 6.03e-03 | +0.3322, -0.3331, -0.3331, -0.3331, -0.6663 | 2.10595e-04, 2.10560e-04 | 6.31839e-05, 6.31681e-05 | 0.956 |
| 0.6 | 20 | 3 | 3, 1431, 4.08e-04, 1.185147e-01 (1.185147e-01), 3/2968/a/1.04e+03 | 0.906, 1.37, 7.25, 32.9, 104, 152 | 1.80e+00, 5.19e-01, 8.31e-02, 1.78e-02, 5.62e-03, 3.84e-03 | +0.3328, -0.3333, -0.3333, -0.3333, -0.6667 | 1.16958e-04, 1.16958e-04 | 3.50872e-05, 3.50873e-05 | 0.734 |
| 0.6 | 27 | 2 | 2, 3002, 2.30e-04, 7.510333e-02 (7.510333e-02), 3/2981/a/1.09e+03 | 0.566, 2.75, 12.5, 57.1, 181, 264 | 5.16e+00, 2.78e-01, 5.63e-02, 1.21e-02, 3.82e-03, 2.61e-03 | +0.3329, -0.3333, -0.3333, -0.3333, -0.6667 | 6.10916e-04, 6.10914e-04 | 1.83275e-04, 1.83274e-04 | 0.994 |
| 0.6 | 29 | 2 | 2, 3956, 1.92e-04, 1.430698e-01 (1.430698e-01), 3/2876/a/1.12e+03 | 0.614, 3.52, 16.2, 74.2, 235, 343 | 4.63e+00, 2.38e-01, 4.78e-02, 1.03e-02, 3.24e-03, 2.22e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 1.34468e-03, 1.34468e-03 | 4.03404e-04, 4.03404e-04 | 0.987 |

The real loop at the masses ε from α = 0 (`loop_nd`, rtol 10⁻⁹) at τ = 400, 1600, 6400: the formed kernels, the dipole pair (the largest |w| and the closest kernel of the other sign: its distance |Δc|/σ + |Δlog σ|, coefficients, its centre of mass's offset from the nearest feature in that feature's width and its net coefficient W against the feature's amplitude, (W − A)/|A|; landed (L) if both within 0.01), the coinciding pairs (within 0.01), the largest |w| over the largest amplitude, the residual:

| Δ ≥ | target | ε | τ | formed | the dipole pair: distance, coefficients, offset (σ*), (W − A)/\|A\| | coinciding pairs | max\|w\|/max\|A\| | residual |
|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 1e-10 | 400 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-10 | 1600 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-10 | 6400 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-09 | 400 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-09 | 1600 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-09 | 6400 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-08 | 400 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 1600 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 6400 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-07 | 400 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-07 | 1600 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-07 | 6400 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-06 | 400 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 1600 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 6400 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-05 | 400 | 4 | 3.12e-01, +3.637/-1.563, 1.3e+00, +2.7e+00 | none | 4.28 | 2.1e-03 |
| 0.4 | 9 | 1e-05 | 1600 | 4 | 2.88e-01, +4.106/-1.818, 1.3e+00, +3.1e+00 | none | 4.83 | 1.9e-03 |
| 0.4 | 9 | 1e-05 | 6400 | 5 | 2.27e-01, +4.290/-1.976, 1.2e+00, +3.1e+00 | none | 5.05 | 6.2e-04 |
| 0.4 | 18 | 1e-10 | 400 | 4 | 6.42e-03, -11.635/+11.130, 7.9e-04, +9.4e-03 L | -11.635/+11.130 | 16.21 | 2.4e-07 |
| 0.4 | 18 | 1e-10 | 1600 | 4 | 4.81e-03, -11.636/+11.129, 3.5e-04, +5.4e-03 L | -11.636/+11.129 | 16.21 | 8.2e-08 |
| 0.4 | 18 | 1e-10 | 6400 | 4 | 3.70e-03, -11.637/+11.128, 1.7e-04, +3.2e-03 L | -11.637/+11.128 | 16.21 | 3.1e-08 |
| 0.4 | 18 | 1e-09 | 400 | 4 | 1.17e-02, -8.025/+7.522, 1.4e-03, +1.4e-02 | none | 11.18 | 5.6e-07 |
| 0.4 | 18 | 1e-09 | 1600 | 4 | 8.79e-03, -8.026/+7.521, 6.5e-04, +8.3e-03 L | -8.026/+7.521 | 11.18 | 1.9e-07 |
| 0.4 | 18 | 1e-09 | 6400 | 4 | 6.75e-03, -8.027/+7.520, 3.1e-04, +5.0e-03 L | -8.027/+7.520 | 11.18 | 7.0e-08 |
| 0.4 | 18 | 1e-08 | 400 | 4 | 2.13e-02, -5.588/+5.089, 2.6e-03, +2.2e-02 | none | 7.79 | 1.4e-06 |
| 0.4 | 18 | 1e-08 | 1600 | 4 | 1.61e-02, -5.591/+5.087, 1.2e-03, +1.3e-02 | none | 7.79 | 4.5e-07 |
| 0.4 | 18 | 1e-08 | 6400 | 4 | 1.24e-02, -5.592/+5.086, 5.9e-04, +7.7e-03 L | none | 7.79 | 1.7e-07 |
| 0.4 | 18 | 1e-07 | 400 | 4 | 1.91e-02, +8.468/-7.882, 3.8e-03, -2.8e-02 | none | 11.80 | 8.1e-06 |
| 0.4 | 18 | 1e-07 | 1600 | 4 | 8.87e-03, +8.474/-7.876, 4.8e-04, -6.7e-03 L | -7.876/+8.474 | 11.81 | 4.7e-07 |
| 0.4 | 18 | 1e-07 | 6400 | 4 | 6.44e-03, +8.475/-7.875, 1.9e-04, -3.6e-03 L | -7.875/+8.475 | 11.81 | 1.4e-07 |
| 0.4 | 18 | 1e-06 | 400 | 4 | 1.29e-01, +5.218/-4.429, 4.7e-01, +3.1e-01 | none | 7.27 | 1.0e-03 |
| 0.4 | 18 | 1e-06 | 1600 | 4 | 2.29e-02, +5.115/-4.522, 1.5e-03, -1.5e-02 | none | 7.13 | 2.3e-06 |
| 0.4 | 18 | 1e-06 | 6400 | 4 | 1.47e-02, +5.118/-4.519, 4.5e-04, -6.4e-03 L | none | 7.13 | 4.3e-07 |
| 0.4 | 18 | 1e-05 | 400 | 4 | 2.23e-01, +3.233/-2.424, 4.8e-01, +3.4e-01 | none | 4.50 | 1.0e-03 |
| 0.4 | 18 | 1e-05 | 1600 | 5 | 1.91e-01, +3.128/-2.515, 2.7e-01, +1.8e-02 | none | 4.36 | 4.2e-04 |
| 0.4 | 18 | 1e-05 | 6400 | 5 | 8.67e-02, +3.091/-2.548, 2.8e-03, -9.9e-02 | none | 4.31 | 6.0e-06 |
| 0.6 | 12 | 1e-10 | 400 | 4 | 6.35e-03, +11.468/-10.642, 1.9e-04, -3.4e-03 L | -10.642/+11.468 | 13.84 | 1.2e-07 |
| 0.6 | 12 | 1e-10 | 1600 | 4 | 4.74e-03, +11.468/-10.642, 8.0e-05, -1.9e-03 L | -10.642/+11.468 | 13.85 | 3.9e-08 |
| 0.6 | 12 | 1e-10 | 6400 | 4 | 3.66e-03, +11.469/-10.641, 3.7e-05, -1.2e-03 L | -10.641/+11.469 | 13.85 | 1.4e-08 |
| 0.6 | 12 | 1e-09 | 400 | 4 | 1.15e-02, +8.265/-7.441, 4.0e-04, -5.6e-03 L | none | 9.98 | 3.3e-07 |
| 0.6 | 12 | 1e-09 | 1600 | 4 | 8.42e-03, +8.266/-7.440, 1.6e-04, -3.0e-03 L | -7.440/+8.266 | 9.98 | 9.6e-08 |
| 0.6 | 12 | 1e-09 | 6400 | 4 | 6.42e-03, +8.266/-7.439, 7.1e-05, -1.8e-03 L | -7.439/+8.266 | 9.98 | 3.3e-08 |
| 0.6 | 12 | 1e-08 | 400 | 4 | 2.14e-02, +6.107/-5.286, 9.5e-04, -9.9e-03 L | none | 7.37 | 1.0e-06 |
| 0.6 | 12 | 1e-08 | 1600 | 4 | 1.51e-02, +6.109/-5.285, 3.4e-04, -5.0e-03 L | none | 7.37 | 2.6e-07 |
| 0.6 | 12 | 1e-08 | 6400 | 4 | 1.12e-02, +6.110/-5.284, 1.4e-04, -2.8e-03 L | none | 7.38 | 8.5e-08 |
| 0.6 | 12 | 1e-07 | 400 | 4 | 4.04e-02, +4.648/-3.835, 2.5e-03, -1.9e-02 | none | 5.61 | 3.7e-06 |
| 0.6 | 12 | 1e-07 | 1600 | 4 | 2.73e-02, +4.652/-3.831, 8.1e-04, -8.9e-03 L | none | 5.62 | 8.3e-07 |
| 0.6 | 12 | 1e-07 | 6400 | 4 | 1.97e-02, +4.654/-3.830, 3.1e-04, -4.7e-03 L | none | 5.62 | 2.4e-07 |
| 0.6 | 12 | 1e-06 | 400 | 4 | 1.55e-01, -5.185/+4.604, 1.2e+00, -1.7e+00 | none | 6.26 | 4.1e-04 |
| 0.6 | 12 | 1e-06 | 1600 | 5 | 1.31e-01, -5.164/+4.637, 1.3e+00, -1.6e+00 | none | 6.23 | 2.2e-04 |
| 0.6 | 12 | 1e-06 | 6400 | 5 | 8.21e-02, -5.103/+4.687, 1.1e+00, +1.5e-01 | none | 6.16 | 1.3e-04 |
| 0.6 | 12 | 1e-05 | 400 | 5 | 4.93e-01, -3.409/+0.988, 1.5e+00, -3.9e+00 | none | 4.12 | 1.0e-03 |
| 0.6 | 12 | 1e-05 | 1600 | 5 | 4.44e-01, -3.347/+1.018, 1.1e+00, -3.8e+00 | none | 4.04 | 2.9e-04 |
| 0.6 | 12 | 1e-05 | 6400 | 5 | 1.53e-01, -3.285/+1.038, 6.1e-01, -3.7e+00 | none | 3.97 | 7.2e-06 |
| 0.6 | 20 | 1e-10 | 400 | 4 | 8.80e-03, -4.838/+4.332, 1.6e-04, +3.0e-03 L | +4.332/-4.838 | 5.84 | 1.2e-07 |
| 0.6 | 20 | 1e-10 | 1600 | 4 | 6.84e-03, -4.838/+4.332, 7.5e-05, +1.8e-03 L | +4.332/-4.838 | 5.84 | 4.3e-08 |
| 0.6 | 20 | 1e-10 | 6400 | 4 | 5.38e-03, -4.838/+4.332, 3.7e-05, +1.1e-03 L | +4.332/-4.838 | 5.84 | 1.7e-08 |
| 0.6 | 20 | 1e-09 | 400 | 4 | 1.61e-02, -3.369/+2.864, 3.0e-04, +4.6e-03 L | none | 4.07 | 2.8e-07 |
| 0.6 | 20 | 1e-09 | 1600 | 4 | 1.24e-02, -3.369/+2.864, 1.4e-04, +2.8e-03 L | none | 4.07 | 1.0e-07 |
| 0.6 | 20 | 1e-09 | 6400 | 4 | 9.70e-03, -3.369/+2.863, 6.7e-05, +1.7e-03 L | +2.863/-3.369 | 4.07 | 3.8e-08 |
| 0.6 | 20 | 1e-08 | 400 | 4 | 3.02e-02, -2.369/+1.865, 6.1e-04, +7.3e-03 L | none | 2.86 | 7.5e-07 |
| 0.6 | 20 | 1e-08 | 1600 | 4 | 2.29e-02, -2.369/+1.865, 2.7e-04, +4.3e-03 L | none | 2.86 | 2.5e-07 |
| 0.6 | 20 | 1e-08 | 6400 | 4 | 1.77e-02, -2.370/+1.864, 1.3e-04, +2.6e-03 L | none | 2.86 | 8.9e-08 |
| 0.6 | 20 | 1e-07 | 400 | 4 | 5.76e-02, -1.694/+1.193, 1.3e-03, +1.2e-02 | none | 2.05 | 2.2e-06 |
| 0.6 | 20 | 1e-07 | 1600 | 4 | 4.30e-02, -1.695/+1.192, 5.5e-04, +6.8e-03 L | none | 2.05 | 6.8e-07 |
| 0.6 | 20 | 1e-07 | 6400 | 4 | 3.27e-02, -1.696/+1.191, 2.4e-04, +4.0e-03 L | none | 2.05 | 2.3e-07 |
| 0.6 | 20 | 1e-06 | 400 | 4 | 3.81e+00, +0.828/-0.227, 1.8e+00, -2.7e-01 | none | 1.00 | 1.2e-05 |
| 0.6 | 20 | 1e-06 | 1600 | 4 | 3.88e+00, +0.828/-0.225, 1.8e+00, -2.7e-01 | none | 1.00 | 3.2e-06 |
| 0.6 | 20 | 1e-06 | 6400 | 4 | 3.93e+00, +0.828/-0.224, 1.9e+00, -2.7e-01 | none | 1.00 | 8.7e-07 |
| 0.6 | 20 | 1e-05 | 400 | 4 | 3.80e+00, +0.826/-0.285, 2.5e+00, -3.5e-01 | none | 1.00 | 4.0e-05 |
| 0.6 | 20 | 1e-05 | 1600 | 4 | 3.86e+00, +0.825/-0.280, 2.5e+00, -3.4e-01 | none | 1.00 | 1.1e-05 |
| 0.6 | 20 | 1e-05 | 6400 | 4 | 3.91e+00, +0.825/-0.278, 2.5e+00, -3.4e-01 | none | 1.00 | 3.2e-06 |
| 0.6 | 27 | 1e-10 | 400 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-10 | 1600 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-10 | 6400 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-09 | 400 | 4 | 3.88e-03, +8.812/-8.330, 1.0e-04, -2.3e-03 L | -8.330/+8.812 | 15.74 | 7.7e-08 |
| 0.6 | 27 | 1e-09 | 1600 | 4 | 2.93e-03, +8.813/-8.330, 4.5e-05, -1.3e-03 L | -8.330/+8.813 | 15.74 | 2.6e-08 |
| 0.6 | 27 | 1e-09 | 6400 | 4 | 2.27e-03, +8.813/-8.330, 2.1e-05, -7.9e-04 L | -8.330/+8.813 | 15.74 | 9.3e-09 |
| 0.6 | 27 | 1e-08 | 400 | 4 | 7.52e-03, +6.043/-5.562, 2.3e-04, -3.9e-03 L | -5.562/+6.043 | 10.79 | 2.2e-07 |
| 0.6 | 27 | 1e-08 | 1600 | 4 | 5.54e-03, +6.044/-5.561, 9.2e-05, -2.1e-03 L | -5.561/+6.044 | 10.79 | 6.8e-08 |
| 0.6 | 27 | 1e-08 | 6400 | 4 | 4.21e-03, +6.044/-5.561, 4.1e-05, -1.2e-03 L | -5.561/+6.044 | 10.80 | 2.3e-08 |
| 0.6 | 27 | 1e-07 | 400 | 4 | 1.53e-02, +4.152/-3.672, 5.7e-04, -7.1e-03 L | none | 7.42 | 7.8e-07 |
| 0.6 | 27 | 1e-07 | 1600 | 4 | 1.09e-02, +4.153/-3.671, 2.1e-04, -3.7e-03 L | none | 7.42 | 2.1e-07 |
| 0.6 | 27 | 1e-07 | 6400 | 4 | 8.07e-03, +4.153/-3.671, 8.7e-05, -2.0e-03 L | -3.671/+4.153 | 7.42 | 6.5e-08 |
| 0.6 | 27 | 1e-06 | 400 | 5 | 2.98e-01, +3.872/-2.075, 2.5e-01, +4.2e+00 | none | 6.92 | 3.8e-05 |
| 0.6 | 27 | 1e-06 | 1600 | 5 | 2.23e-01, +3.898/-2.061, 2.0e-01, +4.3e+00 | none | 6.96 | 8.9e-06 |
| 0.6 | 27 | 1e-06 | 6400 | 5 | 1.79e-01, +3.911/-2.055, 1.6e-01, +4.3e+00 | none | 6.99 | 2.7e-06 |
| 0.6 | 27 | 1e-05 | 400 | 4 | 1.01e-01, +1.745/-1.282, 7.7e-03, -4.3e-02 | none | 3.12 | 1.3e-05 |
| 0.6 | 27 | 1e-05 | 1600 | 4 | 8.25e-02, +1.748/-1.279, 3.8e-03, -3.0e-02 | none | 3.12 | 4.1e-06 |
| 0.6 | 27 | 1e-05 | 6400 | 4 | 6.91e-02, +1.751/-1.278, 1.9e-03, -2.2e-02 | none | 3.13 | 1.5e-06 |
| 0.6 | 29 | 1e-10 | 400 | 4 | 5.07e-03, -7.381/+6.789, 7.0e-05, +1.8e-03 L | +6.789/-7.381 | 8.63 | 6.8e-08 |
| 0.6 | 29 | 1e-10 | 1600 | 4 | 3.90e-03, -7.381/+6.788, 3.2e-05, +1.0e-03 L | +6.788/-7.381 | 8.63 | 2.4e-08 |
| 0.6 | 29 | 1e-10 | 6400 | 4 | 3.05e-03, -7.381/+6.788, 1.6e-05, +6.5e-04 L | +6.788/-7.381 | 8.63 | 9.2e-09 |
| 0.6 | 29 | 1e-09 | 400 | 4 | 9.53e-03, -5.034/+4.442, 1.4e-04, +2.7e-03 L | +4.442/-5.034 | 5.89 | 1.7e-07 |
| 0.6 | 29 | 1e-09 | 1600 | 4 | 7.23e-03, -5.034/+4.442, 6.1e-05, +1.6e-03 L | +4.442/-5.034 | 5.89 | 5.7e-08 |
| 0.6 | 29 | 1e-09 | 6400 | 4 | 5.61e-03, -5.035/+4.442, 2.9e-05, +9.7e-04 L | +4.442/-5.035 | 5.89 | 2.1e-08 |
| 0.6 | 29 | 1e-08 | 400 | 4 | 1.89e-02, -3.389/+2.798, 2.9e-04, +4.6e-03 L | none | 3.96 | 4.7e-07 |
| 0.6 | 29 | 1e-08 | 1600 | 4 | 1.40e-02, -3.390/+2.798, 1.2e-04, +2.5e-03 L | none | 3.96 | 1.5e-07 |
| 0.6 | 29 | 1e-08 | 6400 | 4 | 1.07e-02, -3.390/+2.797, 5.5e-05, +1.5e-03 L | none | 3.96 | 5.0e-08 |
| 0.6 | 29 | 1e-07 | 400 | 4 | 4.07e-02, -2.248/+1.660, 7.0e-04, +8.1e-03 L | none | 2.63 | 1.6e-06 |
| 0.6 | 29 | 1e-07 | 1600 | 4 | 2.92e-02, -2.250/+1.659, 2.7e-04, +4.3e-03 L | none | 2.63 | 4.3e-07 |
| 0.6 | 29 | 1e-07 | 6400 | 4 | 2.17e-02, -2.250/+1.658, 1.1e-04, +2.4e-03 L | none | 2.63 | 1.4e-07 |
| 0.6 | 29 | 1e-06 | 400 | 4 | 8.57e-02, -1.578/+0.993, 1.8e-03, +1.5e-02 | none | 1.85 | 5.8e-06 |
| 0.6 | 29 | 1e-06 | 1600 | 4 | 5.98e-02, -1.580/+0.991, 6.2e-04, +7.4e-03 L | none | 1.85 | 1.5e-06 |
| 0.6 | 29 | 1e-06 | 6400 | 4 | 4.33e-02, -1.581/+0.990, 2.4e-04, +3.9e-03 L | none | 1.85 | 4.1e-07 |
| 0.6 | 29 | 1e-05 | 400 | 4 | 6.32e+00, +0.854/-0.399, 6.1e+00, +1.8e+00 | none | 1.00 | 2.3e-05 |
| 0.6 | 29 | 1e-05 | 1600 | 4 | 6.36e+00, +0.854/-0.395, 6.1e+00, +1.8e+00 | none | 1.00 | 6.0e-06 |
| 0.6 | 29 | 1e-05 | 6400 | 4 | 6.39e+00, +0.854/-0.394, 6.1e+00, +1.8e+00 | none | 1.00 | 1.7e-06 |

(ii)(a)'s width split on the landed pairs (L, 72 rows, 59 distinct to six decimals: frozen states repeat): the observed σ_b − σ_a over (w_b − w_a)Δc²/(6Wσ) 0.9828–1.0458 (median 0.9993), at |w_a − w_b||Δc|/(|W|σ) = 8.74e-03–2.84e-01; its largest deviation 4.6e-02 at the parameter 2.76e-01.

The three phases (`phases`): the dipole block's entry τ_d (the path's) and the next kernel's formation τ₃ in the loop; over the grid times between them the slopes of v³ (v the dipole pair's half-difference) against 3μ₀KΛ(ε) — the block's law in its fast time Λ(ε)τ, K and μ₀ from the flow's last output — (median [range]); the landed v (at τ = 400) against the law's (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}; landed at τ = 400, 1600, 6400; over 1600 → 6400 the exponents in τ of the pair's distance and of the residual:

| Δ ≥ | target | ε | τ_d | τ₃ | d(v³)/dτ / 3μ₀KΛ | v at 400 | the law's v | landed | exponents: distance, residual |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 1e-05 | 17.75 | 22 | 1.360 [1.360, 1.360] | 2.6 | 1.936 | no, no, no | -0.173, -0.794 |
| 0.4 | 9 | 1e-06 | 17.75 | 28 | 1.108 [1.077, 1.223] | 4.244 | 3.81 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-07 | 17.75 | 30 | 1.046 [1.032, 1.123] | 6.414 | 5.935 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-08 | 17.75 | 30 | 1.022 [1.015, 1.063] | 9.467 | 8.711 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-09 | 17.75 | 30 | 1.011 [1.007, 1.031] | 13.92 | 12.79 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-10 | 17.75 | 30 | 1.005 [1.003, 1.015] | 20.43 | 18.77 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 18 | 1e-05 | 25 | 26 | — | 2.829 | 1.057 | no, no, no | -0.571, -3.061 |
| 0.4 | 18 | 1e-06 | 25 | 26 | — | 4.823 | 1.552 | no, no, yes | -0.322, -1.211 |
| 0.4 | 18 | 1e-07 | 25 | 28 | — | 8.175 | 3.288 | no, yes, yes | -0.230, -0.880 |
| 0.4 | 18 | 1e-08 | 25 | 34 | 1.013 [1.010, 1.021] | 5.339 | 6.962 | no, no, yes | -0.190, -0.723 |
| 0.4 | 18 | 1e-09 | 25 | 36 | 1.005 [1.004, 1.010] | 7.773 | 10.93 | no, yes, yes | -0.191, -0.718 |
| 0.4 | 18 | 1e-10 | 25 | 36 | 1.003 [1.002, 1.005] | 11.38 | 16.04 | yes, yes, yes | -0.189, -0.713 |
| 0.6 | 12 | 1e-05 | 24.35 | 26 | — | 2.199 | 0.6147 | no, no, no | -0.771, -2.659 |
| 0.6 | 12 | 1e-06 | 24.35 | 26 | — | 4.894 | 0.9022 | no, no, no | -0.336, -0.374 |
| 0.6 | 12 | 1e-07 | 24.35 | 57.4349 | 1.768 [1.569, 2.688] | 4.241 | 3.594 | no, yes, yes | -0.235, -0.901 |
| 0.6 | 12 | 1e-08 | 24.35 | 57.4349 | 1.427 [1.310, 1.988] | 5.697 | 5.276 | yes, yes, yes | -0.212, -0.821 |
| 0.6 | 12 | 1e-09 | 24.35 | 57.4349 | 1.228 [1.162, 1.558] | 7.853 | 7.744 | yes, yes, yes | -0.196, -0.762 |
| 0.6 | 12 | 1e-10 | 24.35 | 65.9754 | 1.112 [1.073, 1.304] | 11.06 | 12.27 | yes, yes, yes | -0.185, -0.726 |
| 0.6 | 20 | 1e-05 | 32.22 | 34 | — | 0.5552 | 0.5173 | no, no, no | +0.009, -0.886 |
| 0.6 | 20 | 1e-06 | 32.22 | 34 | — | 0.5273 | 0.7593 | no, no, no | +0.010, -0.937 |
| 0.6 | 20 | 1e-07 | 32.22 | 36 | — | 1.443 | 1.433 | no, yes, yes | -0.198, -0.791 |
| 0.6 | 20 | 1e-08 | 32.22 | 38 | 0.984 [0.984, 0.984] | 2.117 | 2.424 | yes, yes, yes | -0.186, -0.741 |
| 0.6 | 20 | 1e-09 | 32.22 | 38 | 0.993 [0.993, 0.993] | 3.116 | 3.557 | yes, yes, yes | -0.178, -0.708 |
| 0.6 | 20 | 1e-10 | 32.22 | 38 | 0.997 [0.997, 0.997] | 4.585 | 5.222 | yes, yes, yes | -0.173, -0.689 |
| 0.6 | 27 | 1e-05 | 31.54 | 32 | — | 1.513 | 0.5708 | no, no, no | -0.128, -0.740 |
| 0.6 | 27 | 1e-06 | 31.54 | 32 | — | 2.974 | 0.8378 | no, no, no | -0.159, -0.846 |
| 0.6 | 27 | 1e-07 | 31.54 | 36 | 1.064 [1.064, 1.064] | 3.912 | 2.628 | yes, yes, yes | -0.219, -0.843 |
| 0.6 | 27 | 1e-08 | 31.54 | 38 | 1.026 [1.018, 1.033] | 5.802 | 4.364 | yes, yes, yes | -0.198, -0.777 |
| 0.6 | 27 | 1e-09 | 31.54 | 46 | 1.006 [1.004, 1.016] | 8.571 | 8.38 | yes, yes, yes | -0.185, -0.730 |
| 0.6 | 27 | 1e-10 | 31.54 | 48 | 1.002 [1.002, 1.008] | 12.86 | 12.84 | yes, yes, yes | -0.000, -0.000 |
| 0.6 | 29 | 1e-05 | 26.75 | 28 | — | 0.6267 | 1.038 | no, no, no | +0.003, -0.925 |
| 0.6 | 29 | 1e-06 | 26.75 | 30 | — | 1.285 | 2.096 | no, yes, yes | -0.233, -0.906 |
| 0.6 | 29 | 1e-07 | 26.75 | 30 | — | 1.954 | 3.076 | yes, yes, yes | -0.213, -0.833 |
| 0.6 | 29 | 1e-08 | 26.75 | 30 | — | 3.094 | 4.515 | yes, yes, yes | -0.194, -0.763 |
| 0.6 | 29 | 1e-09 | 26.75 | 30 | — | 4.738 | 6.628 | yes, yes, yes | -0.183, -0.720 |
| 0.6 | 29 | 1e-10 | 26.75 | 30 | — | 7.085 | 9.728 | yes, yes, yes | -0.177, -0.698 |

The landed v against the mass: v at τ = 400 and its ratio per decade of ε (10^{1/6} = 1.4678 if v ∝ ε^{−1/6}):

| Δ ≥ | target | v at ε = 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 1e-10 | ratios per decade |
|---|---|---|---|
| 0.4 | 9 | 2.6, 4.244, 6.414, 9.467, 13.92, 20.43 | 1.6324, 1.5112, 1.4760, 1.4699, 1.4682 |
| 0.4 | 18 | 2.829, 4.823, 8.175, 5.339, 7.773, 11.38 | 1.7052, 1.6948, 0.6531, 1.4560, 1.4643 |
| 0.6 | 12 | 2.199, 4.894, 4.241, 5.697, 7.853, 11.06 | 2.2262, 0.8666, 1.3431, 1.3785, 1.4078 |
| 0.6 | 20 | 0.5552, 0.5273, 1.443, 2.117, 3.116, 4.585 | 0.9496, 2.7370, 1.4669, 1.4721, 1.4713 |
| 0.6 | 27 | 1.513, 2.974, 3.912, 5.802, 8.571, 12.86 | 1.9648, 1.3156, 1.4832, 1.4772, 1.5002 |
| 0.6 | 29 | 0.6267, 1.285, 1.954, 3.094, 4.738, 7.085 | 2.0507, 1.5203, 1.5832, 1.5316, 1.4953 |

Over the six targets: landed at τ = 6400 in 0 of 6 at ε = 1e-05, 3 of 6 at ε = 1e-06, 6 of 6 at ε = 1e-07, 6 of 6 at ε = 1e-08, 6 of 6 at ε = 1e-09, 6 of 6 at ε = 1e-10; at ε = 1e-10 the dipole phase's slope over the law 0.997–1.112 (the cells' medians, 5 cells with a phase), the landed v over the law's at the grid's τ₃ 0.71–1.09 (on the fine grid κ, below), v's ratio over the last decade 1.4078–1.5002; on the 25 cells landed at both 1600 and 6400, 6 not yet collapsing by τ = 6400 (|exponent| < 10⁻³: landed at a small p, the same law later — GL) (9 at 1e-06, 9 at 1e-07, 9 at 1e-08, 9 at 1e-09, 9 at 1e-10, 27 at 1e-10), the others collapsing with the exponents -0.235 … -0.173 (distance) and -0.906 … -0.689 (residual).


The phases on a fine grid (τ = 10 … 140 by 0.1, with G's times; referee rounds 68, 74): τ₃ the first formation past the block's entry τ_d, τ₄ the next; over (τ_d, τ₃) the slopes of v³ over 3μ₀KΛ(ε) (median); v at the last time before τ₃ over the law's (3μ₀KΛ(ε)(t − τ_d))^{1/3} at that time t; at (τ₃ + τ₄)/2 the formed kernels, the dipole pair's offset from the nearest feature (σ*) and (W − A)/|A|, and over [τ₃, τ₄) the slopes of v³ over 3μ₀K′Λ(ε), K′ the three-kernel block's from the path (below; median); the first landing past τ₃ (centre of mass and W within 0.01): its τ, κ = v_land/v(τ₃⁻) with v(τ₃⁻) the loop's v at the last time before τ₃, κ′ = v_land over the law's v at the path's τ₃ (the loop's v before τ₃ carried by v³ ∝ τ − τ_d: `kappa_loop`, referee rounds 80, 83 — at larger masses not the loop's own κ, whose crossing is the loop's), and p = |w_a − w_b||Δθ|/(|W|σ̄); † another sequence — the loop's formed kernels at its τ₃ or τ₄ not the path's, by index (SQ, referee round 83) —, its κ not the path's:

| target | ε | τ_d | τ₃ | τ₄ | dipole: slope | v(τ₃⁻)/law | intermediate: formed, offset, (W − A)/\|A\|, slope | landing: τ, κ, κ′, p |
|---|---|---|---|---|---|---|---|---|
| 9 | 1e-05 | 17.75 | 22 | 26.2 | 1.2870 | 1.1230 | 3, 0.0358, +0.197, 25.7460 | — |
| 9 | 1e-06 | 17.75 | 27.7 | 81.5 | 1.1006 | 1.0455 | 3, 10.1, -0.723, 1.1058 | 151.572, 1.080, 1.028, 0.017 |
| 9 | 1e-07 | 17.75 | 28.5 | 110.4 | 1.0474 | 1.0231 | 3, 9.56, -0.711, 1.0530 | 128.8, 1.107, 1.082, 0.177 |
| 9 | 1e-08 | 17.75 | 28.9 | 118.8 | 1.0214 | 1.0115 | 3, 9.36, -0.706, 1.0255 | 124.6, 1.112, 1.100, 0.047 |
| 9 | 1e-09 | 17.75 | 29.1 | 122.1 | 1.0102 | 1.0056 | 3, 9.28, -0.704, 1.0120 | 124, 1.113, 1.108, 0.020 |
| 9 | 1e-10 | 17.75 | 29.1 | 123.5 | 1.0048 | 1.0027 | 3, 9.24, -0.703, 1.0056 | 124, 1.117, 1.112, 0.029 |
| 9 | 1e-11 | 17.75 | 29.2 | 124 | 1.0022 | 1.0013 | 3, 9.22, -0.702, 1.0026 | 124.2, 1.115, 1.113, 0.015 |
| 9 | 1e-12 | 17.75 | 29.2 | 124.2 | 1.0010 | 1.0006 | 3, 9.21, -0.702, 1.0012 | 124.3, 1.116, 1.114, 0.015 |
| 18 | 1e-05 | 25 | 25.1 | 35.3 | — | — | 3, 0.483, -0.013, -16.6473 | — |
| 18 | 1e-06 | 25 | 25.1 | 41.3 | — | — | 3, 0.51, -1.655, -21.2937 | 2785.76, —, —†, 0.286 |
| 18 | 1e-07 | 25 | 26.8 | 47.9 | 1.0850 | 1.0347 | 3, 0.997, -1.530, -25.1819 | 918.959, 2.905, 1.652†, 0.275 |
| 18 | 1e-08 | 25 | 34 | 95.2 | 1.0129 | 1.0067 | 3, 2.85, -0.457, 0.9945 | 3200, 0.765, 0.755, 0.284 |
| 18 | 1e-09 | 25 | 34.1 | 97 | 1.0061 | 1.0032 | 3, 2.85, -0.457, 0.9974 | 1055.61, 0.758, 0.752, 0.280 |
| 18 | 1e-10 | 25 | 34.2 | 97.6 | 1.0026 | 1.0015 | 3, 2.85, -0.457, 0.9988 | 400, 0.755, 0.751, 0.276 |
| 18 | 1e-11 | 25 | 34.2 | 97.8 | 1.0012 | 1.0007 | 3, 2.84, -0.457, 0.9995 | 174.11, 0.755, 0.752, 0.284 |
| 18 | 1e-12 | 25 | 34.3 | 97.9 | 1.0006 | 1.0003 | 3, 2.84, -0.457, 0.9997 | 120.7, 0.752, 0.752, 0.285 |
| 12 | 1e-05 | 24.35 | 24.4 | 69 | — | — | 3, 0.88, -0.548, -9.1534 | — |
| 12 | 1e-06 | 24.35 | 24.4 | 91 | — | — | 3, 14.8, +1.068, -10.7375 | — |
| 12 | 1e-07 | 24.35 | 54.8 | 80.2 | 1.7301 | 1.2579 | 3, 20.7, +0.864, 0.9920 | 1392.88, 0.966, 0.929, 0.279 |
| 12 | 1e-08 | 24.35 | 56.5 | 81.3 | 1.3910 | 1.1561 | 3, 20.3, +0.861, 0.9956 | 400, 0.944, 0.925, 0.284 |
| 12 | 1e-09 | 24.35 | 57.4349 | 81.5 | 1.2039 | 1.0898 | 3, 20.1, +0.860, 0.9977 | 174.11, 0.931, 0.921, 0.277 |
| 12 | 1e-10 | 24.35 | 58 | 81.6 | 1.1039 | 1.0492 | 3, 20, +0.859, 0.9988 | 102.9, 0.923, 0.918, 0.285 |
| 12 | 1e-11 | 24.35 | 58.2 | 81.7 | 1.0511 | 1.0258 | 3, 20, +0.859, 0.9994 | 87.3, 0.918, 0.915, 0.285 |
| 12 | 1e-12 | 24.35 | 58.4 | 81.7 | 1.0240 | 1.0130 | 3, 20, +0.858, 0.9997 | 83.2, 0.915, 0.914, 0.284 |
| 20 | 1e-05 | 32.22 | 32.3 | — | — | — | — | — |
| 20 | 1e-06 | 32.22 | 32.3 | — | — | — | — | — |
| 20 | 1e-07 | 32.22 | 36 | — | 0.9468 | 1.0010 | — | 696.44, 1.015, 0.992, 0.281 |
| 20 | 1e-08 | 32.22 | 36.1 | — | 0.9768 | 0.9964 | — | 200, 1.010, 0.995, 0.285 |
| 20 | 1e-09 | 32.22 | 36.2 | — | 0.9897 | 0.9971 | — | 78.5, 1.003, 0.998, 0.286 |
| 20 | 1e-10 | 32.22 | 36.2 | — | 0.9954 | 0.9983 | — | 47.7, 1.005, 0.999, 0.285 |
| 20 | 1e-11 | 32.22 | 36.2 | — | 0.9979 | 0.9991 | — | 39.4, 1.005, 0.999, 0.285 |
| 20 | 1e-12 | 32.22 | 36.2 | — | 0.9990 | 0.9995 | — | 37.1, 1.006, 1.000, 0.285 |
| 27 | 1e-05 | 31.54 | 31.6 | 40.5 | — | — | 3, 1.61, +0.837, 2.9797 | — |
| 27 | 1e-06 | 31.54 | 31.6 | 43.8 | — | — | 3, 0.238, +0.199, 4.6285 | — |
| 27 | 1e-07 | 31.54 | 34.1 | 45.5 | 1.0834 | 1.0327 | 3, 0.202, +0.134, 4.4686 | 229.74, 1.758, 0.949†, 0.280 |
| 27 | 1e-08 | 31.54 | 36.8 | 46.4 | 1.0214 | 1.0112 | 3, 0.2, +0.134, 4.4266 | 86, 1.417, 0.979†, 0.285 |
| 27 | 1e-09 | 31.54 | 44.8 | 46.9 | 1.0055 | 1.0031 | 3, 0.199, +0.134, 4.4036 | 56.3, 1.052, 0.993†, 0.285 |
| 27 | 1e-10 | 31.54 | 47.1 | 51.2 | 1.0024 | 1.0014 | 3, 0.728, +0.163, 1.0000 | 51.2, 1.021, 1.017, 0.009 |
| 27 | 1e-11 | 31.54 | 47.2 | 51.3 | 1.0011 | 1.0006 | 3, 0.727, +0.163, 1.0000 | 51.3, 1.020, 1.018, 0.009 |
| 27 | 1e-12 | 31.54 | 47.2 | 51.3 | 1.0005 | 1.0003 | 3, 0.727, +0.163, 1.0000 | 51.3, 1.021, 1.019, 0.010 |
| 29 | 1e-05 | 26.75 | 27.3 | 34.2 | — | 0.9906 | 3, 7.47, +1.350, 0.7795 | — |
| 29 | 1e-06 | 26.75 | 28.4 | 43.5 | 0.9532 | 0.9941 | 3, 7.47, +1.340, 0.9648 | 918.959, 0.790, 0.643, 0.287 |
| 29 | 1e-07 | 26.75 | 29 | 45.9 | 0.9973 | 0.9994 | 3, 7.49, +1.334, 0.9835 | 303.143, 0.729, 0.662, 0.279 |
| 29 | 1e-08 | 26.75 | 29.3 | 46.7 | 1.0015 | 1.0003 | 3, 7.5, +1.331, 0.9928 | 101.8, 0.753, 0.713, 0.286 |
| 29 | 1e-09 | 26.75 | 29.5 | 46.9 | 1.0014 | 1.0003 | 3, 7.5, +1.330, 0.9969 | 60.5, 0.765, 0.744, 0.285 |
| 29 | 1e-10 | 26.75 | 29.6 | 47 | 1.0008 | 1.0002 | 3, 7.5, +1.330, 0.9986 | 50.5, 0.770, 0.758, 0.284 |
| 29 | 1e-11 | 26.75 | 29.7 | 47 | 1.0004 | 1.0001 | 3, 7.5, +1.330, 0.9994 | 47.9, 0.767, 0.764, 0.284 |
| 29 | 1e-12 | 26.75 | 29.7 | 47 | 1.0002 | 1.0001 | 3, 7.51, +1.329, 0.9997 | 47.3, 0.769, 0.767, 0.270 |

The intermediate phase from the path (`kappa_from_path`, referee rounds 74, 80): K from the third order, −⟨r_dip, ∇³φ:p³⟩/24 at the block's optimal dipole fit, beside G's v²(E − E_inf); τ₃ from the duals; K′/K, K′ the three-kernel block's (the τ₃ entrant added), with that block's closure rms (below 1e-12 the block fits exactly: K′ = 0, no τ₄, κ = 1); τ₄ from the duals moved on by that block's residual; for K′ < 0 the dissolution time τ₃ + K(τ₃ − τ_d)/|K′|, at which v³ would reach 0 (the four stages need τ₄ before it); κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3}; the landing's v by the law v³ = 3μ₀Λ(ε)(K(τ₃ − τ_d) + K′(τ₄ − τ₃)) with the path's τ₃ and τ₄ — each against the loop's at the smallest mass (its κ as in the table above) —; and v at τ = 6400 over v_land at ε = 10⁻¹⁰, minus 1 (v after τ₄):

| target | K: third order (G's) | τ₃: path (loop) | K′/K (the closure's rms) | τ₄: path (loop) | dissolution | κ: path (loop: κ′; its own crossing's bracket) | v_land: law (loop) | v(6400)/v_land − 1 |
|---|---|---|---|---|---|---|---|---|
| 9 | 2.55975e-03 (2.55976e-03) | 29.16 (29.2) | +0.0461 (1.5e-02) | 124.3 (124.2) | — | 1.1145 (1.1139; [1.1126, 1.1159]) | 44.01 (44.01) | +1.5e-08 |
| 18 | 1.77964e-03 (1.77964e-03) | 34.23 (34.3) | -0.0833 (2.0e-02) | 97.89 (97.9) | 144.98 | 0.7520 (0.7517; [0.7498, 0.7525]) | 24.5 (24.5) | +3.1e-06 |
| 12 | 2.10530e-04 (2.10560e-04) | 58.44 (58.4) | -0.3563 (3.0e-02) | 81.64 (81.7) | 154.13 | 0.9116 (0.9135; [0.9139, 0.9148]) | 22.54 (22.89) | +1.2e-05 |
| 20 | 1.16958e-04 (1.16958e-04) | 36.17 (36.2) | 0, exact (4.4e-17) | — (exact) (—) | — | 1.0000 (0.9997; [0.9971, 1.0056]) | 9.906 (9.899) | +2.7e-05 |
| 27 | 6.10913e-04 (6.10914e-04) | 47.19 (47.2) | +0.2278 (5.6e-02) | 51.29 (51.3) | — | 1.0195 (1.0189; [1.0187, 1.0209]) | 27.74 (27.73) | +7.0e-12 |
| 29 | 1.34468e-03 (1.34468e-03) | 29.63 (29.7) | -0.0904 (1.0e-01) | 46.97 (47) | 61.443 | 0.7691 (0.7670; [0.7606, 0.7694]) | 15.48 (15.44) | +1.6e-05 |

The path's κ against the loop's own crossing's bracket at the smallest mass (referee round 83): inside on 9, 18, 20, 27, 29; outside on 12 (the bracket +0.25% … +0.35% off it).

κ′ against the mass (ε = 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 1e-10, 1e-11, 1e-12; † another sequence): 9: —†, 1.028, 1.082, 1.100, 1.108, 1.112, 1.113, 1.114; 18: —†, —†, 1.652†, 0.755, 0.752, 0.751, 0.752, 0.752; 12: —†, —†, 0.929, 0.925, 0.921, 0.918, 0.915, 0.914; 20: —†, —†, 0.992, 0.995, 0.998, 0.999, 0.999, 1.000; 27: —†, —†, 0.949†, 0.979†, 0.993†, 1.017, 1.018, 1.019; 29: —†, 0.643, 0.662, 0.713, 0.744, 0.758, 0.764, 0.767. The grid's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} − 1 that κ′ removes, at ε = 1e-12: 9: +0.18%, 18: +0.10%, 12: +0.13%, 20: +0.58%, 27: +0.19%, 29: +0.32%. Over every fine-grid state with a dipole pair and p < 1, p ≤ 0.28 and the landing agree in 29177 of 29425 (25094 of 25335 at ε ≤ 10⁻⁷; past τ₄ 22800 of 22978).

κ′ − κ (the path's) on the path's sequence (referee round 83; ε as above, the masses on another sequence left out), and |κ′ − κ|'s fall per decade of ε where both are at least 10⁻³ — O(ε^{1/3}) is ×10^{1/3} = ×2.154 —: 9: -0.0864, -0.0329, -0.0143, -0.0063, -0.0028, -0.0013, -0.0006 (×2.22–2.63, as ε^{0.35–0.42}) — over ε^{1/3}: -8.6, -7.1, -6.7, -6.3, -6.1, -5.9, -5.8; 18: +0.0035, +0.0000, -0.0005, -0.0004, -0.0002 — over ε^{1/3}: +1.6, +0.0, -1.1, -1.8, -2.3; 12: +0.0173, +0.0131, +0.0097, +0.0063, +0.0036, +0.0020 (×1.32–1.85, as ε^{0.12–0.27}) — over ε^{1/3}: +3.7, +6.1, +9.7, +13.6, +16.8, +19.6; 20: -0.0085, -0.0046, -0.0024, -0.0012, -0.0006, -0.0003 (×1.86–2.02, as ε^{0.27–0.31}) — over ε^{1/3}: -1.8, -2.1, -2.4, -2.5, -2.6, -2.7; 27: -0.0027, -0.0012, -0.0006 (×2.21, as ε^{0.34}) — over ε^{1/3}: -5.9, -5.7, -5.6; 29: -0.1265, -0.1072, -0.0558, -0.0248, -0.0108, -0.0047, -0.0021 (×1.18–2.30, as ε^{0.07–0.36}) — over ε^{1/3}: -12.6, -23.1, -25.9, -24.8, -23.3, -22.0, -21.1.

The loop's sequence by kernel index (SQ, referee round 83: G's loop rerun at G's outputs, the formed sets at the loop's τ₃ and τ₄ against the path's block at τ_d with its entrants j3 and j4, `kappa_from_path`): the path's 9: {4, 12} + 3 + 11; 18: {10, 16} + 4 + 15; 12: {3, 11} + 2 + 17; 20: {10, 11, 17} + 4; 27: {10, 16} + 17 + 3; 29: {4, 10} + 3 + 17. The loop follows it at 34 of 48 cells (SQ's τ₃ and τ₄ are G's at 48 of 48); elsewhere: 9 at 1e-05: τ₃ 22 {4, 11, 12}, τ₄ 26.2 {3, 4, 11, 12}; 18 at 1e-05: τ₃ 25.1 {9, 10, 16}, τ₄ 35.3 {4, 9, 10, 16}; 18 at 1e-06: τ₃ 25.1 {9, 10, 16}, τ₄ 41.3 {4, 9, 10, 16}; 18 at 1e-07: τ₃ 26.8 {9, 10, 16}, τ₄ 47.9 {4, 9, 10, 16}; 12 at 1e-05: τ₃ 24.4 {3, 10, 11}, τ₄ 69 {2, 3, 10, 11}; 12 at 1e-06: τ₃ 24.4 {3, 10, 11}, τ₄ 91 {3, 10, 11, 17}; 20 at 1e-05: τ₃ 32.3 {10, 11, 16, 17}, τ₄ —; 20 at 1e-06: τ₃ 32.3 {10, 11, 16, 17}, τ₄ —; 27 at 1e-05: τ₃ 31.6 {9, 10, 16}, τ₄ 40.5 {9, 10, 16, 17}; 27 at 1e-06: τ₃ 31.6 {9, 10, 16}, τ₄ 43.8 {9, 10, 16, 17}; 27 at 1e-07: τ₃ 34.1 {9, 10, 16}, τ₄ 45.5 {9, 10, 16, 17}; 27 at 1e-08: τ₃ 36.8 {9, 10, 16}, τ₄ 46.4 {9, 10, 16, 17}; 27 at 1e-09: τ₃ 44.8 {9, 10, 16}, τ₄ 46.9 {9, 10, 16, 17}; 29 at 1e-05: τ₃ 27.3 {3, 4, 10}, τ₄ 34.2 {3, 4, 10, 11}. The intermediate slope's proxy (beyond 1 ± 0.5 times the path's law, round 80) differs at 20 at 1e-05, 20 at 1e-06, 29 at 1e-05.


### GL. The collapse (referee rounds 68, 74): target 20 at ε = 1e-10 to τ = 1e+07

The dipole pair's distance, (W − A)/|A| and the residual, and their local exponents in τ against the law's −1/6, −1/3, −2/3 (E ∝ Δc⁸ with the centre of mass pinned):

| τ | distance | (W − A)/\|A\| | residual | exponents since the previous row: distance, W − A, residual |
|---|---|---|---|---|
| 400 | 8.804e-03 | +3.011e-03 | 1.183e-07 | — |
| 1700 | 6.767e-03 | +1.798e-03 | 4.147e-08 | -0.1818, -0.3565, -0.7244 |
| 7221 | 5.272e-03 | +1.100e-03 | 1.533e-08 | -0.1726, -0.3396, -0.6878 |
| 3.068e+04 | 4.126e-03 | +6.779e-04 | 5.766e-09 | -0.1695, -0.3347, -0.6760 |
| 1.304e+05 | 3.235e-03 | +4.188e-04 | 2.184e-09 | -0.1682, -0.3329, -0.6710 |
| 5.539e+05 | 2.539e-03 | +2.590e-04 | 8.302e-10 | -0.1675, -0.3322, -0.6686 |
| 2.354e+06 | 1.994e-03 | +1.602e-04 | 3.162e-10 | -0.1671, -0.3320, -0.6673 |
| 1e+07 | 1.566e-03 | +9.910e-05 | 1.205e-10 | -0.1669, -0.3320, -0.6667 |

The law's constant (referee round 74) over successive rows: the measured d(Δc⁻⁶)/dτ over 48aGm(W²/(w_a² + w_b²))Λ(ε), a = E/Δc⁸ (the proof's metric: the centres alone), and over it divided by 1 + C, C the slaved coefficients' share of the separation's metric (`gl_constant`); p = |w_a − w_b||Δc|/(|W|σ̄) and C at the interval's end:

| τ | p | ratio | with C | C |
|---|---|---|---|---|
| 574.3 | 0.1456 | 0.9030 | 1.0180 | 0.1190 |
| 2440 | 0.1130 | 0.9418 | 1.0139 | 0.0718 |
| 1.037e+04 | 0.0885 | 0.9656 | 1.0109 | 0.0441 |
| 4.405e+04 | 0.0695 | 0.9802 | 1.0086 | 0.0272 |
| 1.872e+05 | 0.0546 | 0.9890 | 1.0067 | 0.0168 |
| 7.953e+05 | 0.0429 | 0.9943 | 1.0053 | 0.0104 |
| 3.379e+06 | 0.0338 | 0.9975 | 1.0043 | 0.0064 |
| 1e+07 | 0.0282 | 0.9990 | 1.0038 | 0.0045 |

Over the 28 intervals the ratio is 0.9030–0.9990, with C 1.0038–1.0180; the formed kernels' energy over all kernels' 1.000001–1.000787.
The remainder with C against p: (ratio with C − 1)/p 0.123–0.134; a = E/Δc⁸ itself goes 7.893e+06 → 6.981e+06 along the run. E's next order in Δc from the local slope (referee round 80; the formed kernels' energy, the idle ones' floor left out): (d ln E/d ln Δc − 8)/p over the intervals (p their mean) 0.951–0.963, median 0.961 (all kernels' energy: 0.951–1.052, the idle kernels' floor in it), which raises the force over 8aΔc⁷ by 1 + 0.120p; with each interval's own force factor (d ln E/d ln Δc)/8 — E all kernels' energy, as in a = E/Δc⁸ (referee round 83) — the ratio with C is 1 -6.3e-05 … +1.2e-04 (the largest |·| at τ = 574.3); the mean width's share of the metric, p²/36, at most 0.00067.


### GL. The collapse (referee rounds 68, 74): target 9 at ε = 1e-10 to τ = 1e+07

The dipole pair's distance, (W − A)/|A| and the residual, and their local exponents in τ against the law's −1/6, −1/3, −2/3 (E ∝ Δc⁸ with the centre of mass pinned):

| τ | distance | (W − A)/\|A\| | residual | exponents since the previous row: distance, W − A, residual |
|---|---|---|---|---|
| 400 | 2.090e-04 | -2.848e-05 | 9.395e-12 | — |
| 1700 | 2.090e-04 | -2.848e-05 | 9.394e-12 | -0.0000, -0.0000, -0.0000 |
| 7221 | 2.090e-04 | -2.848e-05 | 9.394e-12 | -0.0000, -0.0000, -0.0000 |
| 3.068e+04 | 2.090e-04 | -2.848e-05 | 9.393e-12 | -0.0000, -0.0000, -0.0001 |
| 1.304e+05 | 2.090e-04 | -2.847e-05 | 9.389e-12 | -0.0001, -0.0002, -0.0003 |
| 5.539e+05 | 2.089e-04 | -2.844e-05 | 9.369e-12 | -0.0004, -0.0007, -0.0015 |
| 2.354e+06 | 2.084e-04 | -2.831e-05 | 9.286e-12 | -0.0015, -0.0031, -0.0062 |
| 1e+07 | 2.065e-04 | -2.780e-05 | 8.952e-12 | -0.0063, -0.0127, -0.0253 |

The law's constant (referee round 74) over successive rows: the measured d(Δc⁻⁶)/dτ over 48aGm(W²/(w_a² + w_b²))Λ(ε), a = E/Δc⁸ (the proof's metric: the centres alone), and over it divided by 1 + C, C the slaved coefficients' share of the separation's metric (`gl_constant`); p = |w_a − w_b||Δc|/(|W|σ̄) and C at the interval's end:

| τ | p | ratio | with C | C |
|---|---|---|---|---|
| 574.3 | 0.0151 | 0.9976 | 0.9989 | 0.0013 |
| 2440 | 0.0151 | 0.9978 | 0.9990 | 0.0013 |
| 1.037e+04 | 0.0151 | 0.9976 | 0.9988 | 0.0013 |
| 4.405e+04 | 0.0151 | 0.9951 | 0.9963 | 0.0013 |
| 1.872e+05 | 0.0151 | 0.9974 | 0.9987 | 0.0013 |
| 7.953e+05 | 0.0151 | 0.9974 | 0.9987 | 0.0013 |
| 3.379e+06 | 0.0150 | 0.9973 | 0.9986 | 0.0013 |
| 1e+07 | 0.0149 | 0.9974 | 0.9986 | 0.0013 |

Over the 28 intervals the ratio is 0.9914–0.9979, with C 0.9926–0.9992; the formed kernels' energy over all kernels' 0.999973–0.999995.

## H. Two features of one sign, close

Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 — 0.16 apart, 1.33 times the kernels' starting width σ₀ = 0.12 — unless stated, the moving-width path (as C's: its tolerance 1e-09, its budget 1e+06): the outcome (as C's), the most kernels of the other sign in a block along the path and the most negative coefficient, the end's formed kernels of the other sign and the smallest formed width — counting only kernels visible on the samples (largest φ there at least 1e-10; referee round 74), the end's formed ones of the other sign including the invisible in brackets —, over all the end's active kernels those narrower than the samples' spacing and those of negative width (the widths are not constrained) and how many of them are formed, the formed kernels invisible on the samples, the largest |w|, the residual:

| target | outcome | visible negative kernels along the path (most negative w) | at the end (all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | formed invisible | max\|w\| | residual | the main kernels' errors (w/A − 1) | the others: n, Σw (Σ\|w\|), on the samples, widths, \|c\| max |
|---|---|---|---|---|---|---|---|---|---|---|
| σ* 0.04, equal | exact, null cluster | 10 (-0.427) | 2 (2) | 4.00e-02 | 0, 0 (0, 0) | 0 | 1 | 1.3e-13 | 9.2e-12 (-9.2e-12), 9.2e-12 (-9.2e-12) | 3, -0.000 (0.219), mass +0.000 (window +6.7e-12), 1.0e-11 (rms 3.9e-12), 0.399–0.399, 0.00 |
| σ* 0.06, equal | cancellation | 14 (-0.610) | 7 (7) | 4.96e-02 | 6, 0 (0, 0) | 0 | 3.24 | 4.3e-10 | 1.6e-05 (-1.6e-05), 2.2e+00 (+2.2) | 7, -2.761 (2.761), mass -0.359 (window -3.6e-01), 2.3e+00 (rms 4.2e-01), 0.050–0.059, 0.16 |
| σ* 0.06, A₂ 0.99 | cancellation | 13 (-0.573) | 6 (6) | 5.19e-02 | 5, 1 (0, 0) | 0 | 2.93 | 9.6e-10 | 1.9e+00 (+1.9), 1.2e-05 (-1.2e-05) | 6, -2.287 (2.287), mass -0.306 (window -3.1e-01), 2.0e+00 (rms 3.6e-01), 0.052–0.059, 0.15 |
| σ* 0.06, A₂ 0.95 | cancellation | 14 (-0.587) | 7 (7) | 4.92e-02 | 4, 1 (0, 0) | 0 | 3.14 | 6.0e-10 | 2.1e+00 (+2.1), 1.5e-05 (-1.5e-05) | 7, -2.662 (2.662), mass -0.343 (window -3.4e-01), 2.2e+00 (rms 4.0e-01), 0.049–0.059, 0.16 |
| σ* 0.06, A₂ 0.9 | cancellation | 12 (-0.555) | 7 (7) | 5.12e-02 | 3, 1 (0, 0) | 0 | 2.84 | 8.1e-10 | 1.8e+00 (+1.8), 1.9e-01 (+0.19) | 7, -2.408 (2.408), mass -0.326 (window -3.3e-01), 1.9e+00 (rms 3.6e-01), 0.051–0.060, 0.16 |
| σ* 0.06, c₂ +0.01 | cancellation | 11 (-0.561) | 6 (7) | 5.26e-02 | 5, 1 (1, 0) | 1 | 2.85 | 4.4e-10 | 1.8e+00 (+1.8), 3.5e-06 (-3.5e-06) | 7, -2.175 (2.175), mass -0.291 (window -2.9e-01), 1.9e+00 (rms 3.5e-01), 0.000–0.058, 0.24 |
| σ* 0.06, σ₂ 0.065 | cancellation | 14 (-0.603) | 7 (7) | 4.94e-02 | 5, 2 (0, 0) | 0 | 3.16 | 6.1e-10 | 2.2e+00 (+2.2), 1.1e-05 (-1.1e-05) | 7, -2.663 (2.663), mass -0.346 (window -3.5e-01), 2.2e+00 (rms 4.0e-01), 0.049–0.059, 0.16 |
| σ* 0.04, A₂ 0.6 | exact | 0 (+0.600) | 0 (0) | 4.00e-02 | 0, 0 (0, 0) | 0 | 1 | 5.0e-17 | 2.2e-16 (+2.2e-16), 3.3e-16 (-3.3e-16) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, A₂ 0.6 | exact | 0 (+0.600) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 0 | 1 | 2.5e-16 | 1.2e-15 (+2.2e-16), 3.0e-15 (-3e-15) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |

### H2. The σ* = 0.06 paths to a tolerance of 1e-13

The same paths with the tolerance lowered (the 1e-09 ends were the tolerance's): the outcome, the entries and the last one's τ, the residuals of the last blocks, the formed kernels by the feature nearest each (their count and mass over the feature's, Σw_i|σ_i|/(Aσ*)), and H's columns:

| target | outcome | entries (last τ) | the last blocks' residuals | formed per feature (count: mass) | visible negative along the path (most negative w) | at the end (all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | max\|w\| |
|---|---|---|---|---|---|---|---|---|---|
| σ* 0.06, equal | cancellation | 24 (1.10e+15) | 2.83e-10, 2.32e-10, 2.58e-10, 1.66e-10, 1.12e-10, 8.39e-11 | 1: 1.0000, 8: 1.0000 | 14 (-0.610) | 7 (7) | 5.08e-02 | 6, 0 (0, 0) | 3.26 |
| σ* 0.06, A₂ 0.99 | stopped by the budget | 25 (7.36e+13) | 1.15e-10, 7.94e-11, 5.88e-11, 4.66e-11, 3.73e-11, 2.89e-11 | 7: 1.0000, 1: 1.0000 | 13 (-0.573) | 6 (6) | 5.43e-02 | 5, 1 (0, 0) | 2.96 |
| σ* 0.06, A₂ 0.95 | cancellation | 27 (2.29e+16) | 1.01e-10, 7.49e-11, 6.48e-11, 5.37e-11, 4.97e-11, 4.13e-11 | 8: 1.0000, 1: 1.0000 | 14 (-0.587) | 7 (7) | 5.14e-02 | 4, 1 (0, 0) | 3.18 |
| σ* 0.06, A₂ 0.9 | stopped by the budget | 20 (1.16e+12) | 5.67e-09, 8.13e-10, 5.16e-10, 4.00e-10, 2.92e-10, 2.02e-10 | 7: 1.0002, 2: 0.9998 | 12 (-0.555) | 7 (7) | 5.22e-02 | 3, 1 (0, 0) | 2.85 |
| σ* 0.06, c₂ +0.01 | cancellation | 32 (2.34e+67) | 1.71e-11, 1.70e-11, 1.77e-11, 4.21e-11, 2.60e-11, 1.49e-11 | 7: 1.0000, 2: 0.9999 | 11 (-0.561) | 6 (7) | 5.48e-02 | 6, 2 (1, 0) | 2.88 |
| σ* 0.06, σ₂ 0.065 | cancellation | 24 (2.83e+14) | 1.92e-10, 1.31e-10, 9.62e-11, 7.51e-11, 6.68e-11, 4.99e-11 | 8: 1.0000, 1: 1.0000 | 14 (-0.603) | 7 (7) | 5.15e-02 | 5, 2 (0, 0) | 3.2 |
| σ* 0.06, A₂ 0.6 | exact | 2 (3.09e+02) | 2.28e-01, 4.88e-02, 2.48e-16 | 1: 1.0000, 1: 1.0000 | 0 (+0.600) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 1 |

Why they ended (referee rounds 68, 74): σ* 0.06, equal: the stop dictionary, 21 of 21 kernels active (the dictionary exhausted), 8 of the 17 accepted polishes stopped by their budget (status 0), 0 active kernels of negative width, 3 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, A₂ 0.99: the stop budget, 20 of 21 kernels active, 10 of the 24 accepted polishes stopped by their budget (status 0), 1 active kernels of negative width, 2 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, A₂ 0.95: the stop dictionary, 21 of 21 kernels active (the dictionary exhausted), 12 of the 27 accepted polishes stopped by their budget (status 0), 1 active kernels of negative width, 3 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, A₂ 0.9: the stop budget, 17 of 21 kernels active, 6 of the 19 accepted polishes stopped by their budget (status 0), 1 active kernels of negative width, 3 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, c₂ +0.01: the stop dictionary, 21 of 21 kernels active (the dictionary exhausted), 14 of the 30 accepted polishes stopped by their budget (status 0), 2 active kernels of negative width, 4 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, σ₂ 0.065: the stop dictionary, 21 of 21 kernels active (the dictionary exhausted), 9 of the 24 accepted polishes stopped by their budget (status 0), 2 active kernels of negative width, 4 invisible on the samples (largest φ below 10⁻¹⁰); σ* 0.06, A₂ 0.6: the stop tolerance, 2 of 21 kernels active, 0 of the 2 accepted polishes stopped by their budget (status 0), 0 active kernels of negative width, 0 invisible on the samples (largest φ below 10⁻¹⁰).

### HP. H's path without the polish (referee round 68)

The moving-width path (tolerance 10⁻⁹) with the least-squares polish off, against H's (with it): the outcome, the negative kernels visible on the samples along the path and at the end (all formed at the end in brackets), the largest |w|, the residual, the active kernels and why the path stopped (referee round 74: tolerance, budget, dictionary exhausted), and the fit's split at the end (referee round 80; as HLL's):

| target | outcome (H) | visible negative along the path, at the end (H) | max\|w\| (H) | residual (H) | active | stop | the main kernels' errors (w/A − 1) | the others: n, Σw (Σ\|w\|), on the samples, widths, \|c\| max |
|---|---|---|---|---|---|---|---|---|
| σ* 0.04, equal | not fitted (exact, null cluster) | 10, 0 [0] (10, 2 [2]) | 1 (1) | 7.5e-07 (1.3e-13) | 21 | dictionary | 7.3e-05 (-7.3e-05), 7.4e-05 (-7.4e-05) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, equal | not fitted (cancellation) | 14, 12 [12] (14, 7 [7]) | 0.929 (3.24) | 9.7e-06 (4.3e-10) | 21 | dictionary | 7.1e-02 (-0.071), 7.1e-02 (-0.071) | 13, -1.046 (1.936), mass +0.184 (window +2.8e-02), 8.7e-02 (rms 2.4e-02), 0.069–1.592, 1.69 |
| σ* 0.06, A₂ 0.9 | stopped by the budget (cancellation) | 12, 6 [12] (12, 7 [7]) | 2.39 (2.84) | 1.4e-06 (8.1e-10) | 20 | budget | 1.4e+00 (+1.4), 1.4e-03 (-0.0014) | 12, -2.604 (2.604), mass -0.252 (window -2.5e-01), 1.6e+00 (rms 3.0e-01), 0.000–0.058, 0.28 |
| σ* 0.06, c₂ +0.01 | not fitted (cancellation) | 11, 6 [11] (11, 6 [7]) | 2.26 (2.85) | 9.5e-07 (4.4e-10) | 21 | dictionary | 1.3e+00 (+1.3), 5.9e-04 (-0.00059) | 11, -2.460 (2.460), mass -0.223 (window -2.2e-01), 1.5e+00 (rms 2.7e-01), 0.000–0.056, 0.28 |

### H3. The separation scan

σ* = 0.06, equal amplitudes 1, at ±c (2c in units of σ₀ in brackets; H's ±0.08 is 1.33σ₀), H's columns:

| separation | outcome | visible negative kernels along the path (most negative w) | at the end (all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | formed invisible | max\|w\| | residual | the main kernels' errors (w/A − 1) | the others: n, Σw (Σ\|w\|), on the samples, widths, \|c\| max |
|---|---|---|---|---|---|---|---|---|---|---|
| ±0.05 (0.83σ₀) | stopped by the budget | 10 (-0.240) | 10 (10) | 3.23e-02 | 0, 0 (0, 0) | 0 | 1.42 | 6.5e-06 | 1.1e+00 (+0.42), — | 10, -1.366 (1.366), mass -0.156 (window -1.6e-01), 3.3e-01 (rms 9.0e-02), 0.032–0.062, 0.26 |
| ±0.08 (1.33σ₀) | cancellation | 14 (-0.610) | 7 (7) | 4.96e-02 | 6, 0 (0, 0) | 0 | 3.24 | 4.3e-10 | 1.6e-05 (-1.6e-05), 2.2e+00 (+2.2) | 7, -2.761 (2.761), mass -0.359 (window -3.6e-01), 2.3e+00 (rms 4.2e-01), 0.050–0.059, 0.16 |
| ±0.1 (1.67σ₀) | exact | 0 (+0.000) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 0 | 1 | 3.3e-10 | 3.2e-08 (-3.2e-08), 3.2e-08 (-3.2e-08) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| ±0.12 (2.00σ₀) | exact | 0 (-0.000) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 0 | 1 | 2.5e-12 | 8.0e-11 (-8e-11), 8.0e-11 (-8e-11) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| ±0.15 (2.50σ₀) | coinciding pair | 0 (+0.305) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 0 | 1 | 5.8e-16 | 3.8e-01 (-0.38), 2.2e-16 (+2.2e-16) | 1, +0.378 (0.378), mass +0.057 (window +5.7e-02), 3.8e-01 (rms 6.9e-02), 0.060–0.060, 0.15 |
| ±0.2 (3.33σ₀) | coinciding pair | 0 (+0.231) | 0 (0) | 6.00e-02 | 0, 0 (0, 0) | 0 | 0.631 | 2.3e-15 | 3.7e-01 (-0.37), 3.7e-01 (-0.37) | 2, +0.739 (0.739), mass +0.111 (window +1.1e-01), 3.7e-01 (rms 9.5e-02), 0.060–0.060, 0.20 |

### HL. The real loop on H's targets (referee round 68: the asymmetric ones, to τ = 10⁵)

The free law from α = 0 at the masses ε (`loop_nd`, rtol 10⁻⁹) at the times τ: the formed kernels (|w| above 1% of the smallest amplitude), those visible on the samples (largest φ there at least 1e-10; referee round 74) of the other sign, the invisible formed ones, the smallest visible width, the formed kernels centred outside the samples' window, the largest |w|, the residual, and the fit's split (referee round 80; as HLL's):

| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | max\|w\| | residual | the main kernels' errors (w/A − 1) | the others: n, Σw (Σ\|w\|), on the samples, widths, \|c\| max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| σ* 0.04, equal | 1e-06 | 1600 | 9 | 6 (-0.027) | 0 | 0.036 | 0 | 0.951 | 7.8e-04 | 4.9e-02 (-0.049), 4.9e-02 (-0.049) | 7, -0.063 (0.172), mass +0.016 (window +1.6e-02), 5.5e-02 (rms 1.5e-02), 0.036–0.245, 0.39 |
| σ* 0.04, equal | 1e-06 | 6400 | 2 | 0 | 0 | 0.040 | 0 | 0.995 | 6.7e-05 | 5.0e-03 (-0.005), 5.0e-03 (-0.005) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.04, equal | 1e-06 | 25000 | 2 | 0 | 0 | 0.040 | 0 | 0.999 | 1.7e-05 | 1.3e-03 (-0.0013), 1.3e-03 (-0.0013) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.04, equal | 1e-06 | 100000 | 2 | 0 | 0 | 0.040 | 0 | 1 | 5.3e-06 | 4.3e-04 (-0.00043), 4.3e-04 (-0.00043) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.04, equal | 1e-08 | 1600 | 9 | 8 (-0.408) | 0 | 0.024 | 0 | 0.716 | 7.2e-02 | —, — | 9, -1.909 (3.341), mass +0.207 (window +2.1e-01), 7.2e-01 (rms 2.0e-01), 0.024–0.329, 0.51 |
| σ* 0.04, equal | 1e-08 | 6400 | 2 | 0 | 0 | 0.040 | 0 | 0.998 | 2.0e-05 | 1.8e-03 (-0.0018), 1.8e-03 (-0.0018) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.04, equal | 1e-08 | 25000 | 2 | 0 | 0 | 0.040 | 0 | 1 | 4.6e-06 | 4.3e-04 (-0.00043), 4.3e-04 (-0.00043) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.04, equal | 1e-08 | 100000 | 2 | 0 | 0 | 0.040 | 0 | 1 | 1.5e-06 | 1.4e-04 (-0.00014), 1.4e-04 (-0.00014) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, equal | 1e-06 | 1600 | 7 | 6 (-0.389) | 0 | 0.037 | 0 | 0.979 | 2.7e-02 | —, — | 7, -0.985 (2.943), mass +0.306 (window +3.1e-01), 9.5e-01 (rms 2.7e-01), 0.037–0.252, 0.40 |
| σ* 0.06, equal | 1e-06 | 6400 | 9 | 6 (-0.112) | 0 | 0.047 | 0 | 0.572 | 1.4e-04 | 4.3e-01 (-0.43), 4.3e-01 (-0.43) | 7, -0.043 (1.050), mass +0.159 (window +1.6e-01), 4.9e-01 (rms 1.4e-01), 0.047–0.190, 0.32 |
| σ* 0.06, equal | 1e-06 | 25000 | 11 | 8 (-0.078) | 0 | 0.046 | 0 | 0.638 | 7.0e-05 | 3.6e-01 (-0.36), 3.6e-01 (-0.36) | 9, +0.000 (0.859), mass +0.131 (window +1.3e-01), 4.2e-01 (rms 1.1e-01), 0.046–0.182, 0.39 |
| σ* 0.06, equal | 1e-06 | 100000 | 11 | 8 (-0.052) | 0 | 0.043 | 0 | 0.705 | 2.5e-05 | 3.0e-01 (-0.3), 3.0e-01 (-0.3) | 9, +0.036 (0.676), mass +0.108 (window +1.1e-01), 3.5e-01 (rms 9.4e-02), 0.043–0.168, 0.36 |
| σ* 0.06, equal | 1e-08 | 1600 | 7 | 6 (-0.389) | 0 | 0.037 | 0 | 0.979 | 2.7e-02 | —, — | 7, -0.985 (2.943), mass +0.306 (window +3.1e-01), 9.5e-01 (rms 2.7e-01), 0.037–0.252, 0.40 |
| σ* 0.06, equal | 1e-08 | 6400 | 13 | 12 (-0.456) | 0 | 0.024 | 0 | 0.958 | 2.1e-02 | —, — | 13, -2.981 (4.898), mass +0.304 (window +3.0e-01), 9.3e-01 (rms 2.7e-01), 0.024–0.438, 0.75 |
| σ* 0.06, equal | 1e-08 | 25000 | 17 | 14 (-0.211) | 0 | 0.049 | 2 | 0.799 | 3.5e-05 | 2.0e-01 (-0.2), 2.0e-01 (-0.2) | 15, -1.286 (2.135), mass +0.283 (window +7.4e-02), 2.4e-01 (rms 6.5e-02), 0.049–1.617, 1.68 |
| σ* 0.06, equal | 1e-08 | 100000 | 15 | 12 (-0.208) | 0 | 0.059 | 2 | 0.929 | 9.8e-06 | 7.1e-02 (-0.071), 7.1e-02 (-0.071) | 13, -1.046 (1.936), mass +0.184 (window +2.8e-02), 8.6e-02 (rms 2.4e-02), 0.069–1.592, 1.69 |
| σ* 0.06, A₂ 0.99 | 1e-06 | 1600 | 7 | 6 (-0.392) | 0 | 0.037 | 0 | 0.974 | 2.7e-02 | —, — | 7, -0.980 (2.928), mass +0.304 (window +3.0e-01), 9.5e-01 (rms 2.7e-01), 0.037–0.252, 0.40 |
| σ* 0.06, A₂ 0.99 | 1e-06 | 6400 | 9 | 3 (-0.448) | 3 | 0.044 | 0 | 1.58 | 2.3e-05 | 7.1e-01 (+0.58), 7.6e-03 (-0.0076) | 7, -1.204 (1.345), mass -0.139 (window -1.4e-01), 9.5e-01 (rms 1.7e-01), 0.000–0.062, 0.27 |
| σ* 0.06, A₂ 0.99 | 1e-06 | 25000 | 9 | 3 (-0.438) | 3 | 0.047 | 0 | 1.63 | 8.3e-06 | 6.5e-01 (+0.63), 3.8e-03 (-0.0038) | 7, -1.176 (1.323), mass -0.137 (window -1.4e-01), 9.5e-01 (rms 1.7e-01), 0.000–0.061, 0.27 |
| σ* 0.06, A₂ 0.99 | 1e-06 | 100000 | 9 | 3 (-0.430) | 3 | 0.049 | 0 | 1.67 | 3.3e-06 | 6.7e-01 (+0.67), 2.0e-03 (-0.002) | 7, -1.156 (1.306), mass -0.136 (window -1.4e-01), 9.4e-01 (rms 1.7e-01), 0.000–0.061, 0.27 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 1600 | 7 | 6 (-0.392) | 0 | 0.037 | 0 | 0.974 | 2.7e-02 | —, — | 7, -0.980 (2.928), mass +0.304 (window +3.0e-01), 9.5e-01 (rms 2.7e-01), 0.037–0.252, 0.40 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 6400 | 14 | 13 (-0.472) | 0 | 0.020 | 0 | 0.953 | 2.1e-02 | —, — | 14, -3.274 (5.180), mass +0.302 (window +3.0e-01), 9.2e-01 (rms 2.7e-01), 0.020–0.473, 0.83 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 25000 | 14 | 6 (-0.550) | 6 | 0.038 | 0 | 2.3 | 4.7e-06 | 1.3e+00 (+1.3), 2.7e-03 (-0.0027) | 12, -2.835 (2.835), mass -0.249 (window -2.5e-01), 1.6e+00 (rms 2.9e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 100000 | 14 | 6 (-0.539) | 6 | 0.040 | 0 | 2.38 | 1.8e-06 | 1.4e+00 (+1.4), 1.5e-03 (-0.0015) | 12, -2.784 (2.784), mass -0.253 (window -2.5e-01), 1.6e+00 (rms 3.0e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 1600 | 7 | 6 (-0.407) | 0 | 0.036 | 0 | 0.958 | 2.6e-02 | —, — | 7, -0.958 (2.875), mass +0.298 (window +3.0e-01), 9.3e-01 (rms 2.6e-01), 0.036–0.251, 0.41 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 6400 | 8 | 3 (-0.432) | 3 | 0.050 | 0 | 1.74 | 2.1e-05 | 7.4e-01 (+0.74), 3.8e-03 (-0.0033) | 6, -1.259 (1.259), mass -0.138 (window -1.4e-01), 9.1e-01 (rms 1.6e-01), 0.001–0.057, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 25000 | 8 | 3 (-0.422) | 3 | 0.052 | 0 | 1.78 | 7.3e-06 | 7.8e-01 (+0.78), 1.6e-03 (-0.0015) | 6, -1.235 (1.235), mass -0.139 (window -1.4e-01), 9.2e-01 (rms 1.6e-01), 0.001–0.058, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 100000 | 8 | 3 (-0.416) | 3 | 0.053 | 0 | 1.82 | 2.8e-06 | 8.2e-01 (+0.82), 7.7e-04 (-0.00074) | 6, -1.218 (1.218), mass -0.139 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.001–0.058, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-08 | 1600 | 7 | 6 (-0.407) | 0 | 0.036 | 0 | 0.958 | 2.6e-02 | —, — | 7, -0.958 (2.875), mass +0.298 (window +3.0e-01), 9.3e-01 (rms 2.6e-01), 0.036–0.251, 0.41 |
| σ* 0.06, A₂ 0.95 | 1e-08 | 6400 | 13 | 12 (-0.479) | 0 | 0.023 | 0 | 0.946 | 2.1e-02 | —, — | 13, -2.908 (4.800), mass +0.296 (window +3.0e-01), 9.2e-01 (rms 2.6e-01), 0.023–0.435, 0.76 |
| σ* 0.06, A₂ 0.95 | 1e-08 | 25000 | 14 | 6 (-0.547) | 6 | 0.038 | 0 | 2.29 | 4.6e-06 | 1.3e+00 (+1.3), 2.8e-03 (-0.0028) | 12, -2.622 (2.622), mass -0.247 (window -2.5e-01), 1.6e+00 (rms 2.9e-01), 0.000–0.058, 0.28 |
| σ* 0.06, A₂ 0.95 | 1e-08 | 100000 | 14 | 6 (-0.537) | 6 | 0.040 | 0 | 2.37 | 1.8e-06 | 1.4e+00 (+1.4), 1.6e-03 (-0.0016) | 12, -2.572 (2.572), mass -0.251 (window -2.5e-01), 1.6e+00 (rms 2.9e-01), 0.000–0.058, 0.28 |
| σ* 0.06, A₂ 0.9 | 1e-06 | 1600 | 8 | 4 (-0.469) | 2 | 0.001 | 0 | 1.55 | 2.9e-04 | 6.4e-01 (+0.55), 3.0e-02 (-0.023) | 6, -1.373 (1.373), mass -0.138 (window -1.4e-01), 9.0e-01 (rms 1.6e-01), 0.001–0.057, 0.26 |
| σ* 0.06, A₂ 0.9 | 1e-06 | 6400 | 8 | 3 (-0.432) | 3 | 0.050 | 0 | 1.74 | 2.0e-05 | 7.4e-01 (+0.74), 3.9e-03 (-0.0034) | 6, -1.281 (1.281), mass -0.139 (window -1.4e-01), 9.1e-01 (rms 1.7e-01), 0.001–0.057, 0.26 |
| σ* 0.06, A₂ 0.9 | 1e-06 | 25000 | 8 | 3 (-0.423) | 3 | 0.052 | 0 | 1.79 | 7.2e-06 | 7.9e-01 (+0.79), 1.7e-03 (-0.0016) | 6, -1.257 (1.257), mass -0.139 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.001–0.058, 0.26 |
| σ* 0.06, A₂ 0.9 | 1e-06 | 100000 | 8 | 3 (-0.416) | 3 | 0.053 | 0 | 1.82 | 2.8e-06 | 8.2e-01 (+0.82), 8.2e-04 (-0.00078) | 6, -1.241 (1.241), mass -0.140 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.001–0.058, 0.26 |
| σ* 0.06, A₂ 0.9 | 1e-08 | 1600 | 7 | 6 (-0.427) | 0 | 0.035 | 0 | 0.946 | 2.5e-02 | —, — | 7, -0.927 (2.818), mass +0.290 (window +2.9e-01), 9.2e-01 (rms 2.6e-01), 0.035–0.248, 0.41 |
| σ* 0.06, A₂ 0.9 | 1e-08 | 6400 | 14 | 12 (-0.586) | 0 | 0.010 | 0 | 1.25 | 6.3e-03 | —, 5.4e-01 (-0.54) | 13, -1.997 (4.503), mass +0.247 (window +2.5e-01), 9.8e-01 (rms 2.3e-01), 0.010–0.232, 0.54 |
| σ* 0.06, A₂ 0.9 | 1e-08 | 25000 | 14 | 6 (-0.542) | 6 | 0.038 | 0 | 2.3 | 4.6e-06 | 1.3e+00 (+1.3), 2.9e-03 (-0.0029) | 12, -2.666 (2.666), mass -0.247 (window -2.5e-01), 1.6e+00 (rms 2.9e-01), 0.000–0.058, 0.28 |
| σ* 0.06, A₂ 0.9 | 1e-08 | 100000 | 14 | 6 (-0.532) | 6 | 0.040 | 0 | 2.37 | 1.8e-06 | 1.4e+00 (+1.4), 1.6e-03 (-0.0016) | 12, -2.616 (2.616), mass -0.251 (window -2.5e-01), 1.6e+00 (rms 3.0e-01), 0.000–0.058, 0.28 |
| σ* 0.06, c₂ +0.01 | 1e-06 | 1600 | 2 | 0 | 0 | 0.060 | 0 | 1 | 4.6e-09 | 1.9e-08 (-1.9e-08), 8.2e-09 (-3.4e-09) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, c₂ +0.01 | 1e-06 | 6400 | 2 | 0 | 0 | 0.060 | 0 | 1 | 4.6e-09 | 1.9e-08 (-1.9e-08), 8.2e-09 (-3.4e-09) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, c₂ +0.01 | 1e-06 | 25000 | 2 | 0 | 0 | 0.060 | 0 | 1 | 4.6e-09 | 1.9e-08 (-1.9e-08), 8.2e-09 (-3.4e-09) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, c₂ +0.01 | 1e-06 | 100000 | 2 | 0 | 0 | 0.060 | 0 | 1 | 4.6e-09 | 1.9e-08 (-1.9e-08), 8.2e-09 (-3.4e-09) | 0, +0.000 (0.000), mass +0.000 (window +0.0e+00), 0.0e+00 (rms 0.0e+00), —, — |
| σ* 0.06, c₂ +0.01 | 1e-08 | 1600 | 7 | 6 (-0.391) | 0 | 0.038 | 0 | 0.934 | 3.5e-02 | —, — | 7, -1.073 (2.941), mass +0.307 (window +3.1e-01), 9.2e-01 (rms 2.7e-01), 0.038–0.272, 0.43 |
| σ* 0.06, c₂ +0.01 | 1e-08 | 6400 | 13 | 6 (-0.504) | 5 | 0.035 | 0 | 2.02 | 1.3e-05 | 1.0e+00 (+1), 3.0e-03 (-0.0027) | 11, -2.547 (2.547), mass -0.200 (window -2.0e-01), 1.3e+00 (rms 2.4e-01), 0.000–0.053, 0.28 |
| σ* 0.06, c₂ +0.01 | 1e-08 | 25000 | 13 | 6 (-0.491) | 5 | 0.038 | 0 | 2.12 | 3.9e-06 | 1.1e+00 (+1.1), 1.4e-03 (-0.0013) | 11, -2.488 (2.488), mass -0.208 (window -2.1e-01), 1.4e+00 (rms 2.5e-01), 0.000–0.054, 0.28 |
| σ* 0.06, c₂ +0.01 | 1e-08 | 100000 | 13 | 6 (-0.482) | 5 | 0.040 | 0 | 2.19 | 1.6e-06 | 1.2e+00 (+1.2), 7.8e-04 (-0.00078) | 11, -2.449 (2.449), mass -0.213 (window -2.1e-01), 1.4e+00 (rms 2.5e-01), 0.000–0.055, 0.28 |
| σ* 0.06, σ₂ 0.065 | 1e-06 | 1600 | 6 | 3 (-0.402) | 1 | 0.001 | 0 | 1.51 | 7.4e-05 | 5.1e-01 (+0.51), 6.6e-03 (-0.0045) | 4, -0.832 (0.832), mass -0.097 (window -9.7e-02), 6.3e-01 (rms 1.1e-01), 0.001–0.057, 0.27 |
| σ* 0.06, σ₂ 0.065 | 1e-06 | 6400 | 6 | 2 (-0.387) | 2 | 0.056 | 0 | 1.58 | 1.0e-05 | 5.8e-01 (+0.58), 1.2e-03 (-0.00087) | 4, -0.806 (0.806), mass -0.097 (window -9.7e-02), 6.3e-01 (rms 1.1e-01), 0.001–0.058, 0.27 |
| σ* 0.06, σ₂ 0.065 | 1e-06 | 25000 | 6 | 2 (-0.382) | 2 | 0.057 | 0 | 1.6 | 3.2e-06 | 6.0e-01 (+0.6), 4.4e-04 (-0.00033) | 4, -0.797 (0.797), mass -0.097 (window -9.7e-02), 6.2e-01 (rms 1.1e-01), 0.001–0.059, 0.27 |
| σ* 0.06, σ₂ 0.065 | 1e-06 | 100000 | 6 | 2 (-0.379) | 2 | 0.058 | 0 | 1.62 | 1.2e-06 | 6.2e-01 (+0.62), 1.7e-04 (-0.00014) | 4, -0.791 (0.791), mass -0.097 (window -9.7e-02), 6.2e-01 (rms 1.1e-01), 0.001–0.059, 0.27 |
| σ* 0.06, σ₂ 0.065 | 1e-08 | 1600 | 7 | 6 (-0.403) | 0 | 0.038 | 0 | 1.01 | 2.3e-02 | —, — | 7, -0.915 (2.937), mass +0.318 (window +3.2e-01), 9.6e-01 (rms 2.7e-01), 0.038–0.248, 0.40 |
| σ* 0.06, σ₂ 0.065 | 1e-08 | 6400 | 13 | 12 (-0.474) | 0 | 0.025 | 0 | 0.991 | 1.7e-02 | —, — | 13, -2.891 (4.874), mass +0.316 (window +3.2e-01), 9.4e-01 (rms 2.7e-01), 0.025–0.426, 0.74 |
| σ* 0.06, σ₂ 0.065 | 1e-08 | 25000 | 13 | 7 (-0.645) | 4 | 0.034 | 0 | 2.59 | 5.8e-06 | 1.6e+00 (+1.6), 2.8e-03 (-0.0028) | 11, -3.255 (3.255), mass -0.307 (window -3.1e-01), 1.9e+00 (rms 3.5e-01), 0.000–0.058, 0.28 |
| σ* 0.06, σ₂ 0.065 | 1e-08 | 100000 | 13 | 7 (-0.631) | 4 | 0.036 | 0 | 2.69 | 2.1e-06 | 1.7e+00 (+1.7), 1.5e-03 (-0.0015) | 11, -3.183 (3.183), mass -0.313 (window -3.1e-01), 1.9e+00 (rms 3.5e-01), 0.000–0.058, 0.28 |

### HLL. HL's cells to τ = 10⁹ (referee rounds 74, 80, 83)

HL's columns, and the fit's split (`fit_split`, referee round 83: the main kernels the assignment minimizing the largest error): per feature its main kernel's largest relative error in (c, σ, w) and its signed w/A − 1 (inflated: a cancellation; deficient: a split); the other formed kernels — their number, Σw and Σ\|w\|, their sum's largest value on the samples over max f* (a null cluster where it is small against the residual), their widths and their centres' largest \|c\|:

| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | max\|w\| | residual | the main kernels' errors (w/A − 1) | the others: n, Σw (Σ\|w\|), on the samples, widths, \|c\| max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| σ* 0.06, equal | 1e-06 | 100000 | 11 | 8 (-0.052) | 0 | 0.043 | 0 | 0.705 | 2.5e-05 | 3.0e-01 (-0.3), 3.0e-01 (-0.3) | 9, +0.036 (0.676), mass +0.108 (window +1.1e-01), 3.5e-01 (rms 9.4e-02), 0.043–0.168, 0.36 |
| σ* 0.06, equal | 1e-06 | 1e+06 | 11 | 8 (-0.036) | 0 | 0.041 | 0 | 0.783 | 7.1e-06 | 2.2e-01 (-0.22), 2.2e-01 (-0.22) | 9, +0.030 (0.494), mass +0.082 (window +8.2e-02), 2.5e-01 (rms 7.0e-02), 0.041–0.169, 0.32 |
| σ* 0.06, equal | 1e-06 | 1e+07 | 11 | 8 (-0.025) | 0 | 0.039 | 0 | 0.839 | 2.3e-06 | 1.6e-01 (-0.16), 1.6e-01 (-0.16) | 9, +0.025 (0.362), mass +0.064 (window +6.4e-02), 1.9e-01 (rms 5.2e-02), 0.039–0.175, 0.29 |
| σ* 0.06, equal | 1e-06 | 1e+08 | 11 | 8 (-0.019) | 0 | 0.040 | 0 | 0.88 | 6.5e-07 | 1.2e-01 (-0.12), 1.2e-01 (-0.12) | 9, +0.016 (0.271), mass +0.048 (window +4.8e-02), 1.4e-01 (rms 3.9e-02), 0.040–0.180, 0.32 |
| σ* 0.06, equal | 1e-06 | 1e+09 | 9 | 6 (-0.014) | 0 | 0.043 | 0 | 0.909 | 1.8e-07 | 9.1e-02 (-0.091), 9.1e-02 (-0.091) | 7, +0.033 (0.185), mass +0.040 (window +4.0e-02), 1.1e-01 (rms 3.1e-02), 0.043–0.184, 0.29 |
| σ* 0.06, equal | 1e-08 | 100000 | 15 | 12 (-0.208) | 0 | 0.059 | 2 | 0.929 | 9.8e-06 | 7.1e-02 (-0.071), 7.1e-02 (-0.071) | 13, -1.046 (1.936), mass +0.184 (window +2.8e-02), 8.6e-02 (rms 2.4e-02), 0.069–1.592, 1.69 |
| σ* 0.06, equal | 1e-08 | 1e+06 | 13 | 10 (-0.206) | 0 | 0.060 | 2 | 0.996 | 7.8e-07 | 4.3e-03 (-0.0043), 4.3e-03 (-0.0043) | 11, -0.995 (1.903), mass +0.164 (window +2.1e-03), 5.3e-03 (rms 1.6e-03), 0.204–1.718, 1.84 |
| σ* 0.06, equal | 1e-08 | 1e+07 | 13 | 10 (-0.205) | 0 | 0.060 | 2 | 1 | 1.2e-07 | 5.0e-04 (-0.0005), 5.0e-04 (-0.0005) | 11, -0.984 (1.902), mass +0.169 (window +2.7e-04), 6.3e-04 (rms 1.9e-04), 0.247–1.800, 1.94 |
| σ* 0.06, equal | 1e-08 | 1e+08 | 13 | 10 (-0.205) | 0 | 0.060 | 2 | 1 | 2.3e-08 | 7.6e-05 (-7.6e-05), 7.6e-05 (-7.6e-05) | 11, -0.974 (1.902), mass +0.175 (window +4.6e-05), 1.0e-04 (rms 3.2e-05), 0.282–1.878, 2.04 |
| σ* 0.06, equal | 1e-08 | 1e+09 | 13 | 10 (-0.204) | 0 | 0.060 | 2 | 1 | 5.4e-09 | 1.4e-05 (-1.4e-05), 1.4e-05 (-1.4e-05) | 11, -0.964 (1.902), mass +0.179 (window +1.0e-05), 2.2e-05 (rms 7.0e-06), 0.314–1.953, 2.13 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 100000 | 14 | 6 (-0.539) | 6 | 0.040 | 0 | 2.38 | 1.8e-06 | 1.4e+00 (+1.4), 1.5e-03 (-0.0015) | 12, -2.784 (2.784), mass -0.253 (window -2.5e-01), 1.6e+00 (rms 3.0e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 1e+06 | 14 | 6 (-0.526) | 6 | 0.043 | 0 | 2.48 | 4.4e-07 | 1.5e+00 (+1.5), 6.5e-04 (-0.00065) | 12, -2.722 (2.722), mass -0.259 (window -2.6e-01), 1.7e+00 (rms 3.0e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 1e+07 | 14 | 6 (-0.515) | 6 | 0.045 | 0 | 2.55 | 1.2e-07 | 1.6e+00 (+1.6), 2.8e-04 (-0.00028) | 12, -2.671 (2.671), mass -0.264 (window -2.6e-01), 1.7e+00 (rms 3.1e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 1e+08 | 14 | 6 (-0.506) | 6 | 0.047 | 0 | 2.62 | 3.0e-08 | 1.6e+00 (+1.6), 1.2e-04 (-0.00012) | 12, -2.630 (2.630), mass -0.268 (window -2.7e-01), 1.7e+00 (rms 3.1e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.99 | 1e-08 | 1e+09 | 14 | 6 (-0.498) | 6 | 0.049 | 0 | 2.67 | 7.9e-09 | 1.7e+00 (+1.7), 5.0e-05 (-5e-05) | 12, -2.595 (2.595), mass -0.271 (window -2.7e-01), 1.8e+00 (rms 3.2e-01), 0.000–0.059, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 100000 | 8 | 3 (-0.416) | 3 | 0.053 | 0 | 1.82 | 2.8e-06 | 8.2e-01 (+0.82), 7.7e-04 (-0.00074) | 6, -1.218 (1.218), mass -0.139 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.001–0.058, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 1e+06 | 8 | 3 (-0.407) | 3 | 0.055 | 0 | 1.86 | 6.1e-07 | 8.6e-01 (+0.86), 2.3e-04 (-0.00023) | 6, -1.198 (1.198), mass -0.140 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.000–0.058, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 1e+07 | 8 | 3 (-0.402) | 3 | 0.056 | 0 | 1.89 | 1.3e-07 | 8.9e-01 (+0.89), 6.7e-05 (-6.7e-05) | 6, -1.184 (1.184), mass -0.141 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.000–0.059, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 1e+08 | 8 | 3 (-0.398) | 3 | 0.057 | 0 | 1.91 | 2.9e-08 | 9.1e-01 (+0.91), 1.9e-05 (-1.9e-05) | 6, -1.175 (1.175), mass -0.142 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.000–0.059, 0.26 |
| σ* 0.06, A₂ 0.95 | 1e-06 | 1e+09 | 8 | 3 (-0.395) | 3 | 0.058 | 0 | 1.92 | 6.2e-09 | 9.2e-01 (+0.92), 5.3e-06 (-5.3e-06) | 6, -1.168 (1.168), mass -0.142 (window -1.4e-01), 9.2e-01 (rms 1.7e-01), 0.000–0.059, 0.26 |

The others' sum on the samples (referee round 83): its fall per decade of τ, and at the last τ against the residual's rms, both over max f* (a null cluster where, the main kernels exact, the others' rms is below 1e-08 of max f*: referee round 97): σ* 0.06, equal at 1e-06: ×1.36, 1.36, 1.34, 1.32; at τ = 1e+09 1.1e-01 against 1.8e-07; the others' rms over the residual's 3.94e+03, 1.01e+04, 2.33e+04, 6.22e+04, 1.71e+05 (1% of max f* after 8.5 more decades at the last decade's rate); σ* 0.06, equal at 1e-08: ×16.2, 8.41, 6.14, 4.68; at τ = 1e+09 2.2e-05 against 5.3e-09; the others' rms over the residual's 2.49e+03, 2.06e+03, 1.68e+03, 1.46e+03, 1.33e+03; σ* 0.06, A₂ 0.99 at 1e-08: ×0.975, 0.978, 0.984, 0.985; at τ = 1e+09 1.8e+00 against 7.7e-09; the others' rms over the residual's 1.71e+05, 7.04e+05, 2.76e+06, 1.07e+07, 4.14e+07; σ* 0.06, A₂ 0.95 at 1e-06: ×0.999, 0.996, 1, 0.996; at τ = 1e+09 9.2e-01 against 6.0e-09; the others' rms over the residual's 6.08e+04, 2.81e+05, 1.29e+06, 5.88e+06, 2.78e+07.


The splits' classes (referee rounds 97, 100): over the 94 splits above, the main kernels' (c, σ) errors at most 1.15 (H3 σ* 0.06, ±0.05, where the other feature's only kernel of its sign is this one, left to it by distinctness) and their w/A − 1 -0.543 … +2.24; 30 features unassigned, 29 of them with an unused kernel of their sign, at least 2.81 off in (c, σ) (HL σ* 0.06, σ₂ 0.065 1e-08 τ 1600) — the cap 1.5 on the (c, σ) error alone (w free) lies in that gap, and every split's assignment is the same for every cap in (1.1499, 2.8144) (`cap_interval`, rounded inward; where a same-sign kernel's (c, σ) error is crossed); null clusters (`null_cluster`, `outcome`'s test too: the residual below 1e-08, every main kernel within 1e-06 in c, σ and w, the others' rms below 1e-08 of max f*; referee round 100): H σ* 0.04, equal. On the targets symmetric about 0 the formed kernels' distance from their mirror image (`mirror_asym`, a bottleneck matching; half of it from the reflection's fixed-point subspace; referee rounds 103, 106, 112): H σ* 0.04, equal 2.2e-11 from the mirror image (1.1e-11 from the subspace; all 13 kernels 2.6e-10), the formed kernels' function odd by 1.1e-16 on and between the samples, the whole state's by 1.1e-16 on the samples (the residual's: f* is even) and 1.7e-16 between them; the 13 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.0268) odd by 1.7e-16 on and between the samples (5 formed, the narrowest of them 0.04; 8 idle); the formed threshold's margin: the largest unformed |w| 4.6e-12 (over max A, against the threshold 0.01 of min A), the smallest formed 0.055 (over max A); H σ* 0.06, equal 0.24 from the mirror image (0.12 from the subspace; all 17 kernels 0.48), the formed kernels' function odd by 5.5e-09 on and between the samples, the whole state's by 1.2e-09 on the samples (the residual's: f* is even) and 4.9e-03 between them; the 11 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.041, the widest narrower one 0.001) odd by 1.2e-09 on and between the samples (9 formed, the narrowest of them 0.0496; 2 idle); the formed threshold's margin: the largest unformed |w| 0.0098 (over max A, against the threshold 0.01 of min A), the smallest formed 0.06 (over max A); H3 σ* 0.06, ±0.05 0.0016 from the mirror image (0.00082 from the subspace; all 18 kernels 0.42), the formed kernels' function odd by 1.5e-04 on and between the samples, the whole state's by 1.3e-06 on the samples (the residual's: f* is even) and 1.3e-06 between them; the 18 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.0323) odd by 1.3e-06 on and between the samples (11 formed, the narrowest of them 0.0323; 7 idle); the formed threshold's margin: the largest unformed |w| 0.0072 (over max A, against the threshold 0.01 of min A), the smallest formed 0.036 (over max A); H3 σ* 0.06, ±0.1 6.9e-18 from the mirror image (3.5e-18 from the subspace; all 3 kernels 5.7e-10), the formed kernels' function odd by 1.1e-16 on and between the samples, the whole state's by 1.1e-16 on the samples (the residual's: f* is even) and 2.2e-16 between them; the 3 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.06) odd by 2.2e-16 on and between the samples (2 formed, the narrowest of them 0.06; 1 idle); the formed threshold's margin: the largest unformed |w| 5e-08 (over max A, against the threshold 0.01 of min A), the smallest formed 1 (over max A); H3 σ* 0.06, ±0.12 5.1e-13 from the mirror image (2.6e-13 from the subspace; all 12 kernels 0.00075), the formed kernels' function odd by 2.6e-13 on and between the samples, the whole state's by 2.3e-14 on the samples (the residual's: f* is even) and 2.4e-14 between them; the 12 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.06) odd by 2.4e-14 on and between the samples (2 formed, the narrowest of them 0.06; 10 idle); the formed threshold's margin: the largest unformed |w| 1.8e-07 (over max A, against the threshold 0.01 of min A), the smallest formed 1 (over max A); H3 σ* 0.06, ±0.15 0.3 from the mirror image (0.15 from the subspace; all 3 kernels 0.3), the formed kernels' function odd by 1.4e-15 on and between the samples, the whole state's by 1.4e-15 on the samples (the residual's: f* is even) and 1.4e-15 between them; the 3 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.06) odd by 1.4e-15 on and between the samples (3 formed, the narrowest of them 0.06; 0 idle); H3 σ* 0.06, ±0.2 1.5e-06 from the mirror image (7.6e-07 from the subspace; all 4 kernels 1.5e-06), the formed kernels' function odd by 5.0e-16 on and between the samples, the whole state's by 4.4e-16 on the samples (the residual's: f* is even) and 5.0e-16 between them; the 4 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.06) odd by 5.0e-16 on and between the samples (4 formed, the narrowest of them 0.06; 0 idle); HP σ* 0.04, equal 1.4e-07 from the mirror image (7.2e-08 from the subspace; all 21 kernels 0.012), the formed kernels' function odd by 1.4e-07 on and between the samples, the whole state's by 1.5e-07 on the samples (the residual's: f* is even) and 1.5e-07 between them; the 21 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.028) odd by 1.5e-07 on and between the samples (2 formed, the narrowest of them 0.04; 19 idle); the formed threshold's margin: the largest unformed |w| 7.9e-05 (over max A, against the threshold 0.01 of min A), the smallest formed 1 (over max A); HP σ* 0.06, equal 2.9e-08 from the mirror image (1.4e-08 from the subspace; all 21 kernels 3.8e-05), the formed kernels' function odd by 9.0e-10 on and between the samples, the whole state's by 1.0e-12 on the samples (the residual's: f* is even) and 1.0e-12 between them; the 21 kernels at least as wide as the samples' spacing (σ ≥ h; the narrowest 0.0504) odd by 1.0e-12 on and between the samples (15 formed, the narrowest of them 0.059; 6 idle); the formed threshold's margin: the largest unformed |w| 0.0096 (over max A, against the threshold 0.01 of min A), the smallest formed 0.011 (over max A); the loops' formed kernels (HL, HLL) at most 2.8e-07.

Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.
