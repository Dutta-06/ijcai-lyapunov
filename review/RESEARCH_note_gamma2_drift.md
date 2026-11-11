<!-- Research note (2026-11-05): why Part K of results/chaos/gauge_noise.md saw only part of Katzenberger's O(γ²τ) drift of the masses
under the centres-and-widths noise (D‴.40; Paper 14 §7, §9.5). Computation: experiments/exp_gauge_noise.py Part M (raw key M); every
number below is printed in results/chaos/gauge_noise.md, Part M (§4–5's finite-τ numbers: Part N, raw key N, added 2026-11-05). Not yet
refereed. -->

# D‴.40: the missing O(γ²τ) drift of the deterministic-amplitude loop

Research note, 2026-11-05. Source of every number: `results/chaos/gauge_noise.md`, Part M (sections A–L and raw keys A–L unchanged).
Drifts are per γ²τ and listed per class (W_1, W_2, W_3) unless stated.

## 1. The question

In the variant 'nonlinear only' the centres and widths are dithered and the amplitudes are noise-free (γ = 0.01, τ = 0.0001, the box).
Proposition D‴.40(ii) (Paper 14, Prop. 7.3(a)) says that the masses then have no diffusion. They keep a drift of order γ²τ,
(γ/m)²Σ_lE[φ_l dφ_k/ds], which comes from an antisymmetric correlation of the positions and would vanish for a reversible law. Katzenberger's
limit (`limit_drift`) gives it a size. Over the second half of Part E's runs the drift path is −6.04·10⁻⁵, +5.01·10⁻⁵, +1.03·10⁻⁵.

Part K subtracted from the measured masses two terms:
- (1) the positions' offset from the fit, at first order about the fit at the current masses;
- (2) the Euler chain's O(Δt) effect, from runs at Δt/2 on the same Brownian path.

What was left is far below that path: at Part M(b)'s records, −9.36·10⁻⁶ ± 1.1·10⁻⁵, +7.37·10⁻⁶ ± 3.3·10⁻⁶, +1.98·10⁻⁶ ± 9.6·10⁻⁶,
0.15 ± 0.19, 0.15 ± 0.07, 0.19 ± 0.93 of it. The open question was:
is the drift absent from the loop, or did Part K fail to see it?

## 2. Hypotheses and tests

Each hypothesis below has its test in Part M.

**H1. The prediction is wrong** (Part M(a)). The prediction is computed three independent ways at the masses 1/3:
- from the closed form of the linearised positions' antisymmetric correlation (γ → 0);
- from `limit_drift`;
- from the second differences of the limit map Φ, the noise-free loop integrated by RK4 (`m_limit`, `m_fd`).

| source | γ | drift per γ²τ |
|---|---|---|
| closed form | 0 | −0.01368, +0.01237, +0.00131 |
| `limit_drift` | 10⁻⁵ | −0.01387, +0.01263, +0.00125 |
| `limit_drift` | 0.01 | −0.00803, +0.00662, +0.00141 |
| Φ's second differences, Richardson | 0.01 | −0.00803, +0.00662, +0.00141 |

At γ = 0.01, Φ's second differences agree with `limit_drift` to a largest relative difference of 0.0016. At γ = 0.001 the difference is
0.0242, in the small third class, where the limit map's own error shows. **Excluded:** the prediction stands.

**H2. Reversibility** (Part M(a), (d)). If the positions' law at fixed masses were close to reversible, the drift would be small. The
antisymmetric correlation's share, E[δα_lδα̇_k]/√(E[δα_l²]E[δα̇_k²]), is indeed small: at most 0.0341. But the drift is that small share
times large factors.

In simulation at fixed masses 1/3 (`m_fixed`, two Euler chains on one Brownian path, Richardson), the measured drift against the closed
form is:

| τ | Richardson over the closed form | z against it |
|---|---|---|
| 3·10⁻⁵ | 1.08 ± 0.20, 1.07 ± 0.13 | −0.4, +0.5 |
| 10⁻⁴ | 0.86 ± 0.08, 0.93 ± 0.10 | +1.8, −0.8 |
| 0.0003 | −0.01 ± 0.66, 0.08 ± 0.94 | +1.5, −1.0 |

(First two classes; the third is small and noisy.) The chain's symmetric part at τ = 10⁻⁴ is −0.00177 ± 0.00021 for the Δt chain,
against the closed form's −0.00180. **Not the cause:** the correlation is there at the predicted size.

**H3. No slow mass mode at γ = 0.01** (Part M(c)). If the masses were not a slow mode at γ = 0.01, Katzenberger's limit would not apply.
New runs (48 seeds, 200000 slow units from the exact fit) at γ = 0.01 and 0.003, each as two chains at Δt and Δt/2 on one Brownian path,
measure the rate of Y = Φ(state). The rate per γ²τ at γ = 0.003 over that at γ = 0.01 is:
- Richardson: 1.44 ± 0.02, 1.61 ± 0.04, 0.49 ± 0.66;
- Katzenberger's: 1.43, 1.52, 0.99.

The first transient, Y(20000) − 1/3, is −1.27·10⁻⁶ ± 4.0·10⁻⁷, +1.07·10⁻⁶ ± 3.4·10⁻⁷, +1.99·10⁻⁷ ± 5.7·10⁻⁷ at γ = 0.01. **Excluded:**
the loop's rate follows the finite-γ prediction's γ dependence (the first class closely, the second somewhat above it).

**H4. The estimator** (Part M(b)). Part K's offset (1) is first order about the fit at the current masses, and it ignores the fit's motion
with the masses. When the masses move, the fit moves. The fit's amplitudes then change, and that change moves the masses again. With
M = (γ/m)(I − W1ᵀ)diag(−h_k²/(G_αW_k)) and h_k = a_k/a′_k, the true offset is W − Y = (I − M)⁻¹(1). At the masses 1/3, −diag M =
0.3135, 0.0560, 0.0938: the first class's amplitude is near its saturation.

On Part K's 12 seeds (every 5th record of the second half):

| quantity | value (± s.e.) |
|---|---|
| (1), first order | −1.23·10⁻⁴ ± 3.1·10⁻⁵, +1.46·10⁻⁴ ± 1.4·10⁻⁵, −2.28·10⁻⁵ ± 2.4·10⁻⁵ |
| W − Y, exact | −8.96·10⁻⁵ ± 2.3·10⁻⁵, +1.21·10⁻⁴ ± 1.1·10⁻⁵, −3.18·10⁻⁵ ± 1.9·10⁻⁵ |
| (I − M)⁻¹(1) | −8.82·10⁻⁵ ± 2.3·10⁻⁵, +1.21·10⁻⁴ ± 1.1·10⁻⁵, −3.26·10⁻⁵ ± 1.9·10⁻⁵ |
| Y − 1/3 − (2) | −4.25·10⁻⁵ ± 5.2·10⁻⁶, +3.15·10⁻⁵ ± 2.8·10⁻⁶, +1.10·10⁻⁵ ± 5.5·10⁻⁶ |

Per record, |W − Y − (I − M)⁻¹(1)| is at most 1.4·10⁻⁵, while |W − Y − (1)| reaches 8.6·10⁻⁴. Against Katzenberger's path, Y − 1/3 − (2)
is 0.70, 0.63, 1.07 of it (± 0.09, 0.06, 0.53), with z = +3.4, −6.6, +0.1. **The main cause.**

## 3. The answer

**The loop carries Katzenberger's O(γ²τ) drift, somewhat below its predicted size.** The Richardson-extrapolated rate of the
limit map (Part M(c)) against Katzenberger's drift at 1/3:

| γ | Richardson rate per γ²τ | Katzenberger's | ratio | z |
|---|---|---|---|---|
| 0.01 | −0.00633 ± 0.00065, +0.00538 ± 0.00059, +0.00095 ± 0.00069 | −0.00803, +0.00662, +0.00141 | 0.79 ± 0.08, 0.81 ± 0.09, 0.67 ± 0.49 | +2.6, −2.1, −0.7 |
| 0.003 | −0.00912 ± 0.00100, +0.00866 ± 0.00082, +0.00046 ± 0.00096 | −0.01148, +0.01008, +0.00140 | 0.79 ± 0.09, 0.86 ± 0.08, 0.33 ± 0.68 | +2.4, −1.7, −1.0 |

The two γ use the same seeds and the same Brownian draws, so their agreement is not independent evidence.

At Part E's step Δt, the Euler step's own effect is as large as the drift itself. At γ = 0.01 the Δt chain alone has the rate −0.01363 ± 0.00070,
+0.01609 ± 0.00064, −0.00246 ± 0.00069. Its difference from Richardson, r(Δt) − Richardson, is −0.00730 ± 0.00009, +0.01071 ± 0.00010,
−0.00341 ± 0.00006.

**Part K's shortfall, term by term** (Part M(e)). Katzenberger's path over the second half is the sum of four terms (their sum differs from
it by at most 3.0·10⁻²⁰). As shares of the path, first two classes:

| term | share of the path |
|---|---|
| Part K's secular displacement (at (b)'s records) | 0.15 ± 0.19, 0.15 ± 0.07 |
| the first-order offset's error, (1) − (W − Y) | 0.55 ± 0.14, 0.48 ± 0.07 |
| the excess of Part K's step correction (2) over the loop's own step effect from (c) | 0.09 ± 0.02, 0.08 ± 0.03 |
| the rest | 0.20 ± 0.09, 0.29 ± 0.06 |

In other words:
- About half of the path was hidden by the first-order offset (1). Part K's (2) as a rate is −0.00806, +0.01126, −0.00320, against the
  loop's −0.00730, +0.01071, −0.00341; it over-corrects by 0.09 ± 0.02, 0.08 ± 0.03 of the path.
- The rest is not a fluctuation of Part K's 12 seeds. With (c)'s step effect in place of (2), Y − 1/3 is 0.80 ± 0.09, 0.71 ± 0.06,
  1.22 ± 0.55 of the path. Against (c)'s Richardson ratio at γ = 0.01 this gives z = +0.1, −0.9, +0.7. The rest is the loop's own
  shortfall, the same as in the new runs.

This revises the preliminary reading, which took the rest for a fluctuation.

**Scaling.**
- γ: the loop's rate per γ²τ grows from γ = 0.01 to 0.003 by 1.44 ± 0.02, 1.61 ± 0.04, as Katzenberger's does (1.43, 1.52). The drift is
  O(γ²τ), with the finite-γ factor that `limit_drift` carries.
- τ: at fixed masses the drift per γ²τ is constant in τ within errors at τ = 3·10⁻⁵ and 10⁻⁴ (table in H2). At 0.0003 it is undetermined
  (−0.01 ± 0.66 of the closed form, z = +1.5 against it). [Part N, 2026-11-05: see §4; with 144 more seeds it is 0.19 ± 0.32,
  0.49 ± 0.31.]
- The loop against fixed masses: the loop's ratio to Katzenberger is consistent with the fixed-mass ratio to the closed form at the same
  τ = 10⁻⁴ (0.86, 0.93, 0.19). The z are −0.6, −0.8, +0.4 at γ = 0.01 and −0.5, −0.5, +0.1 at γ = 0.003.

## 4. What remains open

- **The loop's shortfall.** In the first two classes the loop's rate is 0.79 ± 0.08, 0.81 ± 0.09 of Katzenberger's at γ = 0.01
  (z = +2.6, −2.1) and 0.79 ± 0.09, 0.86 ± 0.08 at γ = 0.003 (z = +2.4, −1.7), on shared draws. The fixed-mass correlation at the same τ has a similar, less significant shortfall (0.86 ± 0.08,
  z = +1.8). The likely reading is a finite-τ correction of the positions' law beyond its linearisation. The O(γ²τ) prediction is
  leading order in τ, and the first class's amplitude is near its saturation. That reading is not established:
  - the τ = 3·10⁻⁵ point (1.08 ± 0.20) is consistent with both readings;
  - the τ = 0.0003 point is too noisy.

  Settling it needs more seeds at fixed masses at τ = 10⁻⁴ and 0.0003, or the loop at a smaller τ.

  **[After Part N, 2026-11-05.]** Part N ran `m_fixed` (unchanged) at new seeds: 32 at τ = 0.0001 and 144 at 0.0003, with Part M(d)'s
  run lengths so that the seeds pool, and 16 at 0.001. Ratios are Richardson over the closed form, first two classes.
  - **Not noise.** At τ = 0.0001 the new seeds alone give 0.79 ± 0.09, 0.79 ± 0.07 (against Part M(d)'s seeds z = −0.6, −1.1). Pooled,
    the ratio is 0.82 ± 0.06, 0.86 ± 0.06, with z = −3.0, −2.3 against 1.
  - **At τ = 0.0003** the pooled ratio is 0.19 ± 0.32, 0.49 ± 0.31 (z = −2.5, −1.6 against 1). This point stays wide: the seeds' spread
    of the ratio is 4.06, 3.96 there (100000 slow units), against 0.47, 0.50 at 0.0001 (200000 slow units).
  - **The line.** Weighted least squares through the pooled points at τ = 3·10⁻⁵, 10⁻⁴, 3·10⁻⁴ gives:
    - the slope −0.33 ± 0.14, −0.23 ± 0.12 per 10⁻⁴ of τ (z = −2.3, −1.9 against 0);
    - the limit τ → 0 1.15 ± 0.15, 1.10 ± 0.12 (z = +1.0, +0.8 against 1);
    - χ² 0.02, 0.19 with one degree of freedom.

    A constant ratio gives 0.82 ± 0.06, 0.88 ± 0.05 (z = −3.2, −2.1 against 1), with χ² 5.42, 3.91 on two degrees of freedom
    (p = 0.07, 0.14).
  - **τ = 0.001** (16 seeds, outside the fit) is not informative: 0.67 ± 1.09, −0.27 ± 0.85, against the line's −2.12 ± 1.27,
    −1.21 ± 1.09 (z = +1.7, +0.7).

  **Reading.** The shortfall at τ = 10⁻⁴ is real, not noise. A finite-τ correction is favoured: the ratio falls with τ, and its limit
  τ → 0 is consistent with 1, as it must be, since the closed form is that limit. A constant shortfall would leave that limit 3.2 and
  2.1 errors below 1. The slope is not established, though: its z is −2.3, −1.9, and the constant is not rejected by its own χ²
  (p = 0.07, 0.14). Pinning the slope needs several times Part N's seeds at τ = 0.0003, where the per-seed spread is large, a
  lower-variance estimator, or the loop itself at a smaller τ.
- **Part K's correction.** Part K's (2) model (the response to a forcing held constant from the start) over-corrects by 0.09 ± 0.02,
  0.08 ± 0.03 of the path. That is small, but it is not zero.
- **The masses form.** The 'nonlinear only + barrier' variant was not redone with the exact limit map. There the restoring drift is
  evaluated at the measured masses and cancels the offset, so the question does not arise in the same form.

## 5. Changes Paper 14 needs (not made here)

- **S2.4** (supplement §`sec:more-noise`, the Part K passage).
  - Replace the first-order offset (1) by the exact limit map, or by (I − M)⁻¹(1), with M as above. Give its size: per record
    |W − Y − (1)| up to 8.6·10⁻⁴, against |W − Y − (I − M)⁻¹(1)| at most 1.4·10⁻⁵.
  - Replace "shows only part of Katzenberger's O(γ²τ) drift path … not confirmed at this τ (an open item …)" by Part M's result.
    Y − 1/3 − (2) is 0.70, 0.63 of the path. With the loop's own step effect, it is 0.80 ± 0.09, 0.71 ± 0.06.
  - Add the paired runs at γ = 0.01, 0.003. The Richardson rate is 0.79 ± 0.08, 0.81 ± 0.09 and 0.79 ± 0.09, 0.86 ± 0.08 of
    Katzenberger's.
  - Add the fixed-mass correlation at three τ, with Part N's pooled points and line [after Part N]. The shortfall at τ = 10⁻⁴ is real
    (0.82 ± 0.06, 0.86 ± 0.06 pooled). A finite-τ origin is favoured but not established: the slope is −0.33 ± 0.14, −0.23 ± 0.12 per
    10⁻⁴ of τ, and the limit τ → 0 is 1.15 ± 0.15, 1.10 ± 0.12.
  - Add the decomposition (e).
- **Table S5** (`tab:noiseE`). The 'centres, widths' row's |z| compares the masses W with the SDE of the fit's masses. Add a note
  or row that Katzenberger's SDE describes Y. With the exact limit map, Y − 1/3 − (2) against the path has z = +3.4, −6.6, +0.1, and the
  remainder is the loop's own shortfall (Part M).
- **§9.5** (the full loop paragraph). Replace "in the first, Katzenberger's O(γ²τ) drift is largely absent from the loop" by a statement
  in three parts:
  - the first-order offset misses the fit's motion with the masses;
  - with the exact limit map and the step's effect from paired runs, the loop carries the drift at 0.79 ± 0.08, 0.81 ± 0.09 of
    Katzenberger's at γ = 0.01;
  - the remainder matches the fixed-mass correlation's at τ = 10⁻⁴, which falls with τ (Part N; the slope at z = −2.3, −1.9).
- **Prop. 7.3(a) and Table 1** (row "Prop. 7.3: … the drift's order derived, its size numerical"). Report that the O(γ²τ) drift is seen
  in the loop. Give the Richardson rate −0.00633 ± 0.00065, +0.00538 ± 0.00059 per γ²τ against −0.00803, +0.00662 at γ = 0.01, and the
  fixed-mass antisymmetric correlation at the closed form's size. The antisymmetric share is at most 0.0341, so the law is nearly
  reversible but not enough to suppress the drift.
- **Optionally, §9.5 or S2.4.** At Part E's step Δt the Euler chain's effect on the masses' rate, −0.00730 ± 0.00009, +0.01071 ± 0.00010 per
  γ²τ, is as large as the drift itself. Any measurement of the γ²τ drift needs Richardson's extrapolation in the step.
