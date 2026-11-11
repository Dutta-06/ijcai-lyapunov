# OT σ-modification vs labelled σ-modification (oracle plant, m=6, 5 seeds, σ_m = 0.3)

Same plant, same prior, same geometry and same x₀ in both starts; only the particle labels differ. RMS ‖e‖ over the second half of 40 s; W_end = OT error at T; lab_end = labelled error at T; b = W(μ*, μ⁰); b_lab = labelled bias seen by the labelled leak.

| start | mechanism | RMS ‖e‖ (mean ± 95% CI) | W_end | lab_end | b | b_lab |
|---|---|---|---|---|---|---|
| aligned | projection | 0.0051 ± 0.0010 | 0.2953 | 0.2953 | 0.308 | 0.308 |
| aligned | labelled sigma-mod | 0.0091 ± 0.0018 | 0.3040 | 0.3040 | 0.308 | 0.308 |
| aligned | OT sigma-mod | 0.0091 ± 0.0018 | 0.3040 | 0.3040 | 0.308 | 0.308 |
| relabelled | projection | 0.0051 ± 0.0010 | 0.2953 | 0.2953 | 0.308 | 2.059 |
| relabelled | labelled sigma-mod | 0.0091 ± 0.0018 | 0.3040 | 2.0590 | 0.308 | 2.059 |
| relabelled | OT sigma-mod | 0.0091 ± 0.0018 | 0.3040 | 0.3040 | 0.308 | 2.059 |

## Label invariance (same seed, aligned vs relabelled start): |RMS_aligned − RMS_relabelled|

| mechanism | mean | max | relabelled worse in |
|---|---|---|---|
| projection | 4.86e-18 | 1.13e-17 | 3/5 |
| labelled sigma-mod | 4.05e-06 | 8.51e-06 | 1/5 |
| OT sigma-mod | 2.60e-18 | 5.20e-18 | 0/5 |

Reading: projection and OT σ-modification are exactly label-invariant; the labelled σ-modification is not — with a relabelled start its leak drags every kernel toward a different kernel's prior position, the bias it pays is b_lab instead of b, and tracking degrades accordingly.
