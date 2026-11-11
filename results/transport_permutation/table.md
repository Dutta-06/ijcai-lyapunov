# Labelled vs optimal-transport analysis (permuted-ideal initialisation, mean-field gains)

Means over 3 seeds; 'excess' = |(1/m)Σ Taylor remainders| actually entering V̇; bounds evaluated along the same trajectory.

| m | labelled err (t=0) | OT err W (t=0) | RMS ‖e‖ | mean actual excess | mean OT bound (Thm B) | same bound, identity labelling | ratio labelled/OT | global constant per unit ‖e‖ (Cor. B.3, min-form) | quadratic branch alone at diameter |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 2.675 | 0.1070 | 0.0113 | 1.96e-05 | 8.27e-02 | 2.93e+00 | 35× | 533.2 | 38181 |
| 16 | 2.532 | 0.1044 | 0.0065 | 9.58e-06 | 4.37e-02 | 1.86e+00 | 43× | 533.2 | 38181 |
| 32 | 2.684 | 0.1122 | 0.0079 | 2.69e-05 | 7.07e-02 | 2.11e+00 | 30× | 533.2 | 38181 |
| 64 | 2.764 | 0.1136 | 0.0094 | 1.23e-05 | 6.80e-02 | 2.43e+00 | 36× | 533.2 | 38181 |

OT bound ≥ actual excess at every sample of every run: **True**
Reading: the labelled error never decreases (the kernels converge to a relabelling of the ideal), the OT error is ~1/25 of it, and the identical bound is 30–45× smaller when evaluated with the matched error (the quadratic branch is active for the matched error, the linear branch for the labelled one). Both bounds are still far above the actual excess (the Hessian supremum is global), and the global constant is set by the linear branch, not the quadratic one: the gain of the OT analysis is local, not global. No quantity grows with m.
