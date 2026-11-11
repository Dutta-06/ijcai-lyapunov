# Tangential excitation and label-free parameter convergence (m=3, n=2, permuted start, T = 200 s)

| reference | seed | λ_min(Gramian, first 60 s) | min over 60 s windows | λ_max (first window) | W(0) | W(T)/W(0) | fitted rate of W, t ≥ T/2 (1/s) | labelled(T)/labelled(0) | max ‖e‖ |
|---|---|---|---|---|---|---|---|---|---|
| rich (6 tones) | 0 | 8.60e-03 | 4.12e-03 | 119.1 | 0.0159 | 0.544 | 2.31e-03 | 0.992 | 0.0019 |
| rich (6 tones) | 1 | 5.52e-03 | 1.02e-03 | 72.1 | 0.0250 | 0.409 | 2.32e-03 | 0.962 | 0.0025 |
| rich (6 tones) | 2 | 1.08e-03 | 5.91e-04 | 91.1 | 0.0238 | 0.369 | 3.07e-03 | 1.007 | 0.0022 |
| single sine | 0 | 2.54e-02 | 2.41e-02 | 42.1 | 0.0159 | 0.274 | 4.05e-03 | 0.988 | 0.0015 |
| single sine | 1 | 2.96e-02 | 2.71e-02 | 17.1 | 0.0250 | 0.314 | 4.26e-03 | 0.967 | 0.0023 |
| single sine | 2 | 6.74e-03 | 6.61e-03 | 23.3 | 0.0238 | 0.226 | 2.91e-03 | 1.008 | 0.0036 |

The matched (OT) error decays in every run; the identity-labelled error stays O(1) because the kernels converge to a relabelling of the ideal configuration.
Rank correlation between the uniform excitation level (min over windows) and the fitted rate over the 6 runs: 0.60.  The excitation level is a hypothesis of Theorem E, not a rate predictor (Remark E.1).
