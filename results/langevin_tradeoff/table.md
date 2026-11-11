# Speed–bias trade-off of parameter dither (m = 32, σ = 0.002, k_e = 10, T = 1500 s, 2 seeds)

| τ | L̄ final, collapsed | L̄ final, spread | escape time of the collapsed start (s) | W₂(final collapsed, final spread) | tracking RMS, spread |
|---|---|---|---|---|---|
| 0.0 | 1.66e-01 | 6.96e-03 | ∞, ∞ | 2.89 | 0.0116 |
| 0.0005 | 2.49e-02 | 1.09e-02 | ∞, 1050 | 1.81 | 0.0153 |
| 0.001 | 2.16e-02 | 1.48e-02 | 1480, 1005 | 1.56 | 0.0179 |
| 0.002 | 2.88e-02 | 2.33e-02 | 1115, 800 | 1.33 | 0.0220 |
| 0.004 | 3.92e-02 | 3.62e-02 | 970, 430 | 1.08 | 0.0275 |
| 0.008 | 6.03e-02 | 6.53e-02 | 390, 295 | 1.20 | 0.0363 |
| annealed (hot start) | 3.17e-02 | 3.20e-02 | 990, 530 | 1.11 | 0.0252 |

Escape time = first sample at which the collapsed run's loss ≤ 2 × the spread run's final loss (same seed); ∞ = not within the horizon.  Per-seed values listed.

## Escape-time scaling (τ values at which every seed escapes)

| model | fitted form (seed-averaged) | RMS relative residual, seed-averaged | RMS relative residual, per seed |
|---|---|---|---|
| t_esc = A + B·log(1/τ) | -1693 + 427·log(1/τ) | 0.045 | 0.274 |
| t_esc = a·exp(b/τ) (Arrhenius) | 407·exp(1.25e-03/τ) | 0.253 | 0.330 |

Reading: the seed-averaged fits differ, but on the per-seed values both laws fit poorly, so these data do not identify the escape mechanism (PARADIGM.md Remark H.2).
