# Proposition D‴.51: the transit at a finite mass (`exp_transit_mass.py`)

The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width 0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, masses ε; the inner problem, the overshoot and the predictions from `rich_transit`.

## A. The near-tie law

One bump between the kernels at 0.0 and 0.1, δ from their midpoint; μ* the runner-up's margin when the entrant enters, x = μ*/ε^{1/2}; the runner-up's share |w_r|/(|w_e| + |w_r|) of the bump in the loop (at 1.32 τ_entry) against the ε-free inner problem's at the same x (its limit ε → 0 at fixed x: the bump at the midpoint, the runner-up x behind; the loop's ratio to it → 1 as δ ∝ ε^{1/2}x → 0), the Gibbs split 1/(1 + e^{κx}) (κ = πρ(2c)^{1/2} = 0.7025) and the frozen coefficient K'³ρ⁻³/(|A|(x − Δ)³) (K'³ = 0.0962; Δ from B):

| A | δ | ε | x | share: loop | inner problem | loop/inner | Gibbs split | frozen | the runner-up's centre |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.0001 | 1e-05 | 0.1098 | 0.482 | 0.48 | 1.0037 | 0.481 | — | 0.0501 |
| 0.3 | 0.0001 | 1e-06 | 0.3472 | 0.439 | 0.437 | 1.0042 | 0.439 | — | 0.0501 |
| 0.3 | 0.001 | 1e-05 | 1.096 | 0.324 | 0.309 | 1.0487 | 0.316 | — | 0.0510 |
| 0.3 | 0.0001 | 1e-07 | 1.098 | 0.311 | 0.309 | 1.0051 | 0.316 | — | 0.0501 |
| 0.3 | 0.001 | 1e-06 | 3.466 | 0.0781 | 0.0733 | 1.0663 | 0.0805 | — | 0.0510 |
| 0.3 | 0.0001 | 1e-08 | 3.472 | 0.0735 | 0.073 | 1.0063 | 0.0803 | — | 0.0501 |
| 0.3 | 0.01 | 1e-05 | 10.79 | 0.00138 | 9.99e-4 | 1.3858 | 5.10e-4 | 9.09e-4 | 0.0589 |
| 0.3 | 0.001 | 1e-07 | 10.96 | 9.65e-4 | 9.33e-4 | 1.0341 | 4.53e-4 | 8.47e-4 | 0.0487 |
| 0.3 | 0.0001 | 1e-09 | 10.98 | 9.29e-4 | 9.26e-4 | 1.0028 | 4.47e-4 | 8.41e-4 | 0.0477 |
| 0.3 | 0.01 | 1e-06 | 34.13 | 1.23e-5 | 1.16e-5 | 1.0550 | 3.88e-11 | 1.14e-5 | 0.0089 |
| 0.3 | 0.001 | 1e-08 | 34.66 | 1.11e-5 | 1.10e-5 | 1.0057 | 2.66e-11 | 1.08e-5 | 0.0071 |
| 0.3 | 0.0001 | 1e-10 | 34.72 | 1.10e-5 | 1.10e-5 | 1.0003 | 2.56e-11 | 1.08e-5 | 0.0069 |
| 0.3 | 0.01 | 1e-07 | 107.9 | 2.87e-7 | 2.84e-7 | 1.0120 | 1.19e-33 | 2.83e-7 | 0.0008 |
| 0.3 | 0.001 | 1e-09 | 109.6 | 2.71e-7 | 2.71e-7 | 1.0009 | 3.63e-34 | 2.70e-7 | 0.0006 |
| 0.3 | 0.01 | 1e-08 | 341.3 | 8.37e-9 | 8.34e-9 | 1.0033 | 7.69e-105 | 8.33e-9 | 0.0001 |
| 0.3 | 0.001 | 1e-10 | 346.6 | 7.95e-9 | 7.95e-9 | 1.0001 | 1.79e-106 | 7.95e-9 | 0.0001 |
| 0.3 | 0.01 | 1e-09 | 1079 | 2.58e-10 | 2.58e-10 | 1.0009 | 0 | 2.58e-10 | 0.0000 |
| 0.3 | 0.01 | 1e-10 | 3413 | 8.09e-12 | 8.09e-12 | 1.0003 | 0 | 8.09e-12 | 0.0000 |
| 0.8 | 0.0001 | 1e-05 | 0.1098 | 0.482 | 0.48 | 1.0037 | 0.481 | — | 0.0501 |
| 0.8 | 0.0001 | 1e-06 | 0.3472 | 0.439 | 0.437 | 1.0042 | 0.439 | — | 0.0501 |
| 0.8 | 0.001 | 1e-05 | 1.096 | 0.324 | 0.309 | 1.0487 | 0.316 | — | 0.0510 |
| 0.8 | 0.0001 | 1e-07 | 1.098 | 0.31 | 0.309 | 1.0051 | 0.316 | — | 0.0501 |
| 0.8 | 0.001 | 1e-06 | 3.466 | 0.0779 | 0.073 | 1.0667 | 0.0805 | — | 0.0510 |
| 0.8 | 0.0001 | 1e-08 | 3.472 | 0.0732 | 0.0727 | 1.0063 | 0.0803 | — | 0.0501 |
| 0.8 | 0.01 | 1e-05 | 10.79 | 9.89e-4 | 6.77e-4 | 1.4611 | 5.10e-4 | 6.40e-4 | 0.0599 |
| 0.8 | 0.001 | 1e-07 | 10.96 | 6.50e-4 | 6.25e-4 | 1.0396 | 4.53e-4 | 5.86e-4 | 0.0505 |
| 0.8 | 0.0001 | 1e-09 | 10.98 | 6.22e-4 | 6.20e-4 | 1.0033 | 4.47e-4 | 5.81e-4 | 0.0496 |
| 0.8 | 0.01 | 1e-06 | 34.13 | 5.29e-6 | 5.01e-6 | 1.0553 | 3.88e-11 | 4.90e-6 | 0.0097 |
| 0.8 | 0.001 | 1e-08 | 34.66 | 4.76e-6 | 4.74e-6 | 1.0055 | 2.66e-11 | 4.64e-6 | 0.0077 |
| 0.8 | 0.0001 | 1e-10 | 34.72 | 4.71e-6 | 4.71e-6 | 1.0000 | 2.56e-11 | 4.61e-6 | 0.0075 |
| 0.8 | 0.01 | 1e-07 | 107.9 | 1.12e-7 | 1.11e-7 | 1.0114 | 1.19e-33 | 1.10e-7 | 0.0008 |
| 0.8 | 0.001 | 1e-09 | 109.6 | 1.06e-7 | 1.05e-7 | 1.0009 | 3.63e-34 | 1.05e-7 | 0.0007 |
| 0.8 | 0.01 | 1e-08 | 341.3 | 3.17e-9 | 3.16e-9 | 1.0030 | 7.69e-105 | 3.16e-9 | 0.0001 |
| 0.8 | 0.001 | 1e-10 | 346.6 | 3.02e-9 | 3.02e-9 | 1.0001 | 1.79e-106 | 3.02e-9 | 0.0001 |
| 0.8 | 0.01 | 1e-09 | 1079 | 9.71e-11 | 9.70e-11 | 1.0009 | 0 | 9.70e-11 | 0.0000 |
| 0.8 | 0.01 | 1e-10 | 3413 | 3.04e-12 | 3.04e-12 | 1.0003 | 0 | 3.04e-12 | 0.0000 |
| 2 | 0.0001 | 1e-05 | 0.1098 | 0.482 | 0.48 | 1.0037 | 0.481 | — | 0.0501 |
| 2 | 0.0001 | 1e-06 | 0.3472 | 0.438 | 0.437 | 1.0042 | 0.439 | — | 0.0501 |
| 2 | 0.001 | 1e-05 | 1.096 | 0.324 | 0.309 | 1.0487 | 0.316 | — | 0.0510 |
| 2 | 0.0001 | 1e-07 | 1.098 | 0.31 | 0.309 | 1.0051 | 0.316 | — | 0.0501 |
| 2 | 0.001 | 1e-06 | 3.466 | 0.0778 | 0.0729 | 1.0668 | 0.0805 | — | 0.0510 |
| 2 | 0.0001 | 1e-08 | 3.472 | 0.0731 | 0.0727 | 1.0064 | 0.0803 | — | 0.0501 |
| 2 | 0.01 | 1e-05 | 10.79 | 7.93e-4 | 5.13e-4 | 1.5448 | 5.10e-4 | 5.35e-4 | 0.0600 |
| 2 | 0.001 | 1e-07 | 10.96 | 4.89e-4 | 4.68e-4 | 1.0457 | 4.53e-4 | 4.79e-4 | 0.0510 |
| 2 | 0.0001 | 1e-09 | 10.98 | 4.66e-4 | 4.64e-4 | 1.0038 | 4.47e-4 | 4.74e-4 | 0.0501 |
| 2 | 0.01 | 1e-06 | 34.13 | 2.42e-6 | 2.29e-6 | 1.0560 | 3.88e-11 | 2.23e-6 | 0.0105 |
| 2 | 0.001 | 1e-08 | 34.66 | 2.17e-6 | 2.16e-6 | 1.0056 | 2.66e-11 | 2.11e-6 | 0.0084 |
| 2 | 0.0001 | 1e-10 | 34.72 | 2.15e-6 | 2.15e-6 | 1.0000 | 2.56e-11 | 2.10e-6 | 0.0082 |
| 2 | 0.01 | 1e-07 | 107.9 | 4.64e-8 | 4.59e-8 | 1.0109 | 1.19e-33 | 4.58e-8 | 0.0008 |
| 2 | 0.001 | 1e-09 | 109.6 | 4.38e-8 | 4.37e-8 | 1.0008 | 3.63e-34 | 4.36e-8 | 0.0007 |
| 2 | 0.01 | 1e-08 | 341.3 | 1.28e-9 | 1.28e-9 | 1.0028 | 7.69e-105 | 1.28e-9 | 0.0001 |
| 2 | 0.001 | 1e-10 | 346.6 | 1.22e-9 | 1.22e-9 | 1.0001 | 1.79e-106 | 1.22e-9 | 0.0001 |
| 2 | 0.01 | 1e-09 | 1079 | 3.90e-11 | 3.89e-11 | 1.0008 | 0 | 3.89e-11 | 0.0000 |
| 2 | 0.01 | 1e-10 | 3413 | 1.22e-12 | 1.22e-12 | 1.0002 | 0 | 1.22e-12 | 0.0000 |

The loop's share over the inner problem's: 1.0000–1.0064 at δ = 0.0001; 1.0001–1.0668 at δ = 0.001; 1.0002–1.5448 at δ = 0.01 (x from 0.11 to 3413; the inner problem stationary — every rate below 10⁻⁹|A| — on 51 of 54; the other 3 (x ≈ 34.1–34.7) still drifting, their share changing by at most 9.8e-04 (relative) from s = 4.69e+06 to 1.88e+07).

## B. The entrant's overshoot and the runner-up's pull

The overshoot η_c(α_A)/(Λ₀ρ) (an entered kernel's dual beyond its wall, in units of ε^{1/2}ρ) against the pull Δ the inner problem exerts on a frozen runner-up (x = 200; matched: the runner-up's outer margin at the entrant's arrival less its final ξ) at the exact tie and at δ = 0.001, from the starts xi0 = 400, 800, 1600, 3200 behind the wall (the entrant's far-field creep: the pull converges like a/xi0 + b/xi0²; ∞ the second-order Richardson extrapolation from the last three), and the loop's (δ = 0.001, ε = 1e-08, 1e-09, 1e-10: Δ = x − ξ_r(∞), ξ_r read from the runner-up's coefficient), free and with the runner-up's centre held fixed: the difference is the runner-up's own creep while it waits at its margin, (a²b/(π²c))(g′_r/g_r)²/x in units of ε^{1/2} at leading order (predicted; their ratio in parentheses), and the fixed loop's excess over the inner problem at that offset the rest; I∞ = overshoot − Δ (tie, ∞):

| A | overshoot | Δ: inner, tie (xi0 = 400, 800, 1600, 3200, ∞) | I∞ | Δ: inner, δ (∞) | Δ: loop (x) | Δ: loop, the runner-up's centre fixed | free − fixed: predicted (ratio) | fixed − inner, δ |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 3.6500 | 2.2888, 2.2599, 2.2451, 2.2376, 2.2300 | 1.4200 | 2.2776 | 2.6552 (34.66), 2.3953 (109.6), 2.3134 (346.6) | 2.2800, 2.2783, 2.2768 | 0.3753: 0.3667 (1.023), 0.1170: 0.1159 (1.009), 0.0366: 0.0367 (0.998) | +0.0024, +0.0008, -0.0007 |
| 0.2 | 4.6370 | 3.2337, 3.2048, 3.1900, 3.1825, 3.1749 | 1.4621 | 3.2208 | 3.6091 (34.66), 3.3396 (109.6), 3.2568 (346.6) | 3.2234, 3.2216, 3.2201 | 0.3857: 0.3667 (1.052), 0.1181: 0.1159 (1.018), 0.0367: 0.0367 (1.001) | +0.0026, +0.0008, -0.0007 |
| 0.3 | 5.2143 | 3.7864, 3.7575, 3.7427, 3.7351, 3.7275 | 1.4867 | 3.7725 | 4.1673 (34.66), 3.8920 (109.6), 3.8085 (346.6) | 3.7752, 3.7733, 3.7718 | 0.3921: 0.3667 (1.069), 0.1187: 0.1159 (1.024), 0.0367: 0.0367 (1.002) | +0.0027, +0.0008, -0.0007 |
| 0.5 | 5.9415 | 4.4827, 4.4538, 4.4389, 4.4314, 4.4238 | 1.5177 | 4.4676 | 4.8707 (34.66), 4.5878 (109.6), 4.5036 (346.6) | 4.4703, 4.4683, 4.4668 | 0.4004: 0.3667 (1.092), 0.1195: 0.1159 (1.030), 0.0368: 0.0367 (1.004) | +0.0027, +0.0008, -0.0007 |
| 0.8 | 6.6105 | 5.1233, 5.0944, 5.0796, 5.0720, 5.0644 | 1.5461 | 5.1071 | 5.5182 (34.66), 5.2280 (109.6), 5.1432 (346.6) | 5.1098, 5.1078, 5.1063 | 0.4084: 0.3667 (1.114), 0.1202: 0.1159 (1.037), 0.0369: 0.0367 (1.006) | +0.0028, +0.0008, -0.0007 |
| 1.2 | 7.1877 | 5.6759, 5.6471, 5.6322, 5.6247, 5.6171 | 1.5706 | 5.6588 | 6.0771 (34.66), 5.7804 (109.6), 5.6950 (346.6) | 5.6616, 5.6595, 5.6580 | 0.4155: 0.3667 (1.133), 0.1208: 0.1159 (1.042), 0.0369: 0.0367 (1.008) | +0.0028, +0.0008, -0.0007 |
| 2 | 7.9149 | 6.3722, 6.3433, 6.3285, 6.3210, 6.3134 | 1.6015 | 6.3538 | 6.7815 (34.66), 6.4762 (109.6), 6.3901 (346.6) | 6.3566, 6.3546, 6.3531 | 0.4248: 0.3667 (1.159), 0.1216: 0.1159 (1.049), 0.0370: 0.0367 (1.010) | +0.0028, +0.0008, -0.0007 |
| 4 | 8.9016 | 7.3170, 7.2881, 7.2733, 7.2658, 7.2582 | 1.6435 | 7.2970 | 7.7379 (34.66), 7.4205 (109.6), 7.3333 (346.6) | 7.2998, 7.2977, 7.2962 | 0.4381: 0.3667 (1.195), 0.1228: 0.1159 (1.059), 0.0371: 0.0367 (1.013) | +0.0028, +0.0008, -0.0007 |

The overshoot's growth per unit log A between consecutive amplitudes: 1.4240, 1.4237, 1.4236, 1.4235, 1.4235, 1.4235, 1.4235; the exponential branch's slope a/(4GcΛ₀ρ) = 1.4235 (Λ₀ = 7.0248). The inner profile against its core power law, 𝒲(ξ)/(K'³ξ⁻³) − 1 at ξ = 15, 10, 8.3, 5, 3: -6.7%, -14.0%, -19.1%, -39.5%, -65.2%.

## C. The creep

Target 10 of D‴.48's random targets (centres [-0.752, 0.112, 0.488], amplitudes [0.868, -0.412, -0.514]; the path's kernels [2, 14, 12]). Predicted from the path: the entry threshold 6.12e-7 (kernel 3 at the entry of 2, margin 0.00489), the creep thresholds 15: 4.88e-9, 3: 1.36e-6, 11: 5.56e-5, 13: 2.62e-6: ε_c = 4.88e-9. The loop:

| ε | ε/ε_c | formed | entrants | the path's kernels | the inactive kernels' largest displacement |
|---|---|---|---|---|---|
| 1e-09 | 0.205 | 3 | [2, 14, 12] | yes | 0.0884 |
| 2e-09 | 0.41 | 3 | [2, 14, 12] | yes | 0.0884 |
| 3.5e-09 | 0.717 | 3 | [2, 14, 12] | yes | 0.138 |
| 7e-09 | 1.43 | 3 | [2, 14, 15] | no | 0.288 |
| 1.4e-08 | 2.87 | 3 | [2, 14, 15] | no | 0.288 |
| 3e-08 | 6.14 | 3 | [2, 14, 15] | no | 0.288 |

## D. Prediction from the path (out of sample)

48 random three-bump targets (seed 51); ε_c = min(ε_entry, ε_creep) predicted from the sliding path; the loop recovers when it forms exactly the path's kernels. Predicted (ε < ε_c) against observed:

| ε | recovered | predicted to recover | agree | recovered though not predicted | not recovered though predicted |
|---|---|---|---|---|---|
| 1e-05 | 32 of 48 | 27 | 43 | 5 | 0 |
| 1e-06 | 38 of 48 | 36 | 46 | 2 | 0 |
| 1e-07 | 42 of 48 | 43 | 47 | 0 | 1 |
| 1e-08 | 47 of 48 | 47 | 48 | 0 | 0 |

The disagreements (target, ε, ε_c): (15, 1e-05, 9.22e-07), (27, 1e-05, 6.18e-07), (28, 1e-05, 2.45e-06), (45, 1e-05, 5.73e-06), (46, 1e-05, 1.85e-06), (15, 1e-06, 9.22e-07), (27, 1e-06, 6.18e-07), (20, 1e-07, 1.17e-07); the path itself exact (three kernels, rms residual < 10⁻⁸) on 48 of 48.

## E. The predicted mass over many targets

400 random three-bump targets (seed 52): the fraction predicted not to recover at ε (ε ≥ ε_c), and the fractions with ε ≥ ε_entry (an entry's near tie alone would fail them) and with ε ≥ ε_creep (the creep alone) — not exclusive; the local slope of each in log ε. These are the predictor's fractions; the loop's own were not measured:

| ε | not recovered | by a near tie | by the creep | slope: all | near tie | creep |
|---|---|---|---|---|---|---|
| 0.0001 | 0.7700 | 0.6950 | 0.7025 | — | — | — |
| 3.16e-05 | 0.5975 | 0.4350 | 0.5375 | 0.22 | 0.407 | 0.233 |
| 1e-05 | 0.4650 | 0.2775 | 0.4250 | 0.218 | 0.39 | 0.204 |
| 3.16e-06 | 0.3225 | 0.1550 | 0.2925 | 0.318 | 0.506 | 0.325 |
| 1e-06 | 0.2400 | 0.1050 | 0.2125 | 0.257 | 0.338 | 0.278 |
| 3.16e-07 | 0.1650 | 0.0550 | 0.1525 | 0.325 | 0.562 | 0.288 |
| 1e-07 | 0.1225 | 0.0300 | 0.1175 | 0.259 | 0.526 | 0.226 |
| 3.16e-08 | 0.0900 | 0.0150 | 0.0900 | 0.268 | 0.602 | 0.232 |
| 1e-08 | 0.0700 | 0.0125 | 0.0650 | 0.218 | 0.158 | 0.283 |
| 3.16e-09 | 0.0400 | 0.0050 | 0.0375 | 0.486 | 0.796 | 0.478 |
| 1e-09 | 0.0300 | 0.0025 | 0.0300 | 0.25 | 0.602 | 0.194 |
| 3.16e-10 | 0.0225 | 0.0025 | 0.0225 | 0.25 | 0 | 0.25 |
| 1e-10 | 0.0150 | 0.0025 | 0.0150 | 0.352 | 0 | 0.352 |
| 3.16e-11 | 0.0100 | 0.0025 | 0.0100 | 0.352 | 0 | 0.352 |
| 1e-11 | 0.0075 | 0.0025 | 0.0050 | 0.25 | 0 | 0.602 |
| 3.16e-12 | 0.0025 | 0.0000 | 0.0025 | 0.954 | — | 0.602 |
| 1e-12 | 0.0025 | 0.0000 | 0.0025 | 0 | — | 0 |

Least-squares slopes in log ε where between 10 targets and half of them fail: all 0.294 (ε = 1e-05 … 1e-09); a near tie 0.462 (ε = 3.16e-05 … 1e-07); the creep 0.285 (ε = 1e-05 … 1e-09).

## F. The predicted tail against the loop

Every target of E with 1e-12 ≤ ε_c ≤ 1e-07 (48 targets; 1 more below 1e-12, not run), the loop at 3ε_c (predicted not to recover) and at ε_c/3 (predicted to recover):

| | runs | as predicted | binding: entry | binding: creep |
|---|---|---|---|---|
| 3ε_c | 48 | 43 | 7 | 41 |
| ε_c/3 | 48 | 48 | 7 | 41 |

Not as predicted (target, ε/ε_c, the binding threshold, formed against the path's): (8, 3, creep, [15, 11, 6] vs [15, 11, 6]); (115, 3, creep, [14, 5, 4] vs [14, 5, 4]); (117, 3, creep, [12, 13, 2] vs [12, 2, 13]); (216, 3, creep, [3, 11, 10] vs [3, 11, 10]); (267, 3, creep, [14, 15, 2] vs [14, 15, 2]).

## G. The loop's own threshold against the prediction

D's 48 targets: the largest mass at which the loop still forms exactly the path's kernels, bracketed from ε_c by factors 4 and bisected 4 times in log ε (the bracket ×1.09); ε_loop its geometric mid-point.

| binding | targets | ε_loop/ε_c: median [quartiles] | range | below 1 (the loop fails below ε_c) | within ×1.5 | within ×2 |
|---|---|---|---|---|---|---|
| entry | 32 | 0.648 [0.545, 0.808] | 0.131–43.3 | 30 | 14 | 27 |
| creep | 16 | 0.844 [0.545, 8.24] | 0.324–94.5 | 11 | 6 | 9 |
| all | 48 | 0.771 [0.545, 0.917] | 0.131–94.5 | 41 | 20 | 36 |

Entry-bound targets by the predicted mass (∝ the binding margin squared): ε_c in (0, 2e-06]: 5 targets, median 0.917, range 0.917–43.3; ε_c in (2e-06, 0.0001]: 17 targets, median 0.771, range 0.131–1.09; ε_c in (0.0001, inf]: 10 targets, median 0.545, range 0.273–0.648.

Per target (ε_c, the binding threshold, ε_loop/ε_c): 0: 3.03e-05 (entry) 0.771; 1: 2.14e-05 (entry) 0.917; 2: 8.85e-06 (entry) 0.131; 3: 0.00015 (entry) 0.545; 4: 0.000122 (entry) 0.648; 5: 1.19e-05 (creep) 5.42; 6: 0.000197 (entry) 0.545; 7: 0.000174 (entry) 0.273; 8: 2.47e-05 (entry) 0.771; 9: 5.52e-07 (entry) 0.917; 10: 0.000108 (entry) 0.648; 11: 9.2e-06 (creep) 0.386; 12: 1.9e-05 (entry) 0.771; 13: 2.98e-07 (creep) 0.771; 14: 2.24e-05 (entry) 0.771; 15: 9.22e-07 (creep) 94.5; 16: 5.63e-05 (entry) 0.648; 17: 2.06e-10 (creep) 0.545; 18: 5.07e-05 (entry) 0.648; 19: 0.000262 (entry) 0.545; 20: 1.17e-07 (creep) 0.771; 21: 1.94e-06 (entry) 0.917; 22: 4.38e-06 (entry) 1.09; 23: 1.91e-06 (creep) 0.545; 24: 0.000639 (creep) 0.324; 25: 8.93e-08 (creep) 0.771; 26: 4.74e-05 (entry) 0.545; 27: 6.18e-07 (entry) 43.3; 28: 2.45e-06 (creep) 25.8; 29: 1.57e-05 (entry) 0.771; 30: 2.67e-07 (creep) 0.917; 31: 0.000173 (entry) 0.545; 32: 1.32e-07 (entry) 0.917; 33: 0.000176 (entry) 0.545; 34: 7.99e-05 (entry) 0.648; 35: 1.67e-08 (creep) 0.917; 36: 9.16e-05 (entry) 0.917; 37: 6.34e-05 (entry) 0.648; 38: 4.85e-08 (entry) 0.917; 39: 0.000795 (entry) 0.386; 40: 5.9e-05 (entry) 0.459; 41: 0.000225 (entry) 0.545; 42: 8.23e-08 (creep) 0.917; 43: 1.38e-05 (entry) 0.771; 44: 2.77e-06 (creep) 0.545; 45: 5.73e-06 (creep) 16.7; 46: 1.85e-06 (creep) 19.9; 47: 4.39e-05 (entry) 0.771.

F's targets not as predicted, their loop thresholds bracketed and bisected as in G:

| target (E) | binding | ε_c | ε_loop | ε_loop/ε_c |
|---|---|---|---|---|
| 8 | creep | 3.48e-10 | 3.17e-09 | 9.11 |
| 115 | creep | 3.85e-10 | 4.43e-05 | 1.15e+05 |
| 117 | creep | 6.68e-08 | 0.000156 | 2.33e+03 |
| 216 | creep | 1.04e-08 | 1.04e-07 | 9.93 |
| 267 | creep | 6.21e-09 | 1.4e-06 | 225 |

## H. Near-tie controls

(1) The single bump (A = 0.8) at δ from the midpoint of the kernels at 0.0 and 0.1: the runner-up's formation threshold, x_θ = μ*/ε^{1/2} at the loop's bisected threshold (12 bisections), free and with centres held fixed until the entrant's arrival τ_e (the inner problem's x_θ and the Gibbs value in A):

| δ | μ* | free | the runner-up's centre fixed before τ_e | both centres fixed before τ_e |
|---|---|---|---|---|
| 0.001 | 0.0035 | 6.392 | 5.315 | 6.557 |
| 0.01 | 0.0341 | 7.100 | 5.819 | 6.546 |
| 0.03 | 0.0989 | 8.189 | 6.373 | 6.520 |

The inner problem's own threshold at the exact tie (the limit ε → 0 at fixed x), bisected in x: A = 0.3: 6.404, A = 0.8: 6.300, A = 2: 6.272 (the Gibbs value ln(|A|/θ − 1)/κ = 6.541). The free loop's shift from it against the partners' approach-creep difference at the actual geometry, (a²b/(π²c))((g′_r/g_r)² − (g′_e/g_e)²)/x_θ (g = ⟨fφ⟩, g′ = ⟨f∂_cφ⟩ before the entry):

| δ | μ* | x_θ (free) − inner | creep difference / x_θ | ratio |
|---|---|---|---|---|
| 0.001 | 0.0035 | 0.091 | 0.153 | 0.597 |
| 0.01 | 0.0341 | 0.800 | 1.376 | 0.581 |
| 0.03 | 0.0989 | 1.889 | 3.580 | 0.528 |

The predictor's own threshold for the runner-up (`recovery_threshold`): x = μ*/ε_entry^{1/2} = max(Δ + ξ_θ, x_θκ_ij) — the frozen runner-up's pull plus the formation level's inner coordinate, and the near-tie estimate: δ = 0.001: 6.252 (kernel 10; Δ + ξ_θ = 5.146, x_θκ_ij = 6.541·0.956 = 6.252); δ = 0.01: 6.145 (kernel 10; Δ + ξ_θ = 5.443, x_θκ_ij = 6.541·0.939 = 6.145); δ = 0.03: 5.853 (kernel 10; Δ + ξ_θ = 5.722, x_θκ_ij = 6.541·0.895 = 5.853).

(2) Exact ties made by unequal masses (the farther kernel's mass b set so that both reach their walls together):

| bump | nearer kernel | farther kernel (b) | ε | the farther kernel's share | formed | the path |
|---|---|---|---|---|---|---|
| 0.07625 | 11 | 10 (1.2000) | 1e-08 | 0.8849 | [10, 11] | [10, 11] |
| 0.07625 | 11 | 10 (1.2000) | 1e-10 | 0.8854 | [10, 11] | [10, 11] |
| 0.08 | 11 | 12 (1.6260) | 1e-08 | 0.9972 | [12] | [11, 12] |
| 0.08 | 11 | 12 (1.6260) | 1e-10 | 0.9973 | [12] | [11, 12] |

(3) D's target 2 with its middle bump moved to c₂: the first entry's runner-up and its margin μ*, the predicted ε_c (entry, creep), the loop's threshold (bracketed and bisected as in G), and the kernels it forms just above it (an extra kernel's final centre):

| c₂ | runner-up (μ*) | ε_c: entry, creep | ε_loop | ε_loop/ε_c | ε_loop/μ*² | ε_loop/μ*³ | formed just above (extra: final centre) |
|---|---|---|---|---|---|---|---|
| 0.03 | 11 (0.06205) | 1.08e-05, 1.97e-05 | 8.33e-06 | 0.771 | 0.00216 | 0.0349 | [10, 16, 3, 4] (4: -0.655) |
| 0.035 | 11 (0.04568) | 1.07e-05, 1.95e-05 | 8.25e-06 | 0.771 | 0.00395 | 0.0866 | [10, 16, 3, 4] (4: -0.655) |
| 0.04 | 11 (0.02903) | 1.06e-05, 1.93e-05 | 8.18e-06 | 0.771 | 0.00971 | 0.334 | [10, 16, 3, 4] (4: -0.655) |
| 0.043 | 11 (0.01890) | 9.3e-06, 1.26e-05 | 1.32e-06 | 0.142 | 0.00371 | 0.196 | [10, 11, 16, 3] (11: 0.632) |
| 0.046 | 11 (0.00866) | 1.93e-06, 1.75e-07 | 9.55e-08 | 0.545 | 0.00127 | 0.147 | [10, 11, 3] (11: 0.632) |
| 0.048 | 11 (0.00177) | 8.03e-08, 4.24e-10 | 3.27e-10 | 0.771 | 0.000104 | 0.0587 | [10, 11, 3] (11: 0.632) |

Where the extra kernel is the first entry's runner-up, the loop's threshold against μ*: slope 3.52 in log–log.

(4) E's target 117 — kernel 14 between the third entrant (13) and its bump: its largest coefficient over the run (the formation level θ = 10⁻²max|A*|) and its closest approach to the entrant:

| ε | x = μ/ε^{1/2} | max\|w₁₄\| | θ | min \|c₁₃ − c₁₄\| | formed | the path | the Gibbs split at x |
|---|---|---|---|---|---|---|---|
| 1e-06 | 5.186 | 1e-07 | 0.00896 | 0.0053 | [12, 13, 2] | [12, 2, 13] | 0.026 |
| 1e-05 | 1.640 | 4.15e-06 | 0.00896 | 0.0716 | [12, 13, 2] | [12, 2, 13] | 0.240 |
| 0.0001 | 0.519 | 0.000653 | 0.00896 | 0.0000 | [12, 13, 2] | [12, 2, 13] | 0.410 |

Numbers are produced by `experiments/exp_transit_mass.py`; nothing is edited by hand.
