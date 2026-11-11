# Lemma E.0″: the slow-regime strict Lyapunov function on non-periodic PE regressors (`exp_pe_general.py`)

Exponential averages over a window T on [0, 300] (grid step 0.005, the averages formed from the extended grid); the certificate is the supremum of V̇/V over 400 times and all states; the exact decay rate of (LTV) is −(1/200)log‖Φ(200, 0)‖₂ (top Lyapunov exponent over 200 time units). Constants: 'measured' uses the grid suprema of ‖P₁‖, ‖P₂‖, ‖P₁v‖/‖v‖_{M̄}, ‖P₂v‖/‖v‖_{M̄}, ‖N̄₂‖ (computable constants of the regressor, like the peak factors); 'explicit' uses the closed-form bounds π₁ = 2T(√Λ + ‖M̂‖/√λ) etc.

## quasi-periodic Ψ = (sin t, cos √2 t)

Exact decay rates of (LTV): k = 20: 0.02485, k = 50: 0.00995, k = 100: 0.00497, k = 200: 0.00249, k = 400: 0.00124 (the long-time mean of M has λ_min = 0.5, so λ_min/k = 0.02500, 0.01000, 0.00500, 0.00250, 0.00125).

### window T = 6.283: λ = inf λ_min(M̄_T) = 0.2967, Λ = 0.7276, ρ = 4.00, ρ_D = 6.35, ρ_DD = 10.01, Ψ_M = 1.414, Ψ_D = 1.732, Ψ_DD = 2.236, p₁ = 2.87, p₂ = 8.00, π₁ = 4.71, π₂ = 11.63, ‖N̄₂‖ = 0.706

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | ABCD | 0.119 | 0.00673 | 0.02485 | 3.69 | 0.01483 | -0.02567, -0.01346 | ABD, nan |
| 20 | 0.1 | ABCD | 0.119 | 0.00800 | 0.02485 | 3.10 | 0.01483 | -0.02574, -0.01601 | ABD, nan |
| 50 | 0.2 | ABCD | 0.048 | 0.00408 | 0.00995 | 2.44 | 0.00593 | -0.01112, -0.00815 | ABCD, 0.00308 |
| 50 | 0.1 | ABCD | 0.048 | 0.00464 | 0.00995 | 2.15 | 0.00593 | -0.01115, -0.00927 | ABCD, 0.00364 |
| 100 | 0.2 | ABCD | 0.024 | 0.00222 | 0.00497 | 2.24 | 0.00297 | -0.00572, -0.00444 | ABCD, 0.00209 |
| 100 | 0.1 | ABCD | 0.024 | 0.00251 | 0.00497 | 1.98 | 0.00297 | -0.00574, -0.00502 | ABCD, 0.00238 |
| 200 | 0.2 | ABCD | 0.012 | 0.00115 | 0.00249 | 2.16 | 0.00148 | -0.00291, -0.00230 | ABCD, 0.00113 |
| 200 | 0.1 | ABCD | 0.012 | 0.00130 | 0.00249 | 1.92 | 0.00148 | -0.00291, -0.00259 | ABCD, 0.00128 |
| 400 | 0.2 | ABCD | 0.006 | 0.00058 | 0.00124 | 2.13 | 0.00074 | -0.00146, -0.00117 | ABCD, 0.00058 |
| 400 | 0.1 | ABCD | 0.006 | 0.00066 | 0.00124 | 1.89 | 0.00074 | -0.00147, -0.00132 | ABCD, 0.00066 |

### window T = 12.57: λ = inf λ_min(M̄_T) = 0.3907, Λ = 0.6157, ρ = 4.00, ρ_D = 6.09, ρ_DD = 10.00, Ψ_M = 1.414, Ψ_D = 1.732, Ψ_DD = 2.236, p₁ = 2.92, p₂ = 8.93, π₁ = 4.42, π₂ = 12.72, ‖N̄₂‖ = 0.351

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | ABCD | 0.045 | 0.00999 | 0.02485 | 2.49 | 0.01953 | -0.03393, -0.01998 | ABD, nan |
| 20 | 0.1 | ABCD | 0.045 | 0.01166 | 0.02485 | 2.13 | 0.01953 | -0.03402, -0.02332 | BD, nan |
| 50 | 0.2 | ABCD | 0.018 | 0.00557 | 0.00995 | 1.78 | 0.00781 | -0.01465, -0.01115 | ABCD, 0.00322 |
| 50 | 0.1 | ABCD | 0.018 | 0.00631 | 0.00995 | 1.58 | 0.00781 | -0.01470, -0.01262 | ABCD, 0.00396 |
| 100 | 0.2 | ABCD | 0.009 | 0.00298 | 0.00497 | 1.67 | 0.00391 | -0.00753, -0.00595 | ABCD, 0.00268 |
| 100 | 0.1 | ABCD | 0.009 | 0.00336 | 0.00497 | 1.48 | 0.00391 | -0.00756, -0.00671 | ABCD, 0.00306 |
| 200 | 0.2 | ABCD | 0.004 | 0.00153 | 0.00249 | 1.63 | 0.00195 | -0.00382, -0.00306 | ABCD, 0.00149 |
| 200 | 0.1 | ABCD | 0.004 | 0.00172 | 0.00249 | 1.45 | 0.00195 | -0.00384, -0.00344 | ABCD, 0.00168 |
| 400 | 0.2 | ABCD | 0.002 | 0.00077 | 0.00124 | 1.61 | 0.00098 | -0.00193, -0.00155 | ABCD, 0.00077 |
| 400 | 0.1 | ABCD | 0.002 | 0.00087 | 0.00124 | 1.43 | 0.00098 | -0.00193, -0.00174 | ABCD, 0.00087 |

### window T = 25.13: λ = inf λ_min(M̄_T) = 0.4439, Λ = 0.5577, ρ = 4.00, ρ_D = 6.02, ρ_DD = 10.00, Ψ_M = 1.414, Ψ_D = 1.732, Ψ_DD = 2.236, p₁ = 2.92, p₂ = 9.19, π₁ = 4.26, π₂ = 13.02, ‖N̄₂‖ = 0.174

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | ABCD | 0.020 | 0.01181 | 0.02485 | 2.10 | 0.02219 | -0.03874, -0.02363 | ABD, nan |
| 20 | 0.1 | ABCD | 0.020 | 0.01371 | 0.02485 | 1.81 | 0.02219 | -0.03883, -0.02743 | BD, nan |
| 50 | 0.2 | ABCD | 0.008 | 0.00642 | 0.00995 | 1.55 | 0.00888 | -0.01667, -0.01284 | ABCD, 0.00134 |
| 50 | 0.1 | ABCD | 0.008 | 0.00725 | 0.00995 | 1.37 | 0.00888 | -0.01673, -0.01451 | ABCD, 0.00217 |
| 100 | 0.2 | ABCD | 0.004 | 0.00340 | 0.00497 | 1.46 | 0.00444 | -0.00856, -0.00681 | ABCD, 0.00276 |
| 100 | 0.1 | ABCD | 0.004 | 0.00384 | 0.00497 | 1.30 | 0.00444 | -0.00860, -0.00767 | ABCD, 0.00319 |
| 200 | 0.2 | ABCD | 0.002 | 0.00174 | 0.00249 | 1.43 | 0.00222 | -0.00434, -0.00348 | ABCD, 0.00166 |
| 200 | 0.1 | ABCD | 0.002 | 0.00196 | 0.00249 | 1.27 | 0.00222 | -0.00436, -0.00392 | ABCD, 0.00188 |
| 400 | 0.2 | ABCD | 0.001 | 0.00088 | 0.00124 | 1.41 | 0.00111 | -0.00219, -0.00176 | ABCD, 0.00087 |
| 400 | 0.1 | ABCD | 0.001 | 0.00099 | 0.00124 | 1.26 | 0.00111 | -0.00220, -0.00198 | ABCD, 0.00098 |

## chaotic Ψ = (x₁/20, (x₃ − 25)/20) along a Lorenz orbit traversed 5 times slower than standard time

Exact decay rates of (LTV): k = 20: 0.00759, k = 50: 0.00304, k = 100: 0.00152, k = 200: 0.00076, k = 400: 0.00038.

### window T = 5: λ = inf λ_min(M̄_T) = 0.0152, Λ = 0.5485, ρ = 19.65, ρ_D = 113.28, ρ_DD = 2011.45, Ψ_M = 1.400, Ψ_D = 3.753, Ψ_DD = 15.807, p₁ = 3.60, p₂ = 9.87, π₁ = 13.66, π₂ = 13.66, ‖N̄₂‖ = 1.054

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | BD | 3.459 | nan | 0.00759 | nan | 0.00076 | nan, nan | BD, nan |
| 20 | 0.1 | BD | 3.459 | nan | 0.00759 | nan | 0.00076 | nan, nan | BD, nan |
| 50 | 0.2 | ABD | 1.384 | nan | 0.00304 | nan | 0.00030 | nan, nan | BD, nan |
| 50 | 0.1 | BD | 1.384 | nan | 0.00304 | nan | 0.00030 | nan, nan | BD, nan |
| 100 | 0.2 | ABCD | 0.692 | 0.00001 | 0.00152 | 137.98 | 0.00015 | -0.00030, -0.00002 | ABD, nan |
| 100 | 0.1 | ABCD | 0.692 | 0.00003 | 0.00152 | 59.15 | 0.00015 | -0.00030, -0.00005 | ABD, nan |
| 200 | 0.2 | ABCD | 0.346 | 0.00003 | 0.00076 | 22.81 | 0.00008 | -0.00015, -0.00007 | ABCD, 0.00003 |
| 200 | 0.1 | ABCD | 0.346 | 0.00004 | 0.00076 | 18.63 | 0.00008 | -0.00015, -0.00008 | ABCD, 0.00003 |
| 400 | 0.2 | ABCD | 0.173 | 0.00002 | 0.00038 | 16.12 | 0.00004 | -0.00008, -0.00005 | ABCD, 0.00002 |
| 400 | 0.1 | ABCD | 0.173 | 0.00003 | 0.00038 | 13.90 | 0.00004 | -0.00008, -0.00005 | ABCD, 0.00003 |

### window T = 10: λ = inf λ_min(M̄_T) = 0.0238, Λ = 0.3569, ρ = 26.24, ρ_D = 148.13, ρ_DD = 3071.00, Ψ_M = 1.400, Ψ_D = 3.753, Ψ_DD = 15.807, p₁ = 3.36, p₂ = 8.18, π₁ = 20.66, π₂ = 23.82, ‖N̄₂‖ = 0.341

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | BD | 0.716 | nan | 0.00759 | nan | 0.00119 | nan, nan | BD, nan |
| 20 | 0.1 | BD | 0.716 | nan | 0.00759 | nan | 0.00119 | nan, nan | BD, nan |
| 50 | 0.2 | BCD | 0.287 | nan | 0.00304 | nan | 0.00048 | nan, nan | BD, nan |
| 50 | 0.1 | BCD | 0.287 | nan | 0.00304 | nan | 0.00048 | nan, nan | BD, nan |
| 100 | 0.2 | ABCD | 0.143 | 0.00014 | 0.00152 | 11.03 | 0.00024 | -0.00047, -0.00028 | ABD, nan |
| 100 | 0.1 | ABCD | 0.143 | 0.00016 | 0.00152 | 9.45 | 0.00024 | -0.00047, -0.00032 | ABCD, 0.00001 |
| 200 | 0.2 | ABCD | 0.072 | 0.00008 | 0.00076 | 9.10 | 0.00012 | -0.00024, -0.00017 | ABCD, 0.00006 |
| 200 | 0.1 | ABCD | 0.072 | 0.00010 | 0.00076 | 7.98 | 0.00012 | -0.00024, -0.00019 | ABCD, 0.00008 |
| 400 | 0.2 | ABCD | 0.036 | 0.00004 | 0.00038 | 8.47 | 0.00006 | -0.00012, -0.00009 | ABCD, 0.00004 |
| 400 | 0.1 | ABCD | 0.036 | 0.00005 | 0.00038 | 7.48 | 0.00006 | -0.00012, -0.00010 | ABCD, 0.00005 |

### window T = 20: λ = inf λ_min(M̄_T) = 0.0506, Λ = 0.2424, ρ = 23.78, ρ_D = 158.70, ρ_DD = 3118.95, Ψ_M = 1.400, Ψ_D = 3.753, Ψ_DD = 15.807, p₁ = 5.31, p₂ = 19.40, π₁ = 23.53, π₂ = 49.41, ‖N̄₂‖ = 0.422

| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | BD | 0.417 | nan | 0.00759 | nan | 0.00253 | nan, nan | BD, nan |
| 20 | 0.1 | BD | 0.417 | nan | 0.00759 | nan | 0.00253 | nan, nan | BD, nan |
| 50 | 0.2 | BCD | 0.167 | nan | 0.00304 | nan | 0.00101 | nan, nan | BD, nan |
| 50 | 0.1 | BCD | 0.167 | nan | 0.00304 | nan | 0.00101 | nan, nan | BD, nan |
| 100 | 0.2 | ABCD | 0.083 | 0.00030 | 0.00152 | 4.99 | 0.00051 | -0.00097, -0.00061 | ABD, nan |
| 100 | 0.1 | ABCD | 0.083 | 0.00035 | 0.00152 | 4.31 | 0.00051 | -0.00097, -0.00071 | ABD, nan |
| 200 | 0.2 | ABCD | 0.042 | 0.00018 | 0.00076 | 4.18 | 0.00025 | -0.00050, -0.00036 | ABCD, 0.00010 |
| 200 | 0.1 | ABCD | 0.042 | 0.00021 | 0.00076 | 3.68 | 0.00025 | -0.00050, -0.00041 | ABCD, 0.00013 |
| 400 | 0.2 | ABCD | 0.021 | 0.00010 | 0.00038 | 3.94 | 0.00013 | -0.00025, -0.00019 | ABCD, 0.00009 |
| 400 | 0.1 | ABCD | 0.021 | 0.00011 | 0.00038 | 3.49 | 0.00013 | -0.00025, -0.00022 | ABCD, 0.00010 |

Reading: wherever the lemma's conditions (A)–(D) hold, the certificate sup V̇/V ≤ −γ holds on the grid, and the guaranteed rate is (1 − δ − δ₂)λ_T/k up to the factor 1/(1 + εp₁ + ε²p₂), with λ_T the exponentially weighted PE level of the window and δ₂ the price of the residual second-order term, which decreases like 1/k. The window trades λ_T (closer to the long-time excitation for longer windows) against p₁ and p₂ (which grow with T) and the residual ‖N̄₂‖ (no monotone rule: it decreases with T on the quasi-periodic regressor and is smallest at T = 10 on the chaotic one); the ratio exact/guaranteed is the product of the window loss λ_∞/λ_T, the factor 1/(1 − δ − δ₂) and the O(1/k) factors. Where a condition fails the row says which: (A) fails when the derivative peak factors ρ_D, ρ_DD are large relative to the excitation (a chaotic regressor whose components are nearly collinear has a tiny λ_T and astronomical peak factors, and the threshold in k moves beyond the range), (C) fails when δ₂ = εn₂/λ_T is not small. The explicit closed-form constants are looser than the measured suprema and move the threshold to larger k.

Numbers are produced by `experiments/exp_pe_general.py`; nothing is edited by hand.
