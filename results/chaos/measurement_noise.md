# Proposition D‴.30: measurement noise at perfect tracking (`exp_measurement_noise.py`)

The three-bump instance at σ₀ = 0.15σ_crit = 0.00707, x_m = 2 sin t, the nearest-atom leak toward a fixed (η = 0) or an adaptive prior, at the lift of the exact representation of q = 1/3 with the prior at its atoms; the measurement y = x + ν in the controller and the adaptive law.

## A. The linear theory and the second-order force (no dynamics)

k_e = 25, N₀ = 1. Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩, a = ∂_xf̂ on the orbit: tr Q = 8.7836 = 8.1049 (⟨Ĵ̇ᵀĴ̇⟩) +0.0000 (cross) +0.6786 (⟨a²ĴᵀĴ⟩). Per mode (the eigenbasis of H): h_j, Q_jj, Q_jj/h_j, and the function-level shares of the floor, (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀), for each η.

| j | h_j | Q_jj | Q_jj/h_j | share at η = 0 (fixed) | share at η = 0.00707 | share at η = 0.07074 | share at η = ∞ |
|---|---|---|---|---|---|---|---|
| 1 | 3.0587e-04 | 3.1925e-02 | 104.37 | 2.6462e-05 | 3.2601e-04 | 5.8068e-04 | 6.3850e-04 |
| 2 | 1.9644e-03 | 1.0304e-01 | 52.45 | 4.4787e-04 | 1.1560e-03 | 1.8780e-03 | 2.0607e-03 |
| 3 | 3.1144e-03 | 1.3748e-01 | 44.14 | 8.4049e-04 | 1.6229e-03 | 2.5093e-03 | 2.7497e-03 |
| 4 | 1.5149e-02 | 2.7295e-01 | 18.02 | 3.7213e-03 | 4.1409e-03 | 5.0437e-03 | 5.4591e-03 |
| 5 | 2.1104e-02 | 3.1305e-01 | 14.83 | 4.6891e-03 | 5.0046e-03 | 5.8132e-03 | 6.2610e-03 |
| 6 | 3.2479e-02 | 1.7803e-01 | 5.48 | 2.9237e-03 | 3.0203e-03 | 3.3322e-03 | 3.5605e-03 |
| 7 | 8.8615e-02 | 2.2602e+00 | 25.51 | 4.1863e-02 | 4.2093e-02 | 4.3283e-02 | 4.5205e-02 |
| 8 | 1.3420e-01 | 2.4514e+00 | 18.27 | 4.6574e-02 | 4.6691e-02 | 4.7393e-02 | 4.9029e-02 |
| 9 | 1.9521e-01 | 3.0354e+00 | 15.55 | 5.8585e-02 | 5.8657e-02 | 5.9135e-02 | 6.0708e-02 |

| η | floor (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) | relative to η = ∞ | the weakest mode's parameter variance | its bias-removal rate λ_1^− |
|---|---|---|---|---|
| 0 (fixed) | 1.59670e-01 | 0.9089 | 8.6515e-02 | 0.000e+00 |
| 0.00707 | 1.62711e-01 | 0.9262 | 1.0658e+00 | 1.513e-04 |
| 0.07074 | 1.68968e-01 | 0.9618 | 1.8985e+00 | 2.780e-04 |
| ∞ | 1.75671e-01 | 1.0000 | 2.0875e+00 | 3.059e-04 |

At η = ∞ the floor is (N₀/2k_e)tr Q = 1.75671e-01 whatever H. The second-order force (errors in variables; OU noise of rate B = 100), slow time, per unit of the noise's variance, in the eigenbasis of H: +2.133e-02, +1.880e+00, +6.659e-01, +1.583e-01, +4.332e-01, +9.929e-01, -4.724e+00, -4.150e-02, +5.619e+00. The stationary mean deviation's function-level size per unit of the variance squared: fixed prior Σ_jF_j²h_j/(h_j + σ₀)² = 4.9106e+02, adaptive prior (every η > 0) Σ_jF_j²/h_j = 2.3969e+03 (ratio 4.88).

The fast parts: tr H² = 0.0657 (the function-level floor) and ⟨|Ĵ|⁴⟩ = 0.2991 (the compensation error at the current point), each times N₀/(2k_e). With the OU noise (N₀ = 2Eν²/B) the errors-in-variables bias equals the floor of the parameters at Eν² = 6.50e-06 (std 0.0026) under the fixed prior and 1.47e-06 (std 0.0012) under the adaptive one, and dominates above.

The full second-order mean (`second_order_mean`: the periodic covariance of the linearised loop with the OU noise as a state, contracted with the second derivatives of the true vector field). Its noise part — the noise's own second moments at fixed parameters — equals the quasi-static force above to 1.7e-04 (relative); the parameters' first-order fluctuation adds a force of the same order in Eν² (the compensation's curvature and slope–parameter terms in the error equation, ∇²_ϑΦ̃E[ε_mδϑ] in the law). Function-level biases per unit (Eν²)²: of the mean parameters (‖Ĵȳ‖²_ρ, the averaged response) from the noise part and from the full force, of the mean function (‖Ĵȳ + ½D²_Θf̂[Σ_ΘΘ]‖²_ρ), and of the exact periodic second-order mean of the linearised loop (mean parameters).

| prior | noise part | full force | mean function | exact periodic mean | periodicity |
|---|---|---|---|---|---|
| fixed | 491.112 | 491.11 | 491.114 | 503.72 | 8.9e-13 |
| adaptive (η = ∞; the noise part for every η > 0) | 2397.66 | 2397.68 | 2397.65 | 2652.96 | 5.9e-13 |

The adaptive prior's averaged bias per unit (Eν²)² against k_e at the OU rates B, with its local exponents d log/d log k_e between successive k_e (→ 2 while k_e ≪ B, saturating for k_e ≫ B).

| B | k_e = 12.5 | k_e = 25 | k_e = 50 | k_e = 100 | k_e = 200 | k_e = 400 | local exponents |
|---|---|---|---|---|---|---|---|
| 100 | 742.3 | 2397 | 6649 | 1.495e+04 | 2.657e+04 | 3.825e+04 | 1.69, 1.47, 1.17, 0.83, 0.53 |
| 1000 | 915.6 | 3564 | 1.357e+04 | 4.943e+04 | 1.661e+05 | 4.881e+05 | 1.96, 1.93, 1.87, 1.75, 1.56 |
| 10000 | 936.2 | 3725 | 1.481e+04 | 5.863e+04 | 2.299e+05 | 8.845e+05 | 1.99, 1.99, 1.98, 1.97, 1.94 |

Invariance is not attraction: a separation of class k decays at σ₀/k_e = 2.830e-04 without noise and at σ₀/k_e − Eν²λ_max(S_k) with it (S_k the averaged second-order drift, real time, `separation_drift`). λ_max(S_k) per unit Eν²: 0.0148, 0.0101, 0.2981; the separations stop decaying at the std 0.138, 0.168, 0.031 (classes 1, 2, 3); at Part C's stds the margin is used to 0.9%, 10.5% (the weakest class).

## B. The exactly linearised loop with white measurement noise against the averaged prediction

The periodic stationary covariance P(t) (the monodromy and a discrete Lyapunov equation; N₀ = 1), averaged over a period: the function-level floor E‖f̂ − f‖²_ρ = ⟨tr(HP_yy)⟩ and the compensation error at the plant's current point ⟨ĴP_yyĴᵀ⟩, against the prediction — the slow (parameter) part (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) plus the fast part (N₀/2k_e)tr H² and (N₀/2k_e)⟨|Ĵ|⁴⟩ respectively.

| k_e | η | E‖f̂ − f‖²_ρ: exact | predicted | ratio | at the current point: exact | predicted | ratio | the slow part alone |
|---|---|---|---|---|---|---|---|---|
| 25 | 0 (fixed) | 1.57161e-01 | 1.60985e-01 | 0.9762 | 1.61703e-01 | 1.65653e-01 | 0.9762 | 1.59670e-01 |
| 50 | 0 (fixed) | 7.99990e-02 | 8.04923e-02 | 0.9939 | 8.23123e-02 | 8.28263e-02 | 0.9938 | 7.98351e-02 |
| 100 | 0 (fixed) | 4.01892e-02 | 4.02461e-02 | 0.9986 | 4.13526e-02 | 4.14131e-02 | 0.9985 | 3.99176e-02 |
| 200 | 0 (fixed) | 2.01175e-02 | 2.01231e-02 | 0.9997 | 2.07002e-02 | 2.07066e-02 | 0.9997 | 1.99588e-02 |
| 400 | 0 (fixed) | 1.00612e-02 | 1.00615e-02 | 1.0000 | 1.03528e-02 | 1.03533e-02 | 1.0000 | 9.97939e-03 |
| 800 | 0 (fixed) | 5.03082e-03 | 5.03077e-03 | 1.0000 | 5.17666e-03 | 5.17664e-03 | 1.0000 | 4.98970e-03 |
| 25 | 0.00707 | 1.60323e-01 | 1.64025e-01 | 0.9774 | 1.64871e-01 | 1.68693e-01 | 0.9773 | 1.62711e-01 |
| 50 | 0.00707 | 8.15347e-02 | 8.20126e-02 | 0.9942 | 8.38487e-02 | 8.43466e-02 | 0.9941 | 8.13554e-02 |
| 100 | 0.00707 | 4.09513e-02 | 4.10063e-02 | 0.9987 | 4.21148e-02 | 4.21733e-02 | 0.9986 | 4.06777e-02 |
| 200 | 0.00707 | 2.04978e-02 | 2.05031e-02 | 0.9997 | 2.10806e-02 | 2.10866e-02 | 0.9997 | 2.03389e-02 |
| 400 | 0.00707 | 1.02513e-02 | 1.02516e-02 | 1.0000 | 1.05429e-02 | 1.05433e-02 | 1.0000 | 1.01694e-02 |
| 800 | 0.00707 | 5.12585e-03 | 5.12579e-03 | 1.0000 | 5.27169e-03 | 5.27166e-03 | 1.0000 | 5.08472e-03 |
| 25 | ∞ | 1.73501e-01 | 1.76985e-01 | 0.9803 | 1.78066e-01 | 1.81653e-01 | 0.9803 | 1.75671e-01 |
| 50 | ∞ | 8.80421e-02 | 8.84927e-02 | 0.9949 | 9.03585e-02 | 9.08267e-02 | 0.9948 | 8.78356e-02 |
| 100 | ∞ | 4.41949e-02 | 4.42464e-02 | 0.9988 | 4.53586e-02 | 4.54134e-02 | 0.9988 | 4.39178e-02 |
| 200 | ∞ | 2.21183e-02 | 2.21232e-02 | 0.9998 | 2.27011e-02 | 2.27067e-02 | 0.9998 | 2.19589e-02 |
| 400 | ∞ | 1.10614e-02 | 1.10616e-02 | 1.0000 | 1.13530e-02 | 1.13533e-02 | 1.0000 | 1.09795e-02 |
| 800 | ∞ | 5.53087e-03 | 5.53080e-03 | 1.0000 | 5.67671e-03 | 5.67667e-03 | 1.0000 | 5.48973e-03 |

The weakest mode's variance (period average, N₀ = 1) against its slow part (N₀/2k_e)Q₁₁(h₁ + η)/(h₁(h₁ + η + σ₀)) plus its fast part (N₀/2k_e)h₁, and the ratio of the adaptive (η = ∞) to the fixed prior's variance, exact and averaged.

| k_e | η | exact | slow + fast | ratio |
|---|---|---|---|---|
| 25 | 0 (fixed) | 7.30756e-02 | 8.65208e-02 | 0.8446 |
| 50 | 0 (fixed) | 4.13242e-02 | 4.32604e-02 | 0.9552 |
| 100 | 0 (fixed) | 2.13782e-02 | 2.16302e-02 | 0.9884 |
| 200 | 0 (fixed) | 1.07832e-02 | 1.08151e-02 | 0.9971 |
| 400 | 0 (fixed) | 5.40354e-03 | 5.40755e-03 | 0.9993 |
| 800 | 0 (fixed) | 2.70327e-03 | 2.70377e-03 | 0.9998 |
| 25 | 0.00707 | 1.04118e+00 | 1.06584e+00 | 0.9769 |
| 50 | 0.00707 | 5.29312e-01 | 5.32921e-01 | 0.9932 |
| 100 | 0.00707 | 2.65987e-01 | 2.66460e-01 | 0.9982 |
| 200 | 0.00707 | 1.33170e-01 | 1.33230e-01 | 0.9995 |
| 400 | 0.00707 | 6.66074e-02 | 6.66151e-02 | 0.9999 |
| 800 | 0.00707 | 3.33065e-02 | 3.33075e-02 | 1.0000 |
| 25 | ∞ | 2.04529e+00 | 2.08750e+00 | 0.9798 |
| 50 | ∞ | 1.03756e+00 | 1.04375e+00 | 0.9941 |
| 100 | ∞ | 5.21062e-01 | 5.21876e-01 | 0.9984 |
| 200 | ∞ | 2.60834e-01 | 2.60938e-01 | 0.9996 |
| 400 | ∞ | 1.30456e-01 | 1.30469e-01 | 0.9999 |
| 800 | ∞ | 6.52328e-02 | 6.52345e-02 | 1.0000 |

Adaptive over fixed, the weakest mode: exact 27.99, 25.11, 24.37, 24.19, 24.14, 24.13 (k_e = 25, 50, 100, 200, 400, 800); averaged (slow parts) 24.13.

The deviation of the floor ratio from 1 (η = 0 and ∞) with the noise input −(k_e + a)ν and in the control with −k_eν (the slope a removed, the generator unchanged; its prediction with Q = ⟨Ĵ̇ᵀĴ̇⟩), and the deviation's successive ratios per doubling of k_e (→ 4 for an even expansion in 1/k_e, → 2 for a first-order term).

| η | noise input | ratio − 1 at k_e = 25, 50, 100, 200, 400, 800 | successive ratios |
|---|---|---|---|
| 0 (fixed) | −(k_e + a)ν | -2.375e-02, -6.128e-03, -1.413e-03, -2.783e-04, -3.098e-05, +1.162e-05 | 3.88, 4.34, 5.08, 8.98, -2.67 |
| 0 (fixed) | −k_eν (control) | -2.626e-02, -7.109e-03, -1.817e-03, -4.569e-04, -1.144e-04, -2.861e-05 | 3.69, 3.91, 3.98, 3.99, 4.00 |
| ∞ | −(k_e + a)ν | -1.969e-02, -5.092e-03, -1.163e-03, -2.210e-04, -1.946e-05, +1.309e-05 | 3.87, 4.38, 5.26, 11.36, -1.49 |
| ∞ | −k_eν (control) | -2.191e-02, -5.975e-03, -1.531e-03, -3.850e-04, -9.641e-05, -2.411e-05 | 3.67, 3.90, 3.98, 3.99, 4.00 |

The exact stationary second-order mean of the linearised loop with the OU noise (B = 100, the measured regressor): the function-level bias of its mean parameters over a period against the averaged bias of the quasi-static force, per unit (Eν²)².

| k_e | prior | averaged | exact | exact/averaged | periodicity |
|---|---|---|---|---|---|
| 25 | 0 (fixed) | 491.064 | 503.72 | 1.0258 | 8.9e-13 |
| 50 | 0 (fixed) | 1365.85 | 1374.64 | 1.0064 | 4.4e-14 |
| 100 | 0 (fixed) | 3075.6 | 3080.56 | 1.0016 | 5.1e-14 |
| 200 | 0 (fixed) | 5470.1 | 5472.31 | 1.0004 | 9.5e-15 |
| 25 | ∞ | 2396.95 | 2652.96 | 1.1068 | 5.9e-13 |
| 50 | ∞ | 6648.54 | 6824.75 | 1.0265 | 2.2e-12 |
| 100 | ∞ | 14949.3 | 15048.5 | 1.0066 | 7.2e-13 |
| 200 | ∞ | 26568.3 | 26612.4 | 1.0017 | 1.9e-13 |

## C. The real loop with an Ornstein–Uhlenbeck measurement noise (k_e = 25, B = 100, RK4 dt = 0.001, 300 time units)

E‖f̂ − f‖²_ρ at the records (mean ± standard error over the runs) against the linear prediction: the transient covariance of the exactly linearised loop with the OU noise plus the second-order mean deviation (their sum, and the two parts).

Noise std 0.003 (variance 9.0e-06), η = 0 (fixed), 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 6.6096e-09 ± 1.9e-09 | 7.3368e-09 | 6.7785e-09 | 5.5832e-10 |
| 50 | 1.1271e-08 ± 2.4e-09 | 1.3615e-08 | 1.1720e-08 | 1.8949e-09 |
| 75 | 2.5858e-08 ± 9.1e-09 | 1.8956e-08 | 1.5312e-08 | 3.6437e-09 |
| 100 | 3.2351e-08 ± 8.9e-09 | 2.3547e-08 | 1.7972e-08 | 5.5746e-09 |
| 125 | 3.0619e-08 ± 7.0e-09 | 2.7362e-08 | 1.9816e-08 | 7.5457e-09 |
| 150 | 3.0917e-08 ± 7.9e-09 | 3.0681e-08 | 2.1208e-08 | 9.4728e-09 |
| 175 | 3.3660e-08 ± 6.1e-09 | 3.3732e-08 | 2.2424e-08 | 1.1308e-08 |
| 200 | 3.5042e-08 ± 5.8e-09 | 3.6488e-08 | 2.3461e-08 | 1.3027e-08 |
| 225 | 3.7653e-08 ± 7.7e-09 | 3.8899e-08 | 2.4278e-08 | 1.4621e-08 |
| 250 | 3.5510e-08 ± 7.4e-09 | 4.0995e-08 | 2.4908e-08 | 1.6087e-08 |
| 275 | 4.3328e-08 ± 8.5e-09 | 4.2835e-08 | 2.5405e-08 | 1.7430e-08 |
| 300 | 3.7097e-08 ± 7.7e-09 | 4.4468e-08 | 2.5810e-08 | 1.8657e-08 |

Noise std 0.003 (variance 9.0e-06), η = ∞, 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 6.6504e-09 ± 1.9e-09 | 7.3833e-09 | 6.8211e-09 | 5.6218e-10 |
| 50 | 1.1399e-08 ± 2.4e-09 | 1.3785e-08 | 1.1864e-08 | 1.9204e-09 |
| 75 | 2.6229e-08 ± 9.2e-09 | 1.9301e-08 | 1.5585e-08 | 3.7154e-09 |
| 100 | 3.3026e-08 ± 9.1e-09 | 2.4099e-08 | 1.8382e-08 | 5.7171e-09 |
| 125 | 3.1524e-08 ± 7.2e-09 | 2.8139e-08 | 2.0358e-08 | 7.7809e-09 |
| 150 | 3.1970e-08 ± 8.1e-09 | 3.1692e-08 | 2.1874e-08 | 9.8183e-09 |
| 175 | 3.4940e-08 ± 6.4e-09 | 3.4980e-08 | 2.3202e-08 | 1.1778e-08 |
| 200 | 3.6446e-08 ± 6.0e-09 | 3.7971e-08 | 2.4341e-08 | 1.3630e-08 |
| 225 | 3.9364e-08 ± 8.0e-09 | 4.0614e-08 | 2.5250e-08 | 1.5364e-08 |
| 250 | 3.7337e-08 ± 7.8e-09 | 4.2936e-08 | 2.5961e-08 | 1.6974e-08 |
| 275 | 4.5534e-08 ± 8.9e-09 | 4.4995e-08 | 2.6531e-08 | 1.8464e-08 |
| 300 | 3.9309e-08 ± 8.1e-09 | 4.6839e-08 | 2.7001e-08 | 1.9837e-08 |

Noise std 0.01 (variance 1.0e-04), η = 0 (fixed), 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 1.5887e-07 ± 3.2e-08 | 1.4425e-07 | 7.5317e-08 | 6.8929e-08 |
| 50 | 4.0936e-07 ± 5.1e-08 | 3.6416e-07 | 1.3022e-07 | 2.3394e-07 |
| 75 | 8.1396e-07 ± 1.7e-07 | 6.1997e-07 | 1.7013e-07 | 4.4984e-07 |
| 100 | 1.1820e-06 ± 1.7e-07 | 8.8791e-07 | 1.9969e-07 | 6.8822e-07 |
| 125 | 1.3691e-06 ± 1.7e-07 | 1.1518e-06 | 2.2018e-07 | 9.3157e-07 |
| 150 | 1.5286e-06 ± 1.8e-07 | 1.4051e-06 | 2.3565e-07 | 1.1695e-06 |
| 175 | 1.7691e-06 ± 1.6e-07 | 1.6452e-06 | 2.4915e-07 | 1.3961e-06 |
| 200 | 1.9403e-06 ± 1.8e-07 | 1.8690e-06 | 2.6067e-07 | 1.6083e-06 |
| 225 | 2.1198e-06 ± 2.0e-07 | 2.0748e-06 | 2.6976e-07 | 1.8050e-06 |
| 250 | 2.2769e-06 ± 2.1e-07 | 2.2628e-06 | 2.7675e-07 | 1.9860e-06 |
| 275 | 2.5985e-06 ± 2.4e-07 | 2.4341e-06 | 2.8227e-07 | 2.1519e-06 |
| 300 | 2.5140e-06 ± 2.6e-07 | 2.5902e-06 | 2.8678e-07 | 2.3034e-06 |

Noise std 0.01 (variance 1.0e-04), η = ∞, 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 1.5990e-07 ± 3.2e-08 | 1.4519e-07 | 7.5790e-08 | 6.9404e-08 |
| 50 | 4.1443e-07 ± 5.2e-08 | 3.6891e-07 | 1.3183e-07 | 2.3708e-07 |
| 75 | 8.2804e-07 ± 1.7e-07 | 6.3186e-07 | 1.7317e-07 | 4.5869e-07 |
| 100 | 1.2094e-06 ± 1.7e-07 | 9.1005e-07 | 2.0424e-07 | 7.0581e-07 |
| 125 | 1.4107e-06 ± 1.8e-07 | 1.1868e-06 | 2.2620e-07 | 9.6060e-07 |
| 150 | 1.5839e-06 ± 1.9e-07 | 1.4552e-06 | 2.4305e-07 | 1.2121e-06 |
| 175 | 1.8411e-06 ± 1.7e-07 | 1.7118e-06 | 2.5780e-07 | 1.4540e-06 |
| 200 | 2.0287e-06 ± 1.9e-07 | 1.9532e-06 | 2.7046e-07 | 1.6828e-06 |
| 225 | 2.2266e-06 ± 2.1e-07 | 2.1774e-06 | 2.8056e-07 | 1.8968e-06 |
| 250 | 2.4013e-06 ± 2.2e-07 | 2.3841e-06 | 2.8846e-07 | 2.0956e-06 |
| 275 | 2.7458e-06 ± 2.5e-07 | 2.5743e-06 | 2.9479e-07 | 2.2795e-06 |
| 300 | 2.6739e-06 ± 2.7e-07 | 2.7491e-06 | 3.0001e-07 | 2.4490e-06 |

The adaptive prior's E‖f̂ − f‖²_ρ over the fixed prior's, Monte Carlo (predicted): std 0.003: +0.6% (+0.6%), +3.4% (+3.3%), +6.0% (+5.3%) at t = 25, 150, 300; std 0.01: +0.6% (+0.7%), +3.6% (+3.6%), +6.4% (+6.1%) at t = 25, 150, 300. All four cases run on the same 12 noise paths (the seeds 1000 + run, scaled by the std), so their agreements are not independent.

The lift m = 15 (counts 5, 5, 5) driven by the same noise path as the K₀-atom run (seed 0, std 0.01, η = ∞): the largest difference of the stored f̂ (float32) on the grid over the records 0.0e+00, of E‖f̂ − f‖² 6.0e-19 — the noise is coherent: a replicated start stays replicated (in float64: Part D).

## D. Invariance is not attraction

The lift and the K₀-atom run on one noise path (std 0.01, η = ∞, 25 time units) in float64: the largest difference of f̂ on the grid 5.6e-16 (a mean of five against a weighted atom: rounding), the lift's largest within-class spread 0.0e+00 (identical arithmetic per particle, whether or not the replicated state attracts).

The lift m = 15 with class 3's particles separated by 10⁻⁶·(2, 1, 0, −1, −2) along S₃'s top eigenvector, 8 antithetic pairs of noise paths per std (ν and −ν), 300 time units: the separation's growth rate (least squares on log‖separation‖ over t ≥ 50) — its pair means cancel the first-order part — and λ_eff = (rate + σ₀/k_e)/Eν², against the lift's drift (−σ₀/k_e + λ_max(S₃)Eν²) and against the drift along the second-order transient from the lift (the slope predicted for u̇ = (Eν²S₃(t) − σ₀/k_e)u).

The drift along the transient (Eν² = 1): per 50-unit window, the mean error's rms and the largest eigenvalues of the window-averaged S_k (classes 1, 2, 3); at the lift: 0.0148, 0.0101, 0.2981.

| window | rms E[e] | λ_max(S₁), λ_max(S₂), λ_max(S₃) |
|---|---|---|
| [0, 50) | 0.121 | 0.0095, 0.0013, 0.3069 |
| [50, 100) | 0.273 | 0.0029, 0.0070, 0.3206 |
| [100, 150) | 0.393 | 0.0009, 0.0133, 0.3311 |
| [150, 200) | 0.483 | 0.0002, 0.0143, 0.3391 |
| [200, 250) | 0.550 | -0.0001, 0.0157, 0.3451 |
| [250, 300) | 0.602 | 0.0002, 0.0160, 0.3494 |

The separation's rate vanishes at the std 0.0308 by the lift's drift and at 0.0294 by the transient's over t ∈ [50, 300].

| std | lift's prediction | transient's prediction | measured (pair means ± s.e.) | λ_eff measured | λ_eff: lift; transient |
|---|---|---|---|---|---|
| 0.01 | -2.532e-04 | -2.503e-04 | -2.5027e-04 ± 3.8e-07 | 0.3271 ± 0.0038 | 0.2981; 0.3264 |
| 0.03 | -1.464e-05 | +1.132e-05 | +2.7358e-06 ± 3.3e-06 | 0.3175 ± 0.0037 | 0.2981; 0.3270 |
| 0.05 | +4.624e-04 | +5.372e-04 | +4.5267e-04 ± 8.4e-06 | 0.2943 ± 0.0034 | 0.2981; 0.3281 |

The weighted fit λ_eff = c₀ + c₁Eν²: c₀ = 0.3291 ± 0.0033 (the small-noise limit; the lift's drift 0.2981, the transient's 0.3264 at std 0.01), c₁ = -13.82 ± 2.04 (the next order); the measured rate vanishes at the std 0.0299.

Numbers are produced by `experiments/exp_measurement_noise.py`; nothing is edited by hand.
