# Proposition D‴.31: the adaptive prior is a washout leak; the desired regressor (`exp_eiv_remedy.py`)

The instance of Proposition D‴.30 (σ₀ = 0.00707, x_m = 2 sin t, the lift of the exact representation of q = 1/3 with the prior at its atoms) under the measurement y = x + ν, with an Ornstein–Uhlenbeck noise of rate B. MR: the measured regressor and compensation (Proposition D‴.30); DR: the regressor at the reference, the compensation f̂(y); DR*: both at the reference (the desired compensation adaptation law).

## A. The second-order force and its stationary responses (no dynamics)

Per unit (Eν²)², function level, the averaged stationary response: the bias of the mean parameters under the fixed prior and under the adaptive prior (η = ∞, the noise part's for every η > 0; = FᵀH⁻¹F, the squared norm of the smallest target perturbation exerting F) for the force of the noise's own second moments (quasi-static, `eiv_force`) and for the full second-order force (the parameters' fluctuation included, `second_order_mean`), the bias of the mean function under the adaptive prior with the full force (‖Ĵȳ + ½D²_Θf̂[Σ_ΘΘ]‖²_ρ), the least-squares bound ‖g‖²_ρ (g = k_eE[e] of the noise part; DR, DR*), the rescaled projected mean tracking error at the adaptive prior's stationary point, k_e²‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁ (F₁ the regressor's errors in variables; zero by construction for DR and DR*), and the mean tracking error itself ‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁/k_e². The last column is the noise variance at which the adaptive prior's full bias equals its floor of the parameters (N₀/2k_e)tr Q, N₀ = 2Eν²/B (std above 0.1 is outside the small-noise regime).

| k_e | B | law | fixed: noise part | fixed: full | adaptive: noise part | adaptive: full | adaptive: mean function | ‖g‖²_ρ | k_e²‖P E[e]‖²_ρ | ‖P E[e]‖²_ρ | Eν² at bias = floor (std) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 100 | MR | 4.9106e+02 | 4.9111e+02 | 2.3969e+03 | 2.3977e+03 | 2.3977e+03 | — | 2.394e+03 | 3.830e+00 | 1.47e-06 (0.0012) |
| 25 | 100 | DR | 2.1060e-01 | 2.1469e-01 | 4.7053e-01 | 4.8584e-01 | 4.7591e-01 | 4.7057e-01 | 0.000e+00 | 0.000e+00 | 7.23e-03 (0.0850) |
| 25 | 100 | DR* | 2.3601e-02 | 2.2303e-02 | 5.2621e-02 | 4.8050e-02 | 5.1114e-02 | 5.2663e-02 | 0.000e+00 | 0.000e+00 | 7.31e-02 (0.2704) |
| 25 | 1000 | MR | 7.2955e+02 | 7.2956e+02 | 3.5637e+03 | 3.5638e+03 | 3.5638e+03 | — | 3.557e+03 | 5.691e+00 | 9.86e-08 (0.0003) |
| 25 | 1000 | DR | 5.2899e-01 | 5.2964e-01 | 1.1821e+00 | 1.1845e+00 | 1.1830e+00 | 1.1821e+00 | 0.000e+00 | 0.000e+00 | 2.97e-04 (0.0172) |
| 25 | 1000 | DR* | 3.5056e-04 | 3.3459e-04 | 7.8187e-04 | 7.2490e-04 | 7.6287e-04 | 7.8240e-04 | 0.000e+00 | 0.000e+00 | 4.85e-01 (0.6962) |
| 50 | 100 | MR | 1.3658e+03 | 1.3659e+03 | 6.6485e+03 | 6.6494e+03 | 6.6494e+03 | — | 6.646e+03 | 2.658e+00 | 2.64e-07 (0.0005) |
| 50 | 100 | DR | 6.5041e-02 | 6.6199e-02 | 1.4529e-01 | 1.4971e-01 | 1.4691e-01 | 1.4531e-01 | 0.000e+00 | 0.000e+00 | 1.17e-02 (0.1083) |
| 50 | 100 | DR* | 6.5109e-02 | 6.3987e-02 | 1.4542e-01 | 1.4122e-01 | 1.4391e-01 | 1.4545e-01 | 0.000e+00 | 0.000e+00 | 1.24e-02 (0.1115) |
| 50 | 1000 | MR | 2.7847e+03 | 2.7847e+03 | 1.3569e+04 | 1.3569e+04 | 1.3569e+04 | — | 1.356e+04 | 5.423e+00 | 1.29e-08 (0.0001) |
| 50 | 1000 | DR | 4.7858e-01 | 4.7890e-01 | 1.0695e+00 | 1.0707e+00 | 1.0699e+00 | 1.0695e+00 | 0.000e+00 | 0.000e+00 | 1.64e-04 (0.0128) |
| 50 | 1000 | DR* | 1.3281e-03 | 1.3120e-03 | 2.9666e-03 | 2.9057e-03 | 2.9447e-03 | 2.9671e-03 | 0.000e+00 | 0.000e+00 | 6.05e-02 (0.2459) |
| 100 | 100 | MR | 3.0756e+03 | 3.0756e+03 | 1.4949e+04 | 1.4950e+04 | 1.4950e+04 | — | 1.495e+04 | 1.495e+00 | 5.88e-08 (0.0002) |
| 100 | 100 | DR | 3.1075e-05 | 3.4499e-05 | 4.5092e-05 | 7.0431e-05 | 5.7438e-05 | 5.8406e-05 | 0.000e+00 | 0.000e+00 | 1.25e+01 (3.5315) |
| 100 | 100 | DR* | 1.4623e-01 | 1.4538e-01 | 3.2675e-01 | 3.2351e-01 | 3.2556e-01 | 3.2677e-01 | 0.000e+00 | 0.000e+00 | 2.72e-03 (0.0521) |
| 100 | 1000 | MR | 1.0158e+04 | 1.0158e+04 | 4.9429e+04 | 4.9429e+04 | 4.9429e+04 | — | 4.941e+04 | 4.941e+00 | 1.78e-09 (0.0000) |
| 100 | 1000 | DR | 3.9137e-01 | 3.9151e-01 | 8.7459e-01 | 8.7513e-01 | 8.7479e-01 | 8.7459e-01 | 0.000e+00 | 0.000e+00 | 1.00e-04 (0.0100) |
| 100 | 1000 | DR* | 4.8334e-03 | 4.8178e-03 | 1.0800e-02 | 1.0741e-02 | 1.0778e-02 | 1.0801e-02 | 0.000e+00 | 0.000e+00 | 8.18e-03 (0.0904) |
| 200 | 100 | MR | 5.4701e+03 | 5.4701e+03 | 2.6568e+04 | 2.6569e+04 | 2.6569e+04 | — | 2.657e+04 | 6.643e-01 | 1.65e-08 (0.0001) |
| 200 | 100 | DR | 6.4958e-02 | 6.4671e-02 | 1.4516e-01 | 1.4407e-01 | 1.4476e-01 | 1.4517e-01 | 0.000e+00 | 0.000e+00 | 3.05e-03 (0.0552) |
| 200 | 100 | DR* | 2.5986e-01 | 2.5928e-01 | 5.8069e-01 | 5.7851e-01 | 5.7988e-01 | 5.8070e-01 | 0.000e+00 | 0.000e+00 | 7.59e-04 (0.0276) |
| 200 | 1000 | MR | 3.4156e+04 | 3.4156e+04 | 1.6610e+05 | 1.6610e+05 | 1.6610e+05 | — | 1.661e+05 | 4.152e+00 | 2.64e-10 (0.0000) |
| 200 | 1000 | DR | 2.5984e-01 | 2.5990e-01 | 5.8066e-01 | 5.8088e-01 | 5.8074e-01 | 5.8066e-01 | 0.000e+00 | 0.000e+00 | 7.56e-05 (0.0087) |
| 200 | 1000 | DR* | 1.6240e-02 | 1.6226e-02 | 3.6292e-02 | 3.6237e-02 | 3.6272e-02 | 3.6292e-02 | 0.000e+00 | 0.000e+00 | 1.21e-03 (0.0348) |

The parameters' fluctuation's share of the adaptive prior's bias (full over noise part − 1): MR +0.03%, +0.00%, +0.01%, +0.00%, +0.01%, +0.00%, +0.00%, +0.00%; DR +3.25%, +0.20%, +3.05%, +0.11%, +56.19%, +0.06%, -0.75%, +0.04%; DR* -8.69%, -7.29%, -2.89%, -2.05%, -0.99%, -0.55%, -0.38%, -0.15% (the rows' order); of the mean function's: MR +0.03%, +0.00%, +0.01%, +0.00%, +0.01%, +0.00%, +0.00%, +0.00%; DR +1.14%, +0.07%, +1.12%, +0.04%, +27.38%, +0.02%, -0.28%, +0.01%; DR* -2.86%, -2.43%, -1.03%, -0.74%, -0.36%, -0.20%, -0.14%, -0.06%.

The parameters' part's own function-level response under the adaptive prior (and its learning-law part's), per unit (Eν²)²: MR 3.92e-04 (2.10e-04), 4.01e-06 (2.15e-06), 1.13e-04 (6.34e-05), 1.15e-06 (6.51e-07), 2.93e-05 (1.67e-05), 3.00e-07 (1.72e-07), 7.40e-06 (4.23e-06), 7.58e-08 (4.35e-08); DR 3.92e-04 (2.10e-04), 4.01e-06 (2.15e-06), 1.13e-04 (6.34e-05), 1.15e-06 (6.51e-07), 2.93e-05 (1.67e-05), 3.00e-07 (1.72e-07), 7.40e-06 (4.23e-06), 7.58e-08 (4.35e-08); DR* 3.93e-04 (2.11e-04), 4.02e-06 (2.17e-06), 1.13e-04 (6.35e-05), 1.15e-06 (6.52e-07), 2.93e-05 (1.67e-05), 3.00e-07 (1.72e-07), 7.40e-06 (4.24e-06), 7.58e-08 (4.35e-08) (the rows' order) — set by the parameters' covariance, falling with k_e.

The adaptive prior's full bias, MR over DR: 4.94e+03 (k_e = 25, B = 100), 3.01e+03 (k_e = 25, B = 1000), 4.44e+04 (k_e = 50, B = 100), 1.27e+04 (k_e = 50, B = 1000), 2.12e+08 (k_e = 100, B = 100), 5.65e+04 (k_e = 100, B = 1000), 1.84e+05 (k_e = 200, B = 100), 2.86e+05 (k_e = 200, B = 1000). MR's grows from k_e = 25 to 200 at B = 1000 by 46.6 — the law (k_eB/(k_e + B))² gives 46.7, k_e² 64; the noise parts of DR's and DR*'s stay below ¼‖∂²_xf‖²_ρ = 1.3065e+00 for every k_e and B, and their full biases at most 1.1845 here.

## B. The stationary identity RH + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a)

η = σ₀; R and B from the averaged collective block at zero frequency (k_e = 25, B = 100). h₁ = 3.0587e-04, h_max = 1.9521e-01. The errors-in-variables biases are function-level per unit (Eν²)²; the anchor's biases are function-level for a unit misplacement along the weakest and the strongest mode.

| prior | S₀ (zero-frequency leak) | weakest mode: response r₁ | bias factor b₁ | h₁r₁ + b₁ | EIV bias, MR | EIV bias, DR | anchor bias, weakest | anchor bias, strongest | max|RH + B − I| |
|---|---|---|---|---|---|---|---|---|---|
| adaptive, κ = 0 | 0.0000e+00 | 3.2694e+03 | 0.0000 | 1.000000 | 2.3969e+03 | 4.7053e-01 | 0.000e+00 | 0.000e+00 | 1.1e-14 |
| κ = 0.1η | 6.4313e-04 | 1.0537e+03 | 0.6777 | 1.000000 | 1.5667e+03 | 4.0005e-01 | 1.405e-04 | 2.105e-06 | 2.7e-15 |
| κ = 1η | 3.5372e-03 | 2.6021e+02 | 0.9204 | 1.000000 | 6.8172e+02 | 2.8736e-01 | 2.591e-04 | 6.183e-05 | 1.1e-15 |
| κ = 10η | 6.4313e-03 | 1.4843e+02 | 0.9546 | 1.000000 | 5.1129e+02 | 2.2184e-01 | 2.787e-04 | 1.986e-04 | 8.9e-16 |
| fixed (η = 0) | 7.0744e-03 | 1.3550e+02 | 0.9586 | 1.000000 | 4.9106e+02 | 2.1060e-01 | 2.810e-04 | 2.388e-04 | 4.4e-16 |

Along the family the zero-frequency leak S₀ = σ₀κ/(η + κ) trades the anchor's bias for the response to persistent forces, mode by mode, at h_jr_j + b_j = 1; the trade binds on the modes with h_j < σ₀ (3 of 9 here), where 1/h_j exceeds the fixed prior's 1/σ₀. The desired regressor changes the force instead: MR's bias over DR's (noise parts) along the family 5.09e+03, 3.92e+03, 2.37e+03, 2.3e+03, 2.33e+03.

Priors whose zero-frequency leak S₀ does not commute with H: R and B from the dynamic block of (y, z_f) at zero frequency (`anchored_block`; the block's dimension), against (H + S₀)⁻¹ and (H + S₀)⁻¹S₀ with S₀ = σ₀[E_xE_xᵀ + E_f(ηI + K)⁻¹KE_fᵀ] by elimination of z_f (the largest relative difference of R, the largest difference of B): RH + B = I holds, HR + B = I does not; per mode h_jr_jj + b_jj = 1 holds whatever S₀.

| prior | block dimension | block against the definitions | max|HS₀ − S₀H| | max|RH + B − I| | max|HR + B − I| | max_j|h_jr_jj + b_jj − 1| | EIV bias, MR | EIV bias, DR |
|---|---|---|---|---|---|---|---|---|
| fixed on the amplitudes, η = σ₀ elsewhere | 15 | 1.9e-15 | 8.86e-05 | 4.2e-15 | 0.785 | 5.9e-15 | 1.5926e+03 | 2.0687e-01 |
| fixed on the centres, η = σ₀ elsewhere | 15 | 1.7e-15 | 1.02e-04 | 7.3e-15 | 0.673 | 8.9e-15 | 2.2666e+03 | 4.7051e-01 |
| one active face (atom 1's amplitude), η = ∞ | 9 | 1.1e-16 | 7.57e-05 | 1.4e-14 | 0.575 | 6.0e-15 | 1.9188e+03 | 3.4964e-01 |
| anchor gain K = ηP_amp | 18 | 6.8e-16 | 4.43e-05 | 5.7e-15 | 0.471 | 6.9e-15 | 1.8371e+03 | 2.9146e-01 |
| a random anchor gain, K + Kᵀ ⪰ 0, ‖K‖ ≈ η | 18 | 5.5e-16 | 1.47e-04 | 1.2e-15 | 0.495 | 5.8e-15 | 9.4406e+02 | 3.9894e-01 |

## C. The real loop: MR, DR and DR* driven by the same noise paths (k_e = 25, B = 100, std 0.01, η = ∞, RK4 dt = 0.001)

E‖f̂ − f‖²_ρ at the records (mean ± standard error over 12 runs) against the prediction (covariance part + mean part); MR and DR share the covariance part (their linearisations coincide).

| t | MR: Monte Carlo | MR: predicted | DR: Monte Carlo | DR: predicted | DR*: Monte Carlo | DR*: predicted | covariance part (MR, DR) | mean part: MR | DR | DR* |
|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 1.2049e-07 ± 3.2e-08 | 1.4519e-07 | 8.2906e-08 ± 2.0e-08 | 7.5793e-08 | 8.5242e-08 ± 2.1e-08 | 7.5144e-08 | 7.5790e-08 | 6.9404e-08 | 2.8142e-12 | 3.5996e-13 |
| 50 | 3.6510e-07 ± 8.9e-08 | 3.6891e-07 | 1.3140e-07 ± 3.7e-08 | 1.3184e-07 | 1.3303e-07 ± 3.6e-08 | 1.3086e-07 | 1.3183e-07 | 2.3708e-07 | 1.0366e-11 | 1.3097e-12 |
| 75 | 5.7835e-07 ± 7.4e-08 | 6.3186e-07 | 1.5430e-07 ± 3.0e-08 | 1.7319e-07 | 1.5968e-07 ± 3.3e-08 | 1.7203e-07 | 1.7317e-07 | 4.5869e-07 | 2.1615e-11 | 2.7001e-12 |
| 100 | 7.1395e-07 ± 1.3e-07 | 9.1005e-07 | 1.8124e-07 ± 3.0e-08 | 2.0428e-07 | 1.8119e-07 ± 3.1e-08 | 2.0303e-07 | 2.0424e-07 | 7.0581e-07 | 3.5819e-11 | 4.4282e-12 |
| 125 | 1.0566e-06 ± 1.5e-07 | 1.1868e-06 | 2.0191e-07 ± 3.2e-08 | 2.2626e-07 | 1.9823e-07 ± 3.4e-08 | 2.2504e-07 | 2.2620e-07 | 9.6060e-07 | 5.2439e-11 | 6.4221e-12 |
| 150 | 1.1344e-06 ± 1.6e-07 | 1.4552e-06 | 2.8096e-07 ± 8.1e-08 | 2.4312e-07 | 2.7769e-07 ± 8.2e-08 | 2.4184e-07 | 2.4305e-07 | 1.2121e-06 | 7.1081e-11 | 8.6315e-12 |
| 175 | 1.5203e-06 ± 1.8e-07 | 1.7118e-06 | 3.3507e-07 ± 9.4e-08 | 2.5790e-07 | 3.3152e-07 ± 9.4e-08 | 2.5653e-07 | 2.5780e-07 | 1.4540e-06 | 9.1452e-11 | 1.1021e-11 |
| 200 | 1.6688e-06 ± 2.1e-07 | 1.9532e-06 | 3.0488e-07 ± 5.3e-08 | 2.7057e-07 | 2.9782e-07 ± 5.2e-08 | 2.6917e-07 | 2.7046e-07 | 1.6828e-06 | 1.1333e-10 | 1.3564e-11 |
| 225 | 1.9434e-06 ± 2.1e-07 | 2.1774e-06 | 2.4968e-07 ± 4.6e-08 | 2.8069e-07 | 2.4401e-07 ± 4.4e-08 | 2.7928e-07 | 2.8056e-07 | 1.8968e-06 | 1.3653e-10 | 1.6243e-11 |
| 250 | 2.1220e-06 ± 2.0e-07 | 2.3841e-06 | 2.1120e-07 ± 4.5e-08 | 2.8862e-07 | 2.0796e-07 ± 4.4e-08 | 2.8717e-07 | 2.8846e-07 | 2.0956e-06 | 1.6093e-10 | 1.9042e-11 |
| 275 | 2.3755e-06 ± 2.2e-07 | 2.5743e-06 | 2.8778e-07 ± 6.1e-08 | 2.9497e-07 | 2.8120e-07 ± 5.7e-08 | 2.9348e-07 | 2.9479e-07 | 2.2795e-06 | 1.8640e-10 | 2.1951e-11 |
| 300 | 2.4310e-06 ± 2.4e-07 | 2.7491e-06 | 2.7989e-07 ± 5.5e-08 | 3.0023e-07 | 2.7468e-07 ± 5.3e-08 | 2.9869e-07 | 3.0001e-07 | 2.4490e-06 | 2.1285e-10 | 2.4960e-11 |

The run mean of f̂_MR − f̂_DR on the orbit's grid (common random numbers) against Ĵ(ȳ_MR − ȳ_DR): its mean square (with the run-to-run standard error of the mean function, in the same units), the prediction's, and the relative L²_ρ error of the prediction.

| t | ⟨(run mean)²⟩ | predicted | relative error ‖run mean − Ĵ(ȳ_MR − ȳ_DR)‖/‖Ĵ(ȳ_MR − ȳ_DR)‖ | sampling error of the run mean (rms) | relative error over the relative sampling error |
|---|---|---|---|---|---|
| 25 | 7.0616e-08 | 6.9660e-08 | 0.0527 | 1.02e-05 | 1.37 |
| 50 | 2.3844e-07 | 2.3790e-07 | 0.0283 | 1.15e-05 | 1.20 |
| 75 | 4.7282e-07 | 4.6016e-07 | 0.0337 | 1.36e-05 | 1.68 |
| 100 | 7.3366e-07 | 7.0792e-07 | 0.0297 | 1.38e-05 | 1.81 |
| 125 | 9.9049e-07 | 9.6327e-07 | 0.0215 | 1.47e-05 | 1.43 |
| 150 | 1.2483e-06 | 1.2153e-06 | 0.0229 | 1.38e-05 | 1.83 |
| 175 | 1.4628e-06 | 1.4575e-06 | 0.0152 | 1.41e-05 | 1.30 |
| 200 | 1.6906e-06 | 1.6865e-06 | 0.0146 | 1.33e-05 | 1.43 |
| 225 | 1.9022e-06 | 1.9007e-06 | 0.0123 | 1.42e-05 | 1.19 |
| 250 | 2.1028e-06 | 2.0997e-06 | 0.0052 | 1.78e-05 | 0.43 |
| 275 | 2.2983e-06 | 2.2836e-06 | 0.0070 | 1.42e-05 | 0.74 |
| 300 | 2.4552e-06 | 2.4532e-06 | 0.0079 | 1.36e-05 | 0.90 |

## D. The separations under the desired regressor (Corollary D‴.33)

### D1. At the lift

The separations' second-order drift S_k at the lift (real time, per unit Eν²; k_e = 25, B = 100, the noise's own quasi-static moments): its extreme eigenvalues per class and the noise std at which its largest eigenvalue reaches the separations' stiffness — σ₀/k_e under the adaptive prior, σ₀w_k/k_e under the cohesion leak (Proposition D‴.32, w_k = 1/3); '—' where every eigenvalue is negative.

| law | class 1: [min, max] | class 2 | class 3 | threshold std, adaptive prior | threshold std, cohesion leak |
|---|---|---|---|---|---|
| MR | [-0.37144, +0.01481] | [-0.02570, +0.01009] | [-0.03507, +0.29815] | 0.031 | 0.018 |
| DR | [-0.01122, -0.00050] | [-0.01490, -0.00015] | [-0.00719, -0.00034] | — | — |
| DR* | [+0.00016, +0.00378] | [+0.00005, +0.00498] | [+0.00011, +0.00238] | 0.238 | 0.138 |

### D2. Against k_e

The drift's extreme eigenvalues over the classes (λ_min the smallest, λ_max the largest) against k_e, with the noise's bandwidth scaled with the gain (B = 4k_e) and fixed (B = 100); (ii)'s bound: the largest ratio ‖S_k‖/((1/2k_e)‖∂²_xf‖_∞⟨‖M_k‖⟩) over the classes (M_k = G^{1/2}∇²_ϑΦ̃G^{1/2}; ‖∂²_xf‖_∞ = 5.0155 on the orbit, ⟨‖M_k‖⟩ = 0.4663, 0.5881, 0.3989), and the threshold std per law (adaptive prior). The bound gives the threshold a floor uniform in k_e, std ≥ (2σ₀/(‖∂²_xf‖_∞max_k⟨‖M_k‖⟩))^{1/2} = 0.0693, for DR and DR*.

| B | k_e | MR: λ_max | DR: λ_min, λ_max; k_e·λ_max | DR*: λ_min, λ_max; k_e·λ_max | bound ratio: DR; DR* | threshold std: MR; DR; DR* |
|---|---|---|---|---|---|---|
| 4k_e | 25 | +2.9815e-01 | -1.49e-02, -1.514e-04; -0.0038 | +5.11e-05, +4.979e-03; +0.1245 | 0.253; 0.084 | 0.031; —; 0.238 |
| 4k_e | 50 | +3.0140e-01 | -7.45e-03, -7.585e-05; -0.0038 | +2.54e-05, +2.485e-03; +0.1242 | 0.253; 0.084 | 0.022; —; 0.239 |
| 4k_e | 100 | +3.0305e-01 | -3.73e-03, -3.795e-05; -0.0038 | +1.27e-05, +1.242e-03; +0.1242 | 0.253; 0.084 | 0.015; —; 0.239 |
| 4k_e | 200 | +3.0387e-01 | -1.86e-03, -1.898e-05; -0.0038 | +6.34e-06, +6.208e-04; +0.1242 | 0.253; 0.084 | 0.011; —; 0.239 |
| 4k_e | 400 | +3.0429e-01 | -9.31e-04, -9.494e-06; -0.0038 | +3.17e-06, +3.104e-04; +0.1242 | 0.253; 0.084 | 0.008; —; 0.239 |
| 4k_e | 800 | +3.0450e-01 | -4.66e-04, -4.747e-06; -0.0038 | +1.58e-06, +1.552e-04; +0.1242 | 0.253; 0.084 | 0.005; —; 0.239 |
| 100 | 25 | +2.9815e-01 | -1.49e-02, -1.514e-04; -0.0038 | +5.11e-05, +4.979e-03; +0.1245 | 0.253; 0.084 | 0.031; —; 0.238 |
| 100 | 50 | +2.5211e-01 | -4.14e-03, -4.198e-05; -0.0021 | +4.24e-05, +4.141e-03; +0.2071 | 0.140; 0.140 | 0.024; —; 0.185 |
| 100 | 100 | +1.9045e-01 | -8.31e-06, +1.062e-05; +0.0011 | +3.17e-05, +3.104e-03; +0.3104 | 0.001; 0.211 | 0.019; 2.581; 0.151 |
| 100 | 200 | +1.2742e-01 | +1.06e-05, +1.035e-03; +0.2069 | +2.11e-05, +2.069e-03; +0.4139 | 0.140; 0.281 | 0.017; 0.185; 0.131 |
| 100 | 400 | +7.6592e-02 | +9.50e-06, +9.313e-04; +0.3725 | +1.27e-05, +1.242e-03; +0.4967 | 0.253; 0.337 | 0.015; 0.138; 0.119 |
| 100 | 800 | +4.2589e-02 | +6.16e-06, +6.036e-04; +0.4829 | +7.04e-06, +6.898e-04; +0.5519 | 0.327; 0.374 | 0.014; 0.121; 0.113 |

### D3. Other targets (DR at the lift)

Exactly represented targets with q = 1/K₀ (k_e = 25, B = 100): DR's drift per class (the centre entry, the (γ, γ) entry of S_k, and the extreme eigenvalues), DR*'s and MR's largest eigenvalues.

| target | bumps (A, C, S) | DR: centre entries | DR: [min, max] per class | DR*: λ_max | MR: λ_max |
|---|---|---|---|---|---|
| instance | (1, -1.5, 0.5); (-0.7, 0, 0.4); (0.8, 1.5, 0.5) | -1.12e-02, -1.49e-02, -7.18e-03 | [-1.12e-02, -4.96e-04]; [-1.49e-02, -1.51e-04]; [-7.19e-03, -3.36e-04] | +4.98e-03 | +2.981e-01 |
| two opposite bumps | (1, -1, 0.5); (-0.8, 1, 0.45) | -9.30e-03, -8.19e-03 | [-9.39e-03, -2.40e-04]; [-8.33e-03, -1.45e-04] | +3.21e-03 | +1.397e-02 |
| three same-sign bumps | (0.8, -1.5, 0.5); (0.8, 0, 0.5); (0.8, 1.5, 0.5) | -6.68e-03, -6.93e-03, -6.76e-03 | [-6.85e-03, +8.67e-06]; [-6.93e-03, +2.86e-04]; [-6.94e-03, +1.08e-05] | +2.35e-03 | +2.990e-01 |
| wide bumps (S = 1) | (1, -1.5, 1); (-0.7, 0, 1); (0.8, 1.5, 1) | -2.15e-03, -2.80e-03, -1.44e-03 | [-2.25e-03, +2.27e-04]; [-2.81e-03, +5.55e-04]; [-1.63e-03, +3.40e-04] | +9.34e-04 | +9.451e-02 |
| bumps at the orbit's edge | (1, -2.5, 0.5); (-0.7, 0, 0.4); (0.8, 2.5, 0.5) | -7.03e-04, -1.42e-02, -4.23e-04 | [-7.82e-04, +2.37e-04]; [-1.42e-02, -1.52e-04]; [-5.32e-04, +2.69e-04] | +4.76e-03 | +8.853e-02 |
| two close same-sign bumps | (0.8, -0.4, 0.5); (0.8, 0.4, 0.5) | -3.16e-03, -3.32e-03 | [-4.03e-03, +5.03e-04]; [-4.17e-03, +4.73e-04] | +1.44e-03 | +5.628e-02 |

### D4. At the stationary point

The exact periodic second-order mean of the K₀-atom loop (`second_order_mean(periodic=True)`, the forcing split into the noise part — the noise's own second moments — and the parameters' part), per unit Eν²: the mean error's rms (the lift's; the round-22 averaged noise-only balance (I − P)g/k_e; exact: total, noise part, parameters' part; the tangent part ‖PE[e]‖ and the orthogonal part), and the drift's [min, max] per class from the noise part of E[e], from its parameters' part, from the cross term, and in total.

| law | prior | rms E[e]: lift; averaged (I − P)g/k_e | exact: total; noise; parameters' | tangent; orthogonal | drift, noise part of E[e] | parameters' part of E[e] | cross term | total | threshold std (total) |
|---|---|---|---|---|---|---|---|---|---|
| DR | adaptive | 2.744e-02; 2.335e-04 | 6.697e-04; 2.106e-04; 6.321e-04 | 5.798e-04; 3.351e-04 | [-3.22e-06, +3.76e-06]; [-2.23e-07, +3.61e-07]; [-1.65e-06, +1.35e-06] | [-2.19e-05, +1.66e-05]; [-5.11e-06, +2.59e-04]; [-2.10e-05, +1.36e-05] | [-1.87e-05, +1.69e-05]; [-2.58e-04, +5.27e-06]; [-9.15e-06, +1.32e-05] | [-4.35e-05, +3.68e-05]; [-1.34e-06, +3.00e-06]; [-2.82e-05, +2.44e-05] | 2.773 |
| DR | fixed | 2.744e-02; 2.335e-04 | 1.127e-02; 1.120e-02; 1.469e-04 | 1.093e-02; 2.775e-03 | [-4.52e-03, -1.27e-04]; [-4.87e-03, -3.31e-05]; [-2.04e-03, -7.93e-05] | [-1.15e-05, +4.88e-07]; [-9.66e-05, -7.80e-07]; [-4.20e-06, +5.34e-08] | [-3.07e-05, +1.04e-05]; [-3.22e-04, +4.32e-06]; [-1.90e-05, +6.47e-06] | [-4.55e-03, -1.29e-04]; [-5.29e-03, -3.78e-05]; [-2.06e-03, -8.05e-05] | — |
| DR* | adaptive | 9.179e-03; 2.605e-04 | 6.438e-04; 1.366e-04; 6.273e-04 | 5.813e-04; 2.767e-04 | [-2.00e-06, +2.33e-06]; [-2.25e-07, +3.12e-07]; [-9.72e-07, +8.06e-07] | [-2.25e-05, +1.85e-05]; [-5.06e-06, +2.60e-04]; [-2.09e-05, +1.23e-05] | [-1.77e-05, +1.42e-05]; [-2.59e-04, +5.15e-06]; [-9.69e-06, +1.54e-05] | [-4.17e-05, +3.44e-05]; [-4.90e-07, +2.00e-06]; [-2.89e-05, +2.56e-05] | 2.868 |
| DR* | fixed | 9.179e-03; 2.605e-04 | 3.707e-03; 3.775e-03; 1.175e-04 | 3.589e-03; 9.286e-04 | [+4.24e-05, +1.53e-03]; [+1.18e-05, +1.65e-03]; [+2.63e-05, +6.60e-04] | [-1.04e-05, -2.12e-07]; [-9.56e-05, -8.35e-07]; [-4.88e-06, +6.01e-07] | [-2.97e-05, +8.48e-06]; [-3.23e-04, +4.23e-06]; [-1.95e-05, +7.81e-06] | [+4.06e-05, +1.50e-03]; [+6.20e-06, +1.23e-03]; [+2.50e-05, +6.46e-04] | 0.435 |

### D5. Along the second-order transient from the lift

DR: the exact second-order transient (`long_transient`, 82000 periods ≈ 5.15e+05 time units, 504 grid points per period), per unit Eν²: the period-averaged drift's extreme eigenvalues per class (all second-order terms: the mean error's noise and parameters' parts, the cross term) and the mean error's rms (total; noise; parameters'). The largest eigenvalue first positive at t = 4103, 2947, 2714 (classes 1, 2, 3; with the noise part of E[e] alone at t = 4072, 2721, 2689), its maximum +1.054e-04 (t ≈ 3.1e+04), +4.700e-05 (t ≈ 1.13e+04), +9.518e-05 (t ≈ 2.23e+04) (the noise part alone +1.242e-04 (t ≈ 3.42e+04), +4.248e-05 (t ≈ 1.15e+04), +1.071e-04 (t ≈ 2.48e+04)); at the stationary point (the period maps' limit) [-4.34e-05, +3.67e-05]; [-1.37e-06, +2.98e-06]; [-2.81e-05, +2.43e-05], rms E[e] 6.692e-04.

| t | class 1: [min, max] | class 2 | class 3 | rms E[e]: total; noise; parameters' |
|---|---|---|---|---|
| 0 | [-1.11e-02, -4.72e-04] | [-1.37e-02, -1.58e-04] | [-7.10e-03, -3.11e-04] | 2.644e-02; 2.644e-02; 7.404e-05 |
| 50.27 | [-1.07e-02, -4.59e-04] | [-1.44e-02, -1.60e-04] | [-6.77e-03, -3.23e-04] | 2.615e-02; 2.606e-02; 1.381e-04 |
| 301.6 | [-9.19e-03, -3.61e-04] | [-1.25e-02, -1.41e-04] | [-5.26e-03, -2.61e-04] | 2.253e-02; 2.235e-02; 2.545e-04 |
| 1263 | [-6.11e-03, -2.02e-04] | [-6.97e-03, -5.01e-05] | [-2.54e-03, -1.12e-04] | 1.420e-02; 1.408e-02; 2.065e-04 |
| 2670 | [-3.70e-03, -7.03e-05] | [-3.32e-03, -4.22e-06] | [-1.24e-03, -1.93e-06] | 9.651e-03; 9.590e-03; 1.550e-04 |
| 4128 | [-2.46e-03, +1.03e-06] | [-1.96e-03, +1.49e-05] | [-8.10e-04, +4.04e-05] | 8.221e-03; 8.170e-03; 1.320e-04 |
| 6183 | [-1.64e-03, +4.61e-05] | [-1.49e-03, +3.61e-05] | [-5.92e-04, +6.48e-05] | 7.479e-03; 7.431e-03; 1.148e-04 |
| 1.266e+04 | [-8.01e-04, +8.01e-05] | [-1.77e-03, +4.68e-05] | [-3.74e-04, +8.70e-05] | 6.441e-03; 6.420e-03; 9.594e-05 |
| 3.223e+04 | [-2.76e-04, +1.05e-04] | [-1.91e-03, +3.92e-05] | [-2.20e-04, +9.01e-05] | 4.977e-03; 5.077e-03; 1.508e-04 |
| 9.888e+04 | [-9.50e-05, +4.13e-05] | [-9.70e-04, +1.95e-05] | [-9.87e-05, +3.36e-05] | 2.156e-03; 2.516e-03; 3.826e-04 |
| 2.517e+05 | [-3.12e-05, +1.52e-05] | [-1.92e-04, +4.06e-06] | [-2.12e-05, +5.45e-06] | 2.572e-04; 5.408e-04; 5.819e-04 |
| 5.152e+05 | [-4.26e-05, +3.53e-05] | [-1.15e-05, +1.11e-06] | [-2.76e-05, +2.30e-05] | 6.377e-04; 2.114e-04; 6.289e-04 |

DR*: the exact second-order transient (`long_transient`, 82000 periods ≈ 5.15e+05 time units, 504 grid points per period), per unit Eν²: the period-averaged drift's extreme eigenvalues per class (all second-order terms: the mean error's noise and parameters' parts, the cross term) and the mean error's rms (total; noise; parameters'). The largest eigenvalue first positive at t = 0, 0, 0 (classes 1, 2, 3; with the noise part of E[e] alone at t = 0, 0, 0), its maximum +3.574e-03 (t ≈ 6.28), +3.843e-03 (t ≈ 25.1), +2.245e-03 (t ≈ 6.28) (the noise part alone +3.699e-03 (t ≈ 6.28), +4.669e-03 (t ≈ 6.28), +2.331e-03 (t ≈ 6.28)); at the stationary point (the period maps' limit) [-4.18e-05, +3.44e-05]; [-4.92e-07, +2.01e-06]; [-2.89e-05, +2.57e-05], rms E[e] 6.438e-04.

| t | class 1: [min, max] | class 2 | class 3 | rms E[e]: total; noise; parameters' |
|---|---|---|---|---|
| 0 | [+1.42e-04, +3.57e-03] | [+3.56e-05, +3.29e-03] | [+9.11e-05, +2.21e-03] | 8.719e-03; 8.729e-03; 1.943e-05 |
| 50.27 | [+1.44e-04, +3.46e-03] | [+3.77e-05, +3.84e-03] | [+1.03e-04, +2.15e-03] | 8.628e-03; 8.714e-03; 1.085e-04 |
| 301.6 | [+1.15e-04, +2.98e-03] | [+3.73e-05, +3.43e-03] | [+8.23e-05, +1.69e-03] | 7.283e-03; 7.461e-03; 2.377e-04 |
| 1263 | [+6.51e-05, +2.01e-03] | [+9.41e-06, +1.75e-03] | [+3.65e-05, +8.06e-04] | 4.613e-03; 4.725e-03; 1.862e-04 |
| 2670 | [+2.21e-05, +1.24e-03] | [-3.79e-06, +6.34e-04] | [-2.80e-07, +3.87e-04] | 3.196e-03; 3.255e-03; 1.275e-04 |
| 4128 | [-2.20e-06, +8.41e-04] | [-6.28e-06, +2.32e-04] | [-1.58e-05, +2.53e-04] | 2.733e-03; 2.781e-03; 9.795e-05 |
| 6183 | [-1.87e-05, +5.66e-04] | [-7.82e-06, +1.14e-04] | [-2.65e-05, +1.86e-04] | 2.482e-03; 2.528e-03; 7.245e-05 |
| 1.266e+04 | [-3.67e-05, +2.76e-04] | [-9.85e-06, +2.65e-04] | [-4.07e-05, +1.18e-04] | 2.161e-03; 2.180e-03; 3.680e-05 |
| 3.223e+04 | [-6.46e-05, +1.03e-04] | [-8.24e-06, +3.93e-04] | [-5.11e-05, +7.54e-05] | 1.825e-03; 1.721e-03; 1.241e-04 |
| 9.888e+04 | [-5.73e-05, +6.00e-05] | [-4.01e-06, +2.07e-04] | [-3.99e-05, +5.00e-05] | 1.225e-03; 8.603e-04; 3.736e-04 |
| 2.517e+05 | [-4.48e-05, +3.94e-05] | [-7.31e-07, +4.21e-05] | [-3.09e-05, +3.04e-05] | 7.585e-04; 2.188e-04; 5.765e-04 |
| 5.152e+05 | [-4.19e-05, +3.48e-05] | [-1.61e-07, +3.88e-06] | [-2.90e-05, +2.60e-05] | 6.509e-04; 1.372e-04; 6.241e-04 |

DR over the real loop's 300 time units (`second_order_transient`, dt = 0.0125): per 50-unit window the mean error's rms and the largest eigenvalue of the window-averaged drift per class, with the cross term and without it.

| window | rms E[e] | λ_max(S₁), λ_max(S₂), λ_max(S₃) | without the cross term |
|---|---|---|---|
| [0, 50) | 0.0265 | -4.635e-04, -1.614e-04, -3.248e-04 | -4.598e-04, -1.489e-04, -3.217e-04 |
| [50, 100) | 0.0259 | -4.533e-04, -1.575e-04, -3.189e-04 | -4.508e-04, -1.475e-04, -3.169e-04 |
| [100, 150) | 0.0249 | -4.272e-04, -1.542e-04, -3.057e-04 | -4.260e-04, -1.460e-04, -3.044e-04 |
| [150, 200) | 0.0242 | -3.986e-04, -1.546e-04, -2.929e-04 | -3.984e-04, -1.477e-04, -2.921e-04 |
| [200, 250) | 0.0236 | -3.836e-04, -1.506e-04, -2.806e-04 | -3.830e-04, -1.446e-04, -2.803e-04 |
| [250, 300) | 0.0230 | -3.711e-04, -1.454e-04, -2.689e-04 | -3.707e-04, -1.399e-04, -2.689e-04 |

DR* over [0, 50): the largest eigenvalue +3.6879e-03 with the cross term (+4.4464e-03 without it), the threshold std 0.277 (0.252; at the lift's quasi-static moments 0.238).

### D6. The real loop

The lift m = 15 with class k separated by 10⁻⁶·(2, 1, 0, −1, −2) along a start direction θ — u (MR's top eigenvector of S₃ at the lift) or the law's own top eigenvector of S_k ('top') —, 4 antithetic pairs per family on paths common to the laws, 300 time units, the separation's rate (least squares over t ≥ 50) against −σ₀/k_e = -2.8298e-04 plus Eν² times: the lift's θᵀS_kθ, and the transient's prediction — the least-squares slope of ∫₀ᵗθᵀ(S_k + C_k)θ + E‖P_⊥v₁(t)‖² (the mean error's term along the transient, the cross term, the variance term). The noise-induced part is the rate plus σ₀/k_e; its ratio to the prediction is given with its sampling error.

| law | std | class, θ | measured rate (pair means ± s.e.) | lift prediction | transient prediction (mean-error, cross, variance terms per unit Eν²) | noise part: measured / transient (± s.e.) | runs growing | runs slower than σ₀/k_e |
|---|---|---|---|---|---|---|---|---|
| MR | 0.05 | 3, u | +4.69100e-04 ± 1.3e-05 | +4.62398e-04 | +5.32002e-04 (+3.260e-01, -3.42e-05, +1.70e-05) | 0.923 ± 0.016 | 8 of 8 | 8 of 8 |
| MR | 0.1 | 3, u | +1.24844e-03 ± 3.1e-05 | +2.69852e-03 | +2.97693e-03 (+3.260e-01, -3.42e-05, +1.70e-05) | 0.470 ± 0.009 | 8 of 8 | 8 of 8 |
| DR | 0.05 | 3, u | -2.96505e-04 ± 1.3e-07 | -2.99526e-04 | -2.96816e-04 (-5.519e-03, -3.42e-05, +1.70e-05) | 0.978 ± 0.010 | 0 of 8 | 3 of 8 |
| DR | 0.1 | 3, u | -3.36497e-04 ± 5.3e-07 | -3.49179e-04 | -3.38337e-04 (-5.519e-03, -3.42e-05, +1.70e-05) | 0.967 ± 0.010 | 0 of 8 | 2 of 8 |
| DR | 0.05 | 1, top | -2.83987e-04 ± 2.3e-08 | -2.84217e-04 | -2.83965e-04 (-4.042e-04, -1.50e-06, +9.88e-06) | 1.022 ± 0.023 | 0 of 8 | 3 of 8 |
| DR | 0.1 | 1, top | -2.87007e-04 ± 9.1e-08 | -2.87940e-04 | -2.86934e-04 (-4.042e-04, -1.50e-06, +9.88e-06) | 1.019 ± 0.023 | 0 of 8 | 3 of 8 |
| DR | 0.05 | 2, top | -2.83360e-04 ± 2.5e-08 | -2.83354e-04 | -2.83311e-04 (-1.466e-04, -7.32e-06, +1.96e-05) | 1.146 ± 0.075 | 0 of 8 | 4 of 8 |
| DR | 0.1 | 2, top | -2.84498e-04 ± 1.0e-07 | -2.84490e-04 | -2.84319e-04 (-1.466e-04, -7.32e-06, +1.96e-05) | 1.133 ± 0.075 | 0 of 8 | 4 of 8 |
| DR | 0.05 | 3, top | -2.83698e-04 ± 2.2e-08 | -2.83816e-04 | -2.83699e-04 (-2.949e-04, -1.37e-06, +7.04e-06) | 0.999 ± 0.031 | 0 of 8 | 4 of 8 |
| DR | 0.1 | 3, top | -2.85856e-04 ± 9.2e-08 | -2.86337e-04 | -2.85868e-04 (-2.949e-04, -1.37e-06, +7.04e-06) | 0.996 ± 0.032 | 0 of 8 | 2 of 8 |
| DR* | 0.05 | 2, top | -2.74020e-04 ± 9.3e-07 | -2.70527e-04 | -2.73716e-04 (+4.120e-03, -4.48e-04, +3.12e-05) | 0.967 ± 0.100 | 0 of 8 | 4 of 8 |
| DR* | 0.1 | 2, top | -2.47317e-04 ± 3.7e-06 | -2.33182e-04 | -2.45939e-04 (+4.120e-03, -4.48e-04, +3.12e-05) | 0.963 ± 0.100 | 0 of 8 | 4 of 8 |

Numbers are produced by `experiments/exp_eiv_remedy.py`; nothing is edited by hand.
