# Proposition D‴.58: the transit under width mismatch (`exp_width_mismatch.py`)

Kernels of width σ = 0.12; a feature Aφ_σ*(· − c*), σ* = σ(1 + (τ/σ)²)^{1/2}; t = (σ* − σ)/(σ* + σ), λ = τ²/(σ² + σ*²) (`src/core/wide_feature.py`).

## A. The rules on ℝ

Relative L² errors: the half-squeeze rule (nodes τ(σ/σ*)^{1/2}z_k, amplitudes A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*)))) by the trapezoid rule and by the Mehler series (their relative difference), the series' first term's share of the squared error; moment matching (nodes τz_k, amplitudes A(σ*/σ)w_k); the half-squeeze nodes with least-squares amplitudes; the optimal K-kernel fit (variable projection; not computed below 1e-10, where its solve stops short of the optimum — A3 marks the optimum unconverged where the one-parameter rule beats it). Rows stop where the series falls below 1e-13 (the trapezoid rule's floor). Ratios per kernel e_{K−1}/e_K; the optimum's nodes over the half-squeeze nodes (range over the nodes).

### τ/σ = 0.25 (1/t = 65.9848, 1/λ = 33.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.1431e-02 | 2.1431e-02 (1e-15) | 1.00000 | 2.645e-02 | 2.1430e-02 | 2.1430e-02 | — | — |
| 2 | 3.7504e-04 | 3.7504e-04 (5e-15) | 0.99999 | 7.939e-04 | 3.7504e-04 | 3.7502e-04 | 57.1443, 57.1432, 33.3161 | 0.9999 … 0.9999 |
| 3 | 6.2262e-06 | 6.2262e-06 (1e-12) | 0.99998 | 2.399e-05 | 6.2262e-06 | 6.2260e-06 | 60.2354, 60.2352, 33.0854 | 1.0000 … 1.0000 |
| 4 | 1.0087e-07 | 1.0087e-07 (2e-10) | 0.99998 | 7.263e-07 | 1.0087e-07 | 1.0087e-07 | 61.7230, 61.7229, 33.0349 | 1.0000 … 1.0000 |
| 5 | 1.6114e-09 | 1.6114e-09 (4e-09) | 0.99997 | 2.200e-08 | 1.6114e-09 | 1.6114e-09 | 62.5986, 62.5985, 33.0177 | 1.0000 … 1.0000 |
| 6 | 2.5507e-11 | 2.5507e-11 (4e-07) | 0.99997 | 6.664e-10 | 2.5507e-11 | — | 63.1756, —, 33.0102 | — |
| 7 | 4.0116e-13 | 4.0115e-13 (3e-05) | 0.99997 | 2.019e-11 | 4.0116e-13 | — | 63.5825, —, 33.0064 | — |

### τ/σ = 0.5 (1/t = 17.9443, 1/λ = 9.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 7.8750e-02 | 7.8750e-02 (1e-15) | 0.99999 | 9.924e-02 | 7.8689e-02 | 7.8689e-02 | — | — |
| 2 | 5.0678e-03 | 5.0678e-03 (7e-16) | 0.99989 | 1.105e-02 | 5.0678e-03 | 5.0651e-03 | 15.5394, 15.5354, 8.9774 | 0.9990 … 0.9990 |
| 3 | 3.0939e-04 | 3.0939e-04 (9e-16) | 0.99978 | 1.232e-03 | 3.0939e-04 | 3.0925e-04 | 16.3798, 16.3790, 8.9698 | 0.9994 … 0.9994 |
| 4 | 1.8433e-05 | 1.8433e-05 (5e-13) | 0.99969 | 1.373e-04 | 1.8433e-05 | 1.8425e-05 | 16.7846, 16.7842, 8.9786 | 0.9996 … 0.9996 |
| 5 | 1.0828e-06 | 1.0828e-06 (1e-11) | 0.99962 | 1.528e-05 | 1.0828e-06 | 1.0824e-06 | 17.0229, 17.0227, 8.9849 | 0.9997 … 0.9997 |
| 6 | 6.3030e-08 | 6.3030e-08 (3e-10) | 0.99957 | 1.700e-06 | 6.3030e-08 | 6.3003e-08 | 17.1799, 17.1798, 8.9889 | 0.9997 … 0.9997 |
| 7 | 3.6452e-09 | 3.6452e-09 (4e-09) | 0.99953 | 1.890e-07 | 3.6452e-09 | 3.6437e-09 | 17.2912, 17.2911, 8.9916 | 0.9998 … 0.9998 |
| 8 | 2.0980e-10 | 2.0980e-10 (8e-08) | 0.99950 | 2.102e-08 | 2.0980e-10 | 2.0972e-10 | 17.3742, 17.3738, 8.9934 | 0.9998 … 0.9998 |
| 9 | 1.2031e-11 | 1.2031e-11 (1e-06) | 0.99948 | 2.337e-09 | 1.2031e-11 | — | 17.4385, —, 8.9947 | — |
| 10 | 6.8790e-13 | 6.8790e-13 (1e-05) | 0.99946 | 2.597e-10 | 6.8792e-13 | — | 17.4895, —, 8.9957 | — |

### τ/σ = 1 (1/t = 5.8284, 1/λ = 3.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.4090e-01 | 2.4090e-01 (1e-16) | 0.99946 | 3.237e-01 | 2.3915e-01 | 2.3915e-01 | — | — |
| 2 | 4.7745e-02 | 4.7745e-02 (3e-16) | 0.99867 | 1.112e-01 | 4.7723e-02 | 4.7497e-02 | 5.0456, 5.0350, 2.9107 | 0.9901 … 0.9901 |
| 3 | 8.9780e-03 | 8.9780e-03 (1e-15) | 0.99768 | 3.757e-02 | 8.9771e-03 | 8.9365e-03 | 5.3180, 5.3149, 2.9608 | 0.9939 … 0.9939 |
| 4 | 1.6474e-03 | 1.6474e-03 (4e-15) | 0.99691 | 1.261e-02 | 1.6474e-03 | 1.6402e-03 | 5.4499, 5.4485, 2.9793 | 0.9956 … 0.9957 |
| 5 | 2.9802e-04 | 2.9802e-04 (5e-14) | 0.99634 | 4.221e-03 | 2.9802e-04 | 2.9676e-04 | 5.5277, 5.5269, 2.9873 | 0.9966 … 0.9967 |
| 6 | 5.3418e-05 | 5.3418e-05 (2e-13) | 0.99590 | 1.411e-03 | 5.3418e-05 | 5.3197e-05 | 5.5791, 5.5785, 2.9915 | 0.9972 … 0.9973 |
| 7 | 9.5126e-06 | 9.5126e-06 (1e-12) | 0.99555 | 4.713e-04 | 9.5126e-06 | 9.4739e-06 | 5.6154, 5.6151, 2.9939 | 0.9976 … 0.9977 |
| 8 | 1.6859e-06 | 1.6859e-06 (6e-12) | 0.99527 | 1.573e-04 | 1.6859e-06 | 1.6791e-06 | 5.6426, 5.6423, 2.9954 | 0.9979 … 0.9980 |
| 9 | 2.9767e-07 | 2.9767e-07 (5e-11) | 0.99504 | 5.251e-05 | 2.9767e-07 | 2.9649e-07 | 5.6636, 5.6633, 2.9964 | 0.9982 … 0.9983 |
| 10 | 5.2404e-08 | 5.2404e-08 (1e-10) | 0.99484 | 1.752e-05 | 5.2404e-08 | 5.2197e-08 | 5.6803, 5.6801, 2.9971 | 0.9984 … 0.9984 |

### τ/σ = 2 (1/t = 2.6180, 1/λ = 1.5000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 5.2280e-01 | 5.2280e-01 (0e+00) | 0.98665 | 8.088e-01 | 5.0462e-01 | 5.0462e-01 | — | — |
| 2 | 2.3059e-01 | 2.3059e-01 (2e-16) | 0.98657 | 5.371e-01 | 2.2815e-01 | 2.2364e-01 | 2.2672, 2.2564, 1.5059 | 0.9485 … 0.9485 |
| 3 | 9.6611e-02 | 9.6611e-02 (4e-16) | 0.98401 | 3.600e-01 | 9.6247e-02 | 9.4014e-02 | 2.3868, 2.3787, 1.4920 | 0.9657 … 0.9657 |
| 4 | 3.9497e-02 | 3.9497e-02 (7e-16) | 0.98168 | 2.408e-01 | 3.9438e-02 | 3.8509e-02 | 2.4460, 2.4414, 1.4946 | 0.9746 … 0.9770 |
| 5 | 1.5918e-02 | 1.5918e-02 (9e-16) | 0.97978 | 1.609e-01 | 1.5907e-02 | 1.5539e-02 | 2.4813, 2.4783, 1.4966 | 0.9799 … 0.9818 |
| 6 | 6.3555e-03 | 6.3555e-03 (2e-16) | 0.97824 | 1.074e-01 | 6.3532e-03 | 6.2092e-03 | 2.5046, 2.5025, 1.4977 | 0.9834 … 0.9854 |
| 7 | 2.5208e-03 | 2.5208e-03 (4e-15) | 0.97697 | 7.171e-02 | 2.5203e-03 | 2.4644e-03 | 2.5212, 2.5196, 1.4983 | 0.9858 … 0.9876 |
| 8 | 9.9499e-04 | 9.9499e-04 (5e-15) | 0.97591 | 4.785e-02 | 9.9488e-04 | 9.7318e-04 | 2.5335, 2.5323, 1.4987 | 0.9877 … 0.9894 |
| 9 | 3.9125e-04 | 3.9125e-04 (8e-16) | 0.97502 | 3.192e-02 | 3.9122e-04 | 3.8282e-04 | 2.5431, 2.5421, 1.4990 | 0.9891 … 0.9906 |
| 10 | 1.5339e-04 | 1.5339e-04 (1e-14) | 0.97426 | 2.129e-02 | 1.5338e-04 | 1.5013e-04 | 2.5507, 2.5499, 1.4992 | 0.9902 … 0.9917 |

### A2. The series' tail

D_{2K+2j,K}/D_{2K,K} (integers) for j = 1, …, 4 and K = 1, …, 30; their j-th differences in K, constant if they are polynomials of degree j, against j!·binom(1/2, j)(−2)^j = −1, −1, −3, −15.

| j | K = 1, 2, 3, 4, 5, …, 30 | j-th differences (all K) | j!·binom(1/2, j)(−2)^j |
|---|---|---|---|
| 1 | 0, -1, -2, -3, -4, …, -29 | -1 | -1 |
| 2 | 15, 26, 36, 45, 53, …, -72 | -1 | -1 |
| 3 | 0, 86, 270, 549, 920, …, 31320 | -3 | -3 |
| 4 | 945, 751, 513, 1206, 3790, …, 3183840 | -15 | -15 |

### A3. Scaled rules

The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`) at the half-squeeze variance v_hs = τ²σ/σ*, the natural Fock scaling v_hs/cosh r and the Mercer-zero scaling τ²(1 + 2τ²/σ²)^{−1/2}: their ratios per kernel at the last K against 1/t, (2 − λ)/λ = (1 − t + t²)/t and 1/q, q = (β² − 1)/(β² + 1), β² = (1 + 2τ²/σ²)^{1/2}; the variance of least error (`best_scale`): v/v_hs, its error and the optimal fit's over the rule's, at K = 1 and at the last K with an optimum (in brackets the share of the gap to the optimum that the variance recovers, (e_hs − e_best)/(e_hs − min(e_opt, e_best)); 'optimum unconverged' where e_best < e_opt by more than a relative 1e-09); over all K with an optimum the shares 93.0%–100.0% (1 optima unconverged):

| τ/σ | half-squeeze (1/t) | Fock ((2 − λ)/λ) | Mercer-zero (1/q) | K = 1: v/v_hs, e_best/e_hs, e_opt/e_hs | the last K: K, v/v_hs, e_best/e_hs, e_opt/e_hs |
|---|---|---|---|---|---|
| 0.25 | 63.583 (65.985) | 63.500 (65.000) | 34.018 (33.971) | 0.99954, 0.99994, 0.99994 (100.0%) | 5, 0.99995, 0.99997, 0.99997 (100.0%) |
| 0.5 | 17.489 (17.944) | 17.079 (17.000) | 9.916 (9.899) | 0.99381, 0.99922, 0.99922 (100.0%) | 8, 0.99959, 0.99959, 0.99960 (100.0%, optimum unconverged) |
| 1 | 5.680 (5.828) | 5.034 (5.000) | 3.751 (3.732) | 0.94281, 0.99272, 0.99272 (100.0%) | 10, 0.99678, 0.99608, 0.99605 (99.1%) |
| 2 | 2.551 (2.618) | 2.021 (2.000) | 2.021 (2.000) | 0.74536, 0.96523, 0.96523 (100.0%) | 10, 0.98131, 0.97968, 0.97877 (95.7%) |

### A4. Narrower features: the continued frame

σ* < σ: the nodes imaginary in conjugate pairs, the K-term sum e^{−x²/(2σ²)}Σ_k b_k cos(η_kx/σ²) at the Hermite-zero frequencies η_k = (2σ²|sinh r|)^{1/2}z_k (`narrow_rule`); its relative L² error (trapezoid) against the Mehler series at |t| = (σ − σ*)/(σ + σ*), and the ratios per kernel against 1/|t|:

| σ*/σ | 1/|t| | K: error (relative difference from the series) | ratio per kernel at the last K |
|---|---|---|---|
| 0.5 | 3.0000 | 1: 4.595e-01 (2e-16), 2: 1.769e-01 (0e+00), 3: 6.469e-02 (7e-16), 4: 2.308e-02 (0e+00), 5: 8.115e-03 (2e-16), 6: 2.827e-03 (2e-15), 7: 9.786e-04 (2e-15), 8: 3.370e-04 (2e-15), 9: 1.156e-04 (1e-14), 10: 3.956e-05 (3e-14) | 2.9231 |
| 0.8 | 9.0000 | 1: 1.567e-01 (0e+00), 2: 2.010e-02 (2e-16), 3: 2.447e-03 (1e-15), 4: 2.908e-04 (2e-15), 5: 3.406e-05 (9e-15), 6: 3.953e-06 (6e-13), 7: 4.558e-07 (6e-13), 8: 5.231e-08 (3e-11), 9: 5.981e-09 (7e-11), 10: 6.819e-10 (1e-10) | 8.7717 |

### A5. The frames: which Gauss–Hermite rule a Mehler expansion singles out

The Mehler frames of the kernel e^{−(y−y′)²/(4σ²)} (ratio ρ, Hermite scale ℓ² = 2σ²ρ/(1 − ρ²)): the Gauss–Hermite rule of the frame's own Gaussian N(0, v_ρ), v_ρ = 2σ²τ²/(2σ² + (1 − ρ)τ²), has no defect below 2K in the frame; v′ = v_ρ/ℓ² the measure's variance in the frame's units. At ρ = t (the half-squeeze rule), q (the Mercer expansion) and λ: v_ρ over the half-squeeze variance, and v′; the Mercer frame's rule N(0, 2τ²/(β² + 1)) — not the Mercer-zero rule — and its ratio per kernel at the last three K against 1/t and 1/q:

| τ/σ | ρ = t: v/v_hs, v′ | ρ = q: v/v_hs, v′ | ρ = λ: v/v_hs, v′ | the Mercer frame's rule: ratio per kernel (1/t, 1/q) |
|---|---|---|---|---|
| 0.25 | 1.000000, 2.000000 | 1.000433, 1.029437 | 1.000459, 1.000000 | 62.5377, 63.1014, 63.4950 at K = 5–7 (65.98, 33.97) |
| 0.5 | 1.000000, 2.000000 | 1.005090, 1.101021 | 1.006231, 1.000000 | 17.1134, 17.1485, 17.1718 at K = 8–10 (17.94, 9.899) |
| 1 | 1.000000, 2.000000 | 1.035276, 1.267949 | 1.060660, 1.000000 | 5.3221, 5.3221, 5.3210 at K = 8–10 (5.828, 3.732) |
| 2 | 1.000000, 2.000000 | 1.118034, 1.500000 | 1.341641, 1.000000 | 2.3304, 2.3314, 2.3320 at K = 8–10 (2.618, 2) |

The rules' error series in their frames, exactly (`frame_series`: 110 digits, K = 1 … 40; referee round 81) — the Mercer frame's rule (ρ = q, v′ = 1 + q) and the half-squeeze rule (ρ = t, v′ = 2): the tail after the first term over the first at K = 5, 10, 20, 40 and the first K from which it exceeds it; the ratio per kernel (S(K − 1)/S(K))^{1/2} at K = 10, 20, 30, 40 against 1/t and 1/q; and the series' ratio per kernel against the actual errors' where both are tabulated (the largest relative difference):

| τ/σ | frame | tail/first at K = 5, 10, 20, 40 | the tail above the first from K | ratio per kernel at K = 10, 20, 30, 40 (1/t, 1/q) | against the actual errors |
|---|---|---|---|---|---|
| 0.25 | Mercer | 0.00765, 0.0258, 0.0944, 0.38 | — | 64.1872, 64.9029, 65.0696, 65.1135 (65.98, 33.97) | 3.2e-05 (K ≤ 7) |
| 0.25 | half-squeeze | 2.78e-05, 4.03e-05, 4.81e-05, 5.26e-05 | — | 64.3140, 65.1548, 65.4327, 65.5711 (65.98, 33.97) | 3.1e-05 (K ≤ 7) |
| 0.5 | Mercer | 0.0808, 0.283, 1.21, 8.6 | 19 | 17.1721, 17.2032, 17.1837, 17.1741 (17.94, 9.899) | 2.0e-05 (K ≤ 10) |
| 0.5 | half-squeeze | 0.000378, 0.000545, 0.000651, 0.000711 | — | 17.4897, 17.7185, 17.7941, 17.8318 (17.94, 9.899) | 1.3e-05 (K ≤ 10) |
| 1 | Mercer | 0.448, 1.86, 14.7, 687 | 8 | 5.3210, 5.3115, 5.3100, 5.3097 (5.828, 3.732) | 3.2e-10 (K ≤ 10) |
| 1 | half-squeeze | 0.00367, 0.00518, 0.00617, 0.00674 | — | 5.6803, 5.7549, 5.7796, 5.7918 (5.828, 3.732) | 1.9e-10 (K ≤ 10) |
| 2 | Mercer | 0.967, 4.44, 55.3, 8.35e+03 | 6 | 2.3320, 2.3330, 2.3332, 2.3333 (2.618, 2) | 6.0e-14 (K ≤ 10) |
| 2 | half-squeeze | 0.0206, 0.0264, 0.0307, 0.0335 | — | 2.5507, 2.5848, 2.5960, 2.6015 (2.618, 2) | 1.4e-14 (K ≤ 10) |

The defects below 2K, zero exactly, are at most 1.5e-52 of their parts (the rounding of the 110-digit arithmetic where the rule is exact; it bounds that rounding, not the tail's).

The Mercer frame's rule further out (`a5_far`): its ratio per kernel by K — τ/σ = 2 (140 digits): 60: 2.333302, 80: 2.333316, 100: 2.333322, 120: 2.333326; the successive differences' ratios 0.450, 0.534 (K⁻²'s 0.463, 0.543); the limit by L + c/K² + d/K³ through the last three 2.333333, the formula's (b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{1/2}, 2.333333 (A6), 1/t 2.6180; the last K with 40 more digits the same in double; τ/σ = 1 (140 digits): 60: 5.309511, 80: 5.309459, 100: 5.309437, 120: 5.309425; the successive differences' ratios 0.425, 0.516 (K⁻²'s 0.463, 0.543); the limit by L + c/K² + d/K³ through the last three 5.309401, the formula's (b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{1/2}, 5.309401 (A6), 1/t 5.8284; the last K with 40 more digits the same in double; τ/σ = 0.5 (140 digits): 60: 17.168222, 80: 17.166623, 100: 17.165974, 120: 17.165645; the successive differences' ratios 0.406, 0.508 (K⁻²'s 0.463, 0.543); the limit by L + c/K² + d/K³ through the last three 17.164987, the formula's (b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{1/2}, 17.164966 (A6), 1/t 17.9443; the last K with 40 more digits the same in double; τ/σ = 0.25 (220, 300, 380, 460, 540, 700 digits): 40: 65.113507, 60: 65.111816, 80: 65.092630, 100: 65.078756, 120: 65.070472, 160: 65.062764; the successive differences' ratios 11.348, 0.723, 0.597, 0.930 (K⁻²'s 0.350, 0.463, 0.543, 0.994); the limit by L + c/K² + d/K³ through the last three 65.054019, the formula's (b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{1/2}, 65.055507 (A6), 1/t 65.9848; the last K with 40 more digits the same in double — the low defects at most 3.7e-22. Its limit lies below 1/t at every τ/σ, within 2.3e-05 of the formula's — D‴.61 (iv)'s limit there (the frame's t(v) at most 0.4286, below its branch and bound's r* = 0.8052 already, and with VA (iv) holds at every frame; by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's Appendix A.0 assumed correct, now independent cross-checks; referee rounds 96, 101, 107).

The first tail ratio D_{2K+2}/D_{2K} = −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1) (`gh_defect_v`, exact) at v′ = 1, 3/2, 7/4, 2, 3, K = 1 … 12: exact in all 60 — −(K − 1) at v′ = 2, growing like (v′ − 2)K² at any other v′.


### A6. The scaled rules' rate: one formula (referee round 84; proved by D‴.61 (iv) where a rule has a Mehler frame (by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's Appendix A.0 assumed correct, now independent cross-checks); by D‴.62 (b) where it has none, without computer assistance (referee round 132))

The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`), its relative L² error computed exactly (`a6_error`: the closed-form Gaussian inner products in decimal arithmetic of 2K·log₁₀(1/t) + 60 digits, σ = 1), against the formula t(v) = max{h(τ²/v), h(vσ*²/(σ²τ²))}, h(y) = (y − 1)/(y + 1) — the two arguments are (σ*/σ)(v_hs/v) and (σ*/σ)(v/v_hs), v_hs = τ²σ/σ*, so t(v) = tanh((r + |ln(v/v_hs)|)/2), r = ln(σ*/σ): a mis-scaling of the variance by e^u adds |u| to the squeeze, and t(v) ≥ t with equality only at v_hs, a cusp —: the ratio per kernel e_{K−1}/e_K at K = 20, 40, 60, 80, 120, 160; the successive differences' ratios (as K⁻²: 0.185, 0.350, 0.714, 0.350; as K⁻¹: 0.333, 0.500, 1.000, 0.500); the limit by L + c/K² + d/K³ through the last three K (at v_hs by L + c/K + d/K²) against 1/t(v):

| τ/σ | rule | v/v_hs | 1/t(v) (attained at) | e_{K−1}/e_K at K = 20, 40, 60, 80, 120, 160 | the differences' ratios | the limit | over 1/t(v), minus 1 | e_{K−1}/e_K·t(v) − 1 at K = 160 (its exponent from K = 120) |
|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | 1.030776 | 33.000000 (vσ*²/τ²) | 33.000309, 33.000047, 33.000016, 33.000008, 33.000003, 33.000002 | 0.116, 0.275, 0.591, 0.293 | 33.000000 | +5.3e-11 | +4.6e-08 (2.34) |
| 0.25 | Fock | 0.999541 | 65.000000 (τ²/v) | 64.886538, 65.074448, 65.061894, 65.038777, 65.015220, 65.007463 | -0.067, 1.842, 1.019, 0.329 | 64.998108 | -2.9e-05 | +1.1e-04 (2.48) |
| 0.25 | Mercer-zero | 0.971825 | 33.970563 (τ²/v) | 33.976419, 33.972041, 33.971222, 33.970934, 33.970728, 33.970656 | 0.187, 0.351, 0.716, 0.351 | 33.970563 | +3.0e-10 | +2.7e-06 (2.00) |
| 0.25 | Mercer frame | 1.000433 | 65.055507 (vσ*²/τ²) | 64.902928, 65.113507, 65.111816, 65.092630, 65.070472, 65.062764 | -0.008, 11.348, 1.155, 0.348 | 65.052913 | -4.0e-05 | +1.1e-04 (2.52) |
| 0.25 | 0.5·v_hs | 0.500000 | 2.884033 (τ²/v) | 2.900260, 2.886552, 2.885089, 2.884651, 2.884324, 2.884202 | 0.107, 0.299, 0.746, 0.373 | 2.884035 | +9.3e-07 | +5.9e-05 (1.89) |
| 0.25 | 0.8·v_hs | 0.800000 | 7.933118 (τ²/v) | 7.936888, 7.934118, 7.933574, 7.933378, 7.933235, 7.933184 | 0.196, 0.361, 0.729, 0.356 | 7.933118 | +7.0e-09 | +8.4e-06 (1.97) |
| 0.25 | 1·v_hs | 1.000000 | 65.984845 (both arguments, the cusp) | 65.154799, 65.571139, 65.709331, 65.778318, 65.847233, 65.881663 | 0.332, 0.499, 0.999, 0.500 | 65.984845 | -5.0e-09 | -1.6e-03 (1.00) |
| 0.25 | 1.25·v_hs | 1.250000 | 7.933118 (vσ*²/τ²) | 7.931683, 7.932747, 7.932953, 7.933025, 7.933077, 7.933095 | 0.193, 0.350, 0.714, 0.350 | 7.933118 | +5.8e-11 | -2.9e-06 (2.00) |
| 0.25 | 2·v_hs | 2.000000 | 2.884033 (vσ*²/τ²) | 2.892369, 2.884232, 2.883972, 2.883991, 2.884014, 2.884022 | 0.032, -0.073, 1.208, 0.358 | 2.884033 | +8.7e-08 | -3.6e-06 (2.00) |
| 0.5 | moment matching | 1.118034 | 9.000000 (vσ*²/τ²) | 8.998866, 8.999712, 8.999871, 8.999927, 8.999968, 8.999982 | 0.189, 0.352, 0.717, 0.351 | 9.000000 | +1.9e-10 | -2.0e-06 (2.00) |
| 0.5 | Fock | 0.993808 | 17.000000 (τ²/v) | 17.051366, 17.010770, 17.004247, 17.002287, 17.000981, 17.000543 | 0.161, 0.300, 0.666, 0.335 | 17.000005 | +3.0e-07 | +3.2e-05 (2.05) |
| 0.5 | Mercer-zero | 0.912871 | 9.898979 (τ²/v) | 9.903568, 9.900195, 9.899531, 9.899294, 9.899121, 9.899059 | 0.197, 0.359, 0.726, 0.355 | 9.898980 | +5.7e-09 | +8.1e-06 (1.98) |
| 0.5 | Mercer frame | 1.005090 | 17.164966 (vσ*²/τ²) | 17.203179, 17.174148, 17.168222, 17.166623, 17.165645, 17.165334 | 0.204, 0.270, 0.612, 0.318 | 17.164975 | +5.4e-07 | +2.1e-05 (2.13) |
| 0.5 | 0.5·v_hs | 0.500000 | 2.618034 (τ²/v) | 2.625524, 2.620437, 2.619227, 2.618749, 2.618375, 2.618234 | 0.238, 0.394, 0.782, 0.379 | 2.618038 | +1.6e-06 | +7.6e-05 (1.86) |
| 0.5 | 0.8·v_hs | 0.800000 | 6.030909 (τ²/v) | 6.035311, 6.032123, 6.031469, 6.031230, 6.031054, 6.030992 | 0.205, 0.365, 0.735, 0.358 | 6.030909 | +2.4e-08 | +1.4e-05 (1.96) |
| 0.5 | 1·v_hs | 1.000000 | 17.944272 (both arguments, the cusp) | 17.718495, 17.831753, 17.869341, 17.888104, 17.906847, 17.916211 | 0.332, 0.499, 0.999, 0.500 | 17.944272 | -2.9e-09 | -1.6e-03 (1.00) |
| 0.5 | 1.25·v_hs | 1.250000 | 6.030909 (vσ*²/τ²) | 6.029699, 6.030608, 6.030775, 6.030834, 6.030876, 6.030890 | 0.184, 0.349, 0.713, 0.350 | 6.030909 | +8.5e-11 | -3.1e-06 (2.00) |
| 0.5 | 2·v_hs | 2.000000 | 2.618034 (vσ*²/τ²) | 2.617442, 2.617887, 2.617969, 2.617997, 2.618018, 2.618025 | 0.184, 0.348, 0.712, 0.349 | 2.618034 | +3.3e-11 | -3.5e-06 (2.00) |
| 1 | moment matching | 1.414214 | 3.000000 (vσ*²/τ²) | 2.999321, 2.999835, 2.999927, 2.999959, 2.999982, 2.999990 | 0.180, 0.347, 0.710, 0.348 | 3.000000 | -7.5e-11 | -3.4e-06 (2.01) |
| 1 | Fock | 0.942809 | 5.000000 (τ²/v) | 5.008441, 5.002313, 5.001079, 5.000624, 5.000286, 5.000164 | 0.201, 0.368, 0.743, 0.363 | 5.000001 | +2.1e-07 | +3.3e-05 (1.94) |
| 1 | Mercer-zero | 0.816497 | 3.732051 (τ²/v) | 3.738119, 3.733875, 3.732926, 3.732565, 3.732290, 3.732189 | 0.224, 0.381, 0.760, 0.369 | 3.732052 | +3.2e-07 | +3.7e-05 (1.91) |
| 1 | Mercer frame | 1.035276 | 5.309401 (vσ*²/τ²) | 5.311526, 5.309690, 5.309511, 5.309459, 5.309425, 5.309414 | 0.098, 0.292, 0.644, 0.324 | 5.309401 | +2.5e-08 | +2.5e-06 (2.11) |
| 1 | 0.5·v_hs | 0.500000 | 2.093836 (τ²/v) | 2.102939, 2.097054, 2.095519, 2.094881, 2.094358, 2.094150 | 0.261, 0.416, 0.819, 0.397 | 2.093852 | +7.5e-06 | +1.5e-04 (1.76) |
| 1 | 0.8·v_hs | 0.800000 | 3.604957 (τ²/v) | 3.611070, 3.606807, 3.605847, 3.605481, 3.605201, 3.605098 | 0.225, 0.382, 0.762, 0.370 | 3.604958 | +3.7e-07 | +3.9e-05 (1.91) |
| 1 | 1·v_hs | 1.000000 | 5.828427 (both arguments, the cusp) | 5.754948, 5.791841, 5.804071, 5.810173, 5.816267, 5.819310 | 0.332, 0.499, 0.999, 0.499 | 5.828427 | +2.2e-08 | -1.6e-03 (1.00) |
| 1 | 1.25·v_hs | 1.250000 | 3.604957 (vσ*²/τ²) | 3.604193, 3.604771, 3.604875, 3.604911, 3.604937, 3.604946 | 0.180, 0.346, 0.710, 0.348 | 3.604957 | -2.5e-11 | -3.2e-06 (2.01) |
| 1 | 2·v_hs | 2.000000 | 2.093836 (vσ*²/τ²) | 2.093425, 2.093735, 2.093792, 2.093811, 2.093825, 2.093830 | 0.181, 0.347, 0.711, 0.349 | 2.093836 | +1.8e-11 | -3.0e-06 (2.01) |
| 2 | moment matching | 2.236068 | 1.500000 (vσ*²/τ²) | 1.499815, 1.499955, 1.499980, 1.499989, 1.499995, 1.499997 | 0.180, 0.346, 0.709, 0.348 | 1.500000 | -1.7e-11 | -1.8e-06 (2.01) |
| 2 | Fock | 0.745356 | 2.000000 (τ²/v) | 2.008529, 2.003233, 2.001764, 2.001127, 2.000584, 2.000360 | 0.277, 0.434, 0.852, 0.413 | 2.000028 | +1.4e-05 | +1.8e-04 (1.68) |
| 2 | Mercer-zero | 0.745356 | 2.000000 (τ²/v) | 2.008529, 2.003233, 2.001764, 2.001127, 2.000584, 2.000360 | 0.277, 0.434, 0.852, 0.413 | 2.000028 | +1.4e-05 | +1.8e-04 (1.68) |
| 2 | Mercer frame | 1.118034 | 2.333333 (vσ*²/τ²) | 2.333009, 2.333260, 2.333302, 2.333316, 2.333326, 2.333329 | 0.169, 0.331, 0.690, 0.341 | 2.333333 | -2.1e-09 | -1.8e-06 (2.04) |
| 2 | 0.5·v_hs | 0.500000 | 1.576014 (τ²/v) | 1.586907, 1.580602, 1.578680, 1.577797, 1.577000, 1.576650 | 0.305, 0.460, 0.902, 0.439 | 1.576110 | +6.1e-05 | +4.0e-04 (1.52) |
| 2 | 0.8·v_hs | 0.800000 | 2.114153 (τ²/v) | 2.122458, 2.117220, 2.115809, 2.115203, 2.114692, 2.114483 | 0.270, 0.429, 0.844, 0.409 | 2.114177 | +1.1e-05 | +1.6e-04 (1.70) |
| 2 | 1·v_hs | 1.000000 | 2.618034 (both arguments, the cusp) | 2.584772, 2.601523, 2.607057, 2.609813, 2.612562, 2.613933 | 0.330, 0.498, 0.997, 0.499 | 2.618035 | +3.1e-07 | -1.6e-03 (1.00) |
| 2 | 1.25·v_hs | 1.250000 | 2.114153 (vσ*²/τ²) | 2.113761, 2.114064, 2.114115, 2.114132, 2.114144, 2.114148 | 0.168, 0.338, 0.699, 0.344 | 2.114153 | -1.3e-09 | -2.5e-06 (2.02) |
| 2 | 2·v_hs | 2.000000 | 1.576014 (vσ*²/τ²) | 1.575799, 1.575962, 1.575991, 1.576001, 1.576009, 1.576011 | 0.179, 0.346, 0.708, 0.348 | 1.576014 | -6.9e-11 | -2.0e-06 (2.01) |

Off the cusp the limit lies within 6.1e-05 of the formula's 1/t(v) at all 32 rules (31 distinct) (the farthest: τ/σ = 2, 0.5·v_hs); at v_hs within 3.2e-07 of 1/t. With a K⁻¹ term instead (L + b/K + c/K² through the same K) off the cusp within 1.4e-04 (the farthest: τ/σ = 2, 0.5·v_hs): the second fit's largest deviation; the two fits differ by at most 2.0e-04 (bounds rounded up: referee round 96). Off the cusp, the 18 rules (17 distinct) with a frame, whose limits D‴.61 (iv) proves (the fits' deviations their bias): the first fit within 6.1e-05, the second within 1.4e-04, the two bracketing 1/t(v) at 10 of the 18. Off the cusp, the 14 rules without a frame, whose limits D‴.62 (b) proves by hand (the fits' deviations their bias too; 'the proved rules' below are the frame rules, the calibration predating D‴.62): the first fit within 7.6e-06, the second within 4.5e-05, the two bracketing 1/t(v) at 9 of the 14. By frame, the branch attaining t(v) and nearness to the cusp (far: |ln(v/v_hs)| ≥ 0.05; referee round 96): with a frame, τ²/v, far (8, 7 distinct): within 6.1e-05 and 1.4e-04, bracketing at 8; with a frame, τ²/v, near (3): within 3.0e-05 and 5.4e-06, bracketing at 2; with a frame, vσ*²/τ², far (4): within 2.1e-09 and 1.8e-07, bracketing at 0; with a frame, vσ*²/τ², near (3): within 4.0e-05 and 3.7e-05, bracketing at 0; without, τ²/v, far (5): within 7.6e-06 and 4.5e-05, bracketing at 5; without, vσ*²/τ², far (8): within 8.8e-08 and 3.6e-07, bracketing at 4; without, vσ*²/τ², near (1): within 5.3e-11 and 3.8e-08, bracketing at 0 — the fits' bracket calibrated by the proved limits only where they bracket those. On the fast branch far from the cusp the largest deviations without a frame, both at τ/σ = 0.25, 2·v_hs, are 41.4 and 2.0 times the proved rules' largest (the first fit, the second; referee round 101). The slow branch far from the cusp by t(v) (both fits' deviations): t(v) 0.635 at τ/σ = 2, 0.5·v_hs (a frame): +6.07e-05, -1.39e-04; t(v) 0.500 at τ/σ = 2, Fock (a frame): +1.42e-05, -6.15e-05; t(v) 0.500 at τ/σ = 2, Mercer-zero (a frame): +1.42e-05, -6.15e-05; t(v) 0.478 at τ/σ = 1, 0.5·v_hs (no frame): +7.55e-06, -4.45e-05; t(v) 0.473 at τ/σ = 2, 0.8·v_hs (a frame): +1.10e-05, -5.21e-05; t(v) 0.382 at τ/σ = 0.5, 0.5·v_hs (no frame): +1.59e-06, -1.68e-05; t(v) 0.347 at τ/σ = 0.25, 0.5·v_hs (no frame): +9.29e-07, -1.07e-05; t(v) 0.277 at τ/σ = 1, 0.8·v_hs (a frame): +3.67e-07, -6.37e-06; t(v) 0.268 at τ/σ = 1, Mercer-zero (a frame): +3.21e-07, -5.81e-06; t(v) 0.200 at τ/σ = 1, Fock (a frame): +2.13e-07, -3.55e-06; t(v) 0.166 at τ/σ = 0.5, 0.8·v_hs (no frame): +2.39e-08, -1.04e-06; t(v) 0.126 at τ/σ = 0.25, 0.8·v_hs (no frame): +6.97e-09, -4.56e-07; t(v) 0.101 at τ/σ = 0.5, Mercer-zero (a frame): +5.65e-09, -3.59e-07.

By the argument attaining t(v), away from the cusp (|ln(v/v_hs)| ≥ 0.05): τ²/v (13 rules): the limit over 1/t(v), minus 1, +5.6e-09 … +6.1e-05; the ratio's excess at K = 160 +8.0e-06 … +4.1e-04, its exponent 1.52 … 1.98; with a K⁻¹ term -1.4e-04 … -3.5e-07, the two fits bracketing 1/t(v) at 13 of the 13; vσ*²/τ² (12 rules): the limit over 1/t(v), minus 1, -2.1e-09 … +8.8e-08; the ratio's excess at K = 160 -3.7e-06 … -1.7e-06, its exponent 1.99 … 2.04; with a K⁻¹ term -1.8e-07 … +3.6e-07, the two fits bracketing 1/t(v) at 4 of the 12. Near it: τ/σ = 0.25, moment matching (v/v_hs 1.030776): +4.6e-08, τ/σ = 0.25, Fock (v/v_hs 0.999541): +1.1e-04, τ/σ = 0.25, Mercer-zero (v/v_hs 0.971825): +2.7e-06, τ/σ = 0.25, Mercer frame (v/v_hs 1.000433): +1.1e-04, τ/σ = 0.5, Fock (v/v_hs 0.993808): +3.2e-05, τ/σ = 0.5, Mercer frame (v/v_hs 1.005090): +2.1e-05, τ/σ = 1, Mercer frame (v/v_hs 1.035276): +2.5e-06.

At v_hs the sketch's integral spreads over every angle, so e_K ∝ K^{1/4}t^K and the ratio per kernel is (1 − 1/(4K))/t + O(K⁻²): the ratio over (1 − 1/(4K))/t, minus 1, at K = 20, 40, 60, 80, 120, 160 (its successive ratios; K⁻²'s 0.250, 0.444, 0.562, 0.444, 0.562): τ/σ = 0.25: -8.0e-05, -2.0e-05, -8.8e-06, -4.9e-06, -2.2e-06, -1.2e-06 (0.247, 0.443, 0.561, 0.444, 0.562); τ/σ = 0.5: -8.3e-05, -2.1e-05, -9.1e-06, -5.1e-06, -2.3e-06, -1.3e-06 (0.248, 0.443, 0.562, 0.444, 0.562); τ/σ = 1: -1.1e-04, -2.7e-05, -1.2e-05, -6.9e-06, -3.1e-06, -1.7e-06 (0.253, 0.446, 0.564, 0.445, 0.563); τ/σ = 2: -2.1e-04, -5.7e-05, -2.6e-05, -1.5e-05, -6.8e-06, -3.9e-06 (0.275, 0.461, 0.574, 0.454, 0.569).

A Mehler frame (|ρ| < 1, ρ = 1 + 2σ²/τ² − 2σ²/v: v ∈ (τ²σ²/σ*², τ²); D‴.61) exists for 22 of the 36 rules (21 of 35 distinct; 4 at τ/σ = 0.25, 4 at τ/σ = 0.5, 6 at τ/σ = 1, 8 at τ/σ = 2); ρ = 1 at 4 (moment matching), |ρ| > 1 at 10.

The Fock rule (ρ = 0) and the Mercer-zero rule have a frame at every τ/σ, so with D‴.61's part VA (its (iv) at every frame, r* = 1; by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's Appendix A.0 assumed correct, now independent cross-checks; referee rounds 101, 107) their rates — λ/(2 − λ) and q — are proved at every τ/σ, no threshold left. The earlier tiers (thresholds rounded inward: referee round 96): with part V (r* = 0.9900) while below that r*, σ*/σ < 19.924 (τ/σ < 19.899) and τ/σ < 140.712 — near the frame interval's ends t(v) → λ, above that r* once σ*/σ > 14.107 —; with the branch and bound alone (r* = 0.8052): σ*/σ < 4.187 (τ/σ < 4.066), τ/σ < 6.514, the frame ends above it once σ*/σ > 3.045.

The Mercer frame's rule's exact squared error over its frame series (A5) is one constant per τ/σ: within 2.3e-16 at τ/σ = 0.25 (K = 19, 20, 39, 40), 2.3e-16 at τ/σ = 0.5 (K = 19, 20, 39, 40), 2.3e-16 at τ/σ = 1 (K = 19, 20, 39, 40), 2.3e-16 at τ/σ = 2 (K = 19, 20, 39, 40) — the resolution of the stored doubles, not an agreement measured below it.

A6 beyond K = 160 on the fast branch (A6X, referee rounds 101, 107, 113, 118, 123, 129): the two fits' deviations from 1/t(v) through K = 80, 120, 160 → through 200, 240, 320, and the limits from the 26-digit errors (t(v) exact) — τ/σ = 0.25, 2·v_hs (no frame; 1225 digits): +8.71e-08 → +1.55e-12 and +3.60e-07 → -3.34e-09, the excess times K² at K = 160, 200, 240, 320 -0.0928, -0.0927, -0.0927, -0.0926; by successive orders (four … eight terms) -6.5e-13, +9.3e-11, -2.9e-09, +3.7e-09, +1.2e-08; with the other branch's transient Bq^K, q = (h_other/t(v))² = 0.850630 per kernel, through K = 120 … 320: +6.5e-16 (B 0.2383; B_K = r_K/q^K at K = 60, 80, outside the fit, within 1.2e-04 of B — at the fit's own K it is B by construction —: no power of K beside q^K to that resolution), its next order (with K⁻⁵, from K = 80) +3.1e-15, the residual against the four-term fit falling by 0.85219, 0.85069, 0.85062, 0.85026 per kernel over K = 20 … 120; its carried orders (four … eight terms, Bq^K in place of the top power, on each order's K) +1.5e-12, +6.5e-16, +3.1e-15, -6.0e-13, -7.8e-12; τ/σ = 0.5, moment matching (no frame; 863 digits): +1.93e-10 → +9.03e-12 and +1.91e-08 → +1.81e-09, the excess times K² at K = 160, 200, 240, 320 -0.0518, -0.0519, -0.0519, -0.0519; by successive orders (four … eight terms) -3.9e-14, -2.1e-16, +1.9e-17, -4.8e-19, -7.4e-20; τ/σ = 1, moment matching (no frame; 550 digits): -7.52e-11 → -3.53e-12 and -5.89e-08 → -5.85e-09, the excess times K² at K = 160, 200, 240, 320 -0.0865, -0.0864, -0.0863, -0.0862; by successive orders (four … eight terms) +8.6e-15, -3.3e-16, +2.1e-18, +7.7e-19, -1.2e-19; τ/σ = 2, moment matching (no frame; 328 digits): -1.68e-11 → -9.00e-13 and -3.52e-08 → -3.51e-09, the excess times K² at K = 160, 200, 240, 320 -0.0471, -0.0470, -0.0470, -0.0469; by successive orders (four … eight terms) -3.6e-14, +5.2e-16, -3.6e-18, -2.9e-20, -1.0e-20; τ/σ = 2, 2·v_hs (a frame; 328 digits): -6.94e-11 → -3.35e-12 and -4.49e-08 → -4.46e-09, the excess times K² at K = 160, 200, 240, 320 -0.0517, -0.0516, -0.0515, -0.0515; by successive orders (four … eight terms) -2.0e-14, +4.9e-16, -7.4e-18, +1.9e-19, -3.2e-20; τ/σ = 1, Mercer frame (a frame; 550 digits): +2.55e-08 → +9.76e-10 and +7.55e-07 → +6.20e-08, the excess times K² at K = 160, 200, 240, 320 +0.0631, +0.0619, +0.0611, +0.0602; by successive orders (four … eight terms) -7.0e-11, +8.5e-12, -2.2e-12, +1.2e-12, -1.7e-12; the other branch's transient Bq^K, q = (h_other/t(v))² = 0.674644 per kernel, not carried in the six, seven and eight-term orders, where q^K > 10⁻²⁰ at their smallest K (80, 60, 40): q^K 2.1e-14, 5.6e-11, 1.5e-07 there is below the order's own remainder, its next change (the power the order would drop for Bq^K, at that K, by the next order with Bq^K) 1.3e-08, 3.2e-08, 1.9e-06; on K ≥ 160, where the orders converge (A6C), q^160 = 4.5e-28; τ/σ = 2, Mercer frame (a frame; 328 digits): -2.10e-09 → -9.69e-11 and -1.73e-07 → -1.62e-08, the excess times K² at K = 160, 200, 240, 320 -0.0455, -0.0452, -0.0450, -0.0447; by successive orders (four … eight terms) +1.2e-12, +2.5e-13, -1.3e-13, +7.9e-14, -1.4e-13; the other branch's transient Bq^K, q = (h_other/t(v))² = 0.604938 per kernel, not carried in the six, seven and eight-term orders, where q^K > 10⁻²⁰ at their smallest K (80, 60, 40): q^K 3.4e-18, 8.0e-14, 1.9e-09 there is below the order's own remainder, its next change (the power the order would drop for Bq^K, at that K, by the next order with Bq^K) 4.4e-10, 1.9e-09, 5.4e-08; on K ≥ 160, where the orders converge (A6C), q^160 = 1.2e-35 — at the 4 rules without a frame the deviations shrink by 18.6–5.63e+04 and 10–108 times, in the later window (160, 200, 240 → 200, 240, 320) by 2.47–3.18 and 1.99–2.03; at the 3 proved rules the deviations shrink by 20.7–26.2 and 10–12.2 times, in the later window (160, 200, 240 → 200, 240, 320) by 2.52–2.65 and 2–2.09; a pure K⁻⁴ and K⁻³ remainder gives 21.11 and 10, 2.533 and 2 (window geometry, any rule); every deviation shrinks as the window moves out (a necessary condition); by the orders (four … eight terms, with Bq^K wherever it is carried; an order below its propagated input bound — the 26-digit errors' rounding through the fit, here at most 6.8e-24 — counting as converged), rule by rule at the 6 rules without a transient resolvable at K = 160, (τ/σ, rule) = (0.5, moment matching), (1, moment matching), (2, moment matching), (2, 2·v_hs): decreasing in magnitude through eight terms, the last within ten times its last change (that, the resolution: 4.1e-18, 9.0e-18, 1.9e-19, 2.3e-18) — 1/t(v) to that resolution; (1, Mercer frame), (2, Mercer frame): NOT decreasing in magnitude through eight terms — inconclusive: the limit, the premise, or the window (the expansion not yet asymptotic at the order's smallest K); at (1, Mercer frame), (2, Mercer frame) the window: on K ≥ 160 the orders decrease (A6C), to -1.1e-21, -7.3e-24 at sixteen terms, the last within ten times its last change, the excess's coefficients growing 12.4–17.5, 22.7–24.2 times per order (the last four ratios c_k/c_{k−1}; at most 6.3 at A6C's rules where the check passes), so that at the eight-term order's smallest K, 40, the terms c_k/K^k shrink by only 0.31–0.44, 0.56–0.61 per order: the expansion is not yet asymptotic there; the known cases, the 3 rules proved: the check passes at 1 of 3 — a failure is inconclusive, it does not show the limit off; at (1, Mercer frame), (2, Mercer frame) the other branch's transient, resolvable in the 26-digit errors at some orders' smallest K, lies below the order's own remainder at each: not carried, it does not account for the stall, and where the orders converge, on K ≥ 160 (A6C), q^160 is at most 4.5e-28 (a negative result); the next term (referee round 129, derived from D‴.61 (ii)'s moments by a two-dimensional Laplace step at the mass point; formal at next order — D‴.61 proves only the leading limit, the next-order uniformity and the Laplace remainder are not proved): e_{K−1}/e_K = (1/t(v))(1 + c₂/K² + O(K⁻³)), c₂ = s₁/2, s₁ = 1/16 − ϱ/(4(1 + ϱ)²) − 2ϱ(1 − ϱ)/(1 + ϱ)³ + √q(1 − ϱ)²/(8(1 − √q)(1 + ϱ)²), ϱ = t(v) on the fast branch (−t(v) on the slow), √q = h_other/t(v) signed; per rule c₂ and a fitted one (the limit exact; A6C's where computed, else from the stored errors, the eight-term order's seven powers through the last seven K), relative: (0.25, 2·v_hs) -0.0924369845420337681 and -9.2436984542038e-02 (A6C), 4.4e-14; (0.5, moment matching) -0.052 and -5.1999999999997e-02 (the stored errors), 6.1e-14; (1, moment matching) -0.0859375 and -8.5937499999995e-02 (the stored errors), 5.9e-14; (2, moment matching) -0.04675 and -4.6750000000000e-02 (A6C), 1.2e-41; (2, 2·v_hs) -0.0512432244843701805 and -5.1243224484370e-02 (A6C), 1.6e-42; (1, Mercer frame) 0.0575470734670659269 and +5.7547073467066e-02 (A6C), 1.2e-18; (2, Mercer frame) -0.044 and -4.4000000000000e-02 (A6C), 9.4e-22; with c₂ fixed (R_K − c₂/(t(v)K²) = L + Σ_{k=3..n+1} c_k/K^k through the last n K: the same top powers, one K fewer) the orders (0.5, moment matching) +1.47e-14, +4.56e-17, -2.53e-18, +4.98e-20, +3.84e-21, decreasing, the last within ten times its last change (4.6e-19); (1, moment matching) -3.16e-15, +7.56e-17, -3.58e-19, -6.62e-20, +6.21e-21, decreasing, the last within ten times its last change (7.3e-19); (2, moment matching) +1.33e-14, -1.20e-16, +4.97e-19, +3.34e-21, +5.44e-22, decreasing, the last within ten times its last change (2.8e-20); (2, 2·v_hs) +7.36e-15, -1.14e-16, +9.89e-19, -1.60e-20, +1.64e-21, decreasing, the last within ten times its last change (1.8e-19); (1, Mercer frame) +2.41e-11, -1.71e-12, +2.22e-13, -6.47e-14, +4.49e-14, decreasing, the last within ten times its last change (1.1e-12); (2, Mercer frame) -5.14e-13, -4.39e-14, +1.20e-14, -4.41e-15, +2.75e-15, decreasing, the last within ten times its last change (7.2e-14) — the known cases, the 3 rules proved: with c₂ fixed the check passes at 3 of 3; at the transient's 1 rule ((0.25, 2·v_hs): resolvable at every order's smallest K) the power fits worsen with more terms, B > 0, B_K flat (within 0.001 of B outside the fit), and the transient fit is within 6.6e-16 (rounded up) of 1/t(v), inside its next order's change (or its input bound, if larger) 2.5e-15 (that, the resolution; the transient's own q^{2K} term unmodelled) — the transient's own criterion, not the orders' decrease: its carried orders do NOT decrease in magnitude — numerical, consistent with an expansion in integer powers of 1/K from K⁻² (no K⁻¹ term) and the transient (a non-integer power can pass this check), not a proof; a proved rule's bias calibrates nothing here, since the proved rules' biases differ (by 29.0 times at τ/σ = 2: -3.35e-12 and -9.69e-11, the first fit through K = 200, 240, 320) and at τ/σ = 0.25 no far fast-branch rule has a frame (A6X extends 4 of the 9 fast-branch rules without a frame).

A6X's rules at consecutive K (A6C, referee round 129): the exact errors e_K at every K = 159 … 400 in Arb (python-flint 0.9.0: the Gauss–Hermite nodes from LAPACK, refined by Newton in Arb; e_K² from the Gaussian inner products, independent of A6's decimal arithmetic; ⌈2K·log₂(1/t(v))⌉ + 450 bits, 637–2377; stored to 70 digits, the balls' relative radii at most 3.5e-86), within 2.4e-26 (relative) of A6's and A6X's 26-digit errors at K − 1 and K for K = 160, 200, 240, 320 (their rounding at most 5.0e-26); from the ratios in exact rationals, t(v) exact, the orders by windows on [160, 400] (n + 1 K spread uniformly in 1/K, the limit free; four … sixteen terms), the residual floor (N powers K⁻² … K^{−N−1} of the excess t(v)R_K − 1, the limit exact, through N such K; the largest residual over all 241 K, N = 6, 10, 14, 18, 22, 26), the K⁻¹ coefficient (times t(v); beside 15 powers) and the excess's coefficients c_k (the limit exact; the median of 6 fits, 20, 24, 28 powers on [160, 400] and [200, 400], printed while they agree within 1%) — (2, Mercer frame; a frame; q = 0.604938, q^160 = 1.2e-35): the orders +7.9e-13, +6.4e-14, -8.0e-15, +7.6e-16, -7.5e-17, +7.6e-18, -8.4e-19, +1.0e-19, -1.3e-20, +1.8e-21, -2.7e-22, +4.3e-23, -7.3e-24 — decreasing in magnitude, the last within ten times its last change (5.1e-22; the input bound at most 1.4e-60); the residual floor 6.9e-18, 9.8e-24, 5.3e-29, 6.8e-34, 9.5e-38, 9.4e-41; the K⁻¹ coefficient +3.1e-18; c_2 … c_11 -0.044, -0.235536, -1.09167, -3.56048, +28.2686, +1174.07, +28436.5, +644380, +1.48467e+07, +3.56822e+08, the ratios c_k/c_{k−1} +5.4, +4.6, +3.3, -7.9, +41.5, +24.2, +22.7, +23.0, +24.0; (1, Mercer frame; a frame; q = 0.674644, q^160 = 4.5e-28): the orders -4.4e-11, +2.3e-12, -1.5e-13, +1.1e-14, -1.0e-15, +1.1e-16, -1.4e-17, +2.1e-18, -3.6e-19, +7.0e-20, -1.5e-20, +3.9e-21, -1.1e-21 — decreasing in magnitude, the last within ten times its last change (5.0e-20; the input bound at most 1.2e-60); the residual floor 1.0e-16, 2.1e-22, 5.0e-27, 9.4e-31, 1.9e-34, 2.1e-36; the K⁻¹ coefficient +3.1e-16; c_2 … c_8 +0.0575471, +0.818862, +9.72693, +120.662, +1658.57, +25668.8, +449187, the ratios c_k/c_{k−1} +14.2, +11.9, +12.4, +13.7, +15.5, +17.5; (2, 2·v_hs; a frame; q = 0.0077138, q^160 = 9.2e-339): the orders -1.3e-14, +1.5e-16, -6.7e-19, +3.4e-21, -8.2e-23, +1.7e-24, -2.4e-26, +3.0e-28, -5.6e-30, +1.4e-31, -3.2e-33, +7.1e-35, -1.7e-36 — decreasing in magnitude, the last within ten times its last change (7.3e-34; the input bound at most 1.1e-60); the residual floor 3.0e-23, 2.8e-32, 8.1e-41, 5.6e-49, 1.2e-56, 4.7e-64; the K⁻¹ coefficient +4.5e-30; c_2 … c_15 -0.0512432, -0.0679679, -0.039581, +0.0421205, +0.146609, +0.186194, +0.219848, +1.39776, +8.05706, +29.9115, +92.4443, +387.04, +2376.16, +14282.7, the ratios c_k/c_{k−1} +1.3, +0.6, -1.1, +3.5, +1.3, +1.2, +6.4, +5.8, +3.7, +3.1, +4.2, +6.1, +6.0; (2, moment matching; no frame; q = 0, q^160 = 0): the orders -2.3e-14, +1.6e-16, -3.4e-19, -8.2e-22, -2.7e-23, +6.4e-25, -3.3e-27, -4.5e-29, +3.4e-31, +1.5e-32, -2.6e-34, -2.2e-36, +8.6e-38 — decreasing in magnitude, the last within ten times its last change (2.4e-35; the input bound at most 1.6e-60); the residual floor 7.0e-24, 4.2e-33, 2.6e-42, 1.4e-51, 9.0e-59, 1.2e-65; the K⁻¹ coefficient -1.4e-31; c_2 … c_15 -0.04675, -0.053814, -0.0114805, +0.0777506, +0.156395, +0.098873, -0.0795635, +0.435776, +3.3366, +5.42581, -15.0794, -50.497, +281.765, +1776.77, the ratios c_k/c_{k−1} +1.2, +0.2, -6.8, +2.0, +0.6, -0.8, -5.5, +7.7, +1.6, -2.8, +3.3, -5.6, +6.3; (0.25, 2·v_hs; no frame; q = 0.85063, q^160 = 5.7e-12): on [200, 400] the power orders -1.1e-15, +3.4e-15, -9.2e-15, +2.3e-14, -5.1e-14, +9.7e-14, -1.6e-13, +2.1e-13, -2.2e-13, +1.3e-13, +5.0e-15, -1.1e-13, +8.0e-14; with q^K(B₀ + B₁/K) in the basis +9.6e-17, -1.4e-17, -4.0e-20, +2.3e-22, +1.5e-22, -5.2e-22, +1.6e-21, -4.1e-21, +9.7e-21, -2.0e-20, +2.6e-20, -9.4e-21, -1.1e-20 (B₀ -1.012811, 0.204939, 0.238343, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317, 0.238317; B₁/B₀ -0.0063 at sixteen terms); with B₂/K² too +7.8e-17, -1.1e-17, -3.1e-20, +2.1e-22, +2.4e-24, +1.4e-26, -4.8e-26, +1.6e-25, -4.7e-25, +1.3e-24, -3.8e-24, -9.5e-23, +1.6e-24 (B₀ -59.245390, -0.853179, 0.238960, 0.238317, 0.238316, 0.238316, 0.238316, 0.238316, 0.238316, 0.238316, 0.238316, 0.238316, 0.238316; B₁/B₀ -0.0050: a small 1/K term in B) — the transient's own next term, q^{2K} = 7.9e-29 at K = 200, unmodelled — numerical, not a proof.

The step of the rate's Laplace sketch (A6E): the rule's defect generating function ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] — Σ_n D_n(v′)sⁿ/n! = e^{(v′−1)s²/2}ε̃_K(√v′s), so that the frame series is its Fock norm, (πρv′)⁻¹∫|ε̃_K(ξ)|²e^{−Ax²−By²}dA, A = (1 − ρ(v′ − 1))/(ρv′), B = (1 + ρ(v′ − 1))/(ρv′) — over (K!/(2K)!)ξ^{2K}e^{−ξ²/4}, on ξ = s√K and ξ = is√K for s = 0.1 … 2 (`gen_ratio`, 3K + 60 digits): K = 20, the real axis: 0.804–1.083 (at s = 0.5, 1, 1.5, 2: 0.973, 0.894, 0.807, 1.083); K = 20, the imaginary axis: 1.001–1.408 (at s = 0.5, 1, 1.5, 2: 1.028, 1.110, 1.241, 1.408); K = 40, the real axis: 0.741–1.052 (at s = 0.5, 1, 1.5, 2: 0.971, 0.882, 0.758, 1.052); K = 40, the imaginary axis: 1.001–1.411 (at s = 0.5, 1, 1.5, 2: 1.029, 1.114, 1.246, 1.411); K = 80, the real axis: 0.677–1.021 (at s = 0.5, 1, 1.5, 2: 0.969, 0.874, 0.721, 1.021); K = 80, the imaginary axis: 1.001–1.413 (at s = 0.5, 1, 1.5, 2: 1.030, 1.116, 1.248, 1.413). Bounded above at these K; on the real axis its minimum falls with K toward D‴.61's limit (1 − s²/4)^{1/2}, which vanishes at s = 2 (referee round 88), so its K-th root tends to 1 where s < 2; Laplace's method on the Fock integral — |ξ|^{4K}e^{−(A + ½)x² − (B − ½)y²}, its mass at |ξ| = (2K/min(A + ½, B − ½))^{1/2} — then gives the rate 1/(2 min(A + ½, B − ½)), the formula's max{h(τ²/v), h(vσ*²/(σ²τ²))} (a sketch here; D‴.61 (iii)–(iv) prove it for every rule with a Mehler frame, by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's Appendix A.0 assumed correct, now independent cross-checks, and D‴.62 (b) for every rule without one, by hand).

The tail's share Σ_{n>2K}/(first term) against ₂F₁(−½, −½; 1; t²) − 1 (the share's limit if the ratios' leading coefficients govern the tail):

| τ/σ | t | K = 5 | 10 | 20 | 30 | ₂F₁(−½, −½; 1; t²) − 1 |
|---|---|---|---|---|---|---|
| 0.25 | 0.01515 | 0.00003 | 0.00004 | 0.00005 | 0.00005 | 0.00006 |
| 0.5 | 0.05573 | 0.00038 | 0.00054 | 0.00065 | 0.00069 | 0.00078 |
| 1 | 0.17157 | 0.00367 | 0.00518 | 0.00617 | 0.00655 | 0.00737 |
| 2 | 0.38197 | 0.02064 | 0.02642 | 0.03074 | 0.03251 | 0.03682 |

## B. The transit

The sliding path (polished; a kernel leaves only if its correlation stays inward with its own contribution removed) on f* = 0.8φ_σ*(· − 0.03) sampled at 321 points of [-1.6, 1.6], 21 kernels started at the grid −1, −0.9, …, 1, to 13 entries, the fast flow's chunks up to 10⁵·2^24: after each entry the relaxed block's relative rms error against the optimal fit with as many kernels on the same samples (the better of variable projection from the block and from the half-squeeze nodes: the best local optimum found); a block is optimal within 0.001. The numerical variants: base (polish step 1σ, grid shifted by 0, the polish's budget 5000 evaluations, its gradient test 1e-15), polish 2σ (polish step 2σ, grid shifted by 0, the polish's budget 5000 evaluations, its gradient test 1e-15), grid +0.05 (polish step 1σ, grid shifted by 0.05, the polish's budget 5000 evaluations, its gradient test 1e-15), polish 2σ, 5·10⁴ evaluations (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, its gradient test 1e-15), polish 2σ, 5·10⁴, repeated (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, repeated up to 3 times, its gradient test 1e-15), polish 2σ, 5·10⁴, repeated, no gradient test (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, repeated up to 3 times, no gradient test). 'Relaxed': the block passed the stationarity test. The first non-optimal block's polishes as status/evaluations/accepted (status 0: the budget was spent first; 1: the gradient test), a repeated polish's end over the optimum in brackets. 'none': every block optimal, to the path's cap of 12 kernels (the range censored).

| τ/σ | variant | optimal blocks (kernels) | the first non-optimal: kernels, its optimum's error, block / optimum, relaxed | its polishes | relaxed blocks | chunks redone (stall, failure), failed | polish rejected | exits |
|---|---|---|---|---|---|---|---|---|
| 0.25 | base | 1, 2, 3, 4 | 5, 1.61e-09, 10.85, True | 0/5000/a | 12/12 | 0, 0, 0 | 13 | 1 |
| 0.25 | polish 2σ | 1, 2, 3, 4 | 5, 1.61e-09, 10.95, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 2 |
| 0.25 | grid +0.05 | 1, 2, 3, 4 | 5, 1.61e-09, 10.96, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 2 |
| 0.25 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4 | 5, 1.61e-09, 1.012, True | 1/30240/a | 12/12 | 0, 0, 0 | 2 | 1 |
| 0.25 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4 | 5, 1.61e-09, 1.012, True | 1/30240/a 1/1/r (1.012) | 12/12 | 0, 0, 0 | 2 | 0 |
| 0.25 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5 | 6, 2.55e-11, 4.951, True | 0/50000/a 0/50000/a (8.7) 0/50000/a (6.609) 0/50000/a (4.951) | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.5 | base | 1, 2, 3, 4, 5 | 6, 6.30e-08, 14.34, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 3, 5, 0 | 182 | 0 |
| 0.5 | polish 2σ | 1, 2, 3, 4, 5 | 6, 6.30e-08, 1.46, True | 0/5000/a | 12/12 | 0, 0, 0 | 3 | 1 |
| 0.5 | grid +0.05 | 1, 2, 3, 4, 5 | 6, 6.30e-08, 14.19, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 2, 5, 0 | 143 | 0 |
| 0.5 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 4.605, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.5 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 2.677, True | 0/50000/a 0/50000/a (2.677) 1/1/r (2.677) | 12/12 | 0, 0, 0 | 1 | 0 |
| 0.5 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 1.233, True | 0/50000/a 0/50000/a (2.596) 0/50000/a (1.708) 0/50000/a (1.233) | 12/12 | 1, 0, 0 | 0 | 0 |
| 0.75 | base | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 7.53, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 1, 7, 0 | 138 | 0 |
| 0.75 | polish 2σ | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 1.642, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 8.656, True | 0/5000/r | 12/12 | 1, 3, 0 | 88 | 0 |
| 0.75 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 6.81e-10, 2.305, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 6.81e-10, 1.121, True | 0/50000/a 0/50000/a (1.121) 1/1/r (1.121) | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 | 11, 7.75e-11, 2.685, True | 0/50000/a 0/50000/a (3.932) 0/50000/a (3.203) 0/50000/a (2.685) | 12/12 | 1, 0, 0 | 0 | 0 |
| 1 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 5.22e-08, 1.771, True | 0/5000/a | 12/12 | 0, 0, 0 | 12 | 0 |
| 1 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 5.22e-08, 1.771, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 10 | 9, 2.96e-07, 4.241, True | 1/3635/r 1/3620/r 1/3720/r 1/3725/r 1/3618/r 1/3600/r 1/3648/r 1/3682/r 1/3722/r 1/3604/r 1/3685/r 1/3689/r 1/3684/r 1/3608/r 1/3624/r 1/3594/r 1/3881/r 1/3658/r 1/4129/r 1/3699/r 1/3648/r 1/3718/r 1/3616/r 1/3553/r | 12/12 | 2, 3, 0 | 110 | 0 |
| 1 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 | 12, 1.61e-09, 1.465, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 1, 0 | 71 | 0 |
| 1.5 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 13 | 0 |
| 2 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |

Over every width and variant: the optimal blocks' optimum errors down to 6.81e-10; the non-optimal blocks' from 2.96e-07 down (0 of 116 failed the stationarity test; 53 of those that passed it pass it at any state by Cauchy–Schwarz, the test being absolute — their errors up to 4.0e-10; for the other 63 its threshold over the Cauchy–Schwarz bound — the residual's cosine with the block's directions it allows — is 3.5e-04–0.97 (median 7.7e-03; above 0.1 in 24): a weak test); non-optimal blocks with an optimum above 10⁻⁷: 1; the non-optimal blocks' errors 1.01e-12 … 1.26e-06.

Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in 33 of the 38 blocks whose optimum is at least 1e-10; 1 ends within 1% of it (1.001 times: in its basin, unconverged); 4 end 16.36–24.74 times above it — the blocks' largest node distances from the optimum's 2.14–3.99σ — where variable projection stalls: its two closest nodes 0.000–0.021σ apart, its end not stationary (the largest cosine of its residual with a Jacobian column 2.6e-04–4.7e-03), its stop 3/2554, 0/20000, 3/2166, 3/18995 (status/evaluations: a miss says nothing of the block's basin). The 78 blocks with an optimum below 1e-10, where the solve stops short of the optimum (A), are not counted (0 reached). From the first non-optimal block of each path: 20 of 20; on 2 paths the first non-optimal block's optimum is below 1e-10. There the Jacobian's condition number is 1.3e+06–5e+08 for the joint fit in (w, c) with the centres in x, the polish's parameters (4.5e+06–2.3e+09 with the centres in units of σ) and 5.4e+03–3.8e+05 for variable projection's in c alone (analytic; the same in either unit of the centres); the ratio of their smallest singular values, variable projection's over the joint fit's, is at least 1 in any common scaling of the centres for Kaufman's form of the former, the centre columns projected off the amplitudes' (referee round 77: the joint fit's smallest singular value is at most theirs — a Schur complement), from which the analytic Golub–Pereyra form used here differs by a term proportional to the residual; its size depends on the units: 1–31.7 (median 1.01) with the centres in units of σ, 1–31.7 (median 1.54) in x, 49.2–2.79e+03 (median 570) with every column scaled to unit norm, which scales the two problems' centre columns differently; TRF itself works in Coleman–Li's variables, the bounded amplitudes scaled by v^{1/2}.


### B3. The range ends: the polish restarted from its recorded start, one setting at a time (referee rounds 65, 70, 77)

Each path's first non-optimal block (its error at least 1.001 times the optimum's): the path rerun with every polish's start recorded, and from the block's first polish's start — its error over the optimum's — the polish's own solve (`polish_block`, TRF) with the path's budget and gradient test, with 500000 evaluations and the test, with the path's budget without the test, and with both changed: each its stop (status 0: the budget; 1: the gradient test; 2: the decrease test; 3: the step test)/evaluations, its error over the optimum's and the path's seed rule on its end (the largest move of a kernel near its seed, σ; a: accepted within the path's polish step, r: rejected); the block the path ended with; † the optimum below 1e-10:

| τ/σ | variant | K | optimum | start | the path's block | path | budget | no test | budget, no test | bounded by |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | base | 5 | 1.61e-09 | 62.98 | 10.85 | 0/5000, 10.84750, 0.44σ a | 1/31071, 1.00169, 0.13σ a | 0/5000, 10.84750, 0.44σ a | 3/31740, 1.00000, 0.13σ a | the budget and the gradient test |
| 0.5 | base | 6 | 6.30e-08 | 17.05 | 14.34 | 0/5000, 1.32109, 1.47σ r | 1/6444, 1.00000, 1.51σ r | 0/5000, 1.32109, 1.47σ r | 3/6463, 1.00000, 1.51σ r | the budget (and the radius) |
| 0.75 | base | 8 | 5.22e-08 | 8.567 | 7.53 | 0/5000, 1.59270, 1.50σ r | 1/11030, 1.00000, 1.64σ r | 0/5000, 1.59270, 1.50σ r | 3/11046, 1.00000, 1.64σ r | the budget (and the radius) |
| 1 | base | 10 | 5.22e-08 | 5.767 | 1.771 | 0/5000, 1.77117, 0.90σ a | 1/11923, 1.00000, 1.17σ r | 0/5000, 1.77117, 0.90σ a | 3/11938, 1.00000, 1.17σ r | the budget (and the radius) |
| 0.25 | polish 2σ | 5 | 1.61e-09 | 62.98 | 10.95 | 0/5000, 10.94555, 0.45σ a | 1/30240, 1.01196, 0.13σ a | 0/5000, 10.94555, 0.45σ a | 3/31769, 1.00000, 0.13σ a | the budget and the gradient test |
| 0.5 | polish 2σ | 6 | 6.30e-08 | 17.05 | 1.46 | 0/5000, 1.46063, 1.45σ a | 1/6635, 1.00000, 1.51σ a | 0/5000, 1.46063, 1.45σ a | 3/6654, 1.00000, 1.51σ a | the budget |
| 0.75 | polish 2σ | 8 | 5.22e-08 | 8.567 | 1.642 | 0/5000, 1.64194, 1.49σ a | 1/11196, 1.00000, 1.64σ a | 0/5000, 1.64194, 1.49σ a | 3/11213, 1.00000, 1.64σ a | the budget |
| 1 | polish 2σ | 10 | 5.22e-08 | 5.767 | 1.771 | 0/5000, 1.77117, 0.90σ a | 1/11923, 1.00000, 1.17σ a | 0/5000, 1.77117, 0.90σ a | 3/11938, 1.00000, 1.17σ a | the budget |
| 0.25 | grid +0.05 | 5 | 1.61e-09 | 62.46 | 10.96 | 0/5000, 10.95913, 0.44σ a | 1/30286, 1.01251, 0.12σ a | 0/5000, 10.95913, 0.44σ a | 3/31833, 1.00000, 0.12σ a | the budget and the gradient test |
| 0.5 | grid +0.05 | 6 | 6.30e-08 | 17.14 | 14.19 | 0/5000, 1.41071, 1.04σ r | 1/6567, 1.00000, 1.10σ r | 0/5000, 1.41071, 1.04σ r | 3/6585, 1.00000, 1.10σ r | the budget (and the radius) |
| 0.75 | grid +0.05 | 8 | 5.22e-08 | 8.656 | 8.656 | 0/5000, 1.70316, 1.07σ r | 3/11386, 1.00000, 1.22σ r | 0/5000, 1.70316, 1.07σ r | 3/11386, 1.00000, 1.22σ r | the budget (and the radius) |
| 1 | grid +0.05 | 9 | 2.96e-07 | 5.592 | 4.241 | 1/3635, 1.00000, 1.88σ r | 1/3635, 1.00000, 1.88σ r | 3/3643, 1.00000, 1.88σ r | 3/3643, 1.00000, 1.88σ r | the radius |
| 0.25 | polish 2σ, 5·10⁴ evaluations | 5 | 1.61e-09 | 62.98 | 1.012 | 1/30240, 1.01196, 0.13σ a | 1/30240, 1.01196, 0.13σ a | 3/31769, 1.00000, 0.13σ a | 3/31769, 1.00000, 0.13σ a | the gradient test |
| 0.5 | polish 2σ, 5·10⁴ evaluations | 8 | 2.10e-10 | 445 | 4.605 | 0/50000, 4.60526, 1.20σ a | 1/252146, 1.01569, 0.83σ a | 0/50000, 4.60526, 1.20σ a | 3/275355, 1.00022, 0.82σ a | the budget and the gradient test |
| 0.75 | polish 2σ, 5·10⁴ evaluations | 10 | 6.81e-10 | 150.1 | 2.305 | 0/50000, 2.30523, 1.33σ a | 1/143369, 1.00003, 1.08σ a | 0/50000, 2.30523, 1.33σ a | 3/144109, 1.00000, 1.08σ a | the budget |
| 1 | polish 2σ, 5·10⁴ evaluations | 12 | 1.61e-09 | 57.9 | 1.465 | 0/50000, 1.46508, 1.30σ a | 1/102753, 1.00000, 1.46σ a | 0/50000, 1.46508, 1.30σ a | 3/102930, 1.00000, 1.46σ a | the budget |
| 0.25 | polish 2σ, 5·10⁴, repeated | 5 | 1.61e-09 | 62.98 | 1.012 | 1/30240, 1.01196, 0.13σ a | 1/30240, 1.01196, 0.13σ a | 3/31769, 1.00000, 0.13σ a | 3/31769, 1.00000, 0.13σ a | the gradient test |
| 0.5 | polish 2σ, 5·10⁴, repeated | 8 | 2.10e-10 | 445 | 2.677 | 0/50000, 4.60526, 1.20σ a | 1/252146, 1.01569, 0.83σ a | 0/50000, 4.60526, 1.20σ a | 3/275355, 1.00022, 0.82σ a | the budget and the gradient test |
| 0.75 | polish 2σ, 5·10⁴, repeated | 10 | 6.81e-10 | 150.1 | 1.121 | 0/50000, 2.30517, 1.33σ a | 1/143369, 1.00000, 1.08σ a | 0/50000, 2.30517, 1.33σ a | 3/144109, 0.99998, 1.08σ a | the budget |
| 0.25 | polish 2σ, 5·10⁴, repeated, no gradient test | 6 | 2.55e-11† | 2420 | 4.951 | 0/50000, 13.52667, 0.89σ a | 3/450230, 1.02346, 0.45σ a | 0/50000, 13.52667, 0.89σ a | 3/450230, 1.02346, 0.45σ a | neither (the solve's stop, status 3) |
| 0.5 | polish 2σ, 5·10⁴, repeated, no gradient test | 8 | 2.10e-10 | 445 | 1.233 | 0/50000, 4.57583, 1.19σ a | 3/278779, 1.00016, 0.82σ a | 0/50000, 4.57583, 1.19σ a | 3/278779, 1.00016, 0.82σ a | the budget |
| 0.75 | polish 2σ, 5·10⁴, repeated, no gradient test | 11 | 7.75e-11† | 1351 | 2.685 | 0/50000, 5.11693, 0.86σ a | 0/500000, 1.04673, 1.72σ a | 0/50000, 5.11693, 0.86σ a | 0/500000, 1.04673, 1.72σ a | neither (the budget spent) |

The 'path' restart reproduces the first polish B stored — its stop and evaluations — at 22 of 22 blocks.

Above 1e-10 (20 blocks: τ/σ (variant)) the end is bounded by

- the budget: 7 — 0.5 (polish 2σ); 0.75 (polish 2σ); 1 (polish 2σ); 0.75 (polish 2σ, 5·10⁴ evaluations); 1 (polish 2σ, 5·10⁴ evaluations); 0.75 (polish 2σ, 5·10⁴, repeated); 0.5 (polish 2σ, 5·10⁴, repeated, no gradient test);
- the budget (and the radius): 5 — 0.5 (base); 0.75 (base); 1 (base); 0.5 (grid +0.05); 0.75 (grid +0.05);
- the budget and the gradient test: 5 — 0.25 (base); 0.25 (polish 2σ); 0.25 (grid +0.05); 0.5 (polish 2σ, 5·10⁴ evaluations); 0.5 (polish 2σ, 5·10⁴, repeated);
- the gradient test: 2 — 0.25 (polish 2σ, 5·10⁴ evaluations); 0.25 (polish 2σ, 5·10⁴, repeated);
- the radius: 1 — 1 (grid +0.05).

Below it (2 blocks) the restart with both changed stops at 1.0235 (3/450230: neither (the solve's stop, status 3)), 1.0467 (0/500000: neither (the budget spent)) times the stored optimum, 1.0235, 1.0484 times B4's best end.

The restart with both changed (the polish's TRF): its end over the optimum, minus 1 — below 2e-09 over B4's best end —, by the optimum (one solver from one kind of start; above B4_MAX a negative value is the end below the stored optimum, converged only to that): 2.96e-07: -5.3e-10; 6.3e-08: -9.6e-11 … -5.4e-11; 5.22e-08: -6.3e-09 … -5.3e-10; 1.61e-09: +3.9e-08 … +1.7e-06; 6.81e-10: +2.7e-05; 2.1e-10: +1.7e-04 … +2.3e-04; 7.75e-11: +4.8e-02; 2.55e-11: +2.3e-02.


### B4. The floor: the same solve with the residual in double, in long double and about an anchor, evaluated exactly (referee rounds 70, 77, 81, 84)

B3's blocks whose optimum is below 2e-09, from the path's block and from the stored optimum (its centres with least-squares amplitudes): Levenberg–Marquardt on the joint fit in (w, c) (`lm_fit`: the Jacobian and the step in double; a step accepted when its predicted decrease is positive and the objective, as the solve computes its residual, decreases) with the parameters and the residual f − Φ(c)w by direct subtraction in double or in long double (its ε 1.08e-19), or about an anchor (w₀, c₀) — r = r₀ − Σ_kφ₀_k[(w_k − w₀_k) + w_k·expm1(Δ_k)], Δ_k = (c_k − c₀_k)(2x − c_k − c₀_k)/(2σ²), r₀ by direct subtraction: in double re-anchored at each stall (μ reset) and every 2000 accepted steps, or not at stalls, and in long double from the anchored double's end —, and the polish's own TRF (`polish_block`, 500000 evaluations, no gradient test). Every error is exact (`exact_resid`: the residual in 60-digit decimal arithmetic, every double and long double taken exactly); the start's error over the stored optimum's; the block's best end — the lowest error any solve reached from either start — over it; and each end's excess over that best (formed exactly, `b4_excess`), its iterations/accepted steps and stop (stall: no step the solve computes decreases its objective; flat: the exact rms fell by less than 1e-10 of itself over 2000 accepted steps; cap: 400000 iterations) or TRF's status/evaluations; † the optimum below 1e-10:

| τ/σ | variant | K | optimum | start | its error | the best end | double | long double | double, anchored | double, anchored, no stall re-anchoring | long double, anchored, from the anchored double's end | TRF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | base | 5 | 1.61e-09 | block | 10.8475 | 0.999999986214 | 1.9e-08, 7996/7957, stall | 2.4e-22, 8011/7970, stall | 2.6e-16, 8626/8003, stall | 9.6e-13, 8006/7966, stall | 0.0e+00, 57/1, stall | 1.5e-07, 3/26621 |
| 0.25 | base | 5 | 1.61e-09 | optimum | 1 | 0.999999986214 | 2.4e-08, 28/0, stall | 1.1e-16, 46/11, stall | 9.5e-17, 8011/481, stall | 2.4e-16, 46/11, stall | 4.9e-21, 677/36, stall | 2.4e-08, 3/19 |
| 0.5 | polish 2σ, 5·10⁴ evaluations | 8 | 2.10e-10 | block | 4.60526 | 0.999989242094 | 5.6e-05, 57489/57409, stall | 5.9e-10, 58772/58709, stall | 1.8e-14, 86119/62000, flat | 4.3e-11, 58768/58719, stall | 0.0e+00, 158/12, stall | 2.9e-04, 3/227797 |
| 0.5 | polish 2σ, 5·10⁴ evaluations | 8 | 2.10e-10 | optimum | 1 | 0.999989242094 | 1.1e-05, 28/0, stall | 3.8e-10, 758/720, stall | 1.5e-14, 27856/4000, flat | 8.7e-11, 778/734, stall | 1.5e-19, 95/4, stall | 1.0e-05, 3/29 |
| 0.75 | polish 2σ, 5·10⁴ evaluations | 10 | 6.81e-10 | block | 2.30523 | 0.999976205201 | 1.1e-05, 34530/34451, stall | 5.9e-11, 35202/35152, stall | 1.4e-15, 67613/38000, flat | 1.6e-11, 35198/35155, stall | 0.0e+00, 66/3, stall | 3.9e-05, 3/94801 |
| 0.75 | polish 2σ, 5·10⁴ evaluations | 10 | 6.81e-10 | optimum | 1 | 0.999976205201 | 4.2e-06, 146/95, stall | 2.7e-11, 633/593, stall | 4.5e-15, 39827/4000, flat | 2.0e-13, 660/603, stall | 1.8e-20, 475/24, stall | 4.2e-06, 3/195 |
| 1 | polish 2σ, 5·10⁴ evaluations | 12 | 1.61e-09 | block | 1.46508 | 0.999996043806 | 1.2e-06, 19773/19726, stall | 3.7e-14, 20089/20037, stall | 1.6e-16, 75839/24000, flat | 8.4e-16, 20084/20039, stall | 0.0e+00, 399/19, stall | 3.2e-06, 3/51507 |
| 1 | polish 2σ, 5·10⁴ evaluations | 12 | 1.61e-09 | optimum | 1 | 0.999996043806 | 1.7e-06, 130/87, stall | 2.0e-12, 456/405, stall | 1.1e-16, 51176/4000, flat | 5.3e-13, 439/404, stall | 2.2e-20, 302/15, stall | 2.3e-06, 3/137 |
| 0.25 | polish 2σ, 5·10⁴, repeated, no gradient test | 6 | 2.55e-11† | block | 4.951 | 0.999969127454 | 1.9e-03, 64132/64022, stall | 6.6e-09, 67882/67837, stall | 1.2e-12, 79655/70000, flat | 6.4e-09, 67895/67844, stall | 2.1e-19, 1061/69, stall | 1.5e-02, 3/256857 |
| 0.25 | polish 2σ, 5·10⁴, repeated, no gradient test | 6 | 2.55e-11† | optimum | 1 | 0.999969127454 | 7.0e-07, 91/37, stall | 1.4e-12, 95/51, stall | 1.3e-12, 22046/4000, flat | 2.6e-12, 95/50, stall | 0.0e+00, 1860/149, stall | 1.3e-05, 3/44 |
| 0.75 | polish 2σ, 5·10⁴, repeated, no gradient test | 11 | 7.75e-11† | block | 2.6853 | 0.998373504731 | 2.3e-03, 135961/135766, stall | 1.4e-08, 151823/151738, stall | 3.6e-13, 163853/154000, flat | 1.8e-09, 151940/151880, stall | 0.0e+00, 260/16, stall | 5.3e-03, 3/384829 |
| 0.75 | polish 2σ, 5·10⁴, repeated, no gradient test | 11 | 7.75e-11† | optimum | 1 | 0.998373504731 | 2.3e-05, 961/894, stall | 2.2e-08, 1429/1356, stall | 9.4e-14, 14936/4000, flat | 9.8e-09, 1481/1422, stall | 3.5e-21, 326/19, stall | 7.7e-05, 3/1392 |

The long-double rms by direct subtraction (`lm_rms`, B4's instrument until round 81) against the exact one, minus 1: -2.4e-10 … +1.7e-10 over every end, -2.4e-10 … -2.3e-12 at the long-double solve's ends (12 of 12 below it: that solve minimised the same evaluation).

Each solve's excess over the best end, minus 1, from both starts: double +1.9e-08 … +2.3e-03; long double +2.4e-22 … +2.2e-08; double, anchored +9.5e-17 … +1.3e-12; double, anchored, no stall re-anchoring +2.4e-16 … +9.8e-09; long double, anchored, from the anchored double's end +0.0e+00 … +2.1e-19; TRF +2.4e-08 … +1.5e-02. By the optimum, from the path's blocks — 1.61e-09 (τ/σ = 0.25, K = 5): double 1.9e-08, long double 2.4e-22; 1.61e-09 (τ/σ = 1, K = 12): double 1.2e-06, long double 3.7e-14; 6.81e-10 (τ/σ = 0.75, K = 10): double 1.1e-05, long double 5.9e-11; 2.1e-10 (τ/σ = 0.5, K = 8): double 5.6e-05, long double 5.9e-10; 7.75e-11 (τ/σ = 0.75, K = 11): double 2.3e-03, long double 1.4e-08; 2.55e-11 (τ/σ = 0.25, K = 6): double 1.9e-03, long double 6.6e-09.

The best end is reached by long double, anchored, from the anchored double's end at 6 of the 6 blocks; that solve's ends from the two starts agree within 2.1e-19 (exactly), and a Gauss–Newton step from the best end's exact residual predicts at most 2.1e-19 of further decrease (`b4_gap`: ‖P_Jr‖²/(2‖r‖²), J and r exact). The 55 of the 72 ends whose excess is at least 10³ times that lie at their own Gauss–Newton gap from it: their excess over their gap is 1.000–1.495 (median 1.122) — the best end the minimum the solves approach, each excess a first-order distance to it.

The direct subtraction's floor in double over its floor in long double, from the same start (the two ends' excesses over the best end), where the long double's exceeds that spread (11 of 12 rows): 1.1e+03, 2.8e+04, 9.6e+04, 1.5e+05, 1.7e+05, 1.8e+05, 3e+05, 5e+05, 8.5e+05, 3.3e+07, 2.2e+08 — median 1.83e+05, as (ε_d/ε_ld)^a with a = 0.91–2.52 (median 1.59); a floor linear in the rounding (a = 1) gives ε_d/ε_ld = 2048; at the other 1.61e-09 (τ/σ = 0.25, K = 5) from the block: 2.4e-22.

The anchored double's floor, the rise of the exact rms that a rounding e₀ frozen in its residual predicts, ‖P_Je₀‖²/(2‖r‖²) at its end (`b4_floors`, exactly): with e₀ at the end 1.8e-16 … 8.8e-13; with e₀ at its last anchor, the rounding it froze, 1.4e-16 … 1.5e-12; its gap — the decrease a Gauss–Newton step predicts there — 9.5e-17 … 1.3e-12: the last anchor's prediction 0.50–2.53 times the gap, the end's 0.13–4.43; 2σ² rounded to double (relative 4.66e-17), the same at every anchor, alone predicts 0.10–0.51 of the last anchor's. Its last anchor: 12 stall, 1–28 iterations before the end (9.2e-11 … 4.5e-09 from it at the flat stops); per solve 19–1622 re-anchors at stalls and 0–76 by period; without re-anchoring at stalls, from the blocks, 3–75 by period.

The long double's anchored end of the best end's start with its parameters rounded to double, over the best end, minus 1 (a scale, not the double parameters' floor: B5; the anchored double's excesses in brackets): 1.61e-09 (τ/σ = 0.25, K = 5): 6.1e-17 (2.6e-16, 9.5e-17); 1.61e-09 (τ/σ = 1, K = 12): 1.7e-16 (1.6e-16, 1.1e-16); 6.81e-10 (τ/σ = 0.75, K = 10): 2.8e-16 (1.4e-15, 4.5e-15); 2.1e-10 (τ/σ = 0.5, K = 8): 1.8e-16 (1.8e-14, 1.5e-14); 7.75e-11 (τ/σ = 0.75, K = 11): 3.1e-14 (3.6e-13, 9.4e-14); 2.55e-11 (τ/σ = 0.25, K = 6): 1.1e-12 (1.2e-12, 1.3e-12).

The anchored double's ends from the block and from the stored optimum agree within 2.7e-13; its stops: flat 10, stall 2, after 8011 … 163853 iterations; without re-anchoring at stalls it stops after 46 … 151940. The long double's anchored solve from the anchored double's end gains on it, minus 1, -1.3e-12 … -9.5e-17. The stored optimum over the best end, minus 1: 1.61e-09 (τ/σ = 0.25, K = 5): +2.4e-08; 1.61e-09 (τ/σ = 1, K = 12): +4.0e-06; 6.81e-10 (τ/σ = 0.75, K = 10): +2.4e-05; 2.1e-10 (τ/σ = 0.5, K = 8): +1.1e-05; 7.75e-11 (τ/σ = 0.75, K = 11): +1.6e-03; 2.55e-11 (τ/σ = 0.25, K = 6): +3.1e-05.

The anchored double's error after 100000, 200000, 300000 iterations, over the best end, minus 1: 7.75e-11 (τ/σ = 0.75, K = 11) from the block: +8.7e-02.


### B5. The floors' order: the double's rounding scaled (referee round 84)

B4's solve in double at three blocks with the residual's rounding scaled by κ — r_κ = r_ld + κ(r_d − r_ld), formed in long double from the direct subtractions in double (r_d) and in long double (r_ld) at the same double parameters and rounded to double (`kappa_blend`; κ = 1 the double's residual exactly, κ = 0 the long double's) —, by direct subtraction from five starts (the stored optimum's, θ = 0; the path's block, θ = 1; the centres c_opt + θ(c_block − c_opt) with least-squares amplitudes between) and, from the stored optimum's, about an anchor whose r₀ is r_κ (B4's re-anchoring and flat stop): each end's excess over B4's best end of its block, minus 1 (exact), and the slope of log(excess) against log κ (least squares) over κ in [1/16, 1] and [1/4, 4]:

| τ/σ | K | start θ (its error) | κ = 0 | κ = 1/16 | κ = 1/8 | κ = 1/4 | κ = 1/2 | κ = 1 | κ = 2 | κ = 4 | slope on [1/16, 1] | on [1/4, 4] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 0 (1) | 1.4e-12 | 6.1e-09 | 4.7e-08 | 1.0e-07 | 2.1e-07 | 1.7e-06 | 2.3e-06 | 3.0e-06 | 1.84 | 1.32 |
| 1 | 12 | 0.25 (1.037) | 5.5e-12 | 1.5e-08 | 2.8e-08 | 7.9e-08 | 4.5e-07 | 9.6e-07 | 3.7e-06 | 8.7e-06 | 1.59 | 1.66 |
| 1 | 12 | 0.5 (1.134) | 1.7e-13 | 3.9e-09 | 1.6e-08 | 5.5e-08 | 4.9e-07 | 9.7e-07 | 2.6e-06 | 1.1e-05 | 2.09 | 1.76 |
| 1 | 12 | 0.75 (1.278) | 1.3e-12 | 6.9e-09 | 1.7e-08 | 6.4e-08 | 2.3e-07 | 1.8e-06 | 2.3e-06 | 1.2e-05 | 1.98 | 1.84 |
| 1 | 12 | 1 (1.465) | 1.7e-12 | 6.7e-09 | 3.3e-08 | 2.1e-07 | 4.7e-07 | 1.2e-06 | 4.7e-06 | 8.3e-06 | 1.89 | 1.39 |
| 0.25 | 6 | 0 (1) | 1.6e-12 | 6.2e-09 | 2.7e-08 | 3.2e-08 | 3.4e-07 | 7.0e-07 | 3.1e-05* | 3.1e-05* | 1.73 | 2.64 |
| 0.25 | 6 | 0.25 (87.34) | 7.3e-09 | 4.2e-05 | 6.6e-05 | 2.3e-04 | 7.7e-04 | 1.9e-03 | 4.5e-03 | 1.6e-02 | 1.46 | 1.49 |
| 0.25 | 6 | 0.5 (78.93) | 1.3e-08 | 2.1e-05 | 5.5e-05 | 1.7e-04 | 7.7e-04 | 2.0e-03 | 6.6e-03 | 1.5e-02 | 1.70 | 1.61 |
| 0.25 | 6 | 0.75 (38.23) | 8.0e-09 | 4.4e-05 | 7.1e-05 | 1.6e-04 | 7.3e-04 | 2.0e-03 | 7.3e-03 | 1.6e-02 | 1.44 | 1.67 |
| 0.25 | 6 | 1 (4.951) | 6.4e-09 | 2.7e-05 | 6.2e-05 | 2.2e-04 | 5.1e-04 | 1.9e-03 | 6.1e-03 | 1.6e-02 | 1.53 | 1.59 |
| 0.75 | 11 | 0 (1) | 1.9e-08 | 4.0e-07 | 6.1e-07 | 3.9e-06 | 4.0e-06 | 2.3e-05 | 1.2e-04 | 2.7e-04 | 1.45 | 1.71 |
| 0.75 | 11 | 0.25 (4.995) | 1.8e-08 | 3.3e-05 | 7.5e-05 | 2.9e-04 | 8.6e-04 | 2.4e-03 | 5.9e-03 | 3.1e-02 | 1.59 | 1.62 |
| 0.75 | 11 | 0.5 (4.031) | 2.1e-08 | 2.8e-05 | 8.2e-05 | 3.7e-04 | 1.1e-03 | 2.4e-03 | 6.7e-03 | 2.0e-02 | 1.66 | 1.41 |
| 0.75 | 11 | 0.75 (2.566) | 2.2e-08 | 3.0e-05 | 1.3e-04 | 2.2e-04 | 7.0e-04 | 2.5e-03 | 6.6e-03 | 2.0e-02 | 1.52 | 1.63 |
| 0.75 | 11 | 1 (2.685) | 2.1e-08 | 2.6e-05 | 7.6e-05 | 3.0e-04 | 8.1e-04 | 2.3e-03 | 8.1e-03 | 2.1e-02 | 1.64 | 1.57 |

* at most two accepted steps: the end is the start's, censored (the slopes include it).

The geometric mean over the starts: τ/σ = 1, K = 12 (5 starts): 0: 1.2e-12, 1/16: 7.0e-09, 1/8: 2.6e-08, 1/4: 9.0e-08, 1/2: 3.5e-07, 1: 1.3e-06, 2: 3.0e-06, 4: 7.7e-06 — slopes 1.88 on [1/16, 1], 1.60 on [1/4, 4]; the starts' largest over least at one κ 1.89–32.9; τ/σ = 0.25, K = 6 (5 starts): 0: 1.5e-09, 1/16: 5.8e-06, 1/8: 1.3e-05, 1/4: 3.4e-05, 1/2: 1.5e-04, 1: 4.0e-04, 2: 2.1e-03, 4: 4.6e-03 — slopes 1.57 on [1/16, 1], 1.80 on [1/4, 4]; the starts' largest over least at one κ 233–8.15e+03; τ/σ = 0.75, K = 11 (5 starts): 0: 2.0e-08, 1/16: 1.2e-05, 1/8: 3.2e-05, 1/4: 1.2e-04, 1/2: 2.9e-04, 1: 9.6e-04, 2: 3.0e-03, 4: 9.3e-03 — slopes 1.57 on [1/16, 1], 1.59 on [1/4, 4]; the starts' largest over least at one κ 1.24–275.

Per start (15): the slopes 1.44–2.09 on [1/16, 1] (median 1.64), 1.32–2.64 on [1/4, 4] (median 1.62); the starts' largest over least at one κ at most 8.15e+03 (τ/σ = 0.25, K = 6, κ = 0).

| τ/σ | K | about an anchor, from the stored optimum: κ = 0 | κ = 1/4 | κ = 1/2 | κ = 1 | κ = 2 | κ = 4 | slope on [1/2, 4] | over the best end rounded to double, κ = 0, 1/4, 1/2, 1, 2, 4 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 4.6e-17 (559, stall) | 3.2e-17 (512, stall) | 1.1e-16 (682, stall) | 1.1e-16 (51176, flat) | 1.1e-15 (48338, flat) | 5.3e-15 (47374, flat) | 2.03 | 0.277, 0.192, 0.637, 0.668, 6.81, 31.8 |
| 0.25 | 6 | 1.0e-12 (180, stall) | 1.0e-13 (5915, stall) | 3.0e-13 (22811, flat) | 1.3e-12 (22046, flat) | 4.3e-12 (21936, flat) | 2.8e-11 (22404, flat) | 2.15 | 0.954, 0.0972, 0.275, 1.2, 4.01, 26.2 |
| 0.75 | 11 | 1.8e-13 (1959, stall) | 6.0e-13 (2250, stall) | 4.1e-14 (15449, flat) | 9.4e-14 (14936, flat) | 7.0e-13 (14865, flat) | 2.8e-12 (14659, flat) | 2.11 | 5.77, 19.7, 1.35, 3.09, 22.8, 91.3 |

The anchored solves' slopes on [1/2, 4] 2.03–2.15; their flat stops by κ: 0: 0 of 3, 1/4: 0 of 3, 1/2: 2 of 3, 1: 3 of 3, 2: 3 of 3, 4: 3 of 3.

At κ = 1 the ends are B4's (the double's from the stored optimum and the block, the anchored double's from the stored optimum) in 9 of 9: rms, iterations, accepted steps, stop and parameters.


The entries over the variants that record them, by the block they join — one event per width, starting grid, iteration and kernel whose duals agree within 10⁻⁴ — the same entrant at the same iteration, whether or not the variants' paths have parted (the entries over the variants in brackets; the duals of one event agree to 6.8e-06; referee rounds 65, 70) —: the entrant had the largest current rate |⟨rφ_j⟩|; its rate over the largest; its dual |u_j| at the block's start; whether it starts beyond the block's outermost node (how far) or between its nodes (their duals). The entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩| (`sliding_iss`'s rule), its dual charged during the earlier blocks.

| joining | events (entries) | largest rate | rate / largest: median (range) | dual: median (least) | beyond: entrants, median (range) | between: entrants (duals) |
|---|---|---|---|---|---|---|
| none (the first entry) | 12 (36) | 12 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |
| an optimal block | 107 (316) | 1 | 0.21 (0.001–1.00) | 0.963 (0.570) | 105, 1.17σ (0.11–2.75σ) | 2 (0.989–0.992) |
| a non-optimal block | 102 (116) | 36 | 0.78 (0.004–1.00) | 0.731 (0.020) | 83, 2.56σ (0.17–4.38σ) | 19 (0.237–0.999) |
| all | 215 (468) | 49 | 0.41 (0.001–1.00) | 0.914 (0.000) | 182, 1.61σ (0.11–4.38σ) | 21 (0.237–0.999) |

6 events join an optimal block in some variants and a non-optimal one in others: counted in both rows, once in 'all'.

### τ/σ = 0.25, base (events: 10+, 11+, 8+, 12+, 13+, 7+, 16+, 14+, 4+, 5+, 4−, 17+, 3+, 9+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 18.87 | 2.1430e-02 | 2.1430e-02 | 1.000000 | True | 0.0e+00 | 0.000 | 2.1431e-02 | — |
| 2 | 529.4 | 3.7502e-04 | 3.7502e-04 | 1.000000 | True | 4.2e-11 | 0.000 | 3.7504e-04 | 57.143 |
| 3 | 5.093e+04 | 6.2260e-06 | 6.2260e-06 | 1.000000 | True | 2.2e-09 | 0.000 | 6.2262e-06 | 60.235 |
| 4 | 2.389e+05 | 1.0087e-07 | 1.0087e-07 | 1.000000 | True | 3.9e-07 | 0.000 | 1.0087e-07 | 61.723 |
| 5 | 6.067e+07 | 1.7479e-08 | 1.6114e-09 | 10.847501 | True | 3.1e-01 | 0.314 | 1.6114e-09 | 5.771 |
| 6 | 6.766e+08 | 2.4073e-09 | 2.5507e-11 | 94.379669 | True | 6.5e-01 | 0.649 | 2.5507e-11 | 7.261 |
| 7 | 2.499e+10 | 1.2145e-07 | 4.0115e-13 | 302749.882294 | True | 3.8e+00 | 3.826 | 4.0116e-13 | 0.020 |
| 8 | 2.521e+10 | 3.5753e-10 | 6.2753e-15 | 56974.096023 | True | 1.3e+00 | 1.321 | 6.2929e-15 | 339.685 |
| 9 | 1.021e+11 | 2.8072e-10 | 1.7429e-16 | 1610643.059713 | True | 4.1e+00 | 4.100 | 1.6856e-16 | 1.274 |
| 9 | 2.267e+11 | 2.3861e-10 | 1.7429e-16 | 1369005.248164 | True | 3.3e+00 | 3.290 | 1.6856e-16 | — |
| 10 | 4.159e+11 | 1.7099e-10 | 2.2206e-16 | 770019.312199 | True | 4.4e+00 | 4.355 | 2.9816e-16 | 1.395 |
| 11 | 7.292e+11 | 1.3883e-10 | 1.6488e-16 | 842053.741460 | True | 4.8e+00 | 4.789 | 2.9420e-16 | 1.232 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4): min C/max C within ±1e-07, its peaks beyond the outermost nodes at (left, right) (2.02, 2.02)σ, (2.61, 2.61)σ, (3.08, 3.08)σ, (3.48, 3.48)σ.

### τ/σ = 0.5, base (events: 10+, 11+, 8+, 13+, 7+, 14+, 6+, 17+, 3+, 2+, 19+, 4+, 18+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 18.15 | 7.8689e-02 | 7.8689e-02 | 1.000000 | True | 0.0e+00 | 0.000 | 7.8750e-02 | — |
| 2 | 141 | 5.0651e-03 | 5.0651e-03 | 1.000000 | True | 6.2e-10 | 0.000 | 5.0678e-03 | 15.535 |
| 3 | 3729 | 3.0925e-04 | 3.0925e-04 | 1.000000 | True | 2.7e-09 | 0.001 | 3.0939e-04 | 16.379 |
| 4 | 9700 | 1.8425e-05 | 1.8425e-05 | 1.000000 | True | 1.4e-08 | 0.000 | 1.8433e-05 | 16.784 |
| 5 | 6.597e+05 | 1.0824e-06 | 1.0824e-06 | 1.000000 | True | 2.9e-08 | 0.000 | 1.0828e-06 | 17.023 |
| 6 | 1.07e+07 | 9.0368e-07 | 6.3003e-08 | 14.343496 | True | 2.0e+00 | 1.976 | 6.3030e-08 | 1.198 |
| 7 | 1.518e+07 | 3.0205e-07 | 3.6437e-09 | 82.896419 | True | 1.2e+00 | 1.228 | 3.6452e-09 | 2.992 |
| 8 | 3.371e+08 | 2.9130e-07 | 2.0972e-10 | 1388.997007 | True | 3.6e+00 | 3.623 | 2.0980e-10 | 1.037 |
| 9 | 3.763e+08 | 2.8760e-07 | 1.2031e-11 | 23904.519332 | True | 3.9e+00 | 3.949 | 1.2031e-11 | 1.013 |
| 10 | 7.815e+08 | 2.5457e-07 | 6.8789e-13 | 370068.553314 | True | 4.6e+00 | 4.599 | 6.8790e-13 | 1.130 |
| 11 | 9.679e+08 | 2.4758e-07 | 3.9238e-14 | 6309816.493600 | True | 4.8e+00 | 4.789 | 3.9246e-14 | 1.028 |
| 12 | 1.101e+09 | 2.4274e-07 | 2.2524e-15 | 107768417.576135 | True | 4.6e+00 | 4.639 | 2.2392e-15 | 1.020 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5): min C/max C within ±9e-10, its peaks beyond the outermost nodes at (left, right) (2.06, 2.06)σ, (2.47, 2.47)σ, (2.81, 2.81)σ, (3.09, 3.09)σ, (3.34, 3.34)σ.

### τ/σ = 0.75, base (events: 10+, 11+, 13+, 8+, 7+, 14+, 15+, 5+, 16+, 19+, 4+, 2+, 18+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 17.29 | 1.5617e-01 | 1.5617e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 1.5666e-01 | — |
| 2 | 69.35 | 2.0060e-02 | 2.0060e-02 | 1.000000 | True | 0.0e+00 | 0.003 | 2.0103e-02 | 7.785 |
| 3 | 927.4 | 2.4427e-03 | 2.4427e-03 | 1.000000 | True | 1.2e-09 | 0.003 | 2.4474e-03 | 8.212 |
| 4 | 1281 | 2.9023e-04 | 2.9023e-04 | 1.000000 | True | 1.4e-08 | 0.003 | 2.9075e-04 | 8.416 |
| 5 | 3.709e+04 | 3.3999e-05 | 3.3999e-05 | 1.000000 | True | 0.0e+00 | 0.003 | 3.4058e-05 | 8.537 |
| 6 | 4.95e+05 | 3.9461e-06 | 3.9461e-06 | 1.000000 | True | 1.8e-08 | 0.003 | 3.9529e-06 | 8.616 |
| 7 | 2.59e+06 | 4.5505e-07 | 4.5505e-07 | 1.000000 | True | 2.1e-07 | 0.002 | 4.5582e-07 | 8.672 |
| 8 | 6.861e+06 | 3.9327e-07 | 5.2224e-08 | 7.530460 | True | 1.8e+00 | 1.838 | 5.2311e-08 | 1.157 |
| 9 | 1.494e+07 | 2.9214e-07 | 5.9713e-09 | 48.924738 | True | 1.9e+00 | 1.851 | 5.9811e-09 | 1.346 |
| 10 | 8.66e+08 | 2.8715e-07 | 6.8078e-10 | 421.802203 | True | 4.0e+00 | 3.988 | 6.8186e-10 | 1.017 |
| 11 | 9.907e+08 | 2.8108e-07 | 7.7548e-11 | 3624.569070 | True | 3.8e+00 | 3.765 | 7.7548e-11 | 1.022 |
| 12 | 1.16e+09 | 2.1206e-07 | 8.8020e-12 | 24091.841122 | True | 3.5e+00 | 3.475 | 8.8021e-12 | 1.325 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7): min C/max C within ±4e-08, its peaks beyond the outermost nodes at (left, right) (2.13, 2.12)σ, (2.40, 2.40)σ, (2.62, 2.62)σ, (2.83, 2.83)σ, (3.01, 3.01)σ, (3.18, 3.18)σ, (3.34, 3.34)σ.

### τ/σ = 1, base (events: 10+, 11+, 13+, 7+, 14+, 6+, 15+, 5+, 16+, 4+, 17+, 3+, 2+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 16.51 | 2.3915e-01 | 2.3915e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 2.4090e-01 | — |
| 2 | 44.49 | 4.7497e-02 | 4.7497e-02 | 1.000000 | True | 0.0e+00 | 0.008 | 4.7745e-02 | 5.035 |
| 3 | 374 | 8.9365e-03 | 8.9365e-03 | 1.000000 | True | 2.9e-09 | 0.009 | 8.9780e-03 | 5.315 |
| 4 | 643.1 | 1.6402e-03 | 1.6402e-03 | 1.000000 | True | 0.0e+00 | 0.009 | 1.6474e-03 | 5.449 |
| 5 | 5788 | 2.9676e-04 | 2.9676e-04 | 1.000000 | True | 5.1e-09 | 0.008 | 2.9802e-04 | 5.527 |
| 6 | 5.893e+04 | 5.3197e-05 | 5.3197e-05 | 1.000000 | True | 1.5e-07 | 0.008 | 5.3418e-05 | 5.579 |
| 7 | 1.212e+05 | 9.4739e-06 | 9.4739e-06 | 1.000000 | True | 8.4e-08 | 0.007 | 9.5126e-06 | 5.615 |
| 8 | 4.83e+05 | 1.6791e-06 | 1.6791e-06 | 1.000000 | True | 0.0e+00 | 0.007 | 1.6859e-06 | 5.642 |
| 9 | 5.348e+06 | 2.9649e-07 | 2.9649e-07 | 1.000000 | True | 0.0e+00 | 0.007 | 2.9767e-07 | 5.663 |
| 10 | 5.839e+07 | 9.2441e-08 | 5.2197e-08 | 1.771002 | True | 2.7e-01 | 0.260 | 5.2404e-08 | 3.207 |
| 11 | 8.35e+07 | 1.6136e-08 | 9.1674e-09 | 1.760191 | True | 2.6e-01 | 0.255 | 9.2035e-09 | 5.729 |
| 12 | 5.188e+08 | 6.0433e-08 | 1.6069e-09 | 37.609199 | True | 1.5e+00 | 1.458 | 1.6131e-09 | 0.267 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9): min C/max C within ±1e-08, its peaks beyond the outermost nodes at (left, right) (2.21, 2.21)σ, (2.38, 2.38)σ, (2.53, 2.53)σ, (2.67, 2.67)σ, (2.80, 2.80)σ, (2.92, 2.92)σ, (3.04, 3.04)σ, (3.14, 3.14)σ, (3.25, 3.25)σ.

### τ/σ = 1.5, base (events: 10+, 11+, 7+, 14+, 6+, 15+, 5+, 16+, 3+, 18+, 2+, 19+, 1+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 15.37 | 3.8940e-01 | 3.8940e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 3.9732e-01 | — |
| 2 | 26.96 | 1.2942e-01 | 1.2942e-01 | 1.000000 | True | 0.0e+00 | 0.032 | 1.3149e-01 | 3.009 |
| 3 | 124.3 | 4.0725e-02 | 4.0725e-02 | 1.000000 | True | 0.0e+00 | 0.035 | 4.1300e-02 | 3.178 |
| 4 | 172.1 | 1.2492e-02 | 1.2492e-02 | 1.000000 | True | 0.0e+00 | 0.034 | 1.2657e-02 | 3.260 |
| 5 | 1003 | 3.7762e-03 | 3.7762e-03 | 1.000000 | True | 2.1e-08 | 0.033 | 3.8240e-03 | 3.308 |
| 6 | 3159 | 1.1307e-03 | 1.1307e-03 | 1.000000 | True | 2.5e-08 | 0.031 | 1.1446e-03 | 3.340 |
| 7 | 7097 | 3.3631e-04 | 3.3631e-04 | 1.000000 | True | 1.8e-08 | 0.030 | 3.4036e-04 | 3.362 |
| 8 | 1.718e+04 | 9.9542e-05 | 9.9542e-05 | 1.000000 | True | 9.0e-08 | 0.029 | 1.0072e-04 | 3.379 |
| 9 | 1.939e+05 | 2.9351e-05 | 2.9351e-05 | 1.000000 | True | 1.2e-07 | 0.028 | 2.9693e-05 | 3.391 |
| 10 | 2.05e+05 | 8.6284e-06 | 8.6284e-06 | 1.000000 | True | 2.6e-07 | 0.027 | 8.7280e-06 | 3.402 |
| 11 | 7.793e+05 | 2.5303e-06 | 2.5303e-06 | 1.000000 | True | 1.1e-07 | 0.026 | 2.5593e-06 | 3.410 |
| 12 | 3.246e+06 | 7.4053e-07 | 7.4053e-07 | 1.000000 | True | 1.3e-07 | 0.025 | 7.4893e-07 | 3.417 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12): min C/max C within ±5e-07, its peaks beyond the outermost nodes at (left, right) (2.38, 2.38)σ, (2.41, 2.41)σ, (2.47, 2.47)σ, (2.53, 2.53)σ, (2.60, 2.60)σ, (2.66, 2.66)σ, (2.72, 2.72)σ, (2.77, 2.77)σ, (2.84, 2.84)σ, (2.89, 2.89)σ, (2.94, 2.94)σ, (2.99, 2.99)σ.

### τ/σ = 2, base (events: 10+, 11+, 7+, 14+, 5+, 16+, 4+, 17+, 18+, 3+, 20+, 1+, 0+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 14.69 | 5.0462e-01 | 5.0462e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 5.2280e-01 | — |
| 2 | 20.93 | 2.2364e-01 | 2.2364e-01 | 1.000000 | True | 3.9e-09 | 0.069 | 2.3059e-01 | 2.256 |
| 3 | 65.67 | 9.4014e-02 | 9.4014e-02 | 1.000000 | True | 0.0e+00 | 0.079 | 9.6611e-02 | 2.379 |
| 4 | 69.8 | 3.8509e-02 | 3.8509e-02 | 1.000000 | True | 1.9e-09 | 0.079 | 3.9497e-02 | 2.441 |
| 5 | 426 | 1.5539e-02 | 1.5539e-02 | 1.000000 | True | 8.8e-09 | 0.077 | 1.5918e-02 | 2.478 |
| 6 | 430.4 | 6.2092e-03 | 6.2092e-03 | 1.000000 | True | 5.4e-09 | 0.074 | 6.3555e-03 | 2.502 |
| 7 | 1616 | 2.4644e-03 | 2.4644e-03 | 1.000000 | True | 3.9e-09 | 0.071 | 2.5208e-03 | 2.520 |
| 8 | 5276 | 9.7318e-04 | 9.7318e-04 | 1.000000 | True | 4.0e-08 | 0.068 | 9.9499e-04 | 2.532 |
| 9 | 8671 | 3.8282e-04 | 3.8282e-04 | 1.000000 | True | 9.9e-08 | 0.066 | 3.9125e-04 | 2.542 |
| 10 | 1.156e+04 | 1.5013e-04 | 1.5013e-04 | 1.000000 | True | 1.7e-07 | 0.064 | 1.5339e-04 | 2.550 |
| 11 | 1.214e+05 | 5.8730e-05 | 5.8730e-05 | 1.000000 | True | 2.5e-07 | 0.062 | 5.9987e-05 | 2.556 |
| 12 | 1.231e+05 | 2.2927e-05 | 2.2927e-05 | 1.000000 | True | 1.7e-07 | 0.060 | 2.3413e-05 | 2.562 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12): min C/max C within ±1e-06, its peaks beyond the outermost nodes at (left, right) (2.57, 2.57)σ, (2.51, 2.51)σ, (2.51, 2.51)σ, (2.53, 2.53)σ, (2.56, 2.56)σ, (2.58, 2.58)σ, (2.60, 2.60)σ, (2.63, 2.63)σ, (2.66, 2.66)σ, (2.70, 2.70)σ, (2.72, 2.72)σ, (2.75, 2.75)σ.

## C. Under noise

τ/σ = 1, 10 draws per level, the path to 11 entries: the discrepancy rule's stop (kernels), its prediction from (i) and the predicted error there, the stopped fit's relative error (median [range]), the kernels within 4σ* of the feature at the stop, the largest relative gap's kernel count (in how many draws it is the stop's).

| s | stop: kernels (median [range]) | predicted | error at the stop: median [range] | predicted | kernels at the feature | largest gap at |
|---|---|---|---|---|---|---|
| 0.0001 | 6 [6, 6] (10/10 stopped) | 6 | 8.712e-05 [7.949e-05, 1.053e-04] | 9.533e-05 | 6 [6, 6] | 5 [4, 9] (the stop in 0) |
| 0.001 | 4 [4, 5] (10/10 stopped) | 4 | 1.712e-03 [5.808e-04, 1.813e-03] | 1.769e-03 | 4 [4, 5] | 4 [2, 5] (the stop in 4) |
| 0.01 | 3 [3, 4] (10/10 stopped) | 3 | 1.079e-02 [3.814e-03, 1.274e-02] | 1.057e-02 | 3 [3, 4] | 2 [2, 4] (the stop in 1) |

At the stop (K kernels): the error over the oracle's s(2K/N)^{1/2} (the variance alone, the right support), over the bias–variance prediction (e_K²‖f*‖²_N + 2Ks²/N)^{1/2} — e_K the noiseless optimum's error with K kernels (B's base path at τ/σ = 1) — and over the discrepancy rule's bound s((1 + δ)² − 1 + 4K/N)^{1/2}; paths with a runaway kernel before the stop:

| s | error / oracle: median [range] | / bias–variance | / bound: max | runaway |
|---|---|---|---|---|
| 0.0001 | 1.10 [1.01, 1.33] | 0.915 [0.835, 1.106] | 0.483 | 0 |
| 0.001 | 2.66 [0.81, 2.81] | 0.990 [0.745, 1.348] | 0.871 | 0 |
| 0.01 | 1.93 [0.59, 2.28] | 1.024 [0.573, 1.209] | 0.627 | 0 |

The stop as a margin (`stop_law`): the first-term formula's root K*(s) of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/(s(2δ + δ²)^{1/2}))/log(1/t), whose ceiling is the predicted stop, and the probability of the stop at each K — the rms's expectation s(b_K² + 1 − 2K/N)^{1/2} and its scatter from the noise's chi-square and its cross term with the bias, b_K = e_K‖f*‖_N/s, e_K the half-squeeze rule's error, in the normal approximation (the rms falls along the path) —; the expected number of the draws stopping at the predicted K against the observed:

| s | K*(s) | P(stop at K) | at the predicted K | expected | observed |
|---|---|---|---|---|---|
| 0.0001 | 5.26 | 5: 0.005, 6: 0.993, 7: 0.001 | 0.993 | 9.9 of 10 | 10 |
| 0.001 | 3.91 | 4: 0.789, 5: 0.210, 6: 0.001 | 0.789 | 7.9 of 10 | 7 |
| 0.01 | 2.55 | 3: 0.986, 4: 0.012 | 0.986 | 9.9 of 10 | 9 |

## D. The tail of the exact series (the sharp constant: proved by D‴.61 (iv) at every t < 1 — every σ*/σ —, by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's Appendix A.0 assumed correct, now independent cross-checks; before VA for t < 0.9900 with V, for t < 0.8052 with the branch and bound alone)

The half-squeeze rule's squared error over the series' first term, Σ_j t^{2j}(2K)!/(2K + 2j)!·(D_{2K+2j,K}/D_{2K,K})² (the defects exact, by the remainders modulo He_K), against ₂F₁(−½, −½; 1; t²), the sum of the terms' limits binom(1/2, j)²t^{2j}; K times the gap in parentheses:

| t | ₂F₁ | K = 10 | K = 20 | K = 40 | K = 80 | K = 160 | K = 320 |
|---|---|---|---|---|---|---|---|
| 0.1716 | 1.007373 | 1.005184 (-0.022) | 1.006173 (-0.024) | 1.006743 (-0.025) | 1.007050 (-0.026) | 1.007209 (-0.026) | 1.007291 (-0.026) |
| 0.382 | 1.036820 | 1.026417 (-0.104) | 1.030739 (-0.122) | 1.033484 (-0.133) | 1.035063 (-0.141) | 1.035916 (-0.145) | 1.036361 (-0.147) |
| 0.5 | 1.063544 | 1.047067 (-0.165) | 1.053239 (-0.206) | 1.057594 (-0.238) | 1.060295 (-0.260) | 1.061834 (-0.274) | 1.062665 (-0.281) |
| 0.8 | 1.167810 | 1.162342 (-0.055) | 1.156639 (-0.223) | 1.157221 (-0.424) | 1.159788 (-0.642) | 1.162436 (-0.860) | 1.164502 (-1.058) |

The ratios D_{2K+2j,K}/D_{2K,K} are polynomials in K of degree j (the (j + 1)-th divided difference over K = 1 … j + 2: 0, 0, 0, 0, 0, 0, 0) with the leading coefficients -1, -1/2, -1/2, -5/8, -7/8, -21/16, -33/16 at j = 1 … 7 — −C_{j−1}/2^{j−1} at every j (C the Catalan numbers), binom(1/2, j)(−2)^j.

Numbers are produced by `experiments/exp_width_mismatch.py`; nothing is edited by hand.
