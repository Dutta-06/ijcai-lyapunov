# Proposition D‴.61: the Gauss–Hermite Peano kernel (`exp_gh_peano.py`)

The K-point Gauss–Hermite rule of N(0, 1): ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] = (K!/(2K)!)ξ^{2K}e^{−ξ²/2}ĝ_K(ξ), ĝ_K the Laplace transform of the normalised Peano kernel g_K; conjecture P: ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} for complex |η| < 2.

## P. The defects as weighted paths: the ratios P_j(K) = D(2K + 2j)/D(2K)

D(m) = (E − Q_K)[z^m] as the weighted Dyck paths of length m reaching height K (`defect`; D(2K) = K! on K ≤ 20: yes), each P_j the exact polynomial through K = 1 … 2j + 1, checked on K ≤ 20:

| j | P_j(K) | off the polynomial on K ≤ the check |
|---|---|---|
| 1 | 1 + 1·K + 1·K^2 | 0 |
| 2 | 3 + 13/2·K + 3·K^2 + 2·K^3 + 1/2·K^4 | 0 |
| 3 | 15 + 67/2·K + 223/6·K^2 + 25/2·K^3 + 31/6·K^4 + 3/2·K^5 + 1/6·K^6 | 0 |
| 4 | 105 + 1163/4·K + 6439/24·K^2 + 605/3·K^3 + 347/6·K^4 + 197/12·K^5 + 13/3·K^6 + 2/3·K^7 + 1/24·K^8 | 0 |
| 5 | 945 + 10551/4·K + 139207/40·K^2 + 46025/24·K^3 + 8501/8·K^4 + 3305/12·K^5 + 3629/60·K^6 + 161/12·K^7 + 53/24·K^8 + 5/24·K^9 + 1/120·K^10 | 0 |
| 6 | 10395 + 130551/4·K + 193434/5·K^2 + 7976737/240·K^3 + 4694753/360·K^4 + 132535/24·K^5 + 948773/720·K^6 + 14689/60·K^7 + 5423/120·K^8 + 115/16·K^9 + 283/360·K^10 + 1/20·K^11 + 1/720·K^12 | 0 |

## M. The normalised Peano kernel's cumulants (exact)

g_K = G_K(2K)!/K! (a probability density): its variance exactly, (K² + K + 1)/((K + 1)(2K + 1)) by P_1, and κ_{2j}·K^j against c_j = −(2j − 1)!/4^j (the coefficients of ½log(1 − η²/4): -0.25, -0.375, -1.875, -19.6875, -354.375, -9745.31), with K(κ_{2j}K^j − c_j) in brackets:

| K | variance | K(κ_2 − ½) | κ_4·K^2 | κ_6·K^3 | κ_8·K^4 | κ_10·K^5 | κ_12·K^6 |
|---|---|---|---|---|---|---|---|
| 10 | 0.4805194805 | -0.194805 [+0.5519] | 0.0696649 [+4.447] | 2.51334 [+43.88] | 28.4423 [+481.3] | 215.287 [+5697] | 487.119 [+1.023e+05] |
| 20 | 0.4889663182 | -0.220674 [+0.5865] | -0.107048 [+5.359] | 1.46003 [+66.7] | 32.6903 [+1048] | 590.901 [+1.891e+04] | 9115.61 [+3.772e+05] |
| 40 | 0.4941282746 | -0.234869 [+0.6052] | -0.227061 [+5.918] | 0.220595 [+83.82] | 20.9686 [+1626] | 641.847 [+3.985e+04] | 18990 [+1.149e+06] |
| 80 | 0.4969710912 | -0.242313 [+0.615] | -0.297147 [+6.228] | -0.694112 [+94.47] | 5.98515 [+2054] | 391.966 [+5.971e+04] | 17541.6 [+2.183e+06] |
| 160 | 0.4984617171 | -0.246125 [+0.62] | -0.335048 [+6.392] | -1.24729 [+100.4] | -5.21121 [+2316] | 105.831 [+7.363e+04] | 9347.39 [+3.055e+06] |
| 320 | 0.4992248288 | -0.248055 [+0.6225] | -0.354761 [+6.477] | -1.55128 [+103.6] | -11.9939 [+2462] | -98.3738 [+8.192e+04] | 1585.21 [+3.626e+06] |
| 640 | 0.4996108978 | -0.249025 [+0.6237] | -0.364814 [+6.519] | -1.7106 [+105.2] | -15.7206 [+2539] | -219.301 [+8.645e+04] | -3568.49 [+3.953e+06] |
| 1280 | 0.4998050686 | -0.249512 [+0.6244] | -0.36989 [+6.541] | -1.79215 [+106] | -17.6732 [+2578] | -284.989 [+8.881e+04] | -6519.86 [+4.129e+06] |
| 2560 | 0.4999024391 | -0.249756 [+0.6247] | -0.372441 [+6.552] | -1.83342 [+106.5] | -18.6726 [+2598] | -319.21 [+9.002e+04] | -8097.12 [+4.219e+06] |
| 5120 | 0.4999511957 | -0.249878 [+0.6248] | -0.373719 [+6.557] | -1.85417 [+106.7] | -19.1781 [+2608] | -336.673 [+9.064e+04] | -8912.19 [+4.266e+06] |

At K = 5120 the scaled cumulants lie within 4.9e-04, 3.4e-03, 1.1e-02, 2.6e-02, 5.0e-02, 8.5e-02 of c_1 … c_6 (relative); extrapolated from K = 2560 and 5120 (2s(2K) − s(K), the 1/K term removed) -0.25, -0.374998, -1.87492, -19.6836, -354.136, -9727.27, within 2.5e-07, 5.6e-06, 4.3e-05, 2.0e-04, 6.7e-04, 1.9e-03.

## L. The Laplace transform against conjecture P's limit

ĝ_K(ξ)e^{−ξ²/4} (D‴.58's `gen_ratio`) over (1 − s²/4)^{1/2} on the real axis and (1 + s²/4)^{1/2} on the imaginary, |ξ| = s√K:

| K | axis | s = 0.1 | s = 0.2 | s = 0.3 | s = 0.4 | s = 0.5 | s = 0.6 | s = 0.7 | s = 0.8 | s = 0.9 | s = 1 | s = 1.1 | s = 1.2 | s = 1.3 | s = 1.4 | s = 1.5 | s = 1.6 | s = 1.7 | s = 1.8 | s = 1.9 | largest deviation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | real | 1.0001 | 1.0006 | 1.0014 | 1.0027 | 1.0045 | 1.0070 | 1.0105 | 1.0155 | 1.0226 | 1.0326 | 1.0470 | 1.0680 | 1.0991 | 1.1462 | 1.2198 | 1.3394 | 1.5473 | 1.9524 | 2.9756 | 1.98e+00 |
| 20 | imag | 0.9999 | 0.9994 | 0.9988 | 0.9979 | 0.9970 | 0.9960 | 0.9950 | 0.9942 | 0.9935 | 0.9930 | 0.9926 | 0.9925 | 0.9925 | 0.9927 | 0.9930 | 0.9934 | 0.9938 | 0.9943 | 0.9948 | 7.50e-03 |
| 40 | real | 1.0001 | 1.0003 | 1.0007 | 1.0014 | 1.0023 | 1.0037 | 1.0056 | 1.0084 | 1.0123 | 1.0182 | 1.0268 | 1.0399 | 1.0602 | 1.0924 | 1.1455 | 1.2371 | 1.4063 | 1.7590 | 2.7160 | 1.72e+00 |
| 40 | imag | 0.9999 | 0.9997 | 0.9994 | 0.9989 | 0.9985 | 0.9980 | 0.9975 | 0.9971 | 0.9967 | 0.9965 | 0.9963 | 0.9963 | 0.9963 | 0.9964 | 0.9965 | 0.9967 | 0.9969 | 0.9972 | 0.9974 | 3.74e-03 |
| 80 | real | 1.0000 | 1.0002 | 1.0004 | 1.0007 | 1.0012 | 1.0019 | 1.0029 | 1.0044 | 1.0065 | 1.0097 | 1.0146 | 1.0222 | 1.0344 | 1.0549 | 1.0904 | 1.1555 | 1.2844 | 1.5748 | 2.4382 | 1.44e+00 |
| 80 | imag | 1.0000 | 0.9999 | 0.9997 | 0.9995 | 0.9992 | 0.9990 | 0.9987 | 0.9985 | 0.9984 | 0.9982 | 0.9982 | 0.9981 | 0.9981 | 0.9982 | 0.9983 | 0.9984 | 0.9985 | 0.9986 | 0.9987 | 1.87e-03 |
| 160 | real | 1.0000 | 1.0001 | 1.0002 | 1.0004 | 1.0006 | 1.0010 | 1.0015 | 1.0022 | 1.0033 | 1.0050 | 1.0076 | 1.0118 | 1.0188 | 1.0309 | 1.0529 | 1.0960 | 1.1876 | 1.4127 | 2.1569 | 1.16e+00 |
| 160 | imag | 1.0000 | 0.9999 | 0.9998 | 0.9997 | 0.9996 | 0.9995 | 0.9994 | 0.9993 | 0.9992 | 0.9991 | 0.9991 | 0.9991 | 0.9991 | 0.9991 | 0.9991 | 0.9992 | 0.9992 | 0.9993 | 0.9994 | 9.32e-04 |

The largest deviation by K (its successive ratios): real: 1.98e+00, 1.72e+00, 1.44e+00, 1.16e+00 (0.869, 0.838, 0.804); imag: 7.50e-03, 3.74e-03, 1.87e-03, 9.32e-04 (0.498, 0.499, 0.500); on s ≤ 1.5: real: 2.20e-01, 1.46e-01, 9.04e-02, 5.29e-02 (0.662, 0.621, 0.586); imag: 7.50e-03, 3.74e-03, 1.87e-03, 9.32e-04 (0.498, 0.499, 0.500).

## LC. The Laplace transform off the axes (referee round 90)

ĝ_K(ξ)e^{−ξ²/4} at ξ = s√K·e^{iθ} (`gen_ratio_c`, complex decimal arithmetic) over conjecture P's (1 − ξ²/(4K))^{1/2}: |ratio/limit − 1|:

| K | θ | s = 0.1 | s = 0.2 | s = 0.3 | s = 0.4 | s = 0.5 | s = 0.6 | s = 0.7 | s = 0.8 | s = 0.9 | s = 1 | s = 1.1 | s = 1.2 | s = 1.3 | s = 1.4 | s = 1.5 | s = 1.6 | s = 1.7 | s = 1.8 | s = 1.9 | largest |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 1π/8 | 1.5e-04 | 6.0e-04 | 1.4e-03 | 2.6e-03 | 4.2e-03 | 6.4e-03 | 9.4e-03 | 1.3e-02 | 1.8e-02 | 2.5e-02 | 3.4e-02 | 4.6e-02 | 6.2e-02 | 8.4e-02 | 1.1e-01 | 1.5e-01 | 2.0e-01 | 2.6e-01 | 3.4e-01 | 3.4e-01 |
| 20 | 2π/8 | 1.5e-04 | 5.9e-04 | 1.3e-03 | 2.3e-03 | 3.7e-03 | 5.3e-03 | 7.1e-03 | 9.3e-03 | 1.2e-02 | 1.4e-02 | 1.7e-02 | 2.0e-02 | 2.3e-02 | 2.5e-02 | 2.8e-02 | 3.0e-02 | 3.2e-02 | 3.3e-02 | 3.4e-02 | 3.4e-02 |
| 20 | 3π/8 | 1.5e-04 | 5.7e-04 | 1.3e-03 | 2.2e-03 | 3.2e-03 | 4.3e-03 | 5.5e-03 | 6.6e-03 | 7.6e-03 | 8.5e-03 | 9.2e-03 | 9.6e-03 | 9.8e-03 | 9.8e-03 | 9.6e-03 | 9.3e-03 | 8.8e-03 | 8.2e-03 | 7.6e-03 | 9.8e-03 |
| 40 | 1π/8 | 7.6e-05 | 3.1e-04 | 7.2e-04 | 1.3e-03 | 2.2e-03 | 3.4e-03 | 4.9e-03 | 7.0e-03 | 9.8e-03 | 1.3e-02 | 1.8e-02 | 2.5e-02 | 3.4e-02 | 4.6e-02 | 6.3e-02 | 8.6e-02 | 1.2e-01 | 1.7e-01 | 2.4e-01 | 2.4e-01 |
| 40 | 2π/8 | 7.6e-05 | 3.0e-04 | 6.8e-04 | 1.2e-03 | 1.9e-03 | 2.7e-03 | 3.7e-03 | 4.7e-03 | 5.9e-03 | 7.1e-03 | 8.4e-03 | 9.7e-03 | 1.1e-02 | 1.2e-02 | 1.3e-02 | 1.3e-02 | 1.3e-02 | 1.2e-02 | 1.3e-02 | 1.3e-02 |
| 40 | 3π/8 | 7.5e-05 | 3.0e-04 | 6.5e-04 | 1.1e-03 | 1.6e-03 | 2.2e-03 | 2.8e-03 | 3.3e-03 | 3.8e-03 | 4.2e-03 | 4.6e-03 | 4.8e-03 | 4.9e-03 | 4.8e-03 | 4.7e-03 | 4.6e-03 | 4.3e-03 | 4.1e-03 | 3.8e-03 | 4.9e-03 |
| 80 | 1π/8 | 3.9e-05 | 1.6e-04 | 3.7e-04 | 6.8e-04 | 1.1e-03 | 1.7e-03 | 2.5e-03 | 3.6e-03 | 5.0e-03 | 7.0e-03 | 9.5e-03 | 1.3e-02 | 1.7e-02 | 2.4e-02 | 3.2e-02 | 4.3e-02 | 6.0e-02 | 8.7e-02 | 1.4e-01 | 1.4e-01 |
| 80 | 2π/8 | 3.8e-05 | 1.5e-04 | 3.5e-04 | 6.1e-04 | 9.6e-04 | 1.4e-03 | 1.9e-03 | 2.4e-03 | 3.0e-03 | 3.6e-03 | 4.2e-03 | 4.8e-03 | 5.3e-03 | 5.7e-03 | 6.0e-03 | 6.2e-03 | 6.3e-03 | 6.2e-03 | 6.0e-03 | 6.3e-03 |
| 80 | 3π/8 | 3.8e-05 | 1.5e-04 | 3.3e-04 | 5.6e-04 | 8.3e-04 | 1.1e-03 | 1.4e-03 | 1.7e-03 | 1.9e-03 | 2.1e-03 | 2.3e-03 | 2.4e-03 | 2.4e-03 | 2.4e-03 | 2.4e-03 | 2.3e-03 | 2.2e-03 | 2.0e-03 | 1.9e-03 | 2.4e-03 |
| 160 | 1π/8 | 1.9e-05 | 7.9e-05 | 1.8e-04 | 3.4e-04 | 5.6e-04 | 8.7e-04 | 1.3e-03 | 1.8e-03 | 2.6e-03 | 3.5e-03 | 4.8e-03 | 6.5e-03 | 8.7e-03 | 1.2e-02 | 1.5e-02 | 2.0e-02 | 2.5e-02 | 3.0e-02 | 3.7e-02 | 3.7e-02 |
| 160 | 2π/8 | 1.9e-05 | 7.7e-05 | 1.7e-04 | 3.1e-04 | 4.8e-04 | 6.9e-04 | 9.3e-04 | 1.2e-03 | 1.5e-03 | 1.8e-03 | 2.1e-03 | 2.4e-03 | 2.6e-03 | 2.8e-03 | 3.0e-03 | 3.1e-03 | 3.1e-03 | 3.1e-03 | 3.0e-03 | 3.1e-03 |
| 160 | 3π/8 | 1.9e-05 | 7.6e-05 | 1.7e-04 | 2.8e-04 | 4.2e-04 | 5.6e-04 | 7.1e-04 | 8.4e-04 | 9.6e-04 | 1.1e-03 | 1.1e-03 | 1.2e-03 | 1.2e-03 | 1.2e-03 | 1.2e-03 | 1.1e-03 | 1.1e-03 | 1.0e-03 | 9.4e-04 | 1.2e-03 |

The largest deviation by K (its successive ratios): θ = 1π/8: 3.39e-01, 2.38e-01, 1.39e-01, 3.73e-02 (0.703, 0.582, 0.269); θ = 2π/8: 3.36e-02, 1.29e-02, 6.26e-03, 3.09e-03 (0.383, 0.487, 0.493); θ = 3π/8: 9.83e-03, 4.86e-03, 2.41e-03, 1.20e-03 (0.494, 0.497, 0.499); the largest |arg| 1.5e-01.

## Q. The corner ratio: its limit, (iii)(a)

c_K(u) = Q_K[e^{uz}He_{K−1}(z)]/(u^{K−1}e^{u²/2}) at u = μ√K (`corner_ratio`: the corner entry of e^{uJ_K} over the infinite Jacobi matrix's, complex decimal arithmetic) against its limit 1 − μ² (|μ| < 1; (iii)(a)): |c/(1 − μ²) − 1|:

| K | arg μ | \|μ\| = 0.1 | \|μ\| = 0.3 | \|μ\| = 0.5 | \|μ\| = 0.7 | \|μ\| = 0.9 |
|---|---|---|---|---|---|---|
| 20 | 0π/8 | 4.93e-04 | 5.80e-03 | 2.75e-02 | 1.24e-01 | 8.81e-01 |
| 20 | 2π/8 | 4.76e-04 | 4.29e-03 | 1.21e-02 | 2.47e-02 | 4.82e-02 |
| 20 | 4π/8 | 4.60e-04 | 3.16e-03 | 5.03e-03 | 4.04e-03 | 1.37e-03 |
| 40 | 0π/8 | 2.53e-04 | 3.06e-03 | 1.55e-02 | 7.84e-02 | 6.58e-01 |
| 40 | 2π/8 | 2.44e-04 | 2.19e-03 | 6.02e-03 | 1.12e-02 | 1.44e-02 |
| 40 | 4π/8 | 2.35e-04 | 1.58e-03 | 2.46e-03 | 1.95e-03 | 6.67e-04 |
| 80 | 0π/8 | 1.28e-04 | 1.57e-03 | 8.36e-03 | 4.74e-02 | 4.79e-01 |
| 80 | 2π/8 | 1.23e-04 | 1.11e-03 | 2.98e-03 | 5.23e-03 | 6.41e-03 |
| 80 | 4π/8 | 1.19e-04 | 7.91e-04 | 1.21e-03 | 9.59e-04 | 3.29e-04 |
| 160 | 0π/8 | 6.46e-05 | 8.00e-04 | 4.38e-03 | 2.72e-02 | 3.39e-01 |
| 160 | 2π/8 | 6.21e-05 | 5.56e-04 | 1.48e-03 | 2.54e-03 | 3.13e-03 |
| 160 | 4π/8 | 5.97e-05 | 3.95e-04 | 6.03e-04 | 4.76e-04 | 1.63e-04 |

The largest deviation by K (its successive ratios): arg μ = 0π/8: 8.81e-01, 6.58e-01, 4.79e-01, 3.39e-01 (0.747, 0.728, 0.708); arg μ = 2π/8: 4.82e-02, 1.44e-02, 6.41e-03, 3.13e-03 (0.300, 0.444, 0.489); arg μ = 4π/8: 5.03e-03, 2.46e-03, 1.21e-03, 6.03e-04 (0.488, 0.494, 0.497).

## D. The reduction: Duhamel's identity, the Laplace step and the domination condition (referee round 91)

Duhamel's identity in s: ĝ_K(ξ)e^{−ξ²/4} = K·C(2K, K)∫_0^1 s^{K−1}(1 − s)^K e^{ξ²(1 − 2s)²/4}c_K(ξs)ds, ξ = √Kη, by Gauss–Legendre quadrature in double (c_K by `corner_ratio`) over `gen_ratio_c`, minus 1:

| K | \|η\| | arg η | points | quadrature over the exact, minus 1 |
|---|---|---|---|---|
| 12 | 1.5 | 2π/8 | 56 | 3.1e-15 |
| 12 | 1.9 | 1π/8 | 56 | 3.1e-15 |
| 40 | 1.5 | 2π/8 | 140 | 8.2e-15 |
| 40 | 1.9 | 4π/8 | 140 | 7.8e-15 |
| 40 | 1.95 | 1π/8 | 140 | 8.7e-15 |

The Laplace step at s = ½ (the phase's stationary point for every η): |ĝ_K(ξ)e^{−ξ²/4}(1 − η²/4)^{1/2}/c_K(ξ/2) − 1| (`gen_ratio_c`, `corner_ratio`; in brackets the corner ratio's own |c_K(ξ/2)/(1 − η²/4) − 1|):

| K | arg η | \|η\| = 1 | \|η\| = 1.5 | \|η\| = 1.8 | \|η\| = 1.9 |
|---|---|---|---|---|---|
| 20 | 0π/8 | 5.0e-03 (2.7e-02) | 2.8e-02 (1.9e-01) | 3.8e-02 (8.8e-01) | 5.4e-02 (2.1e+00) |
| 20 | 1π/8 | 3.8e-03 (2.2e-02) | 1.8e-02 (9.5e-02) | 6.9e-02 (2.1e-01) | 1.1e-01 (2.6e-01) |
| 20 | 2π/8 | 2.3e-03 (1.2e-02) | 2.4e-03 (2.9e-02) | 1.7e-02 (4.8e-02) | 3.2e-02 (5.9e-02) |
| 20 | 3π/8 | 2.1e-03 (6.6e-03) | 4.0e-03 (7.4e-03) | 5.1e-03 (6.6e-03) | 5.2e-03 (6.4e-03) |
| 20 | 4π/8 | 2.0e-03 (5.0e-03) | 3.6e-03 (3.4e-03) | 4.3e-03 (1.4e-03) | 4.6e-03 (6.7e-04) |
| 40 | 0π/8 | 2.6e-03 (1.5e-02) | 2.0e-02 (1.2e-01) | 6.1e-02 (6.6e-01) | 1.8e-02 (1.7e+00) |
| 40 | 1π/8 | 1.9e-03 (1.2e-02) | 7.3e-03 (5.7e-02) | 4.4e-02 (1.4e-01) | 8.1e-02 (1.9e-01) |
| 40 | 2π/8 | 1.3e-03 (6.0e-03) | 2.6e-03 (1.3e-02) | 6.8e-03 (1.4e-02) | 9.6e-03 (1.3e-02) |
| 40 | 3π/8 | 1.1e-03 (3.2e-03) | 2.0e-03 (3.5e-03) | 2.4e-03 (3.2e-03) | 2.5e-03 (3.2e-03) |
| 40 | 4π/8 | 1.1e-03 (2.5e-03) | 1.8e-03 (1.7e-03) | 2.2e-03 (6.7e-04) | 2.3e-03 (3.3e-04) |
| 80 | 0π/8 | 1.3e-03 (8.4e-03) | 1.2e-02 (7.7e-02) | 6.5e-02 (4.8e-01) | 6.9e-02 (1.3e+00) |
| 80 | 1π/8 | 9.6e-04 (6.0e-03) | 1.6e-03 (3.0e-02) | 2.4e-02 (8.5e-02) | 5.8e-02 (1.3e-01) |
| 80 | 2π/8 | 6.9e-04 (3.0e-03) | 1.4e-03 (5.7e-03) | 1.8e-03 (6.4e-03) | 1.9e-03 (6.6e-03) |
| 80 | 3π/8 | 5.8e-04 (1.6e-03) | 1.0e-03 (1.7e-03) | 1.2e-03 (1.6e-03) | 1.3e-03 (1.6e-03) |
| 80 | 4π/8 | 5.5e-04 (1.2e-03) | 9.2e-04 (8.2e-04) | 1.1e-03 (3.3e-04) | 1.1e-03 (1.6e-04) |
| 160 | 0π/8 | 6.5e-04 (4.4e-03) | 6.4e-03 (4.6e-02) | 5.5e-02 (3.4e-01) | 9.7e-02 (9.7e-01) |
| 160 | 1π/8 | 4.9e-04 (3.1e-03) | 5.0e-04 (1.5e-02) | 1.1e-02 (3.9e-02) | 3.5e-02 (6.1e-02) |
| 160 | 2π/8 | 3.6e-04 (1.5e-03) | 7.0e-04 (2.8e-03) | 8.7e-04 (3.1e-03) | 9.0e-04 (3.2e-03) |
| 160 | 3π/8 | 3.0e-04 (7.9e-04) | 5.0e-04 (8.4e-04) | 6.0e-04 (7.8e-04) | 6.2e-04 (7.7e-04) |
| 160 | 4π/8 | 2.8e-04 (6.0e-04) | 4.6e-04 (4.1e-04) | 5.4e-04 (1.6e-04) | 5.7e-04 (8.0e-05) |

K times the step's deviation at K = 160 (|η| = 1, 1.5, 1.8, 1.9), by arg η: 0π/8: 0.103, 1.018, 8.811, 15.540; 1π/8: 0.078, 0.080, 1.710, 5.668; 2π/8: 0.057, 0.111, 0.139, 0.144; 3π/8: 0.047, 0.081, 0.096, 0.100; 4π/8: 0.044, 0.074, 0.087, 0.091; the ratio of the last two deviations (K⁻¹'s 80/160 = 0.5): 0π/8: 0.488, 0.523, 0.850, 1.398; 1π/8: 0.508, 0.306, 0.444, 0.611; 2π/8: 0.518, 0.498, 0.477, 0.470; 3π/8: 0.509, 0.501, 0.497, 0.497; 4π/8: 0.507, 0.502, 0.500, 0.500.

The domination condition — sup over s ∈ [1/|η|, 1) of g(s) + Λ⁺(ηs) below g(½), g(s) = log(s(1 − s)) + (Re η²/2)((1 − s)² + s²), Λ⁺(μ) = M(Re μ) − ½ − log|μ| − Re(μ²)/2 (`lam_plus`), the corner ratio's bound beyond Q's disc —, the supremum over s in closed form (`gap_closed`; referee round 93: with w = 1 − s, log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4, η = Re^{iφ}, p = Re η²): it holds at every one of 361 angles in [0, π/2] for |η| < 1.892735 (the least radius, at arg η = 0.345); its gap at |η| = 1.999 on the real axis -2.50e-10, on the imaginary -5.00e-07; violated at |η| = 1.95 for at least arg η ∈ [0.087, 0.615] and at 1.999 for at least arg η ∈ [0.004, 0.764]; the gap's largest value on the grid (120 radii from 1 to c times the angles) -9.14e-03 at |η| ≤ 1.7 (at arg η = 0.000), -3.29e-03 at |η| ≤ 1.78 (at arg η = 0.000), -9.71e-04 at |η| ≤ 1.85 (at arg η = 0.000).

A branch and bound on the closed form (`dom_bb`: on a box of (|η|, arg η) the least radius, the least Re η² and the largest |η|cos(arg η) bound the gap from above; floating point with slack 1e-12, every accepted box rechecked in ball arithmetic — `box_bound_ball`: the same bound with every parameter at its rigorous end, the top edge at the true π/2, the supremum's candidates as balls; referee round 108 —) verifies the condition on 1 < |η| ≤ 1.78, 0 ≤ arg η ≤ π/2 (4935 boxes, the largest bound -3.3e-07, in ball arithmetic -3.2813e-07, 0 unresolved); 1 < |η| ≤ 1.8927, 0 ≤ arg η ≤ π/2 (402664 boxes, the largest bound -1.9e-11, in ball arithmetic -1.9432e-11, 0 unresolved). The violated sector at |η| = 1.95 by bisection on the closed form (`dom_sector`): arg η ∈ (0.08605, 0.61756), and there the supremum's maximiser s has Re(ηs) between 1.0051 and 1.1000 (numerically: floating point at 20001 angles and the ends, the printed digits rounded outward; referee round 131) — beyond 1, where Λ⁺'s M(Re ηs) is its linear piece (`gap_arg`; referee rounds 95, 99).

(iv) needs the limit of (iii) on |η| ≤ c with c² > 4 log u/(u − 1), u = 1/r, r a frame's Laplace rate (`c_need`): R's 22: r ≤ 0.6345, |η| ≤ 1.7773 (τ/σ = 2, 0.5·v_hs); S's 12: r ≤ 0.5279, |η| ≤ 1.6904 (τ/σ = 2, half-squeeze); every frame with r < 0.8052 lies within the verified radius 1.8927 (c_need(0.8052) ≤ 1.89269820 in ball arithmetic, inside it: referee round 110) — at the half-squeeze rule (r = t) σ*/σ < (1 + r)/(1 − r) = 9.26; below r = 0.0966 (|η| ≤ 1, where the condition is void; c_need(0.0966) ≤ 0.99982566 < 1) no check is needed.

Beyond Q's disc on the imaginary axis, c_K(i√Ky) (real there; 4K + 80 digits): over 1 + y², minus 1, at K = 20, 40, 80, 160: y = 1.2: +2.3e-03, +1.1e-03, +5.5e-04, +2.7e-04; y = 1.5: +4.2e-03, +2.1e-03, +1.0e-03, +5.1e-04; y = 1.8: +4.5e-03, +1.0e-02, +1.1e-03, +6.0e-04; beyond y = 1.8589, where the first alias's rate (5/8)y² − log y − ½ − (3/2)log 2 (`alias_rate`; referee round 93) vanishes, c_K − (1 + y²) grows: Prony fits (`prony`, a two-term recurrence) over 13 consecutive K from K = 100, 150 — y = 1.9: the rate 0.07475, 0.07471 (the alias 0.07468), the frequency 2.1318, 2.1318; y = 2: the rate 0.26720, 0.26716 (the alias 0.26713), the frequency 2.2854, 2.2854; y = 2.2: the rate 0.69688, 0.69685 (the alias 0.69682), the frequency 2.5693, 2.5693; y = 2.5: the rate 1.45029, 1.45026 (the alias 1.45024), the frequency 2.9217, 2.9217 — (two-point slopes oscillate with the frequency and do not measure it).

The deformed Duhamel path (`path_ok`): with the bound F(σ) = log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + C(σ) (C = 0 inside the unit disc, Λ⁺ outside) and the ball B of radius rb = ½(1 − |η|/2) around η/2, the two seeds η/2 ∓ 1.2·rb·u on the saddle's steepest-descent line (`path_seeds`: u² ∝ −conj φ″, φ″ = 2 − 8/η²; within 2rb of η/2, so in the unit disc) lie in the components of {F < F(η/2)} ∖ B (8-connected on a grid) that contain 0 and η, those two components differ, and the chord between the seeds stays below F(η/2) but at η/2 — by Taylor's theorem, φ‴(η/2) = 0, where the ratio q = D²/((|η|/2 − D)⁴|φ″|) < 1, D = |b − η/2| (`chord_q`; referee round 99) —, at |η| = 1.9: 19 of 19 angles (grid step 0.004; the seeds' F − F(η/2) -9.8e-05 … -1.9e-03, the chord's ratio q at most 5.8e-03); |η| = 1.95: 19 of 19 angles (grid step 0.002; the seeds' F − F(η/2) -1.2e-05 … -4.6e-04, the chord's ratio q at most 2.6e-03); |η| = 1.98: 19 of 19 angles (grid step 0.001; the seeds' F − F(η/2) -7.3e-07 … -7.3e-05, the chord's ratio q at most 9.5e-04) — 19 angles from 0 to 0.5π, the sector of the straight segment's failure among them; the components a grid check, numerical (referee round 95; part V certifies a valley condition on the annulus: below). Were the condition to hold on the annulus 1.8927 < |η| ≤ 1.98 (the branch and bound covers |η| ≤ 1.8927), (iv) would cover every frame with r < 0.9607 (σ*/σ < 49.89 at the half-squeeze rule).

## V. The valley condition certified (referee round 99's recipe; in ball arithmetic: round 104)

In the saddle's coordinates σ = (η/2)(1 + z) the bound is F_η(σ) − F_η(η/2) = log|1 − z²| + Re(cz²) + C(σ), c = η²/4 (`cert_E`): on the straight segment inside the unit disc, z = x real with (|η|/2)(1 + x) ≤ 1, it is at most log(1 − x²) + x²·Re c < log(1 − x²) + x² < 0 for x ≠ 0 (|c| < 1), so the path from 0 through η/2 to b = η/2 + t·η/2 (0 < t < (2 − |η|)/|η|) needs no check, Laplace's point crossed along the real z-axis (Re(c − 1) < 0); from b to η one path per η-box — the box centre's descent flow of the bound from b at t = 0.9·(2 − R2)/R2, then straight to z = 1 from |1 − z| ≤ 0.02 — with an upper bound of the bound over each segment and the whole box (`cert_seg_ub`: the value at the midpoint and the box centre, plus first-order terms in z and in η from explicit gradient bounds, C's |∇_σC| ≤ h′(|σ|) exactly (referee round 104); in floating point with slack 1e-12), a segment bisected and a box halved until every bound is below −slack: 1.8927 ≤ |η| ≤ 1.98, 0 ≤ arg η ≤ π/2: 760 boxes verified (0 failed), 160111 segment bounds, the largest -2.5e-12, the smallest box side 7.6e-04, the end piece's bound -1.383; 1.98 ≤ |η| ≤ 1.99, 0 ≤ arg η ≤ π/2: 2011 boxes verified (0 failed), 946644 segment bounds, the largest -1.0e-12, the smallest box side 1.9e-04, the end piece's bound -1.357; 1.99 ≤ |η| ≤ 1.995, 0 ≤ arg η ≤ π/2: 10658 boxes verified (0 failed), 9581078 segment bounds, the largest -1.0e-12, the smallest box side 4.8e-05, the end piece's bound -1.344. Every accepted segment rechecked in ball arithmetic (`cert_seg_ball`: Arb at 120 bits, python-flint 0.9.0; the same formula on the exact binary inputs, every monotone term at its rigorous upper end): the largest rigorous bound -2.4767e-12, -1.0258e-12, -1.0004e-12 (by band), at most 2.9e-14, 2.6e-14, 3.1e-14 above the float one; each box's covering disc from its corners rigorously, the top edge at the true π/2 (the float radius short of it at 220, 495, 2584 boxes, the rigorous one used), b inside the unit disc and the path's end within the end piece's radius rigorously, the end piece's rigorous bound -1.382677, -1.357139, -1.344304. The bound is even in η ↦ η̄ and η ↦ −η (z ↦ z̄; σ ↦ −σ), so the quadrant covers the annulus. With the branch and bound on |η| ≤ 1.8927, the limit of (iii)(b) holds on |η| ≤ 1.995 (both certificates in ball arithmetic: referee round 108; since round 148 cross-checks of the hand proof), and (iv) for every frame with r < 0.9900 (c_need(0.9900) ≤ 1.99497695 in ball arithmetic, inside |η| ≤ 1.995: referee round 110) — at the half-squeeze rule σ*/σ < 199.00 (in exact arithmetic from r's floor: referee round 104).

## VC. The corner at η = 2: a blow-up limit (not yet a proof)

With ε = 2 − |η|, θ = arg η = εφ and the saddle's coordinate z = ε(x + iy), the centred bound is ε³Ψ_φ(x, y) + o(ε³), Ψ_φ = −(x² − y²) − 4φxy + (x − ½)₊[(y + φ)² − (x − ½)²/3] (`corner_psi`: the quadratic part's leading order and, outside the unit disc, C's squares differing at order ε³). Ψ_φ is convex in y; its floor y = min(2φx, φ) (`corner_floor`) gives −(1 + 4φ²)x² inside the disc and −x² − φ² − (x − ½)³/3 outside: negative but at the saddle for every φ ≥ 0. Along the floor at x = 0.05, 0.25, 0.5, 1, 2, 4, 8 the true bound (`cert_E`) over ε³Ψ_φ: ε = 0.01: φ = 0 0.998–1.053, φ = 0.5 0.988–1.052, φ = 1 0.971–1.050, φ = 1.44 0.952–1.048, φ = 2 0.913–1.043, φ = 5 0.496–1.104; ε = 0.001: φ = 0 1.000–1.005, φ = 0.5 0.999–1.005, φ = 1 0.997–1.005, φ = 1.44 0.995–1.005, φ = 2 0.991–1.004, φ = 5 0.949–1.010, the bound negative at every point. On the straight segment (y = 0) beyond the disc (x ≥ ½; at the saddle Ψ is 0) Ψ's maximum is -0.2500 at φ = 0 (the real axis' gap −¼ε³, part D's −¼(2 − η)³) and turns positive at φ_dom = 1.4397 (`corner_phi_dom`; the root of 8φ³ − 12φ² + 1); the measured sector of the straight segment near the real axis (`corner_sector`, the domination gap's first violated arg η; floating point, valid for ε ≥ 10⁻³) over ε: |η| = 1.99: 1.4842, |η| = 1.995: 1.4614, |η| = 1.999: 1.4440. Its first correction inside the disc (referee round 115): `cert_E`/ε³ = Ψ_φ + εΨ¹_φ + O(ε²), Ψ¹_φ(w) = Re((¼ − 2φ² − 2iφ)w²) − Re(w⁴)/2 (`corner_psi1`), so that on the floor Ψ¹_φ/Ψ_φ = −(¼ + 5φ² + 8φ⁴ − x²Re((1 + 2iφ)⁴)/2)/(1 + 4φ²) ≈ −2φ²; at x = 0.05 (E/ε³ − Ψ_φ)/(εΨ_φ) in ball arithmetic (Arb, 320 bits) at ε = 10⁻⁶, 10⁻⁸ against it: φ = 0: -0.24875, -0.24875 (-0.24875); φ = 0.5: -1.0025, -1.0025 (-1.0025); φ = 1: -2.65175, -2.65175 (-2.65175); φ = 2: -8.70875, -8.70875 (-8.70875); φ = 5: -50.6287, -50.6287 (-50.6287); φ = 10: -200.257, -200.257 (-200.257) (Ψ_φ in balls; the midpoints kept) — the expansion holds for θ = εφ ≪ √ε only. Far out, at ε = 0.001, the true bound on Ψ_φ's floor and its least over y from 0 to three times the floor's: φ = 20, x = 0.05: -8.07e-10, -2.22e-09; φ = 20, x = 0.25: -2.50e-08, -5.60e-08; φ = 20, x = 0.45: -1.05e-07, -1.86e-07; φ = 25, x = 0.05: +1.54e-09, -2.78e-09; φ = 25, x = 0.25: +2.69e-08, -6.99e-08; φ = 25, x = 0.45: +3.99e-08, -2.30e-07; φ = 30, x = 0.05: +7.16e-09, -3.21e-09; φ = 30, x = 0.25: +1.55e-07, -8.08e-08; φ = 30, x = 0.45: +4.16e-07, -2.65e-07 — the floor positive from φ = 25, the least negative throughout: there another limit governs. At θ = √εψ, z = εx + i√εy, `cert_E`/ε² → Ψ̃_ψ = (1 + 2ψ²)y² − 4ψxy − y⁴/2 + D₊² − P₊², D = x − ½ + y²/2, P = x − ½ − ψ²/2 − ψy (`corner_psi_mid`; the referee's): in ball arithmetic at ε = 10⁻⁴, 10⁻⁶, 10⁻⁸ (ψ, x, y) = (1, 1, 0.5): -0.890684, -0.890626, -0.890625 against -0.890625; (ψ, x, y) = (2, 0.3, 0.4): 0.467173, 0.467200, 0.467200 against 0.467200; (ψ, x, y) = (3, 5, 3): 31.513697, 31.500137, 31.500001 against 31.500000; (ψ, x, y) = (0.2, 3, 0.1), Re σ > 1: -0.007273, -0.005839, -0.005825 against -0.005825. Along y = min(2ψx/(1 + 2ψ²), ψ) (`corner_mid_path`) Ψ̃_ψ is negative at every grid point for 145 ψ in [0.02, 20] and 0.0001 ≤ x ≤ max(200, 2 + 4ψ²) (at most -0.0043·ψ²min(x, 1)²; inside the disc −4ψ²x²/(1 + 2ψ²) − y⁴/2, far out −ψ² − ψ⁴/2 exactly, to 8.9e-13). The limits match: Ψ_φ(x, φỹ)/φ² (φ → ∞) and Ψ̃_ψ(x, ψỹ)/ψ² (ψ → 0) share Ψ_∞ = ỹ² − 4xỹ + (x − ½)₊(ỹ + 1)², its floor ỹ = min(2x, 1) at −4x² and −1 (ε ≪ θ ≪ √ε), and Ψ̃_ψ(x, ŷ/ψ) (ψ → ∞) tends to 2ŷ² − 4xŷ + (x − ½)₊², its floor ŷ = x at −2x² + (x − ½)₊² (θ ≫ √ε, the fixed-θ limit). So the corner has four regimes, each limit with an explicit negative path, all four proved (Ψ̃_ψ's by cases, Lemma VC in D61's (iii)(b) Status: on y = sx, s = 2ψ/(1 + 2ψ²), P ≤ −ψ²/2 and D₊² < 2ψxy + y⁴/2; on y = ψ a convex quadratic in x negative at both ends, then −ψ² − ψ⁴/2; the grid above a check). The four paths are one family (referee round 117): Im z = min(κRe z, θ), κ = Im A/|Re A| ≈ 2θ/(ε + 2θ²), A = η²/4 − 1 — the exact quadratic floor capped at Im z = θ —, in Ψ̃'s coordinates y = min(sx, ψ) to leading order (the cap y = ψ exactly; at ε > 0 the floor's slope is Im A/|Re A|: referee round 120), tending to Ψ_φ's floor as ψ → 0 and to the fixed-θ floor as ψ → ∞; Ψ_φ's own floor, y = min(2ψx, ψ) there, gives Ψ̃ = 4ψ²x²(2ψ² − 1) − 8ψ⁴x⁴ inside the disc, negative near the saddle iff ψ < 1/√2: valid for θ < √(ε/2) (φ = 20 and 25 at ε = 0.001 are ψ = 0.63 and 0.79). Two forms suffice (referee round 120): in the path point's coordinates (u, q) = (ε/X, θ²/X) the four regimes are faces of one square — the ψ chart its rays q/u = ψ², the φ and (ε/θ, θ/√ε) charts the face q = 0, the matching region its corner u = q = 0 (E ≈ −εθ² − θ⁴/2 − X³/3), the (√ε/θ, θ) chart q → ∞, through p = X/θ² —, the junction from b carried by the forms' derivative where their limit vanishes (on the real segment the limits vanish, the margin there V's real-axis bound of order ε³), and part VA certifies the path along them (below); the saddle's degeneration φ″ → 0 is harmless on this route and no overlap with the straight segment is needed (referee rounds 115, 117).

## VA. The last annulus certified in ball arithmetic (referee rounds 120–122, 126, 127, 131)

Along part VC's path family — from b = (1 + ε/4)η/2 (a point of part V's chord, ε = 2 − |η|) straight up to the floor Im z = κRe z, κ = Im A/|Re A|, A = η²/4 − 1, along it to the cap Im z = arg η, along the cap to Re z = 0.05, then a connector and a descent polyline per η-box to |1 − z| ≤ 0.02 and part V's end piece — the bound is certified negative (`gh_lastann`, python-flint at 128 bits) by two scaled forms of E/X² — exact on the path up to two one-sided simplifications ((Re σ − 1)₊ for (|Re σ| − 1)₊ in A, −P₊² dropped in B: upper bounds) —, A in (u, q, ρ) = (ε/X, θ²/X, X) and B in (u, p, ϑ) = (ε/X, X/θ², θ²), each box by the mean value form f(ρ) ≤ f(0) + ρ·sup f′ with f(0) in closed form (Ψ̃'s limit) and sup f′ by forward differentiation in balls: the junction (form A) 40 boxes, 0 failed; the floor and the cap (A) 2011 boxes, 0 failed; the junction (B) 320 boxes (634 outside its domain: ε > ε₀ or Re z > 0.05), 0 failed; the floor up to arg η = π/2 (B) 769 boxes (1319 outside its domain: ε > ε₀ or Re z > 0.05), 0 failed; the outer piece 76 η-boxes and 5193 segments, 0 failed, the largest rigorous bound -1.5e-07; the end piece -1.3314; the joint (`handoff_check`, referee rounds 121, 122), from the pieces certified: at every η-box the connector finite, the inner path's end min(κδ₀, arg η) proved inside it in balls on 484 sub-boxes (and sampled at 1140 η: corners, edges, arg η = 0 and π/2, random), the certified pieces chaining the connector from end to end and the polyline through every vertex from its start δ₀ + i·y_c to its end inside |1 − z| < 0.02 (round 126), y_c inside the connector — 0 failures —, the η-boxes inside the sector with disjoint interiors and exactly its area (in rationals); the forms' value below the certified bound at 240, 18326, 3056, 7434 distinct seeded points and corners (at the expansion variable's upper end and at its lower end where positive; 240, 24770, 3056, 7434 point–branch evaluations) in every accepted box of the four inner pieces (`soundness`), 0 above it, the largest excess -1.5e-13. The forms agree with `cert_E` and with F_η's definition to 7.6e-23 at 1479 path points, the closed forms with the forms at ρ = 0 to 2.1e-59. So the valley condition holds on 1.995 ≤ |η| < 2 (with η ↦ η̄, −η the whole annulus) and, with V (to |η| ≤ 1.995) and D's branch and bound (|η| ≤ 1.8927), both in ball arithmetic and meeting VA at exactly its R1, the limit (iii)(b) holds on the whole open disc |η| < 2 (D's domination is the valley condition on the straight path, so Ω ⊆ Ω′) — computer-assisted (Arb's, python-flint's and CPython's float arithmetic and `math.nextafter` assumed correct), and since referee round 148 an independent cross-check: the valley condition on 1 < |η| < 2 is proved by hand (review/RESEARCH_note_annulus.md, Theorem A). Round 121's referee verified the certificate independently — every identity re-derived, the code read line by line, 444 η with 437 720 path points at 400 bits — and found one defect: in the η-box nearest η = 2 at arg η ≤ π/128, Re A's ball contains 0, κ's ball and the connector's float ends were NaN, and the connector was silently dropped (round 120's 75 η-boxes and 5126 segments, a hole between δ₀ + iY(η) and the polyline); fixed by the referee's patch (a NaN or infinite κ gives the connector [0, θ_hi]; the connector always certified) — that box now split in two, the lower at arg η ≤ 0.01227 with its connector [0, 0.01227] certified; the counts above are the fixed certificate's (a negative result: round 120's certificate was not a proof as submitted). Round 122's referee verified the adopted certificate independently again (1493 η with 1 619 161 path points at 1200 bits, every η-box's corners and edges, 0 defects) and found that the first joint and tiling checks trusted the connector `outer_box` reported and an area sum; they now rebuild the cover from the pieces certified (above).

Part VA again at ε₀ = 0.11 (`part_VA011`, referee round 131, D3dz-3: VA's code, every piece at this ε₀, nothing else changed), from R1 = 2 − ε₀ = 1.89: the junction (form A) 54 boxes, 0 failed; the floor and the cap (A) 2533 boxes, 0 failed; the junction (B) 413 boxes (549 outside its domain), 0 failed; the floor (B) 1189 boxes (1174 outside its domain), 0 failed; the outer piece 220 η-boxes and 16083 segments, 0 failed, the largest rigorous bound -4.2e-09; the end piece -1.3314; the joint proved on 1306 sub-boxes (sampled at 3300 η), the certified pieces chaining connector and polyline, 0 failures, the η-boxes tiling the sector (in rationals); the forms' value below the certified bound at 324, 22258, 3982, 11634 distinct points (324, 30132, 3982, 11634 evaluations), 0 above it; the forms against `cert_E` and F_η's definition to 1.1e-12 at 1425 path points — this ε₀'s bar 1e-10 (`validate_tol`: the forms' midpoints are truncated series, whose tails grow with ε₀; VA's bar 1e-20 is ε₀ = 0.005's) —, the closed forms to 2.1e-59. VA011 certifies the valley condition at every η with 2 − |η| ≤ ε₀ (ε₀ the double 0.11; its outer piece from the double 1.89, at most that), which starts inside the disc |η| ≤ 1.8927 of D's branch and bound in ball arithmetic (both in rationals): so D and VA011 alone cover 1 < |η| < 2 — with |η| ≤ 1, where the domination condition is void, and the reflections, the limit (iii)(b) on the whole open disc without V —; V an independent cross-check on 1.8927 ≤ |η| ≤ 1.995. Since referee round 148 D, V, VA and VA011 are all independent cross-checks: the valley condition on 1 < |η| < 2 is proved by hand (review/RESEARCH_note_annulus.md, Theorem A).

## S. The frame series against conjecture P's Laplace prediction

S_K(ρ, v′) exactly (`frame_series`, 2K·log₁₀(1/r) + 2K + 60 digits) over the prediction — per angle the radius by Laplace's method with |1 − ξ²/(4K)|, the angle by the trapezoid (`log_s_pred`) —, S_K over the closed form (`log_s_simple`: off the cusp r^{2K}(1 ∓ r)/(ρv′(2γ(γ′ − γ))^{1/2}), at it (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/(ρv′)), r = 1/(2 min(γx, γy)), and the last K's series against the same with 40 more digits (the digits they agree to); the half-squeeze frame at v′ = 2 and the Mercer frame are the half-squeeze and the Mercer frame rule's own frames, v′ = 1.6 and 2.5 at ρ = t synthetic:

| τ/σ | frame | ρ | v′ | axis (r) | over the prediction at K = 10, 20, 40, 80 | K(ratio − 1) at the last two | over the closed form at the last K | digits agreed |
|---|---|---|---|---|---|---|---|---|
| 0.5 | half-squeeze | 0.055728 | 2 | cusp (0.055728) | 0.99969, 0.99984, 0.99992, 0.99996 | -0.003, -0.003 | 1.00153 | 419 |
| 0.5 | half-squeeze | 0.055728 | 1.6 | imaginary (0.045085) | 0.99856, 0.99870, 0.99886, 0.99913 | -0.046, -0.070 | 1.09021 | 428 |
| 0.5 | half-squeeze | 0.055728 | 2.5 | real (0.070644) | 1.00200, 1.00220, 1.00217, 1.00178 | +0.087, +0.142 | 1.10577 | 397 |
| 0.5 | Mercer | 0.101021 | 1.10102 | imaginary (0.058258) | 0.99418, 0.99560, 0.99733, 0.99857 | -0.107, -0.114 | 1.01381 | 399 |
| 1 | half-squeeze | 0.171573 | 2 | cusp (0.171573) | 0.99708, 0.99844, 0.99919, 0.99959 | -0.032, -0.033 | 1.00124 | 340 |
| 1 | half-squeeze | 0.171573 | 1.6 | imaginary (0.142136) | 0.98836, 0.99135, 0.99472, 0.99724 | -0.211, -0.221 | 1.01270 | 349 |
| 1 | half-squeeze | 0.171573 | 2.5 | real (0.224077) | 1.02037, 1.02160, 1.01789, 1.01187 | +0.715, +0.949 | 1.05802 | 317 |
| 1 | Mercer | 0.267949 | 1.26795 | imaginary (0.188345) | 0.97649, 0.98723, 0.99358, 0.99680 | -0.257, -0.256 | 1.00157 | 327 |
| 2 | half-squeeze | 0.381966 | 2 | cusp (0.381966) | 0.98634, 0.99233, 0.99588, 0.99786 | -0.165, -0.171 | 0.99987 | 284 |
| 2 | half-squeeze | 0.381966 | 1.6 | imaginary (0.330847) | 0.96806, 0.98360, 0.99219, 0.99620 | -0.313, -0.304 | 1.00029 | 292 |
| 2 | half-squeeze | 0.381966 | 2.5 | real (0.527864) | 1.14137, 1.15246, 1.13228, 1.09843 | +5.291, +7.875 | 1.17622 | 261 |
| 2 | Mercer | 0.500000 | 1.5 | imaginary (0.428571) | 0.96567, 0.98406, 0.99236, 0.99625 | -0.306, -0.300 | 0.99887 | 275 |

At K = 80, by where the mass lies: cusp (3 frames) 0.99786–0.99996, K(ratio − 1) -0.171 … -0.003, over the closed form 0.99987–1.00153; imaginary (6 frames) 0.99620–0.99913, K(ratio − 1) -0.304 … -0.070, over the closed form 0.99887–1.09021; real (3 frames) 1.00178–1.09843, K(ratio − 1) +0.142 … +7.875, over the closed form 1.05802–1.17622; the last K's series and the same with 40 more digits agree to 261–428 digits.

## R. D‴.58's scaled rules: the frame constant and the sharp asymptotics

Each rule of A6 with its own Mehler frame (−1 < ρ < 1: ρ = 1 + 2σ²/τ² − 2σ²/v, ℓ² = 2σ²ρ/(1 − ρ²), v′ = v/ℓ², a = ρv′ = v(1 − ρ²)/(2σ²), b = ρ(v′ − 1); ρ = 0, the Fock scaling, as the limit): the exact squared L² error (`a6_error`, 80 digits) over (σ*/σ)(1 − ρ²)^{1/2}(v/τ²)·S_K (`frame_series_any`), minus 1, at K = 2, 3, and A6's exact e_K² over the constant times conjecture P's prediction of S_K — per angle (`log_s_pred_ab`) at K = 20, 40, 60, 80, 120, 160, and the closed forms (`log_s_simple_ab`) at the last —, by where the mass lies:

| τ/σ | rule | v/v_hs | ρ | a | axis | the constant: deviation | e_K² over the per-angle prediction | over the closed form |
|---|---|---|---|---|---|---|---|---|
| 0.25 | Fock | 0.999541 | +0.000000 | 0.030303 | real | +2.0e-49, +6.6e-49 | 1.00055, 1.00050, 1.00043, 1.00037, 1.00028, 1.00022 | 1.03173 |
| 0.25 | Mercer-zero | 0.971825 | -0.941125 | 0.003367 | real | +7.3e-40, +6.9e-40 | 1.00384, 1.00198, 1.00133, 1.00100, 1.00067, 1.00051 | 1.00088 |
| 0.25 | Mercer frame | 1.000433 | +0.029437 | 0.030304 | imaginary | +1.2e-49, +3.6e-47 | 0.99947, 0.99953, 0.99960, 0.99966, 0.99974, 0.99980 | 1.02997 |
| 0.25 | 1·v_hs | 1.000000 | +0.015155 | 0.030310 | cusp | +3.3e-49, +1.7e-51 | 0.99999, 0.99999, 1.00000, 1.00000, 1.00000, 1.00000 | 1.00078 |
| 0.5 | Fock | 0.993808 | +0.000000 | 0.111111 | real | +8.5e-45, +3.3e-46 | 1.00571, 1.00364, 1.00262, 1.00204, 1.00141, 1.00108 | 1.01004 |
| 0.5 | Mercer-zero | 0.912871 | -0.797959 | 0.037075 | real | +1.5e-40, +1.3e-40 | 1.01659, 1.00868, 1.00589, 1.00445, 1.00299, 1.00225 | 1.00260 |
| 0.5 | Mercer frame | 1.005090 | +0.101021 | 0.111226 | imaginary | +3.5e-45, +3.3e-46 | 0.99560, 0.99733, 0.99814, 0.99857, 0.99903, 0.99927 | 1.00654 |
| 0.5 | 1·v_hs | 1.000000 | +0.055728 | 0.111456 | cusp | +2.1e-46, +1.1e-48 | 0.99984, 0.99992, 0.99994, 0.99996, 0.99997, 0.99998 | 1.00076 |
| 1 | Fock | 0.942809 | +0.000000 | 0.333333 | real | +4.8e-43, +1.5e-43 | 1.03770, 1.02210, 1.01564, 1.01211, 1.00835, 1.00637 | 1.01082 |
| 1 | Mercer-zero | 0.816497 | -0.464102 | 0.226497 | real | +2.8e-41, +2.1e-41 | 1.07563, 1.04210, 1.02929, 1.02248, 1.01536, 1.01167 | 1.01236 |
| 1 | Mercer frame | 1.035276 | +0.267949 | 0.339746 | imaginary | +1.3e-43, +4.2e-41 | 0.98723, 0.99358, 0.99573, 0.99680, 0.99787, 0.99840 | 1.00075 |
| 1 | 0.8·v_hs | 0.800000 | -0.535534 | 0.201724 | real | +4.2e-41, +5.0e-41 | 1.08162, 1.04544, 1.03162, 1.02427, 1.01659, 1.01261 | 1.01310 |
| 1 | 1·v_hs | 1.000000 | +0.171573 | 0.343146 | cusp | +8.9e-44, +3.9e-44 | 0.99844, 0.99919, 0.99945, 0.99959, 0.99972, 0.99979 | 1.00062 |
| 1 | 1.25·v_hs | 1.250000 | +0.737258 | 0.201724 | imaginary | +4.1e-41, +7.0e-41 | 0.98563, 0.99279, 0.99519, 0.99639, 0.99759, 0.99819 | 0.99899 |
| 2 | Fock | 0.745356 | +0.000000 | 0.666667 | real | +1.8e-41, +3.5e-41 | 1.28551, 1.18687, 1.14088, 1.11371, 1.08260, 1.06510 | 1.07129 |
| 2 | Mercer-zero | 0.745356 | +0.000000 | 0.666667 | real | +1.8e-41, +3.5e-41 | 1.28551, 1.18687, 1.14088, 1.11371, 1.08260, 1.06510 | 1.07129 |
| 2 | Mercer frame | 1.118034 | +0.500000 | 0.750000 | imaginary | +1.9e-42, +3.5e-42 | 0.98406, 0.99236, 0.99497, 0.99625, 0.99751, 0.99814 | 0.99944 |
| 2 | 0.5·v_hs | 0.500000 | -0.736068 | 0.204915 | real | +6.5e-40, +6.0e-40 | 1.83295, 1.53654, 1.40673, 1.33117, 1.24478, 1.19578 | 1.19480 |
| 2 | 0.8·v_hs | 0.800000 | +0.102458 | 0.708030 | real | +1.6e-41, +1.0e-41 | 1.22261, 1.14838, 1.11250, 1.09102, 1.06623, 1.05221 | 1.06064 |
| 2 | 1·v_hs | 1.000000 | +0.381966 | 0.763932 | cusp | +2.7e-42, +7.4e-43 | 0.99233, 0.99588, 0.99718, 0.99786, 0.99855, 0.99890 | 0.99991 |
| 2 | 1.25·v_hs | 1.250000 | +0.605573 | 0.708030 | imaginary | +2.5e-41, +1.0e-41 | 0.98507, 0.99264, 0.99511, 0.99634, 0.99757, 0.99818 | 0.99922 |
| 2 | 2·v_hs | 2.000000 | +0.940983 | 0.204915 | imaginary | +7.1e-39, +6.2e-40 | 0.98693, 0.99345, 0.99563, 0.99672, 0.99781, 0.99836 | 0.99936 |

Of A6's 36 rules 22 have a frame (5 with ρ < 0, 5 with ρ = 0, 12 with ρ > 0) and 14 not (|ρ| ≥ 1: τ/σ = 0.25, moment matching (+1.000), τ/σ = 0.25, 0.5·v_hs (-32.970), τ/σ = 0.25, 0.8·v_hs (-8.231), τ/σ = 0.25, 1.25·v_hs (+6.612), τ/σ = 0.25, 2·v_hs (+16.508), τ/σ = 0.5, moment matching (+1.000), τ/σ = 0.5, 0.5·v_hs (-8.889), τ/σ = 0.5, 0.8·v_hs (-2.180), τ/σ = 0.5, 1.25·v_hs (+1.845), τ/σ = 0.5, 2·v_hs (+4.528), τ/σ = 1, moment matching (+1.000), τ/σ = 1, 0.5·v_hs (-2.657), τ/σ = 1, 2·v_hs (+1.586), τ/σ = 2, moment matching (+1.000)). The frame constant holds to 7.1e-39 at all 22 — the 80-digit rounding of a cancelling sum —; at K = 160 the exact e_K² lie over the per-angle prediction 0.99814–1.19578, 16 of the 22 within 0.01 (0.99814–1.00637); the others τ/σ = 1, Mercer-zero (1.01167), τ/σ = 1, 0.8·v_hs (1.01261), τ/σ = 2, Fock (1.06510), τ/σ = 2, Mercer-zero (1.06510), τ/σ = 2, 0.5·v_hs (1.19578), τ/σ = 2, 0.8·v_hs (1.05221); over the closed forms 0.99899–1.19480, 12 of the 22 within 0.01 (0.99899–1.00654); the others τ/σ = 0.25, Fock (1.03173), τ/σ = 0.25, Mercer frame (1.02997), τ/σ = 0.5, Fock (1.01004), τ/σ = 1, Fock (1.01082), τ/σ = 1, Mercer-zero (1.01236), τ/σ = 1, 0.8·v_hs (1.01310), τ/σ = 2, Fock (1.07129), τ/σ = 2, Mercer-zero (1.07129), τ/σ = 2, 0.5·v_hs (1.19480), τ/σ = 2, 0.8·v_hs (1.06064).

Numbers are produced by `experiments/exp_gh_peano.py`; nothing is edited by hand.
