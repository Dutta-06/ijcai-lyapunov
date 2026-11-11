# Lemma E.0: the explicit rate against the exact Floquet rate (`exp_explicit_rate.py`)

The explicit rate lies below the exact rate in every row. The ordering is numerical: Ψ_M, Ψ_D are grid maxima (80–100 samples per period), α a converged quadrature (checked by refinement), ε a finite grid; the margins are at least 17×.

## 1. Remark E.1's instance, x_m = 1.5 sin ωt (m = n = 1, G = I, k_e = 3)

| ω | window T | α(T) | Ψ_M | Ψ_D | exact Floquet rate | explicit rate, closed-form ε | explicit rate, optimised ε | exact / explicit |
|---|---|---|---|---|---|---|---|---|
| 1 | 6.28 | 1.13 | 1.78 | 3.69 | 4.660e-02 | 3.688e-05 | 4.886e-05 | 954 |
| 2 | 3.14 | 0.566 | 1.78 | 7.38 | 2.158e-02 | 7.858e-05 | 1.040e-04 | 208 |
| 4 | 1.57 | 0.283 | 1.78 | 14.8 | 6.834e-03 | 8.696e-05 | 1.152e-04 | 59.3 |
| 8 | 0.785 | 0.142 | 1.78 | 29.5 | 1.831e-03 | 4.736e-05 | 6.301e-05 | 29.1 |
| 16 | 0.393 | 0.0708 | 1.78 | 59 | 4.662e-04 | 1.650e-05 | 2.199e-05 | 21.2 |
| 32 | 0.196 | 0.0354 | 1.78 | 118 | 1.171e-04 | 4.720e-06 | 6.286e-06 | 18.6 |
| 64 | 0.0982 | 0.0177 | 1.78 | 236 | 2.931e-05 | 1.249e-06 | 1.665e-06 | 17.6 |

Log-log slopes over ω ≥ 4: exact -1.97, explicit -1.56 (Remark E.1 reports the exact ω⁻² law); explicit slope over the last octave (ω = 32 → 64): -1.92.

## 2. Slow-adaptation limit: the ideal of the Theorem E tests (m = 3, n = 2, k_e = 6), two-tone reference of period 20π, G = g·diag(2, 2, 1)

| g | Ψ_M | Ψ_D | α(T = 20π) | Ψ_M²/k_e | κ_Ψ = Ψ_M²/λ_min⟨Ψ_GᵀΨ_G⟩ | exact Floquet rate | explicit rate, closed-form ε | explicit rate, optimised ε | exact / explicit (closed form) | exact / explicit (optimised) | 3κ_Ψ(1 + Ψ_D/(k_eΨ_M))² | averaging value λ_min⟨Ψ_GᵀΨ_G⟩/k_e |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 9.08 | 52.3 | 0.0408 | 13.7 | 1.27e+05 | 7.465e-05 | 8.542e-16 | 1.138e-15 | 8.74e+10 | 6.56e+10 | 1.47e+06 | 1.079e-04 |
| 10 | 5.24 | 30.2 | 0.0136 | 4.58 | 1.27e+05 | 2.358e-05 | 2.539e-15 | 3.382e-15 | 9.29e+09 | 6.97e+09 | 1.47e+06 | 3.597e-05 |
| 3 | 2.87 | 16.5 | 0.00408 | 1.37 | 1.27e+05 | 7.130e-06 | 8.204e-15 | 1.094e-14 | 8.69e+08 | 6.52e+08 | 1.47e+06 | 1.079e-05 |
| 1 | 1.66 | 9.55 | 0.00136 | 0.458 | 1.27e+05 | 2.381e-06 | 2.256e-14 | 3.006e-14 | 1.06e+08 | 7.92e+07 | 1.47e+06 | 3.597e-06 |
| 0.3 | 0.908 | 5.23 | 0.000408 | 0.137 | 1.27e+05 | 7.132e-07 | 5.698e-14 | 7.597e-14 | 1.25e+07 | 9.39e+06 | 1.47e+06 | 1.079e-06 |
| 0.1 | 0.524 | 3.02 | 0.000136 | 0.0458 | 1.27e+05 | 2.375e-07 | 9.106e-14 | 1.213e-13 | 2.61e+06 | 1.96e+06 | 1.47e+06 | 3.597e-07 |
| 0.03 | 0.287 | 1.65 | 4.08e-05 | 0.0137 | 1.27e+05 | 7.124e-08 | 8.019e-14 | 1.069e-13 | 8.88e+05 | 6.66e+05 | 1.47e+06 | 1.079e-07 |
| 0.01 | 0.166 | 0.955 | 1.36e-05 | 0.00458 | 1.27e+05 | 2.374e-08 | 4.216e-14 | 5.622e-14 | 5.63e+05 | 4.22e+05 | 1.47e+06 | 3.597e-08 |

As g → 0 (slow adaptation) the ratio exact/explicit decreases towards a constant of order 3κ_Ψ(1 + Ψ_D/(k_eΨ_M))², where κ_Ψ is the condition number of the averaged Gramian (Remark E.0a): the exact rate pays for ill-conditioning once, through λ_min, and the bound pays twice.  The formula is the leading-order value of the closed-form rate kα²/(3c₁²); the measured ratio is smaller because the Young step and the bound p_hi ≤ 3/4 have slack and because ε is optimised.  The ratio grows further with the excitation strength Ψ_M²/k_e (fast adaptation).  α is the PE level over one period (t-independent; 800 Gauss–Legendre nodes, converged).

## 3. The PE window (same ideal, G = diag(3, 3, 1.5))

| T / period | α(T) | explicit rate | ε |
|---|---|---|---|
| 0.25 | 0.000115 | 1.336e-14 | 2.404e-10 |
| 0.5 | 0.000655 | 3.172e-14 | 9.859e-11 |
| 1 | 0.00204 | 2.094e-14 | 2.056e-11 |
| 2 | 0.00408 | 5.468e-15 | 2.654e-12 |

## 4. Lemma E.0′ in the slow regime (the ideal of table 2 at g = 0.01, T = 20π)

λ_min(M̄) = 2.164e-07, κ(M̄) = 7.78e+03, Ψ_M = 0.168, Ψ_D = 0.957, Ψ_DD = 8.39, TΨ_M² = 1.78; peak factors ρ = 69.9, ρ_D = 3721, ρ_DD = 3.67e+05; ‖N₂‖ = 2.9e-21 (the lemma: N₂ = 0); threshold k₀ of conditions (A)–(D): 226 (δ = 0.2), 321 (δ = 0.1).

| k | δ | conditions (A) c₁c₂ ≤ δ/2, (B), (C), (D) | β | c_ζ | guaranteed rate γ/2 | exact Floquet rate | exact / guaranteed | averaged rate λ_min/k | exact / averaged | Lemma E.0 optimised rate | exact / E.0 | sup V̇/V on the grid, −γ | overshoot C |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 0.2 | 2.066: no, yes, yes, yes | 727 | +0.569 | nan | 4.278e-09 | nan | 4.328e-09 | 0.9885 | 2.577e-14 | 1.66e+05 | +nan, +nan | nan |
| 50 | 0.1 | 2.066: no, yes, yes, yes | 1455 | +0.669 | nan | 4.278e-09 | nan | 4.328e-09 | 0.9885 | 2.577e-14 | 1.66e+05 | +nan, +nan | nan |
| 100 | 0.2 | 0.511: no, yes, yes, yes | 713 | +0.743 | nan | 2.157e-09 | nan | 2.164e-09 | 0.9971 | 1.475e-14 | 1.46e+05 | +nan, +nan | nan |
| 100 | 0.1 | 0.511: no, yes, yes, yes | 1426 | +0.843 | nan | 2.157e-09 | nan | 2.164e-09 | 0.9971 | 1.475e-14 | 1.46e+05 | +nan, +nan | nan |
| 200 | 0.2 | 0.127: no, yes, yes, yes | 706 | +0.786 | nan | 1.081e-09 | nan | 1.082e-09 | 0.9993 | 7.917e-15 | 1.37e+05 | +nan, +nan | nan |
| 200 | 0.1 | 0.127: no, yes, yes, yes | 1412 | +0.886 | nan | 1.081e-09 | nan | 1.082e-09 | 0.9993 | 7.917e-15 | 1.37e+05 | +nan, +nan | nan |
| 300 | 0.2 | 0.056: yes, yes, yes, yes | 704 | +0.794 | 5.708e-10 | 7.210e-10 | 1.263 | 7.213e-10 | 0.9997 | 5.402e-15 | 1.33e+05 | -1.436e-09, -1.142e-09 | 26.6 |
| 300 | 0.1 | 0.056: no, yes, yes, yes | 1407 | +0.894 | nan | 7.210e-10 | nan | 7.213e-10 | 0.9997 | 5.402e-15 | 1.33e+05 | +nan, +nan | nan |
| 400 | 0.2 | 0.032: yes, yes, yes, yes | 702 | +0.796 | 4.299e-10 | 5.409e-10 | 1.258 | 5.410e-10 | 0.9998 | 4.099e-15 | 1.32e+05 | -1.078e-09, -8.598e-10 | 26.6 |
| 400 | 0.1 | 0.032: yes, yes, yes, yes | 1405 | +0.896 | 4.839e-10 | 5.409e-10 | 1.118 | 5.410e-10 | 0.9998 | 4.099e-15 | 1.32e+05 | -1.079e-09, -9.677e-10 | 37.6 |
| 600 | 0.2 | 0.014: yes, yes, yes, yes | 701 | +0.798 | 2.875e-10 | 3.606e-10 | 1.254 | 3.606e-10 | 0.9999 | 2.771e-15 | 1.3e+05 | -7.187e-10, -5.750e-10 | 26.5 |
| 600 | 0.1 | 0.014: yes, yes, yes, yes | 1403 | +0.898 | 3.235e-10 | 3.606e-10 | 1.115 | 3.606e-10 | 0.9999 | 2.771e-15 | 1.3e+05 | -7.199e-10, -6.471e-10 | 37.5 |

Reading: above the threshold the guaranteed rate of Lemma E.0′ is the averaged rate λ_min(M̄)/k up to the factor c_ζ/(1 + εp₁ + ε²p₂) ≈ 1 − δ, within a factor 1.3 of the exact rate, where Lemma E.0's optimised rate is five orders of magnitude below it; the certificate sup V̇/V ≤ −γ is verified on the grid wherever the conditions hold. Below the threshold the conditions fail through (A): the derivative peak factor ρ_DD of this regressor is large (the directions that are poorly excited on average are excited in short bursts).

## 5. Corollary E.0′.2: Theorem E in the slow regime on the single-kernel instance (m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = 0.3·I, x_m = 1.5 sin t, k_e = 40)

λ_min(M̄_G) = 5.409e-02, κ(M̄_G) = 6.7, peak factors ρ = 4.8, ρ_D = 30.3, ρ_DD = 461; threshold k₀ = 24.4 (δ = 0.2), 33.5 (δ = 0.1); r = min(ρ, δ) = 0.350, L_J = 4.62 (estimated from samples); exact linear rate 1.348e-03.

| Lyapunov function | explicit ball radius r₀ (w-coordinates) | nonlinear rate λ | overshoot C | exact linear rate / λ |
|---|---|---|---|---|
| Lemma E.0, optimised ε (Theorem E) | 6.146e-06 | 2.628e-05 | 1.43 | 51.3 |
| Lemma E.0′, δ = 0.2 (Corollary E.0′.2) | 1.532e-06 | 4.533e-04 | 11.84 | 3.0 |
| Lemma E.0′, δ = 0.1 (Corollary E.0′.2) | 8.323e-07 | 5.159e-04 | 17.14 | 2.6 |

Reading: on this well-conditioned regressor (κ ≈ 7) Lemma E.0 is only an order of magnitude conservative, and the corollary trades the certified ball for the rate: its nonlinear rate is within a small factor of the exact linear rate (the factor 2 of the Lyapunov argument times 1/(1 − δ)), while its ball radius carries the anisotropy p_lo/p_hi ≈ 1/β of the quadratic form. On ill-conditioned regressors the rate improves by the factor of Remark E.0a (table 4); the radius comparison was not computed there.

## 6. Corollary E.0′.3 — the stationary parameter spread of the dithered (LTV), Σ = √(2·10⁻³)·I

| instance (κ(M̄)) | k | δ | exact period mean of E‖ζ‖² | exact max over the period | averaged k·tr X | explicit bound 2ν̄/(γ lo) | bound / exact | 1/(1 − δ) | ν̄ grid / explicit | exact E‖e‖² | bound on E‖w‖² |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ideal of table 4 (7.78e+03) | 300 | 0.2 | 1.5094e+06 | 1.5094e+06 | 1.5089e+06 | 2.6444e+07 | 17.519 | 1.250 | 0.997 | 5.000e-05 | 3.747e+04 |
| ideal of table 4 (7.78e+03) | 300 | 0.1 | 1.5094e+06 | 1.5094e+06 | — | — (conditions fail) | — | 1.111 | — | 5.000e-05 | — |
| ideal of table 4 (7.78e+03) | 400 | 0.2 | 2.0123e+06 | 2.0123e+06 | 2.0119e+06 | 3.5053e+07 | 17.420 | 1.250 | 0.998 | 3.750e-05 | 4.979e+04 |
| ideal of table 4 (7.78e+03) | 400 | 0.1 | 2.0123e+06 | 2.0123e+06 | 2.0119e+06 | 3.1147e+07 | 15.478 | 1.111 | 0.998 | 3.750e-05 | 2.212e+04 |
| ideal of table 4 (7.78e+03) | 600 | 0.2 | 3.0181e+06 | 3.0181e+06 | 3.0179e+06 | 5.2329e+07 | 17.338 | 1.250 | 0.998 | 2.500e-05 | 7.451e+04 |
| ideal of table 4 (7.78e+03) | 600 | 0.1 | 3.0181e+06 | 3.0181e+06 | 3.0179e+06 | 4.6507e+07 | 15.409 | 1.111 | 0.998 | 2.500e-05 | 3.311e+04 |
| Ψ = (sin t, cos t) (1) | 20 | 0.2 | 8.0250e-02 | 8.0250e-02 | 8.0000e-02 | 2.1599e-01 | 2.691 | 1.250 | 0.816 | 1.000e-04 | 5.167e-03 |
| Ψ = (sin t, cos t) (1) | 20 | 0.1 | 8.0250e-02 | 8.0250e-02 | 8.0000e-02 | 2.0249e-01 | 2.523 | 1.111 | 0.791 | 1.000e-04 | 2.422e-03 |
| Ψ = (sin t, cos t) (1) | 50 | 0.2 | 2.0010e-01 | 2.0010e-01 | 2.0000e-01 | 3.1506e-01 | 1.575 | 1.250 | 0.933 | 4.000e-05 | 1.186e-02 |
| Ψ = (sin t, cos t) (1) | 50 | 0.1 | 2.0010e-01 | 2.0010e-01 | 2.0000e-01 | 2.8203e-01 | 1.409 | 1.111 | 0.928 | 4.000e-05 | 5.307e-03 |
| Ψ = (sin t, cos t) (1) | 200 | 0.2 | 8.0002e-01 | 8.0002e-01 | 8.0000e-01 | 1.0509e+00 | 1.314 | 1.250 | 0.984 | 1.000e-05 | 4.901e-02 |
| Ψ = (sin t, cos t) (1) | 200 | 0.1 | 8.0002e-01 | 8.0002e-01 | 8.0000e-01 | 9.3452e-01 | 1.168 | 1.111 | 0.984 | 1.000e-05 | 2.179e-02 |

Reading: the exact stationary spread is the averaged Langevin value k·tr X (M̄X + XM̄ = ΣΣᵀ) up to O(k⁻²) on both instances (Σ is isotropic, so ΣΣᵀ commutes with M̄; for a non-commuting ΣΣᵀ the deviation is O(1/k)). On the isotropic sinusoid (κ = 1) the explicit bound exceeds the exact value by a factor that decreases toward 1/(1 − δ) as k grows (the Itô term's hi/lo factors and βL_M² are O(1/k)); on the ill-conditioned ideal the bound pays the condition number of M̄ in the trace (it replaces every 1/λ_i by 1/λ_min), as a scalar Lyapunov function must. The bound on E‖w‖² (the error off its slow manifold) is coarse: the exact E‖e‖² is ‖Σ‖²_F/(2k) = tr(M̄X)/k (the error rides on its slow manifold e ≈ Ψζ/k) while the bound carries the factor 1/β.

Numbers are produced by `experiments/exp_explicit_rate.py`; nothing is edited by hand.
