# Corollary J′.3: the Gibbs measures π_m are chaotic with limit μ*_τ (`exp_gibbs_m.py`)

Instance of `exp_dither_stationary.py`; θ-grid of 96 points; μ*_τ on the 201² grid by 59 damped fixed-point iterations (final residual 9.4e-12); π_m by MALA with 4000 replicas, 1500 steps of dt = 0.05 (first 375 discarded), standard errors from the spread of the per-replica averages (independent replicas). The m = 1 row is cross-checked against the quadrature of `dither.GibbsGrid`.

μ*_τ: ∫α dμ* = +0.49632, ∫γ dμ* = +0.18734, ∫α² dμ* = +0.93526, ∫γ² dμ* = +0.56836, L̄(μ*_τ) = 0.13943. Quadrature of π₁: E[α] = +0.47720, E[γ] = +0.25015, E[α²] = +0.53568, E[γ²] = +0.52916, E[L̄] = 0.19640.

| m | E[α₁] ± s.e. | E[γ₁] | E[α₁²] | E[γ₁²] | E[L̄(μ_Θ)] | acceptance | (E[α₁] − ∫α dμ*)·m | (E[γ₁] − ∫γ dμ*)·m | (E[L̄(μ_Θ)] − L̄(μ*))·m |
|---|---|---|---|---|---|---|---|---|---|
| 1 | +0.47941 ± 0.00237 | +0.25117 ± 0.00345 | +0.53596 ± 0.00220 | +0.52858 ± 0.00226 | 0.19603 ± 0.00035 | 0.996 | -0.0169 | +0.0638 | +0.0566 |
| 2 | +0.49547 ± 0.00180 | +0.22468 ± 0.00245 | +0.72923 ± 0.00203 | +0.54595 ± 0.00156 | 0.16812 ± 0.00023 | 0.995 | -0.0017 | +0.0747 | +0.0574 |
| 4 | +0.49746 ± 0.00130 | +0.20401 ± 0.00172 | +0.83244 ± 0.00159 | +0.55560 ± 0.00111 | 0.15403 ± 0.00015 | 0.993 | +0.0046 | +0.0667 | +0.0584 |
| 8 | +0.49531 ± 0.00094 | +0.19363 ± 0.00122 | +0.88097 ± 0.00115 | +0.56205 ± 0.00077 | 0.14701 ± 0.00011 | 0.991 | -0.0081 | +0.0503 | +0.0606 |
| 16 | +0.49594 ± 0.00067 | +0.19051 ± 0.00087 | +0.90880 ± 0.00084 | +0.56501 ± 0.00056 | 0.14318 ± 0.00008 | 0.988 | -0.0061 | +0.0506 | +0.0599 |
| 32 | +0.49756 ± 0.00048 | +0.18996 ± 0.00061 | +0.92216 ± 0.00060 | +0.56651 ± 0.00039 | 0.14126 ± 0.00005 | 0.983 | +0.0396 | +0.0837 | +0.0585 |

Slopes of log|E[·] − (μ*_τ value)| against log m over m = 2–32: -0.05 (α), -1.01 (γ), -0.99 (L̄) (Corollary J′.3: the deviations vanish; the standard rate for mean-field Gibbs measures is 1/m, slope −1). The last three columns are the rescaled deviations; a constant means a clean 1/m law. The Gibbs average of the occupation loss is the cleanest observable: it is a smooth functional of the empirical measure and its 1/m correction is the finite-size (fluctuation) energy of the Gibbs measure.

Numbers are produced by `experiments/exp_gibbs_m.py`; nothing is edited by hand.
