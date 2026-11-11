# Theorem D‴.23: adaptive priors — a learned prior removes the leak's bias (`exp_adaptive_prior.py`)

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = 0.00707; the OT leak toward the prior with centres shifted by +0.15 (masses 1/3); the adaptive prior moves its atoms toward the class barycentres at the rate η. Midpoint rule dt = 0.25 slow units.

## A. The linear theory at the exact representation

The whitened weighted Gauss–Newton eigenvalues h_j of the loss (the collective block without the leak), the fraction σ₀/(h_j + σ₀) of the prior's error that a fixed prior keeps in mode j, and the slower rate λ_j^− of the collective pair for each η (η = ∞: λ_j^− = h_j). The separations decay at σ₀ in every case.

| j | h_j | σ₀/(h_j + σ₀) | λ_j^− at η = 0.00071 | λ_j^− at η = 0.00707 | λ_j^− at η = 0.07074 | λ_j^− at η = ∞ |
|---|---|---|---|---|---|---|
| 1 | 3.0587e-04 | 0.959 | 2.684e-05 | 1.513e-04 | 2.780e-04 | 3.059e-04 |
| 2 | 1.9644e-03 | 0.783 | 1.447e-04 | 9.144e-04 | 1.782e-03 | 1.964e-03 |
| 3 | 3.1144e-03 | 0.694 | 2.061e-04 | 1.388e-03 | 2.821e-03 | 3.114e-03 |
| 4 | 1.5149e-02 | 0.318 | 4.773e-04 | 4.284e-03 | 1.348e-02 | 1.515e-02 |
| 5 | 2.1104e-02 | 0.251 | 5.264e-04 | 4.922e-03 | 1.858e-02 | 2.110e-02 |
| 6 | 3.2479e-02 | 0.179 | 5.790e-04 | 5.600e-03 | 2.788e-02 | 3.248e-02 |
| 7 | 8.8615e-02 | 0.074 | 6.548e-04 | 6.513e-03 | 5.760e-02 | 8.862e-02 |
| 8 | 1.3420e-01 | 0.050 | 6.718e-04 | 6.702e-03 | 6.425e-02 | 1.342e-01 |
| 9 | 1.9521e-01 | 0.035 | 6.826e-04 | 6.818e-03 | 6.704e-02 | 1.952e-01 |

A fixed prior's collective rates are h_j + σ₀ = 7.3803e-03–2.0229e-01.

(v)'s radius. A sampled lower bound on the Lipschitz constant C_N of the collective Jacobian in the whitened weighted norm (the unregularised K₀-atom drift by central differences at Θ*(q) and at 160 random points at whitened distances 0.005–0.05 from the exact representation, where the Jacobian is −H_L to 7.7e-12): C_N ≥ 0.3406. The proof's radius is then at most min{σ₀, λ_1^−}/(4C_N) (the plan-margin term can only lower it), against the norm |||·|||_η of Part C's start (Θ*(q) and the shifted prior; at η = ∞ the prior is the class barycentres and the norm has no ζ part):

| η | λ_1^− | min{σ₀, λ_1^−}/σ₀ | upper bound on ρ | |||start|||_η | start / bound |
|---|---|---|---|---|---|
| 0.00071 | 2.684e-05 | 0.0038 | 1.97e-05 | 1.172e+00 | 59507 |
| 0.00707 | 1.513e-04 | 0.0214 | 1.11e-04 | 3.713e-01 | 3344 |
| 0.07074 | 2.780e-04 | 0.0393 | 2.04e-04 | 1.191e-01 | 584 |
| ∞ | 3.059e-04 | 0.0432 | 2.24e-04 | 2.102e-02 | 94 |

The runs start outside the ball that the proof of (v) covers: Part C's start by the factors of the table; Part D's warm clouds (largest separations of order 5 at s = 0) by 2.5e+05–3.6e+05 at η = 0.00071; 4.5e+04–6.4e+04 at η = 0.00707; 2.4e+04–3.5e+04 at η = 0.07074; 2.2e+04–3.2e+04 at η = ∞; Part F's separated starts by 3344 at η = 0.00707; 131–1127 at η = ∞. They probe the dynamics beyond it.

## B. The fixed prior's bias (shifted prior)

Θ*(q) by Newton (residual 4.3e-17; diagonal plan ε = 0.389 < Δ_ν/2 = 3.162); fit error at Θ*(q) 9.7765e-03 (0 at the exact representation). Mode amplitudes (components on the eigenvectors of H_L, whitened weighted coordinates) of the prior's error x⁰, of the actual bias y = Θ*(q) − exact, and of the first-order prediction σ₀(H_L + σ₀I)⁻¹x⁰; ‖y − prediction‖/‖y‖ = 4.805e-02.

| j | x⁰_j | y_j (actual) | predicted | actual/x⁰ | actual/predicted |
|---|---|---|---|---|---|
| 1 | -1.9429e-03 | -1.5815e-03 | -1.8624e-03 | 0.814 | 0.849 |
| 2 | +1.3012e-02 | +1.0715e-02 | +1.0184e-02 | 0.823 | 1.052 |
| 3 | +3.6625e-03 | +1.8500e-03 | +2.5430e-03 | 0.505 | 0.727 |
| 4 | +4.2598e-03 | +1.2441e-03 | +1.3560e-03 | 0.292 | 0.917 |
| 5 | +2.1436e-02 | +5.6414e-03 | +5.3817e-03 | 0.263 | 1.048 |
| 6 | +2.6124e-02 | +4.7319e-03 | +4.6724e-03 | 0.181 | 1.013 |
| 7 | -1.3160e-01 | -9.9217e-03 | -9.7291e-03 | 0.075 | 1.020 |
| 8 | +1.2244e-01 | +6.3699e-03 | +6.1314e-03 | 0.052 | 1.039 |
| 9 | +3.2208e-01 | +1.1232e-02 | +1.1264e-02 | 0.035 | 0.997 |

The first-order formula against the shift of the prior's centres (Newton residual of Θ*(q); the relative error in norm; actual/predicted in the modes 1–4):

| shift | Newton residual | ‖y − prediction‖/‖y‖ | actual/predicted, modes 1, 2, 3, 4 |
|---|---|---|---|
| 0.15 | 4.3e-17 | 4.805e-02 | 0.849, 1.052, 0.727, 0.917 |
| 0.075 | 3.3e-16 | 2.341e-02 | 0.917, 1.023, 0.871, 0.963 |
| 0.0375 | 1.3e-16 | 1.157e-02 | 0.956, 1.011, 0.937, 0.982 |
| 0.01875 | 6.0e-17 | 5.754e-03 | 0.977, 1.005, 0.969, 0.991 |

## C. The K₀-atom adaptive loop from Θ*(q) with the prior at the shifted centres

Mode amplitudes of the atoms' deviation from the exact representation at s = 0, 1000, 3000, 6000 against the exact solution of the linear pairs (started from the actual Θ*(q) and the actual prior); the fit error and the energy.

| η | mode | amplitude at s = 0, 1000, 3000, 6000 | linear prediction at the same times | ratio at s = 1000, 3000, 6000 | log-slope on [3000, 6000] / λ_j^− | fit error at s = 0, 1000, 3000, 6000 | 𝒥 at s = 6000 |
|---|---|---|---|---|---|---|---|
| 0.00071 | 1 | -1.582e-03, -1.621e-03, -1.570e-03, -1.450e-03 | -1.582e-03, -1.796e-03, -1.702e-03, -1.570e-03 | 0.903, 0.922, 0.923 | 0.987 | 9.776e-03, 5.046e-03, 1.402e-03, 3.612e-04 | 1.658e-07 |
| 0.00071 | 2 | +1.072e-02, +9.278e-03, +6.893e-03, +4.461e-03 | +1.072e-02, +8.974e-03, +6.719e-03, +4.352e-03 | 1.034, 1.026, 1.025 | 1.002 |  |  |
| 0.00071 | 3 | +1.850e-03, +1.708e-03, +1.192e-03, +6.459e-04 | +1.850e-03, +2.083e-03, +1.379e-03, +7.433e-04 | 0.820, 0.864, 0.869 | 0.992 |  |  |
| 0.00071 | 4 | +1.244e-03, +8.167e-04, +3.237e-04, +7.943e-05 | +1.244e-03, +8.589e-04, +3.307e-04, +7.898e-05 | 0.951, 0.979, 1.006 | 0.981 |  |  |
| 0.00071 | 5 | +5.641e-03, +3.343e-03, +1.153e-03, +2.364e-04 | +5.641e-03, +3.240e-03, +1.131e-03, +2.330e-04 | 1.032, 1.020, 1.015 | 1.003 |  |  |
| 0.00071 | 6 | +4.732e-03, +2.670e-03, +8.345e-04, +1.469e-04 | +4.732e-03, +2.658e-03, +8.347e-04, +1.469e-04 | 1.005, 1.000, 1.000 | 1.000 |  |  |
| 0.00071 | 7 | -9.922e-03, -5.155e-03, -1.384e-03, -1.940e-04 | -9.922e-03, -5.090e-03, -1.374e-03, -1.927e-04 | 1.013, 1.007, 1.007 | 1.000 |  |  |
| 0.00071 | 8 | +6.370e-03, +3.226e-03, +8.317e-04, +1.103e-04 | +6.370e-03, +3.147e-03, +8.209e-04, +1.094e-04 | 1.025, 1.013, 1.008 | 1.003 |  |  |
| 0.00071 | 9 | +1.123e-02, +5.699e-03, +1.456e-03, +1.879e-04 | +1.123e-02, +5.711e-03, +1.458e-03, +1.881e-04 | 0.998, 0.999, 0.999 | 1.000 |  |  |
| 0.00707 | 1 | -1.582e-03, -1.421e-03, -1.049e-03, -6.663e-04 | -1.582e-03, -1.500e-03, -1.108e-03, -7.040e-04 | 0.947, 0.946, 0.947 | 1.000 | 9.776e-03, 3.153e-04, 5.717e-05, 1.825e-05 | 7.096e-11 |
| 0.00707 | 2 | +1.072e-02, +4.517e-03, +7.254e-04, +4.670e-05 | +1.072e-02, +4.435e-03, +7.123e-04, +4.585e-05 | 1.019, 1.018, 1.018 | 1.000 |  |  |
| 0.00707 | 3 | +1.850e-03, +5.554e-04, +3.425e-05, +4.782e-07 | +1.850e-03, +6.277e-04, +3.911e-05, +6.082e-07 | 0.885, 0.876, — | — |  |  |
| 0.00707 | 4 | +1.244e-03, +2.490e-05, +1.651e-07, +1.986e-08 | +1.244e-03, +2.234e-05, +4.243e-09, +1.111e-14 | 1.115, —, — | — |  |  |
| 0.00707 | 5 | +5.641e-03, +4.704e-05, +3.469e-08, +5.253e-09 | +5.641e-03, +4.694e-05, +2.489e-09, +9.610e-16 | 1.002, —, — | — |  |  |
| 0.00707 | 6 | +4.732e-03, +2.019e-05, -2.310e-08, -1.649e-08 | +4.732e-03, +2.001e-05, +2.734e-10, +1.381e-17 | 1.009, —, — | — |  |  |
| 0.00707 | 7 | -9.922e-03, -1.599e-05, -2.844e-08, +2.816e-09 | -9.922e-03, -1.549e-05, -3.409e-11, -1.114e-19 | 1.032, —, — | — |  |  |
| 0.00707 | 8 | +6.370e-03, +7.979e-06, -3.637e-08, -8.988e-09 | +6.370e-03, +7.903e-06, +1.192e-11, +2.206e-20 | 1.010, —, — | — |  |  |
| 0.00707 | 9 | +1.123e-02, +1.269e-05, +2.622e-09, +3.392e-09 | +1.123e-02, +1.275e-05, +1.524e-11, +1.993e-20 | 0.996, —, — | — |  |  |
| 0.07074 | 1 | -1.582e-03, -1.196e-03, -6.858e-04, -2.979e-04 | -1.582e-03, -1.222e-03, -7.010e-04, -3.045e-04 | 0.978, 0.978, 0.978 | 1.000 | 9.776e-03, 1.308e-04, 1.884e-05, 8.035e-06 | 1.358e-11 |
| 0.07074 | 2 | +1.072e-02, +1.848e-03, +5.238e-05, +2.564e-07 | +1.072e-02, +1.836e-03, +5.203e-05, +2.483e-07 | 1.007, 1.007, — | — |  |  |
| 0.07074 | 3 | +1.850e-03, +1.128e-04, +3.384e-07, -9.957e-09 | +1.850e-03, +1.199e-04, +4.255e-07, +8.993e-11 | 0.941, —, — | — |  |  |
| 0.07074 | 4 | +1.244e-03, +6.503e-07, +2.122e-08, +3.563e-09 | +1.244e-03, +2.142e-09, +4.167e-21, +1.130e-38 | —, —, — | — |  |  |
| 0.07074 | 5 | +5.641e-03, -6.404e-08, +5.800e-09, +5.237e-10 | +5.641e-03, +6.135e-11, +4.430e-27, +2.717e-51 | —, —, — | — |  |  |
| 0.07074 | 6 | +4.732e-03, +2.316e-08, -1.743e-08, -3.391e-09 | +4.732e-03, +5.551e-15, +3.386e-39, +1.613e-75 | —, —, — | — |  |  |
| 0.07074 | 7 | -9.922e-03, -1.117e-07, +2.846e-09, +8.489e-10 | -9.922e-03, -1.993e-27, -1.844e-77, -1.642e-152 | —, —, — | — |  |  |
| 0.07074 | 8 | +6.370e-03, -5.311e-08, -9.603e-09, -1.622e-09 | +6.370e-03, +1.362e-30, +2.142e-86, +4.221e-170 | —, —, — | — |  |  |
| 0.07074 | 9 | +1.123e-02, -4.462e-09, +3.566e-09, +7.371e-10 | +1.123e-02, +1.277e-31, +7.465e-90, +3.338e-177 | —, —, — | — |  |  |
| ∞ | 1 | -1.582e-03, -1.153e-03, -6.253e-04, -2.498e-04 | -1.582e-03, -1.165e-03, -6.318e-04, -2.524e-04 | 0.990, 0.990, 0.990 | 1.000 | 9.776e-03, 1.078e-04, 1.699e-05, 6.738e-06 | 9.543e-12 |
| ∞ | 2 | +1.072e-02, +1.506e-03, +2.963e-05, +8.633e-08 | +1.072e-02, +1.503e-03, +2.955e-05, +8.150e-08 | 1.003, 1.003, — | — |  |  |
| ∞ | 3 | +1.850e-03, +7.978e-05, +1.088e-07, -7.039e-09 | +1.850e-03, +8.215e-05, +1.620e-07, +1.418e-11 | 0.971, —, — | — |  |  |
| ∞ | 4 | +1.244e-03, +4.577e-07, +1.683e-08, +2.504e-09 | +1.244e-03, +3.280e-10, +2.280e-23, +4.180e-43 | —, —, — | — |  |  |
| ∞ | 5 | +5.641e-03, -1.653e-08, +3.902e-09, +3.653e-10 | +5.641e-03, +3.855e-12, +1.800e-30, +5.745e-58 | —, —, — | — |  |  |
| ∞ | 6 | +4.732e-03, +4.631e-09, -1.467e-08, -2.385e-09 | +4.732e-03, +3.712e-17, +2.284e-45, +1.102e-87 | —, —, — | — |  |  |
| ∞ | 7 | -9.922e-03, -8.316e-08, +2.898e-09, +5.984e-10 | -9.922e-03, -3.247e-41, -3.477e-118, -1.219e-233 | —, —, — | — |  |  |
| ∞ | 8 | +6.370e-03, -4.982e-08, -7.668e-09, -1.139e-09 | +6.370e-03, +3.341e-61, +9.190e-178, +0.000e+00 | —, —, — | — |  |  |
| ∞ | 9 | +1.123e-02, -2.001e-09, +3.072e-09, +5.186e-10 | +1.123e-02, +1.868e-87, +5.166e-257, +0.000e+00 | —, —, — | — |  |  |

Ratio: the amplitude over the linear prediction (constant in time when the rate is the predicted one; its offset from 1 is set in the nonlinear first phase, the shift being finite). Log-slope: the decay rate of the amplitude on [3000, 6000] over λ_j^−. Both only where the amplitudes are above 1e-06: below it the faster modes stop at the level of the second-order remainder driven by the slow modes (the amplitudes below 1e-06 at s ≥ 3000, every mode and η: 5.0e-11–9.7e-07).

Over every mode and η above 1e-06: log-slopes 0.981–1.003 of λ_j^− (13 modes), ratios 0.820–1.115; the largest change of a ratio from s = 1000 to s = 3000 or 6000 (modes with a ratio at s = 1000 and a later one): 5.96% at η = 0.00071, 1.04% at η = 0.00707, 0.02% at η = 0.07074, 0.01% at η = ∞.

## D. m-particle systems from the warm clouds (shifted prior)

Sizes m = 15, 60, 240 and M = 510, 4 seeds, 3000 slow units, plan every 4 steps. Means over the seeds (min–max). ‖ȳ‖: the whitened weighted norm of the class means' deviation from the exact representation (classes by the plan); W²: to the exact representation's measure. Gated: the prior is frozen while mass is in transit (some particle assigned to an atom that is not its nearest), so that the class barycentres it moves to contain no particle of another cluster ('transit'), and in addition, the first time, until the fixed-prior drift is below 10⁻⁵ ('settled': the adaptation starts from the fixed-prior stationary point); t_open is the first time it adapts. In parentheses: the K₀-atom run of Part C at s (ungated) or at s − t_open (gated, mean t_open).

| η | gated | m | t_open (seeds) | fit error at s = 0, 500, 3000 | W²(μ_Θ, μ_exact) at s = 0, 500, 3000 | ‖ȳ‖ at s = 500, 3000 | prior's distance at s = 3000 | final counts |
|---|---|---|---|---|---|---|---|---|
| 0.00000 | no | 15 | — | 2.751e-01, 9.751e-03, 9.776e-03 | 5.236e+00, 6.380e-03, 4.417e-04 | 2.120e-02 (1.95e-02–2.22e-02), 2.102e-02 (2.10e-02–2.10e-02) (Θ*(q): 2.102e-02) | 3.707e-01 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| 0.00000 | no | 60 | — | 1.344e-01, 9.774e-03, 9.776e-03 | 2.754e+00, 3.129e-03, 4.417e-04 | 2.169e-02 (1.89e-02–2.56e-02), 2.102e-02 (2.10e-02–2.10e-02) (Θ*(q): 2.102e-02) | 3.707e-01 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| 0.00000 | no | 240 | — | 6.427e-02, 9.798e-03, 9.776e-03 | 1.166e+00, 2.913e-03, 4.417e-04 | 2.466e-02 (2.06e-02–3.05e-02), 2.102e-02 (2.10e-02–2.10e-02) (Θ*(q): 2.102e-02) | 3.707e-01 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| 0.00000 | no | 510 | — | 5.251e-02, 9.786e-03, 9.776e-03 | 9.937e-01, 2.164e-03, 4.417e-04 | 2.307e-02 (2.07e-02–2.67e-02), 2.102e-02 (2.10e-02–2.10e-02) (Θ*(q): 2.102e-02) | 3.707e-01 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| 0.00707 | no | 15 | — | 2.751e-01, 1.848e-02, 2.969e-03 | 5.236e+00, 3.024e-01, 1.635e-02 | 3.229e-01 (9.05e-02–8.40e-01), 1.029e-01 (2.99e-02–2.27e-01) (7.439e-03, 1.276e-03) | 1.056e-01 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| 0.00707 | no | 60 | — | 1.344e-01, 4.725e-03, 9.864e-04 | 2.754e+00, 1.134e-02, 1.829e-03 | 7.302e-02 (5.60e-03–1.13e-01), 3.573e-02 (9.00e-04–6.54e-02) (7.439e-03, 1.276e-03) | 3.655e-02 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| 0.00707 | no | 240 | — | 6.427e-02, 4.769e-03, 3.637e-04 | 1.166e+00, 1.063e-02, 2.710e-04 | 5.215e-02 (8.52e-03–1.27e-01), 1.283e-02 (2.34e-03–2.77e-02) (7.439e-03, 1.276e-03) | 1.315e-02 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| 0.00707 | no | 510 | — | 5.251e-02, 2.975e-03, 2.627e-04 | 9.937e-01, 4.357e-03, 1.402e-04 | 3.182e-02 (7.55e-03–7.42e-02), 9.380e-03 (2.78e-03–2.10e-02) (7.439e-03, 1.276e-03) | 9.604e-03 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| ∞ | no | 15 | — | 2.751e-01, 2.868e-02, 1.170e-02 | 5.236e+00, 1.336e+00, 1.170e+00 | 7.130e-01 (1.60e-01–1.92e+00), 5.932e-01 (5.48e-02–2.12e+00) (4.259e-03, 6.260e-04) | 5.932e-01 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| ∞ | no | 60 | — | 1.344e-01, 5.917e-03, 1.256e-03 | 2.754e+00, 3.566e-02, 3.105e-03 | 1.309e-01 (3.97e-03–2.11e-01), 4.654e-02 (4.54e-04–8.34e-02) (4.259e-03, 6.260e-04) | 4.654e-02 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| ∞ | no | 240 | — | 6.427e-02, 4.649e-03, 4.742e-04 | 1.166e+00, 2.432e-02, 5.057e-04 | 8.481e-02 (1.71e-02–2.24e-01), 1.757e-02 (3.39e-03–3.83e-02) (4.259e-03, 6.260e-04) | 1.757e-02 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| ∞ | no | 510 | — | 5.251e-02, 2.816e-03, 3.332e-04 | 9.937e-01, 7.369e-03, 2.457e-04 | 4.771e-02 (1.25e-02–1.13e-01), 1.235e-02 (2.76e-03–2.77e-02) (4.259e-03, 6.260e-04) | 1.235e-02 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| 0.00707 | transit | 15 | 112, 91, 88, 149 | 2.751e-01, 6.888e-03, 1.013e-03 | 5.236e+00, 2.042e-02, 1.626e-03 | 8.101e-02 (4.76e-02–1.20e-01), 3.698e-02 (1.47e-02–5.36e-02) (8.238e-03, 1.329e-03) | 3.782e-02 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| 0.00707 | transit | 60 | 94, 93, 114, 112 | 1.344e-01, 3.241e-03, 4.846e-04 | 2.754e+00, 5.308e-03, 4.262e-04 | 4.126e-02 (1.28e-02–5.34e-02), 1.717e-02 (1.23e-03–3.30e-02) (8.238e-03, 1.329e-03) | 1.759e-02 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| 0.00707 | transit | 240 | 117, 102, 96, 102 | 6.427e-02, 4.374e-03, 2.126e-04 | 1.166e+00, 7.691e-03, 7.661e-05 | 4.273e-02 (5.45e-03–9.86e-02), 6.740e-03 (1.53e-03–1.52e-02) (8.238e-03, 1.329e-03) | 6.959e-03 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| 0.00707 | transit | 510 | 109, 99, 92, 98 | 5.251e-02, 2.951e-03, 1.547e-04 | 9.937e-01, 3.779e-03, 3.793e-05 | 2.793e-02 (7.41e-03–6.20e-02), 4.929e-03 (1.55e-03–1.08e-02) (8.238e-03, 1.329e-03) | 5.082e-03 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| ∞ | transit | 15 | 112, 91, 88, 149 | 2.751e-01, 1.207e-02, 1.659e-03 | 5.236e+00, 1.077e-01, 4.690e-03 | 2.107e-01 (8.56e-02–3.98e-01), 6.151e-02 (3.04e-02–1.07e-01) (5.116e-03, 6.458e-04) | 6.151e-02 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| ∞ | transit | 60 | 94, 93, 114, 112 | 1.344e-01, 3.564e-03, 6.098e-04 | 2.754e+00, 1.036e-02, 7.567e-04 | 6.901e-02 (1.14e-02–1.01e-01), 2.260e-02 (4.60e-04–4.40e-02) (5.116e-03, 6.458e-04) | 2.260e-02 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| ∞ | transit | 240 | 117, 102, 96, 102 | 6.427e-02, 3.919e-03, 2.572e-04 | 1.166e+00, 1.547e-02, 1.528e-04 | 6.450e-02 (8.01e-03–1.68e-01), 9.523e-03 (2.39e-03–2.17e-02) (5.116e-03, 6.458e-04) | 9.523e-03 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| ∞ | transit | 510 | 109, 99, 92, 98 | 5.251e-02, 2.461e-03, 1.772e-04 | 9.937e-01, 5.501e-03, 6.875e-05 | 3.727e-02 (6.70e-03–8.85e-02), 6.565e-03 (1.91e-03–1.46e-02) (5.116e-03, 6.458e-04) | 6.565e-03 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| 0.00707 | settled | 15 | 1276, 1235, 1258, 1211 | 2.751e-01, 9.751e-03, 1.594e-04 | 5.236e+00, 6.380e-03, 6.748e-06 | 2.120e-02 (1.95e-02–2.22e-02), 2.597e-03 (2.53e-03–2.66e-03) (2.102e-02, 2.612e-03) | 2.911e-03 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| 0.00707 | settled | 60 | 1232, 1269, 1281, 1281 | 1.344e-01, 9.774e-03, 1.624e-04 | 2.754e+00, 3.129e-03, 6.977e-06 | 2.169e-02 (1.89e-02–2.56e-02), 2.641e-03 (2.58e-03–2.67e-03) (2.102e-02, 2.660e-03) | 2.961e-03 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| 0.00707 | settled | 240 | 1299, 1286, 1241, 1304 | 6.427e-02, 9.798e-03, 1.649e-04 | 1.166e+00, 2.913e-03, 7.161e-06 | 2.466e-02 (2.06e-02–3.05e-02), 2.676e-03 (2.59e-03–2.72e-03) (2.102e-02, 2.660e-03) | 3.001e-03 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| 0.00707 | settled | 510 | 1292, 1291, 1282, 1236 | 5.251e-02, 9.786e-03, 1.638e-04 | 9.937e-01, 2.164e-03, 7.085e-06 | 2.307e-02 (2.07e-02–2.67e-02), 2.661e-03 (2.59e-03–2.69e-03) (2.102e-02, 2.660e-03) | 2.984e-03 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |
| ∞ | settled | 15 | 1276, 1235, 1258, 1211 | 2.751e-01, 9.751e-03, 3.369e-05 | 5.236e+00, 6.380e-03, 9.262e-07 | 2.120e-02 (1.95e-02–2.22e-02), 9.621e-04 (9.29e-04–9.88e-04) (2.102e-02, 9.794e-04) | 9.621e-04 | (5, 5, 5); (5, 5, 5); (5, 5, 5); (5, 5, 5) |
| ∞ | settled | 60 | 1232, 1269, 1281, 1281 | 1.344e-01, 9.774e-03, 3.471e-05 | 2.754e+00, 3.129e-03, 9.665e-07 | 2.169e-02 (1.89e-02–2.56e-02), 9.831e-04 (9.71e-04–9.92e-04) (2.102e-02, 9.922e-04) | 9.831e-04 | (20, 20, 20); (20, 20, 20); (20, 20, 20); (20, 20, 20) |
| ∞ | settled | 240 | 1299, 1286, 1241, 1304 | 6.427e-02, 9.798e-03, 3.544e-05 | 1.166e+00, 2.913e-03, 9.889e-07 | 2.466e-02 (2.06e-02–3.05e-02), 9.944e-04 (9.73e-04–1.01e-03) (2.102e-02, 9.922e-04) | 9.944e-04 | (80, 80, 80); (80, 80, 80); (80, 80, 80); (80, 80, 80) |
| ∞ | settled | 510 | 1292, 1291, 1282, 1236 | 5.251e-02, 9.786e-03, 3.516e-05 | 9.937e-01, 2.164e-03, 9.833e-07 | 2.307e-02 (2.07e-02–2.67e-02), 9.916e-04 (9.73e-04–9.99e-04) (2.102e-02, 9.922e-04) | 9.916e-04 | (170, 170, 170); (170, 170, 170); (170, 170, 170); (170, 170, 170) |

The settled gate against the K₀-atom run of Part C: per seed, ‖ȳ‖ and the fit error at s = 3000 over the K₀-atom run's at s = 3000 − t_open (min–max over the seeds): η = 0.00707, m = 15: ‖ȳ‖ 0.995–0.999, fit error 0.998–1.000; η = 0.00707, m = 60: ‖ȳ‖ 0.999–1.000, fit error 1.000–1.000; η = 0.00707, m = 240: ‖ȳ‖ 1.000–1.000, fit error 1.000–1.000; η = 0.00707, m = 510: ‖ȳ‖ 1.000–1.000, fit error 1.000–1.000; η = ∞, m = 15: ‖ȳ‖ 0.968–0.995, fit error 0.979–0.997; η = ∞, m = 60: ‖ȳ‖ 0.991–1.000, fit error 0.995–1.000; η = ∞, m = 240: ‖ȳ‖ 0.996–0.999, fit error 0.998–1.000; η = ∞, m = 510: ‖ȳ‖ 0.997–1.001, fit error 0.999–1.000.

The settled systems open replicated: at the last record before the opening (32 runs) the rms separation is 1.2e-04–7.9e-04 and the class means lie 7.4e-07–6.7e-05 from Θ*(q) (‖y*(q)‖ = 2.102e-02). Each system is then the K₀-atom system started from its own class means, so the ratios above measure the invariance of the replicated configurations and the opening offset in the weakest modes, not the uniformity of (v); Part F puts the separations in at the start.

System by system at s = 3000 the transit gate ends lower than the settled one in all three of ‖ȳ‖, W² and the fit error in 6 of the 32 systems: the comparison of the table's seed means does not hold for every system.

## E. The ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000, continued

Seed 3, m = 15 (‖ȳ‖ = 2.124 at s = 3000), re-run from its seed (fit error at s = 3000: 4.009106e-02, stored 4.009106e-02) and continued to s = 23000.

| s | fit error | 𝒥 | largest spread of the class centres | largest drift |
|---|---|---|---|---|
| 0 | 5.1580e-01 | 7.826262e-02 | 3.0010 | 9.77e-02 |
| 1000 | 3.7532e-02 | 1.606886e-03 | 0.9566 | 9.78e-04 |
| 2000 | 3.4202e-02 | 1.141301e-03 | 0.7493 | 4.76e-04 |
| 3000 | 4.0091e-02 | 8.551092e-04 | 0.5627 | 3.67e-04 |
| 7000 | 3.4850e-02 | 2.552910e-04 | 0.0000 | 2.37e-05 |
| 11000 | 3.4772e-02 | 2.541452e-04 | 0.0000 | 1.25e-05 |
| 15000 | 3.4719e-02 | 2.533778e-04 | 0.0000 | 1.28e-05 |
| 19000 | 3.4676e-02 | 2.527426e-04 | 0.0000 | 1.34e-05 |
| 23000 | 3.4639e-02 | 2.522042e-04 | 0.0000 | 1.38e-05 |

At s = 23000 the classes are collapsed onto the atoms a = (+3.133, -2.808, +3.037), c = (-1.415, +0.188, +1.224), s = (0.583, 0.694, 0.692) (exact representation: a = (+3.000, -2.100, +2.400), c = (-1.500, +0.000, +1.500), s = (0.500, 0.400, 0.500)); the unregularised weighted loss there is 2.5220e-04, the eigenvalues of its Hessian (whitened, weighted) are -1.84e-05, 7.16e-05, 9.46e-04 (smallest three) — a plateau, not a minimum — and a BFGS minimisation of the unregularised loss from those atoms ends at 1.7e-17 with a = (+3.000, -2.100, +2.400), c = (-1.500, -0.000, +1.500), s = (0.500, 0.400, 0.500).

The collapsed η = ∞ dynamics from the same atoms, the unregularised K₀-atom gradient flow (LSODA, rtol 10⁻¹⁰, flow time counted from the s = 23000 state, the crossings located as events): the loss stays above 0.9 of its 2.5220e-04 until 4.1565e+04, halves at 4.7797e+04, falls below 10⁻¹² at 7.9714e+04 and ends at 7.8e-31 with a = (+3.000, -2.100, +2.400), c = (-1.500, +0.000, +1.500), s = (0.500, 0.400, 0.500): the same basin, crossed on the time scale of the negative curvature (1/1.84e-05 = 5.43e+04).

## F. Separated starts: the structure behind (v) beyond its ball

From the lift of Θ*(q) plus zero-mean per-particle separations (rms in the G⁻¹ norm; 5 draws of their direction, the same for both rms), the prior at the shifted centres, ungated, 1750 slow units: ‖ȳ‖ at the end over the K₀-atom run's at the same time (Part C) − 1, the mean over the draws and their range. A zero-separation system is the K₀-atom system by the invariance of the replicated configurations, so the departure is the separations' second-order effect, which Theorem D‴.18's facts (1)–(2) bound uniformly in m; the starts lie outside (v)'s ball (Part A), so (v) itself asserts nothing here.

| η | rms 0.02, m = 15 | rms 0.02, m = 510 | rms 0.1, m = 15 | rms 0.1, m = 510 |
|---|---|---|---|---|
| 0.00707 | +0.0013 (+0.0007 to +0.0018) | +0.0009 (+0.0007 to +0.0011) | +0.0307 (+0.0175 to +0.0438) | +0.0224 (+0.0177 to +0.0268) |
| ∞ | +0.0077 (+0.0038 to +0.0107) | +0.0063 (+0.0058 to +0.0071) | +0.1852 (+0.0953 to +0.2556) | +0.1552 (+0.1436 to +0.1750) |

From rms 0.02 to 0.1 the departure grows, draw by draw, by the factors 23.77–25.03 (the square of the ratio of the rms, 25, for a second-order effect). At second order it is linear in the separations' class covariances, whose expectation after the normalisation is the same at every m (classes of equal size, isotropic whitened directions): the prediction is a departure of the same expected size at both m, with a spread that shrinks as the covariances concentrate. Across the draws the standard deviation at m = 15 is 3.3–6.0 times that at m = 510, and the means at the two sizes differ by +0.9 to +1.7 combined standard errors (5 draws; the four columns share their draws, so these are not independent tests). The starts' norms |||·|||_η: 0.029–0.371.

## G. The limits of (i)

(a) A non-replicated stationary point. The four-bump target (A, C, S) = (1, -1.5, 0.5), (-0.7, 0, 0.4), (0.4, 1.3, 0.4), (0.4, 2.3, 0.4), which six kernels of weight 1/6 fit exactly with the third class split over the last two bumps (fit error 1.7e-16, G⁻¹-distance of the split pair 3.2892; K₀ = 3, q = 1/3, m = 6). The stationary point of the η = ∞ form (the leak toward the class barycentres) by least squares from that fit (residual 8.7e-17): the split pair at distance 1.6401, the leak force σ₀‖ϑ_i − T̂_π(i)‖ on its particles 0.0058, 0.0058 (largest on the others 3.1e-18), the unregularised gradient there 5.80e-03 (at the exact fit 5.4e-17), fit error 4.246e-02. With the prior at the class barycentres it is stationary for the joint flow — η = 0.00707: the plan the classes True, drift 4.4e-17 (Θ) and 0.0e+00 (Z), largest real part of the joint Jacobian +3.66e-05; η = 0.07074: the plan the classes True, drift 4.4e-17 (Θ) and 0.0e+00 (Z), largest real part of the joint Jacobian +6.64e-05 — a saddle: the leak holds the split class against the loss.

The joint flow from the exact split fit (40000 slow units): the split distance and the fit error at the records; at the end the classes' largest spread, the merged atom (the widest) — its width coordinate κ and width s — and the largest whitened unregularised drift at the class atoms:

| η | split distance at s = 0, 4000, 8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000 | fit error at the same s | at the end: spread, κ, s, drift |
|---|---|---|---|
| ∞ | 3.2892, 0.0076, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 | 1.730e-16, 4.362e-02, 2.964e-02, 2.621e-02, 2.476e-02, 2.394e-02, 2.340e-02, 2.301e-02, 2.271e-02, 2.247e-02, 2.227e-02 | 4.4e-16, -0.245, 1.046, 2.8e-05 |
| 0.00707 | 3.2892, 0.9801, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 | 1.730e-16, 6.012e-02, 3.928e-02, 3.218e-02, 2.891e-02, 2.709e-02, 2.595e-02, 2.518e-02, 2.462e-02, 2.419e-02, 2.385e-02 | 8.9e-16, -0.610, 0.898, 5.4e-05 |

The collapsed η = ∞ dynamics continued from the end state (the unregularised K₀-atom gradient flow, LSODA, rtol 10⁻¹⁰; s counted from the end of the run): the merged atom's width coordinate κ, the gap s_max − s of its width, the fit error.

| η | s | κ | s_max − s | fit error |
|---|---|---|---|---|
| ∞ | 1e+06 | +2.159 | 1.76e-01 | 1.955537e-02 |
| ∞ | 1e+07 | +3.843 | 3.57e-02 | 1.944406e-02 |
| ∞ | 1e+08 | +5.529 | 6.72e-03 | 1.943437e-02 |
| ∞ | 1e+09 | +7.528 | 9.14e-04 | 1.943302e-02 |
| ∞ | 1e+10 | +9.751 | 9.90e-05 | 1.943285e-02 |
| 0.00707 | 1e+06 | +2.144 | 1.78e-01 | 1.955824e-02 |
| 0.00707 | 1e+07 | +3.842 | 3.57e-02 | 1.944408e-02 |
| 0.00707 | 1e+08 | +5.529 | 6.72e-03 | 1.943437e-02 |
| 0.00707 | 1e+09 | +7.528 | 9.14e-04 | 1.943302e-02 |
| 0.00707 | 1e+10 | +9.751 | 9.90e-05 | 1.943285e-02 |

The increments of κ per decade over ln 10: η = ∞: 0.73, 0.73, 0.87, 0.97; η = 0.00707: 0.74, 0.73, 0.87, 0.97 — κ grows like ln s (the logistic saturation of the width) and s tends to s_max like 1/s: the solution leaves every compact set, as in (b). Its physical image approaches the minimiser of the loss on the face where the merged atom's width is s_max (BFGS over the other parameters: fit error 1.943283e-02, gradient 6.5e-10, inside the other bounds: True), where the loss still decreases outward (∂L/∂s = -1.326e-06): a Karush–Kuhn–Tucker point of the box-constrained loss, not a critical point, reached only as s → ∞.

The exact six-kernel fit is not stationary: the joint flow collapses the split class onto a replicated configuration with a larger fit error and then runs off — the K₀-cluster bias of (iv)'s regulariser, and a second instance of the lost a-priori bound, at a target the network fits exactly.

(b) Unbounded solutions. The target (A, C, S) = (1.6, -1.5, 0.5), (-0.7, 0, 0.4), (0.8, 1.5, 0.5): the first bump's weighted amplitude 3·1.6 = 4.8 exceeds a_max = 4, so no exact representation with q = 1/3 exists. The K₀-atom system from the instance's exact atoms with the prior there (LSODA, rtol 10⁻¹⁰): the first atom's amplitude coordinate ϑ_a and its prior atom's z_a (unconstrained), its amplitude a, the fit error and the largest drift.

| η | s | ϑ_a | z_a | a | fit error | largest drift |
|---|---|---|---|---|---|---|
| 0.00000 | 0e+00 | 3.8918 | 3.8918 | 3.000000 | 3.2951e-01 | 6.18e-02 |
| 0.00000 | 1e+02 | 5.6514 | 3.8918 | 3.552385 | 1.8564e-01 | 4.12e-03 |
| 0.00000 | 1e+03 | 5.8692 | 3.8918 | 3.596246 | 1.7498e-01 | 1.84e-08 |
| 0.00000 | 1e+04 | 5.8692 | 3.8918 | 3.596246 | 1.7498e-01 | 7.76e-17 |
| 0.00000 | 1e+05 | 5.8692 | 3.8918 | 3.596246 | 1.7498e-01 | 4.17e-17 |
| 0.00000 | 1e+06 | 5.8692 | 3.8918 | 3.596246 | 1.7498e-01 | 1.73e-17 |
| 0.00707 | 0e+00 | 3.8918 | 3.8918 | 3.000000 | 3.2951e-01 | 6.18e-02 |
| 0.00707 | 1e+02 | 5.7754 | 4.5893 | 3.577871 | 1.7971e-01 | 8.39e-03 |
| 0.00707 | 1e+03 | 7.9940 | 7.7684 | 3.855686 | 1.2530e-01 | 1.60e-03 |
| 0.00707 | 1e+04 | 11.7934 | 11.7668 | 3.978072 | 1.0792e-01 | 1.88e-04 |
| 0.00707 | 1e+05 | 16.2457 | 16.2429 | 3.997627 | 1.0527e-01 | 1.98e-05 |
| 0.00707 | 1e+06 | 20.8276 | 20.8273 | 3.999760 | 1.0498e-01 | 2.00e-06 |
| ∞ | 0e+00 | 3.8918 | — | 3.000000 | 3.2951e-01 | 6.18e-02 |
| ∞ | 1e+02 | 6.1997 | — | 3.655085 | 1.6257e-01 | 1.06e-02 |
| ∞ | 1e+03 | 9.1282 | — | 3.917504 | 1.1623e-01 | 1.51e-03 |
| ∞ | 1e+04 | 13.1310 | — | 3.988751 | 1.0647e-01 | 1.90e-04 |
| ∞ | 1e+05 | 17.6245 | — | 3.998809 | 1.0511e-01 | 1.98e-05 |
| ∞ | 1e+06 | 22.2129 | — | 3.999880 | 1.0497e-01 | 2.00e-06 |

The increments of ϑ_a per decade of s over ln 10, from s = 10³ on: η = 0.00707: 1.650, 1.934, 1.990; η = ∞: 1.738, 1.951, 1.993 — ϑ_a grows like 2 ln s, and a approaches a_max like 1/s: under an adaptive prior the solution leaves every compact set, while the fixed prior's converges.

(c) The relabelling. The exact representation (m = 6) with the prior's atoms permuted cyclically (largest offset 6.225, unconstrained coordinates): the plan maps the classes to the atoms (2, 2, 0, 0, 1, 1), W² = 0.0e+00, drift 1.1e-17 (Θ) and 3.1e-18 (Z) — stationary, with z_σ(k) = ϑ_k for the plan's bijection σ.

## H. Proposition D‴.27: the projected adaptive prior

The box 𝒵 for the prior's atoms (unconstrained coordinates): α ∈ [-5.8889, 5.8889], γ ∈ [-4.4167, 4.4167], κ ∈ [-6, 6] (amplitudes and centres to 0.9 of a_max and c_max). (a) Part G(b)'s target (the first bump's weighted amplitude 4.8 > a_max), the K₀-atom system from the instance's exact atoms with the prior there, the prior projected onto 𝒵 (η = σ₀: the projected midpoint rule — positions projected inside the step — at dt = 0.25 (the rows) and 0.125, to s = 1e+05; η = ∞: LSODA, rtol 10⁻¹⁰, the prior the projected class mean): the first atom's amplitude coordinate ϑ_a, its prior atom's z_a, the amplitude a, the fit error and the largest drift (of the projected flow).

| η | s | ϑ_a | z_a | a | fit error | largest drift |
|---|---|---|---|---|---|---|
| 0.00707 | 0e+00 | 3.8918 | 3.8918 | 3.000000 | 3.2951e-01 | 6.18e-02 |
| 0.00707 | 1e+02 | 5.7753 | 4.5893 | 3.577870 | 1.7971e-01 | 8.39e-03 |
| 0.00707 | 1e+03 | 6.7842 | 5.8889 | 3.739660 | 1.4197e-01 | 4.28e-04 |
| 0.00707 | 1e+04 | 6.7290 | 5.8889 | 3.732612 | 1.4152e-01 | 7.86e-06 |
| 0.00707 | 1e+05 | 6.7278 | 5.8889 | 3.732452 | 1.4155e-01 | 1.75e-15 |
| ∞ | 0e+00 | 3.8918 | 3.8918 | 3.000000 | 3.2951e-01 | 6.18e-02 |
| ∞ | 1e+02 | 6.1736 | 5.8889 | 3.650755 | 1.6330e-01 | 8.82e-03 |
| ∞ | 1e+03 | 6.7547 | 5.8889 | 3.735912 | 1.4141e-01 | 3.37e-04 |
| ∞ | 1e+04 | 6.7278 | 5.8889 | 3.732457 | 1.4154e-01 | 4.03e-07 |
| ∞ | 1e+05 | 6.7278 | 5.8889 | 3.732452 | 1.4155e-01 | 1.18e-16 |
| ∞ | 1e+06 | 6.7278 | 5.8889 | 3.732452 | 1.4155e-01 | 2.50e-17 |

At dt = 0.125 the η = σ₀ run ends within 3.4e-12 of the dt = 0.25 run (atoms and prior), with the projected flow's drift 3.5e-15 at s = 1e+05 (1.8e-15 at dt = 0.25) and the prior's amplitude atom on the face (gaps 0.0e+00, 0.0e+00): the shared fixed point. The scheme's consistency: at s = 300 the trajectory (atoms and prior) differs from a reference at dt = 0.015625 by 6.3e-07 (dt = 0.25), 1.9e-07 (dt = 0.125), 3.8e-08 (dt = 0.0625).

The largest G⁻¹-distance of a particle from 𝒵 along the run: η = 0.00707: 0.9165; η = ∞: 0.8880 (η = σ₀ at dt = 0.125: 0.9165), against the a-priori bound B/σ₀ ≤ 3686.3 (B ≤ 26.08 = (max|f| + a_max)·L, L = 4.658 the closed-form bound of L_ϑ^G; 4.044 by search, `learning_table.md`). Both priors converge to the same bounded stationary point, with the prior's amplitude atom on the boundary of 𝒵; in this run the fit error lies between the fixed prior's and the unprojected adaptive prior's limit (Part G(b)) — (c) shows its dependence on the box.

(b) Part C's K₀-atom loop at the misplaced prior (η = σ₀, 6000 slow units) with and without 𝒵: the largest difference of the atoms and of the prior at the end is 0.0e+00 (the run stays inside 𝒵: True). The box is never approached — along the run the prior keeps a margin of at least 16.9% of the box's width to every face, the particles 16.9% — so the two runs coincide by construction (the projection is the identity inside 𝒵).

(c) The box's size: (a)'s η = ∞ system (LSODA to s = 1e+06) with the amplitude bound of 𝒵 varied (the centre and width bounds as above).

| amplitude bound (× a_max) | the α-face | the start (the exact atoms) inside 𝒵 | ϑ_a | z_a | fit error | largest drift |
|---|---|---|---|---|---|---|
| 0.5 | 2.1972 | no | 5.0616 | 2.1972 | 1.9490e-01 | 5.6e-17 |
| 0.6 | 2.7726 | no | 5.2358 | 2.7726 | 1.8451e-01 | 6.8e-17 |
| 0.65 | 3.1012 | no | 5.3436 | 3.1012 | 1.7921e-01 | 2.3e-17 |
| 0.7 | 3.4692 | no | 5.4730 | 3.4692 | 1.7362e-01 | 5.7e-17 |
| 0.75 | 3.8918 | yes | 5.6352 | 3.8918 | 1.6763e-01 | 1.7e-17 |
| 0.8 | 4.3944 | yes | 5.8511 | 4.3944 | 1.6108e-01 | 5.4e-17 |
| 0.9 | 5.8889 | yes | 6.7278 | 5.8889 | 1.4155e-01 | 2.5e-17 |
| 0.95 | 7.3271 | yes | 7.7789 | 7.3271 | 1.2679e-01 | 9.4e-17 |
| 0.99 | 10.5866 | yes | 10.6806 | 10.5866 | 1.1012e-01 | 2.0e-17 |

The fixed prior's limit is 1.7498e-01 and the unprojected adaptive prior's fit error at s = 1e+06 1.0497e-01 (Part G(b)): for the amplitude bounds 0.7–0.99 the box's fit error lies between them, decreasing as the box grows; the tighter boxes are worse than the fixed prior (0.1949 at 0.5, 0.1845 at 0.6, 0.1792 at 0.65), and the boxes with the bound ≤ 0.7 exclude the instance's exact atoms (the start is projected at once): a box should contain the representation.

(e) At η = ∞ the plan is Lloyd-consistent — optimal for its own class means — not jointly optimal: a configuration of six particles (zero amplitudes, the target 0.02 × the three bumps, the box; found by a seeded search, trial 1) with the prior at the class means of the consistent pairing [(0, 2), (1, 4), (3, 5)]. Along AdaptivePriorLoop's run (plan 'ot', η = ∞, 400 slow units) the selected energy 𝒥(Θ, Z) falls (4.5195e-03 → 4.9483e-05, its largest increment -3.4e-06), while the reduced energy min_Z𝒥 (the minimum over the 15 balanced pairings) rises from 1.2512e-03 to 1.5742e-03 (its largest increment +6.00e-05, at s = 60–70) before the plan becomes optimal (at s = 80): the selected energy is the Lyapunov function, the reduced energy is not.

(d) Part G(a)'s four-bump target with 𝒵 at η = ∞: the joint flow from the exact split fit (40000 slow units; plan classes (0, 0, 1, 1, 2, 2), largest spread 4.4e-16, the merged atom's κ = -0.245, fit error 2.226877e-02), then the collapsed dynamics with 𝒵 (the leak toward the projection of the atoms; LSODA, rtol 10⁻¹⁰; s counted from the end of the run): the merged atom's width coordinate κ, its prior atom's, the fit error and the largest whitened drift.

| s | κ | the prior's κ | fit error | largest drift |
|---|---|---|---|---|
| 1e+06 | +2.159255 | +2.1593 | 1.955537e-02 | 1.6e-06 |
| 1e+07 | +3.843184 | +3.8432 | 1.944406e-02 | 1.4e-07 |
| 1e+08 | +5.529437 | +5.5294 | 1.943437e-02 | 1.6e-08 |
| 1e+09 | +6.000001 | +6.0000 | 1.943376e-02 | 5.3e-17 |
| 1e+10 | +6.000001 | +6.0000 | 1.943376e-02 | 1.8e-16 |

Without 𝒵 (Part G(a)) κ grows like ln s and the solution leaves every compact set; with 𝒵 the merged atom's prior stops on the face κ = 6 and the atom stays within 6.6e-07 of it: a bounded stationary point, its fit error within 9.3e-07 of the face minimiser's of Part G(a).

Numbers are produced by `experiments/exp_adaptive_prior.py`; nothing is edited by hand.
