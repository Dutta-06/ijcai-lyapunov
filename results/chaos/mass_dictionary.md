# Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (`exp_mass_dictionary.py`)

The loop at fixed positions: masses β (Σβ = m = 12, or free), amplitudes a = a_max tanh(α/a_max) (a_max = 1), w = β⊙a/m, G = 1, γ = 0.1; label noise on the measured error; the selection along the exact fits (`mass_dictionary`), in s = (m v/4)t. The paradigm's update: Euler on (α, log β), renormalised; the flat one: Euler on (α, √β), renormalised radially.

## A. Identities

At random exact fits (5 × 12, three points per law): the closed-form drift against Katzenberger's drift of the scheme (`gauge_drift.limit_drift` on a chart of M), largest relative difference:

| masses | the draws | the flat update | the paradigm's update (with its Itô term) |
|---|---|---|---|
| conserved | fresh, h → 0 | 2.4e-06 | 2.4e-06 |
| conserved | AR(1), τ_c = 200, h = 0.5 | 1.4e-06 | 1.4e-06 |
| free | fresh, h → 0 | 1.6e-06 | 1.6e-06 |
| free | AR(1), τ_c = 200, h = 0.5 | 2.0e-06 | 2.0e-06 |

The paradigm's drift against the flat drift without the masses' own selection (the γ-term of D_i differentiated in √β), largest relative difference (the conserved law's leftover is the centring's cross term, `cross_term`; its size along the flows below): free, fresh 1.7e-15; free, AR(1) 1.0e-15; conserved, fresh 2.9e-03; conserved, AR(1) 8.1e-03.

Free masses, the gauge invariant Q_i = β_i − c cosh(2α_i/a_max) (Proposition D‴.45): ∇Q·(Katzenberger's drift) against the closed form (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)) (κ = 1 flat, 2 paradigm), largest relative difference: fresh, flat 8.6e-07; fresh, paradigm 8.6e-07; AR(1), flat 1.5e-06; AR(1), paradigm 1.6e-06; Q's quadratic variation under the paradigm's update relative to its scale: fresh 8.5e-17, AR(1) 6.2e-17.

Loop integrals ∮⟨P⁻¹q, dθ⟩ of the white drifts on M around loops of radius ε (in the metric) at one point per law — zero for a gradient flow in the flow's metric, ε² times a 2-form otherwise:

| masses | the update | ε = 0.01 | ε = 0.005 | ε = 0.0025 | ratios |
|---|---|---|---|---|---|
| conserved | flat | 3.59e-20 | -2.39e-20 | -4.57e-21 | -1.50, 5.24 |
| conserved | the paradigm | 2.40e-09 | 6.00e-10 | 1.50e-10 | 4.00, 4.00 |
| free | flat | -1.03e-19 | 1.36e-20 | 6.65e-22 | -7.54, 20.50 |
| free | the paradigm | 6.25e-09 | 1.56e-09 | 3.91e-10 | 4.00, 4.00 |

K's spectrum against the flow's normal rates (−eigenvalues of its Jacobian in (α, √β)): largest relative difference conserved 3.7e-11, free 2.6e-11.

The excess-mass bound at the conserved law's points (tr K the conserved law's, K_free + its −(γ/N)yyᵀ):

| E/m | Σβa′/m | tr K | tr K_free | (G/(Nm))E²/Σ|φ_i|⁻² |
|---|---|---|---|---|
| 0.4964 | 0.6847 | 1.1594e-01 | 1.1702e-01 | 4.9288e-02 |
| 0.5277 | 0.7008 | 1.2405e-01 | 1.2513e-01 | 5.5699e-02 |
| 0.4966 | 0.6780 | 1.2363e-01 | 1.2472e-01 | 4.9316e-02 |

The gauge-minimal relevance D*(|w|) against (γa_max/m)|w| and the 4/3 law (the threshold γa_max³/(Gm) = 0.008333); its minimiser's t = |a|/a_max (the root of 2Gm|w|(1 − t⁴) = γa_max³t³; D at that t against D*) and the paradigm's gauge equilibrium (Gm|w|(1 − t⁴) = γa_max³t³):

| |w| | D* | over (γa_max/m)|w| | over (3/2)(2G)^{1/3}(γ/m)^{2/3}|w|^{4/3} | D*'s t | D(t)/D* − 1 | the paradigm's t |
|---|---|---|---|---|---|---|
| 1e-05 | 1.6538e-08 | 0.1985 | 0.9882 | 0.1339 | 0.0e+00 | 0.1063 |
| 0.0001 | 3.4139e-07 | 0.4097 | 0.9468 | 0.2878 | 0.0e+00 | 0.2287 |
| 0.001 | 6.1370e-06 | 0.7364 | 0.7900 | 0.5944 | -1.1e-16 | 0.4840 |
| 0.01 | 7.9402e-05 | 0.9528 | 0.4744 | 0.9101 | 0.0e+00 | 0.8420 |
| 0.1 | 8.2904e-04 | 0.9948 | 0.2299 | 0.9897 | -1.1e-16 | 0.9798 |
| 0.3 | 2.4957e-03 | 0.9983 | 0.1600 | 0.9965 | 0.0e+00 | 0.9931 |

The null kernel's potential against the free masses' at six points of its manifold (the real atoms' state an exact fit of the free law to 1.7e-15): white noise, R_null − R_free from -1.082422e-03 to -1.082422e-03 (−(γ/N)|y|² = -1.082422e-03); slow noise (τ = 10000), from -2.1116e-05 to -9.6002e-06.

Along Part B's instance's white flows (the paradigm's update), at s = 0, 100, 1000, 10000, 100000, 1e+06, 1e+07: the Itô term's length against the flat drift's; for the conserved law the cross term's against the net drift's and the pushes on ‖w‖₁ (d‖w‖₁/ds from the cross term, from the rest); for free masses the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms:

| masses | s | ‖w‖₁/a_max | Itô/flat | cross/net | push: cross, rest | gauge ratio (alive) |
|---|---|---|---|---|---|---|
| conserved | 0 | 0.7894 | 0.037 | 0.0239 | 5.77e-06, 4.82e-04 | — |
| conserved | 100 | 0.8330 | 0.035 | 0.0273 | 1.05e-06, 5.47e-04 | — |
| conserved | 1000 | 0.9537 | 0.187 | 0.186 | 2.34e-06, 1.83e-05 | — |
| conserved | 10000 | 0.9766 | 1.304 | 3.16 | 6.59e-07, -4.12e-07 | — |
| conserved | 100000 | 0.9777 | 1.176 | 1.76 | 5.73e-07, -5.70e-07 | — |
| conserved | 1e+06 | 0.9854 | 1.049 | 2.69 | 4.81e-07, -4.80e-07 | — |
| conserved | 1e+07 | 0.9882 | 1.000 | 8.21e+12 | 3.72e-07, -3.72e-07 | — |
| free | 0 | 0.8016 | 0.028 | — | — | 2.6584–467273.6187 (12) |
| free | 100 | 0.9323 | 0.067 | — | — | 2.5220–17.3204 (12) |
| free | 1000 | 0.9694 | 0.673 | — | — | 1.5763–1.9325 (12) |
| free | 10000 | 0.9389 | 0.934 | — | — | 1.0140–1.0233 (12) |
| free | 100000 | 0.6554 | 0.866 | — | — | 1.0000–1.2235 (12) |
| free | 1e+06 | 0.4670 | 1.000 | — | — | 1.0000–1.0000 (2) |
| free | 1e+07 | 0.4670 | 1.000 | — | — | 1.0000–1.0000 (2) |

## B. The paradigm's noisy scheme against the reduced flow, projected (5 × 12, the 2-sparse instance 4, 64 paths)

From the reduced flow's state at s = 1500, a burn-in of 10000 time units, then a window of 450 in s; every path's state at the window's ends mapped onto M by the noiseless flow (`settle`, Katzenberger's limit map); the projected advance's component along the flow's advance over the same window, relative to it (± s.e. over the paths); the raw paths' mean offset from the projected one at the window's end (largest component; in s.e.); the projected and raw means against the flow's state there. Fresh draws at σ = 0.1; AR(1) draws at τ_c = 10³ (in time units) with σ set by σ²τ_c.

| masses | the draws | h | σ²τ_c/h | steps (burn-in + window) | the flow's advance (largest component) | projected advance / the flow's | the raw mean's offset (z) | projected mean − flow | raw mean − flow | start: projected − flow |
|---|---|---|---|---|---|---|---|---|---|---|
| conserved | fresh | 1 | — | 10000 + 75000 | 8.78e-04 | 0.999 ± 0.001 | 2.5e-04 (0.9) | 7.4e-07 | 2.5e-04 | 3.2e-07 |
| conserved | AR(1), σ²τ_c = 0.9 | 1 | 0.9 | 10000 + 833333 | 8.63e-03 | 0.008 ± 0.171 | 3.3e-03 (2.3) | 7.3e-03 | 8.9e-03 | 8.3e-05 |
| conserved | AR(1), σ²τ_c = 0.225 | 1 | 0.225 | 10000 + 3333333 | 8.63e-03 | 0.808 ± 0.077 | 1.2e-03 (1.9) | 1.5e-03 | 2.1e-03 | 1.8e-05 |
| conserved | AR(1), σ²τ_c = 0.05 | 1 | 0.05 | 10000 + 15000000 | 8.63e-03 | 0.970 ± 0.030 | 8.0e-04 (2.0) | 3.3e-04 | 8.0e-04 | 5.8e-06 |
| conserved | AR(1), σ²τ_c = 0.45 | 0.5 | 0.9 | 20000 + 6666667 | 8.61e-03 | 0.113 ± 0.150 | 1.5e-03 (1.3) | 6.3e-03 | 7.8e-03 | 6.1e-05 |
| conserved | AR(1), σ²τ_c = 0.1 | 0.5 | 0.2 | 20000 + 30000000 | 8.61e-03 | 0.818 ± 0.065 | 8.3e-04 (1.1) | 1.4e-03 | 1.7e-03 | 1.3e-05 |
| free | fresh | 1 | — | 10000 + 75000 | 3.80e-04 | 1.001 ± 0.000 | 2.4e-04 (1.0) | 4.6e-07 | 2.4e-04 | 1.6e-07 |
| free | AR(1), σ²τ_c = 0.9 | 1 | 0.9 | 10000 + 833333 | 7.34e-03 | 0.857 ± 0.001 | 2.5e-03 (1.8) | 9.8e-04 | 3.3e-03 | 1.6e-05 |
| free | AR(1), σ²τ_c = 0.225 | 1 | 0.225 | 10000 + 3333333 | 7.39e-03 | 0.963 ± 0.001 | 1.1e-03 (1.6) | 2.6e-04 | 1.2e-03 | 6.4e-07 |
| free | AR(1), σ²τ_c = 0.05 | 1 | 0.05 | 10000 + 15000000 | 7.40e-03 | 0.991 ± 0.000 | 5.3e-04 (2.1) | 6.1e-05 | 5.9e-04 | 2.4e-07 |
| free | AR(1), σ²τ_c = 0.45 | 0.5 | 0.9 | 20000 + 6666667 | 7.44e-03 | 0.926 ± 0.001 | 1.0e-03 (1.1) | 5.2e-04 | 1.5e-03 | 1.7e-06 |
| free | AR(1), σ²τ_c = 0.1 | 0.5 | 0.2 | 20000 + 30000000 | 7.46e-03 | 0.983 ± 0.000 | 7.7e-04 (1.1) | 1.2e-04 | 7.8e-04 | 1.5e-07 |
| null kernel | fresh | 1 | — | 10000 + 75000 | 8.49e-04 | 0.999 ± 0.001 | 2.6e-04 (0.9) | 1.0e-06 | 2.7e-04 | 3.2e-07 |
| null kernel | AR(1), σ²τ_c = 0.9 | 1 | 0.9 | 10000 + 833333 | 3.19e-03 | -0.340 ± 0.162 | 3.4e-03 (2.4) | 6.9e-03 | 7.5e-03 | 7.8e-05 |
| null kernel | AR(1), σ²τ_c = 0.225 | 1 | 0.225 | 10000 + 3333333 | 3.20e-03 | 0.659 ± 0.094 | 1.1e-03 (1.8) | 1.4e-03 | 2.0e-03 | 2.2e-05 |
| null kernel | AR(1), σ²τ_c = 0.05 | 1 | 0.05 | 10000 + 15000000 | 3.20e-03 | 0.904 ± 0.041 | 5.8e-04 (1.9) | 3.1e-04 | 6.3e-04 | 6.4e-06 |
| null kernel | AR(1), σ²τ_c = 0.45 | 0.5 | 0.9 | 20000 + 6666667 | 3.21e-03 | -0.072 ± 0.159 | 1.2e-03 (1.4) | 5.4e-03 | 6.0e-03 | 6.1e-05 |
| null kernel | AR(1), σ²τ_c = 0.1 | 0.5 | 0.2 | 20000 + 30000000 | 3.21e-03 | 0.833 ± 0.099 | 8.2e-04 (1.1) | 8.0e-04 | 1.4e-03 | 1.4e-05 |

## C. The selection's ends on Gaussian designs (5 × 12, 16 instances per k)

The paradigm's update unless marked; reduced flows to s = 1e+08; the start the noiseless flow's end from equal masses and amplitudes N(0, 0.05²). Each row over the instances whose flow finished (n; a stopped flow is listed below with its last state): means of ‖w‖₁/a_max (the truth's over the same instances), the excess E/m over the real atoms (conserved law and null kernel), the atoms used (|w| > 10⁻³; the largest count in brackets), alive (β > 10⁻³); the fractions recovering the truth (within 10⁻⁴·max|w*|) and within 0.001 of basis pursuit's solution, on the instances where basis pursuit is the truth / where it is not; the largest distance to basis pursuit; the null kernel's mass.

| k | masses | the noise | n | ‖w‖₁/a_max (truth's) | E/m | used (max) | alive | the truth: bp ok / bp fails | at basis pursuit: bp ok / bp fails | largest distance to basis pursuit | the null kernel's mass |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | — | basis pursuit | 16 | 0.462 (0.485) | | | | 13/16 | | | |
| 2 | conserved | the start | 16 | 0.668 (0.485) | 0.332 | 12.0 (12) | 12.0 | 0/13, 0/3 | 0/13, 0/3 | 0.200 | — |
| 2 | conserved | white | 16 | 0.988 (0.485) | 0.012 | 5.9 (6) | 10.1 | 0/13, 0/3 | 0/13, 0/3 | 0.363 | — |
| 2 | conserved | slow, τ_c = 1e+03 | 16 | 0.990 (0.485) | 0.010 | 5.9 (6) | 7.9 | 0/13, 0/3 | 0/13, 0/3 | 0.335 | — |
| 2 | conserved | slow, τ_c = 1e+04 | 16 | 0.994 (0.485) | 0.006 | 6.0 (6) | 7.2 | 0/13, 0/3 | 0/13, 0/3 | 0.475 | — |
| 2 | conserved | white, the flat update | 16 | 0.975 (0.485) | 0.025 | 12.0 (12) | 12.0 | 0/13, 0/3 | 0/13, 0/3 | 0.232 | — |
| 2 | free | the start | 16 | 0.676 (0.485) | — | 11.9 (12) | 12.0 | 0/13, 0/3 | 0/13, 0/3 | 0.196 | — |
| 2 | free | white | 16 | 0.468 (0.485) | — | 2.6 (5) | 2.7 | 13/13, 0/3 | 13/13, 2/3 | 0.108 | — |
| 2 | free | slow, τ_c = 1e+03 | 16 | 0.473 (0.485) | — | 2.5 (5) | 2.9 | 12/13, 1/3 | 12/13, 1/3 | 0.270 | — |
| 2 | free | slow, τ_c = 1e+04 | 16 | 0.471 (0.485) | — | 2.4 (5) | 2.6 | 13/13, 1/3 | 13/13, 1/3 | 0.270 | — |
| 2 | free | white, the flat update | 16 | 0.463 (0.485) | — | 3.5 (6) | 6.3 | 2/13, 0/3 | 5/13, 2/3 | 0.015 | — |
| 2 | null kernel | the start | 16 | 0.661 (0.485) | 0.339 | 12.0 (12) | 12.0 | 0/13, 0/3 | 0/13, 0/3 | 0.202 | 0.77 |
| 2 | null kernel | white | 16 | 0.486 (0.485) | 0.514 | 2.8 (5) | 2.9 | 11/13, 1/3 | 11/13, 1/3 | 0.270 | 6.11 |
| 2 | null kernel | slow, τ_c = 1e+03 | 16 | 0.610 (0.485) | 0.390 | 3.7 (5) | 3.8 | 7/13, 0/3 | 7/13, 0/3 | 0.421 | 4.61 |
| 2 | null kernel | slow, τ_c = 1e+04 | 16 | 0.627 (0.485) | 0.373 | 3.6 (5) | 3.9 | 6/13, 1/3 | 6/13, 0/3 | 0.421 | 4.41 |
| 3 | — | basis pursuit | 16 | 0.578 (0.671) | | | | 6/16 | | | |
| 3 | conserved | the start | 16 | 0.731 (0.671) | 0.269 | 11.9 (12) | 12.0 | 0/6, 0/10 | 0/6, 0/10 | 0.315 | — |
| 3 | conserved | white | 16 | 0.988 (0.671) | 0.012 | 5.9 (6) | 10.9 | 0/6, 0/10 | 0/6, 0/10 | 0.462 | — |
| 3 | conserved | slow, τ_c = 1e+03 | 16 | 0.992 (0.671) | 0.008 | 5.9 (6) | 8.2 | 0/6, 0/10 | 0/6, 0/10 | 0.463 | — |
| 3 | conserved | slow, τ_c = 1e+04 | 16 | 0.995 (0.671) | 0.005 | 5.9 (6) | 7.2 | 0/6, 0/10 | 0/6, 0/10 | 0.444 | — |
| 3 | conserved | white, the flat update | 16 | 0.974 (0.671) | 0.026 | 11.9 (12) | 12.0 | 0/6, 0/10 | 0/6, 0/10 | 0.353 | — |
| 3 | free | the start | 16 | 0.763 (0.671) | — | 11.9 (12) | 12.0 | 0/6, 0/10 | 0/6, 0/10 | 0.318 | — |
| 3 | free | white | 16 | 0.581 (0.671) | — | 4.1 (5) | 5.0 | 6/6, 1/10 | 6/6, 3/10 | 0.214 | — |
| 3 | free | slow, τ_c = 1e+03 | 16 | 0.582 (0.671) | — | 4.0 (5) | 5.0 | 6/6, 1/10 | 6/6, 3/10 | 0.214 | — |
| 3 | free | slow, τ_c = 1e+04 | 16 | 0.591 (0.671) | — | 3.9 (5) | 4.9 | 6/6, 2/10 | 6/6, 2/10 | 0.234 | — |
| 3 | free | white, the flat update | 16 | 0.578 (0.671) | — | 5.1 (7) | 7.8 | 0/6, 0/10 | 1/6, 4/10 | 0.010 | — |
| 3 | null kernel | the start | 16 | 0.726 (0.671) | 0.274 | 11.9 (12) | 12.0 | 0/6, 0/10 | 0/6, 0/10 | 0.311 | 0.66 |
| 3 | null kernel | white | 16 | 0.587 (0.671) | 0.413 | 4.2 (5) | 5.4 | 5/6, 1/10 | 5/6, 2/10 | 0.222 | 4.87 |
| 3 | null kernel | slow, τ_c = 1e+03 | 16 | 0.655 (0.671) | 0.345 | 3.8 (5) | 4.7 | 6/6, 4/10 | 6/6, 0/10 | 0.253 | 4.08 |
| 3 | null kernel | slow, τ_c = 1e+04 | 16 | 0.743 (0.671) | 0.257 | 4.1 (5) | 5.0 | 4/6, 3/10 | 4/6, 0/10 | 0.410 | 3.02 |

The null kernel against free masses, instance by instance (the paradigm's update, the same colour): the instances whose ends differ (largest coefficient's distance above 0.001); the instances where the null kernel misses the truth while basis pursuit and free masses recover it:

| k | the noise | ends differ | the null kernel misses the truth where both recover it |
|---|---|---|---|
| 2 | white | 4 of 16 (0, 6, 10, 13) | 2 of 13 (0, 10) |
| 2 | slow, τ_c = 1e+03 | 9 of 16 (0, 3, 6, 8, 10, 11, 13, 14, 15) | 5 of 12 (3, 8, 10, 11, 15) |
| 2 | slow, τ_c = 1e+04 | 9 of 16 (0, 2, 3, 8, 10, 11, 13, 14, 15) | 7 of 13 (0, 2, 3, 8, 10, 11, 15) |
| 3 | white | 5 of 16 (0, 2, 4, 12, 15) | 1 of 6 (0) |
| 3 | slow, τ_c = 1e+03 | 9 of 16 (1, 2, 3, 4, 5, 9, 10, 12, 15) | 0 of 6 (—) |
| 3 | slow, τ_c = 1e+04 | 11 of 16 (1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 15) | 2 of 6 (8, 11) |

Free masses: the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms at the ends (1 at the paradigm's gauge equilibrium, 1/2 at D*'s split), median [range] over the atoms of all finished flows; the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|), that minimiser against basis pursuit's solution (median, largest), and the flat ends on the truth:

| k | white | slow, τ_c = 1e+03 | slow, τ_c = 1e+04 | white, the flat update | flat end − D* minimiser | D* minimiser − basis pursuit | flat ends on the truth |
|---|---|---|---|---|---|---|---|
| 2 | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 0.5000 [0.5000, 0.5000] | 1.1e-04, 3.3e-04 | 1.1e-03, 1.5e-02 | 2 of 16 |
| 3 | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 0.5000 [0.5000, 0.5000] | 1.3e-04, 9.3e-03 | 2.0e-03, 5.7e-03 | 0 of 16 |

## D. The selection's ends on Gaussian kernel dictionaries (16 centres, 8 samples on [0, 1], width 0.08, 8 instances per k)

The paradigm's update unless marked; reduced flows to s = 1e+08; the start the noiseless flow's end from equal masses and amplitudes N(0, 0.05²). Each row over the instances whose flow finished (n; a stopped flow is listed below with its last state): means of ‖w‖₁/a_max (the truth's over the same instances), the excess E/m over the real atoms (conserved law and null kernel), the atoms used (|w| > 10⁻³; the largest count in brackets), alive (β > 10⁻³); the fractions recovering the truth (within 10⁻⁴·max|w*|) and within 0.001 of basis pursuit's solution, on the instances where basis pursuit is the truth / where it is not; the largest distance to basis pursuit; the null kernel's mass.

| k | masses | the noise | n | ‖w‖₁/a_max (truth's) | E/m | used (max) | alive | the truth: bp ok / bp fails | at basis pursuit: bp ok / bp fails | largest distance to basis pursuit | the null kernel's mass |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | — | basis pursuit | 8 | 0.428 (0.453) | | | | 7/8 | | | |
| 2 | conserved | the start | 8 | 0.590 (0.453) | 0.410 | 13.9 (16) | 16.0 | 0/7, 0/1 | 0/7, 0/1 | 0.171 | — |
| 2 | conserved | white | 8 | 0.981 (0.453) | 0.019 | 7.8 (9) | 15.5 | 0/7, 0/1 | 0/7, 0/1 | 0.448 | — |
| 2 | conserved | slow, τ_c = 1e+03 | 8 | 0.983 (0.453) | 0.017 | 7.6 (9) | 15.2 | 0/7, 0/1 | 0/7, 0/1 | 0.449 | — |
| 2 | conserved | slow, τ_c = 1e+04 | 8 | 0.986 (0.453) | 0.014 | 7.5 (9) | 14.1 | 0/7, 0/1 | 0/7, 0/1 | 0.449 | — |
| 2 | conserved | white, the flat update | 8 | 0.966 (0.453) | 0.034 | 16.0 (16) | 16.0 | 0/7, 0/1 | 0/7, 0/1 | 0.192 | — |
| 2 | free | the start | 8 | 0.611 (0.453) | — | 14.1 (16) | 16.0 | 0/7, 0/1 | 0/7, 0/1 | 0.176 | — |
| 2 | free | white | 8 | 0.428 (0.453) | — | 2.2 (4) | 2.8 | 7/7, 0/1 | 7/7, 1/1 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+03 | 8 | 0.428 (0.453) | — | 2.2 (4) | 2.9 | 7/7, 0/1 | 7/7, 1/1 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+04 | 8 | 0.453 (0.453) | — | 2.0 (2) | 2.0 | 7/7, 1/1 | 7/7, 0/1 | 0.298 | — |
| 2 | free | white, the flat update | 8 | 0.428 (0.453) | — | 2.6 (5) | 5.6 | 2/7, 0/1 | 6/7, 1/1 | 0.003 | — |
| 2 | null kernel | the start | 8 | 0.585 (0.453) | 0.415 | 13.9 (16) | 16.0 | 0/7, 0/1 | 0/7, 0/1 | 0.169 | 0.50 |
| 2 | null kernel | white | 8 | 0.428 (0.453) | 0.572 | 2.2 (4) | 2.8 | 7/7, 0/1 | 7/7, 1/1 | 0.000 | 6.82 |
| 2 | null kernel | slow, τ_c = 1e+03 | 8 | 0.453 (0.453) | 0.547 | 2.0 (2) | 2.0 | 7/7, 1/1 | 7/7, 0/1 | 0.298 | 6.53 |
| 2 | null kernel | slow, τ_c = 1e+04 | 8 | 0.453 (0.453) | 0.547 | 2.0 (2) | 2.0 | 7/7, 1/1 | 7/7, 0/1 | 0.298 | 6.53 |
| 3 | — | basis pursuit | 8 | 0.559 (0.654) | | | | 3/8 | | | |
| 3 | conserved | the start | 8 | 0.749 (0.654) | 0.251 | 15.8 (16) | 16.0 | 0/3, 0/5 | 0/3, 0/5 | 0.180 | — |
| 3 | conserved | white | 8 | 0.983 (0.654) | 0.017 | 7.8 (9) | 14.0 | 0/3, 0/5 | 0/3, 0/5 | 0.371 | — |
| 3 | conserved | slow, τ_c = 1e+03 | 8 | 0.986 (0.654) | 0.014 | 7.5 (9) | 13.1 | 0/3, 0/5 | 0/3, 0/5 | 0.362 | — |
| 3 | conserved | slow, τ_c = 1e+04 | 8 | 0.990 (0.654) | 0.010 | 7.5 (9) | 9.9 | 0/3, 0/5 | 0/3, 0/5 | 0.362 | — |
| 3 | conserved | white, the flat update | 8 | 0.966 (0.654) | 0.034 | 16.0 (16) | 16.0 | 0/3, 0/5 | 0/3, 0/5 | 0.200 | — |
| 3 | free | the start | 8 | 0.800 (0.654) | — | 15.8 (16) | 16.0 | 0/3, 0/5 | 0/3, 0/5 | 0.194 | — |
| 3 | free | white | 8 | 0.563 (0.654) | — | 4.4 (7) | 7.0 | 3/3, 1/5 | 3/3, 1/5 | 0.329 | — |
| 3 | free | slow, τ_c = 1e+03 | 8 | 0.566 (0.654) | — | 4.0 (5) | 6.8 | 3/3, 1/5 | 3/3, 2/5 | 0.329 | — |
| 3 | free | slow, τ_c = 1e+04 | 8 | 0.591 (0.654) | — | 3.8 (5) | 6.0 | 3/3, 2/5 | 3/3, 1/5 | 0.329 | — |
| 3 | free | white, the flat update | 7 | 0.529 (0.637) | — | 6.0 (9) | 10.0 | 0/2, 0/5 | 1/2, 2/5 | 0.009 | — |
| 3 | null kernel | the start | 8 | 0.740 (0.654) | 0.260 | 16.0 (16) | 16.0 | 0/3, 0/5 | 0/3, 0/5 | 0.177 | 0.43 |
| 3 | null kernel | white | 8 | 0.571 (0.654) | 0.429 | 4.5 (8) | 6.5 | 3/3, 1/5 | 3/3, 1/5 | 0.329 | 5.05 |
| 3 | null kernel | slow, τ_c = 1e+03 | 8 | 0.599 (0.654) | 0.401 | 4.1 (8) | 5.8 | 3/3, 2/5 | 3/3, 0/5 | 0.329 | 4.73 |
| 3 | null kernel | slow, τ_c = 1e+04 | 8 | 0.599 (0.654) | 0.401 | 4.1 (8) | 6.4 | 3/3, 2/5 | 3/3, 0/5 | 0.329 | 4.73 |

The null kernel against free masses, instance by instance (the paradigm's update, the same colour): the instances whose ends differ (largest coefficient's distance above 0.001); the instances where the null kernel misses the truth while basis pursuit and free masses recover it:

| k | the noise | ends differ | the null kernel misses the truth where both recover it |
|---|---|---|---|
| 2 | white | 0 of 8 (—) | 0 of 7 (—) |
| 2 | slow, τ_c = 1e+03 | 1 of 8 (0) | 0 of 7 (—) |
| 2 | slow, τ_c = 1e+04 | 0 of 8 (—) | 0 of 7 (—) |
| 3 | white | 1 of 8 (3) | 0 of 3 (—) |
| 3 | slow, τ_c = 1e+03 | 2 of 8 (0, 3) | 0 of 3 (—) |
| 3 | slow, τ_c = 1e+04 | 1 of 8 (3) | 0 of 3 (—) |

Free masses: the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms at the ends (1 at the paradigm's gauge equilibrium, 1/2 at D*'s split), median [range] over the atoms of all finished flows; the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|), that minimiser against basis pursuit's solution (median, largest), and the flat ends on the truth:

| k | white | slow, τ_c = 1e+03 | slow, τ_c = 1e+04 | white, the flat update | flat end − D* minimiser | D* minimiser − basis pursuit | flat ends on the truth |
|---|---|---|---|---|---|---|---|
| 2 | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 0.5000 [0.5000, 0.5000] | 5.9e-05, 2.6e-04 | 4.0e-04, 2.8e-03 | 2 of 8 |
| 3 | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] | 0.5000 [0.5000, 0.5000] | 5.7e-04, 7.6e-03 | 1.3e-03, 1.1e-02 | 0 of 7 |

Flows stopped (over the 900-s limit, or failed), left out of the rows: k = 3, instance 1, free, flat: over the limit at s = 2.83e+07 (its last state: ‖w‖₁/a_max 0.776, the truth recovered (distance 4.8e-14), distance to basis pursuit 4.8e-14, 3 atoms used, 3 alive).

## E. The conserved law's white end against the representability boundary (5 × 12, 2-sparse truths, 12 instances; those whose basis-pursuit solution is not representable, ‖w_bp‖₁ ≥ a_max, left out)

(iii)'s bound E/m ≤ (a_max/m)(γmΣ|φ_i|⁻²max|φ_i|²/G)^{1/2} (vacuous at 1 and above); the largest amplitude |a|/a_max at the end (median over the instances); the saturation threshold γa_max³/(Gm) against the median used coefficient |w_i| at the end:

| a_max | instances | ‖w‖₁/a_max at the start | at the end (min) | E/m at the end | the bound | max|a|/a_max | γa_max³/(Gm), median used |w_i| | the atoms used, alive | ‖w_bp‖₁/a_max |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 7 | 0.987 | 0.998 (0.997) | 0.0021 | 0.158 | 0.999 | 0.0010, 0.045 | 6.0, 6.9 | 0.870 |
| 1 | 12 | 0.697 | 0.988 (0.984) | 0.0119 | 0.316 | 0.994 | 0.0083, 0.146 | 5.9, 9.9 | 0.476 |
| 2 | 12 | 0.327 | 0.942 (0.930) | 0.0585 | 0.632 | 0.971 | 0.0667, 0.257 | 6.9, 12.0 | 0.238 |
| 4 | 12 | 0.161 | 0.653 (0.582) | 0.3469 | 1.265 | 0.862 | 0.5333, 0.038 | 11.5, 12.0 | 0.119 |

## F. Katzenberger's range along the flows, Q's sign maps, the null kernel at a finite step (referee round 146)

Along Part B's instance's white flows (the paradigm's update; to s = 1e+08 with Part C's pieces): K's spectrum (the normal rates, decreasing), the atoms alive (β > 10⁻³), the dying atoms' masses (β ≤ 10⁻³: smallest, largest), Q/c on the dying and on the alive atoms (free masses; c = γa_max²/(4G)), the fastest relative mass rate max|d log β_i/ds| over the atoms with β_i > 10⁻¹², and λ_min(K) over that rate per unit time at Part B's noise (σ = 0.1, h = 1: ds/dt = 6.0e-03) — Katzenberger's limit needs it large —, with the σ at which it would be 10 (the ratio ∝ 1/σ² at fixed h):

| masses | s | ‖w‖₁/a_max | K's spectrum | alive | dying masses: min, max | Q/c dying: min, max | Q/c alive: min, max | max\|d log β/ds\| | λ_min(K)/rate at σ = 0.1 | σ for a ratio of 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 0.8016 | 3.3e-02, 2.0e-02, 1.3e-02, 8.7e-03, 3.7e-03 | 12 | — | — | 38.9808, 38.9999 | 3.34e-03 | 187 | 0.43 |
| free | 1000 | 0.9694 | 5.8e-03, 5.3e-03, 4.6e-03, 2.6e-03, 1.6e-03 | 12 | — | — | 13.1523, 20.9737 | 4.08e-05 | 6.71e+03 | 2.6 |
| free | 1500 | 0.9679 | 5.7e-03, 5.3e-03, 4.6e-03, 2.5e-03, 1.6e-03 | 12 | — | — | 10.1605, 16.5789 | 2.98e-05 | 9.04e+03 | 3 |
| free | 3000 | 0.9623 | 5.6e-03, 5.2e-03, 4.5e-03, 2.5e-03, 1.6e-03 | 12 | — | — | 5.3869, 8.9972 | 1.99e-05 | 1.33e+04 | 3.6 |
| free | 10000 | 0.9389 | 5.4e-03, 5.1e-03, 4.4e-03, 2.4e-03, 1.5e-03 | 12 | — | — | 0.4841, 0.8078 | 1.50e-05 | 1.7e+04 | 4.1 |
| free | 30000 | 0.8747 | 5.2e-03, 4.8e-03, 4.0e-03, 2.1e-03, 1.4e-03 | 12 | — | — | -0.0551, -0.0141 | 1.90e-05 | 1.19e+04 | 3.5 |
| free | 100000 | 0.6554 | 5.3e-03, 3.6e-03, 2.6e-03, 9.7e-04, 6.0e-04 | 12 | — | — | -0.5073, -0.0110 | 5.55e-05 | 1.79e+03 | 1.3 |
| free | 200000 | 0.4809 | 6.2e-03, 3.0e-03, 3.4e-04, 4.1e-05, 1.3e-05 | 12 | — | — | -0.8984, -0.0076 | 5.62e-05 | 37.7 | 0.19 |
| free | 300000 | 0.4671 | 6.3e-03, 3.1e-03, 7.8e-07, 2.3e-07, 1.5e-07 | 12 | — | — | -0.9449, -0.0074 | 1.38e-05 | 1.76 | 0.042 |
| free | 500000 | 0.4670 | 6.3e-03, 3.1e-03, 5.7e-08, 2.8e-08, 2.0e-08 | 6 | 7.28e-04, 9.66e-04 | -0.97126, -0.96279 | -0.9603, -0.0074 | 3.72e-06 | 0.874 | 0.03 |
| free | 1e+06 | 0.4670 | 6.3e-03, 3.1e-03, 7.3e-09, 4.3e-09, 3.0e-09 | 2 | 3.27e-04, 3.98e-04 | -0.98692, -0.98460 | -0.0110, -0.0074 | 1.30e-06 | 0.379 | 0.019 |
| free | 3e+06 | 0.4670 | 6.3e-03, 3.1e-03, 5.9e-10, 3.7e-10, 2.5e-10 | 2 | 1.03e-04, 1.12e-04 | -0.99589, -0.99568 | -0.0110, -0.0074 | 3.60e-07 | 0.118 | 0.011 |
| free | 1e+07 | 0.4670 | 6.3e-03, 3.1e-03, 5.3e-11, 3.3e-11, 2.2e-11 | 2 | 3.06e-05, 3.35e-05 | -0.99879, -0.99877 | -0.0110, -0.0074 | 1.02e-07 | 0.0361 | 0.006 |
| free | 3e+07 | 0.4670 | 6.3e-03, 3.1e-03, 7.4e-12, 4.3e-12, 2.7e-12 | 2 | 1.01e-05, 1.26e-05 | -0.99960, -0.99959 | -0.0110, -0.0074 | 3.35e-08 | 0.0137 | 0.0037 |
| free | 1e+08 | 0.4670 | 6.3e-03, 3.1e-03, 1.4e-12, 6.9e-13, 3.9e-13 | 2 | 3.01e-06, 5.64e-06 | -0.99988, -0.99986 | -0.0110, -0.0074 | 9.97e-09 | 0.00652 | 0.0026 |
| conserved | 0 | 0.7894 | 1.8e-02, 9.6e-03, 7.0e-03, 4.3e-03, 1.5e-03 | 12 | — | — | — | 1.04e-03 | 243 | 0.49 |
| conserved | 1000 | 0.9537 | 5.6e-03, 4.7e-03, 3.1e-03, 2.6e-03, 9.5e-04 | 12 | — | — | — | 4.85e-05 | 3.27e+03 | 1.8 |
| conserved | 1500 | 0.9611 | 5.5e-03, 4.7e-03, 3.1e-03, 2.6e-03, 9.1e-04 | 12 | — | — | — | 2.77e-05 | 5.5e+03 | 2.3 |
| conserved | 3000 | 0.9698 | 5.5e-03, 4.6e-03, 3.1e-03, 2.6e-03, 8.7e-04 | 12 | — | — | — | 1.01e-05 | 1.44e+04 | 3.8 |
| conserved | 10000 | 0.9766 | 5.5e-03, 4.7e-03, 3.1e-03, 2.6e-03, 8.4e-04 | 12 | — | — | — | 4.11e-06 | 3.42e+04 | 5.8 |
| conserved | 30000 | 0.9776 | 5.4e-03, 4.7e-03, 3.0e-03, 2.7e-03, 8.2e-04 | 12 | — | — | — | 6.13e-06 | 2.23e+04 | 4.7 |
| conserved | 100000 | 0.9777 | 5.3e-03, 5.0e-03, 2.9e-03, 2.8e-03, 6.9e-04 | 12 | — | — | — | 2.52e-05 | 4.56e+03 | 2.1 |
| conserved | 200000 | 0.9796 | 5.8e-03, 5.3e-03, 3.0e-03, 2.1e-03, 5.2e-04 | 12 | — | — | — | 4.85e-05 | 1.79e+03 | 1.3 |
| conserved | 300000 | 0.9820 | 5.9e-03, 5.5e-03, 2.9e-03, 1.9e-03, 4.7e-04 | 12 | — | — | — | 9.32e-06 | 8.42e+03 | 2.9 |
| conserved | 500000 | 0.9824 | 6.0e-03, 5.6e-03, 2.8e-03, 2.0e-03, 3.7e-04 | 12 | — | — | — | 8.53e-06 | 7.29e+03 | 2.7 |
| conserved | 1e+06 | 0.9854 | 7.4e-03, 5.2e-03, 2.3e-03, 1.8e-03, 2.3e-04 | 12 | — | — | — | 1.01e-05 | 3.84e+03 | 2 |
| conserved | 3e+06 | 0.9882 | 1.0e-02, 3.0e-03, 2.3e-03, 1.3e-03, 2.2e-04 | 12 | — | — | — | 2.38e-09 | 1.55e+07 | 1.2e+02 |
| conserved | 1e+07 | 0.9882 | 1.0e-02, 3.0e-03, 2.3e-03, 1.3e-03, 2.2e-04 | 12 | — | — | — | 1.15e-17 | 3.23e+15 | 1.8e+06 |
| conserved | 3e+07 | 0.9882 | 1.0e-02, 3.0e-03, 2.3e-03, 1.3e-03, 2.2e-04 | 12 | — | — | — | 3.49e-19 | 1.06e+17 | 1e+07 |
| conserved | 1e+08 | 0.9882 | 1.0e-02, 3.0e-03, 2.3e-03, 1.3e-03, 2.2e-04 | 12 | — | — | — | 8.38e-19 | 4.42e+16 | 6.6e+06 |

Free masses: the ratio is 187 at s = 0 and 9.04e+03 at Part B's s = 1500; it falls below 10 at s ≈ 2.4e+05 and below 1 at s ≈ 4.5e+05 (log–log interpolation), and is 0.0065 at s = 1e+08; at s = 1e+08 2 of 12 atoms are alive, the smallest mass 3.01e-06, λ_min(K) = 3.9e-13; the end is Part C's (k = 2, instance 4, free, white) to 0.0e+00.

Conserved masses: the ratio is 243 at s = 0 and 5.5e+03 at Part B's s = 1500; it stays above 1 (at least 243); at s = 1e+08 12 of 12 atoms are alive, the smallest mass 1.40e-03, λ_min(K) = 2.2e-04.

The paradigm's scheme against the free law's reduced flow, projected as in B (64 paths, fresh draws, h = 1, the flow at the step's filter), from the flow's state at s0: a burn-in of 10/λ_min(K) time units, a window of 2000 in s, the limit map (`settle`) run for max(2·10⁵, 100/λ_min(K)); λ_min(K) at s0 and the ratio there at the run's σ:

| σ | s0 | λ_min(K) at s0 | λ_min(K)/rate at s0 | steps (burn-in + window) | the flow's advance (largest component) | projected advance / the flow's | the raw mean's offset (z) | projected mean − flow | raw mean − flow | start: projected − flow |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 100000 | 5.9e-04 | 1.78e+03 | 16810 + 333333 | 1.65e-03 | 1.000 ± 0.000 | 9.8e-04 (2.2) | 9.2e-07 | 9.8e-04 | 1.3e-07 |
| 0.05 | 100000 | 5.9e-04 | 7.11e+03 | 16810 + 1333333 | 1.65e-03 | 1.000 ± 0.000 | 4.5e-04 (2.5) | 2.4e-07 | 4.5e-04 | 2.4e-08 |
| 0.1 | 200000 | 1.2e-05 | 37.3 | 801432 + 333333 | 9.61e-04 | 0.999 ± 0.001 | 8.5e-04 (3.2) | 2.2e-06 | 8.5e-04 | 1.5e-06 |
| 0.05 | 200000 | 1.2e-05 | 149 | 801432 + 1333333 | 1.06e-03 | 1.000 ± 0.000 | 6.3e-04 (2.0) | 3.9e-07 | 6.3e-04 | 1.4e-07 |

Free masses, Q's drift dQ_i/ds = (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)) (`q_drift`, white) on the free law's curves Q/c, at 199 amplitudes t = a/a_max from 0.005 to 0.995 within each curve's range (β > 0): its sign along the curve for the flat update (κ = 1) and the paradigm's (κ = 2), and its zero (bisection on `q_drift`) against the closed form — κ = 2: t = 1/κ′ = ((c + Q)/(c − Q))^{1/2} on −c < Q < 0; κ = 1: t = ((−Q − c)/(c − Q))^{1/4} on a branch:

| Q/c | the curve | the flat update (κ = 1) | its zero (closed form) | the paradigm's (κ = 2) | its zero (closed form) |
|---|---|---|---|---|---|
| 1 | a free curve, Q ≥ 0 | − on [0.005, 0.995] | — | − on [0.005, 0.995] | — |
| 0 | a free curve, Q ≥ 0 | − on [0.005, 0.995] | — | − on [0.005, 0.995] | — |
| -0.5 | a free curve, −c < Q < 0 | − on [0.005, 0.995] | — | − on [0.005, 0.575], + on [0.580, 0.995] | 0.577350 (0.577350) |
| -0.9 | a free curve, −c < Q < 0 | − on [0.005, 0.995] | — | − on [0.005, 0.225], + on [0.230, 0.995] | 0.229416 (0.229416) |
| -1 | the critical curve | − on [0.005, 0.995] | — | + on [0.005, 0.995] | — |
| -1.1 | a branch, Q < −c | + on [0.220, 0.465], − on [0.470, 0.995] | 0.467138 (0.467138) | + on [0.220, 0.995] | — |
| -2 | a branch, Q < −c | + on [0.580, 0.755], − on [0.760, 0.995] | 0.759836 (0.759836) | + on [0.580, 0.995] | — |
| -10 | a branch, Q < −c | + on [0.905, 0.950], − on [0.955, 0.995] | 0.951070 (0.951070) | + on [0.905, 0.995] | — |

Where each update's Q-drift vanishes on the orbit β|a| = m|w|: D*'s split (the flat update, 2Gm|w|(1 − t⁴) = γa_max³t³) and the paradigm's gauge equilibrium (Gm|w|(1 − t⁴) = γa_max³t³) — Q/c on states built there (`gauge_invariant`) against the closed forms −(1 + t⁴)/(1 − t⁴) (a branch for t > 0) and −(1 − t²)/(1 + t²); at the paradigm's equilibrium the Itô term's mass component over the flat drift's before the projection (`ito_term`, `grad`) against (1 + t²)/(3 − t²):

| \|w\| | D*'s split: t | Q/c | −(1 + t⁴)/(1 − t⁴) | the paradigm's: t | Q/c | −(1 − t²)/(1 + t²) | Itô/flat | (1 + t²)/(3 − t²) |
|---|---|---|---|---|---|---|---|---|
| 1e-05 | 0.1339 | -1.0006 | -1.0006 | 0.1063 | -0.9777 | -0.9777 | 0.3384 | 0.3384 |
| 0.0001 | 0.2878 | -1.0138 | -1.0138 | 0.2287 | -0.9006 | -0.9006 | 0.3570 | 0.3570 |
| 0.001 | 0.5944 | -1.2853 | -1.2853 | 0.4840 | -0.6204 | -0.6204 | 0.4463 | 0.4463 |
| 0.01 | 0.9101 | -5.3683 | -5.3683 | 0.8420 | -0.1703 | -0.1703 | 0.7459 | 0.7459 |
| 0.1 | 0.9897 | -48.5077 | -48.5077 | 0.9798 | -0.0204 | -0.0204 | 0.9608 | 0.9608 |
| 0.3 | 0.9965 | -144.5026 | -144.5026 | 0.9931 | -0.0069 | -0.0069 | 0.9864 | 0.9864 |

The null kernel's potential against the free masses' at Part A's six points (Part A's points replayed: its white and slow values reproduced to 0.0e+00): R_null − R_free from -1.082422e-03 to -1.082422e-03 under white noise (spread 6.9e-17); under fresh draws at a finite step, Ψ = −(2/h)log(1 − hκ/2) (Ψ′ = 2/(2 − hκ)), h = 1 from -1.0914e-03 to -1.0885e-03 (spread 2.86e-06); h = 10 from -1.1844e-03 to -1.1513e-03 (spread 3.31e-05); slow noise (τ = 10000) from -2.1116e-05 to -9.6002e-06 (spread 1.15e-05).

Numbers are produced by `experiments/exp_mass_dictionary.py`; nothing is edited by hand.
