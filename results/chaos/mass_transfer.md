# Proposition D‴.34: mass transfer — the replicator weight law with the cohesion leak (`exp_mass_transfer.py`)

The instances of `exp_cohesion_leak.py` (σ₀ = 0.00707439); the weighted network f̂ = (1/m)Σβ_iΦ̃(·;ϑ_i), Σβ = m, the replicator law dβ_i/ds = γβ_i(u_i − λ), u_i = ⟨rΦ̃(·;ϑ_i)⟩, beside the transport law with the cohesion leak (γ = 0: Proposition D‴.32's loop); masses are Σβ/m per cluster (single linkage at 0.25), the classes listed by centre.

## A. Fixed positions (the five-bump atoms, one particle per class, the amplitudes adapted; 40000 slow units)

| γ | amplitude scale | Σ|A_k|/a_max | fit error at s = 4000, 20000, 40000 | masses at the end | needed |A_k|/a_max | |a|/a_max at the end |
|---|---|---|---|---|---|---|
| 0 | 1 | 0.875 | 0.391, 0.388, 0.388 | 0.150, 0.300, 0.200, 0.250, 0.100 | 0.225, 0.150, 0.175, 0.125, 0.200 | 0.999, 0.425, 0.826, 0.383, 0.999 |
| 0.001 | 1 | 0.875 | 0.129, 0.00594, 2.91e-4 | 0.226, 0.215, 0.178, 0.179, 0.201 | 0.225, 0.150, 0.175, 0.125, 0.200 | 0.994, 0.698, 0.980, 0.697, 0.995 |
| 0.01 | 1 | 0.875 | 0.00238, 1.26e-9, 2.99e-14 | 0.234, 0.199, 0.191, 0.168, 0.208 | 0.225, 0.150, 0.175, 0.125, 0.200 | 0.960, 0.755, 0.917, 0.745, 0.961 |
| 0 | 1.25 | 1.094 | 0.496, 0.495, 0.494 | 0.150, 0.300, 0.200, 0.250, 0.100 | 0.281, 0.188, 0.219, 0.156, 0.250 | 1.000, 0.491, 0.991, 0.444, 1.000 |
| 0.001 | 1.25 | 1.094 | 0.181, 0.0739, 0.071 | 0.270, 0.164, 0.194, 0.133, 0.239 | 0.281, 0.188, 0.219, 0.156, 0.250 | 0.998, 0.994, 0.997, 0.992, 0.998 |
| 0.01 | 1.25 | 1.094 | 0.0872, 0.0735, 0.0711 | 0.271, 0.164, 0.194, 0.132, 0.239 | 0.281, 0.188, 0.219, 0.156, 0.250 | 0.997, 0.995, 0.996, 0.994, 0.997 |

## B. The five-bump target from broad starts (m = 60, the annealed Gaussian cohesion leak)

| γ | seed | clusters after R = 2 | fit error | clusters at the end (> 1% of the mass) | fit error | masses | centres | particles | |a|/a_max (max) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 6 | 0.12 | 6 (6) | 0.118 | 0.067, 0.183, 0.067, 0.150, 0.117, 0.417 | -1.879, -1.581, -0.801, +0.015, +0.828, +1.626 | 4, 11, 4, 9, 7, 25 | 0.998 |
| 0 | 1 | 10 | 0.0484 | 9 (9) | 0.045 | 0.133, 0.200, 0.100, 0.183, 0.017, 0.100, 0.017, 0.133, 0.117 | -1.862, -1.498, -0.792, -0.015, +0.290, +0.809, +0.973, +1.505, +1.752 | 8, 12, 6, 11, 1, 6, 1, 8, 7 | 0.995 |
| 0 | 2 | 7 | 0.0589 | 8 (8) | 0.0497 | 0.133, 0.233, 0.100, 0.167, 0.100, 0.017, 0.133, 0.117 | -1.882, -1.517, -0.791, +0.008, +0.802, +1.005, +1.478, +1.781 | 8, 14, 6, 10, 6, 1, 8, 7 | 0.995 |
| 0.001 | 0 | 6 | 0.0367 | 5 (5) | 0.0139 | 0.229, 0.131, 0.174, 0.112, 0.353 | -1.623, -0.811, +0.003, +0.818, +1.621 | 15, 4, 9, 7, 25 | 0.994 |
| 0.001 | 1 | 8 | 0.0177 | 7 (7) | 0.00761 | 0.121, 0.178, 0.131, 0.189, 0.105, 0.045, 0.231 | -1.824, -1.505, -0.796, +0.009, +0.821, +1.024, +1.611 | 8, 12, 6, 10, 6, 3, 15 | 0.987 |
| 0.001 | 2 | 8 | 0.0174 | 6 (6) | 0.00862 | 0.121, 0.208, 0.132, 0.191, 0.115, 0.233 | -1.836, -1.520, -0.795, +0.005, +0.819, +1.617 | 8, 14, 6, 10, 7, 15 | 0.982 |
| 0.01 | 0 | 6 | 0.00651 | 5 (5) | 0.00135 | 0.258, 0.155, 0.210, 0.134, 0.242 | -1.603, -0.799, +0.009, +0.801, +1.593 | 15, 4, 9, 7, 25 | 0.961 |
| 0.01 | 1 | 5 | 0.0043 | 5 (5) | 5.99e-4 | 0.262, 0.157, 0.206, 0.133, 0.242 | -1.600, -0.798, +0.004, +0.798, +1.597 | 20, 6, 10, 6, 18 | 0.958 |
| 0.01 | 2 | 5 | 0.00524 | 5 (5) | 0.00137 | 0.256, 0.164, 0.202, 0.149, 0.228 | -1.600, -0.794, +0.010, +0.793, +1.591 | 22, 5, 8, 10, 15 | 0.929 |

## D. The instance from broad starts (m = 60, the annealed Gaussian cohesion leak)

| γ | seed | clusters after R = 2 | fit error | clusters at the end (> 1% of the mass) | fit error | masses | centres | particles | |a|/a_max (max) |
|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 0 | 3 | 0.00102 | 3 (3) | 2.35e-7 | 0.307, 0.233, 0.460 | -1.500, +0.000, +1.500 | 17, 13, 30 | 0.814 |
| 0.01 | 1 | 3 | 0.00157 | 3 (3) | 9.68e-7 | 0.366, 0.278, 0.356 | -1.500, +0.000, +1.500 | 22, 16, 22 | 0.682 |
| 0.01 | 2 | 3 | 0.0015 | 3 (3) | 8.83e-7 | 0.423, 0.278, 0.300 | -1.500, +0.000, +1.500 | 26, 16, 18 | 0.667 |

## C. The finest resolution stage (the fourth bump, compact R = 1, σ₀/64, from the jittered lift of ϑ^b(q), 100000 slow units)

| γ | seed | clusters at the end (> 1% of the mass) | fit error | masses | centres | particles |
|---|---|---|---|---|---|---|
| 0 | 0 | 8 (8) | 0.00586 | 0.333, 0.033, 0.100, 0.200, 0.233, 0.033, 0.033, 0.033 | -1.499, -0.218, -0.185, +0.192, +1.147, +1.388, +1.717, +1.729 | 10, 1, 3, 6, 7, 1, 1, 1 |
| 0 | 1 | 10 (10) | 0.00718 | 0.333, 0.033, 0.033, 0.067, 0.200, 0.200, 0.033, 0.033, 0.033, 0.033 | -1.499, -0.205, -0.199, -0.187, +0.187, +1.114, +1.263, +1.534, +1.710, +1.721 | 10, 1, 1, 2, 6, 6, 1, 1, 1, 1 |
| 0 | 2 | 10 (10) | 0.0078 | 0.333, 0.033, 0.100, 0.033, 0.167, 0.200, 0.033, 0.033, 0.033, 0.033 | -1.499, -0.207, -0.189, +0.127, +0.193, +1.112, +1.338, +1.518, +1.691, +1.721 | 10, 1, 3, 1, 5, 6, 1, 1, 1, 1 |
| 0.01 | 0 | 8 (8) | 0.00548 | 0.327, 0.034, 0.102, 0.202, 0.232, 0.033, 0.036, 0.034 | -1.499, -0.215, -0.185, +0.196, +1.137, +1.378, +1.717, +1.729 | 10, 1, 3, 6, 7, 1, 1, 1 |
| 0.01 | 1 | 10 (10) | 0.00686 | 0.327, 0.035, 0.034, 0.068, 0.202, 0.199, 0.033, 0.033, 0.034, 0.035 | -1.499, -0.201, -0.200, -0.186, +0.190, +1.105, +1.259, +1.536, +1.711, +1.721 | 10, 1, 1, 2, 6, 6, 1, 1, 1, 1 |
| 0.01 | 2 | 10 (10) | 0.00748 | 0.328, 0.034, 0.102, 0.034, 0.169, 0.199, 0.033, 0.033, 0.033, 0.035 | -1.499, -0.206, -0.189, +0.129, +0.196, +1.103, +1.336, +1.520, +1.692, +1.721 | 10, 1, 3, 1, 5, 6, 1, 1, 1, 1 |

## F. The targets where Proposition D‴.32's count failed (Part K of `exp_cohesion_leak.py`; the annealed schedule, broad starts, m = 60, γ = 0.01)

| target | Σ|A_k|/a_max | seed | clusters at the end (> 1% of the mass) | fit error | masses | centres | needed |A_k|/a_max | |a|/a_max (max) |
|---|---|---|---|---|---|---|---|---|
| three heavy | 0.950 | 0 | 3 (3) | 1.04e-7 | 0.411, 0.309, 0.280 | -1.500, -0.000, +1.500 | 0.400, 0.300, 0.250 | 0.974 |
| three heavy | 0.950 | 1 | 3 (3) | 7.79e-8 | 0.421, 0.308, 0.271 | -1.500, +0.000, +1.500 | 0.400, 0.300, 0.250 | 0.973 |
| three heavy | 0.950 | 2 | 3 (3) | 1.63e-7 | 0.423, 0.313, 0.265 | -1.500, +0.000, +1.500 | 0.400, 0.300, 0.250 | 0.960 |
| two close | 0.400 | 0 | 3 (3) | 0.0683 | 0.059, 0.881, 0.059 | -1.354, +0.000, +1.354 | 0.200, 0.200 | 0.986 |
| two close | 0.400 | 1 | 2 (2) | 0.115 | 0.971, 0.029 | +0.035, +1.360 | 0.200, 0.200 | 0.946 |
| two close | 0.400 | 2 | 1 (1) | 0.13 | 1.000 | -0.000 | 0.200, 0.200 | 0.246 |

## E. The real loop and the linearisation

The real loop (k_e = 25, γ = 2, the positions at the instance's atoms, β(0) = (1.6, 0.7, 0.7), 200 time units): V = ½e² + (1/(γm))Σβ*log(β*/β) from 0.0418 to 2.14e-7, monotone; V̇ (from the vector field) + k_ee² within 1.11e-16 of max k_ee²; Σβ − m ≤ 2.40e-14; β at the end 1.00035, 0.99873, 1.00092.

- The averaged loop's linearisation in (Θ, log β) at the lift of ϑ^e(w), counts (6, 5, 4), γ = 0: the largest drift 2.45e-17; 15 eigenvalues below 10⁻⁹ in modulus; the others from -1.8336e-01 to -3.0798e-04.
- The averaged loop's linearisation in (Θ, log β) at the lift of ϑ^e(w), counts (6, 5, 4), γ = 0.01: the largest drift 2.45e-17; 15 eigenvalues below 10⁻⁹ in modulus; the others from -1.8337e-01 to -3.0812e-04.
- At γ = 0 the 45 nonzero eigenvalues are Proposition D‴.32's predicted spectrum (the separations at −σ₀w_k, the collective modes) to 2.08e-11.

## G. The stationary points under the Gaussian kernel, and Part A's infima

Class level (the classes of a run replicated, equal weights inside each, least squares on the stationarity conditions to rounding), at R = 1: the masses held at the run's, or freed by the replicator's condition u_k = λ.

| target | kernel | masses | stationary masses | fit error | |a|/a_max | largest residual |
|---|---|---|---|---|---|---|
| five (counts 6, 6, 20, 10, 18) | gauss | held | 0.1573, 0.1330, 0.2618, 0.2060, 0.2419 | 3.85e-4 | 0.9582, 0.9514, 0.8596, 0.8525, 0.8285 | 2.50e-17 |
| five (counts 6, 6, 20, 10, 18) | gauss | free (u_k = λ) | 0.1579, 0.1318, 0.2661, 0.2076, 0.2366 | 2.63e-4 | 0.9578, 0.9578, 0.8464, 0.8464, 0.8464 | 5.17e-17 |
| five (counts 6, 6, 20, 10, 18) | compact | held | 0.1573, 0.1330, 0.2618, 0.2060, 0.2419 | 2.46e-16 | 0.9538, 0.9396, 0.8594, 0.8495, 0.8268 | 3.47e-17 |
| five (counts 6, 6, 20, 10, 18) | compact | free (u_k = λ) | 0.1575, 0.1325, 0.2621, 0.2062, 0.2418 | 2.42e-16 | 0.9525, 0.9437, 0.8585, 0.8488, 0.8271 | 6.79e-17 |
| three equal (counts 28, 16, 16) | gauss | held | 0.4667, 0.2667, 0.2667 | 9.39e-5 | 0.4286, 0.7499, 0.7500 | 2.25e-17 |
| three equal (counts 28, 16, 16) | gauss | free (u_k = λ) | 0.5666, 0.2167, 0.2167 | 9.17e-5 | 0.3529, 0.9231, 0.9231 | 1.68e-10 |
| three equal (counts 28, 16, 16) | compact | held | 0.4667, 0.2667, 0.2667 | 7.32e-17 | 0.4286, 0.7500, 0.7500 | 4.25e-18 |
| three equal (counts 28, 16, 16) | compact | free (u_k = λ) | 0.4667, 0.2667, 0.2666 | 1.40e-16 | 0.4285, 0.7500, 0.7502 | 1.94e-17 |

- Part B's seed 1 (γ = 10⁻², the fit error 5.99e-4 at its end) continued: under the Gaussian kernel at R = 1 the fit error 5.99e-4 (s = 0), 4.45e-4 (s = 10000), 3.62e-4 (s = 20000), 3.20e-4 (s = 30000), 3.01e-4 (s = 40000), 2.91e-4 (s = 50000), 2.87e-4 (s = 60000), 2.84e-4 (s = 70000), 2.82e-4 (s = 80000), 2.81e-4 (s = 90000), 2.80e-4 (s = 100000); under the compact kernel at R = 1 from the same state 5.99e-4 (s = 0), 4.46e-4 (s = 10000), 3.59e-4 (s = 20000), 3.11e-4 (s = 30000), 2.85e-4 (s = 40000).
- Part A's target with the amplitudes scaled by 1 at fixed centres and widths: the least fit error with the masses held at (0.15, 0.30, 0.20, 0.25, 0.10) 0.387; over the ℓ¹ ball Σ_k|p_k| ≤ a_max (the masses free) 1.85e-8, at |p_k|/a_max 0.225, 0.150, 0.175, 0.125, 0.200.
- Part A's target with the amplitudes scaled by 1.25 at fixed centres and widths: the least fit error with the masses held at (0.15, 0.30, 0.20, 0.25, 0.10) 0.494; over the ℓ¹ ball Σ_k|p_k| ≤ a_max (the masses free) 0.0684, at |p_k|/a_max 0.270, 0.164, 0.194, 0.132, 0.239.

Numbers are produced by `experiments/exp_mass_transfer.py`; nothing is edited by hand.
