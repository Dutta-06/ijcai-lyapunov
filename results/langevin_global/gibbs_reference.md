# Reference values for Theorem H and Conjecture H.6 (n = 1 experiment)

Computed by `experiments/exp_gibbs_reference.py` (concave dual on a grid; Anderson-mixed fixed point for the exactly averaged loop).  G = diag(4, 1, 1), ϑ⁰ = 0, f = sin 2x + 0.3x, x_m = 1.5 sin t + 0.5 sin 2t.

| quantity | (σ, τ) | k_e | grid (Nγ, Nκ, Nα, N_t) | value | residual |
|---|---|---|---|---|---|
| L̄(μ*_τ) | (0.02, 0.004) | 10 | (241, 241, 1201, 256) | 0.11118 | 2.4e-10 |
| L̄(μ*_τ) | (0.02, 0.004) | 10 | (301, 301, 1601, 320) | 0.11118 | 2.1e-08 |
| L̄(μ*_τ) | (0.002, 0.0005) | 10 | (241, 241, 1201, 256) | 0.00850 | 8.4e-09 |
| L̄(μ*_τ) | (0.002, 0.0005) | 10 | (301, 301, 1601, 320) | 0.00850 | 6.4e-09 |
| L̄(μ*_τ), scaling σ = 0.2/k_e, τ = 0.04/k_e | (0.02, 0.004) | 10 | (161, 161, 801, 256) | 0.11118 | — |
| L̄ of the exactly averaged stationary state, same scaling | (0.02, 0.004) | 10 | (161, 161, 801, 256) | 0.11415 | fixed point 5.7e-11 after 12 iterations |
| ‖f̂(exactly averaged) − f̂(μ*_τ)‖ in L²(ρ), and k_e × it | (0.02, 0.004) | 10 | (161, 161, 801, 256) | 0.00681 | k_e·d = 0.0681 |
| L̄(μ*_τ), scaling σ = 0.2/k_e, τ = 0.04/k_e | (0.01, 0.002) | 20 | (161, 161, 801, 256) | 0.11118 | — |
| L̄ of the exactly averaged stationary state, same scaling | (0.01, 0.002) | 20 | (161, 161, 801, 256) | 0.11148 | fixed point 6.1e-11 after 10 iterations |
| ‖f̂(exactly averaged) − f̂(μ*_τ)‖ in L²(ρ), and k_e × it | (0.01, 0.002) | 20 | (161, 161, 801, 256) | 0.00143 | k_e·d = 0.0287 |
| L̄(μ*_τ), scaling σ = 0.2/k_e, τ = 0.04/k_e | (0.005, 0.001) | 40 | (161, 161, 801, 256) | 0.11118 | — |
| L̄ of the exactly averaged stationary state, same scaling | (0.005, 0.001) | 40 | (161, 161, 801, 256) | 0.11101 | fixed point 3.7e-12 after 11 iterations |
| ‖f̂(exactly averaged) − f̂(μ*_τ)‖ in L²(ρ), and k_e × it | (0.005, 0.001) | 40 | (161, 161, 801, 256) | 0.00097 | k_e·d = 0.0388 |

In the scaling of Conjecture H.6 (σ = σ₀/k_e, τ = τ₀/k_e) μ*_τ does not depend on k_e, so the first column of each pair should be constant.  The offset of the exactly averaged stationary state is measured by the L²(ρ) distance between the two learned functions (the scalar loss difference can change sign); it is the O(1/k_e) offset that makes the fixed-k_e statement false.
