# Theorem J′: the stationary law of the dithered single-kernel loop (`exp_dither_stationary.py`)

Instance: n = 1, m = 1, f(x) = e^{−(x−0.5)²/0.5} + 0.4x, x_m(t) = 1.5[sin t′ + (1/3) sin 2t′] (t′ = ±t), a_max = 3, c_max = 2.5, s = 0.6 frozen, G = diag(2, 2), σ₀ = 0.4, ϑ⁰ = 0, τ₀ = 0.1; σ = σ₀/k_e, τ = τ₀/k_e. Grid 201×201 on [-3.5, 3.5]² (boundary density / peak density = 2.7e-06); Monte Carlo: 20000 replicas started from the first-order law π(1 + εg), Strang splitting with k_e·h = 1/2, burn-in ≥ 10 and averaging ≥ 40 slow-time units (whole numbers of reference periods), standard errors from the spread of the per-replica time averages (one standard error quoted).

## Averaged data

Gibbs expectations: π(alpha) = +0.47720, π(gamma) = +0.25015, π(alpha2) = +0.53568, π(gamma2) = +0.52916, π(alpha_gamma) = +0.29797; 2π(L̄) = 0.39279.

| orientation | rms b_fr | rms b_osc | π(‖b‖²_{G⁻¹}) | non-gradient fraction π(‖b_⊥‖²)/π(‖b‖²) | EPR coefficient τ₀⁻¹π(‖b_⊥‖²_{G⁻¹}) | rms curl(G⁻¹b) / rms sym D(G⁻¹b) | circulation coefficient: ∫Θ∧j₁ + π⟨u∧u̇⟩ |
|---|---|---|---|---|---|---|---|
| +1 | 0.1431 | 0.2247 | 0.0454 | 0.316 | 0.1436 | 0.759 | +0.01165 + -0.01672 = -0.00507 |
| -1 | 0.1431 | 0.2247 | 0.0256 | 0.462 | 0.1181 | 0.925 | -0.01027 + +0.01672 = +0.00645 |

The first-order coefficients ∫φ g dπ of the stationary means (m_ε(φ) = π(φ) + ε∫φg dπ + O(ε²)):

| φ | orientation +1 | orientation −1 | even part | odd part |
|---|---|---|---|---|
| alpha | -0.25517 | +0.16734 | -0.04392 | -0.21126 |
| gamma | -0.42622 | +0.26933 | -0.07845 | -0.34777 |
| alpha2 | -0.12537 | +0.14445 | +0.00954 | -0.13491 |
| gamma2 | -0.02952 | +0.03220 | +0.00134 | -0.03086 |
| alpha_gamma | -0.06144 | +0.11743 | +0.02799 | -0.08944 |

## Table 1 — stationary means of the real loop (orientation +1) against the Gibbs value and the first-order prediction

### φ = alpha

| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |
|---|---|---|---|---|---|---|---|
| 10 | +0.44963 ± 0.00132 | +0.47720 | +0.45169 | -0.02757 | -0.00205 | -0.2757 | -0.2054 |
| 20 | +0.45933 ± 0.00130 | +0.47720 | +0.46445 | -0.01788 | -0.00512 | -0.3575 | -2.0470 |
| 40 | +0.46828 ± 0.00129 | +0.47720 | +0.47082 | -0.00892 | -0.00254 | -0.3568 | -4.0651 |
| 80 | +0.47256 ± 0.00129 | +0.47720 | +0.47401 | -0.00464 | -0.00145 | -0.3714 | -9.2952 |

Slope of log|m_ε − π| against log k_e: -0.87 (first order: −1); slope of log|m_ε − first order|: not resolved beyond the Monte Carlo error (second order: −2).

### φ = gamma

| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |
|---|---|---|---|---|---|---|---|
| 10 | +0.20605 ± 0.00184 | +0.25015 | +0.20753 | -0.04411 | -0.00148 | -0.4411 | -0.1484 |
| 20 | +0.22325 ± 0.00184 | +0.25015 | +0.22884 | -0.02690 | -0.00559 | -0.5381 | -2.2376 |
| 40 | +0.23747 ± 0.00182 | +0.25015 | +0.23950 | -0.01269 | -0.00203 | -0.5074 | -3.2478 |
| 80 | +0.24183 ± 0.00181 | +0.25015 | +0.24483 | -0.00832 | -0.00300 | -0.6659 | -19.1761 |

Slope of log|m_ε − π| against log k_e: -0.83 (first order: −1); slope of log|m_ε − first order|: not resolved beyond the Monte Carlo error (second order: −2).

### φ = alpha2

| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |
|---|---|---|---|---|---|---|---|
| 10 | +0.52765 ± 0.00116 | +0.53568 | +0.52314 | -0.00803 | +0.00451 | -0.0803 | +0.4508 |
| 20 | +0.52663 ± 0.00117 | +0.53568 | +0.52941 | -0.00904 | -0.00278 | -0.1809 | -1.1106 |
| 40 | +0.53101 ± 0.00118 | +0.53568 | +0.53254 | -0.00466 | -0.00153 | -0.1865 | -2.4462 |
| 80 | +0.53380 ± 0.00117 | +0.53568 | +0.53411 | -0.00188 | -0.00031 | -0.1501 | -1.9800 |

Slope of log|m_ε − π| against log k_e: -0.72 (first order: −1); slope of log|m_ε − first order|: -0.70 over the k_e where it is resolved (2 of 4) (second order: −2).

### φ = gamma2

| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |
|---|---|---|---|---|---|---|---|
| 10 | +0.52890 ± 0.00118 | +0.52916 | +0.52621 | -0.00026 | +0.00269 | -0.0026 | +0.2691 |
| 20 | +0.52815 ± 0.00120 | +0.52916 | +0.52769 | -0.00101 | +0.00046 | -0.0202 | +0.1855 |
| 40 | +0.53012 ± 0.00118 | +0.52916 | +0.52843 | +0.00096 | +0.00170 | +0.0383 | +2.7142 |
| 80 | +0.52894 ± 0.00117 | +0.52916 | +0.52880 | -0.00022 | +0.00015 | -0.0177 | +0.9414 |

Slope of log|m_ε − π| against log k_e: -0.08 (first order: −1); slope of log|m_ε − first order|: not resolved beyond the Monte Carlo error (second order: −2).

### φ = alpha_gamma

| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |
|---|---|---|---|---|---|---|---|
| 10 | +0.28996 ± 0.00080 | +0.29797 | +0.29183 | -0.00801 | -0.00186 | -0.0801 | -0.1862 |
| 20 | +0.29253 ± 0.00081 | +0.29797 | +0.29490 | -0.00544 | -0.00237 | -0.1089 | -0.9490 |
| 40 | +0.29628 ± 0.00081 | +0.29797 | +0.29643 | -0.00169 | -0.00016 | -0.0677 | -0.2505 |
| 80 | +0.29702 ± 0.00080 | +0.29797 | +0.29720 | -0.00096 | -0.00019 | -0.0764 | -1.1979 |

Slope of log|m_ε − π| against log k_e: -1.09 (first order: −1); slope of log|m_ε − first order|: 0.35 over the k_e where it is resolved (2 of 4) (second order: −2).

## Table 2 — orientation-odd part of the stationary means, (m₊ − m₋)/2, against its first-order prediction (p₁₊ − p₁₋)/(2k_e)

| k_e | φ | measured odd part ± s.e. | predicted | measured/predicted | ratio·(1 ± s.e./measured) |
|---|---|---|---|---|---|
| 10 | alpha | -0.01946 ± 0.00091 | -0.02113 | 0.921 | ± 0.043 |
| 10 | gamma | -0.03242 ± 0.00127 | -0.03478 | 0.932 | ± 0.037 |
| 10 | alpha2 | -0.01107 ± 0.00083 | -0.01349 | 0.821 | ± 0.062 |
| 10 | gamma2 | -0.00453 ± 0.00083 | -0.00309 | 1.466 | ± 0.271 |
| 10 | alpha_gamma | -0.00799 ± 0.00058 | -0.00894 | 0.893 | ± 0.064 |
| 20 | alpha | -0.01295 ± 0.00091 | -0.01056 | 1.226 | ± 0.086 |
| 20 | gamma | -0.01919 ± 0.00128 | -0.01739 | 1.103 | ± 0.074 |
| 20 | alpha2 | -0.00929 ± 0.00083 | -0.00675 | 1.377 | ± 0.124 |
| 20 | gamma2 | -0.00157 ± 0.00083 | -0.00154 | 1.018 | ± 0.541 |
| 20 | alpha_gamma | -0.00571 ± 0.00058 | -0.00447 | 1.276 | ± 0.129 |
| 40 | alpha | -0.00675 ± 0.00090 | -0.00528 | 1.278 | ± 0.171 |
| 40 | gamma | -0.00999 ± 0.00127 | -0.00869 | 1.149 | ± 0.147 |
| 40 | alpha2 | -0.00410 ± 0.00083 | -0.00337 | 1.216 | ± 0.246 |
| 40 | gamma2 | -0.00012 ± 0.00083 | -0.00077 | 0.160 | ± 1.077 |
| 40 | alpha_gamma | -0.00230 ± 0.00057 | -0.00224 | 1.030 | ± 0.256 |
| 80 | alpha | -0.00329 ± 0.00090 | -0.00264 | 1.246 | ± 0.341 |
| 80 | gamma | -0.00600 ± 0.00128 | -0.00435 | 1.380 | ± 0.293 |
| 80 | alpha2 | -0.00188 ± 0.00083 | -0.00169 | 1.116 | ± 0.492 |
| 80 | gamma2 | -0.00024 ± 0.00082 | -0.00039 | 0.615 | ± 2.138 |
| 80 | alpha_gamma | -0.00125 ± 0.00057 | -0.00112 | 1.114 | ± 0.511 |

## Table 3 — stationary tracking error: k_e²·E‖e‖² against 2π(L̄) (Theorem J′(d): equal to O(1/k_e))

| k_e | k_e²E e² (orientation +1) ± s.e. | (−1) | 2π(L̄) | ratio (+1) | ratio (−1) |
|---|---|---|---|---|---|
| 10 | 0.41416 ± 0.00039 | 0.40175 ± 0.00038 | 0.39279 | 1.0544 | 1.0228 |
| 20 | 0.40586 ± 0.00038 | 0.39842 ± 0.00038 | 0.39279 | 1.0333 | 1.0143 |
| 40 | 0.39959 ± 0.00038 | 0.39577 ± 0.00037 | 0.39279 | 1.0173 | 1.0076 |
| 80 | 0.39641 ± 0.00038 | 0.39450 ± 0.00037 | 0.39279 | 1.0092 | 1.0044 |

## Table 4 — the circulation observable Λ_ε = E[Θ ∧ k_e·drift] against ε[∫Θ∧j₁ + π⟨u∧u̇⟩]

| k_e | orientation | measured Λ_ε ± s.e. | predicted ε·coefficient | (measured − predicted)/s.e. |
|---|---|---|---|---|
| 10 | +1 | +0.00365 ± 0.00051 | -0.00051 | +8.15 |
| 10 | -1 | -0.00057 ± 0.00053 | +0.00065 | -2.29 |
| 20 | +1 | +0.00056 ± 0.00051 | -0.00025 | +1.59 |
| 20 | -1 | +0.00119 ± 0.00052 | +0.00032 | +1.67 |
| 40 | +1 | -0.00041 ± 0.00051 | -0.00013 | -0.55 |
| 40 | -1 | +0.00007 ± 0.00052 | +0.00016 | -0.18 |
| 80 | +1 | +0.00031 ± 0.00051 | -0.00006 | +0.73 |
| 80 | -1 | -0.00034 ± 0.00052 | +0.00008 | -0.81 |

The predicted circulation is below the Monte Carlo resolution of these runs wherever |predicted| < 2 s.e.; the table is a consistency check, not a measurement of the current. The current machinery (Hodge decomposition of b, j₁ = π(b − τ₀G∇g), Proposition J′.2) is checked deterministically on the first-order averaged SDE dΘ = (−G∇J + εb)ds + √(2τ₀)G^{1/2}dW, whose stationary density m_ε and current j_ε are computed on the grid (orientation +1; the discrete current at ε = 0, ∫Θ∧j₀ = -6.56e-05, a residual of the central scheme, is subtracted): Proposition J′.2 predicts ∫Θ∧j_ε = ε∫Θ∧j₁ + O(ε²) with ∫Θ∧j₁ = +0.01165, and m_ε = π(1 + εg) + O(ε²).

| ε | ∫Θ∧(j_ε − j₀)/ε | ratio to ∫Θ∧j₁ | ‖m_ε − m₀(1 + εg)‖₁/ε² |
|---|---|---|---|
| 0.2 | +0.01716 | 1.4722 | 0.307 |
| 0.1 | +0.01406 | 1.2067 | 0.308 |
| 0.05 | +0.01266 | 1.0861 | 0.307 |
| 0.025 | +0.01198 | 1.0283 | 0.305 |

Numbers are produced by `experiments/exp_dither_stationary.py`; nothing is edited by hand.
