# Theorem C′: the barycentric OT leak (`exp_barycentric_leak.py`)

## 1. Pure leak toward a 48-atom prior (m = 6, n = 1, lifted dimension 3, G = diag(2, 1, 0.5), σ_m = 1, dt = 0.01, 3000 steps)

| seed | W²(μ₀, ν) | W²(μ_T, ν) | centroid defect V_bar² at T | within-cell variance Var_γ at T | W² non-increasing | plan switches | best of 20 equal-mass Lloyd runs |
|---|---|---|---|---|---|---|---|
| 0 | 2.5526 | 1.1642 | 1.1e-26 | 1.1642 | True | 3 | 1.1708 |
| 1 | 2.3032 | 1.2128 | 4.8e-26 | 1.2128 | True | 2 | 1.1708 |
| 2 | 2.5807 | 1.2030 | 2.2e-25 | 1.2030 | True | 6 | 1.1708 |

The flow stops exactly where W² = Var_γ (every particle at the barycentre of its cell): an equal-mass centroidal configuration of the prior. Like Lloyd's algorithm it finds a local minimiser of the quantisation energy, not always the global one: the best of 20 Lloyd runs lies between the end points of the three starts (last two columns).
With m₀ = m the same leak gives the exact decay of Theorem C(i): max relative deviation of W(μ_t, μ⁰) from (1 − σ_m dt)^k W₀ over 500 Euler steps: 5.1e-16.

## 2. Closed loop: three leaks of equal gain (σ_m = 0.5) toward three priors; representable plant, m = 4, k_e = 4

| seed | b = W(μ*, cloud) | b_bar (centroid defect of the ideal w.r.t. the cloud) | W(μ*, subset) | labelled bias to the subset | A: RMS‖e‖ [20, 40], W(μ_T, μ*) | B | C | exchange bound on Δ_gap at μ* (cloud) | η(0) = W(μ₀, μ*) | 4η(0)(b + η(0)) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.446 | 0.168 | 0.647 | 1.987 | 0.0000, 0.000 | 0.0035, 0.167 | 0.0431, 0.614 | 0.0251 | 0.848 | 4.385 |
| 1 | 0.328 | 0.097 | 1.027 | 1.479 | 0.0000, 0.000 | 0.0024, 0.097 | 0.0430, 1.016 | 0.0155 | 0.627 | 2.394 |
| 2 | 0.401 | 0.159 | 0.395 | 1.765 | 0.0000, 0.000 | 0.0069, 0.157 | 0.0293, 0.352 | 0.0256 | 0.553 | 2.112 |

Reading: each loop settles at the configuration its leak anchors; for (B) the distance W(μ_T, μ*) equals b_bar, the distance of the ideal from the centroids of the prior's cells, not b = W(μ*, ν) (Theorem C′(i), (v)). The gap condition 4η(b + η) < Δ_gap of Theorem C′(ii)–(iv) fails at t = 0 in every run: 4η(0)(b + η(0)) exceeds the exchange bound on Δ_gap (itself an upper bound on Δ_gap) by the factors 175, 154, 83, so the local estimate (iv) is not what these runs test. The comparison with (C) shows the cost of a prior of the right size at the wrong positions; a labelled leak toward the cloud's own equal-mass 4-quantiser would carry the same bias b_bar (Theorem C′(v)).

Numbers are produced by `experiments/exp_barycentric_leak.py`; nothing is edited by hand.
