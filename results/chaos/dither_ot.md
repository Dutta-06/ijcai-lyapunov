# Remark D‴.16: dither and the mass constraint — the dithered averaged slow loop under the nearest-atom and the OT leak (`exp_dither_ot.py`)

The three-bump instance with the matching prior (three atoms at the exact representation, masses 1/3), σ₀ = 0.15σ_crit = 0.00707, m = 15 particles started at the collapsed configuration with counts (6, 5, 4) (jitter 0.02), the Langevin dynamics dϑ_i = −G∇J_m ds + √(2τ₀)G^{1/2}dW_i in slow time (Euler–Maruyama, ds = 0.05; the leak target recomputed every 10 steps), cluster masses by nearest atom recorded every 10 slow units; horizons 50000, 200000, 50000 slow units at τ₀ = 0.005, 0.01, 0.02; seeds 0, 1. The Gibbs prediction for the nearest-atom leak is the Laplace approximation over the count vectors, P(n) ∝ multinomial(n)·exp(−m𝒢*_{n/m}/τ₀)·Π_k det(σ₀I − R^G_k)^{−(n_k−1)/2}·det(H_n)^{−1/2} (count vectors whose replicated configuration is a saddle carry no Laplace mode and are omitted); 'hop records' counts the records at which some particle's nearest atom changed.

## τ₀ = 0.005 (horizon 50000 slow units; the leak's part of the lowest ridge σ₀Δ_ν²/8 = 0.035 = 6.9τ₀)

The energy–multinomial law over the count vectors, multinomial(n)·exp(−m𝒢*_n/τ₀) (count vectors with n/m ∈ W_adm and the saddles alike; not an order of an expansion: the Laplace determinant factors are of the same order and are omitted), P(n): (5, 5, 5): 0.862, (5, 4, 6): 0.075, (5, 6, 4): 0.043, (6, 4, 5): 0.013, (6, 5, 4): 0.004; with the Laplace determinant factors (saddle count vectors omitted): (5, 5, 5): 0.720, (5, 6, 4): 0.243, (5, 4, 6): 0.035.

| leak | seed | records with a label change / records (undone at the next change) | time fraction by counts (top four) | mean |w − q|₁ | final counts | excursions off the balanced counts after the first balanced record: returns / exchanges; off-balance fraction |
|---|---|---|---|---|---|---|
| nearest | 0 | 0/4999 (0) | (6, 5, 4): 1.000 | 0.133 | (6, 5, 4) | — |
| nearest | 1 | 1/4999 (0) | (6, 5, 4): 0.976, (5, 6, 4): 0.024 | 0.133 | (5, 6, 4) | — |
| ot | 0 | 1/4999 (0) | (5, 5, 5): 0.996, (6, 5, 4): 0.004 | 0.001 | (5, 5, 5) | 0 / 0; 0.00e+00 |
| ot | 1 | 1/4999 (0) | (5, 5, 5): 0.997, (6, 5, 4): 0.003 | 0.000 | (5, 5, 5) | 0 / 0; 0.00e+00 |

## τ₀ = 0.01 (horizon 200000 slow units; the leak's part of the lowest ridge σ₀Δ_ν²/8 = 0.035 = 3.5τ₀)

The energy–multinomial law over the count vectors, multinomial(n)·exp(−m𝒢*_n/τ₀) (count vectors with n/m ∈ W_adm and the saddles alike; not an order of an expansion: the Laplace determinant factors are of the same order and are omitted), P(n): (5, 5, 5): 0.576, (5, 4, 6): 0.156, (5, 6, 4): 0.118, (6, 4, 5): 0.064, (6, 5, 4): 0.037; with the Laplace determinant factors (saddle count vectors omitted): (5, 6, 4): 0.527, (5, 5, 5): 0.383, (5, 4, 6): 0.058.

| leak | seed | records with a label change / records (undone at the next change) | time fraction by counts (top four) | mean |w − q|₁ | final counts | excursions off the balanced counts after the first balanced record: returns / exchanges; off-balance fraction |
|---|---|---|---|---|---|---|
| nearest | 0 | 111/19999 (55) | (5, 5, 5): 0.562, (6, 4, 5): 0.208, (5, 4, 6): 0.097, (6, 5, 4): 0.069 | 0.059 | (5, 4, 6) | — |
| nearest | 1 | 93/19999 (41) | (5, 5, 5): 0.511, (6, 4, 5): 0.203, (6, 5, 4): 0.127, (5, 4, 6): 0.088 | 0.065 | (5, 5, 5) | — |
| ot | 0 | 131/19999 (85) | (5, 5, 5): 0.991, (5, 4, 6): 0.004, (5, 6, 4): 0.002, (6, 4, 5): 0.001 | 0.001 | (5, 5, 5) | 65 / 0; 7.51e-03 |
| ot | 1 | 115/19999 (78) | (5, 5, 5): 0.994, (5, 4, 6): 0.003, (5, 6, 4): 0.001, (6, 5, 4): 0.001 | 0.001 | (5, 5, 5) | 56 / 0; 5.16e-03 |

## τ₀ = 0.02 (horizon 50000 slow units; the leak's part of the lowest ridge σ₀Δ_ν²/8 = 0.035 = 1.7τ₀)

The energy–multinomial law over the count vectors, multinomial(n)·exp(−m𝒢*_n/τ₀) (count vectors with n/m ∈ W_adm and the saddles alike; not an order of an expansion: the Laplace determinant factors are of the same order and are omitted), P(n): (5, 5, 5): 0.332, (5, 4, 6): 0.157, (5, 6, 4): 0.137, (6, 4, 5): 0.101, (6, 5, 4): 0.076; with the Laplace determinant factors (saddle count vectors omitted): (5, 6, 4): 0.520, (5, 5, 5): 0.187, (5, 7, 3): 0.098.

| leak | seed | records with a label change / records (undone at the next change) | time fraction by counts (top four) | mean |w − q|₁ | final counts | excursions off the balanced counts after the first balanced record: returns / exchanges; off-balance fraction |
|---|---|---|---|---|---|---|
| nearest | 0 | 281/4999 (113) | (5, 5, 5): 0.211, (6, 4, 5): 0.163, (5, 4, 6): 0.160, (5, 6, 4): 0.112 | 0.125 | (5, 4, 6) | — |
| nearest | 1 | 275/4999 (146) | (5, 5, 5): 0.287, (6, 4, 5): 0.165, (6, 5, 4): 0.154, (5, 4, 6): 0.152 | 0.105 | (5, 5, 5) | — |
| ot | 0 | 364/4999 (213) | (5, 5, 5): 0.861, (5, 4, 6): 0.057, (5, 6, 4): 0.026, (6, 5, 4): 0.026 | 0.019 | (5, 5, 5) | 150 / 20; 1.37e-01 |
| ot | 1 | 344/4999 (210) | (5, 5, 5): 0.894, (5, 4, 6): 0.037, (6, 4, 5): 0.025, (5, 6, 4): 0.020 | 0.014 | (5, 5, 5) | 148 / 18; 9.77e-02 |

## Reading

Nearest-atom leak: record-level label changes (every 10 slow units, about half undone at the next change — recrossings) at 1.00e-05 per slow unit at τ₀ = 0.005, 5.10e-04 per slow unit at τ₀ = 0.01, 5.56e-03 per slow unit at τ₀ = 0.02 — Arrhenius barriers ΔV = 0.048 from τ₀ = 0.01, 0.02, 0.039 from τ₀ = 0.005, 0.01 (the second pair rests on one event at τ₀ = 0.005); the leak's part of the lowest ridge is σ₀Δ_ν²/8 = 0.035. Time fraction at the balanced counts (5, 5, 5), nearest-atom leak against the OT leak against the energy–multinomial law: τ₀ = 0.005: 0.000 / 0.996 / 0.862; τ₀ = 0.01: 0.537 / 0.992 / 0.576; τ₀ = 0.02: 0.249 / 0.878 / 0.332 (two seeds; the balanced indicator's autocorrelation is long at τ₀ = 0.01, so these fractions carry standard errors of order 0.1 there). OT leak: the plan's masses are q by construction; the counts by nearest atom leave the balanced counts in excursions that end with the crossing particle back in its cluster (returns) or, rarely, with two particles' assignments exchanged — τ₀ = 0.005: 0 returns, 0 exchanges, off-balance fraction 0.00e+00; τ₀ = 0.01: 121 returns, 0 exchanges, off-balance fraction 6.33e-03; τ₀ = 0.02: 298 returns, 38 exchanges, off-balance fraction 1.17e-01 — an effective barrier of 0.058 between τ₀ = 0.01 and 0.02, and label changes as frequent as the nearest-atom leak's (ratios 1.18–1.30): the crossing barrier is the nearest-atom ridge barrier, and what the mass constraint changes is that a crossing is reversible — the particle stays assigned to its own atom and is pulled back. Dither does on the Arrhenius time scale, and only in distribution, what the mass constraint does deterministically in ln 2/σ₀.

Numbers are produced by `experiments/exp_dither_ot.py`; nothing is edited by hand.
