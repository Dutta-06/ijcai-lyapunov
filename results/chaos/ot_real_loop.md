# Corollary D‴.14: the real closed loop with the OT leak — the orbit is the lift of the weighted K₀-atom loop's, the same for every commensurate size (`exp_ot_real_loop.py`)

The real closed loop of Theorem J on the three-bump instance (k_e = 25, σ = σ₀/k_e with σ₀ = 0.15σ_crit = 0.00707, RK4 plant with the control and the reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = 0.001), 100 time units = 15.9 periods = 4.0 slow units; statistics over the last 5 periods on the 400-point grid of the reference orbit. Starts: the clouds at the end of `exp_ot_leak.py`'s averaged runs (contracted where the realised weights are admissible) (commensurate set, m = 15, 60, M = 510; the OT runs' clouds carry the prior's masses q = (1/3, 1/3, 1/3), the nearest-atom runs' clouds their conserved sampled weights), and the collapsed lifts: the three atoms of the exact representation for the OT leak (the weighted 3-atom loop with weights q, itself the m = 3 loop), and the replicated weighted stationary configurations Θ*(w^m) for the nearest-atom leak. Seeds 0, 1; the leak target (the nearest atom, or the barycentric projection of the vertex plan) recomputed every 100 steps.

## 1. The OT leak: every size on the same orbit

| seed | system | cluster weights at the start | within-class spread at the end (rms, whitened) | ‖f̂_m − f̂_M‖²_{L²(orbit)} (M = 510) | ‖f̂_m − f̂_lift‖² | fit error along the orbit | rms tracking error | W²_{G⁻¹}(μ_m, μ_M) at the end | W² to the lift |
|---|---|---|---|---|---|---|---|---|---|
| 0 | m15 | (0.3333, 0.3333, 0.3333) | 3.1e-02 | 2.739e-08 | 3.552e-08 | 2.907e-04 | 7.296e-06 | 1.391e-03 | 9.885e-04 |
| 0 | m60 | (0.3333, 0.3333, 0.3333) | 1.3e-02 | 7.459e-09 | 5.199e-10 | 3.517e-05 | 8.888e-07 | 3.344e-04 | 1.646e-04 |
| 0 | m510 | (0.3333, 0.3333, 0.3333) | 2.0e-02 | 0.000e+00 | 1.168e-08 | 1.666e-04 | 4.205e-06 | 1.533e-15 | 4.157e-04 |
| 0 | lift (m = 3) | (0.3333, 0.3333, 0.3333) | 0.0e+00 | 1.168e-08 | 0.000e+00 | 2.361e-08 | 7.936e-12 | 4.157e-04 | 4.737e-15 |
| 1 | m15 | (0.3333, 0.3333, 0.3333) | 3.3e-02 | 8.582e-09 | 8.595e-09 | 1.430e-04 | 3.486e-06 | 1.188e-03 | 1.088e-03 |
| 1 | m60 | (0.3333, 0.3333, 0.3333) | 2.8e-02 | 8.529e-09 | 4.941e-09 | 1.084e-04 | 2.668e-06 | 8.627e-04 | 7.776e-04 |
| 1 | m510 | (0.3333, 0.3333, 0.3333) | 1.4e-02 | 0.000e+00 | 4.252e-09 | 1.006e-04 | 2.537e-06 | 1.421e-15 | 1.988e-04 |
| 1 | lift (m = 3) | (0.3333, 0.3333, 0.3333) | 0.0e+00 | 4.252e-09 | 0.000e+00 | 2.361e-08 | 7.936e-12 | 1.988e-04 | 4.737e-15 |

## 2. The nearest-atom leak: every size with admissible weights on the orbit of its own weights

| seed | m | conserved weights w^m | splitting margin σ₀ − max κ_R at w^m (Proposition D‴.5(i); ≤ 0: the lift is a saddle orbit) | within-class spread at the end (rms, whitened) | ‖f̂_m − f̂_M‖² | ‖f̂_m − f̂_{lift(w^m)}‖² | ‖f̂_{lift(w^m)} − f̂_{lift(w^M)}‖² | fit error along the orbit (m-system / its lift) | rms tracking error (m-system / its lift) | W²(μ_m, μ_M) | W² to its lift |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 15 | (0.4000, 0.3333, 0.2667) | -0.00348 (saddle) | 1.1e-01 | 1.823e-03 | 3.311e-06 | 1.898e-03 | 4.435e-02 / 4.557e-02 | 1.141e-03 / 1.172e-03 | 4.883e+00 | 1.251e-02 |
| 0 | 60 | (0.3167, 0.3500, 0.3333) | +0.00561 | 5.5e-04 | 1.062e-04 | 1.468e-11 | 1.062e-04 | 1.334e-02 / 1.334e-02 | 3.428e-04 / 3.427e-04 | 1.568e+00 | 3.068e-07 |
| 0 | 510 | (0.3098, 0.3157, 0.3745) | +0.00440 | 1.3e-03 | 0.000e+00 | 2.521e-11 | 0.000e+00 | 2.267e-02 / 2.266e-02 | 5.833e-04 / 5.832e-04 | 1.588e-15 | 1.814e-06 |
| 1 | 15 | (0.4000, 0.2000, 0.4000) | -0.00108 (saddle) | 4.2e-02 | 2.468e-03 | 4.136e-08 | 2.474e-03 | 7.471e-02 / 7.477e-02 | 1.910e-03 / 1.911e-03 | 5.547e+00 | 1.798e-03 |
| 1 | 60 | (0.4000, 0.2500, 0.3500) | -0.00216 (saddle) | 1.0e-01 | 1.075e-03 | 1.831e-06 | 1.119e-03 | 4.455e-02 / 4.541e-02 | 1.141e-03 / 1.163e-03 | 3.735e+00 | 1.072e-02 |
| 1 | 510 | (0.3196, 0.3255, 0.3549) | +0.00544 | 7.2e-04 | 0.000e+00 | 9.607e-12 | 0.000e+00 | 1.245e-02 / 1.245e-02 | 3.203e-04 / 3.203e-04 | 1.588e-15 | 5.172e-07 |

Reading. Under the OT leak the learned functions of m = 15, 60 and M = 510 agree along the orbit to ‖f̂_m − f̂_M‖² = 7.5e-09–2.7e-08 and with the 3-atom lift to 5.2e-10–3.6e-08 (the within-class spread of the clouds, second order); under the nearest-atom leak they disagree by 1.1e-04–2.5e-03. Where both realised weight vectors are admissible this disagreement is the distance between the lifted orbits of the conserved weights to a ratio of 1.0003 (1 pair), and each admissible m-system is within 2.5e-11 of its own lift (within-class spread 5.5e-04–1.3e-03); the 3 systems with inadmissible realised weights — seed 0, m = 15 (margin -0.00348), seed 1, m = 15 (margin -0.00108), seed 1, m = 60 (margin -0.00216) — are not on their lifts, which are saddle orbits: their over-weighted class is split (within-class spread 4.2e-02–1.1e-01), the distance to the lift 4.1e-08–3.3e-06 being the split's second-order effect on f̂, and the lift-to-lift ratio at those pairs is 0.960, 0.998, 0.961. The rms tracking error is 8.9e-07–7.3e-06 under the OT leak (the m-systems' remaining fit error over k_e; the collapsed lift tracks to 7.9e-12, the scheme's floor) against 3.2e-04–1.9e-03 under the nearest-atom leak — the weight bias of the fit is a tracking error, which the mass constraint removes for every size at once.

Numbers are produced by `experiments/exp_ot_real_loop.py`; nothing is edited by hand.
