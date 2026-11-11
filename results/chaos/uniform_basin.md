# Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate (`exp_uniform_basin.py`)

The real closed loop of Theorem J on the three-bump instance: k_e = 25, σ = σ₀/k_e with σ₀ = 0.15σ_crit = 0.00707, RK4 plant with the control and the reference at the stages and the parameters frozen over the step, Euler parameters, dt = 0.005, horizon 8000 time units (1273 periods, 320 slow units); the leak target refreshed every 200 steps; records every 5 time units. Fitted rates: least-squares slope of log d_∞ over the second half of the run; they are compared with the exponents of Part F, which are the test of the rate.

## A. The matching prior at the matching masses (the lift is the exact representation, a constant orbit)

Prediction: λ_sep = 0.00707 (= σ₀ at zero residual), λ_coll = 0.00738, rate min{λ_coll, λ_sep}/k_e = 2.830e-04 per time unit for m ≥ 15; for m = 3 (one particle per class, no separation modes) the slowest mode is the collective one, λ_coll/k_e = 2.952e-04, followed by λ_j(H_q)/k_e = 3.616e-04, 4.076e-04. Every particle is displaced from its class atom by ρ in a random direction of the G⁻¹-unit sphere (the same displacements under both leaks); the Voronoi cells extend at least Δ_ν/2 = 3.13 from the atoms. The rate and C are printed only for runs that end at the counts q with constant counts over the second half.

| leak | ρ | m | final counts by nearest atom | d_∞(0) → d_∞(end) | fitted rate of d_∞ (second half) | C(T) = max_{t≤T} d_∞(t)e^{μt}/d_∞(0) at T = 2000, 4000, 6000, 8000 (μ the prediction for that m) | W²(μ_Θ, ν): start → end | rms tracking error (last 10 periods) |
|---|---|---|---|---|---|---|---|---|
| nearest | 0.3 | 3 | (1, 1, 1) | 0.300 → 1.04e-02 | 3.484e-04 | 1.000, 1.000, 1.000, 1.000 | 9.000e-02 → 8.69e-05 | 1.4e-05 |
| nearest | 0.3 | 15 | (5, 5, 5) | 0.300 → 3.63e-02 | 2.694e-04 | 1.113, 1.114, 1.122, 1.162 (maximum at the last record) | 9.000e-02 → 8.29e-04 | 1.1e-05 |
| nearest | 0.3 | 60 | (20, 20, 20) | 0.300 → 4.01e-02 | 2.802e-04 | 1.291, 1.291, 1.291, 1.291 | 9.000e-02 → 8.71e-04 | 4.5e-06 |
| nearest | 0.3 | 240 | (80, 80, 80) | 0.300 → 3.42e-02 | 2.888e-04 | 1.110, 1.120, 1.120, 1.120 | 9.000e-02 → 9.24e-04 | 1.7e-06 |
| nearest | 1 | 3 | (1, 1, 1) | 1.000 → 6.25e-02 | 2.948e-04 | 1.000, 1.000, 1.000, 1.000 | 1.000e+00 → 2.29e-03 | 5.9e-05 |
| nearest | 1 | 15 | (5, 5, 5) | 1.000 → 1.32e-01 | 2.756e-04 | 1.168, 1.228, 1.257, 1.266 (maximum at the last record) | 1.000e+00 → 9.50e-03 | 2.4e-05 |
| nearest | 1 | 60 | (20, 20, 20) | 1.000 → 1.07e-01 | 2.967e-04 | 1.106, 1.106, 1.106, 1.106 | 1.000e+00 → 9.18e-03 | 1.5e-05 |
| nearest | 1 | 240 | (80, 80, 80) | 1.000 → 1.11e-01 | 2.934e-04 | 1.137, 1.138, 1.138, 1.138 | 1.000e+00 → 8.62e-03 | 9.8e-06 |
| nearest | 2.5 | 3 | (1, 1, 1) | 2.500 → 1.13e-01 | 4.062e-04 | 1.000, 1.000, 1.000, 1.000 | 6.250e+00 → 8.83e-03 | 2.0e-04 |
| nearest | 2.5 | 15 | (5, 5, 5) | 2.500 → 2.51e-01 | 2.900e-04 | 1.026, 1.026, 1.026, 1.026 | 6.250e+00 → 3.89e-02 | 6.2e-05 |
| nearest | 2.5 | 60 | (20, 20, 20) | 2.500 → 2.87e-01 | 2.942e-04 | 1.152, 1.159, 1.159, 1.159 | 6.250e+00 → 3.61e-02 | 2.4e-05 |
| nearest | 2.5 | 240 | (80, 80, 80) | 2.500 → 2.79e-01 | 2.971e-04 | 1.139, 1.143, 1.143, 1.143 | 6.250e+00 → 4.64e-02 | 3.5e-05 |
| nearest | 4.5 | 3 | (2, 0, 1) | 4.500 → 2.50e+00 | — (does not end at q) | — | 2.025e+01 → 1.23e+01 | 9.0e-03 |
| nearest | 4.5 | 15 | (6, 5, 4) | 4.500 → 1.16e+00 | — (does not end at q) | — | 1.866e+01 → 2.81e+00 | 1.2e-03 |
| nearest | 4.5 | 60 | (23, 19, 18) | 4.500 → 1.27e+00 | — (does not end at q) | — | 1.836e+01 → 2.30e+00 | 7.2e-04 |
| nearest | 4.5 | 240 | (87, 70, 83) | 4.500 → 8.93e-01 | — (does not end at q) | — | 1.832e+01 → 1.92e+00 | 5.4e-04 |
| ot | 0.3 | 3 | (1, 1, 1) | 0.300 → 1.04e-02 | 3.484e-04 | 1.000, 1.000, 1.000, 1.000 | 9.000e-02 → 8.69e-05 | 1.4e-05 |
| ot | 0.3 | 15 | (5, 5, 5) | 0.300 → 3.63e-02 | 2.694e-04 | 1.113, 1.114, 1.122, 1.162 (maximum at the last record) | 9.000e-02 → 8.29e-04 | 1.1e-05 |
| ot | 0.3 | 60 | (20, 20, 20) | 0.300 → 4.01e-02 | 2.802e-04 | 1.291, 1.291, 1.291, 1.291 | 9.000e-02 → 8.71e-04 | 4.5e-06 |
| ot | 0.3 | 240 | (80, 80, 80) | 0.300 → 3.42e-02 | 2.888e-04 | 1.110, 1.120, 1.120, 1.120 | 9.000e-02 → 9.24e-04 | 1.7e-06 |
| ot | 1 | 3 | (1, 1, 1) | 1.000 → 6.25e-02 | 2.948e-04 | 1.000, 1.000, 1.000, 1.000 | 1.000e+00 → 2.29e-03 | 5.9e-05 |
| ot | 1 | 15 | (5, 5, 5) | 1.000 → 1.32e-01 | 2.756e-04 | 1.168, 1.228, 1.257, 1.266 (maximum at the last record) | 1.000e+00 → 9.50e-03 | 2.4e-05 |
| ot | 1 | 60 | (20, 20, 20) | 1.000 → 1.07e-01 | 2.967e-04 | 1.106, 1.106, 1.106, 1.106 | 1.000e+00 → 9.18e-03 | 1.5e-05 |
| ot | 1 | 240 | (80, 80, 80) | 1.000 → 1.11e-01 | 2.934e-04 | 1.137, 1.138, 1.138, 1.138 | 1.000e+00 → 8.62e-03 | 9.8e-06 |
| ot | 2.5 | 3 | (1, 1, 1) | 2.500 → 1.13e-01 | 4.062e-04 | 1.000, 1.000, 1.000, 1.000 | 6.250e+00 → 8.83e-03 | 2.0e-04 |
| ot | 2.5 | 15 | (5, 5, 5) | 2.500 → 2.51e-01 | 2.900e-04 | 1.026, 1.026, 1.026, 1.026 | 6.250e+00 → 3.89e-02 | 6.2e-05 |
| ot | 2.5 | 60 | (20, 20, 20) | 2.500 → 2.87e-01 | 2.942e-04 | 1.152, 1.159, 1.159, 1.159 | 6.250e+00 → 3.61e-02 | 2.4e-05 |
| ot | 2.5 | 240 | (80, 80, 80) | 2.500 → 2.79e-01 | 2.971e-04 | 1.139, 1.143, 1.143, 1.143 | 6.250e+00 → 4.64e-02 | 3.5e-05 |
| ot | 4.5 | 3 | (1, 1, 1) | 4.500 → 1.57e-01 | 7.377e-04 | 1.379, 1.963, 1.963, 1.963 | 2.025e+01 → 1.01e-02 | 2.1e-04 |
| ot | 4.5 | 15 | (5, 5, 5) | 4.500 → 9.20e-01 | 3.017e-04 | 1.662, 2.107, 2.107, 2.107 | 1.866e+01 → 2.26e-01 | 2.1e-04 |
| ot | 4.5 | 60 | (20, 20, 20) | 4.500 → 7.55e-01 | 2.690e-04 | 1.470, 1.569, 1.569, 1.611 (maximum at the last record) | 1.836e+01 → 1.65e-01 | 1.2e-04 |
| ot | 4.5 | 240 | (80, 80, 80) | 4.500 → 8.77e-01 | 2.474e-04 | 1.509, 1.642, 1.783, 1.871 (maximum at the last record) | 1.832e+01 → 1.56e-01 | 1.5e-04 |

Inside the cells (ρ = 0.3, 1.0, 2.5; both leaks): the two leaks give identical runs (largest difference of d_∞ 0.0e+00); for m ≥ 15 fitted rates 2.694e-04–2.971e-04 against σ₀/k_e = 2.830e-04 (ratios 0.952–1.050), C(8000) = 1.026–1.291; for m = 3 fitted rates 2.948e-04–4.062e-04 against λ_coll/k_e = 2.952e-04 (ratios 0.999–1.376), C(8000) = 1.000–1.000. C is still growing at the horizon (maximum of d_∞e^{μt} at the last record) for: nearest, ρ = 0.3, m = 15; nearest, ρ = 1, m = 15; ot, ρ = 0.3, m = 15; ot, ρ = 1, m = 15.

m = 3: the final deviation from the exact representation on the three slowest collective modes (amplitudes in the whitened, weighted coordinates √w_k G^{−1/2}u_k; modes at λ_j(H_q)/k_e = 2.952e-04, 3.616e-04, 4.076e-04): ρ = 0.3: 0.0060, 0.0059, 0.0039 (fitted rate 3.484e-04); ρ = 1: 0.0364, 0.0302, 0.0066 (fitted rate 2.948e-04); ρ = 2.5: 0.0287, 0.0237, 0.0860 (fitted rate 4.062e-04). A fitted rate is that of the mode that carries the end of the run: the runs are mixtures of modes, and only a run whose end is carried by the slowest mode shows λ_coll/k_e.

### A′. Across a cell boundary (ρ = 4.5): the cell bookkeeping

Classes by nearest prior atom. For the nearest-atom leak: the weighted stationary point for the realised end weights (Newton by continuation), its residual and splitting margin σ₀ − max κ_R (negative: the lift is a saddle), the largest distance of the final configuration from the lift of that point and the lift's W² to ν. For the OT leak: the particles that end outside the cell of their initial class (the plan relabels) and the log-slope of W² over t > 6000 per slow unit, against the leak's 2σ₀ = 0.01415.

| leak | m | particles outside their class cell at t = 0 | counts at t = 0 → end | count changes (records) / last change at t | within-class spread at the end, per class (rms / max) | end weights: Newton residual, splitting margin | max distance from the lift of the end weights | W²(μ_Θ, ν) end / W² of that lift | particles ending outside their initial class's cell | log-slope of W² over t > 6000 (per slow unit) |
|---|---|---|---|---|---|---|---|---|---|---|
| nearest | 3 | 1 | (2, 0, 1) → (2, 0, 1) | 0 / — | 0.78 / 0.78; 0.00 / 0.00 | no stationary point found (residual 4.4e-02) | — | 12.259 / — | 1 | — |
| nearest | 15 | 3 | (6, 5, 4) → (6, 5, 4) | 0 / — | 0.66 / 0.94; 0.62 / 0.99; 0.24 / 0.40 | 1.9e-16, -0.00348 | 1.211 | 2.809 / 3.213 | 3 | — |
| nearest | 60 | 13 | (21, 22, 17) → (23, 19, 18) | 3 / 420 | 0.55 / 0.81; 0.46 / 1.19; 0.34 / 0.49 | 1.1e-16, -0.00161 | 1.251 | 2.303 / 2.345 | 12 | — |
| nearest | 240 | 55 | (80, 80, 80) → (87, 70, 83) | 12 / 285 | 0.46 / 0.72; 0.39 / 0.84; 0.39 / 0.53 | 2.0e-16, +0.00167 | 0.838 | 1.921 / 1.928 | 56 | — |
| ot | 3 | 1 | (2, 0, 1) → (1, 1, 1) | 1 / 2555 | 0.000 / 0.000; 0.000 / 0.000; 0.000 / 0.000 | — | — | 0.010 / — | 0 | -0.02338 |
| ot | 15 | 3 | (6, 5, 4) → (5, 5, 5) | 2 / 2695 | 0.300 / 0.439; 0.522 / 0.588; 0.514 / 0.916 | — | — | 0.226 / — | 4 | -0.01411 |
| ot | 60 | 13 | (21, 22, 17) → (20, 20, 20) | 8 / 1610 | 0.321 / 0.459; 0.451 / 0.765; 0.394 / 0.669 | — | — | 0.165 / — | 15 | -0.01408 |
| ot | 240 | 55 | (80, 80, 80) → (80, 80, 80) | 21 / 1445 | 0.344 / 0.453; 0.439 / 0.870; 0.374 / 0.515 | — | — | 0.156 / — | 53 | -0.01417 |

## B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (non-zero residual), masses q

Θ*(q) by Newton (residual 4.3e-17); prediction λ_sep = 0.00567, λ_coll = 0.00740, rate 2.270e-04 per time unit (the separation exponents (σ₀ − λ_j(R^G_k))/k_e: class 0: 2.740e-04, 2.856e-04, 3.292e-04; class 1: 2.597e-04, 2.836e-04, 3.193e-04; class 2: 2.270e-04, 2.796e-04, 3.037e-04). The reference is the 3-atom loop from Θ*(q); d_∞ is measured against the lift of the reference trajectory.

| ρ | m | d_∞(0) → d_∞(end) | fitted rate | C(T) at T = 2000, 4000, 6000, 8000 | rms tracking error (m-system / reference) |
|---|---|---|---|---|---|
| 0.3 | 15 | 0.300 → 6.02e-02 | 2.184e-04 | 1.139, 1.189, 1.216, 1.232 (maximum at the last record) | 2.55e-04 / 2.55e-04 |
| 1 | 15 | 1.000 → 1.36e-01 | 2.632e-04 | 1.078, 1.078, 1.078, 1.078 | 2.55e-04 / 2.55e-04 |
| 0.3 | 60 | 0.300 → 5.08e-02 | 2.339e-04 | 1.139, 1.139, 1.139, 1.139 | 2.55e-04 / 2.55e-04 |
| 1 | 60 | 1.000 → 1.48e-01 | 2.290e-04 | 1.079, 1.079, 1.079, 1.079 | 2.56e-04 / 2.55e-04 |
| 0.3 | 240 | 0.300 → 4.82e-02 | 2.242e-04 | 1.040, 1.040, 1.040, 1.040 | 2.55e-04 / 2.55e-04 |
| 1 | 240 | 1.000 → 1.53e-01 | 2.336e-04 | 1.122, 1.122, 1.122, 1.122 | 2.56e-04 / 2.55e-04 |

Fitted rates 2.184e-04–2.632e-04 against the prediction 2.270e-04 (ratios 0.962–1.159): the second-half fits are mixtures of the separation modes, whose gaps (inverse gaps above 10⁴ time units) exceed the fit window; they are consistent with the bound (1 − κ)min{λ_coll, λ_sep}/k_e, and the exponents of Part F are the test of the rate.

## C. One particle across a boundary (Remark D‴.19): the lift of q with one particle of class 0 moved 0.6 of the way to atom 1

| leak | m | mean-square displacement (1/m)Σ‖u_i‖² at t = 0 | counts at the end | W²(μ_Θ, ν): start → end | d_∞ at the end | end weights: Newton residual, splitting margin | max distance from the lift of the end weights | W² of that lift |
|---|---|---|---|---|---|---|---|---|
| nearest | 60 | 0.3114 | (19, 21, 20) | 0.3114 → 7.93e-01 | 3.90e-01 | 2.5e-15, +0.00561 | 0.307 | 0.8655 |
| ot | 60 | 0.3114 | (20, 20, 20) | 0.3114 → 3.03e-03 | 4.11e-01 | — | — | — |
| nearest | 240 | 0.0779 | (79, 81, 80) | 0.0779 → 1.98e-01 | 3.23e-01 | 1.1e-16, +0.00667 | 0.302 | 0.2163 |
| ot | 240 | 0.0779 | (80, 80, 80) | 0.0779 → 8.22e-04 | 4.40e-01 | — | — | — |

## F. The exponents that test the rates, and the proof's radius on this instance

**A (matching prior).** The 3-atom orbit of the real loop by shooting: residual 3.4e-16, period-start offset from Θ*(q) 0.00e+00; its parameter Floquet exponents 2.9355e-04, 3.5579e-04, 3.9986e-04, 8.7324e-04, 1.1091e-03, 1.5731e-03, 3.7043e-03, 5.5228e-03, 7.9337e-03 against λ_j(H_q)/k_e = 2.9521e-04, 3.6155e-04, 4.0755e-04, 8.8892e-04, 1.1271e-03, 1.5821e-03, 3.8276e-03, 5.6508e-03, 8.0914e-03 (relative errors 0.56–3.22%). The ghosts' exponents (separation exponents of every lift, Proposition D‴.15): class 0: 2.8297e-04, 2.8298e-04, 2.8298e-04 (prediction 2.8298e-04, 2.8298e-04, 2.8298e-04; return residual 0.0e+00); class 1: 2.8298e-04, 2.8298e-04, 2.8298e-04 (prediction 2.8298e-04, 2.8298e-04, 2.8298e-04; return residual 5.9e-17); class 2: 2.8297e-04, 2.8298e-04, 2.8298e-04 (prediction 2.8298e-04, 2.8298e-04, 2.8298e-04; return residual 0.0e+00). The slowest exponent overall: 2.8297e-04.

The proof's radius: a sampled lower bound C_N ≥ 0.0925 on the |||·|||-Lipschitz constant of D⟨H⟩ at Θ*_15(q) (displacements of per-particle size 0.05) gives ρ ≤ λ₀/(4C_N) = 0.0191 with λ₀ = 0.00707; the runs start at per-particle radii 0.3, 1, 2.5, i.e. 31–261 times the proved start radius ρ/2. The contraction observed there is evidence that the basin extends to the cell margins, which the theorem does not claim.

**B (shifted prior).** The 3-atom orbit of the real loop by shooting: residual 1.4e-14, period-start offset from Θ*(q) 2.68e-04; its parameter Floquet exponents 2.9422e-04, 3.5373e-04, 4.0064e-04, 8.6899e-04, 1.1136e-03, 1.5758e-03, 3.6478e-03, 5.4843e-03, 8.0248e-03 against λ_j(H_q)/k_e = 2.9591e-04, 3.5944e-04, 4.0853e-04, 8.8481e-04, 1.1316e-03, 1.5851e-03, 3.7665e-03, 5.6261e-03, 8.1826e-03 (relative errors 0.57–3.15%). The ghosts' exponents (separation exponents of every lift, Proposition D‴.15): class 0: 2.7401e-04, 2.8565e-04, 3.2915e-04 (prediction 2.7397e-04, 2.8565e-04, 3.2918e-04; return residual 1.1e-14); class 1: 2.5965e-04, 2.8356e-04, 3.1927e-04 (prediction 2.5965e-04, 2.8357e-04, 3.1926e-04; return residual 8.9e-15); class 2: 2.2701e-04, 2.7959e-04, 3.0367e-04 (prediction 2.2696e-04, 2.7960e-04, 3.0375e-04; return residual 8.0e-15). The slowest exponent overall: 2.2701e-04.

The proof's radius: a sampled lower bound C_N ≥ 0.0922 on the |||·|||-Lipschitz constant of D⟨H⟩ at Θ*_15(q) (displacements of per-particle size 0.05) gives ρ ≤ λ₀/(4C_N) = 0.0154 with λ₀ = 0.00567; the runs start at per-particle radii 0.3, 1, i.e. 39–130 times the proved start radius ρ/2. The contraction observed there is evidence that the basin extends to the cell margins, which the theorem does not claim.

Numbers are produced by `experiments/exp_uniform_basin.py`; nothing is edited by hand.
