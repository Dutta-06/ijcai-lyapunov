# Variable-projection check — configuration 'violating' (m=4, k_e=20, g_a=2, δ=0.1, T=800 s)

Time scales along the run (min–max over 5 samples; rates in 1/s):

| k_e | ω_ref | joint fastest rate (H2a) | fastest weight | slowest weight | geometry (mδ/k_e)‖Γ^{1/2}∇²L*Γ^{1/2}‖ | slope S (H2c) |
|---|---|---|---|---|---|---|
| 20 | 1 | 0.039–0.0656 | 0.0126–0.0163 | 0.000782–0.00147 | 0.0047–0.0143 | 5.44–25 |

(H2a) ω_ref / joint fastest: 15.3–25.6; (H2b) slowest weight / geometry: 0.0993–0.231; (H2c) S ≪ 1: S = 5.44–25.

| quantity | first quarter (median) | last quarter (median) |
|---|---|---|
| G(a): cosine(joint velocity, −(m/k_e)G∇L̄), G⁻¹ metric | 0.9994 | 0.9985 |
| G(a): magnitude ratio | 1.002 | 0.939 |
| G(a), weight block only: cosine | 0.9996 | 0.9985 |
| G(a), geometry block only: cosine | 0.9997 | 0.9986 |
| G(a): share of the residual ‖v − v_pred‖²_(G⁻¹) carried by the weight block | 0.601 | 0.423 |
| G(a): predicted weight-block / geometry-block speed, G⁻¹ norm | 1.76 | 1.02 |
| G(b): cosine(geometry velocity, −(mδ/k_e)Γ∇L*) | -0.118 | 0.385 |
| G(b): magnitude ratio | 4.207 | 0.472 |
| relative gap of the weights to a*(c,s) | 0.799 | 0.206 |

| quantity | first sample | last sample |
|---|---|---|
| profile loss L* = min_a L̄ | 3.8136e-03 | 5.8866e-03 |
| L̄ at the period-mean weights | 2.1892e-01 | 8.0127e-03 |
| ‖∇L*‖ | 2.135e-02 | 3.733e-02 |
| tracking RMS over a period | 0.0329 | 0.0060 |

Sampled once per period (centred differences over ±5 periods): L* non-increasing: **False**; L̄ (period-mean weights) non-increasing: **True** (fraction of decreasing steps 1.00); fraction of samples with G(b) cosine > 0.9: 0.00; with G(a) cosine > 0.99: 1.00
