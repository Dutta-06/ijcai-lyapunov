# Variable-projection check — configuration 'satisfying' (m=2, k_e=20, g_a=13, δ=0.004, T=3000 s)

Time scales along the run (min–max over 5 samples; rates in 1/s):

| k_e | ω_ref | joint fastest rate (H2a) | fastest weight | slowest weight | geometry (mδ/k_e)‖Γ^{1/2}∇²L*Γ^{1/2}‖ | slope S (H2c) |
|---|---|---|---|---|---|---|
| 20 | 1 | 0.138–0.153 | 0.137–0.153 | 0.039–0.04 | 0.000458–0.000494 | 0.00601–0.00602 |

(H2a) ω_ref / joint fastest: 6.53–7.27; (H2b) slowest weight / geometry: 78.8–87.2; (H2c) S ≪ 1: S = 0.00601–0.00602.

| quantity | first quarter (median) | last quarter (median) |
|---|---|---|
| G(a): cosine(joint velocity, −(m/k_e)G∇L̄), G⁻¹ metric | 0.9848 | 0.9624 |
| G(a): magnitude ratio | 0.956 | 0.935 |
| G(a), weight block only: cosine | -0.2128 | 0.5004 |
| G(a), geometry block only: cosine | 0.9992 | 0.9982 |
| G(a): share of the residual ‖v − v_pred‖²_(G⁻¹) carried by the weight block | 0.922 | 0.947 |
| G(a): predicted weight-block / geometry-block speed, G⁻¹ norm | 0.165 | 0.279 |
| G(b): cosine(geometry velocity, −(mδ/k_e)Γ∇L*) | 0.999 | 0.999 |
| G(b): magnitude ratio | 0.969 | 0.970 |
| relative gap of the weights to a*(c,s) | 0.001 | 0.001 |

| quantity | first sample | last sample |
|---|---|---|
| profile loss L* = min_a L̄ | 1.1987e-02 | 4.1073e-03 |
| L̄ at the period-mean weights | 1.1998e-02 | 4.1076e-03 |
| ‖∇L*‖ | 1.405e-01 | 4.417e-02 |
| tracking RMS over a period | 0.0077 | 0.0044 |

Sampled once per period (centred differences over ±20 periods): L* non-increasing: **True**; L̄ (period-mean weights) non-increasing: **True** (fraction of decreasing steps 1.00); fraction of samples with G(b) cosine > 0.9: 1.00; with G(a) cosine > 0.99: 0.00
