# Mean-field consistency in m (Van der Pol, mean-field gains, smooth initial measure quantised at m points)

Theorem D predicts sup_t ‖e^m − e^{2m}‖ + W(μ^m_t, μ^{2m}_t) ≤ e^{Λt} W(μ^m_0, μ^{2m}_0); the initial distance decays like 1/m.

| m | sup_t ‖e^m − e^{2m}‖ | W(μ^m_0, μ^{2m}_0) | W(μ^m_T, μ^{2m}_T) | ratio W_T / W_0 | RMS ‖e‖ (last half) |
|---|---|---|---|---|---|
| 16 vs 32 | 0.0046 | 0.2179 | 0.3658 | 1.68 | 0.2842 |
| 32 vs 64 | 0.0009 | 0.1091 | 0.1788 | 1.64 | 0.2832 |
| 64 vs 128 | 0.0001 | 0.0545 | 0.0920 | 1.69 | 0.2832 |
| 128 vs 256 | 0.0000 | 0.0273 | 0.0459 | 1.68 | 0.2832 |
| 256 | — | — | — | — | 0.2832 |
