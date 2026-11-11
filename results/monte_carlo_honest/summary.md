# Monte Carlo over random initial conditions (Van der Pol, m=16, k_e=5.0, 200 trials)

| check | result |
|---|---|
| all trajectories finite | 200/200 |
| hard parameter bounds held at every recorded step | 200/200 |
| a-posteriori ISS bound ‖e(t)‖ ≤ e^{-kₑ(t-T₁)}‖e(T₁)‖ + sup‖f−f̂‖/kₑ held | 200/200 (max ratio 0.979) |
| steady-state max ‖e‖ (last 25%): median / 95th pct / max | 0.0568 / 0.0938 / 0.1075 |
| steady-state mean ‖e‖: median / 95th pct / max | 0.0302 / 0.0484 / 0.0593 |
| sup‖f − f̂‖ in steady state / kₑ: median / max | 0.0627 / 0.1475 |
| transient max ‖e‖: max over trials | 4.903 |

No a-priori R_UUB threshold is drawn: D_max and C_Ω are not computable for this plant (the 0.5 threshold of the original figure was an arbitrary constant, not a calculated bound).
