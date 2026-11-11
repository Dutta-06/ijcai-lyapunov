# Proposition G(b) needs slow-manifold flatness: closed loop vs averaged predictions

m = 1, n = 1, f(x) = 2exp(−(x−1)²/0.5), x_m = 0.2 sin t, k_e = 20, gains (g_a, δ, δ).  The closed loop and the joint averaged flow start ON the slow manifold of (b′) at ξ₀ = (1.75, 0.9) (a*(ξ₀) = 2.055).  Displacements of the period-mean geometry between the centres of the first and last periods.

| case | horizon (s) | (H2b) ratio | slope S (H2c) | ω_ref / joint fastest rate (H2a) | closed-loop Δξ = (Δc, Δs) | (a) joint averaged Δξ | (b) VarPro Δξ | (b′) induced-metric Δξ |
|---|---|---|---|---|---|---|---|---|
| δ = 0.47 | 232 | 81–178 | 2–2.81 | 12.4–13.0 | (-0.0037, -0.0278) | (-0.0043, -0.0250) | (+0.0161, -0.0645) | (-0.0046, -0.0246) |
| δ = 0.47, all gains / 4 | 928 | 81–178 | 2–2.81 | 49.5–51.8 | (-0.0047, -0.0274) | (-0.0045, -0.0261) | (+0.0168, -0.0672) | (-0.0048, -0.0256) |
| δ = 0.047 | 1238 | 654–1779 | 0.2–0.319 | 31.5–37.7 | (+0.0063, -0.0315) | (+0.0060, -0.0302) | (+0.0093, -0.0367) | (+0.0060, -0.0301) |
| δ = 0.0047 | 3000 | 13162–17790 | 0.02–0.0229 | 36.9–39.8 | (+0.0023, -0.0092) | (+0.0022, -0.0088) | (+0.0023, -0.0090) | (+0.0022, -0.0088) |

| case | ‖closed‖/‖(b)‖ | cos(closed, (b)) | ‖closed‖/‖(b′)‖ | cos(closed, (b′)) | ‖closed‖/‖(a)‖ | cos(closed, (a)) | ‖(a)‖/‖(b′)‖ | max |a| along the run (box 20) |
|---|---|---|---|---|---|---|---|---|
| δ = 0.47 | 0.422 | 0.9303 | 1.121 | 0.9986 | 1.104 | 0.9992 | 1.015 | 2.25 |
| δ = 0.47, all gains / 4 | 0.401 | 0.9156 | 1.065 | 0.9999 | 1.049 | 1.0000 | 1.015 | 2.24 |
| δ = 0.047 | 0.850 | 0.9986 | 1.047 | 1.0000 | 1.044 | 1.0000 | 1.003 | 2.38 |
| δ = 0.0047 | 1.023 | 1.0000 | 1.044 | 1.0000 | 1.044 | 1.0000 | 1.000 | 2.16 |

Reading: compare the cosine and size ratio of the closed loop against (b) and against (b′), row by row, together with S.  (H2b) is satisfied in every row.  ‖(a)‖/‖(b′)‖ isolates the O(1/R) error of (b′); ‖closed‖/‖(a)‖ isolates the averaging and error-channel errors (compare the 'gains / 4' row).
