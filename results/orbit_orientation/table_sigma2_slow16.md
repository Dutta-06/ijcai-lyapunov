# Orientation of the reference orbit changes the learned kernels (Proposition H.7, real closed loop)

n = 2, m = 6, G = diag(4.0, 4.0, 1.0, 1.0, 1.0), σ = 2/k_e, no dither, x_m = (1.2 cos t, ±0.8 sin t) (the same occupation measure for both orientations), same initial parameters, horizon 16·k_e periods.  Distances are labelled G⁻¹ distances between period-mean parameters.

Two linear-response predictions of Θ₊ − Θ₋ at the (+) equilibrium: *frozen* uses only the frozen-measure correction of Proposition H.7 (−2J⁻¹b_odd); *full* adds the second-order averaging term of the intra-period parameter oscillation (−2J⁻¹(b_odd + b_osc)).

| k_e | measured ‖Θ₊ − Θ₋‖ | residual drift, last 10 periods (+, −) | frozen: size ratio measured/predicted | frozen: cosine | full: size ratio measured/predicted | full: cosine | full: vector error ‖meas − pred‖/‖pred‖ |
|---|---|---|---|---|---|---|---|
| 10 | 5.729e-03 | 1.1e-16, 1.0e-16 | 2.670 | 0.7577 | 0.937 | 0.9314 | 0.364 |
| 20 | 2.843e-03 | 3.9e-16, 5.0e-16 | 2.414 | 0.8302 | 0.965 | 0.9809 | 0.195 |
| 40 | 1.403e-03 | 1.8e-16, 1.3e-16 | 2.294 | 0.8631 | 0.969 | 0.9944 | 0.109 |

Fitted exponent of the measured difference against k_e: -1.01 (prediction: −1).
Kernel spread (largest coordinate deviation of a kernel from the mean kernel): 1.05 at the start, at most 5.6e-17 at the end of every run. The 6 kernels converge to a common configuration, so the learned network is effectively a single kernel and the test probes the orientation effect at a collapsed equilibrium of the m = 6 system.
Reading: the occupation measure is identical for the two orientations, so (𝔄) cannot tell them apart. The difference is the orientation-odd part of the first-order correction of the real averaged loop; compare the frozen and full columns to see how much of it the frozen-measure expansion alone accounts for.
