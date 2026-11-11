# Proposition K: orientation-dependent stability of the collapsed configuration (`exp_orientation_stability.py`)

n = 2, m = 3, ellipse (1.2 cos t, ±0.8 sin t), G = diag(4.0, 4.0, 1.0, 1.0, 1.0). Floquet exponents of the collapsed periodic orbit of the real closed loop (Newton shooting; the parameter block of the monodromy matrix restricted to the separation subspace Σu_i = 0 and to the collective subspace). Exponents are −log|λ|/T₀: positive means stable. The leading-order separation exponents are the eigenvalues of (σ₀I − G^{1/2}R₂G^{1/2})/k_e, identical for both orientations. The second-order prediction is −eig of the separation block of the Jacobian of the second-order averaged drift (1/k_e)⟨H⟩ + b_frozen + b_osc at its own equilibrium Θ̄₂ (Newton from Θ*); the column 'at Θ*' uses the Jacobian at the collapsed critical point instead, which omits the third-derivative term ∂[(1/k_e)G∇J_m](Θ̄₂ − Θ*) = O(k_e⁻²). Every Θ̄₂ found is collapsed (kernel spread below 1e-8).

## 1. σ₀ = 0.5 (stable collapse)

| k_e | shooting residual (+, −) | softest separation exponent (+) | (−) | leading order (same for both) | k_e² × (measured − leading) (+, −) | k_e² × (predicted at Θ̄₂ − leading) (+, −) | difference (+) − (−) | k_e² × difference | predicted difference at Θ̄₂ | at Θ* | ‖Θ̄₂ − Θ*‖ (+) | softest collective exponent (+, −) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 1.0e-13, 4.2e-13 | 3.56262e-03 | 2.66009e-03 | 3.88551e-03 | -0.0323, -0.1225 | -0.0824, -0.1601 | +9.025e-04 | +0.0903 | +7.771e-04 | +9.011e-04 | 3.676e-02 | 4.8668e-02, 4.7604e-02 |
| 20 | 4.2e-15, 5.2e-15 | 1.80026e-03 | 1.61408e-03 | 1.94275e-03 | -0.0570, -0.1315 | -0.0829, -0.1506 | +1.862e-04 | +0.0745 | +1.693e-04 | +2.082e-04 | 1.807e-02 | 2.4472e-02, 2.4241e-02 |
| 40 | 2.1e-15, 2.2e-15 | 9.27681e-04 | 8.86260e-04 | 9.71377e-04 | -0.0699, -0.1362 | -0.0831, -0.1459 | +4.142e-05 | +0.0663 | +3.926e-05 | +4.982e-05 | 8.963e-03 | 1.2291e-02, 1.2237e-02 |
| 80 | 1.8e-15, 2.1e-15 | 4.73736e-04 | 4.64029e-04 | 4.85688e-04 | -0.0765, -0.1386 | -0.0831, -0.1435 | +9.708e-06 | +0.0621 | +9.435e-06 | +1.217e-05 | 4.464e-03 | 6.1614e-03, 6.1485e-03 |

Fitted exponent of |difference| against k_e: -2.18 (Proposition K: −2). Largest kernel spread of any Θ̄₂: 1.2e-08.

## 2. σ₀ = σ_crit = 0.460627: the softest separation mode is marginal at leading order

| k_e | softest separation exponent (+) | (−) | sign (+), sign (−) | predicted at Θ̄₂ (+), (−) | predicted at Θ* (+), (−) | second softest (+), (−) |
|---|---|---|---|---|---|---|
| 10 | -4.1307e-04 | -1.2738e-03 | unstable, unstable | -8.7333e-04, -1.6145e-03 | -6.4358e-04, -1.5209e-03 | -4.1307e-04, -1.2738e-03 |
| 20 | -1.6043e-04 | -3.3743e-04 | unstable, unstable | -2.1982e-04, -3.8079e-04 | -1.5392e-04, -3.5689e-04 | -1.6043e-04, -3.3743e-04 |
| 40 | -4.7533e-05 | -8.6854e-05 | unstable, unstable | -5.5061e-05, -9.2334e-05 | -3.7610e-05, -8.6213e-05 | -4.7533e-05, -8.6854e-05 |

Reading: at σ₀ = σ_crit the leading-order separation exponent vanishes and the k_e⁻² terms decide; the orientation-odd part (U_odd, b_osc and the third-derivative term they induce) differs between the two traversals, the orientation-even part does not. The two separation exponents coincide (multiplicity m − 1 = 2, Proposition K(ii)).

## 3. σ₀ designed from the prediction: the orientation decides

σ₀(k_e) solves pred₊(σ₀) + pred₋(σ₀) = 0 for the predicted softest separation exponents at Θ̄₂ (secant from σ_crit), so that the predicted exponents of the two orientations are opposite. The measured exponents are then computed by shooting.

| k_e | σ₀ − σ_crit | leading order (same for both) | predicted (+), (−) | measured (+), (−) | sign (+), sign (−) |
|---|---|---|---|---|---|
| 10 | +0.01251 | +1.2334e-03 | +3.7652e-04, -3.7651e-04 | +8.4970e-04, -2.4994e-05 | stable, unstable |
| 20 | +0.00606 | +2.9896e-04 | +8.1162e-05, -8.1154e-05 | +1.4138e-04, -3.7111e-05 | stable, unstable |

At every k_e of the table one orientation of the same reference orbit (same occupation measure, same averaged system) leaves the collapsed configuration stable and the other makes it unstable.

Numbers are produced by `experiments/exp_orientation_stability.py`; nothing is edited by hand.
