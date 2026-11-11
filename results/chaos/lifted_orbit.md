# Proposition D‴.15: the transverse spectrum of a lifted periodic orbit of the real closed loop (`exp_lifted_orbit.py`)

The real closed loop (CL) of Theorem J on the three-bump instance at k_e = 25, σ = σ₀/k_e with σ₀ = 0.15σ_crit = 0.00707; monodromy matrices by central differences of the Poincaré map (DOP853, rtol 10⁻¹¹); Floquet exponents −log|λ|/T₀ (the error channel's exponent sits at the finite-difference floor and is omitted). The lifts are the collapsed configurations with the K₀-atom orbit's atoms repeated n_k times; the ghost of class k is one particle with weight zero in f̂, driven by the K₀-atom orbit and leaking toward the class target.

## A. The OT leak at the matching masses (prior = exact representation, q = (1/3, 1/3, 1/3))

The 3-atom loop's periodic orbit is the exact representation with zero tracking error (shooting residual 3.4e-16); its parameter Floquet exponents are 2.935e-04, 3.558e-04, 3.999e-04, 8.732e-04, 1.109e-03, 1.573e-03, 3.704e-03, 5.523e-03, 7.934e-03, against Theorem J(b)'s prediction λ_j(H_q)/k_e = 2.952e-04, 3.616e-04, 4.076e-04, 8.889e-04, 1.127e-03, 1.582e-03, 3.828e-03, 5.651e-03, 8.091e-03 with H_q the whitened, weighted collective block of Proposition D‴.5(i) (relative errors 0.56–3.22%, the O(k_e⁻²) of Theorem J(b)). The ghosts' exponents: class 0: 2.830e-04, 2.830e-04, 2.830e-04 (return residual 0.0e+00); class 1: 2.830e-04, 2.830e-04, 2.830e-04 (return residual 5.9e-17); class 2: 2.830e-04, 2.830e-04, 2.830e-04 (return residual 0.0e+00) — σ₀/k_e = 2.830e-04 at zero residual.

| m | classes | return residual of the lift | separation exponents per class (distinct values, rounded to 10⁻⁷; each of multiplicity n_k − 1) | invariance defects | collective exponents (parameters) | collective defect |
|---|---|---|---|---|---|---|
| 6 | (2, 2, 2) | 1.7e-15 | class 0: 2.830e-04 ×1; class 1: 2.830e-04 ×1; class 2: 2.830e-04 ×1 | 9.1e-17, 7.5e-09, 6.4e-09 | 2.935e-04, 3.558e-04, 3.999e-04, 8.733e-04, 1.109e-03, 1.573e-03, 3.704e-03, 5.523e-03, 7.934e-03 | 2.8e-09 |
| 9 | (3, 3, 3) | 1.3e-15 | class 0: 2.830e-04 ×2; class 1: 2.830e-04 ×2; class 2: 2.830e-04 ×2 | 1.8e-08, 1.3e-08, 8.0e-09 | 2.935e-04, 3.558e-04, 3.999e-04, 8.732e-04, 1.109e-03, 1.573e-03, 3.704e-03, 5.523e-03, 7.934e-03 | 2.2e-09 |
| 15 | (5, 5, 5) | 7.1e-16 | class 0: 2.830e-04 ×4; class 1: 2.830e-04 ×4; class 2: 2.830e-04 ×4 | 9.1e-09, 1.5e-08, 1.0e-08 | 2.935e-04, 3.558e-04, 3.998e-04, 8.733e-04, 1.109e-03, 1.573e-03, 3.704e-03, 5.523e-03, 7.934e-03 | 2.1e-09 |

Relative errors over the three lifts: collective exponents against the 3-atom orbit's 0.0079% at most; separation exponents against the ghosts' σ₀/k_e 0.004% at most.

## B. The σ-modification (ϑ⁰ = 0) at unequal weights: the weighted 3-atom orbit, its lift to m = 16, and the saddle

### w = (0.4375, 0.3125, 0.2500), classes (7, 5, 4): splitting margin σ₀ − max κ_R = +0.00257 (admissible)

Weighted 3-atom orbit: shooting residual 1.1e-14, period-start offset from Θ*(w) 1.12e-02 (Theorem J(a): O(1/k_e)); its parameter exponents 2.492e-04, 3.730e-04, 8.157e-04, 8.600e-04, 1.436e-03, 2.387e-03, 3.055e-03, 5.303e-03, 8.308e-03. Lift to m = 16: return residual 4.9e-14; collective exponents 2.492e-04, 3.730e-04, 8.157e-04, 8.600e-04, 1.436e-03, 2.387e-03, 3.055e-03, 5.303e-03, 8.308e-03 (defect 2.1e-08).

| class k | n_k − 1 | separation exponents of the lift | ghost exponents | prediction (σ₀ − λ_j(R^G_k))/k_e | κ_R at the atom | defect |
|---|---|---|---|---|---|---|
| 0 | 6 | 9.989e-05, 4.553e-04, 1.267e-03 (each ×6) | 9.989e-05, 4.553e-04, 1.267e-03 | 1.029e-04, 4.572e-04, 1.275e-03 | +0.00450 | 4.7e-08 |
| 1 | 4 | 4.004e-04, 8.219e-04, 3.245e-03 (each ×4) | 4.004e-04, 8.219e-04, 3.245e-03 | 4.022e-04, 8.238e-04, 3.288e-03 | -0.00298 | 4.0e-08 |
| 2 | 3 | 3.230e-04, 5.542e-04, 1.631e-03 (each ×3) | 3.230e-04, 5.542e-04, 1.631e-03 | 3.261e-04, 5.566e-04, 1.652e-03 | -0.00108 | 5.7e-08 |
Relative errors: separation exponents of the lift against the ghosts' 0.002% at most, against the prediction (σ₀ − λ_j(R^G_k))/k_e 0.23–2.95%; collective exponents of the lift against the weighted 3-atom orbit's 0.0007% at most.

### w = (0.5625, 0.2500, 0.1875), classes (9, 4, 3): splitting margin σ₀ − max κ_R = -0.00219 (a saddle)

Weighted 3-atom orbit: shooting residual 7.1e-15, period-start offset from Θ*(w) 1.60e-02 (Theorem J(a): O(1/k_e)); its parameter exponents 2.131e-04, 5.371e-04, 8.019e-04, 9.409e-04, 1.302e-03, 2.235e-03, 3.171e-03, 5.913e-03, 9.985e-03. Lift to m = 16: return residual 4.5e-14; collective exponents 2.131e-04, 5.371e-04, 8.019e-04, 9.409e-04, 1.302e-03, 2.235e-03, 3.171e-03, 5.913e-03, 9.985e-03 (defect 1.9e-08).

| class k | n_k − 1 | separation exponents of the lift | ghost exponents | prediction (σ₀ − λ_j(R^G_k))/k_e | κ_R at the atom | defect |
|---|---|---|---|---|---|---|
| 0 | 8 | -9.126e-05, 4.059e-04, 1.188e-03 (each ×8) | -9.126e-05, 4.059e-04, 1.188e-03 | -8.769e-05, 4.089e-04, 1.195e-03 | +0.00927 | 5.3e-08 |
| 1 | 3 | 5.030e-04, 9.210e-04, 4.229e-03 (each ×3) | 5.030e-04, 9.210e-04, 4.229e-03 | 5.051e-04, 9.232e-04, 4.292e-03 | -0.00555 | 5.0e-08 |
| 2 | 2 | 5.086e-04, 6.648e-04, 2.276e-03 (each ×2) | 5.086e-04, 6.648e-04, 2.276e-03 | 5.119e-04, 6.676e-04, 2.315e-03 | -0.00572 | 4.6e-08 |
Relative errors: separation exponents of the lift against the ghosts' 0.001% at most, against the prediction (σ₀ − λ_j(R^G_k))/k_e 0.24–4.07%; collective exponents of the lift against the weighted 3-atom orbit's 0.0014% at most.

Reading. The lift of the K₀-atom orbit is a periodic orbit of the m-loop for every size tried (return residuals at the integration tolerance), its monodromy leaves every class's separation subspace invariant (defects at the finite-difference level) and acts there as the ghost's monodromy with multiplicity n_k − 1, while the complement carries the K₀-atom loop's spectrum: the Floquet spectrum of the lift does not depend on m. Theorem J(b) reads the separation exponents off the averaged Hessian, (σ₀ − λ_j(R^G_k))/k_e + O(k_e⁻²): at the matching masses the residual vanishes and every separation exponent is σ₀/k_e; at the saddle weights (0.5625, 0.2500, 0.1875), class 0: -9.126e-05 — the over-weighted class's negative exponent, the orbit-level form of Proposition K's split.

Numbers are produced by `experiments/exp_lifted_orbit.py`; nothing is edited by hand.
