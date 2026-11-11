# Referee report D3z: Proposition D‴.37

Object: Proposition D‴.37 (PARADIGM.md §5d, and its §0/§7 rows, §8 bullet), results/RESULTS_PARADIGM.md §D‴.37, LITERATURE_TRANSPORT.md [S65] and its audit item, PAPERS paper 13, README; src/core/gauge_drift.py, experiments/exp_gauge_drift.py, results/chaos/gauge_drift.{md,_raw.npz}, tests/test_gauge_drift.py. Commit 391b29f. I changed nothing in the repository. My scripts are in scratchpad/ref27/.

## What I checked and how

1. **(i) The identities, re-derived and re-implemented without the eigenbasis.** I differentiated ∂Φ·F ≡ 0 twice at F = 0 and got ∂²Φ[u,Av] + ∂²Φ[v,Au] = −P∂²F[u,v]. From Φ∘γ = γ I got (I−P)∂²γ for the tangent–tangent block. I then wrote an independent drift (`ref27/indep.py`) that uses no eigendecomposition and no division by λ_i+λ_j:
   - P is built from the left kernel of A.
   - The QQ part is ⟨P∂²F, X⟩, with X from the Lyapunov equation (A−P)X + X(A−P)ᵀ = −QΣQᵀ.
   - The PQ part is ⟨P∂²F, Y+Yᵀ⟩, with Y = −PΣQᵀA^{#ᵀ} (A^# the group inverse).
   - The PP part is (I−P)∂²γ[K].
   - The Itô term is added.

   I also wrote my own brute-force limit map (`ref27/bf.py`): BDF to 2·10⁵, second differences at h = 0.02 and 0.01, then Richardson extrapolation. I ran it at two points not used in Part A: (0.27, 0.36, 0.37), which is inside the first class's layer, and (0.29, 0.22, 0.49).
2. **(ii) The slow limit.** I checked the sign and the Gibbs mean, and the closed form of ℓ (below). I checked how drift/(γτ) and D/(γ²τ) behave down to γ = 10⁻⁷. I checked split invariance for the compact kernel and for Gaussian kernels with R = 1 and R = 0.5, with the split t ∈ {0.5, 0.2, 0.05}.
3. **(iii) The layer.** I computed the rates with the actual Γ_11 = 0.2354. I checked the δ* crossing and the zero's sensitivity to the finite-difference steps (h, h2) ∈ {(3e-6, 1e-4), (1e-6, 3e-5), (1e-5, 3e-4), (3e-6, 1e-5)}. Along W_2 = W_3 I computed b_1, D_11 and the zero-flux quantity 2b − ∂_W D. I re-analysed the full-loop records H_rec: α excursions, W below W_min, and time fractions by quarter.
4. **(iv) Part G**, re-analysed per component in units of the standard error.
5. **Numerics.**
   - Part E: an independent zero search on a 2556-point triangulated grid with pad 0.002, using winding numbers per triangle, Newton, and the boundary normal component.
   - Part H: deviation of each mean in standard errors.
   - I regenerated the md from the raw file with `--from-raw` into a scratch OUT. It is byte-identical (`cmp` passes).
   - I checked every number in the text, rows, bullet and RESULTS against the md and the raw file.
6. **Tests.** I ran them (6 passed in 0.7 s). I made 9 mutants in a scratch copy.
7. **Literature.** One web search, on the Li–Wang–Arora isotropic corollary. Everything else is from memory and flagged as such.

## Verdicts

**(i) Correct.** The division by λ_i+λ_j is right. The mixed tangent/non-tangent case is right (it divides by λ_j, which is what the formula gives with λ_i = 0). The (I−P)∂²γ block is right. The Katzenberger drift ½∂²Φ[Σ] with Σ = 2τG is right, and so is the Itô term ½∂²W[PΣPᵀ]. Checks:
- My Lyapunov/group-inverse implementation matches `limit_drift`:

  | W_1 − 0.25 | relative difference |
  |---|---|
  | 1/12 (interior) | 2·10⁻⁶ |
  | 0.02 | 4·10⁻⁵ |
  | 0.005 | 6·10⁻⁴ |
  | 0.002 | 4·10⁻³ (step-limited) |

- My brute force matches the identities to 5.9·10⁻⁵ at (0.27, 0.36) and 3.7·10⁻⁶ at (0.29, 0.22), after Richardson.
- D matches D‴.36's projection to ≤ 4·10⁻⁷.
- Size of each term at W = 1/3 (the full drift's ratio after dropping it):
  - Itô term: 4–5% of the drift.
  - Tangent–tangent term: 14–19%.
  - Mixed PQ term: almost all of the drift (dropping it leaves 0.005 and 0.15 of the drift).

  At (0.27, 0.36), dropping the tangent–tangent term flips the sign of the drift.

**(ii) Correct, and simpler than stated.**
- The signs are right: ⟨−∂_βJ⟩_Gibbs = τ∂_β log Z. Laplace gives −(τ/2)∂ℓ. The replicator form is right, and ℓ is correctly the total β-derivative through Θ*(β).
- Closed form, not in the text: with one particle per class the Hessian at the exact fit is m·Gram, and only the α-rows depend on W (∝ W_k a′_k). So exactly ℓ = 2Σ_k log(W_k − W_min,k²/W_k) + const. I confirmed this against `log_det` to 10⁻⁹: the constant is −19.53724 at all five points. The law is then u_k ∝ −(W_k² + W_min,k²)/(W_k(W_k² − W_min,k²)), which does not depend on the bumps' shapes. The repository's finite-difference `log_det_law` is off from this closed form by 5·10⁻⁵ at W_1 = 0.28 and 4·10⁻⁴ at 0.26.
- drift/(γτ) converges to the closed form with an O(γ) relative error: at W = (0.28, 0.36, 0.36) the W_1 component is −1.274, −1.751, −1.794, −1.7986, −1.7990 for γ = 10⁻³ to 10⁻⁷, against the closed form's −1.7990.
- D/(γ²τ) converges (5.199 at 1/3, 36.02 at 0.28), so D = O(γ²τ) holds. Deterministic motion on the time scale 1/(γτ) with spread O(√γ) is right.
- Split invariance holds: |Δℓ| ≤ 6·10⁻⁸ for the compact kernel and ≤ 5·10⁻⁸ for Gaussian kernels with R = 1 and 0.5. The Gaussian cross-class cohesion is exponentially small on this instance; in general it adds split-dependent cross terms of order σ₀ψ(d_cross) and makes the exact representation non-stationary.
- The metastability point is a wording issue only (ID-7).

**(iii) Partly right; the main conclusion is wrong.**
- What holds: ℓ → −∞ (indeed ℓ ∝ 2 log a′). Both rates re-derive exactly. With the true Γ_11, the amplitude channel 3G_αΓW²a′² crosses the mass channel γa²WΓ at 0.98δ*, so δ* is the right scale. The zero is robust to the finite-difference steps: 0.0383590–0.0383597 at γ = 10⁻², and 0.0032667–0.0032788 at γ = 10⁻⁴ (0.4%). The inside law τG_α∂²_αW is right: it is the Itô term of a drift-free α.
- What fails: the inference that the masses "settle" at the drift's zero, i.e. that the layer "stops" the descent. A zero of an Itô drift with a strongly varying diffusion is not where the mass accumulates (ID-1).

**(iv) First order: correct as a formula. Stated validity: somewhat generous.** Details in ID-5.

**Numerics.**
- Part A agrees, but the text quotes stale numbers (ID-3).
- Part E is confirmed independently. There are exactly 7 triangles with nonzero winding (index sum +1), the 7 zeros match to 4 digits, and the normal drift points inward on all three edges (min b_k = 5.05·10⁻⁸ > 0 at pad 0.002), so Poincaré–Hopf holds.
- Part F is fine.
- Part H overstates (ID-1, ID-4).
- The md regenerates byte-identically.

## Issues

**ID-1 (major): the layer does not stop the masses; the stationary law is not concentrated at the drift's zero.**
- *Location:* Prop. D‴.37 title and (iii) ("stop in a boundary layer", "the masses settle there"); §0 row ("stops it"); §8 bullet ("A boundary layer … stops it"); RESULTS (D) ("stops the descent"); PAPERS paper 13 ("until a boundary layer … stops them"); README ("stopped by a layer"); the [S65] audit ("the boundary layer … that stops it").
- *Analysis:* inside the layer the Proposition's own inside law says α is drift-free (free diffusion at 2τG_α), with W = A/a(α). The zero-flux density in W is then p ∝ |W_α|/D ∝ 1/|W_α| ∝ 1/(W − W_min). That is not integrable, so mass piles up at the saturation (null recurrence in α).
- *Evidence along W_2 = W_3 (γ = 10⁻²):* 2b_1 − ∂_{W_1}D_11 is negative at every one of the 30 points from 0.017δ* to 1.39δ* (−5.7·10⁻⁸ to −6.3·10⁻⁷), including at the zero of b_1. (2b − D′)/D ≈ −1/x near the saturation, consistent with p ∝ 1/x.
- *Evidence in the full-loop records (H_rec, 12 seeds):*
  - Mean over seeds of the largest |α| per record: 4.5, 8.1, 9.3, 11.7 in the four quarters of the run. It is still growing.
  - Fraction of records with |α_k| > 10, by quarter:

    | class | Q1 | Q2 | Q3 | Q4 |
    |---|---|---|---|---|
    | 1 | 0.02 | 0.21 | 0.29 | 0.26 |
    | 3 | 0.00 | 0.05 | 0.12 | 0.28 |

  - Time fraction within 0.004 of the saturation: 0.12 / 0.004 / 0.02 in the first half, against 0.28 / 0.03 / 0.22 in the second half.
  - |α| maxima per seed reach 13–27 in 10 of 12 seeds, against α° ≈ 5.3 at the attracting nodes.
  - W even drops below W_min, i.e. off M: minimum −0.0076, and in the last quarter 11–12% of records for classes 1 and 3.

  This is D‴.36's "unbounded excursion" heuristic, confirmed, not stopped.
- *Fix:*
  - Restate (iii): the layer changes which variable is slow (in α the motion is free), and the zero of b is not a point of accumulation. The relevant quantity is the zero-flux condition 2b = ∂D (multi-dimensional: b − ½∇·D).
  - Re-label E's "attracting nodes" as zeros of the mean drift (the mean path's attractors), not states the masses settle in.
  - Drop "stops it/them" everywhere listed above, and keep "the stationary law of the masses" open, or state its non-normalisability in the 1-D reduction.

**ID-2 (major, literature/novelty): Li–Wang–Arora already have the log-determinant drift for isotropic noise.**
- *Location:* [S64]/[S65] and the D‴.37 audit item; the "New:" sentence of D‴.37; the §0 row.
- *Evidence:* LWA (ICLR 2022), Corollary 5.1: for Σ = I the limiting diffusion on the manifold of minimisers is a Brownian motion on Γ plus the drift −½∇_Γ log|∇²L|₊ (the log pseudo-determinant). This is the "flat minima by log-det" statement in exactly the Katzenberger setting. It is cited only for the framework. Verified by web search on the abstract/summary; I did not open the proof.
- *Also missing (from memory, unverified):*
  - The stiff-vs-rigid constraint correction in molecular dynamics: Fixman 1974; Go–Scheraga 1976; Bornemann–Schütte 1997 on homogenisation with a strong constraining potential. There the soft variables feel ½k_BT log det of the stiff Hessian, the same mechanism as (ii).
  - Label-noise SGD giving tr ∇²L, a contrast showing the selected functional depends on the noise's shape: Damian–Ma–Lee 2021; HaoChen–Wei–Lee–Ma 2021, "Shape matters"; Blanc et al. (already in [S64]).
  - Noise-induced drift in replicator / evolutionary dynamics: Fudenberg–Harris 1992; Foster–Young 1990; Imhof 2005.
  - Temperature-driven weight selection in mixtures: deterministic annealing, Rose 1998.
- *Fix:* cite LWA Cor. 5.1 explicitly in D‴.37. Narrow "New:" to what LWA does not cover: degenerate anisotropic noise (masses unnoised) in a Fisher–Rao/Euclidean product metric; the time-scale-separated γ → 0 limit giving the log-det of the position block only (not the full pseudo-determinant), with the closed form of ID-6; the finite-γ field; the saturation mechanism; the leak's bias. Add the Fixman/Bornemann–Schütte, label-noise and evolutionary-game precedents.

**ID-3 (minor): Part A's numbers are stale.**
- *Location:* D‴.37 Numerics (A) ("to 9.5·10⁻⁴ at W = 1/3 and 1.5·10⁻³ at (0.30, 0.35, 0.35), the finite differences' O(h²)") and the §7 row (same numbers).
- *Evidence:* the md and the raw file give 2.2·10⁻⁵ and 1.1·10⁻⁴ (raw: 2.197e-5, 1.054e-4). My Richardson check shows the h = 0.02 brute force is accurate to about 10⁻⁴, so "O(h²)-limited at 10⁻³" no longer describes the data. The §0 row and RESULTS say "to 10⁻³", which is also stale.
- *Fix:* quote 2.2·10⁻⁵ and 1.1·10⁻⁴, and drop the O(h²) explanation.

**ID-4 (minor): Part H's "means within the 12 seeds' errors at all three times" is an overclaim; the SDE cannot test the saturation.**
- *Location:* D‴.37 (H); RESULTS (H) ("the limit SDE reproducing the means").
- *Evidence:*
  - (full loop − SDE) in standard errors: (1.55, −1.00, −0.82) at 10⁵; (0.36, −0.15, −0.31) at 3·10⁵; (0.12, 1.05, −2.23) at 10⁶.
  - Spreads per component at 10⁶: 0.145 against 0.100 for W_1, and 0.081 against 0.117 for W_3. With 12 seeds the relative standard error of a standard deviation is about 21%, so these are about 2σ and 1.5σ apart. The range comparison "0.10–0.13 against 0.08–0.15" hides this.
  - The SDE is confined to W ≥ W_min + 0.004 with refused steps, which makes the edge a sticky boundary. Its grid spacing is about 0.016, so the layer (δ* = 0.058, zero at 0.038) spans about 3–4 grid cells while D_11 changes by an order of magnitude. By ID-1 the true dynamics put a non-normalisable density exactly where the SDE is cut off. The SDE's near-saturation fractions (0.30, 0.18, 0.26 against 0.50, 0.25, 0.33) are therefore artefacts of confinement.
  - ds = 50 itself is fine: the drift step is about 10⁻⁵, the noise step about 1.2·10⁻³, and rate·ds is about 10⁻³.
- *Fix:* report the deviations (one component at 2.2 SE at 10⁶). State that the comparison tests only the interior. Either run the limit SDE in (α-coordinate) variables that resolve the layer, or drop the near-saturation comparison.

**ID-5 (minor): Part G's "within 1–2 standard errors" hides a consistent sign, and the validity window is narrower than "ε ≳ 3τ".**
- *Location:* D‴.37 (iv) and (G); §0/§7 rows; §8 bullet ("for ε ≳ 3τ").
- *Evidence:*
  - (measured − κ⁻¹b)/SE for W_2 is +2.83, +2.06, +1.68 at ε = 10⁻⁴, 3·10⁻⁴, 10⁻³. Measured/predicted for W_2 is 2.0, 1.7, 1.6, always above the prediction, and still 1.6 at ε = 10τ, where the O(τ/ε) correction should be about 10%.
  - At ε = 3τ the W_1 ratio is 1.40.
  - The three ε use the same seeds (200–215), so the three tests are correlated, not independent confirmations.
  - The reduction also needs ε small against the slowest normal rate (3.06·10⁻⁴), so ε = 10⁻³ is outside it (D‴.36 says so). The honest window at τ = 10⁻⁴ is essentially empty in theory, though the numerics are reasonable.
- *Fix:* state the W_2 excess and its sign. Write the condition as τ ≪ ε ≪ slowest normal rate, and say the instance only approximately satisfies it. Use independent seeds per ε, or say they are shared.

**ID-6 (minor, a strengthening): state the closed form of ℓ.**
- *Location:* (ii)/(iii) and the proof of (iii).
- *Evidence:* ℓ(β) = 2Σ_k log(W_k a′_k) + const = 2Σ_k log(W_k − W_min,k²/W_k) + const, exactly for one particle per class (constant matches to 10⁻⁹). Hence:
  - the slow law is explicit and depends on the bumps only through W_min;
  - ℓ is increasing in each W_k, so its only critical point on the simplex is a maximum (the repeller);
  - "converge toward its minimisers" is empty: ℓ is unbounded below and has no minimiser; the flow simply goes to the boundary.
- *Fix:* add the formula. Replace "minimisers" by "the saturations (ℓ → −∞; no interior minimiser)". The finite-difference `log_det_law` (errors up to 4·10⁻⁴ near the layer) can then be replaced by the exact formula in Part C and in the test.

**ID-7 (minor): the metastability wording.**
- *Location:* (ii) proof ("a reversible diffusion with the Gibbs law", "averaging … over it").
- *Detail:* the global Gibbs law exp(−J/τ) over all positions includes other basins, and for W_k near W_min it includes the flat saturated region (α → ∞, an infinite-measure tail of the Lebesgue reference). The averaging that is actually used is over the local quasi-stationary measure of the exact-fit well.
- *Fix:* say "the quasi-stationary law of the well at the exact representation (metastable on the time scale 1/(γτ))". Add that Laplace is valid only where W_k − W_min ≫ the layer, which ties into ID-1.

**ID-8 (minor): the tests miss the Itô term.**
- *Location:* tests/test_gauge_drift.py.
- *Evidence:* mutants M2 (drop the Itô term) and M9 (double it) pass all 6 tests. The toy tests use linear observables, and on the instance the Itô term is O(γ²) at γ = 10⁻⁵. Yet the Itô term is 4–5% of the drift at γ = 10⁻², the regime of E–H.

  | Mutant | Tests failed (of 6) |
  |---|---|
  | Drop the tangent–tangent block | 2 |
  | Factor 1 instead of ½ | 4 |
  | 2λ_j instead of λ_i+λ_j | 1 (only the integrated-map test) |
  | Sign flip | 3 |
  | Drop P | 2 |
  | P instead of I−P in the tangent block | 3 |
  | Count off-diagonal terms once | 1 |
  | Drop the Itô term | 0 |
  | Double the Itô term | 0 |

- *Fix:* add a toy test with a nonlinear observable (e.g. coords = y₀² on the curved-manifold flow, where the Itô term is ½·2Σ_uu). Add an instance test at γ = 10⁻² against a stored brute-force value; Part A's 2.2·10⁻⁵ agreement would catch a 5% Itô error.

**ID-9 (minor): Katzenberger's regime fails within roughly 0.1δ* of the saturation.**
- *Location:* (i)/(iii) and Status.
- *Evidence:* the full loop's masses fluctuate transversally by about 0.005, reaching W − W_min = −0.0076, below the representable set where no exact fit exists. That is comparable to the distances at which the inside law is checked (0.002–0.0028).
- *Fix:* add to the Status that the W-coordinate description of the limit loses meaning within about 0.1δ* (≈ 0.005) of the saturation. There the description must be in α, with the mass as a fast variable, which is what the inside law says.

**ID-10 (cosmetic): numbers.**
- (B) says "5.35·10⁻⁷"; the md and raw file give 5.343·10⁻⁷, so 5.34.
- RESULTS (B) "the centres' and widths' 10⁴ times smaller": the ratio is about 2.4·10³ at 1/3 and 7·10² at (0.30, 0.35), so "10³ times".
- (D) "0.57δ* as γ → 0": the sequence 0.664, 0.609, 0.585, 0.572, 0.566 is still decreasing (differences roughly halve), so the limit is about 0.56. Say "≈ 0.56–0.57".
- All other numbers in the text, rows, bullet and RESULTS agree with the md: B parts, C, D roots and ratios, E zeros, rates and |a|/a_max = 0.85–0.88, F, G values, H fractions and maxima.

## Summary

| Part | Verdict |
|---|---|
| (i) | Correct. Independently re-derived and re-implemented without the eigenbasis (agreement ≤ 6·10⁻⁶ in the interior). Brute force at two new points, one inside the layer, agrees to 6·10⁻⁵. |
| (ii) | Correct. Convergence to the law is O(γ). D = O(γ²τ) confirmed. Split invariance holds (also for Gaussian kernels on this instance). The exact closed form of ℓ is missing. |
| (iii) | The divergence, both rates, δ* (crossing at 0.98δ*), the zero's robustness to the steps, and the inside law are all right. "The layer stops the masses / they settle at the zero" is wrong: the zero-flux density ∝ 1/(W − W_min) is non-normalisable, and the full loop's α excursions keep growing (ID-1). |
| (iv) | The formula is fine. The stated validity and the "within 1–2 SE" summary are generous (ID-5). |
| Numerics | E confirmed independently (7 zeros, index sum 1, inward boundary). The md regenerates byte-identically. A's quoted numbers are stale. H overclaims. |
