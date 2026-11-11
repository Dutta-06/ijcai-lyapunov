# Referee report D3q (round 21) — PARADIGM §5d: Proposition D‴.30 (statement (i)–(iii), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.30 |", the §8 bullet "Done (Proposition D‴.30)", `results/RESULTS_PARADIGM.md` (Proposition D‴.30), the LITERATURE audit bullet and [S53], `src/core/measurement_noise.py` (the D‴.30 part), `experiments/exp_measurement_noise.py` → `results/chaos/measurement_noise.md` (+ raw), `tests/test_measurement_noise.py`, PAPERS.md paper 12

Branch `claude/festive-bardeen-eopej2` at 637a051, HEAD throughout the review. The concurrent session's uncommitted edits (18 files, among
them PARADIGM.md, RESULTS_PARADIGM.md and LITERATURE_TRANSPORT.md) touch D‴.26–D‴.29, not the D‴.30 text. Every script imports the
637a051 tree, extracted with `git archive` into `ref_D3q/tree`.

Line numbers are 637a051's:
- PARADIGM: the statement l. 3396–3430 ((i) l. 3404–3406, (ii) l. 3408–3417, (iii) l. 3419–3427, "Hence" l. 3428–3429); the proof
  l. 3431–3449; the Status l. 3450–3463; the *Numerics* l. 3465–3483; the §0 row l. 134; the §7 row l. 6021; the §8 bullet l. 6285–6291.
- RESULTS: l. 704–720.
- LITERATURE_TRANSPORT: [S53] l. 73; the audit bullet l. 427–434.
- PAPERS.md: paper 12, l. 165–181.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread (`OMP_NUM_THREADS=1`). pytest ran with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the scratch tree
rewrites `measurement_noise.md` byte-identically (md5 `8f63199d…`, equal to `git show 637a051:results/chaos/measurement_noise.md`).

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3q/scripts/`:

| script | what it does |
|---|---|
| `inst.py` (module) | independent instance: analytic features, ∇_ϑΦ̃, the whitened Ĵ; x- and ϑ-derivatives by the complex step (no repository code) |
| `lin.py` (module) | independent linear-periodic machinery: the collapsed linearised loop (white or OU-augmented), fixed-step RK4, a Kronecker discrete Lyapunov solve, periodic mean responses |
| `check_partA.py` | every Part A number; the force and its split into the regressor's part and the E[e] part; the crossovers |
| `check_partB.py` (+ `check_partB_800res.py`, `derived_B.py`) | Part B's exact periodic floors at k_e = 25–800 against the raw record; per-mode variances; RK4 resolution (N, 2N); fits in 1/k_e |
| `check_parity.py` | the floor with the slope a removed from the noise input (the origin of the odd order) |
| `check_force.py` | (iii) exactly: the stationary second-order mean of the linearised periodic loop with OU noise, from its periodic covariance, with every second-order term (A)–(D), at k_e = 25–200 |
| `transient_second_order.py` | the exact second-order transient from the lift (Part C's setting): covariance and mean ODEs together |
| `sim_batch.py` (module), `run_antithetic.py` | independent batched real loop with the experiment's noise model; the 12 stored paths and their negatives (antithetic pairs) |
| `repro_C.py`, `check_lift.py` | Part C's seed-0 K₀-atom run and m = 15 lift with the committed code, bitwise against the raw record; their float64 comparison |
| `recompute_C.py` | Part C's z-scores, shares and fixed/adaptive differences from the raw record; correlations between cases |
| `check_sep.py`, `sim_sep.py` | the separations of a class under noise: the averaged second-order drift S_k, and finite-time Lyapunov exponents along noisy trajectories |
| `check_misc.py` | Q's cross term on references without the time symmetry; the exact OU floor against the white formula; the bias against k_e and B |
| `check_tests.py` | the tests' asserted quantities; mutants of `measurement_noise.py` (patched in memory) |
| `derived.py` | ratios and crossovers from the outputs above |

Also in `ref_D3q/`: `tests.out` (4 passed, 15.67 s), `from_raw.out`, `md5_before.out`, `md5_after.out`. The experiment was not rerun in
full. Two Part C runs were reproduced bitwise with the committed code, and the 24 stored std-0.01 runs (12 per prior) to ≤ 6.2·10⁻¹³ by the
independent simulator.

## Verdicts

| object | verdict |
|---|---|
| **(i) coherence** | **holds, exactly, as invariance.** Re-derived: with frozen classes, every particle of a class sees the same y, the same ε_m and the same prior atom (the class mean at η = ∞, z_k at finite η). The replicated subspace is invariant for every path, and on it f̂ = Σ_k(N_k/m)Φ̃(ϑ_k) with the unweighted law: the weighted K₀-atom loop at w = N/m, at every order in ν. With recomputed nearest-atom cells this holds while no two classes share a nearest atom (near the lift). The contrast with Theorem J′ is fair: Corollary J′.3 measures the dither's 1/m law. **But invariance is not attraction** (D3q-3): the noise shifts the separations' exponents at second order. On the instance class 2's separations stop decaying at std ≈ 0.03. The "0.0" of the lift is a float32 artefact (D3q-3). |
| **(ii) the floor** | **holds**, as a leading-order formula with an absolute O(k_e⁻²) remainder. Re-derived: the linearisation (δė = −k_eδe − Ĵy − (k_e + a)ν, ẏ = Ĵᵀ(δe + ν) − leak), n = (d/dt + k_e)⁻¹(d/dt − a)ν, the slow diffusion (N₀/k_e)Q, the pair-block Lyapunov solution Σ_yy = q(h + η)/(2h(h + η + σ₀)), and the fast part y_fast ≈ Ĵᵀξ with Var ξ = N₀/(2k_e). The "total derivative" of the proof is the fast part itself, same order in k_e as the slow part (0.74% of the floor here). The white-noise formula applies at B = 4k_e: exact OU/exact white = 0.9942. Two problems. The relative error is O(k_e⁻¹), not O(k_e⁻²): the slope a in the noise input breaks the time-reversal parity, and the ratio crosses 1 between k_e = 400 and 800 (D3q-1). The per-mode statements are not checked: the weakest mode's fixed-prior variance is 15.5% below the formula at k_e = 25 (D3q-2). |
| **(iii) the errors in variables** | **holds for the averaged system; the force is right.** Re-derived: the formulas for E[ε_m²], E[νe], Var(e) and E[e]; the stationary points (H + σ₀I)⁻¹ and H⁻¹, the latter for every η > 0 (the exact periodic loop gives the same 2652.93 at η = σ₀ and ∞). The quasi-static moments are accurate to ±0.32% (E[ε_m²]) and 0.22% (E[e]) along the orbit. The force is confirmed in the real loop mode by mode by antithetic pairs, mode 2 to ≤ 1.7% (that mode carries 75% of the adaptive bias). Three problems. The averaged response understates the exact second-order stationary mean of the linearised loop by 10.7% (adaptive) and 2.6% (fixed) at k_e = 25 (D3q-2). The derivation omits second-order terms, negligible here (≤ 1.8% per mode, 4·10⁻⁶ of the bias), and needs ∂²_xf̂ = ∂²_xf (D3q-4). "Grows like k_e²" holds only for B ≳ 40k_e, where the local exponent is ≥ 1.93 (D3q-4). |
| "Hence" (l. 3428–3429) | **holds** with its qualifier k_e ≪ B. On the instance (B = 100) the bias's local exponent in k_e is 1.69, 1.47 at k_e = 12.5–50 and below 1 beyond 100 (D3q-4). |
| Proof | **correct in structure** (linearisation, slaving, stochastic averaging, Lyapunov equation, second-order moments). It has an omission (D3q-4) and imprecisions (D3q-7). |
| *Status* | **accurate** on what is proved, except "(the ratios tending to 1 like k_e⁻²)" (D3q-1). The novelty is overstated (D3q-6). "The second-order force linear in the noise's variance" is vacuous as a check (D3q-5). |
| *Numerics* | **every number reproduces** (confirmations 1–4). Readings that do not hold: "the deficit falling like k_e⁻²" (D3q-1); "×24" and "491 … 2397 (×4.88)" quoted as the instance's (D3q-2); "which the transient checks" (D3q-2); "to 0.0" (D3q-3); "by the reference's time symmetry", "falls from 104 … to 5.5", "5–6%" (D3q-7). |
| §0 row, §7 row, §8 bullet, RESULTS, PAPERS 12 | **numbers hold.** Wording: "the deficit O(k_e⁻²)" (D3q-1); "×24", "491 and 2397" (D3q-2); "identical", "(12 runs, two noise levels)" (D3q-3, D3q-7); the §8 "Open: the remedies checked (the reference in the regressor …)" is stale after D‴.31; PAPERS' "(iii) … checked against the exact linearised loop" is not done (D3q-2, D3q-7). |
| LITERATURE bullet, [S53] | **bibliographic details right**; the audit omits the closest precedents (D3q-6). |
| `src/core/measurement_noise.py` (D‴.30 part) | **holds.** Against independent code: H and Q to 3.6·10⁻¹⁶; tr H² and ⟨|Ĵ|⁴⟩ exact; F to 7.2·10⁻¹¹; `periodic_covariance` to 1.6·10⁻¹¹–2.1·10⁻⁹; `fixed_generator` likewise; `transient_covariance` to 1·10⁻⁴ at t = 300. `NoisyRealLoop` (defaults) and `simulate_noisy` agree to ≤ 6.2·10⁻¹³ on 24 stored runs. `averaged_covariance`, `floor_closed_form` and `mean_response` hold (test 1, test 4, the transient recomputation). |
| `experiments/exp_measurement_noise.py` | **holds**; regeneration is byte-identical; single runs are reproducible bitwise. `job_C` stores f̂ in float32, which produces the "0.0" (D3q-3). All four Part C cases run on the same 12 normalised noise paths (D3q-7). |
| `results/chaos/measurement_noise.md` | **numbers hold**; readings as above. |
| `tests/test_measurement_noise.py` | **passes (4/4); the force is untested** (D3q-5). Mutants flipping the force's sign, multiplying it by 10, dropping E[e], or setting E[ε_m²] = Eν² pass. Test 3 passes with the noise removed from the adaptive law and with a wrong OU variance. Test 2 runs at η = ∞ only. |

## Confirmed issues

**D3q-1 (minor; Status l. 3452–3453 "(the ratios tending to 1 like k_e⁻²)"; *Numerics* l. 3473–3475 "the deficit falling like k_e⁻²"; §0
row "(k_e = 25–200, the deficit O(k_e⁻²))"; the md Part B).** The floor's relative error has a first-order term of the opposite sign. It
is hidden on k_e = 25–200 by a large second-order one (`check_partB.py`, `derived_B.py`, `check_partB_800res.py`).
- *Ratio exact/predicted − 1, η = ∞, k_e = 25, 50, 100, 200, 400, 800:* −1.969·10⁻², −5.092·10⁻³, −1.163·10⁻³, −2.210·10⁻⁴,
  −1.946·10⁻⁵, **+1.309·10⁻⁵**. Successive ratios: 3.87, 4.38, 5.26, 11.4, then a sign change.
- *η = 0:* −2.375·10⁻², …, −3.098·10⁻⁵, **+1.162·10⁻⁵**.
- *Resolution:* the k_e = 800 value is resolved (1.30928·10⁻⁵ and 1.30935·10⁻⁵ at 2¹⁷ and 2¹⁸ RK4 steps per period).
- *Fit:* ratio − 1 = c₁/k_e + c₂/k_e² + … gives c₁ = +0.029 (η = ∞) and +0.031 (η = 0), with c₂ = −14.7 and −17.5 (max residual ≤ 1.7·10⁻⁷).
  The first-order term dominates beyond k_e ≈ 500.
- *Per mode* the first-order corrections are larger and halve per doubling. Mode 4 (η = ∞, with the fast part included): 1.0226, 1.0106,
  1.0051, 1.0025, 1.0012, 1.0006.

*Mechanism* (`check_parity.py`). The noise input of the error channel, −(k_e + a(t))ν, carries the plant's slope. a is a rate that does
not change sign under time reversal: the parity breaker that Corollary D‴.26's proof names for a in the drift. With the noise input −k_eν
(a removed, generator unchanged), the deviation is even in 1/k_e: −2.191·10⁻², −5.976·10⁻³, …, −2.411·10⁻⁵, successive ratios 3.67,
3.90, 3.98, 3.99, 4.00, c₁ ≈ 2·10⁻⁴. The statement's "+ O(k_e⁻²)" (absolute) is right, and so is the Status's "formal". Only the reading
of the numerics is wrong: the lesson of D3m-1 (the observed order on a limited range).

*Fix:*
- Status, *Numerics*, §0: "the deficit 2.4%, 0.61%, 0.14%, 0.028% at k_e = 25–200 (η = 0), dominated there by an O(k_e⁻²) term
  (−17.5/k_e² at η = 0, −14.7/k_e² at η = ∞). The relative error is O(1/k_e) (+0.031/k_e and +0.029/k_e): the slope a in the noise input
  −(k_e + a)ν breaks the time-reversal parity, and the ratio crosses 1 near k_e = 500."
- Add k_e = 400 and 800 to Part B, or the a-free noise input as a control.

**D3q-2 (minor; *Numerics* l. 3470 "the weakest mode's parameter variance grows from 0.0865 to 2.09, ×24", l. 3471–3473 "the
function-level stationary bias 491 (fixed prior) and 2397 (adaptive, ×4.88) … equals the floor of the parameters at the noise std 0.0026
and 0.0012", l. 3480–3482 "so the stationary amplification rests on the averaged equations, which the transient checks"; the §0 and §7
rows, RESULTS l. 711–712, the md Part A; PAPERS 12 "(ii)–(iii) … checked against the exact linearised loop").** These are the averaged
system's values (k_e → ∞), quoted as the instance's at k_e = 25. Neither the per-mode variances nor the stationary second-order mean was
checked against the exact linearised loop. The exact values differ by more than the floor's 2.4%.
- *Per-mode variances* (`check_partB.py`). At k_e = 25 the weakest mode's exact variance is 0.0731 under the fixed prior (0.845 of
  0.0865: the noise's attenuation by the error filter, O(k_e⁻²), 4.5% at k_e = 50) and 2.045 at η = ∞ (0.980 of 2.0875). The ratio is
  **×28.0, not ×24.1**; it is 25.1, 24.4, 24.2 at k_e = 50, 100, 200. The per-mode formula also omits the fast part (N₀/2k_e)h_j, 1.3% of
  mode 9's variance.
- *The stationary errors-in-variables bias* (`check_force.py`). This is the exact periodic second-order mean of the linearised loop
  with OU noise, all second-order terms, resolved to 6 digits.

  | k_e (B = 100) | adaptive: averaged → exact | excess | fixed: averaged → exact | excess |
  |---|---|---|---|---|
  | 25 | 2396.95 → **2652.93** | +10.7% | 491.06 → **503.71** | +2.6% |
  | 50 | 6648.5 → 6824.7 | +2.65% | 1365.9 → 1374.6 | +0.64% |
  | 100 | 14949 → 15048 | +0.66% | 3075.6 → 3080.6 | +0.16% |
  | 200 | 26568 → 26612 | +0.17% | 5470.1 → 5472.3 | +0.04% |

  The excess falls by 4.0 per doubling, so it is O(k_e⁻²). The fixed/adaptive ratio is **5.27**, not 4.88. Mode 2 is 6.8% above the
  averaged response. The source is the response, not the force: the quasi-static moments are accurate to 0.3%. The weak modes'
  effective stiffness is the filtered Gram H − k_e⁻²⟨Ĵ̇ᵀĴ̇⟩, as in D3m-1. With the exact biases and the exact OU floors, the crossover
  noise levels are std 0.00248 and 0.00114, against 0.00255 and 0.00121.
- *"Which the transient checks."* The transient checks the force; it cannot check the stationary stiffness. In the exact second-order
  transient (`transient_second_order.py`) and in the real loop (confirmation 7), the dominant mode's deviation follows the averaged
  prediction within 1.8% (mode 2) over t ≤ 300. A mode with time constant k_e/h_j ≫ 300 is in its linear growth ȳ_j ≈ k_eF_js
  there, which does not involve h_j.

*Fix:* quote "(averaged; the linearised loop at k_e = 25: ×28.0)", "491 and 2397 (averaged), 504 and 2653 (×5.27) in the exact linearised
loop at k_e = 25, the excess O(k_e⁻²)". Replace "which the transient checks" by "the transient checks the force (the weak modes grow
linearly over 300 time units); the stationary response is checked against the exact linearised loop (+10.7% at k_e = 25)". In PAPERS 12
write "(ii) checked against the exact linearised loop, (iii)'s force in the real loop's transient". Add the per-mode check to Part B,
and the exact stationary mean to the experiment.

**D3q-3 (minor; (i) l. 3404–3406 "a replicated configuration stays replicated under every path of ν … Measurement noise does not average
out with m"; the module docstring "the separations are not excited"; *Numerics* l. 3482–3483 "gives the same f̂ on the grid to 0.0 (the
error's difference 6·10⁻¹⁹): the noise is coherent, as (i) states"; RESULTS l. 716 and the §0 row "identical").**

*Invariance is not attraction.* A separation of class k leaves f̂ unchanged at first order. It obeys δϑ̇_i = [Gε_m(t)∇²_ϑΦ̃(y(t); ϑ_k(t))
− σ₀/k_e]δϑ_i along the noisy class trajectory. Its averaged second-order drift (whitened per particle, `check_sep.py`) is
S_k = ⟨E[ε_m²]G^{1/2}∂_x∇²_ϑΦ̃G^{1/2} + E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩, almost all of it the regressor's part. Per unit Eν², its largest
eigenvalues are 0.0148, 0.0101 and **0.298** in classes 0, 1, 2. Class 2's noiseless margin σ₀/k_e = 2.83·10⁻⁴ is used up at
Eν² = 9.5·10⁻⁴ (std 0.031). The finite-time Lyapunov exponents of class 2's separations along noisy real-loop trajectories (`sim_sep.py`;
4 paths of 300 time units; noiseless −2.830·10⁻⁴) confirm it:

| std | top exponent (measured) | predicted −σ₀/k_e + 0.298Eν² |
|---|---|---|
| 0.01 | −2.490·10⁻⁴ ± 0.031·10⁻⁴ | −2.53·10⁻⁴ |
| 0.03 | **+1.52·10⁻⁵ ± 0.70·10⁻⁵** (all four paths positive) | −1.5·10⁻⁵ |
| 0.05 | **+4.93·10⁻⁴ ± 0.08·10⁻⁴** | +4.62·10⁻⁴ |

At the md's noise levels the margin shrinks by 0.9% (std 0.003) and 10.5% (std 0.01). Above std ≈ 0.03 (1.5% of the reference's
amplitude), the class-2 separations of an m-system started near, rather than at, the replicated configuration grow (e-folding time
≈ 2000 time units at std 0.05). The system then leaves the weighted K₀-atom loop, and D‴.26(a)'s stability of perfect tracking is lost
under noise. The critical variance is (σ₀/k_e)/λ_max(S_k): the
design's small σ₀ makes it low.

*The numerical demonstration cannot see this* (`repro_C.py`, `check_lift.py`). The within-class spread of the lift is exactly 0.0
because every particle's update is the same floating-point arithmetic on the same inputs. That holds whether or not the replicated state
is stable. The lift against the K₀-atom run differs in float64 by 1.1·10⁻¹⁵ in f̂ on the grid (different arithmetic: a mean of 15 against
a weighted sum of 3). The "0.0" of the md is the float32 cast in `job_C` (`fh=….astype(np.float32)`; spacing 6·10⁻⁸ at max|f̂|). The §7
row's "identical to rounding" is the right wording.

*Fix:* add to (i): "Invariance only: under noise the separations' decay rates are σ₀/k_e − Eν²λ(S_k) + …, S_k as above. On the instance the
replicated configuration stops attracting at std ≈ 0.03 (class 2), so the identification with the K₀-atom loop for starts off the
replicated subspace needs the noise below that level." Print S_k's eigenvalues and the critical std in Part A. Store f̂ in float64 and
write "to rounding (1.1·10⁻¹⁵)". Better, run the lift from a perturbed start at two noise levels and print the separations' decay.

**D3q-4 (minor; (iii) l. 3419–3427, proof l. 3444–3449).** The derivation and the claims of (iii) have three gaps.
- *Omitted second-order terms.* The second-order mean of the law has the further term ∇²_ϑΦ̃·E[ε_mδϑ]. The error equation has
  −Ĵ′E[ε_my] and −½Σ_k tr(Ĥ_kΣ_θ,kk), the parameter fluctuations rectified by the network's curvature. The proof's "the second-order part
  of E[ε_m∇_ϑΦ̃(y; ϑ)] at the lift is ∂_x∇_ϑΦ̃·E[ε_m²] + ∇_ϑΦ̃·E[e]" drops them silently. They are O(N₀) = O(Eν²/B), not of higher order
  in ν. At η = ∞ they involve the large weak-mode variances. On the instance they are negligible (`check_force.py`, the stationary
  periodic mean at k_e = 25, B = 100): ≤ 1.8% of a mode's deviation (mode 1), ≤ 0.61% elsewhere, and 4·10⁻⁶ of the bias (2652.92 →
  2652.93). They grow relative to F like 1/B for coloured noise.
- *A hidden hypothesis.* E[e] = −(∂²_xf̂/(2k_e))Eν²(B − k_e − 2a)/(k_e + B) uses one second derivative for f(x) and f̂(x + ν). Without
  ∂²_xf̂ = ∂²_xf on the orbit the error equation keeps ½(∂²_xf − ∂²_xf̂)Var(e), of the same order. The hypothesis states only
  ∂_xf̂ = ∂_xf. The exact representation (an identity of functions) supplies both; the module's `mean_error` docstring says so.
- *"Grows like k_e²(Eν²)² while … k_e ≪ B."* This is right asymptotically but is not realised on the instance (`check_misc.py` (c)). At
  B = 100 the local exponents of the averaged adaptive bias in k_e are 1.69, 1.47, 1.17, 0.83, 0.53, 0.30, 0.17 (k_e = 12.5 → 1600).
  They are 1.96–1.75 at B = 1000 (k_e ≤ 200) and 1.99–1.94 at B = 10⁴. "More slowly beyond" is a saturation: the bias tends to a
  constant for k_e ≫ B.

*Fix:* in the proof, "plus ∇²_ϑΦ̃·E[ε_mδϑ] in the law and −∂_xDf̂[E[ε_mδϑ]] − ½E[D²_ϑf̂[δϑ, δϑ]] in the error equation, of order N₀; neglected
(on the instance ≤ 1.8% of any mode's response)". Add "∂²_xf̂ = ∂²_xf on the orbit" to (iii)'s hypothesis. Add "on the instance
(B = 4k_e) the exponent is 1.5–1.7" and "saturates" after "more slowly beyond".

**D3q-5 (minor; `tests/test_measurement_noise.py`; Status l. 3459–3463; RESULTS l. 718–720).** The second-order force, the object of
(iii), is untested.

| asserted quantity | measured | tolerance |
|---|---|---|
| test 1: closed form vs Lyapunov | 6.0·10⁻¹⁵ | 10⁻¹⁰ |
| test 1: trace identity at η = ∞ | 0.0 | 10⁻¹² |
| test 2: orbit and point floors at k_e = 100, η = ∞ | −1.17·10⁻³, −1.09·10⁻³ | 5·10⁻³ |
| test 2: periodicity | 6.1·10⁻¹⁴ | 10⁻⁶ |
| test 4: max\|F₁\| | 2.3·10⁻⁵ | > 10⁻⁷ |

Mutants (`check_tests.py`, patched in memory):
- *Test 4 passes* with the force's sign flipped, the force ×10, E[e] dropped, or E[ε_m²] = Eν² with no E[e]. "The force linear in the
  noise's variance" holds by construction, since `eiv_force` multiplies everything by `var`. Its first part (`mean_response`) does
  catch a finite-η response (H + σ₀I)⁻¹.
- *The Monte Carlo is no substitute.* The function-level error is even in F, so neither the tests nor Part C's E‖f̂ − f‖² can detect
  the force's sign. The antithetic check of confirmation 7 can.
- *Test 3 passes* with the noise removed from the adaptive law (the true error in the law): coherence holds anyway, and the plant's
  jitter still moves the parameters. It also passes with a wrong OU variance in `simulate_noisy`. It runs at η = ∞, where the prior's
  state is never executed.
- *Test 2 catches* Q without the slope term (+8.2%), the missing fast parts (+0.63%, +3.3%), tr H² in place of ⟨|Ĵ|⁴⟩ (+2.5%), Q ×1.004
  (−0.51%) and the noise input −k_eν (−7.8%). It passes Q ×0.995. It runs at η = ∞ only, so `fixed_generator` and the finite-η loop are
  never executed by a test.

*Fix:*
- Test the force against an independent quadrature of ⟨Ĵ′ᵀE[ε_m²] + ĴᵀE[e]⟩ with analytic Ĵ′, or against the even part of an
  antithetic pair in the real loop over ~25 time units (mode 2 to a few per cent).
- Assert the OU's stationary variance in `simulate_noisy`, and that the law uses the measurement (a run with ν in the controller only
  must differ).
- Add η = 0 and η = σ₀ to test 2.
- Remove "the force linear in the noise's variance" from the Status's Checked list.

**D3q-6 (minor; Status l. 3454–3459 "New are the coherence … the closed form of the floor, in which the adaptive prior's mobility η
enters the Tikhonov filter factor … the loop's high-pass … and the identification of the adaptive prior's robustness loss under noise as
the H⁻¹ response"; the §0 row's last column; the audit bullet l. 427–434; PAPERS 12 "Novelty").** The audit omits the closest
precedents.
- *"Whatever the spectrum."* The η = ∞ floor (N₀/2k_e)tr Q comes from the Lyapunov identity tr(HΣ) = ½tr D for HΣ + ΣH = D. This is the
  mechanism of the LMS misadjustment, μσ_v²tr R/2, which depends on the input's eigenvalues only through tr R (Widrow, McCool,
  Larimore, Johnson, Proc. IEEE 64(8):1151–1162, 1976).
- *The fixed prior's per-mode factor.* The factor h_j/(h_j + σ₀) in a noise floor is the steady state of leaky LMS (Mayyas–Aboulnasr,
  IEEE TSP 45(4):927–934, 1997).
- *The stationary covariance.* A Lyapunov solution for the stationary covariance of a constant-gain adaptive algorithm is its Gaussian
  (diffusion) approximation (Benveniste–Métivier–Priouret, *Adaptive Algorithms and Stochastic Approximations*, Springer 1990;
  Kushner–Yin).
- *Coherence.* That a noise common to all particles does not average out in the mean-field limit is the "common noise" of mean-field
  theory: the conditional McKean–Vlasov limit (Carmona–Delarue–Lacker, Ann. Probab. 44(6):3740–3803, 2016).
- *The H⁻¹ response.* The amplification of a persistent disturbance in weakly excited directions is the classical lack-of-PE drift
  (Anderson, Automatica 21(3):247–258, 1985; Ioannou–Kokotović, [3], [4]).

What remains specific:
- the interpolation (h_j + η)/(h_j + η + σ₀), an elementary 2×2 Lyapunov computation;
- the noise's matrix Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩, which arises because the adaptive law sees the sensitivity-filtered measurement;
- the coherence and the H⁻¹ amplification placed in the transport-leak setting.

The bibliographic details above were checked by a web search. The statements of content are from memory; no full text was fetched.

*Fix:* cite the five precedents in the audit bullet and the Status. Restrict "new" to the η-interpolation, the derivative-Gram form of
the driving matrix, and their combination with the coherence of the K₀-atom identification (with D3q-3's noise threshold).

**D3q-7 (cosmetic; wording and consistency).**
- *"The cross term vanishes by the reference's time symmetry"* (*Numerics* l. 3467–3468). The whole cross matrix
  ⟨a d(ĴᵀĴ)/dt⟩ = (1/T₀)∮a(x)(ĴᵀĴ)′(x)dx vanishes for every periodic reference in one dimension: 2.6·10⁻¹⁶ on 2 sin t + 0.8 sin 2t, which
  has no t → π − t symmetry (`check_misc.py` (a)). *Fix:* "vanishes in one dimension (an exact differential)".
- *"Q_jj/h_j falls from 104 on the weakest mode to 5.5."* It is not monotone: 104.4, 52.5, 44.1, 18.0, 14.8, 5.48 (mode 6), 25.5, 18.3,
  15.5. *Fix:* "ranges from 104 (mode 1) to 5.5 (mode 6)".
- *"The adaptive prior costs at most 10%."* It is 10.02%. *Fix:* "10%".
- *"Over 300 time units the fixed and the adaptive prior differ by 5–6%"* (l. 3480). The Monte Carlo difference grows from 0.6% at
  t = 25 to 6.0% (std 0.003) and 6.4% (std 0.01) at t = 300; predicted 5.3% and 6.1% (`recompute_C.py`). *Fix:* "drift apart to 6.0–6.4%
  at t = 300 (predicted 5.3–6.1%)".
- *"12 runs … per case", "(12 runs, two noise levels)"* (§0). All four cases use the same 12 normalised noise paths (`default_rng(1000 +
  seed)` scaled by the std). Per-run errors correlate at 0.983 between the noise levels (t = 100) and 0.99989 between the priors
  (t = 300). The z-sequences at the two levels correlate at 0.80, and successive records within a run at 0.69. The four agreements are one
  test on 12 paths. *Fix:* "12 common noise paths for the four cases".
- *§8 bullet l. 6290.* "*Open:* … the remedies checked (the reference in the regressor, …)" is done by D‴.31 in the same commit.
  *Fix:* drop it, keep "a projected or partially fixed prior on the weak modes".
- *The proof's notation.*
  - With a time-varying a, "((d/dt − a)/(d/dt + k_e))ν" means (d/dt + k_e)⁻¹[ν̇ − aν].
  - "The total derivative stays bounded" and "y_fast ≈ Ĵᵀξ" describe the same term. The exact decomposition is y = Ĵᵀξ + u with
    ξ = (d/dt + k_e)⁻¹ν and u driven by −(Ĵ̇ + aĴ)ᵀξ.
  - The modes decouple in tr(HΣ_yy) because Q's off-diagonal entries in H's eigenbasis enter only Σ_ij, i ≠ j.
  - The formulas (E[ε_m²], ∂_xĴ, ⟨|Ĵ|⁴⟩, the cross term) are for n = 1, which the statement does not say.
- *"Mode j's parameter variance"* is its slow part. The fast jitter adds (N₀/2k_e)h_j.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration, reproduction and tests.** `--from-raw` is byte-identical (md5 `8f63199d…`). With the committed code, Part C's seed-0
   K₀-atom run (std 0.01, η = ∞) and the m = 15 lift reproduce their stored err series bitwise, and their stored float32 f̂ exactly
   (`repro_C.py`, 100 s and 111 s). The tests pass: 4 in 15.67 s.
2. **Part A, independently** (`check_partA.py`; analytic Ĵ, complex-step Ĵ′, 4096 samples).
   - The representation is exact: max|f̂ − f| = 2.2·10⁻¹⁶ and max|f̂′ − f′| = 6.7·10⁻¹⁶ on [−4, 4]. H and Q equal the raw record to
     3.3·10⁻¹⁶ and 3.6·10⁻¹⁶.
   - tr Q = 8.783561 = 8.104920 (⟨Ĵ̇ᵀĴ̇⟩) + 0.678641 (⟨a²ĴᵀĴ⟩); the cross matrix's largest entry is 6.2·10⁻¹⁷.
   - h_j = 3.0587·10⁻⁴ … 0.19521.
   - The floors per unit N₀: 0.159670, 0.162711, 0.168968, 0.175671 (η = 0, σ₀, 10σ₀, ∞); ratio 1.1002. Modes 7–9 carry 88.2% (η = ∞)
     and 92.1% (fixed) of the floor.
   - The weakest mode's variance: 0.086515 → 2.08750, ×24.13. Its rates λ₁^−: 1.5128·10⁻⁴, 2.7796·10⁻⁴, 3.0587·10⁻⁴.
   - tr H² = 0.065712, ⟨|Ĵ|⁴⟩ = 0.299111; k_e/h₁ = 81 735.
   - F equals the raw record to 7.2·10⁻¹¹. The biases are 491.064 and 2396.95 (ratio 4.881). The crossovers are std 0.00255 and 0.00121.
   - The adaptive bias is the regressor's part's: 2393.5 from ⟨Ĵ′ᵀE[ε_m²]⟩ alone, 0.47 from ⟨ĴᵀE[e]⟩ alone. The fixed prior's is the same
     part's (492.0 alone, 0.21). Mode 2 carries 75.0% of the adaptive bias, modes 7 and 9 carry 44% and 31% of the fixed one.
3. **Part B, independently** (`check_partB.py`; RK4 with 8192–131 072 steps per period, a Kronecker discrete Lyapunov solve).
   - All twelve exact floors equal the raw record to 2.4·10⁻¹⁰–2.1·10⁻⁹ (orbit; 1.6·10⁻¹¹ at k_e = 25 with 2¹⁴ steps) and ≤ 4.2·10⁻⁸
     (point; the experiment's 65-point trapezoid).
   - The ratios 0.976248, 0.977432, 0.980309 (k_e = 25) through 0.999722–0.999779 (k_e = 200) reproduce.
   - The resolution is 2.2·10⁻¹⁰ between N and 2N at k_e = 25.
   - The OU floor: with B = 100, the exact OU floor over the exact white floor (N₀ = 2Eν²/B) is 0.99415 (η = 0) and 0.99421 (η = ∞)
     (`check_misc.py` (b)). The white-noise formula applies at B = 4k_e. The fast part under the OU noise is B/(k_e + B) = 0.8 of the
     white one, < 1% of the floor. The sample-and-hold OU of the scheme has the low-frequency density 1.00083N₀.
4. **Part C from the raw record** (`recompute_C.py`).
   - z ranges from −0.980 to +0.991 (std 0.003, η = 0), −0.989 to +0.986 (η = ∞), −0.295 to +1.753 and −0.277 to +1.750 (std 0.01).
   - The mean part's shares at t = 300 are 0.4196, 0.4235, 0.8893, 0.8909. The Monte Carlo is 8.77 and 8.91 times the covariance part at
     t = 300.
   - My exact transient covariance and the proposition's mean part reproduce the md's prediction columns: 3.0004·10⁻³ against 3.0001·10⁻³
     per unit variance at t = 300; 244.90 and 230.33 per unit variance squared.
5. **(ii)'s derivation**, line by line, including the 2×2 Lyapunov solution (Σ_ζζ = Σ_yζ = ηΣ_yy/(h + η)) and the fast part's variance
   N₀/(2k_e). The closed form is the Lyapunov solution for random H, Q to 6.0·10⁻¹⁵ (test 1). The trace identity at η = ∞ is exact.
6. **(iii)'s formulas and the quasi-static assumption** (`check_force.py`).
   - E[ε_m²], E[νe], Var(e) and E[e] re-derived (spectral integrals).
   - Along the orbit the exact periodic E[ε_m²] lies within [0.99687, 1.00322] of the quasi-static formula (mean 1.000157), and the
     rectified E[e] within 2.2·10⁻³.
   - The stationary point is η-independent for η > 0 in the exact linearised loop too: 2652.93 at η = σ₀ and at η = ∞.
   - The periodic mean is resolved to 6 digits (N = 8192 against 16 384).
7. **The force in the real loop, mode by mode** (`sim_batch.py`, `run_antithetic.py`).
   - *Cross-validation.* The independent simulator, with the experiment's noise model and seeds, reproduces all 12 stored runs at std 0.01
     to ≤ 6.2·10⁻¹³ (η = ∞) and 3.0·10⁻¹³ (η = 0) relative in E‖f̂ − f‖².
   - *Antithetic pairs.* With the negated paths, the even part ½(y(ν) + y(−ν)) isolates the second-order mean. Over 12 pairs it equals the
     exact second-order transient (`transient_second_order.py`) in mode 2 to 0.983, 1.000, 1.001, 1.008 at t = 25, 100, 200, 300
     (s.e. 2.0, 0.9, 0.7, 0.5%). Modes 3, 6, 7, 9 agree within 0.2–3.8% (s.e. 0.6–2.9%) at t ≥ 100; the other modes are at the noise.
   - *Sign and size.* The real loop's mean deviation thus has the proposition's sign and size in the modes that carry both biases.
   - *Function level, t = 300 (η = ∞).* The even part's size is 2.520·10⁻⁶, the exact second-order mean 2.497·10⁻⁶, the proposition
     2.449·10⁻⁶. For η = 0: 2.370, 2.348, 2.303·10⁻⁶.
   - *Pair-averaged E‖f̂ − f‖².* 2.758·10⁻⁶ ± 0.069·10⁻⁶ against the prediction 2.749·10⁻⁶ (η = ∞); 2.595·10⁻⁶ ± 0.065·10⁻⁶ against
     2.590·10⁻⁶ (η = 0). The pairing cuts the standard error by 3.9 against the 12 single runs.
8. **(i)**: the lift's within-class spread is 0.0 and its class means equal the K₀-atom parameters to 4.4·10⁻¹⁶, over 300 time units of
   std-0.01 noise (`check_lift.py`). At std 0.01 the separations' exponents remain negative in every class (−2.83, −2.74, −2.49·10⁻⁴),
   as the averaged drift S_k predicts: class 2 within 1.6%, class 0 within 0.4%, class 1 within its standard error (`sim_sep.py`).

## Not checked

- The full experiment (the 48 Part C runs) was not rerun: 2 runs were reproduced bitwise, and 24 to 6·10⁻¹³ by independent code.
- The stochastic averaging was not proved, and the O(1/k_e) coefficient of the floor (D3q-1) was not derived analytically. Its origin
  was identified by the control experiment of `check_parity.py`.
- The separations' instability above std ≈ 0.03 (D3q-3) was established for the linearised separation equation along noisy
  trajectories (4 paths of 300 time units per level). The nonlinear splitting of a perturbed lift was not simulated.
- The real loop's stationary biases (~10⁵ time units) were not simulated; the exact linearised loop's stationary second-order mean
  was computed instead.
- Only the instance (n = 1, w = q, x_m = 2 sin t) was examined. The weighted path at w ≠ q was not checked, and on a non-symmetric
  reference only Q's cross term was computed, not the floor or the force.
- `sim_sep_005.out` is an earlier run whose top exponent was read from an unconverged QR column. It is superseded by `sim_sep_0.05.out`
  (the singular values of the accumulated R-product).
- The literature statements of D3q-6 are from memory; only the bibliographic details were checked by a web search.
- Round 20 is not in `review/REVIEW_LOG.md` at 637a051 (it ends at round 19). Proposition D‴.31, README's sentence on D‴.30, and the
  concurrent session's uncommitted edits are outside the material.
