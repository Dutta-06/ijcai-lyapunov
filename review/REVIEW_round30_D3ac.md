# Referee report D3ac (round 30): the revised Propositions D‴.41 and D‴.42

Object: Propositions D‴.41 and D‴.42 as revised after round 29 (scratchpad `label/D41_r29.md`, `multi/D42_r29.md` and their pieces,
`label/paper15_r29.md`), `src/core/label_noise.py` (corr_sums, lr_factor, filtered_covariance, lr_potential, field_derivatives),
`src/core/mass_transfer.py` (Heun with AR(1) draws), `src/core/gauge_drift.py` (limit_drift's Schur path), `experiments/exp_gauge_label.py` (A2–G),
`experiments/exp_gauge_classes.py`, the new tests, at 9143f28. The report as delivered (the referee's scratch paths are under the session's scratchpad).


**What I reviewed.**
- The drafts:
  - label/D41_r29.md, read three times: as copied at 02:08, again with (iv)(e) added at 02:14, and again with the revised (ii) "second mechanism" at 02:59.
  - label/pieces_D41_r29.md and label/paper15_r29.md.
  - multi/D42_r29.md and multi/pieces_D42_r29.md.
- The dry runs label/t_report.md and multi/t_classes.md.
- The code. It is now committed as 9143f28; the committed tree equals my copy except for one extra G2 report column (slow/dither_n × 2ττ_c/m).

**Main findings.** The core of (iv), the linear-response filter f_il, survives every test I could make, including hλ_max = 1.56. The revised second-mechanism proof is correct at leading order. The new extensions are overstated:
- **(iv)(e) and Part G.** The slow-noise selection is multistable. At Part C's γ its target stays interior and is not "toward the saturation".
- **(c) and (e).** The trace-to-log-det passage and "a slow measurement noise is a dither" hold only for noise with white innovations (AR(1)/OU, block).
- **D‴.42.** The "internal differences' contribution" is a finite-difference artefact.

## What I checked and how

**Setup.**
- Scratch copy: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review30/repo. I re-copied exp_gauge_label.py after the Part G addition.
- Scripts are in …/review30/work/ and outputs in …/review30/out/.
- I made no changes to the repository. Its HEAD moved only through the coordinator's own commits.
- Every run used OMP_NUM_THREADS=1 and nice 10, with at most 2 processes, each under 15 min, and no EXP_CACHE.

**Re-derived by hand.**
- (iv)'s per-step expansion, term by term at O(σ²): the k ≥ 3 Taylor terms, the O(n²) parts of F and U, the evaluation point of ∂²Φ, the mean shift of n, third moments, and the block phase average.
- f_il and Q_il, and the continuum limits.
- The exact AR(1) form f_ii = 2(1−ρ)/((1−ρa)(1+a)).
- The slow limit at finite h.
- Heun's exact scalar identity Q′ = Q[1 − (h²/4)(c − c̃)² + (h⁴/4)c²c̃²].
- The gauge identity u_k = −g_k∂_{α_k}J.
- D‴.42's Jacobi blocks.

**Numerics.**

| part | runs | script and output |
|---|---|---|
| Scalar model | h·λ = 0.21–1.78; ρ = 0, ±0.9, 0.99, block 4; 4000 paths | scalar_heun |
| Part C, h = 0.25 | ρ = 0.9 and 0.95, 2 seeds × 6·10⁴ slow units; fresh, 2 × 1.6·10⁵ | c_run*.py; c_pred_path.txt |
| Part C, h = 8 (hλ_max = 1.56) | fresh, 8 seeds × 3.2·10⁵; v-scan at 5·10⁻⁵ and 2·10⁻⁴; ρ = −0.5 | c_run*.py; c_pred_path.txt |
| Part C, increments | the increments' modal second moments at h = 8 | c_incr_h8.txt |
| Part E | ρ = 0.99 and 0.999 at full size | e_checks.txt |
| Part F | P_F = 128, T_F = 200; Heun, slow noise, a 300-unit burn-in | f_checks*.txt |
| Part D | at full precision | d_check.txt |
| A2 | γ = 10⁻¹ … 10⁻⁶ | a2_scan.txt |
| Part G checks | the identity Σ = −3AM with M built by hand; slow limits; drifts at γ = 10⁻² … 10⁻⁶ and near the saturation; Gaussian-correlated and block slow noise | g_checks.txt |
| Part G continuation | fine grid (90 steps in τ_c), two guards | g_cont_*.txt |
| Part G zero set | multi-start zero search, and flow integration for the basins | g_multi_*.txt, g_verify.txt |
| D‴.42 | c_m by an explicit unimodular change of variables; Part F's mask rows at three finite-difference steps, eig vs Schur; Schur vs eig at Part B/D/F's points | d42_*.txt |

**Tests.**
- The 17 new tests pass, and so do the 59 gauge tests.
- I ran 45 mutants of the new code (work/mutate.py, out/mutants.txt).

**Literature.** Web search only: arxiv, ncbi and KIAS fail on DNS. Each source is flagged below as snippet or memory.

## Verdicts

### D‴.41

**(ii), the first mechanism.** Correct.
- A = −MH holds exactly, which makes K_hM and K(A)M symmetric. M_flat = blockdiag(G, γ(diag(1/β₁,₂) + 11ᵀ/β₃)) is symmetric.
- Σ + 3AM = 5.5·10⁻¹¹ to 1.2·10⁻¹⁰ of Σ. The cohesion's Jacobian at the m = 3 fit is exactly 0.
- B is antisymmetric on the whole flat tangent space.

**(ii), the second mechanism (revised proof).** Correct at leading order in γ.
- u_k = −g_k∂_{α_k}J holds exactly: u_k = ⟨r a_kφ_k⟩ and ∂_{α_k}J = −β_ka′_k⟨rφ_k⟩.
- Σβ_ju_j is also a total derivative. So the identity integrates along the full flow; no frozen-mass assumption is needed, only convergence to M, which holds near M (the positional Hessian is positive definite, slowest rate 3·10⁻⁴).
- Neither ψ's position part nor the coordinates' term enters. ψ is tangent at Φ(z), so its position part is fixed by its mass part. PΣ_CPᵀ = 0 because ψ = 0 on M.
- Caveats:
  - The constant behaves like g_k ∝ 1/a′, so it blows up at the saturation.
  - The O(γ) law is asymptotic. With γ = 10⁻¹, 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶ in that order (W₀):
    - offset: 6.9·10⁻², 2.1·10⁻², 2.4·10⁻³, 2.4·10⁻⁴, 2.4·10⁻⁵, 2.4·10⁻⁶;
    - cos 2t: 8.2·10⁻², 2.5·10⁻³, 5.4·10⁻⁴, 5.7·10⁻⁵, 5.7·10⁻⁶, 5.7·10⁻⁷;
    - local exponents run 0.35–0.75 between 10⁻¹ and 10⁻² and reach 1.00 only for γ ≤ 10⁻³.
  - sin t is not O(γ²). At W = 1/3 its exponents are 1.01 and 1.03 for γ ≤ 10⁻⁴ (coefficient about 6.5·10⁻³). At W₀ the values fall below the white-noise finite-difference floor (≈ 5·10⁻⁹) from γ = 10⁻⁵ on.
  - D‴.42(iv)'s non-vanishing residual is consistent with the mechanism but does not discriminate. It is computed for white noise, where the antisymmetry also fails.

**(ii), the coloured-noise identity.** It is the right first-order (k = 1) term of tr(K(A)MB) = 0. It holds to 1.3–2.5·10⁻⁸ at every γ, the same finite-difference floor as the white cancellation. So the check is meaningful, but it tests only k = 1.

**(iii).** Correct.
- Direct minimisation of the a-potential gives W_ℓᵃ = (0.29674, 0.45126, 0.25200), the same as w_label(coord='a').
- a_coords is implemented correctly.
- The coord='a' code is untested: mutants A1 and A2 survive.

**(iv), the filter.** Complete at O(v).
- For Gaussian draws, every other term is O(v²); non-Gaussian third moments enter at O(v^{3/2}). The revised proof now says both.
- Block noise: the phase average gives exactly ρ_d = (1 − d/k)₊.
- The cross term's mean vanishes, as above.

**(iv), numbers.** Measured over the path-averaged linear-response prediction, per component (W₁, W₂, W₃):

| run | measured / prediction | the per-step law / prediction |
|---|---|---|
| ρ = 0.9, h = 0.25 (2 seeds × 6·10⁴) | 1.001 ± .034, 0.991 ± .041, 0.974 ± .053 | 1.18, 1.15, 1.09 (rejected at 3.6–5σ) |
| ρ = 0.95, h = 0.25 (2 seeds × 6·10⁴) | 0.995 ± .108, 0.990 ± .081, 0.980 ± .035 | 1.38, 1.32, 1.22 (rejected at 3.5–6σ) |
| fresh, h = 0.25 (2 seeds × 1.6·10⁵) | 1.006, 0.999, 0.986 | — |
| fresh, h = 8 (8 seeds × 3.2·10⁵) | 0.992 ± .006, 1.003 ± .004, 1.036 ± .025 | the law is off by 39% and 30% on W₁, W₂ |

- At h = 8, measured over the Itô drift is 1.64, 1.43, 1.02.
- The v-scan is flat: there is no O(v²) effect.
- The increments' modal second moments match filtered_covariance to 1–2%.
- Per-seed scatter at h = 8 is heavy-tailed (W₃ ranges 0.93–1.16). A 2-seed standard error is meaningless there; my early "7 of 7 seeds above" was a fluctuation.
- ρ = −0.5 at h = 8 (2 seeds) is unresolved (W₁ 0.90 ± 0.03).
- The scalar model: E[c²]/σ² equals lr_factor to 0.1% at hλ up to 1.78.
- Part E: ρ = 0.999 gives 0.977–1.024 and ρ = 0.99 gives 0.988–1.012. The prediction varies by a factor 2.17 across coordinates.
- Part D: ρ = 0.999 gives 0.1916 of the per-step law.
- paper15's "0.19", "2.4%" and "factor 2" are therefore correct.

**(iv), Heun's scheme.** It is a consistent integrator of the ODE driven by the piecewise-constant AR(1) path. That path differs from the OU-driven ODE at O(h/τ_c). Its "zero" holds only as h → 0:
- Heun's selection over Euler's is (hλ)²/4 times the ratio of the E[c²]: 1.1% at hλ = 0.21, 14–24% at 0.85, up to 84% at 1.78 with ρ = 0.99.
- At Part C's steps it is ≤ 6·10⁻⁴.

**(iv)(a)–(d).** The algebra is correct.
- The (1 + 2xx′/(x + x′))/((1 + x)(1 + x′)) limit and Ψ = log(1 + λτ_c)/τ_c check out.
- On the gauge, cos(slow, fresh) = −0.99 at W₀: the colour nearly reverses the selection.
- Damian–Ma–Lee's normalisation (from memory): R = −(1/2η)tr log(1 − η∇²L/2), with tr Ψ = 4R, which is consistent.

**(iv)(e).** Its components:
- **The identity.** Σ = −3AM holds.
- **At finite h.** The limit is m·QM_hQᵀ with a per-mode factor 2/(2 − hλ), not mQMQᵀ. The coordinator's 0.6% residual (my 5.3–6.3·10⁻³) is exactly this factor; with it the residual is 6.9–7.4·10⁻⁴.
- **Near the saturation.** λ_min falls to 3.6·10⁻⁵ at W₁ = 1.02W_min, and the residual there is 2.4% even at τ_c = 10⁶.
- **Direction.** cos(slow, dither) = 1.00000 at γ ≤ 10⁻⁴, including near the saturation. Near the saturation it is 0.76 at γ = 10⁻² and −0.21 at 10⁻³.
- **Magnitude at γ = 10⁻².** The slow noise is 1.9–2.7 times its normal-part equivalent, because the mass block of QMQᵀ adds to it. The dither's own drift is 3.0–3.6 times its normal part.
- **The target at γ = 10⁻² (Part C's γ).** It converges to the interior point (0.2647, 0.5209, 0.2144), with W/W_min − 1 = (0.06, 1.98, 0.07), stable up to τ_c = 10⁵.
- **Multistability.** The multi-start search and flow integration find three stable zeros:
  - at τ_c = 100: (0.2779, 0.4949, 0.2272) [the continued branch], (0.2821, 0.2260, 0.4919) and (0.5419, 0.2266, 0.2315);
  - at τ_c = 10⁴: (0.2649, 0.5204, 0.2146) [the continued branch], (0.2652, 0.1898, 0.545) and (0.5948, 0.190, 0.2152);
  - for τ_c ≥ 10³: W♯ ≈ (0.375, 0.300, 0.325) is an unstable node (exponents +27 to +34), with saddles between.
  - At τ_c = 10⁴, the start (0.40, 0.27, 0.33) next to W₀ ends at (0.5948, 0.19, 0.2152), not on the continued branch.
  - "The target moves continuously" therefore describes one branch of several.
- **At γ = 10⁻⁵.** The branch approaches the corner, with W/W_min − 1 = 1.3·10⁻⁴ and 1.6·10⁻⁴ at τ_c = 10⁵, and stays stable. G1's "no interior zero at 2·10⁴" is an artefact of the guard 1.0005·W_min: with a guard of 1 + 10⁻⁶ the branch continues.
- **The targets sit where the slow noise is not dither-like.** At the targets the softest normal rate × τ_c is 0.005–0.2. So the masses never reach the saturation.
- **AR(1)-specific.** A smooth slow noise (Gaussian step correlations) has f_ii ≈ (λτ_c)⁻² instead of the AR(1)'s ≈ (λτ_c)⁻¹: 2.6·10⁻⁵ against 5.1·10⁻³ at λτ_c = 195. Its Ψ is bounded, and its direction differs:
  - from the dither's: cos 0.977–0.983 at W = 1/3 (γ = 10⁻⁵) and 0.926–0.950 (γ = 10⁻²);
  - from the AR(1)'s: 0.876–0.982.
  - Block noise behaves like the AR(1).
- **Size.** At fixed intensity the effect scales like v·h/τ_c², so it is unobservable on any run's time scale.

**(v).** Fine.

**(vi).** The identities are exact.
- Part F tests fresh noise (1.0000 of the prediction) and Heun (0.3%), but not slow noise. Its noise floor is ±2·10⁻⁶, which is 0.3% of the fresh slope. The predictions are −1.4·10⁻⁶ at ρ = 0.99 and −2.8·10⁻⁸ at ρ = 0.999.
- So "the law holds as stated" on the generic network is supported only for ρ ≤ 0.9.

### D‴.42

**(i).** Verified numerically.
- The barycentre block is ℓ₁ + 9 log 2 (−23.643383); the η block is −52.369544. So c₆ = −46.131219, and c₉ = 9 log 3 − 102.149949 = −92.262438.
- The off-diagonal blocks are 6·10⁻¹⁰ (relative).

**(ii).** Restricted to the isotropic dither, and correct there. The unequal-split O(γτ) effects reproduce. But the equal split's "(the latter to the internal differences' contribution)" is a truncation error:
- For mask (0,1,1) at γ = 10⁻⁵, the W₁ component per γτ is:

  | h2 | m = 6 | per-class law |
  |---|---|---|
  | 10⁻³ | 2.77·10⁻⁷ | −1.73·10⁻⁷ |
  | 2·10⁻³ | 1.31·10⁻⁶ | −4.86·10⁻⁷ |
  | 4·10⁻³ | 5.46·10⁻⁶ | −1.74·10⁻⁶ |
  | Richardson extrapolation | −6.83·10⁻⁸ | −6.84·10⁻⁸ |

- The extrapolated ratio is 0.998 on W₁ and 1.002 on W₂.
- At the experiment's h2 = 10⁻⁴ the predicted truncation errors give 0.934 and 0.949, which is the dry run's 0.9353 and 0.9484. At γ = 10⁻⁶ in GAMMAS_F the error will be about 50%.
- Plausibly the internal differences contribute nothing exactly for a coordinate mask: η's law is Gibbs coordinate by coordinate, and its block does not depend on W.

**(iii).** Correct.

**(iv).**
- The invariant subspace checks out factor by factor: β_i = mW_k/n = KW_k, the cohesion force is 0 at coincidence, the fields and the log-mass field are identical. The class drift over m = 3's is 1.0000.
- The fractions' Itô drift is the full trace (5.515·10⁻⁷ is exactly minus the Wong–Zakai term) plus the residual.
- "A few percent" should be 9.6% on W₁ at (0.2, 0.7, 0.5).

### Code and tests

**limit_drift's Schur path is correct.**
- X_N is supported on the normal space.
- A_D = (A − σP)⁻¹Q is the group inverse, and AA_D = Q.
- The tangent-tangent block is irrelevant because P∂²F[t, t′] = 0.
- At well-conditioned points it agrees with eig to 1.7·10⁻⁸–9.7·10⁻⁶. The auto-switch fires only at n = 3's equal splits (cond 1.5–3.7·10¹⁸), where eig returns up to 4·10³².
- It leaves the eig path unchanged. One behaviour change is untested: cases with |Im λ| between 10⁻¹² and 10⁻⁶ (relative) used to raise and are now silently routed to Schur.

**Mutants.** 38 of 45 are killed. The 7 survivors:
- Untested code: A1 and A2 (the coord='a' closed form and W_ℓᵃ), and H9 (Heun with dither).
- Equivalent mutants: C4, S5 and G6.
- G9 ("Schur everywhere") passes because the two methods agree.
- The experiment-side prediction code (c_prediction, d_prediction, the prediction parts of e_job and f_prediction, g_drift, g1_job, g2_job) has no tests.

### Literature
- **Missing:**
  - Kim, Oh, Kwon, Choe, Baek, Bae, NSPCS2026 abstract (KIAS; snippet only): anti-correlated noise gives a Hessian-dependent implicit regulariser at finite learning rate, with a trade-off set by the correlation time. This is close to (iv).
  - Beneventano, arXiv:2312.16143 (snippet): SGD without replacement penalises a weighted trace with a Hessian-function weight.
  - Lucchi–Proske–Orvieto–Bach–Kersting, NeurIPS 2022, arXiv:2209.09162 (snippet): correlated parameter noise (fPGD).
  - The AOUP and UCNA literature for the filter's continuum limit: the trap's effective temperature T/(1 + kτ) (snippet, arXiv:1910.02667); the det(I + τ∇²U) factor of Maggi et al., Sci. Rep. 2015 (memory); Jung–Hänggi 1987 and Hänggi–Jung 1995 (memory).
- **(e)'s target.** It is Li–Wang–Arora's isotropic-noise log-det drift (slides, snippet).
- **DML's white case** is itself derived from the discrete OU's stationary covariance.
- **The "New" claims** should narrow to the coloured filter for label-noise drift along M.

## Table

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ac-1 | major | (iv)(e) and Part G: the slow-noise selection is multistable (three stable zeros at τ_c = 100 and 10⁴ at γ = 10⁻²; W♯ an unstable node; saddles). The continued branch is not reached from W₀'s neighbour. At Part C's γ the target is interior (0.2647, 0.5209, 0.2144), not toward the saturation. Near the saturation cos(slow, dither) is 0.76 (γ = 10⁻²) and −0.21 (10⁻³). Magnitudes are 1.9–2.7× the equivalent. | D41_r29 (iv)(e), Status, "New": restrict to the slow-replicator limit and drop "the target". exp_gauge_label G: multi-start zero set with stability, the basin of W₀, every branch continued. |
| D3ac-2 | major | (c)/(e) "as the noise slows … log det", "a slow measurement noise is a dither": AR(1)/OU-specific. A smooth slow noise has f ≈ (λτ_c)⁻², a bounded Ψ and a different direction (cos 0.93–0.98). | Name AR(1)/Markov noise in (e)'s title, (c), paper15 "Claim", BULLET, README, ROW0, AUDIT. Add a smooth-correlation row to G2. |
| D3ac-3 | major | D‴.42(ii) "the internal differences' contribution" at the equal split is O(h2²) truncation error. Richardson gives a ratio of 1.00; at γ = 10⁻⁶ the F rows will be about 50% error. | D42_r29 (ii), proof, Status. exp_gauge_classes F: h2-scan or Richardson; flag rows below resolution. |
| D3ac-4 | minor | (e)'s limit only as h → 0 (factor 2/(2 − hλ), which explains the 0.6%) and away from the saturation (λ_min → 0, 2.4% at τ_c = 10⁶). The targets sit where λ_softτ_c ≈ 0.005–0.2. | (e) statement/proof; g2_job's lim_res with the factor. |
| D3ac-5 | minor | G1 stops at γ = 10⁻⁵ because of the 1.0005·W_min guard; with 1 + 10⁻⁶ the branch continues to τ_c = 10⁵. | g1_job guard and its report row. |
| D3ac-6 | minor | (e)'s effect scales as v·h/τ_c² and is unreachable in practice; the drafts do not say so. | (e), "New", paper15. |
| D3ac-7 | minor | Heun's "zero" holds only as h → 0: its O(v) selection is (hλ)²/4 × (E[c²] ratio) of Euler's (up to 84% at hλ = 1.78). | exp_gauge_label C/D–F text, pieces, paper15; the (iv) Status wording. |
| D3ac-8 | minor | Second mechanism: correct, but the constant ∝ 1/a′; O(γ) only for γ ≲ 10⁻³ (6–12% at 10⁻¹); sin t is O(γ), not O(γ²); D‴.42(iv) does not discriminate. | D41 (ii) proof/Status; the A2 report. |
| D3ac-9 | minor | Part F cannot test ρ ≥ 0.99 (floor ±2·10⁻⁶ against −1.4·10⁻⁶ and −2.8·10⁻⁸). | F report: print the floor; (vi) wording. |
| D3ac-10 | minor | c_prediction is evaluated at W₀ only. Over S_C = 2·10⁵ the prediction's W₃ shifts by 0.7–2% (4.4% at h = 8, S = 3.2·10⁵). Seed scatter at large h is heavy-tailed. | summarise_C: average the prediction along the path. Add a large-step scheme with ≥ 8 seeds as the "any step" test. |
| D3ac-11 | minor | Mutants A1, A2, H9 survive; the experiment-side prediction code is untested. | tests/test_gauge_label.py. |
| D3ac-12 | minor | Missing literature (Kim et al. 2026, Beneventano 2023, Lucchi et al. 2022, AOUP/UCNA, Li–Wang–Arora's log-det); the novelty is overstated. | [S69], AUDIT, "New", paper15. |
| D3ac-13 | minor | D‴.42(iv) "a few percent" is 9.6% on W₁. | D42 (iv), pieces. |
| D3ac-14 | cosmetic | "For any step" should read h < 2/λ_max; in (b), ∂²F[X] = −½∇trΨ only after P_T; G2's componentwise ratios; Part D prints 4 decimals; the dry-run C used L_MAP = 5·10³ (a factor e^{−1.5}). | D41 (iv); exp_gauge_label report. |

## Lesson

The derivations held up again: the filter's algebra and the gauge identity are sound, and the filter matches every quantitative test, including hλ_max = 1.56 once enough seeds and a path-averaged prediction were used. The failures come from reading a global statement off a few points. (e) was inferred from one continued branch and three cosines. The drift field itself has several attractors, a repelling W♯, saddles and basins, and its slow limit depends on the noise's short-time roughness, on h and on the distance to the saturation. A sub-percent finite-difference residual was given a mechanism ("the internal differences' contribution") without varying the step. Before naming a target or a contribution, map the whole zero set and its basins. Vary the noise model's spectrum, not only its correlation time. Extrapolate any small numerical effect in the finite-difference step.
