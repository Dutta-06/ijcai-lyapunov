# Referee report D3ag (round 34): Proposition D‴.44, the replicator's own masses on a kernel dictionary under measurement noise

## What I checked and how

**Rules followed**
- I worked only in a scratch copy of the repository. I wrote nothing to the repository and ran no git command that changes anything.
- The D‴.44 files I reviewed are byte-identical to the repository's current versions.
- Every Python run used `OMP_NUM_THREADS=1` and `nice -n 10`. I ran at most one heavy process at a time, each under 15 min, and every
  flow had `max_secs`.

**Re-derived by hand**
- (i):
  - ∂J/∂α_i = −β_ia′_i⟨rφ_i⟩ and ∂J/∂ϱ_i = −2ϱ_ia_i⟨rφ_i⟩, so the mobility is diag(G, γ/4).
  - The replicator's λ is the tangential projection on the sphere.
  - The noise covariance is mv_h·P∇²JP, so K = (m/N)J_θPJ_θᵀ and D_i is as stated.
  - For the conserved mass, K − (γ/N)yyᵀ.
  - The label noise's filtered moment is Σ_f = c·H̃f(H̃), which is diagonal in the normal modes.
  - The Itô term is ϱδℓ²/8, and its cancellation follows.
  - The conserved law's leftover is the cross term −2a_iφ_iᵀf(K)y.
- (ii): D*, its two branches, and the minimiser's equation 2Gm|w|(1−t⁴) = γa_max³t³. The paradigm's own gauge equilibrium is
  Gm|w|(1−t⁴) = γa_max³t³, and it holds for every colour because ω_i cancels.
- (iii): the inequality chain, both Cauchy–Schwarz steps, and the limit at the boundary.
- (iv): K for the null kernel.

**Regeneration and transcription**
- `exp_mass_dictionary.py --from-raw` with `MD_OUT` set to scratch reproduces the md byte for byte.
- Every number in D44.md's numerics paragraph matches the md.
- I also went through the per-instance values in `_raw.npz`.

**Independent checks** (script → what it does)
- `b_katz.py`: the closed-form drift against Katzenberger's drift (`scheme_drift`) at Part B's own parameters (h = 1, fresh and AR(1)
  τ_c = 10³), for all three laws.
- `b_scan.py`, `b_late.py`: the paradigm's AR(1) scheme with every path's state projected onto M by the noiseless flow (Katzenberger's
  limit map Φ), compared with the reduced flow. Run on early windows and on a late window starting from the flow's state at s = 1500,
  for σ = 0.03 and 0.015, h = 1 and 0.5, conserved and free.
- `sat_ratio.py`, `cross_term.py`, `cross_sat.py`: the size of the Itô term and of the centring's cross term along the paradigm's own
  white flow.
- `nongradient.py`: loop integrals of P⁻¹·drift around small loops on M.
- `quick_checks.py`: the bound with the draft's conserved K at 246 random points; R_null − R_free for slow and white Ψ.
- `dstar_check.py`, `dstar_min.py`: D*'s convexity, and the support of the D* minimiser.
- `flow_status.py`: three off-basis-pursuit free ends re-run with solve_ivp's status recorded per piece.
- `timeout_est.py`, `timeout_rerun.py`, `stall_state.py`: Part D's left-out flow.
- `bound_end.py`: the white potential and the bound at flow ends.
- `mutate.py`, `test_vacuity.py`: 34 single-replacement mutants run against `tests/test_mass_dictionary.py`.

**Tests**
- 10 passed in 3.5 s.

**Literature**
- Web search only. ar5iv.arxiv.org and projecteuclid.org could not be resolved in this sandbox. Anything marked "(memory)" is unverified.

## Verdicts

### (i) The selection, the Itô term, the cancellation

**Correct**
- **Flat coordinates, metric and 0-homogeneous extension.** The conserved flow is the gradient flow of J̃ on the sphere.
- **K's nonzero spectrum equals the normal rates.** This is consistent with Part A's 10⁻¹¹ check.
- **The closed form equals Katzenberger's drift**, also at Part B's own h = 1 and τ_c = 10³ (0.96–1.96·10⁻⁶, all three laws).

**Radial renormalisation: right conclusion, wrong reason**
- The proof says the renormalisation is harmless because J̃ is constant along rays. That is not why.
- The unrenormalised scheme does have a radial selection drift (tr Ψ(H̃) is not 0-homogeneous), and the renormalisation adds
  −ϱ tr Σ_ϱϱ/(2m).
- Both terms are P-normal to the sphere. The projection removes them, which is what the code's `F1` does. The proof should say this
  instead.

**"O(γ/G) relative" is false where the selection decides the end**
- Per atom, the Itô term over the amplitude part of the ϱ-drift is γa²/(2Gβa′²), which blows up as a′ → 0.
- Along the paradigm's own white flow (Part B's instance), |Itô|/|flat| is:
  - conserved: 0.010 (s = 0), 0.016 (240), 1.28 (10⁴), 0.95 (10⁵);
  - free: 0.019, 0.13, 0.96, 0.98.
- At a min-norm fit placed at the paradigm's gauge it is 0.98.

**The cancellation for free masses**
- It is an exact identity at order v for every colour. It even holds for any node covariance: in ℓ = log β, the mass's own noise field
  γa_iφ_i·ξ does not depend on ℓ_i.
- But it holds only at order v. The free scheme's projected drift is 0.892 ± 0.001 of the reduced flow's at σ²τ_c = 0.9.

**Conserved mass: "up to the centring's cross term" hides the leading term**
- The leftover is exactly the projected cross term (to 1–5·10⁻¹³).
- Its size relative to the net paradigm drift, along the conserved white flow on two instances:
  - 0.3–1.9% at s = 0;
  - 18–21% at 10³;
  - 2.8–3.2× at 10⁵;
  - 5.8× and 188× at 10⁶.
- At the ends it balances the rest of the drift. Part A's 3.8·10⁻³ and 8.1·10⁻³ come from moderate random points and are
  unrepresentative.

**A consequence the draft does not draw**
- After the cancellation, the paradigm's drift is not a gradient on M. The removed covector Σk_iϱ_ia_i²dϱ_i has
  dν = Σ2k_iϱ_ia_ia′_i dα_i∧dϱ_i ≠ 0 on T_M.
- Loop integrals of P⁻¹·drift around loops of radius ε = 10⁻², 5·10⁻³, 2.5·10⁻³ confirm it:
  - flat drift: −4.5·10⁻¹⁶, −2.8·10⁻¹⁷, −1.8·10⁻¹⁸ (ratio 16, i.e. zero);
  - paradigm drift: 2.1·10⁻⁸, 5.3·10⁻⁹, 1.3·10⁻⁹ (ratio 4.0, i.e. ε² times a nonzero 2-form).
- So no statement in (ii)–(iv) about minimisers of R applies to the paradigm's update — which is the one Parts B–E run.

### (ii) Free masses: the penalty and the support count

**D* and its branches: correct**
- I reproduce 0.9983 and 0.9882.
- D* is within 1% of the ℓ¹ branch only for |w| > 0.051 (6× the threshold), and within 5% of the 4/3 law only for |w| ≲ 10⁻⁴.
- D* is convex, with D*′(0) = 0.

**"Global minimisers are those of Σ|φ_i|²D*" is right for the flat update; "weighted basis pursuit above the threshold" is misleading**
- Because D*′(0) = 0, the minimiser is dense: 9–12 of 12 atoms are above 10⁻⁶, with the off-support ones at 10⁻³–10⁻⁵.
- It sits 1.4·10⁻⁴ to 5·10⁻³ from basis pursuit on six instances, and up to 1.5·10⁻² across the raw data.
- It is never the truth, even where basis pursuit is: C k = 2 i = 0 is 4.4·10⁻³ off. The flat update's truth fraction is 0.12/0,
  against basis pursuit's 0.81/0.38.
- The off-support atoms are always below the threshold, so "above the threshold" never covers them.

**The flat flow does not always end at that minimiser**
- On C k = 2 i = 1 and 4, D k = 2 i = 7, and D k = 3 i = 1, it ends exactly on the sparse truth (5·10⁻¹⁴), 3.3·10⁻⁴ from the dense D*
  minimiser, because dead atoms are absorbing.
- The flow-finishes test asserts exactly this sparse end.

**The paradigm's update does not reach the D* split**
- Its gauge settles at Gm|w|(1−t⁴) = γa³t³: t = 0.484 against D*'s 0.594 at |w| = 10⁻³, and 0.980 against 0.990 at 0.1.
- Round 33 (D3af) already flagged this. D44.md still treats D* as the paradigm's relevance.

**Slow-noise support count: the stated leading term is wrong beyond N**
- Since rank K ≤ min(k, N), the leading term is min(k, N)·log τ_c/τ_c. Every C conserved end has 7–12 alive atoms with N = 5, so they
  all get the same count.
- It counts eigenvalues above 1/τ_c, not atoms with β_i > 0.
- The O(1/τ_c) term contains Σ log κ_j and is not uniform.
- "Whatever the parametrisation" holds only for a gradient selection.

### (iii) The conserved mass inflates the fit

**The chain and Σβa′ ≥ E: correct.**

**tr K ≥ (G/(Nm))E²/Σ|φ_i|⁻²**
- True. With the draft's own conserved K (minus (γ/N)yyᵀ) I found 0 violations at 246 random points, minimum margin 2.4·10⁻³.
- But the proof omits the second Cauchy–Schwarz it needs: (γ/N)|Φw|² ≤ (γ/N)Σ(β_i/m)a_i²|φ_i|² = (γ/N)Σ|φ_i|²|w_i||a_i|.
- Part A and the test check tr K_free, not the stated K.
- The boundary limit is tr K_free's.
- The "fits with ‖w‖₁ = a_max" are limits, not points of M: dead atoms, α = ±∞.

**The bound on E: correct but uninformative**
- It is vacuous for a_max ≥ (Gm/(γΣ|φ|⁻²max|φ|²))^{1/2} = 3.16 here. At Part E's a_max = 4 it gives E/m ≤ 1.27.
- It is 13–26× loose at a_max = 1 (bound 0.316, against 0.025 flat and 0.012 paradigm).
- The bound scales ∝ a_max, while Part E's E/m (0.0021, 0.0119, 0.0585, 0.347) grows ∝ a_max^{2.5}.

**"Every kernel used near saturation — the opposite of pruning" and "the white selection keeps every atom alive"**
- True of the flat update: 12/12 and 16/16 used and alive.
- False for the paradigm's update:
  - in C, 5–6 (= N + 1) of 12 atoms are used at every conserved end, and 6–12 are alive;
  - in D, 5–9 (≤ N + 1 = 9) of 16 are used.
- These are vertex-sized supports of {Φw = y, ‖w‖₁ ≤ a_max}: the paradigm prunes half the dictionary while inflating ‖w‖₁.
- At a_max = 4, max|a|/a_max at the end is only 0.82–0.89. The saturation mechanism needs |w_i| ≫ γa_max³/(Gm) = 0.53, a condition the
  proposition never states.
- "Keeps every atom alive" is a heuristic: it compares two costs, it does not prove anything.

**"Inflate under every colour (numerical)"**
- Only two colours were run, and only as reduced flows (see Part B).

### (iv) The null kernel

**White noise, flat update: correct.** tr K_null − tr K_free = −1.2423·10⁻³ = −(γ/N)|y|² at six points, a constant on M.

**Coloured noise: "proved" is false**
- The null kernel's mass lies on the sphere, so K_null = K_free − (γ/N)yyᵀ; `gram()` computes exactly this.
- By the determinant lemma, R_null − R_free = (1/τ)log(1 − τ(γ/N)yᵀ(I + τK_free)⁻¹y). This is not constant on M: −3.5 to −5.4·10⁻⁵ at
  six points, τ = 10⁴.
- In a one-atom example, the mass part of D is scaled by β_0/m.
- Under the paradigm's update, the null-kernel law also keeps the cross term.

**The numerics contradict the title's "restores the colour's selection"**
- C k = 2, slow noise: the null kernel misses the truth on 6 (τ_c = 10³) and 7 (10⁴) instances where both basis pursuit and free masses
  find it.
- It ends within 10⁻³ of basis pursuit on 0 of 3 (k = 2) and 0 of 10 (k = 3) of the instances where basis pursuit fails.
- White noise: its ends differ from the free ends on 4 of 16 instances (i = 0, 6, 10, 13).

### Numerics

**Transcription is exact. Interpretation is not.**

**"Free masses end at basis pursuit (median < 10⁻³)" and "recovering the truth where basis pursuit does"**
- The medians are taken over instances where basis pursuit equals the truth, which hides the informative ones.
- C k = 3 free white: within 10⁻³ of basis pursuit on 9/16 instances, but on only 3/10 of the instances where basis pursuit fails. The
  largest distance is 0.21.
- C k = 2 at τ_c = 10³: 13/16 within 10⁻³; the largest distance is 0.27.
- The free flows recover the truth where basis pursuit fails: C k = 2 i = 6; C k = 3 i = 5, 6; D k = 2 i = 0; D k = 3 i = 0, 4.
- They miss it where basis pursuit gets it: C k = 2 i = 0 at τ_c = 10³.
- I re-ran three of the off-basis-pursuit ends:
  - all 60 pieces finished with status 0, and the raw values are reproduced;
  - the ends are settled: |Δw| ≤ 7.5·10⁻⁸ from s = 10⁷ to 10⁸;
  - they land on other 5-atom fits, or on the truth.
- Round 33's D3af-1 made the same finding and it has not been absorbed here.

**"The flat update ends at the D* minimiser"**
- True for the median (1.1·10⁻⁴ and 1.3·10⁻⁴). The largest distances are 9.3·10⁻³ (C k = 3) and 7.6·10⁻³ (D k = 3).

**Part D's left-out flow (k = 2, i = 5, free, flat)**
- `report()` puts a 7-instance flow mean next to an 8-instance truth mean (0.453; 0.450 over the same 7) and an 8-instance
  basis-pursuit fraction.
- The flow had already reached the truth (6·10⁻¹⁰) by s = 1.67·10⁶. It then stalls while dead atoms revive: their masses are
  10⁻⁸–10⁻⁵⁶ with log-mass rates up to +5·10⁻³.
- My 720-s rerun advanced only from s = 1.667·10⁶ to 1.677·10⁶.
- Its outcome matters for the row: the truth fraction would be 0.25 at its last state and 0.125 at the D* minimiser, against the 0.14
  printed. The ‖w‖₁ mean would move from 0.422 to 0.428.

**Part E**
- Covered under (iii). In addition, at a_max = 0.5 the start is already at 0.987, so that point shows nothing.

**Part B, fresh noise**
- "Within 0.002" hides deviations of 3–5.6 standard errors (gap/se 4.0, 5.6, −3.3, −3.1). This is an O(v) offset of E[w(x)].
- Projected means agree with the flow to 0.98–0.99 ± 0.005–0.008.

**Code**
- `flow()` never checks `sol.status`. A failed LSODA piece would be silently treated as having reached `hi`. This is latent: none of my
  reruns failed.
- The null kernel's E/m in `summary()` includes the null atom's |w_0|. For C k = 2 white it prints 0.485, where 1 − ‖w‖₁/a_max = 0.514.

### Part B: the AR(1) lag (0.873 ± 0.005 against 0.917)

**What it is not**
- Not an error in the flow. The closed form equals Katzenberger's linear-response drift at B's own h and τ_c to 1.2·10⁻⁶.
- Not an error in the scheme. I re-implemented the update independently and see the same behaviour.
- Not a step-size bias: it gets worse when h shrinks.

**What I measured** — every path projected onto M, on the window s ∈ [1650, 2100] starting from the flow's state at 1500. Projected
advance divided by the flow's advance:

| Law | σ | h | Paths | Ratio |
|---|---|---|---|---|
| Conserved | 0.03 | 1 | 128 | 0.20 ± 0.08 |
| Conserved | 0.015 | 1 | 64 | 0.67 ± 0.07 |
| Conserved | 0.03 | 0.5 | 64 | −0.12 ± 0.18 |
| Free | 0.03 | 1 | 128 | 0.892 ± 0.001 |

- The free law's projections are almost deterministic: standard error 10⁻⁵.
- Early on, the advance up to s = 600 still agrees (0.94 ± 0.06).
- But the paths' mean carries an O(σ²) offset E[w(x)] − w(Φ(x)) of up to 0.009 per component (z ≤ 9.5). That offset is B's early
  "dev" of 0.007–0.011.

**Interpretation: the order-v expansion's limit for the conserved laws**
- The shortfall shrinks with σ and grows as h shrinks.
- That is consistent with an effect of the slow noise acting through the centring that does not depend on h per unit time. It competes
  with the order-v selection, which is a sampled-data effect of size O(h).
- The free law, which has no centring and whose continuous-time loop conserves Q (D‴.45), is barely affected.
- So the slow-noise columns of C–E for the conserved and null-kernel laws are statements about reduced flows, valid only when
  σ²τ_c/h ≪ 1. That condition is untested.

### Tests and mutants

34 single-replacement mutants; 11 killed, 23 survived.

| Group | Mutants | Killed | Why |
|---|---|---|---|
| Cancellation / selection | M01–M08, M25, M26 | 10 of 10 | Katzenberger or manifold test catches them |
| D* | M09–M13, including 2G for G | 0 of 5 | The 4/3-law assertion is vacuous: `np.allclose`'s default atol = 10⁻⁸ against D* values of 1.3·10⁻¹² and 2.8·10⁻¹¹, so it passes even against 0. The ℓ¹ branch is insensitive to G at leading order. |
| Bound / excess | M14–M18 | 0 of 5 | The test is one-sided (≥); the columns have unit norm, so \|φ\|² and \|φ\|⁻² cannot be told apart; a_max = 1 hides a_max factors |
| Noisy scheme | M19–M22 (AR(1) innovation, replicator sign, a′, variance) | 0 of 4 | The only scheme test uses σ = 0 at a point of M, where g = 0 |
| Flow | M23, M24 | 1 of 2 | Mass-rate mutant killed; floor mutant survives |
| Experiment's classification | X01–X08 | 0 of 8 | `summary`, `report`, `c_job` and `filt_of(τ>0)` are never tested |

### Novelty and literature

**The cancellation**
- It is the coordinate dependence of the Itô drift: in log coordinates, noise that enters additively does not select the mass by its own
  noise.
- This is the exponential-learning and stochastic-mirror-descent mechanism: Mertikopoulos–Moustakas, AAP 2010 (the paper's existence is
  search-verified; that it carries this Itô correction is from memory); Raginsky–Bouvrie 2012; Mertikopoulos–Staudigl 2018.
- The draft lists the log-normal correction itself as classical. Its "new" cancellation is a two-line corollary.

**D***
- Its two branches are the label-noise penalty of a quadratic parametrisation (HaoChen–Wei–Lee–Ma 2021; Pillaud-Vivien–Reygner–Flammarion
  2022) and the depth-3 penalty |w|^{2−2/L} with L = 3 (elementary; this is D‴.48's "effective depth-3 unit").

**Support counting**
- Sparse Bayesian learning's log-det (Tipping; Wipf–Rao), and Chen–Kunin–Yamamura–Ganguli's stochastic collapse (NeurIPS 2023), in which
  noise attracts toward vanishing-unit invariant sets.

**Other claims**
- "Every colour" means two colours.
- The null-kernel claim is contradicted (D3ag-4).
- I found no precedent for the conserved-mass inflation.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ag-1 | major | (ii)–(iv)'s minimiser statements describe the flat (√β) update. The paradigm's update — the proposition's subject, and the one run in B–E — is not a gradient on M (loop integrals scale as ε², 2.1·10⁻⁸ at ε = 10⁻², against ≤ 4·10⁻¹⁶ for the flat drift). Its gauge settles at Gm\|w\|(1−t⁴) = γa³t³, not at D*. Its ends differ: ≤ N+1 atoms used, and off basis pursuit on most instances where basis pursuit fails. The numerics are still offered as "the same picture". | D44 (i)–(iv) and Status: say which update each claim is about; state that the paradigm's selection is non-gradient and its ends are numerical; derive its invariant (D‴.45's Q-coordinates) or drop "minimisers" for it. |
| D3ag-2 | major | "O(γ/G) relative" is false at the ends. Along the paradigm's white flow, \|Itô\|/\|flat\| goes from 0.010 to 1.28 (conserved) and from 0.019 to 0.98 (free); it is 0.98 at a saturated gauge. | D44 (i): "O(γa_i²/(Gβ_ia′_i²)), O(1) at the ends". Part A: add saturated points. |
| D3ag-3 | major | The cross term is exactly the conserved leftover, but along the conserved white flow it is 0.3–1.9% of the net drift at s = 0 and 2.8–188× at s ≥ 10⁵, where it balances the rest. Part A's 4–8·10⁻³ is unrepresentative. | D44 (i), numerics (A): report it at the flows' ends; drop the reading "flat drift minus the masses' share" for the conserved law. |
| D3ag-4 | major | (iv) is false for coloured noise and for the paradigm's update. K_null = K_free − (γ/N)yyᵀ, so R_null − R_free is non-constant for slow Ψ (−3.5 to −5.4·10⁻⁵). In C k = 2 under slow noise, the null kernel misses the truth on 6–7 instances where basis pursuit and free masses get it; under white noise its ends differ from the free ends on 4/16. | D44 title, (iv), Status: restrict to white noise and the flat update. Numerics (C): report the disagreements. Add a test of R_null − R_free across M for a slow Ψ. |
| D3ag-5 | major | "Free masses end at basis pursuit" and "recovering the truth where basis pursuit does" rest on medians. Per instance: C k = 3 within 10⁻³ of basis pursuit on 9/16, on only 3/10 of the basis-pursuit-failure instances, max distance 0.21–0.27. The flows recover the truth where basis pursuit fails on 6 instances and miss it where it succeeds on 1. The ends are genuine and settled (rerun: 60/60 pieces status 0, \|Δw\| ≤ 7.5·10⁻⁸). | Numerics (C), (D) and `report()`: per-instance agreement, rows split by `bp_ok`, the fraction within 10⁻³; restate the finding. |
| D3ag-6 | major | "Every kernel used near saturation — the opposite of pruning" and "keeps every atom alive" hold only for the flat update. The paradigm's conserved ends use 5–6 of 12 atoms (5–9 of 16): vertex-sized supports, i.e. pruning. | D44 (iii): restrict to the flat update; describe the paradigm's ends as (N+1)-sparse near-vertices; label "alive" as numerical. |
| D3ag-7 | major | Part B's AR(1) lag is the order-v expansion's limit for the conserved laws, not a bug and not a step-size bias. Projected late-window advance relative to the flow: 0.20 ± 0.08 (σ = 0.03), 0.67 ± 0.07 (σ = 0.015), −0.12 ± 0.18 (h = 0.5); free masses 0.892 ± 0.001. B's "dev" also mixes in an O(σ²) offset of E[w(x)]. The slow-noise claims in C–E rest on reduced flows that are unvalidated at σ²τ_c/h ≈ 1. | Numerics (B) and `b_job`: compare projected means; report the advance ratio with its standard error at σ²τ_c ≤ 0.05 and two values of h; state the validity condition; qualify "inflate under every colour". |
| D3ag-8 | major | Tests kill only 11 of 34 mutants. All D* mutants survive, including the paradigm's 2G (the 4/3-law assertion is vacuous through atol = 10⁻⁸). All bound/excess, all noisy-scheme (σ = 0 test) and all classification mutants survive. | `tests/test_mass_dictionary.py`: use atol=0 and test D*'s minimiser equation; test the bound at saturation with non-unit columns and a_max ≠ 1; test one noisy step's mean and covariance against v_h and AR(1) stationarity; test `summary()` on hand-built states. |
| D3ag-9 | minor | "Weighted basis pursuit above the threshold": the D* minimiser is dense (9–12/12 atoms), 10⁻³–1.5·10⁻² from basis pursuit, and never the truth. The flat flow itself sometimes ends on the sparse truth (C k = 2 i = 1, 4; D k = 2 i = 7; D k = 3 i = 1), 3.3·10⁻⁴ from the minimiser, and the flow test asserts that sparse end. | D44 (ii), numerics (C): call it a dense ℓ^{4/3} perturbation of basis pursuit; report the largest distances. |
| D3ag-10 | minor | The slow count's leading term is min(k, N)·log τ_c/τ_c. It counts eigenvalues above 1/τ_c, not atoms with β_i > 0; the O(1/τ_c) is non-uniform; "whatever the parametrisation" holds only for gradient selections. | D44 (ii), (iii). |
| D3ag-11 | minor | The conserved-K inequality's proof misses the (γ/N)\|y\|² step (the inequality is true). Part A and the test check tr K_free; the boundary limit is tr K_free's; the boundary fits are not points of M. | D44 (iii) proof; label in Part A; test. |
| D3ag-12 | minor | The E bound is vacuous for a_max ≥ 3.16 and 13–26× loose at a_max = 1. Part E's E/m grows ∝ a_max^{2.5}. At a_max = 4 there is no saturation (max\|a\|/a_max 0.82–0.89). The regime condition \|w\| ≫ γa_max³/(Gm) is missing. "Every atom used (11.5 of 12)". | D44 (iii), numerics (E). |
| D3ag-13 | minor | Part D's timed-out flow: the row mixes denominators, and the flow stalls while reviving atoms after reaching the truth. Its fate moves the truth fraction anywhere from 0.125 to 0.25 (0.14 printed). | `report()`: one instance set per row, print n; numerics (D). |
| D3ag-14 | minor | Part B fresh: "within 0.002" hides 3–5.6 standard-error deviations (an O(v) offset). | Numerics (B): projected comparison. |
| D3ag-15 | minor | `flow()` ignores `sol.status`. The null kernel's E/m includes \|w_0\|. | `mass_dictionary.flow`: raise on failure. `exp summary()`: compute exc over the real atoms. |
| D3ag-16 | minor | Novelty is overstated: the cancellation is the classical log-coordinate Itô mechanism; D*'s branches are classical quadratic and depth-3 penalties; "every colour" means two colours; the null-kernel claim is refuted. | D44 Status, "New:". |
| D3ag-17 | cosmetic | The renormalisation argument in (i)'s proof (see Verdicts). "M (with Σβ = m)" is also used for free masses. "The sparsest fit is representable": basis pursuit's representability suffices. The null kernel is compared with "the truth's 0.485" instead of the free end's 0.468. "nan" E/m rows. | D44, `report()`. |

## Lesson

A potential describes one update. Here the proposition proves its minimiser theorems for the flat √β update, but runs its numbers on
the paradigm's log β update. That update's own Itô term removes exactly the share of the potential that kills masses. What remains is a
drift that is not a gradient: for free masses it settles elsewhere, and for the conserved law it is dominated by the centring's cross
term. Medians over instances where every method agrees (basis pursuit equal to the truth) hid the fact that theorems and figures describe
different dynamics. The one simulation of the actual scheme was run at σ²τ_c ≈ 1, where its reduced flow is not the leading term. The
checks are cheap and should come before any claim: a one-second loop integral shows whether the dynamics measured is a gradient flow of
the potential; per-instance outcomes should be reported on the instances where the candidate answers differ; and a reduced flow should be
validated in the limit it claims (σ²τ_c/h → 0), against the projected state Φ(x) rather than the raw mean.

Sources: [Mertikopoulos–Moustakas, arXiv:0906.2094](https://arxiv.org/pdf/0906.2094);
[Mertikopoulos–Staudigl, stochastic mirror descent dynamics, arXiv:1710.01551](https://arxiv.org/pdf/1710.01551);
[Krichene–Bartlett (cites Raginsky–Bouvrie), arXiv:1707.06219](https://arxiv.org/pdf/1707.06219);
[Pillaud-Vivien–Reygner–Flammarion, arXiv:2206.09841](https://arxiv.org/pdf/2206.09841); [Pesme et al., arXiv:2106.09524](https://arxiv.org/pdf/2106.09524);
[Damian–Ma–Lee, arXiv:2106.06530](https://arxiv.org/pdf/2106.06530);
[Chen–Kunin–Yamamura–Ganguli, stochastic collapse, arXiv:2306.04251](https://arxiv.org/html/2306.04251v3).
