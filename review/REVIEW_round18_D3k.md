# Referee report D3k (round 18) — PARADIGM §5d: Proposition D‴.24 ((i)–(vi), proof, Status, *Numerics*) and Remark D‴.25; the §0 and §7 rows "| D‴.24 |", the §8 bullets "Done (Proposition D‴.24(i))" (branch 7) and "Done (Proposition D‴.24)" (branch 8, with its *Open*), the branch-8 edits of Theorem D‴.23's *Open*, the pointer in Remark D‴.10; `results/RESULTS_PARADIGM.md` (Proposition D‴.24); `src/core/floor_bias.py`; the nearest-atom option of `src/core/adaptive_prior.py`; `tests/test_adaptive_floor.py`; `experiments/exp_adaptive_floor.py` → `results/chaos/adaptive_floor.md` (+ `adaptive_floor_raw.npz`); LITERATURE_TRANSPORT [S51] and the D‴.24 audit bullet; README's D‴.24 sentence

Branch `claude/festive-bardeen-eopej2` at d2092a5 (HEAD; `git log --oneline -5`: d2092a5, 4e2e938, b76b3b9, 1aeaa56, 1fcff54). At the start
of the review the working tree had five untracked files of a concurrent session (`experiments/exp_adaptive_floor_2d.py`,
`experiments/exp_adaptive_floquet.py`, `src/core/adaptive_floquet.py`, `src/core/planar_loop.py`, `tests/test_planar_loop.py`) and no change
to tracked files; none of them was reviewed. During the review that session also modified five tracked files: `src/core/adaptive_prior.py`
(an optional `box` for a later Proposition D‴.27; with `box=None` the drift and the integrator are unchanged), `experiments/exp_adaptive_prior.py`,
`results/chaos/adaptive_prior.md` and its raw file, and `tests/test_adaptive_prior.py`. These edits are outside this review, which
covers d2092a5's versions. The material was reviewed from a snapshot of d2092a5 taken with `git archive` (src,
experiments, tests, the raw files and the md's) into the scratch directory, and the committed PARADIGM, RESULTS and LITERATURE were
extracted with `git show d2092a5:…`. Line numbers below are those of d2092a5's `PARADIGM.md` (D‴.24 l. 2799–2917, Remark D‴.25
l. 2919–2937, §0 row l. 120, §7 row l. 5368, §8 bullets l. 5570–5571 and 5609–5614, D‴.10 pointer l. 1761–1762), `RESULTS_PARADIGM.md`
l. 547–570 and `LITERATURE_TRANSPORT.md` l. 71 and 367–376. No repository file was modified: every script writes only to the scratch
directory, Python ran with `PYTHONDONTWRITEBYTECODE=1` (one thread), pytest with `-p no:cacheprovider`, and `--from-raw` ran on the snapshot
with `AF_OUT` pointing to a scratch folder. Scratch scripts and their `.out` files are in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3k/`:

| script | what it checks |
|---|---|
| `from_raw.out` | the `--from-raw` regeneration (byte-identical, md5 `55347aee…`) |
| `inspect_raw.py`, `recompute_md.py` | the raw file's structure; every number of the md, the paragraph, the §0/§7 rows and RESULTS, recomputed independently of `report()`; the counts of all 48 Part B systems and of Part C |
| `mymodel.py`, `check_model.py` | an independent implementation of the instance (own features, analytic first and second derivatives, the weighted stationary map, Newton, exact representations, H_L, R^G); checked against finite differences and against the repository's instance |
| `check_part_i.py` | (i): F by my central differences (three steps), D‴.9's F by the implicit function theorem, the bias form, the law with the full K×K C̃ and on the K − 1 zero-sum coordinates, the mode decomposition, the shares with and without the weight; (iii)'s split |
| `check_sampled.py` | at the 16 sampled weight vectors: Θ*(w), its D‴.5 splitting margins, the bias formula's accuracy; per-seed predictions of ‖f̂_m − f̂_M‖² and of the synchronous coupled distance D against the stored finals |
| `check_series.py`, `check_rates.py` | Part B time series (the split runs' drift; monotonicity and log-slopes of the adaptive runs; the boundary sample's s·h2); h₁(w), λ₁⁻(w) at every sample; when the adaptive distance falls below the fixed floor |
| `check_boundary.py` | (iv): the K₀-atom system at (8, 4, 3)/15 under the adaptive prior (η = ∞, σ₀), LSODA to s = 10⁶ |
| `check_prob.py`, `check_mc_law.py` | (iii): exact multinomial probabilities, Hoeffding's bound, permuted representations, the D‴.5 saddle probability of the fixed prior, the adaptive prior's limiting error at inadmissible samples; Monte Carlo (600 nested samples per m) of the exact floors against the laws |
| `check_W2.py`, `check_ot_D.py` | W² split into cross-cluster transport and within-cluster part; D and W² at the end of Part C for all four designs |
| `check_splits.py`, `check_partC.py`, `check_partC_split.py`, `check_d23_same.py` | the split fixed-prior finals; Part C per seed and its margins; Part C's OT rows against Theorem D‴.23's raw data |
| `check_crossings.py`, `check_empty.py` | (ii) beyond the cells: 𝒥_∘ across crossings, an empty cell |
| `check_partD.py` | Part D: the slaved-error ratio (filter factor, parameter offset); the "exact" run at dt = 10⁻³ and 5·10⁻⁴ (the experiment's own `real_loop`) |
| `check_repro.py` | Part B's records at s = 0, 50 re-run from the seeds; `AdaptivePriorLoop.run` against the lockstep integrator |
| `check_tests.py`, `tests.out`, `collect.out`, `check_fb_dirs.py` | the tests' asserted quantities and stronger variants; 5 passed (0.53 s); 122 tests collected; `bias_floor_law` with explicit directions |

The hour-long experiment was not rerun: every number below comes from the raw file, from re-integrations of single systems, or from the
scripts listed.

## Verdicts

| object | verdict |
|---|---|
| **(i)** bias formula at w | **holds** (first order in x⁰). Re-derived: the stationarity of Θ*(w) linearised at ϑ^e(w) is H_L(w)y + σ₀(y − x⁰(w)) = 0 in y_k = √w_kG^{−1/2}u_k, and f̂ − f = Ĵy; the sign of x⁰ is right (the bias points toward the prior). At the 15 representable sampled weights the formula's relative error is 0.004–0.154 and grows like ≈ 0.08‖x⁰(w)‖ (‖x⁰‖ = 0.048–1.55): the remainder is second order. |
| (i) F at q | **holds, and is exact** (not only first order): x⁰(q) = 0 kills the derivatives of Ĵ, H_L, √w and of the O(‖x⁰‖²) remainder. D‴.9's F by the implicit function theorem equals Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ to 8.8·10⁻¹⁵, 2.7·10⁻¹⁵, 1.1·10⁻¹⁵ (σ₀/σ_crit = 0.05, 0.15, 0.5). My central differences give an O(h²) discrepancy (6.2·10⁻⁵, 6.2·10⁻⁷, 6.2·10⁻⁹ for h = 10⁻³, 10⁻⁴, 10⁻⁵ at 0.05σ_crit), so the quoted 1.1·10⁻⁶–1.2·10⁻⁸ are finite-difference errors (D3k-7). |
| (i) law (B) | **holds as D‴.9's law written on the eigenvectors**: with the full K × K C̃ and raw directions, with the K − 1 zero-sum coordinates (C̃ restricted to w₁, w₂), and as the mode sum, it is 1.355928·10⁻², 6.166158·10⁻², 0.2291306 in all three forms. ‖Ĵv_j‖² = h_j and the Ĵv_j are orthogonal (Gram − diag(h): 4.2·10⁻¹⁶, with the 1/N normalisation of ⟨·⟩_ρ). The weight is ≤ σ₀/4 ((h − σ₀)² ≥ 0). But: the law is conditional on D‴.5-admissibility, which fails for 5 of the 12 sampled m-systems (D3k-1, **major**); its value is far from the expectation at m ≤ 60 (D3k-3); "carried by the modes excited at the level of the leak" is mostly a property of ∂_wx⁰ (D3k-5). |
| **(ii)** | **holds inside the cells.** min_{q′} W²(μ, Σq′_jδ_{z_j}) = ∫min_j‖ϑ − z_j‖²dμ is exact (test 2), and the flow is −mG∇_ϑ𝒥_∘, −(η/(σ₀w_j))G∇_{z_j}𝒥_∘ with the cell masses (test 1: residual 6.4·10⁻¹¹; with q in the metric 8.2·10⁻⁴). Across crossings 𝒥_∘ stays non-increasing (vertex selection; 600 of 600 decreases in four runs). An empty cell has a degenerate metric, and its atom is frozen for good. "Learned by the same energy" means minimised out, and the conservation sentence is a tautology (D3k-8). |
| **(iii)** stationary points | **holds**: at a replicated configuration with z_j = T̂_j the leak vanishes, and ϑ^e(w) with the prior at its atoms is stationary (test 4, drift < 10⁻¹²). The split S_y = H_L(H_L + σ₀I)⁻¹S^e_y is an exact derivative identity at q (my finite differences: 1.45·10⁻⁹, 9.2·10⁻¹⁰, 2.7·10⁻¹⁰). The parameter-level claim (D larger under the adaptive prior) holds dynamically: 1.71–4.82 per sample from the stored finals, and 1.83–2.73 in the Monte Carlo means. The numerics never test this, and the column of Remark D‴.25 mixes D with W² (D3k-2, **major**). |
| (iii) inadmissible samples | **holds asymptotically, overstated on the instance.** Hoeffding's bound is vacuous at m = 15 (1.87). The adaptive prior's expected limiting floor is 1.8 times the fixed prior's law at m = 15. "No exact representation exists" is false under permutations. The parameters run off, so the parameter-level floor is a conditional notion (D3k-4). The W² parenthesis holds: the cross-cluster part is ½‖w − w′‖₁ × 34.5–49.0. |
| **(iv)** | **holds** (the boundary case formal, and confirmed with its constant). With a = a_max tanh(ã/a_max), a_max − a ≈ 2a_max e^{−2ã/a_max} gives a loss ∝ e^{−4ã/a_max}, ã ≈ (a_max/4)ln s and a loss ∝ 1/s. The loss's curvature cancels, so ‖f̂ − f‖² ≈ a_max²w_k(1 + σ₀/η)/(8G_αα s) = 0.4/s (η = ∞) and 0.8/s (η = σ₀). The K₀-atom run gives s·d = 0.4016 and 0.8054 at s = 10⁶, with α₃ − ln s → −0.33 and −1.02. The m = 15 run gives s·h2 = 0.4225 and 0.877 at s = 6000, with local exponents −0.999 and −1.009. The two runs agree to 1% at s = 6000. |
| **(v)** and **Remark D‴.25** | **holds in part.** The bias, transit, collective-rate and a-priori-bound columns and the entropic row are supported by statements and numbers; the entropic row's a-priori bound is D‴.20's T^ε in conv supp ν. The "parameter-level floor" column is inconsistent: under D the OT rows are not 0, and under W² the nearest-atom rows are not D‴.5's law. The Reading's "the mass constraint removes both" therefore fails at the coupling level (D3k-2, **major**). "Resolves Remark D‴.10's trade-off at the function level" is supported at the function level; the price at the parameter level is understated (D3k-9). |
| **(vi)** | **holds.** f̂ = f on the reference orbit (supp ρ of the occupation measure) is all the invariance needs; the instance's representation is exact to 7.8·10⁻¹⁶ on [−2, 2]. At the invariant set the cells and the plan refresh are inert. The 10⁻¹¹ level is RK4 truncation (8.09·10⁻¹² → 5.42·10⁻¹³ when dt halves, ratio 14.9). The fixed prior's ratio is explained to 0.07%: the error dynamics' filter 1/(1 + n²/k_e²) (0.9930, 0.9913) times the real loop's parameter offset from Θ*(w) (1.0089, 1.0050) gives 1.0018 and 0.9962, against 1.0011 and 0.9956 measured. But the m = 15 fixed run sits on a saddle orbit (D3k-1). |
| Proof | **holds** for (i), (ii) (inside the cells), (iii), (iv) (as D‴.23(ii)) and (vi). |
| *Status* | **mostly accurate.** (i)'s exactness is misplaced (D3k-7), the §0 row's "proved" should read "formal" for (B), and the novelty is overstated against Chen–Rotskoff–Bruna–Vanden-Eijnden 2020 and the Tikhonov filter factors (D3k-6). |
| *Numerics* A | **numbers hold** (all rows, margins 0.025–0.10, the boundary sample exactly on the boundary: 3/15 and 0.8/4 are the same double, margin 0.0). "The floor's expectation" mislabels the law's value (D3k-3). |
| *Numerics* B | **numbers hold; the adaptive half confirmed, the fixed half misdescribed.** The counts are conserved in all 48 systems. The adaptive distances fall monotonically: no record-to-record increase, except at the 10⁻¹⁴ level in one run. At η = ∞ the late log-slopes match the rates 2λ₁⁻ of the two systems (mostly between them, at most 16% above); at η = σ₀ they are still pre-asymptotic. The adaptive finals are replicated and near ϑ^e(w). Five of the 12 fixed-prior systems are at D‴.5 saddles and split, and the W² example compares with one of them (D3k-1, **major**). W² is not the gauge's floor (D3k-2). |
| *Numerics* C | **numbers hold.** The OT rows reproduce Theorem D‴.23's settled-gate runs to all digits (12 of 12 fit errors, same t_open). The nearest/fixed row holds the same five D‴.5-inadmissible samples (margins −0.38 to −5.76·10⁻³ at the misplaced prior; 4 split by s = 3000). The ranking depends on the horizon (D3k-9). |
| *Numerics* D | **numbers hold**, with the windows and ratios recomputed; the explanation is under (vi). |
| `src/core/floor_bias.py` | **holds** (Newton with continuation, Ĵ, H_L, the decomposition). The custom-direction branch of `bias_floor_law` is correct only for orthonormal directions; it is latent, since no caller uses it (D3k-10). |
| nearest-atom option of `adaptive_prior.py` | **holds**: cell plan, cell masses, empty cell frozen, one Lloyd step per refresh at η = ∞. `AdaptivePriorLoop.run` equals Part B's lockstep integrator (difference 1.2·10⁻⁷, the float32 rounding of the stored final). |
| `exp_adaptive_floor.py` → md | **holds**: regeneration is byte-identical, and Part B's records at s = 0, 50 reproduce bitwise from the seeds for η = 0, σ₀, ∞. `report()` prints neither the fixed prior's D‴.5 margins nor D (D3k-1, D3k-2). |
| `tests/test_adaptive_floor.py` | **pass (5/5, 0.53 s); none vacuous**, but test 3's mode assertion passes without the weight and test 4's sample is a fixed-prior saddle (D3k-11). README's "122 tests" is right. |
| §0, §7 rows, §8 bullets, D‴.10 pointer, RESULTS, LITERATURE, README | **numbers hold**; wording under D3k-1–D3k-10. |
| Consistency with Theorem D‴.23 / round 17 | **D‴.24 is right, at the level of stationary points.** Round 17 upheld "with learned masses the floor returns" because the free-mass plan conserves the cluster masses. That holds for the parameter level, where D is in fact larger (1.7–4.8×). It fails for the function level: with adaptive atoms the replicated stationary point at every w ∈ W is the exact representation, and the dynamics reaches it at rate 2λ₁⁻(w). The new wording ("does not return at the admissible samples … what remains is the parameter-level floor of the gauge") is exact up to D3k-4: the instance's probabilities, and a parameter-level floor that is conditional on admissibility. It is consistent with D‴.23's revised (i) (replicated stationary points, boundedness a hypothesis) and (v) (sketched; D‴.24's Status does not claim attractivity). |

## Confirmed issues

**D3k-1 (major; D‴.24 *Numerics* Part B l. 2902–2907, Part C l. 2909–2912, Part D l. 2913–2917; Status l. 2888–2890 and test 4;
RESULTS l. 556–565 and the tests line l. 567–570; the experiment's docstring and `report()`; the use of "admissible" in (iii), (v),
Remark D‴.25 and the rows).**
In 5 of the 12 (seed, m) pairs the fixed prior's side of the comparison sits at weights where its replicated stationary point is a saddle
(Proposition D‴.5's W_adm fails). Those systems split, and the text does not say so. Meanwhile "admissible" is used throughout D‴.24 for the
other notion, w ∈ W (an exact representation exists).

(a) *The margins.* I computed σ₀ − max_k λ_max(R^G_k) at Θ*(w), the fixed matching prior, σ₀ = 0.15σ_crit, with my own Newton (residuals
≤ 8·10⁻¹⁶; `check_sampled.py`):
- seed 0, m = 15, (6, 5, 4): −3.48·10⁻³ (round 14's D3g-1 found −0.00348 for this vector);
- seed 1, m = 15, (6, 3, 6): −1.09·10⁻³;
- seed 1, m = 60, (24, 15, 21): −2.16·10⁻³;
- seed 3, m = 15, (8, 4, 3): −6.14·10⁻³;
- seed 3, m = 60, (25, 17, 18): −3.98·10⁻³.

Each negative margin is at class 1, with w₁ ≥ 0.4. The other seven m-systems and all four M = 510 references have margins ≥ +3.79·10⁻³.

(b) *The runs.* Exactly these five Part B fixed-prior systems end split; every other m-system and all four references end replicated
(spread 0.0):
- the maximum within-class spread at s = 6000 is 0.92, 0.42, 0.40, 0.83 and 0.60 (G⁻¹ metric);
- class 1 is broken into sub-clusters (5, 1), (3, 1, 1, 1), (10, 6, 2, 2, 2, 1, 1), (2, 2, 2, 1, 1) and (15, 8, 1, 1) (single linkage at
  0.05);
- their ‖f̂_m − f̂_M‖² is 0.872, 0.973, 0.888, 1.010 and 0.811 times that of the replicated Θ*(w), against 1.0000 at all seven admissible
  samples;
- their h2 drifts until s ≈ 3000: h2(600)/h2(6000) = 1.10, 1.03, 1.08, 1.00, 1.03, against ≤ 1.0006 for the others.

(c) *The paragraph's W² example* ("5.92 and 5.92 at seed 0, m = 15, above the fixed prior's (4.48)"; RESULTS "5.92 against 4.48")
compares the exact representations with a split configuration. The replicated Θ*(w) has W² = 5.21.

(d) *Parts C and D, test 4.*
- Part C (misplaced prior) contains the same five samples, with margins −2.87, −0.38, −1.52, −5.76 and −3.41·10⁻³. Four have split by
  s = 3000 (spread 0.63–1.1); the fifth, at −0.38·10⁻³, has partly split (0.14). All 36 adaptive m-systems of Parts B and C end
  replicated (spread ≤ 1.2·10⁻⁷).
- Part D's m = 15 fixed-prior run and test 4 use (6, 5, 4)/15. The real-loop run sits on a saddle orbit (Proposition D‴.15), held there only
  by the exact replication of its start.
- Test 4's "the fixed prior's stationary point there is biased" refers to a configuration the dynamics leaves.

(e) *This is not rare on the instance.* Summing the multinomial probabilities of the count vectors with a non-positive margin, the fixed
matching prior's replicated point is a saddle with probability 0.376, 0.260 and 0.109 at m = 15, 60 and 240. That is comparable with the
adaptive prior's non-representability: 0.361 + 0.118 (boundary), 0.074 + 0.053, and 0.0020 + 0.0012.
- The law (B), like D‴.9's, holds only on D‴.5-admissible samples.
- The comparison "fixed: the floor (B); adaptive: 0 at admissible w" (v, Remark D‴.25, §0) uses two different admissible sets without
  saying so.
- Conversely, at an exact representation R = 0, so the separations relax at σ₀ for every w ∈ W. The adaptive prior therefore also removes
  this splitting instability, an advantage the text does not claim.

This repeats the lesson the log records after rounds 14 and 16: "every reading of runs at sampled weights must print the margins at those
weights".

*Fix:*
- `report()` prints, per sample, the D‴.5 margin of the fixed prior's Θ*(w) and the end spread.
- The paragraph and RESULTS say "5 of the 12 fixed-prior systems are at weights where Θ*(w) is a saddle (class 1, w₁ ≥ 0.4) and split;
  their floors are 0.81–1.01 of the replicated point's".
- Either quote the W² example against 5.21 or pick an admissible seed.
- Distinguish "D‴.5-admissible (w ∈ W_adm of the fixed prior)" from "representable (w ∈ W)" throughout D‴.24, Remark D‴.25 and the rows.
- Use an admissible sample in test 4, or assert its margin's sign.
- Note in Part D that the m = 15 fixed-prior orbit is a saddle.
- State that at an exact representation the splitting margin is σ₀ for every w ∈ W.

**D3k-2 (major; Remark D‴.25's column "parameter-level floor" and its Reading l. 2924–2937; (v) l. 2852–2855; *Numerics* l. 2905–2907
"W²(μ_m, μ_M) at the end then equals … the parameter-level floor of the gauge"; §0 row "while W² reaches the exact representations' (the
gauge's floor)"; RESULTS l. 558–559; the Part B docstring and md header).** The column mixes two distances:
- the nearest-atom rows quote Proposition D‴.5's law, which is the synchronous coupled distance D of the shared particles;
- the OT rows' "0" is the measure-level W² (at commensurate m the stationary measures coincide).

The column is right in neither distance. I recomputed both from the stored finals (`check_ot_D.py`, `check_W2.py`).

*Under D:*
- The OT leak, fixed and adaptive, ends Part C with D = 2.6–9.4 (m = 15), 2.1–5.4 (m = 60) and 0.19–1.47 (m = 240). The transit relabels
  different particles in the m- and M-systems, which is D‴.12's own *Open*.
- The nearest-atom leak ends with D = 0.13–2.0, 0.055–0.45 and 0.002–0.032 (fixed), and 0.26–5.3, 0.09–0.92 and 0.004–0.066 (adaptive).
- At m = 240 the OT designs have the largest D of the four.

*Under W²:*
- The nearest-atom rows are not D‴.5's law. W² is the cross-cluster transport ½‖w^m − w^M‖₁ × (34.5–49.0), 67–99% of W² in all 23
  computable cases (84–99% fixed, 67–99% adaptive), plus a within-cluster part of 0.002–1.5 (fixed) or 0.004–2.4 (adaptive).
- The leading term is common to both priors and of order ‖δw‖₁ ~ m^{−1/2}.
- So Part B's W² comparison is dominated by that common term: the adaptive/fixed W² ratios are 1.00–1.32 at the representable samples. At
  the boundary sample they are 1.67 and 1.85, because its parameters run off (D3k-4).

Consequences:
- The Reading's "the mass constraint removes both [parts] by fixing the weights … [the adaptive prior] enlarges the second" holds only at
  the measure level. There the enlargement is a second-order, within-cluster effect.
- At the coupling level the mass constraint does not remove the floor; it replaces it by the transit's mismatch.
- The numerics' "parameter-level floor of the gauge" is this W², not the gauge's floor.

*What holds* is (iii)'s claim about D, which the paper never tests dynamically:
- under the adaptive prior, D from the stored finals equals the exact representations' D to 0.95–1.00;
- under the fixed prior it equals the replicated Θ*(w)'s D to 0.94–1.00 (1.0000 at the admissible samples);
- the adaptive/fixed ratio is 1.71–4.82 per sample, and the Monte Carlo means are 1.83, 2.62 and 2.73 at m = 15, 60 and 240, against the
  law's 2.17.

*Fix:* split the column in two and rewrite the Reading accordingly.
- *Measure level (W²):* for the nearest-atom leak under both priors, ½‖δw‖₁Δ² + O(‖δw‖²), the within-cluster gauge part being larger under
  the adaptive prior; for the OT leak at commensurate m, 0.
- *Coupling level (D, D‴.5):* D‴.5's law with H_L(H_L + σ₀I)⁻¹S^e (fixed) or S^e (adaptive); for the OT leak, the transit's mismatch.
- In the numerics, the §0 row and RESULTS, call W² "the measure-level distance, of order ‖δw‖₁ (the cross-cluster transport)", and report
  D per seed as the test of (iii).

**D3k-3 (minor; *Numerics* l. 2898–2899 "the floor's expectation is 3.99·10⁻³, 9.07·10⁻⁴, 1.36·10⁻⁴ at m = 15, 60, 240"; (i)
l. 2817 "the floor is at most (1/m − 1/M)(σ₀/4)tr(…), whatever the excitation spectrum"; the md's Part A header).** These are the
delta-method values. A Monte Carlo over 600 nested multinomial samples (M = 510, the m-system taking the first m labels), with the exact
replicated stationary points, gives E‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖² = 3.08 ± 0.23, 1.41 ± 0.11 and 1.06 ± 0.05 times the law. Conditionally on
D‴.5-admissibility of both samples the factors are 3.50 ± 0.32, 1.43 ± 0.16 and 1.02 ± 0.06.

So the o(1/m) remainder dominates at m = 15 and is about 40% at m = 60. The reason is that the gauge's derivative,
∂α^e/∂w = −A/(w²(1 − (A/(wa_max))²)), varies fast near the saturation: the matching amplitude is 3 against a_max = 4.

Per sample:
- The exact bias at the sampled weights reproduces the observed floor at all seven admissible samples (ratio 1.0000). This is the strongest
  numerical support for "the floor is the bias at the sampled weights", and the paper does not show it.
- The delta method F(w^m − w^M) misses the observed floor by factors 0.86–1.72 at the seven admissible samples (0.42–1.31 at the five
  split ones); at m = 240, seed 0, 8.17·10⁻⁵ is observed against 4.76·10⁻⁵ predicted.
- "The floor is at most …" bounds the law, not the floor: at m = 15 the expectation is 2.3 times the bound.

*Fix:*
- Write "the law's value is", and "the law (B) is at most".
- Add the per-sample comparison (exact bias at the sampled weights, ratio 1.0000) and the Monte Carlo ratios.

**D3k-4 (minor; (iii) l. 2834–2838; (vi) l. 2863; §0 row "no function-level floor at any admissible m"; §8 branch-8 bullet
l. 5610–5613; the audit bullet; the commit message).** The paragraph on inadmissible samples has five problems.

(a) *The bound.* Hoeffding's Σ_k exp(−2m(q_k − |A_k|/a_max)²) is 1.870, 0.602 and 0.036 at m = 15, 60 and 240, so it is vacuous at m = 15.
The exact probabilities are P(inadmissible) = 0.361, 0.074 and 0.0020, plus P(boundary) = 0.118, 0.053 and 0.0012.

(b) *The expected floor on the instance.* At an inadmissible count vector the adaptive prior's limit is the constrained minimum of the
weighted loss with |a| ≤ a_max: for example 4.87·10⁻³ at (3, 6, 6)/15 and 5.14·10⁻⁴ at (14, 23, 23)/60. Summing over the count vectors of
probability > 10⁻⁵ (> 10⁻¹⁰ at m = 240), so as lower bounds, the expected limiting floor is 7.2·10⁻³, 2.0·10⁻⁴ and 6.6·10⁻⁷ at m = 15, 60
and 240:
- 1.8, 0.22 and 0.005 times the fixed prior's law;
- 0.59, 0.16 and 0.005 times its Monte Carlo expectation (D3k-3).

"Exponentially small rather than of order 1/m" is therefore asymptotic, not in force at m ≤ 60 on the instance. It also needs every
representable sample to converge to its exact representation. That is numerical only (23 runs, plus the boundary run at 1/s), and (iii) is
"proved for the stationary points".

(c) *Existence under permutations.* "No exact representation exists" at inadmissible weights is false for 31%, 87% and 100% of the
inadmissible mass at m = 15, 60 and 240, because a permutation of the classes over the bumps is admissible: for (14, 23, 23)/60, class 2 can
take bump 1. What fails is the branch continuous from ϑ^e(q), which the dynamics from the warm clouds follows.

(d) *Runaway parameters.* At inadmissible and boundary samples the amplitude coordinate grows like (a_max/2) ln s and (a_max/4) ln s
respectively, so the parameter-level distance diverges like (ln s)².
- The boundary run has D = 6.38 (η = σ₀) and 7.80 (η = ∞) at s = 6000.
- Its K₀-atom continuation has α₃ = 13.49 at s = 10⁶, against 8.31 at s = 6000.
- Since P(inadmissible ∪ boundary) > 0 at every m, the unconditional expected parameter-level distance has no finite limit at any m. "The
  parameter-level floor of the gauge" is a conditional statement, like D‴.5's law.

(e) *Wording.* "At any admissible m" (§0) and "perfect tracking at every admissible m" ((vi)) make admissibility a property of m. It is a
property of the sampled w.

*Fix:*
- Give the instance's exact probabilities next to Hoeffding's bound.
- Write "if every representable sample converges to its exact representation, the expected limiting floor is exponentially small in m",
  and add the instance's values next to the fixed prior's.
- Write "no exact representation on the branch of ϑ^e(q)".
- Write "the parameter-level floor, conditional on representability (elsewhere the parameters run off)".
- Write "at every admissible sample".

**D3k-5 (minor; (i) l. 2818–2819 "it is carried by the modes excited at the level of the leak"; §0 row; §8 branch-7 bullet
l. 5570–5571; *Numerics* "the modes with h_j/σ₀ between 0.28 and 4.6 carry 99% of it"; Status "carried by the modes within a factor 5 of
σ₀"; test 3's last assertion).** The concentration belongs mostly to the gauge's displacement, not to the weight:
- The shares of z_jᵀC̃z_j alone (the bound's distribution, which does not depend on σ₀) already put 97.2% on modes 2–6.
- With the weight, modes 2–6 hold 99.4%, 99.0% and 97.9% at σ₀/σ_crit = 0.05, 0.15 and 0.5, while their h_j/σ₀ range over 0.83–13.8,
  0.28–4.6 and 0.08–1.38.
- The weight is a broad band-pass: it is ≥ half its peak for h/σ₀ ∈ [3 − 2√2, 3 + 2√2] = [0.17, 5.83], and at 0.15σ_crit it is 0.16–0.87 of
  its maximum on modes 1–6.
- At 0.05σ_crit the modes within a factor 5 of σ₀ (modes 2, 3) carry only 42%.
- The weight's own effect is to move shares toward h ≈ σ₀: modes 2–3 go from 21% unweighted to 42% (0.05σ_crit) and 10% (0.5σ_crit).
- Test 3's "share[near] > 0.8" passes with the weight removed (0.972).

*Fix:*
- Write "the weight is a band-pass centred at h = σ₀, of half-width a factor 5.8; on the instance the floor sits in the modes that carry the
  gauge's displacement, re-weighted toward h ≈ σ₀".
- In test 3, compare the weighted and unweighted shares at two values of σ₀.

**D3k-6 (minor; Status l. 2882–2886 "New is the identification … and the two ways to remove it"; §0 last column; LITERATURE [S51] l. 71
and the audit bullet l. 367–376).**
- *Chen–Rotskoff–Bruna–Vanden-Eijnden 2020.* [S51]'s annotation omits the result closest to (iii). Two independent web searches returned
  the abstract's statement that if the mean-field dynamics converges to a measure that interpolates the training data, the asymptotic
  deviation eventually vanishes in the CLT scaling (I could not open the full text: arxiv.org is blocked by the session's proxy). The
  adaptive prior removes the function-level floor by making the regularised flow's replicated stationary points unregularised interpolants,
  which is the same mechanism. The novelty claim should be delimited against it: what is new is that a leak toward a learned prior restores
  this property under regularisation, together with the bias formulation.
- *The weight and its bound.* The universal weight is the squared Tikhonov filter factor (σ₀/(h + σ₀))² times the mode's visibility h. Its
  bound is the AM–GM bound h/(h + λ)² ≤ 1/(4λ), familiar from ridge regression's filter factors. D‴.23's Status already cites the filter
  factors.
- *The derivation.* F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ is the chain rule applied to D‴.23(iii)'s first-order bias at x⁰(q) = 0.
- *Why the adaptive prior has no function-level floor.* At an exact representation the unregularised fit is weight-invariant by definition
  (f̂ = f at every w ∈ W), so the absence of a function-level floor at the adaptive prior's replicated stationary points follows at once
  from (ii)–(iii). The substantive new content is (i)'s identification and (iii)'s split.

*Fix:*
- Annotate [S51] with Chen et al.'s interpolation result.
- Restrict "new" accordingly.
- Cite the filter-factor bound as classical.

**D3k-7 (cosmetic; Status l. 2877–2879; §0 row "(i) proved … (the law with Proposition D‴.9's delta-method status)"; *Numerics*
l. 2895–2896).** The Status says "the identity F = ∂_w(f̂_{Θ*(w)} − f) is exact, f not depending on w". That identity is trivial. The identity
of interest, F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ at q, is also exact (verdict table), while the bias formula at w ≠ q is first order.
- The numerics' "F equals the bias form to 1.1·10⁻⁶, 1.9·10⁻⁷, 1.2·10⁻⁸" are finite-difference errors, not approximation errors: my central
  differences give O(h²), and the repository uses h = 10⁻⁴.
- D‴.9(iii)'s status is "formal", so "(i) proved" in the §0 row should say "(B) formal".

*Fix:* write "F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ exact at q, the bias formula first order, the law (B) formal (delta method, as D‴.9)", and in the
numerics write "to the finite-difference error".

**D3k-8 (cosmetic; (ii) l. 2821–2829, proof (ii), the `AdaptivePriorLoop` docstring).**

(a) "With the prior's masses learned by the same energy" describes the semi-relaxed minimum over the masses, the infinitely fast limit of a
mass dynamics, not a learning dynamics. Say "minimised out".

(b) "The cell counts are conserved as long as no particle crosses a boundary between the moving cells" is a tautology.

(c) The metric (σ₀w_j/η)G⁻¹ degenerates at an empty cell. The code freezes that atom, and the cell then stays empty for good: from a cloud
near two atoms, with the third atom 6 G⁻¹-units away, both η = σ₀ and η = ∞ end with counts (15, 15, 0), the atom unmoved, fit error
0.73 → 0.59.

(d) Across crossings the flow is the vertex-selection flow of 𝒥_∘ in a branch-dependent metric, and 𝒥_∘ is non-increasing, because the
active branch's derivative bounds the Dini derivative of the minimum. Confirmed in four spread-cloud runs: 600 of 600 decreases between
refreshes, one count change. This is worth stating, as D‴.23(i) does for the OT version.

*Fix:* reword (a) and (b); state (c) and (d).

**D3k-9 (cosmetic; (v) l. 2855–2858; Remark D‴.10's pointer l. 1761–1762; the audit bullet's verdict; Part C's reading).**

*The trade-off.* "This resolves Remark D‴.10's trade-off … at the function level … the price is the rate" is supported at the function level:
- in all 24 adaptive runs of Part B, including the boundary sample, ‖f̂_m − f̂_M‖² is below the fixed prior's final floor by s = 100–150;
- at η = ∞ it then decays at the rates 2λ₁⁻ of the two systems: log-slopes on [5000, 6000] of 5.8–9.6·10⁻⁴, against 2λ₁⁻ = 5.8–8.2·10⁻⁴
  (mostly between the two systems' values; at most 16% above).

But D‴.10's trade-off concerned the floor "at both levels". The adaptive prior:
- enlarges the parameter-level floor (D 1.7–4.8× per sample);
- gives up the a-priori bound;
- has a slowest collective rate λ₁⁻(q) = 1.51·10⁻⁴ (η = σ₀) or 3.06·10⁻⁴ (η = ∞), 49 and 24 times below the fixed prior's
  h₁ + σ₀ = 7.38·10⁻³ and below the σ-modification's κ = 0.00196.

*Part C's ranking depends on the horizon.* Until the settled gate opens (s = 1211–1304) the OT/adaptive row equals the OT/fixed row by
construction (fit error 9.78·10⁻³).
- At s = 1250 the ungated nearest-atom adaptive prior leads, with fit errors 1.5·10⁻³–2.0·10⁻² (seeds 0, 1).
- The OT design leads from s ≈ 2000: 4.0–4.3·10⁻⁴, against 7.6·10⁻⁴–1.3·10⁻² for the nearest-atom adaptive prior.
- At s = 3000 the OT design is ahead despite its ~1750 units of adaptation, so the Reading survives.

*Fix:*
- In the D‴.10 pointer write "at the price of the collective rate and of a larger parameter-level floor".
- In Part C, mention the gate's delay.

**D3k-10 (cosmetic; rows, table entries, code).**
- The §7 row's "the parameter split to 10⁻⁷" should read 1.4·10⁻⁷ (the maximum).
- The §7 row's "e_rms ≤ 9·10⁻¹²": the window maximum is 9.03·10⁻¹²; the instantaneous |e| reaches 1.31·10⁻¹¹.
- The §0 row's "its stationary replicated points are the exact representations of the sampled weights" is wrong as worded: they are the
  critical points of the unregularised weighted loss, the exact representation among them. (iii) states this correctly.
- Remark D‴.25, entropic row: the function-level entry "D‴.5's law at Θ*_ε(q)" should be D‴.9's function-level law there.
- Remark D‴.25, nearest/fixed "transit": "none (counts conserved)" is D‴.9(i)'s hypothesis, not a theorem; round 14 found count changes from
  ρ = 4.5.
- Remark D‴.25, nearest/adaptive "collective rates": "λ_j^−(w) < min{h_j, η}" should read min{h_j(w), η}.
- (vi) needs f̂ = f on the reference orbit, i.e. on supp ρ of the occupation measure, not on a quadrature grid.
- `floor_bias.bias_floor_law(directions=…)` uses VdᵀC̃Vd, which is the coordinates' covariance only for orthonormal directions. Passing the
  default e_j − e_K explicitly gives 0.2490 instead of 0.06166; an orthonormal zero-sum basis gives 0.06166. No caller passes directions.
  Use Vd⁺C̃Vd⁺ᵀ, or document the restriction.

**D3k-11 (cosmetic; `tests/test_adaptive_floor.py`).** All five tests pass in 0.53 s and none is vacuous, but each tests less than its
name says:
- **Test 1** is discriminating (residuals 6.4·10⁻¹¹ and 1.6·10⁻¹¹ against drifts of 3.1·10⁻² and 3.3·10⁻³; 8.2·10⁻⁴ with q in place of the
  cell masses). It runs far from every tie (gap 37 between nearest and second-nearest squared distances), so only the interior of the cells
  is tested.
- **Test 2** uses all four trials, but only random Dirichlet masses (smallest excess 2.53). Perturbing the cell masses by ±1/m would test
  the minimum where it matters.
- **Test 3** has tolerances of 10⁻⁵ against actual errors of 1.9·10⁻⁷ and 3.1·10⁻⁸; "law < bound" is an algebraic identity; the near-σ₀
  share passes without the weight (D3k-5).
- **Test 4**'s sample is a fixed-prior saddle (margin −0.00348; D3k-1).
- **Test 5** is fine (9.2·10⁻⁸; 0.545; 2.17).

*Fix:*
- In test 2, perturb the cell masses by ±1/m.
- In test 3, compare against the unweighted shares.
- In test 4, use a D‴.5-admissible sample, or assert the margin.
- Add a crossing or empty-cell test for (ii).

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and reproducibility.**
   - `--from-raw` rewrites `adaptive_floor.md` byte-identically (md5 `55347aee…`).
   - Part B's records at s = 0 and 50 reproduce bitwise from the seeds for η = 0, σ₀ and ∞.
   - Part C's OT rows reproduce Theorem D‴.23's settled-gate runs to all digits (12 of 12 fit errors, the same t_open).
2. **Every number of the md, the paragraph, the §0 and §7 rows and RESULTS**, from the raw records independently of `report()`:
   - **Part A:** all three rows; the split 1.4·10⁻⁷–2.9·10⁻⁸; the laws 32.27, 21.32 and 8.79 against 46.32 (ratios 1.435, 2.173, 5.269);
     the 99% share of modes 2–6 at 0.15σ_crit (0.990; mode 2 has h/σ₀ = 0.278); the margins 0.025–0.10.
   - **Part B:** every range at s = 6000 (fixed 3.0·10⁻⁴–5.5·10⁻³, 1.1·10⁻⁴–9.9·10⁻⁴, 4.8·10⁻⁶–8.2·10⁻⁵; adaptive admissible maxima
     2.0·10⁻⁶ and 5.4·10⁻⁸ at m = 15); W² 5.92/5.92/4.48; the boundary values 1.5·10⁻⁴ and 7.0·10⁻⁵.
   - **Part C:** all twelve means and min–max ranges.
   - **Part D:** all windows and ratios (1.001, 0.987, 0.996, 0.992).
3. **The instance, independently.** My own implementation (analytic derivatives within 4.4·10⁻¹⁰ and 8.8·10⁻¹⁰ of finite differences)
   equals the repository's: grid and F to 0.0, f̂ to 2.2·10⁻¹⁶, Θ*((6, 5, 4)/15) to 4.4·10⁻¹⁶.
4. **(i) and (B):**
   - F at q is exact (implicit function theorem to 10⁻¹⁵; finite differences O(h²)).
   - (B) is D‴.9's law, identical with the full K × K C̃ and with the zero-sum coordinates.
   - The norms and orthogonality of Ĵv_j hold (4.2·10⁻¹⁶); the weight ≤ σ₀/4.
   - The bias formula at the sampled weights is first order, with relative error ≈ 0.08‖x⁰‖.
   - New identity: the bound equals (σ₀/4) × (the adaptive prior's coupled-distance law), since tr(∂_wx⁰C̃∂_wx⁰ᵀ) = 46.3249 = D_gauge. The
     function-level floor of the fixed prior is thus at most σ₀/4 times the parameter-level law of the adaptive prior, which is worth
     stating.
5. **(iii):**
   - The split is exact at q.
   - D under both priors from the finals: the adaptive prior's equals the exact representations' to 0.95–1.00, the fixed prior's equals the
     replicated points' to 0.94–1.00; the adaptive/fixed ratio is 1.71–4.82.
   - Monte Carlo of the parameter-level laws: fixed E[D]/law = 0.735, 0.863, 0.957; adaptive 0.619, 1.041, 1.203.
   - W² ≈ ½‖δw‖₁ × 34.5–49.0 plus the within-cluster part.
6. **Part B dynamics:**
   - The counts are conserved in all 48 systems.
   - The adaptive distances fall monotonically (one non-monotone run, at 10⁻¹⁴).
   - The late η = ∞ slopes match the rates 2λ₁⁻ of the two systems (at most 16% above). The 24 adaptive m-systems end replicated, and
     the 22 representable ones end near ϑ^e(w) (distances 6.8·10⁻⁴–6.0·10⁻² at η = σ₀, 1.7·10⁻⁴–1.2·10⁻² at η = ∞).
   - The adaptive distance is below the fixed floor by s = 100–150 in all 24 runs.
7. **(iv)** with its constant: K₀-atom integration to 10⁶ (s·d = 0.4016 and 0.8054; α₃ − ln s → −0.33 and −1.02); agreement with the
   m = 15 run to 1% at s = 6000.
8. **(vi) and Part D:** RK4 scaling of the invariant set's level (ratio 14.9 for dt halved); the representation is exact on [−2, 2]
   (7.8·10⁻¹⁶); the slaved-error ratio explained by the filter and the parameter offset to 0.07%.
9. **(ii):** the semi-relaxed identity; the gradient structure with the cell masses; monotonicity of 𝒥_∘ across crossings.
10. **Tests:** 5 passed (0.53 s); every asserted quantity recomputed; README's 122 tests correct (122 collected at d2092a5).
11. **Round 17's *Open* claim, adjudicated:**
    - with learned (cell) masses the parameter-level floor returns, and is larger;
    - the function-level floor does not, at the representable samples, in the stationary limit;
    - D‴.24's revised wording is consistent with D‴.23 as revised in round 17.

## Not checked

- **The hour-long experiment was not rerun.** Re-integrations were limited to:
  - single Part B systems for 50 slow units;
  - the K₀-atom boundary system to s = 10⁶;
  - the real loop's exact case for 100 time units at two values of dt;
  - spread-cloud and empty-cell runs of 600 and 1500 slow units for (ii).
- **Convergence and attractivity.** The convergence of every representable sample to its exact representation was not followed beyond the
  23 runs, and the attractivity of perfect tracking in the real loop was not examined.
- **Inadmissible samples.** Their limiting errors are constrained minima found by L-BFGS-B from the bumps. They were taken as the adaptive
  dynamics' limits, not integrated.
- **Monte Carlo.** The Monte Carlo used 600 samples per m (standard errors quoted).
- **Used as given (rounds 9, 10, 12, 17):** Proposition D‴.5, Proposition D‴.9's delta-method derivation, Theorem J, Proposition D‴.15, and
  Theorem D‴.23(ii) and (v).
- **Literature.**
  - arxiv.org is blocked by the session's proxy, so the full texts of [S51] and [30] were not read.
  - Chen–Rotskoff–Bruna–Vanden-Eijnden's interpolation result is taken from its abstract as returned by two web searches
    ([NeurIPS 2020 abstract page](https://proceedings.neurips.cc/paper/2020/hash/fc5b3186f1cf0daece964f78259b7ba0-Abstract.html),
    [arXiv:2008.09623](https://arxiv.org/abs/2008.09623v3)).
  - The Stochastic Processes and their Applications details of Sirignano–Spiliopoulos
    ([arXiv:1808.09372](https://arxiv.org/abs/1808.09372v2)) were not verified.
- **Concurrent work.** The concurrent session's untracked files (a two-dimensional variant and a Floquet study) are outside the material.
