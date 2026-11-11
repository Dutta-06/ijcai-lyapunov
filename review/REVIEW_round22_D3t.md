# Referee report D3t (round 22) — the round-21 revisions of Propositions D‴.30 and D‴.31: PARADIGM §5d (statements, proofs, Status, *Numerics*), the §0 and §7 rows "| D‴.30 |", "| D‴.31 |", the §8 bullets "Done (Proposition D‴.30)", "Done (Proposition D‴.31)", `results/RESULTS_PARADIGM.md`, LITERATURE [S57], [S58] and the audit bullets, PAPERS.md paper 12, `src/core/measurement_noise.py` (`zero_frequency_operators`, `_whitening`, `second_order_mean`, `mean_biases`, `separation_drift`), `exp_measurement_noise.py` and `exp_eiv_remedy.py` with their md and raw, the new tests

Branch `claude/festive-bardeen-eopej2` at 3e8f8e5. During the review the concurrent session committed a6ed653 (Proposition D‴.32's
code and tests, new files only) and e325236 (the manuscript), neither touching the material. By the end of the review it also had
uncommitted working-tree edits to files of the material: PARADIGM.md (among others around D‴.30(i)–(ii)), `measurement_noise.py`, both
test files, and `eiv_remedy.md` with its raw. Those edits are not reviewed here. Everything here is the committed 3e8f8e5 version,
extracted with `git archive` into `ref_D3t/tree`; the repository's scripts and tests were run from that tree.
`git diff 637a051 9204ea7` also contains round 20 (ac7795c: D‴.26–D‴.29), which is outside the material. The round-21 changes are
ac7795c..9204ea7, plus 3e8f8e5's new Status and proof reference for D‴.30(ii). That reference is Lemma D‴.26a, taken as given.

Line numbers are 3e8f8e5's:
- PARADIGM, D‴.30: l. 3589–3719. (i) l. 3598–3604, (ii) l. 3606–3616, (iii) l. 3618–3628, "Hence" l. 3630–3631; proof l. 3633–3655;
  Status l. 3657–3681; *Numerics* l. 3684–3719 ((A) 3684–3698, (B) 3698–3707, (C) 3707–3715, (D) 3715–3719).
- PARADIGM, D‴.31: l. 3721–3855. (i) l. 3727–3736, (ii) l. 3738–3751, (iii) l. 3753–3757, (iv) l. 3759–3781, "Hence" l. 3783–3786;
  proof l. 3788–3806; Status l. 3808–3828; *Numerics* l. 3831–3855.
- PARADIGM, elsewhere: the §0 rows l. 135–136; the §7 rows l. 6291–6292; the §8 bullets l. 6555–6563 (D‴.31) and 6564–6571 (D‴.30).
- RESULTS: l. 735–761 (D‴.30), 763–787 (D‴.31).
- LITERATURE_TRANSPORT: [S57] l. 77, [S58] l. 78; the audit bullets l. 447–459 (D‴.30) and 460–473 (D‴.31).
- PAPERS: paper 12, l. 170–195.
- `measurement_noise.py`: the docstring l. 1–31; `zero_frequency_operators` l. 108–113; `_whitening` l. 336–348; `second_order_mean`
  l. 351–433; `mean_biases` l. 436–449; `separation_drift` l. 452–478.
- `exp_measurement_noise.py`: `part_A` l. 105–140, `job_B` l. 144–174, `job_B2` l. 177–192, `job_D0` l. 241–258, `job_D` l. 261–283.
- `exp_eiv_remedy.py`: `job_A` l. 65–77, `part_A` l. 80–101, `part_B` l. 105–137.
- `test_measurement_noise.py`: l. 83–102 (the separations), 105–128 (the OU noise and the law), 131–152 (the force).
- `test_eiv_remedy.py`: l. 71–95 (the plan switch), 98–127 (the block, RH + B = I), 156–185 (the desired regressor), 211–231 (the forces
  against the second-order mean).

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
(or `-I`) and one thread (`OMP_NUM_THREADS=1`). pytest ran with `-p no:cacheprovider` on the scratch tree: the two test files pass, 13 tests
in 107.3 s (`ref_D3t/tests.out`). `--from-raw` on the scratch tree rewrites both mds byte-identically: `measurement_noise.md` md5
`188158eb…` and `eiv_remedy.md` md5 `16ee467a…`, each equal to `git show 3e8f8e5:…`.

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3t/scripts/`.
Only `mutants.py` and `test_sep_values.py` run the repository's code. `raw_D.py`, `compare_raw_A.py`, `raw_checks.py`, `scaling.py` and
`fits_B.py` read its records (raw npz or md). All others are independent (no repository import).

| script | what it does |
|---|---|
| `hd.py` (module) | hyper-dual numbers: exact first and second directional derivatives, vectorised over time |
| `inst.py` (module) | independent instance: the saturated network, f, f′, f″, the exact representation at q = 1/3; the field of the weighted K₀-atom loop (MR, DR, DR*; η = ∞ or 0) in physical coordinates with the OU noise as a state; its exact Jacobian and Hessians |
| `check_inst.py` | self-consistency of `inst.py`: the representation, the closed-form gradient, the Jacobian and the Hessians against central differences |
| `lin.py` (module) | RK4 periodic covariance (monodromy, Kronecker discrete Lyapunov); the exact forcing ½∂²F[P], split by covariance blocks; the averaged responses; the exact periodic second-order mean |
| `so_MR.py`, `run_so_all.sh` | the full second-order mean of MR, DR, DR* under both priors at (k_e, B) = (25, 100), (50, 100), (100, 100), (200, 100), (25, 1000), (200, 1000) |
| `compare_raw_A.py` | those 36 cases against the raw records of `eiv_remedy`; the derived shares and ratios |
| `mean_function.py` | the exact, time-resolved second-order bias of the mean function |
| `filtered_gram.py` | the exact periodic mean against the averaged response with the stiffness H − ⟨Ĵ̇ᵀĴ̇⟩/k_e² |
| `law_part.py`, `law_part_ratio.py`, `scaling.py`, `modes.py` | the k_e-scaling of the parameters' part; the per-mode stationary mean deviation and relaxation times |
| `sep.py` | S_k independently: quasi-static at the lift, with the exact periodic moments, with the parameters' cross term, with the stationary mean error |
| `transient_sep.py` | Part D's setting at second order: the transient covariance and mean from the lift, the separation drift along them, the predicted fitted rates |
| `sim_lift.py` | Part D's std-0.05, seed-0 run with an independent simulator (the lift m = 15), with the linearised separation and finite-time exponents along it |
| `sim_batch_sep.py`, `window_threshold.py` | the K₀-atom loop and the linearised separation on 24 antithetic noise pairs per std (new seeds), Part D's estimator and the top finite-time exponent. The run at std 0.01, 0.03, 0.05 is `sim_batch_sep_300.out`. Afterwards a single-path option (argument 4) was added, which leaves the batch path unchanged; it reproduces Part D's seed-0 path (`sim_batch_partD0.out`) |
| `partD0.py` | Part D0 independently: the lift against the K₀-atom loop on one path |
| `partB.py`, `run_partB.sh`, `run_partB2.sh` | Part B: the exact periodic covariance under white noise at k_e = 25–1600, η = 0 and ∞, with the slope and in the control |
| `fits_B.py` | fits of Part B's deviations with and without odd powers; the out-of-sample point k_e = 1600 |
| `check_ii_iii.py` | D‴.31(ii) (the definitions' identity; a dynamic block with a non-commuting anchor gain) and (iii) (random S, a counterexample, the anchor family's S₀) |
| `raw_checks.py`, `raw_D.py` | the raw records without `report()`: ‖P E[e]‖²_ρ, the crossovers, Part B's family, Part C's sampling ratio, Part D's fits |
| `mutants.py`, `test_sep_values.py` | 13 mutants of `measurement_noise.py` (patched in memory) against the round-21 tests; the separation test's measured values |
| `derived.py` | ratios from the outputs above |

## Verdicts

| object | verdict |
|---|---|
| **D‴.30(i) coherence and the threshold** | **Invariance holds, exactly. The threshold formula holds at the lift only (ID-1).** S_k evaluates E[e] as the noise's own mean error at fixed parameters. Along the noisy class trajectory, E[e] also carries the mean tracking error driven by the developing errors-in-variables deviation. That term is of the same order. Over Part D's window it raises the separation drift by 7.5% (9 standard errors at std 0.01). The averaged stationary point has a different drift and a different binding class. Part D's agreement at std 0.05 is a cancellation against the next order. The threshold std over 300 time units from the lift is 0.030, against the stated 0.031. |
| **D‴.30(ii) and Part B** | **Holds.** It is proved for the linearised loop by Lemma D‴.26a (given). Every Part B number reproduces independently at k_e = 25–800, and I added k_e = 1600. The stated order is supported: the increments of k_e²·(ratio − 1) double per doubling, an even fit fails (residual 2.7·10⁻⁵ against 2·10⁻⁷), and the odd fit predicts k_e = 1600 out of sample. The control's evenness is numerical, not proved (ID-6). |
| **D‴.30(iii) with the full second-order mean** | **Holds; the computation is right.** An independent exact computation (hyper-dual Hessians, my own periodic covariance) reproduces `second_order_mean`: the noise parts to ≤ 1.2·10⁻⁸ and the full biases to ≤ 5.3·10⁻⁷ in 36 cases. The exact periodic means 504/2653 and all of B2 reproduce. The excess over the averaged system is the filtered Gram's (2695.5 against 2653.0 at k_e = 25; ≤ 0.17% from k_e = 50). The weighting k_e/c(t) is the slaving's, correct at leading order. The (δe, ν) / rest split is consistent with the text. The covariance blocks are symmetric. |
| D‴.30 proof, Status | The proof is correct in structure, and the round-21 additions (the hypothesis ∂²_xf̂ = ∂²_xf, the parameters' terms, the notation) are right. The Status is accurate except "the attraction threshold formal … and checked on the real loop" (ID-1) and the control's evenness stated as fact (ID-6). |
| **D‴.31(i)** | **Holds.** Re-derived: u := T̂ − z satisfies du = dT̂ − ηu ds for T̂ of bounded variation, with the integrated form for measurable T̂. The split-row term −σ₀(I − M)z has M row-stochastic. |
| **D‴.31(ii)** | **Holds as RH + B = I.** It is an identity of the definitions R = (H + S₀)⁻¹, B = RS₀ (3.6·10⁻¹² for random unstructured S₀). Its content is that a dynamic block has these R and B: I confirmed this for a non-commuting anchor gain to 7.4·10⁻¹⁵. The per-mode identity holds in H's eigenbasis, whatever S₀. The numerics' and the test's "check" for non-commuting leaks is vacuous (ID-5). |
| **D‴.31(iii)** | **Holds.** Re-proved by the congruence. S + Sᵀ ⪰ 0 suffices; the exact condition is S + Sᵀ + SᵀH⁻¹S ⪰ 0, and S = −(h₁/2)v₁v₁ᵀ gives the ratio 4. New: every anchor gain K with K + Kᵀ ⪰ 0 gives S₀ + S₀ᵀ ⪰ 0 (2MMᵀ − η(M + Mᵀ) = η(K + Kᵀ) + 2KKᵀ, M = η + K). So (iii) covers the whole matrix-anchor family. |
| **D‴.31(iv)(c) restated** | **Holds as restated.** The noise part is least-squares with its bound. The parameters' part is named correctly and quantified correctly: all 36 Part A cases reproduce, including +3.25% and +56.19% for DR, −8.69% to −0.15% for DR*, MR ≤ 0.030%, and MR over DR 2.12·10⁸ at k_e = B = 100. "Small on the instance and falling with k_e" (PAPERS, §8, LITERATURE) is overstated (ID-2). The "new" parameters' errors in variables is classical (ID-3). |
| D‴.31 proof, Status | correct; the Status is accurate on what is proved; the novelty (ID-3) |
| *Numerics* (both) | **Every number reproduces** (confirmations 3–9). Readings that do not hold: (D) the measured rates read against the predictions −2.53·10⁻⁴, −1.5·10⁻⁵ ("the threshold") and +4.62·10⁻⁴ (ID-1); (B) "RH + B = I holds to 1.4·10⁻¹⁴ … for zero-frequency leaks that do not commute" (ID-5); "without the slope the expansion is even" (ID-6). |
| `measurement_noise.py` | **Holds.** `second_order_mean`, `mean_biases` (called consistently: `out` from the matching prior's loop, `S` = σ₀ for η = 0), `_whitening` and `zero_frequency_operators` reproduce. `separation_drift` reproduces (0.014809, 0.010085, 0.298149) but computes the lift's drift only (ID-1). The docstring is stale (ID-6). |
| experiments, mds, raw | **Hold.** Both mds regenerate byte-identically. Part D's std-0.05 seed-0 run reproduces to all printed digits with two independent codes, and D0 reproduces (4.4·10⁻¹⁶, spread 0.0). |
| tests (13) | **Pass. They catch the swapped noise and parameters' parts, a wrong c(t), a doubled force and the wrong sign of S_k.** They do not catch S_k ×0.85, a dropped E[e] term in S_k, or `mean_biases` ignoring S. The exact periodic mean, the fixed-prior branch and the mean-function bias are never executed (ID-4). |
| §0, §7 rows, §8 bullets, RESULTS, PAPERS 12, LITERATURE | **Numbers hold** and match the mds. Wording: the threshold and "rates at their predictions" (ID-1); "small on the instance and falling with k_e", "+3.3% of DR's bias" (ID-2); "new … the parameters' own errors in variables" (ID-3); "RH + B = I to 1.4·10⁻¹⁴" (ID-5); "the slope-free control even" (ID-6). |

### Fix audit (round 21)

| round-21 issue | fix in text and code | verdict |
|---|---|---|
| D3q-1 (the order) | Status l. 3660–3663, *Numerics* (B) l. 3698–3703, §0, RESULTS l. 747–749; Part B to k_e = 800 with the control | **applied, correct** (confirmation 5); the control's evenness worded as proved (ID-6) |
| D3q-2 (averaged values as the instance's) | (A) l. 3689–3695, (B) l. 3704–3707, (C) l. 3712–3714; §0, §7, RESULTS, PAPERS l. 181–182; md B2 and the mode table | **applied, correct**: ×28.0, 504/2653 (×5.27) and the excesses 10.7%–0.17% reproduce independently |
| D3q-3 (invariance is not attraction) | (i) l. 3601–3604; `separation_drift`; Part D; the test | **applied, but the new threshold omits a same-order term (ID-1)**; the module docstring l. 13–14, flagged in D3q-3, is unchanged (ID-6) |
| D3q-4 (the omitted terms, ∂²f̂ = ∂²f, "grows like k_e²") | l. 3591, 3621–3623, 3627, proof l. 3648–3651; md exponents | **applied, correct** (local exponents 1.472, 1.169, 0.830 reproduce) |
| D3q-5 (the force untested) | the force test l. 131–152, OU test l. 105–128, test 2 at η = 0, σ₀, ∞ | **applied**; the force mutants are now caught; other blind spots remain (ID-4) |
| D3q-6, D3r-3 (novelty) | [S57], [S58]; Status l. 3666–3675, 3814–3822; audit l. 447–473; PAPERS l. 185–190 | **applied**; bibliographic details of [S58] and Mayyas–Aboulnasr re-checked; one new over-claim (ID-3) |
| D3q-7 (wording) | l. 3687 (exact differential, 104…5.5), 3689 (10%), 3712 (6.0%, 6.4%), 3708 (common paths), §8 l. 6570–6571, proof l. 3638–3641, n = 1, fast part l. 3615–3616 | **applied** |
| D3r-1 (the parameters' fluctuation) | `second_order_mean`, `mean_biases`; md Part A; (iv)(c) l. 3772–3778; Status l. 3812–3814; rows, RESULTS, PAPERS; test l. 211–231 | **applied, correct** (all 36 cases independently, confirmation 3); residual wording (ID-2) |
| D3r-2 (RH + B = I) | (ii) l. 3738–3751, proof l. 3792–3795, `zero_frequency_operators`, md B, test l. 122–127 | **applied, correct**; the non-commuting "check" is the definitions' identity (ID-5) |
| D3r-4 (rescaled projected error) | (iv)(c) l. 3778–3781; md columns; RESULTS l. 772; §7 | **applied, correct**: ‖P E[e]‖²_ρ ∝ (B/(k_e + B))² to 0.1% (5.984–5.979) |
| D3r-5, D3r-6 (BV, split rows; S + Sᵀ ⪰ 0) | (i) l. 3727–3736, proof l. 3788–3791; (iii) l. 3753, proof l. 3795–3796; docstrings l. 311, 316–317 | **applied, correct** |
| D3r-7 (numerics wording) | l. 3838–3855; md; `mean_error` docstring l. 270–276; DR*'s jitter l. 3785; PAPERS | **applied**: 46.6 against the law's 46.7; sampling ratio 0.43–1.83 recomputed |
| D3r-8 (tests) | plan-switch test l. 71–95; active face l. 122–127; regressor's point l. 176–185; the forces l. 211–231 | **applied**; the active-face assertion is the definitions' identity (ID-5) |

No number in the rows, RESULTS or PAPERS disagrees with the regenerated mds (checked item by item).

## Confirmed issues

**ID-1 (major; D‴.30(i) l. 3601–3604 "a separation of class k decays at σ₀/k_e − Eν²λ_max(S_k) to leading order, S_k the period average
of the second-order mean of the law's Jacobian … along the noisy class trajectory (… + E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩), so the replicated
configuration attracts — and an m-system started near it follows the K₀-atom loop — only below the noise level Eν² =
(σ₀/k_e)/max_kλ_max(S_k)". Also: Status l. 3657–3658 "the attraction threshold formal … and checked on the real loop"; *Numerics* (A)
l. 3696–3698 and (D) l. 3715–3719; §0 row l. 135 "the separations' threshold std 0.031, the real loop's separation rates at their
predictions"; §7 row l. 6291; RESULTS l. 750–752; PAPERS l. 174 and 188–189 and LITERATURE l. 458 ("its noise threshold" as new);
`separation_drift` l. 452–478; md l. 47 and Part D.)** The separation drift omits a term of its own order: the mean tracking error that
the parameters' errors-in-variables deviation drives. The stated rate is the rate at the lift. It is not the rate along the noisy
trajectory, and the long-time threshold does not follow from it.

*The omission* (`sep.py`). A separation obeys u̇ = (G^{1/2}ε_m∇²_ϑΦ̃(y; ϑ_k)G^{1/2} − σ₀/k_e)u. Its second-order mean drift is
- E[ε_m]M + E[ε_m²]∂_xM, with M := G^{1/2}∇²_ϑΦ̃G^{1/2}, plus the parameters' cross term Σ_lE[ε_mδϑ_kl]∂_{ϑ_l}M;
- with E[ε_m] = E[e] the mean tracking error along the trajectory.

D‴.31(iv)(c)'s own mean balance gives k_eE[e] = g − Ĵȳ, where ȳ is the mean parameter deviation. `separation_drift` (l. 473) uses
`mean_error`: the noise's own E[e] at fixed parameters (ȳ = 0), whose rms is 0.0274 per unit Eν². The Ĵȳ/k_e term is of the same order
in Eν². It grows as ȳ develops, toward D‴.30(iii)'s errors-in-variables response.

The parts of λ_max(S₃) = 0.298149 at the lift, all reproduced independently (N = 512 and 4096 agree):
- the regressor's part: 0.3048 along the top eigenvector;
- the quasi-static E[e] part: −0.0066;
- with the exact periodic moments and the cross term: 0.298305. The cross term is negligible.

*Over Part D's window* (`transient_sep.py`: the exact second-order transient from the lift, ν stationary, P and m otherwise 0).
- The mean error's rms per unit Eν² is 0.121, 0.275, 0.394, 0.484 and 0.553 in successive 50-unit windows to t = 251.
- The period-averaged drift's λ_max is 0.306, 0.319, 0.329, 0.337 and 0.343 over the same windows. Modes 9 and 7 equilibrate in the window
  (k_e/h_j = 128 and 282) and carry 6.7% and 10.5% of the bias (`modes.py`).
- Part D's estimator is the least-squares slope of log‖u‖ over t ≥ 50, from S₃'s top eigenvector. Under this drift it is predicted at the
  noise-induced rate λ_eff := (rate + σ₀/k_e)/Eν² = 0.326–0.328 at every std (−2.504·10⁻⁴, +1.10·10⁻⁵, +5.36·10⁻⁴ at 0.01, 0.03, 0.05).
  Without any mean error λ_eff is 0.3045. The proposition's value is 0.298.

*Measured* (`sim_batch_sep.py`, `window_threshold.py`). The K₀-atom loop and the linearised separation ran on 24 antithetic noise pairs per
level, with seeds different from Part D's. The pairing cuts the first-order (odd) part of the log-growth. The same code on Part D's
std-0.05 seed-0 path gives +4.369022·10⁻⁴, the stored value.

| std | rate (pair means) | λ_eff | against λ_max(S₃) = 0.2981 | second-order transient prediction |
|---|---|---|---|---|
| 0.01 | −2.5094·10⁻⁴ ± 0.0025·10⁻⁴ | 0.3204 ± 0.0025 | +7.45% (9.0 s.e.) | 0.3261 |
| 0.03 | −0.012·10⁻⁵ ± 0.169·10⁻⁵ | 0.3143 ± 0.0019 | +5.41% (8.6 s.e.) | 0.3267 |
| 0.05 | +4.4532·10⁻⁴ ± 0.040·10⁻⁴ | 0.2913 ± 0.0016 | −2.29% (4.3 s.e.) | 0.3278 |

- A weighted fit λ_eff = 0.3239 − 12.85Eν² has the small-noise limit 0.3239. The transient prediction matches it to 0.7%. The fall with Eν²
  is the next order.
- Part D's agreement at std 0.05 (+4.37·10⁻⁴ against +4.62·10⁻⁴) is therefore the net of two effects: a second-order excess S_k omits
  (0.326/0.298 = +9.4%) and a higher-order deficit. At std 0.01, Part D's four paths (s.e. 0.037·10⁻⁴) cannot resolve the 7.5% excess.
- Over this window the rate changes sign at std 0.0301 (Eν² = 9.06·10⁻⁴), against the lift's 0.0308. The number survives here.
- *Lyapunov exponent or mean.* S_k is a mean drift; Part D measures a finite-time exponent. At this order they differ by the cumulant of
  the first-order multiplicative noise. That cumulant is small because the measured error's noise is high-passed (low-frequency density
  a²N₀/k_e²). The measured limit equals the mean-drift prediction to 0.7%. So the threshold is the sign change of the mean drift, but of
  the time-dependent drift.

*At the averaged stationary point* (`sep.py`: E[e] the exact periodic second-order mean's, rms 2.007 per unit Eν², 73 times the noise
part).
- Adaptive prior: λ_max = 0.0753, 0.0530 and 0.0696 in classes 1–3. Class 1 binds, at the extrapolated std 0.061.
- Fixed prior: 0.0064, 0.0018 and 0.2836 (std 0.0316).
- The expansion does not apply at these noise levels. Mode 2's stationary deviation is 957Eν², i.e. 0.92 in whitened coordinates at std
  0.031 (`modes.py`). So the expansion decides neither the long-time threshold nor which class binds. "Attracts … only below" is a claim
  about the noisy system for all times. It rests on the lift's drift.

*Fix:*
- (i): "S_k at the lift (the parameters at their atoms; E[e] the noise's own mean error). Along the noisy trajectory E[e] also carries the
  mean tracking error −Ĵȳ(t)/k_e of the errors-in-variables deviation, of the same order. So the rate is σ₀/k_e − Eν²λ_max(S_k(t)). On the
  instance λ_max(S₃(t)) rises from 0.306 to 0.343 over 250 time units, and the measured small-noise rate is 0.320. The threshold
  (σ₀/k_e)/max λ_max(S_k) is the lift's (std 0.031; 0.030 over 300 time units from the lift). The long-time threshold is open: at the
  stationary point the expansion's drift is 0.070 (class 3) and 0.075 (class 1, adaptive prior), where the parameters' mean deviation is
  O(1)."
- Status, the rows, RESULTS: "formal at the lift". Replace "at their predictions" with the λ_eff figures. Novelty: "a noise threshold at
  the lift".
- Code: let `separation_drift` take E[e](t) and the cross covariance from `second_order_mean`'s transient or periodic solution.
- Experiment: run Part D on antithetic pairs and print λ_eff with its standard error.

**ID-2 (minor; PAPERS l. 179 "the parameters' fluctuation adds a force that is not least-squares, small on the instance and falling with
k_e"; §8 l. 6561–6562 "the parameters' fluctuation adding a force that is not, small on the instance and falling with k_e"; LITERATURE
l. 471–472 "(small here, +3.3% of DR's bias at k_e = 25)").** The wording holds only for the part's absolute size. The instance's fall with
k_e has a one-dimensional cause that the text does not state (`compare_raw_A.py`, `scaling.py`, `law_part.py`).
- *Share of the bias.* Relative to DR's bias the part is not small: +3.25%, +3.05%, **+56.19%**, −0.75% at k_e = 25, 50, 100, 200
  (B = 100). DR*'s is −8.69% at k_e = 25. "+3.3%" is the most favourable case quoted alone.
- *Falling with k_e.* The part's own response falls like k_e⁻²: 3.92·10⁻⁴, 1.13·10⁻⁴, 2.93·10⁻⁵, 7.40·10⁻⁶ at B = 100, i.e. k_e² times
  it is 0.245, 0.282, 0.293, 0.296. Its share of DR's bias does not fall monotonically.
- *The mechanism.* The law's term k_e⟨∇²_ϑΦ̃E[ε_mδϑ]⟩ falls like 1/k_e (k_e times its norm is 8.4·10⁻³–8.9·10⁻³ for k_e = 25–200).
  - Pointwise, E[ε_mδϑ] = O(N₀/k_e). Its rms times k_e is 6.9·10⁻³, 5.8·10⁻³, 4.4·10⁻³, 2.9·10⁻³, proportional to B/(k_e + B)
    (divided by it: 8.65·10⁻³–8.71·10⁻³, `law_part_ratio.py`).
  - So the slow-time term's leading O(N₀) part is ∝ ⟨M(x_m)∂_xĴ(x_m)ẋ_m⟩ = T₀⁻¹∮M∂_xĴ dx, which vanishes on a closed orbit in ℝ.
  - That cancellation needs n = 1. On a curve in ℝⁿ (n > 1) such a line integral does not vanish in general (an analytic remark, not
    computed). That matters for paper 12's planned second-order plants.

*Fix:*
- PAPERS and §8: "small in absolute terms on the instance (≤ 3.9·10⁻⁴ per (Eν²)²); −8.7% to +3.3% of DR's and DR*'s biases, and +56% of
  DR's at k_e = B = 100; its own response falls like k_e⁻² on a one-dimensional reference, where the law's leading term averages out".
- LITERATURE: add "+56% at k_e = B = 100".

**ID-3 (minor; D‴.31 Status l. 3820–3822 "New are their transfer … and the instance's quantification, including the parameters' own
errors in variables, which the desired regressor does not remove"; the §0 row l. 136, last column; PAPERS l. 188–190; LITERATURE
l. 470–472).** The mechanism listed as new is the classical second-order bias of a constant-gain stochastic approximation from its own
parameter fluctuations: the mean field's curvature in the parameters, contracted with the stationary covariance, plus the noise–parameter
correlation.
- Benveniste–Métivier–Priouret (already [S57]) treat such second-order expansions.
- For constant step-size SGD on non-quadratic objectives, Dieuleveut–Durmus–Bach show that the stationary mean is O(γ) off the optimum and
  remove the bias by Richardson–Romberg extrapolation (*Bridging the gap between constant step size stochastic gradient descent and Markov
  chains*, Ann. Statist. 48(3):1348–1382, 2020).
- The model's curvature in the parameters biasing a least-squares fit at O(σ²) is Box's (*Bias in nonlinear estimation*, JRSS B
  33:171–201, 1971). It is the analogue of the term ½D²_Θf̂[Σ_ΘΘ] and its least-squares force.
- Bibliographic details checked by a web search; the content is from the abstracts and memory.

*Fix:* cite them. Restrict "new" to the instance's quantification and the observation that the desired regressor does not remove the term.

**ID-4 (minor; `tests/test_measurement_noise.py` l. 83–102 and 131–152; `tests/test_eiv_remedy.py` l. 98–127 and 211–231; the Status
"Checked" lists l. 3675–3681 and 3822–3828; RESULTS l. 757–761 and 781–787).** The new tests catch the gross errors but leave the
round-21 numbers' main code paths unexercised.

Mutants (`mutants.py`, patched in memory):

| mutant of `measurement_noise.py` | result |
|---|---|
| `second_order_mean`: noise and parameters' parts swapped | caught (l. 149, l. 224) |
| `second_order_mean`: c(t) = k_e; c(t) = k_e + a | caught by the eiv test (l. 224, through DR*); the MN test passes (c = k_e for MR) |
| `second_order_mean`: the ½ dropped; the noise block ν only | caught |
| `eiv_force`: no E[e] (MR); E[ε_m²] = Eν² | caught |
| `zero_frequency_operators`: B = S₀R | caught (l. 125) |
| `separation_drift`: sign flipped; ×1.15 | caught (l. 91; l. 101) |
| `separation_drift`: **×0.85**; **the E[e] term dropped** | **pass** |
| `mean_biases`: **S ignored** | **passes** (no test calls it with S) |

- *Never executed by a test* (grep): `second_order_mean(periodic=True)`, i.e. the exact periodic mean behind 504, 2653 and all of B2;
  `second_order_mean` at η = 0 (every "fixed" entry); `mean_biases` with `S`, and its `function`/`c_bar`; Part B's control floor.
- *The separation test's tolerance.* It measures −2.5150·10⁻⁴ and +4.5013·10⁻⁴ against −2.5316·10⁻⁴ and +4.6240·10⁻⁴
  (`test_sep_values.py`). The differences, 1.7·10⁻⁶ and 1.2·10⁻⁵, sit 60 and 8 times inside the tolerance 10⁻⁴, which is 35% of σ₀/k_e.
  It cannot see ID-1's 7.5%.

*Fix:*
- Test `mean_biases(…, S = σ₀)` against (H + σ₀I)⁻¹, and the η = 0 branch.
- Test the exact periodic mean against the filtered-Gram averaged response (2%, k_e = 25) or against a transient.
- Test `function` against ‖Ĵȳ + c̄‖² computed independently.
- Test `separation_drift` against an independent quadrature (third derivatives), and the rate with antithetic pairs at a tolerance near
  10⁻⁵.

**ID-5 (cosmetic; *Numerics* (B) l. 3847–3848 "RH + B = I holds to 1.4·10⁻¹⁴ and per mode to 8.5·10⁻¹⁵" for the non-commuting leaks;
the §0 row l. 136 "RH + B = I to 1.4·10⁻¹⁴"; §7 l. 6292; RESULTS l. 776; md l. 56–63; `exp_eiv_remedy.py` l. 129–135; test
l. 122–127).** For these rows R and B come from `zero_frequency_operators`, R = (H + S₀)⁻¹ and B = RS₀. RH + B = R(H + S₀) = I is then an
identity of the definitions (`check_ii_iii.py`: 3.6·10⁻¹² for 200 random, unstructured, even destabilising S₀). The per-mode identity is
its diagonal in H's eigenbasis. The printed residual measures `numpy.linalg.inv`.
- The content of (ii) is that a prior's dynamic block has these R and B. That is checkable. For y′ = −Hy − σ₀(y − z) + F,
  z′ = −η(z − y) − K(z − δ) with K = κP_amp (κ = 10η and 10⁶η, non-commuting), the block's R and B equal (H + S₀)⁻¹ and (H + S₀)⁻¹S₀,
  S₀ = σ₀(η + K)⁻¹K, to 7.4·10⁻¹⁵ and 5.3·10⁻¹⁵. There max|HR + B − I| = 0.736 and 0.785.
- The family rows (`stationary_operators`, from the block) are a genuine check.

*Fix:* derive R and B for the non-commuting rows and in the test from a block, as the family rows do, or label the residual "an identity of
the definitions".

**ID-6 (cosmetic; wording).**
- *The control's evenness stated as fact.* Status l. 3662–3663 "(with the noise input −k_eν the deviation is even in 1/k_e)", *Numerics*
  l. 3703 "without the slope the expansion is even", RESULTS l. 749, and the §0 row "the slope-free control even in 1/k_e". Lemma D‴.26a's
  Status (l. 3272) says "(v) does not prove it". The evenness is numerical (`fits_B.py`): the successive ratios reach 4.002 at
  k_e = 800 → 1600, and an even fit to k_e ≤ 800 predicts d(1600) to 0.06%. *Fix:* "numerically even (k_e ≤ 800)".
- *"adaptive (every η > 0)"* labels the full-force and mean-function columns (measurement_noise.md l. 37 via `exp_measurement_noise.py`
  l. 327; eiv_remedy.md l. 7). These are computed at η = ∞. The parameters' part depends on η through the covariance: the weakest mode's
  variance is 1.0658 at η = σ₀ against 2.0875 at ∞ (md Part A). *Fix:* "(η = ∞; the noise part for every η > 0)".
- *The module docstring* l. 13–14 "so the separations are not excited" holds at first order only. It was flagged in D3q-3 and is
  unchanged. *Fix:* "not excited at first order; their stability needs the noise below a level (`separation_drift`)".
- *(iv)(c) l. 3773* "and its slope–parameter term": DR* has none, since its compensation is at x_m.
- *(iv)(c) l. 3776* "which is not a least-squares force": every F is ⟨Ĵᵀ(ĴH⁻¹F)⟩. The intended content is that no bound in the noise's
  moments applies. *Fix:* "not the learning force of a mean tracking error (no bound in the noise's moments)".

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** Both `--from-raw` runs are byte-identical, and the two test files pass (13 tests, 107.3 s).
2. **The independent instance** (`check_inst.py`):
   - max|f̂ − f| = 4.4·10⁻¹⁶ on [−4, 4];
   - the closed-form gradient equals the hyper-dual one to 8.9·10⁻¹⁶;
   - the exact Jacobian equals central differences to 2.1·10⁻⁹, and the Hessians equal central differences of the exact Jacobian to
     ≤ 2.5·10⁻⁹ (MR, DR, DR*).
3. **The full second-order mean, independently** (`lin.py`, `so_MR.py`, `compare_raw_A.py`).
   - *Method.* Physical coordinates (x, Θ, ν) with the OU noise as a state. Exact hyper-dual Hessians of the true field, not a stencil. RK4
     with N = 4096 and 8192 identical to the printed digits; the covariance periodic to 1.7·10⁻¹³.
   - *Against the raw records,* 36 cases (MR, DR, DR* × both priors × six (k_e, B)). The noise parts agree to ≤ 1.2·10⁻⁸ relative, the full
     biases to ≤ 5.3·10⁻⁷ (the worst DR at k_e = B = 100, where the noise part nearly cancels), the mean function to ≤ 4.9·10⁻⁷, and the
     parameters' part alone to ≤ 4·10⁻⁵ (at 7.6·10⁻⁸). So the 4th-order stencil at h = 2·10⁻³ is accurate.
   - *Code reading.* The covariance blocks are symmetric. The split assigns the cross blocks (δe, ν)–ϑ to the parameters' part, as the text
     says; a mutant moving δe into the parameters' part is caught.
   - *The exact periodic means of MR:* 2652.957, 503.720 (k_e = 25); 6824.751, 1374.643; 15048.456, 3080.562; 26612.417, 5472.311. They
     equal the md. Adaptive over fixed is 5.267, and the excess over the averaged system 10.65% to 0.16%.
   - *New, DR and DR*:* their exact periodic means lie +0.23% and −0.94% from the averaged full response at k_e = 25, and −0.29% and −0.03%
     at k_e = 100.
   - *The weighting.* k_e/c(t) is the slaving m_e = (f_e − Ĵm_y)/c. Making DR*'s stiffness consistent with it (⟨k_eĴᵀĴ/(k_e − a)⟩ for H)
     moves DR*'s bias by −0.40% (k_e = 25).
   - *The mean function* (`mean_function.py`). Time-resolved, with the exact m(t) and P(t), it is 2652.93, 0.47699, 0.05064 against the
     averaged 2397.65, 0.47591, 0.05111. Using the period average Σ̄ is immaterial: ⟨‖c_t‖²⟩ = 8.0990·10⁻⁵ against ‖c̄‖² = 8.0988·10⁻⁵.
4. **The filtered Gram** (`filtered_gram.py`). The averaged response with H − ⟨Ĵ̇ᵀĴ̇⟩/k_e² gives 2695.5 at k_e = 25 against the exact
   2653.0 (0.984); the exact over it is 0.9983, 0.9997, 0.99994 at k_e = 50, 100, 200. Under the fixed prior: 0.9977–0.99998. *Numerics* (B)'s
   explanation holds.
5. **Part B, independently** (`partB.py`, `fits_B.py`).
   - *The table.* All 24 table deviations reproduce: −1.969058·10⁻² … +1.309278·10⁻⁵ (η = ∞) and −2.375209·10⁻² … +1.161518·10⁻⁵
     (η = 0); the control −2.190887·10⁻² … −2.411236·10⁻⁵ and −2.625627·10⁻² … −2.860935·10⁻⁵. 2¹⁸ steps change k_e = 800 by 7·10⁻¹⁰.
   - *Other quantities.* tr Q = 8.783561 = 8.104920 + 0.678641 (the cross term 4·10⁻¹⁷); tr H² = 0.065712. The weakest mode's
     adaptive/fixed ratio is 27.989, 25.108, 24.373, 24.189, 24.143, 24.131, and exact/(slow + fast) is 0.8446 (fixed) and 0.9798 (adaptive)
     at k_e = 25.
   - *The order.* k_e²·d changes by −0.42, +1.10, +2.79, +5.73, +11.49 per doubling (η = ∞); the increments double, the signature of
     c₁k_e. The odd fit gives c₁ = +0.0290 (η = ∞) and +0.0311 (η = 0) with residuals ≤ 2·10⁻⁷; the even fit's residual is 2.7·10⁻⁵.
     Out of sample, the odd fit predicts k_e²d(1600) = 31.62 against the computed 31.38; the even fit predicts −7.02.
   - *The control:* increments ×¼ per doubling, c₁ = 2·10⁻⁴, d(1600) = −6.0253·10⁻⁶ against the even fit's −6.0289·10⁻⁶.
6. **Part D and D0** (`sim_lift.py`, `sim_batch_sep.py`, `partD0.py`, `raw_D.py`).
   - *Std 0.05, seed 0* (Part D's noise path): the independent lift simulator gives the fitted rate +4.369022·10⁻⁴ and the separation norms
     3.250073·10⁻⁶ (t = 50) and 3.630480·10⁻⁶ (t = 300), all equal to the stored record. The linearised separation from the same vector
     gives the same slope; the top finite-time exponent over [0, 300] is +4.336·10⁻⁴.
   - *D0:* the lift and the K₀-atom loop on one path differ by 4.4·10⁻¹⁶ in f̂, and the within-class spread is 0.0.
7. **S_k at the lift** (`sep.py`): λ_max = 0.014809, 0.010085, 0.298149; threshold stds 0.1382, 0.1675, 0.0308. The exact periodic
   E[ε_m²] lies within [0.99687, 1.00322] of the quasi-static formula.
8. **D‴.31(ii) and (iii)** (`check_ii_iii.py`): as in ID-5 and the verdicts.
   - *(iii):* over 2000 random S per class the worst-F ratio is ≤ 1.000000 (symmetric S ⪰ 0, skew S) and 0.990 (non-symmetric with
     S + Sᵀ ⪰ 0); 4.0000 for S = −(h₁/2)v₁v₁ᵀ.
   - *The anchor family:* over 20 000 non-normal K with K + Kᵀ ⪰ 0, λ_min(S₀ + S₀ᵀ) ≥ 3.4·10⁻⁵, as the identity predicts.
   - *Part B's family against the closed forms* (`raw_checks.py`): S₀ = σ₀κ/(η + κ), r₁ = 1/(h₁ + S₀), b₁ = S₀/(h₁ + S₀) to the printed
     digits (3269.38 … 135.497; 0 … 0.958556); n_weak = 3.
9. **The raw records without `report()`** (`raw_checks.py`, `compare_raw_A.py`).
   - ‖P E[e]‖²_ρ = 3.8296, 2.6583, 1.4949, 0.6643 (B = 100), i.e. 5.98 (B/(k_e + B))².
   - The crossover stds are 0.0012, 0.0850, 0.2704.
   - Part C's relative errors are 0.0052–0.0527, at 0.43–1.83 sampling errors; MR/DR at t = 300 is 8.685 (predicted 9.157).
   - MR over DR (full) is 3009 to 2.123·10⁸. The largest DR/DR* full bias is 1.18453, against ¼‖∂²_xf‖²_ρ = 1.30646. DR's noise part over
     ‖g‖²_ρ is 0.99993.
10. **Literature** (web search, bibliographic details only): Yucelen–Haddad, IEEE TAC 58(4):1080–1085, 2013; Mayyas–Aboulnasr, IEEE TSP
    45(4):927–934, 1997; Dieuleveut–Durmus–Bach, Ann. Statist. 48(3):1348–1382, 2020; Box, JRSS B 33:171–201, 1971.

## Not checked

- The Monte Carlo parts were not rerun: Part C's 48 and 36 runs, and 11 of Part D's 12 stored runs. One Part D run and D0 were reproduced
  independently. ID-1's measurement uses 144 new paths.
- The full second-order mean at finite η (the prior in the state, e.g. η = σ₀) was not recomputed; my code covers η = ∞ and 0. Of the
  eight (k_e, B) grid points, (50, 1000) and (100, 1000) were not run.
- The long-time (stationary) separation behaviour beyond the second-order expansion was not simulated (~10⁵ time units). ID-1's stationary
  figures are the expansion's, which does not apply there.
- The first-order noise's cumulant in the separation exponent was estimated (O(N₀/k_e²)), not computed. It is bounded by the 0.7%
  agreement of the measured limit.
- n > 1 (ID-2's remark on the exact differential) was analysed, not computed.
- Lemma D‴.26a, `noise_decoupling` and `decoupled_covariance` were taken as given. Only the two test files of the material were run.
- The literature's content is from abstracts and memory; no full text was fetched.
