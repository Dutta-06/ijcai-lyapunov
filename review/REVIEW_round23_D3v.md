# Referee report D3v (round 23) — Corollary D‴.33 (the desired regressor keeps the replicated configuration attracting under measurement noise): PARADIGM §5d (statement (i)–(iv), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.33 |", the §8 bullet "Done (Corollary D‴.33)", RESULTS "Corollary D‴.33", PAPERS paper 12, the LITERATURE audit bullet, README, `src/core/measurement_noise.py` (`mean_error`, `eiv_force(variant, return_g)`, `separation_drift(variant)`, `second_order_transient`, `separation_drift_series`), `exp_eiv_remedy.py` Part D (`part_D_pred`, `job_D`, the report's section D, `--only-D`) with `results/chaos/eiv_remedy.md` and its raw, and the D‴.33 test

Branch `claude/festive-bardeen-eopej2` at bfdbc0f. I reviewed the committed version only. It was extracted with `git archive bfdbc0f` into
`ref_D3v/tree`, and the repository's scripts and tests were run from that tree. No repository file was modified. The only git commands run
were log, status, archive, diff and show.

During the review the concurrent session committed a9a0e5a (manuscript part 4e). It touches only `manuscript/` and `lyapunov.pdf`, not the
material. It also had uncommitted edits to `experiments/exp_cohesion_leak.py` and the cohesion leak's md and raw. None of this is reviewed
here.
Python ran with `PYTHONDONTWRITEBYTECODE=1` and one thread (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`), pytest with `-p no:cacheprovider`.
Context taken as given: Lemma D‴.26a, Proposition D‴.30(i) as revised in round 22, D‴.31(iv), and the averaging behind them.

Line numbers are bfdbc0f's:
- PARADIGM, D‴.33: l. 3932–3996.
  - Setting l. 3932–3935; (i) l. 3937–3939; (ii) l. 3941–3944; (iii) l. 3946–3948; (iv) l. 3950–3957.
  - Proof l. 3959–3965; Status l. 3967–3975; *Numerics* l. 3977–3996.
- PARADIGM, elsewhere: the §0 code summary l. 37–38; the §0 row l. 138; the §7 row l. 6434; the §8 bullet l. 6698–6705.
- RESULTS l. 815–833; PAPERS paper 12 l. 171–209 (D‴.33 at l. 184–185, 192–193, 202–203); LITERATURE l. 484–490; README l. 65–67.
- `measurement_noise.py`: the docstring l. 25–31; `mean_error` l. 314–326; `eiv_force` l. 328–351; `separation_drift` l. 497–527;
  `second_order_transient` l. 530–574; `separation_drift_series` l. 577–599.
- `exp_eiv_remedy.py`: the docstring l. 28–31; `part_D_pred` l. 209–249 (the stationary point l. 233–239, the transient l. 240–247);
  `job_D` l. 252–276; the report's section D l. 377–426; `--only-Dpred` l. 442, `--only-D` l. 451.
- `eiv_remedy.md` section D: l. 102–146 (the stationary table l. 121–126, the windows l. 128–137, the real loop l. 139–146).
- `tests/test_eiv_remedy.py`: l. 243–268.

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3v/scripts/`.
Only `mutants.py` runs the repository's code, as do the two runs listed last in the table. `raw_D.py` reads the raw npz. Everything else is
independent: no repository import, and the model is re-implemented from its definitions.

| script | what it does |
|---|---|
| `inst.py` (module) | independent instance: the saturated Gaussian kernel, the plant, the exact representation of q = 1/3, x_m = 2 sin t, G, σ₀; closed-form gradient and Hessian in ϑ; a truncated-Taylor jet arithmetic in (x, α, γ, κ) to order 3, for exact Hessians and third derivatives |
| `check_inst.py` | f̂ = f to 2.2·10⁻¹⁶ on [−4, 4]; closed forms = jets to 3.6·10⁻¹⁵; the jets' third derivatives = central differences of the closed forms to 2.2·10⁻⁸; numpy's Generator gives the same draws in one block as one at a time |
| `lift_drift.py` | S_k at the lift for MR, DR, DR* with the quasi-static moments re-derived; the thresholds; MR's top eigenvector u; the k_e scans at B = 4k_e, 100 and 1000 to k_e = 800; the exact periodic noise moments at fixed parameters |
| `raw_D.py` | Part D's stored predictions and its 32 runs, read from `eiv_remedy_raw.npz` |
| `sim.py` (module) | independent real-loop simulator, vectorised over paths. Two modes: the m = 15 lift with Part D's start, and the K₀-atom loop with the exact variational equation of every class's separations. RK4 with the noise held over a step, Part D's noise paths |
| `repro_partD.py`, `check_std01.py` | Part D's DR runs (std 0.05 and 0.1, seeds 0–3, both signs) and MR runs (std 0.05, seeds 0–1) |
| `so_sep.py` (module), `run_so_sep.py`, `windows_parts.py` | the exact second-order theory of a separation's finite-time growth. It uses the covariance of the linearised loop, augmented with the separation's first-order response V_k' = ε_mM_k; the full second-order mean (noise part and parameters' part); the cross term; and E\|P_⊥v₁\|². It predicts Part D's estimator along u and along every eigenvector, and gives the window drifts by parts |
| `dirs_mc.py` | the real loop (K₀ + variational) on 8 antithetic pairs (seeds 0–7) at std 0.05, 0.1, 0.2: Part D's estimator along u and along all eigenvectors of DR's S_k, and the top finite-time exponents |
| `stationary.py` | the exact periodic covariance and second-order mean (noise and parameters' parts); the mean error and the separations' drift with the cross term at the stationary point, under the adaptive and the fixed prior; Part A's biases from the same forcing |
| `longtime.py`, `longtime_parts.py` | the exact second-order transient from the lift over 8.2·10⁴ periods (5.1·10⁵ time units). It is propagated period by period: P = P_per + ΦD_nΦᵀ, plus linear maps of one period. The parts version separates the noise part, the parameters' part and the cross term |
| `targets.py` | DR's drift at the lift for 11 exactly represented targets; (iv)'s centre-entry formula with its c″ term |
| `bound_ii.py` | (ii)'s bound against the drift, for 36 values of (k_e, B) |
| `derived.py` | the ratios, s.e. counts and thresholds quoted in this report, computed from the other scripts' outputs |
| `mutants.py` | 11 mutants of `measurement_noise.py`, patched in memory, against the D‴.33 test and the related tests |
| `from_raw.out`, `tests.out` | `exp_eiv_remedy.py --from-raw` and pytest on the scratch tree |

## Verdicts

| object | verdict |
|---|---|
| **(i) the drift** | **Holds as the mean Jacobian at the lift.** With the regressor at x_m the separation's Jacobian is exactly G^{1/2}ε_m∇²_ϑΦ̃(x_m; ϑ_k)G^{1/2}: on Part D's paths the variational equation and a 10⁻⁶ separation of the m = 15 lift agree to 3.9·10⁻⁸. Its second-order mean is E[e]M_k plus the cross term Σ_lE[ε_mδϑ_kl]∂_{ϑ_l}M_k. There is no compensation or leak term: the separation does not enter f̂ at first order, and the leak is −σ₀/k_e. At the lift DR's E[e] is MR's, since the regressor enters only the parameters' equation. Every printed S_k reproduces. **Two claims are not shown.** The cross term is negligible for MR, but not for DR and DR*: it moves the class-2 eigenvalue by −8% and −17% in the first window. And the exponent differs from the mean drift by E\|P_⊥v₁\|², +13% of DR's weakest direction (ID-5). |
| **(ii) the bound and the k_e-uniformity** | **The bound holds:** ‖S_k‖ is at most 0.40 of it in all 36 cases. It gives a k_e-uniform lower bound on the threshold (std ≥ 0.069 here). **"Does not fall with k_e" needs B ∝ k_e**, the md's B = 4k_e. At the instance's fixed B = 100, DR's drift is positive definite from k_e = 200, and DR's and DR*'s thresholds fall toward the bound's
floor. MR's "falls like 1/k_e in the variance" needs k_e ≪ B (ID-4). |
| **(iii) the stationary point** | **Does not hold at the stated order** (ID-1). The balance ⟨Ĵᵀ(g − Ĵȳ)⟩ = 0 omits the parameters' part of the force. D‴.31(iv)(c) shows that part is of the same order and not least-squares. With it, the exact periodic second-order mean error has rms 6.70·10⁻⁴ per unit Eν² (DR), against the stated 2.34·10⁻⁴. The omitted part contributes 6.32·10⁻⁴, and its tangent component is not removed by I − P. The drift there, with the cross term, is indefinite with largest eigenvalue +3.68·10⁻⁵, against "at most 4.7·10⁻⁶". "Zero when g lies in the tangent space" is false. |
| **(iv) the sign** | **The formula holds exactly**: the centre entry equals −(Eν²G_c/(2k_e))⟨∂²_xf∂²_xΦ̃_k(B − k_e − 2a)/(k_e + B)⟩ − c″G_γ⟨E[e]∂_xΦ̃_k⟩ to all digits in 11 targets. **The centre directions are negative in every target tried. The whole matrix's negative definiteness is the instance's:** it is indefinite for 6 of 10 other targets, among them three separated same-sign bumps (ID-3). **It holds at the lift only.** Along the trajectory DR's drift turns indefinite at t ≈ 2.7·10³, nine times the md's window, and this happens in the corollary's own model (ID-2). DR*'s reversed sign holds (positive definite at the lift). |
| proof | (i), (ii) and (iv) are correct; (iii)'s balance omits the parameters' part (ID-1). (iv)'s sign "proved" in the centre directions also needs the c″ term to be dominated: its sign is not controlled (ID-3). |
| Status | "(i), (iii) and (iv) exact consequences of the second-order expansion" fails for (iii) (ID-1). "(iv)'s sign proved for the centre direction" holds as above. The novelty list includes "the drift vanishing at the stationary point" (ID-1) and an uncited classical framework (ID-8). |
| *Numerics* | **Every printed number reproduces with independent code** (confirmations 3–10). Readings that do not hold: the stationary point (ID-1); "no noise level breaks the replicated configuration's attraction at this order" beyond the lift (ID-2); "DR's rates stay at theirs to 3.6%" (ID-6); "+5.0·10⁻⁵" (ID-10). |
| code | `mean_error`, `eiv_force(variant, return_g)`, `separation_drift(variant)`, `second_order_transient` and `separation_drift_series` are right for what they compute: S_k, g and the transient windows reproduce independently. `part_D_pred`'s stationary point is the averaged noise part only (ID-1). `separation_drift_series` drops the cross term on MR's justification (ID-5) and takes no variant (ID-10). |
| experiment, md, raw | **Hold.** The md regenerates byte-identically. Twenty of Part D's 32 stored runs reproduce with an independent simulator (confirmation 5). The prediction column is the lift's drift, and the start direction is DR's most favourable (ID-6). |
| test | **Passes** (8 tests in 33.5 s). It catches 7 of 11 mutants itself, and 10 of 11 together with the related tests. DR*'s drift doubled passes every test, and the test cannot see ID-1 (ID-7). |
| §0, §7 rows, §8 bullet, RESULTS, PAPERS 12, LITERATURE, README | **Numbers match the md**, with one rounding slip (ID-10). Wording: the stationary point and "the drift vanishing" (ID-1); "the noise stiffens the separations" and the trajectory (ID-2, ID-3); "does not fall with k_e" (ID-4); "at their predictions … to 4%" (ID-6); novelty (ID-8); the citation of D‴.32(vi) (ID-9). |

## Confirmed issues

**ID-1 (major; (iii) l. 3946–3948; proof l. 3962–3963 "the slow force balance is ⟨Ĵᵀ(g − Ĵȳ)⟩ = 0, i.e. Ĵȳ = ĴH⁻¹⟨Ĵᵀg⟩ = Pg";
Status l. 3967–3972 ("(i), (iii) and (iv) exact consequences …", "the drift vanishing at the stationary point" as new); *Numerics*
l. 3985–3987; §0 row l. 138, columns 2–4; §7 row l. 6434; §8 l. 6701–6703; RESULTS l. 824–825; LITERATURE l. 488–489; `exp_eiv_remedy.py`
l. 233–239 and 399–405; md l. 121–126; the test l. 263–268.)** The stationary-point claim omits two terms of its own order, and on the
instance they dominate. Statement (iii) and the numbers quoted for it do not hold at the stated order.

*The omitted terms.*
- The force. The full second-order force is F = ⟨Ĵᵀ(g + g_par)⟩ + k_e⟨f_y,par⟩. The parameters' part comes from the compensation's
  curvature and slope–parameter terms in the error equation, and from the law's ∇²_ϑΦ̃E[ε_mδϑ]. D‴.31(iv)(c) establishes it as of the same
  order and not least-squares. At the adaptive prior's stationary point it gives k_eE[e] = (I − P)(g + g_par) − ĴH⁻¹k_e⟨f_y,par⟩. The last
  term lies in the tangent space, and the projection does not remove it.
- The cross term. The drift there also carries the cross term Σ_lE[ε_mδϑ_kl]∂_{ϑ_l}M_k, which (i)'s proof includes. At the stationary
  covariance it is of the drift's order.

*The numbers* (`stationary.py`: the exact periodic covariance and second-order mean of the linearised K₀-atom loop, OU noise as a state,
k_e = 25, B = 100, adaptive prior, per unit Eν²; Nper = 2048 and 4096 agree to all printed digits).

The forcing is validated against the md: the same forcing gives Part A's biases. For DR the noise part is 0.470472 (md 0.47053), the full
bias 0.485839 (md 0.48584), the parameters' part 3.924·10⁻⁴ and its learning-law part 2.101·10⁻⁴ (md 3.92·10⁻⁴, 2.10·10⁻⁴). For DR* the
full bias is 0.048050 (md 0.04805).

| | DR | DR* |
|---|---|---|
| rms E[e], the md's averaged formula (I − P)g/k_e (my value) | 2.335·10⁻⁴ | 2.605·10⁻⁴ |
| rms E[e], exact: noise part | 2.106·10⁻⁴ | 1.366·10⁻⁴ |
| rms E[e], exact: parameters' part | 6.321·10⁻⁴ | 6.273·10⁻⁴ |
| rms E[e], exact: full | **6.697·10⁻⁴** (2.4% of the lift's 2.744·10⁻²) | **6.438·10⁻⁴** (7.0% of 9.18·10⁻³) |
| full E[e]: tangent part, orthogonal part | 5.798·10⁻⁴, 3.351·10⁻⁴ | 5.813·10⁻⁴, 2.767·10⁻⁴ |
| drift, noise part of E[e] only (the md's quantity) | ≤ 3.8·10⁻⁶ in magnitude | ≤ 2.3·10⁻⁶ |
| drift, parameters' part of E[e] | up to +2.59·10⁻⁴ (class 2) | up to +2.60·10⁻⁴ |
| cross term | down to −2.58·10⁻⁴ (class 2) | down to −2.59·10⁻⁴ |
| **drift, total: [min, max] per class** | [−4.35, +3.68]·10⁻⁵; [−0.13, +0.30]·10⁻⁵; [−2.82, +2.44]·10⁻⁵ | [−4.17, +3.44]·10⁻⁵; [−0.05, +0.20]·10⁻⁵; [−2.89, +2.56]·10⁻⁵ |

- *The tangent part* is √(2.101·10⁻⁴)/k_e = 5.80·10⁻⁴. That is the md's own Part A law-part response, so the omission is visible in the
  repository's numbers.
- *The md's figures, against the exact mean.* The full mean error is 2.9 times (DR) and 2.5 times (DR*) the md's figure. The drift's largest
  eigenvalue is 7.8 times the md's "at most 4.7·10⁻⁶ in magnitude". In class 2 the two omitted parts nearly cancel.
- *DR*'s averaged residual is not resolved even in its noise part.* The averaged slaving uses H and divides by k_e instead of the error
  channel's rate k_e − a. That is an O(|a|/k_e) error in a residual that is only 2.84% of |g| (‖(I − P)g‖/‖g‖; 0.85% for DR): the exact
  noise part is 1.37·10⁻⁴, against the averaged 2.60·10⁻⁴.
- *"Zero when g lies in the tangent space"* fails: the parameters' part remains, its tangent part alone giving rms 5.8·10⁻⁴.
- *The expansion applies here for DR.* The stationary whitened parameter covariance is at most 0.041Eν², and the mean deviation is
  15Eν² (physical units), small at std ≤ 0.05. The failure is therefore not an artefact of a breakdown, unlike MR's in round 22.
- *The fixed prior (the §8 Open item), computed the same way.* The residual mean error keeps 41% of the lift's (rms 1.127·10⁻²). DR's
  drift stays negative definite there (largest −3.8·10⁻⁵), and DR*'s stays positive definite (largest +1.50·10⁻³, threshold std 0.435).
  The cross term is up to −3.2·10⁻⁴.

*Fix:*
- (iii) and its proof. "At the adaptive prior's stationary point ⟨Ĵᵀ(g + g_par − Ĵȳ)⟩ + k_e⟨f_y,par⟩ = 0, so k_eE[e] = (I − P)(g + g_par)
  − ĴH⁻¹k_e⟨f_y,par⟩. The drift there is ⟨E[e]M⟩ plus the cross term. Both parameters' terms are of the order of the noise's residual, and on
  the instance they dominate."
- *Numerics*, the rows, RESULTS. Use the exact periodic mean, which `second_order_mean(periodic=True)` already computes: "rms E[e]
  6.7·10⁻⁴ (DR) and 6.4·10⁻⁴ (DR*) per unit Eν², 2.4% and 7.0% of the lift's, the parameters' part 6.3·10⁻⁴; the drift indefinite,
  |λ| ≤ 4.4·10⁻⁵".
- Drop "zero when g lies in the tangent space" and, from the novelty, "the drift vanishing at the stationary point".
- Code. Compute the stationary drift from `second_order_mean(periodic=True)`, with the cross term from its covariance.
- Test. Assert the stationary mean error from the full forcing (ID-7).

**ID-2 (major; the §0 row l. 138 "under DR the compensation's rectification of the noise makes the drift negative in the separated kernels'
centre directions, so the noise stiffens the separations"; §8 l. 6700–6703 "makes the drift negative definite on the instance (the noise
stiffens the separations), and along the trajectory the drift relaxes toward the residual of g outside the tangent space, zero at the
adaptive prior's stationary point"; *Numerics* l. 3980–3981 "DR's is negative definite in every class, so no noise level breaks the
replicated configuration's attraction at this order", and l. 3987–3990 (the window drifts); the §7 row "DR's window drifts negative definite
over 300 time units"; RESULTS l. 820, 825, 829; PAPERS l. 184–185.)** "The noise stiffens the separations" holds only near
the lift. Along DR's second-order transient the drift turns indefinite after about 2.7·10³ time units, in the corollary's own model.
The md's 300-unit windows are read as if they settled the trajectory, and the §8 description of the trajectory is wrong.

*The computation* (`longtime.py`, `longtime_parts.py`). It is the exact second-order transient from the lift, propagated period by period
over 8.2·10⁴ periods (5.1·10⁵ time units), and it contains every second-order term. The slowest collective mode relaxes in
k_e/h₁ = 8.2·10⁴ time units (h₁ = 3.06·10⁻⁴). Its first periods are consistent with my direct integration and with the md's windows:
- the md's [250, 300): −3.70·10⁻⁴, −1.39·10⁻⁴, −2.68·10⁻⁴;
- mine without the cross term: −3.707, −1.399, −2.690·10⁻⁴;
- the one-period average at t = 251, with the cross term: −3.745, −1.469, −2.724·10⁻⁴.

Largest eigenvalue of DR's period-averaged drift, per unit Eν²:

| t | class 1 | class 2 | class 3 |
|---|---|---|---|
| 0 (the first period) | −4.76·10⁻⁴ | −1.58·10⁻⁴ | −3.14·10⁻⁴ |
| 302 (the md's window ends) | −3.61·10⁻⁴ | −1.41·10⁻⁴ | −2.61·10⁻⁴ |
| 1257 | −2.03·10⁻⁴ | −5.055·10⁻⁵ | −1.13·10⁻⁴ |
| first positive | t = 4097 | t = 2947 | t = 2714 |
| 6283 | +4.75·10⁻⁵ | +3.69·10⁻⁵ | +6.56·10⁻⁵ |
| 32170 | **+1.05·10⁻⁴** | +3.93·10⁻⁵ | +9.02·10⁻⁵ |
| 5.1·10⁵ (≈ stationary) | +3.54·10⁻⁵ | +1.1·10⁻⁶ | +2.31·10⁻⁵ |

*The cause* is the noise part, not the terms of ID-1. With the noise part alone (k_eE[e] = g − Ĵȳ(t), ȳ the response to the noise's own
least-squares force) the drift turns indefinite at the same times: t = 4072, 2727, 2689. Its largest eigenvalue peaks at +1.24·10⁻⁴
(class 1, t ≈ 3.2·10⁴). Near the crossings the other two parts are small: the parameters' part's largest eigenvalue is at most
6.2·10⁻⁶ until t ≈ 6·10³, and the cross term's at most 1.6·10⁻⁵ until t ≈ 1.3·10⁴.

The mechanism follows from (iii)'s own formula.
- The fast modes of ȳ relax first: k_e/h_j = 128–1650 time units for the six modes with h_j ≥ 0.015.
- E[e] therefore passes through (I − P_relaxed)g/k_e, the part of g carried by the three slow modes.
- That residual has no sign relation to ∂²_xf, which is what (iv)'s sign rests on.
- So "k_eE[e] moves from g … to (I − P)g" describes the end points only. The drift between them is not monotone.

*What survives.* The positive eigenvalues reach at most 1.05·10⁻⁴ per unit Eν² (1.24·10⁻⁴ for the noise part alone). They give a formal
threshold std of 1.64 (1.51), far outside the expansion. DR still keeps the replicated configuration attracting at small noise, by a
wide margin over MR. The sign claim does not survive: after t ≈ 2.7·10³ some separations decay slower than σ₀/k_e under the noise.

*Fix:*
- (iv), the §0 row, the §8 bullet, PAPERS, RESULTS: "near the lift (on the instance for t ≲ 2.7·10³) the noise stiffens the
  separations; along the trajectory the drift turns indefinite, its largest eigenvalue reaching 1.05·10⁻⁴ Eν² (t ≈ 3·10⁴), before settling
  at the stationary point's (ID-1)".
- Replace §8's trajectory sentence.
- Part D: carry the transient beyond 300. Period-by-period propagation makes it cheap: the linear maps take about a minute. Print the time
  of the sign change.

**ID-3 (minor; *Numerics* l. 3980–3981; the §0 row l. 138, columns 2–3; (iv) l. 3954 "there a noise fast against the loop
(B > k_e + 2 max|a|) stiffens the separations"; PAPERS l. 185; RESULTS l. 820, 829.)** The whole matrix is negative definite only on the
instance. The text proves the centre direction and checks the rest numerically on one target, then states generally that the noise
stiffens the separations.

*Evidence* (`targets.py`; quasi-static DR drift at the lift, per unit Eν²; 11 exactly represented targets).
- *The centre entry* equals (iv)'s formula with the c″ term to all printed digits, and it is negative in every target, overlapping ones
  included.
- *Negative definite* at (k_e, B) = (25, 100) only for separated alternating-sign kernels with centres well inside the orbit: the
  instance, two opposite bumps, four alternating bumps, narrow bumps (S = 0.32), and unequal weights.
- *Indefinite* for 6 of the 10 other targets, among them alternating-sign kernels that are wide or sit at the orbit's edge, and every
  same-sign target:
  - three separated bumps of the same sign, (0.8, −1.5, 0.5), (0.8, 0, 0.5), (0.8, 1.5, 0.5): largest eigenvalues +8.7·10⁻⁶,
    **+2.86·10⁻⁴**, +1.08·10⁻⁵ (class 2's width entry +1.29·10⁻⁴);
  - wide bumps (S = 1): +2.27·10⁻⁴, +5.55·10⁻⁴, +3.40·10⁻⁴;
  - bumps at the orbit's edge (C = ±2.5): +2.37·10⁻⁴ and +2.69·10⁻⁴ (classes 1 and 3);
  - same-sign overlapping bumps, four same-sign overlapping bumps, and two same-sign close bumps (+5.03·10⁻⁴).
- The same pattern holds at B = 1000 and at (k_e, B) = (100, 400). The positive directions are width and amplitude.
- *The saturation term's sign is not controlled*: −1.98·10⁻⁴ (instance, class 1), +2.43·10⁻⁴ (two bumps), −3.64·10⁻⁴ (edge bumps,
  51.8% of the entry). So "proved in the centre directions" holds where the main term dominates it.
- *At fixed B the sign reverses beyond k_e ≈ B.* At B = 100 DR's drift is positive definite at k_e = 200 (λ_min +1.06·10⁻⁵,
  λ_max +1.03·10⁻³). (iv) states the condition; the §0 row and PAPERS do not.

*Fix:* in the §0 row, PAPERS and RESULTS: "negative definite on the instance (separated alternating-sign kernels); the centre entries negative for
separated kernels whose saturation term is dominated; for same-sign or wide kernels the width directions soften (three separated
same-sign bumps: +2.9·10⁻⁴ per unit Eν²); a noise slower than the loop reverses the sign". Change "no noise level breaks … at this order" to
"at the lift, on the instance". Add a second target to Part D.

**ID-4 (minor; (ii) l. 3943–3944 "the noise level up to which the replicated configuration attracts near the lift does not fall with k_e;
under the measured regressor it falls like 1/k_e (in the variance)"; the §0 row l. 138; Status l. 3971; *Numerics* l. 3983–3985; §8
l. 6700; RESULTS l. 822–823; PAPERS l. 185 and 203; LITERATURE l. 487; README l. 66–67 "whose threshold does not fall with the gain"; md
l. 112.)** Both claims hold only when the noise's bandwidth is scaled with the gain, as in the md's scan at B = 4k_e. Neither holds for a
sensor of fixed bandwidth such as the instance's B = 100. (ii) does not state this hypothesis.

*Evidence* (`lift_drift.py`, `bound_ii.py`).
- *At B = 4k_e* the md reproduces exactly and extends to k_e = 800. k_e·λ_max is −0.0038 (DR) and 0.1242 (DR*), and MR's threshold std
  goes 0.0308 → 0.0054, like k_e^{−1/2}.
- *At B = 100*, threshold std (— where negative definite):

  | k_e | 25 | 50 | 100 | 200 | 400 | 800 |
  |---|---|---|---|---|---|---|
  | DR | — | — | 2.58 (indefinite) | 0.185 | 0.138 | 0.121 |
  | DR* | 0.238 | 0.185 | 0.151 | 0.131 | 0.119 | 0.113 |
  | MR | 0.0308 | 0.0237 | 0.0193 | 0.0167 | 0.0152 | 0.0144 |

  MR's threshold variance falls by 4.6 over a factor 32 in k_e, not like 1/k_e: for k_e ≫ B, E[ε_m²] → Eν²B/k_e, which D‴.30(iii) states for
  the bias.
- *The bound holds.* ‖S_k‖ is at most 0.40 of it in all 36 cases (B = 4k_e, 100, 1000; k_e = 25–800). It gives a k_e-uniform lower bound
  on the threshold, std ≥ 0.069 on the instance for both DR and DR*, rather than invariance.

*Fix:* (ii): "the threshold is bounded below uniformly in k_e (on the instance std ≥ 0.069), and constant when B/k_e is fixed. For a sensor
of fixed bandwidth DR's sign reverses at k_e ≈ B, and the thresholds fall toward the bound's level (B = 100: DR 0.185 → 0.121, DR*
0.238 → 0.113 for k_e = 200 → 800). The measured regressor's falls like 1/k_e in the variance while k_e ≪ B, and saturates for k_e ≫ B."
Add a fixed-B scan to the md.

**ID-5 (minor; (i) l. 3937–3939 "S_k = ⟨E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩"; proof l. 3960–3961 "the parameters' cross term … is as in Proposition
D‴.30(i)", which calls it negligible for MR; the §0 row "the mean error's term … alone"; the thresholds 0.238 and 0.138, l. 3982 and
RESULTS l. 820–821; `separation_drift_series` l. 580–581 "omitted (on the instance below 10⁻³ of the drift, referee round 22)";
`part_D_pred` l. 242–247.)** Two terms are negligible against MR's drift but not against DR's and DR*'s, which are 60–2000 times smaller.

*The cross term* (`so_sep.py`, `windows_parts.py`). It builds up within the first period. Its spectral norm is 7.6·10⁻⁴ per unit Eν² in
class 2 over [0, 50), and 3.7·10⁻⁴ over [250, 300). MR and DR share the linearisation, so it is the same for both, and DR*'s is nearly
the same. "Below 10⁻³ of the drift" is MR's ratio (0.25% of 0.298).
- *DR:* the cross term's norm is 2.6–5.1 times DR's class-2 top eigenvalue. It moves λ_max(S₂) by −8% over [0, 50) (−1.490 → −1.615·10⁻⁴)
  and by −4% over [250, 300).
- *DR*:* it moves the class-2 eigenvalue by −17% over [0, 50) (4.449 → 3.690·10⁻³). DR*'s threshold std near the lift is then 0.277,
  with the transient's E[e] and the cross term over the first 50 time units. The stated value is 0.238 (quasi-static moments; 0.246 with
  the exact periodic moments).

*The variance term.* The exponent is not the mean drift: to second order E log‖u‖ = −st + ∫θᵀE[A₂]θ + E\|P_⊥v₁\|². The last term is
nonnegative. It is +1.96·10⁻⁵ per unit Eν² along DR's weakest direction (class 2's top eigenvector), 13% of DR's drift there with the
opposite sign; over all directions it ranges from +1.3·10⁻⁶ to +3.2·10⁻⁵. D‴.30(i)'s "negligible cumulant" is MR's. On the instance
neither term changes a sign (confirmation 7).

*Fix:*
- (i): "S_k at the lift; once the parameters fluctuate (within a period) the cross term adds, not negligible against DR's and DR*'s drift
  (−8% and −17% of the class-2 eigenvalue over the first 50 time units); the exponent adds the nonnegative E\|P_⊥v₁\|²/t".
- DR*'s threshold: "0.238 at the lift's quasi-static moments, 0.28 over the first 50 time units".
- Code: add the cross term to `separation_drift_series`, since the transient's covariance is available, or restrict its docstring to MR.

**ID-6 (minor; *Numerics* l. 3990–3996 "(0 of 16 runs growing) — faster than σ₀/k_e … DR's rates stay at theirs to 3.6%"; the §0 row "at
their predictions along the start's direction to 4%"; RESULTS l. 827–829; `exp_eiv_remedy.py` l. 411–425 and `job_D` l. 259 (the start
along MR's top eigenvector); md l. 139–146.)** Part D compares DR's rates with the lift's drift, though the drift moves along the run. It
also tests DR in its most favourable direction.

*(a) The prediction column misses by about 18%, and the transient explains it.*
- The column is −σ₀/k_e + Eν²uᵀS₃u, the lift's drift. It misses the measured noise-induced part of DR's rate:

  | std | measured noise part | predicted (lift) | ratio | stored pair means' s.e. |
  |---|---|---|---|---|
  | 0.05 | −1.353·10⁻⁵ | −1.655·10⁻⁵ | 0.82 | 22.5 s.e. off |
  | 0.1 | −5.352·10⁻⁵ | −6.620·10⁻⁵ | 0.81 | 23.7 s.e. off |

- "To 3.6%" divides by the whole rate. At std 0.05, 95% of that rate is the noise-free σ₀/k_e (84% at 0.1).
- *The transient explains the miss* (`run_so_sep.py`). Along u, the actual mean error's term is −5.520·10⁻³ per unit Eν², against the
  lift's uᵀS₃u = −6.620·10⁻³; the cross term adds −3.4·10⁻⁵ and E\|P_⊥v₁\|² +1.7·10⁻⁵. The predicted rates are −2.96817·10⁻⁴ and
  −3.38344·10⁻⁴.
  - Against the stored runs: within 2.3 and 3.5 s.e.
  - Against my own 8-pair runs: within 1.5 and 2.5 s.e. (−2.965827·10⁻⁴ ± 1.6·10⁻⁷, −3.367922·10⁻⁴ ± 6.3·10⁻⁷).
  - The md's runs are thus at the second-order prediction once the transient is used; the lift's prediction is 22.5 s.e. away.
- The miss is (iii)'s own mechanism, the mean error relaxing along the trajectory. The md computes it for its windows but not for the
  estimator.

*(b) The start direction is DR's most favourable one.* The start u is within 18° of DR's most stable eigenvector of S₃ (\|cos\| = 0.952).
uᵀS₃u = −6.62·10⁻³ is 13–44 times DR's top eigenvalues (−4.97·10⁻⁴, −1.51·10⁻⁴, −3.36·10⁻⁴).
- Along DR's top eigenvectors (`dirs_mc.py`, 8 pairs, std 0.05) the rates are −2.839655·10⁻⁴, −2.833422·10⁻⁴, −2.837043·10⁻⁴
  (± 2–3·10⁻⁸).
- These are faster than σ₀/k_e by 0.35%, 0.13% and 0.26%, against 4.8% along u. The noise's "stiffening" was measured where it is
  strongest.

*(c) The md's column "(the rate the separation tends to once it has turned onto the top eigenvector)".* Turning takes
1/(Eν²(λ_top − λ_bottom)) ≈ 5.8·10⁴ time units at std 0.05 (1.5·10⁴ at 0.1). That is beyond the run, and beyond ID-2's sign change.

*(d) "(0 of 16 runs growing) — faster than σ₀/k_e"* holds for the pair means. Individually, 5 of the 16 stored runs decay slower than σ₀/k_e
(e.g. −1.649·10⁻⁴ at std 0.1, seed 2, sign +1). In my runs at std 0.05 it is 3–8 of 16 per direction. None grows.

*(e) The §0 row's "at their predictions … to 4%"* also covers MR, which at std 0.1 grows at less than half its prediction.

*Fix:*
- Compare DR with the transient's prediction along u, print λ_eff with its standard error, and state the ratio 0.82.
- Start DR also along its own top eigenvectors of every class, where the effect is 0.1–0.4% of σ₀/k_e.
- Drop the "turned" column, or state the turning time.
- Write "the pair means decay faster than σ₀/k_e; no run grows".

**ID-7 (minor; `tests/test_eiv_remedy.py` l. 243–268; Status l. 3973–3975; RESULTS l. 831–833; the §7 row l. 6434, column 3.)** The test
catches the gross errors in DR's path. It leaves DR*'s magnitude and every claim beyond the lift untested.

*Mutants* (`mutants.py`). The tree's 8 tests pass in 33.5 s.

| mutant of `measurement_noise.py` | the D‴.33 test | the related tests |
|---|---|---|
| the ∂_x term kept for DR and DR* (M1); DR's E[e] = DR*'s (M2); DR's sign flipped (M3); DR's factor (B − k_e − 2a)/(k_e + B) dropped (M4); DR ×0.85 (M6); DR's factor at B → ∞ (M9); DR's whitening dropped (M10) | caught (l. 257) | — |
| `mean_error`'s shared branch ×0.85 (M7); the sign of a flipped there (M8) | **pass** | caught (quadrature l. 183 of `test_measurement_noise.py`; `second_order_mean` l. 233) |
| DR*'s mean error with k_e + a (M5) | **pass** | caught (l. 233) |
| **DR*'s drift doubled in `separation_drift` (M11)** | **pass** | **pass** |

- The test's identity assertion compares two code paths that share `mean_error` and the same Hessian code. That is why M7 and M8 pass it;
  the related tests catch them because they type E[e] independently.
- For DR* the test checks only its sign and its k_e-ratio, so the threshold 0.238 is untested (M11).
- The stationary assertion (the rms of g − ĴH⁻¹k_eF below 1% of g; the value is 0.85%) tests the averaged noise part. That is a property
  of the instance, and the assertion cannot see ID-1.
- Nothing tests the transient's sign for DR, the k_e-uniformity at fixed B, or Part D's real loop.

*Fix:*
- Test DR*'s drift against an independent quadrature, as `test_measurement_noise.py` l. 159–185 does for MR.
- Test the stationary mean error from `second_order_mean(periodic=True)` with the full forcing.
- Test DR's transient window drifts, and Part D's DR rate against the transient's prediction.

**ID-8 (minor; Status l. 3969–3973; the §0 row, column 4; PAPERS l. 201–203; LITERATURE l. 484–490.)** The audit does not cite the classical
framework for the corollary's question, and two of its "new" items do not hold.
- *The framework.* The attraction of a replicated configuration under a noise common to its members is the transverse stability of a
  noise-synchronized state, and that is classical. Identical uncoupled systems driven by a common noise synchronize when the transverse
  Lyapunov exponent, computed to second order in the noise, is negative. References: J. Teramae, D. Tanaka, PRL 93:204103 (2004);
  D. S. Goldobin, A. Pikovsky, PRE 71:045201 (2005); H. Nakao, K. Arai, Y. Kawamura, PRL 98:184101 (2007). Bibliographic details were
  checked by a web search; the content is from the abstracts. LITERATURE cites Carmona–Delarue–Lacker for the non-averaging of a common
  noise, but not this.
- *The four "new" items:*
  - the drift reduced to the mean error's term is the one-line consequence of a y-independent regressor;
  - the k_e-uniform threshold is the bound's lower bound (ID-4);
  - the stabilising sign holds in the centre directions, and for the whole matrix only on the instance and near the lift (ID-2, ID-3);
  - the drift vanishing at the stationary point does not hold (ID-1).

*Fix:* cite the synchronization works. Restrict "new" to the removal of the regressor's ∂_x term and the instance's quantification.

**ID-9 (minor; the setting l. 3935 "σ₀w_k/k_e under the cohesion leak (Proposition D‴.32(vi))"; Status l. 3972–3973; *Numerics* l. 3982;
§8 l. 6703; LITERATURE l. 490.)** The corollary rests on a statement that is not in the document.
- *The citation is dangling.* PARADIGM.md at bfdbc0f contains no Proposition D‴.32; its only occurrence is this citation. The item "(vi)"
  exists only in `tests/test_cohesion_leak.py` (docstring l. 9: in the noise-free real loop a separation decays by e^{−σ₀w_kT₀/k_e} per
  period). The stiffness itself is in the docstring of `cohesion_leak.py`.
- *The transfer is consistent* (`lift_drift.py`).
  - With the compact kernel and classes at least R apart, the cohesion force vanishes on the replicated configuration and on its class
    means.
  - Its linearisation on a class's separations is −σ₀w_k/k_e, and the noise's S_k does not depend on the leak.
  - So the thresholds scale by √w_k: 0.1376 = 0.2384/√3 and 0.0178 = 0.0308/√3, as printed.
  - The class means feel no leak, so the stationary balance is the adaptive prior's, and ID-1 applies to it.

*Fix:* write D‴.32, or cite the code and its test for the stiffness, or derive σ₀w_k/k_e in a line.

**ID-10 (cosmetic.)**
- DR*'s smallest eigenvalue in class 2, "+5.0·10⁻⁵" (l. 3980, the §7 row, RESULTS l. 819), is 5.11·10⁻⁵ (raw 5.10·10⁻⁵): write
  "+5.1·10⁻⁵".
- (iv)'s "with the metric weight G_c": G acts on γ, so G_c := G_γc′(γ)² needs defining.
- `separation_drift_series` takes no variant. For DR and DR* the caller must pass E[ε_m²] = 0, as `part_D_pred` does (l. 237, 243);
  passing the transient's own Eem2 would put MR's term back. Say so in its docstring.
- md l. 112 "at B = 4k_e (the noise white over the loop's bandwidth)": at B = 4k_e the factor (B − k_e − 2a)/(k_e + B) is about 0.6, not 1.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` on the scratch tree rewrites `eiv_remedy.md` byte-identically: md5 `34623bdf…`, equal to
   `git show bfdbc0f:…`. `tests/test_eiv_remedy.py` passes, 8 tests in 33.5 s.
2. **The instance** (`check_inst.py`): f̂ = f to 2.2·10⁻¹⁶, f̂′ = f′ and f̂″ = f″ to 1.3·10⁻¹⁵. The closed-form gradient and Hessian equal the
   jets to 3.6·10⁻¹⁵, and the jets' third derivatives equal central differences of the closed forms to 2.2·10⁻⁸. The atoms (a, c, s) are
   (3, −1.5, 0.5), (−2.1, 0, 0.4), (2.4, 1.5, 0.5); c′ = 0.75, 1, 0.75 and c″ = 0.25, 0, −0.25.
3. **The moments, re-derived.**
   - E[eν] = −(k_e + a)Eν²/(k_e + B) and Var(e) = (k_e + a)²Eν²/(k_e(k_e + B)).
   - E[ε_m²] = Eν²(a² + k_eB)/(k_e(k_e + B)), by partial fractions of the high-passed OU spectrum.
   - E[e] = −(∂²_xf/(2k_e))(2E[eν] + Eν²) = −(∂²_xf/(2k_e))Eν²(B − k_e − 2a)/(k_e + B), the same for MR and DR at fixed parameters.
   - DR*: E[e] = ½∂²_xf·Var(e)/(k_e − a), with Var(e) = k_e²Eν²/((k_e − a)(k_e − a + B)).
   - All agree with `mean_error`'s formulas.
4. **S_k at the lift** (`lift_drift.py`). Every entry of md table D reproduces: MR [−0.37144, +0.01481]; [−0.02570, +0.01009]; [−0.03507,
   +0.29815]; DR [−0.01122, −0.00050]; [−0.01490, −0.00015]; [−0.00719, −0.00034]; DR* [+0.00016, +0.00378]; [+0.00005, +0.00498];
   [+0.00011, +0.00238].
   - Thresholds: 0.0308, 0.2384, and with the cohesion leak 0.0178, 0.1376.
   - u = (−0.2781, −0.9606, 0.0032); uᵀS₃u is −6.6204·10⁻³ for DR, as the raw `rq`.
   - The B = 4k_e scan: all 12 values. N = 512 and 2048 agree to 10⁻¹⁵.
   - The exact periodic noise moments at fixed parameters give rms E[e] 2.705·10⁻² against the quasi-static 2.744·10⁻²; pointwise they
     differ by up to 19% of the maximum. The eigenvalues move by at most 6% (DR's class-2 top −1.433·10⁻⁴ against −1.514·10⁻⁴; DR*'s
     largest 4.680·10⁻³ against 4.979·10⁻³), inside (ii)'s (1 + O(1/k_e)).
5. **Part D reproduced** with the independent simulator, on Part D's noise paths (`repro_partD.py`, `check_std01.py`).
   - DR at std 0.05, seeds 0–3, both signs (m = 15 lift): all 8 rates to 7 digits, separation norms to ≤ 3.3·10⁻⁹ relative, pair means
     −2.967911, −2.962084, −2.963592, −2.966624·10⁻⁴.
   - MR at std 0.05, seeds 0–1: 4 rates to 7 digits, norms ≤ 1.7·10⁻⁹.
   - DR at std 0.1, seeds 0–3 (linearised): 8 rates to 7 digits, norms ≤ 3.0·10⁻⁸.
   - The K₀-atom loop with the exact variational equation equals the finite separation of the lift to ≤ 3.9·10⁻⁸.
6. **The second-order transient** (`run_so_sep.py`, `windows_parts.py`).
   - DR's window drifts in md l. 132–137 reproduce to ≤ 0.9%: [0, 50) −4.598, −1.490, −3.217·10⁻⁴ against −4.59, −1.48, −3.21·10⁻⁴;
     [250, 300) −3.707, −1.399, −2.690·10⁻⁴ against −3.70, −1.39, −2.68·10⁻⁴.
   - The window rms of E[e] reproduces to ≤ 1%: 0.0265 → 0.0230 against 0.0264 → 0.0228.
   - MR's windows in `measurement_noise.md` reproduce to four digits (0.3066/0.3069 … 0.3494/0.3494; rms 0.1213 … 0.6024).
7. **The real loop along every direction** (`dirs_mc.py`: 8 antithetic pairs at each std 0.05, 0.1, 0.2, so 48 paths; 10 directions per
   path).
   - DR's pair means decay faster than σ₀/k_e in all 10 directions at all three stds. No run grows.
   - They agree with the second-order transient in every direction to ≤ 2.6 s.e. at std 0.05 and 0.1. Examples: class 2's top eigenvector
     −1.467·10⁻⁴ ± 5.9·10⁻⁶ per unit Eν² measured against −1.344·10⁻⁴ predicted; class 1's −3.961·10⁻⁴ ± 1.1·10⁻⁵ against −3.959·10⁻⁴.
   - At std 0.2 the next order appears (along u −5.153·10⁻³ against −5.537·10⁻³). The pair-mean top finite-time exponents are slower than
     σ₀/k_e, through the max over directions of the first-order part, which falls like t^{−1/2}.
   - So DR's sign holds on the instance for 300 time units, as stated.
8. **The averaged stationary point** (`stationary.py`): rms (I − P)g/k_e 2.33523·10⁻⁴ (DR) and 2.60497·10⁻⁴ (DR*), as the md;
   ‖Pg‖/‖g‖ = 0.999964 and 0.999597; H's eigenvalues 3.06·10⁻⁴ … 0.195. Part A's biases also reproduce with the same forcing (ID-1).
9. **(iv)'s formula** (`targets.py`): the centre entry equals the formula with the c″ term exactly in all 11 targets. The B ≫ k_e
   decomposition's diagonal term dominates the cross terms on the instance (6.50, 4.99, 4.16 against +0.024, −0.023, +0.027).
10. **(ii)'s bound** (`bound_ii.py`): ‖S_k‖ ≤ 0.40 of the bound in all 36 cases, with ‖∂²_xf‖_∞ = 5.016 on the orbit and
    ⟨‖M_k‖⟩ = 0.466, 0.588, 0.399.
11. **DR*** along its transient (`longtime.py`) stays positive definite at the lift, its largest eigenvalue decreasing 3.6·10⁻³ → 1.3·10⁻³
    (t ≈ 2500) → 3.5·10⁻⁵ at the stationary point. The stated threshold 0.238 is the lift's, and the conservative end.
12. **Literature** (web search, bibliographic details only): Teramae–Tanaka, PRL 93:204103, 2004 (doi 10.1103/PhysRevLett.93.204103);
    Goldobin–Pikovsky, PRE 71:045201, 2005; Nakao–Arai–Kawamura, PRL 98:184101, 2007.

## Not checked

- The real loop beyond 300 time units. ID-2's sign change is computed in the second-order expansion only. Measuring it means resolving a
  0.1–0.4% change of the rate at t ≥ 3·10³, which needs about 10⁴ time units per path. The stationary point in the real loop (about 10⁵
  time units) was not run either.
- Twelve of Part D's 32 stored runs were not reproduced individually: MR at std 0.05 for seeds 2–3, and MR at std 0.1. Their pair means are
  quoted as stored.
- DR*'s real loop (Part D has no DR* runs), and k_e ≠ 25 in the real loop. The fixed prior and finite η were computed only at the
  stationary point. The cohesion leak was not run with noise.
- The targets of ID-3 were checked at the lift only (no transient, no real loop).
- n > 1, and the content of the synchronization literature beyond its abstracts.
- Lemma D‴.26a, D‴.30 and D‴.31 were taken as given. Only `tests/test_eiv_remedy.py` was run with pytest; the related tests of
  `test_measurement_noise.py` were run inside `mutants.py`.
