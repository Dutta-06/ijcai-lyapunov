# Referee D3af (round 33): report on Proposition D‴.45, the free mass–amplitude law as a mirror flow

## What I checked and how

**Setup and rules followed**
- I worked only in the scratch copy `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review33/repo`, copied at 09:06 with the uncommitted results included. I never touched the repository and ran no git command that changes anything.
- Every run used `OMP_NUM_THREADS=1` and `nice -n 10`. I ran at most one heavy process at a time, never used `pkill -f` or `pgrep -f`, and never set `EXP_CACHE`.
- During the review the author committed D‴.45's results (23b4b80) and filled D45C4RED. `mirror_free.py`, `exp_mirror_free.py` and `test_mirror_free.py` are byte-identical to my copy.
- Scripts and logs are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review33/work/`.

**By hand**
- (i): I re-derived the conservation law dQ = 0 and dw/dα > 0, and worked out the branches.
- (ii): I re-derived the three closed forms, both branches, the κ≈1 identity, dw/dη = mD, ψ* and ψ.
- (iii): I checked the dual flow, its Hessian and the coercivity argument.
- (iv)(b): I re-derived the second-order ΔQ. I also showed its identity with D‴.44's drift: the flat part is (γβ|φ|²/mN)[γa² − 2Gβa′(1+a²/a²_max)] and the Itô part is γ²βa²|φ|²/(mN).
- Also re-derived: (iv)(c)'s Itô term, (v)'s V̇ and (vi).

**Tests**
- 8 passed in 59–71 s. D‴.46's leak test accounts for 54 s; D‴.45's own 7 tests take about 7 s.

**Raw data**
- I read `mirror_free_raw.npz`. The B residuals and the per-instance values behind the medians are stored there but not printed in the report.

**Numerical checks (script → what it does)**
- `eta_edge.py`: q = ±c; the negative branch; a 60-digit Decimal reference up to q = 10¹⁰; the range of `alpha_of_eta`.
- `shape.py`: the ℓ² core, the wall, the ℓ¹ weight and its remainder.
- `c2_exact.py`: the exact E[ΔQ] by Gauss–Hermite at h = 0.5, 0.05 and 0.005.
- `rich_check.py` and `dual_solver.py`: the 7 instances whose unequal-mass rich distances exceed 10⁻³, re-solved with an independent Newton method on the dual objective 𝒟 with Armijo backtracking.
- `branch_flow.py`: the loop itself from three branch starts, run to t = 10⁷.
- `branch_feas.py`: whether an exact fit with the frozen signs exists.
- `branch_cost.py`: the marginal costs η(w) − η(0) of branch and free particles.
- `c4_reduced2.py`: D‴.44's reduced flow to s = 270, called exactly as `part_C4` calls it (D45C4RED).
- `c3_hard.py`: the reduced selection flow to s = 10⁸ on two instances where basis pursuit ≠ truth.
- `gauge_eq.py`: C(3)'s end values of Q against the zero of (iv)(b)'s bracket.
- `flow_time.py`: the loop's time to reach the projection at s = 10⁻², 10⁻³ and 10⁻⁴.
- `d_split.py`: Part D's split error.

**Mutation testing**
- `mutate.py`, `mini_exp.py` and `mut_analyse.py`; log in `mutate_all.log`.
- 36 single-line mutants of `mirror_free.py`, each run against D‴.45's 7 tests and against a 14-second reduced copy of the checks in Parts A–E.

**Literature**
- Web search only (arxiv.org unreachable). Items marked "(memory)" below are not verified.

## Verdicts

### (i) The conservation law
- **Correct as an identity, for every scalar drive.** w is strictly increasing and unbounded along the curve; the branch geometry is right.
- **β = 0 is never reached in finite time.** β(t) = β(0)·exp(γ∫a·g dt) stays positive for a locally bounded drive.
- **Overstated scope.** The paradigm's position vector is ϑ = (α, γ, κ), and the cohesion term −(σ₀/m)Σψ_ij(ϑ_i − ϑ_j) acts on α_i − α_j. That gives dQ_i = −(2c/a_max)·sinh(2α_i/a_max)·(the cohesion's α-force), which is nonzero whenever two particles are within range — always under the Gaussian ψ.
- Moving centres can bring classes into range, so "also when centres and widths move" needs a qualifier. Two particles on one atom (Part D, the pair in (v)) are inside the range by construction.

### (ii) The closed form of η, the metric and the potential's shape
- **The closed forms are correct in all cases, on both branches.** They match quadrature to 1.5·10⁻¹⁴; the negative branch is checked neither by the tests nor by Part A. At q = c exactly, η = α/(2Gc). At q = −c exactly, the code returns NaN.
- dw/dη = mD, ψ* and ψ are correct.
- **The κ≈1 identity is right, but the code is not cancellation-free.** The stable form is used only for |t| < 0.9. For |t| ≥ 0.99 the relative error is 5·10⁻¹¹ at q = 10⁴, 8·10⁻⁷ at q = 10⁸ and 4·10⁻⁵ at q = 10¹⁰.
- **Shape claims hold:**
  - the ℓ² core;
  - the soft wall: η stays within 0.71 (β⁰ = 100c) and 0.095 (β⁰ = 1000c) of (1/γa_max)·log(|w|/w_⊓) for |w| ≥ 1.5·w_⊓;
  - the factor β⁰/(4c);
  - the rich ℓ¹ weight: its approximation is off by 10⁻³ at β⁰ = 10⁻⁴.
- **The stated remainder O(|w|(1+|log|w||)) is false.** Outside the core, R = ψ − λ₁|w| = −(3/2)·K·|w|^{2/3}·(1+o(1)), with K = (a²_max/2Gc)(2c/(a²_max·m))^{1/3}, uniformly in β⁰.
  - At β⁰ = 10⁻⁸, R/(|w|(1+|log|w||)) is −17, −61 and −222 at |w| = 10⁻⁵, 10⁻⁷ and 10⁻⁹.
  - Meanwhile R/(−1.5·K·|w|^{2/3}) stays between 0.994 and 1.000.

### (iii) The implicit bias
- **The convergence proof is sound.** It is the gradient flow of a strictly convex, coercive dual; full row rank suffices, and nearly collinear atoms only slow it down. The limit is an exact fit.
- **The lazy Γ-limit is sound.** The rich one is sound with the corrected remainder, but only when weighted basis pursuit is unique — and the statement omits that hypothesis.
- **The rich limit is reached at scales set by the LP's vertex gap.** The dual solver puts k=3 i=13 at 0.17 from weighted basis pursuit at s = 10⁻⁵ and 10⁻⁷ (excess weighted ℓ¹ 8.3·10⁻⁵), 4.3·10⁻³ at s = 10⁻⁹, and 3·10⁻⁸ only at s = 10⁻¹¹.
- **The time cost grows like 1/β⁰².** The loop reaches the projection to about 10⁻⁸ at t ≈ 10⁶ (s = 10⁻²), 10⁷ (s = 10⁻³) and 10⁹ (s = 10⁻⁴). The experiment never integrates the loop below s = 10⁻².
- **Branches.**
  - The direction of the claim ("no longer basis pursuit") is right.
  - "Their potential an entropy against the others' ℓ¹ weight" is wrong. At β⁰ = 10⁻⁴, branch particles have η − η(0) of 340–480 against 420–1220 for free particles, with λ₁ = 702 — the same order. The limit is a sign-asymmetric weighted ℓ¹ (weights λ₁ ∓ η_i(0)) with sign constraints.
  - The exact fit can disappear: with frozen signs, 2/20 (k=2) and 1/20 (k=3) instances have no exact fit at s = 10⁻⁶; 3/20 and 2/20 at s = 10⁻⁸.
  - B's branch numbers are misreported (D3af-3).

### (iv)(a) Noise in continuous time
- **True pathwise for noise with smooth paths, i.e. Stratonovich white noise:** every colour, every node covariance, every amplitude.
- **It is consistent with D‴.41(ii) and with Pesme–Pillaud-Vivien–Flammarion, but the draft doesn't say why.**
  - D‴.41(ii)'s requirement C ∝ I comes from the conserved law's centring term, which the free law does not have.
  - Pesme et al.'s stochastic gradient flow is an Itô SDE (memory). Its drift correction is exactly the second-order ΔQ of (iv)(b), so their effect is the sampled-data one, not a counterexample.
- **The pieces drop "Stratonovich"** (ROW0, BULLET, README), and as written they are false.
- **Caveat.** What is invariant is the limit map. Under persistent noise the stationary mean E[w] is an exact fit (E[∇𝒟] = 0) displaced from the projection at O(v).
- **C(1) cannot test the interpretation.** A noise that is a smooth ODE input is Stratonovich by construction, so C(1) only checks the integrator.

### (iv)(b) The Euler step and the mirror update
- **The ΔQ formula and its sign are correct.**
- **C(2)'s 14–25 standard errors are the O(h³) term, not an error.** The exact means differ from the closed form by 3.2–3.6% at h = 0.5, 0.29–0.35% at h = 0.05 and 0.03–0.04% at h = 0.005 — linear in h — and the Monte Carlo means agree with them.
- **The identification with D‴.44's drift is exact.** C(3)'s used atoms end exactly at the bracket's zero, γ·a³_max·t³ = G·m·|w|·(1 − t⁴): end values of Q such as −1.863221·10⁻⁴ are reproduced to 7 digits in all six instances.
  - That zero is not D‴.44(ii)'s gauge-minimal relevance D*, which has 2G and would put these atoms at Q ≈ −2 to −3.3.
  - D‴.44's route to basis pursuit is the flat scheme's, so nothing supports basis pursuit under the paradigm's update.
- **"Reaches the rich limit (basis pursuit)" is false in general** (D3af-1).
- **D45C4RED: the reduced flow reproduces the Euler runs, and the early selection moves away from basis pursuit.** ‖w‖₁ is 0.7793/0.7232 against Euler's 0.7782/0.7139, and mean ΔQ −0.338/−0.355 against −0.338/−0.356. The values now filled in agree with my rerun.
  - ‖w‖₁ rises monotonically on s ∈ [0, 270]: 0.706, 0.721, 0.742, 0.759, 0.769, 0.775, 0.777, 0.779.
  - The distance to basis pursuit also rises, 0.196 → 0.215. The early selection moves away from basis pursuit.
- **The mirror update's invariance holds by construction.** β is recomputed from Q, and η moves by h·g ∈ range Φᵀ. It is discrete mirror descent, i.e. exact zero-order-hold integration of each particle.
  - It is implementable in a real controller: closed-form η plus a few Newton steps per particle per sample. But `mirror_step` uses a fixed 8 iterations with no convergence check.
  - Its invariance is classical for stochastic mirror descent (memory).

### (iv)(c) Dither
- Correct; the drift matches the paradigm's dither √(2τG). "Q falls without bound" is shown in expectation only.

### (v) The real loop's Lyapunov function
- **The identity is correct.** Checked: the normalisation f̂ = Σw_iφ_i with w = βa/m, the sign convention e = x − x_m, and V ≥ 0, radially unbounded on the level set. Boundedness and e → 0 (Barbalat) follow.
- It is Boffi–Slotine's Bregman Lyapunov function with ψ = Ψ_Q.
- **The pair statement holds only with σ₀ = 0** (see (i)).
- **Part D's "sensor noise" reaches only the adaptation** — not the feedback, not the regressor φ(x). That is not a sensor noise. The pair's identity survives any common drive, so only the label and V's behaviour change.
- **The split is exact for the current sum.** The reported "split within 0.018" equals 0.668 × the sum's error (0.0272 → 0.0182).
- Parameter convergence is not shown.
- The free law loses D‴.34's a priori bound |f̂| ≤ a_max.

### (vi) The conserved law
- The identity is correct. "Not a mirror flow of a fixed potential" is not proved.
- Under noise, λ ≠ 0 at first order, so "the shift is the transient's" holds only noiselessly. E's numbers check out.

### Novelty and literature
- **(ii)–(iii) are Li–Wang–Lee–Arora's theorem.** In D‴.44's flat coordinates (α/√G, 2√β/√γ) the law is a Euclidean gradient flow of a separable 2→1 parametrisation, so the vector fields commute trivially. Marcotte–Gribonval–Peyré's count gives exactly one conserved function per unit, namely Q. Azulay et al. do the same for u⊙v with δ = u² − v² conserved.
- **The other headline items are classical too:**
  - (iv)(a)'s dual-span invariance — stochastic mirror descent: Gunasekar et al. 2018; Azizan–Hassibi 2019 (memory);
  - (iv)(b)'s mechanism — Pesme et al. and Even–Pesme–Gunasekar–Flammarion 2023, who give the per-step change of the effective initialisation in closed form (memory);
  - (v) — Boffi–Slotine;
  - the ℓ²→entropy shape — Ghai–Hazan–Singer's hypentropy.
- **Genuinely new:** the explicit potential of this saturating unit (the wall and the three regimes), the branches, the identification with D‴.44, and the mirror update as a non-selecting sampled-data law. The power-law rich weights 1/√β⁰ may also be new, but must be checked against depth-≥3 diagonal networks, whose rich-limit weights scale polynomially with the initialisation (memory).
- The AUDIT's own verdict is fair; the draft's "New:" list and ROW0 are not.
- **lit_S72 attribution.** It credits "isotropic noise does not" to Pesme et al.; the theorem is HaoChen–Wei–Lee–Ma (COLT 2021, search-verified).
- My searches found no prior closed form for the tanh-saturated, amplitude-times-mass unit, and no joint weight–amplitude Lyapunov function beyond Boffi–Slotine's general one.

### Numerics: transcription and interpretation
- **Every number in `numerics_D45.md` matches `mirror_free.md`.** The problems are interpretive:
  - silent solver failures inside medians;
  - a bimodal branch median;
  - checks that are true by construction;
  - instances chosen where the claim is automatically true.
- B's flow check is a real, independent test and converged for s ≥ 10⁻² (4.7·10⁻¹¹), but it never exercises the regimes where the solver fails.
- C(1) does exercise colour and node correlation (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes), but cannot discriminate Stratonovich from Itô.

### Tests and mutations
- **9 of 36 mutants pass all 7 tests:** M04 (raw formula everywhere), M09 (branch bracket ε → 10⁻²), M10 (clip 340 → 3), M11 (wrong side when out of range), M12 (Newton Jacobian without m), M13 (undamped Newton), M22 (mirror-step clip 0.05·a_max), M24 (mass law without a), M36 (stable form everywhere).
- **The reduced experiment misses M11, M12, M22, M29 and M36.** It sees M04, M09 and M10 only as hangs in B (timeouts).
- **Missed by both: M11, M12, M22, M36.** M12 and M36 do not change results. M11 and M22 sit exactly where I found real defects (clipping, unchecked Newton).
- **M24 passes the V̇ test** because `real_lyapunov` rebuilds β from (α, q) and never reads the state's β.
- **M13 passes** because the only projection test is at s = 1.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3af-1 | major | (iv)(b) and C(3): "over the selection's time the state reaches the rich limit (basis pursuit on the instances)" holds only because all six C(3) instances have basis pursuit = the 2-sparse truth. Where they differ, D‴.44's reduced flow (white noise, paradigm's update, lazy start, s to 10⁸) stops elsewhere, unchanged from s = 10⁶ to 10⁸: k=2 i=8 ends on the 5-atom fit {1,5,7,8,9} (‖w‖₁ 0.4228 vs basis pursuit's 0.4156, 0.099 from it); k=3 i=0 ends on {3,7,9,10}, 1.2·10⁻³ from basis pursuit. Atoms whose Q reaches −c are dead. No theory backs basis pursuit under this update (D‴.44(ii)'s D* is the flat scheme's). | State "reaches a sparse exact fit; basis pursuit's on 6/6 instances where it is the truth, not in general". Add basis-pursuit-≠-truth instances (k=2: i = 8, 9, 13, 17) to `c3_job` and report the end support (draft (iv)(b); numerics (3); exp `S_C3`/`c3_job`). |
| D3af-2 | major | Part B's medians include silent solver failures. `bregman_projection` uses the residual norm as merit function, accepts a failed line search with an inconsistent (z, w), returns no convergence flag, and inherits `alpha_of_eta`'s silent clip at \|α\| = 340. It diverges to \|w\| ≈ 2.2·10²⁹² in raw ('rich', 10⁻⁷) for k=2 i=15 and k=3 i=12, at s = 10⁻⁹ in all 7 instances tried, and on branch starts k=2 i=15 (residual 1.7·10³, reported distance 1.7·10³) and k=3 i=1 (residual 0.42, reported 0.464). My dual-descent solver converges in all these cases (k=3 i=12 at 8.3·10⁻⁷; k=2 i=15 at 6·10⁻¹¹), and from the two failed branch starts the loop reaches fits 0.318 and 0.026 from basis pursuit (residual 10⁻¹⁶). The 'rich'/'lazy' entries do not store residuals. | Run Newton on 𝒟 with Armijo on 𝒟, reject clipped α, return and record the residual and a flag, exclude and count failures in every aggregate (mirror_free `bregman_projection`; exp `b_job`, `report`). |
| D3af-3 | major | The branch result is misreported and misdescribed. (a) The k=2 median 6.6·10⁻⁴ is the midpoint of a bimodal set: 10 instances ≤ 1.2·10⁻⁴, one at 1.2·10⁻³, 8 at 0.05–0.37, plus the failure (truly 0.318). Half lose basis pursuit because a truth atom sits on a wrong-sign branch (k=2 i=0: 0.37; the loop confirms the projection to 4·10⁻¹²). (b) "Keep their signs in every instance" is true by construction of the solver and the geometry. (c) The "entropy vs ℓ¹" description is wrong (costs 340–480 vs 420–1220, λ₁ = 702). (d) As s → 0 an exact fit can cease to exist (2/20 and 1/20 instances at s = 10⁻⁶). | Rewrite (iii)'s branch sentences. Report the fraction with a wrong-sign truth atom and the full distribution. Drop "signs kept" (draft (iii); numerics (B); exp branch block). |
| D3af-4 | major | ROW0, BULLET and README say continuous-time noise of any colour or covariance "never moves that projection", and the "New:" list claims invariance "under every continuous-time noise". This is false for white noise read in Itô's sense, the standard continuous-time model of SGD noise (Pesme et al.'s stochastic gradient flow). The draft never says that their effect is (iv)(b)'s. | Say "noise with smooth paths (Stratonovich)" everywhere, and add: "Pesme et al.'s stochastic gradient flow is Itô; its drift correction is (iv)(b)'s ΔQ per unit time" (draft (iv)(a), Status; ROW0, BULLET, README, AUDIT). |
| D3af-5 | major | Novelty is overstated: (ii)–(iii) are an instance of Li–Wang–Lee–Arora (with Marcotte–Gribonval–Peyré and Azulay et al.); (iv)(a) is the classical invariance of stochastic mirror descent; (iv)(b)'s mechanism is Pesme et al.'s and Even et al.'s; (v) is Boffi–Slotine's; the shape is the hypentropy. See Novelty above. | Recast as "classical, made explicit for the paradigm's unit"; cite Ghai–Hazan–Singer 2020, Gunasekar et al. 2018, Azizan–Hassibi 2019, HaoChen et al. 2021, Even et al. 2023, Woodworth et al.'s depth-D section (Status; ROW0; lit_S72). |
| D3af-6 | major | (v)'s pair claim and Part D lie outside the stated setting ("one particle per class beyond the cohesion's range"). In the paradigm the cohesion acts on α_i − α_j within a class, which breaks Q and the equal mirror increments, so the split is no longer the Bregman split. | Assume σ₀ = 0, or no cohesion on α, for the pair, or analyse the pair with cohesion (draft (v); exp `part_D` docstring). |
| D3af-7 | minor | (i)'s "also when centres and widths move by any law": moving centres can enter the cohesion's range, and the Gaussian ψ never vanishes. | Add "while the cohesion's α-component vanishes" (draft (i); BULLET "Open"). |
| D3af-8 | minor | "Lowers the masses' scale … not monotonically in ‖w‖₁" understates the early phase. On s ∈ [0, 270] ‖w‖₁ rises monotonically (0.706 → 0.779) and the distance to basis pursuit rises (0.196 → 0.215); C(3) samples only s ≥ 10⁵ and shows no non-monotonicity itself. | Say "first away from basis pursuit (+10% in ‖w‖₁ at s = 270), turning after s ~ 10⁴–10⁵". Add s = 10², 10³, 10⁴ to `S_C3` (draft (iv)(b); numerics (3)–(4)). |
| D3af-9 | minor | The rich remainder O(\|w\|(1+\|log\|w\|\|)) is false; it is −(3/2)K\|w\|^{2/3} (numbers in (ii)). The Γ-limit survives, and the correct form also explains the rate (inactive atoms carry \|w\| ~ β⁰^{3/2}). | Replace with O(\|w\|^{2/3} + \|w\|·log₊\|w\|) in (ii) and in the proof of (iii). |
| D3af-10 | minor | The rich limit needs uniqueness, which only the proof assumes; duplicated atoms (Part D) break it. The limit's scale depends on the LP's vertex gap (k=3 i=13 at 0.17 until s ~ 10⁻¹¹), and B's medians hide 3/20 instances ≥ 3·10⁻³ at s = 10⁻⁷ plus one divergence. | Add "when unique (otherwise the minimiser selected by the next order)". Report maxima and quantiles (draft (iii); numerics (B)). |
| D3af-11 | minor | The time to converge grows like 1/β⁰²: 10⁶, 10⁷ and 10⁹ at s = 10⁻², 10⁻³, 10⁻⁴. "No rate condition" and "starting masses as a sparsity prior" omit this cost. | State the rate; integrate one s = 10⁻³ instance (draft (iii); AUDIT; exp `FLOW_MIN_B`). |
| D3af-12 | minor | Code defects. `eta`: NaN at q = −c; docstring wrong on the negative branch (η → +∞ there); not cancellation-free for \|t\| ≥ 0.9, and with the raw form everywhere (M04) the B lazy projection fails. `alpha_of_eta`: silently returns ±340 (for η = 10⁴ at q = 0.5 it gives η(340) = 6765); on a branch, brentq raises ValueError (η = −10⁵ at q = −2c). `mirror_step`: no convergence check. | Use 1 − t = 2/(1 + e^{2\|α\|/a_max}) with log1p for all t; handle q = −c; raise on clipping; guard the branch bracket; check convergence (mirror_free.py). |
| D3af-13 | minor | Checks that cannot fail: "signs kept"; the mirror update's Q to 1.1·10⁻¹⁶ (β is recomputed from Q); Part D's split (0.668 × the sum's error); A's shift check (recomputes the rhs's own λ); C(1) (an ODE input is Stratonovich by construction). | Replace with discriminating checks: Euler–Maruyama against the mirror update at fixed intensity v as h → 0, and Part D's split against the Bregman split of the current sum (numerics; exp). |
| D3af-14 | minor | "Never selects … at any order" holds for the limit map. Under persistent noise the stationary mean is an exact fit displaced at O(v) in a covariance-dependent direction. | Add one clause (draft (iv)(a)). |
| D3af-15 | minor | The Status overclaims: "(iv) proved" though (iv)(b) contains numerical claims that are false in general (D3af-1); "(vi) proved" though "not a mirror flow" is unproved; (iv)(c) holds in expectation only. | Correct the Status. |
| D3af-16 | minor | Part D's noise does not reach the plant or φ, so it is not a sensor noise (D‴.41's loop puts it in both laws with the desired regressor). "Parameters converge slowly" is not shown. | Rename it, or add n to the feedback and to φ(x + n); show a trend (numerics (D); `part_D`). |
| D3af-17 | minor | lit_S72 credits "isotropic noise does not" to Pesme et al.; precedents are missing (see D3af-5). | Fix lit_S72. |
| D3af-18 | minor | No cross-reference to D‴.41(ii): its C ∝ I requirement is the conserved law's, from the centring in (vi). | One sentence in (iv)(a). |
| D3af-19 | minor | The free law loses D‴.34's bound \|f̂\| ≤ a_max. Under approximation error or disturbance the masses can drift (D‴.46 is the remedy). | One sentence in (v). |
| D3af-20 | minor | The tests miss 9/36 mutants (M04, M09–M13, M22, M24, M36). The V̇ test is blind to the mass law; there is no rich-regime, branch, or `alpha_of_eta` range test. | Make `real_lyapunov` use the state's β or assert Q; add projection tests at s = 10⁻⁴ and 10⁻⁷ with a branch; add η/α round trips at the extremes (tests). |
| D3af-21 | cosmetic | "Median distance … falling 0.14, 0.16" (it rises); "truth recovered as often as basis pursuit" (0.75 vs 0.80); "Σβ = m kept to 10⁻⁶" is print precision; C(2)'s 4·10⁵ Monte Carlo draws could be a 1-D quadrature. | numerics_D45.md; `part_C12`. |
| D3af-22 | cosmetic | The proof of (ii) drops the factor (c − Q) in ∫2a_max dt/(2(…)); (iii) uses β⁰_i for the start mass while (ii) defines β⁰ := Q + c; ROW0's "soft walls for large starting masses" (the walls do not matter in the lazy limit). | Draft; pieces. |

## Lesson
The exact identities were right and cheap to verify. What failed was everything around them that was allowed to pass unexamined:
- a numerical claim generalised from instances on which it holds automatically (basis pursuit = 2-sparse truth);
- medians over solver outputs whose convergence was never recorded, hiding values of 10²⁹² and bimodal splits;
- checks that are true by construction (signs kept, Q under the mirror update, the pair's split);
- summary pieces that drop the one qualifier (Stratonovich) the theorem depends on.

A "numerical" claim must be tested where it can fail, every aggregate must carry its failure count, and every shortened statement must keep its hypotheses.

Sources: [Pesme et al., NeurIPS 2021](https://proceedings.neurips.cc/paper/2021/hash/f4661398cb1a3abd3ffe58600bf11322-Abstract.html); [HaoChen et al., COLT 2021](https://proceedings.mlr.press/v134/haochen21a.html); [Ghai–Hazan–Singer, ALT 2020](https://proceedings.mlr.press/v117/ghai20a.html); [Boffi–Slotine, arXiv 1912.13154](https://arxiv.org/abs/1912.13154v4); [Even et al., NeurIPS 2023](https://papers.neurips.cc/paper_files/paper/2023/hash/5da6ce80e97671b70c01a2e703b868b3-Abstract-Conference.html); [Bai et al., arXiv 2607.12332](https://arxiv.org/pdf/2607.12332).
