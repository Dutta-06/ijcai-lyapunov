# Referee report D3w (round 24) — PARADIGM §5d: Proposition D‴.34 (mass transfer: the replicator weight law; statement (i)–(v), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.34 |", the §8 bullet "Done (Proposition D‴.34)", `results/RESULTS_PARADIGM.md` (Proposition D‴.34), LITERATURE_TRANSPORT [S62] and the audit bullet "Proposition D‴.34", PAPERS.md (paper 13, as far as it concerns D‴.34), README.md (the replicator), `src/core/mass_transfer.py`, `experiments/exp_mass_transfer.py` → `results/chaos/mass_transfer.md` (+ raw), `tests/test_mass_transfer.py`; secondary: Proposition D‴.32 as revised in round 23 ((v), Status, *Numerics* (I)–(L), Parts I–L of `exp_cohesion_leak.py`)

Branch `claude/festive-bardeen-eopej2` at 1ddff8c, HEAD when the review started (`git log --oneline -5`: 1ddff8c, a6249e8, 64cb25e,
a9a0e5a, bfdbc0f).

The tree was extracted with `git archive 1ddff8c` into `ref_D3w/tree`, and everything below was run there. Uncommitted edits in the
shared tree were ignored. No repository file was modified and no git command touching the index or the tree was run. Python ran with
`PYTHONDONTWRITEBYTECODE=1` and one thread (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`), at most two processes at a time; pytest with
`-p no:cacheprovider`. The committed `mass_transfer.md` equals the scratch tree's (md5 `3cd2332c9747569b7170ca7d3adc2700`); `--from-raw` on the scratch tree rewrites it byte-identically (`ref_D3w/md5_before.out`, `md5_after.out`). The four tests pass
(`ref_D3w/tests.out`: 4 passed, 119 s).

Every number below comes from a script of mine in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3w/scripts/`
(each with its `.out`; the runs' outputs are `B_five_1.out`, `F_heavy_0.out`, `partC_0.01_0.out`, `cont_five_gauss.out`,
`cont_five_compact.out`, `gauss_floor_five.out`, `gauss_floor_five_cont.out`, `gauss_floor_3eq.out`). `mi.py` is an independent implementation that imports nothing from the repository: the saturated scalar
network (a = a_max tanh(α/a_max), c = c_max tanh(γ/c_max), s = s_min + (s_max − s_min)σ(κ)), the quadrature on x = 2 sin t (256 points),
the weighted network, its residual, the analytic gradients in the unconstrained coordinates, the fitness u_i, both cohesion kernels and
their energy (my own formulas), the weighted drift, J, the midpoint rule on (Θ, log β) with and without renormalisation, single linkage
(scipy) and the broad-start draw in the experiment's order. The repository's objects are used only where marked as a cross-check
(`check_i.py`), and the repository's tests are run against mutants of its code (`mutants.py`).

| script | what it does |
|---|---|
| `mi.py` (module) | the independent instance above |
| `check_i.py` | (i): the drift against −(metric)⁻¹∇J by central differences of my J (both kernels, random β ∈ [0.3, 2]), the Fisher–Rao projected gradient also by a constrained minimisation (SLSQP), dJ/ds against the metric norm, the √β form, the WFR drift, |f̂| ≤ a_max on 2000 random draws; cross-check of my drift and J against the repository's |
| `check_ii.py` | (ii): the real loop at the instance's atoms (my RHS, DOP853 at 10⁻¹²): V, V̇ from the vector field, Σβ, β(t) → β*; a start off Σβ = m; the identity for n = 3 at random points |
| `check_A.py` | (iv), Part A: the frozen-position runs (γ = 0, 10⁻²; amplitude scales 1 and 1.25, 4·10⁴ slow units), the Gram matrix, the box-constrained infimum (γ = 0) and the ℓ¹-ball infimum (any masses) by L-BFGS-B and SLSQP |
| `check_v.py` | (v), Part E: the Jacobian in (Θ, log β) at the lift of ϑ^e(w), counts (6, 5, 4), compact R = 1, γ = 0, 10⁻², 1; its null space against the predicted redistributions and gauges; the null space on Σδβ = 0; the K-th gauge as an equilibrium off the constraint; the separations and the collective eigenvalues; J's Hessian |
| `broad.py` | broad starts (m = 60, `default_rng(seed)`, the experiment's draw order), the annealed Gaussian schedule R = 2 → 1 (2·10⁴ + 2·10⁴), any γ; saves the end state |
| `gauss_floor.py` | the Gaussian kernel's stationary point at the class level (replicated classes; the positions' balance with the factor β_i and, optionally, u_k = λ), by least squares from a run's end classes; also three equal bumps (D‴.32 (L)) |
| `cont.py` | a saved end state continued at R = 1 (either kernel) |
| `partC.py` | Part C: ϑ^b(q) by my own BFGS continuation in ε, the jittered lift (`default_rng(700 + seed)`), compact R = 1, σ₀/64, 10⁵ slow units |
| `raw_read.py` | the raw record of `mass_transfer_raw.npz` (Part B's last records, Part C's cluster counts along the runs, Parts D and F) |
| `mutants.py` | the repository's four tests against nine in-memory mutants of `mass_transfer.py` |
| `diag_five.py`, `diag_five2.py` | the five-bump end states (after the Gaussian and the compact continuations): class spreads, the collapsed fit, the particle Jacobian (240 × 240), the class-level exact fit nearby, the class-level flow for 2·10⁵ slow units |
| `renorm.py` | Part A (γ = 10⁻²) with and without the renormalisation of Σβ |
| `calc32.py` | D‴.32 (v): the twins' rate at random parameters, the time-scale integral, the dead class's time scale; Part L's decay rates |

Line numbers are 1ddff8c's:
- PARADIGM: the statement l. 4262–4269, (i) l. 4271–4276, (ii) l. 4278–4280, (iii) l. 4282–4289, (iv) l. 4291–4295, (v) l. 4297–4302;
  the proof l. 4304–4317; the Status l. 4319–4327; the *Numerics* l. 4329–4352. The §0 row l. 144, the §7 row l. 6792, the §8 bullet
  l. 7056–7063. D‴.32: (v) l. 3984–3996, its proof l. 4028–4038, the Status l. 4045–4060, *Numerics* (A4) and (I) l. 4074–4086, (J)
  l. 4099–4102, (K), (L) l. 4138–4148.
- RESULTS l. 890–911. LITERATURE_TRANSPORT [S62] l. 82, the audit bullet l. 519–531. PAPERS.md l. 225–243. README l. 71–74.
- Code: `mass_transfer.py` l. 1–23 (docstring), 29–86 (`ReplicatorCohesionLoop`), 89–96 (`class_masses`), 99–116 (the real loop and V);
  `exp_mass_transfer.py` l. 76–90 (A), 94–113 (B, D, F), 117–129 (C), 133–174 (E); `tests/test_mass_transfer.py` l. 39–124.

## Verdicts

| object | verdict |
|---|---|
| **(i) the gradient flow; the bound** | **holds.** With J = mL̄(μ_{Θ,β}) + (σ₀m/2)∫∫W dμ_Θdμ_Θ (my own J and drift) the positions' drift is −G∂_ΘJ and the weights' drift the Fisher–Rao steepest descent on Σδβ = 0, to 9.4·10⁻¹⁰ and 8.3·10⁻¹⁰ (relative, central differences; both kernels, β ∈ [0.3, 2]); the constrained minimisation gives the same weight velocity to 2.9·10⁻⁸ (SLSQP's tolerance); dJ/ds equals minus the metric norm (−0.85950 against −0.85950); Σdβ/ds = −1.3·10⁻¹⁶. With r_i = √β_i the weight law is (γ/4) times the projected Euclidean gradient in r on Σr_i² = m (to 6·10⁻¹¹). The WFR remark is right for the fit term but not complete (ID-4). |f̂| ≤ a_max: max 3.75 < 4 over 2000 random draws. β > 0 and Σβ = m invariant; the run's renormalisation is harmless (Part A without it: Σβ − m = −8.1·10⁻⁹ after 4·10⁴ slow units, β within 3.7·10⁻⁹, the fit error 4.5·10⁻¹⁴ against 3.0·10⁻¹⁴). My drift and J equal the repository's to 2.8·10⁻¹⁷ (cross-check). |
| **(ii) the relative entropy** | **holds** (V̇ = −k_e|e|², re-derived and checked to 1.1·10⁻¹⁶ of max k_ee² along the run). Its hypotheses are partly unneeded and partly unstated (ID-5): V̇ needs only Σβ* = m (not n = 1, not Σβ = m); V ≥ 0 (the Bregman reading) needs Σβ = m, which is invariant in real time but not attracting. β → β* is not claimed; it holds numerically (1.3·10⁻³ at t = 200) and would need persistent excitation. |
| **(iii) stationarity; the transient** | **holds only for the finite range at classes R apart.** The "iff … the transport condition of Proposition D‴.32" misses the factor β_i, and under the Gaussian kernel — the kernel of every broad-start run — the stationary set is not D‴.32's: the masses are fixed by u_k = λ ≠ 0 and the bias changes (ID-1). The transient's rate dW_k/ds = γW_k(1 − W_k)u_k is exact for a coincident class with u_l = 0 elsewhere; "biased for times exponential in α_k" is a heuristic (ID-8). |
| **(iv) representability at fixed positions** | **holds as a dichotomy for critical points**: they exist iff Σ_k|A_k|/a_max < 1 (strict; the boundary case has none), and Łojasiewicz applies (L̄ and the metric are analytic in (log W, α)). Hidden hypotheses (f in the span of the φ_k, one particle per atom), a gap in "leaves every compact set" (true by the ω-limit argument, not by the one given) and summaries that drop "bounded" (ID-3). Part A reproduces; the scale-1.25 run creeps to the ℓ¹-constrained infimum 0.0684, which the text does not give (ID-3). |
| **(v) the zero modes** | **holds exactly.** In (Θ, log β) at the lift (6, 5, 4): 15 zero eigenvalues at γ = 0, 10⁻², 1, the null space spanned by my 12 redistributions and 3 gauges (residual 3·10⁻⁹, J·P = 5.7·10⁻¹¹); 14 on Σδβ = 0; the K-th gauge is an exact equilibrium curve off the constraint (drift ≤ 8.5·10⁻¹⁸ at β ≡ 0.9 and 1.1 with the amplitudes rescaled); J's Hessian PSD (min −4.4·10⁻¹¹, the differences' error; 15 zeros); the 36 separations at −σ₀w_k unchanged by γ; the 9 collective eigenvalues change with γ (the slowest −3.0798·10⁻⁴ → −3.0812·10⁻⁴ → −3.0844·10⁻⁴, the fastest −0.1834 → −0.1834 → −0.6108 at γ = 0, 10⁻², 1). The real loop's "manifold of periodic orbits" is right (e ≡ 0 on exact representations) but not checked numerically. |
| Proof | (i), (ii), (v) **correct**; (iii) **correct for the finite range**, the first sentence wrong as stated (ID-1); (iv) **correct for critical points and bounded trajectories**, the unbounded case argued too briefly (ID-3). |
| *Status* | accurate, except that "(iii)'s stationarity … proved" is read in the §0 row, §8, PAPERS, README and the audit as "the stationary set unchanged" for the runs' kernel (ID-1). |
| *Numerics* | **the printed numbers reproduce**: A (all four runs I made, to every printed digit), B (γ = 10⁻², seed 1), C (γ = 10⁻², seed 0: 8 classes, 0.00548, the cluster counts along the run equal to the raw record), F (three heavy, seed 0: 1.04·10⁻⁷, masses 0.411, 0.309, 0.280), E (the real loop: V 0.0418 → 2.14·10⁻⁷, 1.11·10⁻¹⁶, β at the end 1.00035, 0.99873, 1.00092; the linearisation to all printed digits). Readings overstated: (B)'s "still falling … 3–4% per 1000 slow units" (ID-2); (C)'s "the same splits at the same times" holds except at three records (ID-8). |
| §0 row, §7 row, §8 bullet, RESULTS, PAPERS, README | numbers hold; "the stationary set unchanged / keeps the stationary set / its stationary points are D‴.32's" (ID-1); "at fixed positions the exact fit is reached iff the target is representable" without "bounded" (ID-3); "still falling" (ID-2). |
| LITERATURE [S62], audit | bibliographic details right where checked; the verdict fair; precedents missing (ID-7). |
| `mass_transfer.py` | **holds** (my independent drift, J and real-loop law agree; the module's λ = Σβu/Σβ differs from the text's (1/m)Σβu off the constraint only, ID-5). |
| `exp_mass_transfer.py`, md, raw | **holds**; byte-identical regeneration; every reproduced run matches the raw record. |
| `tests/test_mass_transfer.py` | **pass (4/4)**; the drift mutants are caught only by the gradient-identity test, which compares the code with its own J; the representability and zero-mode tests pass every drift mutant but the sign flip (ID-6). |
| D‴.32's round-23 revisions (D3u-1–D3u-10) | **applied correctly** where checked: the twins' rate re-derived (and checked at random parameters to 4·10⁻¹⁶), the time scales, Parts I–L's readings; one cosmetic slip in (A4) (ID-9). |

## Confirmed issues

**ID-1 (major; (iii) l. 4282–4285 "(Θ, β) is stationary iff the positions satisfy the transport condition of Proposition D‴.32 and u_i = λ
on the support … the stationary points are D‴.32's, for every vector of masses"; the §0 row l. 144 (statement column "the stationary points
are Proposition D‴.32's for every vector of masses", novelty column "the stationary set unchanged"); the Status l. 4326 "its stationary set
unchanged"; the §8 bullet l. 7059 "Its stationary points are D‴.32's"; PAPERS l. 226 "keeps the stationary set"; README l. 72–73 "without
changing the stationary set"; the audit l. 528 "the stationary set unchanged").** True for the finite range at classes at least R apart
(the case the proof treats); false as stated in general and false for the Gaussian kernel, which every broad-start run (B, D, F) uses.
- *The general condition.* The weighted loop's positions are stationary iff β_iG⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ − (σ₀/m)Σ_jψ_ij(ϑ_i − ϑ_j) = 0: the
  fit term carries β_i, the cohesion term does not (the cohesion is the unweighted cloud's). This is not D‴.32's transport condition
  unless the cohesion force vanishes on the configuration (coincident classes R apart under the finite range), where it reduces to
  ⟨r∇Φ̃_k⟩ = 0 and β drops out. Inside a class, particles with different β_i feel different fit forces and the same cohesion, so at a
  biased stationary point a class with unequal weights is not even replicated. In the runs the weights inside a class differ by up to
  2.2× (seed 1 of (B): max/min 1.30, 2.01, 1.45, 1.31, 2.21 per class). At the end of 10⁵ further slow units at Gaussian R = 1 the
  classes' particles lie up to 2.5·10⁻³–6.8·10⁻³ (whitened) from their β-weighted atoms, against 8·10⁻⁶–4.4·10⁻⁴ under the compact kernel
  over 4·10⁴ (`diag_five.out`).
- *The Gaussian kernel* (`gauss_floor.py`, class level: replicated classes, equal weights inside each, the factor β_i, least squares to a
  residual of 10⁻¹⁵). Five bumps, from the end of (B) seed 1, γ = 10⁻², continued 10⁵ slow units at Gaussian R = 1 (the counts 6, 6, 18, 20, 10;
  `gauss_floor_five_cont.out`; from the run's own end, `gauss_floor_five.out`: 3.85·10⁻⁴ at its masses, 2.63·10⁻⁴ with u_k = λ):
  - with the masses held at the run's (0.1581, 0.1322, 0.2408, 0.2623, 0.2066), the stationary point has the fit error 2.44·10⁻⁴;
  - with the replicator's condition u_k = λ added (the masses free), the masses move to (0.1581, 0.1320, 0.2365, 0.2660, 0.2075), the fit
    error is 2.63·10⁻⁴, and all three positive classes sit at |a|/a_max = 0.8468, both negative ones at 0.9567;
  - under the compact kernel R = 1 the same start gives the exact fit (3.0·10⁻¹⁶).
  So under the Gaussian kernel the replicator selects the masses (one stationary point near each class configuration, not one per mass
  vector) and moves the bias. The same holds on D‴.32 (L)'s three equal bumps (counts 16, 16, 28): without the replicator the
  stationary fit error is 9.387·10⁻⁵ (the md's 9.39·10⁻⁵), and with it the masses go from (0.267, 0.267, 0.467) to
  (0.2205, 0.2205, 0.559), the fit error 9.17·10⁻⁵.
- *Why it matters.* The summaries make the stationary set's invariance the proposition's selling point next to the broad-start
  numerics, which run under the Gaussian kernel, where it does not hold.

*Fix:* in (iii) write the stationarity condition with the factor β_i; restrict "the stationary points are D‴.32's for every vector of
masses" to the finite range at classes R apart; add that under the Gaussian kernel the stationary points are fixed by u_k = λ too (the
masses are no longer free), with a bias of the same order as D‴.32 (ii)'s (the five-bump numbers above). Qualify the §0 row, the Status,
the §8 bullet, PAPERS, README and the audit the same way ("at clustered configurations under the finite range").

**ID-2 (minor; *Numerics* (B) l. 4337–4338 "the fit error 6.0·10⁻⁴–1.4·10⁻³ and still falling at the end (by 3–4% per 1000 slow units; the
classes at |a|/a_max up to 0.96)"; RESULTS l. 897 "still falling"; the §0 row's numerics column "6.0·10⁻⁴–1.4·10⁻³, still falling").** The
rate is not 3–4% in every run, and it does not persist: the fit levels off near the Gaussian kernel's stationary point, about half the
printed value, and below that a saturation-slowed collective mode takes over.
- *The rate.* The raw record's last five ratios per 1000 slow units are 1.040–1.042 (seed 0), 1.035–1.039 (seed 1) and 1.019–1.021
  (seed 2) (`raw_read.out`), so 2–4%.
- *The continuation* (`cont.py`, from my reproduction of seed 1, which matches the md: 5.99·10⁻⁴, centres −1.600, −0.798, +0.004, +0.798,
  +1.597). At Gaussian R = 1, γ = 10⁻², the fit error is 5.99, 4.45, 3.62, 3.20, 3.01, 2.91, 2.87, 2.84, 2.82, 2.81, 2.80·10⁻⁴ at
  s = 0, 10⁴, …, 10⁵. The last 10⁴ slow units take off 0.4%. The level is the Gaussian stationary point of ID-1 (2.63·10⁻⁴ at the class
  level; the run's unequal weights inside the classes shift it).
- *The compact kernel.* Under the compact kernel R = 1 the same end state reaches 2.85·10⁻⁴ after 4·10⁴ slow units (5.99, 4.46, 3.59, 3.11,
  2.85·10⁻⁴ at s = 0, 10⁴, …, 4·10⁴), although an exact five-class representation with these counts lies within 0.05–0.09 (whitened) of its classes (class-level
  fit error 3.0·10⁻¹⁶). The drift there is 3·10⁻⁷ and the class-level flow from it falls only to 1.67·10⁻⁴ in 2·10⁵ slow units (2.85 → 2.60·10⁻⁴ over the first 2·10⁴,
  1.76 → 1.67·10⁻⁴ over the last; `diag_five2.out`). The two negative classes sit at |a|/a_max = 0.956 (a′ = 0.086), so the amplitude directions are slow.
  The particle Jacobian there has small positive eigenvalues (up to 1.3·10⁻⁷, `diag_five.out`), which I did not analyse.
- *What limits the fit.* Not time alone: the Gaussian bias at R = 1 (≈ 2.6·10⁻⁴) and, under the compact kernel, slow modes of the nearly
  saturated classes.

*Fix:* replace "still falling at the end (by 3–4% per 1000 slow units)" by the rates per seed (2–4%) and say that the fit levels off at the
Gaussian kernel's stationary point (≈ 2.6·10⁻⁴–2.8·10⁻⁴ for seed 1), and that under the compact kernel it is still slower; report a
continuation if the text keeps "still falling".

**ID-3 (minor; (iv) l. 4291–4295 and the proof l. 4308–4313; the §0 row's statement column "at fixed positions the exact fit is reached iff
Σ_k|A_k|/a_max < 1"; the §8 bullet l. 7060; PAPERS l. 227; README l. 73; RESULTS l. 893–895; *Numerics* (A) l. 4334–4335).**
- *The "if" direction is not proved.* The proposition proves that bounded trajectories converge to a critical point (the Status:
  "boundedness numerical"), but the §0 row's statement column, the §8 bullet, PAPERS and README say "the exact fit is reached iff the
  target is representable". With the masses free a trajectory can drift along the level sets toward W_k → |A_k|/a_max with α_k → ∞ (the
  closure of the exact-fit set). Nothing shown excludes this.
- *Hidden hypotheses.*
  - "the exact fits p = A" presumes f = Σ_kA_kφ_k. In general the critical points are p = p*, the least-squares coefficients, and the
    criterion is Σ_k|p*_k|/a_max < 1.
  - The (W, α) reduction needs one particle per atom (with several, their amplitudes adapt separately).
- *The unbounded case.* The proof's "if there is none it cannot remain bounded" shows only that the trajectory is not bounded. "Every
  trajectory leaves every compact set" needs the standard argument: the dissipation ∫|∇J|² < ∞ makes every ω-limit point critical, so an
  empty critical set means the trajectory eventually leaves each compact set. The argument holds, but it is not the one given.
- *"(the classes saturate)".* At Σ_k|A_k|/a_max ≥ 1 the infimum is the least-squares fit over the ℓ¹ ball Σ_k|p_k| ≤ a_max. There every
  class has |p_k| = W_ka_max, i.e. it saturates or its mass tends to 0 (p_k = 0). The latter is not excluded.
- *Part A's 0.071* (`check_A.out`). The ℓ¹-ball infimum is 0.0684, at the masses (0.270, 0.164, 0.194, 0.132, 0.239). The run's end masses
  already equal these to 10⁻³ (0.271, 0.164, 0.194, 0.132, 0.239), and the fit "creeps down" toward 0.0684 as the amplitudes saturate.
  The box infimum at the fixed masses (γ = 0) is 0.3874 (scale 1, the md's 0.388) and 0.4942 (scale 1.25, the md's 0.494).

*Fix:* write "bounded trajectories reach the exact fit iff …" in the summaries (or prove boundedness); state f ∈ span{φ_k} (or the
criterion with p*) and one particle per atom; give the ω-limit argument; write "(each class saturates or loses its mass)"; print the
infimum 0.068 beside 0.071.

**ID-4 (cosmetic; (i) l. 4273–4276).**
- *The WFR remark.* "(The Wasserstein–Fisher–Rao flow of μ_{Θ,β} differs in the metric on the positions … its positions' drift lacks
  the factor β_i.)" is right for the fit term. In the WFR metric Σβ_i‖δϑ_i‖² the cohesion force is divided by β_i too (checked to
  9·10⁻¹⁰, `check_i.out`). And J is not a functional of μ_{Θ,β}, because its cohesion is the unweighted cloud's, so "the WFR flow of
  μ_{Θ,β}" is defined for the fit part only.
- *The √β form.* "gradient descent on (r, Θ) … the output weights of a two-layer network" holds with the step γ/4 in r (6·10⁻¹¹). But
  the output weights are r_i²/m, not r_i, and r_i clashes with the residual r.

*Fix:* "differs in the metric on the positions (its drift lacks the factor β_i on the fit and divides the cohesion force by β_i; J is not a
function of μ_{Θ,β} alone)"; rename r_i (e.g. ρ_i := √β_i), "the square roots of the output weights", and give the step γ/4.

**ID-5 (cosmetic; (ii) l. 4278–4280, the proof l. 4306–4308, (𝔉) l. 4266; the brief's "n = 1").** The hypotheses of (ii) are not the
ones needed (`check_ii.out`).
- *What V̇ needs.* d/dt Σβ*_i log(β*_i/β_i) = −γeᵀ(mf − mf̂) needs only Σβ* = m and the representation f = (1/m)Σβ*_iΦ̃_i at every state
  the loop visits. Neither n = 1 nor Σβ = m is needed: V̇ + k_e|e|² = 0 at random points for n = 3, and V̇ = −k_ee² to 1.4·10⁻¹⁷ from
  Σβ(0) = 3.3. For n = 1 and analytic features, the representation on the orbit's interval implies it everywhere. In ℝⁿ the averaged
  exact representation (on supp ρ) does not suffice.
- *What V ≥ 0 needs.* The Bregman reading needs Σβ = Σβ*. In real time Σβ = m is invariant but not attracting:
  d(Σβ)/dt = γeᵀf̂(m − Σβ), and from Σβ(0) = 3.3 it ends at 4.70 after 200 time units, where Σβ*log(β*/β) can be negative.
- *The paper and the code.* (𝔉)'s λ := (1/m)Σβ_ju_j and the module's λ = Σβ_ju_j/Σβ_j agree on the constraint only.
- *β → β*.* It is not claimed and needs persistent excitation. Numerically |β − β*|_∞ = 0.058, 0.015, 0.0043, 0.0013 at t = 50, 100,
  150, 200.

*Fix:* state (ii) with Σβ* = m, Σβ(0) = m and the pointwise representation; drop n = 1 or say it is the instance's; note that the
β-normalised λ would conserve Σβ for every start; add "β → β* under persistent excitation of (Φ̃_i(x_m(t)))_i".

**ID-6 (minor; `tests/test_mass_transfer.py` l. 39–124; RESULTS l. 905–911; the §7 row's last column).** Nine mutants of
`mass_transfer.py`, patched in memory (`mutants.py`, `mutants1.out`, `mutants2.out`):

| mutant | gradient identity | real loop V | fixed positions | zero modes |
|---|---|---|---|---|
| λ = 0 (unbalanced) | FAIL | pass | pass | pass |
| λ = the unweighted mean of u | FAIL | pass | pass | pass |
| the factor β_i on the positions dropped | FAIL | pass | pass | pass |
| the factor β_i on the weights dropped (additive law) | FAIL | pass | pass | pass |
| the fitness's sign flipped | FAIL | pass | FAIL | pass |
| real loop: additive weight law | pass | FAIL | pass | pass |
| real loop: λ dropped | pass | FAIL | pass | pass |
| the run without renormalisation | pass | pass | pass | pass (harmless, see (i)) |

- Every drift mutant is caught by one test only, the gradient identity. It compares the code's drift with the code's own J, so a
  consistent error in both (a wrong energy) passes. Here I checked the identity independently against my J.
- The representability test passes with λ = 0, i.e. with Σβ not conserved. It never asserts Σβ = m, and its ">" asserts are loose:
  fits > 0.05 at Σ|A_k|/a_max = 1.094, whose infimum is 0.068.
- The zero-mode test passes every mutant. Its γ = 0 count (7) is trivial, since the log β rows vanish. At γ > 0 it checks the number of
  zeros, not that they are the redistributions and gauges, and of the rest only that it is negative, not that the separations stay at
  −σ₀w_k.

*Fix:* assert Σβ = m along the fixed-position runs; check the null space against the predicted vectors (as in `check_v.py`) and the
36 (here 4) separation eigenvalues at γ > 0; check the drift against a closed form at one configuration (e.g. the class-level formula at
a lift), not only against J.

**ID-7 (minor; LITERATURE_TRANSPORT [S62] l. 82, the audit l. 519–531).** The audit's verdict is fair (classical ingredients; the
combination new), and the cited details are right: Liero–Mielke–Savaré, Invent. Math. 211(3):969–1117, 2018; Chizat, Math. Program.
194:487–532, 2022; Rotskoff–Jelassi–Bruna–Vanden-Eijnden, ICML 2019, PMLR 97:5508–5517, the PMLR title being the one cited (web
search); Boffi–Slotine, Neural Comput. 33(3):590–673, 2021, doi:10.1162/neco_a_01360 (web search); Kivinen–Warmuth and Hofbauer–Sigmund
from memory. Missing, from memory (not checked):
- *The other two origins of the WFR distance:* Kondratyev–Monsaingeon–Vorotnikov (Adv. Differential Equations, 2016) and
  Chizat–Peyré–Schmitzer–Vialard (Found. Comput. Math., 2018).
- *The lift to (√β, Θ):* Chizat–Bach (NeurIPS 2018), gradient descent on the positions and the square-root weights of a two-layer
  network as a WFR-type flow. It is the closest precedent of "(i)'s gradient descent on (√β, Θ)".
- *The metric:* Shahshahani (Mem. Amer. Math. Soc., 1979).
- *Mixture-weight learning:* EM's multiplicative update of mixture weights.
- *The real-time law:* (ii) is the negative-entropy instance, on the simplex, of the mirror-descent adaptive laws with Bregman Lyapunov
  functions that the audit already credits to Boffi–Slotine, so (ii) should be listed as classical, not only "the exponentiated-gradient
  law".

*Fix:* add these to [S62] and name (ii) as an instance of the mirror-descent adaptive law.

**ID-8 (cosmetic; (iii) l. 4283–4287; *Numerics* (C) l. 4340–4341).**
- *"whose amplitudes are not saturated".* The phrase is vacuous: |a| < a_max always holds.
- *"so the fit stays biased for times exponential in α_k".* At fixed positions and equal weights a class short of mass stays biased
  forever (no critical point). What is exponential is the time α_k needs to grow (α ≈ (a_max/2) ln s); the Status calls it numerical,
  and the text should say heuristic.
- *The rate.* dW_k/ds = γW_k(1 − W_k)u_k is exact when the class's particles coincide and u_l = 0 for every other class.
- *(C)'s "from the same splits at the same times".* The raw cluster counts with and without the replicator differ at three of
  51 records (seed 0 at s = 6.8·10⁴: 5 against 6; seed 1 at 7.0·10⁴: 5 against 6 and 8.0·10⁴: 6 against 7). Seed 2 is identical. My
  reproduction of γ = 10⁻², seed 0 has the raw record's counts at every record.

*Fix:* drop "whose amplitudes are not saturated"; "the fit stays biased (at fixed positions, forever)"; "the same splits at nearly the
same times".

**ID-9 (cosmetic; D‴.32 *Numerics* (A4) l. 4076–4077 "A dead class (amplitude 0, at 5.19 from the nearest live class, of mass 4/15)";
D‴.34 (D) l. 4339).**
- *The dead class's mass.* It is 2/15 (`exp_cohesion_leak.py` l. 309–311, the counts (5, 4, 4, 2)). 4/15 is the nearest live class's
  mass, the w of (v)'s time scale. With w = 4/15 and d = 5.19 the time scale is 1.39·10⁷ (R = 1) and 1.76·10⁴ (R = 1.5), the text's
  1.4·10⁷ and 1.8·10⁴; with the dead class's 2/15 it would be twice these (`calc32.out`).
- *D‴.34 (D).* "(without the replicator 1.9·10⁻⁷–2.3·10⁻⁶)" is the range of D‴.32 (C)'s six m = 60 seeds. For the same seeds 0–2 it is
  1.9·10⁻⁷–9.0·10⁻⁷, so "as without" holds even more closely.

*Fix:* "(amplitude 0, mass 2/15, at 5.19 from the nearest live class, of mass 4/15)"; compare (D) with seeds 0–2.

## Independent confirmations

1. **(i)** (`check_i.py`): see the verdict table. The identity, the FR projection, the sphere form and the bound hold for both kernels.
   My drift and J equal the repository's to 2.8·10⁻¹⁷ and 0.
2. **(ii)** (`check_ii.py`, my RHS):
   - V from 0.04181 to 2.140·10⁻⁷, monotone;
   - max|V̇ + k_ee²|/max k_ee² = 1.11·10⁻¹⁶; Σβ − m ≤ 1.4·10⁻¹⁴;
   - β at the end (1.00035, 0.99873, 1.00092), the md's to every printed digit.
3. **(iv), Part A** (`check_A.py`): the four runs reproduce the md's table:
   - γ = 0: 0.391, 0.388, 0.388, |a|/a_max (0.999, 0.425, 0.826, 0.383, 0.999);
   - γ = 10⁻²: 0.00238, 1.26·10⁻⁹, 2.99·10⁻¹⁴, masses (0.234, 0.199, 0.191, 0.168, 0.208); the geometric factor 6.07 per 2000 slow
     units (the text's "about 6");
   - scale 1.25, γ = 0: 0.496, 0.495, 0.494; γ = 10⁻²: 0.0872, 0.0735, 0.0711, |a|/a_max ≥ 0.994.
   The Gram matrix of the five atoms is nonsingular (eigenvalues 0.064–0.192), and Σ|A_k|/a_max = 0.875 and 1.094.
4. **(v), Part E** (`check_v.py`): see the verdict table. The md's γ = 0.01 range [−1.8337·10⁻¹, −3.0812·10⁻⁴] and the 15 zeros
   reproduce, and at γ = 0 the 36 separations equal −σ₀w_k with the 9 collective eigenvalues the unregularised flow's (−0.18336 …
   −3.0798·10⁻⁴).
5. **(B), seed 1, γ = 10⁻²** (`broad.py`, my implementation from `default_rng(1)`): after R = 2 the fit error is 0.0043 with 5 classes;
   at the end 5.99·10⁻⁴, the masses 0.262, 0.157, 0.206, 0.133, 0.242, the centres −1.600, −0.798, +0.004, +0.798, +1.597, the
   particles 20, 6, 10, 6, 18, |a|/a_max ≤ 0.958. All match the md. "Five classes in 3 of 3, the centres the bumps' to 0.01" holds in
   the md (largest deviation 0.010, seed 2's +0.010).
6. **(F), three heavy, seed 0** (`broad.py`): three classes, the fit error 1.04·10⁻⁷, the masses 0.411, 0.309, 0.280 above
   0.400, 0.300, 0.250, all as in the md. The fit error is still falling (1.95, 1.43, 1.04·10⁻⁷ at s = 3.8, 3.9, 4.0·10⁴: a transient, as printed).
7. **(C), γ = 10⁻², seed 0** (`partC.py`, my ϑ^b(q), 3.0·10⁻⁹ from the round-23 referee's independent one, its fit error 0.01695): 8
   classes, the fit error 0.005483, the masses 0.327, 0.034, 0.102, 0.202, 0.232, 0.033, 0.036, 0.034, and the cluster counts along the
   run the raw record's at every one of the 51 records.
   - *"The fitness vanishes at a critical point".* At a critical point of L_w every gradient component vanishes, ⟨r∇Φ̃_k⟩ = 0, so
     u_k = a_k⟨rφ_k⟩ = 0 = λ. The separations change neither u (⟨r∇Φ̃_k⟩ = 0) nor r (first order), and within-class redistributions
     leave r unchanged. So the linearisation's separation block is the replicator-free one at β ≡ 1.
   - The replicator acts only after the split (the end masses 0.327–0.036 against 0.333 and 0.033). The explanation is right.
8. **(F), two close bumps** (raw record): one class of mass 0.881, 0.971 and 1.000 at c = 0.000, +0.035, −0.000, widths 0.888, 0.706,
   0.667, the extra classes single saturated particles. Without the replicator (D‴.32 (K)) the widths are 0.716, 0.690, 0.667 with
   the masses 0.967–1.000. The broad stage merges the bumps in both. The masses of the merged class far exceed the needed 0.2 + 0.2,
   so the failure is the scale's, not the masses', as the text says.
9. **D‴.32's round-23 revisions** (`calc32.py`; the md):
   - *The twins' rate.* σ₀ψ(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²) re-derived from ṗ = 0 (to 4·10⁻¹⁶ at random parameters).
   - *Its time scale.* G_αR²e^{z₀}/(σ₀(w_k + w_l)Δα²), z₀ = Δα²/(2G_αR²), from dz/ds = −2σ₀(w_k + w_l)ze^{−z} for a′ ≈ 1. The exact
     integral exceeds it by 1.29, 1.13, 1.06 at z₀ = 5, 10, 20 (the stated 1 + O(1/z₀)).
   - *The dead class's time scale.* It reproduces the text's numbers with w = 4/15 (ID-9).
   - *Part I.* The 'cx' run (12.7 → 11.5–12.9 over 4·10⁵, the fit error 3.16·10⁻⁴ → 4.54·10⁻³), its annealed version (11.7, 4.99·10⁻³),
     the far twins (5.28 and 8.19 unmoved) and the dead class at R = 1 (unmoved) and 1.5 (joined between s = 2.4·10⁴ and 2.8·10⁴) are as
     the text states.
   - *Parts J and K.* The text's ranges are the md's (J: 4–6 classes, 4.84·10⁻⁴–1.59·10⁻³ at R = 1, 7.82·10⁻⁴–1.40·10⁻³ at R = 2; K:
     two close 0.0983–0.13, widths 0.667–0.716, margins 6.46·10⁻³–7.07·10⁻³; three heavy 0.0577–0.204).
   - *Part L.* The three-equal residual is the Gaussian bias: the stationary point at R = 1 with that run's counts has 9.387·10⁻⁵ (the
     run's 9.39·10⁻⁵ for 2·10⁵ slow units), and the compact kernel removes it. The four-alternating decay is geometric at
     1.49·10⁻⁵ (Gaussian, last 10⁴) and 1.50·10⁻⁵ (compact stage), the weakest collective rate the round-23 report found
     (1.5·10⁻⁵): a slow transient.
10. **The bibliography**: two items checked by a web search (ID-7).

## Not checked

- **Parts B (other seeds and γ = 0, 10⁻³), D and C (γ = 0; seeds 1, 2) and F (seeds 1, 2; two close).** Not rerun; their values come
  from the raw record. Every run I did reproduce (A ×4, B seed 1, C seed 0, F three heavy seed 0, E) matches it to the printed digits.
- **The real loop with transport and replicator**, and the "manifold of periodic orbits" of (v): argued (e ≡ 0 on exact
  representations), not integrated.
- **The positive eigenvalues (≤ 1.3·10⁻⁷)** of the particle Jacobian at the near-stationary five-bump end under the compact kernel
  (ID-2). With unequal weights the separations' margin becomes σ₀n_k/m against β_i times the residual's curvature. I did not identify
  the directions.
- **The weighted splitting margin in general.** D‴.32 (iv)'s σ₀w_k − (q_k/w_k)λ⁺_k assumes β ≡ 1. After a transit with the replicator
  (β ≠ 1) the margin carries the particle fraction n_k/m in the cohesion and β_i in the fit. Neither the text nor I worked it out.
- **Łojasiewicz on the constrained manifold** (iv): checked only as hypotheses (analytic J and metric in (log W, α)).
- **The literature's full texts** (ID-7).
- **The concurrent session's** uncommitted edits.

