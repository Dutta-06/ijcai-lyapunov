# Referee report D3u (round 23) — PARADIGM §5d: Proposition D‴.32 (the cohesion leak: statement (i)–(vi), proof, Status, *Numerics*), the §0 row and the §7 row "| D‴.32 |", the §8 bullet "Done (Proposition D‴.32)", `results/RESULTS_PARADIGM.md` (Proposition D‴.32), LITERATURE_TRANSPORT [S60] and the audit bullet "Proposition D‴.32", PAPERS.md (paper 13), README.md where it mentions the cohesion leak, `src/core/cohesion_leak.py`, `experiments/exp_cohesion_leak.py` → `results/chaos/cohesion_leak.md` (+ raw), `tests/test_cohesion_leak.py`

Branch `claude/festive-bardeen-eopej2` at 64cb25e, HEAD when the review started and when it ended (`git log --oneline -5`: 64cb25e, a9a0e5a,
bfdbc0f, 29f9394, e325236).

The tree was extracted with `git archive 64cb25e` into `ref_D3u/tree`, and everything below was run there. During the review the
concurrent session was running `experiments/exp_eiv_remedy.py --only-D2` (four processes). At the end the shared tree had uncommitted
edits to `experiments/exp_eiv_remedy.py`, `src/core/measurement_noise.py` and `tests/test_eiv_remedy.py`, and untracked
`src/core/mass_transfer.py`, `experiments/exp_mass_transfer.py`, `tests/test_mass_transfer.py` and `results/chaos/mass_transfer_raw.npz`.
None of this is reviewed here.

Every script imports nothing from the repository, with three exceptions. The cross-checks are marked (the repository's object compared
against mine). `check_tests.py` and `tests_values.py` run the repository's tests against mutants of its code. `g_track.py` was run with
the repository's ϑ^b(q) as data, and for seed 0 also with mine, with identical results. My
instance (`inst.py`) is an independent implementation. It has the saturated network, analytic gradients and Hessians of Φ̃, the averaged
learning drift, both kernels and their energies, analytic Jacobians of the particle drift, the whitened weighted Hessian Ã_w, the
Gauss–Newton block H_L(w), single-linkage clustering and a midpoint integrator, all my own.

Line numbers are 64cb25e's:
- PARADIGM: the statement l. 3934–3940, (i) l. 3942–3951, (ii) l. 3953–3960, (iii) l. 3962–3968, (iv) l. 3970–3973, (v) l. 3975–3981,
  (vi) l. 3983–3993; the proof l. 3995–4022 (of (v) l. 4011–4018); the Status l. 4024–4035; the *Numerics* l. 4037–4095 ((A) l. 4038–4052,
  (B) l. 4052–4056, (C) l. 4057–4062, (D) l. 4063–4066, (E) l. 4066–4070, (F) l. 4070–4075, (G) l. 4075–4086, (H) l. 4086–4095).
  The §0 row l. 139, the §7 row l. 6599, the §8 bullet l. 6872–6882.
- RESULTS l. 815–848. LITERATURE_TRANSPORT [S60] l. 80, the audit bullet l. 485–497. PAPERS.md l. 211–235. README l. 67–70.
- Code: `cohesion_leak.py` l. 1–30 (docstring), 37–52 (`psi`, `omega`), 55–70 (`pair_d2`, `cohesion_force`, `interaction_energy`),
  73–78 (`clusters`), 81–125 (`CohesionLoop`), 128–155 (`predicted_spectrum`), 158–176 (`CohesionRealLoop`);
  `exp_cohesion_leak.py` l. 1–71 (header), 142–181 (class drift, Ã_w, Newton), 184–318 (Part A), 322–372 (B), 376–421 (C),
  425–476 (D), 479–509 (E), 513–539 (G), 543–562 (H), 566–604 (F), 607–763 (the report); `tests/test_cohesion_leak.py` l. 52–161.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`); pytest with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the
scratch tree rewrites `cohesion_leak.md` byte-identically (md5 `c7848f6102b148a73b99367526f81d48`, equal to
`git show 64cb25e:results/chaos/cohesion_leak.md`).

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3u/scripts/`:

| script | what it does |
|---|---|
| `inst.py` (module) | independent instance: the saturated network, ∇Φ̃ and ∇²Φ̃ analytic, the averaged drift, both kernels ψ, ω, the cohesion force and J^coh_m, the analytic Jacobian of the particle drift (learning and kernel parts), R_k, Ã_w, H_L(w), single linkage, the midpoint rule; a fast drift for long runs (checked against the slow one to 3·10⁻¹⁶) |
| `check_i.py` | (i): ∇Φ̃, ∇²Φ̃ against differences; ω′ = ψ/2, ∇²W(0) = G⁻¹; the drift against −G∇J^coh_m (my J); the forces' sum; the mean-shift identity; the mean-field velocity at the particles; the barycentre's velocity; the kernels' Lipschitz constants and the compact pair force's maximum; cross-check of my drift and J against the repository's |
| `check_A1.py` | (ii)–(iii) at the 15 count vectors of m = 15: drift, inter-atom distances, the analytic Jacobian's spectrum against −σ₀w_k (separations) and −eig H_L(w) (collective), Ã_w = H_L(w) at r = 0, the leak's Jacobian on the class means, the mixing blocks |
| `check_A2.py` | (ii): the Gaussian kernel's replicated stationary point by Newton (analytic class Jacobian), ‖y‖, Ã⁻¹b, (Ã − Db)⁻¹b, the stated order estimate, the fit error, the particle Jacobian's top eigenvalue; R = 1, 1.5, 2, 3 at (5, 5, 5) and (6, 5, 4) |
| `check_A3.py`, `check_A3b.py`, `cross_bq.py` | (iv): ϑ^b(q) by Newton continuation in ε (cross-check: the repository's BFGS approximant, 5.7·10⁻¹¹ apart); Ã at ϑ^b(q); λ⁺_k; the gauge scaling of R_k; the margins against the analytic Jacobian at (5, 5, 5), (6, 4, 5), (4, 6, 5) and σ₀, σ₀/8; the three split runs reproduced from `default_rng(0)` with the top separation component's growth rate |
| `check_A4.py` | (v): the md's twins (3:1) and dead-class runs reproduced; twins with the shares 0.9 and 0.97; the dead class at R = 1, 1.5, 2 |
| `twins_cx.py`, `twins_hess.py`, `twins_sched.py` | (v): a twin configuration with opposite-sign twins, one at the amplitude coordinate of the class of bump 2; J^coh_m along the gauge curve; the class Hessian's spectrum there; trust-region minimisation of the class energy from it; a 4·10⁵ slow-unit particle run (Gaussian, R = 1); the annealed schedule R = 2 → 1 |
| `dead_modes.py` | (v): the modes of the dead class's mixture (weights w_le^{−α_l²/(2R²)}) on the instance and on a near-equilateral configuration, R = 0.75–3 |
| `shed_energy.py` | (G): the cohesion energy of a class split by shedding k particles at equal Σ‖u_i‖² (k = 1, 2, 3, 5), both kernels |
| `broad_run.py`, `batch.sh`, `list1.txt`–`list4.txt`, `margins_end.py`, `t4rep_bias.py` | (C), (H): broad starts reproduced from the seed (`default_rng(seed)`, the experiment's draw order), any schedule of kernels and scales; the R-window, other targets, other m and seeds; the end states' splitting margins; the Gaussian bias on my four-bump target |
| `check_B.py` | (B): seed 3, compact kernel, the four systems m = 15, 60, 240, 510 reproduced from the seed |
| `check_D.py`, `check_D2.py` | (vi), (D): the real loop's right-hand side, the lift's return, the block structure of the variational generator, the 46-dimensional monodromy and the 10-dimensional collapsed one (matrix ODE, DOP853), the collapsed rates against h_j(w)/k_e; D2's warm-cloud run (m = 15, cohesion) reproduced |
| `check_E.py` | (E): the dither run seed 0, compact kernel, reproduced from `default_rng(100)` |
| `g_track.py`, `g_tally.py`, `raw_G.py` | (G): the σ₀ sweep reproduced from `default_rng(700 + seed)` with every class tracked through its stage (split or not, the pieces' sizes); the raw record read stage by stage |
| `check_tests.py`, `tests_values.py` | the tests: their asserted quantities; 12 mutants of `cohesion_leak.py` patched in memory |

Also in `ref_D3u/`: `tests.out` (5 passed, 4.42 s), `from_raw.out`, `md5_before.out`, `md5_after.out`.

## Verdicts

| object | verdict |
|---|---|
| **(i) the energy and the mean-field form** | **holds.** Re-derived with the repository's normalisation (J = mL̄ + (σ₀/(2m))ΣW, f̂ = (1/m)ΣΦ̃, the per-particle law without 1/m). G∇W(u) = ψ(z)u with ω′ = ψ/2, ∇²W(0) = G⁻¹, and the leak σ₀ψ̄_i(m_i − ϑ_i) = (σ₀R²/2)G∇p̂ with p̂ the kernel density of the shadow profile. All confirmed to ≤ 8·10⁻¹⁰ against differences. The empirical measure solves the mean-field equation exactly (the j = i term vanishes). The pair force is Lipschitz (constant 1, whitened), so Dobrushin applies. The forces sum to zero. Nothing confines the unconstrained coordinates: the counterexample of D3u-1 runs toward the saturation, an unbounded solution, so "every bounded solution converges" is the right hedge. The mean-field wording is imprecise (D3u-8). |
| **(ii) no bias at clustered configurations** | **holds** with the finite range (stationary iff critical for L_w, any K, any w, atoms ≥ R apart). The Gaussian kernel's continuation Ã_w⁻¹b + O(\|b\|²) and every number of Part A2 reproduce. The stated order of the bias is an upper estimate, 7.6–8.3 times the actual on the instance. |
| **(iii) the blocks** | **holds exactly.** Separations −(σ₀w_kI − R^G_k) with multiplicity n_k − 1, the mass factor right (the class's share of the particles in range), class means −Ã_w, no mixing, the leak's Jacobian zero on the class means. ker Df̂ = separations iff H_L(w) ≻ 0 (true at all 15 count vectors). The response Ã_w⁻¹F agrees with D‴.31(ii), whose S₀ = 0 for every η > 0, not only η = ∞. The mass factor is a design property, not a flaw of the proof: light classes are held less, split more readily (√ of the adaptive prior's thresholds) and spread more under noise and dither. |
| **(iv) splitting** | **holds for the linearisation.** The margin σ₀w_k − (q_k/w_k)λ⁺_k is exact at gauge transforms (the scaling to 9·10⁻¹⁴), "never over-split" holds as linear stability. In Part A3 the split is the predicted instability at the predicted rate (1.224·10⁻⁴ against 1.218·10⁻⁴). Part G reproduces in all 21 stages, and its tallies hold: 0 of 48 positive-margin classes split, 8 of 8 below −5·10⁻⁵, 3 of 7 between. Two explanations are missing (D3u-6). |
| **(v) minimality** | **not proved.** The twins argument is wrong: J^coh_m along the gauge curve can have a local minimum away from the merger, and twins need not merge (a counterexample on the instance; far twins do not merge in 4·10⁴ slow units) (D3u-1). The dead-class argument holds only for well-separated live classes and on a time scale e^{d²/(2R²)}/(σ₀w). Its numerical check is at R = 2, outside that regime (D3u-3). |
| **(vi) the real loop** | **holds.** The lift is periodic. The separations' exponents σ₀w_k/k_e are exact, the subspaces invariant, the collapsed block the unregularised K-atom loop's (D‴.26 at η = ∞). Part D1's 45 exponents are confirmed to 5.2·10⁻¹¹ by a variational monodromy. The noise statements are D‴.30's, formal as there. Their threshold is a linear one: under the compact kernel the restoring force on a separation is bounded (the pair force ψ(d²/R²)d is maximal, 0.238R, at d = R/√7 and zero at R), so a separation pushed beyond R is not restored at all. |
| Proof | (i)–(iv), (vi) **correct**; (v) **wrong for twins** (D3u-1), **incomplete for the dead class** (D3u-3). |
| *Status* | accurate except "(v) proved for twins" (D3u-1) and "checked numerically" for the dead class (D3u-3). |
| *Numerics* | **the printed numbers reproduce**: A1–A4, B (seed 3), C (three runs), D1, D2, E (seed 0), G (all 21 stages), H (five bumps, seed 0). The readings overstate in three places: the count selection (D3u-2), (B)'s m = 15 range (D3u-4), and finite-difference precision quoted as agreement (D3u-7). |
| §0 row, §7 row, §8 bullet, RESULTS, PAPERS, README | numbers hold; "twins and dead classes merge" (D3u-1); "selects the number of kernels", "from every broad start" (D3u-2); README combines the compact kernel's "without biasing any clustered critical point" with the Gaussian kernel's merging (D3u-2). |
| LITERATURE [S60], audit bullet | bibliographic details right where checked; precedents omitted (D3u-10). |
| `cohesion_leak.py` | **holds.** `psi`, `omega`, `pair_d2`, `cohesion_force`, `interaction_energy`, `CohesionLoop`, `predicted_spectrum` and `CohesionRealLoop` agree with my implementation (the drift to 4.7·10⁻¹⁶). The docstring's "the force σ₀w_kw_le^{−d²/(2R²)}d" is the class-mass-weighted force (the per-particle one is σ₀w_l…), harmless. |
| `exp_cohesion_leak.py`, md, raw | **holds**; the regeneration is byte-identical; every reproduced run matches the raw record. |
| `tests/test_cohesion_leak.py` | **pass (5/5)**; they catch sign and mass-factor errors but not the kernel's range or scale (D3u-5). |

## Confirmed issues

**D3u-1 (major; (v) l. 3977–3981 "The twins' attraction decreases J^coh_m along the curve all the way to their merger, so a configuration
with twins is not a local minimum … So, under the Gaussian kernel, the stable configurations near exact representations carry no dead or
duplicated classes: the number of classes is the representation's"; the proof l. 4011–4015; the Status l. 4024–4025 "(v) proved for
twins"; the §0 row l. 139 (statement and Status columns); the §8 bullet l. 6877; RESULTS l. 822; PAPERS l. 218–219, 222–223; README
l. 68–69; the audit bullet l. 495).** The twins argument is wrong, and twins do not merge in general (`twins_cx.py`, `twins_sched.py`,
`check_A4.py`).
- *The gap.* The proof compares the twins' own attraction with "the other classes' attraction … only at the order of their
  exponentially small forces". Both are exponentially small: the twins' pull is ∝ e^{−|α_k − α_l|²/(2R²)}, a third class's pull on one
  twin ∝ e^{−d²/(2R²)}. When a third class is closer to one twin than the twins are to each other, the third class wins, and J^coh_m
  restricted to the gauge curve need not be monotone.
- *A counterexample on the instance* (Gaussian kernel, R = 1, σ₀). The counts are (5, 3, 5, 2)/15. Bump 3 is carried by twins k, l
  (masses 5/15 and 2/15, the same centre 1.5 and width 0.5) with the amplitudes +3.8 and −3.5 (|a|/a_max ≤ 0.95). Twin l sits at the
  amplitude coordinate of the class of bump 2 (−3.5): it is 3.63 (whitened) from that class and 12.74 from its twin. The fit error is
  1.7·10⁻¹⁶: an exact, degenerate representation.
  - *J^coh_m along the gauge curve* (exact representations, a_k = (0.8 − w_la_l)/w_k), relative to its value at a_l = −3.85:
    −3.42·10⁻⁶ at a_l = −3.59, −3.06·10⁻⁶ at −3.33, −4.49·10⁻⁷ at −2.81, +9.6·10⁻⁸ at −2.30, +1.69·10⁻⁷ at −1.52, −2.4·10⁻⁸ at −1.00,
    −4.72·10⁻³ at the merger (a_k = a_l = 1.714). The curve has a local minimum near the start, behind a barrier of 3.6·10⁻⁶. The
    sentence "decreases J^coh_m along the curve all the way to their merger" is false here.
  - *Not a local minimum, for another reason.* The class-level Hessian of J^coh_m at the twins has the eigenvalue −1.8·10⁻⁵ (whitened;
    the next ones +5.9·10⁻⁷, +1.5·10⁻⁶). Its eigenvector is the twins' split direction in (γ, β): both centres move toward class 2, twin k
    by 0.37 of twin l's step, which cancels f̂'s first-order change for opposite-sign amplitudes. Along it the fit is only quartic, and
    the attraction from beyond R has negative radial curvature, so the twins are a saddle. A trust-region minimisation from them lowers
    the class energy by 5.4·10⁻⁴. The proof does not consider this direction.
  - *The dynamics.* From the twins with a jitter 10⁻³ the particle flow (m = 15) never merges them. Over 4·10⁵ slow units their distance
    stays 11.5–12.9, four classes throughout. Twin l moves to c = 0.88 and a = −3.87 (toward the saturation), class 2's amplitude to
    −3.87. The fit error rises from 4.1·10⁻⁴ (the jitter) to 4.5·10⁻³. Under the annealed schedule of Part C (R = 2 for 2·10⁴, then
    R = 1 for 2·10⁴) the end is the same: four classes, the fit error 5.0·10⁻³. The leak turns an exact (degenerate) representation into
    a biased four-class fit.
- *Same-sign twins farther apart do not merge on the numerics' time scale* (R = 1, the md's masses 3/15 and 2/15). Bump 3 shared
  0.9 : 0.1 (amplitudes 3.6 and 0.6, distance 5.28): 5.28 → 5.29 over 4·10⁴ slow units, the fit error 1.6·10⁻⁶. Shared 0.97 : 0.03
  (distance 8.19): 8.19 throughout, 5.4·10⁻⁶. The md's twins (3 : 1, distance 2.31) merge by s = 1.8·10⁴, which I reproduce. They are
  the favourable case: twins within about 2R of each other in the amplitude coordinate.
- *What holds.* For twins with no other class nearer to either twin than they are to each other, the twins' pair energy is the
  dominant part of J^coh_m's variation along the curve. The merger is then reachable inside the saturation (verified). Its rate is
  ∝ σ₀w e^{−|α_k−α_l|²/(2R²)}, so even then "merge" holds only on a time scale e^{|α_k−α_l|²/(2R²)}/(σ₀w).

*Fix:* restate (v) for twins as a conditional statement: twins whose own pair force dominates the variation of the other pair energies
along the gauge curve are not local minima, and the flow merges them at the rate above. Drop "(v) proved for twins" and "twins and dead
classes merge" from the Status, the §0 row, the §8 bullet, RESULTS, PAPERS, README and the audit bullet, or qualify them the same way.
Add the counterexample and the far-twin runs to Part A4.

**D3u-2 (major; the title l. 3934 "a leak … that selects the number of kernels"; the §0 row l. 139; the §8 bullet l. 6872 and 6877–6878
"the annealed Gaussian kernel finds the target's three kernels from every broad start"; README l. 67–70 "holds the separations without
biasing any clustered critical point … merges duplicated or dead ones, so the network selects its number of kernels"; PAPERS l. 211–219;
*Numerics* (C) l. 4057–4062 and (H) l. 4093–4095 "So the count is the target's when the transit's masses are representable and the
approximant's below the resolution").** The selection is a property of the annealed schedule, not of the leak, and it is not proved.
The two properties the summaries combine belong to different kernels.
- *What is proved.* The proposition proves no selection. It proves the splitting margin's linearisation (iv) and a merging statement
  (v), which fails as stated (D3u-1). The transit is not proved, as the Status says.
- *The kernel and the scale decide* (`broad_run.py`, broad starts from the seeds of Part C, m = 60):
  - Gaussian at R = 1 alone for 4·10⁴ slow units: 5, 6, 5, 5, 4, 5 classes at seeds 0–5, i.e. 0 of 6 at the target's three. The fit
    errors are 4.8·10⁻⁴–1.6·10⁻³. The extra classes are small fragments (one to five particles, |a|/a_max ≤ 0.42) far from the live
    classes, which the R = 1 kernel no longer reaches.
  - Gaussian at R = 2 alone: three classes, but the fit error is 7.8·10⁻⁴ (the bias).
  - Annealed 1.5 → 1, 2 → 1 (the md's), 3 → 1 and 2 → 0.5: three classes, the fit error 2.4·10⁻⁸, 1.9·10⁻⁷, 5.9·10⁻⁶, 1.9·10⁻⁷
    (seed 0).
  - With the compact kernel the md already reports 17–22 classes.
  So the count is selected by a broad stage at a scale comparable to the classes' distances, i.e. at a biased fit, followed by a fine
  stage. Neither stage is analysed. The md does not report the Gaussian kernel at a fixed scale.
- *Other targets and starts* (the annealed schedule 2 → 1, m = 60 unless stated):
  - *The instance, other starts:* three classes at m = 30 (seeds 0, 1: the fit error 2.4·10⁻⁷, 6.1·10⁻⁶), m = 240 (seed 0: 6.7·10⁻⁷)
    and m = 60, seeds 6, 7 (3.5·10⁻⁷, 5.6·10⁻⁷). On this target the annealed schedule is robust: 14 of 14 runs counting the md's.
  - *Four well-separated bumps of alternating sign* (A = 0.5, −0.45, 0.5, −0.4 at c = ±1.6, ±0.55, s = 0.35; Σ|A_k|/a_max = 0.46):
    four classes in seeds 0, 1, 2. The fit errors are still 1.1·10⁻², 9.6·10⁻³, 1.0·10⁻² at the end. The stationary bias at R = 1 is
    only 1.9·10⁻⁶ (`t4rep_bias.py`), but the weakest collective rate is 1.5·10⁻⁵, so 2·10⁴ slow units are far from enough.
  - *Three equal same-sign bumps* (0.8 at c = −1.5, 0, 1.5, s = 0.45): three classes, 9.4·10⁻⁵ (seed 0).
  - *Three bumps needing 95% of the mass* (A = 1.6, −1.2, 1.0; Σ|A_k|/a_max = 0.95):
    - seed 0: three classes, but the masses 0.467, 0.283, 0.250 do not represent it (bump 1 needs w > 0.40, bump 2 w > 0.30). Two
      classes are saturated, and the fit error is 0.204.
    - seed 1: five classes (bump 3 carried by a pair, a singleton), 0.113.
    Part H's five-bump failure needs only three bumps.
  - *Two same-sign bumps 0.9 apart* (0.8 at c = ±0.45, s = 0.4; Σ|A_k|/a_max = 0.40, representable at any w_k > 0.2):
    - The broad stage merges them into one wide class (mass 0.967, c = 0.000, s = 0.716). Two saturated singletons of negative amplitude
      sit at c = ±1.42. The fit error is 0.098 (seed 0).
    - The merged class's margin is +6.8·10⁻³ with λ_max(R^G) = 6.8·10⁻⁸ (`margins_end.py`): by (iv) it never splits, whatever σ₀. The
      leak's resolution cannot see the missing kernel there.
    - Smaller scales do not help. 1 → 0.5 ends at nine classes (0.036), 0.6 → 0.3 at sixteen (0.021), both with a wide class
      (s ≈ 1.0) at c = 0. Seed 1 (2 → 1) ends at two classes, the target's count, but they are the wide class (mass 0.983, s = 0.69) and
      a saturated singleton, with the fit error 0.117. Without a leak (seed 0, 4·10⁴ slow units) the fit error is 6.3·10⁻³ at 50 classes:
      the leak makes the fit fifteen times worse here.
- *The two kernels.*
  - "No bias at every clustered critical point" holds for the compact kernel, which merges nothing (twins and dead classes persist
    beyond R).
  - The Gaussian kernel merges, but only on the time scale e^{d²/(2R²)}/(σ₀w). It biases the fit by ∝ σ₀‖Ã_w⁻¹‖e^{−d²_kl/(2R²)}, the
    same small parameter.
  - README's sentence gives one leak both properties.

*Fix:* the title and the summaries should say what is shown. Under an annealed Gaussian scale the leak ended at the target's count
on well-separated targets whose transit masses were representable (the md's two- and three-bump targets, my four bumps and three equal
bumps). It failed when two bumps are closer than the broad scale (my two bumps 0.9 apart) and when the transit's masses are not
representable (five bumps; three bumps needing 95% of the mass). Report the fixed-scale Gaussian runs and the scale window in Part C.
"From every broad start" → "in the 9 broad starts of Part C". State in README which kernel has which property.

**D3u-3 (minor; (v) l. 3978–3980 "a dead class moves under the attraction of the live ones, whose local minima — when the live classes are
farther apart than a few R — lie at the live classes, where it becomes a twin and merges"; the proof l. 4015–4018; the Status l. 4025
"checked numerically"; *Numerics* (A) l. 4050–4052).** The regime of the argument and the regime of the check do not overlap
(`check_A4.py`, `dead_modes.py`).
- *The check is outside the hypothesis.* The md's dead class merges at R = 2, which I reproduce (5.19 → 0.33 at s = 4000, to rounding
  by s = 1.6·10⁴, the fit error then 1.64·10⁻³). At R = 2 the live classes are 3.63 apart in (γ, β) (whitened), i.e. 1.8R. The
  components overlap, and the dead class's mixture has a single mode, 0.35 from the class of bump 2.
- *Where the hypothesis holds, nothing moves.* At R = 1 (3.6R apart) the dead class stays at 5.19 from the nearest live class for
  4·10⁴ slow units (the fit error 7.0·10⁻⁷ at the end, the bias). At R = 1.5 it joins at s ≈ 2.8·10⁴. The pull is
  ∝ e^{−d²/(2R²)}: at d = 5.19 and R = 1 it is e^{−13.5}.
- *The weights.* The mixture's weights are w_le^{−α_l²/(2R²)}, not w_l: the dead class sits at α = 0, and the live classes' amplitude
  coordinates enter. At w = q and R = 1 the class of bump 1 (α = 3.89) has no mode, though it is 3.6R from the class of bump 2. So
  "its local minima lie within O(R e^{−d²/(2R²)}) of the components" needs a condition on the weights as well as on the distances.
- *Carreira-Perpiñán–Williams.* An extra mode does not occur on the instance: its triangle (3.63, 3.63, 6.59) gives one mode per
  component at R ≤ 1 and a single mode near the class of bump 2 at R ≥ 1.5. It does occur at the broad stage's scale on a near-equilateral
  target. With c = (−1, 0, 1), s = (0.4, 0.8, 0.4) and A = 0.6 at w = q, the triangle is (4.33, 4.16, 4.33). At R = 2–3 the mixture is
  unimodal, its mode 2.34 (R = 2) to 2.42 (R = 3) from the nearest component, where a dead class is a twin of nothing.

*Fix:* state the time scale e^{d²/(2R²)}/(σ₀w) and the weights' condition; check the dead class at a scale where the live classes are
"farther apart than a few R" (or say that the merger is observed only at R = 2, where the components overlap and the fit is biased).

**D3u-4 (minor; *Numerics* (B) l. 4052–4056 "both kernels remove the floor as the adaptive prior does … 1.9·10⁻¹⁰–7.1·10⁻⁵ … one compact
run at m = 15 ends with a fourth class (a class diameter of 0.81)"; the §7 row (B); RESULTS l. 823–824).** The top of the m = 15 range
and the fourth class come from the one sample outside W, which no leak repairs (`check_B.py`, the raw record).
- *The sample.* Seed 3's m = 15 counts are (8, 4, 3): w₃ = 3/15 = 0.2 = |A₃|/a_max, on the boundary of W, with no exact
  representation.
- *What every leak does there.* At s = 6000, ‖f̂_m − f̂_M‖² is 7.10·10⁻⁵ (compact), 7.04·10⁻⁵ (Gaussian) and 7.04·10⁻⁵ (adaptive
  prior). The fit error is 1.300·10⁻², 1.295·10⁻² and 1.295·10⁻²: no floor is removed.
- *The other three seeds.* At m = 15 the range is 1.9·10⁻¹⁰–5.4·10⁻⁸ under every leak.
- *Reproduction.* My run of seed 3 (compact, the four systems from the seed) gives 7.096·10⁻⁵ and 1.300·10⁻² at m = 15, the class
  diameter 0.81 and four clusters. At m = 60 and 240 it gives 6.451·10⁻¹⁰ and 1.628·10⁻¹⁰, and the class means 8.42·10⁻⁴ and
  1.75·10⁻⁴ from ϑ^e(w^m). All are the raw record's values.

*Fix:* "… at the representable samples (at m = 15: 1.9·10⁻¹⁰–5.4·10⁻⁸); seed 3's m = 15 sample (8, 4, 3) lies on ∂W and keeps the fit
error 1.3·10⁻² under every leak (‖f̂_m − f̂_M‖² 7.0·10⁻⁵–7.1·10⁻⁵), and there the compact kernel ends with a fourth class".

**D3u-5 (minor; `tests/test_cohesion_leak.py` l. 52–161; the Status's "Checked" l. 4030–4035; RESULTS l. 841–845; the §7 row's last
column).** The tests pass (5/5, 4.42 s) and catch sign and mass-factor errors, but no test pins the kernel's range or scale
(`check_tests.py`, `tests_values.py`).
- *Asserted quantities.* Test 2: the drift 6.1·10⁻¹⁸ and 2.8·10⁻¹⁷ (< 10⁻¹⁴), the spectrum error 1.8·10⁻¹¹ and 1.9·10⁻¹¹ (< 10⁻⁹,
  finite differences), −coll[−1]/σ₀ = 0.049 and 0.047 (< 0.2); the classes 6.6 apart. Test 4: the twins start 1.54 apart; after 400 slow
  units at 10σ₀ the Gaussian twins are at 0.0002·d₀ (asserted < 0.9·d₀), the fit error 4.4·10⁻⁴ (< 2·10⁻²).
- *Caught.* The mass factor dropped, as a pure mean-shift step σ₀(m_i − ϑ_i) (tests 1, 2, 5) or consistently in force and energy (1, 2,
  3, 5); the sign flipped, in the force alone or in both (all five); σ₀ in place of σ₀/k_e in the real loop (5); a leak on the class
  means (all five).
- *Not caught.* Consistent mutants of ψ and ω (ω′ = ψ/2 kept) pass all five tests:
  - the compact range 0.5R, √2R or 1.5R;
  - the Gaussian scale R/√2 (e^{−z}) or 2R (e^{−z/8}).
  Test 1 checks only the force–energy consistency. Test 2's classes are 6.6 apart, so any range below 6.6 passes. Test 4's compact twins
  at 1.54 pin the range only below 1.54R. Test 3 has one class (w = 1), so it cannot see the mass factor. `clusters` is untested.
- *Loose.* Test 4 asserts only that the Gaussian twins approach (dg < 0.9d₀ against 2·10⁻⁴ measured) and a 2% fit; it does not test the
  merger's time scale, which is the point of D3u-1.

*Fix:* assert ψ at the range: a pair at 0.99R attracts, one at 1.01R does not, and the Gaussian force ratio is e^{−1/2} at d = R.
Test the twins' merger rate against σ₀w e^{−d²/(2R²)}, and `clusters` on a configuration with a known count.

**D3u-6 (cosmetic; *Numerics* (G) l. 4079–4086 "every class … whose margin was below −5·10⁻⁵ split within its stage (8 of 8), and of
those between −5·10⁻⁵ and 0 … 3 of 7. The split is not a halving: the unstable class sheds one to three single particles … σ₀ sets the
resolution as (iv) says"; the §0 row; the §7 row (G); RESULTS l. 835–837).** The counts hold, but the reading leaves out two
explanations it could give (`g_track.py`, `raw_G.py`, `shed_energy.py`).
- *The counts.* I tracked every class through its stage in all 21 stages. The counts at each stage's start and end equal the raw
  record, and the fit errors agree to the printed digits.
  - No class of two or more particles with a positive margin split (0 of 48).
  - All 8 with a margin below −5·10⁻⁵ split.
  - 3 of the 7 between −5·10⁻⁵ and 0 split: seed 0 at σ₀/32 and σ₀/64, seed 2 at σ₀/64.
  - Every split sheds one to three single particles and none is a halving: the pieces are [8, 1, 1], [9, 1], [8, 1], [7, 1], [9, 1],
    [8, 1], [8, 1, 1], [6, 1, 1], [8, 1, 1], [7, 1, 1, 1], [6, 1, 1].
  The paragraph's tallies are right.
- *The −5·10⁻⁵ criterion.* It is not ad hoc. A class counts as split when a piece passes the single-linkage threshold 0.25, and the
  stage starts from the jitter 10⁻³. That takes a growth by about 250 in 10⁵ slow units, i.e. a margin near −ln(250)/10⁵ = −5.5·10⁻⁵.
  Near that value the outcome depends on the jitter's projection and on the residual's drift during the stage:
  - −4.93·10⁻⁵ and −4.21·10⁻⁵ split (seed 0), −4.13·10⁻⁵ and −3.81·10⁻⁵ did not (seeds 1, 2);
  - −2.75·10⁻⁵ split in the same stage as another class's split (seed 2 at σ₀/64), which changes the residual.
  The paragraph gives the growth factors (1.2–140) but not this criterion.
- *The shedding is explicable.* The unstable eigenspace of a class is (n_k − 1)-fold degenerate: every zero-sum pattern of the
  particles along the top eigenvector grows at the same linear rate. The linear theory therefore does not choose the pattern; the
  kernel's concavity does. ω″ = ψ′/2 < 0, and the compact pair force ψ(d²/R²)d is maximal at d = R/√7 and vanishes at R. At equal
  Σ‖u_i‖² (equal linear fit gain), shedding one particle costs the least cohesion energy.
  - At Σ‖u_i‖² = 0.5 (n = 10, R = 1) the cost is 1.08 for one particle out, against 1.55, 1.74, 1.85 for 2, 3, 5 (compact). For the
    Gaussian it is 2.18 against 2.31, 2.36, 2.38.
  - At Σ‖u_i‖² = 10⁻⁴ the four patterns cost the same to 10⁻⁴ (relative): the linear degeneracy.
- *"σ₀ sets the resolution."* For the first split this is (iv) exactly. For the count it is not: 8–10 classes at σ₀/64 for four
  bumps. The threshold is σ₀w_k² rather than σ₀ (confirmation 5).

*Fix:* add the growth estimate behind −5·10⁻⁵, the degeneracy plus concavity explanation of the shedding, and "σ₀w_k² sets the
threshold of the first split".

**D3u-7 (cosmetic; *Numerics* (A) l. 4040, 4043, 4045 and (D) l. 4064–4065; the §0 row; the §7 row; RESULTS l. 818–821, 828–829).** Several
agreements quoted as properties are the precision of finite differences. With analytic Jacobians and a variational monodromy
(`check_A1.py`, `check_A3.py`, `check_D.py`):
- "(iii)'s spectrum to 2.8·10⁻¹¹" (compact) is 5.6·10⁻¹⁶;
- "the margins the Jacobian's to 2·10⁻¹¹" is ≤ 1.5·10⁻¹⁷ (the spectrum to ≤ 3.1·10⁻¹⁶);
- "the predicted ones to 5.5·10⁻⁶ (relative; median 1.5·10⁻⁷)" is 5.2·10⁻¹¹ (median 1.2·10⁻¹²). The separations' block of the monodromy
  is exp(−σ₀w_kT₀/k_e)I to 1.3·10⁻¹⁵, and the collapsed and separation subspaces are invariant to 6·10⁻¹⁷;
- "the first order Ã_w⁻¹b to 6·10⁻⁸" at R = 1 is a relative difference of 6.3·10⁻⁸ on ‖y‖ = 2.05·10⁻⁸, i.e. 1.3·10⁻¹⁵ absolute:
  rounding, not the second-order term.

The Gaussian kernel's "(iii)'s spectrum to 4.6·10⁻⁹" is real: it is the cross-class terms of the Gaussian kernel, which (iii) excludes.
*Fix:* "exact (to the finite-difference error 3·10⁻¹¹)", "(to the monodromy's finite-difference error)", and "at R = 1 the first order is
exact to rounding".

**D3u-8 (cosmetic; (i) l. 3943–3944 "its mean-field form is the Wasserstein gradient flow of L̄(μ) + (σ₀/2)∫∫W dμdμ — an aggregation
equation with the fit as the external field"; the proof l. 3997–3999; the §0 row; the §8 bullet l. 6873–6874; the audit bullet l. 486–487).**
- *The fit is not an external field.* L̄(μ) = ½⟨f²⟩ − ∫⟨fΦ̃(·;ϑ)⟩dμ + ½∫∫⟨Φ̃(·;ϑ)Φ̃(·;ϑ′)⟩dμdμ is quadratic in μ. It is an external
  potential V(ϑ) = −⟨fΦ̃(·;ϑ)⟩_ρ plus a second interaction K(ϑ, ϑ′) = ⟨Φ̃(·;ϑ)Φ̃(·;ϑ′)⟩_ρ, which is positive semi-definite and not
  translation invariant. The mean-field network equation is that of the standard "potential plus interaction" form.
- *The metric.* The Wasserstein metric is the one of the cost ‖ϑ − ϑ′‖²_{G⁻¹}.
- *What holds.* The velocity field v[μ] = −G∇(δℱ/δμ) is right. At the particles of μ_Θ it equals the drift exactly (0 to rounding,
  `check_i.py`).

*Fix:* "the Wasserstein gradient flow (in the metric G⁻¹) of L̄(μ) + (σ₀/2)∫∫W dμdμ — an aggregation equation in the external potential
−⟨fΦ̃⟩_ρ with two interactions, the fit's ⟨Φ̃Φ̃′⟩_ρ and the cohesion W".

**D3u-9 (cosmetic; (vi) l. 3991–3993 "their stationary spread is τ/(σ₀w_k) per whitened coordinate, against τ/σ₀ under a leak toward a
prior"; *Numerics* (E) l. 4066–4070 "against τ/(σ₀w_k)'s 0.065"; the §0 row; the §7 row (E); RESULTS l. 831–832).**
- *Variance against rms.* τ/(σ₀w_k) is the variance per whitened coordinate (4.24·10⁻³ at w_k = 1/3). The numbers compared with it are
  rms spreads, and 0.065 is √(τ/(σ₀w_k)). The md prints √(·) correctly.
- *The (1 − 1/n_k) factor.* The offsets from the class mean carry the forcing (1 − 1/n_k)·2τ, so the linear prediction is
  √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635 at n_k = 20.
- *The run* (`check_E.py`, seed 0, compact) reproduces the raw record to all digits: 0.06339, 0.07164, 0.07133 at s = 2000, 10⁴, 2·10⁴.
  Averaged over the second half, the spread is 0.0681 (compact) and 0.0638 (Gaussian, raw). The Gaussian kernel matches the linear
  prediction to 0.5%, closer than the snapshots 0.066–0.071 suggest. The compact kernel is 7% above, which is the softening the text
  names.

*Fix:* "their stationary variance is (1 − 1/n_k)τ/(σ₀w_k) per whitened coordinate (rms 0.0635 here)"; compare time averages.

**D3u-10 (minor; the audit bullet l. 485–497; [S60] l. 80; the §0 row's last column; PAPERS l. 226–229).** The audit is fair on the
ingredients it lists. It omits the closest precedents for three of the claims it calls new:
- *A particle dynamics attracted to kernel-localised means.* Bungert, Roith and Wacker, "Polarized consensus-based dynamics for
  optimization and sampling" (Math. Program. 211:125–155, 2025; arXiv:2211.05238). Each particle drifts toward its own kernel-weighted
  local mean; the authors describe it as a bounded-confidence model with a common objective. Clustering of particles under softmax-weighted
  attraction is also analysed in Geshkovski, Letrouit, Polyanskiy and Rigollet, "The emergence of clusters in self-attention dynamics"
  (NeurIPS 2023). Both are the mean-shift attraction of (i) as an interacting-particle dynamics; the cohesion leak adds the kernel-mass
  factor and the gradient-flow structure.
- *The twins' gauge curve (v).* Fukumizu and Amari, "Local minima and plateaus in hierarchical structures of multilayer perceptrons"
  (Neural Networks 13(3):317–327, 2000). Duplicating a hidden unit turns a critical point of the smaller network into lines of critical
  points of the larger one, with segments of local minima and saddles. These are the twins' curves, and the question (v) addresses is
  where the attraction lifts their degeneracy.
- *A prior-free pairwise attraction that selects the number of distinct parameters, with a strength that sets the resolution.*
  - Convex clustering: Hocking, Joulin, Bach and Vert, "Clusterpath" (ICML 2011, 745–752), with Gaussian-kernel weights (as Chi and
    Lange, JCGS 2015, report).
  - The same in identification: Ohlsson, Ljung and Boyd, "Segmentation of ARX-models using sum-of-norms regularization" (Automatica
    46(6):1107–1111, 2010), one regularisation constant trading the fit against the number of segments; Lindsten, Ohlsson and Ljung,
    "Just relax and come clustering! A convexification of k-means clustering" (LiTH-ISY-R-2992, 2011).
  - Deterministic annealing: Rose, Proc. IEEE 86(11):2210–2239, 1998. A cluster splits when the annealing parameter crosses a critical
    value set by the top eigenvalue of the cluster's covariance. That is the structure of (iv)'s "σ₀ sets the resolution", and the
    annealed R of Part C is a scale-space schedule.
- *A leak on the unexcited directions toward a prior* (the D‴.23 analogue rather than this one's): "stabilized forgetting" (Kulhavý and
  Zarrop, "On a general concept of forgetting", UMIST report 757, 1992; J. Milek, ETH dissertation, 1995).

The verdict "classical ingredients" stands. "New" should be restricted to the combination in mean-field adaptive control: the exact
blocks at clustered configurations (with the mass factor σ₀w_k), the gauge-scaled margin σ₀w_k − (q_k/w_k)λ⁺_k, and the real loop's
exponents.

*What was checked.* By web search: the bibliographic details above; CPW's report number EDI-INF-RR-0185 (December 2003) and its
equilateral-triangle result; Lian–Lee–Sudhoff–Żak, IEEE TNN 19(3):460–474, 2008; Milek's dissertation. The contents of Bungert et al.,
Fukumizu–Amari, Hocking et al. and Rose are from abstracts and citing texts, not full texts. Rose's critical-value formula was not seen,
the SSP 2011 version of Lindsten et al. was not found, and what stabilized forgetting does is from memory. The other [S60] entries (Fukunaga–Hostetler, Cheng, Carrillo–McCann–Villani,
Bertozzi–Laurent–Rosado, Hegselmann–Krause, Wu–Wang–Liu, Kulhavý 1987) are from memory and match it.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` rewrites the md byte-identically (md5 `c7848f61…`). The five tests pass in 4.42 s (2.83 s
   of it the real loop's period map).
2. **(i)** (`check_i.py`, 9 particles with pairs inside the range).
   - My ∇Φ̃ and ∇²Φ̃ against differences: 1.4·10⁻¹⁰ and 9.2·10⁻¹¹ (relative). ω′ = ψ/2 on [0, 1.5] to 8·10⁻¹⁰ and ω(0) = 0, for both
     kernels. ∇²W(0) = G⁻¹ to 6·10⁻⁷ (a central second difference, step 10⁻⁴). ψ = (1 − z)³₊ is C² and not C³ at z = 1.
   - The drift equals −G∇J^coh_m to 6.9·10⁻¹⁰ and 8.1·10⁻¹⁰ (relative; compact R = 1.2, Gaussian R = 0.8). The repository's convention
     is the same: J = mL̄ + (σ₀m/2)(1/m²)ΣW, f̂ = (1/m)ΣΦ̃, the per-particle law without 1/m.
   - The forces sum to ≤ 5·10⁻²⁰. They equal σ₀ψ̄_i(m_i − ϑ_i) to 6·10⁻¹⁹.
   - The mean-field velocity at the particles equals the drift (0 and 7·10⁻²¹). The barycentre moves with the mean learning drift alone
     (2.8·10⁻¹⁷).
   - The pair force u ↦ ψ(|u|²/R²)u is Lipschitz with constant 1 (whitened, both kernels). The compact one is maximal at |u| = R/√7
     = 0.378R (0.238R) and decreases to 0 at R.
   - My drift and J equal the repository's to 4.7·10⁻¹⁶ (cross-check).
3. **(ii)–(iii) at the 15 count vectors of m = 15** (`check_A1.py`, analytic Jacobian).
   - The drift is ≤ 3.6·10⁻¹⁷ (compact) and ≤ 6.2·10⁻¹⁰ (Gaussian); the md has 2.79·10⁻¹⁷ and 6.23·10⁻¹⁰.
   - The spectrum is −σ₀w_k ∪ −eig H_L(w) to 5.6·10⁻¹⁶ (compact); the Gaussian kernel's deviation is 4.56·10⁻⁹.
   - At r = 0, Ã_w = H_L(w) to 6·10⁻¹⁵. The leak's Jacobian vanishes on the class means (≤ 6·10⁻²⁰), and the mixing blocks are ≤ 10⁻¹⁷.
   - The smallest inter-atom distance is 5.8035. h₁(w) runs from 2.232·10⁻⁴ to 4.185·10⁻⁴, σ₀w_k from 1.415·10⁻³ to 3.773·10⁻³. At
     every count vector h₁ < σ₀ min_k w_k: the slowest mode is collective, the leak absent there.
   - The kernel of Df̂ is the separations exactly where H_L(w) ≻ 0: on the instance at all 15 count vectors.
4. **The Gaussian kernel's bias** (`check_A2.py`).
   - (5, 5, 5): ‖y‖ = 2.0513·10⁻⁸, 1.2520·10⁻³, 6.8279·10⁻², 0.86024 at R = 1, 1.5, 2, 3. The first order Ã⁻¹b differs by 6.3·10⁻⁸,
     9.06·10⁻⁴, 1.02·10⁻², 0.564; the fit errors are 6.59·10⁻¹⁰, 3.93·10⁻⁵, 2.13·10⁻³, 3.09·10⁻².
   - (6, 5, 4): 1.1354·10⁻⁸ … 0.84177, the relative differences 8.4·10⁻⁸, 7.59·10⁻⁴, 8.42·10⁻³, 0.621.
   - All of these are the md's values; the R = 1 differences are rounding (D3u-7).
   - With b's derivative, (Ã − Db)⁻¹b is closer at R = 1.5 (3.8·10⁻⁴) and farther at R = 2 (1.9·10⁻²), where the O(|b|²) terms take
     over.
   - The stated order σ₀‖Ã⁻¹‖max_kΣ_lw_ld_kle^{−d²/(2R²)} is 1.57·10⁻⁷, 1.04·10⁻², 0.551 at (5, 5, 5): 7.6–8.3 times the actual bias,
     an upper order.
   - The Gaussian stationary points are stable: the particle Jacobian's top eigenvalue is −3.06·10⁻⁴ … −5.26·10⁻⁴.
5. **(iv) and Part A3** (`check_A3.py`, `check_A3b.py`, `cross_bq.py`).
   - My ϑ^b(q) has gradient 2.2·10⁻¹⁷ and fit error 0.016953. It is 5.7·10⁻¹¹ from the repository's (cross-check). Ã ≻ 0 there: its
     eigenvalues run from 2.43·10⁻⁴ to 0.221, so ϑ^b(q) is a strict local minimum.
   - λ⁺_k = 9.9999·10⁻⁵, 2.8607·10⁻⁴, 1.0845·10⁻⁴. The gauge scaling R_k(ϑ^b(w)) = (q_k/w_k)R_k(ϑ^b(q)) holds to 9·10⁻¹⁴ (relative),
     the residual unchanged to 3·10⁻¹⁶.
   - The margins σw_k − (q_k/w_k)λ⁺_k equal the Jacobian's top separation eigenvalues to ≤ 1.5·10⁻¹⁷ at the three count vectors and
     both σ: the md's table to every printed digit. The first split at w = q is at σ = 3λ⁺₂ = 0.1213σ₀.
   - *The mass factor.* At σ₀ the cohesion leak's splitting thresholds are w_k < √(q_kλ⁺_k/σ₀) = 0.069, 0.116, 0.071, against
     w_k < q_kλ⁺_k/σ₀ = 0.0047, 0.013, 0.0051 under the adaptive prior (D‴.28(iii)). At σ₀/8 they are 0.194, 0.328, 0.202. A light class
     is held by σ₀w_k and splits far more easily, which is what (iv)'s w_k² says.
   - *The split runs*, reproduced from `default_rng(0)` (the three draws in order), match the md's diameters to all printed digits:
     (2.19·10⁻³, 1.37·10⁻³, 1.88·10⁻³) → (4.39·10⁻¹², 3.40·10⁻⁴, 1.22·10⁻¹¹) at (5, 5, 5) and σ₀/8, and 2.19·10⁻³ → 1.69·10⁻² →
     0.277 → 1.14 → 1.17 → 1.02 for the middle class at (6, 4, 5) and σ₀/8. The end clusters (0.400, 0.200, 0.067, 0.333) and the fit
     error 0.0148 match too.
   - *The predicted rate.* The middle class's separation component along the top eigenvector of its block grows at 1.224·10⁻⁴ on
     s ∈ [2·10³, 10⁴], against the margin's 1.218·10⁻⁴. It then accelerates (1.25, 1.33, 1.40·10⁻⁴ per 10⁴) as the compact force softens.
     At (5, 5, 5) the decay is 8.76·10⁻⁶ against 8.70·10⁻⁶, and at σ₀ it is 1.529·10⁻³ against 1.5289·10⁻³.
   - *What the split achieves.* The shed particle goes to c = −0.31 with a = −3.32, the other side of the middle bump from the missing
     one (0.8). The fit improves only from 0.0170 to 0.0148.
6. **(v)**: the md's twin and dead-class runs reproduce (`check_A4.py`). The twins: 2.31 → 6.1·10⁻³ at s = 4000, 0 by 2·10⁴, the
   fit error 4.45·10⁻¹⁰; compact 2.31 throughout. The dead class: Gaussian R = 2 to rounding by 1.6·10⁴, the fit error 1.64·10⁻³;
   compact 5.19 throughout. The counterexamples are under D3u-1 and D3u-3.
7. **(vi) and Part D** (`check_D.py`, `check_D2.py`).
   - The lift returns after a period to 7.81·10⁻¹⁶ (the md's 7.81·10⁻¹⁶).
   - The variational generator at the lift has the exact block structure: |QAP| = 8.4·10⁻¹⁷, |PAQ| = 3.5·10⁻¹⁷, the separations' block
     −σ₀w_k/k_e to 4·10⁻²⁰. ∂_xf̂ = ∂_xf on the orbit to 2.8·10⁻¹⁰ (differences).
   - The 46-dimensional monodromy (matrix ODE): its 45 slow rates equal σ₀w_k/k_e (36) and the 10-dimensional collapsed system's
     9 to 5.2·10⁻¹¹.
   - The collapsed rates times k_e are 2.658·10⁻⁴ … 0.1789, against h_j(w) = 3.080·10⁻⁴ … 0.1834: −13.7% for the weakest mode at
     k_e = 25. This is the c_j/k_e² correction of (vi), not an error.
   - D2 (m = 15, cohesion, seed 0) reproduces: the rms 0.00487, 0.00296, 0.00134, 5.8·10⁻⁴ near t = 0, 500, 1500, 2900, the fit error
     0.0218, three clusters.
8. **Part E** (`check_E.py`): seed 0 under the compact kernel reproduces the raw record to all digits (D3u-9).
9. **Part C** (`broad_run.py`, from the seeds, my integrator). Three runs of the annealed Gaussian schedule reproduce the md exactly:
   - m = 60, seed 0: three classes, the fit error 1.925·10⁻⁷, masses 0.500, 0.283, 0.217, |a|/a_max 0.882;
   - m = 60, seed 4: 2.251·10⁻⁶, masses 0.367, 0.333, 0.300;
   - m = 120, seed 2: 1.815·10⁻⁷, masses 0.392, 0.383, 0.225.
10. **Part H** (five bumps, seed 0): six classes after R = 2 (fit error 0.120) and at the end (0.1179). The masses 0.417, 0.183, 0.150,
    0.117, 0.067, 0.067 and the saturated inner classes are the md's.
11. **Part G** (`g_track.py`, from `default_rng(700 + seed)`, all three seeds with the repository's ϑ^b(q); seed 0's first five
    stages also with mine, identical to the printed digits). Every stage's clusters, masses, margins, end clusters and fit errors are
    the raw record's; the fit errors agree to ≤ 3.5·10⁻⁴ (relative, the printed digits). The tally is under D3u-6.
12. **The tests** (`check_tests.py`, `tests_values.py`): the asserted quantities and the mutants of D3u-5.
13. **Bibliography** (web search): see D3u-10.

## Not checked

- **Part F (the cogging servo).** Not rerun; its numbers were read from the md and the raw record only.
- **The other leaks of Parts B–E.** The fixed and adaptive priors, the σ-modification and no leak were not rerun. Neither were Part B's
  seeds 0–2 nor the m = 60 runs of D2. Their values come from the raw record. The reproduced runs (B seed 3, C three runs, D2 m = 15,
  E seed 0, G, H five bumps seed 0) match it to the printed digits, which makes the rest credible.
- **Parts H two and four bumps.** Not rerun.
- **The noise statements of (vi).** They are D‴.30's, formal there, and were not checked. Neither was the separations' noise threshold
  "scaled by w_k".
- **The propagation of chaos and Łojasiewicz's theorem.** Checked only as hypotheses: the field's Lipschitz constant, and the observation
  that solutions can be unbounded (D3u-1's run toward the saturation).
- **c_j(w) of (vi).** Only the collapsed rates at k_e = 25 were compared with h_j(w); the expansion in 1/k_e was not.
- **A stable twin configuration.** The counterexample of D3u-1 is a saddle that the flow leaves toward the saturation, with four
  classes. I did not find a configuration with twins that is a strict local minimum of J^coh_m: the fit is only quartic in the
  twins' split directions, and in my example the attraction of a class beyond R destabilised them. Whether one exists (e.g. with a third class within R)
  is open.
- **The literature.** Full texts were not read (D3u-10).
- **The concurrent session.** Its running experiments (`exp_eiv_remedy.py --only-D2`) and any uncommitted edits are outside this report.
