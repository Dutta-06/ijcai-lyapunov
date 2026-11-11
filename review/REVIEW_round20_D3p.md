# Referee report D3p (round 20) — PARADIGM §5d: Proposition D‴.29 ((i)–(iv), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.29 |", the §8 bullet, RESULTS (Proposition D‴.29), the LITERATURE audit bullet, PAPERS.md paper 11, `src/core/sliding.py`, `tests/test_sliding.py`, `experiments/exp_servo_cogging.py` → `results/servo/servo_cogging.md`; and Proposition D‴.27 as revised in round 19 (compact 𝒵 with Z(0) ∈ 𝒵, the sharper bound, the η = ∞ convergence, the conditional (iv), Status), the projected step of `AdaptivePriorLoop.run`, Part H of `exp_adaptive_prior.py` → `results/chaos/adaptive_prior.md`, test 10 of `tests/test_adaptive_prior.py`

**Setup.** Branch `claude/festive-bardeen-eopej2`, material reviewed at ee2ed27 (`git log --oneline -5` at the start: ee2ed27, 442cd6c,
d2092a5, 4e2e938, b76b3b9).
- The tree.
  - At the start (22:13 UTC) the working tree had no change to tracked files and one untracked file (`src/core/measurement_noise.py`).
  - By 23:20 UTC a concurrent session had modified seven tracked files and added seven untracked ones. After a container restart
    (23:27 UTC) the branch had gained 3609f86 (Proposition D‴.30), and the tree carried uncommitted manuscript edits and D‴.31 files.
  - In 3609f86 every passage reviewed here is textually unchanged; D‴.30 is inserted after D‴.29's *Numerics*. I checked the D‴.27 and
    D‴.29 blocks, both §0 and §7 rows, and both §8 bullets.
  - None of the concurrent material was used or reviewed.
- How the review was run.
  - The material was read from ee2ed27 extracted with `git archive` into `ref_D3p/tree/`. Every script imports that tree, never the
    repository; its key files were re-verified by md5 against `git show ee2ed27:` after the restart.
  - A second copy, `ref_D3p/regen/`, was used for `--from-raw`.
  - No repository file was modified, and no git command touching the index or the tree was run.
  - Python ran with `PYTHONDONTWRITEBYTECODE=1`, one BLAS thread and one job at a time; pytest ran with `-p no:cacheprovider`.
  - Every script had finished before the restart, and its output was intact.

Line numbers are ee2ed27's:
- `PARADIGM.md`:
  - Proposition D‴.29 l. 3317–3398: (i) l. 3325, (ii) l. 3328 (the exact-filter sentence l. 3334), (iii) l. 3337, (iv) l. 3341,
    Proof l. 3348 (the Theorem J sentence l. 3350), *Status* l. 3358, *Numerics* l. 3367 (A at l. 3369, the misplaced reading
    l. 3381–3384, C at l. 3385–3391).
  - Proposition D‴.27 l. 3161–3240: (i) l. 3168, (ii) l. 3173 (the "reduced form" l. 3174), (iii) l. 3181, (iv) l. 3185, Proof l. 3190
    (the D‴.17 reduction l. 3196–3203), *Status* l. 3209, *Numerics* l. 3223 ((b) l. 3232, (c) l. 3234).
  - The rows and bullets: §0 rows l. 128 (D‴.27) and l. 130 (D‴.29); §7 rows l. 5821 and l. 5823; §8 bullets l. 6073–6077 (D‴.28),
    l. 6078–6083 (D‴.29), l. 6084–6089 (D‴.27).
  - Theorem J l. 4964.
- `results/RESULTS_PARADIGM.md`: D‴.27 l. 640–660, D‴.29 l. 684–701.
- `LITERATURE_TRANSPORT.md`: D‴.27 l. 400–409, D‴.29 l. 417–424.
- `PAPERS.md`: paper 11 l. 151–162.

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3p/scripts/`;
each `.out` is one level up, in `ref_D3p/`.

| script | what it checks |
|---|---|
| `from_raw.out`, `tests.out`, `collect.out` | `--from-raw` of both experiments; the 13 tests of `test_sliding.py` and `test_adaptive_prior.py`; the collection count |
| `common.py`, `inspect_raw.py` | shared setup (the archived tree); the servo raw file's structure |
| `mymodel.py`, `check_mymodel.py` | an independent model of the servo instance: the 1-D saturated network, the orbit grid, analytic first and second derivatives, the weighted loss, its Hessian, the residual curvature, Newton; checked against finite differences and the repository |
| `check_partA.py` | Part A independently: the approximant, Ã's spectrum, the rates, λ⁺_k, W_g, the gauge transform, the floor law from finite differences of Newton solutions in w, the floors against the law, the misplaced prior's stationary point and D‴.5 margin |
| `check_best.py`, `check_best_fine.py`, `check_kflow.py`, `check_runoff.py` | is ϑ^b(q) the best approximant? (120 BFGS starts; a 10× finer grid; the loss along the saturation); the K₀-atom adaptive flow from the misplaced atoms and from 40 random starts to s = 10⁵ (four of them to 10⁸) |
| `check_partB.py` | Part B's table from the raw records with my own code; the ratios quoted; two jobs re-run with `run_job` |
| `check_coldstart.py` | Part B's designs from the engineer's model (clouds around the misplaced atoms); the collapsed joint flows to s = 10⁵ |
| `check_noleak.py` | the σ₀ = 0 baseline from Part B's warm clouds (three seeds, m = 16, 64, 128, 256) |
| `myservo.py`, `check_sliding_code.py` | an independent real servo; `src/core/sliding.py` against independent implementations (gradient, drift, residual curvature, Gauss–Newton block, `SlidingRealLoop.rhs`, `SlidingServo.run` against DOP853, `filtered_error` against the closed form) |
| `check_averaging.py`, `check_decompose.py` | D‴.29(ii): the real servo against the averaged loop at k_s = 25, 50, 100, 200; the position-error predictions, and the sources of their error (nonlinearity, the within-period motion of Θ) |
| `check_partC.py`, `check_partC_more.py` | Part C's table from the raw records; one job re-run with `real_job`; dense against 0.2-spaced records; the η = σ₀ prior frozen; the residual's spectrum; the parameter drift over 1000 time units |
| `check_vonmises.py`, `check_vonmises_abs.py`, `check_vonmises_centred.py` | D‴.29(iv): the fibre; the Gauss–Newton zeros on the circle, on ellipses and on a two-harmonic orbit; centred kernels; test 2's discrimination; the full Hessian at a fitted point |
| `apold.py`, `check_projected_step.py` | the new projected step against an event-driven DOP853 reference at six steps, beside the 442cd6c step |
| `check_test10.py` | test 10's quantities at four steps; test 10 with the 442cd6c step patched in |
| `partH_rerun.py`, `partH_stored.py` | Part H re-run with the experiment's `projected_job()` against the stored entry |
| `check_boxscan.py` | Part H(c) extended to amplitude bounds 0.5–0.65 of a_max |
| `search_selection.py`, `check_selection.py`, `check_selection_limit.py` | D‴.27(ii) at η = ∞: Lloyd-consistent non-optimal plans; the reduced energy along the code's flow; a limit from a consistent plan |

The long runs were not rerun in full: 25 of Part B's 27 jobs and 17 of Part C's 18. Every number below comes from the scripts listed.

## Verdicts

### (A) Proposition D‴.29

| object | verdict |
|---|---|
| **(i)** the error channel | **holds exactly.** Re-derived: ë = f(z) + u − ẍ_m = r_Θ(z) − λė − k_s s. So ṡ = −k_s s + r_Θ(z), with r evaluated at the plant's state z = z_m + (e, ė), not at z_m, and ė = −λe + s. The controller needs ẍ_m (x_m ∈ C³ is assumed). `SlidingRealLoop.rhs` satisfies the identity to 1.8·10⁻¹⁵. The law ϑ̇ = G∇(sΦ̃(z;ϑ)) − (σ₀/k_s)(ϑ − T) is Theorem J's (CL) with s for e and k_s for k_e; in τ = t/k_s both give G⟨r∇Φ̃⟩ − σ₀(ϑ − T). |
| **(ii)** the averaged slow loop | **holds numerically; "formal" is the right label, but "Theorem J's argument transferred" overstates (D3p-5).** Independent real servo against an independent averaged loop, same slow horizon (`check_averaging.py`): the relative displacement error is 1.51·10⁻², 4.61·10⁻³, 1.55·10⁻³, 6.31·10⁻⁴ at k_s = 25, 50, 100, 200, i.e. O(1/k_s). s is slaved at the rate k_s. e is not: it is the (d/dt + λ)⁻¹-filtered response, relaxing at the fixed rate λ. That is still fast against the parameters, the separation being λk_s in slow time. The exact-filter claim holds only for small high-frequency residuals (D3p-6). |
| **(iii)** the transfer | **holds for the averaged-loop statements**, which use ρ, the target and the network: the definability of D‴.17 holds for the period average of an analytic orbit, and the constants of Theorem D‴ and D‴.21 are recomputable. The list includes D‴.24(vi), a real-loop statement, and it carries D‴.27(ii)'s flawed η = ∞ convergence into the second-order setting (D3p-5, D3p-3). |
| **(iv)** where the kernels live | **exact as computed.** Along the von Mises fibre f̂ on the orbit changes by ≤ 2.8·10⁻¹⁷. On test 2's circle 5 eigenvalues of 20 are machine zeros (≤ 6.8·10⁻¹⁷ of the largest), and the next is 2.4·10⁻¹⁰. At a fitted critical point the full Hessian has the same zeros (4.2·10⁻¹³–1.4·10⁻¹², the finite-difference level). The fixed/adaptive/projected readings follow (h = 0: rates σ₀, λ⁻ = 0, bounded). But the exactness needs ω₀ = 1 in isotropic units; the weak excitation is generic on other orbits; and "exactly one per kernel" holds only in general position (D3p-7). |
| Proof | (i), (iv) correct; (ii) a sketch whose key sentence is inaccurate (D3p-5); (iii) one line. |
| *Status* | the labels are right; the Theorem J wording, the (iii) exclusion and the novelty need changes (D3p-5, D3p-7). |
| *Numerics* (A) | **every number reproduced by an independent model, but ϑ^b(q) is not the best approximant: the loss's infimum lies at the saturation (D3p-1, major).** |
| *Numerics* (B) | **every md entry recomputed from the raw records, two jobs re-run bitwise.** The misplaced experiment is vacuous at η = ∞ because the network starts at the oracle (D3p-2, major); the classical no-leak baseline is missing (D3p-4); readings (D3p-9). |
| *Numerics* (C) | **every ratio recomputed, one job re-run bitwise** (136 s). The η = σ₀ prior is frozen; 0.2-spaced records; the exact filter's advantage is state-dependent (D3p-5, D3p-6). |
| §0, §7 rows; §8 bullet; RESULTS; LITERATURE; PAPERS 11 | **numbers hold.** The "1.46% / Proposition D‴.28's regime" and the §8 friction contrast (D3p-1); "removed by the adaptive one" (D3p-2); the baselines and the cogging framing (D3p-4); the novelty (D3p-6(d), D3p-7). |
| `src/core/sliding.py` | **correct against independent implementations** (`check_sliding_code.py`): the phase-plane gradient within 3.2·10⁻¹⁰ of finite differences of the forward map; the drift 6.1·10⁻¹²; the residual curvature 3.1·10⁻⁹ of second differences; the Gauss–Newton block 1.5·10⁻¹¹; `SlidingServo.run` (RK4, dt = 2·10⁻³) within 1.6·10⁻⁷ of DOP853 after four periods (displacement 3.3·10⁻³); `filtered_error` within 1.7·10⁻¹⁷ (s) and 2.6·10⁻¹⁸ (e) of the closed form per harmonic. |
| `tests/test_sliding.py` | **passes; test 2 discriminates** the circle: it fails on the ellipse with the velocity axis ×1.001. Test 3 is consistent with O(1/k_s); mine gives 1.5% at k_s = 25 against its 3% bound. |
| `exp_servo_cogging.py`, the md | **reproduce** (`--from-raw` byte-identical). Docstrings: "best" (D3p-1); "e² every 0.1 time units" is 0.2 (D3p-6); `real_job` re-implements `SlidingServo` and freezes the η = σ₀ prior (D3p-5). |

### (B) Proposition D‴.27 as revised in round 19

| object | verdict |
|---|---|
| **(i)** the bound | **holds.** 𝒵 compact and Z(0) ∈ 𝒵^{K₀} are now stated. The sharper bound d(s) ≤ B/σ₀ + (d(0) − B/σ₀)₊e^{−σ₀s} is re-derived: comparison with u′ = B − σ₀u where d > B/σ₀. The split rows' targets lie in 𝒵 by convexity. The unbounded case is remarked. |
| **(ii)** the energy | **finite η: monotonicity holds** (Moreau's decomposition; the plan is a joint minimiser for (Θ, Z), so the vertex selection is active), and the convergence is formal as labelled. **η = ∞: "non-increasing along its reduced form min_Z𝒥" is false, and the reduction to Proposition D‴.17 is invalid (D3p-3, major).** |
| **(iii)** | **holds** as labelled; the local convergence is sketched. |
| **(iv)** | **the conditional holds up to an edge case** (D3p-10). The identification −σ₀(ϑ − Π_𝒵ϑ) is exact (round 19). |
| Proof | (i) correct; (ii) has a gap at η = ∞ (D3p-3). |
| *Status* | "the convergence proved at η = ∞" (D3p-3); the novelty paragraph is as fixed in round 19. |
| the projected step (`run` l. 133–183) | **consistent** (`check_projected_step.py`, η = σ₀, Part H(a)'s system): the prior reaches the face at s = 297.65–297.8 at all six steps (the reference: 297.6488). The error to the event-driven reference at s = 100 falls 3.54·10⁻⁶ → 8.84·10⁻⁷ → 2.21·10⁻⁷ for dt = 0.25 → 0.125 → 0.0625 (second order), and at s = 300 6.30·10⁻⁷ → 1.91·10⁻⁷ → 4.07·10⁻⁸ (order 1.7–2.2). There is no stall; at s = 20000 the gap is 0 and the drift 2.58·10⁻⁷ (the flow's own slow approach), the same at every dt. Its fixed points are the projected flow's when ηdt < 2. At η = ∞ the change is inert (the prior's velocity is 0, so Z_mid = Z), and the plan refresh is untouched. The 442cd6c step stalls at dt = 0.2, 0.125, 0.1, 0.0625 (gaps 2.40·10⁻⁴, 1.17·10⁻⁴; drift 5.94·10⁻³), as round 19 found. |
| Part H (a)–(d), the md | **re-run identical to the stored H** in every printed value (3.4372504842394846·10⁻¹² for (a)'s end states; margins 0.16907, 0.16885; the scan; the G(a) continuation) (`partH_rerun.py`, 210 s). (a) at two steps tests the shared fixed point, not consistency (D3p-10). (b): "at least 17%" is 16.9% (D3p-10). (c): "between" fails below 0.7, and the 0.7 box excludes the exact atom (D3p-8). (d) holds. Part H's end-drift assertion discriminates: the old step leaves 5.94·10⁻³ at dt = 0.125. |
| test 10 | **passes; its new assertion discriminates** (`check_test10.py`): with the 442cd6c step it fails at l. 281, z_a 1.171·10⁻⁴ inside the face at dt = 0.125. Measured with the new step at dt = 0.125, 0.25, 0.2, 0.1: on the face; max ϑ_a 6.805 (< 8); the largest 𝒥 increment between records −4.3·10⁻⁷ to −1.1·10⁻⁶; the projected drift at s = 2000 1.98·10⁻⁴. |
| §0, §7 rows; §8 bullet; RESULTS; LITERATURE | **numbers hold.** "The convergence proved at η = ∞" / "at η = ∞, every solution convergent" / "follows from Proposition D‴.17's argument" (D3p-3); "between" (D3p-8). |

### The round-19 fixes to D‴.27 (D3n-1 … D3n-7), item by item

| id | fixed as logged? | correct and complete? |
|---|---|---|
| D3n-1 | yes: compact 𝒵, Z(0) ∈ 𝒵^{K₀}, the sharper bound, the unbounded remark (statement, rows, §8, branch 7, D‴.23(i), audit) | correct; Part H(c)'s 0.7 row starts outside its box (D3p-8) |
| D3n-2 | yes: positions projected inside the step; Part H at dt = 0.25 and 0.125 with the end drift asserted; test 10's face assertion | correct — consistent, second order before the face, no stall at six steps; "agree to ten digits" is about end states (D3p-10) |
| D3n-3 | yes: (iv) conditional; "between" an observation; Part H(c), (d) | (iv): an edge case (D3p-10); "between" fails below 0.7 (D3p-8); (d) correct |
| D3n-4 | yes | correct |
| D3n-5 | yes | correct (B ≤ 26.08, L ≤ 4.658, 4.044 by search) |
| D3n-6 | yes: margins printed, "by construction", test 10 asserts "passes 11" (l. 275) | correct; the margin is 16.9%, written "at least 17%" (D3p-10) |
| D3n-7 | yes | (a)–(d) correct. **(e) — round 19's own suggestion — is wrong:** the η = ∞ flow selects Lloyd-consistent plans, not joint minimisers, so Proposition D‴.17's argument does not apply as it stands (D3p-3) |
| D3n-8 … D3n-11 | present (spot-checked by `grep`: the five "(restored by projecting the prior, Proposition D‴.27)", "every representable weight", "leads at m = 60 and 240", the branch-minimum wording, the heavy tail) | not re-verified numerically (D‴.24, outside this brief) |

## Confirmed issues

**D3p-1 (major; D‴.29 *Numerics* l. 3369 "(A) The best weighted approximant ϑ^b(q) is a nondegenerate critical point with the relative fit
error 1.46·10⁻² — no exact representation, so Proposition D‴.28 is the analysis"; the §0 row l. 130 "a cogging torque that four kernels
approximate to 1.46% (Proposition D‴.28's regime)"; the §7 row l. 5823; the §8 bullet l. 6082 "velocity-dependent friction, whose best
approximant by a few kernels lies at infinity"; RESULTS l. 688; PAPERS l. 153–156 "without an exact representation (the regime of
Proposition D‴.28)"; `exp_servo_cogging.py` l. 9 and the `approximant` docstring l. 84; the md l. 7).**

ϑ^b(q) is a nondegenerate local minimum, the one BFGS reaches from the experiment's START. It is not the best weighted approximant: the
weighted four-kernel loss's infimum lies at the saturation boundary.

*The better fits are all saturated* (`check_best.py`):
- 120 BFGS starts: 25 end below 0.0146 and none at ϑ^b(q).
- The ten best distinct end points have fit errors 0.00639–0.0101, and every one has at least one amplitude at a_max (|α| ≥ 19).
- The best: a = (0.465, 4.0, −4.0, −3.498), c = (−0.955, 0.498, 0.224, 0.781), s = (0.166, 0.735, 0.298, 0.422), fit error 0.006390.
  It is the same on a ten times finer orbit grid (`check_best_fine.py`), where ϑ^b(q) stays at 0.014620.

*The infimum is approached, not attained* (`check_best_fine.py`):
- Holding the saturated amplitudes at |α| = 6, 8, 10, 14 (|a| = 3.62, 3.86, 3.95, 3.99) and re-optimising the rest gives the fit errors
  0.00714, 0.00665, 0.00649, 0.00640.
- So even with |α| ≤ 6 four kernels fit the torque to 0.71%, half of ϑ^b(q)'s error.

*The adaptive prior's dynamics run off toward it:*
- The collapsed K₀-atom flow at η = ∞ is the unregularised weighted gradient flow. From 4 of 40 random starts (`check_runoff.py`) it has
  fit errors 0.0084–0.0109 at s = 10⁵.
- max|α| then grows 2.5–3.1 → 3.3–3.5 → 4.6 → 7.0 over s = 10⁵, 10⁶, 10⁷, 10⁸, while the fit error falls to 0.0075 — the logarithmic run-off of
  Part G(b).
- None of the 40 starts reaches ϑ^b(q) by s = 10⁵ (`check_kflow.py`).

So, globally, the servo is the case that the D‴.28 bullet lists as open (l. 6076: "critical points at infinity (the best approximant on the
saturation boundary, as for … a cogging function that two kernels fit as a dipole), where only the projected prior keeps the parameters
bounded"). It is not the contrast that l. 6082 draws with friction. Part A's analysis at ϑ^b(q) is correct (see the confirmations), but it
and Parts B–C, all started at ϑ^b(q), are statements about one basin.

*Fix:*
- Replace "the best weighted approximant" by "a nondegenerate local minimum ϑ^b(q) (reached by BFGS from START)" in the paragraph, the md
  and the docstrings.
- Replace "four kernels approximate to 1.46%" by "whose local minimum ϑ^b(q) has the fit error 1.46%; the loss's infimum, ≈ 0.64%, lies at
  the saturation".
- State that globally the servo is Proposition D‴.27's case. Either run the projected adaptive prior from cold starts, or restrict the
  claims to the basin of ϑ^b(q).
- Drop the friction/cogging contrast at l. 6082.

**D3p-2 (major; D‴.29 *Numerics* (B) l. 3381–3384 "while the adaptive prior at η = ∞ is unaffected (… the runs are the matching ones) and
at η = σ₀ removes the bias more slowly"; the §0 row "a misplaced prior's bias kept by the fixed designs … and removed by the adaptive
one"; the §8 bullet l. 6080–6081 "the adaptive prior removes the fixed prior's floor and a misplaced prior's bias"; RESULTS l. 684–701;
PAPERS paper 11 "Claim" l. 153–156 and "Missing" l. 161).**

Every design of Part B starts from clouds of width 0.02 around ϑ^b(q). ϑ^b(q) is computed from the true torque, and the "matching prior" is
ϑ^b(q) itself. So the application starts at the answer.

(a) *At η = ∞ the misplaced-prior experiment tests nothing.*
- The first plan's class means replace the prior.
- The misplaced rows are the matching rows: the stored f̂ records are identical for every seed and m (`check_partB.py`).
- "Removed by the adaptive one" is therefore not demonstrated at η = ∞.

(b) *From the engineer's model the bias decays only slowly.* `check_coldstart.py` uses the experiment's own `design_loop`,
`mixture_start` and `run`, with clouds of the same width around the misplaced atoms, m = 64, seed 0.

| design | fit error at s = 600, 1500, 2000 |
|---|---|
| adaptive, η = ∞ | 0.0328, 0.0252, 0.0230 |
| adaptive, η = σ₀ | 0.0341, 0.0323, 0.0294 |
| nearest atom, fixed | 0.0585, 0.0579, 0.0571 |
| OT, fixed | 0.0498, 0.0484, 0.0483 |

- From the oracle start the adaptive designs reach 0.0148 (η = ∞) and 0.0159 (η = σ₀).
- The collapsed flows from the misplaced atoms give the time scale:
  - η = ∞: 0.0297, 0.0213, 0.0169, 0.0155 at s = 600, 2000, 10⁴, 10⁵, still 0.553 from ϑ^b(q) in the parameters;
  - η = σ₀: 0.0345, 0.0259, 0.0181, 0.0160.
- That is about 10⁵ slow units, i.e. about 2.5·10⁶ time units of the real servo at k_s = 25. By D3p-1, the flow need not end at ϑ^b(q).

(c) *What holds.* From the oracle start at η = σ₀, "removes the bias more slowly" is confirmed on a longer horizon: the K₀-atom joint flow
has the fit errors 0.0240, 0.0183, 0.0154, 0.0149 at s = 600, 2000, 10⁴, 10⁵.

*Fix:*
- In the paragraph, the §0 row, §8, RESULTS and paper 11, say that the network starts at ϑ^b(q) and that the matching prior is ϑ^b(q).
- Replace "removed by the adaptive one" by "at η = ∞ the misplaced prior never acts (the runs are the matching ones); at η = σ₀ its bias
  decays (0.035, 0.020, 0.019 at s = 2000)".
- Add the engineer's-model start to Part B and to paper 11's "Missing".

**D3p-3 (major; D‴.27(ii) l. 3173–3176 "𝒥 … is non-increasing along the projected flow (at η = ∞ along its reduced form,
min_{Z∈𝒵^{K₀}}𝒥) … At η = ∞ every solution converges to a Clarke-critical point of the reduced energy"; Proof l. 3196–3203 "the particles
follow the vertex-selection subgradient flow of a minimum of finitely many C^{1,1} definable functions, so Proposition D‴.17(i)–(ii) apply
as they stand"; *Status* l. 3209–3210 "the convergence proved at η = ∞ (Proposition D‴.17's argument)"; the §0 row l. 128; the §8 bullet
l. 6085–6086; LITERATURE l. 406–407 — the round-19 fix D3n-7(e)).**

At η = ∞ the flow's plan γ is an optimal plan to ν_Z with Z = Π_𝒵T̂(γ): a Lloyd-consistent pair. The code reaches it by alternation at each
refresh (`run` l. 155–160, l. 176). Such a γ is optimal for its own projected barycentres. It need not minimise
E_γ(Θ) := min_{Z∈𝒵^{K₀}}𝒥 over the plans. When it does not, the particles follow −∇E_γ of a piece that is not active in
min_γE_γ = min_Z𝒥. That is not a subgradient trajectory of the reduced energy, and Proposition D‴.17 — whose fixed prior makes every
selected piece active — does not apply as it stands.

*Evidence:*
- `search_selection.py`, pure geometry: of 20000 random 6-point configurations in the plane (K₀ = 3, masses 1/3), 7455 have a
  Lloyd-consistent, non-optimal pairing along whose pure-leak flow the minimal piece increases.
- `check_selection.py`, the repository's code on one of them: the instance's network, m = 6 particles at zero amplitude, f ≡ 0 (so the
  learning force vanishes), Part H's box, `AdaptivePriorLoop(plan='ot', η = ∞)`, the prior at the class means of the consistent pairing.
  - The code keeps that pairing until s ≈ 100.
  - Its own energy 𝒥(Θ, Z) falls from 4.435·10⁻³ (largest increment −2.35·10⁻⁶).
  - The reduced energy min_Z𝒥, enumerated over the 15 balanced pairings, rises from 1.167·10⁻³ at s = 0 to 1.511·10⁻³ at s = 70 (+30%).
  - With a small three-bump target, so that the learning term acts, it rises from 1.251·10⁻³ to 1.574·10⁻³.

So "non-increasing along its reduced form" is false, and the convergence "proved by Proposition D‴.17's argument" is not proved.

*What survives:*
- The selected energy 𝒥(Θ(s), Z(s)) is non-increasing, also at plan switches: the new plan is optimal for the old prior, and the prior then
  re-minimises.
- A Kurdyka–Łojasiewicz argument should give bounded solutions finite length. It needs one desingularising function for the finitely many
  definable pieces that share the limiting value at a limit point.
- The limit is then a critical point of a Lloyd-consistent piece, which need not be Clarke-critical for min_Z𝒥.
- I found no case where it is not: in `check_selection_limit.py` a consistent "wide" pairing became the active one by s = 2000, and the flow
  ended at a replicated two-kernel fit.

The same identification underlies Theorem D‴.23(iv)'s formal "gradient flow of L̄ + (σ₀/2)𝒬_q": as η → ∞ the prior tracks a local, not a
global, minimiser of W²(μ_Θ, ·).

*Fix:*
- In (ii), write that at η = ∞ the selected energy 𝒥(Θ, Π_𝒵T̂(γ)), γ optimal for ν_Z, is non-increasing — not min_Z𝒥, because the plan is
  Lloyd-consistent, not jointly optimal.
- Then either define the η = ∞ flow with a jointly optimal plan, implement it and keep D‴.17's argument, or label the convergence "argued
  (KL on the selected pieces), the limit a critical point of a Lloyd-consistent piece".
- Correct the Status, §0, §8, LITERATURE and D‴.23(iv) accordingly.

**D3p-4 (minor; Part B's designs and readings; the §0 row "σ-modification's bias 0.328"; RESULTS; PAPERS paper 11 "Claim" ("the classical
σ-modification keeps a large bias") and "Missing").**

The application compares the transport leaks with each other, and with the σ-modification toward ϑ⁰ = 0 at the transport leaks' σ₀. It
does not compare them with the classical design without a leak (σ₀ = 0), to which parameter projection adds the a-priori bound; the
projection would be inactive here.

*The no-leak baseline from Part B's own warm clouds* (`check_noleak.py`: the experiment's set-up with σ₀ = 0, three seeds, 2000 slow
units):

| design | fit error at m = 16, 64, 128 | floor at m = 16, 64, 128 |
|---|---|---|
| no leak | 0.0180, 0.01466, 0.01457 | 7.7·10⁻⁶, 2.4·10⁻⁷, 6.6·10⁻⁸ |
| adaptive prior, η = ∞ (the md) | 0.0298, 0.0148, 0.0147 | 4.85·10⁻⁵, 2.88·10⁻⁷, 6.36·10⁻⁸ |

The no-leak runs reach max|α| 5.2, 3.9, 2.6. So in the noise-free simulation no leak does at least as well as every transport design:
- The floors and biases that the adaptive and OT leaks remove are those a fixed leak creates.
- Their purpose, robustness, is never exercised: there is no disturbance, noise or unmodelled dynamics.
- The unprojected adaptive prior has no a-priori bound either (Theorem D‴.23(i)).

The baselines are also weak or missing:
- The σ-modification toward 0 at the transport leaks' σ₀ is the weakest classical choice.
- A Fourier regressor in the cogging harmonics, (sin 3x, cos 3x, sin 7x, cos 7x), represents τ_cog exactly with four linear parameters. So
  "no exact representation" is a property of the Gaussian network, not of the servo.

*Fix:* add the no-leak (and projection-only) baseline to Part B; call the σ-modification "toward ϑ⁰ = 0 at the transport leaks' σ₀"; in
paper 11's "Missing", list the classical baselines, a disturbance or noise for the leaks' robustness to act on, and the comparison with a
Fourier regressor of the known harmonics.

**D3p-5 (minor; D‴.29 *Status* l. 3358 "(ii) and (iii) formal — Theorem J's argument transferred to the fast variable (s, e), not
rewritten"; Proof (ii) l. 3349–3352 "Theorem J's averaging uses only the exponential stability of the error channel and the smoothness of
the regressor"; the §0 row "(Theorem J's argument transferred, not rewritten)"; (ii) l. 3328–3333; (iii) l. 3337–3340 with the Status's
"not the real-loop ones" l. 3359).**

(a) *What Theorem J proves.* Theorem J (l. 4964) is local: a periodic orbit near each hyperbolic critical point of the smooth fixed-prior
energy J_m, with its stability type, basin and location. Its proof uses the smooth leak (Step 3's normal form needs H to be C² in Θ), the
hyperbolicity, and f ∈ C² on a tube.

(ii) claims more: finite-horizon averaging "up to O(1/k_s)" for every §5d target — plan and cell switches, adaptive and projected priors.
- For switching targets, the first-order setting itself has only Theorem D‴.18's locally constant targets.
- The (s, e) channel is not Theorem J's: the e-transient decays at λ, not at k_s, and the error channel's Floquet multipliers are e^{−k_sT₀}
  and e^{−λT₀}. So Steps 1–2 have to be redone.
- The proof's sentence omits the smooth leak and the hyperbolicity.
- Numerically (ii) holds: see the verdict table (1.51·10⁻² → 6.31·10⁻⁴ for k_s = 25 → 200).

(b) *The prior's real-time law is not stated.* At finite η it is dz_j/dt = −(η/k_s)(z_j − T̂_j), and (ii) does not say so.
- Part C's η = σ₀ rows run with the prior frozen at its Part B end state (`real_job`'s `targets` l. 192–198).
- Over the 40 slow units it would move at most 5.5·10⁻³ (G⁻¹ norm; `check_partC.py`).

(c) *(iii)'s list.* It includes Proposition D‴.24 whole. D‴.24(vi) is a real-loop statement: perfect tracking is an invariant set of
Theorem J's loop, with the fixed prior's error "e = (1/k_e)r". For the servo the invariance holds, but the formula is (ii)'s. So the list
contradicts the Status's "not the real-loop ones". The list also carries D‴.27(ii)'s η = ∞ convergence (D3p-3).

*Fix:*
- Status and §0: "(ii) formal: first-order averaging with the fast pair (s, e) (s slaved at k_s, e filtered at the fixed λ). For the smooth
  fixed-prior leak near a hyperbolic critical point, Theorem J's steps carry over once the error-channel steps are redone for (s, e). For
  switching targets not even the first-order loop has more than Theorem D‴.18's locally constant case."
- State the prior's law in (ii), and either integrate the η = σ₀ prior in Part C or say that it is frozen.
- Exclude D‴.24(vi) from (iii), or restate it with (ii)'s filter.

**D3p-6 (minor; (ii) l. 3334–3336 "The exact filter (d/dt + k_s)⁻¹ keeps the error filter's attenuation of the residual's harmonics, which
the slaved value r/k_s drops"; *Numerics* (C) l. 3385–3391 and the md's (C); RESULTS l. 696–698; the §0 row "(the slaved one 0.93–0.996)";
LITERATURE l. 421–422, which lists the exact filter as specific; `real_job`'s docstring l. 184).**

(a) *The exact filter is not uniformly the better prediction.* Both forms are accurate to O(k_s⁻²) in s, i.e. to O(1/k_s) relatively. The
exact filter wins only where its own correction, n²ω₀²/(2k_s²), exceeds the other O(1/k_s) terms: the parameters' within-period motion
(Theorem J(a)'s (1/k)u(Θ, t)) and the nonlinearity r(z) − r(z_m). Evidence from `check_decompose.py` (independent servo, dense output over
the last of 12 periods):
- At test 3's start (fit error 0.040), relative errors of the predictions:

  | k_s | exact filter | slaved r/k_s | frozen-Θ nonlinear periodic error |
  |---|---|---|---|
  | 25 | 1.44% | 0.67% | 1.38% |
  | 50 | 0.93% | 0.76% | 0.90% |
  | 100 | 0.52% | 0.49% | 0.51% |

  The frozen-Θ nonlinear error is almost as far off as the exact filter, so the deviation is the parameters' motion; the nonlinearity
  contributes 0.06%.
- At Part C's matching η = ∞, m = 16 state (fit error 0.0153) the exact filter is off by 0.14%: the nonlinearity contributes 0.11% and the
  motion 0.02%. There the exact filter does win; the slaved form is off by −3.8%.

(b) *Part C's resolution.* The rms is computed from records 0.2 time units apart, 157 per window; `real_job`'s docstring says 0.1.
- For the matching η = ∞, m = 16 design, one more sample moves the rms by 0.23% (1.70437·10⁻⁵ → 1.70051·10⁻⁵).
- The dense value is 1.70310·10⁻⁵, a ratio of 1.00136 against the md's 1.0021 (`check_partC.py`).
- So "1.0000–1.0029" means "within 0.3% at the records' resolution", and its spread is not resolved.

(c) *"The harmonics n ≈ 7 of sin 7x".* At ϑ^b(q) the residual's power is spread over n = 4–13: n = 9, 5, 11, 7, 13 carry 13.2, 11.0, 10.8,
10.2 and 9.9% (`check_partC_more.py`). The e-filter shifts it to n = 2–7, and the slaved deficit computed from this spectrum is 2.7%
(n_eff = 5.8). The reading is right in size, but the residual is the approximation error's broad fine scale, not sin 7x.

(d) *Not a contribution.* The exact linear filter is the closed loop's error transfer function.

*Fix:*
- In (ii): "for a small, high-frequency residual the exact filter is the better prediction; in general both are accurate to O(1/k_s)
  relatively, the within-period motion of Θ contributing at the same order (at test 3's start the slaved form is closer)".
- In Part C, use dense records, correct the docstring, and describe the residual's spectrum as broad.
- In LITERATURE, drop the filter from the "specific" list.

**D3p-7 (minor; (iv) l. 3341–3347; *Status* l. 3362–3364 "the over-parametrisation of phase-plane kernels on circular orbits as an
excitation defect of the measure-space kind" and "exactly one zero Gauss–Newton eigenvalue per kernel, the next at least 10⁴ times larger";
the §0 row's last column; LITERATURE l. 422–423; PAPERS l. 159).**

(iv)'s computation is exact (see the verdict table). But three points are overstated.

(a) *The exact zero depends on units.* The phase-plane orbit of x_m = R sin(ω₀t) is a circle only for ω₀ = 1, in units where the kernel is
isotropic in (x, ẋ).

(b) *The excitation defect it describes is generic.* With test 2's atoms (`check_vonmises_abs.py`), count the whitened Gauss–Newton modes
below 10⁻⁶σ₀ — modes that an adaptive prior leaves practically neutral and a fixed prior pins at σ₀:

| orbit | modes below 10⁻⁶σ₀ | smallest |
|---|---|---|
| circle | 8 (5 exact, 3 more) | machine zero |
| ellipse, velocity axis ×1.1 | 8 | 8.8·10⁻¹⁵σ₀ |
| ellipse, velocity axis ×2 | 6 | 2.7·10⁻¹¹σ₀ |
| x_m = 1.2 sin t + 0.3 sin 2t | 3 | 5.4·10⁻⁹σ₀ |

So "on circular orbits" singles out the exactness, not the defect. The defect is the classical partial excitation of Gaussian networks along
an orbit ([6]–[8] of the reference list: Kurdila–Narcowich–Ward, Farrell, Wang–Hill), here with learned centres and widths.

(c) *"Exactly one per kernel" holds only in general position.*
- With the centres at the origin, 17 of the 20 eigenvalues vanish (`check_vonmises_centred.py`).
- "The next at least 10⁴ times larger" restates test 2's two thresholds (10⁻¹⁵ and 10⁻¹¹); the next eigenvalue is 2.4·10⁻¹⁰ of the largest.
- Test 2 does discriminate the circle: it fails on the ellipse ×1.001 (`check_vonmises.py`).

*Fix:*
- In (iv) and the Status: "on a circular orbit (ω₀ = 1, isotropic units) the restriction is a von Mises function: one exact zero per kernel
  in general position (more at degenerate atoms); on other orbits the same directions are only weakly excited (on test 2's atoms 6–8 modes
  below 10⁻⁶σ₀ on ellipses of aspect 1.1–2)".
- Cite [6]–[8], and restrict "new" to the exact fibre.

**D3p-8 (minor; D‴.27 *Numerics* (c) l. 3234–3236 "monotonically between the fixed prior's 0.1750 and the unprojected adaptive prior's
0.1050"; l. 3223–3225 "the instance's exact atoms are inside"; the §0 row "the box's amplitude bound 0.7–0.99 of a_max trading the fit
error 0.1736–0.1101 between …"; the §7 row; RESULTS l. 650–651; the §8 bullet "the box trading the fit error against the parameters'
bound").**

"Between" holds only on the scanned bounds (`check_boxscan.py`, the experiment's set-up, LSODA to s = 10⁶):
- Below them the projected adaptive prior is worse than the fixed prior: fit errors 0.1792, 0.1845, 0.1949 at amplitude bounds 0.65, 0.6,
  0.5 of a_max, against the fixed prior's 0.174984.
- The scan's own 0.7 row (0.1736) uses a box whose α-face, 3.4692, excludes the instance's exact atom (α = 3.8918, the 0.75 face). That
  contradicts "the exact atoms are inside" and Z(0) ∈ 𝒵^{K₀}; at η = ∞ the start is projected at once.

A box that cuts the representation is a bias of its own, and the design statement should say so.

*Fix:* "for the amplitude bounds 0.7–0.99 (0.7 already excluding the exact atom) the fit error lies between …; tighter boxes are worse than
the fixed prior (0.1792 at 0.65, 0.1949 at 0.5)"; in §8, "the box should contain the representation".

**D3p-9 (cosmetic; D‴.29 *Numerics* and rows).**
- l. 3381–3382 "end at 0.048 (OT; the nearest atom 0.10, 0.060, 0.048), below the replicated point's 0.0545": the nearest atom ends above
  it at m = 16 and 64 (0.1017, 0.0597; `check_partB.py`).
- The §0 row "the adaptive prior's 38–780 times lower with the fit error to the approximant's": at m = 16 the fit error is 0.0298, twice the
  approximant's.
- (A) "two modes far below σ₀": three are below 0.012σ₀ — 2.888·10⁻⁴, 5.396·10⁻⁴, 1.133·10⁻²σ₀ (`check_partA.py`).
- The §7 row "1.7–2.1·10⁻⁵ (adaptive, OT)" omits the misplaced OT design's 7.8·10⁻⁵.
- (ii) "n²/(2k_s²)" has ω₀ = 1 built in; the proof has n²ω₀².
- `real_job` re-implements `SlidingServo`, so the test and the experiment exercise different code. Both agree with my servo.

*Fix:* as indicated.

**D3p-10 (cosmetic; D‴.27).**
- (a) (iv) l. 3185 "If the limit has an atom on ∂𝒵 … the limit keeps a bias" needs T̂_j ∉ 𝒵: an atom on ∂𝒵 with T̂_j = z_j has no bias.
  That case is a stationary point of D‴.23(i) whose prior lies on the boundary. Write "an atom on ∂𝒵 with T̂_j ∉ 𝒵 … — necessarily so when
  no stationary point of D‴.23(i) has its prior in 𝒵".
- (b) *Numerics* (b) l. 3232 and the md "at least 17%": the margins are 16.91% (prior) and 16.88% (particles) (`partH_stored.py`). Write
  "16.9%".
- (c) *Numerics* (a) "the two steps' end states within 3.4·10⁻¹²" shows the shared fixed point, not the scheme's consistency. Print instead
  the trajectory error against a reference: at s = 300 it is 6.3·10⁻⁷ (dt = 0.25), 1.9·10⁻⁷ (0.125), 4.1·10⁻⁸ (0.0625)
  (`check_projected_step.py`).

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration, reproducibility, tests.**
   - `--from-raw` rewrites `servo_cogging.md` (md5 `6af6b7e7…`) and `adaptive_prior.md` (`f2156c5e…`) byte-identically.
   - Part B jobs (nearest-adaptive 0.3/0 and nearest-fixed 0.5/1) re-run bitwise, 53 s each.
   - The Part C job (nearest-adaptive, matching, m = 16) re-runs bitwise.
   - Part H re-runs identical to the stored entry in every printed value.
   - The 13 tests pass in 21.4 s; 140 are collected, as README says.
2. **D‴.29(i)–(ii).**
   - The error channel is re-derived; `SlidingRealLoop.rhs` satisfies it to 1.8·10⁻¹⁵, and its parameter law to 7.1·10⁻¹¹.
   - The averaging error is O(1/k_s): 1.51·10⁻², 4.61·10⁻³, 1.55·10⁻³, 6.31·10⁻⁴ at k_s = 25–200.
   - The position error against the exact filter is 1.0159, 1.0076, 1.0038, 1.0019 at test 3's start, and against the slaved form 1.0087,
     1.0057, 1.0033, 1.0018 (both within test 3's 1% at k_s = 50).
3. **D‴.29(iv).**
   - Along the fibre f̂ on the orbit is unchanged to 2.8·10⁻¹⁷.
   - Five machine zeros for five kernels on test 2's circle, the next at 2.4·10⁻¹⁰ of the largest.
   - At a fitted point the full Hessian has three zeros (4.2·10⁻¹³–1.4·10⁻¹²) where the Gauss–Newton block has three, then 1.77·10⁻⁴.
4. **Part A, by an independent analytic model.** `mymodel` derivatives are within 1.2·10⁻¹⁰ and 2.6·10⁻¹⁰ of finite differences, and the
   orbit grid and F are identical to the repository's.
   - ϑ^b(q) is reproduced (fit 0.014620118708317855 against 0.014620118708317859).
   - Ã's spectrum agrees to 2.7·10⁻¹¹ relative (the smallest eigenvalue to 4.9·10⁻¹⁰).
   - The rates are 1.444·10⁻⁴, 2.698·10⁻⁴, 5.651·10⁻³σ₀ (η = σ₀) and 2.888·10⁻⁴, 5.396·10⁻⁴, 1.133·10⁻²σ₀ (η = ∞); the fixed prior's
     1.0003, 1.0005, 1.0113σ₀.
   - λ⁺_k/σ₀ = 0.0004, 0.0016, 0.0098, 0.0075; the margin at w = q is 9.902·10⁻³; the W_g bounds are 0.0614, 0.0851, 0.1185, 0.0078.
   - The gauge transform at w = (0.3, 0.2, 0.3, 0.2) is stationary to 2.2·10⁻¹⁷ with the same function to 1.7·10⁻¹⁶, and its Ã is
     positive (3.5·10⁻⁴σ₀).
   - The floor law from finite differences of Newton solutions is 1.293966·10⁻² (relative difference 2.3·10⁻⁷); the adaptive law is
     5.7·10⁻²⁷.
   - The observed fixed-prior floors are 2.44, 1.47, 0.89 times the law at m = 16, 64, 128.
   - The misplaced prior: its fit error 0.2000547425245671 is identical; its stationary point has fit error 0.05452039; it is the same point
     from ϑ^b(q) (1.1·10⁻¹⁶); its cell margins are 1.83–2.18; κ_R/σ₀ = 0.474, −0.003, 0.263, 1.397; the margin is −3.965057·10⁻³ (2.5·10⁻¹¹
     from the repository's); the K₀-atom Hessian with the leak is positive (0.969σ₀).
5. **Part B and C numbers.**
   - Every md entry is recomputed from the raw records, including the fit errors from the stored f̂ against my own F.
   - Fixed/adaptive floor ratios: 38.1, 777, 707 (η = ∞) and 23.6, 106, 384 (η = σ₀).
   - The slaved ratios span 0.9316–0.9958 and the exact ones 1.0000–1.0029; my own slaved predictions equal the stored ones.
   - The parameters drift by 2.5·10⁻⁴–2.9·10⁻³ over 1000 time units. The prediction at the averaged end state is within 0.9% (0.998–1.009) of the one
     at the real servo's final parameters, so Part C also holds against the averaged designs.
6. **D‴.27.**
   - (i) is re-derived.
   - The new step is consistent at six steps; the old one stalls at four.
   - Test 10 discriminates: it fails at l. 281 with the old step.
   - Part H(c)'s six values are reproduced; the fixed prior's limit is 0.174984.
   - At finite η the plan is a joint minimiser, so (ii)'s monotonicity holds there.

## Not checked

- **Long runs.** Not rerun in full: 25 of Part B's 27 jobs and 17 of Part C's 18; `exp_adaptive_prior.py` outside Part H. The raw records
  were used.
- **Theory not derived.**
  - Theorem J's steps for the (s, e) channel were not rewritten.
  - The KL argument for the selected-energy flow at η = ∞ (D3p-3) is only sketched. No example of a non-Clarke-critical limit was found.
  - The existence of solutions with plan switching was assumed.
  - Theorems D‴, D‴.21, D‴.23 and Propositions D‴.5, D‴.12, D‴.17, D‴.24, D‴.28 were used as given.
- **Phase-plane dynamics.** No phase-plane dynamics were run beyond the tests: the experiment uses kernels on the position only.
- **Literature.**
  - The full texts of Sanner–Slotine 1992 (IEEE TNN 3(6):837–863; [Semantic Scholar record](https://www.semanticscholar.org/paper/Gaussian-Networks-for-Direct-Adaptive-Control-Sanner-Slotine/ae964fe0a0f394bc0f70a94585347d9118c189ed))
    and of Slotine–Li could not be reached: the proxy refused patents.google.com (403), and dspace.mit.edu failed DNS.
  - Their use of the filtered error with Gaussian networks on the state is from memory and from the paper's own citation.
  - That adaptive cogging compensation is a standard problem rests on search results: RBF networks compensating cogging in a
    transverse-flux machine ([IFAC 2005](https://skoge.folk.ntnu.no/prost/proceedings/ifac2005/Papers/Paper1584.html)), and periodic
    adaptive learning compensation of cogging in PMSM position servos ([Luo, Chen et al.](https://sites.ucmerced.edu/yqchen/publications/high-order-periodic-adaptive-learning-compensator-cogging-effect-pmsm-position-servo)).
- **Concurrent work.** 3609f86 (D‴.30) and the later uncommitted D‴.31 files were not reviewed.
