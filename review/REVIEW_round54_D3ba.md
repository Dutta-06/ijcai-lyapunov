# Referee report D3ba, round 54: Proposition D‴.58, round-48 revision

I reviewed commit 265a133 on branch claude/festive-bardeen-eopej2. During the review, other sessions moved HEAD to 849f5fe. `git diff --stat 265a133 HEAD` is empty for every reviewed file (draft, D3au report, the three source files, the experiment, the md and raw, the four test files), so this review applies unchanged at HEAD.

## Overall verdict

**What holds.**
- **(i) and (v) are exact.** I checked both from scratch at 80 digits, comparing closed-form Gaussian integrals with the Mehler series built on exact integer defects:
  - (i) agrees to 10⁻⁵²–10⁻⁸⁰ at τ/σ = ½, 1, 2 for K ≤ 12.
  - (v) agrees to 10⁻⁴⁴–10⁻⁸⁰ at σ*/σ = 0.5, 0.8 and 0.95; the repository never tests 0.95.
  - The analytic-continuation argument added to the proof of (i) is valid. Three small gaps remain (D3ba-9).
- **(ii)'s rates are numerically right.** At K = 40, τ/σ = 1, the ratios per kernel are 2.99984, 5.7918, 5.0023 and 3.7339, against 3, 5.828, 5 and 3.732. At τ/σ = 2 they are 1.49996, 2.6015, 2.0032 and 2.0032, against 1.5, 2.618, 2 and 2.
- **Reproduction.** The md regenerates from the raw byte for byte (md5 db199a20…), and all tests pass.
- **Numbers.** Every number in (A), (C) and (D) matches the md and the raw. So does most of (B), and "126 of 128", "26 of 30" and "69 of 69" reproduce.
- **C(c) ≥ 0 with double zeros** holds after all 147 optimal blocks of the three variants. The most negative value is −5.6·10⁻⁷ of max C, it sits exactly at a node, and there is no negative lobe between nodes.

**What does not hold.**
- **D3ba-1: the explanation of where the optimal ranges end is refuted.** (iii) attributes it to "the floor of a relaxation in double precision, whose squared error is rounded at ε_mach".
  - In all three polish-2σ paths I instrumented, the range ends exactly at the first block whose polish runs out of its `max_nfev = 5000` evaluations.
  - Raising that cap to 50000, and changing nothing else, extends the ranges by two kernels at τ/σ = ½, ¾ and 1. The optima then reached go down to 3.6·10⁻⁹–9.2·10⁻⁹.
  - The beyond-range numbers the draft quotes are wrong as well.
- **D3ba-2: the mechanism of flank entry is misattributed.**
  - I rebuilt the duals from the raw. The entrant is never the kernel where the current correlation peaks (0 of 128 entries).
  - What decides entry is the time to the wall, (1 − u_j)/C(c_j), with duals charged during earlier blocks.
  - The two "between-node" entrants entered because their duals were already at 0.992 and 0.989.
- **D3ba-3: the round-48 failure path in `sliding_iss` and `sliding_iss_nd` is broken.**
  - If a chunk fails, the path raises UnboundLocalError (first event) or IndexError.
  - With `exit_kkt=True` it instead silently applies the exit rule using the previous block's correlations.
- **D3ba-4: 38 of 56 mutants of the new and changed code survive.** The survivors include every cell of A3, A4, B's summary and C's ratio table that the draft quotes.
- **D3ba-5 to D3ba-11: smaller misquotes and wording.**
  - Misquotes: "96–100%", "(10⁻⁹ at τ/σ = ¼)", "1–3·10⁻⁷", "1·10⁻¹⁰–3·10⁻⁷", "not the stop".
  - (ii) says "q numerically" where the ordering is provable in one line.
  - (iv) states the stop deterministically and misattributes its e_K.
  - Literature points remain open.

## What I checked and how

- **Setup.**
  - I worked in a scratch copy and wrote nothing into the repository. Every process ran with OMP_NUM_THREADS=1 and nice 10, one heavy process at a time. I killed no processes.
  - The venv has no mpmath, so I added it to a scratch library directory. That directory already existed and is shared; I left it as it was.
  - The scratchpad itself turned out to be shared with other agents. Another process deleted my first copy (`scratchpad/repo`) at 12:01. My earlier results ran on a copy whose reviewed files `diff` shows equal to the committed ones. Everything later ran in a fresh, uniquely named copy.
- **Reproduction.**
  - `--from-raw` with WM_OUT in scratch gives a byte-identical md; I did it twice, in both copies.
  - Tests: `test_wide_feature.py` and `test_width_mismatch.py` pass (19). `test_rich_transit.py` passes (18), `test_moving_width.py` (3) and `test_noisy_certificate.py` (10).
- **Part A from the raw.** δ_t(K); the share of the gap recovered by `best_scale`; node ratios; first-term shares; the trapezoid-to-series differences.
- **Part B from the raw.**
  - All 18 paths: optimal ranges; the 69 non-optimal blocks with their errors and optima; ties; exits; node distances in every variant.
  - C(c) on a 32001-point grid after all 147 optimal blocks.
  - The duals along every path, rebuilt from the blocks' residuals and the event times. As a check, each entrant's rebuilt dual reaches its wall to within 4.7·10⁻⁸.
- **The stalls** (all 69 non-optimal blocks):
  - Stationarity measures; `polish_block` and variable-projection (VP) Levenberg–Marquardt (LM) started from each block.
  - Single-block re-relaxations: `sliding_iss`'s inner loop copied verbatim, restarted from the state just after the entry (the previous block plus the event list). With the original settings this reproduces the raw stalls exactly (9.0368·10⁻⁷ and 9.2441·10⁻⁸). I then varied one setting at a time.
  - The polish's termination status and evaluation counts; unbounded trust-region-reflective (TRF), LM and dogbox solvers from each block.
  - Jacobian condition numbers, and the error along the straight segment from block to optimum.
  - Part B's own `b_job` rerun with `max_nfev` raised from 5000 to 50000 for four polish-2σ paths and the base path at τ/σ = 1.
  - An instrumented rerun of three polish-2σ paths logging every polish call's status and evaluation count. It reproduces the committed errors exactly.
- **Part C per draw.**
  - The stop, rms/s at each K, the K of the largest relative gap, and the stop probabilities implied by the rms's draw-to-draw scatter.
  - A surrogate Monte Carlo of the discrepancy stop: 60 seeds per setting, δ = 0.05, 0.1, 0.2, N = 161, 321, 641, s = 10⁻⁴…10⁻². It uses VP fits started from the noiseless optima instead of the transit.
- **(i), (ii), (v)** in 80-digit arithmetic, as above. I re-derived the proof step by step: Bargmann prefactors, the Fock-space Gaussian integral as Mehler's kernel, the domination constant, the radius of the series, the eigenvalue (1 − ϑ)/(2(1 + ϑ)).
- **Failure path.** I monkeypatched `relax_chunk` to report a failed chunk (how = 3) at the first event and at a later one.
- **Mutation testing.** 56 single-replacement mutants, each run against the tests relevant to its code:
  - 16 in `relax_chunk`, `sliding_iss` and `leaves`, against `test_rich_transit.py` and `test_moving_width.py`.
  - 15 in `scaled_rule`, `scalings`, `best_scale`, `narrow_rule` and `narrow_error`, against both D58 test files.
  - 25 in `report()`'s A3, A4, B and C sections, in `part_A` and in `b_job`, against `test_width_mismatch.py`.
- **Literature.** Web-search snippets only. arxiv.org and ams.org are unreachable, and alphaxiv.org did not resolve, so I stopped there. Anything marked [memory] is from memory.

## Verdicts per part

**(i): correct; the proof is valid with small gaps.**
- The 80-digit closed form equals the series, and D_{2K,K} = 2^K K! holds.
- The new argument is sound: domination for small t, then analytic continuation.
  - Left side: the Mehler kernel integrated against Δ̃⊗Δ̃ is analytic near [0, 1), because the Gaussian part's quadratic form has smallest eigenvalue (1 − t)/(4(1 + t)) > 0.
  - Right side: the series has radius ≥ 1.
  - The two sides agree for t < C⁻² with C = 2, since E He_n(√2Z)² ≤ (4/3)4ⁿ n!.
- Gaps (D3ba-9):
  - The atoms Σ w_k|He_n(√2z_k)| need the same Cⁿ(n!)^{1/2} bound (Cramér's inequality); the text bounds only the Gaussian part.
  - "Uspensky's theorem in each variable" does not cover κ_ϑ(√2z, √2z′), which grows along the diagonal like e^{2ϑz²/(1+ϑ)}. The step needs the two-dimensional statement: a positive product rule, dominated by e^{a(z²+z′²)} with a = ϑ/(1 + ϑ) < ½.

**(ii): numerically right (to K = 40); one ordering is provable; wording.**
- "q numerically" understates it. q = (β² − 1)/(β² + 1) and t = (σ*/σ − 1)/(σ*/σ + 1); x ↦ (x − 1)/(x + 1) is increasing; and β² = (1 + 2τ²/σ²)^{1/2} > (1 + τ²/σ²)^{1/2} = σ*/σ. So q > t always.
- "A factor up to 2 per kernel, not 2^K" contradicts itself: a factor of 2 per kernel is 2^K. The intended statement is a gain of ≈ 0.6K^{−1/4}(λ/t)^K. (D3ba-8)

**(iii): the observation holds inside the reported ranges; the range's end and the entry mechanism do not.**
- **Holds:**
  - Inside the ranges the blocks are optimal to 1.000000, and their centres are within 7.6·10⁻⁷σ of the optimum's in every variant.
  - δ_t(K) is as quoted: 0.0776% → 0.0397%, 0.728% → 0.395%, 3.477% → 2.123%. The node deficits are 0.0077%–5.15%.
  - C ≥ 0 after the optimal blocks; interior bumps are at most 17% of the flank peak (median 0.4%).
  - 126 of 128 entrants start beyond the outermost node, counting entries from an optimal block into an optimal block.
  - No entry ties; the 6 exits all happen beyond the ranges.
  - "Best local optimum found": none of my refits from the stalled blocks beat the recorded optima.
- **Does not hold:**
  - D3ba-1: where the range ends.
  - D3ba-2: how kernels enter.
  - D3ba-5: "96–100%".

**(iv): numbers right and the error formula robust; the stop is stated too deterministically.**
- Stops: 6/6/3 reproduce, 26/30 at the predicted K.
- The stopped error is 0.92, 0.99, 1.02 times the bias–variance prediction; the ratio to the bound is ≤ 0.87; the ratio to the oracle is 1.10, 2.66, 1.93.
- In the Monte Carlo the stopped error over the bias–variance prediction has median 0.95–1.01 in all 27 settings.
- The stop itself is a margin question, and the wording should say so (D3ba-6, D3ba-7):
  - The rms's draw-to-draw scatter (≈ 0.04s here) implies P(stop at the predicted K) = 0.99, 0.79, 0.98. That gives 27.6 ± 1.4 of 30 expected, against 26 observed.
  - In the Monte Carlo, the hit rate ranges from 0.43 to 1.00.
  - At δ = 0.05 and 0.2 the misses are often one kernel earlier: 57% at N = 161, s = 10⁻³, δ = 0.05; 40% at N = 161, s = 10⁻⁴, δ = 0.2. So "the others one kernel later" holds only for δ = 0.1.

**(v): correct; the sketch should say where the measure lives.**
- Verified at 80 digits, including σ*/σ = 0.95.
- For t < 0, e^{tz²} = E e^{iWz} with W ~ N(0, 2|t|). N(0, 2t) is not a measure, so the representing measure lies on the imaginary axis, and the Gaussian integral gives Mehler's kernel at ρ = |t| in real variables. The sketch should say this.

**Status.**
- Fair for (i), (v) and the ₂F₁ conjecture.
- "The bias–variance split exact given e_K" overstates: it is first order in the noise (a linearised nonlinear fit).
- (iii)'s double-precision attribution has to go (D3ba-1).

**Numerics paragraph.** Everything matches the md and raw except:
- (iii) and (B) "96–100%";
- (B) "the blocks stalling at 1·10⁻¹⁰–3·10⁻⁷" (the actual range is 1.3·10⁻¹⁰–1.26·10⁻⁶);
- (C) "not the stop";
- (C) "the path to 10 entries" against B's "to 13 entries". Both use `max_entries = K + 1`, so C's should read 11.

**The round-48 fixes (D3au-1 … 13).**

| D3au item | verdict |
|---|---|
| 1 | Partly fixed. K* was dropped and three variants were run. The new "floor" repeats the over-interpretation; gtol and the integrator were never varied, and the polish's cap is shared by all variants. |
| 2 | Partly fixed. `relax_chunk` and its unit test are fine. The path-level failure branch is broken and untested (D3ba-3). |
| 3 | Fixed: 0–2 exits per path, all beyond the ranges, with no re-entry cycles. The path-level wiring is untested (R14, R15). |
| 4 | Fixed: the formula, the sample norm and the numbers. Wording remains (D3ba-6). |
| 5 | (v) itself is fixed and exact. The requested comparison with Derevianko–Plonka is not done (D3ba-10). |
| 6 | Partly fixed (D3ba-10). |
| 7 | Partly fixed. The new code's mutants mostly survive (D3ba-4). |
| 8 | Mostly fixed (D3ba-5). |
| 9 | Mostly fixed (D3ba-9). |
| 10 | Fixed, apart from wording. |
| 11 | The old wording is fixed; a new mechanism sentence is wrong (D3ba-2). |
| 12 | Fixed. |
| 13 | Fixed. |

**Tests: 18 of 56 mutants killed.**

| group | killed | survived |
|---|---|---|
| `relax_chunk` / `sliding_iss` / `leaves` (16) | R01 (LSODA failure ignored), R02, R12, R13 | R03 (Radau failure ignored), R04 (exhausted split reported as a success), R05 (second half started from z, not z1), R06 (first half's state returned), R07 (first-half failure ignored), R08 (RELAX_SPLITS = 0), R09 (n_redo counts stalls), R10 (no break after a failed chunk), R11 (n_fail never incremented), R14 (exit_kkt never applied in the path), R15 (own × nx²), R16 (≤ 0; near-equivalent) |
| `scaled_rule` / `scalings` / `best_scale` / `narrow_*` (15) | W01, W02, W06, W07, W09–W15 (all narrow-rule mutants) | W03 (Fock variance × cosh r), W04 (Mercer variance with 1 + τ²/σ²), W05 (Mercer variance without the root), W08 (`best_scale` xatol 10⁻²) |
| `report()`, `part_A`, `b_job` (25) | E01–E03 | E04 (stall/failure counts swapped), E05–E07 (B's summary line), E08–E10 (flank, ratios, event signs), E11 (stop at (1 + 2δ)s), E12–E14 (bias–variance, bound and oracle ratios), E15 (e_K from the grid variant; equivalent on the committed raw), E16 (gap index), E17 (runaway count), E18 (stopped error), E19–E23 (A3/A4 cells), E24–E25 (code no test runs) |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ba-1 | major | **"The range ends where the optimum's relative error falls below ≈ 10⁻⁷ (10⁻⁹ at τ/σ = ¼): beyond, every block … stays at 1–3·10⁻⁷ (10⁻¹⁰–10⁻⁸ at ¼) …: the floor of a relaxation in double precision, whose squared error is rounded at ε_mach" is refuted.** (a) *The numbers.* The 69 non-optimal blocks lie at 1.3·10⁻¹⁰–1.26·10⁻⁶. Outside τ/σ = ¼, only 19 of 45 lie in 1–3·10⁻⁷; 20 are below 10⁻⁷ (down to 4.0·10⁻¹⁰, all in polish-2σ or at τ/σ = 1) and 6 above 3·10⁻⁷ (up to 1.26·10⁻⁶). At ¼, K = 5's optimum (1.61·10⁻⁹, above 10⁻⁹) is missed in all variants. The grid-shifted path at τ/σ = 1 reaches an optimum at 5.22·10⁻⁸ (K = 10) right after missing one at 2.96·10⁻⁷. The md's own caption calls 10⁻¹⁰ "the double-precision floor". (b) *The blocks are not stationary points of the fit.* They pass a test with max\|g\| ≤ 9.4·10⁻¹⁰·scale (median 3.8·10⁻¹¹), under gtol = 10⁻⁹. Yet double-precision VP-LM from each block's own centres lowers all 69 (median ×47, up to ×5700), and reaches the recorded optimum in 25. (c) *The mechanism.* The polish (`polish_block`: bounded TRF, `max_nfev = 5000`, accepted once, then never repeated) ends with status 0, "max_nfev exceeded", exactly at the first non-optimal block of each polish-2σ path. Its evaluations per block grow about 3× per kernel: 74, 186, 524, 1407, 3796, then 5000 at τ/σ = 1; 2757, then 5000 at ¾; 1125, then 5000 at ½. Raising the cap to 50000, nothing else changed, gives optimal ranges K ≤ 7, 9, 11 instead of 5, 7, 9 at τ/σ = ½, ¾, 1, reaching optima at 3.6·10⁻⁹, 6.0·10⁻⁹ and 9.2·10⁻⁹. At ¼ the K = 5 block then sits 1.3% above its optimum. From the stalled blocks, the same bounded TRF with 50000 evaluations reaches the optimum exactly (5.2197·10⁻⁸ in 7338 evaluations; 5.2224·10⁻⁸ in 6297), as do unbounded TRF, LM and dogbox. The 1σ acceptance radius adds the base variant's stalls: at τ/σ = ½, K = 6 (9.0368·10⁻⁷, 22 rejections) the block reaches 8.32·10⁻⁸ in one chunk once the polish is accepted. The base path at τ/σ = 1 stays at K ≤ 9 even with the cap raised (54 rejections). The flow alone cannot finish: at τ/σ = 1, K = 10, gtol = 10⁻¹³ moves the block only 0.3% of the way to the optimum in the full 1.68·10¹². Nothing here is ε_mach. Across all widths e·cond J(w, c) ≈ 0.09–0.30, so the polish's cost grows like 1/e. That is why the cap binds at a similar e at every width: cond J ≈ 3.5·10⁶ at the τ/σ = 1 stall (VP: 4.0·10⁴), far below 1/√ε_mach = 6.7·10⁷. | (iii), (B), status, novelty: drop "double precision" and "ε_mach". Say the transit's relaxation stops reaching the optimum where its polish exhausts its evaluation budget or is rejected, and the flow (centre speed ∝ w) cannot finish within 1.68·10¹². Quote the beyond-range errors as 1.3·10⁻¹⁰–1.3·10⁻⁶ and variant-dependent; delete "(10⁻⁹ at ¼)"; replace "stationary" with "passed the test". Code: make the polish's `max_nfev` adaptive or larger; re-polish when the test passes, or polish by VP; record the polish's status and evaluation count in the raw. Then rerun B and report where (or whether) the ranges end. |
| D3ba-2 | major | **(iii)'s "so … the next to reach its wall is, as a rule, a flank kernel", and the between-node explanation "where the nodes are spread wider than the kernels' spacing", are not the mechanism.** In 128 entries from an optimal block into an optimal block: the entrant never had the largest current rate C(c_j) (0/128); its rate was a median 0.21 of the largest (range 0.001–0.78); its dual at the block's start was a median 0.96 (minimum 0.57). It sat a median 1.09σ beyond the outermost node, whereas C peaks 2.0–3.5σ out. It was the most-charged dual in only 45/128 cases: what selects it is the time to the wall, (1 − u_j)/C(c_j). The between-node entrants (τ/σ = 2, grid +0.05): kernel 11 at 0.15 sat 0.27σ from the node at 0.182, inside C's double zero, with rate 2.1·10⁻⁴ against a peak of 1.0·10⁻²; kernel 13 at 0.35 had rate 1.7·10⁻⁵ against 5.8·10⁻⁴. Their duals were already 0.992 and 0.989, charged during earlier blocks. | (iii): C ≥ 0 makes every dual grow monotonically. The kernels next to the block, charged while they sat under earlier blocks' flank peaks, reach the wall first; the exceptions are that history. Report each entrant's u at the block's start and its rate against the peak's in Part B. |
| D3ba-3 | major | **The failure branch added in round 48 breaks in `sliding_iss`** (`rich_transit.py` 381–430) and, by code reading, the same way in `sliding_iss_nd` (`rich_transit_nd.py` 270–316). On how = 3 at an iteration's first chunk, the loop breaks before wi, g, h, big and si (and ci in nd) are computed for this block, and the exit rule then uses undefined or stale arrays. Monkeypatched failures: at the first event, UnboundLocalError ('si'). At a later event with `exit_kkt=True` there is no error, but the previous 1-kernel block's arrays broadcast over the 2-kernel block. With `exit_kkt=False`, IndexError. B and C never trigger it (solver_failed = 0), and no test drives a failure through the path. | `sliding_iss`, `sliding_iss_nd`: on how = 3, recompute the correlations at the returned state, or skip exits for an unrelaxed block. Add a path-level test with a failing `relax_chunk` at the first and at a later event. Keep how = 1 vs 2 distinct through the split. |
| D3ba-4 | major | **The new code is mostly untested: 38 of 56 mutants survive** (Tests table above). Covered by no test: `relax_chunk`'s Radau and split path; the failure counts; `exit_kkt`'s wiring into the path; the Fock and Mercer-zero variances; and every cell of B's summary, C's ratio table and A3/A4. These cells carry the quoted "0 of 69 unrelaxed", "chunks redone 0–3/0–7", "0.92, 0.99, 1.02", "0.87 of the bound", "1.10, 2.66, 1.93", "largest gap at 5, 4, 2", A3's 63.50…/34.02… and A4's 2.92/8.77. | Tests: a failing Radau and a split chunk. A path whose kernel cycles without `exit_kkt` but not with it. `scalings` against closed forms (v_Fock = 2τ²σ²/(σ² + σ*²); the Mercer-zero variance as Fasshauer–McCourt's). `report()` on synthetic A3/A4/B/C with nonzero stall and failure counts, known ratios and an unrelaxed non-optimal block, asserting every quoted cell. |
| D3ba-5 | minor | **"96–100% of the gap recovered by retuning the rule's variance alone".** At τ/σ = 2 it is 93.0%–95.7% for K = 2–10, and 96% only at the last K. At τ/σ = ½, K = 8 it is 104.4%, which an optimum cannot exceed. So Part A's "optimal" fit there (2.097·10⁻¹⁰) is 1.75·10⁻⁵ worse, relatively, than the best one-parameter rule, and δ_t(8) ≥ 0.0415%, not 0.040%. | (iii), (A): write 93–100%. Mark Part A's optimum near 2·10⁻¹⁰ as unconverged, or use min(e_opt, e_best). |
| D3ba-6 | minor | **(iv) is stated deterministically, and its e_K is misattributed.** In each draw the stop is decided by the rms's scatter: at s = 10⁻³ the expected rms/s at K = 4 is 1.066 ± 0.043 against the threshold 1.1, so P ≈ 0.79 (7/10 observed). Per Monte Carlo, the misses go one kernel *earlier* at δ = 0.05 and 0.2, so "the others one kernel later" is δ = 0.1-specific. "(e_K the transit's optimal blocks)" is wrong for the stop: `predict_stop` and the "predicted" errors use the half-squeeze rule's L²(ℝ) error, and only the bias–variance ratio uses B's optima. K(s) ≈ log(‖f*‖_N/s)/log(1/t) gives 4.4, 3.1, 1.8 against the observed 6, 4, 3. Keeping the constant gives 5.2, 3.9, 2.6, hence 6, 4, 3: K(s) ≈ log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/(s(2δ + δ²)^{1/2}))/log(1/t). The status's "exact given e_K" is first order. | (iv): the predicted K is the first K whose *expected* rms is within the tolerance, decided per draw when a neighbouring K's margin is within ≈ s(1/(2N))^{1/2}. Name which e_K each prediction uses; give the K(s) formula with its constant. Status: "first order". |
| D3ba-7 | minor | **(C).** "The residual's largest relative gap … not the stop" fails at s = 10⁻³: there the largest gap falls at the stopping K in 4 of 10 draws, and both medians are 4. B says "to 13 entries" while C says "to 10 entries" under the same `max_entries = K + 1` convention; C's is 11. | Numerics (C), the md caption. |
| D3ba-8 | minor | **(ii)** "q numerically", although q > t is a one-line proof (above); and "a factor up to 2 per kernel, not 2^K" contradicts itself. | (ii), status. |
| D3ba-9 | minor | **Proof gaps.** The atoms' Cⁿ(n!)^{1/2} bound is missing. Uspensky is invoked "in each variable" where the two-dimensional dominated statement is needed (a = ϑ/(1 + ϑ) < ½). (v)'s measure lies on the imaginary axis, which gives Mehler's kernel at ρ = \|t\|. | Proof of (i) and of (v). |
| D3ba-10 | minor | **Literature and novelty, from snippets.** (a) D3au-5 asked to compare Derevianko–Plonka's rate with \|t\|; it is still not done. Per the snippets, they work in L²(ℝ, e^{−t²/(2ρ)}) with frequencies at the zeros of a scaled Hermite polynomial. Say whether these equal (v)'s η_k/σ², and whether their error is (v)'s series. Their journal version is Sampl. Theory Signal Process. Data Anal. 24(6), 2026, not "2023". (b) Per the snippets, Kuo–Woźniakowski 2012 analyse the *classical* Gauss–Hermite rule under a length-scale condition, so "at scaled Gauss–Hermite nodes" is misattributed. (c) A snippet reports Karvonen–Särkkä 2019 saying the classical rule converges in this RKHS iff ε² < ½, whereas (ii) has moment matching (the classical rule) converging at λ at every τ/σ. My 80-digit check confirms that convergence at τ/σ = 2 (e₄₀ = 1.1·10⁻⁷). Reconcile the normalisations or the claim. (d) Karvonen–Oates–Girolami 2021's scaling is unknown, so "New: the half-squeeze scaling" is unverified, and "whichever scaling these use" does not settle it. The squeezed-state discretisations of Janszky–Vinogradov and Szabó et al. are the same Bargmann computation [memory] and need checking. | Literature and novelty paragraph: hedge "New" until the papers are read. |
| D3ba-11 | cosmetic | The md's Part A caption ("not computed below 1e-10, the double-precision floor") contradicts (iii)'s ≈ 10⁻⁷. (B) says "all 69 stationary". The ~10⁻¹⁷ trapezoid-to-series difference is attributed to the GH nodes, but rounding of the residual gives the same order, and this was not tested. "The constant of (i) to within 1 − δ_t(K)" is ambiguous. C's docstring says "the first relaxed block", but `report()` never checks `relaxed`; all blocks happen to be relaxed here. | md caption, numerics, the experiment's docstring. |

## Lesson

A ceiling that every variant of a method hits in the same place belongs to what the variants share, not to the arithmetic. The three Part B variants shared the polish's 5000-evaluation cap, its once-only acceptance and an absolute stationarity tolerance. The revision credited the shared ceiling to double precision without the one test that separates the two: another solver, or the same solver with a larger budget, started from the stalled states the raw already held. That test took a minute and moved the ceiling by two orders of magnitude. Likewise, a mechanism offered for a path's decisions should be checked against the state that makes them. The duals could be rebuilt from the raw to 5·10⁻⁸, and they show the entrant is chosen by its time to the wall, not by where the residual's correlation currently peaks. A failure branch added to a solver loop needs a test that drives a failure through the loop, not just through the function that reports it.

Sources (search snippets only):
- [Derevianko–Plonka, arXiv:2307.13587](https://arxiv.org/pdf/2307.13587) and [journal version](https://link.springer.com/article/10.1007/s43670-026-00124-9)
- [Plonka–Riebe–Cuyt, arXiv:2606.03855](https://arxiv.org/pdf/2606.03855)
- [Karvonen–Oates–Girolami, arXiv:2004.12654](https://arxiv.org/abs/2004.12654), [ResearchGate](https://www.researchgate.net/publication/340963491_Integration_in_reproducing_kernel_Hilbert_spaces_of_Gaussian_kernels) and [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/)
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532)
