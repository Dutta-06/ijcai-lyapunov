# Referee report D3bq, round 70: Proposition D‴.58 as revised after round 65

## Scope

**Commit and files.** I reviewed commit 409ec98 on branch claude/festive-bardeen-eopej2. I exported it with `git archive` to `scratchpad/ref70/D3bq_work/repo` and worked only there.
- Files reviewed: `review/drafts/D58.md`, `experiments/exp_width_mismatch.py`, `results/chaos/width_mismatch.md` and `_raw.npz`, `tests/test_width_mismatch.py` (27 tests), `tests/test_rich_transit.py`, `tests/test_wide_feature.py`, `src/core/wide_feature.py`, `src/core/rich_transit.py`, `review/REVIEW_round65_D3bl.md` and the "Round 65" entry of `review/REVIEW_LOG.md`.
- Since 3eb9326 (round 65) only the draft, the experiment (+333 lines), the md, the raw and `test_width_mismatch.py` (+242 lines) changed.
- HEAD moved to 42aa6a3 during the review. `git diff 409ec98 HEAD` is empty for every reviewed file.

**How I worked.**
- All runs used the venv Python with OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1 and `nice -n 10`, with at most two heavy processes of mine at a time.
- WM_OUT pointed at a scratch copy of the raw. Mutants ran in slim copies under `D3bq_work/mut/<id>`.
- Nothing was written to the repository (`git status` is clean). No process of mine is left running.
- Each number is marked [run] (my run), [md] or [raw].

## Summary

**What checks out**
- **Reproduction.**
  - `--from-raw` regenerates the md byte for byte (md5 302648b6…), and the raw is untouched (5b603943…) [run].
  - All 66 tests pass (27 + 27 + 12, in 52.7 s) [run].
  - These reruns with the experiment's own functions are bit-identical to the raw [run]:
    - `b3_job` at (½, 5·10⁴ once, K = 8): 0/50000 2.67693; 0/50000 2.67693; 3/227797 1.000274.
    - `b3_job` at (¾, sixth variant, K = 11): 0/50000 2.28692 twice; 3/384829 1.0036518.
    - B3's 5·10⁵ end at (¼, sixth, K = 6): 3/256857 1.014736573.
    - `part_A5`: all four rows, and the 60/60 check.
    - `b2_job` for two far misses (½ base K = 8; ¾ base K = 10) and for the 31.7 outlier (¾ grid, K = 8).
    - The (½, polish 2σ), (1, 5·10⁴ once) and (¼, base) paths up to their first non-optimal block: blocks and polish logs.
    - B3's "path" end (move and accepted flag) at all 22 rows.
- **Numbers.** Every printed number in (iii), the Status and Numerics A–D matches the md or the raw. Recomputed from the raw [run]:
  - 316 optimal blocks, optima down to 6.81·10⁻¹⁰.
  - 116 non-optimal blocks, all relaxed: optima ≤ 2.96·10⁻⁷, errors 1.01·10⁻¹² … 1.26·10⁻⁶.
  - B3 above the floor: 20 ends, 1530–227797 evaluations, 0.99998–1.00028 times the optimum, 16 reached with 5·10⁴.
  - Rejected "path" polishes: 1.93, 1.71, 1.83, 1.07 and 1.27σ.
  - K = 1 flank peaks: 2.02–2.57σ.
- **Mathematics.** (i), (ii), (iv), (v) and the proof are unchanged. I re-derived the Mehler-frame paragraph by hand:
  - The cross term gives ℓ² = 2σ²ρ/(1 − ρ²) and g(y) ∝ e^{−(1−ρ)y²/(4σ²)}. So g·N(0, τ²) ∝ N(0, v_ρ) with v_ρ = 2σ²τ²/(2σ² + (1 − ρ)τ²), and the Gauss–Hermite rule of N(0, v_ρ) is exact below 2K.
  - D_{2K}(v) = v^K K! and D_{2K+2}(v) = K![−v^{K+1}K(K + 2) + (K + 1)(2K + 1)v^K(v − 1)]. This gives the stated first tail ratio.
  - Its K² coefficient v′ − 2 vanishes only at ρ = t: (1 − ρ)/√ρ = 2σ/τ has the single root ρ = t. Also v′ = 1 if and only if ρ = λ.
  - Fasshauer–McCourt (ε² = 1/(4σ²), α² = 1/(2τ²)) is the frame ρ = q, with ℓ_q² = τ²/β², e^{−δ²y²} = g, its own Gaussian N(0, 2τ²/(β² + 1)) and v′_q = 2β²/(β² + 1).
  - The Mercer-zero measure is N(0, ℓ_q²). q ≠ λ for τ > 0.
  - The √2 bound holds (k ≥ 0, Σw = 1). The relative L² error is (σ*/σ)^{1/2} times the normalised one, and ‖m_μ‖² = σ/σ*.
- **The new (iii) claims hold in substance.**
  - `seed_rule` reproduces `sliding_iss`'s first-polish acceptance: min(start, target) ≤ step holds if and only if either does; same 10⁻² seed set and same residual; the signs agree for these positive amplitudes.
  - The budget binds at every end I restarted from the path's own polish start (3 ends).
  - The seed rule rejects the optimum itself at all five radius ends.
  - The floor is arithmetic. With the residual computed in long double, the same Levenberg–Marquardt solve that stalls in double reaches the optimum from both below-floor blocks.
- **Entries.**
  - 468 entries = 36 first entries + 316 joining an optimal block + 116 joining a non-optimal block.
  - After optimal blocks: 107 events — 1 entrant of the largest current rate, 105 beyond the outermost node (median 1.17σ) and 2 between nodes (duals 0.989–0.992). 6 events fall in two rows. All [run].
  - The 107 is stable for dual thresholds from 10⁻⁵ to 10⁻²; it becomes 110 at 10⁻⁶ and 126 at 10⁻⁸ [run].
- **Round-65 issues.** D3bl-3 and D3bl-5 to D3bl-10 are fixed. D3bl-1 is mostly fixed and D3bl-2 partly (status table).

**What is wrong (majors first)**
1. **D3bq-1 (major): the round's new code is weakly tested.**
   - I ran 91 single-point mutants plus a no-op control. The tests kill 39 and `--from-raw` kills 7 more. 45 survive both, and 25 of those change what a rerun would publish.
   - The B3 test uses only the variant without a gradient test, so passing or dropping the test is invisible (J01, J02).
   - The seed rule's ascent-target branch is untested (S02, S05, S06). One such mutant flips a committed radius verdict.
   - `part_B3` and `part_A5` have no test. A01, which reinstates round 65's Mercer-zero error, survives.
   - B2's far-miss diagnostics are untested.
   - The round-65 "no merge" mutant (X29) still survives.
2. **D3bq-2 (major): B3 does not test what the path did.**
   - It restarts from the block's final state, not from the path's polish start (17–63× above the optimum in my reruns).
   - Its extra columns change budget and gradient test together, and their ends are never judged by the seed rule.
   - At two ends, B3's own "path" column shows the path's settings succeeding. From the recorded starts, I find the budget does bind there. So the claim survives, but B3 does not show it.
3. **D3bq-3 (minor): the floor.**
   - "Censored by the floor" is right (my extended-precision runs) but is not shown by anything in the repository.
   - At ¾ the reference optimum is 0.16% too high.
   - The floor is soft: the step test's excess rises smoothly as the optimum falls.
4. **D3bq-4 (minor):** the claim that the two problems do "not [differ] in their slow modes" holds only with the centres in units of σ.
5. **D3bq-5 (minor):** an "event" also merges entries after the variants have parted.
6. **D3bq-6 (minor):** the Mercer frame's rule is asymptotically slower than the half-squeeze rule, not "nearly as fast".
7. **D3bq-7, D3bq-8 (cosmetic).**

## What I checked and how

- **Reproduction.** As listed above.
- **B3, row by row.**
  - I applied the path's seed rule at the variable-projection (VP) optimum, relative to the final state, at all 22 rows.
  - I wrapped `polish_block` to record its inputs and reran three paths to their first non-optimal block (`max_entries` = block index + 2).
  - From the recorded starts I restarted with 5·10⁵ evaluations, with and without the path's gradient test, and judged each end by the path's acceptance rule (for repolishes, also the 10⁻⁹ improvement test).
  - I reran B3's "path" polish at all 22 rows and scored it with `seed_rule` and seven mutated versions of it.
- **The floor.** From both below-floor blocks:
  - TRF started at the VP optimum.
  - TRF and MINPACK LM restarted from their own ends.
  - Gauss–Newton in double and in long double.
  - A Levenberg–Marquardt solve (Nielsen damping, column scaling) with the residual in double versus in long double, everything else in double.
  - The predicted and actual decrease of small steps, in both precisions.
  - VP from the block.
  - Also, from the raw: the 5·10⁵ end's excess over the optimum against the optimum's level, at all 22 rows.
- **B2.** σ_min ratios and condition numbers at the 20 above-floor first blocks: centres in σ, in length units, and with equilibrated columns. Also Kaufman's projected Jacobian, and whether VP's raw ends are sorted.
- **Entries.** I regrouped the 468 entries by (width, grid, iteration, kernel) for dual thresholds 10⁻⁸ … 10⁻² and listed the straddling events.
- **Scaled rules.** Per-kernel ratios up to K = 30–34 at τ/σ = 1 and 2, for v_hs, the Mercer frame's variance and v_hs·(1 ± 0.05).
- **Mutation testing.** 92 mutants (X00 a no-op control) in:
  - the B3 code: `seed_rule` (S01–S08), `b3_job` (J01–J09), `part_B3` (P01–P02), `b3_lines` (L01–L09);
  - the report code: `polish_cell` (C01–C03), `entry_lines` (E01–E09);
  - B2: `b2_job` (B01–B09), `vp_fit` (V01–V03), `vp_jacobian` (G01–G03);
  - A5: `gh_defect_v` (H01–H04), `frame` (F01–F02), `part_A5`/`a5_lines` (A01–A09);
  - the summaries: `cs_margin`/`cs_vacuous` (M01–M06), B's summary lines (R01–R04), `b2_line` (N01–N04);
  - plumbing: `main`/`part_B`/`merge_B` (X01–X07).
  - Each ran against `test_width_mismatch.py` (with `-x`, deselecting `test_b_paths_reproduce_the_raw`, which calls no mutated code) and against `--from-raw` (byte diff of the md).
- **Literature.** The shared web-search budget was exhausted at my first query, and arxiv.org does not resolve. The literature points rest on memory and on round 65's recorded snippets.

## Verdicts per part

| claim (D58.md line) | verdict | evidence |
|---|---|---|
| (i), proof (:8–14, :82–112) | Unchanged; correct | word-diff 3eb9326→409ec98 |
| (ii) (:15–23) | Unchanged; correct | λ/t = 1 + sech r, the Fock variance, q > t (by hand) |
| (iii) δ_t, 93.0–100%, nodes 0.008–5.2% (:24–27) | Correct | md A, A3 |
| (iii) C ≥ 0, peaks 2.0–2.6σ at K = 1 (:28–30) | Correct | 2.02–2.57σ [run from raw] |
| (iii) entries: 1 of 107, 316, 0.21, 0.963, 105 at 1.17σ, 2 at 0.989–0.992 (:30–34) | Numbers correct; "until they part" is inaccurate | [run]; D3bq-5 |
| (iii) ranges and optima down to 6.8·10⁻¹⁰ (:34–40) | Correct | md B |
| (iii) "the ranges are the solver's (B3)"; above the floor the budget, radius or gradient test binds (:40–46) | Holds where tested; B3 does not show it; the quoted moves (1.07–1.93σ) belong to a non-optimal end | D3bq-2 |
| (iii) gradient-test stops 1.012, 2.677, 1.121; "met at once" (:46–49) | Correct | raw logs; B3 "path" 1/1 [raw] |
| (iii) sixth variant one kernel further; 1.233 (:49–50) | Correct | md |
| (iii) below the floor: 1.0147, 1.0037, "whatever the budget", censored (:50–52) | Numbers reproduce; the floor is arithmetic [run] but the repository does not show it; ¾ reference unconverged; soft | D3bq-3 |
| (iii) VP 20/20; "as does the polish's own solve" (:52–54) | Correct | raw B2, B3 |
| (iii) σ_min 1–31.7 (median 1.01); condition numbers (:54–57) | Numbers correct; the "slow modes" reading depends on units | D3bq-4 |
| (iii) 33/38; 1.001; 4 at 16.36–24.74; coalescing (:57–60) | Correct | md; two far misses rerun identically [run] |
| (iii) 116/116; 53; 63 at 3.5·10⁻⁴–0.97 (:60–62) | Correct | md |
| (iv) (:63–75), (v) (:76–80) | Unchanged; correct | md C, A4 |
| Status: √2 bound, (2σ*/σ)^{1/2} (:134–137) | Correct | by hand |
| Status: Mehler frames (:142–151) | Algebra correct; "nearly as fast" holds only at finite K; one ambiguous clause | D3bq-6 |
| Status: literature and "New" (:128–164) | Hedges unchanged and appropriate | Sources |
| Numerics A–D (:166–202) | Correct | md/raw; A5 rerun identical [run] |

## Status of the round-65 issues

| D3bl | status | note |
|---|---|---|
| 1 | mostly fixed | B3 added. (iii) and the status now set the two below-floor ends apart as censored, and my extended-precision runs confirm the floor is arithmetic. Remaining: the repository shows only TRF's step test; the ¾ reference optimum is 0.16% high; B3's design (D3bq-2, D3bq-3). |
| 2 | partly | Tests now kill analogues of X01, X02, X04/X05/X07, R01, X28 and X30, and the listed boundaries. X29 (no merge) still survives as my X02: the new test reruns both stored paths, so merging and replacing coincide. The new code adds 25 behaviour-changing survivors (D3bq-1). |
| 3 | fixed | Stall wording, gaps 0.000–0.021σ, cosines and stops are reported. Note: at all four far misses VP's raw ends are unsorted, i.e. the coalescing pair actually crossed [run]. |
| 4 | partly | Analytic Golub–Pereyra Jacobian (tested), units stated, σ_min reported, "far better conditioned" dropped. The new slow-mode claim depends on units (D3bq-4). |
| 5 | fixed | Frame algebra right; residual issue in D3bq-6. |
| 6 | fixed | 107 events [run]; wording in D3bq-5. |
| 7 | fixed | Margins quoted [md]. |
| 8 | fixed | "Below 10⁻¹⁵"; the grid row at the optimum is rejected (1.27σ) [md]. |
| 9 | fixed | Status :125–126. |
| 10 | fixed | Sentence replaced. |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bq-1 | major | **Mutation testing of the new code.** 91 mutants plus a control (the control survives both, as it should). The tests kill 39; `--from-raw` kills 7 more (E04, A06, A07, M02, R04, N02, X04); 45 survive both.<br>**20 survivors are equivalent on the committed data:** boundaries never reached, positive amplitudes, and one exact equivalent (H03: p[K] is zero after the reduction).<br>**25 survivors change what a rerun publishes:**<br>• **J01 and J02 (gradient-test plumbing).** The B3 "path" column without the gradient test, or the "no gradient test" columns with it. Both pass because the only B3 test uses the sixth variant, whose gradient test is None. On the data they swap ¼ once/repeated's 1/1, 1.01196 and 3/1530, 1.00000 [raw].<br>• **J06 (seed rule judged at the start).** All five radius rows become "accepted, 0.00σ" [run].<br>• **J08 (move and acceptance attached to the other columns).** `b3_lines` raises an error on a rerun.<br>• **S06 (ascent with the opposite sign).** Flips (¾, grid, K = 8) from 1.07σ rejected to 0.64σ accepted [run].<br>• **S02 and S05 (ascent target).** Change the printed moves. The quoted 1.27σ at the grid row comes from the ascent target; without that branch it would be 1.32σ [run].<br>• **S07 and S08.** The seed set at 10⁻³ instead of `sliding_iss`'s 10⁻², or min instead of max over the seed kernels: moves change (1.71σ → 0.18σ at ¾ base) [run]. The only seed-rule test has a single seed kernel whose ascent target lies farther away than its start.<br>• **P01, P02, X05 (`part_B3`, untested).** The last non-optimal block, ratios over the block's own error, the first six paths only.<br>• **A01–A05, A09 (`part_A5`, untested).** A01 reinstates the Mercer-zero rule — round 65's own error. A rerun would print 34.02, 9.916, 3.751, 2.021 (md A3) instead of 63.49, 17.17, 5.321, 2.332.<br>• **B04, B05, B06, B09, V02, V03 (B2 far-miss diagnostics).** V03 matters because VP's ends at all four far misses come out unsorted (−1.5·10⁻⁵ … −0.021σ) [run].<br>• **X02 (`--variants` without the merge).** | • Test `b3_job` with a variant that has a gradient test; assert the calls (path budget, 10⁻¹⁵) and then (n, None).<br>• Test `seed_rule` with two seed kernels, one decided by its ascent target; pin its wfrac to `sliding_iss`'s default.<br>• Test `part_B3` and `part_A5` on synthetic data, with v_ρ, 2τ²/(β² + 1) and e_mframe checked by hand.<br>• Test `b2_job` on a case where VP's end differs from the block.<br>• Test the order of `vp_fit`'s return values and its sorting.<br>• Run `--variants 5` with stored paths that are not rerun. |
| D3bq-2 | major | **B3 restarts from the block's final state, not from the polish start, and changes two settings at once** (:40–46; md :168–197).<br>**Different state.** From recorded paths [run]: the path polished at 17.05× (½, polish 2σ, K = 6), 57.90× (1, 5·10⁴ once, K = 12) and 62.98× (¼, base, K = 5). Those polishes ended at 1.461×, 1.465× and 10.85×, essentially the final blocks. B3's "path" column is therefore a second polish.<br>**Two ends where B3 shows nothing binding.** There the path's own settings succeed from the final state: 1/1788, 1.00000, accepted; and 0/50000, 1.00005, accepted [raw]. The B3 summary is silent on both.<br>**From the recorded starts, 5·10⁵ evaluations [run]:**<br>• (½, polish 2σ): 1/6635, 1.0000000, accepted (move 1.51σ ≤ 2σ) — the budget binds.<br>• (1, once): 1/102753, 0.999999, accepted — the budget binds.<br>• (¼, base): with the test, 1.0017 (gradient test, 31071 evaluations); without, 1.0000007 — the budget and the gradient test both bind.<br>**Radius.** The quoted 1.07–1.93σ are the moves of the 5000-evaluation end. The seed rule also rejects the optimum itself, at 1.98, 1.84, 1.89, 1.22 and 1.27σ [run]. | Record the polish starts in `b_job` (wrap `polish_block`; these reruns took 4–91 s). Restart from them one setting at a time: 5·10⁵ with the path's test; the test off; the radius off. Judge every end by the path's own acceptance and tabulate which setting binds per end. Relabel the current column as a restart from the final state. |
| D3bq-3 | minor | **The floor is real but not shown, and it is soft** (:41, :50–52, :125).<br>**Arithmetic [run].** The same Levenberg–Marquardt solve from the path's blocks stalls at 1.00224 and 1.00095 with the residual in double, and reaches 0.99997 and 0.99837 with it in long double (67847 and 151803 iterations). At the ¼ block, small steps are pure rounding in double (predicted 9.9·10⁻⁹F, measured 2.5·10⁻⁷F) and agree in long double (9.88·10⁻⁹F).<br>**Other double-precision solves stop elsewhere [run]:** MINPACK LM at 1.00224 and 1.00140 (restarted, 1.000125 at ¾); TRF restarted at 1.01199 and 1.00337; Gauss–Newton makes no progress.<br>**The ¾ reference.** At (¾, K = 11), TRF started at the VP optimum ends at 0.99845×. So the reference optimum is 0.16% too high, and by the experiment's own yardstick LM in double reaches "the optimum".<br>**Soft floor.** The 5·10⁵ excess over the optimum is ≤ 1.3·10⁻⁷ for optima ≥ 1.6·10⁻⁹, 1.5·10⁻⁵ at 6.8·10⁻¹⁰, 2.7·10⁻⁴ at 2.1·10⁻¹⁰, 3.7·10⁻³ at 7.8·10⁻¹¹ and 1.5·10⁻² at 2.6·10⁻¹¹ [raw]. | Add the extended-precision evidence. Note that the ¾ reference is unconverged. Describe 10⁻¹⁰ as the level where the step test's excess crosses 0.1%. Note that "whatever the budget" applies to a single solve. |
| D3bq-4 | minor | **The slow-mode comparison depends on units** (:54–57; md :165).<br>• σ_min(VP)/σ_min(joint), centres in σ: 1.0001–31.75, median 1.011, 17 of 20 below 1.1.<br>• In the polish's own length units: median 1.54, 4 of 20 below 1.1.<br>• With equilibrated columns: 49–2786. All [run].<br>• At the five radius rows the joint fit's smallest direction is a near-seed kernel's centre: Kaufman's projected Jacobian gives exactly 1.000 there. The 31.75 outlier comes from the Golub–Pereyra residual term at a kernel with \|w\|/max = 2.4·10⁻⁷.<br>• The joint condition numbers in the polish's units are 1.3·10⁶–5.0·10⁸, not 4.5·10⁶–2.3·10⁹ [run].<br>• The md's "slow modes as slow in both" is unqualified. | State the units and scale-dependence, or drop the claim. Quote the condition numbers in the polish's own units. |
| D3bq-5 | minor | **"Events" merge entries after the paths part** (:31–32). The 6 two-row events join an optimal block in one variant and a block 1.01–14.3× above its optimum in another. Example: ¼, iteration 5, kernel 7 — the sixth variant's block is optimal, base's is 10.85×, and the duals agree (0.9253085–0.9253116) [run]. | Write "the same entrant at the same iteration, whether or not the paths have parted". |
| D3bq-6 | minor | **The Mercer frame's rule is slower, not "nearly as fast"** (:146–148). Its per-kernel ratio converges to 2.333 at τ/σ = 2 (K = 34) and 5.313 at τ/σ = 1, against 1/t = 2.618 and 5.828. The half-squeeze rule's ratio is still rising at the last K (2.596 at K = 30) [run].<br>In that frame the first term would decay at q·v′/2 = 0.375 < t, so the tail sets the rate — consistent with the (v′ − 2)K² growth.<br>The clause "a frame's own Gaussian only at ρ = λ" can be misread: for τ/σ > 2, N(0, τ²/β²) is the own Gaussian of the frame ρ = (β² − 3)/(β² + 1). | Quote the limiting rate, and reword the clause as "a frame's orthogonality measure is its own Gaussian only at ρ = λ". |
| D3bq-7 | cosmetic | `seed_rule` omits the repolish's improvement test, so B3 prints "accepted" for the 1/1 polishes that B's table logs as "1/1/r". The column "the path's polish" is a restart from the final state. | Add the test, or relabel. |
| D3bq-8 | cosmetic | Numbers no script prints: :116 (j ≤ 10, K ≤ 40), :122 (K = 640, 1280) and :136 (80-digit arithmetic). They are attributed to earlier referee rounds. | Mark them as not reproducible from the repository, or script them. |

## Sources

- I read no full text. The shared web-search budget was exhausted at my first query, and arxiv.org did not resolve (`getaddrinfo ENOTFOUND`). arxiv.org must be allowed in the session's network settings before these can be checked.
- From memory:
  - Karvonen–Oates–Girolami, Math. Comp. 90 (2021): the first class of rules is "an appropriate scaling of the classical Gauss–Hermite rules", with bounds e^{−c₁N}N^{1/4} and e^{−c₂N}N^{−1/4}, c₁ > c₂. This matches :133. Which scaling they use is unchecked; the novelty of the scaling itself depends on it.
  - Karvonen–Särkkä, BIT 59 (2019), "Gaussian kernel quadrature at scaled Gauss–Hermite nodes": consistent with the Mercer-zero scaling.
  - Kuo–Woźniakowski, BIT 52 (2012) 425–436, and Kuo–Sloan–Woźniakowski, Math. Comp. 86 (2017): titles only.
  - Gnewuch–Hefter–Hinrichs–Ritter, J. Complexity 71 (2022): isometry between Gaussian-kernel and Hermite spaces, preserving point evaluations.
  - Irrgeher–Kritzer–Leobacher–Pillichshammer, J. Complexity 31 (2015).
  - Janszky–Vinogradov, PRL 64 (1990) 2771; Szabó et al., PRA 53 (1996); Tang, SIAM J. Sci. Comput. 14 (1993).
- Not checkable here: Derevianko–Plonka (2026), Plonka–Riebe–Cuyt (2026) and Gnewuch–Ritter–Rüßmann (arXiv:2412.05368). The draft's hedges ("unchecked", [snippets]) remain appropriate.

## Lesson

A restart is evidence only about the state it starts from. B3 restarted each range end from the block the path left behind, which is in effect a second polish. At two ends it therefore showed the path's own settings succeeding where the path had not; only the recorded polish starts, 17× and 58× above the optimum, show the budget binding. A floor is a claim about precision, and the clean test changes only the precision: with the residual in long double, the same solve reaches the optimum from both blocks. That is what makes "arithmetic" right; step-test stops alone could not. Finally, a test that pins plumbing must exercise it where it matters. The B3 test uses the one variant whose gradient test is None, so passing or dropping that test leaves the suite green — and the same blind spot hides the seed rule's ascent branch, which decides one committed verdict.
