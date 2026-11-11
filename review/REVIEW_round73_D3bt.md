# Referee report D3bt, round 73: Proposition D‴.52 as revised after round 71

## Scope

- **Commit reviewed:** 1f51d48 on branch claude/festive-bardeen-eopej2 ("D‴.52 round 71: the minimum gains nothing at the tie; the one-sided half at equal size").
  - I exported it with `git archive 1f51d48` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/repo` and worked only there.
  - For diffs, round 71's raw, md, experiment and tests (42aa6a3) are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/old71`.
- **Files reviewed:**
  - `review/drafts/D52.md` (256 lines);
  - `experiments/exp_close_features.py` (1437 lines);
  - `results/chaos/close_features.md` (303 lines) and `close_features_raw.npz`;
  - `tests/test_close_features.py` (43 tests);
  - `src/core/rich_transit.py` (`Bumps`, l. 30–49);
  - `review/REVIEW_round71_D3br.md`, and the Round 67 and Round 71 entries of `review/REVIEW_LOG.md`;
  - the diff 42aa6a3..1f51d48 of the experiment, the tests, D52.md and `review/drafts/README.md`.
- **HEAD moved** to 18bcfa0 ("D‴.59 round 68") while I worked. None of the reviewed files differs between 1f51d48 and HEAD [own]. I reviewed 1f51d48 only.
- **How I ran things:**
  - the repo's venv with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`;
  - `CF_OUT` pointed to a scratch copy of the raw. The one exception: one `--from-raw` ran without `CF_OUT` and rewrote the export's own md, inside my scratch copy. The result is byte-identical to the commit's md [own];
  - at most two heavy processes of mine at a time, on a machine at load 7–10 from other jobs;
  - the background time limit stopped one mutation runner after U71-02 and one idle waiter. I reran the remaining U71 mutants in a fresh runner;
  - nothing of mine is left running, and the repository's `git status` is clean.
- **Tags:** every number is marked [own] (my run), [md] or [raw].
- **I modified nothing in the repository.**

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (`cmp`) [own].
  - The 43 tests pass in 24.8 s [own].
  - Against 42aa6a3's raw, A, B, C, E, F and S are identical. In U only `q2` changed (the new law seed), and `q_full` and `mu0` were added; its `t2_*`, `first_*` and `j2` are identical [own].
  - Rerun with the experiment's own functions and arguments, these are bit-identical to the raw [own]:
    - U at (Δ, s, c̄) = (0.05, 0.03, 0.049), the swapping tie cell, and at (0.0175, 0.003, 0);
    - the E S job at (0.0175, 0.003, 0.045) (1400 draws);
    - F for (0.8, 0.6) under H1 at (c̄, s, Δ) = (0.08, 0.03, 0.05), and the (0.8, −0.4) null at (0.055, 0.003, 0.014).
- **Numbers.**
  - D52.md has 547 decimal numbers. 541 appear verbatim in the md, 5 are roundings of md numbers (0.563, 1.053, 32.3, 0.0254, 0.88), and "90.7%" is the md's 0.907 [own].
  - I checked every round-71 number in context (l. 55–57, 62–64, 70, 221–226, 234–235, 240, 253–256) [md].
  - The derived numbers follow [own]:
    - "≈13% in samples": 1.032⁴ = 1.134;
    - "4–6%": 1.043, and 1.056/1.003 = 1.053.
- **D3br-1 is fixed, and the zero at the tie belongs to the cells, not just to the shared samples.**
  - `q_full` equals E's `q_law` in 15 of 15 cells, and `q2` = `q_full` in the 6 cells at c̄ ≥ 0.045 [raw]. The part is +0.0000 there [md].
  - Population check at those six cells, with 400 000 law samples and 40 000 fresh H1 draws per cell [own]:
    - no sample or draw has the full knot below the threshold while the pair's own knot is above it (0 of 2 400 000; 0 of 240 000);
    - the tie pair sets 100% of the law's 5% tail.
  - So the part is 0 to ≲10⁻⁴ (95% bounds ≈7.5·10⁻⁶ under H0 and 7.5·10⁻⁵ under H1, per cell) [own]. The shared samples only make the zero exact.
- **"The full knot is never later than the pair's own" is true draw by draw.**
  - **Proof.** U's own knot is (1 − sign(v_k)u_k)/|v_k|, built from the path's duals u and rates v, for k = j₂ (or j₁ where j₂ enters first). The path's knot is the minimum of the same expression over every kernel except the noisy first entrant, and k is never that entrant. So full ≤ own.
  - **Law samples:** `knot_of` ≤ `knot_of_kernel` holds exactly in all 6 000 000 samples [own].
  - **Paths:** of 21 000 draws, dt exceeds t2 in 3293, by at most 6.5·10⁻¹² relative. All are draws where the own knot is the minimum, so this is round-off [own]. This also confirms that U is paired with E draw by draw.
- **The equal-size comparison (D3br-3) reproduces and is robust.**
  - **Every printed figure reproduces** in an independent re-implementation [own]:
    - (0.8, 0.6): α′ = 0.055 (size 0.0476 against z's 0.0478), Δ₅₀/NP 1.045 [0.973, 1.226], paired 1.032 [0.973, 1.239], ahead in 4 of 24;
    - (0.8, −0.4): α′ = 0.051, Δ₅₀/NP 1.679, paired 1.652 [1.598, 1.798], ahead in 0 of 24.
  - **The grid's resolution is enough.** The p-values are multiples of 1/M_E = 0.0005, and so is the grid, so no finer α′ changes a rejection set [own].
  - **The size estimate's noise does not change the conclusion** [own]:
    - the paired SE of the two sizes' difference is 0.0019;
    - a 1000-replicate paired bootstrap that re-selects α′ gives α′ in 0.0500–0.0595 and a paired median in [1.015, 1.046], above 1 in every replicate, with the knot ahead in 1–8 cells;
    - for (0.8, −0.4) it gives [1.629, 1.693], with the knot ahead in no replicate.
  - **Other ways of matching the size agree** [own]:

    | α′ chosen | paired median | knot ahead |
    |---|---|---|
    | pooled per pair (as printed) | 1.032 | 4 of 24 |
    | per offset (α′ 0.043–0.066) | 1.041 | 4 of 24 |
    | per cell (α′ 0.033–0.079) | 1.035 | 2 of 24 |
    | the wall's size set to 0.05 (α′ = 0.058) | 1.024 | 4 of 24 |
    | z's threshold moved to the knot's 0.0446 | 1.037 | 3 of 24 |

  - **Where the knot leads.** In 4 cells, by at most 2.7%: (s, c̄) = (0.003, 0.045), (0.003, 0.049), (0.003, 0.055) and (0.01, 0.032) [own].
  - **The other wall is the right one-sided test for the opposite-sign pair.** The fit's-sign wall never reaches power ½ there [md].
- **E′'s margin column is right.**
  - argmax|Φᵀf₀| = argmax|Φᵀf₁| (kernel 10) in all 15 cells, so `mu0` is the null's margin whichever first entrant it uses [own].
  - The column matches Φ⁻¹ of the swap rate [own]:
    - Φ(−μ₀/sd) = 0.2321 against the swap rate 0.233 (μ₀/sd 0.73 against Φ⁻¹(0.767) = 0.729);
    - 0.0190 against 0.016 at the other swapping cell;
    - P(μ̂_u < 0) equals the swap rate in all 15 cells.
- **(iv)'s new remainder O(ϑ²(Δ/σ)² + (Δ/σ)⁴) is right.**
  - In the continuum (60 001 points), (|cos| − (1 + ϑ²/6)^{−1/2})/(ϑ²(Δ/σ)²) tends to −0.0564 for (0.8, 0.6) and to +0.109 for (0.8, −0.4). At ϑ = 0 the remainder is −6.510·10⁻⁴(Δ/σ)⁴ [own].
  - The constants are uniform in the amplitudes [own derivation]. For A₁A₂ < 0, (1 − 3A₁A₂/A²)(Δ/σ)² < ϑ² and ((A₁² + A₂²)/A²)(Δ/σ)² ≤ ϑ², so the fourth- and fifth-moment terms are O(ϑ⁴). Numerically, the remainder over ϑ⁴ is 0.012–0.014 for (0.8, −0.4) and (0.8, −0.6) [own].
- **The loss's first two parts hold.**
  - Off the tie, the own knot and z₂ decide alike on all but 0–2% of H1 draws: paired differences −0.012 … +0.012, SE ≤ 0.008 [own].
  - At c̄ = 0.049 the margin's noise costs −0.035, −0.140 and −0.348 (SE 0.013, 0.025, 0.027) [own].
- **Tests.**
  - Z02, D02, SM01, SM03 and SM05 are now killed by tests that execute the mutated code [own].
  - Eight of twelve `u_job` mutants (seed, q_full, mu0, M, α) are killed. Their reruns change 1–17 md lines [own].
- **Status.** The Liu–Nehorai/Smith sentence is now consistent. Liu–Nehorai's ARL comes from the GLRT's asymptotics and is ∝ √CRB of the separation [snippet].

**What is wrong (majors first)**

No major issue. Four minor issues and loose ends:

1. **D3bt-1 (minor). The zero at the tie is argued from an inequality that only bounds it from below** ((iii) l. 69–72; E′'s fixed sentence, exp l. 866–870).
   - With a common threshold q, "same α quantile" and "never later" give part = P(full < q ≤ own) ≥ 0. The minimum cannot lose; it could still gain.
   - The zero rests on a separate, empirical fact the script never computes: no H1 draw falls between its two knots at q (0 of 2400) [own].
   - Round 71's own report and log state the inequality the wrong way round.
2. **D3bt-2 (minor). The minimum's part off the tie is quoted to four decimals with no Monte Carlo error** (l. 253–254).
   - Per-cell SE: 0.013–0.021 from the 400 paired H1 draws, plus 0.004–0.008 from the law [own].
   - The population values are −0.086 … −0.034, against the md's −0.1225 … −0.0250 [own/md].
3. **D3bt-3 (minor). The round-71 tests execute the mutated code, but often on data where the mutant changes nothing.**
   - Of 69 single-point mutants the tests kill 26. Another 27 are caught only by the md byte diff [own], among them:
     - 6 of 9 FA and 3 of 5 SL mutants that the log credits to the new tests;
     - 7 of the 11 `e_prime_lines` mutants;
     - 3 of the 5 `f_lines` z-range mutants.
   - 16 survive both checks. 14 of them change behaviour on inputs this data does not contain [own].
4. **D3bt-4 (minor). (iii)'s new size wording turns one noisy offset into a pattern** (l. 54–57).
   - c̄ = 0.08 is 0.02 from kernel 11's centre, the mirror of c̄ = 0.02. Its size is 0.0413 [0.0337, 0.0500] [md], yet the text counts it among "below it from c̄ = 0.032 on".
   - For (0.8, 0.6) the mirror offsets agree (Fisher p = 0.34) [own].
5. **D3bt-5 (cosmetic). Loose ends** (see Issues).

## What I checked and how

**Reproducibility**
- `cmp` of the md that `--from-raw` regenerates into a scratch `CF_OUT`.
- `pytest -q -p no:cacheprovider tests/test_close_features.py`; output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/out/pytest_base.txt`.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/cmp_raw71.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/cmp_U71.py`: the raw against 42aa6a3's, recursively and exactly; U field by field.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/rerun_jobs.py`: two U jobs, one E S job and two F jobs with `part_U`'s, `part_S`'s and `f_jobs`' arguments, compared exactly. Outputs: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/out/rerun_S.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/out/rerun_F.txt`.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/numbers_check.py`: every decimal number in D52.md against the md. I checked the round-71 numbers by hand, in context.

**The minimum's part and the tie**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/tie_ineq.py`: E's path knot against U's own knot, draw by draw (21 000 draws), and the draws between the two knots at the threshold.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/pop_part.py` (output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/out/pop_part.txt`):
  - 20 law seeds × 20 000 samples per cell: the law's MC SD, and the inequality `knot_of` ≤ `knot_of_kernel` asserted on every sample;
  - the McNemar SE of E's 400 H1 draws;
  - the population part, from 400 000 law samples and 40 000 fresh H1 draws per cell. These use a vectorised `one_fit` and the knot formulas, validated on E's 1400 draws per cell: centres to 1.5·10⁻¹⁵ and knots to 5·10⁻⁷ relative.
- The knot residual: the path's amplitude differs from the least-squares one by ≤ 2.8·10⁻¹⁴, and the knots are that sensitive (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/w_check.py`).
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/tail_who.py`: which kernel sets the H0 law's 5% tail.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/own_vs_z2.py`: the own knot against z₂ on the same draws, with the paired SE.

**Equal size**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/equal_size_check.py` (output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/out/equal_size.txt`):
  - an independent `equal_size_line`, with the size curve and the H0 overlap of the two tests;
  - α′ per offset, per cell, and for size 0.05; z moved to the knot's size;
  - a 1000-replicate paired bootstrap that re-selects α′.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/equal_size_cells.py`: the cells where the knot leads.

**The margin, the size, the cosine**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/margin_check.py`: `mu0` with j₁ from f₁ or from f₀; Φ(−μ/sd) against the swap rate; P(μ̂_u < 0).
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/size_offsets.py`: the calibrated knot's size per offset, per pair, and by distance to the nearest kernel; the mirror offsets (Fisher).
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/scripts/cos_check.py`: (iv)'s |cos| in the continuum for (0.8, 0.8), (0.8, 0.6), (0.8, −0.4) and (0.8, −0.6), at Δ = 0.05 … 0.00625.
- By hand: the inequality; the identity part = P(full < q ≤ own) at a common q; the remainder's amplitude factors.

**Mutation**
- **Definitions:** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/mutants.py`, 69 single-point mutants:
  - 22 of round 71's, ported verbatim from its definitions: Z02, D02, SM01–SM05, FA01–FA09, SL01–SL05, A06;
  - 47 new ones: `u_job` (12); `e_prime_lines`' margin column and shared-samples summary (11); `f_lines`' z range (5); `equal_size_line` and its call (19).
- **Runner:** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/run_mut.py` runs each mutant against:
  - the full test file;
  - `--from-raw`, with a byte diff of the md;
  - for `u_job`, `s_job`, `score_stats` and `draw_stats` mutants, also a rerun of the 15 U jobs or the 15 E S jobs, with the md regenerated and diffed (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/behav.py`).
  
  The identity mutant passes all three checks [own].
- **Results:** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/results_a.jsonl`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/results_b.jsonl`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/results_c.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref73/D3bt_work/mut/results_d.jsonl`. Outcome [own]:

  | code | mutants | killed by the tests | caught only by the md byte diff | survive both |
  |---|---|---|---|---|
  | round 71's, ported | 22 | 10: Z02, D02, SM01, SM03, SM05, FA02, FA08 (by a KeyError), FA09, SL01, SL04 | 9: FA01, FA03–FA07, SL02, SL03, SL05 | 3: SM02, SM04, A06 |
  | `u_job` | 12 | 8: U71-01, -02, -03, -07, -08, -09, -11, -12 | — (`--from-raw` does not rerun U) | 4: U71-04, -05, -06, -10 |
  | `e_prime_lines` (margin, summary) | 11 | 1: EM08 (by a KeyError) | 7: EM01–EM03, EM07, EM09–EM11 | 3: EM04, EM05, EM06 |
  | `f_lines`' z range | 5 | 0 | 3: FR02, FR04, FR05 | 2: FR01, FR03 |
  | `equal_size_line` and its call | 19 | 7: ES01, ES02, ES06, ES07, ES08, ES11, ES19 | 8: ES04, ES09, ES10, ES14–ES18 | 4: ES03, ES05, ES12, ES13 |

  - **Reruns** [own]: Z02, D02 and SM05 change 16, 15 and 15 md lines, and SM01 changes 2. SM03 changes the raw only. The killed U71 mutants change 1–17 lines; U71-04, -05, -06 and -10 change none.
- **Literature:** web searches (snippets only).

## Verdicts per part

| claim (D52.md line) | verdict | evidence |
|---|---|---|
| Title: the NP bound with the pair's sign known; ≈ 1.09 (1.04–1.18) without it (l. 1–4) | **holds** (unchanged) | z one-sided 1.003 and 1.011; \|z\| 1.085 and 1.094, range 1.039–1.179 [md] |
| (i) second entry, Δ⁻², opposite wall, dipole (l. 8–21) | **holds** (unchanged) | A identical to 42aa6a3's raw [own] |
| (ii) stiffnesses, rates, two regimes (l. 23–41) | **holds** (unchanged) | B identical [own] |
| (iii) the knot's law, LARS = lasso for k ≤ 2, the null across the offset (l. 43–54) | **holds** (unchanged) | C and E identical; the S and F reruns bit-identical [own] |
| (iii) calibrated size: "nominal at and next to a kernel's centre (c̄ = 0, 0.02…), below it from c̄ = 0.032 on (0.0329–0.0413)" (l. 54–57) | **numbers right; the grouping fails at c̄ = 0.08** | c̄ = 0.08 sits 0.02 from kernel 11, with 0.0413 [0.0337, 0.0500] [md]. Mirror offsets agree: p = 0.34; halves p = 0.15 [own] (D3bt-4) |
| (iii) de facto one-sided, 75–100% (l. 57–59) | **holds** (unchanged) | [md] |
| (iii) one-sided: 24 of 24 at the nominal level (0.045 vs 0.048); at equal size 1.045, 3.2% (≈13% in samples), ahead in 4 of 24 (l. 61–63) | **holds, robust** | Re-implemented [own]. Bootstrap [1.015, 1.046], ahead 1–8. Other matchings 1.024–1.041. 1.032⁴ = 1.134 [own] |
| (iii) two-sided; the allocated score; the lead over \|z\| only through the allocation (l. 63–66) | **holds** (unchanged) | F identical [own] |
| (iii) loss parts 1–2: the filter, the margin's noise (l. 66–69) | **holds** | Own knot − z₂ off the tie −0.012 … +0.012 (paired SE ≤ 0.008). At c̄ = 0.049: −0.035, −0.140, −0.348 [own] |
| (iii) part 3 off the tie: a loss "where far kernels' knots share the null's tail" (l. 69–70) | **holds; range noisy; "far" partial** | Population −0.086 … −0.034 [own] against −0.1225 … −0.0250 [md]. The fit's other neighbour sets 9–19% of the tail [own] (D3bt-2) |
| (iii) part 3 at the tie: "none … exactly 0 in 6 of 6 … same α quantile … never later" (l. 70–72) | **holds; the stated reason gives only ≥ 0** | q2 = q_full [raw]. 0 of 2400 H1 draws between the knots; population 0 [own] (D3bt-1) |
| (iii) "the full knot is never later than the pair's own" (l. 71–72) | **holds** (proved; data agree to round-off) | ≤ 6.5·10⁻¹² relative, in 3293 of 21 000 draws [own] |
| (iii) opposite-sign bias; walls at α/2; the moving null (l. 72–78) | **holds** (unchanged) | F identical [own] |
| (iv) z's law, mean and cosine; (a), (b); the proof's remainder (l. 80–92, 128–137) | **holds** | Continuum ratios → −0.0564 and +0.109; −6.510·10⁻⁴(Δ/σ)⁴ at ϑ = 0; uniform in the amplitudes [own] |
| (v) placement orders (l. 94–99) | **holds** (unchanged) | C (D's draws) identical [own] |
| Status: (iii)'s scope (l. 143–146) | **holds; wording** | "the α allocation" appears twice (D3bt-5) |
| Status: Liu–Nehorai / Smith (l. 148–154) | **holds in substance** | LN's ARL ∝ √CRB [snippet]; constants [md]; LN's amplitude model [memory] |
| New (l. 167–174) | **holds as revised** | ≈3% at equal size and "none at it" robust [own]; novelty plausible [memory] |
| Numerics (l. 176–256) | **match the md** | 541 verbatim, 5 roundings, 90.7% [own]. "−0.0075 … +0.0175" (l. 255) is round 67's md (D3bt-5); l. 253–254 has no errors (D3bt-2) |

## Status of the round-71 issues

| D3br | status | evidence |
|---|---|---|
| 1 (major) | **fixed** | `u_job` draws its law from E's samples (exp l. 1024). q_full = q_law in 15 of 15 cells and q2 = q_full in 6 of 6 tie cells [raw]. The part is 0 there [md], and 0 in the population [own]. Residue: the reason given (D3bt-1), inherited from round 71's own wording (REVIEW_round71_D3br.md l. 217; REVIEW_LOG l. 2812). The off-tie values carry no errors (D3bt-2) |
| 2 (minor) | **partly fixed** | Now killed: Z02, D02, SM01, SM03, SM05, FA02, FA08 (a crash), FA09, SL01, SL04 [own]. Still caught only by the md byte diff: FA01, FA03–FA07, SL02, SL03, SL05. Still surviving both: SM02, SM04, A06. The requested "pair of each sign" for `mu_u` was not added (D3bt-3) |
| 3 (minor) | **fixed** | `equal_size_line`; (iii) l. 61–63 and New l. 169–170. Reproduced and robust [own] |
| 4 (cosmetic) | **fixed** | Status l. 148–154 |
| 5 (cosmetic) | **mostly fixed** | Fixed: sign(ŵ)z (l. 64); z's range (l. 256); the two parts labelled (l. 222, 253); the null's margin (E′); the remainder (l. 132); the Status's scope (l. 145). Round 64's log row is corrected by an erratum (REVIEW_LOG l. 2807; l. 2576 unedited). Residue: the new "nominal at 0.02" wording over-reads one offset, as D3br-5 itself did (D3bt-4) |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bt-1 | minor | **The tie's zero is argued from an inequality that only bounds it from below.**<br>**Where:** (iii) l. 69–72; E′'s summary, a fixed string (exp l. 866–870; md l. 210); round 71's report ("equal tails leave no room for a gain", REVIEW_round71_D3br.md l. 217) and log ("as the inequality … requires", REVIEW_LOG l. 2812).<br>**The logic.** At a common threshold q, full ≤ own draw by draw gives part = P(full < q) − P(own < q) = P(full < q ≤ own) ≥ 0. The minimum cannot lose; it can gain. The zero needs a further fact: no H1 draw has its full knot below q while the pair's own knot is above it.<br>**That fact holds** [own]:<br>– 0 of 2400 E draws at the six cells;<br>– 0 of 240 000 fresh H1 draws and 0 of 2 400 000 law samples;<br>– off the tie the same event has probability 0.100–0.208 (c̄ ≤ 0.02) and 0.036–0.048 (c̄ = 0.032) under H1, outweighed there by the lower threshold.<br>**The script never computes the inequality.** Its sentence would print unchanged if a rerun produced a gain. The inequality holds exactly (proof in the Summary), and in the paths up to round-off (≤ 6.5·10⁻¹², 3293 of 21 000 draws) [own]. | **E′ (exp l. 866–870):** print, per tie cell, the H1 draws with dt < q ≤ t2 (0 in 6 of 6) and max(dt/t2 − 1) (6.5·10⁻¹²). Reword: "never later, so the minimum cannot lose; it gains nothing because no H1 draw falls between its two knots at the threshold".<br>**(iii) l. 70–72:** the same, with "0 of 2400".<br>**REVIEW_LOG l. 2812:** "bounds it from below". |
| D3bt-2 | minor | **The minimum's part off the tie is quoted without its Monte Carlo error, and its range is inflated.**<br>**Where:** Numerics l. 253–254 ("−0.1225 … −0.0250"); E′'s per-cell list (md l. 210).<br>**Errors per cell** [own]: E's 400 paired H1 draws give SE 0.013–0.021 (McNemar). The 20 000-sample law adds an SD of 0.004–0.008 (20 seeds).<br>**Population** [own] (400 000 law samples, 40 000 fresh H1 draws per cell):<br>– s = 0.003: −0.0845, −0.0822, −0.0384 at c̄ = 0, 0.02, 0.032;<br>– s = 0.01: −0.0716, −0.0687, −0.0341;<br>– s = 0.03: −0.0856, −0.0832, −0.0357.<br>That is −0.086 … −0.034 against the md's −0.1225 … −0.0250 [md]. The −0.1225 is 2.0 SE from its cell's −0.0845. The loss is real in all nine cells.<br>**Who shares the tail** (l. 69) [own]: in the H0 law's 5% tail, far kernels set 25–63% and the fit's other neighbour 9–19%. At c̄ = 0 that neighbour (kernel 11, 0.16–0.19) is as frequent as the runner-up (kernel 9, 0.18–0.21). | **E′:** print each cell's paired SE beside the part.<br>**Numerics l. 253–254:** "a loss of 0.03–0.09 (SE ≈ 0.02 per cell)".<br>**(iii) l. 69:** "where other kernels — far ones and the fit's other neighbour — share the null's tail". |
| D3bt-3 | minor | **The round-71 tests run the mutated code, but often on data where the mutant changes nothing** (outcome table in "What I checked and how").<br>**What still escapes** [own]:<br>– of the 22 ported mutants, the tests kill only 10;<br>– FA01, FA03–FA07, SL02, SL03 and SL05 are caught only by the md byte diff;<br>– FA08 and EM08 are "killed" only by a KeyError on rows that lack `q_law` or `q_full`;<br>– 10 of the 11 `e_prime_lines` mutants and all 5 z-range mutants escape the tests.<br>**Why** [own]:<br>– `test_f_alloc_lines_on_synthetic_rows` (tests l. 893–910) uses one set of rows for both pairs and both offsets, qA = qB, no H1 rows and no `q_law`, and asserts only the α = 0.05 per-offset figure;<br>– `test_smith_line_…` (l. 913–921) monkeypatches `crb_sep` (so SL05 never runs) and uses exact integer powers;<br>– `test_s_job_margin_…` runs (0.8, 0.6) only;<br>– no test passes U rows with `q_full` or `mu0`, or asserts the range or the equal-size line in `f_lines`.<br>**Behaviour-changing survivors of both checks** (14; equivalent only on this data):<br>– SM02, SM04, U71-04, U71-05, U71-06 (wall −1, ŵ < 0, or a null whose own first entrant differs);<br>– A06, ES03, FR03 (ŵ < 0);<br>– FR01 (an H0 extreme);<br>– ES05 (α′ < 0.05), ES12 (a cell where one Δ₅₀ crosses), ES13 (an incomplete cell);<br>– EM05, EM06 (U's samples not E's).<br>EM04 is equivalent in general, since q2 ≥ q_full always. U71-10 differs only when a third kernel leads: 0 of 300 000 law samples [own]. | **Tests:**<br>– assert `mu_u` in the (−0.8, −0.6) loop of `test_s_job_on_an_e_cell_and_a_c_null` (l. 777–806, which asserts only `mu` at l. 802);<br>– run `u_job` with PAIR_C = (−0.8, −0.6) for `mu0`;<br>– `f_alloc_lines`: rows that differ between pairs and offsets, H1 rows, qA ≠ qB, a `q_law`, and the α = 0.01 figures asserted;<br>– `smith_line`: a non-integer power, and the unpatched `crb_sep` on the orthogonal design of `test_crb_and_np_power_closed_forms`;<br>– `e_prime_lines`: U rows with `q_full`, `mu0`, and a cell where q2 ≠ q_full;<br>– `f_lines`: ŵ of both signs, and the range and equal-size lines asserted;<br>– `equal_size_line`: α′ < 0.05, one non-crossing cell, one incomplete cell. |
| D3bt-4 | minor | **(iii)'s new size wording turns one noisy offset into a pattern** (l. 54–57).<br>**The claim:** "nominal at and next to a kernel's centre (c̄ = 0, 0.02 …), below it from c̄ = 0.032 on (0.0329–0.0413)".<br>**The cell is mirror-symmetric** about the tie at 0.05: kernels at 0 and 0.1, and the samples map onto samples. So c̄ = 0.08 is 0.02 from kernel 11's centre, the mirror of c̄ = 0.02; for (0.8, 0.6) the null sits at c̄ to 10⁻⁴ [own]. Its size is 0.0413 [0.0337, 0.0500] [md], and the text counts it under "below".<br>**The mirrors agree** for (0.8, 0.6) [own]: 62/1200 against 51/1200 at 0.02/0.08 (Fisher p = 0.34); near half 0.0483 against mirrored half 0.0411 (p = 0.15).<br>**By distance to the nearest kernel** (both pairs, mirrors pooled) [own]:<br>– centre: 0.0512 [0.0428, 0.0608];<br>– 0.02: 0.0456 [0.0399, 0.0519];<br>– 0.032: 0.0398 [0.0344, 0.0457];<br>– 0.045: 0.0379 [0.0327, 0.0437];<br>– 0.049: 0.0329 [0.0261, 0.0409].<br>**"The more so at α = 0.01"** fails at 0.02 from the centre (0.0094: 6% below, against 9% at α = 0.05) and at 0.049 (33% against 34%) [own]. | **(iii) l. 54–57:** "nominal at a kernel's centre (0.0512), 0.046 [0.040, 0.052] at 0.02 from it (c̄ = 0.02 and its mirror 0.08), below from 0.03 from it on (0.033–0.040), lowest next to the tie".<br>**F:** print the size by distance to the nearest kernel. |
| D3bt-5 | cosmetic | **Loose ends.**<br>– Status l. 145–146: "the α allocation" appears twice.<br>– Numerics l. 255: "−0.0075 … +0.0175" are round 67's numbers (42aa6a3's md); no current script prints them.<br>– `u_job` takes j₁ from f₁ for the null's margin (exp l. 1023, 1028). argmax\|Φᵀf₀\| is the same kernel in all 15 cells [own], so nothing changes.<br>– `e_prime_lines` falls back silently to H1's `mu_true` without U (exp l. 820), under a header that says "at the null f₀" (exp l. 797).<br>– Status l. 154, "the score, the GLRT's local form": the LRT of one against two components is a supremum of a squared score process (Dacunha-Castelle–Gassiat 1999; Azaïs–Gassiat–Mercadier 2006 [memory/snippet]). The curvature score is its limit along ∂²_cφ.<br>– Not cited [snippet/memory]: Amar–Weiss 2008 (SNR ∝ δ⁻⁴); El Korso–Boyer–Renaux–Marcos (GLRT-based statistical resolution limits); Azaïs–Gassiat–Mercadier 2006 (the mixture LRT's local power: the mixture form of Δ₅₀ ∝ N^{−1/4}). | **Status l. 145–146, 154; Numerics l. 255:** text.<br>**`u_job`:** j₁ from f₀.<br>**`e_prime_lines`:** print "—" without U. |

## Sources

- **Read in the export:**
  - in full: D52.md, `exp_close_features.py`, close_features.md, `test_close_features.py`, REVIEW_round71_D3br.md, and the Round 67 and 71 entries of REVIEW_LOG.md;
  - the diffs 42aa6a3..1f51d48, and `rich_transit.py` l. 30–49.
- **Round 71's scratch mutant definitions** (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/mutants.py`), ported verbatim.
- **Web, snippets only.** arxiv.org and ams.org are unreachable, and fetching profiles.wustl.edu was refused (host not allowed in this session).
  - [Liu–Nehorai 2007 abstract](https://profiles.wustl.edu/en/publications/statistical-angular-resolution-limit-for-point-sources/), via a search snippet: the ARL from the GLRT's asymptotics, ∝ √CRB of the separation.
  - El Korso et al.'s GLRT-based limits from a first-order Taylor model: [hal-00771379](https://hal-supelec.archives-ouvertes.fr/hal-00771379), [hal-00771402](https://hal-centralesupelec.archives-ouvertes.fr/hal-00771402).
  - Amar–Weiss, IEEE TSP 56(11) 5309–5318 (2008), citation only: [reference list](https://hal-supelec.archives-ouvertes.fr/hal-01103527v1/html_references).
  - Azaïs–Gassiat–Mercadier, Bernoulli 12(5) 775–799 (2006): [projecteuclid](https://projecteuclid.org/euclid.bj/1161614946).
- **[memory]:** Harris 1964; Helstrom 1973; Shahram–Milanfar 2004–06; Smith 2005; Neyman–Scott 1966; Kulaitis–Munk–Werner 2021; Fukumizu–Amari 2000; Liu–Wu–Wang 2019; Tibshirani–Taylor–Lockhart–Tibshirani 2016; Azaïs–De Castro–Mourareau 2018; McLachlan 1987; Liu–Zhang; Batenkov–Goldman–Yomdin; Denoyelle–Duval–Peyré 2017; Dacunha-Castelle–Gassiat 1999; Amar–Weiss's δ⁻⁴ law. Whether Liu–Nehorai's model fixes or knows the amplitudes also rests on memory.
- **Assessment.** On these sources the classical attributions look right. "New" is plausible for the delay law, the knot's offset dependence, its allocation, the opposite-sign bias and the loss decomposition.

## Lesson

An inequality settles one side only. "The full knot is never later than the pair's own", together with equal quantiles, shows that the dictionary's minimum cannot lose at the tie. That it gains nothing is a separate, empirical fact about the H1 draws: none falls between its two knots at the threshold. The script should count that fact, not print it as a consequence; round 71 got the zero right and the reason backwards, and the text inherited the reason. The rest of this round is about noise. A Monte Carlo difference quoted to four decimals needs its standard error before its extremes are read: the part's −0.1225 is two standard errors from its cell's population value. On a grid with mirror symmetry, every offset has a twin with an equal claim to describe "next to a kernel's centre": c̄ = 0.02 reads 0.0500 and its mirror c̄ = 0.08 reads 0.0413, and the two agree within noise. Finally, a test that runs the mutated function kills a mutant only if its data reach the changed branch. Identical rows for both pairs and offsets, a single sign, exact integer powers or a monkeypatched helper make most mutants invisible to the test written for them.
