# Referee report D3eb, round 133: Proposition D‴.59 as revised after round 130

## Scope

- **Export and work directory.** I exported commit 58123a7 with `git -C … archive 58123a7 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref133. W/repo/venv is a symlink to the repository's venv (Python 3.13.18, numpy 2.5.4, scipy 1.18.1; x87 long double).
- **Settings.** MW_OUT=W/out and MW_WORKERS=1. Every run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, and only one of my heavy processes ran at a time.
- **Process notes.**
  - I killed no process. I stopped only two of my own Monitor watchers.
  - My first mutant batch was stopped by the 30-minute background time limit after E15. I resumed it from E16 with `--only`, in a fresh slot.
  - Other sessions' processes (ref134; a noise_resonance mutant run) were left alone.
- **The repository is untouched; I used only read-only git.**
  - HEAD moved to c551c89 during the round through other sessions' commits. None of them touches the D59 files: `git diff 58123a7 -- experiments/exp_moving_width.py tests/test_moving_width.py results/chaos/moving_width* review/drafts/D59.md` is empty.
  - The uncommitted manuscript/paper18/paper18.tex in the working tree is another session's.
- **Read in full.**
  - review/drafts/D59.md (461 lines).
  - review/REVIEW_round130_D3dy.md and the Round 130 entry of REVIEW_LOG.md.
  - `git diff 967a0c9 e851e9b` for exp_moving_width.py, test_moving_width.py, moving_width.md and D59.md.
  - tools/mutants/round130_D59_own.py and tools/mutants.py.
  - What the changed code calls: `a2_exit`, `a2_exit_ld`, `_LDBumps`, `BumpsCW.corr`, `critical_rhs_nd`, `critical_jac_nd`, `a2_line_lam`, `a2_tangent`, `a2p_job`, `part_A2P`, `part_A2H`, `a2p_h0`, `a2p_stats`, `a2p_line`, `a2p_shift`, `a2h_line`, `nh_verdict`, and the verdict blocks of `a2_lines`.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my computation.
  - [snippet]: a web-search result I saw (no page fetched).
  - [memory]: not verified this round.

## Summary

**Design.**
1. **Reproduction.**
   - The md byte for byte from `--from-raw`, and the 57 tests.
   - A2P and A2H rows recomputed with the repository's job: rows 0–3 and row 100 of several records.
2. **The statistics, recomputed independently from the raw.**
   - Bootstrap intervals, the heterogeneity between cells, correlations across cells, and an order-blocked sem.
   - The three hypotheses for the 3A rtol 10⁻¹⁰ cell, tested with new runs:
     - **chance:** 400 more fresh orders (seeds 133–136);
     - **atol:** the double run at the long double's atol, 1e-14·σ*, and the long-double run at the double's, 1e-14;
     - **prefactor at equal steps:** a replication of A2H's 5A h₀/16 arm with seed 132.
3. **`nh_verdict`** against the HPS rate conditions, including edge cases and the non-normality of the (α, σ) block.
4. **D59 against the md and raw.**
5. **Mutants.** The repository's E00–E30 and 13 of my own (H01–H13) plus a control. Each ran against `-k "round130 or round125 or round112"` (the md byte test deselected), with the tool's `--from-raw` comparison.
6. **Literature:** four web searches.

**Verdicts.**
- **Reproducibility is exact.**
  - md identical (sha256 d14fdf43…54c41a8; 31.7 s).
  - 57 passed in 239.6 s (one LSODA warning).
  - Sampled A2P and A2H rows: exits, directions and step counts bitwise equal (rows 0–3; row 100 of three records). ln|S| and Λ agree to ≤ 6.8·10⁻⁸; the line's grid depends on the job's largest exit.
- **Every printed statistic is correctly computed.**
  - Pooled exponent 1.0128 ± 0.0092, 1.40 sem from 1. The order-blocked sem is 0.0098 (1.31 sem), because the same orders recur across cells.
- **The cap is applied identically in both arithmetics, and h₀ is well defined.** "Sharing it in part" follows.
- **`nh_verdict` matches the HPS conditions** for the frozen linearization, but for one edge case.
- **D59's numbers are the md's.**
- **All 30 repository mutants are killed semantically.**

**But the reading of the new statistics does not hold.**
- **The 3A rtol 10⁻¹⁰ excess is not chance, nor a defect of the steps' atol or of the paired summation.** It is a real difference between the arithmetics' prefactors, resolved at 3.95 sem on 600 fresh draws.
- **At equal (capped) steps the 5A h₀/16 excess replicates** (pooled +5.1 sem).
- **Consequences:**
  - The 3-sem rule tests equality of the prefactors, which is false. It is "met" only at n = 200, and the code's own rule would withdraw "so the seed scales" with more draws.
  - "Λ compares the arithmetics at equal steps" does not make Λ = 11 ln 2.
  - The exponent's ±0.009 is statistical; the systematic spread across configurations is about ±0.04.
- The scientific claim, that the seed scales with the arithmetic's unit roundoff with O(1) prefactors, survives.

**Findings.**
- **D3eb-1 (minor, borderline major):** the 3A rtol 10⁻¹⁰ cell's excess is real (+3.95 sem on 600 draws, unchanged with matched atol). The 3-sem rule is the wrong test.
- **D3eb-2 (minor):** at equal steps the arithmetics' prefactors still differ (5A h₀/16: Λ 7.91 ± 0.056, +5.1 sem on 200 draws), so "so Λ compares the arithmetics at equal steps" is a non sequitur. The exponent's precision is statistical only.
- **D3eb-3 (cosmetic):** A2H's "per order" shifts are unpaired in effect (r ≈ 0). Brouwer's h^{1/2} is the reference law the shifts should be read against.
- **D3eb-4 (cosmetic):** `nh_verdict(0, tg)` with min tg > 0 returns (True, True). The verdict is the frozen Jacobian's in an adapted metric (the (α, σ) block is strongly non-normal).
- **D3eb-5 (cosmetic):** "extremes … of 200" are extremes over 800; the sem ignores that the same orders recur across cells.
- **D3eb-6 (minor):** test gaps. H07 (the rule at 6 sem) survives, so the threshold is unpinned from above. H01 is equivalent only while K_A2H = K_A2P.

## What I checked and how

Scripts are in W/own/:
- stats.py (raw statistics);
- recomp.py (row recomputation);
- exp133.py (atol and cap arms);
- more133.py (seeds 133–136, default and matched atol);
- nh133.py (verdict edges, log-norms);
- spec133.py (my mutants).

Data: W/own/*.json. Logs: W/logs/. Mutant records: W/mut/res_own.jsonl, res_r130.jsonl, res_r130b.jsonl.

1. **Reproducibility** [own].
   - `--from-raw` takes 31.7 s; `cmp` against the committed md: identical.
   - pytest: 57 passed in 239.6 s. The slowest tests are round106's seed test (84.9 s), a2e (56.5 s), the round87 md test (29.0 s), round125 (12.9 s) and round130 (12.6 s).
   - recomp.py (57 s): `a2p_job` at (3A, 10⁻¹⁰, seed 131, k = 3), (5A, 10⁻¹², 132, 3), (3A, 10⁻¹², 132, 2) and (5A, 10⁻¹⁰, 132, 2), and at both h₀/16 caps (k = 2).
     - Exits, directions and step counts are bitwise those of the raw.
     - max |Δ ln S, ΔΛ| = 6.8·10⁻⁸, from the line's grid (1.05·the job's largest exit).
     - Row 100 of (3A, 10⁻¹², 132), (5A, 10⁻¹⁰, 131) and A2H 5A h₀/16 is bitwise equal.
2. **Statistics** (stats.py, from the raw alone) [own].

   | cell | pooled Λ (t) | per seed (t) | bootstrap 99% | skew, kurtosis | ln(S/uσ*): d, ld |
   |---|---|---|---|---|---|
   | 3A, 10⁻¹⁰ | 8.005 ± 0.143 (+2.65) | +3.55, +0.27 (the seeds differ by 2.3σ) | 7.635 … 8.376 (excludes 7.625) | 0.00, 4.36 | +0.36 ± 0.09, −0.02 ± 0.11 |
   | 3A, 10⁻¹² | 7.756 ± 0.166 (+0.80) | | | | |
   | 5A, 10⁻¹⁰ | 7.513 ± 0.129 (−0.87) | | | | |
   | 5A, 10⁻¹² | 7.616 ± 0.115 (−0.07) | | | | |

   - Heterogeneity between cells: χ² = 7.44 on 3 dof (p ≈ 0.06).
   - Pooled over the 800: 7.7226 ± 0.0699. Stratified: identical. Order-blocked (the four cells of an order averaged): ± 0.0744, i.e. an exponent of 1.0128 ± 0.0098.
   - Correlations across cells of the same order: |r| ≤ 0.20. Double against long double: |r| ≤ 0.04.
   - Under the null, P(max |t| ≥ 3.55 of 8) = 0.25%, and P(max |t| ≥ 2.65 of 4) = 2.7%.
   - One of the 800 orders has the long double earlier (−0.48%, 3A, 10⁻¹⁰, seed 132). One has Λ < 0 (−0.19).
3. **The three hypotheses for the 3A rtol 10⁻¹⁰ cell** [own].
   - **(a) Chance.** Fresh seeds 133–136 (more133.py, the repository's `a2p_job`) give 7.785 ± 0.222, 7.970 ± 0.193, 8.055 ± 0.213 and 7.956 ± 0.226.
     - Pooled over 400: 7.941 ± 0.107 (+2.97 sem).
     - With the raw's 200: 7.963 ± 0.086 (+3.95 sem), an exponent of 1.044 ± 0.011.
     - Round 130's referee's seeds 130 and 2026 (7.63 ± 0.23, 8.11 ± 0.21) agree.
     - Prefactors over the 400: double +0.21 ± 0.08, long double −0.10 ± 0.08 in ln(S/uσ*).
     - **Not chance.**
   - **(b) atol.** The paired runs differ in atol: 10⁻¹⁴ in double against 10⁻¹⁴σ* (3·10⁻¹⁶) in long double. This is why the steps are 60 against 73.
     - With the double run at 10⁻¹⁴σ*: 7.748 … 8.164 on seeds 133–136 (pooled 7.959 ± 0.106, +3.16), and 7.97 ± 0.28 and 7.77 ± 0.20 on seeds 131 and 132.
     - With the long-double run at 10⁻¹⁴ (seed 131): 8.05 ± 0.20.
     - Matching atol makes the double take 68 steps, but it moves the mean step apart: 7.80 against 8.87·10⁻⁴, 12%.
     - **Not the atol.**
   - **(c) The order of summation of the paired runs.** Both runs use the same permutation. Both reduce with `rP.sum(0)/N` in the same code (`BumpsCW.corr` and `_LDBumps.corr`). f* is evaluated elementwise after the permutation, and the ±x terms are bitwise antisymmetric at c = 0.
     - **Not the order** (code-read).
   - **Conclusion:** a real arithmetic-dependent prefactor at DOP853's own steps in this cell, about e^{0.3}.
4. **A2H** (stats.py; exp133.py `cap`) [own], [raw].
   - **h₀ is well defined**: `a2p_h0` = the mean of s_e/n over both seeds' 200 double runs at rtol 10⁻¹⁰, at DOP853's own steps.
   - **The same max_step goes to `a2_exit` and to `a2_exit_ld`.** Under the cap every step sits at it: the double's mean steps are 0.990–0.999 of the cap; the long double's mean steps equal the double's to 0.2%.
   - **Replication at 5A h₀/16 with seed 132** (326 s): Λ = 7.825 ± 0.077 (+2.59). ln(S/uσ*) is −1.31 ± 0.07 in double and −1.51 ± 0.03 in long double.
     - With seed 131's 8.000 ± 0.081: 7.91 ± 0.056, **+5.1 sem**, an exponent of 1.038 ± 0.007.
   - **Shifts.** Own and capped runs of the same order are uncorrelated: r from −0.11 to +0.24, directions agreeing in 43–55 of 100. So "per order" pairing adds nothing; the shift is a difference of independent means.
     - The arithmetics' shift differences are −0.62 ± 0.36, −0.87 ± 0.23 (3A) and −0.11 ± 0.23, +0.43 ± 0.18 (5A). **"Sharing it in part" follows.**
     - Brouwer's injection law (S ∝ h^{1/2}) predicts −0.69 at h₀/4 and −1.39 at h₀/16. Observed: 5A −0.64 and −0.53 at h₀/4, −0.70 and −1.13 at h₀/16; 3A double −0.54 and −0.94, but the 3A long double +0.08 and −0.07.
5. **`nh_verdict`** (nh133.py) [own], [memory].
   - HPS for an unstable normal bundle: λ_u > k·max tg for k = 0 … r. For a stable one: λ_s < k·min tg.
   - The code returns (k = 1 holds, every r holds). That is right for lam > 0, and for lam < 0.
   - **Edge case:** `nh_verdict(0.0, [0.1, 0.2])` returns (True, True). A zero normal rate is not hyperbolic.
   - **Non-normality.** At A2E's states the (α, σ) block is strongly non-normal (off-diagonals −0.028 and −53.5 at 5A, repelling).
     - Its log-norm is +0.047 in the raw coordinates. In the scaled coordinates (α/α, σ/σ) it is +1.39·10⁴, above λ = 1.09·10⁴.
     - So the verdict is the eigenvalues' (an adapted metric, the frozen Jacobian). It is robust in the raw metric.
   - **The plants** in test_round130 flip each branch; the repository's E24–E28 are killed.
6. **D59 against the md and raw.**
   - Every changed number is verbatim the md's:
     - 397/402 of 800; medians 19–23%, IQRs within 15–27%, "one order 0.5% earlier" (true: one of 800);
     - 8.01 ± 0.14 … 7.62 ± 0.12; +2.7, +0.8, −0.9, −0.1; +3.6; 1.013 ± 0.009; geometric means 0.60–1.43;
     - shifts −1.13 … +0.08, 5 of 8 beyond 3 sem; Λ 7.71, 7.46, 7.46, 8.00 (+4.7); 0.979 … 1.049;
     - h₀ 8.71·10⁻⁴ and 6.85·10⁻⁵; the step counts; "apart by … relatively".
   - The Status labels the Remark numerical, which is right.
   - "Extremes −0.19 … 16.27 of 200" and "−0.5% … 47.3% of 200" are extremes over the four cells, i.e. of 800.
7. **Mutants** (W/mut) [own].
   - **The repository's spec:** E00 (control) passes. E01–E30 are all killed by the semantic tests; E02, E09, E10, E14–E17, E19–E23, E26 and E28–E30 also change the md. This matches the log.
   - **Mine:** H00 (control) passes.
     - Killed: H02, H03, H09 (also the md), H11 and H13 (also the md).
     - Survivors, all equivalent on the data:
       - H01 (A2H's base not truncated): equivalent while K_A2H = K_A2P.
       - H04 (the attracting "every r" at min tg ≥ 0) and H05 (lam = 0 counted as repelling): the edges are unpinned.
       - H06 (the largest per-seed offset by signed t): the largest is positive.
       - H08 (the rule at 2.7 sem): the data's largest is 2.65.
       - H10 (> 3 for ≥ 3): no |t| equals 3.
       - H12: an algebraic identity.
       - **H07 (the rule at 6 sem) survives.** The rule's threshold is pinned from below only, by the data, not by a plant.
     - None of the survivors changes a published number now. H07 would keep the conclusion on the more numerous draws where the rule at 3 sem withdraws it. H01 would change A2H's numbers if K_A2H < K_A2P.
8. **Literature:** four searches, no page fetched (Sources).

## Verdicts per part

| claim (D59 unless noted) | verdict | evidence |
|---|---|---|
| md byte for byte; the 57 tests; A2P and A2H rows | **holds** | [own] |
| "Λ … 8.01 ± 0.14, 7.76 ± 0.17, 7.51 ± 0.13, 7.62 ± 0.12 … +2.7, +0.8, −0.9, −0.1 sem; one seed's … +3.6" | **numbers hold** | [raw], [own] |
| "each cell's mean within 3 sem of it … round 125's rule … met by every cell's pooled mean" | **true of these 200, but the rule tests the wrong thing.** The 3A 10⁻¹⁰ excess is real (+3.95 sem on 600; +3.16 with matched atol). The pooled aggregation was chosen after the per-seed record failed (D3eb-1) | [own] |
| "over the 800 draws Λ/(11 ln 2) = 1.013 ± 0.009" | **correctly computed** (order-blocked ± 0.0098). Read as γ, it overreaches: the cells range from 0.985 to 1.044 ± 0.011 at own steps and from 0.979 to 1.038 ± 0.007 capped (D3eb-2) | [own] |
| "each arithmetic's seeds O(1) in units of uσ* at DOP853's steps (… 0.60–1.43)"; "the seed scales with the arithmetic's roundoff" | **holds** | [raw], [own] |
| "its prefactor moves by up to about e^{−1} with the step sequence (… −1.13 … +0.08)"; "sharing it in part" | **holds.** The arithmetics' shift differences: −0.87 ± 0.23 (3A h₀/16), +0.43 ± 0.18 (5A h₀/16) | [md], [own] |
| "— at equal steps Λ/(11 ln 2) is 0.979 … 1.049 —, so Λ compares the arithmetics at equal steps" | **non sequitur.** At equal steps the prefactors differ (5A h₀/16: +5.1 sem replicated). At own steps the sequences are not equal (60/73 steps; atol differs) (D3eb-2) | [own] |
| A2H: the cap identical in both arithmetics; h₀ | **holds**: the same max_step; every step at the cap; h₀ the double's mean step over the 200 | code, [raw] |
| "ln(S/uσ*) moves per order by …" | **numbers hold; "per order" overstates the pairing** (r ≈ 0) (D3eb-3) | [own] |
| `nh_verdict`; "normally hyperbolic, for every r" (repelling), "not normally hyperbolic" (attracting) | **holds** for the frozen Jacobian in an adapted metric. The lam = 0 edge is wrong (D3eb-4) | [own], [memory] |
| "extremes −0.19 … 16.27 of 200", "−0.5% … 47.3% of 200" | **over 800** (D3eb-5) | [raw] |
| "(apart by 2.1·10⁻⁶ relatively)" | **holds**, tested against live exits | [md] |
| proved/numerical labels | **unchanged and correct** | re-read |
| tests | **good**: 30/30 repository mutants killed semantically. The rule's threshold is unpinned from above (H07) (D3eb-6) | [own] |

## Status of the round-130 issues

| D3dy | status | what remains |
|---|---|---|
| 1 (replicated draws) | **applied**: seeds 131 and 132, 100 orders each, the code's order apart; the pooled exponent printed | The +2.7/+3.6 is real (D3eb-1); the rule's aggregation was chosen post hoc |
| 2 (the steps) | **applied**: A2H, step counts, "at DOP853's steps" | "Λ compares … at equal steps" does not follow; the prefactor differs at equal steps too (D3eb-2) |
| 3 (extremes) | **applied**: medians, IQRs, "extremes … of N"; the lateness test on the median | "of 200" for extremes over 800 (D3eb-5) |
| 4 (hard-coded verdicts) | **applied**: `nh_verdict`; plants flip each branch | The lam = 0 edge; the frozen-metric qualifier (D3eb-4) |
| 5 ("apart by") | **resolved** | — |
| 6 (regeneration untested) | **resolved**: the configuration ties; the jobs at (3A, 10⁻¹²) and at the cap; G02, G15, G16, G09 and G12 now killed | H07, H01 (D3eb-6) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3eb-1 | minor (borderline major) | **The 3A rtol 10⁻¹⁰ excess is a real prefactor difference, and the 3-sem rule is the wrong test.** <br>• Fresh seeds 133–136: 7.941 ± 0.107 (+2.97). With the raw's: 7.963 ± 0.086 (+3.95, n = 600) [own]. Round 130's referee's seeds agree. <br>• Not atol: the double at the long double's atol gives 7.959 ± 0.106 (+3.16). Not the paired summation: same permutation, same reduction (code). <br>• The prefactors: double +0.21 ± 0.08, long double −0.10 ± 0.08. <br>• The rule ("within 3 sem of 11 ln 2") tests equal prefactors in the two arithmetics, which is false here. It passes at n = 200 and would withdraw "so the seed scales" with more draws. <br>• Its application to the pooled cell rather than to each seed's record was made after the per-seed record failed (+3.55; P(max of 8) = 0.25% under the null). | **Code (`a2p_line`):** replace the gate by an interval statement. For example: every cell's exponent Λ/(11 ln 2) with its 99% interval inside [0.9, 1.1], and every geometric mean within 0.05–20 uσ*. Print each cell's exponent. **D59 Remark and A2P:** "each cell's exponent 0.985–1.05; at 3A, rtol 10⁻¹⁰ the arithmetics' prefactors differ by about e^{0.3} (+3.9 sem on 600 draws) — the seed scales with u, its prefactor is not the same in the two arithmetics". Optionally K_A2P 200 at 3A, rtol 10⁻¹⁰, or a third seed. |
| D3eb-2 | minor | **At equal steps the prefactors still differ, so "so Λ compares the arithmetics at equal steps" does not follow, and the exponent's ±0.009 is statistical only.** <br>• 5A h₀/16: seed 131 8.000 ± 0.081 (+4.7); my seed 132 7.825 ± 0.077 (+2.6); pooled 7.91 ± 0.056 (+5.1), an exponent of 1.038 ± 0.007 [own]. <br>• Under the caps 3A is near 11 ln 2 (7.71, 7.46), so the steps move the arithmetics' prefactors differently by configuration. <br>• At DOP853's own steps the sequences differ (60 against 73 steps; atol 10⁻¹⁴ against 10⁻¹⁴σ*). <br>• Over the configurations tried the exponent spans 0.979–1.049. | **D59 Remark and A2H, and `a2h_line`'s closing clause:** "the exponent is 1 to within ±0.05 at every step sequence tried (pooled at DOP853's steps 1.013 ± 0.009, statistical), the two arithmetics' prefactors not equal even at equal steps (5A, h₀/16: +5 sem on 200 draws)". Drop "so Λ compares the arithmetics at equal steps". Optionally a second A2H seed. |
| D3eb-3 | cosmetic | **"Moves per order by …" suggests a paired precision that is absent.** Own and capped runs of an order are uncorrelated (r −0.11 … +0.24) [own]. Brouwer's h^{1/2} for derivative rounding predicts −0.69 and −1.39. 5A follows it at h₀/4 (−0.64, −0.53); 3A's long double does not move (+0.08, −0.07). | **D59 and md:** "shifts of the means (the orders' draws independent at different steps)". Cite Brouwer's law (S ∝ h^{1/2} for rounding in the derivative) as the reference law, and say where it fails (3A, long double). |
| D3eb-4 | cosmetic | **`nh_verdict`.** <br>• lam = 0 with min tg > 0 returns (True, True); H04 and H05 survive. <br>• The verdict compares eigenvalues of the frozen Jacobian. The (α, σ) block is strongly non-normal: its log-norm in (α/α, σ/σ) is 1.39·10⁴, above λ = 1.09·10⁴; in the raw coordinates +0.047. So the HPS reading holds in an adapted metric of the frozen linearization [own], [memory]. | **`nh_verdict`:** `if lam == 0: return False, False`, with a plant. **D59:** "(the frozen Jacobian's rates, in an adapted metric)". |
| D3eb-5 | cosmetic | **"Extremes −0.19 … 16.27 of 200" and "−0.5% … 47.3% of 200" are extremes over four cells of 200 each.** The pooled sem treats the 800 as independent, though each seed's orders recur in all four cells: 0.0098 order-blocked against 0.0092 [raw], [own]. | **D59:** "extremes over the 800". **`a2p_line`:** the pooled sem over orders (each order's four cells averaged), or a note. |
| D3eb-6 | minor | **Test gaps.** <br>• H07 (the rule at 6 sem) survives: the "+3 on one seed" plant gives its cell less than 6 sem, so the threshold is pinned from below only by the data. Under D3eb-1's data it decides whether "so the seed scales" is printed. <br>• H01 (A2H's base not truncated) is equivalent only while K_A2H = K_A2P. <br>• H06, H08, H10 and H12 are equivalent. | **Test:** plant a cell 3.5 sem off (shift one record so that the pooled t is 3.5) and assert no conclusion. Plant K_A2H < K_A2P in `a2h_line`. If D3eb-1's interval rule replaces the gate, plant both of its edges. |

## Sources

- [snippet] Lythe, "Noise-controlled dynamic bifurcations" ([adap-org/9707008](https://arxiv.org/abs/adap-org/9707008); also [adap-org/9707007](https://arxiv.org/abs/adap-org/9707007)): if the passage is slow, the dynamics are controlled by additive disturbances "even if they are extremely small"; multiplicative noise does not do this.
- [snippet] Berglund and Gentz, "Pathwise description of dynamic pitchfork bifurcations with additive noise" ([math.PR/0008208](https://arxiv.org/abs/math.PR/0008208); [abstract](https://www.idpoisson.fr/berglund/abs_pitch.html)): the delay and the exit window depend on the noise intensity, through √(ε|log σ|).
- [snippet] "Analytical prediction of delayed Hopf bifurcations in a simplified stochastic model of reed musical instruments" ([arXiv:2512.23269](https://arxiv.org/pdf/2512.23269)): in numerical integration, round-off errors act as an additive noise of very low level and change the observed bifurcation delay. This is the same mechanism as A2P's seed.
- [snippet] Neishtadt, maximal delay of stability loss ([arXiv:2606.07662](https://arxiv.org/pdf/2606.07662)): the deterministic background.
- [snippet] Brouwer's law, as described in [arXiv:1711.06050](https://arxiv.org/pdf/1711.06050): round-off accumulates as a Gaussian random walk with sd ∝ (t/h)^{1/2}.
- [snippet] Kehlet and Logg, a posteriori error analysis of round-off errors in ODE solutions ([Chalmers](https://research.chalmers.se/publication/251849)): the accumulated round-off is inversely proportional to √h. That is the state-rounding regime; A2P's seed, rounding in the derivative of a state near 0, has the opposite sign of exponent, h^{1/2}.
- [snippet] Mosbach and Turner, rounding errors as independent random variables in ODE solution ([arXiv:math/0512364](https://arxiv.org/abs/math/0512364)): near a saddle, a Gaussian-driven transition at a time of about −c·log h. This is the closest literature to A2H's step dependence.
- [snippet] Henrici's caveat that local rounding errors are not random variables, as quoted in a [MathGroup thread](https://forums.wolfram.com/mathgroup/archive/2011/Jul/msg00311.html); Kahan's criticism of probabilistic error analysis ([improber.pdf](https://www.cs.berkeley.edu/%7Ewkahan/improber.pdf)). Relevant to why the two arithmetics' prefactors need not agree.
- [memory] The Hirsch–Pugh–Shub r-normal-hyperbolicity rate conditions.

## Lesson

A significance test of "equal to the scaling's prediction" is not a test of the scaling.
- Λ = 11 ln 2 asks that the two arithmetics' prefactors be equal. With enough draws that fails: here at 3.95 sem on 600 draws at DOP853's steps and at 5.1 sem on 200 at equal capped steps.
- "The seed scales with u" holds throughout.
- State a scaling as an interval with its systematic spread, and replicate an outlier before calling it chance or a defect.
