# Referee report D3bn, round 67: Proposition D‴.52 as revised after round 64

## Scope

- **Commit reviewed.** 88933bf on branch claude/festive-bardeen-eopej2 ("D‴.52 round 64: the knot is de facto one-sided; both walls stay biased; the losses decomposed").
  - I exported it with `git archive 88933bf` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref67/D3bn_work/repo` and worked only there.
  - I extracted round 64's files (8c826c6) beside it for diffs.
- **Files reviewed.**
  - `review/drafts/D52.md` (236 lines).
  - `experiments/exp_close_features.py` (1313 lines: parts A, B, C, C′, D, E, E′/S, U and F, plus `--from-raw`).
  - `results/chaos/close_features.md` and `close_features_raw.npz`.
  - `tests/test_close_features.py` (33 tests).
  - `review/REVIEW_round64_D3bk.md` and the Round 64 entry of `review/REVIEW_LOG.md`.
  - For context: `src/core/rich_transit.py` (`sliding_iss`, `polish_block`), `src/core/mirror_free.py` (`iss_path`), the relative-mode law in `review/drafts/D54.md` and the width law in `D59.md`.
- **HEAD moved** to 908d816 while I worked. None of the files above differ between 88933bf and 908d816. I did not use the uncommitted edits other sessions had in the live tree.
- **How I ran things.**
  - The repo's venv with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`.
  - `CF_OUT` always pointed to a scratch copy of the raw.
  - At most two heavy processes of mine ran at a time; the machine's load from other sessions was 6–9. Nothing of mine is left running.
  - Every number is tagged [own] (my run), [md] or [raw].
- **I modified nothing in the repository.**

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (`cmp`), and the 33 tests pass in 19.5 s [own].
  - Reruns with the experiment's own functions and seeds are bit-identical to the raw [own]:
    - all 15 U jobs (`u_job` with `part_U`'s arguments);
    - all 95 S jobs (`part_S`'s job list);
    - 5 F jobs: (0.8, −0.4) at c̄ = 0.032, s = 0.01, Δ = 0.05 and 0.07 (the non-monotone cells); the (0.8, 0.6) null at (0.049, 0.01); (0.8, 0.6) at (0.049, 0.03, 0.07); (0.8, −0.4) at (0.045, 0.003, 0.04);
    - 3 A rows, among them (0.8, −0.4, 0.0125), which D3bk could reproduce only to 2.4·10⁻¹⁴.
  - Against round 64's raw [own]:
    - B, C, E, F and S are identical;
    - A differs only at round-off (rms ≤ 6·10⁻¹⁶, end coefficients ≤ 6.5·10⁻¹² relative);
    - U is new.
  - I also re-implemented U independently (my own runner-up, projector and law, fresh noise) [own]:
    - the same runner-up in all 15 cells;
    - per-draw statistics equal to 2·10⁻⁹;
    - sizes and powers within the law's Monte Carlo noise (quantile ratio 0.983–1.028; power sd 0.005–0.011).
- **The LARS identity (l. 46–49; proof l. 112–118) is right as now scoped.**
  - I re-derived u_j(τ) = τc_j(1/τ) piece by piece; the intercept ⟨φ_jΦ_A⟩⟨Φ_AᵀΦ_A⟩⁻¹s_A is fixed by continuity at each knot.
  - Test: 400 noisy targets on the experiment's dictionary (one bump, pairs of each sign, three bumps), `iss_path` against my own LARS and lasso homotopies [own]:
    - knots 1–2 agree to 8.8·10⁻¹⁴;
    - every knot inside the predicted range agrees to 2.3·10⁻¹², with 0 counterexamples;
    - every later difference (42, 112, 148 and 190 of 400 at knots 4–7) follows a failure of one of the two conditions;
    - after such a failure the next knot differs in 231 of 235 paths. In the other 4 the next entrant lies ≥ 10σ from the change.
  - My lasso satisfies KKT mid-piece up to its first drop (violations only after it, in 21 of 400 paths).
- **(iv) holds.**
  - For 8 pairs at 4 separations each, including (−0.8, −0.6) and (−0.8, 0.4), the noiseless z has V's sign and ŵz has sign(A₁A₂) [own].
  - I re-derived cos = (1 + ϑ²/6)^{−1/2} (the odd part's norm over the even part's is ϑ/√6).
  - Measured (1 − cos)/(ϑ²/12) is 1.007 for (0.8, 0.6) at Δ = 0.0125. For (0.8, −0.4) at ϑ = 0.31, measured over closed form is 0.986 [own].
  - I re-derived readings (a) and (b). D54's dD/ds = Λ₀(…) is the same law written in its own time.
- **The knot is de facto one-sided (l. 56–57).**
  - 390/430 and 328/364 [md]; 75–100% per offset [md].
  - The knot's law at the truth (20 000 samples) puts 80–100% of its 5% tail at the fit's-sign wall at every offset [own].
- **Both walls at α/2 stay biased (l. 66–68).**
  - Least power 0.010–0.045, below α in 24 of 24 cells; pooled 0.036 and 0.039 [md].
  - The walls test's own components under H1 [own]:
    - fit's-sign wall: 0.0054 and 0.0015 (0.0223 under H0);
    - other wall: 0.0310 and 0.0379 (0.0235 under H0).
  - At α = 0.01 the pooled power is 0.0069 and 0.0075 [own].
- **Part U is defined, computed and calibrated correctly.**
  - `knot_of_kernel` (l. 886–897) and `u_job` (l. 900–935) give the tie pair's second knot, exact in the sliding limit.
  - The convention's corner case, a third kernel entering first, occurs in none of the 21 000 draws [own].
  - The law is the first-order null at the truth's one-bump fit. Sizes are 0.039–0.068 [md]; z₂'s sizes, which the md does not print, are 0.040–0.064 [own].
  - The margin part is confirmed exactly. With the runner-up's margin frozen at its noiseless value, the own knot recovers z₂'s power in every cell: 0.705/0.705, 0.395/0.400 and 0.522/0.530 at c̄ = 0.049 [own].
- **The moving null (l. 68–70).**
  - ĉ − c̄ matches −T₃/(4Aσ²) = 104.2Δ³ for (0.8, −0.4) to within 5% for Δ ≤ 0.02 [own].
  - The quoted quantiles and medians are as printed [md].
- **Smith (Status).**
  - `crb_sep` agrees with my finite-difference Jacobian to 2.2·10⁻⁴ [own].
  - The local slopes are −2 (unknown, unequal amplitudes) and −1 (equal or known amplitudes) [own].
  - 56Δ³, 49.4Δ² and 51.7Δ² check out: 56.2, 49.6 and 52.0 [own].
- **E's new column and the offset-dependent wording (l. 50–52).**
  - The column reads 34–52%, 76–80% and 100% [md].
  - The law agrees: the tie partner sets 20–30% of the null's 5% tail at c̄ ≤ 0.02, 64–65% at 0.032 and 100% at 0.045–0.049 [own].
- **The calibrated knot is conservative (l. 54–56).**
  - F pooled 0.0414 at α = 0.05 and 0.0076 at α = 0.01 [md].
  - Nulls within 0.01 of the tie: 0.0365 and 0.0060; elsewhere 0.0442 and 0.0086 [own].
- **Numbers.** Every decimal number in (iii), the Status and the Numerics appears in the md, verbatim or as one of five roundings [own], except the misquotes in D3bn-4 and D3bn-7.

**What is wrong (majors first)**

1. **D3bn-1 (major): the two-sided half of the round's like-for-like comparison is not like for like.**
   - (iii) l. 59–60 compares "the full knot's 1.06 against the sign-agnostic |z|'s 1.09", and New l. 160 claims "compared like for like".
   - But l. 56–57 itself makes the full knot one-sided.
   - Like for like, the knot is behind:
     - a score test that spends α as the knot does reaches ½ at 1.019 × NP, and the knot (1.064) beats it in only 2 of 24 cells [own];
     - the knot's two-sided version (both walls at α/2) trails |z| in 23 of 24 cells (1.113 against 1.085) [md/own];
     - the one-sided knot trails the one-sided z in 24 of 24 cells, and in 20 of 24 with the sizes matched [own].
   - The full knot's lead over |z| is its α allocation. It pays for it with 2.09 × NP on opposite-sign pairs [md].
2. **D3bn-2 (minor): the decomposition's third part changes sign at the tie.**
   - At c̄ ≥ 0.045 the full knot's power is at or above the runner-up's own knot's in 5 of 6 cells, like for like [own].
   - The draft calls this part a loss and quotes it only at c̄ ≤ 0.032.
3. **D3bn-3 (minor): the Status equates Liu–Nehorai's GLRT limit with Smith's criterion where the two part.**
   - With unknown, unequal amplitudes, Smith's criterion sits at s ∝ Δ³, while the GLRT's local form detects at Δ² ∝ s. For (0.8, −0.4), Smith's Δ is 2.6–3.8 times NP's Δ₅₀ [own].
   - Two constants are local fits quoted as pure powers.
4. **D3bn-4 (minor): E′'s margin-over-noise column is censored at the tie.** It prints 0.89 sd where the margin is 0.72–0.73 sd [own]; the observed swap rate fits 0.73.
5. **D3bn-5 (minor): the tests still miss the data paths behind this round's headlines.**
   - 8 of D3bk's non-equivalent mutants still survive. On a rerun, D01 would print "40 of 430" for "390 of 430", and F04 would scale every Δ₅₀/NP in F by 1.41 [own].
   - Two new survivors are in U: U09, which moves the tie cell's own-knot size from 0.041 to 0.140, and U01.
6. **Cosmetic.**
   - D3bn-6: leading orders are quoted without their ranges.
   - D3bn-7: misquotes and loose ends.

## What I checked and how

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref67/D3bn_work/scripts/`; `common.py` loads the export and the raw. The mutation files are in `…/D3bn_work/mut/`.

- **Reproducibility**
  - `cmp` of the md that `--from-raw` regenerates, and `pytest -q tests/test_close_features.py`.
  - `cmp_old.py`: the raw against 8c826c6's, part by part.
  - `rerun_U.py`: all 15 U jobs.
  - `rerun_S_all.py`: all 95 S jobs.
  - `rerun_SFA.py`: 3 S jobs, 5 F jobs and 3 A rows.
  - `u_indep.py`: U written from scratch.
- **Mathematics**
  - `lars_iss.py`, `lars_dbg.py` and `lasso_kkt.py`: the ISS against LARS and the lasso.
  - `cos_theta.py`: (iv), on the grid and on a 6001-point continuum grid.
  - `crb_sep.py`: the separation's bound, plus an exact-root comparison of Smith's Δ with NP's Δ₅₀.
- **From the raw**
  - `u_decomp.py`: the decomposition, with each knot against its own law at the truth.
  - `margin_mech.py`: the frozen-margin control.
  - `margin_check.py`: E′'s margin column.
  - `tail_kernels.py` and `alloc_law.py`: who sets the null's tail.
  - `f_compare.py`, `size_matched.py`, `z_alloc.py` and `knot_vs_absz_boot.py`: like-for-like Δ₅₀ comparisons with a paired bootstrap.
  - `walls_components.py`, `cal_size_offset.py`, `z_null_shift.py` and `numbers_check.py`.
- **Mutation**
  - `mut/mutants.py` holds 71 single-point mutants of the new and changed code: `knot_of_kernel`, `u_job`, `part_U`, `own_knot`, E′'s U plumbing and margin column, `report`/`main` plumbing, `rejections`' own/other, `f_alloc_lines`, `f_lines`' new lines, `crb_sep`/`smith_line`, E's earliest-5% column and `s_job`.
  - `mut/run_mut.py` runs each mutant three ways:
    - the test file with `-x`;
    - `--from-raw` into a scratch `CF_OUT`, with a byte diff of the md;
    - `mut/behav.py`, a digest of all 15 U jobs, `part_U`'s job list and two reduced S jobs.
  - `mut/mutants_d3bk.py` ports 24 of D3bk's 25 non-equivalent survivors verbatim. J01's string now also matches `u_job`; my S01 is its port.

## Verdicts per part

| claim (D52.md line) | verdict | evidence |
|---|---|---|
| Title: the NP bound with the sign known; 1.04–1.18 × without it (l. 1–4) | **holds** | z one-sided 1.003 [0.952, 1.048] and 1.011 [0.965, 1.058]; \|z\| 1.039–1.179 [md]. The range is the spread around (z₀.₉₇₅/z₀.₉₅)^{1/2} = 1.092 (D3bn-6) |
| (i) second entry, Δ⁻², the T₃ remainder, the opposite wall, the dipole (l. 8–21) | **holds** | Unchanged; 3 A rows rerun bit-identical [own]; l. 168–174 match the md |
| (ii) stiffnesses, flow rates, two regimes (l. 23–41) | **holds** | Unchanged; raw B identical to round 64's |
| (iii) first two knots LARS = lasso; later ones while the LS fit keeps the signs and the lasso drops none (l. 46–49, 112–118) | **holds** | 400 targets, 0 counterexamples [own] |
| (iii) null law, ×10 within a cell, lower quantiles by offset (l. 49–54) | **holds** | E [md]; law-tail shares [own] |
| (iii) "sd 0.14–0.17s" (l. 54) | **low end censored** | 0.148–0.170 uncensored [own] (D3bn-4) |
| (iii) calibrated knot keeps its size, conservatively, more so at α = 0.01 and near a tie (l. 54–56) | **holds, except next to a kernel** | Near the tie 0.0365/0.0060, elsewhere 0.0442/0.0086; at c̄ = 0: 0.0512 and 0.0117 [own] |
| (iii) de facto one-sided, 75–100% (l. 56–57) | **holds** | [md]; the law gives 80–100% [own] |
| (iii) one-sided against one-sided, 1.06 vs 1.00 (l. 59–60) | **holds** | Behind in 24/24 cells; size-matched 1.045 vs 1.003 [own] |
| (iii) full knot 1.06 vs the sign-agnostic \|z\| 1.09 (l. 60) | **not like for like** | D3bn-1 |
| (iii) the loss's three parts (l. 60–64) | **holds off the tie; third part fails at it** | D3bn-2 |
| (iii) opposite-sign bias; walls at α/2 stay biased (l. 64–68) | **holds** | [md]; the walls' components [own] |
| (iii) moving null, ĉ − c̄ ∝ Δ³ (l. 68–70) | **holds to leading order** | 26% and 39% too high where the null crosses the tie (D3bn-6) |
| (iv) z's null law, mean with V's sign, ŵz sign(A₁A₂), cos (l. 72–76, 119–125) | **holds** | 8 pairs × 4 Δ [own]; O((Δ/σ)⁴) term omitted (D3bn-6) |
| (iv)(a)(b) (l. 78–84, 125–129) | **holds** | Re-derived |
| (v) (l. 86–91) | **holds** | Unchanged; D [md] |
| Status: Liu–Nehorai = Smith; constants (l. 140–145) | **exponents right; identification overreaches** | D3bn-3 |
| New (l. 158–164) | **partly** | "Like for like": D3bn-1. "Decomposed": D3bn-2. "Close to the one-sided score": 4–6% in Δ₅₀, about 20% in samples (D3bn-1). The rest plausible as new [memory] |
| Numerics (l. 166–236) | **matches md/raw** | Except 0.89 and 0.137 (D3bn-4), and "24 cells" that are 23 (D3bn-7) |

## Status of the round-64 issues

| D3bk | status | evidence |
|---|---|---|
| 1 (major) | **partly fixed** | Fixed: the title (l. 2–3), the allocation (l. 56–57, F prints 390/430 and 328/364), the one-sided per-wall tests (`rejections` l. 724–726) and the one-sided comparison. Not fixed: the two-sided half still compares the de facto one-sided full knot with \|z\| (D3bn-1) |
| 2 (major) | **fixed** | l. 66–68 "stays biased". F prints least power and pooled CP intervals. The mechanism is quoted with proxies (D3bn-7) |
| 3 | **fixed in substance** | Right exponents and mechanism; `crb_sep`/`smith_line`. Two constants are local fits; the identification overreaches (D3bn-3) |
| 4 | **fixed** | Proof l. 112–118 correct; 400 targets [own] |
| 5 | **fixed, with one latent copy** | l. 74–75 and 122; brackets l. 5–6; D‴.54's time t. A test with (−0.8, −0.6) kills the sign bug in `rejections`. The same bug survives in `f_alloc_lines` l. 878 (D3bn-5) |
| 6 | **fixed** | l. 68–70; F prints ĉ − c̄, q·s and the median ×s; offsets at l. 232–233 |
| 7 | **fixed** | l. 54–56; F prints the α = 0.01 sizes |
| 8 | **fixed, with a new omission** | Wording at l. 50–52, E's column and Part U are in place. The third part's sign at the tie is omitted (D3bn-2); the margin column is censored (D3bn-4) |
| 9 | **partly fixed** | Now killed: W01, W02, Z04, J01–J07, P01, P02, F01, F02 and X02 (D3bf's M08). Still non-equivalent survivors: K01, D01, D02, Z02, S03, F03, F04, F07 (D3bn-5). J08 is equivalent (R̂φ_{j₂} ⟂ f₀: ratio 4·10⁻¹²) and S05 is equivalent on this data (no width escapes in 113 800 rows) [own] |
| 10 | **fixed** | l. 76–77 and 82–84 (Fukumizu–Amari, Liu–Wu–Wang, Shahram–Milanfar; (b) needs a merged block); New l. 162–163 |
| 11 | **mostly fixed** | "Four times" (l. 224); the margin column (now censored, D3bn-4); A bit-identical. z's tested range is still not stated |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bn-1 | major | **The two-sided half of the round's like-for-like comparison is not like for like** (l. 59–60; New l. 160).<br>**What the draft does.** It calls the full knot one-sided (l. 56–57: 90.7% and 90.1% of its null rejections at the fit's-sign wall [md]; 80–100% of its law's tail [own]). It then sets it against the two-sided \|z\|. That is D3bk-1's comparison again.<br>**Same allocation.** A score test that spends α as the knot does rejects when ŵz > 1.692 or −ŵz > 2.601; its H0 size is 0.0457. It reaches ½ at 1.019 [0.970, 1.064] × NP. The calibrated knot (1.064) beats it in 2 of 24 cells; at α = 0.01, 1 of 24 (1.008 against 1.037). For (0.8, −0.4) the same test reaches 1.261 against the knot's 2.093, with the knot ahead in 0 of 23 cells [own, z_alloc.py].<br>**Two-sided.** Both walls at α/2 trail \|z\| in 23 of 24 cells: 1.113 against 1.085 [md], paired median +2.6% [+1.9, +3.6]. At α = 0.01 they trail in 22 of 22 cells that cross (1.108 against 1.057) [own].<br>**One-sided.** The fit's-sign wall trails z in 24 of 24 cells (1.056 against 1.003). With the sizes matched (α′ = 0.0546) it still trails in 20 of 24 (1.045 against 1.003) [own]. Since Δ₅₀ ∝ N^{−1/4}, that is 18–24% more samples.<br>**The apparent lead.** The full knot leads \|z\| in 17 of 24 cells (paired median −1.3% [−2.8, −0.9]). That lead is its allocation, which costs 2.09 × NP against \|z\|'s 1.09 for (0.8, −0.4) [md]. | **(iii) l. 59–60:** "two-sided against two-sided, both walls 1.11 against \|z\|'s 1.09 (behind in 23 of 24 cells); with the knot's own allocation the score reaches ½ at 1.02 against the knot's 1.06".<br>**New l. 160–161:** quantify "close": 4–6% in Δ₅₀, about 20% in samples.<br>**F:** print the allocated score and the paired differences. |
| D3bn-2 | minor | **The decomposition's third part changes sign at the tie** (l. 60–64; E′ summary at exp l. 828–837; Numerics l. 210–215).<br>**At the tie.** With each knot against its own law at the truth (E's `q_law`, U's `q2`), the full knot minus the own knot is +0.005, +0.005, +0.0125, +0.0125 and +0.0175 in 5 of the 6 cells at c̄ = 0.045–0.049 (−0.0075 in the sixth). In those five cells the full knot rejects 2–7 draws the own knot does not, and never the reverse [own]. Against E′'s "larger of two" the gain reaches +0.068 (0.233 against 0.165 at (0.05, 0.03, 0.049)) [md].<br>**Why.** At c̄ ≥ 0.045 the tie partner sets 100% of the null's tail [own], so far kernels do not share it. Under H1 they still reach their wall when the partner's noisy margin fails.<br>**Off the tie.** Like for like, the minimum costs 0.030–0.108 [own]; the md's mixed-calibration figure is 0.012–0.098.<br>**The margin part acts without swaps.** There are 0 swaps at (0.0175, 0.003, 0.049) and at (0.05, 0.03, 0.045), yet the losses to z₂ are 0.040 and 0.089 [md]. Freezing the margin removes them [own]. | **(iii):** "the minimum: a loss where far kernels share the null tail (0.03–0.11 at c̄ ≤ 0.032), a small gain next to a tie; the margin's noise enters the knot's numerator μ̂ (and near a tie swaps the first entrant)".<br>**E′:** print the full knot against its law at the truth, and the third part in all 15 cells. |
| D3bn-3 | minor | **The Status identifies Liu–Nehorai's GLRT limit with Smith's criterion where they part** (l. 141–145).<br>**The contradiction.** "Whose GLRT resolution limit ∝ √CRB is Smith's 2005 Cramér–Rao criterion" and "unknown and unequal … s ∝ Δ³" together put a GLRT's limit at s ∝ Δ³. Yet (iii)–(iv) show that the score, the GLRT's local form, detects such pairs at Δ² ∝ s.<br>**Numbers.** Smith's Δ (exact root of sd(ĉ₂ − ĉ₁) = Δ) over NP's Δ₅₀ at s = 0.003, 0.01, 0.03 is 1.21, 1.00, 0.85 for (0.8, 0.6) and 3.77, 3.08, 2.57 for (0.8, −0.4) [own].<br>**Constants.** "474Δ³" and "109Δ²" are fits that `smith_line` prints as 474Δ^2.99 and 109Δ^1.98. As pure powers the constants are 500 (506 as Δ → 0) and 119 (121) [own]. | **Status:** restrict the identification to equal or known amplitudes. Say that, otherwise, Smith's criterion is an estimation threshold, not a detection limit. Quote 500Δ³ and 119Δ², or the exponents. |
| D3bn-4 | minor | **E′'s margin column is censored at the tie** (exp l. 812, `s_job` l. 513; Numerics l. 213–214; (iii) l. 54).<br>**The cause.** μ̂ = 1 − wall·u_{j₂} is 0 whenever j₂ has the largest correlation.<br>**The effect.** At (0.05, 0.03, 0.049), where 23% of draws swap, sd(μ̂)/s is 0.137, against 0.170 for the uncensored 1 − wall·g_{j₂}/\|g_{j₁}\|. The uncensored value lies in 0.148–0.170 in every cell [own].<br>**Consequence.** The margin is 0.72–0.73 sd, not 0.89. Φ(−0.73) = 0.232 matches the observed 0.233 swaps; Φ(−0.89) = 0.187 does not. | **E′:** use the uncensored sd and the null's margin.<br>**Text:** l. 54 "0.15–0.17s", l. 213 "0.73". |
| D3bn-5 | minor | **The tests miss the data paths behind the round's headlines.** `--from-raw` cannot see a data-path mutant.<br>**My 71 mutants.** 10 are killed by both checks, 28 by the tests only and 25 by `--from-raw` only; 8 survive both.<br>– **U09** (`u_job` l. 928, the dual divided by \|g_{j₁}\| instead of the max): at the tie cell, size 0.041 → 0.140 and power 0.165 → 0.310 [own]. The swap branch is never tested; the test runs at c̄ = 0.032.<br>– **U01** (l. 913, the runner-up forced to j₁ + 1): at c̄ = 0, power 0.637 → 0.593, 0.378 → 0.365, 0.545 → 0.557 [own].<br>– **U10**: U's `first_0/1` fields are read by nothing.<br>– **A06** (l. 878, the sign of ŵ dropped): latent, since all 97 800 fits have ŵ > 0 [own].<br>– **O02, A01, J01, S03**: equivalent on this data.<br>**D3bk's survivors still alive, with their effect on a rerun [own].**<br>– D01 (`draw_stats`' walls swapped): "40 of 430" for "390 of 430", and an 'own' size of 0.006.<br>– K01 (`knot_walls` swapped): 'own' size 0.14–0.48.<br>– F03 (qA at α): the walls test's size 0.0463 → 0.0702.<br>– F04 (NP at s/2): every Δ₅₀/NP × 1.41.<br>– F07: F's nulls move to the second separation.<br>– Z02: E′'s t size → ≈ 0.<br>– D02: E′'s one-sided z_σ power → ≈ 0.<br>– S03 (5% escape window): about 20% of the (0.8, −0.4), s = 0.03 nulls escape.<br>**Also:** `c_prime_lines` is still untested, and `cp_int` is still tested against itself. | **Tests:**<br>– `knot_walls` and `draw_stats`' rows against an independent wall definition; e.g., run `test_the_walls_take_the_pairs_sign` through `f_job` rows.<br>– `f_job`'s qA and np, and `f_jobs`' H0 separation.<br>– `u_job` at a swapping tie cell and at c̄ = 0.<br>– t against its formula; z_σ's sign; the escape window on a (0.8, −0.4) null.<br>– the α = 0.01 sizes on an A < 0 null. |
| D3bn-6 | cosmetic | **Leading orders quoted without their range.**<br>**cos (l. 74, 124).** The formula omits an O((Δ/σ)⁴) term. With ϑ = 0, 1 − cos = 7.7·10⁻⁸ … 3.1·10⁻⁴ at Δ = 0.0125 … 0.1. For (0.8, 0.6) the ratio to ϑ²/12 is 1.007, 1.029, 1.119, 1.498 [own].<br>**ĉ − c̄ ∝ Δ³ (l. 68–69).** It is 26% and 39% too high at Δ = 0.05 and 0.07, where the null crosses the tie (0.0097 against 0.0130; 0.0217 against 0.0357) [own].<br>**Title "1.04–1.18".** This is the spread over the cells around 1.092; the medians are 1.085 and 1.094. | Write "+ O(ϑ⁴ + (Δ/σ)⁴)", "to leading order" and "≈ 1.09 (1.04–1.18)". |
| D3bn-7 | cosmetic | **Smaller misquotes and loose ends.**<br>– l. 218 "24 cells": the opposite-sign knot's median and range use 23 (cell (0.003, 0.045) never crosses).<br>– l. 67 and l. 227–228: the walls mechanism quotes the per-draw one-sided tests at α, not the walls test's own α/2 components (0.0054/0.0015 against 0.0223; 0.0310/0.0379 against 0.0235) [own].<br>– l. 55 "conservatively": not at c̄ = 0 (0.0512 [0.0428, 0.0608]; 0.0117 at α = 0.01) [own].<br>– z's tested range (s/\|ŵ\| ≤ 0.079) is not stated.<br>– E′ omits z₂'s size.<br>– `smith_line` prints inside D (exp l. 1226), not B.<br>– The log's "but for the rms": the end coefficients also differ, by ≤ 6.5·10⁻¹².<br>– Neyman–Scott (l. 146) is usually dated 1966 [memory]. | Correct the text and the log. |

## Sources

- **Read in the export:**
  - in full: D52.md, `exp_close_features.py`, close_features.md, `test_close_features.py`, REVIEW_round64_D3bk.md and the Round 64 log entry;
  - in part: `rich_transit.py` (l. 1–80, 200–510), `mirror_free.py` `iss_path` (l. 635–692), D54.md (l. 1–12, 119–135, 190–200) and D59.md (l. 1–10).
- **Literature.** My one web search was refused (the session's shared search budget was exhausted), and arxiv.org and ams.org are unreachable. Everything below rests on memory or on D3bk's recorded snippets:
  - Neyman & Scott, "On the use of C(α) optimal tests of composite hypotheses", Bull. ISI 41(1):477–497. D3bk's snippet dates it 1965; most citations say 1966 [memory].
  - Fukumizu & Amari, Neural Networks 13 (2000) 317–327; Liu, Wu & Wang, "Splitting steepest descent", NeurIPS 2019 [memory].
  - Shahram & Milanfar: IEEE TIP 2004, TSP 2005, TIT 2006 [memory; D3bk snippets].
  - Liu & Nehorai, IEEE TSP 55(11), 2007: "proportional to the square root of the CRB" [D3bk snippet]. I could not check its amplitude assumptions, so D3bn-3 rests on the internal argument.
  - Smith, IEEE TSP 53(5), 2005 [memory].
  - Tibshirani–Taylor–Lockhart–Tibshirani, JASA 2016; Azaïs–De Castro–Mourareau, Bernoulli 2018, and the grid-less version (ACHA; arXiv 2017); Lockhart et al., Ann. Statist. 2014; McLachlan, JRSS C 1987 [memory].
- **Novelty.** On these sources, "New" stands for the transit's delay law and for the knot's offset dependence, allocation and opposite-sign bias. A fixed-grid spacing test is no substitute here: an off-grid null bump is not a single atom of the grid.

## Lesson

A decomposition and a comparison are claims about every cell, signs included. Part U confirmed the knot's three loss mechanisms where the draft looked, and a frozen-margin control confirms the second one exactly. At the tie, though, the third part flips: the dictionary's other kernels recover power that the runner-up's noisy margin loses. The draft quotes that part only where it is a loss. Round 64's rule was one-sided against one-sided and two-sided against two-sided, but it was applied to only one half of the comparison. The full knot, which the draft itself calls one-sided, is set against the two-sided |z|. Meanwhile the knot's genuinely two-sided version, already in the md, trails |z| in 23 of 24 cells, and a score spending α the knot's way beats it in 22. Before calling a comparison like for like, check that both members spend α the same way. And re-check every fix against the mutants that motivated it: half of D3bk's surviving mutants still survive.
