# Referee report D3br, round 71: Proposition D‴.52 as revised after round 67

## Scope

- **Commit reviewed:** 42aa6a3 on branch claude/festive-bardeen-eopej2 ("D‴.52 round 67: like for like on both halves; the minimum a gain at the tie").
  - I exported it with `git archive 42aa6a3` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/repo` and worked only there.
  - For diffs, I exported round 67's files (766c864) to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/old/r67`.
- **Files reviewed:**
  - `review/drafts/D52.md` (250 lines)
  - `experiments/exp_close_features.py` (1387 lines)
  - `results/chaos/close_features.md` (301 lines) and `close_features_raw.npz`
  - `tests/test_close_features.py` (37 tests)
  - `review/REVIEW_round67_D3bn.md` and the Round 67 entry of `review/REVIEW_LOG.md`
  - For context: `src/core/rich_transit.py` (`sliding_iss` l. 304, `polish_block` l. 202) and `src/core/mirror_free.py` (`iss_path` l. 635).
- **HEAD moved** to 1413849 while I worked; the new commits are the D3bq and D3bp reports. None of the reviewed files differ between 42aa6a3 and HEAD. I reviewed 42aa6a3 only.
- **Container restart.** The container restarted after all my runs had finished. Afterwards I re-checked two things:
  - the export's raw has the same md5 as 42aa6a3's;
  - `--from-raw` still regenerates the md byte for byte.
  Nothing was lost or rerun.
- **How I ran things:**
  - the repo's venv with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`;
  - `CF_OUT` always pointed to a scratch copy of the raw;
  - at most two heavy processes at a time; nothing of mine is left running.
- **Tags:** every number is marked [own] (my run), [md] or [raw].
- **I modified nothing in the repository.**

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (`cmp`) [own].
  - The 37 tests pass in 31.8 s [own].
  - Against round 67's raw (766c864), A, B, C, E, F and U are identical (recursive exact comparison). S is identical except for the 21 000 added `mu_u` fields (15 E cells × 1400 draws) [own].
  - Rerun with the experiment's own functions and arguments, these are bit-identical to the raw, `mu_u` included [own]:
    - all 15 U jobs;
    - all 95 S jobs;
    - two F jobs: (0.8, −0.4) at (c̄, s, Δ) = (0.032, 0.01, 0.07), and the (0.8, 0.6) null at (0.049, 0.01).
- **Numbers.** D52.md has 523 decimal numbers. 517 appear verbatim in the md, 5 are roundings of md numbers (0.563, 1.053, 32.3, 0.0254, 0.88), and "90.7%" is the md's 0.907 [own].
- **The like-for-like comparison (D3bn-1) reproduces and is robust.**
  - An independent re-implementation of the rejections, `crossing` and pairing gives every F figure [own]:
    - allocation 0.907 (390/430) and 0.901 (328/364);
    - allocated score Δ₅₀/NP 1.019 [0.970, 1.064] and 1.257 [1.189, 1.389];
    - paired 1.045 (2 of 24), 1.026 (1 of 24), 1.043 (0 of 24), 0.987 (17 of 24).
  - The conclusion survives every other way of matching the knot's allocation for (0.8, 0.6) [own]:

    | allocation matched | paired median | knot ahead |
    |---|---|---|
    | per offset | 1.041 | 0 of 24 |
    | per cell | 1.040 | 0 of 24 |
    | law-calibrated knot's allocation (0.915 = 422/461) | 1.047 | 2 of 24 |
    | realized per-wall H0 rates matched exactly (knot 0.0406/0.0042; score at ŵz > 1.699, −ŵz > 2.589) | 1.042 | 2 of 24 |

  - F's own allocated score already realizes 0.0418/0.0040 under H0, close to the knot's rates [own].
  - A paired bootstrap (H1 draws resampled per cell) gives a median in [1.034, 1.057], with the knot ahead in 1–5 cells [own].
  - The two-sided half is size-matched: both walls 0.0475 against \|z\| 0.0460 [own].
- **F's new lines recompute exactly** [own]:
  - per-offset sizes 0.0512 [0.0428, 0.0608] … 0.0329 [0.0261, 0.0409] (2400 draws each), and 0.0054–0.0117 at α = 0.01;
  - the walls test's components: H0 0.0223/0.0235 and H1 0.0034/0.0345 for (0.8, −0.4);
  - s/\|ŵ\| from 0.00214 to 0.07346.
- **The uncensored margin is the right quantity.**
  - P(μ̂_u < 0) equals the observed swap rate in all 15 cells (0.233/0.233, 0.016/0.016, 0/0 elsewhere) [own].
  - It moves 2.12 → 2.06 and 0.89 → 0.71 [md]. No conclusion changes.
- **The minimum's part off the tie is a real loss.** It is −0.112 … −0.028 with both laws on common noise (means over 10 seeds) [own], against −0.1075 … −0.0300 in the md [md].
- **(iv) checks out.**
  - At ϑ = 0, 1 − cos = 6.510·10⁻⁴(Δ/σ)⁴ for Δ = 0.1 … 0.0125 (continuum, 40 001 points) [own].
  - For (0.8, 0.6) at Δ = 0.0125, (1 − cos)/(ϑ²/12) = 1.008 [own].
  - I re-derived σ̂ − σ = sz/(ŵσ‖R̂∂²_cφ‖), the continuum constants √π/(2σ³), 3√π/(4σ⁵) and 3/(2σ²), and ĉ − c̄ = −T₃/(4Aσ²) [own].
- **Smith (Status).**
  - `crb_sep` agrees with my finite-difference Jacobian to ≤ 5.2·10⁻⁶ [own].
  - The integer-power constants check out [own]:
    - (0.8, 0.6) unknown amplitudes: 500.5, 504.9, 506 Δ³ at Δ = 0.0125, 0.00625, 0.003125;
    - (0.8, 0.6) known: 52.0Δ²;
    - (0.8, −0.4): 56.2Δ³ (unknown) and 119.2 → 121.4 Δ² (known);
    - (0.8, 0.8) unknown: 49.6Δ².

**What is wrong (majors first)**

1. **D3br-1 (major). The "small gain at the tie" is a Monte Carlo artifact** ((iii) l. 67–69, New l. 168–169, Numerics l. 248–249, E′, README).
   - The full knot and the runner-up's own knot are set against two independently sampled 20 000-sample laws (seeds [seed, 9] and [seed, 7]).
   - At c̄ ≥ 0.045 every decision difference lies between those two thresholds. Their ratios are 0.997–1.032, and the sign of that ratio is the sign of the "part".
   - Drawn from the same samples, the two thresholds coincide exactly, and the part is +0.0000 in 6 of 6 cells (and in 60 of 60 cell-seeds) [own].
   - At the tie the minimum costs nothing and gains nothing.
2. **D3br-2 (minor). The tests the log credits for Z02, D02 and the uncensored margin never reach the mutated code.**
   - Z02 and D02 survive both checks and would change E′ on a rerun:
     - Z02: the t size becomes 0.000 in all 15 cells;
     - D02: the z_σ power becomes 0.000–0.007 [own].
   - `s_job`'s `mu_u` is untested. Reverting it to the censored margin (SM01) survives both checks and would reprint round 67's 0.137 / 0.89 [own].
3. **D3br-3 (minor). The one-sided half compares a conservative knot at the nominal level.**
   - The sizes are 0.0446 (the knot's fit's-sign wall) against 0.0478 (one-sided z) [own].
   - At equal size the knot is behind in 20 of 24 cells, by ≈3% in Δ₅₀ (≈13% in samples) [own].
   - So "behind in 24 of 24" (l. 61) and New's "4–6% … about 20% in samples" (l. 166–167) overstate.
4. **D3br-4 (cosmetic). The Status's Liu–Nehorai sentence attributes the parting to the wrong thing.** A limit ∝ √CRB coincides with Smith's criterion for any amplitudes. What parts is that form from the GLRT's actual detection order.
5. **D3br-5 (cosmetic).** Notation, scope and loose ends (see Issues).

## What I checked and how

**Reproducibility**
- `cmp` of the md that `--from-raw` regenerates, before and after the container restart.
- `pytest -q -p no:cacheprovider tests/test_close_features.py`.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/cmp_all67.py`: the raw against 766c864's, recursively and exactly.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/rerun_SU.py`: all 95 S jobs and all 15 U jobs, with `part_S`/`part_U`'s seeds and fits.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/rerun_F.py`: two F jobs.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/numbers_check.py`: every decimal number in D52.md against the md.

**The allocation and like for like**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/alloc_check.py`: an independent re-implementation, plus allocations per pair, per offset, per cell, from the law-calibrated knot, and from absolute sizes; also size-matched variants.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/alloc_realized.py`: the score's thresholds set to the knot's realized per-wall H0 rates.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/alloc_boot.py`: paired bootstrap of the allocated comparison.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/onesided_boot.py`: paired bootstrap of the one-sided comparison.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/sizes_pair.py`: realized sizes per pair.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/f_newlines_check.py`: F's new lines, the draws' s/\|ŵ\|, the one-sided comparison at matched size, and precise medians.

**E′**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/tie_draws.py`: the full knot (E's Δτ₂) against U's own knot draw by draw, and q_law against q2.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/min_part_seed9.py`: U's law drawn from E's own law samples.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/min_part_crn.py`: 10 further common-noise law seeds, and the null's margin.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/mu_cmp.py`: the censored against the uncensored margin, and swaps.

**Mathematics**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/cos_check.py`: (iv)'s cos in the continuum.
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/smith_check.py`: `crb_sep` against finite differences, the constants, and exact roots of Smith's criterion against NP's Δ₅₀.
- By hand: the knot formula; n's covariance s²ΦᵀRΦ/N²; (ii)'s continuum Schur constants; the odd/even ratio ϑ/√6; (iv)(a); and the allocated score's population bias for opposite signs, which follows from −φ(1.692) + φ(2.601) < 0.

**Mutation**
- Mutant definitions: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/mutants.py`. It holds 70 single-point mutants:
  - round 67's 10 survivors ported verbatim, plus U10 and A06;
  - 58 mutants of the new code: `knot_allocation`, `paired`, `rejections`' alloc, `f_lines`' paired/range lines, `f_alloc_lines`' per-offset/components lines, `e_prime_lines`' new columns and the minimum's part, `s_job`'s `mu_u`, and `smith_line`/`crb_sep`.
- Runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/run_mut.py` runs each mutant three ways:
  - the full test file;
  - `--from-raw` with a byte diff of the md;
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/behav.py`, a digest of all 15 U jobs, `part_U`/`part_S`'s job lists, reduced S jobs (tie and centre E cells, a C null), reduced F jobs (both pairs, H1 and H0), `draw_stats` and `smith_line`.
- Effects of the survivors on a rerun: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/z02_d02_effect.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/scripts/sm_effect.py`.
- Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/results_a.jsonl`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/results_b.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref71/D3br_work/mut/results_c.jsonl`.
- Outcome [own]:
  - **Round 67's set (12):**
    - killed by the tests: D01, K01, F03, F04, F07, S03, U09, U01, U10;
    - surviving both checks and changing a rerun: Z02, D02;
    - latent (equivalent on this data): A06.
  - **New code (58):** 19 killed by the tests; 28 killed only by the byte diff; 11 survive both:
    - SM01, SM05: change the md on a rerun;
    - SM03: changes the raw only;
    - SM02, SM04, PA03, FL04, FL05, EP04, EP10: equivalent on this data;
    - FL01: equivalent, because the ratio cancels.
  - KA06 (p ≤ α) changes 18 md lines because 14 of the 19 200 H0 p-values equal α exactly. The code is consistent: p < α in both `knot_allocation` and `rejections`.

## Verdicts per part

| claim (D52.md line) | verdict | evidence |
|---|---|---|
| Title: the NP bound with the pair's sign known; ≈ 1.09 (1.04–1.18) without it (l. 1–4) | **holds** | z one-sided 1.003 and 1.011 [md]. \|z\| medians 1.0851 and 1.0936, range 1.039–1.179 [own/md] |
| (i) second entry, Δ⁻², opposite wall, dipole (l. 8–21) | **holds** (unchanged) | A identical to 766c864's raw [own] |
| (ii) stiffnesses, rates, two regimes (l. 23–41) | **holds** (unchanged) | B identical [own]; continuum constants re-derived [own] |
| (iii) the knot's law, LARS = lasso for k ≤ 2, the null across the offset (l. 43–54) | **holds** (unchanged) | C, E identical [own] |
| (iii) "sd 0.15–0.17s" (l. 54) | **holds** | Uncensored 0.148–0.170 [md]. P(μ̂_u < 0) equals the swap rate in all 15 cells [own] |
| (iii) the calibrated knot's size; c̄ = 0: 0.0512 [0.0428, 0.0608], 0.0117 (l. 54–56) | **holds; wording** | Per offset recomputed [own]. "Conservatively off a kernel's centre" fails at c̄ = 0.02: 0.0500 [0.0416, 0.0595] [md] (D3br-5) |
| (iii) de facto one-sided, 75–100% (l. 56–58) | **holds** | Per offset 0.754–1.000 [own] |
| (iii) one-sided: 1.06 against 1.00, "behind it in 24 of 24" (l. 60–61) | **holds at the nominal level only** | Sizes 0.0446 against 0.0478. At equal size: paired 1.032, knot ahead in 4 of 24; bootstrap [1.020, 1.041] [own] (D3br-3) |
| (iii) two-sided: both walls 1.11 against \|z\| 1.09, behind in 23 of 24 (l. 61) | **holds** | Sizes 0.0475 and 0.0460; paired 1.026 [own] |
| (iii) allocated score: 1.06 against 1.02, ahead in 2 of 24 (l. 61–63) | **holds, robust** | Per offset: 0 of 24. Realized rates matched: 1.042, 2 of 24. Bootstrap [1.034, 1.057], 1–5 cells [own] |
| (iii) the lead over \|z\| is the allocation; 2.09 against 1.26 (l. 63–64) | **holds** | 0.987 (17 of 24) [md]. With realized rates matched, the allocated score for (0.8, −0.4) is 1.296 [own] |
| (iii) loss parts 1–2: the filter, the margin's noise (l. 64–67) | **holds** | z₂ against the own knot [md]; the margin is right [own] |
| (iii) part 3 off the tie: a loss (l. 67–68) | **holds** | Common-noise laws: −0.112 … −0.028 [own] |
| (iii) part 3 at the tie: "a small gain … other kernels still reach their wall" (l. 68–69) | **fails** | Part = 0.0000 in 6/6 cells with shared law samples; no rejecting draw has full < own [own] (D3br-1) |
| (iii) opposite-sign bias; walls at α/2 biased (l. 69–73) | **holds** | Least power and components recomputed [own] |
| (iii) moving null "to leading order" (l. 73–75) | **holds** | −T₃/(4Aσ²) re-derived [own] |
| (iv) z's law, mean, cos = 1 − ϑ²/12 + O(ϑ⁴ + (Δ/σ)⁴) (l. 77–82; proof l. 124–131) | **holds** | Continuum check [own]. The proof's intermediate step hides ϑ²(Δ/σ)² (D3br-5) |
| (iv)(a)(b) (l. 83–89, 131–134) | **holds** | Re-derived [own] |
| (v) (l. 91–96) | **holds** (unchanged) | D identical [own] |
| Status: LN = Smith only with equal or known amplitudes; an estimation threshold (l. 145–151) | **holds in substance; wording** | Constants checked [own]; D3br-4 |
| Status: losses "decomposed … two pairs, eight offsets (F)" (l. 140–143) | **overstated scope** | The decomposition is E′'s: (0.8, 0.6), c̄ ∈ [0, 0.049] (D3br-5) |
| New (l. 164–171) | **partly** | "A small gain at it" is false (D3br-1). "4–6% … about 20% in samples" is at the nominal level (D3br-3). The rest is plausible as new [memory] |
| Numerics (l. 173–250) | **match the md** | Except "positive in 5 of 6" (l. 249, an artifact: D3br-1), and two figures for the same part (l. 218–219 against 248–249) |

## Status of the round-67 issues

| D3bn | status | evidence |
|---|---|---|
| 1 (major) | **fixed** | (iii) l. 60–64 and New l. 166–167 rewritten. F prints the allocated score and four paired comparisons. All figures reproduce independently and survive per-offset, per-cell, law-knot and realized-rate allocations and a bootstrap [own]. Residue: the one-sided half is compared at unequal sizes (D3br-3) |
| 2 (minor) | **implemented, but the effect is spurious** | E′ prints the part per cell, and (iii)/New/README say "a small gain at the tie". That gain is a law-sampling artifact; my predecessor's evidence rested on the same two independent laws (D3br-1) |
| 3 (minor) | **fixed** | Status l. 146–151 restricted. `smith_line` prints 501/505Δ³ and 119/121Δ² [md]; checked [own]. Wording residue (D3br-4) |
| 4 (minor) | **fixed in the numbers** | `mu_u` gives 0.148–0.170, 2.06, 0.71 [md]. Still untested (D3br-2). μ is still H1's margin (`mu_true`, exp l. 525), not the null's as asked: 0.73 with the null's [own] |
| 5 (minor) | **partly fixed** | Now killed: D01, K01, F03, F04, F07, S03, U09, U01, U10. Still surviving: Z02, D02. New survivors: SM01, SM05, SM03. Unchanged: A06 latent, `c_prime_lines` untested, `cp_int` tested against itself (D3br-2) |
| 6 (cosmetic) | **fixed** | l. 79 and 129: O(ϑ⁴ + (Δ/σ)⁴). l. 74: "to leading order". Title "≈ 1.09 (1.04–1.18)". Residue in the proof (D3br-5) |
| 7 (cosmetic) | **mostly fixed** | 23 cells (l. 226); walls components (l. 237–238); size per offset (l. 56, 246–247); z's range (l. 250); z₂'s size in E′ [md]; Neyman–Scott 1966. Still open: "conservatively" at c̄ = 0.02, and REVIEW_LOG l. 2576 "but for the rms" (D3br-5) |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3br-1 | major | **The minimum's "small gain at the tie" is law-sampling noise.**<br>**Where:** (iii) l. 67–69; New l. 168–169; Numerics l. 248–249; E′ summary (exp l. 855–860); README; log.<br>**How the part is computed.** The full knot is tested against E's `q_law` (`e_job`, rng [seed, 9], exp l. 329). The own knot is tested against U's `q2` (`u_job`, rng [seed, 7], exp l. 977). These are two independent 20 000-sample laws.<br>**At the six cells with c̄ ≥ 0.045:**<br>– q_law/q2 = 1.0083, 1.0153, 0.9968, 1.0269, 1.0091, 1.0321 [own, from raw];<br>– the md's parts are +0.0050, +0.0050, −0.0075, +0.0125, +0.0125, +0.0175 [md]; their signs are the signs of q_law − q2;<br>– every H1 decision difference (2, 2, 3, 5, 5, 7 draws) is a draw whose knot lies between the two thresholds;<br>– draws whose full knot is below the own knot and below either threshold: 0 of 6000 under H0, 0 of 2400 under H1 [own].<br>**With U's law drawn from E's own samples:** q2 = q_law exactly (1162.54, 268.549, 345.725, 55.4552, 102.453, 13.7554), and the part is +0.0000 in all six cells. Off the tie it is −0.1225 … −0.0250 [own].<br>**Ten further shared-noise seeds:** 0.0000 in 60 of 60 cell-seeds, with identical sizes [own].<br>**Why.** The full knot ≤ the own knot draw by draw. At the tie the pair sets the law's whole 5% tail sample by sample, so the quantiles coincide. Under H1 other kernels do reach their wall first in 7–82 draws per cell, but always above both thresholds [own]. | **u_job:** draw the own knot's law from the same noise as E's law (rng [seed, 9]), or compute both in one job. Otherwise print each law's Monte Carlo error: ±0.0015 in level ≈ ±0.006 in power.<br>**(iii) l. 67–69:** "a loss where far kernels' knots share the null's tail (c̄ ≤ 0.032; 0.03–0.12), none next to the tie, where the tie pair sets the whole tail and no other kernel's knot falls below it under H1".<br>**New l. 168–169; Numerics l. 248–249:** drop "a small gain" and "positive in 5 of 6".<br>**README and REVIEW_LOG** round 67. |
| D3br-2 | minor | **Tests named for mutants do not run the mutated code.**<br>**What the new test does.** `test_the_e_prime_report_with_the_runner_ups_knot` (tests l. 634–652) feeds `e_prime_lines` hand-made rows. `score_stats`' t, `draw_stats`' zs and `s_job`'s `mu_u` are never computed in any test. `test_s_job_on_an_e_cell_and_a_c_null` checks `mu` (l. 799), not `mu_u`.<br>**What survives both the tests and the byte diff [own]:**<br>– Z02 (t without √(N−2)): on a rerun E′'s t size is 0.000 in all 15 cells (now 0.039–0.063), and the pooled 0.0506 becomes 0.0000;<br>– D02 (z_σ's sign flipped): E′'s z_σ power is 0.000–0.007 (now 0.315–0.710);<br>– SM01 (`mu_u` censored again, exp l. 514): the tie cells reprint round 67's 0.137/0.89 and 0.164/2.12 (now 0.170/0.71, 0.169/2.06), undoing D3bn-4 silently;<br>– SM05 (noiseless numerator): 0.292/0.42;<br>– SM03 (1 + …): corrupts the raw's `mu_u` while the md is unchanged, since sd(2 − μ̂) = sd(μ̂);<br>– SM02, SM04: equivalent only because every E cell has wall +1 and g_{j₁} > 0.<br>**Covered only by the byte diff.** These new lines have no test at all: `f_alloc_lines`' per-offset and components lines (FA01–FA09), `smith_line`'s constants (SL01–SL05), and `f_lines`' like-for-like/range lines (FL02, FL03, FL06–FL09). | **Tests:**<br>– t and z_σ through `draw_stats` on fixed labels against their formulas: t = z·s/(‖r‖/√(N−2)), and sign z_σ = sign(ŵ)·sign z;<br>– `mu_u` = 1 − wall·g_{j₂}/\|g_{j₁}\| for the pair of each sign, with P(mu_u < 0) = swap at a swapping tie cell (`e_jobs()[14]`);<br>– `f_alloc_lines` and `smith_line` on synthetic F and a closed-form `crb_sep`. |
| D3br-3 | minor | **The one-sided half compares at unequal sizes** (l. 60–61; New l. 166–167).<br>**Sizes.** (0.8, 0.6) H0: the fit's-sign wall 0.0446, one-sided z 0.0478. The other halves are matched: both walls 0.0475 against \|z\| 0.0460; the calibrated knot 0.0448 against the allocated score 0.0457 [own].<br>**At equal size** (the wall at α′ = 0.055, size 0.0476): Δ₅₀/NP 1.045; paired 1.032 [0.973, 1.239], with the knot ahead in 4 of 24. Bootstrap: median [1.020, 1.041], ahead in 2–7 cells [own].<br>**At the nominal level:** [1.031, 1.054], ahead in 0–4 [own].<br>**Consequence.** "Behind it in 24 of 24" and "4–6% in Δ₅₀ (about 20% in samples)" partly measure the knot's conservativeness. At equal size the gap is ≈3% (≈13% in samples). | **l. 60–61 and New:** "behind it at equal size by ≈3% in Δ₅₀ (≈13% in samples, ahead in 4 of 24 cells); in 24 of 24 at the nominal level, where its size is 0.045 against z's 0.048".<br>**F:** print the size-matched comparison. |
| D3br-4 | cosmetic | **Liu–Nehorai and Smith** (Status l. 146–147).<br>**The logic.** "Whose GLRT resolution limit ∝ √CRB coincides with Smith's … only with equal or known amplitudes" is self-contradictory. A limit ∝ √CRB equals Smith's criterion, up to its constant, for any amplitudes. What parts with unknown unequal amplitudes is the √CRB form from the GLRT's actual detection order Δ² ∝ s.<br>**Known amplitudes.** Smith's Δ is 0.50–0.51 × NP's Δ₅₀ for both pairs and all three s [own]: the same order, but a constant ≈4 in s.<br>**Source.** Whether LN's model assumes equal or known amplitudes rests on [memory]. | **Status:** "Liu–Nehorai's √CRB limit, like Smith's criterion, has the detection order Δ² ∝ s only with equal or known amplitudes (Smith's Δ then ≈ ½ NP's Δ₅₀); with unknown unequal amplitudes both are estimation thresholds". |
| D3br-5 | cosmetic | **Loose ends.**<br>– l. 62: "ŵz > 1.692 or −ŵz > 2.601" means sign(ŵ)z. Taken literally, ŵz > 1.692 with ŵ ≈ 1.4 is z > 1.2.<br>– l. 55: "conservatively off a kernel's centre", but c̄ = 0.02 is nominal (0.0500; 0.0092 at α = 0.01) [md]. The knot is conservative from c̄ = 0.032 on.<br>– l. 250: z's range is computed from the truth's fits. The draws' own fits reach s/\|ŵ\| = 0.0826 [own].<br>– l. 218–219 and l. 248–249 give two figures for the same part (0.012–0.098, mixed calibration; 0.0300–0.1075). Both carry the law's MC error; with shared samples it is 0.025–0.123 [own].<br>– E′'s margin column divides H1's μ by H0's sd. The null's μ gives 0.73 [own], which matches Φ⁻¹ of the 0.233 swap rate.<br>– l. 129: "(1 + ϑ²/6)^{−1/2} + O((Δ/σ)⁴)" hides a ϑ²(Δ/σ)² term; for (0.8, −0.4) at Δ = 0.0125 the remainder is −1.11·10⁻⁴ [own]. The final bound O(ϑ⁴ + (Δ/σ)⁴) is right.<br>– Status l. 140–143: the decomposition covers one pair and half a cell, not "two pairs, eight offsets".<br>– REVIEW_LOG l. 2576's "but for the rms" (D3bn-7) is still uncorrected. | Correct the text, the E′ column and the log. |

## Sources

- **Read in the export:**
  - in full: D52.md, `exp_close_features.py`, close_features.md, `test_close_features.py`, REVIEW_round67_D3bn.md, the Round 67 log entry, and the README diff;
  - in part: `rich_transit.py` (its function index) and `mirror_free.py` `iss_path` (l. 635–665).
- **Round 67's mutant strings.** I consulted round 67's own scratch mutant definitions to port the survivors verbatim.
- **Literature.** My one web search was refused (the turn's search budget was exhausted), and arxiv.org and ams.org are unreachable. Everything below is [memory]:
  - Harris, JOSA 54 (1964);
  - Helstrom, IEEE TIT 19 (1973);
  - Shahram–Milanfar, IEEE TIP 2004, TSP 2005, TIT 2006 (the locally optimal test, SNR ∝ Δ⁻⁴);
  - Smith, IEEE TSP 53(5) (2005);
  - Liu–Nehorai, IEEE TSP 55(11) (2007);
  - Neyman–Scott, Bull. ISI 41 (1966);
  - Kulaitis–Munk–Werner, Ann. Statist. 49 (2021);
  - Fukumizu–Amari, Neural Networks 13 (2000);
  - Liu–Wu–Wang, NeurIPS 2019.
  On these, the classical attributions look right. "New" is plausible for the transit's delay law, the knot's offset dependence, its allocation and its opposite-sign bias.

## Lesson

A difference between two tests, each calibrated by its own Monte Carlo law, is only as accurate as the two calibrations agree. E′ set the full knot against one 20 000-sample law and the runner-up's own knot against another, drawn from independent seeds. At the tie the two laws share their whole lower tail, yet their sampled quantiles differed by up to 3%, and that difference, not the dictionary, produced the "small gain" in 5 of 6 cells. Drawing both laws from the same samples is a one-seed change, and it makes the part exactly zero there. The effect entered through round 67's report and was adopted without checking it against the inequality that bounds it: the full knot is never later than the own knot, so equal tails leave no room for a gain. The same round's tests show the companion lesson. A test written to kill a mutant must execute the mutated function; checking the report on hand-made rows cannot see a broken score, width statistic or margin upstream.
