# Referee report D3cd, round 83: Proposition D‴.59 as revised after round 80

## Scope

- **What I reviewed.** I exported commit 3b6aac7 with `git archive` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/repo`. Below, W is the `D3cd_work` directory and W/repo the export. The files were:
  - `review/drafts/D59.md` (301 lines);
  - `experiments/exp_moving_width.py`, `results/chaos/moving_width.md` and `results/chaos/moving_width_raw.npz`, the last loaded with `exp_cohesion_leak.load`;
  - `tests/test_moving_width.py` (43 tests);
  - `src/core/moving_width.py`, and in `src/core/rich_transit_nd.py` the functions `loop_nd`, `sliding_iss_nd` (with `polish_block_nd`), `path_duals_nd`, `critical_rhs_nd` and `critical_jac_nd`;
  - `review/REVIEW_round80_D3ca.md`, the Round 80 entry of `review/REVIEW_LOG.md`, and `git diff cd71612 3b6aac7`.
- **Settings.**
  - Python: the venv's, with `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and `nice -n 10`.
  - pytest: `python -m pytest -q -p no:cacheprovider`, with `-x` for the mutants.
  - Every run of the experiment had `MW_OUT` set to a directory of mine holding a copy of the raw: W/mwout, or W/mut/<id>/out for the mutants.
  - I kept at most two heavy processes of mine at a time.
- **I modified nothing in the repository.** I ran only `git archive`, `log`, `diff` and `status` there.
  - During the review, other sessions committed 8c2c66f and 902fa4a (D‴.52, D‴.58).
  - At the end, D59.md, the script, the md, the raw, the test file and both src modules were byte-identical to my export (cmp), and `git status` was clean.
  - Everything I wrote is under W: scripts in W/scripts, logs in W/logs, pickles in W/out, mutants in W/mut.
- **Nothing of mine is still running** (ps at the end). I signalled no process.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my run, or my computation on the raw; [snippet] a web-search result; [memory] not verified this round.

## Summary

**What checks out**

- **Reproducibility** [own]:
  - `--from-raw` regenerates the md byte for byte.
  - The 43 tests pass in 12.6 s.
  - I reran all 15 A2 cells, all six GT rows, both GL runs and all four HLL cells (to τ = 10⁹). All are bit-identical to the raw in every stored field.
- **D3ca-1 is fixed.**
  - (ii)(c)(2), the proof and the Status now agree: an approach to an infimum only for K′ > 0, a retraction below the optimal dipole fit for K′ < 0.
  - My dissolution times equal the raw's: 144.9773, 154.1288 and 61.4426 for targets 18, 12 and 29.
  - κ³ is the same both ways: 1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)) = (t_diss − τ₄)/(t_diss − τ₃) = 0.42518, 0.75748 and 0.45487.
  - τ₄ comes at 0.675, 0.530 and 0.764 of t_diss. The κ guard is right [own].
- **D3ca-6 is fixed.**
  - CLOSURE_EXACT = 10⁻¹² lies 4.6 decades above target 20's closure rms (4.4·10⁻¹⁷) and 10 decades below the others' (1.5·10⁻² … 1.0·10⁻¹).
  - In the loop, target 20 forms no fourth kernel up to τ = 6400 at any ε.
  - Its κ′ goes 0.9915, 0.9954, 0.9976, 0.9988, 0.9994, 0.9997 (10⁻⁷ … 10⁻¹²).
  - Its landing comes 660, 164, 42, 11.5, 3.2 and 0.9 time units after the loop's τ₃ (≈ ε^0.55). So "lands right after τ₃" holds as ε → 0 [own, md].
- **τ₃ and τ₄ are the path's, by kernel index.** In the loop at small ε, the kernel that forms at τ₃ is the path's j3 (3, 4, 2, 4, 17, 3 for targets 9, 18, 12, 20, 27, 29). The kernel at τ₄ is the path's j4 (11, 15, 17, —, 3, 17). This holds on every path-sequence cell I ran (34 cells, ε = 10⁻⁵ … 10⁻¹²) [own].
- **The † rule flags exactly the cells on another sequence.** Those are 18 at 10⁻⁶ and 10⁻⁷, and 27 at 10⁻⁷ … 10⁻⁹. There kernel 9 forms at τ₃ instead of 4 or 17. The flags are the same for any SEQ_TOL in (0.11, 3.40) [own].
- **D3ca-4's algebra re-derives** [own, by hand]:
  - the centre of mass sits at c* + M₃/(6Wσ²);
  - its motion is −v³Δc²/(W³σ²) per unit Δc, against +v³Δc²/(W³σ²) from the shares at frozen v, using W − A = −v²Δc²/(2Wσ²);
  - the centres' share is (w_a² + w_b²)/W² to O(p⁴);
  - the arithmetic mean width moves by v²Δc²/(6W²σ), a third of the weighted mean's, so its share is p²/36;
  - C ≈ (Gm/(16μ₀))W²p²/(vσ̄²); |W − A|/|A| ≈ p²/8, and the offset is ≈ p³σ/24;
  - GL's p and the landing's differ by p²/72.
- **The parity claim holds.** For target 20's landed pair alone, with E minimised at fixed (Δc, v):
  - E(+Δc)/E(−Δc) − 1 = 4·10⁻¹¹ on a window symmetric about c*;
  - its local slope minus 8 is 0.355p², so there is no O(p) term;
  - on the samples the asymmetry is −2.7·10⁻³, a linear coefficient of about −0.011 against the target's 0.96 [own].
- **D3ca-7 is fixed.** The parity of the pair's sum in 1/v re-derives, and the test asserts the O(v⁻²) relative error.
- **D3ca-2's conclusion holds: on the symmetric pair the cancellation is the polish's.** I reran H's symmetric path [own]:
  - In block 8 (τ = 6968), one accepted polish stopped at its 5000-evaluation cap (status 0).
  - That polish multiplies max|w| by 3.358 and drops the residual from 2.05·10⁻² to 1.98·10⁻⁹. It reproduces H's end (4.26·10⁻¹⁰, |w| 3.239).
  - The loop's main kernels never exceed w = 1 on that pair.
- **D3ca-3's numbers check out.**
  - The plateau E is the same for every w, to six digits: 0.643650, 0.906294, 0.976497.
  - The escapes are as stated, at half-decade resolution.
  - 20A is still on the plateau at 10¹³; its w/A has fallen 20 → 18.6, much as 10A's had (9.97 → 8.45) one output before it escaped.
- **Numbers.** I traced every number this round changed or added in D59.md to the md or the raw and found no misquote.

**What is wrong (majors first)**

There is no major issue. Minors:

1. **D3cd-1: GL's remainder mixes two energies.** a = E/Δc⁸ takes all kernels' E, but the force factor takes the formed kernels'. The published top end, +3.2·10⁻⁴ at τ = 10⁷, is the idle floor (E_formed/E − 1 = 7.9·10⁻⁴ there). With one energy throughout, the remainder is −6.3·10⁻⁵ … +1.2·10⁻⁴.
2. **D3cd-2: `fit_split`'s criterion is fragile and incomplete.**
   - The c-only nearest rule misassigns two HL rows.
   - A one-sign split (H3's coinciding pairs) satisfies the stated cancellation criterion.
   - "Null cluster" is applied where the md's own test fails.
3. **D3cd-3: κ′ is not the loop's own κ, and the O(ε^{1/3}) rate is not shown for target 12.**
   - κ′ uses the path's τ₃, so at larger masses it is not the loop's κ.
   - At 10⁻¹² the loop's own κ is known only to the grid, and target 12's bracket excludes the path's κ.
   - Target 12's κ′ − κ falls only ×1.3–1.8 per decade.
4. **D3cd-4: A2's plateau is misnamed, and its route is roundoff-dependent.**
   - The plateau is the two-sample interpolation, not J_N's pointwise σ → 0 limit.
   - The collapse passes a roundoff-dependent one-sample saddle first.
   - SPIKE_A2 sits 15 % below a non-spike dip.
5. **D3cd-5: four behaviour-changing mutants survive.**
   - Three are in `hll_job` (HJ1–HJ3).
   - One is in `gl_job`'s E_formed (GJe), whose test assertion is vacuous.

Cosmetic: D3cd-6 (the † rule is a proxy for a direct test) and D3cd-7 (wording, rounding, literature).

## What I checked and how

1. **Reproducibility.**
   - `--from-raw` into W/mwout, compared with cmp; pytest on the test file.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/rerun.py` reruns `a2_job` on all 15 cells, `part_GT`, `part_GL` and `hll_job` on all four cells, and compares every field with the raw (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/logs/rerun.log`).
2. **GT (D3ca-1, -6).** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/gt_check.py` computes:
   - the dissolution times and κ³ two ways;
   - K's third order against G's;
   - the loop's κ three ways: uncorrected, κ′ (carried to the path's τ₃), and carried to the loop's own grid crossing.
3. **Sequences (D3ca-5).** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/seq_check.py`:
   - the path's τ₃ and τ₄ entrants (`kappa_from_path`'s internals, with the τ₄ runners-up);
   - the loop on TAUS_FINE, recording which kernel indices form, for 34 (target, ε) cells.
4. **Collapse (D3ca-4).**
   - By hand: the moment algebra above.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/gl_check.py`: the local slope and the remainder computed with E_formed, with E, and mixed.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/parity.py`: target 20's pair alone, E minimised over (W, c_a, σ_a, σ_b) at ±Δc, on the samples and on two windows symmetric about c*.
5. **A2 (D3ca-3).**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/a2_plateau.py`: the kernel's values at the samples on the plateau, against the one- and two-sample energies.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/a2_rows.py`: every stored output.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/a2_early.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/a2_early2.py`: dense output on s ≤ 10, and the same start under LSODA with spans 10, 10³ and 10¹³, and under Radau and BDF.
6. **(ii)(b) (D3ca-2).**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/split_check.py`: `fit_split`'s main kernels against the best-fitting assignment, on H, HP, H3, HL and HLL.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/hl_detail.py`: the kernels of a cell.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/h_anat.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/h_polish.py`: H's polish log, rerun.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/neg_sigma.py`: widths and the signs of the invisible kernels.
7. **Mutants.**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/mutants.py` holds 75 single-point mutants: 61 of this round's new or changed code, 14 of older code.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/run_mutants.py` runs each against the test file and against `--from-raw` (md cmp); results in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/logs/mut_summary.log`.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/surv_check.py` (KF10: GT rerun; HJ3: four HLL cells rerun) and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref83/D3cd_work/scripts/gje_check.py`: whether the survivors change results on the data.
8. **Literature.** Three web searches. I fetched no page.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) J, its single critical point and limits; convergence from E(start) < ½‖f*‖² (8–15; proof 110–118) | holds | unchanged; ∂_σE = ½π^{1/2}w(w − A) at (c*, σ*) re-derived [own] |
| (i) "sufficient, not sharp": 2.5A dips to 0.3σ*; from 3A on the plateau, E 0.644/0.906/0.976 whatever w; escapes at 3·10⁵, 3·10⁷, 1–3·10¹⁰; 20A on it at 10¹³ (15–21; A2 211–215; Status 161–163) | the numbers hold; the plateau is misnamed | 15 cells bit-identical [own]; the plateau is the two-sample interpolation; a one-sample saddle comes first (D3cd-4) |
| (i) J_N's limits on the samples (21–22) | holds | describes the one-sample limit, not A2's plateau (D3cd-4) |
| (ii) C's outcomes 71/3/6 (24–26; 216–221) | holds | [md] |
| (ii)(a) one Gaussian to O(Δc⁴), the width split (27–32; 118–124) | holds | unchanged [md] |
| (ii)(b) the path's cancellation, H2, H3 (33–41; 261–272, 278–280) | holds | [md] |
| (ii)(b) the fit's split tells a cancellation from a null cluster (41–43) | holds for every number D59 quotes; incomplete as a criterion | misassigns two HL rows; a one-sign split meets it (H3) [own] (D3cd-2) |
| (ii)(b) asymmetric pairs: HP; HL at 10⁻⁸ and 10⁻⁶; HLL to 10⁹, undecided (43–51; 272–292) | holds | HLL bit-identical; every number traced [own, md] |
| (ii)(b) symmetric pair: the polish's; decay at 10⁻⁶; null cluster at 10⁻⁸; HP the same state (51–58; 286–290; Status 166–168) | holds; "null cluster" is an extrapolation | one capped polish ×3.358 [own]; the md's null-cluster test is not met at τ ≤ 10⁹ (D3cd-2) |
| (ii)(c) closure, v³ = 3μ₀Ks, K from the third order, O(v⁻⁴) (59–66; 124–134) | holds | parity in 1/v re-derived [own] |
| (ii)(c)(1) τ₃ from the path's duals (66–72; 134–140) | holds | GT bit-identical; the loop's τ₃ entrant is j3 [own] |
| (ii)(c)(2) K′ of either sign, retraction, dissolution, τ₄; the exact closure (72–83; 140–145; 240–244) | holds | t_diss and κ³ [own]; the τ₄ entrant is j4 [own] |
| (ii)(c)(3) κ³ from the path; the loop's κ′ converging; O(ε^{1/3}) a bound; † (83–90; 244–252) | κ_path and † hold; κ′ and the rate are overstated | κ′ uses the path's τ₃; target 12's rate is ε^{0.12–0.26} (D3cd-3) |
| (ii)(c)(4) the collapse law, C, the metric's shares, the force factor, parity (90–102; 145–159) | holds | by hand; parity 4·10⁻¹¹ on a symmetric window [own] |
| GL: exponents, ratios, local slope 0.951–0.963, remainder −4.7·10⁻⁵ … +3.2·10⁻⁴, p²/36 ≤ 0.00067 (253–261) | holds except the remainder's range | bit-identical; the remainder mixes two energies (D3cd-1) |
| (iii), E, F (104–108; 292–301) | holds | unchanged [md] |
| Status (161–176) | fits, with two edits | "to 3.2·10⁻⁴" (D3cd-1); "told from a null cluster by the fit's split" (D3cd-2) |
| Classical, literature (176–198) | fits | [snippet] below |
| New (198–204) | fits; "asymmetric pairs" generalises | 10⁹ reached on two of five asymmetric cells (D3cd-7) |

## Status of the round-80 issues

| D3ca | status | what remains |
|---|---|---|
| 1 (K′ < 0) | fixed | Nothing material. The wording is consistent in (2), the proof and the Status. My t_diss equals the raw's. Tests kill KF5 and KF7; KF4 (the k1 > 0 guard dropped) is equivalent on the data. "145" is 144.98 printed by `.4g` (cosmetic). |
| 2 (symmetric pair) | fixed in substance | The conclusion is confirmed (the capped polish). The new criterion is fragile and incomplete, and "null cluster" is used beyond the md's own test (D3cd-2). |
| 3 (A2's 10A) | fixed | The plateau's identity, the one-sample saddle and SPIKE_A2's margin (D3cd-4). |
| 4 (finite-p) | fixed, one inconsistency | The derivations hold, and parity is verified. GL's remainder mixes two energies (D3cd-1). |
| 5 (the loop's κ) | partly fixed | The grid offset is removed and † is right. κ′ uses the path's τ₃; the grid bracket is missing; target 12's rate (D3cd-3). |
| 6 (target 20) | fixed | Nothing. The threshold is well separated; K′ = 0, no τ₄ and κ = 1 are borne out by the loop. |
| 7 (O(v⁻⁴)) | fixed | Nothing. |
| 8 (tests) | mostly fixed | All ten round-80 survivors are now killed (my GT1, FPf, FPg, FPh, GJd, HS3, HLa, A2b, A2d). HJ1–HJ3 and GJe survive (D3cd-5). |
| 9 (literature) | fixed | The secondary source can be named (D3cd-7). |
| 10 (wording) | fixed | Nothing. |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cd-1 | minor | **GL's remainder mixes two energies (D59.md:257–258; Status :174; `gl_lines`).** <br>• **The mix.** The ratio with C uses a = E/Δc⁸ with E over all kernels (`gl_constant`). The force factor (d ln E/d ln Δc)/8 uses E_formed. <br>• **The floor.** Along the run, E_formed/E − 1 rises from 1.0·10⁻⁶ to 7.87·10⁻⁴. <br>• **What follows.** The published remainder drifts from −4.7·10⁻⁵ (τ = 4.4·10⁴) to +3.22·10⁻⁴ (τ = 10⁷), tracking that floor. With one energy throughout it is −6.3·10⁻⁵ … +1.2·10⁻⁴ (all kernels) or −3.8·10⁻⁴ … +1.2·10⁻⁴ (formed only) [own]. <br>• **The coefficient.** The local slope's 0.951–0.963 is E_formed's; from all kernels' E it is 0.951–1.052. <br>• **What survives.** The conclusion, a remainder ≲ 4·10⁻⁴, holds either way. | **GL and the Status:** <br>• compute the remainder with the total E in both factors (the force is the total energy's slope), keeping E_formed's slope only for the coefficient 0.961; <br>• quote "1 within −6.3·10⁻⁵ … +1.2·10⁻⁴" at :258 and :174. |
| D3cd-2 | minor | **`fit_split`'s criterion (D59.md:41–43, :52, :166–167, :274).** <br>• **(a) The c-only nearest rule misassigns.** In HL, A₂ 0.99 at ε = 10⁻⁶, at τ = 6400 and 2.5·10⁴, it picks a +0.071 kernel at c = −0.042 (0.038 from the feature). The +1.58 kernel at −0.123 (0.043 away) is the one that carries the feature. Main error 0.93 instead of 0.71/0.65; the others' Σ\|w\| 2.86/2.88 instead of 1.35/1.32. At τ = 10⁵ the choice flips (1.31). <br>  – A (c, log σ) rule (mutant FS3) changes these two rows and HL A₂ 0.9 at 10⁻⁸, τ = 6400, but no number D59 quotes [own]. <br>• **(b) A one-sign split meets the criterion.** The stated test ("a main kernel off by O(1), the others summing to O(1)") is also met by H3's coinciding pairs: 0.38 with the others at 0.38 (±0.15); 0.37, 0.37 with 0.37 (±0.2) [md]. <br>  – \|w/A − 1\| discards the sign that tells them apart: +2.2 for H's inflated kernel, −0.38 for a split's deficient one. <br>• **(c) "Null cluster" is used where the md's own test fails.** The test is "small against the residual". <br>  – HP's end: the others sum to 8.7·10⁻² of max f* against a residual of 9.7·10⁻⁶. <br>  – HLL at 10⁻⁸, τ = 10⁹: 2.2·10⁻⁵ against 5.4·10⁻⁹. The others' sum falls by ×16, 8.4, 6.3, 4.5 per decade, decelerating [md]. | **`fit_split`:** <br>• report the signed w/A − 1; <br>• assign main kernels by (c, log σ) or by the smallest error. <br>**(ii)(b):** <br>• state the criterion with the sign: a cancellation has an inflated main kernel plus others of the other sign; a split has a deficient one; <br>• at :52, :167 and :274 write "approaching a null cluster (HLL)", or give the test's numbers. |
| D3cd-3 | minor | **κ′ and the O(ε^{1/3}) rate (D59.md:89–90, :245–252; `kappa_loop`).** <br>• **(a) κ′ is not the loop's κ at larger masses.** κ′ replaces the loop's own crossing by the path's τ₃. Target 12, κ′ against the own crossing (the grid step's midpoint): 0.929 vs 0.965 (10⁻⁷), 0.925 vs 0.943, 0.921 vs 0.931, 0.918 vs 0.922 (10⁻¹⁰) [own]. <br>• **(b) At 10⁻¹² the loop's own κ is known only to the grid.** Bracketed between "v not carried" and "v carried to the grid's τ₃": [1.1126, 1.1159], [0.7498, 0.7525], [0.9139, 0.9148], [0.9971, 1.0056], [1.0187, 1.0209], [0.7606, 0.7694]. <br>  – The path's κ lies inside for five targets. Target 12's 0.9116 is 0.25–0.35 % below its bracket. <br>  – The quoted "−0.27 … +0.21 %" borrows the path's τ₃. <br>• **(c) The rate.** Target 12's κ′ − κ is +0.0173, +0.0131, +0.0097, +0.0063, +0.0036, +0.0020 (10⁻⁷ … 10⁻¹²): ×1.32–1.80 per decade, with (κ′ − κ)/ε^{1/3} rising 3.7 → 20. It is pre-asymptotic. Targets 9 and 29 show ε^{1/3} (×2.15–2.63; ×1.9–2.3). | **(3) and GT:** <br>• call κ′ "v_land over the law's v at the path's τ₃"; <br>• add the loop's own-crossing bracket at 10⁻¹²; <br>• "converges, O(ε^{1/3}) asymptotically; over 10⁻⁷ … 10⁻¹² target 12 converges only as ε^{0.1–0.3}". |
| D3cd-4 | minor | **A2 (D59.md:18–21, :162, :211–215; `a2_spike`).** <br>• **(a) The plateau is the two-sample interpolation.** The kernel equals f* at the two samples nearest c* (to 10⁻⁶) and is ≤ 3·10⁻⁶ elsewhere. So E = 1 − (f*(x_j)² + f*(x_{j±1})²)/Σf*² = 0.643650, 0.906294, 0.976497, which the plateau matches to six digits. <br>  – J_N's σ → 0 limit at a fixed c is the one-sample value, 0.8119, 0.9530, 0.9882; that is the limit :21–22 describes. The plateau is lim_{σ→0} min_c J_N [own]. <br>• **(b) A one-sample saddle comes first.** c* = 0.03 is a sample, and the collapse first passes a symmetric state at c = c*: for 5A at 0.25, σ/σ* = 0.1846 and E = 3.480 at s = 10⁻³. <br>  – Its exit depends on roundoff and solver. LSODA on [0, 10] stays there to s = 10 (E 3.406). LSODA on [0, 10³] or [0, 10¹³] (A2's span), and Radau or BDF on [0, 10], reach the plateau by s = 0.01 [own]. <br>• **(c) SPIKE_A2 has a thin margin.** 0.2 lies 15 % below the 3A/0.25 dip (0.2345). At 0.24 (mutant AS6) that dip reads as a spike with E 1.17–1.20, which is not the plateau. The plateau's widths are ≤ 0.094. <br>• **(d) The dip is at the outputs.** "Dips to 0.295–0.304σ*" is the outputs' minimum; the flow's is 0.293–0.303, reached before s = 1. | **(i), the Status and A2:** <br>• "the two samples nearest c* (lim_{σ→0} min_c J_N)"; <br>• note the one-sample saddle, or start A2's feature between samples; <br>• define the spike by E within 10⁻³ of the two-sample value; <br>• add "(at the outputs)" at :212. |
| D3cd-5 | minor | **Tests (D3ca-8).** Of 75 mutants, the tests kill 48. `--from-raw` alone kills 9 more (KF1, FS3, FS8, GL5, GL6, AS6, GR5, AG1, DP1). 18 survive. <br>• **Equivalent on the data (14):** KF4, KF10 (GT rerun identical), KL2, SQ1, SQ5, FS5, FS10 (no negative widths among HL/HLL kernels), AS1, AS2, AS3, AS5, CF2, GR4, LR1 (all 140 invisible formed kernels negative). <br>• **Behaviour-changing (4):** <br>  – HJ1: `hll_job` runs every cell at EPS_HL[0]. The test's cell has ε = 10⁻⁶ = EPS_HL[0]. <br>  – HJ2: outputs computed at TAUS_HL but labelled TAUS_HLL. The test mocks `loop_nd`. <br>  – HJ3: rtol 10⁻⁶. HLL changes by ≤ 3·10⁻³ (A₂ 0.95's residual at 10⁹: 6.215 vs 6.196·10⁻⁹), below the printed digits. <br>  – GJe: E_formed computed over all kernels. The test's idle kernel at c = 0.5 changes E by 1.9·10⁻⁸ relative, which passes np.isclose's rtol of 10⁻⁵. On the data the local slope would read 0.951–1.052 [own]. | **Tests:** <br>• call `hll_job` at ε = 10⁻⁸ and assert the times and rtol/atol passed to `loop_nd`; <br>• in `test_gl_job_and_its_constant`, put the idle kernel on the feature and compare E_formed with rtol 10⁻¹², atol 0. |
| D3cd-6 | cosmetic | **The † rule is a proxy.** The loop's states hold every kernel, so the τ₃ entrant's index can be compared with the path's j3 directly. <br>• On 34 cells, † marks exactly the cells with a landing on another sequence (kernel 9 instead of 4 or 17) [own]. <br>• 18 at 10⁻⁶ prints "—†" in the table but "—" in the κ′ list. | Compare entrant indices, or report both; print † on "—" in the list too. |
| D3cd-7 | cosmetic | **Wording, rounding, literature.** <br>• "145" is 144.98 printed by `.4g` (:81, :242). <br>• "E even in Δc" holds in the continuum; on the samples it holds to 2.7·10⁻³ (:100, :158). <br>• New's "the loop keeps on asymmetric pairs (to τ = 10⁹)" (:203) rests on two of five asymmetric cells; the others reach 10⁵. <br>• The symmetric "decay" at 10⁻⁶ (×0.73 per decade) needs about eight more decades to reach 1 % of max f*: a trend at the horizon. <br>• Ho–Nguyen's r̄(2) = 4 has a nameable secondary source (arXiv:2403.01684) [snippet]. | Edit the text. |

## Sources

Three web searches; I fetched no page, and arxiv.org is unreachable.

- **Ho, N. and Nguyen, X.**, "Convergence rates of parameter estimation for some weakly identifiable finite mixtures", Ann. Statist. 44(6) (2016) 2726–2755, [Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-44/issue-6/Convergence-rates-of-parameter-estimation-for-some-weakly-identifiable-finite/10.1214/16-AOS1444.full). For location–scale mixtures over-fitted by one component the rate is n^{−1/8} under W₄, and n^{−1/12} under W₆ for two [snippet].
- **"Dendrogram of mixing measures"**, [arXiv:2403.01684](https://arxiv.org/pdf/2403.01684). It attributes r̄(2) = 4, r̄(3) = 6 and r̄(4) ≥ 7 to Ho–Nguyen 2016a [snippet]. The companion report is [arXiv:1501.02497](https://arxiv.org/pdf/1501.02497) [snippet].
- **Gallon, Jentzen and Lindner**, "Blow up phenomena for gradient descent optimization methods in the training of artificial neural networks", [arXiv:2211.15641](https://arxiv.org/pdf/2211.15641) [snippet].
- **[memory], not re-verified:** de Boor 1969; Kammler; Braess; Petersen–Raslan–Voigtlaender 2021; Xu–Du 2023; Dwivedi et al. 2020; Chizat 2022; Pearson; Runnalls; Schiebinger–Robeva–Recht 2018; Duval–Peyré 2017; Golub–Pereyra; Poggio–Girosi; Keriven et al.; Boyd–Schiebinger–Recht; Denoyelle et al.; Kiefer–Wolfowitz; Day; Biernacki–Chrétien; Marr–Hildreth.

## Lesson

This round's fixes are right where they compute something new: the dissolution times, the exact closure, the plateau's escape and the collapse metric all re-derive or rerun exactly. The weak spots are in the comparisons built to test them.
- κ′ measures the loop's κ with the path's τ₃.
- GL's remainder divides an all-kernel constant by a formed-kernel slope.
- `fit_split`'s main kernel is whichever kernel lies nearest in c, and its error drops the sign.

Each choice made a check look tighter or cleaner than the data alone allow. The direct versions were cheap and unambiguous: the loop's own crossing bracket, one energy throughout, the kernels' indices, the sign of w/A − 1. A check should not borrow the quantity it tests.
