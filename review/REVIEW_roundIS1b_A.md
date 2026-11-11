# Referee round IS1b-A: crowded-field astronomy test, transfer test and post hoc kernel study

Repository `/home/user/Lyapunov-Convergence-Kernels-Showcase`, branch `claude/lyapunov-context-setup-yg51ew`, reviewed at commit eee967e (clone in my scratch; every astronomy file I use is identical to HEAD). My scripts and logs are in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_IS1b_A_work/` (`myscore.py`, `rederive1-3.py`, `robust_metric.py`, `align_split.py`, `split3.py`, `adequate_check.py`, `diag_resid.py`, `synth_*.py`, `stack_psf*.py`, `cs_sens.py`, `repro_path.py`, `logs/`).

**Rule compliance.** I did not modify the repository, did not use the network, and posted or pushed nothing. I did not read `astro_data_3`, the *repeat* files or the other referee's directory. Every run used `OMP_NUM_THREADS=1 nice -n 10`, and I killed no process. I kept to two long path runs at a time; short analysis scripts of seconds to about 3 minutes sometimes overlapped, so up to three of my processes ran together on a few occasions. I read the coordinator's scratch run logs (`run_astro_T.log`, `astro_D2.log`, `astro_smoke2.log`) as provenance data. During my review the original working tree gained uncommitted edits by someone else (`HANDOFF.md`, `tests/test_imaging.py`, `tools/mutants/wp2_imaging_own.py`) and more untracked *repeat* files. None of them touches the astronomy files I reviewed.

## 0. Bottom line

- **No blocker.** I found no evidence of test peeking. The freeze, its hashes, the commit order and the run order are consistent.
- **Numbers.** Every headline number reproduces exactly with my own scorer. I re-ran the frozen path on one test patch and got a bit-identical catalogue.
- **Weakest part: the explanation of the M92 failure (claims 3 and 4).** I built synthetic frames from each real patch's own psField stamp and the HUGS list, with an exact PSF, exact noise model and exact truth.
  - The kernel misfit costs the ISS its lead on M92-like data (a tie), but does not make it lose.
  - The bright-star splitting is not cured by the better kernel, and it also occurs on M2.
  - The real M92 and M53 stamps do not match the real stars (central pixel -16 % and +10 %, against +2 % on M2).
- **Lead over crowdsource.** It is not attributable to the inverse-scale-space ordering. A greedy order ties on the one patch where it was run.
- **Smaller issues.**
  - Quoted intervals are conditional on the patches.
  - The frozen configuration sits on the edge of all three grids; only kappa is flagged.
  - The 2500-entry cap censors runs without a flag.
  - F at m < 22 is blind to the faint units.
  - Several wording points.

## A. Findings

### 1. MAJOR (claims). The kernel-misfit explanation does not account for the M92 reversal, and its stated mechanism (split bright stars, "optimistic noise level") is contradicted by control experiments

**Where.** `experiments/exp_astro_psf.py:5-14` (docstring "What explains it ... the cause that the algorithm can address"), `results/astro/psf.md` (intro, bright-star section, "Paired differences" tables with verdicts), claims 3 and 4.

**Method.**
- For each patch I rendered every HUGS source (F606W < 30, box plus 20 px) with that patch's own 51x51 psField stamp, using cubic-spline sub-pixel shifts.
- Fluxes come from the HUGS magnitudes. I added Poisson noise at the patch gain and Gaussian dark noise at its dark variance, and used the pipeline sky.
- Scoring uses the injected sources as truth. I ran the repository's unchanged `run_path` (frozen isotropic kernel, sigma_b 3, r_ex 1.2, kappa 2), `run_path_psf` (adequate kernel) and the crowdsource default.
- The synthetic frames resemble the real ones: median 284 against 323 DN on M92, 5 and 95 percentiles within 15 %.

**Evidence: F at m < 22, real against synthetic (kappa 2, M2-frozen configuration).**

| patch | real iso | real adequate | real crowdsource | synthetic iso | synthetic adequate | synthetic crowdsource |
|---|---|---|---|---|---|---|
| M2 (310,630), control | 0.653 | (same kernel) | 0.377 | 0.589 | (same kernel) | 0.386 |
| M53 (660,420) | 0.199 | 0.330 | 0.157 | 0.197 | 0.291 | 0.156 |
| M92 (1100,220) | 0.148 | 0.221 | 0.323 | 0.306 | 0.370 | 0.359 |
| M92 (1430,550) | 0.104 | 0.316 | 0.160 | 0.232 | 0.399 | 0.185 |
| M92 pooled | 0.122 | 0.277 | 0.226 | 0.262 | 0.387 | 0.256 |

- **Control and M53.** The M2 control reproduces crowdsource (0.386 against 0.377), and the ISS is within 0.06. M53 (660,420) is reproduced for both methods: 0.197 against 0.199 and 0.156 against 0.157.
- **M92 in the exact world.** Kernel misfit takes the ISS from a lead of about +0.20 (M2-like) to a tie (0.262 against 0.256 pooled). In the exact world it does not produce a loss.
- **Real M92.** The real data cost the isotropic ISS a further -0.14 and crowdsource only -0.03. That extra loss is not produced by kernel shape, a wrong noise level or a PSF 6 % too narrow taken alone.
  - Noise variance x1.5 unknown to the methods: ISS 0.223 against 0.232, so the ISS is insensitive.
  - Injected PSF 6 % broader (FWHM +5 to 8 %): ISS 0.195, crowdsource 0.168. Still a win.
- **Combined small systematics.** A PSF 6 % broader, a (0.18,-0.06) px truth offset, 15-30 % of the truth missing and 0.15-0.25 mag truth scatter, applied together, give ISS 0.131-0.156 against crowdsource 0.133-0.146, a tie. This approaches the real -0.056 and -0.176 on the two patches.
  - The ISS has many faint detections near the positional and magnitude tolerance. Its margin in the exact world is thin.
- **Residual scale.** The robust sd of the standardized residual (kappa 2) in the exact world is 2.19 for the isotropic kernel and 1.13 for the adequate kernel. The real values are 2.28 and 1.27. The M2 control gives 0.97 against 0.96, and M53 gives 1.55 and 1.14 against 1.55 and 1.17.
  - The excess is therefore reproduced by kernel misfit with an exact noise model. "Optimistic noise level" is the wrong description: the noise model is exact in the synthetic world and the sd is still 2.2.
- **Splitting is not a kernel symptom.**
  - Definition: bright isolated stars have truth m < 20.5 and no truth source within 3 px brighter than m_star + 2.5. A split means at least 2 units of at least 10 % of the star's flux within 2 px. n = 136, 116, 82 for M2, M53, M92.

  | catalogue | M2 | M53 | M92 |
  |---|---|---|---|
  | real ISS tuned, isotropic | 0.40 | 0.59 | 0.90 |
  | real ISS default, isotropic | 0.22 | 0.41 | 0.55 |
  | real ISS tuned, adequate kernel | - | 0.78 | 0.68 |
  | real ISS default, adequate kernel | - | 0.59 | 0.59 |
  | crowdsource | 0.04 | 0.00 | 0.02 |
  | photutils tuned | 0.10 | 0.05 | 0.01 |

  - Exact-stamp synthetic, M92 (1430,550): isotropic 0.59, adequate 0.73, crowdsource 0.00. M2 control: isotropic 0.36.
  - Splitting exists on M2 with a kernel fitted to 0.1 % of the peak. It is not reduced by the post hoc kernel; the M92 adequate-kernel value (0.68, 0.59) is about the same as the isotropic default (0.55) and the exact-world M92 value (0.59), and on M53 and in the synthetic world it is higher.
  - The parenthesis in `exp_astro_psf.py:295` ("a kernel that does not fit makes extra units") is not supported.
- **M92 deficit is in the flux, not the position.**
  - With position-only matching the pooled M92 gap shrinks from -0.104 to -0.018 (ISS tuned 0.280, crowdsource 0.298).
  - For stars of truth m < 19.5, the median (detection - truth) magnitude is +0.29 and +0.32 for the isotropic ISS on the two patches, +0.38 and +0.31 for the adequate kernel, and +0.09 and +0.13 for crowdsource. The exact synthetic world gives +0.04 and +0.12 for the ISS.

**Fix.**
- Say that the kernel family removes the residual excess and the positional bias and, post hoc, recovers a lead on the same patches. Do not say it explains the transfer failure or the splitting.
- Either include the exact-stamp synthetic control or drop the causal paragraph.
- Name what was not identified: the real M92 data reduce the ISS further, and the measured candidates are in findings 2 and 8.

### 2. MAJOR (data). The M53 and M92 psField stamps were validated as a computation, not against the stars, and they differ from the stars by amounts M2 does not show

**Where.** `results/astro/second_cluster_spec.md` §5.3 (agreement with pydl) and §8 ("a check of the psField PSF against stars of the new frames" is listed as not done). Brief question 4.

**Evidence.**
- I recomputed all eight stamps (M2, M53, M92) with pydl `sdss_psf_recon`. After clipping and unit-sum normalisation they match the stored stamps to at most 2.1e-8.
- I recomputed `dn` for M92 (1100,220), M92 (1430,550) and M53 (660,420) from the frame files and got exact equality. I checked the sha256 of all six second-cluster patch files and of the three truth files against the spec and `freeze.json`.
- The raw psField reconstruction sums to 1.11-1.15 on M53 and about 1.01-1.02 on M92. The clipped negative mass is 2.0-2.5 % and 1.2-1.4 %. All methods share the clipped stamp.
- I stacked isolated bright stars with the median of own-centroid-centred cutouts. Selection: m < 20, no truth source within 4.5 px brighter than m + 1.5. I applied the identical procedure to the synthetic frames built from the stamp, so centroiding blur cancels.

  | frame | stars | central pixel, real minus stamp | note |
  |---|---|---|---|
  | M2 | 18 | +2 % | |
  | M53 | 23 | +10 % | 68 % interval [+1, +11] |
  | M92 | 39 | -16 % | 68 % interval [-16, -10]; FWHM 3.35 against 3.12 px |

  - The M92 stack reproduces the stamp's lumpy asymmetry to within about 4 % per core pixel. The lumps are real; the peak is not.
  - On M92 the 6 %-broadened synthetic reproduces the faint photometric bias of the real catalogue.
  - The psField came from 14 and 16 PSF stars on M53 and M92.

**Consequence.** The second-cluster comparison tests the psField reconstruction as much as the method. All methods share the stamp, but the ISS's rigid analytic kernel with a high-pass is plausibly more exposed (not isolated by my runs).

**Fix.**
- Measure the PSF from the field's own stars, or give each method its own refit.
- State in the limitations that the stamps were not validated against stars, and that on M92 the stars are less peaked than the stamp by about 16 %.
- Crowdsource's own `refit_psf=True` does not help it (finding 9), so the asymmetry is not simply a missing option of the baseline.

### 3. MAJOR (wording). The lead is not shown to come from the inverse-scale-space ordering

**Where.** The design-variant table in `m2.md` (lines 220-231, benchmark patch); `exp_astro_m2.py:250-252`.

**Evidence.**
- The greedy order with the same stop (instantaneous significance) gives F 0.644 [0.608, 0.686], dF against the frozen path -0.009 [-0.029, +0.012]. It took 1152 entrants and 1214 s, against 294 s for the frozen path.
- The dual order with that same stop gives 0.412. The stop rule matters; the order does not.
- No high-pass gives 0.260, below crowdsource's 0.377 on that patch (kappa not re-chosen). No exclusion zone gives 0.603, from a path that hit the 2500 cap. The local background removal at sigma_b = 3 px is the dominant component.
- The replacement of the dropped greedy-with-time-stop row is fair in construction: both orders now share a stop. The dropped row was degenerate by construction, because `sliding_path` honours `tau_max` only in the dual branch (`imaging.py:335-339`), so the greedy path ignored it and ran to the cap. The replacement was made after seeing the first run on a test patch, and the history is recorded only in a code comment, not in `m2.md`.

**Fix.**
- State that the dual order buys speed, about 4x, and not F.
- Attribute the lead to local background removal, the exclusion zone and the time rule.
- Re-run the variants on the development patches with kappa re-chosen per variant.

### 4. MINOR (statistics). The quoted intervals are conditional on the patches

**Where.** Pooled bootstrap, `exp_astro_m2.py:444-451` and the report.

**Evidence.**
- The pooled bootstrap resamples 20 px blocks within patches. `m2.md` already says the patches are not independent trials.
- Per-patch F of the tuned ISS on M2 is (0.653, 0.417, 0.567, 0.492), sd 0.10. The patch-level 95 % interval of the mean is [0.37, 0.69], against the quoted [0.490, 0.535].
- The paired differences are robust. M2 tuned minus crowdsource per patch is (+0.277, +0.276, +0.195, +0.279), patch-level t interval [+0.19, +0.32].
- The M53 small lead is fragile. Tuned minus crowdsource per patch is (+0.028, +0.008, +0.041, +0.041); the verdicts are tie, tie, win, win. The pooled interval is [+0.010, +0.051] and the patch-level t interval is [+0.004, +0.055], with correlated patches. Crowdsource's dev-best variant (finding 9) shrinks the lead to +0.020.
- M92 has two adjacent windows of one frame and one stamp. Both show a loss, but no patch-level inference is possible: the t interval is [-0.88, +0.65].

**Fix.** Call the interval conditional on the patches, or add the patch-level interval, and say "two windows" for M92.

### 5. MINOR (protocol disclosure). The frozen configuration is a boundary optimum on all three grid axes, and it does not transfer

**Where.** `exp_astro_m2.py:55-57`; `m2.md` line 60 flags only kappa (report code lines 518-519).

**Evidence.**
- sigma_b = 3 is the lowest of {3, 5, 8}, r_ex = 1.2 is the highest of {0.8, 1.2}, and kappa = 2 is the lowest of the kappa grid.
- Tuned against default pooled F: M2 test +0.043 (0.511 against 0.468), M53 -0.047 (0.263 against 0.310), M92 -0.028 (0.122 against 0.150). The permissive edge hurts out of distribution.
- The tables "The path at the other significance levels" on the test patches are not labelled descriptive-only.

**Fix.**
- Flag all three edges, or extend the grids.
- Report the a priori default as the primary transfer row.
- Label the kappa tables on the test patches as descriptive.

### 6. MINOR (disclosure). The 2500-entry cap censors runs without a flag

**Where.** `m2.md` prints the cap only for the frozen configuration's run (`exp_astro_m2.py:551-555`); the flag `bic[2]` is in the raw.

**Evidence.**
- 6 of the 12 development runs are capped: (sigma_b 5, r_ex 0.8) and (8, 0.8) and (8, 1.2), on both patches. Their rows in the development table are identical for several kappa. For (5, 0.8), kappa 2 to 4 give 0.469, 0.577, 0.517 and 1871 detections; for sigma_b = 8 all kappa are nearly identical (0.364 for (8, 0.8); 0.402 for (8, 1.2)). They are the last block of the capped path, not F at that kappa. The frozen choice is an uncapped run and is not affected.
- Both M92 a priori default runs are capped (`default|bic` = [293, 2501, 1] and [630, 2501, 1]).
  - The M92 default rows (kappa 5 and BIC) are identical in each patch, so they probably come from blocks far below the cap.
  - The path had a flood of entries at sigma_b = 5 on M92 and was not run to its stopping time.
- The "no exclusion zone" variant (2500 entrants) is capped, so its row is truncated too.

**Fix.** Print the cap flag wherever it is set, and mark censored development rows.

### 7. MINOR (metric). F at m < 22 is blind to faint units and hides over-detection on the wide-PSF frames

**Evidence.**
- Purity counts only detections brighter than 22 mag. In the variant table the no-clean-up row shows dF 0.000 [0.000, 0.000], yet it carries 1073 detections against 1041. The 32 extra units are all at 22 mag or fainter: the sub-22 detection count is 855 in both rows.
- Pooled all-magnitude purity (matched over detected):

  | set | ISS tuned | crowdsource | photutils tuned |
  |---|---|---|---|
  | M2 | 0.533 | 0.477 | 0.468 |
  | M53 | 0.239 | 0.399 | 0.441 |
  | M92 | 0.112 | 0.425 | 0.576 |

  - Unmatched detections on M92: 1899 for the tuned ISS and 351 for crowdsource.
- On M2 the ranking is robust. The ISS is ahead of crowdsource and photutils in every variant I computed:
  - depth 21, 21.5, 22.5 and 23;
  - tolerance 0.75 px, and tolerance 1.0 mag;
  - position-only matching;
  - truth magnitudes shifted by +-0.1 and +-0.25.
- The 22-mag edge has no net asymmetry: pooled matched truth 2485 against matched detections 2480 for the tuned ISS on M2.

**Fix.** Report purity over all detections, or the number of unmatched detections brighter and fainter than 22, beside F. Do not read the ablation rows as "no effect" for faint-unit control.

### 8. MINOR (astrometry). The ISS carries an extra position offset on M53 and M92; it is a partial contributor to the M92 gap, not the explanation

**Evidence.**
- Median offsets (detection minus truth) of bright pairs, truth m < 19.5, position-only matching. The spec's image-based truth-to-image shift on M92 is (+0.134,-0.079) and (+0.182,-0.063), uncorrected.

  | method | M92 (1100,220) | M92 (1430,550) | M2 benchmark patch |
  |---|---|---|---|
  | ISS tuned | (+0.26,-0.26) | (+0.27,-0.21) | (+0.09,+0.02) |
  | crowdsource | (+0.18,-0.11) | (+0.16,-0.07) | |
  | photutils | (+0.15,-0.09) | (+0.15,-0.07) | |

  - On M53 the ISS has an extra y offset of about +0.13. M2 shows none.
- Removing each method's own median offset: M92 F goes 0.122 to 0.167 for the ISS and 0.226 to 0.242 for crowdsource. The gap goes from -0.104 to -0.075.
- Applying the spec's method-independent image-based shift to the truth: M92 F goes 0.122 to 0.147 for the ISS (tuned), 0.150 to 0.185 for the default, 0.226 to 0.240 for crowdsource and 0.166 to 0.179 for photutils. M53 barely moves.
- In the synthetic exact world the isotropic kernel gives its own offset (+0.04,-0.10) and the adequate kernel (-0.03,-0.01). A (0.18,-0.06) px offset costs the isotropic ISS 0.028 F and crowdsource 0.003, because the permissive ISS catalogue is more tolerance-sensitive.
- The stamp's flux centroid (-0.04,+0.04) is not the cause of the ISS offset.

**Fix.** Apply the measured method-independent truth shift and report both, or state that the declared uncorrected offset is worth about 0.01-0.03 of the M92 gap.

### 9. MINOR (baselines). The baselines are close to their own plateaus; some wording is inexact

**Evidence.**
- **Crowdsource "best of 48 post-filters is its unfiltered output".** It is a six-way tie at F = 0.290651. The selected name (the lexicographically largest, `snr5|ff0|qf0|spinf`) drops 2-3 sources per test patch.
- **Crowdsource internal knobs.** I ran 16 variants on the development patches only: threshold 3 and 8, `psfderiv=False`, blendthreshu, psfvalsharpcutfac, maxiter 20, sky-grid npix 6, 10 and 14 (monkeypatched), the ISS's smoothed variance weights, and combinations.
  - Dev F: default 0.291, best (smoothed weights plus sky 10) 0.305, `psfderiv=False` 0.167.
  - That best variant on the test sets changes pooled F by +0.008 (M2, 0.251 to 0.259), +0.010 (M53, 0.233 to 0.243) and -0.005 (M92, 0.226 to 0.221). No ranking changes; the M53 lead of the tuned ISS shrinks from +0.031 to +0.020.
  - Threshold 3 adds 3 detections and changes nothing in F. `refit_psf=True` is worse on the M2 test patches (0.186) and on M92 (0.106) and +0.014 on M53.
- **Photutils.**
  - The best variant is the last extension added (`ext_all_local_i6`, dev F 0.236, against 0.140 for the best of the original 12), so the frontier of the family was not established.
  - The M2 ranking would need roughly 2.4x more F to change.
  - The M92 verdict against photutils depends on the configuration (tuned: loss, default: win, BIC: tie).
- **Inputs.** The baselines get the exact stamp, which is better than the ISS's three-Gaussian fit. The ISS's extra inputs are the smoothed-count variance map (worth about +0.01 to crowdsource when given) and the local high-pass. Both are part of the method.

**Fix.** Say "ties" instead of "is" for crowdsource, and note that the photutils family frontier is open.

### 10. MINOR (wording and disclosure). Items for the draft

- `exp_astro_m2.py:19-20` still says `DEFAULT` was "chosen on no patch". The first smoke run of the path used exactly (5.0, 1.2) on development patch (110,630) (`astro_smoke2.log`, 18:54). IS1a asked for the rewording. The default is also M2-specific in pixel units.
- The replaced test patches are disclosed in the `exp_astro_m2.py` and `astro_patch_survey.py` docstrings and in spec §5.2, not in `m2.md`. They were never scored (no test key exists before 2d1952d).
- Claim 3 says "from the images and the PSF stamps alone (no truth)". The bright-star table in `psf.md` uses truth.
- Numbers in the claims and docstrings:
  - The M2 residual-sd range is 0.94-1.05 in the docstring (`exp_astro_psf.py:11`) and the claims, but 0.96-1.05 in `psf.md` and in my recomputation (0.96, 1.00 on two patches).
  - M53 is 1.18-1.55, not 1.2-1.5.
  - "5 to 15 % of the peak in single core pixels" is, by my own isotropic fit, 1.7 % on M2, 13.5-15.5 % on M53 and 20-23 % on M92.
- `psf.md` prints "win" verdicts for post hoc rows, and the POST HOC label appears once at the top. Label the rows or drop the verdict column.
- The post hoc adequacy threshold is not tuned in effect: any value between 1.1e-3 and 6.4e-3 gives the same decisions on these stamps. The family (up to 10 parts, offsets up to 2.5 px) was chosen after the failure was seen.

### 11. MINOR (freeze). The guard has holes, and the order evidence is metadata

**Evidence.**
- `frozen_code_hash` does not cover `_test_patch_arrays`, `_run_test_patches`, `truth_arrays`, `patch_path` or the test patch files (`exp_astro_m2.py:296-303`). In a copy, touching `imaging.py` or `MAX_ENTRIES` is refused, but a change inside `_test_patch_arrays` is not.
- The last pre-freeze edit (diff b241ab9..d016d32) was a semantics-preserving refactor of exactly that unhashed function (parallelisation).
- The order d016d32 (2026-11-07 23:39:06) before 2d1952d (2026-11-08 00:35:22) rests on author-set git metadata, supported by scratch logs: the development run ended 23:37:30, `run_astro_T.sh` was created at 23:39:36, T finished 23:57:42, T2 00:17:06, V 00:30:06.
- Unlogged runs cannot be excluded. The only pre-freeze path run I can find is on development patch (110,630).

**Fix.** Hash those functions and the patch files, and anchor the freeze commit externally (a tag pushed before the run).

### 12. MINOR (scope)

- The M2 test patches hold 923-2057 scored sources below 22 mag (0.09-0.21 per px^2). The patch survey lists 257 fully covered candidates up to 14,285 scored sources. (360,1040) is the sparsest of them.
- Development patch (210,630) abuts test patch (310,630). (310,730) abuts (310,630). The strip (110,630), (210,630), (310,630) shares one frame.
- Overall there are 3 frames and 10 windows, so crowding claims are bounded to the light end.
- On M53 and M92, 22 mag is a lower S/N than on M2 (spec §5.2: 4.2-5.3 against 6.3-6.5), so completeness at 22 is not like for like.

### 13. MINOR (reproducibility environment)

- The baselines need astropy 8.0.2, photutils 3.1.1 and crowdsource 0.5.9. Crowdsource is installed from a local clone because the PyPI name `crowdsource` is a different package. None of this is in `tools/requirements-pinned.txt`.
- Default data paths are absolute (`/home/user/astro_data`, `/home/user/astro_data_2`) with environment overrides.
- In the repository venv 5 tests skip because astropy is missing.

## B. Verified (re-derived and agrees)

1. **Chronology and freeze integrity.**
   - d016d32 holds `freeze.json`, the development-only raw (2130 keys, parts D and DB only) and the development md. 2d1952d adds the test keys (T, T2, V; 4619 keys).
   - Every common key is identical across d016d32, 2d1952d and eee967e, NaN-aware, except 4 timing arrays of the variant part; the greedy rows were replaced by greedy_snr.
   - `freeze.json` is unchanged since d016d32. `check_freeze` passes at d016d32, 2d1952d and eee967e.
   - Frozen file hashes, the `frozen_code` hash, the `dev_keys` hash and the truth hashes (M2, M53, M92) match.
2. **Choices made on development patches only.**
   - `run_path` and `imaging.py` read no truth.
   - `choose` and `choose_baselines` read only D and DB keys. From the raw I recomputed (sigma_b 3, r_ex 1.2, kappa 2), dev F 0.561546.
   - I recomputed the baseline selection from the stored catalogues with my scorer: 15 photutils variants (best `ext_all_local_i6`, F 0.2365) and 48 crowdsource filters (six-way tie at 0.290651).
3. **Disjointness.** Test patches have zero pixel overlap with the development patches. (310,630) abuts (210,630); see finding 12.
4. **Test-patch selection.** `astro_patch_survey.py` uses only truth positions and magnitudes. I re-ran it (257 candidates). (360,1040) is the only candidate near 923 scored sources; (100,1100) is one of 10 near 1446.
5. **Second clusters run the M2 configuration verbatim.** All ten test patches in the raw carry `cfg` = (3, 1.2, 2); the tuned baseline names are the same everywhere.
6. **Scoring re-derived.** With my own matching and scorer, the four totals (true22, matched22, det22, det22_matched) agree exactly with the stored `tot` arrays in all 60 rows: 10 patches times the 6 rows ISS tuned, ISS default, crowdsource default and tuned, photutils default and tuned.
7. **Pooled F and bootstrap.** My intervals (own seed, B = 4000, stratified by patch, 20 px blocks):
   - M2: ISS tuned 0.511 [0.488, 0.534], crowdsource 0.251, photutils tuned 0.205, ISS default 0.468. dF against crowdsource is +0.260 [+0.239, +0.282].
   - M53: dF +0.031 [+0.010, +0.051] (tuned), +0.077 [+0.055, +0.100] (default).
   - M92: dF -0.104 [-0.133, -0.078] (tuned), -0.076 [-0.107, -0.048] (default); photutils tuned -0.044, default +0.048.
   - The per-patch verdicts (tuned and default against crowdsource and photutils) agree with `m2.md` for all 10 patches.
   - The BIC pooled F of 0.503 is consistent with the stored totals (arithmetic only).
8. **Label-free diagnostics of claim 3.**
   - With my own mixture fit, the stamp-fit rms is 0.00100, 0.00097, 0.00649, 0.00769, 0.00730 and 0.00836 for M2 (310,630), M2 (100,1100), M53 (660,420), M53 (760,590), M92 (1100,220) and M92 (1430,550).
   - My residual sd is 0.96, 1.00, 1.55, 1.18, 1.84 and 2.28 at the tuned configuration, and 1.39, 1.50, 1.94, 1.43, 2.56 and 5.85 at the default. All match `psf.md`.
   - The 5.85 at the default is mostly depth, not misfit: that catalogue is the capped, high-kappa one.
9. **Kernel study code.**
   - `GeneralPixelPSF`: analytic Jacobian against finite differences 8e-10, windows exact, high-pass partner 2e-16, unit integral 1.000000.
   - On all six M2 stamps `fit_adequate` returns the isotropic fit with identical sigma and parts, so the kernel object and values are bitwise identical.
10. **Reproducibility.**
    - `exp_astro_m2.py --from-raw` and `exp_astro_psf.py --from-raw` (scratch `ASTRO_OUT`, `PSF_OUT`, `PSF_M2_RAW`) regenerate `m2.md` and `psf.md` byte-identically.
    - I re-ran the frozen path on M92 (1100,220): 1330 entries in 216 s, detections bit-identical to the stored 1147 rows.
    - I re-ran crowdsource default on all ten test patches (pooled F 0.251, 0.233, 0.226, exactly as stored) and photutils default on two patches (bit-identical).
    - Tests: 83 passed and 5 skipped (141 s, repository venv); with the domain site-packages on `PYTHONPATH`, `test_astro_second.py` gives 24 passed.
11. **Data.** The sha256 of the six second-cluster `patch_arrays.npz` files and the three truth files match the spec and `freeze.json`. Gain and dark constants are the psField camcol-6 values (4.895 and 0.9025 for M53 and M92; 4.6 and 1.0 for M2).
12. **Photometry constant.** Using the M2 `NMGY` for M53 and M92 shifts magnitudes by -0.024 and -0.002, immaterial against the 0.5 mag tolerance.

## C. Not checked

- Whether the ISS was ever run on a test patch before the freeze without logging. The logs show only development patch (110,630).
- Git timestamps are not independently anchored.
- I did not re-run the 12 development path jobs (hours); the development counts were checked only from the stored arrays.
- I did not re-run the photutils extended variants (17 min each) or the shipped DAOPHOT row.
- I did not explore other photutils options beyond the stored family and IS1a's work.
- I did not test the PSF spatial variation across a patch, or PSF wings beyond the 51x51 stamp.
- I did not quantify the HUGS colour term, saturated-star magnitudes or truth incompleteness. The synthetic worlds only show that such effects of the stated size matter.
- I ran the synthetic experiments on 4 patches: M92 x2, M53 (660,420) and the M2 control. I did not test other M53 patches or the a priori default configuration in the synthetic worlds.
- The draft paper text was not available; I reviewed the claims list only.
