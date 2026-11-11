# Plan, round 2: the astronomy baselines given the ISS's front end, with the options that round 1 left out (photutils `xy_bounds`, crowdsource passes and a 2 px high-pass), selected on the development patches only

Declared 2026-11-08 13:39 UTC, before any round-2 catalogue was made or scored. Round 1 (`baselines_hp_plan.md`, `baseline_catalogs_hp.npz`, `baselines_hp_raw.npz`, `baselines_hp.md`) stays as the first output and is read, never
written, by round 2 (one housekeeping change to the round-1 report text, below). SHA-256 of the round-1 files this round reads: plan `944ad6f8…` (as stored with the round-1 selection), raw `84afb2a37dfba139…`,
catalogues `7cf46e410989f4d9…`; the round-1 report after its one-line format fix: `7cf6e7ae2a392789…`.

## 0. What this round is, and what is known beforehand (declared now)

* Origin: referee round BHP (`review/REVIEW_roundBHP.md`, findings 1, 2, 3, 7). Finding 1: the photutils family had no `xy_bounds`; on the development patches the round-1 selected variant scored F 0.192 and 0.263 with
  `xy_bounds=1`. Finding 2: crowdsource passes options (`miniter 8, maxiter 20, fewstars 20`) and a 2 px high-pass raised the development F by 0.017 to 0.030. The coordinator asked for a declared round 2 with these options in the family.
* **The family below was written after the referee's development probes (`probe_phu.json`, `probe_cs.json`) were read; they are scores on the two M2 development patches.** The development F of the round-2 winners is therefore
  optimistic by selection, in the direction of those probes (photutils about 0.26, crowdsource about 0.44 to 0.45). The probes contain no test-patch score (the referee scored none).
* **The round-1 test results are known to the author of this round** (they were reported in round 1 and in the referee's tables): round 2 is not a blind test of a family fixed before any test patch was seen. What stays as in round 1:
  the selection reads the two development patches only, by a rule fixed here, and the test stage runs the two selected variants without any further choice. The test scores of round 2 are results of a declared procedure, not a pristine
  held-out confirmation.
* No round-2 catalogue has been scored against any truth by the author before this file was written. The witness of the order is the hash of this file, stored in the selection (`SEL|plan_sha256`) and checked by the test stage; the author
  commits nothing (the coordinator may commit this file before the first scoring if an external witness is wanted).
* The ISS results are read, as in round 1, only by the report (`m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`), after the selection is frozen.

## 1. The front end and the scoring: those of round 1, unchanged

High-pass of the frame `H (dn - median)` and of the 51 x 51 stamp (zero-padded, interior operator), the ISS's weights, the match at 0.5 px and 0.5 mag, F-measure at mag < 22, the pooled counts of the two M2 development patches
(110, 630) and (210, 630) for the selection; the round-1 functions (`exp_astro_baselines_hp.highpass_frame`, `highpass_stamp`, `iss_weights`, `cs_stamp`, `load_arrays`, `SETS`) are imported, not copied.

## 2. Family 2, photutils (6 variants)

Each variant is a round-1 variant with `xy_bounds` (photutils `IterativePSFPhotometry(xy_bounds=b)`: the largest distance in px that a fitted position may lie from its initial position, in x and in y), everything else as in round 1
(DAOStarFinder at 5 x the noise-model sigma with the round-1 cuts, `ImagePSF(H stamp)`, fit shape 7, aperture radius 3, no local background, sources kept when `flux_fit` is finite and positive).

| tag | round-1 variant | sigma_b | mode | maxiters | grouper | finder cuts | xy_bounds |
|---|---|---|---|---|---|---|---|
| `phu2_hp3_all_i6_b1` | `phu_hp3_all_i6` (selected in round 1) | 3 | all | 6 | `SourceGrouper(min_separation=4.4)` | relaxed | 1 |
| `phu2_hp3_all_i6_b2` | same | 3 | all | 6 | same | relaxed | 2 |
| `phu2_hp3_new_i30_b1` | `phu_hp3_new_i30` | 3 | new | 30 | none | default | 1 |
| `phu2_hp3_new_i30_b2` | same | 3 | new | 30 | none | default | 2 |
| `phu2_hp5_new_i3_b1` | `phu_hp5_new_i3` (the round-1 reference) | 5 | new | 3 | none | default | 1 |
| `phu2_hp5_new_i3_b2` | same | 5 | new | 3 | none | default | 2 |

What is done about the `all` mode: it stays at **6 iterations** (the round-1 plan abandoned 10 for its cost; the referee reports a bounded `all` i10 run at 7.4 CPU minutes on one development patch, which would be 8+ minutes per
patch and per test patch; it is not in this family, and a deeper bounded `all` run is a limitation). Relaxed finder cuts for the `new` variants (the referee's 0.248 with `xy_bounds=1`, against 0.211 with the default cuts) are not in the
family either: the coordinator's instruction fixes the six variants, and the cuts of each variant are those of its round-1 parent.

## 3. Family 2, crowdsource (8 variants)

The round-1 call (`fit_im` on the high-passed image with the high-passed stamp `SimplePSF(., normalize=0)` in the default's flux convention, the ISS weights, `maxstars 1e7, ntilex = ntiley = 1, refit_psf False`) with three things varied:
the width of the high-pass (2 or 3 px), the detection threshold (3 or 5) and the passes.

| tag | high-pass sigma_b | threshold | passes |
|---|---|---|---|
| `cs2_hp2_t3_pd` | 2 | 3 | d |
| `cs2_hp2_t3_pm` | 2 | 3 | m |
| `cs2_hp2_t5_pd` | 2 | 5 | d |
| `cs2_hp2_t5_pm` | 2 | 5 | m |
| `cs2_hp3_t3_pd` | 3 | 3 | d |
| `cs2_hp3_t3_pm` | 3 | 3 | m |
| `cs2_hp3_t5_pd` | 3 | 5 | d |
| `cs2_hp3_t5_pm` | 3 | 5 | m |

Passes `d`: `miniter 4, maxiter 10, fewstars 60` (the round-1 call: the package's `miniter` and `maxiter`, and the `fewstars 60` that the round-1 call passes; the package's own default is 100). Passes `m`: `miniter 8, maxiter 20, fewstars 20`.
`cs2_hp3_t5_pd` is the round-1 selected variant `cs_hp3_t5` again; it is re-run as a check that the new runner reproduces round 1 bit for bit (it also stands as the baseline of the passes and width effects). No post-filter, as in round 1.

## 4. Selection (the round-1 rule, with the round-1 selection as a candidate)

For each family, among the round-2 variants **and the round-1 selected variant of the family** (`phu_hp3_all_i6`, `cs_hp3_t5`; scored again on the two development patches from the round-1 catalogue file, which must reproduce the round-1 development
totals): the largest pooled F at mag < 22 on the two development patches (counts added, `exp_astro_m2.f_measure`), compared after rounding to 6 decimals, a tie going to the larger tag (`max(tags, key=lambda t: (round(F, 6), t))`; the round-1 tags
`phu_hp...` and `cs_hp...` sort above the round-2 tags `phu2_...` and `cs2_...`, so round 2 returns the round-1 variant if no round-2 variant is better; it cannot return something worse on the development patches). A candidate is eligible only
if it has a score on both development patches. The selection is stored with the SHA-256 of this file (`SEL|plan_sha256`), the hash of the development scores and the round-1 selected tags; the test stage refuses unless the stored selection is
what the rule gives, the round-1 selection verifies (`exp_astro_hp.verify_selection`), and this file is unchanged.

## 5. Test stage

The two selected variants only (no reference variant is re-run), on the 16 test patches of round 1: the four M2 test patches, the four M53, the two M92 and the six repeat patches (scored after the flag file of the ISS repeat run, which exists).
If a selected variant is a round-1 variant its catalogues are copied from the round-1 catalogue file (identical by construction) so that the round-2 files are complete. Nothing is selected on a test patch.

## 6. What the report (`baselines_hp2.md`) shows, and how it is read

* The round-2 development table (every candidate: pooled completeness, purity, F, detections < 22 mag, unmatched detections of all magnitudes, CPU seconds) with the selected variants marked and the old baselines for reference.
* For the groups M2, M53, M92, each repeat frame, and the six repeat patches together, and for every test patch: completeness, purity, F (95 % intervals of the round-1 stratified 20 x 20 px block bootstrap, 2000 resamples, seed 0), purity over all
  magnitudes and the unmatched detections of the round-2 selections beside the round-1 selections and the old tuned baselines and the tuned ISS rows of the frozen results (read only).
* Paired differences: tuned ISS minus round-2 baseline; for M53 and M92 (the post hoc stamp kernel) and the repeat patches (the repeat plan's `adq` row, the kernel the rule chooses from the stamp) also ISS-with-the-kernel-from-the-stamp
  minus round-2 baseline; round-2 minus round-1 selection; round-2 minus old tuned.
* A per-patch verdict count (win / tie / loss of the 95 % interval of the paired F difference on the patch's own 25 blocks) generated by the script, for the 16 patches, for tuned ISS minus each round-2 baseline, with two ISS rows:
  the isotropic-kernel row as frozen (all 16 patches) and the row with the kernel from the stamp (M2: the isotropic row, which is that kernel there; M53 and M92: the post hoc stamp-kernel row of `psf_raw.npz`; the repeat patches: the `adq` row
  of `repeat_raw.npz`, equal to the isotropic row where the rule keeps the isotropic kernel).
* Reading: as in round 1 (win: the ISS is better by the 95 % interval, tie: the interval contains 0, loss), pooled intervals are conditional on the patches, the groups are windows of one frame (M53, M92) or of three frames of the M2 windows (repeat).

## 7. Limits (declared now)

* Options not in the family: relaxed finder cuts for `new` variants; `all` mode beyond 6 iterations; `xy_bounds` other than 1 and 2 px (the referee's 0.5 px probe is not reproduced); the grouper radius and the fit shape in `all` mode; crowdsource `blendthreshu`, `titer_thresh`,
  high-pass widths between 2 and 3 px or below 2; a second pass on the residual image.
* Two development patches and 14 candidates per family (6 + 1, 8 + 1): the winners' development F is optimistic, more than in round 1 because the family was written after dev-patch probes (section 0).
* Everything declared as a limit of round 1 stands (interior stamp operator at the frame border, two windows for M92, windows of one frame, shared truth between M2 and the repeat windows, the uncorrected m2r2 offset).

## 8. Housekeeping on round 1 (item 3 of the instruction; no stored number changes)

`experiments/exp_astro_hp.py` printed `95 %%` literally in the 16 per-patch headings of `baselines_hp.md` and now prints `95 %`; the round-1 report was regenerated from the unchanged raw: exactly those 16 lines differ (`Paired differences (first minus
second, same blocks), 95 %% intervals:` becomes `... 95 % intervals:`), the raw is byte-identical (`84afb2a37dfba139…`). The dropped mutant id F11 of `tools/mutants/wp_astrohp_own.py` and its status (equivalent, survived, removed with the redundant condition it mutated)
are recorded in that file's docstring.
