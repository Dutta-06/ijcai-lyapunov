# Plan (pre-registration): does the factorised small-start Gaussian readout count overlapping nuclei better than tuned classical detectors?

Written 2026-11-10 08:27 UTC (first commit 7974836, 08:27:16 UTC, swept in by another agent's work-in-progress commit; see Amendment 1), before any patch was cut, any method was run on any image or any test image was scored. What had been looked at before this plan: the file listing of the
public Kaggle copy of the 2018 Data Science Bowl stage-1 training set and, per image, its size, whether its three colour channels are identical, its median/99th-percentile/mean grey
level, the number of its masks, the median mask area and the mean grey level inside / outside the masks (to see which images are bright-on-dark fluorescence). Nothing else. Later changes are
listed under *Amendments* with their UTC time. Code: `experiments/exp_blobcount.py`; tests: `tests/test_exp_blobcount.py`. This is a SPIKE: one data set, one patch size, one frozen
configuration per method, the test run once.

## Question
The theory (proved in 1-D with fixed widths; implemented for 2-D images in `src/core/imaging.py`) says that a Gaussian-kernel network with a factorised readout trained by gradient descent
from a small start switches units on one at a time, slides them (centre and width) onto features and leaves the unneeded ones dormant, so the number of units (the model order) is chosen by the
training. Earlier tests found a win or a tie on counting overlapping blobs in simulation and in astronomy, and no advantage elsewhere. This spike asks **one** thing on real data with a
ground-truth count: **in fluorescence microscopy patches with overlapping nuclei, is the count that the sliding inverse scale space (ISS) returns more accurate than that of tuned classical
detectors, and does the order of entry (dual) matter against greedy sliding matching pursuit?**

## Data (fixed before cutting)
* 2018 Data Science Bowl, `stage1_train` (670 images with one binary mask per nucleus), as hosted in the public Kaggle dataset `sinjoysaha/nucleiimagesmasksfromdsb2018` (archive SHA-256
  `912a4fa2449ba1ef7745acb6837f0efc91f83df77494df59612324139625d6ce`, fetched inside a private notebook; no image is committed). The competition itself is not attachable to a notebook through
  the available tools; the archive's file list shows the original layout.
* **Eligible = fluorescence**: the three colour channels identical at every pixel (grayscale) and median grey level < 100 (this excludes the bright-background grey tissue images, the polarity
  is the opposite). Mask-free rule. Grey image = the first channel.
* **Split by IMAGE, seeded rule fixed now**: the eligible ids sorted ascending, permuted by `numpy.random.default_rng(1847).permutation`; the first 60 are DEVELOPMENT, the next 150 are TEST, the
  rest are unused.
* **Patches**: 128 x 128, non-overlapping tiles `[128 r, 128 r + 128) x [128 c, 128 c + 128)` of the image (the remainder at the right / bottom border is dropped). Ground truth of a patch: the
  masks of the whole image with area >= 10 pixels whose centroid lies in the patch (centroid = mean of the mask's pixel coordinates); count = their number. A **base sample** is 2 tiles per image
  (all tiles if the image has fewer), drawn by a per-image seeded permutation; a **crowded boost** adds up to 2 further tiles per image that have >= 25 ground-truth nuclei (the stratum would
  otherwise be too small). The **crowded stratum** is every sampled patch (base or boost) with >= 25 ground-truth nuclei; **all patches** means the base sample only (a representative sample;
  boost patches are not in it unless drawn as base). Per patch the share of nuclei touching another mask is recorded and reported (the second description of crowding).
* Every method sees only the 128 x 128 uint8 patch. Predicted count = number of detections whose centre lies inside the patch.

## Methods (each gets the same tuning budget: at most 24 configurations on the development patches; stop rules that read off the same path are configurations of the same path run)
* **Shared label-free front end** (declared, no ground truth): background = median grey level of the patch; pixel variance `a * smooth(g) + b` with `(a, b)` fitted on the patch from the
  robust relation between the squared horizontal first differences and the local mean (a noise model estimated from the image itself); a patch scale estimate `s_hat` (the scale of the
  strongest scale-normalised Laplacian-of-Gaussian responses of the noise-whitened patch). The classical methods may use `s_hat` (through tuned multipliers) or fixed pixel scales;
  ours and greedy use `s_hat` for the start width of every unit.
* **ours (`dual`)**: the repository's sliding path (`imaging_fast.sliding_path_fast`, the unguarded-equivalent fast path used by the astronomy experiments, guard on), Gaussian pixel-integrated
  kernels with a learned isotropic width per unit (centre and width slide, width inside `[lo, hi] * s_hat`), non-negative amplitudes, the whitened patch as data, an exclusion radius `r_ex`,
  and the stop **BIC** (`stop_bic`, penalty `lam` per parameter, 4 parameters per unit by default) followed by the clean-up of units whose own significance is below `sqrt(lam ln N)` (as
  in the astronomy experiment). The path configuration (start width factor, width bounds, `r_ex`, high-pass width) and the stop (`lam`, clean-up on/off) are chosen on development data.
* **greedy**: the same function with `order='greedy'` (matching pursuit with sliding), the same configuration grid and the same stop. This is the control that matters.
* **LoG, DoG** (`skimage.feature.blob_log`, `blob_dog`) on the same whitened background-subtracted patch: threshold, scale range (as multiples of `s_hat` or in pixels), overlap tuned.
* **peak_local_max** of the Gaussian-smoothed patch: smoothing (multiples of `s_hat` or pixels), min distance, relative threshold tuned.
* **Otsu + distance-transform watershed**: smoothing, Otsu or local-adaptive threshold, h-maxima / min-distance of the distance transform, minimum object area tuned.
* **Optional reference**: a pretrained StarDist / Cellpose, only if installable in the notebook and its weights reachable without work-arounds; reported as a reference, not a target.
* Tuning objective for every method (development only): the mean of the count MAE over the development base sample and over the development crowded stratum. The chosen configuration of every
  method is frozen in `results/ml/blobcount_freeze.json` and committed before the test run.

## Metrics
Count error per patch (detected minus true): MAE, relative error (mean |error| / true count over patches with >= 1 nucleus), bias (mean signed error), per stratum (all / crowded) and per
method; detection F1 of centres against mask centroids at a match radius of 6 pixels (optimal one-to-one assignment, `scipy.optimize.linear_sum_assignment`); paired bootstrap over IMAGES
(10 000 resamples of the test images with replacement, seed 0, percentile 95 % interval, ratio-of-sums for the pooled MAE of a stratum).

## Decision rule (single primary verdict, fixed now)
Let `C*` be the classical method (LoG, DoG, peak_local_max, watershed, each in its frozen configuration) with the lowest TEST MAE in the stratum in question (an oracle choice, favouring the
baselines).
* **C1 (crowded)**: MAE(C*) - MAE(ours) on the crowded stratum is at least 5 % of MAE(C*) AND the 95 % paired-bootstrap interval of the difference excludes 0 (lower end > 0).
* **C2 (greedy)**: MAE(ours) <= MAE(greedy) + 0.5 nuclei per patch on the crowded stratum. (Whether ours beats greedy, interval excluding 0, is reported separately as "the dual order contributes";
  it is not needed for GO.)
* **C3 (non-inferior)**: on all test patches, the upper end of the 95 % paired-bootstrap interval of MAE(ours) - MAE(C*) is at most 10 % of MAE(C*).
* **GO if C1, C2 and C3 all hold; KILL otherwise.** Computed by `verdict()` in `results/ml/blobcount.md`.
* If KILL, the cause is read from the same raw data: (i) ours not better than the best classical on crowded patches (which direction: over- or under-counting, bias); (ii) ours not
  non-inferior overall; (iii) greedy better than ours; (iv) the F1 of centres against the counts; (v) the share of patches that hit the entry cap.

## What will be dropped
Nothing after scoring. Patches on which a method returns no detection stay (count 0). A failed run of a method on a patch counts as 0 detections and is recorded. Optional comparators that
cannot be installed or whose weights are not reachable are not included.

## Amendments
1. 2026-11-10 08:28 UTC: the first draft of the header said 08:30; the file was in fact committed at 08:27:16 UTC as part of commit 7974836 (another agent's work-in-progress commit). Header corrected. No rule changed.
2. 2026-11-10 08:40 UTC decided, committed 08:54 UTC (e75e57e); before any method was run on any TEST patch (smoke runs on 18 development patches, and the first full development run, were made after the amendment's decision and are development only). The first bank built with the rule of the plan (crowded = at least 25 ground-truth nuclei, boost 2 tiles per
   image) contained 4 crowded development patches and 9 crowded test patches from 5 test images, too few for a stratum. Counts only (no image content, no method) were then tabulated over ALL tiles of the
   150 test and 60 development images: tiles with >= 25 nuclei: 4 (development) / 14 (test); with >= 15: 27 / 85 (13 / 26 images). Amendment: the **crowded stratum is defined as at least 15
   ground-truth nuclei** (every sampled patch, base or boost), the boost is up to **4** extra tiles per image with >= 15 nuclei, and the stratum with >= 25 nuclei is reported separately as
   "very crowded" (underpowered, descriptive only). Everything else (split, base sample, margins, rule) unchanged. The test patches' ground-truth counts were therefore seen before
   the freeze (as the plan said they would be, for the stratum), their images were not inspected and no method has touched them. The development patches were used freely to design the grids (smoke runs of every method on 18 of them; see Amendment 3).
3. 2026-11-10 08:55 UTC (before any test patch was touched). Grid design on development patches only (an 18-patch smoke sample, then the full development run): (a) the ISS stop `lam` range of the first draft (4-32 per parameter) was far too small: nuclei are not Gaussian and
   their texture is not pixel noise, so the first-difference noise model underestimates the residual variance (weighted mean square residual 2-20 at the end of the path) and the BIC keeps splitting nuclei; the stop therefore has a **noise calibration**: the BIC
   penalty per parameter is `lam` times the weighted mean square residual of the path's last block (at least 1; label-free), with `lam` in {16, 32, 64, 128, 256}, plus one variant with the clean-up at `sqrt(lam ln N)` (`STOP_GRID`); (b) the scale estimate `s_hat` of the front end
   underestimates the nuclei's LoG scale by about 1.5, so the relative multipliers of the start width, the width bounds and the classical detectors' scale ranges were shifted up (`PATH_GRID`, `SCALES`). A quadratic variance function was tried on the 18 patches and dropped (worse). The budget stays at most 24 configurations per method.
4. 2026-11-10 09:15 UTC (committed fd95e32; development only, no test patch touched). **Development round 1** (the grids of the first code version, stored in `results/ml/blobcount_dev_round1_raw.npz`, report `blobcount_dev_round1.md`): the best configuration of ours and of greedy
   sat at the edge of the grid (the largest penalty with clean-up, the larger start width), and several classical optima sat at an edge too (DoG threshold, peak_local_max smoothing and min distance, LoG overlap, watershed h). One extension round was made **for every method alike**
   (all methods again at exactly 24 configurations, same development patches): ours / greedy: stop penalty `lam` in {64, 256, 1024} with and without the clean-up, start widths {1.3, 1.8} x `s_hat`; LoG / DoG: the scale ranges `abs` and `rel_wide` only, three overlaps {0.15, 0.3, 0.6},
   thresholds LoG {0.04, 0.08, 0.15, 0.3}, DoG {0.15, 0.3, 0.6, 1.2}; peak_local_max: a wider absolute grid; watershed: h-maxima heights {0.5, 1, 2} only for the better thresholds. This is round 2; its results define the frozen configurations; no further round.
   Round-1 development numbers (objective = mean of MAE on base and crowded): ours 4.75, greedy 4.16, LoG 2.67, DoG 2.75, peak_local_max 2.85, watershed 2.44.
5. 2026-11-10 09:30 UTC (development only). Development round 2 done (`blobcount_dev.md`, `blobcount_dev_raw.npz`). Objective of the best configuration: ours 4.54, greedy 4.16, LoG 2.67, DoG 2.73, peak_local_max 2.85, watershed 2.44. The best ours / greedy
   configurations are interior to the start-width axis and at the low end of the extended penalty axis; no third round (the plan allows one extension). Frozen configurations: `results/ml/blobcount_freeze.json` (written by `--stage freeze` from the development raw), committed before the test run.
   The test run is made once, in a private notebook, by the committed code.
6. 2026-11-10 09:55 UTC, AFTER the test result (everything below is POST HOC and is labelled so in the reports). Test run made once (`blobcount.md`): verdict **KILL**. Post hoc additions, none of which changes the verdict: (a) `blobcount_posthoc.md`: the frozen paths of ours and greedy
   with an oracle stop (reads the truth, a diagnostic); (b) `blobcount_ref.md`: pretrained Cellpose (nuclei, default diameter) and StarDist (2D_versatile_fluo, default thresholds) on the same test patches, a REFERENCE only (not tuned, not part of the rule; such models are reportedly trained on data
   that include images of this challenge); (c) `blobcount_nocap.md`: the plan fixed an early stop of the sliding path (a path stops a few entries after an entrant of matched-filter significance below 2, never tuned); the oracle diagnostic showed the cut path ending several units short of the true count on crowded patches,
   so the same grids were run without that stop (`ours_nc`, `greedy_nc`), re-tuned on the development patches by the same objective, frozen (`blobcount_nocap_freeze.json`, committed before its test run) and run once on the test patches.
