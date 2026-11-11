# Plan: the astronomy baselines given the ISS's front end (a high-pass) and a wider pipeline search, selected on the development patches only

Declared 2026-11-08 10:02 UTC, before any catalogue of this study was scored against any truth. **Nothing in this study has been scored against a truth list yet**: no truth file has been read, no
catalogue of the new families has been matched to anything. What has been read: the frozen code (`exp_astro_baselines.py`, `exp_astro_baselines_repeat.py`, `exp_astro_m2.py` (code only), `astro_common.py`, `astro_repeat.py`,
`exp_astro_repeat.py`, `exp_astro_psf.py`, `imaging.py`) and the text files `data_access.md`, `repeat_plan.md` and the baseline passages of `repeat_spec.md` (grep: section 7 and the offsets of section 6.2),
the key names of `baseline_catalogs*.npz` and the stored crowdsource default catalogue of the two development patches (compared array by array with a re-run: no truth), the pixels of the two development patches (`dn`, `psf51`, `gain`, `dark_variance`, `sky_dn`; label-free), the packages' sources. **Not opened before this plan, and not to be opened before the selection is frozen**: `repeat_raw.npz`, `repeat.md`, the `T|`, `T2|`, `P|` keys of `m2_raw.npz`
and `psf_raw.npz`, `m2.md` (the ISS's results; they must not influence what is declared or selected here). The ISS's tuned parameters that this plan uses as a front end are the two sigma_b values the ISS itself
uses (3 px tuned, 5 px a priori default; `repeat_plan.md`, `exp_astro_m2.DEFAULT`).

## 1. The question

A referee (astronomy): the baselines are less tuned than the text says; the 48 crowdsource variants are post-filters on one default `fit_im` run and can only delete detections; no pipeline parameter was
searched while the ISS's tuned parameters (sigma_b, r_ex, kappa) are pipeline parameters; the ISS's dominant ingredient, the border-consistent high-pass, is available to no baseline (no baseline was run on a
high-passed patch with a high-passed PSF stamp). **Does a baseline given the same front end, and a wider search of its own pipeline, close the gap to the sliding ISS?** Either answer is a result.

## 2. The front end given to the baselines (exactly the ISS's, `exp_astro_m2.run_path`)

* Data: `y = H (dn - median(dn))`, `H y = y - S_y y S_x^T`, `S = imaging.smoother(n, sigma_b)` (the border-consistent Gaussian smoother of the ISS); computed by `imaging.PixelPSF(ny, nx, ., hp=sigma_b).highpass`, i.e. by
  the ISS's own method. sigma_b in {3, 5} px (the ISS's two values).
* Weights: `sw = 1 / sqrt(max(gaussian_filter(dn, 1.5, mode='reflect'), percentile(dn, 1)) / gain + dark_variance)` (the ISS's smoothed-count variance); the error array of photutils is `1 / sw`, the weights of
  crowdsource are `sw`.
* Model of a star of unit flux: `H` applied to the 51 x 51 psField stamp. The stamp is embedded in a zero frame of side 51 + 2 ceil(4 sigma_b) (the smoother's truncation radius is 4 sigma_b, so the border
  renormalisation never reaches a stamp pixel), `H` is applied with `smoother(side, sigma_b)` and the result is cropped to 51 x 51 (the star's bowl beyond 25 px, at most 4e-6 of the flux, is dropped). The
  stamp is the only model; the border-consistent operator of the frame is exact for the data and is the interior operator for the model, which differs only for stars within about 2 sigma_b of the frame edge
  (planted-star numbers in section 3).

## 3. Label-free checks done before this plan (synthetic images and the pixels of the development patches; no truth)

Environment: `scratchpad/venv_domain` (python 3.13.18, numpy 2.5.4, scipy 1.18.1, astropy 8.0.2, photutils 3.1.1, crowdsource 0.5.9), one process, `nice -n 10`.

1. **The stamp operator.** On the stamp of development patch (110, 630) the embedded high-pass equals `stamp - scipy.ndimage.gaussian_filter(stamp, sigma_b, truncate=4, mode='constant')` of the zero-padded stamp to
   5e-18 (sigma_b = 3 and 5); applying `smoother(51, sigma_b)` to the bare 51 x 51 stamp instead differs by 4.8e-8 (3 px) and 2.0e-6 (5 px), the stamp's peak being 0.128 (3 px) / 0.137 (5 px). The high-passed stamp
   sums to 3.5e-7 (3 px) / 6.7e-5 (5 px), its minimum is -0.0060 / -0.0033.
2. **Planted stars, photutils** (a Gaussian PSF of sigma 0.93 px as stamp and as truth, an isolated star of 20 000 DN at a random sub-pixel position within 0.5 px of (50, 50), sky 300 DN, noise-free, 12 draws,
   `PSFPhotometry` with the high-passed stamp, `ImagePSF(origin=(25, 25))`, fit shape 7, started at the rounded position): flux / input - 1: mean +0.0046, largest |.| 0.0089 (sigma_b = 3); +0.0051, 0.0069 (5);
   position error: mean -0.0007 / -0.0008 px (x / y), largest |.| 0.0042 px (3); mean -0.0001 / +0.0012, 0.0040 (5). Without any high-pass (raw data, raw stamp) the same fit gives +0.0034 / 0.0064 and
   largest position error 0.0033 to 0.0040 px: the +0.3 to +0.5 % is the spline interpolation of the sampled stamp and is not due to the front end. Stars near the frame border (the data operator is border-consistent,
   the stamp is the interior operator): at 25 and 12 px from the border the same numbers as in the interior; at 6 px +0.3 % / -0.2 % (3 / 5 px) and 0.004 to 0.005 px; at 3 px -1.1 % / -1.7 % and 0.026 / 0.015 px.
3. **Planted stars, crowdsource** (same kind of star, 8 draws, `fit_im` on the high-passed image with the high-passed stamp in `SimplePSF(., normalize=0)`, weights 1/8): flux / input - 1: mean +0.0012, largest |.| 0.0025 (3 px),
   +0.0013, 0.0022 (5 px); position error largest |.| 0.0040 px (3) and 0.0026 px (5); one source per image. The default-style call (raw image, stamp normalised, `normalize=19`) gives +0.0007, 0.0012.
   `SimplePSF`'s own normalisation divides by the sum of the stamp's central 19 x 19 pixels, which is about zero for a zero-sum stamp, hence `normalize=0` and the flux convention of the default is imposed on the
   stamp beforehand (section 5).
4. **The detection threshold of photutils is inert at the noise scale.** On the high-passed patches (DAOStarFinder at n_sig x median(1/sw) = n_sig x 9.1 DN, 30 iterations, 'new' mode) n_sig = 3 and n_sig = 5
   returned bit-identical catalogues on both development patches and both sigma_b (1129 / 1411 rows on (110, 630), 1058 / 1422 on (210, 630)). (For the 5 x 5 DAOFIND kernel of the stamp's FWHM the finder
   scales the threshold by 1.089, so n_sig = 3 and 5 are heights of the best-fitting Gaussian of 29.7 and 49.6 DN, i.e. isolated stars of about m = 22.35 and m = 21.8 (height 0.14 of the flux; m = 22 is
   flux 290 DN, height 41 DN): every local maximum of the convolved crowded field is above both.) The brief's definition of the noise scale, 1.4826 MAD of the high-passed image, is 96 DN (sigma_b = 3) and
   123 DN (5) on (110, 630), 11 to 14 times the pixel noise of the variance model (9.1 DN), because the field is crowded (the MAD measures the confusion of the stars, not the noise). A threshold of 3 or
   5 x that is 315 to 670 DN in the height of the best-fitting Gaussian, an isolated star brighter than m = 19.0 to 19.8, 2 to 3 magnitudes above the benchmark's depth; the first-pass candidate list falls from
   94 to 77 / 56 (3 px) and from 95 to 68 / 52 (5 px). That definition is not used; the noise scale of the threshold is the ISS's variance model, `median(1/sw)`.
5. **The iteration count is the pipeline parameter that sets the depth** (photutils 'new' mode, (110, 630), 3 px): new sources per iteration 94, 86, 63, 53, 47, 41, 37, 37, 35, 35, 33, 32, 31, 31, 31, 31, 31, 31, 31, 30,
   31, 30, 29, 29, 29, 28, 28, 29, 28, 28; the catalogue has 243 rows after 3 iterations, 384 after 6, 528 after 10 and 1129 after 30 (1000 with a finite positive flux). The brief's {3, 6} would stop at 243 / 384
   sources; the default of the existing baseline is 3. A run with maxiters = 10 equals the first 10 iterations of a run with maxiters = 30 exactly (maximum absolute difference 0.0, both development patches and
   both sigma_b), so the 'new' variants are nested.
6. **Run times** (one process, shared machine, user CPU): existing best variant of the old family (`photutils_ext_all_local_i6`, raw data) on (110, 630): 93.6 s, 428 sources; the high-passed counterpart (mode 'all',
   grouper 4.4 px, relaxed cuts, maxiters 6, no local background): 87.0 s, 481 rows (398 with a positive flux); 'new' mode on the high-passed data: 1.8 s (3 iterations), 5 to 8 s (10), 12 to 19 s (30);
   crowdsource: 2 to 4 s per call. The 'all' mode costs far more than the square of the number of iterations (every iteration refits every source, in groups of sources closer than 4.4 px that grow with the
   catalogue): 87 s at 6 iterations, more than 24 CPU minutes at 10 (stopped, see Amendment 1 in section 4); the catalogue file also stores the seconds of each run (`|secs`).
7. **crowdsource with a zero-sum PSF is stable** (no exception, finite fields, one source per planted star) but its diagnostics are not calibrated for it: the median `dflux` is 0.6 DN (3 px) / 6 DN (5 px) against
   41 to 47 DN without the high-pass (the fluxes and positions are right, check 3). The old post-filters (`CS_FILTERS`: flux / dflux, fracflux, qf, spread_model) are therefore not part of the new family;
   a cut at flux >= 5 dflux removed 2 of 564 and 3 of 572 sources of the default catalogue (the old tuned crowdsource is nearly the default catalogue, the old tuned row stays in the report).
8. **crowdsource knobs, counts only** (the package's defaults except the one knob; raw image with the ISS's weights, or high-passed 3 px with them): number of sources for threshold 3 / 5 / 8: 599 / 597 / 571 and
   589 / 584 / 543 (raw, (210, 630) and (110, 630)), 942 / 812 and 886 / 752 (3 px, threshold 3 and 5); blend threshold of the first two passes 0.2 instead of 2: +0.8 % and -2.4 % (raw), +1.0 % and -1.6 % (3 px);
   `maxiter` 20 and `fewstars` 20: +1.8 % and +2.6 % (raw), +0.5 % and 0 (3 px); `psfderiv=False`: -16 % and -9 % sources and no centroid refinement;
   sky grid `nskyx = nskyy = 3`: **the installed package raises `ValueError` in `build_sparse_matrix` (line 500, `colnorm[startidx + syloc[i]] == 0` on an array)**, so the sky model cannot be run
   (a package defect, not changed). `exp_astro_baselines.run_crowdsource` reproduces the stored default catalogues of both development patches exactly (every field, 564 and 572 sources).
   Only the detection threshold changes the catalogue by more than a few per cent; it is the knob searched.

## 4. Family 1: photutils on the high-passed patch (8 variants)

Common: data `H (dn - median)`, errors `1 / sw`, PSF `ImagePSF(H stamp, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)`, `DAOStarFinder(threshold=5 x median(1/sw), fwhm=psf_fwhm(stamp))`
(`exp_astro_baselines.psf_fwhm`, 2.22 px on the development stamps), `IterativePSFPhotometry(fit_shape=(7, 7), aperture_radius=3, n_threads=1)` with photutils' default fitter, **no local background, no xy bounds**,
sources kept when `flux_fit` is finite and positive (as the existing runners). n_sig is fixed at 5 (check 4: 3 gives the identical catalogue).

| tag | sigma_b | mode | maxiters | grouper | finder cuts |
|---|---|---|---|---|---|
| `phu_hp5_new_i3` (**reference**: the settings of the default `photutils` baseline on the a priori high-pass) | 5 | new | 3 | none | default (sharpness 0.2 to 1.0, roundness -1 to 1) |
| `phu_hp3_new_i3` | 3 | new | 3 | none | default |
| `phu_hp3_new_i10` | 3 | new | 10 | none | default |
| `phu_hp5_new_i10` | 5 | new | 10 | none | default |
| `phu_hp3_new_i30` | 3 | new | 30 | none | default |
| `phu_hp5_new_i30` | 5 | new | 30 | none | default |
| `phu_hp3_all_i6` | 3 | all | 6 | `SourceGrouper(min_separation=4.4)` | relaxed (sharpness 0 to 2, roundness -2 to 2) |
| `phu_hp5_all_i6` | 5 | all | 6 | `SourceGrouper(min_separation=4.4)` | relaxed |

The last two are the options of the best variant of the old family (`photutils_ext_all_local_i6`: mode 'all', grouper 4.4 px, relaxed finder cuts, maxiters 6) on the high-passed data without its local background.

**Amendment 1 (2026-11-08 10:35 UTC, before any catalogue of this study was scored; no truth was read).** The 'all' variants were first declared with maxiters 10 (`phu_hp3_all_i10`, `phu_hp5_all_i10`). Their runs did
not finish: the run for (110, 630) at 3 px had used 24 CPU minutes and the one for (210, 630) 15 CPU minutes when both were stopped (by their process ids); the cost of the 'all' mode rises steeply with the
iteration count (87 s at 6 iterations, check 6, every iteration refits every source in groups that grow with the number of sources), and 10 iterations on two patches and two widths cannot be had inside the 3 hours
of CPU of the development search, let alone on 16 test patches. The two variants are replaced by maxiters 6, the depth of the best variant of the old family and the 'all' depth of the brief's {3, 6}. The
catalogues of the other 14 variants had been written (not scored) before the amendment and are unchanged by it.

## 5. Family 2: crowdsource (8 variants)

Common: `crowdsource_base.fit_im(image, SimplePSF, weights, dq=zeros(int64), maxstars=1e7, fewstars=60, ntilex=ntiley=1, refit_psf=False, verbose=False)` with the package's other defaults (`psfderiv=True`,
`nskyx=nskyy=0`, `miniter=4`, `maxiter=10`, `blendthreshu=2`, `titer_thresh=2`); the printed centroid lines are discarded; the output is converted as in `run_crowdsource` (x and y swapped, flux, dflux, qf,
fracflux, spread_model, rchi2 stored; the last four are not used). Only `threshold` is varied (check 8). The catalogue of a variant is the `fit_im` output as returned: **no post-filter**.

| tag | image | stamp | weights | threshold |
|---|---|---|---|---|
| `cs_default` (**reference 1**, the current default) | `dn - sky_dn.mean()` | clipped at 1e-10, unit sum, `SimplePSF(., normalize=19)` | `1 / sqrt(max(dn, 1) / gain + dark)` | 5 |
| `cs_w_t5` | `dn - sky_dn.mean()` | as `cs_default` | `sw` | 5 |
| `cs_w_t3` | same | same | `sw` | 3 |
| `cs_w_t8` | same | same | `sw` | 8 |
| `cs_hp3_t5` | `H (dn - median)`, 3 px | `H` of the `cs_default` stamp divided by the sum of its central 19 x 19 pixels (the default's flux convention), `SimplePSF(., normalize=0)` | `sw` | 5 |
| `cs_hp3_t3` | same | same | `sw` | 3 |
| `cs_hp3_t8` | same | same | `sw` | 8 |
| `cs_hp5_t5` (**reference 2**: the default pipeline on the a priori high-pass) | `H (dn - median)`, 5 px | same | `sw` | 5 |

`cs_default` is `exp_astro_baselines.run_crowdsource` called unchanged. The raw-image rows `cs_w_*` isolate the referee's weights remark and give the raw image the same threshold search as the high-passed rows, so that the
effect of the front end and the effect of the pipeline search are separable. Not searched: sky grid (the installed package fails, check 8), `psfderiv`, blend threshold, `maxiter` / `fewstars` (check 8), the old
post-filters (check 7). The new family therefore contains the old default and, up to 2 or 3 faint sources, the old tuned crowdsource.

## 6. Selection rule (exactly `exp_astro_m2.choose_baselines`, applied to the new families)

On the two M2 development patches (110, 630) and (210, 630), the variant of each family with the largest pooled F at mag < 22 (completeness and purity from the counts of the two patches added,
`exp_astro_m2.f_measure` of the summed `tot`), compared after rounding to 6 decimals, **a tie going to the larger tag** (`max(tags, key=lambda t: (round(F, 6), t))`); a variant is eligible only if it has a score on
both development patches. The selection is written to the raw (`SEL|photutils`, `SEL|crowdsource`) together with the SHA-256 of this file; the test stage refuses to run unless the stored selection equals the
rule applied to the stored development scores and this file is unchanged. Nothing is selected on any test patch. Scoring is `astro_common.score` (match 0.5 px and 0.5 mag, one to one, edge margin 2.5 px,
F-measure at mag < 22). The old baselines are scored on the development patches too (from `baseline_catalogs.npz`, by `exp_astro_m2.baseline_scores`) and listed in the development table for reference; they
take no part in the selection.

## 7. What runs where

* Development stage: all 16 variants on the two development patches -> `results/astro/baseline_catalogs_hp.npz` (keys `'x<x0>_y<y0>|<tag>|x'`, `|y`, `|flux`; crowdsource also `|dflux`, `|qf`, `|fracflux`, `|spread`,
  `|rchi2`; patch pixels, 0-based, flux in DN; the 18 patch origins of M2 (6), M53 (4), M92 (2) and the repeat frames (6) are distinct, checked in the script). Written before it is scored.
* Test stage (only after the selection is frozen): for each of the 16 test patches the **selected variants of both families** and the **references** (`phu_hp5_new_i3`, `cs_default`, `cs_hp5_t5`) and nothing else:
  the four M2 test patches (310, 630), (310, 730), (360, 1040), (100, 1100); the four M53 patches (760, 590), (640, 680), (400, 660), (660, 420); the two M92 patches (1100, 220), (1430, 550); the six repeat
  patches m2r1 (316, 124), (316, 224), m2r2 (133, 281), (195, 689), m2r3 (370, 441), (110, 501) (read with `astro_repeat.load_patch`). The repeat patches are scored only after the file
  `scratchpad/repeat_done.flag` exists (the ISS's repeat run is complete). The selection is not revisited after any test result.
* Report `results/astro/baselines_hp.md` from `baselines_hp_raw.npz` (`--from-raw`, byte-identical): the development table (every variant: pooled completeness, purity, F at mag < 22, detections < 22 mag, unmatched
  detections of all magnitudes, run time), the selected variants; for every test patch and pooled per cluster (M2, M53, M92, the three repeat frames and the six repeat patches together): completeness, purity, F
  (95 % intervals of the 20 x 20 px block bootstrap, 2000 resamples, seed 0, stratified within patches, the same resamples for every row of a table), purity over all magnitudes and unmatched detections, of
  the new baselines beside the old tuned (and default) baselines and the tuned ISS rows of the frozen results (read only from `m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`), and the paired differences
  "tuned ISS minus new baseline" and "new baseline minus old tuned baseline" (verdict by the 95 % interval: win / tie / loss, as `exp_astro_m2.verdict`).

## 8. How the answer is read (declared now)

* **The gap is closed** for a family on a cluster if the interval of pooled dF = F(tuned ISS) - F(selected new baseline) contains 0 or lies below 0; **the ISS still leads** if its lower end is above 0. The
  isotropic-kernel ISS rows of the frozen results are the primary ISS row; on M53, M92 and the repeat frames the ISS with the kernel chosen from the stamp (`psf_raw.npz` P rows, `repeat_raw.npz` `adq` rows) is shown
  as well and the same reading is applied to it.
* The effect of the front end and of the search on the baselines is the paired "new minus old tuned" difference per family, and the reference rows separate the front end alone (`phu_hp5_new_i3` against the old
  default photutils, `cs_hp5_t5` against `cs_default`) from the search (selected against reference).
* Whatever the result, it is reported as it is, with the unmatched detections beside F (F at mag < 22 does not see faint false positives).

## 9. Limits (declared now)

* The 'all' mode, which refits every source jointly, was run to 6 iterations only (about 400 sources of the 1000 or so stars of the benchmark's depth) because its cost grows steeply; a deeper 'all' run was not
  possible here and may do better than the variants declared.
* Two development patches, 16 candidates (8 + 8): the winner's development F is optimistic; the test patches are untouched by the selection. The family is small by design (the brief's cap of 8 per family); a
  larger search could do better, and the search is wider than the old family's only in the depth (photutils), the threshold (crowdsource) and the front end.
* The transfer to M53, M92 and the repeat frames is as for the ISS: nothing is re-tuned, the sigma_b and the pipeline settings come from the M2 development patches. The repeat windows share sky and truth with the
  M2 test windows (`repeat_plan.md`); the offset of m2r2 against the truth (up to 0.26 px in y) is not corrected for any method.
* The photutils model is the interior stamp operator (section 2, check 2); crowdsource's diagnostics are uncalibrated with a zero-sum PSF (check 7); the old post-filters are not searched; the sky grid cannot be
  run (check 8); `psfderiv`, blend threshold and `maxiter` are fixed at the package's defaults.
* The ISS models the PSF by a three-Gaussian fit of the stamp (`fit_psf_parts`), the baselines use the stamp itself; the ISS also has a clean-up by the units' own significance that no baseline has.
* Not a held-out comparison of the ISS against these baselines in the sense of the repeat plan: the ISS rows exist and were produced before this study; they are read, not used.

## 10. Deviations from the brief (to be repeated in the report)

1. n_sig in {3, 5} on the MAD of the high-passed image is replaced by n_sig = 5 on the variance model's pixel noise (check 4: the MAD is the confusion of a crowded field and cuts the catalogue at m of about 19 to 20;
   n_sig = 3 and 5 on the noise model give identical catalogues).
2. maxiters in {3, 6} is replaced by {3, 10, 30} for 'new' and kept at 6 for 'all' (check 5: the depth is set by the iteration count and 3 / 6 stop at 243 / 384 of about 1000 stars; the 'all' mode at 10 iterations
   was declared first and abandoned for its cost, Amendment 1).
3. The crowdsource family varies the detection threshold, not the sky grid (the installed package fails), the blend threshold, `maxiter` or `psfderiv` (check 8); the weights and the high-pass are varied as asked.
4. Post-filters are not part of the new crowdsource family (check 7).
5. A second reference per family (the default pipeline on the a priori high-pass) is run on every test patch beside the chosen variants, and `cs_default` is re-run to reproduce the stored default exactly.
