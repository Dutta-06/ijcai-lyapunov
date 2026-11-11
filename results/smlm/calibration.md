# Calibration of MT0.N1.LD from its development frames 0-999

Produced by `experiments/exp_smlm_calibration.py`. The PSF, offset and emitter statistics use the first 1000 frames (the development frames of the protocol); the camera constants (b, c) and the halo profile use the quiet / single-emitter frames of the stack outside the test frames 1000-1499, with the development-only numbers printed next to them; test frames are never used. Frame indices are 0-based stack indices, positions are in nm with the corner convention (see `experiments/smlm_common.py`). Camera, from the simulation parameters published on the dataset pages (https://srm.epfl.ch/srm/dataset/challenge-3D-simulation/MT0.N1.LD/index.html, read 2026-11-07): baseline 100 DN, EM gain 300 (gamma-distributed), 45 e-/ADU, QE 0.90, read-out noise 74.4 e- (Gaussian), spurious (clock-induced) noise Poisson 0.002, wavelength 660 nm, NA 1.49, pixel 100 nm, "autofluorescence" Poisson at 100 nm pixel size; total gain 6.0 DN per photon (photons below are (DN - 100) / 6, i.e. incident photons, QE included).

Data: `sequence-MT0.N1.LD-2D-Exp-as-stack.zip` (sha256 `f2e478af19799ced...`), `activations.csv` (sha256 `eb631ab5417f3cbc...`), `sequence-MT0.N2.LD-2D-Exp-as-stack.zip` (sha256 `0b0df60ddfdfee97...`), `activations.csv` (sha256 `dffb04292a2f0119...`).

## 0. Conventions of the ground truth, checked against the images

Mean excess of the 3 x 3 pixels around the ground-truth position over the frame median (DN), for 560 bright (> 4000 photons) emitters, as a function of the offset between the stack index and (file frame - 1) and of the axis mapping. The correct mapping stands out sharply.

| mapping | offset -2 | offset -1 | offset +0 | offset +1 | offset +2 |
|---|---|---|---|---|---|
| x -> column, y -> row (used) | 470.0 | 1055.4 | 1999.6 | 1108.2 | 507.1 |
| x <-> y swapped | 64.8 | 129.2 | 211.2 | 129.5 | 65.7 |
| x flipped | 22.9 | 44.4 | 83.3 | 56.1 | 30.0 |
| y flipped | 58.5 | 80.0 | 116.0 | 86.9 | 62.2 |

Offset 0 is the stack index equal to (file frame - 1), i.e. the ground truth counts frames from 1. The pixel-centre convention is decided by the model-free centroid of the stacked PSF in section (a) (offset near zero under the corner convention, 50 nm in each axis under the other).

## (a) Empirical PSF: 2D Gaussian + constant on isolated emitters

Selection: ground-truth emitters of frames 0-999 with no other ground-truth emitter of the same frame within 8 px, the full 11 x 11 px fit window inside the image, and at least 2000 photons: **664 emitters** (of 1253 ground-truth activations in these frames, 860 isolated and interior at any brightness). Least-squares fit in DN with amplitude, centre, width and constant free, **fit window 11 x 11 px unless another window is named**; "point" samples the Gaussian at the pixel centres (as Picasso does), "integrated" integrates it over each pixel. Expected from the optics (0.21 lambda / NA, 660 nm, NA 1.49): sigma = 0.93 px = 93 nm.

| quantity | median | mean | std | IQR | robust sd (1.48 MAD) | 5% .. 95% |
|---|---|---|---|---|---|---|
| sigma, point-sampled circular (px, 11 x 11 window) | 1.398 | 1.554 | 0.517 | 0.229 | 0.146 | 1.249 .. 2.531 |
| sigma, pixel-integrated circular (px, 11 x 11 window) | 1.369 | 1.535 | 0.536 | 0.234 | 0.151 | 1.218 .. 2.543 |
| sigma, geometric mean of the elliptical fit (px, 11 x 11 window) | 1.399 | 1.546 | 0.488 | 0.230 | 0.145 | 1.247 .. 2.504 |
| ellipticity |sx - sy| / mean (11 x 11 window) | 0.051 | 0.080 | 0.102 | 0.072 | 0.047 | 0.007 .. 0.244 |

**Median sigma, all |z| (11 x 11 window): 1.398 px (point-sampled, 140 nm), 1.369 px (pixel-integrated, 137 nm)**, against 0.93 px expected from the optics (n = 664 fits). The width is not a constant: the PSF broadens with defocus, and the emitters sit at z = -633 ... 609 nm, so the spread is partly physical (defocus), partly fit noise. The defocus widths depend on the fit window (table below).

**In focus (|z| < 100 nm, 226 emitters, 11 x 11 window): median sigma 1.277 px (pixel-integrated), 1.310 px (point-sampled)**, i.e. 128 nm, 37 % wider than the 0.93 px of the ideal optics.

Pixel-integrated minus point-sampled: median 0.031 px (a point-sampled model must be wider by about sqrt(sigma^2 + 1/12) - sigma = 0.030 px).

Width by axial position |z| of the emitter (11 x 11 window):

| |z| | n | median sigma_int (px) | median sigma_point (px) | IQR sigma_int |
|---|---|---|---|---|
| 0-100 nm | 226 | 1.277 | 1.310 | 1.249 .. 1.321 |
| 100-200 nm | 253 | 1.393 | 1.423 | 1.319 .. 1.446 |
| 200-350 nm | 103 | 1.656 | 1.681 | 1.504 .. 1.786 |
| 350-800 nm | 82 | 2.358 | 2.344 | 1.977 .. 2.944 |

Window dependence of the widths: the SAME 460 emitters (isolated within 14 px, at least 2000 photons, a 19 x 19 window inside the image) fitted with pixel-integrated Gaussians in windows of different size, median sigma in px (a window that is too small for a broad defocused PSF cuts its wings, a large one admits the halo and the neighbours):

| emitters | n | 7 x 7 | 9 x 9 | 11 x 11 | 13 x 13 | 19 x 19 |
|---|---|---|---|---|---|---|
| all |z| | 460 | 1.338 | 1.343 | 1.368 | 1.392 | 1.429 |
| |z| 0-100 nm | 150 | 1.285 | 1.274 | 1.282 | 1.298 | 1.314 |
| |z| 100-200 nm | 203 | 1.364 | 1.371 | 1.394 | 1.417 | 1.449 |
| |z| 200-350 nm | 74 | 1.513 | 1.535 | 1.584 | 1.632 | 1.717 |
| |z| 350-800 nm | 33 | 1.897 | 1.942 | 2.006 | 2.310 | 2.495 |

Across the windows the median width of the in-focus emitters (|z| < 100 nm, n = 150) varies by 3 %, that of the strongly defocused ones (|z| >= 350 nm, n = 33) by 32 %: a defocus width is only meaningful together with the window it was fitted in, an in-focus width hardly depends on it.


Width by brightness (11 x 11 window; to see that the selection does not drive the value):

| photons | n | median sigma_int (px) | std |
|---|---|---|---|
| 2000-3500 photons | 149 | 1.337 | 0.567 |
| 3500-5000 photons | 116 | 1.373 | 0.517 |
| 5000-inf photons | 399 | 1.381 | 0.530 |

Fitted centre minus ground truth (nm, pixel-integrated fit, 11 x 11 window) under the two pixel conventions: the right one has a small median offset, the wrong one is displaced by a further 50 nm in each axis.

| convention | median dx (nm) | median dy (nm) | std dx (nm) | std dy (nm) |
|---|---|---|---|---|
| corner (centre of pixel j at (j + 0.5) * 100 nm; used) | 10.9 | -8.8 | 22.8 | 20.7 |
| centre (centre of pixel j at j * 100 nm) | -39.1 | -58.8 | 22.8 | 20.7 |

### Offset of a symmetric centre against the simulated emitter: it depends on the estimator and on |z|

The simulated PSF is not symmetric, so the centre that a symmetric estimator finds is displaced against the emitter position, and by how much depends on WHICH estimator and on |z|. It is not a constant of the data set: it must not be read as a convention error, and a correction estimated for one estimator (or one z range) does not carry over to another. Estimator by estimator (offset = estimate minus truth, nm, corner convention; isolated emitters of the development frames, at least 2000 photons, a 15 x 15 window inside the image; the background is the quiet-frame mean of section (b)):

| emitters | n | centroid of the stacked PSF, r < 3 px | r < 5 px | r < 6 px (whole PSF) | Gaussian fit to the stacked PSF, r < 5 px (nm; width in px) | per-emitter centroid (median, 4.5 px) | per-emitter integrated Gaussian fit (11 x 11, median) |
|---|---|---|---|---|---|---|---|
| all |z| | 612 | (+6.7, -6.1) | (+2.3, -3.6) | (+0.3, -2.0) | (+10.3, -8.9), sigma 1.361 | (+4.5, -4.0) | (+10.5, -8.4), n 612 |
| |z| 0-100 nm | 220 | (+5.1, -5.6) | (-1.4, -1.3) | (-4.0, +0.4) | (+7.6, -8.2), sigma 1.303 | (-0.2, -2.3) | (+7.2, -7.6), n 220 |
| |z| 100-200 nm | 244 | (+7.3, -5.4) | (+2.6, -2.5) | (+0.6, -0.9) | (+11.4, -7.7), sigma 1.371 | (+6.0, -4.0) | (+11.3, -7.5), n 244 |
| |z| 200-350 nm | 86 | (+8.6, -6.8) | (+7.9, -7.4) | (+7.5, -6.0) | (+15.2, -11.1), sigma 1.519 | (+13.1, -7.7) | (+14.5, -10.2), n 86 |
| |z| 350-800 nm | 62 | (+10.6, -14.0) | (+10.4, -14.8) | (+7.1, -12.2) | (+30.1, -35.3), sigma 1.879 | (+28.7, -23.6) | (+32.0, -27.7), n 62 |

"Stacked PSF": the pixels of all emitters of the row at their position relative to the ground-truth position, each emitter normalised by its ground-truth photons; its centroid within a radius of the ground-truth position is the offset of the mean profile, model free. Pooled over |z| (row "all |z|") it is (+6.7, -6.1) for r < 3 px, (+2.3, -3.6) for r < 5 px and (+0.3, -2.0) for the whole PSF (r < 6 px): the offset of a centre estimate depends on how much of the PSF it uses, and the whole-PSF centroid is within 2 nm of the ground-truth position in each axis, which is the evidence for the corner convention (the other convention would be 50 nm off in each axis). The offsets of the Gaussian fits are those of the PSF core against the centre of mass; per-emitter estimators scatter around them; all of them change with |z| (rows). Picasso's own least-squares fits show their offsets in `picasso_baseline.md`.

Residual offset of the pixel-integrated fit (11 x 11 window, as above) under the corner convention by axial position:

| z | n | median dx (nm) | median dy (nm) |
|---|---|---|---|
| -800 .. -400 nm | 31 | 55.0 | -59.0 |
| -400 .. -150 nm | 56 | 8.1 | -14.6 |
| -150 .. 150 nm | 386 | 9.1 | -7.5 |
| 150 .. 400 nm | 163 | 14.9 | -8.4 |
| 400 .. 800 nm | 28 | 61.8 | -32.2 |

A method that fits a symmetric model therefore carries an offset of about (11, -9) nm (this estimator, all |z|, median) into its RMSE, and the offset grows with |z|; `score` reports `bias_x_nm`, `bias_y_nm` and `rmse_centered_nm`, and the baseline ladder reports every row also with a constant offset estimated on the development frames subtracted (`picasso_baseline.md`).

## (b) Background and noise on the quiet frames, and the PSF halo

**Quiet frames**: frames with no ground-truth emitter in the frame or in either neighbouring frame (so that no emitter, and no signal of one that has just switched on or off, is in the frame). In MT0.N1.LD there are 1964 such frames among the 19496 frames outside the test frames (102 among the 1000 development frames); in MT0.N2.LD 1964 among the 19496 frames outside the test frames. All values in photons, x = (DN - 100) / 6; standard errors from 10 contiguous groups of frames.

| statistic | MT0.N1.LD, quiet frames of the development frames | MT0.N1.LD, quiet frames outside the test frames | MT0.N2.LD, quiet frames outside the test frames |
|---|---|---|---|
| quiet frames | 102 | 1964 | 1964 |
| mean (photons) | 99.932 +- 0.027 | 99.919 +- 0.006 | 9.915 +- 0.002 |
| mean (DN) | 699.59 | 699.51 | 159.49 |
| median (photons) | 99.33 | 99.33 | 9.33 |
| std of the pixel values (photons) | 14.919 | 14.916 | 4.719 |
| skewness | 0.223 | 0.224 | 0.703 |
| temporal variance, mean over pixels (photon^2) | 222.65 +- 0.55 | 222.49 +- 0.09 | 22.27 +- 0.01 |
| variance / mean (photons) | 2.228 | 2.227 | 2.246 |
| sd of the frame means (photons) | 0.230 | 0.233 | 0.073 |
| ... expected from independent pixel noise alone | 0.233 | 0.233 | 0.074 |
| linear trend of the frame mean (photons per 1000 frames) | 0.026 | -0.000 | -0.000 |
| sd of the per-pixel mean map (photons) | 1.455 | 0.340 | 0.106 |
| ... expected from sampling alone | 1.477 | 0.337 | 0.106 |

The background of a quiet frame is **constant in space and in time**: radial means of the quiet-frame mean image (MT0.N1.LD, frames outside the test frames), photons: r 0-8 px 99.914 +- 0.023; r 8-16 px 99.927 +- 0.014; r 16-24 px 99.925 +- 0.011; r 24-32 px 99.914 +- 0.009; r 32-40 px 99.915 +- 0.012; r 40-46 px 99.922 +- 0.032; the sd of the frame means equals that of independent pixel noise (0.233 against 0.233), so there is no common-mode fluctuation from frame to frame, no spatial profile (the sd of the per-pixel mean map equals its sampling error, 0.340 against 0.337) and no trend. The background of MT0.N1.LD is 99.92 photons (DN 699.5), that of MT0.N2.LD 9.92 photons (DN 159.5).

**What the frames with emitters show instead is the PSF halo of the emitters.** In the 13850 frames outside the test frames that have at least one emitter, the pixels farther than 8 px from every ground-truth emitter of the frame ("emitter-free" pixels) have a mean of 100.54 photons, 0.62 above the quiet-frame background, and the sd of the per-frame mean of these pixels is 0.530 (quiet frames: 0.233). Regressing that frame mean on the total ground-truth photons of the frame gives a slope of 7.62e-05 photons per pixel per emitted photon (correlation 0.860); the residual sd after removing it is 0.270 (pixel noise alone: 0.242). Over the 3784 emitter-free pixels per frame this is 0.29 photons beyond 8 px per emitted photon: signal far from the emitter, not illumination (the background of quiet frames is flat) and not a common-mode fluctuation (it is a deterministic function of the photons emitted in the frame).

### Radial excess profile around isolated single-emitter frames

Frames outside the test frames with exactly one ground-truth emitter, brighter than the lower quartile of the whole set's activations (MT0.N1.LD: 1851 photons; MT0.N2.LD: 200), at least 16 px from every border; the image minus the quiet-frame background of its set is binned by the distance of the pixel centre to the ground-truth position. Cumulative fraction of the emitter's ground-truth photons within a radius = sum of the mean excess per pixel of the annuli times their full areas (isotropy assumed beyond the frame edge) divided by the mean ground-truth photons; +- from a leave-one-group-out jackknife over 10 groups of emitters.

| set | emitters | n | mean ground-truth photons | within 3 px | within 5 px | within 10 px | within 25 px | whole image (64 x 64), summed directly |
|---|---|---|---|---|---|---|---|---|
| MT0.N1.LD | all | 2882 | 6389 | 0.740 +- 0.002 | 0.945 +- 0.001 | 1.159 +- 0.002 | 1.424 +- 0.003 | 1.413 |
| MT0.N1.LD | in focus (|z| < 100 nm) | 969 | 6480 | 0.838 +- 0.001 | 1.006 +- 0.001 | 1.189 +- 0.001 | 1.448 +- 0.004 | 1.435 |
| MT0.N1.LD | defocused (300 <= |z| < 700 nm) | 268 | 6267 | 0.510 +- 0.006 | 0.807 +- 0.006 | 1.134 +- 0.006 | 1.448 +- 0.007 | 1.417 |
| MT0.N2.LD | all | 2882 | 691 | 0.741 +- 0.002 | 0.945 +- 0.002 | 1.161 +- 0.001 | 1.436 +- 0.005 | 1.434 |
| MT0.N2.LD | in focus (|z| < 100 nm) | 969 | 700 | 0.842 +- 0.002 | 1.010 +- 0.004 | 1.194 +- 0.005 | 1.465 +- 0.009 | 1.475 |
| MT0.N2.LD | defocused (300 <= |z| < 700 nm) | 268 | 677 | 0.507 +- 0.006 | 0.794 +- 0.007 | 1.120 +- 0.011 | 1.457 +- 0.034 | 1.444 |

Mean excess per pixel per 1000 ground-truth photons in each annulus and the cumulative fraction (all emitters of the row "all" above):

| annulus (px) | MT0.N1.LD excess / px / 1000 photons | MT0.N1.LD cumulative | MT0.N2.LD excess / px / 1000 photons | MT0.N2.LD cumulative |
|---|---|---|---|---|
| 0-1 | 56.6827 | 0.178 | 56.8732 | 0.179 |
| 1-2 | 35.9047 | 0.516 | 35.7958 | 0.516 |
| 2-3 | 14.2549 | 0.740 | 14.3185 | 0.741 |
| 3-4 | 5.1859 | 0.854 | 5.1533 | 0.854 |
| 4-5 | 3.2093 | 0.945 | 3.2011 | 0.945 |
| 5-6 | 1.7869 | 1.007 | 1.7688 | 1.006 |
| 6-7 | 1.0772 | 1.051 | 1.0945 | 1.051 |
| 7-8 | 0.9101 | 1.094 | 0.9145 | 1.094 |
| 8-9 | 0.6253 | 1.127 | 0.6438 | 1.128 |
| 9-10 | 0.5345 | 1.159 | 0.5549 | 1.161 |
| 10-11 | 0.4470 | 1.189 | 0.4434 | 1.190 |
| 11-12 | 0.3469 | 1.214 | 0.3647 | 1.217 |
| 12-14 | 0.2845 | 1.260 | 0.2967 | 1.265 |
| 14-16 | 0.2209 | 1.302 | 0.2309 | 1.309 |
| 16-18 | 0.1808 | 1.340 | 0.1894 | 1.349 |
| 18-20 | 0.1520 | 1.377 | 0.1445 | 1.384 |
| 20-22 | 0.1159 | 1.407 | 0.1336 | 1.419 |
| 22-23 | 0.0683 | 1.417 | 0.0736 | 1.429 |
| 23-24 | 0.0341 | 1.422 | 0.0395 | 1.435 |
| 24-25 | 0.0109 | 1.424 | 0.0051 | 1.436 |
| 25-26 | 0.0050 | 1.424 | 0.0372 | 1.442 |
| 26-28 | -0.0041 | 1.423 | -0.0018 | 1.441 |
| 28-30 | -0.0012 | 1.423 | -0.0102 | 1.438 |
| 30-32 | -0.0039 | 1.421 | 0.0096 | 1.441 |
| 32-36 | 0.0004 | 1.421 | 0.0047 | 1.445 |
| 36-40 | 0.0026 | 1.424 | 0.0066 | 1.452 |
| 40-50 | -0.0022 | 1.418 | -0.0047 | 1.439 |

The ground-truth "intensity" corresponds to the core of the PSF: the cumulative excess within 5 px is 1.01 of it for the emitters in focus (0.95 over all emitters, whose PSFs are broader). But the PSF has a long halo: for the in-focus emitters the cumulative fraction keeps rising to 1.45 at 25 px (44 % of the ground-truth photons arrive beyond 5 px), and then the profile is cut off abruptly: the excess per pixel per 1000 ground-truth photons (all emitters) is 0.116 at 20-22 px, 0.068 at 22-23 px, 0.011 at 24-25 px and -0.004 (consistent with zero) at 26-28 px. In MT0.N2.LD the fractions are 0.74 at 3 px, 0.94 at 5 px, 1.16 at 10 px, 1.44 at 25 px; they differ from those of MT0.N1.LD by at most 0.01, so normalised by the photons the profile is the same in both noise classes: it is a property of the simulated PSF, not of the noise. The cause of the halo was not investigated. Consequences: a "Gaussian + constant" patch model (Picasso's fit) cannot hold the halo, so its photon count is biased low (`picasso_baseline.md`), and any estimate of the background from the neighbourhood of an emitter (as the earlier version of this file did with an 8 px mask: 100.3 photons instead of 99.92) is raised by it.

## (c) Noise model: var = a x + b in photon units, x = (DN - 100) / 6

Nominal EMCCD model (gamma gain noise, excess noise factor^2 = 2, read noise 74.4 e-, quantisation): a = 2 / QE = 2.222, b = (read noise / (QE EM))^2 + 1/(12 G^2) = 0.0782 photon^2.

**One operating point (quiet frames of MT0.N1.LD).** Mean x = 99.919 photons, temporal variance 222.49 +- 0.09 photon^2 (mean over the pixels of the per-pixel variance), variance / mean = 2.227, a_eff = (var - b) / x = 2.226 with the nominal b, against a = 2.222 expected.

**Two operating points (the backgrounds of two noise classes of the same camera).** MT0.N1.LD: x = 99.919 photons, var = 222.49 +- 0.09 photon^2, var / mean = 2.227; MT0.N2.LD: x = 9.915 photons, var = 22.267 +- 0.011 photon^2, var / mean = 2.246. The line through the two points: **a = 2.225 +- 0.001, b = 0.21 +- 0.02 photon^2** (standard errors from 10 paired groups of frames), against the nominal a = 2.222, b = 0.078 photon^2. The fitted b differs from the nominal b by 0.13 photon^2 (7.7 standard errors of the group estimate; a systematic cause such as the spurious noise or an excess-noise factor that is not exactly 2 is not excluded). It is 0.09 % of the variance at the background of MT0.N1.LD (222 photon^2), so it is irrelevant at that level and matters only for very dark pixels (var = a x + b at x = 0).

**Emitter pixels (a rough check).** The 80344 pixels of the 664 isolated-emitter patches (11 x 11 window), residual of the pixel-integrated Gaussian fit divided by (1 - leverage), squared, averaged in bins of the fitted mean x. Weighted regression over the 7 bins (x = 112 .. 781): **a = 2.61 +- 0.20, b = -16 +- 25 photon^2**; with b fixed at the nominal value a = 2.49. The bins are biased by the PSF misfit (upwards: the simulated PSF is not Gaussian, it has the halo of section (b)) and by the imperfect leverage correction of the unweighted fit (downwards at the brightest pixels), so treat them as a consistency check (a of the order of 2-3), not as a measurement.

| mean x (photons) | pixels | variance (photon^2) | nominal a x + b |
|---|---|---|---|
| 111.8 | 39872 | 267.7 | 248.5 |
| 130.0 | 22089 | 344.1 | 289.1 |
| 170.5 | 8193 | 437.2 | 378.9 |
| 240.9 | 5603 | 608.8 | 535.3 |
| 363.8 | 2889 | 884.1 | 808.6 |
| 543.2 | 1405 | 1304.4 | 1207.1 |
| 781.0 | 244 | 1737.8 | 1735.7 |

The regression of per-pixel variance on per-pixel mean across emitter-free pixels is not informative (the background is flat, so the means span only a few tenths of a photon and the sample mean and variance of a skewed distribution are correlated), which is why the two backgrounds, not a regression, fix a and b.

**Result: a = 2.23 (quiet frames; nominal 2.22); two-point a = 2.225, b = 0.21 +- 0.02 photon^2 (nominal 0.078).**


## (d) Ground-truth emitters per frame and photons per emitter (frames 0-999)

1253 ground-truth activations in 1000 frames: mean 1.253 per frame, variance 1.202 (variance / mean = 0.96; 1 for a Poisson count), 28.0 % of the frames are empty, maximum 6.

| emitters in frame | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| frames | 280 | 360 | 226 | 106 | 20 | 5 | 3 |
| Poisson, same mean | 285.6 | 357.9 | 224.2 | 93.7 | 29.3 | 7.4 | 1.5 |

Photons per activation (the ground-truth `intensity`, photons emitted during the frame; partial-frame activations are included; development frames): mean 4693, median 4259, min 81, max 10573. Whole set (24938 activations): median 4081, 25 % quantile 1851.

| 1% | 5% | 25% | 50% | 75% | 95% | 99% |
|---|---|---|---|---|---|---|
| 140 | 386 | 2031 | 4259 | 7180 | 10114 | 10406 |

Fraction of the development activations below 250 / 500 / 1000 / 2000 photons: 3.4 % / 6.8 % / 12.9 % / 24.7 %. A quarter of the activations of the WHOLE set are dimmer than 1851 photons (the quartile used throughout the harness; the development frames alone give 2031); the challenge discarded the dimmest quarter, after pairing, when it assessed software (rule from memory, `smlm_common.py`). A 4081-photon emitter (the median of the whole set) of width 1.37 px has a peak signal of about 347 photons (2081 DN) per pixel above the background of 99.92 photons.

**The illumination of the simulation acts on the emitters, not on the background.** The simulator's summary file (`stats_activation_sumary.csv`) lists an illumination profile (centre of power (2948.5, 3473.7) nm, both coordinates are labelled "X" in the file; distance at half maximum 3520 nm). The brightness of the activations of the whole set falls with the distance of the emitter to that centre (mean photons per activation in the table), while the quiet frames of section (b) are flat:

| distance to the centre | activations | mean photons | median photons |
|---|---|---|---|
| 0-1000 nm | 3848 | 5485 | 5023 |
| 1000-2000 nm | 9397 | 5098 | 4611 |
| 2000-3000 nm | 7939 | 4262 | 3842 |
| 3000-4000 nm | 3754 | 3194 | 2907 |

Axial positions: z from -633 to 609 nm (median |z| 143 nm); 74 % of the activations lie within +-250 nm of the focal plane.

**Persistence of an emitter over consecutive frames.** A fluorophore stays on for several frames (runs of the same coordinates in consecutive frames of the whole ground truth), so neighbouring frames are not independent samples and frame statistics need a block bootstrap (`smlm_common.bootstrap_ci`, blocks of 10 frames):

| set | runs | mean run length (frames) | runs of 3 or more frames | runs of 1 | runs of 2 | runs of 3 | runs of 4 | runs of 5 | runs of 6 | runs of 7 | runs of 8 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MT0.N1.LD | 10989 | 2.269 | 31.0 % | 3582 | 4001 | 1776 | 854 | 392 | 189 | 106 | 51 |
| MT0.N1.HD | 13981 | 2.261 | 31.1 % | 4669 | 4958 | 2283 | 1059 | 527 | 255 | 125 | 61 |
| MT0.N2.LD | 10989 | 2.269 | 31.0 % | 3582 | 4001 | 1776 | 854 | 392 | 189 | 106 | 51 |
| MT0.N2.HD | 13981 | 2.261 | 31.1 % | 4669 | 4958 | 2283 | 1059 | 527 | 255 | 125 | 61 |

