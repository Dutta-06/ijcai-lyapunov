# Picasso baseline ladder on the simulated 2D microtubule sets

Produced by `experiments/exp_smlm_picasso.py` with the harness `experiments/smlm_common.py`. Frame indices are 0-based stack indices, positions in nm (corner convention; Picasso's pixel-index coordinates are converted with nm = (p + 0.5) * 100), matching tolerance 250 nm (Hungarian pairing per frame, inclusive), metrics pooled over frames. Every number below is printed by that script; the file regenerates byte for byte (no timings).

**The baseline is a ladder, not one setting.** A baseline is only as fair as the options it may use (referee finding SM1-1, `review/REVIEW_roundSM1.md`): tuning the minimum net gradient alone leaves Picasso well below what its own documented options reach, and a single operating point is not a baseline (the point of best efficiency and the point of best Jaccard index differ strongly). Rungs, each tuned on development frames only and scored once on held-out test frames:

* **R0**: Picasso's defaults (minimum net gradient 5000, box 7 px, raw frames, no padding, least-squares fit `gausslq`).
* **R1E / R1J**: the minimum net gradient alone tuned (20 values) for efficiency E / Jaccard J (the earlier "tuned (efficiency)" / "tuned (Jaccard)" rows).
* **R2E / R2J**: every documented option tuned together for E / J: identification box, fit box, Picasso's Gaussian filter, median padding of each frame (Picasso cannot identify in a border zone of box // 2 px; the padded frame removes it, positions are shifted back and spots whose centre lies outside the field are dropped), fit method (`gausslq` / `gaussmle`), and post-filters (net gradient, minimum photons, width window |sx|, |sy|). The grids are declared at the top of the script before any run; the search is described in each set's section and the number of configurations tried is reported there, so that another method can be given the same tuning effort.
* **+c** (appended to any row): the same detections with a constant position offset subtracted, estimated on the development frames only as the mean signed error (detection - truth) of the matched true positives of that configuration (the offset of a symmetric fit against the simulated emitter, `calibration.md`; it depends on the estimator and on |z|).

Two scoring settings on the same detections: ALL ground-truth emitters (primary), and the challenge's assessment setting (emitters closer than 450 nm to the border and the dimmest quarter of the whole set's activations dropped after pairing; the rules are from memory except the 250 nm tolerance and the metric definitions, see the docstring of `smlm_common.py`). Intervals are 95 % percentile intervals of a bootstrap over blocks of 10 consecutive frames (an emitter stays on for several consecutive frames, see the run lengths in each set's section and `calibration.md`); the earlier i.i.d. frame bootstrap was too narrow. They cover the sampling of frames only, not other datasets or the simulation, and the development and test frames share the same three microtubules. The tuning effort on the development frames is part of the result: a method that is tuned less is compared against a better-tuned baseline.

**Exact Picasso setting.** Picasso 0.11.3 (`picasso.localize.identify` and `fit`, in-memory stack, single thread; `localize_frames` is the same two calls with one box). Identification: local maxima over a box of 7 px (default) whose net gradient exceeds the minimum net gradient (default 5000, in raw DN, on the frame or, with the Gaussian filter, on the smoothed frame); fit: `gausslq`, Picasso's default least-squares fit (Levenberg-Marquardt, elliptical Gaussian with free sx, sy, amplitude, background, centre; convergence schedule of Picasso: tolerance 0.01, at most 200 iterations), box 7 x 7. Camera passed to Picasso: baseline 100 DN, EM gain 300, sensitivity 45 e-/ADU, QE 0.9 (Picasso converts to photoelectrons (DN - baseline) x sensitivity / gain and does not divide by QE; the harness divides the reported photons by QE; the photon count of its Gaussian fit is biased low against the emitted photons because the Gaussian core does not hold the PSF halo, `calibration.md`; table in each set section). **The Gaussian of Picasso's least-squares fit is evaluated at the pixel centres, not integrated over the pixel** (the integrated-Gaussian fit exists only in its deprecated maximum-likelihood module); with the fitted sigma of about 1.3-1.4 px the pixel-integration correction is a width change of about 0.03 px. R0 and R1 have no post-filtering of the localizations (no sigma, photon or precision filter, no drift correction, no linking of consecutive frames). Cross-checks done once by hand, not by this script: Picasso's own command line (`picasso localize -a lq -b 7 -g 5000 -bl 100 -s 45 -ga 300 -qe 0.9 -px 100`, whose default fit would otherwise be `mle`) on a TIFF of test frames 1000-1199 gave exactly the localizations of this runner (158 spots, largest position difference 2e-4 nm); and the independent pipeline of the SM1 referee (Picasso identification and fit called separately, as here) gave the same R0 numbers and, for the configuration "identification box 7, fit box 9, Gaussian filter 1.5, median pad 5, net gradient > 2000" on the MT0.N1.LD test frames, the same J, precision, recall, RMSE and E as this runner (82.5, 99.4, 82.9, 45.7, 51.1).

Data: MT0.N1.HD/activations.csv sha256 cd335c3029c6...; MT0.N1.HD/sequence-MT0.N1.HD-2D-Exp-as-stack.zip sha256 0b97f4ce2820...; MT0.N1.LD/activations.csv sha256 eb631ab5417f...; MT0.N1.LD/sequence-MT0.N1.LD-2D-Exp-as-stack.zip sha256 f2e478af1979...; MT0.N2.HD/activations.csv sha256 035bd15a12d8...; MT0.N2.HD/sequence-MT0.N2.HD-2D-Exp-as-stack.zip sha256 bebae60f119e...; MT0.N2.LD/activations.csv sha256 dffb04292a2f...; MT0.N2.LD/sequence-MT0.N2.LD-2D-Exp-as-stack.zip sha256 0b0df60ddfdf.... Python 3.13.18, numpy 2.5.4.

## Summary: test frames, every row of the ladder

Point estimate [95 % block-bootstrap interval] on the held-out test frames; the sections of the sets below give every metric, the paired differences, the development-frame search and the ablations.

Jaccard %, all ground-truth emitters (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 71.7 [66.8, 76.3] | 47.8 [46.3, 49.3] | 22.1 [18.4, 25.9] | 18.1 [17.2, 19.1] |
| R0+c | 71.7 [66.8, 76.3] | 47.8 [46.3, 49.3] | 22.1 [18.4, 25.9] | 18.1 [17.2, 19.1] |
| R1E | 64.8 [60.0, 69.6] | 46.9 [45.5, 48.4] | 53.1 [48.1, 58.1] | 35.1 [34.0, 36.5] |
| R1E+c | 64.8 [60.0, 69.6] | 46.9 [45.5, 48.5] | 53.1 [48.1, 58.1] | 35.1 [34.0, 36.5] |
| R1J | 73.4 [68.1, 78.2] | 49.0 [47.6, 50.5] | 55.9 [51.2, 60.5] | 42.9 [41.7, 44.3] |
| R1J+c | 73.4 [68.1, 78.2] | 49.0 [47.6, 50.5] | 55.9 [51.2, 60.5] | 42.9 [41.7, 44.3] |
| R2E | 68.7 [63.8, 73.8] | 51.5 [49.9, 53.0] | 57.3 [52.9, 61.7] | 38.7 [37.4, 39.9] |
| R2E+c | 68.7 [63.8, 73.8] | 51.5 [49.9, 53.0] | 57.3 [52.9, 61.7] | 38.7 [37.4, 39.9] |
| R2J | 82.9 [79.5, 86.3] | 56.4 [54.8, 58.0] | 65.0 [61.6, 68.4] | 48.4 [47.1, 49.8] |
| R2J+c | 83.3 [79.8, 86.6] | 56.5 [55.0, 58.0] | 65.0 [61.6, 68.4] | 48.5 [47.2, 49.9] |

RMSE nm, all ground-truth emitters (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 34.7 [27.8, 41.3] | 49.6 [47.5, 51.8] | 23.4 [20.1, 26.8] | 45.7 [43.6, 48.1] |
| R0+c | 29.0 [22.3, 35.3] | 46.4 [44.3, 48.7] | 19.7 [16.6, 23.3] | 44.5 [42.4, 46.8] |
| R1E | 25.3 [21.4, 29.3] | 48.5 [46.5, 50.7] | 33.5 [30.3, 36.3] | 49.4 [47.4, 51.4] |
| R1E+c | 19.9 [16.0, 23.8] | 45.5 [43.4, 47.7] | 30.7 [28.0, 33.3] | 47.2 [45.2, 49.3] |
| R1J | 41.3 [33.0, 49.0] | 53.4 [51.1, 55.8] | 40.2 [35.5, 45.7] | 60.1 [58.1, 61.9] |
| R1J+c | 36.1 [28.0, 43.8] | 50.1 [47.7, 52.7] | 38.4 [33.1, 44.5] | 58.2 [56.1, 60.1] |
| R2E | 24.2 [20.8, 28.3] | 43.5 [41.7, 45.5] | 33.5 [30.0, 37.5] | 45.7 [44.0, 47.3] |
| R2E+c | 19.1 [15.2, 23.6] | 40.6 [38.7, 42.5] | 31.5 [28.3, 35.3] | 43.9 [42.2, 45.6] |
| R2J | 44.3 [37.3, 51.1] | 53.1 [50.8, 55.5] | 52.2 [44.7, 59.8] | 60.5 [58.4, 62.6] |
| R2J+c | 41.0 [33.7, 47.4] | 50.1 [47.8, 52.4] | 50.3 [42.5, 57.7] | 59.1 [56.9, 61.2] |

efficiency E, all ground-truth emitters (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 55.2 [48.5, 61.3] | 28.0 [25.9, 30.1] | 18.7 [14.9, 22.5] | 6.2 [4.9, 7.5] |
| R0+c | 59.5 [53.3, 65.2] | 30.2 [28.1, 32.3] | 19.7 [16.0, 23.5] | 6.8 [5.6, 8.0] |
| R1E | 56.6 [51.6, 61.4] | 28.1 [25.9, 30.2] | 42.4 [38.0, 46.8] | 18.5 [17.0, 20.1] |
| R1E+c | 59.6 [54.6, 64.4] | 30.1 [28.0, 32.3] | 44.0 [39.6, 48.4] | 19.8 [18.2, 21.3] |
| R1J | 50.8 [43.7, 58.0] | 26.2 [24.0, 28.2] | 40.4 [34.6, 45.9] | 17.1 [15.5, 18.8] |
| R1J+c | 55.1 [48.0, 62.1] | 28.5 [26.4, 30.6] | 41.6 [35.4, 47.5] | 18.5 [16.9, 20.2] |
| R2E | 60.4 [55.2, 65.2] | 34.8 [33.0, 36.7] | 45.7 [41.0, 50.4] | 23.6 [22.0, 25.2] |
| R2E+c | 63.3 [58.0, 68.5] | 36.7 [35.0, 38.6] | 46.9 [42.2, 51.5] | 24.6 [23.0, 26.3] |
| R2J | 52.5 [46.1, 59.1] | 31.3 [29.2, 33.5] | 37.1 [30.0, 44.5] | 20.5 [18.7, 22.4] |
| R2J+c | 55.7 [49.6, 62.4] | 33.7 [31.6, 35.9] | 38.7 [31.6, 46.2] | 21.7 [19.6, 23.6] |

Jaccard %, challenge setting (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 91.8 [88.0, 95.6] | 65.0 [63.1, 66.8] | 33.2 [28.2, 38.1] | 27.0 [25.8, 28.2] |
| R0+c | 92.1 [88.4, 95.7] | 65.0 [63.1, 66.9] | 33.2 [28.2, 38.1] | 27.0 [25.7, 28.2] |
| R1E | 87.1 [82.1, 91.6] | 64.5 [62.6, 66.3] | 78.4 [73.1, 83.5] | 51.7 [50.1, 53.4] |
| R1E+c | 87.3 [82.8, 91.6] | 64.6 [62.6, 66.4] | 78.4 [73.1, 83.5] | 51.6 [50.0, 53.3] |
| R1J | 91.9 [87.9, 95.4] | 64.7 [62.8, 66.4] | 80.3 [75.3, 85.2] | 59.7 [58.1, 61.2] |
| R1J+c | 92.2 [88.4, 95.5] | 64.7 [62.8, 66.5] | 80.3 [75.3, 85.2] | 59.6 [58.0, 61.2] |
| R2E | 83.9 [78.0, 89.5] | 65.8 [63.9, 67.7] | 77.6 [72.5, 82.7] | 52.8 [51.1, 54.4] |
| R2E+c | 83.9 [78.0, 89.5] | 65.9 [64.0, 67.8] | 77.6 [72.5, 82.7] | 52.9 [51.2, 54.6] |
| R2J | 95.8 [94.1, 97.8] | 68.3 [66.4, 70.1] | 83.8 [79.4, 88.3] | 61.1 [59.6, 62.8] |
| R2J+c | 95.8 [94.1, 97.8] | 68.4 [66.5, 70.2] | 84.0 [79.7, 88.5] | 61.2 [59.6, 62.9] |

RMSE nm, challenge setting (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 32.2 [25.4, 38.8] | 49.0 [46.8, 51.3] | 23.7 [20.3, 27.3] | 45.7 [43.5, 48.0] |
| R0+c | 26.3 [19.6, 32.6] | 45.9 [43.6, 48.2] | 20.0 [16.8, 23.6] | 44.4 [42.3, 46.8] |
| R1E | 24.5 [20.4, 28.9] | 48.1 [46.0, 50.4] | 32.6 [29.7, 35.4] | 49.0 [47.1, 51.2] |
| R1E+c | 19.0 [14.6, 23.3] | 45.2 [43.1, 47.4] | 29.6 [27.1, 32.3] | 47.0 [45.0, 49.1] |
| R1J | 38.4 [28.5, 46.9] | 51.3 [48.7, 53.8] | 38.4 [33.4, 44.5] | 57.4 [55.2, 59.6] |
| R1J+c | 32.9 [22.9, 41.4] | 48.0 [45.3, 50.5] | 36.7 [31.0, 43.6] | 55.5 [53.1, 57.8] |
| R2E | 21.0 [18.3, 24.8] | 42.6 [40.4, 44.8] | 29.5 [27.1, 32.4] | 45.0 [43.2, 46.8] |
| R2E+c | 14.6 [11.6, 18.7] | 39.5 [37.3, 41.7] | 27.4 [25.0, 30.3] | 43.2 [41.3, 45.0] |
| R2J | 37.4 [29.2, 44.2] | 49.3 [46.8, 51.7] | 48.6 [40.5, 56.3] | 57.1 [54.8, 59.3] |
| R2J+c | 31.7 [23.7, 38.4] | 45.9 [43.5, 48.2] | 46.7 [38.1, 54.8] | 55.4 [53.0, 57.6] |

efficiency E, challenge setting (test frames):

| row | MT0.N1.LD | MT0.N1.HD | MT0.N2.LD | MT0.N2.HD |
|---|---|---|---|---|
| R0 | 66.8 [59.8, 73.6] | 39.7 [37.3, 42.2] | 29.2 [23.9, 34.2] | 13.9 [12.5, 15.3] |
| R0+c | 72.6 [66.1, 79.3] | 42.3 [39.8, 44.8] | 30.3 [25.2, 35.3] | 14.5 [13.1, 15.9] |
| R1E | 72.3 [67.0, 76.9] | 40.2 [37.7, 42.6] | 60.9 [57.1, 64.5] | 31.2 [29.2, 33.2] |
| R1E+c | 77.2 [71.8, 82.2] | 42.6 [40.1, 45.1] | 63.3 [59.6, 66.9] | 32.5 [30.6, 34.5] |
| R1J | 60.7 [52.0, 70.6] | 37.8 [35.2, 40.4] | 56.8 [50.4, 62.4] | 29.8 [27.8, 31.9] |
| R1J+c | 66.2 [57.6, 76.1] | 40.4 [37.9, 43.2] | 58.4 [50.9, 64.5] | 31.4 [29.4, 33.5] |
| R2E | 73.6 [68.8, 77.8] | 45.3 [43.1, 47.6] | 62.9 [58.9, 66.5] | 34.8 [32.7, 36.7] |
| R2E+c | 78.3 [72.7, 83.1] | 47.8 [45.6, 50.1] | 64.6 [60.7, 68.3] | 36.1 [34.1, 38.0] |
| R2J | 62.4 [55.5, 70.5] | 41.4 [39.0, 43.8] | 48.8 [40.9, 57.1] | 31.0 [28.9, 33.0] |
| R2J+c | 68.0 [61.2, 76.0] | 44.3 [41.9, 46.8] | 50.6 [42.2, 59.4] | 32.4 [30.2, 34.6] |

Number of configurations scored on the development frames (R2): MT0.N1.LD 108 structural settings, 14618 (structure, post-filter) configurations; MT0.N1.HD 108 structural settings, 12474 (structure, post-filter) configurations; MT0.N2.LD 108 structural settings, 13278 (structure, post-filter) configurations; MT0.N2.HD 108 structural settings, 13413 (structure, post-filter) configurations.

## MT0.N1.LD

Development frames 0-999 (1000 frames, 1253 ground-truth activations), test frames 1000-1499 (500 frames, 597 activations). Frame indices are 0-based stack indices. Ground-truth activations of the whole set: 24938; photon threshold of the challenge setting (25 % quantile of the whole set's ground truth): 1851 photons. Persistence: a fluorophore stays on for 2.27 consecutive frames on average (10989 runs, 31 % of them 3 frames or longer), so neighbouring frames are correlated, which the block bootstrap accounts for.

### R1: tuning curve of the minimum net gradient on the development frames (20 values)

| min net gradient | Jaccard % | precision % | recall % | RMSE nm | efficiency | detections / frame |
|---|---|---|---|---|---|---|
| 1000 | 4.0 | 4.1 | 80.6 | 51.0 | -8.7 | 24.71 |
| 1500 | 8.2 | 8.4 | 79.9 | 49.7 | -4.4 | 11.97 |
| 2000 | 19.5 | 20.6 | 79.2 | 48.3 | 6.1 | 4.82 |
| 2500 | 42.8 | 48.7 | 77.9 | 45.5 | 26.9 | 2.00 |
| 3000 | 65.9 | 82.2 | 76.9 | 42.1 | 45.8 | 1.17 |
| 3250 | 71.0 | 91.5 | 76.0 | 40.1 | 50.5 | 1.04 |
| 3500 | 72.8 | 95.8 | 75.2 | 38.2 | 53.1 | 0.98 |
| 3750 | 73.8 | 98.2 | 74.8 | 38.5 | 53.4 | 0.95 |
| **4000** | 73.8 | 99.1 | 74.3 | 37.6 | 54.2 | 0.94 |
| 4500 | 72.9 | 99.7 | 73.0 | 35.1 | 55.6 | 0.92 |
| 5000 | 72.0 | 99.9 | 72.1 | 34.5 | 55.6 | 0.90 |
| 5500 | 70.9 | 100.0 | 70.9 | 34.3 | 55.1 | 0.89 |
| 6000 | 69.9 | 100.0 | 69.9 | 33.3 | 55.1 | 0.88 |
| 7000 | 68.0 | 100.0 | 68.0 | 31.0 | 55.4 | 0.85 |
| **8000** | 66.6 | 100.0 | 66.6 | 29.0 | 55.8 | 0.83 |
| 9000 | 65.4 | 100.0 | 65.4 | 28.4 | 55.2 | 0.82 |
| 10000 | 63.8 | 100.0 | 63.8 | 26.6 | 55.0 | 0.80 |
| 12000 | 59.9 | 100.0 | 59.9 | 24.0 | 53.3 | 0.75 |
| 15000 | 55.5 | 100.0 | 55.5 | 22.4 | 50.2 | 0.70 |
| 20000 | 48.9 | 100.0 | 48.9 | 22.1 | 44.4 | 0.61 |

Bold: the value maximising the efficiency (8000; R1E) and the value maximising the Jaccard index (4000; R1J). The curve is computed from one run at the smallest grid value, filtered at each threshold (identification keeps the maxima with net gradient above the threshold and every spot is fitted independently); this was verified to reproduce direct runs at 8000 and 20000 on the first 60 development frames.

How well defined are these two operating points? On the development frames 8 of the 20 thresholds (from 4500 to 10000) have an efficiency within 1 point of the best (55.8): the efficiency is a plateau and its argmax is arbitrary within it. In 1000 block-bootstrap resamples of the development frames (blocks of 10 frames) the efficiency-optimal threshold is 8000 in 42 %, 4500 in 26 %, 5000 in 16 %; the Jaccard-optimal threshold is 4000 in 59 %, 3750 in 41 %, 3500 in 0 %. A falling RMSE with a rising threshold is largely survivor bias (the hard emitters are lost): compare the column "RMSE common TP" below, where the rows are scored on the same emitters.

### R2: search over all options on the development frames

Declared grid (identical for the sets of one noise class; at the top of `experiments/exp_smlm_picasso.py`): identification box 5/7/9 px, fit box 5/7/9 px, Gaussian filter sigma off/1/1.5 px, median padding 0/4 px, fit gausslq/gaussmle: **108 structural settings**, each run once at net gradient 300 and then post-filtered by a coordinate ascent (at most 3 rounds over, in this order: net gradient (27 values), minimum photons (9), minimum width (5), maximum width (7)) once for the efficiency and once for the Jaccard index. Development frames only.

**Configurations tried on the development frames: 108 structural settings run, 14618 distinct (structure, post-filter) configurations scored** (R1: 20 net-gradient values of the default structure, all contained in them); 5352859 candidate localizations fitted in total. The coordinate ascent converged (a round without change) within the 3 rounds in 189 of 216 runs.

Chosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); the offset is the mean signed error (detection - truth) of the matched pairs, the correction of the "+c" rows:

| rung | configuration | dev J % | dev precision % | dev recall % | dev RMSE nm | dev E | offset (dx, dy) nm | dev TP |
|---|---|---|---|---|---|---|---|---|
| R0 | box 7 px, raw frames, no padding, gausslq, net gradient > 5000 | 72.0 | 99.9 | 72.1 | 34.5 | 55.6 | (15.3, -11.5) | 903 |
| R1E | box 7 px, raw frames, no padding, gausslq, net gradient > 8000 | 66.6 | 100.0 | 66.6 | 29.0 | 55.8 | (14.1, -10.5) | 835 |
| R1J | box 7 px, raw frames, no padding, gausslq, net gradient > 4000 | 73.8 | 99.1 | 74.3 | 37.6 | 54.2 | (16.2, -11.2) | 931 |
| R2E | ident 5 / fit 9 px, Gaussian filter 1 px, median pad 4 px, gausslq, net gradient > 2500, 70 < width < 200 nm | 70.7 | 100.0 | 70.7 | 22.4 | 63.1 | (11.1, -8.7) | 886 |
| R2J | ident 7 / fit 9 px, Gaussian filter 1.5 px, median pad 4 px, gaussmle, net gradient > 1750 | 84.3 | 98.9 | 85.1 | 45.8 | 51.6 | (15.0, -11.1) | 1066 |

Ablation of R2E on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=5|fit_box=9|gauss_sigma=1.0|pad=4|fit=gausslq`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 70.7 | 63.1 | 22.4 |
| identification box -> 7 px | 72.1 | 62.9 | 24.5 |
| fit box -> the identification box | 67.3 | 58.7 | 25.3 |
| Gaussian filter -> off | 71.0 | 62.8 | 23.2 |
| median padding -> 0 | 69.0 | 59.6 | 25.9 |
| fit -> gausslq | (already the default) |  |  |
| photon filter -> off (the other post-filters kept, not re-tuned) | 70.7 | 63.1 | 22.4 |
| width window -> off (the other post-filters kept, not re-tuned) | 78.7 | 57.5 | 36.8 |
| every structural option -> default (post-filters re-tuned) | 70.2 | 58.6 | 28.8 |
| R1: net gradient only, 20 values | 66.6 | 55.8 | 29.0 |
| R0: Picasso defaults | 72.0 | 55.6 | 34.5 |

Ablation of R2J on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=7|fit_box=9|gauss_sigma=1.5|pad=4|fit=gaussmle`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 84.3 | 51.6 | 45.8 |
| identification box -> 7 px | (already the default) |  |  |
| fit box -> the identification box | 83.2 | 48.2 | 49.0 |
| Gaussian filter -> off | 80.5 | 53.6 | 42.1 |
| median padding -> 0 | 75.8 | 53.7 | 39.5 |
| fit -> gausslq | 83.8 | 51.9 | 45.3 |
| photon filter -> off (the other post-filters kept, not re-tuned) | 84.3 | 51.6 | 45.8 |
| width window -> off (the other post-filters kept, not re-tuned) | 84.3 | 51.6 | 45.8 |
| every structural option -> default (post-filters re-tuned) | 73.8 | 54.2 | 37.6 |
| R1: net gradient only, 20 values | 73.8 | 54.2 | 37.6 |
| R0: Picasso defaults | 72.0 | 55.6 | 34.5 |

Each option alone on top of R1 (the default structure with only the net gradient tuned, no other filter; the net gradient is tuned for the objective of the column), development frames:

| setting | dev E (net gradient tuned for E) | dev J % (net gradient tuned for J) |
|---|---|---|
| default structure, net gradient tuned (R1 on the finer grid) | 55.8 | 73.8 |
| box = 5 | 57.2 | 72.4 |
| box = 9 | 51.9 | 70.8 |
| fit_box = 5 | 54.6 | 73.5 |
| fit_box = 9 | 55.3 | 73.3 |
| gauss_sigma = 1.0 | 55.6 | 76.0 |
| gauss_sigma = 1.5 | 54.4 | 76.3 |
| pad = 4 | 58.3 | 80.6 |
| fit = gaussmle | 56.0 | 73.8 |

### Test frames: all ground-truth emitters (primary setting)

Rungs R0 (Picasso defaults), R1E / R1J (net gradient tuned on the development frames for efficiency / Jaccard), R2E / R2J (all options tuned on the development frames for efficiency / Jaccard), each without and with the "+c" offset correction (the constant (dx, dy) estimated on the development frames for that configuration, subtracted from every detection). 95 % intervals of a bootstrap in blocks of 10 consecutive frames (1000 resamples, seed 0; pooled metrics recomputed in every resample). The centred RMSE is the RMSE about the mean error of the matched pairs of the test frames. "RMSE common TP": RMSE over the 359 ground-truth emitters that EVERY one of the ten rows detects, so that no row is flattered by losing its hard emitters.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 71.7 [66.8, 76.3] | 100.0 [100.0, 100.0] | 71.7 [66.8, 76.3] | 34.7 [27.8, 41.3] | 29.1 [22.3, 35.2] | 19.9 [18.7, 21.1] | 55.2 [48.5, 61.3] |
| R0+c | 71.7 [66.8, 76.3] | 100.0 [100.0, 100.0] | 71.7 [66.8, 76.3] | 29.0 [22.3, 35.3] | 29.0 [22.2, 35.2] | 14.3 [13.4, 15.3] | 59.5 [53.3, 65.2] |
| R1E | 64.8 [60.0, 69.6] | 100.0 [100.0, 100.0] | 64.8 [60.0, 69.6] | 25.3 [21.4, 29.3] | 19.9 [15.9, 23.7] | 19.9 [18.7, 21.1] | 56.6 [51.6, 61.4] |
| R1E+c | 64.8 [60.0, 69.6] | 100.0 [100.0, 100.0] | 64.8 [60.0, 69.6] | 19.9 [16.0, 23.8] | 19.8 [15.8, 23.7] | 13.9 [13.0, 14.9] | 59.6 [54.6, 64.4] |
| R1J | 73.4 [68.1, 78.2] | 98.9 [97.9, 99.8] | 74.0 [68.9, 78.8] | 41.3 [33.0, 49.0] | 36.2 [28.0, 43.8] | 19.9 [18.7, 21.1] | 50.8 [43.7, 58.0] |
| R1J+c | 73.4 [68.1, 78.2] | 98.9 [97.9, 99.8] | 74.0 [68.9, 78.8] | 36.1 [28.0, 43.8] | 36.1 [27.9, 43.7] | 14.6 [13.6, 15.6] | 55.1 [48.0, 62.1] |
| R2E | 68.7 [63.8, 73.8] | 100.0 [100.0, 100.0] | 68.7 [63.8, 73.8] | 24.2 [20.8, 28.3] | 19.1 [15.2, 23.5] | 19.6 [18.5, 20.7] | 60.4 [55.2, 65.2] |
| R2E+c | 68.7 [63.8, 73.8] | 100.0 [100.0, 100.0] | 68.7 [63.8, 73.8] | 19.1 [15.2, 23.6] | 19.1 [15.2, 23.5] | 13.4 [12.4, 14.3] | 63.3 [58.0, 68.5] |
| R2J | 82.9 [79.5, 86.3] | 98.6 [97.4, 99.6] | 83.9 [80.9, 87.0] | 44.3 [37.3, 51.1] | 39.6 [32.6, 45.9] | 19.0 [17.6, 20.3] | 52.5 [46.1, 59.1] |
| R2J+c | 83.3 [79.8, 86.6] | 98.8 [97.6, 99.8] | 84.1 [81.0, 87.1] | 41.0 [33.7, 47.4] | 40.9 [33.7, 47.2] | 14.7 [13.6, 15.9] | 55.7 [49.6, 62.4] |

TP / FP / FN: R0 428 / 0 / 169; R1E 387 / 0 / 210; R1J 442 / 5 / 155; R2E 410 / 0 / 187; R2J 501 / 7 / 96. Detected per frame: R0 0.856; R1E 0.774; R1J 0.894; R2E 0.820; R2J 1.016; true per frame: 1.194. Offsets (dx, dy) in nm estimated on the development frames: R0 (15.3, -11.5); R1E (14.1, -10.5); R1J (16.2, -11.2); R2E (11.1, -8.7); R2J (15.0, -11.1).

Paired differences between rows on the same test frames (both evaluated on every bootstrap resample, blocks of 10 frames; 95 % interval; an interval that contains 0 is not distinguishable from frame-sampling noise):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | -6.9 [-9.2, -4.8] | -9.4 [-14.7, -4.6] | +1.5 [-1.5, +4.8] |
| R1J minus R0 | +1.7 [+0.3, +3.4] | +6.6 [+2.1, +11.2] | -4.3 [-7.7, -1.2] |
| R2E minus R0 | -3.0 [-7.3, +1.5] | -10.5 [-17.1, -4.6] | +5.2 [+1.9, +9.3] |
| R2J minus R0 | +11.3 [+7.7, +14.8] | +9.6 [+1.6, +17.4] | -2.7 [-9.3, +4.3] |
| R2E minus R1E | +3.9 [-0.0, +7.9] | -1.1 [-6.5, +4.8] | +3.8 [+0.4, +7.2] |
| R2J minus R1J | +9.5 [+6.0, +13.5] | +3.0 [-5.2, +11.1] | +1.7 [-4.7, +9.0] |
| R1E+c minus R0+c | -6.9 [-9.2, -4.8] | -9.1 [-14.6, -4.4] | +0.1 [-2.5, +3.2] |
| R1J+c minus R0+c | +1.7 [+0.3, +3.4] | +7.1 [+2.6, +11.6] | -4.3 [-7.5, -1.3] |
| R2E+c minus R0+c | -3.0 [-7.3, +1.5] | -9.9 [-16.4, -4.2] | +3.8 [+0.8, +7.5] |
| R2J+c minus R0+c | +11.6 [+8.0, +15.1] | +12.0 [+4.6, +19.2] | -3.7 [-9.5, +2.7] |
| R2E+c minus R1E+c | +3.9 [-0.0, +7.9] | -0.8 [-6.4, +5.4] | +3.7 [+0.6, +7.0] |
| R2J+c minus R1J+c | +9.8 [+6.2, +13.7] | +4.8 [-2.8, +12.7] | +0.6 [-5.7, +7.2] |
| R0+c minus R0 | +0.0 [+0.0, +0.0] | -5.7 [-6.6, -4.9] | +4.3 [+3.4, +5.1] |
| R1E+c minus R1E | +0.0 [+0.0, +0.0] | -5.5 [-6.4, -4.5] | +2.9 [+2.3, +3.6] |
| R1J+c minus R1J | +0.0 [+0.0, +0.0] | -5.2 [-6.1, -4.3] | +4.3 [+3.4, +5.1] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -5.1 [-6.1, -4.3] | +2.9 [+2.4, +3.5] |
| R2J+c minus R2J | +0.3 [+0.0, +0.9] | -3.4 [-5.3, -0.7] | +3.2 [+1.0, +5.0] |

### Test frames: challenge assessment setting (emitters within 450 nm of the border and the dimmest 25 % of the whole set's activations dropped after pairing)

The same ten rows; "RMSE common TP" over the 311 kept emitters that every row detects.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 91.8 [88.0, 95.6] | 100.0 [100.0, 100.0] | 91.8 [88.0, 95.6] | 32.2 [25.4, 38.8] | 26.2 [19.3, 32.3] | 19.0 [17.6, 20.2] | 66.8 [59.8, 73.6] |
| R0+c | 92.1 [88.4, 95.7] | 100.0 [100.0, 100.0] | 92.1 [88.4, 95.7] | 26.3 [19.6, 32.6] | 26.3 [19.3, 32.5] | 12.8 [11.9, 13.7] | 72.6 [66.1, 79.3] |
| R1E | 87.1 [82.1, 91.6] | 100.0 [100.0, 100.0] | 87.1 [82.1, 91.6] | 24.5 [20.4, 28.9] | 18.7 [14.0, 22.9] | 19.0 [17.6, 20.2] | 72.3 [67.0, 76.9] |
| R1E+c | 87.3 [82.8, 91.6] | 100.0 [100.0, 100.0] | 87.3 [82.8, 91.6] | 19.0 [14.6, 23.3] | 18.9 [14.4, 23.2] | 12.3 [11.5, 13.3] | 77.2 [71.8, 82.2] |
| R1J | 91.9 [87.9, 95.4] | 98.6 [97.3, 99.7] | 93.1 [89.5, 96.4] | 38.4 [28.5, 46.9] | 32.8 [22.8, 41.2] | 19.0 [17.6, 20.2] | 60.7 [52.0, 70.6] |
| R1J+c | 92.2 [88.4, 95.5] | 98.6 [97.3, 99.7] | 93.4 [90.0, 96.5] | 32.9 [22.9, 41.4] | 32.9 [22.9, 41.2] | 13.0 [12.2, 14.0] | 66.2 [57.6, 76.1] |
| R2E | 83.9 [78.0, 89.5] | 100.0 [100.0, 100.0] | 83.9 [78.0, 89.5] | 21.0 [18.3, 24.8] | 14.6 [11.6, 18.7] | 18.8 [17.4, 20.0] | 73.6 [68.8, 77.8] |
| R2E+c | 83.9 [78.0, 89.5] | 100.0 [100.0, 100.0] | 83.9 [78.0, 89.5] | 14.6 [11.6, 18.7] | 14.5 [11.6, 18.7] | 11.9 [10.9, 12.8] | 78.3 [72.7, 83.1] |
| R2J | 95.8 [94.1, 97.8] | 99.5 [98.6, 100.0] | 96.3 [94.4, 98.3] | 37.4 [29.2, 44.2] | 31.8 [23.7, 38.3] | 18.2 [16.5, 19.8] | 62.4 [55.5, 70.5] |
| R2J+c | 95.8 [94.1, 97.8] | 99.5 [98.6, 100.0] | 96.3 [94.4, 98.3] | 31.7 [23.7, 38.4] | 31.7 [23.6, 38.3] | 13.5 [12.3, 14.6] | 68.0 [61.2, 76.0] |

Ground-truth emitters counted: 379; detections counted per row: R0 348; R1E 330; R1J 358; R2E 318; R2J 367.

Paired differences (challenge setting):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | -4.7 [-7.3, -2.6] | -7.6 [-12.6, -3.3] | +5.5 [+1.9, +9.2] |
| R1J minus R0 | +0.1 [-1.7, +1.8] | +6.3 [+0.7, +11.3] | -6.1 [-10.8, -1.1] |
| R2E minus R0 | -7.9 [-11.8, -4.2] | -11.2 [-17.6, -4.9] | +6.8 [+2.2, +11.1] |
| R2J minus R0 | +4.0 [+0.8, +7.4] | +5.2 [-1.8, +12.6] | -4.4 [-11.5, +2.2] |
| R2E minus R1E | -3.2 [-6.0, -0.3] | -3.6 [-8.9, +2.1] | +1.3 [-2.7, +5.1] |
| R2J minus R1J | +3.9 [+0.9, +7.6] | -1.0 [-7.5, +4.7] | +1.6 [-3.9, +7.8] |
| R1E+c minus R0+c | -4.7 [-7.3, -2.6] | -7.3 [-11.8, -3.2] | +4.6 [+1.4, +8.0] |
| R1J+c minus R0+c | +0.1 [-1.7, +1.8] | +6.6 [+0.7, +11.8] | -6.4 [-11.0, -1.1] |
| R2E+c minus R0+c | -8.2 [-12.2, -4.4] | -11.7 [-18.2, -5.4] | +5.7 [+1.6, +9.6] |
| R2J+c minus R0+c | +3.7 [+0.6, +6.9] | +5.5 [-1.3, +12.8] | -4.6 [-11.4, +1.8] |
| R2E+c minus R1E+c | -3.4 [-6.3, -0.6] | -4.4 [-10.1, +1.7] | +1.1 [-2.8, +4.6] |
| R2J+c minus R1J+c | +3.6 [+0.7, +6.9] | -1.2 [-8.1, +4.9] | +1.8 [-4.1, +8.5] |
| R0+c minus R0 | +0.3 [+0.0, +1.0] | -5.9 [-6.8, -4.9] | +5.8 [+4.7, +6.7] |
| R1E+c minus R1E | +0.3 [+0.0, +1.0] | -5.5 [-6.7, -4.4] | +4.9 [+3.9, +6.0] |
| R1J+c minus R1J | +0.3 [+0.0, +0.9] | -5.5 [-6.6, -4.6] | +5.5 [+4.5, +6.5] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -6.4 [-7.5, -5.5] | +4.7 [+3.8, +5.8] |
| R2J+c minus R2J | +0.0 [+0.0, +0.0] | -5.7 [-6.6, -4.7] | +5.6 [+4.7, +6.5] |

Pairing precedes the exclusions, so an excluded (dim or border) emitter can take the detection of a kept one. Kept emitters left unmatched although a detection within 250 nm was paired with an excluded emitter, of all kept emitters: R0 2 of 379; R1E 1 of 379; R1J 2 of 379; R2E 0 of 379; R2J 0 of 379.

### Effect of each exclusion on the same detections (Jaccard % / RMSE nm / efficiency, point estimates)

| row | all emitters | border 450 nm dropped | dimmest 25 % dropped | both (challenge setting) |
|---|---|---|---|---|
| R0 | 71.7 / 34.7 / 55.2 | 81.2 / 32.1 / 62.8 | 82.2 / 35.3 / 60.4 | 91.8 / 32.2 / 66.8 |
| R0+c | 71.7 / 29.0 / 59.5 | 81.2 / 26.3 / 67.7 | 82.4 / 29.4 / 65.7 | 92.1 / 26.3 / 72.6 |
| R1E | 64.8 / 25.3 / 56.6 | 73.9 / 25.2 / 63.7 | 77.4 / 24.8 / 66.5 | 87.1 / 24.5 / 72.3 |
| R1E+c | 64.8 / 19.9 / 59.6 | 73.9 / 19.7 / 67.3 | 77.7 / 19.2 / 70.5 | 87.3 / 19.0 / 77.2 |
| R1J | 73.4 / 41.3 / 50.8 | 82.7 / 38.0 / 58.3 | 82.8 / 42.4 / 54.2 | 91.9 / 38.4 / 60.7 |
| R1J+c | 73.4 / 36.1 / 55.1 | 82.7 / 33.0 / 62.8 | 83.0 / 36.7 / 59.6 | 92.2 / 32.9 / 66.2 |
| R2E | 68.7 / 24.2 / 60.4 | 74.7 / 22.4 / 66.1 | 77.9 / 22.7 / 68.3 | 83.9 / 21.0 / 73.6 |
| R2E+c | 68.7 / 19.1 / 63.3 | 74.7 / 16.7 / 69.6 | 77.9 / 17.0 / 72.1 | 83.9 / 14.6 / 78.3 |
| R2J | 82.9 / 44.3 / 52.5 | 87.0 / 37.9 / 60.0 | 92.4 / 45.0 / 54.3 | 95.8 / 37.4 / 62.4 |
| R2J+c | 83.3 / 41.0 / 55.7 | 87.0 / 33.1 / 64.4 | 92.9 / 41.3 / 58.1 | 95.8 / 31.7 / 68.0 |

Picasso never searches the box // 2 pixel rows / columns at the top and left (box // 2 + 1 at the bottom and right: 3 and 4 px for the default 7 px box) when the frame is not padded, so emitters there are missed by construction; the border setting removes most of that zone from the ground truth (450 nm = 4.5 px). Median padding (an R2 option) removes the blind zone.

### Recall by distance to the border, by photons and by |z| (test frames, rungs without the offset correction)

By distance of the emitter to the nearest field edge (all emitters):

| distance to the edge | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-300 nm | 44 | 6.8 | 2.3 | 11.4 | 20.5 | 50.0 |
| 300-450 nm | 48 | 31.2 | 27.1 | 31.2 | 50.0 | 79.2 |
| 450-700 nm | 51 | 58.8 | 47.1 | 62.7 | 41.2 | 76.5 |
| 700-1000 nm | 46 | 63.0 | 52.2 | 69.6 | 60.9 | 87.0 |
| 1000-7000 nm | 408 | 86.0 | 79.7 | 87.7 | 80.4 | 88.7 |

By emitted photons (emitters at least 450 nm from the edge):

| photons | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-500 photons | 27 | 3.7 | 0.0 | 3.7 | 0.0 | 7.4 |
| 500-1000 photons | 38 | 31.6 | 5.3 | 47.4 | 31.6 | 60.5 |
| 1000-2000 photons | 69 | 82.6 | 71.0 | 84.1 | 79.7 | 85.5 |
| 2000-4000 photons | 98 | 88.8 | 80.6 | 88.8 | 82.7 | 93.9 |
| 4000-7000 photons | 128 | 86.7 | 82.0 | 90.6 | 75.8 | 94.5 |
| 7000-inf photons | 145 | 97.9 | 95.2 | 97.9 | 91.0 | 99.3 |

By axial distance |z| from the focal plane (emitters at least 450 nm from the edge and at least 1851 photons; the PSF broadens with |z|, see `calibration.md`):

| |z| | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-150 nm | 217 | 100.0 | 100.0 | 100.0 | 99.5 | 99.1 |
| 150-300 nm | 82 | 96.3 | 93.9 | 96.3 | 91.5 | 93.9 |
| 300-450 nm | 47 | 91.5 | 70.2 | 91.5 | 55.3 | 100.0 |
| 450-1000 nm | 33 | 27.3 | 9.1 | 42.4 | 3.0 | 78.8 |

### Picasso photons against the emitted photons, and what a one-peak-per-blob detector cannot do

Photons reported by Picasso (QE-corrected) divided by the ground-truth photons of the matched emitter, test frames, median and inter-quartile range (the Gaussian core holds only part of the emitted photons, `calibration.md` section on the halo, so the photon count of a Picasso fit is biased low and more so out of focus):

| row | emitters | n | median ratio | IQR |
|---|---|---|---|---|
| R0 | all matched | 428 | 0.78 | 0.67 .. 0.87 |
| R0 | |z| 0-100 nm | 146 | 0.87 | 0.82 .. 0.92 |
| R0 | |z| 100-200 nm | 179 | 0.76 | 0.69 .. 0.83 |
| R0 | |z| 200-350 nm | 79 | 0.59 | 0.50 .. 0.67 |
| R0 | |z| 350-800 nm | 24 | 0.49 | 0.34 .. 0.89 |
| R2E | all matched | 410 | 0.77 | 0.69 .. 0.85 |
| R2E | |z| 0-100 nm | 145 | 0.85 | 0.82 .. 0.87 |
| R2E | |z| 100-200 nm | 192 | 0.76 | 0.70 .. 0.81 |
| R2E | |z| 200-350 nm | 67 | 0.61 | 0.55 .. 0.65 |
| R2E | |z| 350-800 nm | 6 | 0.50 | 0.36 .. 0.54 |

Structural ceiling of a detector that reports one peak per blob: 3.4 % of the 597 ground-truth rows of the test frames have another emitter of the same frame within 250 nm and 3.7 % within 400 nm; on the development frames, keeping EVERY candidate of the chosen structure (R2E: 88.3 per frame, R2J: 42.5 per frame) the recall is at most 92.8 % (R2E structure) / 90.3 % (R2J structure) (pairing with the maximum number of pairs). The development and test frames show the same three microtubules, and `positions.csv` holds the exact emitter sites of the fluorophores (no frame information), which a method that uses the structure as a prior could exploit; Picasso uses neither.

## MT0.N1.HD

Development frames 0-99 (100 frames, 1294 ground-truth activations), test frames 100-499 (400 frames, 5046 activations). Frame indices are 0-based stack indices. Ground-truth activations of the whole set: 31612; photon threshold of the challenge setting (25 % quantile of the whole set's ground truth): 1829 photons. Persistence: a fluorophore stays on for 2.26 consecutive frames on average (13981 runs, 31 % of them 3 frames or longer), so neighbouring frames are correlated, which the block bootstrap accounts for.

### R1: tuning curve of the minimum net gradient on the development frames (20 values)

| min net gradient | Jaccard % | precision % | recall % | RMSE nm | efficiency | detections / frame |
|---|---|---|---|---|---|---|
| 1000 | 24.3 | 31.4 | 51.9 | 58.5 | 4.3 | 21.38 |
| 1500 | 33.5 | 48.9 | 51.5 | 56.6 | 12.6 | 13.62 |
| 2000 | 42.9 | 72.7 | 51.2 | 55.4 | 20.4 | 9.11 |
| 2500 | 48.2 | 90.3 | 50.9 | 55.3 | 24.2 | 7.29 |
| **3000** | 49.6 | 97.3 | 50.3 | 54.0 | 26.2 | 6.69 |
| 3250 | 49.5 | 98.3 | 49.9 | 53.1 | 26.8 | 6.57 |
| 3500 | 49.2 | 98.5 | 49.5 | 52.9 | 26.6 | 6.51 |
| 3750 | 49.0 | 98.9 | 49.3 | 52.4 | 26.9 | 6.45 |
| 4000 | 49.0 | 99.1 | 49.2 | 52.3 | 27.0 | 6.43 |
| 4500 | 48.6 | 99.2 | 48.8 | 51.9 | 26.9 | 6.36 |
| 5000 | 48.3 | 99.5 | 48.4 | 51.8 | 26.8 | 6.29 |
| 5500 | 47.6 | 99.7 | 47.7 | 51.0 | 26.9 | 6.19 |
| **6000** | 47.2 | 99.7 | 47.3 | 49.5 | 27.6 | 6.14 |
| 7000 | 46.7 | 99.8 | 46.8 | 49.1 | 27.5 | 6.06 |
| 8000 | 45.0 | 99.8 | 45.1 | 47.7 | 27.2 | 5.84 |
| 9000 | 44.3 | 100.0 | 44.3 | 46.9 | 27.1 | 5.73 |
| 10000 | 43.4 | 100.0 | 43.4 | 45.7 | 27.2 | 5.61 |
| 12000 | 41.2 | 100.0 | 41.2 | 43.0 | 27.2 | 5.33 |
| 15000 | 38.9 | 100.0 | 38.9 | 39.2 | 27.4 | 5.03 |
| 20000 | 34.2 | 100.0 | 34.2 | 38.9 | 23.6 | 4.43 |

Bold: the value maximising the efficiency (6000; R1E) and the value maximising the Jaccard index (3000; R1J). The curve is computed from one run at the smallest grid value, filtered at each threshold (identification keeps the maxima with net gradient above the threshold and every spot is fitted independently); this was verified to reproduce direct runs at 6000 and 20000 on the first 60 development frames.

How well defined are these two operating points? On the development frames 14 of the 20 thresholds (from 3250 to 15000) have an efficiency within 1 point of the best (27.6): the efficiency is a plateau and its argmax is arbitrary within it. In 1000 block-bootstrap resamples of the development frames (blocks of 10 frames) the efficiency-optimal threshold is 15000 in 28 %, 6000 in 24 %, 10000 in 11 %; the Jaccard-optimal threshold is 3000 in 73 %, 3250 in 26 %, 3500 in 0 %. A falling RMSE with a rising threshold is largely survivor bias (the hard emitters are lost): compare the column "RMSE common TP" below, where the rows are scored on the same emitters.

### R2: search over all options on the development frames

Declared grid (identical for the sets of one noise class; at the top of `experiments/exp_smlm_picasso.py`): identification box 5/7/9 px, fit box 5/7/9 px, Gaussian filter sigma off/1/1.5 px, median padding 0/4 px, fit gausslq/gaussmle: **108 structural settings**, each run once at net gradient 300 and then post-filtered by a coordinate ascent (at most 3 rounds over, in this order: net gradient (27 values), minimum photons (9), minimum width (5), maximum width (7)) once for the efficiency and once for the Jaccard index. Development frames only.

**Configurations tried on the development frames: 108 structural settings run, 12474 distinct (structure, post-filter) configurations scored** (R1: 20 net-gradient values of the default structure, all contained in them); 353483 candidate localizations fitted in total. The coordinate ascent converged (a round without change) within the 3 rounds in 207 of 216 runs.

Chosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); the offset is the mean signed error (detection - truth) of the matched pairs, the correction of the "+c" rows:

| rung | configuration | dev J % | dev precision % | dev recall % | dev RMSE nm | dev E | offset (dx, dy) nm | dev TP |
|---|---|---|---|---|---|---|---|---|
| R0 | box 7 px, raw frames, no padding, gausslq, net gradient > 5000 | 48.3 | 99.5 | 48.4 | 51.8 | 26.8 | (13.9, -9.9) | 626 |
| R1E | box 7 px, raw frames, no padding, gausslq, net gradient > 6000 | 47.2 | 99.7 | 47.3 | 49.5 | 27.6 | (13.1, -10.0) | 612 |
| R1J | box 7 px, raw frames, no padding, gausslq, net gradient > 3000 | 49.6 | 97.3 | 50.3 | 54.0 | 26.2 | (14.2, -9.3) | 651 |
| R2E | box 5 px, raw frames, median pad 4 px, gausslq, net gradient > 5500, 0 < width < 500 nm | 52.0 | 100.0 | 52.0 | 39.6 | 37.8 | (10.7, -9.9) | 673 |
| R2J | box 5 px, raw frames, median pad 4 px, gausslq, net gradient > 3000 | 57.6 | 98.8 | 58.0 | 53.5 | 31.8 | (13.8, -10.6) | 751 |

Ablation of R2E on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=5|fit_box=5|gauss_sigma=None|pad=4|fit=gausslq`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 52.0 | 37.8 | 39.6 |
| identification box -> 7 px | 47.6 | 36.2 | 36.4 |
| fit box -> the identification box | (already the default) |  |  |
| Gaussian filter -> off | (already the default) |  |  |
| median padding -> 0 | 49.4 | 35.9 | 39.4 |
| fit -> gausslq | (already the default) |  |  |
| photon filter -> off (the other post-filters kept, not re-tuned) | 52.0 | 37.8 | 39.6 |
| width window -> off (the other post-filters kept, not re-tuned) | 52.6 | 37.3 | 41.0 |
| every structural option -> default (post-filters re-tuned) | 44.2 | 31.6 | 39.5 |
| R1: net gradient only, 20 values | 47.2 | 27.6 | 49.5 |
| R0: Picasso defaults | 48.3 | 26.8 | 51.8 |

Ablation of R2J on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=5|fit_box=5|gauss_sigma=None|pad=4|fit=gausslq`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 57.6 | 31.8 | 53.5 |
| identification box -> 7 px | 55.7 | 29.5 | 54.8 |
| fit box -> the identification box | (already the default) |  |  |
| Gaussian filter -> off | (already the default) |  |  |
| median padding -> 0 | 53.6 | 31.1 | 50.9 |
| fit -> gausslq | (already the default) |  |  |
| photon filter -> off (the other post-filters kept, not re-tuned) | 57.6 | 31.8 | 53.5 |
| width window -> off (the other post-filters kept, not re-tuned) | 57.6 | 31.8 | 53.5 |
| every structural option -> default (post-filters re-tuned) | 49.7 | 26.3 | 53.8 |
| R1: net gradient only, 20 values | 49.6 | 26.2 | 54.0 |
| R0: Picasso defaults | 48.3 | 26.8 | 51.8 |

Each option alone on top of R1 (the default structure with only the net gradient tuned, no other filter; the net gradient is tuned for the objective of the column), development frames:

| setting | dev E (net gradient tuned for E) | dev J % (net gradient tuned for J) |
|---|---|---|
| default structure, net gradient tuned (R1 on the finer grid) | 27.6 | 49.6 |
| box = 5 | 27.3 | 51.4 |
| box = 9 | 26.4 | 43.9 |
| fit_box = 5 | 33.2 | 49.6 |
| fit_box = 9 | 18.7 | 47.0 |
| gauss_sigma = 1.0 | 27.6 | 46.6 |
| gauss_sigma = 1.5 | 25.7 | 44.8 |
| pad = 4 | 30.1 | 55.8 |
| fit = gaussmle | 25.8 | 49.5 |

### Test frames: all ground-truth emitters (primary setting)

Rungs R0 (Picasso defaults), R1E / R1J (net gradient tuned on the development frames for efficiency / Jaccard), R2E / R2J (all options tuned on the development frames for efficiency / Jaccard), each without and with the "+c" offset correction (the constant (dx, dy) estimated on the development frames for that configuration, subtracted from every detection). 95 % intervals of a bootstrap in blocks of 10 consecutive frames (1000 resamples, seed 0; pooled metrics recomputed in every resample). The centred RMSE is the RMSE about the mean error of the matched pairs of the test frames. "RMSE common TP": RMSE over the 2092 ground-truth emitters that EVERY one of the ten rows detects, so that no row is flattered by losing its hard emitters.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 47.8 [46.3, 49.3] | 99.8 [99.7, 100.0] | 47.8 [46.4, 49.3] | 49.6 [47.5, 51.8] | 46.6 [44.5, 48.8] | 43.6 [41.7, 45.8] | 28.0 [25.9, 30.1] |
| R0+c | 47.8 [46.3, 49.3] | 99.9 [99.7, 100.0] | 47.8 [46.4, 49.4] | 46.4 [44.3, 48.7] | 46.3 [44.2, 48.6] | 40.3 [38.5, 42.3] | 30.2 [28.1, 32.3] |
| R1E | 46.9 [45.5, 48.4] | 99.9 [99.8, 100.0] | 46.9 [45.5, 48.4] | 48.5 [46.5, 50.7] | 45.7 [43.7, 47.8] | 43.6 [41.7, 45.6] | 28.1 [25.9, 30.2] |
| R1E+c | 46.9 [45.5, 48.5] | 100.0 [99.9, 100.0] | 46.9 [45.5, 48.5] | 45.5 [43.4, 47.7] | 45.4 [43.3, 47.6] | 40.3 [38.5, 42.3] | 30.1 [28.0, 32.3] |
| R1J | 49.0 [47.6, 50.5] | 97.7 [97.1, 98.3] | 49.6 [48.1, 51.1] | 53.4 [51.1, 55.8] | 50.3 [48.0, 52.7] | 43.8 [41.8, 45.9] | 26.2 [24.0, 28.2] |
| R1J+c | 49.0 [47.6, 50.5] | 97.8 [97.2, 98.3] | 49.6 [48.1, 51.1] | 50.1 [47.7, 52.7] | 49.9 [47.6, 52.5] | 40.5 [38.5, 42.4] | 28.5 [26.4, 30.6] |
| R2E | 51.5 [49.9, 53.0] | 99.9 [99.8, 100.0] | 51.5 [49.9, 53.0] | 43.5 [41.7, 45.5] | 40.9 [39.1, 42.8] | 37.5 [35.6, 39.3] | 34.8 [33.0, 36.7] |
| R2E+c | 51.5 [49.9, 53.0] | 99.9 [99.8, 100.0] | 51.5 [49.9, 53.0] | 40.6 [38.7, 42.5] | 40.4 [38.6, 42.3] | 33.4 [31.5, 35.3] | 36.7 [35.0, 38.6] |
| R2J | 56.4 [54.8, 58.0] | 98.2 [97.6, 98.8] | 57.0 [55.5, 58.5] | 53.1 [50.8, 55.5] | 50.0 [47.8, 52.3] | 37.5 [35.6, 39.3] | 31.3 [29.2, 33.5] |
| R2J+c | 56.5 [55.0, 58.0] | 98.3 [97.7, 98.9] | 57.1 [55.6, 58.6] | 50.1 [47.8, 52.4] | 50.0 [47.7, 52.2] | 33.2 [31.3, 35.1] | 33.7 [31.6, 35.9] |

TP / FP / FN: R0 2412 / 4 / 2634; R1E 2368 / 2 / 2678; R1J 2501 / 58 / 2545; R2E 2598 / 2 / 2448; R2J 2877 / 52 / 2169. Detected per frame: R0 6.040; R1E 5.925; R1J 6.397; R2E 6.500; R2J 7.322; true per frame: 12.615. Offsets (dx, dy) in nm estimated on the development frames: R0 (13.9, -9.9); R1E (13.1, -10.0); R1J (14.2, -9.3); R2E (10.7, -9.9); R2J (13.8, -10.6).

Paired differences between rows on the same test frames (both evaluated on every bootstrap resample, blocks of 10 frames; 95 % interval; an interval that contains 0 is not distinguishable from frame-sampling noise):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | -0.9 [-1.1, -0.6] | -1.1 [-1.7, -0.5] | +0.1 [-0.2, +0.5] |
| R1J minus R0 | +1.2 [+0.8, +1.7] | +3.8 [+2.7, +5.1] | -1.8 [-2.6, -1.1] |
| R2E minus R0 | +3.7 [+2.5, +4.9] | -6.0 [-8.2, -4.0] | +6.8 [+5.4, +8.3] |
| R2J minus R0 | +8.7 [+7.6, +9.8] | +3.5 [+0.9, +6.0] | +3.3 [+1.6, +5.2] |
| R2E minus R1E | +4.6 [+3.4, +5.7] | -4.9 [-7.1, -2.9] | +6.7 [+5.2, +8.3] |
| R2J minus R1J | +7.4 [+6.2, +8.6] | -0.3 [-3.2, +2.3] | +5.2 [+3.4, +7.1] |
| R1E+c minus R0+c | -0.9 [-1.1, -0.6] | -0.9 [-1.5, -0.4] | -0.0 [-0.3, +0.3] |
| R1J+c minus R0+c | +1.2 [+0.8, +1.7] | +3.7 [+2.6, +4.9] | -1.6 [-2.4, -1.0] |
| R2E+c minus R0+c | +3.7 [+2.5, +4.9] | -5.8 [-8.0, -3.5] | +6.6 [+5.1, +8.1] |
| R2J+c minus R0+c | +8.7 [+7.6, +9.8] | +3.7 [+1.0, +6.2] | +3.5 [+1.9, +5.2] |
| R2E+c minus R1E+c | +4.5 [+3.4, +5.6] | -4.9 [-7.1, -2.8] | +6.6 [+5.2, +8.1] |
| R2J+c minus R1J+c | +7.5 [+6.3, +8.7] | -0.0 [-2.9, +2.8] | +5.2 [+3.4, +7.0] |
| R0+c minus R0 | +0.0 [+0.0, +0.1] | -3.2 [-3.8, -2.5] | +2.2 [+1.7, +2.6] |
| R1E+c minus R1E | +0.0 [+0.0, +0.1] | -3.0 [-3.6, -2.3] | +2.0 [+1.6, +2.4] |
| R1J+c minus R1J | +0.0 [+0.0, +0.1] | -3.3 [-3.9, -2.6] | +2.4 [+1.9, +2.8] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -3.0 [-3.3, -2.6] | +2.0 [+1.7, +2.2] |
| R2J+c minus R2J | +0.1 [+0.0, +0.2] | -3.0 [-3.6, -2.2] | +2.3 [+1.8, +2.8] |

### Test frames: challenge assessment setting (emitters within 450 nm of the border and the dimmest 25 % of the whole set's activations dropped after pairing)

The same ten rows; "RMSE common TP" over the 1944 kept emitters that every row detects.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 65.0 [63.1, 66.8] | 99.8 [99.6, 100.0] | 65.1 [63.2, 66.8] | 49.0 [46.8, 51.3] | 46.1 [43.9, 48.2] | 43.8 [41.8, 45.9] | 39.7 [37.3, 42.2] |
| R0+c | 65.0 [63.1, 66.9] | 99.9 [99.7, 100.0] | 65.1 [63.2, 66.9] | 45.9 [43.6, 48.2] | 45.8 [43.5, 48.0] | 40.5 [38.6, 42.6] | 42.3 [39.8, 44.8] |
| R1E | 64.5 [62.6, 66.3] | 99.9 [99.8, 100.0] | 64.6 [62.6, 66.4] | 48.1 [46.0, 50.4] | 45.3 [43.2, 47.5] | 43.7 [41.7, 45.9] | 40.2 [37.7, 42.6] |
| R1E+c | 64.6 [62.6, 66.4] | 100.0 [99.9, 100.0] | 64.6 [62.6, 66.4] | 45.2 [43.1, 47.4] | 45.1 [43.0, 47.3] | 40.4 [38.5, 42.6] | 42.6 [40.1, 45.1] |
| R1J | 64.7 [62.8, 66.4] | 97.5 [96.9, 98.2] | 65.8 [63.9, 67.6] | 51.3 [48.7, 53.8] | 48.2 [45.7, 50.7] | 43.8 [41.8, 45.9] | 37.8 [35.2, 40.4] |
| R1J+c | 64.7 [62.8, 66.5] | 97.6 [96.9, 98.2] | 65.8 [63.9, 67.7] | 48.0 [45.3, 50.5] | 47.8 [45.2, 50.4] | 40.5 [38.6, 42.6] | 40.4 [37.9, 43.2] |
| R2E | 65.8 [63.9, 67.7] | 99.9 [99.8, 100.0] | 65.8 [64.0, 67.7] | 42.6 [40.4, 44.8] | 39.8 [37.5, 42.0] | 37.5 [35.5, 39.5] | 45.3 [43.1, 47.6] |
| R2E+c | 65.9 [64.0, 67.8] | 99.9 [99.8, 100.0] | 65.9 [64.0, 67.8] | 39.5 [37.3, 41.7] | 39.3 [37.1, 41.4] | 33.4 [31.4, 35.4] | 47.8 [45.6, 50.1] |
| R2J | 68.3 [66.4, 70.1] | 98.7 [98.2, 99.1] | 68.9 [67.0, 70.8] | 49.3 [46.8, 51.7] | 46.1 [43.8, 48.4] | 37.5 [35.5, 39.5] | 41.4 [39.0, 43.8] |
| R2J+c | 68.4 [66.5, 70.2] | 98.7 [98.2, 99.2] | 69.0 [67.1, 70.9] | 45.9 [43.5, 48.2] | 45.8 [43.4, 48.1] | 33.2 [31.2, 35.2] | 44.3 [41.9, 46.8] |

Ground-truth emitters counted: 3314; detections counted per row: R0 2160; R1E 2142; R1J 2235; R2E 2184; R2J 2316.

Paired differences (challenge setting):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | -0.4 [-0.7, -0.2] | -0.9 [-1.5, -0.4] | +0.5 [+0.1, +0.9] |
| R1J minus R0 | -0.3 [-0.7, +0.1] | +2.2 [+1.2, +3.4] | -2.0 [-2.8, -1.2] |
| R2E minus R0 | +0.8 [-0.1, +1.8] | -6.4 [-8.5, -4.4] | +5.6 [+4.1, +7.1] |
| R2J minus R0 | +3.3 [+2.4, +4.3] | +0.3 [-2.1, +2.5] | +1.6 [-0.1, +3.4] |
| R2E minus R1E | +1.3 [+0.4, +2.2] | -5.5 [-7.6, -3.6] | +5.1 [+3.7, +6.6] |
| R2J minus R1J | +3.6 [+2.5, +4.6] | -2.0 [-4.6, +0.4] | +3.6 [+1.8, +5.6] |
| R1E+c minus R0+c | -0.4 [-0.7, -0.2] | -0.7 [-1.3, -0.2] | +0.3 [-0.0, +0.7] |
| R1J+c minus R0+c | -0.3 [-0.7, +0.1] | +2.1 [+1.1, +3.2] | -1.9 [-2.7, -1.2] |
| R2E+c minus R0+c | +0.8 [-0.1, +1.7] | -6.4 [-8.4, -4.3] | +5.5 [+4.0, +6.9] |
| R2J+c minus R0+c | +3.4 [+2.5, +4.3] | +0.0 [-2.1, +2.3] | +2.0 [+0.3, +3.6] |
| R2E+c minus R1E+c | +1.3 [+0.5, +2.1] | -5.7 [-7.8, -3.7] | +5.2 [+3.7, +6.7] |
| R2J+c minus R1J+c | +3.7 [+2.8, +4.6] | -2.1 [-4.6, +0.3] | +3.8 [+2.0, +5.7] |
| R0+c minus R0 | +0.0 [-0.2, +0.3] | -3.2 [-3.8, -2.4] | +2.6 [+2.1, +3.1] |
| R1E+c minus R1E | +0.0 [-0.2, +0.3] | -3.0 [-3.6, -2.2] | +2.4 [+1.9, +2.9] |
| R1J+c minus R1J | +0.0 [-0.2, +0.3] | -3.3 [-3.9, -2.5] | +2.7 [+2.1, +3.2] |
| R2E+c minus R2E | +0.1 [-0.1, +0.2] | -3.1 [-3.5, -2.7] | +2.5 [+2.1, +2.8] |
| R2J+c minus R2J | +0.1 [-0.1, +0.3] | -3.4 [-4.1, -2.6] | +2.9 [+2.2, +3.5] |

Pairing precedes the exclusions, so an excluded (dim or border) emitter can take the detection of a kept one. Kept emitters left unmatched although a detection within 250 nm was paired with an excluded emitter, of all kept emitters: R0 59 of 3314; R1E 58 of 3314; R1J 58 of 3314; R2E 51 of 3314; R2J 55 of 3314.

### Effect of each exclusion on the same detections (Jaccard % / RMSE nm / efficiency, point estimates)

| row | all emitters | border 450 nm dropped | dimmest 25 % dropped | both (challenge setting) |
|---|---|---|---|---|
| R0 | 47.8 / 49.6 / 28.0 | 54.2 / 49.4 / 32.7 | 58.0 / 49.2 / 35.3 | 65.0 / 49.0 / 39.7 |
| R0+c | 47.8 / 46.4 / 30.2 | 54.3 / 46.2 / 35.0 | 58.1 / 45.9 / 37.8 | 65.0 / 45.9 / 42.3 |
| R1E | 46.9 / 48.5 / 28.1 | 53.3 / 48.4 / 32.8 | 57.6 / 48.1 / 35.8 | 64.5 / 48.1 / 40.2 |
| R1E+c | 46.9 / 45.5 / 30.1 | 53.3 / 45.4 / 34.9 | 57.6 / 45.1 / 38.1 | 64.6 / 45.2 / 42.6 |
| R1J | 49.0 / 53.4 / 26.2 | 55.4 / 52.9 / 30.8 | 57.9 / 51.5 / 33.5 | 64.7 / 51.3 / 37.8 |
| R1J+c | 49.0 / 50.1 / 28.5 | 55.4 / 49.7 / 33.2 | 58.0 / 48.1 / 36.1 | 64.7 / 48.0 / 40.4 |
| R2E | 51.5 / 43.5 / 34.8 | 53.5 / 42.6 / 37.0 | 63.9 / 43.7 / 43.3 | 65.8 / 42.6 / 45.3 |
| R2E+c | 51.5 / 40.6 / 36.7 | 53.5 / 39.4 / 39.0 | 63.9 / 40.6 / 45.7 | 65.9 / 39.5 / 47.8 |
| R2J | 56.4 / 53.1 / 31.3 | 57.9 / 49.8 / 34.8 | 67.0 / 52.8 / 37.7 | 68.3 / 49.3 / 41.4 |
| R2J+c | 56.5 / 50.1 / 33.7 | 58.0 / 46.6 / 37.3 | 67.2 / 49.5 / 40.6 | 68.4 / 45.9 / 44.3 |

Picasso never searches the box // 2 pixel rows / columns at the top and left (box // 2 + 1 at the bottom and right: 3 and 4 px for the default 7 px box) when the frame is not padded, so emitters there are missed by construction; the border setting removes most of that zone from the ground truth (450 nm = 4.5 px). Median padding (an R2 option) removes the blind zone.

### Recall by distance to the border, by photons and by |z| (test frames, rungs without the offset correction)

By distance of the emitter to the nearest field edge (all emitters):

| distance to the edge | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-300 nm | 440 | 1.4 | 1.1 | 1.8 | 37.3 | 45.9 |
| 300-450 nm | 287 | 21.3 | 20.6 | 24.7 | 42.2 | 54.0 |
| 450-700 nm | 370 | 42.4 | 40.8 | 46.2 | 38.4 | 46.2 |
| 700-1000 nm | 412 | 51.0 | 49.8 | 54.4 | 45.9 | 53.9 |
| 1000-7000 nm | 3537 | 55.9 | 55.1 | 57.3 | 56.0 | 60.1 |

By emitted photons (emitters at least 450 nm from the edge):

| photons | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-500 photons | 277 | 3.2 | 3.2 | 6.5 | 2.9 | 5.4 |
| 500-1000 photons | 292 | 13.4 | 10.3 | 22.6 | 5.8 | 20.9 |
| 1000-2000 photons | 505 | 34.7 | 31.1 | 38.4 | 26.9 | 39.0 |
| 2000-4000 photons | 943 | 49.5 | 48.8 | 50.6 | 51.0 | 56.1 |
| 4000-7000 photons | 960 | 63.7 | 63.5 | 64.5 | 65.2 | 67.6 |
| 7000-inf photons | 1342 | 77.7 | 77.3 | 78.1 | 77.9 | 79.7 |

By axial distance |z| from the focal plane (emitters at least 450 nm from the edge and at least 1829 photons; the PSF broadens with |z|, see `calibration.md`):

| |z| | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-150 nm | 1877 | 68.8 | 68.7 | 68.8 | 73.6 | 73.7 |
| 150-300 nm | 800 | 68.2 | 68.1 | 68.2 | 69.0 | 71.1 |
| 300-450 nm | 402 | 58.5 | 57.2 | 59.7 | 51.2 | 60.2 |
| 450-1000 nm | 235 | 35.7 | 31.9 | 43.8 | 17.9 | 38.3 |

### Picasso photons against the emitted photons, and what a one-peak-per-blob detector cannot do

Photons reported by Picasso (QE-corrected) divided by the ground-truth photons of the matched emitter, test frames, median and inter-quartile range (the Gaussian core holds only part of the emitted photons, `calibration.md` section on the halo, so the photon count of a Picasso fit is biased low and more so out of focus):

| row | emitters | n | median ratio | IQR |
|---|---|---|---|---|
| R0 | all matched | 2412 | 0.81 | 0.66 .. 1.07 |
| R0 | |z| 0-100 nm | 716 | 0.91 | 0.83 .. 1.23 |
| R0 | |z| 100-200 nm | 952 | 0.79 | 0.70 .. 1.06 |
| R0 | |z| 200-350 nm | 531 | 0.61 | 0.51 .. 0.79 |
| R0 | |z| 350-800 nm | 213 | 0.55 | 0.37 .. 1.05 |
| R2E | all matched | 2598 | 0.91 | 0.67 .. 1.21 |
| R2E | |z| 0-100 nm | 769 | 1.02 | 0.87 .. 1.32 |
| R2E | |z| 100-200 nm | 1144 | 0.89 | 0.72 .. 1.21 |
| R2E | |z| 200-350 nm | 533 | 0.66 | 0.46 .. 1.06 |
| R2E | |z| 350-800 nm | 152 | 0.41 | 0.26 .. 1.01 |

Structural ceiling of a detector that reports one peak per blob: 35.0 % of the 5046 ground-truth rows of the test frames have another emitter of the same frame within 250 nm and 51.6 % within 400 nm; on the development frames, keeping EVERY candidate of the chosen structure (R2E: 89.2 per frame, R2J: 89.2 per frame) the recall is at most 64.8 % (R2E structure) / 64.8 % (R2J structure) (pairing with the maximum number of pairs). The development and test frames show the same three microtubules, and `positions.csv` holds the exact emitter sites of the fluorophores (no frame information), which a method that uses the structure as a prior could exploit; Picasso uses neither.

## MT0.N2.LD

Development frames 0-999 (1000 frames, 1253 ground-truth activations), test frames 1000-1499 (500 frames, 597 activations). Frame indices are 0-based stack indices. Ground-truth activations of the whole set: 24938; photon threshold of the challenge setting (25 % quantile of the whole set's ground truth): 200 photons. Persistence: a fluorophore stays on for 2.27 consecutive frames on average (10989 runs, 31 % of them 3 frames or longer), so neighbouring frames are correlated, which the block bootstrap accounts for.

### R1: tuning curve of the minimum net gradient on the development frames (20 values)

| min net gradient | Jaccard % | precision % | recall % | RMSE nm | efficiency | detections / frame |
|---|---|---|---|---|---|---|
| 100 | 2.0 | 2.0 | 76.6 | 74.2 | -23.0 | 48.77 |
| 150 | 2.1 | 2.2 | 76.5 | 73.9 | -22.6 | 44.50 |
| 200 | 2.4 | 2.4 | 76.1 | 73.2 | -22.0 | 39.77 |
| 250 | 2.7 | 2.7 | 75.7 | 72.8 | -21.5 | 34.86 |
| 300 | 3.1 | 3.1 | 74.9 | 71.6 | -20.5 | 29.91 |
| 400 | 4.4 | 4.5 | 73.5 | 68.3 | -17.5 | 20.65 |
| 500 | 6.7 | 6.9 | 72.0 | 65.9 | -14.2 | 13.04 |
| 600 | 11.1 | 11.6 | 70.8 | 64.5 | -9.8 | 7.62 |
| 800 | 30.6 | 36.0 | 67.0 | 56.4 | 10.6 | 2.33 |
| 1000 | 54.0 | 78.1 | 63.7 | 49.9 | 32.1 | 1.02 |
| **1250** | 58.8 | 98.0 | 59.5 | 45.0 | 39.0 | 0.76 |
| **1500** | 55.3 | 100.0 | 55.3 | 39.3 | 40.5 | 0.69 |
| 2000 | 50.2 | 100.0 | 50.2 | 33.6 | 39.9 | 0.63 |
| 2500 | 44.7 | 100.0 | 44.7 | 32.6 | 35.8 | 0.56 |
| 3000 | 39.1 | 100.0 | 39.1 | 31.3 | 31.5 | 0.49 |
| 3500 | 34.1 | 100.0 | 34.1 | 29.7 | 27.7 | 0.43 |
| 4000 | 29.8 | 100.0 | 29.8 | 29.0 | 24.1 | 0.37 |
| 5000 | 21.5 | 100.0 | 21.5 | 28.4 | 16.5 | 0.27 |
| 6000 | 14.6 | 100.0 | 14.6 | 28.1 | 10.1 | 0.18 |
| 8000 | 5.7 | 100.0 | 5.7 | 33.2 | 0.1 | 0.07 |

Bold: the value maximising the efficiency (1500; R1E) and the value maximising the Jaccard index (1250; R1J). The curve is computed from one run at the smallest grid value, filtered at each threshold (identification keeps the maxima with net gradient above the threshold and every spot is fitted independently); this was verified to reproduce direct runs at 1500 and 8000 on the first 60 development frames.

How well defined are these two operating points? On the development frames 2 of the 20 thresholds (from 1500 to 2000) have an efficiency within 1 point of the best (40.5): the efficiency is a plateau and its argmax is arbitrary within it. In 1000 block-bootstrap resamples of the development frames (blocks of 10 frames) the efficiency-optimal threshold is 1500 in 69 %, 2000 in 28 %, 1250 in 3 %; the Jaccard-optimal threshold is 1250 in 100 %, 100 in 0 %, 150 in 0 %. A falling RMSE with a rising threshold is largely survivor bias (the hard emitters are lost): compare the column "RMSE common TP" below, where the rows are scored on the same emitters.

### R2: search over all options on the development frames

Declared grid (identical for the sets of one noise class; at the top of `experiments/exp_smlm_picasso.py`): identification box 5/7/9 px, fit box 5/7/9 px, Gaussian filter sigma off/1/1.5 px, median padding 0/4 px, fit gausslq/gaussmle: **108 structural settings**, each run once at net gradient 30 and then post-filtered by a coordinate ascent (at most 3 rounds over, in this order: net gradient (28 values), minimum photons (9), minimum width (5), maximum width (7)) once for the efficiency and once for the Jaccard index. Development frames only.

**Configurations tried on the development frames: 108 structural settings run, 13278 distinct (structure, post-filter) configurations scored** (R1: 20 net-gradient values of the default structure, all contained in them); 6291256 candidate localizations fitted in total. The coordinate ascent converged (a round without change) within the 3 rounds in 211 of 216 runs.

Chosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); the offset is the mean signed error (detection - truth) of the matched pairs, the correction of the "+c" rows:

| rung | configuration | dev J % | dev precision % | dev recall % | dev RMSE nm | dev E | offset (dx, dy) nm | dev TP |
|---|---|---|---|---|---|---|---|---|
| R0 | box 7 px, raw frames, no padding, gausslq, net gradient > 5000 | 21.5 | 100.0 | 21.5 | 28.4 | 16.5 | (10.5, -8.4) | 269 |
| R1E | box 7 px, raw frames, no padding, gausslq, net gradient > 1500 | 55.3 | 100.0 | 55.3 | 39.3 | 40.5 | (12.4, -8.6) | 693 |
| R1J | box 7 px, raw frames, no padding, gausslq, net gradient > 1250 | 58.8 | 98.0 | 59.5 | 45.0 | 39.0 | (12.7, -8.7) | 746 |
| R2E | box 7 px, Gaussian filter 1 px, median pad 4 px, gaussmle, net gradient > 1000, photons > 80, 0 < width < 250 nm | 61.3 | 99.7 | 61.4 | 34.6 | 48.1 | (10.2, -7.6) | 769 |
| R2J | box 9 px, Gaussian filter 1.5 px, median pad 4 px, gausslq, net gradient > 700, photons > 60 | 70.0 | 96.2 | 72.0 | 52.0 | 39.9 | (12.7, -8.1) | 902 |

Ablation of R2E on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=7|fit_box=7|gauss_sigma=1.0|pad=4|fit=gaussmle`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 61.3 | 48.1 | 34.6 |
| identification box -> 7 px | (already the default) |  |  |
| fit box -> the identification box | (already the default) |  |  |
| Gaussian filter -> off | 60.3 | 45.7 | 37.1 |
| median padding -> 0 | 56.9 | 45.6 | 33.1 |
| fit -> gausslq | 62.9 | 46.6 | 38.5 |
| photon filter -> off (the other post-filters kept, not re-tuned) | 61.7 | 46.8 | 36.9 |
| width window -> off (the other post-filters kept, not re-tuned) | 64.7 | 46.2 | 40.6 |
| every structural option -> default (post-filters re-tuned) | 54.8 | 41.6 | 37.0 |
| R1: net gradient only, 20 values | 55.3 | 40.5 | 39.3 |
| R0: Picasso defaults | 21.5 | 16.5 | 28.4 |

Ablation of R2J on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=9|fit_box=9|gauss_sigma=1.5|pad=4|fit=gausslq`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 70.0 | 39.9 | 52.0 |
| identification box -> 7 px | 68.6 | 43.8 | 46.7 |
| fit box -> the identification box | (already the default) |  |  |
| Gaussian filter -> off | 65.3 | 42.4 | 45.9 |
| median padding -> 0 | 62.9 | 37.4 | 50.5 |
| fit -> gausslq | (already the default) |  |  |
| photon filter -> off (the other post-filters kept, not re-tuned) | 70.0 | 39.5 | 52.5 |
| width window -> off (the other post-filters kept, not re-tuned) | 70.0 | 39.9 | 52.0 |
| every structural option -> default (post-filters re-tuned) | 58.8 | 39.1 | 44.9 |
| R1: net gradient only, 20 values | 58.8 | 39.0 | 45.0 |
| R0: Picasso defaults | 21.5 | 16.5 | 28.4 |

Each option alone on top of R1 (the default structure with only the net gradient tuned, no other filter; the net gradient is tuned for the objective of the column), development frames:

| setting | dev E (net gradient tuned for E) | dev J % (net gradient tuned for J) |
|---|---|---|
| default structure, net gradient tuned (R1 on the finer grid) | 40.5 | 58.8 |
| box = 5 | 40.5 | 53.9 |
| box = 9 | 38.0 | 59.2 |
| fit_box = 5 | 34.3 | 57.8 |
| fit_box = 9 | 40.9 | 58.6 |
| gauss_sigma = 1.0 | 42.1 | 63.3 |
| gauss_sigma = 1.5 | 42.0 | 63.4 |
| pad = 4 | 43.1 | 63.5 |
| fit = gaussmle | 40.7 | 58.7 |

### Test frames: all ground-truth emitters (primary setting)

Rungs R0 (Picasso defaults), R1E / R1J (net gradient tuned on the development frames for efficiency / Jaccard), R2E / R2J (all options tuned on the development frames for efficiency / Jaccard), each without and with the "+c" offset correction (the constant (dx, dy) estimated on the development frames for that configuration, subtracted from every detection). 95 % intervals of a bootstrap in blocks of 10 consecutive frames (1000 resamples, seed 0; pooled metrics recomputed in every resample). The centred RMSE is the RMSE about the mean error of the matched pairs of the test frames. "RMSE common TP": RMSE over the 130 ground-truth emitters that EVERY one of the ten rows detects, so that no row is flattered by losing its hard emitters.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 22.1 [18.4, 25.9] | 100.0 [100.0, 100.0] | 22.1 [18.4, 25.9] | 23.4 [20.1, 26.8] | 19.7 [16.4, 23.2] | 22.0 [19.7, 24.3] | 18.7 [14.9, 22.5] |
| R0+c | 22.1 [18.4, 25.9] | 100.0 [100.0, 100.0] | 22.1 [18.4, 25.9] | 19.7 [16.6, 23.3] | 19.7 [16.5, 23.2] | 18.2 [16.2, 20.3] | 19.7 [16.0, 23.5] |
| R1E | 53.1 [48.1, 58.1] | 100.0 [100.0, 100.0] | 53.1 [48.1, 58.1] | 33.5 [30.3, 36.3] | 30.4 [27.5, 33.0] | 22.0 [19.7, 24.3] | 42.4 [38.0, 46.8] |
| R1E+c | 53.1 [48.1, 58.1] | 100.0 [100.0, 100.0] | 53.1 [48.1, 58.1] | 30.7 [28.0, 33.3] | 30.4 [27.5, 33.0] | 18.5 [16.5, 20.5] | 44.0 [39.6, 48.4] |
| R1J | 55.9 [51.2, 60.5] | 97.4 [95.7, 98.8] | 56.8 [52.0, 61.5] | 40.2 [35.5, 45.7] | 38.1 [32.7, 44.4] | 22.0 [19.7, 24.3] | 40.4 [34.6, 45.9] |
| R1J+c | 55.9 [51.2, 60.5] | 97.4 [95.7, 98.8] | 56.8 [52.0, 61.5] | 38.4 [33.1, 44.5] | 38.1 [32.7, 44.4] | 18.6 [16.6, 20.6] | 41.6 [35.4, 47.5] |
| R2E | 57.3 [52.9, 61.7] | 100.0 [100.0, 100.0] | 57.3 [52.9, 61.7] | 33.5 [30.0, 37.5] | 31.4 [27.9, 35.2] | 20.5 [18.1, 22.7] | 45.7 [41.0, 50.4] |
| R2E+c | 57.3 [52.9, 61.7] | 100.0 [100.0, 100.0] | 57.3 [52.9, 61.7] | 31.5 [28.3, 35.3] | 31.4 [27.9, 35.3] | 17.2 [14.9, 19.3] | 46.9 [42.2, 51.5] |
| R2J | 65.0 [61.6, 68.4] | 94.2 [92.1, 95.9] | 67.7 [64.5, 71.1] | 52.2 [44.7, 59.8] | 50.2 [42.1, 57.4] | 21.9 [19.5, 24.2] | 37.1 [30.0, 44.5] |
| R2J+c | 65.0 [61.6, 68.4] | 94.2 [92.1, 95.9] | 67.7 [64.5, 71.1] | 50.3 [42.5, 57.7] | 50.2 [42.1, 57.4] | 18.6 [16.4, 20.9] | 38.7 [31.6, 46.2] |

TP / FP / FN: R0 132 / 0 / 465; R1E 317 / 0 / 280; R1J 339 / 9 / 258; R2E 342 / 0 / 255; R2J 404 / 25 / 193. Detected per frame: R0 0.264; R1E 0.634; R1J 0.696; R2E 0.684; R2J 0.858; true per frame: 1.194. Offsets (dx, dy) in nm estimated on the development frames: R0 (10.5, -8.4); R1E (12.4, -8.6); R1J (12.7, -8.7); R2E (10.2, -7.6); R2J (12.7, -8.1).

Paired differences between rows on the same test frames (both evaluated on every bootstrap resample, blocks of 10 frames; 95 % interval; an interval that contains 0 is not distinguishable from frame-sampling noise):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | +31.0 [+27.7, +34.6] | +10.1 [+6.5, +14.1] | +23.7 [+20.4, +26.7] |
| R1J minus R0 | +33.8 [+30.6, +37.0] | +16.8 [+11.2, +23.0] | +21.7 [+17.3, +25.7] |
| R2E minus R0 | +35.2 [+31.7, +38.9] | +10.1 [+6.1, +15.1] | +27.1 [+23.5, +30.7] |
| R2J minus R0 | +42.8 [+39.3, +46.5] | +28.8 [+21.0, +36.6] | +18.4 [+12.7, +24.2] |
| R2E minus R1E | +4.2 [+1.5, +7.1] | -0.0 [-3.7, +4.1] | +3.4 [+0.7, +6.1] |
| R2J minus R1J | +9.0 [+5.7, +12.6] | +12.0 [+5.1, +18.7] | -3.3 [-8.4, +2.3] |
| R1E+c minus R0+c | +31.0 [+27.7, +34.6] | +11.0 [+7.3, +14.7] | +24.3 [+21.0, +27.3] |
| R1J+c minus R0+c | +33.8 [+30.6, +37.0] | +18.7 [+12.7, +25.6] | +21.9 [+17.2, +26.0] |
| R2E+c minus R0+c | +35.2 [+31.7, +38.9] | +11.8 [+7.6, +16.3] | +27.2 [+23.9, +30.7] |
| R2J+c minus R0+c | +42.8 [+39.3, +46.5] | +30.6 [+22.6, +38.6] | +19.1 [+13.3, +24.9] |
| R2E+c minus R1E+c | +4.2 [+1.5, +7.1] | +0.9 [-2.6, +4.8] | +2.9 [+0.4, +5.5] |
| R2J+c minus R1J+c | +9.0 [+5.7, +12.6] | +11.9 [+4.5, +19.0] | -2.9 [-8.6, +3.1] |
| R0+c minus R0 | +0.0 [+0.0, +0.0] | -3.7 [-4.7, -2.6] | +1.0 [+0.7, +1.3] |
| R1E+c minus R1E | +0.0 [+0.0, +0.0] | -2.8 [-3.9, -1.6] | +1.6 [+0.9, +2.2] |
| R1J+c minus R1J | +0.0 [+0.0, +0.0] | -1.8 [-3.1, -0.6] | +1.2 [+0.4, +2.0] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -1.9 [-2.8, -1.0] | +1.2 [+0.6, +1.7] |
| R2J+c minus R2J | +0.0 [+0.0, +0.0] | -1.9 [-3.2, -0.6] | +1.6 [+0.5, +2.6] |

### Test frames: challenge assessment setting (emitters within 450 nm of the border and the dimmest 25 % of the whole set's activations dropped after pairing)

The same ten rows; "RMSE common TP" over the 124 kept emitters that every row detects.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 33.2 [28.2, 38.1] | 100.0 [100.0, 100.0] | 33.2 [28.2, 38.1] | 23.7 [20.3, 27.3] | 20.0 [16.6, 23.6] | 22.3 [20.0, 24.7] | 29.2 [23.9, 34.2] |
| R0+c | 33.2 [28.2, 38.1] | 100.0 [100.0, 100.0] | 33.2 [28.2, 38.1] | 20.0 [16.8, 23.6] | 20.0 [16.7, 23.6] | 18.5 [16.4, 20.6] | 30.3 [25.2, 35.3] |
| R1E | 78.4 [73.1, 83.5] | 100.0 [100.0, 100.0] | 78.4 [73.1, 83.5] | 32.6 [29.7, 35.4] | 29.4 [26.8, 32.1] | 22.3 [20.0, 24.7] | 60.9 [57.1, 64.5] |
| R1E+c | 78.4 [73.1, 83.5] | 100.0 [100.0, 100.0] | 78.4 [73.1, 83.5] | 29.6 [27.1, 32.3] | 29.4 [26.8, 32.1] | 18.7 [16.6, 20.7] | 63.3 [59.6, 66.9] |
| R1J | 80.3 [75.3, 85.2] | 97.8 [96.2, 99.1] | 81.8 [76.9, 86.8] | 38.4 [33.4, 44.5] | 36.5 [30.6, 43.5] | 22.3 [20.0, 24.7] | 56.8 [50.4, 62.4] |
| R1J+c | 80.3 [75.3, 85.2] | 97.8 [96.2, 99.1] | 81.8 [76.9, 86.8] | 36.7 [31.0, 43.6] | 36.5 [30.7, 43.5] | 18.7 [16.7, 20.8] | 58.4 [50.9, 64.5] |
| R2E | 77.6 [72.5, 82.7] | 100.0 [100.0, 100.0] | 77.6 [72.5, 82.7] | 29.5 [27.1, 32.4] | 27.3 [24.8, 30.2] | 20.7 [18.3, 23.1] | 62.9 [58.9, 66.5] |
| R2E+c | 77.6 [72.5, 82.7] | 100.0 [100.0, 100.0] | 77.6 [72.5, 82.7] | 27.4 [25.0, 30.3] | 27.3 [24.8, 30.2] | 17.4 [15.0, 19.7] | 64.6 [60.7, 68.3] |
| R2J | 83.8 [79.4, 88.3] | 94.1 [91.8, 96.1] | 88.4 [84.4, 92.6] | 48.6 [40.5, 56.3] | 46.6 [37.9, 54.7] | 22.1 [19.7, 24.5] | 48.8 [40.9, 57.1] |
| R2J+c | 84.0 [79.7, 88.5] | 94.4 [92.3, 96.2] | 88.4 [84.4, 92.6] | 46.7 [38.1, 54.8] | 46.6 [37.9, 54.7] | 18.8 [16.6, 21.2] | 50.6 [42.2, 59.4] |

Ground-truth emitters counted: 379; detections counted per row: R0 126; R1E 297; R1J 317; R2E 294; R2J 356.

Paired differences (challenge setting):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | +45.1 [+40.6, +50.3] | +8.9 [+5.5, +12.8] | +31.7 [+27.6, +35.8] |
| R1J minus R0 | +47.1 [+42.7, +52.3] | +14.7 [+9.2, +21.8] | +27.7 [+22.2, +32.9] |
| R2E minus R0 | +44.3 [+40.0, +49.5] | +5.8 [+2.7, +9.5] | +33.8 [+30.1, +37.8] |
| R2J minus R0 | +50.5 [+44.7, +56.5] | +24.9 [+16.8, +33.2] | +19.6 [+12.7, +26.6] |
| R2E minus R1E | -0.8 [-2.8, +1.5] | -3.1 [-5.9, -0.2] | +2.1 [+0.1, +4.1] |
| R2J minus R1J | +3.4 [-1.3, +7.9] | +10.1 [+2.8, +18.2] | -8.0 [-14.9, -1.3] |
| R1E+c minus R0+c | +45.1 [+40.6, +50.3] | +9.6 [+6.4, +13.3] | +33.0 [+29.0, +37.2] |
| R1J+c minus R0+c | +47.1 [+42.7, +52.3] | +16.7 [+10.3, +24.6] | +28.0 [+21.8, +33.8] |
| R2E+c minus R0+c | +44.3 [+40.0, +49.5] | +7.4 [+4.1, +11.4] | +34.3 [+30.4, +38.3] |
| R2J+c minus R0+c | +50.7 [+45.0, +56.7] | +26.7 [+18.2, +35.5] | +20.3 [+13.1, +27.7] |
| R2E+c minus R1E+c | -0.8 [-2.8, +1.5] | -2.2 [-4.7, +0.8] | +1.3 [-0.8, +3.1] |
| R2J+c minus R1J+c | +3.6 [-1.0, +8.1] | +10.0 [+1.9, +18.8] | -7.7 [-15.2, -0.1] |
| R0+c minus R0 | +0.0 [+0.0, +0.0] | -3.7 [-4.8, -2.6] | +1.2 [+0.8, +1.5] |
| R1E+c minus R1E | +0.0 [+0.0, +0.0] | -3.0 [-4.1, -1.8] | +2.5 [+1.5, +3.4] |
| R1J+c minus R1J | +0.0 [+0.0, +0.0] | -1.7 [-3.2, -0.4] | +1.5 [+0.4, +2.8] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -2.1 [-3.1, -1.1] | +1.7 [+0.9, +2.4] |
| R2J+c minus R2J | +0.2 [+0.0, +0.8] | -1.9 [-3.1, -0.7] | +1.8 [+0.7, +3.1] |

Pairing precedes the exclusions, so an excluded (dim or border) emitter can take the detection of a kept one. Kept emitters left unmatched although a detection within 250 nm was paired with an excluded emitter, of all kept emitters: R0 0 of 379; R1E 1 of 379; R1J 1 of 379; R2E 1 of 379; R2J 1 of 379.

### Effect of each exclusion on the same detections (Jaccard % / RMSE nm / efficiency, point estimates)

| row | all emitters | border 450 nm dropped | dimmest 25 % dropped | both (challenge setting) |
|---|---|---|---|---|
| R0 | 22.1 / 23.4 / 18.7 | 25.0 / 23.7 / 21.3 | 29.8 / 23.4 / 26.0 | 33.2 / 23.7 / 29.2 |
| R0+c | 22.1 / 19.7 / 19.7 | 25.0 / 20.0 / 22.3 | 29.8 / 19.7 / 27.1 | 33.2 / 20.0 / 30.3 |
| R1E | 53.1 / 33.5 / 42.4 | 60.6 / 33.3 / 48.4 | 69.5 / 32.8 / 55.2 | 78.4 / 32.6 / 60.9 |
| R1E+c | 53.1 / 30.7 / 44.0 | 60.6 / 30.4 / 50.2 | 69.5 / 29.9 / 57.3 | 78.4 / 29.6 / 63.3 |
| R1J | 55.9 / 40.2 / 40.4 | 63.7 / 39.4 / 46.4 | 71.5 / 39.2 / 51.5 | 80.3 / 38.4 / 56.8 |
| R1J+c | 55.9 / 38.4 / 41.6 | 63.7 / 37.8 / 47.6 | 71.5 / 37.3 / 53.0 | 80.3 / 36.7 / 58.4 |
| R2E | 57.3 / 33.5 / 45.7 | 61.2 / 30.7 / 50.5 | 73.4 / 32.5 / 58.0 | 77.6 / 29.5 / 62.9 |
| R2E+c | 57.3 / 31.5 / 46.9 | 61.2 / 28.7 / 51.7 | 73.4 / 30.3 / 59.7 | 77.6 / 27.4 / 64.6 |
| R2J | 65.0 / 52.2 / 37.1 | 68.8 / 50.2 / 40.9 | 79.9 / 50.7 / 45.4 | 83.8 / 48.6 / 48.8 |
| R2J+c | 65.0 / 50.3 / 38.7 | 69.0 / 48.8 / 42.2 | 79.9 / 48.4 / 47.6 | 84.0 / 46.7 / 50.6 |

Picasso never searches the box // 2 pixel rows / columns at the top and left (box // 2 + 1 at the bottom and right: 3 and 4 px for the default 7 px box) when the frame is not padded, so emitters there are missed by construction; the border setting removes most of that zone from the ground truth (450 nm = 4.5 px). Median padding (an R2 option) removes the blind zone.

### Recall by distance to the border, by photons and by |z| (test frames, rungs without the offset correction)

By distance of the emitter to the nearest field edge (all emitters):

| distance to the edge | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-300 nm | 44 | 0.0 | 0.0 | 2.3 | 20.5 | 29.5 |
| 300-450 nm | 48 | 12.5 | 22.9 | 25.0 | 50.0 | 60.4 |
| 450-700 nm | 51 | 3.9 | 25.5 | 29.4 | 23.5 | 47.1 |
| 700-1000 nm | 46 | 19.6 | 43.5 | 43.5 | 45.7 | 50.0 |
| 1000-7000 nm | 408 | 28.2 | 66.9 | 71.3 | 67.6 | 77.2 |

By emitted photons (emitters at least 450 nm from the edge):

| photons | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-50 photons | 25 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 50-100 photons | 38 | 0.0 | 2.6 | 2.6 | 2.6 | 10.5 |
| 100-200 photons | 62 | 0.0 | 12.9 | 24.2 | 22.6 | 37.1 |
| 200-400 photons | 93 | 0.0 | 55.9 | 64.5 | 59.1 | 79.6 |
| 400-700 photons | 124 | 7.3 | 78.2 | 81.5 | 74.2 | 84.7 |
| 700-inf photons | 163 | 71.8 | 90.8 | 91.4 | 90.2 | 95.7 |

By axial distance |z| from the focal plane (emitters at least 450 nm from the edge and at least 200 photons; the PSF broadens with |z|, see `calibration.md`):

| |z| | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-150 nm | 217 | 42.9 | 93.1 | 96.8 | 94.0 | 96.8 |
| 150-300 nm | 82 | 34.1 | 87.8 | 89.0 | 82.9 | 89.0 |
| 300-450 nm | 47 | 10.6 | 48.9 | 57.4 | 46.8 | 91.5 |
| 450-1000 nm | 33 | 0.0 | 0.0 | 0.0 | 0.0 | 27.3 |

### Picasso photons against the emitted photons, and what a one-peak-per-blob detector cannot do

Photons reported by Picasso (QE-corrected) divided by the ground-truth photons of the matched emitter, test frames, median and inter-quartile range (the Gaussian core holds only part of the emitted photons, `calibration.md` section on the halo, so the photon count of a Picasso fit is biased low and more so out of focus):

| row | emitters | n | median ratio | IQR |
|---|---|---|---|---|
| R0 | all matched | 132 | 0.80 | 0.69 .. 0.89 |
| R0 | |z| 0-100 nm | 58 | 0.87 | 0.78 .. 0.96 |
| R0 | |z| 100-200 nm | 60 | 0.76 | 0.68 .. 0.83 |
| R0 | |z| 200-350 nm | 13 | 0.62 | 0.58 .. 0.65 |
| R0 | |z| 350-800 nm | 1 | 0.43 | 0.43 .. 0.43 |
| R2E | all matched | 342 | 0.80 | 0.68 .. 0.94 |
| R2E | |z| 0-100 nm | 124 | 0.89 | 0.80 .. 1.05 |
| R2E | |z| 100-200 nm | 156 | 0.79 | 0.69 .. 0.93 |
| R2E | |z| 200-350 nm | 55 | 0.59 | 0.47 .. 0.71 |
| R2E | |z| 350-800 nm | 7 | 0.47 | 0.19 .. 0.55 |

Structural ceiling of a detector that reports one peak per blob: 3.4 % of the 597 ground-truth rows of the test frames have another emitter of the same frame within 250 nm and 3.7 % within 400 nm; on the development frames, keeping EVERY candidate of the chosen structure (R2E: 63.6 per frame, R2J: 35.2 per frame) the recall is at most 85.6 % (R2E structure) / 83.1 % (R2J structure) (pairing with the maximum number of pairs). The development and test frames show the same three microtubules, and `positions.csv` holds the exact emitter sites of the fluorophores (no frame information), which a method that uses the structure as a prior could exploit; Picasso uses neither.

## MT0.N2.HD

Development frames 0-99 (100 frames, 1294 ground-truth activations), test frames 100-499 (400 frames, 5046 activations). Frame indices are 0-based stack indices. Ground-truth activations of the whole set: 31612; photon threshold of the challenge setting (25 % quantile of the whole set's ground truth): 198 photons. Persistence: a fluorophore stays on for 2.26 consecutive frames on average (13981 runs, 31 % of them 3 frames or longer), so neighbouring frames are correlated, which the block bootstrap accounts for.

### R1: tuning curve of the minimum net gradient on the development frames (20 values)

| min net gradient | Jaccard % | precision % | recall % | RMSE nm | efficiency | detections / frame |
|---|---|---|---|---|---|---|
| 100 | 13.6 | 15.6 | 52.1 | 79.8 | -17.6 | 43.22 |
| 150 | 14.7 | 17.0 | 52.0 | 79.8 | -16.8 | 39.61 |
| 200 | 15.9 | 18.6 | 51.5 | 77.8 | -14.6 | 35.77 |
| 250 | 17.3 | 20.7 | 51.3 | 77.2 | -13.1 | 32.03 |
| 300 | 19.1 | 23.5 | 50.7 | 75.3 | -10.5 | 27.90 |
| 400 | 23.4 | 30.7 | 49.8 | 73.1 | -5.9 | 21.00 |
| 500 | 28.2 | 39.9 | 48.8 | 70.5 | -0.7 | 15.82 |
| 600 | 33.2 | 52.2 | 47.6 | 67.9 | 4.7 | 11.80 |
| 800 | 41.9 | 83.5 | 45.7 | 65.0 | 12.8 | 7.08 |
| **1000** | 42.6 | 95.3 | 43.5 | 59.3 | 17.4 | 5.91 |
| 1250 | 40.9 | 99.3 | 41.0 | 53.8 | 20.1 | 5.35 |
| 1500 | 38.9 | 100.0 | 38.9 | 50.9 | 20.4 | 5.03 |
| **2000** | 35.7 | 100.0 | 35.7 | 46.4 | 20.7 | 4.62 |
| 2500 | 32.6 | 100.0 | 32.6 | 43.7 | 19.7 | 4.22 |
| 3000 | 29.8 | 100.0 | 29.8 | 42.6 | 17.9 | 3.86 |
| 3500 | 27.2 | 100.0 | 27.2 | 42.7 | 15.6 | 3.52 |
| 4000 | 25.0 | 100.0 | 25.0 | 42.7 | 13.7 | 3.24 |
| 5000 | 19.7 | 100.0 | 19.7 | 40.6 | 10.0 | 2.55 |
| 6000 | 15.5 | 100.0 | 15.5 | 41.1 | 6.0 | 2.00 |
| 8000 | 7.3 | 100.0 | 7.3 | 43.6 | -2.5 | 0.94 |

Bold: the value maximising the efficiency (2000; R1E) and the value maximising the Jaccard index (1000; R1J). The curve is computed from one run at the smallest grid value, filtered at each threshold (identification keeps the maxima with net gradient above the threshold and every spot is fitted independently); this was verified to reproduce direct runs at 2000 and 8000 on the first 60 development frames.

How well defined are these two operating points? On the development frames 3 of the 20 thresholds (from 1250 to 2000) have an efficiency within 1 point of the best (20.7): the efficiency is a plateau and its argmax is arbitrary within it. In 1000 block-bootstrap resamples of the development frames (blocks of 10 frames) the efficiency-optimal threshold is 2000 in 56 %, 1500 in 31 %, 1250 in 11 %; the Jaccard-optimal threshold is 1000 in 98 %, 800 in 2 %, 100 in 0 %. A falling RMSE with a rising threshold is largely survivor bias (the hard emitters are lost): compare the column "RMSE common TP" below, where the rows are scored on the same emitters.

### R2: search over all options on the development frames

Declared grid (identical for the sets of one noise class; at the top of `experiments/exp_smlm_picasso.py`): identification box 5/7/9 px, fit box 5/7/9 px, Gaussian filter sigma off/1/1.5 px, median padding 0/4 px, fit gausslq/gaussmle: **108 structural settings**, each run once at net gradient 30 and then post-filtered by a coordinate ascent (at most 3 rounds over, in this order: net gradient (28 values), minimum photons (9), minimum width (5), maximum width (7)) once for the efficiency and once for the Jaccard index. Development frames only.

**Configurations tried on the development frames: 108 structural settings run, 13413 distinct (structure, post-filter) configurations scored** (R1: 20 net-gradient values of the default structure, all contained in them); 516512 candidate localizations fitted in total. The coordinate ascent converged (a round without change) within the 3 rounds in 212 of 216 runs.

Chosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); the offset is the mean signed error (detection - truth) of the matched pairs, the correction of the "+c" rows:

| rung | configuration | dev J % | dev precision % | dev recall % | dev RMSE nm | dev E | offset (dx, dy) nm | dev TP |
|---|---|---|---|---|---|---|---|---|
| R0 | box 7 px, raw frames, no padding, gausslq, net gradient > 5000 | 19.7 | 100.0 | 19.7 | 40.6 | 10.0 | (8.8, -4.7) | 255 |
| R1E | box 7 px, raw frames, no padding, gausslq, net gradient > 2000 | 35.7 | 100.0 | 35.7 | 46.4 | 20.7 | (10.5, -7.8) | 462 |
| R1J | box 7 px, raw frames, no padding, gausslq, net gradient > 1000 | 42.6 | 95.3 | 43.5 | 59.3 | 17.4 | (13.6, -9.4) | 563 |
| R2E | box 5 px, Gaussian filter 1 px, median pad 4 px, gaussmle, net gradient > 1000, photons > 80 | 38.9 | 99.8 | 38.9 | 41.2 | 26.3 | (9.3, -6.0) | 504 |
| R2J | ident 5 / fit 7 px, Gaussian filter 1 px, median pad 4 px, gausslq, net gradient > 500 | 48.7 | 95.4 | 49.8 | 63.4 | 18.5 | (10.8, -10.5) | 645 |

Ablation of R2E on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=5|fit_box=5|gauss_sigma=1.0|pad=4|fit=gaussmle`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 38.9 | 26.3 | 41.2 |
| identification box -> 7 px | 37.2 | 24.4 | 42.1 |
| fit box -> the identification box | (already the default) |  |  |
| Gaussian filter -> off | 39.4 | 25.8 | 42.8 |
| median padding -> 0 | 36.9 | 25.1 | 40.3 |
| fit -> gausslq | 38.9 | 26.1 | 41.7 |
| photon filter -> off (the other post-filters kept, not re-tuned) | 39.0 | 26.3 | 41.4 |
| width window -> off (the other post-filters kept, not re-tuned) | 38.9 | 26.3 | 41.2 |
| every structural option -> default (post-filters re-tuned) | 38.7 | 22.9 | 46.7 |
| R1: net gradient only, 20 values | 35.7 | 20.7 | 46.4 |
| R0: Picasso defaults | 19.7 | 10.0 | 40.6 |

Ablation of R2J on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `box=5|fit_box=7|gauss_sigma=1.0|pad=4|fit=gausslq`):

| configuration | dev J % | dev E | dev RMSE nm |
|---|---|---|---|
| chosen configuration | 48.7 | 18.5 | 63.4 |
| identification box -> 7 px | 47.2 | 18.8 | 61.7 |
| fit box -> the identification box | 48.1 | 20.5 | 60.3 |
| Gaussian filter -> off | 46.8 | 17.2 | 63.5 |
| median padding -> 0 | 45.5 | 17.1 | 62.4 |
| fit -> gausslq | (already the default) |  |  |
| photon filter -> off (the other post-filters kept, not re-tuned) | 48.7 | 18.5 | 63.4 |
| width window -> off (the other post-filters kept, not re-tuned) | 48.7 | 18.5 | 63.4 |
| every structural option -> default (post-filters re-tuned) | 42.6 | 17.4 | 59.3 |
| R1: net gradient only, 20 values | 42.6 | 17.4 | 59.3 |
| R0: Picasso defaults | 19.7 | 10.0 | 40.6 |

Each option alone on top of R1 (the default structure with only the net gradient tuned, no other filter; the net gradient is tuned for the objective of the column), development frames:

| setting | dev E (net gradient tuned for E) | dev J % (net gradient tuned for J) |
|---|---|---|
| default structure, net gradient tuned (R1 on the finer grid) | 20.7 | 42.6 |
| box = 5 | 20.0 | 43.2 |
| box = 9 | 18.5 | 39.5 |
| fit_box = 5 | 21.3 | 42.5 |
| fit_box = 9 | 14.8 | 41.6 |
| gauss_sigma = 1.0 | 20.9 | 42.7 |
| gauss_sigma = 1.5 | 18.4 | 41.4 |
| pad = 4 | 22.9 | 46.9 |
| fit = gaussmle | 18.8 | 42.2 |

### Test frames: all ground-truth emitters (primary setting)

Rungs R0 (Picasso defaults), R1E / R1J (net gradient tuned on the development frames for efficiency / Jaccard), R2E / R2J (all options tuned on the development frames for efficiency / Jaccard), each without and with the "+c" offset correction (the constant (dx, dy) estimated on the development frames for that configuration, subtracted from every detection). 95 % intervals of a bootstrap in blocks of 10 consecutive frames (1000 resamples, seed 0; pooled metrics recomputed in every resample). The centred RMSE is the RMSE about the mean error of the matched pairs of the test frames. "RMSE common TP": RMSE over the 823 ground-truth emitters that EVERY one of the ten rows detects, so that no row is flattered by losing its hard emitters.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 18.1 [17.2, 19.1] | 100.0 [100.0, 100.0] | 18.1 [17.2, 19.1] | 45.7 [43.6, 48.1] | 44.8 [42.6, 47.1] | 42.4 [40.0, 44.8] | 6.2 [4.9, 7.5] |
| R0+c | 18.1 [17.2, 19.1] | 100.0 [100.0, 100.0] | 18.1 [17.2, 19.1] | 44.5 [42.4, 46.8] | 44.5 [42.4, 46.8] | 40.9 [38.5, 43.2] | 6.8 [5.6, 8.0] |
| R1E | 35.1 [34.0, 36.5] | 99.9 [99.8, 100.0] | 35.1 [34.0, 36.5] | 49.4 [47.4, 51.4] | 47.6 [45.6, 49.7] | 42.4 [40.0, 44.8] | 18.5 [17.0, 20.1] |
| R1E+c | 35.1 [34.0, 36.5] | 99.9 [99.8, 100.0] | 35.1 [34.0, 36.5] | 47.2 [45.2, 49.3] | 47.1 [45.1, 49.1] | 40.8 [38.4, 43.2] | 19.8 [18.2, 21.3] |
| R1J | 42.9 [41.7, 44.3] | 95.7 [94.8, 96.6] | 43.7 [42.5, 45.2] | 60.1 [58.1, 61.9] | 58.3 [56.4, 60.2] | 42.4 [40.0, 44.8] | 17.1 [15.5, 18.8] |
| R1J+c | 42.9 [41.7, 44.3] | 95.8 [94.9, 96.6] | 43.8 [42.5, 45.3] | 58.2 [56.1, 60.1] | 58.1 [56.1, 60.0] | 41.1 [38.6, 43.5] | 18.5 [16.9, 20.2] |
| R2E | 38.7 [37.4, 39.9] | 99.8 [99.7, 100.0] | 38.7 [37.5, 40.0] | 45.7 [44.0, 47.3] | 44.0 [42.3, 45.7] | 38.3 [36.1, 40.6] | 23.6 [22.0, 25.2] |
| R2E+c | 38.7 [37.4, 39.9] | 99.8 [99.7, 100.0] | 38.7 [37.5, 40.0] | 43.9 [42.2, 45.6] | 43.8 [42.1, 45.5] | 36.7 [34.5, 38.9] | 24.6 [23.0, 26.3] |
| R2J | 48.4 [47.1, 49.8] | 94.5 [93.7, 95.4] | 49.8 [48.4, 51.3] | 60.5 [58.4, 62.6] | 58.7 [56.6, 60.9] | 43.3 [40.7, 46.0] | 20.5 [18.7, 22.4] |
| R2J+c | 48.5 [47.2, 49.9] | 94.7 [93.8, 95.5] | 49.9 [48.5, 51.4] | 59.1 [56.9, 61.2] | 59.0 [56.9, 61.1] | 42.0 [39.3, 44.8] | 21.7 [19.6, 23.6] |

TP / FP / FN: R0 915 / 0 / 4131; R1E 1773 / 1 / 3273; R1J 2207 / 98 / 2839; R2E 1955 / 3 / 3091; R2J 2514 / 145 / 2532. Detected per frame: R0 2.288; R1E 4.435; R1J 5.763; R2E 4.895; R2J 6.647; true per frame: 12.615. Offsets (dx, dy) in nm estimated on the development frames: R0 (8.8, -4.7); R1E (10.5, -7.8); R1J (13.6, -9.4); R2E (9.3, -6.0); R2J (10.8, -10.5).

Paired differences between rows on the same test frames (both evaluated on every bootstrap resample, blocks of 10 frames; 95 % interval; an interval that contains 0 is not distinguishable from frame-sampling noise):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | +17.0 [+15.8, +18.2] | +3.7 [+1.4, +6.1] | +12.2 [+10.9, +13.7] |
| R1J minus R0 | +24.8 [+23.5, +26.2] | +14.3 [+12.2, +16.4] | +10.9 [+9.3, +12.4] |
| R2E minus R0 | +20.6 [+19.2, +22.0] | -0.1 [-2.5, +2.4] | +17.3 [+15.8, +18.9] |
| R2J minus R0 | +30.3 [+28.7, +32.0] | +14.8 [+12.6, +17.2] | +14.3 [+12.5, +16.2] |
| R2E minus R1E | +3.6 [+2.8, +4.4] | -3.7 [-5.4, -2.1] | +5.1 [+4.0, +6.2] |
| R2J minus R1J | +5.5 [+4.6, +6.5] | +0.4 [-1.4, +2.4] | +3.4 [+1.9, +4.8] |
| R1E+c minus R0+c | +17.0 [+15.8, +18.2] | +2.7 [+0.4, +5.1] | +13.0 [+11.6, +14.4] |
| R1J+c minus R0+c | +24.8 [+23.5, +26.2] | +13.6 [+11.4, +15.8] | +11.7 [+10.2, +13.3] |
| R2E+c minus R0+c | +20.6 [+19.2, +22.0] | -0.7 [-3.2, +1.7] | +17.8 [+16.3, +19.4] |
| R2J+c minus R0+c | +30.4 [+28.8, +32.1] | +14.5 [+12.3, +16.8] | +14.9 [+13.0, +16.8] |
| R2E+c minus R1E+c | +3.6 [+2.8, +4.4] | -3.4 [-5.0, -1.7] | +4.9 [+3.7, +6.0] |
| R2J+c minus R1J+c | +5.6 [+4.7, +6.6] | +0.9 [-1.0, +2.9] | +3.1 [+1.7, +4.7] |
| R0+c minus R0 | +0.0 [+0.0, +0.0] | -1.2 [-1.7, -0.7] | +0.6 [+0.4, +0.8] |
| R1E+c minus R1E | +0.0 [+0.0, +0.0] | -2.2 [-2.6, -1.8] | +1.3 [+1.1, +1.5] |
| R1J+c minus R1J | +0.0 [-0.1, +0.1] | -1.9 [-2.7, -1.1] | +1.4 [+0.9, +1.9] |
| R2E+c minus R2E | +0.0 [+0.0, +0.0] | -1.8 [-2.1, -1.5] | +1.1 [+0.9, +1.2] |
| R2J+c minus R2J | +0.1 [+0.0, +0.2] | -1.5 [-2.1, -0.7] | +1.2 [+0.6, +1.7] |

### Test frames: challenge assessment setting (emitters within 450 nm of the border and the dimmest 25 % of the whole set's activations dropped after pairing)

The same ten rows; "RMSE common TP" over the 806 kept emitters that every row detects.

| row | J % | precision % | recall % | RMSE nm | centred RMSE nm | RMSE common TP nm | efficiency |
|---|---|---|---|---|---|---|---|
| R0 | 27.0 [25.8, 28.2] | 100.0 [100.0, 100.0] | 27.0 [25.8, 28.2] | 45.7 [43.5, 48.0] | 44.6 [42.4, 47.0] | 42.4 [40.0, 44.7] | 13.9 [12.5, 15.3] |
| R0+c | 27.0 [25.7, 28.2] | 100.0 [100.0, 100.0] | 27.0 [25.7, 28.2] | 44.4 [42.3, 46.8] | 44.3 [42.2, 46.7] | 40.8 [38.4, 43.1] | 14.5 [13.1, 15.9] |
| R1E | 51.7 [50.1, 53.4] | 99.9 [99.8, 100.0] | 51.8 [50.1, 53.4] | 49.0 [47.1, 51.2] | 47.2 [45.2, 49.3] | 42.4 [40.0, 44.7] | 31.2 [29.2, 33.2] |
| R1E+c | 51.6 [50.0, 53.3] | 99.9 [99.8, 100.0] | 51.6 [50.0, 53.3] | 47.0 [45.0, 49.1] | 46.9 [44.9, 49.0] | 40.8 [38.4, 43.1] | 32.5 [30.6, 34.5] |
| R1J | 59.7 [58.1, 61.2] | 95.8 [94.8, 96.7] | 61.3 [59.7, 63.0] | 57.4 [55.2, 59.6] | 55.5 [53.2, 57.5] | 42.4 [40.0, 44.7] | 29.8 [27.8, 31.9] |
| R1J+c | 59.6 [58.0, 61.2] | 95.8 [94.9, 96.8] | 61.2 [59.6, 63.0] | 55.5 [53.1, 57.8] | 55.4 [53.1, 57.7] | 41.0 [38.6, 43.4] | 31.4 [29.4, 33.5] |
| R2E | 52.8 [51.1, 54.4] | 99.9 [99.8, 100.0] | 52.8 [51.2, 54.4] | 45.0 [43.2, 46.8] | 43.2 [41.2, 45.0] | 38.3 [36.0, 40.6] | 34.8 [32.7, 36.7] |
| R2E+c | 52.9 [51.2, 54.6] | 99.9 [99.8, 100.0] | 52.9 [51.2, 54.6] | 43.2 [41.3, 45.0] | 43.0 [41.1, 44.8] | 36.7 [34.4, 38.8] | 36.1 [34.1, 38.0] |
| R2J | 61.1 [59.6, 62.8] | 95.0 [94.2, 95.9] | 63.2 [61.5, 64.9] | 57.1 [54.8, 59.3] | 55.1 [52.7, 57.4] | 43.3 [40.7, 46.1] | 31.0 [28.9, 33.0] |
| R2J+c | 61.2 [59.6, 62.9] | 95.1 [94.2, 96.0] | 63.2 [61.5, 64.9] | 55.4 [53.0, 57.6] | 55.3 [52.9, 57.5] | 42.0 [39.3, 44.9] | 32.4 [30.2, 34.6] |

Ground-truth emitters counted: 3314; detections counted per row: R0 896; R1E 1716; R1J 2120; R2E 1752; R2J 2202.

Paired differences (challenge setting):

| difference | Jaccard (points) | RMSE (nm) | efficiency (points) |
|---|---|---|---|
| R1E minus R0 | +24.7 [+22.9, +26.6] | +3.4 [+1.0, +5.7] | +17.3 [+15.4, +19.0] |
| R1J minus R0 | +32.6 [+30.8, +34.6] | +11.8 [+9.6, +14.0] | +15.9 [+13.9, +17.8] |
| R2E minus R0 | +25.8 [+23.9, +27.6] | -0.6 [-3.1, +1.6] | +20.9 [+19.0, +22.7] |
| R2J minus R0 | +34.1 [+32.1, +36.2] | +11.4 [+9.0, +13.9] | +17.0 [+14.9, +19.0] |
| R2E minus R1E | +1.1 [+0.3, +1.8] | -4.0 [-5.6, -2.5] | +3.6 [+2.4, +4.7] |
| R2J minus R1J | +1.5 [+0.6, +2.5] | -0.4 [-2.1, +1.5] | +1.1 [-0.3, +2.4] |
| R1E+c minus R0+c | +24.6 [+22.9, +26.5] | +2.6 [+0.3, +4.9] | +18.0 [+16.2, +19.8] |
| R1J+c minus R0+c | +32.6 [+30.8, +34.6] | +11.1 [+8.7, +13.4] | +16.9 [+15.0, +18.7] |
| R2E+c minus R0+c | +25.9 [+24.0, +27.8] | -1.3 [-3.7, +1.0] | +21.6 [+19.8, +23.4] |
| R2J+c minus R0+c | +34.2 [+32.2, +36.3] | +10.9 [+8.6, +13.3] | +17.9 [+15.7, +19.8] |
| R2E+c minus R1E+c | +1.3 [+0.6, +2.1] | -3.8 [-5.3, -2.3] | +3.6 [+2.5, +4.6] |
| R2J+c minus R1J+c | +1.6 [+0.7, +2.5] | -0.1 [-1.8, +1.6] | +1.0 [-0.3, +2.3] |
| R0+c minus R0 | -0.1 [-0.2, +0.0] | -1.2 [-1.7, -0.8] | +0.6 [+0.4, +0.9] |
| R1E+c minus R1E | -0.2 [-0.3, +0.0] | -2.1 [-2.5, -1.6] | +1.3 [+1.0, +1.7] |
| R1J+c minus R1J | -0.1 [-0.3, +0.2] | -2.0 [-2.9, -0.9] | +1.6 [+0.8, +2.3] |
| R2E+c minus R2E | +0.1 [-0.0, +0.2] | -1.9 [-2.3, -1.4] | +1.3 [+1.1, +1.6] |
| R2J+c minus R2J | +0.1 [-0.1, +0.3] | -1.7 [-2.3, -1.0] | +1.5 [+0.9, +1.9] |

Pairing precedes the exclusions, so an excluded (dim or border) emitter can take the detection of a kept one. Kept emitters left unmatched although a detection within 250 nm was paired with an excluded emitter, of all kept emitters: R0 19 of 3314; R1E 42 of 3314; R1J 53 of 3314; R2E 48 of 3314; R2J 60 of 3314.

### Effect of each exclusion on the same detections (Jaccard % / RMSE nm / efficiency, point estimates)

| row | all emitters | border 450 nm dropped | dimmest 25 % dropped | both (challenge setting) |
|---|---|---|---|---|
| R0 | 18.1 / 45.7 / 6.2 | 20.9 / 45.7 / 8.7 | 23.8 / 45.7 / 11.1 | 27.0 / 45.7 / 13.9 |
| R0+c | 18.1 / 44.5 / 6.8 | 20.9 / 44.4 / 9.3 | 23.7 / 44.5 / 11.7 | 27.0 / 44.4 / 14.5 |
| R1E | 35.1 / 49.4 / 18.5 | 40.3 / 49.1 / 22.7 | 45.7 / 49.2 / 26.8 | 51.7 / 49.0 / 31.2 |
| R1E+c | 35.1 / 47.2 / 19.8 | 40.3 / 47.0 / 24.0 | 45.6 / 47.2 / 28.0 | 51.6 / 47.0 / 32.5 |
| R1J | 42.9 / 60.1 / 17.1 | 48.7 / 59.0 / 21.8 | 53.3 / 58.5 / 25.1 | 59.7 / 57.4 / 29.8 |
| R1J+c | 42.9 / 58.2 / 18.5 | 48.7 / 57.2 / 23.2 | 53.2 / 56.5 / 26.7 | 59.6 / 55.5 / 31.4 |
| R2E | 38.7 / 45.7 / 23.6 | 41.3 / 45.1 / 26.0 | 50.2 / 45.5 / 32.6 | 52.8 / 45.0 / 34.8 |
| R2E+c | 38.7 / 43.9 / 24.6 | 41.3 / 43.2 / 27.1 | 50.3 / 43.8 / 33.8 | 52.9 / 43.2 / 36.1 |
| R2J | 48.4 / 60.5 / 20.5 | 50.4 / 58.8 / 23.0 | 59.2 / 58.8 / 28.4 | 61.1 / 57.1 / 31.0 |
| R2J+c | 48.5 / 59.1 / 21.7 | 50.4 / 57.1 / 24.4 | 59.3 / 57.5 / 29.6 | 61.2 / 55.4 / 32.4 |

Picasso never searches the box // 2 pixel rows / columns at the top and left (box // 2 + 1 at the bottom and right: 3 and 4 px for the default 7 px box) when the frame is not padded, so emitters there are missed by construction; the border setting removes most of that zone from the ground truth (450 nm = 4.5 px). Median padding (an R2 option) removes the blind zone.

### Recall by distance to the border, by photons and by |z| (test frames, rungs without the offset correction)

By distance of the emitter to the nearest field edge (all emitters):

| distance to the edge | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-300 nm | 440 | 0.0 | 0.2 | 1.8 | 21.8 | 36.6 |
| 300-450 nm | 287 | 3.8 | 10.5 | 18.8 | 26.5 | 42.9 |
| 450-700 nm | 370 | 7.8 | 22.2 | 33.5 | 22.7 | 39.2 |
| 700-1000 nm | 412 | 10.7 | 31.6 | 44.2 | 32.3 | 44.9 |
| 1000-7000 nm | 3537 | 23.5 | 43.3 | 52.0 | 44.3 | 53.7 |

By emitted photons (emitters at least 450 nm from the edge):

| photons | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-50 photons | 252 | 0.4 | 0.8 | 2.4 | 1.2 | 4.0 |
| 50-100 photons | 275 | 0.7 | 2.9 | 6.9 | 2.9 | 9.1 |
| 100-200 photons | 484 | 1.2 | 3.7 | 18.6 | 4.5 | 21.7 |
| 200-400 photons | 873 | 1.8 | 19.4 | 40.0 | 21.3 | 45.4 |
| 400-700 photons | 939 | 10.8 | 48.3 | 57.2 | 48.9 | 57.8 |
| 700-inf photons | 1496 | 52.0 | 72.9 | 76.5 | 73.9 | 76.9 |

By axial distance |z| from the focal plane (emitters at least 450 nm from the edge and at least 198 photons; the PSF broadens with |z|, see `calibration.md`):

| |z| | emitters | recall % R0 | recall % R1E | recall % R1J | recall % R2E | recall % R2J |
|---|---|---|---|---|---|---|
| 0-150 nm | 1877 | 39.7 | 62.7 | 68.4 | 63.3 | 68.9 |
| 150-300 nm | 800 | 17.2 | 51.1 | 63.4 | 54.8 | 66.8 |
| 300-450 nm | 402 | 3.0 | 29.6 | 47.5 | 29.1 | 53.0 |
| 450-1000 nm | 235 | 0.0 | 4.3 | 21.3 | 3.4 | 22.6 |

### Picasso photons against the emitted photons, and what a one-peak-per-blob detector cannot do

Photons reported by Picasso (QE-corrected) divided by the ground-truth photons of the matched emitter, test frames, median and inter-quartile range (the Gaussian core holds only part of the emitted photons, `calibration.md` section on the halo, so the photon count of a Picasso fit is biased low and more so out of focus):

| row | emitters | n | median ratio | IQR |
|---|---|---|---|---|
| R0 | all matched | 915 | 0.92 | 0.75 .. 1.35 |
| R0 | |z| 0-100 nm | 398 | 0.95 | 0.81 .. 1.33 |
| R0 | |z| 100-200 nm | 435 | 0.89 | 0.73 .. 1.39 |
| R0 | |z| 200-350 nm | 79 | 0.82 | 0.68 .. 1.24 |
| R0 | |z| 350-800 nm | 3 | 1.10 | 0.84 .. 2.06 |
| R2E | all matched | 1955 | 1.04 | 0.64 .. 2.02 |
| R2E | |z| 0-100 nm | 635 | 1.16 | 0.80 .. 1.87 |
| R2E | |z| 100-200 nm | 896 | 0.99 | 0.63 .. 1.91 |
| R2E | |z| 200-350 nm | 366 | 0.80 | 0.41 .. 2.73 |
| R2E | |z| 350-800 nm | 58 | 0.88 | 0.26 .. 58.10 |

Structural ceiling of a detector that reports one peak per blob: 35.0 % of the 5046 ground-truth rows of the test frames have another emitter of the same frame within 250 nm and 51.6 % within 400 nm; on the development frames, keeping EVERY candidate of the chosen structure (R2E: 86.9 per frame, R2J: 87.9 per frame) the recall is at most 59.0 % (R2E structure) / 59.2 % (R2J structure) (pairing with the maximum number of pairs). The development and test frames show the same three microtubules, and `positions.csv` holds the exact emitter sites of the fluorophores (no frame information), which a method that uses the structure as a prior could exploit; Picasso uses neither.

Detections of every row on the development and test frames: `results/smlm/picasso_ladder_raw.npz` (sets: MT0.N1.LD, MT0.N1.HD, MT0.N2.LD, MT0.N2.HD; sha256 `04ab6c77f1ce8c0c33994826d53307d719e6715791edf3dd1569c65b8cc9ce00`; format in its README entry).

