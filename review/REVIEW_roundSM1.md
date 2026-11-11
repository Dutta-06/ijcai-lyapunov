# Referee report SM1: the SMLM harness, the calibration and the Picasso baseline (WP4)

(Archived verbatim from the referee subagent's final message; the referee worked on commit 40a7252 in its own copy and did not modify the repository.)

## Scope

Commit 40a7252 was exported with `git archive` to W/repo, W = scratchpad/refSM1. Nothing in the repository was touched. The repository venv ran the harness. Picasso 0.11.3 ran in venv_domain with `-I`, one thread, `nice -n 10`, at most two heavy processes at a time. Data: /home/user/smlm_data, read only. My scripts, logs and cached arrays are in W/work and W/logs.

I read in full: smlm_common.py, exp_smlm_calibration.py, exp_smlm_picasso.py, test_smlm_common.py, calibration.md, picasso_baseline.md, fetch.py, checksums.json, then CLAIMS.md section 6 last. The new method was not looked at.

The challenge definitions are [memory] only: Hungarian pairing at 250 nm lateral, J = TP/(TP+FP+FN), RMSE, E = 100 − sqrt((100−J)² + α²RMSE²) with α = 1/nm, the 450 nm border rule and the 25 % dimmest rule. I could not verify them (see SM1-7).

Protocol of my own tuning: every tuning step used development frames only (LD 0–999, HD 0–99). Test frames (LD 1000–1499, HD 100–499) were only scored. I scored about 20 dev-selected configurations on the LD test frames and about 12 on the HD test frames, and I report all of them. The decision to extend the search after the first LD test result (77.4 / 79.4 %) was prompted by that result. Every selection inside a sweep was by dev score alone.

## Summary

Verdicts per stated finding:

- **S1 (conventions): confirmed, with one qualification.** The frame index is 1-based, x is the column axis, and the corner convention holds. The (+11, −9) nm offset is real, but it is not a constant. It depends on the estimator (the centre of mass of the whole PSF has an offset of ≈ 0) and it grows with |z|. It is a property of the simulated, asymmetric PSF seen through a symmetric fit, not a convention error (SM1-4).
- **S2 (calibration): mixed.**
  - Confirmed: the PSF is wider than 0.93 px (in-focus σ 1.27–1.31 px, stable across fit windows); a ≈ 2.22–2.23; 1.25 emitters per frame with 28 % empty frames, 12.6 in HD; photons median 4259 on the dev frames.
  - Refuted: "the background is not spatially flat (illumination profile)" and "frame-level common-mode fluctuation". Both are the PSF halo of the emitters (SM1-3).
  - Qualified: 100.3 photons (true constant background 99.92); the defocus width 2.36 px depends on the fit window (1.86–2.6 px) (SM1-5).
- **S3 (Picasso numbers): reproduced to the digit, but the baseline is far from the best Picasso.** With dev-only tuning I reach on test:
  - LD J = 82.5 % [79.4, 85.6] (reported best 73.4 %), so J ≥ 0.8 is reached.
  - HD J = 56.4 % [55.2, 57.6] (reported best 49.0 %), so J ≥ 0.6 is not reached.
  - LD E = 64.0 [60.6, 67.6] (reported 55.2 / 56.6).
  - HD E = 35.3 [33.7, 36.9] (reported 28.0 / 28.1).
  - Picasso was not given a fair chance (SM1-1, SM1-2).
- **S4 (recall by distance to the border and |z|): reproduced.** Both effects are artefacts of the baseline settings. With padding and a Gaussian pre-filter the recall at 0–300 nm from the border goes from 6.8 % to 50 %, and beyond 450 nm in |z| from 27 % to 73 % (SM1-1).
- **S5 (runner vs Picasso): confirmed.** The CLI gives 158 localizations identical to the runner. Picasso photons are (DN − baseline)·sensitivity/gain with QE not removed. The LS fit is point-sampled (code). The deprecated integrated-Gaussian MLE gives the same J and RMSE.

Numbered findings:

| id | sev. | one line |
|---|---|---|
| SM1-1 | major | Baseline = Picasso with one tuned parameter. Dev-tuned box / Gaussian filter / padding gives +9 J points (LD) and +8.6 (HD). |
| SM1-2 | major | Baseline E is understated by 4–9 points (constant offset, width filter). The "tuned (efficiency)" argmax is arbitrary. The RMSE drop is survivor bias. |
| SM1-3 | major | calibration.md (b) misattributes the PSF halo (≈ 30–44 % of the signal beyond 5 px) to illumination and common mode, and omits the halo. |
| SM1-4 | minor | The (+11, −9) nm offset is estimator- and z-dependent, not constant. |
| SM1-5 | minor | Defocus widths depend on the fit window. |
| SM1-6 | minor | The iid frame bootstrap undercovers (82 % empirical coverage on real data). Paired intervals are iid only. |
| SM1-7 | minor | Challenge rules are unreferenced [memory]. Their implementation has edge-case behaviours (stealing by excluded emitters, NaN). Absolute E depends on the RMSE definition. |
| SM1-8 | minor | fetch.py: a corrupted same-size local file gives a wrong diagnosis, and `--accept-new` then records the corrupted digest. |
| SM1-9 | minor | Omissions: Picasso photon bias, HD structural ceiling, equal tuning effort, shared structure between dev and test. |
| SM1-10 | minor | The two-point noise determination (N2.LD) is available but not used. |
| SM1-11 | cosmetic | Small inconsistencies. |

No error was found in the matching or the metric code. The problems are what the baseline omits and two false sentences in calibration.md.

## What I checked and how

### Item 1: is Picasso given a fair chance?

**Pipeline.** I built W/work/pic.py on `picasso.localize.identify` and `fit`. It allows an identification box different from the fit box, Picasso's own "Gaussian Filter Sigma", padding of the frames, all `gausslq` / `gaussmle` variants, the deprecated integrated-Gaussian MLE, and post-filters on net gradient, photons, sx, sy and bg.

- It reproduces the harness exactly: LD default J 0.7169, 428/0/169.
- float32 and uint16 input give identical results.
- Filtering one low-threshold run equals a direct run [own].

**Search.**
- LD: 112 structural configurations, then 38 with median padding, then about 25 fit and padding variants.
- HD: 128, then 42, then about 25.
- Each configuration got a coordinate-ascent post-filter on dev.
- The camera parameters change nothing in detection or LS positions (only the photon units). The temporal-median window gives no gain (dev J 0.734–0.743 vs 0.738). A gradient-boosting classifier is worse than the rule filter (test J 74.5 vs 79.4 / 82.5).

**Results on the test frames.**

| LD (tuned on dev 0–999) | J % [95 % CI] | P / R % | RMSE nm | E |
|---|---|---|---|---|
| reported default (ng 5000) | 71.7 [68.2, 75.0] | 100 / 71.7 | 34.7 | 55.2 |
| reported "tuned J" (4000) | 73.4 | 98.9 / 74.0 | 41.3 | 50.8 |
| reported "tuned E" (8000) | 64.8 | 100 / 64.8 | 25.3 | 56.6 |
| 7/7 raw + median pad 4, ng > 4000 | 78.4 [75.3, 81.7] | 98.7 / 79.2 | 44.1 | 50.9 |
| 7/7 + Gaussian filter 1.5, no pad, ng > 2000 | 75.5 [72.1, 78.6] | 99.1 / 76.0 | 44.8 | 48.9 |
| GUI-only: box 5 + Gaussian 1.0, ng > 2000, no pad | 76.8 [73.3, 79.9] | 98.7 / 77.6 | 46.9 | 47.7 |
| 7/7 + Gaussian 1.5 + median pad 4, ng > 2000 | 81.6 [78.6, 84.5] | 98.8 / 82.4 | 51.4 | 45.4 |
| **best J: ident 7, fit 9, Gaussian 1.5, median pad 5, ng > 2000 (dev J 82.7)** | **82.5 [79.4, 85.6]** (block-10 [79.4, 85.6]) | 99.4 / 82.9 | 45.7 | 51.1 [45.5, 56.8] |
| **best E: ident 5, fit 7 `gaussmle`, Gaussian 1.0, pad 3, ng > 2500, 0.7 < s < 2.0 px, offset (11.1, −9.1) nm estimated on dev (dev E 64.7)** | 69.0 [65.4, 72.7] | 100 / 69.0 | 18.3 | **64.0 [60.6, 67.6]** (block-10 [58.9, 68.8]) |
| default detections + dev-estimated offset only | 71.7 | 100 / 71.7 | 29.0 | 59.4 [54.6, 63.4] |

Paired (iid, block-10 gives the same sign and similar width):
- Best J vs the reported default: dJ +10.8 [+7.6, +14.0]. Best E vs default: dE +8.8 [+5.6, +12.8]. Offset only: dE +4.3 [+3.5, +4.9].
- Challenge setting: best J gives J 95.3 and E 62.8 (RMSE 36.9) vs the reported 91.8 / 66.8. Best E gives J 85.0 and E 77.4 vs 66.8.
- Recall at precision ≥ 0.99 (threshold chosen on dev): 74.0 % → 82.9 %.
- Robustness on other windows of the stack, same dev-tuned recipe: dJ vs default is +11.3 (frames 5000–5999) and +12.7 (frames 12000–12999), both with block-10 CIs above 0. J of the best recipe is 82.9 and 79.4 there.

| HD (tuned on dev 0–99) | J % [95 % CI] | P / R % | RMSE nm | E |
|---|---|---|---|---|
| reported default | 47.8 [46.7, 48.8] | 99.8 / 47.8 | 49.6 | 28.0 |
| reported "tuned J" | 49.0 | 97.7 / 49.6 | 53.4 | 26.2 |
| reported "tuned E" | 46.9 | 99.9 / 46.9 | 48.5 | 28.1 |
| 7/7 + median pad 4, ng > 3000 | 54.4 [53.4, 55.6] | 97.2 / 55.3 | 57.4 | 26.7 |
| box 5/5 raw, no pad, ng > 3000 | 53.0 [51.8, 54.1] | 98.5 / 53.4 | 51.4 | 30.4 |
| **best J: box 5/5, median pad 3, ng > 3000 (dev J 57.6)** | **56.4 [55.2, 57.6]** | 98.2 / 57.0 | 53.1 | 31.3 |
| the LD recipe on HD (7/9, Gaussian 1.5) | 50.3 | 99.2 / 50.5 | 59.3 | 22.7 |
| **best E: box 5/5 raw, ng > 5000, s < 3.0 px, offset (13.0, −10.4)** | 48.1 | 100 / 48.1 | 38.5 | **35.3 [33.7, 36.9]** |
| default detections + offset only | 46.9 | 100 / 46.9 | 45.4 | 30.2 |

- Recall at precision ≥ 0.99 on HD: 48.6 % → 55.9 %.
- HD J ≥ 0.6 is not reachable. The ceiling, with every candidate kept (about 98 per frame), is a candidate recall of 64 % (dev). On HD, 34.1 % of ground-truth rows have another emitter of the frame within 250 nm and 51.0 % within 400 nm (LD: 4.4 % and 7.1 %) [raw]. A one-peak-per-blob detector cannot split these.
- Ablation: padding is worth +5.0 J (LD) and +5.4 (HD, box 7); the Gaussian filter +2.1 (LD); box 5 +4.0 (HD). The Gaussian filter hurts on HD. The LD best recipe has ONE tuned parameter (ng > 2000).

**Exact recipe for LD J ≥ 0.8**, tuned only on frames 0–999 [snippet]:

```python
P = 5; n = len(frames)                       # frames: (n,64,64) uint16 DN
pad = np.empty((n, 64+2*P, 64+2*P), np.uint16)
pad[:] = np.median(frames.reshape(n, -1), axis=1)[:, None, None]; pad[:, P:-P, P:-P] = frames
movie = picasso.localize._InMemoryMovie(pad, [])
ids, _ = picasso.localize.identify(movie, 2000, 7, threaded=False, gaussian_filter_sigma=1.5)
ids = ids[(ids.x >= 4) & (ids.x < 64+2*P-4) & (ids.y >= 4) & (ids.y < 64+2*P-4)]
locs, _ = picasso.localize.fit(movie=movie, camera_info=dict(Baseline=100, Sensitivity=45, Gain=300, Qe=0.9, Pixelsize=100),
                               identifications=ids, box=9, fitting_method='gausslq', multiprocess=False)
# x_nm = (locs.x - P + 0.5)*100, y_nm likewise; drop centres outside [0, 6400] nm
```

Picasso's blind border zone follows from `_local_maxima` (localize.py:197–232) [snippet]. This is why padding works (3 px top/left, 4 px bottom/right at box 7).

**Recall breakdown on LD test, best-J recipe vs default:**

| bin | default | best-J recipe |
|---|---|---|
| 0–300 nm from the border | 6.8 % (3/44) | 50.0 % |
| 300–450 nm | 31.2 % | 77.1 % |
| 450–700 nm | 58.8 % | 72.5 % |
| \|z\| beyond 450 nm | 27.3 % (9/33, Wilson [15, 44]) | 72.7 % |
| 500–1000 photons | 31.6 % | 57.9 % |

### Item 2: matching and metrics

- **Matcher.** I wrote an exhaustive-search matcher: maximum number of pairs within 250 nm inclusive, then minimum total distance. I ran 4000 random cases on a 50 nm grid (205 with several optimal pairings, exact 3-4-5 distances, empties). There were 0 mismatches against `match` in cardinality, total distance and membership of the optimum set.
- **Pooled scorer.** An independent pooled scorer with border / brightness exclusions applied after pairing agrees on 1200 comparisons (300 random multi-frame cases × 4 exclusion settings). tp, fp, fn, RMSE and E all match.
- **Real detections.** On real Picasso detections (LD ng 5000 and 4000; HD ng 5000 and 3000) the independent scorer equals the harness to the digit, in both settings (the all-emitter and the challenge setting).
- **Edge cases.**
  - No ground-truth emitter lies outside [0, 6400] nm [raw] (x 1.8–6398.8, y 0.8–6398.5).
  - A detection within tolerance of two emitters pairs with one, with no double count.
  - The pairing is made before the exclusions. A detection closer to an excluded dim emitter than to a kept one leaves the kept emitter a false negative although the detection is within 250 nm of it. This is a rule choice, [memory] consistent with the challenge, and rare (about 1 % of rows).
  - NaN / inf detections count as false positives when `border_nm` = 0 but silently vanish when `border_nm` > 0.
- **Definitions.** Pooled counts and the 2D Euclidean RMSE are as I remember the challenge [memory, moderate confidence].

### Item 3: conventions

All figures in this item are [own] from the raw data.

- **Frame offset.**
  - The stack has 19,996 frames, and the largest file frame number is 19,996, so the file is 1-based.
  - Births (first row of a consecutive-frame run of one fluorophore, n = 182): the signal in a 7×7 window is 9 photons at index = file frame − 2 and 4678 at file frame − 1, where 4678 ≈ 0.82 × the mean row photons of 5683.
  - Deaths (n = 178): 4622 at the last row and 25 one frame later.
  - Correlation of the 7×7 signal with the row's ground-truth photons: 0.949 at the claimed index, ≤ 0.27 at ±1, ≤ 0.08 at ±2.
  - Frames 5000–8000: 0.948 vs 0.23 / 0.30. HD: 0.944 vs 0.37 / 0.25.
- **Axes.** Correlation 0.923 for x → column, y → row. The other six dihedral mappings give −0.02 to 0.26. HD: 0.933 vs ≤ 0.31. (A 7×7 window cannot decide a half-pixel shift. Both ±0.5 px give 0.92–0.93.)
- **Pixel convention.** The whole-PSF centre of mass (radius 6 px) of the stack of 296 isolated bright in-focus emitters sits at (−0.2, −1.8) nm from the ground truth (frames 5000–7999: (−0.2, −0.9)). The other convention would give ≈ −50 nm.
- **The offset is estimator-dependent.**

| estimator | offset (nm) |
|---|---|
| centroid within radius 3 px | (+6.6, −5.3) |
| centroid within radius 5 px | (+1.1, −2.7) |
| Gaussian fit to the sub-pixel mean profile (r < 5 px) | (+9.6, −7.9), σ 1.355 |
| per-emitter centroid, median over 377 | (+5.1, −3.8) |
| own integrated fits (reproduces the md) | (10.9, −8.8) |
| Picasso LS 7×7, matched pairs | (+15.2, −11.4) |
| Picasso MLE | (+10.8, −9.6) |

  Per-emitter centroid offsets by z: (+0.6, −2.8) in focus and (+32, −13) for z in 300–800 nm. They show no dependence on x or y sub-pixel phase or on field position.
- **Replication.** The values were replicated on frames 5000–8000.

### Item 4: calibration

All figures are [own].

- **PSF width.**
  - My own LS integrated-Gaussian fits reproduce 664 emitters, median σ 1.369 px, in focus 1.277, |z| ≥ 350 nm 2.360.
  - By window, 7×7 / 11×11 / 19×19: median 1.338 / 1.369 / 1.408; in focus 1.282 / 1.277 / 1.308; |z| ≥ 350 nm 1.86 / 2.36 / 2.44 (peaking at 2.61 for 13×13).
  - The mean radial profile gives a Gaussian-equivalent σ ≈ 1.37 px (the optics value 0.93 is excluded).
- **Background and noise.**
  - Quiet frames (no ground-truth row in the frame or its neighbours): 2024 frames in the whole LD stack have mean 99.919 photons (DN 699.5). The radial bins are flat at 99.91–99.93 ± 0.01. The frame-mean sd is 0.230 against 0.233 for pure pixel noise. var/mean is 2.228 on the dev frames and 2.216 on another window.
  - Frames with emitters: the emitter-free pixel mean is 100.53 and the sd of the frame mean is 0.532. The frame mean regresses on the frame's total emitted photons with slope 7.4e−5 photon/px per photon, correlation 0.862. After removing it, the residual sd is 0.270 against 0.243 expected.
  - N2.LD quiet frames: mean 9.91 photons, var/mean 2.24. The two-point line through the N1 and N2 backgrounds is var = 2.227 x + 0.15 photon² (nominal 2.222 x + 0.078).
- **Halo.**
  - Around isolated single-emitter frames, the mean excess image per ground-truth photon has a cumulative sum of 0.837 within r < 3 px, 1.006 within r < 5 px, and 1.435 at 25 px, where it cuts off abruptly.
  - The same normalised profile appears in N2.LD (1.46), so it is a property of the simulated PSF. The cause is unknown.
  - Beyond r = 5 px: 11 photons/px at r = 5.5, 3 at r = 10, 0.9 at r = 17 for a 7400-photon emitter.
- **Density and photons.** Density and photon statistics reproduce: whole set 24938/19996 = 1.247 per frame; HD 12.645; photons median 4081 and 25 % quantile 1851 over the whole set.
- **Regression claim.** The claim that the pixel-variance-on-mean regression is uninformative holds. The slope of s² on x̄ for flat means is skewness × sd = 3.37, near the reported 3.36.

### Item 5: bootstrap

- **Implementation.** It is correct: frames resampled, pooled metrics recomputed, the paired version uses the same indices, and the circular blocks are right.
- **Toy coverage** (W/work/t10_bootstrap_coverage.py; nominal 95 %; 500 repetitions; 400 resamples):

| setting | iid frames | block 10 |
|---|---|---|
| independent frames, F = 500 (J / RMSE / P / R / E) | 0.952 / 0.930 / 0.942 / 0.958 / 0.948 | 0.956 / 0.914 / 0.944 / 0.964 / 0.952 |
| emitters persisting ≈ 3 frames with one detectability each, F = 500 | 0.636 / 0.930 / 0.882 / 0.610 / 0.690 | 0.906 / 0.934 / 0.930 / 0.904 / 0.902 |
| the same, F = 100 | J 0.632 | J 0.862 |
| paired difference, F = 500, independent | 0.948 / 0.963 / 0.932 (J / R / E) | 0.945 / 0.953 / 0.925 |
| paired difference, F = 500, persistent | 0.662 / 0.588 / 0.755 (J / R / E) | 0.887 / 0.892 / 0.915 |

- **Real data.**
  - Default Picasso J over 39 disjoint 500-frame windows of the LD stack has a window-to-window sd of 2.55, i.e. a half-width of 5.0.
  - The md's iid half-width is 3.4, and only 82 % of the windows fall inside it. The block-10 half-width is 4.7.
  - Paired intervals between thresholds of the same detector hardly change with block length.

### Item 6: reproducibility

- **Tests.** 27 passed.
- **calibration.md.** The rerun is identical except the runtime line ("10 s" → "11 s").
- **picasso_baseline.md.** The full script rerun (1000 resamples) is identical in every metric, CI and count. Only ms/frame, JIT seconds and the total runtime differ.
- **Checksums.** sha256 and size of all 16 files in /home/user/smlm_data match checksums.json.
- **CLI cross-check.** `picasso localize -a lq -b 7 -g 5000 -bl 100 -s 45 -ga 300 -qe 0.9 -px 100` on a TIFF of frames 1000–1199 gives 158 localizations, equal to the runner (max difference 0.0 in my pipeline). The CLI default fit is `mle` (`fit_methods = args.fit_method or ["mle"]`).
- **fetch.py against a local server (W/work/t15_fetch.py).**
  - A truncated file is removed and re-downloaded (OK).
  - One flipped byte at the same size with a record present gives exit "the server content changed. Re-run with --accept-new". The file is NOT repaired. `--accept-new` then overwrites the good digest with the corrupted file's digest and prints "downloaded".
  - A corrupted zip with no record ends in a raw `RuntimeError: CRC error` traceback and writes no record (safe).

### Item 7: overreach and omissions

- **Survivor bias of RMSE.**
  - LD: default RMSE 34.7 (428 TPs). The threshold-8000 setting has 25.3 nm (387 TPs). The default restricted to those 387 common TPs is 25.3 nm.
  - The 41 TPs lost at 8000 have RMSE 80.9 nm, median 2017 photons, median |z| 308 nm. The whole RMSE drop is selection, not a better fit.
  - HD: 44 lost TPs at 89.0 nm.
- **Argmax stability.** Bootstrapping the dev frames (2000 resamples), the E-argmax over the 20-value grid is 8000 in only 45 % (LD; 4500 in 23 %, 5000 in 16 %). On HD, 6000 in 24 % and 15000 in 25 %. E over 4500–10000 spans 0.8 points. The Jaccard argmax is stable (4000 / 3750; HD 3000 / 3250).
- **Photon bias.** Picasso photons (QE-corrected) / ground-truth photons for TPs: median 0.777 (IQR 0.67–0.87). By |z| 0.867, 0.756, 0.588, 0.494 for bins 0–100 / 100–200 / 200–350 / 350–800 nm. HD: 0.812. The Gaussian core captures about 0.8 of the emitted photons (SM1-3).
- **Detections per frame** are reported (0.856 vs 1.194 true per frame).
- **Offset in E.** Offset-corrected E is not reported in the md, only the centred RMSE.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| SM1-1 | major | **"Picasso baseline" = default + tuned min net gradient only (md l.5, S3).** The harness cannot express anything else: `PICASSO_DEFAULTS` (smlm_common.py:392–398), `_merge_params` (423–429) rejects other keys, and `tune_picasso` (510–548) is a ≤ 20-value grid of one parameter. Picasso's own "Gaussian Filter Sigma", frame padding (the 3–4 px blind border, localize.py:197–232), the box, and post-filters are unused. Dev-tuned alternatives give LD test J 82.5 [79.4, 85.6] vs 73.4, and HD 56.4 [55.2, 57.6] vs 49.0 (item 1 tables). S4 (recall 7 % at the border, 27 % beyond 450 nm) becomes 50 % and 73 %. In the challenge setting the gain is smaller (J +3.5 points on LD, +3.3 on HD). | **smlm_common.py:** add `gauss_sigma`, `pad` (median fill), `fit_box` and a post-filter dict (ng, photons, sx, sy) to `PICASSO_DEFAULTS`, `_picasso_worker` (432–458) and `tune_picasso`. **exp_smlm_picasso.py** (89–113): add a rung "Picasso tuned over its documented options on the dev frames" and the ablation ladder. **picasso_baseline.md:** rename "tuned"; change the S3 / S4 wording; state the tuning effort (≈ 150 configurations on dev) so the new method gets equal effort. |
| SM1-2 | major | **Baseline E understated.** (i) The offset is a calibratable constant: the default detections with a dev-estimated offset give LD E 59.4 (+4.3 [+3.5, +4.9]) and HD +2.2. (ii) A width filter (s < 2–3 px) plus offset gives LD 64.0 [60.6, 67.6] and HD 35.3 [33.7, 36.9]. (iii) "tuned (efficiency)" sits on a plateau (argmax 8000 in 45 % of bootstrap resamples, range 0.8 points) and its RMSE gain is survivor bias (41 lost TPs at 80.9 nm). (iv) Picasso's best-J and best-E operating points differ strongly (LD best J: E 51.1; best E: J 69.0), so single points are not a baseline. | **picasso_baseline.md / exp_smlm_picasso.py** (101–113, 125, 136–139): add E with the centred RMSE (already available as `rmse_centered_nm`, smlm_common.py:229–238) and the common-TP RMSE. Report a J–RMSE Pareto set of Picasso operating points rather than three points. |
| SM1-3 | major | **calibration.md (b), l.86 and l.88 (script 295–305): "background not spatially flat (illumination profile)" and "frame-level common-mode fluctuation, sd 0.48 photon" are wrong.** In quiet frames the background is flat (99.91–99.93 over all radii, 2024 frames) and the common mode is absent (0.230 vs 0.233). The effect is the PSF halo (SM1 Item 4): the emitter-free mask at 8 px (script l.31, 270–282) leaves halo (+0.4 photon; "100.3" should be 99.92). The halo itself is omitted: 44 % of the ground-truth photons arrive beyond r = 5 px (cumulative 1.435 at 25 px, abrupt cutoff), same in N2.LD. It biases any "Gaussian + constant" patch model and Picasso's photons (ratio 0.78, 0.49 at \|z\| > 350). The GT "intensity" corresponds to the core (cumulative ≈ 1.0 within 5 px). | **exp_smlm_calibration.py:** estimate the background on frames without ground truth within ±1 frame (99.92 ± 0.01); add a radial excess-profile table around single-emitter frames (cumulative to 25 px, cutoff); delete the illumination and common-mode sentences. **calibration.md:** state the halo as a PSF property. |
| SM1-4 | minor | **(+11, −9) nm is not "a constant offset" (calibration.md l.71).** It is the displacement of the PSF core relative to its centre of mass (≈ 0). It is estimator-dependent (table in item 3; Picasso 7×7 gives (15.2, −11.4)) and z-dependent ((+9, −7.5) in focus up to (+62, −32) beyond ±400 nm, md l.63–69; Picasso TPs at \|z\| 350–800 nm: (+70, −29), RMS 111 nm). | calibration.md l.71 and `score` docs: say "estimator- and z-dependent", name the estimator; keep reporting `bias_*`. |
| SM1-5 | minor | **Defocus widths are window-dependent** (md l.37–44: 2.358 px at 350–700 nm): 1.86 (7×7), 2.08, 2.36 (11×11), 2.61 (13×13), 2.44 (19×19). In focus it is stable (1.27–1.31). | Report the window with every width, or fit only the in-focus range. |
| SM1-6 | minor | **Intervals.** The headline intervals (exp_smlm_picasso.py:119–120) and the paired differences (l.125, `block=1`) assume independent frames. Real emitters persist ≈ 2.3 frames. The iid half-width 3.4 covers 82 % of 39 windows; the toy coverage is 61–69 % under persistence vs 90 % with 10-frame blocks. The block variant is given only for unpaired intervals (l.140–141). | Make `block=10` the default in `bootstrap_ci` / `bootstrap_diff_ci` (smlm_common.py:350, 371) or resample by fluorophore run id (consecutive frames with identical coordinates). Report block intervals for paired differences. State that dataset and simulation variability are not covered (docstring says it). |
| SM1-7 | minor | **Challenge rules are unreferenced.** Docstring l.22–25 (450 nm, 75 % brightest, after pairing) has no source in the repo; CLAIMS §6 says the Metric page was not read. I recall the same rules [memory, moderate confidence] but not the fine print: per row vs per fluorophore, whole-set quantile (1851 vs the dev quartile 2031), whether unpaired detections in the border zone are dropped (`score`, smlm_common.py:296–298), 2D Euclidean vs per-axis RMSE (a per-axis RMSE would add +7 to the LD default E: 62.5 vs 55.2). Edge cases: an excluded dim emitter steals a detection from a kept one; NaN detections vanish when `border_nm` > 0. | Cite the challenge page and section in the docstring; add tests for stealing and NaN; raise or count non-finite detections as FP (smlm_common.py:183–190, 296–298). |
| SM1-8 | minor | **fetch.py:187–204.** A same-size corrupted local file with a record is diagnosed as "the server content changed"; the file is never re-downloaded; `--accept-new` records the corrupted digest and prints "downloaded". A zip CRC error is a raw traceback. | On a recorded-digest mismatch of a local file, re-download to `.part` and compare; reserve `--accept-new` for a fresh download; catch the CRC error. Add a test for a same-size corrupted file. |
| SM1-9 | minor | **Omissions.** (a) Picasso photon bias (0.78 median, 0.49 beyond 350 nm). (b) HD structural ceiling: 34 % of ground-truth rows have a neighbour within 250 nm, 51 % within 400 nm, and the candidate-recall ceiling is 64 %; J ≥ 0.6 needs a multi-emitter method. (c) Dev and test share the same three microtubules, and `positions.csv` (fetched by fetch.py) holds the exact emitter sites, so a learned or prior-using method can exploit it. (d) Equal tuning effort with the new method is not stated. | Add to the md's text and the harness docstring. |
| SM1-10 | minor | **Noise model.** calibration.md l.94–110 fixes `a` at one operating point and calls `b` "not measurable". The N2.LD set already on disk gives a second point: slope 2.227, intercept 0.15 photon². | Add the N2.LD background as a second point. |
| SM1-11 | cosmetic | Dev-frame photon quartile 2031 (calibration.md l.128) vs the whole-set 1851 used by the harness; timings inside an otherwise deterministic md; the runtime sentence prints the JIT time twice (exp_smlm_picasso.py:134–135); "tuned (efficiency)" label. | Tidy. |

## Sources

- **Repository (commit 40a7252):**
  - experiments/smlm_common.py: 22–25, 183–190, 229–238, 241–309, 350–386, 392–458, 510–548.
  - experiments/exp_smlm_picasso.py: 89–141, 180–189.
  - experiments/exp_smlm_calibration.py: 31, 270–305.
  - tests/test_smlm_common.py; data/smlm/fetch.py: 186–204; data/smlm/checksums.json.
  - results/smlm/calibration.md: l.31–44, 54–71, 86, 88, 94–110; results/smlm/picasso_baseline.md: l.5, 38–61, 81–117.
- **Raw data** (/home/user/smlm_data): activations.csv, the stacks, the stats files for MT0.N1.LD / HD and MT0.N2.LD / HD.
- **Picasso 0.11.3 source** (venv_domain/lib/python3.13/site-packages/picasso/): localize.py 197–232 (`_local_maxima`), 1501–1566 (`identify_in_frame`, float32 cast), 1890–2046 (`identify`, Gaussian filter and temporal median), 2661–2680 (`_to_photons`), 3363–3400 (`parse_gauss_code`), 7866+ (`localize_frames`); fitting/gaussfit.py l.15 (pixel-centre models); gaussmle.py 1–35, 1055–1100 (deprecated integrated MLE).
- **The referee's scripts** (refSM1/work): pic.py (pipeline), t4/t18/t21 (sweeps), t14_final.py (test scoring), t5/t5b (independent scorer), t6/t6b (conventions), t7/t7b/t8/t9 (calibration and halo), t10 (bootstrap toy), t13/t22/t25 (breakdowns), t15 (fetch), t19 (argmax), t27/t28 (windows). Logs in refSM1/logs.
- **[memory]:**
  - Sage et al., Nat. Methods 16:387 (2019), Methods: matching, J, RMSE, E, 450 nm border, brightest 75 %.
  - Sage et al., Nat. Methods 12:717 (2015).
  - Smith et al., Nat. Methods 7:373 (2010): MLE with the integrated Gaussian.
  - Zhang et al., Appl. Opt. 46:1819 (2007): sigma = 0.21 lambda/NA.

## Lesson

A baseline is only as fair as the options the harness lets it use. Here the API admitted one parameter, so "tuned Picasso" stayed 9 J points and 4–9 E points below what Picasso's own documented options plus a trivial pre-processing reach on the same dev frames. Build the baseline as a ladder (default, net-gradient-tuned, option-tuned) with an ablation, and give the new method the same tuning budget.

A calibration statement must be tested on the subset where the nuisance is absent. Frames without any emitter would have shown that the "illumination" and "common mode" were the emitters' own halo, and that the background is exactly 100 photons.

Check interval claims against the data's own window-to-window spread, not only a toy: 39 disjoint windows showed the iid interval was too narrow.
