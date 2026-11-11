# Referee report IS1b-S: single-molecule localization experiments, baselines and cost claims

Repository `/home/user/Lyapunov-Convergence-Kernels-Showcase`, branch `claude/lyapunov-context-setup-yg51ew`. The reviewed content is commit de1ec46. Every reviewed file was byte-identical to `git show HEAD:<file>` at the start and again at the end (sha256 in the private copy). HEAD moved to eee967e, an astronomy commit that does not touch these files. Uncommitted edits by other jobs (`HANDOFF.md`, `tests/test_imaging.py` +146 lines, `tools/mutants/wp2_imaging_own.py`, astro files) were not in my copy and were not reviewed.

Rule compliance:
- I did not modify the repository. `git status` shows only the other jobs' changes.
- Work was done in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_IS1b_S_work/`, with the subfolders `repo/` (private copy), `scripts/`, `logs/` and `out/`.
- I used `OMP_NUM_THREADS=1 nice -n 10` and at most 2 busy processes. The machine load average was 4 to 7 on 4 cores.
- I killed nothing, used no network, pushed nothing, and posted nothing.
- Interpreters: repo venv (Python 3.13.18, numpy 2.5.4, scipy 1.18.1), Picasso 0.11.3 in `venv_domain`, torch in `venv_learn`. All Python ran with `-I`.

## Summary

There is no blocker. No final number is test-informed, and every number I re-derived reproduces. The problems are in how the headline is worded (findings 1 and 2) and in disclosures and robustness (findings 3 to 16).

Claim-by-claim verdict:

| claim | verdict |
|---|---|
| 1, tuned part (E 3 wins and 1 tie, J 2 wins and 2 ties, RMSE equal or lower) | supported; re-derived with my own code |
| 1, untuned part ("beats Picasso's defaults on three of four") | literally true, but it measures Picasso's unit-dependent default threshold (finding 1) |
| 2 (2013 sets, 8 of 8) | supported; the tests I added for the model-match advantage do not change it (finding 8) |
| 3 (label-free width costs about nothing) | supported (≤ 0.4 E). The BIC stop is a separate cost (finding 2) |
| 4 (real FRC) | supported against R0. The "equal count" part is understated (finding 9) |
| 5 (CNN, mismatch) | supported as worded. Caveats in finding 11 |
| 6 (cost) | numbers right, range selective (finding 10) |
| 7 (1-D) | supported; "on a par" means indistinguishable with 36 targets (finding 14) |

## Findings

### 1. [major] "The untuned ISS beats Picasso's defaults on three of four" measures a unit-dependent default threshold

- **Code and report at fault:** `experiments/exp_smlm_picasso.py:55` (`DEFAULT_NG = 5000.0`, "for every set"). `results/smlm/iss.md` rows "a priori default ISS vs Picasso defaults (R0)". `results/smlm/picasso_baseline.md` summary tables.
- **Mechanism:** The N2 sets have about 10 times fewer photons (median 441 against 4081). The default net gradient 5000 DN is then far above the signal. R0 reaches E 18.7 (N2.LD) and 6.2 (N2.HD), J 22.1 and 18.1. With only the threshold tuned (R1E), E is 42.4 and 18.5.
- **Evidence:** I computed paired intervals (own scorer, 4000 block resamples, script `scripts/iss0_vs_r1.py`). Entries are the E difference, a priori ISS minus the Picasso rung, in the order N1.LD, N1.HD, N2.LD, N2.HD.
  - Against R0: +1.9 [-4.0, 7.6] tie; +10.3 [8.0, 12.7] win; +26.7 [21.7, 31.2] win; +14.4 [12.3, 16.4] win.
  - Against R1E+c (threshold tuned, offset removed): -2.5 [-8.0, 2.8]; +8.2 [5.9, 10.6]; +1.4 [-2.6, 5.1]; +0.9 [-0.8, 2.6]. That is three ties and one win.
  - Against R2E+c: -6.2 [-11.0, -1.6] loss; +1.6 [-0.7, 3.9]; -1.6 [-4.9, 1.4]; -4.0 [-5.7, -2.3] loss.
  - `iss.md:650` itself says the mean E rise from tuning is +15.9 for Picasso and +7.7 for the ISS.
- **Fix:** Report the untuned ISS against R1E+c as well. State that the R0 gap on the N2 sets is a scale effect of Picasso's default threshold, or give R0 a camera-scaled threshold. Word the claim as "an untuned ISS matches a Picasso whose threshold alone is tuned and is below the fully tuned Picasso on 2 of 4 sets".

### 2. [major] "No regularization parameter, no training data" is not what the evidence supports

- **Code at fault:** `src/core/imaging.py:464-471` (`stop_bic`). `imaging.py:312-313` and `:385-387` (the weights follow the model; `chi2` is the model-weighted mean square). `experiments/exp_smlm_iss.py:60` (`SNR_MIN`, `SNR_EXTRA`).
- **(a) The headline win is by the tuned ISS.**
  - It uses κ, configuration, width bound and offset, all chosen with labelled development frames.
  - The untuned ISS is below R2E+c on 2 of 4 sets (finding 1).
- **(b) The parameter-free variant pays.**
  - Label-free BIC default against the κ=5 default, E difference per set: -4.2 [-7.4, -1.5], -1.2 [-3.0, 0.8], -10.2 [-14.8, -6.2], +2.8 [1.2, 4.4].
  - Against the E-tuned ISS: -15.0, -9.5, -12.5, -6.1.
  - Subtracting point estimates, the BIC default is about 10 E points below R2E+c on N1.LD (53.3 against 63.3) and on N2.LD (35.1 against 46.9). It is about level on the two HD sets.
- **(c) The BIC stop admits satellites.**
  - Script `scripts/bic_analysis.py`, development frames 0-399 of N1.LD, a priori configuration.
  - Repo BIC: 13 false positives, precision 97.1 %, E 56.4.
  - BIC with the log-variance term added: 9 FP, E 58.2.
  - BIC with weights frozen at the end of the path: 6 FP, E 58.6.
  - κ=5: 4 FP, E 60.7.
  - κ=5.77 (the nominal single-unit BIC threshold √(4 ln 4096)): 1 FP, E 59.9.
  - The extra false positives of the repo BIC sit a median 216 nm from a real emitter (92 % within 400 nm) and carry 0.25 of its flux. Their fitted width is a median 2.2 px.
  - The criterion `N·χ²` is not a likelihood. A unit that raises the model also lowers the weights around it, and the `Σ log var` term is missing.
  - The search is truncated by `snr_min=3` and `snr_extra=3`.
- **(d) The inputs are not free.** The method needs the PSF width and family, a halo part, a high-pass, the noise law and about 1000 unlabelled frames for calibration. `iss.md` line 7 lists "sigma0 = 1.40 px" under "Calibration (… no labels)". The docstring (`exp_smlm_iss.py:15-16`) says that value comes from labelled emitters. `b0 = 0.078` is the dataset page's nominal read-out term, not a measurement.
- **Fix:** Reword to "no labelled training data in the label-free variant, which pays X E points". Add the log-variance term or freeze the weights, and re-run parts B and L. Print the BIC deficit next to claim 3.

### 3. [minor] The order of events is not recorded, and "declared before its test run" is weaker than it reads

- **Where:** `exp_smlm_iss.py:11-20` ("nothing is tuned on test frames") and `:93-94` (`DEFAULT`).
- **Evidence (repo history plus the build scripts and logs in the coordinator's scratchpad):**
  - `git show feea3df:experiments/exp_smlm_iss.py` has an `op_table` over (config, κ) only. There is no width bound, no offset and no round 2.
  - `run_iss_test.log` (17:09) shows the first test run with that rule. The dev choices (hp5h κ8, hp3 κ5, hp3 κ6, hp3h κ6) are exactly what I recompute from `iss_raw.npz`. Test E was 63.3, 41.4, 45.5, 24.1 and J was 77.7, 65.0, 57.0, 44.8. These numbers appear nowhere in the repository.
  - Afterwards the width bound, the dev offset and the round-2 variants (D2) were added (`5527d7b`). The Picasso ladder with +c and a σ window followed (`s4_brief.md`, 17:15). The test run was repeated (`run_iss_final.log`, `run_iss_final2.log`).
  - `exp_smlm_iss.py.bak` (19:07) has no `DEFAULT`. T0 ran 19:24-19:35, after the ISS-E and ISS-J test results and the whole ladder's test results existed.
- **Impact, measured (`scripts/round1_check.py`):**
  - D2 does not matter. The base-only E point (hp5h κ6, width ≤ 2.5, offset) re-run on the N1.LD test frames gives E 68.2, J 74.9, RMSE 19.4, identical to the final point. The dev E gain from D2 is +0.42, +0.33, 0.00, +0.26.
  - The rise from 63.3 to 68.2 and so on comes from the offset and the width bound, which mirror Picasso's +c and σ window.
  - `DEFAULT` (hp5h, κ5) is not the dev optimum on 3 of 4 sets.
- **Fix:** Add an "order of events" paragraph like the one in `exp_smlm_real.py:36-50`, and print the first-round numbers.

### 4. [minor] "Tuned like Picasso … with the same objective" is not literally true; the effect is nil

- **Where:** `exp_smlm_iss.py:421-444` (the ISS E point maximises dev E after offset removal). `smlm_common.py:918-927` and `exp_smlm_picasso.py:164-165` (Picasso R2E maximises uncorrected dev E, and +c is applied afterwards).
- **Test:** `scripts/retune_picasso.py` and `scripts/exhaustive_picasso.py`.
  - Full 108-structure search with the corrected objective on N1.LD: same winner and the same filters, dev E(c) 65.93.
  - Exhaustive search of all post-filter combinations for the stored structure: 8505 combinations on N1.LD and 8820 on N2.LD. The optimum equals the ladder's greedy ascent: dev E(c) 65.93 and 49.64, test E 63.3 and 46.9.
  - Re-tune of the stored structure on N1.HD and N2.HD: test E 36.9 (against 36.7) and unchanged 24.6.
- **Fix:** State that the objectives differ and that it does not matter.

### 5. [minor] Calibration information is not equal; I tested the part that could matter

- **Asymmetry:** The ISS gets a labelled σ0, the nominal b0, the noise slope and a background surface. Picasso gets only its camera parameters.
- **Tests:**
  - Giving Picasso the ISS's static background calibration (Tubulins_I, same structure, 2013 grid, corrected objective): test E 63.4 against 63.5 (`scripts/bgsub_picasso_2013.py`).
  - The ISS is insensitive to σ0. A 0.5× or 2× width gives identical output on a 3-emitter frame, and the learned sweep gives dE -0.5 and -1.3.
  - The ISS is sensitive to the noise scale (finding 11).
- **Tuning effort, in counts:**
  - ISS: 168 operating points for 2016 (84 for 2013). Picasso: 12,474 to 14,618 (2016) and 12,762 to 19,063 (2013).
  - Dev-to-test optimism: ISS-E 70.4→68.2, 48.9→46.8, 49.6→47.6, 30.9→29.5. Picasso R2E(c) 65.9→63.3, 39.7→36.7, 49.6→46.9, 27.3→24.6.
  - Neither effort nor information is "equal", but no asymmetry I could test changes a win or a tie.
- **Fix:** Say so, and correct `iss.md` line 7 (finding 2d).

### 6. [minor] The label-free width estimator is biased to in-focus emitters and is not stable across noise levels

- **Where:** `exp_smlm_iss.py:149-187`. Script `scripts/sigma0_analysis.py` (dev frames 0-299 of N1.LD).
- **Evidence:**
  - 178 accepted spots. In the brightest quarter 49 % have |z| < 100 nm and 44 % have 100-200 nm. Only 6.6 % have |z| ≥ 200 nm, against 27.4 % of isolated bright ground-truth emitters.
  - Estimate 1.305 px. In-focus median 1.272. All-z labelled median 1.369 (I reproduced the 664 emitters, medians 1.277, 1.393, 1.657, 2.360).
  - N1 gives 1.305 and N2 gives 1.383 for the same PSF (178 against 86 spots).
- **Impact:** Harmless (finding 5). Label-free against BIC with the labelled width changes E by at most 0.4: -0.4 [-1.4, 0.4], -0.2 [-0.9, 0.5], -0.0, -0.0.
- **Test gap:** `tests/test_exp_smlm_iss.py` plants one width only.
- **Fix:** State the bias, and test with a defocus mixture.

### 7. [minor] The four 2016 sets are two scenes

- **Evidence:**
  - N1.LD and N2.LD have identical ground truth: all 24,938 rows agree in frame, x, y and z. Photons are scaled by about 0.108.
  - N1.HD and N2.HD likewise (31,612 rows).
  - 98.7 % of the LD fluorophore sites also occur in HD.
  - The max photon count is 10,577.
- **Consequence:** "3 of 4 sets" is not three replicates. The 2013 sets are paired the same way, as `smlm2013.md` says.
- **Fix:** Say "2 scenes × 2 noise levels" next to the counts.

### 8. [minor] A statement about the PSF family in `smlm2013.md` is inaccurate; the model-match advantage is not what drives 8 of 8

- **Where:** `results/smlm/smlm2013.md:138`. It says no kernel of either method contains the 2016 halo. All four E-tuned ISS configurations on 2016 are 'h' configurations with a core + halo mixture (25 % at 3× width) and a high-pass. Dev gain of the halo part is about +1 to 2 E (hp5 against hp5h).
- **Test of the 2013 Gaussian advantage (`scripts/refit_picasso_2013.py`):** I refit Picasso's own R2E and R2J detections with a pixel-integrated Gaussian (7×7, dev offset, free σ or σ fixed at σ0) and rescored. No algorithmic step changes. In the list below, "original" is the stored R2E+c or R2J+c.

  | set | rung | refit, σ fixed | refit, σ free | original |
  |---|---|---|---|---|
  | Tubulins_I | R2E | E 64.1 | 62.8 | 63.5 |
  | Tubulins_I | R2J | E 62.6 | 59.2 | 58.8 |
  | Bundled LS | R2E | E 80.7 | 77.8 | 82.1 |
  | Bundled LS | R2J | E 72.9 | 67.0 | 76.2 |

  Picasso's own point-sampled fit is as good as or better than my integrated refit, except for J-tuned Tubulins_I, where it gains at most 3.8 E (ISS still +4.4 ahead). The 8 of 8 pattern is unchanged.
- **Fix:** Correct the sentence and cite this test.

### 9. [minor] Real-data FRC wording understates the equal-count result

- **Evidence (my own FRC, `scripts/frc_indep.py`, reproduces the repo):**
  - Alternating-block resolution: R0 55.1 nm, ISS 40.2, Picasso at net gradient 1000 47.1.
  - Paired ratio against R0: 0.82 [0.71, 0.94]. Against net gradient 1000: 0.96 [0.84, 1.05], a tie.
  - At equal count (42,097): ISS thinned 53.4 [51.8, 59.3] in the repo (55.5 [48.9, 61.0] in my run), Picasso net gradient 1000 thinned 49.5 [49.2, 51.5].
  - `scripts/real_picasso_scan.py` scan of Picasso's threshold (net gradient 5000 → 150): 55.1, 49.1, 47.1, 47.1, 47.1, 47.2, 49.1, 49.1 nm. Picasso never gets below 47.1. It caps near 105 to 114 localizations per frame against the ISS's 160.
- **Reading:** The gain over R0 disappears at equal count and reverses against a lower-threshold Picasso. "Mostly disappears" is too weak.
- **Disclosed in `exp_smlm_real.py`:** the noise-scale revision (after 2 frames, never the FRC) and the headline split changed after the first run.
- **Fix:** "At equal count the ISS is no better than Picasso."

### 10. [minor] Cost reporting is selective and one baseline timing is an outlier

- **Where:** `results/smlm/speed.md`, `iss.md:456`, `real_frc.md`, `smlm2013.md`.
- **Evidence:**
  - `speed.md` uses 12 frames per set and the default configuration. The tuned HD configuration costs 5.3 to 6.5 wall seconds per frame over 400 frames (loaded machine).
  - 2013: 5.3 to 5.7 CPU s per frame on Tubulins and 15.9 on the dense set. Real data: 7.76 CPU s per frame.
  - So "0.36 to 1.8 s" is the best case. The "2 to 4 orders" ratio is conservative.
  - Picasso on N2.LD is 0.0297 s per frame against 0.0005 on the other three sets. That is unexplained (JIT or load) and would give a ratio of 37.
  - My re-measure, process CPU time (`scripts/speed_check.py`):

    | set | ISS full | ISS fast | Picasso R0 | speed-up | ratio to Picasso |
    |---|---|---|---|---|---|
    | N1.LD | 0.354 | 0.110 | 0.00064 | 3.21 | about 555 |
    | N1.HD | 0.892 | 0.444 | 0.00039 | 2.01 | about 2302 |

  - N1.HD 0.89 s on 4 frames against 1.78 s on 12 frames. The frame mix matters.
- **Fix:** Report medians over all test frames for the tuned configuration, the range across all experiments, and repeat the Picasso timing.

### 11. [minor] Mismatch sweep and the CNN stand-in

- **Mismatch injection (`exp_smlm_learned.py:535-547`):**
  - It is not an artefact. The weights use variance `a·(0.5·bg)+B0 ≥ 112.6` photon² (true ≈ 230), with no negative or tiny variance. This is a genuine 2× variance underestimate (κ_eff ≈ 3.5).
  - The reported drop is -18.6 [-27.3, -11.6]. I reproduced it on 80 test frames (`scripts/mismatch_check.py`):

    | condition | E | note |
    |---|---|---|
    | nominal | 61.5 | |
    | bg ×0.5, a unchanged | 41.2 | precision 77.2 %, 1.54 against 1.12 detections per frame, 183 s against 29 s (6× slower) |
    | bg ×0.5, slope recomputed (a/0.5) | 62.3 | |
    | bg ×2 | 60.8 | |
    | bg ×2, slope recomputed (a/2) | 60.7 | |

  - The ISS is insensitive to the background level. It is fragile to the noise scale. The "background" label is misleading, and the calibration (a = median((var−B0)/bg)) would absorb a wrong bg.
  - On the paired CI the E difference between CNN and ISS under this mismatch is a tie: +8.8 [-0.8, 18.6]. The J difference is +13.5 [8.8, 18.3].
  - The CNN is a single seed. The seed-to-seed |ΔE| is up to 4.4. The sweep uses only 200 LD and 100 HD frames. The ISS has no photon-count counterpart.
- **CNN as stand-in (code read: `train_model`, `validation_data`, `infer_model`):**
  - All choices are on held-out development or simulated frames, and test inference happens after `selection.json`. The pilots used dev frames, including the held-out ones.
  - The CNN loses to the E-tuned ISS on all four sets: -7.1 [-12.8, -2.1], -8.5 [-10.2, -6.8], -7.7 [-11.4, -4.4], -3.0 [-4.3, -1.8].
  - The hindsight threshold gives 62.5 against 61.2, so the gap is not a threshold effect.
  - The decoder has one knob. The ISS and Picasso have width and photon filters.
  - Training budget is 2000 × 16 on 800 frames and 289k parameters. The labelled-data table shows E still rising by about +2 per doubling.
- **Fix:** Keep "CPU-scale" in the claim, and call the mismatch "noise scale" for the ISS.

### 12. [minor] The packaged `fit_frame` fails silently on some inputs

- **Where:** `imaging.py:571-578` (`sw_fn`, `y`). Script `scripts/adversarial.py`.
- **Results (12 categories):**
  1. Empty frames: 0 detections in 100 frames at κ=5 and under BIC, 7 at κ=4.
  2. Bright emitters (5e3, 1e5, 1e6 photons): κ=5 gives 1, 1, 4 units; BIC gives 1, 7, 8. The extras are ghosts at the 3.2 px width bound (about 2.7e3 photons each at 1e5). The data max is 10,577 photons, so this does not occur in the benchmark.
  3. Two equal emitters at 1.0, 1.5 and 2.0 px are returned as one unit, also with `r_ex=0`. At 2.5 px and above they are resolved.
  4. Border emitters from -1.5 to 64.6 px are found, centres pulled in by up to 0.67 px.
  5. An exactly symmetric pair is found.
  6. Noise variance assumed ×0.25: 10 units for a 3-emitter frame and 6.7 false units per empty frame (400 in 60 frames). ×4 is unchanged.
  7. A NaN or inf pixel makes the whole frame return 0 detections. Constant and zero frames give 0. An 8×8 frame works.
  8. A bg argument ≤ 0 or a surface with negative pixels gives a negative variance and NaN weights: 0 detections, warning only (`invalid value encountered in sqrt`). The 2013 pipeline avoids this with its own `np.maximum(bg, 0)` (`exp_smlm2013.py:322`). Its background surface has pixels down to -2.2 photons. Assuming bg=0 gives 23 units for one emitter.
  9. Clipping at 600 photons gives 1 unit, with amplitude half and σ 2.12.
  10. 60 emitters (cap 40) give 37 units in 38 s.
  11. A dark spot makes 3 ghost positive units at σ 2.7 px.
  12. A non-square 64×40 frame works.
- **Fix:** Validate finiteness, clamp the variance, and warn.

### 13. [minor] Baseline scope and the "every option" wording

- **Not in the grid:** spherical and rotated fits, `lpx`/`lpy`/`bg` filters, GPU fits, finer σ grids.
- **Test:** Spherical fits on N1.LD (`scripts/retune_fit.py`): dev E(c) 64.7 and 63.7 against 65.9, test E 64.2 and 63.6 against 63.3. No gain. `lpx`/`lpy` were not tried.
- **Scope:**
  - Picasso is the only baseline, and it is a single-emitter localizer. On the dense 2013 bundle its J is 19 to 20 % against the ISS's 68 to 77 %, which is outside its regime. No multi-emitter baseline was run.
  - Picasso 0.11.3 has only point-sampled Gaussian fits (the integrated version is commented out in `gausslq.py:84-88`), so a pixel-integrated Picasso option does not exist.
  - The SFW baseline is not always converged. `scripts/kkt_sfw.py` over 144 runs: 132 flagged converged, 124 have a fine-grid certificate ≤ 1.005, worst 8.8.
- **Fix:** Say "the options of this grid", and name the missing ones.

### 14. [minor] Statistics

- **Few blocks in 2013:** Tubulins test = 10 blocks, dense bundle = 6 blocks. Percentile bootstrap intervals are unreliable at that size.
- **No multiplicity correction** over about 30 comparisons.
- **The 1-D study treats 432 trials as independent** (Wilson intervals). There are 36 targets × 12 trials.
- **My target-clustered paired bootstrap** (exact-order rate, ISS-BIC minus other):

  | other method | difference [95 %] |
  |---|---|
  | OMP-BIC | -0.021 [-0.049, +0.002] |
  | SFW-BIC | +0.044 [-0.025, +0.111] |
  | SFW-discrepancy | -0.009 [-0.069, +0.049] |
  | CPGD-BIC | -0.025 [-0.076, +0.025] |
  | CPGD-discrepancy | +0.025 [-0.044, +0.093] |

  At s=0.3, ISS-BIC is 0.52 and CPGD-BIC is 0.67. "On a par" means indistinguishable, with a power of about ±5 to 7 points.
- **Order-of-entry:** −1.33 [−2.07, −0.67] over all cells. At separation 1.6 widths, flux ratio 1, 2000 photons: +21. At separation 3.0 widths, flux ratio 1: −47 points.

### 15. [minor] Reproducibility

- **Regeneration:** all md files regenerate byte-identically from raw (see Verified), including `smlm2013` with `SMLM_DATA_2013=/nonexistent`.
- **Gaps:**
  - `iss.md --from-raw` is not "from raw alone". It also needs `/home/user/smlm_data` ground truth and the ladder file.
  - `requirements.txt` is unpinned (`>=`) and omits Pillow, Picasso and torch. The environments live outside the repo.
  - Detections depend on the scipy `nnls`/`least_squares` implementation.
  - `calibration.md` and `picasso_baseline.md` have no `--from-raw` and need hours.

### 16. [minor, cosmetic] Arithmetic in the per-set header of `iss.md`

"168 operating points (4 model configurations x 7 x 3)" gives 84, not 168. `exp_smlm_iss.py:640-642` prints `len(CONFIGS)`, but the table contains 8 configurations (the base four plus four round-2 variants).

## Verified (re-derived with my own code unless stated)

- **2016 E, J, RMSE.** Own Hungarian matching at 250 nm per frame, J, RMSE, E, block bootstrap. All four sets, ISS-E, ISS-J, R0, R2E+c, R2J+c (scripts `verify_2016.py`, `indep.py`). The development operating points are recomputed from the `D|…` arrays, and the same points are selected.

  | set | ISS-E J / RMSE / E | R2E+c J / RMSE / E | paired dE | paired dJ |
  |---|---|---|---|---|
  | N1.LD | 74.9 / 19.4 / 68.2 | 68.7 / 19.1 / 63.3 | 4.9 [2.6, 7.4] | 6.2 [3.7, 8.9] |
  | N1.HD | 59.4 / 34.4 / 46.8 | 51.5 / 40.6 / 36.7 | 10.1 [8.4, 11.7] | 8.0 [6.5, 9.4] |
  | N2.LD | 58.5 / 32.0 / 47.6 | 57.3 / 31.5 / 46.9 | 0.7 [-0.9, 2.4] | 1.2 [-0.6, 2.8] |
  | N2.HD | 45.4 / 44.5 / 29.5 | 38.7 / 43.9 / 24.6 | 4.9 [3.8, 6.0] | 6.6 [5.7, 7.6] |

  The J-tuned rows match the md as well: J 85.1 / 68.8 / 65.8 / 54.0 for ISS-J, and paired dJ 1.9 [-0.2, 4.1], 12.3, 0.8 [-1.8, 3.4], 5.5. Intervals from a different RNG (B=4000) differ from the md's by at most 0.3. The win, tie and loss readings are correct: E 3/1/0, J 2/2/0.
- **2013.** All four sets. Tubulins_I ISS-E E 67.8 against R2E+c 63.5 (dE 4.3 [3.5, 5.0]). Tubulins_II 57.1 against 52.9 (4.2 [3.1, 5.4]). Bundled LS 84.0 against 82.1 (1.9 [0.8, 3.0]). Dense 21.2 against 7.0 (14.2 [12.9, 15.5]). J-tuned dJ 4.3, 5.0, 2.0, 57.2. RMSE worse on the dense bundle (+26.0) and Bundled LS E-tuned (+3.0). `scripts/verify_2013.py`.
- **Raw detections come from the code.** A fresh run of `exp_smlm_iss.run_frames` on 12 test frames of N1.LD (hp5h_r25 κ6) and 10 of N2.LD (hp3h κ5) is bitwise equal to the stored arrays. `scripts/rerun_iss.py`.
- **Dev and test ranges are disjoint everywhere**, in `iss_raw`, the ladder raw, the 2013 raw and the learned protocol. 2016 calibrations use development frames only (`exp_smlm_calibration.CAL_EXCLUDE`, my 664-emitter refit). 2013 pre-test steps use dev ranges only (code read). The sha256 of 32 data files and the digests quoted in the md files match.
- **Real data.** Counts 42,099 and 80,047 (`real_raw`). Picasso R0 rerun: 42,097 localizations in the field. My FRC reproduces 55.1, 40.2, 47.1 nm and the thinned values (finding 9).
- **Byte-identical regeneration** of `iss.md`, `smlm2013.md`, `learned.md`, `speed.md`, `noise_order.md`, `order_pairs.md`, `order_pairs_bicc.md`, `order_pairs_snr5.md`, `real_frc.md` (`scripts/run_fromraw.sh`).
- **Tests.** 192 passed (imaging, baselines, exp_smlm_iss, smlm_common, smlm2013, smlm_speed, noise_order, order_pairs), 19 passed (learned localizer, venv_learn), 36 passed (frc, smlm_real).
- **1-D trials.** Six trials re-run: all 16 K-hat columns identical (`scripts/rerun_noise_order.py`).
- **Learned pipeline.** Code read: no test frame enters training, checkpoint, threshold, offset or pilot selection.
- **Asymmetry tests that came out immaterial:** corrected objective (finding 4), exhaustive filter search (finding 4), background calibration (finding 5), integrated-Gaussian refit (finding 8), spherical fits (finding 13), D2 round (finding 3).

## Could not check

- The 2013 challenge rules, and the primary source for the 2016 metric rules (no network); E, the 450 nm border and the 75 % rule remain "from memory".
- That no 2013 test frame was looked at before the protocol was declared. Part-file timestamps and code paths are consistent with it. "100 frames of each set for pixel statistics" is not specified.
- Anything the building agents did without logging. The timeline in finding 3 is reconstructed from scratchpad logs and timestamps.
- Regeneration of `calibration.md` and `picasso_baseline.md`, and a repeat of the full Picasso search on all four sets (several hours each).
- A re-localization of the real sequence by the ISS (about 1 h) and any retraining or re-inference of the CNNs. The CNN and real-data rows rest on the stored raw, whose md I regenerated byte-identically.
- Picasso `lpx`/`lpy`/`bg` filters and other versions.
- Sensitivity to numpy and scipy versions (one environment only).
- Mutation testing of `imaging.py`.
- The `S` and `V` parts of `iss.md` (stop rules, design variants) were not re-derived individually.
- The uncommitted `tests/test_imaging.py` changes of another job.
