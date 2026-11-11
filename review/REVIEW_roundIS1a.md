<!-- Archived verbatim from the referee's final message (round IS1a, pre-test review of the astronomy experiment, 2026-11-07). The coordinator added only this comment. -->

# Referee round IS1a: pre-test review of the astronomy experiment (M2 patch of SDSS against HUGS)

Repository `/home/user/Lyapunov-Convergence-Kernels-Showcase`, branch `claude/lyapunov-context-setup-yg51ew`, commit 1748a01 (the astronomy files equal the commit; `results/astro/baseline_catalogs.npz` is untracked). Review window 20:36 to 21:55 on 2026-11-07. Everything below was done in a private copy (`/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_IS1a/`, scripts in `work/`); the repository was not modified.

**Rule compliance.** I did not run the path (`imaging.py`) and did not score any baseline on a test patch. On the test patches I only read truth positions, pixel statistics and the baseline catalogue arrays (geometry, counts, row counts, array equality), none of it a score. The one permitted exception was used: the shipped DAOPHOT catalogue was scored on the benchmark patch (310, 630), for the shift and for its magnitude calibration only. Every number that involves detections and truth comes from the development patches (110, 630) and (210, 630), synthetic inputs, or the stored dev-patch baseline catalogues.

**Recommendation: do not release the test run as it stands.** One blocker (two of the four test patches are mostly outside the HST footprint, so their purity and F measure footprint geometry), and several majors that are cheap to fix before the run: the protocol-critical selection code has no numeric test; the photutils baseline is under-tuned and half of its tuning space is a no-op; the DAOPHOT row is aligned with test-patch truth; the tuning objective cannot see false positives fainter than 22 mag, so kappa sits at and beyond the lower grid edge; the freeze is not mechanical. Details, evidence and fixes follow. The development-grid run in progress does not need to be restarted for the blocker (it concerns test patches only), see section D.

Index of findings: 1 BLOCKER coverage of the test patches; 2 MAJOR the selection and report code has no numeric test (20 of 24 mutants survive); 3 MAJOR photutils under-tuned, threshold axis inert; 4 MAJOR DAOPHOT row aligned with test-patch truth; 5 MAJOR objective blind to false positives fainter than 22 mag, kappa at and beyond the grid edge, BIC row; 6 MAJOR freeze and provenance; 7 MAJOR (claims) ablation at the frozen kappa; 8 MAJOR (claims) statistics and pooling; 9 MINOR truth and calibration; 10 MINOR algorithm notes. Sections: A findings, B verified, C over-claim risk, D order of work, E remaining degrees of freedom, F not checked, G reproduction.

---

## A. Findings

### 1. BLOCKER. Test patches (410, 630) and (410, 730) lie largely outside the HUGS footprint; every detection there is a false positive by construction

Code at fault: `exp_astro_m2.py:39` (`TESTS`), the docstring ("(410, 630) sparse, (410, 730), (310, 730) dense"), `astro_common.score` (lines 72-113: a detection with no truth counts as a false positive, nothing models coverage).

Evidence (truth positions and pixel statistics only, no scoring). Truth sources brighter than 22 mag per 20 x 20 block (rows are y from 0, columns x from 0):

```
(410,630)             (410,730)             (310,630) [benchmark]   (310,730)
14  0  0  0  0        34 31  0  0  0         all 25 blocks >= 61*    all 25 blocks >= 109*
20  0  0  0  0        52 38  6  0  0         (*: counts of m<24)     (*: counts of m<24)
31  2  0  0  0        58 44 13  0  0
35 15  0  0  0        69 47 39  0  0
35 19  0  0  0        70 53 52  0  0
```

* 17 of 25 blocks of (410, 630) and 11 of 25 of (410, 730) have no truth source brighter than 22 mag (28 of the 100 blocks of the pooled analysis). The truth edge runs diagonally: the largest truth x in the five 20-row bands is 10.8, 15.6, 21.3, 26.1, 27.4 px for (410, 630) and 35.9, 44.4, 48.2, 55.0, 59.1 px for (410, 730) (patch coordinates).
* The pixels in those blocks contain stars: the number of 5-sigma local maxima (3 x 3 maximum, pixel more than 5 sigma above a 4 px Gaussian-smoothed copy, sigma from the Poisson-plus-dark model) per 20 x 20 block is 4 to 14 in the zero-truth blocks and 5 to 15 in the covered ones (same density), and the median DN per block is 186 to 302 in the zero-truth blocks of (410, 630) against 210 to 321 in its covered column. So the missing truth is a coverage edge, not an absence of stars. "Sparse" for (410, 630) is wrong: its 148 scored sources below 22 mag (against 969 in the benchmark patch) are the covered strip only; the label comes from the distance to the cluster core, not from the data.
* Truth quality near the edge: the fraction of truth sources brighter than 22 mag found in a single HUGS exposure (`nf606` <= 1) is 41 % in (410, 630) and 21 % in (410, 730), against 4 % in (310, 630), 7 % in (210, 630), 6 % in (310, 730) and 12 % in the dev patch (110, 630) (whose lower-left blocks (0,0) and (0,1) also have only 15 and 18 truth sources of m<24 against 33 to 101 in their neighbouring blocks: a milder case of the same thing).
* Effect on the metric, simulated on dev patch (110, 630) by keeping truth only for x <= 27 px (the (410, 630) situation; catalogues unchanged): path kappa 3, purity 0.591 -> 0.130, F 0.557 -> 0.204; path kappa 5, F 0.537 -> 0.212; crowdsource F 0.304 -> 0.154; photutils (default) purity 0.488 -> 0.151, F 0.140 -> 0.124. The paired purity difference path minus photutils changes sign (+0.10 -> -0.02), the F gap shrinks from 0.42 to 0.08. With truth kept for x <= 60 the path's F is 0.397 against 0.246 (crowdsource) and 0.149 (photutils). The methods that detect more sources are penalised most. The verdicts on these two patches, the "wins / ties / losses" record and the pooled interval (which weights blocks equally) would reflect where the HST footprint ends.

Fix (decide before the run, from truth geometry only, no detection result involved): replace (410, 630) and (410, 730) by patches whose 25 blocks are all covered, or mask both truth and detections to covered blocks (a footprint polygon from the CAOM `s_region` of the F606W image, or a rule such as "every 20 x 20 block has at least 35 truth sources of m<24") and say so in the protocol. Geometry-only candidates with all blocks covered and a density near the benchmark's (scored sources brighter than 22 mag): (410, 930) 1057, (10, 730) 1512, (310, 1030) 1473, (10, 1130) 862, (410, 1030) 665; the denser (210, 730), (110, 730), (310, 830) have 3100 to 3700 and lie about 150 px from the cluster core. I did not inspect the pixels of these candidates for saturation, spikes or frame edges. Relabel (410, 630) honestly if it is kept.

### 2. MAJOR. The code that decides the frozen configuration and the reported verdicts has no numeric test; 20 of 24 single-point mutants survive

Code: `exp_astro_m2.py:190-216` (`choose`, `choose_baselines`), `:156` (`f_measure`), `:294-331` (`cpf`, `boot_idx`, `boot_cpf`, `boot_diff`, `verdict`), `:55` (`cs_dets`), `:250` (`part_T`); `astro_common.py:42-113`. `tests/test_astro_report.py` runs `report()` on a synthetic raw with no D keys (so `choose` is never called) and checks string needles and the absence of NaN; no test calls `choose`, `choose_baselines`, `f_measure`, `part_D`, `part_T`, `cs_dets` or `baseline_scores`.

Mutants (applied to a copy, run against the fast subset of `test_astro_common.py` + `test_astro_report.py`: 7 tests, 2 s; the slow synthetic path tests only for the one marked):

| id | mutation | result |
|---|---|---|
| A1 | `match`: position AND magnitude -> OR | killed |
| A2 | `score`: detections in the edge zone are scored | killed |
| A3 | `score`: matched truth counted in the detection's block (blk22 col 1) | SURVIVED |
| A4 | `score`: matched detections selected by truth magnitude < 22 instead of detection magnitude | SURVIVED |
| A5 | `match`: any maximum matching, no minimum total distance | SURVIVED |
| A6 | `truth_for_patch`: upper y margin 0.5 px too small | SURVIVED |
| A7 | `score`: non-positive-flux detections not dropped | SURVIVED |
| A8 | `block_of`: last block column merged into the previous | killed |
| A9 | `score`: depth of the truth total 22.5 instead of 22 | SURVIVED |
| E1 | `choose`: selects the smallest F | SURVIVED |
| E2 | `choose`: ties go to the smaller kappa | SURVIVED |
| E3 | `f_measure`: purity denominator = true22 | SURVIVED |
| E4 | `choose`: only the first development patch enters F | SURVIVED |
| E5 | `choose_baselines`: photutils picks the worst variant | SURVIVED |
| E6 | `boot_diff`: interval of (b - a) beside the point estimate (a - b) | SURVIVED |
| E7 | `verdict`: win if the upper bound is positive | SURVIVED |
| E8 | `cs_dets`: S/N filter inverted | SURVIVED |
| E9 | `boot_idx`: last block never resampled | SURVIVED |
| E10 | `boot_cpf`: 90 % interval labelled 95 % | SURVIVED |
| E11 | `cpf`: purity denominator = true22 | SURVIVED |
| E12 | `part_T`: no tuned baselines | SURVIVED |
| E13 | `run_path`: clean-up threshold halved | passes the fast tests; killed only by the slow synthetic path test, by one false positive (6 > 5) |
| E14 | `part_T`: stored default catalogue at the tuned kappa | SURVIVED |
| E15 | `REXS` enlarged without a D part | SURVIVED |

(A10, an equivalent control, survived as it must.) 24 real mutants, 4 killed, 20 survived. A sign or direction error in `choose` would freeze the worst configuration without any test failing; one in `boot_diff` or `verdict` would flip the reported verdicts.

Fix: (a) unit tests with planted raw arrays: `choose` returns the planted maximum, the documented tie rule, and uses both dev patches; `choose_baselines` likewise; `f_measure`/`cpf` on hand-computed totals; `boot_diff` sign (planted a > b gives a positive interval and a "win"), `boot_idx` covers every block, `verdict` truth table, `cs_dets` filter on a planted catalogue; (b) a golden test of `score` with one case per surviving scoring mutant (A3: a matched pair that straddles a block border; A4: a detection of 21.9 mag matched to a truth of 22.2; A5: an assignment where the nearest-first choice differs from the minimum-total-distance one; A6: a truth at y = 97.2; A7: a negative-flux detection; A9: a 22.3 mag truth); (c) a test of `part_T` on a tiny synthetic patch pair.

### 3. MAJOR. The photutils baseline is under-tuned, half of its tuning space is inert, and its fluxes are biased; "the same tuning effort" is not supported

Code: `exp_astro_baselines.py:32` (`PHOT_VARIANTS`), `:54-82` (`run_photutils`: global background = 1st percentile of the patch, line 74; `DAOStarFinder(threshold = nsigma x skysig)` line 76; `IterativePSFPhotometry(..., grouper=None, mode='new', fit_shape=(7, 7))` line 78); `exp_astro_m2.py:46-48`.

Evidence (dev patches (110, 630) + (210, 630) pooled, repo scorer, F of completeness and purity below 22 mag):

* The threshold axis is a no-op. `DAOStarFinder` alone on (110, 630) returns 95 sources at 2, 3, 5 and 8 sigma. In `baseline_catalogs.npz` the catalogues of `*_t3_*` and `*_t5_*` are array-identical (x, y, flux) for five of the six pairs on the three patches I compared (110, 210, 310; array comparison, no scoring) and, for the sixth (local, 6 iterations), have identical F and counts on the dev patches. The 12 variants are 6 distinct catalogues per patch. Reason: the patch minus its 1st percentile is far above the threshold almost everywhere (median 339 and 380 DN against 1st percentiles 201 and 204 DN on the dev patches, i.e. 21 and 27 times the 6.4 DN sky sigma; thresholds of 2 to 8 sigma are 13 to 51 DN), so the threshold never binds and the default sharpness/roundness cuts decide.
* The "tuning" gain is 0.011: default 0.1292, best of the 12 ("tuned", `photutils_local_t5_i6`, tie with `t3_i6`) 0.1398.
* Standard photutils options outside the grid, all on the same dev patches (`work/phot_alt.py`): Background2D + 3 sigma, 3 iterations: 0.131; relaxed finder cuts (sharp 0 to 2, round -2 to 2) + 3 px `SourceGrouper` + Background2D + 2 sigma + 10 iterations + 5 x 5 fit: 0.158; `LocalBackground(6, 10 px)` in the fit + relaxed cuts + Background2D + 3 sigma + 10 iterations: 0.194; a low-percentile-filter background + relaxed cuts + 4.4 px grouper + 6 iterations: 0.199; `mode='all'` (DAOPHOT-style refit of all stars after each iteration, requires a grouper) + `LocalBackground` + relaxed cuts + 6 iterations: **0.244** (C 0.159, P 0.524; 17 min wall on the loaded machine). A grouper alone changes nothing (194/181 sources either way). So a reviewer's photutils scores about 1.7 times the "tuned" row. It stays below crowdsource (0.291) and far below the path (0.54), so the ordering is not at risk; the "same effort" statement and the "tuned photutils" label are.
* Photometric bias: for position-matched pairs the median (detected - HUGS) magnitude is -0.34 and -0.41 for the global-background photutils (35 % and 43 % of its position matches fail the 0.5 mag test) against -0.04 and -0.12 for crowdsource, -0.04 to -0.13 for the shipped DAOPHOT (benchmark patch) and within 0.04 mag of zero for the path at 15 to 21 mag. Plausible causes (not isolated): a background 140 to 175 DN too low (median minus 1st percentile) and neighbour flux in 7 x 7 fits; the Background2D variants shrink the bias to -0.14..+0.02.
* Crowdsource is treated better: the 48 post-filters cannot improve on the unfiltered catalogue (F = 0.2907 for six names; the "tuned" and "default" rows are the same catalogue), but its own knobs do not help either (threshold 3: 0.291; `maxiter` 20: 0.291; `blendthreshu` 1: 0.276; `refit_psf` True: 0.296; sky input pipeline/1st percentile/0: identical output, it fits its own sky; the path's smoothed-expectation weights instead of observed-count weights: 0.303). `nskyx/nskyy` >= 1 crashes under NumPy 2.5.4 (ValueError in `build_sparse_matrix`), so a sky grid could not be tried.
* The path's grid is 3 x 2 x 6 = 36 configurations (6 runs); photutils has 6 distinct catalogues; crowdsource has one catalogue and 48 selections of it.

Fix: either extend the photutils family on the development patches only (add `localbkg_estimator`, relaxed cuts, `mode='all'` with a grouper, more iterations; drop or repair the inert threshold; this needs only the baseline script and a dev re-score, no truth enters the baseline) or state in the protocol that photutils is the weaker, default-like baseline and that crowdsource and the aligned DAOPHOT catalogue are the reference baselines. Report the dev F of every family next to the path's.

### 4. MAJOR. The DAOPHOT row is aligned with the test patch's truth (an oracle shift); as shipped it scores zero

Code: `exp_astro_baselines.py:31` (`DAOPHOT_YSHIFT = -1.0`, docstring lines 12-16: "303 of its 356 sources lie within 0.5 px of a truth source at that shift"). The shift was chosen on the benchmark patch against the HUGS truth that scores it; there is no dev-only way to fix it (the catalogue exists for (310, 630) alone).

Evidence (the permitted exception, DAOPHOT on (310, 630) only): as shipped (dy = 0) 5 of 329 detections below 22 mag match (F 0.008); dy = -0.5: F 0.129; -0.75: 0.307; **-1.00: C 0.233, P 0.689, F 0.348**; -1.06 (the histogram peak, dx -0.04): 0.352; -1.25: 0.322; -1.5: 0.194. Plateau from about -0.95 to -1.10, so -1.0 is neither fragile nor finely tuned. The catalogue is 356 sources, all brighter than 21.4 mag. The unexplained part is a y-only offset of 1 px (a 1-based to 0-based conversion would move x too; `m2_field_spec.md` item 10 lists two other conflicting offsets).

Consequence: the aligned DAOPHOT (F 0.35) is the strongest baseline on the benchmark patch, stronger than the dev F of crowdsource (0.29) and photutils (0.14) (different patches, so indicative only). The alignment favours the baseline, so it is conservative for a "path beats the baselines" claim, but the row is not "the published result": under the same truth the unshifted catalogue scores F 0.008.

Fix: label the row "DAOPHOT catalogue, y shifted by -1 px, shift fitted to HUGS on this patch"; report the unshifted score beside it; the recorded `-1.0` is fine, the docstring should say that the value used the test-patch truth.

### 5. MAJOR. The tuning objective cannot see false positives fainter than 22 mag, so kappa is pushed to and beyond the lower edge of the grid; the BIC row depends on the path length

Code: `exp_astro_m2.py:156-160` (`f_measure`: purity = det22_matched / det22, only detections whose own magnitude is below 22), `:190-200` (`choose`), `:42` (`KAPPAS`), `:109` (`tau_max = 1 / (min(kappas) x scale)`), `:117-122` (BIC readout).

Evidence (path, sigma_b 5, r_ex 1.2, dev patch (110, 630); my own runs with the repo's settings; for the lower kappas I ran the path to the kappa = 2 time, 1516 entries, 496 s):

| kappa | 2 | 2.5 | 3 | 3.5 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|---|
| completeness | 0.563 | 0.548 | 0.526 | 0.504 | 0.492 | 0.475 | 0.449 | 0.387 |
| purity (< 22) | 0.590 | 0.591 | 0.591 | 0.585 | 0.591 | 0.616 | 0.622 | 0.614 |
| **F** | **0.576** | 0.569 | 0.557 | 0.542 | 0.537 | 0.537 | 0.522 | 0.474 |

* F rises monotonically as kappa falls and the lowest value of the grid is 3: the optimum lies outside the grid on (110, 630). On (210, 630), with the path also run to the kappa = 2 time (1627 entries, 546 s), F for kappa 2, 2.5, 3, 3.5, 4, 5, 6, 8 is 0.537, 0.526, 0.530, 0.526, 0.527, 0.527, 0.514, 0.489 (flat within 0.01 for kappa 2 to 5, again highest at the lowest value). Pooled over the two dev patches: 0.5549, 0.5456, 0.5421, 0.5335, 0.5316, 0.5312, 0.5173, 0.4823. So "tuned kappa" will be 3.0 for this slice, with the optimum below the grid.
* The objective counts false detections only among those whose own magnitude is below 22. Detections between 22 and 23.5 mag on (110, 630), matched / total (purity), at kappa 2, 3, 5, 8: 35/195 (0.18), 28/118 (0.24), 9/30 (0.30), 0/0; unmatched detections of all magnitudes inside the scored area: 563, 461, 323, 250. From kappa 5 to kappa 2 the F gains 0.039 while the unmatched detections rise by 240 (74 %). The stopping rule is therefore not tuned for false-positive control: with this objective the optimum is the most permissive stop. The same blindness applies to the tuned baselines. The benchmark's own FDR curve runs to 23.5 mag and will show the cost in the faint bins on the test patch.
* The BIC row depends on `min(KAPPAS)`. On (110, 630) the BIC minimum is interior (block 1291 of 1516 on the path extended to the kappa = 2 time, 1227 units; cost 51011 there, 52428 at the end), essentially the kappa = 3 block (1233 units) with the clean-up at 5.26: F 0.550. On (210, 630) it is not: on the path stopped at the kappa = 3 time (as `run_path` does) the minimum is the last block but one (block 1403 of 1404, 1342 units, cost 495541 against 495549 at the last block); on the path run to the kappa = 2 time it moves to block 1558 of 1627 (1478 units, cost 494785 against 495129 at the end). The F of the row changes little (0.530 against 0.533 on (210, 630); pooled 0.5392 against 0.5408), but the row is "the BIC minimum within the path as truncated at the smallest kappa of the grid", not "nothing to tune"; on a denser patch the truncation will bind harder.

Consequences: the tuned (kappa 3, or 2 if the grid is extended) and default (kappa 5) rows differ by 0.011 (3 against 5) to 0.024 (2 against 5) F pooled on dev; the headline gaps on dev (0.54 to 0.58 against 0.29 against 0.14) do not depend on that choice, but the text must not present kappa as tuned to a false-positive-controlled optimum.

Fix: (a) the dev counts needed for a false-positive-aware objective are already in the D raw (`tot` holds det_all and matched_all, so unmatched detections = det_all - matched_all for every (sigma_b, r_ex, kappa)): decide now, before the D result is read, whether the objective stays F at mag < 22 or gets a term or constraint on the total unmatched detections (or the benchmark's FDR over 15 to 23.5 mag), and apply it identically to the path and to the baselines; report total unmatched detections in every main table; (b) extend `KAPPAS` to 2 or lower, or state that the optimum is on the edge (extending changes D and needs a re-run); (c) run the path past the BIC minimum (for example to the kappa = 2 time) and store whether the minimum is interior (`index < len(hist) - 1`), or call the row "BIC within the path truncated at the kappa = 3 time".

### 6. MAJOR. The freeze is not mechanical and the provenance of what produced D cannot be shown

* The D run started at 20:27:18; `exp_astro_m2.py` was modified at 20:30:55 (the BIC readout and report lines), `exp_astro_baselines.py` at 20:28:05, `baseline_catalogs.npz` regenerated at 20:33:29 and `part_D` reads it only after its path jobs (`pool.map`, line 169; `baseline_scores` line 177). The edits do not change the path outputs the selection reads (I diffed against the 20:17 backup), but nothing records which file text produced D, and T will import different text. `scratchpad/run_queue3.log` shows the baseline step crashing on the sparse patch (`Background2D`, `exclude_percentile`) and the queue going on to D in the same second (20:27:18); the baseline script was then edited (20:28:05) and the catalogues regenerated (20:33:29).
* `baseline_catalogs.npz` is untracked (`git status`: `??`) and has no recorded hash; the test and dev baseline catalogues are inputs of the result.
* `part_D` stores only counts and block counts, not catalogues; any change of the scoring (a coverage mask, finding 1) forces a re-run of all 12 path jobs (505 and 566 s each for my two runs under load 14).
* `part_T` overwrites its keys (`save_part`), nothing enforces "run once", and the raw does not record the commit or file hashes (only `cfg = [hp, rex, kappa]`).
* `DEFAULT = (5.0, 1.2, 5.0)` is described as "chosen on no patch"; the first smoke run of the path (`scratchpad/astro_smoke.py`, log 18:54) used exactly (5.0, 1.2) on dev patch (110, 630) and printed F for kappa 3 to 8. That is a dev patch, so it is not a leak, but "chosen on no patch" cannot be verified; say "set before the grid, first run on (110, 630)".

Fix: before releasing T, write `freeze.json` (git commit, sha256 of `imaging.py`, `exp_astro_m2.py`, `astro_common.py`, `exp_astro_baselines.py`, `baseline_catalogs.npz`, `truth_hugs_meth1.npz`, the chosen configuration and the tuned baseline names) from the D result; make `part_T` refuse to run if a recorded hash differs or if any `T|` key exists; commit the baseline catalogues (or their hash); save the D catalogues (optional: it needs a D re-run); add `n_entries == MAX_ENTRIES` to the raw and the report (the cap, 2500, is likely not binding: entries were 1300 and 1404 at densities 0.114 and 0.139 per px^2, and the densest test patch has 0.228, so perhaps 1700 to 2300, but that is an extrapolation).

### 7. MAJOR (for claims). The design ablation (part V) is run on the benchmark patch at the kappa tuned for the full method

Code: `exp_astro_m2.py:219-247`. `greedy`, `nohp`, `norex`, `snr`, `nocleanup` are evaluated at the frozen (sigma_b, r_ex, kappa) without re-tuning kappa per variant (`_variant_job` passes the frozen `hp`, `rex`, and the report uses the same `kap`). In the SMLM experiments the threshold is re-chosen on development frames for each variant. Removing the high-pass changes the noise scale and hence what kappa means; the `snr` variant uses a kappa chosen for the time rule. A drop in F for a variant measures "this component at a kappa tuned without it", and it is on a single patch. The text must not say "component X is necessary" from it.

Fix: re-choose kappa per variant on the dev patches (cheap in code, 5 more dev paths), or label the table "at the frozen kappa, one patch, descriptive".

### 8. MAJOR (for claims). Statistics: single-method intervals are optimistic, paired intervals are fine, the pooled analysis is dominated by one patch and its units are not independent

Code: `exp_astro_m2.py:302-331` (bootstrap), `:456-490` (pooled analysis).

* Coverage of the 25-block percentile interval, simulated (1500 repetitions, 25 blocks drawn with replacement from the 50 dev-patch blocks of the crowdsource catalogue, B = 1000): 93.2 %, 93.3 %, 93.7 % for completeness, purity, F against the pool value; 92.1 %, 92.9 %, 93.3 % for the paired crowdsource minus photutils differences. That is in the favourable iid-block model.
* The blocks are not independent: lag-1 (rook) correlation of block completeness residuals is 0.36 and 0.40 on the two dev patches. Block-size sensitivity of the bootstrap standard error (the two dev patches, kappa 3): F of the path 0.0123 (10 px), 0.0182 (20 px), 0.0197 (25 px), 0.0275 (50 px, 8 blocks); crowdsource 0.0110, 0.0127, 0.0148, 0.0178; the paired difference path minus crowdsource 0.0140, 0.0160, 0.0157, 0.0154. So single-method intervals at 20 px are too narrow by roughly 10 to 50 %; the paired differences are stable, which is what the verdicts use. The verdicts will not change; the quoted per-method intervals should be widened or computed with 40 px blocks or a patch-level jackknife.
* The blocks assign a truth source and its match to the truth's block and a detection and its match to the detection's block (`astro_common.py:105-110`); a matched pair that straddles a border is split between blocks. That is a valid ratio-of-sums cluster bootstrap but untested (mutant A3).
* The four test patches tile a contiguous 200 x 200 region of one frame, so the "record by patch" (wins / ties / losses) is not four independent trials; the dev patch (210, 630) shares its border with the benchmark patch (310, 630) (no pixel overlap, the 2.5 px margins make a dead strip; the docstring's "non-overlapping" is true, "independent" is not). Scored truth sources below 22 mag per test patch: 969 (benchmark), 148, 541, 2057; the pooled ratio-of-sums is 55 % one patch ((310, 730), a density 1.6 to 2.0 times that of the dev patches, hence an extrapolation, not a held-out replicate of the tuning set) and 26 % the benchmark patch, and the pooled blocks contain the zero-truth blocks of finding 1. Stratify the pooled bootstrap by patch and report per-patch values as the primary.
* Verdicts are on F alone, no multiplicity control (about 8 comparisons per patch); the differences are large so this does not matter, but the text should not read a significance level into them.

### 9. MINOR. Truth and calibration details

* Truth list: all HUGS method-1 sources with a valid F606W (`data/astro/make_truth.py`), no F814W requirement (the shipped metric requires both bands and a margin of 2.5 px in pixel-centre coordinates), no `nf606` cut, WCS mapping only. Counts below 22 mag in the scored area: 969 for (310, 630) against 951 under the shipped selection and the paper's 1,049. Method 1 against method 2 is undocumented as a choice (1,070 against 1,071).
* Common systematic position offset: the median (detected - truth) position of bright matched pairs is (+0.15, +0.02) px on (110, 630) and (+0.10, +0.03) on (210, 630) for the path, (+0.15, +0.06) and (+0.09, +0.03) for photutils, (+0.15, +0.03) and (+0.10, +0.03) for crowdsource (a 14-yr epoch difference between SDSS 2001 and HUGS 2015 is a plausible partial cause, not checked). Shifting the truth by +0.15 px in x changes F on (110, 630) by +0.005 (path kappa 3), +0.009 (kappa 5), +0.012 (crowdsource), +0.005 (photutils): immaterial and common to all methods.
* One calibration constant (`NMGY_PER_COUNT = 0.00546689`) is used for every method; the per-column calibration of the six patches is 0.981 to 0.984 of it (a uniform offset of about 0.02 mag, spread about 0.002 mag), the same for every method.
* `KAPPAS`, `HPS`, `REXS`, `LAT_H = 0.5`, `MARGIN = 2.0`, `FIT_RADIUS = 4`, `MAX_ENTRIES = 2500`, the PHOT and CS grids are unchanged since the first commit of the script (5527d7b); `DEFAULT` appears in 6d4dc32 and is unchanged. The three-part PSF mixture (fit rms 0.0010 of the peak) is fitted per patch from its psField stamp.

### 10. MINOR. Algorithm notes (nothing found that changes a number on this data; each is a statement the text should respect)

* Weights `sw` come from the smoothed data (`exp_astro_m2.py:102-103`); the high-pass leaves the noise slightly correlated with variance sigma^2 (1 - 2 S_ii + sum S_ij^2) = 0.973, 0.990, 0.998 sigma^2 for sigma_b 3, 5, 8: immaterial.
* `unit_significance` (the clean-up at kappa and the BIC clean-up at 5.26) is the unit's amplitude over the standard error it would have standing alone. In this field 40 to 48 % of the m<22 truth sources have an m<22 neighbour within 1.2 px (18 to 25 % within 0.8 px; 70 to 74 % have an m<24 neighbour within 1.2 px; dev patches), so the true standard error of a blended unit is larger. kappa is a threshold of the procedure, not a significance level.
* `run_path` returns the blocks of `sliding_path` as they are (local refits of the units within 4 px of each entrant, `passes=1`, `tol=1e-6`); the global polish `refine` that `fit_frame` uses is not called, so catalogue fluxes are not at the joint least-squares optimum. Consistent between D and T; say so.
* The local windows are +-10 px; the high-pass partner of a kernel has width about 6 px at sigma_b = 5, so its tail beyond 10 px is cut in the local fits: tiny.
* The lattice margin is 2 px; stars 3 to 10 px outside the patch put wings into the border rows and have no kernel, and the scoring drops only 2.5 px. The contamination of the 2.5 to about 6 px ring is the same for all methods (none models outside stars); a larger scored margin would remove it.
* The exclusion zone (0.8 or 1.2 px) forbids a second unit within that distance of an active one, while the truth has an m<22 neighbour within 0.8 px for 18 to 25 % and within 1.2 px for 40 to 48 % of the m<22 sources; with one-to-one matching a detection matches one truth source, and at the SDSS resolution (FWHM 2.2 px) such pairs are in practice unresolvable. A completeness near 0.5 is therefore not by itself a sign of a defective method (I did not quantify the ceiling).
* The path is deterministic (no seeds); my run reproduced the earlier smoke log to all digits (B.10).

---

## B. Verified (re-derived independently and agrees)

0. **Provenance of every choice that enters a test number** (what I could prove from the code, the git history, the scratch logs and the file times):

| choice | fixed on | evidence |
|---|---|---|
| sigma_b, r_ex, kappa of the tuned path | the development patches (110, 630), (210, 630) | `choose` reads only `D|` keys of `raw` (`exp_astro_m2.py:190-200`); no other input |
| form of the time rule, clean-up at kappa, `unit_significance` | design, before D | commits 5527d7b / 7993d4d; exercised on dev patch (110, 630) in the smoke runs only (`scratchpad/astro_smoke.py`) |
| lattice 0.5 px with 2 px margin, fit radius 4 px, `MAX_ENTRIES` 2500, three-part PSF mixture, variance smoothing 1.5 px, 1st-percentile floor, grids `KAPPAS` / `HPS` / `REXS` | design, before D, unchanged since 5527d7b | no log of a path run on a test patch exists in the scratch and queue logs (`run_queue4.sh` holds T behind a release flag); unlogged runs cannot be excluded |
| `DEFAULT` = (5.0, 1.2, 5.0) | no selection; first run on dev (110, 630) | finding 6 |
| photutils / crowdsource variants; which one is "tuned" | variants: design before D; pick: argmax of F on the dev patches (`choose_baselines`, `exp_astro_m2.py:207-216`) | finding 3 |
| baseline catalogues of the test patches | computed 20:33 without any truth; no selection depends on them | `exp_astro_baselines.py` reads no truth |
| crowdsource x/y swap | benchmark-patch truth (`data_access.md`: 380 against 171 of 563), re-derived by me on dev patches only | item 7 |
| DAOPHOT y shift of -1 px | **the test patch's own truth** | finding 4 |
| truth list and scoring rules (0.5 px, 0.5 mag, 2.5 px, depth 22, HUGS method 1, no F814W condition) | the benchmark's definition plus the authors' choices | finding 9 |
| the set of test patches | the authors' guess of crowding level from the distance to the core; coverage not considered | finding 1 |

1. **Scoring.** I re-implemented the benchmark from `data_access.md` and `m2_field_spec.md` without reading `astro_common.score` first (`work/myscore.py`: truth patch box [x0-0.5, x0+99.5), the 2.5 px margin from the box edge, F606W as r, `22.5 - 2.5 log10(flux x 0.00546689)`, 0.5 px and 0.5 mag, one-to-one maximum matching through scipy `maximum_bipartite_matching`). Compared with `astro_common.score` on the 26 stored dev-patch baseline catalogues: the four totals (true22, matched22, det22, det22_matched) agree exactly in 25; in one (crowdsource, (110, 630)) matched22 is 228 against 229, a tie-break between maximum matchings (the repo takes the minimum total distance among them). On the path's own dev catalogues (kappa 3, 5, 8 on both patches): 5 of 6 exact, one 558 against 559 (same cause). On 400 random synthetic cases (random truth, jittered copies, duplicates, fakes, detections outside the patch): true22, det22 and matched22 identical in all 400, det22_matched differs by 1 in 2. No bug in the matching, the margins, the magnitude conversion or the blocks. Detections within 2.5 px of the box edge are dropped from numerator and denominator for every method, as are truth sources there, so the path's 2 px lattice margin brings no advantage.
2. **Many-to-one against one-to-one.** The shipped metric is a ball query, the repo's is one-to-one: for the path on (110, 630) completeness changes by +0.004 and purity by 0 when one-to-one is dropped (kappa 3, 5, 8); for the baselines by at most 0.01. A faint truth source cannot "steal" a detection from a bright one: completeness at 22 mag is identical (540/1027, 488/1027, 229/1027, 84/1027 for path 3, path 5, crowdsource, photutils) whether the assignment uses all truth magnitudes or only those below 22.
3. **Same magnitude for every method.** All methods go through `counts_to_mag(flux in DN)`; PSF stamps are unit sum for photutils and crowdsource; the path's mixture is a unit-sum fit of the same stamp. The F606W-as-r treatment shows no net colour or zero-point offset for the path (median dmag -0.04 to +0.04 for truth 15 to 21 mag, 16 to 84 % range +-0.09 for 15 to 19 mag, -0.11 at 21 to 22 mag on (110, 630)); the offsets of the baselines are the ones in finding 3, not a colour term. The F gap does not hinge on the 0.5 mag tolerance: with tolerance 0.5, 1.0 and unlimited, F on (110, 630) is 0.557, 0.649, 0.675 (path kappa 3), 0.304, 0.355, 0.369 (crowdsource), 0.140, 0.182, 0.199 (photutils).
4. **Geometry.** The patch convention (column = x, row = y, 0-based, centres at integers, box [x0-0.5, x0+99.5)) is confirmed by the data: detections of all three methods sit within 0.15 px of the truth.
5. **High-pass consistency.** For sigma_b 3, 5, 8 and stars inside, on and 2 px outside the border, the high-pass of a star frame equals the high-passed kernels to 3e-14, the windowed local model equals the full model exactly, and the analytic Jacobian matches finite differences to 3e-10.
6. **Protocol reads.** `choose` and `choose_baselines` read only `D|` keys and select by the F of completeness and purity below 22 mag over the two dev patches' pooled counts, the maximum by (rounded F, kappa) or (rounded F, name); they do not look at test numbers. `part_T` calls them before touching a test patch, reads the test patch's own pixels and truth only inside the loop, saves only the default and the chosen baseline variants, and computes but discards the scores of the other variants (`baseline_scores` scores all 60 on the test patch; nothing is printed or stored). The test pixels the path reads are the patch's own (percentile, median, smoothing), no truth.
7. **Patches.** Dev (110, 630), (210, 630) and the four test patches are disjoint (zero pixel overlap) and contiguous: the dev strip x 110-309 continues into the test strip x 310-509 at rows 630-729, and the second test row (rows 730-829) is above it. The crowdsource x/y swap needed in `run_crowdsource` is confirmed on the dev patches (377 of 564 and 415 of 572 sources within 0.5 px of a truth source swapped, 216 and 244 not swapped), so no test-patch truth is needed for it.
8. **Baseline inputs.** Both baselines receive the exact 51 x 51 psField stamp (a better PSF than the path's three-Gaussian fit), the same gain/dark noise (from the observed counts, not the path's smoothed expectation: worth +0.012 F to crowdsource), the same patch and border. Crowdsource ignores the input sky level (identical output for pipeline sky, 1st percentile and zero). The default crowdsource filter drops only non-positive-flux rows (0 on the dev patches, 1 on three test patches). What the baselines do not have is the path's spatially adaptive background (the high-pass): crowdsource fits its own local sky, photutils has a constant (global) or box-median (local) background. That is a difference in capability, not in information, and it is part of what the comparison measures.
9. **Tests.** `pytest tests/test_imaging.py`: 42 passed (86 s); `tests/test_astro_common.py` + `tests/test_astro_report.py`: 9 passed (286 s under load 12).
10. **Reproducibility of the path.** My run of the path on (110, 630) at (5.0, 1.2) reproduces to all printed digits the figures in `scratchpad/astro_smoke2.log` (kappa 3: C 0.526, P 0.591, F 0.557; kappa 8: 0.387, 0.614, 0.474).

---

## C. Over-claim risk: what the paper must not say

* Not "outperforms DAOPHOT / PCAT on the benchmark of Portillo et al.": truth is HUGS (not the ACS-survey catalogue the papers used, which could not be downloaded), the match is one-to-one at 0.5 px (the first paper uses 0.75 px, the follow-up 0.5 px, both with a many-to-one ball query in the shipped code), the truth selection differs (margin 2.5 px from the box edge instead of the script's 2.5 < x < 96.5, no F814W condition), PCAT was not run, and the DAOPHOT row carries a shift fitted to the truth on the same patch.
* Not "beats the standard tools": two off-the-shelf packages, one of which (photutils) is demonstrably under-tuned (dev F 0.14 reported, 0.24 reachable with documented options); no DAOPHOT/ALLFRAME run, no DOLPHOT, no PCAT. Say "against photutils and crowdsource at the settings of the protocol".
* Not "generalises": one globular cluster, one SDSS frame, one band, one PSF family; tuning and testing patches are neighbouring windows of the same frame (the dev patch (210, 630) touches the benchmark patch); four test patches, contiguous, one of them 55 % of the pooled truth, and only two of them, (310, 630) and (310, 730), are fully covered by the truth as the set stands.
* Not "5-sigma detections" or "false-alarm rate": kappa is a threshold of a time rule and of the unit's stand-alone amplitude significance (`unit_significance`, which ignores the covariance with neighbouring units); there is no calibration of the false-detection rate. Purity is not 1 - false-detection rate: the HUGS list has a footprint edge, matches have about 4 % chance coincidence at 22 mag (density of truth within +-0.5 mag times the 0.785 px^2 area), and 20 % (20 to 21 mag) to 32 % (21 to 22 mag) of the path's position matches fail the 0.5 mag test (path, dev (110, 630), kappa 4; at kappa 3 the position-only purity is 0.78 against 0.59).
* Not "false-positive-controlled stopping" or "tuned to an optimum" (finding 5: the objective ignores false detections fainter than 22 mag and the optimum is outside the grid), not "nothing to tune" for the BIC row (its block moves from 1342 to 1478 units on (210, 630) when the path is run longer), not "chosen on no patch" for `DEFAULT`.
* Not "each component is necessary" from part V (finding 7).
* Not a significance level from the intervals (finding 8), and not "four independent test patches".
* Not a speed statement without the cost: the path takes 505 to 566 s for a 100 x 100 patch on one loaded core (1300 to 1400 entries to the kappa = 3 time); photutils and crowdsource take seconds to tens of seconds.
* The faint end: HUGS has 1 to 5 exposures per source and a rollover beyond about 27 mag; at 22 mag it is safe, but 12 % of the dev (110, 630) truth below 22 mag is single-exposure; run the sensitivity with `nf606` >= 2.

---

## D. What can proceed, and the order of work before the release

* The running D grid reads dev patches only; findings 1, 2, 4, 6 and 7 do not require restarting it. Findings that would: an extension of `KAPPAS` below 3 (finding 5b), a coverage mask applied to the dev patches (their blocks (0,0) and (0,1) of (110, 630)), or saving catalogues. The cheapest consistent route is to leave D running, replace the two test patches (no mask needed on the dev patches; document that two dev blocks are partly uncovered), decide the objective (finding 5a) from the D raw before looking at it, add the baseline variants as a separate part scored on the dev patches, and add the tests and the freeze record.
* Order: (1) fix the test set; (2) tests for `choose`/`choose_baselines`/`boot_*`/`verdict`/`score` cases; (3) freeze record and write-once guard; (4) photutils family or relabelling; (5) decide the objective (finding 5a) before D is read; relabel DAOPHOT, DEFAULT, V; (6) kappa extension if the objective stays.

## E. Degrees of freedom a sceptic would still list after the above

(i) The choice of the F-measure of completeness and purity below 22 mag as objective and as the headline (equal weights; the benchmark reports completeness and FDR by bin). (ii) Which patches are tests, and that they were picked before anyone knew the coverage. (iii) The path's constants fixed by design before D (lattice 0.5 px with a 2 px margin, fit radius 4 px, three-part PSF mixture, smoothing of the variance, `rho`, `box`, tolerances, the clean-up at the same kappa, `tau_max` from the kappa grid), partly explored on dev patch (110, 630) in smoke runs; the grids themselves (3 x 2 x 6). (iv) The baseline families and the objective used to tune them (finding 3). (v) The truth list (method 1, no F814W condition, no exposure cut, WCS-only mapping, Vega/AB question for F606W left open, F606W as r). (vi) The 2.5 px edge rule and the symmetric box-edge definition (0.5 px different from the shipped script). (vii) The DAOPHOT shift. (viii) The pooled analysis' unit (blocks of all patches exchangeable). (ix) Presentation: default-against-default, tuned-against-tuned and BIC-against-all rows are all reported, which is good; the verdict is on F alone.

## F. Things I could not check

* The result of the development grid (still running at 21:41 after 1 h 14 min; the four workers get about 13 % CPU each under load 14 to 18). I ran only the (sigma_b 5, r_ex 1.2) slice, on both dev patches, with the path also taken to the kappa = 2 time, so I cannot say which configuration will be frozen, whether the optimum sits on another grid edge in sigma_b or r_ex, or how flat the surface is there.
* Anything about the path on a test patch, by rule; hence the possible effect of the 2500-entry cap on (310, 730), and whether the dense patch behaves like an extrapolation.
* The exact HST footprint (no CAOM `s_region` polygon in the data); the coverage evidence is from truth positions and pixel peaks.
* The ACS-survey catalogue and the shipped metric script on our data; PCAT/Lion; DAOPHOT proper; other versions of photutils (3.1.1) and crowdsource (0.5.9); crowdsource's sky-grid option (it crashes under NumPy 2.5.4).
* Coverage of the bootstrap under the true spatial dependence (only the iid-block simulation and the block-size sensitivity).
* The manuscript text (I found none that describes the astronomy results) and pixel-level properties (saturation, spikes, frame edges) of the suggested replacement patches.
* Mutants of `imaging.py` (not asked; its 42 tests pass).

## G. Reproduction

Scripts (all in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_IS1a/work/`, run with the repository venv, `OMP_NUM_THREADS=1`, `nice -n 10`; the photutils/crowdsource ones with the domain venv): `myscore.py`, `cmp1.py`, `cmp_path.py`, `fuzz.py` (scoring re-implementation and comparisons); `posonly.py`, `posonly2.py`, `steal.py`, `m2o.py`, `offset_sim.py` (matching, calibration, offsets); `phot_alt.py`, `cs_alt.py`, `thr_test.py`, `devF.py` (baseline variants on the dev patches); `path_dev.py`, `path_ext.py`, `analyze_path.py`, `analyze_path_ext.py`, `pooled_kappa.py`, `pooled_kappa_ext.py`, `faint_fp.py` (path runs on dev patches; pickles in the parent directory); `cov_sim.py` (coverage simulation); `boot_cov.py`, `blocksize.py` (bootstrap); `dao_shift.py` (the permitted DAOPHOT check); `mutate.py` (24 mutants, copies in `../mut/`).

