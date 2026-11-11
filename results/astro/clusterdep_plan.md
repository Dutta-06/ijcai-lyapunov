# Plan for the cluster-dependence test: a post hoc diagnosis of the spent patches, one label-free repair, a held-out test on six clusters no earlier test used (declared before any held-out catalogue of the ISS exists and before anything is scored)

**Declaration time: 2026-11-09T14:54:09Z** (`date -u`, taken immediately before this file was committed; the commit time is that of the commit that contains this file). This file, `experiments/astro_clusterdep.py`,
`experiments/exp_astro_clusterdep.py`, `experiments/exp_astro_clusterdep_baselines.py`, `experiments/clusterdep_select.py`, `data/astro/fetch_clusterdep.py`, the selection `clusterdep_selection.txt`, the registry
`clusterdep_patches.json`, the baselines' catalogues `clusterdep_baseline_catalogs.npz` and the freeze record `freeze_clusterdep.json` (hashes of all of them) are committed BEFORE the ISS is run on a held-out patch; the run
refuses to start unless the freeze record matches the files.

## 0. State at the time of writing (what exists, what does not, what I know)

* **Exists.** (i) The spent evidence: 19 test patches in 5 clusters (M2: 4 windows of the benchmark frame and 6 patches of three repeat frames; M53: 4; M92: 2; M3: 1; M15: 2) with the stored guarded ISS catalogues and the front-end
  crowdsource catalogues; the post hoc diagnosis of them (`clusterdep.md` sections 1 and 2, produced by `exp_astro_clusterdep.py --diag` before this plan: its result is described in section 1 below and is POST HOC).
  (ii) The held-out patches as geometry: the selection (`clusterdep_selection.txt`), the registry, the patches' pixels, the HUGS truth lists mapped to the frames, the descriptive rigid offset of the truth against the pixels (a fit of a
  PSF-convolved copy of the truth list to the pixels, no detector output), and the baselines' catalogues (being written by `exp_astro_clusterdep_baselines.py` from pixels and PSF stamps alone, no truth read).
* **Does not exist.** Any ISS run on a held-out patch; any score of any catalogue on a held-out patch; any comparison of a held-out detection with a held-out truth source.
* **What I know about the held-out patches** (stated, because it can bias): their number, density of F606W < 22 truth stars per pixel, stamp FWHM, descriptive rigid offsets (below), and the dry-run funnels (stars brighter than the limit, saturated
  pixels). Nothing else. I know the spent results (section 1) and chose the repair after seeing them.

## 1. The question and what the spent patches say (POST HOC)

**Question.** The lead of the ISS over crowdsource given the same high-pass front end (the strict comparison O2 of `fresh_plan.md`) was a win at M2 and M15, a tie at M3, a loss at M53 (frozen isotropic kernel) and M92. Is the cluster dependence
explained by a label-free property of the frame, and can a label-free repair reduce it? "Reduce" is tested on clusters never scored.

**Diagnosis on the 19 spent patches (`clusterdep.md` sections 1.x; all from stored catalogues, nothing re-run).**
* The gap ISS(isotropic kernel) minus crowdsource is predicted, leave-one-cluster-out, by the label-free misfit of a round Gaussian-mixture fit to the PSF stamp (rank rho = -0.82, LOCO rmse 0.096 against 0.146 for the null, signs right at 5 of 5 clusters),
  by the stamp's ellipticity (-0.84, 0.085, 4 of 5) and by the residual scale of the ISS catalogue (-0.75, 0.082, 5 of 5). Seeing (stamp FWHM), sky, noise, background structure, crowding (density, fraction of stars with a brighter neighbour within 1.5 or 2.5 FWHM)
  and the 21-22 mag completeness do not predict it (LOCO rmse at or above the null, signs 0 to 3 of 5).
* After the kernel rule (the general kernel fitted to the stamp wherever the round fit is inadequate) no covariate predicts the remaining gap in LOCO except the residual scale of the round-kernel catalogue (0.067 against 0.084, 5 of 5): the remainder is
  not explained.
* Where the ISS loses (M53 isotropic, M92): bright stars (< 19 mag) are represented by a doublet or triplet of units (2.0 to 2.1 detections within 1.5 px of a bright truth star against 1.2 for crowdsource), the flux is split so the nearest unit is off in
  magnitude by more than the 0.5 mag tolerance (bright matched fraction 0.65 at M92 against 0.83), and 62 to 66 % of the unmatched detections lie 0.5 to 1.5 px from a truth star. Misses (no detection within 1.5 px) are not the problem (5 to 12 % of the
  truth stars for every method). The positions of the bright stars are within the tolerance (nearest unit 0.36 px, crowdsource 0.37 px). Stopping too late is not the cause (the stored F of the isotropic row is nearly flat in kappa and, on both M92 patches, still below crowdsource's at kappa 8).
* Attempts that did not repair it. A variance floor proportional to the first-pass model flux (a residual-calibrated noise model, `astro_clusterdep.calibrate` and `run_path_weighted`, kept in the module) changed F by +0.00 to +0.02 on the two M92 patches in a
  scratch run (not a script of the repository; a development run of the module on the M92, M53 and M2 patches is made AFTER the held-out run and reported in `clusterdep.md` section 2b, POST HOC) and left the number of units around bright stars unchanged. In scratch checks (also not reproduced by a script) an oracle
  magnitude shift of the M92 rule catalogue changed F by at most 0.016, and an estimate of a PSF stretch from the bright isolated stars of one patch found only 0 to 14 such stars (not enough; the in-situ stamp of `insitu.md` is the stronger form of that
  attempt and was a partial repair).

**The repair declared (one).** *Merge of units that the kernel cannot tell apart*: after the path has stopped (tuned configuration, kappa 2, clean-up by the unit's own significance), two units whose kernels have normalised inner product (coherence)
at least **rho\* = 0.8** are one source: the pair with the largest coherence merges first, the merged unit sits at the flux-weighted centroid and carries the summed flux, coherences with the neighbours are recomputed, until no pair is left
(`astro_clusterdep.merge_coherent`). The coherence is that of the kernel the row uses (the pixel-integrated kernel without high-pass, from `psf_general.fit_adequate`, tabulated on offsets of +-4 px, 0.25 px steps, bilinear), so the radius
scales with the seeing and follows an anisotropic kernel; for the stamp of M2 (FWHM 2.15 px) rho\* = 0.8 is a separation of about 0.9 px, for 3.4 px seeing about 1.4 px. rho\* is the value of the grid {0.6, 0.7, 0.75, 0.8, 0.85, 0.9} with the largest mean over the
19 spent patches of F of the ISS rule row (0.387 against 0.374 unmerged and 0.323 for crowdsource): it was chosen on spent data and is therefore tested here on data never seen; the declared value is final. Its principle: the sparse solution has no
penalty on splitting a source into two units that the dictionary barely distinguishes, so when the kernel does not match the star the optimizer spends the freedom on a doublet; the merge applies the identifiability limit of the dictionary after the fact.
The development gain is small (+0.013 mean F over the 19 spent patches; M2 +0.006, M53 +0.025, M92 +0.037, M3 +0.021, M15 -0.011) and does not remove the M92 loss; the test is whether it is real.

## 2. Data (all fixed before the run)

Universe: the HUGS clusters of the coverage survey (`second_cluster_spec.md`, `fresh_survey.txt`) that no earlier test used: NGC 5053, 5466, M5 (5904), M13 (6205), M10 (6254), NGC 6535, M71 (6838), NGC 6934, M72 (6981) (M2, M53, M92, M3, M15
are excluded). The selection rule and its geometry-only outcome are in `clusterdep_selection.txt`; the code is `experiments/clusterdep_select.py` (the rule of `fresh_rule.md` with the constants below; `astro_second.py` and `astro_fresh.py` imported, not edited).

1. *Version 1* (the relaxed gates of the first dry run, declared before any dry run was read): PSF status any (status 3 = no PSF star, extrapolated PSF: kept and flagged), stamp FWHM at the frame centre in [2.0, 3.8] px (M2's own seeing is 2.1 px, the best frame
   of M15 3.75), at least 30 survey windows, the patch rules of M53 and M92 unchanged (100 x 100 px; 40 px from the frame edge; inside the HUGS footprint; no pixel above 0.8 x the saturation level within 10 px; no truth source brighter than F606W = 14.5
   within 10 px; no HUGS source with F814W < 17 and no F606W within 10 px; density 0.03 to 0.15 per px), frame = the passing frame with most windows, the next passing frame if it leaves no patch. Dry-run outcome: 2 clusters (NGC 6934 with 3 patches, NGC 5466
   with 1); M5, M10, NGC 6535, M71 lose every candidate to the bright-star rule (M5: 227 origins after the saturation rule, 0 after it) and M13 to the saturation rule.
2. *Amendment A (declared 2026-11-09 after the dry run of version 1, before any ISS or baseline catalogue of any held-out patch existed, because the brief asks for at least 4 new clusters)*: the F606W limit of the bright-star rule 14.5 -> **13.0** and the margin
   of the bright-star and saturation rules 10 -> **6 px**; nothing else. I had seen only the counts of surviving origins per cluster for a grid of (limit, margin) values (geometry), and chose the pair for which the first six clusters of the rank order each keep at least
   one patch. The SDSS-saturation rule (pixel maximum) still removes every saturated star; the HST-saturated stars (F606W < 13) are what the limit now excludes.
3. Result (`clusterdep_selection.txt`): six clusters, 8 patches (up to two per cluster, the first two returned by `astro_second.select_patches`, the scarcest density class first):

| cluster | frame (SDSS rerun 301, r) | PSF status / stars | stamp FWHM (px) | patches (x0, y0): F606W < 22 per pixel | descriptive rigid offset of the truth (px) |
|---|---|---|---|---|---|
| M5 (NGC 5904) | run 1458, camcol 4, field 700 | 3 / 0 | 2.51, 2.46 | (640, 240): 0.098; (820, 640): 0.130 | (-0.44, +0.13); (+0.10, +0.21) |
| NGC 6934 | run 1738, camcol 1, field 29 | 2 / 1 | 2.18 | (370, 1230): 0.041; (610, 1220): 0.092 | (-1.00, -0.13); (+0.19, +0.10) |
| M10 (NGC 6254) | run 5391, camcol 5, field 185 | 2 / 5 | 2.81 | (1900, 80): 0.124 | (+1.09, +0.12) |
| NGC 6535 | run 3170, camcol 3, field 49 | 3 / 0 | 3.52 | (150, 1240): 0.040 | (+0.11, +0.23) |
| M71 (NGC 6838) | run 6895, camcol 3, field 52 | 3 / 0 | 2.90 | (1880, 1120): 0.035 | (+0.13, +0.04) |
| NGC 5466 | run 4623, camcol 1, field 324 | 0 / 28 | 3.45 | (1030, 800): 0.033 | (+0.06, +0.35) |

Truth: HUGS F606W used as r, zero point 22.5 and the M2 frame's 0.00546689 nanomaggies per count for every frame, positions by the full asTrans solution, no proper motion, no colour term (as for every earlier cluster). Raw third-party data stay out of git:
`data/astro/fetch_clusterdep.py` and `checksums_clusterdep.json` only.

## 3. Methods, all frozen

* **ISS (frozen).** The configuration of `freeze.json` (sigma_b 3 px, r_ex 1.2 px, kappa 2, the time rule with the units' own-significance clean-up), the guarded path (`imaging_guard.sliding_path_guarded`, installed as in `exp_astro_cyclefix.py`), once per patch:
  (i) *isotropic row*: `exp_astro_m2.run_path`; (ii) *rule row*: the kernel of `psf_general.fit_adequate(stamp)`, `exp_astro_psf.run_path_psf`; where the rule keeps the isotropic fit the two are one run stored under both names. Nothing is re-tuned. **The rule row is the
  ISS of the primary outcomes** (the fresh protocol's primary row). The a priori default and the BIC stop are not run (the fresh test showed they do not decide anything).
* **Repair.** `merge_coherent` at rho\* = 0.8 applied to the kappa = 2 catalogue of each row (the grid of section 1 is also stored, for a sensitivity table that is POST HOC).
* **Comparators (frozen, no re-selection, the catalogues are made by `exp_astro_clusterdep_baselines.py` with the functions of the earlier scripts imported, before the ISS run):** crowdsource with the ISS's front end (`cs2_hp2_t3_pm`, the round-2 selection of
  `baselines_hp2.md`), crowdsource tuned (post-filter `snr5|ff0|qf0|spinf`, `freeze.json`), photutils with the front end (`phu2_hp3_all_i6_b1`), photutils tuned (`photutils_ext_all_local_i6`).
* **Scoring.** `astro_common.score` (match at 0.5 px and 0.5 mag, one to one; truth within 2.5 px of the edge not scored), completeness, purity and F at mag < 22; 95 % percentile intervals of a block bootstrap (20 x 20 px blocks, 2000 resamples, seed 0, stratified
  within patches, the same resamples for every row of a table, so differences are paired). **Primary scoring: the truth shifted by the descriptive rigid offset of the table above (fitted to the pixels before any detection existed)**; secondary scoring: the truth as it is.
  Reason for the choice, declared now: two patches have offsets of 1.0 and 1.1 px, above the 0.5 px match radius, so the truth as it is would score every method near zero there and decide nothing; the offset is a nuisance shared by all methods (it was the
  secondary scoring of the fresh test, whose offsets were at most 0.38 px).

## 4. Outcomes (declared before the run)

Verdict: *win* if the 95 % paired interval of dF (first minus second) is above 0, *loss* if below, *tie* if it contains 0; per cluster, its patches pooled (stratified bootstrap).

* **O2** dF(ISS repaired; crowdsource with the front end); **O2f** the same for the frozen ISS (rule row). **O1/O1f**: against crowdsource tuned. **O3/O3f**: against photutils with the front end. **R** dF(ISS repaired; ISS frozen).
  The isotropic rows (repaired and frozen) are reported beside them, with the same pairs; photutils tuned is reported.
* **Headline tally.** The number of held-out clusters won / tied / lost under O2 (repaired) and under O2f (frozen), with the intervals.
* **Reading of the repair** (decided by the primary scoring only):

| result | reading |
|---|---|
| the pooled R interval (all 8 patches) above 0, and O2 has at most as many losses and at least as many wins as O2f | **the repair reduces the dependence** (stated with n = 6 clusters) |
| the pooled R interval above 0, the tally not weakly better | **partial**: a small real gain that does not change the tally |
| the pooled R interval contains 0 | **no effect**: the dependence is reported as explained (section 1) but not repaired |
| the pooled R interval below 0 | **harmful** |

* **What the paper may say about the lead**, whatever the reading: the tally of the frozen ISS (rule row) over the clusters now scored (5 spent + the held-out ones) with the cluster-level intervals, and the repaired tally over the held-out ones. A single
  loss under O2 at a held-out cluster keeps "cluster dependent" in the statement of limits.
* **Secondary (no claim rests on them):** the same tables with the truth as it is; photutils; the isotropic rows; the sensitivity of R to rho in the stored grid (POST HOC).
* **Dropped:** nothing is dropped after the run. A path that reaches the entry cap (2500) or a guard activation is flagged beside its row and stays. If a baseline fails on a patch the patch stays with the row marked missing (and no claim uses that row).

## 5. Known limits (declared now)

* n = 6 clusters, 8 patches, four clusters with a single patch; three of the six (NGC 5466, NGC 6535, M71) are sparse (0.03 to 0.04 stars per px < 22 mag), where a crowded-field advantage of any method is not expected; the density classes of the selection are the
  scarcest first. The test is of the clusters that the (amended) patch rules let through, not an unbiased sample of globular clusters.
* Three clusters (M5, NGC 6535, M71) have PSF status 3 (no PSF star in the field: the psField stamp is extrapolated). Both methods use the same stamp; a bad stamp is part of what is tested, but it makes the cluster-level verdict a statement about the frame
  as much as about the cluster.
* Amendment A (the bright-star limit 13.0 and the 6 px margins) admits patches with stars of F606W 13 to 14.5 near the box that the earlier clusters' patches do not have; the SDSS-saturation rule is unchanged. Patches with a star of F606W 13.1 to 13.9 inside (M5
  has brightest 13.08, 13.74; NGC 6934 13.59; M10 13.55; NGC 6535 13.96) carry its wings; all methods are scored on them alike.
* The rigid offset is fitted per patch to the HUGS list convolved with the stamp, so the corrected truth is shifted by an estimate with its own error (a few hundredths of a pixel); epochs differ by 9 to 14 years (the proper motion of the cluster is part of the offset,
  the scatter of the stars' motions is not corrected); M2's zero point is used for every frame (0.03 to 0.08 mag); HST-saturated truth stars are missing; blends below the HST resolution.
* rho\* was chosen on spent data (the grid maximum of the mean F over 19 patches); a different value would have been chosen with a different weighting (M92 alone prefers 0.6); the held-out test is of this value.
* The merge is a post-processing of the finished catalogue: it cannot recover a source that the path never entered, and it merges a true close pair of stars that the kernel cannot resolve. The truth also has close pairs: it costs completeness there.
* The residual-calibration prototype (`astro_clusterdep.calibrate`) and the PSF-stretch estimate were attempted and dropped before this plan (section 1); they are not part of the test.

## 6. Run protocol

    python -I experiments/exp_astro_clusterdep_baselines.py            # the baselines' catalogues, once (domain environment; no truth read) -- done before this plan was committed
    python experiments/exp_astro_clusterdep.py --freeze                 # freeze_clusterdep.json: hashes of this file, the code, the registry, the selection, the catalogues
    OMP_NUM_THREADS=1 nice -n 10 python experiments/exp_astro_clusterdep.py --run     # the ISS on the held-out patches, once each; refuses unless the freeze record matches; write-once clusterdep_raw.npz
    python experiments/exp_astro_clusterdep.py --from-raw               # regenerates results/astro/clusterdep.md byte for byte from the raw

One heavy process at a time, `nice -n 10`, `OMP_NUM_THREADS=1`. A bug found after the run is corrected only with the first output kept; the report states both. The raw stores the SHA-256 of the inputs of the diagnosis, the time of each patch's run and
the cap flags and guard records of every path.

## 7. History of this file

* Declared at the time on the first line.
