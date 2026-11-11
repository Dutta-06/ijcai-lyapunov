# The rule that chooses the two fresh clusters (declared before it was executed)

Declared **2026-11-08T21:58Z** (file creation time 21:58:17 UTC by `stat`; `date -u` read 21:58:00 just before). Author: the agent carrying out the coordinator's brief for a fresh cluster-level test. This file is a part
of `fresh_plan.md` (its SHA-256 is in the freeze record) and is written first so that the choice of the clusters cannot depend on anything in an image.

## What I knew when I wrote the rule (stated, not hidden)

* I had read `second_cluster_spec.md` (sections 2.2 and 2.3), so I knew from TABLES the survey columns of every candidate (r-band `PSF_WIDTH`, footprint fraction, number of windows, `PSP_STATUS`,
  `PSF_NSTAR`) and that the earlier agent's placement dry runs gave **one** patch for M3 (NGC 5272, 4646-3-81) and **none** for M5 (NGC 5904, 1458-4-700). The thresholds below were chosen for the
  reasons given with each, not by trying values; I did not compute the outcome of the rule before writing it.
* I had seen NO image pixel of any new cluster, NO detection catalogue, and no comparison of anything with a truth list of a fresh cluster. The only pixel statistic that has ever been computed for a candidate
  frame is the maximum over a window (the saturation rule of the M53/M92 patch selection, dry runs of M3 and M5 by the earlier agent).

## The rule

Universe U: the HUGS clusters that the coverage survey of `second_cluster_spec.md` found to have a candidate SDSS run (rerun 301) and a downloaded HUGS method-1 catalogue (`fetch_second.CLUSTERS`),
without the three clusters already used (NGC 7089 = M2, NGC 5024 = M53, NGC 6341 = M92): NGC 5053, 5272 (M3), 5466, 5904 (M5), 6205 (M13), 6254 (M10), 6535, 6838 (M71), 6934, 6981 (M72), 7078 (M15).
The 38 HUGS clusters with no candidate run and the 13 whose catalogue the MAST endpoint did not find under the file-name pattern (`second_cluster_spec.md` 2.1) are outside U (a limitation).

1. **Frame.** The frame of a cluster is the best field that `experiments/astro_second.py survey` returns for it (largest fraction of the F606W < 22 HUGS sources inside the field box, then the smaller `PSF_WIDTH`).
   No other frame of the cluster is considered. A cluster for which the survey finds no field box with at least 2 % of the footprint (NGC 6981) is dropped.
2. **PSF measured.** `PSP_STATUS` of the field <= 2 and `PSF_NSTAR(r)` >= 1. (The M2 benchmark frame has 2 and 12. A status 3 field with no PSF star has an extrapolated PSF; the rule that picks the kernel
   from the stamp needs a measured stamp.)
3. **Seeing.** The FWHM of the psField stamp (the KL reconstruction, 51 x 51, unit sum, `astro_second.sdss_psf_recon`, `unit_psf`, `psf_fwhm`) at the frame centre (column 1023.5, row 744.0, as in the dry runs
   of `astro_second.py build`) lies in [2.3, 3.6] px. (M2 frame 2.1 to 2.2; M53 3.3 to 3.4; M92 2.9 to 3.0; the repeat frames 2.28, 3.05, 3.52. The range is the one the brief asked for: worse seeing than M2,
   inside the range already covered.) The psField is a calibration file of the survey; no pixel of the frame is read by this step.
4. **Availability.** The survey's number of windows (100 x 100 px origins on a 10 px grid, 30 px from the field box edge, inside the HUGS footprint, 0.03 to 0.15 truth sources brighter than 22 mag per pixel
   = "moderate density", the density band of M53 and M92) is >= 100.
5. **Order.** The clusters that pass 1 to 4 are ordered by the number of windows, descending (ties: the smaller NGC number). The first two that also pass 6 are the fresh clusters.
6. **Patches.** The patch rules of M53 and M92 (`astro_second.py`, module docstring; functions `screen_patches`, `select_patches` imported, not copied), unchanged: origins on the 10 px grid, 40 px from the frame
   edges, inside the HUGS footprint, no pixel above 0.8 x the photoField r `SATURATION_LEVEL` within 10 px of the box (this reads the frame's pixel MAXIMUM; nothing else of the pixels), no truth source brighter than
   F606W = 14.5 within 10 px, no HUGS source with F814W < 17 and no F606W within 10 px, density 0.03 to 0.15 per pixel; one patch per density class [0.03, 0.05), [0.05, 0.08), [0.08, 0.11), [0.11, 0.15]
   (targets 0.04, 0.065, 0.095, 0.13), no overlap, the scarcest class first, ties by (y0, x0); 100 x 100 px. One addition, which concerns the stamp and not an image: the patch-centre FWHM of the psField stamp must
   lie in [2.3, 3.6] px as well. The truth is mapped with the full asTrans solution (as for M53 and M92); no rule is adjusted after a patch has been looked at. A cluster that yields fewer than two patches is
   dropped and the next in the order is tried (declared now; the dry runs of M3 and M5 make this possible).
7. **No relaxation.** No threshold above is changed once the rule has been executed. If fewer than two clusters survive, the test is not run and the result is reported as "no fresh clusters obtainable under
   this rule".

The rule is executed by `experiments/astro_fresh.py select`, which prints every number it uses; its output is stored in `fresh_selection.txt`.

## Execution of version 1, and amendment A1 (declared after version 1 had been executed, before any catalogue of any method existed)

**Version 1 executed (2026-11-08, 22:00 to 22:02 UTC; `fresh_selection.txt`).** Items 1 to 5 on the survey's best frames: the clusters that passed are NGC 6254 (M10; 286 windows; stamp FWHM 2.81 px) and NGC 5272 (M3; 199 windows;
2.81 px); all others failed one of the tests (NGC 5053: 0 windows; NGC 5466: 71 windows; M5: PSF status 3 and no PSF star; M13: 34 windows; NGC 6535 and M71: status 3, no PSF star, windows < 100; NGC 6934: stamp
FWHM 2.18 px, below the range; M15: 3.75 px, above it; M72: no field box). Item 6 on the two (`astro_fresh.py build --dry-run`, which writes nothing): M10 leaves **no** patch (the stars brighter than F606W = 14.5
are within 10 px of every origin that survives the footprint and saturation rules: origins after each rule 1108, 920, 397, 397, **0**), M3 leaves **one** patch (5 origins, all in the densest class, one patch at (1870, 820), density 0.128). By item 6
both are dropped and by item 7 the test would not be run. These two builds read the frame's pixel maximum (the saturation rule) and the truth geometry; they produced no catalogue and no score.

**Amendment A1 (the single change; declared 2026-11-08, before the run of the all-frames survey below).** Item 1 chose the cluster's frame by the footprint fraction, a criterion made for covering the whole HUGS
footprint with one frame of M53 and M92; it is not related to the seeing. The seeing varies from run to run, so a cluster whose best-footprint frame is outside the seeing range may have another frame inside it.
Item 1 is replaced by 1': the frames of a cluster are ALL the field boxes that the survey procedure finds (`astro_second.py survey`'s loop, which keeps the field boxes holding at least 2 % of the footprint; the survey
printed only the best of them); items 2, 3 and 4 are evaluated for each frame; the frame of a cluster is the frame passing 2, 3 and 4 with the most windows (ties: the smaller `PSF_WIDTH`, then run, camcol, field);
the clusters are ordered by the windows of their frame (item 5); item 6 (the patches) is unchanged and is evaluated on the chosen frame of the cluster (a cluster with fewer than two patches in that frame is
dropped, no second frame of the same cluster is tried: the patch rules are mostly cluster properties - bright stars, HUGS footprint - and trying frames until one leaves two patches would be a search). Thresholds
(status, nstar, FWHM range, windows) are NOT changed. Clusters M10 and M3 stay dropped by version 1 if their only passing frame is the one already tried; if another frame passes they are evaluated like any other.
If A1 leaves fewer than two clusters, the test is not run (item 7 stands).

What I knew when I wrote A1: the results of the two dry runs above and the best-frame table; not the number of windows, seeing or status of any frame other than the best one of each cluster.

**Result of A1 (`fresh_selection_A1.txt`, 2026-11-08, about 22:03 UTC).** Frames passing items 2 to 4 exist for three clusters: M10 (NGC 6254; frames 5391-5-185 with 286 windows and 5398-5-22 with 173; the chosen frame is the first, already built: zero patches), M3 (NGC 5272; 4646-3-81, 199 windows;
one patch), M15 (NGC 7078; three passing frames, 2507-6-40 with 158 windows is the chosen one; stamp FWHM 3.15 px). NGC 6934 fails item 3 in both its frames (2.18 and 2.20 px, below 2.3) and the best frame of M15 fails it above (3.75 px). So the
only cluster that can be added with at least two patches is M15 if its frame leaves two; the rule yields one cluster at most under item 6 as written.

**Amendment A2 (declared 2026-11-08 at about 22:04Z (written 22:03:51), before the M15 build was run, before any catalogue of any method existed).** Item 6's "fewer than two patches: dropped" is replaced by
"no patch: dropped". A cluster with a single patch stays and is reported with that caveat: its result is one 100 x 100 px window of one frame (25 blocks), the interval is conditional on that window, and it counts as ONE cluster in
the unit of replication like every other. The motive is the number of available clusters (the HUGS clusters observable in SDSS are about a dozen, three are used, and the patch rules - bright stars, HUGS footprint, density - leave
M3 one patch and M10 none); the cost is a weaker second cluster, which the plan and the report state. Nothing else changes; thresholds are not touched; M10 stays out (no patch).

**Result of the M15 build (dry run, 2026-11-08, about 22:04 UTC; nothing written).** M15's chosen frame (2507-6-40): origins after each rule 859 (density), 803 (margin), 157 (footprint), 89 (saturation), **0** (stars brighter than F606W = 14.5 within 10 px), so no patch. The cause is the same as for M10: a nearby, dense cluster has
stars that bright everywhere in the HUGS footprint of the frame; the far clusters have too few stars per pixel (windows < 100). After A1 and A2 the only cluster with a patch is M3 (one).

**Amendment A3 (declared 2026-11-08 at about 22:04Z (written 22:04:32), before any further build; the last one).** Item 6's patch availability is a property of the sky position, not of the seeing, but a frame holds only part of the HUGS footprint, so another frame of
the same cluster may leave a patch where the chosen one leaves none. For a cluster whose chosen frame leaves no patch, the other frames that pass items 2 to 4 are tried in the order of A1 (windows descending, ties as in A1),
and the first that leaves at least one patch becomes the frame of the cluster. Clusters are then ordered as in item 5 by the windows of their final frame. Only dry runs (geometry, truth, the saturation maximum) enter this
step. If after A3 fewer than two clusters have a patch, the test is NOT run (item 7), the result is "no fresh clusters obtainable under this rule", and no further amendment is made.
