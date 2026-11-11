# Plan: the mechanism on real star-field windows (pre-registration)

Written 2026-11-10; the UTC time of the commit that adds the final version of this file is the registration time (recorded at the end, section 11). No test window (section 3) has been cut, and no
arm has been run on any pixel of a test patch for this study, before that commit. Script: `experiments/exp_mechreal.py` (numerics: `experiments/mechreal_core.py`); freeze record:
`results/astro/mechreal_freeze.json`; report: `results/astro/mechreal.md` (+ `mechreal_raw.npz`, regenerated with `--from-raw`).

## 1. Question

The paper's thesis: small-initialisation gradient training of a factorised Gaussian-kernel network chooses its own number of units under (R) a readout of effective depth >= 3 and (O) an optimizer
that does not normalise each coordinate by its own gradient history. Every test of this MECHANISM so far is synthetic (`results/synthetic/network_imaging.md` trains the network on simulated
frames; the real-data astronomy results use the computed path as a fitter). Question: does the trained network itself, on real SDSS pixels with independent HST truth, (a) count the stars,
(b) lose that ability when (R) or (O) is switched off, as the theory predicts, and (c) form the block of the computed path? This is a test of the mechanism's two conditions. It is NOT a
comparison with any detector: crowdsource, photutils and DAOPHOT are not arms, and no claim of beating a detector will be made from it whatever the outcome.

## 2. Model, loss, units (identical for every arm)

As `exp_network_imaging.py` (tests pin `mechreal_core.py` to that harness: right-hand side, time unit, and a whole trajectory on a simulated frame reproduced to rounding).

* Window: 20 x 20 px cut from a 100 x 100 px SDSS r-band patch. Kernel: the pixel-integrated three-Gaussian fit of the frame's psField stamp (`imaging.fit_psf_parts(stamp, 3, 1.0)`, the
  pipeline's isotropic kernel), high-passed with the border-consistent filter of width 3 px (the window is its own frame). Start lattice: spacing 1 px, margin 1 px, 22 x 22 = 484 units.
* Loss F = 1/2 sum_p W_p (H y - sum_i w_i H phi_i)_p^2, W the pipeline's noise law 1 / (Gaussian-smoothed data (1.5 px) floored at its 1st percentile / gain + dark variance), computed on the
  whole patch and FIXED (no reweighting).
* Units: the data are divided by the flux limit f_lim = 40 / s, s = `imaging.noise_scale` of the window (the median over the lattice of sqrt(sum_p W phi^2)): f_lim is the flux of
  matched-filter significance Z_LIM = 40 (fixed a priori; label-free, from the noise law). The noise scale of the normalised problem is 40 and the time stop of significance kappa is
  tau = 1 / (40 kappa) under the unit wall. A counted star has amplitude >= 1. (The limit is high because the fields are confusion-limited: the development patches hold about 0.1 stars per
  pixel above 5 sigma; above 40 sigma a 15 x 15 px interior holds 2 to 12.)

## 3. Data, development / test split, windows

* Development: the two M2 development patches of `exp_astro_m2.py`, (110, 630) and (210, 630) (tuning data of every earlier astronomy experiment). Truth as it is.
* Test: the held-out patches of `results/astro/clusterdep_patches.json` whose registration did not fail (`clusterdep.md` section 3.3: ngc5466 (1030, 800), ngc5904 (820, 640), ngc6535
  (150, 1240), ngc6838 (1880, 1120), ngc6934 (610, 1220); the three failed fits ngc5904 (640, 240), ngc6254 (1900, 80), ngc6934 (370, 1230) are excluded) and the fresh clusters M3 (ngc5272
  (1870, 820)) and M15 (ngc7078 (640, 1180), (920, 1180)): 8 patches of 7 clusters, none of them M2. Truth shifted by the registry's rigid offset (fitted to the pixels before any detection
  existed). These pixels were used before by the path-as-fitter experiments; no network has been trained on them and nothing of this study has been tuned on them.
* Candidate windows of a patch: the 16 non-overlapping 20 x 20 windows with origins 10, 30, 50, 70 in x and y. A candidate is ELIGIBLE (truth geometry only, no pixel-based outcome): 2 to 12
  counted truth sources in its interior (below), and no HUGS star brighter than 50 f_lim within the window padded by 5 px (unsaturated, moderate dynamic range).
* Selection: per patch a seeded permutation of the candidates (`numpy.random.default_rng([20261110, 0 for development / 1 for test, patch index])`); development: the first 4 eligible of
  each patch (8 windows); test: the first 5 eligible of each patch (up to 40 windows; a patch with fewer eligible candidates contributes what it has). If fewer than 30 test windows result,
  the test is run on what there is and the number is reported; nothing is added. The other eligible development candidates ('devx') are used for arm F only (tuning of the stop).
* Truth of a window: HUGS stars (F606W as r, `astro_common.mag_to_counts`), flux in units of f_lim; stars closer than 0.5 px (single linkage, stars above 0.1 f_lim) are ONE source (summed
  flux, flux-weighted centre); a source is COUNTED when its flux is >= f_lim and its centre lies in the interior (2.5 px from every side, `astro_common.EDGE`). K_true = their number.

## 4. Arms (same dictionary, same loss, same budget rule)

| arm | readout | optimizer | start | clock (t = clock x tau) |
|---|---|---|---|---|
| A (the mechanism) | paper's unit w = nu^2 a_max tanh(alpha / a_max) / m (m 12, G 1, gamma 0.1, a_max 1) | gradient flow (Radau, rtol 1e-8; BDF fallback) | alpha 0, nu^2 = eps | Lambda(eps) |
| B | depth 2, w = u^2 - v^2 (rates 1, centres G m) | gradient flow | u = v = sqrt(eps) | ln(1 / eps) / 4 |
| C | depth 1, w = u (rate G m) | gradient flow | 0 | Lambda(eps): the physical time of arm A (depth 1 has no wall, hence no clock of its own; the convention of `exp_parametrization`) |
| D | paper's unit | Adam, default settings (step 1e-3, betas 0.9 / 0.999, eps_adam 1e-8), raw gradient, one step size for all parameters and centres | as A | none (100000 steps) |
| E | paper's unit | Adam with ONE shared second moment (mean squared gradient after the flow's group rates), step 1e-3 | as A | none (100000 steps) |
| F (reference) | - | the computed path `imaging_guard.sliding_path_guarded` (dual order, global refits, no exclusion zone, fixed widths and weights) | - | tau |

eps = 1e-8 for A to E (the value at which the simulated-frame harness had 18 to 21 of 24). Budget (the same rule for every arm and for development and test): 3600 s of wall clock per run on the machine that runs it, one core; a run that has not reached a snapshot by then
is recorded as NOT CONVERGED: it is scored at its last accepted state (flagged) and counts as a failure for P-path. Heavy runs go to private CPU workers; the same code
(`mechreal_core.work`) runs there. 

## 5. Stop rule (label-free, declared; chosen on development windows only)

* Primary: the time stop of the path, tau = 1 / (kappa x 40), in the clock of each flow (A, B, C) and for F. kappa and the formed threshold theta are chosen from arm F ALONE on the
  eligible development candidates (grid kappa in {40, 30, 20, 15, 10}, theta in {0.25, 0.5, 0.7, 0.85, 1.0}): among the configurations whose mean count error is within one standard error
  of the smallest, the largest kappa (the earliest stop), then the smallest error, then the largest mean F, then the largest theta. The network arms do not enter this choice. The network
  arms are integrated on development windows to kappa = 20 at the latest, so the choice is restricted to kappa >= 20 if the rule returns a smaller value (this would be reported).
* D and E have no clock. Each is trained for a fixed budget of 100000 steps and scored THERE (primary; no stop is tuned for them). Why this and not a tuned step count: the first draft
  of this plan chose for each of them the step count of smallest development count error; development showed that this stops Adam before it has formed anything (5000 to 10000 steps:
  the loss is already below the path's, carried by hundreds of small amplitudes, none above theta), which is not a trained network. On every development window Adam's state no longer
  changes after 50000 steps; E is still slowly forming units at 200000 (it runs on a slow clock), so its reading at a fixed budget is to that extent arbitrary; the equal-fit reading
  below is the clock-free one. The grid of step snapshots is 1000, 2000, 5000, 10000, 20000, 50000, 100000 (development also has 200000; 200000 steps take about 2600 s on the worker).
* Secondary (registered, deciding nothing): (i) "best case": B, C, D, E each at its own best development snapshot AND threshold, i.e. tuned on the truth of the development windows,
  which A is not (for C the snapshots include early times, 1e-6 to 1e-1 of tau(kappa = 20), since depth 1 has no clock; for D and E the step grid up to the budget); (ii) "equal fit"
  for D and E: the first snapshot whose training loss is at most the loss of arm F's block at the time stop (label-free; the last one if none is); (iii) the BIC stop, the same
  label-free rule for A to E: per window the reached snapshot of smallest 2F + 3 K ln N. Also reported for every arm: the number of units with |w| > 0.05 anywhere on the lattice
  ("awake" units; dormancy is the statement that this number stays small).

## 6. Outcomes per window (all arms alike)

Formed units: |w| > theta (either sign); same-sign formed units within R_MERGE = 0.5 px (single linkage) are one unit (summed amplitude, |w|-weighted centre); K_hat = those in the interior.
Count error = |K_hat - K_true| (primary metric); signed error; detection F = 2 TP / (K_hat + K_true) at a match radius of 0.5 px (one-to-one, maximum cardinality then minimum distance).
For arm A: agreement of its formed set with the path's block, taken exactly as the harness does: at the snapshot tau_mid between the last entry before the stop and the next entry
(capped at 4 tau_K), against arm F's block polished to its stationary point, with the harness's formed threshold min(0.01 max|w_block|, 0.1 min|w_block|); the formed SET agrees iff C1
(same count, every formed unit within 0.5 px of a block unit) and C2 (one-to-one within 0.25 px) hold.

## 7. Predictions and decision rules (each reported as HELD or FAILED with its interval)

Intervals: 95 % percentile bootstrap over the test windows (20000 resamples, seed 20261110), paired differences. All at the frozen primary stop, eps = 1e-8, theta, R_MERGE 0.5.

* **P-R (readout depth)** holds iff all of: R1 the mean of (count error of B - count error of A) has its interval above 0; R2 the same for C; R3 C forms many more units than A: the
  interval of mean (K_hat of C - K_hat of A) is above 0 and mean K_hat of C >= 2 x mean K_hat of A. Registered companion R3b (reported, not part of the verdict): no dormancy, the
  same two conditions on the units with |w| > 0.05 anywhere on the lattice. (Recorded before the test: R3 is posed on the thresholded count K_hat, as the task words it; the factor 2 was
  fixed before any run of C. On the development windows the ratio of K_hat is 1.7 and R1 does not hold (B and A do not differ); the thresholds are not changed.)
* **P-O (optimizer)** holds iff both: O1 Adam (D) forms more units and has a larger count error than A: the intervals of mean (K_hat D - A) and of mean (count error D - A) are both above
  0; O2 E is within the margin of A: the interval of mean (count error E - A) lies inside [-1.0, +1.0]. Registered companion O2b (reported, not part of the verdict): E forms no excess
  units, the interval of mean (K_hat E - A) lies below +1.0. (Note recorded before the test: E runs on a slow clock, so whether it under- or over-forms at the fixed budget depends on the window;
  O2 is decided at the budget as registered, and the equal-fit reading is reported next to it.)
* **P-path** holds iff arm A (eps = 1e-8) forms the path's block (section 6) on at least 60 % of the test windows; a run that is missing or did not reach tau_mid is a failure. The exact
  binomial interval is reported. (The harness had 18 to 21 of 24 simulated frames.) Two readings are registered, each HELD or FAILED on its own, same 60 %:
  **P-path** proper, against the block of arm F (the computed path as the pipeline uses it, non-negative amplitudes), and **P-path-S**, against the SIGNED path (Definition 1 of the paper,
  `exp_network_imaging.signed_path`; the network's readout is signed and F is the restriction to non-negative amplitudes). Recorded before the test, from development (section 10): on real
  high-passed pixels the network forms NEGATIVE units, which F cannot have; at eps = 1e-6 arm A formed the block of F on 0 of 8 development windows and the signed block on 2 of 8. The
  threshold of 60 % was fixed from the simulated-frame harness and is not lowered: both readings are expected to be at risk and will be reported as they fall.
* **Informative?** If the mean count error of the reference path F on the test windows exceeds 0.5 x mean K_true, the real PSF / confusion makes every arm fail: this is then reported
  plainly as a negative and the predictions are reported as uninformative.
* What counts against the thesis and will be stated: any of P-R, P-O, P-path FAILED; a best-case or BIC reading in which B, C or D is as good as A; a count error of A clearly above F's.
* Dropped: nothing after the fact. Anything computed after the test results are seen is labelled POST HOC in the report.

## 8. Secondary analyses (registered; none decides a verdict)

Best case, equal fit and BIC stop (section 5); merge radii 0.25, 1.0, 1.5 px; positive units only; windows where every arm reached its snapshot; cluster bootstrap over patches; by cluster; arm A
against the SIGNED path (Definition 1; the network's readout is signed, arm F is the non-negative restriction) and with the harness's full criterion C1 to C4.

## 9. Process

Development (M2 windows) -> `--freeze` (kappa, theta, step counts, best cases, hashes of the code, the plan and the registries) -> commit -> `--select test` (refuses without the freeze
or if a hashed file changed) -> arm F and the network arms on the test windows, ONCE (the raw is write-once for test keys) -> report. A bug found after the test started is fixed, the fix
and its effect are reported, and affected runs are rerun in full (not selectively).

## 10. Development record (everything here is from the two M2 development patches; no test patch is involved)

Tables: section 3 of the report. Values chosen by the rules above (`mechreal_freeze.json`): kappa = 30, theta = 0.85 (arm F on the 24 eligible development candidates: smallest mean count
error 0.88 at kappa 10, theta 1; one standard error 0.18; within it the earliest stop is kappa 30, error 0.96, where theta 0.85 has the larger F). At that stop, on the 8 development windows
(mean K_true 6.5), mean count error and mean number of units with |w| > 0.05 (of 484): A 1.38 / 17; B 1.50 / 39; C 5.62 / 424; D (100000 steps) 21.1 / 441; E (100000 steps) 1.75 / 28;
F 1.12 / 13. Every development run finished (A at most 1078 s, B 139 s, C 374 s; D and E about 1300 to 2700 s for 200000 steps). Arm A formed the block of F on 0 of 8 windows (eps 1e-6 and
1e-8) and the signed block on 2 of 8 (1e-6) and 1 of 8 (1e-8). Best cases (own snapshot and threshold, tuned on the development truth): A 1.00, B 1.00, C 2.00 (at 1e-5 of tau(20), theta
0.5), D 4.25 (1000 steps, theta 0.25), E 1.50.

Changes made during development, before this plan was committed (none after): (1) the tolerance scale of the depth-1 amplitudes (1e-10 absolute instead of 1e-10 sqrt(eps): they start at 0);
(2) the stop of D and E (section 5: a fixed training budget instead of a step count tuned on the development truth); (3) the step grid cut at 200000 (the budget); (4) P-path given its
second reading (signed path) after the negative units were seen; (5) the companion R3b and the "awake units" column; (6) a faster evaluation of the kernel terms (`fast_comps`, pinned to
`imaging.PixelPSF.comps` by a test; stored development runs are reproduced by the final code: a depth-2 run exactly, a shared-moment run to 1e-4, the size of machine-to-machine
rounding amplified by the normalised steps). The development runs of D and E were made on the grid up to 200000 steps; the test runs stop at the budget of 100000.

## 11. Registration

The UTC time of the commit that adds this file together with the freeze record is the registration time. Written at 2026-11-10T20:03Z; at that moment no test window has been selected
(`--select test` refuses to run without the freeze record), and no arm, the reference path included, has been run on any test patch for this study.
