# Plan (pre-registration): does the sliding Gaussian path count pulses in real 1-D signals better than tuned classical counters?

Written 2026-11-10 09:26 UTC (the commit time is the record), before any TEST record was scored. What had been looked at by then, in full: (a) the dataset
listing and the archive's aggregate label statistics over all 2000 records (beat counts, RR intervals; no method was scored on them); (b) a handful of individual
records (ids 1, 2, 3, 500, 1200, 1500) run through prototype fits to measure compute time; (c) the DEVELOPMENT records (150, defined below) through the front-end
choice and the path families (stages `fe`, `dev`, `freeze`; their results are in `results/ml/peakcount_dev.md` and `peakcount_freeze.json`, committed with this plan).
The 500 test records have been loaded by none of the stages. Code: `experiments/exp_peakcount.py`; tests: `tests/test_exp_peakcount.py`. This is a SPIKE.

## Question
The theory (proved in 1-D) says that gradient flow of a Gaussian-unit network with a factorised readout, started small, switches units on one at a time, slides them onto
features and leaves the rest dormant, so the NUMBER of units is chosen by training. The narrow test asked for: on REAL 1-D signals whose number of components has an
independent ground truth, is the number found by the path better than the number found by well-tuned classical counters?

## Data and ground truth
CPSC 2019 challenge training set (Zenodo 23179391, 2000 single-lead ECG records, 10 s at 500 Hz; SHA-256 recorded in `exp_peakcount.py`). Ground truth: the expert-annotated R-peak
positions of each record; target = their number (6 to 38, mean 14.7). Why this dataset and its limits (QRS is not a Gaussian, the beats rarely overlap, label count equals model
count only after a front end): `results/ml/peakcount_scout.md`. The raw archive is not in git; `--fetch` restores it.

## Split (fixed by a seeded rule)
Record ids 1..2000 are permuted with `numpy.random.default_rng(20261110).permutation`; the first 150 are DEVELOPMENT, the next 500 are TEST, the remaining 1350 are not used (reserved).
All choices are made on the development records; the test records are run once.

## Front end (shared by every method)
Band-pass 5-18 Hz (2nd-order Butterworth, zero phase), derivative, square (`sq`) or absolute value (`abs`), 120 ms moving mean, decimation to 50 Hz (500 samples), division by the 95th
percentile. The front end is chosen on development by the lowest mean absolute count error of the tuned `find_peaks` counter (the classical counter picks it, which favours the baselines);
the choice is frozen in `peakcount_fe.json`.

## Methods and equal tuning budgets (each family: 42 development configurations, the one with the lowest development mean absolute count error is frozen; ties go to the first)
* `fp`: `scipy.signal.find_peaks` on the envelope, prominence in {0.03 ... 0.35} (7 values) x minimum distance in {0.15 ... 0.40} s (6).
* `fp_rel`: the same with a prominence relative to the record (q x the 75th percentile of the prominences of a low-threshold pass), q in {0.10 ... 0.60} (7) x distance (6).
* `mf`: a Gaussian matched filter (width in {0.03, 0.05, 0.07} s) then `find_peaks`, prominence (7) x distance in {0.2, 0.3} s.
* `step`: forward selection of grid Gaussians by least squares (the explicit-method stand-in for lmfit forward selection with BIC; `lmfit` is not installed), `baselines.stepwise_bic`.
* `omp`: matching pursuit with sliding (`baselines.omp_slide`), the greedy cousin of the path; this is the "not worse than greedy" comparator.
* `iss`: the sliding inverse scale space path (`baselines.iss_path`, the limit form of the factorised small-start network; not a trained torch network), OURS.
The three path families (`step`, `omp`, `iss`) use Gaussians of one width in {0.03, 0.04, 0.05} s on a start grid of spacing 0.1 s, at most `min(44, int(1.4 P) + 4)` entries where P is the
number of peaks of a loose label-free `find_peaks` pass, and the same two stopping rules: the penalised BIC `N ln(rss/N) + c * 2 k ln N`, c in {1, 2, 4, 8, 16, 32, 64}, or the residual-drop rule
(stop at the first entry that lowers the residual by less than d times the empty-fit residual), d in {0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08}: 3 x 14 = 42 configurations each. The count
is the number of units of the selected block. (Kernel widths above 0.05 s were not run because the sign-constrained polish of the `iss` path then takes 20 to 100 times longer: a compute limit
chosen before any accuracy was seen.) Deep-learning counters (published CPSC 2019 entries) are not run; they would be the ceiling and are outside a CPU spike.

## Metric and verdict (single primary metric: mean absolute count error per test record)
Let `A_m(i) = |K_hat_m(i) - K(i)|` for record i. `best` is the classical family among {`fp`, `fp_rel`, `mf`, `step`} with the lowest DEVELOPMENT mean absolute error (frozen). The paired
bootstrap (10000 resamples of the test records, seed 20261111, percentile 95% interval) is applied to the differences of the per-record absolute errors.
* **Win condition (W)**: the interval of `mean(A_best - A_iss)` lies above 0.
* **Non-inferiority to greedy (N)**: the upper end of the interval of `mean(A_iss - A_omp)` is below 0.10 beats per record.
* **GO if W and N; KILL otherwise.** (Equal error, a loss, or a win over `best` that is not clear of the greedy control are all KILL.) Computed by `verdict()` and printed in `results/ml/peakcount.md`.
* Secondary, reported but not part of the verdict: exact-count share, signed error, position F1 (one-to-one matching, tolerance 0.10 s, after removing the shared envelope-to-annotation
  offset estimated on development), paired differences against every family, errors by true count. Everything not decided here (for example `ours` against `fp` when `fp` is not `best`) is
  descriptive.

## What a result would and would not mean
A GO would mean: on real pulse trains with an independent count, the path's count is better than that of tuned peak counters on a shared front end, and not worse than greedy sliding pursuit.
It would NOT show that the dual order of entry matters (that is the `omp` comparison) and not that Gaussian-unit networks beat modern deep QRS detectors. A KILL, the expected outcome given that
beats rarely overlap (the regime where resolving close components matters), means that on this data model-order selection by the path buys nothing over a tuned peak finder. Errors mix the
model-order error with the front-end error, identically for every method.

## What will be dropped
Nothing after scoring: all 500 test records stay in the raw, a failed fit counts as its count error. Anything done after seeing test results is labelled POST HOC in the report.

## Amendments
1. Before the plan commit (development work only, recorded here for completeness). (a) Kernel widths above 0.05 s were dropped from the path families for compute (see Methods). (b) A third stopping rule for the
   path families, a dormant-unit threshold on the weights of the last block (the analogue of the repository's `threshold_units`), was prototyped on 25 to 40 development records and dropped: units that sit off the
   annotated beats have weights as large as the units on beats (medians 0.56 against 0.76 of the 95th-percentile-normalised envelope), and the rule's error was 11 to 18 beats per record. It is not part of the spike. (c) The
   first full development pass (`peakcount_dev.md`) gave, as mean absolute count error on the 150 development records: matched filter 1.45 (frozen `best`), find_peaks 1.55, relative find_peaks 1.65, forward selection
   3.47, matching pursuit with sliding 3.54, sliding inverse scale space path 3.64. So the expected outcome of the test is KILL: the path families are about two beats per record worse than the peak counters on
   development. The test is run anyway, once, to obtain the held-out verdict; nothing in the rule, the grids or the frozen configurations is changed after this note.
2. Correction of the header time, made before the test run: the header says 09:26 UTC, which was an estimate typed by mistake; the commit of this plan (b532f64) was made at 2026-11-10 09:23 UTC (see `git log`),
   and the test run starts after the commit of this amendment. No test record had been loaded at either time.
