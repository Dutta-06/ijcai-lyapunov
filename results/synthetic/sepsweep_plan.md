# Pre-registration: a separation sweep on the TRAINED network (workstream X1)

Script: `experiments/exp_sepsweep.py` (engine, design, outcomes, statistics, report); tests: `tests/test_exp_sepsweep.py`; report: `results/synthetic/sepsweep.md`; raws: `sepsweep_dev_raw.npz` (development) and `sepsweep_raw.npz` (held-out).  This plan is committed before any held-out target is run; the commit time is recorded below.  Plan committed: **2026-11-10 23:36 UTC** (the git commit time of this file is the authority).  Everything not listed here is POST HOC and is labelled so in the report.

## 1. The question

Two independent simulated reviewers asked for the same experiment: "a separation sweep on the TRAINED network (not the path): K = 2-6 Gaussian features at separations 1-8 sigma and amplitude ratios 1-1e3, noise from 1% to 30% of the amplitude, many targets, eps from 1e-2 to 1e-8: show where the trained network deviates from the path and compare the threshold 4.56 sigma with the empirical transition."  The paper's separation condition S has computed conservative thresholds 4.56-6.70 sigma (Gaussian, K <= 3, ratios 1-10, `theory3e.md` Table A) for the limit path; the path itself is exact on 1154 targets down to 1.5 sigma; the earlier trained-network tests (T1, `theory_t1.md`) were above 4.9 sigma and did not discriminate between "the threshold is the truth" and "the threshold is conservative".  The paper says S is "sufficient, not a transition"; this experiment measures the transition of the network itself and where it parts from the path.

A reduction of the reviewers' grid, fixed here: K in {2, 3, 4} (not 6), Delta up to 7 (not 8), R up to 100 (not 1000), noise up to 20 % (not 30 %), eps in {1e-4, 1e-6, 1e-8} (not 1e-2 and not 1e-8 only): the budget of one core for the whole of it is about 20 core-hours at the measured cost of 9.5 s per target (path plus three flows).

## 2. What was seen before this plan (development data only)

Timing and smoke runs of the engine on throw-away targets (seeds 999 and 5, about 60 targets in all, none of them the development or test seeds): the path and the network both succeed on every K = 2, Delta = 7 noiseless target tried; at eps = 1e-4 the network represented the second feature of a noiseless Delta = 7 pair by two co-located units (merged success, formed set not equal to the path's); one K = 4, R = 10, noise 0.2 target took a path of 4e5 right-hand-side evaluations.  No success rate of any cell was estimated.  The numbers 4.56-8.54 of section 4 come from the repository's own threshold code; the five of them in `theory3e.md` were reproduced.

## 3. Machinery (all existing; nothing is retuned)

* The network: `factorized_rbf.PaperLaw` with the law (m, G, gamma, a_max) = (12, 1, 0.1, 1), the centres' rate Gm, trained by the exact gradient flow `factorized_rbf.train_flow` from the start nu_i^2(0) = eps, alpha_i(0) = 0, integrated by Radau (rtol 1e-9; the LSODA default of the function wastes 1.5e5 evaluations when its stiffness detection fails, and Radau gave the same final coefficients to about 1e-9 on the three targets compared), 601 recorded snapshots on tau in [0, 1.5 tau_K + 1] where tau_K is the path's K-th entry time on the same data.  The kernel matrix is cached between the two calls of one right-hand-side evaluation (`CachedBumps`; bit-identical, 1.8 times faster).
* The dictionary and the data: kernels of width sigma = 0.12 on the lattice of spacing 0.1 (0.833 widths, the benchmark's) over a window of `[0, W/10]`; samples at spacing 1/30; the window is the span of the K features plus 0.5 (4.2 widths) on each side, rounded out to the lattice (translation invariance: the window's place is irrelevant; the random position is the lattice phase).
* The path: `rich_transit_nd.sliding_iss_nd` (polish on, the repository's `reference` setup) on the SAME samples, run to K + 2 entries with a budget of 6e5 right-hand-side evaluations (a path that runs out of budget is a recorded failure of the path and of the targets' agreement; its tau_K is then the last entry time and the flag is recorded).
* Targets: K features with neighbours exactly Delta widths apart (so all pairs are at least Delta), the first at 0.5 + 0.1 u (u uniform: a random lattice phase), amplitudes A_min R^(j/(K-1)), j = 0..K-1, in a random order, all positive, A_min = 0.7 (the scale of theory3b / T1); white noise of sd = (noise) x A_min added to the samples (noise relative to the SMALLEST amplitude, the unit of the theory's gauges); the same unit-variance noise vector for the noise levels of one target (paired).

## 4. Design (fixed)

| factor | levels |
|---|---|
| K | 2, 3, 4 |
| Delta (widths) | 1.5, 2.5, 3.5, 4.5, 5.5, 7 |
| R (largest / smallest amplitude) | 1, 2.4, 10, 100 |
| noise sd / A_min | 0, 0.01, 0.05, 0.2 |
| eps | 1e-4, 1e-6, 1e-8 (every target is trained at all three) |

Drop rules, fixed beforehand: (D1) K = 4 only at R in {1, 10}; (D2) the intermediate noise levels 0.01 and 0.05 only at R in {1, 10} (R = 2.4 and 100: noise 0 and 0.2).  This gives **192 cells** (K = 2: 72, K = 3: 72, K = 4: 48).  Held-out size: **40 targets per cell = 7680 targets** (each: one path and three flows); development: **6 per cell = 1152 targets**, seeds from 7 100 000 (development) and 8 200 000 (test), separate; the position, amplitudes and noise of a target come from `SeedSequence([seed, K, iDelta, iR, idx])`.  If part of the held-out run cannot be completed (the compute is a private CPU notebook, see section 9), cells are dropped WHOLE, in the order of the shards, and the report lists the cells that are missing; no claim is made for them and a (K, R) with a missing noise-0 or noise-0.2 cell is left out of the predictions.

Computed thresholds Delta_entry(K, R) (Proposition G, blockwise spectral norm, grid 0.02, d_0 = 0.4167, Gaussian Gram; `theory3e_kernels.entry_threshold`): (2, 1) 4.56, (2, 2.4) 5.44, (2, 10) 6.70, (3, 1) 5.00, (3, 2.4) 5.80 (these five are the ones printed in `theory3e.md`); (2, 100) 8.22, (3, 10) 6.96, (3, 100) 8.42, (4, 1) 5.22, (4, 2.4) 5.98, (4, 10) 7.10, (4, 100) 8.54 (the same function; not printed in the paper).  The tests re-derive the values.  The thresholds for R = 100 and (4, 10) lie beyond the largest separation tested (7).

## 5. Outcomes per target (all computed from the same snapshot of the network, and of the path, under one configuration)

A *configuration* is (stop rule, formed threshold, merge radius).  A snapshot has coefficients w_i and centres c_i; unit i is *formed* when |w_i| exceeds the threshold; the formed units are merged by single linkage on the line (neighbours within the merge radius join; the merged centre is the |w|-weighted mean, the merged coefficient the sum).

* **(i) success of the network**: after merging there are exactly K units, one per feature (matched by order), each within 0.5 widths of its feature and within 50 % of its amplitude.  The paper's strict exactness (0.05 widths, 5 %) is recorded as `exact` and read for the noiseless cells only.
* **(ii) agreement**: the set of formed dictionary units of the network (unmerged, same threshold) equals the formed set of the path's block at its stopping state.
* **(iii) success of the path**: the same success criterion applied to the path's state (the state after its K-th entry for the stop rules A, the state of least BIC for BIC).
* Stop rules: **A05 / A25 / A50** = the network's snapshot nearest tau = (1 + omega) tau_K with omega = 0.05 / 0.25 / 0.5 (the paper's "the path's time of the K-th entry plus a window"; K is used, as in the paper); the path's state is its state after the K-th entry whatever omega; **BIC** = the snapshot (the path: the state) minimising N ln(rss/N + 1e-10 mean(y^2)) + 2 k ln N, k the number of merged units under the configuration (no use of K or of the truth).
* Formed threshold: **P** = min(0.01 max|A*|, 0.1 min|A*|) (the paper's: `factorized_rbf.formed_threshold`), **D** = 1e-3 max|w| of the snapshot (no use of the truth).
* Merge radius in widths: **0, 0.25, 0.5, 0.75**.
* The 32 configurations = {A05, A25, A50, BIC} x {P, D} x {0, 0.25, 0.5, 0.75}.  Every configuration is scored on every target and stored (the raw keeps the outcome bits, not the snapshots).  Flows that time out (300 s) or fail are failures in the rates and are counted in the report.

## 6. Development rule (development targets only; the configuration is frozen before the test)

The development targets are run once; the **objective** of a configuration is the mean over the development targets at Delta in {5.5, 7} (all K, R, noise) of [network success at eps = 1e-8 + path success]/2 -- the cells where both are expected to succeed, so the transition region is not tuned.  The configuration of maximal objective is chosen; among those within 0.005 of the maximum the first in the listed order (A05, A25, A50, BIC; P, D; radius ascending).  It is written into `exp_sepsweep.FROZEN` and recorded in the freeze record below, committed before any test target is run; `tests/test_exp_sepsweep.py` checks that the frozen configuration is the rule's choice on the stored development raw.  The other 31 configurations are scored on the test and reported as a registered **sensitivity** (does the transition move with the stop rule / threshold / radius?), not as the verdict.

## 7. The single primary metric and the verdict rule

For each (K, R) at noise 0, eps = 1e-8: the empirical **transition** of the network, Delta*_net = the smallest grid separation from which the success rate exceeds 0.9 at every larger grid separation, linearly interpolated from the previous grid point (`<=1.5` if the rate already exceeds 0.9 at 1.5, `>7` if it never does at 7); interval = the 2.5th-97.5th percentile of the transition of 2000 parametric bootstrap resamples of the per-cell success counts (seed 0); the rates carry exact (Clopper-Pearson, 95 %) intervals.  The same quantity for the path (Delta*_path, same targets) and at every eps and noise level.  The comparator is the computed threshold of section 4 for the same (K, R).

* **Verdict "the computed threshold is conservative for the trained network"** iff P1 holds (below); the paper may then say "the trained network forms exactly the features at all separations tested at or above ... widths, below the computed thresholds", with the measured slack.  **Verdict "the threshold is not conservative at ..."** if P1 fails for some (K, R): the paper must say where.
* Whatever happens, the paper states the measured transitions with their intervals, the eps and noise dependence, and the discordance between network and path (section 3 of the report: success rates, agreement, and the per-target tallies of "network succeeds, path does not" and the reverse).

## 8. Registered predictions (stated before any held-out run; held / failed is reported plainly)

* **P1**: at eps = 1e-8 and noise 0 the network's transition is at or below the computed threshold for every (K, R) [point estimates; the number of (K, R) whose upper interval end is also at or below the threshold is reported].  `>7` counts as above any threshold, `<=1.5` as below.  (For R = 100 and (4, 10) the thresholds are above 7, so P1 only asks Delta* <= 7 there.)
* **P2**: at eps = 1e-8, noise 0, the network's transition is within +1 width of the path's, for every (K, R) (network minus path <= 1; both `>7` counts as a difference of 0).
* **P3**: noise at 0.2 of A_min moves the network's transition (eps = 1e-8) above the computed threshold, for every (K, R) whose threshold is at most 7 (7 of the 10 (K, R) of the design).
* **P4**: the larger start eps = 1e-4 moves the transition up relative to eps = 1e-8 at noise 0: not lower in any (K, R) and higher by more than 0.25 widths in at least two thirds of them (7 of 10).
* The success is also expected to be monotone in Delta within noise; non-monotone cells are reported, not smoothed.

## 9. Execution and what is dropped

Compute: the development and held-out runs are large (the held-out run is about 20 core-hours at the machine's measured cost); they run as private notebooks on a CPU service, four processes each, a neutral title, no data and no text posted anywhere public.  A cross-check: the first 12 development targets are also run on the local machine and their stored outcome bits are compared (reported in the freeze record); the platforms' numerical libraries differ, so single targets near a tolerance may differ; the comparison is reported, not assumed.  Dropped, if needed (in this order): cells of the shards that did not finish.  Not run: K = 5, 6; Delta = 8; R = 1000; noise 0.3; eps = 1e-2 (the reviewers' grid beyond section 4); the plane; non-Gaussian Grams; unequal signs.

## 10. What counts as what

* **Held** / **failed** for P1-P4 as defined; there is no "win" or "tie" (the experiment compares the network with the threshold and with the path, not two methods).
* The agreement curves and the discordance are descriptive; they are not predictions.
* The 31 non-primary configurations and the BIC stop are sensitivity, registered here.  Anything else (a different level, an interval from another rule, a subset of cells chosen after seeing) is POST HOC.

## 11. Amendment 1 (before the development run finished, before any held-out target)

* The held-out run is split into three shards of cells (cell index mod 3; one on the local machine, two on private CPU notebooks, which run two at a time and whose cores are about four times slower than the local one), and inside a shard the cells run in **priority tiers**: tier 0 = noise 0 and 0.2 (all the predictions), tier 1 = noise 0.05, tier 2 = noise 0.01; within a tier in the cell order.  If the time runs out, the cells not reached are dropped whole (lowest tier first) and listed in the report.  A (K, R, noise) group enters the tables only if all six of its separations have the full 40 targets (as many as the fullest cell); an incomplete group is listed, not analysed.
* The report carries a section 8 (sensitivity: the verdicts P1-P4 and the median transitions under each of the 32 configurations on the same test targets), registered in section 6; the transition's bootstrap interval is computed with the vectorised crossing (identical values; tested against the row-by-row function).
* The development run (1152 targets, one private notebook) was started with the engine as committed in the plan commit; the later edits of the script touch the report, the shard order and `group_keys`, not the targets, the path, the network or the outcome bits.

## 12. Amendment 2 (development objective; before the development run finished, before any held-out target)

The objective registered in section 6 reads only Delta in {5.5, 7}.  A local development run of the 313 first targets of cells 96-191 (K = 3 and 4, tiers 0, development seeds; a check of the second machine, `ssout3`) showed, noiselessly, network success at eps = 1e-8 of 0.88 / 0.92 / 0.93 under the windows 0.05 / 0.25 / 0.5 (all other choices equal): the network lags the path on the smaller feature of unequal pairs at Delta <= 3.5, so the shortest window is cut off before it has formed.  On the well separated cells all windows tie, and the tie-break of section 6 (the first in order) would then take the shortest.  The objective is therefore replaced by **the mean over ALL development targets (all separations) of [network success at eps = 1e-8 + path success]/2**, the rest of section 6 unchanged (tolerance 0.005, first in order).  This tunes the configuration to the network's own success on development data, which favours the network; the registered sensitivity (section 8 of the report: the 32 configurations on the test targets, including the strict radius 0 and the shortest window) is what shows how much each verdict depends on it.  `DEV_DELTAS` is removed from the script, the tests and the mutants follow.

## 13. Freeze record



(appended before the held-out run)

Written 2026-11-11 01:53 UTC, after the development run (1152 targets, one private CPU notebook, 2 h 10 min; outputs `ckpt_dev_s0.jsonl`) and before any held-out target was run.

* Development objective (amendment 2) over all 1152 development targets: maximum 0.9492 (A50, P, radius 0.75); the configurations within 0.005 of it include A25 / P / radius 0 (the first of them in the listed order, since all A05 configurations are lower).  **Frozen configuration: stop rule A25 (the snapshot at 1.25 tau_K), formed threshold P (min(0.01 max|A*|, 0.1 min|A*|)), merge radius 0 widths (no merging).**  This is the strictest merge radius: a feature carried by two co-located units counts as a failure of the network.  The table of the 12 best configurations is in the report's section 0.
* Platform cross-check: 839 development targets were also run on the local machine (Python 3.13 with numpy 2.5 / scipy 1.18 against the notebook's numpy 2.1 / scipy 1.16): 32 of 839 (3.8 %) differ in the success bit of at least one of the 32 configurations, 3 of 839 (0.4 %) in the frozen configuration; the largest relative difference of tau_K is 4.5e-8.  The platforms agree except on borderline targets.
* The development raw is the notebook's (1152 targets); the local runs are not used in it.  The partial local first half of the development run (cells 0-95) was stopped once the notebook finished.
* Run log: the container was restarted during the development run (the local part was resumed from its checkpoint; the notebook was unaffected).  No held-out job has been launched before this freeze.

## 14. Run log of the held-out run

* Launched 2026-11-11 about 01:55 UTC with the frozen configuration: shard 0 of 3 on the local machine (resumable checkpoint), shards 1 and 2 as private CPU notebooks (`sweep-run-b`, `sweep-run-c`; 4 processes each), each shard the cells with index mod 3 equal to its number, tiers 0, 1, 2 in that order.  If a job is killed before it produces its output, it is rerun with the same configuration and seeds (a restart is not a new look at the data); a restart of the local part resumes from its checkpoint.
