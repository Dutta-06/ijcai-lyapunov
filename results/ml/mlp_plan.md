# Plan (pre-registration): do the two practical messages (depth of the readout, per-coordinate optimizers) hold outside RBF / KAN, in a one-hidden-layer MLP?

Written 2026-11-10 23:31 UTC, before any test target or test run of this experiment (a few development runs on the development target, teacher id 901, were already running; they are
not part of any verdict and the test teachers 1-30 have not been generated or looked at). Committed before the test run; later changes are listed under *Amendments* with their UTC time.
Code: `experiments/exp_mlp_dormancy.py`; tests: `tests/test_exp_mlp_dormancy.py`; mutation spec: `tools/mutants/mlp_dormancy_own.py`. numpy only.

## Question
Two reviewers asked for a test of the paper's two practical messages in a standard ML testbed:
**(R)** a readout of effective depth >= 3 lets small-initialisation gradient training leave unneeded units dormant (depth 2 only logarithmically in eps, depth 1 not at all);
**(O)** per-coordinate-normalising optimizers (Adam, sign descent) erase this, plain gradient descent and a shared-second-moment variant keep it.
A null outside RBF is a valid and informative outcome and will be reported as such.

## Models
* Teacher `f*(x) = sum_{k<=K} b_k tanh(w_k . x)`, `x ~ N(0, I_d)`, `|w_k| = 1`, pairwise `|cos(w_j, w_k)| <= 0.5` (rejection), `b_k = +-U(0.5, 1.5)` with random sign;
  K in {2, 4, 8}, d in {10, 30}; label noise sd 0.05 and 0 (the same teacher and student initialisation at both levels); n = 4000 training rows, 1000 noisy validation rows (early stopping), 20000 clean test rows.
* Student `f(x) = sum_{i<=H} a_i tanh(w_i . x)`, H = 128, no biases. Readout of effective depth p: p = 1: `a = u`; p = 2: `a = u^2 - v^2`; p = 3: `a = u^3 - v^3`. Initial scale `s`, `s^2 = eps in {1e-2, 1e-4, 1e-6}`:
  inner weights `w_i` with entries `N(0, s^2 / d)` (norm about s), `u = v = s` (p >= 2, so a = 0), `u_i = s N(0, 1)` (p = 1). Loss 0.5 * MSE, full batch. float32 forward and backward, float64 parameters.
* Optimizers (all on the flat vector of all parameters): `gd` (one rate); `adam` (default betas 0.9 / 0.999, eps_adam 1e-8); `shared_v` (Adam with ONE scalar second moment, the EMA of the mean of g^2 over all
  parameters); `sign` (sign descent, fixed rate); `adam_big` (Adam with eps_adam = 0.1). Runs: GD at depths 1, 2, 3; every other arm at depth 3 only (7 (arm, depth) cells x 3 eps). Declared: depth is varied for GD only (R is a GD statement); depth 3 is the depth of (O).
* Step budget **80000** full-batch steps; validation every 50 steps; the iterate of lowest validation MSE is kept; armed once the validation MSE is below 0.9 times its initial value; stop when armed and the time since the last
  relative improvement (above 0.1%) exceeds max(5000, the step of that improvement) ("doubling rule"). A run that never arms runs to the budget. Flags stored per run: converged / armed but budget exhausted / never left the start /
  diverged (validation MSE not finite or above 1e6 times its initial value; the best iterate before divergence is scored). **No run is dropped**.

## Development and test (disjoint)
* **Development** (rate tuning, eps sets; nothing else): teacher id 901 in each of the six cells (K, d), noise 0.05, data and initialisation seeds from the range starting at 7,000,000 (`STAGE_SEED['dev']`): 6 runs per configuration.
* **Test**: teacher ids 1-30 in each (K, d) cell (30 per cell, 180 teachers), noise 0.05 and 0 (360 target-noise combinations), data seeds from the range starting at 1,000,000. Teacher ids and seed ranges are disjoint from the development ones;
  the teacher is a function of (K, d, id) only, so no test teacher equals a development teacher. Test size: 360 x 7 x 3 = **7560 runs**.
* **Rates**: 4 values per (arm, depth, eps) (`LR_GRID`: gd depth 1 {0.1, 0.25, 0.6, 1.5}, gd depth 2 and 3 {0.05, 0.15, 0.4, 1.0}, adam {3e-4, 1e-3, 3e-3, 1e-2}, shared_v {3e-3, 1e-2, 3e-2, 1e-1}, sign {1e-4, 3e-4, 1e-3, 3e-3},
  adam_big {0.03, 0.1, 0.3, 1.0}), the same number of values and runs for every arm; chosen by the lowest mean over the development runs of log10(validation MSE / initial validation MSE) (ties: smaller rate).
  If a chosen rate sits at the edge of its grid, that configuration's grid is extended by one value (a factor of about 2.5 / 3) on that side, once, for every arm alike (an amendment with the time of the extension), and re-chosen.
* **Eps sets** (rule, applied on development): `E_O` = the eps for which depth-3 GD at its chosen rate fits (relative clean test MSE < 0.2) in at least half of the development runs; `E_R` = the eps in `E_O` with eps <= 1e-4.
  Eps outside `E_O` are reported in every table (unconditionally) but do not enter the verdicts: where depth-3 GD does not leave its start in the step budget, "dormancy" says nothing about the units that are needed.
* The merge cosine 0.95, the activity threshold 0.05, the sweep grids, the margins below and the budget are fixed here and are not tuned.

## Outcomes
Active neuron: strength `|a_i| |w_i| > tau` (declared tau = 0.05; the teacher's strengths are |b_k| >= 0.5). Merging: single-linkage components of the active neurons' oriented directions `sign(a_i) w_i / |w_i|` with cosine
> 0.95. Count error = |merged count - K|. Relative test MSE = MSE against the clean teacher / teacher variance (test rows). A run "fits" if that is below 0.2. Teacher recovery: a teacher direction (oriented `sign(b_k) w_k`) is recovered if some
active neuron has oriented cosine > 0.9 to it. Sweeps (stored for every run): tau in {0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5} x merge cosine in {0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, no merging}.

## Registered predictions and decision rules
All contrasts are paired over runs (same target, noise and initial inner weights), averaged over the noise levels within a target, with a 95% percentile cluster bootstrap over targets (20000 resamples, stratified by (K, d), cells weighted equally).
A **count-error gain** is established if its interval is above 0 and its mean is at least 0.5 neuron.
* **(R)** For every eps in `E_R`: the count error of depth-3 GD is smaller than that of depth-2 GD and than that of depth-1 GD (both gains established). Held if both hold at every eps in `E_R` (and `E_R` is not empty).
* **(R2)** Depth-1 GD forms about H neurons: the median number of active neurons before merging is at least H / 2 = 64 at both eps <= 1e-4. (This is a stated expectation with a stated threshold; if it fails, it is reported as failed.)
* **(O)** At depth 3, for every eps in `E_O`: count error of Adam minus that of GD, of sign descent minus GD, of Adam minus shared-v, of sign descent minus shared-v: all four gains established (declared tau, merge cosine).
* **(O-sweep)** The same four contrasts over the whole tau x cosine sweep (2000 resamples per cell): at every eps in `E_O` at least half of the 64 cells have the comparator worse (interval above 0). The earlier KAN test showed the Adam gap closing at a larger
  merge radius; the full sweep is printed whatever it shows (both the `+` and `-` shares).
* **(T)** For every eps in `E_O` and each of the six comparators (GD depth 1, GD depth 2, Adam, shared-v, sign, large-eps Adam): the paired difference of the relative test MSE (depth-3 GD minus the comparator) has its interval's upper end at most 0.01
  (1% of the teacher variance): depth-3 GD is not worse by more than that.
* Secondary, not in any verdict: the same contrasts on the runs where both arms fit; count error by K; convergence counts; teacher recovery; adam_big and shared_v descriptively.

## What counts as what
* Each prediction is *held* or *failed* by the rule above; a held (R) and a failed (O) (or the reverse) are reported as such. A failure is reported as a failure of the paper's message in this testbed, with the numbers.
* What will be dropped: nothing. Diverged, never-armed and budget-exhausted runs stay in every table (they are listed in the convergence table). The numbers of the verdicts are those of the declared setting; sweeps are shown in full.
* Anything done after the test run is labelled POST HOC.

## Compute
Development and test runs are made on CPU notebooks (numpy, 4 workers each, shards of the flat job list) and the raw jsonl parts are merged into `results/ml/mlp_dev_raw.npz` / `mlp_raw.npz`; the report regenerates from the raw with `--from-raw`
(byte-identically). The development report is `results/ml/mlp_dev.md`; the frozen choices are `results/ml/mlp_freeze.json` (rates, `E_R`, `E_O`), committed before the test run.

## Amendments
1. Clock correction, 2026-11-10 23:25 UTC. The header time above (23:31) was a mis-estimate written before reading the clock; the plan was committed at 23:24:47 UTC (commit 45e52bb), before any test run. No rule changed.
2. 2026-11-11 03:50 UTC (commit time), after the development runs and before any test run, from development data and from the measured run times only. (a) **Measured cost.** On the CPU notebooks one full-batch step takes about 12 ms
   (K = 8, d = 30; 1.9 ms on the local core): a converged development run took about 100 s, a run that never leaves its start (80000 steps) about 1000 s. The registered test (30 teachers per cell, 7560 runs) would need
   about 80 CPU-hours on the local core. **Test size is therefore reduced to** teacher ids 1-10 per (K, d) cell (10 per cell, 60 teachers) at noise 0.05, noise 0 for the first 3 of them (18 teachers), and `(gd, depth 3, eps 1e-6)` only for teacher 1 of
   each cell at noise 0.05 (6 runs), because on the development targets that configuration never left its start at any of the 4 rates (6 of 6 runs never armed at every rate), so repeating it 60 times adds nothing but compute: **1566 test runs**.
   Contrasts at eps = 1e-6 involving it are computed on the targets run in both arms (`Table.align`). This is below the "30 targets per cell" asked for; the reason is compute, stated here before the test. (b) **Rate extension** (the rule of the Development
   section, applied once): the first development round (4 rates, 504 runs, 3 notebooks) chose a grid-edge rate for sign (all 3 eps), shared_v (eps 1e-2, 1e-4: lowest edge; eps 1e-6: highest edge) and adam_big (eps 1e-6: highest edge); one more rate beyond that edge was run for
   each (`LR_EXTRA`: sign 3e-5, shared_v 1e-3 / 1e-3 / 0.3, adam_big 3.0; 42 runs, local core) and the rate re-chosen: shared_v 0.001 (eps 1e-2, tie with 0.003 within 0.01, the smaller taken), 0.003 (eps 1e-4), 0.1 (eps 1e-6; 0.3 is equal); sign 3e-5 (eps 1e-2, tie with 1e-4), 1e-4 (eps 1e-4 and 1e-6); adam_big keeps 1.0 (3.0 never armed). Edge rates left after the single extension (reported in `mlp_dev.md`): gd depth 3 at eps 1e-6 (every rate the same: never leaves the start), shared_v at eps 1e-2 and sign at eps 1e-2 (the scores of the edge rate and its neighbour are equal to 0.01).
   (c) **Frozen** (`results/ml/mlp_freeze.json`, committed before the test): the rates of `mlp_dev.md`; the eps sets by the registered rule: depth-3 GD fits in 6/6 development runs at eps 1e-2 and 1e-4 and in 0/6 at 1e-6, so **`E_O` = {1e-2, 1e-4}, `E_R` = {1e-4}**. The step budget,
   patience, merge cosine, thresholds and every margin are unchanged. (d) **Run plan:** the 1566 test jobs are cut in 5 shards (job index modulo 5): shards 0-3 on four notebooks (neutral titles), shard 4 on the local core; each shard appends to its own jsonl. **Run log:** amendment 2 was committed 2026-11-11 03:46 UTC (5e1bc17, development raw 8d1cd71), before any test run. Test started 03:46 UTC: shards 0, 1, 2 on notebooks `grid-sweep-1..3` (4 workers each), shard 4 on the local core (1 worker); shard 3 could not be launched (the notebook session cap of 5 was reached by other jobs) and is run as soon as a slot is free, else locally after shard 4. A restart or relaunch of a shard is not a new look at the test data: the same frozen configuration is rerun and the shard's jsonl resumes by job key.
   Run log, continued: shard 4 (local) finished 2026-11-11 08:09 UTC (313 jobs, 4.4 h). The notebook session cap kept shard 3 from launching on a notebook for 4 h (launch retried every 3 min); it was started on the local core at 08:10 UTC (same frozen configuration, same job list, shard 3 of 5). Notebooks `grid-sweep-1..3` (shards 0-2) still running at 08:10 UTC.
