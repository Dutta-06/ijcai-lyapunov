# Plan: interpretable additive models with factorized small-start RBF shape functions (feasibility spike)

Written and committed before any test dataset is scored. UTC time of writing: 2026-11-09 18:24 UTC.
Code: `experiments/exp_gam_spike.py` (shipped to private Kaggle CPU notebooks). Tests: `tests/test_gam_spike.py`.
Report: `results/ml/gam.md` (regenerated from `results/ml/gam_raw.npz` by `--from-raw`). Frozen configurations: `results/ml/gam_frozen.json`.

## 1. Question

Does the paper's mechanism (a Gaussian-kernel network with a factorized readout, weight = mass x saturating amplitude, effective depth >= 3, trained by
plain gradient descent from a small start, switches bumps on one at a time and leaves the unneeded ones dormant) carry to an application where
matching published baselines is feasible: interpretable additive models for tabular data,

    y = b + sum_i f_i(x_i),    f_i(x) = sum_{k=1..K} w_ik exp(-(x - c_ik)^2 / (2 s^2)),    w_ik = nu_ik^2 a_max tanh(alpha_ik / a_max),

alpha_ik(0) = 0, nu_ik(0)^2 = eps, centres c_ik trained (started on a regular grid with a small seeded jitter), one global width s, K = 16 candidate bumps per
feature for EVERY feature and dataset (a ceiling, not a per-feature size)? The number of bumps per feature is then what training leaves active
(|w| > 0.01 in training-sd units of the target; logit units for classification). Primary question: at matched or smaller parameter count, is
it non-inferior to the best neural / boosted interpretable baseline, and better than the same RBF parametrization with a standard random start and
a penalty, with the count chosen automatically?

## 2. Data (fixed by a rule, before any model was run)

Candidate list, in the order written in `CANDIDATES` of the code (public OpenML sets of the kind used in the NAM / EBM / GAMI-Net papers; 12
regression and 12 classification candidates). Rule `choose_datasets`: in that order, per kind, take the first 8 that load and pass the filters
n >= 1000 rows and 3 <= encoded columns <= 40. Encoding (the same for every method): numeric columns as they are (missing -> median, constant
dropped); categorical columns with 2 levels -> one 0/1 column, 3..10 levels -> indicators, more than 10 levels -> dropped; rows without target dropped.
At most 8000 rows per dataset (seeded subsample: sorted first 8000 of `RandomState(777).permutation(n)`; the brief allows 20 000, the cap is for compute).
Per seed s in {0..4} a random 60 / 20 / 20 train / validation / test split (`RandomState(1000 + s)`). The front end fitted on the training rows is shared
by all methods: columns with <= 2 distinct values -> {0, 1}, others -> the training empirical CDF (values in [0, 1]); regression targets standardised by
the training mean and sd (so RMSE is "relative RMSE", in training-sd units); classification: binary, AUC.
Measured candidate sizes (n, encoded width) are in `gam_probe.json` (counts only). Three candidates fail the filters (spambase 57 columns, credit-g 59 columns;
kc1 is not reached). The rule selected:

| dataset | kind | n | encoded d | role |
|---|---|---|---|---|
| california_housing | reg | 20640 | 13 | dev |
| cpu_act | reg | 8192 | 21 | dev |
| MagicTelescope | clf | 19020 | 10 | dev |
| churn | clf | 5000 | 31 | dev |
| Bike_Sharing_Demand | reg | 17379 | 18 | test |
| concrete_compressive_strength | reg | 1030 | 8 | test |
| abalone | reg | 4177 | 10 | test |
| wine_quality | reg | 6497 | 11 | test |
| airfoil_self_noise | reg | 1503 | 5 | test |
| kin8nm | reg | 8192 | 8 | test |
| adult | clf | 48842 | 33 | test |
| bank-marketing | clf | 45211 | 24 | test |
| default-of-credit-card-clients | clf | 30000 | 23 | test |
| higgs | clf | 98050 | 28 | test |
| electricity | clf | 45312 | 14 | test |
| phoneme | clf | 5404 | 5 | test |

Roles: `assign_roles`: within each kind the sorted names are permuted by `RandomState(20261109)`; the first 2 per kind are development datasets, the other 6 per kind
test datasets (4 dev + 12 test). Nothing was chosen after looking at a model's result. (The only model runs before this plan were smoke tests on synthetic toy data.)

## 3. Methods

All fit on train; the validation split is used for stopping / per-run selection only. (EBM, XGBoost and the MLP fit on train + validation and use their
own internal early stopping, so every method sees the same 80 % of the rows.) Test rows are scored once, at the end of the run.

* **ours_gd (the primary row)**: factorized small start, full-batch plain gradient descent, group learning rates (alpha: lr, nu: lr, centres: lr), a_max = 3,
  b fixed (0 for standardized targets; prior logit for classification). Stopping (declared in advance, validation-based): the parameters at the last validation-loss
  improvement of at least 1e-4 (relative), stop when no such improvement came for max(800, t_best) steps (counted only once a bump is active; the path is in logarithmic
  time, an entry can take as long again as everything so far) or at 6000 steps; the number of active bumps and
  stored parameters (2 per active bump: weight and centre) are read at that stop. A step that diverges is halved and restarted (at most 3 times, recorded).
  Development grid: eps in {1e-2, 1e-3, 1e-4} x s in {0.05, 0.08} (lr = 0.5, scaled by 12/d for d > 12 and by 4 for the logistic loss).
* **ours_sgd**: the frozen ours_gd configuration with mini-batches of 512, no tuning of its own. **ours_adam**: Adam on the same parametrization, eps = 1e-4,
  grid s in {0.05, 0.08} x lr in {3e-3, 1e-2}, at most 3000 steps. Both are reported; neither enters the verdict (the conserved quantity of the flow is a
  property of the plain flow; the corpus note says Adam does not select).
* **ctrl_l1 / ctrl_l2 (the key control)**: the SAME shape parametrization (K = 16 bumps per feature, trained centres, one width s) with free weights w_ik from a standard random
  start N(0, 0.1^2), trained by Adam (first draft: plain gradient descent; replaced, see Amendments), with an L1 penalty (proximal step; grid
  {1e-4, 1e-3, 3e-3, 1e-2}) or an L2 penalty (grid {1e-5, 1e-4, 1e-3}) chosen PER RUN on the validation split (more tuning than ours gets). Grid s in {0.05, 0.08}
  on development data. Active bumps and stored parameters as for ours.
* **nam**: neural additive model, one MLP 1 -> h -> h -> 1 (ReLU) per feature, Adam, batch 512, up to 300 epochs, patience 20 on validation, ReLU kinks initialised across [0, 1]; development grid
  h in {16, 48} x lr in {3e-3, 1e-2} x weight decay in {0, 1e-4}. Parameters: all weights, d(h^2 + 4h + 1) + 1. **nam_matched**: the frozen lr / weight decay; single-hidden-layer per-feature networks on the features ours_gd uses, with at least ours' parameter count (see Amendments; this replaces the first draft of this bullet, which used two-layer networks on all features); best of three starts on validation.
* **ebm**: interpret's explainable boosting machine, no interactions, development grid max_bins in {32, 256} x learning rate in {0.01, 0.03}; parameters: values and cut
  points of the step functions after merging equal neighbours.
* **gam_spline**: a spline GAM built from sklearn (cubic B-splines on [0, 1], n_knots in {6, 12, 24} on development data, a ridge / L2-logistic penalty chosen per run on
  validation); pygam is not used (version risk on the current numpy / scipy), recorded as a deviation from the brief's "pygam or similar".
* **xgb**, **mlp** (non-interpretable reference points, not targets): XGBoost (depth in {3, 6}, early stopping on validation), a two-layer sklearn MLP (width in {64, 128}, alpha in {1e-4, 1e-2}).

Tuning: only on the four development datasets, seeds {0, 1}; per method the configuration minimising the mean over datasets of (score / the best configuration's score on
that dataset) (`pick_configs`; the score is relative RMSE or 1 - AUC; ties to the lowest index). It is written to `gam_frozen.json` and committed before the test run. The test
datasets are run once with seeds {0..4}.

## 4. Measures

Test relative RMSE (regression) / AUC (classification); stored parameters; active bumps per feature (and features used); shape complexity (turning points of the shape on a
100-point grid, ignoring reversals under 5 % of the range); stability (importance-weighted mean Pearson correlation of a feature's shape between the five splits);
wall time (Kaggle CPU, noisy, reported only).

## 5. GO / KILL rule (fixed before scoring)

"ours" = ours_gd. Per test dataset (means over the five seeds; score = relative RMSE or 1 - AUC, lower is better):

* **A (accuracy and compactness)**: ours within the margin of the better (by mean test score) of nam and ebm: relative RMSE <= +2 %, AUC >= -0.005 (absolute), AND ours'
  stored parameters <= 0.7 x that comparator's. (Against NAM / EBM at their usual sizes the parameter condition is nearly automatic; this is why B exists.)
* **B (matched parameters)**: ours beats nam_matched (which has at least as many parameters): mean paired improvement over the five seeds > 0 and one-sided paired t-test p < 0.05.
* **C (the control)**: against the more accurate of ctrl_l1 / ctrl_l2 (mean score): ours *wins* if it is more accurate beyond the tie band (0.5 % relative RMSE, 0.002 AUC), or ties
  in accuracy (within the band) and has <= 0.7 x the stored parameters of ctrl_l1; it *loses* if it is less accurate beyond the band; otherwise a tie.

**GO** iff (A holds on >= 8 of the 12 test datasets OR B holds on >= 8 of 12) AND C is a win on >= 7 of 12. Otherwise **KILL**. Rounding: boundaries are inclusive up to 1e-9.
Whatever the verdict, `gam.md` also reports the supplementary table (every method's mean score gap to ours, parameters, time), the Adam and SGD rows, stability and complexity.

If KILL, the report names the cause by the failing clause: A and B fail -> accuracy of the neural / boosted comparators dominates (or ours is not parameter-efficient); C fails ->
the penalised same-parametrization control is as good (the small-start effect adds nothing here); Adam / SGD rows -> whether the effect survives the optimizer.

## 6. What will be dropped / amendments

* A method that errors on a dataset is reported as an error (never silently dropped); a verdict clause needing it is evaluated on the datasets where all needed methods ran and the denominators are stated.
* Any change to this plan after the commit is listed below with its UTC time and is made only before the test run; anything after is labelled POST HOC.
* Known caveats stated now: parameter counts of different model classes are not strictly comparable (the unit is "stored numbers needed to evaluate the shape functions");
  five seeds are five random splits of the same (subsampled) data, so the paired tests are optimistic about between-dataset variation; subsampling to 8000 rows
  favours smooth low-parameter models on the large datasets; one width per method (no per-feature width), centres trained.

### Amendments
* 2026-11-09 18:32 UTC, before any development or test result was read (only synthetic-toy smoke runs had been seen): **nam_matched** changed from "two-layer MLPs on all features, smallest width of the ladder
  with at least ours' parameters" to "single-hidden-layer MLPs (1 -> h -> 1, 3h + 1 numbers) on the features on which ours_gd has at least one active bump (a free feature
  selection for the comparator), h the smallest of {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64} with at least ours' total parameter count (n_used (3h + 1) + 1 >= ours),
  best of three starts on validation". Reason: the smallest two-layer per-feature net has 13 numbers per feature on EVERY feature, so a network "matched" to a handful of bumps
  could not be built and trained poorly (smoke test), which would be a weakened comparator. The tuned `nam` (two-layer, all features) is unchanged.
* 2026-11-09 19:53 UTC, after the development results of the baselines (and a development-only exploration of ours_gd and of the first control on cpu_act and MagicTelescope), before any test result of ours / the controls was read
  (the baselines' test run had finished and was not opened):
  1. **The control is retrained with Adam.** In the exploration the plain-gradient-descent control (random start, proximal L1 / ridge, 4000 steps) was still improving at its step cap and
     reached relative RMSE 0.293 (L1) / 0.302 (L2) on cpu_act, against 0.17-0.25 for ours and 0.14-0.15 for NAM / EBM / spline: an under-trained comparator, which the brief forbids. `ctrl_l1` / `ctrl_l2`
     now use Adam (lr 1e-2, up to 3000 steps, patience 300, betas 0.9 / 0.999), the L1 penalty as a proximal step in Adam's diagonal metric (soft threshold lr x lambda / sqrt(v)),
     the L2 penalty in the gradient; penalty grids L1 {1e-4, 3e-4, 1e-3, 3e-3, 1e-2}, L2 {1e-5, 1e-4, 1e-3, 1e-2}, chosen per run on validation; width grid {0.03, 0.05, 0.08} on development data.
     The first, gradient-descent controls stay in the development records under the names `ctrlgd_l1` / `ctrlgd_l2` (not used in any verdict).
  2. **Second development round for ours_gd** (the plan's grid lets the factorized unit run at most 6000 plain-gradient steps, and the exploration shows accuracy still rising with the step size and
     with finer widths): configurations 6-11 are added, eps in {1e-2, 1e-3} x width in {0.03, 0.05, 0.08} at step size 1.0; configurations 0-5 are the first round. All twelve are in the development
     records and `pick_configs` chooses among all twelve; the step cap (6000) and the stopping rule are unchanged, so the compute of ours is not increased.

### Run log (UTC, from the commits / notebook starts; nothing below was changed after the corresponding test run began)
* 2026-11-09 ~18:55: frozen baseline configurations committed (`gam_frozen_baselines.json`, commit 8b1e68a); test run of nam, ebm, gam_spline, xgb, mlp on the twelve test datasets started ~18:56 and finished ~19:20.
  Their test results were not opened until every method had been frozen.
* 2026-11-09 21:36: frozen controls (Adam) and ours_adam committed (`gam_frozen_controls_adam.json`, commit bf6b560); their test run started ~21:38.
* 2026-11-09 22:00: all development records and `gam_frozen.json` committed (commit 71f79de): ours_gd = eps 1e-2, width 0.05, step 1.0 (configuration 7 of 12; chosen by the accuracy-only rule, which
  prefers the fastest, least selective start); ours_sgd inherits it; test run of ours_gd, ours_sgd, nam_matched started ~22:01 in two notebooks.
* 2026-11-09 22:41 / 23:11: the second half of the test datasets for ours / nam_matched started and finished (the first half ran 22:01-22:40; one notebook at a time per half because of the two-notebook limit);
  the control / Adam test run (started 21:38) finished 2026-11-10 00:18. The test results of ours and the baselines were first opened at ~23:15, after every configuration had been frozen and the first
  seven test jobs finished; no configuration, grid or rule was changed after that. The report's POST HOC section and the development summary were added after the results were read (labelled as such in `gam.md`).
