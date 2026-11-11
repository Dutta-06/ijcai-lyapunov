# REVIEW_roundML: referee report on the four ML feasibility spikes (KAN, additive models, nuclei counting, ECG peak counting)

Archived verbatim from the referee's hand-back (2026-11-10). The referee was read-only for the repository; its scratch code is outside the repository.

REFEREE REPORT, four ML spikes. Scratch code: referee_ML_work/{kan1.py,blob.py,pk.py}. The repo was not modified. Everything was recomputed from the stored raw files. torch is absent, so no model was re-run.

HEADLINE
- KAN: by its own pre-registered rule the verdict is KILL (S met, R not met, 1 of 8 sets). It is not a positive result as a whole.
- The "72% vs 28%" count result does not survive a simple stronger control. The control's surplus units are mostly coincident duplicates, and merging them closes most or all of the gap.
- GAM, nuclei and ECG: the KILLs reproduce and are not artefacts.

1. KAN (kan_plan.md, kan_raw.npz, kan.md)

Findings
1.1 BLOCKER for the wording "72% vs 28%". The control's 28.3% is almost entirely the empty edges.
- 84 of the 300 edges have K=0, and both methods are exact on all 84 (84/300 = 0.28).
- On the 216 non-empty edges the control is exact on 1 (0.5%); ours is exact on 61.6%.
- The control's 8.67 surplus units per run sit entirely on K>0 edges (0.00 on K=0 edges), about 2.4 extra units per non-empty edge.
- Duplicates are the cause. 40.6% of adjacent active-unit gaps in the control are below 0.1 and 54% are below 0.3; for ours, 4.9% and 12.6%. The control's precision at tolerance 0.2 is 0.88 (0.96 at 0.5), so the extra units are copies sitting on true bumps.
- Mechanism: L1 on same-sign coincident units does not penalise splitting an amplitude between them. A larger lambda cannot fix this, which is why the oracle penalty gives 0.287. The missing step is merging, not a stronger penalty.

1.2 Merging the control's coincident units, using the stored centres of the validation-chosen control. Both methods are merged by single linkage at the same radius delta. Gain is ours minus control, cluster bootstrap over the 10 targets, 20,000 resamples.

| delta | control exact | ours exact | gain [95% interval] |
|---|---|---|---|
| 0 (published) | 0.283 | 0.723 | 0.440 [0.363, 0.517] |
| 0.10 | 0.563 | 0.750 | 0.187 [0.10, 0.28] |
| 0.20 | 0.687 | 0.783 | 0.097 [0.03, 0.18] |
| 0.30 | 0.767 | 0.810 | 0.043 [-0.03, 0.12] (t-interval [-0.05, 0.14]) |
| 0.40 | 0.840 | 0.830 | -0.010 |
| 0.50 | 0.910 | 0.863 | -0.047 [-0.10, 0.007] |

- By delta = 0.3 the S criterion (gain at least 0.10 with an interval above 0) fails.
- The generator's true centres are at least 1.0 apart, so any delta up to 0.5 is a legitimate "distinct bump" radius. Initial width and grid spacing are 0.556, so a natural rule gives delta of about 0.3 to 0.55.
- With merging at 0.5 and the best lambda chosen per run by the truth, the control reaches 0.957.
- Caveat in the other direction: the control has no merge step in its protocol, and delta is my choice.
- Not tested, and no torch to run: iterative pruning, group-L1 or Hoyer penalties, pruning with a larger tau. Larger tau does not help: the control is 0.300 at tau = 0.05.

1.3 MAJOR: accuracy and recall go the other way.
- Clean-test RMSE is 0.068 for ours against 0.039 for the control (ours is worse in 55% of runs).
- Recall of true bumps is 0.930 against 0.986; ours under-counts 0.08 units per run.
- The count is bought with a worse fit. Do not describe it as better model selection.

1.4 MAJOR: only 14 of 60 runs (23%) are exact on all 5 edges for ours; the control has 0 of 60.
- Ours is exact on only 46% of K=3 edges.
- Ours' seed sd is 0.099 against 0.005 for the control.
- Ours hit the 12,000-step cap in 32 of 60 runs, so the count partly reflects the step budget.

1.5 MAJOR: hyper-parameters.
- Ours (eps, lr) were chosen on the 3 development targets (12 runs) by the true exact-edge share, which no practitioner has. That is development only; the freeze is committed at 22:53, before the test.
- The control got only 2 to 3 learning rates and a penalty chosen by validation MSE or by the truth.
- Sensitivity of ours on development: exact share 0.23 to 0.70 across the 10 configurations (eps 1e-2 gives 0.28 to 0.50; lr 2 gives 0.47). The count is steered by eps and lr, so "no hyper-parameter to set" overstates.

1.6 MAJOR: free advantage to ours. The truth is a sum of Gaussian bumps with widths 0.25 to 0.45, the same family as the kernel, with a minimum separation of 1.0. The scope is realisable synthetic additive data only. There is no "true count" on real data.

1.7 MINOR.
- tau is applied identically: both go through the same `kan_summary`. Ours at tau 0.01 / 0.02 / 0.05 is 0.68 / 0.72 / 0.81; the control is 0.28 / 0.28 / 0.30.
- The bootstrap is correct. My own 20,000-resample cluster bootstrap over targets gives 0.440 [0.363, 0.517]. All 10 per-target gains are positive.
- The tables kan.md and kan_raw.npz reproduce (ours 0.7233, control 0.2833). 52 of 53 tests pass; 1 skipped for lack of torch.
- The development report's headline "ours 0.500" row is the first configuration, not the frozen one. It says so in the text, which is OK.
- The two-layer and Adam statements are development-only and labelled as such. The pure two-layer small start stayed at the mean predictor (RMSE 0.93 to 1.05) on 3 development sets; the Adam configurations scored 0.23 to 0.28 on 12 runs. Neither is a test claim, and the Adam grid is only 4 configurations.
- The hybrid [d, 4, 1] architecture standardises the second layer and tunes its rate ratio (l2_ratio), so it is not the small-start depth of the theory.

1.8 Amendment timeline matches git: plan ec79f56 at 18:12; amendment 1 at 18:37; amendment 2 at 20:59; amendment 3 at 22:24; freeze e0587f5 at 22:53; amendment 5 at 23:06; kan_sets.json at 00:22; raw at 01:18. Two disclosed gaps:
- The amendment 5 smoke run looked at 2 test jobs: target 1, noise 0.1, seeds 0 and 1, where ours scored 0.90 against 0.60 for the control. Nothing was changed afterwards, and these jobs are inside the 60.
- Amendment 6 swapped treasury and weather_izmir. This is verified against the rule (ascending id order, ranks 3, 6, 9 are development). Treasury was tuned on, so keeping it out of the test is correct. The effect on R is nil: R fails anyway, and the 1-of-8 non-inferior count is the same either way (treasury's hybrid ratio on development was 1.42).

Cheapest fix: add a "merge units within a radius" post-processor to the control, report the gain as a function of delta, and report the exact share on K>0 edges separately.

Allowed paper wording: "On 10 synthetic additive targets with Gaussian-bump truth (60 runs), one-layer small-start gradient descent returns the exact number of active units on 72% of edges without post-processing (62% of non-empty edges). An L1-plus-magnitude-pruning control with the same units returns duplicated units on the same bump and is exact on 28% of edges (the empty ones). After merging units whose centres lie within 0.3 the shares are 81% and 77%, and beyond 0.4 the control is level or ahead. The count is bought with a worse clean-function fit (RMSE 0.068 against 0.039). On 8 real tabular sets the architecture used was non-inferior at matched parameters on 1; the pre-registered verdict is KILL."

Must not claim: that the method selects the model order better than sparsifiers in general; "no hyper-parameter"; superiority on real data; a stronger result for the two-layer or Adam variants; the gain 0.44 without the merge caveat.

2. GAM (additive models)

- The KILL reproduces. Using the module's own `criteria_for_dataset` on the 660 raw records, I get A = 5/12, B = 7/12, C = 3/12 wins. The rule needs (A >= 8 or B >= 8) and C >= 7, so both clauses fail. B is close (7 against the 8 needed).
- The post hoc table (tie bands, control subsets, A margin 1.5x with no parameter condition) never reaches C >= 7; the maximum is 6.
- Not an artefact:
  - The 8,000-row cap, K = 16 bumps and the shared width apply identically to the control. Ours' width grid {0.03, 0.05, 0.08} has its optimum at 0.05, an interior value.
  - 13% of ours' runs hit the 6,000-step cap (median 2,512). That is a minor handicap, and it is small.
- MAJOR design caveat. The frozen "ours" is eps = 1e-2, the largest eps tried and the grid edge. Development accuracy worsens monotonically with smaller eps (ratio 1.01, 1.08, 1.57 for 1e-2, 1e-3, 1e-4). The configuration that wins on accuracy is the least small-start one, so this spike barely tests the dormancy mechanism. eps above 1e-2 was not tried.
- MAJOR: "level with a tuned NAM at 0.3% of its parameters" is true only in a narrow sense.
  - 128 / 38,913 = 0.33%, but the frozen NAM is h = 48 with 46,820 parameters in development, while h = 16 (6,020 parameters) was equally good (ratio 1.0054 against 1.0038). Against h = 16, ours is 2.1% of its size.
  - Against the pre-registered comparator, the better of NAM and EBM (EBM wins 11/12), ours is non-inferior on only 5 of 12.
  - Against NAM alone (my post hoc count with the plan's margins: +2% relative RMSE, -0.005 AUC) ours is inside the margin on 9 of 12. It is outside on abalone (+3.2%), credit-default (+0.008 AUC) and higgs (+0.014 AUC).
  - The only apples-to-apples neural comparator, nam_matched (146 parameters), is beaten on 7 of 12 and is worse by 0.034 on average. A cubic-spline GAM with 203 parameters is level with or slightly better than ours (-0.003 on average; 3 of 12 better beyond the band).
- Parameter counting: ours counts 2 per active bump (weight and centre), without the shared width and bias. NAM counts all weights. EBM counts values plus cut points. The spline counts coefficients. These are not strictly commensurable.
- Control: ctrl_l1 has 1.5x ours' parameters (ratio 0.66) and is better beyond the band on 4 of 12 datasets. Ours is not worse beyond the band on 8 of 12.
- Allowed wording: "Pre-registered verdict KILL (A 5/12, B 7/12, C 3/12). Mean relative accuracy is within 0.001 of our NAM implementation, which uses 300 times more stored numbers than the active bumps of ours. It is 0.026 behind EBM at about 1/26 of its parameters (128 against 3,347), and a spline GAM is level at 1.6 times the parameters." Label the NAM-only margin count as post hoc.
- Must not claim: "level with NAM at 0.3% of parameters" without the h = 16 and spline context; a small-start benefit (the frozen eps is the largest, i.e. no dormancy); "beats the control".

3. Nuclei counting

- Recomputed from blobcount_raw.npz with my own bootstrap by image (344 patches, 150 images; 68 crowded patches from 26 images).
- Crowded MAE: ours 7.853, plm 2.809. The difference (ours minus classical) is +5.04 [1.74, 8.85] against plm, +4.74 [1.57, 8.24] against log, +4.87 [1.56, 8.91] against ws. This matches the report.
- Base sample: ours 4.833, log 2.830. The difference is +2.00 [1.15, 2.81]. The verdict reproduces. Taking image-level means first gives +3.39 [1.21, 6.05].
- Amendment 2 (threshold 25 to 15):
  - The commit history is consistent with the claim: plan 7974836 at 08:27 with the 25 rule; code plus amendment 2 in e75e57e at 08:54; the first dev result in fd95e32 at 09:15; test raw at 09:26.
  - The ground-truth counts of the test tiles were tabulated before the decision (disclosed in the plan); no method had been run.
  - At 25 or more the picture is worse for ours (21.1 against 4.8, 11 patches from 5 images).
  - The commit message says "09:05" in the original e75e57e text; it was corrected later. Cosmetic.
- ISS not handicapped by the cap: the maximum number of detections is 34 against a cap of 120, and the no-cap re-run gives 4.76 / 8.28 (no help). An oracle stop on the path gives 0.153 (all) and 2.29 (crowded). So the path contains the right counts and the BIC stop with its noise calibration fails.
- Count insensitivity: ours returns about 10 nuclei for every true count (mean detections 9.4, 9.4, 10.1, 12.9, 14.7, 11.7 for true-count bins 1-3, 4-6, 7-10, 11-14, 15-24, 25+). The correlation of detected with true count is 0.37 (log 0.79). On the base sample a constant predictor of 4 gets MAE 3.74, better than ours.
- Alternative stops (a count-calibrated or CV stop) were never tried, and the BIC calibration uses the residual of the path's last block, which depends on the early stop. I cannot say whether a better-calibrated stop would flip the result; this is a possible handicap I could not check.
- Allowed wording: "Pre-registered KILL. The sliding-path count is 2.8 times worse than the best tuned classical detector on crowded patches (7.85 against 2.81 MAE, interval of the difference 1.7 to 8.9) and worse on all patches. With an oracle stop that reads the truth the same path reaches MAE 0.15, so the failure is the stopping rule (post hoc)." Mention the reference results: StarDist 1.0 overall; Cellpose 9.9 on crowded.
- Must not claim: that the path "finds the right count" or the oracle-stop numbers as a method result; any advantage on crowded patches.

4. ECG peak counting

- Recomputed from peakcount_raw.npz (E and the path arrays) with my own record-level bootstrap, 20,000 resamples.
- MAE: ours 3.724, matched filter 1.168, find_peaks 1.428, omp 3.668. Matched filter minus ours = -2.556 [-2.94, -2.19]; ours minus omp = +0.056 [-0.26, +0.38] (N requires an upper end below 0.10, so not met). Matches the report. The ids in the raw equal the test split.
- Ours is worse than predicting the median count (constant median MAE 3.27 on the test labels).
- Edge-of-grid: ours' optimal path width is the largest tried (0.05 s; widths above that were dropped for compute before any accuracy was seen). Development MAE falls monotonically 4.29, 3.83, 3.64 with width. A wider grid would help somewhat but cannot close a gap of 2.4.
- No other configuration in the grid beats the frozen one on the test labels (oracle check), so the KILL is not a tuning artefact.
- Timing: the plan was first committed in WIP 3f5b670 (09:09) with a time placeholder; the "09:26" header was a typo (plan commit b532f64 at 09:23, amendment 2 at 09:23:40); test raw at 09:56. Development results were seen before the plan commit (the plan says so). The plan is pre-test, not pre-development.
- Disclosed looks at labels: aggregate beat statistics over all 2,000 records, and prototype fits on records 1, 2, 3, 500, 1200 and 1500. Two of these, 500 and 1200, are TEST records. They were run for timing only and not scored, but the plan does not say they are test records. Effect negligible.
- Allowed wording: "Pre-registered KILL: on 500 held-out ECG records the path count is 2.6 beats worse than a tuned matched filter (3.72 against 1.17 MAE, paired interval of the difference 2.2 to 2.9) and no better than a constant predictor. The path width was at the edge of a compute-limited grid." Do not claim any advantage for the path on ECG data.

Verified (how)
- Reproduced from the stored raw: KAN exact shares (0.723 / 0.283), the gain interval, the duplicate statistics, the per-K shares, tau variants and the all-5-exact counts; the GAM A/B/C counts via the module's `criteria_for_dataset`; blobcount crowded and all MAE with a by-image bootstrap; peakcount MAE and bootstraps by record.
- Checked against git log: plan, amendment and freeze times for all four spikes.
- 52 KAN tests pass.

Could not check
- Any merge, iterative-pruning or group-penalty control re-trained with torch (no torch). The merge test uses stored active centres only, without widths or weights.
- The Kaggle test notebooks were run on the committed code, but I could not confirm that.
- Whether a better-calibrated stop in the nuclei spike, or wider kernels in the ECG spike, would change the ranking.
- Real-data raw sets (not in the repo), and the correctness of the GAM NAM and EBM implementations against published ones.
