# Peak counting spike: test report (500 held-out CPSC 2019 records)

Plan and rule: `results/ml/peakcount_plan.md` (committed before any test record was scored). Data: CPSC 2019 training records, expert R-peak annotations; target = number of annotated R peaks per 10 s record (mean 14.79, range 6 to 38). Front end `abs`. Every configuration below was frozen on the development records (`results/ml/peakcount_freeze.json`).

## Verdict: **KILL**

* Win condition: ours (`iss`) against the best tuned classical counter frozen on development (`mf`): mean absolute count error difference (baseline minus ours) -2.556, 95% bootstrap interval [-2.938, -2.184] -> not met.
* Non-inferiority to the greedy sliding matching pursuit (`omp`): ours minus omp +0.056, 95% interval [-0.254, +0.390], margin 0.10 beats per record -> not met.

## Count error of every family at its frozen configuration (test records)

| family | frozen configuration | mean abs error | exact count | mean signed error | secondary: position F1 (tol 0.10 s) |
|---|---|---|---|---|---|
| find_peaks (prominence x distance) | 0.25, 0.3 | 1.428 | 0.486 | -0.296 | 0.920 |
| find_peaks, prominence relative to the record | 0.3, 0.3 | 1.426 | 0.470 | -0.330 | 0.919 |
| matched filter + find_peaks | 0.07, 0.12, 0.2 | 1.168 | 0.494 | -0.192 | 0.935 |
| forward selection on the grid + stopping rule | width 0.05 s, drop 0.02 | 3.444 | 0.128 | -1.460 | 0.803 |
| matching pursuit with sliding + stopping rule (greedy control) | width 0.05 s, drop 0.01 | 3.668 | 0.236 | +1.572 | 0.854 |
| sliding inverse scale space path + stopping rule (ours) | width 0.05 s, drop 0.01 | 3.724 | 0.224 | -0.444 | 0.829 |

## Paired differences of mean absolute error against ours (baseline minus ours; positive = ours better), 95% bootstrap intervals

| baseline | difference | 95% interval |
|---|---|---|
| find_peaks (prominence x distance) | -2.296 | [-2.690, -1.914] |
| find_peaks, prominence relative to the record | -2.298 | [-2.706, -1.908] |
| matched filter + find_peaks | -2.556 | [-2.938, -2.184] |
| forward selection on the grid + stopping rule | -0.280 | [-0.656, +0.080] |
| matching pursuit with sliding + stopping rule (greedy control) | -0.056 | [-0.390, +0.254] |

## Where the errors are (mean absolute error by true count)

| true beats | records | fp | fp_rel | mf | step | omp | iss |
|---|---|---|---|---|---|---|---|
| 0-10 | 50 | 1.30 | 1.34 | 1.06 | 2.54 | 3.52 | 2.76 |
| 11-14 | 253 | 1.12 | 1.20 | 0.86 | 2.80 | 3.77 | 3.33 |
| 15-19 | 130 | 1.36 | 1.18 | 1.21 | 3.86 | 4.04 | 4.40 |
| 20-99 | 67 | 2.82 | 2.81 | 2.34 | 5.73 | 2.66 | 4.63 |

Notes. Errors mix the model-order error with the front-end error, identically for every method. Positions are scored after removing the shared offset -0.004 s found on development; the F1 is secondary and not part of the verdict. The "ours" family is the limit form of the factorised small-start network (`baselines.iss_path`), not a trained torch network.
