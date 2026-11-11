# Peak counting spike: development summary (150 development records)

Front end: `abs`.  Mean absolute count error of the frozen configuration of each family (lowest over its 42 configurations):

| family | frozen configuration | dev mean abs error | dev mean signed error |
|---|---|---|---|
| find_peaks (prominence x distance) | 0.25, 0.3 | 1.553 | -0.100 |
| find_peaks, prominence relative to the record | 0.3, 0.3 | 1.647 | -0.327 |
| matched filter + find_peaks | 0.07, 0.12, 0.2 | 1.453 | -0.080 |
| forward selection on the grid + stopping rule | width 0.05 s, drop 0.02 | 3.467 | -1.480 |
| matching pursuit with sliding + stopping rule (greedy control) | width 0.05 s, drop 0.01 | 3.540 | +1.740 |
| sliding inverse scale space path + stopping rule (ours) | width 0.05 s, drop 0.01 | 3.640 | +0.027 |

Best tuned classical counter (lowest development error among fp, fp_rel, mf, step): **mf**.  Shared position offset: -0.004 s.
