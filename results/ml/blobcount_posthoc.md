# POST HOC diagnostic (made after the test result): the frozen paths with an oracle stop

The frozen path of `ours` (dual order) and of `greedy` is run on the test patches again; the block of the path whose count is closest to the true count (it reads the ground truth: a diagnostic, not a method) is taken. This separates "the path does not contain the right number of units" from "the stop rule picks the wrong block". Generated from `results/ml/blobcount_posthoc_raw.npz` and `results/ml/blobcount_raw.npz` (`--from-raw`).

| method | stratum | count MAE at the frozen stop | count MAE at the oracle block | recall of true centres at the oracle block | best recall anywhere on the path | mean blocks on the path | MAE of the last block |
|---|---|---|---|---|---|---|---|
| ours | all | 4.833 | 0.153 | 0.743 | 0.841 | 18.7 | 10.547 |
| ours | crowded | 7.853 | 2.294 | 0.666 | 0.714 | 24.4 | 5.588 |
| greedy | all | 4.663 | 0.027 | 0.791 | 0.895 | 42.4 | 30.537 |
| greedy | crowded | 7.559 | 0.324 | 0.811 | 0.908 | 47.0 | 21.779 |

Oracle-stop crowded MAE, ours minus greedy: +1.971, interval [+0.823, +3.350] (paired bootstrap over images).

