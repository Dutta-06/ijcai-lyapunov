# Experiments that were started and stopped

Both were stopped on request before any held-out result existed. They are listed so the gap in the evidence is visible. The plans are in the repository as written.

## 1. Separation sweep on the trained network (`experiments/exp_sepsweep.py`, plan `results/synthetic/sepsweep_plan.md`)
Purpose: the separation condition (S) of the theory has computed conservative thresholds of 4.56-6.70 kernel widths (Gaussian, K <= 3). The limit path itself is exact on 1154 targets down to 1.5 widths, and earlier trained-network checks were all above 4.9 widths, so they could not tell "the threshold is the truth" from "the threshold is conservative". The sweep was to train the actual network at K in {2,3,4}, separations 1.5-7 widths, amplitude ratios 1-100, noise 0-20% of the smallest amplitude, eps in {1e-4,1e-6,1e-8} (192 cells, 40 targets per cell held out) and compare the network's empirical transition with the path and with the thresholds. Present state: engine and development subset only (`results/synthetic/sepsweep_dev_raw.npz`, `sepsweep.md`); no held-out cell was run.

## 2. Teacher-student MLP test (`experiments/exp_mlp_dormancy.py`, plan `results/ml/mlp_plan.md`)
Purpose: test the two practical messages (readout depth >= 3 keeps unneeded units dormant; per-coordinate optimizers such as Adam and sign descent erase it) outside RBF/KAN, in one-hidden-layer tanh MLPs with K in {2,4,8} teacher neurons, d in {10,30}, width 128, depth-p product readouts p in {1,2,3}, six optimizers/arms, 7560 test runs planned. Present state: development runs only (`results/ml/mlp_dev_raw.npz`, `mlp_dev.md`, `mlp_freeze.json`); the test teachers were never generated.

Consequence: the claim that the mechanism extends to standard networks is untested. All evidence is on RBF and one-layer KAN edges.
