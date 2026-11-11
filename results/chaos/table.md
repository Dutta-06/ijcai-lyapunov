# Theorem D′: synchronous coupling of the m-particle dithered loop to an M = 4096 reference (`exp_chaos.py`)

n = 1, k_e = 10, σ = 0.02, G = diag(4.0, 1.0, 1.0), T = 20 s, dt = 0.001, Euler–Maruyama for the particles and RK4 for the plant (as `langevin.simulate`). The m-system shares its initial particles and Brownian motions with the first m particles of the reference (M = 4096). sup_t D_m is averaged over 24 seeds (mean ± standard error of the mean); the median and the spread over seeds are shown because single seeds differ by up to an order of magnitude at fixed m. The shared-particle coupling makes the driving fluctuation proportional to (1/m − 1/M)^{1/2}, so the mean is also divided by that factor (a constant row means the predicted law).

## τ = 0

| m | mean sup_t D_m ± s.e. | median | min, max over seeds | mean / (1/m − 1/M)^{1/2} | mean · m^{1/2} |
|---|---|---|---|---|---|
| 16 | 1.487e-01 ± 1.7e-02 | 1.107e-01 | 2.34e-02, 3.56e-01 | 0.596 | 0.595 |
| 32 | 1.338e-01 ± 1.7e-02 | 1.148e-01 | 2.82e-02, 3.33e-01 | 0.760 | 0.757 |
| 64 | 9.699e-02 ± 1.2e-02 | 8.810e-02 | 1.19e-02, 2.36e-01 | 0.782 | 0.776 |
| 128 | 6.083e-02 ± 6.4e-03 | 5.996e-02 | 6.40e-03, 1.23e-01 | 0.699 | 0.688 |
| 256 | 4.341e-02 ± 6.4e-03 | 3.155e-02 | 1.39e-02, 1.41e-01 | 0.717 | 0.694 |
| 512 | 2.922e-02 ± 4.1e-03 | 2.477e-02 | 8.89e-03, 1.03e-01 | 0.707 | 0.661 |
| 1024 | 2.169e-02 ± 2.3e-03 | 1.936e-02 | 6.82e-03, 5.25e-02 | 0.801 | 0.694 |

Log-log slope of the mean against m over m ≤ 256: -0.47 ± 0.06 (bootstrap s.e. over seeds); of the median: -0.46 (Theorem D′: −0.5). Slope of the mean against (1/m − 1/M)^(1/2) over all m: 0.94 ± 0.08 (predicted 1).

## τ = 0.002

| m | mean sup_t D_m ± s.e. | median | min, max over seeds | mean / (1/m − 1/M)^{1/2} | mean · m^{1/2} |
|---|---|---|---|---|---|
| 16 | 1.524e-01 ± 1.8e-02 | 1.270e-01 | 2.13e-02, 3.63e-01 | 0.611 | 0.610 |
| 32 | 1.333e-01 ± 1.6e-02 | 1.142e-01 | 3.71e-02, 3.14e-01 | 0.757 | 0.754 |
| 64 | 9.818e-02 ± 1.1e-02 | 8.851e-02 | 3.19e-02, 2.45e-01 | 0.792 | 0.785 |
| 128 | 5.857e-02 ± 5.6e-03 | 5.565e-02 | 2.11e-02, 1.24e-01 | 0.673 | 0.663 |
| 256 | 4.602e-02 ± 6.0e-03 | 3.779e-02 | 1.58e-02, 1.30e-01 | 0.760 | 0.736 |
| 512 | 2.955e-02 ± 3.8e-03 | 2.530e-02 | 8.63e-03, 8.26e-02 | 0.715 | 0.669 |
| 1024 | 2.321e-02 ± 2.2e-03 | 1.880e-02 | 9.10e-03, 4.52e-02 | 0.857 | 0.743 |

Log-log slope of the mean against m over m ≤ 256: -0.46 ± 0.05 (bootstrap s.e. over seeds); of the median: -0.45 (Theorem D′: −0.5). Slope of the mean against (1/m − 1/M)^(1/2) over all m: 0.92 ± 0.08 (predicted 1).

Numbers are produced by `experiments/exp_chaos.py` (raw per-seed values in `raw_sup_D.npz`); nothing is edited by hand.
