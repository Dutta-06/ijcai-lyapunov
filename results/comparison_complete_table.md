# Higher-dimensional comparison — complete (nothing excluded)

Single trajectory per cell; steady state = last 40% of 30 s. Step time = mean wall-clock of the method's update call (GP: per update, every 50–100 steps; NTK: per update, every 50 steps at n=1 and every step at n=2, 3, 6; proposed: per step). Every baseline regresses on noise-free samples of the true f(x, t); the proposed method uses only the tracking error. No valid significance test is available (one run per cell).

| benchmark | method | mean SS ‖e‖ | median SS | max SS | rank (mean SS) | mean update time (ms) |
|---|---|---|---|---|---|---|
| n=1 (1-D plant, Lorenz ref) | Proposed (Proj) | 0.2342 | 0.1496 | 1.290 | 3/5 | 0.10 |
| n=1 (1-D plant, Lorenz ref) | Exact GP | 0.2206 | 0.1415 | 1.016 | 1/5 | 273.78 |
| n=1 (1-D plant, Lorenz ref) | Sparse GP | 0.2285 | 0.1914 | 0.874 | 2/5 | 52.92 |
| n=1 (1-D plant, Lorenz ref) | Deep GP | 0.3963 | 0.2994 | 1.436 | 4/5 | 94.34 |
| n=1 (1-D plant, Lorenz ref) | NTK | 0.4286 | 0.3318 | 1.324 | 5/5 | 13.22 |
| n=2 (Van der Pol) | Proposed (Proj) | 0.0372 | 0.0137 | 0.310 | 2/5 | 0.22 |
| n=2 (Van der Pol) | Exact GP | 0.0555 | 0.0370 | 0.278 | 3/5 | 498.50 |
| n=2 (Van der Pol) | Sparse GP | 0.1771 | 0.0929 | 0.886 | 5/5 | 150.30 |
| n=2 (Van der Pol) | Deep GP | 0.1758 | 0.1169 | 1.097 | 4/5 | 415.37 |
| n=2 (Van der Pol) | NTK | 0.0015 | 0.0009 | 0.006 | 1/5 | 1.59 |
| n=3 (Lorenz plant) | Proposed (Proj) | 7.7190 | 5.7197 | 30.476 | 3/5 | 4.11 |
| n=3 (Lorenz plant) | Exact GP | 10.2774 | 8.8437 | 34.467 | 5/5 | 1888.80 |
| n=3 (Lorenz plant) | Sparse GP | 3.2729 | 2.5099 | 19.176 | 1/5 | 264.64 |
| n=3 (Lorenz plant) | Deep GP | 9.5344 | 8.3305 | 33.755 | 4/5 | 802.93 |
| n=3 (Lorenz plant) | NTK | 5.7642 | 4.2944 | 18.240 | 2/5 | 2.35 |
| n=6 (coupled Duffing) | Proposed (Proj) | 0.0044 | 0.0013 | 0.026 | 2/5 | 0.60 |
| n=6 (coupled Duffing) | Exact GP | 0.0841 | 0.0695 | 0.279 | 3/5 | 668.01 |
| n=6 (coupled Duffing) | Sparse GP | 0.1444 | 0.0872 | 0.670 | 5/5 | 523.83 |
| n=6 (coupled Duffing) | Deep GP | 0.0953 | 0.0553 | 0.342 | 4/5 | 1607.48 |
| n=6 (coupled Duffing) | NTK | 0.0013 | 0.0010 | 0.004 | 1/5 | 4.75 |
