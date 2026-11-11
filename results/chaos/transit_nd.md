# Proposition D‴.54: the transit in ℝⁿ (`exp_transit_nd.py`)

The free law with moving centres (m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width 0.15 in the plane on a square lattice of spacing 0.125 (H: a hexagonal one), the samples on a grid of spacing 0.05; masses ε; formed: |w| > 0.01|A|; κ = πρ(2c)^{1/2} = 0.7025.

## A. The k-fold near tie

One feature near an edge (the kernels at (0, 0), (0.125, 0)) or a vertex (and (0.125, 0), (0, 0.125), (0.125, 0.125)) of the lattice — 8 × 8 kernels around it, the samples symmetric about the vertex —, its margins ε^{1/2}x at the entry; the loop to 1.32 τ_entry. The gaps x (the runner-ups', as placed; at the vertex the diagonal's as it follows), the shares of the tie's kernels in the loop against the inner problem's at the gaps of the smallest ε (the feature at the tie; extrapolated in its start from xi0 = 800, 1600, the change between those in the last column), the softmax at those gaps and, at the vertex, the product of the two edges' inner problems posed in the vertex's geometry (the feature H/2 off each edge's axis):

| tie | A | x | shares: loop (ε = 1e-06, 1e-07, 1e-08, 1e-09, 1e-10) | inner problem | softmax | the edges' product | the others' largest \|w\| (ε = 1e-06, 1e-07, 1e-08, 1e-09, 1e-10) | the loop's successive changes | xi0's change |
|---|---|---|---|---|---|---|---|---|---|
| edge | 0.3 | 0.000 | 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5 | 0.5, 0.5 | 0.5, 0.5 | — | 2.50e-08, 7.63e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 1.4e-13, 1.6e-14, 9.5e-13, 5.2e-12 | 1.1e-12 |
| edge | 0.3 | 0.500 | 0.5886, 0.4114; 0.59, 0.41; 0.5905, 0.4095; 0.5906, 0.4094; 0.5906, 0.4094 | 0.5907, 0.4093 | 0.5869, 0.4131 | — | 2.51e-08, 7.64e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 1.4e-03, 4.5e-04, 1.4e-04, 4.6e-05 | 1.6e-06 |
| edge | 0.3 | 1.000 | 0.6718, 0.3282; 0.6744, 0.3256; 0.6752, 0.3248; 0.6754, 0.3246; 0.6755, 0.3245 | 0.6756, 0.3244 | 0.6687, 0.3313 | — | 2.52e-08, 7.64e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.5e-03, 8.1e-04, 2.6e-04, 8.2e-05 | 2.8e-06 |
| edge | 0.3 | 2.000 | 0.8074, 0.1926; 0.8109, 0.1891; 0.812, 0.188; 0.8123, 0.1877; 0.8125, 0.1875 | 0.8125, 0.1875 | 0.803, 0.197 | — | 2.53e-08, 7.65e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 3.5e-03, 1.1e-03, 3.5e-04, 1.1e-04 | 3.9e-06 |
| edge | 0.3 | 4.000 | 0.946, 0.05402; 0.9481, 0.0519; 0.9488, 0.05123; 0.949, 0.05102; 0.949, 0.05095 | 0.9491, 0.05092 | 0.9432, 0.05679 | — | 2.54e-08, 7.66e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.1e-03, 6.7e-04, 2.1e-04, 6.7e-05 | 2.5e-06 |
| edge | 0.3 | 8.000 | 0.9958, 0.004171; 0.9961, 0.003929; 0.9961, 0.003855; 0.9962, 0.003832; 0.9962, 0.003824 | 0.9962, 0.003821 | 0.9964, 0.003612 | — | 2.54e-08, 7.66e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.4e-04, 7.4e-05, 2.3e-05, 7.4e-06 | 3.0e-07 |
| edge | 1 | 0.000 | 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5 | 0.5, 0.5 | 0.5, 0.5 | — | 2.57e-08, 7.69e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 8.0e-15, 6.6e-14, 2.1e-12, 4.3e-13 | 9.3e-13 |
| edge | 1 | 0.500 | 0.5886, 0.4114; 0.59, 0.41; 0.5905, 0.4095; 0.5906, 0.4094; 0.5907, 0.4093 | 0.5907, 0.4093 | 0.5869, 0.4131 | — | 2.58e-08, 7.70e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 1.4e-03, 4.5e-04, 1.4e-04, 4.6e-05 | 1.6e-06 |
| edge | 1 | 1.000 | 0.6719, 0.3281; 0.6744, 0.3256; 0.6752, 0.3248; 0.6755, 0.3245; 0.6756, 0.3244 | 0.6756, 0.3244 | 0.6687, 0.3313 | — | 2.58e-08, 7.71e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 2.5e-03, 8.1e-04, 2.6e-04, 8.2e-05 | 2.8e-06 |
| edge | 1 | 2.000 | 0.8075, 0.1925; 0.811, 0.189; 0.8121, 0.1879; 0.8125, 0.1875; 0.8126, 0.1874 | 0.8126, 0.1874 | 0.803, 0.197 | — | 2.60e-08, 7.72e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 3.5e-03, 1.1e-03, 3.5e-04, 1.1e-04 | 3.9e-06 |
| edge | 1 | 4.000 | 0.9464, 0.05364; 0.9485, 0.0515; 0.9492, 0.05083; 0.9494, 0.05062; 0.9494, 0.05055 | 0.9495, 0.05052 | 0.9432, 0.05679 | — | 2.61e-08, 7.72e-10, 2.40e-11, 7.53e-13, 2.38e-14 | 2.1e-03, 6.7e-04, 2.1e-04, 6.7e-05 | 2.5e-06 |
| edge | 1 | 8.000 | 0.9966, 0.003442; 0.9968, 0.003211; 0.9969, 0.003141; 0.9969, 0.003118; 0.9969, 0.003111 | 0.9969, 0.003108 | 0.9964, 0.003612 | — | 2.61e-08, 7.73e-10, 2.40e-11, 7.53e-13, 2.38e-14 | 2.3e-04, 7.1e-05, 2.2e-05, 7.0e-06 | 2.8e-07 |
| vertex | 0.3 | 0.000, 0.000, 0.000 | 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 3.83e-09, 1.20e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.9e-13, 4.1e-13, 9.0e-13, 2.2e-12 | 7.8e-13 |
| vertex | 0.3 | 0.500, 0.500, 1.000 | 0.3509, 0.2415, 0.2415, 0.1662; 0.3526, 0.2412, 0.2412, 0.165; 0.3532, 0.2411, 0.2411, 0.1646; 0.3534, 0.2411, 0.2411, 0.1645; 0.3534, 0.2411, 0.2411, 0.1645 | 0.3535, 0.241, 0.241, 0.1644 | 0.3445, 0.2424, 0.2424, 0.1706 | 0.3535, 0.2411, 0.2411, 0.1644 | 3.86e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.8e-03, 5.7e-04, 1.8e-04, 5.8e-05 | 3.9e-06 |
| vertex | 0.3 | 1.000, 1.000, 2.000 | 0.4605, 0.2181, 0.2181, 0.1034; 0.4641, 0.2171, 0.2171, 0.1017; 0.4653, 0.2168, 0.2168, 0.1012; 0.4657, 0.2167, 0.2167, 0.101; 0.4658, 0.2166, 0.2166, 0.1009 | 0.4658, 0.2166, 0.2166, 0.1009 | 0.4472, 0.2215, 0.2215, 0.1097 | 0.466, 0.2166, 0.2166, 0.1007 | 3.88e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 3.7e-03, 1.2e-03, 3.7e-04, 1.2e-04 | 8.0e-06 |
| vertex | 0.3 | 1.000, 3.000, 4.000 | 0.613, 0.2905, 0.06527, 0.0313; 0.6175, 0.2889, 0.06347, 0.03012; 0.619, 0.2884, 0.06289, 0.02975; 0.6195, 0.2882, 0.06271, 0.02963; 0.6196, 0.2881, 0.06266, 0.0296 | 0.6197, 0.2881, 0.06263, 0.02958 | 0.5963, 0.2954, 0.07247, 0.0359 | 0.6202, 0.2883, 0.06247, 0.02904 | 3.95e-09, 1.22e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 4.6e-03, 1.5e-03, 4.6e-04, 1.5e-04 | 9.9e-06 |
| vertex | 0.3 | 2.000, 2.000, 4.000 | 0.6668, 0.1496, 0.1496, 0.03399; 0.6728, 0.1473, 0.1473, 0.03274; 0.6747, 0.1465, 0.1465, 0.03234; 0.6753, 0.1463, 0.1463, 0.03221; 0.6755, 0.1462, 0.1462, 0.03218 | 0.6756, 0.1461, 0.1461, 0.03216 | 0.6448, 0.1582, 0.1582, 0.03882 | 0.676, 0.1462, 0.1462, 0.03161 | 3.92e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 6.0e-03, 1.9e-03, 6.1e-04, 1.9e-04 | 1.3e-05 |
| vertex | 0.3 | 2.000, 5.000, 7.000 | 0.7962, 0.179, 0.01962, 0.005226; 0.801, 0.1754, 0.01868, 0.004926; 0.8026, 0.1742, 0.01838, 0.004834; 0.803, 0.1739, 0.01829, 0.004805; 0.8032, 0.1738, 0.01826, 0.004796 | 0.8033, 0.1737, 0.01825, 0.004791 | 0.7797, 0.1913, 0.02326, 0.005707 | 0.8042, 0.1739, 0.01807, 0.003907 | 4.02e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.8e-03, 1.5e-03, 4.9e-04, 1.5e-04 | 1.0e-05 |
| vertex | 0.3 | 4.000, 4.000, 8.000 | 0.9049, 0.04594, 0.04594, 0.003234; 0.9091, 0.04395, 0.04395, 0.003035; 0.9104, 0.04333, 0.04333, 0.002974; 0.9108, 0.04313, 0.04313, 0.002955; 0.9109, 0.04307, 0.04307, 0.002949 | 0.911, 0.04304, 0.04304, 0.002946 | 0.8896, 0.05357, 0.05357, 0.003225 | 0.9119, 0.04304, 0.04304, 0.002032 | 4.00e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.2e-03, 1.3e-03, 4.1e-04, 1.3e-04 | 8.7e-06 |
| vertex | 1 | 0.000, 0.000, 0.000 | 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 3.87e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.3e-13, 5.7e-13, 8.1e-13, 4.4e-12 | 9.6e-13 |
| vertex | 1 | 0.500, 0.500, 1.000 | 0.3509, 0.2415, 0.2415, 0.1661; 0.3527, 0.2412, 0.2412, 0.1649; 0.3533, 0.2411, 0.2411, 0.1645; 0.3535, 0.2411, 0.2411, 0.1644; 0.3535, 0.2411, 0.2411, 0.1644 | 0.3536, 0.241, 0.241, 0.1643 | 0.3445, 0.2424, 0.2424, 0.1706 | 0.3536, 0.241, 0.241, 0.1643 | 3.90e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 1.8e-03, 5.7e-04, 1.8e-04, 5.8e-05 | 3.9e-06 |
| vertex | 1 | 1.000, 1.000, 2.000 | 0.4607, 0.2181, 0.2181, 0.1032; 0.4643, 0.2171, 0.2171, 0.1015; 0.4655, 0.2168, 0.2168, 0.101; 0.4659, 0.2167, 0.2167, 0.1008; 0.466, 0.2166, 0.2166, 0.1007 | 0.4661, 0.2166, 0.2166, 0.1007 | 0.4472, 0.2215, 0.2215, 0.1097 | 0.4661, 0.2166, 0.2166, 0.1007 | 3.93e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 3.7e-03, 1.2e-03, 3.7e-04, 1.2e-04 | 8.0e-06 |
| vertex | 1 | 1.000, 3.000, 4.000 | 0.6136, 0.2907, 0.06498, 0.03078; 0.6182, 0.2891, 0.06317, 0.02959; 0.6196, 0.2886, 0.06259, 0.02921; 0.6201, 0.2884, 0.06241, 0.02909; 0.6202, 0.2883, 0.06235, 0.02905 | 0.6203, 0.2883, 0.06233, 0.02904 | 0.5963, 0.2954, 0.07247, 0.0359 | 0.6204, 0.2883, 0.06231, 0.02896 | 4.00e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.6e-03, 1.5e-03, 4.7e-04, 1.5e-04 | 9.9e-06 |
| vertex | 1 | 2.000, 2.000, 4.000 | 0.6674, 0.1496, 0.1496, 0.03349; 0.6734, 0.1472, 0.1472, 0.03222; 0.6753, 0.1464, 0.1464, 0.03182; 0.6759, 0.1462, 0.1462, 0.0317; 0.6761, 0.1461, 0.1461, 0.03166 | 0.6762, 0.1461, 0.1461, 0.03164 | 0.6448, 0.1582, 0.1582, 0.03882 | 0.6763, 0.1461, 0.1461, 0.03156 | 3.97e-09, 1.22e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 6.0e-03, 1.9e-03, 6.1e-04, 1.9e-04 | 1.3e-05 |
| vertex | 1 | 2.000, 5.000, 7.000 | 0.7975, 0.1791, 0.01897, 0.004466; 0.8023, 0.1755, 0.01802, 0.004174; 0.8038, 0.1744, 0.01772, 0.004085; 0.8043, 0.174, 0.01762, 0.004056; 0.8045, 0.1739, 0.01759, 0.004047 | 0.8046, 0.1738, 0.01758, 0.004043 | 0.7797, 0.1913, 0.02326, 0.005707 | 0.8048, 0.1739, 0.01755, 0.003791 | 4.07e-09, 1.23e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.8e-03, 1.5e-03, 4.9e-04, 1.5e-04 | 1.0e-05 |
| vertex | 1 | 4.000, 4.000, 8.000 | 0.9063, 0.04556, 0.04556, 0.002559; 0.9105, 0.04356, 0.04356, 0.002374; 0.9118, 0.04293, 0.04293, 0.002318; 0.9122, 0.04274, 0.04274, 0.0023; 0.9124, 0.04267, 0.04267, 0.002295 | 0.9124, 0.04264, 0.04264, 0.002292 | 0.8896, 0.05357, 0.05357, 0.003225 | 0.9127, 0.04265, 0.04265, 0.001993 | 4.04e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.2e-03, 1.3e-03, 4.1e-04, 1.3e-04 | 8.8e-06 |

The largest difference of the shares from the inner problem's (at the smallest ε): edge — the loop 5.31e-05; vertex — the loop 9.57e-05; vertex — the edges’ product 9.14e-04; edge — the softmax 9.66e-03; vertex — the softmax 3.14e-02. The loop's successive changes fall by 3.15 per decade of ε (median; √10 = 3.16), their last 1.9e-04 at most; the inner problem's change between its starts at most 1.3e-05. The edges' product against the 4-kernel problem, by amplitude (all shares; the edges' runner-ups only — the rest is the diagonal's frozen tail, which does not factorise, and the entrant's balance): A = 0.3: 9.1e-04 (2.0e-04); A = 1: 3.0e-04 (4.9e-05).

## T. The edge's split along the face

The 2-kernel inner problem (A = 0.8) with the feature on the face at the offset t from the edge's axis: the runner-up's share p(x; t) at x = 2, 4, 6, 8 and the crossing x_θ(t) of 0.01:

| t/H | p(2) | p(4) | p(6) | p(8) | x_θ(t) |
|---|---|---|---|---|---|
| 0 | 0.18740 | 0.05056 | 0.01234 | 0.00319 | 6.300 |
| 0.125 | 0.18679 | 0.05019 | 0.01220 | 0.00315 | 6.284 |
| 0.25 | 0.18498 | 0.04906 | 0.01179 | 0.00302 | 6.233 |
| 0.375 | 0.18195 | 0.04721 | 0.01112 | 0.00282 | 6.149 |
| 0.5 | 0.17771 | 0.04470 | 0.01024 | 0.00256 | 6.034 |

x_θ averaged over the face (t uniform in [−H/2, H/2]): 6.208 (at t = 0: 6.300; at the vertex, t = H/2: 6.034). For a feature uniform in a Voronoi cell the offsets from a face are (2σ²/H)ε^{1/2}x at leading order, so P(N ≥ 2) = 4·(2σ²/H²)·x̄_θ·ε^{1/2} = 71.5·ε^{1/2} (the four faces' layers, less O(ε) where they overlap).

Beyond H/2 (T2, by bisection): x_θ = 5.888 (t = 0.625H), 5.715 (t = 0.75H), 5.517 (t = 0.875H), 5.298 (t = 1H). In the continuum x_θ depends on |t| only, so on hℤⁿ (σ/H as here) a face is the (n − 1)-cube of offsets and E[N − 1] = 2n(2σ²/H²)x̄_θ⁽ⁿ⁾ε^{1/2}: n = 2: x̄_θ = 6.208, 71.5·ε^{1/2}; n = 3: x̄_θ = 6.120, 105.8·ε^{1/2}; n = 4: x̄_θ = 6.033, 139.0·ε^{1/2}; n = 5: x̄_θ = 5.947, 171.3·ε^{1/2} (the n = 2 average by the midpoint rule; the trapezoid rule above).

## B. The formation map at a vertex

A = 0.8; the gaps (x₁₀, x₀₁) of the edge runner-ups on a grid, x₁₀ ≤ x₀₁, the diagonal's x₁₁ = x₁₀ + x₀₁ − ε^{1/2}x₁₀x₀₁ exactly (the tensor grid's correlations multiply); the number of kernels formed by the loop at ε = 1e-06, 1e-08, 1e-10, by the 4-kernel inner problem (B4's), by the independent rule #{j: x_j < x_θ} + 1 with x_θ = 6.034 (the edge's inner problem at the vertex's offset t = H/2 (T)), and by the product rule with the edges' factors at the vertex's offset t = H/2 (T):

| x₁₀ | x₀₁ | x₁₁ | loop | inner problem | the independent rule | the product rule (T) |
|---|---|---|---|---|---|---|
| 0 | 0 | 0.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 1 | 1.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 2 | 2.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 3 | 3.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 4 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 5 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 6 | 6.00 | 2, 2, 2 | 2 | 4 | 2 |
| 0 | 7 | 7.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 8 | 8.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 9 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 10 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 1 | 2.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 2 | 3.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 3 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 4 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 5 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 1 | 6 | 7.00 | 2, 2, 2 | 2 | 3 | 2 |
| 1 | 7 | 8.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 8 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 9 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 10 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 2 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 2 | 3 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 2 | 4 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 2 | 5 | 7.00 | 3, 3, 3 | 3 | 3 | 3 |
| 2 | 6 | 8.00 | 2, 2, 2 | 2 | 3 | 2 |
| 2 | 7 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 8 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 9 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 10 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 3 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 3 | 4 | 7.00 | 3, 3, 3 | 3 | 3 | 3 |
| 3 | 5 | 8.00 | 3, 3, 3 | 3 | 3 | 3 |
| 3 | 6 | 9.00 | 3, 2, 2 | 2 | 3 | 2 |
| 3 | 7 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 8 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 9 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 10 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 4 | 8.00 | 3, 3, 3 | 3 | 3 | 3 |
| 4 | 5 | 9.00 | 3, 3, 3 | 3 | 3 | 3 |
| 4 | 6 | 10.00 | 3, 2, 2 | 2 | 3 | 2 |
| 4 | 7 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 8 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 9 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 10 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 5 | 10.00 | 3, 3, 3 | 3 | 3 | 3 |
| 5 | 6 | 11.00 | 3, 3, 3 | 3 | 3 | 3 |
| 5 | 7 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 8 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 9 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 10 | 15.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 6 | 12.00 | 3, 3, 3 | 3 | 3 | 3 |
| 6 | 7 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 8 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 9 | 15.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 10 | 16.00 | 2, 2, 2 | 2 | 2 | 2 |
| 7 | 7 | 14.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 8 | 15.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 9 | 16.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 10 | 17.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 8 | 16.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 9 | 17.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 10 | 18.00 | 1, 1, 1 | 1 | 1 | 1 |
| 9 | 9 | 18.00 | 1, 1, 1 | 1 | 1 | 1 |
| 9 | 10 | 19.00 | 1, 1, 1 | 1 | 1 | 1 |
| 10 | 10 | 20.00 | 1, 1, 1 | 1 | 1 | 1 |

The loop's count equals the inner problem's at 64 of 66 (ε = 1e-06), 66 of 66 (ε = 1e-08), 66 of 66 (ε = 1e-10); the independent rule at 58 of 66, the product rule at 66 of 66.

## B4. The areas at a vertex

The 4-kernel inner problem on the grid of gaps 0 … 10 in steps of 0.25 (861 problems, x₁₀ ≤ x₀₁): the region N ≥ 3 has the area 33.5 and N = 4 15.5 in the gaps' quarter-plane (trapezoid weights; 648 of 861 problems settled to the shares' tolerance 1e-07; the 213 others rerun with the stop's reason recorded: 209 at the doubling cap (the last drift 1.0e-07 … 3.7e-07), 1 at the clock (the last drift 2.0e-07), 3 settled on rerun, their counts unchanged at 213 of 213, their shares moved by at most 1.1e-06; continued to 24 doublings (B4x): the counts unchanged at 213 of 213, the shares moved after the first cap by at most 1.6e-06, the drifts peaking at doublings 1–13, 211 settled, 0 at the new cap, 2 at the clock); the product rule with the factors at t = H/2 on the same grid: 33.5 and 15.1 (on a grid of step 0.05: 33.5 and 15.3, the coarse grid's error), its count equal to the inner problem's at 1675 of 1681 grid points. Each vertex's sector in a cell maps to the gaps by (2σ²ε^{1/2}/H)·(x₁₀, x₀₁) and a cell has four: P(N ≥ 3) = 4(2σ²/H²)²·33.5·ε = 1112·ε and P(N = 4) = 513·ε.

## C. Several features in the plane

Random three-feature targets (seed 54; centres in [−0.7, 0.7]² at least 3σ apart, amplitudes ±U(0.4, 1)) on the lattice over [−1, 1]² (17² kernels; the samples' grid 41²); the sliding path against the loop to 1.3 times the path's last entry; the never-entering kernels' displacement vectors against the creep's εC_j integrated along the path (`creep_constants_nd`; where |εC_j| is above 10⁻³ of the largest: the largest relative error and the smallest cosine):

| target | the path (entries) | ε | formed = the path's | the never-entering kernels' largest \|w\|·ε^{−3/2} | their largest displacement/ε | creep: relative error, smallest cosine | the block: max \|w − w_path\|, max \|c − c_path\| | rms residual |
|---|---|---|---|---|---|---|---|---|
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-06 | yes | 2862 | 4838 | 4.57e-01, 0.979100 | 3.01e-06, 2.14e-07 | 4.61e-08 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-07 | yes | 2039 | 3849 | 1.21e-01, 0.998353 | 6.94e-08, 4.99e-09 | 1.11e-09 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-08 | yes | 1851 | 3594 | 3.64e-02, 0.999843 | 2.01e-09, 1.43e-10 | 3.23e-11 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-09 | yes | 1797 | 3518 | 1.13e-02, 0.999985 | 6.18e-11, 4.40e-12 | 9.94e-13 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-10 | yes | 1780 | 3495 | 3.56e-03, 0.999998 | 1.94e-12, 1.38e-13 | 3.12e-14 |
| 0 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.508 | 1.006 | | 1.508, 1.509 | |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-06 | yes | 1.369e+05 | 4.749e+04 | 6.02e-01, 1.000000 | 1.36e-04, 3.63e-06 | 1.27e-07 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-07 | yes | 4.676e+04 | 3.729e+04 | 2.58e-01, 1.000000 | 1.41e-06, 1.16e-07 | 1.23e-08 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-08 | yes | 3.576e+04 | 3.193e+04 | 7.71e-02, 1.000000 | 3.41e-08, 2.94e-09 | 3.29e-10 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-09 | yes | 3.305e+04 | 3.035e+04 | 2.38e-02, 1.000000 | 9.95e-10, 8.63e-11 | 9.70e-12 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-10 | yes | 3.225e+04 | 2.987e+04 | 7.48e-03, 1.000000 | 3.07e-11, 2.66e-12 | 3.00e-13 |
| 1 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.522 | 1.014 | | 1.522, 1.522 | |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-06 | yes | 6041 | 5.122e+04 | 1.67e+00, 0.999876 | 4.89e-06, 9.40e-07 | 2.04e-07 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-07 | yes | 3088 | 2.318e+04 | 2.06e-01, 0.999988 | 9.18e-08, 1.17e-08 | 1.76e-09 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-08 | yes | 2736 | 2.027e+04 | 5.49e-02, 0.999999 | 2.60e-09, 3.22e-10 | 4.81e-11 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-09 | yes | 2642 | 1.953e+04 | 1.65e-02, 1.000000 | 7.94e-11, 9.82e-12 | 1.47e-12 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-10 | yes | 2613 | 1.932e+04 | 5.14e-03, 1.000000 | 2.48e-12, 3.07e-13 | 4.60e-14 |
| 2 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.510 | 1.010 | | 1.510, 1.510 | |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-06 | yes | 6.347e+04 | 3.56e+04 | 5.70e-01, 0.999984 | 6.27e-05, 3.88e-06 | 1.65e-07 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-07 | yes | 2.763e+04 | 2.731e+04 | 2.04e-01, 0.999999 | 8.32e-07, 9.93e-08 | 8.09e-09 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-08 | yes | 2.232e+04 | 2.408e+04 | 6.18e-02, 1.000000 | 2.12e-08, 2.62e-09 | 2.22e-10 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-09 | yes | 2.095e+04 | 2.312e+04 | 1.92e-02, 1.000000 | 6.29e-10, 7.79e-11 | 6.63e-12 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-10 | yes | 2.054e+04 | 2.282e+04 | 6.04e-03, 1.000000 | 1.95e-11, 2.42e-12 | 2.06e-13 |
| 3 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.518 | 1.012 | | 1.518, 1.517 | |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-06 | yes | 3.851e+06 | 8.48e+04 | 4.98e-01, 0.999098 | 3.87e-03, 5.69e-06 | 1.91e-07 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-07 | yes | 6.017e+05 | 2.339e+05 | 3.87e-01, 0.999963 | 1.87e-05, 1.51e-06 | 1.49e-07 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-08 | yes | 3.632e+05 | 1.894e+05 | 1.24e-01, 0.999997 | 3.49e-07, 3.69e-08 | 5.02e-09 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-09 | yes | 3.187e+05 | 1.75e+05 | 3.80e-02, 1.000000 | 9.68e-09, 1.04e-09 | 1.45e-10 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-10 | yes | 3.066e+05 | 1.706e+05 | 1.19e-02, 1.000000 | 2.95e-10, 3.17e-11 | 4.43e-12 |
| 4 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.537 | 1.023 | | 1.537, 1.533 | |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-06 | yes | 8.397e+04 | 4.497e+04 | 5.28e-01, 0.999917 | 8.38e-05, 3.77e-06 | 2.17e-07 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-07 | yes | 3.498e+04 | 3.532e+04 | 2.00e-01, 0.999992 | 1.06e-06, 1.03e-07 | 1.28e-08 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-08 | yes | 2.82e+04 | 3.121e+04 | 6.07e-02, 0.999999 | 2.70e-08, 2.73e-09 | 3.54e-10 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-09 | yes | 2.648e+04 | 2.998e+04 | 1.89e-02, 1.000000 | 8.03e-10, 8.13e-11 | 1.06e-11 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-10 | yes | 2.597e+04 | 2.96e+04 | 5.93e-03, 1.000000 | 2.49e-11, 2.52e-12 | 3.29e-13 |
| 5 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.518 | 1.012 | | 1.518, 1.517 | |

## D. The probabilities over uniformly placed features

A = 0.8 at 600 uniform positions in the Voronoi cell of the kernel at (0, 0) (seed 541; the same positions at both masses, so the two rows are not independent); the loop's count N against the product rule at each position's gaps to its two nearest faces and its offsets along them (T), and the fractions against the rule's at the sampled positions and the asymptotics (P(N ≥ 2) = 4w/H − 4w²/H², w = 2σ²x̄_θε^{1/2}/H with the face average x̄_θ = 6.208; P(N ≥ 3) from B4's area):

| ε | N = 1, 2, 3, 4 (loop) | the product rule agrees | P(N ≥ 2): loop (± s.e.) | the rule at the positions | asymptotic | P(N ≥ 3): loop | the rule at the positions | asymptotic |
|---|---|---|---|---|---|---|---|---|
| 1e-05 | 458, 133, 2, 7 | 592 of 600 | 0.2367 ± 0.0174 | 0.2267 | 0.2134 | 0.0150 | 0.0133 | 0.0111 |
| 1e-06 | 548, 49, 1, 2 | 597 of 600 | 0.0867 ± 0.0115 | 0.0833 | 0.0702 | 0.0050 | 0.0033 | 0.0011 |

The misses — ε: the loop's count/the rule's at the gaps (x₁₀, x₀₁) [x_θ at their offsets]: 1e-05: 2/1 at (26.49, 6.28) [6.09, 6.23]; 1e-05: 2/1 at (6.52, 13.35) [6.15, 6.09]; 1e-05: 4/3 at (0.59, 4.85) [6.08, 6.04]; 1e-05: 2/1 at (6.36, 8.41) [6.11, 6.09]; 1e-05: 2/1 at (42.90, 6.49) [6.09, 6.29]; 1e-05: 3/2 at (6.12, 3.45) [6.06, 6.09]; 1e-05: 2/1 at (6.33, 7.32) [6.10, 6.09]; 1e-05: 2/1 at (6.30, 46.00) [6.29, 6.09]; 1e-06: 2/1 at (6.11, 8.25) [6.06, 6.05]; 1e-06: 3/2 at (3.92, 6.03) [6.05, 6.04]; 1e-06: 2/1 at (6.29, 76.61) [6.22, 6.05].

The offset-blind rule (every factor at t = 0, x_θ(0)) agrees with the loop at 596 of 600 (ε = 1e-05), 600 of 600 (ε = 1e-06): at these masses the loop's own threshold is shifted by more than x_θ varies along a face (D2), so D measures the shift, not x_θ(t) (referee round 49).


## H. A hexagonal lattice

Kernels on a hexagonal lattice of spacing 0.125, the feature (A = 0.8) at a Voronoi vertex — three cells meet there, the kernels (0, 0), (H, 0), (H/2, √3H/2) at H/√3 —, placed where the three correlate equally with it on the samples. The 3-kernel inner problem at the exact tie: shares 0.33333, 0.33333, 0.33333; on the gaps' grid 0 … 12 (step 0.5) the region where all three form, N = 3, has the area 36.6 (the square vertex's N ≥ 3: B4). A sector of a hexagonal vertex maps to the gaps with the Jacobian (2σ²ε^{1/2})²/((√3/2)H²) and a cell has six: P(N = 3) = 32σ⁴/H⁴·36.6·ε = 2426·ε. The loop at the exact tie:

| ε | shares | formed | the others' largest \|w\| | the gaps of the runner-ups and of the fourth kernel (units of ε^{1/2}) |
|---|---|---|---|---|
| 1e-08 | 0.33333, 0.33333, 0.33333 | 3 | 2.39e-11 | 5.551e-12, 5.551e-12, 1594 |
| 1e-10 | 0.33333, 0.33333, 0.33333 | 3 | 2.38e-14 | 5.551e-11, 5.551e-11, 1.594e+04 |

The fourth kernel's margin at the vertex is ε-free, 1 − e^{−H²/(4σ²)} in the continuum, so its gap grows like ε^{−1/2}: a fourth kernel never forms at a simple vertex (k ≤ n + 1).

## C2. Designed targets in the plane

Part C's lattice; 'near tie': feature 0 (A = 0.8) near the face between (0, 0.25) and (0.125, 0.25), on their axis, its runner-up's margin 0.0003 when it enters (first), two features far from it (−0.6 at (−0.5, −0.5), 0.7 at (0.5, −0.55)); 'window edge': feature 0 (0.8) at (0.96, 0) with −0.7 at (−0.4, 0.5) and 0.9 at (−0.3, −0.5). The loop to 1.3 times the path's last entry; feature 0's entrant and runner-up (the target's correlations at the start) and its nearest kernel; the runner-up's share in the loop (of the kernels started nearest feature 0) against the edge's 2-kernel inner problem p(x; t) at its gap (T), and its coefficient at the first output after the entry against the fold c(d)ℓ(d) of the same-sign features after it ((v); V's c(d) at A′/A = 0.875, interpolated in d and held at its last value beyond 7σ — there a lower estimate: c(d) grows with d), the largest; the kernels formed for feature 0 (|w| > 1% of its amplitude; by their final positions) against the rule 1 + [x < x_θ(t)], and the formed kernels that ended at another feature than the one they started nearest (captured):

| target | ε | entrant, runner-up, nearest | margin μ, gap x = μ/ε^{1/2} | runner-up's share: loop, p(x; t) | the runner-up after the entry: \|w\| (τ), the fold (d) | formed for feature 0: loop, rule | captured (kernel: start's feature → end's) | formed = the path's | max\|w − w_path\| | rms residual |
|---|---|---|---|---|---|---|---|---|---|---|
| near tie | 1e-06 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 0.3 | 0.4464, 0.4452 | 0.357 (79.3), 0.017 (6.08σ) | 2, 2 | none | no: [72, 146, 163, 208] | 3.57e-01 | 9.87e-08 |
| near tie | 1e-07 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 0.949 | 0.3338, 0.3328 | 0.267 (79.3), 0.017 (6.08σ) | 2, 2 | none | no: [72, 146, 163, 208] | 2.67e-01 | 3.47e-08 |
| near tie | 1e-08 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 3 | 0.1002, 0.09973 | 0.0802 (79.3), 0.017 (6.08σ) | 2, 2 | none | no: [72, 146, 163, 208] | 8.01e-02 | 4.02e-08 |
| near tie | 1e-09 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 9.49 | — (captured), 0.001322 | 0.00106 (79.3), 0.017 (6.08σ) | 1, 1 | (0.125, 0.250): 0 → 2 | no: [72, 146, 163] | 7.00e-01 | 1.95e-13 |
| near tie | 1e-10 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 30 | 8.205e-06, — | 6.48e-06 (79.3), 0.017 (6.08σ) | 1, 1 | none | yes | 6.33e-06 | 4.28e-08 |
| window edge | 1e-06 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 159 | 3.295e-08, — | 1.89e-09 (71.3), 9.1e-07 (9.04σ) | 1, 1 | none | yes | 4.15e-06 | 3.94e-08 |
| window edge | 1e-07 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 504 | 9.693e-10, — | 5.98e-11 (71.3), 9.1e-07 (9.04σ) | 1, 1 | none | yes | 8.72e-08 | 1.03e-09 |
| window edge | 1e-08 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 1.59e+03 | 2.998e-11, — | 1.89e-12 (71.3), 9.1e-07 (9.04σ) | 1, 1 | none | yes | 2.46e-09 | 2.95e-11 |
| window edge | 1e-09 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 5.04e+03 | 9.416e-13, — | 5.97e-14 (71.3), 9.1e-07 (9.04σ) | 1, 1 | none | yes | 7.53e-11 | 9.04e-13 |
| window edge | 1e-10 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 1.59e+04 | 2.971e-14, — | 1.89e-15 (71.3), 9.1e-07 (9.04σ) | 1, 1 | none | yes | 2.35e-12 | 2.83e-14 |

## V. A runner-up held at its feature: the static law of (v)

The feature (A = 0.8) at the origin and an unfitted one (A′ = 0.7) at the distance d (the continuum: a trapezoid of step σ/20 in one dimension, σ/8 in the plane along the diagonal); the entrant (w_e, c_e) and the runner-up's centre free, its coefficient w fixed: the local minimum's displacement Δ of the runner-up toward the unfitted feature (by continuation in w) against (v)'s law Δ² = σ²ℓ(d)·w_e/(w(w_e + w)), ℓ(d) = A′(d/σ)²e^{−d²/(4σ²)} (the pull projected off the entrant's φ, as the binding is: referee round 49); the fold — below it no minimum within 2σ — by bisection, the fold over ℓ(d), c(d), the minimum's displacement at the fold, and the law's there:

| d/σ | Δ/σ (the law) at w = 0.3, 0.1, 0.03, 0.01, 0.003, 0.001 | the fold w_f | ℓ(d) | c(d) = w_f/ℓ | Δ/σ at w_f: the minimum's (the law's) |
|---|---|---|---|---|---|
| 5 | 0.3164 (0.2665), —, —, —, —, — | 0.1424 | 0.0338 | 4.21 | 0.938 (0.442) |
| 5.5 | 0.1641 (0.1516), 0.4503 (0.3228), —, —, —, — | 0.05956 | 0.011 | 5.41 | 0.909 (0.413) |
| 6 | 0.0837 (0.0805), 0.1983 (0.1716), 0.4823 (0.3213), —, —, — | 0.02066 | 0.00311 | 6.64 | 0.860 (0.383) |
| 6.5 | 0.0407 (0.0399), 0.0912 (0.0851), 0.1862 (0.1593), 0.3708 (0.2659), —, — | 0.006018 | 0.000765 | 7.87 | 0.806 (0.355) |
| 7 | 0.0186 (0.0185), 0.0407 (0.0394), 0.0791 (0.0738), 0.1397 (0.1232), 0.3325 (0.2412), — | 0.001495 | 0.000164 | 9.11 | 0.759 (0.331) |

Along V's and V3's distinct one-dimensional branches (the runner-up toward A′; V3's A′/A = 0.875 column repeats V's) the rests' energy E*(w) (the minimum over the other parameters) increases with w at every step (223 steps; the secant dE*/dw from 1.58e-11 to 2.10e-03): with w free the runner-up's coefficient decays toward the fold on this side — the hold of (v) (referee round 49; the far side: VF).

The hold's third-moment law (`e3_law`: E*₃ = ∓(1/6)P₃(M₂*)^{3/2}(1 − w/w_e)(w(1 + w/w_e))^{−1/2}, P₃ = A′√π(d/σ)³e^{−d²/(4σ²)}/(8σ²)): the measured secant over the law's at w ≈ 0.28, 0.1, 0.037 (the displacement there, σ) — toward A′ (V) and away (VF):

| d/σ | side | secant / the law's (displacement) |
|---|---|---|
| 5 | toward | 2.049 (+0.379) |
| 5.5 | toward | 1.410 (+0.191), 3.866 (+0.450) |
| 6 | toward | 1.188 (+0.096), 1.854 (+0.198), 4.328 (+0.414) |
| 6.5 | toward | 1.090 (+0.046), 1.350 (+0.091), 1.862 (+0.168) |
| 7 | toward | 1.042 (+0.021), 1.153 (+0.041), 1.327 (+0.072) |
| 5 | away | 0.791 (-0.276), 0.392 (-0.452), 0.203 (-0.674) |
| 5.5 | away | 0.831 (-0.162), 0.512 (-0.276), 0.313 (-0.428) |
| 6 | away | 0.884 (-0.088), 0.650 (-0.155), 0.464 (-0.249) |
| 6.5 | away | 0.930 (-0.044), 0.779 (-0.080), 0.634 (-0.132) |
| 7 | away | 0.963 (-0.021), 0.878 (-0.038), 0.784 (-0.064) |

The two sides' mean at w ≈ 0.28 (an even next order enters the two ratios with opposite signs): d = 5σ: 1.41973, d = 5.5σ: 1.12067, d = 6σ: 1.03614, d = 6.5σ: 1.00994, d = 7σ: 1.00242.

VF. The far side (the tail started away from A′): the rest's displacement along e (σ; negative: away) at w = 0.3, 0.03, 0.003, 0.001, the fold, and the secant dE*/dw along the branch:

| n | d/σ | displacement at w = 0.3, 0.03, 0.003, 0.001 | fold | dE*/dw: steps < 0 of all, range |
|---|---|---|---|---|
| 1 | 5 | -0.2492, -0.7167, -1.3419, -1.4935 | none down to 0.0001 | 48 of 48, -5.75e-04 … -2.99e-05 |
| 1 | 5.5 | -0.1450, -0.4584, -1.0297, -1.2905 | none down to 0.0001 | 48 of 48, -1.65e-04 … -2.09e-06 |
| 1 | 6 | -0.0783, -0.2680, -0.6732, -0.9328 | none down to 0.0001 | 48 of 48, -3.77e-05 … -1.03e-07 |
| 1 | 6.5 | -0.0393, -0.1430, -0.3907, -0.5715 | none down to 0.0001 | 48 of 48, -5.49e-06 … -3.53e-09 |
| 1 | 7 | -0.0184, -0.0697, -0.2040, -0.3104 | none down to 0.0001 | 48 of 48, -4.37e-07 … -8.38e-11 |
| 2 | 6 | -0.0783, —, -0.6732, — | none down to 0.0001 | 12 of 12, -9.41e-06 … -3.12e-08 |

V3. The fold's constant against the unfitted feature's amplitude: c = w_f/ℓ(d) (one dimension, A = 0.8) at A′/A = 0.437, 0.875, 1.75, 3.5 (where the tail is unbound at w = 0.3 the continuation starts at 0.99A; '>' the tail unbound from there on, the fold above it; '<' bound down to w = 0.0001; 'entrant to A′' where the entrant leaves its feature for A′ — off the branch):

| d/σ | A′/A = 0.437 | A′/A = 0.875 | A′/A = 1.75 | A′/A = 3.5 |
|---|---|---|---|---|
| 5 | 4.502 | 4.215 | 3.761 | 3.135 (fresh start) |
| 6 | 6.731 | 6.644 | 6.480 | 6.182 |
| 7 | 9.116 | 9.106 | 9.086 | 9.047 |

The plane, d = 6σ along the diagonal: w = 0.3: along 0.0837, across 2.3e-08 (the law 0.0805); w = 0.079: along 0.2224, across 2.2e-07 (the law 0.1884); w = 0.0208: along 0.8015, across 7.2e-08 (the law 0.3816); the fold between w = 0.02066 and 0.02066.

## C3. The capture window

C2's near tie (feature 2, 0.7, at 6.08σ from feature 0) and its variant with feature 2 at 5.00σ (−60°), at 15 masses, the loop at the tolerance 1e-12 (captures sit on a saddle): the runner-up's gap x, its coefficient after the entry against the fold at feature 2's distance (V: c(d)ℓ(d), interpolated in d), the kernels formed per feature by their final positions, the captured kernels, and the outcome — bound (feature 0 two kernels, the runner-up formed there), captured (feature 0 one kernel, feature 2 fitted by the migrant), a capture tie (feature 2 two kernels), stranded at the horizon (one kernel per feature, the fit not exact at the end: the tail away from both features then — it can still be recalled to its own), one kernel (the runner-up not formed, the fit exact). The runner-up's coefficient is read at the first output after the entry (on the 5σ variant the runner-up has left by then: the coefficient is the captured kernel's). At a fixed μ the gap x and the mass vary together; C4 separates them:

| target | ε | gap x | the runner-up after the entry | the fold c(d)ℓ(d) | kernels per feature (final) | captured | rms | outcome |
|---|---|---|---|---|---|---|---|---|
| near tie | 1e-08 | 3 | 0.0801 | 0.01668 (FC) | 2, 1, 1 | none | 4.0e-08 | bound |
| near tie | 7.05e-09 | 3.57 | 0.0545 | 0.01668 (FC) | 2, 1, 1 | none | 5.1e-08 | bound |
| near tie | 4.97e-09 | 4.26 | 0.034 | 0.01668 (FC) | 2, 1, 1 | none | 7.1e-08 | bound |
| near tie | 3.5e-09 | 5.07 | 0.0191 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 5.1e-10 | captured |
| near tie | 3e-09 | 5.48 | 0.0143 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 8.6e-13 | captured |
| near tie | 2.47e-09 | 6.04 | 0.00966 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 7.3e-13 | captured |
| near tie | 1.74e-09 | 7.19 | 0.00435 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.2e-13 | captured |
| near tie | 1.22e-09 | 8.59 | 0.00178 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.5e-13 | captured |
| near tie | 8.63e-10 | 10.2 | 0.000722 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.6e-13 | captured |
| near tie | 6.08e-10 | 12.2 | 0.000296 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.3e-13 | captured |
| near tie | 4.29e-10 | 14.5 | 0.000128 | 0.01668 (FC) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.4e-13 | captured |
| near tie | 3.13e-10 | 17 | 6.26e-05 | 0.01668 (FC) | 1, 1, 2 | (0.125, 0.250): 0 → 2 | 1.4e-11 | capture tie |
| near tie | 3.02e-10 | 17.3 | 5.79e-05 | 0.01668 (FC) | 1, 1, 1 | none | 2.9e-06 | unsettled |
| near tie | 2.13e-10 | 20.6 | 2.78e-05 | 0.01668 (FC) | 1, 1, 1 | none | 2.8e-07 | unsettled |
| near tie | 1.5e-10 | 24.5 | 1.39e-05 | 0.01668 (FC) | 1, 1, 1 | none | 9.8e-08 | unsettled |
| near tie 5σ | 1e-08 | 3 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.9e-12 | captured |
| near tie 5σ | 7.05e-09 | 3.57 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.9e-12 | captured |
| near tie 5σ | 4.97e-09 | 4.26 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.7e-12 | captured |
| near tie 5σ | 3.5e-09 | 5.07 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.0e-12 | captured |
| near tie 5σ | 3e-09 | 5.48 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 8.0e-13 | captured |
| near tie 5σ | 2.47e-09 | 6.04 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 6.0e-13 | captured |
| near tie 5σ | 1.74e-09 | 7.19 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 3.5e-13 | captured |
| near tie 5σ | 1.22e-09 | 8.59 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.1e-13 | captured |
| near tie 5σ | 8.63e-10 | 10.2 | 0.339 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.2e-13 | captured |
| near tie 5σ | 6.08e-10 | 12.2 | 0.000303 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 7.3e-14 | captured |
| near tie 5σ | 4.29e-10 | 14.5 | 0.00013 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.3e-14 | captured |
| near tie 5σ | 3.13e-10 | 17 | 6.35e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.7e-14 | captured |
| near tie 5σ | 3.02e-10 | 17.3 | 5.87e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.6e-14 | captured |
| near tie 5σ | 2.13e-10 | 20.6 | 2.82e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.5e-14 | captured |
| near tie 5σ | 1.5e-10 | 24.5 | 1.41e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 9.3e-15 | captured |

## C4. The race: the gap held, the mass varied

C2's near tie with the runner-up's gap x held fixed (its margin μ = xε^{1/2}) at the masses 1e-08, 1e-09, 1e-10, 1e-11, 1e-12 (rtol 1e-12, the absolute tolerance min(10⁻¹⁴, 10⁻⁵ε): `c4_atol`): the runner-up's largest coefficient while within 2σ of feature 0 (after its entry) against the fold c(d)ℓ(d) at feature 2's distance; its departure (2σ from feature 0) in s = (τ − τ_e)/ε^{1/2} against feature 2's own entry s₂ = (τ₂ − τ_e)/ε^{1/2}, its coefficient at the departure; the final kernels per feature and the outcome; at the tolerance 1e-13 (C4b) the outcome, the largest coefficient held and the departure:

| x | ε | side (closest approach, σ) | max w held (/ fold) | departure s | s₂ | w at the departure (/ fold) | kernels per feature (final) | outcome | rtol 1e-13: side, outcome, max w held, departure s |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 1e-08 | toward (3.4e-15) | 0.0802 (4.81) | none | 9.71e+04 | — | 2, 1, 1 | bound | toward, bound, 0.0802, none |
| 3 | 1e-09 | toward (3.0e-15) | 0.0799 (4.79) | none | 3.07e+05 | — | 2, 1, 1 | bound | toward, bound, 0.0799, none |
| 3 | 1e-10 | toward (3.5e-15) | 0.0799 (4.79) | none | 9.7e+05 | — | 2, 1, 1 | bound | toward, bound, 0.0799, none |
| 3 | 1e-11 | toward (3.5e-15) | 0.0798 (4.79) | none | 3.07e+06 | — | 2, 1, 1 | bound | toward, bound, 0.0798, none |
| 3 | 1e-12 | toward (3.5e-15) | 0.0798 (4.79) | none | 9.7e+06 | — | 2, 1, 1 | bound | toward, bound, 0.0798, none |
| 4.9 | 1e-08 | toward (1.5e-14) | 0.0217 (1.30) | none | 9.72e+04 | — | 2, 1, 1 | bound | toward, bound, 0.0217, none |
| 4.9 | 1e-09 | toward (7.4e-15) | 0.0216 (1.29) | none | 3.07e+05 | — | 2, 1, 1 | bound | toward, bound, 0.0216, none |
| 4.9 | 1e-10 | toward (1.2e-14) | 0.0215 (1.29) | 3.52e+05 | 9.7e+05 | 0.0163 (0.977) | 1, 1, 1 | captured | toward, captured, 0.0215, 3.52e+05 (w 0.0163) |
| 4.9 | 1e-11 | toward (1.5e-14) | 0.0215 (1.29) | 3.3e+05 | 3.07e+06 | 0.0172 (1.03) | 1, 1, 1 | captured | toward, captured, 0.0215, 3.3e+05 (w 0.0172) |
| 4.9 | 1e-12 | toward (1.5e-14) | 0.0215 (1.29) | 3.42e+05 | 9.7e+06 | 0.0167 (1.00) | 1, 1, 1 | captured | toward, captured, 0.0215, 3.42e+05 (w 0.0167) |
| 7 | 1e-08 | toward (1.6e-08) | 0.00498 (0.299) | 9.72e+04 | 9.72e+04 | 0.00498 (0.299) | 1, 1, 1 | captured | toward, captured, 0.00498, 9.72e+04 (w 0.00498) |
| 7 | 1e-09 | toward (1.7e-08) | 0.00495 (0.297) | 1.01e+05 | 3.07e+05 | 0.00495 (0.297) | 1, 1, 1 | captured | toward, captured, 0.00495, 1.01e+05 (w 0.00495) |
| 7 | 1e-10 | toward (1.8e-08) | 0.00496 (0.297) | 1.04e+05 | 9.7e+05 | 0.00496 (0.297) | 1, 1, 1 | captured | toward, captured, 0.00496, 1.04e+05 (w 0.00496) |
| 7 | 1e-11 | toward (1.8e-08) | 0.00494 (0.296) | 9.77e+04 | 3.07e+06 | 0.00494 (0.296) | 1, 1, 1 | captured | toward, captured, 0.00494, 9.77e+04 (w 0.00494) |
| 7 | 1e-12 | toward (1.8e-08) | 0.00494 (0.296) | 1.01e+05 | 9.7e+06 | 0.00494 (0.296) | 1, 1, 1 | captured | toward, captured, 0.00494, 1.01e+05 (w 0.00494) |
| 10 | 1e-08 | toward (6.6e-03) | 0.000811 (0.0486) | none | 9.73e+04 | — | 1, 1, 1 | not departed | toward, not departed, 0.000811, none |
| 10 | 1e-09 | toward (6.7e-03) | 0.000858 (0.0514) | 1.67e+05 | 3.07e+05 | 0.000858 (0.0514) | 1, 1, 1 | captured | toward, captured, 0.000858, 1.67e+05 (w 0.000858) |
| 10 | 1e-10 | toward (6.8e-03) | 0.000811 (0.0487) | 1.56e+05 | 9.7e+05 | 0.000811 (0.0487) | 1, 1, 1 | captured | toward, captured, 0.000811, 1.56e+05 (w 0.000811) |
| 10 | 1e-11 | toward (6.8e-03) | 0.000819 (0.0491) | 1.62e+05 | 3.07e+06 | 0.000819 (0.0491) | 1, 1, 1 | captured | toward, captured, 0.000819, 1.62e+05 (w 0.000819) |
| 10 | 1e-12 | toward (6.8e-03) | 0.000808 (0.0484) | 1.52e+05 | 9.7e+06 | 0.000808 (0.0484) | 1, 1, 1 | captured | toward, captured, 0.000808, 1.52e+05 (w 0.000808) |
| 17.3 | 1e-08 | toward (2.2e-01) | 5.78e-05 (0.00346) | none | 9.76e+04 | — | 1, 1, 1 | not departed | toward, not departed, 5.78e-05, none |
| 17.3 | 1e-09 | toward (2.2e-01) | 5.83e-05 (0.00350) | none | 3.07e+05 | — | 1, 1, 1 | not departed | toward, not departed, 5.83e-05, none |
| 17.3 | 1e-10 | toward (2.2e-01) | 6.57e-05 (0.00394) | 5.84e+05 | 9.7e+05 | 6.57e-05 (0.00394) | 1, 1, 1 | captured | toward, captured, 6.57e-05, 5.84e+05 (w 6.57e-05) |
| 17.3 | 1e-11 | toward (2.2e-01) | 6.9e-05 (0.00414) | 6.06e+05 | 3.07e+06 | 6.9e-05 (0.00414) | 1, 1, 1 | captured | toward, captured, 6.9e-05, 6.06e+05 (w 6.9e-05) |
| 17.3 | 1e-12 | toward (2.2e-01) | 7.74e-05 (0.00464) | 6.28e+05 | 9.7e+06 | 7.74e-05 (0.00464) | 1, 1, 1 | captured | toward, captured, 7.74e-05, 6.28e+05 (w 7.74e-05) |
| 24.5 | 1e-08 | toward (3.3e-01) | 1.39e-05 (0.000835) | none | 9.79e+04 | — | 1, 1, 1 | not departed | toward, not departed, 1.39e-05, none |
| 24.5 | 1e-09 | toward (3.3e-01) | 1.39e-05 (0.000836) | none | 3.08e+05 | — | 1, 1, 1 | not departed | toward, not departed, 1.39e-05, none |
| 24.5 | 1e-10 | toward (3.3e-01) | 1.46e-05 (0.000877) | none | 9.71e+05 | — | 1, 1, 1 | not departed | toward, not departed, 1.46e-05, none |
| 24.5 | 1e-11 | toward (3.2e-01) | 2.5e-05 (0.00150) | 1.85e+06 | 3.07e+06 | 2.5e-05 (0.00150) | 1, 1, 1 | captured | toward, captured, 2.5e-05, 1.85e+06 (w 2.5e-05) |
| 24.5 | 1e-12 | toward (3.2e-01) | 2.48e-05 (0.00149) | 1.73e+06 | 9.7e+06 | 2.48e-05 (0.00149) | 1, 1, 1 | captured | toward, captured, 2.48e-05, 1.73e+06 (w 2.48e-05) |

The two tolerances give the same outcome at 30 of 30 points.

The fold divided by: the configuration's own (FC), 0.01668.

C4M. The loop on the mirrored target (features 1 and 2 reflected across the face: feature 2 on the entrant's side, the runner-up displaced away from it), rtol as C4's, the absolute tolerance min(10⁻¹⁴, 10⁻⁵ε) (`c4_atol`: at 10⁻¹⁴ the state before the entry, at the masses' scale, is unresolved at ε = 10⁻¹⁷ and the loop enters thousands of s early), the outputs from s = −100 about the path's entry (s = (τ − τ_e)/ε^{1/2}; a row marked † has the grid of C4, from s = 10⁻⁷s₂, which can begin past the inner layer):

| x | ε | side | the pair's closest approach (σ) | max w held | departure s | s₂ | kernels per feature (final) | the runner-up at the end: w, distance from feature 0 (σ) | outcome |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 1e-10 | far | 1.8e-08 | 0.005422 | none | 9.7e+05 | 1, 1, 1 | 0.00541, 0.016 | held, unformed (far) |
| 7 | 4.5e-12 | far | 1.8e-08 | 0.007412 | none | 4.57e+06 | 1, 1, 1 | 0.007401, 0.006 | held, unformed (far) |
| 7 | 2e-12 | far | 1.8e-08 | 0.008708 | none | 6.86e+06 | 2, 1, 1 | 0.008697, 0.004 | bound |
| 7 | 1e-13 | far | 1.8e-08 | 0.02235 | none | 3.07e+07 | 2, 1, 1 | 0.02234, 0.001 | bound |
| 7 | 6.5e-17 | far | 1.8e-08 | 0.5817 | none | 1.2e+09 | 2, 1, 1 | 0.5812, 0.000 | bound |
| 7 | 3.7e-17 | far | 1.8e-08 | 0.8 | none | 1.59e+09 | 1, 1, 1 | 0.8, 0.000 | captured |
| 7 | 1e-17 | far | 1.8e-08 | 0.8 | none | 3.07e+09 | 1, 1, 1 | 0.8, 0.000 | captured |
| 10 | 6e-13 | far | 3.9e-03 | 0.006553 | none | 1.25e+07 | 1, 1, 1 | 0.006542, 0.004 | held, unformed (far) |
| 10 | 2.5e-13 | far | 2.0e-03 | 0.01046 | none | 1.94e+07 | 2, 1, 1 | 0.01045, 0.002 | bound |

C4W. C4's loop at (x, ε) = (4.9, 10⁻¹⁰) on the window x ∈ [-1, 0.8], y ∈ [-1, 1.5] (the right edge 4.92σ from feature 0, the top 8.33σ): the side the formed runner-up takes, its closest approach, largest held w, departure, the kernels per feature and the outcome, on both targets:

| target | x | ε | side | the pair's closest approach (σ) | max w held | departure s | s₂ | kernels per feature (final) | outcome |
|---|---|---|---|---|---|---|---|---|---|
| near tie | 4.9 | 1e-10 | far | 3.09e-14 | 0.02205 | none | 1.1e+06 | 2, 1, 1 | bound |
| near tie mirrored | 4.9 | 1e-10 | toward | 2.77e-14 | 0.02154 | 3.95e+05 | 1.09e+06 | 1, 1, 1 | captured |

## R. The race without ε (referee round 53)

C2's near tie at μ → 0 (1e-09) with only feature 0's entrant and runner-up active and features 1 and 2 in the residual: the inner problem (the entrant 400 behind its wall, the runner-up x behind it) in s = (τ − τ_e)/ε^{1/2} to 3e+08 — C4 without its horizon s₂ —, rtol 1e-13: the side the runner-up leaves feature 0 to (its displacement along e where it first passes 0.05σ from feature 0: toward feature 2, or away), the pair's closest approach in the inner layer, its coefficient after the layer (s = 10³), its departure s_d (2σ from feature 0) and its coefficient there:

| x | side | closest approach (σ; its share along e) | w after the layer | s_d | w at s_d | at s = 10⁵: w, distance (σ) |
|---|---|---|---|---|---|---|
| 2.5 | toward | 2.4e-15 (+0.88) | 0.1096 | 4.756e+07 | 0.02908 | 0.1102, 0.152 |
| 3 | toward | 3.5e-15 (+0.88) | 0.07927 | 2.346e+07 | 0.02547 | 0.07972, 0.188 |
| 3.5729 | toward | 5.0e-15 (+0.87) | 0.05393 | 9.198e+06 | 0.02133 | 0.05415, 0.245 |
| 4 | toward | 7.4e-15 (+0.88) | 0.04018 | 4.032e+06 | 0.0191 | 0.04024, 0.302 |
| 4.2554 | toward | 9.3e-15 (+0.88) | 0.03363 | 2.275e+06 | 0.01755 | 0.03357, 0.346 |
| 4.5 | toward | 1.1e-14 (+0.89) | 0.02832 | 1.21e+06 | 0.01722 | 0.02814, 0.399 |
| 4.9 | toward | 1.5e-14 (+0.88) | 0.02135 | 3.191e+05 | 0.01674 | 0.02081, 0.534 |
| 5.0709 | toward | 1.7e-14 (+0.88) | 0.01892 | 1.497e+05 | 0.01639 | 0.01797, 0.652 |
| 5.4772 | toward | 2.3e-14 (+0.89) | 0.01419 | 6.083e+04 | 0.0143 | 0.7, 6.079 |
| 6.0363 | toward | 5.8e-14 (+0.97) | 0.00956 | 8.402e+04 | 0.00966 | 0.7, 6.079 |
| 7 | toward | 1.8e-08 (+0.58) | 0.004905 | 9.566e+04 | 0.004949 | 0.7, 6.079 |
| 7.1919 | toward | 1.0e-07 (+0.57) | 0.004309 | 9.797e+04 | 0.004364 | 0.7, 6.079 |
| 8.589 | toward | 2.4e-04 (+0.52) | 0.001764 | 1.2e+05 | 0.00179 | 0.001775, 0.125 |
| 10 | toward | 6.8e-03 (+0.50) | 0.0007993 | 1.516e+05 | 0.0008189 | 0.0008035, 0.091 |
| 10.2121 | toward | 9.3e-03 (+0.50) | 0.0007168 | 1.573e+05 | 0.0007368 | 0.0007205, 0.092 |
| 12.1666 | toward | 5.6e-02 (+0.49) | 0.0002941 | 2.234e+05 | 0.0003293 | 0.0002955, 0.114 |
| 14.4841 | toward | 1.4e-01 (+0.48) | 0.0001271 | 3.458e+05 | 0.000145 | 0.0001277, 0.164 |
| 16.957 | toward | 2.1e-01 (+0.48) | 6.243e-05 | 5.451e+05 | 9.325e-05 | 6.266e-05, 0.223 |
| 17.2631 | toward | 2.2e-01 (+0.48) | 5.774e-05 | 5.754e+05 | 7.706e-05 | 5.796e-05, 0.229 |
| 17.3 | toward | 2.2e-01 (+0.48) | 5.721e-05 | 5.791e+05 | 7.433e-05 | 5.742e-05, 0.230 |
| 20.5557 | toward | 2.8e-01 (+0.48) | 2.772e-05 | 9.881e+05 | 5.606e-05 | 2.781e-05, 0.286 |
| 24.4949 | toward | 3.2e-01 (+0.48) | 1.384e-05 | 1.714e+06 | 2.441e-05 | 1.387e-05, 0.328 |
| 24.5 | toward | 3.2e-01 (+0.48) | 1.382e-05 | 1.715e+06 | 2.428e-05 | 1.386e-05, 0.328 |
| 30 | toward | 3.6e-01 (+0.48) | 6.458e-06 | 3.207e+06 | 1.569e-05 | 6.473e-06, 0.361 |

The rule — captured iff s_d(x) < s₂ = (τ₂ − τ_e)ε^{−1/2} (C3's s₂ from C4's τ₂ − τ_e = 9.699) — against the loops' outcomes: 74 of 74 agree, 1 capture ties:

| loops | x | ε | s₂ | s_d | predicted | the loop's outcome |
|---|---|---|---|---|---|---|
| C3 | 3 | 1e-08 | 9.7e+04 | 2.35e+07 | not departed | bound |
| C3 | 3.57 | 7.05e-09 | 1.16e+05 | 9.2e+06 | not departed | bound |
| C3 | 4.26 | 4.97e-09 | 1.38e+05 | 2.28e+06 | not departed | bound |
| C3 | 5.07 | 3.5e-09 | 1.64e+05 | 1.5e+05 | captured | captured |
| C3 | 5.48 | 3e-09 | 1.77e+05 | 6.08e+04 | captured | captured |
| C3 | 6.04 | 2.47e-09 | 1.95e+05 | 8.4e+04 | captured | captured |
| C3 | 7.19 | 1.74e-09 | 2.33e+05 | 9.8e+04 | captured | captured |
| C3 | 8.59 | 1.22e-09 | 2.78e+05 | 1.2e+05 | captured | captured |
| C3 | 10.2 | 8.63e-10 | 3.3e+05 | 1.57e+05 | captured | captured |
| C3 | 12.2 | 6.08e-10 | 3.93e+05 | 2.23e+05 | captured | captured |
| C3 | 14.5 | 4.29e-10 | 4.68e+05 | 3.46e+05 | captured | captured |
| C3 | 17 | 3.13e-10 | 5.48e+05 | 5.45e+05 | captured | capture tie |
| C3 | 17.3 | 3.02e-10 | 5.58e+05 | 5.75e+05 | not departed | unsettled |
| C3 | 20.6 | 2.13e-10 | 6.65e+05 | 9.88e+05 | not departed | unsettled |
| C3 | 24.5 | 1.5e-10 | 7.92e+05 | 1.71e+06 | not departed | unsettled |
| C4 | 3 | 1e-08 | 9.71e+04 | 2.35e+07 | not departed | bound |
| C4 | 3 | 1e-09 | 3.07e+05 | 2.35e+07 | not departed | bound |
| C4 | 3 | 1e-10 | 9.7e+05 | 2.35e+07 | not departed | bound |
| C4 | 3 | 1e-11 | 3.07e+06 | 2.35e+07 | not departed | bound |
| C4 | 3 | 1e-12 | 9.7e+06 | 2.35e+07 | not departed | bound |
| C4 | 4.9 | 1e-08 | 9.72e+04 | 3.19e+05 | not departed | bound |
| C4 | 4.9 | 1e-09 | 3.07e+05 | 3.19e+05 | not departed | bound |
| C4 | 4.9 | 1e-10 | 9.7e+05 | 3.19e+05 | captured | captured |
| C4 | 4.9 | 1e-11 | 3.07e+06 | 3.19e+05 | captured | captured |
| C4 | 4.9 | 1e-12 | 9.7e+06 | 3.19e+05 | captured | captured |
| C4 | 7 | 1e-08 | 9.72e+04 | 9.57e+04 | captured | captured |
| C4 | 7 | 1e-09 | 3.07e+05 | 9.57e+04 | captured | captured |
| C4 | 7 | 1e-10 | 9.7e+05 | 9.57e+04 | captured | captured |
| C4 | 7 | 1e-11 | 3.07e+06 | 9.57e+04 | captured | captured |
| C4 | 7 | 1e-12 | 9.7e+06 | 9.57e+04 | captured | captured |
| C4 | 10 | 1e-08 | 9.73e+04 | 1.52e+05 | not departed | not departed |
| C4 | 10 | 1e-09 | 3.07e+05 | 1.52e+05 | captured | captured |
| C4 | 10 | 1e-10 | 9.7e+05 | 1.52e+05 | captured | captured |
| C4 | 10 | 1e-11 | 3.07e+06 | 1.52e+05 | captured | captured |
| C4 | 10 | 1e-12 | 9.7e+06 | 1.52e+05 | captured | captured |
| C4 | 17.3 | 1e-08 | 9.76e+04 | 5.79e+05 | not departed | not departed |
| C4 | 17.3 | 1e-09 | 3.07e+05 | 5.79e+05 | not departed | not departed |
| C4 | 17.3 | 1e-10 | 9.7e+05 | 5.79e+05 | captured | captured |
| C4 | 17.3 | 1e-11 | 3.07e+06 | 5.79e+05 | captured | captured |
| C4 | 17.3 | 1e-12 | 9.7e+06 | 5.79e+05 | captured | captured |
| C4 | 24.5 | 1e-08 | 9.79e+04 | 1.72e+06 | not departed | not departed |
| C4 | 24.5 | 1e-09 | 3.08e+05 | 1.72e+06 | not departed | not departed |
| C4 | 24.5 | 1e-10 | 9.71e+05 | 1.72e+06 | not departed | not departed |
| C4 | 24.5 | 1e-11 | 3.07e+06 | 1.72e+06 | captured | captured |
| C4 | 24.5 | 1e-12 | 9.7e+06 | 1.72e+06 | captured | captured |
| C4b | 3 | 1e-08 | 9.71e+04 | 2.35e+07 | not departed | bound |
| C4b | 3 | 1e-09 | 3.07e+05 | 2.35e+07 | not departed | bound |
| C4b | 3 | 1e-10 | 9.7e+05 | 2.35e+07 | not departed | bound |
| C4b | 3 | 1e-11 | 3.07e+06 | 2.35e+07 | not departed | bound |
| C4b | 3 | 1e-12 | 9.7e+06 | 2.35e+07 | not departed | bound |
| C4b | 4.9 | 1e-08 | 9.72e+04 | 3.19e+05 | not departed | bound |
| C4b | 4.9 | 1e-09 | 3.07e+05 | 3.19e+05 | not departed | bound |
| C4b | 4.9 | 1e-10 | 9.7e+05 | 3.19e+05 | captured | captured |
| C4b | 4.9 | 1e-11 | 3.07e+06 | 3.19e+05 | captured | captured |
| C4b | 4.9 | 1e-12 | 9.7e+06 | 3.19e+05 | captured | captured |
| C4b | 7 | 1e-08 | 9.72e+04 | 9.57e+04 | captured | captured |
| C4b | 7 | 1e-09 | 3.07e+05 | 9.57e+04 | captured | captured |
| C4b | 7 | 1e-10 | 9.7e+05 | 9.57e+04 | captured | captured |
| C4b | 7 | 1e-11 | 3.07e+06 | 9.57e+04 | captured | captured |
| C4b | 7 | 1e-12 | 9.7e+06 | 9.57e+04 | captured | captured |
| C4b | 10 | 1e-08 | 9.73e+04 | 1.52e+05 | not departed | not departed |
| C4b | 10 | 1e-09 | 3.07e+05 | 1.52e+05 | captured | captured |
| C4b | 10 | 1e-10 | 9.7e+05 | 1.52e+05 | captured | captured |
| C4b | 10 | 1e-11 | 3.07e+06 | 1.52e+05 | captured | captured |
| C4b | 10 | 1e-12 | 9.7e+06 | 1.52e+05 | captured | captured |
| C4b | 17.3 | 1e-08 | 9.76e+04 | 5.79e+05 | not departed | not departed |
| C4b | 17.3 | 1e-09 | 3.07e+05 | 5.79e+05 | not departed | not departed |
| C4b | 17.3 | 1e-10 | 9.7e+05 | 5.79e+05 | captured | captured |
| C4b | 17.3 | 1e-11 | 3.07e+06 | 5.79e+05 | captured | captured |
| C4b | 17.3 | 1e-12 | 9.7e+06 | 5.79e+05 | captured | captured |
| C4b | 24.5 | 1e-08 | 9.79e+04 | 1.72e+06 | not departed | not departed |
| C4b | 24.5 | 1e-09 | 3.08e+05 | 1.72e+06 | not departed | not departed |
| C4b | 24.5 | 1e-10 | 9.71e+05 | 1.72e+06 | not departed | not departed |
| C4b | 24.5 | 1e-11 | 3.07e+06 | 1.72e+06 | captured | captured |
| C4b | 24.5 | 1e-12 | 9.7e+06 | 1.72e+06 | captured | captured |

The mirrored target (feature 2 across the face: the runner-up held on the far side): the side, the closest approach, the coefficient after the layer and at s = 10000, 100000, 1e+06, 1e+07, 1e+08, 3e+08, the distance from feature 0 at the end, the first s at which a coefficient below θ|A| has grown past it, the departure:

| x | side | closest (σ) | w after the layer | w at the marks | distance at the end (σ) | s at θ\|A\| | departure |
|---|---|---|---|---|---|---|---|
| 3 | toward | 3.3e-15 | 0.0793 | 0.0798, 0.0797, 0.0784, 0.0632, 0.7, 0.7 | 6.079 | — | 2.35e+07 |
| 4.9 | toward | 1.5e-14 | 0.02136 | 0.0215, 0.0208, 0.7, 0.7, 0.7, 0.7 | 6.079 | — | 3.2e+05 |
| 7 | far | 1.8e-08 | 0.004907 | 0.00494, 0.00495, 0.00543, 0.0105, 0.0593, 0.153 | 0.105 | 5.71e+06 | none |
| 10 | far | 6.8e-03 | 0.0007995 | 0.000803, 0.000804, 0.00106, 0.00517, 0.0544, 0.148 | 0.107 | 1.57e+07 | none |
| 17.3 | far | 1.3e-01 | 5.722e-05 | 5.74e-05, 5.74e-05, 6.24e-05, 0.00182, 0.0509, 0.145 | 0.108 | 2.28e+07 | none |

The mirror past R's horizon (referee round 58): the coefficients at the end, the s at which the runner-up outgrows the entrant (the swap), the entrant's departure (2σ from feature 0) and the feature it ends nearest (its distance), the runner-up's distance from feature 0 at the end:

| x | s_end | w: runner-up, entrant | the swap s | the entrant's departure s | the entrant ends at feature (σ) | the runner-up's distance (σ) |
|---|---|---|---|---|---|---|
| 7 | 3e+09 | 0.8, 0.7 | 8.63e+08 | 1.38e+09 | 2 (0.000) | 0.000 |
| 10 | 3e+09 | 0.8, 0.7 | 8.73e+08 | 1.39e+09 | 2 (0.000) | 0.000 |

The far rule's window in ε (referee round 62): s₂ = (τ₂ − τ_e)ε^{−1/2} with τ₂ − τ_e = 9.696 (C4M's path, 9.696–9.699), so s_θ < s₂ < s_{d,e} holds for ((τ₂ − τ_e)/s_{d,e})² < ε < ((τ₂ − τ_e)/s_θ)²: x = 7: 4.9e-17 < ε < 2.9e-12; x = 10: 4.8e-17 < ε < 3.8e-13.

The side against the tolerance (to s = 3e+06): per gap the closest approach at the tightest tolerance and the side at rtol 1e-09, 1e-10, 1e-11, 1e-12, 1e-13 (None: the pair merged — the runner-up never 0.05σ from feature 0):

| x | closest approach (σ) | sides |
|---|---|---|
| 2.5 | 2.0e-15 | far, toward, toward, toward, toward |
| 3 | 3.5e-15 | toward, far, toward, toward, toward |
| 3.5 | 5.2e-15 | far, far, toward, far, toward |
| 4 | 7.4e-15 | far, far, far, toward, toward |
| 4.5 | 1.1e-14 | toward, far, toward, toward, toward |
| 4.9 | 1.5e-14 | far, far, far, toward, toward |
| 5.5 | 2.4e-14 | None, toward, far, toward, toward |
| 6 | 4.3e-14 | far, None, toward, toward, toward |

Runs not toward, per tolerance: 1e-09: 6 of 8; 1e-10: 6 of 8; 1e-11: 3 of 8; 1e-12: 1 of 8; 1e-13: 0 of 8.

FC. The static fold in C4's own geometry (Part C's samples; the entrant and the runner-up's centre free, its coefficient held): toward feature 2: the fold at w = 0.016679, the rest there 0.842σ from feature 0 (+0.840σ along e); away from it: no fold down to w = 0.003 (the rest 0.591σ, -0.578σ along e); toward, feature 1 removed: the fold at w = 0.017196, the rest there 0.853σ from feature 0 (+0.853σ along e).

## RP. The relative mode: the pin (referee rounds 58 and 62)

R's race (the near tie at μ → 0, the runner-up x behind the entrant) with the sums over the samples in float64 or in the 80-bit long double (unit roundoff 5.4·10⁻²⁰; the state float64), on Part C's window with its top edge at 1 (the block 5.00σ below it), nearer or farther, on the window with every edge at least 11.67σ from feature 0 (x ∈ [-2, 2], y ∈ [-1.5, 2]: no edge within reach, the block's rest moving with W only by the samples' own float64 coordinates — the same in the long double and in exact sums — and by the aliasing, RD) and with the right edge 4.92σ from it (x ∈ [-1, 0.8], y ∈ [-1, 1.5]); N the samples (⟨·⟩ averages over them, so a window of N samples runs at N₀/N the speed of Part C's, N₀ = 1681). The relative offset D = c_r − c_e (σ) at s = 500, 700, 1000 against the law integrated from s = 400 with the block's centroid c_b = (w_e c_e + w c_r)/W slaved to its rest — dD/ds = κ[Λ₀(2Gm ww_e/W)⟨r, ∇²_cφ(c̄)⟩D + (2(w − w_e)/W)(ρ′D + (∂c̄*/∂W)dW/ds)], κ = W²/(2(w² + w_e²)), ρ = (w_e − w)/(2W) (the midpoint c̄ = c_b + ρD), ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩, F = ⟨f* − Wφ(c̄), ∇_cφ(c̄)⟩ —; ∂c̄*/∂W (σ per unit W) at s = 700; the least |D| after the wall before the departure and its share along e; the side, the climb's end (|D| back at 0.01σ) and the departure s_d (in brackets times N₀/N: in Part C's time) and |D|'s growth rate on s ∈ [3·10³, 10⁴]:

| target | x | window | N | sums | D at s = 500, 700, 1000 (x, y) | the law | ∂c̄*/∂W (x, y) | least \|D\| (its share along e) | side | climb's end (×N₀/N) | s_d (×N₀/N) | growth rate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| near tie | 3 | top 0.95 (4.67σ) | 1640 | float64 | (+3.24e-16, -1.13e-13); (+1.39e-15, -1.14e-13); (+1.55e-14, -1.24e-13) | (+3.99e-16, -1.13e-13); (+1.31e-15, -1.14e-13); (+1.48e-14, -1.24e-13) | (-3.62e-14, +3.37e-11) | 1.10e-13 (+0.88) | toward | 8679 (8896) | 2.289e+07 (2.346e+07) | 0.0033 |
| near tie | 4.9 | top 0.95 (4.67σ) | 1640 | float64 | (+1.57e-15, -4.54e-13); (+4.76e-15, -4.60e-13); (+3.42e-14, -4.81e-13) | (+1.56e-15, -4.54e-13); (+4.59e-15, -4.60e-13); (+3.33e-14, -4.80e-13) | (-3.59e-14, +3.37e-11) | 4.52e-13 (+0.88) | toward | 2.948e+04 (3.021e+04) | 3.073e+05 (3.15e+05) | 0.000842 |
| near tie | 6.04 | top 0.95 (4.67σ) | 1640 | float64 | (+3.63e-09, -9.92e-11); (+5.09e-13, -1.08e-12); (+9.12e-14, -1.08e-12) | (+3.50e-09, -9.55e-11); (+4.22e-13, -1.07e-12); (+8.05e-14, -1.08e-12) | (-3.61e-14, +3.37e-11) | 1.07e-12 (+0.92) | toward | 6.339e+04 (6.497e+04) | 7.489e+04 (7.676e+04) | 0.000369 |
| near tie | 3 | top 1 (5.00σ) | 1681 | float64 | (+0.00e+00, -4.07e-15); (-4.63e-17, -4.26e-15); (+4.63e-16, -3.89e-15) | (+1.33e-17, -3.81e-15); (+3.91e-17, -3.84e-15); (+3.93e-16, -4.09e-15) | (-1.11e-15, +1.14e-12) | 3.15e-15 (+0.89) | toward | 9900 (9900) | 2.346e+07 (2.346e+07) | 0.00326 |
| near tie | 4.5 | top 1 (5.00σ) | 1681 | float64 | (-4.63e-17, -1.13e-14); (+1.39e-16, -1.11e-14); (+8.79e-16, -1.24e-14) | (+3.87e-17, -1.15e-14); (+1.06e-16, -1.16e-14); (+7.96e-16, -1.21e-14) | (-1.04e-15, +1.14e-12) | 1.09e-14 (+0.88) | toward | 2.619e+04 (2.619e+04) | 1.21e+06 (1.21e+06) | 0.0011 |
| near tie | 4.9 | top 1 (5.00σ) | 1681 | float64 | (+9.25e-17, -1.52e-14); (+2.31e-16, -1.57e-14); (+8.33e-16, -1.65e-14) | (+1.54e-16, -1.53e-14); (+1.38e-16, -1.55e-14); (+9.41e-16, -1.61e-14) | (-1.04e-15, +1.14e-12) | 1.50e-14 (+0.89) | toward | 3.426e+04 (3.426e+04) | 3.191e+05 (3.191e+05) | 0.000829 |
| near tie | 5.07 | top 1 (5.00σ) | 1681 | float64 | (+4.76e-15, -1.74e-14); (+2.78e-16, -1.72e-14); (+9.25e-16, -1.91e-14) | (+4.22e-15, -1.75e-14); (+1.54e-16, -1.76e-14); (+1.01e-15, -1.82e-14) | (-1.04e-15, +1.14e-12) | 1.70e-14 (+0.88) | toward | 3.856e+04 (3.856e+04) | 1.504e+05 (1.504e+05) | 0.000731 |
| near tie | 5.48 | top 1 (5.00σ) | 1681 | float64 | (+6.32e-12, -2.51e-13); (+1.85e-16, -2.35e-14); (+1.34e-15, -2.41e-14) | (+5.84e-12, -2.35e-13); (+2.05e-16, -2.36e-14); (+1.18e-15, -2.43e-14) | (-1.08e-15, +1.14e-12) | 2.23e-14 (+0.94) | toward | 5.09e+04 (5.09e+04) | 6.101e+04 (6.101e+04) | 0.000542 |
| near tie | 6.04 | top 1 (5.00σ) | 1681 | float64 | (+7.06e-09, -1.80e-10); (+7.28e-13, -8.23e-14); (+4.49e-14, -4.22e-14) | (+6.82e-09, -1.74e-10); (+6.07e-13, -7.45e-14); (+3.44e-14, -4.08e-14) | (-1.26e-15, +1.14e-12) | 5.93e-14 (+0.97) | toward | 7.26e+04 (7.26e+04) | 8.439e+04 (8.439e+04) | 0.000362 |
| near tie | 7 | top 1 (5.00σ) | 1681 | float64 | (+1.02e-05, -1.48e-07); (+8.46e-08, -2.94e-09); (+1.98e-08, -1.45e-09) | (+1.01e-05, -1.46e-07); (+7.81e-08, -2.71e-09); (+1.74e-08, -1.27e-09) | (-1.17e-15, +1.14e-12) | 1.81e-08 (+0.58) | toward | 7.483e+04 (7.483e+04) | 9.566e+04 (9.566e+04) | 0.000161 |
| near tie | 4.9 | top 1 (5.00σ) | 1681 | long double | (+9.25e-17, -1.54e-14); (+1.85e-16, -1.55e-14); (+5.55e-16, -1.65e-14) | (+1.55e-16, -1.53e-14); (+1.40e-16, -1.55e-14); (+9.42e-16, -1.61e-14) | (-1.22e-15, +1.14e-12) | 1.50e-14 (+0.88) | toward | 3.441e+04 (3.441e+04) | 3.192e+05 (3.192e+05) | 0.000828 |
| near tie | 3 | top 1.5 (8.33σ) | 2091 | float64 | (+4.63e-17, +3.70e-16); (-4.63e-17, +0.00e+00); (+2.31e-16, +0.00e+00) | (-1.36e-19, +9.18e-20); (-2.22e-19, -4.59e-20); (-3.37e-19, +1.08e-19) | (+6.65e-17, +1.95e-17) | 0.00e+00 (—) | toward | 1.457e+04 (1.172e+04) | 2.919e+07 (2.347e+07) | 0.00263 |
| near tie | 4.9 | top 1.5 (8.33σ) | 2091 | float64 | (+6.70e-11, -1.86e-12); (-4.63e-17, -1.85e-16); (+4.63e-17, -1.85e-16) | (+5.94e-11, -1.66e-12); (-8.60e-19, +4.08e-20); (-4.18e-19, +1.85e-19) | (+6.65e-17, +1.95e-17) | 0.00e+00 (—) | far | 4.787e+04 (3.849e+04) | none | 0.000676 |
| near tie | 6.04 | top 1.5 (8.33σ) | 2091 | float64 | (+2.91e-06, -4.12e-08); (+1.96e-10, -6.90e-12); (+1.18e-13, -1.09e-14) | (+2.86e-06, -4.05e-08); (+1.79e-10, -6.30e-12); (+9.16e-14, -8.34e-15) | (+1.78e-17, +1.17e-17) | 3.70e-14 (+0.64) | toward | 9.295e+04 (7.473e+04) | 1.076e+05 (8.651e+04) | 0.00027 |
| near tie | 4.9 | top 1.5 (8.33σ) | 2091 | long double | (+6.70e-11, -1.86e-12); (+9.25e-17, -3.70e-16); (+0.00e+00, +0.00e+00) | (+5.94e-11, -1.66e-12); (+2.15e-20, -2.04e-21); (-6.64e-21, +5.77e-22) | (+5.37e-19, +2.35e-20) | 0.00e+00 (—) | far | 5.245e+04 (4.216e+04) | none | 0.00046 |
| near tie | 3 | right edge | 1887 | float64 | (-9.07e-15, +1.85e-16); (-9.02e-15, +3.70e-16); (-9.16e-15, +5.55e-16) | (-9.08e-15, +2.67e-17); (-9.10e-15, +5.47e-17); (-9.14e-15, +4.31e-16) | (+2.72e-12, -2.91e-15) | 1.49e-15 (+0.58) | far | 1.094e+04 (9744) | none | 0.0029 |
| near tie | 4.9 | right edge | 1887 | float64 | (+2.11e-13, -9.07e-15); (-3.69e-14, +0.00e+00); (-3.71e-14, +9.25e-16) | (+1.77e-13, -7.90e-15); (-3.68e-14, +2.03e-16); (-3.70e-14, +1.23e-15) | (+2.72e-12, -2.91e-15) | 9.34e-15 (+0.66) | far | 3.783e+04 (3.37e+04) | none | 0.000734 |
| near tie | 3 | far edges | 5751 | float64 | (+9.46e-07, -1.29e-08); (+1.55e-09, -3.25e-11); (+9.25e-17, +0.00e+00) | (+7.58e-07, -1.05e-08); (+6.90e-10, -1.51e-11); (+6.48e-18, +9.49e-20) | (+3.29e-16, -6.89e-17) | 0.00e+00 (—) | far | 3.775e+04 (1.104e+04) | none | 0.000743 |
| near tie | 4.5 | far edges | 5751 | float64 | (+8.05e-04, -5.05e-06); (+6.93e-05, -6.27e-07); (+1.50e-07, -2.45e-09) | (+7.72e-04, -4.89e-06); (+6.06e-05, -5.57e-07); (+1.10e-07, -1.83e-09) | (+2.67e-16, -7.67e-17) | 0.00e+00 (—) | toward | 1.029e+05 (3.006e+04) | 4.154e+06 (1.214e+06) | 0.00026 |
| near tie | 4.9 | far edges | 5751 | float64 | (+2.15e-03, -1.12e-05); (+3.10e-04, -2.29e-06); (+2.68e-06, -3.50e-08) | (+2.09e-03, -1.10e-05); (+2.84e-04, -2.13e-06); (+2.21e-06, -2.94e-08) | (+3.04e-16, -1.14e-16) | 0.00e+00 (—) | toward | 1.362e+05 (3.982e+04) | 1.111e+06 (3.246e+05) | 0.00015 |
| near tie | 5.07 | far edges | 5751 | float64 | (+3.07e-03, -1.48e-05); (+5.30e-04, -3.60e-06); (+7.39e-06, -8.78e-08) | (+3.00e-03, -1.45e-05); (+4.93e-04, -3.39e-06); (+6.32e-06, -7.63e-08) | (+1.21e-16, -4.13e-17) | 0.00e+00 (—) | None | — | none | 0.000154 |
| near tie | 5.48 | far edges | 5751 | float64 | (+6.39e-03, -2.57e-05); (+1.58e-03, -8.84e-06); (+5.66e-05, -5.41e-07) | (+6.29e-03, -2.54e-05); (+1.51e-03, -8.51e-06); (+5.14e-05, -4.97e-07) | (+4.63e-16, -2.12e-16) | 0.00e+00 (—) | far | 2.099e+05 (6.136e+04) | none | -9.95e-06 |
| near tie | 6.04 | far edges | 5751 | float64 | (+1.39e-02, -4.43e-05); (+4.92e-03, -2.14e-05); (+4.45e-04, -3.21e-06) | (+1.38e-02, -4.40e-05); (+4.79e-03, -2.10e-05); (+4.22e-04, -3.07e-06) | (+2.16e-16, -7.86e-17) | 3.67e-14 (+0.63) | toward | 2.556e+05 (7.471e+04) | 2.959e+05 (8.65e+04) | -2.59e-06 |
| near tie | 7 | far edges | 5751 | float64 | (+3.53e-02, -7.66e-05); (+1.84e-02, -5.34e-05); (+4.37e-03, -2.02e-05) | (+3.52e-02, -7.64e-05); (+1.82e-02, -5.30e-05); (+4.28e-03, -1.98e-05) | (+3.97e-16, -8.47e-17) | 1.81e-08 (+0.58) | toward | 2.56e+05 (7.483e+04) | 3.273e+05 (9.566e+04) | -1.87e-06 |
| near tie mirrored | 4.9 | top 1 (5.00σ) | 1681 | float64 | (+0.00e+00, -1.50e-14); (-1.85e-16, -1.55e-14); (-9.25e-16, -1.61e-14) | (+5.45e-17, -1.53e-14); (-1.41e-16, -1.55e-14); (-9.41e-16, -1.61e-14) | (+1.32e-15, +1.14e-12) | 1.50e-14 (+0.88) | toward | 3.431e+04 (3.431e+04) | 3.198e+05 (3.198e+05) | 0.000824 |
| near tie mirrored | 3 | right edge | 1887 | float64 | (-8.88e-15, +0.00e+00); (-9.34e-15, -1.85e-16); (-8.97e-15, -5.55e-16) | (-9.08e-15, -2.64e-17); (-9.10e-15, -5.45e-17); (-9.15e-15, -4.31e-16) | (+2.72e-12, +2.85e-15) | 1.22e-15 (-0.61) | toward | 1.091e+04 (9716) | 2.636e+07 (2.348e+07) | 0.0029 |
| near tie mirrored | 4.9 | right edge | 1887 | float64 | (+2.09e-13, +9.44e-15); (-3.69e-14, -3.70e-16); (-3.70e-14, -7.40e-16) | (+1.76e-13, +8.07e-15); (-3.68e-14, -2.04e-16); (-3.70e-14, -1.23e-15) | (+2.72e-12, +2.98e-15) | 9.13e-15 (-0.63) | toward | 3.777e+04 (3.365e+04) | 3.583e+05 (3.191e+05) | 0.000734 |
| near tie mirrored | 3 | far edges | 5751 | float64 | (+9.44e-07, +1.39e-08); (+1.55e-09, +3.40e-11); (+0.00e+00, -1.85e-16) | (+7.56e-07, +1.13e-08); (+6.87e-10, +1.58e-11); (+6.18e-18, +3.76e-19) | (+2.82e-16, -7.53e-18) | 0.00e+00 (—) | far | 3.648e+04 (1.066e+04) | none | 0.000897 |
| near tie mirrored | 4.5 | far edges | 5751 | float64 | (+8.03e-04, +5.90e-06); (+6.92e-05, +7.00e-07); (+1.49e-07, +2.61e-09) | (+7.71e-04, +5.70e-06); (+6.05e-05, +6.20e-07); (+1.09e-07, +1.94e-09) | (+3.49e-16, -2.33e-16) | 0.00e+00 (—) | far | 1.029e+05 (3.007e+04) | none | 8.09e-05 |
| near tie mirrored | 4.9 | far edges | 5751 | float64 | (+2.15e-03, +1.35e-05); (+3.10e-04, +2.62e-06); (+2.67e-06, +3.78e-08) | (+2.09e-03, +1.32e-05); (+2.84e-04, +2.43e-06); (+2.21e-06, +3.17e-08) | (+4.28e-16, -1.33e-17) | 0.00e+00 (—) | far | 1.454e+05 (4.249e+04) | none | 0.000121 |
| near tie mirrored | 5.07 | far edges | 5751 | float64 | (+3.06e-03, +1.80e-05); (+5.29e-04, +4.16e-06); (+7.37e-06, +9.55e-08) | (+2.99e-03, +1.77e-05); (+4.92e-04, +3.91e-06); (+6.31e-06, +8.28e-08) | (+2.67e-16, +5.90e-17) | 0.00e+00 (—) | far | 1.47e+05 (4.296e+04) | none | 0.000254 |
| near tie mirrored | 5.48 | far edges | 5751 | float64 | (+6.38e-03, +3.25e-05); (+1.58e-03, +1.05e-05); (+5.65e-05, +6.00e-07) | (+6.28e-03, +3.21e-05); (+1.50e-03, +1.01e-05); (+5.13e-05, +5.50e-07) | (+3.42e-16, -7.41e-17) | 0.00e+00 (—) | far | 2.101e+05 (6.141e+04) | none | 0.000259 |
| near tie mirrored | 6.04 | far edges | 5751 | float64 | (+1.39e-02, +5.91e-05); (+4.92e-03, +2.66e-05); (+4.44e-04, +3.67e-06) | (+1.38e-02, +5.86e-05); (+4.79e-03, +2.60e-05); (+4.22e-04, +3.51e-06) | (+3.00e-16, -3.49e-17) | 3.65e-14 (-0.64) | far | 2.558e+05 (7.478e+04) | none | -2.41e-06 |
| near tie mirrored | 7 | far edges | 5751 | float64 | (+3.53e-02, +1.14e-04); (+1.84e-02, +7.30e-05); (+4.37e-03, +2.48e-05) | (+3.51e-02, +1.14e-04); (+1.81e-02, +7.22e-05); (+4.28e-03, +2.43e-05) | (+4.87e-16, -9.74e-17) | 1.80e-08 (-0.58) | far | 2.566e+05 (7.5e+04) | none | -1.78e-06 |

On Part C's window, the law's D_y over the run's at s = 500, 700, 1000 — the slaved centroid's law (round 58's, κ = 1 and ρ′ = 0, in brackets): x = 3: 0.937 (0.937), 0.902 (0.902), 1.053 (1.063); x = 4.5: 1.018 (1.018), 1.045 (1.046), 0.976 (0.988); x = 4.9: 1.010 (1.010), 0.987 (0.988), 0.978 (0.990); x = 5.07: 1.006 (0.995), 1.021 (1.021), 0.954 (0.965); x = 5.48: 0.933 (0.092), 1.004 (1.005), 1.011 (1.023); x = 6.04: 0.965 (0.004), 0.904 (0.427), 0.966 (0.865); x = 7: 0.988 (0.053), 0.920 (0.001), 0.877 (0.000).
The run's D is the difference of float64 centres, so D_y is quantised in steps of the runner-up's c_y's last place, 1.85e-16σ, 1.85e-16σ, 1.85e-16σ, 1.85e-16σ, 1.85e-16σ, 1.85e-16σ, 1.85e-16σ (x = 3, 4.5, 4.9, 5.07, 5.48, 6.04, 7): the pin at x = 3 is 23 such steps, the law there within 2.3 of them at s = 700.

The rate (referee rounds 62 and 66), at every float64 row: the Jacobian's eigenvalues nearest the relative mode's (the full two-kernel system, Λ₀ times `critical_jac_nd`) against the slaved law's κλ_D·eig⟨r, ∇²_cφ⟩ — its error over the largest rate — and round 58's λ_D·eig⟨r, ∇²_cφ⟩ (λ_D = Λ₀·2Gm ww_e/W): near tie, x = 3, top 0.95 (4.67σ), s = 451: -0.5374, -0.5352 against -0.5639, -0.5616 (4.9%; κ = 0.610; ε_b 0.046, 0.046; the (D, y) system -0.5374, -0.5352 (6.2e-16); round 58: -0.9251, -0.9212); near tie, x = 3, top 0.95 (4.67σ), s = 4992: -0.001681, +0.003338 against -0.001681, +0.003338 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.001681, +0.003338 (5.6e-08); round 58: -0.002758, +0.005477); near tie, x = 4.9, top 0.95 (4.67σ), s = 451.2: -0.1409, -0.1404 against -0.143, -0.1425 (1.5%; κ = 0.528; ε_b 0.015, 0.015; the (D, y) system -0.1409, -0.1404 (2.0e-15); round 58: -0.2709, -0.2698); near tie, x = 4.9, top 0.95 (4.67σ), s = 4993: -0.0004239, +0.0008417 against -0.0004239, +0.0008417 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0004239, +0.0008417 (1.6e-11); round 58: -0.0008035, +0.001595); near tie, x = 6.04, top 0.95 (4.67σ), s = 451.2: -0.06416, -0.06391 against -0.06461, -0.06435 (0.7%; κ = 0.513; ε_b 0.007, 0.007; the (D, y) system -0.06416, -0.06391 (1.8e-10); round 58: -0.126, -0.1255); near tie, x = 6.04, top 0.95 (4.67σ), s = 4993: -0.0001865, +0.0003704 against -0.0001865, +0.0003704 (0.0%; κ = 0.512; ε_b -0.000, 0.000; the (D, y) system -0.0001865, +0.0003704 (6.7e-12); round 58: -0.0003642, +0.0007232); near tie, x = 3, top 1 (5.00σ), s = 448.6: -0.5157, -0.5137 against -0.5427, -0.5406 (5.3%; κ = 0.610; ε_b 0.049, 0.049; the (D, y) system -0.5157, -0.5137 (4.3e-16); round 58: -0.8903, -0.8869); near tie, x = 3, top 1 (5.00σ), s = 4994: -0.00164, +0.003257 against -0.00164, +0.003256 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.00164, +0.003257 (1.4e-09); round 58: -0.002691, +0.005343); near tie, x = 4.5, top 1 (5.00σ), s = 448.8: -0.1791, -0.1784 against -0.1828, -0.1821 (2.1%; κ = 0.537; ε_b 0.020, 0.020; the (D, y) system -0.1791, -0.1784 (6.2e-16); round 58: -0.3402, -0.3389); near tie, x = 4.5, top 1 (5.00σ), s = 4994: -0.0005531, +0.001098 against -0.0005531, +0.001098 (0.0%; κ = 0.537; ε_b -0.000, 0.000; the (D, y) system -0.0005531, +0.001098 (1.0e-12); round 58: -0.00103, +0.002046); near tie, x = 4.9, top 1 (5.00σ), s = 448.9: -0.1351, -0.1346 against -0.1373, -0.1368 (1.6%; κ = 0.528; ε_b 0.016, 0.016; the (D, y) system -0.1351, -0.1346 (1.4e-15); round 58: -0.26, -0.259); near tie, x = 4.9, top 1 (5.00σ), s = 4994: -0.0004136, +0.0008212 against -0.0004136, +0.0008212 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0004136, +0.0008212 (4.2e-13); round 58: -0.0007839, +0.001556); near tie, x = 5.07, top 1 (5.00σ), s = 448.9: -0.12, -0.1195 against -0.1217, -0.1212 (1.4%; κ = 0.525; ε_b 0.014, 0.014; the (D, y) system -0.12, -0.1195 (4.2e-15); round 58: -0.2318, -0.2309); near tie, x = 5.07, top 1 (5.00σ), s = 4995: -0.0003656, +0.0007259 against -0.0003656, +0.0007259 (0.0%; κ = 0.524; ε_b -0.000, 0.000; the (D, y) system -0.0003656, +0.0007259 (2.7e-13); round 58: -0.0006972, +0.001384); near tie, x = 5.48, top 1 (5.00σ), s = 448.9: -0.09028, -0.08995 against -0.09125, -0.09091 (1.1%; κ = 0.519; ε_b 0.011, 0.011; the (D, y) system -0.09028, -0.08995 (1.7e-12); round 58: -0.1759, -0.1753); near tie, x = 5.48, top 1 (5.00σ), s = 4995: -0.0002718, +0.0005397 against -0.0002718, +0.0005397 (0.0%; κ = 0.518; ε_b -0.000, 0.000; the (D, y) system -0.0002718, +0.0005397 (9.6e-14); round 58: -0.0005246, +0.001042); near tie, x = 6.04, top 1 (5.00σ), s = 448.9: -0.06173, -0.0615 against -0.06218, -0.06195 (0.7%; κ = 0.513; ε_b 0.007, 0.007; the (D, y) system -0.06173, -0.0615 (3.9e-10); round 58: -0.1213, -0.1208); near tie, x = 6.04, top 1 (5.00σ), s = 4995: -0.000182, +0.0003614 against -0.000182, +0.0003613 (0.0%; κ = 0.512; ε_b -0.000, 0.000; the (D, y) system -0.000182, +0.0003614 (1.8e-12); round 58: -0.0003553, +0.0007056); near tie, x = 7, top 1 (5.00σ), s = 448.9: -0.03317, -0.03304 against -0.0333, -0.03318 (0.4%; κ = 0.507; ε_b 0.004, 0.004; the (D, y) system -0.03317, -0.03304 (1.2e-07); round 58: -0.06571, -0.06546); near tie, x = 7, top 1 (5.00σ), s = 4995: -9.306e-05, +0.0001848 against -9.306e-05, +0.0001848 (0.0%; κ = 0.506; ε_b -0.000, 0.000; the (D, y) system -9.306e-05, +0.0001848 (3.6e-08); round 58: -0.0001838, +0.000365); near tie, x = 3, top 1.5 (8.33σ), s = 450.5: -0.3128, -0.3119 against -0.335, -0.334 (7.1%; κ = 0.610; ε_b 0.064, 0.064; the (D, y) system -0.3128, -0.3119 (1.2e-15); round 58: -0.5492, -0.5476); near tie, x = 3, top 1.5 (8.33σ), s = 5012: -0.001318, +0.002618 against -0.001318, +0.002618 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.001318, +0.002618 (4.3e-14); round 58: -0.002163, +0.004295); near tie, x = 4.9, top 1.5 (8.33σ), s = 450.8: -0.0826, -0.08236 against -0.08435, -0.08411 (2.1%; κ = 0.529; ε_b 0.021, 0.021; the (D, y) system -0.0826, -0.08236 (1.0e-11); round 58: -0.1595, -0.159); near tie, x = 4.9, top 1.5 (8.33σ), s = 5013: -0.0003325, +0.0006602 against -0.0003325, +0.0006602 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0003325, +0.0006602 (5.1e-14); round 58: -0.0006302, +0.001251); near tie, x = 6.04, top 1.5 (8.33σ), s = 450.8: -0.03894, -0.03883 against -0.03934, -0.03923 (1.0%; κ = 0.514; ε_b 0.010, 0.010; the (D, y) system -0.03894, -0.03883 (3.8e-08); round 58: -0.07661, -0.07638); near tie, x = 6.04, top 1.5 (8.33σ), s = 5013: -0.0001463, +0.0002905 against -0.0001463, +0.0002905 (0.0%; κ = 0.512; ε_b -0.000, 0.000; the (D, y) system -0.0001463, +0.0002905 (6.1e-13); round 58: -0.0002857, +0.0005672); near tie, x = 3, right edge, s = 449.5: -0.4068, -0.4055 against -0.4325, -0.4311 (6.3%; κ = 0.610; ε_b 0.058, 0.058; the (D, y) system -0.4068, -0.4055 (9.6e-16); round 58: -0.7094, -0.707); near tie, x = 3, right edge, s = 5003: -0.001461, +0.002901 against -0.001461, +0.002901 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.001461, +0.002901 (4.2e-10); round 58: -0.002397, +0.00476); near tie, x = 4.9, right edge, s = 449.8: -0.1067, -0.1064 against -0.1087, -0.1084 (1.9%; κ = 0.528; ε_b 0.019, 0.019; the (D, y) system -0.1067, -0.1064 (1.6e-13); round 58: -0.2058, -0.2052); near tie, x = 4.9, right edge, s = 5004: -0.0003684, +0.0007316 against -0.0003684, +0.0007315 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0003684, +0.0007316 (6.4e-13); round 58: -0.0006983, +0.001387); near tie, x = 3, far edges, s = 450: -0.02268, -0.02262 against -0.02496, -0.0249 (10.1%; κ = 0.620; ε_b 0.088, 0.088; the (D, y) system -0.02268, -0.02262 (3.8e-09); round 58: -0.04028, -0.04019); near tie, x = 3, far edges, s = 5013: -0.0005367, +0.0008944 against -0.0005367, +0.0008943 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.0005367, +0.0008944 (2.2e-14); round 58: -0.0008807, +0.001467); near tie, x = 4.5, far edges, s = 450.5: -0.009054, -0.009033 against -0.009454, -0.009433 (4.4%; κ = 0.548; ε_b 0.042, 0.042; the (D, y) system -0.009054, -0.009033 (2.0e-06); round 58: -0.01725, -0.01721); near tie, x = 4.5, far edges, s = 5014: -0.0001824, +0.0003002 against -0.0001824, +0.0003002 (0.0%; κ = 0.537; ε_b -0.000, 0.000; the (D, y) system -0.0001824, +0.0003002 (2.1e-13); round 58: -0.0003398, +0.0005592); near tie, x = 4.9, far edges, s = 450.6: -0.007244, -0.007227 against -0.007503, -0.007486 (3.6%; κ = 0.539; ε_b 0.034, 0.034; the (D, y) system -0.007244, -0.007227 (6.5e-06); round 58: -0.01393, -0.0139); near tie, x = 4.9, far edges, s = 5014: -0.0001365, +0.0002243 against -0.0001365, +0.0002243 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0001365, +0.0002243 (2.6e-13); round 58: -0.0002588, +0.0004252); near tie, x = 5.07, far edges, s = 450.7: -0.006607, -0.006592 against -0.006824, -0.006808 (3.3%; κ = 0.535; ε_b 0.032, 0.032; the (D, y) system -0.006607, -0.006592 (1.0e-05); round 58: -0.01275, -0.01272); near tie, x = 5.07, far edges, s = 5014: -0.0001207, +0.0001983 against -0.0001208, +0.0001983 (0.0%; κ = 0.524; ε_b -0.000, 0.000; the (D, y) system -0.0001207, +0.0001983 (2.4e-13); round 58: -0.0002303, +0.0003781); near tie, x = 5.48, far edges, s = 450.7: -0.00533, -0.005318 against -0.005472, -0.005459 (2.7%; κ = 0.528; ε_b 0.026, 0.026; the (D, y) system -0.00533, -0.005318 (2.9e-05); round 58: -0.01036, -0.01033); near tie, x = 5.48, far edges, s = 5014: -8.984e-05, +0.0001473 against -8.984e-05, +0.0001473 (0.0%; κ = 0.518; ε_b -0.000, 0.000; the (D, y) system -8.984e-05, +0.0001473 (2.8e-13); round 58: -0.0001734, +0.0002844); near tie, x = 6.04, far edges, s = 450.8: -0.004036, -0.004027 against -0.004118, -0.004109 (2.0%; κ = 0.522; ε_b 0.020, 0.020; the (D, y) system -0.004036, -0.004027 (1.0e-04); round 58: -0.007897, -0.007879); near tie, x = 6.04, far edges, s = 5014: -6.02e-05, +9.861e-05 against -6.02e-05, +9.86e-05 (0.0%; κ = 0.512; ε_b -0.000, 0.000; the (D, y) system -6.02e-05, +9.861e-05 (1.6e-12); round 58: -0.0001175, +0.0001925); near tie, x = 7, far edges, s = 450.9: -0.002604, -0.002598 against -0.00264, -0.002634 (1.4%; κ = 0.514; ε_b 0.013, 0.013; the (D, y) system -0.002605, -0.002599 (4.7e-04); round 58: -0.005137, -0.005125); near tie, x = 7, far edges, s = 5014: -3.08e-05, +5.04e-05 against -3.08e-05, +5.04e-05 (0.0%; κ = 0.506; ε_b -0.000, 0.000; the (D, y) system -3.08e-05, +5.04e-05 (1.9e-08); round 58: -6.085e-05, +9.956e-05); near tie mirrored, x = 4.9, top 1 (5.00σ), s = 448.8: -0.1352, -0.1347 against -0.1373, -0.1368 (1.6%; κ = 0.528; ε_b 0.016, 0.016; the (D, y) system -0.1352, -0.1347 (6.4e-15); round 58: -0.26, -0.2591); near tie mirrored, x = 4.9, top 1 (5.00σ), s = 4994: -0.0004137, +0.0008215 against -0.0004137, +0.0008215 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0004137, +0.0008215 (3.1e-13); round 58: -0.0007842, +0.001557); near tie mirrored, x = 3, right edge, s = 449.5: -0.4069, -0.4056 against -0.4327, -0.4313 (6.3%; κ = 0.610; ε_b 0.058, 0.058; the (D, y) system -0.4069, -0.4056 (1.4e-15); round 58: -0.7095, -0.7072); near tie mirrored, x = 3, right edge, s = 5003: -0.001461, +0.002902 against -0.001461, +0.002902 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.001461, +0.002902 (4.4e-10); round 58: -0.002398, +0.004761); near tie mirrored, x = 4.9, right edge, s = 449.7: -0.1068, -0.1064 against -0.1088, -0.1084 (1.9%; κ = 0.528; ε_b 0.019, 0.019; the (D, y) system -0.1068, -0.1064 (1.5e-13); round 58: -0.2059, -0.2052); near tie mirrored, x = 4.9, right edge, s = 5004: -0.0003686, +0.0007318 against -0.0003686, +0.0007318 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0003686, +0.0007318 (7.2e-13); round 58: -0.0006986, +0.001387); near tie mirrored, x = 3, far edges, s = 450: -0.02268, -0.02263 against -0.02497, -0.02491 (10.1%; κ = 0.620; ε_b 0.088, 0.088; the (D, y) system -0.02268, -0.02263 (3.8e-09); round 58: -0.04029, -0.0402); near tie mirrored, x = 3, far edges, s = 5013: -0.0005369, +0.0008947 against -0.0005369, +0.0008946 (0.0%; κ = 0.609; ε_b -0.000, 0.000; the (D, y) system -0.0005369, +0.0008947 (9.2e-14); round 58: -0.000881, +0.001468); near tie mirrored, x = 4.5, far edges, s = 450.5: -0.009056, -0.009035 against -0.009456, -0.009434 (4.4%; κ = 0.548; ε_b 0.042, 0.042; the (D, y) system -0.009056, -0.009035 (1.1e-06); round 58: -0.01725, -0.01721); near tie mirrored, x = 4.5, far edges, s = 5014: -0.0001825, +0.0003004 against -0.0001825, +0.0003003 (0.0%; κ = 0.537; ε_b -0.000, 0.000; the (D, y) system -0.0001825, +0.0003004 (1.1e-13); round 58: -0.0003399, +0.0005594); near tie mirrored, x = 4.9, far edges, s = 450.6: -0.007245, -0.007228 against -0.007504, -0.007487 (3.6%; κ = 0.539; ε_b 0.034, 0.034; the (D, y) system -0.007245, -0.007228 (1.2e-06); round 58: -0.01393, -0.0139); near tie mirrored, x = 4.9, far edges, s = 5014: -0.0001366, +0.0002244 against -0.0001366, +0.0002244 (0.0%; κ = 0.528; ε_b -0.000, 0.000; the (D, y) system -0.0001366, +0.0002244 (2.0e-13); round 58: -0.0002589, +0.0004254); near tie mirrored, x = 5.07, far edges, s = 450.6: -0.006608, -0.006593 against -0.006825, -0.006809 (3.3%; κ = 0.535; ε_b 0.032, 0.032; the (D, y) system -0.006608, -0.006593 (2.9e-06); round 58: -0.01275, -0.01272); near tie mirrored, x = 5.07, far edges, s = 5014: -0.0001208, +0.0001983 against -0.0001208, +0.0001983 (0.0%; κ = 0.524; ε_b -0.000, 0.000; the (D, y) system -0.0001208, +0.0001983 (5.7e-14); round 58: -0.0002304, +0.0003782); near tie mirrored, x = 5.48, far edges, s = 450.6: -0.005331, -0.005319 against -0.005473, -0.00546 (2.7%; κ = 0.528; ε_b 0.026, 0.026; the (D, y) system -0.005331, -0.005319 (1.6e-05); round 58: -0.01036, -0.01033); near tie mirrored, x = 5.48, far edges, s = 5014: -8.988e-05, +0.0001474 against -8.988e-05, +0.0001474 (0.0%; κ = 0.518; ε_b -0.000, 0.000; the (D, y) system -8.988e-05, +0.0001474 (1.3e-13); round 58: -0.0001735, +0.0002845); near tie mirrored, x = 6.04, far edges, s = 450.7: -0.004037, -0.004028 against -0.004119, -0.00411 (2.0%; κ = 0.522; ε_b 0.020, 0.020; the (D, y) system -0.004037, -0.004028 (7.6e-05); round 58: -0.007898, -0.00788); near tie mirrored, x = 6.04, far edges, s = 5014: -6.022e-05, +9.864e-05 against -6.022e-05, +9.864e-05 (0.0%; κ = 0.512; ε_b -0.000, 0.000; the (D, y) system -6.022e-05, +9.864e-05 (3.3e-13); round 58: -0.0001176, +0.0001926); near tie mirrored, x = 7, far edges, s = 450.8: -0.002605, -0.002599 against -0.00264, -0.002634 (1.4%; κ = 0.514; ε_b 0.013, 0.013; the (D, y) system -0.002606, -0.0026 (4.0e-04); round 58: -0.005138, -0.005126); near tie mirrored, x = 7, far edges, s = 5014: -3.082e-05, +5.042e-05 against -3.082e-05, +5.042e-05 (0.0%; κ = 0.506; ε_b -0.000, 0.000; the (D, y) system -3.082e-05, +5.042e-05 (1.9e-08); round 58: -6.087e-05, +9.96e-05).
Over these rows κ = 0.506–0.620 and the slaved law's error 0.0%–10.1% of the largest rate: the slaving holds while the centroid's modes are faster than the relative mode, its error growing as they near it. Referee round 72: the slaving drops the centroid's offset y, slaved at first order in D (y = −(ww_e(w − w_e)/(W(w² + w_e²)))R⁻¹QD, Q = ⟨r, ∇²_cφ⟩, R = Q − W⟨∇_cφ∇_cφᵀ⟩), and its first correction is λ ≈ κλ_D(1 − ε_b), ε_b = (w − w_e)²ww_e/(w² + w_e²)²·R⁻¹Q = (1 − κ)λ_D/λ_b (λ_b the centroid's rate): ε_b (the larger) 0.0%–8.8% over these rows against the slaved law's error 0.0%–10.1%; the (D, y) system on the block's moments, y unslaved, within 4.7e-04 of the Jacobian's rates (relative to the largest).
The top edge at 0.95 (4.67σ; N = 1640): ∂c̄*/∂W along y 3.37e-11σ against 1.14e-12σ; x = 3: D_y at s = 700 -1.14e-13 (the law -1.14e-13) against -4.26e-15, the side toward, s_d×N₀/N 2.346e+07 against 2.346e+07; x = 4.9: D_y at s = 700 -4.6e-13 (the law -4.6e-13) against -1.57e-14, the side toward, s_d×N₀/N 3.15e+05 against 3.191e+05; x = 6.04: D_y at s = 700 -1.08e-12 (the law -1.07e-12) against -8.23e-14, the side toward, s_d×N₀/N 7.676e+04 against 8.439e+04.
The top edge at 1.5 (8.33σ; N = 2091): ∂c̄*/∂W along y 1.95e-17σ against 1.14e-12σ; x = 3: D_y at s = 700 0 (the law -4.59e-20) against -4.26e-15, the side toward, s_d×N₀/N 2.347e+07 against 2.346e+07; x = 4.9: D_y at s = 700 -1.85e-16 (the law 4.08e-20) against -1.57e-14, the side far, s_d×N₀/N none against 3.191e+05; x = 6.04: D_y at s = 700 -6.9e-12 (the law -6.3e-12) against -8.23e-14, the side toward, s_d×N₀/N 8.651e+04 against 8.439e+04.
Long double against float64 (x = 4.9, top 1): D at s = 700 (+1.85e-16, -1.55e-14) against (+2.31e-16, -1.57e-14), the side toward against toward, s_d 3.192e+05 against 3.191e+05.
Long double against float64 (x = 4.9, top 1.5): D at s = 700 (+9.25e-17, -3.70e-16) against (-4.63e-17, -1.85e-16), the side far against far, s_d none against none.

## RC. The relative mode in the continuum: continued about a merged block, injected into the race (referee rounds 66, 72)

The continuum here is the far window's Jacobian with the rest's drift dropped (referee round 72: the drift as a forcing, below and RD).

On the window with every edge at least 11.67σ from the block (N = 5751; s in its time, in brackets times N₀/N: Part C's, N₀ = 1681) R's race (`rp_job`'s) to s₁ = 400; there the merged block — both kernels at the centroid c_b less its offset y = (∂F/∂c_b)⁻¹F(c_b) off the merged rest — runs on (`merged_rhs`), and the centres' perturbation (δc_e = y − (w/W)D, δc_r = y + (w_e/W)D) by the race's exact Jacobian's centre block about it, in midpoint matrix exponentials (4 per output, the norm in logarithms: no floor), until D = δc_r − δc_e, past its least, has regrown to 1e-09σ (and to twice its least); there it is injected into the float64 race, run to the departure (2σ) or, on the far side, to its rests (0.05σ beyond feature 0 away from feature 2). The least |D| (log₁₀ σ), its share along e and its s; the contraction's e-folds from s₁; the injection's s; the side; the climb's end (|D| = 0.01σ) and s_d:

| target | x | least \|D\|: log₁₀ σ (share along e), at s | e-folds | injected at s | side | climb's end | s_d |
|---|---|---|---|---|---|---|---|
| near tie | 3 | -98.72 (+0.98), 4000 | 216 | 2.225e+05 | toward | 2.383e+05 (6.965e+04) | 8.048e+07 (2.352e+07) |
| near tie | 4.5 | -35.98 (+0.86), 4086 | 77 | 2e+05 | toward | 2.483e+05 (7.257e+04) | 4.299e+06 (1.257e+06) |
| near tie | 4.9 | -27.62 (+0.81), 4174 | 58 | 1.837e+05 | toward | 2.506e+05 (7.325e+04) | 1.225e+06 (3.58e+05) |
| near tie | 5.07 | -24.72 (+0.78), 4219 | 52 | 1.76e+05 | toward | 2.516e+05 (7.353e+04) | 6.342e+05 (1.854e+05) |
| near tie | 5.48 | -19.00 (+0.72), 4264 | 39 | 1.532e+05 | toward | 2.537e+05 (7.417e+04) | 2.883e+05 (8.428e+04) |
| near tie | 6.04 | -13.43 (+0.66), 4402 | 27 | 1.044e+05 | toward | 2.561e+05 (7.487e+04) | 2.964e+05 (8.665e+04) |
| near tie | 7 | -7.74 (+0.58), 4594 | 15 | 2.398e+04 | toward | 2.56e+05 (7.483e+04) | 3.273e+05 (9.565e+04) |
| near tie mirrored | 3 | -98.75 (-0.98), 4000 | 216 | 2.225e+05 | far | 2.383e+05 (6.966e+04) | — |
| near tie mirrored | 4.5 | -35.99 (-0.86), 4086 | 77 | 2e+05 | far | 2.484e+05 (7.26e+04) | — |
| near tie mirrored | 4.9 | -27.63 (-0.81), 4174 | 58 | 1.837e+05 | far | 2.507e+05 (7.329e+04) | — |
| near tie mirrored | 5.07 | -24.73 (-0.78), 4219 | 52 | 1.76e+05 | far | 2.517e+05 (7.357e+04) | — |
| near tie mirrored | 5.48 | -19.00 (-0.72), 4264 | 39 | 1.532e+05 | far | 2.539e+05 (7.422e+04) | — |
| near tie mirrored | 6.04 | -13.43 (-0.66), 4402 | 27 | 1.044e+05 | far | 2.564e+05 (7.495e+04) | — |
| near tie mirrored | 7 | -7.75 (-0.58), 4594 | 15 | 2.398e+04 | far | 2.566e+05 (7.499e+04) | — |

The rest's drift as a forcing (referee round 72): the merged base is no solution of the race — its kernels' speeds go as their coefficients —, so a block whose rest moves with W forces the relative mode by κk(∂c̄*/∂W)dW/ds, k = 2(w − w_e)/W (δc_e by −w/W of it, δc_r by +w_e/W). The forced D, log₁₀|D| σ (its share along e), at the continuum's least and at its injection (the forced part's size there over the homogeneous one's, inj·σ: the linear response, its sign deciding the side wherever it is the larger), for three drifts: the float64 sums' (`PairSums.block` along the merged run), exact sums on these samples (their float64 coordinates: ∂c̄*/∂W = (+3.151e-16, -8.227e-17)σ per unit W) and on the ideal grid (the aliasing: (+1.249e-37, -6.958e-40)σ) — the exact drifts at s₁ and 2e+06 within 1.6e-09 relative —, against the continuum's least (the first column):

| target | x | least \|D\| (share along e) | forced at the least: float64 sums | exact (these samples) | ideal grid (the aliasing) | at the injection, the forced over the homogeneous (dex; share): float64 sums | exact | ideal grid |
|---|---|---|---|---|---|---|---|---|
| near tie | 3 | -98.72 (+0.98) | -17.72 (-0.95) | -17.85 (-0.85) | -39.32 (-0.72) | +81.5 (-0.99) | +81.4 (-0.99) | +59.8 (-0.99) |
| near tie | 4.5 | -35.98 (+0.86) | -17.32 (-0.92) | -17.42 (-0.80) | -38.88 (-0.63) | +19.0 (-0.99) | +18.8 (-0.99) | -2.7 (-0.99) |
| near tie | 4.9 | -27.62 (+0.81) | -17.19 (-0.91) | -17.30 (-0.79) | -38.75 (-0.62) | +10.5 (-0.99) | +10.3 (-0.99) | -11.2 (-0.99) |
| near tie | 5.07 | -24.72 (+0.78) | -17.16 (-0.90) | -17.25 (-0.78) | -38.70 (-0.61) | +7.6 (-0.99) | +7.5 (-0.99) | -14.0 (-0.99) |
| near tie | 5.48 | -19.00 (+0.72) | -17.03 (-0.92) | -17.13 (-0.77) | -38.58 (-0.59) | +2.2 (-0.99) | +2.0 (-0.99) | -19.5 (-0.99) |
| near tie | 6.04 | -13.43 (+0.66) | -16.88 (-0.89) | -16.97 (-0.76) | -38.41 (-0.58) | -3.3 (-0.99) | -3.5 (-0.99) | -25.0 (-0.99) |
| near tie | 7 | -7.74 (+0.58) | -16.61 (-0.89) | -16.69 (-0.73) | -38.12 (-0.55) | -7.2 (-0.98) | -7.3 (-0.96) | -28.8 (-0.93) |
| near tie mirrored | 3 | -98.75 (-0.98) | -17.93 (+0.72) | -17.98 (+0.49) | -39.33 (+0.72) | +81.3 (+0.99) | +81.1 (+0.99) | +59.9 (+0.99) |
| near tie mirrored | 4.5 | -35.99 (-0.86) | -17.52 (+0.63) | -17.50 (+0.38) | -38.88 (+0.63) | +18.7 (+0.99) | +18.5 (+0.99) | -2.7 (+0.99) |
| near tie mirrored | 4.9 | -27.63 (-0.81) | -17.36 (+0.61) | -17.37 (+0.37) | -38.76 (+0.61) | +10.2 (+0.99) | +10.0 (+0.99) | -11.2 (+0.99) |
| near tie mirrored | 5.07 | -24.73 (-0.78) | -17.31 (+0.59) | -17.32 (+0.36) | -38.70 (+0.61) | +7.4 (+0.99) | +7.2 (+0.99) | -14.0 (+0.99) |
| near tie mirrored | 5.48 | -19.00 (-0.72) | -17.22 (+0.59) | -17.19 (+0.34) | -38.58 (+0.59) | +1.8 (+0.99) | +1.7 (+0.99) | -19.5 (+0.99) |
| near tie mirrored | 6.04 | -13.43 (-0.66) | -17.02 (+0.56) | -17.01 (+0.32) | -38.41 (+0.57) | -3.6 (+0.99) | -3.8 (+0.99) | -25.0 (+0.99) |
| near tie mirrored | 7 | -7.75 (-0.58) | -16.74 (+0.53) | -16.72 (+0.30) | -38.12 (+0.54) | -7.5 (+0.93) | -7.6 (+0.85) | -28.8 (+0.93) |

On the near tie the forced D exceeds the continuum's least at x = 3, 4.5, 4.9, 5.07, 5.48 for the samples' coordinates and at x = 3 for the aliasing.

Checks (the same gap; 16 substeps; y = 0 in the perturbation; the injection at another level): near tie, x = 4.9, 16 substeps, y, 1e-09σ: least -27.62 against -27.62, s_d 1.22494e+06 against 1.22494e+06, the side toward; near tie, x = 4.9, 4 substeps, y = 0, 1e-09σ: least -27.63 against -27.62, s_d 1.2251e+06 against 1.22494e+06, the side toward; near tie, x = 4.9, 4 substeps, y, 1e-07σ: least -27.62 against -27.62, s_d 1.22494e+06 against 1.22494e+06, the side toward; near tie, x = 4.9, 4 substeps, y, 1e-11σ: least -27.62 against -27.62, s_d 1.22489e+06 against 1.22494e+06, the side toward.

Where float64 resolves the offset itself (RP's least |D| on this window above the injection's level) the injected race's s_d against the plain run's: near tie, x = 7: 3.2725e+05 against 3.2727e+05.

s_d is the climb's end plus the hold after it: on Part C's window (RP, top edge at 1) the climb starts from the pin, in the continuum from the layer's offset. In Part C's time — the climb's end, s_d and their shifts, the continuum's less the pin's: x = 3: climb 9900 → 6.965e+04, s_d 2.346e+07 → 2.352e+07 (+0.3%; shifts 5.98e+04 and 5.98e+04); x = 4.5: climb 2.619e+04 → 7.257e+04, s_d 1.21e+06 → 1.257e+06 (+3.8%; shifts 4.64e+04 and 4.64e+04); x = 4.9: climb 3.426e+04 → 7.325e+04, s_d 3.191e+05 → 3.58e+05 (+12.2%; shifts 3.9e+04 and 3.9e+04); x = 5.07: climb 3.856e+04 → 7.353e+04, s_d 1.504e+05 → 1.854e+05 (+23.3%; shifts 3.5e+04 and 3.5e+04); x = 5.48: climb 5.09e+04 → 7.417e+04, s_d 6.101e+04 → 8.428e+04 (+38.1%; shifts 2.33e+04 and 2.33e+04); x = 6.04: climb 7.26e+04 → 7.487e+04, s_d 8.439e+04 → 8.665e+04 (+2.7%; shifts 2.26e+03 and 2.26e+03); x = 7: climb 7.483e+04 → 7.483e+04, s_d 9.566e+04 → 9.565e+04 (-0.0%; shifts -4.04 and -0.0733).

The split is a convention, the shift an identity (referee round 72): s_d less the climb's end at |D| = θσ, θ = 0.001, 0.01, 0.1, on Part C's window and in the continuum (Part C's time), and their relative difference: x = 3: 2.346e+07 and 2.346e+07 (+7.5e-12), 2.345e+07 and 2.345e+07 (+7.5e-12), 2.345e+07 and 2.345e+07 (+7.5e-12); x = 4.5: 1.186e+06 and 1.186e+06 (+5.2e-07), 1.184e+06 and 1.184e+06 (+5.2e-07), 1.182e+06 and 1.182e+06 (+5.2e-07); x = 4.9: 2.876e+05 and 2.876e+05 (+5.1e-06), 2.848e+05 and 2.848e+05 (+5.1e-06), 2.821e+05 and 2.821e+05 (+5.2e-06); x = 5.07: 1.15e+05 and 1.15e+05 (+1.4e-05), 1.118e+05 and 1.118e+05 (+1.5e-05), 1.088e+05 and 1.088e+05 (+1.5e-05); x = 5.48: 1.435e+04 and 1.435e+04 (+2.0e-07), 1.011e+04 and 1.011e+04 (+3.0e-06), 6033 and 6033 (+9.3e-06); x = 6.04: 1.812e+04 and 1.812e+04 (-2.5e-05), 1.179e+04 and 1.179e+04 (-2.2e-05), 5727 and 5727 (-1.3e-05); x = 7: 3.321e+04 and 3.321e+04 (-2.0e-04), 2.083e+04 and 2.082e+04 (-1.9e-04), 9031 and 9029 (-1.7e-04).

The runner-up's coefficient along the injected race toward feature 2 (referee round 72), relative to its value at the injection: at the climb's ends θ = 0.001, 0.01, 0.1, its peak (|D| there; s in Part C's time) and where it is 1% below: near tie, x = 3: +1.38e-08, +1.37e-06, +1.47e-04; +9.64e-04 at |D| = 0.208σ (s = 7.19e+04); 6.556e+05; near tie, x = 4.5: +8.47e-09, +8.46e-07, +8.37e-05; +1.34e-03 at |D| = 0.397σ (s = 7.713e+04); 1.32e+05; near tie, x = 4.9: +7.58e-09, +7.56e-07, +7.46e-05; +1.56e-03 at |D| = 0.491σ (s = 7.898e+04); 1.009e+05; near tie, x = 5.07: +7.31e-09, +7.30e-07, +7.19e-05; +1.69e-03 at |D| = 0.539σ (s = 7.983e+04); 9.314e+04; near tie, x = 5.48: +6.98e-09, +6.97e-07, +6.84e-05; +2.18e-03 at |D| = 0.715σ (s = 8.234e+04); none; near tie, x = 6.04: +7.30e-09, +7.29e-07, +7.13e-05; +7.01e-03 at |D| = 2.008σ (s = 8.665e+04); none; near tie, x = 7: +1.01e-08, +1.01e-06, +9.81e-05; +1.73e-02 at |D| = 2.004σ (s = 9.565e+04); none.

## RD. What moves the block's rest with W without an edge (referee round 72)

RC's merged block (the gap x on the window x ∈ [-2, 2], y ∈ [-1.5, 2]) at s = 400, 4000, 40000: the rest's shift ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩ (σ per unit W; in the continuum 0, ⟨φ, ∇_cφ⟩ = ½∇_c‖φ‖²) with the sums in float64, in the long double and exact (80 digits) on the committed samples (spacing 0.05, `np.arange`'s float64 coordinates, off the ideal grid by up to 3.73e-15), exact on the ideal grid against the aliasing's first harmonic (Poisson summation: ⟨φ, ∂_iφ⟩ = −⟨φ²⟩(2πq/ds)sin(2π(c_i − x₀,ᵢ)/ds), q = e^{−π²σ²/ds²} = 2.650e-39), and the same on a dyadic grid (spacing 0.046875, its float64 coordinates exact: off by 0; q = 1.283e-44):

| target | x | s | W | float64 | long double | exact | ideal grid, exact | the aliasing | dyadic: float64 | long double | exact | the aliasing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| near tie | 3 | 400 | 0.0261 | (+4.362e-16, -1.297e-16) | (+3.152e-16, -8.224e-17) | (+3.151e-16, -8.227e-17) | (+1.249e-37, -6.958e-40) | (+1.249e-37, -6.958e-40) | (+5.000e-17, +1.142e-17) | (-7.864e-20, -1.632e-20) | (-1.015e-45, +5.603e-43) | (-1.015e-45, +5.603e-43) |
| near tie | 3 | 4000 | 0.7989 | (+2.761e-16, -7.081e-17) | (+3.152e-16, -8.234e-17) | (+3.151e-16, -8.227e-17) | (+1.249e-37, -6.958e-40) | (+1.249e-37, -6.958e-40) | (-1.931e-16, +1.583e-18) | (-5.577e-20, +3.509e-20) | (-1.015e-45, +5.603e-43) | (-1.015e-45, +5.603e-43) |
| near tie | 3 | 40000 | 0.8000 | (+2.762e-16, -2.138e-16) | (+3.152e-16, -8.227e-17) | (+3.151e-16, -8.227e-17) | (+1.249e-37, -6.958e-40) | (+1.249e-37, -6.958e-40) | (-1.931e-16, +2.960e-18) | (-5.574e-20, +2.107e-21) | (-1.015e-45, +5.603e-43) | (-1.015e-45, +5.603e-43) |
| near tie mirrored | 3 | 400 | 0.0261 | (+1.243e-16, +2.091e-18) | (+3.150e-16, -8.162e-17) | (+3.149e-16, -8.160e-17) | (+1.249e-37, -4.294e-40) | (+1.249e-37, -4.294e-40) | (-3.605e-17, -9.552e-17) | (-2.261e-20, +2.111e-20) | (+1.015e-45, +5.603e-43) | (+1.015e-45, +5.603e-43) |
| near tie mirrored | 3 | 4000 | 0.7989 | (+1.086e-16, -1.359e-17) | (+3.149e-16, -8.156e-17) | (+3.149e-16, -8.160e-17) | (+1.249e-37, -4.294e-40) | (+1.249e-37, -4.294e-40) | (-1.094e-16, -8.275e-17) | (-4.845e-20, +1.852e-20) | (+1.015e-45, +5.603e-43) | (+1.015e-45, +5.603e-43) |
| near tie mirrored | 3 | 40000 | 0.8000 | (+1.086e-16, +5.992e-18) | (+3.149e-16, -8.160e-17) | (+3.149e-16, -8.160e-17) | (+1.249e-37, -4.294e-40) | (+1.249e-37, -4.294e-40) | (-1.093e-16, +6.002e-17) | (-4.837e-20, +9.429e-20) | (+1.015e-45, +5.603e-43) | (+1.015e-45, +5.603e-43) |

On the committed samples the long double's drift is the exact sums' to 3.2e-04 relative — the drift is the samples' coordinates', in any arithmetic — and float64's differs from it by up to 2.1e-16σ, its sums' rounding; on the ideal grid the exact drift is the aliasing's to 2.4e-14, on the dyadic grid to 4.3e-14, where the long double's (|·| ≤ 9.4e-20σ) and float64's (≤ 1.9e-16σ) are their sums' rounding.

## T3. The product rule at a vertex of hℤ³

The 8 kernels at the corners of the cube {0, H}³ (samples on a tensor grid of spacing 0.05), the feature (A = 0.8) at the cube's centre — each edge's axis at the transverse offset H/√2 ((n − 1)^{1/2}H/2; H/2 only in the plane) —; the corner 1_S behind the entrant by Σ_{i∈S}x_i: the 8-kernel inner problem against the product of the edges' 2-kernel shares at the offset H/√2 (the planar edge problem: in the continuum it depends on |t| only) and at H/2 (the planar offset); the last triples put a corner between x_θ(H/√2) and x_θ(H/2), where only the true offset's product gets the count:

| gaps (x₁, x₂, x₃) | 8-kernel: formed (stop) | the product at H/√2: formed, max \|Δ share\|, the formed corners' relative errors | the unformed corners with shares ≥ 5·10⁻⁴: the product at H/√2's relative errors (corners) | at H/2: formed, max \|Δ share\|, the formed corners' relative errors |
|---|---|---|---|---|
| (1, 1, 1) | 8 (converged) | 8, 8.68e-05, -0.00%, +0.03%, +0.03%, -0.01%, +0.03%, -0.01%, -0.01%, -0.29% | — | 8, 1.03e-02, -3.15%, +0.23%, +0.23%, +3.65%, +0.23%, +3.65%, +3.65%, +6.92% |
| (2, 2, 2) | 7 (converged) | 7, 3.58e-04, +0.06%, +0.09%, +0.09%, -0.47%, +0.09%, -0.47%, -0.47% | -6.5% (1) | 7, 1.95e-02, -3.39%, +3.41%, +3.41%, +10.04%, +3.41%, +10.04%, +10.04% |
| (1, 2, 3) | 6 (converged) | 6, 3.54e-04, +0.07%, -0.01%, +0.09%, +0.11%, -0.54%, -0.01% | -7.2% … -2.4% (2) | 7, 1.58e-02, -3.01%, +7.28%, +3.82%, +0.38%, +10.40%, +7.29% |
| (0.5, 3, 5) | 5 (converged) | 4, 5.56e-04, +0.10%, -1.43%, +0.03%, +0.15%, -0.11% | -35.9% … -3.3% (3) | 5, 9.59e-03, -1.78%, +14.31%, +8.65%, -0.06%, +10.35% |
| (3, 3, 3) | 4 (converged) | 4, 9.28e-04, +0.12%, +0.14%, +0.14%, +0.14% | -37.3% … -4.4% (4) | 4, 1.93e-02, -2.51%, +7.94%, +7.94%, +7.94% |
| (2, 4, 6) | 3 (cap) | 3, 8.79e-04, +0.11%, +0.04%, +0.15% | -50.2% … -1.1% (4) | 3, 1.42e-02, -1.79%, +12.34%, +5.14% |
| (4, 4, 4) | 4 (converged) | 4, 1.16e-03, +0.13%, +0.13%, +0.13%, +0.13% | -20.2% … -20.2% (3) | 4, 1.38e-02, -1.56%, +12.68%, +12.68%, +12.68% |
| (5.9, 9, 9) | 1 (cap) | 1, 2.20e-04, +0.02% | -0.2% … +0.1% (3) | 2, 2.24e-03, -0.23% |
| (5.95, 9, 9) | 1 (clock) | 1, 2.24e-04, +0.02% | -0.2% … +0.1% (3) | 2, 2.19e-03, -0.22% |
| (6, 9, 9) | 1 (cap) | 1, 2.15e-04, +0.02% | -0.2% … +0.1% (3) | 2, 2.13e-03, -0.22% |

## D2. The loop's threshold at a finite mass

The edge loop (Part A's, A = 0.8) with the feature at the offset t along the face: the gap at which its runner-up's share crosses 0.01 (bisection in [4.5, 8]) against the inner problem's x_θ(t) (T); the shift is the finite mass's:

| ε | t/H | the loop's crossing | x_θ(t) | shift |
|---|---|---|---|---|
| 1e-05 | 0.000 | 6.708 | 6.300 | +0.407 |
| 1e-05 | 0.375 | 6.552 | 6.149 | +0.403 |
| 1e-06 | 0.000 | 6.430 | 6.300 | +0.130 |
| 1e-06 | 0.375 | 6.277 | 6.149 | +0.127 |
| 1e-08 | 0.000 | 6.313 | 6.300 | +0.013 |
| 1e-08 | 0.375 | 6.161 | 6.149 | +0.011 |

## W. Where the selection is the Voronoi cell

W1. The samples' window (Part C's grid: the samples on [−1, 1]², the lattice up to ±1). The faces between neighbours on the axis x₂ = 0 (where their correlations with a feature are equal on the samples; the continuum's: the midpoints):

| pair | the midpoint | the samples' face |
|---|---|---|
| 0 \| 0.125 | 0.0625 | 0.0625 |
| 0.125 \| 0.25 | 0.1875 | 0.1875 |
| 0.25 \| 0.375 | 0.3125 | 0.3125 |
| 0.375 \| 0.5 | 0.4375 | 0.4375 |
| 0.5 \| 0.625 | 0.5625 | 0.5625 |
| 0.625 \| 0.75 | 0.6875 | 0.6881 |
| 0.75 \| 0.875 | 0.8125 | 0.8260 |
| 0.875 \| 1 | 0.9375 | 1.1028 |

The largest error of the nearest kernel's neighbours' margins against the continuum's, features on the axis at x₁: 0: 3.3e-16, 0.3: 3.0e-11, 0.5: 4.3e-06, 0.55: 1.4e-05, 0.6: 3.1e-04, 0.65: 7.9e-04, 0.7: 7.2e-03, 0.75: 1.4e-02, 0.8: 2.8e-02, 0.85: 9.2e-02, 0.9: 1.4e-01, 0.95: 2.5e-01; it first exceeds ε^{1/2}x̄_θ (6.208) at x₁ = 0.64 for ε = 1e-08 (2.40σ from the window's edge), 0.57 for ε = 1e-10 (2.87σ from the window's edge); the samples' entrant is not the nearest kernel at 6 of 99 positions, from x₁ = 0.82.

W2. Unequal masses: the kernel at (H, 0) of mass βε, the rest ε (the local lattice, ε = 1e-08); the entrant is the kernel of the least entry time ρ_j/v_j (ρ = b^{−1/2}): a power diagram, the face between (0, 0) and (H, 0) at H/2 − σ²ln β/H in the continuum. Features (A = 0.8) on the axis at the Voronoi face H/2, at the power face p (the samples') and 2 layer widths (w = 2σ²x̄_θε^{1/2}/H = 2.23e-04) to either side:

| β | p: samples, continuum | x₁ | entrant: power, Voronoi | runner-up, its gap | formed | shares of (0, 0), (H, 0) |
|---|---|---|---|---|---|---|
| 1 | 0.06250, 0.06250 | 0.06205 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9997, 0.0003 |
| 1 | 0.06250, 0.06250 | 0.06250 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 1.11e-12 | (+0, +0)H, (+1, +0)H | 0.5000, 0.5000 |
| 1 | 0.06250, 0.06250 | 0.06295 | (+1, +0)H, (+1, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.001 | 0.06232, 0.06232 | 0.06187 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9997, 0.0003 |
| 1.001 | 0.06232, 0.06232 | 0.06232 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 2.22e-12 | (+0, +0)H, (+1, +0)H | 0.4978, 0.5022 |
| 1.001 | 0.06232, 0.06232 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 5 | (+0, +0)H, (+1, +0)H | 0.0251, 0.9749 |
| 1.001 | 0.06232, 0.06232 | 0.06277 | (+1, +0)H, (+1, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.01 | 0.06071, 0.06071 | 0.06026 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9996, 0.0004 |
| 1.01 | 0.06071, 0.06071 | 0.06071 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 0 | (+0, +0)H, (+1, +0)H | 0.4782, 0.5218 |
| 1.01 | 0.06071, 0.06071 | 0.06116 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.01 | 0.06071, 0.06071 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 49.6 | (+1, +0)H | 0.0000, 1.0000 |
| 1.5 | -0.01048, -0.01048 | -0.01093 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9919, 0.0081 |
| 1.5 | -0.01048, -0.01048 | -0.01048 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 5.55e-12 | (+0, +0)H, (+1, +0)H | 0.0205, 0.9795 |
| 1.5 | -0.01048, -0.01048 | -0.01004 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0001, 0.9999 |
| 1.5 | -0.01048, -0.01048 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 1.84e+03 | (+1, +0)H | 0.0000, 1.0000 |

W3. Coarse samples along x₂ (spacing ds, the grid offset by 0.3ds; along x₁ the spacing 0.05): the face between (0, 0) and (0, H) for features at (0, x₂), its shift from H/2 on the samples and from the Poisson sum (harmonics k ≤ 3 of q = e^{−π²σ²/ds²}), the margin's slope there against the continuum's H/(2σ²), and the shift in layer widths at ε = 1e-08:

| ds | ds/σ | q | shift/H: samples, Poisson | slope ratio | shift/w |
|---|---|---|---|---|---|
| 0.05 | 0.333 | 2.650e-39 | -7.772e-16, +0.000e+00 | 1.0000 | 4.35e-13 |
| 0.1 | 0.667 | 2.269e-10 | +2.151e-09, +2.151e-09 | 1.0000 | 1.2e-06 |
| 0.15 | 1 | 5.172e-05 | +3.855e-04, +3.855e-04 | 0.9989 | 0.216 |
| 0.2 | 1.33 | 3.881e-03 | +3.131e-03, +3.131e-03 | 0.9276 | 1.75 |
| 0.3 | 2 | 8.480e-02 | -5.560e-01, -5.560e-01 | 0.7092 | 311 |

W4. Two kernels coinciding at the samples' centre of symmetry, the runner-up x_r behind, the feature (A = 0.8) δ off their centre: the inner problem's runner-up share (extrapolated in xi0) against the profile split's (iii)(b) — off the centre the kernels slide onto the feature at speeds ∝ their coefficients, apart on the way —, and the kernels' largest distance from the feature at the end:

| x_r | δ/σ | inner problem | profile split | difference | at the end (σ) |
|---|---|---|---|---|---|
| 2 | 0 | 0.197048 | 0.197048 | -1.53e-08 | 8.5e-13 |
| 2 | 0.07 | 0.196776 | 0.197048 | -2.72e-04 | 3.5e-13 |
| 2 | 0.2 | 0.194823 | 0.197048 | -2.22e-03 | 4.0e-13 |
| 2 | 0.4 | 0.188132 | 0.197048 | -8.92e-03 | 1.7e-12 |
| 4 | 0 | 0.056862 | 0.056862 | -1.06e-08 | 8.9e-13 |
| 4 | 0.07 | 0.056678 | 0.056862 | -1.84e-04 | 1.8e-12 |
| 4 | 0.2 | 0.055370 | 0.056862 | -1.49e-03 | 3.3e-12 |
| 4 | 0.4 | 0.051028 | 0.056862 | -5.83e-03 | 1.6e-12 |

W5. The 2-fold profile split's runner-up share p(x) against the softmax 1/(1 + e^{κx}) at the formation threshold: p's crossing of 0.01 (x_θ of coinciding kernels at the feature's centre) and the softmax's, the softmax's relative error at p's crossing and at x = 8, and e^{−2α/a} at the runner-up's coefficient (its amplitude coordinate α; the exponential branch's correction, D‴.51(i)):

| A | x_θ: profile split, softmax | the softmax's error: at the crossing, at x = 8 | α/a | e^{−2α/a} |
|---|---|---|---|---|
| 0.3 | 6.688, 6.541 | -0.097, -0.221 | 0.886 | 0.170 |
| 0.8 | 6.579, 6.541 | -0.026, -0.091 | 1.204 | 0.090 |
| 1 | 6.568, 6.541 | -0.018, -0.071 | 1.287 | 0.076 |

Numbers are produced by `experiments/exp_transit_nd.py`; nothing is edited by hand.
