# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)

The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points (the samples' inner product, e_k the starting point of largest correlation with feature k — the nearest in the continuum —, E their set), its margin m = 1 − max_{j∉E}|η(c_j)| and the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below 1e-08, no exits), exact with other kernels (three kernels, the residual below 1e-08), or not exact. The paths by `sliding_iss_nd` (in one dimension on the 1-D samples, with the exact Jacobian), the block polished by its least-squares fit unless a kernel below 10⁻² of the block's largest coefficient would land beyond one kernel width of both its start and its correlation's ascent's end (a rejected polish; the bound does not check single-kernel blocks or the block's largest kernel); the skeletons recomputed from the targets.

## A. One dimension

80 targets per minimal separation (seed 551; σ = 0.12, the grid's spacing 0.1); amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.25 | 80 | 26, 17, 0 | 6, 31, 0 | 0.343 | +0.0738 (0.60), -0.1442 (0.32), +0.1132 (0.26), +0.1949 (0.85), +0.0459 (0.56), +0.1369 (0.47), +0.1346 (0.54), -0.1486 (0.22), -0.0253 (0.13), +0.0880 (0.42), +0.1341 (0.38), +0.0814 (0.34), +0.0777 (0.25), +0.0566 (0.40), +0.1429 (0.42), +0.0832 (0.51), -0.0278 (0.08), -0.0286 (0.15), +0.0279 (0.50), -0.0327 (0.17), +0.1500 (0.61), +0.0897 (0.38), +0.0643 (0.18) |
| 0.3 | 80 | 35, 13, 0 | 6, 26, 0 | 0.174 | -0.0248 (0.11), +0.0351 (0.57), +0.0243 (0.06), +0.0619 (0.55), -0.0789 (0.10), -0.1248 (0.25), -0.0017 (0.03), +0.0343 (0.17), +0.0204 (0.15), +0.0301 (0.21), +0.0737 (0.39), +0.0104 (0.09), -0.0284 (0.09), -0.0583 (0.15), +0.0903 (0.37), +0.0365 (0.31), +0.0803 (0.45), +0.0429 (0.14), +0.1167 (0.31) |
| 0.35 | 80 | 37, 5, 0 | 9, 29, 0 | 0.151 | -0.0180 (0.10), -0.0094 (0.08), +0.0352 (0.12), -0.0078 (0.20), -0.0359 (0.11), -0.1760 (0.36), -0.0751 (0.19), -0.1770 (0.31), +0.0131 (0.11), +0.0233 (0.21), -0.0575 (0.13), -0.0328 (0.10), +0.1083 (0.24), +0.0487 (0.21) |
| 0.45 | 80 | 60, 4, 0 | 5, 11, 0 | 0.048 | +0.0150 (0.11), -0.0003 (0.06), +0.0074 (0.14), +0.0248 (0.09), -0.0077 (0.07), +0.0127 (0.08), -0.0015 (0.06), -0.0214 (0.13), -0.0015 (0.12) |
| 0.6 | 80 | 76, 2, 0 | 2, 0, 0 | 0.008 | -0.0005 (0.01), +0.0016 (0.01), -0.0007 (0.02), +0.0014 (0.01) |

Over all 400 targets the certificate's sign predicts exactness with E on 331; the disagreements' margins |m| ≤ 0.1949, at most 0.7608 times their cross-talk q. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 208 of 211; C = 0.8: 135 of 135; C = 1: 122 of 122. Every path ends at an exact fit with 3 kernels on 400 of 400 targets.

Exactness with E predicted by the sign of each margin (the relative margin min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|; the single-feature margin min_k(1 − max_{j≠e_k}⟨φ_kφ_j⟩/⟨φ_kφ_{e_k}⟩)), of 400: the certificate, 1 − max off E 331, the relative margin 344, the single feature (no signs or sums) 262, the atoms subtracted (no tangent directions) 390, the path's skeleton 397, the fitted skeleton 400, always E 262, E iff the first entrant is in E 327 (the first entrant off E — the outcome decided at the first entry — on 65; the fitted skeleton decided at its first phase on 101); the unpolished flow's outcome the polished path's on 400 of 400 (its blocks relaxed on 400).

## B. The amplitudes' range

40 configurations (seed 552; positions at least 0.35 apart, signs) with |A_k| = 10^{−r·u_k}, r = 0, 1, 2, 3 (the amplitudes' range up to 10^3); the certificate depends on the positions and signs only (21 of 40 with m > 0). The outcome is the same at every range on 31 of 40 configurations, while the order of the entries changes with the range on 35; it changes on 9: configuration 5 (m = -0.1052): other, other, other, nearest; configuration 7 (m = -0.1488): other, other, other, nearest; configuration 11 (m = -0.1783): other, nearest, nearest, nearest; configuration 12 (m = +0.0997): other, nearest, nearest, nearest; configuration 15 (m = -0.1138): other, nearest, nearest, nearest; configuration 23 (m = -0.0277): other, other, nearest, nearest; configuration 24 (m = -0.0007): other, nearest, nearest, nearest; configuration 28 (m = +0.0240): other, nearest, nearest, nearest; configuration 33 (m = -0.0828): other, other, nearest, nearest. The configurations whose outcome changes have |m|/q ≤ 0.67; of those with |m| > 0.8q, 6 of 6 keep their outcome at every range. The ranges realised, max|A|/min|A|, median [max]: r = 0: 1 [1], r = 1: 3.22 [9.08], r = 2: 10.4 [82.4], r = 3: 33.9 [749]; the amplitudes' order is the same for every r > 0 (u_k fixed), so this Part varies the range, not the order (D(c) varies the order).

Exactness with E predicted by the margins' signs, of 160: the certificate, 1 − max off E 119, the relative margin 107, the single feature (no signs or sums) 113, the atoms subtracted (no tangent directions) 140, the path's skeleton 160, the fitted skeleton 160, always E 113, E iff the first entrant is in E 125 (the first entrant off E — the outcome decided at the first entry — on 12; the fitted skeleton decided at its first phase on 34).

## C. The plane

40 targets per minimal separation (seed 553; σ = 0.15, the lattice's spacing 0.125; centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.3 | 40 | 14, 12, 0 | 5, 9, 0 | 0.215 | +0.0380 (0.31), -0.0111 (0.11), +0.0573 (0.22), +0.0736 (0.39), -0.0524 (0.29), +0.0185 (0.61), -0.0199 (0.16), +0.1001 (0.41), -0.0226 (0.24), +0.1289 (0.64), +0.0145 (0.68), +0.1268 (0.31), +0.1368 (0.52), -0.0225 (0.17), +0.0392 (0.46), +0.0192 (0.90), +0.0128 (0.19) |
| 0.45 | 40 | 26, 6, 0 | 4, 4, 0 | 0.107 | -0.0297 (0.10), -0.0907 (0.27), +0.0841 (0.48), +0.0124 (0.06), +0.0150 (0.10), -0.0317 (0.19), +0.0116 (0.13), +0.0142 (0.20), +0.0393 (0.30), -0.1624 (0.35) |
| 0.6 | 40 | 33, 3, 0 | 1, 3, 0 | 0.028 | +0.0189 (0.08), +0.0075 (0.07), +0.0064 (0.07), -0.0049 (0.07) |

Over all 120 targets the certificate's sign predicts exactness with E on 89. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 52 of 52; C = 0.8: 41 of 41; C = 1: 34 of 34. Every path ends at an exact fit with 3 kernels on 120 of 120 targets.

Exactness with E predicted by the margins' signs, of 120: the certificate, 1 − max off E 89, the relative margin 102, the single feature (no signs or sums) 83, the atoms subtracted (no tangent directions) 117, the path's skeleton 120, the fitted skeleton 120, always E 83, E iff the first entrant is in E 98 (the first entrant off E — the outcome decided at the first entry — on 15; the fitted skeleton decided at its first phase on 23).

## D. Designed families (referee round 43)

(a) Two features on starting points (1e-06 off) at −0.1 and −0.1 + D, the first's amplitude 1, the second's ±ρ_A: the margins (amplitude-free) and the outcome at each ρ_A (E exact with E, O exact with others, × not exact, t over 600 s):

| D | signs | m | m/q | relative | single-feature | 0.001 | 0.02 | 0.1 | 0.333 | 1 | 3 | 10 | 50 | 1e+03 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.3 | + + | -0.3400 | -1.622 | -0.1304 | +0.1594 | E | E | O | O | O | O | O | E | E |
| 0.3 | + − | +0.2216 | +1.057 | +0.0119 | +0.1594 | O | O | O | O | E | O | O | O | O |
| 0.4 | + + | -0.0502 | -0.808 | +0.0119 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.4 | + − | +0.1724 | +2.773 | +0.1102 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.5 | + + | +0.0972 | +7.458 | +0.1102 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.5 | + − | +0.1613 | +12.376 | +0.1483 | +0.1594 | E | E | E | E | E | E | E | E | E |

(b) Two features at s·σ, equal amplitudes, the midpoint uniform over a grid cell (20 per row): the outcomes and how often each margin's sign is right:

| s | signs | E, O, ×, t | m > 0 | the certificate right | relative right | single-feature right | atoms subtracted right | the skeleton right | the fitted skeleton right | always E right |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.5 | + + | 5, 15, 0, 0 | 0 | 15 | 20 | 5 | 20 | 20 | 20 | 5 |
| 1.5 | + − | 0, 20, 0, 0 | 20 | 0 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2 | + − | 2, 18, 0, 0 | 20 | 2 | 20 | 2 | 20 | 20 | 20 | 2 |
| 2.5 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2.5 | + − | 1, 19, 0, 0 | 20 | 1 | 20 | 1 | 20 | 20 | 20 | 1 |
| 3 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 3 | + − | 15, 5, 0, 0 | 20 | 15 | 20 | 15 | 20 | 20 | 20 | 15 |
| 3.5 | + + | 12, 8, 0, 0 | 0 | 8 | 19 | 12 | 20 | 20 | 20 | 12 |
| 3.5 | + − | 16, 4, 0, 0 | 20 | 16 | 20 | 16 | 19 | 20 | 20 | 16 |
| 4 | + + | 17, 3, 0, 0 | 15 | 18 | 20 | 17 | 20 | 20 | 20 | 17 |
| 4 | + − | 16, 4, 0, 0 | 20 | 16 | 20 | 16 | 20 | 20 | 20 | 16 |

(c) 30 of A's targets with |m| > 0.8q, the amplitude magnitudes 1, 10^{−1/2}, 10^{−1} in all six orders (the signs kept): the outcome changes with the order on 0 of them.

Over D, exactness with E predicted by the margins' signs: (a) of 54: the certificate, 1 − max off E 33, the relative margin 42, the single feature (no signs or sums) 41, the atoms subtracted (no tangent directions) 34, the path's skeleton 52, the fitted skeleton 54, always E 41, E iff the first entrant is in E 42 (the first entrant off E — the outcome decided at the first entry — on 1; the fitted skeleton decided at its first phase on 1); (b) of 240: the certificate, 1 − max off E 151, the relative margin 239, the single feature (no signs or sums) 84, the atoms subtracted (no tangent directions) 239, the path's skeleton 240, the fitted skeleton 240, always E 84, E iff the first entrant is in E 205 (the first entrant off E — the outcome decided at the first entry — on 121; the fitted skeleton decided at its first phase on 125); the unpolished flow's outcome the polished path's on 240 of 240 (its blocks relaxed on 240); (c) of 180: the certificate, 1 − max off E 180, the relative margin 180, the single feature (no signs or sums) 180, the atoms subtracted (no tangent directions) 180, the path's skeleton 180, the fitted skeleton 180, always E 180, E iff the first entrant is in E 180 (the first entrant off E — the outcome decided at the first entry — on 0; the fitted skeleton decided at its first phase on 29). Paths over 600 s: 0; polishes rejected: 1 on 1 paths.

D(b)'s failures that are not the second feature's own near ties (m₀ > 0.05), per row: 1.5σ + + 11, 1.5σ + − 11, 2σ + + 5, 2σ + − 9, 2.5σ + + 16, 2.5σ + − 15, 3σ + + 9, 3σ + − 0, 3.5σ + + 0, 3.5σ + − 3 (m₀ = 0.054, 0.056, 0.055), 4σ + + 0, 4σ + − 0.

Outside the band, the certificate's sign on the designed targets: |m| > 0.5q: wrong on 47 of 399 (|m|/q up to 1.622, 34 of them with m > 0); |m| > 0.8q: wrong on 25 of 333 (|m|/q up to 1.622, 12 of them with m > 0); |m| > 1q: wrong on 12 of 271 (|m|/q up to 1.622, 8 of them with m > 0).

## E. The two-feature map

The first feature (amplitude 1) on the starting point −0.1 (1e-06 off), the second at the distance D = 0.12 … 0.6 (step 0.01; σ = 0.12) with the amplitude ±10^ℓ, ℓ = -3 … 3 (step 0.125): the fitted skeleton's outcome at every point (E: exact with the nearest kernels, ·: not), shown at every 0.06 in D and every ½ in ℓ:

Same signs:

| D/σ \ ℓ | -3 | -2.5 | -2 | -1.5 | -1 | -0.5 | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 1.50 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.00 | E | E | E | E | E | · | · | · | · | · | · | · | · |
| 2.50 | E | E | E | · | · | · | · | · | · | · | E | E | E |
| 3.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 3.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 5.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |

Opposite signs:

| D/σ \ ℓ | -3 | -2.5 | -2 | -1.5 | -1 | -0.5 | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 1.50 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.50 | · | · | · | · | · | · | E | · | · | · | · | · | · |
| 3.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 3.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 5.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |

The fitted skeleton's threshold D*(ℓ) on this grid (the first feature on a starting point, |ℓ| ≤ 3: not (iii)'s Δ₀) — every D of the grid at or beyond it exact with E —, in units of σ: same signs ℓ = -3: 4.67, ℓ = -2: 4.67, ℓ = -1: 4.67, ℓ = -0.5: 4.67, ℓ = 0: 4.67, ℓ = 0.5: 4.67, ℓ = 1: 4.67, ℓ = 2: 4.67, ℓ = 3: 4.67; opposite signs ℓ = -3: 3.75, ℓ = -2: 3.75, ℓ = -1: 3.75, ℓ = -0.5: 3.75, ℓ = 0: 3.75, ℓ = 0.5: 2.92, ℓ = 1: 2.75, ℓ = 2: 2.75, ℓ = 3: 2.83; over all ℓ at most 4.67σ. Beyond 3.25σ every failure is the second feature's own near tie (its single-feature margin m₀ ≤ 0.05, near a midpoint of the grid: D = 0.44, 0.45, 0.46, 0.55); without them D*(ℓ) is same signs ℓ = -3: 3.17, ℓ = -2: 3.25, ℓ = -1: 3.25, ℓ = -0.5: 3.25, ℓ = 0: 3.33, ℓ = 0.5: 3.17, ℓ = 1: 2.67, ℓ = 2: 2.50, ℓ = 3: 2.50; opposite signs ℓ = -3: 2.83, ℓ = -2: 2.83, ℓ = -1: 2.83, ℓ = -0.5: 2.83, ℓ = 0: 2.83, ℓ = 0.5: 2.67, ℓ = 1: 2.75, ℓ = 2: 2.75, ℓ = 3: 2.83, at most 3.33σ.

The two ends of the ratio ((v)): same signs: constant for ℓ ≤ −2 at 48 and for ℓ ≥ 2 at 49 of 49 separations, the two ends' outcomes different at 7; opposite signs: constant for ℓ ≤ −2 at 49 and for ℓ ≥ 2 at 48 of 49 separations, the two ends' outcomes different at 5.

The margins' signs against the fitted skeleton on the map's 4802 points: the certificate, 1 − max off E 3128, the relative margin 4227, the single feature (no signs or sums) 2754, the atoms subtracted (no tangent directions) 3991, the path's skeleton 4719; exact with E on 2754; the linear skeleton's |m_sk| > 0.01 where its sign differs from the fitted skeleton's at 33 points (|m_sk| up to 0.134).

(ii)'s order: |m_sk − m_fit| against the ratio of the smaller feature's amplitude to the larger's, its log–log slope over 10⁻³ … 10⁻¹ (1: first order in the tail, 2: second): 2.50σ + + second small 1.00, 2.50σ + + first small 1.00, 3.50σ + + second small 0.99, 3.50σ + + first small 0.99, 4.00σ + + second small 0.99, 4.00σ + + first small 0.98, 3.50σ + − second small 0.98, 4.00σ + − second small 0.98.

The path at 240 sampled points (|ℓ| ≤ 2.5; under 300 s each): it agrees with the fitted skeleton at 240 of 240 (0 over the limit); exact with E at 134; polishes rejected 14 on 3 paths.

## F. A wide range

Three features at −s, 0, s (+0.013; s = 0.6, 0.36: 5σ, 3σ), the amplitudes 1, 0.001, 1e-06 in all six orders, the signs + + + and + − +; exits only below 10⁻¹² of the largest coefficient: the outcome per order (E exact with E, O with others, × not, t over the limit), the polishes rejected, the skeletons' signs:

| s | signs | outcomes (orders 012, 021, 102, 120, 201, 210) | polishes rejected | the skeleton right | the fitted skeleton right |
|---|---|---|---|---|---|
| 5σ | + + + | EEEEEE | 0 | 6 of 6 | 6 of 6 |
| 5σ | + − + | EEEEEE | 0 | 6 of 6 | 6 of 6 |
| 3σ | + + + | O×OOOO | 0 | 6 of 6 | 6 of 6 |
| 3σ | + − + | OEOEOO | 0 | 6 of 6 | 6 of 6 |

## G. The polish's audit on every target of A–D (referee round 140)

The unpolished flow (polish off; the blocks relaxed by the fast flow alone) against the polished path, per family (A and D(b) from their Parts, the rest computed here): its outcome the polished path's, its blocks relaxed, and the polished path's outcome as recomputed here the stored one:

| family | targets | the flow's outcome the path's | blocks relaxed | the path recomputed as stored | polishes rejected (paths) |
|---|---|---|---|---|---|
| A | 400 of 400 | 400 | 400 | 400 | 2 (1) |
| B | 160 of 160 | 160 | 160 | 160 | 8 (4) |
| C | 120 of 120 | 120 | 120 | 120 | 0 (0) |
| D(a) | 54 of 54 | 54 | 54 | 54 | 0 (0) |
| D(b) | 240 of 240 | 240 | 240 | 240 | 1 (1) |
| D(c) | 180 of 180 | 180 | 180 | 180 | 0 (0) |

Over A–D: the unpolished flow's outcome the polished path's on 1154 of 1154 targets, its blocks relaxed on 1154.

## H. Baselines off the path: η_V, basis pursuit, continuous OMP (referee round 140)

Per family: the targets; the path's outcomes (exact with E, exact with other kernels — both fit the features exactly with 3 kernels, at their centres —, not exact); η_V nondegenerate (|η_V| < 1 on the grid of spacing 0.002 (one dimension) or 0.01 (the plane) over [-1.2, 1.2] or [-1, 1]² beyond 0.1σ of the features, and s_k∇²η_V(c*_k) negative definite); basis pursuit on the grid of spacing 0.005 (one dimension only; not computed in the plane), its points within half a spacing of a feature replaced by the feature, recovering the features (its least certificate's sup-norm t < 1); continuous OMP exact after 3 steps (the rms residual below 1e-08; its selection grid 0.001 or 0.01, refined), and with the joint slide after each step; and the targets on which the sign of η_V's margin on the starting points off E, ν_V = 1 − max_{j∉E}|η_V(c_j)|, and of the certificate's m predict exactness with E:

| family | targets | path: E, other, not | η_V nondegenerate | basis pursuit exact | OMP exact | OMP with slide exact | ν_V's sign right | m's sign right |
|---|---|---|---|---|---|---|---|---|
| A | 400 | 262, 138, 0 | 400 | 400 | 0 | 400 | 262 | 331 |
| B | 160 | 113, 47, 0 | 160 | 160 | 0 | 160 | 113 | 119 |
| C | 120 | 83, 37, 0 | 120 | — | 0 | 120 | 83 | 89 |
| D(a) | 54 | 41, 13, 0 | 54 | 54 | 0 | 54 | 41 | 33 |
| D(b) | 240 | 84, 156, 0 | 220 | 240 | 0 | 240 | 104 | 151 |
| D(c) | 180 | 180, 0, 0 | 180 | 180 | 0 | 180 | 180 | 180 |
| all | 1154 | 763, 391, 0 | 1134 | 1034 of 1034 (one dimension) | 0 | 1154 | 783 | 903 |

Continuous OMP without the slide: the largest distance of a selected atom from its nearest feature, median 0.0718σ, range 3.19e-07σ–1.66σ; its rms residual after 3 steps median 0.00999; the largest |η_V| beyond 1σ of the features: median 0.7889, range 0.5225–2.5131; basis pursuit's t (one dimension): median 0.9731, the linear program not solved on 0. Where the path is not exact with E but fits the features with other kernels, η_V is nondegenerate on 371 of 391, and basis pursuit exact on 354 of the 354 in one dimension.

## I. BLASSO at a finite regularisation (referee round 140)

The BLASSO on Part H's basis-pursuit grid (spacing 0.005 over [-1.2, 1.2], its points within half a spacing of a feature replaced by the feature), min_x ½⟨(Σ_jx_jφ(c_j) − f*)²⟩ + λ‖x‖₁ (the samples' mean, no noise), at λ = ℓλ_max, λ_max = max_j|⟨φ(c_j)f*⟩| the least λ at which x = 0 is optimal, ℓ = 0.1, 0.01, 0.001, 0.0001, on every one-dimensional target of A, B and D (1034); solved by feature-sign search (exact Newton steps on the support and signs with an exact line search over their zero crossings, warm started from the larger λ) to the KKT conditions: the largest violation over all 4136 solves 4.08e-12 relative to λ. Recovered: every nonzero coefficient within 2 spacings (0.01) of a feature, and each feature's cluster nonempty with its sum of A_k's sign; the support exactly the features; features lost (an empty cluster); the amplitude error max_k|Σ_cluster x − A_k|/|A_k| on the recovered targets, and that error over ℓ; and how recovery splits by the path's outcome (exact with E, exact with other kernels):

| family | ℓ | targets | recovered | support exactly the features | targets with a feature lost | amplitude error median [max] | error/ℓ median | recovered: path E, path other |
|---|---|---|---|---|---|---|---|---|
| A | 0.1 | 400 | 337 | 97 | 37 | 0.158 [0.254] | 1.58 | 249 of 262, 88 of 138 |
| A | 0.01 | 400 | 399 | 97 | 0 | 0.0165 [0.0518] | 1.65 | 262 of 262, 137 of 138 |
| A | 0.001 | 400 | 400 | 97 | 0 | 0.00165 [0.00569] | 1.65 | 262 of 262, 138 of 138 |
| A | 0.0001 | 400 | 400 | 97 | 0 | 1.65e-4 [5.69e-4] | 1.65 | 262 of 262, 138 of 138 |
| B | 0.1 | 160 | 87 | 8 | 66 | 0.166 [0.939] | 1.66 | 63 of 113, 24 of 47 |
| B | 0.01 | 160 | 142 | 10 | 16 | 0.0346 [0.619] | 3.46 | 102 of 113, 40 of 47 |
| B | 0.001 | 160 | 157 | 12 | 3 | 0.00397 [0.375] | 3.97 | 112 of 113, 45 of 47 |
| B | 0.0001 | 160 | 160 | 12 | 0 | 4.19e-4 [0.0749] | 4.19 | 113 of 113, 47 of 47 |
| D(a) | 0.1 | 54 | 12 | 0 | 39 | 0.298 [0.989] | 2.98 | 12 of 41, 0 of 13 |
| D(a) | 0.01 | 54 | 33 | 0 | 21 | 0.0484 [0.507] | 4.84 | 26 of 41, 7 of 13 |
| D(a) | 0.001 | 54 | 43 | 0 | 10 | 0.0096 [0.987] | 9.6 | 33 of 41, 10 of 13 |
| D(a) | 0.0001 | 54 | 52 | 0 | 2 | 0.00105 [0.109] | 10.5 | 41 of 41, 11 of 13 |
| D(b) | 0.1 | 240 | 136 | 0 | 80 | 0.102 [0.122] | 1.02 | 75 of 84, 61 of 156 |
| D(b) | 0.01 | 240 | 220 | 0 | 20 | 0.011 [0.0321] | 1.1 | 84 of 84, 136 of 156 |
| D(b) | 0.001 | 240 | 220 | 0 | 0 | 0.0011 [0.00322] | 1.1 | 84 of 84, 136 of 156 |
| D(b) | 0.0001 | 240 | 220 | 0 | 0 | 1.10e-4 [3.22e-4] | 1.1 | 84 of 84, 136 of 156 |
| D(c) | 0.1 | 180 | 96 | 62 | 84 | 0.999 [1] | 9.99 | 96 of 180, 0 of 0 |
| D(c) | 0.01 | 180 | 180 | 108 | 0 | 0.1 [0.103] | 10 | 180 of 180, 0 of 0 |
| D(c) | 0.001 | 180 | 180 | 108 | 0 | 0.01 [0.0103] | 10 | 180 of 180, 0 of 0 |
| D(c) | 0.0001 | 180 | 180 | 108 | 0 | 1.00e-3 [0.00103] | 10 | 180 of 180, 0 of 0 |
| all | 0.1 | 1034 | 668 | 167 | 306 | 0.153 [1] | 1.53 | 495 of 680, 173 of 354 |
| all | 0.01 | 1034 | 974 | 215 | 57 | 0.0173 [0.619] | 1.73 | 654 of 680, 320 of 354 |
| all | 0.001 | 1034 | 1000 | 217 | 13 | 0.00176 [0.987] | 1.76 | 671 of 680, 329 of 354 |
| all | 0.0001 | 1034 | 1012 | 217 | 2 | 1.77e-4 [0.109] | 1.77 | 680 of 680, 332 of 354 |

At ℓ = 0.0001: recovered on 1012 of the 1014 targets with η_V nondegenerate and on 0 of the 20 without; on 680 of the 680 the fitted skeleton predicts exact with E and on 332 of the 354 it does not; the support's size beyond the number of features (the grid neighbours of Duval–Peyré's extended support), median [max] per ℓ: 0.1: 2 [3], 0.01: 2 [5], 0.001: 2 [4], 0.0001: 2 [4]; feature-sign iterations per solve, median [max]: 2 [58].

## J. Basis pursuit in the plane by a certificate in the features' jets (referee round 140)

Part H's linear program for basis pursuit's least certificate has, in the plane, a p of dimension 1094 (the span of the 6561 grid columns on the 1681 samples, singular values below 10⁻¹³ of the largest dropped) and 6558 pairs of constraints; here p is restricted to the span of the features' second-order jets (φ, ∇_cφ, ∇²_cφ at the features: three dimensions per feature in one dimension, six in the plane), and the same linear program (HiGHS) is solved over every grid point (spacing 0.005 or 0.025, the features replacing their neighbours). Every such p is a certificate, so its t_J bounds Part H's t from above, and t_J < 1 proves that the features are the unique minimum-ℓ¹ fit on the grid; η_V lies in that span, and its largest |η_V| on the same points bounds t_J from above:

| family | targets | t_J < 1 (basis pursuit exact, certified) | t_J median [max] | max\|η_V\| on the grid median [max] | Part H's t median | t_J − t median [range] | path: E, other |
|---|---|---|---|---|---|---|---|
| A | 400 | 400 | 0.99857 [0.99909] | 0.99984 [0.99998] | 0.9740 | 0.0246 [0.0185, 0.1168] | 262, 138 |
| B | 160 | 160 | 0.99858 [0.99892] | 0.99983 [0.99992] | 0.9744 | 0.0242 [0.0204, 0.0619] | 113, 47 |
| C | 120 | 120 | 0.98044 [0.98990] | 0.99781 [0.99912] | — | — | 83, 37 |
| D(a) | 54 | 54 | 0.99839 [0.99854] | 0.99956 [0.99978] | 0.9363 | 0.0622 [0.0293, 0.1191] | 41, 13 |
| D(b) | 240 | 240 | 0.99850 [0.99952] | 0.99979 [2.51306] | 0.9705 | 0.0276 [0.0152, 0.1368] | 84, 156 |
| D(c) | 180 | 180 | 0.99856 [0.99882] | 0.99984 [0.99989] | 0.9740 | 0.0246 [0.0189, 0.1020] | 180, 0 |

In one dimension (1034 targets) t_J < 1 exactly where Part H's t < 1 on 1034, and t_J ≥ t − 10⁻⁹ on 1034. In the plane basis pursuit on the grid is certified exact on 120 of 120 targets, on 37 of the 37 where the path fits the features with other kernels than E; the full linear program in the plane (Part H's t) was not solved here (Part L solves it at three truncations of the span). Wall time per target (one process), median: one dimension 0.0091 s, the plane 0.42 s.

## K. Part I's failures by type (referee round 149)

Every stored solution of Part I classified, without a new solve: each nonzero belongs to its nearest feature, within the tolerance 0.01 to its cluster (Part I's rule), beyond it but within 0.5σ = 0.06 *displaced*, beyond 0.06 *spurious*. A feature is *lost* when no nonzero within 0.06 has it as its nearest feature. A target that Part I does not recover is counted in the first of lost, spurious, displaced, wrong sign (a cluster's sum not of A_k's sign) that applies; the four cover every failure. The last columns: Part I's column "a feature lost" (a cluster empty within the tolerance) restricted to the failures, the displaced nonzeros' largest distance, and the lost features' |A_k|/max|A|:

| family | ℓ | targets | not recovered | a feature lost | else spurious | else displaced | else wrong sign | a cluster empty within the tolerance | displaced: largest distance | lost: \|A_k\|/max\|A\| median [max] |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 0.1 | 400 | 63 | 0 | 0 | 63 | 0 | 37 | 0.0425 | — |
| A | 0.01 | 400 | 1 | 0 | 0 | 1 | 0 | 0 | 0.0101 | — |
| A | 0.001 | 400 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| A | 0.0001 | 400 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| B | 0.1 | 160 | 73 | 49 | 0 | 24 | 0 | 66 | 0.0354 | 0.0305 [0.0969] |
| B | 0.01 | 160 | 18 | 12 | 0 | 6 | 0 | 16 | 0.0354 | 0.0045 [0.00878] |
| B | 0.001 | 160 | 3 | 0 | 0 | 3 | 0 | 3 | 0.0304 | — |
| B | 0.0001 | 160 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| D(a) | 0.1 | 54 | 42 | 30 | 0 | 12 | 0 | 39 | 0.0500 | 0.02 [0.1] |
| D(a) | 0.01 | 54 | 21 | 12 | 0 | 9 | 0 | 21 | 0.0500 | 0.001 [0.001] |
| D(a) | 0.001 | 54 | 11 | 6 | 0 | 5 | 0 | 10 | 0.0450 | 0.001 [0.001] |
| D(a) | 0.0001 | 54 | 2 | 0 | 0 | 2 | 0 | 2 | 0.0150 | — |
| D(b) | 0.1 | 240 | 104 | 0 | 0 | 104 | 0 | 80 | 0.0440 | — |
| D(b) | 0.01 | 240 | 20 | 0 | 20 | 0 | 0 | 20 | 0.0173 | — |
| D(b) | 0.001 | 240 | 20 | 0 | 20 | 0 | 0 | 0 | — | — |
| D(b) | 0.0001 | 240 | 20 | 0 | 20 | 0 | 0 | 0 | — | — |
| D(c) | 0.1 | 180 | 84 | 80 | 0 | 4 | 0 | 84 | 0.0107 | 0.1 [0.1] |
| D(c) | 0.01 | 180 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| D(c) | 0.001 | 180 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| D(c) | 0.0001 | 180 | 0 | 0 | 0 | 0 | 0 | 0 | — | — |
| all | 0.1 | 1034 | 366 | 159 | 0 | 207 | 0 | 306 | 0.0500 | 0.0884 [0.1] |
| all | 0.01 | 1034 | 60 | 24 | 20 | 16 | 0 | 57 | 0.0500 | 0.00161 [0.00878] |
| all | 0.001 | 1034 | 34 | 6 | 20 | 8 | 0 | 13 | 0.0450 | 0.001 [0.001] |
| all | 0.0001 | 1034 | 22 | 0 | 20 | 2 | 0 | 2 | 0.0150 | — |

Over all failures the displaced nonzeros lie between 0.01 and 0.0500 of their features and the spurious ones at least 0.1977 from every feature, so the split does not depend on the radius between those two distances; unclassified failures: 0. The targets with a spurious nonzero at some ℓ: 20, of which η_V is degenerate on 20.

Part G's audit as entry sequences: the unpolished flow enters the same kernels in the same order as the polished path on 1034 of the 1034 one-dimensional targets (640 of A and D(b) from their own rows, the others from Part G's); in the plane, where the unpolished flow's sequence is not stored, its outcome class is the polished path's on 120 of 120.

## L. Basis pursuit's least certificate in the plane, at three truncations of the span (after round 149)

Part H's linear program in the plane, in full: p in the span of the grid's columns (spacing 0.025, 6561 points over the samples' square), every grid point farther than half a spacing from the features constrained. On the tensor samples these columns are G1 ⊗ G1 (G1 the one-dimensional Gaussians), so with G1 = U1S1V1ᵀ the span's orthonormal basis is U1 ⊗ U1 with the singular values s_is_j, truncated at s_is_j > c·s_0² (c = 10⁻¹³ is Part H's truncation; the dimension below); in the coordinates w = (s_is_j)(u_i ⊗ u_j)ᵀp the grid values are η = (V1 ⊗ V1)w. The program, min t subject to η(c*_k) = sign A_k and |η| ≤ t on the grid, is solved by a primal–dual interior-point method (Mehrotra's predictor–corrector) whose normal matrix (V1 ⊗ V1)ᵀD(V1 ⊗ V1) is formed through the tensor factors. Each iterate gives a certified bracket: t, the sup of |η| over the grid of a w that satisfies the equalities, from above, and t_low, the dual iterate made exactly feasible (its stationarity residual removed on the grid and its mass renormalised), from below by weak duality; the iteration stops when the bracket is below 1e-09 or has not narrowed for 8 iterations. t_J is Part J's jet certificate (an upper bound on t, up to the truncation: its p need not lie in the truncated span). |t_dir − t| is the certificate's p, evaluated directly on the samples, against t, and feq its largest |⟨pφ(c*_k)⟩ − sign A_k|: the rounding of p's large components in the directions of small s_is_j.

| truncation c | dimension of p | targets | t < 1 | t median [min, max] | t ≤ t_J | t_J − t median [min, max] | bracket t − t_low max | iterations median [max] | \|t_dir − t\| max | feq max | seconds median |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1e-13 | 1094 | 120 | 120 | 0.62010 [0.56686, 0.78466] | 120 | 0.3601 [0.2050, 0.4087] | 2.30e-6 | 39 [55] | 0.00126 | 3.44e-4 | 2.8 |
| 1e-10 | 839 | 120 | 120 | 0.69615 [0.63928, 0.82913] | 120 | 0.2841 [0.1605, 0.3354] | 2.17e-6 | 39 [51] | 9.22e-7 | 2.21e-7 | 1.9 |
| 1e-07 | 571 | 120 | 120 | 0.78033 [0.73451, 0.87961] | 120 | 0.1991 [0.1100, 0.2402] | 5.43e-6 | 39 [51] | 1.69e-9 | 2.92e-10 | 1.3 |

At Part H's truncation (1e-13) basis pursuit on the grid is exact (t < 1) on 120 of 120 plane targets, on 37 of the 37 where the path fits the features with other kernels than E; t median by the path's outcome (E, other): 0.71195, 0.61285. t grows as the span is truncated more (a smaller span, a larger least t) on 120 of 120. t follows the grid's geometry at the features: where every feature replaced a grid point (one within half a spacing) t has the median [range] 0.60173 [0.56686, 0.62497] (62 targets), where some feature replaced none 0.72392 [0.70496, 0.78466] (58). The points within 10⁻⁷ of t (the program's active set), median per truncation: 1e-13: 1075, 1e-10: 821.5, 1e-07: 558; the equalities' largest residual in w: 2.11e-15.

At three digits (t at Part H's truncation is meaningful to about 10⁻³, see \|t_dir − t\|), the quantities the paper quotes: t at 1e-13 min, median, max 0.567, 0.620, 0.785; t median per truncation 1e-13: 0.620, 1e-10: 0.696, 1e-07: 0.780; t_J − t at 1e-13 min, max 0.205, 0.409; t median where every feature replaced a grid point, where some feature replaced none: 0.602, 0.724; t median by the path's outcome (E, other): 0.712, 0.613.

## M. The BLASSO in the plane (after round 149)

Part I's BLASSO on the plane targets: the same problem on Part H's basis-pursuit grid in the plane (spacing 0.025, its points within half a spacing of a feature replaced by the feature), λ = ℓλ_max for ℓ = 0.1, 0.01, 0.001, 0.0001, feature-sign search warm started from the larger λ, the Gram matrix applied through the samples' tensor factors (never formed); the largest KKT violation over all 480 solves 7.18e-12 relative to λ. Part I's recovery rule with the plane's tolerance: every nonzero within 2 spacings (0.05) of a feature, each feature's cluster nonempty with its sum of A_k's sign. Failures by Part K's types with the radius 0.5σ = 0.075 (σ = 0.15): a feature lost (no nonzero within 0.075 has it as its nearest feature), else a spurious nonzero (beyond 0.075 of every feature), else a displaced one (beyond the tolerance), else a cluster of the wrong sign:

| ℓ | targets | recovered | support exactly the features | a cluster empty within the tolerance | failures: lost, spurious, displaced, wrong sign | amplitude error median [max] | error/ℓ median | recovered: path E, path other | support beyond the features median [max] |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 120 | 118 | 64 | 0 | 0, 0, 2, 0 | 0.157 [0.301] | 1.57 | 83 of 83, 35 of 37 | 0 [6] |
| 0.01 | 120 | 120 | 64 | 0 | 0, 0, 0, 0 | 0.0159 [0.0347] | 1.59 | 83 of 83, 37 of 37 | 0 [6] |
| 0.001 | 120 | 120 | 64 | 0 | 0, 0, 0, 0 | 0.00159 [0.00347] | 1.59 | 83 of 83, 37 of 37 | 0 [6] |
| 0.0001 | 120 | 120 | 64 | 0 | 0, 0, 0, 0 | 1.59e-4 [3.47e-4] | 1.59 | 83 of 83, 37 of 37 | 0 [6] |

At ℓ = 0.0001: recovered on 120 of the 120 targets with η_V nondegenerate (Part H) and on 0 of the 0 without. Over the failures the displaced nonzeros' largest distance from their feature: 0.0526; the spurious ones' smallest distance from every feature: —; unclassified failures: 0. Feature-sign iterations per solve, median [max]: 2 [36]; wall time per target (four λ), median 0.048 s.

Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.
