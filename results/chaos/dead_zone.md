# Proposition D‴.57: the transit's dead zone (`exp_dead_zone.py`)

D‴.56's setting — 21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least 0.45 apart with m > 0.8q; y = f* + ξ, ξ i.i.d. N(0, s²) —, the gate open while the residual's rms exceeds (1 + δ)s, δ = 0.1.

## A. The sliding limit

Per level s: 10 draws on each of the 40 targets; the draws whose noisy path enters E first (no exit before the (K + 1)-th entry), split by whether it keeps the noiseless entry order (D‴.56's classification): the dead zone's stop (the first entry whose relaxed rms is ≤ (1 + δ)s) at the K-th entry, later (none of the first K) or earlier; at the K-th, the energy dead zone's frozen fit against the oracle's (the K-th block relaxed): the errors' ratio, the mean square over its prediction given the stop (± its standard error) and over s²(4K/N + 2δ + δ²) (the unconditional value), the oracle's over s²·2K/N, and the correlation of the oracle's error with the frozen part u = r_c − r_K:

| s | entry order | draws | stop at K, later, earlier | err_dz/err_or: median [5–95%] | mean err_dz²: / prediction given the stop; / s²(4K/N + 2δ + δ²) | mean err_or² / (s²2K/N) | corr(err_or, u): mean ± s.e. |
|---|---|---|---|---|---|---|---|
| 0.01 | noiseless | 395 | 381, 14, 0 | 2.147 [1.134, 3.819] | 1.016 ± 0.026; 1.064 | 0.985 | -0.022 ± 0.021 |
| 0.01 | swapped | 5 | 5, 0, 0 | 2.183 [1.613, 2.620] | 1.534 ± 0.272; 1.607 | 1.419 | -0.315 ± 0.082 |
| 0.03 | noiseless | 384 | 367, 17, 0 | 2.120 [1.124, 4.067] | 0.978 ± 0.026; 1.024 | 0.961 | +0.008 ± 0.020 |
| 0.03 | swapped | 15 | 14, 1, 0 | 2.121 [1.266, 2.759] | 1.047 ± 0.132; 1.097 | 1.034 | -0.034 ± 0.116 |
| 0.1 | noiseless | 306 | 294, 12, 0 | 2.078 [1.055, 3.984] | 0.957 ± 0.028; 1.003 | 0.973 | +0.007 ± 0.023 |
| 0.1 | swapped | 63 | 62, 1, 0 | 2.080 [1.199, 3.887] | 1.239 ± 0.075; 1.297 | 1.248 | -0.172 ± 0.050 |

The frozen part u at the K-th stop: the identity |u|² = (1 + δ)²s²N − |r_K|² (first order) with its largest deviation, the rms of |u|/(s√N) against its prediction given the stop, |cos| of u against the target's first draw with the same last entrant (that draw left out), and the last entrant's amplitude deficit δw = sign(w_last)(w_last of the oracle − w_last frozen), P⊥φ_last its kernel off the rest of the block's tangent at the oracle (the other kernels' coefficients and centres, its own centre) — at first order u = δw·P⊥φ_last:

| s | entry order | at K | \|u\|² / ((1 + δ)²s²N − \|r_K\|²): median [5–95%]; max\|· − 1\| | rms \|u\|/(s√N) (prediction given the stop) | \|cos\| of u, same last entrant: median, 5%, min (draws) | δw/s: median [5–95%] | δw/\|w_last\|: median [5–95%] | δw\|P⊥φ_last\|/\|u\|: median [range]; min \|cos(u, P⊥φ_last)\| |
|---|---|---|---|---|---|---|---|---|
| 0.01 | noiseless | 381 | 1.0000 [1.0000, 1.0000]; 7.1e-05 | 0.567 (0.572) | 1.0000, 0.9998, 0.9996 (341) | 1.76 [0.86, 2.33] | 0.032 [0.015, 0.049] | 0.9999 [0.9987, 1.0000]; 0.9987 |
| 0.01 | swapped | 5 | 1.0000 [1.0000, 1.0000]; 5.5e-07 | 0.618 (0.572) | 0.9999, 0.9999, 0.9999 (3) | 1.92 [0.93, 2.35] | 0.037 [0.019, 0.049] | 1.0000 [1.0000, 1.0000]; 1.0000 |
| 0.03 | noiseless | 367 | 1.0000 [1.0000, 1.0000]; 4.3e-04 | 0.568 (0.572) | 0.9998, 0.9982, 0.9954 (329) | 1.75 [0.79, 2.36] | 0.096 [0.040, 0.146] | 1.0000 [0.9993, 1.0000]; 0.9993 |
| 0.03 | swapped | 14 | 1.0000 [1.0000, 1.0000]; 5.1e-06 | 0.579 (0.572) | 0.9994, 0.9991, 0.9991 (10) | 1.73 [0.94, 2.27] | 0.085 [0.046, 0.125] | 1.0000 [1.0000, 1.0000]; 1.0000 |
| 0.1 | noiseless | 294 | 1.0000 [0.9999, 1.0002]; 7.5e-04 | 0.566 (0.572) | 0.9973, 0.9774, 0.9388 (257) | 1.76 [0.82, 2.37] | 0.324 [0.145, 0.490] | 1.0000 [0.9992, 1.0000]; 0.9997 |
| 0.1 | swapped | 62 | 1.0000 [0.9999, 1.0004]; 5.4e-04 | 0.580 (0.572) | 0.9970, 0.9886, 0.9872 (41) | 1.83 [0.64, 2.37] | 0.295 [0.153, 0.523] | 1.0000 [0.9996, 1.0000]; 0.9998 |

At first order |u|²/N = (1 + δ)²s² − |r_K|²/N, of mean s²(2δ + δ² + 2K/N); without the cross term the errors' rms ratio would be ((4K/N + 2δ + δ²)/(2K/N))^{1/2} = 2.033. The stop at the K-th entry is the event |r_K|² ≤ (1 + δ)²s²N, |r_K|²/s² ~ χ²_{N−2K} (probability 0.954); given it E|r_K|² = s²(N − 2K)P(χ²_{N−2K+2} ≤ c)/P(χ²_{N−2K} ≤ c), c = (1 + δ)²N, = 53.83s² against 55s², so E|u|²/N = 0.3276s² (rms |u|/(s√N) = 0.572) and E err_dz² = s²(2K + c − E[χ²_{N−2K} | ≤ c])/N = 0.4260s² (×1.047 the unconditional value; the rms ratio 2.081). With u = δw·P⊥φ_last, E(δw|P⊥φ_last|)² = s²(c − E[χ²_{N−2K} | ≤ c]) = 19.98s²: both |P⊥φ_last|² and c − E[χ²_{N−2K} | ≤ c] grow like N, so the deficit δw stays O(s) as N grows while the oracle's error falls like s(2K/N)^{1/2}: the rms ratio grows like (2 + N(2δ + δ²)/(2K))^{1/2} (A3).

### A3. The number of samples

The first 20 targets on N equispaced samples of [−1, 1] (the kernels and the features unchanged), s = 0.03, 20 draws per target (the targets whose noiseless path recovers the features exactly); the draws in the noiseless entry order (the swapped ones counted), at the stop at the K-th entry: the mean squared errors against their predictions given the stop (stop_law(N)), the errors' rms ratio, and the last entrant's amplitude deficit:

| N | targets: exact noiseless path | draws in the noiseless order: stop at K, later, earlier; swapped | mean err_dz²/s²: measured (prediction) | mean err_or² / (s²2K/N) | rms err_dz / rms err_or: measured (prediction) | δw/s: median [5–95%] | rms δw\|P⊥φ_last\|/s: measured (prediction (c − E[χ²_{N−2K} \| ≤ c])^{1/2}) |
|---|---|---|---|---|---|---|---|
| 61 | 20/20 | 366, 19, 0; 15 | 0.4328 (0.4260) | 1.013 | 2.084 (2.081) | 1.80 [0.79, 2.34] | 4.475 (4.470) |
| 181 | 20/20 | 386, 8, 0; 6 | 0.2735 (0.2802) | 1.014 | 2.852 (2.907) | 1.49 [0.91, 1.99] | 6.619 (6.687) |
| 541 | 20/20 | 398, 0, 0; 2 | 0.2327 (0.2323) | 0.928 | 4.756 (4.576) | 1.45 [1.04, 1.75] | 10.906 (10.938) |
| 1621 | 20/20 | 399, 0, 0; 1 | 0.2205 (0.2174) | 1.002 | 7.710 (7.664) | 1.43 [1.23, 1.61] | 18.756 (18.612) |

## B. The loop

The loop at ε = 1e-08 from α = 0 on 16 targets × 3 draws at s = 0.01, 0.03, 0.1 (the draws whose noisy path is exactly E; the noise drawn afresh at each level), to τ_K + 3W: on the draws whose sliding limit stops at the K-th entry, the kernels formed (|w| > 0.001) and the fits' rms errors over the oracle's (the limit's K-th block relaxed), the energy dead zone's against the ungated loop's, and the gated loops' fits against the sliding limit's (the energy dead zone's frozen fit f̂_dz, the oracle's f̂_or: the rms of the difference); on the late draws (the stop after the K-th entry) the kernels formed:

| s | draws: stop at K, late | formed at K: ungated, energy, formation (median [range]) | err / err_or (medians): ungated, energy, formation | energy against the ungated (median) | the loop against the limit: max rms(f̂ − f̂_dz) (energy), max rms(f̂ − f̂_or) (formation) | the gate's τ / τ_K (median) | formed, late: ungated, energy, formation |
|---|---|---|---|---|---|---|---|
| 0.01 | 46, 2 | 5 [4, 10]; 3 [3, 3]; 3 [3, 3] | 1.423, 2.236, 1.000 | 1.524 | 3.6e-10, 2.3e-10 | 1.0011 | 7, 7, 7; 4, 4, 4 |
| 0.03 | 45, 3 | 6 [4, 9]; 3 [3, 3]; 3 [3, 3] | 1.499, 2.234, 1.000 | 1.414 | 3.3e-10, 3.3e-10 | 1.0009 | 5, 5, 5; 8, 5, 5; 6, 4, 4 |
| 0.1 | 43, 0 | 5 [4, 8]; 3 [3, 3]; 3 [3, 3] | 1.492, 2.315, 1.000 | 1.578 | 9.9e-10, 4.3e-10 | 1.0007 | — |

## C. The closed loop under measurement noise

D‴.49's Part E (k_e = 10, ω = 0.2, three bumps of width 0.12; 21 kernels from α = 0 at ε = 0.0001 to τ = 600; the noise 16 tones in 10–30 rad/s of rms s, seed 7); the dead zone on the error's energy filtered over T_f = 5 time units, the drive gated smoothly between the rms 1.2s and 1.5s (every kernel, or the kernels with |w| below about the threshold); formed: |w| > 1% of max|A*| (a kernel's first crossing), ejected: |c| > 1.2 at an output time after its formation, outside: |c| > 1.2 at the end; 1201 output times on [0, T]: the fit's rms error on the reference's range [−1, 1] (every kernel) and the formed kernels inside over the last 50% of the run, the true error x − x_m's rms and the measured error's largest value over the last period, the gate's state over time (shut: the filtered energy at or below the lower threshold; fully open: at or above the upper):

| regressor | s | adaptation | formed, ejected | at the end: formed, outside | fit rms on [−1, 1], last half: median [5–95%] | formed inside, last half (median) | true error rms, last period | max\|ẽ\|, last period | the gate: shut, fully open (shares of time) | the gate: first closing, last reopening (τ) |
|---|---|---|---|---|---|---|---|---|---|---|
| x̃ | 0 | ungated | 5, 0 | 5, 0 | 0.0003 [0.0002, 0.0006] | 5 | 0.0000 | 0.0000 | — | — |
| x̃ | 0.01 | ungated | 14, 14 | 14, 14 | 0.3825 [0.3333, 0.4272] | 1 | 0.0372 | 0.1011 | — | — |
| x̃ | 0.01 | the mirror leak λΛ = 0.003 | 5, 0 | 3, 0 | 0.2837 [0.2636, 0.2994] | 3 | 0.0186 | 0.0624 | — | — |
| x̃ | 0.01 | the energy dead zone | 13, 11 | 13, 11 | 0.3815 [0.2451, 0.3834] | 2 | 0.0339 | 0.1010 | 0.01, 0.94 | 187.0, 562.0 |
| x̃ | 0.01 | the formation dead zone | 14, 14 | 14, 14 | 0.3825 [0.3321, 0.4272] | 1 | 0.0372 | 0.1011 | 0.00, 0.98 | —, — |
| x̃ | 0.02 | ungated | 16, 16 | 16, 16 | 0.4274 [0.4126, 0.4276] | 0 | 0.0380 | 0.1170 | — | — |
| x̃ | 0.02 | the mirror leak λΛ = 0.003 | 11, 8 | 3, 0 | 0.4033 [0.3674, 0.4257] | 3 | 0.0358 | 0.1171 | — | — |
| x̃ | 0.02 | the energy dead zone | 16, 16 | 16, 16 | 0.4273 [0.4113, 0.4276] | 0 | 0.0380 | 0.1171 | 0.00, 0.93 | —, — |
| x̃ | 0.02 | the formation dead zone | 16, 16 | 16, 16 | 0.4273 [0.4119, 0.4276] | 0 | 0.0380 | 0.1174 | 0.00, 0.93 | —, — |
| x_m | 0 | ungated | 3, 0 | 3, 0 | 0.0000 [0.0000, 0.0002] | 3 | 0.0000 | 0.0000 | — | — |
| x_m | 0.01 | ungated | 3, 0 | 3, 0 | 0.0254 [0.0153, 0.0707] | 3 | 0.0071 | 0.0345 | — | — |
| x_m | 0.01 | the mirror leak λΛ = 0.003 | 3, 0 | 3, 0 | 0.2800 [0.2613, 0.2981] | 2 | 0.0187 | 0.0624 | — | — |
| x_m | 0.01 | the energy dead zone | 3, 0 | 3, 0 | 0.0613 [0.0565, 0.0700] | 3 | 0.0068 | 0.0339 | 0.55, 0.37 | 174.5, 599.0 |
| x_m | 0.01 | the formation dead zone | 3, 0 | 3, 0 | 0.0254 [0.0153, 0.0708] | 3 | 0.0071 | 0.0345 | 0.55, 0.37 | 174.5, 600.0 |
| x_m | 0.02 | ungated | 12, 11 | 12, 11 | 0.3823 [0.2476, 0.3883] | 1 | 0.0348 | 0.1170 | — | — |
| x_m | 0.02 | the mirror leak λΛ = 0.003 | 3, 0 | 3, 0 | 0.2783 [0.2584, 0.2974] | 2 | 0.0203 | 0.0835 | — | — |
| x_m | 0.02 | the energy dead zone | 7, 3 | 7, 3 | 0.1288 [0.1195, 0.3506] | 4 | 0.0135 | 0.0686 | 0.54, 0.30 | 149.5, 599.0 |
| x_m | 0.02 | the formation dead zone | 12, 11 | 12, 11 | 0.3829 [0.1569, 0.3883] | 1 | 0.0348 | 0.1174 | 0.11, 0.72 | 149.5, 568.5 |

The formed kernels at the end (snapshots of a moving state; the window's statistics above):

| regressor | s | adaptation | coefficients | centres |
|---|---|---|---|---|
| x̃ | 0 | ungated | 0.900, -0.535, -0.065, 0.031, 0.720 | -0.520, 0.070, 0.070, 0.631, 0.609 |
| x̃ | 0.01 | ungated | 0.854, 0.849, 0.785, 0.901, 0.659, 0.816, -0.568, -0.639, 0.727, 0.730, 0.486, 0.659, 0.707, 0.687 | 1.394, 1.394, 1.394, 1.394, 1.394, 1.394, -1.334, -1.333, 1.394, 1.394, 1.393, 1.394, 1.394, 1.394 |
| x̃ | 0.01 | the mirror leak λΛ = 0.003 | 0.434, -0.159, 0.330 | -0.472, 0.170, 0.603 |
| x̃ | 0.01 | the energy dead zone | 0.874, 0.738, 0.906, 0.957, 0.934, -0.531, -0.079, 0.737, 0.611, 0.717, 0.707, 0.719, 0.617 | 1.393, 1.392, 1.393, 1.393, 1.393, 0.050, 0.068, 1.392, 1.392, 1.392, 1.392, 1.392, 1.392 |
| x̃ | 0.01 | the formation dead zone | 0.855, 0.847, 0.784, 0.901, 0.659, 0.818, -0.566, -0.619, 0.669, 0.729, 0.486, 0.659, 0.706, 0.687 | 1.394, 1.394, 1.394, 1.394, 1.393, 1.394, -1.333, -1.333, 1.393, 1.394, 1.393, 1.393, 1.393, 1.393 |
| x̃ | 0.02 | ungated | 0.369, 0.444, 0.444, 0.540, 0.420, 0.398, -0.398, -0.343, -0.340, -0.312, 0.239, 0.284, 0.291, 0.302, 0.270, 0.207 | 1.328, 1.330, 1.330, 1.332, 1.329, 1.329, -1.293, -1.291, -1.291, -1.290, 1.325, 1.326, 1.326, 1.327, 1.326, 1.325 |
| x̃ | 0.02 | the mirror leak λΛ = 0.003 | -0.145, 0.205, 0.028 | 0.142, 1.179, 0.686 |
| x̃ | 0.02 | the energy dead zone | 0.357, 0.346, 0.358, 0.476, 0.464, 0.370, -0.400, -0.337, -0.335, -0.297, 0.221, 0.247, 0.222, 0.221, 0.244, 0.208 | 1.324, 1.324, 1.324, 1.327, 1.327, 1.324, -1.292, -1.290, -1.290, -1.288, 1.321, 1.321, 1.321, 1.321, 1.321, 1.320 |
| x̃ | 0.02 | the formation dead zone | 0.406, 0.353, 0.410, 0.536, 0.408, 0.470, -0.399, -0.336, -0.330, -0.299, 0.222, 0.251, 0.246, 0.247, 0.243, 0.234 | 1.328, 1.326, 1.328, 1.330, 1.328, 1.329, -1.292, -1.290, -1.289, -1.288, 1.323, 1.324, 1.324, 1.324, 1.324, 1.324 |
| x_m | 0 | ungated | 0.898, -0.600, 0.750 | -0.520, 0.070, 0.610 |
| x_m | 0.01 | ungated | 0.930, -0.609, 0.765 | -0.527, 0.058, 0.648 |
| x_m | 0.01 | the mirror leak λΛ = 0.003 | 0.447, -0.164, 0.340 | -0.470, 0.169, 0.564 |
| x_m | 0.01 | the energy dead zone | 0.923, -0.435, 0.762 | -0.524, 0.043, 0.621 |
| x_m | 0.01 | the formation dead zone | 0.930, -0.610, 0.765 | -0.526, 0.058, 0.648 |
| x_m | 0.02 | ungated | 1.002, 1.023, 1.014, 1.036, 1.009, 0.984, -0.645, 1.036, 0.743, 0.813, 0.843, 1.039 | -1.495, 1.557, -1.515, 1.561, 1.552, 1.361, 0.009, -1.476, 1.547, 1.553, 1.525, 1.499 |
| x_m | 0.02 | the mirror leak λΛ = 0.003 | 0.452, -0.165, 0.343 | -0.473, 0.169, 0.568 |
| x_m | 0.02 | the energy dead zone | 0.884, 0.976, 0.995, -0.253, 0.670, 0.843, 0.010 | -0.522, -1.492, -1.473, 0.009, 0.618, 1.467, 0.614 |
| x_m | 0.02 | the formation dead zone | 1.006, 0.998, 1.014, 1.028, 1.017, 1.012, -0.645, 1.035, 1.013, 0.807, 0.583, 0.784 | 1.539, -1.508, -1.522, 1.567, -1.480, 1.363, 0.009, 1.484, 1.562, 1.559, 1.527, 1.524 |

Tests of the mechanism at x_m, s = 0.02 unless stated (referee round 50): the adaptation driven by the true error x − x_m (the control and the gate by the measured one), the centres bounded to [−1.2, 1.2], other gate thresholds and filters, another noise seed; C2: the noise's structure — one tone, two tones far apart, two close ones (the same rms) — and the centres' Wasserstein–Fisher–Rao metric (mobility G·a, bounded, against the paradigm's Gm|w|):

| variant | adaptation | formed, ejected | at the end: formed, outside | fit rms on [−1, 1], last half: median [5–95%] | formed inside, last half (median) | true error rms, last period | max\|ẽ\|, last period | the gate: shut, fully open (shares of time) | the gate: first closing, last reopening (τ) |
|---|---|---|---|---|---|---|---|---|---|
| true error | ungated | 8, 6 | 8, 6 | 0.2943 [0.0327, 0.2997] | 2 | 0.0276 | 0.1169 | — | — |
| box 1.2 | ungated | 12, 11 | 7, 0 | 0.3839 [0.2451, 0.3906] | 8 | 0.0349 | 0.1170 | — | — |
| box 1.2 | the energy dead zone | 7, 3 | 7, 0 | 0.1239 [0.0985, 0.3194] | 7 | 0.0135 | 0.0672 | 0.53, 0.30 | 149.5, 586.5 |
| thresholds 2.0s, 2.5s | the energy dead zone | 1, 0 | 1, 0 | 0.3487 [0.3452, 0.4271] | 1 | 0.0305 | 0.1128 | 0.57, 0.02 | 1.0, 600.0 |
| thresholds 2.0s, 2.5s | the formation dead zone | 5, 3 | 5, 3 | 0.4265 [0.3149, 0.4271] | 0 | 0.0285 | 0.1022 | 0.48, 0.04 | 1.0, 600.0 |
| T_f = 40 | the energy dead zone | 3, 1 | 3, 1 | 0.1969 [0.1962, 0.1990] | 2 | 0.0181 | 0.0875 | 0.41, 0.28 | 194.0, 600.0 |
| T_f = 40 | the formation dead zone | 12, 11 | 12, 11 | 0.3831 [0.2487, 0.3904] | 1 | 0.0348 | 0.1172 | 0.01, 0.84 | 180.0, 382.0 |
| noise seed 8 | ungated | 12, 11 | 12, 11 | 0.3832 [0.2465, 0.3873] | 1 | 0.0351 | 0.1278 | — | — |
| noise seed 8 | the energy dead zone | 7, 3 | 7, 3 | 0.1327 [0.1198, 0.1979] | 4 | 0.0140 | 0.0671 | 0.54, 0.30 | 146.5, 599.0 |
| one tone, 21.37 rad/s, s = 0.02 | ungated | 3, 0 | 3, 0 | 0.0479 [0.0401, 0.1230] | 3 | 0.0112 | 0.0454 | — | — |
| one tone, 21.37 rad/s, s = 0.04 | ungated | 10, 6 | 10, 6 | 0.5405 [0.2258, 0.7014] | 4 | 0.0260 | 0.1169 | — | — |
| two tones, 13.1 and 27.9 rad/s, s = 0.02 | ungated | 4, 0 | 4, 0 | 0.1242 [0.0926, 0.2150] | 3 | 0.0165 | 0.0608 | — | — |
| two tones, 21.37 and 21.44 rad/s, s = 0.02 | ungated | 3, 0 | 3, 0 | 0.0506 [0.0319, 0.1270] | 3 | 0.0131 | 0.0726 | — | — |
| WFR metric, s = 0 | ungated | 15, 0 | 15, 0 | 0.0000 [0.0000, 0.0000] | 15 | 0.0000 | 0.0000 | — | — |
| WFR metric, s = 0.02 | ungated | 15, 0 | 15, 0 | 0.0041 [0.0032, 0.0109] | 15 | 0.0095 | 0.0524 | — | — |
| WFR metric, s = 0.03 | ungated | 15, 0 | 15, 0 | 0.0087 [0.0070, 0.0171] | 15 | 0.0143 | 0.0790 | — | — |

The kernels' timelines at x_m, s = 0.02: kernel: formed at τ (the feature nearest its centre then) → ejected beyond |c| = 1.2 at τ:

| run | timeline |
|---|---|
| base, ungated | 4: 136.4 (-0.52) → 164.5; 16: 163.6 (+0.61) → 197.0; 5: 167.7 (-0.52) → 197.0; 3: 206.5 (-0.52) → 234.5; 15: 207.4 (+0.61) → 234.5; 17: 213.7 (+0.61) → 294.5; 11: 253.8 (+0.07); 6: 267.4 (-0.52) → 294.5; 2: 362.9 (-0.52) → 389.0; 18: 375.6 (+0.61) → 436.0; 14: 441.5 (+0.61) → 496.5; 7: 527.8 (-0.52) → 553.5 |
| base, the mirror leak λΛ = 0.003 | 4: 175.3 (-0.52); 16: 225.0 (+0.61); 11: 467.0 (+0.07) |
| base, the energy dead zone | 4: 136.5 (-0.52) → 186.0; 16: 165.5 (+0.61) → 321.5; 5: 189.2 (-0.52) → 309.0; 3: 318.4 (-0.52); 15: 330.3 (+0.61); 11: 334.3 (+0.07); 17: 354.5 (+0.61) |
| base, the formation dead zone | 4: 136.5 (-0.52) → 164.5; 16: 165.4 (+0.61) → 197.0; 5: 168.6 (-0.52) → 197.0; 3: 207.8 (-0.52) → 233.5; 15: 207.9 (+0.61) → 282.5; 6: 267.6 (-0.52) → 309.0; 17: 268.7 (+0.61) → 282.5; 11: 286.4 (+0.07); 2: 377.7 (-0.52) → 405.0; 18: 378.0 (+0.61) → 405.0; 14: 414.0 (+0.61) → 472.0; 7: 540.9 (-0.52) → 568.5 |
| true error, ungated | 4: 136.4 (-0.52) → 173.5; 16: 163.6 (+0.61); 5: 177.0 (-0.52) → 217.5; 3: 226.7 (-0.52) → 264.0; 11: 253.5 (+0.07); 6: 297.0 (-0.52) → 333.0; 2: 400.4 (-0.52) → 440.0; 7: 550.9 (-0.52) → 588.5 |
| box 1.2, ungated | 4: 136.4 (-0.52) → 168.5; 16: 163.6 (+0.61) → 197.5; 5: 167.7 (-0.52) → 197.5; 3: 206.5 (-0.52) → 236.5; 15: 207.5 (+0.61) → 236.5; 17: 213.9 (+0.61) → 296.0; 11: 253.8 (+0.07); 6: 267.3 (-0.52) → 287.5; 2: 379.4 (-0.52) → 405.5; 18: 381.1 (+0.61) → 405.5; 14: 412.2 (+0.61) → 472.0; 7: 513.1 (-0.52) → 558.0 |
| box 1.2, the energy dead zone | 4: 136.5 (-0.52) → 188.0; 16: 165.5 (+0.61) → 303.0; 5: 189.2 (-0.52) → 310.0; 15: 310.2 (+0.61); 3: 318.2 (-0.52); 11: 325.3 (+0.07); 17: 330.0 (+0.61) |
| thresholds 2.0s, 2.5s, the energy dead zone | 4: 336.7 (-0.52) |
| thresholds 2.0s, 2.5s, the formation dead zone | 4: 335.5 (-0.52) → 364.5; 3: 379.7 (-0.52) → 405.0; 5: 449.7 (-0.52) → 474.0; 17: 561.6 (+0.61); 2: 595.2 (-0.52) |
| T_f = 40, the energy dead zone | 4: 136.4 (-0.52) → 164.5; 16: 163.6 (+0.61); 5: 167.7 (-0.52) |
| T_f = 40, the formation dead zone | 4: 136.4 (-0.52) → 164.5; 16: 163.6 (+0.61) → 197.0; 5: 167.7 (-0.52) → 197.0; 3: 207.9 (-0.52) → 237.0; 15: 208.3 (+0.61) → 237.0; 17: 214.6 (+0.61) → 287.0; 6: 271.0 (-0.52) → 287.0; 11: 297.8 (+0.07); 2: 353.9 (-0.52) → 382.0; 18: 369.6 (+0.61) → 382.0; 14: 407.9 (+0.61) → 472.0; 7: 508.1 (-0.52) → 538.0 |
| noise seed 8, ungated | 4: 136.4 (-0.52) → 161.5; 16: 163.7 (+0.61) → 208.0; 5: 165.0 (-0.52) → 191.5; 15: 208.6 (+0.61) → 236.5; 3: 210.3 (-0.52) → 236.5; 17: 214.9 (+0.61) → 281.0; 11: 253.8 (+0.07); 6: 269.5 (-0.52) → 303.5; 18: 353.9 (+0.61) → 393.0; 2: 371.9 (-0.52) → 393.0; 14: 414.4 (+0.61) → 468.0; 7: 517.3 (-0.52) → 541.5 |
| noise seed 8, the energy dead zone | 4: 136.5 (-0.52) → 198.0; 16: 165.8 (+0.61) → 239.0; 5: 201.1 (-0.52) → 260.0; 15: 248.2 (+0.61); 3: 269.6 (-0.52); 17: 276.7 (+0.61); 11: 330.5 (+0.07) |
| one tone, 21.37 rad/s, s = 0.02, ungated | 4: 136.4 (-0.52); 16: 163.6 (+0.61); 11: 253.5 (+0.07) |
| one tone, 21.37 rad/s, s = 0.04, ungated | 4: 136.8 (-0.52) → 228.5; 16: 164.1 (+0.61) → 228.5; 5: 202.7 (-0.52) → 339.5; 15: 237.0 (+0.61); 17: 242.7 (+0.61); 11: 254.8 (+0.07) → 340.0; 3: 341.9 (-0.52) → 511.5; 10: 359.0 (+0.07) → 511.5; 12: 511.3 (+0.07); 6: 521.5 (-0.52) |
| two tones, 13.1 and 27.9 rad/s, s = 0.02, ungated | 4: 136.4 (-0.52); 16: 163.7 (+0.61); 11: 253.7 (+0.07); 5: 481.8 (-0.52) |
| two tones, 21.37 and 21.44 rad/s, s = 0.02, ungated | 4: 136.4 (-0.52); 16: 163.6 (+0.61); 11: 253.4 (+0.07) |
| WFR metric, s = 0, ungated | 4: 132.6 (-0.52); 5: 132.9 (-0.52); 3: 133.4 (-0.52); 6: 134.8 (-0.52); 2: 136.9 (-0.52); 16: 160.2 (+0.61); 15: 160.8 (+0.61); 17: 161.0 (+0.61); 14: 163.9 (+0.61); 18: 164.1 (+0.61); 11: 250.1 (+0.07); 10: 251.1 (+0.07); 12: 251.5 (+0.07); 9: 255.8 (+0.07); 13: 266.5 (+0.07) |
| WFR metric, s = 0.02, ungated | 4: 132.8 (-0.52); 5: 133.1 (-0.52); 3: 133.6 (-0.52); 6: 135.0 (-0.52); 2: 137.1 (-0.52); 16: 160.4 (+0.61); 15: 161.0 (+0.61); 17: 161.2 (+0.61); 14: 164.1 (+0.61); 18: 164.3 (+0.61); 11: 250.5 (+0.07); 10: 251.5 (+0.07); 12: 251.8 (+0.07); 9: 256.2 (+0.07); 13: 266.8 (+0.07) |
| WFR metric, s = 0.03, ungated | 4: 133.0 (-0.52); 5: 133.3 (-0.52); 3: 133.8 (-0.52); 6: 135.2 (-0.52); 2: 137.3 (-0.52); 16: 160.7 (+0.61); 15: 161.3 (+0.61); 17: 161.5 (+0.61); 14: 164.4 (+0.61); 18: 164.5 (+0.61); 11: 250.9 (+0.07); 10: 251.9 (+0.07); 12: 252.3 (+0.07); 9: 256.6 (+0.07); 13: 267.2 (+0.07) |

Near the fit the loop high-passes the noise in the measured error (ẽ = e + n, e's response to n through −k_eẽ: |iω/(iω + k_e)|), a floor of rms 0.881s with the desired regressor (at x̃ the regressor's noise adds f̂′(x)n): the energy dead zone starts to close when the slow error's rms falls to 1.21s and is shut below 0.81s.
The ungated loops at the end against the noise-free loop with the same regressor (the formed kernels in the order of their centres): x̃, s = 0.01: 14 kernels against 5; x_m, s = 0.01: max|Δw| 3.16e-02, max|Δc| 3.79e-02; x̃, s = 0.02: 16 kernels against 5; x_m, s = 0.02: 12 kernels against 3.

## D. From the fitted state: the formed kernels under measurement noise

C's closed loop (the regressor at x_m, ungated, the 16 tones of rms s) started with each feature fitted by one kernel — the kernel started nearest it, at the feature's centre with the coefficient λA*_k —, for 300 periods of the reference: per kernel the first period at which its period-averaged centre is 0.3 from its feature (— if never) and the coefficient's growth w/(λA*) − 1 at the end; the paradigm's metric (the centres' mobility Gβa = Gm|w|) at λ = 0.5, 0.75, 1, 1.5, the Wasserstein–Fisher–Rao metric (G·a, at most G·a_max) at λ = 1:

| metric | λ | s | feature 1: left at period, w/(λA*) − 1 | feature 2: left at period, w/(λA*) − 1 | feature 3: left at period, w/(λA*) − 1 | the others' max\|w\| |
|---|---|---|---|---|---|---|
| paradigm | 0.5 | 0.005 | —, +0.002 | —, +0.001 | —, +0.001 | 1.3e-11 |
| paradigm | 0.5 | 0.0075 | —, +0.004 | —, +0.002 | —, +0.003 | 2.9e-11 |
| paradigm | 0.5 | 0.01 | —, +0.007 | —, +0.004 | —, +0.005 | 5.1e-11 |
| paradigm | 0.5 | 0.0125 | —, +0.011 | —, +0.006 | —, +0.008 | 8.0e-11 |
| paradigm | 0.5 | 0.015 | —, +0.015 | —, +0.008 | —, +0.012 | 1.2e-10 |
| paradigm | 0.5 | 0.0175 | —, +0.021 | —, +0.011 | —, +0.016 | 1.6e-10 |
| paradigm | 0.5 | 0.02 | —, +0.028 | —, +0.015 | —, +0.021 | 2.1e-10 |
| paradigm | 0.5 | 0.025 | —, +0.044 | —, +0.023 | —, +0.033 | 3.3e-10 |
| paradigm | 0.5 | 0.03 | —, +0.065 | —, +0.033 | —, +0.048 | 4.9e-10 |
| paradigm | 0.5 | 0.04 | —, +0.122 | —, +0.058 | —, +0.086 | 9.5e-10 |
| paradigm | 0.75 | 0.005 | —, +0.003 | —, +0.002 | —, +0.002 | 2.6e-11 |
| paradigm | 0.75 | 0.0075 | —, +0.008 | —, +0.004 | —, +0.006 | 5.9e-11 |
| paradigm | 0.75 | 0.01 | —, +0.014 | —, +0.008 | —, +0.010 | 1.1e-10 |
| paradigm | 0.75 | 0.0125 | —, +0.023 | —, +0.012 | —, +0.016 | 1.7e-10 |
| paradigm | 0.75 | 0.015 | —, +0.034 | —, +0.017 | —, +0.023 | 2.6e-10 |
| paradigm | 0.75 | 0.0175 | —, +0.049 | —, +0.023 | —, +0.032 | 3.6e-10 |
| paradigm | 0.75 | 0.02 | —, +0.068 | —, +0.030 | —, +0.044 | 5.0e-10 |
| paradigm | 0.75 | 0.025 | —, +0.115 | —, +0.047 | —, +0.073 | 8.5e-10 |
| paradigm | 0.75 | 0.03 | —, +0.186 | —, +0.069 | —, +0.109 | 1.4e-09 |
| paradigm | 0.75 | 0.04 | —, +0.341 | —, +0.121 | —, +0.215 | 2.6e-09 |
| paradigm | 1 | 0.005 | —, +0.006 | —, +0.003 | —, +0.004 | 4.9e-11 |
| paradigm | 1 | 0.0075 | —, +0.016 | —, +0.007 | —, +0.010 | 1.2e-10 |
| paradigm | 1 | 0.01 | —, +0.031 | —, +0.013 | —, +0.018 | 2.3e-10 |
| paradigm | 1 | 0.0125 | —, +0.054 | —, +0.021 | —, +0.030 | 4.1e-10 |
| paradigm | 1 | 0.015 | —, +0.083 | —, +0.030 | —, +0.047 | 6.1e-10 |
| paradigm | 1 | 0.0175 | 232, +0.105 | —, +0.042 | —, +0.070 | 1.4e-08 |
| paradigm | 1 | 0.02 | 133, +0.093 | —, +0.056 | —, +0.096 | 3.3e-08 |
| paradigm | 1 | 0.025 | 118, +0.134 | —, +0.090 | 221, +0.248 | 2.2e-08 |
| paradigm | 1 | 0.03 | 147, +0.208 | —, +0.131 | 147, +0.148 | 3.1e-08 |
| paradigm | 1 | 0.04 | 209, +0.450 | —, +0.245 | 210, +0.325 | 2.1e-08 |
| paradigm | 1.5 | 0.005 | —, +0.021 | —, +0.009 | —, +0.013 | 1.5e-10 |
| paradigm | 1.5 | 0.0075 | 205, +0.044 | —, +0.022 | —, +0.029 | 2.9e-08 |
| paradigm | 1.5 | 0.01 | 193, +0.077 | —, +0.037 | —, +0.060 | 3.3e-08 |
| paradigm | 1.5 | 0.0125 | 104, +0.081 | —, +0.062 | 104, +0.050 | 6.0e-08 |
| paradigm | 1.5 | 0.015 | 112, +0.120 | —, +0.090 | 112, +0.079 | 5.8e-08 |
| paradigm | 1.5 | 0.0175 | 132, +0.177 | —, +0.125 | 131, +0.126 | 5.2e-08 |
| paradigm | 1.5 | 0.02 | 154, +0.249 | 291, +0.158 | 125, +0.155 | 4.6e-08 |
| paradigm | 1.5 | 0.025 | 168, +0.371 | 131, +0.128 | 74, +0.141 | 4.9e-08 |
| paradigm | 1.5 | 0.03 | 168, +0.469 | 198, +0.256 | 74, +0.182 | 5.0e-08 |
| paradigm | 1.5 | 0.04 | 41, +0.830 | 198, +0.388 | 190, +0.572 | 2.8e-08 |
| wfr | 1 | 0.02 | —, +0.006 | —, +0.006 | —, +0.007 | 1.0e-10 |
| wfr | 1 | 0.04 | —, +0.025 | —, +0.023 | —, +0.027 | 4.3e-10 |
| wfr | 1 | 0.06 | —, +0.055 | —, +0.050 | —, +0.058 | 9.9e-10 |
| wfr | 1 | 0.08 | —, +0.093 | —, +0.087 | —, +0.098 | 1.8e-09 |
| wfr | 1 | 0.12 | —, +0.186 | —, +0.184 | —, +0.196 | 4.1e-09 |
| wfr | 1 | 0.16 | —, +0.292 | —, +0.318 | —, +0.310 | 6.8e-09 |

The smallest level at which each kernel left its feature: paradigm, |A| = 0.45: kept to s = 0.04; paradigm, |A| = 0.3: kept to s = 0.04; paradigm, |A| = 0.375: kept to s = 0.04; paradigm, |A| = 0.675: kept to s = 0.04; paradigm, |A| = 0.45: kept to s = 0.04; paradigm, |A| = 0.562: kept to s = 0.04; paradigm, |A| = 0.9: 0.0175 (s·|A| = 0.0158); paradigm, |A| = 0.6: kept to s = 0.04; paradigm, |A| = 0.75: 0.025 (s·|A| = 0.0188); paradigm, |A| = 1.35: 0.0075 (s·|A| = 0.0101); paradigm, |A| = 0.9: 0.02 (s·|A| = 0.0180); paradigm, |A| = 1.12: 0.0125 (s·|A| = 0.0141); wfr, |A| = 0.9: kept to s = 0.16; wfr, |A| = 0.6: kept to s = 0.16; wfr, |A| = 0.75: kept to s = 0.16.

The kept kernels' coefficient growth over (|A|s)² under the paradigm's metric: median 373 [5–95% 312, 470] (92 kernels and levels).

Numbers are produced by `experiments/exp_dead_zone.py`; nothing is edited by hand.
