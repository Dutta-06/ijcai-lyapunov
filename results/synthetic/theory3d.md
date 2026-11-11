# Which readouts select: theorems for the readouts of depth p, their checks, and the mismatch margin (round theory3, direction D)

The statements (with the proofs, in the note outside the repository) concern the SCALAR unit: whatever the other units and its own centre do, a unit of depth p enters the network only through its dual eta = int g, g = <r phi(.;c)> the drive, and its coefficient is an explicit function of eta. Everything below tests this against the exact ODE (A), ball arithmetic (B), the repository's network flows with moving centres (C), the depth-1 flow (D), the sliding unit (E), the fixed-centre limit (L), a family of ties at fixed geometry (K), the paths of the benchmark population (P), the stored runs of 200 benchmark targets (F), the per-target scan (G) and the population model (H). Status of each statement is in the last section.

## A. The closed forms against the coordinate ODE (prescribed drive)

Drive g(t) = Amp (sin 0.7t + 0.3), Amp so that max|eta| = 0.9 wall (0.5 for the paper's unit); DOP853 with rtol 1e-12; w from (u, v) against W(eta) with eta from the same integration; row p = 0 is the paper's unit ((alpha, nu) ODE against Lemma 2.1's dual; the last column is the drift of the conserved Q = nu^2 - v_s cosh(2 alpha/a_max) over v_s).

| p | eps | ell | wall | max eta/wall | err eta/wall | max rel err w (p=0: Q drift) | max w |
|---|---|---|---|---|---|---|---|
| 2 | 0.001 | 0.7 | 2.47 | 0.9 | 1.64e-14 | 4.9e-12 | 0.501 |
| 3 | 0.001 | 0.7 | 15.1 | 0.9 | 1.06e-13 | 9.15e-11 | 0.0316 |
| 4 | 0.001 | 0.7 | 179 | 0.9 | 1.57e-13 | 9.88e-11 | 9.97e-05 |
| 5 | 0.001 | 0.7 | 3.01e+03 | 0.9 | 2.19e-13 | 8.74e-11 | 1.46e-06 |
| 6 | 0.001 | 0.7 | 5.95e+04 | 0.9 | 1.5e-13 | 3.19e-11 | 3.12e-08 |
| 0 | 0.0001 | 1 | 702 | 0.5 | 2.33e-12 | 7.36e-15 | 7.42e-07 |

The paper's unit: exact dormant supremum, the cos/sin^3 form and the mu^-3 form of Lemma 2.2 (T1 note), and the leading constant 8/pi^3 a/(m sqrt(2 v_s)) eps^(3/2) mu^-3:

| eps | mu | exact sup | cos/sin^3 form | mu^-3 form | leading | exact <= cos/sin^3 <= mu^-3 |
|---|---|---|---|---|---|---|
| 0.0001 | 0.5 | 7.42e-07 | 7.49e-07 | 3e-06 | 7.69e-07 | yes |
| 0.0001 | 0.2 | 1.16e-05 | 1.23e-05 | 4.76e-05 | 1.2e-05 | yes |
| 0.0001 | 0.05 | 0.000465 | 0.00114 | 0.00442 | 0.000769 | yes |
| 0.0001 | 0.025 | 0.00164 | -0.0207 | -0.0802 | 0.00615 | n/a (s_eps >= 1: the lemma needs cot(pi mu/2) < kappa) |
| 0.0001 | 0.01 | 0.00408 | -0.0135 | -0.0525 | 0.0962 | n/a (s_eps >= 1: the lemma needs cot(pi mu/2) < kappa) |
| 0.0001 | 0.001 | 0.00738 | -0.119 | -0.461 | 96.2 | n/a (s_eps >= 1: the lemma needs cot(pi mu/2) < kappa) |
| 1e-06 | 0.5 | 7.45e-10 | 7.45e-10 | 2.98e-09 | 7.69e-10 | yes |
| 1e-06 | 0.2 | 1.2e-08 | 1.2e-08 | 4.66e-08 | 1.2e-08 | yes |
| 1e-06 | 0.05 | 7.64e-07 | 7.72e-07 | 2.99e-06 | 7.69e-07 | yes |
| 1e-06 | 0.025 | 6e-06 | 6.23e-06 | 2.42e-05 | 6.15e-06 | yes |
| 1e-06 | 0.01 | 8.27e-05 | 0.000105 | 0.000406 | 9.62e-05 | yes |
| 1e-06 | 0.001 | 0.00409 | -0.0135 | -0.0524 | 0.0962 | n/a (s_eps >= 1: the lemma needs cot(pi mu/2) < kappa) |
| 1e-08 | 0.5 | 7.45e-13 | 7.45e-13 | 2.98e-12 | 7.69e-13 | yes |
| 1e-08 | 0.2 | 1.2e-11 | 1.2e-11 | 4.66e-11 | 1.2e-11 | yes |
| 1e-08 | 0.05 | 7.69e-10 | 7.69e-10 | 2.98e-09 | 7.69e-10 | yes |
| 1e-08 | 0.025 | 6.15e-09 | 6.15e-09 | 2.39e-08 | 6.15e-09 | yes |
| 1e-08 | 0.01 | 9.6e-08 | 9.62e-08 | 3.73e-07 | 9.62e-08 | yes |
| 1e-08 | 0.001 | 8.27e-05 | 0.000105 | 0.000406 | 9.62e-05 | yes |
| 1e-10 | 0.5 | 7.45e-16 | 7.45e-16 | 2.98e-15 | 7.69e-16 | yes |
| 1e-10 | 0.2 | 1.2e-14 | 1.2e-14 | 4.66e-14 | 1.2e-14 | yes |
| 1e-10 | 0.05 | 7.69e-13 | 7.69e-13 | 2.98e-12 | 7.69e-13 | yes |
| 1e-10 | 0.025 | 6.15e-12 | 6.15e-12 | 2.39e-11 | 6.15e-12 | yes |
| 1e-10 | 0.01 | 9.62e-11 | 9.62e-11 | 3.73e-10 | 9.62e-11 | yes |
| 1e-10 | 0.001 | 9.6e-08 | 9.62e-08 | 3.73e-07 | 9.62e-08 | yes |

## B. Ball arithmetic (python-flint, 400 bits)

| check | statement | cases | failed |
|---|---|---|---|
| B1 | dormant supremum s^p[mu^-q - (2-mu)^-q] is exact and below s^p mu^-q | 256 | 0 |
| B2 | distance form (ell p (p-2) d)^-q[1 - (x/(2-x))^q] | 256 | 0 |
| B3 | u' = ell p u^(p-1) g, v' = -ell p v^(p-1) g (arb_series derivative of the closed forms) | 256 | 0 |
| B4 | mass identity (rigorous quadrature) and Psi_p <= eps b mu^(-2/(p-2))/(2 ell p) | 256 | 0 |
| B5 | formation distance bracket d^0 (1 - r^q)^(1/q) < d_theta < d^0, d^0 = theta^(-1/q)/(ell p (p-2)) | 128 | 0 |
| B6 | depth 2: exact dormant value, distance form, the bracket of mu_theta and the exact d_theta | 192 | 0 |

Grid: p in 3..6, eps in {1e-4, 1e-6, 1e-8, 1e-12}, b in {1, 3/2}, ell in {1, 1/4}, mu in {1/2, 1/10, 1/100, 1/1000}, theta in {1/100, 1/1000}; inputs exact rationals; "failed" = a comparison not certified.

## C. The repository's network flows (moving centres): the dual is the integral of the drive

Design target (3 features, 21 kernels), plain flow of `factorized_rbf.make_rhs` with the 21 duals eta_j = int g_j dt integrated as extra states (LSODA, rtol 1e-10), horizon tau = 20; w_j(t) is inverted through W^-1 and compared with the integrated eta_j. "ratio" = max over the unformed units of (max_t |w_j|)/(the exact dormant supremum at the margin mu_j = 1 - max_t|eta_j|/wall actually reached): the theorem says = 1 up to the sampling of the time grid (121 points) and never above 1.

| readout | eps | wall | max |eta_int - eta_W^-1|/wall | formed units | smallest dormant margin | ratio min | ratio max | bound-form/exact (max) | max w / bound-form |
|---|---|---|---|---|---|---|---|---|---|
| D2 | 1e-06 | 3.45 | 1.36e-10 | 9 | 0.432 | 1 | 1 | 1 | 1 |
| D3 | 1e-06 | 333 | 2.02e-09 | 3 | 0.0676 | 1 | 1 | 9.55 | 1 |
| D4 | 0.0001 | 1.25e+03 | 5.46e-09 | 3 | 0.0658 | 1 | 1 | 14.1 | 0.999 |
| D5 | 0.01 | 66.7 | 8.54e-09 | 3 | 0.0463 | 1 | 1 | 16.8 | 0.998 |
| P0 | 1e-06 | 7.02e+03 | 2e-09 | 3 | 0.0608 | 1 | 1 | 36.2 | 0.258 |

(bound-form: (eps)^(p/2) mu^(-p/(p-2)) for depth p; eps^mu for depth 2; Lemma 2.2's mu^-3 form for the paper's unit.)

## D. Depth 1 (w = ell eta, no wall)

Design target, 21 kernels at fixed centres: sigma_min(Phi) = 0.00865, sigma_max = 5.13 (full column rank, so the flow converges to the least-squares solution); the LS solution has 19 of 21 coefficients above the formed threshold 0.009 (smallest |w_LS| = 0.0031, largest 0.869). The repository's depth-1 flow with lr_c = 0 against H^+(I - exp(-ell H t))b over t <= 3000: max abs difference 9.28e-11. At t = 3000 19 units are above the threshold. Formation-time certificate (Theorem 4): ||H|| = 0.432, ||g(0)|| = 0.23; units with g_j(0)^2 > 2 theta ||H|| ||g(0)||: 10, all formed by 2 theta/(ell |g_j(0)|): True.

Formation time of the eight units with the largest initial drive (first t with |w_j| > theta, against theta/(ell |g_j(0)|)):

| unit | g_j(0) | t formed (flow) | theta/(ell |g_j(0)|) | ratio |
|---|---|---|---|---|
| 5 | 0.0933 | 0.0082 | 0.00804 | 1.02 |
| 4 | 0.0842 | 0.0091 | 0.00891 | 1.02 |
| 16 | 0.0778 | 0.00985 | 0.00964 | 1.02 |
| 6 | 0.072 | 0.0107 | 0.0104 | 1.02 |
| 17 | 0.0681 | 0.0113 | 0.011 | 1.02 |
| 15 | 0.0611 | 0.0126 | 0.0123 | 1.02 |
| 11 | -0.0608 | 0.0126 | 0.0123 | 1.02 |
| 10 | -0.0567 | 0.0135 | 0.0132 | 1.02 |

## E. The sliding unit: first integral and overtaking margin

Model: eta' = G(c), c' = mu_c W(eta) G'(c) (frozen residual; depth 3, eps = 1e-6, mu_c = 12). First integral Psi*(eta) - Theta(c), Theta' = G/(mu_c G'). Gaussian landscape G(c) = 0.06 exp(-(c - 0.3)^2/(2 sigma^2)), sigma = 0.12, start 1.3 sigma from the maximum, margin x_a of the start:

| x_a | steps | drift of the first integral (relative) | |c - c0|/sigma | margin at the end | w at the end | G(c)/G(c0) |
|---|---|---|---|---|---|---|
| 0.2 | 28 | 5.56e-11 | 0.8 | 0.012 | 0.000574 | 2.05 |
| 0.05 | 18 | 1.71e-10 | 0.8 | 0.0117 | 0.000621 | 2.05 |
| 0.02 | 13 | 1.31e-10 | 0.8 | 0.0103 | 0.000908 | 2.05 |

Linear landscape G = g0 + h c (g0 = 0.02) from the start of the run (eta = 0, c = 0), followed until the frozen dual is at margin 0.02: the exact relation g^2 = g0^2 + 2 kappa Psi*(eta), kappa = mu_c h^2, and the lead (eta - g0 t)/wall against the first order in kappa, kappa eps/(6 g0^2)[1/(1-s) - 1/(1+s) - 2s], and against its leading form kappa a_ro eps/g0^2 (1/x - 1) (a_ro = 1/(6 ell)); the displacement is in widths of the Gaussian kernel (a linear landscape beyond ~1 width is not a landscape of the benchmark; the large-kappa rows show where the first order stops):

| h/g0 (1/length) | kappa/g0^2 | first-integral error (relative, max) | frozen margin | lead (ODE) | first order | leading form | displacement/sigma |
|---|---|---|---|---|---|---|---|
| 0.3 | 1.08 | 4.16e-14 | 0.02 | 8.56e-06 | 8.56e-06 | 8.82e-06 | 0.0125 |
| 1 | 12 | 2.17e-13 | 0.02 | 9.53e-05 | 9.51e-05 | 9.8e-05 | 0.0419 |
| 3 | 108 | 1.41e-12 | 0.02 | 0.000875 | 0.000856 | 0.000882 | 0.133 |
| 10 | 1.2e+03 | 3.31e-11 | 0.02 | 0.0156 | 0.00951 | 0.0098 | 3.01 |

## L. Fixed centres: the network against the inverse scale space of the weighted l1 norm (Corollary 5)

Design target, 21 kernels, centres frozen (lr_c = 0), time tau = t/wall, window [0, 25] (3 entries of the fixed-centre path: kernel 5 at 10.7, kernel 16 at 12.8, kernel 11 at 16.4). Reference: `mirror_free.iss_path` (equal walls). sup dual = max over units and 501 times of |u_net - u_ISS|, u_net = (dual from the coefficient through the closed form)/wall; L2 = rms over time and samples of the fitted function's difference over the rms of the target. Depth 1 has no dual in this sense (its flow is the l2 flow); it is run to the same window of tau at the paper's wall.

| readout | eps | sup dual error | sup/sqrt(eps) (depth 3, paper) or sup x ln(1/eps) (depth 2) | L2 fit error | units formed at the end |
|---|---|---|---|---|---|
| D3 | 0.0001 | 0.0137 | 1.37 | 0.0544 | 5 |
| D3 | 1e-06 | 0.00119 | 1.19 | 0.00963 | 3 |
| D3 | 1e-08 | 0.000119 | 1.19 | 2.35e-05 | 3 |
| D3 | 1e-10 | 1.19e-05 | 1.19 | 2.26e-08 | 3 |
| P0 | 0.0001 | 0.0677 | 6.77 | 0.159 | 5 |
| P0 | 1e-06 | 0.00677 | 6.77 | 0.0519 | 3 |
| P0 | 1e-08 | 0.000677 | 6.77 | 0.0241 | 3 |
| P0 | 1e-10 | 6.77e-05 | 6.77 | 3.79e-07 | 3 |
| D2 | 1e-06 | 0.133 | 1.84 | 0.15 | 9 |
| D2 | 1e-12 | 0.0717 | 1.98 | 0.0984 | 6 |
| D2 | 1e-24 | 0.0306 | 1.69 | 0.0652 | 5 |
| D2 | 1e-48 | 0.00999 | 1.10 | 0.0438 | 5 |
| D2 | 1e-96 | 0.00383 | 0.85 | 0.0299 | 5 |
| D1 | 1 | n/a |  | 0.713 | 19 |

D3: slope of log sup dual error on log eps: 0.509 (1/2 for a sqrt(eps) rate; for depth 2 the rate is 1/ln(1/eps), slope ~ 0).

P0: slope of log sup dual error on log eps: 0.500 (1/2 for a sqrt(eps) rate; for depth 2 the rate is 1/ln(1/eps), slope ~ 0).

D2: slope of log sup dual error on log eps: 0.016 (1/2 for a sqrt(eps) rate; for depth 2 the rate is 1/ln(1/eps), slope ~ 0).

## K. A one-parameter family of ties (one feature between two kernels): the margin against the critical eps at fixed geometry

Feature of amplitude 0.8 at 0.05 + delta (between the kernels at 0 and 0.1 of the 21), delta small; the path enters one kernel and the other is the runner-up. For each delta the critical eps is found by bisection on log10 eps (flip: the formed set or the entry order differs from the path's); "model eps*" is the overtaking model's. A pure sqrt(eps) law at fixed geometry means a slope 1/2 of log margin on log eps*.

D3:

| delta | margin | eps* bracket (from bisection) | geometric mean | model eps* | model/obs |
|---|---|---|---|---|---|
| 0.0003 | 0.00104 | 7.5e-08 - 8.06e-08 | 7.77e-08 | 1.78e-07 | 2.29 |
| 0.0006 | 0.00208 | 2.94e-07 - 3.16e-07 | 3.05e-07 | 7.03e-07 | 2.3 |
| 0.0012 | 0.00416 | 1.15e-06 - 1.24e-06 | 1.2e-06 | 2.75e-06 | 2.29 |
| 0.0025 | 0.00864 | 4.87e-06 - 5.23e-06 | 5.05e-06 | 1.13e-05 | 2.25 |
| 0.005 | 0.0172 | 1.78e-05 - 1.91e-05 | 1.84e-05 | 4.13e-05 | 2.24 |
| 0.01 | 0.0341 | 6.04e-05 - 6.49e-05 | 6.26e-05 | 0.000139 | 2.22 |
| 0.02 | 0.0671 | 0.000191 - 0.000205 | 0.000198 | 0.000408 | 2.06 |
| 0.04 | 0.13 | 0.000562 - 0.000604 | 0.000583 | 0.000989 | 1.7 |

Slope of log10 margin on log10 eps* (observed): 0.537 over 8 members; the same slope with the model's eps*: 0.550. (1/2: pure sqrt(eps).)

P0:

| delta | margin | eps* bracket (from bisection) | geometric mean | model eps* | model/obs |
|---|---|---|---|---|---|
| 0.0006 | 0.00208 | 1.04e-07 - 1.18e-07 | 1.1e-07 | 3.47e-07 | 3.14 |
| 0.0012 | 0.00416 | 4.14e-07 - 4.7e-07 | 4.41e-07 | 1.36e-06 | 3.07 |
| 0.0025 | 0.00864 | 1.65e-06 - 1.88e-06 | 1.76e-06 | 5.6e-06 | 3.17 |
| 0.005 | 0.0172 | 5.83e-06 - 6.61e-06 | 6.21e-06 | 2.04e-05 | 3.28 |
| 0.01 | 0.0341 | 2.05e-05 - 2.33e-05 | 2.19e-05 | 6.85e-05 | 3.13 |
| 0.02 | 0.0671 | 7.23e-05 - 8.2e-05 | 7.7e-05 | 0.000202 | 2.62 |

Slope of log10 margin on log10 eps* (observed): 0.533 over 6 members; the same slope with the model's eps*: 0.542. (1/2: pure sqrt(eps).)

## P. Sliding paths of the benchmark population

1034 one-dimensional targets of `exp_transit_certificate` (families A, B, D), equal masses; per target the stagewise drives and centre-slopes of the dictionary at the starting centres under each block's residual, the duals at the stages' ends, the runner-up margins. Failed or with exits: 0; median stage count 3.

## F. The 200 targets of Part C: what the theorem predicts, what the overtaking model adds

Matched to the population records: 200 of 200. For each cell: *scalar*: the target is predicted to mismatch iff its never-entering margin is below the exact formation margin mu_theta(eps) (theta = the target's formed threshold; depth 3: d_theta/wall, the paper's unit: its Lemma 2.1 distance); *overtaking*: the model of Section E applied to every (stage, competitor) of the path with mu_c = 12 and the readout's constant a_ro (no fitted parameter): a flip is predicted iff some competitor's overtaking margin reaches its margin at the stage's end. 'Mismatch' in the first column is the paper's (the formed set differs from the path's); 'flip' also counts an entry order different from the path's. For the paper's unit mu_theta is negative when theta exceeds the coefficient that the unit carries at its wall (0.0079): such a unit is formed only after its dual has passed the wall, and the scalar prediction is empty for those targets.

| readout | eps | set mismatches | flips | mu_theta (range over targets) | scalar: predicted | scalar vs set mismatch | overtaking: predicted | overtaking vs flip | overtaking vs set mismatch |
|---|---|---|---|---|---|---|---|---|---|
| D3 | 1e-06 | 24 | 25 | 0.00215 - 0.0299 | 18 | TP 12, FP 6, FN 12, TN 170 | 40 | TP 24, FP 16, FN 1, TN 159 | TP 23, FP 17, FN 1, TN 159 |
| D3 | 1e-08 | 5 | 6 | 0.000215 - 0.00299 | 5 | TP 3, FP 2, FN 2, TN 193 | 10 | TP 5, FP 5, FN 1, TN 189 | TP 4, FP 6, FN 1, TN 189 |
| P0 | 1e-06 | 34 | 36 | -0.00365 - 0.0133 | 2 | TP 2, FP 0, FN 32, TN 166 | 46 | TP 31, FP 15, FN 5, TN 149 | TP 29, FP 17, FN 5, TN 149 |
| P0 | 1e-08 | 9 | 10 | -0.000365 - 0.00133 | 1 | TP 1, FP 0, FN 8, TN 191 | 13 | TP 9, FP 4, FN 1, TN 186 | TP 8, FP 5, FN 1, TN 186 |

D3, eps = 1e-06: the scalar theorem predicts a mismatch for 18 targets, whose never-entering margins are 0.00008, 0.00022, 0.00026, 0.00028, 0.00067, 0.00072, 0.00104, 0.00133, 0.00300, 0.00341, 0.00353, 0.00359, 0.00407, 0.00409, 0.00412, 0.00430, 0.00590, 0.00682; the 24 observed set mismatches have margins 0.00008, 0.00022, 0.00026, 0.00028, 0.00067, 0.00104, 0.00133, 0.00300, 0.00359, 0.00412, 0.00517, 0.00568, 0.00590, 0.00682, 0.00929, 0.00950, 0.01055, 0.01116, 0.01243, 0.01459, 0.01479, 0.01621, 0.01806, 0.01853.

D3, eps = 1e-08: the scalar theorem predicts a mismatch for 5 targets, whose never-entering margins are 0.00008, 0.00022, 0.00026, 0.00028, 0.00067; the 5 observed set mismatches have margins 0.00022, 0.00026, 0.00067, 0.00133, 0.00300.

P0, eps = 1e-06: the scalar theorem predicts a mismatch for 2 targets, whose never-entering margins are 0.00067, 0.00590; the 34 observed set mismatches have margins 0.00008, 0.00022, 0.00026, 0.00028, 0.00067, 0.00072, 0.00104, 0.00133, 0.00300, 0.00353, 0.00359, 0.00412, 0.00517, 0.00568, 0.00587, 0.00590, 0.00682, 0.00712, 0.00727, 0.00754, 0.00824, 0.00929, 0.00950, 0.01055, 0.01116, 0.01243, 0.01259, 0.01329, 0.01459, 0.01479, 0.01621, 0.01720, 0.01806, 0.01853.

P0, eps = 1e-08: the scalar theorem predicts a mismatch for 1 targets, whose never-entering margins are 0.00067; the 9 observed set mismatches have margins 0.00008, 0.00022, 0.00026, 0.00028, 0.00067, 0.00133, 0.00300, 0.00517, 0.00568.

Largest runner-up margin among the mismatches (the statistic of the paper), observed and predicted by the overtaking model, and the exponent between 1e-6 and 1e-8:

| readout | max margin obs 1e-6 | obs 1e-8 | exponent obs | max margin pred 1e-6 | pred 1e-8 | exponent pred | mismatches 1e-6 | 1e-8 | count exponent |
|---|---|---|---|---|---|---|---|---|---|
| D3 | 0.0185 | 0.0030 | 0.395 | 0.0185 | 0.0030 | 0.395 | 24 | 5 | 0.341 |
| P0 | 0.0185 | 0.0057 | 0.257 | 0.0226 | 0.0052 | 0.320 | 34 | 9 | 0.289 |

The theorem-level bound (a) on the targets: (max dormant coefficient of the run)/(the exact dormant supremum at the path's never-entering margin), over the targets with mu_never >= 100 sqrt(eps) whose run formed the path's number of units (the actual dual differs from the path's by O(sqrt(eps)), so the ratio is 1 up to ~3 dev sqrt(eps)/mu):

| readout | eps | targets | median ratio | 10th pct | 90th pct | max |
|---|---|---|---|---|---|---|
| D3 | 1e-06 | 40 | 0.975 | 0.941 | 0.985 | 0.99 |
| D3 | 1e-08 | 171 | 0.989 | 0.929 | 0.997 | 1.19 |
| P0 | 1e-06 | 40 | 1.13 | 1.06 | 1.2 | 1.23 |
| P0 | 1e-08 | 171 | 1.02 | 0.975 | 1.1 | 2.04 |

Depth 2 against the exact formation margin mu_theta(eps) = ln(2/(theta + sqrt(theta^2 + 4 eps^2)))/ln(1/eps) (theta per target) and the paper's rule ln(100)/ln(1/eps): predicted outcome "formed set = path's" iff mu_never > the margin, against the stored D2 runs:

| eps | D2 passes (observed) | exact margin | rule |
|---|---|---|---|
| 1e-06 | 0 | TP 0, FP 0, FN 0, TN 200 | TP 0, FP 0, FN 0, TN 200 |
| 1e-12 | 2 | TP 2, FP 7, FN 0, TN 191 | TP 2, FP 1, FN 0, TN 197 |

Part D of the stored runs (D2 on the design and 20 targets at eps = 1e-6 ... 1e-144): outcome predicted correctly by the exact margin in 122 of 126 runs, by the rule in 124.

The design (margin 0.0664, theta = 0.009): the scalar theorem's "path reached iff mu > mu_theta(readout, eps)" against the stored runs (set equal to the path's):

| readout | eps | mu_theta | predicted pass | observed pass | formed units | agree |
|---|---|---|---|---|---|---|
| P0 | 0.01 | -0.0453 | True | False | 15 | NO |
| P0 | 0.0001 | -0.00193 | True | False | 4 | NO |
| P0 | 1e-06 | -0.000191 | True | True | 3 | yes |
| P0 | 1e-08 | -1.91e-05 | True | True | 3 | yes |
| P0 | 1e-10 | -1.91e-06 | True | True | 3 | yes |
| P0 | 0.01 | -0.0453 | True | False | 6 | NO |
| P0 | 0.0001 | -0.00193 | True | True | 3 | yes |
| P0 | 1e-06 | -0.000191 | True | True | 3 | yes |
| P0 | 1e-08 | -1.91e-05 | True | True | 3 | yes |
| P0 | 1e-10 | -1.91e-06 | True | True | 3 | yes |
| D3 | 0.01 | 0.476 | False | False | 12 | yes |
| D3 | 0.0001 | 0.0481 | True | True | 3 | yes |
| D3 | 1e-06 | 0.00481 | True | True | 3 | yes |
| D3 | 1e-08 | 0.000481 | True | True | 3 | yes |
| D3 | 1e-10 | 4.81e-05 | True | True | 3 | yes |
| D4 | 0.01 | 0.105 | False | False | 6 | yes |
| D4 | 0.0001 | 0.00105 | True | True | 3 | yes |
| D4 | 1e-06 | 1.05e-05 | True | True | 3 | yes |
| D4 | 1e-08 | 1.05e-07 | True | True | 3 | yes |
| D4 | 1e-10 | 1.05e-09 | True | True | 3 | yes |
| D5 | 0.01 | 0.0169 | True | True | 3 | yes |
| D5 | 0.0001 | 1.69e-05 | True | True | 3 | yes |
| D5 | 1e-06 | 1.69e-08 | True | True | 3 | yes |
| D5 | 1e-08 | 1.69e-11 | True | True | 3 | yes |
| D5 | 1e-10 | 1.69e-14 | True | False | 1 | yes |
| D2 | 0.01 | 0.905 | False | False | 18 | yes |
| D2 | 0.0001 | 0.511 | False | False | 12 | yes |
| D2 | 1e-06 | 0.341 | False | False | 9 | yes |
| D2 | 1e-08 | 0.256 | False | False | 8 | yes |
| D2 | 1e-10 | 0.205 | False | False | 6 | yes |

The displacement of the dormant centres on the design (stored runs; largest over the unformed units, in kernel widths) against the mass identity: by Theorem 2(iii) it is (centre rate x drive-slope/drive) x Psi_p((1 - mu) wall) with the same target-dependent factor for every depth and eps, Psi_p = eps b mu^(-2/(p-2))/(2 ell p) at the leading order; the ratio disp/Psi_p over the 8 runs (D3, D4, D5 at eps = 1e-6, 1e-8, 1e-10; the failed run is excluded) has median 245 and spread (max/min) 1.06:

| readout | eps | dormant displacement | Psi_p at the design margin | ratio |
|---|---|---|---|---|
| D3 | 1e-06 | 0.00875 | 3.75e-05 | 233 |
| D3 | 1e-08 | 9.11e-05 | 3.75e-07 | 243 |
| D3 | 1e-10 | 9.15e-07 | 3.75e-09 | 244 |
| D4 | 1e-06 | 0.000416 | 1.7e-06 | 245 |
| D4 | 1e-08 | 4.17e-06 | 1.7e-08 | 245 |
| D4 | 1e-10 | 4.17e-08 | 1.7e-10 | 245 |
| D5 | 1e-06 | 0.000117 | 4.74e-07 | 246 |
| D5 | 1e-08 | 1.17e-06 | 4.74e-09 | 246 |

## G. Per-target critical eps (scan of `theory3d_scan.py`) against the overtaking model

D3: 133 targets scanned (grid 1e-11 ... 0.0001, targets of margin below 0.03 and 15 of larger margin from an earlier pass); 42 do not flip even at the largest eps of the grid (critical eps above it), 0 already flip at the smallest (critical eps below it; censored); 91 with a bracketed critical eps. log10(predicted/observed critical eps): median 0.15, quartiles -0.21 / 0.32; Spearman correlation of the two critical eps 0.767; slope of log10(pred) on log10(obs) 0.87; slope of log10(margin) on log10(observed critical eps) 0.24 (1/2 would be a pure sqrt(eps) law, 1/4 the paper's pooled exponent). 

Decomposition: log10 margin = log10 c + (1/2) log10 eps*, c = margin/sqrt(eps*) the observed constant of the target. Var(log10 margin) = 0.405, Var(log10 c) = 0.425, correlation 0.45: the regression slope of log margin on log eps* is (Var a - Cov(a, b))/(2 (Var a + Var b - 2 Cov(a, b))) = 0.239 (identity; if log c were independent of the margin it would be (1/2) Var a/(Var a + Var b) = 0.244): the spread of c (a factor 4.5 per standard deviation) attenuates the slope from 1/2 to the observed value.

P0: 55 targets scanned (grid 1e-09 ... 1e-05, targets of margin below 0.01); 26 do not flip even at the largest eps of the grid (critical eps above it), 4 already flip at the smallest (critical eps below it; censored); 25 with a bracketed critical eps. log10(predicted/observed critical eps): median 0.22, quartiles -0.18 / 0.43; Spearman correlation of the two critical eps 0.489; slope of log10(pred) on log10(obs) 0.45; slope of log10(margin) on log10(observed critical eps) 0.10 (1/2 would be a pure sqrt(eps) law, 1/4 the paper's pooled exponent). 

Decomposition: log10 margin = log10 c + (1/2) log10 eps*, c = margin/sqrt(eps*) the observed constant of the target. Var(log10 margin) = 0.583, Var(log10 c) = 0.741, correlation 0.81: the regression slope of log margin on log eps* is (Var a - Cov(a, b))/(2 (Var a + Var b - 2 Cov(a, b))) = 0.098 (identity; if log c were independent of the margin it would be (1/2) Var a/(Var a + Var b) = 0.220): the spread of c (a factor 7.3 per standard deviation) attenuates the slope from 1/2 to the observed value.

The columns of the table below are, per scanned target: the runner-up margin, the observed bracket (flip at eps_hi, none at eps_lo), the model's critical eps and the ratio of the margin to sqrt(observed eps).

D3, the 15 flipping targets of largest margin:

| target | margin | eps_hi (flip) | eps_lo (none) | model eps* | margin/sqrt(eps_hi) |
|---|---|---|---|---|---|
| 59 | 0.0437 | 1e-05 | 3.16e-06 | 6.01e-06 | 13.8 |
| 6 | 0.0435 | 0.0001 | 3.16e-05 | 0.000195 | 4.3 |
| 11 | 0.0368 | 0.0001 | 3.16e-05 | 0.000169 | 3.7 |
| 47 | 0.0359 | 3.16e-05 | 1e-05 | 5.45e-06 | 6.4 |
| 45 | 0.0355 | 0.0001 | 3.16e-05 | 0.000153 | 3.5 |
| 32 | 0.0350 | 0.0001 | 3.16e-05 | 0.000148 | 3.5 |
| 31 | 0.0324 | 0.0001 | 3.16e-05 | 8.14e-05 | 3.2 |
| 53 | 0.0311 | 0.0001 | 3.16e-05 | 4.85e-05 | 3.1 |
| 191 | 0.0300 | 0.0001 | 3.16e-05 | 1.02e-05 | 3.0 |
| 29 | 0.0295 | 0.0001 | 3.16e-05 | 0.000118 | 2.9 |
| 148 | 0.0286 | 0.0001 | 3.16e-05 | 9.28e-05 | 2.9 |
| 30 | 0.0283 | 0.0001 | 3.16e-05 | 0.000114 | 2.8 |
| 40 | 0.0283 | 0.0001 | 3.16e-05 | 0.000222 | 2.8 |
| 82 | 0.0282 | 1e-05 | 3.16e-06 | 8.46e-06 | 8.9 |
| 85 | 0.0274 | 1e-05 | 3.16e-06 | 3.31e-05 | 8.7 |

P0, the 15 flipping targets of largest margin:

| target | margin | eps_hi (flip) | eps_lo (none) | model eps* | margin/sqrt(eps_hi) |
|---|---|---|---|---|---|
| 50 | 0.0095 | 1e-07 | 1e-08 | 8.48e-08 | 30.0 |
| 0 | 0.0094 | 1e-05 | 1e-06 | 5.47e-06 | 3.0 |
| 57 | 0.0093 | 1e-07 | 1e-08 | 5.23e-08 | 29.4 |
| 36 | 0.0084 | 1e-06 | 1e-07 | 1.37e-07 | 8.4 |
| 152 | 0.0082 | 1e-06 | 1e-07 | 1.73e-06 | 8.2 |
| 67 | 0.0075 | 1e-06 | 1e-07 | 1.3e-06 | 7.5 |
| 162 | 0.0073 | 1e-06 | 1e-07 | 2.94e-06 | 7.3 |
| 73 | 0.0068 | 1e-06 | 1e-07 | 3.2e-07 | 6.8 |
| 74 | 0.0059 | 1e-06 | 1e-07 | 2.08e-07 | 5.9 |
| 151 | 0.0059 | 1e-06 | 1e-07 | 8.16e-07 | 5.9 |
| 7 | 0.0057 | 1e-08 | 1e-09 | 1.28e-08 | 56.8 |
| 168 | 0.0052 | 1e-05 | 1e-06 | 2.03e-06 | 1.6 |
| 102 | 0.0052 | 1e-08 | 1e-09 | 9.66e-09 | 51.7 |
| 80 | 0.0047 | 1e-06 | 1e-07 | 2.66e-07 | 4.7 |
| 101 | 0.0041 | 1e-06 | 1e-07 | 2.41e-07 | 4.1 |

## H. The population: counts and the largest flipped margin of the model as functions of eps, and the exponent

1034 of 1034 population targets with a usable path (no exit, no solver failure). The model's count of flips N(eps) and the largest runner-up margin among them, for the two readouts; in the model N(eps) = #(targets with rho_t < sqrt(a_ro eps)) to leading order, rho_t = min over (stage, competitor) of x/sqrt(q), so the local exponent of N in eps is half the local exponent of the distribution function of rho_t.

D3:

| eps | flips N(eps) | largest margin among flips | its value / sqrt(eps) | 95th pct margin |
|---|---|---|---|---|
| 0.001 | 1007 | 0.1429 | 4.5 | 0.1040 |
| 0.000316 | 841 | 0.1334 | 7.5 | 0.0728 |
| 0.0001 | 668 | 0.1132 | 11.3 | 0.0576 |
| 3.16e-05 | 484 | 0.0725 | 12.9 | 0.0442 |
| 1e-05 | 340 | 0.0725 | 22.9 | 0.0361 |
| 3.16e-06 | 225 | 0.0335 | 18.8 | 0.0213 |
| 1e-06 | 156 | 0.0253 | 25.3 | 0.0171 |
| 3.16e-07 | 104 | 0.0162 | 28.8 | 0.0109 |
| 1e-07 | 71 | 0.0124 | 39.3 | 0.0067 |
| 3.16e-08 | 54 | 0.0112 | 63.1 | 0.0060 |
| 1e-08 | 35 | 0.0112 | 112.2 | 0.0045 |
| 3.16e-09 | 22 | 0.0112 | 199.5 | 0.0063 |
| 1e-09 | 18 | 0.0112 | 354.8 | 0.0073 |
| 3.16e-10 | 12 | 0.0066 | 368.4 | 0.0031 |
| 1e-10 | 12 | 0.0066 | 655.1 | 0.0031 |
| 3.16e-11 | 8 | 0.0001 | 14.6 | 0.0001 |
| 1e-11 | 7 | 0.0001 | 16.2 | 0.0000 |
| 3.16e-12 | 6 | 0.0000 | 3.1 | 0.0000 |
| 1e-12 | 0 | n/a | n/a | n/a |

Bootstrap (2000 samples of 200 targets from the population, the model's flips at 1e-6 and 1e-8, 1997 samples with flips at both): exponent between 1e-6 and 1e-8 of the count N(eps): 5/25/50/75/95th percentiles 0.22, 0.28, 0.33, 0.38, 0.49; of the largest margin among the flips: 0.09, 0.25, 0.38, 0.46, 0.71.

The constant of the flip margin c = sqrt(a_ro mu_c) |h/g| (a competitor flips when its margin is below c sqrt(eps)), over the 4316 (stage, competitor) pairs with margin below 0.1: 10/50/90/99th percentiles 1.3, 3.8, 48.9, 95.3; fraction above 30: 0.256.

Distribution function of rho_t = min over (stage, competitor) of x |g|/(sqrt(mu_c)|h|) (the same for both readouts; the flip condition is rho_t < sqrt(a_ro eps)): count of targets below, local index between consecutive decades (1 would mean a pure sqrt(eps) law for the counts; the first row is saturated: every target has rho_t below 0.01):

| rho range | count below the larger | count below the smaller | local index |
|---|---|---|---|
| 0.01 -> 0.001 | 1034 | 724 | 0.15 |
| 0.001 -> 0.0001 | 724 | 137 | 0.72 |
| 0.0001 -> 1e-05 | 137 | 24 | 0.76 |

P0:

| eps | flips N(eps) | largest margin among flips | its value / sqrt(eps) | 95th pct margin |
|---|---|---|---|---|
| 0.001 | 1034 | 0.1592 | 5.0 | 0.1104 |
| 0.000316 | 954 | 0.1334 | 7.5 | 0.0912 |
| 0.0001 | 765 | 0.1334 | 13.3 | 0.0626 |
| 3.16e-05 | 583 | 0.1101 | 19.6 | 0.0501 |
| 1e-05 | 421 | 0.0725 | 22.9 | 0.0430 |
| 3.16e-06 | 290 | 0.0437 | 24.6 | 0.0326 |
| 1e-06 | 195 | 0.0297 | 29.7 | 0.0208 |
| 3.16e-07 | 130 | 0.0253 | 45.0 | 0.0147 |
| 1e-07 | 88 | 0.0124 | 39.3 | 0.0097 |
| 3.16e-08 | 63 | 0.0112 | 63.1 | 0.0065 |
| 1e-08 | 47 | 0.0112 | 112.2 | 0.0052 |
| 3.16e-09 | 27 | 0.0112 | 199.5 | 0.0056 |
| 1e-09 | 21 | 0.0112 | 354.8 | 0.0066 |
| 3.16e-10 | 15 | 0.0112 | 631.0 | 0.0080 |
| 1e-10 | 12 | 0.0066 | 655.1 | 0.0031 |
| 3.16e-11 | 10 | 0.0066 | 1165.0 | 0.0037 |
| 1e-11 | 8 | 0.0001 | 26.0 | 0.0001 |
| 3.16e-12 | 7 | 0.0001 | 28.8 | 0.0000 |
| 1e-12 | 6 | 0.0000 | 5.5 | 0.0000 |

Bootstrap (2000 samples of 200 targets from the population, the model's flips at 1e-6 and 1e-8, 1999 samples with flips at both): exponent between 1e-6 and 1e-8 of the count N(eps): 5/25/50/75/95th percentiles 0.22, 0.27, 0.31, 0.36, 0.44; of the largest margin among the flips: 0.15, 0.27, 0.34, 0.41, 0.55.

The constant of the flip margin c = sqrt(a_ro mu_c) |h/g| (a competitor flips when its margin is below c sqrt(eps)), over the 4316 (stage, competitor) pairs with margin below 0.1: 10/50/90/99th percentiles 1.8, 5.4, 69.6, 135.7; fraction above 30: 0.337.

## Status of the statements (the note has the proofs)

* PROVED (exact, scalar): the closed forms of the unit of depth p (Theorem 1); the dormant coefficient bound |w| <= (eps b)^(p/2) mu^(-p/(p-2)), sharp, and the same in terms of a drive bound (Theorem 2); the distance form and the formation criterion "formed iff the dual is within d_theta of the wall", d_theta ~ theta^(-(p-2)/p)/(ell p (p-2)) (Theorem 2); depth 2: w = 2 eps b sinh(4 ell eta), formation margin exactly ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/ln(1/(eps b)) and the bracket by ln(1/theta)/ln(1/(eps b)) (Theorem 3); depth 1: w = ell eta, the fixed-centre flow and its formation times (Theorem 4). Verified here: A (ODE, 1e-10), B (ball arithmetic), C (the repository's flows with moving centres), D.
* PROVED up to the cited abstract steps (those of the paper's fixed-centre theorem): the fixed-centre limits (Corollary 5): weighted l1 inverse scale space for depth p >= 3, unweighted l1 for depth 2 (at a logarithmic rate), the l2 flow for depth 1; the epi-convergence of the potentials is proved. Verified in L.
* PROVED (exact, model): the first integral of the frozen-environment sliding unit (E) and the proportionality of the dormant centre's displacement to the mass identity Psi_p (checked on the stored runs of D3, D4, D5: ratio constant to 6%).
* COMPUTER-ASSISTED (a first-order model without a fitted parameter): the competitor's overtaking margin c sqrt(eps), c = sqrt(mu_c a_ro)|h/g|, its agreement with the 200 targets (F), with the per-target critical eps (G) and with the tie family (K, slope 0.53 at fixed geometry), and the explanation of the pooled exponent (the attenuation of the slope by the spread of c; the bootstrap of H). The explanation is not a theorem.
* NOT proved: that the actual duals of the network follow the path's to the accuracy that the thresholds need (the hypothesis of the dormant bound is on the actual dual); the sliding-limit theorem for depth p != the paper's unit; the layer analysis for D_p; the rate of depth 2's dual error (observed to fall faster than 1/ln(1/eps) on the design). The scalar theorem alone does not predict the paper's mismatch margins (it accounts for 12 of the 24 depth-3 set mismatches at 1e-6 and 2 of the 34 of the paper's unit); the dual deviation does.
