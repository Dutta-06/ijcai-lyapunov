# Theory round 3, target F: optimizers (Adam / sign descent) and two-layer small starts

Code: `experiments/theory3f_escape.py`; tests `tests/test_theory3f_escape.py`; report `results/synthetic/theory3f.md`.
Status words: PROVED (argument below, each step checked and run numerically), COMPUTER-ASSISTED, OUTLINE, FAILED.

## 0. Setting and notation

Dictionary Gram matrix `G` (PSD), target correlations `c`, `F(w) = 1/2 (w-w*)^T G (w-w*)`, drive of unit i: `g_i(w) = -dF/dw_i = c_i - (Gw)_i`.
Unit i: `w_i = W(p_i)`, `p_i in R^q`, one map W for all units, equal starts `p_i(0) = p_0` (the equal-mass case of the paper), common step sizes.
D3 unit (repo `Cubic3`): `p = (u, v)`, `W = u^3 - v^3`, `p_0 = (s, s)`, `s^2 = eps b_i`.  `-dF/du = 3u^2 g`, `-dF/dv = -3 v^2 g`.
The paper's unit has the same structure in the escape phase (alpha moves at speed nu^2 g, nu then grows at the rate alpha g; both gradients carry the factor g).
`theta` = a coefficient threshold (the "formed" level); `T_i(theta)` = first step with `|w_i| >= theta`.
`||G_i||_1` = l1 norm of row i of G.  Conventions for Adam: `m_t = b1 m_{t-1} + (1-b1) g_t`, `v_t = b2 v_{t-1} + (1-b2) g_t^2`, hats = division by `1 - b^t`,
update `lr * mhat/(sqrt(vhat)+eps_a)`.  Sign descent = Adam with `b1 = b2 = 0`, `eps_a -> 0`.

Reference (plain gradient flow, frozen drive): `u' = 3 lr g u^2` for `w = u^3`, `u(t) = 1/(1/s - 3 lr g t)`, so the time to reach `u_c` from `s` is
`T^GD = (1/s - 1/u_c)/(3 lr g)`.  (Exact; checked by Euler in the report, part A.)  The ratio of the times of two units is `g_1/g_2`, the gap is
`(1/s - 1/u_c)(1/g_2 - 1/g_1)/(3 lr) -> infinity` as `s -> 0`: this is the ordering/selection mechanism.

## 1. F1a  Sign descent: simultaneous escape (PROVED)

**Lemma 1 (envelope, D3).** Under sign descent each D3 unit has `u(t) = s + lr x_t`, `v(t) = s - lr x_t`, with `x_t = sum_{k<=t} sigma_k`, `sigma_k in {-1,0,+1}`.
Proof: `-dF/du = 3u^2 g` has the sign of `g` (or is 0 if u = 0), `-dF/dv = -3v^2 g` the opposite sign (or 0).  Hence u and v move by `+-lr` in opposite
directions (or stand still if one is exactly 0: a measure-zero event for `s` not a multiple of `lr`; for such coincidences the statement below is
a statement about generic starts and the report uses `s` that are not multiples of the step).  So `w = u^3 - v^3 = 6 s^2 (lr x) + 2 (lr x)^3` and `|w|` is at most
`omega(t) := 6 s^2 lr t + 2 (lr t)^3`, an increasing function of t (`|x_t| <= t`).  Therefore `T_i >= T_env := min{t : omega(t) >= theta}` for every unit.

**Theorem 1.** Let `theta > 0` and let unit i satisfy `|g_i(0)| > theta ||G_i||_1`.  Then `T_i(theta) = T_env(theta, s, lr)` exactly.
*Proof.* For `t < T_env` every unit k has `|w_k(t)| < theta` (Lemma 1), so `|g_i(t) - g_i(0)| = |sum_j G_ij w_j(t)| <= theta ||G_i||_1 < |g_i(0)|` and `g_i(t)` has the sign
of `g_i(0)` for all `t < T_env`.  Then `sigma_k = sign g_i(0)` for all steps up to `T_env`: `x_t = +-t`, `|w_i(t)| = omega(t)`, so `|w_i|` first reaches theta at `T_env`. QED.
Corollaries.  (a) Two units with drives `g_1 > g_2` (both above the hypothesis) cross at the same step, whatever `g_1/g_2`; `T^sign_2 - T^sign_1 = 0` against `T^GD_2 - T^GD_1 -> infinity`.
(b) `T_env` depends on `s` only through `3 s^2 lr t` relative to `(lr t)^3`: with `a_0(s)` the solution of `a^3 + 3 s^2 a = theta/2`,
`T_env = ceil(a_0(s)/lr)`, `0 <= a_0(0) - a_0(s) <= s^2 a_0(0)/a_0(s)^2`, so `T_env(s) - T_env(0) in [-s^2/(a_0 lr) - 1, 0]`: for `s^2 <= a_0 lr` the escape time is independent
of the start up to 1 step.  In the paper's parametrisation (`s^2 = eps b_i`): the escape does not depend on the mass `eps` or on the mass ordering (different `b_i` give steps differing by at most `eps |b_i-b_j|/(a_0 lr) + 1`).
(c) General readouts: the first sentence of the proof of Theorem 1 holds for any map `W` and any equal starts (the path of unit i is determined by `sign(g_i) * sign(grad W)` at each step), so all units whose drive keeps its sign
have identical parameter paths; the paper's unit (alpha, nu) is an example (numerics, report part E).
Numerical check: report part A (six drives spanning 1e5: sign descent window 1.000, Adam 1.000, GD 1e5), part E (the repo's design: sign descent forms the 21 units at the same step, 207 / 1650 / 488 / 4880).

## 2. F1b  Adam (PROVED)

**Prop 2a (scale invariance).** For `eps_a = 0` (with 0/0 := 0) the update `U_t(g_1..g_t) = mhat_t/sqrt(vhat_t)` satisfies `U_t(lambda g) = U_t(g)` for `lambda > 0`.
Hence if the gradient of coordinate i is `lambda_i h(p_i(t))` with a constant `lambda_i > 0` (a frozen drive) then `p_i(t)` does not depend on `lambda_i`: units with different drives have identical trajectories (induction on t).
Numerically (report B): steps to `u = 0.1` are 511 for all drive scales from 1e-8 to 1e8 at `eps_a = 1e-30`.
In the true flow the drive is frozen up to `(1 +- delta)`, `delta = theta ||G_i||_1/|g_i(0)|` (Theorem 1's estimate).

**Prop 2b (window).** One-sided unit `w = u^3` (or the u-coordinate of D3), `eps_a = 0`, gradient of u: `G_t = 3 u_t^2 D_t` with `D_t in [(1-delta) g, (1+delta) g]`, `g > 0`.
(i) For any real gradient sequence `|R_t| <= C_t`, `R_t = mhat_t/sqrt(vhat_t)`, `C_t^2 = (1-b1)^2 (1-b2^t)/((1-b2)(1-b1^t)^2) sum_{j<t} (b1^2/b2)^j`.
*Proof.* Cauchy-Schwarz: `mhat_t = sum_k a_k G_k = sum_k (a_k/sqrt(b_k)) (sqrt(b_k) G_k) <= sqrt(sum a_k^2/b_k) sqrt(sum b_k G_k^2)`, with `a_k = (1-b1) b1^{t-k}/(1-b1^t)`, `b_k = (1-b2) b2^{t-k}/(1-b2^t)` (sum of b_k G_k^2 = vhat_t).  Evaluate the geometric sum. `C = sup_t C_t = 7.27` for `(0.9, 0.999)` (computed in `adam_ratio_bounds`; t* = 37640).
(ii) If the gradient is quasi-monotone, `G_t >= (1-eta) max_{k<=t} G_k`, then `R_t >= (1-b1)(1-eta)`.
*Proof.* `mhat_t >= a_t G_t >= (1-b1) G_t` (all terms nonnegative; `a_t = (1-b1)/(1-b1^t)`), `sqrt(vhat_t) <= max_k G_k <= G_t/(1-eta)`.
In the model above u is increasing (R > 0), so `G_t >= 3 u_t^2 (1-delta) g >= (1-delta)/(1+delta) max_k G_k`: `eta = 2 delta/(1+delta)`.  In the diagonal model (`D = g - u^3`) `u^2 (g - u^3)` is nondecreasing for `u^3 <= 2g/5`, i.e. eta = 0.
(iii) Consequently `u(t) - s in [lr t (1-b1)(1-eta), lr t C]` and for every unit of this kind `T_i(u_c) in [(u_c - s)/(lr C), (u_c - s)/(lr (1-b1)(1-eta))]`: two units differ by a factor of at most `C/((1-b1)(1-eta)) = 72.7/(1-eta)` for the defaults, independent of `g_1/g_2`.
The constants are crude (on the frozen trajectory R ranges in [1, 2.08], report B); the content is independence of the drives.
**Prop 2c (eps_a > 0).** `R^eps = R sqrt(vhat)/(sqrt(vhat) + eps_a)`.  If `eps_a <= kappa sqrt(vhat_t)` for all `t <= T`, then `R/(1+kappa) <= R^eps <= R`, so (iii) holds with `(1-b1)/(1+kappa)`.
For nondecreasing gradients `sqrt(vhat_t) >= G_1 = 3(1-delta) g s^2`: sufficient `eps_a <= kappa 3 (1-delta) g_min s^2`.
(Repo default 1e-8 with `eps = 1e-6`: units with `g_0 < 1e-8/(3 eps)` = 0.003 are slightly off; the two weakest units of the design, g_0 = 0.0017, form last: 110 vs 106-108 steps, report E.)
**Converse (Prop 2d).** If `sqrt(vhat_t) <= eps_a/K` for `t <= T` then `|R^eps - mhat_t/eps_a| <= mhat_t/(eps_a K)`: Adam is the EMA-momentum SGD with rate `lr/eps_a` up to the factor `1 + 1/K`.  Since `sqrt(vhat_t) <= max_{k<=t} G_k <= 3(1+delta) g_max U^2` as long as `u <= U`,
`eps_a >= 3 K (1+delta) g_max U^2` is sufficient for GD-like (ordered) behaviour up to `u = U`.  This is far from necessary: the pair experiments (report D, D2) show the ordering restored from `eps_a* ~ c rho s/(1-rho)` (rho = 0.2: eps*/s = 0.3; 0.5: 0.9-1.6; 0.8: 2.7 or more; c ~ 1-1.5), a heuristic (OUTLINE): at the time unit 1 leaves its GD-like phase unit 2 sits at `u_2 = s/(1-rho)`
(GD ratio of the blow-up) with gradient `3 rho u_2^2 << eps_a`, so it keeps moving at the GD-like speed `lr 3 rho u_2^2/eps_a` and blows up after `eps_a (1-rho)/(3 rho lr s)`; it must not do so before unit 1 has finished the traverse of length `U` at speed `R lr`: `eps_a > 3 rho R U s/(1-rho)`.
The design (report F): D3 is exactly the sliding path's set from eps_adam = 1e-2 on (10 units formed at 1e-4), the pair's crossover being at 1e-4 to 1e-3.

## 3. F1c  The pair and duplicates (sign descent; PROVED for the first part)

Pair: atom 1 is the target (`w* = (1, 0)`), atom 2 has overlap rho: `G = [[1, rho], [rho, 1]]`, `c = (1, rho)`; drives `g_1 = 1 - w_1 - rho w_2`, `g_2 = rho - rho w_1 - w_2`.
Both units follow the same path (Theorem 1's induction needs only that the signs of g_1 and g_2 agree) while `g_2 > 0`, with `w_1 = w_2 = omega(t)`; `g_2 = rho - (1+rho) omega`, so the first step at which the drive of atom 2 is
not positive has `omega in (rho/(1+rho), rho/(1+rho) + omega-increment]`: **peak w_2 >= rho/(1+rho)** (up to one step's increment of order `lr * 6 a^2`).  Gradient descent: peak `~ s^3` scale (stays dormant; report C: 1e-6 at s = 0.013).
Numerics (report C): sign descent peaks 0.167 / 0.334 / 0.445 against 0.167 / 0.333 / 0.444; Adam 0.150 / 0.287 / 0.396 (momentum and the lag of vhat make it smaller); GD, shared-v Adam, normalised GD: below 1e-3.
After the peak, atom 2 slides along `w_2 = rho (1 - w_1)` (sliding mode) and the end state is `(1, 0)` up to a chattering of size lr: this noiseless, full-rank example is not a failure of the limit point, only of the PATH (the formed set is wrong during the transit).

Duplicates (`phi_1 = phi_2`): `g_1(t) = g_2(t)` always, so by Lemma 1's argument `u_i(t) = s_i + lr x_t`, `v_i(t) = s_i - lr x_t` with the same walk `x_t`;
`w_i = 6 s_i^2 a + 2a^3` (`a = lr x_t`), `w_2/w_1 = (3 s_2^2 + a^2)/(3 s_1^2 + a^2)`; at the fit (`w_1 + w_2 = 1`, a ~ 0.63) the ratio is `1 - 3(s_1^2 - s_2^2)/a^2 + ...`, i.e. the amplitude is SHARED equally.
Gradient descent (D3, `u_i' = 3 u_i^2 g`): `1/s_i - 1/u_i` is the same for both, so the one with the larger start escapes first and wins (ratio `(s_2/s_1)^3`-small): selection by the mass.  Adam (no closed form): ratio 1.013 (report C).
This is the mechanism behind "Adam forms all units on the design": with a redundant dictionary (21 smooth kernels for 3 features, Gram nearly singular) every direction in the near null space of G is filled: the path of any normalised
method is the l-infinity-steepest path in parameter space, not the gradient flow.

## 4. F1d  What keeps the mechanism (PROVED; numerics report D, D2, F, H, K)

(0) **Characterisation.** If the update of unit i is `p_i += lr * Phi(g_i(1..t) * grad W)` with `Phi` positively homogeneous of degree 0 in the whole gradient history of that unit (sign, Adam with eps_a = 0, RMSprop, Adagrad-type
rules, per-coordinate normalisation in any coordinates -- in particular in the mirror coordinates eta, where the paper's flow is eta_i' = g_i, a coordinatewise normalisation of eta_i' is the same erasure -- and any per-group step size, which only rescales lr),
then units with proportional gradient histories have identical trajectories (induction as in Prop 2a).  A normalisation that keeps the order must multiply ALL units by the same scalar (degree-1 homogeneity in the
shared history), e.g. a global norm or a shared second moment: relative speeds stay `g_i dW/dp`.

(i) Plain GD flow: ordered (closed form).  (ii) Adam with `eps_a` large against the gradients up to the decision: Prop 2d (momentum SGD with rate `lr/eps_a`): the closed-form ordering of (i) holds up to the lag of the EMA (O(1/(1-b1)) steps).
(iii) **Normalised GD with the global norm**: `theta' = -grad F/||grad F||` has the same orbits as the gradient flow (reparametrise by `d tau = dt/||grad F(theta(t))||`, positive and smooth where the gradient does not vanish);
discrete version = gradient descent with the variable step `lr/||g_t||`.  So the ORDER and the set of entries are those of the paper's flow (the sliding ISS), at the price of non-vanishing steps near stationarity (chatter of size lr).
(iv) **A shared second moment** (one scalar `v_t` for all parameters, `adam_g`): `p_{t+1} = p_t + (lr/sqrt(vhat_t)) mhat_t`; for `b1 = 0` it is gradient descent with the variable positive step `lr/sqrt(vhat_t)`, the same orbits in the limit lr -> 0.
Numerically (report C, D, D2, F, H, K), the answer to "does eps_adam large relative to the dormant gradients restore selection?": YES, but `large` has to be measured against `rho s/(1-rho)`, not against the dormant gradient `3 g s^2`
(for s = 1e-3 these differ by 10^3-10^4):
 * pair: eps* / s = 0.3 (rho 0.2), 0.9-1.6 (0.5), 2.7+ (0.8), report D2;
 * design, D3 (21 kernels, eps = 1e-6, 2e5 steps, lr = 0.05 eps_adam x rates): units formed 21 (1e-8), 10 (1e-4), 8 (3e-4), 5 (1e-3), then exactly the path's 3 units from 3e-3 on (exact, set = path, rms 1e-7);
 * design, P0: 3 units, exact, set = path, at eps_adam = 1e-2 and at 1.0 after 3e6 steps (the wall of P0 is at 7000/drive time units; at 2e5 steps nothing has formed yet);
 * additive KAN [5,1] (numpy replication of the spike: GD exact-edge share 0.67 (torch: 0.70), Adam lr 0.01: 0.27 (torch 0.23-0.28), Adam eps_adam = 1 with GD rates 0.60, eps_adam = 0.01: 0.43, shared v: 0.70);
 * shared second moment (adam_g): restores everywhere tested (pair except a borderline rho = 0.8 / s = 0.013 cell, design D3 and P0 exact in 2e5 steps, KAN 0.70); global-norm GD: restores on the pair and on D3, NOT converged on P0 within 2e5 steps (chatter / step budget:
   the clock runs with the norm of the group-rate-weighted gradient, dominated by the centres of the formed units).  Per-group learning rates: they rescale `lr` of a group by a constant; they do NOT help for per-coordinate normalisation because the normalised step is `lr R`, independent of g.
Mirror-type normalisation (steepest descent in the metric of the mirror map): not tested.

## 5. F2a  Products of D small factors; the diagonal two-layer network (PROVED, quadrature)

`w = prod_{i=1}^D z_i`, `F = 1/2 (w - g)^2`, rate lr for each factor, balanced start `z_i = s`: by symmetry (uniqueness of the ODE) `z_i(t) = z(t)`, `z' = lr z^{D-1}(g - z^D)`.
Time from s to `z_theta` with `z_theta^D < g`: `T = (1/lr) int_s^{z_theta} dz/(z^{D-1}(g - z^D))`; `1/(z^{D-1}(g-z^D)) = (1/g)[z^{-(D-1)} + z/(g - z^D)]`, so
`T = (1/(lr g)) [ (s^{-(D-2)} - z_theta^{-(D-2)})/(D-2) + J ]` (D >= 3), `T = (1/(lr g)) [ ln(z_theta/s) + J ]` (D = 2), `J = int_s^{z_theta} z dz/(g - z^D) in [0, z_theta^2/(2(g - z_theta^D))]`.
D = 2 is the diagonal two-layer linear network `y = sum_j u_j v_j x_j` (balanced `u = v`): `w = uv` obeys `w' = 2 lr w (g - w)`: `w(t) = g/(1 + (g/w_0 - 1) e^{-2 lr g t})`, `t_{1/2} = ln(g/s^2 - 1)/(2 lr g)`.
So: D = 2: `T = ln(1/s)/(lr g) + O(1)` = `ln(1/eps)/(2 lr g) + O(1)` with `s^2 = eps b`: a LOGARITHMIC escape ("saddle to saddle"), time of the j-th saddle `ln(1/eps)/(2 lr g_j)`; the gap between two units `(1/(2lr))(1/g_2 - 1/g_1) ln(1/eps)` is of the
order of ln(1/eps) transits (T of one transit from w = g/10 to g/2 is O(1/(lr g))): selection only on the log scale (the paper's depth-2 rule).  D = 3 (one cubic unit, balanced untied variant): `s^{-1}/(lr g) = eps^{-1/2}`; D = 4: `eps^{-1}`; D = 6: `eps^{-2}`.
(The repo's D3 unit has ONE parameter u cubed, not three tied factors: `u' = 3 lr g u^2`, `T = (1/s - 1/u_c)/(3 lr g)`, the same exponent.)  Exponent of eps: `-(D-2)/2` (D >= 3), logarithm for D = 2.
Checked: report G table (quadrature = ODE to 6+ digits, fitted exponents -1.11, -2.02, -3.00, -4.00 for D = 3..6).
**The brief's "(1/(lr signal)) ln(1/eps)" is the D = 2 law; "polynomial 1/eps" corresponds to D = 4 (cubic unit times a linear second layer) in eps, `eps^{-1/2}` for the single cubic unit.**

## 6. F2b  The paper's start in two layers (PROVED)

Two-layer product `f(x) = sum_k b_k sum_l w_kl phi_l(x)` with `w_kl = W(p_kl)` (first layer) and `b_k = W(q_k)` (second layer), both layers started with ZERO COEFFICIENT (`w = 0`: alpha = 0 for the paper's unit, `(s, s)` for D3), loss `F`, `gamma_l = E[r phi_l]`:
`-dF/dp_kl = nabla W(p_kl) b_k gamma_l`, `-dF/dq_k = nabla W(q_k) sum_l w_kl gamma_l`.  At `b = 0, w = 0` both vanish: the whole set `{b = 0} x {w = 0}` consists of equilibria, for every optimizer whose update vanishes with the gradient (GD, momentum, Adam with eps_a > 0).
RBF second layer `sum_m W2[k,m] psi_m(h_k(x))` (the KAN): at `W1 = 0`, `h = 0` and `-dF/dW2[k,m] = psi_m(0) E[r]` which vanishes when the target is centred (or the bias has absorbed the mean), `-dF/dW1[k,l] = E[r sum_m W2[k,m] psi_m'(h) phi_l] = 0` since `W2 = 0`; the (detached) standardisation of h maps the constant to 0.
This explains the KAN development record exactly (`results/ml/kan_dev.md`: pure small start in both layers of [d,4,1]: 5 active parameters, RMSE 0.93-1.05 = the mean predictor, in 12000 steps, for every configuration and for Adam): not a slow escape but a stationary point.  The numpy [1,3,1] KAN reproduces it (report G: gradient at the start 1.5e-19, nothing moves in 3000 steps of GD or Adam).
Perturbed start (one layer `w = a_0 != 0`): for the paper's unit unsaturated (alpha_1(0) = 0, alpha_2(0) = a_0, nu = s, rates 1): `nu^2 - 2 alpha^2` is conserved in each layer; linearisation `alpha_1' = s^4 c* alpha_2`, `alpha_2' = s^4 c* alpha_1`: eigenvalue `lambda = lr s^4 c*` (`s^4 = eps^2`);
symmetric start: `T = int dalpha/(lr (s^2 + 2 alpha^2)^2 alpha (c* - (s^2+2alpha^2)^2 alpha^2))` exactly (quadrature) `= (ln(s/a_0) + O(1))/(lr c* s^4)`; asymmetric `(0, a_0)`: the pair has the exact first integral `alpha_2^2 - alpha_1^2 = a_0^2` (`alpha_1' = Phi alpha_2`, `alpha_2' = Phi alpha_1` with the common factor `Phi = f(alpha_1) f(alpha_2) (c* - w_1 w_2)`, `f(alpha) = s^2 + 2 alpha^2`), so `alpha_{1,2} = a_0 (cosh, sinh)(tau)` on the clock `d tau = Phi dt`, the sum grows like `a_0 e^tau` (the symmetric start `a_0/2` has `alpha = (a_0/2)e^tau`) and the difference `a_0^2/(alpha_1+alpha_2)` dies out: `T(0, a_0) = T_sym(a_0/2) (1 + O(a_0^2))` (numerically equal to 9 digits at a_0 = 1e-6: 497.9021032057 vs 497.9021032290; a first-order argument, the O(a_0^2) not bounded rigorously).  The time is `eps^{-2} ln` (D_total = 6).
Pure small start (any perturbation) therefore escapes only after `eps^{-2}` times: with the KAN's step budget it never does.

## 7. F2c  The hybrid (PROVED)

First layer: one-sided cubic units `w_kl = u_kl^3` from `u = s` (rate 1), second layer `b_k(0) = b_0 > 0` (rate lr_b), nonnegative drives `gamma_l(t) >= 0`.
**Conservation.** `d(u_kl^2)/dt = 6 b_k w_kl gamma_l`, `d(b_k^2)/dt = 2 lr_b b_k sum_l w_kl gamma_l`: `b_k^2 - kappa^2 sum_l (u_kl^2 - s^2) = b_0^2` for all t, `kappa^2 = lr_b/3`.  (Numerics: residual 1e-10, report G.)
**Theorem (escape).** Assume `|gamma_l(t) - gamma_l(0)| <= delta gamma_l(0)` up to the stop time, `g_* = max_l gamma_l(0)`, `n_+` the number of units of the node with positive drive.  Let `T` be the first time `max u_kl = U`.  Then
`int_s^U du/(3(1+delta) g_* u^2 b_max(u)) <= T <= (1/s - 1/U)/(3(1-delta) g_* b_0)`, `b_max(u)^2 = b_0^2 + kappa^2 n_+ (u^2 - s^2)`.
*Proof.* Upper: the unit with the largest drive has `u_*' = 3 b gamma_* u_*^2 >= 3 b_0 (1-delta) g_* u_*^2` (b is nondecreasing since all `w gamma >= 0`), integrate `-d(1/u)`.  Lower: `u_max = max_l u_l` satisfies the differential inequality `u_max' <= 3 b (1+delta) g_* u_max^2` and `b^2 = b_0^2 + kappa^2 sum_l (u_l^2 - s^2) <= b_max(u_max)^2`
because units with `gamma_l >= 0` have `u_l <= u_max` and at most `n_+` terms are positive.  Integrate.  QED.
Consequences.  (i) If `b_0^2 >= kappa^2 n_+ U^2` then `b_max <= sqrt(2) b_0` and `T = (1/s - 1/U)/(3 g_* b_0) (1 + O(delta))` up to the factor `sqrt 2`: the ONE-LAYER law (exponent `s^{-1} = eps^{-1/2}`) with the effective drive `b_0 g_*`; with `b_0^2 >> kappa^2 n_+ U^2` the factor is `1 + O(kappa^2 n_+ U^2/b_0^2)`.  (ii) Different `b_k(0)` rescale the drives of the nodes: the order of entry is by `b_k(0) g`.
(iii) Pure small start with `b_0 = kappa s` (balanced, one unit): `b = kappa u` exactly, `T = (s^{-2} - U^{-2})/(6 kappa g)`: `eps^{-1}`.
Numerics: 12 cells, all inside the bounds (bounds within a factor 1.1 for b_0 = 1, factor ~3 for the pure start because of the sum over the 21 units of a node).

## 8. Gaps (honest list)

1. Adam: the exact statement is for a frozen drive (2a) and the window theorem has crude constants (72.7 for the defaults; real windows 1.03-1.24 in the design).  The `delta`-perturbed Adam trajectories are not proved to coincide, only to stay inside the window.
2. The window theorem is for the one-sided cubic unit (or the u coordinate); the v coordinate of the two-sided D3 and the (alpha, nu) of the paper's unit are covered by numerics only (the same structure; sign-descent equality holds for any W).
3. `eps_a`: the sufficient condition (SGD regime) is proved; the observed threshold `eps_a* ~ c rho s/(1-rho)` (scan over rho = 0.2, 0.5, 0.8, s = 0.0131, 0.00131; grid factor 1.7) is a heuristic (OUTLINE), c not derived.
4. Mirror-type normalisation and per-group learning rates: argued, not tested beyond the shared-second-moment and global-norm variants.
5. F2: the exact stationarity needs `w = 0` in BOTH layers and a centred target; the saturated paper unit and the Gaussian second-layer edges are covered by the algebra (zero gradients) and by the numpy KAN, not by an escape-time law.
   The perturbed pair law is a quadrature for the symmetric start and a numerical equivalence for the asymmetric one.
6. Not done: the real KAN (torch) was not rerun (no torch here); the statements are compared with `results/ml/kan_dev.md` as recorded.
7. The noise, finite-sample and centre-motion effects are outside (units' centres move in the design runs; they are normalised by Adam too).
