# THEORY3B: a deterministic, non-asymptotic noise theorem for the limit object (sliding ISS), and its Gaussian form

Work note of direction B. Code: `experiments/theory3b_noise.py` (constants, path validation, report), `tests/test_theory3b_noise.py`, mutants `tools/mutants/theory3b_noise_own.py`, report `results/synthetic/theory3b.md`.
Status vocabulary: PROVED (written below, each step checked, numerically sanity-tested), COMPUTER-ASSISTED (finite system evaluated in ball arithmetic on floating-point limit data), OUTLINE, FAILED.

## 0. What is obtained (summary)

Setting of Theorem S-K: the **limit object** (sliding inverse scale space, equal walls, duals piecewise linear, blocks relaxed to critical points of the fit, zero slow time), the **continuum Gram**
`<phi_c phi_c'> = Z Gamma(c-c')`, `n = 1`, `K` Gaussian features of the kernels' width at separation `>= Delta`, amplitude ratio `R`, a finite dictionary. Data `y = f* + xi`, `xi` any element of the
inner-product space. Everything is in **normalised units**: lengths in `sigma`, coefficients in `A_min`, `<phi^2>=1` (so `Z=1`), time in `1/(Z A_min)`, walls `1`.

*Gauges of the noise* (the "dual-norm" hypothesis; `N_r` = union over the K features of `[c*_k - r, c*_k + r]`, `D` = dictionary centres):

    e0 = sup_{c in N_r u D} |<xi, phi_c>|/(Z A_min),   e1 = sup_{N_r} sigma |<xi, d_c phi_c>|/(Z A_min),   e2 = sup_{N_r} sigma^2 |<xi, d_c^2 phi_c>|/(Z A_min),
    ez_j = |<xi, zeta_j>|/(Z A_min),   zeta_j = (I - P_T) phi_j   (j in D; T = tangent space of the K-unit model at the truth).

**Theorem N (deterministic).** There is a finite explicit system `S(e0,e1,e2,ez; r)` of inequalities on the noiseless path's data (Section 5), increasing in the gauges (every bound is monotone), such that, if
`r > d_0^k` (the entrant offsets), the noiseless path is "simple" (H1: unique entrants, signs, dormant duals below the walls) and `S` holds, then the limit object on `y` (with the box relaxation of Prop. S-A,
which converges for the noisy fit too, Prop. N1) satisfies

* (i) the K entries are the noiseless ones in the noiseless order, no other kernel enters before the K-th, no exit; the k-th entry time differs from the noiseless one by at most `eta_k` (explicit, `eta_k = O(e)`);
  the k-th block is within `Theta_k = O(e)` of the noiseless block (Euclidean norm in `(v,d)`);
* (ii) after the K-th entry the block is the K-unit least-squares fit and no other kernel enters for a time `>= W_lo` (crude) and `>= W_lo^z` (refined), explicit, `W = O(1/e)`;
* (iii) in that window `|| f_hat - f* - P_T xi || <= Rem(e)`, `Rem = O(e^2)` explicit, with `P_T xi` the first-order oracle of Prop. noise-full (the K-unit fit with known support).

**Theorem N-G (high probability).** Sub-Gaussian noise, variance proxy `s^2`, N samples: with probability `>= 1 - p` the gauges are `<= E(st,p)`, `st = s/(A_min sqrt(N Z))`, with
`E0 ~ st sqrt(2 ln(4 nev/p))` (explicit, Section 7). If `S(E(st,p))` holds, (i)-(iii) hold with probability `>= 1-p` and the window is `W^z >= min_j (1 - |u_j| - ebar_j)/(st ||zeta_j|| L + Rem)`, `L = sqrt(2 ln(4 nev/p))`: `~ 1/s`.

Status: PROVED conditional on the finite system S; S evaluated COMPUTER-ASSISTED (python-flint balls for every inequality; the noiseless limit-object data in double precision, widened by 1e-9 relative).
Not proved: S for all targets (it is a computation on the target, as (a),(b) of Thm S-K), monotone-in-e is by inspection (every bound increases with the gauges; also tested), the limit data enclosure, the sampled window, `n>1`.

---

## 1. Notation and normalisation

Parameters of a block of `k` units: `theta = (v_1..v_k, d_1..d_k)`, unit `i` assigned to feature `l(i)`: coefficient `w_i = A_min v_i`, centre `c_i = c*_{l(i)} + sigma d_i`.
Model `Phi(theta) = sum_i v_i phi(c*_{l(i)} + d_i)` (normalised: `A_min` and `sigma` scaled out, `phi` has norm 1, `<phi_c phi_c'> = Gamma(c - c')`, `Gamma(x) = e^{-x^2/4}`, `psi(x) = (x/2)Gamma(x)`).
`F^{(k)}_y(theta) = 1/2 || y - Phi(theta) ||^2`, gradient `grad_{v_i} F = -<r phi_i>`, `grad_{d_i} F = -v_i <r d phi_i>`, `r = y - Phi`.
Derivatives of the kernel: `||phi||=1`, `||phi'|| = 1/sqrt2`, `||phi''|| = sqrt3/2`, `||phi'''|| = sqrt(15/8)` (norms `m_j = 1, 1/2, 3/4, 15/8, 105/16`).
`a_l = A_l/A_min`, `R = max |a_l|`. The noiseless blocks `theta_k` (critical point of `F^{(k)}_{f*}`), residual `r_k = f* - Phi(theta_k)`, rates of dictionary kernels `p^{(k)}_j = <r_k phi_j>`, nominal duals `u_j(tau)`.
The noisy blocks `hat theta_k`, rates `p^{y,(k)}_j`, duals `u^y_j`, entry times `tau_k^y`. `T_k = tau_{k+1} - tau_k`.

## 2. Proposition N1 (the entry layer with noise; Prop. S-A with the noise in the couplings)

Prop. S-A: for `K` features at distance `>= Delta`, `ell = Delta - 2r >= sqrt6`, `kappa = (K-1) Gamma(ell) < 1`, `vbar = (R(1+kappa) + Omega)/(1 - kappa)`, `E_v = kappa(R + vbar) + Omega`,
`E_c = (K-1)(R+vbar) psi(ell) + 0.43 Omega`, (L) `Gamma(r) > E_v`, (C) `psi(r) > E_c`, (N) `lam_own - H_c > 0`.
**With noise.** `F_y = F_0 - <xi, Phi(theta)> + const`. Hence in the field of unit `i`:
`-grad_{v_i} F_y = a_i Gamma(d_i) - v_i + P_i^y`, `P_i^y = P_i + <xi phi_i>/(Z A_min)`, `|P_i^y| <= E_v^nu := kappa(R+vbar) + e0`;
`-grad_{d_i} F_y = -v_i [ a_i psi(d_i) - C_i^y ]`, `C^y_i = C_i + sigma <xi d_c phi_i>/(Z A_min)`, `|C^y_i| <= E_c^nu := (K-1)(R+vbar)psi(ell) + e1`, for every unit of the box (`|d_i| <= r`, so `c_i in N_r`).
So Prop. S-A's proof for faces and invariance goes through with `vbar = (R(1+kappa) + e0)/(1-kappa)`, `E_v^nu`, `E_c^nu` ((L): `Gamma(r) > E_v^nu`; (C): `psi(r) > E_c^nu`). The inner problem of the noisy entry layer is the gradient flow of `F_y` with positive mobilities (the network's `eps`-flow too, since the profile is increasing), so Nagumo's invariance gives (i) of Prop. S-A, LaSalle convergence to a critical point of `F_y` in `B`, and the entrant increases monotonically.
**Hessian.** The Hessian of `-<xi,Phi>` is block diagonal over units, own block `[[0, -g1],[-g1, -v g2]]` with `|g1| <= e1`, `|g2| <= e2` (`g1 = sigma<xi d phi_i>/(ZA)`, `g2 = sigma^2<xi d^2 phi_i>/(ZA)`). The own block of `F_0` is `[[1, s],[s, t]]`, `s = a psi(d)`, `t = a v psi'(d)`, `psi'(d) = (1/2 - d^2/4)Gamma(d)`; on `D* = {|d_i| <= d*, v_i in [Gamma(d*) - E_v^nu, R + E_v^nu]}`,
`d* = 2E_c^nu/Gamma(r) < r`: `|s| <= R psi(d*) + e1`, `t >= (Gamma(d*) - E_v^nu)(1/2 - d*^2/4)Gamma(d*) - vbar e2`; the smallest eigenvalue of `[[1,s],[s,t]]` increases with `t` and decreases with `|s|`, so it is `>= lam_own`.
The cross part (other units, unfitted features) is bounded as in S-A by `H_c = max(rowv, rowd)` (Omega = 0). Hence `lam := lam_own - H_c` is a lower bound of the Hessian of `F_y` at every point of the convex set `D*`. **(N)**: `lam > 0`.
Consequences: in `B` every critical point lies in `D*`, there is at most one (strictly convex), the flow converges to it. PROVED (the proof of S-A with the modified constants; the monotone-in-`e` is clear: `vbar, E_v, E_c, d*` increase, `lam` decreases).
(`e0` enters `vbar` and `E_v`, `e1` enters `E_c` and `s`, `e2` enters `t`.)

## 3. Lemma N2 (the noisy block is within `Theta_k` of the noiseless one)

Let `S`'s box conditions hold (so the noiseless block `theta_k in D*_0 subset D*_nu`: also checked on the numbers). Let `hat theta_k` be the critical point of `F^{(k)}_y` in `D*`. Then
`||hat theta_k - theta_k||_2 <= || J(theta_k)^T xi ||/lam <= Theta_k := sqrt( sum_{i<=k} (e0^2 + v_i^2 e1^2) )/lam`.
*Proof.* `grad F_y(theta_k) = grad F_0(theta_k) - J(theta_k)^T xi = -J^T xi` (`theta_k` is a critical point of `F_0`), components `-<xi phi_i>`, `-v_i sigma <xi d phi_i>` bounded by `e0`, `|v_i| e1` (the centre `c_i^k` lies in `N_r`).
`grad F_y(hat theta_k) = 0`. Mean value theorem on the convex `D*`: `< grad F_y(hat theta) - grad F_y(theta_k), hat theta - theta_k > >= lam ||hat theta - theta_k||^2`; Cauchy--Schwarz. QED

*Rate error.* For a dictionary kernel `j`: `p^{y,(k)}_j - p^{(k)}_j = <xi phi_j> - <Phi(hat theta_k) - Phi(theta_k), phi_j>`, and
`<Phi(hat theta) - Phi(theta), phi_j> = sum_i [ Delta v_i Gamma(hat c_i - c_j) + v_i (Gamma(hat c_i - c_j) - Gamma(c_i - c_j)) ]` so
`| p^y_j - p_j | <= eps_j^{(k)} := e0 + Theta_k sqrt( sum_i ( Gbar_ij^2 + v_i^2 psibar_ij^2 ) )` (Cauchy--Schwarz over `(Delta v_i, Delta d_i)`),
with `Gbar_ij = max_{|x| in [max(0,|c_i - c_j| - Theta_k), |c_i - c_j| + Theta_k]} Gamma(x)` and `psibar_ij` the maximum of `psi` on the same interval (`sqrt2` included if inside; `|Gamma'| = psi`). For `k=0`: `eps = e0`.

## 4. Lemma N3 (dual lines) and Theorem N (i)

Between entries the residual is constant, so each inactive dual is **linear**: stage `k` (after `k` entries), noisy `L^y_j(tau) = u^y_j(tau_k^y) + p^y_j(tau - tau_k^y)`, nominal `L^0_j(tau) = u^0_j(tau_k) + p^0_j (tau - tau_k)`.
Induction hypothesis `H_k`: the noisy path has had the entries `j_1..j_k` of the noiseless path in order, its blocks are the `hat theta_i` (Prop. N1: the entry layer relaxes to the unique critical point in the box), `|tau_i^y - tau_i| <= eta_i` (`i <= k`), `|u^y_j(tau_k^y) - u^0_j(tau_k)| <= ebar_{k,j}` for inactive `j`. (`eta_0 = 0`, `ebar_{0,j} = 0`.)
**Lines.** For all `tau`: `|L^y_j - L^0_j|(tau) <= ebar_{k,j} + eps_j |tau - tau_k| + (|p^0_j| + eps_j) eta_k` (expand `L^y - L^0 = Delta u + (p^y - p^0)(tau - tau_k) - p^y (tau_k^y - tau_k)`).
**Entrant** `j*`, nominal approach speed `v = |p^0_{j*}| > eps_{j*}` (so the noisy dual moves toward the same wall, simply and transversally): `L^0(tau_{k+1}) = L^y(tau_{k+1}^y) = sgn`, hence `|L^y - L^0|(tau^y_{k+1}) = v |tau^y_{k+1} - tau_{k+1}|`, and with `|tau^y_{k+1} - tau_k| <= T_k + eta_{k+1}`:
`eta_{k+1} <= [ ebar_{k,j*} + eps_{j*} T_k + (v + eps_{j*}) eta_k ] / (v - eps_{j*})`.
**Competitors** `j != j*` (inactive): on `tau in [tau_k^y, tau^y_{k+1}] subset [tau_k - eta_k, tau_{k+1} + eta_{k+1}]`, `|L^0_j|` is convex piecewise linear so `<= max( |u^0_j(tau_k)| + |p^0_j| eta_k, |u^0_j(tau_{k+1})| + |p^0_j| eta_{k+1} )`; hence
`sup |L^y_j| <= M_{k,j} := max(...) + ebar_{k,j} + eps_j max(eta_k, T_k + eta_{k+1}) + (|p^0_j| + eps_j) eta_k`. If `M_{k,j} < 1` for all `j`, none enters before `j*`. (Order.)
**State error update.** `ebar_{k+1,j} <= ebar_{k,j} + eps_j T_k + (|p^0_j| + eps_j)(eta_k + eta_{k+1})` (compare at the noisy entry time, then shift to the nominal one).
This closes `H_{k+1}` (blocks: Prop. N1 at the new stage; the new block `hat theta_{k+1}` is within `Theta_{k+1}`, Lemma N2). Induction up to `k=K`. The coefficients of the box stay `>= v_lo > 0` (no exit); signs are those of the dual (`sgn p^y = sgn p^0`). PROVED. QED (i)

## 5. The window (ii)

After the K-th entry the noiseless residual is zero and the duals stop; `p^0_j = 0` (we keep the computed nominal value `|p^{(K)}_j| ~ 1e-16`). For `tau >= tau_K^y` and inactive `j`:
`|u^y_j(tau)| <= |u^0_j(tau_K)| + ebar_{K,j} + (eps^{(K)}_j + |p^{(K)}_j|)(tau - tau_K^y)`, so no kernel enters before `tau_K^y + W_lo`,
`W_lo = min_j (1 - |u^0_j(tau_K)| - ebar_{K,j})/(eps_j^{(K)} + |p_j^{(K)}|)` (crude). In particular `W_lo = O(mu/e)`.
**Refined window.** By Section 6, `Phi(hat theta_K) - f* = P_T xi + rho`, `||rho|| <= Rem`. Then `p^y_j = <xi, (I - P_T)phi_j> - <rho, phi_j>/(...)`: `|p^y_j| <= ez_j + Rem` (normalised, `||phi_j|| = 1`). Hence
`W_lo^z = min_j (1 - |u^0_j(tau_K)| - ebar_{K,j})/(ez_j + Rem + |p^{(K)}_j|)` is also a valid lower bound; for small `e` it is asymptotically sharp (numerics: ratio `0.995`).
(The cancellation in `(I-P_T)phi_j` — a kernel near a fitted feature is almost in `T` — is lost by the triangle inequality in the crude bound; the gauge on `zeta_j` keeps it.) PROVED (given Section 6).

## 6. Lemma N5: the oracle remainder (iii)

At the truth `theta* = (a, 0)`, `J* = J(theta*)`, `H* = J*^T J*` (Gauss-Newton Gram, equal to the Hessian of `F_0` there, the residual vanishes), `lmin, lmax` its extreme eigenvalues (explicit, `<phi_k phi_l> = Gamma`, `<phi_k, a_l phi'_l> = a_l (x/2)Gamma`, `<a_k phi'_k, a_l phi'_l> = a_k a_l (1/2 - x^2/4)Gamma`, `x = c_k - c_l`).
`h = hat theta_K - theta*`, `||h|| <= Theta_K` (Lemma N2 with `theta_K = theta*`). Taylor: `Phi(hat theta) - f* = J* h + q`, `||q|| <= (1/2) M_2 ||h||^2`, `M_2 = 1/sqrt2 + (sqrt3/2) vbar` (since `D^2Phi[h,h] = sum_k 2 h_{v_k} h_{d_k} phi' + v_k h_{d_k}^2 phi''`, norms `1/sqrt2`, `sqrt3/2`, `sum |h_v||h_d| <= ||h||^2/2`).
`E := J(hat theta) - J*`: `||E|| <= L_J ||h||`, `L_J = sqrt(1/2 + (1/sqrt2 + (sqrt3/2) R)^2)` (`||phi(c) - phi(c*)|| <= |d|/sqrt2`, `||v phi'(c) - a phi'(c*)|| <= |Delta v|/sqrt2 + a |d| sqrt3/2`), and `||E^T xi|| <= (sqrt2 e1 + R e2) ||h||` (the gauges over the box applied to the same differences: `|<xi, phi(c) - phi(c*)>| <= e1 |d|`, `|<xi, v phi'(c) - a phi'(c*)>| <= |Delta v| e1 + a|d| e2`).
Stationarity at `hat theta`: `J(hat theta)^T (f* + xi - Phi(hat theta)) = 0`, i.e. `(J* + E)^T (xi - J* h - q) = 0`, so `H* h = J*^T xi - J*^T q + E^T (xi - J* h - q)`.
`Phi(hat theta) - f* - P_T xi = J*(h - H*^{-1} J*^T xi) + q = (I - P_T) q + J* H*^{-1} [ E^T xi - E^T (J* h + q) ]`, `P_T = J* H*^{-1} J*^T`, `||J* H*^{-1}|| = lmin^{-1/2}`, `||J*|| = lmax^{1/2}`:

    Rem = (1/2) M_2 Th^2 + lmin^{-1/2} [ (sqrt2 e1 + R e2) Th + L_J Th ( sqrt(lmax) Th + (1/2) M_2 Th^2 ) ],   Th = Theta_K.   PROVED.

`||P_T xi||^2 = (J*^T xi)^T H*^{-1} (J*^T xi)`, so the oracle's squared error is the quadratic form of the `2K` gauges (mean square `2K s^2/N` for white noise): `Rem` is `O(e^2)`, while `||P_T xi|| = O(e)`.

## 7. Lemma G: sub-Gaussian noise (high-probability version)

`xi_s` independent, `E exp(t xi_s) <= exp(t^2 s^2/2)`, `<xi, g> = N^{-1} sum xi_s g(x_s)` sub-Gaussian with proxy `s^2 <g^2>/N`. Normalised gauges `g_j(c) = sigma^j <xi, d_c^j phi_c>/(Z A_min)` have proxy `st^2 m_j`, `st = s/(A_min sqrt(N Z))` (`m_0, m_1, m_2 = 1, 1/2, 3/4`), `<xi, zeta_j>/(ZA_min)` proxy `st^2 ||zeta_j||^2/Z`.
Net: `M = max(10, ceil(4 r sqrt N))`, points `c*_k + i r/M`, `|i| <= M`; `nev = 3K(2M+1) + 2d`; `L = sqrt(2 ln(4 nev/p))`; Gaussian tail: all net values obey `|g_j| <= t_j = st sqrt(m_j) L`, `|<xi,zeta_j>|/(ZA) <= st ||zeta_j||_N L` with probability `>= 1 - p/2` (each event fails with prob `<= 2 exp(-t^2/(2 proxy)) = p/(2 nev)`).
Continuum: for `c` in the box and the net point `c_n`, `|c - c_n| <= m/2`, `m = r/M`; `|g_j(c)| <= |g_j(c_n)| + (m/2) sup |g_{j+1}|` since `d g_j/dc = g_{j+1}/sigma`; the chain is closed by Cauchy--Schwarz,
`sup |g_3| <= ||xi||_N ||sigma^3 d^3 phi||_N/(Z A_min) = st sqrt(N m_3) (||xi||_N/s)`, with `||xi||_N^2 = ||xi||^2/N <= s^2 (1 + 2 sqrt(tau/N) + 2 tau/N)`, `tau = ln(2/p)`, with probability `>= 1 - p/2` (Hsu--Kakade--Zhang tail for sub-Gaussian vectors).
So with probability `>= 1-p`: `E2 = t_2 + (m/2) M3`, `E1 = t_1 + (m/2) E2`, `E0 = t_0 + (m/2) E1` bound `e2, e1, e0` (the dictionary part of `e0` is `<= t_0`). PROVED.
**Theorem N-G.** If `S(E0,E1,E2,Ez; r)` holds with these bounds, then with probability `>= 1-p` (i)-(iii) hold with the bounds computed at `E`. Since `S` is monotone, `st_crit(p)` (largest `st` with `S(E(st,p))` true) is well defined by bisection.
Window: `W_lo^z >= min_j (1 - |u_j^0| - ebar_j)/(st ||zeta_j|| L + Rem)` (+ nominal rate): `~ 1/(s sqrt(ln(nev/p)))`.

## 8. Gaps (honest list)

1. `S` is a finite system on the target, not a closed-form condition in the margins; closed-form corollary only through the recursion (the linear constants are displayed in the report).
2. Limit-object data (`theta_k`, rates, hitting times) are double precision (widened by `1e-9` relative in the ball run; the Newton residual `< 1e-13`); the enclosure of the limit object is not rigorous (as T1b).
3. The suprema over `N_r` of the gauges are over a box of radius `r` about the *true* features; the noisy centres are inside it by Prop. N1 (box invariance).
4. The statement is for the **limit object**, not for the trained network at finite `eps` (Theorem T1 with noise is not proved here). For the network, the noise adds to the dual deviation of Thm S-T1 an additional `O(e)` of the same row-sum structure: OUTLINE only.
5. Continuum Gram only: in the sampled window the Gram departs from `Z Gamma` within two widths of the edge; the validation window is wide so that this is `<1e-12`. `n = 1` only (Prop. N1 is as S-A: the planar case is by hand).
6. The remainder `Rem` is rigorous but 3 orders conservative in the validation (triangle inequalities in `E^T xi`, `||D^2Phi||`); the crude window `W_lo` loses the cancellation `(I-P_T)phi_j` (the refined one does not, up to `ebar`).
7. The `ebar` (dual state errors) are accumulated with the crude `eps` at every earlier stage: a stage-wise orthogonalised version (full Hessian `H_k` incl. residual curvature, `zeta^{(k)}_j = phi_j - J_k H_k^{-1} J_k^T phi_j`) would tighten `eta`, `ebar` at larger `e`; not done.
8. Not covered: unequal walls, correlated or heavy-tailed noise (only sub-Gaussian proxies), data-dependent `r`.

## 9. Remarks and numbers

**Leading order for K = 1 (closed form).** Stage 0 has no block, `eps_j = e0`, `eta_1 = e0 T/(v - e0)`, and the competitor condition reads `mu_j > p_j eta_1 + e0 (T + eta_1)`, i.e. to leading order
`e0 T_1 (1 + p_j/v) < mu_j`; for the runner-up of the same feature `p_j/v ~ Gamma(h) ~ 0.96`, so `e0 <= mu/(2 T_1) (1 + O(e))`. The report finds `c = e0_cert tau_K/mu = 0.50` on every K = 1 target and on every R = 1, K = 2 target
(where the cross-feature terms are not limiting), and `0.10 - 0.38` when the amplitude ratio is 2.4 or K = 3. This is the "explicit multiple of the margin": the entry order survives a noise gauge `e0 < c mu/tau_K`, `c <= 1/2`.
The window after the K-th entry is, to leading order, `W = min_j (1 - |u_j| - ebar_j)/(ez_j + Rem)`, `ebar_j = e0 T (1 + p_j/v)`-type: `W ~ mu_j/(ez_j)` for `e -> 0` (the numbers: the ratio `W_lo^z/actual = 0.995` on one target at `0.01 t_cert`; median over all rows 0.24 at `0.1 t_cert`; the Gaussian-level limit `W st -> c_W` is reached to 4 digits below `1e-4 st_crit`).

**Is the certificate sharp?** The adversarial family (the runner-up kernel pushed with the noise `+- phi_j`) attains `eta_k` to 0.998 and `W_lo^z` to 0.998 of the actual values, and the order is lost at 3 times the certified scale: the first-order part of the constants cannot be improved much for the
adversarial worst case; against white Gaussian noise it is conservative by a median factor of 20 in the scale (the adversarial gauge `e0` is the supremum over 33 + 55 K kernels, white noise realises a small part). The O(e^2) remainder is conservative by `~ 10^3` (triangle inequalities) while the first-order oracle term is exact.

**Why a gauge on `zeta_j = (I - P_T) phi_j`.** After the K-th entry the nominal residual vanishes, so the rate of a dormant kernel is `<xi - (Phi(hat theta) - f*), phi_j> = <(I - P_T) xi, phi_j> - <rho, phi_j>`: the block *absorbs* the noise component in the tangent space, and a kernel close to a fitted feature (`phi_j` almost in `T`) gains almost no rate.
With the plain `e0` gauge (triangle inequality) the window is 3 orders too short; with `ez_j` it is sharp (the first-order proposition's `v_j = <R_K xi, phi_j>`, of the paper). For the earlier stages `k < K` the same orthogonalisation needs the full Hessian `H_k` (with the residual curvature) and a second-order estimate: not done (gap 7).

**What was tried and failed / not attempted.** (a) A sup-norm (Gershgorin) dual-propagation as in S-T1 is *not needed* for the limit object: the duals are exactly linear between entries, so the propagation is the scalar recursion of Section 4 (the row-sum margins of S-T1 are for the network). (b) Using only the dictionary gauge `e0` for all rates (the first version): valid, 500x short on the window; replaced by `ez`.
(c) A stage-wise orthogonalised version, a Lipschitz-in-time bound for the network at finite `eps` (Theorem T1 with noise): not carried out.

**Numbers (results/synthetic/theory3b.md).** 24 noiseless limit paths (K = 1,2,3; separations 6, 8 widths; R = 1, 2.4), 3 noise families x 14 multiples: 288 rows with S true, 0 violations; Gaussian: 7 targets x 4 levels x 60 draws: 0 violations, all gauges below the bound; 600 random adversarial draws: 0 violations; 32 entry-layer flows with noise stay in the box and end at the noisy critical point within `Theta_2`; ball arithmetic certifies S on 72 of 72 at 0.9 `t_cert`.
