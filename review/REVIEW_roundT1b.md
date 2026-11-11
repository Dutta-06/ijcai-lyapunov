# REVIEW_roundT1b: referee report on the widening of Theorem T1 (THEORY_T1.md section 12)

Archived verbatim from the referee's hand-back (2026-11-09). The referee modified nothing in the repository; its own code lives outside the repository.

Scope. I read section 12 of THEORY_T1.md, REVIEW_roundT1.md, and `theory_t1b_wide.py`, `theory_t1b_planar.py`, `theory_t1b_rigorous.py`, `tests/theory_t1b_test.py`, `theory_t1b.md`, `theory_t1b_raw.npz`, and the parts of `exp_theory_t1.py` that matter (`cert`, `box_t1`, `delta_t1`, `follows`). I modified nothing in the repository. My code is in `referee_T1b_work/` (scratch):
- `dT.py`, `d059.py` (box and regime in balls)
- `myrig.py` (my own flint transcription of S(eps)), `cmp1.py`, `cmp2.py`, `chain.py`, `stop9.py`, `cover2.py`, `mut.py`, `inc9.py`
- `yhat.py`
- `lem124.py`
- `plan.py`, `ext.py` (my own planar integrator and checks), with logs `plan.log`, `plan2.log`, `ext.log`

**Bottom line.** I found no mathematical error in the widening. Claims (1), (2) and (4) hold as mathematics. Claim (3) holds as a method, with caveats (finding 2). The problems are wording and scope. The suggested status text in 12.5 overstates the evidence in several places, listed below. There is no blocker.

## Findings

**1. [major, wording] "0 violations in 345 / 81 / 771 pairs" counts vacuous pairs. The violation test is also weaker than the theorem.**
- Location: 12.1, 12.2, 12.5 (the "771 target-eps pairs" sentence), `theory_t1b.md` parts K and P.
- Evidence from `theory_t1b_raw.npz` and `theory_t1_raw.npz`: the number of pairs with `cert` true is
  - K=4..6: 204 of 345
  - planar: 41 of 81
  - K<=3 (the earlier table, 345 pairs): 228 of 345
  - repository network (Part E): 123 of 184

  So 473 of the 771 pairs are informative. A pair with `cert` false cannot be a violation.
- Every informative pair lies at eps <= eps_0. By the paper's own table that is at least 1.5 decades below the measured threshold (ratios 32 to 10^6; 100 to 1000 in the plane). The check therefore has no power near the boundary.
- `follows` and `follows2` test only the formed set at T_fin (`formed == ent`). They do not test item 1 of the theorem, "no dormant unit ever above theta". The ingredient diagnostics cover that, on 30 and 18 trajectories only.
- Cheapest fix: write "no counterexample among 473 pairs with S true (all at least 1.5 decades below the measured threshold; final formed set checked; ingredient inequalities on 48 trajectories)". Drop "771".

**2. [major, wording] "Interval arithmetic for every real eps" is not a proof of S(eps) as stated. Several hypotheses are never inside the verifier.**
- Method: sound as far as I can tell (see Verified). The residual non-rigorous inputs are:
  - the limit-object data: critical points, entry times, duals at breakpoints, rates p_e, margins, and the entrant/feature assignment, all taken from `hybr`/float code as exact numbers (12.3 item 1 says so);
  - `float(ball.upper())` comparisons in all final tests;
  - a `1e-6` flatness threshold in `inv_of`, which is sound.
- Hypotheses that `verify` does not check at all and takes from a prior float `cert` run: r > d_0, the hypothesis that the block critical point lies in D-hat (not (a)/(b)), (N), the sign test, and theta < v_#.
  - I checked that the critical point is in D-hat for all 10 Q targets (`yhat.py`): block offset <= 3.0e-4 against d-hat >= 0.042, coefficients v/a_min within 1e-6 of a_i/a_min. It holds, but nothing in the pipeline asserts it.
  - theta < v_# (part of S5) is checked by neither `cert` nor `verify`. It is trivially true: v_# is about 0.45 a_min and theta <= 0.1 a_min.
- The bootstrap in section 6 needs a strict inequality, E_j(D) < D_j/1.0001. `verify` tests `Etot <= Dv`. In the chain the worst ratio is about 0.94, so this is formal only.
- Cheapest fix:
  - say "S(eps) evaluated in interval arithmetic with floating-point limit-object data (critical points, entry times, assignment) and with the hypotheses r > d_0, block critical point in D-hat, (a), (b) checked in floating point";
  - add assertions for r > d_0, critical point in D-hat and theta < v_# to `verify`;
  - use `E < D/1.0001`.

**3. [major, scope] Planar claim. The hand proof is fine, but the planar object is narrower than the wording suggests.**
- The dictionary is the lattice restricted to 3.5 sigma around the features (`make_target2`, rad=3.5), with 55 to 168 kernels. The statement 12.5 says "a finite dictionary with lattice offsets d_0", which hides this.
- `cert2d` is the source of `cert` rewritten by `inspect`/`replace`, in doubles only. `rgs2` clips its grid to the disc and then uses the full half-diagonal step*sqrt2/2 as the Lipschitz margin. Near the boundary the nearest retained node can be farther than that, so the margin is slightly optimistic. Its effect on D is at most a few percent.
- Truncation test (mine, `ext.py`): K=1, rad 3.5 (J=57) against rad 7 (J=223, 4 times the kernels). Both use the same seed, but the lattice is anchored on the hull, so the phase differs; d_0 is 0.12 and 0.33. `cert2d` stays true at 10^-6.5 and 10^-7.5, and D grows from 76 to 87. The network follows, and the largest dormant coefficient is 1.5e-7, well below theta = 0.007. This is evidence that the truncation is harmless. It is not verified for a full planar window. The K=2 run with rad 6.5 (J=330) did not finish.
- 12.5 quotes Delta_T1 "about 4.9 to 7.5 sigma for K <= 3" for both n=1 and n=2. In the plane d_0 reaches h/sqrt2 = 0.59, and the same box gives 5.2 to 7.9 sigma. My values:

  | R | K=2 | K=3 |
  |---|---|---|
  | 1 | 5.19 | 5.59 |
  | 2.4 | 6.06 | 6.38 |
  | 10 | 7.58 | 7.81 |

  These agree with the author's table P.
- Cheapest fix: add "n=2: square-lattice kernels within 3.5 sigma of a feature; double precision; Delta_T1 5.2 to 7.9 sigma at d_0 = 0.59 sigma".

**4. [minor, wording] The "8 of 10" interval-cover result is a computational-budget outcome, not a failure of S.**
- Evidence for target (2,8,2.4,14), which stopped at 10^-7.09 with `wmin = 1/512` (`stop9.py`): at the stopping point, balls of width 1/2048 or narrower pass with all four candidate scalings. The fragile condition is dc: dcmax goes 0.005, 0.033, 0.43, 12.5 as the ball widens from 1e-5 to 2e-3 decades. This is the window-start gating `wst[i] > q0` becoming undecided.
- With `wmin = 2^-14` (`cover2.py`) the chain continues: 298 balls to 10^-8.0 in 900 s, with minimum width 1/1024 decade.
- The other incomplete target (K=3, R=1.5) is a time-budget stop, as the report says.
- Cheapest fix: lower `wmin` and raise the budget, or state that the two failures are artefacts of the ball width.

**5. [minor] Soundness formalities in `verify` (none changes a verdict at the observed margins).**
- Guards like `ratio.log() if lo(ratio) > 1 else arb(0)` (T2, T3lo, T3hi) and `Y0 if lo(Y0) > 0 else 0` replace an undecided ball by 0 instead of an upper bound. The effect is bounded by log(1+rad)/kd.
- Float edges in `np.linspace(c_lo, c_hi)` for `rig_rgs` and `x_of_s` are not outward-rounded. The `rig_rgs` edges are one ulp.
- The `arb` bug found by the author (the square of a ball containing 0) is fixed by `sq`. I tested flint 0.9.0: `(0 ∪ 1e-30)**2` has a lower end of -5e-61, so the bug is real. The remaining `x*x` uses (`Gam_r`, `psi_r` on straddling subinterval balls) give sound but loose enclosures, with upper end exp(a^2/4) > 1, so I found no similar bug.
- Fix: use `ub/lb` comparisons and `.upper()` throughout. Cosmetic.

**6. [minor] The claimed verify-versus-cert agreement is tested only for K=1.**
- `tests/theory_t1b_test.py::test_ball_arithmetic_verification_agrees_with_the_float_certificate` uses the K=1 target (1,0,1,1). The "about 0.1% on 5 targets" statement for K=2,3 in 12.3 has no test.
- My independent flint transcription (`myrig.py`) is written from sections 2 to 6 of the note. I had read `verify` first, so it is not blind. Results:
  - Point eps = 10^-7.5, target (2,6,1,11): Etot of the largest kernel is 109.41 (author), 109.52 (mine, 1500-cell pre-drift mesh) and 108.80 (mine, 6000 cells). T1 is 44.679 against 44.667, T3 is identical to 12 digits, and Omega and dcmax agree to 0.03%.
  - Over the author's actual chain for (2,6,1,11) (31 balls, down to 10^-12) and (3,8,1,31) (17 balls, to 10^-24), I evaluated 124 and 51 sample points inside the balls. In all of them my point value was at most 0.13% above the author's ball bound, and with the finer pre-drift mesh it was below it.
  - Mutation test (`mut.py`): `verify` with D taken from 10^-7.5 rejects eps = 10^-7.0 and above, through the dc and fixed-point conditions.
- Fix: add a K=2 and a K=3 case to the test.

**7. [minor] Lemma 12.4 is correct and unused by the verified results. The "favourable sign" claim in 12.5 is too strong.**
- The proof text is loose. The correct chain (Lemma 2.2, first form) is:
  - B <= (a/m) s_eps (eps + 2vs s_eps^2/(1 - s_eps^2));
  - with s_eps = cot(theta)/kappa, this is at most (1 - s_eps^2)^-1 (2 vs a/(m kappa^3)) cos(theta)/sin^3(theta);
  - using lambda = pi s/(2 theta) and lambda = pi a kappa/(4 G vs), this equals 4/(G gamma^2 m) s^-3 (1 - s_eps^2)^-1 * theta^3 cos(theta)/sin^3(theta);
  - then theta^3 cos(theta)/sin^3(theta) < 1 by Adamovic-Mitrinovic, and s_eps <= sigma_0/s because cot(theta) <= 1/theta.
- The applicability condition lambda >= 25 means eps <~ 0.03 (lambda(0.03) = 25.7), and the lemma needs eps < 2 vs. It should be stated; "(0, eps_top]" is undefined.
- Numerics: 200-bit arithmetic over about 4000 points per eps, eps from 1e-16 to 0.03, s in [25, 1e6] (and up to lambda). The maximum of B s^3 (1 - 400/s^2)/33.33 is below 1 everywhere (0.99999999878 at 1e-16, 0.395 at 1e-2). The "+6e-8 overshoot" in Part U is float roundoff.
- The lemma bounds an envelope. It does not imply that the exact mass and drift are monotone in eps, and the layer constants are still unsigned. 12.4 says so correctly. The summary table row "mass/drift/(S1)/(S5) have the favourable sign" should read "envelope version".

**8. [minor, wording] Statements about K.**
- 12.1: "No other constant depends on K". In the code, h_dev has a factor (K-1), c_F enters multiplied by K in V_tail, E_c and the row sums carry (K-1), and the old-unit sum has k-1 terms. The proof is K-general; the sentence is not literally true.
- "For any K <= 64 the theorem holds as stated" says only that the box regime (R) exists. S(eps) has been evaluated only for K <= 6. The author's own estimate is +75 in D per old unit, so D is about 5000 for K=64. Say "regime (R) exists for every K <= 64; S not evaluated for K > 6".
- The same-size grid-quantisation caveats for eps_max/eps_0 apply to K=4..6: 3 of 15 targets have near-ties (gap < 0.01), giving ratios 10^3 to 10^6.

**9. [minor, wording] 12.2(b).**
- It says the dormant contributions in row_d are "replaced by G2_sup2". The box constants 0.43 and 1/2 are already sups of the planar operator norms (sup of the norm of the gradient of Gamma is 0.4289; of its Hessian, 1/2 at x=0). `G2_sup2` enters only the dual-deviation part (drift, pre-drift, weights), where it is correct: max(|Gamma''|, Gamma/2).

**10. [minor] The "13 to 56 times for K <= 3 (referee)" in 12.5 comes from three targets in the previous report.**
- It is not a K <= 3 statement. Say "on three targets".

## Verified (re-derived, and how)

- **Delta_T1(K), claim 1.**
  - Own box code (`mybox.py`, arb 80 bits, same conditions as `box_t1`, dormant budget Omega = 0.02, varsigma = 0.01, d_0 = 0.4167). On the author's grid (0.1 in Delta, 0.05 in r) I reproduce every entry I computed (K = 4, 8, 16, 64; R = 1, 2.4, 10; also K = 2, 3):
    - R=1: 5.5, 6.0, 6.3, 6.9
    - R=2.4: 6.3, 6.6, 6.9, 7.4
    - R=10: 7.6, 7.9, 8.1, 8.5
  - A finer r-grid (step 0.005) gives 0.07 to 0.15 lower: K=4: 5.44, 6.19, 7.50; K=8: 5.85, 6.54, 7.77; K=16: 6.18, 6.83, 8.00; K=64: 6.75, 7.34, 8.42.
  - Growth law: Delta^2 = 23.6 + 5.4 ln(K-1) (R=1), 32.6 + 5.2 ln(K-1) (R=2.4), 50.8 + 4.9 ln(K-1) (R=10), over K = 2, 3, 4, 8, 16, 64. Residual in Delta <= 0.03 sigma. So sqrt(log K) is confirmed.
- **K >= 4 proof.**
  - The induction over layers, the triangular recursion, and Lemma 4.1 are K-free except through kappa_c = (K-1) Gamma(ell), the (K-1) factors in the row sums, E_c and h_dev, and the k-1 old-unit summands.
  - Hypothesis check on the 10 Q targets: block critical point in D-hat (see finding 2).
  - K=4..6 certificate against the network: the "0 violations" is correct and the informative count is 204 (finding 1). I did not rerun those integrations.
- **n=2 hand proof, each step (claim 2):**
  - Radial box: rho' = -mu_c v [a psi(rho) - eps <d/rho, C>], with |C|_2 <= sum (a_j+v_j) psi(ell) + 0.43 Omega, since |grad Gamma| = psi. E_v and E_c are unchanged.
  - Hessian of Gamma: Gamma (xx^T/4 - I/2), eigenvalues (|x|^2/4 - 1/2) Gamma and -Gamma/2. The operator norm is |Gamma''|(|x|) for |x| >= 2, and ell >= sqrt6.
  - The cross-block row sums were rederived entry by entry from <e_i, e_j> and agree with row_v, row_d.
  - Own block: a v (I/2 - dd^T/4) Gamma, with the radial eigenvalue smaller than the tangential one. The coupling a grad Gamma has norm <= s_max.
  - Lemma 4.1 with Euclidean blocks and a scalar centre mobility: valid by the argument of 12.2(d). The Dini derivative of |z_d|_2 is well defined where the maximum is positive. The mean-value Hessian lies on the convex D-hat (a product of discs and intervals). V_diam = max(2 d-hat, width) is the Euclidean diameter.
  - The scalar mobility is what makes the argument independent of the direction. With a matrix mobility it would fail, as the note says.
  - Dormant count in the plane: only sum |w_l| enters, and it is computed over the actual list (S1). Kernels beyond 3.5 sigma contribute about 33/lambda^3 each. For a hypothetical 1.4e4 extra kernels at eps = 1e-4 this is about 1.4e-3, far below Omega a_min = 0.014. This is why the truncation seems harmless, but it is not verified.
- **cert2d.** Its construction by `inspect`/`replace` covers all geometry-dependent lines of `cert`:
  - X for dormant kernels
  - X for entrants
  - `resid_grad_sup`, in two places
  - Xj
  - G2_sup

  `ell` uses the Euclidean Delta. `kernel_mass` and `kernel_pieces` are geometry-free.
- **Planar network, my own integrator (`plan.py`, x-formulation, LSODA rtol 1e-9, Gram and gradients written independently of `simulate2`), four targets of my choosing**, with cert from `cert2d`. In every informative case the formed set equals the limit's, no dormant coefficient ever exceeds theta, the level-w1 crossing lies inside [-S^-, S^+] for every entrant, and the final error is <= 3e-9.

  | target (K, sep, R, seed) | Delta, R, mu, gap, J | cert true: eps, D | dual dev / D | notes |
  |---|---|---|---|---|
  | (2, 7.5, 3, 111) | 7.85, 3.0, 0.023, 0.035, 113 | 10^-7.5 (701), 10^-8.5 (674) | 0.14, 0.12 | cert false at 10^-6.5 and 10^-5.5, network still follows |
  | (2, 10, 1, 114) | 12.2, 1.0, 0.019, 0.054, 108 | 10^-7.5, -8.5, -9.5 (D about 125) | 0.31 | d_0 = (0.21, 0.50), a large offset |
  | (3, 9, 1.3, 102) | 9.4, 1.3, 0.022, 0.023, 170 | 10^-8.5 (257), 10^-9.5 | 0.30 | cert false at 10^-7.5 |
  | (2, 10, 1, 104) | 12.5, 1.0, 0.0057, 0.0055, 115 | 10^-10, 10^-11 (D 114) | 0.36 | near-tie, cert false at 10^-8 |
- **Interval arithmetic (claim 3).**
  - The chain is contiguous: consecutive balls share the same float endpoint `10.0**cur`, so the union is an interval. The candidate D is only verified, so soundness does not depend on how it was found.
  - The bisection inversions use three-valued comparisons and are sound. The ball version of s(x) via g(y)/tanh x (series for y < 0.5, alternating bounds) matches s(x) = (a/(2 G vs))(kappa atan(1/(kappa tanh x)) - x).
  - T1 and A_tr: my `acb.integral` gives 44.667 and 6.381. The author's upper Riemann sums give 44.679 and 6.388 (upper, as intended).
  - Pre-drift: the author's mesh is finer and rigorous; my coarser monotone sums converge to it from above.
  - The mass identity Psi/inv, the peak (sum of per-kernel peaks), the drift loop and the order condition are consistent with my transcription. Near the top the eps-ball must be narrow because the window starts and the gated drift terms depend sharply on eps (mass overestimated about 100 times at width 0.06 decade). That is conservative, not unsound.
  - Q-table consistency: the stated tops match the table (seven targets one half-decade below the float eps_0, one equal).
- **Lemma 12.4 (claim 4).** Re-derived (finding 7), Adamovic-Mitrinovic checked (the series of 3 ln(sin t/t) - ln cos t starts t^4/15, and a grid of 2e5 points), and the lemma verified numerically at 200 bits.

## Could not check

- Enclosure of the limit object (Krawczyk/Newton) and its decisions: critical points, entry times, entrant/feature assignment. Not done by the author either.
- K > 6 against the network; K=4..6 and planar `cert` in ball arithmetic.
- The K=2 planar run with the 6.5 sigma dictionary (too slow; killed).
- A full planar window, the sampled-Gram edge of the repository's network, and whether the repository's planar code really uses a scalar centre mobility and the law (12, 1, 0.1, 1).
- The remaining six Q targets: I re-ran only the chains of (2,6,1,11) and (3,8,1,31).
- Whether `hybr` found the intended critical point for K>=4 and in the plane. I checked this for the Q targets only.
- A proof of monotonicity in eps of the layer constants: still open (section 12.4 is honest about it).

## What the paper may honestly claim (exact wording)

> *Theorem 1, computer-assisted statement in a controlled regime. Setting: noiseless data; equal masses; Gaussian features of the kernels' width; the continuum Gram; the P0 unit with (m,G,gamma,a_max) = (12,1,0.1,1), so centre mobility G m/sigma^2 = 833 (scalar); n = 1, or n = 2 with the dictionary restricted to the lattice kernels within 3.5 sigma of a feature; a finite dictionary with entrants at offsets d_0 < r; K features at mutual distance >= Delta_T1(K,R,d_0) (n=1, d_0 = 0.42 sigma: 4.9 to 7.5 sigma for K <= 3, R <= 10, 5.5 to 7.6 sigma at K = 4, 6.9 to 8.5 sigma at K = 64; n=2, d_0 = 0.59 sigma: 5.2 to 7.9 sigma for K <= 3); the limit object satisfies its finite checks (a), (b) and the block critical points lie in the relaxation set (checked in floating point). (i) Proved by hand: if the explicit finite system S(eps) holds, then on [0, lambda(1.05 tau_K + 1)] the flow forms exactly the K units of the limit object, in its order, with entry times within O(eps^{1/2}) in tau = t/lambda, inactive duals within D_j/lambda = O(eps^{1/2}) of the limit's, dormant coefficients bounded by an explicit core bound, and final fit within O(eps^{3/2}). (ii) Computation: S(eps) was evaluated in floating point on 39 targets (K = 1 to 6, n = 1, 2), with no counterexample among the 473 grid pairs for which S held (all at least 1.5 decades below the measured threshold; final formed set checked, and the ingredient inequalities on 48 trajectories); and in interval arithmetic, with floating-point limit-object data, for every real eps in [10^-24, eps_top] (eps_top = 10^-5.5 to 10^-8.5) on 8 targets with K <= 3, n = 1. S is true for fixed data at all sufficiently small eps by a convergence argument, with an uncomputed threshold. (iii) Not proved: validity of S for all eps <= eps_0 (monotonicity), Delta below Delta_T1, unequal masses, noise, the sampled window's edge and a full planar window, eps below 10^-24, enclosure of the limit-object data, and anything beyond the horizon or for matrix centre mobilities.*

Overclaims to remove from 12.5:
- the "771 pairs" count (finding 1);
- "interval arithmetic ... proved" without the float limit data and hypothesis checks (finding 2);
- "a finite dictionary" without the 3.5 sigma truncation, and a single Delta_T1 range for n=1 and n=2 (finding 3);
- "13 to 56 times for K <= 3" (finding 10);
- "favourable sign" for mass/drift from Lemma 12.4 (finding 7);
- "trained network": the object is the continuum-Gram flow.

Still open:
- monotonicity of S in eps
- Delta below Delta_T1
- unequal masses and noise
- the sampled window's edge
- a full planar window
- t to infinity
- enclosure of the limit object
- matrix centre mobility
- K > 6 and n=2 in interval arithmetic
