# REVIEW_roundT1: referee report on THEORY_T1.md (Theorem T1, explicit finite certificate S(eps), eps_0)

Scope. I read THEORY_T1.md line by line against Lemma S1 and the supplement text, THEORY_NOTE/REVIEW_roundTH, experiments/exp_theory_t1.py and results/synthetic/theory_t1.md. I did not use the author's functions for the checks below. I wrote my own code in /tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_T1_work/:
- mybox.py: box and regime (R) in python-flint ball arithmetic.
- unit.py, s1.py, bs3.py: the unit and Lemma S1.
- myt1.py: target, limit path and network integrator.
- mycert.py: my own S(eps), transcribed from Sections 2-6 of the note.
- mono.py, mono2.py: fine-grid monotonicity of S.
- adv.py: random-target test of S against the network.
- run1.py: runs at eps below, at and above eps_0.
- mobil.py, hess.py: mobility sensitivity and Hessian sampling.

**Bottom line.** I found no mathematical error in the conditional statement "S(eps) true implies the network follows", in the stated regime. The weighted-max-norm argument (Lemma 4.1), the exact dormant-core algebra (Lemma 2.2), the scalar-comparison transit bound (Lemma 3.3) and the triangular bootstrap are all correct. My independent S reproduces the author's D and eps_0 on 5 targets. 40 new random targets gave no violation. The problems are what is claimed:
1. "for all eps <= eps_0" is not proved.
2. "Proved" rests on floating point.
3. The theorem is for the continuum Gram, not for the sampled network.
4. S depends on the centre mobility.
5. Several sentences in section 10 are stronger than the proof.

## Findings

**1. [blocker for the text; major as mathematics] "S(eps) holds for all eps <= eps_0(target)" with an explicit eps_0 is not proved.**
- Location: Sec. 0 bullet 2, Sec. 7 Cor. 7.1, Sec. 9 item 3, Sec. 10 (Theorem 1 status and Abstract sentences).
- What is proved: Cor. 7.1 gives the existence of some unspecified eps_1(target) > 0. eps_0 is a grid value with three refinement points. Monotonicity in eps of T_1, A_tr, pre, Env, kappa_*, the weights and the drift is not established (the note admits this in Sec. 9).
- Evidence for the rest of my claims here:
  - Evaluating my own S on a 0.1-decade grid from eps_0 down to 1e-20, on the five targets (K,sep,R,seed) = (1,-,1,2), (1,-,1,1), (2,6,1,11), (2,8,2.4,14), (3,8,1,31), gave no failure. D decreases monotonically and saturates: 83.3 -> 76.0, 93.2 -> 85.4, 110.1 -> 108.9, 426 -> 421.9, 132.2 -> 130.4.
  - So the claim is very probably true, but this is evidence, not proof.
  - Separately, the half-decade grid makes eps_0 conservative by up to a factor 2.4 beyond the true certificate threshold. With steps of 0.125 decade, S is true from 10^-5.125 for target (1,-,1,2) (grid eps_0 is 10^-5.5), from 10^-6.75 for (2,6,1,11) (grid 10^-7), and from 10^-6.75 for (3,8,1.5,15) (grid 10^-7).
  - The headline "factor 32, 100, 316" therefore overstates the gap to the measured eps_max by about 2-3x (13, 18, 56 instead of 32, 32, 100). Both eps_0 and eps_max are grid values, so "median 32" is the grid quantum of three half-decades, not a statistic.
- Cheapest fix: state T1 as "S(eps) true implies ...". Then state separately: "S(eps) was verified at the listed eps (list the grid); a convergence argument gives some eps_1 > 0 below which it holds, not computed". Or prove monotonicity (the constants are explicit integrals of the unit, so this is plausible). Stop saying "explicit eps_0 for all eps <= eps_0".

**2. [major] "Proved" rests on double precision.**
- Location: Sec. 0 not-proved item (5), but Sec. 10 and the Abstract say "proved".
- Items in S that are not rigorous:
  - 0.1% factors on the trapezoid quadratures T_1 and A_tr (8001 points).
  - Inversion of s(x) by np.interp on an 80000-point log table (x_of_s).
  - Sup of |g-hat'| and oscillation of g-hat from 401-point grids plus a Lipschitz margin. The oscillation formula uses the grid max of |g-hat'|, not its sup.
  - The fixed point with the 1.0001 factor; `<=` tests in floating point.
  - Peak of the sum over dormants, sampled in time (in the code, `sumB` at the breakpoints).
- Cheapest fix: label T1 "computer-assisted, floating point, not interval arithmetic" everywhere it is quoted. Optionally recompute the closed-form integrals in flint. Everything in S is a monotone closed-form integral, so this is feasible for K <= 3.
- Reassurance: my independent transcription of Sections 2-6 gives D = 83.3, 110.1, 449.5, 132.2, 100.6 at eps_0 = 10^-5.5, 10^-7, 10^-6.5, 10^-8, 10^-10 on five targets. The author's values are 83, 110, 450, 132, 101. So the code equals the note. This does not prove the note is right.

**3. [major] The theorem is for the continuum Gram. The benchmark network is a sampled window.**
- Location: Sec. 0, Sec. 10 last sentence, Sec. 8.E.
- Part E (184 pairs, 8 interior targets) runs `train_flow` with the continuum constants in `cert`. It is empirical, not covered by T1.
- The Sec. 10 exclusion list ("smaller separations, planar, unequal masses, noise") omits "sampled/finite window".
- REVIEW_roundTH finding 1 (edge effects of the sampled Gram) still applies.
- Fix: add "continuum Gram" to the hypotheses in the Theorem 1 sentence. Call Part E "observed, certificate with continuum constants".

**4. [major] The certificate depends on the centre mobility. "Steady state contains no mobility" is true, but the theorem is not mobility-free.**
- Lemma 4.1: V_infinity indeed contains only beta_v and beta_d. Row-wise dominance makes the max-norm contraction independent of the diagonal mobilities. This part is correct.
- But the dormant drift (Lemma 5.2) and the pre-drift (Lemma 5.3) are proportional to mu_c. T_2 is proportional to 1/mu_c. kappa_* contains M and mu_c.
- Sensitivity of S to mu_c (mobil.py), with mu_c = 833 x factor:

  | Target, eps | factor 0.01 | 0.1 | 1 | 10 | 100 |
  |---|---|---|---|---|---|
  | (1,-,1,2), 1e-6: D | 75.6 | 74.8 | 79.5 | 119.5 | 412 and S fails |
  | (2,6,2.4,12), 1e-7: D | 436 | 434 | 447 | 499.5 | 800 and S fails |

- So the result is for the P0 law with centre rate Gm, mu_c = Gm/sigma^2 = 833. A different bandwidth or other per-group learning rates move eps_0, and above about 10^2 x the certificate fails.
- Fix: put "mu_c = Gm/sigma^2 = 833" and "(m,G,gamma,a_max) = (12,1,0.1,1)" into the theorem's hypotheses. Say that mobility enters only through the layer durations and the dormant drift.

**5. [major, wording] Section 10 sentences stronger than the proof.**
- (a) "dormant coefficients bounded by 4/(G gamma^2 m) s^-3": Lemma 2.2 gives only the mu^-3 forms (constant 0.373 for the crude one, 0.0963 as the leading term). Cor. 2.3 is an asymptotic limit (eps -> 0, s >> 20). The proof uses the exact B(s). Numerically sup_s B(s)s^3 <= 33.3333 for every eps tested (4e-2 to 1e-16, bs3.py). The sup is approached from below as eps -> 0. At eps = 1e-4 the maximum is 32.63 at s = 242; at 1e-12 it is 33.33333. This is probably provable (a one-variable inequality for B_0, with an eps-correction argument) but is not proved. Call it "observed" or "leading order".
- (b) "(H1)-(H4) are consequences": (H4) is only for levels <= a_min v_lo. (H1) needs the finite checks (a),(b).
- (c) "factor 32-316 (median 32)" omits the two near-tie targets (10^4 and 10^6). It is selective, and the figures are grid-quantised (see finding 1).
- (d) "no violation in 529 pairs": these are 15 + 8 targets on one grid, not independent.
- (e) Horizon: all conclusions are on [0, lambda(1.05 tau_K + 1)]. The bound grows like c_N Omega' T_fin, so nothing is claimed as t -> infinity. Say so.
- (f) The listed hypotheses lack: r > d_0^k; (a),(b) of Thm K; the finite dictionary (its dormant count and positions enter D and eps_0); eps_0 depends on the target through mu_min and the entry-time gap, with no uniform bound (eps_0 = 10^-10 for gap 3e-4).

**6. [minor] Lemma 3.1, top face.** "epsilon g <= a_min(R + E_v - vbar) < 0 (definition of vbar)" is "= 0" with vbar = (R(1+kappa_c)+Omega)/(1-kappa_c) exactly. supp_entry.tex has the factor (1+0^+); the note dropped it. Strictness holds at all non-extreme points. At the extreme corner (a_i = R a_min, d_i = 0, all other bounds attained) the Nagumo margin is zero. Fix: restore (1+0^+).

**7. [minor] Proof of Thm T1: unwritten steps.**
- T* must also include the first exit from the ball |c_l - c_l^0| <= d_c. (S2) is computed from bounds that presuppose that ball.
- "Units i > k are still in the pre-phase at t_a^k" is used in D^e_k (the sum i >= k of pre). It needs the first-time argument, which closes by (S4): take the first reaching time among units >= k.
- The claim v-hat in [V_t, v_up] (used in Env) follows from D* (v >= Gamma(d-tilde) - E_v > V_t), not from D-hat. Say so.
- None of these is a gap in substance. There is no circularity: layer k uses only layers < k, plus the inputs D_l.

**8. [minor] v-hat evaluated at the final block in c_oldw and c_oldtail** (vhat_fin[i]), where the stage-(k-1) and stage-k blocks are meant. Immaterial: the stage and final coefficients differ by <= 1e-5 relative on the 15-target set (blocks sit at the features to 1e-4 at Delta >= 6). At Delta near Delta_T1 it is up to about 1e-3. Fix: use the max over stages.

**9. [minor] Lemma 4.2 (Lojasiewicz/PL) is correct but unused.** It is unforced. The transit bound T_1 and the settling bounds T_2, T_3 are scalar comparison arguments (d eta/dt >= a_min v_lo - W(eta); Y' <= -M Y). They are valid, and monotonicity of M in x gives the m_V constants. Sec. 0 (M3) should say "scalar comparison", not PL.

**10. [minor] Regime facts.**
- Gershgorin margins beta_d are small at Delta_T1, so 1/beta_d is large in c_F, c_N and V_infinity: (K,R,Delta) = (2,1,4.9): 0.10 at the best r; (2,2.4,5.7): 0.076; (2,10,7.0): 0.128; (3,1,5.3): 0.107; (3,2.4,6.1): 0.141.
- Sampled true Hessian row margins (hess.py) are 0.34-0.48 against bounds 0.02-0.09, so the bounds are valid and loose.
- Condition (N) never binds: Delta_T1 is identical with and without it. Gershgorin implies (N).

**11. [minor] eps_max/eps_0 for the near-tie targets (10^4, 10^6) is a feature of (S4), not of the network.** The network follows to 1e-4 (re-simulated). The certificate needs disjoint windows. State that eps_0 is only sufficient and not target-uniform.

**12. [minor] Table 8.D "dormant mass ratio 1.00"** is trivially tight at tiny eps (dist - D ~ dist), not saturation. A one-clause explanation would help.

**13. [limits, not errors] n = 1, equal masses b_j = 1, K <= 3.** The row sums and own block were checked for n = 1 only. For n > 1 the dormant count and tangential eigenvalue change. Unequal masses change the walls lambda_j. These are correctly flagged in Sec. 9; keep them in the theorem statement.

## Verified (re-derived by me; how)

- **Lemma S1.** Integrated alpha' = G beta a' g, beta' = gamma beta a g with a prescribed g(t) (DOP853, 1e-12): the integral of g equals eta(alpha) to 1e-11 (eps = 1e-2, 1e-3). beta = eps + 2 v_s sinh^2 x to 1e-15, and w = beta a/m = W(eta).
- **Lemma 2.1 and 2.2.**
  - Algebra by hand: kappa tanh x = cot(vartheta + x/kappa) with the argument in (0, pi/2]; tanh x <= s_eps; the constants 0.0963 = (8/pi^3) a/(m sqrt(2 v_s)) and 0.373; c_u = pi a_max/(2 sqrt2 G sqrt(v_s)) = 7.0248.
  - Numerically in ball arithmetic, 120-bit precision, B(s) by bisection of eta(x): the chain B(mu*lambda) <= first form <= cos/sin^3 form <= mu^-3 form holds in all 36 admissible (eps, mu) combinations with eps from 1e-2 to 1e-12 and mu from 0.01 to 0.9, for every pair with s_eps < 1.
  - B s^3 reaches 33.3333 at eps = 1e-12.
- **Lemma 4.1.**
  - The comparison Phi = V_infinity + c e^{-kappa t} is valid because V_infinity is nondecreasing.
  - The Dini-derivative step is valid because the contraction is row-wise.
  - The row sums (rowv, rowd) agree with my own count of the cross terms of <e_i, e_j>, including the (v_i,d_i) cross entry and the unfitted features (K-1 terms).
  - Finite-difference Hessian sampling of F_0 (1500 corner-biased points for each of 5 (K,R,Delta) cases) gives row margins well above beta_v and beta_d.
- **Lemma 3.3.**
  - T_1: eta' >= a_min v_lo - W(eta); pre via the mass identity (Lemma 2.4).
  - T_2: |d| - d-tilde decays at kappa_d.
  - T_3: Y = V_t - v satisfies Y' <= -M(v) Y, and v >= V_t is preserved; M increasing in x gives m_V.
- **Delta_T1, ball arithmetic, own code.**
  - Smallest Delta with some r > d_0 (r-step 0.005, Omega = 0, varsigma = 0.01):

    | K | R | d_0 | Delta_T1 (mine) | author |
    |---|---|---|---|---|
    | 2 | 1 | 0.4167 | 4.76 | 4.9 |
    | 2 | 2.4 | 0.4167 | 5.61 | 5.7 |
    | 2 | 10 | 0.4167 | 6.85 | 7.0 |
    | 3 | 1 | 0.4167 | 5.17 | 5.3 |
    | 3 | 2.4 | 0.4167 | 5.94 | 6.1 |
    | 3 | 10 | 0.4167 | 7.10 | 7.2 |

  - I get 4.32 (K=2, R=1, d_0=0.2) and 7.27 (K=2, R=10, d_0=0.6). All 18 author entries are >= mine (the conservatism is r-grid and 0.1 rounding).
  - Best r is about 0.422.
- **Certificate reproduced, own code.**

  | Target | eps_0 mine | D mine | author eps_0 / D |
  |---|---|---|---|
  | (1,-,1,2) | 10^-5.5 | 83.3 | 10^-5.5 / 83 |
  | (2,6,1,11) | 10^-7 | 110.1 | 10^-7 / 110 |
  | (2,6,2.4,12) | 10^-6.5 | 449.5 | 10^-6.5 / 450 |
  | (3,8,1,31) | 10^-8 | 132.2 | 10^-8 / 132 |
  | (2,6,1,21) near-tie | 10^-10 | 100.6 | 10^-10 / 101 |

  - T_lay = 117.4-117.8 (author 117).
- **Simulation (own integrator; LSODA, rtol 1e-10, x-formulation; same eps_max as the author on 4 targets).** Formed set == limit, and no dormant ever above theta, down to 1e-7.5 at least:

  | Target (K,sep,R,seed) | follows at | fails at | eps_max | E_act at eps_0 (lambda*mu/E_act) |
  |---|---|---|---|---|
  | (2,6,2.4,12) | 10^-4.0 | 10^-3.5 | 10^-4 | 47.9 at 10^-6.5 (21.6) |
  | (3,8,1.5,15) | 10^-5.0 | 10^-4.5 | 10^-5 | 41.9 at 10^-7 (18.9) |
  | (2,10,1,23) near-tie | 10^-4.0 | 10^-3.5 | 10^-4 | 39.4 at 10^-7 (72) |
  | (1,-,1,1) | 10^-5.5 | 10^-5.0 | 10^-5.5 | 42.8 at 10^-6 (2.5) |

  - E_act is 37-69 over eps from 10^-7.5 to 10^-3.5 (the bound D at eps_0 is 83-450). It grows slowly with eps.
  - At eps_max, lambda*mu/E_act = 0.92-1.2 (e.g. 0.92 and 1.20).
  - Entry shifts, monotone entrants (until a_min v_lo) and offsets <= r hold in all runs below.
- **Adversarial test.**
  - 40 random targets (K in {1,2,2,3}, R in {1,1.5,2.4,3}, sep 6-12, random seeds), my S (the largest grid eps from 10^-5 to 10^-10 with S true), my integrator.
  - Checked: formed set equals the limit, no dormant ever above theta, E/D, entry shift within [-S^-, S^+], entrant monotone until a_min v_lo, offsets <= r, v <= vbar.
  - **0 violations**; max E/D = 0.523 (author 0.527). The factor-2 slack in D is the conservativeness.
- **Fine-grid S** monotone on five targets (finding 1).

## Could not check

- n > 1 and K >= 4; unequal masses; the sampled-window network (Part E: I did not rerun `train_flow`); the repository's tests.
- Rigorous (interval) arithmetic for the quadratures and the drift bound. I re-implemented the drift term from the code, not from an independent derivation of its constants (H^(st), osc, h_dev).
- A formal proof of monotonicity of S in eps; only fine-grid evidence.
- Lemma 4.2 in the forced case (not used).

## What the paper may claim, exact wording

**Theorem T1 (replace Theorem 1's status sentence in Sec. 3, regime stated).**

> "Controlled regime (computer-assisted). Continuum Gram ⟨φ_cφ_{c'}⟩ = ZΓ(c−c'), n = 1, noiseless data, equal masses; K ≤ 3 Gaussian features of the kernels' width, separation Δ ≥ Δ_T1(K,R,d_0) (4.8σ for K = 2 equal features, 5.7σ at ratio 2.4, 7.0σ at ratio 10; 5.2–7.2σ for K = 3), the P0 unit with (m,G,γ,a_max) = (12,1,0.1,1) and centre rate Gm, a finite dictionary, the limit object's finite checks (a),(b), and r > d_0^k. If an explicit finite system of inequalities S(ε) (Thm S-T1) holds, then on [0, λ(1.05τ_K+1)] the trained network forms exactly the limit object's K units and no other, in the limit's order, with entry times within O(ε^{1/2}) and inactive duals within D_j/λ = O(ε^{1/2}) of the limit's (explicit constants), and the final fit is the target to O(ε^{3/2}). S(ε) is evaluated in floating point (not interval arithmetic). For a fixed target it holds for all sufficiently small ε (a convergence argument; the threshold is not computed), and it was verified at every ε of a half-decade grid down to 1e-14 for the 15 test targets, with ε_0 = 10^-8–3·10^-6 (10^-10 for a near-tie)."

**Observed, not proved.**
- eps_max is a factor 13-56 above the true certificate threshold on three targets and 32-316 above the grid eps_0 (two near-ties 10^4 and 10^6).
- The sampled P0 network agrees on 8 interior targets (184 pairs).
- No violation was found by my independent test on 40 random targets.

**Remain outline or open.**
- Sampled/finite-window network.
- n > 1.
- Unequal masses and noise.
- Δ below Δ_T1.
- General (H1)-(H4) and the O(ε^{3/4}) statements of Thm 1 as in sec_theory.tex.
- Explicit eps_0 for all real eps ≤ eps_0.
- t → ∞.

**Remark replacing the K^3 formula.** Keep the Lemma 2.2 inequality with its hypotheses (|η| ≤ (1−μ)λ, s_ε < 1), and add "the distance form 4/(Gγ²m)s^-3 is the eps -> 0, s >> 20 limit". Add "(numerically an upper bound for all s)" only as observed.

**Abstract.** Replace "proved for ε ≤ ε_0 ... one to two orders below" with: "proved, by a floating-point finite certificate, in a controlled continuum regime; the certificate's ε_0 is one to two orders of magnitude below the measured threshold on 15 test targets".

**How conservative.**
- eps_0 is 13-56 below the measured eps_max by the finer threshold (32-316 by the grid).
- That is 3.6-7.5 in sqrt(eps).
- The cause is D versus E_act: E/D <= 0.52 in all 40 + 30 runs. The proof pays D + s_Omega where the network pays E_act, plus grid rounding.
- Gershgorin margins are loose by 5-20x but this barely enters (V_tails are 2-48 of D).
- The sharp prediction is lambda*mu = E_act: at eps_max, lambda*mu/E_act is 0.9-1.9 on all 15 targets.
