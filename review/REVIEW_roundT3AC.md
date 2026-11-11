# REVIEW_roundT3AC: referee report on the small-eps tail (theory3 A) and the enclosed limit object (theory3 C)

Archived verbatim from the referee's hand-back (2026-11-10). The referee left the repository untouched; scratch code and logs are outside it.

REFEREE REPORT, round T3 (A = small-eps tail, C = enclosed limit object). Repository untouched; all scratch code and logs are in referee_T3AC_work/

VERDICT
- Defect (f) is real, and A's repair is correct but incomplete: a second, smaller omission (g) sits in the same term in every verifier.
- C's verifier uses the legacy accounting, so "S(eps) verified, 8 of 8" is not valid as committed.
- When I patch both terms into C's verifier and re-run its full covers, all 8 close on [1e-24, eps_top] with unchanged tops.
- A's tail argument is sound in structure; one stored tail certificate (K=2, 6/1/11) fails once (g) is added and passes with re-solved D_j.
- A and C cannot yet be merged into "every real eps in (0, eps_top] with enclosed limit data".

FINDINGS

1. BLOCKER (for C's claim as committed): C verifies the legacy system, not the S of Lemma 5.2.
- Evidence: `experiments/theory3c_rigorous.py` line 587 has `dcb = muc * (Hst[st] * Ig + hinc)` with `Ig` cumulative from t = 0 (lines 554, 583, 594). This is the same expression A flagged in `theory_t1b_rigorous.verify`.
- I patched it at import (A's three replacements, plus (g)) and re-verified all 364 stored balls of `theory3c_raw.npz` (`exp2all.py`, `exp2all.log`):
  - legacy: 364 of 364 pass;
  - repaired: 362 of 364 pass.
- The two stored balls that fail (`exp9.py`, failing for every candidate factor 1.1 to 4):
  - (3,8,1,31), top ball [10^-8.75, 10^-8.5]: E/D 0.951 becomes 1.069 at factor 1.5;
  - (2,8,1,13), ball [10^-8.1875, 10^-8.0625]: E/D 0.933 becomes 1.003.
- Re-running `R3.cover` from the stored tops with (f)+(g) (`exp2.py ... fullg`, `exp2_fullg.log`): 8 of 8 complete down to 10^-24 with the same tops.

  | target | balls | verify calls | narrowest ball (decades) |
  |---|---|---|---|
  | 0 | 19 | 37 | 1/4 |
  | 1 | 86 | 443 | 1/64 |
  | 2 | 37 | 172 | 1/32 |
  | 3 | 43 | 190 | 1/32 |
  | 4 | 22 | 66 | 1/16 |
  | 5 | 20 | 46 | 1/8 |
  | 6 | 18 | 40 | 1/8 |
  | 8 | 122 | 674 | 1/256 |

- Fix: apply the two patches in `verify_rig`, regenerate `theory3c.md` and its raw file, and add a test that legacy and repaired `dcmax` differ. About 5 CPU minutes.

2. MAJOR: defect (f) is real (re-derived).
- For a non-entrant, c_l' = mu_c w_l h_l with h_l = ghat^(st(t))'(c_l) + d_c R_l. Hence |c_l(t) - c_l^0| <= mu_c ( int_0^t H^(st(t')) B dt' + peak_l int_0^t |d_c R_l| ). The weight H belongs to the time at which the mass was accumulated.
- The legacy form H^(st(t)) * int_0^t B is not a bound. In the last stage H^(K) is about 0.43 a_min Omega_pk, close to 0, so legacy forgets every earlier displacement. S2 then takes the maximum of stage-end values where the lemma gives their sum.
- H is not monotone in the stage: on planar target (2,8,1.5,2) the legacy value is larger (ratio 0.96). So legacy is neither an upper nor a lower bound.
- Size on the float `cert`, all 473 informative pairs (`exp1.py`, `exp1.log`):
  - dcmax at each target's eps_0 changes by a factor 0.96 to 1.90 (median 1.35); the largest over all pairs is 1.98;
  - E_j rises by at most 8e-4 relative.
- Size in the ball verifiers on wide balls near the top: max drift +30 to 40 %, E_j up to +12 %.
- A's note is wrong on two details:
  - "for K = 1 legacy and repaired dcmax coincide": I measure factors 1.58 to 1.90 on the four K=1 targets, because window 1 starts exactly at t_1 and its hd term is absent in stage 0;
  - "dcmax +20..50 %, E +1..4 %" understates both.
- Effect on the earlier float certificate: 2 of 473 pairs flip, both at the target's eps_0 and both on S2 (d_c = 0.3), with no parameter set passing.
  - (2,8,1,13) at 10^-7.5: dcmax 0.240 becomes 0.341; eps_0 moves to 10^-8.0.
  - (6,8,1.5,3) at 10^-8.5: dcmax 0.299 becomes 0.323; eps_0 moves to 10^-9.0.
  - Counts: K<=3 1 of 228; K=4..6 1 of 204; planar 0 of 41. With (g) added the same two flip (`exp1g.log`).
- Consequences: "473 informative pairs" becomes 471; two table entries for eps_0 are half a decade optimistic, and those two ratios eps_max/eps_0 grow by 3.2. "0 violations" is unaffected. Planar (3,8,1.5,2) sits at dcmax 0.2957 against 0.3.
- Fix: patch `cert` (and hence `cert2d`), regenerate the two eps_0 values and the counts.

3. MAJOR: a second omission (g) that A did not find, present in `cert`, `theory_t1b_rigorous.verify`, A's chain, A's tail and C.
- The window gating of hinc counts unit i only after wstart_i = t_i - S_i^-. Before that, the entrant's pre-layer coefficient (up to w_1 = 0.035 at delta_1 = -15) still enters d_c R_l, with int v_i <= pre = 0.60.
- The gated bound therefore drops psi_sup(X_li) * pre for each window not yet started. Lemma 5.2 as written is ungated.
- Size:
  - my independent evaluator: max drift +0.5 to +1.3, worst E/D +0.2 to 0.5 % near the tops;
  - float `cert` at eps_0: E_j up to +4.3 % (K=1 seed 1: 93.22 becomes 97.21).
- A's tail with (g) at the stored D_j (`exp6.py`, `exp6.log`): 7 of 8 still pass, with E/D rising 0.0002 to 0.007. Target (2,6,1,11) fails: E/D 0.9980 becomes 1.0034. `solve_D` at the same eps_tail then passes with worst ratio 0.883.
- A's 64 stored chain balls all pass with (f)+(g) (`exp7.py`).
- Fix: add the psi_sup * pre term for windows not started (one line per verifier) or drop the gating; regenerate.

4. MAJOR: A and C are not compatible as they stand.
- A alone: for every real eps in (0, eps_top] on the 8 targets, a system S with (f) repaired but (g) missing. Limit data are floats taken as exact, hypotheses are checked on floats (final block only), comparisons go through float(ball.upper()), and the chain inherits the T1b ball gating sorted by midpoints. One tail certificate must be regenerated (finding 3).
- C alone, as committed: enclosed limit data and certain comparisons, but the legacy displacement term, on [1e-24, eps_top]. After the patches (my re-run) it gives the repaired S on [1e-24, eps_top], same tops.
- Missing for "(0, eps_top] with enclosed limit data": the interval (0, 1e-24). A's tail reads `prep` floats throughout (`dormant_pieces`, `kernel_margin`, p_e, vhat, `rgs`, `RG.hypotheses`).
- A merged verifier needs:
  - `verify_rig` with (f)+(g);
  - a port of `_verify_tail` to the ball limit object LR (pieces, c_min, p_e, vhat and `rig_rgs` from LR; `hypotheses_rig`; certain comparisons);
  - one tail ball at any eps_hi >= 1e-24, for example 1e-12, where A's max E is 76 to 446 with large slack. No new chain is needed.

5. MINOR: Theorem TAIL's h-ball argument is sound, modulo (g) and the float inputs.
- I found no quantity evaluated at a single eps. `best_cert(10^e)` supplies only pi and the candidate D.
- Checked line by line:
  - h in [0, h_hi] encloses every eps in (0, eps_hi], and the bisections return enclosures of all roots over the ball;
  - d eta/dx is bounded by its eps = 0 form, which is valid because q >= 2 v_s tanh^2 x;
  - env(c) c^3 is decreasing in c, and kap_j is taken at (a_lo, L_hi);
  - the integral of c^-3 over a linear piece is exact;
  - Omega_pk * T_fin = tau_end L^2 sum e_j / a_min;
  - gates use the earliest start tau_i - S_i^- L_hi;
  - the displacement is stage-resolved (Wq);
  - S4 is checked with L_hi.
- The chain and the tail meet at the double 10.0**etail.
- Residual, also in C: `kappa = k1 if up(k1) < up(k2) else k2` selects by upper ends, which is not a lower bound of the minimum if the balls overlap. Here k1 is about 0.06 and k2 about 12, so it is harmless. Fix: take the smaller lower end.

6. MINOR: every verifier uses the final-block coefficient vhat_fin for an old unit i during window k. Lemma 5.4 needs the limit coefficient of the current stage. The stage value exceeds the final one by up to 4.8e-5 relative on (2,6,1,11) and 8.3e-6 on (2,6,2.4,12), where C's cover has max E/D 0.9999. Fix: use the maximum over stages, or a_min v_D^+.

7. MINOR: Lemma E and Lemma M are correct; Proposition N is not a proposition.
- Lemma E, re-derived: kappa tanh x = cot(theta + x/kappa), theta = s/(sigma_0 kappa), the two terms sum to at most A rho^(3/2) (1 - s_eps^2)^-1 cos/sin^3, and A rho^(3/2) kappa^3 = 2 v_s a_max / m.
- The Adamovic-Mitrinovic step is right: 3 cot t + tan t - 3/t = sum 2^(2n) |B_2n| / (2n)! * (2^(2n) - 4) t^(2n-1); the coefficient is 0 at n = 1 and positive for n >= 2 (4/15, 8/63, 4/75, ...).
- The statement should say s <= lambda; beyond it B <= 0 and the bound is trivial.
- Proposition N: only the numerator of pre is proved monotone. The directions of w_1, T_3, Env, kappa_*, m_V are numerical (2 targets, 14 eps). "S(eps_1) => S(eps) is false in that form" is unsupported: no counterexample exists. What is shown is that the ingredient-wise route fails.

8. MINOR: C's Krawczyk implementation is correct.
- K(X) = z0 - C F(z0) + (I - C J(X))(X - z0), with J evaluated on the box, C the midpoint inverse, and a strict interior test in arb.
- The D-hat uniqueness is a single Krawczyk step on the exact product box D-hat, with no subdivision and none needed: my own Jacobian gives ||I - C J(D-hat)||_inf <= 0.015 to 0.073 for every stage of targets 0, 3, 6, 8 at three parameter sets, and my K(X) lies in int X. The test runs inside every `verify_rig` call, so it follows the parameter set.
- Two small points:
  - the entrant candidate is chosen by float midpoint, so a 1-ulp near-tie is refused although decidable (conservative);
  - `inv_of` contains a float gate. It only chooses between two valid bounds, but "no float() of a ball in any decision" is literally untrue.

9. OVERCLAIMS in the proposed wording.
- A:
  - "ball-arithmetic proof ... for every real eps" must carry "floating-point limit data and float comparisons";
  - "(E) defect ... repaired" is incomplete (finding 3);
  - "all 8 covers are re-done" is really 64 balls on [eps_tail, eps_top] for 5 targets;
  - the "K = 1 coincide" and "+20..50 %" statements (finding 2);
  - "Proposition N, the exact obstruction / false in that form" (finding 7).
- C:
  - "Fully rigorous", "8 of 8 targets: S(eps) is verified", "tops unchanged": not for the S of Lemma 5.2 as committed (finding 1);
  - "no float gates" (finding 8);
  - "all hypotheses on enclosures" omits finding 6 and Theorem K(a), which C's own gap 1 concedes.
- Both: "Theorem T1 proved on 8 targets" needs "computer-assisted; the implication from S(eps) is a hand proof; horizon [0, lambda(1.05 tau_K + 1)]".

VERIFIED (how)
- Lemma 5.2 accounting: re-derived by hand; legacy and repaired paths located in `cert`, `RG.verify`, `verify_rig`, `_verify_tail`.
- Flip counts: `cert` and `cert2d` patched at import and run on all 473 informative pairs, with (f) and with (f)+(g).
- C's verifier: all 364 stored balls under three accountings, and full re-covers for 8 targets with (f)+(g).
- A's tail: 8 stored certificates with and without (g); A's 64 stored chain balls with (f)+(g).
- Independent re-evaluation of S(eps) (`indep.py`, about 200 lines written from THEORY_T1.md sections 2 to 6 with scipy quad/brentq, stage-resolved displacement, stage-maximum vhat):
  - targets (2,6,1,11) at eps = 10^-7.501, 10^-9, 10^-13 and (3,8,1,31) at 10^-8.501, 10^-10, 10^-14, using A's stored tail D_j;
  - every condition holds, and my E_j is below A's majorant for every kernel;
  - T_3, Env, pre and kappa_* agree to 6 digits; A's T_1 upper sum is 0.012 above my quadrature.
- Lemma E: 58,120 random (eps, x) pairs with eps in [1e-40, 0.05), in 400-bit arb; 0 not certainly below the envelope; smallest certified gap 3e-26. Series identity checked to 1e-16 with exact Bernoulli coefficients.
- Krawczyk: 12,000 random boxes on three one-unit, two-feature systems; 0 accepted among 4,862 boxes with two or more roots and 5,847 with none. Ties are refused: an exact tie and two mirror-symmetric targets give ok = False; a tie broken by 2^-40 is decided.
- `tests/theory3a_test.py`: 16 passed. `tests/theory3c_test.py`: 38 passed.

COULD NOT CHECK
- The hand proofs of sections 2 to 6 beyond Lemma 5.2's displacement term and Lemmas E and M. I examined only that term, and (g) suggests other transcription gaps may exist.
- The mutant records of A and C; byte-identical regeneration of either report; arb's own correctness.
- A's chain re-generated from scratch with (g): I re-verified the stored balls only.
- An independent re-implementation of C's entry recursion: I ran tie and decision tests only.
- Interval verification for K >= 4 and the plane: floats only.
- Whether the float limit data could change A's verdicts below 1e-24: not rigorous either way.

SENTENCE THE PAPER MAY USE
Allowed once C is regenerated with the two corrected terms (my scratch re-run closes 8 of 8 at the same tops) and A's tail certificate for (2,6,1,11) is regenerated:

"For eight test targets (n = 1, K <= 3, continuum Gram, equal masses, noiseless data, law (12, 1, 0.1, 1)) the finite system S(eps) of Theorem S-T1 has been verified in ball arithmetic: with Krawczyk-enclosed limit-object data (existence and uniqueness of each stage block's critical point in the relaxation set, strict entry decisions) for every real eps in [10^-24, eps_top], and, with floating-point limit-object data taken as exact, for every real eps in (0, eps_top], where eps_top = 10^-5.5 to 10^-8.5. Combined with the hand proof that S(eps) implies the conclusions on the horizon [0, lambda(1.05 tau_K + 1)], Theorem T1 holds in this computer-assisted form on these targets for eps in [10^-24, eps_top]; below 10^-24 it rests on unenclosed limit data. Not verified: that the enclosed critical points are the global minimisers of Theorem K(a), and any monotonicity of S in eps."

Until regeneration, only the second half may be used: "with floating-point limit-object data, for every real eps in (0, eps_top]". Also correct the two eps_0 entries and replace "473" by "471" (finding 2).
