# theory3 N: corrected, merged verifications of Theorem T1 (work note)

Source: review/REVIEW_roundT3AC.md (binding). Code: experiments/theory3n_{kernels,system,tail,run}.py.

## The one corrected system S
Lemma 5.2 (displacement of a dormant kernel l):  |c_l(t)-c_l^0| <= mu_c ( int_0^t H^(st(t')) B_l dt' + peak_l * hinc(t) ),
 hinc(t) = sum_{windows started} (psi_sup(X_li) Ttv_i + G2_sup(X_li) Ttd_i) + sum_{windows NOT started} psi_sup(X_li) * pre + 0.43 * mass_tot.
 (f) stage-resolved first term (W = sum of H^(stage) x mass-of-piece); (g) the 'pre' term for unstarted windows (int v_i <= pre before its window);
 (6) vhat_fin[i] -> max over stages s>=i of |w_i^(s)|; kappa = the min with the smaller LOWER end; inv_of without float gate.
Implementation: anchored source replacement of exp_theory_t1.cert (float, n=1 and planar), theory3c_rigorous.verify_rig (balls); the tail verifier is a port of theory3a_tail._verify_tail
to the ball limit object (pieces, cmin, p_e, vhat, rig_rgs, hypotheses_rig), certain comparisons only.

## Tail ball
h in [0,h_hi], h_hi = (upper bound of sqrt(eps_meet/(2vs-eps_meet))) * (1+2^-40), eps_meet = fdown(10^EHI) = lower end of the chain's last ball.
Covers every eps in (0, eps_of_h(h_hi)], eps_of_h(h_hi) >= eps_meet is tested (certain).
Gates of the windows: earliest start tau_{i+1} - Sm_i L_hi rounded DOWN to a double; a window is 'not started' on a sub-interval only if gate >= the sub-interval's upper end.
Dormant c(tau)=1-|u| is concave on each piece, so the chord (also across an undecided zero crossing) is a lower bound: c^-3 integral by the exact chord formula is an upper bound.

## Status
See results/synthetic/theory3n.md. Gaps: Theorem K(a) (global minimiser), monotonicity of S in eps beyond the chain+tail, soundness of arb, hand proofs of S => theorem, horizon [0, lambda(1.05 tau_K+1)],
K>=4 and planar only floating point, sampled Gram.

## Results (2026-11-10)
* Float cert: 473 -> 471 informative pairs of 771; flips (2,8,1,13) 1e-7.5->1e-8.0 and (6,8,1.5,3) 1e-8.5->1e-9.0 (network re-run: follows, no counterexample); violations 0.
  dcmax factor corrected/old over old-informative pairs: median 1.25 (K<=3), 1.42 (K=4..6), 1.12 (planar); range 0.96-1.98. max E relative change 4.3e-2 (K<=3).
* Stored 364 balls of 3c: legacy 364, corrected 362 (the two the referee found: target 4 (2,8,1,13) and target 6 (3,8,1,31) top balls).
* Enclosed: 8/8 chains complete from the unchanged tops (-7.5,-5.5,-6.5,-7.5,-8.0,-7.0,-8.5,-7.0) to 1e-12, one tail ball (0,1e-12] with strict E<D/1.0001 (worst E/D 0.952-0.964);
  full 3c-style chains to 1e-24 also complete 8/8; tail alone reaches the chain top (-7.5 -> ... see T table).
* Tests 31; mutation specs 3 files; runs in $W/mut.
* Independent check with the referee's indep.py (gated='pre', stage-max vhat): targets 3,6,4: tail D_j at eps_hi=1e-12 evaluated at eps=1e-12,1e-13 and chain D_j at 1e-9:
  every E_j of my tail/chain >= independent E_j, min ratio 1.0012-1.0016 (the 0.1% safety), max 2.75 (small E_j).
