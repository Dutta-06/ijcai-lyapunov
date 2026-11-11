# REVIEW_roundT3N: referee report on the corrected, merged verification of Theorem T1 (workstream N)

Archived verbatim from the referee's hand-back (2026-11-10). The referee left the repository untouched; scratch code and logs are outside it.

REFEREE REPORT, round T3N (merged corrected verification of Theorem T1). Repository untouched. Scratch code and logs are in referee_T3N_work/.

VERDICT: no blocker. All six findings of REVIEW_roundT3AC.md are repaired in every path actually used for the claim, and everything I re-ran reproduces. Remaining issues are one hand-proof/documentation gap (finding 1), stale legacy reports (finding 2) and wording (findings 3-5).

FINDINGS
1. MAJOR (hand proof / documentation, not the code). Every verifier computes a window-GATED displacement term: hinc counts unit i's full totals only after wstart_i = t_i - S_i^-, and uses psi_sup*pre before that (repair (g) is exactly the pre part).
   - THEORY_T1.md Lemma 5.2 is ungated ("totals from 5.4"), and grep finds no statement of the gating anywhere in it.
   - Gating is valid. I re-derived it: before t_a, int v_i <= pre (Lemma 3.3(1)); before t_i the limit coefficient is 0, so the Gamma'' term vanishes; t_a >= t_i - S^- by Lemma 5.3; the entrant's centre stays within r of its feature (Lemma 3.1a). So indep.py's gated='pre' is the correct reading.
   - The literal ungated reading would not close at the tops. My evaluator at D = 1.5 x float E gives:
     * target 6, eps = 10^-8.5: E 131 (gated) vs 254 (ungated), E/D 0.67 vs 1.87;
     * target 4, eps = 10^-8.0: E 108 vs 456, E/D 0.67 vs 3.30.
   - Cheapest fix: add the three-line gated statement and proof to Lemma 5.2. Until then "S implies T1" is a hand proof that does not contain the step the code uses.
2. MINOR. theory3c.md is still headed "Fully rigorous interval verification ... 8 verified targets", and theory3a.md, theory_t1.md and theory_t1b.md carry no "superseded by theory3n.md" banner. Those files use the legacy accounting. Add a banner, or the paper may cite them.
3. MINOR. theory3n.md says the Gram-tail fix of theory3n_kernels.py is "built in". No verifier imports that module (grep: only tests and docstrings), and no number in the report depends on it. Reword to "provided and tested; not used by any result here".
4. MINOR. Repair (6) is present in all four paths but is not verdict-material. The float dcmax is identical to 5 digits with and without it. In the tail, E/D is 0.99927369 vs 0.99927367. No mutation of (6) flips a verdict; only the value tests can see it.
5. MINOR (attribution of the two float flips).
   - (2,8,1,13) at 1e-7.5 is caused by (f) alone: dcmax goes from 0.2441 to 0.3472 at the best parameter set, against d_c = 0.3.
   - (6,8,1.5,3) at 1e-8.5 needs both (f) and (g). Reverting (f) alone leaves 0.3015, reverting (g) alone leaves 0.3228, and reverting both gives 0.2991, which passes.
   - So (g) is not cosmetic in the float grid.
6. NOTE. At eps <= 1e-12 the tail verdict is insensitive to (f) and (g). dcmax <= 2e-6, and the legacy and corrected E agree to 1e-13 relative. The tail ball at 1e-12 therefore has large slack and does not depend on the repairs. They matter near the tops.

(1) TRACE AND MUTATIONS. I traced the code and mutated it myself.
- Float cert: cert_n and cert2d_n are compiled from the anchored source (each anchor occurs exactly once; build_cert(legacy_*) reproduces the old function). Compiled code has Wprev, hpre-style terms and _stage_max_vhat, and the globals are bound correctly (rgs2_cached and G2_sup2 for the planar one).
- Chain: verify_rig_n has Wprev, hpre, _stage_max_vhat_rig, _min_lower and _inv_of. Gate arithmetic checked by hand (window st can start inside stage st; windows i<st have started; windows i>st start after the stage by the 'order' condition).
- Tail: no flags. (f), (g) and (6) are hard-coded: Wq stage-resolved, hpre where gate >= the sub-interval end, vhat_fin stage maximum. I wrote mutants tail_mg.py, tail_mf.py, tail_m6.py and tail_mall.py.
- Mutation results:
  * Chain, stored failing balls: target 6 [-8.75,-8.5] and target 4 [-8.1875,-8.0625] are rejected at every factor. Reverting (f) accepts both at factor 1.5 (E/D 0.954 and 0.938). Reverting (g) accepts target 4 at 1.5 (0.987). Reverting (6) changes nothing. Dcmax moves from 0.0925 to 0.0786 and from 0.1716 to 0.1456.
  * Tail, target 3 at eps_hi = 10^-7.5, same D scaled to 0.99, 0.98 and 0.96: the corrected tail rejects at every one. The (g)-mutant rejects at 0.99 and below, the (f)-mutant accepts down to 0.96, the all-legacy mutant down to 0.96. At unscaled D the E/D values are corrected 0.9993, (g) reverted 0.9939, (f) reverted 0.9712, all reverted 0.9658.
  * Float, 2 flips: corrected certificate false at all 16 parameter sets; legacy true (see finding 5).
  * Planar (2,8,1,1) at 1e-6 and 1e-7: dcmax and Dnew move with (f) and (g), so the planar path has them.

(2) INDEPENDENT RE-EVALUATION. I ran the earlier referee's indep.py (gated='pre', stage-maximum vhat, stage-resolved displacement), whose gating I checked in finding 1.
- Targets 8 (2,6,2.4,12) and 6 (3,8,1,31). Evaluation points:
  * 1e-14: inside the tail ball;
  * 1e-12: the join, evaluated with the stored tail D;
  * 1e-9 or 1e-10: chain, with my own fixed-point D, compared with verify_rig_n on a narrow ball;
  * the top (target 6 at 10^-8.5, target 8 at 10^-7) on the same footing.
- All conditions hold in my evaluator. My E_j is below the repo's majorant for every kernel at every point (ratios 0.9908-0.9997 in the tail, 0.74-0.9996 in the chain).
- At 10^-8.5, target 6, the repo majorant for one kernel is 21 % above mine (0.826). That is conservative, not a violation.
- I did not run a K=3 evaluation of target 6 at the top with the repo's stored D (it was not needed for the comparison).

(3) JOIN AND TAIL INPUTS.
- Chain contiguity from the stored raw (E_segs), all 8 targets: exponents are contiguous, the top ball upper end >= 10^top, the bottom ball lower end is fdown(10^-12) = eps_meet = 1e-12 exactly, consecutive balls overlap as doubles.
- Tail: h_hi = upper(sqrt(eps_meet/(2vs - eps_meet)))(1+2^-40). eps_of_h(h_hi) >= eps_meet is tested with lo(...) >= arb(eps_meet) and is True for 8 of 8. The union of the chain and the tail ball is therefore contiguous on (0, eps_top].
- Tail float usage: float()/np appear only in prep['K'] and prep['As'] (target data), the quadrature grids (outward partitions), and info diagnostics. The verdict is arb: Etot*(10001/10000) < A(Dv[j]), strict and certain. There are no .mid() or float(ball) calls in any decision. Limit data come from LR.
- Re-run on my side: the tail re-verified for all 8 targets from the stored D, with all conditions and the join true. Chains for targets 0, 4, 5, 6 regenerated with the same ball and call counts as stored (7/25, 10/54, 8/34, 6/28).
- Tests: tests/theory3n_test.py gives 32 passed.

(4) RECOUNT from the stored raw (my own script):
- informative pairs, old to corrected: K<=3 228 to 227; K=4..6 204 to 203; planar 41 to 41; total 473 to 471 of 771 scanned.
- true to false: 2. The pairs are C (2,8,1,13) at 1e-7.5 and K (6,8,1.5,3) at 1e-8.5. False to true: 0. Violations: 0.
- eps_0: 1e-7.5 to 1e-8.0 and 1e-8.5 to 1e-9.0.
- The "old" counts equal the earlier raw files' counts.

(5) OVERCLAIM CHECK of the md. The bold line "S(eps) is verified ... every real eps in (0, 10^TOP]" is acceptable only together with the Status block. That block correctly lists:
- Theorem K(a) global minimality as not verified;
- the hand proof S => T1 and the horizon;
- K>=4 and planar as floating point only;
- ingredient tables not regenerated.
Add to the paper: no monotonicity claim; continuum Gram only; each ball has its own (pi, D), which is fine because S(eps; pi) is per eps. Findings 1-3 are the remaining wording issues.

COULD NOT CHECK
- Hand proofs of sections 2-6 beyond Lemma 5.2, Lemmas E and M (from the previous round), and the gating step in finding 1.
- Full regeneration of the chains for targets 1, 2, 3, 8 (about 74-110 balls); I checked their stored contiguity and ran spot balls instead.
- Whether the float certificate at the other 769 grid pairs is correct. I checked the two flips, the counts from the raw, and the mutation sensitivity, not every pair.
- Soundness of arb; K>=4 and planar rigour (floating point only).
- Mutant records of the repo (own spec reads sensibly; not re-run).

SENTENCE FOR THEOREM T1's STATUS (use verbatim):
"For eight test targets (n = 1, K <= 3, continuum Gram, equal masses, noiseless data, law (12, 1, 0.1, 1)) the finite system S(eps) of Theorem S-T1, with the corrected displacement accounting (stage-resolved, window-gated with the pre-layer term), has been verified in ball arithmetic with Krawczyk-enclosed limit-object data and certain comparisons, for every real eps in (0, eps_top], eps_top = 10^-7.5, 10^-5.5, 10^-6.5, 10^-7.5, 10^-8.0, 10^-7.0, 10^-8.5, 10^-7.0 (a chain of eps-balls on [10^-12, eps_top] and one tail ball in h = 1/kappa on (0, 10^-12]). Combined with the hand proof that S(eps) implies the conclusions of Theorem T1 on the horizon [0, lambda(1.05 tau_K + 1)], this gives Theorem T1 on these targets in a computer-assisted form. Not verified: that the enclosed block critical points are the global constrained minimisers of Theorem K(a), the hand proof by machine, any monotonicity of S in eps beyond these covers, sampled Gram, noise."
(Use only after finding 1, the gated Lemma 5.2, is written into the hand proof.)

SENTENCE FOR THE FLOAT GRID (use verbatim):
"With the corrected accounting, the floating-point certificate holds on 471 of the 771 grid pairs (227 of 345 for K <= 3, 203 of 345 for K = 4..6, 41 of 81 in the plane); on all 471 the integrated network forms exactly the limit's units (0 violations). Two pairs of the earlier 473 flipped, and eps_0 of (2,8,1,13) and (6,8,1.5,3) moves by half a decade. This is floating-point arithmetic without enclosures (K >= 4 and the plane carry no interval verification), and the earlier ingredient tables were not regenerated."
