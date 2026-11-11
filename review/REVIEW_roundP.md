# REVIEW_roundP: referee report on the pre-registered optimizer (P1) and depth (P2) tests

Archived verbatim from the referee's hand-back (2026-11-10). The referee left the repository untouched; scratch code is outside it.

REFEREE REPORT, P1 (optimizer) and P2 (depth) tests. No major or blocker problem in the numbers: every reported value reproduces. Two paper claims go beyond the evidence (findings 1 and 2), and the P1 gap depends on the activity threshold and learning rate in a way the text does not say (finding 3). Scratch: referee_P_work/ (an1..an9.py, rerun.py, tau.py, flow3.py, scal.py, rerun.jsonl, tau.jsonl). Repo untouched (git status clean, bytecode off).

FINDINGS
1. MAJOR (abstract, overclaim). "A merged penalised control is level with gradient descent only when its penalty is chosen by the truth" is true only at the declared radius 0.3.
   - My paired contrast from the raw records, with the validation one-SE control minus GD (control / GD):
     - 0.3: 0.654 / 0.726 (GD ahead by 0.072 [0.004, 0.140]).
     - 0.4: 0.760 / 0.758 (level).
     - 0.5: 0.840 / 0.786 (control ahead by 0.053 [-0.001, 0.110]).
     - 0.6: 0.889 / 0.773 (control ahead by 0.116 [0.068, 0.167]).
   - Supp S4.7 admits "ahead of GD at 0.4 and 0.5"; the abstract and sec_exp do not.
   - Fix: add "at the declared radius 0.3; at radii 0.4 and above the validation-chosen control is level or ahead (post hoc)". Do not claim that GD selects better than the control. P1 supports the optimizer separation only.
2. MAJOR (P2 wording).
   - (i) "The depth exponents p/2 are confirmed on the scalar unit" overclaims. By construction the margin is fixed relative to lam, so the exponent is p/2 identically. I re-integrated in eta with DOP853, rtol 1e-13: exponents 1.5000000, 2.0000000, 2.5000000, relative deviation from the closed form at most 3e-13. This checks the integrator and the repository's `scale(eps)`, not the theory. The supp itself says "near-identity".
   - (ii) The abstract omits that the pre-registered primary arm (unequal masses) was INCONCLUSIVE; only the secondary equal-mass arm held. The "exactly 1.000 once the margin factor is divided out" statements, and the "identity check" column, are circular. The margin is inferred from w with the same closed form. They show nothing about the network.
   - Fix: say "confirmed as a check of the closed form and the numerics (scalar); on the flow, equal-mass arm held, unequal-mass primary inconclusive".
3. MAJOR (robustness, text silent). The Adam and sign gap at radius 0.3 is robust to the learning rate and the weak tuning, but it closes at larger radius and activity threshold.
   - Own rerun: 10 test runs (targets 411, 415, 420, 425, 430; noise 0.1 and 0.4; seed 42), K>0 exact share at radius 0.3, selected on the test truth:

     | arm | shared-v | GD | Adam, one lr (8 lrs 1e-4..0.3) | Adam, GD group rates (7 lrs) | sign, one lr (7 lrs) | sign, group rates (5 lrs) |
     |---|---|---|---|---|---|---|
     | share | 0.842 | 0.783 | best single lr 0.383 | per-run oracle over all 15 Adam configs 0.475 | best single lr 0.358 | per-run oracle over 12 sign configs 0.508 |

     (The 0.383 is the best single lr, reached by the group-rate variant at lr 1.0; the best uniform lr gives 0.233.)
     - So an Adam or sign lr tuned on the test truth does not reach GD at radius 0.3.
     - The registered Adam arm used one global lr with no group rates, while GD carries the spike's lr 1 and group rates (nu 0.05, bias 0.1). The group-rate variant I added does not rescue Adam at 0.3. The dev grids (Adam lr 0.003–0.1, sign 1e-4–1e-2) were reasonable; Adam's lower bound 0.003 versus my 1e-4 changes nothing.
     - The Adam lr curve is non-monotone and noisy (12 dev runs).
   - At larger radius a tuned Adam does catch up. In the same 10 runs: radius 0.5, best single lr 0.633 against GD 0.683 (frozen lr 0.03: 0.458); radius 0.6, 0.708 against 0.708; per-run oracle 0.858 and 0.883 against GD 0.683 and 0.708.
   - Activity-threshold sweep on the same 10 runs, refitting GD and the frozen Adam, sign and shared-v arms and varying tau (radius 0.5):

     | tau | GD | Adam | sign |
     |---|---|---|---|
     | 0.02 | 0.683 | 0.458 | 0.425 |
     | 0.05 | 0.733 | 0.658 | 0.608 |
     | 0.1 | 0.733 | 0.767 | 0.808 |

     At radius 0.3 and tau 0.1, sign is 0.608 against GD 0.833, so the gap there is 0.225, still above 0.20. Adam is 0.492, so its gap is 0.34. At tau 0.1 and radius 0.5 Adam and sign reach or pass GD.
   - Only 10 runs on one seed, so noisy. The point is qualitative: "Adam is not rescued by merging" holds only at tau 0.02, the registered lr, and radius up to 0.99. With the frozen lr it holds for all radii on the full 120 runs: the (a) gap is at least 0.297 at 0.5 and at least 0.30 from 0.6 to 0.99, with lower ends of about 0.19 and above.
   - Fix: qualify the supp "(a) holds everywhere" with "at tau 0.02 and the learning rates frozen at radius 0.3". Add that Adam's duplicated units are removable by a larger merge radius together with a larger activity threshold. That is what theory O says: Adam loses the ordering but not the units.
4. MINOR (P2.3 not fully held-out).
   - The D2 and D3 columns (D2 0/200 at 1e-6; D3 176/200 and 195/200) are bit-identical replications of results/chaos/parametrization.md lines 115-118 on the same 200 targets and code path. Only D4, D5 and the extra eps are new.
   - The plan's "the author did not read the fractions" is disclosed, but the paper headings say "held-out".
   - Fix: write "re-run on the readouts study's 200 targets, replicating its D2 and D3 numbers" for those columns.
5. MINOR (D5 integrations).
   - The 4 failures are real: DOP853 "Required step size is less than spacing between numbers", targets 120, 122, 123, 124, D5 at 1e-8. At the same cells D2, D3 and D4 succeeded. At 1e-10, 6 of 7 unequal-mass and 3 of 7 equal-mass cells also failed, so the p=5 exponent is fitted on 3 eps values (1e-4, 1e-6, 1e-8) for most targets. The paper should say so.
   - The verdict arithmetic is right: f5(1e-8) = 196/200 = 0.980 < f4 = 1.000 and < f5(1e-6) = 1.000. FAILED is the correct label because the plan says a failed integration is a failure.
   - INCONCLUSIVE would not be rule-compliant; the plan has no such state for P2.3. The registered design was fragile: with dev saturated and tied at 1.000, one integrator failure flips the verdict.
   - Post hoc, among completed integrations, 196/196 at 1e-8, so the ordering holds with ties.
6. MINOR (P2.2).
   - The sensitivity (no margin exclusion) at p=5 is 2.448, which is 0.052 from 2.5, i.e. outside the 0.05 band the registered rule uses. p=3 is 1.470 and p=4 is 2.000.
   - Excluded targets had exponents 1.234 (s149, p=3), 2.124 (design, p=5) and 2.556 (s134, p=5).
   - Unequal-mass p=3 identity-check exponents are 1.532, 1.557 and 1.543 for three targets, because the maximising unit changes with eps (b_j varies).
   - The equal-mass exponents are within 0.001 of p/2 mainly because the realised margin is eps-independent there, which is the "observed law".
7. MINOR (referee-number overlap).
   - The referee's numbers (GD 0.72, Adam 0.25, shared-v 0.81, paired gap +0.48) were known, and the test reproduces them (unmerged K>0: referee 0.62, 0.012, 0.75; test 0.627, 0.026, 0.664). P1 (a) and (b) are a replication on fresh targets, not a discovery. The plan and supp already say this. Sec_exp "registered before any held-out run" is fine.

VERIFIED (how)
- Timestamps, git log %ci: plan 4b9754e at 13:55:02 UTC; development, freeze and amendment 1 at 5e03c0d, 14:33:54; amendment 2 at ad40ecb, 14:36:58; raw and report at 18:49:21. Between 5e03c0d and 3153764 the script changed only by a post hoc table of successful integrations and a default path for the freeze file. The post hoc table is labelled POST HOC.
- Summed fit seconds in the stored records (held-out about 11,500 s) fit inside the 14:34–18:49 window (15,300 s) with no sign of an earlier held-out run. The development sum (2,370 s) is about 40 s over the 2,332 s window between plan and freeze, which suggests a short overlap or parallelism; immaterial.
- Freshness: `exp_optdepth.py --check-unused` passes (ids 401-403 and 411-430, generator seeds, 230 data-draw seeds, no identical bumps). Only `exp_kan_spike` and `exp_optdepth` use these ids. The plan's used-id list includes the referee's 101-103 and 301-306 with seeds 1-3, the spike 1-10, and theory3f.
- P1 numbers recomputed with my own merge and cluster bootstrap (20,000 resamples). All table values, all four gaps, (b), (d), the radius table and the exact-on-all-5-edges counts agree. GD 35, adam_big 43, shared-v 54, Adam 0, sign 0; control 23 (validation) and 33 (oracle); oracle share 0.764.
  - My intervals differ from the report's by 0.001 at most. Minimum true bump spacing is 1.007.
  - Adam's share by radius: 0.49 at 0.5, 0.47 at 0.6, 0.40–0.43 from 0.7 to 0.99. It never reaches GD (0.707–0.786).
  - The report regenerates identically from the raw (diff empty).
- Scalar P2.1: independently integrated, as above. The repository's lam matches s^-(p-2)/(ell p (p-2)).
- P2.2 exponents recomputed from raw: all 12 means match (included, all-target and equal-mass), for example unequal masses p=3 1.510, p=4 1.983, p=5 2.492, and equal masses 1.502, 2.001, 2.501. I re-ran the repository flow for the design target and sample target s143 (equal masses, my own slope code): design p=3 1.486, p=4 2.001, p=5 2.501 (3 eps values; eps = 1e-10 failed); s143 1.501, 2.001, 2.500.
- P2.3: counts recomputed from the raw. D2 0/200 at all four eps; D3 7, 89, 176, 195; D4 148, 182, 197, 200; D5 178, 193, 200, 196.

COULD NOT CHECK
- That the held-out run really started after 14:34. Commit times are written by the author, and there is no external clock or raw-file mtime.
- Runs of the full 53-test suite or the mutation record.
- A reimplemented flow beyond two targets.
- The tau and radius sweeps beyond 10 runs; there are no per-unit rms values in the stored records.
- Adam variants beyond lr, group rates and eps (betas, warmup).

WORDING ALLOWED
- "Pre-registered, 20 fresh targets, 120 runs per arm: at the declared radius 0.3 and activity threshold 0.02, GD (0.726) and shared-second-moment Adam (0.788) exceed default Adam (0.240) and sign descent (0.184) by 0.49–0.60 in K>0 exact-edge share, and are within 0.10 of each other. This replicates a prior referee rerun."
- "The gap shrinks as the merge radius grows (0.30 at 0.5) because merging repairs Adam's duplicated units; an Adam or sign descent learning rate tuned on the test truth does not close it at radius 0.3 (10-run check, post hoc)."
- "P2.1 holds as a check of the closed form and the numerics. P2.2 is INCONCLUSIVE for the primary unequal-mass arm (p=4: 2 of 7 targets met the margin rule) and HELD for the secondary equal-mass arm."
- "P2.3 FAILED by the registered rule, solely through 4 D5 integration failures at eps 1e-8 (a failed integration counts as a failure). Post hoc, among completed integrations the ordering holds, with D4 and D5 tied at 1.000."

MUST NOT CLAIM
- GD or shared-v select better than the merged control. Control minus GD is +0.116 at radius 0.6 and level at 0.4.
- Adam or sign is not rescued by merging, without "at the registered lr and tau 0.02".
- Exponents p/2 "confirmed" on the flow without the equal-mass qualifier, or the scalar result as evidence for the theory.
- P2.3 as evidence against the depth ordering, or "held-out" for its D2 and D3 columns.
- The D5 exponent as a 4-point fit.
- Statements for the fresh data that the referee's numbers were unknown.
