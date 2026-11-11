# Referee report PZ1: the generality study (WP1)

## Scope

Reviewed state: commit **5083a49** on `claude/lyapunov-context-setup-yg51ew`, all ten parts A-J. It supersedes 7d45e0c (named in the brief) and 6c87a4e (Part F rerun at h = 0.02, Part J added). Part C, the rerun of Part I and the Radau fallback arrived in 5083a49. I reviewed the three states in turn and re-checked every change. All verdicts are on 5083a49 and nothing is pending on Part C any more.

Method:
- `git archive` exports of the three commits under W = /tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/refPZ1 (repo = 7d45e0c, repo2 = 6c87a4e, repo3 = 5083a49).
- The repository was never modified. Its venv was used read-only (Python 3.13.18, numpy 2.5.4, scipy 1.18.1). Every run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. Nothing was written outside W.
- One heavy process at a time, except three overlaps with a second process (a convergence job of about 10 min, a coherence job of about 5 min, and two scripts under 10 s). They were stopped when noticed.
- Read in full: the files named in the brief at all three commits, with diffs; the paper passages named; `rich_transit_nd.sliding_iss_nd/BumpsND/path_duals_nd`; `rich_transit.moving_rhs/critical_rhs/transit_scale/relax_chunk`; `mirror_free.FreeLaw`. CLAUDE-style summary file CLAIMS.md was read last, and again after it was updated for Part C.

Not done: Part C was not rerun in full (96 cells = 16 targets x 6, and 84 cells = the 14 targets with horizon > 1000 or amplitude ratio > 100). Of Part J I reran 4 of 36 cells. Literature: one full text (Woodworth et al.) and three abstracts.

## Summary

### Design of the study
One Gaussian-kernel network (D'''.48's dictionary: 21 kernels of width 0.12 at spacing 0.1, 61 samples, centre rate Gm = 12) is trained from s^2 = eps with factorised readouts: D1 w = p, D2 u^2 - v^2, D3 u^3 - v^3, P0 the paper's unit, P0u without saturation, P4 with equal rates. Training is by gradient flow (LSODA, now with a Radau fallback), fixed-step descent and Adam. Each network is compared with the sliding inverse scale space path (`sliding_iss_nd`, polished) through:
- the formed set (|w| > min(0.01 max|A*|, 0.1 min|A*|)), the entry order and entry times;
- the gap between formed and dormant coefficients, and the dormant displacement;
- an "exact" test (K units, 0.05 sigma, 5%).

Ten parts, 1,702 cells, all with status ok. Strengths: failures are tabulated (max columns, failing-cell list); the md regenerates byte for byte from the raw; cells recompute bit for bit; the unstable descent run of 7d45e0c was found and corrected by the authors.

### Verdicts on F1-F4

**F1: CONFIRMED.**
- Derivation: with w = nu^2 a(alpha)/m, dF/dalpha = -g beta a'/m, dF/dnu = -2 nu g a/m, dF/dc = -w<r grad_c phi>. The paper's eqs (free), (centres) (paper17.tex l.157-165) then give exactly the rates (Gm, gamma m/4, Gm), and v_s/a_max^2 = gamma/(4G) = lr_nu/lr_alpha [own].
- Finite differences at three constant sets with G != 1 and a_max != 1, (m,G,gamma,a_max) = (12,1,0.1,1), (5,2.3,0.37,0.7), (0.8,0.45,3.1,1.9): paper law = -lr * grad F to 1.5e-8, 1.8e-9, 1.0e-9. Wrong-rate controls (lr_nu = gamma m/2, /8, gamma m, 4 gamma m; lr_c = Gm/2) fail by 0.5 to 15 (relative) [own].
- Repository `PaperLaw` vs `moving_rhs`: 5.5e-16, 4.2e-16, 4.6e-16 at the same sets [own]. The md's 8.21e-16 is the same check at the paper's constants [md].
- Tests are not vacuous for rates: three rate mutants are killed by the identity test (rtol 1e-12). They never vary G or a_max, so four mutants that drop G or add a spurious a_max^2 or G survive; the G/a_max dependence is verified by me, not by the suite.
- Scope: "exactly Euclidean" holds for the paradigm centre law, not the WFR variant (the study says so).

**F2: PARTLY CONFIRMED; overstated as worded.**
- Design: confirmed and reproduced independently (see section 1-6 below).
- 200 targets [md C, 1,200 cells ok]:
  - K units on the K features: P0 187/200 (eps 1e-6) and 197 (1e-8); D3 192 and 198.
  - The path's own kernels: P0 166 and 191; D3 176 and 195.
  - Median gap 5.4/8.4 decades (P0), 4.5/7.3 (D3); minimum 0.74 (P0) and 0.79 (D3) at 1e-6.
  - Dormant displacement: median 0.036/0.0006 sigma (P0), 0.019/0.0003 (D3); maximum 2.57/0.78 (P0) and 1.2/1.45 sigma (D3).
- So "exactly the path's units ... at 1e-6 and below", "dormant centres staying put" and "several decades" hold for the median target and for a stated class of targets, not for every target. The class is a margin condition [own]; see PZ1-1.

**F3: CONFIRMED on the benchmark's generic (off-grid, coherent) targets; the constants and "no practical epsilon" need their conditions.**
- Reproduced with my own code on 21 targets and on the 200 targets [md C].
- No D2 variant I tried passes at eps in [1e-12, 1e-6] on this dictionary. The rule behind it is an inequality on the runner-up margin (PZ1-2).

**F4: CONFIRMED.**
- Depth-1 control: 12 units at every eps on the design [md B].
- Fixed-step descent reproduces the flow once stable. The authors' stability estimate is right (Part F, below).
- Adam as reported. One CLAIMS statement about Adam at 1e-6 is wrong (PZ1-10).

### Numbered findings (detail in the Issues table)
- PZ1-1 (major): F2's set identity, order and "static dormant units" are conditional on the runner-up margin. 98/98 hold above margin 0.02; the failures are all below.
- PZ1-2 (minor): F3's constant 4.74 is ln(1/theta); "no practical eps" is a margin inequality, ln(1/theta)/ln(1/eps) against the margin.
- PZ1-3 (minor): Radau fallback bookkeeping. The claim "A, B, D-H, J cannot trigger it" is false for E and J; the raw holds two versions of one computation; usage is not recorded.
- PZ1-4 (minor): the wall-clock limit was lowered after the runs. The committed code does not reproduce J[8] (and F[2], E[3]) as ok on this machine.
- PZ1-5 (minor): non-ok cells are not counted as successes, but they drop out of denominators silently in D, E, G, I, J.
- PZ1-6 (minor): Part C presentation: a truncated failing-cell list, and a "max tau err." column that is a record-grid artefact.
- PZ1-7 (minor): the test suite. 23 of 77 mutants survive (G and a_max never varied; tau_rel, order and Adam semantics unpinned).
- PZ1-8 (minor): "depth >= 3" is tested only at depth 3 in the study. My D4 and D5 pass.
- PZ1-9 (minor): Part G's amplitude rows do not test the paper's logarithmic regime. My branch-pair test does.
- PZ1-10 (cosmetic): md and CLAIMS wording.
- PZ1-11 (minor, paper text): l.127 "not a Euclidean one".

## What I checked and how

**1. F1.** Derived independently of the code; finite-difference, wrong-rate controls and repository-vs-`moving_rhs` at three constant sets (see F1 above). Mutants of the learning rates, saturation and centre law (section 9).

**2. Trying to make depth 2 pass, and depth 3 fail** [own, design target unless stated; own code `myfr.py` with LSODA, Radau, DOP853 and reduced mirror coordinates].
- D2 baseline at eps 1e-6: 9 formed. Counts at theta = 0.01/0.03/0.1/0.3/0.5 of max|A*|: 9/8/6/3/2.
- Threshold-free gap (K-th over (K+1)-th largest |w|): D2 0.17 at 1e-6, 0.37 at 1e-8, 0.68 at 1e-12, 1.55 at 1e-24, 3.2 at 1e-48; D3 5.3 at 1e-6. D2's widest gap sits after 12 units, D3's after 3.
- D2 variants at eps 1e-6 and 1e-12:
  - centre rate x0.01, x0.1, x10, x100; lr 0.1 and 10 (only the ratio lr_c/lr matters, as it should: lr 10 equals centre x0.1);
  - horizon x3 and x10; asymmetric starts u0^2/v0^2 = 10^+-6; rtol 1e-12.
  - None passes. At theta = 0.01, at least 6 units form; the three largest amplitudes are off by 31-160% (1e-6) and 14-54% (1e-12).
  - At theta = 0.3 the set is right at 1e-6 and 1e-8, but the amplitudes are off by 39% and 29%, and the rms is 2.8e-3 against D3's 1.3e-7.
- D2 on Part D's 21 targets, threshold sweep: the path's set at theta = 0.01 / 0.1 / 0.3 / 0.5 is
  - eps 1e-6: 0 / 0 / 7 / 1
  - eps 1e-12: 0 / 3 / 12 / 5
  - eps 1e-24: 4 / 7 / 17 / 9
  - eps 1e-48: 9 / 11 / 17 / 12
  - eps 1e-96: 12 / 17 / 18 / 15
  - eps 1e-144: 17 / 18 / 18 / 15
  - "Top-K exact" (K largest units within 0.05 sigma and 5%) is 0, 1, 10, 13, 20, 21 of 21. Median top-K amplitude error is 0.49, 0.34, 0.09, 0.001, 0, 0.
  - D3 on the same 21 targets: exact 11 (1e-4), 21 (1e-6), 21 (1e-8); gap 2.3, 4.4, 7.3 decades.
- Margin rule: D2 selects iff every not-yet-active kernel stays at margin mu > ln(1/theta)/ln(1/eps). The threshold is 0.341 at 1e-6, 0.170 at 1e-12, 0.085 at 1e-24, 0.043 at 1e-48, 0.021 at 1e-96, 0.014 at 1e-144 (theta = 0.009). With no fitted parameter, it predicts the outcome on 21/21, 21/21, 20/21, 18/21, 20/21, 20/21 of Part D's targets. Benchmark margins: median 0.019, maximum 0.149.
- Dictionary spacing (off-grid design features, repository engine, own reference path per dictionary): D2 forms 9, 8, 6, 5, 4, 4 units at eps 1e-6 for 21, 16, 13, 11, 9, 7 kernels (spacing 0.1 to 0.33); at 1e-12 it passes only at spacing 0.33. D3 and P0 pass everywhere at 1e-6 (gap 5-8 decades). With features exactly on kernel centres and neighbour overlap <= 0.5 (spacing >= 0.2) D2 does pass at 1e-6 (own `myfr`, gap 2.75-6.2 decades). That is the non-generic case, and it shows the verdict is about margins, not about depth alone.
- Entry-time "error": D2's formation times are all early by the same factor. Signed errors -0.3405, -0.3399, -0.3391 at 1e-6 and -0.1704, -0.1702, -0.1689 at 1e-12. Entry-time ratios therefore match the path to 0.1%. 4.74 = ln(1/0.009) + 0.03 (4.71). At theta = 0.5 max|A*| the error is 0.052 at 1e-12 instead of 0.170.
- Equal initial coefficient growth rather than equal eps (D3 coefficient grows as eps^2 * eta, D2's as eps * eta): D3 at eps 1e-3 forms 6 units (gap 0.94); at 1e-4 it passes. D2 at 1e-8 forms 8.
- Trying to make depth 3 fail (design):
  - centre rate x0.01 and x100, Radau rtol 1e-11, LSODA rtol 1e-12, eps 1e-6 and 1e-8, D3 and P0: all 3 units, set, order, exact.
  - D4 (u^4 - v^4) and D5: pass at 1e-4 and 1e-6 (D4 also 1e-8), gaps D4 5.4/9.4/13.4, D5 7.8/12.8 decades. D5 at 1e-8 crashed the 6c87a4e engine (see 10) and finishes under 5083a49.
  - K = 2, tie and near-tie targets (samples 106, 110, 117, 129, 143, 120): D3 and P0 at 1e-6 give set, order and exact.
  - Genuine failures of "same set as the path" found on benchmark targets: PZ1-1.

**3. Measures.** The formed threshold is the paper's convention (paper17.tex l.526) and is not biased against D2. Raising it to 0.3 helps D2's set but not its amplitudes. The study's gap is measured at the threshold and is more generous to D2 than the threshold-free K-th/(K+1)-th gap (0.76 vs 0.17 at 1e-6). Dormant displacement is measured on unformed units only, which also favours D2, and it is still 0.3-0.8 sigma. The tolerances (0.05 sigma, 5%) decide the exact flag in 10 of 1,054 cells with K units formed: 6 fixed-centre rows, 2 P0 1e-4 descent rows, and 2 D2 cells (targets 134, 143 at 1e-12, errors 0.127 sigma/12% and 0.042 sigma/5.1%). The header's "lowered where a feature is 10^3 times smaller" is wrong: the second term of the min binds from amplitude ratio 10 (20 of the 200 targets have ratio > 10).

**4. Integration.**
- Independent integrators (own RHS written from the paper's formulas): LSODA, Radau, DOP853 and reduced mirror coordinates (D2, D3) agree on D3 1e-4 and 1e-6, D2 1e-6 and 1e-24, and P0 1e-6 (LSODA and Radau). Final W differs by <= 4e-10 (D3 1e-4), 1e-12 (D3 1e-6), 1.3e-9 (D2), 1.4e-13 (P0); every printed measure is identical.
- Fixed-step descent with halved steps (D3, eps 1e-6, benchmark target 1): h = 0.04, 0.02, 0.01 give |dW|/max|W| of 4.1e-6, 2.0e-6, 1.0e-6 and |dC| of 5.2e-4, 2.6e-4, 1.3e-4 sigma against the flow (first order in h). h = 0.05 reproduces the old failure (rms 8.9e-2, |dC| 0.40 sigma).
- Flow rms 6.99e-6 at 1x horizon, 6.82e-6 at 3x, 6.28e-6 at 10x.
- P0 at eps 1e-4 on targets 1-4, h = 0.02: |dW| <= 1.8e-5, |dC| <= 3.5e-5 sigma, formation times <= 1e-5, even where neither flow nor descent is clean.
- Statuses: the code never counts a non-ok cell as success (PZ1-5). The committed raw has no non-ok cell.

**5. Sample and reference path.**
- The sample is stratified, not uniform: A 70 of 400 (nearest/other 35/35), B 30 of 160, D(a) 20 of 54, D(b) 50 of 240, D(c) 30 of 180. The "other-kernels" class is over-weighted (82/200 = 41% vs 354/1034 = 34%); amplitude ratio > 100 on 11/200 vs 24/1034. This makes it harder than the population.
- The reference path reproduces the stored entry order on 200/200 (own run, with relaxed = True on every block, no solver failure, no exits, final rms <= 2.3e-14, one polish rejected) and the stored class (nearest = E on 118, other on 82). The J reference reproduces the stored order on 12/12.
- K = 2 targets are handled correctly: `analyse` uses len(A*); P0 1e-6 gives K units on 65/70 K = 2 targets and 122/130 K = 3 targets.
- Hazards: PZ1-6. The family-A-only subsets (Parts D, E, G, I) are the easy family.

**6. Reproducibility.**
- The md regenerated from the raw with `--from-raw`, twice, is byte-identical to the committed md (md5 34de2a94fb41a32530d8058d6702f393) at 5083a49 (also at 7d45e0c and 6c87a4e).
- Tests: 45 pass at 5083a49 (33 + 12; 42 at 7d45e0c).
- Cells recomputed bit for bit (every stored field except wall-clock secs): 42 cells of B, D, E, G, H, I, F under both the 6c87a4e and the 5083a49 engine (16 groups, all identical); 96/96 Part C cells (targets 0, 1, 8, 36, 42, 59, 70, 75, 100, 108, 119, 124, 140, 150, 171, 190); 84/84 cells of the 14 Radau-prone C targets; J[4], J[2], J[11] P0.
- Raw parts A, B, D-H, J are identical between 6c87a4e and 5083a49 apart from secs; in I 130 of 132 cells equal the previous raw.

**7. CLAIMS.md against the md.** Section 5 numbers are verbatim-correct, including the new Part C row (187/192, 197/198 exact; 166/176, 191/195 set; gaps; displacements; D2 0/200 and 2/200), Part F (h*lambda_max 0.44-0.87; 43 -> h < 0.046), Part I (11/11), Part J. Exceptions:
- "Adam: at 1e-6 nothing has formed" is true only for P0. D3 and D2 have 21 formed (rms 0.39, 0.28; md H shows "21 (21/21/21/0)").
- The claim map lists test_factorized_rbf.py (31); it now has 33 tests.
- "What is new ... the quantification (the 4.74/ln(1/eps) law)": the 1/ln(1/eps) scaling is in Woodworth et al. for limit points, and 4.74 is ln(1/theta) (PZ1-2).
- "dormant units static" in the Reading is supported by medians only.
- P0u and P4 appear only on the design, while F2 lists P0u among the readouts "on the 200 targets".
- The statement "the flow is exact on [target 1] ... at 1x, 3x, 10x" is correct (rms 7.0e-6, 6.8e-6, 6.3e-6).
- The new "Limits" sentence ("an unsigned u.v, asymmetric starts or other learning-rate splits were not tried") is answered by section 2: asymmetric starts and rate splits fail. An unsigned u.v is one-sided and cannot fit the design's negative amplitude, and the signed u.v is D2 with a rescaled time.

**8. Literature.** [snippet] Woodworth et al. 2020 (arXiv:2002.09277), full text:
- D = 2 diagonal model: (1/log(1/alpha^2)) Q_alpha(beta) = ||beta||_1 + O(1/log(1/alpha^2)).
- Theorem 2: alpha <= exp(-d/(eps ||beta*||_1)) (exponentially small) is needed to approximate the l1 solution.
- Theorem 3 (D >= 3) and the text: "even at order D = 3, the scale of alpha needed to approximate l1 is polynomial rather than exponential"; increasing D hastens the transition.
- The model w+^D - w-^D with w+(0) = w-(0) = alpha is the study's `Cubic3` at D = 3.
- These are statements about limit points of interpolating gradient flow on a fixed linear design, not trajectories or moving centres.
- Gissin et al. (arXiv:1909.12051, abstract): shallow models need exponentially small initialisation for incremental learning, deeper models polynomial.
- Pesme and Flammarion (arXiv:2304.00488, abstract): depth-2 diagonal linear networks from vanishing initialisation jump saddle to saddle (LARS-like) to the minimum-l1 solution. The abstract does not state the log time rescaling; paper17.tex l.123 does.
- Berthier (arXiv:2208.14673, search-result summary only): incremental learning in the small-initialisation limit, restricted to underparametrised, anti-correlated features. This does not cover the study's positively correlated dictionary.
- Nothing contradicts depth-3 sufficiency or the depth-2 logarithmic law; both are what the literature predicts. The study's contribution is the moving-centre RBF setting and the comparison with the sliding path, not the contrast itself.

**9. Mutants.** Single-point mutants of `src/core/factorized_rbf.py` run against the two test files in my own copy.
- 7d45e0c engine (= 6c87a4e): 66 mutants, 50 killed, 16 survived. Plus 4 constant-dependence mutants (X01-X04), all of which survive.
- 5083a49: 77 mutants (the 66, 7 on the fallback and failure padding, 4 constants), 54 killed, 23 survived. M28 (rtol 1e-9 -> 1e-4) is now killed by the new Radau test. N02, N03, N04 are killed.
- Survivors on 5083a49:
  - Equivalent or immaterial: M33 (ceil -> floor; one digit of one descent cell), M42 (>= -> >), M48 (formed >= thr).
  - Failure bookkeeping without published effect: M29 (failed solve reported ok), M31 and N07 (padding with zeros), N05 (nfev not accumulated), N01 and N06 (the default budget disabled or x1000: tests monkeypatch STIFF_BUDGET = 30, so the default value is untested).
  - Changing printed values or flags when the sentinel cells (B, D, F, G, H) are recomputed:
    - M40 (Adam with the group rates): H P0 rms 7.4e-8 -> 2.0e-5.
    - M53 (dormant max|w| at the end only): 5 of 12 rows, D1 row 0.0125 -> 0.00018.
    - M59 (tau_rel mean instead of max): 9 of 12 rows (P0 1e-4 0.0025 -> 0.00123).
    - M64 (order_equal counts spurious units): D2, D1 and P0 1e-4 rows "yes" -> "no".
    - M66 (sweep relative to min|A*|): 5 of 12 rows.
    - M67 (rms -> mean|r|): 12 of 12 rows.
    - M63 (order_equal ignoring order), M46, M55, M61 (amplitude tolerance variants): flip exactly one stored flag of 1,702 (F target 2, P0).
  - X01-X04: invisible at G = 1 and a_max = 1 by construction.

**10. The coordinator's updates.**
- *Part F stability.* The single-unit estimate lambda = Gm A^2 <(d_c phi)^2> is right to 7%. The full-flow Jacobian at the end state has lambda_max = 35.32, 43.33, 22.78, 34.49, 34.78 (D3; P0 1e-4 gives 35.08, 43.29, 22.75, 33.50, 22.15) for targets 0-4, against the study's 35.30, 43.32, 21.90, 32.33, 32.61 (up to 6.7% low where nearby units couple). Critical steps 2/lambda: 0.0566, 0.0462, 0.0878, 0.0580, 0.0575. Empirically h = 0.04 is stable and accurate, h = 0.05 (h*lambda = 2.17) is unstable. All rows in the table have h*lambda < 2 (maximum 1.77 for the design at h = 0.05, 13% below the limit). Transients with overshoot are not covered by the estimate, but the descent matches the flow throughout.
- *Radau fallback* (`train_flow`: LSODA with a 150,000-evaluation budget per piece; over budget or status < 0, the piece is redone by Radau from its start; both failing records "failed"; `cell()` turns any exception into status "error: ..." with exact = set_equal = False, n_formed = 0).
  - Claim "A, B, D-H, J cannot trigger it": TRUE for A, B, D, F, G, H (flow cells <= 41,637 evaluations, or descent/Adam). FALSE for E: E[3] P0 used 4,288,854 evaluations (360 s) and is exactly the stall the fallback exists for. FALSE for J: J[8] P0 piece 4 exceeds the budget and Radau finishes it in 5,296 evaluations (J[4], J[2], J[11] stay under 115,000 and reproduce bit for bit).
  - Does it change results? Cells that do not trigger it: bit-identical (above). Triggering cells: I[3] P0 (same computation as E[3]) 4,288,854 -> 21,492 evaluations, 360 s -> 14.6 s, seven stored fields differ by <= 5.2e-4 relative (pos_err and amp_err of order 7e-8, tau_rel 4e-5, rms 1.5e-6), nothing at 3 digits; the committed raw holds both versions. J[8] P0 with the time limit lifted: <= 1.3e-5 relative, same flags. I[2] P0 x10: a 1,110 s "timeout" (13.5M evaluations, rms 0.154) became ok in 14 s with 3 formed units, exact. Part C uses it: 14 of the 84 cells I reran contain Radau pieces (budget on targets 82, 85, 93, 100, 103; LSODA "Unexpected istate" on 108, 112, 119, up to 7 of 8 pieces). They agree with all-Radau rtol 1e-11 (five cells): no flag or count differs, gap/disp/tau <= 2e-6 relative (rms of noise-level residuals up to 8e-4).
  - The crash: the 6c87a4e engine raised AttributeError ('list' object has no attribute 'shape') when LSODA failed with an empty output (I reproduced it on D5 at eps 1e-8, design); 5083a49 fixes it (`Y = np.asarray(sol.y)`) and D5 1e-8 finishes in 4 s. LSODA's failure there is the same "Unexpected istate" (status -1, 21 evaluations, `sol.y` a list) seen on target 108.
  - Can an error or timeout be counted as a success? No: see PZ1-5.

## Verdicts per part

**A (reformulation): CONFIRMED.** Identity 8.21e-16; raw vs reduced flow <= 1.43e-7 of max|A*| and <= 9.1e-9 sigma; conserved quantities <= 7.2e-8 [md]; independent derivation and checks at three constant sets [own].

**B (design): CONFIRMED, reproduced independently.** Own D2, D3, P0 reproduce the table (D3 1e-6: 3 formed, gap 5.27, rms 1.26e-7; D2 1e-6: 9, 0.76, 2.8e-3); four integrators agree; nine cells bit-identical. P0 1e-4's fourth unit is kernel 10, the near-tie neighbour of kernel 11 [raw form_order].

**C (200 benchmark targets): F2 conditionally confirmed, F3 confirmed.** All 1,200 cells ok; reference order = stored on 200/200.
- D2: set 0/200 at 1e-6 and 2/200 at 1e-12, exact 0/200, median 8 (3-12) and 5 (2-8) units, gap 0.8 and 1.1. The two D2 "set = path" cells are not exact (errors 0.127 sigma/12% and 0.042 sigma/5.1%).
- P0 and D3: margin table, own (mu = smallest distance 1-|u_j| of a not-yet-active kernel at any entry; `path_duals_nd`):

| margin bin (n) | set mismatches P0 1e-6 / 1e-8, D3 1e-6 / 1e-8 | non-exact (same order) |
|---|---|---|
| < 0.002 (16) | 8 / 6, 7 / 4 | 5 / 3, 4 / 2 |
| 0.002-0.005 (19) | 5 / 1, 3 / 1 | 3 / 0, 1 / 0 |
| 0.005-0.01 (20) | 10 / 2, 6 / 0 | 5 / 0, 2 / 0 |
| 0.01-0.02 (47) | 11 / 0, 8 / 0 | 0 / 0, 1 / 0 |
| 0.02-0.05 (45) | 0 / 0, 0 / 0 | 0 / 0, 0 / 0 |
| >= 0.05 (53) | 0 / 0, 0 / 0 | 0 / 0, 0 / 0 |
| all 200 | 34 / 9, 24 / 5 | 13 / 3, 8 / 2 |

- Order = path also holds 98/98 for margin >= 0.02, and the two order misses with set equal (targets 36 and 42) come from entries closer than the first-order shift (target 42: entries 15.464 and 15.654; tau_rel 0.018 at 1e-8, and the order is right at 1e-10 and 1e-12 [own]).
- Non-exact cells always have an extra formed unit (13, 8 at 1e-6; 3, 2 at 1e-8); none has fewer than K. 34 distinct targets have a set mismatch in at least one P0/D3 cell; the md's "72" counts cells.
- tau err.: median 0.0007/0.0003 (P0), 0.0056/0.0008 (D3). The maxima (0.476, 0.602, 0.408, 0.567, all target 119) are a record-grid artefact; on resolved targets the 95th percentile is 0.24 (P0 1e-6), 0.045 (1e-8), a first-order tail.
- Set mismatches that persist to 1e-8: P0 on target 7 (kernel 15 for 12), D3 on target 9 (kernel 4 for 8); fixed at 1e-10 [own].

**D (depth 2 in eps): CONFIRMED and replicated.** My code reproduces set 0, 0, 4, 9, 12, 17 and exact 0, 0, 6, 13, 15, 20 of 21 [md/own]; qualifications in PZ1-2.

**E (learning rates): exactness CONFIRMED (11/11 at all rates); set identity 8/11 or 9/11.** The mismatching targets are 2, 8, 10 (margins 0.013, 0.006, 0.001). The table prints medians only: the D3 lr = 1 maximum dormant displacement is 0.89 sigma (target 2). E[3] holds the pre-fallback value of a computation also in I[3] (PZ1-3).

**F (fixed-step descent): CONFIRMED after the rerun** (see 10 above). D3 at 1e-6 is exact and has the path's set on the design and targets 1-4 at h = 0.02; P0 1e-6 on the design at h = 0.05 is identical to the flow. The P0 1e-4 rows have no flow comparator in the raw; I computed them (<= 1.8e-5). F[2] (234 s) would time out under the committed limit (PZ1-4).

**G (random starts): masses CONFIRMED (P0 11/11 exact; D3 10/11); amplitude rows are not a depth test.** [own] P0 with amp0 up to 0.3 N(0,1) at 1e-6 and 1e-10 forms exactly 3 units because random weights 1/(gamma|alpha_i(0)|) differ by decades. The paper's depth-2 regime is two kernels per atom with alpha(0) = +-a0 of equal magnitude. [own] P0 on such branches: a0 = 0.3, eps 1e-6 forms 6 units (gap 0.71 decades; D2-like), eps 1e-12 forms 3 (gap 2.2); a0 = 0.1, eps 1e-6 forms 3 with gap 1.9 (D3: 5.3). This is a direct confirmation of the thesis that depth, not the parametrisation, decides.

**H (Adam): as reported** ("reported, not claimed"): all 21 units form at step sizes 1e-3..1e-5. The CLAIMS statement for 1e-6 is wrong for D3 and D2. Cells recompute bit for bit (D3 1e-3).

**I (centres' rate): rerun CONFIRMED.** Fixed centres: 0/11 exact (by construction); D3 11/11 at x0.1, x1, x10; P0 10/11, 11/11, 11/11 (the earlier "10 of 10 that finished" is now 11/11); D2 0/11 at every rate. Only the P0 x10 and P0 x1 target-3 cells changed.

**J (the plane): CONFIRMED with two failures reported honestly.** P0 and D3: set 10/12, exact 11/12 (median gap 4.8 and 4.1). Target 8: a different kernel but exact; target 11: 4 formed. Dormant displacement up to 0.56 sigma (P0) and 0.36 sigma (D3). D2: median 19.5 units, exact 0/12. The reference order is reproduced 12/12 and the sample has 6 of 12 "other" targets (population 31%). J[8] P0 is not reproduced by the committed code on this machine (PZ1-4).

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| PZ1-1 | major | F2 as worded ("exactly the path's units ... at eps 1e-6 and below", "dormant centres staying put") holds only where the path's runner-up margin is large. [md] C: set = path P0 166/200 (1e-6), 191 (1e-8); D3 176, 195; K units on features 187, 197 / 192, 198; max dormant displacement 2.57 and 0.78 sigma (P0), 1.2 and 1.45 (D3); max gap deficit: minimum gap 0.74 and 0.79 decades at 1e-6. [own] 98 targets whose smallest distance to its wall of any not-yet-active kernel at an entry is >= 0.02 (138 if only never-entering kernels count): set, order and exactness 98/98 (138/138) for P0 and D3 at both eps. All 34 (P0) and 24 (D3) mismatches at 1e-6 are below 0.02 (about 20 eps^1/2); the 9 and 5 at 1e-8 are below 0.0057 (about 57 eps^1/2). The substituted kernel is the neighbour (target 8: kernel 15 for 12 for P0 at 1e-8; target 10: kernel 4 for 8 for D3 at 1e-8), still an exact fit, and the set is right at 1e-10. The paper's own finite-eps evidence uses margins 0.063-0.151 (l.553) and says the limit is not uniform in the margin (l.455, 517). Target 59 (margin 0.044): a dormant kernel moves 2.57 sigma at eps 1e-6 (theory: eps * mu^-3). | Add the margin (`path_duals_nd`) to `reference()` and a margin-stratified table to part C; compare with `rich_transit.recovery_threshold`. Reword the brief's F2, CLAIMS section 2 Remark and section 5 "Reading": "K units on the features; the path's kernels when the runner-up margin exceeds about 20-60 eps^1/2; dormant units static to O(eps mu^-3)". |
| PZ1-2 | minor | F3's constants and conditions. [own] The 4.74 equals ln(1/0.009) + 0.03; D2's formation times are all early by the same factor, so entry-time ratios match the path to 0.1% and the "error" is a time-axis normalisation (wall at abs(w) = 1 vs crossing at theta); at theta = 0.5 max abs A* it is 0.052 at 1e-12 (0.170 at theta 0.009). "No practical eps" is the inequality mu > ln(1/theta)/ln(1/eps) (0.34 at 1e-6, 0.17 at 1e-12, 0.043 at 1e-48; no benchmark target above 0.15); it predicts 21/21, 21/21, 20/21, 18/21, 20/21, 20/21 of Part D. On-grid features with overlap <= 0.5 pass at 1e-6; off-grid ones fail at 1e-6 for spacings 0.1-0.33 and pass at 1e-12 only at 0.33. Part D's 21 targets are family A "nearest" only (the easy family); 4/21 already reach the set at 1e-24, 17/21 at 1e-144. The scaling is Woodworth et al.'s hypentropy 1/ln(1/alpha^2) [snippet]. | State the rule and its conditions in the md Part D text and CLAIMS; print the threshold-free gap and top-K exactness next to the study's; add a spacing sweep (cheap). Replace "quantification" by "entry-time analogue of Woodworth et al.". |
| PZ1-3 | minor | Fallback bookkeeping. (a) "A, B, D-H, J cannot trigger it" is false for E[3] (4,288,854 evaluations) and J[8] (see 10 above). (b) The raw holds two versions of one computation: E[3] (LSODA) and I[3] (Radau); seven fields differ by <= 5.2e-4 relative, none at 3 digits. (c) Part C uses it on at least 14 cells (targets 82, 85, 93, 100, 103 by the budget; 108, 112, 119 by LSODA status -1), unrecorded; `nfev` omits the aborted LSODA work (target 82 D3 1e-6 shows < 40,000 evaluations although each Radau piece followed 150,000). (d) The five Radau-computed cells I checked agree with all-Radau rtol 1e-11. | Store `radau_pieces` and the aborted-evaluation count per cell (`train_flow` return and `cell()`); rerun E[3] with 5083a49; put STIFF_BUDGET and MAX_SECS in the md header. |
| PZ1-4 | minor | Wall-clock limit. MAX_SECS was lowered from 400 to 120 s after the runs. Recorded times above 120 s: E[3] 360 s, F[2] 234 s, ten J cells (121-171 s). Rerun with the committed constants on this machine, J[8] P0 ends "timeout": the aborted LSODA piece costs about 90 s, then the check fires (n_formed 2, rms 5.2e-2, against the committed ok, exact, 2.7e-8); with the limit lifted it takes 274 s and agrees to <= 1.3e-5. F[2] would time out at its first record point after 120 s (`train_gd` checks at every record). Status therefore depends on machine speed; the check also fires only between pieces (the old 1,110 s vs 400 s overshoot). | Use evaluation budgets (deterministic) or a limit above the slowest cell; rerun F[2], J[8] with the final constants. |
| PZ1-5 | minor | A non-ok cell is never counted as a success: `cell()` returns exact = set_equal = False, n_formed = 0; every table filters status = ok; the cell record is tested. But D, E, G, I, J drop non-ok cells from the denominators; only C and J print "runs finished"; E, G, I do not; an error cell has secs = 0.0 in the median. No test builds a table from a non-ok cell (a report counting all statuses is not exercised). The committed raw has none. | Count non-ok runs as failures, or print "runs finished" in every table; test a report with one error cell. |
| PZ1-6 | minor | Part C presentation. (a) The failing-cell table lists 40 of 72 rows, all family A and ending at target 68, without saying so; its heading "Targets ...: 72" counts cells (34 distinct targets). (b) "max tau err." is the record grid: first entry tau = 0.00956 against a step of 0.0093 (target 119); an ideal step record gives 0.87 on that grid; 16 targets have a floor above 1%. (c) The reference's `relaxed`, `solver_failed`, `out_of_budget` are discarded (all fine on the 200). (d) P0u, P4, D1 are not in C. | Print all rows or "40 of 72" and the distinct-target count; log-spaced or event-refined record grid and tau err. in units of eps^1/2 on resolved entries; store the reference's flags; add P0u to C or drop it from F2. |
| PZ1-7 | minor | Test suite (section 9): 23 of 77 mutants survive on 5083a49. Rates of P0 are pinned; G and a_max never vary (X01-X04); tau_rel as max, order_equal both ways, dormant max abs w over time, Adam on P0, a failing-solve record, the default STIFF_BUDGET are all unpinned. | Parametrise the F1 test over (m, G, gamma, a_max); synthetic records pinning tau_rel as max, order_equal negative and spurious-unit cases, a time-varying dormant unit; Adam test on P0; a failing-solve test with empty `sol.y`; a test that the default budget fires on a known stall. |
| PZ1-8 | minor | "Depth >= 3" is tested at depth 3 only (P0, P0u and P4 on the design, D3, P0 and D3 in C). [own] D4 and D5 on the design pass (3 units, set, exact; gaps D4 5.4/9.4/13.4, D5 7.8/12.8). | Add D4, D5 to part B. |
| PZ1-9 | minor | Part G's random alpha_i(0) rows are not the paper's logarithmic-regime test; [own] branch pairs (alpha(0) = +-a0): D2-like at a0 = 0.3, eps 1e-6 (6 units, gap 0.71). | Replace the amplitude rows by +-a0 pairs. |
| PZ1-10 | cosmetic | md and CLAIMS wording: header "10^3" should read "10"; Part D text "later or earlier" should read "earlier"; "order = path" ignores spurious units (D2 192/200 yes); CLAIMS: Adam at 1e-6 (only P0 has nothing formed), test count 31 vs 33, "quantification new" (PZ1-2), "dormant units static" medians only. | Edit the header sentence in `report()` and CLAIMS section 5. |
| PZ1-11 | minor (paper text) | paper17.tex l.127: "a metric of Fisher-Rao type in the mass, not a Euclidean one". In (alpha, nu = sqrt(beta), c) the metric is Euclidean (dbeta^2/(gamma m beta) = 4 dnu^2/(gamma m)); the stated reason for lying outside AGF's setting fails (the other reasons given there remain). | Add this remark to Lemma (o) and revise the sentence. |

## Sources

- [md] results/chaos/parametrization.md at 7d45e0c, 6c87a4e, 5083a49.
- [raw] results/chaos/parametrization_raw.npz at the same commits; results/chaos/transit_certificate_raw.npz (families A, B, C, D: 400, 160, 120, 474 targets; 1-D 1,034); transit_limit.md l.7 and transit_certificate.md l.145-159 (CLAIMS section 0 spot-checks: 12/48, 763/391 hold).
- Code: src/core/factorized_rbf.py, experiments/exp_parametrization.py, tests/test_factorized_rbf.py, tests/test_parametrization.py (all three commits); rich_transit.py, rich_transit_nd.py, mirror_free.py.
- Paper: manuscript/paper17/paper17.tex l.100-135, 150-220, 258-300, 395-455, 517-553 (the lemma (o), Prop. unit, Prop. branch, Thm sliding limit).
- CLAIMS.md at /tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/paper/ (read last; re-read at its 11:33 update).
- [own] my scripts in W: myfr.py (independent engine), f1_own.py, f1_general.py, d2attempts.py, d2_21.py, d2_margin.py, coherence.py, coherence2.py, branch_test.py, stab.py, gd_check.py, integr_check.py, bitident.py, c_bitident.py, c_radau_count.py, radau_validate.py, fallback_J.py, j8_nobudget.py, margins.py, margins_all.py, margins_runnerup.py, analyze_C.py, ref200.py, grid_floor.py, mutate3.py (+ mutate.py, mutate_sentinel.py), nfev_scan.py, cmp_raw23.py.
- [snippet] arXiv:2002.09277 (PDF text: eq. before Theorem 2, Theorems 2 and 3, the paragraph after Theorem 3); arXiv:1909.12051 (abstract via fetch); arXiv:2304.00488 (abstract via fetch); arXiv:2208.14673 (search-result summary; the JMLR page was unreachable from the session).
- [memory]: explicit Euler is stable for h < 2/lambda_max; the Fisher-Rao metric on (0, infinity) is flat in sqrt(beta); LSODA's automatic stiffness switching can fail (the repository documents it in `relax_chunk`).

## Lesson

A limit statement checked on a random sample must be reported against the parameter that measures the distance from the limit: here the runner-up margin against eps^1/2 for depth 3 and against ln(1/theta)/ln(1/eps) for depth 2. A bare fraction (166 of 200) hides a clean 98 of 98, and "depth 2 fails" hides that the verdict is the same inequality with a logarithmic instead of a power-law right-hand side. Every budget that can change a cell's status (wall-clock, evaluations, a solver switch) has to be fixed before the runs and recorded in each cell. Here a fallback, a lowered time limit and a stale cell made the committed raw a mixture that no single command reproduces, although no conclusion moved.
