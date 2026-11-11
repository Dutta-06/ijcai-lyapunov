# REVIEW_roundT3E: referee report on the widening of the entry-layer theory to general Grams (theory round 3, direction E)

Archived verbatim from the referee's hand-back (2026-11-10). The referee was read-only for the repository; its scratch code is outside the repository.

REFEREE REPORT, theory3 direction E (general Gram). Work files: referee_T3E_work/ (a_gram.py ... m_nd.py, pytest_core.log). Nothing was modified in the repository.

BOTTOM LINE. The mathematics of Prop. G, Cor. G1 and Thm G holds as stated, for the continuum (L2(R^n)) network. I found no step that silently used a Gaussian property. The numbers reproduce independently. The problems are scope and wording, one real (small) code defect in the T1 generalisation, and a few false or inconsistent sentences.

FINDINGS (ranked)

1. [major, wording/scope] Claim (iv) holds for the continuum inner product only. Gamma = autocorrelation of exp(-|x|) is the continuum Gram. Under sampled data the Gram of the Laplace basis has a cusp in c at every data point, and the sampled single-unit optimum is displaced from c*. Evidence (f_sampled.py), Laplace n=1, data spacing 0.28 sigma_eq (the repo's 61-point window): the slope jump of the sampled own term at a data point is 0.20, against the smooth restoring force psi(0.417) = 0.22. The argmax of the sampled own term sits at c = 0.073 for c* = 0 (phase 0.37). A Gaussian is spectrally accurate under sampling; the Laplace basis is only O(h). No non-Gaussian network exists in the repo (grep finds none), so the "network" evidence for these Grams is the continuum simulator `theory3e_path.simulate_network`. The paper's own earlier referee (REVIEW_roundTH, finding 1) already flagged continuum versus sampled.
Cheapest fix: put "continuum inner product, data spacing small against the kernel's correlation length" in the hypothesis of every non-Gaussian statement. Add: "for sampled data the Laplace basis behaves like a cusped Gram near data points and Prop. G does not apply".
Also say "Laplace BASIS function exp(-|x|/l)". If a reader takes "Laplace kernel" to mean the Gram exp(-d), the case is Prop. G'' (cusp), not Prop. G.

2. [moderate, code; verdicts unchanged] `theory3e_t1gen` underestimates the Gram constants at arguments below L_MIN = 0.4. `Gram.tail(order, l <= 0.4)` returns the sup over [0.4, inf), not over [l, inf). `exp_theory_t1` evaluates `Gam(X)` and `G2_sup(X)` at X = 0 for the entrants and their lattice neighbours (X = max(0, dist - r - dc), and X[ent,i] = 0), which the replacement maps to TAIL0 and `tail(2, .)`.
Measured at argument 0 (g_t1.py): TAIL0 = 0.961 (Gaussian), 0.972 (Laplace), 0.975 (Matern-3/2), where the true value is 1. G2_sup = 0.442 / 0.270 / 0.415, where the true value is 0.5 (Laplace 46% low). This probably also explains "D within 3-4%" in the Gaussian regression.
I patched `tail` in scratch (`tail(l <= 0.4)` = max(table, 1, M1, M2)) and re-ran section H for rows 0, 3, 6, 15. The certificate strings are unchanged (..YYY, ...YY, ...YY, ..YYY), so section H survives.
Fix: apply that patch in `Gram.tail` and add a test at l < 0.4. The existing tail test only uses l >= 0.7.

3. [minor, false sentence] "every kernel is below 1e-8 beyond 40 sigma" is false for the rational-quadratic Gram. At 40 sigma_eq, |G|, |G'|, |G''| are 2.5e-5, 2.5e-6, 3.1e-7 (beta = 2) and 4.1e-7, 6.1e-8, 1.1e-8 (beta = 3). It is harmless numerically: the table is a suffix maximum and the thresholds use l = Delta - 2r of about 4-12, so only monotonicity beyond 40 matters. Monotonicity is only numerically checked on [5, 45] (c_struct.py). The note has a hand proof only for half-integer Matern. Gaussian, RQ and the integer-nu Bessel Matern are "easy but not written".
Fix: replace "<1e-8" by "assumed non-increasing; hand-proved for ...". Do not call the tabulated rows "certified" without that clause.

4. [minor, answers attack (b)] The dormant mass Omega enters kernel-free (|Gamma| <= 1, |Gamma'| <= M_1), so a heavy tail never helps or hurts through Omega. The sums S_j = (K-1) m_j(Delta - 2r) are valid for K features at pairwise distance >= Delta. For a lattice of features the Babel sum converges only if the algebraic exponent s > n. Cor. G1's lattice remark and the decay form "(C(K-1)R^2/c_*)^(1/s)" need the hypothesis 2*beta > n for the RQ case, or K-1 must be read literally.
T1 (section H) was run only for four n=1 Grams (Gaussian, Laplace, Matern-3/2, Wendland). No algebraic-decay Gram was run in T1, and RQ appears only in the box and path parts. T1 "for general Grams" is therefore not established for slowly decaying kernels. The Omega = 0.02 of the certificate is a parameter checked on the target, not a derived constant.

5. [moderate, wording] The "network following" evidence does not test certification.
- Section G has 31 targets: 16 certified and 15 uncertified.
- Formed exactly the recursion's units: at eps 1e-4, 4/16 certified vs 3/15 uncertified; at 1e-6, 16/16 vs 13/15; at 1e-8, 15/15 vs 14/15. The n=2 target has no 1e-8 run.
- Uncertified targets follow as well as certified ones. The only failures (Laplace R=2.4, margin 0.0065, fails at 1e-6 and follows at 1e-8; Wendland, margin 0.0000, never follows) are small dual margins, i.e. hypothesis (b). The box threshold Delta_entry is sufficient and conservative, not discriminating.
- Section H has 4 Grams, one target each, 5 decades. Only about 10 cells have the certificate true, so "0 violations" is a statement about those cells.
- There is one planar target and none in n=3.
- Hypothesis (b) margins (0.006-0.12) are not tested for likelihood.

6. [minor, internal inconsistencies]
- NOTE 9.1 says the network follows from "1e-6 for Wendland". The report's table H gives eps_max = 1e-5 for Wendland.
- NOTE 10 says "1 of 8 failures at Delta_entry - 2 (Gaussian)". The F table also shows Wendland at -2: offsets in box 6/8, coefficients 15/16, monotone 7/8.
- Laplace n=2 leaves the box at only Delta_entry - 1 (7/8), so "the box is conservative by 1-2 sigma" is not uniform. The sentence "the box is conservative, not a transition" is unsupported for n=2.

7. [minor, units]
- "sigma-equivalent" (curvature-matched, -Gamma''(0) = 1/2) is not a width of anything for non-Gaussian bases. State it once and define d_0 = 0.4167 as "half the dictionary spacing 0.833 in these units". The thresholds are conditional on that spacing.
- The cusp row uses half-height sigma (x 2.402) and d_0 = 0.4167/2.402. Comparing "cusp 3.8 < Matern-3/2 6.72" mixes conventions. Converted, 6.72 eq is about 4.7 half-height sigma, and d_0 differs (0.29 vs 0.417), so "smaller because (N) is absent" is not like-for-like.
- Cor. G2 (anisotropic Gaussian): d_0 is taken as 0.4167 in whitened coordinates for every cond(S). For a fixed Euclidean lattice, |d'| <= |d|/sqrt(lambda_min(S)), so the table is for a lattice whitened with S.

8. [minor] The "rational-quadratic" and "Wendland" rows are Gram classes. RQ has no natural basis function (it is defined through sqrt of its spectral density). Wendland is the basis, and its Gram is the support-2 autocorrelation. Do not call them "kernels" without that.

9. [minor] Prop. G does not use Gamma''' (only continuity of Gamma''). It needs Gamma in C^2, so the ODE field is C^1. For the Laplace basis in n=1, Gamma''' jumps (from 2 to -2), so Gamma is C^{2,1} but not C^3 (numerical G''' at 0+ = 1.9955). Say "C^2, Gamma'' Lipschitz" and do not claim smoother.

VERIFIED (how)
- (a) Gram identification. I computed autocorrelations by direct quadrature (a_gram.py) and compared them with the closed forms:
 - exp(-|x|), n=1: (1+d)e^{-d}, agreement 1e-16.
 - exp(-|x|), n=2: Matern-2, c d^2 K_2(d), agreement 3e-12.
 - exp(-|x|), n=3: (1+d+d^2/3)e^{-d}, agreement 1e-16.
 - Matern-3/2 basis, n=1: Matern-7/2, agreement 1e-16.
 - Wendland n=1: 12 digits, my Gauss-Legendre autocorrelation vs `fG`, `fdG`, `fd2G`.
 - Wendland n=3: bipolar quadrature vs `rawf`, 1e-13.
 - `family()` reproduces these in eq units. Finite-difference -Gamma''(0) = 0.5 for all n.
- (i) Gamma in C^2 is equivalent to phi in H^1 (Bochner/Plancherel, second moment of |phi_hat|^2). The "-Gamma''(0) > 0" clause is automatic for phi != 0.
- (d) Thresholds, independent transcription from the NOTE formulas with my own Gram definitions and a finer r grid (d_thr.py, m_nd.py, j_wend.py):
 - Gaussian 4.56 / 5.44 / 6.68 / 4.98 / 5.78. With the row-sum Hc: 4.66 / 5.50 / 6.70 / 5.08 / 5.84, the paper's table exactly.
 - Laplace n=1 6.72 / 8.58 / 12.68 / 7.96 / 9.78; n=2 5.70 / 7.26 / 10.44; n=3 5.34 / 6.74 / 9.52.
 - Matern-3/2 basis n=1 5.04 / 6.30 / 8.66 / 5.74 / 6.94; Matern-5/2 4.82 / 5.96 / 7.96; RQ beta=3 4.96 / 6.34 / 9.54.
 - Wendland n=1 (own exact autocorrelation) 4.48 / 5.28 / 6.10.
 - All agree with the report to within the r-grid resolution (0.02).
- Cor. G1: I re-derived the constants 1.7R, 2.7, 4.59, 7.29, 11.88 by hand and they are correct. c_* recomputed (i_cstar.py): Gaussian R=1: 0.0349; R=2.4: 0.0311; R=10: 0.0165. The note's "0.017" for R=10 is rounded up.
- Sign structure (c_struct.py): psi > 0 on (0, 1.35] for all 17 Grams. Gamma'' < 0 up to 1.265 (RQ beta=2) ... 1.44 (Wendland), so rho_min > 0 on the whole range used (r <= 1.35, dtilde < 0.28 in the cases checked). No sign-change or positivity assumption is violated for Wendland, Matern or RQ. Gamma >= 0 and Gamma' <= 0 on the whole axis for every Gram checked.
- Adversarial check of the proof's steps on F itself (e_adv.py, n=1, K=2,3, random signs and amplitude placement, box faces and corners). At the threshold Delta with the optimal r (r = 0.417 in all 8 cases checked), over 600 random trials per case, for Gaussian, Laplace, Matern-3/2, Matern-5/2 and RQ:
 - max|P_i| <= E_v and max|C_i| <= E_c with slack factors 2-6 (e.g. Laplace E_v = 0.175, attained 0.031).
 - min eigenvalue of the full Hessian over D* is positive (0.23-0.48, refined by Nelder-Mead), so (N) is satisfied with margin.
- Hessian block structure (own block [[1, s],[s, t]], tangential stiffness psi/d, Schur-test bound) re-derived: correct. `--hc norm2` is a valid bound.
- Tests: `tests/test_theory3e_core.py`: 16 passed (pytest_core.log).
- Anisotropic case: the condition psi(r) > cond(S) E_c follows from d/dt of |d'|^2/2 with mobility S^{-1}; correct.

COULD NOT CHECK
- Full `tests/test_theory3e.py` (not run; heavy). Ball-arithmetic enclosures were re-derived only by sampled bounds, not audited line by line.
- Monotonicity of the envelopes beyond 40 sigma (hand proofs).
- T1 proof "uses the kernel only through the replaced quantities": not re-read line by line. I only checked the residual Gaussian constants in `exp_theory_t1.py` and found the item 2 defect.
- Prop. G'' (cusp) and n=2 Wendland (not computed anywhere); the sampled network with any non-Gaussian kernel (none exists in the repo); the likelihood of hypotheses (a), (b).

WHAT THE PAPER MAY CLAIM (with hypotheses)
"For the continuum network (L2(R^n) inner product, equal walls, noiseless), the fit depends on the basis function only through the Gram Gamma = autocorrelation. If the basis function is in H^1 (equivalently Gamma in C^2, radial), then Prop. G and Thm G hold for conditions (L), (C), (N) evaluated from the Gram's tail envelopes. This applies to the Gaussian, the Matern (including the exponential basis function exp(-|x|), whose Gram is Matern-(1+n/2)), the Wendland basis function and the RQ Gram, in all cases as computer-assisted ball-arithmetic sufficient thresholds (in curvature-matched sigma-equivalents, starting offset d_0 = 0.4167), conditional on the tail envelopes being non-increasing beyond 40 sigma_eq. Hypotheses (a) and (b) of Thm G are target-dependent finite checks. Cor. G1: if (K-1) R^2 max(m_0, m_1, m_2)(Delta-2r) <= c_*(r) (c_* = 0.035 at R=1 for the Gaussian), then (L), (C), (N) hold, for K features at pairwise distance >= Delta. For features that are not pairwise equidistant, read it with the Babel sums."
Also allowed: "The T1 certificate evaluates and holds on one n=1 target each for the Gaussian, Laplace-basis, Matern-3/2 and Wendland Grams, with 0 violations on the 10 cells where it is true (computer-assisted; proof outline only)."

MUST NOT CLAIM
- That Laplace or Matern networks on sampled data are covered; that the cusped Gram exp(-d) is covered by Prop. G (only the Filippov statement G'' for the idealised limit object); that "the Laplace kernel is smooth".
- "Theorem T1 for general kernels" (only a computed certificate instance for four n=1 Grams, with item 2 corrected).
- "Certified targets follow the path" as a validation (uncertified follow equally; the network results do not test the threshold).
- That every Gram is below 1e-8 beyond 40 sigma; that "sigma" is a width for non-Gaussian bases; that Delta_entry is a transition separation (it is a sufficient, conservative bound).
- "Any radial Gamma": the statement is for Gamma in C^2 with the (L), (C), (N) numbers actually evaluated.
