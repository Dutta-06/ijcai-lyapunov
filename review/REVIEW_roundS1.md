# REVIEW_roundS1: referee report on the supplement's proof sketches (S1) and the main text's theory section

Archived verbatim from the referee's hand-back (2026-11-10). The referee left the repository untouched; scratch code is outside it. File names refer to the private paper drafts (not in any repository).

REFEREE REPORT on S1 (supp_*.tex) and sec_theory.tex. Repo untouched. Scratch: referee_S1_work/ (sign.py cex.py adam.py hyb.py pair.py unit.py curve.py dorm.py cusp.py gauss.py).

Verdict: no blocker. The T1 gating is correct and present. I found 9 major wording/gap items, none of them an arithmetic failure of a theorem, plus minors. The F and G'' proofs re-derive, with the gaps listed below.

MAJOR
1. supp_optim.tex, Prop F1b/F1d status line. It says (iii)/(iv) are "proved in continuous time, or in discrete time for lr<<s". REVIEW_roundT3BDF F1 allows continuous time only.
 - I ran global-norm GD on the referee's 3-unit counterexample. The flow gives T=331.8, 664.3, 696.4 (matches the text). lr=1e-3 gives order (1,3,2); lr=2e-5=s/50 gives (1,2,3).
 - That is one test, not a proof.
 - (iv) "same orbits" is exact only for beta1=0 (note F1d(iv)); with beta1>0 it holds only as lr->0.
 - Fix: "continuous time (lr->0); discrete lr<<s observed on one 3-unit case".
2. supp_optim F1d(0) and sec_theory Prop 3 ("erases the ordering by the drive"). The proof is an induction that gives identical trajectories only for exactly proportional histories (frozen drive). In the coupled network it is exact for sign descent (Prop F1a) only. For Adam and RMSprop it is numerical.
 - The window theorem is vacuous (72.7 vs spread); observed windows are 1.0-1.6.
 - RMSprop needs eps=0 (not stated).
 - "mirror coordinates" is an untested illustration (note F gap 4).
 - Fix: "exact for frozen drives; for the coupled network proved for sign descent, observed for Adam".
3. supp_optim Prop F2(b), listed under "Status: proved".
 - The start (0,a0) being "the symmetric start a0/2 up to O(a0^2)" is the note's gap 5 (first order/numerical).
 - c_* is never defined (main text too). The hypothesis that both layers are identical units with identical rates is dropped. The result needs the normalisation beta=s^2+2alpha^2, w=beta*alpha, rates 1.
 - Cheapest fix, which I verified in pair.py: alpha1'=Phi*alpha2, alpha2'=Phi*alpha1 with Phi=beta1*beta2*(c*-w1*w2), so alpha=(a0*sinh tau, a0*cosh tau) and T=int dtau/Phi exactly. alpha2^2-alpha1^2 is conserved to 1e-15.
 - Numerics at s=0.1-0.3, a0=1e-6-1e-9: T(0,a0)=T_sym(a0/2) to 1e-11, and T*s^4*c* - ln(s/a0) = -0.153 for every s, a0. So the O(1) is uniform. T s^4/ln = 0.987-0.992, which matches the referee's 0.96-0.99.
4. supp_optim F1a proof sketch. "sgn(g_i) is constant until the first crossing (the drive keeps its sign by the hypothesis)" skips the key step. The step is the envelope lemma: for every dictionary unit k, |w_k(t)|<=omega(t)=6s^2*lr*t+2(lr*t)^3, hence |g_i(t)-g_i(0)|<=theta*||G_i||_1<|g_i(0)| for t<T_env. All dictionary units must be sign-descent D3 units from the common start, which should be said.
 - Numerics (sign.py, three (s,lr,theta) settings, 8-unit random Gram) confirm the proposition. Every unit meeting the hypothesis crosses exactly at T_env (170, 80, 1641). The others cross later (414, 514, 4079).
 - The mass-spread bound eps|b_i-b_j|/(a_0 lr)+1 also holds (spread 2 vs bound 3.2). Derivation: |da/dx|<=1/a_0 on a^3+3xa=theta/2.
5. supp_optim Prop F2(c) and sec_theory Prop escape (c). The note's hypotheses were dropped.
 - Needed: nonnegative drives gamma_l>=0 (b must stay nondecreasing; n_+ counts positive drives), and |gamma_l(t)-gamma_l(0)|<=delta*gamma_l(0) for ALL units, not "its drive".
 - Needed: one-sided cubic first layer and rates 1 (kappa^2=lr_b/3).
 - Conservation holds for any drive signs. Bounds re-derived by hand and checked on 3 ODE cases (hyb.py): T lies in [lower, upper], the conserved quantity is good to 1e-11, and the lower bound is tight for equal drives.
6. supp_readouts_thy Cor. fixed. The text says the depth-2 dual error falls "as 1/ln(1/eps) (0.133, 0.0717, 0.0306, 0.00999 at 1e-6, 1e-12, 1e-24, 1e-48)".
 - Those are repository numbers (theory3d.md Part L). Their error x ln(1/eps) is 1.84, 1.98, 1.69, 1.10, 0.85: faster than 1/ln (note D gap 4, "not explained").
 - The parenthesis attributes them to the referee. The referee's own run gave 0.155, 0.078, 0.039, 0.020 (exactly 1/ln). The 0.509/0.500 slopes are also repository values (the referee's slope was 0.49).
 - Fix: quote one set and say "bounded by", or give both sets.
7. supp_readouts_thy Thm depth2, last sentence: "dual at which the unit reaches amplitude A is lambda_2(1+ln(A/b)/ln(1/(eps b)))" is wrong for b!=1. With lambda_2=ln(1/(eps b))/(4*ell) it should be lambda_2(1+ln A/ln(1/(eps b))). The note's form is in b=1-wall units: ln(A/b)/ln(1/eps).
8. supp_t1.tex status and cross-references. "Status of the implication: proved by hand."
 - REVIEW_roundT3N lists as could-not-check the hand proofs of sections 2-6 other than the displacement term, Lemmas E and M and the gating. The supplement itself concedes this in its last paragraph, so the first status line contradicts it.
 - The gated Lemma cites Lemma 3.3(1), 5.3, 5.4, 3.1(a), "pre" and S_i^-, none of which are defined or proved in the supplement (they live in THEORY_T1.md).
 - The theorem hypothesis block still carries "n=2" and "Delta_T1 ... K=64, 6.9-8.5 sigma", taken from the SUPERSEDED theory_t1b.md. T3N and T3AC reviewed only the n=1, K<=3, eight-target system.
 - Fix: "proved by hand (referee-re-derived: displacement lemma, E, M, gating)", define or inline the cited lemmas, and mark the n=2 and K>3 content as float/outline.
9. sec_theory Thm path.
 - "a basis function in H^1, so Gamma in C^2 with Gamma'' Lipschitz": H^1 gives only C^2. Lipschitz Gamma'' is extra (REVIEW_roundT3E finding 9). Counterexample: |phi_hat|^2~|xi|^-3.5 in n=1. Prop G uses only continuity of Gamma''.
 - The status omits the referee-required clause "tail envelopes assumed non-increasing beyond 40 sigma_eq". The clause is present only in the tab:gram caption and the Not-claimed paragraph. It matters for the RQ row; it is hand-proved for half-integer Matern and trivial for Gaussian and Wendland.
 - "6.72 (Laplace basis ... Gram Matern-3/2)" holds for n=1 only (n=2: 5.70, n=3: 5.34).

MINOR
10. S1.2 setting says depth 1 starts at u=v=s, but the Lemma gives w=ell*eta, which needs u(0)=0 (the repository's Linear has w(0)=0). The same text is in sec_theory.
11. supp_optim F2a sketch writes "z'=D z^{D-1} g(...)". For D independent factors it is z'=lr*z^{D-1}(g-z^D). The tied-parameter form gives T=s^{-(D-2)}/(D(D-2)g), not the stated constant.
 - The hypothesis z_theta^D<g and the model F=1/2(w-g)^2 are not stated ("drive g" is the target value).
 - Add the note's remark that the repository's D3 unit (u cubed) has the same exponent but constant 1/(3sg).
12. supp_gram Prop G'' (cusp).
 - I re-derived the proof: faces by Nagumo, sliding because |C_i|<=E_c<psi_min<=kappa_0, and finite-time locking. It holds.
 - Gaps: v_min is undefined in "<= r/(mu_c v_min(psi_min-E_c))". It is v_0 for the entrant, which tends to 0 as v_0 tends to 0+. Say the bound holds after the entrant passes v_lo.
 - The cusp thresholds need d_0=0.4167/2.402=0.1735 (correlation-length units). That is not stated.
 - Recomputed (cusp.py): 1.565, 2.255, 3.56, 2.26, 2.945 vs the stated 1.58, 2.26, 3.58, 2.28, 2.96. These agree to the 0.02 grid.
13. supp_gram: E_cv is used in D but never defined (it equals E_c).
14. supp_optim "What restores selection": the heading says "Pair, rho=0.5" but lists rho=0.2, 0.5, 0.8.
 - The design sentence "shared second moment and normalised GD form the path's 3 units" is false for P0 normalised GD in tab:optF (1 unit, rms 0.31, not converged).
 - P0 Adam 1e-4: "3 units" has rms 0.194, not converged.
15. sec_theory Sec. 3.4 example. "Six units with drives spanning 1e5: window 1 for sign descent and Adam" is the orthonormal-atoms diagonal model, threshold 1e-6, lr 1e-4, frozen-drive regime, Adam with eps_adam=1e-30. Default Adam is 1.02-1.61 and fails at tiny drives (adam.py: eps_adam=1e-8 gives -1 at scale 1e-8). Add the qualifier.
16. sec_theory Thm noise and N-G: "S(E) holds up to s/A_min=..." means for all smaller s. That needs monotonicity of S in the gauges; the referee only tested it (504 random pairs). B6 said to drop the monotonicity claim. The theorem statement in supp_noise still says "every bound increasing in the gauges". Reword to "at s_crit".
17. sec_intro tabular, last row says "depth test: exponents confirmed". REVIEW_roundP "must not claim" bars this without the equal-mass qualifier and without "check of the closed form".
 - sec_exp is already fixed.
 - The abstract's "(O) ... two layers from a pure small start escape only after eps^-2 ln(1/eps)" lacks "gradient descent, exact zeros" (F4).
 - The abstract's "(proved in the continuum)" omits "thresholds computer-assisted".
18. supp_optim F1a text: "T->ceil((theta/2)^{1/3}/lr)" is correct (checked). The note's t*=37640 for C=7.27 differs from my argmax 35803 (same sup; irrelevant).

VERIFIED (how)
- Exact unit (scripts unit.py, curve.py). ODE for p=2,3,4,5, ell=1, 12, with a sign-changing drive: relative error <=3e-11 vs W_p; lambda=s^{-(p-2)}/(ell*p*(p-2)), q=p/(p-2), depth-2 sinh form, bounds s^p[mu^-q-(2-mu)^-q], 0.963/0.70 ratios.
- Paper-unit Lemma curve (i)-(iii) to 1e-15; wall 7.0248/sqrt(eps).
- Paper-unit dormant constant 8/pi^3*a_max/(m*sqrt(2*v_s)): ratio 0.9996-1.0 for mu<=0.2, 0.969 at mu=0.5. It is a mu->0 statement.
- Sign descent T_i(theta)=min{t:6s^2*lr*t+2(lr*t)^3>=theta}: derived (u=s+a, v=s-a, w=2(3s^2a+a^3)) and simulated. The six-unit table values 78 and 80 reproduce.
- Pair peak rho/(1+rho) re-derived from g_2=rho-(1+rho)omega.
- Adam eps_adam=0: 511 steps for scale 1e-8..1e8. C=7.2703, C/(1-beta1)=72.7. Cauchy-Schwarz and the quasi-monotone lower bound (1-b1)(1-eta) re-derived.
- Referee's 3-unit counterexample: flow T=331.8/664.3/696.4, global-norm GD (1,3,2) at 1e-3 and (1,2,3) at s/50.
- Hybrid conservation and bounds (hyb.py).
- Two-layer pair conservation and law (pair.py).
- F2a quadrature by hand.
- Depth-1 formation bound re-derived. Drift bound Psi_p re-derived.
- Gaussian thresholds, independent float implementation of (L),(C),(N) (gauss.py): 4.56, 5.44, 6.70, 5.00, 5.80 exactly, plus n=2,3 offsets 4.92 and 5.18 (text 4.90, 5.18). Cusp row (cusp.py).
- Cor. G1 constants 1.7R, 2.7, 4.59, 7.29, 11.88 re-derived.
- Table tab_gram rows match the referee's recomputations to the 0.02 grid.
- T1 gating: the supplement contains the gated statement (hinc with "started"/"not started" windows, wstart_i=t_i-S_i^-) and the referee's three-step argument (pre before t_a; limit coefficient 0 before t_i; t_a>=t_i-S^-, centre within r). It matches REVIEW_roundT3N finding 1 and the note. The numbers 131/254, 108/456, 0.67/1.87/3.30 and eps_top, 471/771 (227/203/41), 362/364 match the referee verbatim sentences.
- Noise (supp_noise vs note B and REVIEW_roundT3BDF B): statement, statuses, 24 targets, 288 rows, ratios 0.998/0.641/0.998/0.0407, the tab_noiseN values; B3-B6 are reflected except item 16.
- D "may use" sentences are reproduced faithfully.

COULD NOT CHECK
- T1 hand proofs, sections 2-6 (only the gating).
- Ball arithmetic, Krawczyk and python-flint soundness.
- Lemma E, the noise hand proofs and the N-G net argument; torch KAN; design runs with large eps_adam.
- The shared-v Adam counterexample (only the global-norm variant replicated).
- Wendland, Matern and RQ thresholds (only Gaussian and cusp recomputed); the T1 certificate for other Grams.
- The tab:t1grid float grid.
- supp_sep, supp_entry and S1.5 numbers (skimmed only).
- Whether supp.pdf renders; Section S4 and the P/Q experimental numbers (only the allowed/must-not wording checked against sec_intro, abstract and sec_exp).
