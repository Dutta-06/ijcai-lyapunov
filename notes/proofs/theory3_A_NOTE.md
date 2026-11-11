# THEORY3A: the small-eps tail of Theorem T1 (closing "S(eps) holds for all eps <= eps_0")

Working note, round theory3, direction A.  Setting, notation and the system S(eps) are those of THEORY_T1.md (sections 2-6, 12); code in the repository:
`experiments/theory3a_tail.py` (the tail verifier), `theory3a_chain.py` (the displacement accounting made sound), `theory3a_run.py` (driver, report `results/synthetic/theory3a.md`,
raw `theory3a_raw.npz`, `--from-raw` byte for byte), tests `tests/theory3a_test.py`, mutants `tools/mutants/theory3a_tail_own.py`.

Status words: PROVED = argument written here, each step checked, sanity-tested numerically; COMPUTER-ASSISTED = a finite computation in python-flint balls (100 bits) with the limit object's data
(critical points, entry times, duals) as floating-point inputs taken as exact (as in round T1b: this is NOT enclosed here either); OUTLINE; FAILED.

## 0. What is obtained

* **(A) Lemma E (PROVED).**  For every eps in (0, 2 v_s) and every distance s > sigma_0 = a_max/(2 G v_s) = 20 from the wall,  B_eps(s) <= C/(s (s^2 - sigma_0^2)),  C = 4/(G gamma^2 m) = 100/3.
  This is Lemma 12.4 of THEORY_T1.md, re-proved here from Lemma 2.2 (section 2), with the Adamovic-Mitrinovic inequality proved (section 2.2) rather than quoted.
* **(B) Lemma M (PROVED): the signs.**  For fixed x > 0 the functions q_eps, w_eps, M_eps (x fixed) increase with eps and s_eps decreases; hence the entry positions x_1(eps), x_#(eps) DEcrease as eps increases.
  Consequence (**Proposition N, the exact obstruction**): the ingredient pre(eps) = Psi(x_1(eps))/(v_lo - w_1(eps)) has the numerator Psi(x_1(eps)) INCREASING as eps decreases, so the
  value at eps_1 does NOT bound the value at smaller eps; numerically pre, T_3, Env, kappa_*, m_V all move in the adverse direction by a relative amount ~ c*eps (c = 1.6 for pre).  So
  "S(eps_1) => S(eps) for eps < eps_1" cannot be obtained ingredient by ingredient from the values at eps_1.  (The adverse excess is tiny, 1.6 eps, but it is not signed in the right direction,
  and the proof's strict test E_j < D_j/1.0001 has no room for it unless D_j is chosen with room.)  T_1, A_tr, pre_drift move favourably (larger at larger eps), numerically.
* **(C) Theorem TAIL (COMPUTER-ASSISTED).**  For each of the 8 targets of Part Q of theory_t1b.md there is eps_tail (table of theory3a.md, Part T: 10^-6.25 ... 10^-8.5) and one set of bounds D_j
  such that the system S(eps) holds with these D_j for EVERY real eps in (0, eps_tail].  The majorant is the supremum over eps <= eps_tail of each ingredient, computed by evaluating the ingredient on
  the ball h in [0, h_tail], h = 1/kappa = sqrt(eps/(2 v_s - eps)); the dormant side is bounded through Lemma E in powers of L = 1/lambda.  This is the "monotone majorant" form of the statement and it covers
  eps < 10^-24 (and every eps in (0, 10^-24] simultaneously, no grid, no chain).
* **(D) Corollary (COMPUTER-ASSISTED).**  With the chain of Part C (balls from eps_top down to eps_tail, the verifier of round T1b with the displacement accounting corrected, see (E)), S(eps) holds for every real
  eps in (0, eps_top] on those 8 targets: the interval cover of round T1b ended at 10^-24, the tail removes that end.
* **(E) A defect of the verifier of round T1b found and repaired (PROVED defect, repaired in `theory3a_chain.py`).**  Lemma 5.2 bounds the displacement of a dormant kernel by
  mu_c ( int_0^t H^(st(t')) B(t') dt' + peak [...] ) with the stage-resolved H.  `cert` and `rigorous.verify` use mu_c ( H^(st(t)) I(t) + hinc ) with the current stage's H times the cumulative mass I(t),
  which forgets the weight H^(s) of the mass accumulated in earlier stages s < st(t) (H^(s) decreases with s).  The effect (K = 2 targets, eps near the top): dcmax +20..+50 %, E_j + 1..4 %; all 8 covers are
  re-done with the repaired accounting and close.  The stored covers of theory_t1b.md were computed with the legacy accounting.
* **What is NOT obtained:** a proof of S(eps_1) => S(eps) with the SAME constants (false in that form, Proposition N); anything for targets outside the 8; enclosures of the limit object's data;
  the hand proofs of sections 2-6 of THEORY_T1.md are not re-checked here beyond Lemma 2.2, which is used.

## 1. The parametrisation

Law (m, G, gamma, a_max) = (12, 1, 0.1, 1), v_s = gamma a_max^2/(4 G) = 1/40, kappa = sqrt((2 v_s - eps)/eps), lambda = pi a_max kappa/(4 G v_s).  Put
h := 1/kappa = sqrt(eps/(2 v_s - eps)), so  eps = 2 v_s h^2/(1 + h^2),  2 v_s - eps = 2 v_s/(1 + h^2),  L := 1/lambda = 4 G v_s h/(pi a_max).
The unit (Lemma 2.1) becomes:  q(x) = eps + (2 v_s - eps) tanh^2 x = 2 v_s (h^2 + tanh^2 x)/(1 + h^2);  w = (a_max/m) tanh x (eps + 2 v_s sinh^2 x);
M = (G/m) q (q + 4 v_s sinh^2 x);  Psi = a_max^2 sinh^2 x/(2 G m);  and the distance to the wall  s(x) = lambda - eta(x) = sigma_0 ( g(y)/tanh x - x ),  y = h/tanh x,  g(y) = arctan(y)/y,  sigma_0 = a_max/(2 G v_s) = 20
(Lemma 2.1: s = sigma_0 kappa (arctan(1/(kappa tanh x)) - x/kappa) = sigma_0 (arctan(y)/h - x), and arctan(y)/h = (arctan(y)/y)/tanh x).
`tests/theory3a_test.py::test_parametrisation...` compares every function with the existing `RU`/`Unit` at eps = 1e-4, 1e-9, 1e-14.

Facts used (all elementary): (i) g is decreasing on [0, inf), g(0) = 1, so s_eps(x) = sigma_0 (g(h/t)/t - x) is decreasing in h, hence in eps; (ii) q_eps(x) = eps sech^2 x + 2 v_s tanh^2 x is increasing in eps, and
q >= 2 v_s tanh^2 x for every eps (because (h^2 + t^2)/(1 + h^2) >= t^2 for t <= 1): d eta/dx = a_max/(G q) <= a_max/(2 G v_s tanh^2 x), the eps = 0 value;
(iii) w_eps(x) is increasing in eps; M_eps(x) is increasing in eps (product of positive increasing factors); M is increasing in x.

## 2. Lemma E and its proof

**2.1 Statement.**  For 0 < eps < 2 v_s and s > sigma_0:  B_eps(s) := W_eps(lambda - s) <= C / (s (s^2 - sigma_0^2)),  C = 2 v_s a_max sigma_0^3/m = 4/(G gamma^2 m) = 100/3.

*Proof.*  Lemma 2.2 (THEORY_T1.md): with theta := pi s/(2 lambda) = s/(sigma_0 kappa) (since lambda = pi a_max kappa/(4 G v_s)) one has the exact identity kappa tanh x = cot(theta + x/kappa), the argument in (0, pi/2].
(It follows from s = sigma_0 kappa (pi/2 - arctan(kappa tanh x) - x/kappa), i.e. theta = pi/2 - arctan(kappa tanh x) - x/kappa.)  Since x/kappa >= 0, cot is decreasing: tanh x <= cot(theta)/kappa =: s_eps.
Because cot(theta) <= 1/theta: s_eps <= sigma_0 / s < 1 for s > sigma_0.  w = (a_max/m) tanh x (eps + 2 v_s sinh^2 x), sinh^2 = tanh^2/(1 - tanh^2), and w increases with tanh x, so
w <= (a_max/m) s_eps (eps + 2 v_s s_eps^2/(1 - s_eps^2)).  Insert s_eps = cot(theta) sqrt(eps/(2 v_s - eps)); with A := a_max eps^(3/2)/(m sqrt(2 v_s)) and rho := 2 v_s/(2 v_s - eps) the two terms are
A rho^(1/2) cot(theta) and A rho^(3/2) cot^3(theta)/(1 - s_eps^2); their sum is <= A rho^(3/2) (1 - s_eps^2)^(-1) cot(theta)(1 + cot^2 theta) = A rho^(3/2) (1 - s_eps^2)^(-1) cos(theta)/sin^3(theta).
Adamovic-Mitrinovic (2.2): cos(theta)/sin^3(theta) <= 1/theta^3 = (sigma_0 kappa/s)^3 on (0, pi/2).  And rho^(3/2) eps^(3/2) kappa^3 = (2 v_s)^(3/2).  So
w <= (1 - s_eps^2)^(-1) (a_max (2 v_s) sigma_0^3/m) s^(-3) <= C s^(-3) (1 - sigma_0^2/s^2)^(-1) = C/(s(s^2 - sigma_0^2)).  []   (numerically: 33.333 s^-3 = 4/(G gamma^2 m) s^-3, as in Cor. 2.3.)

**2.2 Adamovic-Mitrinovic: (sin t/t)^3 > cos t on (0, pi/2).**  Put g(t) = 3 ln(sin t) - 3 ln t - ln cos t; g(0+) = 0 and g'(t) = 3 cot t + tan t - 3/t.  With cot t = 1/t - sum_{n>=1} c_n t^(2n-1), c_n = 2^(2n)|B_2n|/(2n)!,
and tan t = sum_{n>=1} d_n t^(2n-1), d_n = 2^(2n)(2^(2n) - 1)|B_2n|/(2n)!  (both series converge for |t| < pi/2 and have nonnegative coefficients), g'(t) = sum_n (d_n - 3 c_n) t^(2n-1) = sum_n 2^(2n)|B_2n|/(2n)! (2^(2n) - 4) t^(2n-1).
The coefficient vanishes for n = 1 and is positive for n >= 2, so g' > 0 on (0, pi/2), g > 0.  []

**2.3 Ball arithmetic check** (Part U of theory3a.md): for eps = 1e-2 ... 1e-24 and about 500 values of x per eps, the upper end of the ball for w(x) is below the lower end of C/(s(s^2 - 400)) at the upper end of the ball for s(x) (all points with s >= 25).
Largest ratio 0.9715 at eps = 1e-4 ... (see report).

## 3. Lemma M and Proposition N (the obstruction)

**3.1 Lemma M (PROVED).**  Fix x > 0.  (a) s_eps(x) is strictly decreasing in eps (Section 1 (i)).  (b) w_eps(x), q_eps(x), M_eps(x) are increasing in eps.  (c) x_1(eps) (the root of s_eps = delta_1) and
x_#(eps) (the root of w_eps = v_#) are decreasing in eps: at larger eps, s_eps is lower pointwise and w_eps higher, and s_eps is decreasing in x, w_eps increasing in x.  (d) lambda(eps) is decreasing in eps.
(e) The coefficient mobility at a given COEFFICIENT (not at a given x), m_V(eps) = M_eps(x_#(eps)), has no sign from (b), (c): M_eps(x) increases in eps but x_#(eps) decreases and M increases in x.  Numerically it is decreasing as eps decreases (Part S).

**3.2 Proposition N.**  pre(eps) = Psi(x_1(eps))/(v_lo - w_1(eps)), w_1(eps) = w_eps(x_1(eps)).  By (c) Psi(x_1(eps)) increases as eps decreases; the sign of w_1(eps) cannot be derived from (b),(c) (w_eps(x) increases in eps, x_1(eps) decreases);
numerically w_1 also increases as eps decreases, so pre does.  Hence the exact obstruction to "monotone in eps": the ingredients pre, T_3 = ln(.)/m_V, Env, kappa_*, and through them De_k, S_k, W_k, T^v_k, T^d_k, V_tail are all
**non-decreasing as eps -> 0** (adverse), with a relative gap to their eps = 0 limit of c*eps, c = 1.56 (pre), ~1.6e-4 eps (T_3, Env), ~4e-5 eps (kappa) on the K = 1 and K = 2 targets of Part S;
T_1, A_tr and the entrants' pre-drift are favourable (larger at larger eps; the Riemann sums used have a discretisation error of 1e-6, so their last digits carry no meaning).
Consequence: for fixed D_j, "E_j(eps_1) < D_j/1.0001" does not imply "E_j(eps) < D_j/1.0001"; but the eps = 0 value is the supremum of the adverse ones, and the ball [0, h_1] (the supremum over eps <= eps_1)
is what Theorem TAIL certifies.  The strongest TRUE monotone statement obtained is therefore Theorem TAIL (monotone majorant), not a statement about S(eps) itself.

## 4. Theorem TAIL: the verifier and the proof of its soundness

Input: the data of Theorem T1 (prep: the limit object of Thm K as floating-point data, assumed exact), the parameters pi = (r, Omega, delta_1, dhat-fraction, varsigma, d_c) of `exp_theory_t1.PARS`, the candidate bounds D_j (j non-entrants), h_hi.
Output: True iff every inequality of S(eps; pi) holds with these D_j for every h in [0, h_hi] (equivalently eps in (0, eps_hi]).  Every number below is a python-flint ball (100 bits, outward rounding); a quantity that is needed as
an upper (lower) bound is replaced by the upper (lower) end of its ball before it is compared.

1. *Box and regime (R).*  eps-free: `rig_box` of round T1b (Prop. A with the dormant budget Omega), unchanged.
2. *Layer constants over the ball.*  x_1 = root of s_h(x) = delta_1, x_# = root of w_h(x) = v_# are obtained by bisection of the monotone functions with three-valued comparisons for a ball h (so the enclosure contains every root over the ball).
   T_1 = int_{x_1}^{x_#} (d eta/dx)/(v_lo - w), A_tr, pre, T_2, T_3, Env, T_lay, kappa_*, the weights, the entrants' pre-drift: as in round T1b, evaluated by upper Riemann sums on balls of subintervals of the enclosure [lo x_1, up x_#].
   For x close to 0 (the pre-drift integral runs from 0): d eta/dx <= a_max/(2 G v_s tanh^2 x) (eps-free, valid at every x), and on [0, 1e-9] the product Psi (d eta/dx) <= a_max^3 cosh^2 x/(4 G^2 m v_s) removes the singularity.
   Nothing here needs eps > 0: y = h/tanh x is at most h_hi/tanh(x_1) ~ 1e-4 on [x_1, x_#].
3. *Dormant kernels.*  For a non-entrant j let c_j(tau) = 1 - |u_j(tau)| >= c_j^min on [0, tau_end] (the limit's dual is piecewise linear, |u| is convex on each piece, so the maximum of |u| is at a breakpoint or at the horizon),
   dist - D_j = (c - D_j L)/L =: s >= (c_min - D_j L_hi)/L_hi, required >= 26.  By Lemma E, B(dist - D_j) <= C L^3 kappa_j c(tau)^(-3), kappa_j := c_min^3/(a_lo (a_lo^2 - sigma_0^2 L_hi^2)) >= 1, a_lo := lo(c_min - D_j L)
   (env(c) c^3 is decreasing in c: c/(c - a) and c^2/((c-a)^2 - b^2) are).  For a linear piece, int c^(-3) dtau = dtau (c_a + c_b)/(2 c_a^2 c_b^2) exactly.  Hence
   mass_j = lambda int B dtau <= C L^2 kappa_j Q_j, Q_j = sum of the pieces; peak_j <= C L^3 kappa_j / c_min^3 =: L^3 e_j;  Omega' <= sum peak_j/a_min;  omega_pk Tfin <= tau_end L^2 sum e_j /a_min (all enter V_tail and (S1)).
4. *Drift.*  |c_l - c_l^0| <= mu_c ( W + peak [sum_{windows started} hd_i + 0.43 mass_tot] ) with W = sum over sub-intervals of H^(stage) x (mass in the sub-interval) -- each earlier stage with ITS H (this is the corrected accounting, (E)),
   and the dual deviation from the displacement is <= int H^(st) dcb dt <= sum over sub-intervals Delta tau H^(st) mu_c (C L kappa_j W_end + L^2 e_j (sum hd_i + 0.43 mass_tot)); the window starts are gated by their earliest possible
   start tau_i - S_i^- L_hi (a superset of the windows that have started).  No cap by the oscillation (an upper bound).  All terms are O(L) or smaller.
5. *Conditions.*  (S1) omega_pk <= Omega; (S2) max dcb <= d_c; (S3) w_1 < v_#; (S4) lambda(tau_{k+1} - tau_k) >= X iff tau_{k+1} - tau_k >= X L, checked with L_hi; (S5) peaks < theta < a_min (Gamma(dhat) - E_v - varsigma), theta < v_#; the hypotheses r > d_0, critical point in D-hat (floats).
   The strict bootstrap test E_j < D_j/1.0001.  A NaN anywhere is an error, never an accepted value (`bmax_up`, finiteness checks).
6. *Why this proves S(eps) for all eps in (0, eps_hi].*  For each fixed eps in that range the quantities of S(eps) in `exp_theory_t1.cert` are: functions of the unit at that eps (T_1 ... pre_drift, kappa_*, windows), which lie in the enclosures of step 2 because h(eps) is in the ball; the true B is
   below env (Lemma E, as dist - D >= 26 > sigma_0), the true piecewise-linear dist(t) is lambda c(tau) with c >= c_min, the mass identity and the peak are upper-bounded as in step 3, the drift by step 4; (6.1) is monotone in each of these
   quantities, so the majorants satisfy (6.1) whenever the true ones do, and the verified inequalities with the majorants imply the true ones.

Remark on the first-order scaling.  The terms of E_j that depend on the dormant side are O(L) = O(eps^(1/2)) (drift), O(L^2) (mass); all the others converge to their eps = 0 values with a relative O(eps) gap (Section 3).  Hence
D_j^infty (the sum of the layer terms at eps = 0) + the dormant terms is the proved bound, and eps_tail is determined by (i) the margin (c_min - D_j L_hi) > 26 L_hi (lambda mu_j - D_j > 26) and (ii) the sum of the dormant terms staying below
the slack D_j - (layer terms), (iii) d_c.

## 5. The finding on the verifier of round T1b (E)

In `cert` (and `rigorous.verify`) the displacement is  dcb = mu_c ( Hst[st] * I + hinc )  with I the cumulative mass up to the current time and Hst[st] the current stage's H; in Lemma 5.2 it is  mu_c ( int^t H^(st(t')) B(t') dt' + ... ).
H^(s) = sup |ghat^(s)'| is decreasing in the stage (ghat^(K) = 0), so the mass accumulated in stage s < st is weighted by H^(st) < H^(s): the legacy bound can be smaller than the lemma's.  `theory3a_chain.verify_fixed` replaces
Hst[st]*Ig by Wprev + Hst[st]*(Ig - Icum) (Wprev = sum over the earlier sub-intervals of H^(stage) x their mass), three anchored textual replacements in `rigorous.verify` (each anchor must occur exactly once).  Effect at the top of the covers (K = 2): see the
test `test_displacement_accounting...` and report Part C; for K = 1 the legacy and repaired values of dcmax coincide at the maximum (stage 1 has H ~ 0).

## 6. Gaps (explicit)

1. The limit object's data are floating-point inputs (critical points, entry times, duals, rates, margins; also the sign/assignment decisions).  Not enclosed (no Krawczyk).
2. The hand proofs of sections 2-6 of THEORY_T1.md (Nagumo, Dini, Lemma 4.1, ...) are not formalised; only Lemma 2.2 is re-derived here.
3. Only the 8 targets of Part Q (n = 1, K <= 3, continuum Gram, equal masses, noiseless).  No claim for K >= 4, n = 2 (the planar certificate is in doubles).
4. The envelope-based dormant terms are cruder than the exact ones (no cap by the oscillation, c_min instead of the time-resolved margin inside kappa_j): the tail starts one to two decades below the top of the cover; the chain covers the rest.
5. Rounding: comparisons through float(ball.upper()) as in round T1b.
6. The Riemann sums are upper sums on a fixed number of cells (1500 / 3000): rigorous as bounds, not as asymptotics.
7. H^(s) is evaluated by `rig_rgs` on a fixed subdivision (60/200 pieces) with ball bounds of g-hat'; correct as enclosures.

## 7. Files, commits, reproduction (final)

Repository (local commits, nothing pushed): df99cda (code, tests, mutant spec, report), 8f05a64 (raw).  `python experiments/theory3a_run.py --run` (about 5 min, one core) writes results/synthetic/theory3a_raw.npz and theory3a.md;
`--from-raw` regenerates the report byte for byte (checked).  `pytest tests/theory3a_test.py`: 16 tests, 15 s.  Mutants: 51 + control, 45 killed, the 6 survivors equivalent on the data (header of tools/mutants/theory3a_tail_own.py).
Numbers (Part T/C of theory3a.md): eps_tail (log10) = -8.50, -6.25, -7.25, -7.50, -8.25, -7.25, -8.50, -6.75 for the 8 targets (tops of the stored covers -7.5, -5.5, -6.5, -7.5, -8.0, -7.0, -8.5, -7.0);
3 of 8 tails reach eps_top without any chain; the other 5 chains (2 to 41 balls, 16 to 259 verify calls) close with the repaired verifier.  Max E_j/D_j 0.91 .. 0.998 (the strict test E < D/1.0001 passes).
At eps = 1e-24 the tail is re-verified separately (max E = 85.476 on target 0 against 86.514 at 1e-12 and 92-93 near 1e-7: the deviation bound converges to its eps = 0 value 85.476).
