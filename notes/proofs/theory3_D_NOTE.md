# NOTE (theory3, direction D): which readouts select -- theorems for the readouts of depth p

Status vocabulary: PROVED (full argument below, each step checked, numerically sanity-tested), COMPUTER-ASSISTED (what is computed, in what arithmetic), OUTLINE, FAILED.
Code: `experiments/theory3d_unit.py` (closed forms, ball arithmetic), `theory3d_scan.py` (per-target scan), `theory3d_readouts.py` (parts A-H, report `results/synthetic/theory3d.md`).
Conventions are those of the repository (`src/core/factorized_rbf.py`): all parameters of a unit follow the gradient flow of F = 1/2 <r^2>, r = y - sum_i w_i phi(.;c_i), with the learning rate ell
(ell = 1 for D2..D5, ell = G m = 12 for D1), from u(0) = v(0) = s, s^2 = eps b (b > 0 the "mass"; b = 1 equal masses).

## 0. What is claimed (read first)

Let a *unit* have the coefficient w and the parameters (u, v) (depth >= 2) or u (depth 1). Its *drive* is g(t) = -dF/dw = <r phi(.;c(t))>; its *dual* is eta(t) = int_0^t g. Whatever the other
units, their coefficients and the unit's own centre do, the unit sees the network only through g.

* **Theorem 1 (the unit, exactly).** depth 1: w = ell eta. depth 2: w = 2 eps b sinh(4 ell eta). depth p >= 3: w = s^p [(1-eta/lam)^-q - (1+eta/lam)^-q], q = p/(p-2),
  lam = s^-(p-2)/(ell p (p-2)); the solution exists exactly while |eta| < lam (a wall). PROVED (section 1).
* **Theorem 2 (depth p >= 3: dormant units are polynomially small).** If |eta| <= (1-mu) lam on [0,T] then |w| <= s^p [mu^-q - (2-mu)^-q] <= (eps b)^(p/2) mu^(-p/(p-2)), sharp. In the distance d = lam - |eta|:
  |w| = (ell p (p-2) d)^-q [1 - (x/(2-x))^q] (x = d/lam): eps-free to leading order. The unit is formed (|w| > theta) iff d < d_theta, with theta^(-(p-2)/p)/(ell p (p-2)) (1 - O(eps^(p/2)...)) <= d_theta <= theta^(-(p-2)/p)/(ell p (p-2));
  in margin units mu_theta = d_theta/lam = theta^(-(p-2)/p) (eps b)^((p-2)/2) (1 - tiny). The centre of a dormant unit with drive >= g_min moves by at most
  ell_c H_max eps b/(2 ell p g_min) mu^(-2/(p-2)). PROVED (section 2). For p = 3: eps^(3/2) mu^-3 (the paper), formation margin theta^(-1/3) eps^(1/2); p = 4: eps^2 mu^-2, margin theta^(-1/2) eps; p = 5: eps^(5/2) mu^(-5/3), margin theta^(-3/5) eps^(3/2).
  The paper's own unit has the same form with other constants (Lemma 2.2 of the T1 note: PROVED there; here: the same chain re-verified numerically, Part A).
* **Theorem 3 (depth 2: exactly logarithmic).** w = 2 eps b sinh(4 ell eta); with the wall lam_2 = ln(1/(eps b))/(4 ell): sup |w| over |eta| <= (1-mu) lam_2 equals (eps b)^mu - (eps b)^(2-mu) (so "eps^mu"); the unit is
  formed iff d < d_theta = ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/(4 ell), i.e. iff mu < mu_theta = ln(2/(theta+sqrt(theta^2+4(eps b)^2)))/ln(1/(eps b)), and ln(1/theta)/ln(1/(eps b)) - delta <= mu_theta <= ln(1/theta)/ln(1/(eps b)),
  delta = -ln(1 - (eps b)^(2-2 mu_0))/ln(1/(eps b)) (negligible). A depth-2 network keeps a unit dormant at the (actual) margin mu only if ln(1/eps) > ln(1/theta)/mu: PROVED (section 3); the rule of the paper is this statement with the path's margin in place of the actual one.
  Also: the dual at which the unit reaches the amplitude A is lam_2 (1 + ln(A/b)/ln(1/eps)): the depth-2 walls depend on the amplitudes at the relative order 1/ln(1/eps) (a non-uniformity; PROVED, the consequences for the path are not tested).
* **Theorem 4 (depth 1: no selection).** w = ell eta exactly. For fixed centres with independent columns: g(t) = exp(-ell H t) g(0), w(t) = ell int g -> the least-squares solution, a gradient flow of the quadratic loss (the ridge-like
  path), and |w_j(t) - ell t g_j(0)| <= ell^2 ||H|| ||g(0)|| t^2/2; every unit with g_j(0)^2 > 2 theta ||H|| ||g(0)|| is formed by t_j = 2 theta/(ell |g_j(0)|), independent of eps. PROVED (section 4).
* **Corollary 5 (limit objects).** For fixed centres and ell, b fixed, the mirror flows of p >= 3 and of the paper's unit converge to the inverse scale space of sum_j b_j^(-(p-2)/2) |w_j| (weights b^-1/2 for p = 3); depth 2
  converges to the unweighted l1 inverse scale space, but only logarithmically (the walls are amplitude dependent at the relative order 1/ln(1/eps); on the design the dual error falls like a power of 1/ln(1/eps), L); depth 1 is the l2 (ridge-like) flow, never sparse. The convergence statement reduces to the same Mosco/Attouch/Brezis-Pazy argument
  as the paper's fixed-centre theorem; what is new is the verification of epi-convergence for D_p and D2 (PROVED, section 5). The sliding limit (moving centres) is not proved for p != the paper's unit (OUTLINE): its dormant part is Theorem 2.
* **Proposition 6 (the sliding unit; the mismatch margin).** Frozen-environment unit: eta' = G(c), c' = mu_c W(eta) G'(c) has the first integral Psi*(eta) - Theta(c), Theta' = G/(mu_c G'); for a linear G: g^2 = g0^2 + 2 kappa Psi*(eta) (PROVED, section 6).
  This gives the competitor's lead and an overtaking margin sqrt(kappa a_ro eps)/|g0| (a model; COMPUTER-ASSISTED comparison with the 200 targets, section 7), which explains why the mismatch margin is not ~ eps^(1/2) with a universal constant and how the observed
  pooled exponent ~ 0.3 (paper: 0.26) arises.

## 1. Theorem 1: the exact unit

Setting. The network is a system of ODEs in (u_i, v_i, c_i); unit i's parameter equations are u' = ell dw/du g, v' = ell dw/dv g with g(t) = <r(t) phi(.;c_i(t))> continuous along any solution (F depends on unit i only through w_i and c_i;
-dF/dw_i = <r phi_i> = g_i). So for the unit, g is a given continuous function of t on the interval of existence of the solution; nothing else enters.

*Depth p >= 3.* u' = ell p u^(p-1) g, v' = -ell p v^(p-1) g, u(0) = v(0) = s > 0. u and v stay positive: u = 0 is an equilibrium of the (locally Lipschitz in u) equation, so a solution that starts at s > 0 cannot reach it.
Put A = u^-(p-2), B = v^-(p-2): A' = -(p-2) u^-(p-1) u' = -ell p (p-2) g, B' = +ell p (p-2) g, so
A = s^-(p-2) - ell p (p-2) eta = s^-(p-2) (1 - eta/lam), B = s^-(p-2) (1 + eta/lam), lam := s^-(p-2)/(ell p (p-2)).
u, v finite and positive means A, B > 0, i.e. |eta| < lam; conversely, while |eta| < lam the formulas are bounded, so the (u, v)-part cannot blow up before |eta| reaches lam (the full system continues as long as the other units and the centres stay finite). Hence u = s (1-eta/lam)^(-1/(p-2)), v = s (1+eta/lam)^(-1/(p-2)), and
w = u^p - v^p = s^p [(1-eta/lam)^-q - (1+eta/lam)^-q], q = p/(p-2). As eta -> lam, u -> infinity: |w| -> infinity (the second term is bounded by 1 times s^p). The wall is lam; lam = s^-(p-2)/(ell p (p-2)) = (eps b)^(-(p-2)/2)/(ell p (p-2)); the
repository's `PowerN.scale` and `Cubic3.scale = 1/(3 ell sqrt(eps))` are this lam (b = 1).
W_p(eta) := that function is odd, W_p' = (s^p q/lam)[(1-sg)^(-q-1) + (1+sg)^(-q-1)] > 0 (sg = eta/lam), W_p -> +-infinity at +-lam: strictly increasing and onto R.
Mass identity: Psi_p(eta) := int_0^eta W_p = s^p lam/(q-1) [(1-sg)^-(q-1) + (1+sg)^-(q-1) - 2], q - 1 = 2/(p-2) (differentiate).
Distance form. Put d = lam - eta, x = d/lam in (0,1]. Then 1 - sg = x, 1 + sg = 2 - x, and s^p x^-q = s^p lam^q d^-q with lam^q = s^-p (ell p (p-2))^-q, so
W_p = (ell p (p-2) d)^-q [1 - (x/(2-x))^q]. (D)
*Depth 2.* u' = 2 ell u g, v' = -2 ell v g: u = s e^(2 ell eta), v = s e^(-2 ell eta), w = s^2 (e^(4 ell eta) - e^(-4 ell eta)) = 2 eps b sinh(4 ell eta). Psi_2 = eps b (cosh(4 ell eta) - 1)/(2 ell).
*Depth 1.* w = u, u' = ell g, w(0) = 0 (`Linear`): w = ell eta.
Numerical checks: Part A (ODE vs closed form, errors <= 1e-10 relative), Part B (the derivative identities du/deta = ell p u^(p-1) etc. in arb_series), Part C (the repository's network flows, moving centres: |eta_integrated - W^-1(w)|/wall <= 1e-8 for D2..D5 and the paper's unit). QED.

Remark (the paper's unit). For w = nu^2 a_max tanh(alpha/a_max)/m with the rates (Gm, gamma m/4), Lemma curve (ii) of the supplement gives eta' = g exactly and w = W^eps(eta), with wall lam = pi a_max kappa/(4 G v_s) ~ 7.0248 eps^(-1/2); it behaves as depth 3
because near its critical curve w ~ u^3. THEORY_T1 Lemma 2.2 gives the dormant bound with the explicit constant.

## 2. Theorem 2 (depth p >= 3)

**(i) Dormant bound.** Let T > 0, mu in (0,1], |eta(t)| <= (1-mu) lam for t <= T. W_p is odd and increasing, so |w(t)| <= W_p((1-mu) lam) = s^p [mu^-q - (2-mu)^-q] <= s^p mu^-q = (eps b)^(p/2) mu^(-p/(p-2)).
Equality in the first inequality iff |eta| attains (1-mu) lam. QED. (Ball arithmetic B1; Part C: ratio max_t|w|/exact-bound = 1 to 1e-4.)
For p = 3, b = 1: eps^(3/2) mu^-3; the paper's Thm (a) has C eps^(3/2) mu^-3 with C = 0.373 (the supplement's constant) for its unit.

**(i') In terms of the drive.** If sup_[0,T] |g| <= gamma (a dormant unit's drive is bounded by its path rate), then |eta(t)| <= gamma t, so the hypothesis of (i) holds up to the time T_mu = (1-mu) lam/gamma = (1-mu) lam/gamma: for t <= T_mu, |w(t)| <= (eps b)^(p/2) mu^(-p/(p-2)).
In the path's time tau = t/Lambda (Lambda = lam for b = 1) the unit stays at the margin mu until tau = (1-mu)/gamma: exactly the statement 'the coefficient stays below C eps^(p/2) mu^(-p/(p-2)) up to the time the dual reaches its wall'. (With the actual drive replaced by the path's, the hypothesis is the dual-deviation hypothesis of 7.1.)

**(ii) Distance form and universality.** By (D), for d = lam - |eta| the coefficient is (ell p (p-2) d)^-q times a factor in (1 - O((d/lam)^q), 1]. Hence |w| <= (ell p (p-2) d)^-q whatever eps: the *absolute* distance of the dual to the wall decides, and eps enters only
through the wall height lam ~ eps^-(p-2)/2. The formation distance: F(d) := W_p(lam-d) is strictly decreasing in d, F(0+) = infinity, F(lam) = 0, so |w| > theta iff d < d_theta (unique root). With d_0 := theta^(-1/q)/(ell p (p-2)) (solves (ell p (p-2) d)^-q = theta):
F(d_0) = theta (1 - r(x_0)^q) < theta, r(x) = x/(2-x), x_0 = d_0/lam < 1 (true for eps small), so d_theta < d_0; and with d_lo := d_0 (1 - r(x_0)^q)^(1/q): F(d_lo) = theta (1 - r(x_lo)^q)/(1 - r(x_0)^q) >= theta (x_lo <= x_0),
so d_lo <= d_theta. In margin units d_0/lam = theta^(-(p-2)/p) (eps b)^((p-2)/2) (ell cancels: d_0/lam = theta^(-1/q) s^(p-2)), and 1 - d_lo/d_0 = 1 - (1 - r(x_0)^q)^(1/q) <= r(x_0)^q ~ (x_0/2)^q (since (1-y)^(1/q) >= 1-y for y in [0,1], q >= 1): negligible
(< 1e-8 for eps <= 1e-6, p = 3, theta = 0.01; 5e-5 at eps = 1e-4, theta = 1e-3).
**Consequence ("what selects")**: if the actual dual of a dormant unit stays at margin >= mu_theta(1+tiny) the unit is below theta forever; if it comes within d_theta it is formed. mu_theta ~ theta^(-(p-2)/p) eps^((p-2)/2) b^((p-2)/2): 4.64 eps^(1/2) (p = 3, theta = 0.01),
10 eps (p = 4), 15.8 eps^(3/2) (p = 5). The polynomial rate in eps of the selection margin is (p-2)/2: deeper readouts select at smaller margins. (B5 in ball arithmetic.)

**(iii) Drift of a dormant centre.** The centre obeys c' = ell_c w <r d_c phi> (the repository: ell_c = G m = 12). If g has a constant sign with |g| >= g_min on [0,T] and |<r d_c phi>| <= H_max, then
|c(T) - c(0)| <= ell_c H_max int_0^T |w| dt = ell_c H_max int_0^{|eta(T)|} W(eta)/|g| d eta <= ell_c H_max Psi_p(|eta(T)|)/g_min <= ell_c H_max (eps b)/(2 ell p g_min) mu^(-2/(p-2)).
Last step: Psi_p((1-mu) lam) = s^p lam/(q-1) [(1-sg)^-(q-1) + (1+sg)^-(q-1) - 2] with (1+sg)^-(q-1) < 1, so the bracket is < mu^-(q-1); and s^p lam/(q-1) = (eps b)^(p/2) (eps b)^(-(p-2)/2)/(ell p (p-2)) (p-2)/2 = eps b/(2 ell p). (Ball arithmetic B4.)
For p = 3: eps mu^-2 (the paper's "O(eps)"), p = 4: eps mu^-1, p = 5: eps mu^(-2/3). **Test on the stored runs of the design**: the largest dormant displacement divided by Psi_p((1-mu) lam) at the design margin is 245 +- 3% over D3, D4, D5 at eps = 1e-6, 1e-8, 1e-10 (Part F): a displacement proportional to the mass identity,
for three depths and three eps, as the first integral of section 6 says. The bound needs the *lower* bound g_min: section 7 shows that a dormant unit with a small drive (a nearly stalled dual) is the weak point.

## 3. Theorem 3 (depth 2)

From Theorem 1: w = 2 eps b sinh(4 ell eta). Define lam_2 by e^(4 ell lam_2) = 1/(eps b) (the repository's `Product2.scale = ln(1/eps)/(4 ell)` for b = 1: the exact wall where |w| = 1 is lam_2 + O((eps b)^2)).
(i) sup over |eta| <= (1-mu) lam_2 of |w| = 2 eps b sinh((1-mu) ln(1/(eps b))) = (eps b)^mu - (eps b)^(2-mu) =: f(mu), strictly decreasing from f(0)=1-(eps b)^2 to f(1) = 0 (f' = ln(eps b)(...) < 0), attained.
(ii) distance form: with d = lam_2 - |eta|: |w| = e^(-4 ell d) - (eps b)^2 e^(4 ell d). |w| >= theta iff y := e^(-4 ell d) satisfies y - (eps b)^2/y >= theta iff y >= (theta + sqrt(theta^2 + 4 (eps b)^2))/2: so exactly
d_theta = ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/(4 ell), mu_theta = d_theta/lam_2 = ln(2/(theta + sqrt(theta^2+4(eps b)^2)))/ln(1/(eps b)).
(iii) The two descriptions agree: with d = mu lam_2, e^(-4 ell d) = e^(-mu L) and (eps b)^2 e^(4 ell d) = e^(-(2-mu) L), L = ln(1/(eps b)), so |w| = f(mu) := e^(-mu L) - e^(-(2-mu) L) and mu_theta is the root of f(mu) = theta. f is strictly decreasing on [0,1]
(f' = -L (e^(-mu L) + e^(-(2-mu)L)) < 0). With mu_0 = ln(1/theta)/L (< 1 iff eps b < theta): f(mu_0) = theta - e^(-(2-mu_0) L) < theta, so mu_theta < mu_0. With delta := -ln(1 - e^(-(2-2 mu_0) L))/L >= 0 (so e^(delta L) = 1/(1 - e^(-(2-2mu_0)L))):
f(mu_0 - delta) = e^(-(mu_0-delta) L) (1 - e^(-(2 - 2 mu_0 + 2 delta) L)) >= theta e^(delta L) (1 - e^(-(2-2mu_0) L)) = theta, so mu_theta >= mu_0 - delta. Hence mu_0 - delta <= mu_theta <= mu_0 with delta ~ e^(-(2-2mu_0)L)/L (< 1e-8 for eps <= 1e-6, theta = 0.01); ball arithmetic B6.
Consequence: a depth-2 unit whose actual dual stays at the margin mu of lam_2 is not formed iff mu > mu_theta ~ ln(1/theta)/ln(1/eps): "exp(-mu ln(1/eps))". For the margin mu of a target to be admissible: ln(1/eps) > ln(1/theta)/mu,
e.g. design: mu = 0.0664, theta = 0.009: eps < 1e-30.8 (Part D of the stored runs: D2 fails at 1e-24, passes at 1e-48).
(iv) Amplitude dependence of the wall: |w| = A at eta_A with 2 eps b sinh(4 ell eta_A) = A: eta_A = lam_2 (1 + ln(A/b)/ln(1/eps)) + O((eps b/A)^2): the entry of a depth-2 unit onto the amplitude A happens at the
relative dual 1 + ln(A/b)/ln(1/eps) (in units of the wall of b = 1): at eps = 1e-6, amplitudes spanning 10^3 move walls by +- 50%. (Not tested here: the path of D2 at finite eps is a weighted-l1 path with these walls.)

## 4. Theorem 4 (depth 1)

(a) w = ell eta (Theorem 1) for any centre motion: there is no wall, no eps, no mass; |w| > theta iff |eta| > theta/ell.
(b) Fixed centres: Phi (N x d) with independent columns, H = Phi^T Phi/N > 0, b0 = Phi^T y/N, g = (g_j) = Phi^T r/N = b0 - H w. Then g' = -H w' = -ell H g, so g(t) = e^(-ell H t) g(0), w(t) = ell int_0^t g = H^-1 (I - e^(-ell H t)) b0 -> H^-1 b0 = w_LS as t -> infinity (all eigenvalues of H positive).
The flow is the gradient flow of 1/2 N^-1 |y - Phi w|^2 from 0: for t finite it is the ridge-type regularisation path (Tikhonov with parameter ~ 1/(ell t) in the spectral sense), for t = infinity the LS fit; neither is sparse. Generic density: the set of y with (w_LS)_j = 0 is a hyperplane.
(c) Formation. ||e^(-A) - I|| <= ||A|| for A >= 0 gives |w_j(t) - ell t g_j(0)| <= ell int_0^t ell s ||H|| ||g(0)|| ds = ell^2 ||H|| ||g(0)|| t^2/2. Hence for t <= T_j := |g_j(0)|/(ell ||H|| ||g(0)||): |w_j(t)| >= ell t |g_j(0)|/2.
At t_j = 2 theta/(ell |g_j(0)|) this is theta, and t_j <= T_j iff g_j(0)^2 >= 2 theta ||H|| ||g(0)||. So every such unit is formed by t_j, a time with no eps in it (while the first depth-p entry takes the time of order lam ~ eps^-(p-2)/2 or ln(1/eps)).
(d) Moving centres: the depth-1 centre moves with speed ell_c w <r d_c phi>, proportionally to a coefficient that is not small, so all units whose drive is not tiny also slide: no dormancy to protect them (the stored runs of Part B: 12 units formed on the design at every eps from 1e-2 to 1e-10).
Checks (Part D): flow vs matrix exponential 1e-10; sigma_min(Phi) > 0; LS solution has 19 of 21 coefficients above the formed threshold; formation times vs the certificate.

## 5. Corollary 5 (limit objects) -- proof

Fixed centres. Write tau = t/Lambda, u_j = eta_j/Lambda, Lambda := lam(eps, b = 1) (p >= 3) or lam_2(eps, b=1). Then u_j' = <r phi_j> and w_j = W_j(Lambda u_j) with wall rho_j = lam_j/Lambda = b_j^(-(p-2)/2) (p >= 3);
rho_j = ln(1/(eps b_j))/ln(1/eps) -> 1 (p = 2). As in the supplement (Thm iss), u = Phi^T zeta and zeta' = -grad D_eps(zeta)/N with D_eps(zeta) = sum_j g_{eps,j}((Phi^T zeta)_j) - y^T zeta/N... where g_{eps,j}(v) = Psi_j(Lambda v)/Lambda is convex, even, >= 0, with
 p >= 3: g_{eps,j}(v) = b_j eps^(p/2)/(q-1) [(1-v/rho_j)^-(q-1) + (1+v/rho_j)^-(q-1) - 2] for |v| < rho_j and +infinity otherwise (barrier);
 p = 2:  g_{eps,j}(v) = (2 eps b_j/L)(cosh(L v) - 1), L = ln(1/eps).
(i) p >= 3: for |v| <= (1-mu) rho_j: 0 <= g <= b_j eps^(p/2) mu^-(q-1)/(q-1) -> 0 uniformly on compacts of (-rho_j, rho_j); for |v| >= rho_j g = +infinity. So g_eps -> chi_{[-rho_j, rho_j]} (the indicator): lim inf: if v_eps -> v with |v| > rho_j then
g(v_eps) = +infinity eventually; if |v| <= rho_j then g >= 0 = chi(v). Recovery: v_eps = (1 - delta_eps) v with delta_eps -> 0 slowly. Hence epi-convergence (and Mosco convergence in finite dimension).
(ii) p = 2: g(v) = (2 eps b_j/L)(cosh(L v) - 1) <= (eps b_j/L) eps^-|v| = b_j eps^(1-|v|)/L -> 0 for |v| < 1 (uniformly on compacts), -> +infinity for |v| > 1 (>= b_j eps^(1-|v|)/(2L)(1 - ...)) and g(+-1) <= b_j/L -> 0. So g_eps -> chi_{[-1,1]} again, epi-convergence, with the rate in eps only logarithmic (g ~ eps^(1-|v|), all margins
smaller than ~ ln(1/theta)/L matter).
(iii) The rest of the argument (Attouch's theorem: graphical convergence of subdifferentials; Brezis-Pazy: convergence of the gradient flows from zeta(0) = 0 in the interior of the limit's domain; the energy identity for the strong convergence of Phi w) is that of the supplement's Theorem (iss) verbatim,
since it uses only (a) Mosco convergence of D_eps to chi_P - y^T zeta/N, P = {|(Phi^T zeta)_j| <= rho_j}, and (b) D_eps(0) = 0. Hence the limit is the inverse scale space flow of R(w) = sum_j rho_j |w_j|. Status: as the paper's Theorem (iss): the abstract steps are cited, the verification of their hypotheses for D_p and D2 is proved here.
(iv) p = 1: w = ell eta, Psi_1(eta) = ell eta^2/2: the mirror map is quadratic, R(w) = w^2/(2 ell), the flow is the gradient flow of the least squares (the 'inverse scale space' of the Euclidean norm: the Bregman distance is the squared distance): no wall, no thresholding, no l1.
Moving centres: the dormant part is Theorem 2 for every p >= 3; the layer analysis (H3), the tracking and the bootstrap of the T1 note were carried out for the paper's unit only; for D_p the numerical evidence is Part B/C of the stored runs (the formed set is the path's on 166 to 195 of 200 targets for P0 and D3 at 1e-6, 1e-8; D4, D5 on the design only). OUTLINE.

## 6. Proposition 6 (the sliding unit)

Model: a dormant unit in a frozen environment: the residual r is fixed, g(c) := <r phi(.;c)>, so g_c := dg/dc = <r d_c phi>. With c' = mu_c w g_c (mu_c = ell_c) and w = W(eta), eta' = g(c):
 eta' = G(c), c' = mu_c W(eta) G'(c).
(1) First integral. On an interval where G' != 0 define Theta(c) = int_{c0}^c G/(mu_c G') dz. Then d/dt [Psi*(eta) - Theta(c)] = W G - (G/(mu_c G')) mu_c W G' = 0. QED. (Psi* = int W, Theorem 1; any readout.) So the displacement of the centre is a function of the dual alone.
 Linear G = g0 + h (c - c0): Theta - Theta(c0) = ((g0 + h dc)^2 - g0^2)/(2 mu_c h^2), so g^2 = g0^2 + 2 kappa Psi*(eta), kappa = mu_c h^2 (with Psi*(eta) measured from the start of the run). [Part E: first-integral error <= 1e-10 relative.]
(2) Lead. To first order in kappa, g - g0 = kappa Psi*(eta)/g0, and the dual gains, over the frozen dual eta_0(t) = g0 t, eta - eta_0 = (kappa/g0^2) int_0^{eta_0} Psi* d eta' = (kappa/g0^2) lam (eps/(6 ell)) [1/(1-sg) - 1/(1+sg) - 2 sg] (depth 3; Psi*_3 from Theorem 1),
 sg = eta_0/lam; at margin x = 1 - sg << 1 the bracket is 1/x, and accrued while the margin goes from x_a to x: kappa a_ro eps/g0^2 (1/x - 1/x_a), a_ro = 1/(6 ell) (depth 3), a_ro = 2/(G gamma^2 m c_u^2) = 0.3377 (paper's unit, from B_0 s^3 -> 4/(G gamma^2 m)).
 [Part E: agrees with the ODE to 0.02% (kappa/g0^2 = 1), 2% (108); fails by 40% at kappa/g0^2 = 1200 where the displacement is 3 widths.]
(3) Overtaking margin. A competitor j with margin x_a at the start and x at the end of the stage in which the entrant arrives (the entrant has margin 0 there) has gained the lead above; it overtakes the entrant iff lead >= x. With K_c := kappa a_ro eps/g0^2, the root of x = K_c (1/x - 1/x_a) is
 x* = (-K_c/x_a + sqrt(K_c^2/x_a^2 + 4 K_c))/2 = sqrt(K_c)(1 - O(sqrt(K_c)/x_a)): so x* ~ sqrt(mu_c a_ro eps) |h|/|g0|. The competitor flips the order (or forms) iff x < x*: the margin below which a competitor intrudes is ~ c_j eps^(1/2) with c_j = sqrt(mu_c a_ro) |h_j/g_j|, the *logarithmic derivative* of its drive landscape:
 large when a dormant unit has a nearly vanishing drive on a steep flank. No fitted parameter (mu_c = 12; a_ro as above). [Overtaking and the first integral are for a competitor that is never captured; capture by the feature is the end of the model.]

## 7. The mismatch margin of the 200 targets and the exponent 0.26 (numbers: results/synthetic/theory3d.md, Parts F, G, H, K, L)

**7.1 What the scalar theorem predicts.** Hypothesis of Theorem 2: the *actual* dual of the unit. In path terms: if |eta_j|/lam <= u_j^path + delta on the whole run (delta the dual deviation, = eps^(1/2) dev_j in the paper's Thm 1(d)) then the actual margin is >= mu - delta and
|w_j| <= (eps b)^(p/2) (mu - delta)^(-q). Formation: a never-entering unit is formed iff its actual margin is < mu_theta = d_theta/lam. For depth 3 (b = 1): mu_theta = theta^(-1/3) eps^(1/2) = 4.64e-3 (eps = 1e-6) and 4.64e-4 (1e-8) for theta = 0.01;
with the per-target formed threshold theta_t = min(0.01 max|A|, 0.1 min|A|) the range over the 200 targets is 0.0022-0.030 at 1e-6 (small features have small theta_t). So the theorem predicts (dual = path dual): D3 mismatch iff margin_never < mu_theta,t(eps): 18 targets at 1e-6 (observed
set mismatches: 24; 12 of the 18 are among them) and 5 at 1e-8 (observed 5; 3 among them). For the paper's unit the same statement is empty: its coefficient at the wall is already 0.0079 < theta = 0.009 (a P0 unit is formed only after its dual has passed the wall: mu_theta < 0 for theta >= 0.0079; positive only for
the 2 targets with theta_t < 0.0079): all of P0's mismatches (34 and 9) are the effect of the dual deviation of the *actual* flow from the path. Hence the theorem cannot, by itself, give the paper's mismatch margins: they are set by the deviation delta, i.e. by the overtaking mechanism below.

**7.2 Mechanism: competitor drift (Proposition 6) and the overtaking margin.** A dormant competitor j in a stage with drive g (positive) and slope h = <r d_c phi> has coefficient w = W(eta) > 0 hence centre speed mu_c w h: its drive grows, dg/dt = mu_c w h^2, with the exact first integral g^2 = g0^2 + 2 kappa Psi*(eta).
The dual gains a lead kappa a_ro eps/g0^2 (1/x - 1/x_a); the competitor intrudes if its lead exceeds its margin x at the entry of the stage's entrant, i.e. if x < x* ~ c_j eps^(1/2), c_j = sqrt(mu_c a_ro) |h/g|. Applied to every (stage, competitor) of the sliding path of each target (no fitted parameter; mu_c = 12; a_ro = 1/6 for D3, 0.3377 for the paper's unit):
 * D3, eps = 1e-6: 40 targets predicted, 24 of the 25 flips caught, 16 false alarms; eps = 1e-8: 5 of 6 caught, 5 false alarms. Paper's unit: 31 of 36 (1e-6) and 9 of 10 (1e-8) caught, 15 and 4 false alarms. The model is conservative (it ignores the entrant's own lead): its critical eps is 2.2 (D3) to 3.1 (paper) times the measured one on the tie family (K), and 1.4 (median, D3) over the scanned targets.
 * the largest margin of a predicted mismatch: D3 0.0185 / 0.0030 at 1e-6 / 1e-8 (observed: 0.0185 / 0.0030, identical); paper's unit 0.0226 / 0.0052 (observed 0.0185 / 0.0057).

 * per-target critical eps (scan, D3, 91 bracketed targets): model and observed critical eps have Spearman correlation 0.77, log-log slope 0.87 and a median ratio 10^0.15 = 1.4 (G); paper's unit (25 bracketed): Spearman 0.49.
**7.3 The exponent.** Two statements.
 (i) At fixed geometry the law is eps^(1/2): along a one-parameter family of ties (K; one feature between two kernels, margin varied from 1e-3 to 0.13) the margin versus the bisected critical eps has slope 0.537 (D3) and 0.533 (paper) -- the sqrt law with small corrections;
 the model gives 0.55 / 0.54. The fixed-centre duals of the network converge to the ISS duals at sqrt(eps) (L: slope 0.509, 0.500) and, for depth 2, at a logarithmic rate.
 (ii) The pooled exponents have a simple origin: heterogeneity of the constant. Per (target, competitor) margin = c eps*^(1/2) (i). Across targets, log10 margin = log10 c + (1/2) log10 eps*, so the regression slope of log margin on log eps* is
 (Var a - Cov(a,b))/(2 (Var a + Var b - 2 Cov(a,b))), a = log10 margin, b = log10 c; for log c independent of the margin it is (1/2) Var a/(Var a + Var b): an errors-in-variables attenuation of 1/2. Measured on the scan (D3, 91 targets with a bracketed critical eps, margins 3e-4 to 0.06): Var a = 0.405, Var b = 0.425 (c varies by a factor 4.5 per standard deviation: c = sqrt(mu_c a_ro)|h/g| is the inverse logarithmic-derivative scale of the competitor's drive: median 3.8, 90th percentile 49, 99th percentile 95 over the competitor pairs near their wall),
 correlation 0.45, the slope is 0.239 -- the paper's 0.26 -- and the independent-c formula gives 0.244 (Spearman correlation of model and observed critical eps 0.77, log-log slope 0.87). For the paper's unit the scan is smaller (55 targets of margin < 0.01, a grid of five decades, 25 bracketed): Spearman 0.49 and a slope of 0.10 in the truncated margin range (correlation of c with the margin 0.81): not a test of the independent formula. The same fact in the form of the sample's counts:
 N(eps) = n F_rho(sqrt(a eps)), F_rho the distribution function of rho_t = min over (stage, competitor) of x/sqrt(mu_c a_ro)(|g/h|); for a product of two variables with positive density at 0 the distribution function is r(1 + ln(1/r)), local index 1 - 1/(1 + ln(1/r)) (0.9 at r = 1e-4);
 the measured local index of F_rho on the 1034-target population is 0.72-0.76 on [1e-5, 1e-3] (heavier than the product law because of the tail of c), so N(eps) falls as eps^0.36-0.38 locally, as eps^0.31-0.33 between 1e-6 and 1e-8 (model, population) and 0.29-0.34 observed on the 200 targets (Poisson uncertainty +-0.1).
 The largest margin of the sample is an extreme-value statistic: model bootstrap of 200-target samples: exponent 0.38 (D3; 5-95%: 0.09-0.71) and 0.34 (paper; 0.15-0.55), observed 0.40 and 0.26. So "0.26, not 1/2" is neither a law nor an anomaly: each (target, competitor) obeys the sqrt law with its own constant c_j (spread over two orders of magnitude), and the sample statistic mixes them.
 (A candidate explanation through the exponents in mu -- the dormant coefficient eps^(3/2) mu^-3 times a time of order 1/mu would suggest mu^4 ~ eps, i.e. eps^(1/4) -- does not arise in the first-integral analysis: the accumulated mass is Psi* ~ eps mu^-2 (not eps^(3/2) mu^-3 x time), so a unit with margin mu needs mu^2 ~ kappa a_ro eps/g0^2, i.e. mu ~ eps^(1/2) with the constant above, as the tie family confirms.)
This explanation is COMPUTER-ASSISTED (a first-order model with no fitted parameter, checked against 1200 stored runs and the scan), not a theorem.

## 8. Wording the paper may use

* "For a unit of depth p >= 3 trained by gradient flow from s^2 = eps b, w = W_p(eta) exactly, eta = int g (Lemma); at the fraction mu of the wall |w| <= (eps b)^(p/2) mu^(-p/(p-2)), and the unit is formed iff its dual is within d_theta ~ theta^(-(p-2)/p)/(ell p (p-2)) of the wall, i.e. at a margin below theta^(-(p-2)/p) (eps b)^((p-2)/2)." (PROVED)
* "Depth 2: w = 2 eps b sinh(4 ell eta); a unit at the margin mu of the wall carries (eps b)^mu - (eps b)^(2-mu), and is formed iff mu < ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/ln(1/(eps b)): the logarithmic scale is exact for the scalar unit." (PROVED)
* "Depth 1: w = ell eta, no wall: in fixed centres the flow is the least-squares gradient flow (dense), every unit with initial drive above ... is formed by t = 2 theta/(ell |g|), at no cost in eps." (PROVED)
* "The fixed-centre limit is the weighted l1 inverse scale space for p >= 3 (weights b^-(p-2)/2), unweighted l1 at a logarithmic rate for p = 2, the l2 flow for p = 1." (PROVED up to the cited Attouch/Brezis-Pazy steps; verified in Part L.)
* "The mismatch margin of a competitor is c eps^(1/2) with c = sqrt(mu_c a_ro)|h/g|; the pooled 0.26 is the local exponent of a sample statistic." (model, COMPUTER-ASSISTED)
Must NOT claim: that the depth-p sliding limit is proved for p != the paper's unit; that the scalar bound alone predicts the paper's mismatch margins (it explains 12 of 24 D3 mismatches and 2 of 34 for the paper's unit); that the overtaking model is exact (60% precision); that the depth-2 rule's accuracy (124 of 126) is a theorem (the exact margin of the scalar unit is; the path's margin enters through the dual deviation, of relative order 1/ln(1/eps)).

## 9. Gaps

1. The hypotheses of Theorems 2 and 3 are on the actual dual, not on the path's. Bridging them is the tracking estimate (T1 note) for the paper's unit only; for D_p it is numerical (Parts C, F, L).
2. The sliding limit with moving centres for D_p (p != the paper's unit): layer analysis (H3) and bootstrap not done. OUTLINE.
3. The overtaking model: first-order in the drift, frozen residual, entrant's own lead ignored; its precision is 60% and it overpredicts the critical eps by 1.7-3x; the exponent explanation (7.3) inherits this. The product-form heuristic for the local index (1 - 1/ln(1/r)) is a heuristic.
4. Corollary 5 for D2: the rate of the dual error is only bounded by the walls' amplitude dependence (relative ln A/ln(1/eps)); the numerics (L) show a faster decay on the design, not explained.
5. Variants of the overtaking model tried (200 targets, flips, D3 at 1e-6): the lead measured from the stage's start (fresh centre: x* = sqrt(K_c) x_a/(x_a + sqrt(K_c))): TP 20, FP 11, FN 5; the pure sqrt(K_c): TP 25, FP 51, FN 0; the model used (displacement accumulated since the start of the run, lead since the stage's start): TP 24, FP 16, FN 1. The choice was made after seeing the counts (the model with the history is the one that matches the creep law of the repository's D'''.48/D'''.53, which integrates the displacement along all stages); no parameter was fitted.
6. The T1 constants (H1)-(H4) hypotheses for D_p were not checked. The planar/noisy/other-kernel cases untouched.

## 10. Reproduce

* `venv/bin/python experiments/theory3d_scan.py --run` (resumable; D3 on 133 targets x 15 eps, P0 on 55 x 5; about 35 min on a loaded core) -> `results/synthetic/theory3d_scan_raw.npz`.
* `venv/bin/python experiments/theory3d_readouts.py --run` (A-E, L, K, P in about 10 min, then F, G, H in seconds; `--analysis` recomputes F, G, H only; `--from-raw` rewrites `results/synthetic/theory3d.md` byte-identically).
* `venv/bin/python -m pytest tests/theory3d_test.py` (25 tests, about 20 s). Mutants: `tools/mutants/theory3d_unit_own.py` (50 + control), `theory3d_readouts_own.py` (43 + control), commands in their headers.
* Statement -> code: Thm 1: `theory3d_unit.W_p, W_2, uv_p, Psi_p`, `part_A`, `part_C`; Thm 2: `dormant_bound_p`, `w_of_distance_p`, `d_theta_p`, `mu_theta_p`, `ball_checks` B1-B5; Thm 3: `W_2, dormant_bound_2, d_theta_2, mu_theta_2, ball_checks` B6;
  Thm 4: `d1_fixed_centres, d1_formation_certificate`, `part_D`; Cor 5: `part_L`; Prop 6: `integrate_sliding_unit, drift_lead, lead_first_order_D3, overtake_margin`, `part_E`; the model on targets: `path_record, stage_pairs, overtake_score, critical_eps`, `part_F, part_G, part_H, part_K`.
* Commits (local, not pushed): 29c7483 (text, code, tests, mutation specs, report), 45dfc30 (raw data).
