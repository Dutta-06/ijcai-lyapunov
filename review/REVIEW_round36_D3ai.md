# Referee D3ai (round 36): Proposition D‴.47 (new), and the revisions of D‴.45 and D‴.46

## What I checked and how

**Rules followed**
- I worked only in a scratch copy of the repository at 7421bcf with a clean `git status`. The repository itself was never touched, and no
  git command that changes anything was run.
- Every process ran with `OMP_NUM_THREADS=1` and `nice -n 10`. I ran at most one heavy process at a time, and no single process ran longer
  than about 12 minutes. I never used `pkill -f` or `pgrep -f` and never set `EXP_CACHE`.
- Locations (in the round's scratch directory): scripts in `work/`, logs in `logs/`, mutant trees in `mut/`, mutant outputs in `mutres/`
  (with `summary.txt`).

**Reproduction**
- `tests/test_mirror_free.py`: 18 passed in 11 s.
- A full rerun of `exp_mirror_branch.py` (about 3 min) reproduces `results/chaos/mirror_branch.md` byte for byte.

**Transcription**
- I checked every number in numerics_D47/D45/D46.md against `mirror_branch.md`, `mirror_free.md` and `mirror_leak.md`.
- All are present and correctly transcribed, to the stated rounding: for example "0.5–10%" stands for 0.45–10.1%, and "below 10⁻¹⁰³" for
  7.22e-104.
- The problems are interpretive, and they are listed below.

**By hand**
- D‴.47(i)–(vi):
  - (i): the slope, the temperature, ψ(0), B(0, w) = β/(γm), the branch criterion, and the critical curve's −K|w|^{−1/3};
  - (ii): the proof of (a)–(c), including the zero-sum argument with Z′;
  - (iii): the rate derivation;
  - (iv): the conjugate and ψ_pair;
  - (v): the polytope and the bound;
  - (vi).
- The changed parts of D‴.45(ii)–(iv) and of D‴.46(i)–(iii).

**Numerical probes (scripts in `work/`)**
- `amax_check.py`: ParticleSet identities at a_max = 0.3 and 2.5. Every test and experiment uses a_max = 1. All identities hold: slope
  3e-9, log_mD 7e-15, `_dstate` exact, round trips 5e-14.
- `mixed.py`, `mixed_trace.py`: the mixed start, including a sign pattern that cannot fit y, and a solver trace.
- `limit_tests.py`: mixed starts with 2 free particles and inexact fits; a rank-4 design; ± pairs.
- `rank4*.py`, `partB_lsq.py`: `sign_ls_set` and lsq_linear.
- `dup_inexact*.py`: inexact fits where S is not a single point.
- `degenerate.py`: vertices with |S| < N.
- `cert_scan.py`, `cert_pct.py`, `cert_exact.py`: the nonnegative-certificate margin over all 3-supports, and 80-digit Chebyshev
  certificates.
- `e_generic.py`: the loop to t = 10⁹ on typical targets.
- `float_vs_exact.py`: where the float Newton fails.

**Mutation testing**
- 54 mutants of `mirror_free.py`:
  - D‴.47: ParticleSet, `sign_ls_set`, `dual_newton`/`_exact_newton_step`, `rich_leading_order`, `certificate_slacks`, `core_tail`;
  - the shared helpers `weighted_bp` and `nonneg_certificate`;
  - D‴.46: `filtered_gram`, `irrepresentable`.
- Each mutant was run against the 18 tests and against a reduced `exp_mirror_branch.py`, diffed against the unmutated run. The reduced
  experiment covers A, C, F; D at u = 10⁻⁴, 10⁻¹⁶, 10⁻⁶⁴, 10⁻¹²⁸; B seeds 0, 1, 3; and E's 21-atom case. The unmutated run matches the
  committed rows exactly.
- 5 report-code mutants, regenerated with `--from-raw`.
- Rounds 33/35's surviving mutants, re-created on the new code and run against the tests only.

**Literature**
- Web search only; arxiv.org and most mirrors were unreachable. Items marked "memory" were not verified.

## Verdicts

### Proposition D‴.47

**(i) Correct.**
- dη/dlog|w| = w/(mD) = a/(γa² + Gβa′²). It is signed: −1/(γa₀) on negative branches.
- The temperature 1/(γa₀), the expansion η = s log|w|/(γa₀) + C + O(|w|), the branch criterion, and K = 3.2183 on the critical curve all
  check.
- B(0, w) = β/(γm) is convention-free. ψ(0) = (Q + c)/(γm) is only a convention, because on a branch it depends on how ψ* is pinned.

**Overstated: the parenthesis "D‴.45(iii)'s rich limit needs α_i(0) = 0 to that precision".**
- Having a free sign is not enough. With θ = |α(0)|/(β(0)/2c)^{1/2} < 1:
  - β⁰ = (1 − θ²)β(0);
  - η(0) = (2/π)·arcsin θ·λ₁, which I verified to 4 digits at θ = 0.1, 0.5, 0.9 and ε = 10⁻⁶ to 10⁻¹⁰.
- So the rich weights become λ₁(1 ∓ (2/π)arcsin θ). This is sign-asymmetric unless θ → 0 (D3ai-12).

**(ii)(a)–(b) The proof is airtight.**
- The zero-sum identity Σ(η_i(t) − η_i(0))(w̄_i − w_∞,i) = zᵀΦ(w̄ − w_∞) = 0, with the −∞ terms on Z′, forces w̄_{Z′} = 0.
- After that, KKT holds on the face, and on the face the problem is strictly convex.
- Neither full row rank nor a bounded S is used:
  - sublevel sets of B(w̄, ·) are bounded by superlinearity;
  - range(Φ_Rᵀ) is closed;
  - KKT is sufficient under linear constraints.
- I confirmed the limit numerically, each time with a correct ŷ:
  - a rank-4 design with S unbounded: 9.5e-11;
  - duplicated atoms with an inexact fit and S not a single point: 4.7e-11 and 5.4e-9 (flow to t = 3·10⁹);
  - mixed starts with 2 free and 10 branch particles: ≤ 4.7e-10.

**What the numerics never test**
- In every inexact instance of Part B and of the unit test, |Z| = 7–8 and the active atoms are linearly independent. So S is a single
  point, and any convergent descent flow ends there.
- The non-trivial content of (ii)(b) — Bregman selection together with extinction — is untested (D3ai-4).

**The code fails exactly where the brief asked.** `sign_ls_set` silently returns a wrong ŷ on rank-deficient or duplicated designs
(D3ai-3).

**(ii)(c) Correct, with two misstatements.**
- "No other coefficient vanishes in the limit" is false for free particles, for example a zero column started at α = 0.
- The exponential rate does not apply on the critical curve (a₀ = 0): extinction there is algebraic, about t⁻³.
- Part B's fitted rates (0.5–10%) are consistent with the claim.

**(iii) The limit (Γ-convergence) and the rate derivation hold at an N-atom vertex.**
- z = log(1/ε)ξ + O(1), and log|w_j| = −μ_j log(1/ε) + O(1).
- "O(1) on S" is valid only on compacts: S is unbounded for pairs.

**At degenerate vertices (|S| < N)**
- On two instances with |S| = 3 < 8, the approach is still a power: ε^{0.13} and ε^{0.26}. The exponent is the slack (0.148) of the dual
  direction z/log(1/ε) that the entropy selects.
- `certificate_slacks` and `rich_leading_order` cannot compute this; they call `solve` on a non-square Φ_S.
- The hypothesis is stated in the draft. D‴.45's Part F and D‴.47's Part C exclude these instances (D3ai-6).

**The mixed-start claim "the free particles keeping O(ε² log(1/ε))" is false whenever the branch atoms' frozen signs cannot fit y.**
- On Part C's design with such an odd-atom sign pattern (the first one tried), Σ|w_free| = 1.36 at ε = 10⁻⁴ and 1.10 at ε = 10⁻⁸.
- Their η reaches λ₁ (7.03e4), i.e. a weighted ℓ¹ at weight ∝ ε^{−1/2}.
- Part C's own check holds by construction: the truth's atoms [1, 5, 9] are all odd, and they carry the truth's signs (`odd_feasible` =
  True).
- Where the claim does apply, Σ|w_free|/(ε² log(1/ε)) = 27, 43, 28, 27 over ε = 10⁻⁴ to 10⁻¹⁶. This is consistent with the bound (D3ai-2).

**The leading-order checks are honest but only asymptotic.**
- In Part D, |W − W₁| = 0.00995 against a distance of 0.0111 at u = 10⁻⁶⁴. The leading order beats the vertex only from u ≈ 10⁻¹²⁸.
- In Part C at ε = 10⁻¹⁶, the leading order is about 70% accurate.
- These "masses" (10⁻³² to 10⁻²⁵⁶) have no physical meaning (D3ai-14).

**(iv) Correct as convex analysis.** ψ_pair′ = asinh(W/2u)/(γa₀) and the curvature 1/(2γa₀u) both check. But it drops D‴.45(v)'s qualifier
"without cohesion between them" — a regression of D3af-6 (D3ai-8).

**(v) The proof is correct.**
- The polytope, S = {w*}, and the bound Σ_off w_i ≤ (|v|/μ)(2N·B/t)^{1/2} all hold.
- I checked the existence of a certificate for Gaussian translates with N = 2k + 1 in 80-digit arithmetic: it holds even where double
  precision returns margin 0.

**The claim that the loop "converges to w* from every positive-branch start, at every mass scale", and Part E's support for it, are
misleading.**

Margins over all 3-supports of Part E's designs:

| design | supports | median margin | margin ≤ 1e-9 | Part E's target | supports with a larger margin |
|---|---|---|---|---|---|
| 21 atoms, width 0.15, 7 samples | 1330 | 4.0e-6 | 338 (25%) | μ = 0.0146 | 10.8% |
| 41 atoms, width 0.15, 7 samples | 10660 | 3.8e-9 | 45% | μ = 0.003 | 5.6% |
| 41 atoms, width 0.1, 11 samples | — | — | — | μ = 0.027 | 8.3% |

- For the numerically degenerate supports, the exact certificate's margin is 1e-32 to 1e-34.
- In double precision, an LP finds nonnegative fits of the data carrying 0.41–1.39 of mass off the support. Uniqueness has no numerical
  meaning there.
- The loop at t = 10⁹ (Part E's law and integrator):

| target | margin μ | max\|w − w*\| | off-target mass |
|---|---|---|---|
| (3, 9, 14) | 1.9e-6 | 0.63 / 0.60 | 0.93 / 1.2, on 11 atoms |
| (0, 2, 4) | ≈ 0 | 0.49 / 0.45 | 0.60 / 0.66 |
| (2, 8, 17) | 3.4e-4 | 0.078 | — |
| Part E's own target | 0.0146 | 5.4e-5 | — |

The decay is about t^{−0.88} (not 1/t) for 41/0.15 (D3ai-1).

**(vi) Correct.** Part F's V̇ = −k_ee² check is a genuine finite-difference check.

**Novelty: little beyond the precedents.**
- (i) is immediate from β̇ = γβag.
- (ii) is Alvarez–Bolte–Brahic's Hessian–Riemannian convergence combined with the classical range argument (Gunasekar et al.). The boundary
  limits are handled as in ABB's linear-programming section (memory). The inconsistent case with extinction is a modest corollary; I found
  no paper stating it in exactly this form.
- Chou–Maly–Stöger, "How to induce regularization in linear models: a guide to reparametrizing gradient flow" (arXiv 2308.04921),
  characterises the implicit bias of reparametrised linear models by a Bregman divergence. It is not cited (search; abstract only).
- (iii) is Cominetti–San Martín and Weed at temperature 1/log(1/ε).
- (iv) is exactly Woodworth et al.'s / Ghai–Hazan–Singer's u² − v².
- (v)'s identifiability is Eftekhari et al.'s / Schiebinger et al.'s. Their stability constants, not uniqueness alone, govern the margins
  above.
- lit_S74 wrongly puts D‴.45's α = 0 rate under the exponential penalty "at temperature ≈ ε^{1/2}". That rate is the power ε^{3/2} from the
  core's tangent; the penalty at that temperature would give exp(−μ/ε^{1/2}) (D3ai-10, D3ai-11).

**Code and tests**
- The decimal fallback works:
  - the float Newton fails from ε = 10⁻³² (residual 0.645 at 10⁻³², 1.05–1.5 below);
  - the exact step converges in 36–55 iterations;
  - the two are bitwise identical for ε ≤ 10⁻²⁴.
- Its precision rule is untested. Mutant X2 (precision 40 digits) passes all 18 tests, and Parts C/D then silently tabulate failed
  projections (residuals 1.05–2.36, distance 0.75) at ε ≤ 10⁻⁶⁴. No part of the D‴.47 experiment asserts or counts residuals (D3ai-5).
- In an infeasible mixed start at ε = 10⁻¹² the masses span e^{−3·10⁵}. Each 3000-digit step takes about 10 s; after 50 iterations
  (286 s) the residual is still 0.96. The stall rule never fires (D3ai-13).

**Mutants that survive the tests (28 of 68)**

| target | mutant | what happens |
|---|---|---|
| `_dstate` | S1, sign dropped | tests pass; Part B hangs |
| `_dstate` | S3, δ dropped from β/δ | tests pass; Part B: flow 0.33 from the projection, B rises by 0.22; F: V̇ off by 9.7e-4 |
| `_dstate` | S4, speed ×(1 + δ) | tests pass; Part B: flow 0.071 from the projection; F: V̇ off by 7.4e-4 |
| `_eta_branch_log` | E4, O(δ) error | only A1's slope check moves, 2.3e-9 → 1.6e-5 |
| `sign_ls_set` | Z2, Z tolerance 1e-16 | \|Z\| changes; every check still passes |
| `log_mD` | L4, Hessian without Gβa′² | iterations and floor only |
| `dual_newton` | R2, radius never shrinks | iterations only |
| `dual_newton` | T3, no floor break | iterations only (34 → 3999) |
| `dual_newton` | T4, no fallback on failure | dead path in practice |
| `_exact_newton_step` | X2, precision 40 digits | tests pass; C/D projections fail silently at ε ≤ 10⁻⁶⁴ |
| `_exact_newton_step` | X4, X5 | redundant |
| `filtered_gram` | FG1, anti-causal filter | tests pass |
| `irrepresentable` | IR1 | tests pass |
| `nonneg_certificate` | NC1, margin capped at 1e-6 | tests pass |

- Equivalent or benign: P4 (Q is dead data on branches), E1, E2, E3, S2, B3, L3, Z1.
- Report code: all 5 report-code mutants survive (nothing checks the md).
- Killed by the tests: the δ(0) solve (P1–P3), bdiv B1/B2/B4, L1/L2/L5, Z3, N1–N4, R1, T1/T2/T5, X1/X3, O1–O3, C1–C3, K1–K4, WB1.

### Proposition D‴.45: the D3af items
- **D3af-1:** fixed. (iv)(b) now reads "settles on a sparse exact fit", and C(3) includes basis-pursuit≠truth instances. Leftover: "then
  atoms die", while the table's "dead" column is 0 throughout (Q → −c only asymptotically).
- **D3af-2:** fixed in D‴.45 (residuals returned, failures counted, 0 failures). Not carried over to D‴.47's experiment (D3ai-5).
- **D3af-3:** largely fixed (distribution, wrong-sign fraction, exact-fit count). numerics_D45 still quotes the bimodal s = 10⁻⁴ median
  6.6·10⁻⁴ without its spread (10/20 at ≤ 10⁻³, largest 0.37).
- **D3af-4:** fixed (Stratonovich everywhere; Pesme's Itô reading).
- **D3af-5, D3af-17:** fixed. Chou–Maly–Stöger is still missing.
- **D3af-6, D3af-7:** fixed in D‴.45 (re-broken in D‴.47(iv)).
- **D3af-8 to D3af-12, D3af-14, D3af-15, D3af-18, D3af-19:** fixed.
  - D3af-9: the −(3/2)K|w|^{2/3} remainder is correct.
  - D3af-12: the η code is checked at q = −c and on branches.
  - Re-created mutants: M10n, M11n, M24n are killed. M12n (Newton Hessian without m) passes in 292 s instead of 11 s; M22n and M36n are
    (near-)equivalent.
- **D3af-11:** fixed in the draft. ROW0/BULLET drop its qualifier: 1/β⁰² holds only for |S| < N.
- **D3af-13:** partly fixed. C(1)/C(4) still report Q kept by the mirror update as a result.
- **D3af-16:** "Parameters converge slowly" is still not shown.
- **D3af-20:** mostly fixed (V̇ test now sees the mass law, projections tested at 10⁻¹⁰), apart from M12n.
- **D3af-21:** not fixed: "falling 0.14, 0.16" and "as often as basis pursuit (0.75 against 0.80)".
- **D3af-22:** fixed.
- **New: the companion leading order (iii) is right asymptotically, but its constants are large.**
  - The relative error is about 27·β⁰^{1/2}/μ in Part C.
  - In Part F the largest closed-form errors are 113% and 138% at ε = 10⁻⁷ (among μ_min ≥ 0.01). The distance ratios range 185–2200.
  - numerics_D45 quotes medians only, and "O(β⁰^{3/2})" hides a μ⁻³ factor (D3ai-18).

### Proposition D‴.46: the D3ah items
- **D3ah-1:** reworded ("creeps toward"), and a case without persistent excitation was added. The σ = 0 run is still 1.3τ, not ≥ 5τ. The
  last quarter's 0.341 exceeds the stated limit 0.337 without comment.
- **D3ah-2, D3ah-5:** fixed. The draft now has H_L and d̄_L, the averaging is cited, and σ → 0 gives w* + H_L⁻¹d̄_L. The Bregman selection
  is correct because null(H_L) = null(H_Lᵀ) for H_L = k_e⟨φL[φ]ᵀ⟩.
- **D3ah-3:** fixed (runs start at w_σ^L, period means, lag split out).
- **D3ah-4:** fixed in the draft, the pieces and D(2). The AUDIT verdict still says "sparsity-selecting". The irrepresentable condition is
  stated unweighted although the weights λ₁(β⁰_i) differ for unequal masses (D3ai-20).
- **D3ah-6, D3ah-7, D3ah-8, D3ah-9, D3ah-10, D3ah-13, D3ah-14, D3ah-16:** fixed.
- **D3ah-11:** partly fixed. The test checks the bound where it binds. The experiment's assertion holds trivially at t = 0, so its printed
  "1.0000" is uninformative. By the reported numbers (C/λ = 1.74, 0.30 against late V = 0.085, 0.13), three of the five single bound
  mutants would still pass.
- **D3ah-12:** mostly fixed. R6, R7, K1, K2, P1, T1 are now killed. D2 (an Armijo step accepting a 10× larger residual) and L1 (`lasso` at
  k_e ≠ 1, still untested) survive.
- **D3ah-15:** fixed for continuous time (Q is read off the state). No Euler run.
- **D3ah-17:** partly fixed. Part D's ridge "support 6–12 atoms" is still a threshold artefact.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ai-1 | major | (v)'s "every start, every mass scale" rests on top-decile targets. Margins over all 3-supports: median 4.0e-6 (21 atoms) and 3.8e-9 (41 atoms); 25% and 45% numerically 0; exact margins 1e-32 to 1e-34; double-precision nonnegative fits with 0.4–1.4 off-support mass. At t = 10⁹ typical targets stay 0.45–0.63 away with 0.6–1.2 off-target mass. Part E's targets beat 89–94% of supports. | State the margin dependence in (v), BULLET, ROW0 and README; cite stability bounds, not uniqueness alone; report margin distributions and run Part E on random or median-margin targets (D47 (v); `e_case`). |
| D3ai-2 | major | (iii)'s "the free particles keeping O(ε² log(1/ε))" is false when the branches' frozen signs cannot fit y: free particles carry 1.36 and 1.10 at ε = 10⁻⁴ and 10⁻⁸, at the ℓ¹ wall λ₁. Part C's check holds by construction (the truth's atoms 1, 5, 9 are odd, with the truth's signs). | Add the hypothesis, or state the limit as weighted basis pursuit over the free atoms; add an infeasible-sign mixed row (D47 (iii); `part_C`). |
| D3ai-3 | major | `sign_ls_set` never checks lsq_linear's status. Rank-4 design, S unbounded: status 0 at max_iter, x ≈ 1e14, \|y − ŷ\| 0.103 instead of 0.520; the projection then runs 4000 iterations to residual 0.55. Duplicated atoms: status −1; \|Z\| 0 instead of 14, so the "projection" is 0.082 from the flow. With nnls's ŷ, flow = projection to 1e-10. | Use nnls/bvls or an active-set QP, check the status, verify KKT for ŷ, raise on failure (`mirror_free.sign_ls_set`). |
| D3ai-4 | major | The Bregman selection with extinction is never tested: every inexact instance (Part B and the unit test) has a single-point S, and exact-fit instances have Z = ∅. Mutants S1 (sign), S3 and S4 (branch speed) pass all 18 tests, and no test integrates a negative branch. | Add a test and a Part B row with non-singleton S and Z ≠ ∅ (two particles per atom, inexact) and negative branches; assert flow = projection (tests; `b_case`). |
| D3ai-5 | minor | The D47 experiment neither asserts nor counts projection residuals. X2 (40-digit precision) passes all tests and silently tabulates residuals 1.05–2.36 at ε ≤ 10⁻⁶⁴. | Count and exclude failures as in D45's Part B; add a 10⁻⁶⁴ fallback test (`exp_mirror_branch.py` C/D/B; tests). |
| D3ai-6 | minor | Degenerate vertices (\|S\| < N): the approach is ε^{μ*} with the entropy-selected dual (0.13, 0.26 observed), which the code cannot compute; these instances are excluded silently. | State what happens there; compute μ* from z/log(1/ε) (D47 (iii); D45 (iii); `certificate_slacks`). |
| D3ai-7 | minor | (ii)(c): "no other coefficient vanishes" is false for free particles, and the exponential rate fails on the critical curve (a₀ = 0, extinction about t⁻³). | Say "no other branch coefficient"; exclude Q = −c (D47 (ii)(c)). |
| D3ai-8 | minor | (iv) drops "without cohesion between them" (D3af-6 regression). | Restore the qualifier (D47 (iv)). |
| D3ai-9 | minor | BULLET's "masses decaying exponentially" is unconditional, against (v)'s algebraic decay. ROW0 calls (iii)'s rate "from the exponential-penalty asymptotics" while the Status says "proved". | Hedge and align the wording (pieces BULLET, ROW0). |
| D3ai-10 | minor | lit_S74 (and lit_S72's Cominetti/Weed entry) put D‴.45's α = 0 rate under the exponential penalty "at ε^{1/2}". It is the power ε^{3/2}. | Correct both entries (lit_S74, lit_S72). |
| D3ai-11 | minor | Novelty is overstated. Chou–Maly–Stöger 2023 is missing; (ii) is ABB plus the range argument, (iv) is Woodworth et al. / Ghai–Hazan–Singer, (v)'s identifiability is Eftekhari et al.'s. | Narrow the "New:" list and add the citation (Status, AUDIT, lit_S74). |
| D3ai-12 | minor | "α_i(0) = 0 to that precision" is not enough: η(0) = (2/π)arcsin θ·λ₁ and β⁰ = (1 − θ²)β(0), so the weights become λ₁(1 ∓ (2/π)arcsin θ). | Require θ → 0 and state the asymmetric limit (D47 intro; D45 (iii)). |
| D3ai-13 | minor | The decimal fallback can run for hours: about 10 s per 3000-digit step; 50 steps leave residual 0.96; the stall rule never fires. "37–55 steps" is instance-specific. | Add a time/iteration budget and a returned flag; report the cost (`dual_newton`). |
| D3ai-14 | minor | The rich-limit checks pass only at u, ε ≤ 10⁻¹²⁸ (Part D) and ≤ 10⁻¹⁶ (Part C, 70% there). The numerics omit the rows where the leading order fails, and the scales are unphysical. | Report all rows and state the reachable scales (numerics_D47 (C)–(D)). |
| D3ai-15 | cosmetic | A2 ("ten digits") and A4 are true by construction. E's columns labelled 10³, 10⁵, 10⁷ are at t = 1.22e3, 7.7e4, 8.8e6. "≈ 1/t" is t^{−0.88} for 41/0.15. "O(1) on S" holds only on compacts. Φ_S injectivity is implicit in (v). | Edit numerics_D47, the report labels and the draft. |
| D3ai-16 | minor | Untested solver and helper rules: R2, T3, T4, X4, X5, L4, E4, Z2, FG1, IR1, NC1 survive the tests; 5/5 report mutants survive. | Test causality, the irrepresentable condition and the margin; drop redundant rules; check the md (tests; `mirror_free`; experiments). |
| D3ai-17 | minor | D45 numerics leftovers: "falling 0.14, 0.16"; "as often as" (0.75 against 0.80); the bimodal median; "converge slowly" unshown; "atoms die" against a dead column of 0; Q kept by construction reported as a result. | Edit numerics_D45. |
| D3ai-18 | minor | D45 (iii)'s closed-form constants are large (about 27·β⁰^{1/2}/μ; 113–138% at 10⁻⁷) and O(β⁰^{3/2}) hides μ⁻³; numerics_D45 quotes medians only. ROW0 and BULLET give the 1/β⁰² time unqualified. | State the μ⁻³ factor, report largest values and ranges, qualify "\|S\| < N" (D45 (iii); numerics_D45 (F); pieces). |
| D3ai-19 | minor | M12n (Newton Hessian without m) passes, only 25× slower. | Bound the iterations or time in the tests (tests). |
| D3ai-20 | minor | D46 (iii): the irrepresentable condition must be weighted when the β⁰_i differ; "iff < 1" should be "< 1 sufficient, ≤ 1 necessary". | Fix the statement (D46 (iii)). |
| D3ai-21 | cosmetic | D46 leftovers: AUDIT "sparsity-selecting"; the ridge "support"; σ = 0 run at 1.3τ, with 0.341 > 0.337 unexplained; "V/bound 1.0000" is the t = 0 value; three of five bound mutants pass; D2 and L1 survive. | Edit numerics_D46 and AUDIT; print max_{t>0} V/bound; test `lasso` at k_e ≠ 1. |

## Lesson

The proofs are now sound: the zero-sum argument is airtight, and (iii)'s rate holds at N-atom vertices. Every failure is in what the
numerics were chosen to show:
- every inexact test instance has a single-point S, so "the flow lands on the Bregman projection" cannot fail there, and the dynamics'
  metric can be broken without any test noticing;
- the mixed start is checked only where the truth lies on the branch atoms;
- the positive-kernel targets come from the top decile of certificate margins, while the typical margin is 10⁻⁶ or numerically zero;
- the independent reference solver was exercised only where lsq_linear happens to converge, and its status is never read.

Before a number is called confirmation, ask what a wrong answer would look like on that very instance: whether S is a point, how the
margin is distributed over targets, whether the hypothesis holds by construction, and whether the reference solver reported success.

Sources cited by the referee: Chou–Maly–Stöger, arXiv 2308.04921; Chou–Maly–Verdun et al., AISTATS 2025; Eftekhari et al., ACHA 2021;
Alvarez–Bolte–Brahic (HAL hal-01928141).
