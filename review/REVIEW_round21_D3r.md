# Referee report D3r (round 21) — PARADIGM §5d: Proposition D‴.31 (statement (i)–(iv) and the closing paragraph, proof, Status, *Numerics*), the §0 and §7 rows "| D‴.31 |", the §8 bullet "Done (Proposition D‴.31)", `results/RESULTS_PARADIGM.md` (Proposition D‴.31), the LITERATURE audit bullet and [S54], `src/core/measurement_noise.py` (the D‴.31 part), `experiments/exp_eiv_remedy.py` → `results/chaos/eiv_remedy.md` (+ raw), `tests/test_eiv_remedy.py`, PAPERS.md paper 12

Branch `claude/festive-bardeen-eopej2` at 637a051, HEAD throughout the review. The concurrent session edited tracked files during the
review. At the start these were only files outside the material (`exp_gauge_floor.py`, `exp_servo_cogging.py`, the servo md and raw).
By the end they included PARADIGM.md, LITERATURE_TRANSPORT.md and RESULTS (other propositions' round-20 edits). No changed line
mentions D‴.31, `eiv_remedy`, the washout or the desired regressor, and `measurement_noise.py`, the experiment, its md and raw, the test
file and PAPERS are unchanged. Everything here is the committed 637a051 version: every script imports the 637a051 tree, extracted with
`git archive` into the scratch directory. Line numbers are 637a051's:
- PARADIGM: D‴.31 l. 3485–3586 — (i) l. 3491–3497, (ii) l. 3499–3508, (iii) l. 3510–3514, (iv) l. 3516–3532, the closing paragraph
  l. 3533–3535, the proof l. 3537–3551, the Status l. 3553–3567, the *Numerics* l. 3570–3586; the §0 row l. 135; the §7 row l. 6022; the
  §8 bullet l. 6277–6284.
- RESULTS l. 722–739. LITERATURE_TRANSPORT: [S53] l. 73, [S54] l. 74, the audit bullet l. 436–446. PAPERS paper 12 l. 165–182.
- `src/core/measurement_noise.py`: the docstring l. 18–23; `collective_block` … `washout_transfer` l. 62–105; `collapsed_noise_input`,
  `reference_compensation` l. 151–172; `NoisyRealLoop` l. 186–216; `mean_error`, `eiv_force`, `equivalent_perturbation`, `response_size`
  l. 254–303.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread; pytest ran with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the scratch tree rewrites `eiv_remedy.md`
byte-identically (md5 `fee61b56…`, equal to `git show 637a051:results/chaos/eiv_remedy.md`).

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3r/scripts/`:

| script | what it does |
|---|---|
| `inst.py` (module) | independent re-implementation of the instance: the saturated network in unconstrained coordinates with analytic derivatives, the exact representation at q = 1/3, Ĵ, ∂_xĴ, H, Q on 4096 orbit points, and the K₀-atom real loop with MR, DR and DR* (no repository code) |
| `recompute_partA.py` | every entry of Part A (24 rows) and Part B (5 priors) of `eiv_remedy.md` from `inst.py`, against the raw records |
| `check_identity.py` | (ii) beyond the scalar prior family (priors fixed on coordinate subspaces, the projected prior's active face, random and non-symmetric S₀, a block with a matrix anchor gain); (iii) for non-symmetric S; the minimum-norm characterisation |
| `check_washout.py` | (i) along a path with jumps; the OT leak with split rows (the repository's `plan_to` at incommensurate m) |
| `second_order_exact.py` | the exact second-order small-noise mean of the real K₀-atom loop at the lift: all quadratic terms of the true vector field contracted with the stationary periodic covariance of the linearised loop (the OU noise as a state), split exactly into the noise blocks and the parameter blocks of the covariance; the exact periodic second-order mean; MR, DR, DR* at k_e = 25 and 100, B = 100, η = ∞ and 0, and DR at η = σ₀ (the prior in the state) |
| `transient_diff.py` | Part C's transient at second order, exactly: the transient covariance from the lift, the second-order forcings, the means; the MR − DR difference and the curvature of the mean function at the records (`transient_exact.py`, a joint integration, timed out and is superseded) |
| `recompute_partC.py` | every number of Part C, the paragraph and the rows from the raw records, independently of `report()`; the sampling error of the MR − DR test; the decomposition of MR's Monte Carlo deficit; Proposition D‴.30's runs of the same loop |
| `reproduce_C.py` | Part C's runs of seed 0 (MR, DR, DR*) re-run with the committed `job_C` against the raw records |
| `check_tests.py` | the asserted quantities of `tests/test_eiv_remedy.py` with their margins; mutants (a regressor at the true state; a wrong prior rate; a wrong prior target) |
| `derived.py` | ratios quoted below, from the raw records and the saved second-order results |

Also in `ref_D3r/`: `tests.out` (5 passed, 13.77 s), `from_raw.out`, `md5_before.out`, `md5_after.out`. The experiment was not rerun in
full; Part A and B were recomputed independently, and Part C's seed 0 was reproduced bitwise for all three laws (confirmation 10).

## Verdicts

| object | verdict |
|---|---|
| **(i) the washout identity** | **holds**: it is the variation-of-constants formula of dz/ds = −η(z − T̂), exact along a path with four jumps (6.7·10⁻¹⁶). "Along every solution" needs T̂ of bounded variation, which the proof asserts ("it jumps where the plan does") but a sliding (chattering) solution of the vertex-selection flow need not have; the integrated form holds for every measurable T̂ (D3r-5). |
| (i) the leak decomposition | **holds for integral plans** (D‴.23's commensurate m, the nearest-atom cells). For the OT leak toward the plan-weighted prior with split rows (incommensurate m, which `plan_to` and D‴.27's proof handle) the barycentres also feel −σ₀(I − M)z, M_jl = (m/q_j)Σ_iγ_ijγ_il, which is not washed out: 2.49·10⁻³ at z = T̂ in an m = 7 example (D3r-5). |
| (i) the stiffness h_j + σ₀iω/(iω + η) | **holds per mode, for every η**: the prior and the leak are scalars in the whitened weighted coordinates, so they commute with H. It is a statement at the exact representation; "at a stationary point" is too general (D3r-5). The washout analogy ([S54]) holds for the barycentres only: the separations keep the full leak, so only the replicated equilibria are preserved. |
| **(ii) HR + B = I** | **holds for the prior family and for every S₀ that commutes with H; false as a matrix identity for a general affine time-invariant prior.** With R and B as the proof defines them the identity is RH + B = I. For non-commuting S₀ (a prior fixed on some coordinates, the linearised projected prior with an active face, a random S₀) max\|HR + B − I\| = 0.21–0.79, while RH + B = I holds to 10⁻¹⁴ (D3r-2). The per-mode statement h_jr_jj + b_jj = 1 holds for every S₀. The identity is the stationary equation itself; "no prior gives both" binds only on the modes with h_j < σ₀ (3 of 9 here). |
| **(iii)** | **holds** (re-proved by the congruence). Symmetry of S is not needed: S + Sᵀ ⪰ 0 suffices, and the exact condition is S + Sᵀ + SᵀH⁻¹S ⪰ 0 (D3r-6). The minimum-norm characterisation and the projection hold. Classical (D3r-3). |
| **(iv)(a)** | **holds.** The lift is a solution under DR and DR*. DR's linearisation is Corollary D‴.26's exactly; DR*'s has −(k_e − a) in the error channel. The averaged loops coincide at leading order: Theorem J's averaged field H(Θ, t) is already evaluated at x_m. Its proof carries over (Step 1 with r_Θ(x) replaced by f(x) − f̂_Θ(x_m) for DR*). Sadegh–Horowitz's gain condition concerns global stability, not this local statement. |
| **(iv)(b)** | **holds.** MR − DR is second order path by path: the test's difference ratio is 100.2 for a factor 10 in the noise's std. DR*'s covariance part is 0.9955 of MR's at t = 300. |
| **(iv)(c) the formulas** | **hold for the truncated force.** g for DR and DR* re-derived. The exact second-order computation reproduces the md's biases from the noise part of the covariance to ≤ 1.8·10⁻³ relative; the quasi-static E[νe], Var(e) hold to 0.31 % and 1.7 % at k_e = 25. The pointwise bounds hold, with the constant 1 nearly attained (0.954 at k_e = 25, B = 1000); DR* needs k_e > max a. |
| **(iv)(c) the structure** | **incomplete — the second-order force omits the parameters' own fluctuation** (D3r-1). With it the force is not a mean-tracking-error force: DR's adaptive-prior bias exceeds ‖g‖²_ρ by 3.3 % at k_e = 25, B = 100, and by 20.5 % at k_e = B = 100, where the md's 3.3·10⁸ is an artefact (2.1·10⁸). "The projected mean tracking error … grows like k_e²" is k_e times it; the mean tracking error itself falls with k_e (D3r-4). |
| closing paragraph | **holds for the truncated force**, and survives D3r-1: the omitted part falls with k_e (its own response 3.9·10⁻⁴ at k_e = 25, 2.9·10⁻⁵ at 100). The jitter formula is MR's and DR's; DR*'s is k_e²Eν²/((k_e − a)(k_e − a + B)) (D3r-7). |
| Proof | (i)–(iii) correct under the BV and commutation hypotheses. The first step of (iv)(c) drops the parameters' fluctuation (D3r-1). |
| *Status* | accurate in calling (iv)(c) formal and its consequences "exact for the averaged system". Overstated in: "(ii) proved" (false in its generality, D3r-2); "F₁ … the only part … that the adaptive prior amplifies" (D3r-1); the novelty (D3r-3). |
| *Numerics* (D‴.31) | **every number reproduces** (confirmations 2, 8, 10). Readings: "3.3·10⁸ … where DR's g nearly cancels" and "‖g‖²_ρ to 0.01 %" are properties of the truncation (D3r-1). "MR's grows like k_e²" is ×46.6 for ×64 (D3r-7). "0.5–5.3 % … the sampling error" holds within 0.43–1.83 sampling errors, but up to 2.8 % of it is the averaged prediction's own O(1/k_e) error, as the exact second-order transient shows (D3r-7, confirmation 9). The curvature of the mean function is negligible (confirmation 9). |
| `src/core/measurement_noise.py` | **holds**: H, Q, the forces, the biases and the prior family agree with an independent implementation to ≤ 4.6·10⁻¹⁰ relative. `stationary_operators` covers scalar κ only. `mean_error`'s docstring misnames a term (D3r-7). |
| `experiments/exp_eiv_remedy.py`, `eiv_remedy.md` | **hold**: regeneration byte-identical, seed 0 of all three laws reproduced bitwise. DR*'s prediction uses its own generator and noise input, as (iv)(a)–(b) say, but H for its mean (negligible). |
| `tests/test_eiv_remedy.py` | **passes (5/5); weaker than its names** (D3r-8). Test 1 is the variation-of-constants formula of the implemented equation, with frozen classes. Test 2 cannot reach D3r-2. Test 4 passes a regressor at the true state. Test 5 checks the force against its own definition. |
| §0 and §7 rows, §8 bullet, RESULTS, audit bullet, PAPERS 12 | **numbers hold.** Wording: D3r-1 ("least-squares", "uniformly in … the spectrum", "3.3·10⁸", "to 0.01 %"); D3r-2 ("HR + B = I", "no affine time-invariant prior"); D3r-3 (novelty); D3r-4 ("projected mean error 2394"). PAPERS omits (iv)'s status and says "(ii) proved" (D3r-7). |

## Confirmed issues

**D3r-1 (major; (iv)(c) l. 3524–3532; proof l. 3546–3549 "The mean of ε_m∇_ϑΦ̃(x_m; ϑ) is ∇_ϑΦ̃(x_m)E[e], since Eν = 0. For DR, E[e] is
Proposition D‴.30's …"; Status l. 3556–3562; the §0 row "its bias at most ¼‖∂²_xf‖²_ρ(Eν²)² uniformly in k_e, B and the spectrum"; the §8
bullet "turning the noise's persistent force into a least-squares one, its bias bounded uniformly in k_e and in the spectrum"; Numerics
l. 3572–3577; RESULTS l. 726–729; PAPERS l. 172–173; md Part A; `eiv_force`; test 5).** The second-order force omits the parameters' own
fluctuation. With it, the force is not a mean-tracking-error force, and the bound ‖g‖²_ρ fails on the instance.

*The omission.* The proof's step evaluates the regressor at fixed parameters. At the lift the parameters fluctuate at first order:
D‴.30(ii) gives mode j the stationary variance (N₀/2k_e)Q_jj(h_j + η)/(h_j(h_j + η + σ₀)). At second order this adds two terms:
- in the learning law, ∇²_ϑΦ̃(x_m)·E[ε_mδϑ]. This is the correlation of the regressor's dependence on the (noisy) parameters with the
  measured error. It is an errors-in-variables term in the parameters, and evaluating the regressor at the reference does not remove it.
- in the mean error, −(1/k_e)(½∂²_Θf̂[Σ_ΘΘ] + ∂_x∂_Θf̂·E[(δe + ν)δΘ]): the compensation's curvature in the parameters times their variance,
  and its slope–parameter cross term.

Both are O(Eν²), the order of g. Their size relative to g is set by the parameters' stationary variance, not by k_e and B alone.

*The computation* (`second_order_exact.py`). The exact second-order small-noise mean of the real K₀-atom loop at the lift:
- all quadratic terms of the true vector field (complex-step and central differences), contracted with the stationary periodic covariance
  of the linearised loop, the OU noise being a state;
- no averaging and no quasi-static step;
- the forcing is linear in the covariance, so it splits exactly into the part from its noise blocks (δe, ν) and the part from its
  parameter blocks.

The noise part reproduces the md to ≤ 1.8·10⁻³ (confirmation 4). With the parameter part, per (Eν²)², η = ∞:

| case | ‖g‖²_ρ (g from the noise part) | bias from the noise part (the md: 0.47053, 0.052621, 4.509·10⁻⁵, 0.32675, 2396.9 and 14949) | full second order | change | ‖g_full‖²_ρ (mean error incl. the parameter terms) | the learning-law term's own response |
|---|---|---|---|---|---|---|
| DR, k_e = 25, B = 100 | 0.47051 | 0.47047 | 0.48584 | +3.3 % | 0.48036 | 2.10·10⁻⁴ |
| DR*, k_e = 25, B = 100 | 0.05273 | 0.05269 | 0.04805 | −8.8 % | 0.04956 | 2.11·10⁻⁴ |
| DR, k_e = B = 100 | 5.843·10⁻⁵ | 4.501·10⁻⁵ | 7.043·10⁻⁵ | +56 % | 6.441·10⁻⁵ | 1.67·10⁻⁵ |
| DR*, k_e = B = 100 | — | 0.32678 | 0.32351 | −1.0 % | — | — |
| MR, k_e = 25 and 100 | — | 2397.66, 14950.1 | 2397.68, 14950.1 | < 10⁻⁵ | — | — |

Under the fixed prior (η = 0, k_e = 25, B = 100) the bias moves from 0.2106 to 0.2147 (+2.0 %, DR) and from 0.02363 to 0.02230 (−5.6 %,
DR*). With the prior in the state at η = σ₀ (DR, k_e = 25, B = 100) the full bias is 0.48597 (+3.3 %; ‖g_full‖²_ρ = 0.47932; the
learning-law term's own response 2.63·10⁻⁴; the weakest mode's variance 0.0206Eν²), so "the same for every η > 0" holds to 0.03 % here
(0.48597 against 0.48584 at η = ∞). The exact periodic second-order mean agrees with the averaged response to the full force to 0.2–0.3 %
for DR, so the excess is not an averaging artefact.

*Consequences.*
- *"By (iii) the bias is at most ‖g‖²_ρ" fails.* DR's full bias exceeds ‖g‖²_ρ by 3.3 % at k_e = 25, and by 20.5 % at k_e = B = 100. It
  exceeds even the least-squares bound of the full mean tracking error, ‖g_full‖²_ρ, by 1.1 % and 9.4 %: the learning-law term is not of
  the form ⟨Ĵᵀ(k_eE[e])⟩. (DR*'s full bias happens to lie below its ‖g‖²_ρ, the parameter part lowering it.)
- *The extremes of Part A are artefacts of the truncation.* At k_e = B = 100 DR's g "nearly cancels", and the omitted part raises DR's
  bias there by 56 % (it is 36 % of the full bias), so MR over DR is 1.4950·10⁴/7.043·10⁻⁵ = 2.1·10⁸, not 3.3·10⁸. Likewise "DR's bias attains ‖g‖²_ρ to 0.01 %" is a property
  of the truncated force; the parameters' fluctuation moves it by 3.3 %.
- *F₁ is not "the only part of the noise's persistent force that the adaptive prior amplifies"* (Status). The learning-law term is answered
  by H⁻¹ as F₁ is, and its equivalent perturbation ĴH⁻¹F_y has no bound in terms of the noise's moments.
- *"Uniformly in k_e, in B and in the spectrum of H"* holds for the truncated force only. The omitted part is linear in the parameters'
  stationary covariance, whose mode j is (N₀/2k_e)Q_jj(h_j + η)/(h_j(h_j + η + σ₀)). That has no bound uniform in the spectrum: it
  diverges as h_j → 0 at η > 0. Here the weakest mode's variance is 0.0404Eν² at η = ∞ against 0.0014Eν² under the fixed prior. On this
  instance the strong modes set the parameter part: the parameter force has about the same size under the fixed prior
  (FᵀH⁻¹F = 4.4·10⁻⁴ against 3.9·10⁻⁴).

*What survives.* On the instance the conclusions stand:
- In the computed cases DR's full bias stays ≤ 0.49, far below ¼‖∂²_xf‖²_ρ = 1.31, and three to four orders of magnitude below MR's.
- The omitted part falls with k_e (its own response 3.9·10⁻⁴ at k_e = 25, 2.9·10⁻⁵ at 100), so "the trade-off in k_e disappears" holds.
- MR's numbers are unchanged.

*Fix:*
- State (iv)(c) for "the part of the second-order force carried by the noise's own second moments", and add: "the parameters' first-order
  fluctuation adds a second-order force — the curvature of f̂ in the parameters and the correlation E[ε_mδϑ] of the regressor's parameter
  dependence with the measured error — of relative size set by the stationary parameter variance; it is not a mean-tracking-error force and
  is answered by H⁻¹ like F₁".
- Print it in Part A. It is cheap: D‴.30(ii)'s periodic covariance and one Hessian contraction. Here it is +3.3 % (DR) and −8.8 % (DR*) at
  k_e = 25, B = 100, and +56 % at k_e = B = 100.
- Replace "3.3·10⁸" by 2.1·10⁸, or drop the extreme. Write "to 0.01 % for the truncated force".
- In the Status, write "(iv)(c) formal and truncated (the parameters' fluctuation omitted)", and drop "the only part" and "uniformly … in the
  spectrum".
- For the D‴.30 referee: D‴.30(iii)'s F omits the same term. It is negligible for MR (2397.66 → 2397.68).

**D3r-2 (minor; (ii) l. 3499–3508 "y = RF + Bδ with R = (H + S₀)⁻¹ and B = (H + S₀)⁻¹S₀, so HR + B = I"; proof l. 3539–3541; the §0 row
"along every affine time-invariant prior … HR + B = I"; the §8 bullet; PAPERS l. 171–172 "no affine time-invariant prior … (HR + B = I)";
`stationary_operators`; test 2; md Part B).** HR + B = I is false for a general affine time-invariant prior; the identity is RH + B = I.
- *The algebra.* With R and B as defined, RH + B = (H + S₀)⁻¹(H + S₀) = I always. HR + B = I holds iff H commutes with S₀; for symmetric
  S₀, HR + Bᵀ = I.
- *The counterexamples* (`check_identity.py`; RH + B = I to ≤ 1.7·10⁻¹⁴ and h_jr_jj + b_jj = 1 to ≤ 10⁻¹⁴ in every case):

  | S₀ | max\|HR + B − I\| |
  |---|---|
  | the scalar family (κ/η = 0, 0.1, 1, 10) | ≤ 2.9·10⁻¹⁵ |
  | σ₀ on the weak eigenmodes (commuting) | 1.1·10⁻¹⁵ |
  | a prior fixed on the amplitude coordinates | 0.785 |
  | one active face of the projected prior (atom 1's amplitude) | 0.575 |
  | fixed on the centre coordinates | 0.673 |
  | random symmetric S₀ ⪰ 0 | 0.415 |
  | non-symmetric S₀ (a prior gain coupling the atoms) | 0.210 (and \|HR + Bᵀ − I\| = 0.24) |

  The linearised projected prior of D‴.27 with an active face is the natural case: there S₀ is σ₀ times a coordinate projection. A block
  with the matrix anchor gain dz/ds = −η(z − T̂) − σ₀P_amp(z − z_a) reproduces the proof's R and B (to 7.3·10⁻¹² and 2.1·10⁻¹⁵) and gives
  max\|HR + B − I\| = 0.471.
- *What holds.* "On every mode a vanishing bias forces the response 1/h_j" holds for every S₀, as the diagonal v_jᵀ(RH + B)v_j =
  h_jr_jj + b_jj = 1 (r_jj := v_jᵀRv_j, b_jj := v_jᵀBv_j).
- *Where the trade-off binds.* "No such prior both removes the leak's bias and keeps its bound" binds only where 1/h_j exceeds the fixed
  prior's 1/σ₀, i.e. on the modes with h_j < σ₀: 3 of the 9 here (h₁, h₂, h₃ = 3.06·10⁻⁴, 1.96·10⁻³, 3.11·10⁻³ against σ₀ = 7.07·10⁻³).
- *A hidden normalisation.* The anchor enters through I − 𝒦: the filter must return the anchor when the barycentre sits there. A prior
  whose anchor gain ℒ(0) differs from I − 𝒦(0) gives B = σ₀(H + S₀)⁻¹ℒ(0), and the identity fails. So does a constant offset when
  𝒦(0) = I, which no anchor can absorb.
- *Content.* The identity is the stationary equation rewritten: per mode, the Tikhonov filter factor h_j/(h_j + s) plus the residual factor
  s/(h_j + s) is 1 (D3r-3).

*Fix:*
- Write "RH + B = I, i.e. R = (I − B)H⁻¹: the response is the unregularised one times the share not taken by the anchor", or restrict to
  𝒦(0) commuting with H.
- State the per-mode form with r_jj and b_jj, and "the trade-off binds on the modes with h_j < σ₀".
- State the anchor's normalisation.
- Add a non-commuting S₀ (the active-face projection) to test 2, asserting RH + B = I.

**D3r-3 (minor; Status l. 3557–3562; the §0 row's last column; the audit bullet l. 436–446; PAPERS l. 177–179).** The novelty is
overstated.
- *The washout reading.* Yucelen and Haddad's low-frequency learning (IEEE TAC 58(4):1080–1085, 2013) has the same structure. Its update
  law carries a modification term that keeps the weight estimate close to a low-pass-filtered copy of itself — a leak toward a filtered
  estimate. It is introduced as filtering the high-frequency content of the adaptation while preserving the asymptotic stability of the
  error dynamics. That is the washout reading: full leak above the filter's bandwidth, none at zero frequency.
  - In elastic averaging the centre is "a moving average, taken in time and also in space" of the workers (Zhang–Choromanska–LeCun,
    NeurIPS 2015, cited by D‴.23 as [S49]); the elastic force is then the washout structure.
  - The bibliographic details were checked by a web search; the content is from the abstracts and the search snippets. The full texts
    could not be fetched (DNS and proxy), and the exact form of Yucelen–Haddad's law is from memory.
- *HR + B = I.* The audit bullet itself calls it elementary, while the §0 row, the Status and PAPERS list it as new. Per mode it is the
  Tikhonov filter factor plus the residual factor (D‴.23's own Status). In matrix form it is the Bayesian/Kalman decomposition of a
  regularised estimate: posterior mean = (I − K)·prior + K·least-squares estimate, with K = RH and I − K = B.
- *(iii)* is the minimum-norm solution of ⟨Ĵᵀg⟩ = F together with the contraction of the regularised fit (filter factors ≤ 1).
- *The split F₁ + ⟨Ĵᵀg⟩* is the classical errors-in-variables distinction between noise in the regressor, which biases, and noise in the
  output, of which only the mean biases. The audit bullet names this as the desired compensation law's classical motivation.
- *DR\* is not quite Sadegh and Horowitz's law.* Their desired compensation adaptation law carries a nonlinear damping term, or a gain
  condition, for global stability; DR* without it is a local statement. This is fine here, but it should be said.

*Fix:* cite Yucelen–Haddad 2013 (and EASGD's moving-average centre) for the washout reading. Label the identity, the bound and the split
classical. Restrict "new" to their transfer to the mass-constrained transport leaks (uniform in m) and the instance's quantification.

**D3r-4 (minor; (iv)(c) l. 3530–3532 "the adaptive prior's stationary point keeps the projected mean tracking error −ĴH⁻¹F₁, of squared
norm F₁ᵀH⁻¹F₁, which grows like k_e² while k_e ≪ B"; Numerics l. 3576–3577 "its projected stationary mean tracking error vanishes, against
MR's 2394"; RESULTS l. 728; the §7 row "the projected mean error 2394 for MR, 0 for DR and DR*"; md Part A's column).** This is k_e times
the projected mean tracking error.
- At the adaptive prior's stationary point the mean of the error equation gives k_eE[e] = g − Ĵȳ at second order. So
  P(g − Ĵȳ) = −ĴH⁻¹F₁ is the projection of k_eE[e]: the projected mean tracking error itself is −ĴH⁻¹F₁/k_e.
- Its squared norm is F₁ᵀH⁻¹F₁/k_e² ∝ (B/(k_e + B))², and it falls with k_e. Per (Eν²)², for k_e = 25, 50, 100, 200 (`derived.py`):
  - 5.69, 5.42, 4.94, 4.15 at B = 1000;
  - 3.83, 2.66, 1.49, 0.66 at B = 100.

  What grows like k_e² is the rescaled quantity, which is also the function-level bias.
- The md's 0.000e+00 for DR and DR* is zero by construction (F₁ := Fs − Fl, with Fl computed for the same variant), not a computed
  stationary point.

*Fix:* "the rescaled projected mean tracking error k_eP E[e] = −ĴH⁻¹F₁ … grows like k_e²; the mean tracking error itself stays O(1) in
k_e". Label the md's column "k_e²‖P E[e]‖²_ρ (zero by construction for DR and DR*)".

**D3r-5 (cosmetic; (i) l. 3491–3497 and its proof l. 3537–3539).**
- *Bounded variation.* "T̂_j has bounded variation in s: it jumps where the plan does" assumes finitely many plan switches on bounded
  intervals. A sliding solution of the vertex-selection flow (D‴.17), whose plan chatters, need not have them. The integrated form
  T̂ − z = T̂ − e^{−ηs}z(0) − η∫₀ˢe^{−η(s−τ)}T̂dτ holds for every measurable T̂.
- *Split rows.* The decomposition of the leak needs integral plans. Theorem D‴.23 assumes commensurate m, but `plan_to` and D‴.27's proof
  handle split rows, whose leak target is mΣ_jγ_ijz_j. Then the barycentres feel −σ₀(T̂ − Mz) = −σ₀(T̂ − z) − σ₀(I − M)z, and the second
  part is not washed out (`check_washout.py`: 2.41·10⁻³, equal to −σ₀(I − M)z to 1.3·10⁻¹⁸; 2.49·10⁻³ with z = T̂; 4.9·10⁻¹⁹ for an integral
  plan).
- *"Linearised at a stationary point".* H and h_j are defined at the exact representation. At a replicated stationary point the block is
  the whitened Hessian, and at a non-replicated one the separations couple to the collective modes.
- *"No leak opposes a persistent force"* should read "no leak opposes a persistent force on the class barycentres": the separations keep σ₀.
- *The washout analogy* holds for the replicated equilibria only. The non-replicated stationary points of D‴.23(i), held by the leak, are
  not equilibria of the unregularised flow.

*Fix:* add "for integral plans (commensurate m)" and "T̂ of bounded variation (otherwise in integrated form)"; "at the exact
representation"; "on the class barycentres".

**D3r-6 (cosmetic; (iii) l. 3510 "For every symmetric S ≥ 0").** Symmetry is not needed.
- The response is Fᵀ(H + S)⁻ᵀH(H + S)⁻¹F, and the bound holds iff S + Sᵀ + SᵀH⁻¹S ⪰ 0, so S + Sᵀ ⪰ 0 suffices.
- `check_identity.py`, 2000 draws each: the ratio to FᵀH⁻¹F is at most 0.421 (symmetric S ⪰ 0), 0.392 (skew S) and 0.428 (symmetric
  part ⪰ 0). It held in 445 of 445 draws with an indefinite symmetric part satisfying the exact condition. It fails, at 4.0, for
  S = −(h₁/2)v₁v₁ᵀ.
- H must be invertible (D‴.24's hypothesis).

*Fix:* "for every S with S + Sᵀ ⪰ 0" (this covers the non-symmetric zero-frequency leaks of D3r-2), H ≻ 0.

**D3r-7 (cosmetic; Numerics l. 3570–3586, md, rows, RESULTS, PAPERS, docstrings).**
- *Averaged system against the exact linearised loop.* The md's 4–5-digit biases are averaged-system values. The exact periodic
  second-order mean is 10.6 % above FᵀH⁻¹F for MR at k_e = 25 and 0.66 % at k_e = 100: O(k_e⁻²), the error filter of D‴.26(c). For DR it
  is 0.2 %. The Status says "exact for the averaged system"; the paragraph should too.
- *"MR's grows like k_e² for B ≫ k_e (3564 to 1.66·10⁵ from k_e = 25 to 200 at B = 1000)".* That is ×46.6 for ×64 in k_e². The law is
  (k_eB/(k_e + B))² (×46.7); at B/k_e = 5 it is not yet k_e².
- *The md's "¼‖∂²_xf‖²_ρ(1 + O(1/k_e)) = 1.3065e+00".* The number is ¼‖∂²_xf‖²_ρ alone.
- *The md's crossover std 4.41 (DR, k_e = B = 100).* It is far outside the small-noise regime, and it rests on D3r-1's artefact.
- *"Three to four orders of magnitude below at every prior".* The MR/DR ratios along the family are 2.3·10³–5.1·10³.
- *"… to a relative L²_ρ error of 0.5–5.3 % at every record, the sampling error of the run mean".* The ratio to the run mean's sampling
  error is 0.43–1.83. Up to 2.8 % of the error is the averaged prediction's own O(1/k_e) error: the exact second-order transient lies
  0.25–2.8 % above it along its direction, the most at t ≤ 125. The run mean is within 0.59–1.67 sampling errors of the exact transient
  (confirmation 9).
- *`mean_error`'s docstring, "and the plant's of its own jitter".* The term −f̂″E[eν] is the compensation's slope at the jittering state,
  correlated with the noise; the plant's own jitter cancels against the compensation's (f ≡ f̂).
- *The closing paragraph's jitter Eν²(k_e + a)²/(k_e(k_e + B))* is MR's and DR's; DR*'s is k_e²Eν²/((k_e − a)(k_e − a + B)). "For all k_e
  and B" needs k_e > max a for DR* (max\|a\| = 1.56 here).
- *PAPERS l. 175 "Proposition D‴.31 ((i)–(iii) proved; the remedy checked in the real loop on common noise paths)".* It omits (iv)'s
  status, and (ii) is not proved as stated (D3r-2). Part C checks F₁ through MR − DR, and that DR and DR* sit at their covariance part
  over 300 time units. It checks neither DR's own force (its predicted mean part is 7.1·10⁻⁴ of the covariance part) nor any stationary
  bias.

*Fix:* in the order above: "(the averaged system; the exact linearised loop is 10.6 % higher for MR at k_e = 25, O(k_e⁻²))";
"like (k_eB/(k_e + B))²"; drop "(1 + O(1/k_e))" from the md line; drop or flag the 4.41; "about three and a half orders"; "0.5–5.3 %
(0.4–1.8 sampling errors; the averaged prediction itself up to 2.8 % low at t ≤ 125, O(1/k_e))"; reword the docstring; give DR*'s
jitter; PAPERS as in D3r-1 and D3r-2.

**D3r-8 (cosmetic; `tests/test_eiv_remedy.py`, the Status's "Checked", RESULTS l. 735–739, the §7 row's test list).** The tests pass, are
not vacuous as implementation tests, and are weaker than their names (`check_tests.py`).

| test | asserted quantity | measured | tolerance | what it cannot see |
|---|---|---|---|---|
| 1 (the washout) | the identity at t = 3T₀ | 8.7·10⁻¹⁵ | 10⁻⁹ | It is the variation-of-constants formula of the implemented prior equation. It catches a doubled rate (2.97·10⁻³) and a wrong target, the weighted class mean (6.8·10⁻⁴). The classes are frozen: no plan jump, no split row, so "every path of the plan" is untested. |
| 2 (block, HR + B = I) | the scalar family | passes (≤ 2.9·10⁻¹⁵ with the independent H) | 10⁻¹² | Non-commuting S₀, where the claim fails (D3r-2). |
| 3 (the bound) | random F, S ⪰ 0 | passes | 10⁻¹⁰ relative | Nothing beyond the theorem. |
| 4 (DR keeps the first order) | dev ratio / dif ratio / dif/dev / DR* | 9.14 / 100.2 / 0.027 / 0.035 | (7, 14) / (70, 140) / 0.05 / 0.1 | A mutant with the regressor at the true state x passes all of them (9.14, 100.1, 0.026; it differs from DR by 2.1·10⁻⁶ against MR − DR's 5.2·10⁻⁵ at std 10⁻²). The test checks a noise-free regressor, not the reference. |
| 5 (least-squares forces) | Fr = Ĵᵀg/k_e; \|g\| ≤ bound; FᵀH⁻¹F ≤ ‖g‖²; F₁ above DR | passes; max \|g\|/bound 0.60 / 0.21 (DR / DR*, k_e = 25, B = 100), 0.82 / 0.092 (k_e = 100, B = 1000); F₁'s response over DR's bias 5.1·10³ and 5.6·10⁴ | 10⁻¹², 1, 10⁻⁹, 10³ | Fr and g come from the same Ee in `eiv_force`, so "least-squares" holds by construction and the third check is Bessel's inequality. The force is never compared with the dynamics, which give a force that is not least-squares (D3r-1). |

*Fix:*
- Test 1: a nearest-atom run with a crossing, the barycentre computed independently.
- Test 2: the active-face S₀, asserting RH + B = I.
- Test 4: a unit check of `rhs_noisy`'s regressor point (ν = 0, x ≠ x_m: DR's law is G∇_ϑ(eΦ̃(x_m; ϑ))).
- Test 5: the force against the exact second-order mean (the noise part to 10⁻³) or a Monte Carlo at the lift.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` rewrites `eiv_remedy.md` byte-identically (md5 `fee61b56…`). The test file passes: 5 tests in
   13.77 s.
2. **Part A and Part B from an independent implementation** (`recompute_partA.py`; analytic derivatives, 4096 orbit points).
   - All 24 rows: fixed, adaptive, ‖g‖²_ρ, the projected error and the crossover, to ≤ 4.6·10⁻¹⁰ relative; h_j to 1.2·10⁻¹⁴.
   - MR/DR: 5094, 3015, 4.58·10⁴, 1.27·10⁴, 3.32·10⁸, 5.65·10⁴, 1.83·10⁵, 2.86·10⁵ (truncated force).
   - The largest DR/DR* bias is 1.1821, and ¼‖∂²_xf‖²_ρ = 1.30646.
   - DR's bias over ‖g‖²_ρ is 0.999928 at k_e = 25, B = 100, 1.000000 at B = 1000, and 0.772 at k_e = B = 100.
   - The crossover stds are 0.00121 (MR), 0.0864 (DR), 0.2584 (DR*).
   - Part B: S₀, r₁ = 3269.4 … 135.50, b₁ = 0 … 0.9586, h₁r₁ + b₁ = 1 to 10⁻¹², MR's 2397 → 491, DR's 0.4705 → 0.2106, the anchor's biases —
     all to the printed digits.
3. **The tangent-space reading.** For a bump, ∂²_xφ = (1/s)∂_sφ − (A/s²)∂_Aφ exactly, so the f″ part of g lies in the span of Ĵ. Only the
   a·f″ part leaves it: 0.007 % of ‖g‖²_ρ at k_e = 25, B = 100, and 23 % at k_e = B = 100, where the f″ part cancels.
4. **The quasi-static formulas against the exact second order** (`second_order_exact.py`; periodicity of P and m₂ to 10⁻¹¹–10⁻¹⁵).
   - The noise part of the slow force in H's eigenbasis is (0.0215, 1.880, 0.666, 0.158, 0.434, 0.993, −4.724, −0.0406, 5.619) for MR,
     against D‴.30's (0.02133, 1.880, 0.6659, 0.1583, 0.4332, 0.9929, −4.724, −0.0415, 5.619).
   - The adaptive biases from it: 2397.66 (MR), 0.470472 (DR), 0.052686 (DR*), against 2396.95, 0.470534, 0.052621. At k_e = B = 100:
     14950.1, 4.501·10⁻⁵, 0.32678, against 14949.3, 4.509·10⁻⁵, 0.32675.
   - Along the orbit, E[νe] and Var(e) agree with the quasi-static formulas to 0.31 % and 1.7 % (MR/DR, k_e = 25), 0.07 % and 1.0 %
     (DR*), and 0.05 % and 0.15 % (k_e = 100). ‖g‖²_ρ is 0.470507 against 0.470568.
   - The weakest mode's stationary variance is 0.0404Eν², against (N₀/2k_e)Q₁₁/h₁ = 0.0417Eν².
5. **(i)** (`check_washout.py`): the identity with jumps to 6.7·10⁻¹⁶; the split-row remainder as in D3r-5.
6. **(ii) and (iii)** (`check_identity.py`): as in D3r-2 and D3r-6. The minimum-norm g = ĴH⁻¹F exerts F to 10⁻¹⁴, with ‖g‖²/FᵀH⁻¹F = 1 to
   10⁻¹².
7. **(iv)(a)–(b).**
   - The lift stays exact under DR and DR* (test 4's clean run, < 10⁻¹²).
   - MR − DR scales like the variance (100.2 for a factor 10 in std) while the deviations scale like the std (9.14).
   - DR*'s covariance part is 0.99549 of MR's at t = 300; its prediction uses −(k_e − a) and −k_e as (iv)(a)–(b) say.
8. **Part C from the raw records** (`recompute_partC.py`).
   - The z-scores run from −2.04 to −0.04 (MR), −1.71 to 0.82 (DR) and −1.80 to 0.79 (DR*).
   - At t = 300: 2.4310·10⁻⁶, 2.7989·10⁻⁷, 2.7468·10⁻⁷; ratio 8.685 (predicted 9.157); DR's and DR*'s mean parts are 7.10·10⁻⁴ and
     8.36·10⁻⁵ of their covariance parts.
   - The MR − DR test: relative errors 0.0052–0.0527, against the run mean's relative sampling error 0.0087–0.0386 — ratios 0.43–1.83.
   - *MR's deficit is a sampling effect of the common paths.* E‖f̂_MR − f‖² = mean‖Δ‖² + 2mean⟨Δ, f̂_DR − f⟩ + mean‖f̂_DR − f‖², with
     Δ = f̂_MR − f̂_DR. At t = 300 MC − prediction = −3.18·10⁻⁷, of which the cross term is −3.06·10⁻⁷: the 12 paths' mean first-order
     deviation along Δ. D‴.30's 12 runs on other seeds give ΣMC/Σpred = 1.0699 ("7 % above"); pooled, the 24 runs give 0.974 with
     \|z\| ≤ 1.1.
9. **Part C's transient at second order, exactly** (`transient_diff.py`: the transient covariance, nfev 4.9·10⁵; the forcings on a grid of
   step 0.05; the means).
   - The covariance part tr(HP_yy) equals the experiment's at every record to the printed digits (7.5790·10⁻⁸ … 3.0001·10⁻⁷).
   - The exact MR − DR difference departs from the averaged prediction Ĵ(ȳ_MR − ȳ_DR) by 0.60–4.97 % (relative L²_ρ). Its amplitude along
     the prediction is 1.0283, 1.0260, 1.0253, 1.0250, 1.0212, 1.0141 at t = 25–150, and 1.0025–1.0097 later: the averaged prediction is
     O(1/k_e) low early on.
   - The run mean's error against the exact transient is 0.73–2.55 %, i.e. 0.59–1.67 sampling errors (against the averaged prediction
     0.52–5.27 %, 0.43–1.83).
   - The curvature ½∂²_Θf̂[Σ_ΘΘ] of the mean function at t = 300 changes DR's mean part from 2.129·10⁻¹⁰ to 2.000·10⁻¹⁰ (4·10⁻⁵ of its
     covariance part) and MR's by < 10⁻⁴. In MR − DR it cancels at second order, the first-order fluctuations being common.
   - MR's exact mean part at t = 300 is 2.4968·10⁻⁶, 2.0 % above the averaged 2.4490·10⁻⁶.
10. **Reproducibility** (`reproduce_C.py`). Seed 0 of MR, DR and DR* re-run with the committed `job_C` reproduces the raw records
    bitwise: max difference 0.0 in err and dev for all three, 95–97 s per run; E‖f̂ − f‖² at t = 300 of 2.859157·10⁻⁶, 4.206586·10⁻⁸ and
    4.119774·10⁻⁸.
11. **The literature's bibliographic details** (web search): [S54] (Hassouneh, Lee, Abed, Proc. ACC 2004, paper FrA02.2; the abstract:
    advantages and limitations, generalized washout filters with intercoupled auxiliary states); Yucelen–Haddad, IEEE TAC 58(4):1080–1085,
    2013; Zhang–Choromanska–LeCun, NeurIPS 2015.

## Not checked

- The experiment was not rerun in full (36 runs). Seed 0 was reproduced bitwise for all three laws, the regeneration is byte-identical,
  and Parts A–B were recomputed independently.
- The exact second order was computed at k_e = 25 (η = ∞, σ₀, 0) and at k_e = 100 (η = ∞), both at B = 100; B = 1000 and k_e = 50, 200
  were not run.
- No rigorous treatment of the stochastic averaging or of the second-order expansion was attempted. The computation of D3r-1 is the
  small-noise expansion of the real loop to second order, exact at that order; it is not a Monte Carlo of the stationary biases.
- Theorem J's proof was checked for its structure (Steps 1–2 with DR's and DR*'s error equations), not rewritten.
- The full texts of [S54], Yucelen–Haddad 2013 and the EASGD paper could not be fetched (DNS and proxy); their content is from the
  abstracts, the search snippets and memory.
- `review/REVIEW_LOG.md` at 637a051 has no round-20 entry; rounds 18–19 were read.
- Proposition D‴.30 is reviewed separately. Its formulas were taken as given except where they enter D‴.31; there D3r-1's omission
  applies to D‴.30(iii) as well.
