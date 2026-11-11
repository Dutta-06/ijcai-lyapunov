# Referee report D3eo, round 146: Proposition D‴.44 as revised after round 34

## Scope
I reviewed `review/drafts/D44.md` at b9bd244 with its supporting files: `experiments/exp_mass_dictionary.py`, `src/core/mass_dictionary.py`, `results/chaos/mass_dictionary.md` (+ `_raw.npz`) and `tests/test_mass_dictionary.py`. I read it against:
- round 34 (`review/REVIEW_round34_D3ag.md`; REVIEW_LOG l. 2235);
- the round-141 correction of D‴.45 (iv)(b) (PARADIGM l. 5600–5620);
- the round-142 correction of D‴.41's statuses (PARADIGM l. 153, 5273);
- Paper 16's Proposition (Euler) and Corollary 5.3 (`manuscript/paper16/paper16.tex` l. 410–460).

The repository was not modified. I worked in W = …/scratchpad/ref146 (`repo/` is the archive of b9bd244 with the venv symlinked; my scripts are in `own/`; downloads are in `dl_abs/` and `dl_pdf/`). All runs used OMP_NUM_THREADS=1, nice -n 10 and one process at a time.

## Summary
**Design.** I did six things:
- re-derived (i)–(iv) by hand;
- checked every closed form with my own scripts;
- measured K's spectrum, the masses and Q along Part B's free white reduced flow out to s = 10⁸, which tests Katzenberger's range;
- regenerated the md from its raw file and ran the tests;
- ran 10 mutants of my own;
- read Li–Wang–Arora and Pillaud-Vivien–Reygner–Flammarion in full, and seven abstracts.

**Verdicts.**
- Every identity in (i)–(iv) is correct:
  - the gradient structure and tr Ψ(K);
  - the replicator drift and its exact cancellation for free masses;
  - dν ≠ 0;
  - the Euler-lag reading and Q's drift with κ = 2 and κ = 1;
  - D*, ζ and their asymptotics;
  - the excess bound;
  - the null-kernel lemma.
- The md regenerates byte for byte, and every number in the draft is the md's: 160 of the 173 decimals verbatim, and 13 re-rounded by hand from md values, correctly.
- The round-34 issues are resolved.

Three gaps remain. The first is major: (ii)(b) is not restricted to the range where it holds.

**Findings.**
- **D3eo-1 (major).** The proposition drops Katzenberger's range, and (ii)(b) describes "ends" that lie outside it. The relevant facts:
  - The proved conditions hold at rest points with every β_i > 0. Paper 16's Corollary 5.3 was narrowed to exactly that case.
  - The ends observed are sparse limits in which the unused masses fall like 1/s. There K loses rank: 3 of its 5 eigenvalues fall from ~10⁻³ to 4·10⁻¹³ by s = 10⁸.
  - At Part B's noise level, the normal rate falls below the selection rate from s ≈ 3·10⁵ on.
  - Part B validates the reduced flow only at s = 1500.
- **D3eo-2 (moderate).** D‴.41's round-142 status change is not propagated. "(i) proved at order v given D‴.41(iv)" should be split:
  - proved (given Katzenberger/LWA) only for fresh draws as h → 0, on compact parts of M where rank K = N;
  - otherwise derived under an assumed fixed-step analogue of Katzenberger's theorem.
- **D3eo-3 (moderate).** Q's drift is stated without its sign, and the two updates have opposite sign structures.
  - The paradigm's update (κ = 2) has the corrected D‴.45 sign map: no particle crosses onto a branch.
  - The flat update (κ = 1) lowers Q on every free curve and on the critical curve. Its zero, which is D*'s split, lies only on branches, at Q = −c(1 + t⁴)/(1 − t⁴) < −c. The flat update therefore drives each alive atom onto a branch, where its sign is frozen. The draft never says this.
- **D3eo-4 (minor).** Three phrases claim more than the text proves:
  - "local minima of a concave penalty reached": the reduced flow is non-gradient by (i)(c), and P is convex below the threshold;
  - "minimising sequences inflate … toward the boundary": the bound only gives E/m ≤ 0.316·a_max;
  - (ii)(b)'s proof says "at an end the drift vanishes" without saying that this is for a rest point inside M.
- **D3eo-5 (minor).** "White" means h → 0 (f ≡ 1). At a finite step, fresh draws have f = 2/(2 − hκ), and then (iv)'s constant is not constant: the spread is 5.8·10⁻⁶ at h = 1, against the slow noise's 1.15·10⁻⁵. Part B's fresh rows use the finite-step filter, while C's white flows use f ≡ 1.
- **D3eo-6 (minor).** (i)(b) says the Itô term is "small only away from the saturation". That is incomplete: the ratio γa²/(2Gβa′² + γa²) also tends to 1 for light atoms. At the paradigm's equilibrium it equals (1 + t²)/(3 − t²) ∈ [1/3, 1).
- **D3eo-7 (minor).** Literature is missing:
  - Kunin et al., Even et al., Pesme et al. and Nacson et al. on finite steps breaking conservation laws. Paper 16 already cites them; the draft still lists "Q's drift for every colour" as new.
  - Li–Wang–Arora §6.1, where the reduced flow's stationary points likewise sit on the boundary of the manifold.
  - Pillaud-Vivien–Reygner–Flammarion eq. (7)–(10): their log-β Itô drift δh is exactly the "masses' own selection" that (i)(b) cancels.
- **D3eo-8 (cosmetic).** ζ denotes both the drive and the mirror rate; "per unit v" and "per unit s" are mixed; "before the projection" is missing in (i)(b); the null kernel is still compared with "the truth's 0.485"; 13 hand-rounded decimals.

## What I checked and how
- **[read]** D44.md in full. Round 34's report and its log entry. The PARADIGM rows and statements of D‴.41, D‴.43 and D‴.45, with the round-141 and round-142 corrections. Paper 16's Proposition (Euler) and Corollary 5.3. The experiment, the library and the tests.
- **[own] Re-derived by hand:**
  - J's gradient in θ = (α, ϱ) with mobility P = diag(G, γ/4);
  - ∂D_i/∂ϱ_i = (2ϱ_i/m²)(2Gβa′² + γa²);
  - the flat ϱ-drift −(γϱω/(2Nm))(2Gβa′² + γa²) and the Itô ϱ-drift γ²a²ωϱ/(2Nm), hence the ratio and the exact cancellation;
  - ν and dν = Σ(4γω_i/(Nm))ϱ_ia_ia′_i dα_i∧dϱ_i;
  - the paradigm's drift −(2ω/(Nm))(Gβa′(γa + Gβa″), γGβa′²) = −(2ω/(Nm))dv/dζ;
  - Q's drift as (2ω/(Nm))vᵀ∇²Qv, using a′²cosh(2α/a_max) = 1 − t⁴ and the second-order coefficient of ΔQ in both charts: (g²γβ/4)(κγa² − 2Gβ(1 − t⁴)), κ = 2 for (α, log β), κ = 1 for (α, √β);
  - D*, its split, both branches, the envelope slope (2γa_max/m)t/(1 + t²) and convexity;
  - ζ = (2γa_max/(Nm))t(1 − t²)/(1 + t²) and its asymptotics;
  - (iii)'s chain, both Cauchy–Schwarz steps and the boundary limit;
  - (iv)'s determinant lemma.
- **[own] `own/closed_forms.py`:**
  - sign maps of the κ = 1 and κ = 2 brackets on Q = ±0.5c, 0, −c and −2c, −10c (branches);
  - the zero sets Q = −c(1 − t²)/(1 + t²) for κ = 2 and −c(1 + t⁴)/(1 − t⁴) for κ = 1, matched to machine precision;
  - D*'s argmin against the root of its split (to 10⁻⁶), and D*′ by finite differences against the closed form (7 digits);
  - Q/c at the flat split (−1.014, −5.37, −144.5 at |w| = 10⁻⁴, 10⁻², 0.3) and at the paradigm's split (−0.90, −0.17, −0.0069);
  - the Itô/flat ratio at the paradigm's equilibrium (0.357, 0.746, 0.986).
- **[own] `own/katz_range.py` (output in `katz_range.out`):**
  - Part B's instance, free masses, white noise, the paradigm's reduced flow to s = 10⁸;
  - K's spectrum, the masses, Q/c, and λ_min(K) divided by the fastest relative mass rate at Part B's noise (σ = 0.1, h = 1).
- **[own] `own/null_finite_h.py`:** R_null − R_free with Damian–Ma–Lee's Ψ = −(2/h)log(1 − hκ/2), at six points of M_null.
- **[own] `own/mutants.py`:** 10 single-replacement mutants of the library (κ, the gauge k, ζ's factor 1 − t², the bound's Σ|φ|⁻², the cross term's factor, the Euler lag's a″, the Itô /4, the gauge ratio's t⁴, the relevance), each run against the test file in a scratch copy.
- **[own] `own/numcheck.py`:** every decimal in the numerics paragraph against the md.
- **[md]/[raw]:**
  - `MD_OUT=W/out python experiments/exp_mass_dictionary.py --from-raw` with the raw copied to W/out; `cmp` reports the md identical.
  - `pytest tests/test_mass_dictionary.py`: 21 passed in 2.5 s.
- **[paper]** Li–Wang–Arora (arXiv:2110.06914): Assumption 3.1, Theorem 4.6's stopping times, §6.1. Pillaud-Vivien–Reygner–Flammarion (arXiv:2206.09841): §2.4, Theorem 3. Both via pdftotext.
- **[snippet]** arXiv abstract pages: HaoChen–Wei–Lee–Ma 2006.08680, Shalova–Schlichting–Peletier 2404.12293, Kunin et al. 2012.04728, Even–Pesme–Gunasekar–Flammarion 2302.08982, Pesme–Pillaud-Vivien–Flammarion 2106.09524. export.arxiv.org was refused by the proxy; arxiv.org/abs and /pdf worked.

## Verdicts per part
**(i)(a)** Correct.
- The flow is a gradient flow in |dα|²/G + (4/γ)|dϱ|², of J or of its 0-homogeneous extension.
- The renormalisation argument is now the right one: the radial drifts are P-normal and the projection removes them.
- K's nonzero spectrum equals the normal rates (Part A: 3.7·10⁻¹¹ and 2.6·10⁻¹¹).
- Status: see D3eo-2. The selection is D‴.41(iv)(b)'s, which round 142 relabelled "derived … under an assumed fixed-step analogue". Only the fresh, h → 0 case is covered by a theorem, Katzenberger's (for this gradient flow, LWA Theorem 4.6), and only under Assumption 3.1 (rank ∇²L = codim on M) on compact parts of M. The draft states that hypothesis nowhere. D‴.43's status had "where rank K = N and λ_min(K) ≫ v_h"; D44 lost it.

**(i)(b)** Correct.
- The ratio γa²/(2Gβa′² + γa²) is right.
- The cancellation is exact for free masses and every colour (Part A: 1.7·10⁻¹⁵ and 1.0·10⁻¹⁵).
- The conserved leftover is the cross term −(γ²/(2N))a_iφ_iᵀf(K)y.
- The reading "small only away from the saturation" is incomplete (D3eo-6).

**(i)(c)** Correct. ν and dν are as stated; dν ≠ 0 on T M when d > N; the loop integrals scale with ratio 4.00, against ≤ 10⁻¹⁹ for the flat drift.

**(i)(d)** Correct.
- The Euler-lag identity holds.
- Q has no martingale part, because ∇Q·v ≡ 0 exactly.
- Q's drift is (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)); it agrees with D‴.45 and Paper 16 (b) for κ = 2 with ω_i = |φ_i|².
- What is missing is the sign map (D3eo-3). For κ = 2 the bracket is 2·[−G(1 − t²)(β⁰ − (2c − β⁰)t²)]:
  - negative on Q ≥ 0 and in the core of −c < Q < 0;
  - zero on |κ′t| = 1;
  - positive beyond, on the critical curve and on every branch.

  This is the corrected D‴.45 (iv)(b), and my own run confirms it: the dead atoms' Q/c rises toward −1 from above (−0.9846, −0.9957, −0.9988, −0.99986 at s = 10⁶ … 10⁸), never crossing it.
- For κ = 1 the bracket is −2G[Q(1 − t⁴) + c(1 + t⁴)]:
  - negative on every curve with Q ≥ −c (on the critical curve it is −4Gct⁴);
  - zero only on branches.

  The flat update's rest split, which is D*'s, therefore sits on a branch (Q/c = −1.014 to −144.5 above), with the atom's sign frozen (D‴.47). This is the opposite of the paradigm's "no crossing", and it is the content the round-141 correction asks to be checked.

**(ii)(a)** Correct:
- D*, its split, ≈ (γa_max/m)|w|, the 4/3 law and its constant;
- D*′ = (2γa_max/m)t/(1 + t²);
- convexity, D*′(0) = 0 and the dense ℓ^{4/3}-perturbation;
- the zero-coefficient condition φ_iᵀμ = 0.

Add the branch location (D3eo-3).

**(ii)(b)** The algebra is correct:
- the gauge equilibrium Gβ(1 − t⁴) = γa², i.e. "D*'s split with G halved";
- ζ and its asymptotics;
- the logarithmic P above the threshold.

But:
1. **The range (D3eo-1).** The proof's "at an end the drift vanishes" makes sense only at a rest point inside M, with every β_i > 0. The ends the draft talks about, and the only ones the numerics find, are sparse limits:
   - 2.4–2.6 atoms used in C (k = 2), 2 alive from s = 10⁶ in A;
   - the unused masses decay like 1/s (3.3·10⁻⁴, 3.1·10⁻⁵, 3.0·10⁻⁶ at s = 10⁶, 10⁷, 10⁸);
   - (Q, α) of the dead atoms tends to the corner (−c, 0), where the mirror chart degenerates (λ₁ ∝ β⁰^{−1/2});
   - K loses rank: its three smallest eigenvalues fall to 1.5·10⁻⁷ (s = 3·10⁵), 3.0·10⁻⁹ (10⁶) and 3.9·10⁻¹³ (10⁸), while the two alive atoms' modes stay at 6.3·10⁻³ and 3.1·10⁻³.

   On a k < N truth the projection drops these directions in the code itself (RCOND = 10⁻¹⁰). λ_min(K) over the selection's relative rate, at Part B's noise, is:

   | s | ratio |
   |---|---|
   | 10³ | 6.7·10³ |
   | 3·10⁴ | 1.2·10⁴ |
   | 10⁵ | 1.8·10³ |
   | 3·10⁵ | 1.76 |
   | 10⁶ | 0.38 |
   | 10⁷ | 0.036 |
   | 10⁸ | 0.0065 |

   So:
   - Extinction, and the choice of which sparse fit is reached, happen where the reduced flow is not the scheme's limit at that noise. For any fixed v the separation fails as s → ∞.
   - This is the situation of LWA §6.1: the stationary points lie on the boundary of Γ, "where the limiting flow itself is not even defined", and LWA Theorem 6.1 controls it with an ε–T statement.
   - Part B validates the reduced flow only at s = 1500, where the ratio is about 10⁴.
   - Paper 16's Corollary 5.3 states exactly the restricted version: rest points with every β_i > 0, −c < Q_i < 0, no zero coefficient, the sparse ends outside the theorem's range. D44, its source, does not.
2. **The penalty's role (D3eo-4).** "Local minima of a concave penalty reached according to the path" is not shown. The flow is non-gradient by (i)(c), and Σω_iP is not shown to be a Lyapunov function, so the rest points are critical points of Σω_iP, of unknown stability. P is convex (the 4/3 law) below the threshold and concave only above it.

**(ii)(c)** Correct (min(k, N) log τ_c/τ_c; rank K; the non-uniform O(1/τ_c) term).

**(iii)** The bound is correct (both Cauchy–Schwarz steps; it holds with the conserved K; it is vacuous for a_max ≥ 3.16), and so are the boundary limit and the threshold. But "the white potential's minimising sequences inflate the coefficients toward the representability boundary" says more than the bound gives: E/m ≤ 0.316 at a_max = 1, while the flat ends, presumably interior minimisers, have E/m = 0.025–0.034. The flows' inflation and the near-vertex ends are labelled numerical; that is correct.

**(iv)** Correct for f ≡ 1 and the flat update. The slow-noise lemma is right. The restriction to white noise holds only as h → 0 (D3eo-5). With Ψ = −(2/h)log(1 − hκ/2), R_null − R_free ranges over [−1.0936, −1.0877]·10⁻³ at h = 1, and [−1.2580, −1.1585]·10⁻³ at h = 10, against −1.082422·10⁻³.

**Numerics** All reproduced:
- Part A as printed;
- Part B's projected advances;
- C–E per instance, as the md prints them.

## Status of the round-34 issues
| id | status |
|---|---|
| D3ag-1 | Resolved in substance by (i)(c), (i)(d) and (ii)(b), and the loop integrals confirm it. The "ends" still need restricting: D3eo-1. |
| D3ag-2 | Resolved: the ratio and its O(1) at the ends; along the flows 0.037 → 1.000 (at most 1.304). Wording: D3eo-6. |
| D3ag-3 | Resolved: the cross fitness; the pushes balance from s = 10⁵. |
| D3ag-4 | Resolved: white noise and the flat update; the slow lemma; per-instance disagreements. The finite-step caveat is new: D3eo-5. |
| D3ag-5 | Resolved: rows split by bp_ok, with n and the largest distance. |
| D3ag-6 | Resolved: flat ends use 12.0 of 12; the paradigm's ≤ N + 1 (5.9–6.0, at most 6). |
| D3ag-7 | Resolved: projected comparison; σ²τ_c/h → 0 stated; the ratios tend to 1. The validation covers one early window only: D3eo-1. |
| D3ag-8 | Resolved: 21 tests; my own 10 of 10 mutants killed. |
| D3ag-9 | Resolved ("dense ℓ^{4/3}-perturbation"; the distances printed). |
| D3ag-10, -11, -12 | Resolved. |
| D3ag-13 | Resolved (n per row; the stopped flow listed with its last state). |
| D3ag-14 | Resolved. |
| D3ag-15 | Resolved: FlowFailed after Radau; null-kernel E/m 0.514 = 1 − 0.486. |
| D3ag-16 | Partly resolved: the cancellation and the branches of D* are credited, but the finite-step conservation-law literature is not: D3eo-7. |
| D3ag-17 | Resolved, except the null kernel's "(the truth's 0.485)" remnant: D3eo-8. |

## Overall assessment
D‴.44 now gets its mathematics right, its numerics regenerate exactly, and its tests are sharp. One thing is still needed before Paper 16 can lean on it: say where the reduced flow is the scheme's limit.
- Restrict (ii)(b) to rest points with every β_i > 0, as Corollary 5.3 does.
- State Katzenberger's (LWA's) hypothesis: rank K = N on compact parts of M.
- Label the sparse ends as reduced-flow statements outside that range.
- Propagate D‴.41's corrected status.
- Give Q's sign map for both updates.

I recommend a minor revision, with D3eo-1 as the one major item. It needs restricting and labelling, not new runs. A run of the scheme from the reduced flow's state at s ≈ 10⁵ at a smaller σ would be welcome but is optional.

## Issues
| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3eo-1 | major | Katzenberger's range is dropped. (ii)(b)'s "at an end every alive atom rests…" and Status "(ii)(b) the ends' conditions proved (necessary)" cover sparse limits. Along Part B's free white flow, K's three smallest eigenvalues fall from ~10⁻³ to 3.9·10⁻¹³ (s = 10⁸), the dead masses fall like 1/s, and λ_min(K) over the selection rate drops from 6.7·10³ (s = 10³) to 1.76 (3·10⁵) and 0.0065 (10⁸) at σ = 0.1, h = 1 [own `katz_range.out`]. Part B validates only at s = 1500. LWA §6.1 meets the same boundary stationary points. Paper 16's Cor. 5.3 is already restricted. | D44 (i) and Status: add "on compact parts of M where rank K = N (LWA's Assumption 3.1)". (ii)(b): "at a rest point with every β_i > 0 … (then −c < Q_i < 0, |κ′_it_i| = 1, no zero coefficient)", and say that sparse ends lie where the normal rates vanish, outside the theorem's range. Status: the ends' fits are those of the reduced flow, numerical, not the scheme's. Optional: B at s ≈ 10⁵ at smaller σ, or an LWA Thm 6.1–style ε–T statement. |
| D3eo-2 | moderate | Round 142 changed D‴.41(iv)'s status to "derived … under an assumed fixed-step analogue of Katzenberger". D44's Status still reads "(i) proved at order v given Proposition D‴.41(iv)". | Status (i): proved for fresh draws as h → 0, given Katzenberger/LWA within the range of D3eo-1; for a finite step and for coloured noise (Part A's AR(1), C–E's slow columns), derived under the assumed fixed-step analogue. The algebraic parts — the cancellation, ν, the Euler lag, Q's zero quadratic variation — stay proved. |
| D3eo-3 | moderate | Q's drift is given without its sign. κ = 2 has the D‴.45-corrected map: positive on the critical curve and on branches, so no crossing (own run: Q/c of the dead atoms → −1⁺). κ = 1 is negative on every curve with Q ≥ −c and vanishes only on branches at Q = −c(1 + t⁴)/(1 − t⁴). D*'s split, where the flat ends rest (gauge ratio 0.5000), has Q/c = −1.014 … −144.5, a branch with frozen sign [own `closed_forms.py`]. | (i)(d): the factored brackets −2G(1 − t²)(β⁰ − (2c − β⁰)t²) for κ = 2 and −2G[Q(1 − t⁴) + c(1 + t⁴)] for κ = 1, with their sign regions. (ii)(a): the flat update pushes alive atoms across the critical curve and ends on branches (D‴.47). (ii)(b): the equilibrium has −c < Q < 0. |
| D3eo-4 | minor | (ii)(b)'s "local minima of a concave penalty reached": the flow is non-gradient ((i)(c)), there is no Lyapunov function, and P is convex below the threshold. (iii)'s "minimising sequences inflate … toward the representability boundary": the bound gives only E/m ≤ 0.316·a_max, while the flat ends have E/m = 0.025–0.034. | "Critical points of Σω_iP (concave above the threshold), reached according to the path (numerical)". (iii): "the minimisers' excess is at most …". |
| D3eo-5 | minor | f ≡ 1 is the h → 0 limit. At a finite step fresh draws have f = 2/(2 − hκ) (D‴.43(ii), Damian–Ma–Lee), and (iv)'s difference varies: spread 5.8·10⁻⁶ at h = 1, 1.0·10⁻⁴ at h = 10 [own `null_finite_h.py`], against 1.15·10⁻⁵ for the slow noise. B's fresh rows use the finite-step filter; C's white flows use f ≡ 1. | Setting: "white noise as hλ_max → 0 (f ≡ 1); at a finite step f = 2/(2 − hκ)". (iv): add "as h → 0". Numerics: say which filter each part uses. |
| D3eo-6 | minor | (i)(b) says "small only away from the saturation". The ratio γa²/(2Gβa′² + γa²) also → 1 for light atoms at fixed amplitude; at the paradigm's equilibrium it is exactly (1 + t²)/(3 − t²) ∈ [1/3, 1). | (i)(b): "of order one wherever Gβa′² ≲ γa² (saturation or small mass); (1 + t²)/(3 − t²) at the gauge equilibrium"; "before the projection". |
| D3eo-7 | minor | Literature. The finite-step breaking of conservation laws (Kunin et al. 2012.04728; Even et al. 2302.08982; Nacson et al.) and the Itô shrinking of the effective initialisation (Pesme et al. 2106.09524) are cited in Paper 16 but not here, while D44 calls "Q's drift for every colour" new. LWA §6.1 meets boundary stationary points. PVRF eq. (7)–(10): the log-β Itô drift δh (the weighted-Lasso term) is the masses' own selection that (i)(b) removes. | Classical list: add them. "New:" narrowed to the sign map in both charts, the gauge equilibrium, the logarithmic penalty, the conserved-mass inflation and the null kernel. Cite PVRF as the quadratic-case precedent of the cancellation's mechanism. |
| D3eo-8 | cosmetic | ζ denotes both the drive and the mirror rate. The cross fitness is "per unit v" while the drifts are per unit s. The null kernel is compared with "the truth's 0.485" (the free end is 0.468). 13 decimals are hand-rounded rather than printed by the md (0.049, 0.056, 0.116, 0.124, 0.497, 0.528, 0.685, 0.701, 4.8·10⁻⁴, 6.6·10⁻⁷, −4.1·10⁻⁷, −2.11·10⁻⁵, −9.60·10⁻⁶; all correct roundings) [own `numcheck.py`]. | Rename the drive. One unit throughout. Compare the null kernel with the free end. Print the rounded values in the md, or quote the md's digits. |

## Sources
- [paper] Z. Li, T. Wang, S. Arora, *What happens after SGD reaches zero loss? A mathematical framework*, arXiv:2110.06914: Assumption 3.1 (rank ∇²L = M on Γ), Theorem 4.6 (stopping times, compact K), §6.1 (stationary points on the boundary of Γ, where the flow is undefined), Theorem 6.1 (ε–T).
- [paper] L. Pillaud-Vivien, J. Reygner, N. Flammarion, *Label noise (stochastic) gradient descent implicitly solves the Lasso for quadratic parametrisation*, arXiv:2206.09841: §2.4 eq. (7)–(10), Theorem 3.
- [snippet] J. Z. HaoChen, C. Wei, J. D. Lee, T. Ma, arXiv:2006.08680; A. Shalova, A. Schlichting, M. Peletier, arXiv:2404.12293; D. Kunin et al., arXiv:2012.04728; M. Even, S. Pesme, S. Gunasekar, N. Flammarion, arXiv:2302.08982; S. Pesme, L. Pillaud-Vivien, N. Flammarion, arXiv:2106.09524.
- [read] PARADIGM.md (D‴.41 row l. 153, round-142 block l. 5273, D‴.43, D‴.45 (iv)(b) l. 5600–5620); `manuscript/paper16/paper16.tex` l. 410–460; `review/REVIEW_round34_D3ag.md`; REVIEW_LOG rounds 34, 141, 142, 145.
- [md]/[raw] `results/chaos/mass_dictionary.md`, regenerated byte-identically into W/out.
- [own] Scripts in W/own (W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref146): `closed_forms.py`, `katz_range.py` (output in `katz_range.out`), `null_finite_h.py`, `mutants.py`, `numcheck.py`, and `dq_debug.py` (a finite-difference check of ΔQ that roundoff ruins; superseded by the analytic second derivative).

## Lesson
Revising a proposition into the right algebra does not finish it: the range of the limit theorem has to travel with the formulas. The quantity that decides this, λ_min(K) against the selection rate, takes one run along the flow to print. Here it shows that the sparse ends, the outcome the proposition is about, lie where its own reduced flow has stopped being the scheme's limit. A sister manuscript (Paper 16) had already restricted its corollary for exactly this reason, and a source that is still unreviewed should be checked against such downstream corrections, both the range and the sign of a drift on every curve, before it is cited.
