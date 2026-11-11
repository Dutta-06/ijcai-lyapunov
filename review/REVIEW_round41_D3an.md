# Referee report D3an (round 41): Theorem D‴.53 (the sliding limit) and Lemma D‴.53a

## What I checked and how

**Setup and rules.**
- All work was done in a scratch copy, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review41/repo`. Its D‴.53 objects are byte-identical to 7f5d1e4, and the repository's later HEAD (aab77ca) has not changed them.
- Every Python process ran with venv/bin/python, `OMP_NUM_THREADS=1` and `nice -n 10`, one heavy process at a time. Nothing in the repository was written to.
- The container restarted mid-review. I reran the interrupted jobs (the smaller-margin targets).
- Procedural slip: one wait loop I launched contained `pgrep -x -f`. It only listed processes and killed nothing. I stopped it within seconds and used `ps` from then on.

**Tests and reproduction.**
- `tests/test_transit_limit.py`: 2 passed (15 s).
- `test_transit_mass.py` and `test_rich_transit.py`: 20 passed.
- A full rerun of the experiment (`TL_WORKERS=1`, `TL_OUT` in scratch) reproduces `results/chaos/transit_limit.md` byte for byte.
- I checked every number in the numerics paragraph against `transit_limit.md` and `transit_limit_raw.npz`.

**Lemma D‴.53a.**
- Re-derived the gradient-flow identity by hand.
- On the real grid (61 samples of [−1, 1], σ = 0.12): found the critical points of the single-kernel fit (zeros of d log(⟨φ*φ(c)⟩²/⟨φ(c)²⟩)/dc on [−4, 4]) and E's infimum as c leaves the samples.
- Integrated the isolated-entrant inner problem in 5 cases near the ends of the interval.

**Step 1.** Computed the exact profile W^ε(ρ(1 − ε^{1/2}ξ)) from the closed-form mirror map (`FreeLaw.eta` on Q = ε − c) and compared it with 𝒲(ξ), for ε = 10⁻⁴ … 10⁻⁸.

**`loop_with_duals`.**
- Independent integration: (α, c) only, one segment, no dual integration, η recovered from α by the closed form. Run at ε = 10⁻⁶ and 10⁻⁸.
- Its own tolerance convergence: rtol 10⁻¹¹ → 10⁻¹² → 10⁻¹³ at ε = 10⁻⁹ … 10⁻¹¹, on all five targets.

**`entry_pull`.**
- Convergence of the pulls and of s_θ in xi0 (200 … 12800) and in tolerance (rtol 10⁻⁹ → 10⁻¹¹).
- Re-predicted all five targets with converged values (tight tolerance, Richardson from xi0 = 1600/3200) and re-compared.

**Conclusion (b) as stated.** On target 39 at ε = 10⁻⁷ … 10⁻¹⁰, the block's distance from S_k at four points: τ_k + 400ε^{1/2}, τ_k + 20ε^{1/2}log(1/ε), the midpoint, and τ_{k+1} − ε^{1/4}.

**Target selection.**
- For all 48 targets: eligibility, margins, approach speeds v_k, slides, s_θ, σ_k and pulls.
- The full comparison (converged predictions, loop at rtol 10⁻¹²) on seven more targets (31, 4, 16, 0, 8, 14, 5), and target 12 at the committed tolerances, at ε = 10⁻⁷ … 10⁻¹¹.

**Mutation testing.** 27 mutants of `loop_with_duals`, `creep_constants`, `entry_pull` and the experiment's recursion. Each was run against the D‴.53 tests and against the experiment's metrics (target 39 at 10⁻¹⁰, target 41 at 10⁻⁹). 22 were killed and 5 survived.

**Literature.** Web search only. arxiv, ar5iv and infosciences.epfl.ch were not reachable (the latter two would need to be allowed in the session's network settings). Comparisons marked [memory] are from memory and not verified.

## Verdicts

### Lemma D‴.53a
- **The identity** (the inner problem is a gradient flow of E with mobilities |𝒲′(ξ_i)|/ρ_i for w and Λ₀Gm for c): correct and exact.
- **Isolated entrant in the continuum** (unique critical point with w ≠ 0, compact sublevel sets, LaSalle): correct.
- **On the finite grid and near the interval's ends: not proved.**
  - The finite critical point is unique numerically (c* = 0, 0.5, 0.8, 0.95), and the inner problem converges in all 5 boundary cases I ran.
  - But the compactness step fails, because ⟨φ(c)²⟩ → 0 off the samples (D3an-3).
- **Coupled case: not proved** (D3an-3).

### Theorem D‴.53

**Hypotheses.** (H1)–(H3) are the substance. Also used but not stated:
- no exits (only implied through (H2));
- θ must be a regular value of the entrant's coefficient (needed for (c));
- the inner problem's limit must be identified with the path's polished S_k;
- the code assumes equal masses.

More importantly, every constant depends on quantities missing from the statement: the approach speeds v_k, the entrants' slides, the gaps between entries, the later entrants' margins at earlier entries, and powers of μ₀ beyond μ₀⁻³ (D3an-1, D3an-10).

**Proof, step by step.**
- **Step 1:** the profile is accurate to O(ε^{1/2}) in relative terms, not O(ε) (D3an-6).
- **Step 2:** correct in substance.
  - The fast equilibrium is exponentially stable uniformly (linearisation −DH_k).
  - The feedback from the slow (inactive) kernels is contained in p.
  - The bound O(ε^{3/2}μ₀⁻³ + ε^{3/4}) is right, with a constant proportional to v_{k+1}⁻³. The ε^{3/4} term is actually attained at τ_{k+1} − ε^{1/4}.
- **Step 3:** the "bounded fundamental matrix" claim is false. The conclusion can be repaired (D3an-2).
  - Λε^{1/2} = Λ₀ exactly.
  - The O(ε^{1/2}) perturbation comes from the profile.
- **Step 4:** the matching and the 1/X tail are right. Numerically the tail is 1/xi0 with a 1/xi0² second order, so Richardson extrapolation works.
  - The shift term σ_k(a_{k−1,j} − a_{k,j}) is **derived, not fitted**: it is a first-order Taylor expansion in the arrival time.
  - It is confirmed: flipping its sign, halving it or dropping it is each caught by the test.
  - Its tail coefficient grows like g′²/v³, which is what breaks the expansion for slow entrants (D3an-1).
- **Step 5:** εD_j is formally O(ε); its interaction with the pulls is O(ε^{3/2}μ₀⁻⁴). But D_j ranges from about 2.7·10² to 1.2·10⁴, and the induction drops Step 3's error (D3an-7, D3an-11).

**Conclusions.**
- **(a):** the bound is right, and the coefficients follow the core's profile to O(ε^{1/2}) (1.1–1.9·10⁻³ at 10⁻¹⁰). The printed formulas for C_j and D_j are wrong (D3an-5). The displacement ratio converges like ε^{1/2} as claimed (3.1–3.3 per decade).
- **(b):** right, but tested only at the midpoints (D3an-12).
- **(c) and (d):** asymptotically correct on every target I ran, and better supported than the committed numerics show (D3an-4). But the range of ε where the expansion is useful varies by orders of magnitude from target to target (D3an-1).
  - "+ εD_j + o(ε^{1/2})" is empty, because εD_j is itself o(ε^{1/2}).
  - "The next order in u is O(ε)" is shown only numerically (D3an-7).
- **Status line:** overstated. "Proved under (H1)–(H3)" rests on Step 3's false claim and Step 1's wrong rate. "(H3) proved for isolated entrants" holds only in the continuum with one kernel. "(H3) verified numerically at each entry" is not in the committed code (D3an-9).

### Numerics
- **Independence:** the comparison is genuinely independent; the loop never calls `entry_pull`, and I found no common-mode error.
  - `loop_with_duals` matches the independent integration to 2.7·10⁻⁷ at ε = 10⁻⁶ and 4.6·10⁻⁶ at 10⁻⁸ (units of ε^{1/2}).
  - Tightening its rtol to 10⁻¹³ moves the observed values by at most 3.3·10⁻⁵ at 10⁻¹⁰ and 2.0·10⁻⁴ at 10⁻¹¹.
  - The per-component atol scaling matters: the unscaled mutant fails at 6.2·10⁻³.
- **The residual's decay is real.** The weak part is the predictions, not the loop (D3an-4).
- **Signs:** correct, and recorded per entry. That choice is not protected by the test (D3an-8).
- **Selection is biased** (D3an-1).

### Novelty
The method is classical. The structural picture has close precedents that D‴.53 does not cite (D3an-14).

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3an-1 | major | **The selected targets hide constants that are not uniform.** The remainders depend on the entrants' approach speeds and slides, the gaps between entries, the later entrants' margins at earlier entries, and powers of μ₀. None of these appear in the statement. The inner problem's far-field creep tail has coefficient ∝ g′²/v³: about 2.5·10² for target 41's 2nd entry (v = 0.05) and about 7·10⁶ for target 12's 3rd (v = 0.0013). The eligibility filter excludes none of the 48 targets. The five chosen are the top margins, and all have v ≥ 0.044 and slides ≤ 0.06. On 8 of 48 targets (5, 12, 15, 27, 28, 35, 45, 46; v = 0.0013–0.030, slides 0.06–0.30), \|s_θ+σ\| is 6.9·10²–9.8·10⁴ and the pulls are 28–3.7·10³. The top five have ≤ 135 and ≤ 5.9. **Target 12** (μ₀ = 0.025): the 3rd entrant slides 0.30, and s_θ ≈ −9.42·10⁴ (the committed recipe at xi0 = 400/800 gives −9.22·10⁴). The loop crosses at 16.41, 15.83 and 14.44 for ε = 10⁻¹¹, 10⁻¹⁰, 10⁻⁹, where the path has τ₃ = 16.71 and the prediction is 16.40, 15.73, 13.60. The residual decays like ε^{1/2} (3.4·10³ε^{1/2}, then 1.0·10⁴ε^{1/2}). The tracers sit −3.7·10³ε^{1/2} off the path (predicted −3.9·10³). **Target 8** (entries 0.014 apart): the block is 0.76, 0.27 and 3·10⁻⁵ off S_k at 10⁻⁷, 10⁻⁸, 10⁻⁹. The 2nd crossing is off by 5.0, 3.6, 1.6, 0.56, 0.18 (units of ε^{1/2}) for 10⁻⁷ … 10⁻¹¹. | Statement and Steps 3–4: make the constants explicit, C(μ₀, min v_k, slides, gaps, later margins, λ), and state a regime, e.g. ε^{1/2}\|s_θ,k+σ_k\| ≪ gaps and ε^{1/4} ≪ v_k³/g′_k². Numerics: report all 48 targets (v_k, slide, s_θ) and include targets 12, 5 and 8 in the experiment. |
| D3an-2 | major | **Step 3's "bounded fundamental matrix on (−∞, S]" is false.** The entrant's ξ_e depends on its centre through the O(1) coefficient s⟨r_{k−1}∂_cφ_e⟩/ρ. So a perturbation δc_e(s′) produces δξ_e(s) ≈ (s−s′)g′_eδc_e, which grows linearly. On top of that come the neutral tracer and translation directions. "Its centre's coupling is O(ξ⁻³)" bounds the c_e equation, not ξ_e's dependence on c_e. The conclusion O(ε^{1/2}(X+S)) survives because the perturbations of the c_e equation carry the factor w_e = O(ξ⁻³), and the mismatch of c_e at the cut is O(X⁻²), giving O(1/X) in ξ_e. That constant is g′²/v³, the one in D3an-1. | Step 3: replace the claim with a weighted Gronwall estimate in (ξ_e, c_e) with explicit constants, or use a blow-up of the wall. |
| D3an-3 | major | **The lemma's coupled case and "(H3) proved for isolated entrants" are not proved.** (i) "{w_e=0} lies above the sublevel set by ½A²G(0)−O(δ)" needs S_{k−1} to be (up to O(δ)) the global minimiser of the block's sign-constrained fit; (H2) gives only a local one. The block kernels' own degenerate sets {w_i=0} are not addressed (an exit, or the entrant taking over a block kernel's feature). The implicit function theorem only locates critical points; convergence to the one near (S_{k−1}, (A, c*)) needs a shadowing or basin argument. The limit is never identified with the path's polished S_k. (ii) On the grid the compactness step fails, because ⟨φ(c)²⟩ → 0 off the samples. For c* = 0.8, E → 0.990E₀ as c → ∞, and E = 0.973E₀ at c = 1.5 with w = 2.3·10³. For c* = 0.95, E = 0.68E₀ at c = 1.5. The lemma's conclusion still holds numerically. | Lemma and status: say "proved in the continuum, one kernel". Add the hypothesis inf E_block > E(S_{k−1}) − ½A²G(0) + O(δ) and a shadowing proof, or call the coupled case a conjecture. On the grid, require E(s₁) < liminf_{\|c\|→∞} E. |
| D3an-4 | major | **The predictions in the numerics are not converged, and the paragraph misreads floors and trends.** Richardson from xi0 = 400/800 with `entry_pull` at rtol 10⁻⁹ is off by up to 0.014 in s_θ: target 41's 2nd entry gives 103.385 against 103.371 (from 6400/12800 at rtol 10⁻¹¹), target 24's 3rd gives 76.376 against 76.362. The pulls are off by up to 7·10⁻⁴ (4.8149 against 4.8142). With converged predictions and the loop at rtol 10⁻¹³, the residuals at 10⁻¹¹ are 3.1·10⁻⁵, 5.8·10⁻⁴, 3.1·10⁻⁵, 2.8·10⁻⁵ and 5.3·10⁻⁵, not (2–8)·10⁻⁴. Four of the five targets fall by 2.8–3.3 per decade (ε^{1/2}) down to 10⁻¹¹; target 6 is irregular. The crossings at 10⁻¹¹ are ≤ 1.8·10⁻³, except target 24's 3rd entry at 8.5·10⁻³ (0.094, 0.029, 0.0085 over 10⁻⁹ … 10⁻¹¹). So: "(2–8)·10⁻⁴ … (the integration's floor)" is misattributed, since the loop's own error there is ≤ 2·10⁻⁴. "Falls like ε^{1/2} (slopes 0.32–0.54)" is not what the committed slopes show (three of the five are 0.32–0.36). "(c) … falling like ε^{1/2} to 0.003–0.017" is contradicted by the committed data: target 41's 2nd crossing goes 0.0013, −0.0089, −0.0126, −0.0166 (10⁻⁸ … 10⁻¹¹), and target 6's goes 0.0047, 0.0071, 0.0089. Target 24's −0.0063 at 10⁻¹¹ is an accidental cancellation; the true value is +0.0085. The theorem comes out better than reported. | `exp_transit_limit.py`: XI0 = (1600, 3200) with `entry_pull` at rtol 10⁻¹¹, report the spread of the extrapolation, and run the loop at rtol 10⁻¹³ for ε ≤ 10⁻¹⁰. Rewrite the (c) and (d) sentences of the numerics paragraph. |
| D3an-5 | major | **Statement (a)'s formulas are wrong as printed.** C_j = ∫κ(μ_j)s_j a_{k,j}dτ′ and D_j = ∫C_j a_{k,j}dτ′ use a_{k,j} = ⟨r_kφ(c_j⁰)⟩, the rate. The creep law and `creep_constants` use ⟨r_k∂_cφ(c_j⁰)⟩. The printed D_j, run as a mutant, fails the test (2.8·10⁻³ > 2·10⁻³) and puts target 41 at 10⁻⁹ off by 1.2·10⁻² instead of 5·10⁻⁴. | (a): replace a_{k(τ′),j} with ⟨r_{k(τ′)}∂_cφ(c_j⁰)⟩ in both integrals. |
| D3an-6 | minor | **Step 1's "W^ε(sρ(1−ε^{1/2}ξ)) = s𝒲(ξ) + O(ε)" is false.** The relative error is O(ε^{1/2}): W^ε/𝒲 − 1 = 6.86·10⁻², 6.76·10⁻³, 6.76·10⁻⁴ at ξ = 0 for ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, and the same for ξ from −8 to 10. It comes from D‴.51's second integral, ∫(D_ε−D₀)dα = O(ε^{1/2}), which shifts ξ by O(ε^{1/2}) (O(ε) in u). Step 3 already uses O(ε^{1/2}), so the orders survive. This offset feeds the O(ε) remainders in (c) and (d). | Step 1: write "+O(ε^{1/2}) (O(ε) in u)". |
| D3an-7 | minor | **(d)'s "ε^{1/2}dev+εD_j+o(ε^{1/2})" is empty in its εD_j term, and implies εD_j is the O(ε) correction.** With converged predictions, the residual left after εD_j is still O(ε) (about 3.2 per decade in units of ε^{1/2}), and amounts to 1–31% of εD_j on targets 0, 4, 14, 16 and 39. So other O(ε) terms are missing: D3an-6's offset, the creep-modified pulls εC_j∂_cΔ_kj, and others. "The next order in u is O(ε)" is a numerical observation only. D_j is about 2.7·10² at μ₀ = 0.15 and 1.2·10⁴ at 0.034. On target 0 at 10⁻⁷ the creep term (3.7) rivals the pulls (≤ 5.5), in units of ε^{1/2}. | (d): state ε^{1/2}dev + O(ε^{3/4}), which is what is proved. Present εD_j as the observed dominant part of the O(ε) term. |
| D3an-8 | minor | **The test is weak.** It uses one target (39), which has no tracer changing sign and almost no sliding. Its thresholds are 20× (duals), 15× (crossings) and 2700× (block) above the observed values. **Surviving mutants:** 100× looser loop tolerances (target 39 at 10⁻¹⁰: residual 9.9·10⁻⁵ → 8.3·10⁻⁴, crossings 0.0013 → 0.011); a per-kernel sign instead of per-entry (target 41: 5·10⁻⁴ → 4.2, the very point the paragraph stresses); no Richardson extrapolation (target 41 crossings 0.0145 → 0.34); b^{3/2} → b^{1/2} in κ (untestable at b = 1); the path's rate not zeroed on the walls (harmless). **Killed only narrowly:** the printed D_j (2.8·10⁻³), D_j = 0 (2.7·10⁻³), and Λ₀ halved or doubled in the inner centre equation (2.3·10⁻³ and 3.1·10⁻³), all against a threshold of 2·10⁻³. | `test_transit_limit.py`: add target 41 at 10⁻⁹ with thresholds about 3× the observed values, and check `w_end` and the end centres against S_k. |
| D3an-9 | minor | **"(H3) verified numerically at each entry (the inner problem's end at the path's state)" is not in the committed experiment.** `a_job` ignores `stationary` and `w_end`, and `entry_pull` does not return the end centres. (H3) does hold where I looked: \|w_end − w^k\| ≤ 1.5·10⁻¹⁰ on every entry I integrated, and the loop's block matches S_k to ≤ 7·10⁻¹⁴ at 10⁻¹¹. | Experiment: record max\|w_end−w^k\| and max\|c_end−c^k\| per entry in `transit_limit.md`. |
| D3an-10 | minor | **Implicit hypotheses and definitions.** No exits is only implied through (H2). θ must be a regular value for (c). The inner problem's limit must equal the polished S_k. v_k is \|a_{k−1,j_k}\| in (c) but s_k⟨r_{k−1}φ_{j_k}⟩/ρ in Step 3; these differ by ρ. "The k-th entry arrives" and "the outer intervals" are not defined. The statement allows general masses b_i, but the code is equal-mass only: `path_duals` pins duals at ±1, `creep_constants` uses μ = 1−\|u\|, and `entry_pull` takes a scalar ρ. | Statement: list these. Code: assert b = 1 or implement ρ_i. |
| D3an-11 | minor | **Step 5 and the constant in (b).** Step 5's "starts within O(e^{−λS}) of S_k" omits Step 3's error, which is O(ε^{1/2}) in the stable directions and O(ε^{1/4}) in the neutral ones. The C in (b) is left free although the layers relax slowly. On target 39, at τ_k + 400ε^{1/2} the block is still 1–3·10⁻⁹ off S_k at every ε (λ ≈ 0.05, so C ≳ 30). At τ_k + 20ε^{1/2}log(1/ε) the distance goes from 1.3·10⁻⁷ to 3·10⁻¹¹. | Step 5 and (b): carry the Step 3 error; require C ≥ 3/(2λ_min). |
| D3an-12 | minor | **(b) is tested only at the midpoints**, which probe the interior O(ε^{3/2}), not the stated window. At τ_{k+1} − ε^{1/4} the ε^{3/4} term is attained: 3.9·10⁻⁶, 1.1·10⁻⁶, 2.5·10⁻⁷, 5.1·10⁻⁸ for ε = 10⁻⁷ … 10⁻¹⁰ (target 39). | Experiment: add the window's end points. |
| D3an-13 | minor | **Numerics paragraph details.** "(b) falling 31–32 times per decade … (one target reaching the floor)" holds from 10⁻⁸ to 10⁻¹⁰. In the last decade the factors are 1.5× (target 39), 12.6× (19) and 38× (41), so two targets hit a floor, not one. "Whose path is exact with simple entries and no exits" excludes none of the 48 targets. All other numbers match `transit_limit.md`. | Draft numerics paragraph. |
| D3an-14 | minor | **Novelty.** AGF (Kunin et al., NeurIPS 2025) alternates utility maximisation for dormant units with cost minimisation over active ones, and is proved for diagonal linear networks. It is the closest precedent for "between entries the loop is slaved to the block's fit". D‴.48 cites it; D‴.53 does not. Sliding Frank–Wolfe/ADCG (Denoyelle–Duval–Peyré–Soubies) is the discrete analogue of the sliding path [memory]. Krupa–Szmolyan blow-up is the standard tool for the wall, where the mobility ∝ \|w\|^{4/3} vanishes [memory]. The ε^{1/2} rate of u already appears in D‴.48's numerics. Pesme–Flammarion resolve the jumps by an arc-length reparametrisation (search snippet only). | Classical/New paragraph: add these citations, and scope "new" to the moving-centre limit for this law, its explicit inner problem and the first-order correction. |
| D3an-15 | cosmetic | **Two code nits.** D‴.51's `near_tie_inner` docstring says the xi0 error is ∝ xi0⁻², while the D‴.51 and D‴.53 experiments and the data say ∝ 1/xi0. `loop_with_duals` never checks `sol.status`, so a failed segment would be indexed silently. | `rich_transit.py`. |

## Lesson

A matched-asymptotics theorem is only as useful as its untracked constants are small, and the numerics decide what the reader believes about them. Here the five targets with the largest margins also happen to have fast entrants that barely slide. Among the other 43, eight have an inner-problem tail ∝ g′²/v³ that makes the "first-order correction" 10²–10³ times larger, so that the path is a poor guide down to ε = 10⁻¹¹. Meanwhile the eligibility filter, which should have flagged a narrow selection, excludes nothing. Conversely, a residual can only be read for its rate when the prediction is more accurate than the residual. The Richardson extrapolation at xi0 = 400/800 and rtol 10⁻⁹ was not, and it produced a "floor" that the paragraph blames on the integrator, and a crossing residual that grows where the text says it falls. Converge every predicted constant, and report its convergence, before reading slopes off a residual; and choose test cases by the quantities the constants depend on, not only by the one the hypotheses name.

Sources: [AGF, NeurIPS 2025](https://proceedings.nips.cc/paper_files/paper/2025/hash/06cbd2e81dfbd3bb4cb0abce95b32584-Abstract-Conference.html); [Pesme–Flammarion, NeurIPS 2023](https://proceedings.neurips.cc/paper_files/paper/2023/hash/17a9ab4190289f0e1504bbb98d1d111a-Abstract.html); [Berthier, JMLR 2023](https://jmlr.org/papers/v24/22-1395.html).

Scripts and logs are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review41/scripts` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review41/out`.
