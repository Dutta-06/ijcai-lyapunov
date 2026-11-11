# Referee report D3i (round 16) — PARADIGM §5d: Theorem D‴.21 (statement, proof, Status, *Numerics*), Remark D‴.22; the §0 and §7 rows "| D‴.21 |", the §8 bullet "Done (Theorem D‴.21)" and the branch-8 sentence; Remark D‴.19's last sentence; `results/RESULTS_PARADIGM.md` (Theorem D‴.21); `src/core/ot_basin.py`; `tests/test_ot_basin.py`; `experiments/exp_ot_basin.py` → `results/chaos/ot_basin.md` (+ `ot_basin_raw.npz`); `LITERATURE_TRANSPORT.md` [S47], [S48] and the D‴.21 audit bullet

## Scope and method

**Branch and concurrent edits.**
- The review is of branch `claude/festive-bardeen-eopej2` at c29c34a (HEAD, the D‴.21 commit). Line numbers are c29c34a's.
- During the review a concurrent session had uncommitted edits in the working tree: the D‴.20 block of PARADIGM.md (l. 2381–2424), its audit bullet in LITERATURE_TRANSPORT.md (which also moves D‴.20's stray verdict back, see D3i-6(d)), the entropic-leak files, and new adaptive-prior files.
- None of these touches the material. The D‴.21 block (statement through *Numerics*, 121 lines) is byte-identical in the working tree, shifted by 12 lines. The §0 and §7 rows, the §8 bullet, the branch-8 sentence, RESULTS and the four D‴.21 files are unchanged.

**No repository file was modified.**
- Every script writes only to the scratch directory.
- Python ran with `PYTHONDONTWRITEBYTECODE=1`, and pytest with `-p no:cacheprovider`.
- `--from-raw` ran on a scratch copy of the script with `src` symlinked. It rewrites `ot_basin.md` byte-identically (md5 `693fdcb7…`, 5 s).
- The 25-minute experiment was not rerun. Instead, 34 single systems were re-integrated with the experiment's own integrator; their records are bitwise equal to the stored ones.

**Scratch scripts and their `.out` files**, in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3i/`:
- `my_model.py`: an independent implementation, with no import of `src/`. It covers the saturated features with analytic gradient and Hessian, J^π, ∇J^π, R[μ](ϑ), the weighted stationary atoms, the Hungarian and LP plans, and the 48-point Gauss–Legendre segment integrals of (E) and (P).
- `check_model.py`: `my_model` against finite differences and against the repository.
- `check_identities.py`: (E) and (P) with genuine integral remainders; non-commensurate m.
- `check_PL.py`: the 9-point grid against true segment suprema; the chain of (ii) at 58 configurations.
- `check_sharpness.py`: (v), the whitened Hessians and the PL ratios along directions.
- `recompute_md.py` and `alpha_at_splus.py`: every number of the md, the paragraph, the rows and RESULTS from the raw records.
- `rerun_dense.py`, with outputs `rerun_dense_small.out`, `rerun_dense_240.out`, `rerun_dense_510.out`: the plan at every plan update, α at every record, and the trajectory box.
- `check_apriori.py`, `check_Hstar_global.py`, `check_box_reach.py`, `check_basin_shapes.py`: the a-priori constants, (H*) and the sublevel set.
- `check_spurious.py`, `check_nearest_finals.py`, `check_counts.py`: the nearest-atom contrast and admissibility.
- `check_other_minimum.py`: a local minimum of J^OT other than μ* under (H*).
- `check_moved.py`: the measure-level table and the one-particle runs at a fine step and with true suprema.
- `check_tests.py` and `tests.out`: the tests.
- `derived_numbers.py`, `from_raw.out`, `inspect_raw.py`.

**Source of the numbers.** Every number below comes from these scripts.

## Verdicts

| object | verdict |
|---|---|
| **(i)** identities (E), (P) | **hold.** Re-derived by hand: signs, the (1 − t) and t weights, the stationarity σ₀G⁻¹(ϑ*_j − ϑ⁰_j) = ⟨r*∇Φ̃(ϑ*_j)⟩, and R_{Δf̂} = R[μ*] − R[μ_Θ]. Confirmed with integral remainders to ≤ 3.2·10⁻¹⁵ at 50 configurations, for non-optimal as well as optimal assignments. The diagonal plan's optimality enters only through 𝒥*. For (E) the counts matter only as column marginals, and (E) holds for fractional couplings too. (P) needs an assignment, and at incommensurate m the PL inequality relative to 𝒥* is false (D3i-8(c)). |
| **(ii)** PL inequality | **holds.** The chain was checked with true suprema at 58 configurations. The lower bound on E does not need α > 0 (D3i-8(a)). The constant is sharp only at a matching prior: at the shifted prior it is 0.69 of the local PL constant (D3i-7(e)). |
| **(iii)** a-posteriori decay | **holds for vertex-selection solutions**, whose existence is assumed (§4b; D3i-8(d)). The Dini argument, Grönwall with E < 0, and the distance bounds (which need κ* < σ₀, as stated) are correct. |
| **(iv)** a-priori basin | **holds for solutions that remain in B.** The statement needs B to be forward-invariant; the box of the numerics is not (D3i-3). The uniqueness is of the global minimiser only (D3i-4). The constants α ≥ κ/2 and λ ≥ λ₀ are consistent with (ii), and L₂^G = 11.111 is the exact supremum. At the shifted prior (H*) actually holds on the whole space (κ = 4.7·10⁻⁴). |
| **(v)** matching prior, sharpness | **holds.** The separations attain 2σ₀; the slowest collective mode gives 2λ_coll = 0.01476; no direction does better near μ*. |
| **(vi)** nearest-atom contrast | **holds as stated ("when admissible"), with two overstatements** (D3i-4): "for every count vector N" holds only near q, and "no minimum other than μ* is left" is false for local minima (a strict local minimum of J^OT₃ at the matching prior). |
| Remark D‴.22 | The hand-over inequality max_i‖u_i‖² ≤ mD_π and the time O(λ⁻¹log m) are correct. "Propagation of chaos uniform in time" is overstated (D3i-5). |
| *Status*, novelty, [S47]/[S48] | The plan segments are a generalised geodesic based at ν (AGS §9.2), and PL inequalities through optimal maps have precedents. The reason given for [S48] not applying is wrong, and the "metric slope" wording is imprecise (D3i-6). The bibliographic data are correct. |
| *Numerics*, OT runs | **numbers hold.** No relabelling holds at every plan update, not only at the records. α > 0 holds at every record after s₊. Minor wording issues: D3i-7. |
| *Numerics*, nearest-atom contrast | **numbers hold, reading fails** in 10 of 32 runs, which sit at saddles (D3i-1, **major**). |
| *Numerics*, a-priori basin and measure-level basin | **numbers hold** up to a 1 % underestimate of the box supremum. The "fixed fraction" claims are vacuous at the proven E₀ (D3i-2, D3i-3). |
| §0 row, §7 row, §8 bullet, RESULTS, Remark D‴.19 | Inherit D3i-1, D3i-2, D3i-4 and D3i-7(a). |
| `src/core/ot_basin.py` | **correct.** The identities are in exact-remainder form. The 9-point grid is an upper estimate of λ, by ≤ 2.7 % at the starts checked and 10 % at the moved-particle configuration (disclosed in the docstring). |
| `experiments/exp_ot_basin.py`, `ot_basin.md` | **reproducible** (bitwise re-integration, byte-identical regeneration). Docstring and comment mismatches (D3i-7(g)); no admissibility column (D3i-1). |
| `tests/test_ot_basin.py` | **pass (4/4, 1.32 s); none vacuous; some weaker than their names** (D3i-9). |

## Confirmed issues

**D3i-1 (major; D‴.21 *Numerics* l. 2550–2553 and l. 2559–2560; §0 row l. 115 "the nearest-atom runs stop at the energy of the replicated configuration of their sampled weights"; §7 row l. 4989 "the nearest-atom E stays at …, the spurious minima's energy (within 0.3 % at m ≥ 240)"; RESULTS l. 490–491 and l. 496–497; `ot_basin.md`, the nearest-atom tables).**

The quoted passage reads: "Under the nearest-atom leak E stays at 3.5·10⁻³–1.9·10⁻² …: at or below the excess energy of the replicated configuration of the sampled weights, within 0.3 % at m ≥ 240. That configuration is a critical point of the nearest-atom energy … and not of the OT energy …: the floor is a spurious minimum." Every number reproduces. But at m = 15 and 60 the sampled weights of 5 of the 8 runs per prior are inadmissible: the replicated configuration is a saddle, not a minimum, and the runs are leaving it.

(i) *Admissibility of the sampled counts.*
- The final Voronoi counts equal the initial ones; no particle changes cell in any nearest-atom run (`check_counts.py`).
- For every count vector, the weighted stationary point exists with its atoms in their own cells, and the weighted K-atom Hessian is positive (smallest eigenvalue ≥ 0.00731).
- The class-0 separation margin σ₀ − λ_max R^G is negative for (6, 5, 4), (6, 3, 6), (8, 4, 3)/15 and (24, 15, 21), (25, 17, 18)/60. The margins are −0.00348, −0.00108, −0.00614, −0.00216, −0.00398 at the matching prior and −0.00287, −0.00038, −0.00576, −0.00152, −0.00341 at the shifted prior.
- Only (5, 4, 6)/15, (19, 21, 20)/60, (19, 18, 23)/60 and all eight vectors at m = 240 and 510 are admissible (`check_spurious.py`).
- Under multinomial(15, 1/3) sampling, an inadmissible count vector has probability 0.37 at the matching prior (`derived_numbers.py`).

(ii) *The runs at inadmissible counts are not at that configuration.*
- E(600) lies below the configuration's energy by 8.8, 0.58, 15.8 % (m = 15) and 2.3, 11.0 % (m = 60) at the matching prior, and by 6.4, 0.46, 16.5 % and 2.0, 11.1 % at the shifted prior.
- E is still falling. At the matching prior, (6, 5, 4)/15 has E = 1.143·10⁻², 1.116·10⁻², 1.076·10⁻² at s = 300, 450, 600; (25, 17, 18)/60 has 1.294·10⁻², 1.259·10⁻², 1.240·10⁻².
- In the stored finals class 0 has split. Its rms spread is 0.067–0.62 (matching) and 0.043–0.62 (shifted). The largest distance of a particle from ϑ*_0(N/m) is 0.08–1.19 (`check_nearest_finals.py`).
- At the admissible counts the clouds are collapsed: spreads ≤ 0.006, largest distances ≤ 0.018. There E sits 0.02–0.24 % below the spurious minimum's energy and rises slowly toward it.

(iii) *The one-particle runs.* "E rises to the spurious minima's 5.7·10⁻³ and 5.6·10⁻⁴": these are the values at s = 300.
- The spurious minima of (9, 11, 10)/30 and (99, 101, 100)/300 have E = 6.31·10⁻³ and 6.13·10⁻⁴. Both are admissible.
- At s = 300 both runs are at 0.906 of these values. Continued, they reach 6.304·10⁻³ and 6.128·10⁻⁴ at s = 1000, and 6.3126·10⁻³ and 6.1328·10⁻⁴ at s = 3000 (`check_moved.py`).

*What holds:*
- In all 16 runs at m ≥ 240 and in the 6 admissible runs at m ≤ 60 (22 of 32), the end state is the spurious minimum. It is critical for the nearest-atom energy (drift ≤ 5.6·10⁻¹⁴) and not for the OT energy (drift 2.3–4.4·10⁻²).
- (vi) itself says "when they are admissible" and is not affected.
- This is round 14's D3g-1 failure mode again: inadmissible realised weights have saddle lifts.

*Fix:*
- Have `job()` print, at the sampled counts, the separation margins (`PD.splitting_margin_prior`) and the smallest eigenvalue of the K-atom Hessian, and mark the inadmissible rows.
- Numerics: "at the admissible counts (every run at m ≥ 240, 3 of 8 at m ≤ 60) E sits at the spurious minimum's energy to 0.25 %. At the inadmissible ones (5 of 8 at m ≤ 60) the replicated configuration is a saddle (class-0 margin −0.0004 to −0.0061), class 0 splits, and E ends 0.5–16.5 % below it, still falling."
- §0: "the nearest-atom runs at admissible counts stop at …". §7 and RESULTS likewise.
- One-particle runs: "E rises toward the spurious minimum's 6.3·10⁻³ and 6.1·10⁻⁴ (0.91 of it at s = 300, reached by s ≈ 3000)".

**D3i-2 (minor).**

Locations:
- §0 row l. 115: "the excess energy decays exponentially from an m-independent sublevel set that contains configurations with a fixed fraction of the particles in other cells"; the same in the §8 bullet, l. 5137–5139.
- Remark D‴.19, l. 2320–2322: "the basin at the level of the empirical measure is larger: … proves an m-independent energy basin, which contains configurations with a fixed fraction of the particles in other cells".
- *Numerics* "*The measure-level basin.*", l. 2556; §7 row "+0.0024 with 1/30 of the particles moved into another cell".

The sublevel set {E ≤ E₀} proved by (iv) contains a misplaced fraction only at a size no realistic m resolves.

(i) *How large a fraction fits.*
- (iv) gives E ≥ (κ/2)D_π. A particle outside the cell of its class atom is at G⁻¹-distance ≥ 3.1315 (matching) or ≥ 2.9471 (shifted) from that atom.
- So E ≥ 0.0347p (matching, κ = σ₀) and E ≥ 0.0118p (shifted, κ = 0.0027) for a misplaced fraction p.
- With E₀ = 2.03·10⁻⁷ and 3.0·10⁻⁸, the set {E ≤ E₀} holds p ≤ 5.8·10⁻⁶ and 2.5·10⁻⁶. A single misplaced particle therefore needs m ≥ 1.7·10⁵ and 3.9·10⁵ (`check_apriori.py`).
- Directly: one particle just across the cell-0/1 boundary has E = 12 336 E₀ at m = 30 and 803 E₀ at m = 300.
- The numerics' p = 1/30 configuration has E = 3.35·10⁻³ = 16 500 E₀.

(ii) *The warm-cloud starts.* E(0) ≈ 0.097/√m (matching) and 0.095/√m (shifted). A warm i.i.d. start therefore lies in {E ≤ E₀} only for m ≈ 2·10¹¹ and 10¹³.

(iii) *"Larger".* The two sets are not nested.
- Theorem D‴.18's per-particle ball (proof radius ≤ 0.019, round 14) is not in {E ≤ E₀}. A separation moving every particle by 0.019 (respectively ρ/2 = 0.0095) has E = 6.30 E₀ (1.58 E₀) at m = 30 and at m = 300.
- Conversely, {E ≤ E₀} allows a single particle to move by up to 0.042 (m = 30) or 0.131 (m = 300) (`check_basin_shapes.py`).
- The comparison is of parameter sets only: D‴.18 is about the real loop and D‴.21 about the averaged loop.

*What holds:*
- The p = 1/30 example is covered a posteriori by (iii). Along the computed run α stays positive: λ ≥ 0.00148 at every 0.5 slow units with true suprema. So the bound E(s) ≤ E(0)exp(−∫λ₊) applies along that path.
- It is not covered by the a-priori basin (iv).

*Fix:*
- Replace the three sentences with: "an m-independent sublevel set (on the instance, E₀ = 2·10⁻⁷ holds a fraction p ≤ 6·10⁻⁶ of the particles in other cells, so none below m ≈ 2·10⁵); configurations with 1/30 of the particles in another cell are covered a posteriori by (iii), α staying positive along the run."
- In Remark D‴.19, "is larger" → "is of a different shape".

**D3i-3 (minor; (iv) l. 2476–2482 and proof l. 2508–2511; *Numerics* l. 2553–2556; `exp_ot_basin.py` docstring l. 17–18 and `a_priori_constants`).**

(a) *Invariance of B.*
- The rate part of (iv) holds for solutions that stay in B^m. The proof's "E is non-increasing, so this persists" keeps α ≥ κ/2 only while Θ(s) ∈ B^m, since κ* ≤ σ₀ − κ needs the segments in B.
- The example given, N_{R_Θ}(conv supp ν), is forward-invariant: T̄_i ∈ conv supp ν, and Step 0's projection argument gives D⁺dist ≤ L_ϑ^G(F_ρ + K_Φ) − σ₀ dist.
- But (iv) does not require invariance, and the numerics apply it to a box that is not invariant. The box's faces are 0.227 (G⁻¹) from atom 1. One particle placed just outside the box, with all others at the atoms, has E = 2.12, 1.06, 0.85, 0.21 E₀ at m = 3000, 6000, 7500, 30 000 (`check_box_reach.py`). So for m ≳ 6400 the set {E ≤ E₀} reaches outside the box, and nothing keeps a solution inside it.
- In the runs themselves the trajectories do stay inside: the whole trajectories of all 34 re-integrated systems lie in the box.

(b) *The box supremum.*
- Multi-start L-BFGS-B over the box gives sup λ_max R^G[μ*] = 4.390·10⁻³, against the reported 4.349·10⁻³ from 20 000 random points (0.95 % low). A 61×61×31 grid gives 4.371·10⁻³.
- The supremum is attained on the box's boundary (α and κ at their lower faces), so it depends on the 10 % padding.
- The corrected values are κ = 0.00268 and E₀ = 2.92·10⁻⁸.

(c) *(H*) holds on the whole parameter space at the shifted prior.*
- A grid over the full saturation ranges (484 169 points) plus Nelder–Mead gives sup over ℝ³ of λ_max R^G[μ*] = 6.601·10⁻³ < σ₀ = 7.074·10⁻³.
- The supremum is approached at the saturation corner a → a_max, s → s_min, c = 1.82, where only the centre–centre entry survives; a 1-D maximisation there confirms it (`check_Hstar_global.py`).
- So κ_global = 4.73·10⁻⁴. Then μ* is the unique minimiser among all m-point empirical measures, not only those in the box, and the invariance question disappears (E₀ = 9.1·10⁻¹⁰).

*Fix:*
- In (iv), "a convex set B, forward-invariant for (𝔄_m^{OT}), …" (or "every vertex-selection solution that remains in B").
- In the numerics, report (H*) on the whole space (κ = 4.7·10⁻⁴, E₀ = 9·10⁻¹⁰), or the optimised box value with the invariance caveat.
- Docstring: the box is that of the end configurations, priors and atoms (the trajectories happen to lie inside), not "a box containing every trajectory".

**D3i-4 (minor; (vi) l. 2487–2492 and proof l. 2513–2515; §0 l. 115 and §8 l. 5139–5140 "the unique minimum among the empirical measures").**

(a) *"Under (H*) no minimum other than μ* is left" (l. 2492).*
- (iv) proves that μ* is the unique *global* minimiser among the m-point empirical measures in B. (ii) excludes other critical points only where α > 0. A critical point with E > 0 needs κ^Θ > 2σ₀ on its segments, which (H*) does not prevent.
- Counterexample at the matching prior, where (H*) holds on the whole space with κ = σ₀, at m = 3 (commensurate):
  - the second local minimum of 𝒢_q (round 13), recomputed independently: G_q = 0.107139, a = (3.148, −3.176, −1.417), c = (−0.445, −0.135, 0.169), s = (1.409, 0.497, 0.520), fit error 0.366;
  - its diagonal plan is the unique optimal plan (cost gap 29.1 to the next permutation), and it is critical for J^{OT}_3 (|∇| = 9.3·10⁻¹⁰);
  - the whitened Hessian of J^{OT}_3 there has smallest eigenvalue 0.00861 > 0, so it is a strict local minimum with E = 0.107 ≈ 5·10⁵ E₀;
  - α = −0.0210 there (κ^Θ = 0.0561 > 2σ₀), as (ii) requires of any critical point with E > 0 (`check_other_minimum.py`).
- At m = 15 the same configuration replicates to a saddle (round 13), so the claim may hold for larger m, but that is not proved.

(b) *"for every count vector N".*
- The weighted stationary point exists by the implicit function theorem on a neighbourhood of q (Proposition D‴.5(i)).
- At m = 15, continuation from q finds it, with its atoms in their own cells, for 90 of the 91 count vectors with no empty class. It is not found for (1, 13, 1). Of the 90, 55 are admissible and 35 are saddles (`check_spurious.py`).
- Count vectors with an empty class are outside Proposition D‴.5.

*Fix:*
- (iv): "… μ* is the unique global minimiser".
- (vi): "under (H*) μ* is the unique global minimiser in B, and the only critical point where α > 0 (in particular in {E ≤ E₀})"; "for every count vector N in the domain of w ↦ Θ*(w) (a neighbourhood of q) whose atoms lie in their cells".
- §0 and §8: "the unique global minimum".

**D3i-5 (minor; Remark D‴.22, l. 2517–2526).** The remark reads: "in the basin both the particle system and the mean-field flow converge to μ* exponentially in W₂, at the same m-independent rate: W₂(μ_{Θ(s)}, μ_s) ≤ W₂(μ_{Θ(s)}, μ*) + W₂(μ_s, μ*) → 0 — propagation of chaos uniform in time with no floor".
- The displayed bound shows that the long-time limits coincide. At any fixed s it is not small in m, so it does not give propagation of chaos, which requires sup_s E W₂(μ_{Θ(s)}, μ_s) → 0 as m → ∞.
- At s = 0 the bound is W₂(μ_{Θ(0)}, μ*) + W₂(μ_0, μ*). For the warm clouds W²(μ_{Θ(0)}, μ*) = 0.99–5.25 (means over the seeds), and an i.i.d. warm start lies in the a-priori basin only for m ≈ 10¹¹ (D3i-2).
- The remark itself disclaims "smallness in m at finite times — the transit itself". That smallness is exactly what uniform-in-time propagation of chaos needs.
- Genuine uniform-in-time propagation of chaos in mean-field training is known for entropy-regularised (mean-field Langevin) dynamics under a uniform log-Sobolev inequality: F. Chen, Z. Ren, S. Wang, arXiv:2212.03050 (2022); T. Suzuki, A. Nitanda, D. Wu, ICLR 2023 (an O(1/N) error uniform in time). These are the natural comparison for the remark and for §8.

*Fix:* "… → 0: no floor. For every commensurate m the particle system and the mean-field flow have the same limit μ*. Uniform-in-time propagation of chaos would also need smallness in m during the transit, which is not claimed."

The rest of the remark checks out:
- max_i‖u_i‖² ≤ mD_π ≤ 2mE/(σ₀ − κ*) gives the hand-over condition.
- That condition matches D‴.18(b)'s start set: the per-particle ball of radius ρ/2 around Θ*_m(q), with the plan as class map.
- The time to reach it is λ₀⁻¹log(2mE(0)/((σ₀ − κ*)(ρ/2)²)) = O(λ⁻¹log m).

**D3i-6 (minor; novelty and literature: §0 novelty column l. 115; *Status* l. 2531–2536; LITERATURE l. 328–338).**

(a) *The plan segments are a generalised geodesic based at the prior.*
- With the diagonal plan ν → μ* (optimal by hypothesis) and the assignment plan ν → μ_Θ, the counts make "mass 1/m on (ϑ⁰_{π(i)}, ϑ*_{π(i)}, ϑ_i)" a three-plan with these two-marginals. Its interpolation (1 − t)ϑ*_{π(i)} + tϑ_i is the curve of the plan segments.
- That curve is the generalised geodesic of Ambrosio–Gigli–Savaré §9.2 (Definition 9.2.2), along which ½W²(·, ν) is 1-convex. This is where the (σ₀/2)D_π of (E) comes from.
- Passing from convexity along a curve that ends at the point to a slope (PL) bound — the derivative at the endpoint divided by the L²(μ_Θ)-norm √D_π of the curve's velocity there — is the standard argument (AGS Ch. 2.4 and 4).
- So "the plan as the test direction" is AGS's device. Closer precedents for PL inequalities in Wasserstein space obtained through optimal maps, for functionals that are not geodesically convex:
  - S. Chewi, T. Maunu, P. Rigollet, A. J. Stromme, *Gradient descent algorithms for Bures–Wasserstein barycenters*, COLT 2020 (a PL inequality for the barycentre functional);
  - Q. Mérigot, F. Santambrogio, C. Sarrazin, *Non-asymptotic convergence bounds for Wasserstein approximation using point clouds*, NeurIPS 2021 (a modified PL inequality for ½W²(μ_Y, ρ) over uniform point clouds — the Lloyd/power-cell, capacity-constrained case — with an error term set by the inter-atom distances).
- I found no statement of a PL inequality for this energy (a mean-field loss plus a W² penalty to a discrete prior). What is specific here:
  - the flat-convex loss, linearised exactly in f̂ with the curvature taken at the stationary residual;
  - the explicit local constants;
  - the observation that commensurability makes the generalised geodesic a particle interpolation, with constants free of m.

(b) *[S48].*
- Blanchet–Bolte also treat λ-convex functionals (their abstract).
- L̄ is semiconvex along generalised geodesics, with constant −L₂^G·sup⟨|r|⟩ (the saturated features are bounded). Hence 𝒥 is (σ₀ − K)-convex along generalised geodesics based at ν.
- Here K > σ₀ globally, and D‴.21 is the local form of this structure.
- So "flat-convex but not displacement-convex, so [S48] does not apply" gives the wrong reason. The right one is "not semiconvex with a positive constant globally".

(c) *The leak alone.*
- The metric slope of ½W²(·, ν) equals W(μ, ν) at every μ: the geodesic toward ν gives ≥, and the triangle inequality gives ≤. So "its metric slope being W(·, ν) where the optimal plan is a map" puts the qualifier in the wrong place.
- What needs a map is that a particle velocity field realises the slope. At a split plan the barycentric gradient is shorter; the m = 16 stationary configuration of D3i-8(c) has ‖∇J‖ = 0 with W > 0.
- "[S47]" follows an AGS statement that it does not support.

(d) *Double verdict.* At c29c34a the D‴.21 audit bullet was inserted between D‴.20's text and D‴.20's verdict. D‴.21 therefore carries two verdicts (l. 337–338) and D‴.20 none. The concurrent session's uncommitted edit already moves the verdict back.

The bibliographic data are correct as far as checked: Blanchet–Bolte, J. Funct. Anal. 275(7):1650–1673, 2018 (one listing says vol. 25; 275 is right). Polyak 1963 and Karimi–Nutini–Schmidt (ECML PKDD 2016) were not web-checked.

*Fix:*
- Cite AGS §9.2 and the two precedents.
- Recast the novelty as: "the PL inequality along the generalised geodesic based at the prior, with the loss linearised at the stationary residual; the counts make it a particle interpolation and the constant m-independent".
- Correct (b), (c) and (d).

**D3i-7 (cosmetic; D‴.21 *Numerics*, §7 row, md, experiment script).**

(a) *"the smallest ratio/λ is … 1.58–1.72 at the shifted one" (also RESULTS l. 486).*
- The per-size minima of ratio/λ are 1.6064, 1.6738, 1.6753, 1.7237, i.e. 1.61–1.72. The 1.58 is the minimum of Q²/(D_πE)/λ.
- The §7 row's "≥ 1.58 (shifted)" likewise mixes the two.
- At the matching prior the unrounded minimum of ratio/λ is 1.00035, so "at least λ" holds there.

(b) *"at least λ at every record".* The report checks the 25 segment records, where λ is computed. My re-integration of 26 OT systems computes α and λ at every 5-unit record. It confirms the claim: the smallest ratio/λ after s₊ is 1.0011 (matching) and 1.674 (shifted), and α > 0 at every record after s₊.

(c) *"W²(μ_Θ, μ*) from 1.0–5.2"* should read 0.99–5.25.

(d) *λ from 9-point grids.*
- The grids give a lower estimate of the segment suprema, hence an upper estimate of λ.
- Against the true suprema, the overestimate is at most +2.7 % at the 26 starts checked (shifted prior, m = 60, seed 0) and 0 at the 32 finals.
- At the moved-particle configuration it is +10 %: λ = 0.00148 and α = +0.0023, against the md's 0.00163 and +0.0024.
- The "ratio ≥ λ" checks are conservative in this direction. "+0.0024 at p = 1/30" should read +0.0023.

(e) *"the curvature along the segments … costs a factor of about 0.55".*
- The factor is λ/(2σ₀) = 0.552, which splits as 0.802 × 0.688.
- 0.802 = 2λ_sep/(2σ₀) is the true cost of the curvature: the local PL constant 2λ_min of the whitened Hessian at μ* is 2λ_sep = 0.01135.
- 0.688 = λ/(2λ_sep) is the slack of (ii), which bounds E with σ₀ + κ*_− and Q with (κ* + κ^Θ)/2 separately.
- So (ii) is sharp only at a matching prior; say so.

(f) *"identically at m = 30 and m = 300, the margin being a function of the empirical measure".* The two configurations are the same measure, so the identity holds by construction.

(g) *`exp_ot_basin.py`.*
- Docstring l. 17–18 says "over a box containing every trajectory"; the code boxes the end configurations, the priors and the atoms.
- The comment at l. 271 says "λ piecewise constant between the segment records"; the code interpolates linearly (`np.interp`) and integrates by the trapezoid rule.
- The docstring at l. 175, "λ held at the smaller endpoint over every 5 slow units" (with n_t = 5), is not a rigorous evaluation of ∫λ₊. It is harmless here: sampling λ every 0.5 slow units with true suprema gives bounds 9.90·10⁻⁵ (m = 30) and 3.69·10⁻⁶ (m = 300), against the md's 1.017·10⁻⁴ and 3.705·10⁻⁶, with E(300) = 2.546·10⁻⁵ and 3.093·10⁻⁶.
- The midpoint step ds = 0.5 is accurate: ds = 0.05 gives the same E(300) to five digits.

**D3i-8 (cosmetic; the statement and the proof).**

(a) *The lower bound on E in (ii)* (l. 2463–2464): "If π is optimal at Θ and the plan-segment margin α … is positive, then ½‖Δf̂‖² + ((σ₀ − κ*)/2)D_π ≤ E".
- This bound follows from (E) for every assignment with the counts, with no condition on α and without κ* ≤ σ₀. It is informative (it gives E ≥ 0) only when κ* ≤ σ₀.
- α > 0 does not imply κ* < σ₀, since κ^Θ can be negative.
- *Fix:* state the bound before the hypothesis.

(b) *"every optimal plan is then an assignment"* (l. 2445–2446). Only every *vertex* optimal plan is; at a tie, convex combinations are optimal too. J_m^{OT} = min_π J^π_m over assignments is correct.

(c) *"The counts enter through Δf̂ …"* (l. 2461).
- They also enter the constant term of the leak.
- (E) in fact holds for every coupling of μ_Θ with the masses q, fractional plans included: in coupling form the residual is 3·10⁻¹⁶ at m = 16 and 1·10⁻¹⁵ at m = 17. For (E) the counts matter only as the column marginals.
- With an assignment of counts N ≠ mq, (E) fails by exactly Σ_k(N_k/m − q_k)c_k, where c_k = −⟨r*, Φ̃(ϑ*_k)⟩ + (σ₀/2)‖ϑ*_k − ϑ⁰_k‖². At the shifted prior c = (0.85, 3.38, 7.94)·10⁻⁴; at a matching prior c = 0.
- With such an assignment, (P) fails by exactly −⟨Δf̂, f̂[μ_N] − f̂[μ*]⟩, a relative error of 0.75–1.2.
- So it is (P), with one test direction per particle, that needs the assignment. At incommensurate m, (P) with the barycentric test direction fails by 52–83 %.
- The OT flow at m = 16 converges to a configuration with two split particles, where ‖∇J‖²_G = 4·10⁻²⁷ (matching) and 5·10⁻²⁴ (shifted) while E = 5.4·10⁻³ and 5.9·10⁻³. So no PL inequality relative to 𝒥* holds there (`check_identities.py`).
- *Fix:* say this in one sentence.

(d) *Vertex-selection solutions.* (iii), (iv) and the *Status* say "proved for vertex-selection solutions", but their existence is assumed (§4b: "their existence for all time is assumed"). Say so, as Theorem C′'s §0 row does. For this leak a two-particle tie is never attracting:
- with h := ⟨ϑ_i − ϑ_{i′}, ϑ⁰_k − ϑ⁰_j⟩_{G⁻¹}, one has J^{π₁} − J^{π₂} = σ₀h;
- on the two sides of the tie dh/ds = ℓ ∓ σ₀Δ², with the same loss term ℓ;
- so no sliding motion is possible. Multi-way ties are not covered.

**D3i-9 (cosmetic; `tests/test_ot_basin.py`).** All 4 tests pass in 1.32 s, and none is vacuous, but three are weaker than their names.
- *Test 1.*
  - At the matching prior ‖r*‖ = 9.5·10⁻¹⁷, so (E) reduces identically to E = ½‖Δf̂‖² + (σ₀/2)D and the check is trivial there. Only the shifted prior (‖r*‖ = 6.3·10⁻³) exercises the identities. The residuals are ≤ 1.9·10⁻¹⁵ against a tolerance of 10⁻¹².
  - `BasinDiagnostics.terms` asserts counts = mq, so the role of the counts, which the test's name refers to, is never exercised.
- *Test 2.* Part 2a is sharp and meaningful: ratio/(2σ₀) − 1 = 1.2·10⁻⁶ and λ/(2σ₀) − 1 = 3.5·10⁻⁹, against a tolerance of 10⁻⁴. In part 2b the chain holds with wide margins: ratio/(Q²/DE) = 2.3–5.5 and (Q²/DE)/λ = 1.5–7.2.
- *Test 3.* The docstring's "the margin is positive" is not asserted (α = +0.0024 at the start). The bound is asserted with a 1.01 slack, where the actual margin is a factor 4.
- *Test 4.* "admissible weights: a strict minimum" asserts only the separation margins (smallest 0.00087). The K-atom Hessian is not checked (its smallest eigenvalue is 0.00734).

*Fix:* add, at the shifted prior, an assertion with counts N ≠ mq against the two mismatch formulas of D3i-8(c); assert α > 0 in test 3; check the K-atom Hessian in test 4.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_ot_basin.py --from-raw` on a scratch copy rewrites `ot_basin.md` byte-identically (md5 `693fdcb7…`).
2. **An independent implementation.** `my_model.py` uses analytic derivatives of the saturated features. They match finite differences to 4·10⁻¹⁰ (gradient) and 8·10⁻¹⁰ (Hessian), and ∇J^π matches to 6·10⁻¹¹. It reproduces the repository's stationary atoms to 4·10⁻¹⁶, 𝒥* (4.6447085·10⁻⁴ at the shifted prior) and J^π to 12 digits.
3. **(E) and (P) with genuine segment integrals** (48-point Gauss–Legendre on the analytic Hessian).
   - Relative residual ≤ 3.2·10⁻¹⁵ over 50 configurations: m = 12, 30, 60; random (non-optimal) and optimal assignments; displacement amplitudes 0.1–2.5 (D_π up to 25.7); six particles moved across cells; both priors.
   - The hand derivation agrees in every sign and factor.
   - The diagonal plan's optimality enters only through 𝒥*. If the diagonal plan were not optimal, (E) would hold with the diagonal energy in place of 𝒥*. (At both priors the diagonal plan is optimal, so this hypothesis is not exercised.)
4. **(ii).** The chain ‖∇J‖²/(mE) ≥ Q²/(D_πE) ≥ λ and the bounds ½‖Δf̂‖² + ((σ₀ − κ*)/2)D ≤ E ≤ ½‖Δf̂‖² + ((σ₀ + κ*_−)/2)D and Q ≥ ‖Δf̂‖² + αD hold, with the true segment suprema, at all 58 configurations checked. These are the 26 initial configurations rebuilt from the seeds (their E(0) and segment curvatures equal the stored records) and the 32 OT finals.
5. **(iii) and (iv), re-derived.**
   - The Dini argument holds, with equality at points of differentiability.
   - Grönwall works when E < 0: E is non-increasing, so E(s) ≤ E(0) ≤ E(0)e^{−∫λ₊}.
   - The distance bounds and ‖R^G[μ_Θ] − R^G[μ*]‖ ≤ L₂^G‖Δf̂‖_{L¹} hold.
   - α ≥ κ − L₂^G(E/2)^{1/2} ≥ κ/2 for E ≤ E₀, and λ ≥ λ₀.
   - The equality case of the uniqueness is exactly D_π = 0.
   - L₂^G = 11.111111 is the exact supremum: 400 Nelder–Mead starts reach G_c a_max/s_min² at the saturation corner.
6. **(v).**
   - At the matching prior the whitened Hessian of J^π at Θ*_m (m = 6, 30) has separation eigenvalue exactly σ₀ = 0.007074 and collective spectrum 0.00738–0.2023.
   - The PL ratio tends to 2σ₀ along separations (0.014149 at rms 10⁻³) and to 2λ_coll = 0.01476 along the slowest collective mode. Along mixed and random collective directions it is ≥ 0.154.
   - The theorem's λ tends to 2σ₀, with ratio/λ = 1.0000: sharp.
   - At the shifted prior the separation spectrum is 0.005674–0.00823 and the collective spectrum 0.007398–0.2046.
7. **Every number of the md, the paragraph, the §0/§7 rows and RESULTS**, recomputed from the raw records:
   - identity residuals 5.4·10⁻¹⁵ / 4.1·10⁻¹⁵ (matching) and 5.6·10⁻¹³ / 4.2·10⁻¹⁵ (shifted);
   - no relabelling at the records in any of the 32 OT runs;
   - s₊ ≤ 300, and s₊ ≤ 25 in 19 of 32 runs (13 at 0); α(s₊) ≥ 1.6·10⁻⁴;
   - E(0) 4.11·10⁻³–3.98·10⁻², E(600) 8.29·10⁻⁷–8.08·10⁻⁶;
   - W² 0.993–5.254 → 2.34·10⁻⁴–2.28·10⁻³;
   - the a-posteriori bound lies above E(600) in all 32 runs, and E never increases between records;
   - the late ratio 0.014163–0.014173 equals −d log E/ds of the records to 10⁻⁶;
   - the late λ 0.013962–0.014056 lies within 0.76–1.49 % of the observed ratio;
   - at the shifted prior, λ 0.00766–0.00779, observed 0.01367–0.01389, κ* and κ^Θ 1.40–1.48·10⁻³;
   - nearest-atom E(600) 3.46·10⁻³–1.86·10⁻², W² 1.01–4.72, within 0.25 % of the spurious minimum at m ≥ 240 and always below;
   - drifts ≤ 5.6·10⁻¹⁴ (nearest-atom) and 2.33–4.40·10⁻² (OT);
   - E₀ = 2.03·10⁻⁷ and 3.01·10⁻⁸.
8. **Dense re-integration** with the experiment's integrator (records bitwise equal). Covered: 26 OT systems (m = 15, 60, 240 for all four seeds at both priors, and m = 510 for seed 0 at both priors), plus the 8 shifted-prior nearest-atom systems at m = 15, 60.
   - No plan change at any of the 3000 plan updates of an OT run, not only at the records.
   - α > 0 at every 5-unit record after s₊, and positive from at most s₊ (for example from 145 against s₊ = 150).
   - Every trajectory, of all 34 systems, lies inside the (H*) box.
9. **The measure-level table, recomputed with true suprema.**
   - κ^Θ = 0.00957, 0.0191, 0.0287, 0.0479 (the grid gives 0.00934, 0.0187, 0.0280, 0.0467), and α = +0.0023, −0.0025, −0.0073, −0.0169.
   - The moved particles are in cell 1 and kept in class 0.
   - The one-particle OT runs reproduce E(300) = 2.546·10⁻⁵ and 3.093·10⁻⁶. The particle is back in cell 0, the ratio at the end is 0.01415 = 2σ₀, E is non-increasing, and λ ≥ 0.00148 throughout.
   - The nearest-atom runs reproduce 5.72·10⁻³ and 5.56·10⁻⁴, with the particle in cell 1.
10. **The spurious configurations.**
    - Newton residuals ≤ 2.6·10⁻¹⁶ for all 19 count vectors used, with every atom in its own cell.
    - E_spur equals the stored values to 5 digits and is positive in every case, as (iv) requires: (H*) holds globally at both priors.
    - (iv) applies to these configurations with their own OT plan, which reassigns ½‖N − mq‖₁ particles.
11. **Tests.** 4 passed in 1.32 s. README's "108 tests" equals the number of test functions at c29c34a.
12. **Literature.** Abstracts checked by web search: Blanchet–Bolte (J. Funct. Anal. 2018: λ-convex and displacement-convex Łojasiewicz/Talagrand equivalences), Mérigot–Santambrogio–Sarrazin (NeurIPS 2021), Chewi–Maunu–Rigollet–Stromme (COLT 2020), Suzuki–Nitanda–Wu (ICLR 2023).

## Not checked

- **The experiment as a whole.** The 25-minute run was not repeated. Six of the eight m = 510 OT systems were not re-integrated at plan-update resolution.
- **Relabellings inside a plan-update interval** (10 steps) are not checkable: the dynamics freezes the plan there by construction.
- **Solution concepts.** The existence of vertex-selection solutions (assumed by the authors) and multi-way ties were not analysed.
- **Formal parts of Remark D‴.22** were not checked beyond the inequality and the time estimate: the existence of the mean-field (Laguerre-cell) OT flow, its chain rule, and the transfer to the real loop.
- **Citations from memory.** The AGS section and lemma numbers (§9.2, Definition 9.2.2; Ch. 2.4 and 4; §10.4 as cited by the authors) and the Polyak and Karimi et al. data come from memory or secondary sources.
- **Context results.** Propositions D‴.5, D‴.9, D‴.12, D‴.17 and Theorem D‴.18 were used as given (rounds 9–14).
- **Concurrent edits.** The concurrent session's uncommitted edits (the D‴.20 block, the entropic-leak and adaptive-prior files) were not reviewed.
