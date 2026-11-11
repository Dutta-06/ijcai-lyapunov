# Referee report D3aj (round 37): Proposition D‴.48 (the rich transit), objects at 5995696

## What I checked and how
- **Setup.** I copied the repository to SCRATCH/repo at commit 5995696 and worked only there. Nothing in the working repository was
  touched; I ran only read-only git commands against it. Its HEAD has since moved to 1dc7594 through other sessions. Those commits add
  `K.filt` and `moving_rhs_widths` but leave the logic of `sliding_iss`, `moving_rhs`, `wfr_inactive_rhs` and `iss_path` unchanged, so
  every finding below still applies at HEAD.
- **Tests.** The D48 tests pass (6 of 6). The full tests/test_mirror_free.py passes (18 of 18).
- **Reproduction.** I re-ran `exp_mirror_transit.py` with MT_OUT set to SCRATCH/out (176 s). The regenerated md is byte-identical to the
  committed `results/chaos/mirror_transit.md`. I checked every number in numerics_D48.md against the md and against the raw npz.
- **Mutation testing.** 53 mutants, run through the harness `work/mutate.py` with a 60 s per-test alarm. For each surviving mutant I
  measured its effect on the experiment's own paths (`work/mut_effect*.py`).
- **Independent numerics** (scripts and logs in the round's scratch directory, `work/`):
  - iss_path's defining properties, and the loop against the path down to ε = 10⁻¹², including a Q-preserving integrator: `A_check.py`,
    `A5*.py`.
  - The instantaneous creep speed against both formulas: `C6_check.py`.
  - Counterexamples (ties, near ties, close and off-grid bumps, coarse grids, edges, amplitudes): `cex.py`, `close_diag.py`, `edge_*.py`.
  - The WFR limit along trajectories, with ρ ≠ 1, and the formed-kernel count down to ε = 10⁻¹⁰: `wfr_traj.py`, `wfr_small.py`,
    `wfr_rand.py`.
  - The single-unit transit scale for α(0) = θε^p: `o_check.py`.
  - The sensitivity of `sliding_iss` to its parameters: `slide_sens.py`, `slide_tf.py`, `seed_dep.py`, `wider.py`.
- **Literature.** I used web search only; arXiv and author homepages could not be fetched.
  - Berthier–Pillaud-Vivien, arXiv:2606.23198 exists. Its abstract matches lit_S75: a general convex lsc potential started near its
    domain's boundary, a limiting flow whose potential is the domain's indicator, and an explicit link to the inverse scale space flow.
  - The statements of Rockafellar–Wets Thm 7.17, Attouch, Brezis–Pazy and Pesme–Flammarion are **from memory**.

## Verdicts

**(o)** The expansions are correct: w = (2c/(ma_max²))α³(1+O(α²)) and α̇ = (2Gc/a_max²)α²g(1+O(α²)). The "dichotomy … set by its start's
amplitude" is overstated. I tested one unit (Φ = [1], y = 1, mass ε) and measured the time T for w to reach ½:

| start | regime | measured |
|---|---|---|
| α(0) = 0 | power scale | T·ε^{1/2} → 7.03 |
| α(0) = 0.4 | log scale | T/log(1/ε) = 21.5 → 24.7 |
| α(0) = 0.4ε^{1/4} (on a branch for ε ≤ 10⁻⁶) | intermediate, ε^{−1/4}log(1/ε) | T·ε^{1/4}/log(1/ε) = 6.7, 7.9, 8.8, 9.4 (ε = 10⁻⁶ … 10⁻¹²) |
| α(0) = 10ε^{1/2} (on a branch) | power scale | T·ε^{1/2} = 3.23 |

The last row contradicts "on a branch … the depth-2 behaviour". What selects the regime is α(0)/ε^{1/2}, and there is a continuum of
scales in between. The weights 1/(γ|a(0)|) in (ii) are not uniform in the start.

**(i)** The proof is essentially sound.
- Inside the box: g_ε ≤ (εb/(γmΛ))cot²(πδ/2)(1+O(ε)) = O(ε^{3/2}δ^{−2}) uniformly on |v| ≤ (1−δ)ρ.
- Beyond the wall: α ≥ 2Gc(η−λ₁) holds exactly (arctan ≤ π/2), so the "−O(1)" is unnecessary and the exponential lower bound is uniform.
- RW 7.17 is applicable, but not by the 1-D statement alone. Compacts of ℝᴺ∖P can meet the hyperplanes (Φᵀζ)_j = ±ρ_j, where g_ε,j does
  not converge uniformly. The missing step uses g ≥ 0 together with max_i|(Φᵀζ)_i|/ρ_i ≥ 1+δ on such a compact, or the dense-set
  criterion.
- Attouch and Brezis–Pazy are used correctly (0 ∈ int P). The minimal-norm section is indeed the least-squares element of 𝒩(u).
- The energy identity has the wrong factor: dζ/dτ = −∇𝒟_ε/N gives N∫|ζ̇|², not N⁻¹.
- The path's structure checks out. On designs 3–5, max(|u|/ρ) − 1 ≤ 1.5·10⁻¹⁴, w lives only on the wall, |v| on the support ≤ 3·10⁻¹⁶,
  τ₁ = min ρN/|Φᵀy| holds to machine precision, and exits occur only at entry times.
- "w^ε → w in L² where the active columns are independent" is right on compact subintervals between events. It needs uniform margins of
  the inactive atoms, which the proof does not state.
- The ε^{1/2} rate is confirmed and is cleaner than reported. With a Q-preserving integrator the loop follows the path to ε = 10⁻¹²:
  sup|u−u_ISS| = 6.7·10⁻⁶, 9.0·10⁻⁶, 2.0·10⁻⁵ on designs 3, 4, 5. Design 5's slope is 0.4998; the reported 0.471 is an artefact (D3aj-8).

**(ii)** The proof route is fine, and sup|u−u_ISS| is exactly ∝ 1/log(1/ε). The sign constraint does bind in Part B: the two-sided path
would enter atom 9 at its wrong wall at τ = 68.5, so B genuinely tests the one-sided box. The unit test's one-sided case, however, is
identical to the two-sided path. The 1/log rate leaves the limit practically invisible:
- at ε = 10⁻⁴ the error is 40–80% of the weights (which are 18.8–40.5, not "17–50"), and the entry-time error (55) exceeds the first entry
  (21.8);
- at ε = 10⁻¹⁶ the entry-time error is still 14;
- for ± pairs the first step is a factor 2.65, not 2.

Given D‴.47(i), (ii) is Berthier–Pillaud-Vivien's entropic case and is not new.

**(iii)(a)** The coefficient bound is correct: the exact w is 0.25–0.26 of the μ^{−3} bound, and the cos/sin³ profile matches the curve to
< 2% for μ ≫ √ε. The creep speed is wrong. The ODE's dc/dτ equals (πa_max²/(4c))ε b_i^{3/2}·cos(πμ/2)/sin³(πμ/2)·|⟨r∂_cφ⟩| to 3–4 digits:

| case | ODE | stated formula | exact profile |
|---|---|---|---|
| kernel 15, ε = 10⁻⁸, τ = 17.70 | 9.747·10⁻³ | 3.792·10⁻² | 9.785·10⁻³ |
| ε = 10⁻¹⁰, τ = 20.50 | 3.212·10⁻⁴ | 1.245·10⁻³ | 3.212·10⁻⁴ |

The stated formula is too large by π³/8 ≈ 3.9, and it is missing b_i^{3/2} and sign(u_i).

**(iii)(b)** Formal, as labelled. With t_fast ≥ 10⁶ the relaxation was seed-independent on every case I tried (seeds 10⁻³ to 10⁻⁹). The
claim "k well-separated bumps are formed by k kernels at their centres" fails at ties:
- **Exact tie.** One bump at 0.05, midway between grid kernels: the loop forms **two** kernels at 0.05 with 0.40 each at ε = 10⁻⁵, 10⁻⁶
  and 10⁻⁸. By symmetry this holds for every ε. `sliding_iss` enters one kernel, by its argmin tie-break.
- **Near tie.** A bump at 0.051 gives two kernels at ε ≤ 10⁻⁶; only at 10⁻⁸ does one kernel form it.
- **Three bumps with one at a midpoint:** four kernels at 10⁻⁶.

No failure came from coarse grids (spacing 0.3 and 0.5), close mixed-sign bumps, bumps at ±0.96/0.97, or 1.0/0.1 amplitudes. The
non-uniformity in the margins is generic, though: on close off-grid bumps the loop's second entry is at τ = 16.2, 29.6, 33.3 (ε = 10⁻⁶,
10⁻⁸, 10⁻¹⁰) against the path's 33.5, because the next kernel creeps from 0.10 to 0.28 at a margin ≈ 0.04.

**(iv)** The tan formula and the 1/ρ factor are right. Along the trajectory the loop matches the limit to ≤ 1.4·10⁻⁸ with ρ = 1, and with
b ∈ [¼, 4] to ≤ 6·10⁻⁹ before the first entry, where wrong powers of ρ are off by 0.03–0.3. The repository never tests this (D3aj-10).
Condensation of the positions is real. "Several kernels share each feature" is not a limit property:
- On the C target the WFR loop forms 15, 12, 8, 4, 3, 3, 3 kernels at ε = 10⁻⁴ … 10⁻¹⁰, with the same entrants (5, 16, 11) as under the
  paradigm's metric.
- On the 16 random targets it forms exactly three kernels at the bumps on 10/16 at ε = 10⁻⁸ and on 16/16 at 10⁻¹⁰.

**Novelty.**
- (i) uses Berthier–Pillaud-Vivien's route, applied to an ε-dependent family of potentials.
- (ii) is their entropic case.
- The headline "metric dichotomy" (one kernel per feature against condensation) is refuted in the limit; what survives is only where the
  massless kernels end up, plus a slower convergence in ε under WFR.
- Missing references: Kunin et al. 2024, "Get rich quick" (conserved quantities and layer-specific learning rates set lazy/rich,
  unbalanced inits freeze or accelerate features); Labarrière et al. 2024 (deep diagonal nets as mirror flows); Maennel–Bousquet–Gelly
  2018 (feature quantization at small init).

**Code and tests.** 27 of 53 mutants survive, and about 20 of those are not equivalent.
- `iss_path`: NNLS replaced by unconstrained LS, the sign taken from Φᵀy, the exit rule disabled or flipped, the hitting time computed via
  |u|, and the exit tolerance raised to 10⁻² all survive. They change the experiment's design 4/5 event lists, and the final w misses
  weighted basis pursuit by up to 0.37. The test's design has no exits.
- `moving_rhs`: dropping the metric's factor m (centres 12× slower), dropping a′, and scaling the WFR speed ×1.2 all survive.
- `wfr_inactive_rhs`: k×2, k×0.8, tan → linear, ρ → 1/ρ, and dropping 1/ρ all survive.
- `sliding_iss`: changing the exit tolerance to 10⁻² or 0, the seed to 10⁻³, using correlations at c0, and the |u| hitting time all
  survive.
- `ParticleSet`: the free-particle drive without a′ survives all 23 tests.
- The experiment has no tests at all.
- Killed: `transit_scale`, `critical_w`, `critical_rhs`'s curve, the centre's sign, the seed's sign, seed 10⁻¹¹, t_fast 10⁴ or 10³, and
  `core_tail` (the last only by D‴.45's test).

**Numerics paragraph.** Every number is transcribed correctly. These are overstated or wrong:
- "slopes … 0.471" (an integration artefact, D3aj-8);
- "(≈ ε)": the displacements are ≈ 2.3–3.3·10³ε, i.e. ∝ ε;
- "17–50": the actual weights are 18.8–40.5;
- "14 of 16": it is 15 of 16 (D3aj-9);
- the C5 "56 events": 25 of them are immediate re-entries of one kernel at τ = 23.95;
- C6's "predicted speed 0.039 moves it by 0.16": the speed formula is 3.9× too large, and `disp` is measured from c0, not over the window.
  Kernel 3 was already at −0.7038 at τ = 15.5, and its in-window motion is ≈ 2·10⁻⁵;
- "two bumps 0.2 apart are resolved": that target sits on grid kernels (±0.1), so no sliding is tested.

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3aj-1 | major | (iv) "several kernels share each feature" and the metric dichotomy are false as ε → 0. WFR forms 3 kernels for ε ≤ 10⁻⁸ on the C target (15, 12, 8, 4, 3 for ε = 10⁻⁴ … 10⁻⁸) and on 16/16 random targets at 10⁻¹⁰. | D48_draft (iv) and New; pieces ROW0, BULLET, README, AUDIT: restate as condensation of the massless kernels' positions plus slower convergence in ε. Add an ε-sweep to 10⁻¹⁰ to C2 and C4. |
| D3aj-2 | major | (iii)(b) fails at ties: one bump midway between grid kernels is formed by 2 kernels with 0.40 each at every ε. Near ties fail non-uniformly. `sliding_iss`'s argmin tie-break is not the limit. | Add a genericity hypothesis to (iii)(b) and state the equal split at ties. Implement simultaneous entries in `sliding_iss`. Add the tie targets to C5. |
| D3aj-3 | major | The creep speed (πa_max²/(4c))εμ^{−3}⟨r∂_cφ⟩ is too large by π³/8 and lacks b_i^{3/2} and sign(u_i); the true speed uses cos/sin³(πμ/2) (checked to 3–4 digits). C6's "validation" compares the total displacement since τ = 0 with speed × window. | (iii)(a) and its proof; experiment C6: report the instantaneous dc/dτ against both formulas and cc[kb]−cc[k2]; numerics (6). |
| D3aj-4 | major | `sliding_iss` is neither faithful nor robust. Exits fire without an inward-velocity check, t_fast is fixed with no convergence test, and the seed is coupled to t_fast. Results: t_fast = 10⁴ recovers 0/17 targets; seed 10⁻⁹ makes target 3 fail (1 kernel, 14 exits); seed 10⁻¹¹ fails the paper's own test target. The defaults produce spurious exit/re-entry cycles (target 3; wider bumps; close bumps; spurious late entries at τ ≈ 5·10⁴). On edge bumps it does not return within 14 min. | `rich_transit.sliding_iss`: exit only if sign(u_j)⟨rφ_j⟩ < 0 at the relaxed state; stop the fast flow on a gradient tolerance or double t_fast until converged; scale t_fast ∝ seed^{−1/3}/\|v_j\|; check seed → 0; handle ties. Regenerate C1 and C4–C6. |
| D3aj-5 | major | The tests do not discriminate the claims: 27/53 mutants survive, including `iss_path`'s exit rule, NNLS and hitting time, `moving_rhs`'s m and a′, the WFR limit's k, tan and 1/ρ, and `sliding_iss`'s exit and seed parameters. The test design has no exits and its one-sided path equals the two-sided one. The frozen-bound test uses a bound that is 4× loose. The experiment is untested. | tests/test_mirror_free.py: add designs 4/5 and a binding sign constraint, and compare event lists with an independent ISS. tests/test_rich_transit.py: test dc/dτ against Gβa⟨r∂_cφ⟩ and the creep formula, the WFR limit at τ = 3, 5, 7 with ρ ≠ 1, a tie, and an exit. |
| D3aj-6 | minor | (o) and New: "depth set by the start's amplitude" is a continuum, not a dichotomy. Branch starts with α(0) ≲ √ε are on the power scale; α(0) ∝ ε^p (p < ½) gives ε^{−p}log(1/ε). | Restate (o) as two limits with the crossover at α(0) ~ √ε (D48_draft (o), New, AUDIT). |
| D3aj-7 | minor | The near-wall creep is generic on close bumps. Second entries are at 16.2/29.6 vs 33.5, 18.0/23.7 vs 24.7, and 34.7/54.6 vs 58.5, so the limit shows only near ε ≈ 10⁻¹⁰. "Well-separated" is not quantified. C5's 0.2-apart target is on-grid. | (iii)(b): state ε ≪ μ_k³ for the margins left at each entry. C5: use off-grid close bumps. |
| D3aj-8 | minor | The (α, log β) integration drifts off Q. The relative drift of the start mass reaches 2.3% at ε = 10⁻⁸, 4 at 10⁻¹⁰ and 259 at 10⁻¹², so exited atoms land on wrong curves. This produces design 5's slope 0.471 (true 0.4998) and loop errors of 0.05–0.35 at 10⁻¹⁰–10⁻¹². | experiments/exp_mirror_transit.py (A, and C at ε ≤ 10⁻⁸): integrate α only, with β = εb + 2c sinh²(α/a_max). |
| D3aj-9 | minor | "Entering the path's kernels on 14 of 16" is miscounted. Target 3's sliding order [2, 3, 3, 12] contains a spurious re-entry; with t_fast = 10⁶ it is [2, 3, 12], matching the loop, so the count is 15/16. | part_C / report: compare after removing exit/re-entry pairs; numerics (4). |
| D3aj-10 | minor | The WFR limit is checked only at τ = 10, at the attractor, where wrong limits also agree to 3·10⁻⁹. ρ = 1 throughout. The formula itself is right (my trajectory and ρ ≠ 1 check). | C2 and test: max over τ < τ_first of \|c−c_lim\| with unequal masses. |
| D3aj-11 | minor | Proof slips: the energy identity should read N∫\|ζ̇\|²; the RW 7.17 step for the composition with Φᵀ is missing; the "−O(1)" is superfluous; the inactive margins needed for w^ε → w are not stated. | D48_draft, proof of (i). |
| D3aj-12 | minor | Novelty is overstated: (ii) is Berthier–Pillaud-Vivien's entropic case, (i) uses their route, and the metric dichotomy is refuted. Kunin et al. 2024, Labarrière et al. 2024 and Maennel et al. 2018 are not cited. | New paragraph; lit_S75; AUDIT. |
| D3aj-13 | minor | (ii) is practically invisible: errors are 40–80% of the weights at 10⁻⁴, entry errors exceed the first entry, and the ± pairs' first ratio is 2.65. | Qualify (ii) as asymptotic; numerics (B). |
| D3aj-14 | minor | `ParticleSet`'s free branch is untested (the a′ mutant survives everything). The η's α-term is invisible to the D48 tests. The unit test's one-sided case is degenerate. | tests/test_mirror_free.py. |
| D3aj-15 | cosmetic | "(≈ ε)" should be ∝ ε; "17–50" should be 18.8–40.5; "(iii)(a)'s bound" is called a speed. The (iii)(b) formula uses c_j(0) where the code uses the slid c (and the c0 mutant survives). The margin is relative to Λρ in the statement but to λ₁ in the code. | numerics_D48; D48_draft (iii); `sliding_iss` docstring. |

## Lesson
The soft-analysis limit (Mosco, then Brezis–Pazy) was the part that held up. What failed were the claims resting on finite-ε pictures and
on a discretisation tuned to the showcase targets. Four cheap moves exposed them:
- sweeping ε two to four decades past the experiment, which turned "several kernels per feature" into one;
- measuring instantaneous rates instead of accumulated displacements, which exposed a factor of π³/8;
- feeding symmetric or degenerate inputs (ties), which broke "one kernel per feature";
- checking that a test can tell the claimed formula from a wrong one. A check made at an attractor, a design with no exits, or a single
  favourable parameter setting validates nothing about the trajectory or the algorithm.
