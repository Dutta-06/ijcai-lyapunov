# Referee report D3j (round 17) — PARADIGM §5d: Theorem D‴.23 (statement (i)–(vi), proof, Status, *Numerics*); the §0 and §7 rows "| D‴.23 |", the §8 bullets "Done (Theorem D‴.23(ii))" (branch 7) and "Done (Theorem D‴.23)" (branch 8, with its *Open*); `results/RESULTS_PARADIGM.md` (Theorem D‴.23); `src/core/adaptive_prior.py`; `tests/test_adaptive_prior.py`; `experiments/exp_adaptive_prior.py` → `results/chaos/adaptive_prior.md` (+ `adaptive_prior_raw.npz`); LITERATURE [S49] and the D‴.23 audit bullet; README's D‴.23 sentence

Branch `claude/festive-bardeen-eopej2` at 1aeaa56 (HEAD when the review started). During the review the concurrent session committed
b76b3b9 (rounds 15–16 applied to D‴.20 and D‴.21) and left uncommitted edits in the working tree: `src/core/adaptive_prior.py`, the
manuscript, and new untracked files `exp_adaptive_floor.py`, `floor_bias.py`, `test_adaptive_floor.py`. The commit does not touch the reviewed
material; of the uncommitted edits, only the one to `adaptive_prior.py` does, and the review ignores it. The D‴.23 block (statement through *Numerics*), both rows, both §8 bullets, the RESULTS section, [S49], the audit bullet and the five D‴.23
files are identical at 1aeaa56 and b76b3b9 (text diff of the extracted blocks, `git diff --quiet` per file). The README sentence is unchanged;
only "rounds 7–14" became "7–16". The uncommitted `adaptive_prior.py` was not used: every script imports the 1aeaa56 tree, extracted with
`git archive` into the scratch directory.

Line numbers are 1aeaa56's. At b76b3b9:
- add 61 in §5d (Theorem D‴.23 starts at l. 2626) and in §7 (the row is at l. 5166);
- add 62 in §8 (l. 5362 and l. 5397);
- the §0 row stays at l. 117;
- RESULTS l. 503 → 516; LITERATURE [S49] l. 70 → 71 and the audit bullet l. 340 → 350.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread; pytest ran with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the scratch tree rewrites `adaptive_prior.md`
byte-identically (md5 `2578bd0f…`).

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3j/`:

| script (in `scripts/` unless noted) | what it does |
|---|---|
| `recompute_md.py` (top level) | every number of the md, the PARADIGM paragraph, the rows and RESULTS, recomputed from the raw records independently of `report()`: λ_j^− by `np.roots`, pair solutions by `expm`, all nine modes |
| `check_partD_raw.py` | the settled gate's state before opening; the mode content of ȳ at s = 3000; 𝒥 between records in all 112 m-runs |
| `check_partD_rerun.py` | bitwise re-runs of stored systems; a test of (v) with separations at the start |
| `check_open_state.py` | the state at the opening refresh itself; K₀-atom runs from it; the semi-relaxed OT identity |
| `check_partC_floor.py` | Part C at dt = 0.25, 0.125, 0.0625; the plan along the runs; slaved fast amplitudes from the nonlinear drift |
| `check_partE.py` | the plateau's Hessian at five finite-difference steps, with Richardson extrapolation; a nearby critical point; the unregularised gradient flow to s = 2·10⁶ |
| `check_linear.py` | the pair algebra; the joint Jacobian at m = 12, 30 in separation, class-mean and prior coordinates; the first-order bias for four shifts |
| `check_unbounded.py` | boundedness under fixed and adaptive priors; the relabelling converse |
| `check_stationary.py`, `check_stationary2.py` | a non-replicated stationary point of the joint flow |
| `check_iv_gate.py` | (iv)'s drift by finite differences; 𝒥 at every plan refresh across the gates |
| `check_drag_ref.py` | the K₀-atom run from the exact representation |
| `check_v_radius.py` | (v)'s radius bound and the starts' norms |
| `check_tests.py` | the tests' asserted quantities and tolerances |
| `check_misc.py` | the leak as an internal force; transit against settled at s = 3000 |

Also in that directory: `tests.out` (7 passed, 2.09 s), `from_raw.out` and `misc_repo.out`. The 30-minute experiment was not rerun. Single
systems were reproduced bitwise, see confirmation 3.

## Verdicts

| object | verdict |
|---|---|
| **Theorem D‴.23 (i)** | **holds for replicated configurations**. The gradient structure is right: with the metric (1/m)G⁻¹ ⊕ (σ₀q_j/η)G⁻¹ the joint drift is −M⁻¹∇𝒥 (finite-difference error 2.6·10⁻¹¹; M·Jacobian symmetric to 2.6·10⁻¹¹). Proposition D‴.17's vertex-selection/KL argument carries over, because 𝒥 is a minimum over the vertices of the fixed transportation polytope of exp-polynomial functions of (Θ, Z). Three problems remain. The heading and the rows say "stationary points" where only replicated ones are characterised, and non-replicated stationary points with a non-zero leak force exist (D3j-1). The converse holds only up to relabelling of equal-mass atoms (D3j-1). Boundedness is a genuine hypothesis: Step 0's a-priori bound is lost, and solutions can diverge (D3j-2). |
| **(ii)** | **holds.** Re-derived by hand, and checked on the joint Jacobian at m = 12 and 30. The separation block is −σ₀I to 10⁻¹². Separations, class means and prior offsets do not couple (≤ 2.2·10⁻¹¹). The prior rows are −ηI on ζ and independent of the separations. The collective spectrum equals the pair rates to eight digits. Self-adjoint in the metric of (i). The lower bound has no violation on a grid of (h, η). λ⁻ increases in h (dλ⁻/dh = (η − λ⁻)/(λ⁺ − λ⁻) > 0), but this is not stated. The displayed strict "<" fails exactly when h_jη = 0, the case the next sentence treats (D3j-6). |
| **(iii)** | **holds to first order.** The relative error halves with the shift (4.80, 2.34, 1.16, 0.58 % at 0.15, 0.075, 0.0375, 0.01875), and L̄ = ½yᵀH_Ly + O(‖y‖³). At the experiment's shift the norm agreement hides 15–27 % errors in the weak modes (D3j-8). "The fit is restored faster than the parameters" is a correct qualitative reading. |
| **(iv)** | **holds as stated (formal).** The η = ∞ drift equals −mG∇(L̄ + (σ₀/2)Σγ‖ϑ − T̂(γ, Θ)‖²) with the plan held to 5.9·10⁻¹¹ (envelope theorem). "Capacity-constrained K₀-means distortion" is apt: size-constrained K-means with cluster sizes mq_j at commensurate m, an LP relaxation otherwise. "Formal" is adequate; only the η → ∞ limit through plan switches is unproved. Nothing is claimed as proved beyond the drift identity. |
| **(v)** | **holds as a sketch; "proved" is overstated** (D3j-3). The ingredients are right and were verified. Missing: the radius bookkeeping, in particular the plan margin for the prior's offsets, which the |||·|||_η-ball bounds only by ρ(η/(σ₀q_k))^{1/2}; the degenerate limits η → 0 and η → ∞; and any statement that the runs start 94–6·10⁴ times outside the proved ball. The m-particle numerics offered for (v) do not test it (D3j-4). A direct test with separations at the start supports it (confirmation 9). |
| **(vi)** | **holds in substance.** The drag is the transit's: the K₀-atom run from the exact representation has ‖ȳ‖(500) = 4.4·10⁻³ at η = σ₀ and 1.6·10⁻¹⁷ at η = ∞, against the m-systems' 0.032–0.71. The residue sits in the weak modes (projection). The η = ∞ plateau is a plateau: the negative eigenvalue is resolved, the gradient flow leaves it into the exact representation's basin, and no critical point lies near it. The remedy works but costs a delay (D3j-11). "Gating keeps 𝒥 non-increasing" is true of every η-schedule, gated or not. At η = ∞ the opening is a downward jump, not a flow (D3j-9). |
| Proof | **correct where it is a proof.** (v) is a pointer (D3j-3). "Conversely ż = 0 forces z_k = T̂_k = ϑ_k" needs "up to the plan's labelling" (D3j-1). |
| *Status* | the (i)–(iv) and (vi) labels are accurate; "(v) proved" is overstated (D3j-3). The novelty claim, "new are … the exact linear theory … with its mode-by-mode trade-off", is overstated: the per-mode formulas are the classical quadratic analyses (D3j-5). The four [S49] citations are bibliographically correct. The description of the tests is accurate, apart from test 7's name (D3j-10). |
| *Numerics* paragraph | **every number reproduces** (confirmation 2). Three readings fail or overreach: "the uniformity in m of (v)" (D3j-4), "every mode amplitude above 10⁻⁶" (only five of nine modes computed; D3j-7), and "describes it to 4.8 %" (D3j-8). The floor explanation, the drag, the residue in the weak modes and the plateau all hold. |
| `src/core/adaptive_prior.py` | **holds.** Drift, energy, gates, H_L, pair rates and pair solutions are correct. The η = ∞ refresh lowers 𝒥 (−44 % and −96 % at the openings). The docstring repeats the overgeneralisation of D3j-1. |
| `experiments/exp_adaptive_prior.py` | **holds.** Regeneration is byte-identical. The ungated raw entries, which lack the `gate`/`t_open` keys and so came from an earlier `m_job`, are reproduced bitwise by the committed code. `report()` computes slopes and ratios for modes 1–4 and 9 only (D3j-7). Part E rests on BFGS alone (D3j-12). |
| `results/chaos/adaptive_prior.md` | **numbers hold**; wording as in D3j-7, D3j-12. |
| `tests/test_adaptive_prior.py` | **pass (7/7, 2.09 s); none vacuous**, and the tolerances discriminate: test 6 rejects a doubled prior weight (asymmetry 7.8·10⁻² against 2.6·10⁻¹¹). Test 7's name is not what its tolerance 10⁻³ checks (D3j-10). |
| §0 row, §7 row, §8 bullets, RESULTS, LITERATURE, README | **numbers hold** (115 tests at 1aeaa56 is right). Wording: D3j-1 (§0, §8 branch 8), D3j-2 (§8 branch 7), D3j-3 (§0), D3j-4 (§0 numerics clause), D3j-5 (§0 last column, audit bullet), D3j-6 (§0, §8 branch 7, audit), D3j-7 (§7, RESULTS), D3j-11 (§0 "must", README). The §8 *Open* claim on learned masses **holds**: min_q W²(μ, ν_{Z,q}) is the nearest-atom distortion (semi-relaxed OT, to 1.8·10⁻¹⁵), and with free masses the plan becomes the nearest-atom assignment, which conserves the cluster masses (Proposition D‴.12(i)), so the floor returns. |

## Confirmed issues

**D3j-1 (minor; (i) heading l. 2573 "a gradient flow whose stationary points carry no leak bias", l. 2576–2580, proof l. 2628–2629; §0
row l. 117 "its stationary points carry no leak bias (the critical points of the unregularised K₀-atom loss)"; §8 branch-8 bullet
l. 5336 "its stationary points are the unregularised K₀-atom critical points (no leak bias)"; module docstring of `adaptive_prior.py`
"the stationary points are the critical points of the unregularised K₀-atom loss").** The body characterises *replicated* stationary
configurations only, and the heading and summaries drop that restriction. Three defects follow.

(a) *Non-replicated stationary points carry a leak force.* At a stationary point every particle satisfies G⟨r∇Φ̃(ϑ_i)⟩ =
σ₀(ϑ_i − T̂_{π(i)}) with z = T̂. Summing over a class gives only that the class-average learning force vanishes; each particle may still be
held by the leak against the loss. Example (`check_stationary.py`): K₀ = 3, q = 1/3, m = 6, on a four-bump target that a six-particle
network fits exactly with one class split (fit error 1.8·10⁻¹⁶).
- The joint flow has a stationary point with that class split at G⁻¹-distance 1.64, against 3.29 at the unregularised critical point. Its
  drift is ≤ 3.4·10⁻¹⁷ in Θ and 0 in Z, for η = σ₀ and 10σ₀, and its optimal plan is the class partition.
- The leak force on each split particle is 0.0058 (whitened), and the unregularised loss gradient there is 5.8·10⁻³: this is not a
  critical point of the unregularised loss. Fit error 4.2·10⁻².
- It is a saddle of 𝒥: the largest real part of the joint Jacobian is +3.7·10⁻⁵.
- From the exact split fit, the joint flow collapses the class by s ≈ 8000 and converges to a replicated K₀-atom critical point with fit
  error 1.9·10⁻² (`check_stationary2.py`). The exact fit is not stationary: this is the K₀-cluster bias of (iv)'s regulariser.

(b) *The converse holds only up to relabelling.* ż = 0 forces z_{σ(k)} = ϑ_k for the plan's class-to-atom bijection σ, not z_k = ϑ_k. At
the exact representation with the prior cyclically permuted (z_k differs from ϑ_k by up to 6.2 in the unconstrained coordinates), the plan
maps the classes to (2, 0, 1), W² = 0, and the joint drift is ≤ 1.1·10⁻¹⁷ (`check_unbounded.py`).

(c) *"The stationary set does not depend on the initial prior" is a tautology.* Z is a state variable, so the initial prior never enters
the stationary equations. Which stationary point is *reached* does depend on it: Part E's drag leaves the system on a plateau for more
than 2·10⁴ slow units.

*Fix:*
- Heading: "(… whose replicated stationary points carry no leak bias)". "Iff z_{σ(k)} = ϑ_k for the plan's bijection σ (equal masses)".
- Add: "non-replicated stationary points — a class held together by the leak toward its barycentre — carry the K₀-means bias of (iv)",
  with the example.
- Replace the tautology by "the stationary equations do not involve the prior's initial position (the limit reached may: (vi))".
- §0, §8 and the docstring accordingly.

**D3j-2 (minor; (i) l. 2575–2576 "bounded solutions converge to Clarke-critical points (Proposition D‴.17's argument …)", proof l. 2626;
§8 branch-7 bullet l. 5300–5304 "the measure-space form of the classical trade-off between a leak's robustness and persistent
excitation").** Proposition D‴.17(ii) bounded every solution by Step 0's a-priori estimate, using conv supp ν of the *fixed* prior. With an
adaptive prior the leak is an internal force between a class and its atom. With the plan held,

    d/ds[(1/m)Σ_{i∈k}ϑ_i + (σ₀q_k/η)z_k] = (1/m)Σ_{i∈k}G⟨r∇Φ̃(ϑ_i)⟩,

to 10⁻¹⁸ (`check_misc.py`). Nothing confines the pair, the sublevel sets of 𝒥 are unbounded (W² = 0 with Z following Θ, L̄ bounded), and
solutions can diverge. Example (`check_unbounded.py`): the instance with the first bump's weighted amplitude 3A₀ = 4.8 > a_max = 4.

| prior | ϑ_{a,0} at s = 10³, 10⁴, 10⁵, 10⁶ | outcome |
|---|---|---|
| fixed (η = 0) | 5.869 from s = 10³ on | converges, drift 1.7·10⁻¹⁷ at s = 10⁶, fit error 0.175 |
| adaptive, η = σ₀ | 7.99, 11.79, 16.25, 20.83 (≈ 2 ln s) | diverges; z_{a,0} follows to 3·10⁻⁴, a₀ → a_max, drift ∝ 1/s, infinite length |
| adaptive, η = ∞ | 22.21 at s = 10⁶ | diverges |

On the instance itself, a least-squares search for ∇L = 0 from the plateau atoms runs off to a critical point at the saturation boundary
(a₃ = a_max, L = 2.5726·10⁻⁴). The theorem is conditional, but it never says that boundedness is now an assumption; the confinement it
loses is precisely the "robustness" half of the classical trade-off. The §8 bullet invokes that trade-off through the rates only.

*Fix:* in (i), "bounded solutions (boundedness is no longer automatic: the leak toward an adaptive prior does not confine the parameters;
example …)". In §8 branch 7, add that the adaptive prior gives up the leak's a-priori bound, which is the robustness side of the trade-off.

**D3j-3 (minor; (v) l. 2607–2614, proof l. 2637–2638, *Status* l. 2641 "(v) proved, as Theorem D‴.18's contraction argument …", §0 row
"(i)–(iii) and (v) proved").** The proof of (v) is two sentences pointing to Theorem D‴.18. The structural part is correct and verified
(`check_linear.py`, test 6 of `check_tests.py`): separations decay at −σ₀I, decoupled from (ū, ζ), and the collective block is
self-adjoint in the metric with eigenvalues −λ_j^±. The radius bookkeeping is missing, in three respects.

(a) *The plan margin for the prior.* A |||·|||_η-ball of radius ρ bounds the prior's offsets only by ‖ζ_k‖_{G⁻¹} ≤ ρ(η/(σ₀q_k))^{1/2}. The
diagonal plan stays optimal only if ρ ≲ (σ₀q_min/η)^{1/2}·(margin)/(C(1 + q_min^{−1/2})), so for finite η the radius shrinks like
η^{−1/2}. At η = ∞ the prior is slaved and the norm reduces to D‴.18's, which needs its own sentence.

(b) *η → 0 and the constants.* As η → 0, λ_1^− ≤ η → 0, so the rate and the radius κλ_1^−/C_N vanish, and the ζ-extent of the ball
vanishes faster still. C_N must be the |||·|||_η-Lipschitz constant of the Jacobian (round 14's D3g-2 fix). Since the prior's part is
linear, it equals the particles', and a sampled lower bound from the collective block alone is 0.34 (`check_v_radius.py`).

(c) *The runs are far outside the proved ball.* With ρ ≤ min{σ₀, λ_1^−}/(4C):

| η | bound on ρ | |||·|||_η of the starts (Θ*(q), shifted prior) | factor outside |
|---|---|---|---|
| σ₀/10 | 2.0·10⁻⁵ | 1.17 | 6·10⁴ |
| σ₀ | 1.1·10⁻⁴ | 0.37 | 3344 |
| 10σ₀ | 2.0·10⁻⁴ | 0.12 | 584 |
| ∞ | 2.2·10⁻⁴ | 0.021 | 94 |

The radius is λ_1^−/σ₀ = 0.0038–0.043 times D‴.18's at the same C_N. Neither the md nor the paragraph says so; round 14 required exactly
this for D‴.18 (D3g-4).

*Fix:*
- "(v) sketched as an adaptation of Theorem D‴.18's contraction argument".
- State ρ(η) = min{κ·min{σ₀, λ_1^−(η)}/C_N, (σ₀q_min/η)^{1/2}·margin/(C(1 + q_min^{−1/2}))}, with C_N the Jacobian-Lipschitz constant,
  and treat η = ∞ in the ζ-free norm.
- Quote the proof's radius on the instance and present the runs as probing beyond it.

**D3j-4 (minor; *Numerics* l. 2672–2676 "at s = 3000 ‖ȳ‖ and the fit error of every system equal the K₀-atom run's at s − t_open —
per-seed ratios 0.995–1.000 (η = σ₀) and 0.968–1.001 (η = ∞) for every m from 15 to 510, the uniformity in m of (v)"; the md's last
Part D paragraph; §0 row "started from the fixed-prior stationary point every size m = 15–510 follows the K₀-atom run …"; §7 row;
RESULTS l. 515–517).** The comparison does not test (v). By the opening, the clouds have been leaking at σ₀ for 1211–1304 slow units.

- *The systems open replicated.* Over all 32 settled runs, at the last record before the opening, the rms separations are
  1.2·10⁻⁴–7.9·10⁻⁴ and the class means lie 7.4·10⁻⁷–6.7·10⁻⁵ from Θ*(q) (`check_partD_raw.py`). Captured at the opening refresh itself: seed 3, m = 15 has rms 7.4·10⁻⁴ (max
  1.5·10⁻³) and an offset of 6.1·10⁻⁵; seed 1, m = 510 has 1.8·10⁻⁴ and 4.2·10⁻⁶ (`check_open_state.py`).
- *Each m-system is then the K₀-atom system started from its own class means.* That run reproduces the m-system at s = 3000 to 10⁻⁵ in
  ‖ȳ‖ and in fit error (ratios 1.00001/1.00000 and 1.00000/1.00000), by the invariance of the replicated configurations — for any m,
  independently of (v).
- *The departures from 1 are the opening offset, not a finite-m effect.* For seed 3, m = 15 the offset in the weakest mode is 5.7·10⁻⁵,
  and the K₀-atom run from Θ*(q) gives 0.968/0.979 (η = ∞) and 0.995/0.998 (η = σ₀) — exactly the md's per-seed values. The departures
  are larger at small m because the transit is.

A test of (v) needs separations at the start (`check_partD_rerun.py`, Part 3). Start from the lift of Θ*(q) plus zero-mean per-particle
separations, with the prior at the shifted centres, and compare ‖ȳ‖ at s = 1750 with the K₀-atom run:

| η | rms separation 0.02, m = 15 / 510 | rms separation 0.1, m = 15 / 510 |
|---|---|---|
| σ₀ | +0.15 % / +0.11 % | +3.7 % / +2.7 % |
| ∞ | +1.0 % / +0.7 % | +24 % / +17.5 % |

The effect is second order in the separations (factor 25 for a factor 5) and nearly m-independent, which is what (v) and facts (1)–(2)
predict.

*Fix:* replace "the uniformity in m of (v)" with: "as the invariance of the replicated configurations requires: the systems open
replicated to 10⁻³, and the K₀-atom run started from their class means at the opening reproduces every system to 10⁻⁵; the departures
from 1 are the opening offset in the weakest mode". Add the separated-start runs as the test of (v), and have `report()` print the
opening state.

**D3j-5 (minor; *Status* l. 2645–2648 "New are … the exact linear theory at an exact representation with its mode-by-mode trade-off
between the leak's bias σ₀/(h_j + σ₀) and the rate λ_j^− …"; §0 row last column "new are … the mode-by-mode trade-off between the leak's
bias and the rate set by the excitation spectrum"; audit bullet l. 340–350 "Verdict: … new mechanism (the bias–excitation trade-off per
mode, uniform in m)").** At an exact representation the collective dynamics is the linearisation of a quadratic with Hessian H_L. Every
per-mode formula is classical:

- σ₀/(h_j + σ₀) is the Tikhonov (ridge) residual factor in the Gauss–Newton eigenbasis.
- Iterated Tikhonov, i.e. the proximal point method on a quadratic (x_l = σ(AᵀA + σI)⁻¹x_{l−1} + (AᵀA + σI)⁻¹Aᵀf, Hanke–Groetsch, JOTA
  98 (1998) 37–53), contracts mode j by σ₀/(h_j + σ₀) per step. The slow-prior limit of the lower bound, ηh_j/(h_j + σ₀), is its
  continuous-time form: the gradient flow of the Moreau envelope, whose curvature is h_jσ₀/(h_j + σ₀).
- EASGD's own analysis is the stability and convergence analysis of the 2 × 2 worker–centre linear map on a quadratic, in terms of the
  learning rate and the moving rate. This is the discrete-time version of the pair matrix [[−(h + σ₀), σ₀], [η, −η]].
- Empirical Bayes / type-II ML for the mean of a Gaussian prior is the same proximal iteration.

What is specific to D‴.23:
- the identification of H_L as the measure-space excitation spectrum;
- the decoupling of the separations at σ₀ (hence the m-uniformity, by D‴.18's machinery);
- the mass constraint;
- the transit's drag and the gate.

The "classical trade-off" sentence covers only the rate half (see D3j-2).

*Fix:* cite the quadratic analyses (Tikhonov/iterated-Tikhonov filter factors, Hanke–Groetsch; EASGD's quadratic stability analysis) as
the source of the per-mode formulas, and restrict "new" to the transfer to the mass-constrained particle system, uniformly in m, and to the
drag and its gate. The literature check is limited: arxiv.org and proceedings.nips.cc were unreachable from this session (403/DNS), so the
EASGD and Hanke–Groetsch statements rest on search-result summaries and memory.

**D3j-6 (cosmetic; (ii) l. 2589–2592, §0 row, §8 branch-7 bullet l. 5302, audit bullet; the hypotheses "η ≥ 0" l. 2566 and (i)).**
"ηh_j/(h_j + σ₀ + η) ≤ λ_j^− < min{h_j, η}" fails as displayed exactly when h_jη = 0, where λ_j^− = 0 = min{h_j, η}. The sentence that
follows treats h_j = 0. For h_j, η > 0 both inequalities are strict; the lower bound is an equality only when h_jη = 0, and is tight as
η → 0 (ratio 0.9895 for h₁ at η = σ₀). At η = ∞, λ_j^− = h_j. The preamble allows η ≥ 0, but (i)'s metric (σ₀q_j/η)G⁻¹ and its "iff"
need η > 0: at η = 0, ż ≡ 0. The monotonicity of λ^− in h, used for "the slowest rate is λ_1^−", is not stated.

*Fix:* "for h_j, η > 0: ηh_j/(h_j + σ₀ + η) < λ_j^− < min{h_j, η}; λ_j^− = 0 if h_jη = 0 and λ_j^− = h_j at η = ∞; λ^− increases with h
(dλ^−/dh = (η − λ^−)/(λ^+ − λ^−))"; (i) for η > 0.

**D3j-7 (cosmetic; *Numerics* l. 2660–2662 "every mode amplitude above 10⁻⁶ decays at its predicted rate (log-slopes on [3000, 6000]
between 0.981 and 1.002 times λ_j^−)"; md Part C; §7 row l. 5105; RESULTS l. 510).** `report()` computes slopes and ratios for modes 1–4
and 9 only. At η = σ₀/10, modes 5–8 are above the floor on [3000, 6000], with log-slopes 1.0031, 1.0000, 1.0004 and 1.0026. Over every
mode above 10⁻⁶ the range is 0.981–1.003; the ratio range 0.820–1.115 is unchanged (`recompute_md.py`). Separately, the §7 row's "the
transit gate leaves a residue … (‖ȳ‖ 4.9·10⁻³–3.7·10⁻² against 1.3·10⁻³)" omits "(η = σ₀)".

*Fix:* print all nine modes, or say "modes 1–4 and 9", and quote 0.981–1.003; add "(η = σ₀)" to the §7 row.

**D3j-8 (cosmetic; *Numerics* l. 2658–2660 "the first-order formula y* = σ₀(H_L + σ₀I)⁻¹y⁰ describes it to 4.8 %"; md Part B).** The
4.8 % is a norm, dominated by the strong modes. Per mode, actual/predicted is 0.849 (mode 1), 1.052, 0.727 (mode 3), 0.917 (mode 4) and
then within 5 % (`recompute_md.py`). These are the weak modes the trade-off is about. The first-order law itself is right: halving the
shift from 0.15 to 0.01875 takes mode 3 through 0.727, 0.871, 0.937, 0.969 and mode 1 through 0.849, 0.917, 0.956, 0.977, and the norm
error through 4.80, 2.34, 1.16, 0.58 % (`check_linear.py`). The paragraph quotes mode 1's kept fraction only (0.814 against 0.959). It also
writes y⁰ for the theorem's x⁰.

*Fix:* add "(per mode: 0.85 and 0.73 of the prediction in modes 1 and 3, converging linearly in the shift)", and write x⁰.

**D3j-9 (cosmetic; (vi) l. 2619–2621 "Freezing the prior while mass is in transit … keeps 𝒥 non-increasing, both phases being gradient
flows of 𝒥"; proof (vi) l. 2638–2639; *Status* l. 2642 "(vi) the Lyapunov property under either gate proved"; §0 "the gates' Lyapunov
property").** 𝒥 does not depend on η, so every schedule η(s) ≥ 0, gated or ungated, keeps it non-increasing. It is not a property of the
gate. At η = ∞ the opening is a jump of Z to the class barycentres, not a flow. 𝒥 drops at it: by 44 % at the transit gate and 96 % at
the settled gate (seed 0, m = 15), against 1.5 % per refresh at η = σ₀. Between refreshes 𝒥 is strictly decreasing in every gated and
ungated run checked (`check_iv_gate.py`), and between records in all 112 m-runs (largest relative increase 2.8·10⁻¹⁴, rounding;
`check_partD_raw.py`).

*Fix:* "𝒥 is non-increasing under any schedule of η, and at η = ∞ the opening jump lowers it. The gate is there against the drag, not for
monotonicity."

**D3j-10 (cosmetic; `tests/test_adaptive_prior.py::test_the_settled_gate_opens_at_the_fixed_prior_stationary_point`, tolerance 10⁻³).**
In the test's setting, the gate opens at s = 441 with tolerance 10⁻³ (`check_tests.py`). At that moment the rms separation is 0.042,
twice the bias ‖y*(q)‖ = 0.021, and the class means are 1.8·10⁻³ from y*(q). Only at the experiment's 10⁻⁵ is it at the stationary
point (s = 1103: 4.1·10⁻⁴ and 2.5·10⁻⁶; at 10⁻⁴, s = 773: 4.1·10⁻³ and 6.6·10⁻⁵). What the test checks is the gate's logic: it opens
when the fixed-prior drift first falls below the tolerance, after the transit, with the prior frozen until then.

The other tests are meaningful.

| test | quantity | measured | tolerance / comparison |
|---|---|---|---|
| 1 | finite-difference error of the gradient structure | 2.6·10⁻¹¹ | 10⁻⁷; plan swap gap 6.4 |
| 3 | rate error | 6.5·10⁻¹² | 10⁻⁶ absolute (0.66 % of the smallest rate) |
| 4 | first-order error | 1.55 % | 5 % |
| 5 | opening time | t_open = 15 | — |
| 6 | asymmetry of M·Jacobian | 2.6·10⁻¹¹ | 10⁻⁷; a doubled prior weight gives 7.8·10⁻² |

*Fix:* rename test 7 ("… opens when the fixed-prior drift first falls below its tolerance, after the transit"), or run it at 10⁻⁵ with an
assertion on the opening state's distance from Θ*_m(q).

**D3j-11 (cosmetic; (vi) l. 2622–2624 "The adaptation should therefore start from the fixed-prior stationary point"; §0 "so the
adaptation must start from the fixed-prior stationary point"; README l. 51–52 "provided the adaptation starts after the transit").** The
settled gate spends ≈ 1200 of the 3000 slow units waiting. At s = 3000 and η = σ₀ (`check_misc.py`):
- it wins everywhere in ‖ȳ‖ and W²;
- in fit error it wins at m ≤ 240 (1.59·10⁻⁴ against 1.01·10⁻³ at m = 15, 1.65·10⁻⁴ against 2.13·10⁻⁴ at m = 240) and at η = ∞ for every
  m;
- at m = 510 the transit gate's fit error is lower (1.547·10⁻⁴ against 1.638·10⁻⁴).

"Must" is too strong. The README's remedy, "after the transit", is the transit gate, not the theorem's: after the transit *and* after the
fixed-prior relaxation.

*Fix:* §0 "should"; (vi) "when the parameters, not only the fit, are wanted at a given horizon"; README "after the transit and the
fixed-prior relaxation".

**D3j-12 (cosmetic; md Part E and *Numerics* l. 2669–2671 "and the loss descends from there to the exact representation").** In the md
this rests on a BFGS run, which is not the flow and could change basins. The gradient flow confirms it (`check_partE.py`). The
unregularised K₀-atom flow, which is the collapsed η = ∞ dynamics, started from the stored s = 23000 atoms (LSODA, rtol 10⁻¹⁰):
- stays on the plateau (L = 2.52·10⁻⁴ → 2.42·10⁻⁴) for 3.6·10⁴ more slow units;
- halves L by ≈ 5·10⁴;
- reaches the exact representation (L = 1.9·10⁻¹⁷) by ≈ 9.7·10⁴.

This is the same basin, on the time scale 1/1.84·10⁻⁵ = 5.4·10⁴ of the negative curvature, which bears out "crosses at that loss's own
curvature". *Fix:* quote the flow and its crossing time in place of, or next to, BFGS.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_adaptive_prior.py --from-raw` on the scratch tree rewrites `adaptive_prior.md` byte-identically (md5
   `2578bd0f…`).
2. **Every number of the md, the PARADIGM paragraph, the §0/§7 rows and RESULTS**, recomputed from the raw records independently of
   `report()` (`recompute_md.py`).
   - Part A: h_j = 3.0587·10⁻⁴ … 0.19521; σ₀/(h_j + σ₀) = 0.959 … 0.035; h_j + σ₀ = 7.3803·10⁻³ … 0.20229. λ_j^− by `np.roots` equals the
     stored rates to 2·10⁻¹⁶; λ_1^− = 2.684·10⁻⁵, 1.513·10⁻⁴, 2.780·10⁻⁴, 3.059·10⁻⁴.
   - Part B: Newton residual 4.3·10⁻¹⁷; fit error 9.7765·10⁻³; ‖y − pred‖/‖y‖ = 4.805 %; all nine mode rows.
   - Part C: all 20 rows, with ratios to the pair solutions recomputed by `expm`: 0.820–1.115; the change from s = 1000 on is ≤ 1.05 % for
     η ≥ σ₀ and up to 5.7 % at σ₀/10. Fit errors 3.612·10⁻⁴, 1.825·10⁻⁵, 8.035·10⁻⁶, 6.738·10⁻⁶. The floor amplitudes are
     5.0·10⁻¹¹–9.7·10⁻⁷ (2522 of them).
   - Part D: all 28 rows. t_open = 88–149 for the transit gate and 1211–1304 for the settled gate. Per-seed ratios 0.9953–1.0000 (η = σ₀)
     and 0.9676–1.0005 (η = ∞). Counts are q in all 112 runs.
   - Part E: all records; ‖ȳ‖ = 2.124 at s = 3000; L = 2.5220·10⁻⁴; Hessian eigenvalues −1.84·10⁻⁵, 7.16·10⁻⁵, 9.46·10⁻⁴.
   - Paragraph: "0.032–0.32 / 0.048–0.71", "0.0094–0.10 / 0.012–0.59", "4.9·10⁻³–3.7·10⁻²", "6.7–7.2·10⁻⁶", "9.3–9.9·10⁻⁷" and "3.4–3.5·10⁻⁵"
     are all exact seed-mean ranges.
3. **Reproducibility.** Single systems were re-integrated with the committed code from their seeds (`check_partC_floor.py`,
   `check_partD_rerun.py`, `check_open_state.py`); the stored records are reproduced bitwise (max difference 0.0).
   - Part C at dt = 0.25, η = σ₀ and ∞.
   - Four ungated Part D systems whose raw entries come from an earlier `m_job` (seed 0: η = σ₀ at m = 15 and 60, η = 0 at m = 15; seed 3:
     η = ∞ at m = 15).
   - The settled seed 3, m = 15 system (t_open = 1211).
4. **(i) and (iv).** Test 1's configuration has a unique plan (swap gap 6.4), and the drift is −M⁻¹∇𝒥 to 2.6·10⁻¹¹ (Θ) and 7.8·10⁻¹¹
   (Z). The η = ∞ drift is the envelope gradient of L̄ + (σ₀/2)𝒬_q to 5.9·10⁻¹¹.
5. **(ii).** The algebra:
   - no lower-bound violation on a grid of (h, η);
   - the strict upper bound fails in exactly the 12 of 42 grid cases with hη = 0;
   - λ^− increasing in h on [10⁻⁶, 10];
   - the symmetrised form with off-diagonal (σ₀η)^{1/2}.
   The joint Jacobian at m = 12 and 30 for η = σ₀ and 10σ₀ is block-diagonal as stated (confirmation of the verdict).
6. **Part C's floor is the second-order remainder.**
   - The fast-mode amplitudes at s = 3000 and 6000 are identical, to the four digits printed, for dt = 0.25, 0.125 and 0.0625, so they
     are not integrator error.
   - The plan is the diagonal at every step of every run, so they are not plan effects.
   - With modes 1–3 held at their recorded values, the fast amplitudes slaved by the full nonlinear drift reproduce the recorded floor.
     At η = σ₀, s = 6000: 1.88·10⁻⁸ against 1.99·10⁻⁸ (mode 4) and 4.90·10⁻⁹ against 5.25·10⁻⁹ (mode 5). At η = ∞, s = 6000:
     2.40·10⁻⁹ against 2.50·10⁻⁹. At η = σ₀, s = 3000, mode 4 is the worst: 1.31·10⁻⁷ against 1.65·10⁻⁷.
7. **The transit gate's residue lies in the weak modes.** At η = σ₀ and s = 3000, modes 1–3 carry all but ≤ 6.1·10⁻⁶ of ‖ȳ‖² (mode 1:
   0.85–0.99 in seed means). The amplitudes reach 50 and 43 times the K₀-atom run's in modes 1 and 3, and 0.08–3.8 times in mode 2. At
   η = ∞, mode 1 carries ≥ 0.990.
8. **The drag.** The K₀-atom run from the exact representation, the clouds' centre, with the shifted prior has ‖ȳ‖(500) = 4.4·10⁻³ at
   η = σ₀, below the md's reference 7.4·10⁻³, and stays at 1.6·10⁻¹⁷ at η = ∞. The ungated excess (0.032–0.71) is the transit's.
9. **(v) with separations** (D3j-4). The second-order, nearly m-independent effect of per-particle separations on the collective
   trajectory supports (v)'s uniformity.
10. **The plateau.**
    - The negative eigenvalue is resolved: −1.841425·10⁻⁵ for finite-difference steps 10⁻⁷–10⁻⁵, −1.841422·10⁻⁵ at 10⁻⁴ and
      −1.84106·10⁻⁵ at 10⁻³; Richardson gives −1.841425·10⁻⁵; the asymmetry of the finite-difference Jacobian is ≤ 2.7·10⁻¹⁰.
    - The whitened gradient there is 1.12·10⁻⁵ (component 9.1·10⁻⁶ on the negative direction).
    - No critical point is nearby: least squares ends at a saturated one at whitened distance 14 with a higher loss.
    - The drift passes through a minimum (1.24·10⁻⁵ at s = 12000) and then grows. "Plateau, not a minimum" is the right word, not
      "saddle". The flow reaches the exact representation (D3j-12).
11. **§8 *Open*.** min_q W²(μ, ν_{Z,q}) = (1/m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹} to 1.8·10⁻¹⁵ on 20 random instances (LP with free column sums).
12. **Tests.** 7 passed in 2.09 s; every asserted quantity recomputed (D3j-10). README's "115 tests" is right at 1aeaa56.
13. **Citations.** Rockafellar, SIAM J. Control Optim. 14(5):877–898, 1976; Zhang–Choromanska–LeCun, NeurIPS 2015; Nowlan–Hinton, Neural
    Computation 4(4):473–493, 1992; Ullrich–Meeds–Welling, ICLR 2017. All four are correct. The descriptions of the methods are fair,
    except for the novelty claim built on them (D3j-5).

## Not checked

- The 30-minute experiment was not rerun. Single systems were reproduced bitwise and the regeneration is byte-identical.
- (v)'s constants (κ, C, the plan margins) were not derived; only a sampled lower bound for the collective block was computed.
- The KL/definability argument of Proposition D‴.17 was taken as given, as were Theorem D‴.18's facts (round 14).
- A stable non-replicated stationary point on the instance itself was not searched for. D3j-1's example is a saddle on a four-bump target.
- Incommensurate sizes (split plans) under an adaptive prior, the real loop with an adaptive prior and learned masses are outside the
  material or open by the authors' account.
- The literature check is limited: the full texts of EASGD, Lookahead and Hanke–Groetsch could not be fetched (arxiv.org 403 at the
  proxy, DNS failure for proceedings.nips.cc). The statements in D3j-5 rest on search-result summaries and memory.
- The concurrent session's commit b76b3b9 and its uncommitted edits (including `src/core/adaptive_prior.py`, `floor_bias.py` and
  `exp_adaptive_floor.py`) were not reviewed.
