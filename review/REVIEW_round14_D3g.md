# Referee report D3g (round 14) — PARADIGM §5d: Theorem D‴.18 (statement, proof, Status), Remark D‴.19, the D‴.18 *Numerics*; the §0 and §7 rows "| D‴.18 |", the §8 branch-3 bullet "Done (Theorem D‴.18)", `results/RESULTS_PARADIGM.md` (Theorem D‴.18); `experiments/exp_uniform_basin.py` → `results/chaos/uniform_basin.md` (+ `uniform_basin_raw.npz`); `tests/test_uniform_basin.py`

Branch `claude/festive-bardeen-eopej2` at 5ff2e0a (the D‴.18 commit, HEAD when the review started). During the review a concurrent session
advanced the branch to 80ab38d (cffb7be entropic-leak code and tests; d6c1396 round-13 fixes to D‴.16/D‴.17; 80ab38d Proposition D‴.20,
inserted after the D‴.18 numerics). None of them touches the material: the D‴.18 block (statement through numerics), its two table rows,
the §8 bullet, the RESULTS section and the four D‴.18 files are identical at 5ff2e0a and 80ab38d (`git diff`, text diff of the extracted
blocks). Line numbers below are 5ff2e0a's; at 80ab38d add 22 in §5d (Theorem D‴.18 starts at l. 2233), and the §0 row is l. 111, the §7
row l. 4825, the §8 bullet l. 4970. No repository file was modified by this review:
every script writes only to the scratch directory, Python ran with `PYTHONDONTWRITEBYTECODE=1`, the `--from-raw` regeneration ran on a scratch
copy of the script with `src` symlinked, and the copies of `exp_uniform_basin.py` and `exp_ot_real_loop.py` used were checked equal to 5ff2e0a's (d6c1396 later changes only
report strings of `exp_ot_real_loop.py`, not `real_run`).
The regeneration rewrites `uniform_basin.md` byte-identically (md5 `4d427b58…`, 0.7 s). Scratch scripts and their `.out` files in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3g/`:
`recompute_md.py` (every entry of both tables and every range of the md and the PARADIGM paragraph from the raw records; identity of the two
leaks; counts along the run; quarter-window rates and fit residuals; the horizon dependence of C; the W² time series at ρ = 4.5),
`check_initial.py` (the initial configurations rebuilt from the seeds: particles outside their class cell, cell margins, the fraction of each
sphere in another cell), `check_rho45.py`, `check_rho45_lifts.py`, `check_rho45_crossings.py` (the ρ = 4.5 finals: within-class spreads,
distance to the lifts of the realised weights, splitting margins, drift, cell changes during the run, W² of the ideal lifts),
`check_integrator.py` (the batched step against `exp_ot_real_loop.real_run` on single systems; one system alone against inside a batch),
`check_predictions.py` (λ_sep, λ_coll and the full spectra; the 3-atom orbits by shooting and their ghost exponents at k_e = 25),
`check_m3_modes.py` (the m = 3 finals on the collective modes; the class-mean / separation split of the m ≥ 15 finals), `check_partB.py`,
`check_structure.py` (the three structural facts on the averaged loop at three replicated configurations, m ≤ 160), `check_jac_fixedpattern.py`,
`check_misc.py` (norm equivalence, fact (1), a fixed one-particle remainder, the averaged contraction at a non-zero residual),
`check_poincare_remainder.py` (the real loop's one-period map: second-order remainder in |||·||| for m = 6, 30, 150), `rel_err_collective.py`, `derived_numbers.py` and `check_C_argmax.py` (the relative errors, derived quantities and the location of C quoted below), `check_tests.py` (the
tests' asserted quantities, tolerances and stronger variants), `tests.out` (3 passed, 1.04 s), `from_raw.out`, `collect.out`. The 12-minute
experiment was not rerun; every number below comes from the raw file, from re-integrations of single systems, or from the scripts listed.

## Verdicts

| object | verdict |
|---|---|
| **Theorem D‴.18 (a)** | **holds with a radius fix.** The fixed point is unique where P contracts, and the lift is a fixed point there (its separations vanish and |||ϑ^Γ − ϑ*(w)||| = ‖Θ^Γ − Θ*(w)‖_w = O(1/k_e), the weighted K₀-atom distance, m-independent). But the contraction is proved on the |||·|||-ball of radius ρ, while the stated uniqueness set is the per-particle ball of radius ρ, which only lies in the |||·|||-ball of radius (2 − 2/n_max)ρ (sharp): uniqueness is proved in the per-particle ball of radius ρ/2 (D3g-2). |
| **(b)** | **holds**, modulo D3g-2: rate (1 − κ)min{λ_coll, λ_sep}/k_e and C independent of m; the ‖x(0) − x_p(0)‖/k_e is justified (Step 1: Θ moves by O(ε)‖e(0) − e_p(0)‖ during the transient); ρ/2 is the right start radius. The cell/plan margins must bound the per-particle reach (1 + w_min^{−1/2})Cρ of the |||·|||-ball, not ρ. |
| **(c)** | **holds** (Proposition D‴.15 and Theorem J(b) for (CL_w)). On the moving-orbit instance the true ghost exponents at k_e = 25 equal (σ₀ − λ_j(R^G_k))/k_e to ≤ 3·10⁻⁴ relative, the collective Floquet exponents λ_j(H_q)/k_e to 0.57–3.15 %. |
| "In particular …" (proviso void; D‴.7 settled) | **holds in substance**: D‴.18 is an m-uniform replacement of Corollary D‴.14(ii)'s proviso, and with per-particle starts Corollary D‴.7's transients become m-uniform ((1/m)Σ‖u_i‖² ≤ max‖u_i‖² ≤ (1 + w_min^{−1/2})²|||u|||²). Wording: Theorem J's own orbit in its unnormalised ball still needs k_e ≥ √m·C_q/ρ₀ (D3g-3, D3g-8). |
| Proof, structural fact (1) | **holds**: the first-order term of f̂ in the separations vanishes (remainder 0.057 s²), the remainder after the second-order term is 0.0068 s³. "Steps 1–2 involve Θ only through r_Θ" is imprecise (the ζ-forcing contains Θ̇; per particle O(ε), so the constants stay m-independent) (D3g-8). |
| Proof, structural fact (2) | **holds**: the off-diagonal blocks of the drift Jacobian are exactly O(1/m) (m·max entry 0.502, 0.510, 0.262 at three configurations, the same for every m); for a fixed one-particle displacement m·N_i (i ≠ 0) is m-independent to seven digits; C_N for random per-particle displacements is bounded (≤ 0.068) with no growth in m; the real loop's one-period map has a second-order remainder |||R|||/|||u|||² ≤ 0.0102 at m = 6, 30, 150, largest at m = 6 — the "inherit" claim at the level of the full map. "∂²_{ϑ_jϑ_l} vanishes" holds for ⟨H⟩ only (D3g-8). |
| Proof, structural fact (3) | **holds**: at the zero-residual matching prior, at the shifted prior (R ≠ 0) and at w = (7, 5, 4)/16 under the σ-modification (non-uniform weights), for m up to 150/160: class means and separations decouple to the finite-difference floor (≤ 3.6·10⁻¹⁰), the separation action is I ⊗ (GR_k − σ₀I) (to 5·10⁻¹⁰), exactly G⁻¹-self-adjoint, the collective compression has H_w's spectrum, and |||I + (T₀/k_e)D⟨H⟩||| = 1 − (T₀/k_e)min{λ_coll, λ_sep} exactly (0.998222, 0.998574, 0.999353), attained by sampling. D⟨H⟩ is m-uniformly Lipschitz in |||·|||, but with a constant 2–5 times the remainder constant C_N, not C_N (D3g-2). |
| Steps 4 and 6 as adapted | **hold up to the radius bookkeeping** (D3g-2): ρ := min{λ₀/(4C_N), margins} conflates two constants and two norms. The m-uniformity itself is sound: |||·||| is equivalent to the per-particle max-norm with the constants 2 and 1 + w_min^{−1/2}, and every map built in Steps 1–4 is of mean-field type. |
| *Status* "proved as an adaptation" | **acceptable once D3g-2 is fixed** — a sketch at Theorem J's level of detail, whose m-uniformity rests on the norm equivalence and the mean-field structure (worth saying). The description of the tests is accurate (with D3g-7's caveats). |
| **Remark D‴.19** | **holds in part.** ρ₀ is m-independent in the unnormalised norm (the whitened Hessian's Lipschitz constant there is 0.004–0.093, m-independent at a non-zero residual), but Theorem J's threshold in that norm grows like √m, which the remark omits. The no-mean-square-ball argument holds for m large; at m = 15 the shifted weights are a saddle. "The basin at the level of the empirical measure is larger" is unproved (D3g-3). |
| *Numerics* (A), inside the cells | **numbers hold, reading partly overstated.** All 24 rows and the ranges reproduce. The leaks are identical to 0.0; no particle leaves its cell (confirmation 3). The m = 3 explanation fits one of three runs, C depends on the horizon for two m = 15 runs, and the runs are 12–640 times outside the proof's ball (D3g-4, D3g-5). |
| *Numerics* (A), ρ = 4.5 | **numbers hold, reading fails** (D3g-1, **major**; D3g-6). |
| *Numerics* (B) | **numbers hold**; the fitted rates are pre-asymptotic mixtures and the Floquet exponents are the real test of the rate (D3g-5). The reference is as described. |
| `experiments/exp_uniform_basin.py` | **holds.** The batched step equals `exp_ot_real_loop.real_run` (same dt, same target refresh) to ≤ 4.4·10⁻¹⁶ over 10⁴ steps in four cases, including target switches across a boundary. One system alone and the same system inside a batch of five are bitwise equal. Records and statistics are as described; regeneration is byte-identical. Fixes: D3g-1, D3g-5. |
| `results/chaos/uniform_basin.md` | **holds with cosmetic fixes**: the summary's "ratios 0.952–1.435" compares m = 3 with σ₀/k_e (D3g-5); the ρ = 4.5 rows print a "rate" and "C" for distances that do not tend to zero (D3g-8). |
| `tests/test_uniform_basin.py` | **pass (3/3, 1.04 s); none vacuous**, but weaker than their names. Test 1 is meaningful (the finite-difference linear part is accurate), but at zero residual the per-particle quadratic coefficient vanishes and its collective part is the lift's invariance. Test 2's "same at every m" holds by construction at zero residual. Test 3 would pass at 0.8σ₀ (D3g-7). |
| §0 row, §7 row, §8 bullet, RESULTS | **numbers hold.** §0's "attracting every solution in it at the rate min{·}/k_e" and its unnormalised clause (D3g-3, D3g-8); §7's "constants 1.00–1.29" for m ≥ 15 is 1.03–1.29; the §8 bullet's "with admissible weights" means margins λ₀ on a compact 𝒲; RESULTS is accurate (it does not make the "lifts" claim). |

## Confirmed issues

**D3g-1 (major; D‴.18 *Numerics* (A), PARADIGM l. 2303–2306).** "At ρ = 4.5 some particles start in another cell: the nearest-atom leak keeps
them there and ends with the counts (2, 0, 1), (6, 5, 4), (23, 19, 18), (87, 70, 83) and W²(μ_Θ, ν) = 12.3, 2.8, 2.3, 1.9 — the lifts of
other weights, as Remark D‴.19 says". Three of its clauses are false.

(i) *"keeps them there".* The initial configurations, rebuilt from the seeds, have 1, 3, 13 and 55 particles outside their class cell
(20–26 % of each G⁻¹-sphere of radius 4.5 lies in another cell; none at ρ ≤ 3.13). At m = 60 and 240 the counts change during the run:
3 and 12 changes between records, the last at t = 420 and 285. 3 and 14 particles end in a different cell from their starting one,
and 1 and 7 of these had started in their own cell. At m = 240 the starting counts are exactly (80, 80, 80) = q (the 55 misplaced
particles balance), so the end counts (87, 70, 83) come entirely from crossings during the transient. At m = 3 and 15 the counts are
constant.

(ii) *"the lifts".* At t = 7995 no cloud has collapsed:
- within-class spread (max): 0.78 (m = 3, class 0); 0.94/0.99/0.40 (m = 15); 0.81/1.19/0.49 (m = 60); 0.72/0.84/0.53 (m = 240); rms 0.24–0.66 for m ≥ 15;
- the particles sit up to 1.21, 1.25 and 0.84 (m = 15, 60, 240) from the lifts Θ*_m(w′) of their realised weights;
- the averaged drift is 2.5·10⁻³–6.9·10⁻³ per slow unit.

Two realised weight vectors are inadmissible: (6, 5, 4)/15 with σ₀ − max κ_R = −0.00348 and (23, 19, 18)/60 with −0.00161. Their lifts
are saddle orbits (Proposition D‴.15), which the systems cannot approach. This is round 12's D3e-1 again. (2, 0, 1)/3 empties a class:
the two class-0 particles sit 0.78 from their mean, the fit error is 0.35, and Newton finds no collapsed two-atom stationary point with
weights (2/3, 1/3) near the prior atoms (residual 5.9·10⁻²). Only (87, 70, 83)/240 is admissible (+0.00167).

(iii) *"as Remark D‴.19 says".* The remark's mechanism is one particle across a boundary at vanishing mean-square displacement. These runs
displace every particle by 4.5 (mean square 20.25), and at m = 240 the change of weights comes from crossings during the run, not from
the initial misplacement.

*What holds:* the end counts differ from q and are frozen after t = 420 (m = 60) and t = 285 (m = 240), and at m = 3 and 15 from the
start. The W² at the end equals the W² of the lifts of the realised weights to within 2 % at m = 60 and 240 (2.303 vs 2.345, 1.921 vs
1.928) and to 13 % at m = 15 (2.809 vs 3.213): the cross-cluster transport of the weight mismatch. Theorem D‴.18 is not affected.

*Fix:* "ends at other weights (counts …, frozen after t ≈ 400, W² that of the lifts of these weights), with clouds not yet collapsed
(spread …), and for m = 15 and 60 at inadmissible weights whose lifts are saddles". Have `report()` print, per system, the start and end
cell counts, the number of cell changes, the within-class spread and the splitting margin of the end weights (as round 12 asked of
`exp_ot_real_loop.py`). To illustrate Remark D‴.19, add a run in which one particle is pushed across a boundary (m = 60, 240).

**D3g-2 (minor; D‴.18 (a), proof l. 2266–2274).** The radius bookkeeping of the adapted Steps 4 and 6 is inconsistent in three places.
None of them touches m-independence.

(a) *Uniqueness ball.* The contraction is shown on the |||·|||-ball of radius ρ, while (a) claims uniqueness in the per-particle ball
of radius ρ. Since |||u||| ≤ (2 − 2/n_max)·max_i‖u_i‖, sharp (one particle of the largest class at +e, the others at −e: ratios 1.60,
1.67, 1.95, 1.98 at n_max = 5, 6, 40, 96), that ball lies only in the |||·|||-ball of radius ≈ 2ρ. (b)'s ρ/2 is right for this reason,
and (a) should also read ρ/2, or the contraction should be run on 2ρ. The sentence "the norm equivalence max_i‖u_i‖ ≤ (1 + w_min^{−1/2})|||u|||
translates the |||·|||-balls into per-particle balls" uses the wrong direction for the start and uniqueness sets: they need the converse
|||u||| ≤ 2max_i‖u_i‖, which is not stated.

(b) *Cell and plan margins.* With ρ := min{λ₀/(4C_N), the cell and plan margins}, the |||·|||-ball of radius ρ reaches the per-particle
distance (1 + w_min^{−1/2})ρ (sharp: 2.73, 3.24, 3.45, 8.07 at w_min = 1/3, 0.2, 1/6, 0.02), and the solutions of (b) reach C times that.
For the leak's targets to stay locally constant, the margins must bound C(1 + w_min^{−1/2})ρ, not ρ.

(c) *Two constants called C_N.* C_N is used both for the remainder bound ‖N_i(u)‖ ≤ C_N|||u|||² (fact (2)) and for "D⟨H⟩(Y) differs
from D⟨H⟩(Θ*_m(w)) by at most C_N|||Y − Θ*_m(w)|||" (fact (3)). For a quadratic remainder the second constant is at least twice the
first (DN(u)u = 2N(u)), and up to four times with the polarisation constant of a max-type norm. On the instance, for the same
displacement, the |||·|||-norm of D⟨H⟩(Θ* + u) − D⟨H⟩(Θ*) per |||u||| is 0.13–0.23 (sampled lower and block upper bound) against a
remainder ratio of 0.068 (matching prior, m = 6, per-particle 0.05), and 0.076–0.136 against 0.029 at the shifted prior. With
ρ = λ₀/(4C_N) and C_N the remainder constant, the Jacobian perturbation can reach λ₀/2 to λ₀ on the ball, so κ < 1 is not guaranteed.

*Fix:* let C_N denote the |||·|||-Lipschitz constant of D⟨H⟩ on the ball (then ‖N_i‖ ≤ (C_N/2)|||u|||² follows by integration); set
ρ := min{λ₀/(4C_N), margins/(C(1 + w_min^{−1/2}))}; state uniqueness in the per-particle ball of radius ρ/2; state both directions of the
norm equivalence.

**D3g-3 (minor; Remark D‴.19 l. 2277–2284; §0 row clause "in the unnormalised norm of Theorem J the basin is m-independent too but shrinks per particle like m^{−1/2}").**

(a) *k₀ in the unnormalised norm.* "In the unnormalised norm ‖·‖_{G_m⁻¹} of Theorem J the basin radius ρ₀ is m-independent too" is true
of ρ₀. The whitened Hessian of J_m is Lipschitz in that norm with an O(1) constant: 0.004–0.093 for one-particle displacements and
0.005–0.061 for spread ones, and m-independent at a non-zero residual (0.016, 0.012, 0.013 at m = 16, 48, 160). Theorem J's threshold is
not m-independent in that norm. Every periodic orbit near Θ*_m sits √m·O(1/k_e) away in that norm (the lift: √m·sup_t‖Θ^Γ − Θ*(w)‖_w),
and Ã_ε(Θ*) = O(ε√m) there, so the Lipschitz implicit-function step needs k_e ≥ C√m/ρ₀. This is exactly Corollary D‴.14(ii)'s proviso,
and the actual reason the per-particle norm is needed. The remark and the §0 clause omit it.

(b) *Admissibility of the shifted weights.* "The solution converges to the lift of different weights however small its mean-square
displacement" holds for m large, where the shifted weights are admissible. For small m they need not be: at m = 15 on the matching prior,
moving one particle out of a class of five gives (6, 5, 4)/15, margin −0.00348, a saddle lift. The solution then splits instead; it
still does not converge to Γ_q, which is all the remark needs.

(c) *Empirical-measure basin.* "The basin at the level of the empirical measure is larger" is stated as a fact but is unproved; the §8
bullet correctly lists it as open.

*Fix:* add "but its threshold k₀ grows like √m (the orbit itself is √m/k_e away in that norm)"; add "for m large enough that the
shifted weights are admissible"; replace "is larger" with "appears larger in the numerics (open, §8)".

**D3g-4 (minor; D‴.18 *Numerics* (A), experiment docstring "The prediction of Theorem D‴.18(b): every d_∞ inside the per-particle ball decays at the rate …").**
The numerics probe a ball far outside the proved one, and the text does not say so. With the tests' remainder constant C_N = 0.036, the
proof's ρ = λ₀/(4C_N) is 0.050 in |||·|||; with the Jacobian constant of D3g-2(c) (0.13–0.23) it is 0.008–0.014, and the start ball of
(b) has per-particle radius ρ/2. The runs use per-particle ρ = 0.3–2.5 (|||u(0)||| up to ≈ 2ρ), which is 12–640 times the proved
radius. The clean contraction observed there is evidence that the basin extends to the cell margins, a statement the theorem does not
make. *Fix:* give the proof's ρ on the instance, and present ρ = 0.3–2.5 as probing beyond it.

**D3g-5 (minor; D‴.18 *Numerics* (A) "at 2.95·10⁻⁴–4.06·10⁻⁴ for m = 3 … (the slowest mode is then the collective one …)", "one rate and one constant across a factor 80 in m", (B) "ratios 0.96–1.16"; `uniform_basin.md` summary line).**

(a) *m = 3.* The explanation fits one run of three. On the eigenvectors of the collective block, the stored final deviation at ρ = 2.5
carries amplitude 0.086 on the third collective mode (rate 4.08·10⁻⁴ averaged, 4.00·10⁻⁴ Floquet) and 0.029 on the slowest. The fitted
4.06·10⁻⁴ is therefore the third mode's rate (quarter-window rates 1.21, 1.41, 1.44, 1.43μ: a clean exponential, but of the wrong mode).
At ρ = 0.3 the three slowest modes have comparable amplitudes (0.0060, 0.0059, 0.0039), and the decay over the second half is not
exponential (quarter rates 1.88, 1.56, 1.43, 1.09μ; rms log-residual of the fit 2.8·10⁻², against ≤ 7.3·10⁻³ in every other inside-cell
run). Only ρ = 1.0's 2.948·10⁻⁴ is the slowest collective rate (Floquet 2.9355·10⁻⁴); the slowest mode would take over about 10⁴ time
units after the horizon. The md's "ratios 0.952–1.435" compares m = 3 with σ₀/k_e, which is not m = 3's prediction: without a
separation block it is λ_coll/k_e = 2.952·10⁻⁴.

(b) *C depends on the horizon* where the fitted rate is below μ. In the other 12 runs with m ≥ 15 (Parts A and B) the maximum of d_∞e^{μt} is attained at t = 530–3365 (`check_C_argmax.py`). For m = 15, C(T) at T = 2000, 4000, 6000, 8000 is:

| run | C(T) | last-quarter rate | maximum of d_∞e^{μt} |
|---|---|---|---|
| (A) ρ = 0.3 | 1.113, 1.114, 1.122, 1.162 | 0.938μ | at the last record |
| (A) ρ = 1.0 | 1.168, 1.228, 1.257, 1.266 | 0.987μ | at the last record |
| (B) ρ = 0.3 | 1.139, 1.189, 1.216, 1.232 | 0.971μ_B | at the last record |

(c) *Part B fits.* The fitted ratios 0.96–1.16 are mixtures that have not reached the asymptotic regime. The true separation exponents at
k_e = 25 (ghosts along the shot 3-atom orbit) are 2.270·10⁻⁴ (class 2), 2.597·10⁻⁴ (class 1), 2.740·10⁻⁴ (class 0), then 2.796–3.29·10⁻⁴,
all equal to (σ₀ − λ_j(R^G_k))/k_e to ≤ 3·10⁻⁴ relative. At the end of every run d_∞ is carried by class 2, but class 1 is close (m = 15,
ρ = 1: 0.136 against 0.130) and the next directions of class 2 decay at 2.80 and 3.04·10⁻⁴. So 1.16 is a mixture, and 0.96 (m = 15,
ρ = 0.3, quarter rates 0.90, 0.91, 0.95, 0.97μ) lies below every linear exponent, as the bound (1 − κ)μ allows. The inverse gaps
(≥ 3·10⁴ time units) exceed the 4000-unit fit window.

*Fix:* quote the Floquet exponents as the test of the rate; present the fits as consistent with (1 − κ)μ; describe the m = 3 runs as
mode mixtures; give C with its horizon; in the md, compare m = 3 with λ_coll/k_e.

**D3g-6 (cosmetic; D‴.18 *Numerics* "W² = 0.010–0.23, still falling (the transit)", l. 2305–2306; RESULTS "(still falling)").** Under
the OT leak at ρ = 4.5 the counts reach q at t = 2555 (m = 3) and by t = 2695, 1610, 1445 (last change) for m = 15, 60, 240. Over
t > 6000, W² decays at exactly 2σ₀ per slow unit for m ≥ 15 (log-slopes −0.01411, −0.01408, −0.01417 against 2σ₀ = 0.01415). This is the
leak's contraction of every particle's remaining displacement (ρ²e^{−2σ₀s} = 0.219 at s = 320, against the measured 0.156–0.226), not the
transit, which has ended. The OT finals are also the lift of q under a different class map: 4, 15 and 53 particles end outside the cell
of their initial class, because the plan relabels. "Restores q" holds for the measure, not for π₀. *Fix:* "still falling at the leak's
rate 2σ₀ (the transit having ended by t ≈ 2700)".

**D3g-7 (cosmetic; `tests/test_uniform_basin.py`, the Status, RESULTS).** All three tests pass in 1.04 s and none is vacuous, but each is
weaker than its name.

- **Test 1.** The finite-difference linear part (h = 10⁻⁵) is accurate, so the remainder test is meaningful: with the same directions,
  C_N converges as the displacement shrinks from 0.2 to 0.05 to 0.01 (0.0121, 0.0109, 0.0106 at m = 30; 0.00659, 0.00658, 0.00657 at
  m = 150), so the remainder is genuinely quadratic. But at zero residual the O(1) per-particle coefficient ∂²_{ϑ_iϑ_i} ∝ ⟨r∇³Φ̃⟩ vanishes. For
  a fixed one-particle displacement, |N_0| = 9.252·10⁻⁴·6/m exactly, and m|N_i| = 1.677·10⁻³ for the other particles of its class, the
  same to seven digits for m = 6–600. Only the mean-field terms of N are exercised: C_N = 0.0357, 0.0121, 0.0066, 0.0041, decreasing
  because random class means are O(m^{−1/2}); slack 2.8. The collective part (0.0636, equal across m to 3·10⁻⁹) is exact by Corollary
  D‴.14(i): it checks the lift's invariance, not D‴.18.
- **Test 2.** The "one particle" separation perturbs two particles (±v). At zero residual the response is the leak alone (−σ₀u₀ to
  10⁻¹³), so "the same at every m" holds by construction. The tolerances 10⁻⁹ and 10⁻⁸ face exact zeros and 10⁻¹⁴. At the shifted prior
  the response (GR₀ − σ₀I)u₀ differs from −σ₀u₀ by 8·10⁻⁵ and is m-independent to 2·10⁻¹⁴: that is the test the name promises.
- **Test 3.** The constants are 0.999, 0.999, 0.994 against a tolerance of 1.5, and |||u(s)||| decays at 1.000σ₀ for every m. With μ = σ₀
  the constants are still ≤ 0.9996, so the factor 0.9 is pure slack, and a pure exponential at 0.8σ₀ would pass (e^{0.1σ₀·400} = 1.33 < 1.5).
  At the shifted prior (λ_sep = 0.80σ₀) the averaged loop contracts at 1.00–1.15λ_sep with C ≤ 0.998 for m = 6, 30, 150.

*Fix:* run tests 1–3 at the shifted prior, where R ≠ 0 and λ_sep < σ₀; tighten test 3 to μ = λ_sep with C ≤ 1.05.

**D3g-8 (cosmetic; wording in the statement, proof, rows and md).**

- **(a) Fact (1).** "Steps 1–2 … involve Θ only through r_Θ" is imprecise: the ζ-forcing −∂_Θe_Θ·Θ̇ and the remainder R = ∂_eF·ζ involve
  the particles' drifts. These are O(ε) per particle, so the constants are m-independent, but by (2), not (1); likewise "Lip_ζP_Y = O(ε) …
  m-independent by (1)".
- **(b) Fact (2).** "∂²_{ϑ_jϑ_l} vanishes for distinct j, l ≠ i because r is linear in the Φ̃'s" holds for ⟨H⟩ but not for H₂, which
  contains ∂_xr_Θr_Θ, quadratic in f̂, with entries O(1/m²). The bound survives.
- **(c) Dependence list.** "the plant, the reference, the features, σ₀, G, 𝒲 and λ₀" omits the leak's prior (its cell and plan margins
  enter ρ) and R; the Status's "the only inputs are those of Theorem J and the margins λ₀" also omits D‴.14(i) and D‴.15.
- **(d) "The lift of Γ_q is Theorem J's orbit at Θ*_m(q) for every commensurate m".** D‴.18 is the m-uniform replacement. Theorem J's
  own orbit, unique in its unnormalised ball, coincides with the lift wherever Theorem J applies with ρ₀(m) ≤ ρ/2, but in its own norm it
  still needs k_e ≥ √m·C_q/ρ₀ (D3g-3). Say "the proviso is replaced by".
- **(e) Rows and bullet.**
  - §0 row: "unique in a per-particle ball of radius ρ, attracting every solution in it at the rate min{λ_coll, λ_sep}/k_e" — but (b)
    attracts the ρ/2-ball, at (1 − κ)min{·}/k_e.
  - §7 row: "for m ≥ 15 … with constants 1.00–1.29" — for m ≥ 15 they are 1.03–1.29; 1.00 is m = 3.
  - §8 bullet: "with admissible weights" means margins λ₀ on a compact 𝒲.
- **(f) md table.** The ρ = 4.5 rows print a "fitted rate" and "C" for distances that do not tend to zero (nearest-atom d_∞ at the end:
  2.50, 1.16, 1.27, 0.89).

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_uniform_basin.py --from-raw` on a scratch copy rewrites `uniform_basin.md` byte-identically (md5 `4d427b58…`).
2. **Every number of the md and the paragraph** from the 1600 records per system:
   - **Part A**, all 32 rows: final counts, d_∞(0) → d_∞(end), second-half rates, C, W² start → end, rms tracking errors.
   - **Inside the cells:**
     - rates for m ≥ 15: 2.6940–2.9710·10⁻⁴ (ratios 0.952–1.050 to σ₀/k_e = 2.8298·10⁻⁴);
     - rates for m = 3: 2.9479–4.0619·10⁻⁴;
     - C: 1.000–1.291 (m ≥ 15: 1.026–1.291).
   - **ρ = 4.5:** W² at the end 12.259, 2.809, 2.303, 1.921 (nearest-atom) and 0.0101, 0.2260, 0.1647, 0.1561 (OT).
   - **Part B**, all 6 rows: rates 2.1840–2.6316·10⁻⁴ (ratios 0.962–1.159 to 2.2696·10⁻⁴), C 1.040–1.232, tracking errors 2.548–2.556·10⁻⁴
     against the reference's 2.553·10⁻⁴; horizon 8000 = 1273 periods = 320 slow units.
3. **"Inside the cells" holds along the whole run.** At ρ ≤ 2.5 no particle starts outside its class cell (0 of 318 at each ρ; initial
   half-gap margins ≥ 0.657 at ρ = 2.5), the counts are constant over all 1600 records, max_t d_∞ ≤ 2.502 < Δ_ν/2 = 3.1315, and the
   nearest-atom and OT runs are identical to 0.0 in d_∞, W², the final configuration and e_rms.
4. **The integrator.** `Batch.step` equals `exp_ot_real_loop.real_run` (with its dt and target refresh set to 5·10⁻³ and 200) to
   ≤ 4.4·10⁻¹⁶ after 10⁴ steps for (A) m = 15, ρ = 1 OT; (A) m = 15, ρ = 4.5 with both leaks (target switches); and (B) m = 15, ρ = 1. A
   system alone and the same system inside a batch of five agree bitwise. The hand-written f̂ and ∇_ϑ(eΦ̃) equal `SaturatedLiftedNetwork`'s
   formulas.
5. **Predictions and the true spectra.**
   - **Part A:** the three separation blocks are σ₀I exactly (residual 7·10⁻¹⁷), λ_coll = 0.007380, the ghost exponents equal σ₀/k_e to 10⁻⁸
     (the finite-difference floor), and the 3-atom collective Floquet exponents are 2.9355, 3.5579, 3.9986, … ·10⁻⁴.
   - **Part B:** Θ*(q) residual 4.3·10⁻¹⁷; λ_min of the separation blocks 0.006849, 0.006491, 0.005674 (classes 0, 1, 2), λ_coll = 0.007398;
     the 3-atom orbit by shooting has residual 1.4·10⁻¹⁴ and period-start offset 5.6·10⁻⁴; its ghost exponents equal the prediction to
     ≤ 3·10⁻⁴ relative and its collective exponents λ_j(H_q)/k_e to 0.57–3.15 % (Part A: 0.56–3.22 %, as round 12 found).
6. **Structural fact (3)** at three replicated configurations and m ≤ 160: decoupling, separation action, G⁻¹-self-adjointness, collective
   spectrum and the exact |||·|||-norm 1 − (T₀/k_e)min{λ} of I + (T₀/k_e)D⟨H⟩ (see the verdict). The perturbation D⟨H⟩(Θ* + u) − D⟨H⟩(Θ*)
   is m-uniformly Lipschitz: per |||u||| it is 0.007–0.13 (sampled) to 0.07–0.23 (block bound), and with a fixed per-particle pattern its
   per-particle ∞-bound is 0.099, 0.100, 0.106 (and 0.112, 0.113, 0.114 with a class-mean shift) at m = 12, 48, 192.
7. **Structural facts (1) and (2)** as in the verdicts, and the norm equivalence in both directions with sharp constants.
8. **The full one-period map (DOP853, rtol 10⁻¹¹)** at the lifted orbits, for m = 6, 30, 150:
   - lift residuals ≤ 1.7·10⁻¹³;
   - second-order remainder |||R|||/|||u|||² = 0.0040/0.0066, 0.0010/0.0034, 0.0007/0.0006 (matching prior, per-particle 0.05/0.2) and
     0.0102/0.0077, 0.0013/0.0011, 0.0015/0.0018 (shifted prior);
   - linear part |||DPu|||/|||u||| = 0.998224 = e^{−T₀σ₀/k_e} at every m (matching prior) and 0.9982–0.9984 ≤ 1 − (T₀/k_e)λ_sep = 0.998574
     (shifted prior).
9. **Contraction of the averaged loop at a non-zero residual** (shifted prior, λ_sep = 0.80σ₀): |||u(s)||| decays at 1.15, 1.006, 1.002·λ_sep
   (second half) for m = 6, 30, 150, with C = 0.998, 0.998, 0.994. This is the theorem's rate at R ≠ 0, the same for every m.
10. **The finals.** For m ≥ 15 inside the cells the |||·||| of the final deviation is dominated by the separations (0.03–0.28) over the class
    means (0.0008–0.06), so d_∞ measures the separation (leak) rate. In Part B every system keeps its classes (counts constant) and the
    reference ends within 4.1·10⁻⁴ of Θ*(q); d_∞ is a distance between two solutions, as the text says.
11. **The ρ = 4.5 bookkeeping** behind D3g-1 and D3g-6 (initial misplacements, cell changes, spreads, margins, W² of the ideal lifts, the
    2σ₀ decay of the OT W²).
12. **Tests.** 3 passed in 1.04 s; every asserted quantity recomputed (D3g-7). README's "100 tests" is right at 5ff2e0a (104 at cffb7be,
    with the 4 entropic-leak tests).

## Not checked

The 12-minute experiment was not rerun. Single systems were re-integrated only for 10⁴ steps (the integrator identity). The long-time fate
of the ρ = 4.5 nearest-atom systems at saddle weights, the convergence of C(T) for the m = 15 runs, and the take-over by the slowest mode
at m = 3 are inferred from the stored finals and the linear modes, not followed beyond the horizon. The proof's constants (k₀, κ, C) were
not derived. The uniqueness statement (a) is not testable numerically. The structural facts were checked on the averaged loop at three configurations
and on the one-period map at two, with m ≤ 160 — evidence, not proof. Theorem J's own proof, Proposition D‴.5(i), Corollary D‴.14(i) and
Proposition D‴.15 were used as given (rounds 5, 9, 12). The measure-level basin of the OT leak (open by the authors' account), dither, and
incommensurate sizes are outside the material. The concurrent session's commit cffb7be, its uncommitted edits and Proposition D‴.20 were
not reviewed.
