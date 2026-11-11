# Referee report D3c — PARADIGM §5d, the material new since round 9: Proposition D‴.9 (i)–(iii) and its Status, Remark D‴.10, the D‴.9 *Numerics* block (tables 1–4, Reading), Remark D‴.11; `src/core/prior_design.py`; `experiments/exp_prior_design.py`; `tests/test_prior_design.py`

Branch `claude/festive-bardeen-eopej2` at 0d6bb36. No repository file was modified (the `--from-raw` regeneration rewrote
`results/chaos/prior_design.md` byte-identically, md5 `ce7fd1db…` before and after). Scratch scripts and their outputs (`.out`) in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3c/`:
`check_gradJ.py` (grad_J by finite differences at random points and on/near a Voronoi boundary, single-atom equivalence, drift = −G∇J,
cell_margin against the exact distance to the boundary, weighted_energy_grad_prior, the decomposition for a pair straddling a boundary),
`recompute_from_npz.py` (every number of tables 1–4 and the prose from `prior_design_raw.npz` / `learning_raw.npz`),
`check_signs_admissibility.py` (Remark D‴.11 against the residual profiles, the admissible set of the matching prior, the margin at the
reference weights, the K-atom Hessian, the exact stationary distance per pair and the plateau extrapolation, Proposition K's R₂, test 2's
slack), `relaxation_compare.py` (observed relaxation of D in the two experiments), `saddle_runs.py` and `saddle_followup.py` (the two
inadmissible m-systems re-simulated autonomously to s = 9000 and 6000), `from_raw.out`. The full experiment (1.5 h) was not rerun; the
m-systems of (𝔄_m^ν) are autonomous (the reference enters only the diagnostics), so they were reproduced exactly from their seeds.

## Verdicts

| object | verdict |
|---|---|
| **Proposition D‴.9 (i)** | **holds with fixes.** The leak enters Theorem D‴'s proof only through Step 0's a priori bound (which goes through with dist(·, supp ν) in place of ‖· − ϑ⁰‖: the upper Dini derivative of the min of quadratics along the flow is bounded by the active branch's, so the stated R_Θ is right) and through Step 2's curvature σ₀G⁻¹ on the segments [ϑ̄_i, ϑ_i]; the latter needs the segment inside one cell, which the stated hypothesis guarantees only for s ≥ s₁ while (2.1)/(5.2) are also used on [0, s₁) — a hidden hypothesis (D3c-1, major, trivial fix). The "equivalently" is a "sufficient" (D3c-7). Theorem J applies at the replicated configuration with ρ₀ below the cell distance (unstated, cosmetic). |
| **Proposition D‴.9 (ii)** | **holds.** r ≡ 0 gives R ≡ 0 at every ϑ, hence κ_R[μ*](A) = 0 for every A and margin σ₀; the exact representation is stationary; its K-atom Hessian is Gram + σ₀diag(q)G⁻¹ ≻ 0 (whitened/weighted eigenvalues 0.0074–0.20, all ≥ σ₀ on the instance), so the IFT applies and the O(δ) statements follow; the matching prior's stationary measure is a strict minimiser (separation blocks σ₀G⁻¹). Cosmetic: the restriction "ϱ′ below half the cell margin" is (i)'s Voronoi hypothesis, not (H); the masses q_j play no role in the leak (D3c-11). |
| **Proposition D‴.9 (iii)** | **holds with fixes.** The function-level law is the delta method with the right normalisation (‖·‖²_{L²(ρ)} = (1/N)Σ on the uniform N-point grid; the shared-label covariance (1/m − 1/M)C̃ re-derived); F and S confirmed by finite differences; "neither floor vanishes" is argued correctly (at (7,5,4)/16 the stationary amplitude of the over-weighted atom is 2.446, between the fit's 2.286 and the prior's 3.0). "Both sensitivities grow with the amplitudes" is a trend of the amplitude rows, not a law (D3c-6); E D_∞ quoted as an equality (D3c-12). |
| *Remark D‴.10* | **holds with fixes.** "margin above σ₀ … all tested weights admissible" for the shrunk prior: admissible yes, above σ₀ only at (7,5,4)/16 (D3c-2); "3.6 times faster" is the ratio of the theorem's 1/κ, the observed e-folding ratio is ≈ 3 (D3c-10). |
| *Numerics* block and *Reading* | **holds with fixes.** Every number of tables 1–4 and the prose reproduces from the raw data (see confirmations). Over-readings: "margins at or above σ₀" (D3c-2), the saddle explanation of the ratios 0.91/0.81 — the plateaus are not left in D (D3c-4), seed-mean statistics quoted as maxima/values (D3c-5), "3.6 times sooner" (D3c-10), table 4's Jensen explanation (D3c-13). |
| *Remark D‴.11* | **holds with fixes.** ∂²_{cc}Φ̃ = −(a/s²)(1 − (x − c)²/s²)φ is right, and for the matching prior the mechanism is exactly as stated (the top eigenvector of R^G at the over-weighted atom is the pure centre direction). The two other attributions are wrong on the instance: the σ-modification's inadmissible seed has R_cc < 0 at its over-weighted atom (it still undershoots at the centre; the unstable direction is amplitude–width), and Proposition K's R₂ is not carried by the region "between the bumps" (D3c-3). "The shrunk prior undershoots everywhere" is false at three of the four weights (D3c-2). "Loses definiteness" should be "loses margin" (D3c-3). |
| `src/core/prior_design.py` | **holds.** grad_J = FD of J to 5.6·10⁻¹⁰ relative away from boundaries; drift = −G∇J exactly; the single-atom loop equals `AveragedSlowLoop` to 8.9·10⁻¹⁶ (also for ϑ⁰ ≠ 0); J = mL̄ + σ₀ΣV_ν with V_ν = ½min‖·‖²_{G⁻¹}; `weighted_energy_grad_prior` = ∇𝒢_w/w_k to 6·10⁻¹²; `function_sensitivity` = FD of f̂_{Θ*(w)}; `cell_margin`/2 ≤ exact distance to the boundary (equality at the atoms). On a Voronoi boundary grad_J is the one-sided derivative of the active branch (as documented). Docstrings: the module says the matching prior "lowers S" (opposite of the table, D3c-8) and `cell_margin` says the gap criterion is the hypothesis (D3c-7). |
| `experiments/exp_prior_design.py` | **holds with minor fixes.** `report()` reproduced from the raw data (regeneration byte-identical, 1 s); `design_row`, the W_TEST block, `coupled_run`, `window_stats`, `slope_se`, the per-seed laws and section 4 do what they say. Fixes: the "max_s κ_R(r_ref)" and the fit/drift endpoints are statistics of the seed-mean curve (D3c-5); the per-seed check should also report the exact law, which it already has the ingredients for (D3c-4). |
| `tests/test_prior_design.py` | **passes (3/3 in 0.7 s) and tests what its docstring says**, except that test 2's T_B assertion is vacuous (T_B < 0 against a positive bound) and test 3's `fm < 4·fs` is loose (D3c-9). Test 1's FD rows are 0.6–1.75 whitened units from the nearest boundary (safe). |

## Confirmed issues

**D3c-1 (major; statement and Status of (i)).** The no-boundary hypothesis is imposed "for s ≥ s₁", but the proof of Theorem D‴ uses Step 2's
curvature bound (2.1) and the per-particle bound (5.2) on all of [0, S], with a(s) ≥ σ₀ − L₂^G(F_ρ + K_Φ) on [0, s₁) feeding Λ₁. With the
nearest-atom leak, a pair (ϑ_i, ϑ̄_i) in two different cells adds to ½dD/ds the term (σ₀/m)⟨δ_i, ϑ⁰_{j(ϑ_i)} − ϑ⁰_{j(ϑ̄_i)}⟩_{G⁻¹}, which is
≥ 0 (add the two cell inequalities ‖ϑ_i − ϑ⁰_j‖² ≤ ‖ϑ_i − ϑ⁰_{j̄}‖², ‖ϑ̄_i − ϑ⁰_{j̄}‖² ≤ ‖ϑ̄_i − ϑ⁰_j‖²): the leak of a min of quadratics is
*expanding* across the boundary (V_ν is semiconcave; its distributional Hessian has a negative singular part on the boundaries). Confirmed
numerically: for a straddling pair, `decomposition` returns half_dD = −3.061·10⁻² (the proof's convention) while the true ½dD/ds from the
actual drifts is −7.261·10⁻³; the difference +2.334·10⁻² equals the extra term to all printed digits. Its size is up to σ₀·diam(supp ν)·D^{1/2}
(σ₀‖ϑ⁰_1 − ϑ⁰_2‖_{G⁻¹} = 0.051 on the instance) — of the D^{1/2} type of T_C but O(1) instead of O(m^{−1/2}) — so it cannot be absorbed in
a(s), in Λ₁ or in the sampling term, and "the same constants" fails on [0, s₁) unless the segments stay in one cell there too. The same
hypothesis is needed for the well-posedness of (𝔄₀^ν) (the vector field is discontinuous on the boundaries, so Dobrushin's argument and
the homeomorphism property of the characteristic flow need supp μ_s to avoid them for all s). *Fix:* require "N_{ϱ′}(supp μ_s) meets no
Voronoi boundary for every s ≥ 0" (the warm-start numerics have s₁ = 0 anyway), state in Status (i) that this is where the local constancy
is used, and add that Theorem J's basin radius ρ₀ at the replicated configuration must be taken below min_k dist_{G⁻¹}(ϑ*_k, ∂cell_k)
(= half the gap at the atoms, 3.13 on the instance) so that 𝒦 = {‖Θ − Θ*‖ ≤ 2ρ₀} meets no boundary.

**D3c-2 (minor; Numerics ×2, Remark D‴.10, Remark D‴.11, `results/chaos/prior_design.md`).** "the shrunk prior keeps all four admissible with
margins at or above σ₀" is false for three of the four rows: the margins are +0.00733, +0.00570, +0.00536, +0.00663 against σ₀ = 0.00707;
only (7, 5, 4)/16 is above σ₀. At (9, 4, 3)/16, (10, 3, 3)/16 and (0.5, 0.3, 0.2) the shrunk prior's over-weighted first atom *overshoots*:
f̂(c₁) = 1.0185, 1.0271, 1.0059 against f(c₁) = 0.9989, r(c₁) = −0.020, −0.028, −0.007, κ_R(atom 1) = +0.0014, +0.0017, +0.0004. Hence
Remark D‴.11's "the shrunk prior undershoots everywhere" is also false (it holds at w̄ and at (7, 5, 4)/16 only), and Remark D‴.10's
"amplitudes ×0.5–0.7: margin above σ₀ … all tested weights admissible" conflates the design-table margin at w̄ (which is above σ₀) with
the splitting margins at the tested weights. *Fix:* "all four admissible, with margins +0.0054 to +0.0073 — above σ₀ at (7, 5, 4)/16, slightly
below it at the other three, where the over-weighted first atom overshoots by 1–3 % and κ_R ≤ +0.0017".

**D3c-3 (minor; Remark D‴.11, attributions).** The centre-overshoot mechanism is exact for the matching prior: at the four table-2 weights
the over-weighted atom 1 has f̂(c₁) = 1.069–1.100 against f(c₁) = 0.999, R_cc(phys) = +0.078 to +0.095, and the top eigenvector of R^G is
(0.01, 1.00, −0.02) in (a, c, s) — the pure centre direction; under-weighted atoms have R_cc < 0. The two other attributions fail on the
instance. (a) Theorem D‴'s inadmissible σ-modification seed (9, 4, 3)/16: the over-weighted kernel still *undershoots* at its centre
(f̂(c₁) = 0.743 against f(c₁) = 0.953, r(c₁) = +0.21, R_cc(phys) = −0.084 < 0, as at all four tested weights of the σ-modification); its
positive eigenvalue +0.0093 lives in a mixed amplitude–width direction (eigenvector (0.67, −0.37, −0.64)) — the broad kernel (s = 0.875
against the bump's 0.5) overshoots in its tails, through the cross-curvature ⟨r φ(x − c)²/s³⟩, not at the centre. (b) Proposition K's
collapsed atom: R_cc(phys) = +0.0553, of which the region "between the bumps" (0.5 ≤ |x − c| < s) contributes +0.0020; the centre overshoot
over the negative bump (|x − c| < 0.5: +0.0276, r(0) = −1.01) and the *under*shoot of the outer bumps beyond the inflection points
(|x − c| ≥ s, where 1 − (x−c)²/s² < 0: +0.0257) carry it, and the top eigenvector (0.71, −0.32, 0.63) is again mixed. (c) "the separation
block σ₀G⁻¹ − R loses definiteness" overstates: an over-weighted atom gets κ_R > 0 (loses *margin*), splitting needs κ_R > σ₀. On the
matching prior this happens only for the first (largest-amplitude) bump: excess ε on atom 1 alone becomes inadmissible between ε = 0.03
(+0.00111) and 0.04 (−0.00031), while ε = 0.08 on atom 2 or 3 alone stays admissible (+0.0021, +0.0032); the realised (6, 4, 6)/16 has
κ_R = +0.0065, +0.0019 at two over-weighted atoms and is admissible (+0.00056). *Fix:* keep the mechanism for the matching prior (where it is
verified), drop or correct the attributions (a), (b), and say "the admissible set is smallest in the direction of the largest bump".

**D3c-4 (minor; Numerics, the inadmissible pairs).** "the two inadmissible pairs are among them (0.91 and 0.81), their plateaus being saddles'
that a longer horizon would see left" (and `report()`'s "would be left on a longer horizon"). Re-simulated autonomously from the same
clouds: at s = 300 both m-systems sit on the replicated saddle (whitened Hessian of J_m^ν: 25 eigenvalues exactly −0.00296 at m = 64, 5 of
−0.00052 at m = 16, i.e. σ₀ − κ_R^max, confirming D‴.5(i)'s block splitting), and they do leave it along the unstable eigenspace at exactly
the predicted rate (the m = 16 projection grows by e^{0.00052·600} = 1.37 between s = 300 and 900, measured 0.0720 → 0.0985). But the escape
is the local splitting of Proposition K(i), not a migration between basins: the over-weighted cluster separates into two sub-clusters along
the centre direction inside its own Voronoi cell (m = 64: 15 particles at c ≈ −1.62 and 11 at c ≈ −1.32 around the saddle's −1.49, 3.1
whitened units from the nearest boundary), the cluster weights by nearest atom stay (26, 17, 21) and (6, 5, 5) up to s = 6000 and 9000, the
rms whitened displacement from the saddle saturates at 0.25 (m = 64, by s ≈ 2000) and 0.11 (m = 16, by s ≈ 9000), J_m is lower by 2.1 % and
0.4 %, and the m-system's own κ_R falls from 0.0099 to 0.0073 (m = 64). The coupled distance does not leave its plateau: D against the
reference's replicated configuration Θ*(w^M) goes 0.4356 → 0.4300 (−1.3 %) at m = 64 between s = 300 and 6000 (details in confirmation 12).
Proposition D‴.5(i)'s "the m-system leaves the configuration by migrating particles between basins" and Remark D‴.2(iii)'s floor are
therefore not what happens at these two saddles; the (9, 4, 3)/16 σ-modification seed of round 9 did migrate.
The ratios 0.91 and 0.81 are the linearisation remainder at the two largest weight fluctuations (|w^m − w^M|_∞ = 0.057, 0.084): the exact
form of D‴.5(ii), Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖², gives measured/exact 0.995 and 0.994 at s = 300 (0.987–0.998 over all nine pairs;
0.995–1.000 for the geometrically extrapolated plateaus). This contrasts with the σ-modification seed of round 9 ((9, 4, 3)/16, D from 1.29 to
5.13 by basin migration). *Fix:* report the exact law for the nine pairs (as the D‴ numerics already do), and describe the saddles' fate
as a soft split inside the basin with D essentially unchanged; `report()` can compute it with `weighted_stationary_atoms_prior` at w^m and w^M.

**D3c-5 (minor; Numerics and `report()`).** "with the reference's own residual at most 0.00488 < σ₀" is the maximum of the seed-*averaged*
κ_R curve (`kR_own = np.mean(...)` then `.max()`); the per-seed maxima are 0.00634, 0.00674, 0.00155 (all at the first record, s = 0.02), so the
margin along the run is 0.0003 for seed 1 at the start. Likewise "fit error falls from 0.029 to 0.008" and "rms drift 3.6·10⁻³ → 3.0·10⁻⁵" are
seed means (per seed 0.041, 0.036, 0.011 → 0.012, 0.009, 0.003). *Fix:* take the maximum over seeds and say "seed mean" where one is quoted.

**D3c-6 (minor; statement (iii)).** "both sensitivities grow with the amplitudes of the stationary atoms" is stated as part of the proposition
but is a trend of the amplitude-scaled rows only: at the same first-atom amplitude a₁ ≈ 3.0, ‖S₁‖ is 13.01 (matching, s₁ = 0.50), 11.13
(widths ×1.5, s₁ = 0.67) and 10.39 (widths ×2, s₁ = 0.77), and widths ×0.7 has a₁ = 3.19 > 3.0 with ‖S₁‖ = 13.00. *Fix:* move it to the
numerical findings ("on the amplitude rows both floors rise with the amplitudes") or qualify it.

**D3c-7 (cosmetic; (i), `cell_margin` docstring, table header).** "(equivalently the gap between the distances to the two nearest atoms exceeds
2ϱ′)" — the gap criterion is sufficient, not equivalent: dist_{G⁻¹}(ϑ, ∂cell) = min_{j≠j₁}(d_j² − d₁²)/(2‖ϑ⁰_{j₁} − ϑ⁰_j‖) ≥ (d₂ − d₁)/2, with
equality only on the segment between the atoms (on 400 sampled points the exact distance exceeds gap/2 by up to a factor 1.97; at the atoms
they coincide, 3.34, 3.13, 3.13). *Fix:* "in particular when"; same in the docstring ("the gap must exceed 2ϱ′") and the table header.

**D3c-8 (cosmetic; `prior_design.py` module docstring).** "a prior whose atoms sit near the plant's structure … lowers the weight sensitivity S
(the floor of Proposition D‴.5)" — the opposite of table 1 (‖S_k‖ 13.01, 6.18, 8.91 against 8.38, 4.66, 8.32; floors 1.61× and 2.18× higher) and
of the test file's docstring. Stale text from the work-in-progress commit.

**D3c-9 (cosmetic; tests).** Test 2's `assert d['T_B'] <= 11.12 ** 2 * d['D'] ** 2 / 16` is vacuous: T_B = −1.301·10⁻⁴ (≈ −‖h‖² = −1.307·10⁻⁴)
against a bound of +3.636·10⁻⁴; the non-trivial content of (3.1) is |T_B + ‖h‖²| ≤ (L₂^G/2)‖h‖D (6·10⁻⁷ against 4.3·10⁻⁴ here). Test 3's
`fm < 4 * fs` (2.18 < 4) is loose. The identity `total == half_dD` is an algebraic tautology of the three brackets under the same-cell
convention (fine as a code check, but it is not a test of the leak).

**D3c-10 (cosmetic; Numerics, Remark D‴.10, Reading).** "3.6 times faster / sooner" and "relaxation time 141 against 511 (512 in Remark D‴.10
and in Theorem D‴'s numerics)" are ratios of the theorem's 1/κ, not observed relaxation times; the observed late-time e-folding times of D are
48–82 slow units (matching prior) against 110–307 (σ-modification warm runs; medians 59 against 190, ratio ≈ 3.2), and at s = 300 the
σ-modification runs were at 82–99.7 % of their plateaus against ≥ 99 % for the matching prior — the qualitative claim holds, the numbers quoted
are not the observed ones. The collective block's smallest eigenvalue (0.00738, 136 slow units) is the relevant rate for the matching prior,
not σ₀ alone. *Fix:* "the theorem's relaxation time 1/κ"; 511 vs 512 (1/0.001954 = 511.8).

**D3c-11 (cosmetic; statement).** ν carries masses q_j, but the nearest-atom leak depends on supp ν only (q enters only through the exact
representation in (ii)); one sentence would help, since the contrast with the OT leak's mass constraint is the point of Remark D‴.10. In (ii),
"κ = σ₀ for every ϱ′ below half the cell margin": κ_R[μ*](A) = 0 for *every* A (R ≡ 0 since r ≡ 0); the restriction is (i)'s Voronoi hypothesis.
The non-degeneracy used by the IFT (Gram + σ₀diag(q)G⁻¹ ≻ 0) could be stated.

**D3c-12 (cosmetic; (iii)).** E D_∞ = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) is quoted as an equality; D‴.5(ii) has it conditional on
admissibility with O(m^{−3/2}) + O(m⁻¹M^{−1/2}) (round-9 fix D3b-5).

**D3c-13 (cosmetic; table 4).** "(a mean of norms over four seeds lies below the root-mean-square)" explains only part of the ratio 0.48 at
m = 256: the RMS over seeds gives 0.84, 0.79, 0.52. The per-seed function-level law from the realised weights gives 0.83–1.24 (m = 16),
0.60–1.57 (64), 1.01–1.18 (256) — the sharper, supported statement (the m = 256 seeds have small realised fluctuations, as Theorem D‴'s
numerics already say).

## Independent confirmations (what was recomputed and agrees)

1. **Code.** `MixturePriorLoop` with one atom at 0 (and at (0.7, −0.4, 0.9)) equals `AveragedSlowLoop` in drift, grad_J and J to 8.9·10⁻¹⁶
   on cold and warm clouds; grad_J against central differences of J at six random configurations: relative error ≤ 5.6·10⁻¹⁰ (minimum
   distance to a Voronoi boundary 0.022); drift + G·grad_J = 0 exactly. On the 1|2 boundary the FD derivative is one-sided (−0.0208 / +0.0469
   along the normal, analytic −0.0208 = the active branch); 10⁻³ inside the cell the agreement is 6·10⁻¹². `weighted_energy_grad_prior`
   against FD of 𝒢_w: 6.4·10⁻¹². `cell_margin`/2 ≤ exact boundary distance on 400 points (min excess 6.4·10⁻⁴, equality at the atoms).
2. **Table 1** (all 12 rows: fit, κ_R at atoms, κ_R over N_{0.1}, margin, cell margin, ‖S_k‖, D₆₄, f₆₄) from the stored rows; D₆₄, f₆₄ and ‖S_k‖
   re-derived from the stored S, F. Margin ratio 3.616; 1/margin 141.4 vs 511.2; D₆₄ ratio 1.606; f₆₄ ratio 2.177; cell margins 1.22–6.26;
   closed-window rows exactly those named; ×0.5/×0.7 margins +0.00740/+0.00731 with κ_R negative at every atom at w̄.
3. **Table 2** recomputed from scratch (`splitting_margin_prior`): all 12 margins and all κ_R-at-atoms agree with `prior_design.md`; the nine
   realised-weight margins (+0.00056 to +0.00609; −0.00052 and −0.00296) agree; the reference weights are admissible (margins +0.00595,
   +0.00591, +0.00646 for seeds 0–2).
4. **Table 3**: window maxima and standard errors for m = 16, 64, 256, last-window means, h means; last/first 2.188, 2.067, 1.856; last/previous
   1.046, 1.036, 1.042; slope −1.228 ± 0.402 (per seed −0.65, −2.00, −1.03), function level −0.626 ± 0.212; predicted D and f columns.
   Saturation: D(300)/D(225) = 1.007–1.018 in all nine pairs; geometric extrapolation gives plateaus within 0.2–0.5 % of D(300) and e-folding
   times 48–82 slow units.
5. **Per-seed laws**: ratios 0.909, 1.002, 1.001, 1.114, 0.805, 0.831, 0.815, 0.958, 1.038 (parameter) and 0.795, 1.081, 1.121, 1.691, 0.612, 0.927,
   0.944, 1.060, 1.217 (function) — ranges 0.81–1.11 and 0.61–1.69 as quoted; the inadmissible pairs are 0.91 and 0.81.
6. **Exact law** (not in the document): measured D(300)/Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖² = 0.987–0.998 in all nine pairs.
7. **Table 4**: predicted RMS 0.04191, 0.02071, 0.009839; measured 0.03073 ± 0.0097, 0.01542 ± 0.0033, 0.004674 ± 0.0012; ratios 0.733, 0.745, 0.475.
8. **Hypothesis along the run**: seed-mean max κ_R(r_ref) = 0.00488, fixed-residual κ_R = 5.7·10⁻¹⁸; seed-mean fit 0.0294 → 0.0082, drift
   3.56·10⁻³ → 3.04·10⁻⁵ (per-seed values in D3c-5).
9. **Remark D‴.11's formula** and its unconstrained-coordinate version: R_γγ(code) = R_cc(phys)(dc/dγ)² + ⟨r∂_cΦ̃⟩d²c/dγ² to all printed digits
   at every atom and weight tested (e.g. −0.08709, +0.05006, −0.30030).
10. **K-atom Hessian** at the exact representation under the matching prior (FD, whitened and weighted): eigenvalues 0.00738, 0.00904, 0.01019,
    0.02222, 0.02818, 0.03955, 0.09569, 0.14127, 0.20229, all ≥ σ₀ — the IFT's non-degeneracy, and the slowest collective time 136 slow units.
11. **Admissibility of the matching prior**: margin vs excess ε on atom 1 alone: +0.00477 (0.01), +0.00279 (0.02), +0.00111 (0.03), −0.00031 (0.04),
    −0.00149 (0.05), −0.00392 (0.08); on atom 2 or 3 alone it stays admissible to ε = 0.08 (+0.00208, +0.00316). Patterns: (6,4,6)/16 +0.00056,
    (6,5,5)/16 −0.00052, (5,6,5)/16 +0.00372, (5,5,6)/16 +0.00441, (6,6,4)/16 −0.00133, (4,6,6)/16 +0.00375, (7,4,5)/16 −0.00446.
12. **Saddle runs** (D3c-4): m = 64 seed 1 and m = 16 seed 0 reproduced (initial weights equal the saved final ones); Hessian of J_m^ν at the
    replicated configuration has exactly n₁ − 1 negative eigenvalues equal to σ₀ − κ_R^max (−0.00296 ×25; −0.00052 ×5); unstable projection
    grows at the predicted rate; no label change up to s = 6000 / 9000; D against Θ*(w^M) 0.4356 → 0.4300 (m = 64, s = 6000). Final state of
    the m = 64 system at s = 6000: the over-weighted cluster (26 particles) has split into two sub-clusters along the centre direction,
    15 at c ≈ −1.62 and 11 at c ≈ −1.32 (one straggler at −1.43 still moving), with a ≈ 2.67–2.68 and s ≈ 0.47–0.48 (saddle: c = −1.494,
    a = 2.605, s = 0.483); the two other clusters are unchanged to three digits; whitened Hessian of J_m^ν there: one eigenvalue −0.0002, the
    rest positive, |∇J|_∞ = 6.8·10⁻⁵ (not yet converged; the split is Proposition K(i)'s separation inside the cell, 3.1 whitened units from
    the nearest Voronoi boundary); fit error 0.0409 → 0.0379, J_m 0.09704 → 0.09497, ‖f̂ − f̂(saddle)‖_{L²(ρ)} = 6.6·10⁻³.
13. **Regeneration**: `exp_prior_design.py --from-raw` rewrites `prior_design.md` byte-identically (1 s).
14. **Tests**: 3 passed in 0.70 s; test 1's FD rows are 0.60, 1.25, 0.88, 1.75 whitened units from the nearest boundary.

## Not checked

The full experiment was not rerun (the design rows and the M = 2048 reference are taken from the npz; the reference's final particles are not
stored, so D of the continued saddle runs was measured against the replicated Θ*(w^M), which matched the saved D(300) to 0.2–0.6 %).
Theorem C′'s relation to the nearest-atom leak and §8's open items were not reviewed. Theorem D‴'s proof Steps 0–9 were not re-refereed
(rounds 8–9) beyond locating where the leak enters.
