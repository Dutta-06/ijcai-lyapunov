# Referee report D3e (round 12) — PARADIGM §5d: Corollary D‴.14 (i)–(iii) with Status and *Numerics*, Proposition D‴.15 with proof, Status and *Numerics*; the §0 summary-row clauses, the §7 rows D‴.14/D‴.15, the §8 branch-3 bullet, `results/RESULTS_PARADIGM.md` (D‴.14, D‴.15); `src/core/lifted_orbit.py`; `experiments/exp_ot_real_loop.py` → `results/chaos/ot_real_loop.md`; `experiments/exp_lifted_orbit.py` → `results/chaos/lifted_orbit.md`; the last three tests of `tests/test_ot_leak.py`

Branch `claude/festive-bardeen-eopej2` at 2c0e404 (f0cfc940…). No repository file was modified: the `--from-raw` regenerations were run on
copies of the two experiment scripts in a scratch tree (`regen/`, with `src` symlinked) and the outputs diffed against the repository —
both byte-identical (`ot_real_loop.md` md5 `ff893214…`, 3 s; `lifted_orbit.md` md5 `50c437b6…`, 1 s). Scratch scripts and their `.out`
files in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3e/`:
`recompute_ot_real_loop.py` (every entry of tables 1–2 and every range of the Reading and of the D‴.14 prose from `ot_real_loop_raw.npz`;
the within-class spreads, class-mean offsets and the 'second order' decomposition of the OT clouds; the splitting margins of the six
realised weight vectors of the nearest-atom clouds and their spreads at the start and the end of the real loop),
`recompute_lifted_orbit.py` (every number of `lifted_orbit.md` and the D‴.15 prose from `lifted_orbit_raw.npz`; the relative errors
against Theorem J(b); the whitening of the collective prediction), `check_code.py` (`ClassLeakLoop.rhs` against `ClosedLoop.rhs`, `lift`,
the two bases, the invariance of S_k and of the collective subspace under the linearised vector field at a collapsed state, the m = 6
monodromy's block structure beyond the exponents, the FD floor, the third test's quantities against its tolerances),
`check_scheme_floor.py` (the integration scheme at the exact representation for dt = 2·10⁻³, 10⁻³, 5·10⁻⁴, against the continuous loop
and against a scheme without the control hold; the nearest-atom lift of seed 0, m = 15 against the true periodic orbit Γ_w),
`check_ghost_averaging.py` (the ghost monodromy by the matrix ODE along the orbit and by the period-averaged ghost Jacobian, against
the FD values and the prediction; the same at k_e = 50 for the saddle weights), `check_saddle_clouds.py` (the three inadmissible
nearest-atom clouds: drift, the split along the unstable eigenvector, the function-level decomposition), `check_plan_chain.py` (the
½‖w − q‖₁ claim on the instance and on a collinear prior), `check_tests_vacuity.py`, `tests.out` (the three tests, 52 s),
`from_raw_*.out`. The two experiments (7 and 10 minutes) were not rerun; every quantity below was recomputed from the raw files or
from independent integrations.

## Verdicts

| object | verdict |
|---|---|
| **Corollary D‴.14 (i)** | **holds.** The particle law depends on i only through (ϑ_i, T_i); on 𝒞_m(π₀) the targets agree within a class (nearest-atom leak inside the cells; OT leak: with w = q and each class atom nearest to its own prior atom the diagonal plan is the unique optimal plan — the row-sum argument needs only the Voronoi condition, so the stated Δ_ν/2 bound is a sufficient condition), the vector field is tangent and equals (CL_w)'s, the plant sees Θ only through f̂ = Σ_kw_kΦ̃. Confirmed: the linearised field at a collapsed state leaks out of S_k by 2·10⁻¹⁶ and out of the collective subspace by 7·10⁻¹⁸; the collapsed 15-system and the 3-atom loop agree exactly (0.0) after 300 scheme steps. The non-invariance count "½‖w − q‖₁ … those particles leave" is a lower bound in general (D3e-10). |
| **Corollary D‴.14 (ii)** | **holds with fixes.** Theorem J's proof transfers to (CL_w) (its averaged drift is −(1/k_e)(diag(w)⁻¹⊗G)∇𝒢_w, similar to a symmetric matrix through (diag(w)⊗G⁻¹)^{1/2}, which is all Steps 4–6 use; "−G∇𝒢_w" is imprecise, D3e-13); the √m scaling of the G_m⁻¹-distance is right; the lift is Theorem J's orbit whenever it lies in the uniqueness ball, and "exactly, not up to O(1/k_e)" is justified (one K₀-atom orbit). But the "in particular for k_e ≥ √m·sup_t‖Θ^q_p − Θ*(q)‖_w/ρ₀(m)" is garbled by a factor k_e (D3e-2), (i)'s cell condition along Γ_q is a hidden hypothesis — W_adm does not contain it (D3e-3), and at the matching masses the ghost exponent is exactly σ₀/k_e (D3e-13). |
| **Corollary D‴.14 (iii)** | **holds with fixes.** D_orb(t) = Σ_kw^m_k‖ϑ^Γ_k(t; w^m) − ϑ^Γ_k(t; w^M)‖²_{G⁻¹} is exact under the shared-label convention of the coupling once both Theorem J orbits are the lifts; it needs Θ*(w^m), Θ*(w^M) hyperbolic for the identification and w^m, w^M ∈ W_adm (Corollary D‴.7's hypothesis) for the coupled distance to converge to it — unstated (D3e-3), and not academic: three of the six realised weight vectors of the numerics are inadmissible (D3e-1). |
| *Status* of D‴.14 | **holds with fixes.** The dependence on Theorem J's uniqueness ball is stated honestly, with the √m factor and the mean-norm basin as the open point; "(ii), (iii) proved … under the stated proviso" also needs the cell and W_adm hypotheses (D3e-3). |
| *Numerics* of D‴.14 (`ot_real_loop.md`, prose, §7 row, RESULTS) | **numbers hold, reading fails in part.** Every entry of tables 1–2 and every range of the Reading and the prose reproduces from the raw data (confirmations 2–3). Three of the six nearest-atom systems (seed 0 m = 15; seed 1 m = 15 and 60) have inadmissible realised weights: their lifts are saddle orbits and the m-systems sit on Proposition K's split of the over-weighted class, not on the lift — "every size on the orbit of its own weights" and "contracted clouds" are false for them, and the ratios 0.960 are the split (D3e-1, **major**). "m = 60 and 510 within 3·10⁻¹¹" is false for seed 1 (D3e-6). The OT-leak tracking error 3.9·10⁻⁵ and fit error 2.1·10⁻⁴ are the zero-order-hold floor of the scheme (∝ dt; the continuous orbit has e ≡ 0), so "8–50 times the OT leak's" compares against the floor (D3e-7). |
| **Proposition D‴.15** (statement, proof) | **holds with fixes.** The linearisation (δe = δx, the plant row −(1/m)Σ∂_ϑΦ̃δ_i, the particle blocks) is right; S_k is invariant with 𝕄|_{S_k} = M_k ⊗ I_{n_k−1}, the complementary subspace (plant coordinates and class means) is itself invariant — not only its quotient — and carries (CL_w)'s linearisation; confirmed at the field level (BᵀJB = I ⊗ A_k, off-pattern 0) and at the monodromy level (‖BᵀDP B − M_k‖ ≤ 2.1·10⁻⁷, collective block similar to the 3-atom monodromy through diag(1, √n_k)). The stability "iff" and the m-independence follow. The last claim (M_k's exponents = (σ₀I − R^G_k)/k_e + O(k_e⁻²)) is true, but its proof through Theorem J(b)'s multiset ("must therefore be") is incomplete and the Status ties it to the m-proviso unnecessarily: first-order averaging of the linear ghost equation proves it whenever Γ_K is Theorem J's orbit of (CL_w) (D3e-4). "targets constant along Γ" should be "locally constant near Γ" (D3e-3). |
| *Numerics* of D‴.15 (`lifted_orbit.md`, prose, §7 row, RESULTS) | **numbers hold with two wrong ranges.** Residuals, defects, exponents, multiplicities, predictions and the saddle's −9.13·10⁻⁵ ×8 all reproduce. "within 0.6–2.4%" is 0.56–3.22% and "within 0.3–3%" is 0.23–4.07% (D3e-5); the collective prediction is mislabelled by the weight whitening — the quoted numbers are λ_j(diag(q)^{−1/2}G^{1/2}∇²𝒢_qG^{1/2}diag(q)^{−1/2})/k_e, three times the stated matrix's (D3e-8); "to four digits" is 10⁻⁴ relative, and the md's "9.990e-05" is a rounding-order artefact (D3e-9). |
| `src/core/lifted_orbit.py` | **holds.** `rhs` equals `ClosedLoop.rhs` exactly for repeated targets (0.0), the weighted f̂ equals the replicated uniform one (4·10⁻¹⁷); `lift` exact; `class_separation_bases` orthonormal (≤ 7·10⁻¹⁶) with zero plant component, zero off-class entries and zero class sums, dimension (n_k − 1)p; `collective_basis` orthonormal, orthogonal to every S_k (≤ 2·10⁻¹⁶), total rank n + mp; `compression_spectrum`'s defect is zero iff span B is invariant (orthonormal B); `ghost_monodromy` co-integrates the ghost with the K-atom orbit (exact x_p, e_p), g_star is the class atom at the period start (return residual 0 to 6·10⁻¹⁷), the FD agrees with the matrix-ODE monodromy to 2.6·10⁻⁹; `monodromy_fd` is a plain central-difference Jacobian with a relative step. |
| `experiments/exp_ot_real_loop.py` | **holds with fixes.** Does what its docstring says (RK4 plant with the control held, Euler parameters, targets every 100 steps, window statistics on the 400-point orbit grid, float32 records 1.2·10⁻⁷ relative — harmless at the printed digits); regeneration byte-identical. Fixes: report the splitting margin of Θ*(w^m) and the within-class spread per system (D3e-1); label the scheme's tracking/fit floor or integrate the control inside the stages (D3e-7). |
| `experiments/exp_lifted_orbit.py` | **holds with cosmetic fixes.** `hess_weighted` is D‴.5(i)'s weighted whitened block (right; the md label is not, D3e-8); the `np.unique(np.round(ex, 7))` display loses a significant digit below 10⁻⁴ (D3e-9); regeneration byte-identical. |
| the three tests of `tests/test_ot_leak.py` | **pass (3/3 in 52 s) and test what their names say**, with one tautological assertion ("the atoms moved" compares particles with particle 0: 6.62 at t = 0) and one trivially true one (W² > 10⁻³ holds at t = 0 with W² = 3.37) in the two D‴.14 tests; test 3's tolerances have slack 40–10⁵ but bite at the 0.3 % level of the exponents (D3e-12). |
| §0 row, §7 rows, §8 bullet, RESULTS, README | **hold with fixes.** The §0 row and the Status state the proviso honestly; the §8 bullet's "Theorem J's orbit at a replicated configuration is the lift" is unconditional (D3e-11); the §7/RESULTS rows carry the wrong ranges (D3e-5), the "m = 60 and 510" clause lives in the prose only (D3e-6); §0's "within 10⁻⁸" is 7.5·10⁻⁹–2.7·10⁻⁸ in squared norm (D3e-13). |

## Confirmed issues

**D3e-1 (major; D‴.14 *Numerics* prose, `ot_real_loop.md` §2 title and Reading, §7 row, RESULTS).** Three of the six nearest-atom systems have
inadmissible realised weights, and the texts never mention admissibility — the central caveat of Proposition D‴.5(i),(iv), Corollary D‴.7
(w^m, w^M ∈ W_adm) and rounds 9–10. `splitting_margin_prior` at the stored weights: seed 0, m = 15, w = (6, 5, 4)/15: σ₀ − max κ_R = −0.00348
(κ_R = +0.01056 at atom 1); seed 1, m = 15, (6, 3, 6)/15: −0.00108 (+0.00816); seed 1, m = 60, (24, 15, 21)/60: −0.00216 (+0.00923); the
other three are admissible (+0.00561, +0.00440, +0.00544). For the three saddles the lifted orbit Γ_{w^m} is the unstable orbit of the
(9, 4, 3)/16 type, and the m-systems are not on it: at the start of the real loop (end of the averaged run) the over-weighted class has
whitened deviations from ϑ*_k(w^m) of −0.19 … +0.33 (seed 0, m = 15), −0.13 … +0.08 (seed 1, m = 15), −0.31 … +0.27 (seed 1, m = 60), with both
signs along the softest eigenvector of σ₀I − R^G_k (eigenvalue −0.0035, −0.0011, −0.0022; the centre direction, (0.004, 1.00, −0.001)) —
Proposition K's split, in progress: the averaged drift at the clouds is 2.0·10⁻³, 2.6·10⁻⁴, 7.9·10⁻⁴ against 1.5·10⁻⁵–5.3·10⁻⁵ for the three
admissible clouds and ≤ 2.5·10⁻¹⁵ at the lifts; the spread² of the split class grows 3.00·10⁻² → 3.07·10⁻² over the 4 slow units of the real
loop (e^{2·0.0035·4} = 1.03), the class means sit 0.024, 0.005, 0.023 (G⁻¹) from the lift's atom, against ≤ 1.6·10⁻⁴ and spread² ≤ 4.6·10⁻⁶
for the admissible systems. "Contracted clouds" is therefore false for these three (rms spread 0.17, 0.07, 0.16 whitened units in the
over-weighted class, growing), "every size on the orbit of its own weights" describes a saddle transient, and the agreement "within
3.3·10⁻⁶ of its own lift" is the split's second-order effect in f̂ (seed 0, m = 15 at t_end: ‖f̂_m − f̂_lift‖² = 3.4·10⁻⁶ while the
class-means configuration alone is 5.6·10⁻⁶ off). The two ratios 0.960 of "0.960–1.000" are exactly the two large splits (seed 0 m = 15:
0.9600; seed 1 m = 60: 0.9604); the small split gives 0.9976 and the admissible pair 1.0003. Corollary D‴.14(iii)'s coupled distance does not
converge to D_orb for these sizes (the orbit is a saddle), so the numerics support (iii) in one pair of four, not four. *Fix:* report the
splitting margin of Θ*(w^m) and the within-class spread per system (as `exp_chaos_learning.py` and `exp_prior_design.py` already do), say
which systems sit on saddles and that their lifts are the unstable orbits whose distance to the m-system is the split, keep the
function-level comparison with that reading, and state w^m, w^M ∈ W_adm in (iii) (hyperbolicity suffices for the identification,
admissibility is needed for the basin and the coupled distance).

**D3e-2 (minor; (ii) statement).** "in particular for k_e ≥ max{k₀(m), √m·sup_t‖Θ^q_p(t) − Θ*(q)‖_w/ρ₀(m)}" is dimensionally garbled: the
uniqueness-ball condition is √m·sup_t‖Θ^q_p(t) − Θ*(q)‖_w ≤ ρ₀(m), which carries no k_e; with Theorem J(a) for (CL_q), sup_t‖·‖_w ≤ C_q/k_e, it
becomes k_e ≥ √m·C_q/ρ₀(m). As written (sup‖·‖_w in place of C_q = k_e·sup‖·‖_w) the condition is weaker by the factor k_e = 25 on the
instance. ‖·‖_w is used without definition (Σ_kw_k‖·‖²_{G⁻¹}). *Fix:* write the condition on the distance, then the k_e form with C_q.

**D3e-3 (minor; (ii), (iii), Status; D‴.15 hypothesis).** Hidden hypotheses. (a) (ii) needs (i)'s cell condition along Γ_q — sup_t max_k
‖ϑ^Γ_k(t) − ϑ⁰_k‖_{G⁻¹} < Δ_ν/2, i.e. ε(q) + C_q/k_e < Δ_ν/2 — for the lift to be an orbit of (CL) with the OT leak; q ∈ W_adm does not imply it
(D‴.12(iii) states "if ε < Δ_ν/2" explicitly). (b) (iii) needs Θ*(w^m), Θ*(w^M) hyperbolic (Theorem J) and both weights admissible for the
coupled distance to converge to D_orb (Corollary D‴.7's hypothesis), and the Voronoi condition along Γ_{w^m}; on the instance three of six
weights are inadmissible (D3e-1). (c) Proposition D‴.15's "any leak whose targets are constant along Γ within each class" should read
"locally constant in Θ in a neighbourhood of Γ" (a unique vertex plan; atoms off the Voronoi boundaries): that is what makes the Poincaré map
differentiable with the stated linearisation, and what `ClassLeakLoop`'s frozen targets assume. On the instance all three hold (ε = 0 at
the matching masses; the Θ*(w^m) atoms are within 1.1 of the prior atoms against Δ_ν/2 = 3.13). *Fix:* state them.

**D3e-4 (minor; D‴.15 proof, last sentence, and Status).** "Theorem J(b) at Θ*_m(w) … whose separation eigenvalues are those of σ₀I − R^G_k with
multiplicity n_k − 1 (Proposition D‴.5(i)) and must therefore be the exponents of M_k" is not a proof: Theorem J(b) gives the *multiset* of
the mp parameter exponents up to O(k_e⁻²), and attributing the (n_k − 1)-fold values to class k needs an argument (at the matching masses
all classes have the same separation exponent, so multiplicities cannot attribute; in general classes' exponents may coincide). The
direct route proves the claim and removes the m-proviso: along Theorem J's orbit of (CL_w), e_p = εr_{Θ_p}(x_m) + O(ε²) (Step 1), so the
ghost's coefficient is A_k(t) = ε[G∇²_ϑΦ̃(x_m(t);ϑ_{p,k}(t))r_{Θ_p(t)}(x_m(t)) − σ₀I] + O(ε²), with period mean ε(GR[μ_{Θ_p}](ϑ_{p,k}) − σ₀I) +
O(ε²) = ε(GR_k − σ₀I) + O(ε²) (Θ_p = Θ* + O(ε)); hence M_k = I + εT₀(GR_k − σ₀I) + O(ε²), diagonalisable (GR_k ∼ R^G_k), and Bauer–Fike as in
Step 5 gives the exponents (1/k_e)λ_j(σ₀I − R^G_k) + O(k_e⁻²) whenever Γ_K is Theorem J's orbit of (CL_w) — no uniqueness ball in m is
needed, so the Status' "holds when the lift is Theorem J's orbit (Corollary D‴.14(ii)'s proviso) and is formal otherwise" is
over-conditional. Confirmed numerically: the FD ghost monodromy equals the matrix-ODE monodromy Ṃ = A_k(t)M along the stored orbits to
2.6·10⁻⁹ in the exponents, and the eigenvalues of the period-averaged −Ā_k equal the ODE's to ≤ 7.3·10⁻⁹ (the commutator term is
negligible), while the prediction (σ₀I − R^G_k)/k_e differs by 7·10⁻⁶–6·10⁻⁵ — the whole O(k_e⁻²) error is R evaluated at Θ*(w) instead of
along Θ_p. At k_e = 50 (shooting residual 2.7·10⁻¹¹, period-start offset 5.6·10⁻³ against 1.6·10⁻² at k_e = 25) the errors fall by 5.7–12.4×
(≥ the k_e⁻² factor 4): the saddle's negative exponent −4.422·10⁻⁵ against the prediction −4.384·10⁻⁵, relative error 4.1 % → 0.86 %. *Fix:*
replace the last proof sentence by the averaging argument and drop the proviso from the Status.

**D3e-5 (minor; D‴.15 *Numerics* prose ×2, §7 row, RESULTS, commit message).** "the collective exponents … within 0.6–2.4% of Theorem J(b)'s
prediction": the relative errors |ex − pred|/pred are 0.56, 1.59, 1.89, 1.76, 1.60, 0.57, **3.22**, 2.27, 1.95 % (the seventh exponent 3.704·10⁻³
against 3.828·10⁻³) — range 0.6–3.2%. "the prediction (σ₀ − λ_j(R^G_k))/k_e within 0.3–3%": over the 18 Part-B exponents the relative errors
are 0.23 … 4.07 %, the largest at the saddle's negative exponent (−9.126·10⁻⁵ against −8.769·10⁻⁵ — quoted in the same sentence as
"−9.13·10⁻⁵ (predicted −8.77·10⁻⁵)", which is 4.1 %), the smallest at class 1 of (7, 5, 4)/16 (8.219·10⁻⁴ against 8.238·10⁻⁴). *Fix:* 0.6–3.2 %
and 0.2–4.1 %, with Theorem J's own remark that an O(k_e⁻²) absolute error (here 2·10⁻⁶–7·10⁻⁵) is a large fraction of a small exponent.

**D3e-6 (minor; D‴.14 *Numerics* prose).** "(each m-system within 3.3·10⁻⁶ of its own lift, m = 60 and 510 within 3·10⁻¹¹)": seed 1's m = 60 is
1.84·10⁻⁶ from its lift (a saddle, D3e-1); only seed 0's m = 60 (1.47·10⁻¹¹) and both m = 510 (2.52·10⁻¹¹, 9.61·10⁻¹²) are at 3·10⁻¹¹. *Fix:*
"the three admissible systems within 3·10⁻¹¹".

**D3e-7 (minor; D‴.14 *Numerics* prose, `ot_real_loop.md` Reading, §0 row, §7 row, RESULTS).** The OT-leak "rms tracking error 3.9·10⁻⁵–4.0·10⁻⁵
for every size" and the lift's fit error 2.1·10⁻⁴ are the zero-order-hold floor of the integration scheme, not a property of the loop. The
continuous (CL) at the exact representation has e ≡ 0 (DOP853 from (x_m(0), exact): max|e| = 2.2·10⁻¹³ over a period, max|Θ − exact| =
6·10⁻¹⁵; `exp_lifted_orbit`'s shooting residual 3.4·10⁻¹⁶); the scheme of `real_run` (RK4 with u = ẋ_m − k_ee − f̂ held over the step) gives
rms e = 7.83·10⁻⁵, 3.91·10⁻⁵, 1.96·10⁻⁵ and fit error 4.22·10⁻⁴, 2.11·10⁻⁴, 1.06·10⁻⁴ at dt = 2·10⁻³, 10⁻³, 5·10⁻⁴ — exactly ∝ dt — and
re-evaluating the control in every RK4 stage gives rms e = 8·10⁻¹². The scheme is first-order accurate in the control despite the RK4
plant. So "the weight bias of the fit is a tracking error, 8–50 times the OT leak's" compares the nearest-atom errors against the scheme's
floor (the true OT-leak tracking error at the matching masses is zero), and the §0/§7/RESULTS "rms tracking error 4·10⁻⁵ for every size"
is a discretisation number. The nearest-atom values are unaffected (seed 0 m = 15's lift: scheme 1.193·10⁻³ against the true orbit's
1.182·10⁻³). *Fix:* "at the scheme's floor 3.9·10⁻⁵ (∝ dt; the continuous orbit at the matching masses has e ≡ 0)", or integrate the control
inside the stages.

**D3e-8 (cosmetic; D‴.15 *Numerics*, `lifted_orbit.md`, `exp_lifted_orbit.py` docstring/report).** Theorem J(b)'s collective prediction is
labelled "λ_j(G^{1/2}∇²𝒢_qG^{1/2})/k_e"; the stored and quoted numbers (2.952·10⁻⁴, …, 8.091·10⁻³) are the eigenvalues of Proposition
D‴.5(i)'s weighted whitened block diag(q)^{−1/2}G^{1/2}∇²𝒢_qG^{1/2}diag(q)^{−1/2}/k_e (= G_3^{1/2}∇²J_3G_3^{1/2}/k_e at q = 1/3), which `hess_weighted`
computes correctly; the matrix as labelled has eigenvalues three times smaller (9.84·10⁻⁵, …, 2.70·10⁻³). *Fix:* the label.

**D3e-9 (cosmetic; D‴.15 *Numerics*, §7 row, RESULTS, `lifted_orbit.md`).** "the same for m = 3, 6, 9 and 15 to four digits": the m = 15
collective exponents differ from m = 3's in the fourth digit for two of nine (3.9984·10⁻⁴ vs 3.9986·10⁻⁴, 8.7329·10⁻⁴ vs 8.7324·10⁻⁴; m = 6 for
one) — 8·10⁻⁵ relative. The md's lift separation exponent "9.990e-05" against the ghost's "9.989e-05" is a display artefact:
`np.unique(np.round(ex, 7))` rounds to 10⁻⁷ before the `.3e` format, losing a digit below 10⁻⁴ (the actual values 9.98934·10⁻⁵ and
9.98916·10⁻⁵ agree to 4·10⁻⁹). *Fix:* "to 10⁻⁴ relative"; format before rounding.

**D3e-10 (cosmetic; (i) statement).** "the plan assigns ½‖w − q‖₁ of the mass of the over-weighted classes to other atoms, and those particles
leave their class": the reassigned mass is at least ½‖w − q‖₁; with the squared cost the optimal vertex plan can route the excess through
an intermediate class. Collinear prior at 0, 1, 2 (G = I), q uniform, counts (3, 2, 1): the chain 0→1, 1→2 costs 2 against the direct 4, the
plan moves two particles for ½‖w − q‖₁·m = 1. On the instance the direct route is cheaper (44.7 against 91.1 through the middle atom) and
exactly one particle moves for (6, 5, 4)/15, as the test asserts. *Fix:* "at least ½‖w − q‖₁".

**D3e-11 (cosmetic; §8 bullet).** "Theorem J's orbit at a replicated configuration is the lift of one K₀-atom orbit" is stated unconditionally,
where the corollary has it under the uniqueness-ball proviso — whose √m factor is the point: the ball is in the sum norm, so a within-class
spread δ per particle is √m·δ from Θ*_m and both the identification and Theorem J(c)'s basin tighten with m unless ρ₀(m) grows like √m (the
"basin in the mean norm" the bullet lists as open). *Fix:* "under Corollary D‴.14(ii)'s proviso".

**D3e-12 (cosmetic; tests).** Test 1's "the atoms moved (a genuine orbit, not a fixed point)" asserts `np.abs(Th15 − Th15[0]).max() > 1e-3`, the
spread between particle 0 and the others — 6.62 at t = 0 by construction — not motion (the actual motion over the 300 steps is 2.3·10⁻³;
the collapsed agreement is exactly 0.0 against the 10⁻¹² tolerance, the scheme being permutation-exact). Test 2's final `W2 > 1e-3` holds at
t = 0 (W² = 3.37 between the (6, 5, 4)/15 and (5, 5, 5)/15 lifts of the same atoms). The docstrings cite "Corollary D‴.14(a)/(b)" for (i)/(ii).
Test 3's slack: defects 10⁻⁸ against 10⁻⁶, |ex − exg| 2.5·10⁻⁸ against 10⁻⁶, |ex − σ₀/k_e| 2.5·10⁻⁸ against 10⁻⁵, residual 1.7·10⁻¹⁵ against
10⁻¹⁰ — loose, but 10⁻⁶ on exponents of 2.8·10⁻⁴ is a 0.3 % test. *Fix:* compare with the initial configuration; drop or sharpen the W² line.

**D3e-13 (cosmetic; (ii) wording, §0 row).** At the matching masses the ghost exponent is exactly σ₀/k_e, not "σ₀/k_e + O(k_e⁻²)" (e_p ≡ 0
makes A_k ≡ −σI; ‖M_k − e^{−σT₀}I‖ ≤ 8·10⁻⁸, the FD floor); "the averaged gradient structure −G∇𝒢_w" is −(G/w_k)∇_{ϑ_k}𝒢_w, the gradient flow in
the metric diag(w)⊗G⁻¹ — the similarity to a symmetric matrix that Steps 5–6 need; §0's "within 10⁻⁸ of each other" is 7.5·10⁻⁹–2.7·10⁻⁸ in
squared L²(orbit) norm.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_ot_real_loop.py --from-raw` and `exp_lifted_orbit.py --from-raw` (on scratch copies) rewrite both md files
   byte-identically (md5 `ff893214…`, `50c437b6…`); 3 s and 1 s.
2. **`ot_real_loop.md` table 1** (all 8 rows × 6 columns) and **table 2** (6 rows × 9 columns) from the raw records (1000 records at 0.1 time
   units, 314 in the window t > 68.58; float32 f̂ records). Prose ranges: OT ‖f̂_m − f̂_M‖² 7.452·10⁻⁹–2.743·10⁻⁸; to the lift 5.198·10⁻¹⁰–3.557·10⁻⁸;
   fit 2.111·10⁻⁴ (lift) / 2.124·10⁻⁴–3.748·10⁻⁴; rms e 3.927·10⁻⁵–4.002·10⁻⁵ (lift 3.912·10⁻⁵); W²(μ_m, μ_M) 3.346·10⁻⁴–1.390·10⁻³; nearest-atom
   ‖f̂_m − f̂_M‖² 1.064·10⁻⁴–2.467·10⁻³, ratios to the lifts' distance 0.9600, 1.0003, 0.9976, 0.9604, each m-system to its lift ≤ 3.32·10⁻⁶, fit
   0.0125–0.0747, rms e 3.033·10⁻⁴–1.914·10⁻³ (7.6–48.7 times the OT values), W² 1.569–5.547; horizon 15.92 periods = 4.0 slow units.
3. **The clouds.** Nearest-atom weights conserved from the averaged run to the end of the real loop (all six). OT clouds carry q = 1/3; their
   within-class spread² is 10⁻⁷ in the classes that kept their particles and 2·10⁻⁴–3.4·10⁻³ in the classes that received migrants (rms 0.014–
   0.058 whitened units — the transit's arrivals, still contracting at σ₀), the class means 2–5 × spread² off the lift's atoms (e.g. 0.0093
   against 3.0·10⁻³), and ‖f̂_m − f̂_lift‖² at t_end (3.5·10⁻⁸, 5.1·10⁻¹⁰, 1.2·10⁻⁸, 8.5·10⁻⁹, 4.9·10⁻⁹, 4.2·10⁻⁹) equals the window values within 2 %:
   "the within-class spread, at second order" holds.
4. **`lifted_orbit.md`**: 3-atom residual 3.4·10⁻¹⁶ (zK = (x_m(0), exact) exactly); parameter exponents 2.9355·10⁻⁴ … 7.9337·10⁻³; ghosts 2.82975·10⁻⁴
   within 8.5·10⁻⁹ of σ₀/k_e = 2.829754·10⁻⁴; lifts m = 6, 9, 15: residuals 1.7·10⁻¹⁵, 1.3·10⁻¹⁵, 7.1·10⁻¹⁶, separation exponents all within
   3.1·10⁻⁸ of σ₀/k_e with the right counts (3, 6, 12 per class), defects 9.1·10⁻¹⁷–1.8·10⁻⁸, collective defects 2.1·10⁻⁹–2.8·10⁻⁹, collective
   exponents within 8.7·10⁻⁸ of m = 3's. Part B: margins +0.00257 / −0.00219, κ_R (+0.00450, −0.00298, −0.00108) and (+0.00927, −0.00555,
   −0.00572), shooting residuals 1.1·10⁻¹⁴ / 7.1·10⁻¹⁵, offsets 1.12·10⁻² / 1.60·10⁻², lift residuals 4.9·10⁻¹⁴ / 4.5·10⁻¹⁴, lift separation
   exponents equal to the ghosts' within 4.8·10⁻⁹ with multiplicities 6, 4, 3 and 8, 3, 2, defects 4.0·10⁻⁸–5.7·10⁻⁸, collective exponents
   equal to the K-orbit's within 4.8·10⁻⁹; exactly 8 negative exponents −9.1256·10⁻⁵ … −9.1259·10⁻⁵ in the saddle lift's spectrum.
5. **Code.** `ClassLeakLoop.rhs` − `ClosedLoop.rhs` = 0.0 (repeated targets, weights None or 1/m); weighted f̂ against the replicated uniform
   f̂ 4·10⁻¹⁷; `lift` exact; bases: orthonormality ≤ 6.7·10⁻¹⁶, plant components and off-class entries 0, class sums ≤ 4.4·10⁻¹⁶, Bcᵀ B_k ≤
   2.2·10⁻¹⁶, dimensions (n_k − 1)p and n + K₀p, total rank n + mp (for (2,2,2), (7,5,4), (9,4,3), (1,2,3) — a singleton class has no S_k).
6. **The block structure at the field level** (FD Jacobian of `rhs` at a collapsed (3, 2, 2) state off the orbit): relative leakage out of S_k
   2.0·10⁻¹⁶–4.7·10⁻¹⁶, out of the collective subspace 7.1·10⁻¹⁸, plant row on S_k 1.9·10⁻¹⁷, BkᵀJBk = I_{n_k−1} ⊗ A_k with off-pattern norm 0
   (class 0: 1.5·10⁻¹⁹) — Proposition D‴.15's algebra verbatim.
7. **The m = 6 monodromy beyond the exponents**: ‖BᵀDP B − M_k‖ = 6.9·10⁻⁸, 2.1·10⁻⁷, 1.9·10⁻⁷ (‖M_k‖ = 1.73) and ‖M_k − e^{−σT₀}I‖ ≤ 7.9·10⁻⁸;
   collective compression equals diag(1, √n_k)·DPK·diag(1, √n_k)⁻¹ to 1.8·10⁻⁶ (dominated by the plant column, where the FD step meets the
   fast decay), eigenvalues to 4.2·10⁻⁸. FD floor: the smallest multiplier of DP6 is 2.3·10⁻⁹ against the true error-channel multiplier
   e^{−k_eT₀} ≈ 6·10⁻⁶⁹ — the error channel is invisible, as the md says; the defects of 10⁻⁸ are the central-difference roundoff floor
   (ε_mach/fd ≈ 10⁻⁹ per entry), and the 9·10⁻¹⁷ defect of class 0 at m = 6 is the floating-point commutativity of the first two summands
   of f̂, not a better invariance.
8. **The ghost by three routes** (D3e-4): FD = matrix ODE to 2.6·10⁻⁹, matrix ODE = e^{T₀Ā_k} to 7.3·10⁻⁹, prediction off by 7.5·10⁻⁶–6.3·10⁻⁵;
   k_e = 50: error ratios 9.4, 12.4, 11.9 (class 0), 6.8, 5.7, 7.0 (class 1), 7.2, 6.1, 5.7 (class 2) — O(k_e⁻²) or better.
9. **The scheme** (D3e-7): e_rms and fit error ∝ dt at the exact representation; continuous e ≡ 0 to 2·10⁻¹³; no-hold scheme 8·10⁻¹². The
   nearest-atom lift of seed 0, m = 15 through the scheme reproduces the table (1.193·10⁻³, 4.549·10⁻²); the true orbit Γ_w (shooting
   residual 1.2·10⁻¹³, offset 3.1·10⁻³ from Θ*(w)) has rms e 1.182·10⁻³ and fit 4.596·10⁻², and the scheme's state over the window is
   1.5·10⁻⁷ from it in ‖f̂‖²_{L²(orbit)} and 2.9·10⁻³ in ‖·‖_w — the slowest collective exponent 2.93·10⁻⁴ per time unit (relaxation 3400 time
   units) means no system has relaxed to its orbit within the 100-unit horizon, but the function-level comparisons are unaffected (10⁻⁷).
10. **The saddles** (D3e-1): margins, κ_R, unstable eigenvectors, the two-signed projections, the drifts, the growth of the spread, and the
    decomposition ‖f̂_m − f̂_lift‖² = 3.4·10⁻⁶ / ‖f̂_m − f̂_means‖² = 1.4·10⁻⁵ / ‖f̂_means − f̂_lift‖² = 5.6·10⁻⁶ (seed 0, m = 15); the admissible
    systems' class deviations ≤ 6·10⁻³ (max) with projections ≤ 10⁻³ on the softest direction.
11. **The plan** (D3e-10): on the instance (6,5,4), (6,4,5), (5,6,4)/15 move exactly one particle (0→2, 0→1, 1→2), (7,4,4)/15 two; the pairwise
    prior costs 51.9, 44.7, 39.2; the collinear counterexample moves two.
12. **Algebra.** ‖Θ − Θ*_m‖²_{G_m⁻¹} = Σ_kn_k‖ϑ_k − ϑ*_k‖² = m‖·‖²_w (the √m); the lift at the matching masses sits at Θ*_m(q) exactly, so it is
    Theorem J's orbit for every commensurate m as soon as k_e ≥ k₀(m); the stability "iff" (the lift's multipliers are Γ_K's and M_k's);
    exp(−σT₀) exactly for the ghost at zero residual.
13. **Tests**: 3 passed in 52.4 s; `git status` clean.

## Not checked

The two experiments were not rerun from scratch (deterministic; everything above was recomputed from the raw files or re-integrated for
the systems concerned). Theorem J's proof, Proposition D‴.12(iii)'s lemma and D‴.5(i)'s block splitting were used as given (rounds 5, 9, 11).
The O(k_e⁻²) scaling of D3e-4 was checked at k_e = 25 and 50 for the saddle weights only. Whether the three saddle m-systems of D3e-1
eventually split locally or migrate (round 10 saw both) was not followed beyond the stored horizon. The basin radius ρ₀(m) in the mean
norm is open by the authors' own account and was not examined. Incommensurate sizes in the real loop are not part of the material.
