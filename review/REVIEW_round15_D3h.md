# Referee report D3h (round 15) — PARADIGM §5d: Proposition D‴.20 (the entropic OT leak) with its Status and *Numerics*; the §0 row, the §7 row and the §8 clause; `results/RESULTS_PARADIGM.md` (D‴.20); `src/core/entropic_leak.py`; `experiments/exp_entropic_leak.py` → `results/chaos/entropic_leak.md` (+ `entropic_leak_raw.npz`); `tests/test_entropic_leak.py`; [S45], [S46] and the D‴.20 audit bullet of `LITERATURE_TRANSPORT.md`

Branch `claude/festive-bardeen-eopej2`, reviewed at 80ab38d (HEAD when the review started). While it ran, two commits from other sessions
landed: ebd4c03 (referee round 14 applied, Theorem D‴.18) and c29c34a (Theorem D‴.21, inserted after D‴.20's *Numerics*, with its §0/§7 rows,
§8 text, RESULTS section and audit bullet). Neither changes a line of the reviewed material — the D‴.20 block of PARADIGM.md, its §0 and §7 rows
and the §8 clause are byte-identical at c29c34a, only shifted; the module, script, raw file, md and tests are untouched — except that c29c34a
inserted the D‴.21 audit bullet between the D‴.20 bullet and its Verdict line (D3h-11). Every finding applies to c29c34a, and **line numbers
below are those of c29c34a** (PARADIGM.md: D‴.20 at 2374–2440, i.e. 80ab38d's + 41; the §0 row at 114, the §7 row at 4988, the §8 clause at
5210; RESULTS_PARADIGM.md 467–478; LITERATURE_TRANSPORT.md [S45], [S46] at 66–67, the D‴.20 bullet at 323–327, its Verdict at 338).
No repository file was modified: the `--from-raw` regeneration ran on a copy of the script in a scratch tree (`regen/`, `src` symlinked); every
other script imports the repository read-only (`PYTHONDONTWRITEBYTECODE=1`; pytest with `-p no:cacheprovider`); single runs were reproduced
from their seeds with the script's own scheme. Scratch scripts and their `.out` files in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3h/`:
`indep.py` (an independent implementation of J^ε_m — features, quadrature, entropic OT by Newton on the gauge-fixed semi-dual, OT_ε as primal
value, the K-atom energy, R and Cov^G — no repository code), `check_i_ii_iii.py` ((i)–(iii) by finite differences of `indep.J`, m = 20,
N = 256), `misc_checks.py` (collective block, semi-dual kernel, replicated OT_ε), `crit_repo_check.py` (the critical-excess routine of the
script re-run with the Newton residual logged), `crit_continuation.py` (the same with continuation in ε), `crit_indep.py` (the table by
`indep.py`), `fd_delta02.py` (δ = 0.2 by the finite-difference Hessian of `indep.J` at m = 30), `critical_excess_eps.py` (δ_c(ε), p_c(ε), the
w = q statement), `check_A_posdef.py` (σ₀I − R^G_k along ε at six weights), `inspect_raw.py`, `recompute_runs.py` (every dynamics number per
seed), `seed_margins.py` (the margins at every run's start and end weights), `final_drift.py` (the stored finals' drift and spread),
`long_run.py` (+ `long_run_{1_4_240,3_16_60,0_16_15,3_64_15}.out/.npz`: four runs reproduced from their seeds and continued),
`endstate_margins.py`, `bias_fit.py`, `solver_audit.py` (+ `solver_audit_2_1_510.out`), `floors_compare.py`, `floor_law.py`,
`incomm_floor.py`, `check_tests.py`, `test4_A.py`, `tests.out` (4 passed, 1.25 s), `from_raw.out`, `regen/`.

## Verdicts

| object | verdict |
|---|---|
| **D‴.20 (i)** | **holds** (cosmetic wording, D3h-9). (𝔄_m^ε) is exactly the G-gradient flow of J^ε_m: the finite-difference gradient of an independent J^ε_m equals ∇V_L + σ₀G⁻¹(ϑ_i − T^ε_i) to 2–3·10⁻¹¹ (m = 7, 12; ε = 0.5, 2, 30), and ∂OT_ε/∂ϑ_i = (2/m)G⁻¹(ϑ_i − T^ε_i) to 5·10⁻¹⁰ (the factor 2 from ∂C_{ij}/∂ϑ_i = 2G⁻¹(ϑ_i − ϑ⁰_j), 1/m from the row mass, σ₀m/2 · 2/m = σ₀). The KL reference does not depend on Θ: on labels it is a⊗b with a = 1/m; read on the measure, the unique plan is symmetric among coincident particles and OT_ε of a replicated configuration equals the K-atom value to 14 digits. The semi-dual is strictly concave modulo constants (Hessian eigenvalues −0.0099, −0.0067, 0 with kernel 1), so the IFT applies on the gauge g₁ = 0; L̄ is real-analytic (tanh/logistic saturations, Gaussians with s ≥ s_min > 0, a finite quadrature or a period average of an analytic integrand); Łojasiewicz applies, and boundedness is in fact automatic (T^ε_i ∈ conv supp ν, bounded fit force). Limits: \|T^ε − T̄\| = 1.26, 0.78, 0.27, 7·10⁻⁴ at ε = 1, 0.3, 0.1, 0.03 (unique LP plan), \|T^ε − Σq_jϑ⁰_j\| = 0.70, 0.067, 0.0066 at ε = 10², 10³, 10⁴. |
| **D‴.20 (ii)** | **holds.** Class symmetry π^ε_{ij} = π^K_{kj}/n_k to 2·10⁻¹⁷, the class targets to 2·10⁻¹⁵, the replicated weighted stationary point critical for J^ε_m to 2·10⁻¹⁶, J^ε_m = m𝒢^ε_w to 12 digits (m = 20, counts (8, 7, 5), ε = 12 and 30; the energy identity at ε = 12). |
| **D‴.20 (iii)** | **holds** (derived independently: dT_i = (1/ε)Cov_{P_i}(ϑ⁰, dg − dC_{i·}); a zero-sum class perturbation leaves the column sums stationary, so dg = 0 and ∂T_i/∂ϑ_i = (2/ε)Cov G⁻¹; the Gauss–Newton term of mL̄ cancels on zero-sum perturbations). At m = 20, N = 256 (not the test's size or quadrature), ε = 12 (class 1 unstable) and 30, the finite-difference Hessian compressed to each class's separation subspace equals I_{n_k−1} ⊗ S^ε_k **entrywise** to 1.4·10⁻¹², the subspace is invariant (coupling to its complement ≤ 2.9·10⁻¹²: decoupled from the collective block and from the other classes), and the collective block equals the whitened weighted K-atom Hessian to 7.9·10⁻¹² (D‴.5(i)'s form). R^G_k is D‴.5's object, the residual-weighted curvature at the configuration's own (non-zero) residual. |
| **D‴.20 (iv)** | **the rank-one algebra holds** (S = A − (2σ₀/ε)p(1−p)Δ̂Δ̂ᵀ, A ≻ 0 ⇒ negative eigenvalue iff (2σ₀/ε)p(1−p)Δ̂ᵀA⁻¹Δ̂ > 1); the w = q sentence holds (entropic term 10⁻¹⁶ and 3·10⁻⁴–1·10⁻³ of σ₀ at ε = 1, 4). **"since σ₀I − R^G_k ≻ 0" is a hidden hypothesis** — it fails at the test's own weight (6, 5, 4)/15 and at (8, 4, 3)/15 for ε ≤ 4 (D3h-5); ε_c is implicit (p, A depend on ε); "stable" needs the collective block. **"it splits, and the split-off part carries the excess to atom j"** is a reading, contradicted on the instance (D3h-4). |
| **D‴.20 (v)** | **formal, holds as a scaling** (p_c ≈ ε/(2‖Δ̂‖²) confirmed for the conditional fraction within 8–26% at ε ≤ 8; δ_c(ε) > 0 at every ε tried), with false or loose wording: "only the hard constraint … for every m", "instead of conserved", p used for two quantities (D3h-6, D3h-7). |
| Status | **holds** for (i)–(iii) and (v) ("formal"); (iv) "the rank-one criterion proved, its two-atom form an approximation" is fair, but the two-atom form is never realised in the numerics (D3h-3). The deterministic-annealing gloss is inaccurate (D3h-8). |
| *Numerics*, the critical excess | **δ ≤ 0.1 hold** (reproduced independently to 4–5 digits: ε_c = 5.1066, 13.5077, 23.1045, p, the rank-one values, every grid margin). **The δ = 0.2 row fails**: ε_c = 136, p = 0.62, the rank-one 88 and three of the five grid margins come from unconverged Newton solves; the true values are ε_c = 65.0, p = 0.574, rank-one 81.9 (D3h-1, **major**). **"charges mainly two atoms" is false** at every δ (D3h-3). |
| *Numerics*, dynamics | **numbers hold** (the md regenerates byte-identically; every per-seed number recomputed; four runs reproduced bit for bit at s = 600). **Readings fail**: the ε = 4 "residual excesses … the critical excess at that ε" are 600-unit snapshots of unstable systems (D3h-2, **major**); "keep their sampled weights" (ε = 16), "no excess is equalised" (ε = 64) and the split-equalises gloss fail on continued runs (D3h-4); "the floor is the nearest-atom leak's" (§0) is 0.23–0.84 of it (D3h-6). |
| `src/core/entropic_leak.py` | **holds as a solver** (along four reproduced and continued runs, to s = 3000–6000, the plan's column marginals are exact to ≤ 1·10⁻¹³ at every refresh; 3006 solves of the hardest run, ε = 1, m = 510: ≤ 7 Newton iterations, no `max_iter`, no `lstsq` fallback). **`weighted_stationary_atoms_entropic` fails silently** (returns a non-converged point; the caller must check) — the cause of D3h-1. Cosmetic: the docstring's "log-domain Sinkhorn plan". |
| `experiments/exp_entropic_leak.py` | **the integrator holds** (midpoint rule, ds = 0.02, the target lagged ≤ 0.2 slow units; the clouds are D‴.12's — the start weights coincide in all 12 (seed, m) pairs); **`critical_row` ignores the Newton residual** and restarts Newton from the exact representation at every ε (D3h-1); the report prints max-over-seeds columns that mix seeds and no margins (D3h-2); its header names the wrong solver (D3h-9). |
| `tests/test_entropic_leak.py` | **pass (4 in 1.25 s) and test what they say**; no assertion is vacuous (test 3's control σ₀I − R^G differs from the FD spectrum by 1.3·10⁻⁴–7.6·10⁻³ against the tolerance 10⁻⁶), but slacks are 10³–10¹⁷, test 3 compares spectra only and test 4 ignores the Newton residual (D3h-10). The (6, 5, 4)/15 threshold, ε_c ∈ (17, 18) at n_quad = 96 and 17.87 at N = 256, is consistent with the table. |
| [S45], [S46], audit bullet | **adequate as citations** (Cuturi 2013; Peyré–Cuturi 2019; Cominetti–San Martín 1994 for ε → 0 — the limit is the maximal-entropy optimal plan, so "wherever the OT plan is unique" is the right proviso; Courty et al. 2017 for barycentric maps; Rose–Gurewitz–Fox 1990 and Rose 1998 for T_c = 2λ_max(C_{x\|y})); "with the roles exchanged" is not an accurate description of the analogy, and the ε → 0/ε → ∞ interpolation has a precedent not cited (D3h-8); at c29c34a the bullet has lost its Verdict line to the D‴.21 bullet (D3h-11). |
| §0 row (114), §7 row (4988), §8 clause (5210), RESULTS (467–478) | carry D3h-1 (136, 88), D3h-2 ("stop at the residual excess(es)"), D3h-4 ("keep their sampled weights", "and the masses equalise"), D3h-6 ("the floor is the nearest-atom leak's", "floor-free for every m", "equalises … down to a critical excess"). |

## Confirmed issues

**D3h-1 (major; PARADIGM lines 2429–2431, §7 row line 4988, RESULTS line 470, `entropic_leak.md` line 12, `exp_entropic_leak.py`
`critical_row` lines 72–105, `entropic_leak.py` `weighted_stationary_atoms_entropic` lines 176–197; the 80ab38d commit message "ε_c = 5.1–136").**
The δ = 0.2 row of the critical-excess table — "ε_c = … 136 for δ = … 0.2; at ε_c the class sends the fraction p = … 0.62 …, and the rank-one
formula gives … 88 — … off at δ = 0.2, where it charges all three" — is an artefact of unconverged Newton solves.
(a) *The raw file says so.* The `grid` entries of `entropic_leak_raw.npz` store the Newton residual: for δ = 0.2 it is 4.1·10⁻², 3.7·10⁻²,
4.9·10⁻² at ε = 1, 16, 64 (≤ 7.4·10⁻¹³ for every other entry of the table), so the printed margins −0.17299, −0.04808, −0.36601 are evaluated off
the stationary point (`inspect_raw.py`). `margin(eps)` restarts Newton from the exact representation at every ε, and
`weighted_stationary_atoms_entropic` returns whatever it reached (its line search accepts the last trial once λ < 10⁻⁸); `critical_row`
never looks at `res`.
(b) *The bisection found Newton's convergence boundary, not a sign change.* Re-running the routine with the residual logged
(`crit_repo_check.py`): 15 of the 30 bisection evaluations are unconverged (residuals up to 5.3·10⁻², margins ≈ −0.36), the converged ones in
(136.6, 181) are all positive (+0.00090 to +0.00108), and at the reported ε_c = 136.14 the converged margin is **+8.96·10⁻⁴**; there
2σ₀λ_max(A^{−1/2}Cov^G A^{−1/2}) = 69.3, not 136.
(c) *The true values.* Continuing Θ*_ε(w_0.2) in ε from 256 down with the residual ≤ 10⁻¹² at every step (`crit_continuation.py` with the
repository's functions; `crit_indep.py` with the independent implementation): the class-1 margin is monotone along the branch, −0.1274,
−0.0311, −0.00655, −0.00003, +0.00123 at ε = 1, 4, 16, 64, 256 (the md: −0.17299, −0.03114, −0.04808, −0.36601, +0.00123), with one sign change at
**ε_c = 64.99** (65.02 by interpolation with the repository's functions); there **p = 0.574** (conditional plan (0.426, 0.307, 0.267)) and the
rank-one formula gives **81.9 = 1.26 ε_c**; A = σ₀I − R^G_1 ≻ 0 (min eigenvalue 0.0020) and the collective block is positive (0.0078). By the
finite-difference Hessian of the independent J^ε_m at m = 30 with counts (16, 7, 7) = 30·w_0.2 (`fd_delta02.py`): class-1 separation λ_min =
+0.000895 at ε = 136, +0.000132 at 70, −0.000029 at 64, −0.000156 at 60, −0.006547 at 16 (the formula to six digits; the full Hessian's λ_min is
the same number). So the formula is not "off" at δ = 0.2: 1.22, 1.23, 1.31, 1.26 times ε_c at the four δ (but see D3h-3).
*Fix:* in `critical_row` continue the stationary point in ε (from 256 down, and inside the bisection from the nearest converged point), assert
`res < 1e-9`, print the residuals; let `weighted_stationary_atoms_entropic` warn or raise on non-convergence; regenerate; text, §7 and RESULTS:
"ε_c = 5.1, 13.5, 23.1, 65 for δ = 0.02, 0.05, 0.1, 0.2; p = 0.057, 0.15, 0.33, 0.57; rank-one formula 6.2, 16.6, 30.3, 82 (1.22–1.31 times the
exact threshold)", and drop "off at δ = 0.2, where it charges all three".

**D3h-2 (major; PARADIGM lines 2434–2435, §0 row line 114 "at ε = 4 the larger sizes stop at the residual excess", §7 row line 4988
"residual excesses 0.017–0.029 at ε = 4", RESULTS line 472, `entropic_leak.md` table 2; the commit message).** "At ε = 4 m = 15 equalises,
while m = 60, 240 and 510 stop at the residual excesses 0.017, 0.029, 0.028 (from 0.083, 0.054, 0.041) — the critical excess at that ε" is
false: two of the three numbers are not stopping points and are not the critical excess.
(a) *The critical excess at ε = 4* (`critical_excess_eps.py`, continuation in δ at fixed ε): the over-weighted class's margin at Θ*_4(w)
vanishes at \|w − q\|_∞ = **0.0160** for the w_δ family and **0.0117** for w = q + δ(e₁ − e₃); the md's own table has ε_c(0.02) = 5.1 > 4
(margin −0.00154 at ε = 4).
(b) *Per seed* (`recompute_runs.py`, `seed_margins.py`, `final_drift.py`). The 0.029 is m = 240 seed 1, whose counts (77, 76, 87) never
changed in 600 slow units; the 0.028 (14/510 = 0.0275) is m = 510 seed 3, (168, 184, 158), which never changed either. Both are **unstable
weights** — class margins −0.00136 and −0.00146 at Θ*_4(w), e-folding times 735 and 685 slow units, longer than the run — and both are
splitting at s = 600 (max whitened drift 8.9·10⁻³ and 1.8·10⁻², a particle 0.92 and 1.24 G⁻¹-units from its class mean, against
1.5–3.1·10⁻⁵ for the three larger systems at rest). The "(from 0.083, 0.054, 0.041)" are the start maxima of *other* seeds for m = 240, 510
(seed 0 in both), and those two are still changing counts at s = 534 and 590 (m = 510 seed 0 ends at (165, 162, 183), margin −0.00042);
m = 60 seed 2 equalises only at s = 567. At s = 600, five of the twelve larger systems sit at unstable weights ((19, 21, 20)/60: −0.00035,
e-folding 2900 slow units; (77, 76, 87)/240; (165, 162, 183), (163, 166, 181), (168, 184, 158)/510: −0.00042, −0.00037, −0.00146), eight of twelve
have whitened drifts 8.7·10⁻⁴–2.1·10⁻², and the per-seed function-level distance changes by factors 0.39–2.16 between the last two 150-unit
windows.
(c) *Continued.* m = 240 seed 1 reproduced from its seed (equal to the stored final bit for bit at s = 600) and run on (`long_run.py`): the
first particle leaves class 3 at s ≈ 710, the counts reach (78, 79, 83) at s ≈ 1360 and stay there (drift 1.9·10⁻⁵ at s = 3000; the final
within 8.4·10⁻³ of Θ*_240(w), margins +0.0071, +0.0069, +0.0024 — stable): the end excess is **0.0125**, at the critical excess, which the
600-unit 0.029 is not. Only m = 60's 0.017 (one particle) is near it (seed 1 ends at the stable (20, 19, 21)/60, margin +0.00081).
*Fix:* integrate the ε = 4 (and ε = 16, 64: D3h-4) systems until counts and drift settle (several thousand slow units at margins of −0.0004 to
−0.0015), or report the s = 600 state as such ("at s = 600 the larger systems are still equalising: …"); have `report()` print per seed the end
counts, the separation margins at Θ*_ε(w_end) and the final drift (the md's max-over-seeds columns mix seeds); compare the settled excesses with
the critical excess at that ε (0.012–0.016 at ε = 4, depending on the deficit atoms); §0, §7, RESULTS accordingly.

**D3h-3 (minor; PARADIGM lines 2430–2431, §0 row line 114 "its two-atom form within 22–31% of the exact threshold", `entropic_leak.md` table 1
header, `exp_entropic_leak.py` lines 11, 97–105).** "within 22–31% for δ ≤ 0.1, where the conditional plan charges mainly two atoms". At ε_c the
class-1 conditional plans are (0.943, 0.028, 0.029), (0.849, 0.076, 0.075), (0.665, 0.185, 0.150), (0.426, 0.307, 0.267) for δ = 0.02, 0.05, 0.1,
0.2 (`crit_repo_check.py`, `crit_indep.py`): the outflow p splits 51/49, 50/50, 55/45, 54/46 between the two under-weighted atoms — as w_δ's
marginals force (atoms 2 and 3 are each short by δ/2) — so the plan charges three atoms at every δ and Cov^G_1 has rank two (eigenvalues 0.55
and 2.06, 1.47 and 4.97, 3.14 and 8.79 at δ ≤ 0.1). The script's rank-one value puts all of p on the atom with the larger share (`argmax`: atom 3
at δ = 0.02 by 0.0289 against 0.0277, atom 2 otherwise); with the other atom it gives 6.23, 15.3, 19.8, 45.4 = 1.22, 1.13, 0.86, 0.70 times ε_c.
The 22–31% is a property of the tie-break, not of an approximation the plan supports. The exact threshold is ε_c = 2σ₀λ_max(A^{−1/2}Cov^G_1A^{−1/2})
(equal to ε_c to four digits at all four δ). A two-atom plan does occur for a two-atom excess at small ε: w = q + δ(e₁ − e₃) at ε = 1 has
(0.9897, 1·10⁻¹¹, 0.0103). *Fix:* "the conditional plan sends its outflow about equally to the two under-weighted atoms, so Cov^G_1 has rank
two; the rank-one formula with all of p sent to the atom of the larger share is 1.22–1.31 times the exact threshold (0.70–1.22 with the other
atom)"; print the generalised eigenvalue; if the two-atom form is to be tested, use w = q + δ(e₁ − e₃).

**D3h-4 (minor; PARADIGM lines 2406–2407 "it splits, and the split-off part carries the excess to atom j as the OT leak's transit does",
2435–2438 "the larger sizes keep their sampled weights", "At ε = 64 no excess is equalised", §0 row line 114 "so an over-weighted class splits —
and the masses equalise — only below …", "at ε = 16 they keep their sampled weights", §7 row line 4988, RESULTS lines 472–474,
`entropic_leak.md` line 14 "the split-off particles carry the excess to the under-weighted atom").** The end-state readings at ε = 16, 64 and
the split-equalises gloss are 600-unit observations, and two are contradicted on the instance (`seed_margins.py`, `long_run.py`,
`endstate_margins.py`, `bias_fit.py`; all continued runs equal the stored finals bit for bit at s = 600).
(a) At s = 600, three of the sixteen ε = 16 systems and one of the sixteen ε = 64 systems sit at unstable weights: (6, 5, 4)/15 (seed 0,
margin −0.00083), (6, 6, 3)/15 (seed 3's end, −0.00147), (25, 17, 18)/60 (seed 3, −0.00135) at ε = 16; (8, 4, 3)/15 (seed 3, −0.00029) at
ε = 64. Continued: m = 60 seed 3 at ε = 16 sends two particles from class 1 to atom 2 (s ≈ 1200, 1490) and settles at (23, 19, 18) (margins
+0.00084, +0.0051, +0.0063; drift 3.1·10⁻⁵ at s = 4000) — not every m ≥ 60 system keeps its sampled weights; m = 15 seed 3 at ε = 64 moves one
particle (s ≈ 1360) and settles at (7, 5, 3) (within 1.8·10⁻⁵ of Θ*_15(w), drift 3.6·10⁻⁸ at s = 6000) — an excess is equalised at ε = 64.
(b) *A split that equalises nothing.* The test's weight (6, 5, 4)/15 at ε = 16, below its ε_c = 17.87 (seed 0): the class-1 spread grows from
0.06 at s = 600 to 1.2 at s = 3600, and at s ≈ 3665 the split-off particle joins **atom 2, which was at q**; the end state (5, 6, 4)/15 is stable
(margins +0.0058, +0.0027, +0.0072) and ½‖w − q‖₁ = 1/15 before and after. The conditional plan of class 1 charges atoms 2 and 3 (0.112 and
0.122 at ε_c), the unstable direction mixes both, and nothing in (iv) selects the under-weighted atom: a split moves the excess, it need not
reduce it. *Fix:* (iv) "Below ε_c the over-weighted collapsed class is a saddle and splits; its split-off particles leave toward an atom its
conditional plan charges, which need not be under-weighted ((6, 5, 4)/15 at ε = 16 ends at (5, 6, 4)/15), so a split need not reduce the
excess"; §0 without "— and the masses equalise —"; the numerics "over 600 slow units", with the slow saddles and the continued runs reported
(as D‴.5(i) already says, a saddle may also split locally without changing the weights).

**D3h-5 (minor; PARADIGM line 2401 "since σ₀I − R^G_k ≻ 0", lines 2404, 2407 "Above ε_c the class is stable", §0 row line 114).** The
positivity of A = σ₀I − R^G_k is asserted, not assumed, and it fails on the instance (`check_A_posdef.py`, independent implementation,
continuation from ε = 4096): λ_min(A_1) at Θ*_ε(w) is −0.00087 (ε = 4) and −0.00105 (ε = 1) at the test's own weight (6, 5, 4)/15 and −0.0019 at
(8, 4, 3)/15 (ε ≤ 4; seed 3's m = 15 start), positive at ε ≥ 16 there and at every ε for (7, 4, 4), (6, 4, 5), (6, 3, 6)/15 and (25, 17, 18)/60.
Where it fails, λ_min(S^ε_k) ≤ λ_min(A) < 0 whatever ε (S = A − (2σ₀/ε)Cov^G ⪯ A), so the "iff" has no content and the class is a saddle for
the fit's reason (D‴.5's splitting condition), not the entropic one; test 4's ε = 4 point is such a case (−0.00087 with the test's own loop,
n_quad = 96, as at N = 256; `test4_A.py`). At the four w_δ points A ≻ 0 holds
(min eigenvalues 0.0054, 0.0050, 0.0051, 0.0020 at ε_c), so the table is unaffected. Further: p and A depend on ε through Θ*_ε(w) and the plan,
so "ε < ε_c" is an implicit condition (on the instance the margin changes sign once along each of the five branches computed); and "the class
is stable" concerns its separation block — the replicated configuration is a local minimum only if the collective block and the other classes'
blocks are positive too (they are at the w_δ points: collective λ_min 0.0073–0.0078, other classes 0.0057–0.0091). *Fix:* "if σ₀I − R^G_k ≻ 0 (the
splitting condition of Proposition D‴.5 at Θ*_ε(w), which fails at large excesses and small ε, e.g. (6, 5, 4)/15 and (8, 4, 3)/15 at ε ≤ 4)";
"ε_c is the root of λ_min(S^ε_k) = 0 along Θ*_ε(w) (one root per branch on the instance)"; "and the configuration is a local minimum when, in
addition, the collective block is positive".

**D3h-6 (minor; PARADIGM lines 2411–2413, 2439–2440, §0 row line 114 "the floor is the nearest-atom leak's", "only the hard mass constraint is
floor-free for every m", §7 row line 4988, the D‴.20 audit bullet of `LITERATURE_TRANSPORT.md` (its Verdict, line 338), §8 line 5210).** The floor statements.
(a) "Only the hard constraint ε → 0 removes the floor for every m" / "floor-free for every m": the hard OT leak is floor-free on the commensurate
sizes (D‴.12(iii)); on the incommensurate sizes of D‴.12's own runs every plan keeps 1–2 split particles and the function-level distance on the
last window is 3.4·10⁻⁴, 3.1·10⁻⁵, 1.1·10⁻⁶ at m = 16, 64, 256 against M = 512 (`incomm_floor.py`; nearest-atom 2.5·10⁻³, 5.0·10⁻⁴, 4.4·10⁻⁵) —
a smaller, steeper floor, but a floor.
(b) "the floor is the nearest-atom leak's" (§0; "of the nearest-atom leak's order" in the numerics is fair): on the same clouds (start weights
identical in all twelve (seed, m) pairs) the ε = 16 floors are 0.23–0.84 of the nearest-atom ones per seed, 0.61, 0.43, 0.60 in seed means
(`floors_compare.py`). What returns is D‴.5's law at the *entropic* stationary point: (1/m − 1/M)tr(FC̃Fᵀ)/N with F = ∂f̂_w/∂w along zero-sum
directions at Θ*_ε(q) predicts 2.2·10⁻³, 5.0·10⁻⁴, 7.5·10⁻⁵ at ε = 16 — 0.55 of the nearest-atom law's 4.0·10⁻³, 9.1·10⁻⁴, 1.4·10⁻⁴, which the same
computation reproduces — and the observed/predicted ratios 0.71, 0.45, 0.47 are those of the nearest-atom leak on the same clouds (0.63, 0.58,
0.43; `floor_law.py`); at ε = 64 the law is 0.74 of the nearest-atom one.
(c) "with the excess absorbed by the class biases instead of conserved": at large m the entropic leak conserves the cluster masses (in the
600 slow units every m ≥ 60 run at ε = 16 and 64 keeps its counts) — the excess is conserved, as under the nearest-atom leak, and the atoms are
biased. *Fix:* "only the hard constraint removes the floor (on the commensurate sizes; elsewhere its K₀ − 1 bridge particles leave an
O(m⁻²) floor)"; "the floor is Proposition D‴.5's law at the entropic stationary point (0.55 of the nearest-atom leak's at ε = 16)"; "the excess
is conserved and partly absorbed by the class biases".

**D3h-7 (cosmetic; PARADIGM lines 2400–2401, 2407, 2410–2412).** In (v), p denotes first the critical *conditional fraction*
(p_c ≈ ε/(2‖Δ̂‖²)) and then the *mass* excess of the sample ("p ~ m^{−1/2}"). At small ε the plan sends exactly the excess, so the critical mass
excess is w_kp_c (`critical_excess_eps.py`: δ_c = 0.0042, 0.0083, 0.0160, 0.0296 = w₁p_c at ε = 1, 2, 4, 8; p_c/(ε/(2‖Δ̂₁₃‖²)) = 1.12, 1.08,
1.02, 0.92 for w_δ and 0.92, 0.85, 0.76, 0.74 for the two-atom family, where R ≠ 0 lowers it), and with the multinomial s.d. √(q(1−q)/m) the
threshold is m ≫ ((1 − q)/q)(2‖Δ̂‖²/ε)² = 2(2‖Δ̂‖²/ε)² — the scaling is right, the constant and the symbol are not. "T_k = (1 − p)ϑ⁰_k + pϑ⁰_j" is
the two-atom form, never realised in the numerics (D3h-3). *Fix:* distinguish p_c from the excess δ_c = w_kp_c.

**D3h-8 (cosmetic; PARADIGM lines 2420–2422, §0 row line 114, LITERATURE lines 66–67 and 323–324).** (a) "with the roles exchanged: here the
particles are softly assigned to the prior's atoms, and what splits is a cluster of particles" — in deterministic annealing the moving points
(codevectors) are the ones that split, and the covariance is that of the fixed points (data) under the conditional law given the moving one
(C_{x\|y}); here the moving points (particles) split and the covariance is that of the fixed points (atoms) under π(·\|k). The correspondence is
direct (particles ↔ codevectors, atoms ↔ data, ε ↔ T, Cov_{π(·\|k)}(ϑ⁰) ↔ C_{x\|y}); what differs is the two-sided mass constraint (the
conditional covariance of a collapsed class comes from an excess, not from a cluster's data spread) and the fit's curvature. (b) "a critical
regularisation equal to twice the largest eigenvalue of a conditional covariance" holds only at R = 0: in general ε_c = 2σ₀λ_max(A^{−1/2}Cov^G A^{−1/2});
on the instance 2λ_max(Cov^G_1) = 4.11, 9.95, 17.6 against ε_c = 5.11, 13.5, 23.1 (δ ≤ 0.1). (c) The ε → 0 / ε → ∞ interpolation of (i)
(OT to a barycentre/MMD-type term) has a precedent in the Sinkhorn-divergence literature (Ramdas–García Trillos–Cuturi 2017; Genevay–Peyré–Cuturi
2018; Feydy et al. 2019, "Interpolating between optimal transport and MMD using Sinkhorn divergences"), which also gives the standard removal of
the barycentric bias that dominates at ε = 64; and size-constrained soft assignment by Sinkhorn (equipartition, e.g. Asano–Rupprecht–Vedaldi
2020) is related prior art for the soft mass constraint (references from memory). *Fix:* "with the same roles (particles as codevectors, atoms
as data); what is new is the two-sided mass constraint, which makes the conditional covariance of a collapsed class an effect of its excess";
"twice the largest generalised eigenvalue of Cov^G relative to I − R^G/σ₀"; add the Sinkhorn-divergence references to [S45].

**D3h-9 (cosmetic).** `entropic_leak.md` line 3 says "the entropic plan by log-domain Sinkhorn (tolerance 10⁻¹³)" and the module docstring
(line 12) "implements the log-domain Sinkhorn plan" — for three columns the code uses Newton on the semi-dual, as PARADIGM and RESULTS say;
"0.028" is 14/510 = 0.0275 (0.027); "the strictly concave semi-dual" is strictly concave modulo constants (fix the gauge); "every bounded
trajectory" — boundedness is automatic (one line); the ε = 1 "function-level distances 1.0–5.6·10⁻⁷" are still falling (last/previous
150-unit window 0.008–0.15 per seed, 0.02–0.07 in seed means; final drifts 5.8·10⁻⁴–1.3·10⁻²) — a stopping-time value, as in D3f-8; the
experiment docstring's "the receiving atom" (two receive, D3h-3); the md's dynamics table reports max-over-seeds at the end and at the start,
which come from different seeds (D3h-2).

**D3h-10 (cosmetic; `tests/test_entropic_leak.py`).** The tests pass and discriminate, but: test 2's ε → 0 assertion holds with slack 10¹⁷
(6.5·10⁻²¹ against 10⁻³; it would still pass at ε = 0.2, 4·10⁻⁵, and fails at 0.5, 4.8·10⁻²) and its "unique plan here" is a comment, not a
check; test 3 compares sorted spectra (tolerance 10⁻⁶ against 1.7·10⁻¹² achieved) rather than the blocks and the invariance of the separation
subspace (both hold to 10⁻¹², `check_i_ii_iii.py`); test 4 ignores the Newton residuals (here 1.5·10⁻¹³ at worst — harmless, but the same
unchecked path produced D3h-1) and its ε = 4 point has σ₀I − R^G_1 indefinite (−0.00087, D3h-5). *Fix:* assert `res < 1e-9`; compare
BᵀH_wB with I ⊗ S_k entrywise and bound ‖(I − P)H_wB‖; test the ε → 0 limit at an ε where the error is near the tolerance.

**D3h-11 (cosmetic; `LITERATURE_TRANSPORT.md` lines 323–338 at c29c34a).** The commit that added Theorem D‴.21 (c29c34a, landed during this
review) inserted its audit bullet between the D‴.20 bullet (lines 323–327) and the D‴.20 Verdict line, which now sits at line 338 under the
D‴.21 bullet — D‴.21 has two Verdict lines and D‴.20 none. *Fix:* move "Verdict: standard ingredients, new object and new consequence
(only the hard constraint is floor-free for every m)" back to line 328 (and amend it per D3h-6).

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_entropic_leak.py --from-raw` (scratch copy) rewrites `entropic_leak.md` byte-identically (md5 d34b8ec7…).
2. **Every number of `entropic_leak.md` and of the prose, §0, §7, RESULTS** (`recompute_runs.py`, `inspect_raw.py`): ε_c 5.107, 13.508,
   23.105, 136.138; p 0.0566, 0.1515, 0.3345, 0.6233; rank-one 6.218, 16.584, 30.275, 88.091 (ratios 1.218, 1.228, 1.310 → "22–31%"); every grid
   margin; the 16 dynamics rows (max end/start excesses, fit errors, W², function-level distances on (450, 600]); the prose ranges 1.7–3.0·10⁻⁴,
   1.0–5.6·10⁻⁷, 8.8–9.6·10⁻³, 4.5–7.7·10⁻⁵, 0.20 → 0.13, 1.5·10⁻³/2.2·10⁻⁴/3.5·10⁻⁵, 0.033–0.069, 0.12–0.14; the D‴.12 nearest-atom floors 2.52·10⁻³,
   5.22·10⁻⁴, 5.86·10⁻⁵ on the same clouds.
3. **The runs are what the script says.** Four systems reproduced from their seeds with the script's scheme equal the stored float32 finals bit for
   bit at s = 600 (m = 240 seed 1 ε = 4; m = 60 seed 3 ε = 16; m = 15 seeds 0 and 3 at ε = 16, 64), and a fifth (m = 510 seed 2, ε = 1) as well.
4. **The solver.** Column marginals ≤ 1·10⁻¹³ at every plan refresh of the four continued runs (to s = 3000–6000); 3006 solves of m = 510 at ε = 1:
   median 3, max 7 Newton iterations, no `max_iter`, no `lstsq` fallback; cold-started plans at all 64 stored finals: column error ≤ 7·10⁻¹⁴.
5. **(i)–(iii)** as in the verdict table (gradient 2–3·10⁻¹¹, envelope 5·10⁻¹⁰, primal = dual, the limits, the semi-dual kernel, class symmetry
   2·10⁻¹⁷, the separation blocks entrywise to 1.4·10⁻¹² with invariance 2.9·10⁻¹², the collective block = the K-atom Hessian to 7.9·10⁻¹²).
6. **The critical excess for δ ≤ 0.1** by the independent implementation (`crit_indep.py`): ε_c = 5.1066, 13.5077, 23.1045, p, rank-one values and
   grid margins to the printed digits; A ≻ 0, collective block and other classes positive at ε_c; one sign change per branch; the exact threshold
   2σ₀λ_max(A^{−1/2}Cov^GA^{−1/2}) equals ε_c to four digits. The test weight (6, 5, 4)/15: ε_c = 17.87 (N = 256), ∈ (17, 18) at n_quad = 96.
7. **(iv)'s first sentence.** At w = q the entropic term is 3·10⁻¹⁸–7·10⁻¹⁶ of σ₀ at ε = 1 and 3·10⁻⁴–1.1·10⁻³ at ε = 4 (off-diagonal conditional
   masses 1.6–6.9·10⁻⁵ ≈ e^{−Δ²/ε}); at ε = 16 (not ≪ Δ_ν² = 39.2) it is 41–45%.
8. **(v)'s scaling.** δ_c(ε) = 0.0042, 0.0083, 0.016, 0.030, 0.062, 0.14, 0.20 at ε = 1, 2, 4, 8, 16, 32, 64 (w_δ family); p_c ≈ ε/(2‖Δ̂‖²) within
   8–26% at ε ≤ 8.
9. **"the bias toward the barycentre dominates" at ε = 64.** The balanced Θ*_ε(q) alone has fit error 0.120 at ε = 64 (runs 0.12–0.14), 0.029 at
   ε = 16, 1.4·10⁻⁵ at ε = 4.
10. **The ε = 1 statements** (every one of the 16 systems ends at q; last count change s = 117–279, 556 for m = 510 seed 2) and "only m = 15 moves" at
    ε = 16, "no count changes" at ε = 64, over 600 slow units.
11. **The tests**: 4 passed in 1.25 s; the actual values behind every assertion (`check_tests.py`): marginals 10⁻¹⁷/10⁻¹⁴, solver agreement 3·10⁻¹⁵,
    gradient 8·10⁻¹¹, limits 6.5·10⁻²¹ and 6.6·10⁻⁴, criticality 3·10⁻¹⁷, symmetry 0, blocks 1.7·10⁻¹², margins at ε = 4, 16, 20, 64: −0.0224,
    −0.00083, +0.00076, +0.0056.

## Not checked

The full experiment was not rerun (deterministic; the 16 jobs at m ≤ 510): five systems were reproduced from their seeds and four continued
(to s = 3000–6000); the other systems at unstable weights at s = 600 (m = 60 seed 0 and the m = 510 systems at ε = 4, m = 15 seed 3 at ε = 16) and
the M = 510 references were not continued, so the floors after the slow splits were not recomputed. The entropic weight-sensitivity law was
evaluated at w̄ = q, not per realised weight. σ₀I − R^G_k was scanned at six weights and seven values of ε. The literature statements
(Rose et al.'s T_c = 2λ_max(C_{x\|y}); Cominetti–San Martín's limit as the maximal-entropy optimal plan; Courty et al.'s barycentric mapping;
the Sinkhorn-divergence and equipartition references) were checked against memory, not against the papers (no web access was used), and the
novelty claim was not searched. The analyticity of the period average for a continuous occupation measure was argued, not computed (the
numerics use the 256-point quadrature).
