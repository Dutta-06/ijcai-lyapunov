# Referee report D3d — PARADIGM §5d, the material new since round 10: Proposition D‴.12 (i)–(v) and its Status, Remark D‴.13, the D‴.12 *Numerics* block (tables 1–4, Reading), the §0/§7/§8 clauses on D‴.12, `results/RESULTS_PARADIGM.md` (D‴.12), `LITERATURE_TRANSPORT.md` [S40]–[S42] and the audit bullet; `src/core/ot_leak.py`; `experiments/exp_ot_leak.py` → `results/chaos/ot_leak.md`, `ot_leak_raw.npz`; `tests/test_ot_leak.py`

Branch `claude/festive-bardeen-eopej2`. The review is pinned to commit **c2e1ac6** ("Proposition D‴.12, part 2 — numerics"), the HEAD when
the brief was issued (`git log --oneline -3`: c2e1ac6, 056047d, ff2854f). While the review ran, the branch moved to fb23968 (Corollary D‴.14,
part 1: `experiments/exp_ot_real_loop.py`, two D‴.14 tests appended to `tests/test_ot_leak.py`) and the working tree acquired further
uncommitted edits by the author (25 lines appended to `src/core/ot_leak.py` — two helper functions after `bridge_report`, nothing in the
reviewed code changed — a third appended test, and untracked files); none of these is reviewed here, and none was made by the referee. No
repository file was modified by the review: the `--from-raw` regeneration was run with the output directory redirected to the scratch
directory (md5 of the regenerated `ot_leak.md` `b537eab2…` = the repository's), pytest was run with `-p no:cacheprovider` and
`PYTHONDONTWRITEBYTECODE=1`, and the c2e1ac6 versions of the five scope files are snapshotted in `…/ref_D3d/c2e1ac6/`. Scratch scripts and
their outputs (`.out`) are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3d/`:
`inspect_npz.py` (the raw data's structure), `recompute_from_npz.py` (every number of tables 1–4, the time courses, the bridge list and the
prose, with per-seed values, decay rates and the D(OT) time course), `check_code.py` (LSA-with-replication against the transportation LP,
the forest structure of the LP plans, Danskin by finite differences, drift = −G∇J, `W2_empirical` by three routes, the Hessian of J₁₅^{OT}
at Θ*₁₅(q), `transit_mask` on tied rows, the balanced-sample probability), `check_dynamics.py` (drift at the 64 stored final
configurations, within-class spreads, the exact cross-cluster W² and D laws of the nearest-atom runs, the D(OT) plateau from the mismatched
pairs, W²(μ_m, ν) and its m-law, the plan at s = 0 from the regenerated clouds), `check_followup.py` (migrated particles at the commensurate
finals, the floors against μ*(q) = ν with their m-slopes, the constants of (iv)'s bounds, the design table from scratch), `check_fnlaw.py`
(the exact per-seed function-level law of the nearest-atom runs; the m = 16 transit series), `fromraw_wrapper.py` (the regeneration),
`continue_m16.py` → `continue_seed{0,1,2,3}.out/.npz` (the incommensurate OT runs of m = 16, 64 against M = 512 reproduced from their seeds
and continued to s = 2400) with `analyse_continue.py` (their plateaus against the reported last-window values), `pytest_head.out`,
`pytest_c2e1ac6.out`. The full experiment (1.8 h in the author's log) was not rerun; the
m-systems and the references are autonomous under either leak, so the runs were reproduced from their seeds (`mixture_start`).

## Verdicts

| object | verdict |
|---|---|
| **Proposition D‴.12 (i)** | **holds.** On 𝒢_S under D‴.9(i)'s hypothesis (for every s) no particle changes cell, so w^m(s) = w^m(0); the numerics confirm max_s‖w(s) − w(0)‖_∞ = 0 in all 32 nearest-atom systems. |
| **Proposition D‴.12 (ii)** | **holds as a conditional statement, but its hypothesis is incompatible with Theorem D‴'s i.i.d. start except on an event of probability ≈ 0.83/m (D3d-1, major).** Total unimodularity and the stability of a unique optimal vertex are used correctly (confirmed: the LSA-with-replication plan equals the LP optimum to 5·10⁻¹⁵ and is integral on 9 commensurate configurations; the LP plans are forests in all 28 configurations tested); Danskin holds (grad_J against finite differences of J to 3·10⁻¹⁰–3·10⁻⁹ relative; drift = −G∇J exactly); the leak enters Theorem D‴'s proof only through the common target of ϑ_i and ϑ̄_i (Steps 1–2), the a priori bound (Step 0, R_Θ with the hull of supp ν) and V_max — exactly what the stated hypothesis provides. But the hypothesis "for every s ≥ 0 the plan assigns each particle to the atom of its reference's Laguerre cell" forces w^m(0) = q: in all 32 OT systems of the numerics the plan at the first record assigns 1–21 particles away from their nearest atom, while the mean-field flow of the mixture start has masses q from the start, and under the shared-particle coupling the mismatch is permanent. The Status's open item ("the coupling *through* the transit") understates this. The non-uniqueness of the Laguerre cells for a disconnected support is a cosmetic gap (D3d-15). |
| **Proposition D‴.12 (iii)** | **holds.** The cycle argument is right (every feasible γ′ − γ_D is a circulation; a conformal cycle decomposition has as many off-diagonal as diagonal edges per cycle, each off-diagonal pair costing ≥ (Δ_ν − ε)² > ε² ≥ each diagonal pair), uniqueness survives coincident particles (the matrix γ is unique even though swapping identical rows is a symmetry), and the Hessian is Proposition D‴.5's with the plan frozen: at Θ*₁₅(q) the FD Hessian of J₁₅^{OT} equals the nearest-atom J₁₅^ν's to 1.7·10⁻¹⁸, with whitened eigenvalues σ₀ (multiplicity 36 = 3 classes × (5 − 1) × 3 coordinates) and the nine collective eigenvalues 0.00738…0.20229 of round 10, |drift| = 1.3·10⁻¹⁷. "Strict local minimiser iff q ∈ W_adm" has the usual degenerate-boundary caveat (as D‴.5(i)). The "W₂-distance ½‖w^m − w^M‖₁·O(Δ_ν²)" is the *squared* distance (D3d-7). |
| **Proposition D‴.12 (iv)** | **holds (the count proved, the bounds formal, as stated).** Split rows = Σ_i(deg_i − 1) ≤ (m + K₀ − 1) − m = K₀ − 1 is right; the two bounds are valid upper bounds on a specific coupling and on a direct function-level estimate, with constants vacuous on the instance (4(K₀ − 1)R_Θ²/m = 4·10⁶ at m = 16 against the measured 1.0–1.4, the function-level 1.00 against 0.018, D3d-12). "The W₂ floor keeps the order 1/m" compares with the wrong law — the nearest-atom leak's W² floor is O(m^{−1/2}) by (iii) (D3d-7). The orders are confirmed by the right estimators: W²(μ*_m, ν) = 19.3/m, 20.9/m, 20.0/m (slope −0.98 ± 0.03) and ‖f̂*_m − f‖/‖f‖ = 0.45/m–0.60/m (squared slope −1.81 ± 0.05; −1.91 ± 0.07 with M = 512). |
| **Proposition D‴.12 (v)** | **holds with one false clause.** The transit (ln 2/σ₀ = 98 against s_eq = 93–164, median 103.5), D = O(Δ_ν²‖w^m − w^M‖₁) and the orders of W² and the function level during the transit are as stated; "after the transit all three vanish on the commensurate sizes" is false for D, which converges from below to the permanent-mismatch plateau Σp_{kk′}Δ²_{kk′} (D3d-2). The mean-field transit claim is not checked by the numerics (the mixture start has w = q in the mean field; the M = 510/512 references transit because they are finite samples). |
| *Status* | **holds with fixes.** The proof sketches of (ii)–(iv) are correct; "(iii) proved", "(iv) the count proved, the two bounds formal", "(v) formal" are right; the open item should say that the synchronous coupling never recovers from the transit and that (ii) applies to balanced samples (D3d-1). |
| *Remark D‴.13* | **holds with a nuance.** "The memory is erased" holds for the stationary measure; the synchronous coupling keeps the memory for ever (D(OT) at the end is 0.97 of the plateau set by the particles reassigned in one system and not in the other, D3d-2), so the memory moves from the measure to the coupling. (a)–(c) are fine. |
| *Numerics* block, tables, *Reading* | **hold with fixes.** Every number of tables 1–4, the time courses, the bridge list and the design table reproduces from the raw data, the design table also from scratch, and `--from-raw` is byte-identical. Over- or mis-statements: "one to twelve particles" (1–21, D3d-3); the fit-error explanation by the jitter's spread 0.02e^{−σ₀·600} ≈ 3·10⁻⁴ — the residual spread is the migrated particles' Δ_ν e^{−σ₀·600} ≈ 0.09 (D3d-4); the m-laws of (iv) quoted from distances to a reference that carries its own bridge particles (D3d-5); "stationary" for the m = 16 function-level value (0.82 per window; its plateau, reached by s ≈ 700–960 in the continued runs, is 21% below the quoted value, D3d-6); "within 5–14%" (5–16% of the formula, D3d-11); "about 0.5δ" (0.49–0.69δ, D3d-10); "the two or three particles to move are the split ones" (D3d-9); the W₂/W² wording (D3d-7). |
| §0 row, §7 row, §8 bullets, RESULTS, README (c2e1ac6) | **consistent with the data** (checked clause by clause); the §0 clause "the frozen-plan and stationary statements are proved" and the §8 open item inherit D3d-1; "O(1/m) in W₂" means W² (D3d-7); "one to twelve" is repeated in RESULTS (D3d-3). |
| LITERATURE [S40]–[S42], audit bullet | **plausible and honestly hedged** ("details from memory; not re-verified"); the attributions (Hoffman–Kruskal, Klee–Witzgall, Villani Thm 5.10, Aurenhammer–Hoffmann–Aronov, Kitagawa–Mérigot–Thibert) match the uses; bibliographic details not verified here. |
| `src/core/ot_leak.py` | **holds.** `plan` (LSA on replicated atoms = LP optimum, integral; LP vertices are forests), `leak_target`, `drift` = −G∇J (Danskin confirmed), `J`, `split_mask`, `W2_empirical` (replication = LP = lcm-assignment to 5·10⁻⁹, symmetric to 6·10⁻¹⁰), `diagonal_plan_margin`, `stationary_atoms` (its ε < Δ_ν/2 output is the cell check, left to the caller), `coupled_run_ot` (vertex-selection with switches at multiples of 0.2 slow units; records as documented), `equalisation_time`, `bridge_report` do what they say. `transit_mask` is undefined on a tied split row and counts the three-way-split bridge particle as "in transit" for ever (D3d-8). |
| `experiments/exp_ot_leak.py` | **holds with minor fixes.** `report()` reproduced (byte-identical); the prose it feeds should quote the per-seed transit counts (D3d-3), the exact per-seed laws it has the ingredients for (D3d-11, D3d-14), and the distances to ν for the incommensurate m-laws (D3d-5). |
| `tests/test_ot_leak.py` (c2e1ac6) | **passes (5/5 in 0.94 s) and tests what its docstring says.** Test 2's FD tolerance has slack 1.4·10⁴; test 1's W² identity is a tautology of `plan`; `CAB` in test 5 is dead code (D3d-13). The working-tree file (8 tests, D‴.14 included) passes in 54 s; not reviewed. |

## Confirmed issues

**D3d-1 (major; statement (ii), Status "(ii) proved" and "Open", §0 row "the frozen-plan … statements are proved", §8 bullet).** The
hypothesis of (ii) — "for every s ≥ 0 the plan assigns each particle to the atom of its reference's Laguerre cell L_·[μ_s]" — must hold at
s = 0. The reference's cells are those of the mean-field flow μ_s, which for the mixture start has masses exactly q from s = 0 (so no mass
migrates in the mean field and the cells are the clusters' own), while the m-system's vertex plan at s = 0 must move ½‖w^m(0) − q‖₁ of mass
to atoms other than the nearest. Hence the hypothesis holds at s = 0 only on {w^m(0) = q}, whose probability under Theorem D‴'s i.i.d. start
is the multinomial P(m/3, m/3, m/3) = 0.0527, 0.0136, 0.00344, 0.00162 at m = 15, 60, 240, 510 (≈ 3√3/(2πm) → 0). On the numerics' own
clouds the plan at the first record assigns 1, 1, 13, 21 | 1, 3, 11, 21 | 2, 5, 7, 11 | 2, 4, 4, 11 | 1, 3, 4, 3 | 2, 3, 5, 3 | 3, 5, 2, 14 |
4, 6, 4, 14 particles (m = 15/16, 60/64, 240/256, 510/512; eight (seed, set) pairs) away from their nearest atom — the hypothesis fails in all
32 systems — and under the shared-particle coupling the pairs reassigned in one system and not in the other never re-match: at s = 600 the
m-systems have 2, 3, 8 | 2, 5, 1 | 1, 3, 6 | 3, 7, 7 (commensurate; seeds 0–3) mismatched pairs out of m, and D sits at 0.970–0.979 of the
plateau Σ_{mismatched}Δ²_{kk′}/m in all twelve commensurate pairs (D3d-2). So "Theorem D‴, Proposition D‴.5 and Corollary D‴.7 hold for
(𝔄_m^{OT}) with the same constants" is, for i.i.d. initial data, a statement about an event of probability O(1/m): vacuous in the limit that
propagation of chaos is about. Nothing in (ii), the Status or the §0 row says so; the Status's open item ("the uniform-in-time coupling
*through* the transit") suggests the hypothesis fails only during the transit. *Fix:* state in (ii) that the hypothesis entails w^m(0) = q
and the nearest-atom assignment at s = 0 — (ii) is Theorem D‴ for a *balanced* sample, which Theorem D‴'s Step 8 (McDiarmid for i.i.d.
particles) does not cover as written (a stratified start with mq_k particles i.i.d. within each class would need the bounded-differences step
redone class by class) — give P(w^m(0) = q), and say in Status/§0/§8 that for the generic unbalanced sample the synchronous coupling carries
the permanent mismatch and what is proved about the OT dynamics is the stationary statement (iii), the dynamical statement being open for
any coupling.

**D3d-2 (minor; statement (v), Remark D‴.13 nuance).** "after the transit all three vanish on the commensurate sizes ((iii))" is false for
the synchronous D, which the Numerics paragraph itself calls "the one distance the OT leak does not improve": the OT window means of D are
0.766, 2.078, 4.163, 5.344, 5.765 (m = 15), 0.351, 1.186, 2.322, 3.036, 3.316 (m = 60), 0.070, 0.290, 0.629, 0.855, 0.943 (m = 240) over the
five windows — still rising by 8–10% in the last window (last/previous 1.079, 1.092, 1.104) — against the stationary nearest-atom values
0.871, 0.252, 0.017. D converges from below to the plateau (1/m)Σ_{i : label_m(i) ≠ label_M(i)}Δ²_{kk′} set by the permanently mismatched
pairs: end/plateau = 0.970–0.979 in all twelve commensurate pairs (the migrated particles are still Δe^{−σ₀·600} ≈ 0.09 from their atoms,
D3d-4). The memory of the initial masses that Remark D‴.13 says the OT leak erases is erased from the stationary *measure* and kept by the
*coupling*. *Fix:* "the two measure-level distances vanish ((iii)); D converges to Σ_{k≠k′}p_{kk′}Δ²_{kk′}, the mass reassigned in one system
and not in the other (0.97 of it at s = 600)".

**D3d-3 (cosmetic; Numerics ×2, Reading, `RESULTS_PARADIGM.md`).** "(m = 15, 60, 240, 510: one to twelve particles)" and "moves one to
twelve particles across": the mass in transit at s = 0 is, per system, 1, 2, 1, 3 (m = 15), 1, 5, 3, 5 (60), 13, 7, 4, 2 (240), 21, 11, 3, 14
(510) particles — one to 21; "twelve" is the seed mean at M = 510 (0.0240·510 = 12.2). Likewise the fit errors "0.067, 0.029, 0.017, 0.012
… to 2.8·10⁻⁴, 1.2·10⁻⁴, 1.2·10⁻⁴, 0.8·10⁻⁴" are seed means quoted as values (per seed 0.029–0.119 at m = 15 and 8.9·10⁻⁵–5.8·10⁻⁴ under
the OT leak; the table says "mean over seeds", the prose does not). *Fix:* "one to 21 particles (seed means 1.75, 3.5, 6.5, 12.2)"; "seed
means".

**D3d-4 (minor; Numerics, the fit-error explanation; the h² decay gloss).** "the residual within-class spread (0.02·e^{−σ₀·600} ≈
3·10⁻⁴)" is off by two orders of magnitude. At s = 600 the migrated particles — 1–21 per system — still sit 0.064–0.157 whitened units from
their new atom (0.088–0.106 in 55 of 64 cases), which is Δ_{kk′}e^{−σ₀·600} = 0.090, 0.096, 0.103 for the three atom pairs: the pure leak's
approach, as (v)'s ln 2/σ₀ says; and the non-migrated particles are displaced by 10⁻³–10⁻¹ (rms 0.001–0.058 per class; seed 3, m = 15:
0.058) because the fit force makes the class compensate the misplaced kernel. The within-class spread of the end configurations is therefore
≈ 0.09, not 3·10⁻⁴; the fit error is 10⁻⁴ because the fit force cancels the first-order effect of the displaced particles (the text's own
"class means agreeing to first order"), leaving a second-order residual. Consistently, the OT finals are not stationary: max whitened |drift|
5.7·10⁻⁴–1.1·10⁻³ in all 16 commensurate OT systems against 7.9·10⁻⁶–1.0·10⁻³ (15 of 16 below 7·10⁻⁴) under the nearest-atom leak. The decay
gloss "e^{−4σ₀·150} = 0.014: second order in the within-class spread" is right in mechanism and approximately in rate (h² decays at 3.8σ₀,
3.0σ₀, 3.0σ₀ on (300, 600], factors 0.018, 0.043, 0.042 per 150; W² at 2.00σ₀, factors 0.119–0.122 ✓). *Fix:* "the migrated particles'
remaining approach Δ_ν e^{−σ₀·600} ≈ 0.09, whose first-order effect the fit cancels".

**D3d-5 (minor; Numerics incommensurate, §0 row, §7 row, RESULTS: the m-laws of (iv)).** The slopes quoted as "the O(1/m²) of (iv)"
(−2.07 ± 0.05) and "the order 1/m of (iv)" (−1.21 ± 0.04, five standard errors from −1; per seed −1.20, −1.31, −1.15, −1.17) are measured
against the M = 512 reference, which carries its own two bridge particles (W²(μ_M, ν) = 0.037–0.040, fit error 0.0011). At m = 256 the
reference's remainder is comparable to the m-system's and the two systems' bridge particles partly match: W²(μ₂₅₆, μ₅₁₂) = 0.036–0.056 is
*below* W²(μ₂₅₆, ν) = 0.074–0.085, which steepens the apparent law. The quantities (iv) is about — the distances to the mean-field stationary
measure μ*(q), which on this instance is ν itself — are in the raw data (`W2nu`, `fit`) and give clean laws: W²(μ*_m, ν) = 1.204, 0.326,
0.0780, 0.0376 (m = 16, 64, 256, 512), i.e. m·W² = 19.3, 20.9, 20.0, 19.2 ≈ 0.5Δ_ν² with slope −0.98 ± 0.03 (−1.00 ± 0.03 with M); and
‖f̂*_m − f‖/‖f‖ = fit error = 0.0282, 0.0093, 0.0023, 0.0011, i.e. m·fit = 0.45, 0.60, 0.60, 0.54, squared slope −1.81 ± 0.05 (−1.91 ± 0.07 with
M). *Fix:* report these as the test of (iv) and keep the reference-based distances for the comparison with the nearest-atom leak.

**D3d-6 (minor; Numerics incommensurate, "The stationary distances (last/previous window 0.82–1.00)").** At m = 16 the OT function-level
distance is not stationary: 0.819 per window, log-slope −0.00108/slow unit on (450, 600] (e-folding ≈ 900 slow units), seed-mean 3.17·10⁻⁴
at s = 600; the bridge particles' whitened drift is 1.6–1.9·10⁻³ at m = 16 (2–7·10⁻⁴ at m = 256). Continued to s = 2400 from the seeds
(`continue_m16.py`, all four seeds, m = 16 and 64 against M = 512): the m = 16 function-level distance settles at 2.455, 3.132, 3.132,
2.554·10⁻⁴ (seeds 0–3; within 2% of its plateau from s = 710, 961, 961, 709), mean 2.82·10⁻⁴ — the reported 3.41·10⁻⁴ is 21% above it, and
the ratio OT/nearest at m = 16 is 0.114, not 0.14; at m = 64 the plateau is 2.85·10⁻⁵ (reported 3.12·10⁻⁵, 9% high). W² is stationary
(plateaus 1.163 and 0.293 against the reported 1.157 and 0.285), the fit errors settle at 0.02677 (m = 16) and 0.00915 (64) against 0.0282
and 0.0093, and D(OT) keeps rising to 5.11 and 2.92 against 4.65 and 2.80 (D3d-2). With the plateaus the seed-mean h² slope is −2.00
instead of −2.06 (W²: −1.21 either way), so (iv)'s order stands. At s = 2400 the bridge particles are stationary (whitened drift
5·10⁻⁷–3·10⁻⁶ at m = 16, ≤ 8·10⁻⁸ at m = 64, ≤ 4·10⁻⁹ at M = 512) at 0.254–0.551 Δ_ν from the nearest atom and 0.011–0.557 from their target.
*Fix:* quote the m = 16 entry as still relaxing (or extend its horizon) and correct the ratio.

**D3d-7 (cosmetic; (iii), (iv), §0 row, Reading: W₂ against W²).** (iii) "at W₂-distance ½‖w^m − w^M‖₁·O(Δ_ν²) = O(Δ_ν²m^{−1/2})" is the
squared distance (the tables' W²; W₂ itself is O(Δ_νm^{−1/4})); (iv) "the W₂ floor keeps the order 1/m, with a deterministic quantisation
constant in place of the fluctuation's" compares with the D-law, but the nearest-atom leak's W² floor is O(m^{−1/2}) by (iii) (measured slopes
−0.67 ± 0.13, −0.79 ± 0.12), so the OT leak gains half an order in W² (−0.98) and a full order at function level; "O(1/m) in W₂" (§0) and "a
W₂ floor of order m^{−1/2}" (Reading) also mean W². *Fix:* write W² (or "squared W₂") throughout and state the comparison m^{−1/2} → m^{−1}.

**D3d-8 (cosmetic; `transit_mask`, table 1 incommensurate).** A three-way-split row ((1/3, 1/3, 1/3)/m — the single bridge particle of
seeds 0 and 3 at m = 16, which sits near the centroid, 0.58 Δ_ν from the nearest atom) has no assigned atom; `np.argmax` picks index 2 by
round-off while the nearest atom is index 1, so the particle counts as "in transit" at every record (seed 0: transit ≡ 1 particle = nsplit
throughout). The table's "mass in transit … at the end 0.0312" at m = 16 (= 1 particle in two of four seeds) and three of the four s_tr = 0
entries are this artefact; and the column header "mass in transit at s = 0 (= ½|w(0) − q|₁)" is an identity on the commensurate sizes only
(exact in all 16 systems there) — on the incommensurate sizes the first-record transit is 0.1406 against ½|w(0) − q|₁ = 0.1094 at m = 16
(per system 1, 2, 2, 4 particles against 0.67, 1.33, 1.33, 3.67), as (v)'s "(up to the split particles)" says and the header does not.
*Fix:* exclude split rows from `transit_mask` (or label the column "unsplit particles in transit"), and qualify the header.

**D3d-9 (cosmetic; Numerics incommensurate).** "at m = 16 the two or three particles to move are the split ones": per seed the particles in
transit at the first record are 1, 2, 2, 4 and the split ones 1, 2, 2, 2; in seeds 0–2 the one or two transit particles are the split ones
(s_tr = 0.02, by D3d-8's convention), in seed 3 (w(0) = (9, 4, 3)/16) two unsplit particles migrate: transit 4 → 3 (s = 90) → 2 (154) → 1
(162), nsplit 2 → 1 at s = 110.

**D3d-10 (cosmetic; Numerics table 4).** "the bias of the mass constraint, about 0.5δ": fit/δ = 0.69, 0.66, 0.63, 0.57, 0.49 at δ = 0.01,
0.02, 0.03, 0.05, 0.1 — "0.5–0.7δ, decreasing with δ".

**D3d-11 (cosmetic; Numerics commensurate, §7 row, `report()`).** "½‖w^m − w^M‖₁Δ_ν² = 4.9, 2.5, 0.70, within 5–14%": measured/formula =
1.054, 1.079, 1.158 (5–16% of the formula; 5–14% of the measured). The exact cross-cluster transport — the 3 × 3 optimal plan between the two
weighted stationary measures Σ_kw^m_kδ_{ϑ*_k(w^m)} and Σ_kw^M_kδ_{ϑ*_k(w^M)} at the realised weights — reproduces the measured last-window W²
with ratios 0.993–1.001 in 17 of the 24 (seed, m) pairs, 0.950–0.974 in three, and 0.846, 0.859, 0.906, 0.911, 0.933 at the saddle weights
(8, 4, 3)/15, (9, 4, 3)/16, (26, 17, 21)/64, (27, 19, 18)/64, (25, 17, 18)/60, where the m-system is not at the replicated configuration;
the exact D law gives 0.996–1.002 (20 pairs) and 0.948–0.993 at the saddles — round 10's confirmation 6 for the new runs. `report()` has the
ingredients (`weighted_stationary_atoms_prior`, `transport_plan`). *Fix:* print the exact law and quote "5–16%".

**D3d-12 (cosmetic; (iv)'s constants).** The two formal bounds are vacuous on the instance: R_Θ = max{R₀, L_ϑ^G(F_ρ + K_Φ)/σ₀} + max_j‖ϑ⁰_j −
ϑ⁰_1‖ = 7.26 ∨ 2858 + 7.2 = 2865, so 4(K₀ − 1)R_Θ²/m = 4.1·10⁶ at m = 16 against the measured W²(μ*₁₆, ν) = 1.03–1.37, and 2(K₀ − 1)K_Φ/m = 1.00
against ‖f̂*₁₆ − f‖_{L²(ρ)} = 0.0282·0.648 = 0.018 (factor 55). The measured constants are m·W² ≈ 20 = 0.51Δ_ν² (the bridge particles' plan
shares at 0.28–0.58 Δ_ν) and m·‖f̂*_m − f‖/‖f‖ ≈ 0.55. Worth one sentence, as for C₁ and m₀.

**D3d-13 (cosmetic; tests).** Test 2's `|g − FD| < 10⁻⁶·max(1, |g|)` has slack 1.4·10⁴ (actual error 7.4·10⁻¹¹ at its m = 10
configuration); test 1's `|W2 − Σγ_ijC_ij| < 10⁻¹²` is an identity of `plan` (which returns that sum); `CAB` in test 5 is unused. The
round-off assertions (`drift < 10⁻¹¹` against 1.3·10⁻¹⁷, `W2_empirical < 10⁻¹²`) are fine as zero-checks.

**D3d-14 (cosmetic; table 2's "predicted (1/m − 1/M)tr(FC̃Fᵀ)/N" column).** The w̄-law 3.99·10⁻³, 9.07·10⁻⁴, 1.36·10⁻⁴ is 1.6–2.3 times
the seed means (ratios 0.63, 0.58, 0.43) and 0.02–2.0 times the individual seeds; the exact function-level law from the realised weights,
‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖², reproduces the measured h² with ratios 0.999–1.009 in 19 of 24 pairs, 0.969–0.990 in two, and 0.818–0.906 at
the saddle weights (D3c-13's point, again).

**D3d-15 (cosmetic; (ii)'s mean-field form, the reference as stand-in).** For a disconnected supp μ the dual weights ψ, hence the Laguerre
cells, are not unique (any boundary in the gaps between clusters gives the masses q), so "N_ϱ′(supp μ_s) meets no Laguerre boundary" should
read "for some choice of the dual weights"; the map T_μ on supp μ is unique. And since the mixture start has w = q in the mean field, the
mean-field flow has no transit while the M = 510/512 references move 3–21 particles: during the transit the reference is a poor stand-in,
and the measured h² and W² against it carry the reference's own finite-M transit and, on the incommensurate sizes, its two bridge
particles (D3d-5).

## Independent confirmations (what was recomputed and agrees)

1. **Tables 1–4 and the prose** from `ot_leak_raw.npz`, independently of `report()`: all entries of the four cluster-mass tables (mean
   |w(0) − q|₁, ½|w^m(0) − w^M(0)|₁, max_s|w(s) − w(0)|_∞ = 0 in 32 systems, |w(end) − q|_∞ = (2/3)/m exactly on the incommensurate sizes,
   s_eq and s_tr per seed, transit at s = 0 and at the end, split counts), the fit errors with their standard errors, the last-window
   distances (h², W², D under both leaks, the predicted column, the cross-cluster column), the ratios OT/nearest (1.95·10⁻⁴, 9.4·10⁻⁵,
   1.27·10⁻⁴; 1.54·10⁻³, 8.6·10⁻⁴, 8.5·10⁻⁴), last/previous (0.997, 0.973, 0.996 / 0.018, 0.035, 0.028; 0.995, 0.995, 1.002 / 0.118,
   0.120, 0.121; incommensurate 0.819, 0.898, 0.971 / 0.973, 1.002, 0.977), the slopes with their s.e. (OT h² −1.31 ± 0.60 / −2.07 ± 0.05,
   W² −0.77 ± 0.17 / −1.21 ± 0.04, D −0.71 ± 0.21 / −0.72 ± 0.19; nearest h² −1.36 ± 0.37 / −1.59 ± 0.27, W² −0.67 ± 0.13 / −0.79 ± 0.12, D
   −1.40 ± 0.26 / −1.60 ± 0.26), the ten window means per m, the 30 bridge entries (dist/Δ_ν 0.276–0.583, distance to target 0.025–0.595,
   plan rows ∈ {(1/3, 1/3, 1/3), (2/3, 1/3, 0) and permutations}), (K₀ − 1)Δ_ν²/m = 4.90, 1.23, 0.31 with ratios 0.236, 0.233, 0.134, and the
   design table (fit, κ_R at the atoms, κ_R over N_{0.1}, margin, splitting margin, ε < 3.131, atoms) to all printed digits. s_eq over the
   16 commensurate systems: 93–164, median 103.5; ln 2/σ₀ = 97.98; Δ_ν = 6.263 (atom pairs 7.205, 6.686, 6.263); 600σ₀ = 4.245.
2. **Design table from scratch** (`mass_design_row` re-executed): identical to the stored rows to all digits (atoms to 0.0); the δ = 0.02
   splitting margin +0.00279 equals Proposition D‴.9's "excess on atom 1 alone" value — the same computation (w = w̄ + δ(1.5e₁ − 0.5)).
3. **Regeneration**: `report()` from the raw data, output redirected, byte-identical (md5 `b537eab2530b54d96384357ac7c12986`).
4. **The plan.** LSA on the replicated atoms against the transportation LP: equal cost to ≤ 5.3·10⁻¹⁵ on 9 commensurate configurations
   (m = 15, 60, 240; jitter 0.02, 0.4 and cold), all integral; on 15 incommensurate configurations (m = 16, 64, 256, 7, 23) and the 16 stored
   final plans the LP solution is a vertex (support a forest: edges = m + K₀ − components, one component), with 2 split rows and entries in
   {1, 2, 3}/(3m); nsplit ≥ 1 on every incommensurate record (min over the runs 1 or 2), ≤ 2 always.
5. **Danskin.** grad_J against central differences of J: relative error 3.4·10⁻¹⁰–2.7·10⁻⁹ at eight configurations (q uniform with m = 15,
   16; q = (0.5, 0.3, 0.2) with m = 10, 7; jitter 0.02 and 0.4); drift + G·grad_J = 0 exactly.
6. **The Hessian at Θ*₁₅(q)**: FD Hessians of J₁₅^{OT} and of the nearest-atom J₁₅^ν agree to 1.7·10⁻¹⁸ (scale 0.16); whitened eigenvalues
   σ₀ = 0.007074 with multiplicity 36 and 0.00738, 0.00904, 0.01019, 0.02222, 0.02818, 0.03955, 0.09569, 0.14127, 0.20229 (round 10's K-atom
   Hessian); |drift| = 1.3·10⁻¹⁷; the plan is diagonal.
7. **`W2_empirical`**: replication route against the LP at 16|512 (diff 4.8·10⁻⁹, the LP's tolerance) and 15|510 (4·10⁻⁹); the LP route at
   60|510 against the assignment problem on the lcm replication (1020 × 1020): 5.0·10⁻⁹; symmetric to 6.4·10⁻¹⁰.
8. **Reproducibility**: the clouds regenerated from the seeds give the stored first-record weights in all eight (seed, set) pairs; the
   four incommensurate OT runs (m = 16, 64 and the M = 512 reference) re-integrated from their seeds reproduce the stored records at
   s = 599.02 to all printed digits (D, h², w, fit error, W²(μ, ν)) in all four seeds — the systems are autonomous and deterministic, as the
   brief says.
8a. **The stationary state at incommensurate m is set by the plan's split pattern, not by the seed** ((iv)'s "deterministic quantisation
   constant"): at s = 2400 seeds 1 and 2 (m = 16; bridge rows (0, 2/3, 1/3) and (1/3, 2/3, 0)) have identical plateaus — h² 3.1321·10⁻⁴,
   W² 1.3394, fit 0.028266, bridge particles at 0.254 and 0.320 Δ_ν — and seeds 0 and 3 (one (1/3, 1/3, 1/3) particle at 0.551 Δ_ν, the
   centroid) the identical fit 0.025283, their h² (2.455 against 2.554·10⁻⁴) differing only through the reference, whose own stationary state
   is the same for seeds 0–2 (fit 1.2260·10⁻³, bridge rows (0, 2/3, 1/3) + (2/3, 0, 1/3)) and different for seed 3 (6.772·10⁻⁴, rows
   (0, 1/3, 2/3) + (2/3, 1/3, 0)); at m = 64 seeds 0, 1 (fit 0.0092864) and 2, 3 (0.0090041) coincide likewise.
9. **Stationarity of the stored finals** (max whitened |drift|): nearest-atom 7.9·10⁻⁶–1.0·10⁻³ (32 systems); OT commensurate
   5.7·10⁻⁴–1.1·10⁻³ (the migrated particles at Δe^{−σ₀·600}); OT incommensurate 5.9·10⁻⁴–1.9·10⁻³, the bridge particles at 1.6–1.9·10⁻³
   for m = 16 and 2.1·10⁻⁴–1.7·10⁻³ for m ≥ 64.
10. **The exact laws of the nearest-atom runs** (realised weights, weighted stationary points by continuation): W² ratios measured/exact
    0.993–1.001 (17 pairs), 0.950–0.974 (3), 0.846–0.938 (saddle weights); D ratios 0.996–1.002 (20), 0.948–0.993 (saddles); function level
    0.999–1.009 (19), 0.969–0.990 (2), 0.818–0.906 (saddles).
11. **The D(OT) plateau**: end/[(1/m)Σ_{mismatched}Δ²_{kk′}] = 0.970–0.979 in the twelve commensurate pairs; 0.47–1.14 in the incommensurate
    ones (bridge particles and the reference's own bridge particles intervene).
12. **The transit**: particles assigned away from their nearest atom at the first record 1–21 in all 32 OT systems (= m·½|w(0) − q|₁ exactly
    for the commensurate sizes); migrated particles at the end at 0.064–0.157 from their atom (Δ_{kk′}e^{−σ₀·600} = 0.090–0.103).
13. **The laws against ν** (D3d-5): W²(μ*_m, ν) slope −0.98 ± 0.03, m·W² = 19.3, 20.9, 20.0, 19.2; fit² slope −1.81 ± 0.05, m·fit = 0.45, 0.60,
    0.60, 0.54.
14. **Decay rates on (300, 600]** (commensurate OT): W² at 2.01σ₀, 2.00σ₀, 1.99σ₀ (factors 0.119, 0.120, 0.122 per 150, as the text's
    e^{−2σ₀·150} = 0.120); h² at 3.77σ₀, 2.98σ₀, 3.00σ₀ (0.018, 0.043, 0.042; text 0.02–0.035, e^{−4σ₀·150} = 0.014).
15. **P(w^m(0) = q)** = 0.0527, 0.0136, 0.00344, 0.00162 (m = 15, 60, 240, 510).
16. **Tests**: c2e1ac6 file 5 passed in 0.94 s; working tree (8 tests) 8 passed in 54 s.

## Not checked

The full experiment was not rerun (1.8 h); the m = 256 incommensurate runs and the commensurate runs were not continued beyond s = 600
(their last/previous ratios 0.97–1.00 and 0.02–0.12 say what they do); the mean-field transit of (v) was not simulated (the mixture start
has no mean-field transit);
the bibliographic details of [S40]–[S42] were not verified beyond recognising the results cited; Theorem C′'s proof and Theorem D‴'s Steps
0–9 were not re-refereed beyond locating where the leak enters; the D‴.14 material of fb23968 and the working-tree edits were outside the
brief; the 60|510 and 240|510 W² values of the tables were taken from the raw data (one 60|510 final pair was recomputed by two routes).
