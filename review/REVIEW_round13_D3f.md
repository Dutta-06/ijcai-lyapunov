# Referee report D3f (round 13) — PARADIGM §5d: Remark D‴.16 and its *Numerics*; the paragraph *Numerics in two dimensions*; Proposition D‴.17 with Status; the §0 clauses, the §7 rows D‴.16, D‴.14 (n = 2), D‴.17, the §8 bullets, `results/RESULTS_PARADIGM.md` (D‴.16, D‴.14 in two dimensions); `experiments/exp_dither_ot.py` → `results/chaos/dither_ot.md`; `experiments/exp_ot_2d.py` → `results/chaos/ot_2d.md`; the last test of `tests/test_ot_leak.py`; [S43], [S44] and the D‴.17 audit bullet of `LITERATURE_TRANSPORT.md`; verification of the fixes of rounds 11 and 12

Branch `claude/festive-bardeen-eopej2`, reviewed at ed2ee81 (HEAD when the review started). While it ran, three commits from other sessions
landed (e4fbed6 manuscript, 5ff2e0a Theorem D‴.18, cffb7be Proposition D‴.20 part 1); `git diff ed2ee81 HEAD` touches none of the reviewed
passages, scripts, raw files, tests, LITERATURE entries or review files (only the header lists and one unrelated §8 line), so every finding
applies to the current HEAD. Line numbers below are those at cffb7be (PARADIGM.md: ed2ee81's + 1 in §5d, + 101 in §7, + 106 in §8;
RESULTS_PARADIGM.md: + 14 from its D‴.16 section on).
No repository file was modified: the `--from-raw` regenerations ran on copies of the two scripts in a scratch tree (`regen/`, `src` symlinked),
every other script imports the repository read-only (`PYTHONDONTWRITEBYTECODE=1`, pytest with `-p no:cacheprovider`); the first rerun
imported `exp_ot_2d.py` before that flag was set and left the git-ignored `experiments/__pycache__/exp_ot_2d.cpython-313.pyc`, removed at the end.
Scratch scripts and their `.out` files in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3f/`
(the folder also holds files of another session — `fr-*.png`, `pg-*.png`, `satellites.py`, `stretch*`, `worddiff.py`, … — not part of this review):
`recompute_dither.py` (every number of `dither_ot.md` and the D‴.16 prose from `dither_ot_raw.npz`; count-vector transitions, per-particle
classification of the OT excursions, IATs, Poisson intervals, the whole count distribution against the laws, the OT off-balance fractions),
`dither_batch_means.py` (batch-means errors, dwell times, the restricted energy's curvature against the multinomial entropy),
`gibbs_variants.py` (energy-only / energy × multinomial / Laplace laws, TV distances, seed-to-seed TV), `reproduce_dither_prefix.py` (the
first 3000 slow units of two runs from their seeds), `ot_dither_plan_events.py`, `ot_dither_plan_changes_detail.py` (OT runs with the
plan recorded at every refresh),
`ridge_barrier.py`, `swap_barrier.py`, `swap_path.py` (the ridge minima, the plan's tie-set minimum and the exchange path on the instance),
`split_minima.py` (the restricted minima at saddle count vectors), `recompute_ot2d.py` (every number of `ot_2d.md` and of the 2-d paragraph),
`rerun_ot2d_seed0.py` (seed 0's OT run reproduced from its seed and continued to t = 5000), `check_test_D317.py` (the D‴.17 test instrumented,
the tie ranked, a control off the tie), `check_pure_leak_plan.py`, `check_pure_leak_plan_instance.py`, `pure_leak_counterexample.py`
(the pure-leak plan stability for split plans), `check_plan_switches.py` (cold starts at σ₀ and 10σ₀), `other_critical_points.py`,
`second_collapsed_point.py` (D‴.17(iii)), `check_D3d1.py`, `check_D3d4.py`, `check_D3d5.py`, `check_D3e1.py`, `check_D3e7.py`, `check_D3e12_and_misc.py`,
`misc_numbers.py`, `tests.out` (9 passed, 54.6 s), `from_raw_*.out` under `regen/`.

## Verdicts

| object | verdict |
|---|---|
| **Remark D‴.16**, the nearest-atom half | **holds (formal) with fixes.** The Gibbs measures of the min-of-quadratics potentials J^ν_m, J^{OT}_m are well defined (continuous, bounded below, quadratic growth from the leak; the a.e. gradient is the drift; non-degenerate noise), and ΔV ≥ σ₀Δ_ν²/8 is rigorous for a hop out of the balanced state (J = 0 there, L̄ ≥ 0, the other leaks ≥ 0): the computed ridge minima are 0.040 (1\|2 ridge, leak part 0.035) to 0.074. From an unbalanced state the barrier toward balance is lower by the start's energy (the same ridge set). "min_Θ J^ν_m(Θ; w) = m𝒢*_w" needs w ∈ W_adm, and the Gaussian variance τ₀/(m∇²𝒢*) omits the multinomial entropy (D3f-2). |
| **Remark D‴.16**, the OT half | **fails.** A nearest-atom count changes under the OT leak at the nearest-atom ridge barrier — the plan keeps a particle on a ridge assigned to its own atom, J^{OT} = J^ν = 0.0402 at the lowest ridge point — not at σ₀Δ_ν²/2 = 0.139; the plan switch (an exchange) costs at most 0.083; the measured off-balance probabilities are 10²–10⁴ times "exp(−σ₀Δ_ν²/(2τ₀))"; and "a hop is answered by a swap of assignments and the transit of one particle back" is false: the crossings are reversible excursions of one particle, pulled back by its own atom (D3f-1, **major**). |
| *Numerics* of D‴.16 (`dither_ot.md`, prose, §0 clause, §7 row, RESULTS) | **numbers hold, readings fail in part.** Every number reproduces from the raw data; the md regenerates byte-identically; the first 3000 slow units of two runs reproduce the stored labels exactly. The OT reading is false (D3f-1); the "leading-order Gibbs law" is not an order, agrees with the balanced fraction within error bars of ≈ 0.1 at τ₀ = 0.01 and misses it by 2.1 SE at 0.02, and mis-predicts the count distribution by factors 2–5 (D3f-2); the barrier 0.039 rests on one event, and the "hop rates" count record-level label changes, half of them immediately reversed (D3f-3). "frozen for 350/σ₀" is fair (0.43 hops expected from the Arrhenius extrapolation, one seen). |
| `experiments/exp_dither_ot.py` | **holds.** Euler–Maruyama of dϑ_i = −G∇_{ϑ_i}J ds + √(2τ₀)G^{1/2}dW_i (the drift is `loop.drift` = −G∇J^ν or −G∇J^{OT} with the leak target lagged at most 10 steps = 0.5 slow units; the noise √(2τ₀ds)G^{1/2}ξ) — a consistent discretisation (lag ∝ ds); the Laplace table is the right Laplace expansion (multinomial(n) × e^{−m𝒢*/τ₀} × the separation blocks with multiplicity n_k − 1 × the whitened weighted collective block; no further label-symmetry factor). Deterministic from the seeds. Cosmetic: np.int64 reprs in the md; the docstring repeats the false mechanism (D3f-1). |
| *Numerics in two dimensions* (`ot_2d.md`, prose, §0 clause, §7 row, RESULTS) | **holds with cosmetic fixes.** All 12 rows, the equalisation times and every range reproduce; regeneration byte-identical; seed 0's OT run reproduced from its seed matches the stored configurations at t = 2500 bit for bit; the clouds are one array for both leaks (nested sizes); the exact representation is inside the box (|4A_k| ≤ 0.8a_max) and reproduces f to 3.3·10⁻¹⁶; the control and the reference are evaluated at the RK4 stages; the OT tracking error is the fit error over k_e (e_rms·k_e/‖r‖_rms = 1.05–1.07). "Three hundred to a thousand" is 285–480; the OT distance at 5/σ₀ is a snapshot of a steady decay (D3f-8). |
| **Proposition D‴.17 (i)** | **holds.** J^{OT}_m = min over the vertices of the smooth J^γ_m (−G∇_{ϑ_i}J^γ = G⟨r∇Φ̃⟩ − σ₀(ϑ_i − T̄_i) since Σ_jγ_ij = 1/m); quadratic minus a maximum of linear functions, so semiconcave + smooth; ∂_C of a minimum of C¹ functions is the convex hull of the active gradients (−J is a regular maximum); vertex-selection velocities lie in −G_m∂J; the one-sided derivative along an active drift is min_{γ′ active}⟨∇J^{γ′}, −G∇J^γ⟩ ≤ −‖∇J^γ‖²_G, so J is non-increasing without any definability. |
| **D‴.17 (ii)** | **holds.** Definable: for the finite quadrature in ℝ_exp on all of ℝ^{mp}; for the period average restricted-analytic on the ball (ℝ_an suffices); the generalised (optimal-face) velocities are −G_m∇J^γ with γ in the optimal face, hence in −G_m∂J, and Filippov's set is inside Clarke's; the chain rule (DDKL; or upper regularity) gives dJ/ds = −‖Θ′‖²_{G_m⁻¹}; KL + boundedness (Step 0 with conv supp ν) gives finite length and convergence. The metric G_m needs the linear change Y = G_m^{−1/2}Θ, unstated (D3f-9). |
| **D‴.17 (iii)** | **holds as stated** (a unique optimal plan is a vertex and stays optimal nearby, so the limit is a critical point of J^γ); **the gloss classifying the limits is false** — a second collapsed critical point with a unique diagonal plan exists on the instance, a saddle although q ∈ W_adm (D3f-5). |
| Status of D‴.17 | **holds except one sentence:** "under the pure leak a plan never switches once optimal" is proved for integral plans only and false for split plans; the empirical extension ("the fit force did not overturn it … the switches of Theorem C′'s numerics belong to priors with many atoms") holds at 10σ₀ only — at the experiments' σ₀ cold starts switch the plan (D3f-4). |
| the last test of `tests/test_ot_leak.py` | **passes and tests what it says.** The tie is genuine (the two cheapest (2,2,2)-assignments are π_a, π_b at 7.520474, the third at 7.632397) and the assertion discriminates (off the tie by 10⁻³ it fails for the non-optimal plan); the run's monotonicity check covers a smooth flow only — the plan never switches (D3f-4). |
| [S43], [S44], the audit bullet | **adequate** (BDL 2007, BDLS 2007, DDKL 2020, Kurdyka 1998, vdDMM 1994 are the right sources for KL, the nonsmooth KL for definable functions, the chain rule and the o-minimality; the authors flag the details as from memory); "(smooth: Łojasiewicz [S44])" cites Kurdyka/vdDMM for Łojasiewicz's inequality (D3f-9). |
| fixes of rounds 11 and 12 | **in and correct, with one fix missing and leftovers.** D3d-1 (0.053, 0.014, 0.0034, 0.0016 ≈ 0.83/m), D3d-4, D3d-5 (m·W² = 19.4–21.0, −0.98 ± 0.03; m·fit 0.45–0.60, −1.81 ± 0.05), D3e-1 (margins, spreads, ratios), D3e-4 (the averaging proof is correct) and D3e-7 (the OT tracking error is the fit error over k_e: ratio 0.94–0.97) verified; D3e-12 is logged as applied but is not in the repository (D3f-6); leftovers of D3e-1, D3e-9 and D3d-2 (D3f-7). |

## Confirmed issues

**D3f-1 (major; Remark D‴.16 lines 2159–2163, its *Numerics* lines 2172–2174, `dither_ot.md` Reading, `exp_dither_ot.py` docstring lines
13–14, RESULTS line 460, §8 bullet lines 4883–4884).** The OT half of the Remark — "J_m^{OT} charges a particle sitting in one cluster while
the plan assigns it to another atom with … ≈ σ₀Δ_ν²/2 — a barrier of 0.14 … — so … the masses by nearest atom are q up to events of
probability of order exp(−σ₀Δ_ν²/(2τ₀)): a hop across a ridge is answered by a swap of assignments and the transit of one particle back" — is
false on both counts, and the Remark's own numerics show it.
(a) *The barrier.* A count by nearest atom changes as soon as one particle crosses a Voronoi ridge, and the plan does not follow it: with
integral plans the exchange of the assignments of a ∈ class j and b ∈ class k lowers the cost iff ⟨ϑ_a − ϑ_b, ϑ⁰_k − ϑ⁰_j⟩_{G⁻¹} > 0, i.e. only
once a has passed b, so on the ridge the crossing particle is still assigned to its own atom and J^{OT} = J^ν there — the barrier of a count
change is the nearest-atom ridge barrier. On the instance (`ridge_barrier.py`, `swap_barrier.py`): the lowest point of the 1\|2 ridge (one
class-1 particle on the ridge, the 14 others relaxed; |∇| ≤ 7·10⁻⁹) has J^{OT} = J^ν = 0.0402 with the diagonal plan, and nudged across it the
counts are (5, 4, 6); the other ridges give 0.0445, 0.0523, 0.0678, 0.0743 — against σ₀Δ_ν²/2 = 0.1387 (0.1581, 0.1836 for the other pairs).
The plan does switch on its tie set: the lowest tie point for an exchange between classes 1 and 2 is 0.0831 (0.1079 for 0↔1, 0.1244 for 0↔2;
there the two particles are at equal distances 6.36, 3.66, 2.62 from the three atoms), and descending from it on either side reaches the balanced minimum (J = 5.7·10⁻⁸) with the particles
exchanged or not, the maximum along both descents 0.0831 (`swap_path.py`): an exchange costs at most 0.083, not 0.139.
(b) *The data* (`recompute_dither.py`). After the initial transit the OT counts are off balance in a fraction 0 (τ₀ = 0.005, both seeds),
6.3·10⁻³ (τ₀ = 0.01; 7.5·10⁻³, 5.2·10⁻³) and 0.117 (τ₀ = 0.02; 0.137, 0.098) of the records, against exp(−σ₀Δ_ν²/(2τ₀)) = 8.9·10⁻¹³, 9.4·10⁻⁷,
9.7·10⁻⁴ — 6700 and 121 times the stated order; the effective barrier of the off-balance probability between τ₀ = 0.01 and 0.02 is 0.058,
inside the ridge range 0.040–0.074, not 0.139. The text's own "the label changes are as frequent as the nearest-atom leak's hops" (records
131/111, 115/93, 364/281, 344/275 = 1.18–1.30) is what equal barriers predict; with 0.139 against 0.048 they would be e⁻⁹ as frequent at τ₀ = 0.01.
(c) *The mechanism.* Classifying every off-balance excursion by the per-particle labels of the balanced records before and after it: at
τ₀ = 0.01 all 121 (65 + 56) end with the same particle back in its cluster and none with an exchange; at τ₀ = 0.02, 298 of 336 (89%) are
returns and 38 exchanges; at τ₀ = 0.005 the one label change of each run is the initial transit (records 20 and 16; no excursion in 50 000 slow
units). Median durations 1 and 2 records (10–20 slow units). Re-running seed 0 with the plan recorded at every refresh (0.5 slow units; 20 000
slow units after the transit; same trajectory — off-balance 0.164 against the stored records' 0.161 at τ₀ = 0.02): at τ₀ = 0.01 the plan never
changes (0 changes, 21 excursions); at τ₀ = 0.02 it changes 110 times, all inside excursions and 77 of them undone at the next change, and
only 9 of the 314 excursions end with the plan changed (`ot_dither_plan_events.py`, `ot_dither_plan_changes_detail.py`).
"88% … the excursions being the transits that follow a swap of assignments" (md), "(swaps followed by transits)" (RESULTS), "a hop is answered
by a swap of assignments and the transit of one particle" (Numerics), the docstring's "a hop is answered by a swap of assignments and a
transit", and the §8 bullet's "keeps the plan's masses at q up to a barrier σ₀Δ_ν²/2 per misplaced particle" (the plan's masses are q by
construction) all fall with it. The 99.6% at τ₀ = 0.005 is 1 − transit/horizon (200 and 160 slow units of 50 000).
*Fix:* "Under the OT leak the plan keeps a particle that crosses a Voronoi ridge assigned to its own atom (an exchange needs two particles to
pass each other: at most 0.083 on the instance, at least σ₀Δ²/4 in leak), so the counts by nearest atom leave q at the nearest-atom leak's
ridge barrier (≥ σ₀Δ_ν²/8; ≈ 0.040 on the instance) but only for a reversible excursion — beyond the ridge the leak toward the particle's own
atom keeps rising and pulls it back — of stationary probability of order exp(−ΔV_ridge/τ₀); the plan's masses are q always, the masses by
nearest atom q except during these excursions"; in the numerics "balanced 99.6% (the initial transit only), 99.2% and 88% of the time: the
excursions are returns of the crossing particle (all 121 at τ₀ = 0.01, 89% at 0.02, the rest exchanges)"; RESULTS, §8 and the docstring accordingly.

**D3f-2 (minor; Numerics lines 2169–2171, `dither_ot.md` τ₀ sections and Reading, §7 row line 4732, RESULTS line 458; Remark lines
2155–2157).** The "leading-order Gibbs law multinomial(n)·exp(−m𝒢*_n/τ₀)" and the agreement read off it.
(a) *Not an order.* The multinomial and the Laplace determinant factors are both τ₀-independent prefactors: as τ₀ → 0 the leading order is
e^{−m𝒢*_n/τ₀} alone, and the next order contains both. On the instance the multinomial factor of each of the six neighbours of (5, 5, 5) is
0.83 relative to it, while the dropped determinant factors are e^{−0.99} … e^{+1.90} for the four admissible neighbours and undefined at
(6, 4, 5), (6, 5, 4) (saddles; 37 of the 91 count vectors are); the energy-only law fits the data as well (TV 0.202 against 0.203 at τ₀ = 0.01,
0.140 against 0.160 at 0.02; `gibbs_variants.py`).
(b) *The agreement "0.54, 0.25 against 0.58, 0.33".* The balanced indicator's integrated autocorrelation time at τ₀ = 0.01 is 2306 and 744
records (23 000 and 7 400 slow units; N_eff = 8.7 and 26.9; batch means give per-seed SEs 0.08–0.13 and 0.06–0.10, with 20-batch values from 0.00
to 1.00): the seed mean 0.537 carries an SE ≈ 0.1 and 0.576 is 0.4 SE away — the comparison cannot discriminate; at τ₀ = 0.02 (N_eff = 59 per
seed, SE 0.040) the law is 2.1 SE above the data (0.332 against 0.249).
(c) *The distribution.* TV(data, law) = 0.20 at τ₀ = 0.01 against a seed-to-seed TV of 0.074: the law under-weights the saddle count vectors
(6, 4, 5) and (6, 5, 4) by 3.2× and 2.7× and over-weights (5, 6, 4), (4, 6, 5), (4, 5, 6) (data/law 0.49, 0.21, 0.19); at τ₀ = 0.02 (7, 4, 4) is
under-weighted 4.1× (TV 0.16, but the seeds differ by 0.18 there — the horizon does not resolve the tail). The split below a saddle gains
little energy (0.0004, 0.0009, 0.0034 for (6, 4, 5), (6, 5, 4), (7, 4, 4): factors 1.04, 1.09, 1.40 at τ₀ = 0.01; `split_minima.py`), so the
excess is the entropy of the soft and unstable separation modes — the very factors the "leading order" drops.
(d) In the Remark, "w ↦ min_Θ J_m^ν(Θ; w) = m𝒢*_w" holds only for w ∈ W_adm (at a saddle weight the split lies lower, by the amounts in (c)), and
"Gaussian fluctuations of variance τ₀/(m∇²_w𝒢*_w)" omits the multinomial entropy: along the neighbours' directions the entropy precision
3m = 45 is 7–16% (τ₀ = 0.01) and 14–32% (τ₀ = 0.02) of the energy precision m·uᵀ∇²𝒢*u/τ₀ = 285–666 and 142–333, and the cubic part of 𝒢* is
11–20% of the quadratic at one count. *Fix:* "the energy–multinomial law (the Laplace factors, which matter as much, are omitted)", report the
TV and the SEs with the comparison, say that the soft/unstable separation modes carry the excess at the saddle count vectors; "variance
τ₀/(m(∇²𝒢* + τ₀ diag(1/q)))", "for w ∈ W_adm".

**D3f-3 (minor; Numerics line 2167, `dither_ot.md` Reading, §7 row line 4732, RESULTS line 457).** "Arrhenius barriers 0.039–0.048 against
the leak's ridge height σ₀Δ_ν²/8 = 0.035". The 0.039 rests on one event (τ₀ = 0.005: 0 + 1 hops in 100 000 slow units; 95% Poisson interval
of the rate [2.5·10⁻⁷, 5.6·10⁻⁵], hence a barrier in [0.021, 0.078]); 0.048 rests on 204 and 556 record events ([0.043, 0.052]). The "rates"
count records (every 10 slow units) at which some label changed, and 55/111, 41/93, 113/281, 146/275 of them are undone by the same particle
at the next change — recrossings; they are sampling-dependent crossing rates, not hop rates. The comparison with σ₀Δ_ν²/8 compares with the
leak part of the lowest ridge only; the computed single-particle barriers out of (5, 5, 5) are 0.040–0.074 and 0.048 lies among them.
*Fix:* "0.048 (0.043–0.052 at 95%, from τ₀ = 0.01, 0.02; τ₀ = 0.005 has one event)", "record-level label changes, about half of them undone",
"against the ridge minima 0.040–0.074 (leak part 0.035)".

**D3f-4 (minor; D‴.17 Status lines 2205–2209; §7 row line 4734).** "Under the pure leak a plan never switches once optimal: … the cost
difference … is the convex combination of its value at Θ(0), nonnegative, and at the targets, positive for distinct atoms." The value at
the targets is C_{γ′}(T̄) − C_γ(T̄): positive for every γ′ ≠ γ when γ is integral (T̄_i an atom, C_γ(T̄) = 0), i.e. for commensurate m, but not
for a split plan, whose barycentric targets are not atoms — such a plan need not be optimal at its own targets, and then the pure leak
switches it. Counterexample (`pure_leak_counterexample.py`; planar, G = I, atoms (−0.3905, 1.6273), (−1.1755, 0.1601), (−2.1378, −0.0016),
q = 1/3, m = 4, Θ(0) = (1.3493, −0.355), (−0.944, 0.3473), (1.0502, 0.9955), (2.9587, 0.3138)): the vertex plan sends particles 1, 2, 4 to
atoms 2, 3, 1 and splits particle 3 in thirds; at its own targets it costs 0.0022 more than the optimum, and along the exact pure-leak segment
it stops being optimal at σ₀s = 1.818. Random incommensurate instances switch in 216 of 400 (K₀ = 3–5, random q, G; `check_pure_leak_plan.py`)
and in 62/300, 48/300, 58/300 planar ones with q = 1/3 at m = 4, 5, 7; commensurate ones in 0 of 200; on the instance's prior 0 of 400 clouds
(`check_pure_leak_plan_instance.py`). The empirical extension "the fit force did not overturn it either …; the switches of Theorem C′'s numerics
belong to priors with many atoms" is specific to the ten-times leak of the test: at the experiments' σ₀ cold starts on the instance switch the
plan in 2 of 6 runs at m = 15 and in all 6 at m = 16 (1–4 switches in 100 slow units), at 10σ₀ in none (`check_plan_switches.py`); and in Theorem C′'s
numerics every row of the plan is split (m = 6 particles share 48 atoms), which is the case the argument does not cover. The test's run never switches (0 in 300 steps, the plan's cost gap
growing 0.056 → 3.5; `check_test_D317.py`), so its monotonicity assertion checks a smooth gradient flow and the nonsmooth part of (i) is tested
at the static tie only. *Fix:* "for commensurate m (integral plans); a split plan can switch even under the pure leak", replace the
attribution to many atoms, and, if wanted, add a run at σ₀ with a switch to the test.

**D3f-5 (minor; D‴.17(iii) lines 2195–2199).** "for commensurate m the plan is integral and the limit is a critical point of J_m with every
particle leaking toward its assigned atom — the replicated weighted configuration Θ*_m(q) when the classes are collapsed, or a split
configuration of Proposition K's type when q lies outside W_adm —" reads as a classification of the unique-plan limits, which (iii) does not
prove: a collapsed limit is the replication of *any* critical point of 𝒢_q at which the diagonal plan is the unique optimal plan, and a split
can occur for q ∈ W_adm. On the instance, Newton from 300 random starts finds 14 critical points of 𝒢_q (`other_critical_points.py`); one is a
second local minimum of 𝒢_q (energy 0.107, collective eigenvalues positive; amplitudes 3.15, −3.18, −1.42) whose replication to m = 15 has the
diagonal plan as its unique optimal plan (cost gap 1.94) and zero drift (1.7·10⁻¹⁷): a Clarke-critical point of J^{OT}_15 with J = 1.607 (fit
error 0.37) other than Θ*_15(q), and a saddle — its class-0 separation block has the eigenvalue −0.0261 although q ∈ W_adm
(`second_collapsed_point.py`); a collapsed start near it stays collapsed and converges to it. "a critical point of J_m" should read J^π_m (J_m is
the single-prior energy). *Fix:* "a critical point of the frozen energy J^π_m — Θ*_m(q), the replication of another critical point of 𝒢_q, or a
split configuration —; a trajectory that enters the basin of Θ*_m(q) converges to it".

**D3f-6 (minor; `review/REVIEW_LOG.md` D3e-12; `tests/test_ot_leak.py` lines 159, 174, 191, 212; D‴.14 Status line 2063 "the atoms
moving").** The round-12 fix "motion measured against the start; the two lifts' displacements compared" is logged as applied but is not in the
repository (56beabb changed only the round-11 test items): test 1 still asserts `np.abs(Th15 − Th15[0]).max() > 1e-3` — 6.6247 at t = 0 and
6.6223 after the 300 steps, while the motion against the start is 2.3·10⁻³; test 2 still ends with `W2 > 1e-3` — 3.3687 at t = 0, 3.3686
after 200 steps, both lifts moving by 3.4·10⁻⁴ and their class-0 atoms 5.7·10⁻⁸ apart; the docstrings still cite "Corollary D‴.14(a)/(b)"
(`check_D3e12_and_misc.py`). *Fix:* apply the logged fix (compare with the starting configuration; compare the two lifts' displacements or
f̂ against their starts) and the docstrings, or correct the log.

**D3f-7 (cosmetic; leftovers of D3e-1, D3e-9, D3d-2).** "the contracted clouds" for every start (lines 2071–2072, §7 row line 4730, RESULTS
line 406, `ot_real_loop.md`, `exp_ot_real_loop.py` docstring, comments and report) although the three saddle clouds have within-class spread
0.04–0.11 (D3e-1); `ot_real_loop.md` §2's title "every size on the orbit of its own weights" (false for three of six) and "a ratio of
1.0003–1.0003 (1 pairs)"; RESULTS line 433 "identical across m to four digits" (D3e-9: 10⁻⁴ relative; PARADIGM was corrected);
`ot_leak.md`'s incommensurate line "D/plateau = 0.449–1.116 … D converges to the plateau from below" (a ratio above 1).

**D3f-8 (cosmetic; 2-d paragraph line 2097, D‴.14 numerics line 2079, `exp_ot_2d.py` docstring line 7, `ot_2d.md` header, `dither_ot.md`
last column).** "three hundred to a thousand times smaller": the per-system nearest/OT tracking ratios are 285–480 (the cross-range
116–1737); the 1-d "40 to 2000 times the OT leak's" is the cross-range (per system 126–548). The OT distance "5.0·10⁻⁸–1.6·10⁻⁷" at 5/σ₀ is a
snapshot of a steady decay: seed 0's run, reproduced from its seed and continued (`rerun_ot2d_seed0.py`), gives for m = 8 5.6·10⁻⁷, 1.6·10⁻⁷,
4.6·10⁻⁸, 1.4·10⁻⁸ at s = 90, 100, 110, 120 (×3.4–3.7 per 10 slow units) and 1.8·10⁻¹² at s = 200 = 10/σ₀ (m = 32: 7.3·10⁻¹³; still ×2.9 per
10 slow units, no floor), while the spread contracts exactly as e^{−σ₀s} (×0.606 per 10 slow units) — the number is the stopping time's and
the "five orders of magnitude" grows with the horizon. The docstring says "RK4 plant with the control held over the
step" (the code evaluates it at the stages); `ot_2d.md`'s header prints `np.float64(…)` reprs and `dither_ot.md`'s final counts `np.int64(…)`.
*Fix:* "290–480 times (per system)", "130–550 times", "at 5/σ₀, still falling ×3.5 per 10 slow units", the docstring, the formatting.

**D3f-9 (cosmetic; D‴.17 (ii) lines 2185–2192, Remark D‴.16 line 2159, the D‴.17 test docstring).** The KL and finite-length theorems of
[S43] are stated for ẋ ∈ −∂f(x); (𝔄_m^{OT}) is dΘ/ds ∈ −G_m∂J, which the change Y = G_m^{−1/2}Θ maps onto ẏ ∈ −∂(J∘G_m^{1/2})(y) — one line worth
writing. For the finite quadrature J^{OT}_m is definable in ℝ_exp on all of ℝ^{mp}, and the period average is restricted-analytic on the
ball (ℝ_an suffices). "(smooth: Łojasiewicz [S44])" attributes Łojasiewicz's gradient inequality to Kurdyka/vdDMM. The Remark's "Lipschitz,
piecewise-smooth potential" is locally Lipschitz (the leak grows quadratically). The test docstring's "(the Clarke-subgradient inclusion)"
tests the descent along the active drifts, not the inclusion; "the plan's masses never leave q" is the constraint, not a finding.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration.** `exp_dither_ot.py --from-raw` and `exp_ot_2d.py --from-raw` (scratch copies) rewrite `dither_ot.md` (md5 8ae62b7e…) and
   `ot_2d.md` (md5 48f71d7a…) byte-identically.
2. **`dither_ot.md`** (`recompute_dither.py`): hop records 0, 1, 1, 1 / 111, 93, 131, 115 / 281, 275, 364, 344; every time fraction, mean
   |w − q|₁ and final count; the leading-order and Laplace predictions (0.8616, 0.0755, …; 0.7198, 0.2434, …); the Reading's rates 1.00·10⁻⁵,
   5.10·10⁻⁴, 5.56·10⁻³, barriers 0.0478 and 0.0393, ridge estimate 0.0347, balanced fractions 0.000/0.9964/0.8616, 0.5365/0.9925/0.5759,
   0.2486/0.8778/0.3317; the prose's 350/σ₀ (50 000σ₀ = 353.7), 0.139 = 28τ₀, 14τ₀, 7τ₀, 19.8 × 0.007; §7, RESULTS and §0 agree.
3. **The runs are what the script says.** The first 3000 slow units of (nearest, τ₀ = 0.02, seed 0) and (OT, 0.02, 0) reproduce the stored
   300 label records exactly (`reproduce_dither_prefix.py`).
4. **The Laplace table and the saddles.** Splitting margins +0.00707 (5,5,5), +0.00395 (5,4,6), +0.00230 (5,6,4), −0.00247 (6,4,5), −0.00348
   (6,5,4) — the start, D3e-1's saddle weight —, +0.00321 (4,6,5), +0.00313 (4,5,6), −0.00563 (7,4,4), −0.00108 (6,3,6), matching the table's
   saddle flags; the collective block diag(w)^{−1/2}G^{1/2}∇²𝒢_wG^{1/2}diag(w)^{−1/2} is the class-mean block of G_m^{1/2}∇²J_mG_m^{1/2}.
5. **The nearest-atom barrier.** J = 6·10⁻³² at the balanced replicated configuration; the ridge minima 0.0402 (Δ₁₂ = Δ_ν = 6.263), 0.0445,
   0.0523, 0.0678, 0.0743, all above σ₀Δ²_{jk}/8 (0.0347, 0.0395, 0.0459): ΔV ≥ σ₀Δ_ν²/8 out of balance holds. The Arrhenius extrapolation from
   τ₀ = 0.01, 0.02 predicts 4.3·10⁻⁶ hops per slow unit at τ₀ = 0.005, 0.43 in the two runs (P(≥ 1) = 0.35); one was seen.
6. **`ot_2d.md`** (`recompute_ot2d.py`): the 12 rows (counts, spreads 2.8·10⁻⁴–1.8·10⁻², e_rms, fits, distances, W²), equalisation at
   13.44–13.84 slow units (ln 2/σ₀ = 13.86), every range of the prose, §7, RESULTS and §0; float32 storage harmless (f32 and f64 agree to the
   printed digits); amplitudes (3.2, 1.2), (−2, 2.4), (1.6, −2.8), (−2.4, −1.6) ≤ 0.8a_max, f reproduced to 3.3·10⁻¹⁶, Δ_ν = 6.137; per-pair
   distance ratios 1.1–2.1·10⁵; e_rms·k_e/‖r‖_rms = 1.053–1.074 (OT), 0.988–0.997 (nearest): the weight bias is a tracking error in 2-d too.
   Seed 0's OT run reproduced from its seed gives at t = 2500 the stored distances (1.5846·10⁻⁷, 6.6248·10⁻⁸), fits and spreads to every digit
   and the stored final configurations of all three sizes bit for bit (max |ΔΘ| = 0.0); continued to t = 5000 it shows no floor (D3f-8).
7. **The D‴.17 test** (`check_test_D317.py`): the initial plan reassigns 2 particles (counts (5, 7, 3)); the decrease per step ≥ 2.8·10⁻³
   against a tolerance of 10⁻⁹ relative; |dJ/ds + ‖drift‖²| = 2.3·10⁻¹⁰ against 10⁻⁶; J^{OT} = J^{π_a} = J^{π_b} = 1.8999566796 at the tie; each
   drift lowers J by ε‖d‖² (ratio 1.0000) and enters its own plan's region; 9 tests of the file pass in 54.6 s.
8. **D‴.17 (i)–(ii)**: no counterexample to the algebra (the gradients of J^γ, the semiconcave form, Clarke's subdifferential of a minimum,
   generalised and Filippov velocities in −G_m∂J, the chain rule, finite length under KL with bounded trajectories, definability).
9. **The OT ridge and exchange** (`swap_barrier.py`, `swap_path.py`): J^{OT} = J^ν = 0.0402 at the ridge minimiser with the diagonal plan;
   the tie-set minimum 0.0831 joins the balanced minimum to its exchanged copy (descents to 5.7·10⁻⁸ on both sides).
10. **Round 11.** D3d-1: P(balanced) = 0.05274, 0.01363, 0.00344, 0.00162 at m = 15, 60, 240, 510, m·P = 0.791–0.826 → 3^{3/2}/(2π) = 0.827;
    `mixture_start` gives 0.0512 and 0.0143 in 20 000 draws at m = 15, 60 — "0.053, 0.014, 0.0034, 0.0016, ≈ 0.83/m" holds. D3d-4: the migrated
    particles sit 0.064–0.157 from their new atoms at s = 600 against Δ_{kk′}e^{−600σ₀} = 0.090, 0.096, 0.103 (the "max whitened drift 5–8·10⁻⁴"
    is the md's largest whitened component; in the G⁻¹ norm 5.7·10⁻⁴–1.1·10⁻³). D3d-5: m·W²(μ_m, ν) = 19.40, 20.99, 20.10, 19.49, slope
    −0.984 ± 0.029 (−0.998 ± 0.025 with 512), squared-fit slope −1.808 ± 0.053 (−1.909 ± 0.069), m·fit = 0.452–0.597, m·W² = 0.510Δ_ν².
    The other round-11 fixes are in the text as logged (spot checks: the plateau ratio, the per-system counts, the transit mask, the
    tolerance 10⁻⁸, the LP check of the plan's cost, "−7% to +28%", the vacuous constants, Laguerre non-uniqueness).
11. **Round 12.** D3e-1: margins −0.00348, −0.00108, −0.00216 / +0.00561, +0.00440, +0.00544, spreads 0.11, 0.042, 0.10 against 5.5·10⁻⁴–1.3·10⁻³,
    distances to the lifts 4.1·10⁻⁸–3.3·10⁻⁶ against 9.6·10⁻¹²–2.5·10⁻¹¹, lift-to-lift ratios 0.9600, 0.9976, 0.9605 and 1.0003 (`check_D3e1.py`, from the raw records) — prose, §7 and
    RESULTS agree. D3e-4: the averaging proof is correct (ė = r − k_e e gives e_p = r_{Θ*}(x_m)/k_e + O(k_e⁻²); A_k = (1/k_e)[G∇²(rΦ̃) − σ₀I] +
    O(k_e⁻²) with period mean (1/k_e)(GR_k − σ₀I); the Magnus term is O(k_e⁻²); GR_k is similar to R^G_k, so Bauer–Fike applies) and needs no
    condition on m. D3e-7: the OT tracking errors are the fit error over k_e — e_rms·k_e/‖r‖_rms = 0.940–0.975 for all six m-systems
    (‖f‖_rms = 0.648; below 1 by the low-pass factor 1/√(k_e² + n²ω²)); the lift's 7.9·10⁻¹² is the scheme floor (its ‖r‖_rms/k_e = 6.1·10⁻¹⁰
    is the float32 storage of f̂). The other round-12 fixes are in (the proviso, W_adm in (iii), "locally constant near Γ", 0.6–3.2% and 0.2–4.1%,
    H_q, "at least ½‖w − q‖₁", the §8 proviso, "exactly σ₀/k_e", the display 9.989e-05), except D3e-12 (D3f-6) and the leftovers (D3f-7).

## Not checked

The full experiments were not rerun (deterministic; exp_dither_ot ≈ 45 min, exp_ot_2d ≈ 40 min): two 3000-slow-unit prefixes of the dither
runs, two instrumented OT dither runs of 20 000 slow units and seed 0's 2-d OT run (to t = 5000) were; the 2-d nearest-atom runs and seed 1 were
not. The ridge and tie-set minima are local minimisations (L-BFGS / SLSQP from symmetric starts): upper bounds on the lowest ridge and tie
points (the ridge's lower bound σ₀Δ²/8 is rigorous), which is all D3f-1 needs. The statements of [S43]/[S44] were checked against memory, not
against the papers (no web access was used); Theorem J′'s extension to the non-smooth leaks, the dithered real loop and incommensurate sizes
with dither are outside the material; whether the plan switches seen from cold starts at σ₀ change the limits was not followed.
