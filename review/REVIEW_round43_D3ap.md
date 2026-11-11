# Referee report D3ap (round 43): Theorem D‴.55 (the transit's certificate)

The theorem's empirical claims fail. The band "|m| ≤ 0.8q" (Status/numerics) and the claim that the amplitudes' range does not matter outside it both break, on both sides, at |m|/q of 1.06 and −1.62. Opposite-sign features are not "recovered closer together", which (iv) claims. "Every path ends at an exact fit" depends on a numerical parameter of the polish. The asymptotic statement in (iii) (q₀ and C exist) cannot be refuted numerically, because the constants are never given. It is a sketch rather than a proof, and its O(q) accounting leaves out a screening term.

## What I checked and how
- **Setup.** I copied the repository to `SCRATCH/repo` (commit 6666e70) and ran everything there with `venv/bin/python`, `OMP_NUM_THREADS=1` and `nice -n 10`. I ran one heavy process at a time, each under 14 minutes, with per-path alarms. I killed only my own PIDs. I wrote nothing to the repository.
  - The live repository has since moved to b1aa7ca, and `sliding_iss` gained a `max_entries` argument. The polish, exit and tolerance logic I tested look unchanged in that diff.
  - Scripts and JSONL outputs are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review43/work`.
- **Reading.** D55_draft, the D51/D53/D54 drafts, `exp_transit_certificate.py`, its tests, `rich_transit.sliding_iss` and `polish_block`, and `rich_transit_nd`.
- **Tests and numbers.** The tests pass (4 in about 3 s). Rebuilding the md from `_raw.npz` with the committed `report()` reproduces it byte for byte. I checked every number in the numerics paragraph against it.
- **Mutations.** 13 mutants. For each I ran the tests and regenerated the report from the raw data.
- **Reruns.**
  - Part A (400 targets) and Part C (120) with the polish off.
  - Part A with `t_fast=1e3`, and with the polish tolerances loosened to 1e-8.
  - Baselines and alternative certificates on all 520 targets.
  - Accounting check: the final and running duals against η.
- **New targeted experiments.**
  - Same-sign and opposite-sign pairs in 1-D (E1, 560 paths) and 2-D (E5, 128).
  - Pairs with amplitude ratios (E4); a lattice scan (E6); same-sign pairs at lattice points (E7).
  - Every ordering of the amplitudes at ranges 10³ and 10⁶ on targets outside the band (E2b, 137 paths).
  - K = 4 and 5 features, features near the edge, features exactly at lattice points (E3).
  - One feature at the edge; a feature exactly on a kernel.
- **Not tested separately:** unequal masses (out of the statement's scope, ρ = 1) and features at 2-D lattice vertices.
- **Literature:** from memory, plus three web checks (cited below).

## Verdicts

**(i) One feature.**
- Correct in the continuum (Lemma D‴.53a).
- False on samples, which the theorem explicitly allows:
  - c* = 0.955–0.98: the nearest kernel (1.0) is unique and the continuum margin is +0.017 to +0.099. The samples' certificate gives m = −0.18 to −0.24, and the path takes the kernel at 0.9. The same happens at c* = −0.96.
  - Near the edge ⟨φ_c∂_cφ_c⟩ ≠ 0 and Z is not constant, so the one-feature formulas do not hold there.
- "ẇ has A's sign at w = 0" is wrong. On the critical curve ẇ = 0 at w = 0, which is an equilibrium; the sign is kept by invariance, and the code needs a seed of 1e-7 to leave it.

**(ii) The accounting.**
- The order of entries is supported: as predicted on 253 of 262 (A) and 89 of 109 (B) exact-with-E paths. Every violation is a near-tie between consecutive keys (ratio ≤ 1.12).
- The per-phase cancellation of the amplitude (time ∝ 1/|A|) is real.
- But "duals at the end = η + O(q), uniform in the amplitudes" is not borne out:
  - On the committed paths the largest |u_j − η_j| is 0.40 = 0.947q.
  - With large amplitude ratios it grows. For f* = φ(·;0) + Aφ(·;0.3), with both features on starting points: at A ≤ 1/50 the kernels at 0.1 and 0.2 end at 0.975 and 0.9993, against η = 1.34. That is a deviation of 1.6–1.73q. At D = 0.4 the deviation is 1.55q.
- **Mechanism (screening).** The early-fitted large block absorbs part of the later feature's residual. The kernels next to that block then collect O(G(Δ − h − R)/G(R)) less than the certificate assigns. This is not O(G(Δ − 2R)/G(R)), and the excess over q grows like Δ when R → 0.
- The sketch's "block's error O(|A_m′|G(Δ − 2R))" leaves this term out.

**(iii) Sufficiency.**
- Plausible as an asymptotic statement, but unproved. q₀ and C are never given.
- The trapping region T has radius 2R, so it degenerates at R = 0 (all features on lattice points). Its forward invariance then needs C ∝ σ/R.
- The step "every j ∉ E stays below |η| + O(q)" is asserted, not shown.
- Empirically the band constant must be at least 1.057:
  - f* = φ(·;−0.1) − Aφ(·;0.2), both features on starting points: m = 0.2216, q = 0.2096.
  - The path is exact with E at A = 1, but not at A = 0.5, 1/3 or 10⁻³, nor with the roles reversed. The kernel at 0.3 slides onto the feature sitting exactly on the kernel at 0.2.
  - Confirmed with the polish off. An offset version (m/q = 1.047) took 527 s without polish and gave the same result.
- Equal-amplitude opposite-sign pairs at 2.5σ give five failures outside the 0.8q band, with m/q from 0.80 to 0.96.

**(iii) Necessity.**
- The contradiction argument is logically sound only given the unproved sufficiency machinery.
- It fails empirically:
  - The same-sign lattice pair above has m = −1.62q, yet it is exact with E for amplitude ratios ≥ 50 (and 'other' at 10–30).
  - A committed target with m = −0.90q ('other' as recorded) becomes exact with E on 9 of 11 completed reorderings at ranges 10³ and 10⁶.
  - Lattice triples are exact with E at m = −0.81q.

**(iii) Finite mass.**
- "ε small against (m/x_θ)²" hides a dependence on the amplitudes, contradicting the title:
  - x_θ = ln(|A|/θ − 1)/κ, and θ must be below min|A|, so x_θ grows like ln(range).
  - The inactive coefficients are Cε^{3/2}μ⁻³ (D‴.53a). These must stay below θ, which forces ε ≲ (min|A|·μ³)^{2/3}.
- The margins at entries are μ = 1 − |η| + O(q) only up to the screening above.

**(iv) Where the certificate fails.**
- The 3.3σ algebra is right, but only at a midpoint that is itself a starting point, with R₁ = R₂.
- Same-sign pairs fail well beyond 3.3σ (E1, 40 random offsets per distance): at 3.5σ the certificate gives m ≤ 0 on 40 of 40 and the path is 'other' on 16; at 4σ, m ≤ 0 on 15 and 'other' on 12.
- "Opposite signs … are recovered closer together" is contradicted. Exact-with-E counts:

| Separation | 1-D, equal amplitudes (of 40) | 2-D (of 16) |
|---|---|---|
| 1.5σ | 0 | – |
| 2σ | 14 | 2 |
| 2.5σ | 2 | 1 |
| 3σ | 25 | 11 |

  - The signed certificate says m > 0 on all 280 one-dimensional opposite pairs and is right on only 134.
  - On Part C the unsigned certificate predicts better than the signed one (96 vs 89 of 120).

**Numerics.**
- **Numbers.** Every number in the paragraph matches the md, except that "≤ 0.76q" rounds down a maximum of 0.7608q.
- **"135 of 135" is weak evidence.**
  - 131 of the 135 have m > 0 and are 'nearest', so the trivial rule "always nearest" scores 131/135.
  - The single-feature (Voronoi) margin, with no signs and no sums, scores 137/138 outside its own 0.8q band and 125/125 outside the q band.
  - Outside the q band there is no target with m ≤ 0, so necessity is untested there.
  - Signs and sums change the prediction on 125 targets; 121 of them lie inside the 0.8q band.
- **Where the 135 come from.** By minimal separation Δ = 0.25, 0.3, 0.35, 0.45, 0.6 they number 10, 14, 10, 35, 66. At Δ ≥ 0.25, 70 of 80 targets lie inside the band. The 80 targets at Δ ≥ 0.6 are one configuration with ±0.09 jitter.
- **The band constant.**
  - Chosen after the fact: 0.761 observed → 0.8 quoted (0.461 → 0.5 in 2-D).
  - It depends on the normalisation: without the −2R in q it must be 1.77, and "135 of 135" becomes "265 of 282".
  - With K = 4 or 5 features, 0 of 40 targets lie outside the band.
- **Part B.**
  - The certificate's invariance is a tautology.
  - The realised range at r = 3 has median 34 and maximum 749; the six configurations outside the band have ranges 1.24–376, not "up to 10³".
  - The amplitude order is fixed for r > 0, so "order changes on 36" is mostly r = 0 against r > 0 (10 configurations change among r = 1–3).
  - The sentence "keep their outcome at every range" is hard-coded in `report()`. I fed it a Part B in which one of the six (|m|/q = 12) changes outcome, and it still printed that sentence.
- **"Every path exact" and the polish.**
  - The flow alone reproduces the outcome and order on 397 of 400 (A) and 120 of 120 (C).
  - The three that do not finish (in 45 s, or in about 10 min for target 114) rely on a polish jump. Example: target 114's third entrant has w ≈ 0 and crawls from 0.598 to 0.587 by t = 4·10⁵; the polish moves it 0.45 (3.8σ) in one solve.
  - The polish tolerances do not matter (loosened to 1e-8: 0 changes), but the amount of flow before the polish does. With `t_fast=1e3`, 11 of 400 flip from 'nearest' to 'other', four of them at m/q 1.5–3.0. For target 267 the polish sends the new entrant from −0.2 to the feature at 0.6 (6.7σ away) instead of its own feature 0.027 away.

**Code at large ranges and edge cases.**
- At the full range 10³ with every ordering, 41 of 89 (positive) and 15 of 48 (negative) out-of-band paths do not finish in 20 s.
- At 10⁶, `exit_tol = 1e-6` makes the smallest kernel exit and re-enter: three paths end on E's kernels but count as 'other', and one chatters until `max_events = 60`.
- The stationarity tolerances are absolute, so they mean nothing for features much smaller than 1.
- With f* = φ(·;0), a single feature exactly on a kernel, `sliding_iss` hangs in LSODA (no result in 40 s). Offsets of 0.001 or 0.5 take 0.1 s.

**Tests and mutations.** 8 of the 13 mutants survive the tests:

| Mutant | Tests | Effect on the regenerated report |
|---|---|---|
| Normalise by the max-correlation kernel | survives | none (equivalent on this data) |
| Accept exits | survives | none (no exits in the data) |
| RMS_TOL 1e-3 | survives | none |
| RMS_TOL 1e-14 | survives | Part B: 29 → 28 configurations keep their outcome |
| One-sided band | survives | still "131 of 131" |
| Band with q² | survives | claim fails: 285 of 329 |
| q without −2R | survives | claim fails: 265 of 282, constant 1.77 |
| q without /G(R) | survives | 138 of 138 |
| Drop the sign | 1 test fails | 1-D 314/400, 132 of 135 outside the band; 2-D improves to 96/120 |
| Normalise by the entrant's self-correlation | 1 test fails | constant 5.2 |
| Normalise by the entrant's total correlation (amplitude-dependent) | 3 tests fail | predicts better: 382/400 and 108/120 |
| Include E in the max | 2 tests fail | 187/400 |
| Any exact fit counts as 'nearest' | 1 test fails | — |

Tests 1–2 test the definition against itself. Nothing tests `band()`, the exit logic, q for K ≥ 2, the 2-D job, or the md numbers.

**Novelty.** The literature points below are from memory unless marked web-checked.
- An amplitude-free certificate built from positions and signs is classical:
  - Fuchs' condition and the irrepresentable condition (Zhao–Yu; Wainwright).
  - Duval–Peyré's η_V.
  - Inverse scale space / linearized Bregman sign consistency, Osher–Ruan–Xiong–Yao–Yin 2016 (exists, web-checked: https://arxiv.org/abs/1406.7728).
  - The homotopy/LARS k-step property: any k-sparse solution with any amplitudes and signs in k steps when k ≤ (1 + 1/μ)/2 (Donoho–Tsaig 2008; restated in https://arxiv.org/pdf/1104.3792).
  - Sliding Frank–Wolfe terminates in finitely many steps under non-degeneracy (Thm 3 of https://arxiv.org/pdf/1811.06416).
- D55's η is a Gram-free approximation of such certificates restricted to the starting set.
- What is new (the sliding version on the starting set) is exactly the part left at an unspecified O(q).

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ap-1 | major | The band "\|m\| ≤ 0.8q" (and even ≤ q) is false on both sides. Sufficiency fails at m/q = 1.057 (opposite-sign pair on lattice points, ratio ≥ 2; polish off agrees) and at 0.80–0.96 (equal-amplitude opposite pairs at 2.5σ). Necessity fails at −1.62 (same-sign lattice pair, ratio ≥ 50), −0.90 (reordered Part A target) and −0.81 (lattice triples). | D55 Status and numerics: drop 0.8; report the largest \|m\|/q observed (≥ 1.62) as a lower bound; add these designed families to `exp_transit_certificate.py`. |
| D3ap-2 | major | "The amplitudes' range does not enter" fails outside the band. The same certificate (m = 1.057q, and −1.62q) gives 'nearest' at one amplitude ratio and 'other' at another. The final duals' deviation from η grows from ≤ 0.95q (comparable amplitudes) to 1.73q (ratio ≥ 50). | Title, (ii), (iii): say the band's constant depends on the order of entries, i.e. on the amplitudes. |
| D3ap-3 | major | (ii) omits screening: an earlier-fitted block absorbs part of a later feature, so the duals next to it fall short by O(G(Δ − h − R)/G(R)), which exceeds q when R ≪ h. T (radius 2R) degenerates at R = 0. q₀ and C are never given. | Proof (ii)/(iii): add the projection term, define q with the competitors' distance (Δ − h − R), give T a radius of at least cσ, state the constants or call it a conjecture. |
| D3ap-4 | major | (iv) "opposite signs cancel in η and are recovered closer together" is contradicted. Exact with E on 0, 14, 2 of 40 at 1.5, 2, 2.5σ (1-D) and 2, 1 of 16 at 2, 2.5σ (2-D), while the certificate says m > 0 on all. The signed certificate is right on 134 of 280; the unsigned one beats it in 2-D (96 vs 89/120). | D55 (iv): delete the claim. |
| D3ap-5 | major | The margin compares \|η(c_j)\| with 1 although η(c_{e_k}) = 1 ± O(q), so a "margin" of 0.22 hides a 0.01 near-tie (η = 0.7904 at E, 0.7784 at the intruder). The relative margin min_k s_kη(e_k) − max_{j∉E}\|η(c_j)\| predicts 344/400 (vs 331), 102/120 (vs 89), and 276/280 opposite pairs (vs 134). | Use the relative margin (or the projected certificate) in `certificate` and in D55's definition. |
| D3ap-6 | major | The out-of-band evidence is trivial. 131 of the 135 are m > 0 'nearest' ("always nearest" scores 131/135). The Voronoi margin without signs or sums scores 137/138 and 125/125. Necessity is untested outside q. 66 of the 135 are a single configuration. | `report()` and D55 numerics: add baselines; design targets with negative m outside the band; deduplicate the Δ ≥ 0.6 set. |
| D3ap-7 | major | The band constant is fitted after the fact (0.761 → 0.8; 0.461 → 0.5) and depends on the normalisation (without −2R: 1.77, 265 of 282). With K = 4 or 5, 0 of 40 targets lie outside the band, so the statement is vacuous there. | Fix C beforehand; report the needed C per normalisation and per K. |
| D3ap-8 | major | "Every path exact" depends on the polish. With `t_fast=1e3`, 11 of 400 outcomes flip (four at m/q 1.5–3.0) because `polish_block` jumps an entrant still at its seed (6.7σ in one case). Three committed paths depend on a 3.8σ polish jump that the flow alone has not completed (one not after about 10 min). | `rich_transit.sliding_iss`: polish only after the entrant has grown, or bound the polish step; cross-check against `polish=False` in the experiment. |
| D3ap-9 | major | Large ranges and lattice points are not computable with this code. At range 10³ with every ordering, 41 of 89 and 15 of 48 paths time out. At 10⁶, `exit_tol=1e-6` causes exit/re-entry and chattering. The `gtol` and tolerance checks are absolute. f* = φ(·;0) hangs LSODA. | Use relative tolerances, an exit threshold relative to the block's own amplitude, and a stiff integrator with the Jacobian (as `rich_transit_nd` does). |
| D3ap-10 | major | The finite-mass claim depends on the amplitudes: x_θ ∝ ln(range), the inactive coefficients force ε ≲ (min\|A\|μ³)^{2/3}, and the entries' margins are not 1 − \|η\| + O(q) (screening). | D55 (iii), last sentence: state the dependence on min\|A\| and max\|A\|; use the path's running margins. |
| D3ap-11 | minor | (i) is false on samples near the edge (c* = 0.955–0.98: m = −0.18 to −0.24, the path takes the kernel at 0.9). "ẇ has A's sign at w = 0" is wrong. | D55 (i): restrict to the continuum, or take e_k as the kernel of largest correlation; fix the wording. |
| D3ap-12 | minor | The 3.3σ threshold is a midpoint-only statement: same-sign pairs fail at 3.5σ (m ≤ 0 on 40/40, 'other' on 16/40) and at 4σ. "\|η\| ≥ 1 at a runner-up" cannot happen for one feature in the continuum. | D55 (iv). |
| D3ap-13 | minor | Part B is weak: realised ranges have median 34 and maximum 749; the six outside the band have ranges 1.24–376; the amplitude order is fixed for r > 0; the invariance of the certificate is a tautology. | `part_B`: set ranges explicitly and permute the orders. |
| D3ap-14 | minor | `report()` hard-codes "keep their outcome at every range" (demonstrated with a modified Part B). | `exp_transit_certificate.report`: compute the sentence. |
| D3ap-15 | minor | Tests are weak: 8 of 13 mutants survive, including ones that break the band claim (q without −2R, band with q²). | Add tests for `band()`, `outcome()` with synthetic paths (exits, tolerance), q for K ≥ 2, the lattice-pair counterexample, and md regeneration. |
| D3ap-16 | minor | `outcome()`: absolute `RMS_TOL`; 'other' also counts paths that end on E after exits; no check that e_k ends on feature k; `unique` checks distinctness, not ties in the nearest point. | `exp_transit_certificate.outcome` and `certificate`. |
| D3ap-17 | minor | The target generator avoids edges ([−0.7, 0.7]), uses a narrow amplitude range (factor 2.5), and puts most targets at q ≈ 0.15–1, beyond any plausible q₀. | `targets()`: add edge, ratio and lattice-point designs. |
| D3ap-18 | minor | Novelty is overstated: irrepresentable/Fuchs conditions, inverse scale space sign consistency, the homotopy k-step property, sliding Frank–Wolfe finite termination and η_V are not addressed. | D55 "Classical"/"New". |
| D3ap-19 | minor | Status "proved up to constants" and the title's "iff" overstate a sketch that excludes a band of unknown width. | Title and Status. |
| D3ap-20 | cosmetic | "≤ 0.76q" (actual maximum 0.7608q); η(e_k) should read sign(A_k)(1 + O(q)); `q` is reused as the configuration index in `part_B`. | D55 and code. |

## Lesson
A certificate claimed only up to an unspecified O(q) band, then "confirmed" with a constant fitted to the same random targets, cannot be falsified. Those random targets also put almost all the out-of-band mass on configurations where "nearest kernels always win" is already right. Designed families break both the band's constant and the amplitude-free claim within minutes: features exactly on starting points, opposite-sign neighbours, and amplitude ratios applied in every order. The polish shortcut also needs auditing, because it can decide the outcome it is supposed to compute. New certificates should be tested against trivial baselines and adversarial families, every ordering of the amplitudes should be varied, and the reported constants should be fixed before looking at the outcomes.
