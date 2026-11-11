# Referee report D3bl, round 65: Proposition D‴.58 as revised after round 60

## Scope

**Commit and files.** I reviewed commit 3eb9326 on branch claude/festive-bardeen-eopej2. During the review, other sessions moved HEAD to 6d48964. `git diff --stat 3eb9326 6d48964` is empty for every reviewed file, and the working tree is clean for them, so the review applies unchanged at HEAD. The files reviewed were:
- `review/drafts/D58.md`
- `experiments/exp_width_mismatch.py`
- `src/core/rich_transit.py` and `src/core/wide_feature.py`
- `results/chaos/width_mismatch.md` and `results/chaos/width_mismatch_raw.npz`
- `tests/test_width_mismatch.py`, `tests/test_rich_transit.py` and `tests/test_wide_feature.py`
- `review/REVIEW_round60_D3bg.md` and the "Round 60" entry of `review/REVIEW_LOG.md`

**How I worked.** Everything ran in a `git archive` export under `scratchpad/ref65/D3bl_work`; mutants ran in slim copies under `D3bl_work/mut/<id>`. Nothing was written to the repository.
- All runs used the venv Python with OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1 and `nice -n 10`, with at most two heavy processes of mine at a time.
- No process of mine is left running.
- Each number below is marked as coming from my run, the md or the raw.

## Summary

**What checks out**

- **Reproduction.**
  - `--from-raw` with WM_OUT pointing at a copy of the raw regenerates the md byte for byte (md5 3000ba35…). The raw is unchanged (79346a90…).
  - All 56 tests pass (52.5 s).
  - Four full B paths rerun bit-identically (all keys but `secs`): (2, sixth), (1, sixth), (¼, base) and (1, polish 2σ).
  - The sixth variant's three paths at ¼, ½ and ¾ also rerun bit-identically up to their first non-optimal block (6, 8 and 11 blocks). This covers events, `entry_info` and polish logs.
  - 14 of 14 B2 jobs rerun bit-identically.
- **Numbers.** Every number in (iii), (iv), the status and numerics A–D matches the md or raw, with the exceptions listed below.
  - From the raw: δ_t 0.0776→0.0397%, 0.728→0.395%, 3.477→2.123%; node deviations 0.0077% and 5.147%.
  - C(c) ≥ 0 after all 316 optimal blocks, to rounding (min C ≥ −9.5·10⁻¹⁵ absolute).
- **Counts.**
  - `cs_vacuous` is exactly the 53 non-optimal blocks with e ≤ 3.98·10⁻¹⁰; the smallest non-vacuous error is 4.32·10⁻¹⁰.
  - B2 above the floor: 33 of 38 reached, 1 near miss, 4 far misses, and 20 of 20 first blocks reached.
- **Sixth variant.** Its ranges 5, 7, 10, 12 are right. At ½, K = 8 the budget does bind: from the same state, 5·10⁵ evaluations reach 1.00027× (final state) and 1.00016× (first start) of the optimum.
- **Maths.** (i), (ii) and (v) are unchanged since round 60 (word-diff shows reflow only).
  - The Mercer ratio ε²/(α²+δ²+ε²) = (β²−1)/(β²+1) is correct. By quadrature the eigenvalue ratios are 0.029437, 0.101021, 0.267949, 0.500000, equal to q.
  - The Derevianko–Plonka (DP) identity f*−f_K = e^{−x²/(2σ²)}(g − Σb_k cos(η_kx/σ²)) is correct.

**What is wrong (majors first)**

1. **D3bl-1 (major).** "Where [a range] ends, the polish's budget, its acceptance radius or its gradient test binds" is false at the two range ends the sixth variant created.
   - Both lie below E_FLOOR: optima 2.55·10⁻¹¹ (¼, K = 6) and 7.75·10⁻¹¹ (¾, K = 11).
   - From the same states, 10–40× the budget leaves the polish 1.47% and 0.37% above the optimum. It stops on its step test (xtol); with xtol off it stops on the decrease test (ftol) at the same point. That is the double-precision floor, not any setting in the list.
2. **D3bl-2 (major).** 19 of 51 non-control single-point mutants survive both the tests and `--from-raw`.
   - The most serious is X01: `b_job` silently dropping the sixth variant's gradient setting. Its rerun at τ/σ = ¼ is bit-identical to the fifth variant's path (optimal K ≤ 4 instead of 5).
   - Also surviving: B2's condition-number and node-distance code, and the `--variants` CLI.
3. **Minor issues.**
   - B2's far misses are VP stalls at coalescing node pairs, not evidence about basins (D3bl-3).
   - The condition numbers depend on the length unit, the VP one is underestimated (3.8·10⁵, not 3.6·10⁵), and they do not explain VP's success: the smallest singular values are within 1.0–5.4× (D3bl-4).
   - The Mercer-frame sentence misidentifies the rule (D3bl-5).
   - The entry statistics count one event several times across variants: 316 entries are 107 distinct events, and "5" largest-rate entrants are 1 (D3bl-6).
   - The stationarity test is weak beyond the 53 vacuous blocks (D3bl-7).
   - The status-1 wording is incomplete. At the "one exception" block, all 24 polishes reach the optimum and the 1σ radius discards them (D3bl-8).
   - The status needs a floor qualifier (D3bl-9), plus cosmetics (D3bl-10).

## What I checked and how

- **Reproduction.** As listed in the Summary.
- **Range ends: same solver, same state, one setting changed at a time.**
  - I wrapped `polish_block` to record every call while rerunning traced paths (each prefix bit-identical to the raw).
  - From the recorded first-polish starts and from the blocks' final states, I reran TRF changing only one of: `max_nfev` (×10, ×40), `gtol`, `xtol` (None), `ftol` (None), `x_scale='jac'`.
  - Ends checked:
    - base ¼, K = 5
    - polish 2σ, τ/σ = 1, K = 10
    - sixth variant: ¼, K = 6; ½, K = 8; ¾, K = 11
    - the grid exception (τ/σ = 1, K = 9): all 24 of its polishes, with their end errors, gradients and the seeds' moves against the radius
  - I also reran the VP references with 10× the budget.
  - The other ends are bracketed by the published variants, which differ in one setting at a time.
- **B2.**
  - Recomputed both condition numbers.
  - Built an analytic Golub–Pereyra VP Jacobian (QR) to compare with the finite differences, plus σ_min, units and column equilibration.
  - Far misses: VP-LM rerun with status, residual–Jacobian cosines, node gaps and 10× the budget.
- **Entries.** Recomputed the regime split, then deduplicated by (width, grid, iteration, kernel).
- **Stationarity.** Per-block Cauchy–Schwarz thresholds, effective cosine thresholds, and margins at the first non-optimal blocks.
- **Literature algebra.**
  - Mercer eigenfunctions by quadrature.
  - Mercer-basis defects of three rules.
  - The family of Hermite-space isometries f ↦ e^{cx²}f(sx), derived by hand.
  - The first tail ratio for a general variance, checked numerically.
  - The √2 bound and the DP identity.
- **Mutation testing.** 52 single-point mutants (R01–R20 in `rich_transit.py`, X01–X32 in the experiment, X31 a no-op control).
  - Each was run against `test_rich_transit.py` + `test_width_mismatch.py` (R-mutants) or `test_width_mismatch.py` (X-mutants), plus a `--from-raw` byte diff.
  - Survivors were followed up with path reruns (X01, R12) and recomputed B2 output (X02, X04, X05, X07).
- **Literature search.** Four web searches, after which the shared search budget was exhausted. WebFetch to link.springer.com and ar5iv.labs.arxiv.org failed (DNS, outside this session's network). Everything literature-related below rests on search snippets or memory.

## Verdicts per part

| claim (D58.md line) | verdict | evidence |
|---|---|---|
| (i), proof (:8–14, :72–102) | Unchanged; correct as verified in round 60 | word-diff 3eb9326 vs fef9c70 |
| (ii) Mercer variance τ²/β², q, q > t (:15–23) | Correct | quadrature eigen-ratios = q to 6 digits (my run) |
| (iii) δ_t, 93.0–100%, nodes 0.008–5.2% (:25–28) | Correct | raw |
| (iii) C ≥ 0; flank peaks 2.0–2.6σ at K = 1 (:28–30) | Correct | 316 optimal blocks, min C ≥ −9.5·10⁻¹⁵ (my run) |
| (iii) entries: 5/316, 0.21, 0.959, 314/316 at 1.22σ, 2 at 0.989–0.992 (:31–34) | Numbers match; but they come from only 107 distinct events | D3bl-6 |
| (iii) ranges 4,5,7,9,12,12 / 4,7,9,11 / 4,7,9,12 / 5,7,10,12; optima to 6.8·10⁻¹⁰ (:36–40) | Correct | md; reruns bit-identical |
| (iii) "where one ends, budget, radius or gradient test binds" (:40–41) | Holds above E_FLOOR; **fails** at the sixth variant's two ends below it | D3bl-1 |
| (iii) gradient-test stops 1.012, 2.677, 1.121 (:41–43) | Correct; "below any fixed absolute tolerance" is overstated | md; D3bl-8 |
| (iii) sixth variant: one kernel further at ¼, ¾; ½ K = 8 at 1.233 (:43–44) | Correct; budget binds at ½ | +5·10⁵ → 1.00027 (my run) |
| (iii) VP 20 of 20 at the first non-optimal blocks above the floor (:44–46) | Correct | raw; B2 reruns |
| (iii) condition numbers 5.4·10³–3.6·10⁵ / 1.3·10⁶–5·10⁸; "far better conditioned" (:46–47) | As coded; analytic VP maximum is 3.8·10⁵; unit-dependent; not the mechanism | D3bl-4 |
| (iii) 33/38; 1.001; 4 at 16.36–24.74×, nodes 2.14–3.99σ (:47–50) | Numbers match; the four are VP stalls | D3bl-3 |
| (iii) 116/116 pass; 53 vacuous below 4.0·10⁻¹⁰; errors 1.0·10⁻¹²…1.3·10⁻⁶ (:50–52) | Correct | the other 63 are weakly tested (D3bl-7) |
| (iv) (:53–65) | Unchanged but for K*(s); numbers match | md |
| (v) (:66–70) | Unchanged; DP form verified | algebra |
| Status (:104–118) | Fair; needs the floor qualifiers | D3bl-1, -9 |
| Literature (:118–148) | Mercer ratio correct; Mercer-frame link wrong; √2 normalisation; Kuo–Woźniakowski (KW) unresolved | D3bl-5 |
| Numerics A, C, D (:150–162, :176–183) | All match | md |
| Numerics B (:163–176) | Match, except the status-1 enumeration | D3bl-8 |

## Status of the round-60 issues

| D3bg | status | note |
|---|---|---|
| 1 | partly | The mechanism is now the absolute gtol, and the sixth variant demonstrates it. But the closed list fails at the sixth variant's two below-floor ends (D3bl-1), and the replacement "conditioning" explanation does not hold up (D3bl-4). |
| 2 | mostly | The floor is applied: 33/38 recomputed. But the 4 far misses are VP stalls (D3bl-3). |
| 3 | partly | Round 60's surviving `sliding_iss` mutants are now killed: non-improving repolish accepted, extra repolish, `polish_nfev` ignored, no chunk reset, `entry_info`'s u0/rate_max, the accepted flags. But 19 new survivors remain, including X01 (D3bl-2). |
| 4 | fixed | "93.0–100%". |
| 5 | mostly | Split by regime; the counts still repeat the same events across variants (D3bl-6). |
| 6 | fixed | The 53 vacuous blocks are exact. The rest are covered in D3bl-7. |
| 7 | mostly | (b)–(e) fixed; (a) is still incomplete (D3bl-8). |
| 8 | partly | Gnewuch–Hefter–Hinrichs–Ritter (GHHR), Irrgeher–Kritzer–Leobacher–Pillichshammer (IKLP), Gnewuch–Ritter–Rüßmann (GRR), the √2 bound and DP's weighted problem were added. The Mercer sentence is wrong (D3bl-5); Karvonen–Oates–Girolami (KOG) is still unchecked. |
| 9 | fixed | K*(s) vs ⌈K*(s)⌉; "times 1 − δ_t(K)"; `.4g` formats. B2 still recomputes B's two VP candidates, which is harmless. |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bl-1 | major | **The sixth variant's two new range ends are the arithmetic floor's, not the budget's** (:40–41; status :114–115).<br>**The blocks:** optima 2.551·10⁻¹¹ (¼, K = 6) and 7.755·10⁻¹¹ (¾, K = 11), both below E_FLOOR. (iii) never says so; only the md's B2 line notes two such paths. The path's four 5·10⁴ polishes give 2420→13.5→8.70→6.61→4.951× and 1351→5.12→3.93→3.20→2.685×.<br>**One-setting reruns from the final state (gtol off), my runs:**<br>• +5·10⁴: 3.77× and 2.29× (budget).<br>• +5·10⁵ and +2·10⁶: identically **1.01474× (status 3, 256857 evaluations)** and **1.00365× (status 3, 384829)**.<br>• xtol off: ftol stops at the same points (+246 and +249 evaluations).<br>• `x_scale='jac'`: 1.0217× and 1.0144×.<br>• From the first polish's start with 5·10⁵: 1.0235× (xtol) and 1.0467× (budget; the newest kernel moved 2.08σ, beyond the 2σ radius).<br>• VP from those blocks also stops short (1.0069×, 1.1137×); VP from the half-squeeze nodes reaches the reference in 21 and 26 evaluations.<br>**Where the list does hold (my runs):**<br>• ½, K = 8: +5·10⁵ → 1.00027 (xtol, 74011).<br>• base ¼, K = 5: 5·10⁴ → 1.0017 (gtol, 31071).<br>• τ/σ = 1, polish 2σ, K = 10: 5·10⁴ → 1.000000.<br>• grid τ/σ = 1, K = 9: radius (D3bl-8). | (iii) and the status: "above E_FLOOR the budget, radius or gradient test binds; the sixth variant's ends at τ/σ = ¼, ¾ lie below it, where the joint solve stops on its step/decrease tests 0.4–1.5% above the optimum whatever the budget — those two ranges are censored by the floor". Mark below-floor first non-optimal blocks in B's table. |
| D3bl-2 | major | **Mutation testing: 52 mutants; the tests kill 31 and `--from-raw` one more (X25); 19 non-control mutants survive both.**<br>**Nine change behaviour:**<br>• **X01** (:153, `polish_gtol=gt`→1e-15): rerun (¼, sixth) is bit-identical to the committed fifth-variant path, optimal K ≤ 4 not 5. `test_b_paths_reproduce_the_raw` runs only base and the fifth variant at τ/σ = 2, where gtol is irrelevant.<br>• **X02** (:214, dc taken from the block's own VP end): far misses become 1.62–3.87σ instead of 2.14–3.99σ.<br>• **X04, X07** (:215, joint J without amplitudes / with unit amplitudes): 5.5·10⁵–1.6·10⁹ instead of 1.3·10⁶–5·10⁸.<br>• **X05** (:212, h = 10⁻⁴σ): VP 3.8·10⁵.<br>• **R01** (rich_transit.py:202, default gtol None): no effect on D‴.58, but `exp_transit_certificate.py:127` and `exp_noisy_certificate.py:358` call `polish_block` without gtol.<br>• **X28–X30** (:665, :668, :178, the `--variants` CLI): off by one, no merge, inverted filter — the path that built the sixth variant's raw.<br>**Ten are equivalent on the committed data, but their boundaries are untested:**<br>• R12: (1, sixth) and (½, fifth) rerun bit-identical.<br>• R17: `t_min` is unused.<br>• R18: no ties in 36 paths.<br>• X03: centres are already stored sorted.<br>• X09, X14, X18, X23, X24: boundary changes (the 1.01/1.02 miss limits, the CS equality, d = 0, the 10⁻¹² top tolerance) that never occur in the data.<br>• X20: all 116 blocks relaxed.<br>**Killed by tests:** R02–R11, R13–R16, R19, R20, X06, X08, X10–X13, X15–X17, X19, X21, X22, X26, X27, X32. | Tests:<br>• `b_job` passes `VARIANTS_B[i][5]` (monkeypatched `sliding_iss`), or the sixth variant at ¼ to 6 entries.<br>• `b2_job`'s condition numbers against an analytic VP Jacobian on a 2–3-kernel case.<br>• dc where the two starts end apart.<br>• `main --variants` with a stub `b_job`.<br>• Pin `polish_block`'s default.<br>• The listed boundaries on synthetic data. |
| D3bl-3 | minor | **B2's far misses are VP stalls** (:50; md :151). From all four blocks, VP-LM ends with two nodes coalescing (gaps 0.000, 0.010, 0.021, 0.012σ). The endpoints are not stationary: the largest residual–Jacobian cosine is 2.6·10⁻⁴–4.7·10⁻³ against MINPACK's 10⁻¹⁵ cosine test. Three stop on xtol; ¾ base K = 10 stops on its 2·10⁴ budget (2·10⁵ evaluations: 24.74 → 20.14×, then xtol). By the draft's own rule, a miss says nothing of a basin. | Say "VP stalls at a coalescing pair", or rerun VP with a separation guard or perturbed starts before quoting 16–25× and the distances (:50, md). |
| D3bl-4 | minor | **The condition numbers** (:46–47).<br>**(a) FD accuracy.** The analytic Golub–Pereyra VP Jacobian agrees with the finite differences to 8·10⁻⁸–6·10⁻⁶ (Frobenius) but gives 3.83·10⁵, not 3.48·10⁵, at (½, 5·10⁴ once, K = 8); below the floor 1.71·10⁶ vs 7.40·10⁵.<br>**(b) Units.** With c measured in σ the joint range is 4.5·10⁶–2.3·10⁹. With equilibrated columns: joint 1.7·10⁵–4.2·10⁸, VP 1.6·10³–2.5·10⁵.<br>**(c) Not the mechanism.** σ_min(VP)/σ_min(joint) is 1.00–5.35 at 19 of 20 blocks (median 1.54; 31.8 at grid ¾, K = 8): the slow modes are equally slow, and the gap lies at the stiff end. A Jacobian-scaled TRF from the same starts is no faster: 11.45 vs 10.85 (¼ base, 5000); 1.0017 vs 1.0017 (5·10⁴); 1.92 vs 1.77 (τ/σ = 1, K = 10); 17.3 vs 4.58 (½ sixth, 5·10⁴). | Quote 3.8·10⁵ (analytic) with the units; drop "far better conditioned" as the explanation, or support it with σ_min and a scale-free comparison (:46–47, md). |
| D3bl-5 | minor | **The Mercer sentence** (:128–133) [snippets for the papers; algebra mine].<br>**(a) Wrong rule.** In the Mercer basis the Mercer-zero rule's defects do not vanish below 2K: at τ/σ = 1, K = 3 they are 1.5·10⁻³, 3.3·10⁻², 1.8·10⁻¹ (relative, n = 0, 2, 4). The Gauss–Hermite (GH) rule whose defects do vanish there (≤ 5·10⁻¹⁶) is that of N(0, 2τ²/(β²+1)) = 1.0004–1.118·v_hs. It falls at 63.1, 17.15, 5.32, 2.33 per kernel (τ/σ = ¼…2), near 1/t, not 1/q = 33.97, 9.90, 3.73, 2.00.<br>**(b) Wrong measure.** The Mercer basis's Hermite polynomials are orthonormal under N(0, v_M), not under the transformed integration measure. The frame where the two coincide has Mehler ratio λ, not q.<br>**(c) What actually singles out the half-squeeze.** In every frame of the family the tail ratio D_{2K+2}/D_{2K} = −v′K(K+2) + (v′−1)(K+1)(2K+1) (checked numerically) grows like (v′−2)K² unless v′ = 2. Only the half-squeeze frame has degree-j ratios.<br>**(d) Normalisation.** e < √2 holds for the normalised RKHS error; in relative L² it is (2σ*/σ)^{1/2}.<br>**(e) KW unresolved.** "2αγ² > 1" appears in none of the snippets I found. KOG's restatement is "converges … if ℓ > α"; Karvonen–Särkkä's is "converges iff ε² < 1/2" — convergence thresholds that the √2 bound does not address. | Rewrite :128–133 per (a)–(c); mark KW as unresolved. |
| D3bl-6 | minor | **Entries counted across variants** (:31–34; md :158–161). The 316 entries after optimal blocks are 107 distinct events: 49 occur in five variants, with duals equal to 6.4·10⁻⁶. The "5 of 316" largest-rate entrants are one event (τ/σ = 2, iteration 12, kernel 0), i.e. 1 of 107. On distinct events the medians are 0.208, 0.963 and 1.17σ (105 beyond, 2 between). Likewise first entries 36 = 12 distinct and non-optimal 116 = 87. | Count distinct events or report per variant (`entry_lines`; :32). |
| D3bl-7 | minor | **The stationarity test is weak beyond the 53 vacuous blocks** (:51). For the other 63 its effective cosine threshold is 3.5·10⁻⁴–0.97 (median 7.7·10⁻³; 24 above 0.1). Where every polish was rejected (base ½, ¾; grid ½, ¾, 1), the blocks passed at 0.39–0.77 of the threshold while 4–14× above the optimum. | Quote the relative margins, or make the test relative (:51; `rich_transit.py:431`). |
| D3bl-8 | minor | **Status-1 wording** (:172–173; (iii) :43).<br>• Status 1 also ends the once-variant's polish at ¼ (1/30240/a, 1.012) and the third repeated polish at ½ and ¾ (1/1/r).<br>• At τ/σ = 1, grid +0.05, K = 9 (the "exception") **all 24 polishes reach the optimum**: e₁/e_opt = 1.00000000, max\|g\| ≤ 7·10⁻¹⁶. Each is rejected because the newest kernel must move 1.88→1.32σ, beyond the 1σ radius. The block then ends 4.24× above the optimum.<br>• "Below any fixed absolute tolerance" should read "below 10⁻¹⁵". | Fix the wording. In B's table, distinguish status 1 at the optimum from status 1 above it. |
| D3bl-9 | minor | **Status floor qualifier** (:115–116). "The first non-optimal one in VP's basin" holds on 20 of 22 paths. At the two below the floor, VP from the block stops at 1.0069× and 1.1137×. | Add "whose optimum is ≥ 10⁻¹⁰". |
| D3bl-10 | cosmetic | "The block's nodes up to 2.14–3.99σ" reads as a range of maxima. `b2_job` recomputes B's own two VP candidates (identical to B's e_opt). | Wording. |

## Sources

All literature points rest on search snippets or memory; I read no full text.
- KOG, *Integration in RKHS of Gaussian kernels*, [arXiv:2004.12654](https://arxiv.org/pdf/2004.12654) and [html](https://arxiv.org/html/2004.12654). Snippets: "KW 2012 and KSW 2017 analyse integration … under the constraint ℓ > α"; "the Gauss–Hermite rule converges with an exponential rate if ℓ > α"; "the matching lower bound claimed in KW 2012 is erroneous".
- Karvonen–Särkkä, [arXiv:1803.09532](https://arxiv.org/pdf/1803.09532). Fragmentary snippet: "converges iff ε² < 1/2".
- Derevianko–Plonka, *Differential approximation of the Gaussian by short cosine sums with exponential error decay*, [STSIP 24, art. 6](https://link.springer.com/article/10.1007/s43670-026-00124-9) and [arXiv:2307.13587](https://arxiv.org/abs/2307.13587). Snippet: weighted L² with weight e^{−t²/(2ρ)}; frequencies at the zeros of a scaled Hermite polynomial.
- GHHR, [arXiv:2110.05778](https://ar5iv.labs.arxiv.org/html/2110.05778) (title only).
- Kuo–Woźniakowski, BIT 52(2):425–436 (2012): reference entry only.
- link.springer.com and ar5iv.labs.arxiv.org were unreachable (DNS); the shared web-search budget ran out after four searches.

## Lesson

A list of causes is itself a claim, and it has to be checked at every end it covers. The revision replaced one unchecked mechanism with a list of three, then added a variant whose two new ends sit below the experiment's own floor. There no budget helps, and the step and decrease tests stop the solve 0.4–1.5% short; checking that took two reruns from recorded states. A diagnostic borrowed from a referee needs the same scrutiny: the condition numbers differ only at the stiff end, while the slow modes are equally slow in both formulations. And code that produces the round's main new evidence needs a test that fails when that code is removed.
