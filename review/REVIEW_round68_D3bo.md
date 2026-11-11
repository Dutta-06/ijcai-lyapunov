# Referee report D3bo, round 68: Proposition D‴.59 as revised after round 61

## Scope

I reviewed commit 2391685 from a `git archive` export at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref68/D3bo_work/repo`. The files were:
- `review/drafts/D59.md`;
- `experiments/exp_moving_width.py`, `results/chaos/moving_width.md` and `moving_width_raw.npz`;
- `tests/test_moving_width.py` (17 tests);
- `src/core/moving_width.py`, `src/core/rich_transit_nd.py` and `src/core/rich_transit.py`;
- `review/REVIEW_round61_D3bh.md` and the Round 61 entry of `review/REVIEW_LOG.md`.

HEAD moved to 766c864 while I worked. `git diff 2391685 HEAD` is empty for every file above, so the findings also hold at HEAD.

I modified nothing in the repository. Everything I wrote is under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref68/D3bo_work/`: scripts in `scripts/`, logs in `logs/`, loop states in `out/`, and the `--from-raw` output in `mwout/`.

Every run used the venv python with `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and `nice -n 10`. I aimed to run at most two heavy processes at a time; three may have overlapped for a minute or two at one hand-off. I signalled no process. Nothing of mine is still running.

Every number below is tagged [md] (the committed md), [raw] (the committed raw) or [own] (my run).

## Summary

**What checks out**
- **Reproducibility.**
  - `--from-raw` regenerates the md byte for byte [own].
  - The 17 tests pass in 9.6 s [own].
  - These reruns are bit-identical to the raw in every field [own]:
    - H2 "σ* 0.06, equal" and H2 "σ* 0.06, c₂ +0.01" at tolerance 10⁻¹³;
    - H3 at ±0.1 and ±0.12;
    - F target 2.
  - I reran the G loop for all six targets at ε = 10⁻⁷, 10⁻⁸, 10⁻⁹ and 10⁻¹⁰ (24 cells × 56 grid times). Every grid state agrees with the raw to ≤ 4.3·10⁻¹⁴ [own].
- **Numbers.** Every number in the statement, Status and Numerics matches the md or the raw. I found no misquote.
- **(i) as corrected.**
  - E = J + ½‖φ‖²(w − w*)² holds.
  - E ≤ ½‖f*‖² − η gives (w − w*)² ≤ w*² − 2η/‖φ‖². So w keeps w*'s sign and stays in (0, 2|w*|).
  - The compactness and Łojasiewicz steps hold.
  - J, ∂_cJ, ∂_σJ and the single critical point re-derive.
- **(ii)(a).**
  - The k-th moment (w_iw_j/W^k)((−1)^k w_j^{k−1} + w_i^{k−1})Δc^k re-derives. For a cancelling pair it is ≈ (k−1)|W|(vΔc/(|W|σ))^k, so the stated expansion parameter is right.
  - The split (w_j − w_i)Δc²/(6Wσ) re-derives.
  - The 0.9828–1.0458 check matches the raw (72 rows, of which 62 are distinct states).
- **(ii)(c), the block's law.**
  - v³ = 3μ₀Ks, K = −⟨r_dip, ∇³φ:p³⟩/24, and the vanishing second-order term (from stationarity in θ̄ and p) all re-derive.
  - Λ(ε) is the right fast time: it is the loop's own clock, since `loop_nd` integrates t = Λτ.
  - Comparing slopes with K and μ₀ from the path's flow and τ_d from the path's entry is fair.
- **The dipole phase in the loop holds better than the draft shows.** On a 0.1-spaced τ grid, v at the true τ₃ is the following fraction of (3μ₀KΛ(τ₃ − τ_d))^{1/3} [own]:
  - 0.998–1.049 at ε = 10⁻¹⁰;
  - 0.9995–1.013 at ε = 10⁻¹².
- **The ε^{−1/6} scaling holds.** I extended the run to ε = 10⁻¹¹ and 10⁻¹². The ratio of v at τ = 400 over the last decade is 1.4677, 1.4469, 1.4675, 1.4689, 1.4683 and 1.4728, against 10^{1/6} = 1.4678 [own].
- **τ₃ can be predicted from the path alone.** Advance the path's duals at τ_d with the residual of the block flow at s = 10¹¹. This gives τ₃ = 29.16, 58.44, 34.23, 36.17, 47.19 and 29.63. The loop at 10⁻¹² gives 29.2, 58.4, 34.3, 36.2, 47.2 and 29.7 [own]. So "τ₃ − τ_d settles" can be proved rather than assumed.
- **(iii), E and F** are as stated.

**What is wrong (majors first)**
1. **"The cancellation is the path's" and "which the loop does not keep" are contradicted (D3bo-1).**
   - The real loop at ε = 10⁻⁸ forms and keeps the cancellation to τ = 10⁵ on all five asymmetric H targets: 11–12 negative kernels, max|w| 2.19–2.69, sub-sample widths [own].
   - HL tested only the two symmetric targets. On "σ* 0.06, equal" the path's 1 + 8 cluster comes from its least-squares polish, which hit its 5000-evaluation cap.
   - With the polish off, the same path stays symmetric with max|w| 0.929, as the loop does [own].
2. **The pair does not land at τ₃ (D3bo-2).**
   - On 5 of 6 targets, τ₃ opens an intermediate phase with three formed kernels in which the pair sits on no feature. The pair lands only after a fourth kernel forms.
   - v changes on the way by a factor κ = 0.75–1.12. This factor depends on the target, not on ε.
   - The md's "0.71–1.09" mixes κ with a τ₃ read from a grid that is up to 8 units late.
   - "Landed" (within 1%) is the threshold p = |w_a − w_b||Δc|/(|W|σ) ≲ 0.28 on (a)'s kinematics: |W − A| ≈ p²|A|/8 and offset ≈ p³σ/24. It is not a separate phase.
3. **H2's "keeps entering kernels at exploding times while the residual creeps down" is a set of numerical artefacts (D3bo-3).**
   - 4 of 6 paths stop because all 21 kernels are active.
   - The entry at τ = 2.34·10⁶⁷ is the re-entry of a spike of negative width whose largest φ on the samples is 2.7·10⁻⁵³.
   - The polishes stop at max_nfev = 5000 and are accepted anyway.
   - With a converged polish, the same 17 kernels of H's end reach rms 5.6·10⁻¹² instead of 4.3·10⁻¹⁰.

Minors:
- **The collapse (D3bo-4).** The asymptotic law is τ^{−1/6} for the distance and τ^{−2/3} for the residual, because E ∝ Δc⁸ with the centre of mass pinned. "Frozen" is the same law: target 9 starts collapsing at τ ≈ 10⁸.
- **Mutants (D3bo-5).** 23 of 73 survive; 13 of them are not equivalent. The H2/H3/HL job code is still never run by a test or by `--from-raw`.
- **(i)'s "from any start" (D3bo-6)** is overstated.
- **Literature (D3bo-7).** The collapse belongs to a known class: over-parameterized fits converge algebraically.
- **Leftovers (D3bo-8).** C's table still prints pre-block states.
- **Cosmetic (D3bo-9).**

## What I checked and how

1. **Reproducibility.**
   - `--from-raw` with `MW_OUT=mwout/`, compared with `cmp`.
   - `pytest -p no:cacheprovider tests/test_moving_width.py`.
   - `rerun_a.py` (h_job on H3 ±0.1 and ±0.12; f_job on target 2) and `h2_hist.py` (H2 "equal" and "c₂ +0.01"). Each was compared field by field with the raw, keeping each path's events and history.
   - `loop_fine.py`: `loop_nd` with TAUS_G plus a 0.1-spaced grid on [10, 80]. The end stays at 6400, so the solver takes the same steps as G. 36 cells: six targets × ε = 10⁻⁷…10⁻¹².
2. **Numbers.** Every number in D59.md:7–163 against the md and the raw (`raw_overview.py`, `raw_h.py`, `raw_c.py`, `raw_hl.py`).
3. **Theory, re-derived by hand.**
   - (i): E = J + ½‖φ‖²(w − w*)², ∂_σE at (c*, σ*).
   - (ii)(a): the k-th moment.
   - (ii)(c): the block law; the second-order term; |W − A| ≈ p²|A|/8 and offset ≈ p³σ/24 from (a)'s expansion.
   - The collapse: E ∝ Δc⁸ with mobility Gm·W²/(w_a² + w_b²) for the separation at a pinned centre of mass.
4. **(ii)(c) in the loop.**
   - `an_phases.py`, `an_mid.py`: the true τ₃, v(τ₃⁻), κ, τ₄ and the landing time.
   - `t3_pred.py`: τ₃ from the path's duals.
   - `land_param.py`: the landing flag against p on 1116 stored states.
   - `collapse_energy.py`: E(Δc) at frozen v, everything else optimised by least squares.
   - `loop_long.py`: target 20 to τ = 10⁷ and target 9 to τ = 10¹⁰ at ε = 10⁻¹⁰.
   - `frozen_tol.py`: target 9 with rtol 10⁻¹¹ and atol 10⁻¹⁴.
   - `frozen_rhs.py`: the right-hand side in 40-digit decimal arithmetic.
5. **(ii)(b).**
   - `an_h2*.py`: the H2 events, the polish logs, and the kernels behind the last entries.
   - `polish_end.py`: H's stored "equal" end polished at 5·10³, 5·10⁴ and up to 4·10⁵ evaluations.
   - `h_variant.py`: H's path with `polish=False` on "equal", "c₂ +0.01", "A₂ 0.9" and "σ* 0.04, equal".
   - `hl_long.py`: the loop to τ = 10⁵ on "equal" (10⁻⁶), "A₂ 0.9" (10⁻⁶, 10⁻⁸), "A₂ 0.95" (10⁻⁶, 10⁻⁸), "A₂ 0.99", "σ₂ 0.065" and "c₂ +0.01" (10⁻⁶, 10⁻⁸).
6. **Mutants.** `mutants.py` and `run_mutants.py`: 73 single-point mutants. They cover `phases`, `g_phase_summary`, the landed-pair split line, `pair_fits`, `outcome`, `pair_index`, `h_stats`, `clusters`, the H2/H3/HL report and job code, `oracle_error`, `landing`, `dipole_pair` and the L flag. Each mutant ran against the test file and against `--from-raw` with an md byte diff (`logs/mutants.log`).
7. **Literature.** No source could be fetched (see Sources). Everything is from memory.

## Verdicts per part

| claim (D59.md line) | verdict | evidence |
|---|---|---|
| (i) J, its gradient, one critical point, the limits (:7–12, :59–61) | holds | re-derived; B [md] |
| (i) the block converges from a seed start; the corrected sign bound (:62–67) | holds | re-derived |
| (i) "from any start whose coefficient has the feature's sign" (:13) | overstated | the proof needs E(start) < ½‖f*‖² (D3bo-6) |
| (ii) C: 71 exact / 3 pairs / 6 dipoles (:18–19) | holds | [md]; `outcome` recomputed on the stored ends |
| (ii)(a) expansion, split, validity, k-th moment (:24–26, :67–73) | holds | re-derived |
| (ii)(a) the landed pairs' split 0.9828–1.0458 (:139) | holds | [raw]; 62 distinct states among the 72 rows |
| (ii)(b) the path's ends: 11–14 negatives along the path, 6–7 at the end, \|w\| ≤ 3.24 (:27–30, :142–147) | holds as a description of the polished path | [md]; on "equal" the cluster is the polish's (D3bo-1) |
| (ii)(b) "does not end exact … exploding times … creeps down" (:29–31, :147–150; Status :90–91) | not supported | dictionary exhausted; an invisible negative-width entrant; capped polishes (D3bo-3) |
| (ii)(b) the separation window (:31–32, H3) | holds for the polished path | ±0.1 and ±0.12 bit-identical [own] |
| (ii)(b) "The cancellation is the path's" (:33–34); New "which the loop does not keep" (:110–111) | contradicted | the loop keeps it on all five asymmetric targets (D3bo-1) |
| (ii)(c) non-closedness; v³ = 3μ₀Ks (:35–41, :73–82) | holds | re-derived; G [md] |
| (ii)(c) the loop's phase v³ = 3μ₀KΛ(τ − τ_d) (:42–45, :82–84) | holds, better than shown | 0.998–1.049 (10⁻¹⁰), 0.9995–1.013 (10⁻¹²) at the true τ₃ [own] |
| (ii)(c) "until the next kernel forms, at τ₃; the pair then lands … v near (…)^{1/3}" (:45–48) | wrong on 5 of 6 targets | intermediate phase; landing after τ₄; κ = 0.75–1.12 (D3bo-2) |
| (ii)(c) v ∝ ε^{−1/6} "while τ₃ − τ_d settles" (:48, :84) | holds, and provable | τ₃ predicted from the path; ratios 1.4469–1.4728 at 10⁻¹¹→10⁻¹² [own] |
| (ii)(c) "≈ τ^{−0.2}, residual ≈ τ^{−0.7…−0.9}, or frozen"; "exact only as τ → ∞" (:48–50) | pre-asymptotic; "frozen" is the same law; "only as τ → ∞" is right | −1/6 and −2/3 (D3bo-4) |
| (ii)(c) at ε ≥ 10⁻⁶ no landing or no dipole (:50–51) | holds | [md] |
| (iii), E, F (:53–57, :154–163) | holds | F target 2 bit-identical [own] |
| Status and New (:87–111) | overclaim | D3bo-1, D3bo-2, D3bo-3, D3bo-7 |

## Status of the round-61 issues

| D3bh | round-61 issue | status | what remains |
|---|---|---|---|
| 1 | (ii)(b) "exact" was the tolerance; H's ends; widths; σ* 0.04's triple; pair rule | partly fixed | Fixed: wording, H columns, `outcome`'s null-cluster classes, H2 added. Remaining: H2 is dictionary exhaustion, an invisible entrant and capped polishes (D3bo-3); the "equal" cluster and σ* 0.04's null triple are the polish's (D3bo-1) |
| 2 | the loop sentence; the missing ε-law | partly fixed | Fixed: phases, the ε^{−1/6} law, centre of mass and W − A, pair column at report time. Remaining: the landing follows τ₄ with factor κ, "landed" is p ≲ 0.28, the grid τ₃ (D3bo-2); the collapse law and "frozen" (D3bo-4) |
| 3 | the cancellation's condition; (b) the path's | condition fixed; conclusion wrong | 1.33σ₀ / within 10% and the H3 scan are fine. HL ran only symmetric targets, and the loop keeps the cancellation on the asymmetric ones (D3bo-1) |
| 4 | (a)'s check and its constant | fixed | 10 repeated rows (D3bo-9) |
| 5 | the pairs' mechanism | fixed | — |
| 6 | F's drift; (iii)'s noise kernel | fixed | — |
| 7 | mutants; the job functions untested | partly fixed | 23 of 73 survive, 13 not equivalent; h_job, hl_job and h3_targets are never run (D3bo-5) |
| 8 | raw provenance; C's table | mostly fixed | A–D regenerated (the `stalled` key is present); `oc` is `outcome`; polished \|w\| marked. C's table still prints the pre-block state (D3bo-8) |
| 9 | (i)'s proof | fixed | "from any start" (D3bo-6) |
| 10 | σ unconstrained | stated and counted | its consequences now enter the claims: the 10⁶⁷ entrant and invisible "formed" spikes (D3bo-3) |
| 11 | literature and "New" | partly fixed | New's cancellation clause is contradicted; the collapse class is known; one citation unverified (D3bo-7) |
| 12 | misquotes | fixed | — |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bo-1 | major | **"The cancellation is the path's" (D59.md:33–34) and New's "which the loop does not keep" (:110–111) are contradicted.** <br>• **The loop on asymmetric targets.** At ε = 10⁻⁸, to τ = 10⁵ [own]: <br>  – c₂ +0.01: 13 formed / 11 negative; max\|w\| 2.02 at τ = 6400 → 2.19. <br>  – A₂ 0.99: 14/12; 2.38. <br>  – A₂ 0.95: 14/12; 2.37. <br>  – A₂ 0.9: 14/12; 2.21 at τ = 10⁴ → 2.37. <br>  – σ₂ 0.065: 13/11; 2.69. <br>  Every one has formed widths below the sample spacing, and rms 1.6–2.1·10⁻⁶. At 10⁻⁶, A₂ 0.9 and A₂ 0.95 keep 6 negatives and max\|w\| 1.55 → 1.82 (τ = 1600 → 10⁵), with three spikes of width 0.001. Only c₂ +0.01 at 10⁻⁶ is exact with two kernels (rms 4.6·10⁻⁹). <br>• **HL tested only symmetric targets.** σ* 0.04 is not a cancellation target: its path end is "exact, null cluster". On σ* 0.06, equal, the path's 1 + 8 cluster is the polish's [own]: <br>  – At τ = 6968 kernel 8 enters without its mirror 12, which enters at τ = 5.0·10¹¹. <br>  – That block's polish stops at max_nfev = 5000 (status 0), is accepted, and multiplies max\|w\| by 3.358. <br>  – With `polish=False`, 8 and 12 enter together and the path stays symmetric: ±0.080 at w 0.929, wide negatives (σ ≥ 0.05, two centred at ±1.694, outside the window), max\|w\| 0.929, rms 9.7·10⁻⁶ with all 21 kernels. <br>  – The loop does the same: at 10⁻⁶ it stays symmetric to τ = 10⁵ (±0.081 at w 0.572 → 0.705, negatives shrinking to −0.052). <br>• σ* 0.04's null triple (−0.0546/+0.1093/−0.0546 [raw]) is also the polish's: without it the path ends with two kernels of w 1.000 [own]. | (b) and New: say that the cancellation is the path's and the loop's on asymmetric pairs, and that on the exactly symmetric pair the polish breaks the symmetry. HL: run the asymmetric H targets at 10⁻⁶ and 10⁻⁸ to τ ≥ 10⁵, and report max\|w\| against τ — it keeps growing, so say whether this is a slow run-off. `sliding_iss_nd`: flag or reject polishes that stop at max_nfev or grow max\|w\| beyond a set factor. H: add a `polish=False` control. |
| D3bo-2 | major | **The landing (D59.md:45–48, :92–93, :134–138).** <br>• **An intermediate phase.** The five targets whose dipole block has two kernels pass through a phase with three formed kernels between τ₃ and τ₄. In it the pair sits on no feature. At 10⁻¹⁰ [own]: <br>  – 18: offset 2.84σ*, (W − A)/\|A\| = −0.457 over τ = 34.2–80; v 15.1 → 12.7. <br>  – 9: 9.24σ*, −0.70; v 18.3 → 20.0. <br>  – 12: 20σ*, +0.86. <br>  – 29: 7.5σ*, +1.33. <br>  – 27: 0.73σ*, +0.16. <br>  The pair lands only after a fourth kernel forms (τ₄ ≈ 47–132). Only target 20, whose block already has three kernels, lands right after τ₃. <br>• **v changes on the way.** κ = v_land/v(τ₃⁻) is 1.116, 0.915, 0.753, 1.006, 1.021 and 0.769 for targets 9, 12, 18, 20, 27 and 29 at 10⁻¹². It moves by ≤ 0.6% between 10⁻¹⁰ and 10⁻¹² [own]. <br>• **The grid τ₃ is late.** The md's grid τ₃ at 10⁻¹⁰ is 30, 65.98, 36, 38, 48 and 30 [md]; the true values are 29.1, 58.0, 34.2, 36.2, 47.1 and 29.6 [own]. This biases the law's v upward by up to 13%. The md's "0.71–1.09" is κ divided by this bias. <br>• **"Landed" is a threshold.** In (a)'s regime \|W − A\|/\|A\| ≈ p²/8 and the offset ≈ p³σ/24; the stored states give 0.123p² and 0.0415p³ [raw]. So "landed within 1%" means p ≲ 0.283. The landed rows have p ≤ 0.287, and the rule p ≤ 0.28 reproduces the flag on 928 of 982 rows with p < 1 [raw]. At 10⁻¹⁰ the first landing falls at τ = 47.7–400; at 10⁻⁸ it falls as late as τ = 3200 [own]. | (c) and the proof: describe four stages — the dipole law; the intermediate phase τ₃ → τ₄; (a)'s regime from τ₄; the collapse. Derive τ₃ from the path's duals (my prediction matches to ≤ 0.07). State v_land = κ(3μ₀KΛ(τ₃ − τ_d))^{1/3} with κ measured. Define "landed" by p, or drop it as a phase. G: locate formation times by events or on a fine grid, and add ε = 10⁻¹¹ and 10⁻¹². |
| D3bo-3 | major | **H2 and "does not end exact … keeps entering kernels at exploding times while the residual creeps down" (D59.md:29–31, :90–91, :147–150).** <br>• **The paths run out of kernels.** 4 of the 6 H2 paths stop because all 21 kernels are active, so no entry time is finite [raw: n_act 21]. "equal" makes 24 entries and 3 exits in 18 blocks; "c₂ +0.01" makes 32 and 11 in 31 blocks [own]. The md labels these ends "cancellation". <br>• **The τ = 2.3412·10⁶⁷ entry is meaningless.** It is kernel 18 re-entering at c = 0.7754 with σ = −2.94·10⁻⁴. Its largest φ on the samples is 2.7·10⁻⁵³ and its correlation 8.5·10⁻⁶⁸. It had left at τ = 3.16·10¹⁴ after four exits and re-entries at that one τ, its width going 0.086 → 0.0086 → −3.0·10⁻⁴ [own]. <br>• **"Formed" counts invisible spikes.** In "c₂ +0.01" the right feature's "2 kernels" are w = 1.000 plus a spike of −0.018 and width 1.75·10⁻⁴, whose largest φ on the samples is 6·10⁻⁶⁹ [own]. <br>• **The polishes are capped.** From the cancellation's block on, 8 of 9 ("equal") and 14 of 18 ("c₂") polishes stop at max_nfev = 5000 (status 0) and are accepted [own]. <br>• **A converged polish goes much lower with the same kernels.** Polishing H's stored "equal" end with the same sign-constrained least squares takes rms 4.26·10⁻¹⁰ → 2.42·10⁻¹⁰ (5000 evaluations) → 2.51·10⁻¹¹ (5·10⁴) → 5.63·10⁻¹² (converged after 155 778). The 17 kernels and the 1 + 8 clusters stay the same, and max\|w\| goes 3.249 → 3.300 [own]. So the 10⁻⁹ ends and the 10⁻¹³ creep measure the polish's cap. Whether the cluster reaches zero, or runs off with growing \|w\|, is open. | H2: report the stop reason. Drop "exploding times" and the 10⁶⁷ entry, or constrain σ > 0 (log σ). Exclude kernels with negligible φ on the samples from "formed". Raise or adapt max_nfev and flag status-0 polishes. Status: "at the tolerance the residual is the polish's". |
| D3bo-4 | minor | **The collapse: "≈ τ^{−0.2}, residual ≈ τ^{−0.7…−0.9} … or frozen" (D59.md:49, :140–141).** <br>• **E ∝ Δc⁸.** E(Δc) at frozen v, everything else optimised, has local slope d ln E/d ln Δc = 8.008, 8.001, 7.99 and 7.97 (targets 20, 29, 9, 18 at 10⁻¹⁰) [own]. <br>• **The law.** The centre mobility is Gm, and the centre of mass is pinned, so the separation's mobility is Gm·W²/(w_a² + w_b²). This gives Δc⁻⁶ = Δc_L⁻⁶ + 48Gm(W²/(w_a² + w_b²))a·Λ(ε)(τ − τ_L): the distance falls as τ^{−1/6} and the residual as τ^{−2/3}. <br>• **A collapsing cell.** Target 20 at 10⁻¹⁰, run to τ = 10⁷: the distance exponent goes −0.268 (τ ≈ 110) → −0.194 (400) → −0.1668; the residual's → −0.6666; that of W − A → −0.332 [own]. <br>• **A "frozen" cell.** Target 9 at 10⁻¹⁰ collapses later: its exponent is 0.0000 to τ ≈ 2·10⁴, then −0.026 (2.7·10⁷), −0.10 (2.2·10⁸) and −0.162 (5·10⁹) [own]. The formula predicts a relative change of 7.8·10⁻⁶ over τ = 200 → 6400; I observe 8.0·10⁻⁶. It predicts Δc = 2.68·10⁻⁵ at τ = 1.3·10⁸; I observe 2.79·10⁻⁵ at 7.7·10⁷ and 2.56·10⁻⁵ at 2.2·10⁸. rtol 10⁻¹¹ / atol 10⁻¹⁴ gives the same. <br>• The frozen cells simply landed at p ≈ 0.009–0.015. | (c): state the law and its −1/6 and −2/3 exponents. Drop "frozen", or define it by Δc_L. Fit the exponents from τ_L. |
| D3bo-5 | minor | **Mutants.** 73 mutants: 50 killed (16 by the tests only, 21 by `--from-raw` only, 13 by both) and 23 survive [own]. <br>• Equivalent on the stored data or at the printed precision (9): P2, S2, L5, O2, O7, H5, C2, R4, R5. E1 is probably equivalent. <br>• Not equivalent (13): <br>  – F2: `pair_fits` sign-blind in W (a pair of the wrong sign fits); <br>  – F4: W = 0 fits; <br>  – F5: midpoint instead of centre of mass; <br>  – I1: `pair_index` without \|Δlog σ\| (D's d = 0 widths 0.3/0.06 would pair); <br>  – H1: `neg_path` counting x < 0 (the test's "−0.005 is not one" never decides the maximum); <br>  – H7: the threshold at 1% of the largest amplitude; <br>  – C1: nearest feature by centre only; <br>  – R7: H2's last-block slice; <br>  – R8: hl_job's formed threshold; <br>  – R9: `last_tau` taken from the first event; <br>  – R10: H3's second amplitude 0.9; <br>  – E2: oracle error against the data; <br>  – E3: oracle fitted to f + ½ξ. Both E2 and E3 pass the oracle test's 10⁻⁵ < e < 10⁻³. <br>`test_the_loop_phases` tests only the "lands at τ₃" scenario, which D3bo-2 shows is the exception. | tests: `pair_fits` with W of the wrong sign, W = 0 and an off-midpoint pair; `pair_index` at one centre with two widths; `h_stats` with a block whose only negative is −0.5%; h_job, hl_job and h3_targets on a cheap case; the oracle against a direct fit of f + ξ; `phases` with an intermediate phase. |
| D3bo-6 | minor | **(i) "from any start whose coefficient has the feature's sign" (D59.md:13).** The proof uses E(start) < ½‖f*‖², which a seed coefficient gives. At (c*, σ*), ∂_σE = ½π^{1/2}w(w − A), which is > 0 for w > A > 0. A large starting coefficient therefore drives σ down, and the compactness step fails. | Restrict to E(start) < ½‖f*‖², for instance a seed. |
| D3bo-7 | minor | **Literature and New (D59.md:93–111).** <br>• The algebraic convergence of a fit with one component too many is a known class [memory]: over-parameterized student–teacher (Xu–Du, COLT 2023: Θ(T⁻³) for one neuron) and over-specified Gaussian mixtures under EM (Dwivedi et al., Ann. Statist. 2020). <br>• The cancellation clause is contradicted (D3bo-1). <br>• I could not verify the Kranz–Gallon–Dereich–Jentzen (NeurIPS 2025) entry. | Cite these. Claim the (−1/6, −2/3) collapse law only once it is derived (D3bo-4). Narrow New to the ε^{−1/6} landing law with τ₃ from the path. |
| D3bo-8 | minor | C's table still prints the pre-block state of the budget stops: formed kernels [1, 4] and maximum residual 2.7·10⁻¹ [md]. The six dipole ends have one formed kernel and rms 0.16–0.27 [raw] (left over from D3bh-8). | Print the stalled block in C's table. |
| D3bo-9 | cosmetic | <br>• "An idle null cluster" (D59.md:146): the triple's coefficients are 5–11% of the amplitude, so they are formed, not idle [raw]. <br>• "The five cells with a phase" (:136) is the grid's doing: on a fine grid all six have one, target 29 included with 22 slopes, median 1.0005 [own]. <br>• The 72 split rows contain 10 repeats [raw]. | Edit as stated. |

## Sources

I could not consult any source:
- this turn's web-search budget was exhausted on my first query;
- papers.neurips.cc does not resolve from this session (it would have to be allowed in the session's network settings);
- arxiv.org and ams.org are unreachable, as the brief says.

Everything below is therefore [memory] and unverified:
- **Gallon–Jentzen–Lindner 2022**, "Blow up phenomena for gradient descent optimization methods in the training of artificial neural networks" (arXiv 2211.15641): GF risk converging to the infimum while the parameters diverge. Fits the citation.
- **Kranz–Gallon–Dereich–Jentzen, NeurIPS 2025**: not verified, title unknown to me.
- **Chizat**, Math. Program. 194 (2022), "Sparse optimization on measures with over-parameterized gradient descent": conic particle GD with a multiplicative (Fisher–Rao) weight mobility. Fits.
- **Pearson**, Phil. Trans. A 185 (1894): moments of a two-normal mixture. Fits.
- **Runnalls**, IEEE TAES 43 (2007): the moment-preserving merge. Fits.
- **Schiebinger–Robeva–Recht**, Inf. Inference 7 (2018), "Superresolution without separation". Fits.
- **Denoyelle–Duval–Peyré–Soubies**, Inverse Problems 36 (2020), "The sliding Frank–Wolfe algorithm". Fits.
- **Suggested additions**: Xu–Du, COLT 2023; Dwivedi–Ho–Khamaru–Wainwright–Jordan–Yu, Ann. Statist. 48 (2020).

## Lesson

A conclusion of the form "X is the path's", "the pair lands" or "the cell is frozen" needs the control that could have made it false. The revision compared the loop with the path only on the symmetric pair. That is the one case where the path's own flow, with its polish off, disagrees with the path. It timed the landing on a grid coarser than the event, and it read "frozen" from a 6400-unit window. Running the controls overturned three conclusions:
- the polish-free path and the loop on the asymmetric pairs both keep the cancellation;
- the landing follows a later formation and rescales v by an O(1) factor;
- the "frozen" cell collapses by the same τ^{−1/6} law once the horizon is long enough.

The same controls strengthened the main new result: the dipole phase obeys the block's law to 10⁻³, τ₃ follows from the path alone, and v ∝ ε^{−1/6} holds down to 10⁻¹².
