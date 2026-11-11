# Referee report D3bu, round 74: Proposition D‴.59 as revised after round 68

## Scope

- I reviewed commit 18bcfa0, exported with `git archive` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref74/D3bu_work/repo`. The files were:
  - `review/drafts/D59.md` (202 lines);
  - `experiments/exp_moving_width.py`, `results/chaos/moving_width.md` and `results/chaos/moving_width_raw.npz`;
  - `tests/test_moving_width.py` (26 tests);
  - `src/core/moving_width.py`, and in `src/core/rich_transit_nd.py` the functions `sliding_iss_nd`, `loop_nd`, `path_duals_nd`, `critical_rhs_nd`/`critical_jac_nd` and `polish_block_nd`;
  - `review/REVIEW_round68_D3bo.md` and the Round 68 entry of `review/REVIEW_LOG.md`.
- When I finished, HEAD was still 18bcfa0 and none of these files had changed.
- I modified nothing in the repository. Everything I wrote is under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref74/D3bu_work/`:
  - scripts in `scripts/`, logs in `logs/` and loop states in `out/`;
  - the `--from-raw` output in `mwout/`;
  - the mutation workspaces in `mut/`.
- The experiment writes to `$MW_OUT`, or to `results/chaos` if that is unset. I always set `MW_OUT` to `mwout/` (or `mut/out*` for the mutants), each holding a copy of the raw.
- Every run used the venv python with `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and `nice -n 10`. I kept at most two heavy processes of mine at a time; other sessions held the load average at 5–10. I signalled no process, and nothing of mine is still running.
- Tags: [md] the committed md, [raw] the committed raw, [own] my run or my computation on the raw.

## Summary

**What checks out**
- **Reproducibility** [own].
  - `--from-raw` regenerates the md byte for byte.
  - The 26 tests pass in 11.9 s.
  - These reruns are bit-identical to the raw:
    - G's fine-phase cells at ε = 10⁻¹² for all six targets. Every field of `fine_phases` matches (τ₃, τ₄, v(τ₃⁻), the intermediate phase, the landing, the agreement counts) and so do the stored slope points.
    - HL "σ* 0.06, σ₂ 0.065" at 10⁻⁸: max |Δw| = max |Δθ| = 0 at all four times.
    - HP "σ* 0.06, c₂ +0.01": every field.
    - GL: every row.
    - `t3_from_duals` on the six paths: 29.161347, 34.227368, 58.436588, 36.168160, 47.188430 and 29.627607, which are GT's values.
- **Numbers.** Every number in D59.md matches the md or the raw. I found no misquote. Where I flag a number below, its meaning is wrong, not its transcription.
- **τ₃ from the path.**
  - `t3_from_duals` is right. It advances the duals from the last history time to τ_d with that interval's residual, then moves them by the block's residual at s = 10¹¹.
  - All 8 of my mutants of `t3_from_duals` are killed.
  - The limit argument holds. Its rate follows (D3bu-7).
- **(a)'s landing algebra re-derives.** |W − A|/|A| ≈ p²/8 comes from the second moment and offset ≈ p³σ/24 from the third, so the 1% landing is p ≤ √0.08 = 0.283. Checked on the loop's states [own from raw]:
  - target 18 at 10⁻¹⁰, τ = 400: p = 0.276 and (W − A)/|A| = 0.0094, against p²/8 = 0.0095;
  - target 20: p = 0.156 and 0.0030 against 0.0030; offset 1.58·10⁻⁴ against p³/24 = 1.57·10⁻⁴.
- **The dipole phase in the loop fits better than D59 shows.** Evaluated at the time of v(τ₃⁻), the law gives 1.0006, 1.0003, 1.0130, 0.9995, 1.0003 and 1.0001 at 10⁻¹² for targets 9, 18, 12, 20, 27 and 29 [own from raw]. D59 quotes 0.9886–1.0120 [md].
- **K's third-order formula.** K = −⟨r_dip, ∇³φ:p³⟩/24 reproduces G's K = v²(E − E_inf) to within 5·10⁻⁴ on all six blocks, and to 2·10⁻⁵ on four [own].
- **The collapse law re-derives.** E ≈ aΔc⁸; the separation's mobility GmW²/(w_a² + w_b²) is the paradigm's metric induced on dc̄ = 0; the exponents are −1/6, −1/3 and −2/3. Its constant holds as p → 0 [own]:
  - The loop's state sits at the minimum of E over everything but (Δc, v): to 5.3·10⁻⁵ for target 20 and 8.2·10⁻⁴ for target 9.
  - v is frozen: it changes by 10⁻⁵ over τ = 400 → 10⁷.
  - Target 9 is a "not yet collapsing" cell (landed at p = 0.015). It collapses at 0.997–0.999 of the law's rate over τ = 735 … 10¹⁰.
  - Target 20 reaches 0.999 of the law's rate at τ = 10⁷.
- **(i).** These re-derive: J, its single critical point, its limits, ∂_σE = ½π^{1/2}w(w − A) at (c*, σ*), and the proof under E(start) < ½‖f*‖².
- **C's table** now counts the stalled blocks: formed kernels 2–4, the largest residual 0.12 and 0.14 [md]. D3bo-8 is fixed.
- **H2's stop reasons** establish that its 10⁻¹³ ends come from the exhausted dictionary, the budget and the capped polishes.
- **(iii), E and F** are unchanged and hold as stated.

**What is wrong (majors first)**
1. **The intermediate phase is a second dipole phase, and κ follows from the path (D3bu-1, major).**
   - Between τ₃ and τ₄, v³ grows linearly in τ with slope 3μ₀K′Λ(ε). K′ is the three-kernel block's dipole constant: K′/K = +0.046, −0.083, −0.356, +0.228, −0.090 [own].
   - From the path alone, κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3}. This uses the closure fit of the block plus the τ₃ entrant, its K′, and the duals' τ₄. It gives 1.1145, 0.7519, 0.9116, 1.0195 and 0.7690, against the loop's 1.1159, 0.7525, 0.9148, 1.0209 and 0.7694 [own].
   - D59 calls the phase and κ "numerical", and κ "an O(1) factor of the target, ε-free to 0.4%". Over the md's own cells κ = 0.729–2.905 [md]; it only converges as ε → 0.
2. **(ii)(b)'s evidence from the loop and the polish-free path counts dead spikes (D3bu-2, major).**
   - At 10⁻⁸ and τ = 10⁵ the loop keeps 11–12 negative kernels. Of these, 4–6 are spikes invisible on the samples (largest φ ≤ 1.4·10⁻¹¹), and they are exactly D59's "widths below the samples' spacing". The visible kernels are 6–7 negatives, all at least 0.036 wide [own from raw].
   - HP's 12 and 11 negatives include 6 and 5 such spikes [raw].
   - The symmetric loop also keeps 8–12 visible negatives: 8 at 10⁻⁶ out to τ = 10⁷, with |w| growing to 0.839 [own]. Only the inflated cluster is "the polish's".
3. **Minor and cosmetic issues:**
   - the collapse law's constant fails at the landing: the slaved coefficients' drag is 10% at p = 0.15 and ≈44% at p = 0.28; a is undefined; τ_L is misplaced (D3bu-3);
   - the agreement statistic for p ≲ 0.28 (D3bu-4);
   - the argument for (i)'s "a large one is not enough" (D3bu-5);
   - 13 behaviour-changing surviving mutants (D3bu-6);
   - the τ₃ proof does not state the residual's limit, and the dipole-phase check has a grid offset (D3bu-7);
   - literature and New (D3bu-8);
   - cosmetic (D3bu-9).

## What I checked and how
1. **Reproducibility.**
   - `--from-raw` with `MW_OUT=mwout/`, compared with `cmp`.
   - `pytest -p no:cacheprovider tests/test_moving_width.py`.
   - `loop_states.py`: G's loop on `union(TAUS_G, TAUS_FINE)` at 10⁻¹² for each target, then `fine_phases`, compared field by field (11–32 s per loop).
   - `rerun_hl.py`: `hl_job` with `EPS_HL = (1e-8,)`.
   - `rerun_hp.py`: `h_job(…, None, False)`.
   - `gl_constant.py`: GL's loop, rows compared.
   - `t3_from_duals` on the six paths.
2. **Numbers.** Every number in D59.md:1–202 against the md and the raw, with `raw_overview.py`, `raw_fine.py`, `inter_phase.py`, `visible.py` and inline checks.
3. **Theory, re-derived by hand.**
   - (i)'s ∂_σE.
   - (a)'s second and third moments, giving p²/8 and p³σ/24.
   - The dipole law and K's third-order formula.
   - The separation's mobility on dc̄ = 0, and the slaved coefficients' share of that metric.
   - The duals' convergence rate.
4. **The four stages.**
   - `inter_phase.py`: G's coarse states at every ε (v, p, offset, (W − A)/|A|, slopes of v³).
   - `kappa_pred.py`: the 10⁻¹² fine states (v³ on (τ₃, τ₄), and K′ from `dipole_law`).
   - `t4_pred.py` and `kappa_path.py`: K, K′, τ₃, τ₄ and κ from the path alone, by a least-squares closure fit and the third-order formula.
   - `agree_phase.py`: where p ≤ 0.28 and the landing disagree, sorted by phase.
5. **The collapse.**
   - `gl_constant.py`, on target 20 to τ = 10⁷ and target 9 to 10¹⁰, both at 10⁻¹⁰. It compares E at the loop's state with its minimum at fixed (Δc, v), computes a, and compares the measured d(Δc⁻⁶)/dτ with the law.
   - `gl_metric.py`: the rate against the metric induced on every moving parameter.
6. **(ii)(b).**
   - `visible.py`: the largest φ on the samples of every formed kernel in H, HP, H2 and HL.
   - `hl_long.py`: "A₂ 0.99" at 10⁻⁸ and the equal pair at 10⁻⁶, both to τ = 10⁷.
7. **(i).** `one_kernel_large_w.py`: one kernel's flow on the critical curve from (c*, σ*) with w = 1.5A … 10A, at σ*/σ₀ = 0.25, 1 and 4.
8. **Mutants.**
   - `mutants.py` and `run_mutants.py`: 71 single-point mutants of `t3_from_duals`, `gt_job`, `part_GT`, `fine_phases`, `P_LAND`, `g_fine_lines`, `g_phase_summary`, `moving_state`, `stop_anatomy`, `hl_job`, `h_job`, `gl_job` and `gl_lines`. Each ran against the test file and against `--from-raw` (md byte diff).
   - `fp_equiv.py`: the surviving `fine_phases` mutants applied to the real 10⁻¹² states.
   - `run_mutants68.py`: round 68's 13 non-equivalent survivors, run against the current code.
9. **Literature.** Four web searches; the results are in Sources.

## Verdicts per part

| claim (D59.md line) | verdict | evidence |
|---|---|---|
| (i) J, its single critical point, its limits (:7–12, :75–77) | holds | re-derived; B [md] |
| (i) convergence from E(start) < ½‖f*‖² (:12–14, :78–83) | holds | re-derived |
| (i) "a large one is not enough: … ∂_σE > 0 … drives the width down" (:14–15) | the formula is right; the inference does not follow | w ≤ 5A converges exactly; w = 10A ends as a spike [own] (D3bu-5) |
| (ii) C's outcomes 71 / 3 / 6 (:19–21) | holds | [md] |
| (ii)(a) expansion, split, validity, k-th moment (:22–27, :83–89) | holds | re-derived; the landed pairs' split 0.9828–1.0458 [md] |
| (ii)(b) the path's ends (:28–31, :174–178) | the numbers hold; "widths below the samples' spacing and of negative sign among them" are idle or dead kernels | [own from raw] (D3bu-2) |
| (ii)(b) H2: the dictionary's and the polish's (:31–34, :179–184) | holds | stop reasons [md] |
| (ii)(b) separation window (:34–35, :186–188) | holds | [md] |
| (ii)(b) asymmetric pairs: "the path's and the loop's" (:35–39, :184–193) | holds qualitatively; the counts include 4–6 invisible spikes, which are the "sub-sample widths" | reruns bit-identical; visible counts [own] (D3bu-2) |
| (ii)(b) "\|w\| … still growing" (:37, :192) | holds; run-off or convergence undecided | 2.38 → 2.55 by τ = 10⁷ [own] |
| (ii)(b) symmetric pair: "the polish's" (:39–41) | partly | the loop keeps 8–12 visible negatives, 8 at τ = 10⁷ [own] (D3bu-2) |
| (ii)(c) closure and the block's law v³ = 3μ₀Ks (:42–48, :89–98) | holds | re-derived; K's formula matches G's K [own] |
| (ii)(c)(1) dipole phase, τ₃ from the duals, ε^{−1/6} (:49–55, :98–103, :165–168) | holds, better than stated | 0.9995–1.0006 at matched τ (target 12: 1.013) [own]; GT reproduced [own] (D3bu-7) |
| (ii)(c)(2) intermediate phase (:55–57) | the numbers hold; the mechanism is missing | it is a second dipole phase with K′ [own] (D3bu-1) |
| (ii)(c)(3) landing at p ≲ 0.28, κ (:57–60, :169–172) | the threshold holds in (a)'s regime; κ is derivable and is not ε-free | κ from the path [own]; κ = 0.729–2.905 over all cells [md] (D3bu-1, D3bu-4) |
| (ii)(c)(4) collapse law (:60–66, :104–107, :172–174) | exponents hold; the constant holds only as p → 0; τ_L is misplaced | 0.903 → 0.999 [own] (D3bu-3) |
| (ii)(c) "a pair landed at a small p collapses the same way, later" (:66) | holds | target 9 follows the law's constant to 10¹⁰ [own] |
| (ii)(c) at larger masses the stages overlap (:66–67) | holds | [md] |
| (iii), E, F (:69–73, :193–202) | holds | [md]; unchanged |
| Status (:109–116) | understates (c), overstates (b) | D3bu-1, D3bu-2, D3bu-4 |
| Classical and New (:116–137) | the citations fit; New should be narrowed | D3bu-8 |

## Status of the round-68 issues

| D3bo | round-68 issue | status | what remains |
|---|---|---|---|
| 1 | The cancellation's attribution (major) | partly fixed | **Fixed:** HL on the five asymmetric targets to τ = 10⁵; HP (no polish); the wording. <br>**Remaining:** <br>• HL's and HP's counts include invisible spikes, and these are the "sub-sample widths". <br>• The symmetric loop keeps both signs. <br>• `sliding_iss_nd` still accepts polishes that stop at max_nfev. <br>• Run-off or convergence is undecided (D3bu-2). |
| 2 | The landing (major) | mostly fixed | **Fixed:** four stages; a fine grid down to 10⁻¹²; τ₃ from the duals (within 0.073); κ measured; landing as p ≲ 0.28. <br>**Remaining:** <br>• the intermediate phase and κ are derivable (D3bu-1); <br>• κ is ε-free only in the limit (D3bu-1); <br>• the agreement statistic (D3bu-4); <br>• the 0.1 grid offset (D3bu-7). |
| 3 | H2's artefacts (major) | mostly fixed | **Fixed:** H2 prints its stop reasons; the wording. <br>**Remaining:** <br>• invisible kernels are still counted as formed in H, HP and HL (D3bu-2); <br>• `polish_block_nd`'s max_nfev = 5000 is unchanged. |
| 4 | The collapse law (minor) | fixed, with gaps | **Fixed:** the law and GL. <br>**Remaining:** <br>• a is undefined; <br>• GL checks exponents only; <br>• the metric omits the coefficients; <br>• τ_L (D3bu-3). |
| 5 | Mutants; untested job code (minor) | partly fixed | 11 of round 68's 13 non-equivalent survivors are now killed [own]. <br>C1 survives; R8 survives and is behaviour-changing. <br>There are 13 new behaviour-changing survivors (D3bu-6). |
| 6 | (i)'s start (minor) | fixed | The added justification is a non sequitur (D3bu-5). |
| 7 | Literature (minor) | fixed | Ho–Nguyen is missing; New is too broad (D3bu-8). |
| 8 | C's table (minor) | fixed | — (test added; mutants MS2 and MS3 killed) |
| 9 | Cosmetic | partly fixed | "an idle null cluster" (:178) is still there; the repeated split rows remain (D3bu-9). |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bu-1 | major | **The intermediate phase is a second dipole phase; τ₄, κ and v_land follow from the path (D59.md:55–60, :103, :116, :134–135, :169–171).** <br>• **v³ is linear on (τ₃, τ₄).** At 10⁻¹² a straight line fits v³ to within 7.7·10⁻⁶ (relative), with slope 3μ₀K′Λ(ε) [own]. <br>  – K′/K = +0.0461, −0.0833, −0.3562, +0.2278 and −0.0904 for targets 9, 18, 12, 27 and 29. <br>  – `dipole_law` at an intermediate state gives the same values to four digits. <br>  – So the pair grows (targets 9, 27) or retracts (18, 12, 29) by the block's own law. <br>  – After τ₄, v freezes: it changes by less than 10⁻⁴ out to τ = 6400 [own from raw]. <br>• **From the path alone** (`kappa_path.py`): add the τ₃ entrant to the block's state at s = 10¹¹; fit the closure model (other kernels + Wφ(θ̄) + p·∇_θφ(θ̄)); set K′ = −⟨r′_dip, ∇³φ:p′³⟩/24; move the duals by r′_dip [own]. <br>  – τ₄ = 124.31, 97.89, 81.64, 51.29 and 46.97, against the loop's 124.2, 97.9, 81.7, 51.3 and 47.0 [raw]. <br>  – κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3} = 1.1145, 0.7519, 0.9116, 1.0195 and 0.7690, against the loop's 1.1159, 0.7525, 0.9148, 1.0209 and 0.7694 [raw]. <br>  – Hence v_land³ = 3μ₀Λ(ε)[K(τ₃ − τ_d) + K′(τ₄ − τ₃)]. This is a law, not "an O(1) factor of the target". <br>• **"ε-free to 0.4%" holds only between 10⁻¹¹ and 10⁻¹².** <br>  – Over all the md's cells κ = 0.729–2.905 [md]. Target 27 has 1.758, 1.417, 1.052 and 1.021 at 10⁻⁷ … 10⁻¹⁰; target 18 has 2.905 at 10⁻⁷. <br>  – κ inherits the ε^{1/3} convergence of τ₃ and τ₄ (D3bu-7). | **(c)(2)–(3), the proof (:98–104), Status (:116), New (:134–135):** <br>• state the intermediate phase as the block's law with the three-kernel block's K′ (either sign); <br>• τ₄ from the duals; <br>• κ³ = 1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)); <br>• "κ converges as ε → 0", not "ε-free". <br>**G/GT:** print K′, the predicted τ₄ and κ beside the loop's, with a `kappa_from_path` next to `t3_from_duals`. |
| D3bu-2 | major | **(ii)(b)'s evidence counts dead spikes, and "the polish's" overstates (D59.md:28–41, :175–177, :184–193; Status :112–113).** <br>• **The loop.** At 10⁻⁸ and τ = 10⁵ the five asymmetric targets have 12, 12, 12, 11 and 11 negative kernels [md]. <br>  – Of these, 6, 6, 6, 5 and 4 are spikes invisible on the samples (largest φ ≤ 1.4·10⁻¹¹) [own from raw]. <br>  – These spikes are all of the "widths below the samples' spacing". <br>  – Visible negatives: 6, 6, 6, 6 and 7. Every visible formed kernel is at least 0.036 wide [own from raw]. <br>  – At 10⁻⁶, 2–3 of the 4–6 negatives are invisible, leaving 2–3 visible [own from raw]. <br>• **HP.** 6 of A₂ 0.9's 12 negatives and 5 of c₂ +0.01's 11 are invisible; HP's own raw field `invisible` reads 6 and 5. That leaves 6 visible in each [raw]. <br>  – None of HP's four paths reaches the tolerance: three exhaust the dictionary (21 kernels active) and A₂ 0.9 exhausts its budget. <br>  – Their residuals are 7.5·10⁻⁷ … 9.7·10⁻⁶ [raw]. <br>• **The path (H).** <br>  – The 3–6 sub-sample and 0–2 negative-width kernels at its ends (:30, :176) are idle: \|w\| below 1% of the amplitude and largest φ ≤ 5.3·10⁻⁶. <br>  – The one exception is c₂'s spike, 1.75·10⁻⁴ wide (w = −0.018, largest φ 6·10⁻⁶⁹). <br>  – The formed kernels are at least 0.049 wide [own from raw]. <br>  – So H, HP and the loop at 10⁻⁸ all keep 6–7 visible negatives. This supports the attribution, but "11–12 … widths below the samples' spacing" does not describe the fit. <br>• **The symmetric pair.** <br>  – At τ = 10⁵ the loop keeps 8 visible negatives at 10⁻⁶ and 12 at 10⁻⁸ [raw]. <br>  – At 10⁻⁶, run to τ = 10⁷, it still has 8. Max\|w\| goes 0.705 → 0.782 → 0.839 at τ = 10⁵, 10⁶ and 10⁷, and the residual 2.5·10⁻⁵ → 2.3·10⁻⁶ [own]. <br>  – So kernels of both signs are the loop's here too. What the polish adds is the inflated 1 + 8 cluster (\|w\| 3.24). <br>• **"Still growing".** A₂ 0.99 at 10⁻⁸, run to τ = 10⁷ [own]: <br>  – max\|w\| goes 2.381 → 2.476 → 2.554, with per-decade increments of 0.095 and 0.078; <br>  – the residual falls as τ^{−0.59}; the 6 visible negatives are unchanged. <br>  – Run-off or convergence is undecided. If the increments keep shrinking geometrically, the limit is ≈ 2.9. | **(b) and Numerics (H, HP, HL):** <br>• count only kernels visible on the samples (largest φ ≥ 10⁻¹⁰): apply D3bo-3's filter in `h_stats`, `hl_job` and the HL/HP tables, not only in H2's line; <br>• drop "widths below the samples' spacing" from the description of the fit; <br>• say the symmetric loop also fits with kernels of both signs at \|w\| < A; <br>• give HP's stop reasons; <br>• add the τ = 10⁷ run. |
| D3bu-3 | minor | **The collapse law's constant holds only as p → 0; a is undefined; τ_L is misplaced (D59.md:60–66, :104–107, :172–174).** <br>• **a** is not defined in D59. It can be computed as E_min/Δc⁸ at fixed (Δc, v). <br>  – For target 20 at 10⁻¹⁰ it goes 7.89·10⁶ → 6.99·10⁶ over τ = 400 → 10⁷ [own]. <br>  – Across ε it scales as v⁸/W⁶: a·W⁶/v⁸ = 0.64, 0.66, 0.69, 0.77 (target 20) and 0.956, 0.951, 0.945, 0.983 (target 12) at 10⁻¹⁰ … 10⁻⁷ [own from raw]. <br>• **The constant.** For target 20, the ratio of d(Δc⁻⁶)/dτ to 48aGmW²/(w_a² + w_b²)Λ is 0.903, 0.961, 0.987 and 0.999 at τ = 574, 7.2·10³, 1.3·10⁵ and 10⁷ (p = 0.151, 0.097, 0.060, 0.029) [own]. GL shows only the exponents. <br>• **Why the gap.** <br>  – The slaved W − A ∝ Δc² moves at the coefficients' mobility μ₀\|w\| = 0.07–0.48, far below the centres' Gm = 12. <br>  – That motion adds (Gm/(16μ₀))W²p²/(vσ²) to the metric the proof uses: 0.129 at p = 0.151, against a measured 0.127. <br>  – With every moving parameter in the metric, the predicted rate is 0.992–0.994 of the measured rate at every τ [own]. <br>  – At the landing's p = 0.28 the correction is ≈ 0.44 for target 20 at 10⁻¹⁰. It vanishes as p → 0 or v → ∞. <br>• **τ_L** must lie past the fast relaxation after τ₄. <br>  – Target 9 at 10⁻⁷ meets the 1% landing at τ = 128.8 with p = 0.137 [md]. <br>  – At τ = 132 its p = 0.016, unchanged to τ = 6400 [own from raw]. <br>• **The small-p claim (:66)** rested on round 68's run. It holds, with the constant: target 9 at 10⁻¹⁰ collapses at 0.9966–0.9988 of the law over τ = 735 … 10¹⁰. Its Δc goes 3.006 → 1.459·10⁻⁵, with local exponent −0.164 at the end [own]. | **(c)(4) and the proof:** <br>• define a; <br>• state the law as the p → 0 limit, or add the coefficients' term to the metric; <br>• take τ_L after the relaxation that follows τ₄. <br>**GL:** report the constant; add target 9. |
| D3bu-4 | minor | **"landed ⇔ p ≲ 0.28" and its statistic (D59.md:58–59, :103, :116, :171–172).** <br>• The equivalence holds for a pair on a feature in (a)'s regime. <br>• The 36535 − 34067 = 2468 disagreements [md] break down as follows [own from raw]: <br>  – 1548 are at ε ≥ 10⁻⁶, where the stages overlap. <br>  – 719 are target 27's at ε ≤ 10⁻⁷, almost all in its σ-dipole phase (\|p_σ\|/\|p\| = 0.994 [md]). There p built from Δc alone is 0.27, although the pair is not a single Gaussian; with \|Δθ\| = (Δc² + Δσ²)^{1/2} it is 2.43. <br>• With \|Δθ\| in p, the coarse grid's disagreements drop from 96 to 15 out of 1291 [own]. <br>• :103 and :116 call the threshold "numerical", but it is derived. | Use \|Δθ\| in p, or count only states past τ₄. Say which phases the statistic covers. Call the threshold derived. |
| D3bu-5 | minor | **(i)'s "a large one is not enough" (D59.md:14–15).** <br>• ∂_σE = ½π^{1/2}w(w − A) re-derives. But it is also positive for A < w < 2A, where E(start) < ½‖f*‖² and the theorem applies. So it shows a transient, not a failure. <br>• One kernel's flow on the critical curve from (c*, σ*), at σ*/σ₀ = 0.25, 1 and 4 [own]: <br>  – w = 2.5A, 3A and 5A (E(start) = 2.25, 4 and 16 times ½‖f*‖²) converge to the exact fit, as do starts at 0.5σ*, 2σ* and c* + σ*: w/A = 1, σ/σ* = 1, E ≈ 10⁻³². On the way the width dips as low as 0.006σ*. <br>  – w = 10A ends as a spike: σ = 0.005–0.08σ*, w ≈ 10A, and E = 0.64–0.98 of ½‖f*‖² at s = 10⁹. <br>• So the energy condition is sufficient but not sharp. | **(i):** "the proof needs E(start) < ½‖f*‖²; from (c*, σ*) the flow converges up to w = 5A and, on the samples, shrinks to a spike at 10A". |
| D3bu-6 | minor | **Mutants: 29 of 71 survive, 13 of them behaviour-changing [own].** The tests kill 22; `--from-raw` alone kills 20 more. <br>• **Equivalent on the data or below the printed precision (15):** <br>  – GT2: no exits occur before τ_d. <br>  – PG2. <br>  – FP2, FP4, FP5, FP9 and FP10: identical output on the six 10⁻¹² cells (`fp_equiv.py`). <br>  – GS2, MS1, MS4, HL1. <br>  – HL2: `smin` is not printed. <br>  – HJ4, GJ2 (probably), GL5. <br>• **Undecidable: SA2.** The raw does not store polish statuses. <br>• **Behaviour-changing (13):** <br>  – `part_GT` PG1 (`status != 0`): GT comes out empty, and `g_fine_lines` then silently prints G's stale t3_pred (D3bu-9). <br>  – `fine_phases` FP6 (p from σ_a only), FP8 (p < 2), FP11 (`dip[::5]`) and FP13 (\|dW\|). `--from-raw` cannot see `fine_phases`, and its one synthetic test misses these. <br>  – `stop_anatomy` SA4 (invisible means φ below 10⁻⁵): H's "equal" end has active kernels with largest φ 3·10⁻⁷ … 5·10⁻⁶. <br>  – `h_job` HJ1 and HJ2 (HP runs with the polish), HJ5 (H2 runs at 10⁻⁹) and HJ6 (clusters divided by the swapped amplitudes). <br>  – `gl_job` GJ1, GJ3 and GJ4: no test runs `gl_job`. <br>• **Round 68's 13 non-equivalent survivors:** <br>  – 11 are now killed: F2, F4, F5, I1, H1, H7, R7, R9, R10, E2, E3. <br>  – C1 survives. <br>  – R8 (HL's formed threshold at 0.1%) survives and is behaviour-changing: at τ = 6400 the equal pair at 10⁻⁶ has 2 kernels between 0.1% and 1% [own]. | **Tests to add:** <br>• `fine_phases` with a σ-dipole, a state with 0.28 < p < 1, the slope sampling, and the sign of dW; <br>• `part_GT` on a status-0 entry; <br>• `h_job` with `(…, None, False)` and with a tolerance, checking (with `sliding_iss_nd` monkeypatched) that `polish` and `tol` reach it; <br>• `clusters` through `h_job` with unequal amplitudes; <br>• `gl_job` on a short horizon against `loop_nd`; <br>• `stop_anatomy` with a kernel of largest φ 10⁻⁶; <br>• `hl_job` with a kernel between 0.1% and 1%. |
| D3bu-7 | minor | **The τ₃ proof and the dipole phase's check (D59.md:100–103, :167–168).** <br>• **The state diverges, the residual converges.** <br>  – The block's state has no limit, since \|w\| ∝ s^{1/3}. Only its residual converges: r(s) = r_dip + O(v⁻²) = r_dip + O(s^{−2/3}), given the numerically observed approach to the infimum. <br>  – The duals' error is then (1/Λ)∫₀^{Λ(τ−τ_d)}O(s^{−2/3})ds = O(ε^{1/3}). <br>  – Target 12's shortfall in τ₃ is 3.64, 1.94, 1.00, 0.44 and 0.24 at 10⁻⁷ … 10⁻¹¹. It falls by about 2 per decade, against 10^{1/3} = 2.15 [own from raw]. <br>• **Grid offset.** <br>  – v(τ₃⁻) is v at τ₃ − 0.1, but the law is evaluated at τ₃. <br>  – The resulting bias, (1 − 0.1/(τ₃ − τ_d))^{1/3}, is −1.1% for target 29 and −0.9% for target 20. <br>  – At matched τ the ratios are 0.9995–1.0006 (target 12: 1.013) [own]. κ carries the same shift. | **Proof:** "its residual converges (O(s^{−2/3})) though its state runs off", plus the rate. <br>**`fine_phases`/`g_fine_lines`:** evaluate the law at the time of v. |
| D3bu-8 | minor | **Literature and New (D59.md:116–137).** <br>• (a)'s O(Δc⁴) result, and hence E ∝ Δc⁸, is Ho–Nguyen's r̄ = 4 for a location–scale Gaussian mixture over-fitted by one component (Ann. Statist. 44 (2016) 2726–2755). The venue rests on [snippet] and r̄ = 4 on [memory]. This is the moment algebra the draft credits to Pearson. <br>• The collapse exponents are the generic algebraic rate of a gradient flow at a degenerate minimum: E ∝ x⁸ gives x ∝ t^{−1/6}. Xu–Du's Θ(T⁻³) lower bound uses a potential on the students' pairwise distances [snippet], the same mechanism. <br>• What is specific here is the stage sequence (K, τ₃, K′, τ₄ and κ, all from the path) and the law's constant in the paradigm's metric, including the coefficient term. | Cite Ho–Nguyen beside Pearson. Narrow New to the stage sequence and the constant, not the exponents. |
| D3bu-9 | cosmetic | • "an idle null cluster" (:178) is unchanged, although the round-68 log says it was fixed. The triple's −0.0546/+0.1093/−0.0546 are 5–11% of A [raw]. <br>• :163's "72 rows" repeat frozen states: only 60 split ratios are distinct at six decimals [own from raw]. <br>• The raw's G entries still carry the pre-fix t3_pred (39.17, 37.56, 110.54, 40.92, 76.52, 51.14) [raw], which the current `g_job` would not produce. <br>• :103 and :116 call the threshold "numerical". | Edit the text. Drop `t3_pred` from `g_job`, or regenerate G. |

## Sources

Web search answered four queries. projecteuclid.org does not resolve from this session; it would have to be allowed in the session's network settings. arxiv.org and ams.org were unreachable, as the brief says.

- **Ho, N. and Nguyen, X.**, "Convergence rates of parameter estimation for some weakly identifiable finite mixtures", Ann. Statist. 44(6) (2016) 2726–2755. [Project Euclid page](https://projecteuclid.org/journals/annals-of-statistics/volume-44/issue-6/Convergence-rates-of-parameter-estimation-for-some-weakly-identifiable-finite/10.1214/16-AOS1444.full) (search result; the page itself did not load). Title, venue and scope are [snippet]; r̄(1) = 4 is [memory].
- **Xu, W. and Du, S. S.**, "Over-Parameterization Exponentially Slows Down Gradient Descent for Learning a Single Neuron", COLT 2023, [PMLR 195:1155–1198](https://proceedings.mlr.press/v195/xu23a.html). It proves O(T⁻³) with a matching Ω(T⁻³), the lower bound via a potential on pairwise distances [snippet]. D59's "Θ(T⁻³) for one neuron" is right.
- **Gallon, D., Jentzen, A. and Lindner, F.**, "Blow up phenomena for gradient descent optimization methods in the training of artificial neural networks", [arXiv:2211.15641](https://arxiv.org/abs/2211.15641). Gradient-flow trajectories with risk below a threshold diverge [snippet]. Fits the citation.
- **Dwivedi, Ho, Khamaru, Wainwright, Jordan and Yu**, "Singularity, misspecification and the convergence rate of EM", Ann. Statist. 48(6) (2020) 3161–3182, [doi:10.1214/19-AOS1924](https://www.doi.org/10.1214/19-AOS1924). For over-specified fits, EM is exponentially slower with balanced weights [snippet]. Fits the citation.
- **[memory], not verified this round:**
  - de Boor 1969, "On the approximation by γ-polynomials" (in Schoenberg, ed., *Approximations with Special Emphasis on Spline Functions*);
  - Kammler, J. Approx. Theory 1973 and 1976;
  - Braess, *Nonlinear Approximation Theory*, Springer 1986;
  - Petersen–Raslan–Voigtlaender, Found. Comput. Math. 21 (2021);
  - Chizat, Math. Program. 194 (2022);
  - Pearson, Phil. Trans. A 185 (1894);
  - Runnalls, IEEE TAES 43 (2007);
  - Schiebinger–Robeva–Recht, Inf. Inference 7 (2018);
  - Duval–Peyré, Inverse Problems 33 (2017).

  Each fits the way D59.md uses it.

## Lesson

Round 68 asked for controls. This round's findings come from asking two more questions of what those controls produced. First, is a newly described stage really new, or is it the previous stage's law with new data? The intermediate phase was presented as an O(1) factor "of the target". It is in fact the block's dipole law with the three-kernel block's constant, so τ₄, κ and the landing velocity follow from the path just as τ₃ already did, to 0.4%. Second, are the kernels a count counts visible to the fit? Half the negative kernels in (ii)(b)'s loop and polish-free evidence, and all of its sub-sample widths, are spikes invisible on the samples. This is the same artefact round 68 flagged in H2, and the revision filtered it only there. Finally, a derived law should be checked by its constant, not only by its exponents. The collapse law's constant holds as p → 0 but not at the landing, because the slaved coefficients still move, and the metric must charge them for it.
