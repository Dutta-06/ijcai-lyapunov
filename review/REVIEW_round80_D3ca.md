# Referee report D3ca, round 80: Proposition D‴.59 as revised after round 74

## Scope

- **What I reviewed.** Commit 74e57ea, exported with `git archive` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/repo` (called W/repo below; W is the `D3ca_work` directory). The files were:
  - `review/drafts/D59.md` (258 lines);
  - `experiments/exp_moving_width.py`, `results/chaos/moving_width.md` and `results/chaos/moving_width_raw.npz`;
  - `tests/test_moving_width.py` (38 tests);
  - `src/core/moving_width.py`, and in `src/core/rich_transit_nd.py` the functions `_curve_terms`, `critical_rhs_nd`, `critical_jac_nd`, `loop_nd`, `polish_block_nd`, `sliding_iss_nd` and `path_duals_nd`;
  - `review/REVIEW_round74_D3bu.md` and the Round 74 entry of `review/REVIEW_LOG.md`.
- **Settings.** Every run used the venv python with `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and `nice -n 10`.
  - Every run of the experiment had `MW_OUT` set to a directory of mine holding a copy of the raw: W/mwout, or `W/mut/<id>/out` for the mutants.
  - I kept at most two heavy processes of mine at a time.
- **I modified nothing in the repository**, and I ran no git command that writes.
  - During the review, other sessions' commits moved HEAD from 74e57ea to dd88e96 (ed55ab1: D‴.52, round-78 report; dd88e96: D‴.58, round 77). Neither touches the reviewed files: D59.md, the md and the raw are byte-identical to my export.
  - Everything I wrote is under W: scripts in `W/scripts`, logs in `W/logs`, slow-manifold arrays in `W/out`, mutants in `W/mut`.
- **Nothing of mine is still running.** The python processes left on the machine belong to another referee's workspace (`ref79/D3bz_work`, checked via /proc). I signalled no process.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my run, or my computation on the raw; [snippet] a web-search result; [memory] not verified this round.

## Summary

**What checks out**

- **Reproducibility** [own]:
  - `--from-raw` regenerates the md byte for byte (cmp).
  - The 38 tests pass in 9.6 s.
  - Re-running GT through `kappa_from_path` on all six paths gives bit-identical results.
  - GL is bit-identical (both targets, all 29 rows); HLL is bit-identical (both cells); four A2 cells are bit-identical.
- **Numbers.** I checked every number in D59.md that I traced against the md or the raw and found no misquote. The one exception is a rounding: "1.115" stands for 1.1145, which rounds to 1.114.
- **The second phase's law holds for both signs of K′.**
  - Along G's stored loop states at 10⁻¹⁰ in (τ₃, τ₄), a closure fit started from each state gives (E − E′_dip)v² = K′ to ≤ 4·10⁻⁴ relative, on all five targets [own].
  - K′ does not depend on the entrant's seed: the KP7 mutant reproduces it to six digits.
  - τ₃ and τ₄ come from the duals, and κ³ = 1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)) follows. Once the fine grid's offset is removed, the loop's κ at 10⁻¹² lies within −0.27 … +0.21 % of the path's [own].
- **The O(ε^{1/3}) rate holds as a bound.**
  - τ₃'s shortfall falls by ×1.9–2.5 per decade (targets 12, 29, 9) — that is ε^{1/3} [own from raw].
  - τ₄'s falls faster, by ×2.4–4.2 per decade.
  - The bound needs E − E_inf = K/v² + O(v⁻⁴), which is true (D3ca-7).
- **The collapse law and C re-derive.** I did this independently:
  - I computed target 20's slow manifold (E minimised at fixed (Δc, v)), and the induced metric of every coordinate.
  - The coefficients' share is 0.12338 against C's 0.12299 at p = 0.148.
  - The prediction (force factor over metric factor) reproduces GL's measured ratio within 8·10⁻⁴ (3·10⁻⁴ at τ = 10⁷).
  - It reproduces target 9's ratio including C: 0.9989 predicted against 0.9986–0.9990 measured [own].
  - Attributing the remainder to E's next order is right — with the caveats of D3ca-4.
- **(a) holds for signed pairs.** For a two-component pair with w_i = Wx, the fourth-order coefficient is −Wx(1 − x)(1 − x + x²)Δ⁴/36. It vanishes for no real x ∉ {0, 1} [own]. So a pair fits one Gaussian only to O(Δc⁴) for either sign: Ho–Nguyen's r̄ = 4 extends to cancelling pairs.
- **Visibility.** The visible-negative counts in H, H2, HP, HL and HLL do not change for any threshold between 10⁻¹⁰ and 10⁻³ [own from raw].
  - The invisible spikes are 3–6·10⁻⁴ wide and sit between samples at c ≈ 0.22–0.28, with largest φ between 7·10⁻¹⁵ and 2·10⁻¹¹.
  - HP's stop reasons are correct.
  - "Undecided" for the asymmetric clusters is fair. Run to τ = 10⁹, A₂ 0.99 at 10⁻⁸ keeps 6 negatives (Σ −2.30 → −2.11) while max|w| goes 2.38 → 2.67 by shrinking increments; A₂ 0.95 at 10⁻⁶ keeps 3 [own].
- **The rest holds:** (i)'s proof and ∂_σE = ½π^{1/2}w(w − A); the moments behind the landing (p²/8, p³σ/24); K's third-order formula; Ho–Nguyen's venue and r̄(2) = 4 [snippet]; Xu–Du [snippet]; (iii), E and F (unchanged).

**What is wrong (majors first)**

1. **D3ca-1 (major): for K′ < 0 the second dipole phase is not an approach to an infimum.** For targets 18, 12 and 29 the three-kernel block lies *below* its optimal dipole fit by |K′|/v², and the pair retracts away from it. "That block's asymptotic residual" is likewise not r′_dip. Had τ₄ not intervened, the dipole would have dissolved at τ = 145.0, 154.1 and 61.5.
2. **D3ca-2 (major): on the symmetric pair the loop does not keep a cancellation.**
   - Run to τ = 10⁹, at 10⁻⁶ the extra kernels decay by ×0.73 per decade.
   - At 10⁻⁸, and at HP's end (the same state), the two features are fitted exactly by one kernel each (error 1.4·10⁻⁵ at τ = 10⁹). The negatives belong to a wide null cluster leaving the window.
   - So on the symmetric pair the cancellation is the polish's. Round 74 reversed this claim, and the reversal is wrong.
3. **D3ca-3 (major): A2's 10A "spike end" comes from the horizon.**
   - The 10A kernel converges exactly at s ≈ 2–6·10¹⁰. The spike is J_N's σ → 0 plateau, and 3A and 5A sit on it too.
   - The round-74 test asserts the spike.
4. **Minors:**
   - the collapse's finite-p correction: the mean width's share is p²/36, not p²/4; the centre of mass moves; E's p⁹ term comes from the neighbouring kernels, and its sign depends on the target (D3ca-4);
   - the loop's κ carries the fine grid's offset, and the loop switches branches at larger ε (D3ca-5);
   - target 20's GT row is roundoff, and nothing checks the sign of κ³ (D3ca-6);
   - the dipole energy's next order is O(v⁻⁴), not O(v⁻³), and the rate claim depends on it (D3ca-7);
   - 10 behaviour-changing surviving mutants, among them round 68's R8 (D3ca-8).
5. **Cosmetic:** literature (D3ca-9) and wording (D3ca-10).

## What I checked and how

1. **Reproducibility.**
   - `--from-raw` with `MW_OUT=W/mwout`, compared with cmp; pytest on the test file.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/gt_detail.py`: `kappa_from_path` on the six paths, compared with the raw, plus its internals (closure residuals, moments, the order of wall arrivals).
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/gl_rerun.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/hll_rerun.py`: GL, HLL and A2 against the raw.
2. **Second phase.** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/phase2.py`:
   - For every stored loop state at 10⁻¹⁰ in (τ₃, τ₄): E, the closure fit's E′_dip, (E − E′_dip)v² and K₃ at the state.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/t18.py`: the loop's kernels for targets 18 and 27 around τ₃.
3. **κ.** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/kappa_corr.py`:
   - the law's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} and the corrected κ;
   - τ₃'s and τ₄'s shortfalls by mass.
4. **Collapse.**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/slowman.py` (and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/slowman2.py` for target 9):
     - E_min(Δc) at fixed v by Levenberg–Marquardt with an analytic Jacobian, over W, c_a, σ_a, σ_b and the other formed kernels' (w, c, σ), on 61 values of Δc;
     - the local slope d ln E/d ln Δc;
     - the induced metric (dq/dΔc)ᵀ Mob⁻¹ (dq/dΔc), split by coordinate;
     - run four ways: on the samples, on a wide window [−8, 8], with the pair alone, and with the pair alone on the wide window.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/gl_compare.py`: the prediction against GL's intervals, averaged as GL averages.
   - By hand: C; the mean width's share; the midpoint's motion.
5. **Moment algebra.** By hand, the 4th-order coefficient for signed weights (the sympy script failed: sympy is not installed); and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/kv.py`: v²(E − E_inf)/K − 1 against v.
6. **(ii)(b).**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/visible.py`: counts at three visibility thresholds.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/hll_long.py`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/hll_null.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/hl_long_any.py`: `hll_job` to τ = 10⁹ for the symmetric pair at 10⁻⁶ and 10⁻⁸, A₂ 0.99 at 10⁻⁸ and A₂ 0.95 at 10⁻⁶. They report the main kernels' error, and the other kernels' Σ|w|, their sum on the samples, widths and centres.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/hp_eq.py`: HP, H and HL ends.
7. **(i).** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/a2_long.py`: `a2_job` with S_A2 to 10¹³ for 3A, 10A and 20A at the three ratios; and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/a2_detail.py` on the raw.
8. **Mutants.**
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/mutants.py` holds 84 single-point mutants of round 74's code.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/run_mutants.py` runs each against the test file and against `--from-raw`, comparing the md byte for byte.
   - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/gt_mut.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref80/D3ca_work/scripts/surv_check.py`: whether the survivors change results on the data.
9. **Literature.** Four web searches (results in Sources). I fetched no page; arxiv.org and ams.org are unreachable.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) J, its single critical point, its limits (8–13; proof 94–96) | holds | re-derived, including ∂_σE = ½π^{1/2}w(w − A) at (c*, σ*); B [md] |
| (i) convergence from E(start) < ½‖f*‖² (13–15; 97–102) | holds | re-derived (compactness, sign of w, Łojasiewicz) |
| (i) "sufficient, not sharp": 2.5A, 3A, 5A converge, "dipping"; "from 10A it ends on the samples as a spike" (15–18; A2 181–183; Status 137–138) | the first half holds; the 10A end does not | 10A converges at s ≈ 2–6·10¹⁰. The "dip" is a plateau lasting to s ≈ 3·10⁵ (3A) and ≈ 10⁷ (5A) [own] (D3ca-3) |
| (i) on the samples, J_N's limits (19–20) | holds | A2's plateau is J_N's σ → 0 limit: E = 0.644, 0.906, 0.976 of ½⟨f*²⟩ for every w [own] |
| (ii) C's outcomes 71 / 3 / 6 (22–24; 184–189) | holds | [md] |
| (ii)(a) one Gaussian to O(Δc⁴), the split, validity (25–30; 102–108; 186–192, 202–203) | holds, cancelling pairs included | fourth-order coefficient −Wx(1−x)(1−x+x²)Δ⁴/36 ≠ 0 [own]; splits [md] |
| (ii)(b) the path's cancellation: 11–14 / 6–7 visible negatives, every formed visible kernel ≥ 0.049 wide (31–34; 224–230) | holds | counts unchanged for visibility thresholds 10⁻¹⁰ … 10⁻³ [own from raw] |
| (ii)(b) H2 (34–37; 230–235) and the separation window, H3 (37–39; 239–241) | holds | [md] |
| (ii)(b) asymmetric pairs: HP; the loop at 10⁻⁸ and 10⁻⁶; HLL "undecided" (39–44; 235–249) | holds | to τ = 10⁹: A₂ 0.99 at 10⁻⁸ max\|w\| 2.38 → 2.67 (increments 0.096, 0.077, 0.064, 0.053), residual 7.9·10⁻⁹; A₂ 0.95 at 10⁻⁶ keeps 3 negatives [own] |
| (ii)(b) symmetric pair: "the loop too fits with kernels of both signs at \|w\| < A … and so does the path without its polish … what the polish adds is the inflated cluster" (44–47; 247–249; Status 141–142) | does not hold as an outcome | at 10⁻⁶ the extra kernels decay; at 10⁻⁸ and in HP an exact one-sign fit plus a null cluster [own] (D3ca-2) |
| (ii)(c) closure, the block's law, K from the third order (48–55; 108–117; G 193–198) | holds | [md]; but E − E_inf = K/v² − 8.68K/v⁴, not "+ O(v⁻³)" [own] (D3ca-7) |
| (ii)(c)(1) dipole phase, τ₃ from the duals, ε^{−1/6} (55–61; 117–123; 204–207) | holds | GT bit-identical [own]; τ₃'s shortfall ×1.9–2.5 per decade |
| (ii)(c)(2) second phase: the law with K′ of either sign, τ₄ from the duals (61–68; 123–127; 207–213) | the law and τ₄ hold; "approaches its own dipole infimum" and "that block's asymptotic residual" are false for K′ < 0 | (E − E′_dip)v² = K′ < 0 along the loop's states for 18, 12, 29 [own] (D3ca-1) |
| (ii)(c)(2) target 20: K′ = 0, lands right after τ₃ (68; 210) | holds in substance; the printed K′ and τ₄ are roundoff | closure rms 4.4·10⁻¹⁷ (D3ca-6) |
| (ii)(c)(3) landing at p ≲ 0.28; κ³; "κ converges as ε → 0, at τ₃'s and τ₄'s rate" (69–75; 127–128; 210–217) | holds once the loop follows the path's sequence; the md's loop κ is biased by the grid | corrected: −0.27 … +0.21 % at 10⁻¹²; target 27 switches branch at ε ≥ 10⁻⁹ [own] (D3ca-5) |
| (ii)(c)(4) collapse law, a, the metric, C (75–85; 128–135; GL 218–222) | holds; the mean width's share is p²/36, and c̄ moves | induced metric over every coordinate [own] (D3ca-4) |
| (ii)(c)(4) "the rate is the law's over 1 + C, to E's next order" (84; 220–222) | holds per target; the next order comes from the neighbouring kernels, with either sign | pair alone: no p⁹ term; target 20 +0.97p, target 9 −0.58p [own] (D3ca-4) |
| (ii)(c)(4) a pair landed at small p collapses the same way, later (85; 222–224) | holds | target 9: 0.9989 predicted with C, 0.9986–0.9990 measured [own, md] |
| (ii)(c) at larger masses the stages overlap (85–86) | holds | [md] |
| (iii), E, F (88–92; 249–258) | holds | unchanged; [md] |
| Status (137–148) | overstates (b) and (c), and (i)'s 10A | D3ca-1, D3ca-2, D3ca-3 |
| Classical, literature (148–169) | fits | Ho–Nguyen's systems have positive weights; "to fourth order" is ambiguous (D3ca-9) |
| New (169–174) | fits | its last clause, "which the loop keeps on asymmetric pairs", is consistent with HLL to 10⁹ but remains undecided |

## Status of the round-74 issues

| D3bu | round-74 issue | status | what remains |
|---|---|---|---|
| 1 | The intermediate phase is a second dipole phase; κ from the path (major) | fixed, with a new error | **Fixed:** K′, τ₄, κ and v_land from the path (bit-identical); the slopes over 3μ₀K′Λ (0.9997–1.0012); "κ converges". <br>**Remaining:** <br>• for K′ < 0 the block is not approaching an infimum (D3ca-1); <br>• the loop's κ carries the grid offset, and branches switch (D3ca-5); <br>• target 20 is handled through roundoff (D3ca-6). |
| 2 | Dead spikes counted; the symmetric loop (major) | partly fixed | **Fixed:** the visibility filter (robust to the threshold), HP's stop reasons, HLL. <br>**Remaining:** the symmetric-loop claim adopted from D3bu-2 is wrong as an outcome (D3ca-2). |
| 3 | The collapse law's constant (minor) | fixed, with errors | **Fixed:** a defined; C; τ_L; the GL constant table; target 9 added. <br>**Remaining:** <br>• the mean width's share is p²/36, not p²/4; <br>• c̄ moves on the manifold; <br>• the remainder comes from the neighbouring kernels and has either sign (D3ca-4). |
| 4 | The landing statistic (minor) | fixed | p uses \|Δθ\|; the threshold is called derived. Only GL's p, which uses Δc alone, is left (cosmetic, D3ca-10). |
| 5 | (i)'s large coefficient (minor) | fixed, with a new error | "Sufficient, not sharp" holds. The 10A "spike end" comes from the horizon (D3ca-3). |
| 6 | Surviving mutants (minor) | partly fixed | Of round 74's code: 56 of 84 mutants are killed by the tests and 9 by `--from-raw`. 19 survive, 10 of them behaviour-changing, among them round 68's R8 (HL's 0.1 % threshold) (D3ca-8). |
| 7 | The τ₃ proof; the grid offset (minor) | mostly fixed | **Fixed:** the residual's limit and its rate in the proof; the law at v's own time in the v(τ₃⁻)/law column. <br>**Remaining:** <br>• κ still uses v at the grid time (D3ca-5); <br>• "O(v⁻³)" should read "O(v⁻⁴)" (D3ca-7). |
| 8 | Literature and New (minor) | fixed | Ho–Nguyen's systems have positive weights; one tag can be upgraded (D3ca-9). |
| 9 | Cosmetic | fixed | "idle" is dropped; "72 rows, 59 distinct"; the stale t3_pred is gone from G's raw. |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ca-1 | major | **For K′ < 0 the second dipole phase is not an approach to an infimum (D59.md:61–68, :111–112, :123–127, Status :145–148).** <br>• **The data.** I took G's stored loop states at 10⁻¹⁰ in (τ₃, τ₄), with three formed kernels, and ran `closure_fit` from each state. (E − E′_dip)v² is: <br>  – −1.4816·10⁻⁴ … −1.4813·10⁻⁴ for target 18 (K′ = −1.4822·10⁻⁴); <br>  – −7.496·10⁻⁵ for target 12 (K′ = −7.500·10⁻⁵); <br>  – −1.2151·10⁻⁴ … −1.2147·10⁻⁴ for target 29 (K′ = −1.2157·10⁻⁴) [own]. <br>  For example, target 18 at τ = 36 (v = 15.03) has E = 2.04886·10⁻⁴ < E′_dip = 2.05542·10⁻⁴. <br>  For the positive cases it is +1.1834·10⁻⁴ (target 9; K′ = 1.1796·10⁻⁴) and +1.3920·10⁻⁴ (target 27; K′ = 1.3918·10⁻⁴). <br>• **What this means.** For K′ < 0 the block lies below its optimal dipole fit. That fit is not the infimum of the block's family; the energy falls *away* from it as v retracts. <br>• **The residual.** Likewise, the s → ∞ residual of that block is not r′_dip. What moves the duals is the residual at v = ∞, which the loop holds only because v ∝ ε^{−1/6}. <br>• **Without τ₄ the dipole would dissolve.** v³ = 0 at τ₃ + K(τ₃ − τ_d)/\|K′\| = 145.0, 154.1 and 61.5 (targets 18, 12, 29). τ₄ comes first, at 97.9, 81.6 and 47.0. The pair loses 58 %, 24 % and 55 % of v³ (κ³ = 0.425, 0.757, 0.455). <br>• **What still holds.** The law itself: slopes 0.9997–1.0012 [md]. <br>• **Code.** `kappa_from_path` takes `np.cbrt` without checking the sign. Under mutant GT1 it prints κ(12) = −3.30 silently. | **(ii)(c)(2) and the proof:** <br>• For K′ > 0 (targets 9, 27): "approaches its dipole infimum from above". <br>• For K′ < 0 (18, 12, 29): "the optimal dipole fit is not an infimum; the block lies below it by \|K′\|/v² and the pair retracts. The law holds while K(τ₃ − τ_d) + K′(τ − τ₃) > 0, i.e. if τ₄ < τ₃ + K(τ₃ − τ_d)/\|K′\| (here 97.9 < 145.0, 81.6 < 154.1, 47.0 < 61.5)". <br>• Replace "asymptotic residual" with "the residual of its optimal dipole fit". <br>**Status:** restrict "the flows approach the dipole infima" to K, K′ > 0. <br>**Code:** `kappa_from_path` returns None if the bracket is ≤ 0; GT prints the dissolution time. |
| D3ca-2 | major | **On the symmetric pair the loop's both-sign kernels are transient, or a null cluster beside the exact fit (D59.md:44–47, :247–249, Status :141–142).** I ran `hll_job` to τ = 10⁹ [own]; at 10⁵–10⁷ it agrees with the raw to 2·10⁻⁹. <br>• **ε = 10⁻⁶.** <br>  – The 9 kernels other than the two nearest the features have Σ\|w\| = 0.676, 0.494, 0.362, 0.271, 0.185 at τ = 10⁵ … 10⁹ (×0.73 per decade). <br>  – The central kernel goes 0.356 → 0.109. <br>  – The main kernels' largest error in (c, σ, w) goes 0.295 → 0.091 (w 0.705 → 0.909). <br>  – Visible negatives go 8 → 6. <br>• **ε = 10⁻⁸, and HP's end** (the same state at τ = 10⁵, to 10⁻³). <br>  – The kernels at ±0.080 fit the features to 7.1·10⁻², 4.3·10⁻³, 5.0·10⁻⁴, 7.6·10⁻⁵, 1.4·10⁻⁵ (τ = 10⁵ … 10⁹). <br>  – The other 11 kernels (10 negative, Σ\|w\| 1.90 fixed) are a wide null cluster: widths 0.2–1.95 (the central one σ 1.59 → 1.95, w +0.45), centres out to \|c\| 1.69 → 2.13, beyond the window. <br>  – Their sum on the samples falls from 8.9·10⁻² to 2.3·10⁻⁵. <br>• **Conclusion.** Neither the loop nor the polish-free path fits the symmetric pair with a cancellation. On that pair the cancellation (H: a +3.24 kernel with 7 negatives) is the polish's. <br>• **Asymmetric pairs.** These do keep their clusters to τ = 10⁹: A₂ 0.99 at 10⁻⁸ has 6 negatives (Σ −2.30 → −2.11); A₂ 0.95 at 10⁻⁶ has 3 (Σ −1.02 → −0.97). | **(ii)(b), the Status, HLL:** <br>• "on the symmetric pair the cancellation is the polish's: the loop's negatives decay (10⁻⁶) or form a null cluster beside the exact one-sign fit (10⁻⁸, HP)"; <br>• run HLL's symmetric cells to τ = 10⁹ and report the extras' Σ\|w\| and the main kernels' error; <br>• in H, HP and HL, tell a null cluster (its sum on the samples below the residual) apart from a cancellation. |
| D3ca-3 | major | **A2's 10A spike is a plateau, not an end (D59.md:15–18, :181–183, Status :137–138; `test_a2_one_kernel_from_large_coefficients`).** I ran `a2_job` with S_A2 to 10¹³ [own]. <br>• **10A.** The kernel sits on the spike (σ/σ* = 0.0772, 0.0194, 0.0049) until s ≈ 3·10¹⁰ (ratio 0.25) or 1.8·10¹⁰ (ratios 1, 4). Meanwhile w/A falls: 9.972, 9.910, 9.690, 8.449 at s = 10⁹, 3·10⁹, 10¹⁰, 3·10¹⁰. It then converges exactly (w/A = σ/σ* = 1) by s = 5.6·10¹⁰ or 3.2·10¹⁰. <br>• **5A.** The "dip to 0.0058σ*" is the same spike, held from s ≤ 1 to 10⁷ (raw: w/A 4.93, 4.89, 4.89 at 10⁷); exact by 10⁸. <br>• **3A.** At ratios 1 and 4 the kernel sits on a spike (σ/σ* 0.028, 0.0070) until 3·10⁵. <br>• **The plateau.** Its E is the same for every w: 0.644, 0.906 and 0.976 of ½⟨f*²⟩. It is J_N's σ → 0 limit (line 19), a plateau rather than a minimum. Its duration grows steeply with w: 20A is still on it at 10¹³. <br>• **The test** asserts the 10A spike at s = 10⁹. | **(i) and A2:** "from 3A on, the width collapses at once onto the nearest samples (J_N's σ → 0 limit, E 0.64–0.98 of ½⟨f*²⟩). The kernel stays there for s ≈ 3·10⁵ (3A), 10⁷ (5A), 2–6·10¹⁰ (10A) and > 10¹³ (20A), then converges exactly". <br>**Code:** S_A2 to 10¹³, with the escape time as a column. <br>**Test:** assert 10A's escape, not the spike. |
| D3ca-4 | minor | **The collapse's finite-p correction (D59.md:82–85, :129–135, :220–224).** I computed target 20's slow manifold at 10⁻¹⁰: E minimised at fixed (Δc, v) over all other coordinates of the four formed kernels. The loop's τ = 400 state lies on it to 7.6·10⁻⁹ [own]. <br>• **(a) C holds.** The coefficients' share is 0.12338 against 0.12299 at p = 0.148; the other kernels' shares are below 10⁻⁵. <br>• **(b) The mean width's share is p²/36, not p²/4:** 6.3·10⁻⁴ at p = 0.148, against p²/36 = 6.1·10⁻⁴ and p²/4 = 5.5·10⁻³. The reason: (a)'s mean-width shift is the shift of the coefficient-weighted mean Σw_kσ_k/W = σ̄ + (v/W)(σ_a − σ_b). With the split set by the third moment, the arithmetic mean moves only a third as far, (1/6)v²Δc²/(W²σ). <br>• **(c) c̄ moves on the manifold.** It shifts by δc̄ = m₃/(6Wσ²), which is −p²/4 of the separation's leading motion. W's motion contributes +p²/4 in the midpoint, so the centres' share stays (w_a² + w_b²)/W² to O(p⁴) (0.99987 at p = 0.148). The proof's "restricted to dc̄ = 0" gets the right metric through this cancellation without stating it. <br>• **(d) E's next order.** <br>  – The local slope is d ln E/d ln Δc − 8 = 0.955p, so the force grows by 1 + 0.119p; D59 has 1 + 0.131p. With (a)–(c), the prediction matches GL's ratio within 8·10⁻⁴ on 27 of 28 intervals; the first lies at my grid's edge. <br>  – But for the pair alone, on the same samples or on a wide window, the slope minus 8 is 0.03p + 0.33p² and −0.0005p + 0.36p²: there is no p⁹ term, because the order-k terms have parity (−1)^k about c̄. <br>  – So the p⁹ term comes from the other kernels absorbing part of the fourth-order residual, and its sign depends on the target. Target 9 has −0.58p + 0.23p², giving 0.9989 including C against the measured 0.9986–0.9990 [md]. | **Proof and GL:** <br>• the mean width's share is p²/36; <br>• state that c̄ moves and that W cancels it in the midpoint; <br>• "the remainder is E's next order — O(p) only through the neighbouring kernels (O(p²) for an isolated pair), with a coefficient set by the target (+0.97 for 20, −0.58 for 9)"; <br>• print the local slope instead of the linear fit's k = 1.048. |
| D3ca-5 | minor | **The loop's κ carries the fine grid's offset; κ converges only once the loop follows the path's sequence (D59.md:74, :211–212, :215; `g_fine_lines`).** <br>• **The offset.** κ_loop = v_land/v(t_v3), where t_v3 is the grid time before τ₃. At 10⁻¹² the law's offset is +0.18, +0.10, +0.13, +0.58, +0.19 and +0.32 % (targets 9, 18, 12, 20, 27, 29), and all six loop κ lie above the path's. <br>• **Corrected,** the loop's κ differs from the path's by −0.05, −0.03, +0.21, −0.03, −0.06, −0.27 % at 10⁻¹² and by −0.11, −0.05, +0.40, −0.06, −0.12, −0.62 % at 10⁻¹¹. <br>• **Target 20.** Its 1.005–1.006 at 10⁻¹⁰ … 10⁻¹² is just its +0.58 % offset; corrected, κ = 0.9988, 0.9994, 0.9997. <br>• **Branch switches.** Target 27 at ε ≥ 10⁻⁹ and target 18 at 10⁻⁷ follow another sequence. At 10⁻⁹, target 27's τ₃ = 44.8 and the pair relocates to the feature at −0.684 (at τ = 48: w −8.34/+8.81 at c −0.661/−0.662). The intermediate slope is 4.40 times the path's law, and κ = 1.758, 1.417, 1.052. <br>• **Rates.** τ₄'s shortfalls fall ×2.4–4.2 per decade, faster than ε^{1/3}. | **`g_fine_lines`:** take v(τ₃⁻) from the law at the crossing (or at the path's τ₃), or quote v_land only against the law. <br>**(3):** "converges once the loop follows the path's sequence (here ε ≤ 10⁻¹⁰)"; O(ε^{1/3}) is a bound. |
| D3ca-6 | minor | **Target 20's GT row is roundoff (D59.md:68, :210; `kappa_from_path`).** <br>• The four-kernel closure fits exactly: rms 4.4·10⁻¹⁷, \|p′\| 2.2·10⁻⁹. So K′ = 2.0·10⁻⁴³ and τ₄ = 1.41·10¹⁷ are noise. <br>• Seeding the entrant at half width (KP7) gives 1.2·10⁻²⁷ and 3.5·10¹³ instead, with κ = 1 either way. κ = 1 only because K′(τ₄ − τ₃) ∝ \|p′\|³. <br>• The md prints "1.414e+17 (—)". | **`kappa_from_path`:** if the closure's residual is below a tolerance, set K′ = 0, τ₄ = None and κ = 1 ("lands at τ₃"). <br>**GT:** print "— (exact)". |
| D3ca-7 | minor | **E − E_inf = K/v² + O(v⁻⁴), not O(v⁻³) (D59.md:111; the docstring of `test_d3phi_and_the_third_order_k`).** <br>• **Why.** Expanded about the optimal dipole fit, the pair's sum has only even powers of 1/v: W multiplies the even orders of Δθ = p/v, and w_a − w_b = 2v the odd ones. <br>• **Check.** In the test's own setup, v²(E − E_inf)/K − 1 = −8.68/v² for v = 10 … 320 (−8.68·10⁻² … −8.57·10⁻⁵) [own]. The test says the error falls "like 1/v". <br>• **Why it matters.** An O(v⁻³) term would give v³ = 3μ₀Ks(1 + O(s^{−1/3})), a law error of O(ε^{1/6}). That would dominate the O(ε^{1/3}) claimed for τ₃, τ₄ and κ. | **Proof (:111):** "O(v⁻⁴)", with the parity reason. <br>**Test:** fix the docstring. |
| D3ca-8 | minor | **Tests' strength.** 84 single-point mutants of round 74's code: 56 are killed by the tests and 9 more only by `--from-raw` (GF4, GF5, GF7, GC3, GC5, GL1, GL3, GL5, AL1). 19 survive. <br>• **Equivalent on the data (9):** PP3, FW3, KP5, KP7, FPa, GJa (no kernel between 0.1 % and 1 % in GL's runs; idle \|w\| ≤ 3.3·10⁻¹²), GJc, LR2, LR3. <br>• **Behaviour-changing (10):** <br>  – GT1: `gt_job` is untested; with w scaled by 0.9, κ(12) → −3.30 and κ(27) → 23.5; <br>  – FPf, FPg, FPh in `fine_phases`: a landing before τ₃; the intermediate state taken at τ₃; τ₄ at an exit; <br>  – GJd: σ̄ from one kernel in `gl_job`; <br>  – HS3: H's "formed invisible" would count the 3–4 unformed invisible active kernels at H's ends; <br>  – HLa: HL's formed threshold at 0.1 %, i.e. round 68's R8, still alive (the round-74 test covers `hll_job`, not `hl_job`); <br>  – HLLb: ε × 10 in `hll_job` (the test patches `loop_nd` and never checks ε); <br>  – A2b: σ/σ₀ in place of σ/σ* (the test uses ratio 1 only); <br>  – A2d: A2 at rtol 10⁻⁶. <br>• The A2 test also asserts the artefact of D3ca-3. | **Tests to add:** <br>• `gt_job` against `kappa_from_path` on a stored block; <br>• `fine_phases` with a landed state before τ₃, an exit before τ₄, and an intermediate state that is not the midpoint; <br>• `gl_job` with unequal widths; <br>• `h_stats` with an unformed invisible kernel; <br>• `hl_job` with a kernel between 0.1 % and 1 %; <br>• `hll_job` checking that ε reaches `loop_nd`; <br>• `a2_job` at a ratio ≠ 1 and at a tight tolerance. |
| D3ca-9 | cosmetic | **Literature (D59.md:159–163).** <br>• Ho–Nguyen's r̄(2) = 4 for one extra component, rate n^{−1/8}, is now [snippet] (secondary), not [memory]. <br>• Their systems are for mixtures, with positive weights c_j². That r̄ = 4 also holds for the cancelling pair of the collapse is (a)'s own algebra: the coefficient −Wx(1−x)(1−x+x²)Δ⁴/36 [own]. <br>• "The pair one Gaussian to fourth order" (:161) can be read as an O(Δc⁵) error. | Tag [snippet]; add "(positive weights; for signed pairs, (a))"; write "to third order, the fourth-order term non-vanishing". |
| D3ca-10 | cosmetic | <br>• **Rounding.** "0.752–1.115" (:74): the path's 1.1145 rounds to 1.114. <br>• **A negative width.** Target 18's phase-2 closure puts the entrant at σ = −0.0970 (the loop's kernel has +0.0971). This is harmless, since φ is even in σ. <br>• **Two definitions of p.** GL's p uses \|Δc\| (`gl_constant` and the GL header), while (3) defines p with \|Δθ\|. Here they differ by 0.03 %. <br>• **A missing target.** "Past τ₄ 22800 of 22978" (:217) leaves out target 20, which has no τ₄. | Edit the text; take \|σ\| in `closure_fit`'s output; say which p each statistic uses. |

## Sources

Four web searches. arxiv.org and ams.org are unreachable, and I fetched no page.

- **Ho, N. and Nguyen, X.**, "Convergence rates of parameter estimation for some weakly identifiable finite mixtures", Ann. Statist. 44(6) (2016) 2726–2755, [Project Euclid PDF](https://projecteuclid.org/journals/annals-of-statistics/volume-44/issue-6/Convergence-rates-of-parameter-estimation-for-some-weakly-identifiable-finite/10.1214/16-AOS1444.pdf). Venue and pages are [snippet], as is the technical report's abstract ([arXiv:1501.02497](https://arxiv.org/pdf/1501.02497): rates set by the order of a solvable polynomial system).
- **"r̄(2) = 4, r̄(3) = 6, r̄(4) ≥ 7"**, with n^{−1/8} for one extra component and n^{−1/12} for two, attributed to Ho–Nguyen 2016a. This is a [snippet] from a secondary source that the search did not identify precisely; the candidates are [arXiv:2002.05818](https://arxiv.org/pdf/2002.05818) and [arXiv:2601.11465](https://arxiv.org/pdf/2601.11465).
- **Xu, W. and Du, S. S.**, "Over-Parameterization Exponentially Slows Down Gradient Descent for Learning a Single Neuron", COLT 2023, [PMLR 195:1155–1198](https://proceedings.mlr.press/v195/xu23a.html). O(T⁻³) and a matching Ω(T⁻³) for the loss, via a potential on the students' pairwise distances [snippet].
- **[memory], not re-verified:** de Boor 1969; Kammler 1971–76; Braess 1986; Petersen–Raslan–Voigtlaender 2021; Gallon–Jentzen–Lindner 2022; Dwivedi et al., Ann. Statist. 48(6) 2020; Chizat 2022; Pearson 1894; Runnalls 2007; Schiebinger–Robeva–Recht 2018; Duval–Peyré 2017. Each fits the way D59.md uses it.

## Lesson

Two of this round's three majors are claims about where a flow ends that the flow, run longer, does not bear out. A2's 10A "spike end" at s = 10⁹ and the symmetric loop's "both-sign fit" at τ = 10⁷ were states at the horizon. Two to four more decades show the kernel escaping to the exact fit, and the cancellation decaying or turning into a null cluster that leaves the window. An end needs a trend that has stopped, not one that has slowed. The third major is about naming: the same v³ law covers both signs of K′, but only K′ > 0 is an approach to an infimum. For K′ < 0 the block sits below its dipole fit and retracts, and the four stages survive only because τ₄ arrives before the dipole dissolves. Finally, the isolated pair is a cheap control for an attribution. E's p⁹ term, which explains GL's remainder, vanishes for a pair alone by parity. It belongs to the neighbouring kernels, and target 9's opposite sign said so.
