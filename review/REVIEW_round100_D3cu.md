# Referee report D3cu, round 100: Proposition D‴.59 as revised after round 97

## Scope

- **Export.** I exported commit d0ec681 (branch claude/vigilant-clarke-hf939b) with `git archive d0ec681 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref100.
- **Base for the diff.** The base is 0ea58aa, the commit the round-97 referee reviewed. The diffs run 0ea58aa → d0ec681.
- **The raw has not changed since 66cc9e1 (round 87).** This round changes only code that runs when the report is built (`a2_invariance`, `fit_split`/`split_summary`/`split_cells`, the A2E_ROUND filter) and text.
- **Read in full:**
  - review/drafts/D59.md (354 lines).
  - Every line of exp_moving_width.py that changed since 0ea58aa: A2E_ROUND, `a2_invariance`, the symmetric-line and exact-start lines of `a2_lines`, the roundoff filter and its count, SPLIT_CS, NULL_TOL, `fit_split` (ecs/er, main_cs, near_cs, rest_win), `split_summary`, `split_cells`, and the HLL null line in `report`.
  - The code these call: `a2e_job`, `outcome`, `match`, `end_split`, `critical_rhs_nd`/`critical_jac_nd`, `BumpsCW`.
  - results/chaos/moving_width.md (619 lines) and tests/test_moving_width.py (1125 lines, 49 tests).
  - review/REVIEW_round97_D3cr.md, the Round 97 entry of review/REVIEW_LOG.md, tools/mutants/round97_D59_own.py, the D59 row of review/drafts/README.md, and the diffs of all of these.
- **Settings.**
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, in the repository's venv (numpy 2.5.4, scipy 1.18.1).
  - MW_OUT pointed to W/mwout (holding a copy of the raw) or to each mutant slot's `out` directory.
  - pytest used `--basetemp` inside W. My copy of the mutant tool adds `--basetemp` per slot.
  - At most two heavy processes of mine ran at a time, and every run had a timeout.
- **I modified nothing in the repository.** I used only read-only git (`archive`, `log`, `diff`, `status`). At the end, D59.md, the script, the md, the raw, the tests and the round-97 spec are byte-identical to the export, and nothing of mine is running.
  - Meanwhile another session committed 0cb96bd and 49d78c8 (D‴.58, D‴.61) and is editing D‴.60's files. None of these is a D‴.59 file, and I touched none of them.
- **One wrong step.** A wait loop of mine (`until ! pgrep -f "spec_r100.py"`) matched its own command line and would never have ended. I confirmed through /proc/PID/cmdline that it was mine, then stopped it by its task ID (TaskStop). I did not use pkill, and no result depended on it.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my own run or computation; [snippet] a web-search result; [memory] not verified this round.

## Summary

**Design.** I worked in this order:
1. Reproduced the md from the raw (`--from-raw`) and ran the tests.
2. Reran a sample of jobs (all 15 A2, all 3 A2E, one H, one HL) against the raw.
3. Tested the invariance claim three ways:
   - the exact-arithmetic drift of ċ on c = c* from the six unpaired samples, computed in log space;
   - the float drift on the line;
   - LSODA's steps of the exact-start 3A run all the way to s = 10 (the code records steps only to s = 10⁻²).
4. Replayed A2E's BDF steps with the float drift included.
5. Derived the symmetric-sample lemma and checked its decomposition against the Jacobian.
6. Rebuilt all 94 splits by brute force:
   - every cap from 0.5 to 4 in steps of 0.01;
   - near_cs against every kernel;
   - the window mass against the continuum mass;
   - the null test against `outcome`, with planted ends.
7. Recomputed κ′ − κ and its ratio to ε^{1/3} from the raw's fine phases.
8. Checked SQ's event lists.
9. Ran the repository's 26 round-97 mutants and 40 of my own, plus job-level reruns of A2E inside mutated copies.
10. Ran four web searches.

**Verdicts.**
- **Reproducibility is exact.** The md matches byte for byte (0.8 s), 49 of 49 tests pass, and all 20 rerun jobs are bit-identical to the raw.
- **The numbers quoted in the draft are the md's**, apart from rounding (checked line by line).
- **The κ′ − κ rates, SQ's counts and type, the window masses, the gap 1.15/2.81 and the null-cluster row all reproduce** with my own code.
- **Round 97's issues are fixed in wording**, but the premise of D3cr-1 was wrong and the new text inherits it:
  - The line c = c* is invariant only for samples exactly symmetric about c*.
  - For the problem as posed, both the exact flow and the float flow leave the saddle in finite s.
  - The runs that "stay" are held by BDF damping, at the implicit step's fixed point of the samples' drift.
- **Two order-recording survivors of round 97 (AE2, AE3) still survive** and still rewrite a published line on a rerun.

**Findings.**
1. **D3cu-1 (minor; the principal finding): "invariant for the problem as posed" is wrong, and so are its consequences.**
   - The six unpaired samples carry f*/A ≈ 5·10⁻⁶⁰³ at σ*/σ₀ = 0.25 (it is 0 only in double) and 2.3·10⁻³⁸ at 1.
   - The saddle amplifies this. In exact arithmetic, the flow from c* reaches an offset of 10⁻²σ* at s ≈ 4.2 (3A) and 0.26 (5A) at 0.25, and at 0.55 and 0.021 at σ*/σ₀ = 1.
   - In double the line drifts at ċ ≈ −1.0·10⁻¹¹σ* per unit s.
   - LSODA's held runs sit at BDF's fixed point −d/λ.
   - A2E's step left out "into roundoff" is that drift.
2. **D3cu-2 (minor): the invariance and the saddle can be proved, not just evaluated.** For samples symmetric about c*, a lemma gives both: the line is the fixed set of a reflection, the transverse rate is −Gm∂²_cE_N, and the saddle comes from the kernel overshooting its sample.
3. **D3cu-3 (minor): NULL_TOL is not `outcome`'s test.** Planted ends make the two disagree in both directions. `split_summary`'s null cluster can be a fit with residual 1.1·10⁻⁷.
4. **D3cu-4 (minor): the tests.**
   - AE2 and AE3 survive and on a rerun rewrite the A2E line. The log's "kills AE2, AE3" is wrong.
   - The ratio-4 checks of `a2_invariance` sit inside a comment (test line 1099).
   - The null test's main-kernel tolerance is untested (R19).
5. **D3cu-5 (cosmetic): the gap's justification.** The gap is right on the data: by brute force no assignment changes for any cap in (1.1498, 2.8144). But near_cs alone does not establish it, 29 of the 30 unassigned features (not 30) have an unused kernel, and the 1.15 maximum is a budget stop's shared kernel.
6. **D3cu-6 (cosmetic):** "f* is 0" at 0.25; "the residual's threshold"; and the README's round-83 phrase "settling on … 29", still unchanged.

## What I checked and how

All paths are under W.

1. **Reproducibility.**
   - `--from-raw` into W/mwout, then `cmp` with the committed md: identical.
   - pytest: 49 passed in 14.4 s (logs/pytest_full.log).
   - scripts/rerun.py compares every field except timings, bit for bit (logs/rerun_*.log). Identical: all 15 `a2_job`, all 3 `a2e_job`, `h_job` (σ* 0.04, equal) and `hl_job` (σ* 0.06, σ₂ 0.065, 69 s).
2. **Invariance in exact arithmetic** (scripts/exact_seed.py, logs/exact_seed.log).
   - I integrated the symmetric trajectory with c pinned at c* (LSODA on the repository's samples).
   - Along it I computed:
     - the exact-arithmetic drift from the unpaired samples x = −1.6 … −1.55 (the paired ones cancel exactly), in log space;
     - the transverse rate λ(s) = ∂ċ/∂c from `critical_jac_nd`;
     - the linearised offset δ′ = λδ + d from δ(0) = 0, as a log-space quadrature with |d| as the envelope.
   - Results [own]:

     | σ*/σ₀, w | log₁₀\|ċ\| at s = 0 | λ on the saddle (max) | δ reaches 10⁻²σ* at s ≈ |
     |---|---|---|---|
     | 0.25, 3A | −1202.2 | 665 | 4.23 |
     | 0.25, 5A | −1201.7 | 1.089·10⁴ | 0.256 |
     | 1, 3A | −74.0 | 333 | 0.554 |
     | 1, 5A | −73.5 | 1.008·10⁴ | 0.021 |

   - Each decade of cancellation in the seed would shift these times only by ln 10/λ: 0.0035 (3A) and 0.0002 (5A) at 0.25.
3. **LSODA's held 3A run** (scripts/a2e_3a_steps.py, logs/a2e_3a_steps.log), stepped as `a2e_job` does but to s = 10.
   - Adams up to s ≈ 0.01.
   - BDF from there on: orders 1–2 at hλ 0.70–15.6 on (0.01, 0.1], order 2 at hλ 15.6–111 on (0.1, 1], order 3 at hλ 296–1.25·10³ on (1, 10], with λ ≈ 629–665.
   - The offset rises to 9.4·10⁻¹³σ*, is then damped, and ends at +1.04·10⁻¹⁵σ*.
4. **The float drift and BDF's fixed point** (scripts/float_drift.py, scripts/bdf_drift.py, logs/bdf_drift.log).
   - On the line at A2E's held states [own]:
     - 5A: ċ(c*) = −2.95·10⁻¹³ (−1.02·10⁻¹¹σ* per unit s), λ = 1.05·10⁴, so −d/λ = +9.37·10⁻¹⁶σ*. The held offset is +9.25·10⁻¹⁶σ*.
     - 3A: ċ(c*) = −1.90·10⁻¹⁴, λ = 629, so −d/λ = +1.01·10⁻¹⁵σ*. The held offset is +1.04·10⁻¹⁵σ*.
   - Backward Euler with the drift, δ′ = (δ + hd)/(1 − hλ), reproduces all six order-1 steps of the raw:

     | hλ | raw multiplier | replayed with the drift | backward Euler alone |
     |---|---|---|---|
     | 3.986 | ×−0.3348 | ×−0.3348 | ×−0.3349 |
     | 3.986 | ×−0.3350 | ×−0.3351 | ×−0.3349 |
     | 7.556 | ×−0.1521 | ×−0.1521 | ×−0.1525 |
     | 7.556 | ×−0.1554 | ×−0.1555 | ×−0.1525 |
     | 49.53 ("into roundoff") | ×−0.0040 | ×−0.0040 | ×−0.0206 |
     | 49.53 | ×−4.00 | ×−4.17 | ×−0.0206 |

5. **The lemma** (scripts/saddle_lemma.py, logs/saddle_lemma.log). ∂ċ/∂c = Gm·w·⟨r∂²_cφ − w(∂_cφ)²⟩, split sample by sample, on the repository's samples and on exactly symmetric ones (identical to the digits shown) [own]:

   | state | centre sample | k = ±1 | k = ±2 | total | Jacobian |
   |---|---|---|---|---|---|
   | 5A, σ/σ* 0.184, h/σ 1.81 | +15704 | −4919 | +106 | 10890 | 10890 |
   | 3A, σ/σ* 0.2345, h/σ 1.42 | +2900 | −2489 | +248 | 663.5 | 663.46 |

   The centre sample's term equals the leading order (Gm w/(Nσ²))(w − f*(c*)).
6. **The splits** (scripts/split_check.py, scripts/split_extra.py, logs/split_check.log). All 94 rows of H, H3, HP, HL and HLL, rebuilt from the raw with my own assignment [own]:
   - `fit_split` reproduces max main_cs 1.1498 (H3 ±0.05), w/A − 1 from −0.5434 to +2.2390, 30 unassigned features, and min near_cs 2.8144.
   - No assignment changes for any cap in (1.14, 2.82) (steps of 0.01). The first changes above are at 2.82 and 2.87.
   - No excluded kernel would beat a main kernel's combined error.
   - The smallest excluded same-sign kernel, unused, sits 1.805 off an assigned feature (HL σ* 0.06, equal, 10⁻⁶, τ 10⁵).
   - H3 ±0.05's second feature has no unused kernel of its sign. Its one positive kernel, at (0, 0.129) with w 1.42, is 1.15 off both features and is the first feature's main kernel.
   - Without that row the main kernels are at most 0.818 off.
   - The window mass and the continuum mass agree in sign on 94 of 94. No inflated main kernel sits beside a positive window mass, and no deficient one beside a negative one.
   - The only null-cluster row is H σ* 0.04, equal, under both `outcome` and `split_summary`.
7. **Planted ends** (scripts/null_planted.py, logs/null_planted.log). See D3cu-3.
8. **κ′ − κ** (scripts/kappa_rate.py, logs/kappa_rate.log). κ′ = v_land/(v₃·((τ₃^path − τ_d)/(t_v3 − τ_d))^{1/3}) from G's fine phases and GT [own]. Every printed value, per-decade factor and ratio over ε^{1/3} reproduces:
   - 9: ×2.22–2.63, ratio −8.6 → −5.8;
   - 20: ×1.86–2.02, ratio −1.8 → −2.7;
   - 27: ×2.21, ratio −5.9 → −5.6;
   - 29: ×1.18–2.30, ratio −12.6 → −25.9 → −21.1;
   - 12: ×1.32–1.85, ratio +3.7 → +19.6;
   - 18: +0.0035 → −0.0002.
9. **SQ** (scripts/sq_types.py, scripts/td2.py), from the raw's event lists:
   - Target 20, kernel 16 forms at τ = 24.2 (10⁻⁵) and 29.5 (10⁻⁶), before τ_d = 32.224. Kernel 4 does not form there; it forms at 36.0 at 10⁻⁷.
   - Target 18, kernel 9 forms at 23.5 and 25.0, before τ_d = 25.003.
   - Target 27, kernel 9 forms at 25.3 and 25.8, before 31.54.
   - Target 12, kernel 10 forms at 15.2 and 23.2, before 24.35.
   - Five of the kernel-9 cells land, four of them with a κ.
10. **Numbers.** scripts/numcheck_all.py, run over the whole draft against the md. The only misses are rounding (−2.46 for −2.460, 8.4 for 8.41, and so on), counts (71 of 80), percentages and years.
11. **Mutants.**
    - Tool: scripts/mutants_tool.py, which is tools/mutants.py with REPO set to the export, PY set to the venv, and `--basetemp` added.
    - Slot a: the repository's round97_D59_own.py.
    - Slots b, c, d: my scripts/spec_r100.py (34), spec_r100b.py (4) and spec_r100c.py (2).
    - Results in mut/res_{a,b,c,d}.jsonl and logs/mut_*.log.
    - scripts/jobmut.py reruns `part_A2E` inside a mutated copy and diffs the rendered md. Its no-op control changes 0 lines (logs/jobmut.log).
12. **Literature.** Four web searches. I fetched no page.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) the landscape, convergence from E(start) < ½‖f*‖², the proof (8–15, 134–142) | **holds**; proved in the continuum; unchanged | re-read |
| (i) "started just off the line c = c* (at (c* + 10⁻⁹σ*, σ*))" (16) | **holds** | [md], [raw], own reruns |
| (i) "at σ*/σ₀ = 0.25 and 1 the line c = c* is invariant … f* is 0 and 2.3·10⁻³⁸ … from c* exactly the flow stays on that saddle, as LSODA's run of 3A at 0.25 does, and a solver leaves it only by amplifying roundoff" (20–23) | **does not hold** (D3cu-1) | [own] exact drift, escape times, BDF steps, fixed point |
| A2: durations, escapes, smallest σ, E 2.23–2.76·10⁻³², 20A ends (236–240) | **holds** | [md]; 15 of 15 reruns bit-identical [own] |
| A2E: symmetric state, ∂ċ/∂c, Adams growth, plateau from 0.0056, BDF multipliers ×−2.54 … −1.67 and ×−0.335 … −0.152 (241–247) | **numbers hold**; "one further step, into roundoff" is the line's drift, not roundoff of the offset (D3cu-1) | [raw]; [own] replay |
| A2E: "the line c = c* is invariant there … the runs that stay … follow the exact flow … those that leave … amplify roundoff" (248–252) | **numbers hold**; the conclusion does not (D3cu-1) | `a2_invariance` reproduces; [own] |
| (ii)(b): cap 1.5 on (c, σ), gap 1.15 against 2.81 (45–49) | **holds on the data**; the justification is incomplete (D3cu-5) | [own] brute force |
| (ii)(b): window mass as the side (49–50, 59, 65, 69–71) | **holds**; the sign agrees with the continuum mass on 94 of 94; the classes are consistent | [own] |
| (ii)(b): null cluster by the others' rms below 10⁻⁸ of max f*, "the residual's threshold" (52–54) | **H's row holds**; the test is not `outcome`'s (D3cu-3) | [own] planted ends |
| (3): κ′ − κ per decade and over ε^{1/3}; the target-18 sentence (107–113) | **holds** | [own] = [md] |
| SQ: 34 of 48; five landings, four with a κ; target 20's type (289–298) | **holds** | [raw] event lists |
| Status (186–188) | **holds**, except "from just off the invariant line" (D3cu-1) | — |
| Mathematics | (i) proved in the continuum; (ii)(a) derived; (ii)(b) numerical; (c)'s law derived on the slow manifold, the approach numerical; (iii) numerical. The invariance and the saddle's sign condition are provable for symmetric samples (D3cu-2); the saddle's size at A2E's states is numerical; the exact-start escape is my linearised estimate | — |
| Literature (201–223) | unchanged; the lemma needs a reference for the invariance of fixed-point sets of equivariant flows | Sources |

## Status of the round-97 issues

| D3cr | status | what remains |
|---|---|---|
| 1 (the exact start on an invariant line) | **applied as proposed, but the proposal was wrong** | The line is not invariant for the problem as posed, and the held runs are BDF's. Superseded by D3cu-1 and D3cu-2. |
| 2 (mass beyond the window) | **fixed** | rest_win is printed and quoted, and the in-window masses reproduce. |
| 3 (the null test) | **partly fixed** | The tolerance is now relative to max f*, and H's row is reconciled. The test is not `outcome`'s (D3cu-3). |
| 4 (SPLIT_MAX's justification) | **fixed** | The (c, σ) cap and the printed gap are confirmed as the robustness interval. The justification is incomplete (D3cu-5). |
| 5 (order-recording survivors) | **partly fixed** | FM2 is now killed [own], and the raw's orders are asserted. AE2 and AE3 still survive and rewrite the A2E line on a rerun; the log's "kills AE2, AE3" is wrong (D3cu-4). |
| 6 (the target-18 slip) | **fixed** | — |
| 7 (stale end energies) | **fixed** (2.23–2.76·10⁻³²) | — |
| 8 (SQ by type) | **fixed**, verified from the event lists | — |
| 9 (wording) | **mostly fixed** | The −0.00401 is gone, but its filter's label "into roundoff" is D3cu-1. 29's peak and the GL comment are done. The README's round-83 phrase "settling on targets 9, 27 and 29" is unchanged (D3cu-6). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cu-1 | minor | **The line c = c* is not invariant for the problem as posed; the held runs are held by BDF.** (D59.md:20–23, :187, :245–252; md:45 and :55, :58; code `a2_lines` 424–434, 467–476, 490–494; README row.) <br>• **Exact arithmetic.** The six samples x ≤ −1.55 have no mirror image. Their f*/A is e^{−1386.9} ≈ 5·10⁻⁶⁰³ at 0.25 (0 only in double) and 2.3·10⁻³⁸ at 1. They give the line a drift of 10^{−1202} (0.25) and 10^{−74} (1) at the start. The saddle (λ up to 665 and 1.09·10⁴ at 0.25; 333 and 1.01·10⁴ at 1) amplifies it. Linearised, the flow from c* reaches 10⁻²σ* at s ≈ 4.2 (3A) and 0.26 (5A) at 0.25, and 0.55 and 0.021 at 1 [own]. That is inside A2E's own span-10 window, so "the exact flow from c* never leaves" is false. <br>• **Double precision.** linspace's asymmetry gives ċ(c*) = −2.95·10⁻¹³ (5A) and −1.90·10⁻¹⁴ (3A) on the line [own]. "To the samples' roundoff" therefore names a drift, not an invariance. <br>• **The held runs are BDF's.** The 3A exact-start run turns to BDF at s ≈ 0.01, with hλ up to 1.25·10³, beyond 2, 4, 6.67, and its offset of 9.4·10⁻¹³σ* is damped. Both held offsets are BDF's fixed point −d/λ: +1.04·10⁻¹⁵ against +1.01·10⁻¹⁵σ* (3A), and +9.25·10⁻¹⁶ against +9.37·10⁻¹⁶σ* (5A) [own]. <br>• **The step "into roundoff" is the drift.** Backward Euler with d reproduces it: ×−0.0040 raw and replayed, against −0.0206 for backward Euler alone. The drift also reproduces the other order-1 steps to 3–4 digits [own]. <br>• **Consequence.** Only samples exactly symmetric about c* keep the line. The exact start of the problem as posed is not a special case: it leaves too, only later. A2's offset of 10⁻⁹σ* exceeds both seeds, so A2's route is the generic flow, independent of the solver. No printed number changes. | **(i):** "on samples symmetric about c* the line c = c* is invariant (Lemma D‴.59a); the problem as posed is symmetric only up to six unpaired samples (f*/A ≈ 5·10⁻⁶⁰³ and 2.3·10⁻³⁸) and, in double, linspace's roundoff (ċ ≈ −10⁻¹¹σ* per unit s on the line), which the saddle amplifies within s ≈ 0.02–4.2: from c* exactly the flow leaves too; a solver stays only by damping the unstable direction (BDF at hλ ≫ 2–6.7, at its fixed point −d/λ); A2 starts off the line, the generic flow". <br>**A2E numerics, md:45 and :58, README:** say the same. <br>**Code:** record the 3A run's steps to s = 10 (S_A2E_STEPS for A2E_MORE) and print its BDF hλ and −d/λ. Replace A2E_ROUND's "into roundoff" with "set by the line's drift d" and print the backward-Euler replay with d. <br>**Status:** "from just off the line". |
| D3cu-2 | minor | **The symmetric state's invariance and saddle are stated from evaluations, but are provable.** <br>• If X = 2c* − X and f* is symmetric about c* on X, then E_N(α, 2c* − c, σ) = E_N(α, c, σ). The block's flow (α̇ = M(α)⟨rφ⟩, ċ = −Gm∂_cE_N, σ̇ = −Gm∂_σE_N) commutes with the reflection, so {c = c*} is invariant: it is the fixed set of an equivariant flow. <br>• On the line ∂ċ/∂α = ∂ċ/∂σ = ∂α̇/∂c = ∂σ̇/∂c = 0 (even × odd sums), so c decouples. Its rate is ∂ċ/∂c = −Gm∂²_cE_N = Gm·w·⟨r∂²_cφ − w(∂_cφ)²⟩. <br>• For c* a sample, the centre sample contributes (Gm w/(Nσ²))(w − f*(c*)): the saddle is the kernel overshooting the sample it sits on. The neighbours restore, at O((h/σ)²e^{−h²/(2σ²)}). <br>• Check: at 5A the centre gives +15704 and k = ±1 give −4919, total 10890 (the Jacobian). At 3A they give +2900 and −2489, total 663.5 [own]. | Add **Lemma D‴.59a** to (i)'s proof, with this decomposition as its numerical check. Print the centre and neighbour terms in A2E. Cite the fixed-point-subspace invariance (Sources). Then state the problem as posed as a perturbation of the lemma (D3cu-1). |
| D3cu-3 | minor | **NULL_TOL is not `outcome`'s test** (code comment at 1249: "(`outcome`'s residual)"; D59.md:52–54, "the residual's threshold"). <br>• `outcome` calls an end "exact, null cluster" when `match`'s nearest kernels are within 10⁻⁶ (absolute) and the residual is below 10⁻⁸. `split_summary` asks that `fit_split`'s main kernels be within 10⁻⁶ (relative) and the others' rms be below 10⁻⁸ of max f*. <br>• Planted ends [own]: <br>  – (1) a main kernel deficient by 5·10⁻⁷, made up by a formed ±0.3 pair at its θ: `outcome` says "exact, null cluster"; the summary says no (others' rms 7.4·10⁻⁸). <br>  – (2) main kernels off by 5·10⁻⁷ and a far cancelling pair: `outcome` says "not fitted" (residual 1.1·10⁻⁷); the summary says null cluster. <br>• On the data the two agree on their only row (H σ* 0.04). | Use **one predicate** (main kernels within EXACT_TOL, residual below 10⁻⁸, others' rms below NULL_TOL·max f*) in both `outcome` and `split_summary`. Add tests with the two planted ends. Fix the comment and "the residual's threshold". |
| D3cu-4 | minor | **The tests.** <br>• **Totals.** The repository's 26 round-97 mutants: all killed. Mine: 40, of which 30 killed and 10 survive. <br>• **Behaviour-changing survivors:** <br>  – AE2 (`iw[13]` → `iw[14]`) and AE3 (→ `iw[18]`) in `a2e_job` survive. On a rerun they rewrite the published A2E line: "order-1 steps at hλ 2.41–7.56 … −2.43 … −0.152", and "order-2 steps at hλ 2.41–7.56 … −2.54 … −0.152", the damping order-1 steps gone (jobmut [own]). <br>  – The new assert checks the raw, not the code, so the log's "(kills AE2, AE3)" is wrong. <br>• **Test line 1099.** The ratio-4 checks "and cr > 1e-5 and cy < 1e-14" sit inside the comment. <br>• **Equivalent on the data:** <br>  – R19 (the null test's main tolerance 10⁻³): untested. <br>  – R17 (max f* without abs): equivalent on H's positive targets. <br>  – R01, R07, R09, NX1 (boundaries) and R25, R26 (guards). | In `test_a2e_the_early_collapse`, which already runs `a2e_job(10.0)`: assert its symmetric-state BDF orders == [2, 2, 2, 1, 1, 1, 1, 1, 1]. That kills AE2 and AE3. <br>Move the ratio-4 conditions out of the comment. <br>Add a planted end with main error 5·10⁻⁴ and others' rms below 10⁻⁸ (kills R19). <br>Correct the log entry. |
| D3cu-5 | cosmetic | **The cap's gap is right but incompletely justified** (md:617; D59.md:47–49). <br>• Brute force: no assignment changes for any cap in (1.1498, 2.8144) on 94 of 94 [own]. <br>• near_cs looks only at unassigned features. A same-sign kernel 1.805 off an assigned feature lies inside (1.5, 2.81) and is excluded by the cap; it changes nothing only because the min–max objective prefers the main kernel. <br>• Only 29 of the 30 unassigned features have an unused kernel. H3 ±0.05's second feature is unassigned by distinctness: its sign's only kernel, (0, 0.129), w 1.42, is the first feature's main kernel, 1.15 from both. <br>• That budget stop is the 1.15 maximum. Otherwise the maximum is 0.818. | Print the brute-force interval of caps over which no assignment changes. Write "29 of the 30 (one left by distinctness)". Name the 1.15 row as a budget stop's kernel shared by both features. |
| D3cu-6 | cosmetic | **Smaller wording points.** <br>• "whose f* is 0" (D59.md:21; md "0.0e+00") is the double's underflow; the exact value is about 5·10⁻⁶⁰³·A. <br>• "the residual's threshold" compares an absolute rms with an rms over max f*. <br>• The README's round-83 text still says "(κ′ − κ)/ε^{1/3} settling on targets 9, 27 and 29" (D3cr-9) and repeats "κ′ − κ (κ′ − κ)/ε^{1/3}". | Write "≈ 5·10⁻⁶⁰³ (0 in double)" and "the same number". In the README, mark the round-83 phrase as revised ("29 falling in size since 10⁻⁸"). |

## Sources

- **[snippet]** Invariance of fixed-point sets under equivariant flows:
  - [arXiv:2206.04442](https://arxiv.org/pdf/2206.04442) states that Fix(Γ) is invariant for Γ-equivariant systems.
  - [ChaosBook lecture notes (Kohler)](https://www.chaosbook.org/projects/Kohler/workspace.pdf): "the action of the flow leaves them invariant".
  - **[memory]** Golubitsky–Stewart–Schaeffer, *Singularities and Groups in Bifurcation Theory* II, as the standard reference. The searches did not surface it.
- **[snippet]** Backward Euler's region |ζ − 1| ≥ 1 and its factor 1/(1 − hλ), so a positive eigenvalue is damped only for hλ ≥ 2. Search summary over [Driscoll, FNC, "Absolute stability"](https://tobydriscoll.net/fnc-julia/diffusion/absstab-diffusion.html) and [Purdue lecture 4](https://www.math.purdue.edu/~qi117/Spring26/docs/lecture4_stability.pdf).
- **[snippet]** Gear–Kevrekidis, "Reverse integration for computing stationary points of unstable stiff systems", [arXiv:nlin/0302055](https://arxiv.org/pdf/nlin/0302055): strongly damped components coupled with mildly growing ones.
- **[snippet]** Round-off-level symmetry errors amplified into broken symmetric solutions:
  - [arXiv:2106.01872](https://arxiv.org/pdf/2106.01872), which traces them to the non-associativity of floating-point arithmetic.
  - [Math. Comput. Simul. 37(4), 323–339 (1994), IDEAS record](https://ideas.repec.org/a/eee/matcom/v37y1994i4p323-339.html): roundoff breaking a symmetry and unleashing a symmetry-breaking instability.
- **[memory]** The draft's other references are unchanged this round. Round 97 confirmed Ho–Nguyen and Xu–Du [snippet in that round].

## Lesson

"Invariant up to 10⁻³⁸" is not invariance when the state is a saddle. The asymmetry's logarithm divided by λ is a finite time, here s ≈ 0.02–4.2. A symmetry argument proves invariance only for the exactly symmetric problem, and the problem as posed should be stated as a perturbation of it. Before calling a held state "the exact flow", check two things:
- whether the solver resolves the unstable direction (hλ against BDF-k's bound);
- whether the held offset is the implicit step's fixed point −d/λ of a drift. Here it was, to two digits.

Two related lessons:
- A residual a filter calls "roundoff" may be deterministic; replaying the step with the drift separates the two.
- An assertion on stored data does not kill a mutant of the code that writes that data. Test the job's output when the job is cheap.
