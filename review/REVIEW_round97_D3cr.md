# Referee report D3cr, round 97: Proposition D‴.59 as revised after round 87

## Scope

- **Export.** Commit 0ea58aa (branch claude/vigilant-clarke-hf939b), exported with `git archive 0ea58aa | tar -x` into W/repo. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref97.
- **Base for the diff.** The round-87 referee reviewed commit 7be2cb0 ("D‴.59 round 83 …"). It is the last commit before 66cc9e1 that touches D‴.59's files; 32a1144, the round-87 report, touches only the report. All diffs are 7be2cb0 → 0ea58aa, saved in W/diffs/.
- **Read in full:**
  - review/drafts/D59.md (340 lines).
  - The 81 changed lines of experiments/exp_moving_width.py, line by line: A2_OFF, A2E_MORE, `a2_job`, `a2e_job`, `part_A2E`, `bdf_bound`, `a2_lines`, SPLIT_MAX/`fit_split`/`split_cells`, HLL's ratio line, and `g_fine_lines`' over-ε^{1/3} line. Also the code these call or feed: `sq_job`, `sq_args`, `sq_check`, `g_fine_lines`, `gl_lines`/`gl_constant`, the H/HP/HL/HLL parts of `report`, `critical_rhs_nd`/`critical_jac_nd`, and `BumpsCW`.
  - results/chaos/moving_width.md (614 lines) and tests/test_moving_width.py (1074 lines, 48 tests).
  - review/REVIEW_round87_D3ch.md, the Round 87 entry of review/REVIEW_LOG.md, and the D59 row of review/drafts/README.md.
- **Settings.**
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` and the repository's venv (numpy 2.5.4, scipy 1.18.1, matching tools/requirements-pinned.txt).
  - MW_OUT pointed to W/mwout, or to the mutant slots' own `out` directories, each holding a copy of the raw.
  - pytest ran with `--basetemp` W/pyt. No test uses a temporary-directory fixture, so the mutant runner's pytest, which has no `--basetemp`, wrote nothing outside W.
  - At most two heavy processes of mine ran at a time.
- **I modified nothing in the repository.**
  - I used only read-only git: `archive`, `log`, `show`, `diff`, `status`, `rev-parse`.
  - During the review another session committed 99edc3d and bc0c89b (D‴.60 round 94, handoff). They touch none of D‴.59's files.
  - At the end D59.md, the script, the md, the raw and the tests are byte-identical to the export, and nothing of mine is running.
  - One wrong step: waiting loops of mine expected 51 mutants instead of 50, so they would never have ended. I stopped them through the task tool, by task ID.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my run or computation; [snippet] a web-search result; [memory] not verified this round.

## Summary

**Design.** I worked in this order:

1. Reproduced the md from the raw and ran the tests.
2. Reran jobs against the raw.
3. Recomputed A2 with my own right-hand side, using five integrators and starts from 0 to 10⁻⁶σ*.
4. Tested whether the symmetric line c = c* is invariant.
5. Read A2E's BDF steps from the raw and computed the BDF-k roots myself.
6. Rebuilt `fit_split` by brute force on all 94 rows, and computed the others' mass, both over ℝ and over the samples' window.
7. Recomputed κ′ − κ and (κ′ − κ)/ε^{1/3}, and rebuilt SQ's 14 off-sequence cells from the event lists.
8. Ran 50 single-point mutants of this round's code, plus 6 analogues of round 87's survivors.
9. Did six web searches.

**Verdicts.**
- Reproducibility is exact: the md byte for byte, 48 of 48 tests, and every rerun identical to the raw.
- D3ch-1's numbers are fixed: A2's 3A at σ*/σ₀ = 0.25 now reaches the plateau to 3·10⁵, escapes by 10⁶ and ends exact. My own LSODA and Radau from ±10⁻⁹σ*, 10⁻¹²σ* and 10⁻⁶σ*, and DOP853 and RK45, all confirm it.
- The BDF-k bound Σ_{j≤k}2^j/j is correct; I checked it for k ≤ 6.
- (κ′ − κ)/ε^{1/3} and the SQ counts reproduce exactly.
- What remains is interpretation and criteria, not numbers.

**Findings.**
1. **D3cr-1 (minor): the exact start lies on an invariant line, so the labels are inverted.** At σ*/σ₀ = 0.25 and 1 the line c = c* is invariant for the problem as posed. The state held there is the exact flow from (c*, σ*), not "a solver artefact". The integrators that leave it are amplifying roundoff of about 10⁻¹⁶.
2. **D3cr-2 (minor): the others' "signed mass" counts mass outside the samples' window.** In the symmetric 10⁻⁸ cell the mass +0.184 → +0.179 lies almost entirely outside [−1.6, 1.6]. Inside the window the others' mass falls +2.8·10⁻² → +1.0·10⁻⁵.
3. **D3cr-3 (minor): the one-norm null-cluster test contradicts H's own "exact, null cluster".** For σ* 0.04 the others' rms is 29.5 times the residual's, and the test has no tolerance (D3ch-3(c) remains open).
4. **D3cr-4 (minor): SPLIT_MAX's stated justification is wrong.** The stated window is "≤ 2" for cancellations and "3–7" for the dropped kernels. The data show cancellations up to 2.24 and dropped kernels from 2.81. A cap on the (c, σ) error alone would separate the classes with a wide margin.
5. **D3cr-5 (minor): two order-recording mutants pass all 48 tests.** On a rerun each rewrites the published BDF-order line, which is this round's D3ch-2 fix. Four further survivors are equivalent on the data.
6. **D3cr-6 (minor): the sentence on target 18 carries target 12's mass range and bracket** (an editing slip).
7. **D3cr-7 (cosmetic):** "E 1.8–2.6·10⁻³²" is stale; the current raw gives 2.2–2.8·10⁻³².
8. **D3cr-8 (cosmetic):** SQ has five landed cells among the kernel-9 cells, not four; target 20's kernel 16 forms before τ_d.
9. **D3cr-9 (cosmetic):** the "−0.00401" multiplier is a roundoff ratio; 29's ratio is quoted by its endpoints only; the README's "settling" is generous for 29; a test comment claims more than its assert.

## What I checked and how

All paths are under W.

1. **Reproducibility.**
   - `--from-raw` into W/mwout, then `cmp` with the committed md: identical (0.6 s).
   - pytest: 48 passed in 14.75 s (W/logs/pytest_full.log).
   - W/scripts/rerun_a2.py: all 15 `a2_job` and all 3 `a2e_job` reruns are bit-identical to the raw in every stored field (W/logs/rerun_a2.log).
   - W/scripts/rerun_other.py: `sq_job` for target 20 (45 s) and `hl_job` for "σ* 0.06, A₂ 0.95" (38 s) are identical too.
   - W/scripts/explore.py: against 7be2cb0's raw (from `git show`), only A2 and A2E changed.
   - The other 14 A2 cells changed beyond roundoff in a few fields. The plateau side flips (dc ±0.165σ* at 0.25, ±0.083σ* at 1), and w moves by up to 10⁻⁴ near the escapes. The printed durations, escapes and smallest σ are unchanged [own].
2. **A2 recomputed with my own code.**
   - W/scripts/own_a2.py: my own one-kernel right-hand side on the critical curve, with a finite-difference Jacobian.
   - For 3A at 0.25 (W/logs/own_a2_3a.log):
     - LSODA from ±10⁻⁹σ*, 10⁻¹²σ* and 10⁻⁶σ*: the plateau 1 … 3.16·10⁵, escape 10⁶, smallest σ 0.1111σ*, exact end.
     - LSODA and scipy BDF from 0 and from 10⁻¹⁵σ*: stay symmetric (smallest σ 0.2346) and end exact.
     - Radau stays from 10⁻¹⁵σ* but leaves from 0 (its asymmetric finite-difference Jacobian).
     - DOP853 (rtol 10⁻¹²) and RK45 from 0: leave by s = 0.1 (dc −0.1646, σ 0.1111σ*, E 0.6436).
     - BDF from ≥ 10⁻¹²σ* reaches the plateau but fails at the escape (status −1).
   - All 15 cells under Radau from +10⁻⁹σ*, −10⁻⁹σ* and +10⁻¹²σ* (W/logs/own_a2_all.log): the md's durations, escapes, smallest σ and 20A ends (w 18.5905, 17.5957, 17.4894). Radau fails at the 10A escape at all three ratios.
3. **Invariance.** W/scripts/invariance.py (W/logs/invariance.log) evaluates ċ at c = c* on the repository's samples and on samples exactly symmetric about c* (see D3cr-1).
4. **A2E's steps and the BDF bound.**
   - W/scripts/a2e_steps.py (W/logs/a2e_steps.log): each step's method, order, hλ and multiplier, from the raw.
   - W/scripts/bdf_bound.py (W/logs/bdf_bound.log): roots of BDF-k's characteristic polynomial on z ∈ (0, 30].
     - For k = 1 … 6 the unstable set is exactly (0, Σ_{j≤k}2^j/j] = 2, 4, 6.667, 10.667, 17.067, 27.733. Just beyond it the largest root has modulus < 1.
     - BDF2's roots at z = 2.41 are −2.424 and 0.227.
5. **Split, mass and norms.**
   - W/scripts/split_check.py (W/logs/split_check.log): brute-force assignment with and without the cap on all 94 rows of H, HP, H3, HL and HLL. `fit_split` equals my assignment on 94 of 94, and `rest_mass`/`rest_rms` match.
   - The same script gives the others' mass over ℝ, over the window, on the samples (Σ·Δx) and at the features.
   - W/scripts/sym_rest.py: the symmetric 10⁻⁸ cell kernel by kernel.
   - W/scripts/cap_parts.py (W/logs/cap_parts.log): each main kernel's error split into its (c, σ) part and its w part.
6. **κ′ and SQ.**
   - W/scripts/kappa_rate.py: κ′ = v_land/(v(τ₃⁻)((τ₃^path − τ_d)/(t_v3 − τ_d))^{1/3}) from the raw's fine phases. It reproduces every printed κ′ − κ, ratio and over-ε^{1/3} value.
   - W/scripts/sq_types.py (W/logs/sq_types.log): first-formation times per kernel for the 14 off-sequence cells, against τ_d, j3 and j4.
7. **Mutants.**
   - Tool: W/scripts/mutants_tool.py, which is tools/mutants.py with REPO set to the export and PY to the venv.
   - W/scripts/spec_r97.py: 50 mutants of this round's code, run in 2 slots (W/mut/res_a.jsonl, W/mut/res_b.jsonl, W/logs/mut_batch1_summary.log).
   - W/scripts/spec_r97b.py: 6 analogues of round 87's SA3, SJ4, AE12, GL1, FS12 and SC3 (W/logs/mut_batch2_summary.log).
   - W/scripts/jobmut.py: reruns A2E inside a mutated copy and diffs the rendered md. A no-op control gives 0 changed lines.
8. **Numbers.** W/scripts/numcheck_own.py (tools/numcheck.py on the saved diff): every number on the added lines is in the md. A whole-draft pass plus a manual check found the stale range in D3cr-7.
9. **Literature.** Six web searches; no page was fetched.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) landscape, convergence from E(start) < ½‖f*‖², proof (8–15, 126–134) | holds, proved in the continuum, unchanged | re-read: the energy bound keeps (c, σ) compact and w away from 0; Łojasiewicz applies for the analytic E in a non-degenerate metric |
| (i) ∂_σE = ½π^{1/2}w(w − A) at (c*, σ*) (17) | holds | own derivation |
| (i) "from 3A on, at every width ratio, it collapses at once …; from c* exactly LSODA holds 3A at 0.25 there, so A2 starts at c* + 10⁻⁹σ*" (16–21) | numbers hold; **the attribution does not** (D3cr-1) | [own] invariance; integrators |
| A2: from (c* + 10⁻⁹σ*, σ*); 2.5A dips to 0.295–0.304σ*; plateau from the first output, σ to 0.00485–0.111σ*; E₂, E₁; durations; escapes; 20A (228–232) | holds | [md], [raw], own Radau/LSODA/DOP853/RK45 |
| "all twelve end exact (E 1.8–2.6·10⁻³²)" (231–232) | **stale** (D3cr-7) | [raw]: 2.233–2.759·10⁻³² |
| A2E: the symmetric state, ∂ċ/∂c, Adams growth 1.092·10⁴, plateau from 0.0056 (232–235) | holds | [raw] |
| A2E: order-2 steps at hλ 2.41 ×−2.54 … −1.67; order-1 at hλ 3.99–49.5 ×−0.335 … −0.00401 (235–237) | holds as read off the raw; the last value is roundoff (D3cr-9) | [raw]: −0.3348, −0.335, −0.1521, −0.1554, then −0.004008 at an offset of 5.8·10⁻¹⁴ → −2.3·10⁻¹⁶ (backward Euler's factor there is −0.0206) |
| BDF-k damps a positive real eigenvalue only beyond 2, 4, 6.667 (237) | holds for constant steps | [own], k ≤ 6. LSODA's variable-step BDF is close (observed −2.43 against the root −2.42) |
| A2E's 3A exact-start run: \|c − c*\| ≤ 6.3·10⁻¹³σ* to s = 10, σ to 0.2345σ*; "A2's route is the flow's, the span-10 run's a solver artefact" (238–240) | numbers hold; **the labels are inverted** (D3cr-1) | [raw], [own] |
| (ii)(b): the split's main kernel within 2.5; the others' signed mass; the null cluster's rms (43–48) | computation holds (94/94); **the mass is out of window** (D3cr-2); **the null test contradicts H** (D3cr-3); **the cap's justification is wrong** (D3cr-4) | [own] |
| asymmetric cells' masses −0.253 → −0.271 etc.; HP, HL, HLL numbers (49–56, 305–332) | holds | [md]; in-window mass equals continuum mass here |
| symmetric 10⁻⁶: mass +0.108 → +0.040, "a split" (58–61) | holds | in-window, [own] |
| symmetric 10⁻⁸: "the mass +0.184 → +0.179 … by their mass a split too" (62–63) | sign holds; **the magnitude is outside the window** (D3cr-2) | in-window +2.8·10⁻² → +1.0·10⁻⁵ [own] |
| "the others' rms over the residual's 2.49·10³ … 1.33·10³, falling and decelerating: undecided at 10⁹" (64–66) | holds | [md], [own] |
| (3) κ′ − κ per decade and over ε^{1/3} for 9, 20, 27, 29 (101–104) | holds; 29's quote hides the peak (D3cr-9) | [own] = [md] |
| "on 18 κ′ − κ changes sign over 10⁻⁷ … 10⁻¹², where the loop's own crossing still puts its κ 0.25–0.35% above the path's" (105–106) | **misattributed** (D3cr-6) | [md] |
| SQ: 34/48; the 14 by type; "four … with a landing"; "nine … no landing"; proxy misses three (277–284) | counts hold; **"four" should be five**, and 20's type is mislabeled (D3cr-8) | [own] sq_types |
| Status (178–192) | holds, but inherits D3cr-1 and D3cr-2 | — |
| README D59 row | fits, except "the symmetric 10⁻⁸ cell a split" (D3cr-2) and "settling on … 29" | — |
| mathematics: what is proved | (i) proved in the continuum; (ii)(a) derived; (ii)(b) numerical; (c)'s law derived on the slow manifold, with the approach numerical; the four stages derived from the path and measured; (iii) numerical; A2/A2E numerical; the BDF bound is a constant-step fact (the root −1 at z = Σ2^j/j), its real-axis extent checked numerically, not proved | — |
| literature (193–215) | unchanged; Ho–Nguyen and Xu–Du confirmed [snippet] | Sources |

## Status of the round-87 issues

| D3ch | status | what remains |
|---|---|---|
| 1 (A2's 3A at 0.25 on the saddle) | **fixed** (numbers, AS6 inverted, A2E_MORE) | The draft's "LSODA holds … a solver artefact" mislabels an invariant line (D3cr-1). The stale end range (D3cr-7). |
| 2 (A2E's damping) | **fixed** (orders recorded, multipliers printed, `bdf_bound` correct) | The order-recording mutants survive (D3cr-5). The −0.00401 value (D3cr-9). |
| 3 (the split criterion) | **partly fixed** | (a) The classes by the mass: the sign is right, but the mass counts out-of-window mass (D3cr-2). (b) Σw's flips: moot. (c) The null test now uses one norm but has no tolerance, and contradicts H's label (D3cr-3). (d) Broad kernels: fixed by SPLIT_MAX, but its justification is misquoted (D3cr-4). |
| 4 ("approaching a null cluster") | **fixed** | — |
| 5 (tests) | **fixed** for the named gaps: the SA3, SJ4, AE12, GL1, FS12 and SC3 analogues are all killed, and the full-md test exists | New gaps (D3cr-5). The GL test comment overclaims (D3cr-9). |
| 6 (SQ's 14 by type) | **mostly fixed** | Five landed, not four; 20's kernel 16 forms before τ_d (D3cr-8). |
| 7 (rate wording) | **fixed** | The 18 sentence slip (D3cr-6); 29's endpoints (D3cr-9). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cr-1 | minor | **The exact start is on an invariant line: the held state is the exact flow, not a solver artefact (D59.md:16–21, :238–240; md A2E's last line; README).** <br>• **Invariance at 0.25 and 1.** At σ*/σ₀ = 0.25 and 1 the samples are symmetric about c* within the kernels' reach: the unpaired samples carry f*/A = 0 (0.25) and 2.3·10⁻³⁸ (1). So c = c* is invariant for the problem as posed. <br>• **ċ on the line is roundoff.** On the repository's samples \|ċ\| ≤ 1.8·10⁻¹² (linspace is symmetric about c* only to 4.4·10⁻¹⁶). On exactly symmetric samples \|ċ\| ≤ 6.4·10⁻¹⁶ [own]. <br>• **Ratio 4 differs.** The window's unpaired samples (f*/A 4.4·10⁻³) drift c, with ċ = +8.6·10⁻⁵ at the start, so there the exact start is generic. <br>• **What the solvers show.** From c* exactly, LSODA and scipy BDF (my code) keep the symmetric solution, as does Radau from 10⁻¹⁵σ*. DOP853 and RK45 leave by s = 0.1, by amplifying roundoff through the saddle. <br>• **Consequence.** The plateau is the behaviour of every start off c = c*, from 10⁻¹²σ* in every integrator I tried. It is not the behaviour of (c*, σ*) itself at 0.25 and 1. Nothing printed changes: the exact start also ends exact for 3A. | **(i):** "started at (c*, σ*)" → "started off the line c = c* (A2: c* + 10⁻⁹σ*)". <br>**(i) and A2E:** state that c = c* is invariant at σ*/σ₀ = 0.25 and 1, and a saddle in c on the symmetric state; LSODA's held run is the exact solution from that measure-zero start; DOP853/RK45 leave by amplifying roundoff. <br>**A2E:** replace "A2's route is the flow's, the span-10 run's a solver artefact" with "A2's route is the generic flow's". <br>**md and README:** "the symmetric state held by LSODA" → "the invariant symmetric state, which LSODA follows". |
| D3cr-2 | minor | **The others' mass counts mass outside the window (D59.md:45, :62–63; `fit_split` `rest_mass`; README "the symmetric 10⁻⁸ cell a split").** <br>• **The symmetric 10⁻⁸ cell (HLL).** The others' continuum mass is +0.184, +0.164, +0.169, +0.175, +0.179 at τ = 10⁵ … 10⁹. Inside [−1.6, 1.6] it is +2.8·10⁻², 2.1·10⁻³, 2.7·10⁻⁴, 4.6·10⁻⁵, 1.0·10⁻⁵ (the samples' Σ·Δx agrees). <br>• **What +0.18 is made of.** It is the difference between one broad kernel at c = 0 (σ 1.59 → 1.95; mass +1.78 → +2.30, only +1.22 → +1.35 inside) and ten negative ones (−1.59 → −2.12), drifting to \|c\| 2.13. <br>• **HP's equal end** is the same: +0.184 over ℝ, +0.0285 inside. <br>• **No class flips:** the sign agrees on all 94 rows [own]. But "the mass +0.184 → +0.179 … a split too" shows a steady O(0.2) split. On the samples the others' mass falls ×7.5–13 per decade with the main kernels' deficit. | **`fit_split`:** compute the mass over the window (the erf form), or as Σ_x oth(x)Δx, and report both. <br>**(ii)(b) and the HLL lines:** quote the in-window mass (+2.8·10⁻² → +1.0·10⁻⁵, a split of vanishing size). <br>**README:** adjust the row. |
| D3cr-3 | minor | **The null cluster's one-norm test contradicts H's "exact, null cluster" (D59.md:47–48, :302–303; md H row σ* 0.04).** <br>• **The H row.** The others (three kernels at θ ≈ (0, 0.399), Σw 0, Σ\|w\| 0.219) have rms 3.9·10⁻¹² of max f*, against the residual's 1.3·10⁻¹³: ratio 29.5 [own]. By (ii)(b)'s "others' rms … below the residual's" it is not a null cluster, yet `outcome` and the Numerics call it one. <br>• **Why the test fails.** The main kernels absorb the cluster's sum, so the residual can lie below the others' rms even for an exactly cancelling cluster. <br>• **No tolerance was added** (D3ch-3(c)). | **(ii)(b):** a test relative to max f* with a tolerance — e.g. others' rms ≤ 10⁻⁸ of max f*, or ≤ 10× the residual's rms — or the same test as `outcome`'s. <br>**md:** print `rest_rms` in the H/HL tables. <br>**H, (ii)(b):** reconcile the σ* 0.04 label. |
| D3cr-4 | minor | **SPLIT_MAX's justification is wrong and its margin thin (code comment at `SPLIT_MAX`; log D3ch-3: "above a cancellation's inflation (≤ 2 here), below the transient broad kernels' errors (3–7)"; D59.md:44).** <br>• **Kept side.** The cancellations' main kernels have uncapped errors up to 2.239 (H equal), 2.163 (σ₂ 0.065) and 2.144 (A₂ 0.95). <br>• **Dropped side.** The newly unassigned rows start at 2.814 (HL σ₂ 0.065, 10⁻⁸, τ 1600) and 2.867 (HL A₂ 0.9, 10⁻⁸, 6400; a displaced kernel, w/A − 1 = +0.25), then 3.13–7.23. So the cap sits in the data gap (2.24, 2.81). <br>• **The md cannot see it.** Mutant SM1 (cap 2.3) leaves the md unchanged; only the literal `SPLIT_MAX == 2.5` assert kills it. <br>• **Clean separation by kind.** Every cancellation's error lies in w, with its (c, σ) part ≤ 0.82 (1.15 with H3's budget stop). Every dropped row's error lies in (c, σ) (≥ 2.81), with \|w/A − 1\| ≤ 0.28 [own]. <br>• **Risk.** A growing cancellation (A₂ 0.99 +0.05–0.1 per decade; H's polish ends already at +2.2) would lose its main kernel beyond 2.5. | **`fit_split`:** cap the (c, σ) error only (e.g. 1.5), leaving w free, and state the observed gap (≤ 1.15 against ≥ 2.81). Or keep 2.5 and quote the true gap (2.24, 2.81). <br>**Comment and log:** fix them. <br>**Tests:** a test that an inflated main kernel (w/A − 1 = 3) with exact (c, σ) stays the main kernel. |
| D3cr-5 | minor | **Tests: order-recording mutants survive and change published results on a rerun.** <br>• **Totals.** 50 mutants of this round's code: 44 killed by the tests, none caught by the md alone, 6 survive both. <br>• **Behaviour-changing:** <br>  – AE2 (`iw[13]` → `iw[14]`, NQCUR for NQU). On a rerun the md's A2E line becomes "BDF order-1 steps at hλ 2.41–49.5 multiply the offset by −2.43 … −0.00401": an amplifying step counted as order 1. <br>  – AE3 (`iw[13]` → `iw[18]`, the method for the order). The rerun md reads "order-2 steps at hλ 2.41–49.5 … −2.54 … −0.00401", and the damping order-1 steps vanish (W/scripts/jobmut.py). <br>  – The tests check only 1 ≤ order ≤ 12 and BDF ≤ 5. <br>• **Equivalent on the data:** AE4 (`off` recorded as 0 in `a2e_job`: every A2E run has off = 0), AL11 (E₁ for E₂ in the exact-start line: its E stays ≥ 1.19), SM3 (`<` at the cap), FM2 (signed σ in the mass). <br>• **Round 87's analogues** (SA3, SJ4, AE12, GL1, FS12, SC3): all killed now. | **`test_a2e_the_early_collapse`:** assert the raw's span-10 orders on the symmetric state, i.e. (2, 2, 2, 1, 1, 1, 1, 1, 1) at s 0.00264 … 0.0143, or that a rerun's steps equal the raw's. <br>**Tests:** add one with an "other" kernel of negative width (FM2). |
| D3cr-6 | minor | **Editing slip in (3) (D59.md:104–106).** "and on 18 κ′ − κ changes sign over 10⁻⁷ … 10⁻¹², where the loop's own crossing still puts its κ 0.25–0.35% above the path's (GT)". <br>• 18's on-sequence masses are 10⁻⁸ … 10⁻¹² (10⁻⁷ is †). <br>• 18's bracket [0.7498, 0.7525] holds κ = 0.7520 [md]. <br>• The range and the 0.25–0.35% are target 12's: at 7be2cb0 the clause followed "on target 12 only ×1.32–1.85". | Move "over 10⁻⁷ … 10⁻¹², where the loop's own crossing … (GT)" back to target 12. Write "on 18 κ′ − κ changes sign over 10⁻⁸ … 10⁻¹² (+0.0035 → −0.0002)". |
| D3cr-7 | cosmetic | **Stale end energies (D59.md:231–232).** "all twelve end exact (E 1.8–2.6·10⁻³²)" is round 83's raw (1.781–2.627·10⁻³²). The current raw and md give 2.233–2.759·10⁻³². numcheck missed it because the line is unchanged. | Write "2.2–2.8·10⁻³²". Run numcheck over every line of the draft that quotes a changed part. |
| D3cr-8 | cosmetic | **SQ by type (D59.md:279–284).** <br>• **Five landings, not four.** Of the kernel-9 cells, five have a landing: 18 at 10⁻⁶ lands at τ = 2785.76 (intermediate slope −21.3, κ "—") besides the four quoted. Only five matches "nine of the 14 have no landing". <br>• **Target 20's type.** At 10⁻⁵ and 10⁻⁶ kernel 16 forms at 24.2 and 29.5, before τ_d = 32.22, and the path's j3 = 4 never forms. These are the same type as 18, 27 and 12 at 10⁻⁵ and 10⁻⁶ (the loop's τ₃ only the first output past τ_d), not just "an extra kernel" [own]. | "at the five with a landing (four with a κ) …"; "kernel 16 in place of the path's 4, formed before τ_d (20 at 10⁻⁵, 10⁻⁶)". |
| D3cr-9 | cosmetic | **Smaller wording points.** <br>• **A2E's last multiplier (D59.md:237).** "−0.00401" is the ratio as the offset drops into roundoff (5.8·10⁻¹⁴ → 2.3·10⁻¹⁶); backward Euler's factor at hλ 49.5 is −0.0206. <br>• **Target 29 (D59.md:103).** "−12.6 → −21.1" hides the peak −25.9 at 10⁻⁸. The evidence for O(ε^{1/3}) on 29 is the fall since then (×0.94–0.96 per decade). <br>• **README.** "(κ′ − κ)/ε^{1/3} settling on … 29" is generous. <br>• **GL test comment.** "≥ 1% of E at every row" holds only for the last row, the only one asserted: the 10⁻²¹ floor is 0.013%, 0.34% and 48% of E at the three rows. GL1 is killed anyway. | "×−0.335 … −0.155 (then roundoff)"; "−12.6 → −25.9 (10⁻⁸) → −21.1"; "near-constant on 9 and 27, drifting −5% per decade on 29"; "at the last row" in the test comment. |

## Sources

- **[own]** The BDF-k unstable set on the positive real axis is (0, Σ_{j≤k}2^j/j] for k = 1 … 6. These are numerical root computations, not a proof that the boundary locus crosses the real axis only at ζ = ±1. I found no source stating the formula; the search results covered general BDF stability regions only [snippet]:
  - [arXiv:2412.19959](https://arxiv.org/pdf/2412.19959)
  - [kfupm 377.pdf](https://faculty.kfupm.edu.sa/MATH/kabir/OLD/377.pdf)
  - [NCSU MA583 ch. 5](https://mtchu.math.ncsu.edu/Teaching/Lectures/MA583/chapter5.pdf)
  - [Akrivis–Katsoprinakis](https://www.cs.uoi.gr/~akrivis/AK2a.pdf)
  - [Lóczi](https://www.math.u-szeged.hu/12QTDE/abspdf.php?tid=69)
  - A result list cited Hairer–Wanner, "On the instability of the BDF formulas" (SIAM J. Numer. Anal. 1983), only as a reference (BDF unstable for k ≥ 7) [snippet].
- **[snippet]** Backward Euler's amplification factor 1/(1 − λΔt): [AM111 session 8](https://am111.readthedocs.io/en/latest/session8_stabiilty.html), [LibreTexts 10.03](https://phys.libretexts.org/Bookshelves/Mathematical_Physics_and_Pedagogy/Computational_Physics_(Chong)/10%3A_Numerical_Integration_of_ODEs/10.03%3A_Backward_Euler_Method). I found no source on implicit solvers holding unstable symmetric states.
- **[snippet]** Ho–Nguyen, Ann. Statist. 44(6), 2726–2755 (2016), DOI 10.1214/16-AOS1444: [Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-44/issue-6/Convergence-rates-of-parameter-estimation-for-some-weakly-identifiable-finite/10.1214/16-AOS1444.full).
- **[snippet]** Xu–Du, COLT 2023, PMLR 195:1155–1198, Θ(t⁻³) for an over-parameterised student of one ReLU neuron, via a potential on the students' pairwise distances: [PMLR](https://proceedings.mlr.press/v195/xu23a.html), [arXiv:2302.10034](https://arxiv.org/pdf/2302.10034).
- **[memory]** The draft's other references (de Boor, Kammler, Braess, Petersen–Raslan–Voigtlaender, Gallon–Jentzen–Lindner, Chizat, Dwivedi et al., Schiebinger–Robeva–Recht, Denoyelle–Duval–Peyré–Soubies, Duval–Peyré, Golub–Pereyra, Kiefer–Wolfowitz, Day, Biernacki–Chrétien, Marr–Hildreth, Runnalls, Poggio–Girosi), unchanged this round.

## Lesson

Before calling one solver run an artefact and another the flow, ask whether the disputed state is invariant for the exact problem. Here the cheapest test was one evaluation of ċ at the start, on samples exactly symmetric about c*. It shows that at σ*/σ₀ = 0.25 the "held" run is the exact solution from a measure-zero start, and that the runs that leave are amplifying roundoff. Both rounds 87 and 97 rightly moved A2 off the line. The draft should say that it reports the generic flow, not that LSODA misbehaved. Two related lessons:
- A classification quantity should be computed where the fit lives. The others' mass over ℝ is dominated by kernel tails outside the samples' window; over the window it vanishes with the deficit.
- A cap should act on the coordinate that tells the classes apart. Cancellations differ from transient kernels in w against (c, σ). A combined error capped at 2.5 separates them only by the narrow data gap (2.24, 2.81), while a (c, σ) cap separates them by 0.82 against 2.81.
