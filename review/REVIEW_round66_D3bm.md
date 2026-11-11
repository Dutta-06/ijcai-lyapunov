# Referee report D3bm, round 66: Proposition D‴.54 as revised after round 62

## Scope

**What I reviewed.** Commit 6d48964 on branch claude/festive-bardeen-eopej2:
- `review/drafts/D54.md`: (iv)'s exception; (v); the proof's relative-mode paragraph; the Status; Classical/New; the Numerics paragraphs C4, C4W, C4M, R and RP.
- `experiments/exp_transit_nd.py`: `relative_law`, `centroid_velocity`, `rp_job`, `rp_lines`, `c4_job`, `c4_atol`, `c4m_lines`, `c4w_lines`, `window`/`edges`, `race_job` and `--from-raw`.
- `results/chaos/transit_nd.md` and `transit_nd_raw.npz`.
- `tests/test_transit_nd.py` (46 tests).
- `review/REVIEW_round62_D3bi.md` and the "Round 62" entry of `review/REVIEW_LOG.md`.
- Supporting code: `src/core/rich_transit_nd.py` (`loop_nd`, `critical_jac_nd`, `sliding_iss_nd`) and `src/core/rich_transit.py`.
- For comparison, round 62's raw, md and script from 86fbc02, read with `git show`.

**Working copy and settings.**
- Copy made with `git archive 6d48964 | tar -x` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref66/repo`.
- Scripts are in `…/ref66/scripts`, outputs in `…/ref66/out`, logs in `…/ref66/logs`, mutant copies in `…/ref66/mut/<id>`.
- Every run used the repository's `venv/bin/python` (numpy 2.5.4, scipy 1.18.1) with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`. `TN_OUT` pointed to scratch.
- At most two heavy processes of mine ran at a time. I killed one process, a duplicate of my own run, after checking `/proc/PID/cmdline`. None of my processes is still running.

**What I ran.**
- `--from-raw` and the test suite.
- Bitwise reruns of all 28 RP rows, both C4W rows and the C4M rows at 6.5·10⁻¹⁷, 3.7·10⁻¹⁷ and 10⁻¹⁷.
- A field-by-field comparison of the C4M rows with round 62's raw.
- An independent continuation of the relative mode below the floor: the exact first-order dynamics about a merged block. I injected its result back into the full race to get the continuum's s_d.
- Diagnostics of RP's "meas" forcing and of the slaving.
- 13 loops at other absolute tolerances, plus 3 loops with no window edge nearby.
- 34 single-point mutants, with rerun checks of the survivors.

**Repository untouched.** I wrote nothing under `/home/user/Lyapunov-Convergence-Kernels-Showcase`. I ran only read-only git commands (`archive`, `show`, `log`, `diff`, `status`).
- The live HEAD has since moved to 2391685, but `git diff 6d48964 HEAD` is empty for every D54 file.
- The working tree's uncommitted edits (exp_noise_resonance, exp_width_mismatch, REVIEW_LOG, D60 and their tests) are not mine.

**Literature.** The session's web-search budget was exhausted at my first query. The literature points rest on memory and are flagged as such.

## Summary

**What checks out**

*Reproducibility*
- `--from-raw` reproduces the md byte for byte (md5 959600ecc93f56783146845631f42944).
- All 46 tests pass (66 s).
- Fresh reruns equal the raw bitwise in every field except `secs`:
  - all 28 RP rows, including the far-edges rows at x = 5.48 on both targets;
  - both C4W rows;
  - C4M at (7, 6.5·10⁻¹⁷), (7, 3.7·10⁻¹⁷) and (7, 10⁻¹⁷).
- The C4M rows at ε ≥ 10⁻¹⁴ equal round 62's raw in every field except three:
  - `closest` and `d_held` depend on the output grid, which now starts at s = −100. For example, (7, 10⁻¹³) has `closest` 1.804958·10⁻⁸ against 1.804949·10⁻⁸.
  - The keys `early` and `win` are new.
- The C4/C4b raw rows no longer carry `e_t2`/`w_t2`.

*Mathematics.* I re-derived the following and they hold:
- **The identity.** dD/ds = Λ₀(2Gm ww_e/W)⟨r, ∇²φ(c̄)⟩D + (2(w − w_e)/W)dc̄/ds + O(|D|³).
- **The slaved form.** With c̄ = c_b + ρD, 1 − kρ = 1 + (w − w_e)²/W² = 1/κ.
- **The mobility matrix** in (c_b, D) is Λ₀Gm[[(w² + w_e²)/W², (w − w_e)/W], [(w − w_e)/W, 2]] (det = 1). Eliminating a fast centroid quasi-statically gives exactly κ(A + kρ′)D + κk·ċ_b*. So the κ form is an identity; slaving enters only through dc_b/ds.
- **The pin.** D* = −(A + kρ′)⁻¹k·ċ_b*, and κ cancels.
- **The pin's direction.** ∂c_b*/∂W = (∂F/∂c_b)⁻¹·½∇‖φ‖², where (∂F/∂c_b)⁻¹ is negative definite and ½∇‖φ‖² points inward. So the rest drifts toward the nearest edge, and with k < 0, D* points against the drift: "behind the drift, away from the nearest edge" is right. It matches RP:
  - top edge: drift +y, pin −y;
  - right edge: drift +x, pin −x.

*The continuum's side.* My independent check confirms it. It uses the exact c-block Jacobian, with no slaving and no extrapolated forcing.
- The result is toward on C2's target and far on the mirror, at x = 3, 4.9, 5.48, 6.04 and 7.
- The share along e never changes sign.
- It reproduces the run's D to 0.000 dex at s = 500, 700 and 900.

*C4W and C4M.*
- C4W supports "the nearest edge decides". It is bitwise reproducible, and unchanged at atol 10⁻¹⁵.
- C4M's far-rule window is (9.696/1.385·10⁹)² = 4.90·10⁻¹⁷ < ε < (9.696/5.71·10⁶)² = 2.88·10⁻¹² at x = 7. Its two boundary loops sit at s₂/s_{d,e} = 0.87 and 1.15.
- I confirmed the atol diagnosis at 10⁻¹⁷ independently:
  - at atol 10⁻¹⁴ the entrant is already fully formed (w = 0.795) at s = −5000, and the closest approach is 0.1085σ;
  - at atol 10⁻²² the entrant has w = 3.6·10⁻⁷ at s = −5000, the layer runs at s ≈ 0–1000, and the closest approach is 1.81·10⁻⁸σ.

*Numbers.* Every revised number I checked in (iv), (v), the Status, C4, C4W, C4M, R and RP matches the md and raw. The exceptions are in D3bm-2, D3bm-3 and D3bm-9.

**What is wrong, majors first**
- **D3bm-1 (major).** A held runner-up's s_d is not pin-independent: the hold starts when the climb ends, so s_d shifts one for one with the climb's end. In the continuum s_d is +12% at x = 4.9 (3.58·10⁵), +23% at 5.07 (1.85·10⁵, which exceeds C3's s₂ = 1.64·10⁵, so C3's capture is the window's) and +38% at 5.48. The draft's own top-0.95 row shows the effect.
- **D3bm-2 (major).** RP's continuation is right about the sides but wrong in its numbers:
  - The climb times use the run's base after that base has left the merged state. With a merged base they are 7.0–7.5·10⁴, not 2.99–6.86·10⁴.
  - The "measured relaxation" is the pair's own feedback, extrapolated. Without the extrapolation 10^{−45} becomes 10^{−103}.
  - The least |D| at x = 3 is 10^{−98.7}σ, not 10^{−103.5}σ.
- **D3bm-3 (minor).** R's rule's only disagreement, C4's (4.9, 10⁻¹¹), is an absolute-tolerance artefact. At atol ≤ 10⁻¹⁵ it is captured toward, which gives 74 of 74.
- **D3bm-4 (minor).** `test_c4_atol_below_the_masses` is vacuous because `np.isclose` uses a default atol of 10⁻⁸. 17 of 34 mutants survive, and 11 of those change the reruns.
- **D3bm-5 (minor).** The rate check is printed only where it is best. κ reaches 0.61 at x = 3, where the slaved rate is 5–9% fast. The slaving condition is not stated.
- **D3bm-6 (minor).** "In the continuum" holds for the ε-free inner problem in exact arithmetic, not for any float64 loop.
- **D3bm-7 (minor).** RP's table presentation:
  - the "unforced" rows on pinned windows inherit the pin's side;
  - "exact zero" hides two continuations that climbed and a NaN.
- **D3bm-8 (minor).** The overlap-singularity literature is missing, and "New" is overstated.
- **D3bm-9 (cosmetic).** Small numeric slips in the text.

## What I checked and how

**Reruns** (`rerun.py`, `cmp_c4m.py`)
- `rp_job`, `c4_job` (C4W and C4M) and the committed arguments are compared field by field with the raw, bitwise, ignoring `secs`.
- The C4M rows are also compared with 86fbc02's raw.

**Atol** (`atol_loop.py`, `run_atol_batch.sh`, `c4w_atol.py`)
- These run `c4_job`'s loop at a chosen atol with outputs from s = −5000, −2000 or −200.
- Cases:
  - (7, 10⁻¹⁷) mirror at atol 10⁻¹⁴ and 10⁻²²;
  - C4's (4.9, 10⁻¹¹) at atol 10⁻¹⁴ (control), 10⁻¹⁵ and 10⁻¹⁶;
  - (4.9, 10⁻¹²) at 10⁻¹⁴ and 10⁻¹⁷; (4.9, 10⁻¹⁰) at 10⁻¹⁵; (3, 10⁻¹²) at 10⁻¹⁷; (3, 10⁻¹¹) at 10⁻¹⁶;
  - C4M (7, 10⁻¹³) at 10⁻¹⁸ and (10, 2.5·10⁻¹³) at 2.5·10⁻¹⁸;
  - C4W (near tie) at 10⁻¹⁵.
- All at rtol 10⁻¹².

**"meas" and the slaving** (`capture.py`, `decomp.py`, `farwin.py`, `tangent_dbg.py`, `meas_noext.py`)
- I ran `rp_job` with its internals captured.
- Along the WIN_FAR run I compared:
  - the measured dc_b/ds against −ρ′D;
  - the position part v_b from the right-hand side;
  - the merged block's own force h(c_b);
  - the cross term;
  - the full Jacobian's eigenvalues and those of its c-block against the slaved law;
  - "meas" with the measured velocity used only where resolved.

**The continuum** (`merged.py`, `slaved_merged.py`, `inject.py`)
- *Merged base.* A block with both kernels at c_b and its own coefficient flow (LSODA, rtol 10⁻¹³). It starts from the run's state at s₁ = 400 or 600 on WIN_FAR, with the centroid's D-induced offset y = (∂h/∂c_b)⁻¹h(c_b) removed.
- *Continuation about it.* δċ = J_cc(z_m(s))δc, where J_cc is the 4×4 c-block of `critical_jac_nd`. The α couplings vanish on a merged block at rest (⟨φ∇φ⟩ = 0, h ≈ 0). It is integrated log-scaled with midpoint matrix exponentials, 4 substeps per output interval (16 for a convergence check). Initial vector: δc_e = y − (w/W)D, δc_r = y + (w_e/W)D, or the same with y = 0.
- *Injection.* The continuation is carried to |D| = 10⁻⁹σ (also 10⁻⁷ and 10⁻¹¹σ) and injected into the full float64 race (PairSums, rtol 10⁻¹³), which runs on to the departure.
- *Comparison.* `slaved_merged.py` reruns the committed `relative_law` on the run's base and on the merged base.

**Loops with no near edge** (`creep_test.py`). C4's loop at (4.9, 10⁻⁸) on the window with the top at 1.5, at rtol 10⁻¹² and 10⁻¹³, and on the mirror.

**Mutants** (`mutants.py`, `mutdriver.py`, `survivor_check.py`)
- 34 single-point mutants of the new or changed code.
- Each mutant was run through `--from-raw` (with an md diff) and through `pytest -x`.
- Each surviving job-function mutant was rerun on five RP rows: Part C's x = 4.9 and 6.04; WIN_FAR x = 3 and 5.48; top 1.5, x = 3.

**Mathematics.** By hand: the identity, the mobility matrix, κ, the pin and its direction, the quasi-static elimination of the centroid.

## Verdicts per part

| claim (D54.md line) | verdict | evidence |
|---|---|---|
| (iv) l. 65–74: sides; −1, +1, capture tie; the formed far-side −1; 4.9·10⁻¹⁷ < ε < 2.9·10⁻¹² at x = 7 | holds | C4M rows bitwise; interval recomputed 4.90·10⁻¹⁷ / 2.88·10⁻¹²; boundary loops at s₂/s_{d,e} = 0.87, 1.15 |
| (iv) l. 67 "a formed runner-up's the same in the continuum" | holds for the ε-free race in exact arithmetic only (D3bm-6) | my exact continuation: same sides at all gaps; a float64 loop with no near edge (top 1.5, ε = 10⁻⁸) merges the pair to 0 / 1.9·10⁻¹⁶σ and goes toward on both targets |
| (v) l. 92–95: collapsed for x ≲ 8.6, held for x ≲ 5.07 | holds | R's table |
| (v) l. 119–125: identity, slaved law, κ | holds as algebra; the slaved law is approximate (D3bm-5) | re-derived; mobility matrix with det = 1 |
| (v) l. 127 "κ (here 0.51–0.53)" | inaccurate | raw: κ = 0.610 (x = 3, Part C), 0.620 (WIN_FAR) |
| (v) l. 128–130: a tail's side is the layer's | holds | R; C4M |
| (v) l. 130–133: in the continuum the side is the layer's | holds | exact continuation: e-shares +0.98, +0.80, +0.72, +0.66, +0.58 (C2's target), negatives on the mirror, no sign change; robust to s₁ = 400/600, y, 4 or 16 substeps |
| (v) l. 131 "the only forcing is the centroid's relaxation" | true in the two-kernel race; but negligible, and not what "meas" measures (D3bm-2b) | the merged block's own D-creation sits at float64's floor (10⁻¹⁷–2·10⁻¹⁵σ per unit s) from s = 400 |
| (v) l. 132–133: offsets 10^{−19} (5.48), 10^{−25}–10^{−28} (4.9), 10^{−45}–10^{−104} (3) | partly wrong (D3bm-2) | exact: 10^{−19.0}, 10^{−27.6}, 10^{−98.7}σ |
| (v) l. 133–136: pin D*, κ cancels, behind the drift, away from the nearest edge | holds | derivation; RP's top-edge and right-edge rows |
| (v) l. 137–138: a formed runner-up's side belongs to the sampled problem | holds | C4W (bitwise; robust to atol 10⁻¹⁵); RP |
| (v) l. 138 "s_d … the hold where there is one (x ≲ 5.07: pin-independent)" | wrong (D3bm-1) | continuum s_d +3.8% (4.5), +12.2% (4.9), +23.4% (5.07) |
| (v) l. 139 "x ≈ 5.5–6: from Part C's pin earlier" | holds qualitatively | +38% at 5.48, +2.7% at 6.04 |
| proof l. 192–199 | algebra holds; the slaving condition is unstated | see D3bm-5 |
| Status l. 211–213 "the relative mode's law exact to first order …; slaved rate checked …; the continuum's side by its continuation" | identity exact; rate check partial; side confirmed | D3bm-5 |
| New l. 239–241 | overstated (D3bm-8); offsets wrong (D3bm-2) | |
| C4 l. 310–314: (4.9, 10⁻¹¹) "the one outcome on which the two tolerances differ", closest 7.0·10⁻¹⁵σ | an atol artefact (D3bm-3) | captured toward at atol 10⁻¹⁵ and 10⁻¹⁶ |
| C4W l. 314–317 | holds | bitwise; 3.95·10⁵ × 1681/1887 = 3.52·10⁵ |
| C4M l. 317–326 | holds | bitwise; atol diagnosis confirmed |
| R l. 337–339 "73 of 74 … the one disagreement … at rtol 10⁻¹²" | becomes 74 of 74 at atol ≤ 10⁻¹⁵ | D3bm-3 |
| RP l. 348–369: pins, law ratios, Jacobian at 4.9/7, long double, top edges, right edge, far window's float64 sides | hold | bitwise reruns; the Jacobian check is selective (D3bm-5) |
| RP l. 370–377: continuation least, "meas", climbs, "later by 1.7·10⁴" | wrong numbers, sides right (D3bm-2) | see D3bm-2 |

## Status of the round-62 issues

| D3bi | status | evidence |
|---|---|---|
| -1 (major: continuum's side) | partly fixed | The sides are now stated and I confirm them independently. The numbers are off (10^{−98.7} not 10^{−103.5}; "meas" is an artefact; climbs 7.0–7.5·10⁴). D3bi-1's s_d consequences (+12% at 4.9, +38% at 5.48), which I confirm at +12.2% and +38%, were rejected as "pin-independent" (D3bm-1). |
| -2 (major: rate) | fixed | κ is in the law and the proof. The Jacobian agrees to 0.4–1.6% at x ≥ 4.9 on Part C's window; residue at x = 3 (D3bm-5). |
| -3 (major: s_d) | rescaling fixed; attribution not | N₀/N is validated: x = 7 gives 9.566·10⁴ on both windows. "Pin-independent" is false (D3bm-1). |
| -4 | fixed | C4W and RP right-edge rows, bitwise |
| -5 | fixed | decider, formed far-side −1 and the interval all in (iv); qualifier still needed (D3bm-6) |
| -6 | fixed | collapsed x ≲ 8.6, held x ≲ 5.07 |
| -7 | fixed | quantisation printed; "within 2 steps" is 2.3 (D3bm-9) |
| -8 | fixed | 2.346·10⁷; Λ₀; the garbled ending; s_θ to two digits |
| -9 | fixed, with residue | early grid, boundary loops and atol 10⁻⁵ε below 10⁻¹⁴ are all correct. The 10⁻¹⁴ floor kept at ε ≥ 10⁻¹⁴ is unsafe at 10⁻¹¹ (D3bm-3), and its test is vacuous (D3bm-4). |
| -10 | partly fixed | 17 of 34 new mutants survive (D3bm-4) |
| -11 | fixed | WFR contrast and clustered spikes cited; residue in D3bm-8 |
| -12 | fixed | no `e_t2`/`w_t2` in C4/C4b |

## Issues

### D3bm-1 (major): a held runner-up's s_d shifts one for one with the end of the climb

**Where.** D54.md:138 says "its s_d to the hold where there is one (x ≲ 5.07: pin-independent)".

**Method.** Continuum s_d from the exact continuation about the merged block on WIN_FAR, injected into the full race at |D| = 10⁻⁹σ. Injecting at 10⁻⁷ or 10⁻¹¹σ gives the same s_d to 5 digits.

**Validation.** At the two gaps where float64 itself resolves the offset, the injected s_d matches the plain run:
- x = 6.04: 8.665·10⁴ against 8.650·10⁴;
- x = 7: 9.5655·10⁴ against 9.5659·10⁴.

**Results** (Part C's time):

| x | s_d on Part C's window | s_d in the continuum | change | climb's end on Part C's window | climb's end in the continuum |
|---|---|---|---|---|---|
| 4.5 | 1.210·10⁶ | 1.257·10⁶ | +3.8% | 2.63·10⁴ | 7.32·10⁴ |
| 4.9 | 3.191·10⁵ | 3.580·10⁵ | +12.2% | 3.46·10⁴ | 7.40·10⁴ |
| 5.07 | 1.497·10⁵ | 1.847·10⁵ | +23.4% | 3.86·10⁴ | 7.40·10⁴ |
| 5.48 | 6.10·10⁴ | 8.43·10⁴ | +38% | 5.13·10⁴ | 7.49·10⁴ |
| 6.04 | 8.44·10⁴ | 8.67·10⁴ | +2.7% | 7.27·10⁴ | 7.57·10⁴ |

**The shifts equal the climb's shifts.**

| x | s_d shift | climb shift |
|---|---|---|
| 4.5 | 4.64·10⁴ | 4.69·10⁴ |
| 4.9 | 3.90·10⁴ | 3.94·10⁴ |
| 5.07 | 3.50·10⁴ | 3.55·10⁴ |
| 5.48 | 2.33·10⁴ | 2.36·10⁴ |

**The same holds inside one window.** On WIN_FAR at x = 4.9 the float64 run climbs from the rounding (ending at 4.01·10⁴) and departs at 3.246·10⁵. Shift against the continuum: 3.34·10⁴ against 3.39·10⁴.

**The draft's own data show it.** At top 0.95 (a pin 30 times larger), x = 4.9's s_d moves −1.3% beyond the rescaling. That is −4.1·10³, which equals ln 30/0.000829.

**Mechanism.** Coinciding kernels receive the same correlation ((iii)(b)), so their split is frozen while the pair is merged. The hold's decay starts only once the runner-up reaches its rest.

**Consequence.** C3's capture at x = 5.07 (ε = 3.5·10⁻⁹, s₂ = 1.639·10⁵) belongs to Part C's window. In the continuum s_d = 1.847·10⁵ > s₂, so the runner-up stays bound.

**Proposed fix.**
- (v) l. 137–139: "s_d = (the climb's end) + (the hold after it). The latter is pin-independent, the former is not. From Part C's pin the climb ends at 1.0·10⁴, 3.5·10⁴, 3.9·10⁴, 5.1·10⁴ (x = 3, 4.9, 5.07, 5.48); from the continuum's offset at 7.0–7.5·10⁴. So the continuum's s_d is 3.58·10⁵, 1.85·10⁵, 8.4·10⁴."
- R and C3: mark the x = 5.07 capture as the window's.
- RP: add the injected continuum s_d. It costs about 5 s per gap.

### D3bm-2 (major): RP's continuation below the floor: the sides are right, the numbers are artefacts

**(a) The climb times.**
- **The cause.** `relative_law` evaluates H, w, w_e and ρ′ on the run's base. On WIN_FAR at x = 3, 4.9 and 5.48 the float64 pair leaves the merged state from the rounding, reaching |D_run| = 0.01σ at 1.11·10⁴, 4.01·10⁴ and 6.15·10⁴ (Part C's time). That is long before the continuation's own climb.
- **On a merged base** the same law gives climbs of 7.32·10⁴, 7.49·10⁴ and 7.57·10⁴, against the committed 2.99·10⁴, 4.58·10⁴ and 6.86·10⁴.
- **The exact continuation** gives 7.01·10⁴, 7.40·10⁴ and 7.49·10⁴.
- **A mutant confirms it.** Mutant rl9 evaluates the block at the entrant instead of the midpoint, which is equivalent to first order. It moves the committed climb at x = 3 from 1.02·10⁵ to 1.73·10⁵ (far-window s).
- **Consequence.** "Later by 1.7·10⁴" at x = 5.48 (l. 374–375) should be 2.4·10⁴. It is 3.9·10⁴ at x = 4.9 and 6.0·10⁴ at x = 3.

**(b) "meas" is not the centroid's relaxation.**
- **What is measured.** On WIN_FAR at x = 3 the measured dc_b/ds has cosine −1.000 with D for s = 350–900. The merged block's own force is parallel to −D to within 10⁻³.
- **What it is.** It is the D-induced offset of the slaved centroid: 2.2·10⁻⁹ predicted against 2.86·10⁻⁹ measured at s = 400.7.
- **How it is extrapolated.** `centroid_velocity` extrapolates it from s = 713.6 at a fixed decay of 0.0364/s. Meanwhile D's own contraction accelerates to 0.075–0.15/s.
- **Without the extrapolation** "meas" gives 10^{−103.1}, 10^{−28.0} and 10^{−19.2}σ at x = 3, 4.9 and 5.48, against the committed 10^{−45.0}, 10^{−25.4} and 10^{−19.0}.
- **It also depends on the code's tuning.** Changing the resolution threshold (cv1) gives 10^{−53.7}; changing the fit window (cv5) gives 10^{−42.4}.

**(c) The least |D| at x = 3.**
- **Why the slaved law is wrong there.** On WIN_FAR at s = 400.7 the full Jacobian's relative-mode eigenvalue is −0.0203, equal to the run's rate (−0.0202). The slaved law gives −0.0224. The centroid modes are at −0.169, so the separation is only about 8.
- **Consequence for D.** The law falls to 0.444 of the run at s = 700 and to 0.191 at s = 900.
- **The exact continuation's least |D|:**

| x | least \|D\| (σ) | share along e |
|---|---|---|
| 3 | 10^{−98.7} | +0.98 |
| 4.9 | 10^{−27.6} | +0.80 |
| 5.48 | 10^{−19.0} | +0.72 |
| 6.04 | 10^{−13.4} | +0.66 |
| 7 | 10^{−7.7} | +0.58 |

  These are for C2's target; the mirror gives the negatives. The continuation matches the run to 0.000 dex while the run resolves D. Using 16 substeps changes nothing.

**Proposed fix.**
- RP: continue about a merged block. Alternatively, stop using the run's base once |D_run| > 10⁻⁶σ.
- Drop "meas", or replace it with the merged block's own D-creation.
- Requote l. 132 (x = 3: ≈10^{−99}σ; x = 4.9: 10^{−27.6}σ), l. 241 (10⁻¹⁹ … 10⁻⁹⁹σ) and l. 370–375.
- Test: the climb on WIN_FAR at x = 3 within 5% of 7.0·10⁴.

### D3bm-3 (minor): C4's (4.9, 10⁻¹¹) is an absolute-tolerance artefact

**Evidence.** At rtol 10⁻¹²:

| atol | side | outcome | closest approach (σ) | w_max | departure s |
|---|---|---|---|---|---|
| 10⁻¹⁴ (committed; my control reproduces it) | far | bound, [2, 1, 1] | 1.54·10⁻¹⁴ | 0.02322 | none |
| 10⁻¹⁵ | toward | captured | 1.52·10⁻¹⁴ | 0.02153 | 3.507·10⁵ |
| 10⁻¹⁶ | toward | captured | 1.54·10⁻¹⁴ | 0.02153 | 3.507·10⁵ |

The tighter tolerances give C4b's outcome.
- **The mechanism differs from 10⁻¹⁷.** At 10⁻¹¹ the entry is not early. Instead, atol = 10⁻³ε sets the global step on the 287 dormant kernels' α, and the step error then exceeds the 1.5·10⁻¹⁴σ pin.
- **The draft's tolerance argument misses this.** C4 and R put the difference down to rtol ("the one outcome on which the two tolerances differ"), and `c4_atol` keeps 10⁻¹⁴ "so the earlier rows are unchanged".
- **Unchanged at atol 10⁻⁵ε.** Outcome, side and w_max stay the same (w_max to 4 digits) for:
  - (4.9, 10⁻¹²), (4.9, 10⁻¹⁰), (3, 10⁻¹¹) and (3, 10⁻¹²);
  - C4M (7, 10⁻¹³) and (10, 2.5·10⁻¹³);
  - C4W.

**Proposed fix.** Set `c4_atol = 10⁻⁵ε` for every ε. Rerun C4/C4b at ε ≤ 10⁻¹¹ (about 40 s each). Change C4, R and (v) to read 74 of 74, and drop the 7.0·10⁻¹⁵σ closest approach.

### D3bm-4 (minor): tests

**The c4_atol test is vacuous.** Its second assertion is `np.isclose(c4_atol(1e-17), 1e-22, rtol=1e-12)`. numpy's default atol is 10⁻⁸, so `np.isclose(1e-14, 1e-22, rtol=1e-12)` is True.
- As a result at2 (atol 10⁻⁴ε) survives.
- So does at3 (atol 10⁻¹⁴ always), which restores the round-62 bug.

**Mutation results.** 17 of 34 mutants killed: all 7 md-changing ones by the golden test, and 10 job-function ones by the new RP tests.

**The 17 survivors:**

| mutant | effect on the reruns |
|---|---|
| rl2, rl3 (ρ's sign; ρ′ dropped) | ≤ 2% on the leasts; near-equivalent |
| rl9 (block at the entrant, not the midpoint) | climb at x = 3: 1.02·10⁵ → 1.73·10⁵ |
| rl10 (κ = W²/(2w_e²)) | least at far x = 3: 10^{−103.5} → 10^{−104.8} |
| cv1, cv4, cv5 (threshold; swapped weights; fit window) | "meas" at far x = 3: 10^{−53.7}, 10^{−48.1}, 10^{−42.4} |
| rp1 | the top-1.5 x = 3 "meas" side flips |
| rp3 | s_climb 2.92·10⁷ instead of None |
| rp7 (s₀ = 600) | Part C x = 4.9 "none" least: 10^{−28.0} → 10^{−20.3} |
| rp8 (growth window) | growth rates −2% to −25%; far 5.48: −2.2·10⁻³ against −10⁻⁵ |
| at2, at3 | change the low-mass C4M loops (at3 restores the round-62 bug) |
| at1, cw1, rp4, rp5 | none on the committed data (equivalent) |

**Proposed fix.**
- Use exact equality, or `atol=0`, in the c4_atol test.
- A golden test of `rp_job` at far x = 3 (the leasts to 1 dex, s_climb to 1%).
- A `centroid_velocity` test with w ≠ w_e.
- A κ test at x = 3 (κ = 0.610).
- A test of the growth window.

### D3bm-5 (minor): the rate check is printed only where it is best, κ's range is understated, and the slaving condition is not stated

**κ.** D54.md:127 says κ is "here 0.51–0.53". The raw has κ = 0.610 at x = 3 on Part C's window and 0.620 on WIN_FAR.

**The rate check.** `rp_lines` prints the Jacobian only for x ∈ {4.9, 7} on Part C's window, where the slaved rate is off by −1.6% and −0.4%. The raw's other rows at s ≈ 450:

| row | Jacobian against the slaved law |
|---|---|
| x = 3, Part C's window | −5.0% |
| x = 3, top 0.95 | −4.7% |
| x = 3, top 1.5 | −6.6% |
| x = 3, right edge | −6.0% |
| x = 3, WIN_FAR | −9.2% |
| x = 4.9, WIN_FAR | −3.5% |
| x = 5.48, WIN_FAR | −2.6% |

The slaved rate errs by about |κA|/λ_b.

**The text.** The proof (l. 197–198) and the Status (l. 211–212, "exact to first order") conflate the exact identity with the slaved law.

**Proposed fix.** Print all the Jacobian rows. Give κ as 0.51–0.62. State the slaving condition and its error in (v), the proof and the Status.

### D3bm-6 (minor): "in the continuum" means the ε-free inner problem in exact arithmetic

**Where.** D54.md:67, l. 130–133 and l. 241.

**Why it matters.** The offsets (10⁻¹⁹ … 10⁻⁹⁹σ) lie below float64's floor.
- C4's loop at (4.9, 10⁻⁸) on the top-1.5 window (nearest edge 6.25σ) merges the pair to 0σ (rtol 10⁻¹² and 10⁻¹³) and to 1.9·10⁻¹⁶σ (mirror).
- It goes toward on both targets, although the mirror's continuum side is far.
- Nothing shows that a finite-ε loop's own perturbations stay below these offsets.

**Proposed fix.** Write "in the ε-free inner problem in exact arithmetic (no float64 loop realises it below x ≈ 6)".

### D3bm-7 (minor): RP's table presentation

**Unforced rows on pinned windows.** On Part C's and the right-edge windows, D(400) is already pinned, so the "unforced" continuation inherits the edge's side. For example, on the right edge at x = 3 it is far on C2's target (10^{−90.9}) and toward on the mirror: the opposite of the continuum. The table labels every row "unforced".

**"Exact zero" (top 1.5, x = 3).** This label hides three things:
- "rest" climbs to 0.01σ at s = 1.60·10⁴ (far);
- "meas" climbs at 2.03·10⁴ (toward);
- "none" is NaN at s = 4.5·10⁶, because `relative_law`'s eigen-step has no guard.

**The stop rule.** It is sound (rl5 and rl6 are killed). But if |D| never falls below 0.01σ it would report the run's end as s_climb.

**Proposed fix.** Label which rows are the continuum and which start pinned. Print NaN as NaN. Guard the eigen-step. Report None when the rule never arms.

### D3bm-8 (minor): literature, from memory

**Missing classical relatives.**
- Two coinciding units under plain gradient descent, with a relative-mode rate ∝ (product of the output weights) × the residual's curvature, a drift along the singular line, and plateaus. This is the overlap-singularity literature:
  - Fukumizu–Amari, Neural Networks 13 (2000);
  - Amari–Park–Ozeki, Neural Comput. 18 (2006);
  - Wei–Zhang–Cousseau–Ozeki–Amari, Neural Comput. 20 (2008);
  - Wei–Amari, Neural Networks 21 (2008), for Gaussian RBF networks.
- An escape time ∝ log(1/asymmetry) is Saad–Solla's plateau (Phys. Rev. E 52, 1995). It is D3bm-1's mechanism.
- Related: duplicated-neuron critical sets (Şimşek et al., ICML 2021; Zhang et al., the embedding principle, NeurIPS 2021).
- Position mobility ∝ coefficient is ordinary gradient descent in the (a, θ) parameterisation (balancedness: Du–Hu–Lee, NeurIPS 2018).

**"New" is overstated.** The window-edge pin is a property of the sampled setup, not a dynamical law. It belongs in a caveat.

**Proposed fix.** Cite these after verifying them. Recast l. 239–241.

### D3bm-9 (cosmetic)

- l. 353 says "x = 3 within 2 steps". It is 2.3 steps: the law gives −3.84·10⁻¹⁵ and the run −4.26·10⁻¹⁵σ, with a step of 1.85·10⁻¹⁶σ.
- l. 130 says "contracts through tens of e-folds". At x = 3 it is about 225.
- C4's "closest approach 7.0·10⁻¹⁵σ" at (4.9, 10⁻¹¹) is the artefact's value (D3bm-3).

## Sources

- **No web source was consulted this round.** The session's web-search budget was exhausted at my first query (one attempt, unanswered), and I did not work around it.
- **The WFR and conic-particle points** rest on round 62's search results as recorded in `review/REVIEW_round62_D3bi.md`. I did not re-fetch them. Among them:
  - https://projecteuclid.org/journals/annals-of-statistics/volume-52/issue-4/Learning-Gaussian-mixtures-using-the-WassersteinFisherRao-gradient-flow/10.1214/24-AOS2416.full
  - https://arxiv.org/pdf/1907.10300 (a search summary only)
- **From memory, not verified this round:**
  - Fukumizu–Amari (Neural Networks 13:317–327, 2000);
  - Amari–Park–Ozeki (Neural Comput. 18:1007–1065, 2006);
  - Wei et al. (Neural Comput. 20:813–843, 2008);
  - Wei–Amari (Neural Networks 21:989–1005, 2008);
  - Saad–Solla (Phys. Rev. E 52:4225, 1995);
  - Şimşek et al. (ICML 2021);
  - Zhang–Zhang–Luo–Xu (NeurIPS 2021);
  - Du–Hu–Lee (NeurIPS 2018);
  - Denoyelle–Duval–Peyré (JFAA 23, 2017);
  - Batenkov–Goldman–Yomdin (Inf. Inference 10, 2021).

## Lesson

A continuation below the floor is only as good as the base it is linearised about. RP's law has the right form and was validated at the pin. It was then continued over the run's own trajectory after that trajectory had left the merged state. Its "measured relaxation" was the pair's own feedback, frozen and extrapolated. Together these moved the continuum's numbers by up to 54 decades and by up to 4·10⁴ in s, while the sides happened to survive. Two things would have caught it:
- Linearise the exact Jacobian about a merged base, and inject the result back into the full race. That gives the least |D|, the side, the climb and s_d with a validation built in. It also shows that a hold begins only when the climb ends, so s_d inherits the pin one for one.
- Remember that a tolerance "below the masses" is not enough when the outcome is decided at 10⁻¹⁴σ. The absolute tolerance sets the step for every component, including the dormant ones. A test of such a tolerance must not use `np.isclose`'s default atol.
