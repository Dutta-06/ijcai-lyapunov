# Referee report D3bp, round 69: Proposition D‴.60 as revised after round 63

## Scope

- **Commit reviewed:** e65d6a0, exported with `git archive` to W/repo. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref69/D3bp_work.
  - HEAD has since moved to 42aa6a3. `git diff e65d6a0 HEAD` is empty over every file below.
- **Files read:**
  - review/drafts/D60.md (142 lines);
  - experiments/exp_noise_resonance.py (1222 lines) and its imports: exp_closed_transit.py (LAW, constants), exp_dead_zone.py (`fitted_state`, `noise_of`), exp_cohesion_leak.py (`load`/`save`), src/core/closed_loop.py and closed_loop.c, src/core/mirror_free.py (via `law`);
  - results/chaos/noise_resonance.md (345 lines) and noise_resonance_raw.npz;
  - tests/test_noise_resonance.py (32 tests);
  - review/REVIEW_round63_D3bj.md;
  - the "Round 63" entry of review/REVIEW_LOG.md (lines 2629–2670).
- **I modified nothing in the repository.** All work was in W:
  - W/repo: the export;
  - W/out: the raw copy and my outputs;
  - W/scripts: my scripts;
  - W/mut: the mutants;
  - W/logs: run logs.
- **Run settings:** `nice -n 10`, OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1, LCK_BUILD=W/build, pytest `-p no:cacheprovider`.
- **Load:**
  - At most two heavy processes of mine ran at a time. To run short checks I paused one mutation run with SIGSTOP and resumed it afterwards.
  - No process of mine is left.
- **Web search was unavailable.** The session's shared search budget was already used up at my first query. Sources are marked [memory] or [D3bj snippet, not re-checked].
- Every number is marked (md), (raw) or (mine).

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (md5 f39245037eed3c07fe604f89a0f8797a).
  - The 32 tests pass in 121 s.
  - Cells rerun with the experiment's own functions equal the raw in every field (largest relative difference 0):
    - `f_job`, paradigm and WFR at m = 12, including the Floquet sums and saturations;
    - `b3_job` (paradigm, 13.05) and `b5_job` (paradigm, 13.15);
    - `c3_job` at s = 0.02 for seed 10 (leave 114) and for the draw's pair (leave 218);
    - all nine `a2_rows`;
    - one C cell and one C2 band cell.
- **Every number in D60.md is printed by the md.** The only exception is 146.5 = 29.3/0.2, which is arithmetic. Some ranges are pooled over kernels (D3bp-8).
- **(i) re-derives.**
  - The oscillator law and √K_i = 15.758, 10.505, 13.132 check, as does the bound 100/20.69 = 4.83 (mine).
  - At m = 12 the quasi-static peak sits about 10% below √K_i, and the cause is the sweep, broadened by k_e (D3bp-5).
  - The linearised loop's broad maxima are at 0.84–0.90·√K_i (mine).
- **F's monodromy is right, checked independently.**
  - I linearised the C loop itself by central differences over one period: 43 states, rtol 10⁻¹².
  - The formed block has the multipliers 0.993999, 0.992446, 0.992370, 0.202106, 0.079050, 0, 0 (mine).
  - `monodromy` equals that linearisation to 2·10⁻⁸.
  - The slow eigenvectors are the coefficients w2, w1, w3, with centre components −0.26, −0.07 and +0.09; cond(V) = 1.30.
- **(ii)'s Floquet algebra is right.**
  - The recursion, the geometric sum, the bound, the alternation and the offset all re-derive.
  - The proof's period-mean form ⟨z⟩_j = Σ_k(Pv_k)a_kj + Re(e^{iΩjT}H) reproduces `floquet_means` to 5·10⁻¹¹ (mine).
  - The single tones' even parts approach saturation with time constants of 131–140, 166–167 and 131–132 periods (raw), which are F's three τ.
- **At O(s²) the averaged loop is the constant part's response.**
  - At 13.15, `envelope_means` over B5's constant part is:
    - 0.998–1.005 at s = 0.00125 (mine);
    - 0.993–1.022 at s = 0.0025 (mine);
    - 0.964–1.085 at s = 0.005 (raw).
  - The gap therefore shrinks like s².
  - 13.1's averaging failure is its commensurate 2ω part (B5: +69.84 against −1.02 and −1.584).
  - The 16 tones' averaged-loop error is also O(s²): 0.737, 0.989, 1.233 at s = 0.00125 (mine) against 0.779, 0.998, 1.269 at s = 0.005 (md).
- **C3's chronology holds.**
  - The first kernel's centre stays within ~0.05 of its feature until one or two periods before the leave, and its coefficient does not change through the jump (raw).
  - The rise at the leave is the loop's even part at that level: in the 8 runs whose −n twin had not yet left, the ±n half-sum at s = 0.02 is 0.971–1.021 of the rise (mine).
- **A memoryless escape from the start is rejected.** Under an exponential fitted to C3's 11 kernel-1 leaves at s = 0.02, no leave before period 114 has probability 6·10⁻⁴; with s = 0.03 pooled it is 3·10⁻⁷ (mine).
- **Single tones never jump.** This holds also at s = 0.12 and 0.16 (mine: 8 runs, largest distance 0.34–0.51).
- **The jump-level trends are real.** The jump level falls with m, rises with k_e and is ~16–21× higher under WFR. A basin edge explains all three (D3bp-2).
- **Round 63:** D3bj-1, -2, -4, -7, -9, -10, -13 and -14 are fixed, and D3bj-3 is fixed at O(s²).

**What is wrong (majors first)**

1. **(iii)'s "by about (ii)'s even part at its saturation" is not a like-for-like comparison, and done like-for-like it fails (D3bp-1).**
   - Measured on its own noise draw at its own leave period, s²⟨w₂⟩ is +0.052 … +0.107, at 64–101% of its period-300 value.
   - The rise at the leave is 1.02–1.33 times that.
   - The printed agreement comes from two errors that cancel.
2. **(iii) has no criterion for the jump, and its "drift, then a jump" cannot explain single tones (D3bp-2).**
   - In the draft's own raw, single tones raise the first coefficient by +0.47 … +0.96 and keep the kernel at its feature. The bands jump after +0.06 … +0.29.
   - The missing piece is the fit's basin edge. Without noise the loop ejects kernel 1 within two periods once its coefficient exceeds the feature by +0.075, and keeps it at +0.070 (mine).
   - That edge moves with m, k_e and the metric. "The drift reaches the edge" predicts the bands' first-jump levels within the level grid, and within 1.7× at k_e = 20 and for the 25–75 band.
   - Noise raises the edge: bands raise it slightly, single tones beyond +0.5.
3. **Minors:**
   - the statement's Floquet formula describes the mode coordinate, not the period means (H is missing), and "exponentially stable" holds for the formed block only (D3bp-3);
   - the "resonance factor" is the adaptation loop's sensitivity, which exceeds 1 above √(K/2) for every k_e (D3bp-4);
   - the 10% peak shift is the quasi-static law's, and only at m = 12 (D3bp-5);
   - B5 at s = 0.005 mixes in the averaged loop's O(s⁴) (D3bp-6);
   - the increments do not fall monotonically with ω (D3bp-7);
   - some quotes pool kernels, and the first-passage dichotomy is not tested (D3bp-8);
   - B3's list of exceptions is incomplete (D3bp-9);
   - 18 behaviour-changing mutants survive, and 6 of 12 replicas of round 63's survivors still survive (D3bp-10);
   - bursting and slow-drift instability are missing from the sources (D3bp-11);
   - cosmetic items (D3bp-12).

## What I checked and how

**Reproducibility**
- I copied the raw to W/out and ran `NR_OUT=W/out … exp_noise_resonance.py --from-raw`, then `cmp`.
- `pytest tests/test_noise_resonance.py`.
- `rerun_quick.py` and `rerun_slow.py` reran the jobs listed above and compared them with the raw field by field.

**Numbers**
- `num_check.py` searched every D60.md number in the md.
- `raw_analysis.py` and `c3_traj.py` recomputed the per-kernel ranges, the per-(m, k_e) jump levels, B3's remainders and the trajectories around each leave.

**Mathematics**
- By hand: (i), the proof of (ii), the Floquet algebra and the sensitivity identity.
- `fd_monodromy.py`: the C loop's one-period map at the fitted state, all 43 states, central differences with h = 10⁻⁶, compared with `monodromy`.
- `floquet_check.py`:
  - the slow modes' increments ℓ_kᵀG per unit rms at 13 tones;
  - the period means split into mode sums and Re(e^{iΩjT}H).
- `a2_peaks.py`: the quasi-static peak computed four ways — without the sweep (K and f̂′ frozen at the kernel's K peak), with f̂′ = 0, with both, and with k_e = 1…20; also \|S\| at the bound.
- `var_peaks.py`: the linearised loop (`variational`) on a fine grid near the peaks.

**(ii) numerics**
- `b5_small_s.py`: `envelope_means` at 13.15 for s = 0.00125 and 0.0025.
- `b_small_s.py`: the 16 tones' `envelope_means` at s = 0.00125 and 0.0025.
- `mo3_effect.py`: what the surviving rotation mutants would print.

**(iii) numerics** (C's loop and its 8-sample period means unless stated)
- `seed_w2.py`: for each C3 seed 1–11:
  - `sim_means` under ±s·n at s = 0.00125 over 300 periods, giving each draw's own ⟨w₂⟩;
  - the −n twin at s = 0.02 and at 0.03.
  The +n run reproduces C3's W exactly.
- `single_tones.py`: 13.1 and 15.0 at s = 0.02, 0.04, 0.12, 0.16, seeds 7 and 1.
- `basin.py`, `basin2.py`, `basin3.py`: the noise-free loop from the fit with one coefficient raised (the plant unchanged), scanned against m, k_e and the metric.
- `prerise.py`: the bands and a single tone started with the first coefficient pre-raised.
- `edge_vs_drift.py`: ⟨w₂⟩ for seed 7's draw at m = 3 and 8, and at k_e = 5 and 20.
- `band_w2.py`: ⟨w₂⟩ for the 10–30 and the 25–75 band.
- `pair_removal.py`: each band without its most nearly commensurate pair, and without a control pair.
- `first_passage.py`: exponential-from-zero tests on the leave periods.
- From the raw: C2's end coefficients (`w_end`), and the draw's pair's W in C3.

**Mutants**
- `mutants.py` defines 65 single-point mutants: 53 in the new or changed code and 12 replicas of round 63's survivors.
- `mutdriver.py` applies each in a private copy (experiments/ and tests/ copied, src linked), runs the test file with `-x`, then `--from-raw` and `cmp`.
- `mut_targeted.py` runs the replicas against the only tests that execute the mutated function.

## Verdicts per part

| claim (D60.md line) | verdict | evidence |
|---|---|---|
| (i) first order = the linearisation; y″ + k_eẏ + Ky = K(ṅ − f̂′n), quasi-static (8–13, 50–55) | verified | Re-derived. A: loop over linearisation 0.9899–1.0072; quasi-static over linearisation 0.999–1.082 at ω ≥ 10 (md) |
| (i) √K_i = (Gm)^{1/2}\|w_i\|/(σe^{1/2}); the bound on m (15–17, 55–56) | verified | 15.758, 10.505, 13.132; bound 4.83 (mine) |
| (i) "the response peaks about 10% below √K_i (the sweep and the damping k_e)" (16) | verified for the quasi-static law at m = 12; the cause is the sweep | 0.899, 0.872, 0.886 (mine). With the sweep removed: 0.982–0.992. With no sweep and f̂′ = 0: exactly 1.000. At m = 3, f̂′ dominates (D3bp-5) |
| (i) "reduces the resonance, without removing it: the half-power width is k_e"; 1.50 → 1.23 (17–19) | numbers verified (first kernel); the reason is the K = 0 reference | the factor is the loop's sensitivity, > 1 for ω > √(K/2) whatever k_e (D3bp-4) |
| (ii) even/odd split, finite horizons (21–24, 57–60) | proved; numerics mostly agree | B3 doublings mostly ×16 and ×8, but not for 13.15's or 13.05's first kernel (D3bp-9) |
| (ii) slow multipliers 0.99400, 0.99245, 0.99237 = the coefficient modes, the rest < 0.21 (24–26) | verified independently | finite differences of the C loop (mine) |
| (ii) "the linearised loop is … exponentially stable"; "all \|ρ\| < 1" (24–25, 64) | true of the formed block only | 36 further multipliers within 8·10⁻¹⁰ of 1 (the unformed kernels; mine) |
| (ii) "the period means follow m_{j+1} = ρm_j + e^{iΩjT}δ", the closed form, bounded, alternating, offset with τ (26–33) | right for the mode coordinate (proof 60–64); the period means lack the H term | H on the coefficients is 1.36× the largest mean at 13.1 (mine; D3bp-3) |
| (ii) "increments falling steeply with ω" (34) | true overall, not monotone | D3bp-7 |
| (ii) constant part = the averaged loop's response, an offset reached with τ, "not a drift without bound" (35–37) | verified at O(s²); consistent with the finite-horizon expansion | 13.15: 0.998–1.005 at s = 0.00125 (mine); the saturations above. The bounded offset is what exponential stability gives. It is τ-amplified: at s = 0.02 the band's even part is 1.23–1.33× s²⟨w₂⟩, and 13.1's is 0.55–0.59× at period 100 (mine) |
| (ii) the averaged loop's error = the dropped parts' Floquet sums; 13.1 against 13.05/13.15/29.3 (37–39, 64–67) | holds at O(s²) (argued, now checked); the Σ\|…\| formula is a bound; 29.3 inferred | B5 (md); my small-s runs (D3bp-6) |
| (iii) the coefficient drifts up while the centre stays (41–42) | verified as chronology | C3 trajectories (raw) |
| (iii) "by about (ii)'s even part at its saturation, s²⟨w₂⟩" (42; 140–141) | refuted as stated | D3bp-1 |
| (iii) "then it jumps out of the orbit's range within a period" (42–43) | verified for kernel 1 | 3 of 38 C3 jumps end inside [−1, 1], at −0.50 … −0.54 (raw) |
| (iii) "after about one slow-mode time constant — … not a memoryless first passage" (43–44) | timing verified; the dichotomy is untested and false | median 164 periods against τ = 131–166. Memorylessness from the start is rejected; a delayed first passage fits (mine) |
| (iii) single tones never jump; pairs and bands do (44–46) | "never" verified (also at 0.12 and 0.16); "their centre too" holds only at 0.08; pairs: 2 of 4 classes | C2 (md); mine |
| (iii) the mechanism "a drift on the slow modes, then a jump" (1–2, 41, 87) | incomplete | no criterion; the basin edge (D3bp-2) |
| (iii) the jump level falls smoothly with Gmw (46) | verified (4 seeds, 5 levels); explained by the edge | C2 m-scan 0.08 … 0.03, and 0.015 at m = 12 (md) |
| (iii) the jump level rises with k_e (46–47) | verified (4 seeds) | kernel 1's first jump: 0.0075, 0.015, 0.06, 0.08 (1/4) (raw) |
| (iii) a band above √K_max raises it, "mostly through the roll-off" (47) | level verified; attribution inferred from the first-order jitter | 0.06 against 0.015 (md); the drift's band ratio is 5.38, not 2.92² (mine; D3bp-4) |
| (iii) WFR "an order of magnitude higher" (47–48) | verified | ×16 at k_e = 5 (0.12 against 0.0075), ×21 at k_e = 10 (0.32 against 0.015) (md); WFR's edge is (0.3, 0.4] (mine) |
| Status (69–75) | mostly accurate | "averaged loop's error located in the commensurate 2ω part" holds for 13.1 only (D3bp-6) |
| Classical / New (76–88) | partly | bursting uncited; "a drift … to (ii)'s offset" is unsupported (D3bp-1, -2, -11) |
| Numerics (90–142) | all printed by the md | pooled ranges at 138–140 (D3bp-8); B3 exceptions at 107–108 (D3bp-9) |

## Status of the round-63 issues

| id | status | what remains |
|---|---|---|
| D3bj-1 | fixed | C2 now runs single tones (6 × 0.04–0.08 × 4 seeds), pairs (4 classes × 0.02–0.08 × 4 seeds) and D‴.49's own pair, with jumps reported apart. Residue: "pairs … do" holds for 2 of the 4 classes (D3bp-8) |
| D3bj-2 | fixed | Floquet sums; F's multipliers (independently verified); the commensurate saturation (−0.5230 against the loop's −0.521). Residue: the stated period-mean formula lacks H; "exponentially stable" holds for the formed block (D3bp-3) |
| D3bj-3 | fixed at O(s²) | B5 puts 13.1's failure in its 2ω part, and my small-s run confirms the identity. Residue: B5 at s = 0.005 includes O(s⁴); 29.3 is inferred; the 16 tones' error is not located (D3bp-6) |
| D3bj-4 | fixed | B4 runs to 1200 periods with bands over periods 601–1200; the saturations follow F's τ (raw) |
| D3bj-5 | partly fixed | A2's split, "reduces" and the m-scan under four seeds are in. Residue: the factor is the sensitivity; the ejection attribution rests on the first-order jitter (D3bp-4) |
| D3bj-6 | partly fixed | Classical sources cited; "New" limited to the instance. Residue: bursting and slow-drift instability for (iii); the topic of KRP 1985 (D3bp-11) |
| D3bj-7 | fixed | C2 runs under seeds 7, 1, 2, 3; the draw's pair is one run per level by design |
| D3bj-8 | partly fixed | C3 measures the anatomy. Residue: the size claim fails like-for-like (D3bp-1); the criterion is missing (D3bp-2); "not a memoryless first passage" is untested (D3bp-8) |
| D3bj-9 | fixed | WFR to s = 0.32 |
| D3bj-10 | fixed | leaves and jumps are reported apart in C, C2 and C3 |
| D3bj-11 | mostly fixed | the remainders are printed. Residue: the list of exceptions (D3bp-9) |
| D3bj-12 | partly fixed | 6 of my 12 replicas of round 63's survivors are now killed and 6 still survive; 18 new behaviour-changing survivors (D3bp-10) |
| D3bj-13 | fixed | 0.898, 0.871, 0.887 are printed. Residue: the cause and the m-dependence (D3bp-5) |
| D3bj-14 | fixed | the rectification's sign is dropped; "slow modes" is now plural |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bp-1 | major | **(iii)'s "by about (ii)'s even part at its saturation, s²⟨w₂⟩" (D60.md:42, 140–141; `c3_lines`' last paragraph) is not like-for-like, and like-for-like it is wrong in both qualifiers.**<br>• What is printed: B4's band for seed 7's draw (periods 601–1200, s = 0.00125 scaled by s²) against C3's leaves under seeds 1–11. That is ten other draws, at periods 114–274, at 16–24 times the level.<br>• Like-for-like (each run's own draw at its own leave period; mine): s²⟨w₂⟩ = +0.052 … +0.107 at s = 0.02. That is 0.64–1.01 (median 0.78) of its period-300 value, so the leaves come before saturation.<br>• The rise equals the loop's even part at the level: the ±n half-sum is 0.971–1.021 of it in the 8 seeds whose −n twin had not left (0.97–1.03 at s = 0.03, 5 seeds).<br>• The rise is 1.23–1.33× the own s²⟨w₂⟩ in 7 of those 8 seeds (1.02 in seed 11); 0.88–1.36 at s = 0.03. Across draws it tracks it (correlation 0.82, 11 seeds).<br>• So the printed overlap (+0.094 … +0.116 against +0.064 … +0.143) comes from errors that cancel: the later horizon adds about as much as the missing higher orders take away. At the threshold level the overlap is close to automatic anyway, since the threshold is where the saturated drift first reaches the edge (D3bp-2). | (iii) and C3 (`c3_job`, `c3_lines`): compare each run with its own draw's even part at its own leave period and level (the −n twin costs one run). Say "the even part at the level, ~25% above s²⟨w₂⟩ at s = 0.02, before saturation". Drop "at its saturation" (D60.md:42, 140) and "to (ii)'s offset" (87) |
| D3bp-2 | major | **(iii) states no criterion for the jump, and "a drift on the slow modes, then a jump" (title; D60.md:41–43, 87) cannot explain why single tones, which drift further, never jump. The missing piece is the fit's basin edge in the coefficient, which noise moves.**<br>• Drift without a jump, in the draft's own raw:<br>&nbsp;&nbsp;– all 24 single-tone C2 runs at s = 0.04 rise +0.47 … +0.70 by period 300, with the centre within 0.139, and none leaves;<br>&nbsp;&nbsp;– 21.3 at s = 0.08 rises +0.96 with the centre within 0.23;<br>&nbsp;&nbsp;– D‴.49's pair rises +0.48 and +0.74 at s = 0.04 and 0.06 without leaving, though at s = 0.02 it jumps after +0.12 (C3);<br>&nbsp;&nbsp;– C's seed 1 at s = 0.04 keeps kernel 1 (distance ≤ 0.296) at +0.556.<br>&nbsp;&nbsp;Mine: 13.1 and 15.0 at s = 0.02 rise +0.31 … +0.35, with the centre ≤ 0.045.<br>• The edge (mine; noise-free loop from the fit, one coefficient raised, the plant unchanged):<br>&nbsp;&nbsp;– kernel 1 returns from +0.070 and is ejected within 2 periods from +0.075, at rtol 10⁻¹² too; kernel 3's edge is in (0.06, 0.08];<br>&nbsp;&nbsp;– the edge widens as m falls: (0.13, 0.16] at m = 3, (0.08, 0.10] at m = 8;<br>&nbsp;&nbsp;– it widens with k_e: (0.04, 0.06] at k_e = 5, (0.10, 0.13] at 20, (0.13, 0.16] at 40;<br>&nbsp;&nbsp;– under WFR it is (0.3, 0.4].<br>• "The drift reaches the edge", s* = (edge/⟨w₂⟩₃₀₀)^{1/2} with seed 7's ⟨w₂⟩ (mine), against the observed first jumps (md):<br>&nbsp;&nbsp;– m = 12: (0.016, 0.019] against 0.015 (2/4) and 0.02 (all);<br>&nbsp;&nbsp;– m = 8: (0.028, 0.031] against 0.03;<br>&nbsp;&nbsp;– m = 3: (0.066, 0.074] against 0.08;<br>&nbsp;&nbsp;– k_e = 5: (0.008, 0.010] against 0.0075;<br>&nbsp;&nbsp;– k_e = 20: (0.035, 0.040] against 0.06 (none at 0.04);<br>&nbsp;&nbsp;– the 25–75 band: (0.038, 0.044] against 0.06 (none at 0.04).<br>• Noise moves the edge:<br>&nbsp;&nbsp;– bands leave at +0.064 … +0.143 at s = 0.02, straddling the noise-free 0.072, and at +0.18 … +0.29 at s = 0.03;<br>&nbsp;&nbsp;– a single tone holds the kernel far beyond it: pre-raised by +0.31 at s = 0.02, or by +0.50 at s = 0.04, it never leaves (mine).<br>• Pre-raising a band run's coefficient advances its leave: seed 7 from 133 to 60 with +0.084, seed 3 from 157 to 28 with +0.112, seed 10 from 114 to 14 with +0.15. But not always: seed 10 with +0.064 still leaves at 114. So the drift is the clock and a noise event the trigger (mine).<br>• Removing the most nearly commensurate pair from a band delays the leave no more than removing a median pair (seeds 7, 10, 3, 1; mine), so no single pair triggers it. | State (iii) with the criterion: the drift carries the coefficient to the basin's edge, ≈ +8% at m = 12, k_e = 10. Noise raises that edge, much more for single tones (a vibrational stabilisation: Kapitza, Meerkov, already cited). Print the edge against m, k_e and the metric, test s* = (Δ/⟨w₂⟩)^{1/2}, and only then name the mechanism in the title and in "New" |
| D3bp-3 | minor | **(ii)'s statement of the Floquet sums (D60.md:24–30, 37) describes the mode coordinate, not the period means; its stability claim holds for the formed block only.**<br>• The proof has the period means right (D60.md:60–64). Per mode they are (Pv_k)a_kj + Re(e^{iΩjT}H), which reproduces `floquet_means` to 5·10⁻¹¹ (mine).<br>• "The period means follow m_{j+1} = ρm_j + e^{iΩjT}δ" holds only with m₀ = h ≠ 0. The closed form from rest misses a rotating term, which on the coefficients is 1.36× the largest period mean at 13.1, 0.23× at 21.37 and 0.003× at 21.4 (mine). The stated bound fails at 13.1 without + \|H\|.<br>• "Σ_Ω\|δ_Ω\|\|ρ^j − e^{iΩjT}\|/\|ρ − e^{iΩT}\|" (D60.md:37) sums moduli over parts and modes, without P and H. It is a bound, not "the averaged loop's error".<br>• "Exponentially stable", "all \|ρ\| < 1" (D60.md:24–25, 64) holds for the 7-state formed block. The full linearisation has 36 more multipliers within 8·10⁻¹⁰ of 1 (mine): the unformed kernels, whose centres are neutral (v = 0 at α = 0) and whose α are driven at O(ε). | Write the period means as the mode sums plus Re(e^{iΩjT}H) and add \|H\| to the bound; write "≤" at line 37; say "the formed kernels' linearisation" and note the neutral, O(ε)-coupled rest |
| D3bp-4 | minor | **(i)'s "reduces the resonance, without removing it: the oscillator's half-power width is k_e" (D60.md:17–19): the reason is the K = 0 reference, not k_e. (iii)'s "mostly through the band's roll-off" (47) is inferred from the first-order jitter.**<br>• Per tone, A2's factor is the adaptation loop's sensitivity \|S\| = \|iω(iω + k_e)\|/\|K − ω² + ik_eω\| (f̂′ cancels). From \|S\|² − 1 = K(2ω² − K)/\|K − ω² + ik_eω\|², \|S\| > 1 for every ω > √(K/2) and every k_e. So no m > 0 removes it (Bode's waterbed); k_e only sets its size, \|S(√K)\| = (1 + K/k_e²)^{1/2}.<br>• At the bound the band sees \|S\| = 1.41, 1.39, 1.24, 1.11 at ω = 10, 15, 20, 30 (mine). The 25–75 band at m = 12, wholly above √K_max, still has factors 1.144, 1.056, 1.092 (md). This is the sensitivity's skirt, not a peak.<br>• As a factorisation (closed-loop over open-loop) the decomposition is well defined, and its roll-off part is K-independent (2.24, 2.21, 2.22 at m = 3 and 12; md). But it applies to the jitter. The drift — (iii)'s own mechanism — has no K = 0 counterpart. Its band ratio is 5.38, 3.29, 3.91 (mine), not 2.92², 2.61², 2.79² = 8.5, 6.8, 7.8. | Call the factor the adaptation loop's sensitivity and give \|S\| at the bound. For (iii), compare the drift (or the edge, D3bp-2) between the bands, or drop "mostly through the roll-off" |
| D3bp-5 | minor | **"The response peaks about 10% below √K_i (the sweep and the damping k_e)" (D60.md:16) is the quasi-static law's, at m = 12.**<br>• Mine: without the sweep the peak is at 0.992, 0.982, 0.988·√K_i. With no sweep and f̂′ = 0 it is exactly 1.000: the damping alone does not move a velocity response's peak. With the sweep the shift grows with k_e: 0.983 at k_e = 1, 0.899 at 10, 0.849 at 20.<br>• At m = 3 the plant's slope dominates: 0.738, 0.381, 0.626, against 0.888, 0.884, 0.883 with f̂′ = 0. A2's "peaks" at 4.00 for m ≤ 4 sit on the grid's edge (D60.md:98 "lower fractions at smaller m").<br>• The linearised loop (mine): kernel 1 is flat within 1% over ω = 12–14.5 (max 1.142 at 13.5 = 0.86·√K₁); kernel 2 peaks at 9.5 (0.90); kernel 3 at 11.0 (0.84). The quasi-static law overshoots there by up to 8%. | Say "the quasi-static response at m = 12; the sweep, broadened by k_e". Extend OM_A2 below 4, or mark edge maxima |
| D3bp-6 | minor | **B5's evidence (D60.md:113–116) and the Status's "located in the commensurate 2ω part" (73–74).**<br>• At s = 0.005, B5's averaged loop over the constant part (0.959–1.085; kernel 3 at 1.07–1.085) is mostly the averaged loop's own higher order. That error is larger than the 2ω parts it is meant to separate at the neighbours (0.1% of kernel 3's ⟨w₂⟩).<br>• At 13.15 the ratio is 0.998–1.005 at s = 0.00125 and 0.993–1.022 at 0.0025 (mine). The identity holds at O(s²), but the printed B5 does not show it.<br>• 29.3's "increment is small" is inferred from B3's 1.0003–1.0019; B5 has no 29.3 row.<br>• The 16 tones' averaged-loop error is O(s²) (0.737, 0.989, 1.233 at s = 0.00125; mine) but is not located. | Run B5 at s = 0.00125; add 29.3; split the 16 tones' error (drop each near-commensurate pair's cross average); restrict "located" to 13.1 |
| D3bp-7 | minor | **"Their increments falling steeply with ω" (D60.md:34) and the forcing's frequencies (34–35, 64–65).**<br>• The slow modes' increments per unit rms fall overall, from ≈ 0.02 at ω = 10 to ≈ 10⁻⁴ at 29.4, but not monotonically (mine):<br>&nbsp;&nbsp;– kernel 1's mode: 0.0027 at 15.0, 0.0101 at 18.9, 0.0043 at 21.4;<br>&nbsp;&nbsp;– kernel 2's mode: 4.6·10⁻⁵ at 25.0, 2.1·10⁻³ at 29.3, 8·10⁻⁶ at 29.4.<br>&nbsp;&nbsp;The orbit's two passes per period interfere (fringes ≈ 0.42 rad/s for kernel 2).<br>• "Products … carry 0, 2ω_k, ω_j ± ω_k" holds modulo nω_ref. It also omits the first order's transient on the slow modes (τ = 130–166 periods). B5's φ-split is exact as φ-dependence, so its half-sum includes that transient. | "fall on average, with interference"; "modulo nω_ref, the transient included" |
| D3bp-8 | minor | **C3 quotes and (iii)'s wording (D60.md:43–46, 138–140).**<br>• "Its centre stays … −0.016 … +0.020" pools kernels 1 and 3; kernel 1 alone gives −0.016 … −0.011, always outward. "At 0.03 likewise (1.092–1.322)" pools kernels 1–3; kernel 1 alone gives 1.201–1.322 (raw).<br>• "Not a memoryless first passage" has no test. Memorylessness from the start is rejected (P = 6·10⁻⁴), but a delayed first passage fits (delay 114 + mean 55 periods). The 11 runs are 11 different draws (mine).<br>• "Earlier at higher levels": C's kernel-1 medians are 160, 154, 152, 209 at s = 0.02–0.04, and clearly earlier only from 0.06 (75) (raw).<br>• "Tone pairs … do": only 2 of the 4 classes, plus the draw's pair (md).<br>• "Out of the orbit's range": 3 of 38 C3 jumps and 11 of 141 C2 jumps end inside [−1, 1], mostly kernel 3 at −0.50 … −0.64, the vacated first feature (raw). | Per-kernel ranges; print a leave-time test, or say "none before ~0.7τ"; "some pairs"; "out of its feature" |
| D3bp-9 | minor | **B3's "mostly ×16 and ×8 but not on the 16 tones' first kernel" (D60.md:107–108) omits the resonant neighbours.**<br>• 13.15's first kernel: odd ×7.61 then ×0.092; even ×12.8 then ×3.98; its relative even error is −1.04%, −3.33%, −3.31%.<br>• 13.05's first kernel: odd ×8.93 then ×15.0 (raw).<br>• 21.37's third kernel (×76) and WFR 29.3's first (×42) are at noise level. | List these; note that the expansion is not polynomial near √K₃ = 13.13 by s = 0.005 |
| D3bp-10 | minor | **18 of the 53 single-point mutants of the new code survive both the tests and `--from-raw`, all behaviour-changing. A 19th, CL6, is equivalent on this raw.**<br>• Killed: 15 by the tests only, 16 by `--from-raw` only, 3 by both. Every report-function mutant but CL6 changes the md.<br>• Survivors:<br>&nbsp;&nbsp;– `monodromy` and `floquet_means`: MO3, MO4, FM1–FM3 (Im G, Im H, or the rotation conjugated or dropped). The tests use only 21.4 and 21.3, where e^{iωT} = ±1. MO3 would turn F's 21.37 line from 2.40·10⁻² into 1.65 (mine).<br>&nbsp;&nbsp;– `f_job`: FJ4 (phase dropped), FJ5 (×1/√2), FJ6 (two periods).<br>&nbsp;&nbsp;– `b5_job`: B5b and B5c (the averaged loop at the wrong s, or per s); B5d (c2_shift).<br>&nbsp;&nbsp;– all five `c3_job` mutants (no test calls it).<br>&nbsp;&nbsp;– `leave_stats`' window (LS1).<br>&nbsp;&nbsp;– `a2_rows`' K = 0 reference without f̂′ (AR1).<br>• Round 63's replicas: 6 of 12 are now killed (b4_job: π/2, /s, WFR; c_noise's seed; c_run's metric and k_e). Six survive: b4_job's centre offset (below the test's 10% tolerance) and its odd_c_per column; c_noise's band seed; part_C with one seed; C2's WFR case under the paradigm; b3_job's even_c offset (T3). | tests/test_noise_resonance.py:<br>• Floquet means at 21.37 over 3 periods against `second_order`.<br>• `f_job`'s sums against `floquet_means` with phase and rms, and its 3-period table.<br>• `b5_job`'s w_env against `envelope_means`/s².<br>• `c3_job` with a fake `simulate`.<br>• `leave_stats` with a varying offset.<br>• `a2_rows`' jit0.<br>• odd_c_per's column, the band seed, C2's WFR metric, b3's even_c |
| D3bp-11 | minor | **Sources (D60.md:76–88).**<br>• (iii)'s "drift, then a jump" is the phenomenology of slow-drift instability (Fu–Sastry 1987, cited) and of bursting (Anderson 1985 [memory], uncited); see also parameter drift without persistent excitation (Egardt 1979; Ioannou–Kokotović 1984 [memory]).<br>• Kokotović–Riedle–Praly 1985 is a stability criterion for slow adaptation, not disturbance-driven drift [memory].<br>• "Without removing it" is Bode's sensitivity integral [memory].<br>• With the formed block exponentially stable, the expansion could be uniform in time (Hoppensteadt 1966; Khalil §10.1 [memory]). | Cite these in the Status; keep the phenomenology out of "New" |
| D3bp-12 | cosmetic | • "1.50 at m = 12 to 1.23 at the bound" (18–19) is the first kernel's; the others go 1.248 → 1.095 and 1.370 → 1.152 (md).<br>• C2's caption "> 0.6 (out of the orbit's range)": for kernel 2 (c = 0.07), a distance of 0.6 stays inside [−1, 1].<br>• The round-63 log credits tests to "c3_job's statistics"; they cover `leave_stats` and `c3_lines` only. | Wording |

## Sources

- B. D. Riedle, B. Cyr, P. V. Kokotović, "Disturbance instabilities in an adaptive system", IEEE TAC 29(9):822–824, 1984 — [memory; D3bj snippet].
- P. V. Kokotović, B. D. Riedle, L. Praly, "On a stability criterion for continuous slow adaptation", Systems & Control Letters 6(1):7–14, 1985 — [memory; D3bj snippet https://cas.minesparis.psl.eu/~praly/Telechargement/Journaux/1985-SCL-Kokotovic-Riedle-Praly.pdf].
- L.-C. Fu, S. Sastry, "Slow-drift instability in model reference adaptive systems — an averaging analysis", Int. J. Control 45(2):503–527, 1987 — [D3bj snippet https://scholars.lib.ntu.edu.tw/handle/123456789/332237].
- B. D. O. Anderson, "Adaptive systems, lack of persistency of excitation and bursting phenomena", Automatica 21(3):247–258, 1985 — [memory].
- Anderson, Bitmead, Johnson, Kokotović, Kosut, Mareels, Praly, Riedle, *Stability of Adaptive Systems: Passivity and Averaging Analysis*, MIT Press 1986 — [memory].
- B. Egardt, *Stability of Adaptive Controllers*, Springer LNCIS 20, 1979; P. A. Ioannou, P. V. Kokotović, Automatica 20(5):583–594, 1984 — [memory].
- Rohrs, Valavani, Athans, Stein, IEEE TAC 30(9):881–889, 1985; Sastry–Bodson 1989 — [memory].
- B. Widrow et al., Proc. IEEE 63(12):1692–1716, 1975; J. R. Glover, IEEE Trans. ASSP 25(6):484–491, 1977 — [D3bj snippets].
- H. W. Bode, *Network Analysis and Feedback Amplifier Design*, 1945; J. S. Freudenberg, D. P. Looze, IEEE TAC 30(6):555–565, 1985 — [memory].
- V. A. Yakubovich, V. M. Starzhinskii, *Linear Differential Equations with Periodic Coefficients*, 1975; N. M. Wereley, PhD thesis, MIT 1991 — [memory].
- F. C. Hoppensteadt, Trans. AMS 123:521–535, 1966; H. K. Khalil, *Nonlinear Systems*, 3rd ed., §10.1 — [memory].
- A. I. Neishtadt, Russian Math. Surveys, 2014 — [D3bj snippet https://www.mathnet.ru/eng/rm9603]; LMS non-Wiener effects — [D3bj snippet https://vtechworks.lib.vt.edu/handle/10919/26772].

## Lesson

To test a mechanism written as a sequence — "a drift on the slow modes, then a jump" — look at the cases that show the first step and not the second. The draft's own raw has them: every single tone at s = 0.04, and D‴.49's pair at 0.04–0.06, push the same coefficient further with the centre at its feature and never leave. Asking why leads to what the draft lacks: the fit's basin edge. Without noise the loop ejects an over-amplitude kernel at +8%, and that edge reproduces the jump levels' dependence on m, k_e and the metric. A theory-against-event check must also fix the draw, the time and the level. Comparing one draw's saturated second order with ten other draws' leaves, one to two time constants earlier and at sixteen times the level, matched only because the later horizon's surplus cancelled the missing higher orders.
