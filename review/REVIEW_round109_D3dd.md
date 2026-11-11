# Referee report D3dd, round 109: Proposition D‴.59 as revised after round 106

## Scope

- **Export and work directory.** I exported commit 319e8a5 (branch claude/vigilant-clarke-hf939b) with `git -C … archive 319e8a5 | tar -x -C W/repo`. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref109.
  - W/repo/venv is a symlink to the repository's venv (Python 3.13.18, numpy 2.5.4, scipy 1.18.1, pytest 9.1.1).
  - The venv has no mpmath, so I did the high-precision work with Python's `decimal`.
- **Diffs (5ee87f7 → 319e8a5).**
  - The raw `moving_width_raw.npz` is not in the diff. Every number the revision added is computed by `--from-raw`.
  - experiments/exp_moving_width.py gains `a2_ulps_exact`, `a2_exit`'s `method/rtol/off/signed`, `a2_exits`, `a2_lemma_signs`, the bottleneck `mirror_asym`, `mirror_odd`, the changed `a2_lines` lines (symmetric line, exact-start clause, held line, exit line, lemma edges), `split_summary`'s mirror block and `report`'s C sentence.
  - The md changes at the same places.
  - tests/test_moving_width.py gains two tests (now 54).
  - tools/mutants/round106_D59_own.py is new.
  - review/drafts/D59.md changes at (i), the lemma, the Fenichel clause, (ii)(b) and the A2E numerics. Its README row changes too.
- **Read in full:**
  - D59.md (406 lines);
  - every changed part of exp_moving_width.py, and the code it calls: `a2_state`, `a2_ulps`, `a2_fixed_point`, `a2_lemma`, `a2_lemma_scan`, `a2e_job`, `critical_rhs_nd`/`critical_jac_nd`, `critical_w`, `_alpha_of`, `BumpsCW.corr/synth`, `h_job`, `split_summary`;
  - results/chaos/moving_width.md; tests/test_moving_width.py; tools/mutants/round106_D59_own.py; tools/mutants.py;
  - the Round 106 entry of review/REVIEW_LOG.md; review/REVIEW_round106_D3da.md; the git diff of these files.
- **Settings.**
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` under `timeout`. MW_OUT was W/out (a copy of the raw) or a mutant slot's `out`.
  - pytest ran with `--basetemp` inside W. My copy of the mutant tool, W/own/mutants_tool.py, sets REPO to W/repo and adds `--basetemp` per slot.
  - At most two heavy processes of mine ran at a time.
- **One of my own processes stopped.** A background wait loop of mine used `pgrep -f "rerun.py h"`, a pattern that matches its own command line. I stopped it by its task id (TaskStop) before it mattered. No process I did not start was touched.
- **The repository is untouched.** I used only read-only git (archive, diff, log, status, rev-parse).
  - Meanwhile another session committed ce7cfac (D‴.58) and holds uncommitted D60 work. Neither touches a D59 file: `git diff 319e8a5` on D59.md, the script, the md, the raw, the tests and the spec is empty at the end.
- **Tags.** [md] the committed md or draft; [raw] the committed raw; [own] my computation; [snippet] a web-search result I saw; [memory] not verified this round.

## Summary

**Design.**
1. I regenerated the md from the raw, ran the 54 tests, and reran 10 jobs against the raw.
2. I recomputed the zero of ċ at the two held states with my own decimal code: the full prefactor Gm·w/(Nσ²), w computed from α in decimal, at 30, 50 and 80 digits. I checked the formula against `critical_rhs_nd` at large offsets, and checked linearity, σ sensitivity and f* as stored against f* exact.
3. I reproduced `a2_exits`. Then I reran the exit from c* with the right-hand side's arithmetic changed and the data kept: numpy sum, `math.fsum`, reversed sample order, and long double with a δ = c − c* state. I did this on A's samples and on exactly symmetric ones, plus the symmetric problem from several offsets.
4. I followed the exit's seed along the route: the data's drift d(s) = ċ at c = c* in 50 digits, λ(s) = ∂ċ/∂c, and the linear seed ∫e^{−Λ}d. I compared it with the offset start's seed and with the symmetric problem's summation residue.
5. I located the lemma's sign edges by root finding (Brent), and checked the totals against finite differences of the flow.
6. I tested `mirror_asym` against an independent bottleneck (linear threshold scan with augmenting paths) and against brute force over all permutations and over involutions: 400 adversarial inputs, n ≤ 21, ties, duplicates, negative widths, {a, a, Ra}. For `mirror_odd` I used actual sample pairs, a grid that resolves every kernel, and the formed kernels' function alone.
7. Mutants: the repository's 15, plus 26 of my own, each set with a control.
8. Web searches.

**Verdicts.**
- **Reproducibility is exact.** The md matches byte for byte. 54 tests pass. All 10 rerun jobs match the raw bit for bit.
- **Every number the revision added is reproduced** [own].
- **"On A's samples the exit's seed is the data's" holds, more sharply than the md shows.** A's exit keeps its time (to 0.4%) and its direction under every change of the arithmetic. The symmetric problem never leaves with exact summation, and leaves the other way when its samples are summed in reverse order. But the printed conclusion is conditioned on directions only, and the offset-start test is a logarithmic magnitude check (D3dd-1).
- **`a2_ulps_exact` is right**: formula, prefactor, kernel, w, interpolation and precision.
- **`mirror_asym` is a correct bottleneck.** Its value is twice the sup-norm distance to the fixed-point subspace, not the distance itself.
- **`mirror_odd`'s on-sample value is the residual's odd part** (equal in all nine ends). Its between-samples value is set by invisible unformed spikes and is grid-limited.
- **(ii)(b)'s "a symmetry breaking of the decomposition, not of the fit" is true**, but what supports it is the formed kernels' function, symmetric to 5.5·10⁻⁹ on and between the samples, and that is not printed (D3dd-2).
- **The lemma's revised statement is right.** The edges are grid floors of roots at 0.8437, 0.4696, 0.7746, 0.6195 and 0.5583.
- **The tests are strong.** All 15 repository mutants are killed, and 19 of my 26. The 7 survivors are equivalent on the data; none changes a published result.
- **No major issue.**

**Findings.**
1. **D3dd-1 (minor): the exit's seed.** The conclusion holds, but its printed evidence and condition are weak.
   - The "data's" clause is printed whenever every A's-samples exit and the offset start go downward. A missing offset start counts as downward. Times are not compared. The symmetric problem's "its seed is the arithmetic's" is printed unconditionally.
   - The offset start uses the held state's double zero (9.9 ulps), not the 50-digit one (9.57).
   - The data's actual seed is 1.9× (5A) and 2.3× (3A) that offset, accrued along a route on which d changes sign twice and λ crosses zero. The times agree only logarithmically.
   - The symmetric problem's "upward" is the summation order's: reversing the samples gives the same times with the opposite sign.
2. **D3dd-2 (minor): the mirror evidence.**
   - "0.24 from the fixed-point subspace" is the bottleneck distance from the mirror image. The sup-norm distance to the subspace is 0.12.
   - The all-kernel odd part on the samples (1.18·10⁻⁹) equals the residual's odd part (1.18·10⁻⁹), and does so at every symmetric end.
   - Between the samples, 4.5·10⁻³ is a 5·10⁻⁴ grid's value. Resolved it is 4.90·10⁻³, from six unformed sub-spacing spikes on one side, of |w| ≤ 0.0098 and φ ≤ 5·10⁻⁶ on the samples.
   - The formed kernels' own function is symmetric to 5.5·10⁻⁹ everywhere: that is the support, and it is unprinted.
3. **D3dd-3 (cosmetic): "(not at the linearisation's −d/λ)" misleads.**
   - ċ is linear in c to 4·10⁻¹⁴ over 16 ulps. With d in 50 digits, −d/λ is the zero (9.5705). The 8.72 is the double sum's 9% error in d.
   - Also, "the centre sample's displacement drifts c": the other pairs' errors offset 30% (5A) and 74% (3A) of it.
4. **D3dd-4 (cosmetic): stale or loose text.**
   - md:54 still says "left once roundoff has seeded the offset", as does `a2e_job`'s docstring.
   - D59 (i) says "downward for every integrator tried" and "scipy's Radau leaving", but omits that scipy's BDF holds 3A.
5. **D3dd-5 (cosmetic): the lemma's edges are grid floors.** "The ±1 term below 0.46" has its root at 0.4696.
6. **D3dd-6 (cosmetic): tests.**
   - Seven equivalent survivors. R09 (A's own exits dropped from the conclusion's condition) is the notable one.
   - 3A's exit is untested.
   - The 0.6-ulp tolerance passes 3A with a 0.03-ulp margin.
   - The exact-start clause and the held list use different "held" criteria (the maximum over outputs against the end state).

## What I checked and how

Scripts are in W/own/ and logs in W/logs/; below they are named by basename only.

1. **Reproducibility** [own].
   - `--from-raw` into W/out took 3.8 s; `cmp` with the committed md: identical (sha256 dd76694b…).
   - pytest: 54 passed in 26.9 s; the md in W/out is unchanged afterwards.
   - rerun.py compares every field except `secs`, bit for bit. Identical:
     - `part_A2E` (all 3 entries);
     - `a2_job` (0.25, 3A), (0.25, 5A), (4, 20A);
     - `h_job` H σ* 0.06 equal and σ* 0.04 equal (the mirror inputs);
     - H3 ±0.05 (905 s) and ±0.15.
2. **Zeros of ċ** (zeros.py, zeros.log) [own]. "Zero" is in ulps of c* above c*.

   | held state | own decimal, 30/50/80 digits | repo `a2_ulps_exact` | f* exact | double (`a2_ulps`) | −d/λ, d in double | −d/λ, d in 50 digits | held |
   |---|---|---|---|---|---|---|---|
   | 5A, span 10 | 8.11834 | 8.118335 | 8.12427 | 8.11835 | 8.10 | 8.118 | 8 |
   | 3A, span 10¹³ | 9.57046 | 9.570462 | 9.64264 | 9.93209 | 8.72 | 9.5705 | 9 (last 10 steps 9 and 10) |

   - **Formula against `critical_rhs_nd`.** At k = 10⁵ ulps: 3.6418844·10⁻⁹ (double) against 3.6418831·10⁻⁹ (decimal); the difference is the double's roundoff, 10⁻¹⁵ absolute.
   - **ċ at c* itself:** −2.9516·10⁻¹³ against −2.9568·10⁻¹³ (5A); −1.9018·10⁻¹⁴ against −2.0881·10⁻¹⁴ (3A: the double sum is 9% off).
   - **w:** w from α in decimal equals `RT.critical_w` (3.9604942091396564).
   - **Linearity:** |second difference|/|first difference| ≤ 4·10⁻¹⁴ over k = 0…15, so the linear interpolation is exact.
   - **Sensitivity:** moving σ by one ulp leaves the zero unchanged to 5 digits.
   - The prefactor Gm·w/(Nσ²) > 0 and the mean's 1/N do not move the zero. The kernel exp(−(x − c)²/(2σ²)) and the residual f* − wφ are `BumpsCW`'s.
3. **Exits from c*** (exits.py, revsym.py; exits_*.log, revsym.log) [own].
   - **`a2_exits` reproduced exactly:**
     - 3A: 0.0483 (−), 0.0478 (−), Radau 0.0567 (−), BDF holds; symmetric 0.0507 (+), 0.0514 (+); offset start 0.0489 (−).
     - 5A: 0.0035 (−) twice, Radau 0.0039 (−), BDF 0.0039 (−); symmetric 0.0040 (+), 0.0039 (+); offset start 0.0035 (−).
     - Note: `kzs` uses `a2_ulps`' double zeros (9.93 and 8.12).
   - **Arithmetic varied, data kept** (DOP853 at rtol 10⁻¹⁰ / 10⁻¹²):

     | right-hand side | A's 5A | A's 3A | symmetric 5A | symmetric 3A |
     |---|---|---|---|---|
     | numpy sum | 0.003500 / 0.003480 (−, −) | 0.04820 / 0.04845 (−, −) | 0.00393 (+) / 0.00367 (−) | 0.05244 (−) / 0.05353 (+) |
     | `math.fsum` | 0.003481 / 0.003481 (−, −) | 0.04840 / 0.04841 (−, −) | **no exit by s = 1** | **no exit by s = 1** |
     | samples reversed | 0.003472 / 0.003493 (−, −) | 0.04821 / 0.04836 (−, −) | 0.00393 (−) / 0.00367 (+) | 0.05244 (+) / 0.05353 (−) |
     | long double, δ state | 0.003482 / 0.003482 (−, −) | 0.04830 / 0.04830 (−, −) | 0.00438 (−) / 0.00447 (−) | 0.0647 (+) / 0.0635 (−) |

   - **The repository's own right-hand side on reversed symmetric samples:** the same times (0.05072, 0.05135, 0.00397, 0.00386), every direction (−) where the code's order gives (+).
   - **Symmetric samples with `fsum`, from offsets:**
     - −8.118 ulps → 0.003529 (−); +8.118 → 0.003529 (+); −81 and −100 → 0.00333 and 0.00332 (−).
     - But −0.81 and −1.0 ulp → 0.00365 and 0.00364 (**+**).
     - 3A: ∓9.570 → 0.04939 (∓); −1.0 → 0.0524 (+).
     - So within about an ulp the exit's direction does not follow the offset's sign.
4. **The seed along the route** (seed.py, seed.log) [own]. I integrated A's problem from c* (DOP853, rtol 10⁻¹², dense output) and took d(s) in 50 digits and λ(s) from the exact Jacobian on 800 states.

   | | 5A | 3A |
   |---|---|---|
   | λ crosses 0 at | s = 6.3·10⁻⁴ (σ/σ* 0.285) | s = 3.3·10⁻³ (σ/σ* 0.271) |
   | min ∫λ | −0.215 | −0.619 |
   | d > 0 (a short interval) | on s ∈ [3.5, 5.7]·10⁻⁴ (σ/σ* 0.56 → 0.33), d changing sign twice | on s ∈ [1.06, 2.68]·10⁻³ |
   | seed ∫e^{−Λ}d | −5.47·10⁻¹⁷ | −7.76·10⁻¹⁷ |
   | offset start's −kz·ulp | −2.82·10⁻¹⁷ (×1.94) | −3.32·10⁻¹⁷ (×2.34) |
   | symmetric problem's summation-residue seed | +3.7·10⁻¹⁹ (data/arithmetic 148) | +1.9·10⁻¹⁹ (419) |
   | stable phase's share of the seed | 35% | 17% |
   | −d/λ at the exit | 8.11 ulps | 9.46 ulps |
   | linear prediction against actual offset at 0.7 s_e | −3.80 against −3.55·10⁻⁹ | −1.62 against −2.27·10⁻⁸ |

   At the λ crossing, |d| (6.8·10⁻¹⁵ at 5A) is no larger than the symmetric problem's summation drift (up to 8.5·10⁻¹⁵). The data dominates through coherence along the route, not through pointwise size.
5. **The drift's make-up** (split_d.py) [own]. At the held states the centre sample's term is 1.43× (5A) and 3.78× (3A) the total drift d; the other pairs' errors offset the rest. Check: 8 × (2840.7/628.9) × (1/3.776) = 9.57, the 3A zero.
6. **The lemma** (lemma.py, lemma.log) [own].
   - My closed form equals central differences of `critical_rhs_nd` to 6 digits (−213.725, +193.723, +14761 at w = 4A; −74.7662, −20.9593, +852.271 at the least-squares w).
   - There is one sign change on [0.15, 1.2] for each predicate:
     - w = 4A: the total at σ/h = 0.8437, the ±1 term at 0.4696; the ±1 term never rises above the centre's;
     - least-squares w: the total at 0.7746, the ±1 term at 0.6195, the ±1 term above the centre's at 0.5583.
   - `a2_lemma_signs` prints the grid floors 0.84, 0.46, —; 0.77, 0.61, 0.55 [md].
   - XS is symmetric to 4.44·10⁻¹⁶, C0 to 2.2·10⁻¹⁶.
7. **Mirror** (mirror.py, mirror2.py, mirror3.py; mirror*.log) [own].
   - **Random tests.** 400 adversarial inputs (n 1–21, quantised ties, duplicates, exact mirrors with one kernel perturbed, negative widths, {a, a, Ra} types): 0 mismatches against my bottleneck. For n ≤ 7, 0 mismatches against brute force over all permutations, and over involutions too.
   - **Why involutions suffice.** D_ii = 2|c_i − c0|, and D_ij ≥ |c_i + c_j − 2c0|. In an odd cycle some kernel has |c − c0| ≤ b/2, so cycles reduce to involutions.
   - **Factor 2.** With the midpoint construction, the sup-norm distance to the fixed-point subspace is b/2: {a, a, Ra} gives b = 0.2 against a distance of 0.1, and H σ* 0.06's 0.24 corresponds to 0.12.
   - NaN input returns NaN (not on the data).

   | symmetric end | bottleneck formed / all | all kernels' odd part: on paired samples, against the residual's odd part, resolved between samples | formed kernels' function, resolved |
   |---|---|---|---|
   | H σ* 0.04 | 2.16·10⁻¹¹ / 2.64·10⁻¹⁰ | 1.6·10⁻¹⁵, 2.2·10⁻¹⁶, 2.2·10⁻¹⁶ | 1.1·10⁻¹⁶ |
   | **H σ* 0.06** | **0.240 / 0.483** | **1.18·10⁻⁹, 1.18·10⁻⁹, 4.90·10⁻³** (grid: 4.51·10⁻³) | **5.5·10⁻⁹** |
   | H3 ±0.05 | 1.65·10⁻³ / 0.423 | 1.25·10⁻⁶, 1.25·10⁻⁶, 1.25·10⁻⁶ | 1.5·10⁻⁴ |
   | H3 ±0.12 | 5.1·10⁻¹³ / 7.55·10⁻⁴ | 2.35·10⁻¹⁴, 2.35·10⁻¹⁴, — | 2.6·10⁻¹³ |
   | H3 ±0.15 | 0.30 / 0.30 | 1.6·10⁻¹⁵, 1.4·10⁻¹⁵, 1.5·10⁻¹⁵ | 1.5·10⁻¹⁵ |
   | HP σ* 0.04 / 0.06 | 1.43·10⁻⁷, 2.85·10⁻⁸ / 0.0122, 3.8·10⁻⁵ | the residual's (1.54·10⁻⁷, 1.02·10⁻¹²) | 1.4·10⁻⁷, 9.0·10⁻¹⁰ |

   - **H σ* 0.06's spikes.** Six active kernels narrower than the spacing (σ 8.9·10⁻⁵ … 1.0·10⁻³, w −0.0001 … −0.0098) sit at c −0.236 … −0.395, all on one side. Their largest φ on the samples is ≤ 5.3·10⁻⁶. The 4.90·10⁻³ maximum is the −0.0098 spike of width 3.3·10⁻⁴ at x = 0.2471.
   - **H σ* 0.06's formed fit.** One kernel at −0.08 (w 1.0000) against a cluster of eight at +0.04 … +0.16 (one +3.239, seven negative). Yet its function equals its mirror image to 5.5·10⁻⁹.
   - **Pairing.** `mirror_odd(xs=XS[XS >= 0])` evaluates at −x, not at the samples. With the actual pairs the printed 1.1·10⁻¹⁶ (σ* 0.04) becomes 1.6·10⁻¹⁵; the difference is immaterial.
   - **The loops' formed kernels (HL):** at most 2.81·10⁻⁷, equal to the md.
8. **Numbers.** Every number added to D59.md and the md since 5ee87f7 is reproduced [own]:
   - the ulps: 8, 8.12, 9.57, 8.1, 9.9, "9, 10";
   - the exits: 0.0483, 0.0478, 0.0035 twice, 0.0567, 0.0039 twice, BDF's hold of 3A, 0.0507, 0.0514, 0.0040, 0.0039, 0.0489, 0.0035;
   - XS's 4.4e-16;
   - the lemma: 0.84, 0.46, 0.77, 0.61, 0.55;
   - the mirror: 0.24, 0.48, 17, 1.2·10⁻⁹, 4.5·10⁻³ (grid), 2.2·10⁻¹¹, 0.3, 1.4·10⁻¹⁵, 1.4·10⁻⁷, 21, 0.012, 2.8·10⁻⁷, and the md's per-target list.
9. **Mutants** (W/mut/res_{r106,r109}.jsonl; mut_*.log).
   - **The repository's spec:** 15 mutants, all killed; the control passes. Z04 and Z13 are killed by the tests alone (planted inputs); the other 13 by the tests and the md.
   - **spec109.py (26 + control): 19 killed.**
   - Killed by the md test, and each one showing that a published number depends on a choice:
     - R04, the offset start from the 50-digit zero;
     - R13, a 10× finer odd-part grid;
     - R15, actual sample pairs;
     - R18, edge step 0.001;
     - R25, reversed symmetric samples: the directions flip.
   - Killed by the tests alone: R05 and R11.
   - **The 7 survivors** (tests pass, md unchanged), all equivalent on the data:
     - R06: a missing offset start fails, instead of passing;
     - R08: implicit exits dropped from the condition;
     - **R09: A's own exits dropped from the condition**;
     - R10: `perm_type='row'`;
     - R12: signed widths in the bottleneck;
     - R14: the odd part without x = 0;
     - R19: edge scan from 0.9.
10. **Literature.** Eight web searches; no page fetched. Sources lists only what I saw.

## Verdicts per part

| claim (D59.md unless noted) | verdict | evidence |
|---|---|---|
| (i) continuum landscape and convergence; the proof's (i) | **holds**; unchanged | re-read |
| (i) "in double the centre sample lies 8 ulps of c* above c*" | **holds** | [own] |
| (i) "the centre sample's displacement drifts c on the line" (23) | **loose**: the other pairs offset 30% / 74% of it (D3dd-3) | [own] item 5 |
| (i) "a repelling near-symmetric state where ċ vanishes, 8–10 ulps above c*" | **holds** (8.118, 9.570 at 30–80 digits) | [own] item 2 |
| (i) "downward for every integrator tried, as does the symmetric problem started at the posed problem's offset: the exit's seed the data's" (25–26) | **conclusion holds**, confirmed sharply (item 3); the cited test is a magnitude check on a log scale; "every integrator" omits that scipy's BDF and LSODA hold 3A (D3dd-1, D3dd-4) | [own] items 3–4 |
| (i) "(c* lies below the problem's repelling state)" | **a simplification**: true in the unstable phase (83% / 65% of the seed); d changes sign twice along the route | [own] item 4 |
| (i) "on exactly symmetric samples from their centre upward …: the arithmetic's" | **the arithmetic's: holds** (exact sums: no exit); **"upward" is the summation order's** | [own] item 3 |
| (i) "the location the problem's, the stability LSODA's step sequence's" | **holds**; scipy's BDF also holds 3A | [md], [own] |
| Lemma D‴.59a, revised: expansion, O((h/σ)²φ₂), the saddle conditions, the least-squares overshoot | **right** | [own] item 6; round 106 |
| Lemma's numbers (−213.7; 0.84, 0.46; 0.77, 0.61, 0.55) | **hold as grid floors** (D3dd-5) | [own] |
| Fenichel clause ("compact pieces where strictly repelling") | **adequate** | — |
| A2E numerics: held offsets, 50-digit zeros, exits | **hold** | [own] |
| md held line "(not at the linearisation's −d/λ)" | **misleading** (D3dd-3) | [own] item 2 |
| (ii)(b) "formed kernels lie 0.24 from the reflection's fixed-point subspace (… 0.48 with all 17)" | **numbers hold; 0.24 is the distance from the mirror image** (0.12 from the subspace) | [own] item 7 |
| (ii)(b) "its function symmetric to 1.2·10⁻⁹ on the samples (4.5·10⁻³ between them, at spikes …)" | **numbers reproduce**; the first is the residual's odd part, the second grid-limited (4.90·10⁻³) | [own] |
| (ii)(b) "a symmetry breaking of the decomposition, not of the fit" | **holds**, for the formed kernels' function (5.5·10⁻⁹ on ℝ); not shown by the printed values (D3dd-2) | [own] |
| (ii)(b) σ* 0.04 2.2·10⁻¹¹; H3 ±0.15 gauge; HP 1.4·10⁻⁷ (all 0.012); loops 2.8·10⁻⁷ | **hold** | [own] |
| `mirror_asym` (the bottleneck) | **correct** (0 mismatches in 400) | [own] |
| `a2_ulps_exact` | **correct** | [own] |
| exact-start clause from the data | **applied**; its criterion differs from the held list's (D3dd-6) | [md], code |
| tests | **strong**; 7 equivalent survivors | [own] item 9 |

## Status of the round-106 issues

| D3da | status | what remains |
|---|---|---|
| 1 (the exit's seed) | **applied**: signed exits, two tolerances, Radau and BDF, the offset start, the wording | The condition is directions-only. The offset uses the double zero 9.9, not 9.57. "Upward" is the summation order's. The mechanism parenthetical is simplified (D3dd-1). md:54 is still stale (D3dd-4). |
| 2 (the hold) | **applied** | D59 (i) omits that scipy's BDF holds 3A (D3dd-4). |
| 3 (the mirror) | **applied**: the bottleneck (verified), per target, all kernels, the odd part | Factor 2. The on-sample odd part is the residual's. The off-sample value is grid-limited. The formed kernels' function is unprinted (D3dd-2). |
| 4 (the lemma's wording) | **applied**, right | The edges are grid floors (D3dd-5). |
| 5 (ulps) | **applied** (8.12, 9.57) | "(not at the linearisation's −d/λ)", taken over from round 106's proposed text, misleads (D3dd-3). |
| 6 (stale text) | **applied** (the C sentence, the saddle qualifier, Fenichel, XS's 4.4e-16) | md:54's "left once roundoff has seeded the offset"; `a2e_job`'s docstring (D3dd-4). |
| 7 (tests) | **applied**: two tests, 15 mutants killed (my rerun) | R06, R08, R09 survive. 3A's exit is untested. 3A's 0.03-ulp margin (D3dd-6). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dd-1 | minor | **The exit's seed: the conclusion is right, but its condition and evidence are weak.** (`a2_lines` exit line; md:61; D59 (i) 25–27, A2E 295–299.) <br>• `dn` checks signs only. `r['symoff'] or (0.0, -1)` passes a missing offset start (R06 survives). Dropping A's own exits (R09) or the implicit ones (R08) survives. "At nearly A's times" is unchecked. "Its seed is the arithmetic's" is printed for "either way" too. <br>• The offset start uses `a2_ulps`' double zero (9.93 at 3A) while the md cites 9.57. <br>• The data's seed ∫e^{−Λ}d is 1.94× / 2.34× the offset used. d changes sign twice along the route, and λ crosses 0 at σ/σ* 0.285 / 0.271. The time test bounds the seed only within a factor of about 2–4 [own]. <br>• Sharper evidence [own]: A's exit is 0.003472–0.003500 (5A) and 0.04820–0.04845 (3A), always (−), over numpy, `fsum`, reversed order and long double at two tolerances. The symmetric problem with `fsum` never leaves. Reversed, it leaves at the same times with the opposite sign. The data's seed is 148× / 419× the summation residue's. | **Code (`a2_exits`, `a2_lines`):** add per w an exact-sum (`math.fsum`) and a reversed-order right-hand side for both problems. Condition the clause on: A's time within 1% and its direction kept; the symmetric problem not leaving with exact sums; `symoff` present with \|t_symoff/t_A − 1\| < 0.05. Take kz from `a2_ulps_exact`. <br>**md / D59:** "on exactly symmetric samples from their centre upward in the code's summation order (downward with the order reversed; with exact sums it stays: the arithmetic's)"; replace "(c* lies below the problem's repelling state)" by "(the data's drift integrated along the route seeds it: 1.9 and 2.3 times the held state's offset, 150–420 times the summation residue)". <br>**Tests:** a planted fake with A's exits upward (kills R09) and `symoff=None` (kills R06). |
| D3dd-2 | minor | **The mirror evidence for "a symmetry breaking of the decomposition, not of the fit".** (D59 (ii)(b) 73–78; md:621; `mirror_asym`/`mirror_odd` docstrings.) <br>• `mirror_asym` is min over π of max_i ‖θ_i − Rθ_{π(i)}‖∞, the distance from the mirror image. The sup-norm distance to Fix(R∘π) is half: 0.12 for H σ* 0.06 [own, proof and brute force]. <br>• The all-kernel odd part on the samples equals the residual's odd part at all nine ends (1.18·10⁻⁹ = 1.18·10⁻⁹, 1.54·10⁻⁷ = 1.54·10⁻⁷, …), so it carries no information beyond the residual. <br>• "4.5·10⁻³ between them" is a 5·10⁻⁴ grid's value; resolved, 4.90·10⁻³. It comes from six unformed one-sided spikes (\|w\| ≤ 0.0098, φ ≤ 5.3·10⁻⁶ on the samples). <br>• The support exists but is unprinted: the formed kernels (0.24 off) give a function symmetric to 5.5·10⁻⁹ on and between the samples. | **D59 (ii)(b):** "the polished end's formed kernels 0.24 from their mirror image (a bottleneck matching; 0.12 from the fixed-point subspace in the sup norm), their function symmetric to 5.5·10⁻⁹ on and between the samples — a symmetry breaking of the decomposition, not of the fit; between the samples six unformed one-sided spikes (\|w\| ≤ 0.0098, invisible on the samples) make the whole state's function 4.9·10⁻³ odd". <br>**Code:** `mirror_odd` also at every kernel's centre and mirror (or adaptive); print the formed kernels' odd part and the residual's beside the on-sample value; pair actual samples (f(XS) against f(XS[::-1])); the docstring says "distance from the mirror image (twice the distance to the subspace)". |
| D3dd-3 | cosmetic | **"(not at the linearisation's −d/λ)" and "the centre sample's displacement".** (md held line; D59 (i) 23, A2E 289–293.) <br>• ċ is linear in c to 4·10⁻¹⁴ over 16 ulps. With d in 50 digits, −d/λ = 9.5705 = the zero. The printed −d/λ (8.72 ulps) comes from d summed in double (−1.902 against −2.088·10⁻¹⁴) [own]. <br>• The centre sample's term is 1.43× / 3.78× the drift; the other pairs' errors offset 30% / 74% [own]. | "within an ulp of the zero of ċ (−d/λ with d in double, 8.7 ulps at 3A: the double sum's error in d; with d in 50 digits −d/λ is the zero)"; D59 (i): "in double the samples' roundoff — the centre sample's 8 ulps, partly offset by the pairs' — drifts c". |
| D3dd-4 | cosmetic | **Stale or loose text.** <br>• md:54 and `a2_lines`:662 still say "a saddle in c … left once roundoff has seeded the offset", which contradicts the revised seed; `a2e_job`'s docstring says the same. <br>• D59 (i) 25–28 says "downward for every integrator tried … scipy's Radau leaving" while scipy's BDF holds 3A (md:61) [md, own]. <br>• `a2_exits`' docstring calls kz "the zero of ċ" (the double one). | Replace the phrase with "left once the data's drift has seeded the offset (in double; on exact symmetric data the summation's)"; D59 (i): "every exit downward (LSODA and scipy's BDF hold 3A, Radau leaves)". |
| D3dd-5 | cosmetic | **The lemma's edges are grid floors.** The roots are 0.8437, 0.4696 (w = 4A) and 0.7746, 0.6195, 0.5583 (least-squares w), one sign change each on [0.15, 1.2]. "The ±1 term below 0.46" understates by 0.01 [own]. | Print the roots (`brentq` between grid points) or the brackets (0.84, 0.85), etc.; D59 lemma and A2E. |
| D3dd-6 | cosmetic | **Tests.** <br>• Survivors R06, R08, R09, R10, R12, R14, R19 are equivalent on the data; none changes a published result. <br>• The exit test covers 5A only. <br>• `abs(a − b) < 0.6` passes 3A at 0.570 (a 0.03-ulp margin). <br>• The exact-start clause uses the maximum over outputs (`mx < 1e-12`) while the held list uses the end state, so a held run with a transient above 10⁻¹²σ* would print "leaving". | Plant inputs for R06/R09 (D3dd-1); add 3A's exit (downward, 0.047 < t < 0.050); widen to 0.75 ulp or test against the alternation's mean; use one criterion for "held". |

## Sources

- [snippet] SciPy, `scipy.sparse.csgraph.maximum_bipartite_matching` ([docs 1.13.1](https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.sparse.csgraph.maximum_bipartite_matching.html)). Unmatched vertices are −1; `perm_type='column'` gives each row's column. This fits R10's equivalence.
- [snippet] Threshold (binary-search) algorithms for bottleneck matching: [arXiv:1405.0900](https://arxiv.org/pdf/1405.0900) ("at each stage, the edges with bigger weight than the threshold are ignored"); Panagiotas et al., ESA 2023 ([LIPIcs.ESA.2023.87](https://drops.dagstuhl.de/storage/00lipics/lipics-vol274-esa2023/LIPIcs.ESA.2023.87/LIPIcs.ESA.2023.87.pdf)). Only the search snippets were seen.
- [snippet] Engel and Jardón-Kojakhmetov, "Extended and Symmetric Loss of Stability for Canards in Planar Fast-Slow Maps" ([arXiv:1912.10286](https://arxiv.org/pdf/1912.10286); [SIAM J. Appl. Dyn. Syst., 10.1137/20M1313611](https://epubs.siam.org/doi/10.1137/20M1313611)). An explicit Runge–Kutta discretisation can prolong the delay of loss of stability at transcritical and pitchfork singularities: a discretisation changing an exit, related to A2E's holds.
- [snippet] Neishtadt and Treschev, survey on stability loss delay, Russian Math. Surveys 76:5 (2021) 883–926 ([mathnet rm10023](https://www.mathnet.ru/eng/rm10023)). Background for an exit after λ crosses 0 along the route (item 4).
- [snippet] Backward Euler's stability region \|1 − z\| ≥ 1 ([Wikipedia](https://en.wikipedia.org/wiki/Backward_Euler_method); [Wisconsin Math 714, lecture 10](https://web.math.wisc.edu/math714/schedule/714_lecture10.pdf)): the md's "BDF-1 damps a positive real eigenvalue beyond hλ = 2".
- [snippet] Petzold, SIAM J. Sci. Stat. Comput. 4 (1983) 136–148, method switching, as cited by [SciPy's LSODA](https://docs.SciPy.org/doc/scipy/reference/generated/scipy.integrate.LSODA.html).
- [snippet] `math.fsum` as Shewchuk's exactly rounded summation, per [PreciseSums (CRAN)](https://search.r-project.org/CRAN/packages/PreciseSums/PreciseSums.pdf) and [accupy](https://pypi.org/pypi/accupy); numpy's pairwise summation only along the fast axis ([numpy PR #3685](https://github.com/numpy/numpy/pull/3685)). This explains why the code's 2-D `K.corr` sums and my 1-D sums leave in different directions on symmetric samples.
- [memory] Golubitsky–Stewart–Schaeffer on fixed-point subspaces; Fenichel 1971 (cited by D59 with round 103's snippets; not re-searched).

## Lesson

A seed belongs to the data when the exit survives a change of the arithmetic that leaves the data alone: exact summation, reversed summation order, wider precision. It belongs to the arithmetic when such a change removes the exit or reverses it.
- Here A's exit kept its time to 0.4% and always went down.
- The symmetric problem stayed with exact sums and flipped with the order.
- A start at a matched offset only checks the order of magnitude on a logarithmic scale; here the seeds it matched differ by a factor of two.

A symmetry claim about a fit needs the function of the kernels that make the fit, on a grid that resolves them:
- on the samples, a fit's odd part is its residual's;
- between them, it is whatever invisible spikes the state carries;
- a bottleneck distance from the mirror image is twice the distance to the fixed-point subspace.
