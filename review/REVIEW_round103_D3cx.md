# Referee report D3cx, round 103: Proposition D‴.59 as revised after round 100

## Scope

- **Export.** I exported commit 7b92874 (branch claude/vigilant-clarke-hf939b) with `git archive 7b92874 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref103.
- **Base for the diffs.** d0ec681, the commit the round-100 referee reviewed. The diffs run d0ec681 → 7b92874.
- **The raw.** The only raw change since d0ec681 is A2E: the 3A exact-start run now records its steps to s = 10 (201 steps instead of 169, the first 169 identical), and every A2E entry gains a `steps_to` field. All other parts are field-identical [own, scripts/rawdiff.py].
- **Read in full:**
  - review/drafts/D59.md (375 lines).
  - Every line of exp_moving_width.py that changed since d0ec681: `a2_invariance`'s logarithm, `a2_state`, `a2_fixed_point`, `a2_lemma`, the symmetric-line, lemma, held-run and exact-start step lines of `a2_lines`, `a2e_job`'s `steps_to`, A2E_MORE, `null_cluster` and its uses in `outcome` and `split_summary`, `cap_interval`, and `fit_split`'s `cap` and `main`.
  - The code these call: `match`, `outcome`, `pair_index`/`pair_fits`, `end_split`, `h_job`, `h_table`, `report`'s `oc`, `critical_rhs_nd`/`critical_jac_nd`, `BumpsCW.corr`.
  - results/chaos/moving_width.md, tests/test_moving_width.py (51 tests), review/REVIEW_round100_D3cu.md, the Round 100 entry of review/REVIEW_LOG.md, tools/mutants/round100_D59_own.py, tools/mutants.py, the D59 row of review/drafts/README.md, and the diffs of all of these.
- **Settings.**
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, in the repository's venv (numpy 2.5.4, scipy 1.18.1).
  - MW_OUT pointed to W/mwout (a copy of the raw) or to each mutant slot's `out` directory.
  - pytest used `--basetemp` inside W. My copy of the mutant tool (W/scripts/mutants_tool.py) adds `--basetemp` per slot.
  - At most two heavy processes of mine ran at a time, each under `timeout`.
  - My wait loops tested `/proc/PID` of PIDs I had verified as mine through /proc/PID/cmdline; none matched its own command line. I stopped no process.
- **I modified nothing in the repository.** I used only read-only git (`archive`, `show`, `log`, `diff`, `status`). At the end, D59.md, the script, the md, the raw, the tests and the round-100 spec are byte-identical to the commit.
  - Meanwhile another session committed ba0b3d1 … fbefb75 and is editing D‴.60's files, and other sessions' mutant runs (ref101, ref102) were running. None of this touches a D‴.59 file (`git diff 7b92874 HEAD` on them is empty), and I touched none of it.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my own run or computation; [snippet] a web-search result; [memory] not verified this round.

## Summary

**Design.** I worked in this order:
1. Reproduced the md from the raw and ran the tests.
2. Reran all 3 A2E jobs (including the 3A run stepped to s = 10), all 15 A2 jobs, H 0 and H 7, and HL 0 against the raw.
3. Recomputed the round's new quantities with my own code:
   - the lemma's split, by finite differences and an independent formula;
   - the sign of the neighbour terms as σ/h falls;
   - the line's drift at the held states: ċ at c* + k ulps, the same drift in 60-digit decimal arithmetic on the double data, and its share by sample group;
   - the exact-start flow integrated with an explicit method (DOP853) on A's samples and on exactly symmetric ones;
   - the 3A run's BDF steps with their hλ and multipliers.
4. Rebuilt the 94 splits and brute-forced `cap_interval`.
5. Tested `outcome` and `null_cluster` on every stored end and on planted ends.
6. Checked the mirror symmetry of the ends on H's symmetric targets (the multi-kernel question).
7. Ran mutants:
   - the repository's 22 round-100 mutants;
   - round 97's AE2 and AE3;
   - 39 of my own;
   - job-level reruns of the survivors inside mutated copies.
8. Ran six web searches.

**Verdicts.**
- **Reproducibility is exact.** The md matches byte for byte. All 51 tests pass (20.3 s). All 21 rerun jobs are bit-identical to the raw, the rerun A2E included.
- **Every number the revision added to the draft is the md's verbatim.** The md lines changed since d0ec681 are 45, 55, 58–60 and 619. All six of round 100's issues are applied.
- **The cap interval is exact.** The brute-force interval is [1.149818, 2.814436), so the printed (1.1499, 2.8144) is right and rounded inward.
- **The null predicate agrees with `outcome` on all 25 stored path ends (H, H2, H3, HP) and on C's 80.**
- **The lemma's invariance, decoupling, Jacobian formula and centre term are proved and check to 3·10⁻⁸ against finite differences.**
- **What does not hold:**
  - The interpretation of the held runs.
  - Part of the lemma's saddle argument.
  - The claim of one null predicate, in the code's terms.

**Findings.**
1. **D3cx-1 (minor; the principal finding): the held runs sit on the computed problem's own symmetric state, not on "the solver's" state.**
   - In double, the sample that "is c*" is XS[163] = C_A + 2.78·10⁻¹⁷, 8 ulps above C_A.
   - The held offsets are exactly 8 and 9 ulps (+9.252·10⁻¹⁶σ* and +1.041·10⁻¹⁵σ*).
   - The double ċ along the line vanishes there: between 8 and 9 ulps for 5A, and between 9 and 10 ulps for 3A, whose last 13 steps alternate 9 and 10 ulps.
   - The drift d is the data's, not summation roundoff: 60-digit arithmetic on the double XS and f* gives −9.856·10⁻¹² against −9.839·10⁻¹² (5A). About 98.5% of it is the centre sample's 8-ulp displacement.
   - −d/λ is the zero of the flow's linearised ċ, and so the fixed point of every consistent method. BDF supplies only its stability.
2. **D3cx-2 (minor): "the flow from c* exactly leaves the saddle after a time set by the asymmetry's logarithm over λ" is not quantified, and in double it does not separate the two problems.**
   - The md prints no exit time. Only the 5A Adams run leaving is evidence.
   - With DOP853 [own], the flow from c* reaches 10⁻²σ* at s = 0.0035 (5A) and 0.048 (3A) on A's samples.
   - On exactly symmetric samples the times are 0.0040 and 0.051.
   - So in double the exit is set by roundoff, and the lemma's symmetric problem leaves too.
3. **D3cx-3 (minor): Lemma D‴.59a's saddle half is argued with a wrong sign, and the Status overclaims it.**
   - "The neighbours restore … so for σ well below h the symmetric state is a saddle" is wrong in that regime. At w = 4A the ±1 term turns positive for σ/h ≤ 0.45.
   - At the least-squares w, the ±1 term is positive for σ/h ≤ 0.6 and dominates the centre term (+196 against +19.5 at σ/h 0.3).
   - Overshoot is not sufficient near σ ≈ h: at w = 4A and σ/h = 1 the total is −214, which is stable.
   - A2E's states are at σ/h 0.55 and 0.70, so their saddle is numerical, not proved.
4. **D3cx-4 (minor): "one predicate, `outcome`'s and the splits'" holds on the data but not in the code.**
   - `outcome` keeps a second label, `'null cluster'` (a coinciding pair that fits no feature), which `null_cluster` does not govern.
   - It also gates `'exact, null cluster'` on `match`'s nearest kernel.
   - Round 100's planted end 1 is labelled `'null cluster'` by `outcome`, while `null_cluster` says False [own].
5. **D3cx-5 (cosmetic): wording.**
   - "Far beyond BDF-k's bounds" covers three BDF-2 steps at hλ 0.696 < 4. Those steps amplify ×2.17–2.21 and produce the 9.4·10⁻¹³σ* maximum.
   - The |ċ| 2.9·10⁻⁴ at σ*/σ₀ = 4 comes from the unpaired samples, not from "linspace's roundoff".
   - The held offsets are integer ulps, printed to four digits.
   - "The flow stays on the symmetric state" sits next to "the held state is the solver's".
6. **D3cx-6 (cosmetic): the lemma's scope and literature.**
   - The invariance extends to several kernels (a reflection combined with a kernel permutation). The scalar decoupling does not.
   - On H's symmetric targets the polished path ends leave the fixed-point subspace (mirror asymmetry 0.24), while the loops stay on it to 10⁻⁸ [own]. The draft does not say this.
   - D3cu-2's citation for invariance of fixed-point subspaces was not added.
   - H3 ±0.05's shared kernel goes to feature 1 by a tie.
7. **D3cx-7 (cosmetic): the tests.**
   - The repository's 22 round-100 mutants: 21 killed, and V20 equivalent, as logged.
   - AE2 and AE3 are now killed.
   - Of my 39, 24 are killed and 15 survive. None changes a published result: the md from the raw, and job reruns of A2E and H 0, are unchanged.
   - Untested logic remains: the summary's null test on residuals and on HL/HLL rows, the start of `cap_interval`'s scan, and `outcome`'s formed filter.

## What I checked and how

All paths are under W.

1. **Reproducibility.**
   - `--from-raw` into W/mwout (2.2 s), then `cmp` with the committed md: identical.
   - pytest: 51 passed in 20.3 s (logs/pytest_full.log).
   - scripts/rerun.py compares every field except `secs`, bit for bit:
     - `part_A2E`: all 3 entries identical (5A spans 10 and 10¹³; 3A to s = 10 with 201 steps and `steps_to` 10.0).
     - `a2_job`: all 15 identical.
     - `h_job` H 0 (σ* 0.04, equal; its outcome `exact, null cluster` from the new code) and H 7: identical.
     - `hl_job` HL 0: identical (25 s).
     - Logs: logs/rerun_*.log.
2. **Stored outcomes against the current `outcome`** (scripts/outcomes.py). All 25 stored H, H2, H3 and HP outcomes equal the recomputed ones. `null_cluster` is True only on H σ* 0.04, equal. C's 80 give exact 71, coinciding pair 3, dipole run-off 6. The H table prints r['outcome'] stored at job time, so the new predicate reaches it only on a rerun; today the two are identical.
3. **Planted ends** (scripts/planted.py) [own]:

   | end | residual | `outcome` | `null_cluster` | summary |
   |---|---|---|---|---|
   | P1, round 100's #1: a main kernel deficient by 5·10⁻⁷ plus a ±0.3 pair at its θ | 2.7·10⁻¹⁷ | **'null cluster'** (pair branch) | False (others' rms 8.8·10⁻⁸) | none |
   | P2: a ±0.3 pair at a feature's θ, the main kernel 10⁻¹⁰ off (`match` picks the pair) | 2.2·10⁻¹⁰ | **'null cluster'** | True | lists it |
   | P3, P4: exact main kernels plus a far pair cancelling to ≤ 2·10⁻⁸ | ≤ 3.3·10⁻⁹ | 'exact, null cluster' | True | lists it |

4. **The line's drift at the held states** (scripts/drift.py, logs/drift.log) [own].
   - XS[163] − C_A = 2.776·10⁻¹⁷, which is exactly 8 ulps of 0.03 (1 ulp/σ* = 1.157·10⁻¹⁶). All 315 pairs x_k + x_{326−k} − 2c* are nonzero, up to 3.9·10⁻¹⁶.
   - ċ/σ* in double at C_A + k ulps:

     | run | k = 0 | k = 8 (= XS[163]) | k = 9 | k = 10 | zero of ċ | −d/λ | held |
     |---|---|---|---|---|---|---|---|
     | 5A, span 10 | −9.84·10⁻¹² | −1.43·10⁻¹³ | +1.06·10⁻¹² | +2.31·10⁻¹² | 8.1 ulp | 8.10 ulp | 8 ulp |
     | 3A, span 10¹³ | −6.34·10⁻¹³ | −1.19·10⁻¹³ | −7.6·10⁻¹⁴ | +5.6·10⁻¹⁵ | 9.9 ulp | 8.72 ulp | 9 ulp (last 13 steps 9/10) |

   - 60-digit decimal arithmetic on the double XS and f* gives ċ/σ* = −9.856·10⁻¹² (5A) and −6.960·10⁻¹³ (3A). The six unpaired samples contribute 0, and the paired samples carry all of it. So d is the data's, not the summation's.
   - At every ratio [own, scripts/inv4.py], a2_invariance's |ċ| from all samples against paired samples only:
     - σ*/σ₀ = 0.25: 1.50·10⁻¹² and 1.50·10⁻¹² (1.12·10⁻¹³ with c at XS[163]);
     - σ*/σ₀ = 1: 7.2·10⁻¹⁴ and 7.2·10⁻¹⁴;
     - σ*/σ₀ = 4: **2.86·10⁻⁴** and **9.2·10⁻¹⁴**, so at 4 it is the unpaired samples.
5. **The exact-start flow with an explicit integrator** (scripts/explicit.py, logs/explicit.log; DOP853, rtol 10⁻¹⁰, atol 10⁻¹⁴, σ*/σ₀ = 0.25) [own]. The value of s at which |c − c₀|/σ* first passes each threshold:

   | samples, start | w | 10⁻¹² | 10⁻⁶ | 10⁻² |
   |---|---|---|---|---|
   | A's, C_A | 5A | 0.00137 | 0.00264 | 0.00348 |
   | A's, f* and c at XS[163] | 5A | 0.00148 | 0.00275 | 0.00359 |
   | exactly symmetric (hk), 0 | 5A | 0.00186 | 0.00313 | 0.00397 |
   | A's, C_A | 3A | 0.0136 | 0.0344 | 0.0483 |
   | A's, f* and c at XS[163] | 3A | 0.0164 | 0.0371 | 0.0510 |
   | exactly symmetric, 0 | 3A | 0.0161 | 0.0369 | 0.0507 |

   - These match ln(10⁻²/10⁻¹⁵)/λ: 0.0027 and 0.045 after the saddle is reached.
   - In exact arithmetic with exact sample positions, the seed is about 10^{−1205} [own estimate]. That gives s ≈ ln(10^{1203})/λ ≈ 0.25 (5A) and 4.2 (3A), agreeing with round 100's 0.256 and 4.23. The md prints none of these times.
6. **The 3A run's BDF steps** (scripts/steps3a.py) [raw]:
   - Adams up to s = 0.01171.
   - Three BDF-2 steps at hλ 0.696 (BDF-2's bound is 4) grow the offset ×2.17, ×2.21 and ×2.21, from −8.8·10⁻¹⁴ to −9.36·10⁻¹³ (the printed maximum).
   - BDF-1 steps at hλ 6.96–9.44 multiply it by −0.169, −0.159 and −0.165, then +0.139. The last of these matches backward Euler with the drift: (δ + hd)/(1 − hλ) = +6.1·10⁻¹⁶ against +5.8·10⁻¹⁶.
   - From s = 0.21, BDF-2 at hλ 15.6–111 and BDF-3 at hλ 295–1.24·10³ keep the offset at 8–10 ulps.
7. **The lemma** (scripts/saddle.py, logs/saddle.log) [own]:
   - At 5A's largest output, central differences of `critical_rhs_nd` give ∂ċ/∂c = 10888.7217. My sum gives 10888.7221, the same as `a2_lemma` and the Jacobian.
   - The split reproduces: +15702.86, −4919.80, +105.61, +0.052.
   - The centre term equals (Gm w/(Nσ²))(w − A) to all digits. The same holds at 3A (+2901.56, −2488.46, +248.25; total 664.74).
   - On exactly symmetric samples, J[1,0], J[1,2], J[0,1] and J[2,1] are at most 5·10⁻¹⁴: the decoupling holds.
   - The neighbour terms as σ/h falls, on symmetric samples at σ*/σ₀ = 0.25 (format: centre / ±1 / total):

     | σ/h | w = 4A | least-squares w on the line |
     |---|---|---|
     | 1.0 | +2871 / −2817 / **−214** | (1.342A) +110 / −317 / −75 |
     | 0.8 | +4486 / −4600 / +194 | +184 / −369 / −21 |
     | 0.6 | +7975 / −3794 / +4341 | (1.315A) +275 / **+82** / +411 |
     | 0.45 | +14180 / **+576** / +14760 | +194 / +656 / +852 |
     | 0.3 | +31900 / **+759** / +32660 | (1.007A) +19.5 / **+196** / +216 |
     | 0.2 | +71780 / +4.0 / +71780 | (1.0000A) +0.042 / +1.01 / +1.05 |

8. **Splits and caps** (scripts/caps.py, logs/caps.log) [own]:
   - `cap_interval` returns lo = 1.1498178473 and hi = 2.8144359576.
   - Brute force over caps 0.5–4 in steps of 0.005, plus caps at lo·(1 ± 10⁻⁹) and hi·(1 ± 10⁻⁹): the assignments are unchanged exactly on [lo, hi). Just below lo, H3 ±0.05 changes; at hi, HL σ₂ 0.065 10⁻⁸ τ 1600 changes.
   - Over the 94 splits: the largest main (c, σ) error is 1.1498; 30 features are unassigned, 29 of them with an unused same-sign kernel, the nearest 2.8144 off.
   - H3 ±0.05's only positive kernel, (0, 0.129) with w 1.419, is 1.1498 from both features. It goes to feature 1 by enumeration order, a tie.
9. **Symmetric targets** (scripts/symH.py) [own]:
   - XS and C0 are symmetric about 0 to 4.4·10⁻¹⁶.
   - Mirror asymmetry of the formed kernels:
     - H σ* 0.06, equal: **0.24** (the fit's odd part 3·10⁻¹⁰);
     - H2, the same target: 0.23;
     - H3 ±0.15: **0.30** (a pair on one feature);
     - HP (no polish): 3·10⁻⁸;
     - every HL and HLL cell: at most 2.8·10⁻⁷, mostly 10⁻⁸ to 10⁻¹⁵.
10. **Numbers.** I compared every number in the passages changed since d0ec681 with the md: D59.md lines 20–27, 48–60, 148–154, 199–202 and 255–273. All are verbatim, including 10^{−602.3}, 2.3·10⁻³⁸, 4.2·10⁻¹⁵, the drifts, λ, the fixed points and held offsets, the lemma split, the hλ ranges, (1.1499, 2.8144), 29 of 30, 1.15 and 3.9·10⁻¹².
11. **Mutants.** Tool: scripts/mutants_tool.py, which is tools/mutants.py with REPO set to the export, PY set to the venv, and `--basetemp` added. Results in mut/res_{a,b,c}.jsonl and logs/mut_*.log. scripts/jobmut.py reruns `part_A2E` and H 0's `h_job` inside a mutated copy and re-renders the md; its control changes 0 lines (logs/jobmut.log).
12. **Literature.** Six web searches. I fetched no page.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) the landscape and convergence from E(start) < ½‖f*‖² (8–15, 140–148) | **holds**; proved in the continuum; unchanged | re-read |
| (i) A2's durations, escapes and ends (16–27, 250–254) | **holds** | [md]; 15 of 15 bit-identical reruns [own] |
| (i) "on samples symmetric about c* the line … invariant and the state on it a saddle where the kernel overshoots the sample it sits on" (20–21) | **invariance holds**; "a saddle where the kernel overshoots" needs σ well below h and an O(1) overshoot. At σ = h, w = 4A it is stable (D3cx-3) | [own] table |
| (i) "symmetric … but for six … 10^{−602.3} and 2.3·10⁻³⁸ …, and in double linspace's roundoff drifts c" (21–23) | **holds**. The dominant double asymmetry is the centre sample lying 8 ulps above c* (D3cx-1) | [own] drift |
| (i) "LSODA's run of 3A … held at BDF's fixed point of the drift" (24–25) | **the numbers hold**; the fixed point is the flow's own zero of ċ, and BDF's part is its stability (D3cx-1) | [raw], [own] |
| Lemma D‴.59a (148–154): reflection equivariance, invariance of the fixed set, decoupling, ∂ċ/∂c = Gm·w⟨r∂²_cφ − w(∂_cφ)²⟩, the centre term | **proved** (standard equivariance; checked to 3·10⁻⁸ by finite differences) | [own] |
| Lemma: "the neighbours restore at O((h/σ)²e^{−h²/(2σ²)}), so for σ well below h … a saddle" | **the order is right, the sign is not.** The conclusion holds under w(w − f*(c*)) ≫ (h/σ)²e^{−h²/(2σ²)}\|r₁\| or r₁ > 0. At A2E's states it is numerical (D3cx-3) | [own] |
| A2E (255–273): the symmetric state, the Adams growth, the BDF multipliers, the 3A steps, the held-run line, the lemma split | **the numbers hold** (verbatim, rerun bit-identical) | [md], [raw] |
| A2E: "far beyond BDF-k's bounds" (263) | **not for the three BDF-2 steps at hλ 0.696**, which amplify (D3cx-5) | [raw] |
| A2E: "the held state is the solver's, BDF balancing the drift, not the flow's" (269–270) | **does not hold as stated** (D3cx-1) | [own] |
| A2E: "the flow from c* exactly leaves the saddle after a time set by the asymmetry's logarithm over λ" (270–271) | **holds in substance**: in double, shown by the 5A Adams run [md] and DOP853 [own]. The time is not printed, the double exit is roundoff's, and the exactly symmetric samples leave too (D3cx-2) | [own] |
| (ii)(b) the null cluster, "one predicate, `outcome`'s and the splits'" (57–60) | **holds on the data** (one row, H σ* 0.04); **not in the code**: `outcome`'s pair-branch 'null cluster' (D3cx-4) | [own] planted ends |
| (ii)(b) the cap interval (1.1499, 2.8144) by crossing the errors (52–53) | **holds exactly**: [1.149818, 2.814436) | [own] brute force |
| (ii)(b) "at most 1.15 … a budget stop's kernel, the other feature's only one of its sign, left to it by distinctness"; "29 of the 30" (50–52) | **holds**. The kernel is equally far from both features; its assignment is a tie (D3cx-6) | [own] |
| Status (i) (199–202): "Lemma D‴.59a — the symmetric line's invariance and its saddle … proved" | **the invariance is proved; the saddle is not** (D3cx-3) | — |
| (ii)(a), (ii)(c), (iii) and their numerics | **unchanged since round 100**; their md lines are unchanged | diff |

**Mathematics.**
- (i)'s continuum statement is proved.
- The lemma's invariance part is a standard fixed-point-subspace result. It extends to the coupled multi-kernel block: on symmetric (X, f*), the flow is equivariant under c_i ↦ 2c* − c_i combined with a kernel permutation, so the mirror-symmetric configurations, Fix(Σ), are invariant (arXiv:2308.10448 [snippet]).
- The decoupling to one scalar rate needs a single kernel on the axis. With several kernels the odd isotypic block is a matrix: it couples the centres of the self-mirrored kernels with the mirrored pairs' antisymmetric (α, c + c′, σ) modes. The saddle condition then becomes that block's spectrum.
- For the problem as posed, the line is not invariant. The repelling line persists as a nearby repelling invariant set, within O(asymmetry/λ) of it, by persistence of normally hyperbolic manifolds (Fenichel [snippet]); the held run sits on it.
- The saddle's sign at A2E's states, and every exit time, are numerical.

## Status of the round-100 issues

| D3cu | status | what remains |
|---|---|---|
| 1 (not invariant for the problem as posed; held runs BDF's) | **applied**: withdrawn, the exact logarithm printed, steps recorded to s = 10, −d/λ printed, "into roundoff" relabelled | The fixed point is the flow's own symmetric state, 8–10 ulps above C_A on the centre sample; BDF supplies only its stability (D3cx-1). The exit time is not printed (D3cx-2). The ratio-4 attribution is wrong (D3cx-5). |
| 2 (the lemma) | **applied**; the split is printed and reproduces | The sign of the neighbour terms and the hypotheses (D3cx-3). The citation for fixed-point-subspace invariance is not added (D3cx-6). |
| 3 (one null predicate) | **applied** for 'exact, null cluster' and the summary; planted ends and R19 tested | `outcome`'s second 'null cluster' label and its `match` gate (D3cx-4). |
| 4 (tests) | **fixed** | AE2 and AE3 are killed by the job-level NQU assert [own]. The ratio-4 checks are out of the comment. R19 is covered. The log is corrected. |
| 5 (the gap's justification) | **fixed** | (1.1499, 2.8144) is exact. "29 of the 30" and the 1.15 shared kernel are right [own]. |
| 6 (wording) | **fixed** | "f*/A 10^{−602.3} (0 in double)", "the number of the residual's test", and the README's round-83 phrase. |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cx-1 | minor | **The held state is the computed problem's own symmetric state; the solver supplies its stability, not its location.** (D59.md:24–25, :268–271; md:60 "the held state is the solver's, not the flow's"; `a2_fixed_point`'s docstring; README.) <br>• In double the centre sample is XS[163] = C_A + 2.776·10⁻¹⁷ = C_A + 8 ulps. <br>• The held offsets +9.252·10⁻¹⁶σ* and +1.041·10⁻¹⁵σ* are exactly 8 and 9 ulps. <br>• The double ċ on the line changes sign between 8 and 9 ulps (5A: −1.43·10⁻¹³ to +1.06·10⁻¹²) and between 9 and 10 ulps (3A: −7.6·10⁻¹⁴ to +5.6·10⁻¹⁵). The 3A run's last 13 BDF steps alternate 9 and 10 ulps. <br>• d is the data's: 60-digit arithmetic on the double XS and f* gives −9.856·10⁻¹² and −6.960·10⁻¹³ against −9.839·10⁻¹² and −6.339·10⁻¹³ in double. The six unpaired samples contribute 0. With c at XS[163], 5A's ċ is 69× smaller. <br>• −d/λ is the zero of the linearised ċ. Any consistent method's steady state solves ċ = 0. What is BDF's is \|R(hλ)\| < 1 beyond 2, 4 and 6.67. <br>• Read so, the problem as posed has a repelling near-symmetric state (persisting the line, Fenichel). The flow from C_A starts 8 ulps off it and leaves; BDF makes it attracting. | **D59 (i) and A2E:** "held on the computed problem's own symmetric state — where ċ vanishes, 8–10 ulps above c*, essentially on the centre sample XS[163] = c* + 2.8·10⁻¹⁷ — a repelling state that the flow from c* leaves and that BDF, at hλ beyond its bounds, makes attracting". <br>**Code (`a2_lines`):** print XS[163] − C_A in ulps, the zero of ċ (bisect over ulps), and the held offsets in ulps. Consider measuring the offset from XS[163]. <br>Drop "not the flow's". |
| D3cx-2 | minor | **The exit "after a time set by the asymmetry's logarithm over λ" is unquantified, and in double it does not distinguish the symmetric problem from the posed one.** (D59.md:23, :270–271; md:45.) <br>• The md prints no exit time. Only the 5A Adams run (growth 1.092·10⁴ = λ; plateau from s = 0.0056) is evidence. <br>• DOP853 from c* [own]: 10⁻²σ* at s = 0.0035 (5A) and 0.048 (3A) on A's samples; 0.0040 and 0.051 on exactly symmetric samples; 0.0036 and 0.051 from XS[163]. <br>• In double the seed is roundoff (about 10⁻¹⁵σ*) for both problems. Only exact arithmetic separates them: s ≈ 0.25 and 4.2 for the posed problem [own estimate, round 100's 0.256 and 4.23], and never for the symmetric one. | **Code:** add an explicit-integrator (non-damping) exit time for 3A and 5A, on A's samples and on symmetric samples, to the A2E lines. <br>**Draft:** "in double the flow from c* leaves at s ≈ 0.0035 (5A) and 0.05 (3A) on A's samples and, by roundoff, on exactly symmetric ones alike; in exact arithmetic the posed problem's 10^{−1205} seed gives s ≈ 0.26 and 4.2 (an estimate)". |
| D3cx-3 | minor | **Lemma D‴.59a's saddle half is argued with a wrong sign, and the Status overclaims it.** (D59.md:20–21, :152–154, :199–200; README "proves … its saddle".) <br>• The ±1 term is 2(Gm w/N)[r₁φ₁(h²/σ⁴ − 1/σ²) − wφ₁²h²/σ⁴]. For σ ≪ h its leading part has the sign of w·r₁: destabilising when the kernel undershoots its neighbours. <br>• [own], σ*/σ₀ = 0.25: at w = 4A the ±1 term is positive for σ/h ≤ 0.45. At the least-squares w it is positive for σ/h ≤ 0.6 and dominates the centre (+196 against +19.5 at σ/h 0.3). <br>• Overshoot is not sufficient near σ ≈ h: at w = 4A and σ/h = 1 the total is −214, stable. <br>• At a stationary w on the line, w − f*(c*) = 2φ₁r₁ + … is exponentially small, so the "overshoot" mechanism does not carry the saddle there. <br>• A2E's states are at σ/h 0.55 and 0.70, not "well below h": their sign is the Jacobian's. <br>• "w > f*(c*)" assumes w > 0. | **Lemma:** "∂ċ/∂c = (Gm w/(Nσ²))[(w − f*(c*)) + 2φ₁r₁((h/σ)² − 1) − 2wφ₁²(h/σ)² + O(φ₂)], φ₁ = e^{−h²/(2σ²)}; for σ ≪ h positive if w(w − f*(c*)) ≫ (h/σ)²φ₁\|r₁\| (a transient overshoot) or, at the stationary w, if r₁ > 0 (then the neighbours' term leads); at A2E's states (σ/h 0.55, 0.70) the sign is the Jacobian's". <br>**(i):** add "for σ well below h". <br>**Status and README:** "invariance and decoupling proved; the saddle's sign asymptotic under these hypotheses, at A2E's states numerical". |
| D3cx-4 | minor | **`outcome` keeps a second null-cluster label outside the shared predicate.** (D59.md:58–59; `outcome` and its docstring.) <br>• `'null cluster'` is still returned for a residual below 10⁻⁸ with a coinciding pair that fits no feature, and `null_cluster` does not govern it. <br>• `'exact, null cluster'` additionally requires `match`'s nearest kernel within EXACT_TOL. <br>• Planted [own]: P1 (round 100's end 1) gives `outcome` 'null cluster' with `null_cluster` False and summary none. P2 gives `null_cluster` True (the summary lists it) with `outcome` 'null cluster', not 'exact, null cluster'. <br>• On the data no end takes the pair-branch label. The 25 stored outcomes equal the recomputed ones, and H's table prints outcomes stored at job time. | Rename the pair branch (e.g. 'coinciding pair off the features'), or decide it with `null_cluster`. Update the docstring. Test `outcome` on P1 and P2. Or narrow the draft to "one predicate for 'exact, null cluster' and the splits". |
| D3cx-5 | cosmetic | **Wording.** <br>• D59.md:263 and md:59, "far beyond BDF-k's bounds": the three BDF-2 steps at hλ 0.696 (< 4) amplify ×2.17–2.21, producing the 9.4·10⁻¹³σ* maximum. Damping starts with the BDF-1 steps at hλ 6.96. <br>• D59.md:267 and md:45, "in double linspace's roundoff drifts c …: \|ċ\| … 1.5·10⁻¹², 7.2·10⁻¹⁴ and 2.9·10⁻⁴": at σ*/σ₀ = 4 it is the unpaired samples (9.2·10⁻¹⁴ without them) [own]. <br>• The held offsets are integer ulps quoted to four digits. <br>• D59.md:261, "the flow stays on the symmetric state to s = 10", sits beside "the held state is the solver's". <br>• "c* a sample": in double the sample is 8 ulps off. | "Three order-2 steps at hλ 0.70 follow the growth, then all beyond the bounds"; "the unpaired samples at 4"; offsets in ulps; "the solution stays"; "c* a sample (in double 8 ulps above)". |
| D3cx-6 | cosmetic | **The lemma's scope, and the literature.** <br>• Invariance extends to several kernels (a reflection with a kernel permutation); the scalar decoupling needs one kernel. <br>• On H's targets, symmetric to 4.4·10⁻¹⁶, the polished path ends are mirror-asymmetric: 0.24 (H σ* 0.06, equal), 0.23 (H2), 0.30 (H3 ±0.15). HP and the HL/HLL loops stay on the fixed-point subspace to 10⁻⁷–10⁻¹⁵ [own]. "The cancellation is the polish's" is thus a symmetry breaking the draft does not name. <br>• D3cu-2's citation for invariance of fixed-point subspaces is missing. <br>• H3 ±0.05's shared kernel is 1.1498 from both features and is feature 1's by enumeration order. | Add to the proof: "the invariance holds for any number of kernels on Fix of the reflection with the kernel permutation; the decoupling needs one". Note the symmetric targets' symmetry-broken polished ends in (ii)(b). Cite Golubitsky–Stewart–Schaeffer II [memory] or arXiv:2308.10448 [snippet], and Fenichel 1971 for the persistence. Say "a tie" for H3 ±0.05. |
| D3cx-7 | cosmetic | **The tests.** <br>• Repository's round-100 spec: 22 mutants, 21 killed, V20 equivalent, as logged [own]. <br>• Round 97's AE2 and AE3 are now killed. <br>• My 39 (spec103 and W51): 24 killed, 15 survive: <br>  – W03, W11: equivalent (w > 0; z[1] = C_A); <br>  – W16, W18, W25: boundaries; <br>  – W23: positive widths; <br>  – W26: `outcome`'s formed filter dropped; <br>  – W29: the summary's H null test without the residual; <br>  – W30, W31: HL and HLL never null; <br>  – W45: `null_cluster` at cap 10; <br>  – W32: `cap_interval` given all active kernels; <br>  – W36: hi scanned from 2.5; <br>  – W37: strict cap; <br>  – W47: lo default. <br>• Job-level reruns (W03, W16, W18, W23, W26) change 0 md lines, and H 0 stays 'exact, null cluster'. **No survivor changes a published result.** | Add planted tests: an H row with exact main kernels, small others and residual 1.1·10⁻⁷ through `split_summary` (W29); an HL and an HLL null row (W30, W31); a `cap_interval` case whose assignment changes at a cap in (1.5, 2.5) (W36); an unformed kernel in `outcome`'s null path (W26). |

## Sources

- **[snippet]** Fixed-point subspaces are flow-invariant for equivariant vector fields:
  - "A Bifurcation Lemma for Invariant Subspaces", [arXiv:2308.10448](https://ar5iv.labs.arxiv.org/html/2308.10448): proves that Fix(Σ) is F-invariant when F is Γ-equivariant.
  - [arXiv:2206.04442](https://arxiv.org/pdf/2206.04442): states Fix(Γ) invariance as a proposition.
  - **[memory]** Golubitsky–Stewart–Schaeffer, *Singularities and Groups in Bifurcation Theory* II, as the standard text. Not surfaced by the searches.
- **[snippet]** Persistence of normally hyperbolic invariant manifolds: N. Fenichel, "Persistence and smoothness of invariant manifolds for flows", Indiana Univ. Math. J. 21 (1971/72) 193–226, via the [Wikipedia summary](https://en.wikipedia.org/wiki/Normally_hyperbolic_invariant_manifold) ("NHIMs and their stable and unstable manifolds persist under small perturbations").
- **[snippet]** Implicit methods altering stability at large steps: Varsakelis–Anagnostidis, "On the susceptibility of numerical methods to computational chaos and superstability", Commun. Nonlinear Sci. Numer. Simul. 33 (2016) 118–132 ([record](https://research.dial.uclouvain.be/entities/publication/0e8a9d84-6b73-4a25-81ce-d616a85934a2)). Its maps for forward and backward Euler bifurcate at large δt, the solutions are spurious, and backward Euler suppresses chaotic behaviour. I saw the abstract only. Backward Euler's factor 1/(1 − hλ) is from round 100's snippet sources.
- **[snippet]** Roundoff breaking symmetric solutions:
  - [arXiv:2106.01872](https://arxiv.org/pdf/2106.01872): floating-point non-associativity; symmetry cannot be preserved unless operations at symmetric positions are identical.
  - [arXiv:1003.0724](https://arxiv.org/pdf/1003.0724) (Fyris Alpha): floating-point noise breaks the perfect symmetry of the initial conditions.
- **[memory]** The draft's other references are unchanged since round 100, which confirmed Ho–Nguyen and Xu–Du [snippet in that round].

## Lesson

A state that a stiff solver holds is usually the problem's own. Before calling it "the solver's", find where ċ vanishes on the double grid. Here that zero was the centre sample, which in double lies 8 ulps above the nominal c*, and that displacement was almost the whole "drift". The implicit method contributed only the stability of that state.

Two related lessons:
- An asymptotic sign argument ("for σ well below h") does not cover states at σ ≈ 0.55–0.7h. Tabulate each term's sign across σ/h before relying on it.
- In floating point, "invariant on exactly symmetric samples" does not keep a computed flow on a saddle either: summation roundoff seeds both problems alike.
