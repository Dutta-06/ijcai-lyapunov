# Referee report D3dt, round 125: Proposition D‴.59 as revised after round 116

## Scope

- **Export and work directory.**
  - I exported commit 783d55b (branch claude/vigilant-clarke-hf939b) with `git -C … archive 783d55b | tar -x -C W/repo`.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref125.
  - W/repo/venv is a symlink to the repository's venv: Python 3.13.18, numpy 2.5.4, scipy 1.18.1. `np.longdouble` is x87 extended precision (float128, 18 digits).
- **Diff eb22ec4 → 783d55b on the D59 files.**
  - In exp_moving_width.py, `a2_exit_ld` now takes f* as the stored doubles cast to long double. Also changed:
    - `_route_rates` is new (d and λ along a run);
    - `a2_sym_lag` is new;
    - `a2_seed` also returns Λ(s_e);
    - `a2_exits` adds BDF at rtol 10⁻¹², `lag`, `Lam` and `ss`;
    - `a2_lines` changes the S and ratio printout, the linear response, the "within 3·10⁻⁵" clause, `flipa`, `dev` (the largest double deviation), Λ against 11 ln 2, and the symmetric exits to 3 digits;
    - `split_summary` prints the formed/idle split and the margin over min A.
  - The md changes on two lines (the exit line and the mirror block).
  - In D59.md, Remark D‴.59b moves after the proof's ∎ and is rewritten. (ii)(b), the Fenichel clause, the Status and the A2E numerics also change.
  - In the tests, the round-106 exit-clause test gains plants and `test_round112_…` is retightened. There are 55 tests in all.
  - tools/mutants/round116_D59_own.py (W00–W14) is new.
- **Read in full:**
  - D59.md (427 lines);
  - every changed hunk and what it calls: `_LDBumps`, `a2_exit_ld`, `_route_rates`, `a2_sym_lag`, `a2_seed`, `a2_exit`, `a2_exits`, `a2_lines`' exit block, `split_summary`'s mirror block, `critical_rhs_nd`, `critical_w`, `curve_mobility`, `n_formed`;
  - the md, the two tests line by line, the mutant spec and tools/mutants.py;
  - review/REVIEW_round116_D3dk.md and the Round 116 entry of review/REVIEW_LOG.md;
  - `git diff eb22ec4 783d55b` of these files.
- **Settings.**
  - Every heavy process ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with the repository's venv and MW_OUT in W.
  - pytest ran with `--basetemp` in W.
  - At most two heavy processes of mine ran at a time.
  - My copy of the mutant tool (W/own/mutants_tool.py) sets REPO to W/repo and adds `--basetemp` per slot.
- **Process notes.**
  - I stopped two of my own runs with `kill PID`, after reading /proc/PID/cmdline:
    - PID 13227 (`own/resol.py`): its atol of 10⁻³⁰σ* on δ made DOP853 thrash; I reran it at 10⁻²⁰ and 10⁻²²σ*;
    - PID 13556 (`own/perm.py`): its float32 arm was too slow.
  - A float32 timing probe ended by `timeout`. I dropped the float32 arithmetic.
  - I killed nothing by name. Some of my monitors expired or were stopped with TaskStop.
- **The repository is untouched.** I used only read-only git.
  - HEAD has moved to f346010.
  - The working tree holds another session's uncommitted exp_noise_resonance.py, exp_width_mismatch.py and noise_resonance.md.
  - `git diff 783d55b` on the D59 files is empty.
- **Tags.** [md] the committed md or draft; [raw] the committed raw; [own] my computation; [snippet] a web-search result I saw (no page fetched); [memory] not verified this round.

## Summary

**Design.**
1. **Reproduction.** `--from-raw` byte for byte (it reruns `a2_exits` live, including the long-double exits, `a2_seed` and `a2_sym_lag`); the 55 tests; every new number at full precision.
2. **The data's seed S in 40-digit decimal** on the current long-double route (stored f*, Simpson on 2001 points, my closed-form λ). Also the pointwise drift d(s) against the repository's long double.
3. **The double run's own seed, computed independently.** The symmetric runs' effective seeds S(s) = c(s)·e^{−Λ(0,s)} at the solver's steps, in double and in long double. Their log-ratio against `a2_sym_lag`'s Λ.
4. **The seed's distribution within each arithmetic.** Ten arithmetic-equivalent summation orders (permuted samples) in double and in long double: Λ's mean against 11 ln 2, the seeds in units of uσ*, and the exit directions.
5. **What limits "within 3·10⁻⁵".** A's exit and the seeded run at six solver settings (rtol 10⁻¹⁰ … 10⁻¹³, atol on δ down to 10⁻²²σ*), and planted seeds S(1 ± q).
6. **My own long-double right-hand side**, written from the law.
7. **Normal hyperbolicity** of the line at A2E's states and along the route (the exact Jacobian).
8. **Mutants:** the repository's 15 (W00–W14) and 21 of my own (Y00–Y20). All ran with the md byte test deselected, so semantic kills are told apart from md-only ones; the tool's `--from-raw` comparison tells whether a published number changes.
9. **Web searches.**

**Verdicts.**
- **Reproducibility is exact.**
  - md identical (sha256 7577bccb…89b0);
  - 55 passed in 224.1 s (the warning is LSODA's notice);
  - every new number matches the live run [own].
- **The long-double exit is now on the stored data.** `a2_exit_ld`'s f* is the cast of exactly the doubles `a2_exit` and `feature_sum` use (bit-identical) [own]. D3dk-1 is resolved.
- **"Its seed is the data's" holds, decided** [own]:
  - 40-digit S agrees with `a2_seed`'s to 7.4·10⁻⁶ (3A) and 7.4·10⁻⁷ (5A);
  - e^{Λ(s_e)}|S| = 1.0036 and 1.0031 times the exit level;
  - my own implementation gives A's exit 0.0482972/0.0482975 and 0.00348205/0.00348204, and the seeded run within 3.3·10⁻⁶ and 2.8·10⁻⁶.
- **Λ = ∫λ ds between the two exits is the right measure.** It equals the log-ratio of the two runs' own effective seeds to 0.005–0.008 [own]. So "the two seeds' ratio e^Λ" follows.
- **But each Λ is one draw.** Over ten equivalent summation orders Λ ranges over 5.6–12.6, and the exit direction is a coin in either arithmetic. The scaling claim is right only by the mean: 8.98 ± 0.63 and 7.58 ± 0.57 against 7.62, with each seed O(uσ*) (D3dt-1).
- **"Within 3·10⁻⁵" is a valid bound** (observed ≤ 5.1·10⁻⁶ at every setting). "The integrator's resolution" is the wrong cause: no integrator setting moves the floor (D3dt-3).
- **The Remark** is now complete and correctly placed.
- **One latent crash** in `a2_seed` (D3dt-2).
- **No major issue.**

**Findings.**
1. **D3dt-1 (minor): Λ against 11 ln 2, the arithmetic's flip and "13–26%" are single draws of random variables.**
   - The printed range (5.91 … 8.59), the direction flip and the test (one draw within ±2.5) are realisations.
   - Over ten summation orders Λ spans 5.61–12.57 (3A) and 5.32–10.79 (5A), sd 1.98 and 1.79. The direction is upward in 5/10 and 6/10 (double) and 5/10 and 3/10 (long double).
   - Three of ten orders fall outside the test's window at 3A.
   - My own implementation does not flip at 3A.
   - The conclusion holds by the means and by the seeds' O(1) prefactors in uσ*.
2. **D3dt-2 (minor): `a2_seed`'s no-exit branch returns three values; `a2_exits` unpacks four.** The ValueError is demonstrated [own]. The None handling in `a2_lines` is unreachable, and the docstring is stale.
3. **D3dt-3 (cosmetic): "(the integrator's resolution)" names the wrong cause.**
   - The floor (±5·10⁻⁶) is the long-double rounding of ċ, realised differently by each step sequence. Tightening rtol to 10⁻¹³ or δ's atol to 10⁻²²σ* does not move it.
   - Planted S(1 ± 10⁻³) moves the time by the predicted 2.7–3.1·10⁻⁵: the bound resolves S to 10⁻³.
   - A's exit "0.0482972, 0.0482973" carries a noise digit: it is 0.0482971–0.0482974 over the six settings.
4. **D3dt-4 (cosmetic): "the drift at c* in 40-digit arithmetic agrees … to 10⁻⁵".** It is the integrated seed that agrees to 10⁻⁵. Pointwise, d agrees to 5.2·10⁻⁵ and 1.6·10⁻⁵ of max|d|, and d changes sign twice along the route.
5. **D3dt-5 (cosmetic): the Fenichel clause.**
   - Domination holds where the Remark uses it: in the repelling phase, at A2E's states.
   - It fails over the whole attracting phase, where σ contracts faster than c, not only at the sign change.
   - A compact piece persists as a locally invariant manifold.
6. **D3dt-6 (cosmetic): wording and code nits.**
   - The md's "(scipy, rtol 1e-10)" label also covers the BDF run at 10⁻¹².
   - Lemma D‴.59a's ∎ falls mid-paragraph, with (ii)(a) continuing after it.
   - "Over min A" is redundant on equal-amplitude targets.
   - The arithmetic-flip printout reads only rtol 10⁻¹⁰'s directions.
   - "The samples' roundoff" omits f*'s rounding.
7. **D3dt-7 (cosmetic): tests and mutants.**
   - Of the 36 mutants (two of them controls), 16 are killed by semantic tests and 10 only by the md byte test.
   - The 8 survivors are all equivalent on the data; none changes a published number.
   - The round-116 log does not record W09 and W12–W14 as md-only kills.

## What I checked and how

Scripts are in W/own/: vals.py, seeds.py, dprof.py, symseed.py, perm.py, resol.py, own_ld.py, nh.py, spec125.py, mutants_tool.py. Logs are in W/logs/; mutant records in W/mut/r116/res_a.jsonl and W/mut/own/res_b.jsonl.

1. **Reproducibility** [own].
   - `--from-raw`: 30.3 s, then `cmp` against the committed md: identical.
   - pytest: 55 passed in 224.1 s.
   - Live values (vals.py):

     | | 3A | 5A |
     |---|---|---|
     | A, long double, rtol 10⁻¹⁰ / 10⁻¹² | 0.04829717 / 0.04829727 (−, −) | 0.003482048 / 0.003482051 (−, −) |
     | symmetric, long double | 0.0608993 (−) / 0.0642705 (−) | 0.00451514 (−) / 0.00456995 (−) |
     | S | −7.76217·10⁻¹⁷ = −2.587·10⁻¹⁵σ* = −22.37 ulp(c*) | −5.47329·10⁻¹⁷ = −1.824·10⁻¹⁵σ* = −15.78 ulp |
     | S / (−kz·ulp) | 2.3377 | 1.9433 |
     | Λ(s_e); e^{Λ}\|S\| / level | 28.9865; 1.00357 | 29.3354; 1.00305 |
     | seeded run's time against A's | 3.21·10⁻⁶ | 4.11·10⁻⁶ |
     | `a2_sym_lag`, rtol 10⁻¹⁰ / 10⁻¹² | 6.7753 / 8.5914 | 5.9143 / 7.7271 |
     | double, A: DOP853 10⁻¹⁰, 10⁻¹², fsum, rev | 0.0482766, 0.0477874, 0.0481417, 0.0480398 | 0.0034833, 0.0034883, 0.0034579, 0.0034812 |
     | double, symmetric | 0.0507161 (+), 0.0513544 (+) | 0.0039725 (+), 0.0038608 (+) |

   - The largest double deviation is −1.056% (DOP853 at 10⁻¹²), printed "1.1%".
   - "Later": 13.66, 18.37, 20.08 and 25.15%, printed "13–26%".
   - Every md and D59 number matches.
2. **The seed in 40 digits** (seeds.py, dprof.py) [own].
   - **S (Simpson, 2001 points):** −7.762113·10⁻¹⁷ and −5.4732845·10⁻¹⁷. `a2_seed` agrees to +7.4·10⁻⁶ and +7.4·10⁻⁷.
   - **λ:** the repository's central difference is within 1.0·10⁻⁸ and 2.7·10⁻⁸ of my closed form, λ = Gm w/(Nσ²)Σ[rφ(u²/σ² − 1) − wφ²u²/σ²].
   - λ changes sign at s = 3.21·10⁻³ and 6.29·10⁻⁴, and λ_e·s_e = 32.60 and 38.32.
   - **Pointwise d:**
     - |d_LD − d_40| ≤ 1.65·10⁻¹⁸ (5.2·10⁻⁵ of max|d|) at 3A and ≤ 5.1·10⁻¹⁸ (1.6·10⁻⁵) at 5A;
     - the median relative error is 1.6·10⁻⁵ and 1.5·10⁻⁶;
     - d changes sign at s ≈ 9.7·10⁻⁴ and 2.66·10⁻³ (3A) and at 3.3·10⁻⁴ and 5.7·10⁻⁴ (5A), and the largest relative errors (0.16%, 1.5%) sit next to these zeros;
     - the positive part of ∫e^{−Λ}d is 9% (3A) and 2% (5A) of |S|.
3. **The symmetric runs' own seeds** (symseed.py) [own].
   - S(s) = c(s)e^{−Λ(0,s)} at the steps is constant to 10⁻⁵ once |c| ≥ 10⁻¹⁰ of the level.
   - In both arithmetics c stays exactly 0 up to s ≈ 4.8·10⁻⁴ (3A) and 1.3–2.0·10⁻⁴ (5A).

   | run | S_d (double) | S_ld (long double) | ln\|S_d/S_ld\| | `a2_sym_lag` |
   |---|---|---|---|---|
   | 3A, rtol 10⁻¹⁰ | +1.5553·10⁻¹⁷ | −1.7896·10⁻²⁰ | 6.767 | 6.775 |
   | 3A, rtol 10⁻¹² | +1.0176·10⁻¹⁷ | −1.9046·10⁻²¹ | 8.584 | 8.591 |
   | 5A, rtol 10⁻¹⁰ | +2.6249·10⁻¹⁹ | −7.1271·10⁻²² | 5.909 | 5.914 |
   | 5A, rtol 10⁻¹² | +8.8551·10⁻¹⁹ | −3.9238·10⁻²² | 7.722 | 7.727 |

   - The 0.005–0.008 gap is the exit level's nonlinearity: the same 0.4% as the linear response's 1.004.
   - On the line the couplings ∂ċ/∂(α, σ) and ∂(α̇, σ̇)/∂c are ≤ 2.4·10⁻¹³ (nh.py), so the double and long-double routes coincide.
   - So λ along the long-double route is the right integrand, and e^Λ is the seeds' ratio.
4. **Ten summation orders** (perm.py: the samples permuted, the same data, rtol 10⁻¹²; ln S = ln(level) − Λ(0, t_exit)) [own].

   | | ln S, double: mean (sd) | ln S, long double: mean (sd) | paired Λ: range; mean ± sem | upward: double, long double |
   |---|---|---|---|---|
   | 3A | −40.06 (1.05) | −49.04 (1.59) | 5.61–12.57; 8.98 ± 0.63 | 5/10, 5/10 |
   | 5A | −40.76 (1.21) | −48.35 (1.19) | 5.32–10.79; 7.58 ± 0.57 | 6/10, 3/10 |

   - The geometric-mean seeds are 1.20 and 0.59 u_dσ* in double and 0.31 and 0.62 u_ldσ* in long double. Each is O(1) in units of uσ*, so the seed does scale with the unit roundoff.
   - The ±1 spread in ln S matches the ln|N| of a Gaussian amplitude (sd π/√8 = 1.11) [memory].
   - "Later" over the orders: 15.4–37.2% (3A) and 12.6–27.3% (5A).
   - The test's window |Λ − 7.62| < 2.5 fails for 3 of 10 paired orders at 3A (10.38, 12.57, 10.81).
5. **The seeded run's floor** (resol.py) [own].

   | setting | 3A: A's exit; seeded − A | 5A: A's exit; seeded − A |
   |---|---|---|
   | repository (rtol 10⁻¹², atol 10⁻¹⁴σ*) | 0.04829727; +3.2·10⁻⁶ | 0.003482051; −4.1·10⁻⁶ |
   | rtol 10⁻¹², atol(δ) 10⁻²⁰σ* | 0.04829740; −5.1·10⁻⁶ | 0.003482041; +4.4·10⁻⁶ |
   | rtol 10⁻¹², atol(δ) 10⁻²²σ* (6360 steps) | 0.04829729; −1.1·10⁻⁶ | 0.003482038; +2.4·10⁻⁶ |
   | rtol 10⁻¹³, atol(δ) 10⁻²⁰σ* | 0.04829737; −1.2·10⁻⁶ | 0.003482048; −0.5·10⁻⁶ |
   | rtol 10⁻¹¹, atol(δ) 10⁻²⁰σ* | 0.04829711; −1.7·10⁻⁶ | 0.003482045; +1.8·10⁻⁶ |
   | rtol 10⁻¹⁰ | 0.04829717; −5.0·10⁻⁶ | 0.003482048; −4.9·10⁻⁶ |

   - Planted S(1 ± q) with q = 10⁻⁶, 10⁻⁵ or 10⁻⁴ moves the time non-monotonically by at most 4.4·10⁻⁶.
   - q = 10⁻³ moves it by −3.06/+2.86·10⁻⁵ (3A) and −2.80/+2.66·10⁻⁵ (5A), against the predicted ∓q/(λ_e s_e) = 3.07·10⁻⁵ and 2.61·10⁻⁵.
   - So the floor is the long-double rounding of ċ (the same size as item 2's pointwise error), not the integrator's.
   - 3·10⁻⁵ is the response to a 10⁻³ error in S. That still separates the exact from the stored f* (1.6% and 0.4%).
6. **My own long double** (own_ld.py, everything but the solver in np.longdouble) [own].
   - A: 0.0482972/0.0482975 and 0.00348205/0.00348204.
   - Seeded from my 40-digit S: within −3.3·10⁻⁶ and +2.8·10⁻⁶.
   - Symmetric: 0.0641 (+) / 0.0639 (+) at 3A and 0.00438 (+) / 0.00437 (−) at 5A. Against the repository's double runs this is 10–27% later, with no arithmetic flip at 3A.
7. **The double arithmetic on A's samples** [own].
   - The four DOP853 double exits lie 0.04–1.06% before the long-double one: Δln S = +0.014 … +0.34 (λ ≈ 665). The double arithmetic adds 1–40% to the data's seed.
   - On the symmetric samples its seed is 2–23% of |S| over the ten orders.
   - So in double the data's seed dominates by a factor of 2.5–70. That is the 1.1%, and the long double settles it (rounding ~10⁻⁴ of S).
8. **Normal hyperbolicity** (nh.py, the exact Jacobian on symmetric samples) [own].
   - Normal ∂ċ/∂c against the tangential eigenvalues of the (α, σ) block:
     - A2E's held states: 628.9 against −0.0006 and −2875 (3A); 10498 against −0.0012 and −14459 (5A);
     - 5A over the long span: 63173 against +5517 and −0.0012;
     - route, repelling phase: 664–675 against −0.0006 and −2960 (3A); 9493–11006 against −0.0012 and −14757 (5A).
   - In the attracting phase the normal rate is −212 … −229 against tangential −320 … −529 (3A), and −356 … −382 against −541 … −836 (5A). The σ direction contracts faster than c: not dominated.
9. **D59 against the md.** Every new number is the md's and verbatim:
   - (ii)(b): 9 formed, 0.0496; 2 idle; 0.0098, 2% under 0.01·A;
   - the Remark: 1.1%; 0.0482972, 0.0482973; 0.00348205; −2.59 and −1.82·10⁻¹⁵σ*; 2.34, 1.94; 1.004, 1.003; 3·10⁻⁵; 13–26%; 5.91 … 8.59; 7.62; 0.0516;
   - A2E: 0.0516, 0.0041; 0.0609, 0.0643, 0.00452, 0.00457; 0.0482974, 0.00348204.
10. **Mutants** (`-k 'not test_round87_split_mass_bdf_and_the_md'`).
    - Both controls pass (54 passed, 1 deselected).
    - **The repository's spec:**
      - killed by semantic tests: W01–W08 and W11;
      - killed only by the md (published output): W09, W12, W13, W14;
      - survivor: W10 (equivalent: `sym` requires equal amplitudes).
    - **My spec125.py:**
      - killed by semantic tests: Y02 (Λ's double exit on A's samples), Y03 (its long-double route on A's samples), Y04 (∫d instead of ∫λ), Y09 (`flipa` any), Y11 (the deviation over every integrator), Y19 (d at δ = hd), Y20 (the deviation rounded down);
      - killed only by the md: Y08 (Λ at one tolerance twice), Y12 (the resolution at 3·10⁻⁶), Y13 (Λ one point early), Y14 (BDF at 10⁻¹¹), Y15 (the symmetric exits to 6 digits), Y16 (idle counted over all kernels).
    - **Survivors**, all equivalent on the data with the md unchanged:
      - Y01: λ's step 10⁻⁶σ*;
      - Y05 and Y06: Λ on 11 points, and by a left Riemann sum (λ is 664–675 between the exits);
      - Y07: f* left in double (numpy promotes it to long double);
      - Y10: the deviation without the sums' variants (the maximum is DOP853's);
      - Y17: the flip printed from rtol 10⁻¹²'s run;
      - Y18: the seeded run against rtol 10⁻¹⁰'s exit (5.3·10⁻⁶ < 3·10⁻⁵).
    - No survivor changes a published result.
11. **Literature.** Three web searches; no page fetched.

## Verdicts per part

| claim (D59 unless noted) | verdict | evidence |
|---|---|---|
| Remark: "on the stored data with only the arithmetic raised to long double" | **holds**: f* is bit-identical to the double runs' | [own] |
| "A's exit moves by at most 1.1% from double (DOP853, either tolerance and the sums' order or rounding)" | **holds** (1.056%). It is the double arithmetic's own seed: 1–40% of S on A's samples | [own] |
| "the same at two tolerances (0.0482972, 0.0482973; 0.00348205)" | **holds** (2·10⁻⁶). The sixth digit is noise: 0.0482971–0.0482974 over six settings (D3dt-3) | [own] |
| S = −2.59·10⁻¹⁵σ*, −1.82·10⁻¹⁵σ*; 2.34, 1.94 times c*'s distance below the zero of ċ | **holds** | [own] |
| "its linear response e^{Λ(s_e)}\|S\| 1.004 and 1.003 times the exit's level" | **holds** (40 digits: 1.00357, 1.00305). It is a real test: the 0.4% is the exit level's nonlinearity | [own] |
| "planted … gives A's exit within 3·10⁻⁵, the integrator's resolution" | **the bound holds** (≤ 5.1·10⁻⁶). **The cause is misnamed**: it is the long-double rounding, which no setting moves (D3dt-3) | [own] |
| "round 116's drift at c* in 40-digit arithmetic agrees with the long-double one to 10⁻⁵" | **holds for S** (7.4·10⁻⁶, 7.4·10⁻⁷ this round). Pointwise d agrees to 5.2·10⁻⁵ of max\|d\| (D3dt-4) | [own] |
| "exact sums keep c at the centre, reversing the sums mirrors the run" | **identities**, correctly stated | md |
| "the arithmetic flips its direction (upward in double, downward in long double)" | **true of these runs; a coin** (5/10, 6/10 and 5/10, 3/10 upward; my implementation does not flip at 3A) (D3dt-1) | [own] |
| "in long double it comes 13–26% later, by Λ … 5.91 … 8.59 against 11 ln 2" | **numbers right; one draw.** Orders give 12.6–37.2% and Λ 5.3–12.6 (D3dt-1) | [own] |
| "the two seeds' ratio e^Λ" (md) | **follows**: ln\|S_d/S_ld\| equals Λ to 0.008 | [own] |
| "the seed scales with the arithmetic's roundoff" | **holds by the mean** (8.98 ± 0.63, 7.58 ± 0.57 against 7.62; seeds 0.3–1.2 uσ*) | [own] |
| "at 10⁻¹² BDF leaves, at 0.0516" | **holds**, now in the md | [md] |
| Remark after the proof; Status names it | **holds** | re-read |
| Fenichel: "normally hyperbolic — the normal rate nonzero and dominating the tangential ones; it changes sign along the route" | **holds at A2E's states**. It fails over the whole attracting phase, and the persisting object is a locally invariant manifold (D3dt-5) | [own], [snippet] |
| (ii)(b) "9 formed, the narrowest 0.0496, and 2 idle, the narrowest 0.041 …; 0.0098 of A, 2% under 0.01·A" | **holds** | [md], round 116 [own] |
| What is proved and what numerical | **unchanged and correctly labelled**: (i) proved in the continuum; Lemma D‴.59a proved; the Remark numerical | re-read |
| tests | **semantic coverage good**. The Λ test tests one draw (fails for 3 of 10 orders); 10 mutants are guarded only by the md byte test (D3dt-1, -7) | [own] |

## Status of the round-116 issues

| D3dk | status | what remains |
|---|---|---|
| 1 (the data changed) | **resolved**: f* is the stored doubles, bit-identical [own]; the numbers updated; the test windows tightened (W01 is killed by the semantic test) | — |
| 2 (noise digits) | **applied**: "within 3·10⁻⁵"; Λ(s_e) and the linear response; `a2_sym_lag`; the symmetric exits to 3 digits; the tolerance flip only from the data | Λ, "13–26%" and the arithmetic's flip are still single draws (D3dt-1). A's six-digit exit keeps a noise digit, and the resolution's cause is misnamed (D3dt-3) |
| 3 (the ratio's name and units) | **resolved** | — |
| 4 (the decisive step) | **resolved**: 1.1% printed and stated beside the symmetric samples' Λ; the 40-digit check cited | the citation's wording (D3dt-4) |
| 5 ((ii)(b)) | **resolved** | — |
| 6 (structure) | **applied**: the Remark moved; the Status names it; BDF at 10⁻¹² in the md; the domination wording | Lemma D‴.59a's ∎ mid-paragraph; the md's "(scipy, rtol 1e-10)" label; Fenichel's attracting phase (D3dt-5, -6) |
| 7 (tests) | **applied**: W04, W06, W07, W08 and W11 are killed by the new plants; the window is 3·10⁻⁵; V01 recorded | W09 and W12–W14 are killed only by the md byte test and are not recorded as such (D3dt-7) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dt-1 | minor | **Single draws stand where a distribution is needed.** <br>• e^Λ = \|S_d/S_ld\| is confirmed: the runs' own seeds agree with `a2_sym_lag` to 0.005–0.008 [own]. <br>• But over 10 equivalent summation orders, Λ is 5.61–12.57 (3A, sd 1.98) and 5.32–10.79 (5A, sd 1.79). <br>• The exit direction is upward in 5/10 and 6/10 (double) and 5/10 and 3/10 (long double); my implementation keeps 3A upward in long double. <br>• "Later" ranges over 12.6–37.2% across the orders [own]. <br>• The test's \|Λ − 7.62\| < 2.5 fails for 3 of the 10 orders at 3A. <br>• The non-random quantities agree with the claim: the mean Λ is 8.98 ± 0.63 and 7.58 ± 0.57; the seeds are 1.20 and 0.59 u_dσ* and 0.31 and 0.62 u_ldσ* [own]. | **Code (`a2_sym_lag`):** run K (≈ 10) summation orders, permuting the samples (the same data). Return the paired Λ's mean, its sem and the seeds over uσ*. **md (`a2_lines`):** print "Λ over K orders: mean ± sem against 11 ln 2; the seeds (0.3–1.2)uσ*; the direction upward in k of K". Drop "the arithmetic flips it" as evidence, or print it as the coin. **D59 Remark and A2E:** the same. **Test:** assert the mean within 3 sem of 11 ln 2, or the seeds' prefactor in (0.05, 20), not one draw. |
| D3dt-2 | minor | **A latent crash.** `a2_seed` returns `None, ex, (None, 0)` when the long-double run does not leave, while `a2_exits` unpacks `S, exl, syS, Le`: ValueError (demonstrated with `a2_exit_ld` monkeypatched) [own]. The `seed is None` and `.get('Lam') is None` branches in `a2_lines` are unreachable. The docstring still lists three returns. | **Code (`a2_seed`):** `return None, ex, (None, 0), None`, and update the docstring. **Test:** a plant with `a2_exit_ld` returning no exit, asserting the line prints without the seed clause. |
| D3dt-3 | cosmetic | **"(the integrator's resolution)" and a noise digit.** <br>• The seeded run's mismatch is −5.1 … +4.4·10⁻⁶ at rtol 10⁻¹⁰ … 10⁻¹³ and atol(δ) 10⁻¹⁴ … 10⁻²²σ*: tightening does not move it. <br>• S(1 ± q) with q ≤ 10⁻⁴ moves the time non-monotonically by ≤ 4.4·10⁻⁶; q = 10⁻³ moves it by the predicted 2.7–3.1·10⁻⁵ [own]. <br>• So the floor is the long-double rounding of ċ, realised by each step sequence. <br>• A's exit is 0.0482971–0.0482974 and 0.00348204–0.00348205 over the six settings: the sixth digit of "0.0482972, 0.0482973" is noise. | **md and D59:** "within 3·10⁻⁵ (a 10⁻³ change of S; the floor ±5·10⁻⁶ is the long double's own rounding, which no tolerance moves)". Print A's long-double exits to 5 digits (0.048297, 0.0034820) or with the settings' spread. |
| D3dt-4 | cosmetic | **The 40-digit citation.** What agrees to 10⁻⁵ is the seed S: round 116 on the old route; this round 7.4·10⁻⁶ and 7.4·10⁻⁷ on the new one [own]. The drift itself agrees to 5.2·10⁻⁵ and 1.6·10⁻⁵ of max\|d\|, and it changes sign twice along the route; near those zeros the relative errors reach 0.16% and 1.5% [own]. | **D59 Remark:** "the seed S in 40-digit arithmetic agrees with `a2_seed`'s to 10⁻⁵ (referee rounds 116, 125)". |
| D3dt-5 | cosmetic | **Fenichel.** <br>• In the repelling phase, and at A2E's held states, the normal rate dominates: 628.9 against −0.0006 and −2875; 10498 against −0.0012 and −14459 [own]. <br>• In the whole attracting phase (s < 3.2·10⁻³ / 6.3·10⁻⁴), σ contracts faster than c: −212 … −229 against −320 … −529 (3A), −356 … −382 against −541 … −836 (5A) [own]. That is not "only where it changes sign". <br>• A compact piece with boundary persists as a nearby *locally* invariant manifold [snippet]. | **D59 proof:** "… where it is normally hyperbolic (the repelling phase: the normal rate 10²–10⁴ above the tangential ones; in the attracting phase the width contracts faster than c) as a nearby locally invariant manifold (Fenichel 1971)". |
| D3dt-6 | cosmetic | **Wording and code nits.** <br>• (a) The md's "0.0516 (−) (BDF at rtol 1e-12) (scipy, rtol 1e-10)": the label also covers the 10⁻¹² run. <br>• (b) Lemma D‴.59a's ∎ falls mid-paragraph and (ii)(a) continues after it (D3dk-6's second option was not taken). <br>• (c) "Over min A" is printed beside "over max A" on targets that `sym` restricts to equal amplitudes (hence W10). <br>• (d) The arithmetic-flip printout reads only rtol 10⁻¹⁰'s directions (Y17); opposite flips at the two tolerances would be misreported. <br>• (e) The Remark's "the samples' roundoff" omits f*'s double rounding, which is 1.6% and 0.4% of S (round 116). | (a) `a2_lines`: "(scipy, rtol 1e-10 unless stated)". (b) Start "(ii)(a)" on a new paragraph. (c) Drop "over min A", or print it only when the amplitudes differ. (d) Print both tolerances' directions in `flipa`'s clause. (e) "the samples' and f*'s roundoff". Optionally add item 7's "the double arithmetic adds 1–40% to S on A's samples" to the Remark. |
| D3dt-7 | cosmetic | **Tests and mutants.** <br>• 16 kills by semantic tests, 10 by the md alone (W09, W12, W13, W14, Y08, Y12–Y16), 8 equivalent survivors (W10, Y01, Y05, Y06, Y07, Y10, Y17, Y18) [own]. <br>• The md-only kills include deterministic output errors: W13 (S over σ₀), W14 (the level over σ₀), Y14 (BDF's tolerance), Y16 (the idle count), W09 (the narrowest formed kernel). Only the md byte test guards these. <br>• The round-116 log does not mark them. | **Log:** record W09 and W12–W14 as md-only kills. **Test:** cheap semantic asserts in `test_round112_…`: the printed σ* unit (S/(0.25·S0)), the response's level and the BDF entry's rtol. |

## Sources

- [snippet] Fenichel, "Persistence and smoothness of invariant manifolds for flows", Indiana Univ. Math. J. 21 (1971), 193–225, as cited in Jones–Shkoller ([arXiv:math/9807090](https://arxiv.org/pdf/math/9807090)). Persistence of overflowing invariant manifolds.
- [snippet] Statements of Fenichel's theorem for compact manifolds with boundary, in which the persisting slow manifold is *locally* invariant: "Slow Invariant Manifolds of Slow-Fast Dynamical Systems" ([arXiv:2012.06770](https://arxiv.org/pdf/2012.06770)) and Verhulst's notes ([Utrecht](https://webspace.science.uu.nl/~verhu101/SPT07FV.pdf)). These are the basis of D3dt-5.
- [snippet] Kifer, "The exit problem for small random perturbations of dynamical systems with a hyperbolic fixed point", Israel J. Math. 40 (1981) 74–96 ([HUJI](https://cris.huji.ac.il/en/publications/the-exit-problem-for-small-random-perturbations-of-dynamical-syst/)): the exit time ~ λ₊⁻¹ ln ε⁻¹. Almada Monter–Bakhtin, "Normal forms approach to diffusion near hyperbolic equilibria" ([arXiv:1006.3000](https://arxiv.org/pdf/1006.3000)), and [arXiv:1710.03293](https://arxiv.org/pdf/1710.03293): the error term τ_ε − λ⁻¹ ln ε⁻¹ has a limiting *distribution*. This is the theory behind D3dt-1: Λ = ln(ε_d/ε_ld) plus a random O(1) term. It would be a citation for the Remark.
- [snippet] Forest, Goedde and Sinha, Math. Comput. Simul. 37 (1994) ([RePEc](https://ideas.repec.org/a/eee/matcom/v37y1994i4p323-339.html)): roundoff at machine level breaks a symmetry of the data and equations and unleashes an instability absent in the symmetric subspace.
- [memory] Var ln|N| = π²/8 for a standard Gaussian N. Golubitsky–Stewart–Schaeffer, as cited by D59 from earlier rounds; not re-searched.

## Lesson

A seed made by rounding is a random variable.
- Before comparing it with a scaling law, draw it again under an arithmetic-equivalent change: here ten orders of the same sums.
- Then compare the mean, and report the spread.
- One draw of Λ (8.59) passed a ±2.5 test that three of ten equivalent orders fail. The mean (7.58–8.98) and the seeds in units of uσ* (0.3–1.2) show the scaling without depending on the realisation.

Name a resolution after what moves it.
- No integrator setting moved the seeded run's ±5·10⁻⁶; a 10⁻³ change of the seed did.
- So the bound measures the seed, at a floor set by the long double's own rounding.
