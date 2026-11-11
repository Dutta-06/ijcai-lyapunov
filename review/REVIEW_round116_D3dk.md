# Referee report D3dk, round 116: Proposition D‴.59 as revised after round 112

## Scope

- **Export and work directory.** I exported commit eb22ec4 (branch claude/vigilant-clarke-hf939b) with `git -C … archive eb22ec4 | tar -x -C W/repo`.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref116.
  - W/repo/venv is a symlink to the repository's venv: Python 3.13.18, numpy 2.5.4, scipy 1.18.1.
  - The long double is x87 extended: `np.finfo(np.longdouble)` gives eps 1.084·10⁻¹⁹ and nmant 63 (a 64-bit significand).
- **Diff cac7a71 → eb22ec4 on the D59 files.**
  - exp_moving_width.py (+132/−):
    - new: `_LDBumps`, `a2_exit_ld`, `a2_seed`;
    - `a2_exits` gains an `ld` block;
    - the exit line gets new `data`/`arith`/`flip`/`late` conditions and a long-double sentence;
    - the mirror sentence gains the width split and the threshold margin;
    - also changed: `mirror_asym`'s docstring, `edge(pred)`, and the mirror block's lead.
  - The md changes at two lines (the exit line and the mirror block).
  - D59.md changes in four hunks: (i), (ii)(b), Remark D‴.59b, and the A2E numerics.
  - tests/test_moving_width.py: the round-106 exit-clause test is extended, and `test_round112_…` is new (55 tests in all).
  - tools/mutants/round112_D59_own.py is new.
- **Read in full:**
  - D59.md (418 lines);
  - every changed hunk of the script and the code it calls: `critical_rhs_nd`, `critical_w`, `curve_mobility`, `_alpha_of`, `coef_crit`, `BumpsCW`, `a2_exit`, `a2_ulps_exact`, `split_summary`;
  - the md, the two tests line by line, the mutant spec, tools/mutants.py;
  - review/REVIEW_round112_D3dg.md and the Round 112 entry of review/REVIEW_LOG.md;
  - `git diff cac7a71 eb22ec4` of these files.
- **Settings.**
  - Every heavy process ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, with MW_OUT = W/out (a copy of the raw) or a mutant slot's `out`.
  - pytest ran with `--basetemp` inside W.
  - At most two heavy processes of mine ran at a time.
  - My copy of the mutant tool (W/own/mutants_tool.py) sets REPO to W/repo and adds `--basetemp` per slot.
- **Process notes.**
  - I stopped two of my own waiter shells, PIDs 4499 and 4583: their `until ! pgrep -f "mutants_tool.py"` matched their own command lines. I checked each via /proc/PID/cmdline and used `kill PID`.
  - Some other waiters of mine reached the harness's time limit after their jobs had finished.
  - One light `python -c` of mine ran without PYTHONDONTWRITEBYTECODE and wrote `__pycache__` into W/repo (my export only); I deleted it.
  - The repository's experiments/__pycache__ entries dated 21:23 belong to another session's gh_peano work.
- **The repository is untouched.** I used only read-only git.
  - HEAD moved to 1ad8abb during the round, and another session holds uncommitted gh_peano files.
  - `git diff eb22ec4` on the D59 files and src/core is empty.
- **Tags.** [md] the committed md or draft; [raw] the committed raw; [own] my computation; [snippet] a web-search result I saw (no page fetched); [memory] not verified this round.

## Summary

**Design.**
1. Reproduction:
   - `--from-raw` byte for byte. It recomputes `a2_exits`, so the long-double exits and the seed are rerun live.
   - The 55 tests.
   - The long-double exits and `a2_seed` at full precision, and the seed at n = 2001 … 16001.
2. An independent seed in 40-digit decimal arithmetic on the long-double route:
   - d(s) and λ(s) in closed form, summed over A's samples;
   - Simpson rule on 4001 and 2001 points;
   - f* taken two ways: the exact Gaussian on the double samples (what `a2_exit_ld` computes), and the stored double f* (the posed problem, and `a2_ulps_exact`'s data).
3. My own long-double right-hand side, written from the law and not from the repository's classes:
   - f* exact or stored, at four (rtol, atol) settings;
   - on A's samples and on the symmetric ones.
4. The repository's long-double code with the stored f*. Cross-planting each data choice's seed into the other.
5. A negative control: the planted test applied to the symmetric problem's own (arithmetic) seed.
6. Λ between the double and long-double symmetric exits, against 11 ln 2.
7. scipy's BDF at rtol 10⁻¹². The 50-digit zero of ċ with both f*.
8. The accuracy of numpy's long-double `exp`; whether x − c* is exact in double.
9. Mirror: my own resolved odd part, by width and by threshold.
10. Mutants:
    - the repository's 12 plus control;
    - 21 of my own plus control;
    - the mutants killed only via the md were rerun with the md byte test deselected.
11. Web searches.

**Verdicts.**
- **Reproducibility is exact.**
  - The md is identical (sha256 63168657…7dd).
  - 55 tests pass.
  - The live long-double numbers match the md to every printed digit [own].
- **"Its seed is the data's", decided in long double, holds** [own]:
  - My 40-digit seed equals `a2_seed`'s to 1.9·10⁻⁶ (3A) and 1.0·10⁻⁵ (5A).
  - The linear response accounts for the exit: e^{Λ(s_e)}|S| = 1.0035 and 1.0030 times the exit level.
  - Two independent long-double implementations agree on A's exit to 2·10⁻⁵ when they use the same f*.
- **`_LDBumps` is long double where it matters.** The Gaussians, products and sums are in x87 extended precision, and `exp` is accurate to about 4·10⁻²⁰. The state, w, the mobility, the returned correlations and the solver are double, which is harmless.
- **But `a2_exit_ld` recomputes f* in long double.** The data change, against the docstring (D3dk-1).
- **λ and the 8001-point trapezoid are more than accurate enough.**
- **The new printed comparisons carry noise digits (D3dk-2):**
  - the seeded run "to 1.2·10⁻⁶ and 1.8·10⁻⁶";
  - "15–25% later";
  - "(3A in long double)".
- **The width split holds:** 11 kernels, 1.19·10⁻⁹.
- **The Remark's move is clean** but placed inside the proof (D3dk-6).
- **Tests are strong on the semantics.** Two of the round's "killed" mutants (V01, V10) are killed only by noise digits of the md, and V01 is not in double at all (D3dk-7).
- **No major issue.**

**Findings.**
1. **D3dk-1 (minor): the long-double run changes the data.**
   - `a2_exit_ld` computes f* as the exact Gaussian in long double. The double problem, the stored raw and `a2_ulps_exact` (the 9.57 and 8.12 ulps) use the double-rounded f*.
   - The rounding of f* is −1.55% (3A) and −0.39% (5A) of the seed.
   - With the posed problem's f*, the long-double exit is 0.0482972/0.0482973 (3A) and 0.00348205 (5A), not 0.0482742/0.0482736 and 0.00348169. S is −7.762·10⁻¹⁷ and −5.473·10⁻¹⁷.
   - The printed "2.37" divides an exact-f* seed by a stored-f* zero. Consistent ratios are 2.36 (exact f*) or 2.34 (stored).
   - The planted test resolves the difference: the cross-planted seed misses A's time by 4.9·10⁻⁴ at 3A.
2. **D3dk-2 (minor): noise digits published as results.**
   - The seeded run's "to 1.2·10⁻⁶ / 1.8·10⁻⁶" is the integrator's floor. It is 1.2–6.4·10⁻⁶ as n varies, while S moves by at most 7·10⁻⁶ relative, which would shift the time by at most 2·10⁻⁷.
   - The symmetric samples' six-digit exits, "15–25% later" and "the tolerance flips (3A in long double)" are one realisation of rounding noise.
   - Arithmetic-equivalent changes give 9.6–30.5% and other flips: the mutant V01 (σ² in double), my mutant X05 (absolute atol), the stored f*, my implementation.
   - The noise-free quantities are ΔΛ = 6.3–8.8 against 11 ln 2 = 7.62, and e^{Λ}|S|/level = 1.003.
3. **D3dk-3 (cosmetic): "2.37 and 1.95 times the held offset".**
   - The denominator is c*'s distance below the zero of ċ (9.57 and 8.12 ulps), not the held offset (9.00 and 8.00 ulps, which would give 2.53 and 1.98).
   - The sign is opposite.
   - S is printed in units of c, while the held offsets are printed in σ*.
4. **D3dk-4 (cosmetic): the Remark's argument leaves its decisive step implicit.**
   - It cites the tolerance-independence and the planted run. It does not state that A's exit moves by ≤ 1% from double to long double while the symmetric samples' moves by ΔΛ ≈ 11 ln 2.
   - Nothing in the code shows that the long-double drift is the data's (my 40-digit seed does).
5. **D3dk-5 (cosmetic): (ii)(b).**
   - The margin that matters is 0.0098 against the threshold 0.01 (2%), not against 0.06.
   - The "11 kernels" are the 9 formed ones plus 2 idle ones (|w| 1.1·10⁻¹¹ and 1.3·10⁻⁸); the "narrowest 0.041" is idle.
6. **D3dk-6 (cosmetic): structure and wording.**
   - Remark D‴.59b sits inside the Proof, and (ii)(a)'s proof continues in the same paragraph.
   - The Status does not name the Remark.
   - "At 10⁻¹² BDF leaves" is not in the md (I verified it: 0.05156, downward).
   - The Fenichel clause equates normal hyperbolicity with a nonzero normal rate.
7. **D3dk-7 (cosmetic): tests and mutants.**
   - V01 "_LDBumps in double" stays in long double (numpy promotes to float128). It and V10 survive once the md byte test is deselected.
   - Of my 21 mutants, 8 survive. X10 (the long-double exits' direction unchecked) has no planted input. The other seven are equivalent on the data.
   - The test's 10⁻⁵ window for the seeded run has a margin under 2×.

## What I checked and how

Scripts are in W/own/: ld_repo.py, seed_dec.py, own_ld.py, misc.py, late.py, kz_exact.py, mirror_chk.py, v01.py, lin.py, spec116.py and mutants_tool.py. Logs are in W/logs/ and mutant results in W/mut/*/res_*.jsonl (W as in Scope).

1. **Reproducibility** [own].
   - `--from-raw` took 19.8 s; `cmp` against the committed md: identical.
   - pytest: 55 passed in 140.4 s. The one warning is LSODA's convergence notice, as before.
   - The live long-double numbers at full precision (ld_repo.py):

     | | 3A | 5A |
     |---|---|---|
     | A's exit, rtol 10⁻¹⁰ / 10⁻¹² | 0.04827419 / 0.04827362 (−, −) | 0.0034816933 / 0.0034816939 (−, −) |
     | symmetric samples | 0.0633521 (+) / 0.0624181 (−) | 0.00477875 (−) / 0.00444005 (−) |
     | S | −7.884721·10⁻¹⁷ (−2.628·10⁻¹⁵σ*) | −5.494592·10⁻¹⁷ (−1.832·10⁻¹⁵σ*) |
     | S / (−kz·ulp) | 2.3746 | 1.9508 |
     | seeded run's time against A's | 1.22·10⁻⁶ | 1.76·10⁻⁶ |

     Every printed md and D59 number matches.
   - n = 2001, 4001, 16001 at 3A: S moves by −6.7·10⁻⁶, −1.3·10⁻⁶ and +1.5·10⁻⁶ relative. The seeded run's mismatch is 1.21·10⁻⁶, 6.39·10⁻⁶ and 2.36·10⁻⁶.
2. **`_LDBumps` end to end** [own reading and runs].
   - **In long double:**
     - xr = x − c* (exact);
     - u = xr − δ, s and P = exp(−u²/2s²); numpy's float128 `exp` is within 4·10⁻²⁰ relative of 40 digits, against double's 1–6·10⁻¹⁷;
     - `synth` (matmul in float128), the residual r = f* − wφ, rP, a, b, and all sums.
   - **In double:**
     - the state (α, δ, σ) and the solver's arithmetic;
     - w = `critical_w(α)` and `curve_mobility(α)`;
     - g and h, rounded on return;
     - the event.
   - These are symmetric or relative roundings that make no asymmetric drift. ċ is Gm·w·h with h computed before rounding, so the drift keeps its long-double accuracy (about 10⁻⁵ relative here).
   - The δ state does not change the problem: it is the same ODE in another coordinate. It only removes c's quantisation at ulp(c*) = 3.47·10⁻¹⁸.
   - The tolerance on δ is loose (atol 3·10⁻¹⁶ ≫ δ early on), but harmless: my runs at atol 10⁻²⁰σ* and 10⁻²²σ* agree to 2·10⁻⁵.
   - **f* is recomputed:** `f = LD(A_A)·exp(−xr²/(2ss²))`, not the stored values (item 4).
   - In double, x − c* is exact for |x − c*| < 0.26. The 271 inexact samples lie beyond 8.7σ*, so X01 is equivalent.
3. **The seed in 40 digits** (seed_dec.py, on `a2_exit_ld`'s rtol-10⁻¹² route).

   | | 3A | 5A |
   |---|---|---|
   | S, exact f*, Simpson 4001 (2001 agrees to 6·10⁻¹⁰) | −7.8847066·10⁻¹⁷ | −5.4945369·10⁻¹⁷ |
   | `a2_seed` against it | +1.9·10⁻⁶ | +1.0·10⁻⁵ |
   | S, stored f* | −7.7621184·10⁻¹⁷ (−1.555%) | −5.4732864·10⁻¹⁷ (−0.387%) |
   | Λ(s_e) against ln(level/\|S\|) | 28.97081 / 28.96729 | 29.33150 / 29.32846 |
   | e^{Λ}\|S\|/level | 1.0035 | 1.0030 |
   | λ changes sign at | s = 3.198·10⁻³ | 6.276·10⁻⁴ |
   | min Λ | −0.619 | −0.215 |
   | stable phase's share | 0.171 | 0.347 |

   - The central difference (10⁻⁹σ*) and the trapezoid are accurate to about 10⁻⁵. A 10⁻⁴ test on the time needs S only to about 3·10⁻³ (λ_e·s_e = 32.6 and 38.3).
   - The stored-f* ratios, 2.338 and 1.943, are round 112's numbers.
   - The 50-digit zero of ċ at the held states (kz_exact.py): 9.5705 and 8.1183 ulps with the stored f*, 9.6426 and 8.1243 with the exact one. The consistent ratios are therefore 2.3568/1.9494 (exact) or 2.3377/1.9432 (stored); the printed 2.37 mixes the two.
4. **The data choice** (ld_repo.py dblf, misc.py cross, own_ld.py).

   | | A, LD f* exact (the md) | A, LD f* stored | symmetric, f* stored | S, stored | cross-planted |
   |---|---|---|---|---|---|
   | 3A | 0.0482742 / 0.0482736 | **0.0482972 / 0.0482973** | 0.0608993 (−) / 0.0642705 (−) | −7.762·10⁻¹⁷ (2.338) | off A's time by 4.87·10⁻⁴ |
   | 5A | 0.00348169 / 0.00348169 | **0.00348205 / 0.00348205** | 0.00451514 (−) / 0.00456995 (−) | −5.473·10⁻¹⁷ (1.943) | off by 9.95·10⁻⁵ |

   - My own long-double right-hand side gives the same pairs:
     - A, exact f*: 0.0482736–0.0482741 and 0.00348169–0.00348170;
     - A, stored f*: 0.0482972–0.0482975 and 0.00348204–0.00348206.
   - With the stored f*, the md would print no flip and "13–26% later".
5. **What the planted test shows** (misc.py neg) [own].
   - On A's samples it confirms the linear response to the drift at c*: e^{Λ}|S| = 1.003·level.
   - In the symmetric problem's own (arithmetic) case:
     - the drift at δ = 0 integrates to S_sym = −2.5·10⁻²³ (3A) and −2.4·10⁻²³ (5A);
     - that is ×260 and ×68 below its effective seed e^{−Λ(s_e)}·0.01σ* (6.5·10⁻²¹, 1.6·10⁻²¹);
     - so the arithmetic's seed arises off the line, and planting S_sym moves the exit by 1.4·10⁻⁴ and 7.8·10⁻³.
   - Had A's seed been the arithmetic's, S would miss the exit by orders of magnitude, and the planted symmetric run would leave at its own noise time.
   - The residual 1–6·10⁻⁶ matches the symmetric run's own long-double seed, about 10⁻⁴ of S: ln(1 + 10⁻⁴)/32.6 ≈ 3·10⁻⁶.
6. **The arithmetic's lateness** (late.py) [own]. Λ along the symmetric route between the double and the long-double exits is 8.405 and 7.360 (3A, rtol 10⁻¹⁰ and 10⁻¹²) and 8.785 and 6.313 (5A), against 11 ln 2 = 7.625.

   | variant | 3A later by | 5A later by | direction (3A) |
   |---|---|---|---|
   | md | 24.9%, 21.5% | 20.3%, 15.0% | +, − |
   | stored f* | 20.1%, 25.2% | 13.7%, 18.4% | −, − |
   | V01 (σ² in double; u stays float128) | 30.2%, 30.5% | 9.6%, 14.3% | −, + |
   | X05 (absolute atol) | 27.8%, 22.4% | 10.7%, 12.9% | −, + |
   | my implementation | 21–25% | 10–19% | mixed |
   | round 112 | 25–34% | 10–15% | mixed |

7. **BDF** (misc.py bdf) [own]. 3A is held to s = 1 at rtol 10⁻¹⁰ and leaves at 0.0515552 (−) at 10⁻¹². 5A leaves at 0.0039016 (−) and 0.0041300 (−).
8. **Mirror** (mirror_chk.py; a 2·10⁻⁵ grid plus 0.01σ grids over ±8σ of every kernel, then bounded maximisation) [own].
   - **H σ* 0.06, equal:**
     - The kernels with σ ≥ h are 11: the 9 formed ones, plus idle kernels at σ 0.0410 (|w| 1.1·10⁻¹¹) and σ 0.1183 (|w| 1.28·10⁻⁸).
     - The widest narrower kernel is 1.01·10⁻³.
     - The 11 are odd by 1.1935·10⁻⁹ (maximum at x 0.357) and 1.1802·10⁻⁹ on the samples. The set is the same for σ ≥ 2h and σ ≥ h/2.
     - Formed at 1%: 5.4589·10⁻⁹. At 0.95%: 4.9015·10⁻³, at x 0.247133.
     - The largest unformed |w| is 0.009803 (1.97% under 0.01); the smallest formed is 0.060378.
   - **σ* 0.04:** 2.2·10⁻¹⁶.
9. **D59 against the md.** Every new number in (i), Remark D‴.59b, (ii)(b) and A2E is verbatim the md's, or round 112's (BDF at 10⁻¹², "several per cent"), and reproduced [own].
10. **Mutants.**
    - **The repository's spec** (W/mut/r112, W/mut/r112k): the control passes. V01–V08, V10 and V11 are killed; V09 and V12 survive (equivalent), as the log says.
      - With the md byte test (`test_round87_split_mass_bdf_and_the_md`) deselected, V01 and V10 survive: they change only noise digits.
      - V01 changes A's long-double exit by 1.4·10⁻⁵ and the symmetric exits to 0.0660 (−) / 0.0670 (+).
    - **My spec116.py (21 plus control).**
      - Killed by the semantic tests (rerun with the md byte test deselected): X02 (the stored f*), X03 (the c-correlation summed in double), X08 (the drift at the run's own offset), X09 (planted on the other side).
      - Killed by the planted tests: X13 (floor → round).
      - Killed in the full run, by the md byte test (no deselected rerun): X12, X14, X15, X16, X20, X21.
      - Killed only by md noise digits: X05 (absolute atol) and X07 (2001 points: 5A's "1.8e-06" becomes 2.0e-06).
      - Survivors:
        - X01 (x − c* in double: exact within ±0.26);
        - X04 (g in double: α only);
        - X06 (λ step 10⁻⁶σ*);
        - **X10 (the long-double exits' direction unchecked: no planted input)**;
        - X11 (the seeded run against the 10⁻¹⁰ run: 1.2·10⁻⁵ < 10⁻⁴);
        - X17 (strict width test);
        - X18 (all kernels: every w ≠ 0);
        - X19 (the margin over min A: equal pairs).
      - All survivors are equivalent on the data. None changes a published result.
11. **Literature.** Four web searches; no page fetched.

## Verdicts per part

| claim (D59 unless noted) | verdict | evidence |
|---|---|---|
| (i) one clause, "which the problem as posed perturbs: Remark D‴.59b" | **holds**; the move is clean | re-read |
| Remark: "in double the centre sample lies 8 ulps …; 8–10 ulps above c*" | **holds** (9.57/8.12 stored f*, 9.64/8.12 exact) | [own] |
| Remark: "(LSODA, and scipy's BDF at rtol 10⁻¹⁰, hold 3A; at 10⁻¹² BDF leaves)" | **holds** (0.05156, −); not in the md | [own] |
| Remark: "its seed is the data's" | **holds, decided**: 40-digit S equals the long-double S; e^{Λ}\|S\|/level = 1.003; two implementations agree | [own] |
| "in long double the exit is the same at two tolerances (0.0482736, 0.00348169)" | **holds** (1.2·10⁻⁵), **for the exact-f* problem**; the posed problem gives 0.0482973, 0.00348205 (D3dk-1) | [own] |
| route seed S = ∫e^{−Λ}ċ\|c* ds; −7.885·10⁻¹⁷, −5.495·10⁻¹⁷ | **right**, accurate to 2·10⁻⁶ and 10⁻⁵ (central difference and trapezoid adequate) | [own] |
| "2.37 and 1.95 times the held offset" | **mixed data** (2.36/1.95 exact, 2.34/1.94 stored); the denominator is the zero of ċ, not the held offset (D3dk-1, -3) | [own] |
| "planted … gives A's exit to 1.2·10⁻⁶ and 1.8·10⁻⁶" | **holds to 10⁻⁵**; the digits are the integrator's floor (1.2–6.4·10⁻⁶) (D3dk-2) | [own] |
| planting S at s = 0 as the test | **right**: it measures the linear response; the negative control fails it by orders of magnitude in seed | [own] |
| "exact sums keep c at the centre, reversing the sums mirrors the run" | **true, stated as identities** | md, round 112 |
| "the tolerance can flip it" | **holds** (flips occur in every variant) | [own] |
| "in long double it comes 15–25% later" | **qualitatively right** (ΔΛ ≈ 11 ln 2); **the range is one realisation** (9.6–30.5%) (D3dk-2) | [own] |
| md exit line's `data`/`arith`/`flip`/`late` conditions | **print correctly on the data**; `data` decided in long double as stated | [own], mutants |
| `_LDBumps` "only the arithmetic changes" | **sums right; f* changed in `a2_exit_ld`** (D3dk-1) | code, [own] |
| (ii)(b) "the 11 kernels at least as wide … odd by 1.2·10⁻⁹ … (the narrowest 0.041, the widest narrower one 0.001)" | **holds** (1.19·10⁻⁹; the gap ×40) | [own] |
| (ii)(b) "the formed kernels' 5.5·10⁻⁹ rests on the formed threshold's margin, 0.0098 against 0.06" | **number right; the margin is 2% against 0.01** (D3dk-5) | [own] |
| md mirror lead and `mirror_asym` docstring | **fixed** | [md] |
| Fenichel "normally hyperbolic, strictly attracting or repelling in c" | **loose**: normal hyperbolicity also needs domination over the tangential rates (D3dk-6) | [snippet] |
| tests | **strong on semantics**; V01 and V10 are killed only by md noise digits; X10 is unplanted (D3dk-7) | [own] |

## Status of the round-112 issues

| D3dg | status | what remains |
|---|---|---|
| 1 (the exit's tests) | **applied**: `_LDBumps`, `a2_exit_ld` and `a2_seed` added; the data clause is decided in long double; the identities are stated; the flip and lateness are computed | f* is recomputed (D3dk-1). Noise digits are printed (D3dk-2). |
| 2 (the formed threshold) | **applied**: the width split and the margin are printed; D59 leads with the width split | The margin's phrasing (D3dk-5). |
| 3 (D59 wording) | **applied** (the BDF clause, "the arithmetic's own scatter", "follows the sums' order" replaced) | The BDF 10⁻¹² fact is cited but not in the md (D3dk-6). |
| 4 (stale or loose text) | **applied** (lead and docstring, `edge`, Fenichel, Remark D‴.59b) | The Remark sits inside the proof; the Status omits it; the normal-hyperbolicity wording (D3dk-6). |
| 5 (tests) | **applied** (planted inputs T01, T03 and the long-double ones; the new test; 12 mutants) | V01 is not in double; V01 and V10 are killed only by md digits; X10 is unplanted (D3dk-7). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dk-1 | minor | **The long-double run changes the data.** `a2_exit_ld` sets `f = LD(A_A)·exp(−xr²/(2ss²))`; the docstring says "the data are the double samples; only the arithmetic changes". The double runs, the raw and `a2_ulps_exact` use the stored (double-rounded) f*. <br>• f*'s rounding is −1.555% (3A) and −0.387% (5A) of S (40 digits) [own]. <br>• Posed-problem long-double exits are 0.0482972/0.0482973 and 0.00348205, against the printed 0.0482742/0.0482736 and 0.00348169. S is −7.762·10⁻¹⁷ and −5.473·10⁻¹⁷ (round 112's numbers) [own]. <br>• "2.37" = S(exact f*)/kz(stored f*); consistent: 2.36/1.95 (exact f*, kz 9.64/8.12) or 2.34/1.94 (stored) [own]. <br>• The planted test resolves the choice: cross-planting gives 4.9·10⁻⁴ (3A) [own]. The conclusion is unaffected. | **Code (`a2_exit_ld`):** use the stored f* (`feature_sum(...)` or `(A_A·np.exp(−(xs − c0)²/(2ss²))).astype(LD)`), the posed problem and the same data as `a2_ulps_exact`. Or keep the exact Gaussian, say so in the docstrings, and compute kz on it. **md / D59 (Remark, A2E):** the numbers follow (0.0482973, 0.00348205; S −7.762·10⁻¹⁷, −5.473·10⁻¹⁷; 2.34, 1.94). Optionally add "f*'s double rounding is 1.6% and 0.4% of the seed". **Test:** tighten the window 2.3–2.45 to the chosen data. |
| D3dk-2 | minor | **Noise digits published as results.** <br>• "A's time to 1.2·10⁻⁶ and 1.8·10⁻⁶": 1.21, 6.39, 1.22, 2.36·10⁻⁶ for n = 2001 … 16001, while S moves by ≤ 7·10⁻⁶ (a time shift ≤ 2·10⁻⁷). X07 turns 5A's figure into 2.0·10⁻⁶, X05 into 1.9 and 6.2·10⁻⁶ [own]. <br>• The symmetric samples' six-digit exits, "15–25% later" and "(3A in long double)" are noise. V01 (σ² in double) gives 30.2/30.5% and 9.6/14.3%; X05 gives 10–28% and flips both directions; the stored f* gives 13–26% and no flip; my implementation 10–25% [own]. <br>• Noise-free equivalents: ΔΛ between the double and long-double symmetric exits is 6.3–8.8 against 11 ln 2 = 7.62; e^{Λ(s_e)}\|S\|/level = 1.0035 and 1.0030 [own]. | **md (`a2_lines`):** print the seeded run as "within 10⁻⁵ of A's time (the integrator's resolution)". Have `a2_seed` return Λ(s_e) and print it against ln(level/\|S\|). For the symmetric samples print ΔΛ between the double and long-double exits against 11 ln 2 ("the seed scales with the unit roundoff"), and the times to 3 digits. **D59 (Remark, A2E):** "to 10⁻⁵", "later by Λ ≈ 11 ln 2 (15–25% here)", "the tolerance can flip it". |
| D3dk-3 | cosmetic | **"2.37 and 1.95 times the held offset."** The denominator is −kz·ulp, c*'s distance below the zero of ċ (9.57 and 8.12 ulps). The held offsets are 9.00 and 8.00 ulps above c* (ratios 2.53, 1.98), of opposite sign. S is printed in units of c, the held offsets in σ* (S = −2.63·10⁻¹⁵σ*) [own]. | md and D59: "2.37 and 1.95 times c*'s distance below the zero of ċ (9.57, 8.12 ulps)"; "S = −7.885·10⁻¹⁷ (−2.63·10⁻¹⁵σ*)". |
| D3dk-4 | cosmetic | **The Remark's decisive step is implicit.** It rests "the data's" on the tolerance-independence (not decisive alone: the symmetric run's time also moves by only 1.5%) and on the planted run. The precision-invariance of A's exit is not stated: double 0.0483/0.0478, long double 0.04827 (≤ 1%), against the symmetric problem's ΔΛ ≈ 11 ln 2. Nothing in the repository shows that the long-double drift is the data's; my 40-digit S does (2·10⁻⁶) [own]. | D59 Remark: "from double to long double A's exit moves by ≤ 1%, the symmetric samples' by Λ ≈ 11 ln 2; the drift at c* in 40 digits equals the long-double one". Optionally add a decimal check of d(s) at a few route points to `a2_seed`. |
| D3dk-5 | cosmetic | **(ii)(b)'s margin and the "11 kernels".** The formed set's 5.5·10⁻⁹ rests on 0.0098 against the **threshold 0.01** (2%), not against the smallest formed 0.06. The 11 wide kernels are the 9 formed ones plus 2 idle ones (\|w\| 1.1·10⁻¹¹ at σ 0.041, 1.3·10⁻⁸ at σ 0.118); the "narrowest 0.041" is idle [own]. | D59 (ii)(b) and the md: "the largest unformed \|w\| 0.0098, 2% under the threshold 0.01·A"; "(the 9 formed and 2 idle ones; the narrowest formed 0.050)". |
| D3dk-6 | cosmetic | **Structure and wording.** <br>• Remark D‴.59b sits inside the Proof, and (ii)(a)'s proof continues in the same paragraph. <br>• The Status does not name the Remark, whose content is numerical. <br>• "At 10⁻¹² BDF leaves" is a referee's number, not the md's (verified: 0.05156, −) [own]. <br>• "Normally hyperbolic, strictly attracting or repelling in c": normal hyperbolicity also requires the normal rate to dominate the tangential ones, and it fails near λ's sign change (s ≈ 3.2·10⁻³ / 6.3·10⁻⁴) [snippet]. | Put the Remark after the proof's ∎, or start (ii)(a) on a new paragraph. Status: "Remark D‴.59b numerical (in long double; the route seed)". Add BDF at rtol 10⁻¹² to `a2_exits`' `imp`. Fenichel: "normally hyperbolic (its normal rate nonzero and dominating the tangential ones)". |
| D3dk-7 | cosmetic | **Tests and mutants.** <br>• V01 "_LDBumps in double" keeps u in float128 (numpy promotion); only σ² is in double. It and V10 survive once the md byte test is deselected: killed by noise digits [own]. <br>• X10 (the long-double exits' direction dropped from `data`) survives: no planted input. <br>• The 10⁻⁵ window on the seeded run has a margin under 2× (6.4·10⁻⁶ at n = 4001). | Replace V01 by X03 (the c-correlation summed in double, killed by `test_round112`) or by an xr-and-sums-in-double mutant. Plant an upward long-double exit (`ld.a = [(0.05, 1), (0.05, 1)]` → no "data's" clause). Widen to 3·10⁻⁵ or print as in D3dk-2. Record mutants killed only via md noise digits as such. |

## Sources

- [snippet] NumPy, "Data types" ([numpy.org](https://numpy.org/doc/stable/user/basics.types.html)): `longdouble` is the C compiler's long double — 80-bit on most x86, padded to 96 or 128 bits (hence `float128`); the extra precision is easily lost through Python floats; check with `np.finfo(np.longdouble)`. This fits eps 1.084·10⁻¹⁹ here.
- [snippet] Wikipedia, [Extended precision](https://en.wikipedia.org/wiki/Extended_precision) and [Long double](https://en.wikipedia.org/wiki/Long_double): the x87 80-bit format has a 64-bit significand.
- [snippet] "Breakdown of normal hyperbolicity … generalized Lyapunov-type numbers" ([arXiv:0912.5419](https://arxiv.org/pdf/0912.5419)), and lecture notes on invariant manifolds ([Utrecht](https://webspace.science.uu.nl/~verhu101/InvManSPT2004.pdf)): normal hyperbolicity requires the normal expansion or contraction to dominate the tangential flow (Lyapunov-type numbers below critical values). This is the basis of D3dk-6.
- [snippet] Forest, Goedde and Sinha, Math. Comput. Simul. 37 (1994) ([RePEc](https://ideas.repec.org/a/eee/matcom/v37y1994i4p323-339.html)): roundoff of the order of machine arithmetic can break a symmetry of the data and equations and unleash an instability absent in the symmetric subspace.
- [memory] The linear escape time (1/λ)·ln(level/seed), behind ΔΛ ≈ ln(ε_double/ε_ld) = 11 ln 2. Also Fenichel 1971 and Golubitsky–Stewart–Schaeffer, as cited by D59 from earlier rounds; not re-searched.

## Lesson

Raise the precision and keep every stored input bit for bit.
- `a2_exit_ld` changed the arithmetic, but it also re-evaluated the target. That moved the seed by 1.6%, and the planted test meant to decide the question resolves 1.6%.
- A precision test separates the data's seed from the arithmetic's only if the data are the stored ones. Every comparison quoted next to it (the double runs, the 50-digit zero) must use the same inputs.

Print an instrument at the resolution it has.
- A seeded time agrees to 10⁻⁵, not to the 1.2·10⁻⁶ the integrator happened to give.
- A noise-seeded exit is characterised by the gap in Λ (≈ 11 ln 2), not by its six-digit times or a single "15–25%".
