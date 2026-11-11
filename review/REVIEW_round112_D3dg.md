# Referee report D3dg, round 112: Proposition D‴.59 as revised after round 109

## Scope

- **Export and work directory.** I exported commit cac7a71 (branch claude/vigilant-clarke-hf939b) with `git -C … archive cac7a71 | tar -x -C W/repo`.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112.
  - W/repo/venv is a symlink to the repository's venv: Python 3.13.18, numpy 2.5.4, scipy 1.18.1.
  - The machine's long double is x87 80-bit (64-bit mantissa, eps 1.08·10⁻¹⁹).
- **Diff 319e8a5 → cac7a71 (the D59 files).**
  - The raw is not in the diff; every new number comes from `--from-raw`.
  - exp_moving_width.py gains `_FsumBumps`, `a2_exit(summ=…)`, `a2_exits`' `summ` block, the bisected `a2_lemma_signs`, `mirror_odd`'s kernel grid, the halved bottleneck in `split_summary`, and the exit line's `data`/`arith` conditions. It also takes kz from `a2_ulps_exact`, uses one "held" criterion, and changes the a2e_job docstring.
  - The md changes at lines 54, 61–63 and 621.
  - tests/test_moving_width.py extends two tests (still 54 tests). tools/mutants/round109_D59_own.py is new.
  - D59.md changes in (i), the lemma's edges, (ii)(b) and the A2E numerics. Its README row changes too.
- **Read in full:**
  - D59.md (409 lines);
  - every changed part of exp_moving_width.py and the code it calls: `a2e_job`, `a2_state`, `a2_fixed_point`, `a2_lemma`, `a2_ulps`, `a2_ulps_exact`, `a2_lemma_scan`, `a2_lines`, `split_summary`'s mirror block, `BumpsCW` (`_parts`, `corr`, `synth`), `critical_rhs_nd`/`critical_jac_nd`, `critical_w`, `curve_mobility`, `_alpha_of`, `coef_crit`;
  - results/chaos/moving_width.md, tests/test_moving_width.py, tools/mutants.py, tools/mutants/round109_D59_own.py;
  - the Round 109 entry of review/REVIEW_LOG.md, review/REVIEW_round109_D3dd.md, and the git diff of these files.
- **Settings.**
  - Every process ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` and `timeout`, with MW_OUT = W/out (a copy of the raw) or a mutant slot's `out`.
  - pytest ran with `--basetemp` inside W.
  - My copy of the mutant tool, /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112/own/mutants_tool.py, sets REPO to W/repo and adds `--basetemp` per slot.
  - At most two heavy processes of mine ran at a time.
- **Process and file notes.**
  - I stopped four processes of my own: two slow `seed.py` runs, twice (PIDs 5108/5109, then 5658/5659). I checked each one via /proc/PID/cmdline and stopped it with `kill PID`.
  - The other sessions' processes (D60 mutants, width_mismatch tests) were not touched.
  - One rerun used `python -I`, which ignores PYTHONDONTWRITEBYTECODE. It wrote `__pycache__` into W/repo, my export only, and I deleted it before the later runs.
- **The repository is untouched.** I used only read-only git (archive, diff, log, status, rev-parse).
  - Meanwhile HEAD moved to dafe1d0 (rounds 107, 110, 111 reports; part VC), and other sessions hold uncommitted D58/D60 work.
  - `git diff cac7a71` on D59.md, the script, the md, the raw, the tests, the mutant spec and src/core/moving_width.py is empty at the end.
- **Tags.** [md] the committed md or draft; [raw] the committed raw; [own] my computation; [snippet] a web-search result I saw (no page fetched); [memory] not verified this round.

## Summary

**Design.**
1. I regenerated the md from the raw, ran the 54 tests, and reran 11 jobs against the raw bit for bit.
2. I reproduced `a2_exits` at full precision.
   - I then wrote my own independent right-hand side and varied the arithmetic with the data kept: numpy sums, `math.fsum`, reversed order, an offset (δ = c − c*) state in double, long double with a double c, and long double with a δ state. I ran it on A's samples and on exactly symmetric ones, at two tolerances.
   - With the repository's own right-hand side I also varied the tolerance (10⁻⁹ … 10⁻¹³) and the start's α by ±1, ±2 ulps.
3. I integrated the data's seed along the route in long double: d(s) = ċ at c = c*, λ(s) in closed form, and S = ∫e^{−Λ}d. I then started the symmetric problem at S.
4. Lemma: my own implementation of the terms, a scan of step 10⁻⁴ over σ/h ∈ [0.03, 1.2], and Brent roots.
5. Mirror:
   - my own resolved odd part (a 10⁻⁵ grid, ±8σ at 0.01σ around every kernel, then bounded maximisation);
   - actual sample pairs;
   - a scan of the formed threshold;
   - brute force over all involutions (and over all 9! permutations for H σ* 0.06).
6. Mutants: the repository's 10 plus a control, and 16 of my own plus a control.
7. Web searches.

**Verdicts.**
- **Reproducibility is exact.** The md matches byte for byte (sha256 464dd02d…). 54 tests pass. All 11 rerun jobs match the raw bit for bit.
- **Every number the revision added reproduces** [own].
- **`_FsumBumps` changes only the summation.** The products, the residual, `synth`, the 1/N and the shapes are untouched. "Exactly rounded" is right for the sums of the double-rounded products.
- **"On A's samples the exit's seed is the data's" holds**, and more sharply than printed:
  - In long double with an offset state, A's exit is 0.0482973 (3A) and 0.0034820 (5A). It agrees to 10⁻⁶ across tolerances.
  - The data's drift integrated along the route is −7.762·10⁻¹⁷ (3A) and −5.473·10⁻¹⁷ (5A), 2.338 and 1.943 times the held offset.
  - The symmetric problem started at that seed leaves at A's long-double time to 5·10⁻⁶ [own].
- **But the printed conditions are weak tests (D3dg-1).**
  - The arithmetic's two conditions are identities on exactly symmetric data: exact sums keep c at the centre exactly, and reversed sums mirror the run bit for bit.
  - The symmetric problem's direction also changes with the tolerance and the start's last bits.
  - The offset start's 5% window lies inside the double arithmetic's scatter: my independent double implementation would print no "data's" clause.
- **`mirror_odd` resolves every kernel on the data**, and `mirror_asym`/2 is the sup-norm distance (0.119964 by brute force).
- **(ii)(b)'s "the formed kernels' function odd by 5.5·10⁻⁹ … between the samples" rests on a 2% margin at the formed threshold (D3dg-2).** A threshold-free version holds (1.2·10⁻⁹).
- **The lemma's edges are the roots**, with exactly one sign change per predicate on [0.03, 1.2].
- **Tests are strong.** All 10 repository mutants are killed. 8 of my 16 are killed, and the 8 survivors are equivalent on the data.
- **The story of (i) is now consistent** from statement to lemma to perturbation to numerics. Its quantitative mechanism (the route seed) is still missing, and three phrases overstate (D3dg-3).
- **No major issue.**

**Findings.**
1. **D3dg-1 (minor): the exit's printed tests.**
   - "No exit by s = 1 with exact sums" and "the reversed run flips" cannot fail on exactly symmetric samples. The fsum run keeps c = 0.0 exactly over all 490 / 2337 steps. The reversed run's times equal the forward ones bit for bit.
   - "Its direction follows the sums' order" is incomplete. In the code's own order, rtol 10⁻¹³ (3A) or 10⁻¹¹ (5A), or α moved by −2 ulps (5A), flips the symmetric problem's exit; A's never flips.
   - The offset start's 5% window lies within double scatter. From −9.57 ulps the symmetric exit is 0.04899–0.05195 across double variants (up to +7.6% of A's), and 0.050688 (+4.99%) with the repository's own `_FsumBumps`.
   - The decisive tests (long double; the route seed) are not printed.
2. **D3dg-2 (minor): the formed threshold in (ii)(b).**
   - The largest kernel below the threshold is the invisible −0.0098 spike (σ 3.25·10⁻⁴, φ ≤ 1.3·10⁻¹⁷ on the samples), 2% under 0.01·A.
   - With a threshold of 0.0095·A or lower, "the formed kernels' function" is odd by 4.90·10⁻³. My mutant T14 (0.9%) changes the published 5.5e-09.
   - The 11 kernels at least as wide as the spacing give 1.19·10⁻⁹ on and between the samples, with a ×40 gap in width.
3. **D3dg-3 (cosmetic): D59 wording.**
   - "At the same time to 2% … as does the symmetric problem started at the posed problem's offset": that start is +2.6% at 3A (the md says 5%).
   - "(LSODA and scipy's BDF hold 3A)": scipy's BDF holds at rtol 10⁻¹⁰ (to s = 10 too) but leaves at 10⁻¹² (0.0516, downward).
   - "Its direction follows the sums' order" (see D3dg-1).
4. **D3dg-4 (cosmetic): stale or loose text.**
   - md:621 and the code at 1692 introduce the mirror block as "the formed kernels' distance from the reflection's fixed-point subspace", while each entry says "from the mirror image (… from the subspace)". `mirror_asym`'s docstring opens the same way.
   - `edge(pred, k)`'s k is unused.
   - The Fenichel clause restricts to pieces "where it is strictly repelling".
   - Statement (i) carries a 12-line numerical parenthetical.
5. **D3dg-5 (cosmetic): tests.**
   - Survivors T01 (the offset start's time unchecked) and T03 (A's fsum/rev directions dropped) show two printed conditions with no planted input.
   - The arithmetic's plants exercise identities.

## What I checked and how

Scripts: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112/own/rerun.py, …/own/exits.py, …/own/repo_exits.py, …/own/symdir.py, …/own/adir.py, …/own/seed2.py, …/own/lemma.py, …/own/mirror.py, …/own/invol.py, …/own/spec112.py; logs in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112/logs/ (…/own/ abbreviates /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112/own/).

1. **Reproducibility** [own].
   - `--from-raw` took 11.2 s; `cmp` against the committed md: identical.
   - pytest: 54 passed in 61.2 s; the md in W/out is unchanged afterwards.
   - rerun.py compares every field except `secs`. Identical:
     - `part_A2E` (3 entries);
     - `a2_job` (0.25, 3A), (0.25, 5A), (1, 10A), (4, 20A);
     - `h_job` H σ* 0.06 equal and σ* 0.04 equal (the mirror inputs);
     - H3 ±0.15 and ±0.12.
2. **`_FsumBumps`** [own reading and runs].
   - `corr` builds the same `rP = r[:, None]·P` and the same factors `a = u/s²`, `b = u²/s³`, and the same elementwise products. Only `.sum(0)` becomes `math.fsum` per column, before the same `/N`. The shapes are (d,) and (d, 2).
   - `synth` for one kernel is a single product. `second` is never requested (DOP853 needs no Jacobian).
   - On exactly symmetric samples the products are exactly antisymmetric, so fsum returns 0.0 and c stays 0.0 at every step: max |c| = 0 over 490 (3A) and 2337 (5A) steps to s = 1.
3. **The exits, the repository's code at full precision** (repo_exits.py) [own].

   | | A's samples, DOP853 1e-10 / 1e-12 | fsum / rev (A) | symmetric 1e-10 / 1e-12 | fsum / rev (sym) | offset start −kz (50-digit) |
   |---|---|---|---|---|---|
   | 3A | 0.048277 (−) / 0.047787 (−) | 0.048142 (−0.28%) / 0.048040 (−0.49%), (−) | 0.050716 (+) / 0.051354 (+) | none / 0.050716 (−) | 0.049528 (−), +2.59% |
   | 5A | 0.0034833 (−) / 0.0034883 (−) | 0.0034579 (−0.73%) / 0.0034812 (−0.06%), (−) | 0.0039725 (+) / 0.0038608 (+) | none / 0.0039725 (−) | 0.0035200 (−), +1.05% |

   - Radau is +17.4% (3A) and +11.9% (5A); BDF holds 3A to s = 1 and leaves 5A at 0.0039016.
   - The reversed symmetric run's times equal the forward run's bit for bit (0.050716099718251716; also at 1e-12).
   - Extra starts, all with the code's sums except the last:
     - from the double zero (9.93 ulps): 0.048925;
     - from +kz: 0.048658 (+);
     - with `_FsumBumps` from −kz: 0.050688 (−), +4.99% of A's;
     - 5A, `_FsumBumps` from −kz: 0.0035349.
4. **Own right-hand side, arithmetic varied** (exits.py; DOP853, rtol 1e-10 / 1e-12) [own].

   | arithmetic | A 3A | A 5A | sym 3A | sym 5A |
   |---|---|---|---|---|
   | numpy (1-D) | 0.047871 / 0.048634 (−,−) | 0.003507 / 0.003488 (−,−) | 0.051860 (−) / 0.054681 (+) | 0.003930 (+) / 0.003669 (−) |
   | fsum | 0.048590 / 0.048367 (−,−) | 0.003484 / 0.003490 | none / none | none / none |
   | reversed | 0.048222 / 0.048696 (−,−) | 0.003484 / 0.003488 | mirror of numpy | mirror of numpy |
   | δ-state, double | 0.048549 / 0.047760 | 0.003466 / 0.003488 | = numpy | = numpy |
   | long double, double c | 0.048236 / 0.048401 | 0.003481 / 0.003488 | — | — |
   | **long double, δ-state** | **0.0482973 / 0.0482973** | **0.0034820 / 0.0034821** | **0.065077 (−) / 0.068778 (+)** | **0.004362 (+) / 0.004430 (−)** |

   - On A's samples the double runs scatter ±1% around the long-double value. Their tolerance dependence comes from the double representation of the offset: in long double the two tolerances agree to 10⁻⁶.
   - On symmetric samples long double delays the exit by 25–34% (3A) and 10–15% (5A), and the direction changes with the tolerance.
5. **Sensitivity with the repository's own right-hand side** (symdir.py, adir.py; DOP853 rtol 1e-10 unless stated) [own].
   - Symmetric, 3A:
     - tolerance 1e-9 … 1e-13: 0.04945 (+), 0.05072 (+), 0.05152 (+), 0.05135 (+), **0.05287 (−)**;
     - α at −2, −1, +1, +2 ulps: 0.05777, 0.05523, 0.05358, 0.05162 (all +).
   - Symmetric, 5A:
     - tolerance: 0.004186 (+), 0.003972 (+), **0.003725 (−)**, 0.003861 (+), 0.003768 (+);
     - α at −2 ulps **0.003920 (−)**; −1, +1, +2: 0.003691, 0.003916, 0.003674 (+).
   - A's samples under the same nine changes: 3A 0.04779–0.04836, 5A 0.003458–0.003519, every one (−).
6. **The seed along the route** (seed2.py; long double, δ-state, rtol 1e-12, 40001 points) [own].

   | | 5A | 3A |
   |---|---|---|
   | λ crosses 0 at | s = 6.28·10⁻⁴ (σ/σ* 0.286) | s = 3.21·10⁻³ (0.275) |
   | min Λ | −0.2153 | −0.6191 |
   | d changes sign at | 3.47·10⁻⁴, 5.79·10⁻⁴ | 1.06·10⁻³, 2.74·10⁻³ |
   | seed S = ∫e^{−Λ}d | −5.473·10⁻¹⁷ | −7.762·10⁻¹⁷ |
   | against −kz·ulp | −2.817·10⁻¹⁷ (×1.943) | −3.320·10⁻¹⁷ (×2.338) |
   | stable phase's share | 0.346 | 0.165 |
   | e^Λ·S against δ at 0.5 / 0.7 / 0.9 s_e | to 4 digits | to 4 digits |
   | λ at the exit; ln(S/kz·ulp)/λ | 1.101·10⁴; 6.0·10⁻⁵ of s | 675.1; 1.26·10⁻³ of s |
   | **symmetric problem started at S (long double)** | **0.00348205 (both tolerances)** against A's 0.00348204/0.00348206 | **0.0482972 / 0.0482970** against A's 0.0482973 / 0.0482973 |
   | symmetric started at −kz (long double) | 0.0035431 (+1.75%) | 0.049574 (+2.64%) |
   | symmetric started at −kz, my double variants | 0.003524–0.003596 | 0.04899–0.05195 |

   Round 109's figures (×1.94, ×2.34; λ crossing at σ/σ* 0.285, 0.271) are reproduced independently.
7. **The lemma** (lemma.py) [own].
   - My terms equal `a2_lemma_scan`'s to all printed digits: +2871.03, −2816.51, −213.725 at w = 4A, σ/h = 1; +109.608, −316.792, −74.7662 at the least-squares w.
   - On [0.03, 1.2] at step 10⁻⁴ there is exactly one sign change per predicate:
     - w = 4A: the total at 0.843696, the ±1 term at 0.469602, the ±1 term above the centre's nowhere;
     - least-squares w: 0.774625, 0.619540, 0.558315.
   - The repository's bisection gives 0.8436964, 0.4696022, None; 0.7746253, 0.6195403, 0.5583148.
8. **Mirror** (mirror.py, invol.py) [own].

   | end | repo formed / own | repo all, between / own | on the samples = residual's | threshold margin (largest unformed / smallest formed \|w\|) |
   |---|---|---|---|---|
   | H σ* 0.04 | 1.1·10⁻¹⁶ / 1.7·10⁻¹⁶ | 1.7·10⁻¹⁶ / 2.8·10⁻¹⁶ | roundoff | 4.6·10⁻¹² / 0.055 |
   | **H σ* 0.06** | **5.4588·10⁻⁹ / 5.4589·10⁻⁹** | **4.9015·10⁻³ / 4.9015·10⁻³** (x = 0.247133) | 1.1802·10⁻⁹ = 1.1802·10⁻⁹ | **0.009803 / 0.0604** |
   | H3 ±0.05 | 1.5080·10⁻⁴ / 1.5080·10⁻⁴ | 1.2503·10⁻⁶ / 1.2504·10⁻⁶ | 1.2501·10⁻⁶ | 0.0072 / 0.0362 |
   | H3 ±0.12 | 2.643·10⁻¹³ / 2.643·10⁻¹³ | 2.35·10⁻¹⁴ / 2.35·10⁻¹⁴ | 2.35·10⁻¹⁴ | 1.8·10⁻⁷ / 1 |
   | HP σ* 0.04 | 1.397·10⁻⁷ / 1.397·10⁻⁷ | 1.536·10⁻⁷ / 1.536·10⁻⁷ | 1.536·10⁻⁷ | 7.9·10⁻⁵ / 1.0 |
   | HP σ* 0.06 | 8.96·10⁻¹⁰ / 8.96·10⁻¹⁰ | 1.02·10⁻¹² / 1.02·10⁻¹² | 1.02·10⁻¹² | 0.00955 / 0.01115 |

   - **Threshold scan, H σ* 0.06** (the formed kernels' odd part, resolved): 5.46·10⁻⁹ at thresholds of 0.0099–0.05·A; **4.90·10⁻³ at 0.0095·A and below** (10 kernels, the −0.0098 spike included).
   - That spike sits at c −0.24713 with σ 3.25·10⁻⁴; its largest φ on the samples is 1.3·10⁻¹⁷, so it is invisible.
   - The 11 kernels with σ ≥ h give 1.19·10⁻⁹. The 14 visible ones give 1.75·10⁻³ (narrow visible spikes, |w| ≤ 0.0035).
   - **HP σ* 0.06:** 2.1·10⁻¹¹ at ≤ 0.0095·A and 2.6·10⁻⁹ at 0.02·A (immaterial).
   - **Bottleneck:** brute force over involutions equals `mirror_asym` on all nine formed sets (2.16·10⁻¹¹, 0.24, 0.00165, 6.94·10⁻¹⁸, 5.15·10⁻¹³, 0.3, 1.52·10⁻⁶, 1.43·10⁻⁷, 2.85·10⁻⁸). Over all 9! permutations, H σ* 0.06 gives 0.239927, so the sup-norm distance is 0.119964.
9. **Numbers** [own].
   - Every number added to the md and D59 since 319e8a5 is reproduced:
     - the exits: 0.0481, 0.0480, 0.0035 ×2, "no exit by s = 1", 0.0507 (−), 0.0040 (−), 9.57 and 8.12 ulps, 0.0495, 0.0035;
     - the lemma: 0.8437, 0.4696, 0.7746, 0.6195, 0.5583;
     - the mirror block, all nine ends: from the mirror image, from the subspace, all kernels, formed odd, on the samples, between.
   - Also reproduced: the 0.75-ulp test margins, which are 0.57 (3A) and 0.12 (5A).
10. **Mutants** (res files in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref112/mut/).
    - **The repository's spec:** 10 mutants, all killed; the control passes. S04, S05, S06 and S10 are killed by the tests alone; the other six by the tests and the md.
    - **spec112.py (16 + control): 8 killed.**
      - By the tests only: T02 (2% → 25%) and T10 (10 bisection steps, through the ±10⁻⁶ sign probe).
      - By the tests and the md: T07 (numpy for the σ-correlation), T08 (numpy for the c-correlation), T11 (bisection reversed), T13 (formed odd part on the samples only), **T14 (formed threshold 0.9%: the published 5.5e-09 becomes 4.9e-03)**, T16 (only the symmetric problem reversed).
    - **8 survivors, all equivalent on the data:**
      - **T01: the offset start's time unchecked**;
      - **T03: A's fsum/rev directions dropped**;
      - T04: implicit exits dropped;
      - T05: the centre run's exit unchecked (implied by the flip check through the mirror identity);
      - T06: numpy for g;
      - T09: odd-part grid at the centres only;
      - T12: held at 10⁻⁶;
      - T15: reference time a[1].
    - None changes a published result.
11. **Literature.** Six web searches; no page fetched. Sources lists only what I saw.

## Verdicts per part

| claim (D59.md unless noted) | verdict | evidence |
|---|---|---|
| (i) continuum landscape, convergence, Lemma D‴.59a's invariance and decoupling | **hold**; unchanged | re-read |
| (i) "the samples' roundoff — the centre sample's 8 ulps, partly offset by the other pairs' — drifts c" | **holds** (round 109's split) | [own] round 109, re-read |
| (i) "a repelling near-symmetric state … 8–10 ulps above c* at A2E's states" | **holds** (8.118, 9.570) | [own] |
| (i) "at s = 0.0483 (3A) and 0.0035 (5A) by DOP853, downward with every integrator that leaves" | **holds** | [own] item 3 |
| (i) "(LSODA and scipy's BDF hold 3A)"; A2E "its BDF (which holds 3A)" | **tolerance-dependent**: BDF holds at rtol 10⁻¹⁰ (to s = 10 too), leaves at 10⁻¹² (0.0516, −) (D3dg-3) | [own] |
| (i) "at the same time to 2% with the samples' sums exactly rounded … or reversed" | **holds** (≤ 1.0%) | [own] items 3–5 |
| (i) "… as does the symmetric problem started at the posed problem's offset" | **overstated**: +2.6% at 3A (the md says 5%); in double the time is arithmetic-sensitive (+1.5 … +7.6%) (D3dg-1, -3) | [own] item 6 |
| (i) "the exit's seed the data's" | **holds, sharply**: the long-double exit is tolerance-free; the symmetric problem seeded with the route seed leaves at A's time to 5·10⁻⁶ | [own] items 4, 6 |
| (i) "on exactly symmetric samples from their centre its direction follows the sums' order" | **incomplete**: reversal mirrors exactly, but the tolerance or the start's last bits also flip it (D3dg-1) | [own] item 5 |
| (i) "with exactly rounded sums there is no exit by s = 1: the arithmetic's" | **true, an identity** (c ≡ 0.0) | [own] item 2 |
| md "the samples' sums exactly rounded" | **right** for the sums of the rounded products | [own], [snippet] |
| md exit line's printed conditions (2%, 5%, exact-sum and reversed symmetric runs) | **print correctly on the data**; the arithmetic's two conditions are identities; the 5% window is within double scatter (D3dg-1) | [own] |
| Lemma's edges "each a root of its term" | **right**: one sign change per predicate on [0.03, 1.2] | [own] item 7 |
| `a2_lemma_signs`' bisection | **right**; `k` unused | code, [own] |
| `mirror_asym` "twice the sup-norm distance" | **right** (0.119964) | [own] item 8 |
| md:621 lead, `mirror_asym`'s opening sentence ("distance from the … subspace") | **stale** (D3dg-4) | [md] |
| `mirror_odd`'s resolving grid | **resolves every kernel on the data** (own maximiser to 4–5 digits; T09 equivalent) | [own] |
| (ii)(b) "0.24 from their mirror image (… 0.12 …; 0.48 with all 17)" | **holds** | [own] |
| (ii)(b) "their function is odd by 5.5·10⁻⁹ only, on and between the samples" | **number right, threshold-fragile** (2% margin; 4.9·10⁻³ at 0.95%) (D3dg-2) | [own] item 8 |
| (ii)(b) "the whole state's 1.2·10⁻⁹ … the residual's, its 4.9·10⁻³ … the unformed kernels' spikes" | **holds** | [own] |
| (ii)(b) "a symmetry breaking of the decomposition, not of the fit" | **holds**; threshold-free support: the kernels wider than the spacing, 1.19·10⁻⁹ | [own] |
| Fenichel clause "on its compact pieces where it is strictly repelling" | **too narrow**: normally hyperbolic pieces include strictly attracting ones (the stable phase) (D3dg-4) | [memory] |
| Story of (i), statement → lemma → perturbation → numerics | **consistent**; the quantitative mechanism (the route seed) is unstated; the statement's parenthetical is overlong | re-read, [own] |
| tests | **strong**; T01, T03 untested conditions | [own] item 10 |

## Status of the round-109 issues

| D3dd | status | what remains |
|---|---|---|
| 1 (the exit's seed) | **applied**: fsum and reversed runs, conditions on directions, times and the offset start, kz from the 50-digit zero, plants | The arithmetic's conditions are identities. The symmetric direction also follows the tolerance and the start. The 5% window is within double scatter. Long double and the route seed are not adopted; the log says "not reproduced", and they are reproduced here (D3dg-1). |
| 2 (the mirror evidence) | **applied**: halving, resolving grid, per-target sentence, D59 | The formed threshold's 2% margin (D3dg-2). The stale lead and docstring opening (D3dg-4). |
| 3 (−d/λ, the centre sample) | **applied** | — |
| 4 (stale text) | **applied** (md:54, docstring, "LSODA and scipy's BDF hold 3A") | BDF's hold is tolerance-dependent (D3dg-3). |
| 5 (grid floors) | **applied**: roots, verified one per bracket | `k` unused (D3dg-4). |
| 6 (tests) | **applied**: plants for R06, R09, the exact-sum run leaving, A's time, the flip; 3A's exit; 0.75 ulp; one "held" criterion; 10 mutants all killed (my rerun) | T01, T03 survive (D3dg-5). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dg-1 | minor | **The exit's printed tests are weak; the decisive ones are missing.** (`a2_exits`, `a2_lines` exit line; md:62; D59 (i) and A2E.) <br>• On exactly symmetric samples, fsum makes ċ = 0.0 exactly (c ≡ 0.0 over 490 / 2337 steps), and reversal is the mirror map (times equal bit for bit). `arith`'s two conditions cannot fail [own]. <br>• The symmetric problem's direction is not set by the order alone. With the code's own sums, rtol 10⁻¹³ (3A), 10⁻¹¹ (5A) or α −2 ulps (5A) flip it, and its times scatter ±8%. A's: ±1%, always (−) [own]. <br>• Offset start: the long-double value is +2.64% / +1.75%, because the route seed is 2.34× / 1.94× the held offset. In double the same start gives 0.04899–0.05195 (3A) across my variants, and 0.050688 (+4.99%) with the repository's `_FsumBumps`. My independent double implementation would not print the "data's" clause (+8.5%) [own]. <br>• Decisive [own]: A's long-double δ-state exit is 0.0482973 / 0.0034820, tolerance-free to 10⁻⁶. The symmetric problem's exit moves 10–34% later in long double. Seeded with S = −7.762·10⁻¹⁷ / −5.473·10⁻¹⁷, the symmetric problem leaves at A's time to 5·10⁻⁶. | **Code:** add `a2_exit(arith='ld')`: a long-double right-hand side with the offset δ = c − c* as the state. Add `a2_seed`: S = ∫e^{−Λ}d in long double on that run's dense output (d and λ in closed form). Print A's and the symmetric problem's long-double exits and the symmetric problem started at S. Condition "the data's" on \|t_sym(S)/t_A − 1\| < 10⁻⁴ (long double) and on A's double runs lying within 2% of the long-double time. State fsum and reversal as identities. <br>**md / D59 (i), A2E:** "on exactly symmetric samples the exit is the rounding's: exact sums keep c at the centre exactly, reversing the sums mirrors the run, the tolerance or the start's last bits can flip it, and in long double it leaves 10–34% later; on A's samples the data's drift integrated along the route seeds it (−7.8·10⁻¹⁷, −5.5·10⁻¹⁷: 2.3 and 1.9 times the held offset), and the symmetric problem started there leaves at A's time to 5·10⁻⁶". <br>**Tests:** a plant with `symoff` at +10% (kills T01) and one with A's fsum/rev exits upward (kills T03). |
| D3dg-2 | minor | **"The formed kernels' function odd by 5.5·10⁻⁹ … between the samples" rests on the formed threshold.** (D59 (ii)(b); md:621.) <br>• The largest kernel below 0.01·A is −0.009803, an invisible spike: σ 3.25·10⁻⁴, φ ≤ 1.3·10⁻¹⁷ on the samples, 2% under the threshold. <br>• At thresholds of 0.0095·A or lower the formed set's function is odd by 4.90·10⁻³. T14 (0.9%) changes the published number [own]. <br>• Threshold-free: the 11 of 17 kernels with σ ≥ h (the narrowest 0.041, the widest spike 1.0·10⁻³: a ×40 gap) are odd by 1.19·10⁻⁹ on and between the samples [own]. | **D59 (ii)(b) and md:** "the kernels at least as wide as the samples' spacing (11 of 17) form a function odd by 1.2·10⁻⁹ on and between the samples; the narrower ones, invisible or nearly, make the whole state's 4.9·10⁻³ odd between them". Or keep "formed" and print the margin (the largest unformed \|w\| 0.0098, the smallest formed 0.060). **Code:** `split_summary`'s mirror block: add the σ ≥ h set. |
| D3dg-3 | cosmetic | **D59 wording.** <br>• (i): "at the same time to 2% … as does the symmetric problem started at the posed problem's offset": that start is +2.6% at 3A (md: 5%). <br>• (i) "(LSODA and scipy's BDF hold 3A)" and A2E "(which holds 3A)": scipy's BDF holds at rtol 10⁻¹⁰ (also to s = 10), and leaves at 10⁻¹², 0.0516 (−) [own]. <br>• "Its direction follows the sums' order" (D3dg-1). | "… or reversed (0.0480, 0.0035), and the symmetric problem started at the posed problem's offset within 5% (0.0495, 0.0035)"; "(LSODA, and scipy's BDF at rtol 10⁻¹⁰, hold 3A; at 10⁻¹² BDF leaves at 0.0516, downward)". |
| D3dg-4 | cosmetic | **Stale or loose text.** <br>• md:621 and the code at 1692: "the formed kernels' distance from the reflection's fixed-point subspace (`mirror_asym` …)", while the entries say "from the mirror image (… from the subspace)". `mirror_asym`'s docstring opens with "the distance of a configuration from the reflection's fixed-point subspace". <br>• `edge(pred, k)`: `k` unused. <br>• Fenichel clause "where it is strictly repelling": normally hyperbolic pieces include the strictly attracting stable phase (s < 6.3·10⁻⁴ / 3.2·10⁻³) [memory]. <br>• Statement (i) holds a 12-line numerical parenthetical on the symmetric state. | "the formed kernels' distance from their mirror image (`mirror_asym` …; half of it from the fixed-point subspace)"; the same in the docstring; drop `k`; "on its compact pieces where it is normally hyperbolic (strictly attracting or repelling in c)"; move the exit and hold numerics from (i) to a remark after Lemma D‴.59a. |
| D3dg-5 | cosmetic | **Tests.** <br>• 8 of my 16 mutants survive, all equivalent on the data. T01 (the offset start's 5% unchecked) and T03 (A's fsum/rev directions dropped) are printed conditions with no plant. <br>• The plants for S05/S06 exercise identities. <br>• The lemma's sign probe covers only the totals. | Plants as in D3dg-1. A long-double exit test (3A: 0.04829 < t < 0.04831 at both tolerances). The ±1-term roots in the sign probe. |

## Sources

- [snippet] Python `math.fsum` docs: "Return an accurate floating point sum of values in the iterable. Avoids loss of precision by tracking multiple intermediate partial sums". Accuracy depends on IEEE-754 and half-even rounding; extended-precision builds may double-round ([Python 2 docs](https://docs.python.org//2//library/math.html), [3.4 source](https://docs.python.org/3.4/_sources/library/math.txt)). Shewchuk's algorithm per [accupy (PyPI)](https://pypi.org/project/accupy) and [PreciseSums (CRAN)](https://search.r-project.org/CRAN/packages/PreciseSums/PreciseSums.pdf). This fits "exactly rounded" here (x86-64 SSE doubles); that CPython's result is correctly rounded is [memory].
- [snippet] NumPy `numpy.sum` notes: partial pairwise summation, used along an axis only when it is the fast (contiguous) axis; `math.fsum` is slower but more precise ([numpy 2.0 docs](https://numpy.org/doc/2.0/reference/generated/numpy.sum.html)). This fits the different orders of the code's (N, 1) `.sum(0)` and my 1-D sums.
- [snippet] Berglund and Gentz, "Pathwise description of dynamic pitchfork bifurcations with additive noise" ([arXiv:math/0008208](https://arxiv.org/abs/math.PR/0008208); [HAL](https://hal.archives-ouvertes.fr/hal-00130541)): the exit from the unstable branch depends on the noise intensity through log σ. This is why time windows resolve seeds only logarithmically and why a precision change is the stronger lever.
- [snippet] U. Mutze, mp_arc 09-113 ([abstract](https://web.ma.utexas.edu/mp_arc-bin/mpa?yn=09-113)): rounding errors destroy a symmetric behaviour in simulation, shown by computing at 16, 20, 30 and 40 digits. This is the precision test proposed in D3dg-1.
- [snippet] Forest, Goedde and Sinha, Math. Comput. Simul. 37 (1994) ([RePEc](https://ideas.repec.org/a/eee/matcom/v37y1994i4p323-339.html)): errors of the order of machine arithmetic and roundoff can break a symmetry of the data and the equations.
- [snippet] Bubbling transition, Phys. Rev. E (1996) ([OSTI](https://www.osti.gov/biblio/286590)): with noise or asymmetry, an invariant subspace's attractor gives way to transients or bursts. Background only.
- [memory] Fenichel 1971 on normally hyperbolic invariant manifolds (attracting, repelling or saddle type); Golubitsky–Stewart–Schaeffer on fixed-point subspaces. Both are cited by D59 from earlier rounds and were not re-searched.

## Lesson

A change of the arithmetic tests a seed only if the change could have removed it.
- On exactly symmetric data, exact sums and reversed order are identities. Exact sums leave c exactly at the centre, and reversal mirrors the run bit for bit, so they confirm the symmetry but cannot discriminate.
- The levers that separate the data's seed from the arithmetic's are two:
  - **Precision.** In long double the arithmetic's exit moved 10–34% later and changed direction with the tolerance, while A's did not move (to 10⁻⁶).
  - **Planting the data's seed.** The drift integrated along the route reproduces the exit exactly when planted in the symmetric problem.
- A matched-offset start in double cannot do this: its time carries the double arithmetic's own scatter of several percent.

A threshold that splits kernels needs a gap.
- A spike 2% below "formed" makes "the formed kernels' function" another function at a threshold of 0.95% of A.
- Split by width against the samples' spacing instead: the data fix the function of the wide kernels, and here the gap is a factor 40.
