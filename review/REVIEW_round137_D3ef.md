# Referee report D3ef, round 137: Proposition D‴.59 as revised after round 133

## Scope

- **Export and work directory.** I exported commit e53ec0d with `git -C … archive e53ec0d | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref137. W/repo/venv is a symlink to the repository's venv.
- **Settings.** MW_OUT = W/out and MW_WORKERS = 1. Every run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, one heavy process at a time, run as chained background jobs.
- **Process notes.**
  - I killed no process by name. With TaskStop I stopped only my own tasks: one float32 adaptive run that stalled (see "What I checked", item 6), one chain I restarted after a bug, and two redundant waiters.
- **The repository is untouched.** I used only read-only git.
  - HEAD moved to 91cb3fe through other sessions' commits. `git diff --stat e53ec0d HEAD` over the D59 files (exp_moving_width.py, test_moving_width.py, moving_width.md and its raw, D59.md) is empty.
  - The untracked manuscript/paper16/ belongs to another session.
- **Read in full.**
  - review/drafts/D59.md (482 lines).
  - review/REVIEW_round133_D3eb.md and the Round 133 entry of REVIEW_LOG.md.
  - tools/mutants/round133_D59_own.py.
  - `git diff 80205c8 a80663c` for the D59 files. Only a80663c touches them in that range.
  - The code the new parts call: `a2p_job`, `a2_exit`, `a2_exit_ld`, `_LDBumps`, `BumpsCW.corr`, `critical_rhs_nd`, `a2p_stats`, `a2p_exp`, `a2p_blocked`, `a2p_line`, `a2p_shift`, `a2h_line`, `a2r_jobs`, `part_A2R`, `a2r_line`, `nh_verdict`, `a2_line_lam`.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my computation.
  - [snippet]: a web-search result I saw (no page fetched).
  - [memory]: not verified this round.

## Summary

**Design.**
1. **Reproduction.** The md byte for byte from `--from-raw`. The 58 tests. A2R rows recomputed with the repository's `a2p_job` on four sampled jobs.
2. **Statistics from the raw, recomputed independently.**
   - The interval rule's stability: per seed, seed pairs, disjoint replicates, and a bootstrap of the pass probability against n.
   - Its sensitivity to the band.
   - Gating on replication alone against the pooled data.
   - Heterogeneity across configurations (DerSimonian–Laird).
   - Robust statistics (median, trimmed mean).
3. **Two own experiments on the mechanism.**
   - **(a) A third arithmetic.** The repository's right-hand side in float32 as well as double and x87 long double (bitwise equal to `BumpsCW` and `_LDBumps` at those precisions). Integrated by fixed-step RK4, so every arithmetic takes the identical step sequence. Same summation orders across arithmetics. Fresh seeds 137 and 138.
   - **(b) A scan of the computed ċ(δ) near the symmetric line**, at A2E-like states at 3A and 5A, and the repository's own DOP853 3A and 5A runs traced in ulps.
4. **D59's changed claims** against the md and raw.
5. **Mutants.** The authors' M00–M24, and 15 of my own (R01–R15) plus a control, each against `-k "round133 or round130"` and the `--from-raw` md.
6. **Literature.** Seven web searches.

**Verdicts.**
- **Reproducibility is exact.**
  - The md is identical (sha256 9ce33dfd…84ec, 32.6 s).
  - 58 passed in 300.9 s, with the known LSODA warning.
  - A2R's sampled rows are bitwise equal in exits, directions and step counts. ln S and Λ agree to ≤ 1.1·10⁻⁷, from the line's grid.
- **Every new number in D59 is the md's**, and the md's numbers are the raw's.
- **The authors' 24 mutants are all killed.**
- **The revision answers round 133's findings as asked.** But four things do not hold as printed:
  - the interval rule's verdict is a near coin-flip at the design's n;
  - the "not chance" and "differ at equal steps too" verdicts lean on the discovery draws;
  - the "systematic spread" is a range of noisy estimates;
  - two arithmetics cannot separate the exponent from the prefactor.
- **The decisive new result is mechanistic** [own].
  - **5A: the scaling is settled.** With a third arithmetic at identical steps, the exponent is 1 to within 1% (0.9926 ± 0.0044 and 1.0059 ± 0.0029, three-point fits over a 40 ln 2 lever arm). The prefactors differ by up to e^{±0.26}.
  - **3A: not a rounding-noise seed at all.** The computed ċ(δ) has attracting zeros at about ±0.5 ulp of the ±1 samples' coordinate. Under fixed RK4 steps of h₀/4 neither the double nor the long double leaves by s = 0.5 (0 of 30 each). At h₀/2 and h₀ they leave, with effective seeds falling steeply with the step. So the 3A cell's "prefactor difference" (round 133's D3eb-1) is real, but it is the integrator's escape from ulp-scale traps.

**Findings.**
- **D3ef-1 (major):** at 3A the symmetric exit is a step-sequence escape from rounding traps, not a noise-seeded linear escape. Its existence and its seed are the integrator's. The Remark's "on exactly symmetric samples from their centre the exit is the rounding's", and the reading of the 3A excess as an arithmetic prefactor, need qualifying.
- **D3ef-2 (minor, borderline major):** the interval rule is precision-limited at n = 200 and gated on A2P alone. A fresh 200-draw replicate fails it with probability ≈ 0.38, and one disjoint replicate does fail it. The band [0.9, 1.1] was set after the data, and [0.95, 1.05] fails at every n. The band is inconsistent with the geometric-mean band.
- **D3ef-3 (minor):** two arithmetics cannot identify an exponent once the prefactors are declared unequal. The "exponent" Λ/(11 ln 2) assumes the very equality the text then denies. A third arithmetic settles it at 5A (D3ef-1's runs).
- **D3ef-4 (minor):** "not chance", "not the atol" and "differ at equal steps too" are gated on pools that include the selected discovery draws. On the replications alone the t values are 2.97 and 2.59. The code's own 3-sem gate would print "not resolved" (my mutants R01, R12, R14 are killed: the gates flip). "+3.0 sem" is printed for 2.97.
- **D3ef-5 (cosmetic):** "its systematic spread 0.979 … 1.050" is the range of 8 estimates with sems of 0.011–0.037. Noise alone gives a range at least this wide with p = 0.24. The heterogeneity is real (Q = 27.7 on 7 dof), and its size is τ ≈ 0.03.
- **D3ef-6 (cosmetic):** the Brouwer wording "3A's long double not at all" contradicts the md's own count, since at h₀/4 it lies 2.96 sem from the law.
- **D3ef-7 (cosmetic):** test gaps. R06, R09 and R10 survive, all equivalent on the data. The interval rule is never exercised on the pooled 3A cell.

## What I checked and how

Scripts are in W/own/:
- stats137.py, het137.py, robust137.py (statistics from the raw);
- recompR.py (row recomputation);
- f32_137.py (a dtype-generic `_LDBumps`, and the bitwise checks);
- rk3_137.py (three arithmetics at identical fixed steps);
- trap.py and trap5.py (the scan of ċ(δ));
- trap_dop.py (the repository's DOP853 runs traced in ulps);
- an_rk3.py (analysis);
- spec137.py (my mutants).

Logs are in W/logs/, data in W/own/rk3_*.json, and mutant records in W/mut/res_r133.jsonl and res_own.jsonl.

1. **Reproducibility** [own].
   - `--from-raw` takes 32.6 s; `cmp` against the committed md: identical.
   - pytest: 58 passed in 300.9 s. The round133 test takes 51.4 s.
   - recompR.py runs `a2p_job` at:
     - (3A, 10⁻¹⁰, seed 135, k = 3);
     - the matched arm on seed 136 and on seed 131 (atol_d = 10⁻¹⁴σ*, k = 2);
     - the cap on seed 132 (5A, h₀/16, k = 2).
   - Exits, directions and step counts are bitwise those of the raw. max |Δ ln S, ΔΛ| is 2.5·10⁻⁸, 7.1·10⁻⁸, 1.1·10⁻⁷ and 2.7·10⁻¹⁰.
2. **The interval rule** [own], [raw].
   - **At 3A, rtol 10⁻¹⁰**, each seed's exponent with its 99% interval:

     | seed | exponent | 99% interval | rule |
     |---|---|---|---|
     | 131 | 1.093 | 1.026 … 1.160 | fails |
     | 132 | 1.007 | 0.939 … 1.075 | passes |
     | 133 | 1.021 | 0.946 … 1.096 | passes |
     | 134 | 1.045 | 0.980 … 1.111 | fails |
     | 135 | 1.056 | 0.984 … 1.129 | fails |
     | 136 | 1.044 | 0.967 … 1.120 | fails |

   - **At A2P's design (n = 200):**
     - 5 of the 15 seed pairs fail.
     - Of the disjoint replicates, (131, 132) passes at 1.098, (133, 134) passes at 1.083, and (135, 136) fails at 1.102.
     - Bootstrapping from the 600, the rule fails with probability 0.71 at n = 100, 0.38 at n = 200, 0.071 at n = 400 and 0.009 at n = 600.
   - **On the pooled 600 it passes with a margin** (1.015 … 1.073). Yet `ok` is computed on A2P's 200 per cell only.
   - **With the band [0.95, 1.05] the rule fails at every n** (the 600's upper end is 1.073).
   - **The band's meaning.** The band ±0.1 is |ln(C_d/C_ld)| ≤ 0.76, a prefactor ratio within a factor of 2.1. The geometric-mean band 0.05–20 uσ* admits a ratio of up to e^{6}, i.e. an "exponent" of 0.21–1.79.
   - **Seed 131 at 3A, rtol 10⁻¹² also fails per seed** (0.977 … 1.138).
3. **Gating on replication** [own], [raw].
   - **The fresh 400 alone:** 7.9415 ± 0.1067, t = 2.97 (printed "+3.0"), exponent 1.042 (1.006 … 1.078). With the 200 selected as the farthest of four cells: 7.963 ± 0.086, t = 3.95.
   - **The matched-atol arm:** 7.928 ± 0.091, t = 3.33, exponent 1.040 ± 0.012. Its long-double runs are bitwise those of the default arm, so it is not independent evidence of a prefactor difference.
   - **The 5A h₀/16 cap on seed 132 alone:** t = 2.59 (exponent 1.026, 1.000 … 1.052). Pooled with the selected seed 131: +5.1.
   - **The two replications combined** (Stouffer): 3.93. So the effect is likely real. The verdicts as coded, however, rest on the selected draws.
   - **Robust statistics agree.** The 3A, rtol 10⁻¹⁰ median exponent is 1.045 on the 600 (bootstrap 99% 1.015 … 1.076), the 10%-trimmed exponent 1.045. One to two draws per cell lie beyond 4 robust sd.
4. **Spread and heterogeneity** [own].
   - The md's 8 configurations: Q = 27.7 on 7 dof (p = 2.5·10⁻⁴); DerSimonian–Laird τ = 0.029; random-effects mean 1.009 ± 0.012 (99%: 0.978 … 1.040).
   - A2P's 4 cells alone: Q = 7.3 (p = 0.06), τ = 0.021. A2H's 4 caps: Q = 20.0 (p = 1.7·10⁻⁴), τ = 0.040.
   - Under homogeneity, the range of 8 estimates with these sems has a median of 0.056 and P(≥ 0.071) = 0.24.
5. **Three arithmetics at identical steps** [own].
   - **Bitwise checks.** My dtype-generic class at longdouble is `_LDBumps` bitwise. At float64 its right-hand side is `BumpsCW`'s bitwise, at c₀ = 0 on three states.
   - **5A, h₀/4 (1.71·10⁻⁵), 200 orders:**
     - ln(S/uσ*) means: float32 −2.59 ± 0.08, double −2.59 ± 0.29, long double −2.33 ± 0.05.
     - The double carries one order at −59, an exact cancellation. Medians: −2.22, −2.03, −2.10.
     - Three-point slope per order: 0.9926 ± 0.0044; the median slope 0.9975 (bootstrap 99% 0.992 … 1.001).
     - The float32/long-double exponent is 0.9907 ± 0.0032, with ln(C_32/C_ld) = −0.26 ± 0.09.
   - **5A, h₀/16, 100 orders:**
     - Slope 1.0059 ± 0.0029 (median 1.0043).
     - Pairwise: f32/d 1.0064 ± 0.0046, d/ld 1.0038 ± 0.0147.
     - ln(C_32/C_ld) = +0.16 ± 0.09.
   - So the exponent is 1 to within 1%, while the prefactors are arithmetic- and step-dependent at the ±0.3 level. Two arithmetics alone resolve the exponent only to about ±0.04.
   - **3A, fixed RK4:**
     - At h₀/4 (2.18·10⁻⁴, 30 orders), neither double nor long double leaves by s = 0.5 (0 of 30 each). float32 leaves at mean s = 0.145, effective seed e^{−77} uσ*.
     - At h₀/2 (10 orders), the double and long double leave at s = 0.061 and 0.072, with median seeds e^{−5.2} and e^{−4.4} uσ*.
     - At h₀ (10 orders) they leave at s = 0.055 and 0.067, with seeds e^{−0.6} and e^{−1.3}. This is DOP853-like.
     - The route of (α, σ) agrees with DOP853 to 6 digits; the difference is only in δ.
6. **The mechanism** [own].
   - **At the 3A line state** (σ/σ* 0.2345, s ≈ 0.06, one order), the computed ċ(δ) on |δ| ≤ 4 ulp(0.01) is a sawtooth.
     - Its slope below half an ulp is the centre sample's alone (about 4.4λ).
     - At |δ| = ½ ulp the ±1 samples' coordinate x − δ rounds and ċ jumps down.
     - The result is attracting zeros at δ ≈ ±0.51 ulp, in long double, double and float32 alike (and ±1.5 ulp in long double, +2.5 in float32).
     - The rms of ċ − λδ is 0.6–0.8 λuσ*.
     - ulp(0.01) = 0.521 uσ* in every arithmetic, since 0.01 lies in the same binade.
   - **At 5A** (σ/σ* 0.1846) there are no attracting zeros, and the rms is 0.07 λuσ*.
   - **The repository's own DOP853 3A runs** (trap_dop.py, exits and step counts bitwise the raw's, rows 1–3 of seed 131) hold |δ| at 0.3–2.6 ulp during the transient. They spend up to 9 steps with 0.3 < |δ|/ulp < 0.7 and leave the sub-ulp region at s = 0.0010, 0.0045 and 0.0072. The 5A runs leave it by s ≤ 0.001.
   - **A side result:** DOP853 at rtol 10⁻¹⁰ with a float32 right-hand side stalls (the error control chases the rounding). So the third arithmetic needs fixed steps, or an atol scaled with u.
7. **Mutants** [own].
   - **The authors' spec:** the control passes. M01–M24 are all killed by the tests; M01, 02, 05, 07–12, 19, 21, 22 and 24 also change the md.
   - **Mine (spec137.py):** the control passes.
     - Killed: R01, R02, R03, R04, R05, R07, R08, R11, R12, R13, R14 and R15.
       - R01: "not the atol" on the fresh 400's t.
       - R12: the cap's verdict on the replication seed alone.
       - R14: "not chance" on the fresh 400.
       - R15: the gate on one seed per cell.
     - Survive, all equivalent on the data:
       - R06: the own-atol steps over the fresh seeds only (60.30 against 60.28; both print 60.3).
       - R09: the largest offset by signed t (the largest is positive).
       - R10: the geometric-mean band closed.
8. **D59 against the md** (word diff of a80663c).
   - Every new number is verbatim the md's: 1.050/1.017/0.985/0.999; the 99% intervals; ±0.0098; 7.963 ± 0.086 (+3.9); 7.941 ± 0.107 (+3.0); 1.044 (1.015 … 1.073); +0.26/−0.08/+0.34; 7.928 ± 0.091 (+3.3); 67.7/60.3/73.1; 7.825 ± 0.077 (+2.6); 7.912 ± 0.056 (+5.1); 1.038; 0.979 … 1.050; the Brouwer −0.69/−1.39; 2 of 8.
   - The Status labels Remark D‴.59b numerical, which is right.

## Verdicts per part

| claim (D59 unless noted) | verdict | evidence |
|---|---|---|
| md byte for byte; 58 tests; A2R rows | **holds** | [own] |
| "each cell's exponent … 1.050, 1.017, 0.985 and 0.999, its 99% interval inside [0.9, 1.1]" | **numbers hold.** As a gate it is a near coin-flip at n = 200 (P(fail) ≈ 0.38; one disjoint replicate fails); it is computed without A2R's 600, at which it passes robustly (D3ef-2) | [raw], [own] |
| "the seed scales with the arithmetic's roundoff" | **holds at 5A, more strongly than shown**: three arithmetics at equal steps, slope 0.993 ± 0.004 and 1.006 ± 0.003. **At 3A it holds only through the ulp structure**, with the exit's existence and prefactor set by the steps (D3ef-1, D3ef-3) | [own] |
| "its prefactor not the same in the two arithmetics … about e^{0.34}, not by chance and not through the atol" | **real but mislocated.** At 3A it is the trap × step sequence, not the arithmetic. "Not chance" rests on the selected draws (2.97 on the replication alone) (D3ef-1, D3ef-4) | [own] |
| "1.013 ± 0.009, statistical (0.0098 over the orders)" | **correctly computed.** Under the measured heterogeneity (τ ≈ 0.03) the exponent's uncertainty is ±0.012 (random effects) (D3ef-5) | [own] |
| "its systematic spread … 0.979 … 1.050" | **mislabelled**: a range of noisy estimates (P = 0.24 under pure noise). The heterogeneity is real, τ = 0.029 (D3ef-5) | [own] |
| Brouwer: "5A follows at h₀/4 and 3A's long double not at all, its mean not moving" | **partly.** 3A's long double is 2.96 sem from the law at h₀/4, inside the md's own 3-sem count; 3A's double follows at both caps; 5A's double is off at h₀/16 (+4.15) (D3ef-6) | [md] |
| "their prefactors unequal at equal steps too, 5A at h₀/16 +4.7 sem, on a second seed +2.6, +5.1 over both" | **numbers hold.** The replication alone is +2.6; the verdict relies on the selected seed (D3ef-4). My equal-step three-arithmetic runs at 5A, h₀/16 give d/ld 1.004 ± 0.015 with RK4 | [md], [own] |
| "on exactly symmetric samples from their centre the exit is the rounding's" | **needs qualifying at 3A**: under fixed RK4 steps at h₀/4 there is no exit in double or long double by s = 0.5 (D3ef-1) | [own] |
| proved/numerical labels | **correct** | re-read |
| tests | **good.** All 24 authors' mutants and 12 of my 15 are killed; 3 equivalent survivors (D3ef-7) | [own] |

## Status of the round-133 issues

| D3eb | status | what remains |
|---|---|---|
| 1 (3-sem rule) | **applied**: the interval rule, A2R, matched atol | The rule is precision-limited and gated on A2P only; the band is post hoc (D3ef-2). The 3A excess is the trap × steps (D3ef-1) |
| 2 ("Λ compares at equal steps") | **applied**: dropped; spread printed | "Systematic spread" is a range (D3ef-5); the equal-step verdict leans on the selected seed (D3ef-4) |
| 3 (unpaired shifts, Brouwer) | **applied** | Wording (D3ef-6). The 3A non-Brouwer behaviour is now explained by the traps (D3ef-1) |
| 4 (`nh_verdict` edge, metric) | **resolved** (R03 killed) | — |
| 5 (extremes, blocked sem) | **resolved** | — |
| 6 (test gaps) | **resolved** (M01–M24 killed) | R06, R09, R10 equivalent; the rule is untested on the pooled cell (D3ef-7) |

## Overall assessment

(i)–(iii) of D‴.59 are untouched by round 133 and stand as before.

**Remark D‴.59b (the symmetric exit) is settled at 5A.**
- With the three arithmetics at identical steps, the seed scales as u^γ with γ = 1 to within 1%.
- The prefactors are O(1) and arithmetic- and step-dependent at e^{±0.3}.
- This is a better statement than any two-arithmetic interval can give.

**It is not settled at 3A, and the revision's new reading of the 3A cell is wrong in mechanism.**
- There the computed field has ulp-scale attracting zeros. Whether and when the flow leaves the symmetric line is decided by the integrator's steps:
  - DOP853 at its own or capped steps leaves;
  - RK4 at h₀/2 and h₀ leaves, with seeds e^{−5} to e^{−1};
  - RK4 at h₀/4 does not leave by s = 0.5.
- S ∝ u holds there only because ulp(x₁) = 0.521 uσ* in every arithmetic.

**What would still change the verdict:**
1. The trap scan at more orders and states, and on A's samples (the posed problem, where c* is not 0).
2. A third arithmetic in the repository, with fixed or u-scaled atol.
3. Random rounding (CESTAC-style) instead of summation permutations, as a principled draw of the rounding.
4. Another explicit integrator family at 3A.
5. Whether the 5A result persists at DOP853's own adaptive steps with a third arithmetic.

**Recommendation:** minor revision with one major item (D3ef-1). The Remark should say that at 3A the exit is a trap escape whose existence and prefactor are the integrator's, and it should report the 5A three-arithmetic exponent.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ef-1 | major | **At 3A the symmetric exit is an escape from rounding traps decided by the step sequence, not a rounding-seeded linear escape.** <br>• ċ(δ) at the 3A line state has attracting zeros at δ ≈ ±0.51 ulp(0.01) in float32, double and long double (rms of ċ − λδ 0.6–0.8 λuσ*); none at 5A (0.07) [own]. <br>• Fixed RK4: at h₀/4, 0 of 30 exits in double and in long double by s = 0.5; at h₀/2 seeds e^{−5.2}, e^{−4.4}; at h₀ e^{−0.6}, e^{−1.3} [own]. <br>• The repository's DOP853 3A runs sit at 0.3–2.6 ulp for up to 9 steps before leaving (bitwise the raw's exits) [own]. <br>• This explains why 3A's long double ignores Brouwer and why 3A's "prefactor" differs. | **D59 Remark D‴.59b and Status:** "at 3A the computed field has ulp-scale attracting zeros (the ±1 samples' coordinate rounding at ½ ulp); the exit's existence and seed are the integrator's (fixed RK4 at h₀/4 holds the line), S ∝ u through ulp(x₁) = 0.52 uσ*". Reword A2R's "the two arithmetics' prefactors differ" at 3A as the step sequence's. **New part:** a ċ(δ) scan per cell (`a2_trap`) and a fixed-step arm. |
| D3ef-2 | minor (borderline major) | **The interval rule is precision-limited and gated on A2P alone.** <br>• The 3A cell's upper end is 1.098. A fresh 200-draw replicate fails with P ≈ 0.38 (bootstrap); 5 of 15 seed pairs fail; the disjoint (135, 136) fails at 1.102; 4 of 6 seeds fail at n = 100 [own]. <br>• It passes robustly on the 600 (P(fail) 0.009). <br>• The band [0.9, 1.1] was set after the data; [0.95, 1.05] fails at every n. It means a prefactor ratio within e^{±0.76}, while the geometric-mean band admits e^{±6}. | **`a2p_line`:** gate each cell on all its draws (A2P with A2R). State the band as a prefactor ratio and justify it from the alternatives to be excluded (for example γ = ½ or 3/2), or drop the gate and print the interval only. **D59:** "the rule's pass at n = 200 is marginal at 3A". |
| D3ef-3 | minor | **With two arithmetics the exponent and the prefactor ratio are not identifiable.** "Λ/(11 ln 2)" assumes the very equality the text then denies [own]. <br>• A third arithmetic at equal steps gives 0.9926 ± 0.0044 (h₀/4, n = 200) and 1.0059 ± 0.0029 (h₀/16, n = 99) at 5A, with ln(C_32/C_ld) −0.26 ± 0.09 and +0.16 ± 0.09. | **New part (A2T):** float32 (or any third precision) at fixed steps, or DOP853 with atol ∝ u. A three-point slope per order. **D59:** report γ from three arithmetics, and the prefactors as a separate O(1) scatter. |
| D3ef-4 | minor | **The verdicts rest on the selected discovery draws.** <br>• "Not chance" uses 600 draws including the cell chosen as the farthest of 4. The fresh 400 alone give t = 2.97 (printed "+3.0"). <br>• "Differ at equal steps too" relies on the farthest cap of 4; the replication seed alone gives +2.59. <br>• R01, R12 and R14 (gates on the replication alone) are all killed: the printed verdicts flip [own]. <br>• Stouffer over the two replications: 3.93, so the effect is likely real. | **`a2r_line`:** gate on the replication alone, or print both. Give t to two decimals. Name the selection. **D59:** "on the further draws alone +2.97 and +2.59 sem (combined 3.9)". |
| D3ef-5 | cosmetic | **"Its systematic spread 0.979 … 1.050" is a range of 8 estimates with sems of 0.011–0.037.** <br>• Noise alone gives a range of at least 0.071 with P = 0.24. <br>• The heterogeneity is real (Q = 27.7 on 7 dof, p = 2.5·10⁻⁴): DL τ = 0.029, random-effects mean 1.009 ± 0.012 [own]. | **`a2h_line` and D59:** "between-configuration sd τ ≈ 0.03 (DerSimonian–Laird), the exponent 1.009 ± 0.012 over configurations", beside the range. |
| D3ef-6 | cosmetic | **Brouwer wording.** "3A's long double not at all": at h₀/4 it is 2.96 sem from the law, inside the md's own count, where only 2 of 8 are beyond 3 sem, both at h₀/16. The wording also omits that 3A's double follows at both caps and that 5A's double departs at h₀/16 (+4.15) [md]. | **D59:** "2 of 8 beyond 3 sem of Brouwer's law: 3A's long double and 5A's double at h₀/16", plus the trap explanation (D3ef-1). |
| D3ef-7 | cosmetic | **Test gaps.** R06 (the step count over the fresh seeds, 60.30 against 60.28), R09 and R10 survive, all equivalent on the data. No test exercises the rule's verdict on the pooled 3A cell, or a plant at the 3A cell's own edge [own]. | **Test:** print the steps to two decimals or assert the base; plant the interval edges on the 3A cell; with D3ef-2's fix, assert the gate's input is A2P with A2R. |

## Sources

- [snippet] Bergeot, Almeida, Vergez, Gazengel, "Prediction of the dynamic oscillation threshold in a clarinet model with a linearly increasing blowing pressure" ([arXiv:1207.4636](https://arxiv.org/pdf/1207.4636); Nonlinear Dyn. 74, 2013): the bifurcation delay is very sensitive to numerical noise. It vanishes at 7 digits, matches theory at about 5000 digits (mpmath), and grows with precision in between. This is the precedent for a many-precision design (D3ef-3).
- [snippet] Bergeot et al., the comparison paper ([arXiv:1407.3782](https://arxiv.org/pdf/1407.3782)): in double precision the predicted thresholds are far from the simulated ones; finite-precision thresholds are predicted following Baesens. The noise paper ([arXiv:1304.2294](https://arxiv.org/pdf/1304.2294)) models round-off as additive white noise tied to the number of digits.
- [snippet] "Slow passage through a Hopf-like bifurcation in piecewise linear systems" ([arXiv:2111.12416](https://arxiv.org/pdf/2111.12416)): the way-in/way-out saturation matches the integrator's precision.
- [snippet] Baer, Erneux, Rinzel, SIAM J. Appl. Math. 49 (1989): delay and memory in slow passage. Bilinsky and Baer ([NIMBioS abstract](https://legacy.nimbios.org/SMB2012/pdfs/36_Lydia_Bilinsky.pdf)): round-off is a serious issue with slow ramps.
- [snippet] Baesens, "Slow sweep through a period-doubling cascade", Physica D 53 (1991) ([Warwick](https://wrap.warwick.ac.uk/22337/)): small noise plays a key role in delayed bifurcations.
- [snippet] A 1994 Math. Comput. Simul. paper ([IDEAS](https://ideas.repec.org/a/eee/matcom/v37y1994i4p323-339.html)): machine-arithmetic errors break a symmetry of the data and open an instability absent in the symmetric simulation. Also [arXiv:1910.11976](https://arxiv.org/pdf/1910.11976) (double precision loses a spatial symmetry at t ≈ 580) and [arXiv:2505.18696](https://arxiv.org/pdf/2505.18696) (double against quadruple precision, qualitatively different routes).
- [snippet] CESTAC/CADNA (La Porte and Vignes 1974; Vignes, Math. Comput. Simul. 35, 1993) ([Jézéquel slides](https://perso.lip6.fr/Fabienne.Jezequel/ARTICLES/CADNA_pres.pdf)): random rounding as a principled draw of the rounding, an alternative to summation-order permutations.
- [snippet, round 133] Brouwer's law as described in arXiv:1711.06050; Mosbach–Turner (arXiv:math/0512364).
- [memory] DerSimonian–Laird random effects; Filippov sliding at a discontinuity of a vector field (the traps' exact-flow analogue).

## Lesson

A scaling law measured from two arithmetics is one equation in two unknowns. And a rounding "seed" need not be noise at all: where the computed field has ulp-scale zeros, the integrator decides whether there is an exit. Before reading a two-precision ratio as an exponent or a prefactor, add a third precision, and look at the computed field at the scale of the unit roundoff.
