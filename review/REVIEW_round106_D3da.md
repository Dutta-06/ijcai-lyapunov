# Referee report D3da, round 106: Proposition D‴.59 as revised after round 103

## Scope

- **Export and work directory.** I exported commit 5ee87f7 (branch claude/vigilant-clarke-hf939b) with `git archive 5ee87f7 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref106. W/repo/venv is a symlink to the repository's venv (numpy 2.5.4, scipy 1.18.1, pytest 9.1.1).
- **Diffs.** They run from 7b92874 to 5ee87f7. The raw is unchanged over that range (moving_width_raw.npz is not in the diff). The md changed at the symmetric-line paragraph, the held-run line, the new exit and lemma-scan lines, and the splits' mirror sentence.
- **Read in full:**
  - review/drafts/D59.md (395 lines);
  - every changed part of experiments/exp_moving_width.py (`a2_ulps`, `a2_exit`, `a2_lemma_scan`, `mirror_asym`, the changed `a2_lines` lines, `split_summary`'s mirror block, `outcome`), and the code they call (`a2_state`, `a2_fixed_point`, `a2_lemma`, `a2_invariance`, `critical_rhs_nd`/`critical_jac_nd`, `BumpsCW`, `null_cluster`, `match`, `h_job`/`hl_job`, `sliding_iss_nd`'s ties, `polish_block_nd`);
  - results/chaos/moving_width.md; tests/test_moving_width.py (52 tests); tools/mutants/round103_D59_own.py; tools/mutants.py;
  - the Round 103 entry of review/REVIEW_LOG.md; review/REVIEW_round103_D3cx.md; the D59 row of review/drafts/README.md.
- **Settings.**
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, under `timeout`, with MW_OUT=W/out (a copy of the raw) or a mutant slot's `out`.
  - pytest ran with `--basetemp` inside W. My copy of the mutant tool, W/own/mutants_tool.py, sets REPO to W/repo and adds `--basetemp` per slot.
  - At most two heavy processes of mine ran at a time.
  - I stopped two of my own processes, PIDs 18229 and 20590: method-comparison scripts stuck at atol 10⁻²⁵ and 10⁻²⁰ on symmetric samples. I checked /proc/PID/cmdline before stopping each. No wait loop matched its own command line.
- **The repository is untouched.** I used only read-only git (archive, diff, log, status). `git diff 5ee87f7` on the D59 files is empty at the end. Another session committed 04bca94 (D‴.58) meanwhile; it touches no D59 file.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my computation; [snippet] a web-search result; [memory] not verified this round.

## Summary

**Design.**
1. I regenerated the md from the raw, ran the tests, and reran jobs against the raw.
2. I located the zero of ċ along c = c* + k·ulp at the two held states, in double and in 50-digit decimal arithmetic on the double data.
3. I integrated the flow from c* with DOP853, RK45, LSODA, scipy BDF and Radau at several tolerances, and with a fixed-step RK4 in float64 and in x87 long double:
   - on A's samples and on exactly symmetric samples;
   - from several seeds;
   - with the Jacobian's spectrum along the route.
4. I recomputed Lemma D‴.59a's terms by finite differences, compared the truncated expansion with the full sum, and scanned σ/h on a fine grid.
5. I recomputed `mirror_asym` per target and compared it with:
   - a bijective (bottleneck) matching;
   - the same measure on all active kernels;
   - the fitted function's odd part;
   - H3 without polish;
   - HL's loops in their whole state.
6. I ran mutants: the repository's 15, 25 of my own, and 3 job-level ones, each set with a control.
7. I ran web searches.

**Verdicts.**
- **Reproducibility is exact.** The md matches byte for byte. All 52 tests pass. All 13 rerun jobs match the raw bit for bit.
- **Every number the revision added is the md's verbatim** and is reproduced by my own code.
- **"The location the problem's" holds, to the ulp, wherever a solver holds the state.** The 50-digit zeros of ċ are 8.12 ulps (5A) and 9.57 ulps (3A). LSODA holds at 8 and 9 (alternating 9 and 10). scipy BDF holds 3A at 10.
- **"The stability the solver's" is narrower than written.** It is LSODA's step sequence (and scipy BDF's at rtol 10⁻¹⁰ on 3A), not implicit methods in general. A small-step explicit method also holds the state in double, by absorption.
- **The exit claim is misattributed.** On A's samples the exit is seeded by the problem's own data, not by roundoff.
- **The lemma's expansion is right**, with small slips. The weight "asymptotic … numerical at A2E's states" is right.
- **The mirror sentence overgeneralises**, and its measure is not a distance to the fixed-point subspace.
- **The tests are strong.** An md-from-raw equality test kills every printing mutant. All 15 repository mutants and all 3 job-level mutants are killed. Of my 25, 21 are killed and 4 survive; the 4 are equivalent on the data. No survivor changes a published result.
- **No major issue.**

**Findings.**
1. **D3da-1 (minor; the principal finding): the exit from c* is seeded by the data on A's samples, by the arithmetic only on symmetric samples.**
   - "Roundoff seeding both", "alike", and the md's "in double the exit is seeded by roundoff alike for the symmetric problem and the posed one" do not hold.
   - On A's samples the exit time does not move with precision, integrator or tolerance. Every exit goes the same way, downward (c* lies below the repelling state).
   - On exactly symmetric samples the exit time moves with precision, tolerance, integrator and summation; its direction flips; with some integrators it never comes.
   - In exact arithmetic the symmetric problem never exits. Round 103's exact-arithmetic estimate for the exact-linspace problem (s ≈ 0.26, 4.2) was not carried over.
2. **D3da-2 (minor): the hold.**
   - Correct for LSODA's runs.
   - Not true of implicit methods in general: scipy BDF and Radau leave.
   - Not exclusive to implicit methods: a fixed-step explicit RK4 in double holds c bit-exactly at c*, 8 ulps off the zero of ċ, because its increments are below half an ulp.
   - So the md's "a solver stays on it only by damping" and "(DOP853, which neither damps nor holds …)" generalise beyond the evidence.
3. **D3da-3 (minor): the mirror sentence.**
   - Only 2 of the 7 symmetric polished ends leave the fixed-point subspace.
   - H3's 0.3 is a coinciding pair's gauge; the unpolished H3 ±0.15 path leaves too (0.048).
   - "The unpolished path (1.4·10⁻⁷) and the loops stay on it" holds for their formed kernels only. HP's whole state is 1.2·10⁻² off; HL σ* 0.04 at 10⁻⁸ is 1.1·10⁻³ off and growing.
   - `mirror_asym` is a nearest-image lower bound and vanishes off the subspace.
4. **D3da-4 (cosmetic): the lemma's wording.**
   - A dimensional slip in the overshoot condition.
   - The remainder is O((h/σ)²φ₂), not O(φ₂).
   - At the least-squares w the overshoot is of the same order as the neighbours' term, not "exponentially small" relative to it.
   - "Carries the sign of w·r₁" holds only when |r₁| ≫ |w|φ₁.
   - The sign boundaries sit on the scan's grid.
5. **D3da-5 (cosmetic):** the "zero of ċ at 9.9 ulps" is the double sum's (in 50 digits, 9.57), and the md's "zero of the linearised ċ" is −d/λ = 8.72 ulps.
6. **D3da-6 (cosmetic): stale or loose text.**
   - The md's C-outcome sentence still says "a null cluster where it fits none".
   - The symmetric-line paragraph still calls the state "a saddle in c" unconditionally.
   - The Fenichel sentence applies to a plane that is not normally hyperbolic along the route.
   - "Exactly symmetric pair": XS is symmetric only to 4.4·10⁻¹⁶.
7. **D3da-7 (cosmetic): tests.**
   - Four equivalent survivors.
   - The 3A held-offset assertion passes with a 0.93-ulp margin against a 1-ulp tolerance.
   - `a2_lines`' exact-start sentence is hard-coded and would contradict itself on a run that leaves.

## What I checked and how

Scripts are in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref106/own/ and logs in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref106/logs/. Below they are named by basename only.

1. **Reproducibility** [own].
   - `--from-raw` into W/out took 2.5 s; `cmp` with the committed md: identical (sha256 cff5cdfa…).
   - pytest: 52 passed in 17.7 s (W/pt/tests.log).
   - rerun.py compares every field except `secs`, bit for bit. Identical:
     - `part_A2E` (all 3 entries);
     - `a2_job` (0.25, 3A), (0.25, 5A), (1, 10A), (4, 20A);
     - `h_job` H σ* 0.06 equal;
     - H3 ±0.12 and ±0.15;
     - HP σ* 0.06 equal (151 s);
     - `hl_job` HL σ* 0.06 equal.
2. **Zeros of ċ** (ulps.py, ulps.log) [own].
   - XS[163] − C_A = 2.776·10⁻¹⁷ = 8 ulps of 0.03; 1 ulp/σ* = 1.156·10⁻¹⁶.
   - XS is symmetric about 0 to 4.4·10⁻¹⁶, C0 to 2.2·10⁻¹⁶. The sample pairs about C_A, or about XS[163], are off by up to 4.4·10⁻¹⁶.

   | held state (A2E, s = 10) | σ/h | λ | held | double zero (`a2_ulps`) | −d/λ (double) | 50-digit zero on the double data (f* as stored / f* exact) |
   |---|---|---|---|---|---|---|
   | 5A, span 10 | 0.554 | 1.0498·10⁴ | 8 ulps | 8.118 | 8.104 | 8.118 / 8.124 |
   | 3A, span 10¹³ | 0.706 | 628.85 | 9 ulps (last 15 BDF steps 9/10) | 9.932 | 8.717 | 9.570 / 9.643 |

   - At 3A the double ċ is noisy at the sub-ulp level. Its per-ulp increments are 2.1·10⁻¹⁴ … 1.3·10⁻¹³ against λ·ulp = 7.3·10⁻¹⁴. At 5A they are 1.17–1.24·10⁻¹² against 1.21·10⁻¹².
3. **Exit and hold from c* exactly** (exits.py, methods2.py, methods.log, methods2.log) [own]. Values are the s at which |c − c₀| ≥ 10⁻²σ*; the sign is the exit's direction.

   | | A's samples | exactly symmetric samples |
   |---|---|---|
   | 5A, DOP853, rtol 10⁻⁸ / 10⁻¹⁰ / 10⁻¹² / 10⁻¹³ | 0.00351 / 0.00348 / 0.00347 / 0.00350, all − | 0.00364 (−) / 0.00397 (+) / 0.00377 (+) |
   | 5A, RK45, rtol 10⁻¹⁰ / 10⁻¹² | 0.00352 / 0.00371 (−) | 0.00391 (−) |
   | 5A, LSODA, span 10 | **held** (+9.252·10⁻¹⁶σ*) | **0.00350 (+)** |
   | 5A, LSODA, span 10¹³ (no event) | 0.00352 | — |
   | 5A, scipy BDF, rtol 10⁻¹⁰ / 10⁻¹² | 0.00390 / 0.00418 (−) | 0.00395 (+) |
   | 5A, Radau, rtol 10⁻¹⁰ / 10⁻¹² | 0.00390 / 0.00424 (−) | 0.00392 (+) |
   | 3A, DOP853, rtol 10⁻⁸ / 10⁻¹⁰ / 10⁻¹² | 0.04807 / 0.04828 / 0.04799 (−) | 0.05616 / 0.05072 (+) |
   | 3A, RK45, rtol 10⁻¹⁰ | 0.04949 (−) | 0.06285 (−) |
   | 3A, LSODA, span 10 | **held** (9 ulps) | **held** (+9.7·10⁻¹⁷σ*) |
   | 3A, LSODA, span 10¹³ | **held** to the exact fit (max 2.2·10⁻¹⁰σ*) | — |
   | 3A, scipy BDF, rtol 10⁻¹⁰ / 10⁻¹² | **held** (10 ulps) / 0.05114 (−) | **held** (−8.4·10⁻¹⁸σ*) |
   | 3A, Radau, rtol 10⁻¹⁰ / 10⁻¹² | 0.05572 / 0.05330 (−) | 0.07530 (+) |

   - On A's samples all 18 exits are downward. On symmetric samples both directions occur.
4. **Precision and seeds** (ldexit.py, lddiag.py, lddiag2.py, the seed one-liner; ldexit.log, lddiag*.log, seeds.log) [own].
   - Fixed-step RK4: h = 2·10⁻⁶ at 5A and 10⁻⁵ at 3A; halving h changes nothing.
   - **On A's double samples, long double gives 0.00348 (5A) and 0.04830 (3A).** With f* recomputed in long double: 0.00348 and 0.04828. These are DOP853's double values: the seed is the data.
   - **The same RK4 in float64 never leaves.** c stays bit-exact at C_A to s = 0.006, because its increments (6.0·10⁻¹⁹) are below half an ulp of c* (1.73·10⁻¹⁸).
   - **On symmetric samples, 5A moves with the precision:** 0.00401 in float64, 0.00465 in long double. The shift 6.4·10⁻⁴ matches the predicted ln 2¹¹/λ = 7.0·10⁻⁴.
   - **On symmetric samples, 3A never exits by s = 0.3 in either precision.** c is held at −8.8·10⁻¹⁹ in float64 (½ ulp of the ±1 samples) and at −4.2·10⁻²² in long double.
   - **Symmetric samples started at the posed problem's offset:** −8.1 ulp(0.03) gives 0.00355 (downward), +8.1 ulp gives 0.00357 (upward), −6 ulp gives 0.00360; 3A at −9.6 ulp gives 0.0495. These reproduce A's 0.00348 and 0.0483.
   - Exact arithmetic, from this:
     - symmetric problem: never;
     - the computed problem (double XS): ≈ 0.0035 and 0.048;
     - the exact-linspace problem: ≈ ln(10^{1201})/λ ≈ 0.25 and 4.2, round 103's estimate, which my order count reproduces.
5. **Spectrum along the symmetric route** (exits.py `spectrum`; spectrum.log) [own].
   - 5A: ∂ċ/∂c goes from −356 … −380 (σ/h 2.96 → 1.84) to +1.07–1.09·10⁴ (σ/h 0.55). The tangential (α, σ) eigenvalues there are −0.0012 and −1.48·10⁴.
   - 3A: from −212 … −93 (σ/h 2.99 → 0.88) to +665 (0.70).
   - The decoupling holds: |J[1,0]| and |J[1,2]| are ≤ 5·10⁻¹³ on the line.
6. **The lemma** (lemma.py, lemma.log) [own].
   - Central differences of `critical_rhs_nd` equal `a2_lemma_scan`'s totals to six digits at all 12 states. At A2E's four line states, the finite differences, the Jacobian and the split agree to 8 digits.
   - The 3-term expansion against the exact total, at w = 4A:

     | σ/h | 1 | 0.8 | 0.6 | 0.45 | 0.3 |
     |---|---|---|---|---|---|
     | expansion | +54.5 | −113.9 | +4180.8 | +14753.7 | +32659.7 |
     | exact | −213.7 | +193.7 | +4340.7 | +14761 | +32659.7 |

     So the truncated expansion has the wrong sign at σ/h = 1 and 0.8.
   - The remainder over (Gm w/(Nσ²))·φ₂(h/σ)² is −1.7 … +4.9 for σ/h ≥ 0.3.
   - At the least-squares w, the overshoot w − A is 7.284·10⁻³ against 2φ₁f*(x₁) = 7.314·10⁻³ (σ/h 0.3), and 7.05·10⁻⁶ for both at 0.2.
   - The ratio (±1 term)/(centre term) is 10.07 against (h/σ)² − 1 = 10.11 at 0.3, 24.0 against 24 at 0.2, and 43.44 against 43.44 at 0.15.
   - At σ/h 0.8 the ±1 term is −369 although r₁ = +0.26.
   - Fine grid (step 0.01):
     - w = 4A: the total is positive below σ/h 0.85, and the ±1 term positive below 0.47;
     - least-squares w: the total is positive below 0.78, the ±1 term positive below 0.62, and it leads the centre only below 0.56.
   - At 3A's held state the split is +2840.7, −2465.1, +249.7 and +3.6 (total 628.9): the ±2 term is 40% of the total.
7. **Mirror symmetry** (mirror.py, h3np.py, hpdump.py, loopsym.py; mirror.log, h3np.log, loopsym.log) [own].

   | symmetric end | `mirror_asym`, formed | all active | bottleneck matching | the fit's odd part / max f* |
   |---|---|---|---|---|
   | H σ* 0.04 (polished) | 2.2·10⁻¹¹ | 2.6·10⁻¹⁰ | 2.2·10⁻¹¹ | 3.5·10⁻¹⁷ |
   | H σ* 0.06 (polished) | **0.224** | 0.419 | **0.240** | 6.1·10⁻¹⁰ |
   | H2 σ* 0.06 | 0.218 | 0.32 | 0.231 | 1.1·10⁻¹⁰ |
   | H3 ±0.05 | 1.7·10⁻³ | 0.096 | 1.7·10⁻³ | 6.3·10⁻⁷ |
   | H3 ±0.1, ±0.12, ±0.2 | 7·10⁻¹⁸, 5·10⁻¹³, 1.5·10⁻⁶ | — | same | — |
   | H3 ±0.15 (polished: 0.378 + 0.622 on −0.15, 1.0 on +0.15) | **0.30** | 0.30 | 0.30 | **7.3·10⁻¹⁶** |
   | H3 ±0.15 *unpolished* (shares 0.272/0.728 against 0.776/0.224) | **0.048** | 0.048 | — | — |
   | HP σ* 0.04 / 0.06 | 1.4·10⁻⁷ / 2.9·10⁻⁸ | **1.2·10⁻²** / 3.8·10⁻⁵ | same as `mirror_asym` | 1.4·10⁻⁷ / 8.9·10⁻¹³ |

   - **Counterexample:** `mirror_asym({a, a, Ra}) = 0`, while the bottleneck matching distance is 0.2.
   - **HL rerun with all 21 kernels** (i ↔ 20 − i, α included):
     - σ* 0.06: the whole state equals the formed value (≤ 1.1·10⁻⁸);
     - σ* 0.04 at ε = 10⁻⁸: the whole state is 2.3·10⁻⁴ → 6.4·10⁻⁴ → 1.06·10⁻³ off (τ 6400 → 10⁵), against formed values of 2.8·10⁻⁷ → 9.2·10⁻⁸;
     - σ* 0.04 at 10⁻⁶: 1.6·10⁻⁶ off, against 6·10⁻¹⁰ formed.
   - HLL (formed kernels) stays bounded: 1.1 → 1.7·10⁻⁸ over τ 10⁵ → 10⁹ at 10⁻⁸.
8. **Numbers.** Every number added to D59.md since 7b92874 is the md's verbatim [md]: 8; 8.1, 9.9; 0.0483, 0.0035, 0.0507, 0.0040; +2871, −2817, −213.7, +575.8, +196.3, +19.5; 0.22, 0.3, 1.4·10⁻⁷, 2.8·10⁻⁷. My own code reproduces each [own].
9. **Mutants** (W/mut/res_{r103,z106,j106}.jsonl; mut_*.log).
   - Repository spec: 15 mutants, all killed. Y10, Y13 and Y14 are killed by the tests alone; the other 12 by the tests and the md (Y11 also crashes `--from-raw`).
   - spec106.py (25): 21 killed. Survivors: Z02 (kmax 12), Z03 (only rising sign changes), Z15 (signed σ) and Z23 (the held ulps printed from −d/λ, which round to the same 8 and 9). All four are equivalent on the data.
   - spec106b.py, job level: J01 (A2E_MORE's 3A from c* + 10⁻¹²σ*), J02 (`hl_job` threshold 10⁻³) and J03 (`a2e_job` rtol 10⁻⁹): all 3 killed.
   - All controls pass.
   - A job-level rerun of J01's input, outside the tool, shows the md's exact-start sentence would read "|c − c*|/σ* at most 1.6e-01 … on the plateau from s = 0.056 — held on the symmetric state by LSODA's damping". The wording is hard-coded.
10. **Literature.** Seven web searches. The one fetch I tried (arxiv.org) failed: the host is not in this session's allowed network.

## Verdicts per part

| claim (D59.md) | verdict | evidence |
|---|---|---|
| (i) continuum landscape and convergence (8–15, 144–152) | **holds**; unchanged | re-read |
| (i) "in double the centre sample lies 8 ulps of c* above c*" (20) | **holds** | [own] |
| (i) "for widths well below the spacing with an O(1) overshoot … a saddle in c" (21–22) | **holds asymptotically** | [own] fd, expansion |
| (i) "a repelling near-symmetric state where ċ vanishes, 8–10 ulps above c* at A2E's states" (24) | **holds**: 50-digit zeros at 8.12 and 9.57 | [own] |
| (i) "leaves it … at s = 0.0483 and 0.0035, alike on exactly symmetric samples, roundoff seeding both" (25–26) | **times hold; the seed attribution does not** (D3da-1) | [own] items 3–4 |
| (i) "LSODA's run … stays on it because BDF's steps beyond their bounds make it attracting: the location the problem's, the stability the solver's" (26–27) | **location holds** to the ulp; **stability is LSODA's step sequence**, not implicit methods' (D3da-2) | [own] items 2–4 |
| Lemma D‴.59a: equivariance, invariance, decoupling, ∂ċ/∂c formula (152–156) | **proved** (standard); checked | [own] |
| Lemma: the expansion and its saddle conditions (156–161) | **right in substance**; remainder O((h/σ)²φ₂), wording slips (D3da-4) | [own] item 6 |
| "the saddle's sign asymptotic under these hypotheses and numerical at A2E's states" (161–162; Status 212) | **right weight** | — |
| multi-kernel invariance on Fix(reflection ∘ permutation); citation (162–165) | **holds** for flows (the loops confirm it in σ* 0.06's whole state); citation adequate | [own], [snippet] |
| Fenichel persistence of "the repelling line" (165–166) | **loose** (D3da-6) | [own] item 5, [snippet] |
| A2E numerics: held-run, exit and scan lines (282–292) | **numbers hold**; "zero at 9.9" is the double sum's (D3da-5) | [md], [own] |
| (ii)(b) "the polished ends leave the fixed-point subspace … while the unpolished path and the loops stay on it" (71–74) | **overgeneralised; measure inadequate** (D3da-3) | [own] item 7 |
| `outcome`'s label 'coinciding pair off the features' | **applied**; the md prose is stale (D3da-6) | [md] line 89 |
| tests | **strong**; 4 equivalent survivors (D3da-7) | [own] item 9 |

## Status of the round-103 issues

| D3cx | status | what remains |
|---|---|---|
| 1 (the held state the problem's) | **applied** (`a2_ulps`, the wording; round 100's phrase withdrawn) | Holds to the ulp. The "zero at 9.9" is the double sum's. "The stability the solver's" is LSODA's (D3da-2, D3da-5). |
| 2 (exit unquantified) | **partly applied**: the DOP853 times are printed | The seed reading is wrong for A's samples (D3da-1). Round 103's exact-arithmetic sentence (≈ 0.26, 4.2) was dropped without note. |
| 3 (the lemma's saddle sign) | **applied** (expansion, `a2_lemma_scan`, Status) | Slips taken over from the referee's own proposed formula (D3da-4). |
| 4 (one null predicate) | **applied** (the rename; the tests follow) | The md's C sentence still says "a null cluster where it fits none" (D3da-6). |
| 5 (wording) | **applied** | — |
| 6 (scope and literature) | **applied** (several kernels, citations, mirror line, the tie) | The mirror claim overreaches (D3da-3); Fenichel loose (D3da-6). |
| 7 (tests) | **applied**: W26, W29, W30, W31, W36 planted; 15 own mutants, all killed [own rerun] | Equivalent survivors; the 3A ulps margin; hard-coded exact-start text (D3da-7). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3da-1 | minor | **The exit's seed on A's samples is the data's, not roundoff.** (D59.md:25–26 "alike on exactly symmetric samples, roundoff seeding both"; :288 "in double roundoff seeds the exit alike"; md exit line; README.) <br>• A's samples: long-double RK4 gives 0.00348 and 0.04830 (DOP853 in double: 0.00348, 0.0483). DOP853 over rtol 10⁻⁸ … 10⁻¹³ gives 0.00347–0.00351 and 0.04799–0.04828. All 18 exits over 5 integrators are downward [own]. <br>• Symmetric samples: 5A moves 0.00401 → 0.00465 with precision (ln 2¹¹/λ predicted). DOP853 gives 0.00364/0.00397/0.00377 with flipping sign. 3A spans 0.0507–0.0753, or never: LSODA, BDF and RK4 hold it to s = 0.3–10 [own]. <br>• Symmetric samples started at −8.1 ulp(0.03) reproduce A's exit (0.00355) [own]. <br>• The seeds differ by ~e^{λΔs} ≈ 200 (5A); the times are "alike" only because the seeds enter logarithmically. <br>• In exact arithmetic: the symmetric problem never exits; the computed problem ≈ the double times; exact linspace ≈ 0.25 and 4.2. | **D59 (i) and the A2E numerics:** "on A's samples the exit is the data's — c* lies 8 ulps below the problem's repelling near-symmetric state: the same time in long double and for every explicit tolerance, always downward —; on exactly symmetric samples it is the arithmetic's (it moves with precision, tolerance and summation, and need not occur); in exact arithmetic the symmetric problem stays, the computed one leaves at these times, the exact-linspace one at ≈ 0.26 (5A) and 4.2 (3A), an estimate". <br>**Code (`a2_exit` line):** also print the exit's direction, and the symmetric exit at a second tolerance. <br>**Test:** the A's-samples exit is downward and tolerance-stable. |
| D3da-2 | minor | **The hold is LSODA's step sequence, neither generic to implicit methods nor exclusive to them.** (D59.md:26–27, :285–286; md's symmetric-line "a solver stays on it only by damping the unstable direction"; md exit line "(DOP853, which neither damps nor holds …)".) <br>• scipy BDF and Radau leave 5A (0.0039–0.0042); Radau leaves 3A (0.0533–0.0557); BDF at rtol 10⁻¹² leaves 3A (0.0511). <br>• LSODA itself leaves 5A on symmetric samples (0.00350) and on span 10¹³ (0.00352). <br>• scipy BDF at rtol 10⁻¹⁰ holds 3A at 10 ulps, LSODA at 9; both bracket the exact zero 9.57, which confirms the location. <br>• A fixed-step RK4 in double holds c = c* bit-exactly (increments 6·10⁻¹⁹ < ½ ulp 1.7·10⁻¹⁸). That held state is 8 ulps off the zero: an explicit method holding by absorption [own]. | "LSODA's steps (BDF beyond its bounds after its switch) make it attracting; scipy's BDF and Radau leave at rtol 10⁻¹⁰–10⁻¹² (5A); every implicit hold sits within an ulp of the zero of ċ". Replace "only by damping" with "by damping (LSODA), or by absorption in a small-step explicit method in double". Drop "neither … holds" or qualify it by step size. |
| D3da-3 | minor | **The mirror sentence overgeneralises, and `mirror_asym` is not a distance to Fix.** (D59.md:71–74; md splits line; `mirror_asym`.) <br>• Polished symmetric ends: H σ* 0.04 2.2·10⁻¹¹; H σ* 0.06 0.224; H3 ±0.05 1.7·10⁻³; ±0.1, ±0.12, ±0.2 ≤ 1.5·10⁻⁶; ±0.15 0.30. Only 2 of 7 leave. <br>• H3's 0.3 is a coinciding pair's gauge with a kernel count of 3; its fitted function is symmetric to 7·10⁻¹⁶. The unpolished H3 ±0.15 also leaves (0.048, gauge shares). <br>• HP "1.4·10⁻⁷" is its formed kernels; all 21 active are 1.2·10⁻² off. HL σ* 0.04 at 10⁻⁸ is 1.1·10⁻³ off in its whole state and growing. <br>• By the gauge-invariant odd part, HP σ* 0.04 (1.4·10⁻⁷) is less symmetric than the polished ends (6·10⁻¹⁰, 7·10⁻¹⁶). <br>• `mirror_asym` uses the nearest image, not a bijection: {a, a, Ra} gives 0 against a bottleneck distance of 0.2, and H gives 0.224 against 0.240 (round 103's 0.24). | Use a bottleneck matching (e.g. `maximum_bipartite_matching` on thresholds) in `mirror_asym`. Print per-target values, the all-active value and the fit's odd part. **D59 (ii)(b):** "on the equal σ* 0.06 pair the polished end's cancellation cluster sits on one feature (its kernels 0.24 from their mirror images, its function symmetric to 6·10⁻¹⁰): a symmetry breaking of the decomposition; the unpolished path's formed kernels and the loops' formed kernels stay within 3·10⁻⁷ (the HP end's weak kernels 10⁻² off, HL σ* 0.04's 10⁻³)". Drop H3 ±0.15 or label it a gauge. |
| D3da-4 | cosmetic | **The lemma's wording.** (D59.md:156–161.) <br>• "w(w − f*(c*)) ≫ (h/σ)²φ₁\|r₁\|" mixes w² with w. It should read (w − f*(c*)) ≫ (h/σ)²φ₁\|r₁\| with w > 0. This came from round 103's proposed text. <br>• The remainder is O((h/σ)²φ₂): the truncated sum is +54.5 against −213.7 at σ/h 1 and −113.9 against +193.7 at 0.8. <br>• At the least-squares w the overshoot is 2φ₁f*(x₁)(1 + O(φ₁)), the same order as the neighbours' term, which leads by the factor (h/σ)² − 1; the saddle then holds since the total is ≈ 2φ₁f*(x₁)(h/σ)² > 0. <br>• "The neighbours' term carries the sign of w·r₁" only for \|r₁\| ≫ \|w\|φ₁ (−369 at σ/h 0.8 with r₁ = +0.26). <br>• Fine-grid boundaries: 0.85 and 0.47 (w = 4A); 0.78, 0.62 and 0.56 (least-squares w). <br>• "(A2E: …)" labels `a2_lemma_scan`'s hypothetical states. | Rewrite as stated. Write "+ O((h/σ)²φ₂)". Say "the neighbours lead by (h/σ)² − 1". Write "(the scan, not A2E's states)". |
| D3da-5 | cosmetic | **Ulps precision.** (md held line; D59.md:284.) <br>• "The zero of ċ at 9.9 ulps" is the double sum's (per-ulp noise ±0.6 ulp); in 50 digits on the same data it is 9.57 (9.64 with f* exact). <br>• "The offset sits at the zero of the linearised ċ" names −d/λ = 8.72 ulps, a different zero. <br>• The held "9" is the last step of a 9/10 alternation. | "Held at 8 and 9–10 ulps, the zero of ċ at 8.12 and 9.57 ulps (exact arithmetic on the double data; 8.1 and 9.9 in double)". |
| D3da-6 | cosmetic | **Stale or loose text.** <br>• md:89 (`report`'s C sentence) still says "a null cluster where it fits none". <br>• The md's symmetric-line paragraph says "and the symmetric state a saddle in c" unconditionally. <br>• Fenichel: the invariant plane c = c* has a normal rate that changes sign along the route (−93 → +665 at 3A; −380 → +1.09·10⁴ at 5A). The classical result is for compact normally hyperbolic manifolds; the noncompact case needs Eldering's bounded geometry [snippet]. So persistence holds on its strictly repelling pieces, and "the state where ċ vanishes" is the zero set of one component. <br>• "The exactly symmetric pair": XS is symmetric to 4.4·10⁻¹⁶ and C0 to 2.2·10⁻¹⁶. | Fix the C sentence's wording. Write "a saddle in c for σ well below h (Lemma)". Write "Fenichel, on compact pieces where it is strictly repelling". Write "symmetric to 4·10⁻¹⁶ in double". |
| D3da-7 | cosmetic | **Tests.** <br>• Survivors Z02, Z03, Z15 and Z23 are equivalent on the data. <br>• `abs(kz − held) < 1` passes 3A at 0.93. <br>• `a2_lines`' exact-start and held texts are hard-coded: a run that leaves would print "… on the plateau from s = 0.056 — held on the symmetric state by LSODA's damping" [own rerun of J01's input]. A2E_MORE's pin and the raw assertions make this latent. | Condition the exact-start and held sentences on the data. Test the held offsets against the 50-digit zeros within 0.6 ulp. |

## Sources

- **[snippet]** Neuberger, Sieben and Swift, "A Bifurcation Lemma for Invariant Subspaces", arXiv:2308.10448 ([ar5iv](https://ar5iv.labs.arxiv.org/html/2308.10448); SIAM J. Appl. Dyn. Syst. 2024, [DOI 10.1137/23M1595540](https://epubs.siam.org/doi/10.1137/23M1595540)). It describes itself as generalising the equivariant branching lemma from fixed-point subspaces to invariant subspaces.
- **[snippet]** Golubitsky, Stewart and Schaeffer, *Singularities and Groups in Bifurcation Theory II*, Springer AMS 69, 1988 ([Springer](https://link.springer.com/book/9780387966526)). Only the bibliography surfaced. The invariance of fixed-point subspaces is [memory].
- **[snippet]** Fenichel, Indiana Univ. Math. J. 21 (1971) 193–225 ([citation listing](https://modeldb.science/citations/81561)). Eldering, "Persistence of noncompact normally hyperbolic invariant manifolds in bounded geometry" ([arXiv:1205.2753](https://arxiv.org/pdf/1205.2753), [arXiv:1204.1310](https://arxiv.org/pdf/1204.1310)): it generalises Fenichel and Hirsch–Pugh–Shub from compact to noncompact manifolds.
- **[snippet]** Compensated summation for small-step increments in RK schemes ([arXiv:1702.03354](https://arxiv.org/pdf/1702.03354)): "for small step lengths, the increments are smaller than the solution components". This is related to the absorption hold.
- **[snippet]** Backward Euler's factor 1/(1 − hλ) for λ > 0, from a search summary of course notes. I did not attribute it to a specific page. **[memory]** BDF-k's bounds Σ_{j≤k} 2^j/j, as printed.
- **Not reachable.** arxiv.org could not be fetched: it has to be allowed in the session's network settings.

## Lesson

A seed that survives a change of arithmetic belongs to the problem. One that moves with the precision, the tolerance or the summation order belongs to the computation. Before calling an exit "roundoff's", rerun it in long double and with a second integrator, and record its direction:
- here A's exit stayed put and always went down, which marks the data's 8-ulp offset;
- the symmetric problem's exit wandered, flipped sign, or never came.

A held state likewise needs more than one solver:
- LSODA's hold is its step sequence, not a property of implicit methods;
- a small-step explicit method holds by absorption, at a point that is not the zero of ċ.

A fixed-point-subspace claim needs three things:
- a bijective matching;
- all of the state, not only the formed kernels;
- a gauge-invariant check (the fitted function's odd part), because a coinciding pair is asymmetric in its parameters while its function is symmetric to 10⁻¹⁵.
