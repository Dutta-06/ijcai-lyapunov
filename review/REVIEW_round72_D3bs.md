# Referee report D3bs, round 72: Proposition D‴.54 as revised after round 66

## Scope

**What I reviewed.** Commit 37918e1 on `claude/festive-bardeen-eopej2` ("D‴.54 round 66: s_d the climb's end plus the hold; the continuum about a merged block"). I exported it with `git archive 37918e1 | tar -x` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref72/D3bs_work/repo` and reviewed:
- **`review/drafts/D54.md`:**
  - (iv) l. 61–89 and (v) l. 91–159;
  - the proof's relative-mode paragraph, l. 200–211;
  - the Status, Classical and New paragraphs, l. 214–266;
  - Numerics: C3 (l. 323), C4 (l. 329), C4W (l. 340), C4M (l. 343), R (l. 351), RP (l. 374) and RC (l. 397–410).
- **`experiments/exp_transit_nd.py`:** `c4_atol`, `c4_job`, `race_job`, `PairSums` (`__call__`, `moments`, `block`, `force`), `relative_law`, `rp_job`, `growth_rate`, `merged_rhs`, `rc_job`, `part_RC`, `rp_lines`, `rc_lines`, `r_lines`, `main` / `--from-raw`.
- **Results and tests:** `results/chaos/transit_nd.md`, `results/chaos/transit_nd_raw.npz` and `tests/test_transit_nd.py` (49 tests).
- **Imported code:** `src/core/rich_transit_nd.py` (`BumpsGrid2`, `critical_rhs_nd`, `critical_jac_nd`, `loop_nd`) and `src/core/rich_transit.py` (`critical_w`, `curve_mobility`, `transit_scale`, `eta_crit`, `_alpha_of`).
- **History:** `review/REVIEW_round66_D3bm.md`, the "Round 66" entry of `review/REVIEW_LOG.md`, and `git diff 6d48964 37918e1`.

**Settings.**
- I used the repository's `venv/bin/python` (Python 3.13.16, numpy 2.5.4, scipy 1.18.1) with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`.
- `exp_transit_nd.py` writes to `TN_OUT` (its default is `results/chaos`). I set `TN_OUT` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref72/D3bs_work/out`, which holds a copy of the raw.
- Everything is under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref72/D3bs_work`: scripts in `/scripts`, logs in `/logs`, outputs in `/out`, mutant workspaces in `/mut/ws1`–`/mut/ws3`.
- At most two of my heavy processes ran at once. I once started a third, my own `thresh.py` (PID 2822). I checked `/proc/2822/cmdline`, killed it and reran it later.
- Nothing of mine is still running. The machine was shared, with load averages up to 9.5.

**Repository untouched.** In `/home/user/Lyapunov-Convergence-Kernels-Showcase` I ran only read-only git commands (`archive`, `log`, `diff`, `show --stat`, `status`, `rev-parse`). The live HEAD has since moved to 18bcfa0, but `git diff 37918e1 HEAD` is empty for every D54 file.

**Literature.**
- WebSearch worked: 11 queries, results read as snippets only.
- WebFetch could not resolve `jmlr2020.csail.mit.edu` or `www.jmlr.org`, so I stopped fetching.
- I did not use arxiv.org or ams.org.
- Anything resting on snippets or memory is marked as such.

## Summary

**What checks out**

*Reproducibility*
- `--from-raw` regenerates the md byte for byte (md5 8b1748c8ea9c175dc7c2615328581e63; own run).
- All 49 tests pass, in 55.6 s (own run).
- I reran a sample of jobs with the experiment's own functions. Each equals the raw in every field except `secs` (own run):
  - `rc_job` ('near tie', 3) and ('near tie mirrored', 4.9);
  - the RC check rows ('near tie', 4.9, 4, 'zero', 10⁻⁹) and ('near tie', 4.9, 16, 'y', 10⁻⁹);
  - `rp_job` at x = 3, on Part C's window and on the far window;
  - C4 at (4.9, 10⁻¹¹): captured, toward, departing at s = 3.297·10⁵, atol 10⁻¹⁶;
  - C4M at (7, 10⁻¹³) on the early grid: bound, far.

*RC (the continuum's relative mode) is right about the continuum*
- **An independent check.** I wrote the ε-free race from scratch with closed-form Gaussian integrals: no samples, so ⟨φ, ∇φ⟩ = 0 exactly. It runs in 130-digit `decimal` arithmetic and is fully nonlinear: the coefficients are free, there is no merged base and no linearisation. I integrated it from rc_job's state at s₁ = 400. It reproduces RC's least |D| (own run against raw):

| x | decimal race | RC (raw) | difference |
|---|---|---|---|
| 3 | 10^{−98.71610} (share along e +0.981574) | 10^{−98.71606} (+0.981574) | 4.6·10⁻⁵ dex |
| 4.9 | 10^{−27.61785} (+0.805332) | 10^{−27.61786} (+0.805332) | 7.5·10⁻⁶ dex |
| 7 | 10^{−7.74341} (+0.582587) | 10^{−7.74468} (+0.582602) | 1.3·10⁻³ dex |

  At x = 3, from s = 400 to 8000, the two trajectories agree to 1.6·10⁻⁴ dex. Steps h = 0.2 and h = 0.1 agree to 3·10⁻⁹ dex.
- **What RC neglects is negligible in the continuum.**
  - The coefficients' coupling is ∝ F and ⟨φ, ∇φ⟩, both zero at the continuum's rest.
  - Second order in D costs ≤ 1.3·10⁻³ dex, at x = 7, where |D(400)| = 0.045σ. Restarting at s₁ = 1000 gives −7.74341 and −27.61784 (own run).
  - The base's own motion is included.
- **The bookkeeping is sound.**
  - 16 substeps move the least by 1.3·10⁻⁵ dex (raw).
  - The centroid part of the normalised 4-vector is slaved to D (|y|/|D| = 0.04–0.1), so D never loses precision.
  - Including the coefficients naively in the same normalised vector floors at 10^{−38.2}σ, because float64 seeds δα (own run). So the centre-block restriction is necessary for this bookkeeping, and it is justified by the continuum.
- **The injection is placed right.**
  - Injecting at 10⁻⁷ or 10⁻¹¹σ changes s_d by +9.5·10⁻⁷ or −4.2·10⁻⁵ (raw).
  - At x = 7 the injected s_d is 3.27250·10⁵ against the plain run's 3.27266·10⁵ (raw).
  - At x = 6.04 RC's least, 10^{−13.430}, matches the float64 run's 10^{−13.436} (raw). The draft does not cite this agreement.

*s_d and its mechanism*
- **The numbers.** The continuum's s_d values (2.352·10⁷, 1.257·10⁶, 3.58·10⁵, 1.854·10⁵, 8.43·10⁴ at x = 3, 4.5, 4.9, 5.07, 5.48) are the raw's times N₀/N.
- **The shift identity holds.** (Δs_d/Δclimb) − 1 = 2.9·10⁻⁹, 1.3·10⁻⁵, 3.8·10⁻⁵, 4.7·10⁻⁵, 1.3·10⁻⁶ at those x, and −1.1·10⁻⁴ at x = 6.04 (raw). It is unchanged at climb thresholds of 10⁻³ and 10⁻¹σ (own run).
- **C3's capture at x = 5.07 belongs to the window.** s₂ = 1.640·10⁵ lies between Part C's s_d (1.497·10⁵ in R, 1.504·10⁵ in RP) and the continuum's 1.854·10⁵ (raw).
- **The mechanism is right.**
  - On the merged base the coefficients stay at (0.72026204, 0.0797716) from s = 2·10⁴ to 2.5·10⁵ (far window, x = 3; own run).
  - The split is a zero mode: its eigenvalue is −1.4·10⁻¹⁸ … 2.7·10⁻¹², against −4.84, −2.8·10⁻³, −4.8·10⁻⁴ and +9.5·10⁻⁴ for the other modes (own run).
  - The runner-up's coefficient moves only ∝ |D|²: +7.6·10⁻⁹ at |D| = 10⁻³σ and +7.6·10⁻⁷ at 10⁻²σ (x = 4.9; own run).

*Slaving, tolerances, mathematics*
- **κ and the slaving error are as stated.** κ = 0.5062–0.6197. The slaved law's error is 0.40–10.08% at s ≈ 450 and ≤ 0.013% at s ≈ 5000 (raw).
- **The error has a closed form.** It is the next order of eliminating the centroid, ε_b = (1 − κ)·λ_D/λ_b. That gives 8.8% at x = 3 and 3.5% at x = 4.9 on the far window, against 10.1% and 3.6% measured (md).
- **`c4_atol`** = min(10⁻¹⁴, 10⁻⁵ε) is tested exactly. C4 and C4b agree at 30 of 30 points, and R's rule holds at 74 of 74 (md).
- **I re-derived:**
  - the relative-mode identity and 1 − kρ = 1/κ;
  - the merged-block system in (c_b, D) from `critical_jac_nd`'s centre block, whose Schur complement is κA;
  - the pin, with κ cancelling;
  - the direction of the drift and the pin. My forcing model reproduces RP's measured pins: 3.8–3.9·10⁻¹⁵σ against 3.9–4.3·10⁻¹⁵σ at x = 3, and 1.54–1.58·10⁻¹⁴σ against 1.52–1.65·10⁻¹⁴σ at x = 4.9, share +0.88 in both (own run against md).
- **Numbers.** Every number I checked in (iv), (v), the proof, the Status, New and the seven Numerics paragraphs matches the md or the raw. The exceptions are the slips in D3bs-6.

**What is wrong (majors first)**
- **No major issue this round.**
- **D3bs-1 (minor).** "In the continuum" holds only when ⟨·⟩ is the integral. On samples, even in exact arithmetic:
  - **The far window's drift is not rounding.** The draft calls its ∂c*/∂W "the sums' rounding". It is 3.2·10⁻¹⁶σ per unit W, comes from the float64 coordinates of the samples themselves, and is the same in the long double. It decides x ≤ 5.48.
  - **Aliasing decides x = 3.** On an ideal grid the aliasing (1.2·10⁻³⁷σ per unit W) sends x = 3 to the far side. The draft never names aliasing as a decider.
- **D3bs-2 (minor).** "s_d = the climb's end + the hold": the shift is an identity, but the split is a convention.
  - The parts move by ln 10/λ for each decade of the 0.01σ threshold.
  - The hold's decay starts only at the rest: |D| ≈ 0.21σ at x = 3 and 0.53σ at x = 4.9.
  - Before that, the coefficient first rises by about 10⁻³.
- **D3bs-3 (minor).** The proof's y is the slaved centroid's displacement, which is first order in D. It is not a forcing: in the continuum nothing forces D. Setting y to 0 moves the least by 0.016–0.048 dex, not 0.01.
- **D3bs-4 (minor).** Of 62 mutants of the new and changed code, 24 survive and 16 of those change behaviour. The worst three make RC's checks vacuous without any test noticing.
- **D3bs-5 (minor).** Literature and New:
  - missing: Guo et al. 2018 (Gaussian RBF, elimination singularities), Cousseau–Ozeki–Amari 2008, Biehl–Riegler–Wöhler 1996;
  - New still claims the generic trailing;
  - the Du–Hu–Lee citation does not fit Gaussian centres.
- **D3bs-6 (cosmetic).** Slips:
  - "0.01 dex";
  - x = 7 counted as formed;
  - 6.08 against 6.10·10⁴ at "5.48";
  - "8.33σ" (the nearest edge is at 6.25σ);
  - "±0.63";
  - and others listed in the Issues table.

## What I checked and how

**Reproducibility** (`scripts/rerun.py`).
- `--from-raw` with `TN_OUT` on a scratch copy, compared with `cmp` and md5; the test file run with `pytest -p no:cacheprovider`.
- Reruns of rc_job (×4, two of them check rows), rp_job (×2), C4 at (4.9, 10⁻¹¹) and C4M at (7, 10⁻¹³), each compared with the raw field by field.

**The continuum, independently** (`decrace.py`, `rhscheck.py`, `cmp_dec.py`).
- The two-kernel ε-free race uses closed-form inner products, normalised as the far window's mean: ⟨φ(a), φ(b)⟩ = (πσ²/(N ds²))e^{−\|a−b\|²/(4σ²)} with N = 5751.
- Python `decimal` at 130 digits, classical RK4 with h = 0.2 and 0.1, started from rc_job's float64 state at s₁ = 400 and recorded on RC's output grid.
- At that state its right-hand side agrees with PairSums in the long double to ≤ 6·10⁻¹⁸ absolute. The residual is the grid's coordinates, see D3bs-1.

**RC's linearisation** (`rcbase.py`, `rcvar.py`, `spectrum.py`).
- RC's loop reimplemented, with restarts at s₁ = 1000, y = 0 at x = 3, and a naive 6×6 variant.
- The full and centre-block spectra along the merged base.

**The decomposition** (`thresh.py`, `holdstart.py`).
- rc_job and rp_job (Part C's window) with D_CLIMB = 10⁻³, 10⁻² and 10⁻¹σ at x = 3, 4.9 and 5.48.
- The injected race's runner-up coefficient against |D|.

**What moves the rest on samples** (`alias.py`, `gridsrc.py`, `forced.py`).
- Exact (80-digit) sums of ⟨φ, ∇_cφ⟩ and ∂c*/∂W at the merged block on three grids: the committed float64 grid, the same grid stretched but unrounded, and the ideal grid. Done for the far window, Part C's window and the top-1.5 window, and compared with PairSums in float64 and in the long double.
- RC's linear system with the drift added as a forcing of D, κk(∂c*/∂W)dW/ds, validated against RP's measured top-edge pins.

**Mutants** (`mutspec.py`, `mutdriver.py`, `survivor_check.py`, `survivors.sh`).
- 62 single-point mutants of `c4_atol`, `PairSums`, `relative_law`, `growth_rate`, `merged_rhs`, `rc_job`, `rc_lines`, `rp_lines`, `part_RC` and `RC_CHECKS`.
- Each was run through `--from-raw` (md byte diff) and `pytest -x` on the test file.
- Each surviving job-function mutant was rerun on RP and RC rows and compared with the raw.

**Mathematics.** By hand:
- the identity and κ;
- the (c_b, D) system on a merged block, its quasi-static elimination and the next order;
- y's order in D;
- the pin;
- the coefficient couplings ∂α̇/∂c and ∂ċ/∂α, which are ∝ F and ⟨φ, ∇φ⟩.

**Numbers.** I read the md and computed derived quantities from the raw with `numbers.py`.

**Literature.** WebSearch (11 queries).

## Verdicts per part

| claim (D54.md line) | verdict | evidence |
|---|---|---|
| (iv) l. 66–68: a formed runner-up's side is the layer's "in the ε-free inner problem in exact arithmetic (RC; no float64 loop resolves its offset below x ≈ 6), but on samples the window's nearest edge's and, without one, the rounding's" | holds for ⟨·⟩ the integral; fails as a statement about samples (D3bs-1) | decimal continuum race reproduces RC (own run). Counterexamples (own run): on an ideal grid in exact arithmetic, aliasing forces D to 4.7·10⁻⁴⁰σ at x = 3, against the layer's 1.9·10⁻⁹⁹σ, with share −0.99 (far); on the committed far grid the coordinates' drift decides x ≤ 5.48 in any arithmetic |
| (v) l. 119–126: identity, slaved law, κ = W²/(2(w² + w_e²)) | holds; the slaved law is approximate, as now stated | re-derived; the Schur complement of the (c_b, D) system is κA |
| (v) l. 127–129: κ 0.51–0.62; error 0.4–10.1% at s ≈ 450, largest at x = 3 | holds | raw: κ 0.5062–0.6197; error 0.40% (x = 7, top 1) to 10.08% (x = 3, far window); next-order ε_b reproduces it (D3bs-3) |
| (v) l. 130–132: a tail's side is the layer's; closest approach ≥ 1.8·10⁻⁸σ at x ≥ 7 | holds | md (R, C4M) |
| (v) l. 132: "A formed runner-up's relative mode contracts through 15–216 e-folds (x = 7 … 3)" | numbers hold; x = 7 is a tail | raw e-folds 14.7–215.6; md R: w after the layer at x = 7 is 0.0049 < θ\|A\| = 0.008 (D3bs-6) |
| (v) l. 132–136: continuum leasts 10^{−7.7} … 10^{−98.7}σ, the layer's sign | holds, for the integral | raw; decimal race (own run) |
| (v) l. 136–140: the pin, κ cancels, behind the drift; "without an edge, in floating point, the rounding decides" | pin holds; last clause incomplete (D3bs-1) | re-derived; RP (md); exact sums (own run) |
| (v) l. 140–146: s_d = climb's end (0.01σ) + hold; climbs 9.9·10³ … 5.09·10⁴ and 6.97–7.49·10⁴; continuum s_d; "to 10⁻⁴" | numbers and shift identity hold; the split is a convention (D3bs-2) | raw; own run at thresholds 10⁻³ and 10⁻¹σ |
| (v) l. 142–143: merged → same correlation → split frozen; hold starts when the climb ends | mechanism holds; the hold starts at the rest, not at 0.01σ | own run: w frozen, split a zero mode; decay begins at \|D\| = 0.21σ (x = 3) and 0.53σ (x = 4.9) |
| (v) l. 151–155: s_d least at x = 5.48 on both windows; C3's x = 5.07 capture is the window's | holds | md/raw: 6.08·10⁴ (R, x = 5.4772), 8.43·10⁴ (RC); s₂ = 1.640·10⁵ against 1.497·10⁵ and 1.854·10⁵ |
| proof l. 200–208: identity; 1 − kρ = 1/κ; the slaving's error is of the order of the relative rate over the centroid's | holds | re-derived; ε_b = (1 − κ)λ_D/λ_b |
| proof l. 203–204 and 208–211: "D is forced only by the merged block's own offset from its rest, which moves the least \|D\| by 0.01 dex" | imprecise (D3bs-3) | y ∝ D (y/\|D₁\| = 0.105 at x = 3, 0.037 at 4.9; raw); 0.016 dex (raw), 0.048 dex at x = 3 (own run) |
| Status l. 221–228: identity exact; slaved law approximate; RC checked and validated at x = 7 | holds; the checks cover x = 4.9 only, and nothing in the repository checks the deep continuation | raw; my decimal race supplies the missing check (own run) |
| Status l. 243–250 and New l. 257–266 | the sources exist (partly verified); New still claims generic structure (D3bs-5) | search snippets |
| C3 l. 323–329 | holds | md |
| C4 l. 329–340 (30 of 30; departures; closest approach 3.0·10⁻¹⁵ … 1.5·10⁻¹⁴σ; atol) | holds | md, raw; (4.9, 10⁻¹¹) rerun bitwise (own run). At (4.9, 10⁻⁹) and rtol 10⁻¹² the closest approach is 7.45·10⁻¹⁵σ, half the pin (raw; D3bs-6) |
| C4W l. 340–343 | holds | md: 3.95·10⁵ × 1681/1887 = 3.52·10⁵ |
| C4M l. 343–351 | holds | md; (7, 10⁻¹³) rerun bitwise (own run) |
| R l. 351–373 (74 of 74; s_d list; RC cross-references) | holds | md; "6.08·10⁴ (x = 5.48)" is R's x = 5.4772 (D3bs-6) |
| RP l. 374–397 | numbers hold; "∂c_b*/∂W ≈ 3·10⁻¹⁶σ, the sums' rounding" is wrong (D3bs-1); "8.33σ" mis-describes the window (D3bs-6) | md; exact sums (own run) |
| RC l. 397–410 | holds; "the least within 0.01 dex" is 0.016 | raw/md; reruns; decimal race (own run) |

## Status of the round-66 issues

| D3bm | status | evidence |
|---|---|---|
| -1 (major: s_d inherits the climb) | fixed | (v) l. 140–146 now gives s_d = climb's end + hold, with the continuum's s_d 2.352·10⁷ … 8.43·10⁴ and the shift identity to ≤ 4.7·10⁻⁵ at x ≤ 5.48 (raw). C3's x = 5.07 is marked as the window's (l. 154–155, 325–326, 364–365). Residue: the split's convention (D3bs-2). |
| -2 (major: RP's continuation) | fixed | RC (merged base, exact centre block, log norm, injection) replaces it. Its leasts are D3bm-2's and are reproduced by an independent 130-digit nonlinear race (≤ 1.3·10⁻³ dex; own run). "meas" and `centroid_velocity` are removed. Residue: RC's numbers are the continuum's, not the far window's (D3bs-1). |
| -3 (minor: atol artefact) | fixed | `c4_atol` = min(10⁻¹⁴, 10⁻⁵ε). C4 (4.9, 10⁻¹¹) is now captured toward (rerun bitwise; own run). 30 of 30 and 74 of 74 (md). |
| -4 (minor: tests) | partly fixed | The exact c4_atol test kills at1–at3; new tests cover PairSums, merged_rhs, the growth window and a golden rc_job. Still 24 of 62 mutants survive and 16 change behaviour (D3bs-4); round 66 had 17 of 34 and 11. |
| -5 (minor: rate check, κ, slaving) | fixed | Every float64 RP row is printed (md); κ 0.51–0.62; the slaving condition is in (v) l. 127–129, the proof l. 206–208 and the Status l. 223–225. Residue: y and ε_b (D3bs-3). |
| -6 (minor: "in the continuum") | partly fixed | The wording is added (l. 67, 133, 263). But samples in exact arithmetic are still not the continuum: aliasing and the samples' coordinates (D3bs-1). |
| -7 (minor: RP presentation) | fixed | The continuation rows are gone; merged pairs show side None and least "0.00e+00 (—)" (md, far window x = 5.07). |
| -8 (minor: literature, New) | partly fixed | Overlap singularities, Saad–Solla, Şimşek, Zhang and Du–Hu–Lee are added [memory]; New is recast. Residue: D3bs-5. |
| -9 (cosmetic) | fixed | "2.3 steps" (l. 378), "15–216 e-folds" (l. 132), C4's 7.0·10⁻¹⁵ dropped. New slips: D3bs-6. |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bs-1 | minor | **On samples, in exact arithmetic, aliasing and the samples' own float64 coordinates move the block's rest with W, and they decide sides the draft gives to the layer or to "the rounding".** (iv) l. 67–68 names only "the window's nearest edge's and, without one, the rounding's"; (v) l. 139–140 says "without an edge, in floating point, the rounding decides". RP l. 393 calls the far window's ∂c_b*/∂W ≈ 3·10⁻¹⁶σ "the sums' rounding", and the md's RP caption (`rp_lines`, script l. 1801) calls that window "the continuum, where the block's rest does not move with W". **The far window's drift is the coordinates (own run).** Exact (80-digit) sums at the merged block (x = 3, s = 4000) on the committed far grid give ∂c*/∂W = (+3.151·10⁻¹⁶, −8.23·10⁻¹⁷)σ per unit W. The long-double PairSums gives the same, (+3.152·10⁻¹⁶, −8.23·10⁻¹⁷); float64 gives (+2.76·10⁻¹⁶, −7.08·10⁻¹⁷). The source is the samples' coordinates: `np.arange(-2, 2+1e-9, 0.05)` steps by 0.0500000000000000444 and rounds each point, up to 3.7·10⁻¹⁵ off. On that grid S₁/S₀ = −8.40·10⁻¹⁶ per unit length; on the same grid stretched but unrounded it is −3.33·10⁻³⁷. **On an ideal grid the drift is aliasing (own run).** ∂c*/∂W = (1.249·10⁻³⁷, −7·10⁻⁴⁰)σ = (2σ²/W)(2πq/ds), with q = e^{−π²σ²/ds²} = 2.65·10⁻³⁹, (ii)(c)'s own harmonic. **What it does to RC's linear system** (drift added as the forcing κk(∂c*/∂W)dW/ds; own run): on the ideal grid at x = 3 the forced D is 4.7·10⁻⁴⁰σ at RC's least, where the homogeneous D is 1.9·10⁻⁹⁹σ. At s = 8083 it is 1.6·10⁻³⁸σ with share −0.99 along e, against +0.99: the far side on C2's target and toward on the mirror, the continuum's sides reversed. At x = 4.5 the layer still wins: 1.3·10⁻³⁹ against 1.05·10⁻³⁶σ. On the committed grid the forced D is 1.4·10⁻¹⁸ (x = 3), 5.0·10⁻¹⁸ (4.9), 7.4·10⁻¹⁸ (5.48) and 1.1·10⁻¹⁷σ (6.04), share −0.76 … −0.85. That is above the continuum's least for x ≤ 5.48 and below it at 6.04 (3.7·10⁻¹⁴σ). The same method reproduces RP's top-edge pins (summary). **No committed outcome changes:** Part C's top edge dominates there (1.14039·10⁻¹²σ, exact on both grids). | (iv) l. 67–68, (v) l. 132–140 and New l. 263–264: write "in the continuum (⟨·⟩ the integral over ℝⁿ) in exact arithmetic"; and "on samples the nearest edge's, else the aliasing's (≈ 10⁻³⁷σ per unit W here, deciding x ≲ 4), or the samples' float64 coordinates' (3·10⁻¹⁶σ on the far window, in any arithmetic, x ≲ 5.5)". Replace "no float64 loop" with "no loop on these samples, in any arithmetic". RP l. 393 and `rp_lines`' caption: name the coordinates (the long double gives the same). RC (D54 l. 397 and the md caption): "the continuum, with the far window's Jacobian and the rest's drift dropped". Optional: dyadic samples (e.g. ds = 3/64) leave only the aliasing. |
| D3bs-2 | minor | **"s_d = the climb's end + the hold": the shift is an identity, the split a convention.** (v) l. 141–143 says the hold is "pin-independent but starting only when the climb ends (\|D\| = 0.01σ)". **The identity holds and does not depend on the threshold:** 2.9·10⁻⁹ … 4.7·10⁻⁵ at x = 3 … 5.48 (raw); 3.74, 3.75 and 3.76·10⁻⁵ at x = 4.9 for thresholds 10⁻³, 10⁻² and 10⁻¹σ (own run). **The parts move by ln 10/λ per decade** (own run). At x = 4.9 the climb ends at 7.05, 7.33 and 7.60·10⁴ in the continuum and at 3.15, 3.43 and 3.70·10⁴ on Part C's window. "The hold" is therefore 2.876, 2.848 or 2.821·10⁵. At x = 5.48, where there is no hold, "the hold" is 1.44·10⁴, 1.01·10⁴ or 6.0·10³. **The decay does not start at 0.01σ** (own run, injected continuum race). Below the rest the second-moment pull −½M₂P favours a larger M₂ = w(1 + w/w_e)\|δ\|², so the runner-up's coefficient first rises ∝ \|D\|²: +2.3·10⁻⁶ at 0.01σ, peaking at +9.6·10⁻⁴ at \|D\| = 0.21σ (s = 7.20·10⁴; x = 3) and at +1.6·10⁻³ at 0.49σ (s = 7.90·10⁴; x = 4.9). It decays only from the rest onward: down 1% at s = 6.6·10⁵ (x = 3) and 1.01·10⁵ (x = 4.9). | (v) l. 141–146 and RC l. 403–408: "s_d − s_climb(θ) is the same on every window for any θ in the linear range (10⁻³–10⁻¹σ checked): the departure is a time-shift of one trajectory, set by the offset the climb starts from; the coefficient's decay begins at the rest (\|D\| ≈ 0.21σ, 0.53σ), after a rise ∝ \|D\|²". |
| D3bs-3 | minor | **The proof's y is the slaved centroid's first-order response, not a forcing.** Proof l. 203–211 says "F(c_b) = 0 at D = 0" and "there D is forced only by the merged block's own offset from its rest, which moves the least \|D\| by 0.01 dex". **Derivation.** From `critical_jac_nd`'s centre block on a merged block, with μ = Λ₀Gm, Q = ⟨r, ∇²φ⟩ and H = ⟨∇φ∇φᵀ⟩: Ḋ = μ[(2ww_e/W)QD + (w − w_e)(Q − WH)y] and ẏ = μ[((w² + w_e²)/W)(Q − WH)y + (ww_e(w − w_e)/W²)QD]. The slaved y = −(Q − WH)⁻¹(ww_e(w − w_e)/(W(w² + w_e²)))QD is first order in D. **Numbers.** In the layer, \|y\|/\|D\| ≈ ((A − W)/A)\|ww_e(w − w_e)\|/(W(w² + w_e²)) = 0.092 at x = 3, s₁ = 400; measured 0.105 (raw y and lg_D1), 0.037 at x = 4.9 (raw) and 0.025 at s₁ = 1000 (own run). So in the continuum nothing forces D: y is part of the homogeneous solution. Setting y to 0 moves the least by 0.016 dex at x = 4.9 (raw: −27.634 against −27.618) and by 0.048 dex at x = 3 (own run: −98.764). **The κ-law's error.** The κ-law drops ẏ; at next order λ ≈ κA/(1 + ε_b), with ε_b = (1 − κ)λ_D/λ_b = (w − w_e)²ww_e Q/((w² + w_e²)²(Q − WH)). That gives 8.8% (x = 3) and 3.5% (x = 4.9) on the far window at s ≈ 450, against 10.1% and 3.6% measured (md). | Proof l. 203–211: "the slaved centroid sits at c_b*(W) + y(D), y first order in D; the κ-law drops ẏ, an error ε_b = (1 − κ)λ_D/λ_b"; "in the continuum nothing forces D; y is part of the homogeneous solution, and y = 0 moves the least by 0.016 (x = 4.9) to 0.05 dex (x = 3)". l. 405: "within 0.02 dex". |
| D3bs-4 | minor | **Tests: 24 of 62 single-point mutants survive both the test file and `--from-raw`, and 16 of them change behaviour** (own run, survivors rerun against the raw). **Killed:** all of `c4_atol` (3) and `merged_rhs` (3); 11 of 12 `PairSums`; rc_job's main path (10, by the golden x = 5.48 row); md-changing `rc_lines`/`rp_lines` (5). **relative_law:** rl1 (ρ's sign), rl2 (ρ′ = 0), rl3 (block at the entrant) and rl6 (κ = W²/(2w_e²)) move RP's law by up to 6%, 3%, 4% (dcdW 145%) and 106% on the far window at x = 3; the tests check the law only at x = 4.9 and 6.04, within 3–10%. **growth_rate:** gr4 (no D > 0 filter) turns `grow` into NaN on far-window rows x = 3 and 5.07. **rc_job:** rc6 (no "twice the least") moves x = 7's injection from s = 2.398·10⁴ to 4.64·10³; rc7 (\|share\|) flips the mirror's e_least to +0.806; rc8 (e-folds in ln 2) gives 17.5 for 58.1; rc10 and rc11 (the far-side stop) run the mirror to the end (status 0 for 1). **rc12, rc13, rc14 (ymode, sub, inj ignored):** these make the y = 0, 16-substep and injection-level checks copies of the main row (y = 0: least −27.6179 for −27.6343), and `test_the_continuum_rc_on_the_raw`'s check assertions (0.02 dex, 2·10⁻⁴) would still pass on such a rerun. **part_RC / RC_CHECKS:** pr1 (no checks), pr2 (no mirror) and pr3 (16-substep check at 4) — no test touches them. **Equivalent on the data:** rl10, mo3 (3·10⁻⁹), gr1, gr3, gr5 (≤ 2·10⁻¹⁴), rc17, rcl4, rpl3. | `tests/test_transit_nd.py`: golden `rc_job` rows for a mirrored gap (side, e_least < 0, status, s_climb), for x = 7 (s_inj) and for each check row, asserting that each check differs from the main row (e.g. y = 0 by > 0.01 dex); `efolds` in the golden comparison; relative_law's law at far-window x = 3 to 1%; rp_job's `grow` against the raw; a `part_RC` job-list test (14 main rows plus 4 checks with their arguments). |
| D3bs-5 | minor | **Literature and New.** (a) **Verified by search snippets:** Amari–Park–Ozeki, Neural Comput. 18(5):1007–1065 (2006); Fukumizu–Amari, Neural Networks 13:317–327 (2000); Şimşek et al., ICML 2021 (PMLR 139:9722–9732); Zhang et al. NeurIPS 2021; Du–Hu–Lee NeurIPS 2018 (homogeneous models). Wei et al. 2008 and Wei–Amari 2008 are confirmed only through citing papers (Guo et al. 2018), with venue and pages unconfirmed. I did not find Saad–Solla, PRE 52 (1995). (b) **Missing:** Guo, Wei, Ong, Rubio Hervas, Zhao, Wang, Zhang, JMLR 19 (2018), on Gaussian RBF networks, where plateaus are "mainly caused by elimination singularities" — the w → 0 line where (v)'s tails wait "only because its mobility is ∝ w" (l. 95); Cousseau–Ozeki–Amari, IEEE TNN 19(8):1313–1328 (2008); Biehl–Riegler–Wöhler, J. Phys. A 29:4769 (1996), for the plateau length. Barber–Saad–Sollich (EPL 1996) show that more symmetric starts are trapped longer (snippet). (c) **New** l. 261–263 still claims "the relative mode's law in this setting … the lighter kernel trails the block … forced by the centroid's motion". All three ingredients are generic: the trailing, Ḋ ⊃ (w − w_e)·(common force), is immediate in any (a, θ) gradient flow; κ = 1/(1 + z²), with z = (w_e − w)/W, is the Schur complement of its mobility, and Wei et al. 2008 reduce the overlap singularity to an (h, z) system [snippet: Orhan–Pitkow], probably with this factor (unverified); and the rest moves with W for any input measure that is not translation-invariant, of which the window's edge is one case. (d) Du–Hu–Lee's balancedness is for positively homogeneous units, and Gaussian centres are not. (e) The frozen split is the neutral direction along Fukumizu–Amari's line of critical points of a split unit. | Status l. 243–250: add Guo et al. 2018, Cousseau–Ozeki–Amari 2008 and Biehl–Riegler–Wöhler 1996; cite Fukumizu–Amari at l. 142–143; restrict Du–Hu–Lee to "(a, θ) gradient descent". New l. 261–263: claim the race's consequences in this paradigm (the side, the pin, the climb and s_d), not the law; check Wei et al. 2008's reduced equations for the z-factor. |
| D3bs-6 | cosmetic | **Slips.** (1) l. 405 "the least within 0.01 dex" and l. 211 "0.01 dex": the y = 0 check moves it 0.0164 dex (raw). (2) l. 132 ("A formed runner-up's … 15–216 e-folds (x = 7 … 3)"), l. 133–135 and New l. 264 ("10^{−7.7} …") include x = 7, which is a tail; the formed range is 27–216 e-folds and 10^{−13.4} … 10^{−98.7}σ. (3) l. 144 says "6.97–7.49·10⁴" but l. 403 says "6.97·10⁴ … 7.50·10⁴"; the latter includes the mirror's 7.499·10⁴ (raw). (4) l. 355's "6.08·10⁴ (x = 5.48)" is R's x = 5.4772, while l. 408 quotes RP's x = 5.48 as 6.10·10⁴ (md). (5) Status l. 224 says "within 0.4–10.1% … at every RP row"; that range is at s ≈ 450, and at s ≈ 5000 the error is ≤ 0.013% (raw). (6) l. 388 "at 1.5 (8.33σ…) … the rounding level (1.95·10⁻¹⁷σ)": the nearest edge is the right one at 6.25σ, whose exact drift is 4.7·10⁻¹⁹σ along x (own run; long double 5.4·10⁻¹⁹, md), and float64's larger component is x's 6.65·10⁻¹⁷ (md). (7) l. 396 "±0.63": the mirror's share is −0.64 (md). (8) Status l. 226: the checks are at x = 4.9 only; cite x = 6.04's agreement (RC 10^{−13.430} against float64 10^{−13.436}; raw). (9) `c4_atol`'s docstring says "at ε = 10⁻¹¹ and 10⁻¹⁴"; drop 10⁻¹⁴. (10) C4 at (4.9, 10⁻⁹) and rtol 10⁻¹² has its closest approach at 7.45·10⁻¹⁵σ, half the pin (C4b 1.50·10⁻¹⁴; raw), so there "the window's pin" is resolved only to a factor 2. | as listed |

## Sources

**Repository** (commit 37918e1, exported): `review/drafts/D54.md`, `experiments/exp_transit_nd.py`, `results/chaos/transit_nd.md`, `results/chaos/transit_nd_raw.npz`, `tests/test_transit_nd.py`, `src/core/rich_transit_nd.py`, `src/core/rich_transit.py`, `review/REVIEW_round66_D3bm.md`, `review/REVIEW_LOG.md` (Round 66), and `git diff 6d48964 37918e1`.

**Web.** WebSearch on 2026-10-28, search result snippets only; I read no full texts.
- Amari, Park, Ozeki 2006: https://research.knu.ac.kr/en/publications/singularities-affect-dynamics-of-learning-in-neuromanifolds/
- Fukumizu, Amari 2000: https://www.ism.ac.jp/~fukumizu/localmin3.pdf and https://link.springer.com/chapter/10.1007/978-3-7908-2656-2_3
- Guo et al. 2018, with its account of Wei–Amari 2008 and Wei et al. 2008: https://mlanthology.org/jmlr/2018/guo2018jmlr-numerical (the JMLR pages were listed but could not be fetched)
- Wei et al.'s (h, z) system as described by Orhan–Pitkow: https://arxiv.org/pdf/1701.09175 (snippet only)
- Cousseau, Ozeki, Amari 2008: https://pure.teikyo.jp/en/publications/dynamics-of-learning-in-multilayer-perceptrons-near-singularities/
- Şimşek et al. 2021: https://proceedings.mlr.press/v139/simsek21a.html
- Zhang et al. 2021: https://neurips.cc/virtual/2021/poster/28572
- Du, Hu, Lee 2018: https://proceedings.neurips.cc/paper/2018/hash/fe131d7f5a6b38b23cc967316c13dae2-Abstract.html
- Biehl, Riegler, Wöhler 1996: https://mailman.srv.cs.cmu.edu/pipermail/connectionists/1996-March/016988.html
- Barber, Saad, Sollich 1996: https://publications.aston.ac.uk/id/eprint/679/1/0295-5075_34_2_151.pdf
- Plateau length against the log of the initial variance: https://arxiv.org/pdf/2104.14546 (snippet)
- Wei–Amari IJCNN 2006 precursor: https://dblp1.uni-trier.de/rec/conf/ijcnn/WeiA06.html

**Failed.** WebFetch could not resolve jmlr2020.csail.mit.edu or www.jmlr.org in this session's network settings, so I stopped fetching.

**From memory, not verified this round:**
- the volumes and pages of Wei–Zhang–Cousseau–Ozeki–Amari (Neural Comput. 20, 2008) and Wei–Amari (Neural Networks 21, 2008);
- Saad–Solla, PRE 52:4225 (1995), and whether it states the log(1/asymmetry) escape;
- the content of Wei et al.'s z-dependence.

## Lesson

A claim made in the continuum but computed on samples needs an inventory of what the samples add — the window's edges, aliasing, and the samples' own float64 coordinates — measured at the depth the claim reaches. RC is right about the continuum: an independent 130-digit nonlinear race in closed-form integrals reproduces its least to 5·10⁻⁵ dex at x = 3. It is right precisely because it drops the drift of the block's rest. At 10⁻⁹⁹σ, though, that drift is large on every sampled problem. On an ideal grid, aliasing alone (10⁻³⁷σ per unit W) forces D 59 decades above the continuum's least and onto the other side. `np.arange`'s coordinates (3·10⁻¹⁶σ per unit W, the same in the long double) do the same for every gap up to 5.48. The second lesson is the same idea on a smaller scale: test a decomposition at a threshold by moving the threshold. Here the shift survived the move and the parts did not — the climb's end slid by ln 10/λ per decade, and the hold turned out to begin at the rest, after a rise of the coefficient.
