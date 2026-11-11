# Referee report D3bz, round 79: Proposition D‴.54 as revised after round 72

## Scope

**What I reviewed.** Commit 7d2e76d on `claude/festive-bardeen-eopej2` ("D‴.54 round 72: what moves the rest on samples; the shift an identity; the unslaved rate"). I exported it with `git archive 7d2e76d | tar -x` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref79/D3bz_work/repo`. Below, W stands for `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref79/D3bz_work`.
- **`review/drafts/D54.md`.** I read all of it, concentrating on:
  - (iv) l. 61–90 and (v) l. 93–173;
  - the proof, l. 175–231;
  - Status, Classical and New, l. 233–293;
  - the Numerics paragraphs C4, R, RP (l. 402–427), RC (l. 427–444) and RD (l. 444–448).
- **`experiments/exp_transit_nd.py`.** `sample_coords`, `exact_drift`, `alias_drift`, `PairSums`, `relative_law`, `rp_job`, `growth_rate`, `merged_rhs`, `merged_block`, `rc_job`, `rd_job`, `part_RC`, `part_RD`, `rp_lines`, `rc_lines`, `rd_lines`, `main`, `c4_atol`.
- **Results and tests.** `results/chaos/transit_nd.md`, `transit_nd_raw.npz`, and `tests/test_transit_nd.py` (55 tests).
- **History.** `review/REVIEW_round72_D3bs.md`; the Round 72 entry of `review/REVIEW_LOG.md`; `git diff c22b8f1 7d2e76d`; and only the Issues tables of D3bm and D3bi.

**Settings.**
- The repository's venv: Python 3.13.16, numpy 2.5.4, scipy 1.18.1.
- Every run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`.
- The experiment always ran with `TN_OUT=W/out`, which holds a copy of the raw.
- Mutants were applied only in copies under `W/mut/ws1`–`ws3`.
- At most two heavy processes of mine ran at a time. The machine was shared, with load up to about 6.
- Scripts are in `W/scripts`, logs in `W/logs`, and decimal-race outputs in `W/out/dec` and `W/out/vec`.

**Repository untouched.** In `/home/user/Lyapunov-Convergence-Kernels-Showcase` I ran only read-only git commands (`archive`, `log`, `show --stat`, `diff`, `status`, `rev-parse`).
- `git status` is clean.
- HEAD has since moved to 5201021 (D‴.52), but `git diff 7d2e76d HEAD` is empty for every D54 file.
- Nothing of mine is still running (checked with `ps`).

**Literature.** 7 WebSearch queries, read as snippets only. No WebFetch. arxiv.org and ams.org were not accessed.

## Summary

### What checks out

**Reproducibility** [own]
- `--from-raw` regenerates the md byte for byte (md5 6da44019f2ab8e3803eacad6a3984493).
- All 55 tests pass, in 62.9 s.
- These jobs rerun bit-identically to the raw in every field except `secs`: `rd_job` (both targets), and `rc_job` at (near tie, 3), (near tie, 6.04), (mirror, 3) and the y = 0 check.

**RC's forced response is right, confirmed independently** [own]. I wrote a fully nonlinear two-kernel race from scratch (`decrace.py`, not `PairSums`) in Python `decimal`:
- exact sums over the far window's samples, with the committed np.arange coordinates converted exactly (40 digits) or on the ideal grid (60 digits);
- classical RK4 with h = 0.2 (h = 0.1 changes \|D\| by ≤ 2.8·10⁻⁸);
- started from the float64 race's state at s = 400;
- its right-hand side matches the long-double PairSums to ≤ 6·10⁻¹¹ relative.

At every RC output from s = 495 to 6000, its D equals RC's linear D_h + D_f in direction (0.00°) and in size to ≤ 7·10⁻⁴:

| run | decimal race | RC prediction | outcome |
|---|---|---|---|
| committed grid, x = 3 | 10^{−17.848} (share −0.853) at s = 3999.7 | 10^{−17.848} (−0.853) | −0.984 at s = 6000: far |
| mirror, committed, x = 3 | +0.486 at s = 4000 | same | +0.970 at 6000: toward feature 2 |
| committed, x = 5.48 | 10^{−17.139} (−0.771) | D_h 10^{−18.998} + D_f 10^{−17.133} | coordinates decide |
| ideal grid, x = 3 | 10^{−39.324} (−0.724) | 10^{−39.325} (−0.724) | aliasing decides |
| ideal grid, x = 4.5 | 10^{−35.978} (+0.861) | the layer's; D_f(alias) is 10^{−38.878} | layer decides |

So (iv) and (v)'s "the coordinates decide x ≤ 5.48, the aliasing x = 3 and not 4.5, far on C2's target and toward on the mirror" holds in exact arithmetic on these samples. Round 72's derivation of the forcing also holds:
- The merged base's own equation gives F = Wċ_m/(μ(w² + w_e²)) exactly, so Ḋ = κk·ċ_m.
- Each kernel's share, ∓w·(w − w_e)ċ*/(w² + w_e²), is independent of W. That is why my mutant rc19 is equivalent.
- The augmented exponential is the exact variation of constants over a substep.

**RD** [own]
- My own 80-digit sums at x = 3, s = 4000 give the raw to 16 digits:
  - committed grid: (3.151289·10⁻¹⁶, −8.227459·10⁻¹⁷)σ;
  - ideal grid: (1.248587·10⁻³⁷, −6.958·10⁻⁴⁰)σ.
- I re-derived (2σ²/A)(2πq/ds)sin(·): the Poisson sum, the sign, and ∂F/∂c = −A‖φ‖²/(2σ²)I at the feature. It gives 1.2486·10⁻³⁷σ.
- On the dyadic grid, float64's and the long double's values are in the ratio 2048.6, the ratio of their unit roundoffs (2048).

**The (D, y) system and ε_b.**
- Re-derived from ċ_i = Gm w_i⟨r∇_cφ_i⟩. The determinant term is (ac − bd) = μ²ww_e. The slaved law is κλ_DQ exactly, whatever R is. ε_b = (1 − κ)λ_D/λ_b = (w − w_e)²ww_e/(w² + w_e²)²·q/r.
- The first-order law κλ_D(1 − ε_b) is within 0.42% of the Jacobian where the κ-law errs by 10.1% (far window, x = 3, s ≈ 450). The corresponding pairs are 0.21%/7.1% (top 1.5), 0.17%/6.3% (right edge) and 0.12%/5.3% (Part C's) [own, from the raw].

**The shift identity, the climbs and the coefficient** [own, from the raw]
- s_d − s_climb(θ) agrees between Part C's window and every other window of C2's target (top 0.95, top 1.5, the long double, the far window) to ≤ 1.2·10⁻⁹. Against RC it agrees to ≤ 2.5·10⁻⁵ for formed gaps and 2·10⁻⁴ at x = 7. This is better than the md shows.
- The hold moves by about ln 10/λ per decade of θ.
- The coefficient rises ∝ \|D\|² (ratios 99.6 and 106.7 per decade of \|D\|); its peaks and 1% decays are as stated.

**Literature additions** (all [snippet])
- Guo et al., JMLR 19:1–39 (2018): "the plateau phenomena are mainly caused by the elimination singularities".
- Cousseau–Ozeki–Amari, IEEE TNN 19(8):1313–1328 (2008).
- Biehl–Riegler–Wöhler, J. Phys. A 29(16):4769 (1996), on plateau length.
- Yan–Wang–Rigollet, Ann. Statist. 52(4):1774–1795 (2024).
- Du–Hu–Lee is now correctly restricted to positively homogeneous units [memory].

**Numbers.** Every number I checked in (iv), (v), the proof, the Status, New and Numerics RP/RC/RD matches the md or the raw, apart from the items in D3bz-3 and D3bz-7.

### What is wrong (majors first)
- **No major issue.**
- **D3bz-1 (minor). Test gaps.** No test reruns `rd_job` or `part_RD`, or the (D, y)/ε_b code in `rp_job`. And `_same_round72`'s default atol (10⁻⁸) makes its d_coord, d_alias and w_climbs[0] comparisons vacuous. Of my 60 mutants, 31 are killed, 9 survive and are equivalent, and 20 survive while changing an output.
- **D3bz-2 (minor). "Deciding x ≤ 5.48, in any arithmetic" overreaches.** The forced D (10^{−17.85} … 10^{−17.13}σ) lies below a float64 state's quantum (1.85·10⁻¹⁶σ). No run in the repository can follow it, and the far-window float64 runs take the predicted side at only 2 of 9 gaps. RC's "float64 sums" column predicts the float64 race wrongly at 7 of 9.
- **D3bz-3 (minor). The coordinates' mechanism is misdescribed.** The quoted 3.7·10⁻¹⁵ is a harmless uniform stretch. The drift comes from ±2.2·10⁻¹⁶ roundings at five points, and only on the far window's axis from −2. On the top-1.5 window the exact drift is the right edge's (4.71·10⁻¹⁹σ), not "the rounding level".
- **D3bz-4 (minor). The (D, y) "agreement within 5·10⁻⁴" is a change of coordinates.** The first-order residual (≤ 0.42%) is the real check.
- **D3bz-5, -6, -7 (cosmetic).**
  - The criterion: the least against the growing mode; and the 10^{72}σ linear size behind "+81.4 dex".
  - "Mobility ∝ coefficient" against (iii)'s constant Λ₀Gm.
  - Slips.

## What I checked and how

**Reproducibility** (`rerun.py`, `cmputil.py`).
- `--from-raw` regenerated the md, compared with `cmp` and md5.
- I ran pytest and reran the jobs listed above, comparing each with the raw field by field.

**Fully nonlinear exact race** (`rcvec.py`, `decrace.py`, `check_decrhs.py`, `cmpdec.py`).
- `rcvec.py` copies rc_job's continuation but keeps the vectors D_h(s), D_f,coord(s), D_f,alias(s) and the race's state at s₁.
- `decrace.py` integrates the race in decimal from that state, on the committed grid (40 digits) and the ideal grid (60 digits), using exact Gaussians by recurrence. Its right-hand side was checked against PairSums.
- Runs: x = 3 (both grids, plus the mirror on the committed grid), x = 5.48 (committed), x = 4.5 (ideal), and an h = 0.1 repeat. Each takes 4–8 min.

**RD independently** (`mydrift.py`).
- My own 80-digit ∂c*/∂W on: the committed grid; the ideal grid; the stretched-only grid; and the grid with rounding kept only within 3σ of the block.
- The same on Part C's, the top-1.5, the right-edge and the far window.
- np.arange's coordinate errors split into the stretch and per-point rounding.

**The (D, y) system and ε_b.** By hand, then checked numerically against every RP Jacobian row (`holds.py`).

**Holds and the coefficient** (`rawnums.py`, `holds.py`). s_d − s_climb(θ) on every RP and RC row, compared across windows. The coefficient peaks compared with the leading-order rest.

**Criterion** (`rcvec.py` outputs). The forced/homogeneous ratio and the shares, at the least and at s = 2·10⁴.

**Mutation testing** (`mutspec.py`, `mutdriver.py`, `survivor_check.py`, `classify.sh`).
- 60 single-point mutants of round 72's code. Each went through `--from-raw` (md byte comparison) and `pytest -x` on the test file.
- Survivors were rerun against the raw (rd_job ×2; rp_job far window at x = 3 and 4.9; rc_job ×4 or ×14).

**Literature.** WebSearch, 7 queries.

## Verdicts per part

| claim (D54.md lines) | verdict | evidence |
|---|---|---|
| (iv) l. 65–70: a formed runner-up's side is the layer's in the continuum; on samples the edge's, else the aliasing's (deciding x = 3, not 4.5) or the coordinates' ("in any arithmetic, deciding x ≤ 5.48") | holds in exact sums; "in any arithmetic" fails for float64 (D3bz-2) | decimal races (own); RP far-window float64 sides (raw): 2 of 9 on the predicted side |
| (v) l. 121–131: identity, κ-law, κ 0.51–0.62, slaving error 0.4–10.1% | holds | re-derived; raw: κ 0.506–0.620 |
| (v) l. 132–139: tails' side; 27–216 e-folds (x = 7: 15); continuum offsets 10^{−13.4} … 10^{−98.7}σ | holds | raw: 27.03 … 215.58 e-folds and 14.74; my ideal-grid race at x = 4.5 gives 10^{−35.978} against RC's −35.977 |
| (v) l. 139–142: the merged block is no solution; the forcing κk(∂c_b*/∂W)dW/ds; the pin; κ cancels; behind the drift | holds | exact for the merged base with ċ_m; the decimal race shows the pin plateau at 10^{−17.96}σ for s = 1.1–2.7·10³ |
| (v) l. 142–143: D away from the nearest edge, whatever the arithmetic | holds | RP: −1.55 (long double) and −1.57·10⁻¹⁴σ (float64) |
| (v) l. 143–149: the aliasing formula, 1.25·10⁻³⁷σ; the coordinates' 3.15·10⁻¹⁶σ "off that grid by up to 3.7·10⁻¹⁵"; forced D 10^{−17.9} … 10^{−17.1}σ; the sides; float64 rounding ≤ 1.9·10⁻¹⁶σ; quantum 1.85·10⁻¹⁶σ | numbers and sides hold; the mechanism is misattributed (D3bz-3) | own sums; np.arange decomposition (own) |
| (v) l. 149–156: s_d − s_climb(θ) the same on every window; climbs 9.9·10³ … 5.09·10⁴ and 6.97–7.49·10⁴; continuum s_d; ln 10/λ (2.876, 2.848, 2.821·10⁵) | holds, and is understated | all six windows of C2's target agree to ≤ 1.2·10⁻⁹ (own, from the raw) |
| (v) l. 156–160: frozen split; rise ∝ \|D\|² (+1.4·10⁻⁸, +1.4·10⁻⁶, +1.5·10⁻⁴); peak at the rest (0.21σ, 0.49σ); 1% decay at 6.56·10⁵ and 1.01·10⁵ | holds; "at the rest" plausible but unchecked | raw; the peaks are 1.1–1.4× the leading-order rest (0.188 and 0.348σ) (own) |
| proof l. 215–220: identity; 1 − kρ = 1/κ | holds | re-derived |
| proof l. 221–225: the (D, y) system, y, the κ-law, ε_b, "within 5·10⁻⁴" | the formulas hold; the 5·10⁻⁴ is an identity (D3bz-4) | re-derived; raw (4.3·10⁻¹⁶ for merged pairs) |
| proof l. 225–230: the drift; ∂F/∂c = −A‖φ‖²/(2σ²)I; y = 0 moves the least 0.016 dex; the aliasing formula | holds | re-derived; raw 0.0164 dex |
| Status l. 233–254 | holds, except that the forced response is presented without any nonlinear check (D3bz-2), and D3bz-4 | — |
| Literature l. 254–282 | the sources exist; wording (D3bz-6) | snippets |
| New l. 282–293 | appropriately narrowed; holds | — |
| RP l. 402–427 | numbers hold; l. 416–417 mislabels the long double's 5.4·10⁻¹⁹σ (D3bz-3) | md; own sums on the top-1.5 grid |
| RC l. 427–444 | holds; the forced columns are confirmed nonlinearly (own); the md's "float64 sums" side claim fails (D3bz-2) | decimal races; RP far rows |
| RD l. 444–448 | holds | own 80-digit sums |
| tests (55) | pass, but rd_job and the (D, y) code are untested (D3bz-1) | mutants |

## Status of the round-72 issues

| D3bs | status | what remains |
|---|---|---|
| -1 (what moves the rest on samples) | fixed | RD numbers reproduced exactly (own); RC's forced response confirmed by an exact nonlinear race (own). Residue: "in any arithmetic, deciding", the f64 column (D3bz-2); the mechanism and the top-1.5 edge (D3bz-3) |
| -2 (s_d's split) | fixed | the identity holds on all C2 windows (≤ 1.2·10⁻⁹), but the md prints only two; w_peak at x = 6.04 and 7 is the departure, not a peak (D3bz-7) |
| -3 (the proof's y; the κ-law's error) | fixed | the formulas are right; the "5·10⁻⁴" is an identity; quote the first-order residual (≤ 0.42%) instead (D3bz-4) |
| -4 (surviving mutants) | partly fixed | 31 of my 60 mutants are killed, including every rc_job forcing mutant. rd_job/part_RD and the (D, y)/ε_b code are never rerun; the default atol makes three comparisons vacuous; 20 output-changing survivors (D3bz-1) |
| -5 (literature, New) | fixed | additions verified by snippets. Residue: the "mobility" wording; Saad–Solla PRE 52 still not found; Wei et al.'s factor still [unverified] (D3bz-6) |
| -6 (slips) | fixed (all ten) | c4_atol's docstring is now redundant; l. 416–417 still calls the long double's value "the rounding level" (D3bz-3, -7) |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3bz-1 | minor | **Tests.** No test reruns `rd_job` or `part_RD`, or rp_job's (D, y)/ε_b code; the tests read only their committed output. `_same_round72` compares d_coord and d_alias (values of order 10⁻¹⁶ and 10⁻³⁷) and w_climbs[0] (7·10⁻⁹) with numpy's default atol = 10⁻⁸, so those checks pass whatever the values. **Results** (60 single-point mutants; own run): 31 killed. 9 survive and are equivalent on the data: sc3; ed2 (that term is O(drift²)); ed7; mb2; rc10; rc19 (W cancels in κk/W); rd3; rl1; rdl1. 20 survive and change an output: rd1 (dy_ld 4·10³× off); rd2 ('ideal' 2.5·10²¹× off); rd4 (dy_alias 1.2·10⁸× off); rd5 ('ld' 38–103% off); rd6, rd7, rd8; ad5 (dyadic aliasing 549× off — 2lo/ds is an integer only for ds = 0.05); rp1 (ε_b 2×, law_y 1.85×); rp2 (ε_b 44×); rp3 (law_y 2.6% at s ≈ 450 and 700× at 5000); rp4 (ε_b 7.9%); rp6 (law_y 7.5×); rp7 (≤ 5·10⁻⁶); rc12 (argmax\|rel\| gives w_peak −80% at \|D\| = 2σ and removes s_wdecay at x = 3–5.07 — no golden row has a hold); rc8 and rc16 (d_change); rc20, rc21, rc22 (d_coord ×2, d_alias's sign flipped, w_climbs[0] +72%, all passing) | `tests/test_transit_nd.py`: a golden `rd_job(('near tie', 3.0))` compared with atol = 0 (it runs in 1 s); law_y and eps_b added to the far-window rp_job rerun (rtol 10⁻⁹); a golden rc_job row with a hold (x = 4.9: w_peak, s_wdecay); atol = 0 in `_same_round72` for d_coord, d_alias, w_climbs and w_peak; d_change in the golden comparison |
| D3bz-2 | minor | **What decides x ≤ 5.48, and in which arithmetic.** (iv) l. 68–70 says the coordinates' 3·10⁻¹⁶σ decides x ≤ 5.48 "in any arithmetic". RC's md caption says each forced column's sign decides "wherever it is the larger". **In exact sums the prediction is right**: my nonlinear race reproduces D_h + D_f to 7·10⁻⁴ (see the Summary). **No run in the repository resolves it.** The forced D at the least (10^{−17.85} … 10^{−17.13}σ) is 30–120× below a float64 state's quantum (1.85·10⁻¹⁶σ). That quantum also applies to RP's long-double runs, whose state is float64. **The far window's float64 runs:** C2's target goes far, toward, toward, merged, far at x = 3, 4.5, 4.9, 5.07, 5.48; the mirror goes far at all five. The prediction is far on C2's target and toward on the mirror, so they agree at 2 of 9 decided gaps. RP l. 423–425 itself calls these sides "the integrator's". **The "float64 sums" forced column** predicts far on C2's target and toward on the mirror (+81.5 … +1.8 dex at the injection), so it is wrong for the float64 race at 7 of 9. Its forcing is noise in the float64 force at the merged base's ulp-quantized centre, not a drift | (iv) l. 68–70 and (v) l. 145–148: "the same in exact sums and in the long double; deciding x ≤ 5.48 in a loop resolving 10⁻¹⁸σ (an exact race reproduces RC's forced D to 7·10⁻⁴, referee round 79), not in float64, whose state quantises D at 1.85·10⁻¹⁶σ". RC md: relabel the "float64 sums" column as a loop with float64 sums and an exact state, which no run has, or drop it. Status l. 248–249: say that the forced response is linear and checked nonlinearly. Optionally add the decimal race to RC (about 5 min per gap) |
| D3bz-3 | minor | **What "the samples' coordinates" are.** (v) l. 146 and RD's caption: "off that grid by up to 3.7·10⁻¹⁵". **That offset is a harmless stretch.** It is np.arange's uniform stretch k·δ, with δ = fl(lo + ds) − lo = 0.05 + 4.44·10⁻¹⁷ (3.55·10⁻¹⁵ at k = 80). A stretched grid only aliases: my exact sums on the stretched-only far grid give 1.2486·10⁻³⁷σ, the ideal grid's value [own]. **The drift is the rounding of k·δ ≥ 2 at ulp 4.4·10⁻¹⁶.** On the far window's x-axis (from −2), the points x = 0.05, 0.15, 0.25, 0.35, 0.45 sit ∓2.22·10⁻¹⁶ off the stretched grid; on its y-axis, y = 0.55 and 0.65. The rounding within 3σ of the block alone gives (3.04, −0.98)·10⁻¹⁶σ of the full (3.15, −0.82)·10⁻¹⁶σ. **The axes from −1 add nothing.** On Part C's, the top-1.5 and the right-edge windows no point within 3σ is off the stretched grid; exact committed equals ideal there to 10⁻⁸. **So on the top-1.5 window** the exact drift, (4.71·10⁻¹⁹, −5·10⁻²²)σ on both grids, is the right edge's (6.25σ). RP l. 416–417 still calls the long double's 5.4·10⁻¹⁹σ "the rounding level", the residue of D3bs-6(6). Scaling RC's response linearly [own estimate], that edge would force 10^{−20.1}σ at x = 4.9 and 10^{−20.7}σ at x = 3, far above the continuum's offset, so in exact sums both would go far. **Docstring:** `sample_coords` says np.arange rounds "lo + k·fl(ds)" | (v) l. 146 and the RD caption: "np.arange's uniform stretch (3.6·10⁻¹⁵ at the ends, which only aliases) plus its rounding of k·δ at five points within 3σ (±2.2·10⁻¹⁶), specific to the far window's axis from −2". RP l. 416–417: "float64's rounding (6.65·10⁻¹⁷σ); in exact sums the right edge's 4.71·10⁻¹⁹σ (the long double's 5.4·10⁻¹⁹σ)". `sample_coords` docstring: δ = fl(lo + ds) − lo |
| D3bz-4 | minor | **The (D, y) "agreement within 5·10⁻⁴" (proof l. 225, Status l. 244, REVIEW_LOG lesson) is a coordinate change.** For a merged pair, the (D, y) matrix is the Jacobian's centre block written in (D, c_b). So the agreement is 4–6·10⁻¹⁶ wherever \|D\| ≤ 10⁻¹⁴σ: e.g. 4.3·10⁻¹⁶ at x = 3 on Part C's window, s = 448.6 [md]. The 4.7·10⁻⁴ is the far window at x = 7, where \|D\| = 0.035σ at s ≈ 450; at 6.04 (\|D\| = 0.014σ) it is 1.0·10⁻⁴, so it scales ∝ \|D\|². The informative check is the first-order law: κλ_D(1 − ε_b) is within 0.42%, 0.21%, 0.17% and 0.12% of the Jacobian where the κ-law errs by 10.1%, 7.1%, 6.3% and 5.3% [own, from the raw] | proof l. 225, Status l. 244, rp_lines: "κλ_D(1 − ε_b) within 0.42% of the Jacobian where the κ-law errs by 10.1% (the (D, y) matrix is the Jacobian's centre block for a merged pair)" |
| D3bz-5 | cosmetic | **The criterion and the injection ratios.** What decides the side is each response's amplitude on the growing mode. RC compares norms at the homogeneous least, where the two shares differ (+0.98 against −0.85 at x = 3). This changes no verdict here: from the least to s = 2·10⁴ (both shares ±0.99) the ratio moves by ≤ 0.22 dex, and by ≤ 0.6 dex up to the injection. The closest margins are 1.6 dex (mirror, x = 5.48, coordinates) and 2.9 dex (x = 4.5, aliasing) [own]. "+81.4 dex" is a ratio of two linear responses; the linear forced D there is 10^{72.4}σ [raw] | RC md caption: "the ratio of the two linear responses on the growing mode (the forced one's linear size is not physical)" |
| D3bz-6 | cosmetic | **Literature wording.** l. 260 and 273 use "mobility ∝ weight/coefficient" in the WFR sense: the factor on ⟨r∇_cφ⟩. But (iii) l. 41–42 gives the centres the constant mobility Λ₀Gm; it is the force −∂E/∂c_i = w_i⟨r∇_cφ_i⟩ that carries w. Saad–Solla, PRE 52:4225 (1995): search finds only its NeurIPS 1995 version [snippet]; the PRE reference is [memory]. Wei et al. 2008's factor remains unverified; the snippets describe only universal trajectories | l. 260 and 273: "a position's speed ∝ its coefficient (the force carries it; the mobility Λ₀Gm is constant)"; keep the [unverified] tag |
| D3bz-7 | cosmetic | **Slips.** (1) RC md: at x = 6.04 and 7 the coefficient's "peak" is the departure (\|D\| = 2.008σ and 2.004σ, s = s_d), not a peak. (2) `c4_atol` docstring: "…at ε = 10⁻¹¹ … ; referee round 72: at ε = 10⁻¹¹)", which is redundant. (3) RD md: "the drift is the samples' coordinates', in any arithmetic", yet float64 differs from it by up to 2.06·10⁻¹⁶σ, 65% of it (mirror, s ≥ 4000: (1.09, −0.14) against (3.15, −0.82)·10⁻¹⁶σ) [raw]. (4) Status l. 249: "each arithmetic's rounding — computed exactly". (5) "The same on every window": the md shows two windows, though all six of C2's agree to ≤ 1.2·10⁻⁹. The mirrored target differs from C2's by 0.25% at x = 4.9, as a different target should. (6) "Peaks as D reaches its rest": the peaks (0.208, 0.397, 0.491σ) are 1.1–1.4× the leading-order rest (0.188, 0.304, 0.348σ); this is consistent with V's leading order being low near the fold, but it is unchecked | as listed; print RP's other windows in the identity paragraph |

## Sources

**Repository** (commit 7d2e76d, exported): `review/drafts/D54.md`, `experiments/exp_transit_nd.py`, `results/chaos/transit_nd.md`, `results/chaos/transit_nd_raw.npz`, `tests/test_transit_nd.py`, `src/core/rich_transit.py` (critical_w, curve_mobility, transit_scale), `review/REVIEW_round72_D3bs.md`, `review/REVIEW_LOG.md` (Round 72), the Issues tables of D3bm and D3bi, and `git diff c22b8f1 7d2e76d`.

**My scripts and outputs** (under `W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref79/D3bz_work`):
- Scripts in `W/scripts`: rerun.py, cmputil.py, rcvec.py, decrace.py, check_decrhs.py, cmpdec.py, rawnums.py, holds.py, mydrift.py, mutspec.py, mutdriver.py, survivor_check.py, classify.sh.
- Logs in `W/logs`: the mutation results in mut1.jsonl, mut2.jsonl and mut3.jsonl, and the survivor classification in classify.log.

**Web** (WebSearch 2026-10-29, snippets only):
- Guo et al. 2018: https://www.jmlr.org/papers/v19/16-210.html
- Cousseau–Ozeki–Amari 2008: https://pure.teikyo.jp/en/publications/dynamics-of-learning-in-multilayer-perceptrons-near-singularities/
- Biehl–Riegler–Wöhler 1996: https://mailman.srv.cs.cmu.edu/pipermail/connectionists/1996-March/016988.html; plateau-length usage: https://www.esann.org/sites/default/files/proceedings/legacy/es2019-173.pdf
- Yan–Wang–Rigollet 2024: https://projecteuclid.org/journals/annals-of-statistics/volume-52/issue-4/Learning-Gaussian-mixtures-using-the-WassersteinFisherRao-gradient-flow/10.1214/24-AOS2416.full
- Saad–Solla (NeurIPS 1995 version): https://papers.neurips.cc/paper_files/paper/1995/file/a1519de5b5d44b31a01de013b9b51a80-Paper.pdf
- Barber–Saad–Sollich 1996: https://publications.aston.ac.uk/id/eprint/679/1/0295-5075_34_2_151.pdf
- Wei et al. 2008 and Wei–Amari 2008, as described in citing works (search-result snippets; not fetched): https://arxiv.org/pdf/1701.09175 and https://jmlr2020.csail.mit.edu/papers/volume19/16-210/16-210.pdf

**From memory, not verified:** Saad–Solla, PRE 52:4225; Du–Hu–Lee NeurIPS 2018 (positively homogeneous models); Fukumizu–Amari, Neural Networks 13:317 (2000), the line of critical points.

## Lesson

A linear response at 10⁻¹⁸σ can only be tested by a computation that resolves 10⁻¹⁸σ. RC's forced response was right, but nothing in the repository could show it: its float64 runs quantise D at 1.85·10⁻¹⁶σ and take the integrator's side. A from-scratch race in exact sums confirmed it to 7·10⁻⁴ in about five minutes per gap. The size of an error also says little about its effect: np.arange's 3.7·10⁻¹⁵ offset is a harmless stretch, and the drift comes from a ±2.2·10⁻¹⁶ rounding at five points. And a test proves nothing about a job if it only reads the job's committed output, or compares 10⁻¹⁶-sized numbers under numpy's default atol of 10⁻⁸.
