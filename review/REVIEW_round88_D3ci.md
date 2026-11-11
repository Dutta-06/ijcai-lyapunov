# Referee report D3ci, round 88: Proposition D‴.58 as revised after round 84

## Scope

- **Subject.** Commit 4bc7b07, branch claude/vigilant-clarke-hf939b. I exported it with `git archive` into W/repo, where W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref88`.
  - W/wmout is WM_OUT and holds a copy of the raw.
  - My scripts are in W/scripts, logs in W/logs, mutant slots and results in W/mut, probes in W/probe.
- **Read in full.** D58.md, `exp_width_mismatch.py` (2139 lines), `width_mismatch.md`, `test_width_mismatch.py` (53 tests), REVIEW_round84_D3ce.md, the Round 84 entry of REVIEW_LOG.md (:3348–3388), D61.md, both round-84 mutant specs and `tools/mutants.py`.
- **Read in part.** `git diff 902fa4a 4bc7b07` for these files (the hunks, plus a word diff of D58).
- **Run settings.**
  - Python 3.13.18, numpy 2.5.4, scipy 1.18.1; long double is x87 (ε 1.08·10⁻¹⁹).
  - Every process ran with OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1 and `nice -n 15`, one heavy process at a time.
  - pytest ran with `--basetemp` under W, also inside the mutant runner (my copy, `W/scripts/mutants88.py`).
- **Nothing in the repository was modified.** I ran only `git archive/log/diff/status`, and `git status` is clean.
  - Meanwhile another session moved HEAD to a02bbba (D‴.60, D‴.61, HANDOFF). None of those commits touches a D‴.58 file.
  - All my processes have exited. One background *waiting* loop (an `until pgrep` loop that matched its own command line) was stopped at its time limit; it did no computation and lost nothing.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my computations; [snippet] web-search results only (arxiv.org, ams.org, ar5iv and eprints.ncl.ac.uk unreachable); [memory] recalled facts.

## Summary

**Design.**
- I reproduced the md, the tests and nine sets of jobs.
- I recomputed A6's exact errors with my own code, and B5's anchored ends and the rounded best end with my own exact evaluator.
- I refitted all 36 A6 limits with seven ansätze, computed each A6 rule's Mehler-frame ratio ρ, and checked the closed-form identities.
- I ran 46 new single-point mutants of this round's code, plus a control, and re-ran 15 of the round's own and the round-84 referee's (13 survivors and the two controls).

**Verdicts.**
- **Reproducibility holds** [own]:
  - `--from-raw` regenerates the md byte for byte.
  - 53/53 tests pass in 51.9 s.
  - Every sampled job is bit-identical to the raw in every field.
  - My independent A6 code agrees with the raw to the last bit of the ratio, and to 25 digits of the errors, at 6 cells up to K = 160.
- **Every number** of this round's new D58 text traces to the md and the raw.
- **The mathematics of (ii)'s formula is right**: its identities, the cusp slope (1 − t²)/2, and the Laplace sketch's algebra (consistent with D‴.61).
- **The cusp claim is strongly supported.** The 1/K coefficient of ratio·t − 1 is −0.250000 at τ/σ ≤ 1 and −0.250003 at 2, so e_K ∝ K^{1/4}t^K.
- **The vσ*²/τ² branch is robust.** Its excess decays cleanly as K⁻², and the limits sit within 10⁻⁷ of t(v) by any matching ansatz.
- **Round 84's issues are fixed.** I re-ran the round's survivor accounting and confirmed it.
- **I found no major issue.**

**Findings.**
1. **D3ci-1 (minor): the τ²/v branch's limits depend on the ansatz.** The quoted "within 6.1·10⁻⁵" and "+5.7·10⁻⁹ … +6.1·10⁻⁵ above" come from forcing L + c/K² + d/K³.
   - Other reasonable ansätze give −1.4·10⁻⁴ … +1.2·10⁻⁴ at τ/σ = 2, v/v_hs = ½.
   - Across the branch's 13 rules the range is −1.4·10⁻⁴ … +1.2·10⁻⁴; the uniform positive sign is an artefact of the K⁻² ansatz.
   - The conjecture holds there to about 1.4·10⁻⁴, not 6.1·10⁻⁵.
2. **D3ci-2 (minor): "the double parameters' own floor" is not a floor.** The best end rounded to double (1.07·10⁻¹² at τ/σ = ¼, K = 6) is undercut by B5's own anchored double end at κ = ¼ (1.04·10⁻¹³, 0.0972 times it). At τ/σ = 1 it is undercut at every κ ≤ 1 (0.192–0.668 times it). So "at the lowest optimum the double parameters cost as much" (D58 :83, Status :185) is contradicted by the round's own data, and the κ = 0 stalls are left unexplained.
3. **D3ci-3 (minor): D58 never refers to D‴.61.** Its Laplace sketch is still "the referee's" and "unproved", without:
   - D‴.61's conditional proof;
   - the sketch's domain: only rules with a Mehler frame, |ρ| < 1, which is 22 of A6's 36 rules;
   - the fact that A6E's real-axis ratio is not bounded below uniformly on |ξ| ≤ 2√K. Its conjectured limit (1 − s²/4)^{1/2} vanishes at s = 2, and the minimum falls with K (0.804, 0.741, 0.677).
4. **D3ci-4 (minor): 13 of 46 new mutants survive both checks; 3 change published results on a rerun.**
   - P01 (`a6_error`'s nodes at `gh_nodes`' default 110 digits): A6 at τ/σ = ¼, K = 60 reads 64.6435 instead of 65.0619.
   - P02 (`gen_ratio`'s nodes likewise): A6E at K = 80, s = 0.1 on the real axis reads 2.06·10⁶⁴ instead of 0.99879.
   - R01: three cusp labels flip.
   - The real-job tests run only where 110 digits suffice.
5. **D3ci-5 (cosmetic): seven wording and printing items**, listed in the Issues table.

## What I checked and how

**Reproduction.**
- `--from-raw` into W/wmout, then `cmp` with the committed md: identical (md5 5f40713…; raw 959dc055…).
- The test file: 53 passed in 51.93 s.
  - This includes `test_b4_job_on_three_cells_against_the_raw` (B4 at τ/σ = ¼, K = 5 from both starts, and τ/σ = ¾, K = 11 from the optimum: every field) and `test_b5_job_against_the_raw`.
- `b5_rerun.py`, compared with `==` in every field:
  - B5 anchored at (¾, 11): 6 κ, 16 s;
  - B5 direct at (1, 12, θ = ¼): 8 κ, 35 s. Neither job is covered by the tests.
- `a6_rerun.py`:
  - A6 at (2, ½·v_hs) to K = 160 (6 s), (1, Mercer-zero) to K = 160 (18 s), (½, Fock) to K ≤ 80 (3 s) and (¼, Fock) to K = 160 at 642 digits (127 s): header and rows bit-identical;
  - A6E at K = 40 on the real axis (all 20 s values) and K = 80 on the imaginary axis (s = 0.1, 1, 2): bit-identical;
  - `part_A5(series=False, far=False)` rows and checks, `a5_series` at all four τ/σ, and `a5_far_row` at (2, 120), (1, 60), (½, 100), (¼, 60): bit-identical.

**Independent recomputation** [own].
- **A6 (`my_a6.py`).**
  - Nodes by Newton on the orthonormal Hermite recurrence from numpy's nodes. Weights by the Christoffel sum 1/Σ_{j<K}p_j(z)², which is independent of `gh_nodes`' formula.
  - Error from my own derivation: ‖f*‖² = √πσ*, ⟨f*, φ(·−y)⟩ = σ*√π(2/(2 + τ²))^{1/2}e^{−y²/(2(2+τ²))}, Gram √πe^{−(y−y′)²/4}.
  - At (2, ½·v_hs) K = 80 and 160, (2, v_hs) 160, (1, ½·v_hs) 160, (¼, Mercer frame) 80 and (½, 2·v_hs) 120, the ratio agrees with the raw to relative 0.0 and e₁ to 25 digits.
- **B5 (`b5_check.py`).** From the raw's rms strings I recomputed every excess over B4's best end, the per-start slopes (1.44–2.09 on [1/16, 1], median 1.64), the geometric means (1.88, 1.57, 1.57), the spreads (maximum 8.15·10³ at κ = 0) and the anchored slopes (2.03, 2.15, 2.11). All identical to the md.
- **Exact evaluator (`dblfloor.py`).** Every double taken exactly, exp at 90 digits. It reproduces the stored best-end rms to ≤ 6.6·10⁻⁵², and every B5 anchored end and `rms_d` to the printed digits (D3ci-2).
- **`kappa_blend`.** At κ = 1 it equals r_d bit for bit at 6 ends, and no sample lies outside Sterbenz's range there.
- **Identities**, to double rounding:
  - tanh((r + ln cosh r)/2) = λ/(2 − λ) = t/(1 − t + t²);
  - tanh((r + |r − ln β²|)/2) = q, and tanh r = λ;
  - the cusp's slope is (1 − t²)/2, and λ/t = 1 + sech r;
  - the Mercer frame's rule's rate b(b − 1)/(b² + b + 2), b = β², is attained at vσ*²/τ² (algebra).
- **Frames.** ρ = 1 + 2σ²/τ² − 2σ²/v for all 36 A6 rules: |ρ| < 1 at 22 (4, 4, 6, 8 at τ/σ = ¼, ½, 1, 2), ρ = 1 at moment matching, and |ρ| > 1 at the other 10 (D3ci-3).

**Limits (`a6_refit.py`).** For every A6 rule I fitted seven ansätze:
- the published L + c/K² + d/K³ through K = 80, 120, 160;
- L + b/K + c/K² through the same three K;
- L + b/K + c/K² + d/K³ through four K;
- L + aK^{−3/2} + bK^{−2} + cK^{−5/2} through four K;
- L + cK^{−p} with p free, through three K;
- two-point K⁻² and K^{−3/2}.

I also computed the excess's local exponents between successive K, and the cusp's 1/K coefficient by 3- and 4-point fits.

**Mutation testing** (`W/mut/spec88.py`).
- 46 mutants of this round's new code that are in neither round-84 list, plus a control. They cover the precision bookkeeping, `a6_lines`, `rate_conj`, `a6_scalings`, `a6_job`, `gen_ratio`, `dcos`, `a6e_lines`, B5's lines, jobs and job, `kappa_blend`, `lm_fit`'s κ and anchor log, B4's new summaries and floors, and A5's far line.
- Each was run by a private copy of `tools/mutants.py` (`pytest -x` plus `--from-raw` and a byte comparison).
- I also re-ran the round's own survivors BF4, B53, L54, DD1, BE1, BL2 and FS2 (with X00), and the round-84 referee's BJ05, LF01, BB02, FA06, SL02 and SL03 (with CTRL).
- Behaviour-changing survivors were probed by rerunning the mutated `a6_job` and `gen_ratio` against the raw (`probe.py`).

**Literature.** Four web searches, all reported as [snippet] (Sources).

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i), (iv), (v), proofs (:8–14, :116–165) | Unchanged this round; not re-derived | word diff |
| (ii) one curve t(v) = tanh((r + \|ln(v/v_hs)\|)/2); λ at u = r, λ/(2 − λ) at u = −ln cosh r, q at u = r − ln β²; the cusp t + ((1 − t²)/2)\|u\| (:23–26) | Correct as stated (conjectural) | identities [own] |
| (iii) the direct floor "superlinear": 1.1·10³ … 2.2·10⁸ (a = 0.91–2.52, median 1.59); B5 κ^{1.44–2.09}, median 1.64; geometric means 1.88, 1.57, 1.57 (:87–91) | Correct. The 8.15·10³ spread is at κ = 0, outside the "κ = 1/16 … 1" just named (7.15·10³ within) | [md], b5_check [own]; D3ci-5b |
| (iii) the anchored floor κ^{2.03–2.15} over ½ … 4; flat stops 2 of 3 at ½, all at ≥ 1; at τ/σ = 1 0.277, 0.192, 0.637, 0.668 (:76–78) | Numbers correct. "Quadratic" is a least-squares slope over per-doubling exponents of 0.07–3.35 | [md], [own]; D3ci-5e |
| (iii) last-anchor floor 1.4·10⁻¹⁶ … 1.5·10⁻¹², ×0.50–2.53 of the gap; end ×0.13–4.43; 2σ² term 0.10–0.51 (:78–82) | Correct | [md] |
| (iii) "the double parameters cost as much": 1.1·10⁻¹², κ = 0 at 1.0·10⁻¹² (0.954), 3.1·10⁻¹⁴ and 1.8·10⁻¹³ (5.77) (:83–85) | Numbers correct; conclusion contradicted by B5's own κ = ¼ row | D3ci-2 |
| (iii) best end; Gauss–Newton gaps 2.1·10⁻¹⁹, 55 of 72 at 1.000–1.495 (:62–66) | Correct | [md] |
| A5/A6 paragraph: 32 rules within 6.1·10⁻⁵; branches −2.1·10⁻⁹ … +8.7·10⁻⁸ (K⁻², 2.00–2.04) and +5.7·10⁻⁹ … +6.1·10⁻⁵ (K^{−1.52 … −1.98}); "holds there to the precision the run reaches" (:221–230) | Numbers trace to the md. The τ²/v-branch limits are ansatz-dependent to 1.4·10⁻⁴, and the "above" sign is an artefact | D3ci-1 |
| Cusp: K⁻¹, within 3.1·10⁻⁷; e_K ∝ K^{1/4}t^K; remainder ratios 0.247 … 0.569 (:230–238) | Correct and strongly supported (1/K coefficient −0.250000 … −0.250003) | [md], [own] |
| Laplace sketch; A6E 0.677–1.021 and 1.001–1.413 at K = 80, "bounded … converging in K" (:231–236) | Numbers correct; "bounded" holds only above on the real axis; domain unstated; D‴.61 not cited | D3ci-3 |
| Status (ii), (iii) (:176–187) | As above | D3ci-1, -2, -3 |
| Numerics A5, A6, A6E, B5 (:268–283, :305–308) | Trace to the md. ":281 at vσ*²/τ² within 8.7·10⁻⁸" omits "away from the cusp" | D3ci-5a |
| Literature (:190–255) | Karvonen–Oates–Girolami's non-matching bounds and Karvonen–Särkkä's condition match my snippets. No one-formula rate found | Sources |

## Status of the round-84 issues

| D3ce | status | note |
|---|---|---|
| 1 | fixed | "Linear" is withdrawn. B5 measures the order: direct κ^{1.44–2.09}, anchored κ^{2.03–2.15}; "superlinear … path-dependent" in (iii) and the Status. I reproduced it exactly from the raw, and re-ran two B5 jobs bit-identically. |
| 2 | fixed for its list | The round's accounting holds [own]: BF4, B53, L54 and BJ05 are now killed by the tests; DD1, BE1, BL2, FS2, LF01, BB02, SL02, SL03 and FA06 survive as recorded (RS01 not re-run). New survivors in this round's code: D3ci-4. |
| 3 | fixed | A6 and A6E state the formula as a conjecture; τ/σ = ¼'s turnaround and its limit are in D58. Residues: D3ci-1, D3ci-3. |
| 4 | fixed, with a residue | The last anchor's floor, the 2σ² term and the rounded best end are in place. Residue: the rounded best end is read as a floor, which B5 contradicts (D3ci-2). |
| 5 | fixed | Exact excesses; 9.5·10⁻¹⁷. |
| 6 | fixed | `frame_series` stops only on a nonzero small term; a K = 1 test exists; the md's 1.5·10⁻¹¹ is gone (1.4·10⁻¹⁴). |
| 7 | fixed | (a)–(g) all done. |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ci-1 | minor | **The τ²/v branch's limits are ansatz-dependent.** Where: Status :177–178; A5/A6 paragraph :224–230; Numerics :280–281; md :171, :173.<br>**The published fit and the alternatives [own, from the raw].** The published fit forces the excess to fall as K⁻². Other fits of the same ratios:<br>• τ/σ = 2, v/v_hs = ½: published +6.1·10⁻⁵; L + b/K + c/K² −1.4·10⁻⁴; four-point with K⁻¹ −8.9·10⁻⁵; free exponent −5.7·10⁻⁵ (p = 1.37); half-integer powers −2.9·10⁻⁵; two-point K⁻² +1.2·10⁻⁴.<br>• Over the branch's 13 rules: −1.4·10⁻⁴ … +1.2·10⁻⁴ (τ/σ = 1, ½·v_hs: −4.4·10⁻⁵ … +2.2·10⁻⁵; τ/σ = 2, Fock: −6.1·10⁻⁵ … +3.6·10⁻⁵).<br>**No power law yet.** The local exponents still drift: 1.25, 1.34, 1.40, 1.46, 1.52 at τ/σ = 2, ½·v_hs. No ansatz is justified, and "all 13 above" is the K⁻² ansatz's bias.<br>**The robust parts.** The vσ*²/τ² branch: exponents 2.00–2.02, within 7.5·10⁻⁶ by every ansatz and ≲ 10⁻⁷ by those matching K⁻². The cusp: b = −0.25000.<br>**Same rule, two limits.** τ/σ = ¼, Mercer frame gets 65.054019 (A5, K = 100, 120, 160) and 65.052913 (A6, K = 80, 120, 160), both quoted. Its three-point ansätze give −3.7·10⁻⁵ … −4.1·10⁻⁵, while the K = 160 ratio is +1.1·10⁻⁴ above and falling. | Quantify "to the precision the run reaches" as ≈ 1.4·10⁻⁴ (τ/σ = 2, ½·v_hs) and ≲ 5·10⁻⁵ elsewhere on that branch, and drop "above it". Better: have `a6_lines` print a second ansatz with a K⁻¹ term or a free exponent beside the published one, plus the bracket "the K = 160 ratio +4.0·10⁻⁴ above, still falling". |
| D3ci-2 | minor | **The rounded best end is a scale, not a floor.** Where: D58 :83–85, Status :185; md :311 label, B5's last column.<br>**Exactly evaluated [own].** At (¼, 6), B5's anchored double end at κ = ¼ (stall, 953 accepted steps) lies 1.04·10⁻¹³ above the best end. That is 0.0972 times the rounded best end's 1.07·10⁻¹², and 12× below κ = 1's 1.28·10⁻¹².<br>**At (1, 12)** the κ ≤ 1 ends are 0.192–0.668 times the rounded best end.<br>**So the double parameters' floor at the lowest optimum is ≤ 1.0·10⁻¹³.** The κ = 1 end fits the frozen-rounding law there (κ = ½ → 1: ×4.3). The κ = 0 stalls are not explained by the parameters: 1.02·10⁻¹² after 64 accepted steps at (¼, 6), and 1.76·10⁻¹³ at (¾, 11), above κ = 1's 9.45·10⁻¹⁴. | Reword: "the best end rounded to double lies 1.1·10⁻¹² above it — a scale, not the parameters' floor: B5's anchored double end at κ = ¼ lies 1.0·10⁻¹³ above". In the Status, replace "as large as the double parameters' own" with what κ = 0 and ¼ show. State that the κ < ½ stalls are unexplained. Relabel md :311. |
| D3ci-3 | minor | **D58 is out of step with D‴.61.**<br>**No cross-reference.** D58 never names D‴.61 (I grepped: D‴.52, .53, .58, .59 only), though the round-84 log says the mechanism is taken up there.<br>**The sketch's status.** D58 keeps the sketch as "the referee's … Unproved". D‴.61 (iv) proves the asymptotics conditionally on conjecture P, controls the tails unconditionally, and gives (i)'s sharp constant conditionally.<br>**The domain.** The sketch needs a frame (\|ρ\| < 1, v ∈ (τ²σ²/σ*², τ²)): 22 of A6's 36 rules [own]. Moment matching (ρ = 1) and 13 others (\|ρ\| > 1) rest on A6's numbers alone.<br>**A6E on the real axis.** "Stays bounded … converging in K" holds only above there: the minimum falls 0.804, 0.741, 0.677 at K = 20, 40, 80 (s = 1.6–1.7). The ratio to D‴.61's limit (1 − s²/4)^{1/2} at s = 1.9 is 2.98, 2.72, 2.44, and the limit vanishes at s = 2 [raw, own]. The mass sits at s = 2√t(v) (1.59 at τ/σ = 2, ½·v_hs), where convergence is slowest — D‴.61's account of the slow branch. | One sentence in the A5/A6 paragraph and the Status pointing to D‴.61 (iii)–(v): the conditional proof, its domain v ∈ (τ²σ²/σ*², τ²) (22 of 36 rules), and the rules outside it as numerical only. Restate A6E as "bounded above on \|ξ\| ≤ 2√K, converging pointwise for s < 2 (to D‴.61's (1 ∓ s²/4)^{1/2})" in D58 :233–235 and md :179. |
| D3ci-4 | minor | **13 of 46 new mutants survive both checks.** Tests kill 26 and the md alone 7 (L03, L05, L06, L08, F04, Q01, A02); the control passes both.<br>**Survivors that change results on a rerun (probed):**<br>• P01, `a6_error` with `gh_nodes(K)` (110 digits): A6 at τ/σ = ¼, Fock, K = 60 gives 64.64355 instead of 65.06189 (e₁ 3.575·10⁻¹⁰⁹ instead of 3.552·10⁻¹⁰⁹). It is harmless where the cancellation is below 110 digits, e.g. (½, ½·v_hs) at K = 120.<br>• P02, `gen_ratio` likewise: K = 80, s = 0.1 on the real axis gives 2.06·10⁶⁴ instead of 0.99879.<br>• R01, `rate_conj` with `a > b`: flips "τ²/v, the cusp" at τ/σ = ¼, ½, 2, where a == b exactly in double.<br>**Cause.** The real-job tests run `a6_job` only at (2, Mercer frame), K ≤ 40 (77–94 digits), and `gen_ratio` at K = 20 (120 digits), where 110-digit nodes suffice.<br>**Equivalent on the data:** L02, J03, G02 (probed: identical), F03, F07, Q02, Q07, Q08, A01, A03. | Add `a6_job` at τ/σ = ¼ with A6_K = (60,) (~10 s), or `a6_error` at K = 59, 60 with 279 digits, against the raw's e₀/e₁. Add `gen_ratio(80, 0.1, 'real', 300)` against the raw (~1 s). Label the cusp "(both arguments)". |
| D3ci-5 | cosmetic | **(a)** Numerics :281 "at vσ*²/τ² within 8.7·10⁻⁸": add "away from the cusp" (near it that branch reaches −4.0·10⁻⁵).<br>**(b)** :91 "up to 8.15·10³ apart" is at κ = 0, outside "κ = 1/16 … 1" (7.15·10³ within).<br>**(c)** md B5: the θ = 0 start at (¼, 6) accepts 2 steps at κ = 2 and 4, so its excess is the start's 3.1·10⁻⁵ both times. Its [¼, 4] slope 2.64, the per-start maximum, is censored.<br>**(d)** md :140, :149, :158, :167: the cusp's "attained at" is a floating-point tie (at τ/σ = 1, b exceeds a by 1 ulp).<br>**(e)** :76 "κ^{2.03–2.15}" summarises per-doubling exponents of 0.07–3.35; at (1, 12) the excess is flat from κ = ½ to 1 (1.06·10⁻¹⁶, 1.12·10⁻¹⁶).<br>**(f)** md :177, :282 "to 2.2·10⁻¹⁶" is one double ulp: e₀ and e₁ are taken as floats, so the comparison's resolution, not an agreement measured below it.<br>**(g)** Seven print lines are pinned by the md comparison only, not by the tests: L03, L05, L06, L08, F04, Q01, A02. | Reword as indicated; (c) mark censored rows (accepted steps ≤ 2); (f) compare the Decimal strings. |

## Sources

- [Karvonen–Särkkä, Gaussian kernel quadrature at scaled Gauss–Hermite nodes, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532); [Aalto record](https://aaltodoc.aalto.fi/items/0c4e072c-7344-4a4f-ac45-7a6e2cca9121) [snippet]
  - Scaled nodes with approximate weights; exponential rate under a weight-growth assumption.
  - States that Kuo–Woźniakowski's Gauss–Hermite rule converges only for ε² < ½.
- [Karvonen–Oates–Girolami, Integration in RKHS of Gaussian kernels, arXiv:2004.12654](https://arxiv.org/pdf/2004.12654); [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/) [snippet]
  - Scaled Gauss–Hermite rules converge exponentially for every length-scale.
  - Bounds exp(−c₁N^{1/d})N^{1/(4d)} ≲ · ≲ exp(−c₂N^{1/d})N^{−1/(4d)}, c₁ > c₂.
  - Notes that Kuo–Woźniakowski 2012's lower bound is erroneous.
  - No explicit dependence on the scale in the snippets.
- [Scaling Optimized Hermite Approximation Methods, arXiv:2412.08044](https://arxiv.org/pdf/2412.08044) [snippet]
  - The arXiv number and title match D58's Hu–Yu citation; authors and DOI not seen.
- [Wang–Zhang, root-exponential convergence of Gauss–Hermite quadrature, arXiv:2312.07940](https://arxiv.org/abs/2312.07940) [snippet]
- [Multi- and infinite-variate integration … Gaussian kernels, arXiv:2412.05368](https://arxiv.org/pdf/2412.05368) [snippet]
- [How sharp are error bounds? Lower bounds on quadrature worst-case errors for analytic functions](https://research.manchester.ac.uk/en/publications/how-sharp-are-error-bounds-lower-bounds-on-quadrature-worst-case-/) [snippet]
- **One formula for the rates.** No snippet gives a closed form for a scaled Gauss–Hermite rule's rate in a Gaussian-kernel space; the known results are non-matching two-sided bounds. As far as these sources show, t(v) = tanh((r + |ln(v/v_hs)|)/2) is new.
- **Not reachable.** eprints.ncl.ac.uk (KOG's repository record) is outside this session's network; arxiv.org, ams.org and ar5iv were not fetched.
- **[memory].** The Christoffel-sum weights; the Gaussian inner products; Sterbenz's lemma.

## Lesson

An extrapolated limit is only as precise as its ansatz. Where the convergence exponent is still drifting, the fit's form sets both the size and the sign of the deviation, so the spread over two or three reasonable ansätze should be printed beside the limit.

Likewise, a "floor" must be a lower bound. One rounded point is a scale: B5's own κ-scan already held a double point ten times closer.

And a precision argument the tests exercise only at small K is not pinned. The two mutants that corrupt A6 and A6E by orders of magnitude survive because every real-job test runs where the default 110 digits suffice.
