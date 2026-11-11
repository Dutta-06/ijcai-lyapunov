# Referee report D3bh, round 61: Proposition D‴.59 as revised after round 57

**Scope.** I reviewed commit a1f5d85 from a `git archive` export. The files were:
- `review/drafts/D59.md`;
- `review/REVIEW_round57_D3bd.md` and the round-57 entry of `review/REVIEW_LOG.md` (l. 2155–2191);
- `experiments/exp_moving_width.py`;
- `results/chaos/moving_width.md` and `moving_width_raw.npz`;
- `src/core/moving_width.py`, `src/core/rich_transit_nd.py` and `src/core/rich_transit.py`;
- `tests/test_moving_width.py`.

Other sessions moved HEAD to ada0e1c during the review. None of these files changed, so everything below also holds at HEAD.

I modified nothing in the repository. All runs were in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref61/D3bh_work/`:
- scripts in `scripts/`, logs in `logs/`;
- loop states in `out/` and `out2/`;
- mutants in `mut/`.

Every run used the venv python with `OMP_NUM_THREADS=1`, `nice -n 10` and `PYTHONDONTWRITEBYTECODE=1`. I ran at most two heavy processes at a time and signalled only processes I had started.

## Summary

**What checks out (verified).**
- **Reproducibility.**
  - `--from-raw` into my own `MW_OUT` is byte-identical to the committed md. The script reads the raw from `MW_OUT`, so I copied the raw there first.
  - The 9 tests pass in 13 s.
  - Reruns with the current code are bit-identical to the raw:
    - A at σ*/σ₀ = 3, both paths;
    - D at d = 0.6;
    - H "σ* 0.04, equal";
    - E target 2 at s = 10⁻³.
  - All 36 rows of G's loop table reproduce their midpoint offsets to three digits. Target 18 at ε = 10⁻⁸ matches to |Δd| ≤ 5·10⁻¹⁰ and |Δrms| ≤ 1.5·10⁻¹³.
- **Numbers.** Every number in the numerics and status paragraphs matches the md or raw, except "to 2·10⁻⁴ on target 12", which is 2.5·10⁻⁴.
- **(i) in the continuum.**
  - J, ∂_cJ, ∂_σJ, the single critical point and the limits re-derive.
  - The convergence argument holds: E ≥ J, the sublevel set is compact in ℝ × (0, ∞), |w| is bounded away from 0 and ∞, the metric is non-degenerate there, and Łojasiewicz applies. There is one slip.
  - The restriction to the continuum and to a start of the feature's sign is stated correctly.
- **(ii)(a).**
  - The expansion and the identity ∂³_cφ = σ⁻¹∂_c∂_σφ − σ⁻²∂_cφ are right.
  - The width split σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ) takes the one-Gaussian residual from O(Δc³) to O(Δc⁴). Per halving of Δc the residual falls by 8.0 with equal widths and by 16.0 with the split. A cancelling (−0.786, 0.357) pair gives 16.2.
  - The split matches C's targets 1 and 7 and D's d = 0.6 to 5·10⁻⁵, 7·10⁻⁵ and 6·10⁻⁴.
  - It also holds in the real loop's landed pairs, to 4·10⁻⁴ on target 29 at ε = 10⁻¹⁰. The draft does not use this.
- **(ii)(c).**
  - dv/ds = −¼(μ_a + μ_b)E′(v) and v³ = 3μ₀K·s re-derive. μ = Gβa′w′ equals μ₀|w| to 2·10⁻⁵ at |w| = 0.57 (μ₀ = 0.1).
  - I re-fitted the optimal dipole and computed K₃ = −⟨r_dip, ∇³φ:p³⟩/24 from an exact series. It equals G's measured K = v²(E − E_inf) to 0.9·10⁻⁷–3.1·10⁻⁶ on five blocks, and to 1.4·10⁻⁴ on target 12.
  - The second-order term is ≤ 1.7·10⁻⁹, against 24K₃ = 2.8·10⁻³–6.1·10⁻².
  - The law also governs the dipole phase **inside the real loop**. On target 9 at ε = 10⁻⁸, over τ = 18–28, d(v³)/dτ = 57.3, 55.6, 55.2, 54.9, 54.8 against 3μ₀KΛ(ε) = 53.9.
- **(iii)/E.** The exceptions are stated honestly.
  - On 18 of 19 stops the error equals the oracle fit's to 7·10⁻⁸, better than the quoted 10⁻⁴.
  - Target 2 is at 1.0486 and 1.2670 times the oracle fit's error.
  - The dipole target never stops.

**What is wrong.**

Majors:
1. **(ii)(b)'s "the fit exact" is the path's stopping tolerance, and H's ends are misdescribed.**
   - The six σ* = 0.06 residuals, 4.3–9.6·10⁻¹⁰, all sit just under `sliding_iss_nd`'s tol = 10⁻⁹.
   - With tol = 10⁻¹³ the path keeps entering kernels up to τ = 1.1·10¹⁵, and the residual only creeps to 8.4·10⁻¹¹.
   - Every σ* = 0.06 end carries 3–6 unreported sub-sample kernels, and five of the six carry a kernel of *negative* width.
   - σ* = 0.04's "coinciding pair" is a zero-sum null triple beside an exact fit of both features.
2. **(ii)(c)'s real-loop sentence.**
   - "Migrates toward a feature" is a midpoint artefact. In 26 of the md's 36 rows the pair's centre of mass already sits on a feature.
   - "Lands … an exact fit" holds only as τ → ∞. The landed pair is (ii)(a)'s pair a little apart, collapsing as d ∝ τ^{≈−0.2} with the residual ∝ Δc⁴.
   - "Lands" versus "in transit" is set by PAIR_TOL.
   - At ε ≥ 10⁻⁶, 8 of 12 cells do not land by τ = 6400 or form no dipole.
   - The ε-dependence obeys a law the draft never states: v_land ∝ ε^{−1/6}, the block's runaway law integrated over an ε-independent window.
3. **The cancellation's stated condition is false, and (b) describes the path, not the dynamics.**
   - "Closer than the kernels' starting width" contradicts its only evidence, which is at 1.33σ₀.
   - My scan finds negative kernels at 0.83σ₀ and 1.33σ₀ and none at 1.67σ₀ or beyond.
   - In the real loop σ* = 0.04 ends with one-sign kernels and an exact fit.

Minors:
- (ii)(a)'s check omits the one stored pair where the prediction fails. Its O(Δc⁴) constant is not uniform in |w_i − w_j|/|W|.
- The mechanism of C's target 7 is misdescribed.
- F's pair does not "stay"; it drifts apart.
- (iii)'s "a coinciding pair adding a noise kernel": under noise no pair forms.
- 20 of 49 mutants survive the tests and `--from-raw`; 17 of them are not equivalent.
- The raw mixes code versions, and C's table still mixes pre-block states.
- (i)'s proof has an O(·) slip.
- The code allows σ ≤ 0.
- "New" still overclaims, and the gradient-flow divergence literature is missing.

Cosmetic: one misquote, and B's nine rows are one computation.

## What I checked and how
1. **Reproducibility.**
   - `--from-raw` into `MW_OUT`, compared with `cmp`.
   - `pytest -p no:cacheprovider tests/test_moving_width.py`.
   - `scripts/rerun_cells.py`: `a_job(3.0)`, `d_job(0.6)`, `h_job(H_TARGETS[0])`, `e_job((2, 1e-3, …))`, each compared field by field with the raw. The only difference is that A's and D's raw lack the round-57 `stalled` key.
   - `scripts/loop_probe.py` and `loop_analyze.py`: G's loop for target 18 at 10⁻⁸ to τ = 25 600, and target 9 at 10⁻⁶ with early snapshots.
   - `scripts/com_vs_mid.py`: all 36 G loop rows.
2. **Raw dumps** (`scripts/raw1.py`–`raw6.py`): A, C, D, E, G and H in full precision, including the negative and sub-sample widths.
3. **Theory.**
   - (i), (ii)(a) and (ii)(c) re-derived by hand.
   - `scripts/split_check.py`: the residual's order and the optimal split, for same-sign and cancelling pairs and for the loop's large-|w| pairs.
   - `scripts/dipole_K3.py`: the optimal-dipole re-fit, K₃ from an exact series for the third directional derivative (checked against differences), and the second-order term.
   - `scripts/mobility.py`: μ/(μ₀|w|).
4. **Path histories.**
   - `scripts/c_pairs_hist.py`: C's three pair targets, block by block.
   - `scripts/h_tol.py`: H "σ* 0.06, equal", block residuals at tol 10⁻⁹ and 10⁻¹³.
5. **One-sign separation scan.** `scripts/h_sep.py` and `h_close.py`: σ* = 0.06, equal amplitudes, at ±0.05, ±0.1, ±0.12, ±0.15 and ±0.2, using `h_job` with the experiment's budget.
6. **Real loop.**
   - `scripts/loop_batch.py`: the six G targets at ε = 10⁻⁵, 10⁻⁷, 10⁻⁹ and 10⁻¹⁰ to τ = 6400, plus the md's 10⁻⁶ and 10⁻⁸. That is 36 cells; each reports the pair's centre of mass, its net W against the amplitudes, and v.
   - `scripts/loop_split.py`: the width split in the loop's pairs.
   - `scripts/loop_dipole_phase.py`: v³ against τ in the dipole phase.
   - `scripts/vland.py` and `exps.py`: the ε-scaling and the collapse exponents.
   - `scripts/h_loop.py`: the loop on H's targets to τ = 4.1·10⁵.
7. **E and F.** `scripts/ef_check.py`: E target 2's stopped blocks, and F target 2 to 150 times the last entry at ε = 10⁻⁶ and 10⁻⁸.
8. **Mutants.**
   - `scripts/mutants.py` and `run_mutants.py`: 49 single-point mutations of `width_split`, `pairs_of`/PAIR_TOL, `closest_pair`, `outcome`, `dipole_law`, the G/H/E report code, `g_job`, `h_job`, `oracle_error`, `oc`/`split_text` and `n_formed`. Each was run against the tests and against `--from-raw` (log in `logs/mutants.log`).
   - `--from-raw` with PAIR_TOL set to 0.005, 0.02 and 0.05.
9. **Literature.**
   - WebSearch snippets only. WebFetch to arxiv.org and www.semanticscholar.org was refused by this session's network allowlist, so no full texts were read.
   - Items marked [memory] were not re-checked.

## Verdicts per part

| part / claim | verdict | evidence |
|---|---|---|
| (i) J, its gradient, one critical point, limits (D59.md:7–12, 45–48) | holds | re-derived; `one_kernel_grad`; B |
| (i) the block's convergence in the continuum (:48–52) | holds, one slip | E ≥ J; compact sublevel set; w bounded away from 0 and ∞; non-degenerate metric; Łojasiewicz. "+ O(E − J)" should be O((E − J)^{1/2}) (D3bh-9) |
| (i) the restriction to the continuum and to a start of the feature's sign (:13–17, :68) | stated correctly | the sign parenthetical holds on the critical curve only (D3bh-9) |
| (i) on the samples: one kernel, exact at nine ratios (A) | holds | raw; the rerun at ratio 3 is bit-identical |
| B "one local minimum at every ratio" | true, but one computation | the grid scales with σ*, and J/σ* is one function (D3bh-12) |
| (ii) C: 71 exact, 3 pairs, 6 dipoles of 80 | holds | raw; C/D outcomes unchanged for PAIR_TOL 0.005, 0.02, 0.05 |
| (ii)(a) expansion, ∂³ identity, split, O(Δc⁴) "for any split" (:52–57) | holds asymptotically | 8 → 16 per halving; the constant grows with \|w_i − w_j\|/\|W\| (D3bh-4) |
| (ii)(a) "checked on the stored pairs" (:70) | 3 of 5 | matches C 1, C 7, D 0.6; C Δ ≥ 0.6 target 2 is off by 8 orders and omitted; D d = 0 is vacuous (D3bh-4) |
| (ii)(a) the mechanism (:22–24) | partly | targets 1 and 2: two early kernels collapse together onto the widest feature; target 7: the last entrant cancels an early kernel's excess (D3bh-5) |
| (ii)(b) close one-sign features make negative kernels on the path | holds | H; also at 0.83σ₀ in my scan |
| (ii)(b) "the fit exact" (:28, :111) | not supported | the stopping tolerance (D3bh-1) |
| (ii)(b) σ* = 0.04 "a coinciding pair" (:111–112) | mislabelled | a null triple beside an exact fit (D3bh-1) |
| (ii) "a fourth … appears when features are closer than the kernels' starting width" (:20–21, :27) | false as stated | H is at 1.33σ₀; none at ≥ 1.67σ₀; 0.83σ₀ is a budget stall (D3bh-3) |
| (ii)(b) "A target of one sign does not keep its kernels of one sign" (:28–29) | fair for the path, not for every loop | the loop on σ* = 0.04 ends one-sign and exact (D3bh-3) |
| (ii)(c) non-closedness, dipole and higher confluent limits | classical, correct | de Boor |
| (ii)(c) dv/ds, μ → μ₀\|w\|, K from the third order, v³ = 3μ₀Ks (:57–66) | holds | K₃ = K to ≤ 3.1·10⁻⁶ (five blocks) and 1.4·10⁻⁴ (target 12); second-order term ≤ 1.7·10⁻⁹ |
| (ii)(c) G's exponents and constants (:98–104) | match the raw | one misquote (D3bh-12) |
| loop: "the dipole is transient" (:37) | holds | dipole phase at τ ≈ 17–29 obeys the law; the third entry ends it |
| loop: "the pair migrates toward a feature" (:38) | not supported | midpoint artefact (D3bh-2) |
| loop: "lands … an exact fit with one kernel too many" (:38–39) | only as τ → ∞ | algebraic collapse; 6 of 28 landed cells frozen (D3bh-2) |
| loop: "when, and on which feature, depends on ε" (:39–40) | holds; the law is missing | v_land ∝ ε^{−1/6}; 8 of 12 cells at ε ≥ 10⁻⁶ do not land (D3bh-2) |
| (iii)/E (:41–42, :112–116) | holds as stated | "a coinciding pair adding a noise kernel" is misworded (D3bh-6) |
| (iii)/F (:43, :116–119) | numbers hold; target 2 misdescribed | the pair drifts apart (D3bh-6) |
| "New" (:81–83) | overclaims | D3bh-11 |

**Status of the round-57 issues (D3bd).**

| D3bd | round-57 issue | status | what remains |
|---|---|---|---|
| 1 | one-sign targets; "separated" undefined | partly | (ii) is restricted to C's distribution, and (b) and H are added; H's ends are misdescribed and its condition is false (D3bh-1, -3) |
| 2 | the loop's "frozen exact" landing; PAIR_TOL | partly | τ to 6400 at two masses, PAIR_TOL = 0.01, "transient" stated; "migrates" and "lands/in transit" are artefacts, the landing is asymptotic, the ε-law is missing (D3bh-2) |
| 3 | the width split | fixed | one stored pair omitted; the constant is not uniform (D3bh-4) |
| 4 | (iii)'s exceptions | fixed | F target 2 "stays" (D3bh-6) |
| 5 | the law's prefactor; κ's order | fixed | K₃ verified independently |
| 6 | the pairs' mechanism; "a gauge" | partly | target 7; "already fitted" (D3bh-5) |
| 7 | (i)'s scope | fixed | B is nine copies of one computation (D3bh-12) |
| 8 | proof gaps and wording | mostly fixed | the O(E − J) slip; the sign parenthetical (D3bh-9) |
| 9 | stalls reported at the wrong state; circular classification | partly | `stalled`/`polish_log` are recorded and used. C's raw was not regenerated (no `stalled`; `oc()` patches it with G's data). C's table still shows pre-block states. The polish that multiplies \|w\| by 592–1125 is recorded, not rejected (D3bh-8) |
| 10 | 25 of 55 mutants survive | partly | 20 of 49 survive, 17 of them not equivalent (D3bh-7) |
| 11 | literature and novelty | partly | references added; the divergence literature, SFW, positivity and the moment algebra are missing (D3bh-11) |
| 12 | misquotes and docstrings | fixed | one new misquote (D3bh-12) |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bh-1 | major | **(ii)(b)'s "the fit exact" (D59.md:28, :111) is the path's stopping tolerance, and H's ends are misdescribed.** (1) The six σ* = 0.06 residuals (4.3, 9.6, 6.0, 8.1, 4.4, 6.1·10⁻¹⁰) all lie under `sliding_iss_nd`'s tol = 10⁻⁹ (rich_transit_nd.py:212, 252–253). In "σ* 0.06, equal", blocks 2–8 add symmetric negative pairs with the residual stuck at 3.9→2.0·10⁻². Block 9 (τ = 6968) jumps to 1.98·10⁻⁹. Two entries at τ ≈ 2.2·10¹⁰ bring it to 4.26·10⁻¹⁰, and the path stops. (2) Rerun with tol = 10⁻¹³, the path continues with 24 entries, 21 kernels and entries up to τ = 1.1·10¹⁵. The residual creeps 2.83→2.32→2.58→1.66→1.12→0.84·10⁻¹⁰ and max\|w\| goes 3.21→3.26. The left feature is fitted by one kernel; the right one is approximated by eight kernels (+3.24 at σ 0.063 and seven negatives at σ 0.050–0.059). It is not an exact fit. (3) Unreported: each σ* = 0.06 end carries 3–6 sub-sample kernels (\|σ\| 8.9·10⁻⁵–1.0·10⁻³, \|w\| 1.1·10⁻⁴–1.8·10⁻²). Five ends hold one or two kernels of *negative* σ (−1.6·10⁻⁴ … −5.8·10⁻⁴), mostly invisible on the samples. The md's "smallest width" counts formed kernels only (exp_moving_width.py:472; mutant H4 survives). (4) σ* = 0.04's "coinciding pair" is a zero-sum null triple (−0.0546/+0.1093/−0.0546 at (0, 0.399)), beside an exact one-kernel fit of each feature (errors ≤ 9.2·10⁻¹²). `outcome` (exp_moving_width.py:179–180) labels any formed pair within PAIR_TOL, whether or not it shares a feature. | (b) and H: say "the path reaches its tolerance 10⁻⁹ with one feature approximated by an eight-kernel cluster plus 3–6 sub-sample kernels". Report the tol = 10⁻¹³ continuation. Add H columns for sub-sample and negative widths over all kernels. `outcome`: count a pair only if its net coefficient and centre of mass fit a feature, else a new class "null cluster". σ* = 0.04: "exact, with an idle null triple". |
| D3bh-2 | major | **(ii)(c)'s real-loop sentence (D59.md:37–40) and numerics (:104–108).** (1) *"Migrates toward a feature".* In 26 of the md's 36 rows the opposite-sign pair's centre of mass is already within 4·10⁻⁸ … 2.6·10⁻³σ* of a feature, with W within 2.2% of its amplitude, by τ = 400. The md's "offset" is the midpoint (exp_moving_width.py:428–430), and it equals the kinematic \|w_a − w_b\|\|Δc\|/(2\|W\|σ*) to ≤ 0.003. For target 18 at 10⁻⁸ the md has 0.210/0.161/0.125, against centre-of-mass offsets 2.6·10⁻³/1.2·10⁻³/5.9·10⁻⁴. (2) *"Lands … an exact fit".* After landing, the pair is (ii)(a)'s pair a little apart: its split is within 0.04–10% of (w_j − w_i)Δc²/(6Wσ). It collapses algebraically, d ∝ τ^{−0.17…−0.24}, with the residual ∝ \|Δc\|^{3.9–4.1}. Target 18 at 10⁻⁸ is still at d = 9.5·10⁻³ and rms 6.2·10⁻⁸ at τ = 25 600. Over ε = 10⁻⁵…10⁻¹⁰, only 6 of 28 landed cells are frozen: target 9 at ε ≤ 10⁻⁶ and 27 at 10⁻¹⁰, all landed at d ≤ 1.1·10⁻³. Only 8 of 28 are below 10⁻⁸ at τ = 6400. "Lands" for 27 at 10⁻⁸ versus "in transit" for 12, 18 and 29 at 10⁻⁸ is PAIR_TOL = 0.01; 18 crosses it at τ ≈ 2·10⁴. (3) At ε ∈ {10⁻⁵, 10⁻⁶}, 8 of 12 cells do not land by τ = 6400 or form no dipole. Targets 9, 12 and 18 at 10⁻⁵, and 12 and 27 at 10⁻⁶, have five formed kernels with W matching no amplitude (18 at 10⁻⁵ is 10% off) and rms 2.7·10⁻⁶–6.2·10⁻⁴. Targets 20 at 10⁻⁵ and 10⁻⁶, and 29 at 10⁻⁵, form no dipole. At ε ≤ 10⁻⁷ all 24 cells land by τ ≈ 100–400. (4) *Missing law.* The landed half-difference grows as ε^{−1/6}. On target 9, v = 4.24, 6.41, 9.47, 13.9, 20.4 at ε = 10⁻⁶…10⁻¹⁰, i.e. per decade 1.51, 1.48, 1.47, 1.47; on target 20 per decade 1.467, 1.472, 1.471. The reference is 10^{1/6} = 1.468. This is because the loop's dipole phase follows the block's law in t = Λ(ε)τ (see the Summary). | (ii)(c): describe three phases. First, a dipole phase with v³ = 3μ₀KΛ(ε)(τ − τ_d) until the next entry. Second, the landing by τ ≈ 100–400, with the centre of mass on a feature and W → A. Third, an algebraic collapse in (ii)(a)'s regime (d ∝ τ^{≈−0.2}, rms ∝ d⁴), exact only as τ → ∞. State v_land ∝ ε^{−1/6} and test it at three or more masses. List the non-landing cells at ε ≥ 10⁻⁶. G: report the centre-of-mass offset and W − A instead of the midpoint, and compute the pair column at report time. |
| D3bh-3 | major | **The cancellation's condition (D59.md:20–21, :27; md:143) is contradicted, and (b) is the path's.** The only evidence (H) is at ±0.08, i.e. 0.16 = 1.33σ₀ apart (:109), not "closer than σ₀". My scan at σ* = 0.06 with equal amplitudes: <br>• 0.83σ₀: 22 entries and 10 negative kernels (down to −0.24); stopped by the 10⁶ budget at rms 6.5·10⁻⁶, a non-dipole stall. <br>• 1.33σ₀: H. <br>• 1.67σ₀ and 2σ₀: exact, with no formed negative. <br>• 2.5σ₀ and 3.33σ₀: coinciding pairs, no negative. <br>With A₂ = 0.6, H itself is exact. In the real loop: <br>• σ* = 0.04 ends with the two exact one-sign kernels by τ = 6400 (w 0.995 and 0.998 at ε = 10⁻⁶ and 10⁻⁸); the 6–8 negatives present at τ = 1600 are gone. <br>• σ* = 0.06, equal, keeps negatives at τ = 4.1·10⁵, but as wide cancelling kernels around one positive kernel per feature (w 0.988 at 10⁻⁸), with centres out to ±1.80, outside the samples' window ±1.6. The rms is 1.9·10⁻⁶ and still falling. <br>• A₂ = 0.9 keeps the path-like cluster and five formed sub-sample spikes (\|w\| ≤ 0.115). | (ii) and (b): replace the condition with the observed window: separations ≈ 0.8–1.3σ₀ and amplitude ratio ≥ 0.9 give negatives; none at ≥ 1.67σ₀ or at ratio 0.6. Add the scan to H. State (b) as the path's outcome and add the loop's ends. |
| D3bh-4 | minor | **(ii)(a): selective check, non-uniform constant.** "Checked on the stored pairs" (:70; numerics :92–94) quotes C's targets 1 and 7. It omits C Δ ≥ 0.6 target 2 (md:48: observed +3.2834·10⁻⁸, predicted +1.8958·10⁻¹⁶). That pair is a σ-split at Δc = 2.8·10⁻⁸, unresolved at the residual's 10⁻¹⁵ floor. D d = 0 has equal shares, so the check is vacuous there. "For any split" (:56) is asymptotic, with a constant that grows with \|w_i − w_j\|/\|W\|. For (−9.18, 9.75) at σ 0.144 (target 9's landed pair at 10⁻⁸), the predicted split worsens the fit at Δc = 0.02 (rms 0.30 against 0.10). The 16-per-halving rate appears only below Δc ≈ 0.005, where \|w_i − w_j\|Δc/(\|W\|σ) ≈ 1.2. Above that, the optimal split deviates 8–11% from the formula. REVIEW_LOG.md:2177 says "to 4–5 digits"; D 0.6 agrees to 6·10⁻⁴. | (a): "on the three pairs a little apart; at coincidence the split is unresolved". State the validity condition \|w_i − w_j\|Δc ≪ \|W\|σ; it is what separates (a) from the dipole (c). Add the loop's landed pairs as a further check. |
| D3bh-5 | minor | **The pairs' mechanism (D59.md:22–24).** "A kernel of an earlier block collapsing at the last block onto a feature already fitted": <br>• Target 1 (Δ ≥ 0.4): kernels 13 and 14 sit at (0.212, 0.095) and (0.482, 0.182) through blocks 2–3. Together they approximate the widest feature (0.351, 0.273), which neither fits. They collapse onto it at block 4. Target 2 (Δ ≥ 0.6) is the same, with kernels 10 and 11. <br>• Target 7: kernel 3 (block 1) lands on the −0.43 feature with −0.786. The *last entrant* (12, +0.357, entered at τ = 105.4 from 0.2) cancels the excess. That feature was not fitted before. | Describe what was observed: "two early kernels approximating one wide feature collapse onto it together, or the last entrant, of the other sign, cancels an early kernel's excess". |
| D3bh-6 | minor | **F and (iii) wording.** <br>• F target 2 "stays 7.0–8.8·10⁻³ apart" (:118–119): the md's own numbers grow from 7.0 to 8.8·10⁻³. At ε = 10⁻⁶ the growth continues: 1.17·10⁻² at 30 and 1.57·10⁻² at 150 times the last entry (\|Δc\| 1.5→4.5·10⁻³). The shares are frozen at 0.329/0.489 against the path's 0.231/0.587. The residual, 5.2→4.5·10⁻⁷, is the idle kernels' floor (their \|w\| ≈ 5·10⁻⁶). At 10⁻⁸ the pair sits at 2.57·10⁻³ with residual 2.1·10⁻¹⁰. <br>• (iii)'s "a coinciding pair adding a noise kernel" (:42): under noise target 2 forms no pair. Its fourth kernel is noise: (−0.129, 0.264) with w 0.022 at s = 10⁻³, and (−1.315, 0.015) with w 0.013 at 10⁻². | F: "drifts apart at 10⁻⁶ (residual set by the idle kernels); 2.6·10⁻³ at 10⁻⁸"; add 10⁻⁸ and longer times. (iii): "on the coinciding-pair target the stop adds a noise kernel (no pair forms)". |
| D3bh-7 | minor | **Mutants.** 49 mutants; 29 killed (10 by the tests only, 8 by `--from-raw` only, 11 by both). 20 survive both. P1 (`<=` PAIR_TOL), C1 (`> 0` for `>= 0`) and E1 (10⁻² in E's count) are equivalent on the stored data. The other 17: <br>• C4: `closest_pair` without its width term. <br>• `outcome`: O1 (dipole threshold 10→100·max\|A\|), O2 (distance 0.05→0.5), O3 (stalled pair of either sign), O4 (pair rms 10⁻⁸→10⁻⁶), O7 (formed at 10⁻³), O8 (pairs among all kernels). <br>• `dipole_law`: D4 (E_inf without ½), D5 (b the farthest opposite kernel). <br>• `g_job`: G6 (offset in S0 units), G7 (wmax over min\|A\|), G8 (stalled pair of either sign). <br>• `h_job`: H1 (negative threshold 10⁻¹), H2 (`max` for `min`), H3 (counts positives), H4 (smin over all kernels). <br>• E4: `oracle_error` started from 1.1·A. <br>The job functions are never run by a test or by `--from-raw`. `outcome`'s stalled branch never sees real data, because C's raw has no `stalled`. | tests/test_moving_width.py: `outcome`'s stalled branch at its three thresholds; its pair rule with formed and unformed kernels and with a null cluster; `closest_pair` with unequal widths; `dipole_law` with two opposite kernels and a known nonzero E_inf; `h_job`, `g_job` and `oracle_error` on cheap synthetic cases (e.g. a 2·10⁴ budget, S_G = (1e3, 1e5), TAUS_G = (25,)). |
| D3bh-8 | minor | **Raw provenance and C's table.** A–D in the raw predate round 57's `match` (no `stalled` key). C's dipole class therefore comes from the report's `oc()` patch with G's data (exp_moving_width.py:509–517). C's table still prints pre-block states: [1, 4] formed kernels and a 2.7·10⁻¹ maximum residual (md:43–44). G's "coinciding pairs" column is computed at run time, while the md header prints the current PAIR_TOL: with 0.005, 0.02 or 0.05, only the header text changes. The stalled blocks' \|w\| (1431–3964) are the accepted polish's, which multiplied \|w\| by 592–1125 (md:93–98). That state is solver-dependent. My reruns of A at 3 and D at 0.6 are bit-identical, so no numbers should move. | Regenerate A–D with the current code and drop `oc()`. Print the stalled block in C's table. Compute G's pair column in `report`. Mark the polished \|w\| as solver-dependent, or reject polishes that grow \|w\| more than about 10 times. |
| D3bh-9 | minor | **(i)'s proof.** "\|w\| = \|⟨f*, φ⟩\|/‖φ‖² + O(E − J)" (:50) is wrong in order. E = J + ½‖φ‖²(w − w*)², so \|w − w*\| = (2(E − J)/‖φ‖²)^{1/2}. The lower bound follows from E ≤ ½‖f*‖² − η, which gives \|w − w*\|² ≤ w*² − 2η/‖φ‖². "On the free law's curve a coefficient keeps its sign" (:14) holds on the critical curve (b0 = 0). On the loop's curves (b0 = ε), α̇ = Gε⟨rφ⟩ ≠ 0 at α = 0. | Correct the order and add the one-line bound. Say "on the critical curve". |
| D3bh-10 | minor | **σ is unconstrained.** `BumpsCW` has no σ > 0 (moving_width.py:28–32), and `polish_block_nd` bounds only the coefficients (rich_transit_nd.py:195–196). Five H ends hold negative widths. `match`, `pairs_of`, `width_split` and `outcome` use \|σ\| without saying so. | Use log σ coordinates or bound σ in the polish; at minimum count and report σ ≤ 0 and σ below the sample spacing. |
| D3bh-11 | minor | **Literature and "New" (:72–83).** <br>• Missing: gradient flows that diverge toward unattained infima (Gallon–Jentzen–Lindner 2022; Kranz–Gallon–Dereich–Jentzen, NeurIPS 2025). <br>• The exponent 1/3 comes from the coefficient's mobility μ ∝ \|w\|. This is the Fisher–Rao/mirror-descent weight mobility of conic particle descent (Chizat 2022; [memory] for the mobility form). A Euclidean mobility gives v⁴ = 4μKs; generally v^{a+2−b} ∝ s for E − E_inf ∝ v^{−a} and μ ∝ v^b. So the law is a two-line consequence once E − E_inf = K/v². <br>• The width split is the vanishing third central moment of the pair seen as a signed two-component mixture. That is Pearson's moment algebra; Runnalls' moment-preserving merge is its second-order part. <br>• One-sign cancellations are what positivity removes: positive sources with a Gaussian PSF are recoverable without separation (Schiebinger–Robeva–Recht 2018). <br>• Sliding Frank–Wolfe (Denoyelle–Duval–Peyré–Soubies 2020) belongs beside CL-OMPR and ADCG. | Add these. Narrow "New" to the loop's three phases and the ε^{−1/6} landing law (D3bh-2). Claim the width split and the cancellations only as observed in these flows. |
| D3bh-12 | cosmetic | <br>• "To 2·10⁻⁴ on target 12" (:103): Δv³/Δs is 6.31839·10⁻⁵ against 3μ₀K = 6.31681·10⁻⁵, i.e. 2.5·10⁻⁴. <br>• B's nine ratios (:90–91) are one computation, because the grid scales with σ*. <br>• E's "to 10⁻⁴" (:115) understates: the agreement is 7·10⁻⁸. <br>• md:143 repeats "closer than the kernels' starting width". | Edit as stated. |

Sources (search snippets; full texts not fetched):
- de Boor, γ-polynomials: https://ftp1.cs.wisc.edu/Approx/gamma.pdf
- Petersen–Raslan–Voigtlaender, FoCM 2021: https://link.springer.com/article/10.1007/s10208-020-09461-0
- Kranz–Gallon–Dereich–Jentzen, NeurIPS 2025: https://papers.neurips.cc/paper_files/paper/2025/hash/05fb0f4e645cad23e0ab59d6b9901428-Abstract-Conference.html
- Gallon–Jentzen–Lindner 2022: https://arxiv.org/pdf/2211.15641 and https://thesis.unipd.it/handle/20.500.12608/10206
- Chizat, Math. Program. 2022: https://www.deepai.org/publication/sparse-optimization-on-measures-with-over-parameterized-gradient-descent
- Schiebinger–Robeva–Recht: https://ar5iv.arxiv.org/html/1506.03144
- Sliding Frank–Wolfe: https://bigwww.epfl.ch/publications/denoyelle2001.html
- HyperBF: https://en.wikipedia.org/wiki/Hyper_basis_function_network and https://dspace.mit.edu/handle/1721.1/6014
- Pearson's two-Gaussian moments: https://blog.mrtz.org/2014/04/22/pearsons-polynomial.html
- Runnalls 2007: https://www.cs.kent.ac.uk/pubs/2007/2797
- Varah 1985 (bibliographic only; whether it treats coalescing rates was not checked): https://epubs.siam.org/doi/10.1137/0906003
- Not re-checked [memory]: Kammler, Braess, Marr–Hildreth, CL-OMPR, ADCG, Kiefer–Wolfowitz, Day, Biernacki–Chrétien, Duval–Peyré, Golub–Pereyra.

## Lesson
Each of the draft's new words was produced by a threshold, not by the dynamics.
- "Exact" was the path's stopping tolerance.
- "Coinciding pair" was PAIR_TOL applied to a null triple.
- "In transit" versus "lands" was PAIR_TOL applied to a distance still falling as τ^{−0.2}.
- "Migrates" was a midpoint that a pair of unequal coefficients displaces kinematically.

Before naming a configuration, measure its invariants: the centre of mass and net coefficient against the feature, how the residual scales with the separation, and whether the end survives a smaller tolerance or a longer horizon. Then ask whether the law derived for one stage, the block's runaway, predicts what is seen downstream. Here it does: the loop's dipole phase follows v³ = 3μ₀KΛ(ε)τ, and the landed coefficient scales as ε^{−1/6}. That is the quantitative content the proposition is missing.
