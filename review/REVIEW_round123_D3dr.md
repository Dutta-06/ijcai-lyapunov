# Referee report D3dr, round 123: Proposition D‴.58 as revised after round 118 and D‴.61's part VA

## Scope

- **Subject.** Commit 783d55b on claude/vigilant-clarke-hf939b, diffed against a2f3b77. That range holds two D‴.58 commits:
  - b2088e8, round 118 applied: A6X's verdict by successive orders and their last change, the transient derived and tested by B_K, its next order, the resolvability trigger, the premise.
  - 783d55b, the tiers after D‴.61's part VA.
- **Work directory.** W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref123`.
  - Export: `git archive 783d55b | tar -x -C W/repo`, with `W/repo/venv` a symlink to the repository's venv.
  - WM_OUT = W/out and GP_OUT = W/gpout, each holding a copy of its raw.
  - Scripts are in W/own, logs in W/logs, mutant slots and results in W/mut.
- **Read in full.**
  - D58.md, and `git diff a2f3b77 783d55b` for D58.md, exp_width_mismatch.py, width_mismatch.md, the tests and the two mutant specs.
  - In exp_width_mismatch.py: `a6x_lines`, `a6x_resolvable`, `a5_series_lines`, `a6_lines`, `a6e_lines`, `d61_step`, and the constants D61_RSTAR, _V, _VA and D61_STATUS.
  - The tests of rounds 101, 107 and 118: `test_round101_*` (3), `test_round107_d61_status_pinned_to_its_raw_and_the_a6x_verdict`, `test_round107_a6x_verdict_plants`.
  - tools/mutants/round118_D58_own.py and round121_D58_tiers.py.
  - REVIEW_round118_D3dm.md; REVIEW_LOG Rounds 118, 120 and 121; REVIEW_round121_D3dp's log entry.
  - D61.md (i)–(v), its proof of (iv) and its Status; gh_peano.md's sections on (iv), V and VA; HANDOFF §5.
- **Run settings.**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv, at most two heavy processes of mine at a time, pytest with `--basetemp` inside W.
  - My PIDs: 13088 (`a6_job` for (2, Mercer frame)), 17300 (`a6_job` for (1, Mercer frame)), 14551/14552 (mutant slots), 17787/18102 (mutant pass 2).
  - I killed nothing. 17787 and 18102 were still running at handback; they write only into W/mut and exit on their own.
- **The repository was not modified.** I ran only `git archive/log/diff/status/show`.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my code, or the repository's functions on inputs I changed.
  - [snippet]: web-search results only; arxiv.org and www.ams.org do not resolve in this session.
  - [memory]: recalled facts.

## Summary

**Design.**
- Regenerated the md and ran the 63 tests. Reran D‴.61's part VA from the repository's code and compared it with the raw.
- Recomputed every A6X quantity independently: t(v) in Arb at 600 bits, the fits in my own exact-rational elimination.
- Bisected the new verdict's acceptance windows through the repository's own `a6x_lines`.
- Calibrated the verdict on known cases. I ran `a6_job` at K = 200, 240, 320 for the two proved Mercer-frame rules (τ/σ = 1 and 2; 413 s and 108 s).
- Fed the verdict 30 exact mock models, including models that should fail. Among them: offsets, K⁻¹, K^{−5/2} and log K/K² terms, unmodelled geometric components, a transient with a K^γ prefactor, B < 0, and interference at √q.
- Checked the tier claims against D‴.61 (iv) word by word, recomputed every threshold in closed form, and tested "every t(v) < 1" numerically at a frame rule with t(v) = 0.9945, beyond V's tier.
- Mutants: the repository's 10 (T01–T06, V01–V04) plus 36 of my own (X01–X36), with a second pass that excludes the md regression test.
- Five web searches; two fetches, both blocked.

**Verdicts.**
- **Reproducibility is exact.**
  - The md is byte-identical (md5 74245f51…). The raw is unchanged (f1079628…).
  - 63/63 tests pass in 51.6 s.
  - Part VA rerun from `gh_lastann` (2.8 s) gives the raw's 'VA' record for record, accepted boxes included.
- **Every new number in D58 and the md is right** to the digits printed [own].
- **The tiers (b) follow from D‴.61 (iv) exactly as D‴.61 states it.**
  - Every frame has t(v) < λ(σ*/σ) < 1.
  - The half-squeeze frame ρ = t exists at every σ*/σ.
  - (iv) needs uniform convergence only on |η| ≤ c < 2, a compact subset.
  - Nothing claims uniformity in t; the frameless rules (ρ = 1, |ρ| > 1) and A6X stay numerical/conjectural.
  - Numerically [own]: at τ/σ = 20, v = 0.9τ² (ρ = 0.99944, t(v) = 0.99447) the ratio is 1.005639 at K = 160 against 1/t(v) = 1.005557, and falling.
- **Round 118's fixes are in.** The trigger no longer misfires on the proved Mercer frame (D3dm-4).
- **A6X's new verdict is better founded but still weakly calibrated.**
  - It is global over rules.
  - It fails both proved Mercer-frame rules.
  - It cannot tell a premise failure from a wrong limit.
  - A misspecified power can pass it while hiding an offset 290 times its printed resolution.
  - Its trigger and its orders look at different K.
- **No major issue.**

**Findings.**
1. **D3dr-1 (minor).** The transient trigger and the orders use different K. `a6x_resolvable` tests q^160, but the seven- and eight-term orders reach K = 60 and 40. A geometric component with 0.3 ≲ q < 0.75 is neither modelled nor kept out of them, and "decrease" fails.
2. **D3dr-2 (minor).** On known cases the verdict fails 2 of the 3 proved fast-branch rules extended to K = 320, and one rule's failure prints the failure for all rules. "Checked rule by rule" is not what the code does.
3. **D3dr-3 (minor).** The orders' check has limited power and misattributes its failures.
   - A K^{−5/2} term passes and hides a 7.6·10⁻¹² offset behind a printed 2.6·10⁻¹⁴.
   - A K⁻¹ term with the limit exact prints "a limit off 1/t(v)".
   - Exact data at the rounding floor fail "decrease": there is no noise floor.
4. **D3dr-4 (minor).** The transient checks need restating and tightening.
   - B_K's flatness "over K = 60 … 160" is two points: B_120 = B_160 = B identically, since they are nodes of the transient fit.
   - B > 0 is printed, not tested.
   - `tr_ok` accepts outlier offsets of ±6–7·10⁻¹⁵, against the stated resolution of about 2–3·10⁻¹⁵ (D3dm-2's request).
5. **D3dr-5 (cosmetic).** The tier upgrade rests on D‴.61's adoption of round 121's fix, whose own round (122) is in flight. Two attributions are missing: round 120 for VA, and "Arb's correctness assumed".
6. **D3dr-6 (cosmetic).** Residues: "25-digit" in the md and docstring; a stale `a6x_lines` docstring; the md heading without round 118; "calibrated" in a test docstring; two roundings of one quantity; "about 2·10⁻¹⁵" not marked as unreproduced; `d61_step` for t ≥ 1.
7. **D3dr-7 (cosmetic).** Tests: the repository's 10 mutants are all killed. 15 of my 36 survive, none changing the md. The thresholds of `mono`, `cons`, `flat`, `tr_ok` and the trigger are unpinned.

## What I checked and how

**1. Reproduction** [md, raw, own]
- `--from-raw` into W/out (1.3 s), then `cmp` with results/chaos/width_mismatch.md: identical.
- Full test file: 63 passed in 51.56 s.
- `W/own/rerun_va.py` calls `GP.part_VA()` and compares with the raw's 'VA'. All 12 records are the same, the accepted η-boxes included: validate; A-J 40; A-FC 2011; B-J 320 (634 skipped); B-F 769 (1319 skipped); O 76 boxes and 5193 segments, worst −1.4993·10⁻⁷; end −1.3314; handoff 1140 points, tiled; soundness 80/5358/640/1538 points, none failed. The counts are round 121's patched ones.

**2. A6X recomputed independently** (`W/own/recompute.py`: t(v) and q in Arb at 600 bits from my own v's, fits in exact rationals) [own, raw]
- **The outlier (¼, 2·v_hs).**
  - Orders (four … eight terms): −6.5193·10⁻¹³, +9.3071·10⁻¹¹, −2.8792·10⁻⁹, +3.7491·10⁻⁹, +1.1725·10⁻⁸.
  - Transient fit +6.5243·10⁻¹⁶, B = 0.238302; next order (K⁻⁵, from K = 80) +3.1066·10⁻¹⁵.
  - B_K/B − 1 at K = 60 and 80: −1.18·10⁻⁴ and −2.5·10⁻⁵. At K = 120 and 160 it is exactly 0, because the fit passes through those K. At 40 and 20: −1.65·10⁻³ and −3.8%.
  - q = 0.85062954.
- **(½, MM):** −3.861·10⁻¹⁴, −2.092·10⁻¹⁶, +1.909·10⁻¹⁷, −4.829·10⁻¹⁹, −7.406·10⁻²⁰.
- **(1, MM):** +8.612·10⁻¹⁵, −3.267·10⁻¹⁶, +2.084·10⁻¹⁸, +7.741·10⁻¹⁹, −1.183·10⁻¹⁹.
- **(2, MM):** −3.568·10⁻¹⁴, +5.185·10⁻¹⁶, −3.641·10⁻¹⁸, −2.913·10⁻²⁰, −1.049·10⁻²⁰.
- **(2, 2·v_hs), proved:** −1.985·10⁻¹⁴, +4.947·10⁻¹⁶, −7.429·10⁻¹⁸, +1.885·10⁻¹⁹, −3.160·10⁻²⁰.
- **Input propagation.** From the 26-digit errors (8.3·10⁻²⁶ each, Σ|w| = 32, 61, 89, 117, 139) the orders carry at most 2.7·10⁻²⁴ … 1.2·10⁻²³, 900 times below the smallest |o₈|.
- **Margins on the data.** The last-pair ratios |o₈/o₇| are 0.15, 0.15, 0.36, 0.17. `cons`'s |o₈|/(10|o₈ − o₇|) is 0.018, 0.013, 0.056, 0.014.
- Every printed value matches: 6.5/6.6·10⁻¹⁶, 3.1/3.2·10⁻¹⁵, 1.2·10⁻⁴, 1.2·10⁻¹⁹, the orders, 18.6–5.63·10⁴, 21.11/10/2.533/2.
- **The derivation in D58 :260–262 is right.** With ε_K = (b_K/a_K)² = Cq^K(1 + O(1/K)),
  R = R^{(a)}[1 + (ε_{K−1} − ε_K)/2 + O(ε²)] = R^{(a)} + (C(1 − q)/(2tq))·q^K(1 + O(1/K)) + O(q^{2K}),
  so B > 0. The data give C = 2tqB/(1 − q) = 0.94.

**3. The verdict's acceptance windows** (`W/own/windows.py`, `scan.py`: the repository's `a6x_lines`, a relative offset d planted in one rule's ratios and e₀, as the tests' `pl` does; scanned on a grid first, so the accepted set is an interval) [own]

| rule | `mono` ∧ `cons` accept d in | `mono` alone |
|---|---|---|
| (½, MM) | −4.0·10⁻¹⁸ … +2.8·10⁻¹⁹ | −9.3·10⁻¹⁸ … +2.8·10⁻¹⁹ |
| (1, MM) | −3.3·10⁻¹⁹ … +9.0·10⁻¹⁸ | … +1.6·10⁻¹⁶ |
| (2, MM) | −1.8·10⁻¹⁹ … **+2.0·10⁻²⁰** | −2.6·10⁻¹⁶ … |
| (2, 2·v_hs) | −7.8·10⁻²⁰ … +2.2·10⁻¹⁸ | … +3.6·10⁻¹⁸ |
| outlier (`tr_ok`) | −7.2·10⁻¹⁵ … +5.8·10⁻¹⁵ | — |

- Each edge sits at that rule's own |o₈| scale, and the printed "at most" moves with the offset. The criterion is self-consistent per rule, not post hoc as round 107's and 113's were.
- The windows are asymmetric by up to 10⁴.
- `cons` accepts a limit off by up to 10 times the last change. A mock offset of 10⁻¹⁷ passes with a last change of 4·10⁻¹⁸.

**4. Known cases: the proved Mercer frames** (`W/own/extra_job.py` = the repository's `a6_job`; `known.py`) [own]
- **(2, Mercer frame), limit 7/3.**
  - Ratios 2.3333306963519576, 2.3333315103923598, 2.333332313711865.
  - Orders +1.2·10⁻¹², +2.5·10⁻¹³, −1.3·10⁻¹³, +7.9·10⁻¹⁴, −1.4·10⁻¹³. This reproduces round 118's stall.
  - Through 200–320 it deviates by −9.69·10⁻¹¹, which is 28.9 times (2, 2·v_hs)'s −3.35·10⁻¹². That reproduces round 113's "29", which the md marks as not reproduced.
- **(1, Mercer frame), limit 5.309401** (proved, D‴.61).
  - Ratios 5.3094092931583425, 5.309406712202298, 5.309404198357432.
  - Orders −7.0·10⁻¹¹, +8.5·10⁻¹², −2.2·10⁻¹², +1.2·10⁻¹², −1.7·10⁻¹²: a stall near 10⁻¹².
  - Excess × K² drifts: +0.0631, +0.0619, +0.0611, +0.0602.
- **What the verdict prints** with either rule added: "the orders do NOT decrease in magnitude through eight terms at the 5 [6] rules without a transient", beside the hardcoded caveat "(at the proved τ/σ = 2, Mercer frame they stall …: not reproduced by the repository's scripts)".
- With q = 0.605 and 0.675 the trigger does not fire. Round 118's misfire (B = −1.1·10¹⁵) is gone.

**5. Models that should, and should not, pass** (`mocks.py`, `mocks2.py`, `mocks3.py`: exact 60-digit mock errors through `a6x_lines`) [own]
- **The base** 4(1 + 1/(K + 1)²): MONO CONS (last +1.1·10⁻¹⁹).
- **Offsets.** +10⁻¹⁴ is caught (CONS fails). +10⁻¹⁷ passes (last +1.0·10⁻¹⁷). ±10⁻²⁰ pass.
- **A K⁻¹ term** 10⁻⁶, 10⁻⁹, 10⁻¹², limit exact: the orders stall at 8.2·10⁻¹⁰, ·10⁻¹³, ·10⁻¹⁶. MONO without CONS prints "a limit off 1/t(v)", for an exact limit.
- **A K^{−5/2} term** 10⁻³, limit exact: −1.8·10⁻¹¹ … −7.6·10⁻¹², MONO CONS.
  - **With a +7.6·10⁻¹² offset on top:** −1.1·10⁻¹¹, −5.3·10⁻¹², −2.5·10⁻¹², −8.5·10⁻¹³, **+2.6·10⁻¹⁴, MONO CONS.** The verdict would print "1/t(v) to the orders' resolution" at 2.6·10⁻¹⁴ for a limit 7.6·10⁻¹² off.
  - The same term at 10⁻⁶ is not mono.
- **log K/K²** at 10⁻³ and 10⁻⁶: MONO without CONS ("a limit off").
- **A finite power expansion, exact:** the orders sit at the 10⁻⁶¹ rounding and the verdict prints not mono. It has no noise floor.
- **An unmodelled 0.24·q^K, limit exact:** q = 0.30, 0.45, 0.55, 0.60, 0.675, 0.74 give o₈ = −5.7·10⁻²⁷, −6.3·10⁻²⁰, −1.9·10⁻¹⁶, −6.2·10⁻¹⁵, −6.7·10⁻¹³, −2.0·10⁻¹¹, all not mono. The trigger fires only for q > 0.7499.
- **The outlier's geometry** (t(v) and q exact, a smooth power part, and Bq^K):

  | model | `tr_ok` |
  |---|---|
  | γ = 0 | passes (dev 2·10⁻⁶⁰) |
  | γ = 10⁻⁴, 10⁻³ | passes (dev −1.9·10⁻¹⁷, −1.9·10⁻¹⁶) |
  | γ = 0.01, 0.1 | fails |
  | B < 0 | **passes** |
  | interference c·√q^K, c = 10⁻⁸ | fails |
  | interference, c = 10⁻¹⁰ | passes |
  | a −Bq^{2K} term | passes |
  | offset 10⁻¹⁵, 5·10⁻¹⁵ | passes |
  | offset 10⁻¹⁴ | fails |
  | (1 + 5/K) prefactor | fails |
  | a second geometric component 0.1·0.80^K | fails |

**6. The tiers against D‴.61 (iv)** [md, own]
- **What D61 (iv) states.** "For every r < 1 with VA's … the rate t(v) and the sharp asymptotics hold for every scaled rule with a frame … with VA for every t(v) < 1", and the half-squeeze constant "for every σ*/σ with VA's".
- **Every frame is covered.** |ρ| < 1 ⟺ v ∈ (τ²σ²/σ*², τ²). On that interval t(v) ∈ [t, λ) with λ = tanh r < 1, so r = t(v) < 1. D61's tails need c² > 4r ln(1/r)/(1 − r) and its mass needs 4r: c = 0.99983, 1.89270, 1.99498 and 1.99950 at r = 0.0966, 0.8052, 0.99 and 0.999, all inside the open disc. (iv) uses (iii)(b) only on |η| ≤ c < 2, a compact set.
- **The half-squeeze frame.** ρ = t, v′ = 2 and γx = γy = 1/(2t) (the cusp) hold at every σ*/σ.
- **The Fock and Mercer-zero rules.** ρ = 0 exactly, and ρ = 1 − 2(√(1 + 2x) − 1)/x ∈ (−1, 1), so both have a frame at every τ/σ.
- **Uniformity.** The claims are pointwise in t. D58 claims nothing uniform as t → 1, rightly, and none is available.
- **Thresholds, recomputed in closed form.** Fock 19.924859/19.899749, Mercer-zero 140.712473, frame end 14.106736, branch and bound 4.187348/4.066188, 6.514453 and 3.044165 (printed rounded up, 3.045), half-squeeze 199.00 and 9.2669. The md's 19.924, 19.899, 140.712, 14.107, 4.187, 4.066, 6.514, 3.045 and 199 are all right.
- **What stays unproved.** D58 keeps the frameless rules (4 moment matching, 10 with |ρ| > 1) conjectural and A6X numerical.
- **Numerically at t(v) → 1** (`near1.py`, the repository's `a6_error` at 120 digits) [own].
  - τ/σ = 20, v = 0.9τ² (t(v) = 0.994474, VA's tier only): ratio 1.031623, 1.017261, 1.010085, 1.006690, 1.005639 at K = 10 … 160, against 1.005557.
  - The half-squeeze rule at t = 0.904875: ratio over (1 − 1/(4K))/t, minus 1, falls to +1.6·10⁻⁶; e_K over the sharp prediction is 1.00448, 1.00047, 0.99916 at K = 40, 80, 160 (pre-asymptotic, as D58 says at t = 0.8).

**7. Tests, line by line** [md, own]
- **`test_round107_d61_status_…`.**
  - It pins VA's raw: every piece with failed = 0, O's 76 boxes from R1 = 1.995, the end, the handoff tiled, four soundness pieces.
  - It pins D61_RSTAR_VA = 1, the full verdict string, B_K's "1.2e-04", "+3.1e-15", and `a6x_resolvable` at 0.850630 (yes) and 0.605 (no).
  - Plants: 10⁻¹⁵ at (½, MM) breaks "decrease"; K = 60's e₀ moved alone breaks flatness.
  - Its docstring still reads "calibrated on the raw … is NOT calibrated".
- **`test_round107_a6x_verdict_plants`.** The 40-digit mock passes; a 10⁻¹² offset gives "NOT each within ten times"; a frame rule at t(v) = 0.995 is "the proved rule". It has no transient mock and no threshold edges.
- **`test_round101_a5_label_compared`.** `d61_step` at 0.5, 0.9 and 0.995; the tiers lowered by monkeypatch; a frame at t(v) = 0.995 proved, ρ = 3 − 2/1.2 not.
- **`test_round101_outward_…`.** The A6 header, the Fock/Mercer-zero sentence, and the a5 sentence.
- **The c_need check** is a grid evaluation of an analytic inequality (ln x < x − 1). It is harmless, but it is not a check on the raw.

**8. Mutants** (`W/own/spec123.py`; runner `W/own/mutants123.py` = tools/mutants.py with REPO, PY and a per-slot `--basetemp`; `--extra gh_peano_raw.npz`. Pass 2 reruns the md-changing kills with `-k "not round88_precision_and_the_md"`; the md regression is test 54 of 63, so "1 failed, 53 passed" means the md regression killed it) [own]
- **The repository's specs.** T01–T06 and V01–V04 are all killed; the controls T00 and V00 pass. T05, T06, V03 and V04 were first killed by the md regression; their pass 2 had not finished at handback.
- **Mine, killed (21).**
  - By unit tests, md unchanged: X02, X04, X08, X21, X25.
  - Changing the md and killed again with the md regression excluded: X03, X05, X06, X09, X11, X17, X29.
  - Changing the md, pass 2 not finished: X12, X18, X26, X28, X30, X31, X32, X33, X34.
- **Mine, surviving all 63 tests with the md unchanged (15).**
  - X01 (`d61_step`'s ≤) and X07 (`d61_step` on the unrounded t(v)): equivalent unless a t(v) falls on a tier.
  - X10 (the frame tolerance dropped): equivalent at dyadic τ, where moment matching's float ρ is exactly 1. Latent for other τ: moment matching could print as proved.
  - X13 (≥), X14 (10⁵, q > 0.739), X15 (K = 120, q > 0.68), X16 (q > 0 dropped): the trigger is pinned only at 0.605 and 0.85.
  - X19 (`mono` with ≤), X20 (`mono`'s last pair unchecked), X22 (`cons` factor 100), X23 (`cons` factor 1), X24 (`cons` against the change two back): the verdict's edges are unpinned. X23 survives because the data's margins are ≤ 0.056.
  - X27 (`flat` < 10⁻²) and X35 (dev/10): unpinned.
  - X36 (B > 0 asserted): an improvement, equivalent on the data.
- **Which change published results.** None of the survivors changes the md. X10 and X13–X16 would change verdicts on rules the raw does not hold.

**9. Literature** (five searches; arxiv.org and www.ams.org fetches blocked) [snippet]
- **Karvonen–Oates–Girolami (Math. Comp. 90, 2021).** The scaled Gauss–Hermite class has bounds exp(−c₁N^{1/d})N^{1/(4d)} and exp(−c₂N^{1/d})N^{−1/(4d)}. The exp(−c₃N^{1/(2d)}) in the abstract is the second, nested class. D58's one-dimensional citation matches.
- **Kjelgaard Mikkelsen–López-Villellas, "How accurate is Richardson's error estimate?"**, and **Stetter (1968).** The estimate is reliable only when the neglected terms are dominated, and Richardson's fraction diagnoses the leading order. This is what D3dr-3's K^{−5/2} mock exhibits.
- **Derevianko–Plonka, arXiv:2307.13587** (v2 May 2025). The search did not confirm D58's venue, Sampl. Theory Signal Process. Data Anal. 24, article 6, 2026.
- **arXiv:2606.03855**, "Approximation by short exponential sums with geometric error decay based on Gauss quadrature". This is probably D58's "Plonka–Riebe–Cuyt 2026"; the authors were not shown.
- **Nothing found** on the Gauss–Hermite error for e^{xz} through a Laplace transform or a semicircle law.

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) :13–17, Status :172–175: the sharp constant "proved by D‴.61 (iv) at every t, every σ*/σ"; earlier tiers as history | **Right.** It follows from D61 (iv) as stated: ρ = t at every σ*/σ, pointwise in t. Corroborated at t = 0.905 to 8·10⁻⁴ at K = 160 | §6 |
| (ii) :26–30, Status :185–188, :193–196, :289–291: t(v) proved for every rule with a frame "at every t(v) < 1"; Fock and Mercer-zero at every τ/σ; frameless conjectural | **Right**, and not overstated. Every frame has t(v) < λ < 1; the thresholds are right; numerically confirmed at t(v) = 0.9945 | §6 |
| :234–235, md a5 line: "t(v) at most 0.4286, below the branch and bound's r* = 0.8052 already" | Right (b(b − 1)/(b² + b + 2) = 3/7 at τ/σ = 2) | §6 |
| :258–264: the transient q = 0.850630, the derivation, "B_K flat to 1.2·10⁻⁴ over K = 60 … 160", B > 0 | **Derivation right; numbers right.** "Over 60 … 160" is two points (60, 80); B > 0 is untested by the code | §2, §5; D3dr-4 |
| :265–267: the power fits worsen; within 6.5·10⁻¹⁶; next order +3.1·10⁻¹⁵; "about 2·10⁻¹⁵ across transient models" | Numbers right. The last figure is round 118's (it needed extra K) and is not reproduced | §2; D3dr-6 |
| :268–272: the orders decrease at the three MM rules and the proved one, each within ten times its last change; "checked rule by rule by the orders" | **Numbers right; the check overstated.** It fails 2 of 3 proved fast-branch rules, a K^{−5/2} misspecification passes it, and the code's verdict is global | §3–5; D3dr-1, -2, -3 |
| :273–276: window geometry; "29 times … not reproduced" | Right. "29" is reproduced by my run (28.9) | §4 |
| Numerics (A6X) :343–346 | Right | §2 |
| md A6X line and verdict | Numbers right; latent defects | D3dr-1–4 |
| (iii)–(v), proofs, B, C, D | Unchanged since a2f3b77 except the tier phrases; tests pass | [md] |

**Status at 783d55b.**
- **Proved, computer-assisted:** t(v) and the sharp asymptotics for every Mehler-frame rule, and the half-squeeze constant at every σ*/σ. This rests on D61's three ball-arithmetic steps, Arb's correctness assumed, and on its adoption pending round 122.
- **Numerical:** A6X's limits, at four frameless fast-branch rules.
- **Conjectural:** ρ = 1 and |ρ| > 1.
- **New observation:** the integer-power premise fails at both proved Mercer-frame rules (stalls at 10⁻¹² and 10⁻¹³) and holds at (2, 2·v_hs).

## Status of the round-118 issues

| D3dm | status | note |
|---|---|---|
| 1 | **fixed in substance**, with residues | The 10× threshold is replaced by a monotone decrease through eight terms plus last-change consistency; "geometrically" is dropped; the Mercer stall is quoted. The criterion is per-rule consistent (§3). Residues: the verdict is global; the known case is not in A6X, and with it (or (1, Mercer)) the md would print a failure for all rules (D3dr-2); its power is limited (D3dr-3) |
| 2 | **partly fixed** | The next order (+3.1·10⁻¹⁵) is computed and printed; D58 says "about 2·10⁻¹⁵". `tr_ok`'s window is still ±6–7·10⁻¹⁵ (D3dr-4); "about 2·10⁻¹⁵" is unmarked (D3dr-6) |
| 3 | **largely fixed** | The derivation is in D58 and correct; B_K replaces the contaminated decays. Residues: K = 120 and 160 are fit nodes; B > 0 is not tested (D3dr-4) |
| 4 | **fixed**, with a residue | `a6x_resolvable`: no misfire at q = 0.605 or 0.675 (§4). New gap: the trigger's K and the orders' K differ (D3dr-1) |
| 5 | **fixed** in wording | "From K⁻² (no K⁻¹ term), checked rule by rule": the check is weaker than the word (D3dr-3) |
| 6 | **partly fixed** | "29" marked; "26-digit" and "at one frameless rule" in D58; the test renamed. Still open: "25-digit" in the md and docstring, "calibrated" in the docstring, the stale `a6x_lines` docstring (D3dr-6) |
| 7 | **partly fixed** | The plants are rewritten; T01–T06 are all killed (reproduced). The transient mock and the threshold edges are still deferred; 15 of my 36 mutants survive (D3dr-7) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dr-1 | minor | **The trigger and the orders look at different K.**<br>• `a6x_resolvable` fires iff q^160 > 10⁻²⁰ (q > 0.7499), but the seven- and eight-term orders use K = 60 and 40.<br>• A component 0.24·q^K with q = 0.55, 0.60, 0.675, 0.74 spoils o₈ by −1.9·10⁻¹⁶, −6.2·10⁻¹⁵, −6.7·10⁻¹³, −2.0·10⁻¹¹, and the verdict prints "do NOT decrease" with the limit exact [own mocks].<br>• Rules in that window, q from the formula: frameless (¼, 1.25·v_hs) 0.581 and (½, 2·v_hs) 0.548 (two of the five unextended frameless fast-branch rules); proved (1, Mercer) 0.675, (2, Mercer) 0.605, (½, Mercer) 0.834 and (¼, Mercer) 0.944 [own] | **`a6x_lines`:** judge resolvability at each order's smallest K (q^{K_min} > 10⁻²⁰ gives q > 0.32 at K = 40), and carry Bq^K in every order's basis where it is resolvable there, or drop from an order the K at which q^K exceeds the input noise. Add a mock with q = 0.6 |
| D3dr-2 | minor | **On known cases the verdict fails 2 of 3, and globally.**<br>• Proved (2, Mercer) at K = 200–320 (my `a6_job`, 108 s): orders +1.2·10⁻¹² … −1.4·10⁻¹³.<br>• Proved (1, Mercer) (413 s): −7.0·10⁻¹¹, +8.5·10⁻¹², −2.2·10⁻¹², +1.2·10⁻¹², −1.7·10⁻¹².<br>• Either added to A6X, the md prints "the orders do NOT decrease … at the 5 [6] rules without a transient", the four passing rules included, beside a hardcoded caveat calling the stall unreproduced.<br>• D58 :271 says "checked rule by rule".<br>• My run also reproduces round 113's "29" (28.9) [own] | **`a6x_lines`:** a per-rule flag ("decreasing: (½, MM), (1, MM), (2, MM), (2, 2·v_hs); not: (2, Mercer) …"). Add (2, Mercer) and (1, Mercer) to A6X_RULES (about 9 min) and print their stalls and the 29 from the run. **D58 :270–272:** "the check fails at 2 of the 3 proved fast-branch rules extended (a non-power component), so a failure says nothing of the limit" |
| D3dr-3 | minor | **The orders' check has limited power and misattributes its failures** [own mocks, 60-digit exact].<br>• A K^{−5/2} term (10⁻³) passes MONO CONS; with a +7.6·10⁻¹² limit offset on top the last order is +2.6·10⁻¹⁴, and the md would read "1/t(v) to the orders' resolution", 290 times too fine.<br>• A K⁻¹ term or log K/K², limit exact, prints "a limit off 1/t(v)".<br>• An exact finite expansion at the 10⁻⁶¹ floor prints "do NOT decrease": there is no noise floor (the propagated input bound is ≤ 1.2·10⁻²³, §2).<br>• `cons` accepts offsets up to 10 times the last change (a 10⁻¹⁷ offset passes with a 4·10⁻¹⁸ change) | **D58 :270–272 and the md's closing clause:** "consistent with an expansion … (a non-integer power can pass this check)". **`a6x_lines`:** the NOT branch "the limit or the premise off"; orders below the propagated input bound count as converged; state the resolution as 10 times the last change, not "\|o₈\|" |
| D3dr-4 | minor | **The transient checks need restating and tightening.**<br>• The transient fit has 5 unknowns through 5 K (120 … 320), so B_120 = B_160 = B exactly. "Flat to 1.2·10⁻⁴ over K = 60 … 160" is K = 60 (−1.18·10⁻⁴) and 80 (−2.5·10⁻⁵).<br>• B > 0 is claimed (D58 :263) but not tested: a mock with B < 0 passes `tr_ok`, and X36 survives.<br>• `tr_ok` accepts outlier offsets of −7.2·10⁻¹⁵ … +5.8·10⁻¹⁵ against the stated resolution of about 2·10⁻¹⁵ and next order 3.1·10⁻¹⁵ (D3dm-2's request).<br>• `flat < 10⁻³` admits γ ≤ 10⁻³ but rejects a (1 + b₁/K) factor at \|b₁\| ≳ 0.1, which the derivation allows [own] | **D58 :263 and the md:** "B_K at K = 60, 80, outside the fit (at 120 and 160 it is B by construction), within 1.2·10⁻⁴". **`tr_ok`:** add B > 0, and bound \|dev\| by \|dev5\| (or 10⁻¹⁴·2) rather than \|o₄\|/100. Add the transient mock (the outlier's t(v), q and smooth power part; γ = 0 passes, γ = 0.01 and B < 0 fail) |
| D3dr-5 | cosmetic | **Dependencies and attributions.**<br>• The tier upgrade rests on D61's adoption of round 121's fix (04d9ead); HANDOFF §5 lists round 122 on that adoption as in flight. My rerun of `part_VA` reproduces the raw's 'VA' exactly, and its counts are round 121's patched ones, so the risk is low.<br>• D58 omits round 120, the certificate's author, and "Arb's correctness assumed", which gh_peano.md states [own, md] | **D58 (i), Status, :287–289:** "VA (referee rounds 120, 121; its adoption under review in round 122), Arb's correctness assumed" |
| D3dr-6 | cosmetic | **Residues.**<br>• The md line and the `a6x_lines` docstring say "25-digit errors"; D58 says "26-digit".<br>• The `a6x_lines` docstring describes n = 3, 4, 5 and the trigger "q > ½".<br>• The md heading reads "(A6X, referee rounds 101, 107, 113)".<br>• The test docstring keeps "calibrated on the raw … is NOT calibrated".<br>• D58 rounds one quantity two ways: 6.5 and 6.6·10⁻¹⁶ (:266, :345), 3.1 and 3.2·10⁻¹⁵.<br>• "About 2·10⁻¹⁵ across transient models (referee round 118)" is not marked as unreproduced.<br>• `d61_step` returns VA for t ≥ 1, and the a5 line rounds t(v) up, so t ≥ 0.99995 would print "1.0000, below VA's r* = 1" | Fix the strings; "(%.1e; rounded up %s)" once; mark "about 2·10⁻¹⁵" "not reproduced by the repository's scripts"; assert t < 1 in `d61_step` |
| D3dr-7 | cosmetic | **Test coverage.** The repository's 10 mutants are all killed. 15 of my 36 survive, all with the md unchanged: X01, X07 (tier boundary and rounding); X10 (the frame tolerance, equivalent at dyadic τ, latent otherwise); X13–X16 (the trigger's threshold, pinned only at 0.605 and 0.85); X19, X20, X22–X24 (`mono`'s and `cons`'s edges, the last pair, the factor); X27, X35 (`flat`, `tr_ok`'s dev); X36 (equivalent). Several md-changing kills were caught first by the md regression; pass 2 confirms unit-test kills for X03, X05, X06, X09, X11, X17 and X29 [own] | Mocks at the edges: an offset just past the last pair (X20); a factor-1 vs 10 case (X22, X23); q = 0.70 and 0.80 against the trigger (X13–X16); flat = 2·10⁻³ (X27); a non-dyadic τ/σ moment-matching rule (X10) |

## Sources

- [Karvonen, Oates, Girolami, Integration in RKHS of Gaussian kernels, arXiv:2004.12654](https://arxiv.org/pdf/2004.12654) [snippet]; [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/) [snippet; host not reachable]; [Newcastle eprint](https://eprints.ncl.ac.uk/276230) [snippet].
- [Karvonen–Särkkä, Gaussian kernel quadrature at scaled Gauss–Hermite nodes, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532) [snippet].
- [Kjelgaard Mikkelsen, López-Villellas, How accurate is Richardson's error estimate?](https://zaguan.unizar.es/record/168390) [snippet]; [The need for accuracy and smoothness in numerical simulations, arXiv:2406.08257](https://arxiv.org/pdf/2406.08257) [snippet]; [Stetter, Apl. Mat. 13 (1968)](https://dml.cz/bitstream/handle/10338.dmlcz/103154/AplMat_13-1968-2_13.pdf) [snippet].
- [Derevianko–Plonka, Differential approximation of the Gaussian by short cosine sums, arXiv:2307.13587](https://arxiv.org/pdf/2307.13587) [snippet]; [Approximation by short exponential sums with geometric error decay based on Gauss quadrature, arXiv:2606.03855](https://arxiv.org/pdf/2606.03855) [snippet].
- [Convergence analysis of Hermite approximations for analytic functions, arXiv:2312.07940](https://arxiv.org/pdf/2312.07940) [snippet].
- Repository [md, raw]: D58.md, D61.md, exp_width_mismatch.py, exp_gh_peano.py, gh_lastann.py, width_mismatch.md and its raw, gh_peano.md and its raw, tests/test_width_mismatch.py, tools/mutants/round118_D58_own.py and round121_D58_tiers.py, REVIEW_round118_D3dm.md, REVIEW_LOG Rounds 118, 120 and 121, HANDOFF §5; `git diff a2f3b77 783d55b`.
- My scripts, all in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref123/own/`:
  - recompute.py, windows.py, scan.py, known.py, known1.py, mocks.py, mocks2.py, mocks3.py, near1.py, rerun_va.py, extra_job.py, spec123.py, mutants123.py;
  - outputs x_mercer1.json and x_mercer2.json.
- Logs: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref123/logs/`. Mutant results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref123/mut/res_*.jsonl`.
- [memory]: an extrapolated limit's error is set by the first neglected term; a basis missing a component stalls or biases successive orders without necessarily breaking their monotonicity.

## Lesson

**Calibrate a convergence diagnostic on every known case you can afford, not on one.** Here the one proved fast-branch rule in A6X passes. The two other proved fast-branch rules, about nine minutes of computing, both fail it. So a failure of "the orders decrease" says nothing about the limit, and a pass says only that this premise check did not object.

**Then try models built to fool it.** A K^{−5/2} term with an offset on top passes, and the printed resolution comes out 290 times too fine. A diagnostic's window, its trigger and its fits must also look at the same K: a transient that is "unresolvable at K = 160" still dominates an order that reaches K = 40.

**When an upstream result is upgraded, check the claim's form, not just its tier.** D‴.61 (iv) is pointwise in r on compact subsets of the disc, and that is exactly what D‴.58 needs. The upgrade was right, and saying "pointwise, not uniform as t → 1" costs nothing.
