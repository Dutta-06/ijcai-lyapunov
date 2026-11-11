# Referee report D3dm, round 118: Proposition D‴.58 as revised after round 113

## Scope

- **Subject.** Commit a2f3b77 on claude/vigilant-clarke-hf939b, diffed against 5be79c1. In that range, round 113's changes are in `a6x_lines` of `experiments/exp_width_mismatch.py`, the A6X paragraph and Numerics of `review/drafts/D58.md`, one line of `results/chaos/width_mismatch.md`, two tests, and the new `tools/mutants/round113_D58_own.py`.
- **Work directory.** W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref118`.
  - Export: `git archive a2f3b77 | tar -x -C W/repo`.
  - WM_OUT = W/out and GP_OUT = W/gpout, each holding a copy of its raw.
  - My scripts are in W/own, logs in W/logs, mutant slots in W/mut.
- **Read in full.**
  - D58.md, and its diff from 5be79c1.
  - Every changed hunk of `exp_width_mismatch.py`. `a6x_lines` line by line: `solve`, `pw`/`one`/`efit`, `ratio_hp`, `geom`, t(v) from `rate_conj` on `a6_scalings`' 60-digit Decimal variance, `hs`/`ho`/q, the transient fit, `res`/`dec`, `conv`, `tr_ok`, `nlab`, `nfast`, and the printed verdict.
  - width_mismatch.md.
  - Both round-107 tests as rewritten; `tools/mutants/round113_D58_own.py`; REVIEW_round113_D3dh.md; the Round 113 entry of REVIEW_LOG.md; HANDOFF §5.
- **Run settings.**
  - Every run used the repository's venv with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. pytest always ran with `--basetemp` inside W.
  - At most two heavy processes of mine ran at a time. No foreground sleep loops.
  - All my processes (PIDs 27898, 27899, 29239, 30166, 30167, 30673, 31430) have exited; I stopped none.
  - My mutant runner `W/own/mutants118.py` is `tools/mutants.py` with three lines changed: REPO = W/repo, PY = the venv's python, and a per-slot `--basetemp`.
- **The repository was not modified.** I ran only `git archive/log/diff/status/rev-parse`. HEAD was still a2f3b77 at the end, with a clean status.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my code, or the repository's functions on inputs I changed.
  - [snippet]: web-search results only. Two hosts did not resolve: asidi.cswp.cs.technion.ac.il and ar5iv.arxiv.org must be allowed in the session's network settings to read the papers.
  - [memory]: recalled facts.

## Summary

**Design.**
- Regenerated the md from the raw and ran the 63 tests.
- Recomputed every new A6X quantity independently: in Arb ball arithmetic (python-flint 0.9.0) and in my own exact-rational solver, with t(v) derived from scratch.
- Bounded how the 26-digit inputs propagate into each fit.
- Bisected the verdict's acceptance windows through the repository's own `a6x_lines`.
- Derived the transient's form and tested it four ways: q free, the prefactor exponent γ, its sign, and its small-K corrections.
- Ran the repository's `a6_job` at new K:
  - the outlier (¼, 2·v_hs) at K = 100, 140, 180 (260 s), for an out-of-sample test of the transient model;
  - moment matching at τ/σ = 2 at K = 400, 480 (559 s), to test the six-term extrapolation beyond K = 320;
  - the proved (2, Mercer frame) at K = 200, 240, 320 (99 s), and its errors at every K = 99…131.
- Tested the no-K⁻¹ assumption, and computed q for all 36 A6 rules.
- Mutants: round 113's spec (U00–U12) plus 33 of my own (V00–V30, W01, W02).
- Five web searches; two fetches, both blocked.

**Verdicts.**
- **Reproducibility is exact.**
  - The md is byte-identical (md5 8b3ba8c5…). The raw (f1079628…) is unchanged since ce7cfac.
  - 63/63 tests pass in 52.4 s.
  - My `a6_error` at K = 119 and 120 reproduces the raw's errors to all 26 digits.
- **Every new number in D58 and the md is right** to the digits printed [own]. That includes −6.5·10⁻¹³, +9.3·10⁻¹¹, −2.9·10⁻⁹; q = 0.850630; the decays 0.85219, 0.85069, 0.85062, 0.85026; +6.5·10⁻¹⁶ (6.6·10⁻¹⁶ rounded up), B 0.2383; the orders at the four rules without a transient (the last within 2.0·10⁻¹⁷); 21.11, 10, 2.533, 2; and "4 of the 9".
  - The fits are exact well past the 26-digit inputs. The linear error propagation is at most Σ|w|·8.3·10⁻²⁶ ≤ 7.4·10⁻²⁴.
- **Round 107's "calibrated by the proved rule" is withdrawn** in D58, the md and the tests.
- **The transient's form Bq^K is right**, and better supported than D58 says (D3dm-3).
- **The new pass/fail criteria are the weak point.**
  - "Each at least 10 times smaller" was set with the data's 10.95 in hand. It flips under a 1.7·10⁻¹⁸ offset, and it fails at a proved rule (D3dm-1).
  - "Within 6.6·10⁻¹⁶" is one model's value (D3dm-2).
  - The transient trigger misfires on a proved rule (D3dm-4).
  - The integer-power premise is not generic: it fails visibly at a proved rule (D3dm-5).
- **No major issue.**

**Findings.**
1. **D3dm-1 (minor).** The convergence criterion "successive orders each at least 10 times smaller" is post hoc, resolves finer than the evidence, and rejects a proved limit.
2. **D3dm-2 (minor).** The outlier's "within 6.6·10⁻¹⁶ with its transient" is the minimal model's value. Its next order moves it to 3.1·10⁻¹⁵, and the verdict's own window is ±6·10⁻¹⁵. The data support ≈ 2·10⁻¹⁵.
3. **D3dm-3 (minor).** The transient's form is right, and I derive it here. D58 states it without derivation. The code tests the decay against a reference that absorbs part of the transient, and tests neither the prefactor's power of K nor the omitted q^{2K} term.
4. **D3dm-4 (minor).** The trigger q > ½ misfires where 0.5 < q ≲ 0.7. Added to A6X, the proved (2, Mercer frame) prints "the transient fit does NOT account for the residual", with B = −1.1·10¹⁵.
5. **D3dm-5 (minor).** "Under an expansion in integer powers of 1/K": the bases start at K⁻², and the premise is rule-specific. At the proved (2, Mercer frame) the orders stall at ≈ 10⁻¹³ through eight terms.
6. **D3dm-6 (cosmetic).** Residues: the hardcoded, unreproduced "29 times"; a stale test name and docstring; "converge geometrically"; "25-digit"; "beyond the frames" rests on one rule.
7. **D3dm-7 (cosmetic).** Tests: 17 of my 30 non-control mutants survive, all equivalent on the data. Three are caught only by the md regression. Round 113's Z06 and Z08 are still unconstrained.

## What I checked and how

**1. Reproduction** [md, raw, own]
- `--from-raw` into W/out (1.6 s), then `cmp`: identical.
- Full test file: 63 passed in 52.4 s.
- `W/own/consec.py` computes e_K at every K = 99…131 for (2, Mercer frame). Its e₁₁₉ and e₁₂₀ match the raw's strings to all 26 significant digits.

**2. A6X recomputed independently** (`W/own/recompute.py` in Arb at 600 bits; `W/own/variants.py` in exact rationals) [own, raw]

| rule | t(v) float − exact | four terms (160–320) | five terms (120–320) | six terms (80–320) | ratios |
|---|---|---|---|---|---|
| (¼, 2·v_hs) | −6.75·10⁻¹⁸ | −6.5193·10⁻¹³ | +9.3071·10⁻¹¹ | −2.8792·10⁻⁹ | worsening |
| (½, MM) | −6.17·10⁻¹⁸ | −3.861·10⁻¹⁴ | −2.092·10⁻¹⁶ | +1.91·10⁻¹⁷ | 184.6, 10.95 |
| (1, MM) | −1.85·10⁻¹⁷ | +8.612·10⁻¹⁵ | −3.267·10⁻¹⁶ | +2.08·10⁻¹⁸ | 26.4, 157 |
| (2, MM) | −3.70·10⁻¹⁷ | −3.568·10⁻¹⁴ | +5.185·10⁻¹⁶ | −3.64·10⁻¹⁸ | 68.8, 142 |
| (2, 2·v_hs), proved | +2.39·10⁻¹⁷ | −1.985·10⁻¹⁴ | +4.947·10⁻¹⁶ | −7.43·10⁻¹⁸ | 40.1, 66.6 |

- **The outlier's transient fit.** +6.5243·10⁻¹⁶, B = 0.23830161. The limit's weights are 0.0022, −1.416, 9.134, −14.16, 7.442 on K = 120…320: the fit is essentially a four-term fit through 160–320, corrected by the transient that B fitted at K = 120.
- **The outlier's decays.** 0.852194, 0.850693, 0.850622, 0.850262 (dec/q − 1 = +1.84·10⁻³, +7.5·10⁻⁵, −8.3·10⁻⁶, −4.32·10⁻⁴). q = 0.8506295388.
- **The four-term deviation is the transient's.** Predicted −6.526·10⁻¹³ from the transient against the fitted −6.519·10⁻¹³.
- **Geometry:** 21.1111, 10, 2.53333, 2.
- **Input precision.** Each stored error has 26 significant digits ('.25e'), so each ratio carries a relative error of at most 8.3·10⁻²⁶.
  - The limit weights' Σ|w| is 31.96, 60.69 and 88.75 for four, five and six terms, so the linear bound is ≤ 7.4·10⁻²⁴.
  - Arb's interval-LU enclosure is rigorous but loose: 2·10⁻²² … 1.1·10⁻¹⁸.
  - q rounded to a double moves the transient fit by about 3·10⁻²⁷.
  - So the exact-rational fits are exact far below every published value, and "t(v) exact" (60 digits) is more than enough. The float t(v)'s error of 6·10⁻¹⁸ … 4·10⁻¹⁷ was indeed at the six-term level.

**3. The verdict's acceptance windows** (`W/own/windows.py`; the repository's `a6x_lines`, e₀ and the ratio planted by (1 + d), as the test's `pl` does) [own]

| rule | `conv` holds for d in | `tr_ok` holds for d in |
|---|---|---|
| (½, MM) | (−4.447·10⁻¹⁷, **+1.658·10⁻¹⁸**) | — |
| (1, MM) | (−3.862·10⁻¹⁷, +2.784·10⁻¹⁷) | — |
| (2, MM) | (−4.387·10⁻¹⁷, +6.175·10⁻¹⁷) | — |
| (2, 2·v_hs) | (−3.827·10⁻¹⁷, +6.330·10⁻¹⁷) | — |
| (¼, 2·v_hs) | — | (−7.248·10⁻¹⁵, +5.814·10⁻¹⁵) |

- At (½, MM), a relative offset of +1.7·10⁻¹⁸ prints "successive orders NOT each 10 times smaller". That is 1/12 of the published "within 2.0·10⁻¹⁷", and below the six-term fit's own bias (its seven- and eight-term successors are −4.8·10⁻¹⁹ and −7.4·10⁻²⁰).
- `tr_ok` accepts limits about 10 times farther off than the "within 6.6·10⁻¹⁶" it prints.
- `tr_ok`'s margins on the data:
  - |dev| < |o₄|/100: the ratio on the data is 999.
  - The decays within 2·10⁻³ of q: 4.6× margin. The excluded first decay (1.84·10⁻³) would pass too.
  - The worsening condition: 143× and 31×.

**4. Successive orders beyond six terms** (own; through the last n + 1 of the nine K)

| rule | 4 … 8 terms | ratios |
|---|---|---|
| (½, MM) | −3.86·10⁻¹⁴, −2.09·10⁻¹⁶, +1.91·10⁻¹⁷, −4.83·10⁻¹⁹, −7.41·10⁻²⁰ | 185, 11.0, 39.5, 6.5 |
| (1, MM) | +8.61·10⁻¹⁵, −3.27·10⁻¹⁶, +2.08·10⁻¹⁸, +7.74·10⁻¹⁹, −1.18·10⁻¹⁹ | 26, 157, **2.7**, 6.6 |
| (2, MM) | −3.57·10⁻¹⁴, +5.18·10⁻¹⁶, −3.64·10⁻¹⁸, −2.91·10⁻²⁰, −1.05·10⁻²⁰ | 69, 142, 125, **2.8** |
| (2, 2·v_hs) | −1.98·10⁻¹⁴, +4.95·10⁻¹⁶, −7.43·10⁻¹⁸, +1.88·10⁻¹⁹, −3.16·10⁻²⁰ | 40, 67, 40, 5.9 |
| **(2, Mercer frame), proved** (my K = 200–320 run) | +1.20·10⁻¹², +2.52·10⁻¹³, −1.25·10⁻¹³, +7.87·10⁻¹⁴, −1.38·10⁻¹³ | **4.75, 2.01, 1.60, 0.57** |

- At the four published rules the magnitudes decrease monotonically through eight terms. That is a sound statement; "each ≥ 10×" holds only for the three orders chosen.
- At the proved Mercer-frame rule (limit 7/3, D‴.61 and A5) the orders stall at about 10⁻¹³ with alternating signs. Round 113's report had already tabulated its first three.
- **What the Mercer frame's stall is** (`W/own/consec_an*.py`, `mercer_rate.py`):
  - excess × K² is smooth over K = 100…131, with no (−1)^K alternation;
  - with the limit fixed at 7/3, the residual against K⁻²…K⁻⁵ falls at about 0.90 per kernel;
  - a single fitted rate (0.912) misses K = 80 and 60 by 10⁻¹¹ and 4·10⁻¹⁰;
  - no fixed q′ makes the free-limit orders converge.
  - It is not (h_other/t(v))² = 0.605, nor its root 0.778. I leave it unidentified [own].

**5. Extrapolation out of sample, moment matching at τ/σ = 2** (my `a6_job` at K = 400, 480; numpy's weight-overflow warning is harmless, since only its nodes seed Newton) [own]
- The ratios are 1.499999560456823 and 1.4999996949078018.
- Fitted on K ≤ 320, the five-, six- and seven-term fits predict R₄₀₀ to −3.1·10⁻¹⁷, +1.9·10⁻¹⁹ and +1.3·10⁻²¹ (relative).
- With K = 400 and 480 included, the four- … eight-term orders are −3.75·10⁻¹⁵, +3.27·10⁻¹⁷, −1.14·10⁻¹⁹, −5.88·10⁻²², −6.49·10⁻²³, down to the input noise.
- So the moment-matching limit λ is numerically secure there to about 10⁻²², far beyond the published claim. The evidence is strong; the pass/fail criterion is what is fragile.

**6. The transient: its form, derived** [own, memory]
- **Setting.** Suppose e_K² = a_K² + b_K² with
  - a_K = α t^K K^p (1 + O(1/K)), t = t(v);
  - b_K = β |h|^K K^{p′} (1 + O(1/K)), |h| < t, h = h_other.
- **The ratio.** Put ε_K = (b_K/a_K)² = (β/α)² K^γ q^K (1 + O(1/K)), with q = (h/t)² and γ = 2(p′ − p). Then

  R_K = e_{K−1}/e_K = (a_{K−1}/a_K)·[(1 + ε_{K−1})/(1 + ε_K)]^{1/2} = R_K^{(a)} + B K^γ q^K (1 + b₁/K + …) + O(q^{2K}),

  with B = (β/α)²(1 − q)/(2tq) > 0 and R_K^{(a)} = t⁻¹(1 − p/K + c₂/K² + …).
- **Consequences.**
  - Bq^K is the leading form iff γ = 0, i.e. both axes carry the same power of K. For the Laplace sketch's two axis maxima I expect both to carry K⁰ off the cusp.
  - B > 0 iff the contributions add in e_K². An interference term 2a_Kb_K would decay at √q = 0.922 instead.
  - The omitted terms are a relative O(1/K) and an absolute O(q^{2K}), whose relative size is ε_K = 2tqBq^K/(1 − q) = 0.94q^K.
  - The power part has no K⁻¹ term iff p = 0.
  - Sidi's GREP^(m) for A_n ∼ A + Σ_k ζ_k^n × (a power series in n) is this form [snippet].
- **The data** (`variants.py`, `prefactor.py`, with my K = 100, 140, 180):
  - **q free:** 0.8506301 fitted from K = 80 and 0.8506312 from K = 60, against the predicted 0.8506295.
  - **γ:** r_K/q^K = 0.238274, 0.238296, 0.238300, 0.238302, 0.238298, 0.238302 at K = 60…160, so the local γ is within ±1·10⁻⁴ over K = 80–160.
  - **Sign:** B = +0.2383.
  - **The q^{2K} term:** r_K/q^K is 3.8% low at K = 20 and 0.16% low at K = 40, against ε_K = 3.7% and 0.15%.
  - **b₁:** fitted at −0.006.
- **So the two-axis reading is well supported at this rule:** "as if the two arguments' contributions to e_K² added".
- **The decay 0.85026 (80→120)** is 4.3·10⁻⁴ off q only because the four-term reference fit passes through K = 160, where the transient is still 1.4·10⁻¹².

**7. The outlier's limit under alternative models** (exact rationals; my K = 100, 140, 180 out of sample) [own]
- **Transient plus n power terms, on the stored data:** +6.52·10⁻¹⁶ (n = 3, K 120–320), +3.11·10⁻¹⁵ (n = 4, 80–320), −6.01·10⁻¹³ (n = 5, 60–320), −7.82·10⁻¹² (n = 6, 40–320). The same diagnostic that indicts the power fits indicts the transient model below K ≈ 100, through its unmodelled q^{2K} term.
- **Out of sample at K = 100, 140, 180**, (R − prediction)·t:

  | model | K = 100 | K = 140 | K = 180 |
  |---|---|---|---|
  | four power terms | +7.8·10⁻⁹ | +1.0·10⁻¹¹ | −8.2·10⁻¹⁴ |
  | minimal transient | −3.9·10⁻¹⁴ | −1.8·10⁻¹⁶ | +1.4·10⁻¹⁷ |
  | transient + five powers (80–320) | +1.2·10⁻¹³ | −9.6·10⁻¹⁶ | +6.0·10⁻¹⁷ |

  The minimal transient model is validated.
- **On the clean set K ≥ 120** (`clean_orders.py`, which adds my K = 140 and 180):
  - four powers + Bq^K over 20 windows: +5.6·10⁻¹⁶ … +2.0·10⁻¹⁵;
  - five powers + Bq^K over 7 windows: −8.8·10⁻¹⁷ … +8.1·10⁻¹⁶;
  - six powers + Bq^K: −6.3·10⁻¹⁶.
- **The limit's sensitivity to b₁** is 1.35·10⁻¹⁵ per unit.
- **Verdict:** the limit is 1/t(v) to about 2·10⁻¹⁵, 300 times better than the four-term power fit. "Within 6.6·10⁻¹⁶" is one model's value.

**8. Other assumptions** (`kterm_slow.py`) [own, raw]
- **The K⁻¹ term.** Adding it to the bases gives limits of +2.9·10⁻¹³, +3.0·10⁻¹⁵, −3.7·10⁻¹⁶ at (½, MM). Its fitted coefficient over L is −3.1·10⁻¹⁰, −3.3·10⁻¹², +4.3·10⁻¹³, shrinking with the order. So excluding it is consistent with the data, but it is an assumption (p = 0).
- **q over A6.**
  - Only the outlier and its slow-branch mirror (¼, 0.5·v_hs) have q = 0.850630.
  - At the mirror, the transient changes A6's limit only from 9.3·10⁻⁷ to 7.9·10⁻⁷.
  - The slow branch's slow excess decay (exponent 1.52 at (2, 0.5·v_hs), where q = 0.0077) is not a transient.
  - Among the frameless fast-branch rules not yet extended, (½, 2·v_hs) has q = 0.548 and (¼, 1.25·v_hs) has q = 0.581. The proved (2, Mercer frame) has q = 0.605.
- **The frame census.** The frame condition is |ln(v/v_hs)| < r. At τ/σ = ¼, r = 0.0303 < A6_NEAR, so "no far fast-branch rule has a frame there" holds exactly.

**9. The latent trigger** (`latent.py`: the repository's `a6x_lines` on the raw's A6X plus my (2, Mercer frame) run) [own]
- The run itself matches round 113's last-window values: −9.69·10⁻¹¹ and −1.62·10⁻⁸.
- The output for that rule:
  - the transient path triggers (q = 0.604938);
  - B = −1108942274738012;
  - the transient fit equals the four-term fit, +1.2·10⁻¹²;
  - the decays read 0.70455, 0.92863, 0.92073.
- The verdict then prints "at the transient's rules the transient fit does NOT account for the residual", for a proved limit.
- A q^K term is resolvable at K ≥ 160 in 26-digit data only if q^{160} ≳ 10⁻²⁵, i.e. q ≳ 0.70.

**10. Tests, line by line** [md, own]
- **`test_round107_d61_status_pinned_to_its_raw_and_the_a6x_verdict_calibrates`.**
  - The D‴.61 pins are unchanged; V's contiguity with DB is now asserted ✓.
  - `report()` is asserted on the full verdict string, the geometry, q and the decays, and "calibrated by it" absent ✓.
  - Plants go through e₀ as well as the ratio:
    - 5·10⁻⁸ at the outlier gives "does NOT account" ✓;
    - 10⁻¹⁵ at (½, MM) gives "NOT each 10" ✓;
    - the unplanted data give no "NOT" ✓;
    - the K = 40 plant (U07) gives "does NOT account" while the orders still pass ✓.
  - The name and the docstring still describe the withdrawn calibration.
- **`test_round107_a6x_verdict_plants`.**
  - Its mock ratio 4(1 + 1/(K + 1)²) is exact to 40 digits ✓; conv passes ✓; "1 of the 1" ✓; a 10⁻¹² offset gives NOT ✓; the t(v) = 0.995 frame rule gives "the 2 unproved rules" and "2 of the 2" ✓.
  - The docstring's "25-digit errors" and "every order of 1/K present" are loose: the mock has 40 digits and no K⁻¹ term.
  - No mock carries a transient, so `tr_ok` is pinned only by the raw.

**11. Mutants** (`W/own/spec118.py`, combined with round 113's spec in `W/own/spec_all.py`) [own]
- **Method.**
  - Pass 1: `-k "round101 or round107 or round88 or a6"`, with the md compared.
  - The kills were then classified without `round88_precision_and_the_md`.
  - The survivors were run against all 63 tests. The controls U00 and V00 pass.
- **Round 113's spec.** U01–U04 and U07–U12 are killed by unit tests (10). U05 and U06 survive the full file, equivalent on the data. This reproduces the log.
- **Mine.**
  - Killed by unit tests (12): V01, V04, V06, V09, V14 (the second factor at 11, since 10.95 < 11), V16 (/1000, since 999 < 1000), V19, V24, V26, V27, V29, V30.
  - Changing the md but caught only by the md regression (3): V05 (the printed K₀), V07 (B rescaled by q), V08 (B from the e/K⁴ coefficient).
  - Survivors, all 63 tests passing and the md unchanged (17):
    - V02 (t(v) at 30 digits);
    - V03 (the other h always hs[0]; every extended rule is on the fast branch);
    - V10, V11, V28 (equivalent by construction);
    - V12 (`bool(nofit)`);
    - V13 (the first conv condition dropped);
    - V15 (factor 9), V17 (/10), V18 (5·10⁻⁴);
    - V20 (the decay over the first window only), V21 (over all windows);
    - V22 (absolute tolerance), V23 (worsening on the first pair only);
    - V25 (min over one transient rule);
    - W01 = round 113's Z06 (proved at D61_RSTAR);
    - W02 = Z08 (the necessary condition on fit 1 only).
  - V20 survives because the K = 40 plant also moves the 20→40 decay, to 0.85541.
- **Which change published results.** None of the survivors does. V05, V07 and V08 change the md under passing unit tests. V03, W01, W02, V25 and the threshold survivors would change verdicts on data the tests do not plant: a slow-branch transient rule, a frame rule with 0.8052 ≤ t(v) < 0.99, a fit-2-only growth, or several transient rules.

**12. Literature** (five searches; two fetches blocked) [snippet]
- **Sidi, GREP^(m) (Math. Comp. 2010).** Convergence and stability for A_n ∼ A + Σ_k ζ_k^n × (a power series in n). This is the transient-plus-power form derived in §6, including the (1 + b₁/K + …) factor that the repository's single Bq^K omits.
- **Ford–Sidi's W^(m)-algorithm (SINUM 1987)** and **Sidi's GREP (1979)**: Richardson generalised to known basis functions.
- **Sidi's vectorized GREP (arXiv:1605.02630).**
- **A 2026 extrapolation survey (arXiv:2602.01301)**: Shanks' transformation annihilates pure exponential components; Richardson handles algebraic ones.
- **Laplace-method notes**: subdominant maxima contribute their own terms.
- I could not read the full texts: asidi.cswp.cs.technion.ac.il and ar5iv.arxiv.org must be allowed in the session's network settings.
- D58's literature passages are unchanged since 5be79c1.

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| :256–259 "what remains is a geometric transient of the formula's other branch … 0.85069, 0.85062, 0.85026 … the predicted q = (h(τ²/v)/t(v))² = 0.850630, as if the two arguments' contributions to e_K² added" | **Right, and stronger than stated:** q free to 7·10⁻⁷, γ = 0 ± 10⁻⁴, B > 0, the predicted O(q^{2K}) size at K = 20 and 40. 0.85026's offset of 4.3·10⁻⁴ is the reference fit's contamination | §2, §6; D3dm-3 |
| :259 "evidence for the two-axis picture beyond the frames" | Right at one rule (ρ = 16.5). Its slow-branch mirror and the other frameless rules are untested | §8; D3dm-6 |
| :260–261 "power fits worsen … −6.5·10⁻¹³, +9.3·10⁻¹¹, −2.9·10⁻⁹ … from the 25-digit errors, t(v) exact" | Right (26 significant digits; t(v) to 60 digits) | §2 |
| :261 "with the transient Bq^K (B 0.2383) through K = 120 … 320 its limit lies within 6.6·10⁻¹⁶ of 1/t(v)" | **Number right, the reading over-precise:** about 2·10⁻¹⁵ across validated models | §7; D3dm-2 |
| :263–265 "successive orders converge geometrically (… −3.9·10⁻¹⁴, −2.1·10⁻¹⁶, +1.9·10⁻¹⁷ …; the last within 2.0·10⁻¹⁷ at all four) — numerical evidence … under an expansion in integer powers of 1/K, not a proof" | **Numbers right, the status right.** "Geometrically" is loose (ratios 185 and 11). The expansion starts at K⁻² and is rule-specific; it fails at the proved Mercer frame | §2, §4, §5; D3dm-1, -5 |
| :265–267 "window geometry — … 21.11 and 10, 2.533 and 2, as the proved rule's 20.7, 10, 2.52, 2" | Right | §2 |
| :267–269 "a proved rule's bias calibrates nothing: … differ by up to 29 times …; at τ/σ = 1/4 no far fast-branch rule has a frame …; round 107's 'calibrated by the proved rule' withdrawn" | Right. "29" is round 113's number, hardcoded in the md and not reproduced. τ/σ = ¼ is exact | §8; D3dm-6 |
| Numerics :334–336 (A6X) | Right ("up to +9.03·10⁻¹² and −5.85·10⁻⁹", signed values) | §2 |
| md A6X line and verdict | Numbers right; the 10× criterion is fragile and the transient trigger latent | §3, §9; D3dm-1, -4 |
| (i)–(v), proofs, Status, B, C, D, A6E | Unchanged since 5be79c1; tests pass | [md] |

**What is proved, numerical or conjectural, at a2f3b77.** Unchanged from round 113, with these updates:
- **Numerical, now stronger:**
  - The moment-matching limits λ: at τ/σ = 2, about 10⁻²² with my K = 400 and 480.
  - The outlier's limit, under the derived transient form: about 2·10⁻¹⁵.
- **Numerical:** the transient's rate q = (h_other/t(v))² (to 7·10⁻⁷) and the additivity of the two axes in e_K² at one frameless rule.
- **Conjectural:** the two-axis picture for |ρ| > 1 in general.
- **New observation:** the integer-power premise fails at the proved (2, Mercer frame) at the 10⁻¹³ level through K = 320, so that premise is checked rule by rule, not assumed.

## Status of the round-113 issues

| D3dh | status | note |
|---|---|---|
| 1 | **fixed**, with residues | "Calibrated by it" is withdrawn in D58, the md and the tests; geometry is printed beside a pure remainder's; the τ/σ = ¼ impossibility is stated. Residues: the replacement criterion (10×) repeats the post-hoc-threshold pattern and fails round 113's own relabel-a-proved-rule test (D3dm-1). "29 times" is hardcoded and not reproduced (D3dm-6) |
| 2 | **fixed**, with residues | The transient is fitted at fixed q and printed; the decays, B and the worsening are computed; "(the proved 2.0e-14)" is dropped; D58 is rewritten. Residues: over-precision (D3dm-2); no derivation and weak acceptance tests (D3dm-3); the trigger (D3dm-4) |
| 3 | **fixed** | The verdict reads computed conditions only; plants go through e₀; `nfast` counts the same population ("2 of the 2 unproved"); the cal-digit point is moot. New latent defect in the transient path: D3dm-4 |
| 4 | **largely fixed** | Mock rules with exact errors, contiguity, full strings; U07 is killed by a plant. Residues: Z06 and Z08 (my W01, W02) still survive; V05, V07, V08 are md-only kills; the thresholds are unpinned (D3dm-7) |
| 5 | **fixed** | Signed values in the Numerics; "the 2 unproved rules"; HANDOFF §5 updated; R19 and R21 recorded in the log |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dm-1 | minor | **The 10× convergence criterion is post hoc, resolves finer than the evidence, and rejects a proved limit.**<br>• The threshold is 10 with the data's smallest ratio 10.95 in hand ((½, MM), five → six terms; margin 1.095). Round 107's 3 had the data's 2.70 in hand (margin 1.11).<br>• Through the repository's `a6x_lines`, a relative offset of +1.66·10⁻¹⁸ at (½, MM) prints "NOT each 10 times smaller" (window −4.45·10⁻¹⁷ … +1.66·10⁻¹⁸). That is 1/12 of the printed "within 2.0·10⁻¹⁷", and below the six-term fit's own bias (its seven-term successor is −4.8·10⁻¹⁹) [own].<br>• One order further, the criterion fails at (1, MM) (2.7×); two orders further, at all four (6.5, 6.6, 2.8, 5.9×).<br>• At the proved (2, Mercer frame), limit 7/3, the orders are +1.20·10⁻¹², +2.52·10⁻¹³, −1.25·10⁻¹³, +7.87·10⁻¹⁴, −1.38·10⁻¹³ (ratios 4.75, 2.01, …), so the criterion rejects a true limit. Round 113 had tabulated this rule's first three orders [own run, raw].<br>• The evidence itself is strong: at (2, MM), K = 400 and 480 are predicted to 2·10⁻¹⁹, and the eight-term limit is within 6.5·10⁻²³ [own]. | **`a6x_lines` and D58 :263–265, Numerics :335–336:**<br>• Replace the pass/fail by what is robust: print the orders through eight terms and state "decreasing in magnitude through eight terms at the four rules (the last ≤ 1.2·10⁻¹⁹)"; state the limit's agreement at the last order's size; drop "geometrically".<br>• Say that the diagnostic is rule-specific: at the proved (2, Mercer frame) it stalls at about 10⁻¹³ (add that rule to `A6X_RULES`, 99 s).<br>• If a criterion is kept, define it relative to the next order's change, not to 1/t(v), and test it on that proved rule |
| D3dm-2 | minor | **The outlier's "within 6.6·10⁻¹⁶" is one model's value.**<br>• The next transient order (five powers through K = 80) gives +3.11·10⁻¹⁵, then −6.0·10⁻¹³ (from 60) and −7.8·10⁻¹² (from 40). The cause is the transient's own O(q^{2K}) term, outside the basis.<br>• The verdict accepts offsets of −7.2·10⁻¹⁵ … +5.8·10⁻¹⁵, ten times the printed bound.<br>• With my K = 140 and 180 (the repository's `a6_job`), every transient model on K ≥ 120 lies in −6.3·10⁻¹⁶ … +2.0·10⁻¹⁵. The minimal model predicts R at K = 140 and 180 to 1.8·10⁻¹⁶ and 1.4·10⁻¹⁷ [own]. | **D58 :261, Numerics :336 and the md:** "within 6.5·10⁻¹⁶ by the minimal transient fit; about 2·10⁻¹⁵ across transient models (its next order through K = 80 is 3.1·10⁻¹⁵, the q^{2K} term unmodelled)".<br>Optionally add K = 140 and 180 to A6X for the outlier (≈ 4 min) and print the spread over the K ≥ 120 windows. Make `tr_ok`'s dev condition match the stated resolution |
| D3dm-3 | minor | **The transient's form is right but underived, and the code tests the wrong things.**<br>• Derivation (§6): if e_K² = a_K² + b_K², then R_K = R^{(a)} + BK^γq^K(1 + b₁/K + …) + O(q^{2K}), with B > 0 and γ = 2(p′ − p).<br>• The data give q free 0.8506301 (predicted 0.8506295); r_K/q^K flat at 0.238298–0.238302 over K = 60–160 (γ within ±10⁻⁴); B > 0, so additive, with no √q interference; 3.8% and 0.16% low at K = 20 and 40, against ε_K = 3.7% and 0.15%.<br>• `tr_ok` checks the decays against the four-term fit, which absorbs the transient at K = 160 (hence 0.85026, 4.3·10⁻⁴ off). It tests neither γ nor the q^{2K} correction [own]. | **D58 :257–259:** add the two-line derivation and its four checks; this is stronger evidence for "as if the contributions to e_K² added".<br>**`a6x_lines`:** test r_K/q^K = B against the transient fit's power part over K = 60–160 instead of the contaminated decays. Optionally add the q^{2K} term, so that K ≤ 80 can enter |
| D3dm-4 | minor | **The transient trigger q > ½ misfires.**<br>• At K ≥ 160, q^K is resolvable in 26-digit data only for q ≳ 0.70.<br>• With the proved (2, Mercer frame) (q = 0.605) added to A6X, the repository's `a6x_lines` prints B = −1108942274738012, a transient fit equal to the four-term fit (+1.2·10⁻¹²), and "the transient fit does NOT account for the residual", for a proved limit [own].<br>• The same would hit (½, 2·v_hs) (q = 0.548) and (¼, 1.25·v_hs) (q = 0.581), two of the five unextended frameless fast-branch rules.<br>No published number changes. | **`a6x_lines`:** trigger on the transient's size in the fitting window (B q^{160} above about 10⁻²⁰ relative, i.e. q ≳ 0.75), not on q > ½. Add a mock test with q = 0.6 |
| D3dm-5 | minor | **"Under an expansion in integer powers of 1/K" is right as a status, but stated too generally.**<br>• The bases start at K⁻²: the K⁻¹ term is excluded by assumption. That is consistent with the data: fitted, its coefficient over L is ≲ 4·10⁻¹³ and falling.<br>• The premise is rule-specific. At the proved (2, Mercer frame), a smooth non-power component (about 0.90 per kernel over K = 100–131, no (−1)^K alternation, not one geometric term) stalls the orders at about 10⁻¹³ [own consecutive-K run]. | **D58 :264–265 and the md's closing clause:** "in integer powers of 1/K from K⁻² (no K⁻¹ term), checked per rule by the orders' convergence — which fails at the proved (2, Mercer frame), whose limit rests on D‴.61" |
| D3dm-6 | cosmetic | **Residues.**<br>• The md hardcodes "(by up to 29 times at τ/σ = 2: referee round 113)"; the repository's raw holds one proved fast-branch rule at A6X's K.<br>• The test name `…_a6x_verdict_calibrates` and its docstring describe the withdrawn calibration; the mock docstring says "25-digit … every order of 1/K".<br>• "25-digit errors": they have 26 significant digits ('.25e').<br>• "Evidence … beyond the frames" rests on one rule; its slow-branch mirror (¼, 0.5·v_hs), with the same q, is untested beyond K = 160. | Mark "29" as "not reproduced by the repository's scripts" (as D58 does elsewhere), or add (2, Mercer frame) and (2, 1.25·v_hs) to `A6X_RULES` (≈ 100 s each). Rename the test and fix both docstrings. "26-digit". "At one frameless rule" |
| D3dm-7 | cosmetic | **Test coverage.**<br>• 17 of my 30 non-control mutants survive, all equivalent on the data: V03 (the slow-branch "other h"), V13 and V15/V17/V18/V20–V23 (`conv`'s and `tr_ok`'s thresholds, windows and conditions), V12, V25, V02, V10/V11/V28 (by construction), and W01/W02 (round 113's Z06/Z08, still unconstrained).<br>• V05, V07 and V08 (the printed K₀ and B) are caught only by the md regression.<br>• No mock carries a transient. | A mock rule with an exact transient (R = 4(1 + 1/(K + 1)²) + Bq^K at q = 0.85, on the slow branch for V03) asserting q, B and K₀ in the line. Mocks at each threshold's edge. A frame rule with 0.8052 ≤ t(v) < 0.99 (W01). A fit-2-only growth (W02) |

## Sources

- [Sidi, GREP^(m) for linearly convergent sequences, Math. Comp. (2010), P107](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P107_Math.Comp_.GREP^m.pdf) [snippet; host not reachable]. A_n ∼ A + Σ_k ζ_k^n × (a power series in n): the form derived in §6.
- [Ford–Sidi, the W^(m)-algorithm, SINUM 1987, P034](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P034_SINUM.Wm-alg.pdf) and [Sidi, the W-algorithm, Numer. Math. 1982, P015](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P015_Numer.Math_.W-alg.pdf) [snippet]. Recursive implementations of GREP.
- [Sidi, GREP, JIMA 1979, P006](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P006_JIMA.GREP_.pdf) [snippet]. Richardson extrapolation for A + Σ c_kφ_k(y) with known φ_k.
- [Sidi, On a vectorized version of a generalized Richardson extrapolation process, arXiv:1605.02630](https://arxiv.org/pdf/1605.02630) [snippet].
- [A survey of scalar and vector extrapolation, arXiv:2602.01301](https://arxiv.org/html/2602.01301v1) [snippet]. Shanks' transformation annihilates exponential components; Richardson handles algebraic ones.
- [Estimagic, Richardson extrapolation](https://estimagic.readthedocs.io/en/stable/explanations/differentiation/richardson_extrapolation.html) [snippet].
- [Ivrii, APM textbook, Expansion of Integrals 1](https://www.math.toronto.edu/ivrii/APM-textbook/Chapter2/L2.1.html) [snippet]. Laplace's method with several maxima, each contributing.
- Repository [md, raw]:
  - D58.md, `exp_width_mismatch.py`, width_mismatch.md and its raw, `tests/test_width_mismatch.py`, and `tools/mutants/round113_D58_own.py`, all at a2f3b77;
  - REVIEW_round113_D3dh.md; REVIEW_LOG.md Round 113; HANDOFF §5;
  - `git diff 5be79c1 a2f3b77`.
- My scripts and outputs (all in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref118/own/`):
  - recompute.py, variants.py, windows.py, kterm_slow.py, oos.py, oos_mm.py, prefactor.py, clean_orders.py, latent.py, mercer_orders.py, consec.py, consec_an.py, consec_an2.py, mercer_rate.py, extra_job.py;
  - spec118.py, spec_all.py, mutants118.py;
  - x_outlier.json, x_mm2.json, x_mercer2.json, consec_mercer2.json.

  Logs are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref118/logs/`; mutant results in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref118/mut/res_p1.jsonl`, `res_p2.jsonl` and `res_p3.jsonl`.
- [memory]: an extrapolated limit's error is set by the first neglected term, and successive orders stall when a component lies outside the basis or the series is only asymptotic.

## Lesson

**A convergence test needs the same check a calibration does: run it on a known case.** Here the known case was already on record. Round 113 had tabulated the proved (2, Mercer frame)'s orders, +1.2·10⁻¹², +2.5·10⁻¹³, −1.3·10⁻¹³, and on them "each at least 10 times smaller" fails. A threshold set just below the data's own 10.95, and broken by a 1.7·10⁻¹⁸ offset, is the round-107 threshold again in a new form. Report the orders and their monotone decrease; do not report a pass.

**When a new term in the basis fixes the fit, derive its form before trusting its limit, and apply to the enlarged model the same diagnostic that exposed the gap.** Here the derivation predicts Bq^K with γ = 0, B > 0, and an O(q^{2K}) next term, and all three are confirmed. The same derivation explains why the transient model, too, worsens when K = 80 enters. Its limit is good to about 2·10⁻¹⁵, not 6.6·10⁻¹⁶.
