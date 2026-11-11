# Referee report D3dx, round 129: Proposition D‴.58 as revised after round 123

## Scope

- **Subject.** Commit d3db5ac on claude/vigilant-clarke-hf939b (round 123 applied to D‴.58), diffed against 783d55b.
- **Work directory.** W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129`.
  - Export: `git archive d3db5ac | tar -x -C W/repo`, with W/repo/venv a symlink to the repository's venv.
  - WM_OUT = W/out and GP_OUT = W/gpout, each holding a copy of its raw.
  - pytest ran with `--basetemp` inside W.
  - My scripts are in W/own, logs in W/logs, mutant slots and results in W/mut.
- **Read in full.**
  - review/drafts/D58.md.
  - In exp_width_mismatch.py: every hunk changed since 783d55b, and what those hunks call — `a6x_lines`, `a6x_decreasing`, `a6x_resolution`, `a6x_transient_ok`, `a6x_resolvable`, `part_A6X` (with `--rules`, and `main`'s parsing of it) and `d61_step`. Also `a6_job`, `a6_error`, `a6_scalings`, `rate_conj` and `gh_nodes`.
  - results/chaos/width_mismatch.md (its A6X line in full).
  - tests/test_width_mismatch.py (the diff, and the three new round-123 tests with `_a6x_mock`).
  - tools/mutants/round123_D58_own.py and tools/mutants.py.
  - review/REVIEW_round123_D3dr.md and the Round 123 entry of REVIEW_LOG.md.
  - `git diff 783d55b d3db5ac` for these files.
  - D61.md (i)–(v), as the input to the derivation in §5.
- **Run settings.**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv.
  - Heavy PIDs of mine:
    - 2099 and 2101: the two `a6_job` reruns.
    - 2395: the pytest run.
    - 2749, 2854, 3423 and 4074: my exact-error scans.
    - 4484, 4790 and 9159: three mutant slots.
    - The outlier and (1, MM) scans ran under a wrapper shell.
  - **One breach of the limit.** Analysis script PID 6029 (`floor3.py`, exact-rational transient fits) turned out heavier than I expected. It ran about 2 minutes beside the two mutant slots, so three of my heavy processes overlapped. I checked `/proc/6029/cmdline` and killed it; it was the only kill. It was redone in arb, where it took 1.3 s.
  - Every process I started had exited at handback.
- **The repository was not modified.** I ran only `git archive/diff/log/status`. The working tree's uncommitted changes (noise_resonance, D60) belong to other sessions.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my code, or the repository's functions run on inputs I chose.
  - [snippet]: web-search results only; arxiv.org does not resolve in this session.
  - [memory]: recalled facts.

## Summary

**Design.**
- Reproduced the md, the 66 tests, and the two new A6X records with `a6_job`.
- Recomputed every A6X number independently (flint `fmpq_mat` exact solves, t(v) in Arb at 800 bits).
- Attacked the open question with my own exact errors: python-flint arb, independent nodes, ~88 good digits, at every K = 2…400. Rules covered:
  - both proved Mercer frames;
  - three controls: (2, 2·v_hs), (2, MM) and (1, MM);
  - the frameless outlier;
  - two slow-branch frame rules, (2, 0.5·v_hs) and (2, 0.8·v_hs).
- Derived the next term of the ratio from D‴.61 (iv)'s Fock integral and tested it against my data and, in exact arithmetic, against the raw's stored errors.
- Ran the repository's 43 round-123 mutants and 25 of my own, with a second pass on my md-changing ones.
- Six web searches; one fetch, blocked.

**Verdicts.**
- **Reproducibility is exact.**
  - The md is byte-identical; 66 tests pass in 61.2 s.
  - Both new A6X records are bit-identical to the raw, every field: (2, Mercer) in 90 s, (1, Mercer) in 396 s.
  - My independent errors match the raw's 26-digit strings to 3.4·10⁻²⁶.
- **Every number D58 and the md add is right** [own, raw].
- **The open question has an answer.** The Mercer frames' stall comes from the frame's own (iv) correction series, not from a component outside the basis.
  - The excess's integer-power expansion from K⁻² holds there exactly. Its K⁻² coefficient is the one derived from D‴.61 (iv), to 10⁻²³ and 2·10⁻¹⁹.
  - But its coefficients grow about 12–25 times per order, against ≤ 7 at the four passing rules. So the orders whose smallest K is 80, 60 and 40 are limited by the truncation remainder at that K: the omitted powers reproduce each order with its sign.
  - At K ≥ 160 the orders decrease monotonically through 16 terms, to 7.3·10⁻²⁴ and 1.1·10⁻²¹ of 1/t(v). No K^{−1/2}, oscillatory or second geometric term is present there above 10⁻³⁶–10⁻⁴⁰.
  - The verdict's "the limit or the premise off" is therefore wrong at both known failures: the window is off.
- **No major issue.**

**Findings.**
1. **D3dx-1 (minor).** The stall's cause, and the verdict's misattribution of it (§4).
2. **D3dx-2 (minor).** The "negative result" is not a test. The carried-transient fits run down to K = 20. Resolvability is judged against the input digits, not against the power series' own remainder, which at K ≤ 60 is the transient's size and entangled with it. The printed B values are artefacts. The conclusion survives for another reason: the transient is below 10⁻²⁷ where the orders converge.
3. **D3dx-3 (minor).** The next term, derived (§5).
   - c₂ = s₁/2, with s₁ = 1/16 − ϱ/(4(1+ϱ)²) − 2ϱ(1−ϱ)/(1+ϱ)³ + √q(1−ϱ)²/(8(1−√q)(1+ϱ)²). Here ϱ = ±t(v), the sign by branch, and √q = h_other/t(v) taken signed.
   - Confirmed at both branches. It also holds, formally, at every frameless fast-branch rule, moment matching included (exact rationals −13/250, −11/128, −187/4000), and at the outlier.
   - With c₂ fixed, the stored data pass the check at all six non-transient rules, the 3 known cases among them.
4. **D3dx-4 (minor).** The outlier is judged by a different, weaker criterion than the other rules. Its transient-carried higher orders (computed, not printed) do not decrease. At consecutive K, with q^K(B₀ + B₁/K + B₂/K²) carried, the orders do decrease, to ~10⁻²⁴; B₀ = 0.238316 > 0 and B₁/B₀ ≈ −0.005.
5. **D3dx-5 (cosmetic).** Residues:
   - the md's D61_STATUS still lacks round 120 and "Arb's correctness assumed";
   - "changing sign, ranging over a factor 3.35·10⁴" merges two rules' facts;
   - the transient is "second branch's" in one place and "other branch's" in another;
   - `ct5` is dead code.
6. **D3dx-6 (cosmetic).** Coverage.
   - The repository's 42 round-123 mutants are all killed; the control survives.
   - 16 of my 24 are killed. My 9 md-changing kills are also killed with the md regression deselected.
   - 8 survive, all equivalent on the data. Two are latent: the bound index in `a6x_decreasing`, and a transient gate that would crash for 0.68 < q < 0.75.

## What I checked and how

**1. Reproduction** [md, raw, own]
- **The md.** `--from-raw` into W/out: `cmp`-identical to results/chaos/width_mismatch.md (md5 402859c3…; raw ae4a0d69…).
- **The tests.** The full file: 66 passed in 61.24 s.
- **The new records.** I ran `W/own/extra_job.py` = `WM.a6_job((r, 3, WM.A6X_K))` at τ/σ = 1 and 2, then compared every field with the raw's A6X record: r, name, v, v_hs, rate, arg, t, and per K the prec (550/328), e0, e1 (26-digit strings) and ratio. Both are BIT-IDENTICAL.
- **Independent errors.** `W/own/arberr.py` is my own implementation:
  - nodes by LAPACK then Newton in arb;
  - e_K² = 1 − 2√2(1+σ*²)^{−1/2}Σa_k e^{−y_k²/(2(1+σ*²))} + σ*⁻¹Σa_ja_k e^{−(y_j−y_k)²/4}, re-derived from the Gaussian inner products.
  - It matches all 18 stored errors per Mercer rule to ≤ 3.4·10⁻²⁶ relative (the strings' rounding).

**2. A6X recomputed independently** (`W/own/recompute.py`) [own, raw]
- **Orders.** Every printed order of the seven rules reproduces, and so do the carried-transient orders (−2.17·10⁻¹², +1.16·10⁻¹², −6.35·10⁻¹² and −1.25·10⁻¹³, +7.83·10⁻¹⁴, −1.34·10⁻¹³) with their B (107.9, 1.417, −2.249 and 1976, 3.496, 0.05898; factor 3.35·10⁴).
- **Input bound.** At most 6.70·10⁻²⁴, printed 6.8·10⁻²⁴ (rounded up).
- **Resolutions.** 4.088·10⁻¹⁸, 8.924·10⁻¹⁸, 1.864·10⁻¹⁹ and 2.201·10⁻¹⁸, printed 4.1, 9.0, 1.9 and 2.3 (rounded up).
- **The outlier.**
  - dev +6.5243·10⁻¹⁶ and dev5 +3.1066·10⁻¹⁵, so the resolution is 2.45·10⁻¹⁵, printed 2.5·10⁻¹⁵.
  - B_K/B − 1 at K = 20, 40, 60, 80: −3.8·10⁻², −1.66·10⁻³, −1.18·10⁻⁴, −2.50·10⁻⁵.
- **The bias ratio** is 29.0.

**3. The new checks' code** [md, own]
- **What I read.** `a6x_decreasing` (bound of the later order), `a6x_resolution`, `a6x_transient_ok` (B > 0, flat < 10⁻³, |dev| ≤ max(|dev5 − dev|, bnd), power fits worsening), the trigger per order (`rz`), the merge in `part_A6X`, `main`'s `--rules`, and `d61_step`'s assert.
- **Correctness.** All are correct as written.
- **Observations.**
  - `ct5` is assigned but unused; the check uses `ow[1]`.
  - The transient-carried fit uses one K more than the order whose smallest K triggered it, down to K = 20 for the eight-term order (D3dx-2).
  - At the outlier, `ow[2:]` (six to eight terms with the transient) are computed but neither printed nor checked (D3dx-4).
  - The md's D section prints D61_STATUS "VA verified by referee round 121" five times; D58 now says "rounds 120, 121 … Arb's correctness assumed" (D3dx-5).

**4. The open question: why the Mercer frames' orders stall** (`W/own/scan.py`, `analyze1–3.py`, `coefs.py`, `predict_orders.py`, `floor*.py`, `orders_mm.py`) [own]

*(a) The premise holds at K ≥ 160.* From my consecutive-K errors, windows on [160, 400], limit free:

| rule | orders, 4 … 16 terms | residual floor, N = 6 … 26 terms (max over all 241 K) | K⁻¹ coefficient |
|---|---|---|---|
| (2, Mercer) | +7.9·10⁻¹³ … −7.3·10⁻²⁴, monotone | 6.9·10⁻¹⁸ → 9.4·10⁻⁴¹ | 3·10⁻¹⁸ |
| (1, Mercer) | −4.4·10⁻¹¹ … −1.1·10⁻²¹, monotone | 1.0·10⁻¹⁶ → 2.1·10⁻³⁶ | 3·10⁻¹⁶ |
| (2, 2·v_hs), proved control | … 1.7·10⁻³⁶ | — | — |
| (2, MM), frameless | … 8.6·10⁻³⁸ | — | — |
| (1, MM), frameless | … 2.5·10⁻³⁵ | — | 4·10⁻²⁹ |

- The limits agree with the proved 1/t(v) to that resolution.
- No K^{−1/2}, oscillatory or other non-power component exists there above those floors.
- The other branch's transient is q¹⁶⁰ = 1.2·10⁻³⁵ and 4.5·10⁻²⁸.

*(b) The stall is the windows' small K, not the sparse grid.*
- My data at A6X's K reproduce it: +1.20·10⁻¹², +2.52·10⁻¹³, −1.25·10⁻¹³, +7.87·10⁻¹⁴, −1.38·10⁻¹³ and −6.99·10⁻¹¹, +8.53·10⁻¹², −2.18·10⁻¹², +1.21·10⁻¹², −1.74·10⁻¹².
- The four- to seven-term orders already decrease at both rules. Only the eight-term order, the one through K = 40, breaks the decrease.
- With all consecutive K on [40, 400] the Mercer orders still stall through 31 terms, near 10⁻¹⁵ (τ/σ = 2) and 2·10⁻¹³ (τ/σ = 1). On [80, 400] they stall near 10⁻²⁰ and 10⁻¹⁷.
- Moment matching on [40, 400] decreases to 1.6·10⁻³¹ (τ/σ = 2) and 1.7·10⁻²⁷ (τ/σ = 1).

*(c) The coefficients are large.* c_k of the excess ε_K = t(v)R_K − 1, limit exact, nine fits on [140–200, 400]:

| rule | coefficients, k = 2 … | reliability | ratio c_k/c_{k−1} |
|---|---|---|---|
| (2, Mercer) | −0.044, −0.235536, −1.091672, −3.560478, +28.26863, +1174.068, +2.84365·10⁴, +6.4438·10⁵, +1.4847·10⁷, +3.6·10⁸ | to c₁₀; c₁₁ to ±9% | 41.5, 24.2, 22.7, 23.0, 24.0 at k = 7–11 |
| (1, Mercer) | +0.0575471, +0.818862, +9.726926, +120.6617, +1658.57, +2.5669·10⁴, +4.49·10⁵ | to c₈ (±1.3%); c₉ ±50% | 14.2, 11.9, 12.4, 13.7, 15.5, 17.5, (19.8) |
| (2, 2·v_hs), control | \|c_k\| ≤ 2.4·10³ through k = 14 | — | 0.6–6.4 |
| (2, MM), control | \|c_k\| ≤ 2.8·10² through k = 14 | — | — |

- At K = 40 the Mercer terms c_k/40^k shrink only about 0.6× (τ/σ = 2) and 0.3–0.5× (τ/σ = 1) per order.

*(d) The anatomy of each order.*
- An n-term order's error is Σ_{k>n+1} c_k g_k, with g_k = Σ_j w_jK_j^{−k} the window's response.
- Partial sums over the known c_k reproduce every order with its sign:

| order | (2, Mercer): actual | predicted | (1, Mercer): actual | predicted |
|---|---|---|---|---|
| four-term | +1.20·10⁻¹² | +1.20·10⁻¹² | −6.99·10⁻¹¹ | −6.99·10⁻¹¹ |
| five-term | +2.52·10⁻¹³ | +2.52·10⁻¹³ | +8.53·10⁻¹² | +8.47·10⁻¹² |
| six-term | −1.25·10⁻¹³ | −1.23·10⁻¹³ | −2.18·10⁻¹² | −1.96·10⁻¹² |
| seven-term | +7.87·10⁻¹⁴ | +6.6·10⁻¹⁴, still rising | +1.21·10⁻¹² | +5.7·10⁻¹³, rising |
| eight-term | −1.38·10⁻¹³ | −4.7·10⁻¹⁴, rising | −1.74·10⁻¹² | −1.3·10⁻¹³, rising |

- The eight-term order's remainder at K = 40 is T₄₀ = 6.5·10⁻⁹ and 8.9·10⁻⁸, against a weight of −7.8·10⁻⁵ at that K.

*(e) The transient.*
- At K = 40, q⁴⁰ = 1.9·10⁻⁹ and 1.5·10⁻⁷, the size of T₄₀.
- Carrying q^K(B₀ + … + B₄/K⁴) on [40, 400] or [60, 400] lowers the floor by 1–2 decades. But B₀ wanders over 0.06 … 1.9 and B₁ over −8 … −324 across orders, so it is not a resolvable smooth-prefactor term there.

*(f) A model for the growth.*
- I took the leading-order angular (Legendre) model of D‴.61 (iv)'s integral, ∝ (K!/(2K)!)²Γ(2K+1)∫(γx cos² + γy sin²)^{−(2K+1)}dθ, exact via P_{2K}, K ≤ 2400.
- Its ratio's coefficients are Gevrey-1: the increments of c_k/c_{k−1} converge to 1.98–1.99 and 2.53–2.54 at k = 15–22, against 1/(−ln q) = 1.9895 and 2.5408.
- So the Borel singularity sits at the other axis's exponent, the transient's q^K, consistent with Dingle/Berry–Howls resurgence [snippet, memory].
- The true τ/σ = 1 increments, 1.3, 1.8, 2.0, 2.3 at k = 5–9, approach 2.54. Suggestive, not established.

**5. The next term derived from D‴.61 (iv), and tested** (`cumulants.py`, `nextterm.py`, `closed.py`, `c2*.py`) [own, raw]

*The 1/K term of ĝ.*
- From D‴.61 (ii)'s exact moments I built P_j(K) = D(2K+2j)/K! as polynomials. They were interpolated on K = j+1 … 3j+1 and checked on two more K; leading coefficient 1/j!, and P₁ = K² + K + 1.
- Series arithmetic in x = 1/K then gives the cumulants: κ_{2j}(S) = c_j x^j + d_j x^{j+1} + …, with c_j = −(2j−1)!/4^j (D‴.61 (iii)(c), confirmed) and **d_j = −c_j·j(8j² + j + 1)/4, exact for j ≤ 22**.
- Hence log ĝ_K(√Kη) − Kη²/4 = ½log(1 − y) + F₁(y)/K + …, with **F₁(y) = y(1+y)/(1−y)³ + y/(8(1−y)²) + y/(8(1−y))**, y = η²/4.

*The Laplace step.*
- Two-dimensional Laplace to O(1/K) at the mass point |η*|² = 4t(v): Hessian diag(2(γp−γs), 4γs); third and fourth derivatives of 2log|η|²; the curvature of |1 − η²/4|; F₁; Stirling's 1/(12K).
- This gives S_K = S_K⁽⁰⁾(1 + s₁/K + …) and e_{K−1}/e_K = (1/t(v))(1 + c₂/K² + O(K⁻³)) with c₂ = s₁/2. The closed form is in finding 3; p = σ_p²/|η*|² = √q/(4(1−√q)), the angular width, which grows without bound at the cusp.
- At (2, Mercer): 1/12 + 5/12 − 9/40 − 27/400 − 591/2000 = −11/125, so c₂ = −11/250.
- No K⁻¹ term arises off the cusp.
- **This is formal at next order.** D‴.61 proves only the leading limit of ĝ; the next-order uniformity and the Laplace remainder are not proved.

*Tests against my data* (relative agreement of c₂):

| rule | c₂ derived | agreement |
|---|---|---|
| (2, Mercer) | −0.044 | 1.1·10⁻²³ |
| (1, Mercer) | 0.0575470734670659269 | 1.7·10⁻¹⁹ |
| (2, 2·v_hs) | −0.051243224484370181 | < 10⁻³⁰ |
| (2, 0.8·v_hs), slow branch | 5.7308947080 | 9.8·10⁻⁹ |
| (2, 0.5·v_hs), slow branch | 21.988001779 | 6.8·10⁻⁵, fits rising toward it |
| (1, MM), frameless | −11/128 | exact |
| (2, MM), frameless | −187/4000 | exact |
| outlier, formal | −0.0924369845420 | 4.9·10⁻¹⁰ (fit −0.0924369845874; B's 1/K term unmodelled) |

*Tests against the stored errors* (exact rationals, limit exact, 2–4-term fits):
- Frame rules: within 1.2·10⁻⁹ … 1.6·10⁻⁵ at 7 of the 18. These are the fast branch off the cusp and the Mercer-zero rules at τ/σ = ¼ and ½. At the other 11 (the slow branch, where c₃/c₂ = −82 … −263, and the rules near the cusp) the fits move toward it, 7·10⁻⁴ … 0.27 off.
- Frameless rules: within 2.7·10⁻¹⁰ … 2.4·10⁻⁹ at 5 of the 6 fast-branch rules, (¼, 2·v_hs) included; (¼, MM) is within 5.0·10⁻⁷.
- This large c₃/c₂ on the slow branch is why D58's slow-branch "K^{−1.52 … −1.98} at K = 160" is pre-asymptotic.

*The check with c₂ known* (the stored data; the same top powers, one K fewer):

| rule | orders | resolution | the code's resolution |
|---|---|---|---|
| (2, Mercer) | −5.14·10⁻¹³, −4.39·10⁻¹⁴, +1.20·10⁻¹⁴, −4.41·10⁻¹⁵, +2.75·10⁻¹⁵ | 7.2·10⁻¹⁴ | fails |
| (1, Mercer) | +2.41·10⁻¹¹, −1.71·10⁻¹², +2.22·10⁻¹³, −6.47·10⁻¹⁴, +4.49·10⁻¹⁴ | 1.1·10⁻¹² | fails |
| (½, MM) | last 3.8·10⁻²¹ | 4.6·10⁻¹⁹ | 4.1·10⁻¹⁸ |
| (1, MM) | last 6.2·10⁻²¹ | 7.2·10⁻¹⁹ | 9.0·10⁻¹⁸ |
| (2, MM) | last 5.4·10⁻²² | 2.8·10⁻²⁰ | 1.9·10⁻¹⁹ |
| (2, 2·v_hs) | last 1.6·10⁻²¹ | 1.8·10⁻¹⁹ | 2.3·10⁻¹⁸ |

- Every order decreases, so the check passes at 3 of 3 known cases.

**6. The outlier at consecutive K** (`outlier.py`; frameless, q = 0.8506295) [own]
- Power-only orders on [160, 400] stall near 10⁻¹², as they must with q¹⁶⁰ = 5.6·10⁻¹².
- With q^K(B₀ + B₁/K + B₂/K²) on [200, 400] the orders fall to ~10⁻²⁴–10⁻²⁶ (scatter to 10⁻²²): 1/t(v) to ~10⁻²² at a frameless rule.
- B₀ = 0.238316–0.238317, stable and > 0. B₁/B₀ = −0.005 … −0.006: a small 1/K prefactor. The code's four-term B is 0.238302.

**7. Tests and mutants** (`W/mut/mutants129.py` = tools/mutants.py with REPO = W/repo and a per-slot `--basetemp`; `--extra gh_peano_raw.npz`) [own]
- **The repository's spec, round123_D58_own.py.**
  - J00 survives, as a control should; J01–J42 are all killed.
  - 23 of them changed the md, and the md regression killed those first. I did not re-run those 23 without it.
- **My spec, W/mut/spec129.py (Y00–Y24).**
  - Y00 survives, as a control should.
  - Killed: Y02–Y06, Y11–Y14, Y16, Y18, Y19, Y21–Y24.
  - The 9 md-changing kills (Y02, Y04, Y05, Y06, Y12, Y13, Y16, Y19, Y24) are killed again with `-k "not round88_precision_and_the_md"`.
  - Survivors, all equivalent on the data:
    - Y01: the earlier order's bound; latent, the tests' bound vectors are equal at the compared positions.
    - Y07: the transient analysis gated at the five-term order's K; for 0.68 < q < 0.75 it would raise IndexError at `ct[4]`, an untested path.
    - Y08: `ct5`, dead code.
    - Y09: `bnd` from the next order.
    - Y10 and Y15: strict versus non-strict comparisons.
    - Y17: the bias ratio's minimum; only one τ/σ has two proved rules.
    - Y20: equivalent by construction, a flaw of mine.

**8. Literature** (six searches; arxiv.org unreachable) [snippet]
- Berry–Howls, hyperasymptotics for integrals with saddles: exponentially small contributions of other saddles, the remainder through "adjacent" saddles, Dingle's resurgence relating late to early terms.
- Berry's lecture: in resurgence, low orders reappear in high orders.
- A paper on factorially divergent inverse-power series favours an error ansatz over a power-series remainder.
- Stetter (1968) treats Richardson extrapolation for well-behaved expansions only.
- Nemes and Wojdylo give coefficients of one-dimensional Laplace's method; Fukuda–Kagaya–Ueda and Katsevich give multidimensional expansions. I did not see a formula in any of them.
- Nothing found on the Gauss–Hermite rule's frame series at next order.

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) :14–18, (ii) :27–29, Status :174–177, :298–301: VA attributed to rounds 120, 121, its adoption to round 122, Arb's correctness assumed | **Right** (D3dr-5 fixed in D58; the md string is a residue) | §3; D3dx-5 |
| :271–276: A6X extended by 3 known cases; at MM ½, 1, 2 and (2, 2·v_hs) the orders decrease within ten times the last change (≤ 9.0·10⁻¹⁸), input bound ≤ 6.8·10⁻²⁴ | **Right.** Confirmed far finer at τ/σ = 1, 2 by consecutive K (to 10⁻³⁵–10⁻³⁸) | §2, §4a |
| :277–279: the check fails at 2 of 3 proved rules; the Mercer stalls (values) | **Numbers right**; "a failure says nothing of the limit" right | §2, §4b |
| :279: "the limit or the premise off" | **Wrong at both known cases.** Limit and premise hold; the windows' smallest K is outside the expansion's range | §4, §5; D3dx-1 |
| :280–282: the transient resolvable, carried, the orders still do not decrease, B values; "a negative result; its cause is open" | **Numbers right.** The conclusion is right, but the test is not one; the cause is found | §4d–f; D3dx-1, -2 |
| :283–287: window geometry; the 29.0 ratio from the raw | **Right** | §2 |
| :262–271: the outlier's B_K at 60, 80 within 1.2·10⁻⁴, B > 0, within 6.6·10⁻¹⁶ inside 2.5·10⁻¹⁵ | **Numbers right.** A different criterion than the other rules'; B carries a small 1/K term (B₁/B₀ ≈ −0.005); its limit holds to ~10⁻²² at consecutive K | §2, §6; D3dx-4 |
| Numerics (A6X) :355–360 | **Right** | §2 |
| (ii)/Status: frameless rules conjectural | Still conjectural. New support: the next coefficient c₂ also matches there (MM exactly) | §5; D3dx-3 |
| (iii)–(v), proofs, B, C, D | Unchanged except the tier wording; tests pass; 22 unchanged numbers spot-checked in the md | §1 |

## Status of the round-123 issues

| D3dr | status | note |
|---|---|---|
| 1 | **fixed**, with a residue | Each order's trigger is at its smallest K; the transient is carried. Residue: the carried fit reaches one K lower (to K = 20), and "resolvable" is measured against the input digits, not the series' remainder (D3dx-2) |
| 2 | **fixed**, with a residue | (1, Mercer) and (2, Mercer) are added and reproduce bit-identically; the verdict is per rule; the ratio 29.0 comes from the raw; the caveats are computed. Residue: the failure's attribution (D3dx-1) |
| 3 | **fixed** | The input bound (my 6.70·10⁻²⁴), the resolution as ten times the last change, the non-integer-power caveat |
| 4 | **fixed** | B > 0; \|dev\| ≤ \|dev5 − dev\| (accepts offsets −3.1 … +1.8·10⁻¹⁵); B_K at 60 and 80, outside the fit. Residues in D3dx-4 |
| 5 | **fixed in D58** | The md's D61_STATUS is unchanged (D3dx-5) |
| 6 | **fixed** | "26-digit" throughout, the docstrings, the heading, one rounding per quantity, "not reproduced" marked, `d61_step` asserts t < 1 |
| 7 | **fixed** | Plants at the edges, the q = 0.6 and transient mocks, the non-dyadic moment-matching rule, the D header. The repository's 42 mutants are killed; my 8 survivors are equivalent on the data (D3dx-6) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dx-1 | minor | **The open question has an answer, and the NOT verdict misattributes it.** Both known failures are the windows' small K, not the limit or the premise [own]:<br>• On [160, 400], 16 terms, the orders decrease to 7.3·10⁻²⁴ and 1.1·10⁻²¹. The residual floor is ≤ 9·10⁻⁴¹ and 2·10⁻³⁶, the K⁻¹ coefficient ≤ 3·10⁻¹⁸ and 3·10⁻¹⁶, and c₂ is as derived.<br>• The coefficients grow 12–25 times per order (≤ 7 at the passing rules), so the omitted powers reproduce every A6X order with its sign.<br>• Only the eight-term order, through K = 40, breaks the decrease. With consecutive K on [40, 400] the stall persists (≈ 10⁻¹⁵, 2·10⁻¹³); moment matching there converges to 10⁻³¹.<br>• The leading-order model's coefficients are Gevrey with slope 1/(−ln q). | **D58 :277–283:** replace "(a negative result; its cause is open)" with the cause and these numbers.<br>**`a6x_lines`, the NOT branch:** "the limit, the premise, or the window — the expansion not yet asymptotic at the order's smallest K".<br>**A6X:** more K beyond 160 (python-flint, already in the venv, gives K = 2 … 400 in 2–7 min per rule [own]) with windows K ≥ 160, or c₂ from D3dx-3 |
| D3dx-2 | minor | **The negative result's test is not a test.**<br>• The carried fits pass through K = 20, 40, 60. At K = 20 the terms c_k/20^k do not decrease; at K = 40 the power remainder T₄₀ = 6.5·10⁻⁹ and 8.9·10⁻⁸ is the size of q⁴⁰ = 1.9·10⁻⁹ and 1.5·10⁻⁷.<br>• Even with B(1 + … + B₄/K⁴) the B's wander (B₀ 0.06–1.9, B₁ −8 … −324). The printed B (108, 1.42, −2.25; 1976, 3.5, 0.059) are artefacts.<br>• "Resolvable" (q^K > 10⁻²⁰) is measured against the input digits, not the competitor.<br>• The conclusion stands because at K ≥ 160, where the orders converge, q^K < 10⁻²⁷ | **`a6x_lines`:** carry Bq^K only where q^K exceeds the order's own remainder (its next change at that K), and never below the power order's smallest K.<br>**D58 :280–282 and the md:** restate the negative result with the K ≥ 160 evidence |
| D3dx-3 | minor | **The next term, derived and confirmed, but not used.**<br>• c₂ = s₁/2 with the closed form above (from d_j = −c_j j(8j²+j+1)/4, F₁ closed, the 2-D Laplace step); formal at next order.<br>• Own data: 1.1·10⁻²³, 1.7·10⁻¹⁹ and < 10⁻³⁰ at the fast-branch frames; 9.8·10⁻⁹ on the slow branch. Raw: 1.2·10⁻⁹ … 1.6·10⁻⁵ at 7 frame rules; frameless 2.7·10⁻¹⁰ … 2.4·10⁻⁹ at 5 of 6 fast-branch rules, moment matching exact.<br>• With c₂ fixed the stored data pass the check at all 3 known cases and sharpen the moment-matching resolutions 7–13 times.<br>• The slow branch's c₃/c₂ = −82 … −263 explains the "K^{−1.52 …}" | **D58 Status (ii) and the A6X passage:** state the derived next term, "formal at next order", and that it holds at moment matching (−13/250, −11/128, −187/4000) and the outlier.<br>**exp_width_mismatch.py:** `a6x_c2`, and an optional known K⁻² term in the orders.<br>**Tests:** pin c₂(2, Mercer) = −11/250 and c₂(½, MM) = −13/250 |
| D3dx-4 | minor | **The outlier is judged by a weaker criterion.**<br>• Its transient-carried orders are +6.5·10⁻¹⁶, +3.1·10⁻¹⁵, −6.0·10⁻¹³, −7.8·10⁻¹², −7.8·10⁻¹² (`ow`, the last three neither printed nor checked); the rule passes on \|dev\| ≤ \|dev5 − dev\| alone.<br>• "No power of K beside q^K" holds only at leading order: B₁/B₀ ≈ −0.005 [own, raw] | Print the outlier's carried orders and name its criterion; or use windows K ≥ 160 with q^K(B₀ + B₁/K), where my data reach ~10⁻²⁴ with B₀ = 0.238316 > 0 |
| D3dx-5 | cosmetic | **Residues.**<br>• The md's D61_STATUS reads "VA verified by referee round 121", five times, without round 120 or "Arb's correctness assumed".<br>• "its B across those orders changing sign, ranging over a factor 3.35e+04" merges two rules (sign change at τ/σ = 1, factor at 2).<br>• "second branch's" and "other branch's" transient are both used.<br>• "a failure says nothing of the limit — the limit or the premise off" offers the limit as a cause in the same breath.<br>• `ct5` is dead | Fix the strings (the D61_STATUS test pins it; update both); attribute B's behaviour per rule; drop `ct5` |
| D3dx-6 | cosmetic | **Coverage.** The repository's 42 mutants are all killed. 8 of my 24 survive, equivalent on the data:<br>• latent: Y01 (bound index) and Y07 (IndexError for 0.68 < q < 0.75);<br>• equivalent: Y08 (dead code), Y09, Y10, Y15, Y17, Y20 (equivalent by construction) | Plants: a bound vector unequal across orders (Y01); a rule with q = 0.70 through `a6x_lines` (Y07) |

## Sources

- **Berry–Howls hyperasymptotics** [snippet]. These search results were summarised together; I cannot attribute a statement to a single link:
  - https://eprints.soton.ac.uk/496796
  - https://ir.cwi.nl/pub/2293/2293D.pdf
  - https://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/60642/1/0968-2.pdf
  - https://library.wolfram.com/infocenter/Articles/3153/
  - https://arxiv.org/pdf/2205.07163v1
- Berry, "Divergent series: from Thomas Bayes's bewilderment to today's resurgence via the rainbow", https://www.talks.cam.ac.uk/talk/index/54642 [snippet].
- Exponential asymptotics and problems with coalescing singularities, https://arxiv.org/pdf/1403.7182 [snippet].
- Asymptotic Approximations to Truncation Errors of Series Representations for Special Functions, https://arxiv.org/pdf/math/0511074 [snippet].
- Stetter, Apl. Mat. 13 (1968), https://dml.cz/bitstream/handle/10338.dmlcz/103154/AplMat_13-1968-2_13.pdf [snippet].
- Sidi, BIT, https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P059_BIT.gen_Rich.extrap.further_results.pdf [snippet].
- Nemes, An explicit formula for the coefficients in Laplace's method, https://arxiv.org/pdf/1207.5222 [snippet].
- Wojdylo (J. Comput. Appl. Math. 196; SIAM Rev. 48, 2006), https://link.springer.com/article/10.1007/s00365-013-9202-6 [snippet].
- Fukuda–Kagaya–Ueda, https://arxiv.org/html/2504.01310v3 [snippet; fetch blocked].
- Katsevich, https://arxiv.org/html/2602.23151v1 [snippet].
- **Repository** [md, raw]: D58.md, D61.md, exp_width_mismatch.py, width_mismatch.md and its raw, tests/test_width_mismatch.py, tools/mutants/round123_D58_own.py and tools/mutants.py, REVIEW_round123_D3dr.md, REVIEW_LOG Round 123.
- **My files.**
  - Scripts, all in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/own/`:
    - exact errors: arberr.py, scan.py and its outputs scan_*.txt;
    - reproduction: extra_job.py, cmp_raw.py, recompute.py;
    - the stall: analyze1.py, analyze2.py, analyze3.py, coefs.py, predict_orders.py, floor.py, floor2.py, floor3.py, orders_mm.py, legmodel.py;
    - the next term: cumulants.py, f1check.py, nextterm.py, closed.py, c2raw.py, c2raw_all.py, c2formal.py, c2scan.py, c2fixed_all.py, c2one.py;
    - the outlier: outlier.py, bk.py;
    - outputs x_mercer1.json and x_mercer2.json.
  - Mutant spec and runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/mut/spec129.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/mut/mutants129.py`.
  - Mutant results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/mut/res_J.jsonl`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/mut/res_Y.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/mut/res_Y2.jsonl`.
  - Logs: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref129/logs/`.
- [memory]: the late terms of a Laplace expansion grow like Γ(k + β)/F^k, with F the exponent gap to the adjacent critical point (Dingle, Berry–Howls).

## Lesson

**Move the window before inventing a component.** A stalled sequence of extrapolations says something about the window as well as the function. Here the same fits on K ≥ 160 converged through sixteen orders. So the "stall" was the windows reaching K = 40, where this expansion is not yet asymptotic: its coefficients grow twentyfold per order.

**Exact arithmetic was cheap enough to settle the question with data, not hypotheses.** With python-flint, every K up to 400 took minutes, against seven minutes per three K in the repository's Decimal code.

**Ask the theory for the next term.** D‴.61's Fock integral, carried one order further, predicted the K⁻² coefficient to 23 digits. That calibrates the check better than any known case can: with that coefficient fixed, all three proved rules pass. It also carried over, formally, to the unproved frameless rules.
