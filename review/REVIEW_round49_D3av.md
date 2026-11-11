# Referee report D3av (round 49): Proposition D‴.54, revised after round 47 (D3at)

**Summary.** The fixes to (ii)–(iv) and to the bookkeeping hold. The md regenerates byte for byte from the raw, all 17 tests pass, and
every number in the numerics paragraph matches the md.

Part (v), rebuilt as a fold law, is still wrong in two ways:
- **The leading-order law has the wrong factor.** The binding term is projected off φ but the pull term is not. This is D3at-1's failure
  mode again: two balanced terms expanded under different constraints.
- **The fold is not the capture threshold.** Hold the coefficient w fixed and the runner-up has a rest point. Free w (as the loop does)
  and that rest point is not stationary: ∂E/∂w > 0 along the whole branch. The held runner-up's coefficient decays, on an inner time that
  does not depend on ε, until it reaches the fold. Capture is then a race against the unfitted feature's own entry, so at a fixed gap the
  outcome flips with ε:
  - x = 4.9: bound at ε = 3.75·10⁻⁹ and 10⁻⁹, captured at 10⁻¹⁰.
  - C3's "stranded" gaps 20.6 and 24.5: captured at 10⁻¹¹ and 10⁻¹².

  So C3's "window in x" is one line (μ = 3·10⁻⁴) through the (x, ε) plane. (iv)'s restriction, "every runner-up's coefficient above
  w_f", fails; C3's own x = 5.07 row (captured at 1.11 w_f) already shows it.

Other findings:
- At the fold the minimum's displacement is 0.76–0.94σ, not 0.3–0.4σ.
- c(d) depends on A′/A.
- "Stranded" tails actually turn back toward their own feature.
- The "788 doublings" argument is not valid, although the counts it protects do survive.
- Part D's stated reason is wrong: the finite-ε shift dominates.
- T3's count cannot tell the offset H/√2 from H/2.
- 18 of my 26 mutants survive.

## 1. What I checked and how
- **Setup.** A scratch copy with its own venv; separate script and mutation directories. Every Python process ran with
  OMP_NUM_THREADS=1 under nice -n 10, one heavy process at a time, the longest about 7 minutes.
- **Repository untouched.** Nothing in it was modified; `cmp` confirms the seven D54 files there equal the copy I reviewed. I ran one
  `git status` in the repository (read-only); all other git use was in the copy. I killed no processes.
- **Regeneration and tests.** `--from-raw` reproduces the md byte for byte; 17 of 17 tests pass in 10.4 s.
- **Numerics paragraph.** Every number checked against the md and the raw.
- **(v), by hand and by an independent solver** (`vstatic.py`, `v2d.py`). The solver uses closed-form Gaussian inner products, Newton's
  method, and the saddle-node system ∇E = 0, det H = 0. With it I computed:
  - the fold for d = 5–7σ and A′/A = 0.44–3.5;
  - the A′ → 0 limit of the rest point against both laws;
  - ∂E/∂w along the branch toward the unfitted feature and the branch away from it, at d = 5, 6, 8, 10σ;
  - the saddle-node structure of the fold;
  - the planar Hessian with the unfitted feature at 0°, 17° and 45°.
- **C3 dynamics** (`c3run.py`, `c3runmu.py`: c2_job's loop, with the runner-up's trajectory recorded), all at rtol 10⁻¹² unless stated:
  - at μ = 3·10⁻⁴, the gaps x = 4.5, 4.74, 4.9 inside C3's unresolved interval;
  - at fixed x across ε: x = 4.9 at 10⁻⁹ and 10⁻¹⁰, x = 4.5 and 4.0 at 10⁻¹⁰;
  - the upper edge: x = 20.6 at 10⁻¹¹, x = 24.5 at 10⁻¹²;
  - three boundary cases repeated at rtol 10⁻¹³;
  - x = 17.3 continued to 5.1 τ_e.
- **B4u.** Ten unsettled problems (the four tightest margins and six whose drifts were still growing) continued to 24 doublings
  (`b4ext.py`, a replica of `inner_k` that records shares per doubling; it reproduces the raw at doubling 12).
- **Part D.** The loop's edge share at gaps 5.8–6.7, offsets 0 and 3H/8, ε = 10⁻⁵, 10⁻⁶, 10⁻⁸ (`dshift.py`). The offset rule against the
  offset-blind rule on the raw positions.
- **T3.** The 8-kernel shares against both products, corner by corner.
- **Context.** D‴.51 and D‴.53 in PARADIGM.md, to check (i)'s citations.
- **Mutants** (`mutate.py`). 26 mutants of `c3_outcome`, `w_unbind`, `formed_by_feature`, `b4u_note`, `inner_k`'s stop logic and the
  report's V and C3 formulas, each run against the tests and the md regeneration.

## 2. Verdicts
**(v), the static law.**
- *What is right.*
  - The first-order terms cancel when the first moment is zero. With c_e free the first moment shifts by an amount ∝ A′e^{−d²/4σ²}: at
    the fold the entrant sits +0.056σ, +0.005σ, +0.0002σ from the compensating position for d = 5, 6, 7σ. That is higher order.
  - B = π^{n/2}σ^{n−4}/2 is right.
- *What is wrong: the pull.* The slaved coefficient projects the residual off φ, so the pull is ⟨r′,(∂²_eφ)_⊥⟩ =
  A′π^{n/2}σ^{n−2}(d/σ)²e^{−d²/4σ²}/4. That gives M₂ = σ²ℓ(d) and |δ|² = σ²ℓ(d)·w_e/(w(w_e + w)), with the same ℓ as the fold law. The
  transverse pull is 0, not negative; "rest on e" still holds.
- *Verified numerically.* As A′ → 0, the exact rest point divided by the draft's law tends to (d²/(d²−2))^{1/2} = 1.0426, 1.0290, 1.0211
  for d = 5, 6, 7σ (at A′ = 7·10⁻⁴: 1.0459, 1.0300, 1.0213). Divided by the corrected law it tends to 1 (1.0032, 1.0010, 1.0003).
- At w = 0.3 the corrected law is low by 15.8, 7.6, 3.7, 1.8, 0.8%, not 19–3%.

**(v), the fold.**
- *Reproduced.* w_f = 0.14239, 0.059565, 0.020663, 0.0060185, 0.0014952, so c = 4.215, 5.414, 6.644, 7.867, 9.110.
- *It is a saddle-node of the fixed-w problem.* The Hessian has a zero eigenvalue there. At w = 1.001 w_f (d = 6σ) there is a minimum at
  0.838σ (eigenvalue +4·10⁻⁶) and a saddle at 0.885σ (eigenvalue −4·10⁻⁶).
- *Wrong:* at the fold the displacement is 0.94, 0.91, 0.86, 0.81, 0.76σ. The quoted "≈0.3–0.4σ" is the law's extrapolation.
- *c depends on A′/A.* By homogeneity c = c(d, A′/A):

  | d | A′/A = 0.44 | 0.875 | 1.75 | 3.5 |
  |---|---|---|---|---|
  | 5σ | 4.50 | 4.21 | 3.76 | 3.14 |
  | 6σ | 6.73 | 6.64 | 6.48 | 6.18 |
  | 7σ | 9.12 | 9.11 | 9.09 | 9.05 |

- *The away branch never folds:* the mirror minimum on the far side persists down to w = 10⁻⁷.
- *Dimension-freeness* holds exactly for the static continuum problem: the Gaussian separates, and the planar Hessian stays positive up to
  the fold at 0°, 17° and 45°. V's planar run checks this identity; it is not an independent test.

**(v), dynamics: the fold is not the capture threshold.**
- *No stationary two-kernel state with w free.* On the toward branch ∂E/∂w > 0 at every w, in units of E/K:
  - d = 5σ: 1.2·10⁻⁴ … 4.4·10⁻³;
  - d = 6σ: 3.7·10⁻⁷ … 4.6·10⁻⁴;
  - d = 8σ: 5·10⁻¹⁴ … 4·10⁻⁷.

  On the away branch it is negative.
- *What the loop does:*
  - After the inner problem the runner-up coincides with the entrant, 4.7·10⁻⁴σ from the feature.
  - At s ≈ 3·10⁴ it departs to the static rest point. Held displacements in the loop are 0.343/0.393/0.452/0.518σ against the static
    0.356/0.410/0.477/0.537σ, and 25–40% above the leading-order law.
  - While held, its coefficient falls. At x = 4.9: 0.0216 → 0.0201 at ε = 3.75·10⁻⁹, and 0.0215 → 0.0169 (below w_f) at 10⁻⁹.
  - It is released when feature 2's own entrant enters (τ = 1.125 τ_e).
- *The decay does not depend on ε in s, so at fixed x a smaller ε captures.* For example, x = 4.9 is captured at 10⁻¹⁰ with w = 1.25 w_f.
- *Inference, not measured at small x:* at fixed x the capture boundary in x tends to 0 as ε → 0.

**C3.**
- *Classification.* `c3_outcome` gives the right final counts on all 30 rows, but:
  - "bound" is a transient hold;
  - "stranded" is a horizon snapshot. Continued to 5.1 τ_e, the x = 17.3 tail turns back toward feature 0 (1.21σ → 0.70σ,
    w 7.4·10⁻⁵ → 5.5·10⁻⁵);
  - the rms cut of 10⁻⁸ is untested, and C2's ε = 10⁻¹⁰ residual (4.3·10⁻⁸) equals the bound loops' finite-mass residual.
- *Tolerance:* converged; rtol 10⁻¹³ agrees on x = 5.07 and on x = 4.9 at both masses.
- *Mass:* not converged. Both edges move with ε at fixed x (x = 4.9; x = 20.6 and 24.5).
- *The 4.26–5.07 gap:* at μ = 3·10⁻⁴, x = 4.5, 4.74 and 4.9 are bound and 5.07 is captured. So the edge lies in (4.90, 5.07), with the
  runner-up's coefficient in (1.11, 1.25)·w_f, not at the fold's 5.22. The fold does not explain it, and "captures from 1.1 times the
  fold" contradicts "above w_f ⇒ (iv)".
- *5σ rows:* the md column shows 0.7 on 9 of the 15 rows (also 0.6998, 0.3395). That is the coefficient of the kernel already captured at
  1.005 τ_e, not the tail's.

**(iv)'s restriction does not follow.** Its definition through the fold is falsified (x = 5.07 at 1.11 w_f; x = 4.9 at 1.25 w_f with
ε = 10⁻¹⁰), and at fixed x the capture range grows as ε → 0.

**T3.**
- H/√2 = (n − 1)^{1/2}h/2 is the right offset for every edge at a vertex (the cube's centre).
- The shares to 1.2·10⁻³ and the count of 6 of 7 reproduce. But H/2 also gets 6 of 7 (it misses (1, 2, 3) instead of (0.5, 3, 5)), so
  only the shares favour H/√2.
- *The seventh triple:* the (0,0,1) corner has share 0.01012 in the 8-kernel problem against 0.00997 in the product. The formed set stays
  down-closed, so the hierarchy is not violated.
- But the 8-kernel shares exceed the product on every corner with gap ≥ 5, in all seven triples: +1.5% at gap 5, +4.6–7% at gap 6,
  +40–56% at gaps 8–8.5, ×4 at gap 12. The product's formation boundaries are biased low near θ systematically, the same direction as
  B4's three planar disagreements (4 against 3). The miss is not "a share at θ" by accident.

**B4 and B4u.**
- *Areas reproduced:* 33.516 and 15.453, giving 1112ε and 512.7ε.
- *Grid:* the product rule at steps 0.25 / 0.05 / 0.025 gives N ≥ 3: 33.516 / 33.518 / 33.511 and N = 4: 15.078 / 15.261 / 15.270
  (−1.3% at 0.25). The 4-kernel problem's own grid error is not measured.
- *"788 doublings":* margin divided by drift is not a bound while the drifts grow. At (6, 10) the drift doubled every doubling, from
  6.8·10⁻¹⁰ to 3.0·10⁻⁷ at the cap; at that rate the 2.4·10⁻⁴ margin lasts about 10 doublings.
- *Continuation:* the drifts peak at doublings 11–13, and the total change after the cap is ≤ 1.7·10⁻⁶; no count changes. The
  conclusion holds; the argument does not.
- *`remaining`* is the geometric (Aitken) tail of the last two max-norm drifts: an extrapolation, not a bound. It is infinite while drifts
  grow and overestimates by 4.7–6.2× where checked: (6, 8) gives 1.74·10⁻⁶ against 2.8·10⁻⁷ actual; (6, 7.25) gives 8.6·10⁻⁷ against
  1.8·10⁻⁷.

**Part D.**
- The numbers match the md.
- On the 7 positions where the two rules differ (4 at 10⁻⁵, 3 at 10⁻⁶), the loop sides with the offset-blind rule every time.
- The cause is the finite-ε threshold shift, measured at both offsets: +0.38–0.43 at 10⁻⁵, +0.12–0.14 at 10⁻⁶, +0.01 at 10⁻⁸. The
  variation of x_θ along a face is 0.27, so D measures the shift; "too few positions" is not the reason.
- D does not count against (iii)(c): at 10⁻⁸ the loop's crossings are about 6.31 (t = 0) and 6.16 (t = 3H/8), against x_θ = 6.300 and
  6.149.

**(i)–(iv) as revised.**
- (ii)(a) and (b) are right; I re-derived δ_j and μ_j, and the examples (ratio 1.01, ratio 2 = one spacing) check.
- (ii)(c)'s formulas and W1/W3 are right. The bound 8qσ²/|c_j − c_e| holds only for neighbours that differ in one coordinate.
- (iii)(d)'s multiplicativity holds, and "equal masses" is now stated.
- E[N − 1] (71.5, 105.8, 139.0, 171.3) and the face averages are right; so are the codimension law and the hexagonal factor 32σ⁴/H⁴.
- (i) attributes "the rate ε^{1/2}" to D‴.53(a)–(b); it belongs to D‴.53(c)–(d). Its list of hypotheses omits (H3).
- The log's "(ii)(c): separation at the entries is not enough" does not appear in the draft.

**Code and tests.**
- `inner_k`'s stop reasons work; `remaining` is as described under B4u.
- 8 of 26 mutants are killed. Survivors:
  - `b4u_note`: 0 of 5 killed, 4 of them change the md;
  - `c3_outcome`: the rms cut ×100 survives and changes 2 md rows;
  - the report's formulas: `v_law` with the projected pull, C3's fold column with d for d², and V's ℓ(d) with d all survive (6, 30 and 5
    md lines changed);
  - `inner_k`: the post-chunk clock check, the zero-drift branch and ≤ vs < all survive;
  - equivalent on the stored data: `formed_by_feature`'s amplitude and captured-set choices, and two `c3_outcome` reorderings.

**Round-47 issues.**

| issue | status | note |
|---|---|---|
| D3at-1 | partly fixed | the pull is unprojected (D3av-2); the rest point drifts (D3av-1) |
| D3at-2 | partly fixed | the restriction is falsified (D3av-1); (ii)(c) unchanged (D3av-11) |
| D3at-3 | fixed | the "stranded" branch is misdescribed (D3av-5) |
| D3at-4 | partly fixed | tolerance done; the ε-scan runs along one μ only; A′ not scanned |
| D3at-5 | cited | the reconciliation with D‴.51(iii) is missing (D3av-13) |
| D3at-6 | fixed | count non-discriminating (D3av-9) |
| D3at-7 | fixed | |
| D3at-8 | cause fixed | the new 788 argument is invalid (D3av-7) |
| D3at-9 | fixed | `remaining` is a heuristic; branches untested |
| D3at-10 | fixed | rate attribution and (H3) (D3av-10) |
| D3at-11 | reworded | wrong reason given (D3av-8) |
| D3at-12 | partly fixed | D3av-12 |
| D3at-13 | mostly fixed | W5's caption still says "mirror coordinate α"; C2's stale level (D3av-6) |

## 3. Issues
| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3av-1 | major | **The fold is not the capture threshold.**<br>- With w free, the fixed-w rest has ∂E/∂w > 0 on the whole toward branch (d = 5, 6, 8σ).<br>- The held runner-up's coefficient decays, on an ε-free inner time, to the fold. Capture ⇔ the hold and the travel end before the unfitted feature's own entry.<br>- At fixed x: x = 4.9 is bound at 3.75·10⁻⁹ and 10⁻⁹ but captured at 10⁻¹⁰; x = 20.6 and 24.5 are captured at 10⁻¹¹ and 10⁻¹² (rtol 10⁻¹² and 10⁻¹³).<br>- C3's x = 5.07 is captured at 1.11 w_f.<br>- So "every runner-up above w_f ⇒ (iv)" is false, and C3's window in x is a fixed-μ artefact. | (v): describe the hold, the drift of w and the race; say w_f is the threshold only in the zero-time limit. (iv): restrict to features with no same-sign feature still unfitted at their entry, within a reach that grows as ε → 0. C3: scan fixed x across ε (10⁻⁸…10⁻¹²) and report the hold time in s. Numerics: drop "captures from 1.1 times the fold" as support. |
| D3av-2 | major | **The leading-order law is wrong.**<br>- The pull term lacks the projection off φ, which the binding term has.<br>- Correct: P = A′π^{n/2}σ^{n−2}(d/σ)²e^{−d²/4σ²}/4, so \|δ\|² = σ²ℓ(d)w_e/(w(w_e+w)), and the transverse pull is 0.<br>- The A′ → 0 limit confirms it to 3·10⁻⁴.<br>- The error in \|δ\| is 2–4%; the derivation is false. | (v): statement and proof; `v_law` and V's caption; numerics: "low by 0.8–16% at w = 0.3". |
| D3av-3 | minor | **The fold's displacement is misstated.**<br>- "Disappears when \|δ\| ≈ 0.3–0.4σ": the minimum is at 0.94–0.76σ at the fold.<br>- At the coefficients that decide capture (1.2–2 w_f) the law is 25–40% low. | (v); V: report the minimum's \|δ\| at the fold. |
| D3av-4 | minor | **c depends on A′/A**: 3.14–4.50 at 5σ, 6.18–6.73 at 6σ (A′/A = 3.5–0.44). | (v): write c(d, A′/A) with this table; V: scan A′. |
| D3av-5 | minor | **"Stranded on the way, coefficient left to decay" is wrong.**<br>- The tail is recalled toward its own feature (x = 17.3, run to 5.1 τ_e).<br>- "Stranded" rests on the horizon and on an arbitrary, untested rms cut of 10⁻⁸.<br>- C2's ε = 10⁻¹⁰ residual cannot establish stranding. | (v), C3's caption, `c3_outcome`: classify by the tail's distance at growing horizons. |
| D3av-6 | minor | **Wrong md columns.**<br>- C3's 5σ rows print the captured kernel's 0.7 as "the runner-up after the entry".<br>- C2 still prints round 42's refuted "unbinding level" 0.00041, which the captured tail (0.00106) exceeds. | `c2_job`: record the tail's coefficient before it leaves; `report()`: replace with c(d)ℓ(d). |
| D3av-7 | minor | **The 788-doubling argument is invalid** while the drifts grow (×2 per doubling up to the cap). The counts survive by continuation (peak at doublings 11–13; post-cap change ≤ 1.7·10⁻⁶). `remaining` is an extrapolation. | Numerics and `b4u_note`: cite the continuation, or run all unsettled problems to 24 doublings; call `remaining` an estimate. |
| D3av-8 | minor | **Part D's reason is wrong.**<br>- The finite-ε shift (+0.4 at 10⁻⁵, +0.13 at 10⁻⁶) exceeds the face variation (0.27).<br>- The loop sides with the blind rule 7 of 7.<br>- At 10⁻⁸ the loop resolves x_θ(t). | Numerics (D), `report()`; rerun D, or face-stratified positions, at 10⁻⁸. |
| D3av-9 | minor | **T3.**<br>- The count does not discriminate: H/2 also gets 6 of 7.<br>- The product underestimates every far corner (+1.5–7% at gaps 5–6): the miss is systematic, not "a share at θ". | (iii)(d), (iv), numerics; add triples with corners near θ. |
| D3av-10 | minor | **(i)'s citations.** It cites D‴.53(a)–(b) for "the rate ε^{1/2}", which is (c)–(d)'s; it omits (H3), which for k ≥ 2 is only numerical in ℝⁿ. | (i) and the status. |
| D3av-11 | minor | **Log and draft disagree.** The log's D3at-2 fix "(ii)(c): separation at the entries is not enough" is absent from the draft. | (ii)(c), or correct the log. |
| D3av-12 | minor | **Tests miss the new code.** 18 of 26 mutants survive; 8 of them change the md:<br>- `b4u_note` untested;<br>- the rms cut;<br>- `v_law`, C3's fold column and V's ℓ;<br>- `inner_k`'s clock between chunks and zero drift. | Tests: `b4u_note` on a synthetic B4u; `c3_outcome` at the rms boundary; `v_law` against the A′ → 0 limit; the C3 and V report formulas on a synthetic raw; `inner_k` with max_secs > 0 and a zero drift; `formed_by_feature` with a small kernel that changes feature. |
| D3av-13 | minor | **"Its scale D‴.51(iii)'s" is asserted, never compared**: no comparison of C3's edges with D‴.51(iii)'s migration scale; D3at-5's request is unmet. | Status, (v). |
| D3av-14 | cosmetic | Wording and notation:<br>- c(k) (codimension) and c(d) (fold) clash with c and c_i;<br>- W5's caption still says "mirror coordinate α";<br>- "depends on δ only through M₂" holds only along e;<br>- "its own feature keeps one kernel fewer" holds only for formed runner-ups (C2 at 10⁻⁹: 1 = 1);<br>- the aliasing bound doubles for oblique (hexagonal) neighbours. | Statement, md. |

## 4. Lesson
A static fold computed with one variable held fixed is a threshold only if that variable is slaved. Here the coefficient's own gradient
is nonzero at every held state. The "rest" therefore drifts to the fold on an inner time that does not depend on ε, and the outcome is a
race whose clock runs ε^{−1/2} faster in τ. A scan along one line of the parameter plane (fixed μ) cannot tell a threshold in x from such
a race; holding x fixed and varying ε can, and two runs at x = 4.9 decided it. For the second round running, two balanced terms were
expanded under different constraints (the binding projected, the pull not). When the coefficients are slaved, every term must be
projected. The cheapest test is the limit in which the leading order must be exact (A′ → 0): there the draft's law missed by a fixed
factor of 1.03–1.04.
