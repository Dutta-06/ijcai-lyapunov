# REVIEW_roundTH: referee report on THEORY_NOTE.md (coupled entry layer, Prop. A, Theorem K, network invariance)

Scope. I read the note line by line against the paper's theory sections. I also read exp_theory_h3.py, theory_h3.md and the path code `sliding_iss_nd`. I re-implemented (L), (C), (N) in ball arithmetic (python-flint 0.9.0, repo venv) without importing the author's functions. I then ran my own simulations. Work files are in /tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_TH_work/. The main ones are mycond.py (ball arithmetic), faces.py (face signs), hesstest.py (Hessian), sim.py, below.py, adv.py (flows), gram.py, edge.py and edge2.py (window geometry), and misc.py (Ω, d_0).

**Bottom line.** The mathematics of Lemma 1.1, Prop. A and Theorem K is correct. This holds for the continuum, translation-invariant Gram ⟨φ_cφ_c'⟩ = ZΓ(c−c'), for any n, noiseless. I found no sign error, no face or corner gap, and no circularity. The numbers reproduce. The proof is also conservative by 1–3σ. The problems are in what is claimed and about whom:

- The "certified" counts for the experiment are not covered by the theorem.
- The suggested replacement sentences over-claim.
- The network remark has unsupported constants.
- The path-versus-certificate agreement is not an independent test of Prop. A.

No blocker for the mathematics. Finding 1 is a blocker for the text if the Section 8 sentences are pasted unqualified.

## Findings

**1. [blocker for the text / major] The theorem is for the continuum Gram, but 300 of the 568 "certified" targets live where that Gram is false.**

- Location: Sec. 0 (lines 13, 16), Sec. 4C, Sec. 8 ("covers all 360 …", "120 of 120 three-feature").
- Evidence:
  - The experiment uses 61 samples on [−1,1], σ = 0.12. `certify` computes (a) and (b) and the block in the sample inner product. It checks (c) with the continuum constants (`box_exists(K,R,sep,d0)`).
  - Part A (window ±0.7) is fine: the Gram is translation-invariant to 2e-4 apart from the constant factor 60/61.
  - Part C and the K=3 targets at 6 and 7σ draw centres in ±0.95. By `edge.py`, a feature lies at |c| > 0.8 (< 1.7σ from the edge) in 52–63% of Part C targets at 10 and 12σ and in 100% at 14σ. For K=3 at 6σ the fraction is 95%, and at 7σ it is 100%. A feature lies at |c| > 0.9 in up to 80%.
  - `gram.py`: ⟨φ_c²⟩/Z falls from 0.984 (interior) to 0.902 at c = 0.9 and 0.773 at 0.95. ⟨φ_c, ∂_cφ_c⟩/Z, which is 0 in the continuum, is −0.89 and −1.70 per unit length at c = 0.9 and 0.95. In σ units that is −0.11 and −0.20.
  - Prop. A uses ⟨φ,∂φ⟩ = 0 and ⟨φ²⟩ = Z in Steps 0, 2 and 4. The extra term has the size of the (C) margin (ψ(r) − E_c ≈ 0.09 at Δ = 4.64).
  - In `edge2.py` the sampled centre-face field for a single feature at c* = 0.95 still points inward, but its margin is 0.098 against 0.196 in the interior. So the claim is probably true at the edge, but it is not covered by the proof.
  - Count: 268 of the 568 certified targets are interior (Part A 237, K=3 Part B 31). The other 300 are edge targets (Part C 180, K=3 at 6 and 7σ 120).
- Cheapest fix:
  - Say "the continuum conditions (c) hold and (a),(b) hold in the samples", never "the theorem covers N targets". Or restrict the count to targets with every feature ≥ 3σ from the window edge, or add a window-correction term to E_c and λ_own.
  - The mixed continuum/sample test must be labelled as such in `certify`'s docstring and in theory_h3.md.

**2. [major] Section 8 and Sec. 0 replacement sentences over-claim.**

- "For features at least about 5σ apart the path itself is proved to be the recursion":
  - 5σ is below the computed thresholds: 5.49 (K=2, R=2.4), 5.07 (K=3, R=1), 5.83 (K=3, R=2.4), 6.7 (R=10), and the Part D maximum is 5.83. At 5σ the certificate holds on 57/60 targets, and the Part D median is 5.10σ. So "about 5σ" is false for roughly half of the local data.
  - The sentence also drops these hypotheses: continuum Gram, entrant offset d_0 ≤ 0.42σ, genericity (a),(b), equal walls ρ = 1, and noiseless data.
- "4.3–5.8σ over the local data": this is the range over 150 draws of a certificate computed with the continuum inner product, K=2, R ≤ 2.5. It is not a bound over local data. The (a),(b)-threshold is not uniform: it grows like √ln(1/margin).
- "Δ_0 is computed": what is computed is a sufficient Δ_entry. My simulations show the dynamics are benign far below it (finding 8). It is not a transition separation.
- "Theorem 2" is silently replaced by Theorem K, which has a different conclusion (the recursion, not the decoupled entries). The paper's Theorem 2 text (best kernel per feature, in the order of T_k) is still not proved at any Δ, and the note says so only in Sec. 0.
- Proposed wording is in the last section.

**3. [major] The "568/568" and "1019/1020" agreement is not an independent test of Prop. A or Theorem K's (c).**

- `sliding_iss_nd` (rich_transit_nd.py ~L270–295) runs the flow for t_fast = 1e5 and then replaces it by a least-squares polish. The polish is accepted if the seed units land within 1 width of their start or ascent target.
- The certificate's block is a Newton root started at the features (`Inner.block`).
- Both therefore pick "the critical point near the features". Agreement tests the dual recursion (a),(b) and the arg-min. It does not test that the entry flow converges to that point.
- Convergence is tested only in Part B: continuum, 2 units, one mobility law, 108 runs. The 12 failures there are cases where the first unit has no critical point near the weak feature, so the premise fails. The note handles this honestly.
- Fix: say so in Sec. 4C. Optionally rerun a few targets with `polish=False`.

**4. [major] Remark 2.2 (network invariance): unsupported constants, inconsistent numbers, never tested.**

- The face argument itself is valid. Lemma 1.1 holds for the ε-flow, because W^ε is increasing and η̇ = g, so ẇ = −W^ε′ ∂_wF. The bounds on P_i and C_i need only |Σ_dormant w_j| ≤ Ω, with Γ ≤ 1 and |Γ′| ≤ 0.4289. Invariance needs only (L) and (C).
- "costs nothing in ε (it is much weaker than the margin condition)" is unsubstantiated.
  - Take the core bound W = a_max(εb)^{3/2}/(m√(2v_s)) · cosθ/sin³θ. With m = 12, v_s = 0.025, a_max = 1, b = 1 this gives |w_j| ≈ 0.096 ε^{3/2} μ^{−3}.
  - With about 20 dormant kernels, Ω ≤ Ω_max ⇔ μ ≳ (20 · 0.096/Ω_max)^{1/3} ε^{1/2}. That is about 2 ε^{1/2} for Ω_max = 0.26, a_min = 1, and about 2.6 ε^{1/2} for a_min = 0.4.
  - This is the same ε^{1/2} scaling with a constant comparable to ξ_θ ≈ 2 in the margin condition. Example: ε = 1e-5, μ = 0.005 gives Σ|w| ≈ 0.49, which exceeds Ω_max.
  - The constant C of (M2) is not carried out, and Sec. 5 shows the mismatches occur exactly at small margins.
- Remark 2.2 says Ω_max ≈ 0.05 at Δ = 6, R = 1. Sec. 0 and the report say 0.26. I get 0.2608. 0.05 is stale.
- Because only (L) and (C) are needed for invariance, the correct Ω_max is much larger: my values at (R, Δ) = (1, 6), (2.4, 6), (10, 7) are 0.69, 0.58, 0.71. The note uses the (N)-limited values 0.26, 0.12, 0.037, which are needlessly restrictive.
- No run anywhere checks box invariance on the network.

**5. [minor–moderate] The entrant "at ξ = +∞" versus v_0 ∈ (0, v_lo] is not argued (Thm K proof, Prop. A statement).**

- M(v) ~ v^{4/3} makes v = 0 an equilibrium line. The entering orbit leaves it only asymptotically, so "v_0 ↓ 0" needs justification.
- What is needed: the actual orbit reaches v_e = v_0 inside B with |d_e| < r. That means r > d_0 strictly, plus a far-field displacement bound for d_e and for the block. The block sits inside D*, which has slack in B (about 0.03 in v and 0.2 in d at Δ = 4.64).
- Observation in the note's favour: the proof of (i) uses only "initial point ∈ B". The hypothesis "block at a critical point of its own fit" in Prop. A is superfluous. Say so, because it makes the far-field step easy.
- Theorem K (c) says "r ≥ d_0". Prop. A is fine with equality, but the passage needs strict inequality, or a sentence on continuity of the entering orbit.

**6. [minor] d_0 = 0.41 is not the benchmark's worst case.**

- h/(2σ) = 0.1/0.24 = 0.4167. The nearest-kernel offset can reach 0.4167. The paper's own Sec. 4A text says "0.41σ at h = 0.83σ", which is incorrect rounding down.
- At d_0 = 0.4167 my thresholds are 4.642, 5.490, 5.068 and 5.831 for (K, R) = (2, 1), (2, 2.4), (3, 1), (3, 2.4). The 4.64 and 5.48 figures hold only for 0.41.
- In 391 of 1154 paths the entrant is not the nearest kernel, so d_0 can be larger. Theorem K covers this by a per-target check. State it.
- Use d_0 = 0.42 in the headline numbers.

**7. [minor] (H4) is covered only up to v_lo ≈ 0.90–0.96 a_min.**

- The paper's (H4) concerns all levels below the final coefficient (≈ R a_min for a strong entrant). The note says "levels ≤ v_lo", which is honest.
- Add one line: this suffices if the formation threshold θ is below v_lo a_min.

**8. [minor] Sharpness and "just below the threshold".**

- Simulation (`below.py`, `sim.py`) of the 2-unit flow with the critical-curve-type mobility, unit 1 at its critical point and the entrant at d_0 = 0.4167σ:
  - For R = 1 with equal or opposite signs, and for (2.4, ±1) and (10, 1) with the strong feature first, the flow ends at the target (error ≤ 1e-9) for every Δ down to 1.5σ.
  - The only failure is a premise failure with the weak feature first: (1, 2.4) at 2.5σ, where unit 1 has no critical point near its feature and ends 11σ away. Part B reports the same.
  - There is no sharp transition at the thresholds. The inequalities are conservative by about 1–3σ.
- Do not describe Δ_entry as where the structure appears or breaks.
- The tables are conservative at d_0 = 0, because the author's r-grid is coarse. For d_0 = 0 the table exceeds my threshold by up to 0.078σ. All 50 entries are ≥ my value (see Verified).

**9. [minor] Theorem K assumptions that "nothing depends on…" hides.**

- "Nothing depends on mobility, mass, walls, slide length" is true for Prop. A: the proof uses only positivity of M_i and μ_c. It is not true for Theorem K.
  - K needs equal walls (the recursion Δt is written for ρ = 1).
  - (H2)'s λ_k = (λ_own − H_c) min(M, μ_c) depends on the mobility.
  - Δ_entry should be the max over the entries' d_0^k, not a single d_0.
- Dimension n > 1: I checked the proof by hand and it carries over. The radial 2×2 reduction of the own block works because the tangential eigenvalue is larger. Schur block norms give the same ψ and |Γ″| bounds for ℓ ≥ √6, where radial is greater than tangential. But it has not been verified numerically, and no planar certificate has been run. Keep it labelled "same proof, unchecked".

**10. [minor] Internal inconsistencies and loose wording.**

- Sec. 4C says times "≤ 1e-13 except one 2σ target". The table has 2.4e-13 (K=2, 3σ).
- "31/180 at 3–5.8σ" is really min-separation ≥ 3, 4, 5.
- "Main experiment" includes Part C, which was added after the review.
- "Δ_0 grows like 2σ√ln R" is only qualitative. The increments are 0.64, 1.68, 3.2, 4.36 versus 1.87, 3.03, 4.29, 5.26.
- Corollary A's "δ_0 is the largest δ with (L), (C), (N)" is loose, because the box and δ_0 = Γ(Δ − 2r) are tied by r.
- Part B's "in box" test divides w_min by a[0] rather than a_min and does not test the entrant's lower bound. This is harmless (conservative) but should be corrected.
- Sec. 5: the slope of the largest mismatching margin is an extreme-value statistic of 200 targets. Its 0.28 should not be read as a scaling exponent. The note's suggested sentence ("falls as ε^0.28 over these four values") is acceptable only if it is called descriptive.

## Verified (independently re-derived)

- **Algebra (Step 0 and Step 1).**
  - F(v,d) expansion, −∂_vF = aΓ(d) − v + εP, and −∂_cF = −v[aψ(d) − εC]. The sign convention Γ′ = −ψ is correct.
  - Own Hessian entries 1, s = aψ(d), t = va(½ − d²/4)Γ(d).
  - Monotonicity: ψ increases on [0,√2] and decreases beyond; |Γ″| decreases beyond √6; (½ − x²/4)Γ decreases on (0,√6); sup|Γ′| = 0.4289 < 0.43.
  - The v̄ formula makes R + E_v = v̄ at the "1+0⁺" limit, as required.
- **Faces and corners.**
  - Every face, including corners, is handled by uniform bounds, so the Nagumo / Bony–Brezis argument is valid.
  - Numerically, `faces.py` (15,000 random points per case on the surface and corners, including opposite signs and random spacings) found all inward margins > 0 for (K,R,Δ) = (2,1,4.64), (2,1,4.70), (2,2.4,5.48), (3,1,5.06), (3,2.4,5.82), (2,10,6.68). The smallest was 0.0019 on the lower v face at R = 10.
  - Entrant lower face at v_0: margin ≈ 0.48.
- **Hessian (N).** `hesstest.py` sampled D* for the same cases, 4,000 points each. The true minimum eigenvalue is 0.39–0.50 against the certified lower bound λ_own − H_c of 0.0026–0.0095. The bound is valid and very loose.
- **Box invariance with adversarial mobilities.** `adv.py`: 128 runs (8 configurations × 16) from random starting points on the faces and corners of B. Per-unit scales 0.1–10, centre mobility 0.1–10, M = const, v, v^{4/3}, random amplitudes and signs. The (K, R, Δ, units) cases are (2, 1, 4.64, [0,1] and [1]), (2, 2.4, 5.48, [0,1] and [1]), (3, 1, 5.06, [0,1,2] and [0,2]), (3, 2.4, 5.82, [0,1,2]) and (2, 10, 6.68, [0,1]). Result: 0 exits from the box, the entrant monotone below v_lo in all, and every run converged to the target.
- **Thresholds.**
  - Δ_entry(K=2, R=1, d_0=0.41) = 4.629 by bisection, so the table's 4.64 is the 0.02-grid ceiling.
  - At Δ = 4.64 the feasible radius window is r ∈ [0.410, 0.4155]. Γ(Δ − 2r) = 0.0260–0.0262 ≤ 0.0265. The margin (N) is 0.0095.
  - (L),(C) alone: 4.232.
  - Other entries: 5.477 (R = 2.4); K=3: 5.055 and 5.818; R = 10: 6.678.
  - Whole tables (K = 2, 3; R = 1…1000; d_0 = 0…0.8): all 50 entries are ≥ my values (no unsound entry). The maximum excess is 0.078σ, at d_0 = 0 (coarse r-grid).
  - Existence of the box is monotone in Δ on [3.5, 14] for (K,R) = (2, 1) and (3, 2.4) (a single flip).
- **Ω_max.** 0.2608 (R = 1, Δ = 6), 0.1245 (R = 2.4, Δ = 6), 0.1158 (Δ = 5), 0.0372 (R = 10, Δ = 7). These match the report and not Remark 2.2's 0.05.
- **Coverage counts.** The counts in Sec. 0 are consistent with theory_h3.md: 57 + 180 = 237 of 540, plus 180 for Part C; 1 + 4 + 26 = 31 of 180, plus 120; total 568; 540 + 180 + 180 + 120 = 1020. Part B's 8 + 16 + 18 = 42/42.
- **Flow checks.** Equal, opposite and unequal signs, R ∈ {1, 2.4, 10}, d_0 = ±0.41 and Δ ∈ {3.5, 4.3, 4.7, 5.5, 8}: for the cases where unit 1 has a critical point near its feature (all strong-first and equal-amplitude cases; the exceptions are in finding 8), the flow stays within |d| ≤ 0.41, the entrant rises monotonically, and the flow ends at the target to ~2e-11. This includes Δ = 4.3 (only (L),(C) hold) and Δ = 3.5 (no box).
- **Logic.** Induction in Theorem K is well-founded: the block at step k is the previous layer's output; (b) is implied by (a) plus convexity of |u| on piecewise-linear pieces; the sign condition is in (a). The only non-trivial passage is finding 5.

## Could not check

- The path-vs-certificate agreement itself (1019/1020, 568/568): I did not rerun the repository's path (cost). I only analysed the code (finding 3).
- Parts D–F numbers (150-draw thresholds, ε_max table, mismatch margins). Internal consistency only, plus the arithmetic of the slopes and the units of ε^{1/2}.
- Lemma 1.1's profile facts for M (W′ > 0, M ~ v^{4/3} at 0⁺, local Lipschitz). I took the critical-curve formulas from the supplement; Lemma S1(ii) itself is outside this note.
- Theorem 1 itself (M1–M4); the claim for n > 1 numerically; the network at ε > 0 (no run); the planar certificate; the decoupled-prediction threshold Δ_dec.
- Whether the sampled-Gram edge effect actually breaks the face conditions for the Part C and K=3 targets. I showed only that the proof does not cover them and that one single-feature edge case keeps its sign with half the margin.

## What the paper may claim, with wording

**Proposition S-A (proved).**

> Let Γ(d) = e^{−|d|²/4σ²} and let the inner product be ⟨φ(c)φ(c′)⟩ = ZΓ(c−c′) on ℝⁿ (continuum, noiseless). Let f* be K Gaussian features of the kernels' width, pairwise ≥ Δ apart, with amplitude ratio R. Let a block of units sit at their fit's critical point and let an entrant start at offset d_0 ≤ r. Suppose conditions (L), (C), (N) hold for some r ≥ d_0 (e.g. Δ ≥ 4.64σ for K=2, R=1, d_0=0.41σ; ≥5.49σ for R=2.4, d_0=0.4167σ). Then for every positive coefficient mobility and every positive centre mobility the box is positively invariant, the flow converges to the unique critical point in the box (nondegenerate), and the entrant's coefficient rises strictly through all levels ≤ v_lo. With all K features fitted the limit is the target.

State the n = 1 numerics and call the "any n" statement "same proof, checked by hand".

**Theorem K (proved, conditional).**

> For the limit object of Definition 1 with equal walls, in the continuum Gram, if (a), (b) (finite generic checks) hold and Δ ≥ max_k Δ_entry(K, R, d_0^k), the path is the recursion with exactly K entries, no exit, and the final block is the target. (H1)–(H3) hold, and (H4) holds at levels ≤ v_lo a_min.

**Not claimed.**

- Anything about the network's ε → 0 behaviour (Theorem 1 stays an outline).
- "Δ_0 exceeds every tested separation" is false for Theorem K. Replace by: "In the continuum model the entry-layer condition holds from ≈ 5.5σ (K=2, R ≤ 2.4) and ≈ 5.8σ (K=3); the sampled experiment's targets are only partially covered (268 of 568 certified targets have interior geometry)."
- No part of the theorem covers 2–4σ. The statements at those separations remain numerical.

**Abstract sentence (replace the Section 8 version).**

> "For features separated by at least about 6σ (two features, amplitude ratio ≤ 2.5; three features similar) the limit path is proved, in the continuum, to be the recursion of Theorem K: one unit per feature, exact fit, no unit leaves. The separation is computed (4.6–5.8σ across our benchmark's local data, by a conservative sufficient condition), and the trained network is observed to form one unit per feature on 1154 targets down to 1.5σ."

**Remark 2.2 sentence.**

> "The same box is invariant for the trained network while the dormant coefficients sum to less than Ω_max(Δ, R) (0.69 a_min at Δ = 6, R = 1 for (L),(C) alone). We have not quantified how this condition scales in ε."
