# Referee report D3er, round 149: the Paper 17 manuscript, second review

## Scope
- **What I reviewed.** manuscript/paper17/paper17.tex (31 pp.) and supplement.tex (16 pp.), at commit 7f5d4b3, branch claude/vigilant-clarke-hf939b. Title: "The small-mass transit of self-adaptive kernel networks: the sliding limit and a predictor of recovery". The venue considered is SIAM J. Math. Data Sci. or Found. Comput. Math.
- **Against what.** The first review (review/REVIEW_round140_D3ei.md, 22 issues) and the revision record: the README's "Revision after round 140" and "Runs and shortening after round 140", REVIEW_LOG's Round 140 entry and its "Paper 17: the runs declined after round 140" entry, and PARADIGM.md's block "Corrections after referee round 140".
- **Workspace.** W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref149 (repo export in W/repo; my scripts in W/own). I did not modify the repository.

## Status of the round-140 issues

| D3ei id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 σ_k's 1/ρ_{j_k} | yes | Re-derived in Step 4 [read][own]: s_k(u + ε^{1/2}dev) = ρ with s_k u growing at ρv_k gives σ_k = −s_k dev/(ρv_k) = −s_k dev/\|a_{k−1,j_k}\|. Part U [md]: 0.26581/(0.9601·0.05411) = 5.116 against the printed 5.1162; without the ρ, 4.912 against 4.9121. The crossings fall like ε^{1/2} (0.1095 → 0.03466 → 0.01097, ratio 3.16 per decade). Without the factor they do not fall (0.1931). `test_unequal_masses_the_shift_needs_one_over_rho` exists and passes [own]. | — |
| 2 Thm 5.3(b), left end | yes | New Step 2 [read]. It fixes s₂ and a radius r, puts the block within r of S_k at s₂ (Step 3 and (H3)), and uses V = (y−S)ᵀH(y−S), which decays at exactly 2λ_k because M^{1/2}HM^{1/2} ≥ λ_k ⇔ H⁻¹ ≤ M/λ_k [own]. The O(ε^{1/4}) of the old Step 5 is absorbed into "within r", so this is correct. Part W [md]: slopes 1.424–1.597 at C_k = 3/(2λ_k). The old C = 30 start follows ε^{30λ_k} (0.797 against 0.794). | Two technical points: D3er-3 (the window's end lies outside J_k^ε) and D3er-4 (C₁ = Ce^{λ_k s₂}r). |
| 3 (H3), the identity t = Λ₀s | yes | Lemma 5.2(i′) re-derived [own]: η_c′ = coth²/(2Gς) = 1/(Gβ_c a′), so dα/ds = Λ₀Gβ_c a′⟨rφ⟩. (H3) is now rephrased as the relaxation's convergence, and "another metric" is gone. (iii) is proved in the continuum, with far-field control in which ξ_e serves as time; the bootstrap and the IFT/LaSalle end are sound [own]. | D3er-5: the splitting F_B/F_× is inconsistent with S_{k−1} (cosmetic to the argument). |
| 4 Prop. 3.4's sharp correction | yes | Checked by exact 400-bit (python-flint arb) bisection of η = (1−μ)λ₁ [own]. (rel. error − stated correction)·μ⁴/ε² is bounded: 0.45–0.557 for b = 1 and 2.8–3.48 for b = 2.5, at ε = 10⁻⁶…10⁻¹², μ = 0.3…0.003. So the form −(εb/ς)cot²(πμ/2) − εb/(4ς) + O(ε²μ⁻⁴) is confirmed. The supplement's two expansion lines are re-derived, and the dual-margin term (πεb/(8ς))(1−μ)(tan + 3cot) checks. With the runs' constants (ς = 0.025, ε = 10⁻⁸) the exact ratio is 0.98956, 0.99147, 0.99465, 0.99613; the md gives 0.9896, 0.9915, 0.9946, 0.9961. | D3er-8 (wording). |
| 5 the headline | yes | Title, abstract, Table 5 caption and §7.5 now present 1154/1154 as a check of the basin hypothesis. Part G covers all 1154 targets. D(c)'s selection is disclosed. Baselines: η_V, basis pursuit, continuous OMP, BLASSO (Part I), the plane (Part J). | D3er-6 and D3er-7 (statements about I and J). D3er-11: Part G is stronger than stated. |
| 6 literature | mostly | EGSH21, SRR18, JGSHG21, TA20 added. Chizat is compared on its main result. BPV Theorems 3.8/4.2 are used. GSD20 is corrected. | AGF is under-used: D3er-1. |
| 7 self-containedness and length | yes for content, no for placement | Motivation via Lemma 2.1(o); the companion paper is no longer needed (Lemma 2.1(v) and §S6.3). Figure 1 added. 31 pp. | The main theorems' proofs now live only in the supplement: D3er-2. |
| 8 Prop. 4.2 statuses | yes | (b): the symmetric case re-derived (⟨r∂_cφ_e⟩ = (Zd/2σ²)(−AΓ(d) + 2wΓ(2d)) forces Γ(2d) = 1). (c): stated as the exact ratio with the Gibbs limit. (d): s_A added, and I∞ > 0 marked numerical. Table 1 matches. | — |
| 9 dual potentials add | yes | §S6.3 [read]: the indicators sum to the box; the primal side is an infimal convolution. The branch bounds are re-derived. | — |
| 10 Prop. 6.3(d) | yes | Now reads r_k = f* − Σ_fitted A_mφ_m. | — |
| 11 Lemma 2.1(iv) | yes | Re-derived from arctan y = π/2 − 1/y + O(y⁻³). | — |
| 12 notation | largely | Table 2 rewritten; A_k and 𝒜_k now separate. | D3er-12 (q, p, δ, λ, η still overloaded). |
| 13 dimension | yes | §6 opening. | — |
| 14 §7.2 pairs, staircase | yes | §S7.2. | — |
| 15 Remarks 6.5, 6.6 | yes | x_θ, ξ_θ and the midpoint inequality re-derived. | — |
| 16 bisection resolution, drifting runs | yes | Table S1 caption; §7.3 reports 51 of 54. | — |
| 17 energy identity | yes | §S6.2: limsup ∫\|ζ̇^ε\|² ≤ −D_∞(ζ(T)) = ∫\|ζ̇\|², so the convergence is strong [own]. | — |
| 18 μ³ vs μ⁴; FS 3.5 symmetry | yes | Paragraph after FS 3.5. | — |
| 19 Prop. 3.1 | yes | u̇ = h·3u²g re-derived; the scales marked heuristic. | — |
| 20 Step 3's term | sources named | The O(ε) remainder is listed as open problem (2). | Acceptable. |
| 21 margins, C = 30 | yes | §7.4, Part W. | — |
| 22 bibliography | partly | 30 % CHECK items remain (out of scope). | The AGF venue can now be cleared: D3er-13. |

## New findings

### Major

**D3er-1. Literature and novelty: AGF.** The sliding inverse scale space is an instance of AGF, and Theorem 5.3 proves a case of AGF's open conjecture (§1.3, p. 4, l. 125; claim (ii), l. 139).

I read Kunin et al. (arXiv 2506.06489v4) in full for §§1–3 and 6 [read]. AGF alternates two steps:
- dormant units accumulate a utility S_i = ∫κU_i, with U_i = ⟨f_i, r⟩, and activate when S_i reaches a threshold c_i ∝ ‖θ_i(0)‖^{2−κ} (κ > 2);
- active units then minimise the loss to a local minimum, and collapsed units return to dormancy.

The paper's unit maps onto this directly:
- It is depth-3, so κ = 3. The threshold ∝ ‖θ(0)‖⁻¹ ∝ (εb_i)^{−1/2} is exactly the wall Λρ_i, with ρ_i = b_i^{−1/2}.
- The accumulated utility of a kernel with a frozen centre is ∫⟨rφ(·;c_i)⟩, i.e. the dual u_i.
- The cost step is the block relaxation of Definition 2.3, and exits correspond to AGF's collapse.
- The one difference is AGF's directional step (dormant units rotate to maximise utility). Under (2.2) it is switched off, because the centre's drift carries the mass.

So Definition 2.3 is AGF's Algorithm 1 for this unit with the utility step frozen. AGF's own §6 states that whether the AGF ↔ gradient-flow correspondence holds as the initialisation vanishes "remains open" beyond diagonal linear networks. Theorem 5.3 therefore proves that correspondence for a nonlinear model with moving positions, with the rate ε^{1/2} and the first-order correction.

This is both a fairer and a stronger positioning than the present "closest precedent of 'slaved to the fit'". As written, claim (ii) ("the sliding inverse scale space as the small-mass limit") reads as if the limit object were new; it is AGF's.

**D3er-2. Placement: the main theorems' proofs are only in the supplement** (§1.1 l. 89, §§3–6 "Proof. Supplement"). The shortening moved every proof except Lemma 2.1 and Proposition 6.2, including those of Theorems 3.2 and 5.3 and Lemma 5.2. The main text keeps outlines.

The moved proofs are complete, the xr-hyper references resolve (rebuild: 0 undefined, no "??", text identical to the shipped PDFs), and each outline is accurate. But for FoCM, and arguably for SIMODS, the proof of the central theorem must be refereed as part of the paper.

Proposed fix:
- Return §S6.2 (Theorem 3.2), §S6.8 (Lemma 5.2) and §S6.9 (Theorem 5.3) to an appendix of the paper.
- Move instead secondary numerics to the supplement: §7.2's creep paragraph, §7.5's BLASSO/plane paragraph, or Table 3/4.

### Minor

**D3er-3. Theorem 5.3(b): the window's right end lies outside the interval of the proof** (p. 18, (b); §S6.9 Step 2).
- J_k^ε ends at τ_{k+1} − ε^{1/2}X/v_{k+1} = τ_{k+1} − ε^{1/4}/v_{k+1}.
- (b)'s window ends at τ_{k+1} − ε^{1/4}, which is later whenever v_{k+1} < 1. That is every entry computed (v_k ∈ 0.00132–0.125).
- The Lyapunov argument extends there, since |p| ≤ Cε^{3/4}v⁻³ stays small, but as written (b) is not inside J_k^ε.
- Fix: end the window at τ_{k+1} − ε^{1/4}/v_{k+1}, or say Step 2 holds wherever v_{k+1}(τ_{k+1} − τ) ≥ ε^{1/4}.

**D3er-4. Theorem 5.3(b): the constant C₁ hides the layer's position** (Step 2, "C₁ = Ce^{λ_k s₂}r").
- s₂ must exceed the translated layer, s^θ_k + σ_k + O(1/λ_k). Over the 48 targets |s^θ + σ| reaches 9.856·10⁴ with λ_k ≈ 0.02, so C₁ ≈ e^{2000}.
- The statement is correct as an asymptotic one, but the window start τ_k + Cε^{1/2}log(1/ε) is informative only once log(1/ε) ≫ λ_k s₂/C.
- On Part W's eight targets this is harmless, since C_k log(1/ε) ≈ 300 > s^θ ≈ 70–130 at 10⁻⁶.
- Fix: state the start as τ_k + ε^{1/2}(s^θ_k + σ_k)₊ + Cε^{1/2}log(1/ε), or add a sentence to Remark 5.4.

**D3er-5. Lemma 5.2(iii) proof: the splitting does not match S_{k−1}** (§S6.8).
- F_B is "the block's fit of the residual without the entrant's feature", but S_{k−1} is critical for the full F at w_e = 0, i.e. for F_B + F_×|_{w_e=0}, not for F_B.
- "V = F_B − F_B(S_{k−1})" and "F_B ≥ F_B(S_{k−1})" should use F_B + F_×(·, w_e = 0). Equivalently, put the block's static cross term with the entrant's feature into F_B; then the w_e-dependent cross term is O(δ|w_e|), which is what the far-field estimate O(δξ⁻³) uses.
- The conclusion is unaffected.

**D3er-6. Part I's failure description is wrong in two places** (§7.5 p. 28, "Its failures at the larger ℓ are features lost below the regularisation, in the wide amplitude ranges of B and D(a)"; §S7.5 "a feature lost … accounts for every failure at ℓ = 10⁻¹ but the 20 targets of D(b)").
- I recounted from the stored solutions [own]. At ℓ = 10⁻¹ there are 366 failures:
  - 306 have a lost feature, spread over A 37, B 66, D(a) 39, D(b) 80, D(c) 84;
  - 60 have every feature present but an atom 0.01–0.012 from it, beyond the two-spacing tolerance: A 26, B 7, D(a) 3, D(b) 24. Example: an atom at 0.405 for a feature at 0.393. This is the BLASSO's O(λ) support shift, not a loss.
- The supplement's sentence is true at ℓ = 10⁻⁴ (22 failures = 2 D(a) losses + the 20 η_V-degenerate D(b) targets). "10⁻¹" there looks like a typo for 10⁻⁴.
- Fix: correct both sentences, and tabulate the failures by type (lost / displaced / spurious) against the 0.01 tolerance.

**D3er-7. Part J: "below one exactly where the full one is" is vacuous** (§7.5 p. 28; §S7.5).
- On all 1034 one-dimensional targets both t (≤ 0.9843) and t_J (≤ 0.99952) are below 1. So the equivalence was never tested on a failing case.
- Say "both below one on all 1034; t_J − t ∈ [0.0152, 0.137]".
- The restriction itself is valid, and t_J ≥ t holds (D3er's verification below).
- Also say plainly that basis pursuit "with the features added" is the oracle on-grid problem.

**D3er-8. §7.2 creep paragraph: two wording points** (p. 23).
- "the size of the next order (ε/μ²)²" is dimensionally (εb/(ςμ²))². With ς = 0.025 that is ≈ 1.1·10⁻⁴, whereas (ε/μ²)² ≈ 4·10⁻⁷.
- The ratio 0.5515 at μ = 4.49·10⁻⁴, described as "not small, kernel entering", is reproduced by the exact closed form of §S6.4 (0.5513 [own]). One sentence saying the exact core formula still holds there would strengthen the paragraph.

**D3er-9. Theorem 5.3, Step 3 is a careful sketch rather than an estimate** (§S6.9).
- The claim that sources (2) and (3) "act, to this order, along the translation in s" is argued from autonomy and the block's equilibrium, not bounded.
- Table 1's "proved under (H1)–(H4)" is defensible, but a FoCM referee will want Step 3 written as a Gronwall estimate in the variational equation along the inner solution.

**D3er-10. README inconsistencies** (manuscript/paper17/README.md).
- l. 9 says "PARADIGM.md was not changed", but PARADIGM.md l. 6549 carries the dated corrections block.
- The "Numbers" section gives two incompatible check results: 180/156 tokens with supplement 168 (current; my run agrees) against 252/221 with supplement 71 (stale).

**D3er-11. Part G understates its own result** (§7.5 p. 26).
- The stored entry sequences of the unpolished flow equal the polished path's on all 1034 one-dimensional targets (640 in A/D(b) + 394 in G) [own]. In the plane, the outcome class agrees on 120 of 120.
- "Its outcome is the polished path's" can be strengthened to "the same entries in the same order (one dimension)".

**D3er-12. Remaining notation overloads.**
- q: sample index (§2.1) and cross-talk (§6).
- p: certificate p(c), Step 2's forcing p(τ), and §S7.5's dual vector p.
- δ: 5 meanings (offset §7.3, overlap Lemma 5.2(iii), Thm 3.2 proof, δ_k twice).
- λ: wall, block rate, BLASSO weight.
- η: mirror coordinate, η_V, η₀ in §S6.8.

**D3er-13. Cosmetic: the AGF % CHECK can be cleared.** arXiv v4 of 2506.06489 prints "39th Conference on Neural Information Processing Systems (NeurIPS 2025)" [read].

## Verification I did
- **[own] Exact arithmetic for Prop. 3.4.** W/own/prop34.py and creep_md.py, 400-bit arb. Results as in the table: the sharp correction and its remainder; the md's four creep ratios; 0.5513 at μ = 4.49·10⁻⁴.
- **[own] Derivations.**
  - Lemma 2.1(o), (i), (iv), (v).
  - Prop. 3.1's h.
  - Prop. 3.3's bounds and the dual sum.
  - The constants κ = 0.7025, k₃ = 0.0962, Λ₀ = 7.0248.
  - Prop. 4.3's thresholds.
  - Lemma 5.2(i), (i′).
  - Thm 5.3: Step 2's Lyapunov rate, Step 4's σ_k and the dev recursion.
  - Prop. 4.2(b)'s symmetric case.
  - Prop. 6.2.
  - Remarks 6.5, 6.6 and S3.2 (μ = 1 − e^{−h(h−2|δ|)/(4σ²)}, P ≈ 4σ²y/h²).
  - Part U's σ₃ = 5.116 with the factor, 4.912 without.
- **[own] Reproduction.** transit_certificate.md and transit_limit.md were regenerated with `--from-raw` from copied raw files (TC_OUT, TL_OUT): both byte-identical (`cmp`).
- **[own] Part I.**
  - Independent KKT check of the stored BLASSO solutions on 15 targets: worst violation 2.5·10⁻¹² relative to λ.
  - An independent recount of the recovery rule over all 1034 × 4 solves: 668, 974, 1000, 1012, as printed.
  - The failure-type breakdown in D3er-6.
  - The D(b) degenerate cases fail on spurious atoms of size ±3.5·10⁻⁴ ∝ λ at about ±0.2, consistent with η_V's degeneracy.
- **[own] Part J.**
  - The jet certificate re-solved with my own finite-difference jets in the raw (non-orthonormalised) basis with scipy HiGHS, on 8 one-dimensional and 3 plane targets.
  - t_J matches the stored values to 6 digits. The certificate's η satisfies the equality constraints to ≤ 3·10⁻¹⁰, and its sup over the grid equals t_J.
  - Validity: any p in the jet span with ⟨pφ(c*_k)⟩ = sign A_k and |η| < 1 elsewhere on the grid, together with injective feature columns, certifies unique ℓ¹ recovery on the grid. Restricting p can only raise the minimum, so t_J ≥ t. Observed min(t_J − t) = 0.0152 over 1034.
  - η_V lies in the span, so t_J ≤ sup|η_V|.
- **[own] Part G.** Entry sequences, polished against unpolished: 1034 of 1034 identical (one dimension); class 120 of 120 in the plane.
- **[own] Build, number check, tests.**
  - Build (paper/supplement alternately ×3): 31 + 16 pp., 0 undefined references, text identical to the shipped PDFs.
  - numcheck17: paper 180/156, with only 0.9895 and 0.9914 not in an md; supplement 168/159, 0 not found.
  - pytest tests/test_transit_limit.py and tests/test_transit_certificate.py: 47 passed (148 s, nice, one process).
- **[md]** Parts G, H, I, J, U, W and Table 5 against the text: all match, except the descriptions in D3er-6 and D3er-7.
- **[read]** README; REVIEW_LOG (Round 140 and the runs entry); PARADIGM.md's correction block; every proof in the supplement.

## Literature
- **AGF (read in full: §§1–3, 6).** Its framework covers the sliding path's structure. Its convergence theorem is for diagonal linear networks only, via Pesme–Flammarion, and it calls general correspondence open. See D3er-1. The paper's description "for diagonal linear networks they prove convergence to gradient flow" is fair, since AGF → PF algorithm and PF → gradient flow.
- **The other comparisons are now adequate:**
  - Chizat (regularised, non-vanishing masses, global);
  - EGSH21 (continuous OMP; completely monotone functions of |·|^p);
  - SRR18 (positive spikes without separation);
  - TA20 (basins);
  - JGSHG21; BPV26 (Thms 3.8/4.2, the weighted against unweighted point);
  - DDPS20 and BSR17 (discrete-time sliding analogues);
  - PF23 and Ber23.
- **What would help a SIMODS reader.** Part H's "continuous OMP with a joint slide" is exact on 1154 of 1154; one sentence relating it to DDPS20's sliding Frank–Wolfe (which selects by the certificate's maximum, not by duals) would be useful.

## Overall assessment
**Minor revision.**

The round-140 majors are resolved correctly:
- σ_k now carries 1/ρ_{j_k} and is checked on unequal masses;
- (b)'s left end is now proved with the block's own stable field;
- (H3) is identified as the path's relaxation;
- Prop. 3.4's sharp correction is confirmed to high precision;
- the recovery headline is honest, with real baselines (η_V, basis pursuit, BLASSO, continuous OMP).

The new numerics are valid. The BLASSO solver meets its KKT conditions, and the jet certificate is a valid upper bound that proves exactness on the grid. All reports reproduce byte for byte.

What remains is positioning and presentation:
- D3er-1: the AGF correspondence should be stated. It changes novelty claim (ii), and in the paper's favour.
- D3er-2: the main theorems' proofs belong in the paper, at least as an appendix, especially for FoCM.
- D3er-3 to D3er-9: technical precision and two incorrect descriptive sentences about Part I.

If D3er-2 is done, I see no obstacle for SIMODS. For FoCM, the rigorous core is still Theorem 3.2 plus Theorem 5.3 under (H3), and Step 3 should be tightened (D3er-9).

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3er-1 | major | Def. 2.3 is AGF's Algorithm 1 for κ = 3 with the utility step frozen: threshold ∝ (εb)^{−1/2} = wall Λρ, accumulated utility = u_i. AGF §6 calls the gradient-flow correspondence open. Thm 5.3 proves it here; the text only says "closest precedent" | §1.3 l. 125 and claim (ii) l. 139: state the correspondence; reframe (ii) as a proof of the AGF limit, with rate and correction, for a model with moving positions |
| D3er-2 | major (FoCM) / minor (SIMODS) | Proofs of Thms 3.2, 5.3 and Lemma 5.2 only in the supplement | Restore §§S6.2, S6.8, S6.9 as an appendix; move secondary §7 material to the supplement |
| D3er-3 | minor | (b)'s window ends at τ_{k+1} − ε^{1/4}, but J_k^ε ends at τ_{k+1} − ε^{1/4}/v_{k+1}, earlier on every entry (v < 1) | Thm 5.3(b) or Step 2: adjust the end, or extend Step 2 explicitly |
| D3er-4 | minor | C₁ = Ce^{λ_k s₂}r with s₂ ≳ s^θ + σ (up to ≈ 10⁵): astronomically large on slow entries | (b): window start + ε^{1/2}(s^θ_k + σ_k)₊, or Remark 5.4 |
| D3er-5 | minor | §S6.8: S_{k−1} is critical for F_B + F_×(w_e = 0), not for F_B | Redefine the split |
| D3er-6 | minor | Part I failures misdescribed: at ℓ = 0.1, 60 of 366 failures are displaced atoms (A 26, B 7, D(a) 3, D(b) 24), and losses occur in every family. §S7.5's "ℓ = 10⁻¹" should read 10⁻⁴ | §7.5 p. 28 and §S7.5: correct; table the failure types and the tolerance |
| D3er-7 | minor | Part J's "exactly where the full one is" is vacuous (all 1034 below 1); the grid with features inserted is an oracle grid | §7.5 and §S7.5: rephrase; state the oracle grid |
| D3er-8 | minor | "next order (ε/μ²)²" should be (εb/(ςμ²))²; 0.5515 is the exact closed form's 0.5513 | §7.2 p. 23 |
| D3er-9 | minor | Thm 5.3 Step 3: the "along the translation" step is argued, not estimated | §S6.9: a Gronwall estimate in the variational equation |
| D3er-10 | minor | README l. 9 contradicts PARADIGM.md's correction block; two different number-check results | README |
| D3er-11 | minor | Part G is stronger than stated: identical entry sequences on 1034 of 1034 one-dimensional targets | §7.5 p. 26 |
| D3er-12 | minor | q, p, δ, λ, η still overloaded | Table 2, passim |
| D3er-13 | cosmetic | The AGF venue is confirmed (NeurIPS 2025, printed on arXiv v4) | Bibliography: clear the % CHECK |

## Sources
- **Repository.** All at 7f5d4b3, files under W/repo:
  - manuscript/paper17/{paper17.tex, supplement.tex, README.md, numcheck17.py};
  - results/chaos/{transit_certificate, transit_limit}.md and their _raw.npz;
  - experiments/exp_transit_certificate.py (lasso_fs, blasso_verdict, jets, bp_jet_certificate, part_G, job_1d) and experiments/exp_transit_limit.py (entry_shift, Part W);
  - tests/test_transit_{limit,certificate}.py;
  - PARADIGM.md l. 6549ff; review/REVIEW_LOG.md ll. 5191ff and 5483ff; review/REVIEW_round140_D3ei.md.
- **Literature [read].** Kunin et al., arXiv 2506.06489v4: PDF fetched into W/dl_agf and read via pdftotext only.
- **My scripts.** W/own/{prop34.py, creep_md.py, check_IJ.py, check_G.py, reasons_I.py}. Build in W/build; reports regenerated in W/tc and W/tl.

## Lesson
- Count failures by type before describing them. A recovery rule with a hard tolerance turns the BLASSO's O(λ) support shift into "failures" that a one-line summary then attributes to the wrong cause; an exponent typo (10⁻¹ for 10⁻⁴) hides in the same sentence.
- A claim of the form "A < 1 exactly where B < 1" needs at least one case where B ≥ 1, or it is a tautology.
- When a limit object is defined algorithmically, map it onto the nearest algorithmic framework (here AGF) term by term. The correspondence can turn a "new object" claim into a stronger "first proof of the framework's conjecture for this model".
