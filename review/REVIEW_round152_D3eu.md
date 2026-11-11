# Referee report D3eu, round 152: the Paper 17 manuscript, third review

## Scope
- **What I reviewed.** manuscript/paper17/paper17.tex (34 pp.) and supplement.tex (17 pp.), the shipped PDFs, README.md and numcheck17.py at commit 868e774 (branch claude/vigilant-clarke-hf939b). This is the version revised after round 149 (review/REVIEW_round149_D3er.md).
- **Revision records read.** The README's "Revision after round 149" and "The plane's least certificate and BLASSO", and REVIEW_LOG's Round 149 entry and its plane entry.
- **Venue considered.** SIAM J. Math. Data Sci. or Found. Comput. Math.
- **Workspace.** W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref152. It holds the repo export (W/repo, venv symlinked), my own scripts (W/own), the regenerated report (W/tc), the rebuild (W/build) and the AGF download (W/dl_agf, read with pdftotext only). The repository was not modified.
- **Out of scope.** Bibliographic placeholders and the `% CHECK` items.

## Status of the round-149 issues

| D3er id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 AGF | **yes, but stated too broadly** | §1.3 paragraph (l. 125), claim (ii) (l. 139), §1.2 item 3 (l. 85) and the abstract now state the correspondence. I re-derived the threshold: 3hk·a_max² = 2Gς holds exactly, so 1/(3h‖θ₀‖) = (2/π)λ₁ [own]. | Scope of the claim: D3eu-1 |
| 2 proofs only in the supplement | yes | Appendix A covers Thm 3.2, Lemma 5.2 and Thm 5.3 on pp. 26–32. Sections 1–8 run to p. 26. Rebuild: 34 + 17 pp., 0 undefined references, no "??". | Abstract length: D3eu-6 |
| 3 (b)'s window end | yes | Step 2 now runs to τ_{k+1} − ε^{1/2}X′/v_{k+1}. (b) adds "for ε ≤ v_{k+1}⁴". The arithmetic checks: X′ = v·ε^{−1/4} ≥ 1 ⇔ ε ≤ v⁴. | One imprecision: the next entrant's shift σ_{k+1} is not accounted for (D3eu-3) |
| 4 C₁ and the layer's position | partly | C₁ = Cr·e^{λ_k s_r}·e^{λ_k(s^θ+σ)₊} is correct, and Remark 5.4 is updated. | The same exponential dependence now sits in Step 3's near-field Gronwall factor: D3eu-4 |
| 5 the split F_B / F_× | yes | I re-derived it: F_B = ½⟨(f_B*−f̂_B)²⟩ + ⟨(f_B*−f̂_B)Aφ*⟩, F_× = −w_e⟨(f_B*−f̂_B)φ_e⟩. Then dV/ds ≤ −cV + Cδ\|w_e\|V^{1/2} holds. | The far field's "integrable difference" omits the Gronwall step (cosmetic, D3eu-9) |
| 6 Part I's failure description | yes | Part K recounted independently (see Verification). §7.5 and §S7.5 are correct, and Table S5 is correct. | — |
| 7 Part J's vacuous sentence | yes | Removed from the paper. The oracle grid is named for the line. | Name the oracle grid for the plane too (D3eu-8). The md's Part J still says the plane LP "was not solved" (cosmetic) |
| 8 §7.2's order | yes | §S7.2: (εb/(ςμ²))² = 6.78·10⁻⁴ at μ = 0.00392; εb/(ςμ²) = 1.98 at 4.49·10⁻⁴ [own]. The closed-form values agree with D3er's 400-bit values. | — |
| 9 Step 3 as an estimate | **yes** | The Gronwall estimate is sound (verdict (b) below). | D3eu-3, 4 and 5 |
| 10 README | yes | l. 9 and "Numbers" are consistent. My numcheck run: paper 115/98 tokens, 0 not found; supplement 278/265, 9 not found, all derived (as stated). | README says "6558–6559" constraint pairs; the stored range is 6558–6560 (D3eu-7) |
| 11 Part G understated | yes | Recount: 1034/1034 sequences identical; 120/120 classes in the plane [own]. | — |
| 12 notation | partly | q, p, λ and η were handled. | Appendix A.3 introduces new collisions: D3eu-5 |
| 13 AGF venue | yes | arXiv v4 (24 Dec 2025) is still the latest version and prints "NeurIPS 2025" [read]. | — |

## The new content (verdicts per item)

**(a) The AGF correspondence.** I read Kunin et al., arXiv 2506.06489v4 (the latest version), §§1–3 and §6.

What the paper gets right:
- The paraphrase of AGF matches the source:
  - S_i = ∫κU_i;
  - directional rate ‖θ‖^{κ−2};
  - threshold c_i = (‖θ₀‖^{2−κ} − 1)/(κ − 2);
  - acceleration η = E[c_i];
  - collapsed neurons return to dormancy.
- The diagonal-network proof goes through Pesme–Flammarion (Thm 3.1).
- AGF calls the general correspondence open at the end of §2.2 and in §6.
- The threshold identification holds: threshold ∝ (εb_i)^{−1/2} ∝ Λρ_i, and the factor is exactly 2/π at the stated effective norm.
- The directional-rate comparison holds: ∝ w ∝ ‖θ‖³ under (2.2), and ∝ a ∝ ‖θ‖ under WFR.

Where the claim goes too far:
1. **The model is outside AGF's class.** AGF's conjecture concerns Euclidean gradient flow of a two-layer network f_i = a_iσ(w_iᵀg_i(x)) from a small random start. Equation (2) relies on f_i being homogeneous in θ_i.
   - Here the metric is Fisher–Rao in β.
   - The centre is a non-homogeneous parameter (AGF's g_i is fixed).
   - The amplitude saturates.
   - The start is θ = 0, regularised by the mass, so AGF's c_i(‖θ_i(0)‖) is infinite at the actual start.
   - The "up to π/2" is therefore a matching at a chosen effective norm, not AGF's threshold. Started from the norm the core actually has (at η = λ₁/2), the power law would give 1.137λ₁. What is exact is the ratio of the thresholds, ρ_i.
   - The relative thresholds, which are what determine the path, are exact.
2. **The sign.** With the scalar θ_i = u_i of degree 3, AGF's dormant law u̇ = 3hu²g cannot change sign from θ₀ ≠ 0. On the path, however, a dual may change sign between entries (Thm 5.3(d)), and the entrant takes its dual's sign. So S_i = 3s_iΛu_i holds with the sign chosen at activation, as in AGF's κ = 2 Pesme–Flammarion case (log cosh → \|·\|), not as AGF's κ_A = 3 radial law.
3. **What the theorem covers.** Theorem 5.3 assumes (H1): simple entries and no exits. AGF's collapses, which §1.3 maps to exits, are therefore excluded. It also assumes the cost step converges (H3).
4. **Internal inconsistency.** §1.3 says "Definition 2.3 is this algorithm for our unit, term by term" and, three sentences later, "The one difference is the directional step".

The abstract's "a case of that framework's conjectured correspondence" is acceptable only with these qualifications. §1.3's "Theorem 5.3 proves it" and claim (ii)'s "a proof of that framework's open correspondence" overstate them. See D3eu-1. The round-149 report (D3er-1) proposed exactly this stronger wording; I correct it here.

**(b) Appendix A.**
- **A.1 (Theorem 3.2).** Unchanged in substance. I re-checked the core bound g_ε = O(ε^{3/2}μ^{−2}), the bound α ≥ 2Gς(η − λ₁) beyond the wall, and the energy argument. Correct.
- **A.2 (Lemma 5.2).**
  - (i) and (i′) are correct; I re-derived η_c′ = 1/(Gβ_c a′).
  - (ii) is correct.
  - (iii) is correct with the new split. The far field is at sketch level (D3eu-9). The bootstrap, LaSalle and the implicit function theorem are sound.
- **A.3 (Theorem 5.3).**
  - Step 1: correct, including e′ = O(εξ²) from εb/β_c with β_c ∝ ξ⁻².
  - Step 2:
    - V = (y − S)ᵀH(y − S) decays at exactly 2λ_k.
    - The C₁ dependence is correct.
    - The extension to X′ is correct up to the next entrant's shift (D3eu-3).
  - Step 3, the key inequality \|D(ξ)\|_H ≤ M(\|D(X)\|_H + ∫\|b\|), M = e^{C/(2ξ₁²)}, is correct:
    - The block's linear part (1/v)(−MH) is dissipative in the H-norm, with the H of the coordinates of y (D3eu-5).
    - Every coupling carries w_e, 𝒲′ or dc_e/dξ, so it is O(ξ⁻³), which is integrable.
    - The dependence of the denominator \|dξ/ds\| on W multiplies dW/ds = O(ξ⁻³).
    - Source (1) is bounded generously.
    - Source (2) integrates to ε·log X.
    - With X = ε^{−1/4}: \|D(X)\| = O(ε^{3/4} + g′ε^{1/2}/v), giving \|D\| = O(ε^{1/2}).
    - The arrival integral gives σ̂ = O((1 + g′²/v³)ε^{1/4}) + O(εX³/v), as stated.
    - Re-synchronising ξ_e at ξ₁ removes the old "argued translation": the near-field comparison starts with \|D(ξ₁)\| + \|e(ξ₁)\| = O(ε^{1/2}).
  - Weakness: the near field uses a plain Lipschitz Gronwall over [s₀(ξ₁), s₂]. Its length contains (s^θ + σ)₊ (up to 9.856·10⁴), so the threshold ε₀ below which (b)–(d) hold is exponentially small on slow entries (D3eu-4). Table 1's "proved under (H1)–(H4), with constants that are not uniform" remains accurate.

**(c) Part K.** The classification is sound.
- No failure has a nonzero between 0.0500 and 0.1977 of its feature, so the classes are unchanged for every radius in that interval; I checked this independently.
- The order lost > spurious > displaced > wrong sign is applied consistently.
- The text uses it correctly:
  - 159 lost (B 49, D(a) 30, D(c) 80);
  - 207 displaced, in every family;
  - 147 = 306 − 159;
  - at ℓ = 10⁻⁴: 2 displaced + 20 spurious, all 20 on targets where η_V is degenerate.

**(d) Part L.**
- The tensor reduction is exact. ⟨u_i⊗u_j, (G₁⊗G₁)e_{ab}⟩ = s_is_j V_{ai}V_{bj}, so the grid values are (V₁⊗V₁)w exactly. The feature values are the projection, which equals ⟨zφ(c*)⟩ exactly for z in the span.
- The dual bound is valid. Adding δ = B(BᵀB)⁻¹(−r) split into its parts, then renormalising, gives exact dual feasibility; weak duality gives t ≥ −s·ν. It is certified up to floating point.
- Wording problem: "computed exactly" (§7.5 l. 573, p. 25) and the five printed digits overstate the result (D3eu-2):
  - The "least t" is the optimum over the span truncated at 10⁻¹³.
  - t rises monotonically, on 120 of 120 targets, as the truncation coarsens: medians 0.620 → 0.696 → 0.780.
  - The 10⁻¹³ value is therefore only an upper bound for the least t over the full 1681-dimensional span.
  - The bracket certifies the truncated program, not the full one.
  - Rounding moves the direct sup by up to 1.26·10⁻³.
  - Exactness on the grid was already proved by t_J ≤ 0.98990 < 1 (Part J). Part L adds only the value, and that value is mostly a grid-geometry effect (0.602 against 0.724), as the supplement says.

**(d′) Part M.** Correct. The KKT conditions hold under my own Gram matrix, and the counts recount exactly.

## New findings
1. **AGF scope** (abstract l. 56; §1.2 l. 85; §1.3 l. 125, p. 5; claim (ii) l. 139). See (a); this is D3eu-1.
2. **Part L "computed exactly"** (§7.5 l. 573, p. 25; §S7.5 l. 336–338; md Part L title). Truncation-dependent; see (d). This is D3eu-2.
3. **Step 2's next-entrant bound** (App. A.3 Step 2, l. 655, p. 29).
   - The bound uses the path's margin v_{k+1}(τ_{k+1} − τ).
   - The loop's next entrant arrives at τ_{k+1} + ε^{1/2}σ_{k+1} + O(ε^{3/4}), so its inner coordinate at the interval's end is X′ + v_{k+1}σ_{k+1}, not X′.
   - "X′ ≥ 1" should read X′ ≥ 1 + v_{k+1}(σ_{k+1})₋, with the loop's margin in 𝒫. This is harmless asymptotically, since X′ → ∞.
4. **Remark 5.4 / Step 3 near field** (l. 451, p. 18; l. 683, p. 31).
   - The near-field Gronwall bound is e^{L(s₂ − s₀(ξ₁))}·O(ε^{1/2}), with s₂ ⊇ (s^θ + σ)₊.
   - So "(b) holds from τ_k + ε^{1/2}(s^θ+σ)₊ + Cε^{1/2}log(1/ε) with a constant independent of the layer's position" is true of the prefactor. It is not true of the ε₀ below which the statement applies, which is still exponentially small in the layer's length on slow entries.
5. **New notation collisions in App. A.3** (l. 661–690).
   - Z (the continuum normaliser of Table 2) is reused as the state.
   - W (also used for W_i^ε and 𝒲_ρ) is reused as the far-field state.
   - 𝒜(ξ) clashes with the active sets 𝒜_k (the collision D3ei-12 separated).
   - M (Gronwall constant) clashes with M_k (mobilities).
   - b(ξ) clashes with the masses b_i.
   - e(ξ) (profile offset) and the subscript e (entrant) clash with e_k and E of §6.
   - s₁ means two different times in Lemma 5.2 (ii) and (iii).
   - H_{k−1} in \|D\|_H must be the Hessian in y's (α, c) coordinates; (H2) defines it in (w, c). They agree only at critical points and after the change of variables.
6. **Abstract length.** About 380–420 words against SIAM's 250-word limit; FoCM also expects a short abstract.
7. **Constraint pairs.** §S7.5 says "6558 pairs" and the README says "6558–6559"; the stored n_con is 6558 (62 targets), 6559 (44) and 6560 (14). Also, the md's Part J sentence "the full linear program in the plane … was not solved" is now superseded by Part L, with no pointer.
8. **Oracle grid in the plane.** §7.5's "In the plane basis pursuit on the grid is exact" should say the features were added (oracle grid), as it does for the line.
9. **Lemma 5.2(iii)'s far field** (l. 648). "Rates differ by O(δ\|w_e\|), integrable" skips the Lipschitz/Gronwall step, which needs the state difference too. Citing Step 3's estimate (dissipative block, O(ξ⁻³) couplings) would close it.

No error was found in Theorem 3.2, Lemma 5.2(i)–(iii), the orders of Theorem 5.3, or any number in the paper or the supplement. Table 1's statuses are accurate.

## Verification I did
- **[own] Reproduction.**
  - Copied the raw into W/tc (TC_OUT) and ran `--from-raw`: transit_certificate.md is byte-identical (`cmp`), and the raw is unchanged.
  - pytest tests/test_transit_certificate.py: 54 passed (94 s, nice, one process).
- **[own] Part K** (W/own/check_KM.py: own rule from §S7.5's text, no import of the experiment code).
  - "all" rows at ℓ = 10⁻¹…10⁻⁴:
    - not recovered: 366, 60, 34, 22;
    - lost: 159, 24, 6, 0;
    - spurious: 0, 20, 20, 20;
    - displaced: 207, 16, 8, 2;
    - wrong sign: 0;
    - empty within tolerance: 306, 57, 13, 2.
  - Displaced maximum 0.0500, spurious minimum 0.1977.
  - At ℓ = 0.1: A displaced 63; B lost 49 / displaced 24; D lost 110 / displaced 120. All exactly as printed.
  - Sequences: 1034/1034 identical; plane classes 120/120.
- **[own] Part M** (same script). I built my own samples (41², σ = 0.15), my own grid (spacing 0.025, features added) and my own correlations.
  - λ_max reproduced to 10⁻¹².
  - Worst KKT violation 6.75·10⁻¹² relative to λ (printed: 7.18·10⁻¹²).
  - Recovered 118/120/120/120; 2 displaced at ℓ = 0.1; support exactly the features on 64 at every ℓ; at most 6 atoms beyond the features.
  - η_V nondegenerate on 120/120 (raw H).
- **[own] Part L** (W/own/check_L.py: a dense SVD of the 1681 × 6561 grid matrix, not the Kronecker factors; HiGHS, not the paper's interior-point method).

  | target | truncation | dimension | constraint pairs | t, HiGHS (mine) | t, stored | t_low, stored |
  |---|---|---|---|---|---|---|
  | C 0 | 10⁻⁷ | 571 | 6558 | 0.7574492598 | 0.7574492603 | — |
  | C 63 | 10⁻⁷ | 571 | 6560 | 0.8515063804 | 0.8515063810 | — |
  | C 7 | 10⁻¹⁰ | 839 | 6558 | 0.6501637577 | 0.6501630185 | 0.6501622463 |

  - The dimensions and constraint counts match the stored values.
  - On C 7, HiGHS lies 7·10⁻⁷ above the stored value, within its feasibility tolerance (my equality residual is 1.2·10⁻⁷).
  - Each certificate evaluated directly on the samples gives t_dir equal to the stored t.
  - W/own/geom_L.py: 62/58 split with medians 0.60173 and 0.72392 and the stated ranges; overall range 0.56686–0.78466, median 0.62010; t_J − t ∈ [0.2050, 0.4087], median 0.3601; monotone in the truncation on 120/120; brackets ≤ 5.43·10⁻⁶.
- **[own] AGF.**
  - Threshold algebra symbolic and numeric (ε = 10⁻⁸: 44 721 against λ₁ = 70 248, ratio 0.6366 = 2/π).
  - The h, k identity holds.
  - The dormant rates ‖θ‖³ (under (2.2)) and ‖θ‖ (under WFR) check.
- **[own] Derivations.** F_B/F_× split; V's rate 2λ_k; e′ = O(εξ²); Step 3's \|D(X)\|, ∫\|b\| and σ̂ orders; C₁; X′ = vε^{−1/4}.
- **[own] Build.** Paper and supplement built alternately, three times each: 34 + 17 pp.; 0 undefined references; 0 overfull boxes; no "??"; text identical to the shipped PDFs. Appendix A is on pp. 26–32 and the references on pp. 32–34.
- **[own] numcheck17.** Paper 115/98, 0 not found. Supplement 278/265, 9 not found, all derived in §S7.2.
- **[read]** README (both new sections), REVIEW_LOG (Round 149 and the plane entry), the code of Parts K, L and M.

## Overall assessment
**Minor revision.**
- All 13 round-149 issues are addressed.
- The new Appendix A is correct. Step 3 is now a real Gronwall estimate, and the translation is fixed rather than argued.
- Parts K, L and M reproduce byte for byte and agree with independent recomputation.

What remains is precision:
- the AGF claim should be scoped (the unit is outside AGF's class, the sign is chosen at activation, there are no exits under (H1));
- "exact t" in the plane should say "at the truncation 10⁻¹³";
- two effectivity remarks on Theorem 5.3;
- notation;
- an abstract over SIAM's limit.

For SIMODS, acceptable after these minor edits. For FoCM, the rigorous core (Theorems 3.2 and 5.3 under (H1)–(H4), and Lemma 5.2) now stands in the paper itself. Its non-uniform constants are honestly flagged.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3eu-1 | minor (positioning) | The AGF claim overstates its scope. The unit is outside AGF's class: non-Euclidean metric, non-homogeneous centre, saturating amplitude, θ₀ = 0 start, so "π/2" is a matching at a chosen norm. A degree-3 scalar AGF unit keeps its sign, while the path's duals may change sign (Thm 5.3(d)). (H1) excludes exits/collapses. §1.3 says both "term by term" and "the one difference". | §1.3 l. 125, claim (ii) l. 139, §1.2 l. 85, abstract: "an AGF-type correspondence (thresholds ∝ ρ_i, accumulated utility = the dual with its sign chosen at activation as in AGF's κ = 2 case, cost step = relaxation; directional step frozen), proved for paths with simple entries and no exits under (H2)–(H4); an instance beyond AGF's stated setting, not a case of its conjecture as stated"; drop "term by term" |
| D3eu-2 | minor | Plane t "computed exactly" with 5 digits. t depends on the truncation (medians 0.620, 0.696, 0.780 at 10⁻¹³, 10⁻¹⁰, 10⁻⁷, monotone on 120/120), so the 10⁻¹³ value only bounds the full-span least t from above. The bracket certifies the truncated program. Rounding is 1.26·10⁻³. | §7.5 l. 573: "solved to optimality at the line's truncation 10⁻¹³"; 3 digits; one sentence on the monotone truncation dependence; same in §S7.5 and the md L title ("exactly") |
| D3eu-3 | minor | Step 2's next-entrant bound uses the path's margin. The loop's next entrant arrives ε^{1/2}σ_{k+1} later, so its coordinate is X′ + v_{k+1}σ_{k+1}. | App. A.3 Step 2: require X′ ≥ 1 + v_{k+1}(σ_{k+1})₋ and use the loop's margin in 𝒫 |
| D3eu-4 | minor | The layer dependence moved from C₁ to the near-field Gronwall factor e^{L(s₂ − s₀(ξ₁))}, with s₂ ⊇ (s^θ+σ)₊ ≤ 9.856·10⁴. ε₀ remains exponentially small on slow entries. | Remark 5.4: say so; optionally, in Step 3's near field, use the block's dissipativity until the entrant grows |
| D3eu-5 | minor/cosmetic | New notation collisions in App. A.3: Z, W, 𝒜(ξ) vs 𝒜_k, M vs M_k, b(ξ) vs b_i, e(ξ)/subscript e vs e_k and E, s₁ twice in Lemma 5.2. H_{k−1} must be the Hessian in y's (α, c) coordinates. | Rename (e.g. 𝐙, 𝐘, 𝔅(ξ), M_G, β(ξ), o(ξ)); one line on H's coordinates |
| D3eu-6 | minor (venue) | The abstract is about 380–420 words, against SIAM's limit of 250. | Cut to at most 250 words: one sentence each for the theorem, AGF, the finite mass and recovery |
| D3eu-7 | cosmetic | "6558 pairs of constraints" (§S7.5 l. 336) and README's "6558–6559": the stored range is 6558–6560 (62/44/14). The md's Part J says the plane LP "was not solved", with no pointer to L. | §S7.5, README; one line in md section L |
| D3eu-8 | cosmetic | §7.5: the plane's basis pursuit grid is an oracle grid (features added), not stated there. | §7.5 l. 573 |
| D3eu-9 | cosmetic | Lemma 5.2(iii)'s far field skips the Gronwall step on the state difference. | App. A.2: cite the argument of Step 3 (dissipative block, O(ξ⁻³) couplings) |

## Sources
- **Repository** (868e774, under W/repo):
  - manuscript/paper17/{paper17.tex, supplement.tex, README.md, numcheck17.py, PDFs};
  - results/chaos/transit_certificate.md and transit_certificate_raw.npz;
  - experiments/exp_transit_certificate.py (k_classify, part_K, _kron_basis, _dual_lower, bp_certificate_kron, KronGram, blasso_kron, lasso_fs, blasso_verdict, bp_certificate);
  - experiments/exp_transit_nd.py (full);
  - tests/test_transit_certificate.py;
  - review/REVIEW_round149_D3er.md;
  - review/REVIEW_LOG.md ll. 5597ff and 5720ff.
- **Literature [read].** Kunin et al., arXiv:2506.06489v4 (24 Dec 2025, the latest; the history lists v1–v4): §2 eqs. (1)–(3), Algorithm 1, §2.2's "a general proof remains open", §3 Thm 3.1, §6. The PDF is in W/dl_agf, read via pdftotext.
- **My scripts.** W/own/check_KM.py, check_L.py, geom_L.py; outputs in W/own/check_L.out and check_L2.out.

## Lesson
- When a referee proposes that a limit "is" a known framework's algorithm, check the framework's model class and its sign and initialisation conventions before adopting "proves the framework's conjecture". A correspondence at the level of thresholds and steps can hold while the model lies outside the conjecture's stated setting.
- Report an "exact" optimum of an ill-conditioned linear program with its truncation and its monotone dependence on it. A certified bracket certifies the program as posed, not the problem one has in mind.
- When a constant is made explicit in one step of an asymptotic proof (here C₁), follow where the same dependence reappears (here the near-field Gronwall factor).
