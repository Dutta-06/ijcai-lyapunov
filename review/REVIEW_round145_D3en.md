# Referee report D3en, round 145: the Paper 15 manuscript, second review

## Scope
This is the second review of `manuscript/paper15/paper15.tex` (26 pp.; §1–7 end on p. 23) and `supplement.tex` (11 pp.) at commit f8ba0a1, together with the PDFs, `README.md` ("Revision after round 142"), `make_figs.py`, `numcheck15.py`, the Round 142 entry of `review/REVIEW_LOG.md` and the first report `review/REVIEW_round142_D3ek.md`.

- The commit was exported to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/repo` with the venv symlinked. The repository was not modified.
- Bibliographic placeholders and `% CHECK` items are out of scope. The one exception is Katzenberger's hypotheses for non-gradient fields, which I was asked to check; I treat it as a question about the mathematics.
- Evidence tags: [read] read in the files; [md]/[raw] from the reports or the raw npz; [own] my own computation; [paper] a work I read in full text; [snippet] web search snippets only; [memory] my background knowledge, not checked.
- Line numbers are those of `paper15.tex` unless an "S" is prefixed (`supplement.tex`).

## Status of the round-142 issues

| D3ek id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 (formal labelled proved) | **Yes, largely.** The new proof of Thm 3.4(b) is correct. | Table 1 relabelled. Assumption 2.5 is stated, and every Euler-update result (Prop 4.1, Thm 4.2, Cor 4.3, Thm 5.1, Prop S1.1) is "derived under" it. Thm 3.4(b): I re-derived ½∂²_ηηΦ̃[Σ̃] = ⟨∇ψ·𝒢(A)U⟩ by hand and checked it numerically on a generic 3-parameter, 2-node non-gradient-integrable model with OU noise [own] (Verification). Per-node Katzenberger terms of the augmented limit map ≈ 9.1·10⁻⁴; they match ϖ_x∇ψ_x·𝒢(A)U_x to 5·10⁻⁶ relative; their sum cancels to 5·10⁻⁶ of one term; ∂_ηΦ̃ ≈ 10⁻⁹ (no diffusion). | (i) The C² assumption on Φ̃ is very likely removable (D3en-5). (ii) Katzenberger's non-gradient hypothesis can now be supported from sources (D3en-5). (iii) Lemma 3.3 is stated for analytic 𝗀, but (c) uses arbitrary ϱ (D3en-9). |
| 2 (independence over nodes) | **Yes.** | §1.1 l. 63 and §2.3 l. 159 state what the hypothesis means. Prop 3.7's formula follows from the proofs with C replacing δ/ϖ. The S3 identity B_C + B_Cᵀ = −(B̃ + B̃ᵀ) is correct. The OU-offset expansion ½(1 − τA)⁻¹ ≈ −τ⁻¹A_D/2 is correct. On my generic model the common offset's drift is 2.6·10⁻⁴, 28% of a per-node term, so it is clearly nonzero [own]. | none |
| 3 (AR(1) closed form) | **Yes.** | Exact-rational check against Thm 4.2's f_il for ρ ∈ {−0.9, …, 0.99} and x up to 1.9 [own]. f(x*) = 8ρ(1−ρ)/(1+ρ)² is exact. The grid minimiser matches x*. c₂ = (7ρ²−4ρ+1)/(4(1−ρ)²) agrees by finite differences. A brute-force double series at ρ = 0.7, x = 0.8 agrees to 4·10⁻¹⁶. Fig. 1(c) marks ½ and 17/18 correctly. | The abstract's "Ψ is concave only below hλ = (3ρ−1)/(2ρ)" holds for ρ > 1/3 only (minor). |
| 4 (Prop 5.6's reach) | **Partly.** | Prop 5.7(b), (c) are correct; I checked Q_∞ and the certificate against my own implementation on all 100 K2 ends to 1.3·10⁻¹² [own]. The restored "basins" sentence and Part K2 are good. | (a)'s general-position hypothesis is self-contradictory as written (D3en-2). "The failures are basin failures" omits the τ_c = 100 column, where 38% of the k = 6 ends have lower R than the truth (D3en-3). The abstract's "competing fits … remain local minima" holds only for Q_∞ < 1 (D3en-8). |
| 5 (§6 at saturations) | **Yes.** | §6 moved to S1. The saturation is stated in S1.1 and Remark S1.4, Prop S1.2's domain caveat is added, and the 1/a′ caveat is restored in Prop S1.3. | Supplement results have no status line (D3en-12). |
| 6 (cosines vs magnitudes) | **Yes.** | §6.2 l. 485 and S2 numerics print 1.8284–3.2693 and −6.9116 [md, matched]. | none |
| 7 (literature) | **Partly.** | Ziyin–Wang–Li–Wu, Egardt, Ioannou–Kokotović, Anderson and the others are added, and the old l. 99 sentence is gone. | The SSP paragraph (l. 99) is still written from the abstract only. arxiv.org is now reachable, and the full comparison is in Literature (D3en-1). The "isotropic noise does not select" sentence is inaccurate (D3en-10). |
| 8 (statistics) | **Yes.** | `colour_selection.md` regenerates byte-identically from its raw file. My exact McNemar p-values: 6.10·10⁻⁵, 1.53·10⁻⁵, 0.125, 0.6875, 0.03125, 0.0078. My Clopper–Pearson intervals match the printed ones (e.g. 30/40 → [0.588, 0.873]; 22/40 → [0.385, 0.707]) [own]. | none (0.12 for 0.125 is the md's rounding; disclosed) |
| 9 (α = 2 wording, horizon) | **Yes.** | l. 494, l. 500, S4 | none |
| 10 (length) | **Partly** (author declined the 20-page target). | 23 pp. main, 11 pp. supplement | Acceptable for SIMODS. The abstract is ≈346 words against SIAM's 250-word limit (D3en-13). |
| 11 (sign convention) | **Yes.** | S2.2(a) proof: r_{n+1} = (1−hλ)r_n − hλξ_n; θ = u² − v² | none |
| 12 (factor m) | **Yes.** | l. 356 | none |
| 13 (strict convexity) | **Yes.** | l. 412; S5 | none |
| 14 (order v; placement; λτ_c ≫ 1) | **Yes.** | l. 256, 351, 432 | Still overclaimed at l. 61 and l. 101 ("whatever/every temporal spectrum"; D3en-7). |
| 15 (Σ_f in Prop 4.1(b)) | **Yes.** | l. 306; I checked the sign of the added drift | none |
| 16 (room) | **Yes.** | l. 480 | none |
| 17 (ρ = 0.95, at-W₀) | **Yes.** | l. 483 (0.9919–1.0219, 0.9141, 1.4573, 1.0043–1.0265 all in `gauge_label.md`) | none |
| 18 (Heun) | **Mostly.** | Remark 4.5, S2.2(d) | "A quarter of it at hλ = 1" extrapolates a leading-order formula; the first referee's MC gave 12.9% at hλ_max ≈ 0.97 (D3en-14). |
| 19 (Lévy area; zero-set index) | **Yes.** | l. 297; S1.1 "Completeness"; S3 | none |
| 20 (notation) | **Mostly.** | Table 2 completed | New overloads: Π (stationary covariance in Thm 3.4(b), Π_D, Π^⊥), 𝖼 vs certificate c vs c_k, 𝖻 vs drift b, η (augmented state, η(W)), ν₀/ν_j/ν_k (D3en-15). |
| 21 (Fig. 2 SBL markers) | **Yes.** | Rendered fig2: stars on top, with interval bars | none |
| 22 (P_i ≠ 0) | **Yes.** | l. 442, l. 450 | none |

## New findings

### Major
None. The new mathematics is correct, and I verified it independently (below).

### Moderate
1. **The SSP comparison is outdated (l. 99; README "Open items").**
   - The paragraph still says the full text "was not reachable" and that the authors "cannot say whether their framework yields a Stratonovich/Itô distinction …".
   - The arXiv text (2404.12293v1, 49 pp.) is now readable, and I read it in full [paper]. It contains none of the three; details in Literature.
   - The paragraph must be rewritten on the basis of the text, and the priority statement made definite. This is easy, and the result favours the paper.
2. **Prop. 5.7's general-position hypothesis is self-contradictory as written (l. 456).**
   - It demands that, for all supports S, S′, span X_S ∩ span X_{S′} have dimension max(0, |S|+|S′|−N).
   - This fails for every pair of overlapping or nested supports. For example, with N = 10 and |S| = |S′| = 3 sharing one column, the intersection has dimension ≥ 1 but the formula gives 0. So no design satisfies the hypothesis once d ≥ 5, and the proposition is vacuous as stated.
   - The proof (l. 470) uses only "every N columns independent". From that, dim = |S|+|S′| − min(N, |S∪S′|), and the proof's case split is then correct; I re-derived (a) that way.
   - Fix: delete the second clause, or replace it by the correct formula as a consequence.
3. **"The failures are basin failures" (l. 496) is overstated for the flow actually run (τ_c = 100).**
   - K2's own table [md], whose figures I recomputed from the raw [own], shows that at τ_c = 100 the converged end has a lower R than the truth in 0.09 (k = 5) and 0.38 (k = 6) of the failures. At k = 6, 15 of the 40 failures end at a fit the τ_c = 100 potential prefers to the truth. There the truth is not the global minimiser, and the failure is the potential's, not the basin's.
   - The paper prints the τ_c = 1, 10 and ≥ 10⁴ columns and omits τ_c = 100.
   - Related: Prop 5.7(c)'s O(1/τ_c) remainder has a large constant. On K2's first k = 3 instance, τ_c·[τ_c(R−R*) − prediction] ≈ 3.4·10³, 9.1·10³ and 9.3·10³ at τ_c = 10³, 10⁵ and 10⁷, so the remainder is ≈ 3.4 at τ_c = 10³ [own]. The expansion therefore says nothing at τ_c = 100.
   - Fix: print the τ_c = 100 column. Say that the failures are basin failures for τ_c ≥ 10⁴ evaluation, and partly potential failures at the run's τ_c. State that (c) is asymptotic, with a constant of order 1/(τ_c λ⁺_min).
4. **The fixed-intensity scaling of the coloured selection is not said in the main text (Cor. 4.3(c), §1.1 l. 61).**
   - The log-det selection runs at v_h = v(1−ρ)/(1+ρ) ≈ vh/(2τ_c). At a fixed physical intensity v it is O(h) and vanishes as h → 0, which is how Cor 4.3(c) ("finely sampled, h → 0") is consistent with Theorem 3.4 (no selection in continuous time).
   - Only S2.1's last line says this (τ_eff ≈ mvh/(4τ_c²)). The paper's central dichotomy needs it stated where the reader meets both results.
   - Fix: one sentence after Cor 4.3(c) and in §1.1.

### Minor
5. **Thm 3.4(b): the hypotheses (l. 253, Remark 2.3, Table 1 row 1).**
   - *Katzenberger for non-gradient fields.* Parsons–Rogers (arXiv:1510.07031, App. B and §II) restate Katzenberger for a general drift f with a manifold of equilibria that attracts normally [paper]. Durrett–Popovic (arXiv:0902.4780, p. 11–12) apply Katzenberger to a non-gradient system after checking that the normal linearisation has all eigenvalues of negative real part [paper]. Remark 2.3's description agrees with both. I could not read Katzenberger 1991 itself (projecteuclid.org returns 403). SSP and LWA state the theorem only for −∇L [paper].
   - *C² of Φ̃.* LWA obtain C² of Φ from Falconer 1983, Thm 5.1, applied to the time-1 map φ(·,1) (Lemma B.1), which is not a gradient-specific argument [paper]. SSP cite [Kat91, Cor. 5.1] for C² under C³ data [paper]. The augmented field is C³ and Γ×{0} is normally attracting (spectrum {−λ_i} ∪ spec E). So the "assumed" C² very likely follows from [Fal83] or [Kat91]; check this and drop the assumption.
   - The scaling is right: the augmented noise √v ϖ^{-1/2}𝖻 dW is O(√v) against an O(1) drift on η, and the stationary η has size √v.
6. **Prop 3.1's no-diffusion is known for M = I (l. 224).** LWA Cor 5.2 (deterministic label-noise limit) and SSP Cor 5.11 ("∇f ∈ Range ∇²L … the noise term vanishes") contain it [paper]. Credit them, and claim the general mobility and the drift cancellation as new.
7. **"Whatever the temporal spectrum" (l. 61) and "for every temporal spectrum" (l. 101).** These are proved for white and rational spectra and formal otherwise. Qualify them. Likewise Thm 3.4's closing sentence (l. 256) covers (c), which is formal.
8. **Abstract, "competing fits on as many atoms as samples remain local minima".** This holds only for those with Q_∞ < 1: 0.30–0.56 of the observed ends; the rest lose the certificate as τ_c grows (K2: 1.00 → 0.30 at k = 3). Say "some (those with Q_∞ < 1)".
9. **Lemma 3.3 (l. 242) assumes 𝗀 analytic,** but Thm 3.4(c) applies it to "any correlation function ϱ". For a non-exponentially decaying ϱ, 𝗀(λ) = ∫ϱe^{uλ} is not analytic at λ = 0. The proof needs only 𝗀 defined on the (real, finite) spectrum of the diagonalisable A. State it so.
10. **Isotropic noise (l. 101): "isotropic parameter noise does not select among the fits, parameter-dependent noise does [HWLM21]".** LWA Cor 5.1 [paper] gives, for isotropic noise, a Brownian motion on Γ plus the drift −½∂Φ∇ ln|∇²L|₊: a selection. Reword to what HaoChen et al. show; my recollection [memory] is that Langevin noise does not give the sparse bias.
11. **Table 1 omits Lemma 3.3** (new, proved).
12. **The supplement's results (Props S1.1–S1.3, S2.1, S2.2) have no status line** since the rows left Table 1. Prop S1.2's averaging and Laplace steps are only sketched, but its statement is unqualified. Add a short status table to the supplement.
13. **The abstract is ≈346 words;** SIAM journals limit abstracts to 250.
14. **"A quarter of it at hλ = 1" (Remark 4.5, S2.2(d)) extrapolates the leading-order (hλ)²/4.** The h⁴c²c̃²/4 term and c − c̃ = hλc + h²c²uv matter there; the first referee's MC gave 12.9% at hλ_max ≈ 0.97. Say "(hλ)²/4 for hλ ≪ 1".
15. **Notation overloads introduced by the revision:**
    - Π: stationary covariance (l. 253), Π_D (l. 398), Π^⊥ (l. 426);
    - 𝖼 vs the certificate c vs c_k (S1);
    - 𝖻 vs the drift b (Thm 2.2);
    - η: augmented state vs η(W) (S1.1 proof);
    - ν₀ vs ν_j vs ν_k;
    - E: stable matrix vs 𝔼.
16. **Remark S1.4: "This is the multistability of the concave selection of Proposition 5.5 on this network".** Prop 5.5 concerns factorised coefficients, so this is an analogy; say so.
17. **Wording on SSP's limit (l. 99): "a gradient flow of an explicit functional".** That is SSP's non-degenerate case (Thm 4.1). Their degenerate case (Thm 4.7) is a constrained SDE, deterministic only for label noise.

## Verification I did
- **Reproduction [own].** I copied `colour_selection_raw.npz` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/out_cs`, then ran `CS_OUT=… exp_colour_selection.py --from-raw`. `cmp`: the md is byte-identical, and the raw is unchanged.
- **Tests and checks [own].**
  - `tests/test_colour_selection.py`: 11 passed.
  - `numcheck15.py results/chaos paper15.tex supplement.tex`: 397 tokens, 1 not in the md (5.4, derived and declared), as logged.
  - I spot-checked 18 new §6.2/S1 numbers in `gauge_label.md`; all found.
- **Part K [own, `scripts/k2_check.py`].**
  - I rebuilt the 100 K2 instances from their seeds. With my own ω, c and Q_∞, the certificate ratios at τ_c = 10², 10³, 10⁴ and 10⁵ and Q_∞ match the raw to 1.3·10⁻¹².
  - The certificate at τ_c = 10⁷ tends to Q_∞ (e.g. 1.54310 vs 1.54331).
  - Fractions recomputed: all 100 ends are on 10 atoms and certified at τ_c = 100; Q_∞ < 1 in 0.30, 0.5625, 0.382, 0.475; R(end) < R(truth) at τ_c = 100 in 0, 0, 0.088, 0.375.
  - τ_c(R − R*) at 10⁴ reproduces the raw (27.960 etc.).
  - McNemar p and Clopper–Pearson intervals as in the status table.
- **Thm 3.4(b) [own, `scripts/aug_check.py`].** Model: f̂₁ = z₁ + 0.3z₂z₃ + 0.2z₁², f̂₂ = z₂ + 0.4 sin z₁ z₃ + 0.1z₃²; M = I; ϖ = ½; normal rates 0.409 and 0.773; OU noise with τ_c = 2 (𝖻 = 1/τ_c, so ∫ϱ = 1). The augmented flow was integrated with DOP853 (rtol 10⁻¹³) and second η-derivatives taken by central differences. Results as in the D3ek-1 row: per-node terms equal the formula, the sum cancels, there is no diffusion, and the common offset gives a nonzero drift.
- **Prop 4.4 [own, `scripts/ar1_num.py`]:** as in the D3ek-3 row (sympy is not in the venv, so I used exact `Fraction` arithmetic).
- **By hand [own].** The S3 identity for B_C; the OU-offset expansion; S2.1(b)'s Lyapunov solution (P^⊥A = A, AMP^⊥ᵀ = AM); Prop 4.1(b)'s sign; Prop 5.7(b)'s limits and (a)'s proof with the corrected hypothesis.
- **Build [own].** In a scratch copy: three passes each, 0 errors, 0 undefined references, 0 overfull boxes; 26 + 11 pp. The text of both PDFs is identical to the committed ones. I rendered Figs 1 and 2 and checked them visually.

## Literature
**Shalova–Schlichting–Peletier, arXiv:2404.12293v1 (the only arXiv version) [paper].** I read it in full. JMLR (jmlr.org) was unreachable, so possible changes in the journal version are unchecked.

Their setting is discrete GD, w_{k+1} = w_k − α∇_wL̂(w_k, η_k), with η_{k,i} i.i.d. across steps and coordinates, analysed as α → 0 through Katzenberger's theorem, which they state for −∇L only (Thm 3.15).

| SSP result | content | relation to this paper |
|---|---|---|
| Thm 4.1 (non-degenerate, time scale ασ²) | constrained gradient flow of Reg = ½Δ_ηL̂(w,0) (dropout, DropConnect, L(w+η)) | Noise not through the error. No counterpart here. |
| Thm 4.7 (degenerate, time scale α²σ²) | constrained SDE | Contains LWA's label-noise result as Cor 5.11, Reg = Δ_wL/(2N) as α → 0 with fixed σ allowed: this paper's Prop 4.1(a) at h → 0 (Itô). |
| Cor 5.11 (label noise) | the noise term vanishes by orthogonality | Prop 3.1's no-diffusion for M = I (D3en-6). |
| Cor 5.10 (minibatch) | trivial evolution on both time scales | A different mechanism: the per-sample gradients vanish on Γ, so the noise vanishes there. Here the noise is nonzero, normal, and its drift cancels. |
| §6.2 "Correlated noise" | correlation across the coordinates of one η_k, fixed per step: a reweighted Hessian trace | Not temporal correlation; no filter. |
| Remark 4.9 | Itô corrections | No Stratonovich reading, no Wong–Zakai, no smooth noise, no fixed step, no temporal colour anywhere in SSP. |

Conclusion: SSP contains no no-selection result of Theorem 3.4's type for any noise class, no Stratonovich/Itô distinction, and no coloured filter. The paper can claim Theorem 3.4, Theorem 4.2's filter, Prop 4.4 and Prop 3.7 without qualification relative to SSP, while crediting LWA and SSP for the label-noise flow and the no-diffusion (D3en-1, D3en-6).

**Other works.**
- **Li–Wang–Arora [paper].** Cited correctly: Katzenberger for SGD near zero loss, arbitrary covariance, η → 0. Their label-noise flow −¼∂Φ∇tr[c∇²L] matches Prop 4.1(a). Their Cor 5.1 contradicts the "isotropic noise does not select" sentence (D3en-10). C² of Φ comes via Falconer on the time-1 map (D3en-5).
- **Ziyin–Wang–Li–Wu (arXiv:2402.07193v3, NeurIPS 2024) [paper].** "Gradient noise creates a systematic motion … along the degenerate direction to a unique … fixed point" is cited correctly. Their model is a continuous-time Itô SGD, with a discrete-time extension in App. D. Worth one clause: their Noether flow is an Itô effect, which Thm 3.4(a) removes in Stratonovich's reading of label noise.
- **Katzenberger 1991.** Not readable (Project Euclid blocked). For the non-gradient form I relied on Parsons–Rogers and Durrett–Popovic [paper], as above.

## Overall assessment
**Minor revision.**

- **The revision is thorough and honest.** Every round-142 item was acted on, and the closed forms and the number provenance are exemplary.
- **The new proof of Thm 3.4(b) is correct, and I confirmed it independently on a model with no symmetry.** I re-derived it analytically and checked the augmented Katzenberger terms to 5·10⁻⁶. Its two stated hypotheses look dischargeable from [Kat91]/[Fal83] and secondary sources.
- **Props 3.7, 4.4 and 5.7(b), (c) and Part K are right and reproducible.** I reproduced the md byte for byte and recomputed the certificate ratios to 10⁻¹².
- **What remains is narrow:**
  - rewrite the SSP paragraph from the now-available text;
  - correct Prop 5.7's hypothesis;
  - report the τ_c = 100 column before calling the failures basin failures;
  - make the fixed-intensity O(h) scaling explicit;
  - fix a handful of wording and literature points.

None of these changes a result.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3en-1 | moderate | The SSP paragraph is still from the abstract. The full text [paper] has no Stratonovich/Itô distinction, no temporal colour, no fixed step, no no-selection of this type. Its label noise = LWA (Cor 5.11); its minibatch triviality has a different mechanism. | l. 99, README Open items: rewrite per Literature; make the priority claim definite. |
| D3en-2 | moderate | Prop 5.7's "spans meet in dimension max(0, \|S\|+\|S′\|−N) for all S, S′" fails for any overlapping or nested supports, so the hypothesis is vacuous. The proof needs only "every N columns independent". | l. 456: drop the clause, or state dim = \|S\|+\|S′\|−min(N, \|S∪S′\|). |
| D3en-3 | moderate | "The failures are basin failures": at the run's τ_c = 100, R(end) < R(truth) in 0.09/0.38 of the k = 5/6 failures [md K2, own]. Prop 5.7(c)'s O(1/τ_c) remainder is ≈ 3.4 at τ_c = 10³ on an instance [own]. | l. 496: print the τ_c = 100 column; qualify. l. 464: the remainder's constant. |
| D3en-4 | moderate | At fixed intensity the coloured selection runs at v_h ≈ vh/(2τ_c), O(h): this is how Cor 4.3(c) is consistent with Thm 3.4. It is stated only in S2.1. | After Cor 4.3(c); §1.1 l. 61. |
| D3en-5 | minor | Thm 3.4(b)'s hypotheses: Katzenberger's non-gradient form is supported by Parsons–Rogers and Durrett–Popovic [paper]. C² of Φ̃ likely follows from Falconer Thm 5.1 on the time-1 map (as LWA Lemma B.1) or [Kat91, Cor 5.1] (as SSP). | Remark 2.3, Thm 3.4(b), Table 1: cite; try to drop "assumed". |
| D3en-6 | minor | Prop 3.1's no-diffusion is known for M = I (LWA Cor 5.2; SSP Cor 5.11) [paper]. | l. 224 or §1.3: credit. |
| D3en-7 | minor | "Whatever/every temporal spectrum" overclaims (proved for white and rational spectra; formal otherwise); Thm 3.4's closing sentence covers (c). | l. 61, 101, 256. |
| D3en-8 | minor | Abstract: N-atom competitors "remain local minima" holds only for Q_∞ < 1 (0.30–0.56 of the ends). | Abstract. |
| D3en-9 | minor | Lemma 3.3 assumes analytic 𝗀, but (c) uses any ϱ. | l. 242: "any function on the spectrum". |
| D3en-10 | minor | "Isotropic parameter noise does not select" contradicts LWA Cor 5.1 (−½∂Φ∇ln\|∇²L\|₊) [paper]. | l. 101. |
| D3en-11 | minor | Lemma 3.3 is missing from Table 1. | Table 1. |
| D3en-12 | minor | The supplement's propositions have no status labels; S1.2 is partly sketched. | Supplement: a status table. |
| D3en-13 | minor | The abstract is ≈346 words (SIAM: 250). | Abstract. |
| D3en-14 | minor | "A quarter at hλ = 1" extrapolates the leading order (MC: 12.9% at 0.97). | Remark 4.5; S2.2(d): "for hλ ≪ 1". |
| D3en-15 | minor | New overloads: Π ×3, c ×3, b ×2, η ×2, ν ×3, E/𝔼. | Rename; Table 2. |
| D3en-16 | minor | Remark S1.4 calls the network's multistability "of Prop 5.5", which is an analogy. | S1, Remark S1.4. |
| D3en-17 | minor | SSP's limit is called "a gradient flow of an explicit functional"; that is only their non-degenerate case. | l. 99. |

## Sources
- [read] `manuscript/paper15/paper15.tex`, `supplement.tex`, `README.md`, `numcheck15.py`; `review/REVIEW_round142_D3ek.md`; `review/REVIEW_LOG.md` Round 142; `experiments/exp_colour_selection.py` (constants, `k_stats`, `k_job`, `vertex_limit`, `main`); `src/core/colour_selection.py` (`weights`, `potential`).
- [md]/[raw] `results/chaos/colour_selection.md` (Parts K1, K2) and `colour_selection_raw.npz`; `results/chaos/gauge_label.md`.
- [own] Scripts in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/scripts/` (`k2_check.py`, `aug_check.py`, `ar1_num.py`). The regenerated report is `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/out_cs/colour_selection.md`. The build and rendered figures are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/build/`.
- [paper] Downloaded to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref145/dl/`, text extracted with pdftotext:
  - [SSP, arXiv:2404.12293](https://arxiv.org/abs/2404.12293);
  - [Li–Wang–Arora, arXiv:2110.06914](https://arxiv.org/abs/2110.06914);
  - [Ziyin–Wang–Li–Wu, arXiv:2402.07193](https://arxiv.org/abs/2402.07193);
  - [Parsons–Rogers, arXiv:1510.07031](https://arxiv.org/pdf/1510.07031) (Katzenberger restated for a general drift);
  - [Durrett–Popovic, arXiv:0902.4780](https://arxiv.org/pdf/0902.4780) (Katzenberger applied to a non-gradient system).
- Not reachable: Katzenberger 1991 (projecteuclid.org, 403) and jmlr.org.

## Lesson
A hypothesis written as a closed formula for "all supports" should be tested on the simplest pair, two overlapping supports, before it is printed; here that one test shows the clause is never satisfiable.

When a table reports a quantity at several parameter values, the value at which the experiment was actually run must be among those printed. Prop 5.7's "basin failures" holds at τ_c ≥ 10⁴ but only partly at the τ_c = 100 the flow used.

A full text that was once unreachable should be read and compared result by result as soon as it can be. Here it settles the priority question in the paper's favour.
