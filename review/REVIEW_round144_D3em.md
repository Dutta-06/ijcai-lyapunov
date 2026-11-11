# Referee report D3em, round 144: the Paper 16 manuscript, second review

## Scope
I reviewed `manuscript/paper16/paper16.tex` (681 lines, 25 pp.), `supplement.tex` (7 pp.), `README.md` ("Revision after round 141"), `make_figs.py` and `numcheck16.py` at commit f8ba0a1. I also read the round-141 report and the Round 141 entry of `REVIEW_LOG.md`.

The commit was exported to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref144/repo` with the venv symlinked. The live repository was not touched. Line numbers below are those of `paper16.tex` unless "S" (supplement) is given.

## Status of the round-141 issues

| D3ej id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 (Prop. 4.9) | yes | [own] The proof now uses the certificate to show only that every nonnegative fit vanishes off T. Uniqueness comes from the rank of Φ_T; without it the limit is the Bregman projection onto the face (Thm 4.1, since min_K J = 0). The counterexample (l. 394) checks: h = (0, 0, 1.5) and Φᵀ(1, 1) = (2, 2, 2.5). The Gaussian-translate argument is right: Φ_T = D_x[e^{x_k c_i/ω²}]D_c, and e^{xc} is strictly totally positive. | "the convergence is algebraic" (l. 387) is asserted without the line that proves it. Under the rank condition, \|w_T − w*_T\| ≤ ‖Φ_T⁺‖(\|Φw − y\| + \|Φ_{Tᶜ}w_{Tᶜ}\|) = O(t^{−1/2}). Add this line and say "at least O(t^{−1/2})". |
| 2 (sign of ΔQ) | yes | [own] By hand: bracket = G(2ςt² − β⁰ − (2ς − β⁰)t⁴) = −G(1 − t²)(β⁰ − (2ς − β⁰)t²). The stated region is correct (Q ≥ 0; the core \|κ′t\| < 1; zero set ⇔ γa_max³\|t\|³ = Gm\|w\|(1 − t⁴); positive on the critical curve and on branches). The core drift dβ⁰/dτ ≈ −(2\|φ\|²γG/(mN))(β⁰)² checks. Numerically (`own_euler.py`): the exact Euler ΔQ matches the formula with O(h) relative error at 46 points, and the signs are as stated. Near the corner the h³ term is O(hα) relative to the h² term (∝ α⁵ against α⁴), so "to leading order in h" is safe. | none |
| 3 (support count) | yes | [raw→md] `mirror_leak.md` regenerates byte-identically from its raw file (ML_OUT, `--from-raw`, `cmp`). Every count in §7.5 l. 567 matches md E (416; 1760/1760; 437; 415; 1717; 15 = 8 + 7; IR 0.9930–1.0035; 0.00197 / 2.78·10⁻⁴). [own, partial] My independent dual-Newton equilibrium and FISTA + active-set LASSO, on my own copy of the orbit Gram (eigenvalues 1.417·10⁻⁴ to 0.5533, as in md D), give an equilibrium–LASSO distance of 3.68·10⁻⁵ for the D(1) truth at μ = 10⁻³. That is the md's value; KKT residual 8.9·10⁻¹⁶. My full 1760-pattern rerun did not finish in the time available (see Verification). | Two points: see D3em-3 and D3em-4. |
| 4 (literature) | yes | All the missing works are now cited, and the novelty claim at l. 111 has been narrowed. | The description of PVRF22 is slightly inaccurate (D3em-5). |
| 5 (Cor. 5.3 range) | yes | [own] Rest point ⇒ bracket = 0 ⇒ (2ς − β⁰)t² = β⁰, which needs 0 < β⁰ < ς and t ≠ 0. Katzenberger's range is stated (l. 432, l. 446), and so is Table 1. | D44.md is still unreviewed, which the paper says. |
| 6 (Itô chart) | yes | [own] In (α, β): ½Q_αα(dα)² = −½γGβ²(1 − t⁴) d⟨η⟩. The Euler step in (α, β) gives the same drift, checked numerically to O(h). | none |
| 7 (Thm 4.1(c) on the critical curve) | yes | [own] η ≈ −s·a_max²/(2Gς\|α\|) and β ≈ 2ςα²/a_max² give β ≈ a_max²/(2G²ςg²t²). | none |
| 8 (Derived st. 4.4 error) | yes | [own] The log-derivative is −(π/2)tan − (3π/2)cot, so the relative error is O(ε^{1/2}/(μ(1 − μ))). I re-derived the profile itself from κ′t = cot(πμ/2). | README l. 33 still says "O(ε^{1/2}/μ)". |
| 9 (wording) | yes | l. 358 | none |
| 10 (Prop. 6.4 scalings) | yes | [own] (a): ρ = k_eσm/(G(β⁰)²) gives (H + ρI)⁻¹Hw*. (b): μ̄ = σΛ(ε) is fixed and σ → 0; the convergence is uniform on compacts by Prop. 3.3(d). | none |
| 11 (Prop. 6.2) | yes | [own] I re-derived the V̇ identity, the three-point identity, V̇ ≤ −λV + C, the noise term n(w − w*)ᵀφ ≤ σB(w, w_p) + C_n (finite because B(·, w_p) grows like \|w\|log\|w\|), and d(t) → 0 in the example (growth ratios 3.30, 3.19, 3.16 ≈ 100^{1/4}). | φ̄ must be a sup-norm bound on the components to pair with \|·\|₁ (D3em-9). |
| 12 (Prop. 3.4(e)) | yes | [own] The relative error is θ⁴β(0)/(6ς(1 − θ²)); I checked the arcsin tilt. | none |
| 13 (cosmetic) | yes | — | none |
| 14 (Fig. 3, SRR) | yes | Fig. 3 rendered: the guide line is clipped to log₁₀(1/ε) ∈ [16, 64]. | At ε = 10⁻⁶⁴ the guide lies about 13× above the data (local slope 0.115–0.127 against μ_min = 0.1039). Optional caption note: "asymptotic slope". |
| 15 (length, motivation) | mostly | 25 pp.; Sections 1–8 end on p. 22; Table 1 marks the ten supplement proofs with †. The S6 references (Prop. 3.2, (3.1)–(3.3), Remark 4.5/4.7/6.5, (4.1)) all match my build's .aux. | The new §1.1 text has an inaccuracy (D3em-1). The supplement carries review history and an inaccurate "no statement uses" (D3em-7). |

## New findings

### Major
None.

### Minor

1. **l. 66 (new §1.1): "in the coordinates (α, 2β^{1/2}) the pair is the plain gradient flow".** [own] With u = 2β^{1/2}, α̇ = −G∂_αL and u̇ = −γ∂_uL. That is a diagonally preconditioned gradient flow; it is plain only if G = γ, and the defaults have G = 1, γ = 0.1. §2 (l. 123) states the metric correctly.
2. **l. 66: "kernels whose frozen sign cannot help lose their mass exponentially".** The decay is only algebraic on the critical curve (Thm 4.1(c), as revised), and no rate is given when g_i(∞) = 0. Add "on branches, when g_i(∞) ≠ 0".
3. **l. 567: the 8 mismatches attributed to "the core's profile, of order (β⁰)^{3/2}μ_j^{−3}, Derived statement 4.4".** The md does not support this.
   - md E reports the tangent profile at the LASSO's slacks as up to 1.75 of the largest coefficient, against the equilibrium's actual ≤ 0.00372, with a largest relative error of 204.
   - [own, analytic] With the md's constant a_max(β⁰)^{3/2}/(m(2ς)^{1/2}) = 3.73·10⁻¹⁰, the profile at slack 0.00134 is about 0.04 absolute (~0.15 relative). The equilibrium has 0.00108–0.00348.
   - The reason is feedback. An inactive coordinate's own core coefficient moves its slack by about H_jj·w_j/μ, with μ = 10⁻³. So the formula at the LASSO's slack applies only when μ_j⁴ ≫ (core·H_jj/μ), that is μ_j ≳ 10⁻². The 8 patterns (slacks 0.00134–0.007) are outside that range.
   - Fix: say that the coefficient is the core's, set self-consistently, and drop "of order (β⁰)^{3/2}μ_j^{−3}". Alternatively, evaluate the profile at the equilibrium's own slack.
4. **l. 567: "the direct check of Proposition 6.4(b)".** This is a single β⁰ = 10⁻⁶, so it is a consistency check, not a check of a limit.
   - The largest distance, 0.00197, is about a quarter of the LASSO's own shrinkage at this penalty (max\|w − w*\| = 0.00752 for the D(1) truth).
   - Either report the distance at two or more values of β⁰ (D(1) already shows a linear dependence on μ), or reword as "consistent with, at β⁰ = 10⁻⁶".
5. **l. 111, PVRF22 [paper].** Three inaccuracies:
   - Their model is β = θ⊙θ, a nonnegative quadratic parametrisation, not u⊙v (signed u² − v² is only remarked on).
   - The selected object is a nonnegative Lasso weighted by h = diag(XᵀX), whose penalty δ is proportional to step × noise variance.
   - The flow does not converge to it. Off the support, coordinates vanish exponentially at fixed noise; on the support, the flow fluctuates within O(√δ).
   Two consequences for the paper. First, PVRF's analysis is non-asymptotic at fixed noise and covers exactly the extinction regime that the paper says lies outside the reduced-flow theorem; it should be cited in Open problem (4). Second, since Prop. 5.2(b) also carries \|φ_i\|², the paper's Euler ends should be compared with the \|φ_i\|²-weighted basis pursuit (the analogue of PVRF's weights) before being called "a reweighting".
6. **l. 111, KSGYT21 [paper].** Their Appendix F also derives, by Itô calculus, the noise forcing of a conserved quantity (for scale symmetry, "analogous" for rescale). Cite it for Prop. 5.2(b) as well, not only for finite steps. Neither PVRF nor Kunin et al. contains the paper's claimed new items: the sign of ΔQ on every curve for the saturating unit (no linear symmetry group), the chart dependence, and the exact on-curve update. Prop. 5.2(c) is stochastic mirror descent on Ψ_Q, and its invariance is the cited GLSS18/AH19 result; say so explicitly.
7. **Supplement: review history and an inaccurate disclaimer.**
   - S l. 75 reads "round 141 moved them here for length", which is review history in a submission (README l. 23 claims there is none). S l. 5 contains a literal `\n` inside a comment.
   - Paper l. 528 and S l. 32 say that no statement uses S1–S5. Yet §7 rests claims on them: l. 539 (unequal masses, S2), l. 544 ("Where S is not a point…", S3), l. 554 (S3), l. 557 (Heun, S4), l. 559 (S4), and Remark 6.5 l. 523 (S5).
   - Fix: say "S1–S5 hold the runs behind the statements marked (S·)".
8. **Remark 4.7 (l. 361) against S3 and `mirror_branch.md`.** The remark predicts free coefficients of order ε²log(1/ε). The md's Σ\|w_free\|/ε² is 251, 784, 782 at ε = 10⁻⁴, 10⁻⁸, 10⁻¹². It is flat between 10⁻⁸ and 10⁻¹², where log(1/ε) grows by a factor 1.5, and S3 says "about ε²". Either the log factor cancels (that is, (Φᵀξ)_j = 0 for the free atoms, or z_⊥ compensates when the branch support has fewer than N atoms) or the remark is wrong. I did not resolve this numerically (a script is prepared, not run).
9. **Prop. 6.2 (l. 485):** state that \|φ(x(t))\|_∞ ≤ φ̄.
10. **README l. 33:** the stale "O(ε^{1/2}/μ)".

## Verification I did
- **[own] Re-derivations by hand:**
  - Prop. 4.9 with the rank condition, its counterexample, and the total-positivity argument.
  - Prop. 5.2(a): expansion, factorisation, sign region, zero set; (b) and its τ-scaling; the core 1/τ law.
  - Cor. 5.3.
  - Itô in (α, β) and in (α, log β); the Euler step in (α, β).
  - Thm 4.1(c) on the critical curve; existence via ψ*.
  - Derived st. 4.4: the profile and its error.
  - Prop. 6.2: identity, comparison bound, bounded noise, example.
  - Prop. 6.3: H_L + H_Lᵀ, the 1/k_e expansion, the root; Prop. 6.4(a) and (b); Remark 5.4; Remark 6.5.
  - S6 proofs of Props. 3.2 (all readings, including coth), 3.3(c) (the log(2β⁰/ς) term), 3.4(b), (c), (e), 4.6(a), (b) and 4.8(b): complete and correctly referenced.
- **[own] `own_euler.py`:** the exact Euler ΔQ in both charts against the formulas; signs on Q/ς ∈ {40, 1, 0, −0.5, −0.999, −1, −2, −10}.
- **[raw→md]** `mirror_leak.md` regenerated from `mirror_leak_raw.npz` into ML_OUT and compared with `cmp`: IDENTICAL.
- **[own, partial] `own_partE.py`:** my own Gram, a dual-Newton equilibrium on D(η) = Σψ* − w*ᵀη + (σ/2)ηᵀH⁻¹η, and FISTA + active-set LASSO.
  - λ₁ = 7024.74 and σ = 1.4235·10⁻⁷.
  - The D(1) truth reproduces 3.68·10⁻⁵ (the md's value).
  - The 1760-pattern rerun was stopped unfinished: the first version's solver tolerance was too tight, and the patched version was at pattern 0–100 when handback was required. Both were my own PIDs, verified via /proc/PID/cmdline.
  - So the 1760-pattern counts are verified only through [raw→md], not recomputed.
- **[own run] `numcheck16.py`:** 292 tokens, 1 missing (10.4, derived in the text). Note that it checks only decimals: the integers 1760, 416, 437, 415, 1717, 15, 8 and 7 were checked by hand against md D and E. The cross-reference masking hides only references.
- **Build:** pdflatex ×3 for the paper and ×2 for the supplement, no errors and no overfull boxes. 25 + 7 pp. `pdftotext` is identical to the committed PDFs; the .aux numbering matches the supplement's references.
- **Fig. 3:** rendered and inspected.

## Literature
- **[paper] Pillaud-Vivien–Reygner–Flammarion, COLT 2022 (arXiv:2206.09841), read in full text** (§§1–6, abstract, theorems 2–4). It is the closest prior result. It does not contain the paper's results (no conserved Q, no saturating unit, no rest-point condition). Its description in the paper needs the corrections in D3em-5, and it should be pointed to for the extinction regime.
- **[paper] Kunin et al., ICLR 2021 (arXiv:2012.04728), read in full text** (abstract, §§1–6, Apps E–F). It gives exact broken-conservation dynamics for translation, scale and rescale symmetries. For rescale it notes that the discretisation force can have either sign, and it treats noise forcing by Itô calculus. It does not cover non-group conserved quantities such as Q, so the paper's narrowed novelty claim (sign on every curve, chart dependence, exact update) is fair. Add the App. F credit.
- The other comparisons I did not re-read.

## Overall assessment
**Recommendation: minor revision.**

All five round-141 majors are resolved correctly, and I re-derived every changed statement. The proofs moved to S6 are complete and correctly referenced, Table 1 is honest, the build and the number check are clean, and the md reproduces from its raw file. What remains is wording and the reading of numbers:
- §1.1's "plain gradient flow";
- the core-profile attribution and the "direct check" in §7.5;
- the precise description of PVRF;
- Remark 4.7's log factor;
- the supplement's review history and disclaimer;
- one missing proof line in Prop. 4.9.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3em-1 | minor | "plain gradient flow in (α, 2β^{1/2})" holds only for G = γ [own] | Say "the gradient flow with gains G, γ" (l. 66). |
| D3em-2 | minor | "lose their mass exponentially": the decay is algebraic on the critical curve, and no rate is given if g(∞) = 0 | Qualify the claim (l. 66). |
| D3em-3 | minor | The 8 mismatches attributed to the profile "(β⁰)^{3/2}μ_j^{−3}" at the LASSO's slacks: the md shows that profile off by up to 204×, because of self-consistent feedback [md E; own analytic] | Reword, or use the equilibrium's own slack (l. 567). |
| D3em-4 | minor | A single β⁰ is called "the direct check" of a limit; the distance (0.00197) is about ¼ of the LASSO's shrinkage | Add a β⁰ sequence or say "consistent with" (l. 567). |
| D3em-5 | minor | PVRF: θ⊙θ nonnegative; an h-weighted Lasso with penalty ∝ step × noise; O(√δ) fluctuation; it treats extinction at fixed noise [paper] | Fix l. 111; cite in §8 (4); compare the Euler ends with \|φ_i\|²-weighted basis pursuit (§7.4). |
| D3em-6 | minor | Kunin App. F also has the Itô forcing of conserved quantities; Prop. 5.2(c) is SMD invariance | Credit both at l. 111. |
| D3em-7 | minor | Supplement l. 75 "round 141"; l. 5 literal `\n`; "no statement uses S1–S5" is false (l. 528, 539, 544, 554, 557, 559, 523) | Remove the history; restate the disclaimer. |
| D3em-8 | minor | Remark 4.7's ε²log(1/ε) against the md's flat Σ\|w_free\|/ε² (784 → 782 from 10⁻⁸ to 10⁻¹²) [md mirror_branch] | Re-derive; correct the remark or S3. |
| D3em-9 | cosmetic | φ̄'s norm is unspecified in Prop. 6.2 | Use \|φ\|_∞ (l. 485). |
| D3em-10 | minor | "the convergence is algebraic" is unproved in the text | Add the O(t^{−1/2}) line via Φ_T⁺ (l. 387–391). |
| D3em-11 | cosmetic | README l. 33's stale error order; Fig. 3's guide sits about 13× above the data at 10⁻⁶⁴ | Edit; add "asymptotic slope" to the caption. |

## Sources
- **[read]** `paper16.tex`, `supplement.tex`, `README.md`, `numcheck16.py`, `make_figs.py` (Fig. 3 code), `exp_mirror_leak.py` (header, Parts D and E, report), `tests/test_mirror_leak.py`, `exp_mirror_branch.py` ll. 182–240, `REVIEW_round141_D3ej.md`, `REVIEW_LOG.md` Round 141.
- **[md]** `mirror_leak.md` D and E; `mirror_branch.md` (mixed starts); `mirror_free.md` C(3).
- **[raw]** `mirror_leak_raw.npz`.
- **[own]** In `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref144/`: `own_euler.py`, `own_partE.py`, `time1.py`, `own_mixed.py` (not run), and the build in `build/`.
- **[paper]** arXiv:2206.09841 and arXiv:2012.04728. The downloads are in `dl_pvrf/` and `dl_kunin/` in the same scratchpad directory; text extracted only, nothing executed.

## Lesson
When a revision explains a residual by a closed-form asymptotic profile, evaluate the profile at the actual data before citing it. Here the md's own column showed that the profile at the LASSO's slacks was off by up to 204×, because the coefficient feeds back into its own slack. Also, give a long job a progress line and a convergence criterion at the solver's floor, so it can finish within the review.
